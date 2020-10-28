/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "noise_model.h"
#include "noise_util.h"
#include "mathutils.h"
#include "EbLog.h"

#define kLowPolyNumParams 3

static const int32_t k_max_lag = 4;

void *svt_aom_memalign(size_t align, size_t size);
void  svt_aom_free(void *memblk);

void un_pack2d(uint16_t *in16_bit_buffer, uint32_t in_stride, uint8_t *out8_bit_buffer,
               uint32_t out8_stride, uint8_t *outn_bit_buffer, uint32_t outn_stride, uint32_t width,
               uint32_t height);

void pack2d_src(uint8_t *in8_bit_buffer, uint32_t in8_stride, uint8_t *inn_bit_buffer,
                uint32_t inn_stride, uint16_t *out16_bit_buffer, uint32_t out_stride,
                uint32_t width, uint32_t height);

// Defines a function that can be used to obtain the mean of a block for the
// provided data type (uint8_t, or uint16_t)
#define GET_BLOCK_MEAN(INT_TYPE, suffix)                          \
    static double get_block_mean_##suffix(const INT_TYPE *data,   \
                                          int32_t         w,      \
                                          int32_t         h,      \
                                          int32_t         stride, \
                                          int32_t         x_o,    \
                                          int32_t         y_o,    \
                                          int32_t         block_size) {   \
        const int32_t max_h      = AOMMIN(h - y_o, block_size);   \
        const int32_t max_w      = AOMMIN(w - x_o, block_size);   \
        double        block_mean = 0;                             \
        for (int32_t y = 0; y < max_h; ++y) {                     \
            for (int32_t x = 0; x < max_w; ++x) {                 \
                block_mean += data[(y_o + y) * stride + x_o + x]; \
            }                                                     \
        }                                                         \
        return block_mean / (max_w * max_h);                      \
    }

GET_BLOCK_MEAN(uint8_t, lowbd);
GET_BLOCK_MEAN(uint16_t, highbd);

static INLINE double get_block_mean(const uint8_t *data, int32_t w, int32_t h, int32_t stride,
                                    int32_t x_o, int32_t y_o, int32_t block_size,
                                    int32_t use_highbd) {
    if (use_highbd)
        return get_block_mean_highbd((const uint16_t *)data, w, h, stride, x_o, y_o, block_size);
    return get_block_mean_lowbd(data, w, h, stride, x_o, y_o, block_size);
}

// Defines a function that can be used to obtain the variance of a block
// for the provided data type (uint8_t, or uint16_t)
#define GET_NOISE_VAR(INT_TYPE, suffix)                                     \
    static double get_noise_var_##suffix(const INT_TYPE *data,              \
                                         const INT_TYPE *denoised,          \
                                         int32_t         stride,            \
                                         int32_t         w,                 \
                                         int32_t         h,                 \
                                         int32_t         x_o,               \
                                         int32_t         y_o,               \
                                         int32_t         block_size_x,      \
                                         int32_t         block_size_y) {            \
        const int32_t max_h      = AOMMIN(h - y_o, block_size_y);           \
        const int32_t max_w      = AOMMIN(w - x_o, block_size_x);           \
        double        noise_var  = 0;                                       \
        double        noise_mean = 0;                                       \
        for (int32_t y = 0; y < max_h; ++y) {                               \
            for (int32_t x = 0; x < max_w; ++x) {                           \
                double noise = (double)data[(y_o + y) * stride + x_o + x] - \
                               denoised[(y_o + y) * stride + x_o + x];      \
                noise_mean += noise;                                        \
                noise_var += noise * noise;                                 \
            }                                                               \
        }                                                                   \
        noise_mean /= (max_w * max_h);                                      \
        return noise_var / (max_w * max_h) - noise_mean * noise_mean;       \
    }

GET_NOISE_VAR(uint8_t, lowbd);
GET_NOISE_VAR(uint16_t, highbd);

static INLINE double get_noise_var(const uint8_t *data, const uint8_t *denoised, int32_t w,
                                   int32_t h, int32_t stride, int32_t x_o, int32_t y_o,
                                   int32_t block_size_x, int32_t block_size_y, int32_t use_highbd) {
    if (use_highbd)
        return get_noise_var_highbd((const uint16_t *)data,
                                    (const uint16_t *)denoised,
                                    w,
                                    h,
                                    stride,
                                    x_o,
                                    y_o,
                                    block_size_x,
                                    block_size_y);
    return get_noise_var_lowbd(data, denoised, w, h, stride, x_o, y_o, block_size_x, block_size_y);
}

static void equation_system_free(AomEquationSystem *eqns) {
    if (!eqns) return;
    free(eqns->A);
    eqns->A = NULL;
    free(eqns->b);
    eqns->b = NULL;
    free(eqns->x);
    eqns->x = NULL;
    eqns->n = 0;
}

static void equation_system_clear(AomEquationSystem *eqns) {
    const int32_t n = eqns->n;
    memset(eqns->A, 0, sizeof(*eqns->A) * n * n);
    memset(eqns->x, 0, sizeof(*eqns->x) * n);
    memset(eqns->b, 0, sizeof(*eqns->b) * n);
}

static void equation_system_copy(AomEquationSystem *dst, const AomEquationSystem *src) {
    const int32_t n = dst->n;
    if (svt_memcpy != NULL)
    {
        svt_memcpy(dst->A, src->A, sizeof(*dst->A) * n * n);
        svt_memcpy(dst->x, src->x, sizeof(*dst->x) * n);
        svt_memcpy(dst->b, src->b, sizeof(*dst->b) * n);
    }
    else
    {
        svt_memcpy_c(dst->A, src->A, sizeof(*dst->A) * n * n);
        svt_memcpy_c(dst->x, src->x, sizeof(*dst->x) * n);
        svt_memcpy_c(dst->b, src->b, sizeof(*dst->b) * n);
    }
}

static int32_t equation_system_init(AomEquationSystem *eqns, int32_t n) {
    eqns->A = (double *)malloc(sizeof(*eqns->A) * n * n);
    eqns->b = (double *)malloc(sizeof(*eqns->b) * n);
    eqns->x = (double *)malloc(sizeof(*eqns->x) * n);
    eqns->n = n;
    if (!eqns->A || !eqns->b || !eqns->x) {
        SVT_ERROR("Failed to allocate system of equations of size %d\n", n);
        equation_system_free(eqns);
        return 0;
    }
    equation_system_clear(eqns);
    return 1;
}

static int32_t equation_system_solve(AomEquationSystem *eqns) {
    const int32_t n   = eqns->n;
    double *      b   = (double *)malloc(sizeof(*b) * n);
    double *      A   = (double *)malloc(sizeof(*A) * n * n);
    int32_t       ret = 0;
    if (A == NULL || b == NULL) {
        SVT_ERROR("Unable to allocate temp values of size %dx%d\n", n, n);
        free(b);
        free(A);
        return 0;
    }
    if (svt_memcpy != NULL)
    {
        svt_memcpy(A, eqns->A, sizeof(*eqns->A) * n * n);
        svt_memcpy(b, eqns->b, sizeof(*eqns->b) * n);
    }
    else
    {
        svt_memcpy_c(A, eqns->A, sizeof(*eqns->A) * n * n);
        svt_memcpy_c(b, eqns->b, sizeof(*eqns->b) * n);
    }
    ret = linsolve(n, A, eqns->n, b, eqns->x);
    free(b);
    free(A);

    if (ret == 0) return 0;
    return 1;
}
/*
static void equation_system_add(AomEquationSystem *dest,
                                AomEquationSystem *src) {
  const int32_t n = dest->n;
  int32_t i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j)
      dest->A[i * n + j] += src->A[i * n + j];
    dest->b[i] += src->b[i];
  }
}
*/

static void noise_strength_solver_clear(AomNoiseStrengthSolver *solver) {
    equation_system_clear(&solver->eqns);
    solver->num_equations = 0;
    solver->total         = 0;
}

/*
static void noise_strength_solver_add(AomNoiseStrengthSolver *dest,
                                      AomNoiseStrengthSolver *src) {
  equation_system_add(&dest->eqns, &src->eqns);
  dest->num_equations += src->num_equations;
  dest->total += src->total;
}
*/

static void noise_strength_solver_copy(AomNoiseStrengthSolver *dest, AomNoiseStrengthSolver *src) {
    equation_system_copy(&dest->eqns, &src->eqns);
    dest->num_equations = src->num_equations;
    dest->total         = src->total;
}

// Return the number of coefficients required for the given parameters
static int32_t num_coeffs(const AomNoiseModelParams params) {
    const int32_t n = 2 * params.lag + 1;
    switch (params.shape) {
    case AOM_NOISE_SHAPE_DIAMOND: return params.lag * (params.lag + 1);
    case AOM_NOISE_SHAPE_SQUARE: return (n * n) / 2;
    }
    return 0;
}

static int32_t noise_state_init(AomNoiseState *state, int32_t n, int32_t bit_depth) {
    const int32_t k_num_bins = 20;
    if (!equation_system_init(&state->eqns, n)) {
        SVT_ERROR("Failed initialization noise state with size %d\n", n);
        return 0;
    }
    state->ar_gain          = 1.0;
    state->num_observations = 0;
    return svt_aom_noise_strength_solver_init(&state->strength_solver, k_num_bins, bit_depth);
}

static void set_chroma_coefficient_fallback_soln(AomEquationSystem *eqns) {
    const double  k_tolerance = 1e-6;
    const int32_t last        = eqns->n - 1;
    // Set all of the AR coefficients to zero, but try to solve for correlation
    // with the luma channel
    memset(eqns->x, 0, sizeof(*eqns->x) * eqns->n);
    if (fabs(eqns->A[last * eqns->n + last]) > k_tolerance)
        eqns->x[last] = eqns->b[last] / eqns->A[last * eqns->n + last];
}

int32_t svt_aom_noise_strength_lut_init(AomNoiseStrengthLut *lut, int32_t num_points) {
    if (!lut) return 0;
    lut->points = (double(*)[2])malloc(num_points * sizeof(*lut->points));
    if (!lut->points) return 0;
    lut->num_points = num_points;
    memset(lut->points, 0, sizeof(*lut->points) * num_points);
    return 1;
}

void svt_aom_noise_strength_lut_free(AomNoiseStrengthLut *lut) {
    if (!lut) return;
    free(lut->points);
    memset(lut, 0, sizeof(*lut));
}

static double noise_strength_solver_get_bin_index(const AomNoiseStrengthSolver *solver,
                                                  double                        value) {
    const double val   = fclamp(value, solver->min_intensity, solver->max_intensity);
    const double range = solver->max_intensity - solver->min_intensity;
    return (solver->num_bins - 1) * (val - solver->min_intensity) / range;
}

static double noise_strength_solver_get_value(const AomNoiseStrengthSolver *solver, double x) {
    const double  bin    = noise_strength_solver_get_bin_index(solver, x);
    const int32_t bin_i0 = (int32_t)floor(bin);
    const int32_t bin_i1 = AOMMIN(solver->num_bins - 1, bin_i0 + 1);
    const double  a      = bin - bin_i0;
    return (1.0 - a) * solver->eqns.x[bin_i0] + a * solver->eqns.x[bin_i1];
}

void svt_aom_noise_strength_solver_add_measurement(AomNoiseStrengthSolver *solver, double block_mean,
                                                  double noise_std) {
    const double  bin    = noise_strength_solver_get_bin_index(solver, block_mean);
    const int32_t bin_i0 = (int32_t)floor(bin);
    const int32_t bin_i1 = AOMMIN(solver->num_bins - 1, bin_i0 + 1);
    const double  a      = bin - bin_i0;
    const int32_t n      = solver->num_bins;
    solver->eqns.A[bin_i0 * n + bin_i0] += (1.0 - a) * (1.0 - a);
    solver->eqns.A[bin_i1 * n + bin_i0] += a * (1.0 - a);
    solver->eqns.A[bin_i1 * n + bin_i1] += a * a;
    solver->eqns.A[bin_i0 * n + bin_i1] += a * (1.0 - a);
    solver->eqns.b[bin_i0] += (1.0 - a) * noise_std;
    solver->eqns.b[bin_i1] += a * noise_std;
    solver->total += noise_std;
    solver->num_equations++;
}

int32_t svt_aom_noise_strength_solver_solve(AomNoiseStrengthSolver *solver) {
    // Add regularization proportional to the number of constraints
    const int32_t n       = solver->num_bins;
    const double  k_alpha = 2.0 * (double)(solver->num_equations) / n;
    int32_t       result  = 0;
    double        mean    = 0;

    // Do this in a non-destructive manner so it is not confusing to the caller
    double *old_a = solver->eqns.A;
    double *A     = (double *)malloc(sizeof(*A) * n * n);
    if (!A) {
        SVT_ERROR("Unable to allocate copy of A\n");
        return 0;
    }
    if (svt_memcpy != NULL)
        svt_memcpy(A, old_a, sizeof(*A) * n * n);
    else
        svt_memcpy_c(A, old_a, sizeof(*A) * n * n);

    for (int32_t i = 0; i < n; ++i) {
        const int32_t i_lo = AOMMAX(0, i - 1);
        const int32_t i_hi = AOMMIN(n - 1, i + 1);
        A[i * n + i_lo] -= k_alpha;
        A[i * n + i] += 2 * k_alpha;
        A[i * n + i_hi] -= k_alpha;
    }

    // Small regularization to give average noise strength
    mean = solver->total / solver->num_equations;
    for (int32_t i = 0; i < n; ++i) {
        A[i * n + i] += 1.0 / 8192.;
        solver->eqns.b[i] += mean / 8192.;
    }
    solver->eqns.A = A;
    result         = equation_system_solve(&solver->eqns);
    solver->eqns.A = old_a;

    free(A);
    return result;
}

int32_t svt_aom_noise_strength_solver_init(AomNoiseStrengthSolver *solver, int32_t num_bins,
                                          int32_t bit_depth) {
    if (!solver) return 0;
    memset(solver, 0, sizeof(*solver));
    solver->num_bins      = num_bins;
    solver->min_intensity = 0;
    solver->max_intensity = (1 << bit_depth) - 1;
    solver->total         = 0;
    solver->num_equations = 0;
    return equation_system_init(&solver->eqns, num_bins);
}

double svt_aom_noise_strength_solver_get_center(const AomNoiseStrengthSolver *solver, int32_t i) {
    const double  range = solver->max_intensity - solver->min_intensity;
    const int32_t n     = solver->num_bins;
    return ((double)i) / (n - 1) * range + solver->min_intensity;
}

// Computes the residual if a point were to be removed from the lut. This is
// calculated as the area between the output of the solver and the line segment
// that would be formed between [x_{i - 1}, x_{i + 1}).
static void update_piecewise_linear_residual(const AomNoiseStrengthSolver *  solver,
                                             const AomNoiseStrengthLut *lut, double *residual,
                                             int32_t start, int32_t end) {
    const double dx = 255. / solver->num_bins;
    for (int32_t i = AOMMAX(start, 1); i < AOMMIN(end, lut->num_points - 1); ++i) {
        const int32_t lower = AOMMAX(
            0, (int32_t)floor(noise_strength_solver_get_bin_index(solver, lut->points[i - 1][0])));
        const int32_t upper = AOMMIN(
            solver->num_bins - 1,
            (int32_t)ceil(noise_strength_solver_get_bin_index(solver, lut->points[i + 1][0])));
        double r = 0;
        for (int32_t j = lower; j <= upper; ++j) {
            const double x = svt_aom_noise_strength_solver_get_center(solver, j);
            if (x < lut->points[i - 1][0]) continue;
            if (x >= lut->points[i + 1][0]) continue;
            const double y = solver->eqns.x[j];
            const double a =
                (x - lut->points[i - 1][0]) / (lut->points[i + 1][0] - lut->points[i - 1][0]);
            const double estimate_y = lut->points[i - 1][1] * (1.0 - a) + lut->points[i + 1][1] * a;
            r += fabs(y - estimate_y);
        }
        residual[i] = r * dx;
    }
}

int32_t svt_aom_noise_strength_solver_fit_piecewise(const AomNoiseStrengthSolver *solver,
                                                   int32_t                       max_output_points,
                                                   AomNoiseStrengthLut *    lut) {
    // The tolerance is normalized to be give consistent results between
    // different bit-depths.
    const double k_tolerance = solver->max_intensity * 0.00625 / 255.0;
    if (!svt_aom_noise_strength_lut_init(lut, solver->num_bins)) {
        SVT_ERROR("Failed to init lut\n");
        return 0;
    }
    for (int32_t i = 0; i < solver->num_bins; ++i) {
        lut->points[i][0] = svt_aom_noise_strength_solver_get_center(solver, i);
        lut->points[i][1] = solver->eqns.x[i];
    }
    if (max_output_points < 0) max_output_points = solver->num_bins;
    double *residual = malloc(solver->num_bins * sizeof(*residual));
    ASSERT(residual != NULL);
    memset(residual, 0, sizeof(*residual) * solver->num_bins);

    update_piecewise_linear_residual(solver, lut, residual, 0, solver->num_bins);

    // Greedily remove points if there are too many or if it doesn't hurt local
    // approximation (never remove the end points)
    while (lut->num_points > 2) {
        int32_t min_index = 1;
        for (int32_t j = 1; j < lut->num_points - 1; ++j) {
            if (residual[j] < residual[min_index]) min_index = j;
        }
        const double dx           = lut->points[min_index + 1][0] - lut->points[min_index - 1][0];
        const double avg_residual = residual[min_index] / dx;
        if (lut->num_points <= max_output_points && avg_residual > k_tolerance) break;
        const int32_t num_remaining = lut->num_points - min_index - 1;
        memmove(lut->points + min_index,
                lut->points + min_index + 1,
                sizeof(lut->points[0]) * num_remaining);
        lut->num_points--;

        update_piecewise_linear_residual(solver, lut, residual, min_index - 1, min_index + 1);
    }
    free(residual);
    return 1;
}

int32_t svt_aom_flat_block_finder_init(AomFlatBlockFinder *block_finder, int32_t block_size,
                                      int32_t bit_depth, int32_t use_highbd) {
    const int32_t     n = block_size * block_size;
    AomEquationSystem eqns;
    if (!equation_system_init(&eqns, kLowPolyNumParams)) {
        SVT_ERROR("Failed to init equation system for block_size=%d\n", block_size);
        return 0;
    }

    double *at_a_inv = (double *)malloc(kLowPolyNumParams * kLowPolyNumParams * sizeof(*at_a_inv));
    double *A        = (double *)malloc(kLowPolyNumParams * n * sizeof(*A));
    if (at_a_inv == NULL || A == NULL) {
        SVT_ERROR("Failed to alloc A or at_a_inv for block_size=%d\n", block_size);
        free(at_a_inv);
        free(A);
        equation_system_free(&eqns);
        return 0;
    }

    block_finder->A             = A;
    block_finder->at_a_inv      = at_a_inv;
    block_finder->block_size    = block_size;
    block_finder->normalization = (1 << bit_depth) - 1;
    block_finder->use_highbd    = use_highbd;

    for (int32_t y = 0; y < block_size; ++y) {
        const double yd = ((double)y - block_size / 2.) / (block_size / 2.);
        for (int32_t x = 0; x < block_size; ++x) {
            const double  xd               = ((double)x - block_size / 2.) / (block_size / 2.);
            const double  coords[3]        = {yd, xd, 1};
            const int32_t row              = y * block_size + x;
            A[kLowPolyNumParams * row + 0] = yd;
            A[kLowPolyNumParams * row + 1] = xd;
            A[kLowPolyNumParams * row + 2] = 1;

            for (int i = 0; i < kLowPolyNumParams; ++i)
                for (int j = 0; j < kLowPolyNumParams; ++j)
                    eqns.A[kLowPolyNumParams * i + j] += coords[i] * coords[j];
        }
    }

    // Lazy inverse using existing equation solver.
    for (int i = 0; i < kLowPolyNumParams; ++i) {
        memset(eqns.b, 0, sizeof(*eqns.b) * kLowPolyNumParams);
        eqns.b[i] = 1;
        equation_system_solve(&eqns);

        for (int j = 0; j < kLowPolyNumParams; ++j) at_a_inv[j * kLowPolyNumParams + i] = eqns.x[j];
    }
    equation_system_free(&eqns);
    return 1;
}

void svt_aom_flat_block_finder_free(AomFlatBlockFinder *block_finder) {
    if (!block_finder) return;
    free(block_finder->A);
    free(block_finder->at_a_inv);
    memset(block_finder, 0, sizeof(*block_finder));
}

void svt_aom_flat_block_finder_extract_block(const AomFlatBlockFinder *block_finder,
                                            const uint8_t *const data, int32_t w, int32_t h,
                                            int32_t stride, int32_t offsx, int32_t offsy,
                                            double *plane, double *block) {
    const int32_t block_size = block_finder->block_size;
    const int32_t n          = block_size * block_size;
    const double *A          = block_finder->A;
    const double *at_a_inv   = block_finder->at_a_inv;
    double        plane_coords[kLowPolyNumParams];
    double        at_a_inv__b[kLowPolyNumParams];
    int32_t       xi, yi, i;

    if (block_finder->use_highbd) {
        const uint16_t *const data16 = (const uint16_t *const)data;
        for (yi = 0; yi < block_size; ++yi) {
            const int32_t y = clamp(offsy + yi, 0, h - 1);
            for (xi = 0; xi < block_size; ++xi) {
                const int32_t x = clamp(offsx + xi, 0, w - 1);
                block[yi * block_size + xi] =
                    ((double)data16[y * stride + x]) / block_finder->normalization;
            }
        }
    } else {
        for (yi = 0; yi < block_size; ++yi) {
            const int32_t y = clamp(offsy + yi, 0, h - 1);
            for (xi = 0; xi < block_size; ++xi) {
                const int32_t x = clamp(offsx + xi, 0, w - 1);
                block[yi * block_size + xi] =
                    ((double)data[y * stride + x]) / block_finder->normalization;
            }
        }
    }
    multiply_mat(block, A, at_a_inv__b, 1, n, kLowPolyNumParams);
    multiply_mat(at_a_inv, at_a_inv__b, plane_coords, kLowPolyNumParams, kLowPolyNumParams, 1);
    multiply_mat(A, plane_coords, plane, n, kLowPolyNumParams, 1);

    for (i = 0; i < n; ++i) block[i] -= plane[i];
}

typedef struct {
    int32_t index;
    float   score;
} IndexAndscore;

static int compare_scores(const void *a, const void *b) {
    const float diff = ((IndexAndscore *)a)->score - ((IndexAndscore *)b)->score;
    return diff < 0 ? -1 : diff > 0;
}

int32_t svt_aom_flat_block_finder_run(const AomFlatBlockFinder *block_finder,
                                     const uint8_t *const data, int32_t w, int32_t h,
                                     int32_t stride, uint8_t *flat_blocks) {
    // The gradient-based features used in this code are based on:
    //  A. Kokaram, D. Kelly, H. Denman and A. Crawford, "Measuring noise
    //  correlation for improved video denoising," 2012 19th, ICIP.
    // The thresholds are more lenient to allow for correct grain modeling
    // if extreme cases.
    const int32_t  block_size        = block_finder->block_size;
    const int32_t  n                 = block_size * block_size;
    const double   k_trace_threshold = 0.15 / (32 * 32);
    const double   k_ratio_threshold = 1.25;
    const double   k_norm_threshold  = 0.08 / (32 * 32);
    const double   k_var_threshold   = 0.005 / (double)n;
    const int32_t  num_blocks_w      = (w + block_size - 1) / block_size;
    const int32_t  num_blocks_h      = (h + block_size - 1) / block_size;
    int32_t        num_flat          = 0;
    double *       plane  = (double *)malloc(n * sizeof(*plane));
    double *       block  = (double *)malloc(n * sizeof(*block));
    IndexAndscore *scores = (IndexAndscore *)malloc(num_blocks_w * num_blocks_h * sizeof(*scores));
    if (plane == NULL || block == NULL || scores == NULL) {
        SVT_ERROR("Failed to allocate memory for block of size %d\n", n);
        free(plane);
        free(block);
        free(scores);
        return -1;
    }

#ifdef NOISE_MODEL_LOG_SCORE
    SVT_ERROR("score = [");
#endif
    for (int32_t by = 0; by < num_blocks_h; ++by) {
        for (int32_t bx = 0; bx < num_blocks_w; ++bx) {
            // Compute gradient covariance matrix.
            double  g_xx = 0, g_xy = 0, g_yy = 0;
            double  var  = 0;
            double  mean = 0;
            svt_aom_flat_block_finder_extract_block(
                block_finder, data, w, h, stride, bx * block_size, by * block_size, plane, block);

            for (int32_t yi = 1; yi < block_size - 1; ++yi) {
                for (int32_t xi = 1; xi < block_size - 1; ++xi) {
                    const double gx =
                        (block[yi * block_size + xi + 1] - block[yi * block_size + xi - 1]) / 2;
                    const double gy = (block[yi * block_size + xi + block_size] -
                                       block[yi * block_size + xi - block_size]) /
                                      2;
                    g_xx += gx * gx;
                    g_xy += gx * gy;
                    g_yy += gy * gy;

                    mean += block[yi * block_size + xi];
                    var += block[yi * block_size + xi] * block[yi * block_size + xi];
                }
            }
            mean /= (block_size - 2) * (block_size - 2);

            // Normalize gradients by BlockSize.
            g_xx /= ((block_size - 2) * (block_size - 2));
            g_xy /= ((block_size - 2) * (block_size - 2));
            g_yy /= ((block_size - 2) * (block_size - 2));
            var = var / ((block_size - 2) * (block_size - 2)) - mean * mean;

            {
                const double  trace   = g_xx + g_yy;
                const double  det     = g_xx * g_yy - g_xy * g_xy;
                const double  e1      = (trace + sqrt(trace * trace - 4 * det)) / 2.;
                const double  e2      = (trace - sqrt(trace * trace - 4 * det)) / 2.;
                const double  norm    = e1; // Spectral norm
                const double  ratio   = (e1 / AOMMAX(e2, 1e-6));
                const int32_t is_flat = (trace < k_trace_threshold) &&
                                        (ratio < k_ratio_threshold) && (norm < k_norm_threshold) &&
                                        (var > k_var_threshold);
                // The following weights are used to combine the above features to give
                // a sigmoid score for flatness. If the input was normalized to [0,100]
                // the magnitude of these values would be close to 1 (e.g., weights
                // corresponding to variance would be a factor of 10000x smaller).
                // The weights are given in the following order:
                //    [{var}, {ratio}, {trace}, {norm}, offset]
                // with one of the most discriminative being simply the variance.
                const double weights[5] = {-6682, -0.2056, 13087, -12434, 2.5694};
                const float  score =
                    (float)(1.0 /
                            (1 + exp(-(weights[0] * var + weights[1] * ratio + weights[2] * trace +
                                       weights[3] * norm + weights[4]))));
                flat_blocks[by * num_blocks_w + bx]  = is_flat ? 255 : 0;
                scores[by * num_blocks_w + bx].score = var > k_var_threshold ? score : 0;
                scores[by * num_blocks_w + bx].index = by * num_blocks_w + bx;
#ifdef NOISE_MODEL_LOG_SCORE
                SVT_ERROR("%g %g %g %g %g %d ", score, var, ratio, trace, norm, is_flat);
#endif
                num_flat += is_flat;
            }
        }
#ifdef NOISE_MODEL_LOG_SCORE
        SVT_ERROR("\n");
#endif
    }
#ifdef NOISE_MODEL_LOG_SCORE
    SVT_ERROR("];\n");
#endif
    // Find the top-scored blocks (most likely to be flat) and set the flat blocks
    // be the union of the thresholded results and the top 10th percentile of the
    // scored results.
    qsort(scores, num_blocks_w * num_blocks_h, sizeof(*scores), &compare_scores);
    const int32_t top_nth_percentile = num_blocks_w * num_blocks_h * 90 / 100;
    const float   score_threshold    = scores[top_nth_percentile].score;
    for (int32_t i = 0; i < num_blocks_w * num_blocks_h; ++i) {
        if (scores[i].score >= score_threshold) {
            num_flat += flat_blocks[scores[i].index] == 0;
            flat_blocks[scores[i].index] |= 1;
        }
    }
    free(block);
    free(plane);
    free(scores);
    return num_flat;
}

int32_t svt_aom_noise_model_init(AomNoiseModel *model, const AomNoiseModelParams params) {
    const int32_t n         = num_coeffs(params);
    const int32_t lag       = params.lag;
    const int32_t bit_depth = params.bit_depth;
    int32_t       i         = 0;

    memset(model, 0, sizeof(*model));
    if (params.lag < 1) {
        SVT_ERROR("Invalid noise param: lag = %d must be >= 1\n", params.lag);
        return 0;
    }
    if (params.lag > k_max_lag) {
        SVT_ERROR("Invalid noise param: lag = %d must be <= %d\n", params.lag, k_max_lag);
        return 0;
    }
    if (svt_memcpy != NULL)
        svt_memcpy(&model->params, &params, sizeof(params));
    else
        svt_memcpy_c(&model->params, &params, sizeof(params));

    for (int c = 0; c < 3; ++c) {
        if (!noise_state_init(&model->combined_state[c], n + (c > 0), bit_depth)) {
            SVT_ERROR("Failed to allocate noise state for channel %d\n", c);
            svt_aom_noise_model_free(model);
            return 0;
        }
        if (!noise_state_init(&model->latest_state[c], n + (c > 0), bit_depth)) {
            SVT_ERROR("Failed to allocate noise state for channel %d\n", c);
            svt_aom_noise_model_free(model);
            return 0;
        }
    }
    model->n      = n;
    model->coords = (int32_t(*)[2])malloc(sizeof(*model->coords) * n);
    if (!model->coords) {
        SVT_ERROR("Failed to allocate memory for coords\n");
        svt_aom_noise_model_free(model);
        return 0;
    }
    for (int32_t y = -lag; y <= 0; ++y) {
        const int32_t max_x = y == 0 ? -1 : lag;
        for (int32_t x = -lag; x <= max_x; ++x) {
            switch (params.shape) {
            case AOM_NOISE_SHAPE_DIAMOND:
                if (abs(x) <= y + lag) {
                    model->coords[i][0] = x;
                    model->coords[i][1] = y;
                    ++i;
                }
                break;
            case AOM_NOISE_SHAPE_SQUARE:
                model->coords[i][0] = x;
                model->coords[i][1] = y;
                ++i;
                break;
            default:
                SVT_ERROR("Invalid shape\n");
                svt_aom_noise_model_free(model);
                return 0;
            }
        }
    }
    assert(i == n);
    return 1;
}

void svt_aom_noise_model_free(AomNoiseModel *model) {
    int32_t c = 0;
    if (!model) return;

    free(model->coords);
    for (c = 0; c < 3; ++c) {
        equation_system_free(&model->latest_state[c].eqns);
        equation_system_free(&model->combined_state[c].eqns);

        equation_system_free(&model->latest_state[c].strength_solver.eqns);
        equation_system_free(&model->combined_state[c].strength_solver.eqns);
    }
    memset(model, 0, sizeof(*model));
}

// Extracts the neighborhood defined by coords around point (x, y) from
// the difference between the data and denoised images. Also extracts the
// entry (possibly downsampled) for (x, y) in the alt_data (e.g., luma).
#define EXTRACT_AR_ROW(INT_TYPE, suffix)                                                 \
    static double extract_ar_row_##suffix(int32_t(*coords)[2],                           \
                                          int32_t               num_coords,              \
                                          const INT_TYPE *const data,                    \
                                          const INT_TYPE *const denoised,                \
                                          int32_t               stride,                  \
                                          int32_t               sub_log2[2],             \
                                          const INT_TYPE *const alt_data,                \
                                          const INT_TYPE *const alt_denoised,            \
                                          int32_t               alt_stride,              \
                                          int32_t               x,                       \
                                          int32_t               y,                       \
                                          double *              buffer) {                              \
        for (int32_t i = 0; i < num_coords; ++i) {                                       \
            const int32_t x_i = x + coords[i][0], y_i = y + coords[i][1];                \
            buffer[i] = (double)data[y_i * stride + x_i] - denoised[y_i * stride + x_i]; \
        }                                                                                \
        const double val = (double)data[y * stride + x] - denoised[y * stride + x];      \
                                                                                         \
        if (alt_data && alt_denoised) {                                                  \
            double  avg_data = 0, avg_denoised = 0;                                      \
            int32_t num_samples = 0;                                                     \
            for (int32_t dy_i = 0; dy_i < (1 << sub_log2[1]); dy_i++) {                  \
                const int32_t y_up = (y << sub_log2[1]) + dy_i;                          \
                for (int32_t dx_i = 0; dx_i < (1 << sub_log2[0]); dx_i++) {              \
                    const int32_t x_up = (x << sub_log2[0]) + dx_i;                      \
                    avg_data += alt_data[y_up * alt_stride + x_up];                      \
                    avg_denoised += alt_denoised[y_up * alt_stride + x_up];              \
                    num_samples++;                                                       \
                }                                                                        \
            }                                                                            \
            assert(num_samples > 0);                                                     \
            buffer[num_coords] = (avg_data - avg_denoised) / num_samples;                \
        }                                                                                \
        return val;                                                                      \
    }

EXTRACT_AR_ROW(uint8_t, lowbd);
EXTRACT_AR_ROW(uint16_t, highbd);

static int32_t add_block_observations(AomNoiseModel *noise_model, int32_t c,
                                      const uint8_t *const data, const uint8_t *const denoised,
                                      int32_t w, int32_t h, int32_t stride, int32_t sub_log2[2],
                                      const uint8_t *const alt_data,
                                      const uint8_t *const alt_denoised, int32_t alt_stride,
                                      const uint8_t *const flat_blocks, int32_t block_size,
                                      int32_t num_blocks_w, int32_t num_blocks_h) {
    const int32_t lag           = noise_model->params.lag;
    const int32_t num_coords    = noise_model->n;
    const double  normalization = (1 << noise_model->params.bit_depth) - 1;
    double *      A             = noise_model->latest_state[c].eqns.A;
    double *      b             = noise_model->latest_state[c].eqns.b;
    double *      buffer        = (double *)malloc(sizeof(*buffer) * (num_coords + 1));
    const int32_t n             = noise_model->latest_state[c].eqns.n;

    if (!buffer) {
        SVT_ERROR("Unable to allocate buffer of size %d\n", num_coords + 1);
        return 0;
    }
    for (int32_t by = 0; by < num_blocks_h; ++by) {
        const int32_t y_o = by * (block_size >> sub_log2[1]);
        for (int32_t bx = 0; bx < num_blocks_w; ++bx) {
            const int32_t x_o = bx * (block_size >> sub_log2[0]);
            if (!flat_blocks[by * num_blocks_w + bx]) continue;
            int32_t y_start = (by > 0 && flat_blocks[(by - 1) * num_blocks_w + bx]) ? 0 : lag;
            int32_t x_start = (bx > 0 && flat_blocks[by * num_blocks_w + bx - 1]) ? 0 : lag;
            int32_t y_end   = AOMMIN((h >> sub_log2[1]) - by * (block_size >> sub_log2[1]),
                                   block_size >> sub_log2[1]);
            int32_t x_end =
                AOMMIN((w >> sub_log2[0]) - bx * (block_size >> sub_log2[0]) - lag,
                       (bx + 1 < num_blocks_w && flat_blocks[by * num_blocks_w + bx + 1])
                           ? (block_size >> sub_log2[0])
                           : ((block_size >> sub_log2[0]) - lag));
            for (int32_t y = y_start; y < y_end; ++y) {
                for (int32_t x = x_start; x < x_end; ++x) {
                    const double val =
                        noise_model->params.use_highbd
                            ? extract_ar_row_highbd(noise_model->coords,
                                                    num_coords,
                                                    (const uint16_t *const)data,
                                                    (const uint16_t *const)denoised,
                                                    stride,
                                                    sub_log2,
                                                    (const uint16_t *const)alt_data,
                                                    (const uint16_t *const)alt_denoised,
                                                    alt_stride,
                                                    x + x_o,
                                                    y + y_o,
                                                    buffer)
                            : extract_ar_row_lowbd(noise_model->coords,
                                                   num_coords,
                                                   data,
                                                   denoised,
                                                   stride,
                                                   sub_log2,
                                                   alt_data,
                                                   alt_denoised,
                                                   alt_stride,
                                                   x + x_o,
                                                   y + y_o,
                                                   buffer);
                    for (int32_t i = 0; i < n; ++i) {
                        for (int32_t j = 0; j < n; ++j) {
                            A[i * n + j] +=
                                (buffer[i] * buffer[j]) / (normalization * normalization);
                        }
                        b[i] += (buffer[i] * val) / (normalization * normalization);
                    }
                    noise_model->latest_state[c].num_observations++;
                }
            }
        }
    }
    free(buffer);
    return 1;
}

static void add_noise_std_observations(AomNoiseModel *noise_model, int32_t c, const double *coeffs,
                                       const uint8_t *const data, const uint8_t *const denoised,
                                       int32_t w, int32_t h, int32_t stride, int32_t sub_log2[2],
                                       const uint8_t *const alt_data, int32_t alt_stride,
                                       const uint8_t *const flat_blocks, int32_t block_size,
                                       int32_t num_blocks_w, int32_t num_blocks_h) {
    const int32_t           num_coords            = noise_model->n;
    AomNoiseStrengthSolver *noise_strength_solver = &noise_model->latest_state[c].strength_solver;

    const AomNoiseStrengthSolver *noise_strength_luma =
        &noise_model->latest_state[0].strength_solver;
    const double luma_gain  = noise_model->latest_state[0].ar_gain;
    const double noise_gain = noise_model->latest_state[c].ar_gain;
    for (int32_t by = 0; by < num_blocks_h; ++by) {
        const int32_t y_o = by * (block_size >> sub_log2[1]);
        for (int32_t bx = 0; bx < num_blocks_w; ++bx) {
            const int32_t x_o = bx * (block_size >> sub_log2[0]);
            if (!flat_blocks[by * num_blocks_w + bx]) continue;
            const int32_t num_samples_h = AOMMIN(
                (h >> sub_log2[1]) - by * (block_size >> sub_log2[1]), block_size >> sub_log2[1]);
            const int32_t num_samples_w = AOMMIN(
                (w >> sub_log2[0]) - bx * (block_size >> sub_log2[0]), (block_size >> sub_log2[0]));
            // Make sure that we have a reasonable amount of samples to consider the
            // block
            if (num_samples_w * num_samples_h > block_size) {
                const double block_mean = get_block_mean(alt_data ? alt_data : data,
                                                         w,
                                                         h,
                                                         alt_data ? alt_stride : stride,
                                                         x_o << sub_log2[0],
                                                         y_o << sub_log2[1],
                                                         block_size,
                                                         noise_model->params.use_highbd);
                const double noise_var  = get_noise_var(data,
                                                       denoised,
                                                       stride,
                                                       w >> sub_log2[0],
                                                       h >> sub_log2[1],
                                                       x_o,
                                                       y_o,
                                                       block_size >> sub_log2[0],
                                                       block_size >> sub_log2[1],
                                                       noise_model->params.use_highbd);
                // We want to remove the part of the noise that came from being
                // correlated with luma. Note that the noise solver for luma must
                // have already been run.
                const double luma_strength = c > 0
                                                 ? luma_gain * noise_strength_solver_get_value(
                                                                   noise_strength_luma, block_mean)
                                                 : 0;
                const double corr = c > 0 ? coeffs[num_coords] : 0;
                // Chroma noise:
                //    N(0, noise_var) = N(0, uncorr_var) + corr * N(0, luma_strength^2)
                // The uncorrelated component:
                //   uncorr_var = noise_var - (corr * luma_strength)^2
                // But don't allow fully correlated noise (hence the max), since the
                // synthesis cannot model it.
                const double uncorr_std =
                    sqrt(AOMMAX(noise_var / 16, noise_var - pow(corr * luma_strength, 2)));
                // After we've removed correlation with luma, undo the gain that will
                // come from running the IIR filter.
                const double adjusted_strength = uncorr_std / noise_gain;
                svt_aom_noise_strength_solver_add_measurement(
                    noise_strength_solver, block_mean, adjusted_strength);
            }
        }
    }
}

static int32_t ar_equation_system_solve(AomNoiseState *state, int32_t is_chroma) {
    const int32_t ret = equation_system_solve(&state->eqns);
    state->ar_gain    = 1.0;
    if (!ret) return ret;

    // Update the AR gain from the equation system as it will be used to fit
    // the noise strength as a function of intensity.  In the Yule-Walker
    // equations, the diagonal should be the variance of the correlated noise.
    // In the case of the least squares estimate, there will be some variability
    // in the diagonal. So use the mean of the diagonal as the estimate of
    // overall variance (this works for least squares or Yule-Walker formulation).
    double        var = 0;
    const int32_t n   = state->eqns.n;
    for (int32_t i = 0; i < (state->eqns.n - is_chroma); ++i)
        var += state->eqns.A[i * n + i] / state->num_observations;
    var /= (n - is_chroma);

    // Keep track of E(Y^2) = <b, x> + E(X^2)
    // In the case that we are using chroma and have an estimate of correlation
    // with luma we adjust that estimate slightly to remove the correlated bits by
    // subtracting out the last column of a scaled by our correlation estimate
    // from b. E(y^2) = <b - A(:, end)*x(end), x>
    double sum_covar = 0;
    for (int32_t i = 0; i < state->eqns.n - is_chroma; ++i) {
        double bi = state->eqns.b[i];
        if (is_chroma) bi -= state->eqns.A[i * n + (n - 1)] * state->eqns.x[n - 1];
        sum_covar += (bi * state->eqns.x[i]) / state->num_observations;
    }
    // Now, get an estimate of the variance of uncorrelated noise signal and use
    // it to determine the gain of the AR filter.
    const double noise_var = AOMMAX(var - sum_covar, 1e-6);
    state->ar_gain         = AOMMAX(1, sqrt(AOMMAX(var / noise_var, 1e-6)));
    return ret;
}

AomNoiseStatus svt_aom_noise_model_update(AomNoiseModel *const noise_model,
                                         const uint8_t *const data[3],
                                         const uint8_t *const denoised[3], int32_t w, int32_t h,
                                         int32_t stride[3], int32_t chroma_sub_log2[2],
                                         const uint8_t *const flat_blocks, int32_t block_size) {
    const int32_t num_blocks_w = (w + block_size - 1) / block_size;
    const int32_t num_blocks_h = (h + block_size - 1) / block_size;
    //  int32_t y_model_different = 0;
    int32_t num_blocks = 0;
    int32_t i = 0, channel = 0;

    if (block_size <= 1) {
        SVT_ERROR("BlockSize = %d must be > 1\n", block_size);
        return AOM_NOISE_STATUS_INVALID_ARGUMENT;
    }

    if (block_size < noise_model->params.lag * 2 + 1) {
        SVT_ERROR("BlockSize = %d must be >= %d\n", block_size, noise_model->params.lag * 2 + 1);
        return AOM_NOISE_STATUS_INVALID_ARGUMENT;
    }

    // Clear the latest equation system
    for (i = 0; i < 3; ++i) {
        equation_system_clear(&noise_model->latest_state[i].eqns);
        noise_model->latest_state[i].num_observations = 0;
        noise_strength_solver_clear(&noise_model->latest_state[i].strength_solver);
    }

    // Check that we have enough flat blocks
    for (i = 0; i < num_blocks_h * num_blocks_w; ++i) {
        if (flat_blocks[i]) num_blocks++;
    }

    if (num_blocks <= 1) {
        SVT_ERROR("Not enough flat blocks to update noise estimate\n");
        return AOM_NOISE_STATUS_INSUFFICIENT_FLAT_BLOCKS;
    }

    for (channel = 0; channel < 3; ++channel) {
        int32_t        no_subsampling[2] = {0, 0};
        const uint8_t *alt_data          = channel > 0 ? data[0] : 0;
        const uint8_t *alt_denoised      = channel > 0 ? denoised[0] : 0;
        int32_t *      sub               = channel > 0 ? chroma_sub_log2 : no_subsampling;
        const int32_t  is_chroma         = channel != 0;
        if (!data[channel] || !denoised[channel]) break;
        if (!add_block_observations(noise_model,
                                    channel,
                                    data[channel],
                                    denoised[channel],
                                    w,
                                    h,
                                    stride[channel],
                                    sub,
                                    alt_data,
                                    alt_denoised,
                                    stride[0],
                                    flat_blocks,
                                    block_size,
                                    num_blocks_w,
                                    num_blocks_h)) {
            SVT_ERROR("Adding block observation failed\n");
            return AOM_NOISE_STATUS_INTERNAL_ERROR;
        }

        if (!ar_equation_system_solve(&noise_model->latest_state[channel], is_chroma)) {
            if (is_chroma) {
                set_chroma_coefficient_fallback_soln(&noise_model->latest_state[channel].eqns);
            } else {
                SVT_ERROR("Solving latest noise equation system failed %d!\n", channel);
                return AOM_NOISE_STATUS_INTERNAL_ERROR;
            }
        }

        add_noise_std_observations(noise_model,
                                   channel,
                                   noise_model->latest_state[channel].eqns.x,
                                   data[channel],
                                   denoised[channel],
                                   w,
                                   h,
                                   stride[channel],
                                   sub,
                                   alt_data,
                                   stride[0],
                                   flat_blocks,
                                   block_size,
                                   num_blocks_w,
                                   num_blocks_h);

        if (!svt_aom_noise_strength_solver_solve(
                &noise_model->latest_state[channel].strength_solver)) {
            SVT_ERROR("Solving latest noise strength failed!\n");
            return AOM_NOISE_STATUS_INTERNAL_ERROR;
        }

        // Check noise characteristics and return if error.
        //    if (channel == 0 &&
        //        noise_model->combined_state[channel].strength_solver.num_equations >
        //            0 &&
        //        is_noise_model_different(noise_model)) {
        //      y_model_different = 1;
        //    }

        noise_model->combined_state[channel].num_observations =
            noise_model->latest_state[channel].num_observations;
        equation_system_copy(&noise_model->combined_state[channel].eqns,
                             &noise_model->latest_state[channel].eqns);
        if (!ar_equation_system_solve(&noise_model->combined_state[channel], is_chroma)) {
            if (is_chroma) {
                set_chroma_coefficient_fallback_soln(&noise_model->combined_state[channel].eqns);
            } else {
                SVT_ERROR("Solving combined noise equation system failed %d!\n", channel);
                return AOM_NOISE_STATUS_INTERNAL_ERROR;
            }
        }

        noise_strength_solver_copy(&noise_model->combined_state[channel].strength_solver,
                                   &noise_model->latest_state[channel].strength_solver);

        if (!svt_aom_noise_strength_solver_solve(
                &noise_model->combined_state[channel].strength_solver)) {
            SVT_ERROR("Solving combined noise strength failed!\n");
            return AOM_NOISE_STATUS_INTERNAL_ERROR;
        }
    }

    return AOM_NOISE_STATUS_OK;
}

void svt_aom_noise_model_save_latest(AomNoiseModel *noise_model) {
    for (int32_t c = 0; c < 3; c++) {
        equation_system_copy(&noise_model->combined_state[c].eqns,
                             &noise_model->latest_state[c].eqns);
        equation_system_copy(&noise_model->combined_state[c].strength_solver.eqns,
                             &noise_model->latest_state[c].strength_solver.eqns);
        noise_model->combined_state[c].strength_solver.num_equations =
            noise_model->latest_state[c].strength_solver.num_equations;
        noise_model->combined_state[c].num_observations =
            noise_model->latest_state[c].num_observations;
        noise_model->combined_state[c].ar_gain = noise_model->latest_state[c].ar_gain;
    }
}

int32_t svt_aom_noise_model_get_grain_parameters(AomNoiseModel *const noise_model,
                                                AomFilmGrain *       film_grain) {
    if (noise_model->params.lag > 3) {
        SVT_ERROR("params.lag = %d > 3\n", noise_model->params.lag);
        return 0;
    }
    uint16_t random_seed = film_grain->random_seed;
    memset(film_grain, 0, sizeof(*film_grain));
    film_grain->random_seed = random_seed;

    film_grain->apply_grain       = 1;
    film_grain->update_parameters = 1;

    film_grain->ar_coeff_lag = noise_model->params.lag;

    // Convert the scaling functions to 8 bit values
    AomNoiseStrengthLut scaling_points[3];
    svt_aom_noise_strength_solver_fit_piecewise(
        &noise_model->combined_state[0].strength_solver, 14, scaling_points + 0);
    svt_aom_noise_strength_solver_fit_piecewise(
        &noise_model->combined_state[1].strength_solver, 10, scaling_points + 1);
    svt_aom_noise_strength_solver_fit_piecewise(
        &noise_model->combined_state[2].strength_solver, 10, scaling_points + 2);

    // Both the domain and the range of the scaling functions in the film_grain
    // are normalized to 8-bit (e.g., they are implicitly scaled during grain
    // synthesis).
    const double strength_divisor  = 1 << (noise_model->params.bit_depth - 8);
    double       max_scaling_value = 1e-4;
    for (int32_t c = 0; c < 3; ++c) {
        for (int32_t i = 0; i < scaling_points[c].num_points; ++i) {
            scaling_points[c].points[i][0] =
                AOMMIN(255, scaling_points[c].points[i][0] / strength_divisor);
            scaling_points[c].points[i][1] =
                AOMMIN(255, scaling_points[c].points[i][1] / strength_divisor);
            max_scaling_value = AOMMAX(scaling_points[c].points[i][1], max_scaling_value);
        }
    }

    // Scaling_shift values are in the range [8,11]
    const int32_t max_scaling_value_log2 = clamp((int32_t)floor(log2(max_scaling_value) + 1), 2, 5);
    film_grain->scaling_shift            = 5 + (8 - max_scaling_value_log2);

    const double scale_factor = 1 << (8 - max_scaling_value_log2);
    film_grain->num_y_points  = scaling_points[0].num_points;
    film_grain->num_cb_points = scaling_points[1].num_points;
    film_grain->num_cr_points = scaling_points[2].num_points;

    int32_t(*film_grain_scaling[3])[2] = {
        film_grain->scaling_points_y,
        film_grain->scaling_points_cb,
        film_grain->scaling_points_cr,
    };
    for (int32_t c = 0; c < 3; c++) {
        for (int32_t i = 0; i < scaling_points[c].num_points; ++i) {
            film_grain_scaling[c][i][0] = (int32_t)(scaling_points[c].points[i][0] + 0.5);
            film_grain_scaling[c][i][1] =
                clamp((int32_t)(scale_factor * scaling_points[c].points[i][1] + 0.5), 0, 255);
        }
    }
    svt_aom_noise_strength_lut_free(scaling_points + 0);
    svt_aom_noise_strength_lut_free(scaling_points + 1);
    svt_aom_noise_strength_lut_free(scaling_points + 2);

    // Convert the ar_coeffs into 8-bit values
    const int32_t n_coeff   = noise_model->combined_state[0].eqns.n;
    double        max_coeff = 1e-4, min_coeff = -1e-4;
    double        y_corr[2]         = {0, 0};
    double        avg_luma_strength = 0;
    for (int32_t c = 0; c < 3; c++) {
        AomEquationSystem *eqns = &noise_model->combined_state[c].eqns;
        for (int32_t i = 0; i < n_coeff; ++i) {
            max_coeff = AOMMAX(max_coeff, eqns->x[i]);
            min_coeff = AOMMIN(min_coeff, eqns->x[i]);
        }
        // Since the correlation between luma/chroma was computed in an already
        // scaled space, we adjust it in the un-scaled space.
        AomNoiseStrengthSolver *solver = &noise_model->combined_state[c].strength_solver;
        // Compute a weighted average of the strength for the channel.
        double average_strength = 0, total_weight = 0;
        for (int32_t i = 0; i < solver->eqns.n; ++i) {
            double w = 0;
            for (int32_t j = 0; j < solver->eqns.n; ++j)
                w += solver->eqns.A[i * solver->eqns.n + j];
            w = sqrt(w);
            average_strength += solver->eqns.x[i] * w;
            total_weight += w;
        }
        if (total_weight == 0)
            average_strength = 1;
        else
            average_strength /= total_weight;
        if (c == 0)
            avg_luma_strength = average_strength;
        else {
            y_corr[c - 1] = avg_luma_strength * eqns->x[n_coeff] / average_strength;
            max_coeff     = AOMMAX(max_coeff, y_corr[c - 1]);
            min_coeff     = AOMMIN(min_coeff, y_corr[c - 1]);
        }
    }
    // Shift value: AR coeffs range (values 6-9)
    // 6: [-2, 2),  7: [-1, 1), 8: [-0.5, 0.5), 9: [-0.25, 0.25)
    film_grain->ar_coeff_shift =
        clamp(7 - (int32_t)AOMMAX(1 + floor(log2(max_coeff)), ceil(log2(-min_coeff))), 6, 9);
    double   scale_ar_coeff = 1 << film_grain->ar_coeff_shift;
    int32_t *ar_coeffs[3]   = {
        film_grain->ar_coeffs_y,
        film_grain->ar_coeffs_cb,
        film_grain->ar_coeffs_cr,
    };
    for (int32_t c = 0; c < 3; ++c) {
        AomEquationSystem *eqns = &noise_model->combined_state[c].eqns;
        for (int32_t i = 0; i < n_coeff; ++i) {
            ar_coeffs[c][i] = clamp((int32_t)round(scale_ar_coeff * eqns->x[i]), -128, 127);
        }
        if (c > 0) {
            ar_coeffs[c][n_coeff] =
                clamp((int32_t)round(scale_ar_coeff * y_corr[c - 1]), -128, 127);
        }
    }

    // At the moment, the noise modeling code assumes that the chroma scaling
    // functions are a function of luma.
    film_grain->cb_mult      = 128; // 8 bits
    film_grain->cb_luma_mult = 192; // 8 bits
    film_grain->cb_offset    = 256; // 9 bits

    film_grain->cr_mult      = 128; // 8 bits
    film_grain->cr_luma_mult = 192; // 8 bits
    film_grain->cr_offset    = 256; // 9 bits

    film_grain->chroma_scaling_from_luma = 0;
    film_grain->grain_scale_shift        = 0;
    film_grain->overlap_flag             = 1;
    return 1;
}

static void pointwise_multiply(const float *a, float *b, int32_t n) {
    for (int32_t i = 0; i < n; ++i) b[i] *= a[i];
}

static float *get_half_cos_window(int32_t block_size) {
    float *window_function = (float *)malloc(block_size * block_size * sizeof(*window_function));
    ASSERT(window_function);
    for (int32_t y = 0; y < block_size; ++y) {
        const double cos_yd = cos((.5 + y) * PI / block_size - PI / 2);
        for (int32_t x = 0; x < block_size; ++x) {
            const double cos_xd                 = cos((.5 + x) * PI / block_size - PI / 2);
            window_function[y * block_size + x] = (float)(cos_yd * cos_xd);
        }
    }
    return window_function;
}

#define DITHER_AND_QUANTIZE(INT_TYPE, suffix)                                                      \
    static void dither_and_quantize_##suffix(float *   result,                                     \
                                             int32_t   result_stride,                              \
                                             INT_TYPE *denoised,                                   \
                                             int32_t   w,                                          \
                                             int32_t   h,                                          \
                                             int32_t   stride,                                     \
                                             int32_t   chroma_sub_w,                               \
                                             int32_t   chroma_sub_h,                               \
                                             int32_t   block_size,                                 \
                                             float     block_normalization) {                          \
        for (int32_t y = 0; y < (h >> chroma_sub_h); ++y) {                                        \
            for (int32_t x = 0; x < (w >> chroma_sub_w); ++x) {                                    \
                const int32_t result_idx = (y + (block_size >> chroma_sub_h)) * result_stride +    \
                                           x + (block_size >> chroma_sub_w);                       \
                INT_TYPE new_val =                                                                 \
                    (INT_TYPE)AOMMIN(AOMMAX(result[result_idx] * block_normalization + 0.5f, 0),   \
                                     block_normalization);                                         \
                const float err = -(((float)new_val) / block_normalization - result[result_idx]);  \
                denoised[y * stride + x] = new_val;                                                \
                if (x + 1 < (w >> chroma_sub_w)) { result[result_idx + 1] += err * 7.0f / 16.0f; } \
                if (y + 1 < (h >> chroma_sub_h)) {                                                 \
                    if (x > 0) { result[result_idx + result_stride - 1] += err * 3.0f / 16.0f; }   \
                    result[result_idx + result_stride] += err * 5.0f / 16.0f;                      \
                    if (x + 1 < (w >> chroma_sub_w)) {                                             \
                        result[result_idx + result_stride + 1] += err * 1.0f / 16.0f;              \
                    }                                                                              \
                }                                                                                  \
            }                                                                                      \
        }                                                                                          \
    }

DITHER_AND_QUANTIZE(uint8_t, lowbd);
DITHER_AND_QUANTIZE(uint16_t, highbd);

int32_t svt_aom_wiener_denoise_2d(const uint8_t *const data[3], uint8_t *denoised[3], int32_t w,
                                 int32_t h, int32_t stride[3], int32_t chroma_sub[2],
                                 float *noise_psd[3], int32_t block_size, int32_t bit_depth,
                                 int32_t use_highbd) {
    float *plane = NULL, *window_full = NULL, *window_chroma = NULL;
    DECLARE_ALIGNED(32, float, *block);
    block                          = NULL;
    double *               block_d = NULL, *plane_d = NULL;
    struct aom_noise_tx_t *tx_full       = NULL;
    struct aom_noise_tx_t *tx_chroma     = NULL;
    const int32_t          num_blocks_w  = (w + block_size - 1) / block_size;
    const int32_t          num_blocks_h  = (h + block_size - 1) / block_size;
    const int32_t          result_stride = (num_blocks_w + 2) * block_size;
    const int32_t          result_height = (num_blocks_h + 2) * block_size;
    float *                result        = NULL;
    int32_t                init_success  = 1;
    AomFlatBlockFinder     block_finder_full;
    AomFlatBlockFinder     block_finder_chroma;
    const float            k_block_normalization = (float)((1 << bit_depth) - 1);
    if (chroma_sub[0] != chroma_sub[1]) {
        SVT_ERROR(
            "svt_aom_wiener_denoise_2d doesn't handle different chroma "
            "subsampling");
        return 0;
    }
    init_success &=
        svt_aom_flat_block_finder_init(&block_finder_full, block_size, bit_depth, use_highbd);
    result  = (float *)malloc((num_blocks_h + 2) * block_size * result_stride * sizeof(*result));
    plane   = (float *)malloc(block_size * block_size * sizeof(*plane));
    block   = (float *)svt_aom_memalign(32, 2 * block_size * block_size * sizeof(*block));
    block_d = (double *)malloc(block_size * block_size * sizeof(*block_d));
    plane_d = (double *)malloc(block_size * block_size * sizeof(*plane_d));
    window_full = get_half_cos_window(block_size);
    tx_full     = svt_aom_noise_tx_malloc(block_size);

    if (chroma_sub[0] != 0) {
        init_success &= svt_aom_flat_block_finder_init(
            &block_finder_chroma, block_size >> chroma_sub[0], bit_depth, use_highbd);
        window_chroma = get_half_cos_window(block_size >> chroma_sub[0]);
        tx_chroma     = svt_aom_noise_tx_malloc(block_size >> chroma_sub[0]);
    } else {
        window_chroma = window_full;
        tx_chroma     = tx_full;
    }

    init_success &= (int32_t)((tx_full != NULL) && (tx_chroma != NULL) && (plane != NULL) &&
                              (plane_d != NULL) && (block != NULL) && (block_d != NULL) &&
                              (window_full != NULL) && (window_chroma != NULL) && (result != NULL));
    for (int32_t c = init_success ? 0 : 3; c < 3; ++c) {
        float *                window_function = c == 0 ? window_full : window_chroma;
        AomFlatBlockFinder *   block_finder    = &block_finder_full;
        const int32_t          chroma_sub_h    = c > 0 ? chroma_sub[1] : 0;
        const int32_t          chroma_sub_w    = c > 0 ? chroma_sub[0] : 0;
        struct aom_noise_tx_t *tx              = (c > 0 && chroma_sub[0] > 0) ? tx_chroma : tx_full;
        if (!data[c] || !denoised[c]) continue;
        if (c > 0 && chroma_sub[0] != 0) block_finder = &block_finder_chroma;
        memset(result, 0, sizeof(*result) * result_stride * result_height);
        // Do overlapped block processing (half overlapped). The block rows can
        // easily be done in parallel
        for (int32_t offsy = 0; offsy < (block_size >> chroma_sub_h);
             offsy += (block_size >> chroma_sub_h) / 2) {
            for (int32_t offsx = 0; offsx < (block_size >> chroma_sub_w);
                 offsx += (block_size >> chroma_sub_w) / 2) {
                // Pad the boundary when processing each block-set.
                for (int32_t by = -1; by < num_blocks_h; ++by) {
                    for (int32_t bx = -1; bx < num_blocks_w; ++bx) {
                        const int32_t pixels_per_block =
                            (block_size >> chroma_sub_w) * (block_size >> chroma_sub_h);
                        svt_aom_flat_block_finder_extract_block(
                            block_finder,
                            data[c],
                            w >> chroma_sub_w,
                            h >> chroma_sub_h,
                            stride[c],
                            bx * (block_size >> chroma_sub_w) + offsx,
                            by * (block_size >> chroma_sub_h) + offsy,
                            plane_d,
                            block_d);
                        for (int32_t j = 0; j < pixels_per_block; ++j) {
                            block[j] = (float)block_d[j];
                            plane[j] = (float)plane_d[j];
                        }
                        pointwise_multiply(window_function, block, pixels_per_block);
                        svt_aom_noise_tx_forward(tx, block);
                        svt_aom_noise_tx_filter(tx, noise_psd[c]);
                        svt_aom_noise_tx_inverse(tx, block);

                        // Apply window function to the plane approximation (we will apply
                        // it to the sum of plane + block when composing the results).
                        pointwise_multiply(window_function, plane, pixels_per_block);

                        for (int32_t y = 0; y < (block_size >> chroma_sub_h); ++y) {
                            const int32_t y_result =
                                y + (by + 1) * (block_size >> chroma_sub_h) + offsy;
                            for (int32_t x = 0; x < (block_size >> chroma_sub_w); ++x) {
                                const int32_t x_result =
                                    x + (bx + 1) * (block_size >> chroma_sub_w) + offsx;
                                result[y_result * result_stride + x_result] +=
                                    (block[y * (block_size >> chroma_sub_w) + x] +
                                     plane[y * (block_size >> chroma_sub_w) + x]) *
                                    window_function[y * (block_size >> chroma_sub_w) + x];
                            }
                        }
                    }
                }
            }
        }
        if (use_highbd) {
            dither_and_quantize_highbd(result,
                                       result_stride,
                                       (uint16_t *)denoised[c],
                                       w,
                                       h,
                                       stride[c],
                                       chroma_sub_w,
                                       chroma_sub_h,
                                       block_size,
                                       k_block_normalization);
        } else {
            dither_and_quantize_lowbd(result,
                                      result_stride,
                                      denoised[c],
                                      w,
                                      h,
                                      stride[c],
                                      chroma_sub_w,
                                      chroma_sub_h,
                                      block_size,
                                      k_block_normalization);
        }
    }
    free(result);
    free(plane);
    svt_aom_free(block);
    free(plane_d);
    free(block_d);
    free(window_full);

    svt_aom_noise_tx_free(tx_full);

    svt_aom_flat_block_finder_free(&block_finder_full);
    if (chroma_sub[0] != 0) {
        svt_aom_flat_block_finder_free(&block_finder_chroma);
        free(window_chroma);
        svt_aom_noise_tx_free(tx_chroma);
    }
    return init_success;
}

EbErrorType svt_aom_denoise_and_model_alloc(AomDenoiseAndModel *ctx, int32_t bit_depth,
                                           int32_t block_size, float noise_level) {
    ctx->block_size  = block_size;
    ctx->noise_level = noise_level;
    ctx->bit_depth   = bit_depth;

    EB_MALLOC_ARRAY(ctx->noise_psd[0], block_size * block_size);
    EB_MALLOC_ARRAY(ctx->noise_psd[1], block_size * block_size);
    EB_MALLOC_ARRAY(ctx->noise_psd[2], block_size * block_size);
    return EB_ErrorNone;
}

static void denoise_and_model_dctor(EbPtr p) {
    AomDenoiseAndModel *obj = (AomDenoiseAndModel *)p;

    free(obj->flat_blocks);
    for (int32_t i = 0; i < 3; ++i) {
        EB_FREE_ARRAY(obj->denoised[i]);
        EB_FREE_ARRAY(obj->noise_psd[i]);
        EB_FREE_ARRAY(obj->packed[i]);
    }
    svt_aom_noise_model_free(&obj->noise_model);
    svt_aom_flat_block_finder_free(&obj->flat_block_finder);
}

EbErrorType denoise_and_model_ctor(AomDenoiseAndModel *object_ptr, EbPtr object_init_data_ptr) {
    DenoiseAndModelInitData *init_data_ptr = (DenoiseAndModelInitData *)object_init_data_ptr;
    EbErrorType              return_error  = EB_ErrorNone;
    uint32_t                 use_highbd    = init_data_ptr->encoder_bit_depth > EB_8BIT ? 10 : 8;

    int32_t chroma_sub_log2[2] = {1, 1}; //todo: send chroma subsampling
    chroma_sub_log2[0]         = (init_data_ptr->encoder_color_format == EB_YUV444 ? 1 : 2) - 1;
    chroma_sub_log2[1]         = (init_data_ptr->encoder_color_format >= EB_YUV422 ? 1 : 2) - 1;

    object_ptr->dctor = denoise_and_model_dctor;

    return_error =
        svt_aom_denoise_and_model_alloc(object_ptr,
                                       init_data_ptr->encoder_bit_depth > EB_8BIT ? 10 : 8,
                                       DENOISING_BlockSize,
                                       (float)(init_data_ptr->noise_level / 10.0));
    if (return_error != EB_ErrorNone) return return_error;
    object_ptr->width     = init_data_ptr->width;
    object_ptr->height    = init_data_ptr->height;
    object_ptr->y_stride  = init_data_ptr->stride_y;
    object_ptr->uv_stride = init_data_ptr->stride_cb;

    //todo: consider replacing with EbPictureBuffersDesc

    EB_MALLOC_ARRAY(object_ptr->denoised[0],
                    (object_ptr->y_stride * object_ptr->height) << use_highbd);
    EB_MALLOC_ARRAY(object_ptr->denoised[1],
                    (object_ptr->uv_stride * (object_ptr->height >> chroma_sub_log2[0]))
                        << use_highbd);
    EB_MALLOC_ARRAY(object_ptr->denoised[2],
                    (object_ptr->uv_stride * (object_ptr->height >> chroma_sub_log2[0]))
                        << use_highbd);

    if (use_highbd) {
        EB_MALLOC_ARRAY(object_ptr->packed[0], (object_ptr->y_stride * object_ptr->height));
        EB_MALLOC_ARRAY(object_ptr->packed[1],
                        (object_ptr->uv_stride * (object_ptr->height >> chroma_sub_log2[0])));
        EB_MALLOC_ARRAY(object_ptr->packed[2],
                        (object_ptr->uv_stride * (object_ptr->height >> chroma_sub_log2[0])));
    }

    return return_error;
}

static int32_t denoise_and_model_realloc_if_necessary(struct AomDenoiseAndModel *ctx,
                                                      EbPictureBufferDesc *sd, int32_t use_highbd) {
    int32_t chroma_sub_log2[2] = {1, 1}; //todo: send chroma subsampling

    const int32_t block_size = ctx->block_size;

    free(ctx->flat_blocks);
    ctx->flat_blocks = NULL;

    ctx->num_blocks_w = (sd->width + ctx->block_size - 1) / ctx->block_size;
    ctx->num_blocks_h = (sd->height + ctx->block_size - 1) / ctx->block_size;
    ctx->flat_blocks  = malloc(ctx->num_blocks_w * ctx->num_blocks_h);

    svt_aom_flat_block_finder_free(&ctx->flat_block_finder);
    if (!svt_aom_flat_block_finder_init(
            &ctx->flat_block_finder, ctx->block_size, ctx->bit_depth, use_highbd)) {
        SVT_ERROR("Unable to init flat block finder\n");
        return 0;
    }

    const AomNoiseModelParams params = {AOM_NOISE_SHAPE_SQUARE, 3, ctx->bit_depth, use_highbd};
    //  svt_aom_noise_model_free(&ctx->noise_model);
    if (!svt_aom_noise_model_init(&ctx->noise_model, params)) {
        SVT_ERROR("Unable to init noise model\n");
        return 0;
    }

    // Simply use a flat PSD (although we could use the flat blocks to estimate
    // PSD) those to estimate an actual noise PSD)
    const float y_noise_level =
        svt_aom_noise_psd_get_default_value(ctx->block_size, ctx->noise_level);
    const float uv_noise_level =
        svt_aom_noise_psd_get_default_value(ctx->block_size >> chroma_sub_log2[1], ctx->noise_level);
    for (int32_t i = 0; i < block_size * block_size; ++i) {
        ctx->noise_psd[0][i] = y_noise_level;
        ctx->noise_psd[1][i] = ctx->noise_psd[2][i] = uv_noise_level;
    }
    return 1;
}

static void pack_2d_pic(EbPictureBufferDesc *input_picture, uint16_t *packed[3]) {
    const uint32_t input_luma_offset =
        ((input_picture->origin_y) * input_picture->stride_y) + (input_picture->origin_x);
    const uint32_t input_bit_inc_luma_offset =
        ((input_picture->origin_y) * input_picture->stride_bit_inc_y) + (input_picture->origin_x);
    const uint32_t input_cb_offset = (((input_picture->origin_y) >> 1) * input_picture->stride_cb) +
                                     ((input_picture->origin_x) >> 1);
    const uint32_t input_bit_inc_cb_offset =
        (((input_picture->origin_y) >> 1) * input_picture->stride_bit_inc_cb) +
        ((input_picture->origin_x) >> 1);
    const uint32_t input_cr_offset = (((input_picture->origin_y) >> 1) * input_picture->stride_cr) +
                                     ((input_picture->origin_x) >> 1);
    const uint32_t input_bit_inc_cr_offset =
        (((input_picture->origin_y) >> 1) * input_picture->stride_bit_inc_cr) +
        ((input_picture->origin_x) >> 1);

    pack2d_src(input_picture->buffer_y + input_luma_offset,
               input_picture->stride_y,
               input_picture->buffer_bit_inc_y + input_bit_inc_luma_offset,
               input_picture->stride_bit_inc_y,
               (uint16_t *)packed[0],
               input_picture->stride_y,
               input_picture->width,
               input_picture->height);

    pack2d_src(input_picture->buffer_cb + input_cb_offset,
               input_picture->stride_cr,
               input_picture->buffer_bit_inc_cb + input_bit_inc_cb_offset,
               input_picture->stride_bit_inc_cr,
               (uint16_t *)packed[1],
               input_picture->stride_cr,
               input_picture->width >> 1,
               input_picture->height >> 1);

    pack2d_src(input_picture->buffer_cr + input_cr_offset,
               input_picture->stride_cr,
               input_picture->buffer_bit_inc_cr + input_bit_inc_cr_offset,
               input_picture->stride_bit_inc_cr,
               (uint16_t *)packed[2],
               input_picture->stride_cr,
               input_picture->width >> 1,
               input_picture->height >> 1);
}

static void unpack_2d_pic(uint8_t *packed[3], EbPictureBufferDesc *outputPicturePtr) {
    uint32_t luma_buffer_offset =
        ((outputPicturePtr->origin_y) * outputPicturePtr->stride_y) + (outputPicturePtr->origin_x);
    uint32_t chroma_buffer_offset =
        (((outputPicturePtr->origin_y) >> 1) * outputPicturePtr->stride_cb) +
        ((outputPicturePtr->origin_x) >> 1);
    uint16_t luma_width    = (uint16_t)(outputPicturePtr->width);
    uint16_t chroma_width  = luma_width >> 1;
    uint16_t luma_height   = (uint16_t)(outputPicturePtr->height);
    uint16_t chroma_height = luma_height >> 1;

    un_pack2d((uint16_t *)(packed[0]),
              outputPicturePtr->stride_y,
              outputPicturePtr->buffer_y + luma_buffer_offset,
              outputPicturePtr->stride_y,
              outputPicturePtr->buffer_bit_inc_y + luma_buffer_offset,
              outputPicturePtr->stride_bit_inc_y,
              luma_width,
              luma_height);

    un_pack2d((uint16_t *)(packed[1]),
              outputPicturePtr->stride_cb,
              outputPicturePtr->buffer_cb + chroma_buffer_offset,
              outputPicturePtr->stride_cb,
              outputPicturePtr->buffer_bit_inc_cb + chroma_buffer_offset,
              outputPicturePtr->stride_bit_inc_cb,
              chroma_width,
              chroma_height);

    un_pack2d((uint16_t *)(packed[2]),
              outputPicturePtr->stride_cr,
              outputPicturePtr->buffer_cr + chroma_buffer_offset,
              outputPicturePtr->stride_cr,
              outputPicturePtr->buffer_bit_inc_cr + chroma_buffer_offset,
              outputPicturePtr->stride_bit_inc_cr,
              chroma_width,
              chroma_height);
}

int32_t svt_aom_denoise_and_model_run(struct AomDenoiseAndModel *ctx, EbPictureBufferDesc *sd,
                                     AomFilmGrain *film_grain, int32_t use_highbd) {
    const int32_t block_size = ctx->block_size;
    uint8_t *     raw_data[3];
    int32_t       chroma_sub_log2[2] = {1, 1}; //todo: send chroma subsampling
    int32_t       strides[3]         = {sd->stride_y, sd->stride_cb, sd->stride_cr};

    if (!denoise_and_model_realloc_if_necessary(ctx, sd, use_highbd)) {
        SVT_ERROR("Unable to realloc buffers\n");
        return 0;
    }

    if (!use_highbd) { // 8 bits input
        raw_data[0] = sd->buffer_y + sd->origin_y * sd->stride_y + sd->origin_x;
        raw_data[1] = sd->buffer_cb + sd->stride_cb * (sd->origin_y >> chroma_sub_log2[0]) +
                      (sd->origin_x >> chroma_sub_log2[1]);
        raw_data[2] = sd->buffer_cr + sd->stride_cr * (sd->origin_y >> chroma_sub_log2[0]) +
                      (sd->origin_x >> chroma_sub_log2[1]);
    } else { // 10 bits input
        pack_2d_pic(sd, ctx->packed);

        raw_data[0] = (uint8_t *)(ctx->packed[0]);
        raw_data[1] = (uint8_t *)(ctx->packed[1]);
        raw_data[2] = (uint8_t *)(ctx->packed[2]);
    }

    const uint8_t *const data[3] = {raw_data[0], raw_data[1], raw_data[2]};

    svt_aom_flat_block_finder_run(
        &ctx->flat_block_finder, data[0], sd->width, sd->height, strides[0], ctx->flat_blocks);

    if (!svt_aom_wiener_denoise_2d(data,
                                  ctx->denoised,
                                  sd->width,
                                  sd->height,
                                  strides,
                                  chroma_sub_log2,
                                  ctx->noise_psd,
                                  block_size,
                                  ctx->bit_depth,
                                  use_highbd)) {
        SVT_ERROR("Unable to denoise image\n");
        return 0;
    }

    const AomNoiseStatus status = svt_aom_noise_model_update(&ctx->noise_model,
                                                            data,
                                                            (const uint8_t *const *)ctx->denoised,
                                                            sd->width,
                                                            sd->height,
                                                            strides,
                                                            chroma_sub_log2,
                                                            ctx->flat_blocks,
                                                            block_size);

    int32_t have_noise_estimate = 0;
    if (status == AOM_NOISE_STATUS_OK || status == AOM_NOISE_STATUS_DIFFERENT_NOISE_TYPE) {
        svt_aom_noise_model_save_latest(&ctx->noise_model);
        have_noise_estimate = 1;
    }

    film_grain->apply_grain = 0;
    if (have_noise_estimate) {
        if (!svt_aom_noise_model_get_grain_parameters(&ctx->noise_model, film_grain)) {
            SVT_ERROR("Unable to get grain parameters.\n");
            return 0;
        }
        film_grain->apply_grain = 1;

        if (!use_highbd) {
            if (svt_memcpy != NULL)
            {
                svt_memcpy(raw_data[0], ctx->denoised[0], (strides[0] * sd->height) << use_highbd);
                svt_memcpy(raw_data[1],
                    ctx->denoised[1],
                    (strides[1] * (sd->height >> chroma_sub_log2[0])) << use_highbd);
                svt_memcpy(raw_data[2],
                    ctx->denoised[2],
                    (strides[2] * (sd->height >> chroma_sub_log2[0])) << use_highbd);
            }
            else
            {
                svt_memcpy_c(raw_data[0], ctx->denoised[0], (strides[0] * sd->height) << use_highbd);
                svt_memcpy_c(raw_data[1],
                    ctx->denoised[1],
                    (strides[1] * (sd->height >> chroma_sub_log2[0])) << use_highbd);
                svt_memcpy_c(raw_data[2],
                    ctx->denoised[2],
                    (strides[2] * (sd->height >> chroma_sub_log2[0])) << use_highbd);
            }
        } else
            unpack_2d_pic(ctx->denoised, sd);
    }
    svt_aom_flat_block_finder_free(&ctx->flat_block_finder);
    svt_aom_noise_model_free(&ctx->noise_model);
    free(ctx->flat_blocks);
    ctx->flat_blocks = NULL;

    return 1;
}
