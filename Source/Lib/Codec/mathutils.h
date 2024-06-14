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

#ifndef AOM_AV1_ENCODER_MATHUTILS_H_
#define AOM_AV1_ENCODER_MATHUTILS_H_

#include <memory.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// Solves Ax = b, where x and b are column vectors of size nx1 and A is nxn
static INLINE int32_t linsolve(int32_t n, double *A, int32_t stride, double *b, double *x) {
    const double tiny_near_zero = 1.0E-16;
    int32_t      i, j, k;
    double       c;
    // Forward elimination
    for (k = 0; k < n - 1; k++) {
        // Bring the largest magnitude to the diagonal position
        for (i = n - 1; i > k; i--) {
            if (fabs(A[(i - 1) * stride + k]) < fabs(A[i * stride + k])) {
                for (j = 0; j < n; j++) {
                    c                       = A[i * stride + j];
                    A[i * stride + j]       = A[(i - 1) * stride + j];
                    A[(i - 1) * stride + j] = c;
                }
                c        = b[i];
                b[i]     = b[i - 1];
                b[i - 1] = c;
            }
        }
        for (i = k; i < n - 1; i++) {
            if (fabs(A[k * stride + k]) < tiny_near_zero)
                return 0;
            c = A[(i + 1) * stride + k] / A[k * stride + k];
            for (j = 0; j < n; j++) A[(i + 1) * stride + j] -= c * A[k * stride + j];
            b[i + 1] -= c * b[k];
        }
    }
    // Backward substitution
    for (i = n - 1; i >= 0; i--) {
        if (fabs(A[i * stride + i]) < tiny_near_zero)
            return 0;
        c = 0;
        for (j = i + 1; j <= n - 1; j++) c += A[i * stride + j] * x[j];
        x[i] = (b[i] - c) / A[i * stride + i];
    }

    return 1;
}

////////////////////////////////////////////////////////////////////////////////
// Least-squares
// Solves for n-dim x in a least squares sense to minimize |Ax - b|^2
// The solution is simply x = (A'A)^-1 A'b or simply the solution for
// the system: A'A x = A'b
static INLINE int32_t least_squares(int32_t n, double *A, int32_t rows, int32_t stride, double *b, double *scratch,
                                    double *x) {
    int32_t i, j, k;
    double *scratch_ = NULL;
    double *at_a, *atb;
    if (!scratch) {
        scratch_ = (double *)malloc(sizeof(*scratch) * n * (n + 1));
        scratch  = scratch_;
    }
    at_a = scratch;
    atb  = scratch + n * n;
    assert(at_a);
    for (i = 0; i < n; ++i) {
        for (j = i; j < n; ++j) {
            at_a[i * n + j] = 0.0;
            for (k = 0; k < rows; ++k) at_a[i * n + j] += A[k * stride + i] * A[k * stride + j];
            at_a[j * n + i] = at_a[i * n + j];
        }
        atb[i] = 0;
        for (k = 0; k < rows; ++k) atb[i] += A[k * stride + i] * b[k];
    }
    int32_t ret = linsolve(n, at_a, n, atb, x);
    if (scratch_)
        free(scratch_);
    return ret;
}

// Matrix multiply
static INLINE void multiply_mat(const double *m1, const double *m2, double *res, const int32_t m1_rows,
                                const int32_t inner_dim, const int32_t m2_cols) {
    double sum;

    int32_t row, col, inner;
    for (row = 0; row < m1_rows; ++row) {
        for (col = 0; col < m2_cols; ++col) {
            sum = 0;
            for (inner = 0; inner < inner_dim; ++inner) sum += m1[row * inner_dim + inner] * m2[inner * m2_cols + col];
            *(res++) = sum;
        }
    }
}

#endif // AOM_AV1_ENCODER_MATHUTILS_H_
