/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include "convolve.h"
#include "common_dsp_rtcd.h"

// Note: Fixed size intermediate buffers, place limits on parameters
// of some functions. 2d filtering proceeds in 2 steps:
//   (1) Interpolate horizontally into an intermediate buffer, temp.
//   (2) Interpolate temp vertically to derive the sub-pixel result.
// Deriving the maximum number of rows in the temp buffer (135):
// --Smallest scaling factor is x1/2 ==> y_step_q4 = 32 (Normative).
// --Largest block size is 128x128 pixels.
// --128 rows in the downscaled frame span a distance of (128 - 1) * 32 in the
//   original frame (in 1/16th pixel units).
// --Must round-up because block may be located at sub-pixel position.
// --Require an additional SUBPEL_TAPS rows for the 8-tap filter tails.
// --((128 - 1) * 32 + 15) >> 4 + 8 = 263.
#define WIENER_MAX_EXT_SIZE 263

static INLINE int32_t horz_scalar_product(const uint8_t *a, const int16_t *b) {
    int32_t sum = 0;
    for (int32_t k = 0; k < SUBPEL_TAPS; ++k) sum += a[k] * b[k];
    return sum;
}

static INLINE int32_t highbd_horz_scalar_product(const uint16_t *a, const int16_t *b) {
    int32_t sum = 0;
    for (int32_t k = 0; k < SUBPEL_TAPS; ++k) sum += a[k] * b[k];
    return sum;
}

static INLINE int32_t highbd_vert_scalar_product(const uint16_t *a, ptrdiff_t a_stride,
                                                 const int16_t *b) {
    int32_t sum = 0;
    for (int32_t k = 0; k < SUBPEL_TAPS; ++k) sum += a[k * a_stride] * b[k];
    return sum;
}

static const InterpKernel *get_filter_base(const int16_t *filter) {
    // NOTE: This assumes that the filter table is 256-byte aligned.
    return (const InterpKernel *)(((intptr_t)filter) & ~((intptr_t)0xFF));
}

static int32_t get_filter_offset(const int16_t *f, const InterpKernel *base) {
    return (int32_t)((const InterpKernel *)(intptr_t)f - base);
}

static void convolve_add_src_horiz_hip(const uint8_t *src, ptrdiff_t src_stride, uint16_t *dst,
                                       ptrdiff_t dst_stride, const InterpKernel *x_filters,
                                       int32_t x0_q4, int32_t x_step_q4, int32_t w, int32_t h,
                                       int32_t round0_bits) {
    const int32_t bd = 8;
    src -= SUBPEL_TAPS / 2 - 1;
    for (int32_t y = 0; y < h; ++y) {
        int32_t x_q4 = x0_q4;
        for (int32_t x = 0; x < w; ++x) {
            const uint8_t *const src_x    = &src[x_q4 >> SUBPEL_BITS];
            const int16_t *const x_filter = x_filters[x_q4 & SUBPEL_MASK];
            const int32_t        rounding = ((int32_t)src_x[SUBPEL_TAPS / 2 - 1] << FILTER_BITS) +
                                     (1 << (bd + FILTER_BITS - 1));
            const int32_t sum = horz_scalar_product(src_x, x_filter) + rounding;
            dst[x]            = (uint16_t)clamp(
                ROUND_POWER_OF_TWO(sum, round0_bits), 0, WIENER_CLAMP_LIMIT(round0_bits, bd) - 1);
            x_q4 += x_step_q4;
        }
        src += src_stride;
        dst += dst_stride;
    }
}

static void convolve_add_src_vert_hip(const uint16_t *src, ptrdiff_t src_stride, uint8_t *dst,
                                      ptrdiff_t dst_stride, const InterpKernel *y_filters,
                                      int32_t y0_q4, int32_t y_step_q4, int32_t w, int32_t h,
                                      int32_t round1_bits) {
    const int32_t bd = 8;
    src -= src_stride * (SUBPEL_TAPS / 2 - 1);

    for (int32_t x = 0; x < w; ++x) {
        int32_t y_q4 = y0_q4;
        for (int32_t y = 0; y < h; ++y) {
            const uint16_t *     src_y    = &src[(y_q4 >> SUBPEL_BITS) * src_stride];
            const int16_t *const y_filter = y_filters[y_q4 & SUBPEL_MASK];
            const int32_t        rounding =
                ((int32_t)src_y[(SUBPEL_TAPS / 2 - 1) * src_stride] << FILTER_BITS) -
                (1 << (bd + round1_bits - 1));
            const int32_t sum = highbd_vert_scalar_product(src_y, src_stride, y_filter) + rounding;
            dst[y * dst_stride] = clip_pixel(ROUND_POWER_OF_TWO(sum, round1_bits));
            y_q4 += y_step_q4;
        }
        ++src;
        ++dst;
    }
}

void eb_av1_wiener_convolve_add_src_c(const uint8_t *const src, const ptrdiff_t src_stride,
                                      uint8_t *const dst, const ptrdiff_t dst_stride,
                                      const int16_t *const filter_x, const int16_t *const filter_y,
                                      const int32_t w, const int32_t h,
                                      const ConvolveParams *const conv_params) {
    const int32_t             x_step_q4 = 16;
    const int32_t             y_step_q4 = 16;
    const InterpKernel *const filters_x = get_filter_base(filter_x);
    const int32_t             x0_q4     = get_filter_offset(filter_x, filters_x);

    const InterpKernel *const filters_y = get_filter_base(filter_y);
    const int32_t             y0_q4     = get_filter_offset(filter_y, filters_y);

    uint16_t      temp[WIENER_MAX_EXT_SIZE * MAX_SB_SIZE];
    const int32_t intermediate_height =
        (((h - 1) * y_step_q4 + y0_q4) >> SUBPEL_BITS) + SUBPEL_TAPS;

    assert(w <= MAX_SB_SIZE);
    assert(h <= MAX_SB_SIZE);
    assert(y_step_q4 <= 32);
    assert(x_step_q4 <= 32);

    convolve_add_src_horiz_hip(src - src_stride * (SUBPEL_TAPS / 2 - 1),
                               src_stride,
                               temp,
                               MAX_SB_SIZE,
                               filters_x,
                               x0_q4,
                               x_step_q4,
                               w,
                               intermediate_height,
                               conv_params->round_0);
    convolve_add_src_vert_hip(temp + MAX_SB_SIZE * (SUBPEL_TAPS / 2 - 1),
                              MAX_SB_SIZE,
                              dst,
                              dst_stride,
                              filters_y,
                              y0_q4,
                              y_step_q4,
                              w,
                              h,
                              conv_params->round_1);
}

static void highbd_convolve_add_src_horiz_hip(const uint8_t *src8, ptrdiff_t src_stride,
                                              uint16_t *dst, ptrdiff_t dst_stride,
                                              const InterpKernel *x_filters, int32_t x0_q4,
                                              int32_t x_step_q4, int32_t w, int32_t h,
                                              int32_t round0_bits, int32_t bd) {
    const int32_t extraprec_clamp_limit = WIENER_CLAMP_LIMIT(round0_bits, bd);
    uint16_t *    src                   = CONVERT_TO_SHORTPTR(src8);
    src -= SUBPEL_TAPS / 2 - 1;
    for (int32_t y = 0; y < h; ++y) {
        int32_t x_q4 = x0_q4;
        for (int32_t x = 0; x < w; ++x) {
            const uint16_t *const src_x    = &src[x_q4 >> SUBPEL_BITS];
            const int16_t *const  x_filter = x_filters[x_q4 & SUBPEL_MASK];
            const int32_t         rounding = ((int32_t)src_x[SUBPEL_TAPS / 2 - 1] << FILTER_BITS) +
                                     (1 << (bd + FILTER_BITS - 1));
            const int32_t sum = highbd_horz_scalar_product(src_x, x_filter) + rounding;
            dst[x] =
                (uint16_t)clamp(ROUND_POWER_OF_TWO(sum, round0_bits), 0, extraprec_clamp_limit - 1);
            x_q4 += x_step_q4;
        }
        src += src_stride;
        dst += dst_stride;
    }
}

static void highbd_convolve_add_src_vert_hip(const uint16_t *src, ptrdiff_t src_stride,
                                             uint8_t *dst8, ptrdiff_t dst_stride,
                                             const InterpKernel *y_filters, int32_t y0_q4,
                                             int32_t y_step_q4, int32_t w, int32_t h,
                                             int32_t round1_bits, int32_t bd) {
    uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
    src -= src_stride * (SUBPEL_TAPS / 2 - 1);
    for (int32_t x = 0; x < w; ++x) {
        int32_t y_q4 = y0_q4;
        for (int32_t y = 0; y < h; ++y) {
            const uint16_t *     src_y    = &src[(y_q4 >> SUBPEL_BITS) * src_stride];
            const int16_t *const y_filter = y_filters[y_q4 & SUBPEL_MASK];
            const int32_t        rounding =
                ((int32_t)src_y[(SUBPEL_TAPS / 2 - 1) * src_stride] << FILTER_BITS) -
                (1 << (bd + round1_bits - 1));
            const int32_t sum = highbd_vert_scalar_product(src_y, src_stride, y_filter) + rounding;
            dst[y * dst_stride] = clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, round1_bits), bd);
            y_q4 += y_step_q4;
        }
        ++src;
        ++dst;
    }
}

void eb_av1_highbd_wiener_convolve_add_src_c(
    const uint8_t *const src, const ptrdiff_t src_stride, uint8_t *const dst,
    const ptrdiff_t dst_stride, const int16_t *const filter_x, const int16_t *const filter_y,
    const int32_t w, const int32_t h, const ConvolveParams *const conv_params, const int32_t bd) {
    const int32_t             x_step_q4 = 16;
    const int32_t             y_step_q4 = 16;
    const InterpKernel *const filters_x = get_filter_base(filter_x);
    const int32_t             x0_q4     = get_filter_offset(filter_x, filters_x);

    const InterpKernel *const filters_y = get_filter_base(filter_y);
    const int32_t             y0_q4     = get_filter_offset(filter_y, filters_y);

    uint16_t      temp[WIENER_MAX_EXT_SIZE * MAX_SB_SIZE];
    const int32_t intermediate_height =
        (((h - 1) * y_step_q4 + y0_q4) >> SUBPEL_BITS) + SUBPEL_TAPS;

    assert(w <= MAX_SB_SIZE);
    assert(h <= MAX_SB_SIZE);
    assert(y_step_q4 <= 32);
    assert(x_step_q4 <= 32);
    assert(bd + FILTER_BITS - conv_params->round_0 + 2 <= 16);

    highbd_convolve_add_src_horiz_hip(src - src_stride * (SUBPEL_TAPS / 2 - 1),
                                      src_stride,
                                      temp,
                                      MAX_SB_SIZE,
                                      filters_x,
                                      x0_q4,
                                      x_step_q4,
                                      w,
                                      intermediate_height,
                                      conv_params->round_0,
                                      bd);
    highbd_convolve_add_src_vert_hip(temp + MAX_SB_SIZE * (SUBPEL_TAPS / 2 - 1),
                                     MAX_SB_SIZE,
                                     dst,
                                     dst_stride,
                                     filters_y,
                                     y0_q4,
                                     y_step_q4,
                                     w,
                                     h,
                                     conv_params->round_1,
                                     bd);
}
static INLINE int vert_scalar_product(const uint8_t *a, ptrdiff_t a_stride, const int16_t *b) {
    int sum = 0;
    for (int k = 0; k < SUBPEL_TAPS; ++k) sum += a[k * a_stride] * b[k];
    return sum;
}

static void convolve_horiz(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                           ptrdiff_t dst_stride, const InterpKernel *x_filters, int x0_q4,
                           int x_step_q4, int w, int h) {
    src -= SUBPEL_TAPS / 2 - 1;
    for (int y = 0; y < h; ++y) {
        int x_q4 = x0_q4;
        for (int x = 0; x < w; ++x) {
            const uint8_t *const src_x    = &src[x_q4 >> SUBPEL_BITS];
            const int16_t *const x_filter = x_filters[x_q4 & SUBPEL_MASK];
            const int            sum      = horz_scalar_product(src_x, x_filter);
            dst[x]                        = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
            x_q4 += x_step_q4;
        }
        src += src_stride;
        dst += dst_stride;
    }
}

static void convolve_vert(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                          ptrdiff_t dst_stride, const InterpKernel *y_filters, int y0_q4,
                          int y_step_q4, int w, int h) {
    src -= src_stride * (SUBPEL_TAPS / 2 - 1);

    for (int x = 0; x < w; ++x) {
        int y_q4 = y0_q4;
        for (int y = 0; y < h; ++y) {
            const unsigned char *src_y    = &src[(y_q4 >> SUBPEL_BITS) * src_stride];
            const int16_t *const y_filter = y_filters[y_q4 & SUBPEL_MASK];
            const int            sum      = vert_scalar_product(src_y, src_stride, y_filter);
            dst[y * dst_stride]           = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
            y_q4 += y_step_q4;
        }
        ++src;
        ++dst;
    }
}

void aom_convolve8_horiz_c(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                           ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4,
                           const int16_t *filter_y, int y_step_q4, int w, int h) {
    const InterpKernel *const filters_x = get_filter_base(filter_x);
    const int                 x0_q4     = get_filter_offset(filter_x, filters_x);

    (void)filter_y;
    (void)y_step_q4;

    convolve_horiz(src, src_stride, dst, dst_stride, filters_x, x0_q4, x_step_q4, w, h);
}

void aom_convolve8_vert_c(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                          ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4,
                          const int16_t *filter_y, int y_step_q4, int w, int h) {
    const InterpKernel *const filters_y = get_filter_base(filter_y);
    const int                 y0_q4     = get_filter_offset(filter_y, filters_y);

    (void)filter_x;
    (void)x_step_q4;

    convolve_vert(src, src_stride, dst, dst_stride, filters_y, y0_q4, y_step_q4, w, h);
}
static INLINE const int16_t *av1_get_interp_filter_subpel_kernel(
    const InterpFilterParams filter_params, const int32_t subpel);
