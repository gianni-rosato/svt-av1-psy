/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>

#include "aom_convolve8_neon.h"
#include "common_dsp_rtcd.h"
#include "definitions.h"
#include "filter.h"
#include "mem_neon.h"
#include "transpose_neon.h"

// Filter values always sum to 128.
#define FILTER_WEIGHT 128

// clang-format off
DECLARE_ALIGNED(16, static const uint8_t, kDotProdPermuteTbl[48]) = {
    0, 1, 2,  3,  1, 2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6,
    4, 5, 6,  7,  5, 6,  7,  8,  6,  7,  8,  9,  7,  8,  9,  10,
    8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14
};
// clang-format on

static INLINE int16x4_t convolve8_4_h(const uint8x16_t samples, const int8x8_t filters,
                                      const uint8x16x2_t permute_tbl) {
    // Transform sample range to [-128, 127] for 8-bit signed dot product.
    int8x16_t samples_128 = vreinterpretq_s8_u8(vsubq_u8(samples, vdupq_n_u8(128)));

    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    // { 4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
    int8x16_t perm_samples[2] = {vqtbl1q_s8(samples_128, permute_tbl.val[0]),
                                 vqtbl1q_s8(samples_128, permute_tbl.val[1])};

    // Accumulate into 128 * FILTER_WEIGHT to account for range transform.
    int32x4_t acc = vdupq_n_s32(128 * FILTER_WEIGHT);
    int32x4_t sum = vdotq_lane_s32(acc, perm_samples[0], filters, 0);
    sum           = vdotq_lane_s32(sum, perm_samples[1], filters, 1);

    // Further narrowing and packing is performed by the caller.
    return vqmovn_s32(sum);
}

static INLINE uint8x8_t convolve8_8_h(const uint8x16_t samples, const int8x8_t filters,
                                      const uint8x16x3_t permute_tbl) {
    // Transform sample range to [-128, 127] for 8-bit signed dot product.
    int8x16_t samples_128 = vreinterpretq_s8_u8(vsubq_u8(samples, vdupq_n_u8(128)));

    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    // { 4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
    // { 8,  9, 10, 11,  9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14 }
    int8x16_t perm_samples[3] = {vqtbl1q_s8(samples_128, permute_tbl.val[0]),
                                 vqtbl1q_s8(samples_128, permute_tbl.val[1]),
                                 vqtbl1q_s8(samples_128, permute_tbl.val[2])};

    // Accumulate into 128 * FILTER_WEIGHT to account for range transform.
    int32x4_t acc = vdupq_n_s32(128 * FILTER_WEIGHT);
    // First 4 output values.
    int32x4_t sum0 = vdotq_lane_s32(acc, perm_samples[0], filters, 0);
    sum0           = vdotq_lane_s32(sum0, perm_samples[1], filters, 1);
    // Second 4 output values.
    int32x4_t sum1 = vdotq_lane_s32(acc, perm_samples[1], filters, 0);
    sum1           = vdotq_lane_s32(sum1, perm_samples[2], filters, 1);

    // Narrow and re-pack.
    int16x8_t sum = vcombine_s16(vqmovn_s32(sum0), vqmovn_s32(sum1));
    return vqrshrun_n_s16(sum, FILTER_BITS);
}

static INLINE void convolve8_horiz_8tap_neon_dotprod(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                                                     ptrdiff_t dst_stride, const int16_t *filter_x, int w, int h) {
    const int8x8_t filter = vmovn_s16(vld1q_s16(filter_x));

    if (w == 4) {
        const uint8x16x2_t perm_tbl = vld1q_u8_x2(kDotProdPermuteTbl);
        do {
            uint8x16_t s0, s1, s2, s3;
            load_u8_16x4(src, src_stride, &s0, &s1, &s2, &s3);

            int16x4_t d0  = convolve8_4_h(s0, filter, perm_tbl);
            int16x4_t d1  = convolve8_4_h(s1, filter, perm_tbl);
            int16x4_t d2  = convolve8_4_h(s2, filter, perm_tbl);
            int16x4_t d3  = convolve8_4_h(s3, filter, perm_tbl);
            uint8x8_t d01 = vqrshrun_n_s16(vcombine_s16(d0, d1), FILTER_BITS);
            uint8x8_t d23 = vqrshrun_n_s16(vcombine_s16(d2, d3), FILTER_BITS);

            store_u8x4_strided_x2(dst + 0 * dst_stride, dst_stride, d01);
            store_u8x4_strided_x2(dst + 2 * dst_stride, dst_stride, d23);

            src += 4 * src_stride;
            dst += 4 * dst_stride;
            h -= 4;
        } while (h > 0);
    } else {
        const uint8x16x3_t perm_tbl = vld1q_u8_x3(kDotProdPermuteTbl);

        do {
            int            width = w;
            const uint8_t *s     = src;
            uint8_t       *d     = dst;
            do {
                uint8x16_t s0, s1, s2, s3;
                load_u8_16x4(s, src_stride, &s0, &s1, &s2, &s3);

                uint8x8_t d0 = convolve8_8_h(s0, filter, perm_tbl);
                uint8x8_t d1 = convolve8_8_h(s1, filter, perm_tbl);
                uint8x8_t d2 = convolve8_8_h(s2, filter, perm_tbl);
                uint8x8_t d3 = convolve8_8_h(s3, filter, perm_tbl);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s += 8;
                d += 8;
                width -= 8;
            } while (width != 0);
            src += 4 * src_stride;
            dst += 4 * dst_stride;
            h -= 4;
        } while (h > 0);
    }
}

static INLINE int16x4_t convolve4_4_h(const uint8x16_t samples, const int8x8_t filters, const uint8x16_t permute_tbl) {
    // Transform sample range to [-128, 127] for 8-bit signed dot product.
    int8x16_t samples_128 = vreinterpretq_s8_u8(vsubq_u8(samples, vdupq_n_u8(128)));

    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    int8x16_t perm_samples = vqtbl1q_s8(samples_128, permute_tbl);

    // Accumulate into 128 * FILTER_WEIGHT to account for range transform.
    // (Divide by 2 since we halved the filter values.)
    int32x4_t acc = vdupq_n_s32(128 * FILTER_WEIGHT / 2);
    int32x4_t sum = vdotq_lane_s32(acc, perm_samples, filters, 0);

    // Further narrowing and packing is performed by the caller.
    return vmovn_s32(sum);
}

static INLINE uint8x8_t convolve4_8_h(const uint8x16_t samples, const int8x8_t filters,
                                      const uint8x16x2_t permute_tbl) {
    // Transform sample range to [-128, 127] for 8-bit signed dot product.
    int8x16_t samples_128 = vreinterpretq_s8_u8(vsubq_u8(samples, vdupq_n_u8(128)));

    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    // { 4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
    int8x16_t perm_samples[2] = {vqtbl1q_s8(samples_128, permute_tbl.val[0]),
                                 vqtbl1q_s8(samples_128, permute_tbl.val[1])};

    // Accumulate into 128 * FILTER_WEIGHT to account for range transform.
    // (Divide by 2 since we halved the filter values.)
    int32x4_t acc = vdupq_n_s32(128 * FILTER_WEIGHT / 2);
    // First 4 output values.
    int32x4_t sum0 = vdotq_lane_s32(acc, perm_samples[0], filters, 0);
    // Second 4 output values.
    int32x4_t sum1 = vdotq_lane_s32(acc, perm_samples[1], filters, 0);

    // Narrow and re-pack.
    int16x8_t sum = vcombine_s16(vmovn_s32(sum0), vmovn_s32(sum1));
    // We halved the filter values so -1 from right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

static INLINE void convolve8_horiz_4tap_neon_dotprod(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                                                     ptrdiff_t dst_stride, const int16_t *filter_x, int width,
                                                     int height) {
    const int16x4_t x_filter = vld1_s16(filter_x + 2);
    // All 4-tap and bilinear filter values are even, so halve them to reduce
    // intermediate precision requirements.
    const int8x8_t filter = vshrn_n_s16(vcombine_s16(x_filter, vdup_n_s16(0)), 1);

    if (width == 4) {
        const uint8x16_t permute_tbl = vld1q_u8(kDotProdPermuteTbl);

        do {
            uint8x16_t s0, s1, s2, s3;
            load_u8_16x4(src, src_stride, &s0, &s1, &s2, &s3);

            int16x4_t t0 = convolve4_4_h(s0, filter, permute_tbl);
            int16x4_t t1 = convolve4_4_h(s1, filter, permute_tbl);
            int16x4_t t2 = convolve4_4_h(s2, filter, permute_tbl);
            int16x4_t t3 = convolve4_4_h(s3, filter, permute_tbl);
            // We halved the filter values so -1 from right shift.
            uint8x8_t d01 = vqrshrun_n_s16(vcombine_s16(t0, t1), FILTER_BITS - 1);
            uint8x8_t d23 = vqrshrun_n_s16(vcombine_s16(t2, t3), FILTER_BITS - 1);

            store_u8x4_strided_x2(dst + 0 * dst_stride, dst_stride, d01);
            store_u8x4_strided_x2(dst + 2 * dst_stride, dst_stride, d23);

            src += 4 * src_stride;
            dst += 4 * dst_stride;
            height -= 4;
        } while (height > 0);
    } else {
        const uint8x16x2_t permute_tbl = vld1q_u8_x2(kDotProdPermuteTbl);

        do {
            const uint8_t *s = src;
            uint8_t       *d = dst;
            int            w = width;

            do {
                uint8x16_t s0, s1, s2, s3;
                load_u8_16x4(s, src_stride, &s0, &s1, &s2, &s3);

                uint8x8_t d0 = convolve4_8_h(s0, filter, permute_tbl);
                uint8x8_t d1 = convolve4_8_h(s1, filter, permute_tbl);
                uint8x8_t d2 = convolve4_8_h(s2, filter, permute_tbl);
                uint8x8_t d3 = convolve4_8_h(s3, filter, permute_tbl);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s += 8;
                d += 8;
                w -= 8;
            } while (w != 0);
            src += 4 * src_stride;
            dst += 4 * dst_stride;
            height -= 4;
        } while (height > 0);
    }
}

void svt_aom_convolve8_horiz_neon_dotprod(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride,
                                          const int16_t *filter_x, int x_step_q4, const int16_t *filter_y,
                                          int y_step_q4, int w, int h) {
    assert((intptr_t)dst % 4 == 0);
    assert(dst_stride % 4 == 0);

    (void)x_step_q4;
    (void)filter_y;
    (void)y_step_q4;

    src -= (SUBPEL_TAPS / 2) - 1;

    int filter_taps = get_filter_taps_convolve8(filter_x);

    if (filter_taps == 2) {
        convolve8_horiz_2tap_neon(src + 3, src_stride, dst, dst_stride, filter_x, w, h);
    } else if (filter_taps == 4) {
        convolve8_horiz_4tap_neon_dotprod(src + 2, src_stride, dst, dst_stride, filter_x, w, h);
    } else {
        convolve8_horiz_8tap_neon_dotprod(src, src_stride, dst, dst_stride, filter_x, w, h);
    }
}
