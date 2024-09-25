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

// clang-format off
DECLARE_ALIGNED(16, static const uint8_t, kMatMulPermuteTbl[32]) = {
    0,  1,  2,  3,  4,  5,  6,  7,  2,  3,  4,  5,  6,  7,  8,  9,
    4,  5,  6,  7,  8,  9, 10, 11,  6,  7,  8,  9, 10, 11, 12, 13
};

DECLARE_ALIGNED(16, static const uint8_t, kDotProdPermuteTbl[48]) = {
    0, 1, 2,  3,  1, 2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6,
    4, 5, 6,  7,  5, 6,  7,  8,  6,  7,  8,  9,  7,  8,  9,  10,
    8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14
};
// clang-format on

static INLINE int16x4_t convolve8_4_h(const uint8x16_t samples, const int8x8_t filters,
                                      const uint8x16x2_t permute_tbl) {
    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    // { 4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
    uint8x16_t permuted_samples[2] = {vqtbl1q_u8(samples, permute_tbl.val[0]), vqtbl1q_u8(samples, permute_tbl.val[1])};

    int32x4_t sum = vusdotq_lane_s32(vdupq_n_s32(0), permuted_samples[0], filters, 0);
    sum           = vusdotq_lane_s32(sum, permuted_samples[1], filters, 1);

    // Further narrowing and packing is performed by the caller.
    return vqmovn_s32(sum);
}

static INLINE uint8x8_t convolve8_8_h(const uint8x16_t samples, const int8x8_t filters,
                                      const uint8x16x3_t permute_tbl) {
    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    // { 4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
    // { 8,  9, 10, 11,  9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14 }
    uint8x16_t permuted_samples[3] = {vqtbl1q_u8(samples, permute_tbl.val[0]),
                                      vqtbl1q_u8(samples, permute_tbl.val[1]),
                                      vqtbl1q_u8(samples, permute_tbl.val[2])};

    // First 4 output values.
    int32x4_t sum0 = vusdotq_lane_s32(vdupq_n_s32(0), permuted_samples[0], filters, 0);
    sum0           = vusdotq_lane_s32(sum0, permuted_samples[1], filters, 1);
    // Second 4 output values.
    int32x4_t sum1 = vusdotq_lane_s32(vdupq_n_s32(0), permuted_samples[1], filters, 0);
    sum1           = vusdotq_lane_s32(sum1, permuted_samples[2], filters, 1);

    // Narrow and re-pack.
    int16x8_t sum = vcombine_s16(vqmovn_s32(sum0), vqmovn_s32(sum1));
    return vqrshrun_n_s16(sum, FILTER_BITS);
}

static INLINE void convolve8_horiz_8tap_neon_i8mm(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
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

static INLINE int16x4_t convolve6_4_h(const uint8x16_t samples, const int8x16_t filter, const uint8x16_t permute_tbl) {
    // Permute samples ready for matrix multiply.
    // { 0,  1,  2,  3,  4,  5,  6,  7,  2,  3,  4,  5,  6,  7,  8,  9 }
    uint8x16_t perm_samples = vqtbl1q_u8(samples, permute_tbl);

    // These instructions multiply a 2x8 matrix (samples) by an 8x2 matrix
    // (filter), destructively accumulating into the destination register.
    int32x4_t sum = vusmmlaq_s32(vdupq_n_s32(0), perm_samples, filter);

    // Further narrowing and packing is performed by the caller.
    return vmovn_s32(sum);
}

static INLINE uint8x8_t convolve6_8_h(const uint8x16_t samples, const int8x16_t filter,
                                      const uint8x16x2_t permute_tbl) {
    // Permute samples ready for matrix multiply.
    // { 0,  1,  2,  3,  4,  5,  6,  7,  2,  3,  4,  5,  6,  7,  8,  9 }
    // { 4,  5,  6,  7,  8,  9, 10, 11,  6,  7,  8,  9, 10, 11, 12, 13 }
    uint8x16_t perm_samples[2] = {vqtbl1q_u8(samples, permute_tbl.val[0]), vqtbl1q_u8(samples, permute_tbl.val[1])};

    // These instructions multiply a 2x8 matrix (samples) by an 8x2 matrix
    // (filter), destructively accumulating into the destination register.
    int32x4_t sum0123 = vusmmlaq_s32(vdupq_n_s32(0), perm_samples[0], filter);
    int32x4_t sum4567 = vusmmlaq_s32(vdupq_n_s32(0), perm_samples[1], filter);

    // Narrow and re-pack.
    int16x8_t sum = vcombine_s16(vmovn_s32(sum0123), vmovn_s32(sum4567));
    // We halved the filter values so -1 from right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

static INLINE void convolve8_horiz_6tap_neon_i8mm(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                                                  ptrdiff_t dst_stride, const int16_t *filter_x, int width,
                                                  int height) {
    // Filter values are even, so halve to reduce intermediate precision reqs.
    const int8x8_t x_filter = vshrn_n_s16(vld1q_s16(filter_x), 1);
    // Stagger the filter for use with the matrix multiply instructions.
    // { f0, f1, f2, f3, f4, f5,  0,  0,  0, f0, f1, f2, f3, f4, f5,  0 }
    const int8x16_t filter = vcombine_s8(vext_s8(x_filter, x_filter, 1), x_filter);

    if (width == 4) {
        const uint8x16_t perm_tbl = vld1q_u8(kMatMulPermuteTbl);
        do {
            uint8x16_t s0, s1, s2, s3;
            load_u8_16x4(src, src_stride, &s0, &s1, &s2, &s3);

            int16x4_t t0 = convolve6_4_h(s0, filter, perm_tbl);
            int16x4_t t1 = convolve6_4_h(s1, filter, perm_tbl);
            int16x4_t t2 = convolve6_4_h(s2, filter, perm_tbl);
            int16x4_t t3 = convolve6_4_h(s3, filter, perm_tbl);
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
        const uint8x16x2_t perm_tbl = vld1q_u8_x2(kMatMulPermuteTbl);

        do {
            int            w = width;
            const uint8_t *s = src;
            uint8_t       *d = dst;
            do {
                uint8x16_t s0, s1, s2, s3;
                load_u8_16x4(s, src_stride, &s0, &s1, &s2, &s3);

                uint8x8_t d0 = convolve6_8_h(s0, filter, perm_tbl);
                uint8x8_t d1 = convolve6_8_h(s1, filter, perm_tbl);
                uint8x8_t d2 = convolve6_8_h(s2, filter, perm_tbl);
                uint8x8_t d3 = convolve6_8_h(s3, filter, perm_tbl);

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

void svt_aom_convolve8_horiz_neon_i8mm(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride,
                                       const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4,
                                       int w, int h) {
    assert((intptr_t)dst % 4 == 0);
    assert(dst_stride % 4 == 0);

    (void)x_step_q4;
    (void)filter_y;
    (void)y_step_q4;

    src -= (SUBPEL_TAPS / 2) - 1;

    int filter_taps = get_filter_taps_convolve8(filter_x);

    if (filter_taps == 2) {
        convolve8_horiz_2tap_neon(src + 3, src_stride, dst, dst_stride, filter_x, w, h);
    } else if (filter_taps <= 6) {
        convolve8_horiz_6tap_neon_i8mm(src + 1, src_stride, dst, dst_stride, filter_x, w, h);
    } else {
        convolve8_horiz_8tap_neon_i8mm(src, src_stride, dst, dst_stride, filter_x, w, h);
    }
}
