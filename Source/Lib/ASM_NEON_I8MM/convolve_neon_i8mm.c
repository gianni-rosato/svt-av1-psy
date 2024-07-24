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

#include "common_dsp_rtcd.h"
#include "convolve_neon.h"
#include "definitions.h"
#include "sum_neon.h"
#include "mem_neon.h"

DECLARE_ALIGNED(16, static const uint8_t, kMatMulPermuteTbl[32]) = {
    // clang-format off
    0,  1,  2,  3,  4,  5,  6,  7,  2,  3,  4,  5,  6,  7,  8,  9,
    4,  5,  6,  7,  8,  9, 10, 11,  6,  7,  8,  9, 10, 11, 12, 13
    // clang-format on
};

static INLINE uint8x8_t convolve8_8_x(uint8x16_t samples, const int8x8_t filter, const uint8x16x3_t permute_tbl,
                                      const int32x4_t horiz_const) {
    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    // { 4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
    // { 8,  9, 10, 11,  9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14 }
    uint8x16_t perm_samples[3] = {vqtbl1q_u8(samples, permute_tbl.val[0]),
                                  vqtbl1q_u8(samples, permute_tbl.val[1]),
                                  vqtbl1q_u8(samples, permute_tbl.val[2])};

    int32x4_t sum0123 = vusdotq_lane_s32(horiz_const, perm_samples[0], filter, 0);
    sum0123           = vusdotq_lane_s32(sum0123, perm_samples[1], filter, 1);

    int32x4_t sum4567 = vusdotq_lane_s32(horiz_const, perm_samples[1], filter, 0);
    sum4567           = vusdotq_lane_s32(sum4567, perm_samples[2], filter, 1);

    int16x8_t sum_s16 = vcombine_s16(vmovn_s32(sum0123), vmovn_s32(sum4567));
    // We halved the convolution filter values so - 1 from the right shift.
    return vqrshrun_n_s16(sum_s16, FILTER_BITS - 1);
}

static INLINE void convolve_x_sr_8tap_neon_i8mm(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                                                ptrdiff_t dst_stride, int width, int height, const int16_t *filter_x,
                                                const int32x4_t horiz_const) {
    // Filter values are even, so halve to reduce intermediate precision reqs.
    const int8x8_t     x_filter    = vshrn_n_s16(vld1q_s16(filter_x), 1);
    const uint8x16x3_t permute_tbl = vld1q_u8_x3(kDotProdPermuteTbl);

    do {
        const uint8_t *s = src;
        uint8_t       *d = dst;
        int            w = width;

        do {
            uint8x16_t s0, s1, s2, s3;
            load_u8_16x4(s, src_stride, &s0, &s1, &s2, &s3);

            uint8x8_t d0 = convolve8_8_x(s0, x_filter, permute_tbl, horiz_const);
            uint8x8_t d1 = convolve8_8_x(s1, x_filter, permute_tbl, horiz_const);
            uint8x8_t d2 = convolve8_8_x(s2, x_filter, permute_tbl, horiz_const);
            uint8x8_t d3 = convolve8_8_x(s3, x_filter, permute_tbl, horiz_const);

            store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

            s += 8;
            d += 8;
            w -= 8;
        } while (w != 0);
        src += 4 * src_stride;
        dst += 4 * dst_stride;
        height -= 4;
    } while (height != 0);
}

static INLINE int16x4_t convolve6_4_x(uint8x16_t samples, const int8x16_t filter, const uint8x16_t permute_tbl,
                                      const int32x4_t horiz_const) {
    // Permute samples ready for matrix multiply.
    // { 0,  1,  2,  3,  4,  5,  6,  7,  2,  3,  4,  5,  6,  7,  8,  9 }
    uint8x16_t perm_samples = vqtbl1q_u8(samples, permute_tbl);

    // These instructions multiply a 2x8 matrix (samples) by an 8x2 matrix
    // (filter), destructively accumulating into the destination register.
    int32x4_t sum = vusmmlaq_s32(horiz_const, perm_samples, filter);

    // Further narrowing and packing is performed by the caller.
    return vmovn_s32(sum);
}

static INLINE uint8x8_t convolve6_8_x(uint8x16_t samples, const int8x16_t filter, const uint8x16x2_t permute_tbl,
                                      const int32x4_t horiz_const) {
    // Permute samples ready for matrix multiply.
    // { 0,  1,  2,  3,  4,  5,  6,  7,  2,  3,  4,  5,  6,  7,  8,  9 }
    // { 4,  5,  6,  7,  8,  9, 10, 11,  6,  7,  8,  9, 10, 11, 12, 13 }
    uint8x16_t perm_samples[2] = {vqtbl1q_u8(samples, permute_tbl.val[0]), vqtbl1q_u8(samples, permute_tbl.val[1])};

    // These instructions multiply a 2x8 matrix (samples) by an 8x2 matrix
    // (filter), destructively accumulating into the destination register.
    int32x4_t sum0123 = vusmmlaq_s32(horiz_const, perm_samples[0], filter);
    int32x4_t sum4567 = vusmmlaq_s32(horiz_const, perm_samples[1], filter);

    int16x8_t sum = vcombine_s16(vmovn_s32(sum0123), vmovn_s32(sum4567));
    // We halved the convolution filter values so - 1 from the right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

static INLINE void convolve_x_sr_6tap_neon_i8mm(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                                                ptrdiff_t dst_stride, int width, int height, const int16_t *filter_x,
                                                const int32x4_t horiz_const) {
    // Filter values are even, so halve to reduce intermediate precision reqs.
    const int8x8_t x_filter_s8 = vshrn_n_s16(vld1q_s16(filter_x), 1);
    // Stagger the filter for use with the matrix multiply instructions.
    // { f0, f1, f2, f3, f4, f5,  0,  0,  0, f0, f1, f2, f3, f4, f5,  0 }
    const int8x16_t x_filter = vcombine_s8(vext_s8(x_filter_s8, x_filter_s8, 1), x_filter_s8);

    if (width == 4) {
        const uint8x16_t permute_tbl = vld1q_u8(kMatMulPermuteTbl);
        do {
            uint8x16_t s0, s1, s2, s3;
            load_u8_16x4(src, src_stride, &s0, &s1, &s2, &s3);

            int16x4_t t0 = convolve6_4_x(s0, x_filter, permute_tbl, horiz_const);
            int16x4_t t1 = convolve6_4_x(s1, x_filter, permute_tbl, horiz_const);
            int16x4_t t2 = convolve6_4_x(s2, x_filter, permute_tbl, horiz_const);
            int16x4_t t3 = convolve6_4_x(s3, x_filter, permute_tbl, horiz_const);
            // We halved the filter values so -1 from right shift.
            uint8x8_t d01 = vqrshrun_n_s16(vcombine_s16(t0, t1), FILTER_BITS - 1);
            uint8x8_t d23 = vqrshrun_n_s16(vcombine_s16(t2, t3), FILTER_BITS - 1);

            store_u8x4_strided_x2(dst + 0 * dst_stride, dst_stride, d01);
            store_u8x4_strided_x2(dst + 2 * dst_stride, dst_stride, d23);

            src += 4 * src_stride;
            dst += 4 * dst_stride;
            height -= 4;
        } while (height != 0);
    } else {
        const uint8x16x2_t permute_tbl = vld1q_u8_x2(kMatMulPermuteTbl);
        do {
            const uint8_t *s = src;
            uint8_t       *d = dst;
            int            w = width;

            do {
                uint8x16_t s0, s1, s2, s3;
                load_u8_16x4(s, src_stride, &s0, &s1, &s2, &s3);

                uint8x8_t d0 = convolve6_8_x(s0, x_filter, permute_tbl, horiz_const);
                uint8x8_t d1 = convolve6_8_x(s1, x_filter, permute_tbl, horiz_const);
                uint8x8_t d2 = convolve6_8_x(s2, x_filter, permute_tbl, horiz_const);
                uint8x8_t d3 = convolve6_8_x(s3, x_filter, permute_tbl, horiz_const);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s += 8;
                d += 8;
                w -= 8;
            } while (w != 0);
            src += 4 * src_stride;
            dst += 4 * dst_stride;
            height -= 4;
        } while (height != 0);
    }
}

void svt_av1_convolve_x_sr_neon_i8mm(const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w, int h,
                                     InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y,
                                     const int subpel_x_qn, const int subpel_y_qn, ConvolveParams *conv_params) {
    if (w == 2 || h == 2) {
        svt_av1_convolve_x_sr_c(src,
                                src_stride,
                                dst,
                                dst_stride,
                                w,
                                h,
                                filter_params_x,
                                filter_params_y,
                                subpel_x_qn,
                                subpel_y_qn,
                                conv_params);
        return;
    }

    const uint8_t horiz_offset = filter_params_x->taps / 2 - 1;
    src -= horiz_offset;

    const int16_t *x_filter_ptr = av1_get_interp_filter_subpel_kernel(*filter_params_x, subpel_x_qn & SUBPEL_MASK);

    int filter_taps = get_filter_tap(filter_params_x, subpel_x_qn & SUBPEL_MASK);

    // A shim of 1 << (ROUND0_BITS - 1) enables us to simplify computation in the
    // convolution kernels: Adding this shim enables us to use a single rounding
    // right shift by FILTER_BITS instead of two rounding right shifts: first by
    // ROUND0_BITS, and then subsequently by FILTER_BITS - ROUND0_BITS.
    // Halve the total because we will halve the filter values.
    const int32x4_t horiz_const = vdupq_n_s32((1 << ((ROUND0_BITS - 1)) / 2));

    if (filter_taps <= 6) {
        convolve_x_sr_6tap_neon_i8mm(src + 1, src_stride, dst, dst_stride, w, h, x_filter_ptr, horiz_const);
        return;
    }

    convolve_x_sr_8tap_neon_i8mm(src, src_stride, dst, dst_stride, w, h, x_filter_ptr, horiz_const);
}
