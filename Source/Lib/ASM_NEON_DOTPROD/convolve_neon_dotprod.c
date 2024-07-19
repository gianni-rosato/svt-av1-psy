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
#include "definitions.h"
#include "filter.h"
#include "mem_neon.h"
#include "utility.h"

DECLARE_ALIGNED(16, static const uint8_t, kDotProdPermuteTbl[48]) = {
    0, 1, 2, 3, 1, 2, 3, 4,  2, 3, 4,  5,  3, 4,  5,  6,  4,  5,  6,  7,  5,  6,  7,  8,
    6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14};

static INLINE int16x4_t convolve4_4_x(const uint8x16_t samples, const int8x8_t filters, const uint8x16_t permute_tbl) {
    // Transform sample range to [-128, 127] for 8-bit signed dot product.
    int8x16_t samples_128 = vreinterpretq_s8_u8(vsubq_u8(samples, vdupq_n_u8(128)));

    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    int8x16_t perm_samples = vqtbl1q_s8(samples_128, permute_tbl);

    // Dot product constants:
    // Accumulate into 128 << FILTER_BITS to account for range transform.
    // Adding a shim of 1 << (ROUND0_BITS - 1) enables us to use a single rounding
    // right shift by FILTER_BITS - instead of a first rounding right shift by
    // ROUND0_BITS, followed by second rounding right shift by FILTER_BITS -
    // ROUND0_BITS. Halve the total because we halved the filter values.
    int32x4_t acc = vdupq_n_s32(((128 << FILTER_BITS) + (1 << ((ROUND0_BITS - 1)))) / 2);
    int32x4_t sum = vdotq_lane_s32(acc, perm_samples, filters, 0);

    // Further narrowing and packing is performed by the caller.
    return vmovn_s32(sum);
}

static INLINE uint8x8_t convolve4_8_x(const uint8x16_t samples, const int8x8_t filters,
                                      const uint8x16x2_t permute_tbl) {
    // Transform sample range to [-128, 127] for 8-bit signed dot product.
    int8x16_t samples_128 = vreinterpretq_s8_u8(vsubq_u8(samples, vdupq_n_u8(128)));

    // Permute samples ready for dot product.
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    // { 4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
    int8x16_t perm_samples[2] = {vqtbl1q_s8(samples_128, permute_tbl.val[0]),
                                 vqtbl1q_s8(samples_128, permute_tbl.val[1])};

    // Dot product constants:
    // Accumulate into 128 << FILTER_BITS to account for range transform.
    // Adding a shim of 1 << (ROUND0_BITS - 1) enables us to use a single rounding
    // right shift by FILTER_BITS - instead of a first rounding right shift by
    // ROUND0_BITS, followed by second rounding right shift by FILTER_BITS -
    // ROUND0_BITS. Halve the total because we halved the filter values.
    int32x4_t acc = vdupq_n_s32(((128 << FILTER_BITS) + (1 << ((ROUND0_BITS - 1)))) / 2);

    int32x4_t sum0123 = vdotq_lane_s32(acc, perm_samples[0], filters, 0);
    int32x4_t sum4567 = vdotq_lane_s32(acc, perm_samples[1], filters, 0);

    // Narrow and re-pack.
    int16x8_t sum = vcombine_s16(vmovn_s32(sum0123), vmovn_s32(sum4567));
    // We halved the filter values so -1 from right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

static INLINE void convolve_x_sr_4tap_neon_dotprod(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                                                   ptrdiff_t dst_stride, int width, int height,
                                                   const int16_t *filter_x) {
    const int16x4_t x_filter = vld1_s16(filter_x + 2);
    // All 4-tap and bilinear filter values are even, so halve them to reduce
    // intermediate precision requirements.
    const int8x8_t filter = vshrn_n_s16(vcombine_s16(x_filter, vdup_n_s16(0)), 1);

    if (width == 4) {
        const uint8x16_t permute_tbl = vld1q_u8(kDotProdPermuteTbl);

        do {
            uint8x16_t s0, s1, s2, s3;
            load_u8_16x4(src, src_stride, &s0, &s1, &s2, &s3);

            int16x4_t t0 = convolve4_4_x(s0, filter, permute_tbl);
            int16x4_t t1 = convolve4_4_x(s1, filter, permute_tbl);
            int16x4_t t2 = convolve4_4_x(s2, filter, permute_tbl);
            int16x4_t t3 = convolve4_4_x(s3, filter, permute_tbl);
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
        const uint8x16x2_t permute_tbl = vld1q_u8_x2(kDotProdPermuteTbl);

        do {
            const uint8_t *s = src;
            uint8_t       *d = dst;
            int            w = width;

            do {
                uint8x16_t s0, s1, s2, s3;
                load_u8_16x4(s, src_stride, &s0, &s1, &s2, &s3);

                uint8x8_t d0 = convolve4_8_x(s0, filter, permute_tbl);
                uint8x8_t d1 = convolve4_8_x(s1, filter, permute_tbl);
                uint8x8_t d2 = convolve4_8_x(s2, filter, permute_tbl);
                uint8x8_t d3 = convolve4_8_x(s3, filter, permute_tbl);

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

static INLINE uint8x8_t convolve8_8_x(uint8x16_t samples, const int8x8_t filter, const uint8x16x3_t permute_tbl) {
    // Transform sample range to [-128, 127] for 8-bit signed dot product.
    int8x16_t samples_128 = vreinterpretq_s8_u8(vsubq_u8(samples, vdupq_n_u8(128)));

    // Permute samples ready for dot product. */
    // { 0,  1,  2,  3,  1,  2,  3,  4,  2,  3,  4,  5,  3,  4,  5,  6 }
    // { 4,  5,  6,  7,  5,  6,  7,  8,  6,  7,  8,  9,  7,  8,  9, 10 }
    // { 8,  9, 10, 11,  9, 10, 11, 12, 10, 11, 12, 13, 11, 12, 13, 14 }
    int8x16_t perm_samples[3] = {vqtbl1q_s8(samples_128, permute_tbl.val[0]),
                                 vqtbl1q_s8(samples_128, permute_tbl.val[1]),
                                 vqtbl1q_s8(samples_128, permute_tbl.val[2])};

    // Dot product constants:
    // Accumulate into 128 << FILTER_BITS to account for range transform.
    // Adding a shim of 1 << (ROUND0_BITS - 1) enables us to use a single rounding
    // right shift by FILTER_BITS - instead of a first rounding right shift by
    // ROUND0_BITS, followed by second rounding right shift by FILTER_BITS -
    // ROUND0_BITS. Halve the total because we halved the filter values.
    int32x4_t acc = vdupq_n_s32(((128 << FILTER_BITS) + (1 << ((ROUND0_BITS - 1)))) / 2);

    int32x4_t sum0123 = vdotq_lane_s32(acc, perm_samples[0], filter, 0);
    sum0123           = vdotq_lane_s32(sum0123, perm_samples[1], filter, 1);

    int32x4_t sum4567 = vdotq_lane_s32(acc, perm_samples[1], filter, 0);
    sum4567           = vdotq_lane_s32(sum4567, perm_samples[2], filter, 1);

    // Narrow and re-pack.
    int16x8_t sum_s16 = vcombine_s16(vmovn_s32(sum0123), vmovn_s32(sum4567));
    // We halved the convolution filter values so - 1 from the right shift.
    return vqrshrun_n_s16(sum_s16, FILTER_BITS - 1);
}

void svt_av1_convolve_x_sr_neon_dotprod(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride,
                                        int32_t w, int32_t h, InterpFilterParams *filter_params_x,
                                        InterpFilterParams *filter_params_y, const int32_t subpel_x_qn,
                                        const int32_t subpel_y_qn, ConvolveParams *conv_params) {
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

    if (filter_taps <= 4) {
        convolve_x_sr_4tap_neon_dotprod(src + 2, src_stride, dst, dst_stride, w, h, x_filter_ptr);
        return;
    }

    const int16x8_t x_filter_s16 = vld1q_s16(x_filter_ptr);

    const uint8x16x3_t permute_tbl = vld1q_u8_x3(kDotProdPermuteTbl);
    // Filter values are even, so halve to reduce intermediate precision reqs.
    const int8x8_t x_filter = vshrn_n_s16(x_filter_s16, 1);

    do {
        int            width = w;
        const uint8_t *s     = src;
        uint8_t       *d     = dst;

        do {
            uint8x16_t s0, s1, s2, s3;
            load_u8_16x4(s, src_stride, &s0, &s1, &s2, &s3);

            uint8x8_t d0 = convolve8_8_x(s0, x_filter, permute_tbl);
            uint8x8_t d1 = convolve8_8_x(s1, x_filter, permute_tbl);
            uint8x8_t d2 = convolve8_8_x(s2, x_filter, permute_tbl);
            uint8x8_t d3 = convolve8_8_x(s3, x_filter, permute_tbl);

            store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

            s += 8;
            d += 8;
            width -= 8;
        } while (width != 0);
        src += 4 * src_stride;
        dst += 4 * dst_stride;
        h -= 4;
    } while (h != 0);
}
