/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_CONVOLVE8_NEON_H_
#define AOM_CONVOLVE8_NEON_H_

#include <arm_neon.h>

#include "common_dsp_rtcd.h"
#include "definitions.h"
#include "filter.h"
#include "mem_neon.h"
#include "transpose_neon.h"

static INLINE int get_filter_taps_convolve8(const int16_t *filter) {
    if (filter[0] | filter[7]) {
        return 8;
    }
    if (filter[1] | filter[6]) {
        return 6;
    }
    if (filter[2] | filter[5]) {
        return 4;
    }
    return 2;
}

static INLINE int16x4_t convolve8_4_x(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
                                      const int16x4_t s4, const int16x4_t s5, const int16x4_t s6, const int16x4_t s7,
                                      const int16x8_t filter) {
    const int16x4_t filter_lo = vget_low_s16(filter);
    const int16x4_t filter_hi = vget_high_s16(filter);

    int16x4_t sum = vmul_lane_s16(s0, filter_lo, 0);
    sum           = vmla_lane_s16(sum, s1, filter_lo, 1);
    sum           = vmla_lane_s16(sum, s2, filter_lo, 2);
    sum           = vmla_lane_s16(sum, s3, filter_lo, 3);
    sum           = vmla_lane_s16(sum, s4, filter_hi, 0);
    sum           = vmla_lane_s16(sum, s5, filter_hi, 1);
    sum           = vmla_lane_s16(sum, s6, filter_hi, 2);
    sum           = vmla_lane_s16(sum, s7, filter_hi, 3);

    return sum;
}

static INLINE uint8x8_t convolve8_8_x(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
                                      const int16x8_t s4, const int16x8_t s5, const int16x8_t s6, const int16x8_t s7,
                                      const int16x8_t filter) {
    const int16x4_t filter_lo = vget_low_s16(filter);
    const int16x4_t filter_hi = vget_high_s16(filter);

    int16x8_t sum = vmulq_lane_s16(s0, filter_lo, 0);
    sum           = vmlaq_lane_s16(sum, s1, filter_lo, 1);
    sum           = vmlaq_lane_s16(sum, s2, filter_lo, 2);
    sum           = vmlaq_lane_s16(sum, s3, filter_lo, 3);
    sum           = vmlaq_lane_s16(sum, s4, filter_hi, 0);
    sum           = vmlaq_lane_s16(sum, s5, filter_hi, 1);
    sum           = vmlaq_lane_s16(sum, s6, filter_hi, 2);
    sum           = vmlaq_lane_s16(sum, s7, filter_hi, 3);

    // We halved the filter values so -1 from right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

static INLINE void convolve8_horiz_2tap_neon(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst,
                                             ptrdiff_t dst_stride, const int16_t *filter_x, int w, int h) {
    // Bilinear filter values are all positive.
    const uint8x8_t f0 = vdup_n_u8((uint8_t)filter_x[3]);
    const uint8x8_t f1 = vdup_n_u8((uint8_t)filter_x[4]);

    if (w == 4) {
        do {
            uint8x8_t s0 = load_unaligned_u8(src + 0 * src_stride + 0, src_stride);
            uint8x8_t s1 = load_unaligned_u8(src + 0 * src_stride + 1, src_stride);
            uint8x8_t s2 = load_unaligned_u8(src + 2 * src_stride + 0, src_stride);
            uint8x8_t s3 = load_unaligned_u8(src + 2 * src_stride + 1, src_stride);

            uint16x8_t sum0 = vmull_u8(s0, f0);
            sum0            = vmlal_u8(sum0, s1, f1);
            uint16x8_t sum1 = vmull_u8(s2, f0);
            sum1            = vmlal_u8(sum1, s3, f1);

            uint8x8_t d0 = vqrshrn_n_u16(sum0, FILTER_BITS);
            uint8x8_t d1 = vqrshrn_n_u16(sum1, FILTER_BITS);

            store_u8x4_strided_x2(dst + 0 * dst_stride, dst_stride, d0);
            store_u8x4_strided_x2(dst + 2 * dst_stride, dst_stride, d1);

            src += 4 * src_stride;
            dst += 4 * dst_stride;
            h -= 4;
        } while (h > 0);
    } else if (w == 8) {
        do {
            uint8x8_t s0 = vld1_u8(src + 0 * src_stride + 0);
            uint8x8_t s1 = vld1_u8(src + 0 * src_stride + 1);
            uint8x8_t s2 = vld1_u8(src + 1 * src_stride + 0);
            uint8x8_t s3 = vld1_u8(src + 1 * src_stride + 1);

            uint16x8_t sum0 = vmull_u8(s0, f0);
            sum0            = vmlal_u8(sum0, s1, f1);
            uint16x8_t sum1 = vmull_u8(s2, f0);
            sum1            = vmlal_u8(sum1, s3, f1);

            uint8x8_t d0 = vqrshrn_n_u16(sum0, FILTER_BITS);
            uint8x8_t d1 = vqrshrn_n_u16(sum1, FILTER_BITS);

            vst1_u8(dst + 0 * dst_stride, d0);
            vst1_u8(dst + 1 * dst_stride, d1);

            src += 2 * src_stride;
            dst += 2 * dst_stride;
            h -= 2;
        } while (h > 0);
    } else {
        do {
            int            width = w;
            const uint8_t *s     = src;
            uint8_t       *d     = dst;

            do {
                uint8x16_t s0 = vld1q_u8(s + 0);
                uint8x16_t s1 = vld1q_u8(s + 1);

                uint16x8_t sum0 = vmull_u8(vget_low_u8(s0), f0);
                sum0            = vmlal_u8(sum0, vget_low_u8(s1), f1);
                uint16x8_t sum1 = vmull_u8(vget_high_u8(s0), f0);
                sum1            = vmlal_u8(sum1, vget_high_u8(s1), f1);

                uint8x8_t d0 = vqrshrn_n_u16(sum0, FILTER_BITS);
                uint8x8_t d1 = vqrshrn_n_u16(sum1, FILTER_BITS);

                vst1q_u8(d, vcombine_u8(d0, d1));

                s += 16;
                d += 16;
                width -= 16;
            } while (width != 0);
            src += src_stride;
            dst += dst_stride;
        } while (--h > 0);
    }
}

static INLINE uint8x8_t convolve4_8(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
                                    const int16x4_t filter) {
    int16x8_t sum = vmulq_lane_s16(s0, filter, 0);
    sum           = vmlaq_lane_s16(sum, s1, filter, 1);
    sum           = vmlaq_lane_s16(sum, s2, filter, 2);
    sum           = vmlaq_lane_s16(sum, s3, filter, 3);

    // We halved the filter values so -1 from right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

#endif // AOM_CONVOLVE8_NEON_H_
