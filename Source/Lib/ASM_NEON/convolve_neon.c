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

#include "common_dsp_rtcd.h"
#include "definitions.h"
#include "sum_neon.h"
#include "mem_neon.h"
#include "transpose_neon.h"
#include "filter.h"
#include "inter_prediction.h"
#include "utility.h"

#define ROUND0_BITS 3

static INLINE int32x4_t convolve12_4_2d_v(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2,
                                          const int16x4_t s3, const int16x4_t s4, const int16x4_t s5,
                                          const int16x4_t s6, const int16x4_t s7, const int16x4_t s8,
                                          const int16x4_t s9, const int16x4_t s10, const int16x4_t s11,
                                          const int16x8_t y_filter_0_7, const int16x4_t y_filter_8_11) {
    const int16x4_t y_filter_0_3 = vget_low_s16(y_filter_0_7);
    const int16x4_t y_filter_4_7 = vget_high_s16(y_filter_0_7);

    int32x4_t sum = vmull_lane_s16(s0, y_filter_0_3, 0);
    sum           = vmlal_lane_s16(sum, s1, y_filter_0_3, 1);
    sum           = vmlal_lane_s16(sum, s2, y_filter_0_3, 2);
    sum           = vmlal_lane_s16(sum, s3, y_filter_0_3, 3);
    sum           = vmlal_lane_s16(sum, s4, y_filter_4_7, 0);
    sum           = vmlal_lane_s16(sum, s5, y_filter_4_7, 1);
    sum           = vmlal_lane_s16(sum, s6, y_filter_4_7, 2);
    sum           = vmlal_lane_s16(sum, s7, y_filter_4_7, 3);
    sum           = vmlal_lane_s16(sum, s8, y_filter_8_11, 0);
    sum           = vmlal_lane_s16(sum, s9, y_filter_8_11, 1);
    sum           = vmlal_lane_s16(sum, s10, y_filter_8_11, 2);
    sum           = vmlal_lane_s16(sum, s11, y_filter_8_11, 3);

    return sum;
}

static INLINE uint8x8_t convolve12_8_2d_v(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
                                          const int16x8_t s3, const int16x8_t s4, const int16x8_t s5,
                                          const int16x8_t s6, const int16x8_t s7, const int16x8_t s8,
                                          const int16x8_t s9, const int16x8_t s10, const int16x8_t s11,
                                          const int16x8_t y_filter_0_7, const int16x4_t y_filter_8_11,
                                          const int16x8_t sub_const) {
    const int16x4_t y_filter_0_3 = vget_low_s16(y_filter_0_7);
    const int16x4_t y_filter_4_7 = vget_high_s16(y_filter_0_7);

    int32x4_t sum0 = vmull_lane_s16(vget_low_s16(s0), y_filter_0_3, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s1), y_filter_0_3, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s2), y_filter_0_3, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s3), y_filter_0_3, 3);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s4), y_filter_4_7, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s5), y_filter_4_7, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s6), y_filter_4_7, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s7), y_filter_4_7, 3);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s8), y_filter_8_11, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s9), y_filter_8_11, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s10), y_filter_8_11, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s11), y_filter_8_11, 3);

    int32x4_t sum1 = vmull_lane_s16(vget_high_s16(s0), y_filter_0_3, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s1), y_filter_0_3, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s2), y_filter_0_3, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s3), y_filter_0_3, 3);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s4), y_filter_4_7, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s5), y_filter_4_7, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s6), y_filter_4_7, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s7), y_filter_4_7, 3);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s8), y_filter_8_11, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s9), y_filter_8_11, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s10), y_filter_8_11, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s11), y_filter_8_11, 3);

    int16x8_t res = vcombine_s16(vqrshrn_n_s32(sum0, 2 * FILTER_BITS - ROUND0_BITS),
                                 vqrshrn_n_s32(sum1, 2 * FILTER_BITS - ROUND0_BITS));
    res           = vsubq_s16(res, sub_const);

    return vqmovun_s16(res);
}

static INLINE void convolve_2d_sr_vert_12tap_neon(int16_t *src_ptr, int src_stride, uint8_t *dst_ptr, int dst_stride,
                                                  int w, int h, const int16x8_t y_filter_0_7,
                                                  const int16x4_t y_filter_8_11) {
    const int       bd        = 8;
    const int16x8_t sub_const = vdupq_n_s16(1 << (bd - 1));

    if (w <= 4) {
        int16x4_t s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10;
        load_s16_4x11(src_ptr, src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7, &s8, &s9, &s10);
        src_ptr += 11 * src_stride;

        do {
            int16x4_t s11, s12, s13, s14;
            load_s16_4x4(src_ptr, src_stride, &s11, &s12, &s13, &s14);

            int32x4_t d0 = convolve12_4_2d_v(
                s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, y_filter_0_7, y_filter_8_11);
            int32x4_t d1 = convolve12_4_2d_v(
                s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, y_filter_0_7, y_filter_8_11);
            int32x4_t d2 = convolve12_4_2d_v(
                s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, y_filter_0_7, y_filter_8_11);
            int32x4_t d3 = convolve12_4_2d_v(
                s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, y_filter_0_7, y_filter_8_11);

            int16x8_t dd01 = vcombine_s16(vqrshrn_n_s32(d0, 2 * FILTER_BITS - ROUND0_BITS),
                                          vqrshrn_n_s32(d1, 2 * FILTER_BITS - ROUND0_BITS));
            int16x8_t dd23 = vcombine_s16(vqrshrn_n_s32(d2, 2 * FILTER_BITS - ROUND0_BITS),
                                          vqrshrn_n_s32(d3, 2 * FILTER_BITS - ROUND0_BITS));

            dd01 = vsubq_s16(dd01, sub_const);
            dd23 = vsubq_s16(dd23, sub_const);

            uint8x8_t d01 = vqmovun_s16(dd01);
            uint8x8_t d23 = vqmovun_s16(dd23);

            store_u8_4x1(dst_ptr + 0 * dst_stride, d01, 0);
            store_u8_4x1(dst_ptr + 1 * dst_stride, d01, 1);
            store_u8_4x1(dst_ptr + 2 * dst_stride, d23, 0);
            store_u8_4x1(dst_ptr + 3 * dst_stride, d23, 1);

            s0  = s4;
            s1  = s5;
            s2  = s6;
            s3  = s7;
            s4  = s8;
            s5  = s9;
            s6  = s10;
            s7  = s11;
            s8  = s12;
            s9  = s13;
            s10 = s14;
            src_ptr += 4 * src_stride;
            dst_ptr += 4 * dst_stride;
            h -= 4;
        } while (h != 0);

    } else {
        do {
            int      height = h;
            int16_t *s      = src_ptr;
            uint8_t *d      = dst_ptr;

            int16x8_t s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10;
            load_s16_8x11(s, src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7, &s8, &s9, &s10);
            s += 11 * src_stride;

            do {
                int16x8_t s11, s12, s13, s14;
                load_s16_8x4(s, src_stride, &s11, &s12, &s13, &s14);

                uint8x8_t d0 = convolve12_8_2d_v(
                    s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, y_filter_0_7, y_filter_8_11, sub_const);
                uint8x8_t d1 = convolve12_8_2d_v(
                    s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, y_filter_0_7, y_filter_8_11, sub_const);
                uint8x8_t d2 = convolve12_8_2d_v(
                    s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, y_filter_0_7, y_filter_8_11, sub_const);
                uint8x8_t d3 = convolve12_8_2d_v(
                    s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, y_filter_0_7, y_filter_8_11, sub_const);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s0  = s4;
                s1  = s5;
                s2  = s6;
                s3  = s7;
                s4  = s8;
                s5  = s9;
                s6  = s10;
                s7  = s11;
                s8  = s12;
                s9  = s13;
                s10 = s14;
                s += 4 * src_stride;
                d += 4 * dst_stride;
                height -= 4;
            } while (height != 0);
            src_ptr += 8;
            dst_ptr += 8;
            w -= 8;
        } while (w != 0);
    }
}

static INLINE int16x4_t convolve8_4_2d_v(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
                                         const int16x4_t s4, const int16x4_t s5, const int16x4_t s6, const int16x4_t s7,
                                         const int16x8_t y_filter) {
    const int16x4_t y_filter_lo = vget_low_s16(y_filter);
    const int16x4_t y_filter_hi = vget_high_s16(y_filter);

    int32x4_t sum = vmull_lane_s16(s0, y_filter_lo, 0);
    sum           = vmlal_lane_s16(sum, s1, y_filter_lo, 1);
    sum           = vmlal_lane_s16(sum, s2, y_filter_lo, 2);
    sum           = vmlal_lane_s16(sum, s3, y_filter_lo, 3);
    sum           = vmlal_lane_s16(sum, s4, y_filter_hi, 0);
    sum           = vmlal_lane_s16(sum, s5, y_filter_hi, 1);
    sum           = vmlal_lane_s16(sum, s6, y_filter_hi, 2);
    sum           = vmlal_lane_s16(sum, s7, y_filter_hi, 3);

    return vqrshrn_n_s32(sum, 2 * FILTER_BITS - ROUND0_BITS);
}

static INLINE uint8x8_t convolve8_8_2d_v(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
                                         const int16x8_t s4, const int16x8_t s5, const int16x8_t s6, const int16x8_t s7,
                                         const int16x8_t y_filter, const int16x8_t sub_const) {
    const int16x4_t y_filter_lo = vget_low_s16(y_filter);
    const int16x4_t y_filter_hi = vget_high_s16(y_filter);

    int32x4_t sum0 = vmull_lane_s16(vget_low_s16(s0), y_filter_lo, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s1), y_filter_lo, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s2), y_filter_lo, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s3), y_filter_lo, 3);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s4), y_filter_hi, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s5), y_filter_hi, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s6), y_filter_hi, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s7), y_filter_hi, 3);

    int32x4_t sum1 = vmull_lane_s16(vget_high_s16(s0), y_filter_lo, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s1), y_filter_lo, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s2), y_filter_lo, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s3), y_filter_lo, 3);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s4), y_filter_hi, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s5), y_filter_hi, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s6), y_filter_hi, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s7), y_filter_hi, 3);

    int16x8_t res = vcombine_s16(vqrshrn_n_s32(sum0, 2 * FILTER_BITS - ROUND0_BITS),
                                 vqrshrn_n_s32(sum1, 2 * FILTER_BITS - ROUND0_BITS));
    res           = vsubq_s16(res, sub_const);

    return vqmovun_s16(res);
}

static INLINE void convolve_2d_sr_vert_8tap_neon(int16_t *src_ptr, int src_stride, uint8_t *dst_ptr, int dst_stride,
                                                 int w, int h, const int16x8_t y_filter) {
    const int       bd        = 8;
    const int16x8_t sub_const = vdupq_n_s16(1 << (bd - 1));

    if (w <= 4) {
        int16x4_t s0, s1, s2, s3, s4, s5, s6;
        load_s16_4x7(src_ptr, src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6);
        src_ptr += 7 * src_stride;

        do {
            int16x4_t s7, s8, s9, s10;
            load_s16_4x4(src_ptr, src_stride, &s7, &s8, &s9, &s10);

            int16x4_t d0 = convolve8_4_2d_v(s0, s1, s2, s3, s4, s5, s6, s7, y_filter);
            int16x4_t d1 = convolve8_4_2d_v(s1, s2, s3, s4, s5, s6, s7, s8, y_filter);
            int16x4_t d2 = convolve8_4_2d_v(s2, s3, s4, s5, s6, s7, s8, s9, y_filter);
            int16x4_t d3 = convolve8_4_2d_v(s3, s4, s5, s6, s7, s8, s9, s10, y_filter);

            uint8x8_t d01 = vqmovun_s16(vsubq_s16(vcombine_s16(d0, d1), sub_const));
            uint8x8_t d23 = vqmovun_s16(vsubq_s16(vcombine_s16(d2, d3), sub_const));

            store_u8_4x1(dst_ptr + 0 * dst_stride, d01, 0);
            store_u8_4x1(dst_ptr + 1 * dst_stride, d01, 1);
            store_u8_4x1(dst_ptr + 2 * dst_stride, d23, 0);
            store_u8_4x1(dst_ptr + 3 * dst_stride, d23, 1);

            s0 = s4;
            s1 = s5;
            s2 = s6;
            s3 = s7;
            s4 = s8;
            s5 = s9;
            s6 = s10;
            src_ptr += 4 * src_stride;
            dst_ptr += 4 * dst_stride;
            h -= 4;
        } while (h != 0);
    } else {
        // Width is a multiple of 8 and height is a multiple of 4.
        do {
            int      height = h;
            int16_t *s      = src_ptr;
            uint8_t *d      = dst_ptr;

            int16x8_t s0, s1, s2, s3, s4, s5, s6;
            load_s16_8x7(s, src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6);
            s += 7 * src_stride;

            do {
                int16x8_t s7, s8, s9, s10;
                load_s16_8x4(s, src_stride, &s7, &s8, &s9, &s10);

                uint8x8_t d0 = convolve8_8_2d_v(s0, s1, s2, s3, s4, s5, s6, s7, y_filter, sub_const);
                uint8x8_t d1 = convolve8_8_2d_v(s1, s2, s3, s4, s5, s6, s7, s8, y_filter, sub_const);
                uint8x8_t d2 = convolve8_8_2d_v(s2, s3, s4, s5, s6, s7, s8, s9, y_filter, sub_const);
                uint8x8_t d3 = convolve8_8_2d_v(s3, s4, s5, s6, s7, s8, s9, s10, y_filter, sub_const);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s0 = s4;
                s1 = s5;
                s2 = s6;
                s3 = s7;
                s4 = s8;
                s5 = s9;
                s6 = s10;
                s += 4 * src_stride;
                d += 4 * dst_stride;
                height -= 4;
            } while (height != 0);
            src_ptr += 8;
            dst_ptr += 8;
            w -= 8;
        } while (w != 0);
    }
}

static INLINE int16x4_t convolve6_4_2d_v(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
                                         const int16x4_t s4, const int16x4_t s5, const int16x8_t y_filter) {
    const int16x4_t y_filter_lo = vget_low_s16(y_filter);
    const int16x4_t y_filter_hi = vget_high_s16(y_filter);

    int32x4_t sum = vmull_lane_s16(s0, y_filter_lo, 1);
    sum           = vmlal_lane_s16(sum, s1, y_filter_lo, 2);
    sum           = vmlal_lane_s16(sum, s2, y_filter_lo, 3);
    sum           = vmlal_lane_s16(sum, s3, y_filter_hi, 0);
    sum           = vmlal_lane_s16(sum, s4, y_filter_hi, 1);
    sum           = vmlal_lane_s16(sum, s5, y_filter_hi, 2);

    return vqrshrn_n_s32(sum, 2 * FILTER_BITS - ROUND0_BITS);
}

static INLINE uint8x8_t convolve6_8_2d_v(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
                                         const int16x8_t s4, const int16x8_t s5, const int16x8_t y_filter,
                                         const int16x8_t sub_const) {
    const int16x4_t y_filter_lo = vget_low_s16(y_filter);
    const int16x4_t y_filter_hi = vget_high_s16(y_filter);

    int32x4_t sum0 = vmull_lane_s16(vget_low_s16(s0), y_filter_lo, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s1), y_filter_lo, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s2), y_filter_lo, 3);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s3), y_filter_hi, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s4), y_filter_hi, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s5), y_filter_hi, 2);

    int32x4_t sum1 = vmull_lane_s16(vget_high_s16(s0), y_filter_lo, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s1), y_filter_lo, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s2), y_filter_lo, 3);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s3), y_filter_hi, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s4), y_filter_hi, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s5), y_filter_hi, 2);

    int16x8_t res = vcombine_s16(vqrshrn_n_s32(sum0, 2 * FILTER_BITS - ROUND0_BITS),
                                 vqrshrn_n_s32(sum1, 2 * FILTER_BITS - ROUND0_BITS));
    res           = vsubq_s16(res, sub_const);

    return vqmovun_s16(res);
}

static INLINE void convolve_2d_sr_vert_6tap_neon(int16_t *src_ptr, int src_stride, uint8_t *dst_ptr, int dst_stride,
                                                 int w, int h, const int16x8_t y_filter) {
    const int       bd        = 8;
    const int16x8_t sub_const = vdupq_n_s16(1 << (bd - 1));

    if (w <= 4) {
        int16x4_t s0, s1, s2, s3, s4;
        load_s16_4x5(src_ptr, src_stride, &s0, &s1, &s2, &s3, &s4);
        src_ptr += 5 * src_stride;

        do {
            int16x4_t s5, s6, s7, s8;
            load_s16_4x4(src_ptr, src_stride, &s5, &s6, &s7, &s8);

            int16x4_t d0 = convolve6_4_2d_v(s0, s1, s2, s3, s4, s5, y_filter);
            int16x4_t d1 = convolve6_4_2d_v(s1, s2, s3, s4, s5, s6, y_filter);
            int16x4_t d2 = convolve6_4_2d_v(s2, s3, s4, s5, s6, s7, y_filter);
            int16x4_t d3 = convolve6_4_2d_v(s3, s4, s5, s6, s7, s8, y_filter);

            uint8x8_t d01 = vqmovun_s16(vsubq_s16(vcombine_s16(d0, d1), sub_const));
            uint8x8_t d23 = vqmovun_s16(vsubq_s16(vcombine_s16(d2, d3), sub_const));

            store_u8_4x1(dst_ptr + 0 * dst_stride, d01, 0);
            store_u8_4x1(dst_ptr + 1 * dst_stride, d01, 1);
            store_u8_4x1(dst_ptr + 2 * dst_stride, d23, 0);
            store_u8_4x1(dst_ptr + 3 * dst_stride, d23, 1);

            s0 = s4;
            s1 = s5;
            s2 = s6;
            s3 = s7;
            s4 = s8;
            src_ptr += 4 * src_stride;
            dst_ptr += 4 * dst_stride;
            h -= 4;
        } while (h != 0);
    } else {
        // Width is a multiple of 8 and height is a multiple of 4.
        do {
            int      height = h;
            int16_t *s      = src_ptr;
            uint8_t *d      = dst_ptr;

            int16x8_t s0, s1, s2, s3, s4;
            load_s16_8x5(s, src_stride, &s0, &s1, &s2, &s3, &s4);
            s += 5 * src_stride;

            do {
                int16x8_t s5, s6, s7, s8;
                load_s16_8x4(s, src_stride, &s5, &s6, &s7, &s8);

                uint8x8_t d0 = convolve6_8_2d_v(s0, s1, s2, s3, s4, s5, y_filter, sub_const);
                uint8x8_t d1 = convolve6_8_2d_v(s1, s2, s3, s4, s5, s6, y_filter, sub_const);
                uint8x8_t d2 = convolve6_8_2d_v(s2, s3, s4, s5, s6, s7, y_filter, sub_const);
                uint8x8_t d3 = convolve6_8_2d_v(s3, s4, s5, s6, s7, s8, y_filter, sub_const);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s0 = s4;
                s1 = s5;
                s2 = s6;
                s3 = s7;
                s4 = s8;
                s += 4 * src_stride;
                d += 4 * dst_stride;
                height -= 4;
            } while (height != 0);
            src_ptr += 8;
            dst_ptr += 8;
            w -= 8;
        } while (w != 0);
    }
}

static INLINE int16x4_t convolve12_4_x(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
                                       const int16x4_t s4, const int16x4_t s5, const int16x4_t s6, const int16x4_t s7,
                                       const int16x4_t s8, const int16x4_t s9, const int16x4_t s10, const int16x4_t s11,
                                       const int16x8_t x_filter_0_7, const int16x4_t x_filter_8_11,
                                       const int32x4_t horiz_const) {
    const int16x4_t x_filter_0_3 = vget_low_s16(x_filter_0_7);
    const int16x4_t x_filter_4_7 = vget_high_s16(x_filter_0_7);

    int32x4_t sum = horiz_const;
    sum           = vmlal_lane_s16(sum, s0, x_filter_0_3, 0);
    sum           = vmlal_lane_s16(sum, s1, x_filter_0_3, 1);
    sum           = vmlal_lane_s16(sum, s2, x_filter_0_3, 2);
    sum           = vmlal_lane_s16(sum, s3, x_filter_0_3, 3);
    sum           = vmlal_lane_s16(sum, s4, x_filter_4_7, 0);
    sum           = vmlal_lane_s16(sum, s5, x_filter_4_7, 1);
    sum           = vmlal_lane_s16(sum, s6, x_filter_4_7, 2);
    sum           = vmlal_lane_s16(sum, s7, x_filter_4_7, 3);
    sum           = vmlal_lane_s16(sum, s8, x_filter_8_11, 0);
    sum           = vmlal_lane_s16(sum, s9, x_filter_8_11, 1);
    sum           = vmlal_lane_s16(sum, s10, x_filter_8_11, 2);
    sum           = vmlal_lane_s16(sum, s11, x_filter_8_11, 3);

    return vqrshrn_n_s32(sum, FILTER_BITS);
}

static INLINE void convolve_x_sr_12tap_neon(const uint8_t *src_ptr, int src_stride, uint8_t *dst_ptr,
                                            const int dst_stride, int w, int h, const int16_t *x_filter_ptr) {
    const int16x8_t x_filter_0_7  = vld1q_s16(x_filter_ptr);
    const int16x4_t x_filter_8_11 = vld1_s16(x_filter_ptr + 8);

    // A shim of 1 << (ROUND0_BITS - 1) enables us to use a single rounding right
    // shift by FILTER_BITS - instead of a first rounding right shift by
    // ROUND0_BITS, followed by second rounding right shift by FILTER_BITS -
    // ROUND0_BITS.
    const int32x4_t horiz_const = vdupq_n_s32(1 << (ROUND0_BITS - 1));

    do {
        const uint8_t *s     = src_ptr;
        uint8_t       *d     = dst_ptr;
        int            width = w;

        uint8x8_t t0, t1, t2, t3;
        load_u8_8x4(s, src_stride, &t0, &t1, &t2, &t3);
        transpose_elems_inplace_u8_8x4(&t0, &t1, &t2, &t3);

        int16x4_t s0 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
        int16x4_t s1 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
        int16x4_t s2 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));
        int16x4_t s3 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t3)));
        int16x4_t s4 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
        int16x4_t s5 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
        int16x4_t s6 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));
        int16x4_t s7 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t3)));

        load_u8_8x4(s + 8, src_stride, &t0, &t1, &t2, &t3);
        transpose_elems_inplace_u8_8x4(&t0, &t1, &t2, &t3);

        int16x4_t s8  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
        int16x4_t s9  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
        int16x4_t s10 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));

        s += 11;

        do {
            load_u8_8x4(s, src_stride, &t0, &t1, &t2, &t3);
            transpose_elems_inplace_u8_8x4(&t0, &t1, &t2, &t3);

            int16x4_t s11 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
            int16x4_t s12 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
            int16x4_t s13 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));
            int16x4_t s14 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t3)));

            int16x4_t d0 = convolve12_4_x(
                s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, x_filter_0_7, x_filter_8_11, horiz_const);
            int16x4_t d1 = convolve12_4_x(
                s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, x_filter_0_7, x_filter_8_11, horiz_const);
            int16x4_t d2 = convolve12_4_x(
                s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, x_filter_0_7, x_filter_8_11, horiz_const);
            int16x4_t d3 = convolve12_4_x(
                s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, x_filter_0_7, x_filter_8_11, horiz_const);

            transpose_elems_inplace_s16_4x4(&d0, &d1, &d2, &d3);

            uint8x8_t d01 = vqmovun_s16(vcombine_s16(d0, d1));
            uint8x8_t d23 = vqmovun_s16(vcombine_s16(d2, d3));

            store_u8_4x1(d + 0 * dst_stride, d01, 0);
            store_u8_4x1(d + 1 * dst_stride, d01, 1);
            store_u8_4x1(d + 2 * dst_stride, d23, 0);
            store_u8_4x1(d + 3 * dst_stride, d23, 1);

            s0  = s4;
            s1  = s5;
            s2  = s6;
            s3  = s7;
            s4  = s8;
            s5  = s9;
            s6  = s10;
            s7  = s11;
            s8  = s12;
            s9  = s13;
            s10 = s14;
            s += 4;
            d += 4;
            width -= 4;
        } while (width != 0);
        src_ptr += 4 * src_stride;
        dst_ptr += 4 * dst_stride;
        h -= 4;
    } while (h != 0);
}

static INLINE uint8x8_t convolve4_4_x(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
                                      const int16x4_t filter, const int16x4_t horiz_const) {
    int16x4_t sum = horiz_const;
    sum           = vmla_lane_s16(sum, s0, filter, 0);
    sum           = vmla_lane_s16(sum, s1, filter, 1);
    sum           = vmla_lane_s16(sum, s2, filter, 2);
    sum           = vmla_lane_s16(sum, s3, filter, 3);

    // We halved the convolution filter values so - 1 from the right shift.
    return vqrshrun_n_s16(vcombine_s16(sum, vdup_n_s16(0)), FILTER_BITS - 1);
}

static INLINE uint8x8_t convolve8_8_x(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
                                      const int16x8_t s4, const int16x8_t s5, const int16x8_t s6, const int16x8_t s7,
                                      const int16x8_t filter, const int16x8_t horiz_const) {
    const int16x4_t filter_lo = vget_low_s16(filter);
    const int16x4_t filter_hi = vget_high_s16(filter);

    int16x8_t sum = horiz_const;
    sum           = vmlaq_lane_s16(sum, s0, filter_lo, 0);
    sum           = vmlaq_lane_s16(sum, s1, filter_lo, 1);
    sum           = vmlaq_lane_s16(sum, s2, filter_lo, 2);
    sum           = vmlaq_lane_s16(sum, s3, filter_lo, 3);
    sum           = vmlaq_lane_s16(sum, s4, filter_hi, 0);
    sum           = vmlaq_lane_s16(sum, s5, filter_hi, 1);
    sum           = vmlaq_lane_s16(sum, s6, filter_hi, 2);
    sum           = vmlaq_lane_s16(sum, s7, filter_hi, 3);

    // We halved the convolution filter values so - 1 from the right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

void svt_av1_convolve_x_sr_neon(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w,
                                int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y,
                                const int32_t subpel_x_qn, const int32_t subpel_y_qn, ConvolveParams *conv_params) {
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

    if (filter_params_x->taps > 8) {
        convolve_x_sr_12tap_neon(src, src_stride, dst, dst_stride, w, h, x_filter_ptr);
        return;
    }

    // This shim of 1 << ((ROUND0_BITS - 1) - 1) enables us to use a single
    // rounding right shift by FILTER_BITS - instead of a first rounding right
    // shift by ROUND0_BITS, followed by second rounding right shift by
    // FILTER_BITS - ROUND0_BITS.
    // The outermost -1 is needed because we will halve the filter values.
    const int16x8_t horiz_const = vdupq_n_s16(1 << ((ROUND0_BITS - 1) - 1));

    if (w <= 4) {
        // 4-tap filters are used for blocks having width <= 4.
        // Filter values are even, so halve to reduce intermediate precision reqs.
        const int16x4_t x_filter = vshr_n_s16(vld1_s16(x_filter_ptr + 2), 1);

        src += 2;

        do {
            uint8x8_t t0 = vld1_u8(src); // a0 a1 a2 a3 a4 a5 a6 a7
            int16x4_t s0 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
            int16x4_t s4 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));

            int16x4_t s1 = vext_s16(s0, s4, 1); // a1 a2 a3 a4
            int16x4_t s2 = vext_s16(s0, s4, 2); // a2 a3 a4 a5
            int16x4_t s3 = vext_s16(s0, s4, 3); // a3 a4 a5 a6

            uint8x8_t d0 = convolve4_4_x(s0, s1, s2, s3, x_filter, vget_low_s16(horiz_const));

            store_u8_4x1(dst, d0, 0);

            src += src_stride;
            dst += dst_stride;
        } while (--h != 0);
    } else {
        // Filter values are even so halve to reduce precision requirements.
        const int16x8_t x_filter = vshrq_n_s16(vld1q_s16(x_filter_ptr), 1);

        while (h >= 8) {
            uint8x8_t t0, t1, t2, t3, t4, t5, t6, t7;
            load_u8_8x8(src, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);

            transpose_elems_inplace_u8_8x8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
            int16x8_t s0 = vreinterpretq_s16_u16(vmovl_u8(t0));
            int16x8_t s1 = vreinterpretq_s16_u16(vmovl_u8(t1));
            int16x8_t s2 = vreinterpretq_s16_u16(vmovl_u8(t2));
            int16x8_t s3 = vreinterpretq_s16_u16(vmovl_u8(t3));
            int16x8_t s4 = vreinterpretq_s16_u16(vmovl_u8(t4));
            int16x8_t s5 = vreinterpretq_s16_u16(vmovl_u8(t5));
            int16x8_t s6 = vreinterpretq_s16_u16(vmovl_u8(t6));

            int            width = w;
            const uint8_t *s     = src + 7;
            uint8_t       *d     = dst;

            __builtin_prefetch(d + 0 * dst_stride);
            __builtin_prefetch(d + 1 * dst_stride);
            __builtin_prefetch(d + 2 * dst_stride);
            __builtin_prefetch(d + 3 * dst_stride);
            __builtin_prefetch(d + 4 * dst_stride);
            __builtin_prefetch(d + 5 * dst_stride);
            __builtin_prefetch(d + 6 * dst_stride);
            __builtin_prefetch(d + 7 * dst_stride);

            do {
                uint8x8_t t8, t9, t10, t11, t12, t13, t14;
                load_u8_8x8(s, src_stride, &t7, &t8, &t9, &t10, &t11, &t12, &t13, &t14);

                transpose_elems_inplace_u8_8x8(&t7, &t8, &t9, &t10, &t11, &t12, &t13, &t14);
                int16x8_t s7  = vreinterpretq_s16_u16(vmovl_u8(t7));
                int16x8_t s8  = vreinterpretq_s16_u16(vmovl_u8(t8));
                int16x8_t s9  = vreinterpretq_s16_u16(vmovl_u8(t9));
                int16x8_t s10 = vreinterpretq_s16_u16(vmovl_u8(t10));
                int16x8_t s11 = vreinterpretq_s16_u16(vmovl_u8(t11));
                int16x8_t s12 = vreinterpretq_s16_u16(vmovl_u8(t12));
                int16x8_t s13 = vreinterpretq_s16_u16(vmovl_u8(t13));
                int16x8_t s14 = vreinterpretq_s16_u16(vmovl_u8(t14));

                uint8x8_t d0 = convolve8_8_x(s0, s1, s2, s3, s4, s5, s6, s7, x_filter, horiz_const);
                uint8x8_t d1 = convolve8_8_x(s1, s2, s3, s4, s5, s6, s7, s8, x_filter, horiz_const);
                uint8x8_t d2 = convolve8_8_x(s2, s3, s4, s5, s6, s7, s8, s9, x_filter, horiz_const);
                uint8x8_t d3 = convolve8_8_x(s3, s4, s5, s6, s7, s8, s9, s10, x_filter, horiz_const);
                uint8x8_t d4 = convolve8_8_x(s4, s5, s6, s7, s8, s9, s10, s11, x_filter, horiz_const);
                uint8x8_t d5 = convolve8_8_x(s5, s6, s7, s8, s9, s10, s11, s12, x_filter, horiz_const);
                uint8x8_t d6 = convolve8_8_x(s6, s7, s8, s9, s10, s11, s12, s13, x_filter, horiz_const);
                uint8x8_t d7 = convolve8_8_x(s7, s8, s9, s10, s11, s12, s13, s14, x_filter, horiz_const);

                transpose_elems_inplace_u8_8x8(&d0, &d1, &d2, &d3, &d4, &d5, &d6, &d7);

                store_u8_8x8(d, dst_stride, d0, d1, d2, d3, d4, d5, d6, d7);

                s0 = s8;
                s1 = s9;
                s2 = s10;
                s3 = s11;
                s4 = s12;
                s5 = s13;
                s6 = s14;
                s += 8;
                d += 8;
                width -= 8;
            } while (width != 0);
            src += 8 * src_stride;
            dst += 8 * dst_stride;
            h -= 8;
        }
        while (h-- != 0) {
            uint8x8_t t0 = vld1_u8(src); // a0 a1 a2 a3 a4 a5 a6 a7
            int16x8_t s0 = vreinterpretq_s16_u16(vmovl_u8(t0));

            int            width = w;
            const uint8_t *s     = src + 8;
            uint8_t       *d     = dst;

            __builtin_prefetch(d);

            do {
                uint8x8_t t8 = vld1_u8(s); // a8 a9 a10 a11 a12 a13 a14 a15
                int16x8_t s8 = vreinterpretq_s16_u16(vmovl_u8(t8));

                int16x8_t s1 = vextq_s16(s0, s8, 1); // a1 a2 a3 a4 a5 a6 a7 a8
                int16x8_t s2 = vextq_s16(s0, s8, 2); // a2 a3 a4 a5 a6 a7 a8 a9
                int16x8_t s3 = vextq_s16(s0, s8, 3); // a3 a4 a5 a6 a7 a8 a9 a10
                int16x8_t s4 = vextq_s16(s0, s8, 4); // a4 a5 a6 a7 a8 a9 a10 a11
                int16x8_t s5 = vextq_s16(s0, s8, 5); // a5 a6 a7 a8 a9 a10 a11 a12
                int16x8_t s6 = vextq_s16(s0, s8, 6); // a6 a7 a8 a9 a10 a11 a12 a13
                int16x8_t s7 = vextq_s16(s0, s8, 7); // a7 a8 a9 a10 a11 a12 a13 a14

                uint8x8_t d0 = convolve8_8_x(s0, s1, s2, s3, s4, s5, s6, s7, x_filter, horiz_const);

                vst1_u8(d, d0);

                s0 = s8;
                s += 8;
                d += 8;
                width -= 8;
            } while (width != 0);
            src += src_stride;
            dst += dst_stride;
        }
    }
}

static INLINE int16x4_t convolve6_4_y(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
                                      const int16x4_t s4, const int16x4_t s5, const int16x8_t y_filter_0_7) {
    const int16x4_t y_filter_0_3 = vget_low_s16(y_filter_0_7);
    const int16x4_t y_filter_4_7 = vget_high_s16(y_filter_0_7);

    // Filter values at indices 0 and 7 are 0.
    int16x4_t sum = vmul_lane_s16(s0, y_filter_0_3, 1);
    sum           = vmla_lane_s16(sum, s1, y_filter_0_3, 2);
    sum           = vmla_lane_s16(sum, s2, y_filter_0_3, 3);
    sum           = vmla_lane_s16(sum, s3, y_filter_4_7, 0);
    sum           = vmla_lane_s16(sum, s4, y_filter_4_7, 1);
    sum           = vmla_lane_s16(sum, s5, y_filter_4_7, 2);

    return sum;
}

static INLINE uint8x8_t convolve6_8_y(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
                                      const int16x8_t s4, const int16x8_t s5, const int16x8_t y_filters) {
    const int16x4_t y_filter_lo = vget_low_s16(y_filters);
    const int16x4_t y_filter_hi = vget_high_s16(y_filters);

    // Filter values at indices 0 and 7 are 0.
    int16x8_t sum = vmulq_lane_s16(s0, y_filter_lo, 1);
    sum           = vmlaq_lane_s16(sum, s1, y_filter_lo, 2);
    sum           = vmlaq_lane_s16(sum, s2, y_filter_lo, 3);
    sum           = vmlaq_lane_s16(sum, s3, y_filter_hi, 0);
    sum           = vmlaq_lane_s16(sum, s4, y_filter_hi, 1);
    sum           = vmlaq_lane_s16(sum, s5, y_filter_hi, 2);
    // We halved the convolution filter values so -1 from the right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

static INLINE void convolve_y_sr_6tap_neon(const uint8_t *src_ptr, int src_stride, uint8_t *dst_ptr,
                                           const int dst_stride, int w, int h, const int16x8_t y_filter) {
    if (w <= 4) {
        uint8x8_t t0 = load_unaligned_u8_4x1(src_ptr + 0 * src_stride);
        uint8x8_t t1 = load_unaligned_u8_4x1(src_ptr + 1 * src_stride);
        uint8x8_t t2 = load_unaligned_u8_4x1(src_ptr + 2 * src_stride);
        uint8x8_t t3 = load_unaligned_u8_4x1(src_ptr + 3 * src_stride);
        uint8x8_t t4 = load_unaligned_u8_4x1(src_ptr + 4 * src_stride);

        int16x4_t s0 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
        int16x4_t s1 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
        int16x4_t s2 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));
        int16x4_t s3 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t3)));
        int16x4_t s4 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t4)));

        src_ptr += 5 * src_stride;

        do {
            uint8x8_t t5 = load_unaligned_u8_4x1(src_ptr + 0 * src_stride);
            uint8x8_t t6 = load_unaligned_u8_4x1(src_ptr + 1 * src_stride);
            uint8x8_t t7 = load_unaligned_u8_4x1(src_ptr + 2 * src_stride);
            uint8x8_t t8 = load_unaligned_u8_4x1(src_ptr + 3 * src_stride);

            int16x4_t s5 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t5)));
            int16x4_t s6 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t6)));
            int16x4_t s7 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t7)));
            int16x4_t s8 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t8)));

            int16x4_t d0 = convolve6_4_y(s0, s1, s2, s3, s4, s5, y_filter);
            int16x4_t d1 = convolve6_4_y(s1, s2, s3, s4, s5, s6, y_filter);
            int16x4_t d2 = convolve6_4_y(s2, s3, s4, s5, s6, s7, y_filter);
            int16x4_t d3 = convolve6_4_y(s3, s4, s5, s6, s7, s8, y_filter);

            // We halved the convolution filter values so -1 from the right shift.
            uint8x8_t d01 = vqrshrun_n_s16(vcombine_s16(d0, d1), FILTER_BITS - 1);
            uint8x8_t d23 = vqrshrun_n_s16(vcombine_s16(d2, d3), FILTER_BITS - 1);

            store_u8_4x1(dst_ptr + 0 * dst_stride, d01, 0);
            store_u8_4x1(dst_ptr + 1 * dst_stride, d01, 1);
            store_u8_4x1(dst_ptr + 2 * dst_stride, d23, 0);
            store_u8_4x1(dst_ptr + 3 * dst_stride, d23, 1);

            s0 = s4;
            s1 = s5;
            s2 = s6;
            s3 = s7;
            s4 = s8;
            src_ptr += 4 * src_stride;
            dst_ptr += 4 * dst_stride;
            h -= 4;
        } while (h != 0);

    } else {
        do {
            const uint8_t *s      = src_ptr;
            uint8_t       *d      = dst_ptr;
            int            height = h;

            uint8x8_t t0, t1, t2, t3, t4;
            load_u8_8x5(s, src_stride, &t0, &t1, &t2, &t3, &t4);

            int16x8_t s0 = vreinterpretq_s16_u16(vmovl_u8(t0));
            int16x8_t s1 = vreinterpretq_s16_u16(vmovl_u8(t1));
            int16x8_t s2 = vreinterpretq_s16_u16(vmovl_u8(t2));
            int16x8_t s3 = vreinterpretq_s16_u16(vmovl_u8(t3));
            int16x8_t s4 = vreinterpretq_s16_u16(vmovl_u8(t4));

            s += 5 * src_stride;

            do {
                uint8x8_t t5, t6, t7, t8;
                load_u8_8x4(s, src_stride, &t5, &t6, &t7, &t8);

                int16x8_t s5 = vreinterpretq_s16_u16(vmovl_u8(t5));
                int16x8_t s6 = vreinterpretq_s16_u16(vmovl_u8(t6));
                int16x8_t s7 = vreinterpretq_s16_u16(vmovl_u8(t7));
                int16x8_t s8 = vreinterpretq_s16_u16(vmovl_u8(t8));

                uint8x8_t d0 = convolve6_8_y(s0, s1, s2, s3, s4, s5, y_filter);
                uint8x8_t d1 = convolve6_8_y(s1, s2, s3, s4, s5, s6, y_filter);
                uint8x8_t d2 = convolve6_8_y(s2, s3, s4, s5, s6, s7, y_filter);
                uint8x8_t d3 = convolve6_8_y(s3, s4, s5, s6, s7, s8, y_filter);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s0 = s4;
                s1 = s5;
                s2 = s6;
                s3 = s7;
                s4 = s8;
                s += 4 * src_stride;
                d += 4 * dst_stride;
                height -= 4;
            } while (height != 0);
            src_ptr += 8;
            dst_ptr += 8;
            w -= 8;
        } while (w != 0);
    }
}

static INLINE int16x4_t convolve8_4_y(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
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

static INLINE uint8x8_t convolve8_8_y(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
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

    // We halved the convolution filter values so -1 from the right shift.
    return vqrshrun_n_s16(sum, FILTER_BITS - 1);
}

static INLINE void convolve_y_sr_8tap_neon(const uint8_t *src_ptr, int src_stride, uint8_t *dst_ptr,
                                           const int dst_stride, int w, int h, const int16x8_t y_filter) {
    if (w <= 4) {
        uint8x8_t t0 = load_unaligned_u8_4x1(src_ptr + 0 * src_stride);
        uint8x8_t t1 = load_unaligned_u8_4x1(src_ptr + 1 * src_stride);
        uint8x8_t t2 = load_unaligned_u8_4x1(src_ptr + 2 * src_stride);
        uint8x8_t t3 = load_unaligned_u8_4x1(src_ptr + 3 * src_stride);
        uint8x8_t t4 = load_unaligned_u8_4x1(src_ptr + 4 * src_stride);
        uint8x8_t t5 = load_unaligned_u8_4x1(src_ptr + 5 * src_stride);
        uint8x8_t t6 = load_unaligned_u8_4x1(src_ptr + 6 * src_stride);

        int16x4_t s0 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t0)));
        int16x4_t s1 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t1)));
        int16x4_t s2 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t2)));
        int16x4_t s3 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t3)));
        int16x4_t s4 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t4)));
        int16x4_t s5 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t5)));
        int16x4_t s6 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t6)));

        src_ptr += 7 * src_stride;

        do {
            uint8x8_t t7  = load_unaligned_u8_4x1(src_ptr + 0 * src_stride);
            uint8x8_t t8  = load_unaligned_u8_4x1(src_ptr + 1 * src_stride);
            uint8x8_t t9  = load_unaligned_u8_4x1(src_ptr + 2 * src_stride);
            uint8x8_t t10 = load_unaligned_u8_4x1(src_ptr + 3 * src_stride);

            int16x4_t s7  = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t7)));
            int16x4_t s8  = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t8)));
            int16x4_t s9  = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t9)));
            int16x4_t s10 = vreinterpret_s16_u16(vget_low_u16(vmovl_u8(t10)));

            int16x4_t d0 = convolve8_4_y(s0, s1, s2, s3, s4, s5, s6, s7, y_filter);
            int16x4_t d1 = convolve8_4_y(s1, s2, s3, s4, s5, s6, s7, s8, y_filter);
            int16x4_t d2 = convolve8_4_y(s2, s3, s4, s5, s6, s7, s8, s9, y_filter);
            int16x4_t d3 = convolve8_4_y(s3, s4, s5, s6, s7, s8, s9, s10, y_filter);

            // We halved the convolution filter values so -1 from the right shift.
            uint8x8_t d01 = vqrshrun_n_s16(vcombine_s16(d0, d1), FILTER_BITS - 1);
            uint8x8_t d23 = vqrshrun_n_s16(vcombine_s16(d2, d3), FILTER_BITS - 1);

            store_u8_4x1(dst_ptr + 0 * dst_stride, d01, 0);
            store_u8_4x1(dst_ptr + 1 * dst_stride, d01, 1);
            store_u8_4x1(dst_ptr + 2 * dst_stride, d23, 0);
            store_u8_4x1(dst_ptr + 3 * dst_stride, d23, 1);

            s0 = s4;
            s1 = s5;
            s2 = s6;
            s3 = s7;
            s4 = s8;
            s5 = s9;
            s6 = s10;
            src_ptr += 4 * src_stride;
            dst_ptr += 4 * dst_stride;
            h -= 4;
        } while (h != 0);
    } else {
        do {
            const uint8_t *s      = src_ptr;
            uint8_t       *d      = dst_ptr;
            int            height = h;

            uint8x8_t t0, t1, t2, t3, t4, t5, t6;
            load_u8_8x7(s, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6);

            int16x8_t s0 = vreinterpretq_s16_u16(vmovl_u8(t0));
            int16x8_t s1 = vreinterpretq_s16_u16(vmovl_u8(t1));
            int16x8_t s2 = vreinterpretq_s16_u16(vmovl_u8(t2));
            int16x8_t s3 = vreinterpretq_s16_u16(vmovl_u8(t3));
            int16x8_t s4 = vreinterpretq_s16_u16(vmovl_u8(t4));
            int16x8_t s5 = vreinterpretq_s16_u16(vmovl_u8(t5));
            int16x8_t s6 = vreinterpretq_s16_u16(vmovl_u8(t6));

            s += 7 * src_stride;

            do {
                uint8x8_t t7, t8, t9, t10;
                load_u8_8x4(s, src_stride, &t7, &t8, &t9, &t10);

                int16x8_t s7  = vreinterpretq_s16_u16(vmovl_u8(t7));
                int16x8_t s8  = vreinterpretq_s16_u16(vmovl_u8(t8));
                int16x8_t s9  = vreinterpretq_s16_u16(vmovl_u8(t9));
                int16x8_t s10 = vreinterpretq_s16_u16(vmovl_u8(t10));

                uint8x8_t d0 = convolve8_8_y(s0, s1, s2, s3, s4, s5, s6, s7, y_filter);
                uint8x8_t d1 = convolve8_8_y(s1, s2, s3, s4, s5, s6, s7, s8, y_filter);
                uint8x8_t d2 = convolve8_8_y(s2, s3, s4, s5, s6, s7, s8, s9, y_filter);
                uint8x8_t d3 = convolve8_8_y(s3, s4, s5, s6, s7, s8, s9, s10, y_filter);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s0 = s4;
                s1 = s5;
                s2 = s6;
                s3 = s7;
                s4 = s8;
                s5 = s9;
                s6 = s10;
                s += 4 * src_stride;
                d += 4 * dst_stride;
                height -= 4;
            } while (height != 0);
            src_ptr += 8;
            dst_ptr += 8;
            w -= 8;
        } while (w != 0);
    }
}

static INLINE int16x4_t convolve12_4_y(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
                                       const int16x4_t s4, const int16x4_t s5, const int16x4_t s6, const int16x4_t s7,
                                       const int16x4_t s8, const int16x4_t s9, const int16x4_t s10, const int16x4_t s11,
                                       const int16x8_t y_filter_0_7, const int16x4_t y_filter_8_11) {
    const int16x4_t y_filter_0_3 = vget_low_s16(y_filter_0_7);
    const int16x4_t y_filter_4_7 = vget_high_s16(y_filter_0_7);
    int16x4_t       sum;

    sum = vmul_lane_s16(s0, y_filter_0_3, 0);
    sum = vmla_lane_s16(sum, s1, y_filter_0_3, 1);
    sum = vmla_lane_s16(sum, s2, y_filter_0_3, 2);
    sum = vmla_lane_s16(sum, s3, y_filter_0_3, 3);
    sum = vmla_lane_s16(sum, s4, y_filter_4_7, 0);

    sum = vmla_lane_s16(sum, s7, y_filter_4_7, 3);
    sum = vmla_lane_s16(sum, s8, y_filter_8_11, 0);
    sum = vmla_lane_s16(sum, s9, y_filter_8_11, 1);
    sum = vmla_lane_s16(sum, s10, y_filter_8_11, 2);
    sum = vmla_lane_s16(sum, s11, y_filter_8_11, 3);

    // Saturating addition is required for the largest filter taps to avoid
    // overflow (while staying in 16-bit elements.)
    sum = vqadd_s16(sum, vmul_lane_s16(s5, y_filter_4_7, 1));
    sum = vqadd_s16(sum, vmul_lane_s16(s6, y_filter_4_7, 2));

    return sum;
}

static INLINE uint8x8_t convolve12_8_y(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
                                       const int16x8_t s4, const int16x8_t s5, const int16x8_t s6, const int16x8_t s7,
                                       const int16x8_t s8, const int16x8_t s9, const int16x8_t s10, const int16x8_t s11,
                                       const int16x8_t y_filter_0_7, const int16x4_t y_filter_8_11) {
    const int16x4_t y_filter_0_3 = vget_low_s16(y_filter_0_7);
    const int16x4_t y_filter_4_7 = vget_high_s16(y_filter_0_7);
    int16x8_t       sum;

    sum = vmulq_lane_s16(s0, y_filter_0_3, 0);
    sum = vmlaq_lane_s16(sum, s1, y_filter_0_3, 1);
    sum = vmlaq_lane_s16(sum, s2, y_filter_0_3, 2);
    sum = vmlaq_lane_s16(sum, s3, y_filter_0_3, 3);
    sum = vmlaq_lane_s16(sum, s4, y_filter_4_7, 0);

    sum = vmlaq_lane_s16(sum, s7, y_filter_4_7, 3);
    sum = vmlaq_lane_s16(sum, s8, y_filter_8_11, 0);
    sum = vmlaq_lane_s16(sum, s9, y_filter_8_11, 1);
    sum = vmlaq_lane_s16(sum, s10, y_filter_8_11, 2);
    sum = vmlaq_lane_s16(sum, s11, y_filter_8_11, 3);

    // Saturating addition is required for the largest filter taps to avoid
    // overflow (while staying in 16-bit elements.)
    sum = vqaddq_s16(sum, vmulq_lane_s16(s5, y_filter_4_7, 1));
    sum = vqaddq_s16(sum, vmulq_lane_s16(s6, y_filter_4_7, 2));

    return vqrshrun_n_s16(sum, FILTER_BITS);
}

static INLINE void convolve_y_sr_12tap_neon(const uint8_t *src_ptr, int src_stride, uint8_t *dst_ptr, int dst_stride,
                                            int w, int h, const int16_t *y_filter_ptr) {
    const int16x8_t y_filter_0_7  = vld1q_s16(y_filter_ptr);
    const int16x4_t y_filter_8_11 = vld1_s16(y_filter_ptr + 8);

    if (w <= 4) {
        uint8x8_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        load_u8_8x11(src_ptr, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, &t10);
        int16x4_t s0  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
        int16x4_t s1  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
        int16x4_t s2  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));
        int16x4_t s3  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t3)));
        int16x4_t s4  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t4)));
        int16x4_t s5  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t5)));
        int16x4_t s6  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t6)));
        int16x4_t s7  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t7)));
        int16x4_t s8  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t8)));
        int16x4_t s9  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t9)));
        int16x4_t s10 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t10)));

        src_ptr += 11 * src_stride;

        do {
            uint8x8_t t11, t12, t13, t14;
            load_u8_8x4(src_ptr, src_stride, &t11, &t12, &t13, &t14);

            int16x4_t s11 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t11)));
            int16x4_t s12 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t12)));
            int16x4_t s13 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t13)));
            int16x4_t s14 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t14)));

            int16x4_t d0 = convolve12_4_y(
                s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, y_filter_0_7, y_filter_8_11);
            int16x4_t d1 = convolve12_4_y(
                s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, y_filter_0_7, y_filter_8_11);
            int16x4_t d2 = convolve12_4_y(
                s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, y_filter_0_7, y_filter_8_11);
            int16x4_t d3 = convolve12_4_y(
                s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, y_filter_0_7, y_filter_8_11);

            uint8x8_t d01 = vqrshrun_n_s16(vcombine_s16(d0, d1), FILTER_BITS);
            uint8x8_t d23 = vqrshrun_n_s16(vcombine_s16(d2, d3), FILTER_BITS);

            store_u8_4x1(dst_ptr + 0 * dst_stride, d01, 0);
            store_u8_4x1(dst_ptr + 1 * dst_stride, d01, 1);
            store_u8_4x1(dst_ptr + 2 * dst_stride, d23, 0);
            store_u8_4x1(dst_ptr + 3 * dst_stride, d23, 1);

            s0  = s4;
            s1  = s5;
            s2  = s6;
            s3  = s7;
            s4  = s8;
            s5  = s9;
            s6  = s10;
            s7  = s11;
            s8  = s12;
            s9  = s13;
            s10 = s14;
            src_ptr += 4 * src_stride;
            dst_ptr += 4 * dst_stride;
            h -= 4;
        } while (h != 0);

    } else {
        do {
            const uint8_t *s      = src_ptr;
            uint8_t       *d      = dst_ptr;
            int            height = h;

            uint8x8_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
            load_u8_8x11(s, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, &t10);
            int16x8_t s0  = vreinterpretq_s16_u16(vmovl_u8(t0));
            int16x8_t s1  = vreinterpretq_s16_u16(vmovl_u8(t1));
            int16x8_t s2  = vreinterpretq_s16_u16(vmovl_u8(t2));
            int16x8_t s3  = vreinterpretq_s16_u16(vmovl_u8(t3));
            int16x8_t s4  = vreinterpretq_s16_u16(vmovl_u8(t4));
            int16x8_t s5  = vreinterpretq_s16_u16(vmovl_u8(t5));
            int16x8_t s6  = vreinterpretq_s16_u16(vmovl_u8(t6));
            int16x8_t s7  = vreinterpretq_s16_u16(vmovl_u8(t7));
            int16x8_t s8  = vreinterpretq_s16_u16(vmovl_u8(t8));
            int16x8_t s9  = vreinterpretq_s16_u16(vmovl_u8(t9));
            int16x8_t s10 = vreinterpretq_s16_u16(vmovl_u8(t10));

            s += 11 * src_stride;

            do {
                uint8x8_t t11, t12, t13, t14;
                load_u8_8x4(s, src_stride, &t11, &t12, &t13, &t14);

                int16x8_t s11 = vreinterpretq_s16_u16(vmovl_u8(t11));
                int16x8_t s12 = vreinterpretq_s16_u16(vmovl_u8(t12));
                int16x8_t s13 = vreinterpretq_s16_u16(vmovl_u8(t13));
                int16x8_t s14 = vreinterpretq_s16_u16(vmovl_u8(t14));

                uint8x8_t d0 = convolve12_8_y(
                    s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, y_filter_0_7, y_filter_8_11);
                uint8x8_t d1 = convolve12_8_y(
                    s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, y_filter_0_7, y_filter_8_11);
                uint8x8_t d2 = convolve12_8_y(
                    s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, y_filter_0_7, y_filter_8_11);
                uint8x8_t d3 = convolve12_8_y(
                    s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, y_filter_0_7, y_filter_8_11);

                store_u8_8x4(d, dst_stride, d0, d1, d2, d3);

                s0  = s4;
                s1  = s5;
                s2  = s6;
                s3  = s7;
                s4  = s8;
                s5  = s9;
                s6  = s10;
                s7  = s11;
                s8  = s12;
                s9  = s13;
                s10 = s14;
                s += 4 * src_stride;
                d += 4 * dst_stride;
                height -= 4;
            } while (height != 0);
            src_ptr += 8;
            dst_ptr += 8;
            w -= 8;
        } while (w != 0);
    }
}

void svt_av1_convolve_y_sr_neon(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w,
                                int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y,
                                const int32_t subpel_x_qn, const int32_t subpel_y_qn, ConvolveParams *conv_params) {
    if (w == 2 || h == 2) {
        svt_av1_convolve_y_sr_c(src,
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

    const int y_filter_taps  = get_filter_tap(filter_params_y, subpel_y_qn);
    const int clamped_y_taps = y_filter_taps < 6 ? 6 : y_filter_taps;
    const int vert_offset    = clamped_y_taps / 2 - 1;

    src -= vert_offset * src_stride;

    const int16_t *y_filter_ptr = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_qn & SUBPEL_MASK);

    if (y_filter_taps > 8) {
        convolve_y_sr_12tap_neon(src, src_stride, dst, dst_stride, w, h, y_filter_ptr);
        return;
    }

    // Filter values are even so halve to reduce precision requirements.
    const int16x8_t y_filter = vshrq_n_s16(vld1q_s16(y_filter_ptr), 1);

    if (y_filter_taps < 8) {
        convolve_y_sr_6tap_neon(src, src_stride, dst, dst_stride, w, h, y_filter);
    } else {
        convolve_y_sr_8tap_neon(src, src_stride, dst, dst_stride, w, h, y_filter);
    }
}

static INLINE int16x4_t convolve12_4_2d_h(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2,
                                          const int16x4_t s3, const int16x4_t s4, const int16x4_t s5,
                                          const int16x4_t s6, const int16x4_t s7, const int16x4_t s8,
                                          const int16x4_t s9, const int16x4_t s10, const int16x4_t s11,
                                          const int16x8_t x_filter_0_7, const int16x4_t x_filter_8_11,
                                          const int32x4_t horiz_const) {
    const int16x4_t x_filter_0_3 = vget_low_s16(x_filter_0_7);
    const int16x4_t x_filter_4_7 = vget_high_s16(x_filter_0_7);

    int32x4_t sum = horiz_const;
    sum           = vmlal_lane_s16(sum, s0, x_filter_0_3, 0);
    sum           = vmlal_lane_s16(sum, s1, x_filter_0_3, 1);
    sum           = vmlal_lane_s16(sum, s2, x_filter_0_3, 2);
    sum           = vmlal_lane_s16(sum, s3, x_filter_0_3, 3);
    sum           = vmlal_lane_s16(sum, s4, x_filter_4_7, 0);
    sum           = vmlal_lane_s16(sum, s5, x_filter_4_7, 1);
    sum           = vmlal_lane_s16(sum, s6, x_filter_4_7, 2);
    sum           = vmlal_lane_s16(sum, s7, x_filter_4_7, 3);
    sum           = vmlal_lane_s16(sum, s8, x_filter_8_11, 0);
    sum           = vmlal_lane_s16(sum, s9, x_filter_8_11, 1);
    sum           = vmlal_lane_s16(sum, s10, x_filter_8_11, 2);
    sum           = vmlal_lane_s16(sum, s11, x_filter_8_11, 3);

    return vshrn_n_s32(sum, ROUND0_BITS);
}

static INLINE void convolve_2d_sr_horiz_12tap_neon(const uint8_t *src_ptr, int src_stride, int16_t *dst_ptr,
                                                   const int dst_stride, int w, int h, const int16x8_t x_filter_0_7,
                                                   const int16x4_t x_filter_8_11) {
    const int bd = 8;
    // A shim of 1 << (ROUND0_BITS - 1) enables us to use non-rounding shifts -
    // which are generally faster than rounding shifts on modern CPUs.
    const int32x4_t horiz_const = vdupq_n_s32((1 << (bd + FILTER_BITS - 1)) + (1 << (ROUND0_BITS - 1)));

    do {
        const uint8_t *s     = src_ptr;
        int16_t       *d     = dst_ptr;
        int            width = w;

        uint8x8_t t0, t1, t2, t3;
        load_u8_8x4(s, src_stride, &t0, &t1, &t2, &t3);
        transpose_elems_inplace_u8_8x4(&t0, &t1, &t2, &t3);

        int16x4_t s0 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
        int16x4_t s1 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
        int16x4_t s2 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));
        int16x4_t s3 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t3)));
        int16x4_t s4 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
        int16x4_t s5 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
        int16x4_t s6 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));
        int16x4_t s7 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t3)));

        load_u8_8x4(s + 8, src_stride, &t0, &t1, &t2, &t3);
        transpose_elems_inplace_u8_8x4(&t0, &t1, &t2, &t3);

        int16x4_t s8  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
        int16x4_t s9  = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
        int16x4_t s10 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));

        s += 11;

        do {
            load_u8_8x4(s, src_stride, &t0, &t1, &t2, &t3);
            transpose_elems_inplace_u8_8x4(&t0, &t1, &t2, &t3);

            int16x4_t s11 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
            int16x4_t s12 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t1)));
            int16x4_t s13 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t2)));
            int16x4_t s14 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t3)));

            int16x4_t d0 = convolve12_4_2d_h(
                s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, x_filter_0_7, x_filter_8_11, horiz_const);
            int16x4_t d1 = convolve12_4_2d_h(
                s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, x_filter_0_7, x_filter_8_11, horiz_const);
            int16x4_t d2 = convolve12_4_2d_h(
                s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, x_filter_0_7, x_filter_8_11, horiz_const);
            int16x4_t d3 = convolve12_4_2d_h(
                s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, x_filter_0_7, x_filter_8_11, horiz_const);

            transpose_elems_inplace_s16_4x4(&d0, &d1, &d2, &d3);
            store_s16_4x4(d, dst_stride, d0, d1, d2, d3);

            s0  = s4;
            s1  = s5;
            s2  = s6;
            s3  = s7;
            s4  = s8;
            s5  = s9;
            s6  = s10;
            s7  = s11;
            s8  = s12;
            s9  = s13;
            s10 = s14;
            s += 4;
            d += 4;
            width -= 4;
        } while (width != 0);
        src_ptr += 4 * src_stride;
        dst_ptr += 4 * dst_stride;
        h -= 4;
    } while (h > 4);

    do {
        const uint8_t *s     = src_ptr;
        int16_t       *d     = dst_ptr;
        int            width = w;

        do {
            uint8x16_t t0  = vld1q_u8(s);
            int16x8_t  tt0 = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(t0)));
            int16x8_t  tt1 = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(t0)));

            int16x4_t s0  = vget_low_s16(tt0);
            int16x4_t s4  = vget_high_s16(tt0);
            int16x4_t s8  = vget_low_s16(tt1);
            int16x4_t s12 = vget_high_s16(tt1);

            int16x4_t s1  = vext_s16(s0, s4, 1); //  a1  a2  a3  a4
            int16x4_t s2  = vext_s16(s0, s4, 2); //  a2  a3  a4  a5
            int16x4_t s3  = vext_s16(s0, s4, 3); //  a3  a4  a5  a6
            int16x4_t s5  = vext_s16(s4, s8, 1); //  a5  a6  a7  a8
            int16x4_t s6  = vext_s16(s4, s8, 2); //  a6  a7  a8  a9
            int16x4_t s7  = vext_s16(s4, s8, 3); //  a7  a8  a9 a10
            int16x4_t s9  = vext_s16(s8, s12, 1); //  a9 a10 a11 a12
            int16x4_t s10 = vext_s16(s8, s12, 2); // a10 a11 a12 a13
            int16x4_t s11 = vext_s16(s8, s12, 3); // a11 a12 a13 a14

            int16x4_t d0 = convolve12_4_2d_h(
                s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, x_filter_0_7, x_filter_8_11, horiz_const);
            vst1_s16(d, d0);

            s += 4;
            d += 4;
            width -= 4;
        } while (width != 0);
        src_ptr += src_stride;
        dst_ptr += dst_stride;
    } while (--h != 0);
}

static INLINE int16x4_t convolve4_4_2d_h(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2, const int16x4_t s3,
                                         const int16x4_t filter, const int16x4_t horiz_const) {
    int16x4_t sum = horiz_const;
    sum           = vmla_lane_s16(sum, s0, filter, 0);
    sum           = vmla_lane_s16(sum, s1, filter, 1);
    sum           = vmla_lane_s16(sum, s2, filter, 2);
    sum           = vmla_lane_s16(sum, s3, filter, 3);

    // We halved the convolution filter values so -1 from the right shift.
    return vshr_n_s16(sum, ROUND0_BITS - 1);
}

static INLINE int16x8_t convolve8_8_2d_h(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2, const int16x8_t s3,
                                         const int16x8_t s4, const int16x8_t s5, const int16x8_t s6, const int16x8_t s7,
                                         const int16x8_t filter, const int16x8_t horiz_const) {
    const int16x4_t filter_lo = vget_low_s16(filter);
    const int16x4_t filter_hi = vget_high_s16(filter);

    int16x8_t sum = horiz_const;
    sum           = vmlaq_lane_s16(sum, s0, filter_lo, 0);
    sum           = vmlaq_lane_s16(sum, s1, filter_lo, 1);
    sum           = vmlaq_lane_s16(sum, s2, filter_lo, 2);
    sum           = vmlaq_lane_s16(sum, s3, filter_lo, 3);
    sum           = vmlaq_lane_s16(sum, s4, filter_hi, 0);
    sum           = vmlaq_lane_s16(sum, s5, filter_hi, 1);
    sum           = vmlaq_lane_s16(sum, s6, filter_hi, 2);
    sum           = vmlaq_lane_s16(sum, s7, filter_hi, 3);

    // We halved the convolution filter values so -1 from the right shift.
    return vshrq_n_s16(sum, ROUND0_BITS - 1);
}

static INLINE void convolve_2d_sr_horiz_neon(const uint8_t *src, int src_stride, int16_t *im_block, int im_stride,
                                             int w, int im_h, const int16_t *x_filter_ptr) {
    const int bd = 8;

    const uint8_t *src_ptr    = src;
    int16_t       *dst_ptr    = im_block;
    int            dst_stride = im_stride;
    int            height     = im_h;

    if (w <= 4) {
        // A shim of 1 << ((ROUND0_BITS - 1) - 1) enables us to use non-rounding
        // shifts - which are generally faster than rounding shifts on modern CPUs.
        // (The extra -1 is needed because we halved the filter values.)
        const int16x4_t horiz_const = vdup_n_s16((1 << (bd + FILTER_BITS - 2)) + (1 << ((ROUND0_BITS - 1) - 1)));
        // 4-tap filters are used for blocks having width <= 4.
        // Filter values are even, so halve to reduce intermediate precision reqs.
        const int16x4_t x_filter = vshr_n_s16(vld1_s16(x_filter_ptr + 2), 1);

        src_ptr += 2;

        do {
            uint8x8_t t0 = vld1_u8(src_ptr); // a0 a1 a2 a3 a4 a5 a6 a7
            int16x4_t s0 = vget_low_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));
            int16x4_t s4 = vget_high_s16(vreinterpretq_s16_u16(vmovl_u8(t0)));

            int16x4_t s1 = vext_s16(s0, s4, 1); // a1 a2 a3 a4
            int16x4_t s2 = vext_s16(s0, s4, 2); // a2 a3 a4 a5
            int16x4_t s3 = vext_s16(s0, s4, 3); // a3 a4 a5 a6

            int16x4_t d0 = convolve4_4_2d_h(s0, s1, s2, s3, x_filter, horiz_const);

            vst1_s16(dst_ptr, d0);

            src_ptr += src_stride;
            dst_ptr += dst_stride;
        } while (--height != 0);
    } else {
        // A shim of 1 << ((ROUND0_BITS - 1) - 1) enables us to use non-rounding
        // shifts - which are generally faster than rounding shifts on modern CPUs.
        // (The extra -1 is needed because we halved the filter values.)
        const int16x8_t horiz_const = vdupq_n_s16((1 << (bd + FILTER_BITS - 2)) + (1 << ((ROUND0_BITS - 1) - 1)));
        // Filter values are even, so halve to reduce intermediate precision reqs.
        const int16x8_t x_filter = vshrq_n_s16(vld1q_s16(x_filter_ptr), 1);

        while (height > 8) {
            const uint8_t *s     = src_ptr;
            int16_t       *d     = dst_ptr;
            int            width = w;

            uint8x8_t t0, t1, t2, t3, t4, t5, t6, t7;
            load_u8_8x8(s, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
            transpose_elems_inplace_u8_8x8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);

            int16x8_t s0 = vreinterpretq_s16_u16(vmovl_u8(t0));
            int16x8_t s1 = vreinterpretq_s16_u16(vmovl_u8(t1));
            int16x8_t s2 = vreinterpretq_s16_u16(vmovl_u8(t2));
            int16x8_t s3 = vreinterpretq_s16_u16(vmovl_u8(t3));
            int16x8_t s4 = vreinterpretq_s16_u16(vmovl_u8(t4));
            int16x8_t s5 = vreinterpretq_s16_u16(vmovl_u8(t5));
            int16x8_t s6 = vreinterpretq_s16_u16(vmovl_u8(t6));

            s += 7;

            do {
                load_u8_8x8(s, src_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);

                transpose_elems_inplace_u8_8x8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);

                int16x8_t s7  = vreinterpretq_s16_u16(vmovl_u8(t0));
                int16x8_t s8  = vreinterpretq_s16_u16(vmovl_u8(t1));
                int16x8_t s9  = vreinterpretq_s16_u16(vmovl_u8(t2));
                int16x8_t s10 = vreinterpretq_s16_u16(vmovl_u8(t3));
                int16x8_t s11 = vreinterpretq_s16_u16(vmovl_u8(t4));
                int16x8_t s12 = vreinterpretq_s16_u16(vmovl_u8(t5));
                int16x8_t s13 = vreinterpretq_s16_u16(vmovl_u8(t6));
                int16x8_t s14 = vreinterpretq_s16_u16(vmovl_u8(t7));

                int16x8_t d0 = convolve8_8_2d_h(s0, s1, s2, s3, s4, s5, s6, s7, x_filter, horiz_const);
                int16x8_t d1 = convolve8_8_2d_h(s1, s2, s3, s4, s5, s6, s7, s8, x_filter, horiz_const);
                int16x8_t d2 = convolve8_8_2d_h(s2, s3, s4, s5, s6, s7, s8, s9, x_filter, horiz_const);
                int16x8_t d3 = convolve8_8_2d_h(s3, s4, s5, s6, s7, s8, s9, s10, x_filter, horiz_const);
                int16x8_t d4 = convolve8_8_2d_h(s4, s5, s6, s7, s8, s9, s10, s11, x_filter, horiz_const);
                int16x8_t d5 = convolve8_8_2d_h(s5, s6, s7, s8, s9, s10, s11, s12, x_filter, horiz_const);
                int16x8_t d6 = convolve8_8_2d_h(s6, s7, s8, s9, s10, s11, s12, s13, x_filter, horiz_const);
                int16x8_t d7 = convolve8_8_2d_h(s7, s8, s9, s10, s11, s12, s13, s14, x_filter, horiz_const);

                transpose_elems_inplace_s16_8x8(&d0, &d1, &d2, &d3, &d4, &d5, &d6, &d7);

                store_s16_8x8(d, dst_stride, d0, d1, d2, d3, d4, d5, d6, d7);

                s0 = s8;
                s1 = s9;
                s2 = s10;
                s3 = s11;
                s4 = s12;
                s5 = s13;
                s6 = s14;
                s += 8;
                d += 8;
                width -= 8;
            } while (width != 0);
            src_ptr += 8 * src_stride;
            dst_ptr += 8 * dst_stride;
            height -= 8;
        }

        do {
            const uint8_t *s     = src_ptr;
            int16_t       *d     = dst_ptr;
            int            width = w;

            uint8x8_t t0 = vld1_u8(s); // a0 a1 a2 a3 a4 a5 a6 a7
            int16x8_t s0 = vreinterpretq_s16_u16(vmovl_u8(t0));

            do {
                uint8x8_t t1 = vld1_u8(s + 8); // a8 a9 a10 a11 a12 a13 a14 a15
                int16x8_t s8 = vreinterpretq_s16_u16(vmovl_u8(t1));

                int16x8_t s1 = vextq_s16(s0, s8, 1); // a1 a2 a3 a4 a5 a6 a7 a8
                int16x8_t s2 = vextq_s16(s0, s8, 2); // a2 a3 a4 a5 a6 a7 a8 a9
                int16x8_t s3 = vextq_s16(s0, s8, 3); // a3 a4 a5 a6 a7 a8 a9 a10
                int16x8_t s4 = vextq_s16(s0, s8, 4); // a4 a5 a6 a7 a8 a9 a10 a11
                int16x8_t s5 = vextq_s16(s0, s8, 5); // a5 a6 a7 a8 a9 a10 a11 a12
                int16x8_t s6 = vextq_s16(s0, s8, 6); // a6 a7 a8 a9 a10 a11 a12 a13
                int16x8_t s7 = vextq_s16(s0, s8, 7); // a7 a8 a9 a10 a11 a12 a13 a14

                int16x8_t d0 = convolve8_8_2d_h(s0, s1, s2, s3, s4, s5, s6, s7, x_filter, horiz_const);

                vst1q_s16(d, d0);

                s0 = s8;
                s += 8;
                d += 8;
                width -= 8;
            } while (width != 0);
            src_ptr += src_stride;
            dst_ptr += dst_stride;
        } while (--height != 0);
    }
}

void svt_av1_convolve_2d_sr_neon(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w,
                                 int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y,
                                 const int32_t subpel_x_qn, const int32_t subpel_y_qn, ConvolveParams *conv_params) {
    if (w == 2 || h == 2) {
        svt_av1_convolve_2d_sr_c(src,
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

    const int      y_filter_taps  = get_filter_tap(filter_params_y, subpel_y_qn);
    const int      clamped_y_taps = y_filter_taps < 6 ? 6 : y_filter_taps;
    const int      im_h           = h + clamped_y_taps - 1;
    const int      im_stride      = MAX_SB_SIZE;
    const int      vert_offset    = clamped_y_taps / 2 - 1;
    const int      horiz_offset   = filter_params_x->taps / 2 - 1;
    const uint8_t *src_ptr        = src - vert_offset * src_stride - horiz_offset;

    const int16_t *x_filter_ptr = av1_get_interp_filter_subpel_kernel(*filter_params_x, subpel_x_qn & SUBPEL_MASK);
    const int16_t *y_filter_ptr = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_qn & SUBPEL_MASK);

    if (filter_params_x->taps > 8) {
        DECLARE_ALIGNED(16, int16_t, im_block[(MAX_SB_SIZE + MAX_FILTER_TAP - 1) * MAX_SB_SIZE]);

        const int16x8_t x_filter_0_7  = vld1q_s16(x_filter_ptr);
        const int16x4_t x_filter_8_11 = vld1_s16(x_filter_ptr + 8);
        const int16x8_t y_filter_0_7  = vld1q_s16(y_filter_ptr);
        const int16x4_t y_filter_8_11 = vld1_s16(y_filter_ptr + 8);

        convolve_2d_sr_horiz_12tap_neon(src_ptr, src_stride, im_block, im_stride, w, im_h, x_filter_0_7, x_filter_8_11);

        convolve_2d_sr_vert_12tap_neon(im_block, im_stride, dst, dst_stride, w, h, y_filter_0_7, y_filter_8_11);
    } else {
        DECLARE_ALIGNED(16, int16_t, im_block[(MAX_SB_SIZE + SUBPEL_TAPS - 1) * MAX_SB_SIZE]);

        convolve_2d_sr_horiz_neon(src_ptr, src_stride, im_block, im_stride, w, im_h, x_filter_ptr);

        const int16x8_t y_filter = vld1q_s16(y_filter_ptr);

        if (clamped_y_taps <= 6) {
            convolve_2d_sr_vert_6tap_neon(im_block, im_stride, dst, dst_stride, w, h, y_filter);
        } else {
            convolve_2d_sr_vert_8tap_neon(im_block, im_stride, dst, dst_stride, w, h, y_filter);
        }
    }
}

#define MAX_MASK_VALUE (1 << WEDGE_WEIGHT_BITS)

/**
 * See svt_av1_wedge_sse_from_residuals_c for details of the parameters and
 * computation.
 */
uint64_t svt_av1_wedge_sse_from_residuals_neon(const int16_t *r1, const int16_t *d, const uint8_t *m, int N) {
    assert(N % 64 == 0);

    uint64x2_t v_csse[2] = {vdupq_n_u64(0), vdupq_n_u64(0)};

    int i = 0;
    do {
        int32x4_t sum[4];
        int32x4_t sse[2];
        int16x4_t sum_s16[4];

        const int16x8_t r1_l = vld1q_s16(r1 + i);
        const int16x8_t r1_h = vld1q_s16(r1 + i + 8);
        const int16x8_t d_l  = vld1q_s16(d + i);
        const int16x8_t d_h  = vld1q_s16(d + i + 8);
        // The following three lines are a bit inelegant compared to using a pair
        // of vmovl_u8()... but it forces the compiler to emit a ZIP1, ZIP2 pair -
        // which can be executed in parallel with the subsequent SSHL instructions.
        // (SSHL can only be executed on half of the Neon pipes in modern Arm
        // cores, whereas ZIP1/2 can be executed on all of them.)
        const uint8x16x2_t m_u16 = vzipq_u8(vld1q_u8(m + i), vdupq_n_u8(0));
        const int16x8_t    m_l   = vreinterpretq_s16_u8(m_u16.val[0]);
        const int16x8_t    m_h   = vreinterpretq_s16_u8(m_u16.val[1]);

        sum[0] = vshll_n_s16(vget_low_s16(r1_l), WEDGE_WEIGHT_BITS);
        sum[1] = vshll_n_s16(vget_high_s16(r1_l), WEDGE_WEIGHT_BITS);
        sum[2] = vshll_n_s16(vget_low_s16(r1_h), WEDGE_WEIGHT_BITS);
        sum[3] = vshll_n_s16(vget_high_s16(r1_h), WEDGE_WEIGHT_BITS);

        sum[0] = vmlal_s16(sum[0], vget_low_s16(m_l), vget_low_s16(d_l));
        sum[1] = vmlal_s16(sum[1], vget_high_s16(m_l), vget_high_s16(d_l));
        sum[2] = vmlal_s16(sum[2], vget_low_s16(m_h), vget_low_s16(d_h));
        sum[3] = vmlal_s16(sum[3], vget_high_s16(m_h), vget_high_s16(d_h));

        sum_s16[0] = vqmovn_s32(sum[0]);
        sum_s16[1] = vqmovn_s32(sum[1]);
        sum_s16[2] = vqmovn_s32(sum[2]);
        sum_s16[3] = vqmovn_s32(sum[3]);

        sse[0] = vmull_s16(sum_s16[0], sum_s16[0]);
        sse[1] = vmull_s16(sum_s16[2], sum_s16[2]);
        sse[0] = vmlal_s16(sse[0], sum_s16[1], sum_s16[1]);
        sse[1] = vmlal_s16(sse[1], sum_s16[3], sum_s16[3]);

        v_csse[0] = vpadalq_u32(v_csse[0], vreinterpretq_u32_s32(sse[0]));
        v_csse[1] = vpadalq_u32(v_csse[1], vreinterpretq_u32_s32(sse[1]));

        i += 16;
    } while (i < N);

    uint64_t csse = horizontal_add_u64x2(vaddq_u64(v_csse[0], v_csse[1]));
    return ROUND_POWER_OF_TWO(csse, 2 * WEDGE_WEIGHT_BITS);
}

int svt_aom_satd_neon(const TranLow *coeff, int length) {
    const int32x4_t zero  = vdupq_n_s32(0);
    int32x4_t       accum = zero;
    do {
        const int32x4_t src0  = vld1q_s32(&coeff[0]);
        const int32x4_t src8  = vld1q_s32(&coeff[4]);
        const int32x4_t src16 = vld1q_s32(&coeff[8]);
        const int32x4_t src24 = vld1q_s32(&coeff[12]);
        accum                 = vabaq_s32(accum, src0, zero);
        accum                 = vabaq_s32(accum, src8, zero);
        accum                 = vabaq_s32(accum, src16, zero);
        accum                 = vabaq_s32(accum, src24, zero);
        length -= 16;
        coeff += 16;
    } while (length != 0);

    // satd: 26 bits, dynamic range [-32640 * 1024, 32640 * 1024]
    return horizontal_add_s32x4(accum);
}

void svt_aom_subtract_block_neon(int rows, int cols, int16_t *diff, ptrdiff_t diff_stride, const uint8_t *src,
                                 ptrdiff_t src_stride, const uint8_t *pred, ptrdiff_t pred_stride) {
    int r, c;
    if (cols > 16) {
        r = rows;
        do {
            c = 0;
            do {
                const uint8x16_t v_src_00     = vld1q_u8(&src[c + 0]);
                const uint8x16_t v_src_16     = vld1q_u8(&src[c + 16]);
                const uint8x16_t v_pred_00    = vld1q_u8(&pred[c + 0]);
                const uint8x16_t v_pred_16    = vld1q_u8(&pred[c + 16]);
                const uint16x8_t v_diff_lo_00 = vsubl_u8(vget_low_u8(v_src_00), vget_low_u8(v_pred_00));
                const uint16x8_t v_diff_hi_00 = vsubl_u8(vget_high_u8(v_src_00), vget_high_u8(v_pred_00));
                const uint16x8_t v_diff_lo_16 = vsubl_u8(vget_low_u8(v_src_16), vget_low_u8(v_pred_16));
                const uint16x8_t v_diff_hi_16 = vsubl_u8(vget_high_u8(v_src_16), vget_high_u8(v_pred_16));
                vst1q_s16(&diff[c + 0], vreinterpretq_s16_u16(v_diff_lo_00));
                vst1q_s16(&diff[c + 8], vreinterpretq_s16_u16(v_diff_hi_00));
                vst1q_s16(&diff[c + 16], vreinterpretq_s16_u16(v_diff_lo_16));
                vst1q_s16(&diff[c + 24], vreinterpretq_s16_u16(v_diff_hi_16));
                c += 32;
            } while (c < cols);
            diff += diff_stride;
            pred += pred_stride;
            src += src_stride;
        } while (--r != 0);
    } else if (cols > 8) {
        r = rows;
        do {
            const uint8x16_t v_src     = vld1q_u8(&src[0]);
            const uint8x16_t v_pred    = vld1q_u8(&pred[0]);
            const uint16x8_t v_diff_lo = vsubl_u8(vget_low_u8(v_src), vget_low_u8(v_pred));
            const uint16x8_t v_diff_hi = vsubl_u8(vget_high_u8(v_src), vget_high_u8(v_pred));
            vst1q_s16(&diff[0], vreinterpretq_s16_u16(v_diff_lo));
            vst1q_s16(&diff[8], vreinterpretq_s16_u16(v_diff_hi));
            diff += diff_stride;
            pred += pred_stride;
            src += src_stride;
        } while (--r != 0);
    } else if (cols > 4) {
        r = rows;
        do {
            const uint8x8_t  v_src  = vld1_u8(&src[0]);
            const uint8x8_t  v_pred = vld1_u8(&pred[0]);
            const uint16x8_t v_diff = vsubl_u8(v_src, v_pred);
            vst1q_s16(&diff[0], vreinterpretq_s16_u16(v_diff));
            diff += diff_stride;
            pred += pred_stride;
            src += src_stride;
        } while (--r != 0);
    } else {
        r = rows;
        do {
            c = 0;
            do { diff[c] = src[c] - pred[c]; } while (++c < cols);
            diff += diff_stride;
            pred += pred_stride;
            src += src_stride;
        } while (--r != 0);
    }
}

#if defined(__ARM_FEATURE_DOTPROD)

static INLINE void sse_16x1_neon(const uint8_t *src, const uint8_t *ref, uint32x4_t *sse) {
    uint8x16_t s = vld1q_u8(src);
    uint8x16_t r = vld1q_u8(ref);

    uint8x16_t abs_diff = vabdq_u8(s, r);

    *sse = vdotq_u32(*sse, abs_diff, abs_diff);
}

static INLINE void sse_8x1_neon(const uint8_t *src, const uint8_t *ref, uint32x2_t *sse) {
    uint8x8_t s = vld1_u8(src);
    uint8x8_t r = vld1_u8(ref);

    uint8x8_t abs_diff = vabd_u8(s, r);

    *sse = vdot_u32(*sse, abs_diff, abs_diff);
}

static INLINE void sse_4x2_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                uint32x2_t *sse) {
    uint8x8_t s = load_unaligned_u8(src, src_stride);
    uint8x8_t r = load_unaligned_u8(ref, ref_stride);

    uint8x8_t abs_diff = vabd_u8(s, r);

    *sse = vdot_u32(*sse, abs_diff, abs_diff);
}

static INLINE uint32_t sse_8xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                    int height) {
    uint32x2_t sse[2] = {vdup_n_u32(0), vdup_n_u32(0)};

    int i = height;
    do {
        sse_8x1_neon(src, ref, &sse[0]);
        src += src_stride;
        ref += ref_stride;
        sse_8x1_neon(src, ref, &sse[1]);
        src += src_stride;
        ref += ref_stride;
        i -= 2;
    } while (i != 0);

    return horizontal_add_u32x4(vcombine_u32(sse[0], sse[1]));
}

static INLINE uint32_t sse_4xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                    int height) {
    uint32x2_t sse = vdup_n_u32(0);

    int i = height;
    do {
        sse_4x2_neon(src, src_stride, ref, ref_stride, &sse);

        src += 2 * src_stride;
        ref += 2 * ref_stride;
        i -= 2;
    } while (i != 0);

    return horizontal_add_u32x2(sse);
}

static INLINE uint32_t sse_wxh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int width,
                                    int height) {
    uint32x2_t sse[2] = {vdup_n_u32(0), vdup_n_u32(0)};
    int        i, j;
    if ((width & 0x07) && ((width & 0x07) < 5)) {
        i = height;
        do {
            j = 0;
            do {
                sse_8x1_neon(src + j, ref + j, &sse[0]);
                sse_8x1_neon(src + j + src_stride, ref + j + ref_stride, &sse[1]);
                j += 8;
            } while (j + 4 < width);

            sse_4x2_neon(src + j, src_stride, ref + j, ref_stride, &sse[0]);
            src += 2 * src_stride;
            ref += 2 * ref_stride;
            i -= 2;
        } while (i != 0);
    } else {
        i = height;
        do {
            j = 0;
            do {
                sse_8x1_neon(src + j, ref + j, &sse[0]);
                sse_8x1_neon(src + j + src_stride, ref + j + ref_stride, &sse[1]);
                j += 8;
            } while (j < width);

            src += 2 * src_stride;
            ref += 2 * ref_stride;
            i -= 2;
        } while (i != 0);
    }
    return horizontal_add_u32x4(vcombine_u32(sse[0], sse[1]));
}

#else // !defined(__ARM_FEATURE_DOTPROD)

static INLINE void sse_16x1_neon(const uint8_t *src, const uint8_t *ref, uint32x4_t *sse) {
    uint8x16_t s = vld1q_u8(src);
    uint8x16_t r = vld1q_u8(ref);

    uint8x16_t abs_diff    = vabdq_u8(s, r);
    uint8x8_t  abs_diff_lo = vget_low_u8(abs_diff);
    uint8x8_t  abs_diff_hi = vget_high_u8(abs_diff);

    *sse = vpadalq_u16(*sse, vmull_u8(abs_diff_lo, abs_diff_lo));
    *sse = vpadalq_u16(*sse, vmull_u8(abs_diff_hi, abs_diff_hi));
}

static INLINE void sse_8x1_neon(const uint8_t *src, const uint8_t *ref, uint32x4_t *sse) {
    uint8x8_t s = vld1_u8(src);
    uint8x8_t r = vld1_u8(ref);

    uint8x8_t abs_diff = vabd_u8(s, r);

    *sse = vpadalq_u16(*sse, vmull_u8(abs_diff, abs_diff));
}

static INLINE void sse_4x2_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                uint32x4_t *sse) {
    uint8x8_t s = load_unaligned_u8(src, src_stride);
    uint8x8_t r = load_unaligned_u8(ref, ref_stride);

    uint8x8_t abs_diff = vabd_u8(s, r);

    *sse = vpadalq_u16(*sse, vmull_u8(abs_diff, abs_diff));
}

static INLINE uint32_t sse_8xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                    int height) {
    uint32x4_t sse = vdupq_n_u32(0);

    int i = height;
    do {
        sse_8x1_neon(src, ref, &sse);

        src += src_stride;
        ref += ref_stride;
    } while (--i != 0);

    return horizontal_add_u32x4(sse);
}

static INLINE uint32_t sse_4xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                    int height) {
    uint32x4_t sse = vdupq_n_u32(0);

    int i = height;
    do {
        sse_4x2_neon(src, src_stride, ref, ref_stride, &sse);

        src += 2 * src_stride;
        ref += 2 * ref_stride;
        i -= 2;
    } while (i != 0);

    return horizontal_add_u32x4(sse);
}

static INLINE uint32_t sse_wxh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int width,
                                    int height) {
    uint32x4_t sse = vdupq_n_u32(0);
    int        i, j;
    if ((width & 0x07) && ((width & 0x07) < 5)) {
        i = height;
        do {
            j = 0;
            do {
                sse_8x1_neon(src + j, ref + j, &sse);
                sse_8x1_neon(src + j + src_stride, ref + j + ref_stride, &sse);
                j += 8;
            } while (j + 4 < width);

            sse_4x2_neon(src + j, src_stride, ref + j, ref_stride, &sse);
            src += 2 * src_stride;
            ref += 2 * ref_stride;
            i -= 2;
        } while (i != 0);
    } else {
        i = height;
        do {
            j = 0;
            do {
                sse_8x1_neon(src + j, ref + j, &sse);
                j += 8;
            } while (j < width);

            src += src_stride;
            ref += ref_stride;
        } while (--i != 0);
    }
    return horizontal_add_u32x4(sse);
}

#endif // defined(__ARM_FEATURE_DOTPROD)

static INLINE uint32_t sse_128xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                      int height) {
    uint32x4_t sse[2] = {vdupq_n_u32(0), vdupq_n_u32(0)};

    int i = height;
    do {
        sse_16x1_neon(src, ref, &sse[0]);
        sse_16x1_neon(src + 16, ref + 16, &sse[1]);
        sse_16x1_neon(src + 32, ref + 32, &sse[0]);
        sse_16x1_neon(src + 48, ref + 48, &sse[1]);
        sse_16x1_neon(src + 64, ref + 64, &sse[0]);
        sse_16x1_neon(src + 80, ref + 80, &sse[1]);
        sse_16x1_neon(src + 96, ref + 96, &sse[0]);
        sse_16x1_neon(src + 112, ref + 112, &sse[1]);

        src += src_stride;
        ref += ref_stride;
    } while (--i != 0);

    return horizontal_add_u32x4(vaddq_u32(sse[0], sse[1]));
}

static INLINE uint32_t sse_64xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                     int height) {
    uint32x4_t sse[2] = {vdupq_n_u32(0), vdupq_n_u32(0)};

    int i = height;
    do {
        sse_16x1_neon(src, ref, &sse[0]);
        sse_16x1_neon(src + 16, ref + 16, &sse[1]);
        sse_16x1_neon(src + 32, ref + 32, &sse[0]);
        sse_16x1_neon(src + 48, ref + 48, &sse[1]);

        src += src_stride;
        ref += ref_stride;
    } while (--i != 0);

    return horizontal_add_u32x4(vaddq_u32(sse[0], sse[1]));
}

static INLINE uint32_t sse_32xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                     int height) {
    uint32x4_t sse[2] = {vdupq_n_u32(0), vdupq_n_u32(0)};

    int i = height;
    do {
        sse_16x1_neon(src, ref, &sse[0]);
        sse_16x1_neon(src + 16, ref + 16, &sse[1]);

        src += src_stride;
        ref += ref_stride;
    } while (--i != 0);

    return horizontal_add_u32x4(vaddq_u32(sse[0], sse[1]));
}

static INLINE uint32_t sse_16xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                     int height) {
    uint32x4_t sse[2] = {vdupq_n_u32(0), vdupq_n_u32(0)};

    int i = height;
    do {
        sse_16x1_neon(src, ref, &sse[0]);
        src += src_stride;
        ref += ref_stride;
        sse_16x1_neon(src, ref, &sse[1]);
        src += src_stride;
        ref += ref_stride;
        i -= 2;
    } while (i != 0);

    return horizontal_add_u32x4(vaddq_u32(sse[0], sse[1]));
}

int64_t svt_aom_sse_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int width,
                         int height) {
    switch (width) {
    case 4: {
        return sse_4xh_neon(src, src_stride, ref, ref_stride, height);
    }
    case 8: {
        return sse_8xh_neon(src, src_stride, ref, ref_stride, height);
    }
    case 16: {
        return sse_16xh_neon(src, src_stride, ref, ref_stride, height);
    }
    case 32: {
        return sse_32xh_neon(src, src_stride, ref, ref_stride, height);
    }
    case 64: {
        return sse_64xh_neon(src, src_stride, ref, ref_stride, height);
    }
    case 128: {
        return sse_128xh_neon(src, src_stride, ref, ref_stride, height);
    }
    default: {
        return sse_wxh_neon(src, src_stride, ref, ref_stride, width, height);
    }
    }
}

int64_t svt_av1_block_error_neon(const TranLow *coeff, const TranLow *dqcoeff, intptr_t block_size, int64_t *ssz) {
    int64x2_t error   = vdupq_n_s64(0);
    int64x2_t sqcoeff = vdupq_n_s64(0);

    assert(block_size >= 8);
    assert((block_size % 8) == 0);

    do {
        const int16x8_t c       = load_tran_low_to_s16q(coeff);
        const int16x8_t d       = load_tran_low_to_s16q(dqcoeff);
        const int16x8_t diff    = vsubq_s16(c, d);
        const int16x4_t diff_lo = vget_low_s16(diff);
        const int16x4_t diff_hi = vget_high_s16(diff);
        // diff is 15-bits, the squares 30, so we can store 2 in 31-bits before
        // accumulating them in 64-bits.
        const int32x4_t err0 = vmull_s16(diff_lo, diff_lo);
        const int32x4_t err1 = vmlal_s16(err0, diff_hi, diff_hi);
        const int64x2_t err2 = vaddl_s32(vget_low_s32(err1), vget_high_s32(err1));
        error                = vaddq_s64(error, err2);

        const int16x4_t coeff_lo = vget_low_s16(c);
        const int16x4_t coeff_hi = vget_high_s16(c);
        const int32x4_t sqcoeff0 = vmull_s16(coeff_lo, coeff_lo);
        const int32x4_t sqcoeff1 = vmlal_s16(sqcoeff0, coeff_hi, coeff_hi);
        const int64x2_t sqcoeff2 = vaddl_s32(vget_low_s32(sqcoeff1), vget_high_s32(sqcoeff1));
        sqcoeff                  = vaddq_s64(sqcoeff, sqcoeff2);

        coeff += 8;
        dqcoeff += 8;
        block_size -= 8;
    } while (block_size != 0);

    *ssz = vaddvq_s64(sqcoeff);
    return vaddvq_s64(error);
}

int8_t svt_av1_wedge_sign_from_residuals_neon(const int16_t *ds, const uint8_t *m, int N, int64_t limit) {
    int64x2_t sum;
    int32x4_t acc[4] = {vdupq_n_s32(0), vdupq_n_s32(0), vdupq_n_s32(0), vdupq_n_s32(0)};

    do {
        int16x8_t ds_l = vld1q_s16(ds);
        int16x8_t ds_h = vld1q_s16(ds + 8);

        int8x16_t m_s8 = vreinterpretq_s8_u8(vld1q_u8(m));
        int16x8_t m_l  = vmovl_s8(vget_low_s8(m_s8));
        int16x8_t m_h  = vmovl_s8(vget_high_s8(m_s8));

        acc[0] = vmlal_s16(acc[0], vget_low_s16(ds_l), vget_low_s16(m_l));
        acc[1] = vmlal_s16(acc[1], vget_high_s16(ds_l), vget_high_s16(m_l));
        acc[2] = vmlal_s16(acc[2], vget_low_s16(ds_h), vget_low_s16(m_h));
        acc[3] = vmlal_s16(acc[3], vget_high_s16(ds_h), vget_high_s16(m_h));

        ds += 16;
        m += 16;
        N -= 16;
    } while (N != 0);

    sum = vpaddlq_s32(acc[0]);
    sum = vpadalq_s32(sum, acc[1]);
    sum = vpadalq_s32(sum, acc[2]);
    sum = vpadalq_s32(sum, acc[3]);

    return (horizontal_add_s64x2(sum) > limit);
}

uint8_t svt_av1_compute_cul_level_neon(const int16_t *const scan, const int32_t *const quant_coeff, uint16_t *eob) {
    if (*eob == 1) {
        if (quant_coeff[0] > 0)
            return (AOMMIN(COEFF_CONTEXT_MASK, quant_coeff[0]) + (2 << COEFF_CONTEXT_BITS));
        if (quant_coeff[0] < 0) {
            return (AOMMIN(COEFF_CONTEXT_MASK, ABS(quant_coeff[0])) | (1 << COEFF_CONTEXT_BITS));
        }
        return 0;
    }

    int32x4_t sum_256_0 = vdupq_n_s32(0);
    int32x4_t sum_256_1 = vdupq_n_s32(0);

    for (int32_t c = 0; c < *eob; c += 8) {
        const int16x4_t scan_64_0 = vld1_s16(scan + c + 0);
        const int16x4_t scan_64_1 = vld1_s16(scan + c + 4);

        const int32x4_t scan_128_0 = vmovl_s16(scan_64_0);
        const int32x4_t scan_128_1 = vmovl_s16(scan_64_1);

        // no gather operation, unfortunately
        const int32_t quant_coeff_0 = *(quant_coeff + vgetq_lane_s32(scan_128_0, 0) + 4 * 8);
        const int32_t quant_coeff_1 = *(quant_coeff + vgetq_lane_s32(scan_128_0, 1) + 4 * 8);
        const int32_t quant_coeff_2 = *(quant_coeff + vgetq_lane_s32(scan_128_0, 2) + 4 * 8);
        const int32_t quant_coeff_3 = *(quant_coeff + vgetq_lane_s32(scan_128_0, 3) + 4 * 8);

        const int32_t quant_coeff_4 = *(quant_coeff + vgetq_lane_s32(scan_128_1, 0) + 4 * 8);
        const int32_t quant_coeff_5 = *(quant_coeff + vgetq_lane_s32(scan_128_1, 1) + 4 * 8);
        const int32_t quant_coeff_6 = *(quant_coeff + vgetq_lane_s32(scan_128_1, 2) + 4 * 8);
        const int32_t quant_coeff_7 = *(quant_coeff + vgetq_lane_s32(scan_128_1, 3) + 4 * 8);

        int32x4_t quant_coeff_128_0 = vcombine_s32(
            vcreate_s32((((uint64_t)quant_coeff_1) << 32) + (uint64_t)quant_coeff_0),
            vcreate_s32((((uint64_t)quant_coeff_3) << 32) + (uint64_t)quant_coeff_2));
        int32x4_t quant_coeff_128_1 = vcombine_s32(
            vcreate_s32((((uint64_t)quant_coeff_5) << 32) + (uint64_t)quant_coeff_4),
            vcreate_s32((((uint64_t)quant_coeff_7) << 32) + (uint64_t)quant_coeff_6));

        quant_coeff_128_0 = vabsq_s32(quant_coeff_128_0);
        quant_coeff_128_1 = vabsq_s32(quant_coeff_128_1);

        sum_256_0 = vaddq_s32(sum_256_0, quant_coeff_128_0);
        sum_256_1 = vaddq_s32(sum_256_1, quant_coeff_128_1);
    }

    int32x4_t partial_sums = vaddq_s32(sum_256_0, sum_256_1);
    partial_sums           = vpaddq_s32(partial_sums, partial_sums);
    partial_sums           = vpaddq_s32(partial_sums, partial_sums);

    const int32_t cul_level = AOMMIN(COEFF_CONTEXT_MASK, vgetq_lane_s32(partial_sums, 0));

    // DC value, calculation from set_dc_sign()
    if (quant_coeff[0] < 0)
        return (cul_level | (1 << COEFF_CONTEXT_BITS));
    if (quant_coeff[0] > 0)
        return (cul_level + (2 << COEFF_CONTEXT_BITS));
    return (uint8_t)cul_level;
}
