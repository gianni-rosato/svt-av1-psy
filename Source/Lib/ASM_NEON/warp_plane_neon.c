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

#include <assert.h>
#include <arm_neon.h>
#include <memory.h>
#include <math.h>

#include "definitions.h"
#include "transpose_neon.h"
#include "warped_motion.h"
#include "convolve.h"
#include "sum_neon.h"

static INLINE void load_filters_4(int16x8_t out[], int offset, int stride) {
    out[0] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 0 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[1] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 1 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[2] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 2 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[3] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 3 * stride) >> WARPEDDIFF_PREC_BITS)));
}

static INLINE void load_filters_8(int16x8_t out[], int offset, int stride) {
    out[0] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 0 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[1] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 1 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[2] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 2 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[3] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 3 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[4] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 4 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[5] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 5 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[6] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 6 * stride) >> WARPEDDIFF_PREC_BITS)));
    out[7] = vld1q_s16((int16_t *)(svt_aom_warped_filter + ((offset + 7 * stride) >> WARPEDDIFF_PREC_BITS)));
}

static INLINE int16x8_t horizontal_filter_4x1_f4(const uint8x16_t in, int sx, int alpha) {
    const int32x4_t add_const = vdupq_n_s32(1 << (8 + FILTER_BITS - 1));

    // Loading the 8 filter taps
    int16x8_t f[4];
    load_filters_4(f, sx, alpha);

    int16x8_t in16_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(in)));
    int16x8_t in16_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(in)));

    int16x8_t m0 = vmulq_s16(f[0], in16_lo);
    int16x8_t m1 = vmulq_s16(f[1], vextq_s16(in16_lo, in16_hi, 1));
    int16x8_t m2 = vmulq_s16(f[2], vextq_s16(in16_lo, in16_hi, 2));
    int16x8_t m3 = vmulq_s16(f[3], vextq_s16(in16_lo, in16_hi, 3));

    int32x4_t m0123_pairs[] = {vpaddlq_s16(m0), vpaddlq_s16(m1), vpaddlq_s16(m2), vpaddlq_s16(m3)};

    int32x4_t tmp_res_low = horizontal_add_4d_s32x4(m0123_pairs);

    tmp_res_low = vaddq_s32(tmp_res_low, add_const);

    uint16x8_t res = vcombine_u16(vqrshrun_n_s32(tmp_res_low, ROUND0_BITS), vdup_n_u16(0));
    return vreinterpretq_s16_u16(res);
}

static INLINE int16x8_t horizontal_filter_8x1_f8(const uint8x16_t in, int sx, int alpha) {
    const int32x4_t add_const = vdupq_n_s32(1 << (8 + FILTER_BITS - 1));

    // Loading the 8 filter taps
    int16x8_t f[8];
    load_filters_8(f, sx, alpha);

    int16x8_t in16_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(in)));
    int16x8_t in16_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(in)));

    int16x8_t m0 = vmulq_s16(f[0], in16_lo);
    int16x8_t m1 = vmulq_s16(f[1], vextq_s16(in16_lo, in16_hi, 1));
    int16x8_t m2 = vmulq_s16(f[2], vextq_s16(in16_lo, in16_hi, 2));
    int16x8_t m3 = vmulq_s16(f[3], vextq_s16(in16_lo, in16_hi, 3));
    int16x8_t m4 = vmulq_s16(f[4], vextq_s16(in16_lo, in16_hi, 4));
    int16x8_t m5 = vmulq_s16(f[5], vextq_s16(in16_lo, in16_hi, 5));
    int16x8_t m6 = vmulq_s16(f[6], vextq_s16(in16_lo, in16_hi, 6));
    int16x8_t m7 = vmulq_s16(f[7], vextq_s16(in16_lo, in16_hi, 7));

    int32x4_t m0123_pairs[] = {vpaddlq_s16(m0), vpaddlq_s16(m1), vpaddlq_s16(m2), vpaddlq_s16(m3)};
    int32x4_t m4567_pairs[] = {vpaddlq_s16(m4), vpaddlq_s16(m5), vpaddlq_s16(m6), vpaddlq_s16(m7)};

    int32x4_t tmp_res_low  = horizontal_add_4d_s32x4(m0123_pairs);
    int32x4_t tmp_res_high = horizontal_add_4d_s32x4(m4567_pairs);

    tmp_res_low  = vaddq_s32(tmp_res_low, add_const);
    tmp_res_high = vaddq_s32(tmp_res_high, add_const);

    uint16x8_t res = vcombine_u16(vqrshrun_n_s32(tmp_res_low, ROUND0_BITS), vqrshrun_n_s32(tmp_res_high, ROUND0_BITS));
    return vreinterpretq_s16_u16(res);
}

static INLINE int16x8_t horizontal_filter_4x1_f1(const uint8x16_t in, int sx) {
    const int32x4_t add_const = vdupq_n_s32(1 << (8 + FILTER_BITS - 1));

    int16x8_t f_s16 = vld1q_s16((int16_t *)(svt_aom_warped_filter + (sx >> WARPEDDIFF_PREC_BITS)));

    int16x8_t in16_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(in)));
    int16x8_t in16_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(in)));

    int16x8_t m0 = vmulq_s16(f_s16, in16_lo);
    int16x8_t m1 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 1));
    int16x8_t m2 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 2));
    int16x8_t m3 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 3));

    int32x4_t m0123_pairs[] = {vpaddlq_s16(m0), vpaddlq_s16(m1), vpaddlq_s16(m2), vpaddlq_s16(m3)};

    int32x4_t tmp_res_low = horizontal_add_4d_s32x4(m0123_pairs);

    tmp_res_low = vaddq_s32(tmp_res_low, add_const);

    uint16x8_t res = vcombine_u16(vqrshrun_n_s32(tmp_res_low, ROUND0_BITS), vdup_n_u16(0));
    return vreinterpretq_s16_u16(res);
}

static INLINE int16x8_t horizontal_filter_8x1_f1(const uint8x16_t in, int sx) {
    const int32x4_t add_const = vdupq_n_s32(1 << (8 + FILTER_BITS - 1));

    int16x8_t f_s16 = vld1q_s16((int16_t *)(svt_aom_warped_filter + (sx >> WARPEDDIFF_PREC_BITS)));

    int16x8_t in16_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(in)));
    int16x8_t in16_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(in)));

    int16x8_t m0 = vmulq_s16(f_s16, in16_lo);
    int16x8_t m1 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 1));
    int16x8_t m2 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 2));
    int16x8_t m3 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 3));
    int16x8_t m4 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 4));
    int16x8_t m5 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 5));
    int16x8_t m6 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 6));
    int16x8_t m7 = vmulq_s16(f_s16, vextq_s16(in16_lo, in16_hi, 7));

    int32x4_t m0123_pairs[] = {vpaddlq_s16(m0), vpaddlq_s16(m1), vpaddlq_s16(m2), vpaddlq_s16(m3)};
    int32x4_t m4567_pairs[] = {vpaddlq_s16(m4), vpaddlq_s16(m5), vpaddlq_s16(m6), vpaddlq_s16(m7)};

    int32x4_t tmp_res_low  = horizontal_add_4d_s32x4(m0123_pairs);
    int32x4_t tmp_res_high = horizontal_add_4d_s32x4(m4567_pairs);

    tmp_res_low  = vaddq_s32(tmp_res_low, add_const);
    tmp_res_high = vaddq_s32(tmp_res_high, add_const);

    uint16x8_t res = vcombine_u16(vqrshrun_n_s32(tmp_res_low, ROUND0_BITS), vqrshrun_n_s32(tmp_res_high, ROUND0_BITS));
    return vreinterpretq_s16_u16(res);
}

static INLINE void vertical_filter_4x1_f1(const int16x8_t *src, int32x4_t *res, int sy) {
    int16x4_t s0 = vget_low_s16(src[0]);
    int16x4_t s1 = vget_low_s16(src[1]);
    int16x4_t s2 = vget_low_s16(src[2]);
    int16x4_t s3 = vget_low_s16(src[3]);
    int16x4_t s4 = vget_low_s16(src[4]);
    int16x4_t s5 = vget_low_s16(src[5]);
    int16x4_t s6 = vget_low_s16(src[6]);
    int16x4_t s7 = vget_low_s16(src[7]);

    int16x8_t f = vld1q_s16((int16_t *)(svt_aom_warped_filter + (sy >> WARPEDDIFF_PREC_BITS)));

    int32x4_t m0123 = vmull_lane_s16(s0, vget_low_s16(f), 0);
    m0123           = vmlal_lane_s16(m0123, s1, vget_low_s16(f), 1);
    m0123           = vmlal_lane_s16(m0123, s2, vget_low_s16(f), 2);
    m0123           = vmlal_lane_s16(m0123, s3, vget_low_s16(f), 3);
    m0123           = vmlal_lane_s16(m0123, s4, vget_high_s16(f), 0);
    m0123           = vmlal_lane_s16(m0123, s5, vget_high_s16(f), 1);
    m0123           = vmlal_lane_s16(m0123, s6, vget_high_s16(f), 2);
    m0123           = vmlal_lane_s16(m0123, s7, vget_high_s16(f), 3);

    *res = m0123;
}

static INLINE void vertical_filter_4x1_f4(const int16x8_t *src, int32x4_t *res, int sy, int gamma) {
    int16x8_t s0, s1, s2, s3;
    transpose_elems_s16_4x8(vget_low_s16(src[0]),
                            vget_low_s16(src[1]),
                            vget_low_s16(src[2]),
                            vget_low_s16(src[3]),
                            vget_low_s16(src[4]),
                            vget_low_s16(src[5]),
                            vget_low_s16(src[6]),
                            vget_low_s16(src[7]),
                            &s0,
                            &s1,
                            &s2,
                            &s3);

    int16x8_t f[4];
    load_filters_4(f, sy, gamma);

    int32x4_t m0 = vmull_s16(vget_low_s16(s0), vget_low_s16(f[0]));
    m0           = vmlal_s16(m0, vget_high_s16(s0), vget_high_s16(f[0]));
    int32x4_t m1 = vmull_s16(vget_low_s16(s1), vget_low_s16(f[1]));
    m1           = vmlal_s16(m1, vget_high_s16(s1), vget_high_s16(f[1]));
    int32x4_t m2 = vmull_s16(vget_low_s16(s2), vget_low_s16(f[2]));
    m2           = vmlal_s16(m2, vget_high_s16(s2), vget_high_s16(f[2]));
    int32x4_t m3 = vmull_s16(vget_low_s16(s3), vget_low_s16(f[3]));
    m3           = vmlal_s16(m3, vget_high_s16(s3), vget_high_s16(f[3]));

    int32x4_t m0123_pairs[] = {m0, m1, m2, m3};

    *res = horizontal_add_4d_s32x4(m0123_pairs);
}

static INLINE void vertical_filter_8x1_f1(const int16x8_t *src, int32x4_t *res_low, int32x4_t *res_high, int sy) {
    int16x8_t s0 = src[0];
    int16x8_t s1 = src[1];
    int16x8_t s2 = src[2];
    int16x8_t s3 = src[3];
    int16x8_t s4 = src[4];
    int16x8_t s5 = src[5];
    int16x8_t s6 = src[6];
    int16x8_t s7 = src[7];

    int16x8_t f = vld1q_s16((int16_t *)(svt_aom_warped_filter + (sy >> WARPEDDIFF_PREC_BITS)));

    int32x4_t m0123 = vmull_lane_s16(vget_low_s16(s0), vget_low_s16(f), 0);
    m0123           = vmlal_lane_s16(m0123, vget_low_s16(s1), vget_low_s16(f), 1);
    m0123           = vmlal_lane_s16(m0123, vget_low_s16(s2), vget_low_s16(f), 2);
    m0123           = vmlal_lane_s16(m0123, vget_low_s16(s3), vget_low_s16(f), 3);
    m0123           = vmlal_lane_s16(m0123, vget_low_s16(s4), vget_high_s16(f), 0);
    m0123           = vmlal_lane_s16(m0123, vget_low_s16(s5), vget_high_s16(f), 1);
    m0123           = vmlal_lane_s16(m0123, vget_low_s16(s6), vget_high_s16(f), 2);
    m0123           = vmlal_lane_s16(m0123, vget_low_s16(s7), vget_high_s16(f), 3);

    int32x4_t m4567 = vmull_lane_s16(vget_high_s16(s0), vget_low_s16(f), 0);
    m4567           = vmlal_lane_s16(m4567, vget_high_s16(s1), vget_low_s16(f), 1);
    m4567           = vmlal_lane_s16(m4567, vget_high_s16(s2), vget_low_s16(f), 2);
    m4567           = vmlal_lane_s16(m4567, vget_high_s16(s3), vget_low_s16(f), 3);
    m4567           = vmlal_lane_s16(m4567, vget_high_s16(s4), vget_high_s16(f), 0);
    m4567           = vmlal_lane_s16(m4567, vget_high_s16(s5), vget_high_s16(f), 1);
    m4567           = vmlal_lane_s16(m4567, vget_high_s16(s6), vget_high_s16(f), 2);
    m4567           = vmlal_lane_s16(m4567, vget_high_s16(s7), vget_high_s16(f), 3);

    *res_low  = m0123;
    *res_high = m4567;
}

static INLINE void vertical_filter_8x1_f8(const int16x8_t *src, int32x4_t *res_low, int32x4_t *res_high, int sy,
                                          int gamma) {
    int16x8_t s0 = src[0];
    int16x8_t s1 = src[1];
    int16x8_t s2 = src[2];
    int16x8_t s3 = src[3];
    int16x8_t s4 = src[4];
    int16x8_t s5 = src[5];
    int16x8_t s6 = src[6];
    int16x8_t s7 = src[7];
    transpose_elems_inplace_s16_8x8(&s0, &s1, &s2, &s3, &s4, &s5, &s6, &s7);

    int16x8_t f[8];
    load_filters_8(f, sy, gamma);

    int32x4_t m0 = vmull_s16(vget_low_s16(s0), vget_low_s16(f[0]));
    m0           = vmlal_s16(m0, vget_high_s16(s0), vget_high_s16(f[0]));
    int32x4_t m1 = vmull_s16(vget_low_s16(s1), vget_low_s16(f[1]));
    m1           = vmlal_s16(m1, vget_high_s16(s1), vget_high_s16(f[1]));
    int32x4_t m2 = vmull_s16(vget_low_s16(s2), vget_low_s16(f[2]));
    m2           = vmlal_s16(m2, vget_high_s16(s2), vget_high_s16(f[2]));
    int32x4_t m3 = vmull_s16(vget_low_s16(s3), vget_low_s16(f[3]));
    m3           = vmlal_s16(m3, vget_high_s16(s3), vget_high_s16(f[3]));
    int32x4_t m4 = vmull_s16(vget_low_s16(s4), vget_low_s16(f[4]));
    m4           = vmlal_s16(m4, vget_high_s16(s4), vget_high_s16(f[4]));
    int32x4_t m5 = vmull_s16(vget_low_s16(s5), vget_low_s16(f[5]));
    m5           = vmlal_s16(m5, vget_high_s16(s5), vget_high_s16(f[5]));
    int32x4_t m6 = vmull_s16(vget_low_s16(s6), vget_low_s16(f[6]));
    m6           = vmlal_s16(m6, vget_high_s16(s6), vget_high_s16(f[6]));
    int32x4_t m7 = vmull_s16(vget_low_s16(s7), vget_low_s16(f[7]));
    m7           = vmlal_s16(m7, vget_high_s16(s7), vget_high_s16(f[7]));

    int32x4_t m0123_pairs[] = {m0, m1, m2, m3};
    int32x4_t m4567_pairs[] = {m4, m5, m6, m7};

    *res_low  = horizontal_add_4d_s32x4(m0123_pairs);
    *res_high = horizontal_add_4d_s32x4(m4567_pairs);
}

static INLINE int clamp_iy(int iy, int height) { return clamp(iy, 0, height - 1); }

static INLINE void warp_affine_horizontal(const uint8_t *ref, int width, int height, int stride, int p_width,
                                          int p_height, int16_t alpha, int16_t beta, const int64_t x4, const int64_t y4,
                                          const int i, int16x8_t tmp[], const uint8x16_t indx_vec) {
    const int bd                = 8;
    const int reduce_bits_horiz = ROUND0_BITS;
    const int height_limit      = AOMMIN(8, p_height - i) + 7;

    int32_t ix4 = (int32_t)(x4 >> WARPEDMODEL_PREC_BITS);
    int32_t iy4 = (int32_t)(y4 >> WARPEDMODEL_PREC_BITS);

    int32_t sx4 = x4 & ((1 << WARPEDMODEL_PREC_BITS) - 1);
    sx4 += alpha * (-4) + beta * (-4) + (1 << (WARPEDDIFF_PREC_BITS - 1)) +
        (WARPEDPIXEL_PREC_SHIFTS << WARPEDDIFF_PREC_BITS);
    sx4 &= ~((1 << WARP_PARAM_REDUCE_BITS) - 1);

    if (ix4 <= -7) {
        for (int k = 0; k < height_limit; ++k) {
            int     iy      = clamp_iy(iy4 + k - 7, height);
            int16_t dup_val = (1 << (bd + FILTER_BITS - reduce_bits_horiz - 1)) +
                ref[iy * stride] * (1 << (FILTER_BITS - reduce_bits_horiz));
            tmp[k] = vdupq_n_s16(dup_val);
        }
        return;
    } else if (ix4 >= width + 6) {
        for (int k = 0; k < height_limit; ++k) {
            int     iy      = clamp_iy(iy4 + k - 7, height);
            int16_t dup_val = (1 << (bd + FILTER_BITS - reduce_bits_horiz - 1)) +
                ref[iy * stride + (width - 1)] * (1 << (FILTER_BITS - reduce_bits_horiz));
            tmp[k] = vdupq_n_s16(dup_val);
        }
        return;
    }

    uint8x16_t in[15];
    if (((ix4 - 7) < 0) || ((ix4 + 9) > width)) {
        const int out_of_boundary_left  = -(ix4 - 6);
        const int out_of_boundary_right = (ix4 + 8) - width;

        for (int k = 0; k < height_limit; ++k) {
            const int      iy    = clamp_iy(iy4 + k - 7, height);
            const uint8_t *src   = ref + iy * stride + ix4 - 7;
            uint8x16_t     src_1 = vld1q_u8(src);

            if (out_of_boundary_left >= 0) {
                int        limit    = out_of_boundary_left + 1;
                uint8x16_t cmp_vec  = vdupq_n_u8(out_of_boundary_left);
                uint8x16_t vec_dup  = vdupq_n_u8(*(src + limit));
                uint8x16_t mask_val = vcleq_u8(indx_vec, cmp_vec);
                src_1               = vbslq_u8(mask_val, vec_dup, src_1);
            }
            if (out_of_boundary_right >= 0) {
                int        limit    = 15 - (out_of_boundary_right + 1);
                uint8x16_t cmp_vec  = vdupq_n_u8(15 - out_of_boundary_right);
                uint8x16_t vec_dup  = vdupq_n_u8(*(src + limit));
                uint8x16_t mask_val = vcgeq_u8(indx_vec, cmp_vec);
                src_1               = vbslq_u8(mask_val, vec_dup, src_1);
            }
            in[k] = src_1;
        }
    } else {
        for (int k = 0; k < height_limit; ++k) {
            const int      iy  = clamp_iy(iy4 + k - 7, height);
            const uint8_t *src = ref + iy * stride + ix4 - 7;
            in[k]              = vld1q_u8(src);
        }
    }

    if (p_width == 4) {
        if (beta == 0) {
            if (alpha == 0) {
                for (int k = 0; k < height_limit; ++k) { tmp[k] = horizontal_filter_4x1_f1(in[k], sx4); }
            } else {
                for (int k = 0; k < height_limit; ++k) { tmp[k] = horizontal_filter_4x1_f4(in[k], sx4, alpha); }
            }
        } else {
            if (alpha == 0) {
                for (int k = 0; k < height_limit; ++k) {
                    const int sx = sx4 + beta * (k - 3);
                    tmp[k]       = horizontal_filter_4x1_f1(in[k], sx);
                }
            } else {
                for (int k = 0; k < height_limit; ++k) {
                    const int sx = sx4 + beta * (k - 3);
                    tmp[k]       = horizontal_filter_4x1_f4(in[k], sx, alpha);
                }
            }
        }
    } else {
        if (beta == 0) {
            if (alpha == 0) {
                for (int k = 0; k < height_limit; ++k) { tmp[k] = horizontal_filter_8x1_f1(in[k], sx4); }
            } else {
                for (int k = 0; k < height_limit; ++k) { tmp[k] = horizontal_filter_8x1_f8(in[k], sx4, alpha); }
            }
        } else {
            if (alpha == 0) {
                for (int k = 0; k < height_limit; ++k) {
                    const int sx = sx4 + beta * (k - 3);
                    tmp[k]       = horizontal_filter_8x1_f1(in[k], sx);
                }
            } else {
                for (int k = 0; k < height_limit; ++k) {
                    const int sx = sx4 + beta * (k - 3);
                    tmp[k]       = horizontal_filter_8x1_f8(in[k], sx, alpha);
                }
            }
        }
    }
}

static INLINE void warp_affine_vertical(uint8_t *pred, int p_width, int p_height, int p_stride, int is_compound,
                                        uint16_t *dst, int dst_stride, int do_average, int use_dist_wtd_comp_avg,
                                        int16_t gamma, int16_t delta, const int64_t y4, const int i, const int j,
                                        int16x8_t tmp[], const int fwd, const int bwd) {
    const int bd                = 8;
    const int reduce_bits_horiz = ROUND0_BITS;
    const int offset_bits_vert  = bd + 2 * FILTER_BITS - reduce_bits_horiz;
    int       add_const_vert;
    if (is_compound) {
        add_const_vert = (1 << offset_bits_vert) + (1 << (COMPOUND_ROUND1_BITS - 1));
    } else {
        add_const_vert = (1 << offset_bits_vert) + (1 << (2 * FILTER_BITS - ROUND0_BITS - 1));
    }
    const int sub_constant = (1 << (bd - 1)) + (1 << bd);

    const int offset_bits   = bd + 2 * FILTER_BITS - ROUND0_BITS;
    const int res_sub_const = (1 << (2 * FILTER_BITS - ROUND0_BITS - COMPOUND_ROUND1_BITS - 1)) -
        (1 << (offset_bits - COMPOUND_ROUND1_BITS)) - (1 << (offset_bits - COMPOUND_ROUND1_BITS - 1));

    int32_t sy4 = y4 & ((1 << WARPEDMODEL_PREC_BITS) - 1);
    sy4 += gamma * (-4) + delta * (-4) + (1 << (WARPEDDIFF_PREC_BITS - 1)) +
        (WARPEDPIXEL_PREC_SHIFTS << WARPEDDIFF_PREC_BITS);
    sy4 &= ~((1 << WARP_PARAM_REDUCE_BITS) - 1);

    if (p_width > 4) {
        for (int k = -4; k < AOMMIN(4, p_height - i - 4); ++k) {
            int              sy    = sy4 + delta * (k + 4);
            const int16x8_t *v_src = tmp + (k + 4);

            int32x4_t res_lo, res_hi;
            if (gamma == 0) {
                vertical_filter_8x1_f1(v_src, &res_lo, &res_hi, sy);
            } else {
                vertical_filter_8x1_f8(v_src, &res_lo, &res_hi, sy, gamma);
            }

            res_lo = vaddq_s32(res_lo, vdupq_n_s32(add_const_vert));
            res_hi = vaddq_s32(res_hi, vdupq_n_s32(add_const_vert));

            if (is_compound) {
                uint16_t *const p       = (uint16_t *)&dst[(i + k + 4) * dst_stride + j];
                int16x8_t       res_s16 = vcombine_s16(vshrn_n_s32(res_lo, COMPOUND_ROUND1_BITS),
                                                 vshrn_n_s32(res_hi, COMPOUND_ROUND1_BITS));
                if (do_average) {
                    int16x8_t tmp16 = vreinterpretq_s16_u16(vld1q_u16(p));
                    if (use_dist_wtd_comp_avg) {
                        int32x4_t tmp32_lo = vmull_n_s16(vget_low_s16(tmp16), fwd);
                        int32x4_t tmp32_hi = vmull_n_s16(vget_high_s16(tmp16), fwd);
                        tmp32_lo           = vmlal_n_s16(tmp32_lo, vget_low_s16(res_s16), bwd);
                        tmp32_hi           = vmlal_n_s16(tmp32_hi, vget_high_s16(res_s16), bwd);
                        tmp16              = vcombine_s16(vshrn_n_s32(tmp32_lo, DIST_PRECISION_BITS),
                                             vshrn_n_s32(tmp32_hi, DIST_PRECISION_BITS));
                    } else {
                        tmp16 = vhaddq_s16(tmp16, res_s16);
                    }
                    int16x8_t res  = vaddq_s16(tmp16, vdupq_n_s16(res_sub_const));
                    uint8x8_t res8 = vqshrun_n_s16(res, 2 * FILTER_BITS - ROUND0_BITS - COMPOUND_ROUND1_BITS);
                    vst1_u8(&pred[(i + k + 4) * p_stride + j], res8);
                } else {
                    vst1q_u16(p, vreinterpretq_u16_s16(res_s16));
                }
            } else {
                int16x8_t res16 = vcombine_s16(vshrn_n_s32(res_lo, 2 * FILTER_BITS - ROUND0_BITS),
                                               vshrn_n_s32(res_hi, 2 * FILTER_BITS - ROUND0_BITS));
                res16           = vsubq_s16(res16, vdupq_n_s16(sub_constant));

                uint8_t *const p = (uint8_t *)&pred[(i + k + 4) * p_stride + j];
                vst1_u8(p, vqmovun_s16(res16));
            }
        }
    } else {
        // p_width == 4
        for (int k = -4; k < AOMMIN(4, p_height - i - 4); ++k) {
            int              sy    = sy4 + delta * (k + 4);
            const int16x8_t *v_src = tmp + (k + 4);

            int32x4_t res_lo;
            if (gamma == 0) {
                vertical_filter_4x1_f1(v_src, &res_lo, sy);
            } else {
                vertical_filter_4x1_f4(v_src, &res_lo, sy, gamma);
            }

            res_lo = vaddq_s32(res_lo, vdupq_n_s32(add_const_vert));

            if (is_compound) {
                uint16_t *const p = (uint16_t *)&dst[(i + k + 4) * dst_stride + j];

                int16x4_t res_lo_s16 = vshrn_n_s32(res_lo, COMPOUND_ROUND1_BITS);
                if (do_average) {
                    uint8_t *const dst8     = &pred[(i + k + 4) * p_stride + j];
                    int16x4_t      tmp16_lo = vreinterpret_s16_u16(vld1_u16(p));
                    if (use_dist_wtd_comp_avg) {
                        int32x4_t tmp32_lo = vmull_n_s16(tmp16_lo, fwd);
                        tmp32_lo           = vmlal_n_s16(tmp32_lo, res_lo_s16, bwd);
                        tmp16_lo           = vshrn_n_s32(tmp32_lo, DIST_PRECISION_BITS);
                    } else {
                        tmp16_lo = vhadd_s16(tmp16_lo, res_lo_s16);
                    }
                    int16x4_t res  = vadd_s16(tmp16_lo, vdup_n_s16(res_sub_const));
                    uint8x8_t res8 = vqshrun_n_s16(vcombine_s16(res, vdup_n_s16(0)),
                                                   2 * FILTER_BITS - ROUND0_BITS - COMPOUND_ROUND1_BITS);
                    vst1_lane_u32((uint32_t *)dst8, vreinterpret_u32_u8(res8), 0);
                } else {
                    uint16x4_t res_u16_low = vreinterpret_u16_s16(res_lo_s16);
                    vst1_u16(p, res_u16_low);
                }
            } else {
                int16x4_t res16 = vshrn_n_s32(res_lo, 2 * FILTER_BITS - ROUND0_BITS);
                res16           = vsub_s16(res16, vdup_n_s16(sub_constant));

                uint8_t *const p   = (uint8_t *)&pred[(i + k + 4) * p_stride + j];
                uint8x8_t      val = vqmovun_s16(vcombine_s16(res16, vdup_n_s16(0)));
                vst1_lane_u32((uint32_t *)p, vreinterpret_u32_u8(val), 0);
            }
        }
    }
}

void svt_av1_warp_affine_neon(const int32_t *mat, const uint8_t *ref, int width, int height, int stride, uint8_t *pred,
                              int p_col, int p_row, int p_width, int p_height, int p_stride, int subsampling_x,
                              int subsampling_y, ConvolveParams *conv_params, int16_t alpha, int16_t beta,
                              int16_t gamma, int16_t delta) {
    const int       w0                    = conv_params->fwd_offset;
    const int       w1                    = conv_params->bck_offset;
    const int       is_compound           = conv_params->is_compound;
    uint16_t *const dst                   = conv_params->dst;
    const int       dst_stride            = conv_params->dst_stride;
    const int       do_average            = conv_params->do_average;
    const int       use_dist_wtd_comp_avg = conv_params->use_dist_wtd_comp_avg;

    static const uint8_t k0To15[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    const uint8x16_t     indx_vec   = vld1q_u8(k0To15);

    assert(IMPLIES(is_compound, dst != NULL));
    assert(IMPLIES(do_average, is_compound));
    assert(6 == WARPEDPIXEL_PREC_BITS);

    for (int i = 0; i < p_height; i += 8) {
        for (int j = 0; j < p_width; j += 8) {
            const int32_t src_x = (p_col + j + 4) << subsampling_x;
            const int32_t src_y = (p_row + i + 4) << subsampling_y;
            const int64_t dst_x = (int64_t)mat[2] * src_x + (int64_t)mat[3] * src_y + (int64_t)mat[0];
            const int64_t dst_y = (int64_t)mat[4] * src_x + (int64_t)mat[5] * src_y + (int64_t)mat[1];

            const int64_t x4 = dst_x >> subsampling_x;
            const int64_t y4 = dst_y >> subsampling_y;

            int16x8_t tmp[15];
            warp_affine_horizontal(
                ref, width, height, stride, p_width, p_height, alpha, beta, x4, y4, i, tmp, indx_vec);
            warp_affine_vertical(pred,
                                 p_width,
                                 p_height,
                                 p_stride,
                                 is_compound,
                                 dst,
                                 dst_stride,
                                 do_average,
                                 use_dist_wtd_comp_avg,
                                 gamma,
                                 delta,
                                 y4,
                                 i,
                                 j,
                                 tmp,
                                 w0,
                                 w1);
        }
    }
}
