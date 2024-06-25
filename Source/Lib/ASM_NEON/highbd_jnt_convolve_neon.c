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

#include <arm_neon.h>
#include <assert.h>

#include "definitions.h"
#include "common_dsp_rtcd.h"
#include "convolve.h"

static INLINE int32x4_t svt_aom_convolve(const uint16x8_t *const s, const int16x8_t *const coeffs) {
    const int32x4_t res_0 = vpaddq_s32(vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[0])), vget_low_s16(coeffs[0])),
                                       vmull_s16(vget_high_s16(vreinterpretq_s16_u16(s[0])), vget_high_s16(coeffs[0])));
    const int32x4_t res_1 = vpaddq_s32(vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[1])), vget_low_s16(coeffs[1])),
                                       vmull_s16(vget_high_s16(vreinterpretq_s16_u16(s[1])), vget_high_s16(coeffs[1])));
    const int32x4_t res_2 = vpaddq_s32(vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[2])), vget_low_s16(coeffs[2])),
                                       vmull_s16(vget_high_s16(vreinterpretq_s16_u16(s[2])), vget_high_s16(coeffs[2])));
    const int32x4_t res_3 = vpaddq_s32(vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[3])), vget_low_s16(coeffs[3])),
                                       vmull_s16(vget_high_s16(vreinterpretq_s16_u16(s[3])), vget_high_s16(coeffs[3])));

    const int32x4_t res = vaddq_s32(vaddq_s32(res_0, res_1), vaddq_s32(res_2, res_3));

    return res;
}

static INLINE void prepare_coeffs(const int16_t *const filter, int16x8_t *const coeffs /* [4] */) {
    const int16x8_t coeff = vld1q_s16(filter);

    // coeffs 0 1 0 1 0 1 0 1
    coeffs[0] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(vreinterpretq_s32_s16(coeff), 0)));
    // coeffs 2 3 2 3 2 3 2 3
    coeffs[1] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(vreinterpretq_s32_s16(coeff), 1)));
    // coeffs 4 5 4 5 4 5 4 5
    coeffs[2] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(vreinterpretq_s32_s16(coeff), 2)));
    // coeffs 6 7 6 7 6 7 6 7
    coeffs[3] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(vreinterpretq_s32_s16(coeff), 3)));
}

static INLINE int32x4_t highbd_comp_avg_neon(const uint32x4_t *const data_ref_0, const uint32x4_t *const res_unsigned,
                                             const int32x4_t *const wt0, const int32x4_t *const wt1,
                                             const int use_dist_wtd_avg) {
    int32x4_t res;
    if (use_dist_wtd_avg) {
        const int32x4_t wt0_res = vmulq_s32(vreinterpretq_s32_u32(*data_ref_0), *wt0);
        const int32x4_t wt1_res = vmulq_s32(vreinterpretq_s32_u32(*res_unsigned), *wt1);

        const int32x4_t wt_res = vaddq_s32(wt0_res, wt1_res);
        res                    = vshrq_n_s32(wt_res, DIST_PRECISION_BITS);
    } else {
        const int32x4_t wt_res = vaddq_s32(vreinterpretq_s32_u32(*data_ref_0), vreinterpretq_s32_u32(*res_unsigned));
        res                    = vshrq_n_s32(wt_res, 1);
    }
    return res;
}

static INLINE int32x4_t highbd_convolve_rounding_neon(const int32x4_t *const res_unsigned,
                                                      const int32x4_t *const offset_const, const int round_shift) {
    const int32x4_t res_signed = vsubq_s32(*res_unsigned, *offset_const);

    return vrshlq_s32(res_signed, vdupq_n_s32(-round_shift));
}

void svt_av1_highbd_jnt_convolve_x_neon(const uint16_t *src, int32_t src_stride, uint16_t *dst0, int32_t dst_stride0,
                                        int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
                                        const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
                                        const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd) {
    if (w <= 4) {
        svt_av1_highbd_jnt_convolve_x_c(src,
                                        src_stride,
                                        dst0,
                                        dst_stride0,
                                        w,
                                        h,
                                        filter_params_x,
                                        filter_params_y,
                                        subpel_x_q4,
                                        subpel_y_q4,
                                        conv_params,
                                        bd);
        return;
    }
    CONV_BUF_TYPE        *dst        = conv_params->dst;
    int                   dst_stride = conv_params->dst_stride;
    const int             fo_horiz   = filter_params_x->taps / 2 - 1;
    const uint16_t *const src_ptr    = src - fo_horiz;

    int        i, j;
    uint16x8_t s[4];
    int16x8_t  coeffs_x[4];

    const int        do_average       = conv_params->do_average;
    const int        use_jnt_comp_avg = conv_params->use_jnt_comp_avg;
    const int        w0               = conv_params->fwd_offset;
    const int        w1               = conv_params->bck_offset;
    const int32x4_t  wt0              = vdupq_n_s32(w0);
    const int32x4_t  wt1              = vdupq_n_s32(w1);
    const int32x4_t  round_const_x    = vdupq_n_s32(((1 << conv_params->round_0) >> 1));
    const int        offset_0         = bd + 2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    const int        offset           = (1 << offset_0) + (1 << (offset_0 - 1));
    const int32x4_t  offset_const     = vdupq_n_s32(offset);
    const int        rounding_shift   = 2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    const uint16x8_t clip_pixel_to_bd = vdupq_n_u16(bd == 10 ? 1023 : (bd == 12 ? 4095 : 255));

    const int16_t *filter_x = av1_get_interp_filter_subpel_kernel(*filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    prepare_coeffs(filter_x, coeffs_x);

    for (j = 0; j < w; j += 8) {
        /* Horizontal filter */
        for (i = 0; i < h; i += 1) {
            const uint16x8_t row00 = vld1q_u16(&src_ptr[i * src_stride + j]);
            const uint16x8_t row01 = vld1q_u16(&src_ptr[i * src_stride + (j + 8)]);

            // even pixels
            s[0] = vextq_u16(row00, row01, 0);
            s[1] = vextq_u16(row00, row01, 2);
            s[2] = vextq_u16(row00, row01, 4);
            s[3] = vextq_u16(row00, row01, 6);

            int32x4_t res_even = svt_aom_convolve(s, coeffs_x);
            res_even           = vshlq_s32(vaddq_s32(res_even, round_const_x), vdupq_n_s32(-conv_params->round_0));

            // odd pixels
            s[0] = vextq_u16(row00, row01, 1);
            s[1] = vextq_u16(row00, row01, 3);
            s[2] = vextq_u16(row00, row01, 5);
            s[3] = vextq_u16(row00, row01, 7);

            int32x4_t res_odd = svt_aom_convolve(s, coeffs_x);
            res_odd           = vshlq_s32(vaddq_s32(res_odd, round_const_x), vdupq_n_s32(-conv_params->round_0));

            res_even = vreinterpretq_s32_u32(
                vshlq_u32(vreinterpretq_u32_s32(res_even), vdupq_n_s32(FILTER_BITS - conv_params->round_1)));
            res_odd = vreinterpretq_s32_u32(
                vshlq_u32(vreinterpretq_u32_s32(res_odd), vdupq_n_s32(FILTER_BITS - conv_params->round_1)));

            int32x4_t  res1            = vzip1q_s32(res_even, res_odd);
            uint32x4_t res_unsigned_lo = vreinterpretq_u32_s32(vaddq_s32(res1, offset_const));
            if (w - j < 8) {
                if (do_average) {
                    const uint16x4_t data_0     = vld1_u16(&dst[i * dst_stride + j]);
                    const uint32x4_t data_ref_0 = vmovl_u16(data_0);

                    const int32x4_t comp_avg_res = highbd_comp_avg_neon(
                        &data_ref_0, &res_unsigned_lo, &wt0, &wt1, use_jnt_comp_avg);
                    const int32x4_t round_result = highbd_convolve_rounding_neon(
                        &comp_avg_res, &offset_const, rounding_shift);

                    const uint16x8_t res_16b  = vcombine_u16(vqmovun_s32(round_result), vqmovun_s32(round_result));
                    const uint16x8_t res_clip = vminq_u16(res_16b, clip_pixel_to_bd);
                    vst1_u16(&dst0[i * dst_stride0 + j], vget_low_u16(res_clip));
                } else {
                    uint16x8_t res_16b = vcombine_u16(vqmovun_s32(vreinterpretq_s32_u32(res_unsigned_lo)),
                                                      vqmovun_s32(vreinterpretq_s32_u32(res_unsigned_lo)));
                    vst1_u16(&dst[i * dst_stride + j], vget_low_u16(res_16b));
                }
            } else {
                int32x4_t  res2            = vzip2q_s32(res_even, res_odd);
                uint32x4_t res_unsigned_hi = vreinterpretq_u32_s32(vaddq_s32(res2, offset_const));
                if (do_average) {
                    const uint16x8_t data_0        = vld1q_u16(&dst[i * dst_stride + j]);
                    const uint32x4_t data_ref_0_lo = vmovl_u16(vget_low_u16(data_0));
                    const uint32x4_t data_ref_0_hi = vmovl_u16(vget_high_u16(data_0));

                    const int32x4_t comp_avg_res_lo = highbd_comp_avg_neon(
                        &data_ref_0_lo, &res_unsigned_lo, &wt0, &wt1, use_jnt_comp_avg);
                    const int32x4_t comp_avg_res_hi = highbd_comp_avg_neon(
                        &data_ref_0_hi, &res_unsigned_hi, &wt0, &wt1, use_jnt_comp_avg);

                    const int32x4_t round_result_lo = highbd_convolve_rounding_neon(
                        &comp_avg_res_lo, &offset_const, rounding_shift);
                    const int32x4_t round_result_hi = highbd_convolve_rounding_neon(
                        &comp_avg_res_hi, &offset_const, rounding_shift);

                    const uint16x8_t res_16b = vcombine_u16(vqmovun_s32(round_result_lo), vqmovun_s32(round_result_hi));
                    const uint16x8_t res_clip = vminq_u16(res_16b, clip_pixel_to_bd);
                    vst1q_u16(&dst0[i * dst_stride0 + j], res_clip);
                } else {
                    uint16x8_t res_16b = vcombine_u16(vqmovun_s32(vreinterpretq_s32_u32(res_unsigned_lo)),
                                                      vqmovun_s32(vreinterpretq_s32_u32(res_unsigned_hi)));
                    vst1q_u16(&dst[i * dst_stride + j], res_16b);
                }
            }
        }
    }
}

void svt_av1_highbd_jnt_convolve_y_neon(const uint16_t *src, int32_t src_stride, uint16_t *dst0, int32_t dst_stride0,
                                        int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
                                        const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
                                        const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd) {
    if (w <= 4) {
        svt_av1_highbd_jnt_convolve_y_c(src,
                                        src_stride,
                                        dst0,
                                        dst_stride0,
                                        w,
                                        h,
                                        filter_params_x,
                                        filter_params_y,
                                        subpel_x_q4,
                                        subpel_y_q4,
                                        conv_params,
                                        bd);
        return;
    }
    CONV_BUF_TYPE        *dst        = conv_params->dst;
    const int             dst_stride = conv_params->dst_stride;
    const int             fo_vert    = filter_params_y->taps / 2 - 1;
    const uint16_t *const src_ptr    = src - fo_vert * src_stride;

    int       i, j;
    const int do_average       = conv_params->do_average;
    const int use_jnt_comp_avg = conv_params->use_jnt_comp_avg;

    const int       w0            = conv_params->fwd_offset;
    const int       w1            = conv_params->bck_offset;
    const int32x4_t wt0           = vdupq_n_s32(w0);
    const int32x4_t wt1           = vdupq_n_s32(w1);
    const int32x4_t round_const_y = vdupq_n_s32(((1 << conv_params->round_1) >> 1));

    const int        offset_0         = bd + 2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    const int        offset           = (1 << offset_0) + (1 << (offset_0 - 1));
    const int32x4_t  offset_const     = vdupq_n_s32(offset);
    const int        rounding_shift   = 2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    const uint16x8_t clip_pixel_to_bd = vdupq_n_u16(bd == 10 ? 1023 : (bd == 12 ? 4095 : 255));
    uint16x8_t       s[16];
    int16x8_t        coeffs_y[4];

    const int16_t *filter_y = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_q4 & SUBPEL_MASK);
    prepare_coeffs(filter_y, coeffs_y);

    for (j = 0; j < w; j += 8) {
        const uint16_t *data = &src_ptr[j];
        /* Vertical filter */
        {
            uint16x8_t s0 = vld1q_u16(data + 0 * src_stride);
            uint16x8_t s1 = vld1q_u16(data + 1 * src_stride);
            uint16x8_t s2 = vld1q_u16(data + 2 * src_stride);
            uint16x8_t s3 = vld1q_u16(data + 3 * src_stride);
            uint16x8_t s4 = vld1q_u16(data + 4 * src_stride);
            uint16x8_t s5 = vld1q_u16(data + 5 * src_stride);
            uint16x8_t s6 = vld1q_u16(data + 6 * src_stride);

            s[0] = vzip1q_u16(s0, s1);
            s[1] = vzip1q_u16(s2, s3);
            s[2] = vzip1q_u16(s4, s5);

            s[4] = vzip2q_u16(s0, s1);
            s[5] = vzip2q_u16(s2, s3);
            s[6] = vzip2q_u16(s4, s5);

            s[0 + 8] = vzip1q_u16(s1, s2);
            s[1 + 8] = vzip1q_u16(s3, s4);
            s[2 + 8] = vzip1q_u16(s5, s6);

            s[4 + 8] = vzip2q_u16(s1, s2);
            s[5 + 8] = vzip2q_u16(s3, s4);
            s[6 + 8] = vzip2q_u16(s5, s6);

            for (i = 0; i < h; i += 2) {
                data = &src_ptr[i * src_stride + j];

                const uint16x8_t s7 = vld1q_u16(data + 7 * src_stride);
                const uint16x8_t s8 = vld1q_u16(data + 8 * src_stride);

                s[3] = vzip1q_u16(s6, s7);
                s[7] = vzip2q_u16(s6, s7);

                s[3 + 8] = vzip1q_u16(s7, s8);
                s[7 + 8] = vzip2q_u16(s7, s8);

                const int32x4_t res_a0       = svt_aom_convolve(s, coeffs_y);
                int32x4_t       res_a_round0 = vshlq_s32(res_a0, vdupq_n_s32(FILTER_BITS - conv_params->round_0));
                res_a_round0 = vshlq_s32(vaddq_s32(res_a_round0, round_const_y), vdupq_n_s32(-conv_params->round_1));

                const int32x4_t res_a1       = svt_aom_convolve(s + 8, coeffs_y);
                int32x4_t       res_a_round1 = vshlq_s32(res_a1, vdupq_n_s32(FILTER_BITS - conv_params->round_0));
                res_a_round1 = vshlq_s32(vaddq_s32(res_a_round1, round_const_y), vdupq_n_s32(-conv_params->round_1));

                const uint32x4_t res_unsigned_lo_0 = vreinterpretq_u32_s32(vaddq_s32(res_a_round0, offset_const));
                const uint32x4_t res_unsigned_lo_1 = vreinterpretq_u32_s32(vaddq_s32(res_a_round1, offset_const));

                if (w - j < 8) {
                    if (do_average) {
                        const uint16x4_t data_0 = vld1_u16(&dst[i * dst_stride + j]);
                        const uint16x4_t data_1 = vld1_u16(&dst[i * dst_stride + j + dst_stride]);

                        const uint32x4_t data_ref_0 = vmovl_u16(data_0);
                        const uint32x4_t data_ref_1 = vmovl_u16(data_1);

                        const int32x4_t comp_avg_res_0 = highbd_comp_avg_neon(
                            &data_ref_0, &res_unsigned_lo_0, &wt0, &wt1, use_jnt_comp_avg);
                        const int32x4_t comp_avg_res_1 = highbd_comp_avg_neon(
                            &data_ref_1, &res_unsigned_lo_1, &wt0, &wt1, use_jnt_comp_avg);

                        const int32x4_t round_result_0 = highbd_convolve_rounding_neon(
                            &comp_avg_res_0, &offset_const, rounding_shift);
                        const int32x4_t round_result_1 = highbd_convolve_rounding_neon(
                            &comp_avg_res_1, &offset_const, rounding_shift);

                        const uint16x8_t res_16b_0  = vqmovun_high_s32(vqmovun_s32(round_result_0), round_result_0);
                        const uint16x8_t res_clip_0 = vminq_u16(res_16b_0, clip_pixel_to_bd);
                        const uint16x8_t res_16b_1  = vqmovun_high_s32(vqmovun_s32(round_result_1), round_result_1);
                        const uint16x8_t res_clip_1 = vminq_u16(res_16b_1, clip_pixel_to_bd);

                        vst1_u16(&dst0[i * dst_stride0 + j], vget_low_u16(res_clip_0));
                        vst1_u16(&dst0[i * dst_stride0 + j + dst_stride0], vget_low_u16(res_clip_1));

                    } else {
                        const uint16x8_t res_16b_0 = vqmovun_high_s32(
                            vqmovun_s32(vreinterpretq_s32_u32(res_unsigned_lo_0)),
                            vreinterpretq_s32_u32(res_unsigned_lo_0));

                        const uint16x8_t res_16b_1 = vqmovun_high_s32(
                            vqmovun_s32(vreinterpretq_s32_u32(res_unsigned_lo_1)),
                            vreinterpretq_s32_u32(res_unsigned_lo_1));

                        vst1_u16(&dst[i * dst_stride + j], vget_low_u16(res_16b_0));
                        vst1_u16(&dst[i * dst_stride + j + dst_stride], vget_low_u16(res_16b_1));
                    }
                } else {
                    const int32x4_t res_b0       = svt_aom_convolve(s + 4, coeffs_y);
                    int32x4_t       res_b_round0 = vshlq_s32(res_b0, vdupq_n_s32(FILTER_BITS - conv_params->round_0));
                    res_b_round0                 = vshlq_s32(vaddq_s32(res_b_round0, round_const_y),
                                             vdupq_n_s32(-conv_params->round_1));

                    const int32x4_t res_b1       = svt_aom_convolve(s + 4 + 8, coeffs_y);
                    int32x4_t       res_b_round1 = vshlq_s32(res_b1, vdupq_n_s32(FILTER_BITS - conv_params->round_0));
                    res_b_round1                 = vshlq_s32(vaddq_s32(res_b_round1, round_const_y),
                                             vdupq_n_s32(-conv_params->round_1));

                    uint32x4_t res_unsigned_hi_0 = vreinterpretq_u32_s32(vaddq_s32(res_b_round0, offset_const));
                    uint32x4_t res_unsigned_hi_1 = vreinterpretq_u32_s32(vaddq_s32(res_b_round1, offset_const));

                    if (do_average) {
                        const uint16x8_t data_0 = vld1q_u16(&dst[i * dst_stride + j]);
                        const uint16x8_t data_1 = vld1q_u16(&dst[i * dst_stride + j + dst_stride]);

                        const uint32x4_t data_ref_0_lo_0 = vmovl_u16(vget_low_u16(data_0));
                        const uint32x4_t data_ref_0_lo_1 = vmovl_u16(vget_low_u16(data_1));

                        const uint32x4_t data_ref_0_hi_0 = vmovl_high_u16(data_0);
                        const uint32x4_t data_ref_0_hi_1 = vmovl_high_u16(data_1);

                        const int32x4_t comp_avg_res_lo_0 = highbd_comp_avg_neon(
                            &data_ref_0_lo_0, &res_unsigned_lo_0, &wt0, &wt1, use_jnt_comp_avg);
                        const int32x4_t comp_avg_res_lo_1 = highbd_comp_avg_neon(
                            &data_ref_0_lo_1, &res_unsigned_lo_1, &wt0, &wt1, use_jnt_comp_avg);
                        const int32x4_t comp_avg_res_hi_0 = highbd_comp_avg_neon(
                            &data_ref_0_hi_0, &res_unsigned_hi_0, &wt0, &wt1, use_jnt_comp_avg);
                        const int32x4_t comp_avg_res_hi_1 = highbd_comp_avg_neon(
                            &data_ref_0_hi_1, &res_unsigned_hi_1, &wt0, &wt1, use_jnt_comp_avg);

                        const int32x4_t round_result_lo_0 = highbd_convolve_rounding_neon(
                            &comp_avg_res_lo_0, &offset_const, rounding_shift);
                        const int32x4_t round_result_lo_1 = highbd_convolve_rounding_neon(
                            &comp_avg_res_lo_1, &offset_const, rounding_shift);
                        const int32x4_t round_result_hi_0 = highbd_convolve_rounding_neon(
                            &comp_avg_res_hi_0, &offset_const, rounding_shift);
                        const int32x4_t round_result_hi_1 = highbd_convolve_rounding_neon(
                            &comp_avg_res_hi_1, &offset_const, rounding_shift);

                        const uint16x8_t res_16b_0  = vqmovun_high_s32(vqmovun_s32(round_result_lo_0),
                                                                      round_result_hi_0);
                        const uint16x8_t res_clip_0 = vminq_u16(res_16b_0, clip_pixel_to_bd);

                        const uint16x8_t res_16b_1  = vqmovun_high_s32(vqmovun_s32(round_result_lo_1),
                                                                      round_result_hi_1);
                        const uint16x8_t res_clip_1 = vminq_u16(res_16b_1, clip_pixel_to_bd);

                        vst1q_u16(&dst0[i * dst_stride0 + j], res_clip_0);
                        vst1q_u16(&dst0[i * dst_stride0 + j + dst_stride0], res_clip_1);
                    } else {
                        const uint16x8_t res_16bit0 = vqmovun_high_s32(
                            vqmovun_s32(vreinterpretq_s32_u32(res_unsigned_lo_0)),
                            vreinterpretq_s32_u32(res_unsigned_hi_0));
                        const uint16x8_t res_16bit1 = vqmovun_high_s32(
                            vqmovun_s32(vreinterpretq_s32_u32(res_unsigned_lo_1)),
                            vreinterpretq_s32_u32(res_unsigned_hi_1));
                        vst1q_u16(&dst[i * dst_stride + j], res_16bit0);
                        vst1q_u16(&dst[i * dst_stride + j + dst_stride], res_16bit1);
                    }
                }
                s[0] = s[1];
                s[1] = s[2];
                s[2] = s[3];

                s[4] = s[5];
                s[5] = s[6];
                s[6] = s[7];

                s[0 + 8] = s[1 + 8];
                s[1 + 8] = s[2 + 8];
                s[2 + 8] = s[3 + 8];

                s[4 + 8] = s[5 + 8];
                s[5 + 8] = s[6 + 8];
                s[6 + 8] = s[7 + 8];

                s6 = s8;
            }
        }
    }
}
