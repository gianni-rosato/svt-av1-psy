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

#include <assert.h>

#include <arm_neon.h>

#include "definitions.h"
#include "common_dsp_rtcd.h"
#include "convolve.h"

static INLINE int32x4_t highbd_comp_avg_neon(const int32x4_t *const data_ref_0, const int32x4_t *const res_unsigned,
                                             const int32x4_t *const wt0, const int32x4_t *const wt1,
                                             const int use_dist_wtd_avg) {
    int32x4_t res;
    if (use_dist_wtd_avg) {
        const int32x4_t wt0_res = vmulq_s32(*data_ref_0, *wt0);
        const int32x4_t wt1_res = vmulq_s32(*res_unsigned, *wt1);

        const int32x4_t wt_res = vaddq_s32(wt0_res, wt1_res);
        res                    = vshrq_n_s32(wt_res, DIST_PRECISION_BITS);
    } else {
        const int32x4_t wt_res = vaddq_s32(*data_ref_0, *res_unsigned);
        res                    = vshrq_n_s32(wt_res, 1);
    }
    return res;
}

static INLINE int32x4_t highbd_convolve_rounding_neon(const int32x4_t *const res_unsigned,
                                                      const int32x4_t *const offset_const,
                                                      const int32x4_t *const round_const, const int round_shift) {
    const int32x4_t res_signed = vsubq_s32(*res_unsigned, *offset_const);
    const int32x4_t res_round  = vshlq_s32(vaddq_s32(res_signed, *round_const), vdupq_n_s32(-round_shift));

    return res_round;
}

static INLINE int32x4_t multiply_then_pairwise_add(const int16x8_t a, const int16x8_t b) {
    int32x4_t a_even = vmovl_s16(vget_low_s16(vuzp1q_s16(a, a)));
    int32x4_t a_odd  = vmovl_s16(vget_low_s16(vuzp2q_s16(a, a)));

    int32x4_t b_even = vmovl_s16(vget_low_s16(vuzp1q_s16(b, b)));
    int32x4_t b_odd  = vmovl_s16(vget_low_s16(vuzp2q_s16(b, b)));

    int32x4_t res_even = vmulq_s32(a_even, b_even);
    int32x4_t res_odd  = vmulq_s32(a_odd, b_odd);

    return vaddq_s32(res_even, res_odd);
}

void svt_av1_highbd_jnt_convolve_2d_neon(const uint16_t *src, int32_t src_stride, uint16_t *dst16, int32_t dst16_stride,
                                         int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
                                         const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
                                         const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd) {
    if (w < 8) {
        svt_av1_highbd_jnt_convolve_2d_c(src,
                                         src_stride,
                                         dst16,
                                         dst16_stride,
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

    DECLARE_ALIGNED(16, int16_t, im_block[(MAX_SB_SIZE + MAX_FILTER_TAP - 1) * MAX_SB_SIZE]);
    ConvBufType          *dst        = conv_params->dst;
    int                   dst_stride = conv_params->dst_stride;
    int                   im_h       = h + filter_params_y->taps - 1;
    int                   im_stride  = MAX_SB_SIZE;
    int                   i, j;
    const int             do_average       = conv_params->do_average;
    const int             use_jnt_comp_avg = conv_params->use_jnt_comp_avg;
    const int             fo_vert          = filter_params_y->taps / 2 - 1;
    const int             fo_horiz         = filter_params_x->taps / 2 - 1;
    const uint16_t *const src_ptr          = src - fo_vert * src_stride - fo_horiz;

    const int w0 = conv_params->fwd_offset;
    const int w1 = conv_params->bck_offset;

    const int offset_0       = bd + 2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    const int offset         = (1 << offset_0) + (1 << (offset_0 - 1));
    const int rounding_shift = 2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;

    const int32x4_t wt0 = vdupq_n_s32(w0);
    const int32x4_t wt1 = vdupq_n_s32(w1);

    const int32x4_t offset_const   = vdupq_n_s32(offset);
    const int32x4_t rounding_const = vdupq_n_s32((1 << rounding_shift) >> 1);

    const uint16x8_t clip_pixel_to_bd_128 = vdupq_n_u16(bd == 10 ? 1023 : (bd == 12 ? 4095 : 255));

    // Check that, even with 12-bit input, the intermediate values will fit
    // into an unsigned 16-bit intermediate array.
    assert(bd + FILTER_BITS + 2 - conv_params->round_0 <= 16);

    /* Horizontal filter */
    {
        const int16_t  *x_filter = av1_get_interp_filter_subpel_kernel(*filter_params_x, subpel_x_q4 & SUBPEL_MASK);
        const int16x8_t coeffs_x = vld1q_s16(x_filter);

        // coeffs 0 1 0 1 2 3 2 3
        const int16x8_t tmp_0 = vreinterpretq_s16_s32(
            vzip1q_s32(vreinterpretq_s32_s16(coeffs_x), vreinterpretq_s32_s16(coeffs_x)));
        // coeffs 4 5 4 5 6 7 6 7
        const int16x8_t tmp_1 = vreinterpretq_s16_s32(
            vzip2q_s32(vreinterpretq_s32_s16(coeffs_x), vreinterpretq_s32_s16(coeffs_x)));

        // coeffs 0 1 0 1 0 1 0 1
        const int16x8_t coeff_01 = vreinterpretq_s16_s64(
            vzip1q_s64(vreinterpretq_s64_s16(tmp_0), vreinterpretq_s64_s16(tmp_0)));
        // coeffs 2 3 2 3 2 3 2 3
        const int16x8_t coeff_23 = vreinterpretq_s16_s64(
            vzip2q_s64(vreinterpretq_s64_s16(tmp_0), vreinterpretq_s64_s16(tmp_0)));
        // coeffs 4 5 4 5 4 5 4 5
        const int16x8_t coeff_45 = vreinterpretq_s16_s64(
            vzip1q_s64(vreinterpretq_s64_s16(tmp_1), vreinterpretq_s64_s16(tmp_1)));
        // coeffs 6 7 6 7 6 7 6 7
        const int16x8_t coeff_67 = vreinterpretq_s16_s64(
            vzip2q_s64(vreinterpretq_s64_s16(tmp_1), vreinterpretq_s64_s16(tmp_1)));

        const int32x4_t round_const  = vdupq_n_s32(((1 << conv_params->round_0) >> 1) + (1 << (bd + FILTER_BITS - 1)));
        const int32_t   round_shift  = conv_params->round_0;
        const int32x4_t vround_shift = vdupq_n_s32(round_shift);

        for (i = 0; i < im_h; ++i) {
            for (j = 0; j < w; j += 8) {
                const int16x8_t data  = vld1q_s16((int16_t *)&src_ptr[i * src_stride + j]);
                const int16x8_t data2 = vld1q_s16((int16_t *)&src_ptr[i * src_stride + j + 8]);

                // Filter even-index pixels
                const int32x4_t res_0 = multiply_then_pairwise_add(data, coeff_01);
                const int32x4_t res_2 = multiply_then_pairwise_add(vextq_s16(data, data2, 4 / 2), coeff_23);
                const int32x4_t res_4 = multiply_then_pairwise_add(vextq_s16(data, data2, 8 / 2), coeff_45);
                const int32x4_t res_6 = multiply_then_pairwise_add(vextq_s16(data, data2, 12 / 2), coeff_67);

                int32x4_t res_even = vaddq_s32(vaddq_s32(res_0, res_4), vaddq_s32(res_2, res_6));
                res_even           = vshlq_s32(vaddq_s32(res_even, round_const), -vround_shift);

                // Filter odd-index pixels
                const int32x4_t res_1 = multiply_then_pairwise_add(vextq_s16(data, data2, 2 / 2), coeff_01);
                const int32x4_t res_3 = multiply_then_pairwise_add(vextq_s16(data, data2, 6 / 2), coeff_23);
                const int32x4_t res_5 = multiply_then_pairwise_add(vextq_s16(data, data2, 10 / 2), coeff_45);
                const int32x4_t res_7 = multiply_then_pairwise_add(vextq_s16(data, data2, 14 / 2), coeff_67);

                int32x4_t res_odd = vaddq_s32(vaddq_s32(res_1, res_5), vaddq_s32(res_3, res_7));
                res_odd           = vshlq_s32(vaddq_s32(res_odd, round_const), -vround_shift);

                // Pack in the column order 0, 2, 4, 6, 1, 3, 5, 7
                int16x8_t res = vcombine_s16(vqmovn_s32(res_even), vqmovn_s32(res_odd));
                vst1q_s16(&im_block[i * im_stride + j], res);
            }
        }
    }

    /* Vertical filter */
    {
        const int16_t  *y_filter = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_q4 & SUBPEL_MASK);
        const int16x8_t coeffs_y = vld1q_s16(y_filter);

        // coeffs 0 1 0 1 2 3 2 3
        const int16x8_t tmp_0 = vreinterpretq_s16_s32(
            vzip1q_s32(vreinterpretq_s32_s16(coeffs_y), vreinterpretq_s32_s16(coeffs_y)));
        // coeffs 4 5 4 5 6 7 6 7
        const int16x8_t tmp_1 = vreinterpretq_s16_s32(
            vzip2q_s32(vreinterpretq_s32_s16(coeffs_y), vreinterpretq_s32_s16(coeffs_y)));

        // coeffs 0 1 0 1 0 1 0 1
        const int16x8_t coeff_01 = vreinterpretq_s16_s64(
            vzip1q_s64(vreinterpretq_s64_s16(tmp_0), vreinterpretq_s64_s16(tmp_0)));
        // coeffs 2 3 2 3 2 3 2 3
        const int16x8_t coeff_23 = vreinterpretq_s16_s64(
            vzip2q_s64(vreinterpretq_s64_s16(tmp_0), vreinterpretq_s64_s16(tmp_0)));
        // coeffs 4 5 4 5 4 5 4 5
        const int16x8_t coeff_45 = vreinterpretq_s16_s64(
            vzip1q_s64(vreinterpretq_s64_s16(tmp_1), vreinterpretq_s64_s16(tmp_1)));
        // coeffs 6 7 6 7 6 7 6 7
        const int16x8_t coeff_67 = vreinterpretq_s16_s64(
            vzip2q_s64(vreinterpretq_s64_s16(tmp_1), vreinterpretq_s64_s16(tmp_1)));

        const int32x4_t round_const  = vdupq_n_s32(((1 << conv_params->round_1) >> 1) -
                                                  (1 << (bd + 2 * FILTER_BITS - conv_params->round_0 - 1)));
        const int32_t   round_shift  = conv_params->round_1;
        const int32x4_t vround_shift = vdupq_n_s32(round_shift);

        for (i = 0; i < h; ++i) {
            for (j = 0; j < w; j += 8) {
                const int16_t *data = &im_block[i * im_stride + j];

                // Filter even-index pixels
                const int16x8_t data_0 = vld1q_s16(data + 0 * im_stride);
                const int16x8_t data_1 = vld1q_s16(data + 1 * im_stride);
                const int16x8_t data_2 = vld1q_s16(data + 2 * im_stride);
                const int16x8_t data_3 = vld1q_s16(data + 3 * im_stride);
                const int16x8_t data_4 = vld1q_s16(data + 4 * im_stride);
                const int16x8_t data_5 = vld1q_s16(data + 5 * im_stride);
                const int16x8_t data_6 = vld1q_s16(data + 6 * im_stride);
                const int16x8_t data_7 = vld1q_s16(data + 7 * im_stride);

                const int16x8_t src_0 = vzip1q_s16(data_0, data_1);
                const int16x8_t src_2 = vzip1q_s16(data_2, data_3);
                const int16x8_t src_4 = vzip1q_s16(data_4, data_5);
                const int16x8_t src_6 = vzip1q_s16(data_6, data_7);

                const int32x4_t res_0 = multiply_then_pairwise_add(src_0, coeff_01);
                const int32x4_t res_2 = multiply_then_pairwise_add(src_2, coeff_23);
                const int32x4_t res_4 = multiply_then_pairwise_add(src_4, coeff_45);
                const int32x4_t res_6 = multiply_then_pairwise_add(src_6, coeff_67);

                const int32x4_t res_even = vaddq_s32(vaddq_s32(res_0, res_2), vaddq_s32(res_4, res_6));

                // Filter odd-index pixels
                const int16x8_t src_1 = vzip2q_s16(data_0, data_1);
                const int16x8_t src_3 = vzip2q_s16(data_2, data_3);
                const int16x8_t src_5 = vzip2q_s16(data_4, data_5);
                const int16x8_t src_7 = vzip2q_s16(data_6, data_7);

                const int32x4_t res_1 = multiply_then_pairwise_add(src_1, coeff_01);
                const int32x4_t res_3 = multiply_then_pairwise_add(src_3, coeff_23);
                const int32x4_t res_5 = multiply_then_pairwise_add(src_5, coeff_45);
                const int32x4_t res_7 = multiply_then_pairwise_add(src_7, coeff_67);

                const int32x4_t res_odd = vaddq_s32(vaddq_s32(res_1, res_3), vaddq_s32(res_5, res_7));

                // Rearrange pixels back into the order 0 ... 7
                const int32x4_t res_lo = vzip1q_s32(res_even, res_odd);
                const int32x4_t res_hi = vzip2q_s32(res_even, res_odd);

                const int32x4_t res_lo_round    = vshlq_s32(vaddq_s32(res_lo, round_const), -vround_shift);
                const int32x4_t res_unsigned_lo = vaddq_s32(res_lo_round, offset_const);

                const int32x4_t res_hi_round    = vshlq_s32(vaddq_s32(res_hi, round_const), -vround_shift);
                const int32x4_t res_unsigned_hi = vaddq_s32(res_hi_round, offset_const);

                if (do_average) {
                    const uint16x4_t data_lo = vld1_u16(&dst[i * dst_stride + j + 0]);
                    const uint16x4_t data_hi = vld1_u16(&dst[i * dst_stride + j + 4]);

                    const int32x4_t data_ref_0_lo = vreinterpretq_s32_u32(vmovl_u16(data_lo));
                    const int32x4_t data_ref_0_hi = vreinterpretq_s32_u32(vmovl_u16(data_hi));

                    const int32x4_t comp_avg_res_lo = highbd_comp_avg_neon(
                        &data_ref_0_lo, &res_unsigned_lo, &wt0, &wt1, use_jnt_comp_avg);
                    const int32x4_t comp_avg_res_hi = highbd_comp_avg_neon(
                        &data_ref_0_hi, &res_unsigned_hi, &wt0, &wt1, use_jnt_comp_avg);

                    const int32x4_t round_result_lo = highbd_convolve_rounding_neon(
                        &comp_avg_res_lo, &offset_const, &rounding_const, rounding_shift);
                    const int32x4_t round_result_hi = highbd_convolve_rounding_neon(
                        &comp_avg_res_hi, &offset_const, &rounding_const, rounding_shift);

                    const uint16x8_t res_16b = vcombine_u16(vqmovun_s32(round_result_lo), vqmovun_s32(round_result_hi));
                    const uint16x8_t res_clip = vminq_u16(res_16b, clip_pixel_to_bd_128);

                    vst1q_u16(&dst16[i * dst16_stride + j], res_clip);
                } else {
                    const uint16x8_t res_16b = vcombine_u16(vqmovun_s32(res_unsigned_lo), vqmovun_s32(res_unsigned_hi));

                    vst1q_u16(&dst[i * dst_stride + j], res_16b);
                }
            }
        }
    }
}

static INLINE void svt_prepare_coeffs_12tap(const int16_t *const filter, int16x8_t *coeffs /* [6] */) {
    int32x4_t coeffs_y  = vld1q_s32((int32_t const *)filter);
    int32x4_t coeffs_y2 = vld1q_s32((int32_t const *)(filter + 8));

    coeffs[0] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(coeffs_y, 0))); // coeffs 0 1 0 1 0 1 0 1
    coeffs[1] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(coeffs_y, 1))); // coeffs 2 3 2 3 2 3 2 3
    coeffs[2] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(coeffs_y, 2))); // coeffs 4 5 4 5 4 5 4 5
    coeffs[3] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(coeffs_y, 3))); // coeffs 6 7 6 7 6 7 6 7

    coeffs[4] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(coeffs_y2, 0))); // coeffs 8 9 8 9 8 9 8 9
    coeffs[5] = vreinterpretq_s16_s32(vdupq_n_s32(vgetq_lane_s32(coeffs_y2, 1))); // coeffs 10 11 10 11 10 11 10 11
}

static INLINE int32x4_t convolve_12tap(const int16x8_t *s, const int16x8_t *coeffs) {
    const int32x4_t d0     = vmull_s16(vget_low_s16(s[0]), vget_low_s16(coeffs[0]));
    const int32x4_t d1     = vmull_s16(vget_low_s16(s[1]), vget_low_s16(coeffs[1]));
    const int32x4_t d2     = vmull_s16(vget_low_s16(s[2]), vget_low_s16(coeffs[2]));
    const int32x4_t d3     = vmull_s16(vget_low_s16(s[3]), vget_low_s16(coeffs[3]));
    const int32x4_t d4     = vmull_s16(vget_low_s16(s[4]), vget_low_s16(coeffs[4]));
    const int32x4_t d5     = vmull_s16(vget_low_s16(s[5]), vget_low_s16(coeffs[5]));
    const int32x4_t d_0123 = vaddq_s32(vaddq_s32(d0, d1), vaddq_s32(d2, d3));
    const int32x4_t d      = vaddq_s32(vaddq_s32(d4, d5), d_0123);
    return d;
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

static INLINE int32x4_t svt_aom_convolve(const int16x8_t *const s, const int16x8_t *const coeffs) {
    const int32x4_t res_0 = vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(s[0])), vmovl_s16(vget_low_s16(coeffs[0]))),
                                       vmulq_s32(vmovl_s16(vget_high_s16(s[0])), vmovl_s16(vget_high_s16(coeffs[0]))));
    const int32x4_t res_1 = vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(s[1])), vmovl_s16(vget_low_s16(coeffs[1]))),
                                       vmulq_s32(vmovl_s16(vget_high_s16(s[1])), vmovl_s16(vget_high_s16(coeffs[1]))));
    const int32x4_t res_2 = vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(s[2])), vmovl_s16(vget_low_s16(coeffs[2]))),
                                       vmulq_s32(vmovl_s16(vget_high_s16(s[2])), vmovl_s16(vget_high_s16(coeffs[2]))));
    const int32x4_t res_3 = vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(s[3])), vmovl_s16(vget_low_s16(coeffs[3]))),
                                       vmulq_s32(vmovl_s16(vget_high_s16(s[3])), vmovl_s16(vget_high_s16(coeffs[3]))));

    const int32x4_t res = vaddq_s32(vaddq_s32(res_0, res_1), vaddq_s32(res_2, res_3));

    return res;
}

void svt_av1_highbd_convolve_2d_sr_neon(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
                                        int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
                                        const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
                                        const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd) {
    DECLARE_ALIGNED(32, int16_t, im_block[(MAX_SB_SIZE + MAX_FILTER_TAP) * 8]);
    int                   im_h      = h + filter_params_y->taps - 1;
    int                   im_stride = 8;
    int                   i, j;
    const int             fo_vert  = filter_params_y->taps / 2 - 1;
    const int             fo_horiz = filter_params_x->taps / 2 - 1;
    const uint16_t *const src_ptr  = src - fo_vert * src_stride - fo_horiz;

    // Check that, even with 12-bit input, the intermediate values will fit
    // into an unsigned 16-bit intermediate array.
    assert(bd + FILTER_BITS + 2 - conv_params->round_0 <= 16);

    const int32x4_t round_const_x = vdupq_n_s32(((1 << conv_params->round_0) >> 1) + (1 << (bd + FILTER_BITS - 1)));
    const int32x4_t round_shift_x = vsetq_lane_s32(conv_params->round_0, vdupq_n_s32(0), 0);

    const int32x4_t round_const_y = vdupq_n_s32(((1 << conv_params->round_1) >> 1) -
                                                (1 << (bd + 2 * FILTER_BITS - conv_params->round_0 - 1)));
    const int32x4_t round_shift_y = vsetq_lane_s32(conv_params->round_1, vdupq_n_s32(0), 0);

    const int       bits             = FILTER_BITS * 2 - conv_params->round_0 - conv_params->round_1;
    const int32x4_t round_shift_bits = vsetq_lane_s32(bits, vdupq_n_s32(0), 0);
    const int32x4_t round_const_bits = vdupq_n_s32((1 << bits) >> 1);
    const int16x8_t clip_pixel       = vdupq_n_s16(bd == 10 ? 1023 : (bd == 12 ? 4095 : 255));
    const int16x8_t zero             = vdupq_n_s16(0);

    const int16_t *const x_filter = av1_get_interp_filter_subpel_kernel(*filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    const int16_t *const y_filter = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_q4 & SUBPEL_MASK);

    if (filter_params_x->taps == 12) {
        int16x8_t coeffs_x[6], coeffs_y[6], s[24];
        svt_prepare_coeffs_12tap(x_filter, coeffs_x);
        svt_prepare_coeffs_12tap(y_filter, coeffs_y);

        for (j = 0; j < w; j += 8) {
            /* Horizontal filter */
            {
                for (i = 0; i < im_h; i += 1) {
                    const int16x8_t row00 = vld1q_s16((const int16_t *)&src_ptr[i * src_stride + j]);
                    const int16x8_t row01 = vld1q_s16((const int16_t *)&src_ptr[i * src_stride + (j + 8)]);
                    const int16x8_t row02 = vld1q_s16((const int16_t *)&src_ptr[i * src_stride + (j + 16)]);

                    // even pixels
                    s[0] = vextq_s16(row00, row01, 0);
                    s[1] = vextq_s16(row00, row01, 2);
                    s[2] = vextq_s16(row00, row01, 4);
                    s[3] = vextq_s16(row00, row01, 6);
                    s[4] = vextq_s16(row01, row02, 0);
                    s[5] = vextq_s16(row01, row02, 2);

                    int32x4_t res_even = convolve_12tap(s, coeffs_x);

                    res_even = vshlq_s32(vaddq_s32(res_even, round_const_x), -round_shift_x);

                    // odd pixels
                    s[0] = vextq_s16(row00, row01, 1);
                    s[1] = vextq_s16(row00, row01, 4);
                    s[2] = vextq_s16(row00, row01, 5);
                    s[3] = vextq_s16(row00, row01, 7);
                    s[4] = vextq_s16(row01, row02, 1);
                    s[5] = vextq_s16(row01, row02, 3);

                    int32x4_t res_odd = convolve_12tap(s, coeffs_x);
                    res_odd           = vshlq_s32(vaddq_s32(res_odd, round_const_x), -round_shift_x);

                    int16x8_t res_even1 = vqmovn_high_s32(vqmovn_s32(res_even), res_even);
                    int16x8_t res_odd1  = vqmovn_high_s32(vqmovn_s32(res_odd), res_odd);
                    int16x8_t res       = vzip1q_s16(res_even1, res_odd1);

                    vst1q_s16(&im_block[i * im_stride], res);
                }
            }

            /* Vertical filter */
            {
                int16x8_t s0  = vld1q_s16(im_block + 0 * im_stride);
                int16x8_t s1  = vld1q_s16(im_block + 1 * im_stride);
                int16x8_t s2  = vld1q_s16(im_block + 2 * im_stride);
                int16x8_t s3  = vld1q_s16(im_block + 3 * im_stride);
                int16x8_t s4  = vld1q_s16(im_block + 4 * im_stride);
                int16x8_t s5  = vld1q_s16(im_block + 5 * im_stride);
                int16x8_t s6  = vld1q_s16(im_block + 6 * im_stride);
                int16x8_t s7  = vld1q_s16(im_block + 7 * im_stride);
                int16x8_t s8  = vld1q_s16(im_block + 8 * im_stride);
                int16x8_t s9  = vld1q_s16(im_block + 9 * im_stride);
                int16x8_t s10 = vld1q_s16(im_block + 10 * im_stride);

                s[0] = vzip1q_s16(s0, s1);
                s[1] = vzip1q_s16(s2, s3);
                s[2] = vzip1q_s16(s4, s5);
                s[3] = vzip1q_s16(s6, s7);
                s[4] = vzip1q_s16(s8, s9);

                s[6]  = vzip2q_s16(s0, s1);
                s[7]  = vzip2q_s16(s2, s3);
                s[8]  = vzip2q_s16(s4, s5);
                s[9]  = vzip2q_s16(s6, s7);
                s[10] = vzip2q_s16(s8, s9);

                s[12] = vzip1q_s16(s1, s2);
                s[13] = vzip1q_s16(s3, s4);
                s[14] = vzip1q_s16(s5, s6);
                s[15] = vzip1q_s16(s7, s8);
                s[16] = vzip1q_s16(s9, s10);

                s[18] = vzip2q_s16(s1, s2);
                s[19] = vzip2q_s16(s3, s4);
                s[20] = vzip2q_s16(s5, s6);
                s[21] = vzip2q_s16(s7, s8);
                s[22] = vzip2q_s16(s9, s10);

                for (i = 0; i < h; i += 2) {
                    const int16_t *data = &im_block[i * im_stride];

                    int16x8_t s11 = vld1q_s16(data + 11 * im_stride);
                    int16x8_t s12 = vld1q_s16(data + 12 * im_stride);

                    s[5]  = vzip1q_s16(s10, s11);
                    s[11] = vzip2q_s16(s10, s11);

                    s[17] = vzip1q_s16(s11, s12);
                    s[23] = vzip2q_s16(s11, s12);

                    const int32x4_t res_a0       = convolve_12tap(s, coeffs_y);
                    int32x4_t       res_a_round0 = vshlq_s32(vaddq_s32(res_a0, round_const_y), -round_shift_y);
                    res_a_round0 = vshlq_s32(vaddq_s32(res_a_round0, round_const_bits), -round_shift_bits);

                    const int32x4_t res_a1       = convolve_12tap(s + 12, coeffs_y);
                    int32x4_t       res_a_round1 = vshlq_s32(vaddq_s32(res_a1, round_const_y), -round_shift_y);
                    res_a_round1 = vshlq_s32(vaddq_s32(res_a_round1, round_const_bits), -round_shift_bits);

                    if (w - j > 4) {
                        const int32x4_t res_b0       = convolve_12tap(s + 6, coeffs_y);
                        int32x4_t       res_b_round0 = vshlq_s32(vaddq_s32(res_b0, round_const_y), -round_shift_y);
                        res_b_round0 = vshlq_s32(vaddq_s32(res_b_round0, round_const_bits), -round_shift_bits);

                        const int32x4_t res_b1       = convolve_12tap(s + 18, coeffs_y);
                        int32x4_t       res_b_round1 = vshlq_s32(vaddq_s32(res_b1, round_const_y), -round_shift_y);
                        res_b_round1 = vshlq_s32(vaddq_s32(res_b_round1, round_const_bits), -round_shift_bits);

                        int16x8_t res_16bit0 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_b_round0);
                        res_16bit0           = vminq_s16(res_16bit0, clip_pixel);
                        res_16bit0           = vmaxq_s16(res_16bit0, zero);

                        int16x8_t res_16bit1 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_b_round1);
                        res_16bit1           = vminq_s16(res_16bit1, clip_pixel);
                        res_16bit1           = vmaxq_s16(res_16bit1, zero);

                        vst1q_s16((int16_t *)&dst[i * dst_stride + j], res_16bit0);
                        vst1q_s16((int16_t *)&dst[i * dst_stride + j + dst_stride], res_16bit1);
                    } else if (w == 4) {
                        int16x8_t res_a_round0_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_a_round0);
                        res_a_round0_s16           = vminq_s16(res_a_round0_s16, clip_pixel);
                        res_a_round0_s16           = vmaxq_s16(res_a_round0_s16, zero);

                        int16x8_t res_a_round1_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_a_round1);
                        res_a_round1_s16           = vminq_s16(res_a_round1_s16, clip_pixel);
                        res_a_round1_s16           = vmaxq_s16(res_a_round1_s16, zero);

                        vst1_s64((int64_t *)&dst[i * dst_stride + j],
                                 vreinterpret_s64_s16(vget_low_s16(res_a_round0_s16)));
                        vst1_s64((int64_t *)&dst[i * dst_stride + j + dst_stride],
                                 vreinterpret_s64_s16(vget_low_s16(res_a_round1_s16)));
                    } else {
                        int16x8_t res_a_round0_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_a_round0);
                        res_a_round0_s16           = vminq_s16(res_a_round0_s16, clip_pixel);
                        res_a_round0_s16           = vmaxq_s16(res_a_round0_s16, zero);

                        int16x8_t res_a_round1_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_a_round1);
                        res_a_round1_s16           = vminq_s16(res_a_round1_s16, clip_pixel);
                        res_a_round1_s16           = vmaxq_s16(res_a_round1_s16, zero);

                        *((uint32_t *)(&dst[i * dst_stride + j])) = vgetq_lane_s32(
                            vreinterpretq_s32_s16(res_a_round0_s16), 0);

                        *((uint32_t *)(&dst[i * dst_stride + j + dst_stride])) = vgetq_lane_s32(
                            vreinterpretq_s32_s16(res_a_round1_s16), 0);
                    }
                    s[0] = s[1];
                    s[1] = s[2];
                    s[2] = s[3];
                    s[3] = s[4];
                    s[4] = s[5];

                    s[6]  = s[7];
                    s[7]  = s[8];
                    s[8]  = s[9];
                    s[9]  = s[10];
                    s[10] = s[11];

                    s[12] = s[13];
                    s[13] = s[14];
                    s[14] = s[15];
                    s[15] = s[16];
                    s[16] = s[17];

                    s[18] = s[19];
                    s[19] = s[20];
                    s[20] = s[21];
                    s[21] = s[22];
                    s[22] = s[23];

                    s10 = s12;
                }
            }
        }
    } else {
        int16x8_t coeffs_x[4], coeffs_y[4], s[16];
        prepare_coeffs(x_filter, coeffs_x);
        prepare_coeffs(y_filter, coeffs_y);

        for (j = 0; j < w; j += 8) {
            /* Horizontal filter */
            {
                for (i = 0; i < im_h; i += 1) {
                    const int16x8_t row00 = vld1q_s16((const int16_t *)&src_ptr[i * src_stride + j]);
                    const int16x8_t row01 = vld1q_s16((const int16_t *)&src_ptr[i * src_stride + (j + 8)]);

                    // even pixels
                    s[0] = vextq_s16(row00, row01, 0);
                    s[1] = vextq_s16(row00, row01, 2);
                    s[2] = vextq_s16(row00, row01, 4);
                    s[3] = vextq_s16(row00, row01, 6);

                    int32x4_t res_even = svt_aom_convolve(s, coeffs_x);
                    res_even           = vshlq_s32(vaddq_s32(res_even, round_const_x), vdupq_n_s32(-round_shift_x[0]));

                    // odd pixels
                    s[0] = vextq_s16(row00, row01, 1);
                    s[1] = vextq_s16(row00, row01, 3);
                    s[2] = vextq_s16(row00, row01, 5);
                    s[3] = vextq_s16(row00, row01, 7);

                    int32x4_t res_odd = svt_aom_convolve(s, coeffs_x);
                    res_odd           = vshlq_s32(vaddq_s32(res_odd, round_const_x), vdupq_n_s32(-round_shift_x[0]));

                    int16x8_t res_even1 = vqmovn_high_s32(vqmovn_s32(res_even), res_even);
                    int16x8_t res_odd1  = vqmovn_high_s32(vqmovn_s32(res_odd), res_odd);
                    int16x8_t res       = vzip1q_s16(res_even1, res_odd1);

                    vst1q_s16(&im_block[i * im_stride], res);
                }
            }

            /* Vertical filter */
            {
                int16x8_t s0 = vld1q_s16(im_block + 0 * im_stride);
                int16x8_t s1 = vld1q_s16(im_block + 1 * im_stride);
                int16x8_t s2 = vld1q_s16(im_block + 2 * im_stride);
                int16x8_t s3 = vld1q_s16(im_block + 3 * im_stride);
                int16x8_t s4 = vld1q_s16(im_block + 4 * im_stride);
                int16x8_t s5 = vld1q_s16(im_block + 5 * im_stride);
                int16x8_t s6 = vld1q_s16(im_block + 6 * im_stride);

                s[0] = vzip1q_s16(s0, s1);
                s[1] = vzip1q_s16(s2, s3);
                s[2] = vzip1q_s16(s4, s5);

                s[4] = vzip2q_s16(s0, s1);
                s[5] = vzip2q_s16(s2, s3);
                s[6] = vzip2q_s16(s4, s5);

                s[0 + 8] = vzip1q_s16(s1, s2);
                s[1 + 8] = vzip1q_s16(s3, s4);
                s[2 + 8] = vzip1q_s16(s5, s6);

                s[4 + 8] = vzip2q_s16(s1, s2);
                s[5 + 8] = vzip2q_s16(s3, s4);
                s[6 + 8] = vzip2q_s16(s5, s6);

                for (i = 0; i < h; i += 2) {
                    const int16_t *data = &im_block[i * im_stride];

                    int16x8_t s7 = vld1q_s16(data + 7 * im_stride);
                    int16x8_t s8 = vld1q_s16(data + 8 * im_stride);

                    s[3] = vzip1q_s16(s6, s7);
                    s[7] = vzip2q_s16(s6, s7);

                    s[3 + 8] = vzip1q_s16(s7, s8);
                    s[7 + 8] = vzip2q_s16(s7, s8);

                    const int32x4_t res_a0       = svt_aom_convolve(s, coeffs_y);
                    int32x4_t       res_a_round0 = vshlq_s32(vaddq_s32(res_a0, round_const_y),
                                                       vdupq_n_s32(-round_shift_y[0]));
                    res_a_round0                 = vshlq_s32(vaddq_s32(res_a_round0, round_const_bits),
                                             vdupq_n_s32(-round_shift_bits[0]));

                    const int32x4_t res_a1       = svt_aom_convolve(s + 8, coeffs_y);
                    int32x4_t       res_a_round1 = vshlq_s32(vaddq_s32(res_a1, round_const_y),
                                                       vdupq_n_s32(-round_shift_y[0]));
                    res_a_round1                 = vshlq_s32(vaddq_s32(res_a_round1, round_const_bits),
                                             vdupq_n_s32(-round_shift_bits[0]));

                    if (w - j > 4) {
                        const int32x4_t res_b0       = svt_aom_convolve(s + 4, coeffs_y);
                        int32x4_t       res_b_round0 = vshlq_s32(vaddq_s32(res_b0, round_const_y),
                                                           vdupq_n_s32(-round_shift_y[0]));
                        res_b_round0                 = vshlq_s32(vaddq_s32(res_b_round0, round_const_bits),
                                                 vdupq_n_s32(-round_shift_bits[0]));

                        const int32x4_t res_b1       = svt_aom_convolve(s + 4 + 8, coeffs_y);
                        int32x4_t       res_b_round1 = vshlq_s32(vaddq_s32(res_b1, round_const_y),
                                                           vdupq_n_s32(-round_shift_y[0]));
                        res_b_round1                 = vshlq_s32(vaddq_s32(res_b_round1, round_const_bits),
                                                 vdupq_n_s32(-round_shift_bits[0]));

                        int16x8_t res_16bit0 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_b_round0);
                        res_16bit0           = vminq_s16(res_16bit0, clip_pixel);
                        res_16bit0           = vmaxq_s16(res_16bit0, zero);

                        int16x8_t res_16bit1 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_b_round1);
                        res_16bit1           = vminq_s16(res_16bit1, clip_pixel);
                        res_16bit1           = vmaxq_s16(res_16bit1, zero);

                        vst1q_s16((int16_t *)&dst[i * dst_stride + j], res_16bit0);
                        vst1q_s16((int16_t *)&dst[i * dst_stride + j + dst_stride], res_16bit1);
                    } else if (w == 4) {
                        int16x8_t res_a_round0_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_a_round0);
                        res_a_round0_s16           = vminq_s16(res_a_round0_s16, clip_pixel);
                        res_a_round0_s16           = vmaxq_s16(res_a_round0_s16, zero);

                        int16x8_t res_a_round1_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_a_round1);
                        res_a_round1_s16           = vminq_s16(res_a_round1_s16, clip_pixel);
                        res_a_round1_s16           = vmaxq_s16(res_a_round1_s16, zero);

                        vst1_s64((int64_t *)&dst[i * dst_stride + j],
                                 vreinterpret_s64_s16(vget_low_s16(res_a_round0_s16)));
                        vst1_s64((int64_t *)&dst[i * dst_stride + j + dst_stride],
                                 vreinterpret_s64_s16(vget_low_s16(res_a_round1_s16)));
                    } else {
                        int16x8_t res_a_round0_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_a_round0);
                        res_a_round0_s16           = vminq_s16(res_a_round0_s16, clip_pixel);
                        res_a_round0_s16           = vmaxq_s16(res_a_round0_s16, zero);

                        int16x8_t res_a_round1_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_a_round1);
                        res_a_round1_s16           = vminq_s16(res_a_round1_s16, clip_pixel);
                        res_a_round1_s16           = vmaxq_s16(res_a_round1_s16, zero);

                        *((uint32_t *)(&dst[i * dst_stride + j])) = vgetq_lane_s32(
                            vreinterpretq_s32_s16(res_a_round0_s16), 0);

                        *((uint32_t *)(&dst[i * dst_stride + j + dst_stride])) = vgetq_lane_s32(
                            vreinterpretq_s32_s16(res_a_round1_s16), 0);
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
}
