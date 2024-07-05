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

#include "common_dsp_rtcd.h"
#include "convolve.h"
#include "definitions.h"
#include "mem_neon.h"

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

static INLINE int32x4_t convolve_12tap(const uint16x8_t *s, const int16x8_t *coeffs) {
    const int32x4_t d0     = vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[0])), vget_low_s16(coeffs[0]));
    const int32x4_t d1     = vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[1])), vget_low_s16(coeffs[1]));
    const int32x4_t d2     = vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[2])), vget_low_s16(coeffs[2]));
    const int32x4_t d3     = vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[3])), vget_low_s16(coeffs[3]));
    const int32x4_t d4     = vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[4])), vget_low_s16(coeffs[4]));
    const int32x4_t d5     = vmull_s16(vget_low_s16(vreinterpretq_s16_u16(s[5])), vget_low_s16(coeffs[5]));
    const int32x4_t d_0123 = vaddq_s32(vaddq_s32(d0, d1), vaddq_s32(d2, d3));
    const int32x4_t d      = vaddq_s32(vaddq_s32(d4, d5), d_0123);
    return d;
}

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

void svt_av1_highbd_convolve_y_sr_neon(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
                                       int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
                                       const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
                                       const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd) {
    (void)filter_params_x;
    (void)subpel_x_q4;
    (void)conv_params;
    int                   i, j;
    const int             fo_vert = filter_params_y->taps / 2 - 1;
    const uint16_t *const src_ptr = src - fo_vert * src_stride;
    const int             bits    = FILTER_BITS;

    const int32x4_t round_shift_bits = vsetq_lane_s32(bits, vdupq_n_s32(0), 0);
    const int32x4_t round_const_bits = vdupq_n_s32((1 << bits) >> 1);
    const int16x8_t clip_pixel       = vdupq_n_s16(bd == 10 ? 1023 : (bd == 12 ? 4095 : 255));
    const int16x8_t zero             = vdupq_n_s16(0);

    if (filter_params_y->taps == 12) {
        uint16x8_t           s[24];
        int16x8_t            coeffs_y[6];
        const int16_t *const y_filter = av1_get_interp_filter_subpel_kernel(*filter_params_y,
                                                                            subpel_y_q4 & SUBPEL_MASK);
        svt_prepare_coeffs_12tap(y_filter, coeffs_y);

        for (j = 0; j < w; j += 8) {
            const uint16_t *data = &src_ptr[j];
            /* Vertical filter */
            uint16x8_t s0  = vld1q_u16(data + 0 * src_stride);
            uint16x8_t s1  = vld1q_u16(data + 1 * src_stride);
            uint16x8_t s2  = vld1q_u16(data + 2 * src_stride);
            uint16x8_t s3  = vld1q_u16(data + 3 * src_stride);
            uint16x8_t s4  = vld1q_u16(data + 4 * src_stride);
            uint16x8_t s5  = vld1q_u16(data + 5 * src_stride);
            uint16x8_t s6  = vld1q_u16(data + 6 * src_stride);
            uint16x8_t s7  = vld1q_u16(data + 7 * src_stride);
            uint16x8_t s8  = vld1q_u16(data + 8 * src_stride);
            uint16x8_t s9  = vld1q_u16(data + 9 * src_stride);
            uint16x8_t s10 = vld1q_u16(data + 10 * src_stride);

            s[0] = vzip1q_u16(s0, s1);
            s[1] = vzip1q_u16(s2, s3);
            s[2] = vzip1q_u16(s4, s5);
            s[3] = vzip1q_u16(s6, s7);
            s[4] = vzip1q_u16(s8, s9);

            s[6]  = vzip2q_u16(s0, s1);
            s[7]  = vzip2q_u16(s2, s3);
            s[8]  = vzip2q_u16(s4, s5);
            s[9]  = vzip2q_u16(s6, s7);
            s[10] = vzip2q_u16(s8, s9);

            s[12] = vzip1q_u16(s1, s2);
            s[13] = vzip1q_u16(s3, s4);
            s[14] = vzip1q_u16(s5, s6);
            s[15] = vzip1q_u16(s7, s8);
            s[16] = vzip1q_u16(s9, s10);

            s[18] = vzip2q_u16(s1, s2);
            s[19] = vzip2q_u16(s3, s4);
            s[20] = vzip2q_u16(s5, s6);
            s[21] = vzip2q_u16(s7, s8);
            s[22] = vzip2q_u16(s9, s10);

            for (i = 0; i < h; i += 2) {
                data = &src_ptr[i * src_stride + j];

                uint16x8_t s11 = vld1q_u16(data + 11 * src_stride);
                uint16x8_t s12 = vld1q_u16(data + 12 * src_stride);

                s[5]  = vzip1q_u16(s10, s11);
                s[11] = vzip2q_u16(s10, s11);

                s[17] = vzip1q_u16(s11, s12);
                s[23] = vzip2q_u16(s11, s12);

                const int32x4_t res_a0       = convolve_12tap(s, coeffs_y);
                int32x4_t       res_a_round0 = vshlq_s32(vaddq_s32(res_a0, round_const_bits),
                                                   vdupq_n_s32(-round_shift_bits[0]));

                const int32x4_t res_a1       = convolve_12tap(s + 12, coeffs_y);
                int32x4_t       res_a_round1 = vshlq_s32(vaddq_s32(res_a1, round_const_bits),
                                                   vdupq_n_s32(-round_shift_bits[0]));

                if (w - j > 4) {
                    const int32x4_t res_b0       = convolve_12tap(s + 6, coeffs_y);
                    int32x4_t       res_b_round0 = vshlq_s32(vaddq_s32(res_b0, round_const_bits),
                                                       vdupq_n_s32(-round_shift_bits[0]));

                    const int32x4_t res_b1       = convolve_12tap(s + 18, coeffs_y);
                    int32x4_t       res_b_round1 = vshlq_s32(vaddq_s32(res_b1, round_const_bits),
                                                       vdupq_n_s32(-round_shift_bits[0]));

                    int16x8_t res_16bit0 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_b_round0);
                    res_16bit0           = vminq_s16(res_16bit0, clip_pixel);
                    res_16bit0           = vmaxq_s16(res_16bit0, zero);

                    int16x8_t res_16bit1 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_b_round1);
                    res_16bit1           = vminq_s16(res_16bit1, clip_pixel);
                    res_16bit1           = vmaxq_s16(res_16bit1, zero);

                    vst1q_u16(&dst[i * dst_stride + j], vreinterpretq_u16_s16(res_16bit0));
                    vst1q_u16(&dst[i * dst_stride + j + dst_stride], vreinterpretq_u16_s16(res_16bit1));
                } else if (w == 4) {
                    int16x8_t res_a_round0_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_a_round0);
                    res_a_round0_s16           = vminq_s16(res_a_round0_s16, clip_pixel);
                    res_a_round0_s16           = vmaxq_s16(res_a_round0_s16, zero);

                    int16x8_t res_a_round1_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_a_round1);
                    res_a_round1_s16           = vminq_s16(res_a_round1_s16, clip_pixel);
                    res_a_round1_s16           = vmaxq_s16(res_a_round1_s16, zero);

                    vst1_u16(&dst[i * dst_stride + j], vget_low_u16(vreinterpretq_u16_s16(res_a_round0_s16)));
                    vst1_u16(&dst[i * dst_stride + j + dst_stride],
                             vget_low_u16(vreinterpretq_u16_s16(res_a_round1_s16)));
                } else {
                    int16x8_t res_a_round0_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_a_round0);
                    res_a_round0_s16           = vminq_s16(res_a_round0_s16, clip_pixel);
                    res_a_round0_s16           = vmaxq_s16(res_a_round0_s16, zero);

                    int16x8_t res_a_round1_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_a_round1);
                    res_a_round1_s16           = vminq_s16(res_a_round1_s16, clip_pixel);
                    res_a_round1_s16           = vmaxq_s16(res_a_round1_s16, zero);

                    *((uint32_t *)(&dst[i * dst_stride + j])) = vgetq_lane_u32(vreinterpretq_u32_s16(res_a_round0_s16),
                                                                               0);
                    *((uint32_t *)(&dst[i * dst_stride + j + dst_stride])) = vgetq_lane_u32(
                        vreinterpretq_u32_s16(res_a_round1_s16), 0);
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
    } else {
        uint16x8_t     s[16];
        int16x8_t      coeffs_y[4];
        const int16_t *filter = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_q4 & SUBPEL_MASK);
        prepare_coeffs(filter, coeffs_y);

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

                    uint16x8_t s7 = vld1q_u16(data + 7 * src_stride);
                    uint16x8_t s8 = vld1q_u16(data + 8 * src_stride);

                    s[3] = vzip1q_u16(s6, s7);
                    s[7] = vzip2q_u16(s6, s7);

                    s[3 + 8] = vzip1q_u16(s7, s8);
                    s[7 + 8] = vzip2q_u16(s7, s8);

                    const int32x4_t res_a0       = svt_aom_convolve(s, coeffs_y);
                    int32x4_t       res_a_round0 = vshlq_s32(vaddq_s32(res_a0, round_const_bits),
                                                       vdupq_n_s32(-round_shift_bits[0]));

                    const int32x4_t res_a1       = svt_aom_convolve(s + 8, coeffs_y);
                    int32x4_t       res_a_round1 = vshlq_s32(vaddq_s32(res_a1, round_const_bits),
                                                       vdupq_n_s32(-round_shift_bits[0]));

                    if (w - j > 4) {
                        const int32x4_t res_b0       = svt_aom_convolve(s + 4, coeffs_y);
                        int32x4_t       res_b_round0 = vshlq_s32(vaddq_s32(res_b0, round_const_bits),
                                                           vdupq_n_s32(-round_shift_bits[0]));

                        const int32x4_t res_b1       = svt_aom_convolve(s + 4 + 8, coeffs_y);
                        int32x4_t       res_b_round1 = vshlq_s32(vaddq_s32(res_b1, round_const_bits),
                                                           vdupq_n_s32(-round_shift_bits[0]));
                        int16x8_t       res_16bit0   = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_b_round0);
                        res_16bit0                   = vminq_s16(res_16bit0, clip_pixel);
                        res_16bit0                   = vmaxq_s16(res_16bit0, zero);

                        int16x8_t res_16bit1 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_b_round1);
                        res_16bit1           = vminq_s16(res_16bit1, clip_pixel);
                        res_16bit1           = vmaxq_s16(res_16bit1, zero);

                        vst1q_u16(&dst[i * dst_stride + j], vreinterpretq_u16_s16(res_16bit0));
                        vst1q_u16(&dst[i * dst_stride + j + dst_stride], vreinterpretq_u16_s16(res_16bit1));
                    } else if (w == 4) {
                        int16x8_t res_a_round0_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_a_round0);
                        res_a_round0_s16           = vminq_s16(res_a_round0_s16, clip_pixel);
                        res_a_round0_s16           = vmaxq_s16(res_a_round0_s16, zero);

                        int16x8_t res_a_round1_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_a_round1);
                        res_a_round1_s16           = vminq_s16(res_a_round1_s16, clip_pixel);
                        res_a_round1_s16           = vmaxq_s16(res_a_round1_s16, zero);

                        vst1_u16(&dst[i * dst_stride + j], vget_low_u16(vreinterpretq_u16_s16(res_a_round0_s16)));
                        vst1_u16(&dst[i * dst_stride + j + dst_stride],
                                 vget_low_u16(vreinterpretq_u16_s16(res_a_round1_s16)));
                    } else {
                        int16x8_t res_a_round0_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round0), res_a_round0);
                        res_a_round0_s16           = vminq_s16(res_a_round0_s16, clip_pixel);
                        res_a_round0_s16           = vmaxq_s16(res_a_round0_s16, zero);

                        int16x8_t res_a_round1_s16 = vqmovn_high_s32(vqmovn_s32(res_a_round1), res_a_round1);
                        res_a_round1_s16           = vminq_s16(res_a_round1_s16, clip_pixel);
                        res_a_round1_s16           = vmaxq_s16(res_a_round1_s16, zero);

                        *((uint32_t *)(&dst[i * dst_stride + j])) = vgetq_lane_u32(
                            vreinterpretq_u32_s16(res_a_round0_s16), 0);
                        *((uint32_t *)(&dst[i * dst_stride + j + dst_stride])) = vgetq_lane_u32(
                            vreinterpretq_u32_s16(res_a_round1_s16), 0);
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

static INLINE uint16x4_t highbd_convolve6_4_2d_v(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2,
                                                 const int16x4_t s3, const int16x4_t s4, const int16x4_t s5,
                                                 const int16x8_t y_filter, const int32x4_t round_shift,
                                                 const int32x4_t offset, const uint16x4_t max) {
    // Values at indices 0 and 7 of y_filter are zero.
    const int16x4_t y_filter_0_3 = vget_low_s16(y_filter);
    const int16x4_t y_filter_4_7 = vget_high_s16(y_filter);

    int32x4_t sum = vmlal_lane_s16(offset, s0, y_filter_0_3, 1);
    sum           = vmlal_lane_s16(sum, s1, y_filter_0_3, 2);
    sum           = vmlal_lane_s16(sum, s2, y_filter_0_3, 3);
    sum           = vmlal_lane_s16(sum, s3, y_filter_4_7, 0);
    sum           = vmlal_lane_s16(sum, s4, y_filter_4_7, 1);
    sum           = vmlal_lane_s16(sum, s5, y_filter_4_7, 2);

    sum            = vshlq_s32(sum, round_shift);
    uint16x4_t res = vqmovun_s32(sum);
    return vmin_u16(res, max);
}

static INLINE uint16x8_t highbd_convolve6_8_2d_v(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
                                                 const int16x8_t s3, const int16x8_t s4, const int16x8_t s5,
                                                 const int16x8_t y_filter, const int32x4_t round_shift,
                                                 const int32x4_t offset, const uint16x8_t max) {
    // Values at indices 0 and 7 of y_filter are zero.
    const int16x4_t y_filter_0_3 = vget_low_s16(y_filter);
    const int16x4_t y_filter_4_7 = vget_high_s16(y_filter);

    int32x4_t sum0 = vmlal_lane_s16(offset, vget_low_s16(s0), y_filter_0_3, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s1), y_filter_0_3, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s2), y_filter_0_3, 3);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s3), y_filter_4_7, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s4), y_filter_4_7, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s5), y_filter_4_7, 2);

    int32x4_t sum1 = vmlal_lane_s16(offset, vget_high_s16(s0), y_filter_0_3, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s1), y_filter_0_3, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s2), y_filter_0_3, 3);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s3), y_filter_4_7, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s4), y_filter_4_7, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s5), y_filter_4_7, 2);

    sum0 = vshlq_s32(sum0, round_shift);
    sum1 = vshlq_s32(sum1, round_shift);

    uint16x8_t res = vcombine_u16(vqmovun_s32(sum0), vqmovun_s32(sum1));
    return vminq_u16(res, max);
}

static INLINE void highbd_convolve_2d_sr_vert_6tap_neon(const uint16_t *src_ptr, int src_stride, uint16_t *dst_ptr,
                                                        int dst_stride, int w, int h, const int16_t *y_filter_ptr,
                                                        ConvolveParams *conv_params, int bd, const int offset) {
    const int16x8_t y_filter         = vld1q_s16(y_filter_ptr);
    const int32x4_t offset_s32       = vdupq_n_s32(offset);
    const int       round1_shift     = conv_params->round_1;
    const int32x4_t round1_shift_s32 = vdupq_n_s32(-round1_shift);

    if (w == 4) {
        const uint16x4_t max = vdup_n_u16((1 << bd) - 1);
        const int16_t   *s   = (const int16_t *)src_ptr;
        uint16_t        *d   = dst_ptr;
        int16x4_t        s0, s1, s2, s3, s4;
        load_s16_4x5(s, src_stride, &s0, &s1, &s2, &s3, &s4);
        s += 5 * src_stride;

        do {
            int16x4_t s5, s6, s7, s8;
            load_s16_4x4(s, src_stride, &s5, &s6, &s7, &s8);

            uint16x4_t d0 = highbd_convolve6_4_2d_v(
                s0, s1, s2, s3, s4, s5, y_filter, round1_shift_s32, offset_s32, max);
            uint16x4_t d1 = highbd_convolve6_4_2d_v(
                s1, s2, s3, s4, s5, s6, y_filter, round1_shift_s32, offset_s32, max);
            uint16x4_t d2 = highbd_convolve6_4_2d_v(
                s2, s3, s4, s5, s6, s7, y_filter, round1_shift_s32, offset_s32, max);
            uint16x4_t d3 = highbd_convolve6_4_2d_v(
                s3, s4, s5, s6, s7, s8, y_filter, round1_shift_s32, offset_s32, max);

            store_u16_4x4(d, dst_stride, d0, d1, d2, d3);

            s0 = s4;
            s1 = s5;
            s2 = s6;
            s3 = s7;
            s4 = s8;
            s += 4 * src_stride;
            d += 4 * dst_stride;
            h -= 4;
        } while (h != 0);
    } else {
        const uint16x8_t max = vdupq_n_u16((1 << bd) - 1);

        do {
            int            height = h;
            const int16_t *s      = (const int16_t *)src_ptr;
            uint16_t      *d      = dst_ptr;
            int16x8_t      s0, s1, s2, s3, s4;
            load_s16_8x5(s, src_stride, &s0, &s1, &s2, &s3, &s4);
            s += 5 * src_stride;

            do {
                int16x8_t s5, s6, s7, s8;
                load_s16_8x4(s, src_stride, &s5, &s6, &s7, &s8);

                uint16x8_t d0 = highbd_convolve6_8_2d_v(
                    s0, s1, s2, s3, s4, s5, y_filter, round1_shift_s32, offset_s32, max);
                uint16x8_t d1 = highbd_convolve6_8_2d_v(
                    s1, s2, s3, s4, s5, s6, y_filter, round1_shift_s32, offset_s32, max);
                uint16x8_t d2 = highbd_convolve6_8_2d_v(
                    s2, s3, s4, s5, s6, s7, y_filter, round1_shift_s32, offset_s32, max);
                uint16x8_t d3 = highbd_convolve6_8_2d_v(
                    s3, s4, s5, s6, s7, s8, y_filter, round1_shift_s32, offset_s32, max);

                store_u16_8x4(d, dst_stride, d0, d1, d2, d3);

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

static INLINE uint16x4_t highbd_convolve8_4_2d_v(const int16x4_t s0, const int16x4_t s1, const int16x4_t s2,
                                                 const int16x4_t s3, const int16x4_t s4, const int16x4_t s5,
                                                 const int16x4_t s6, const int16x4_t s7, const int16x8_t y_filter,
                                                 const int32x4_t round_shift, const int32x4_t offset,
                                                 const uint16x4_t max) {
    const int16x4_t y_filter_lo = vget_low_s16(y_filter);
    const int16x4_t y_filter_hi = vget_high_s16(y_filter);

    int32x4_t sum = vmlal_lane_s16(offset, s0, y_filter_lo, 0);
    sum           = vmlal_lane_s16(sum, s1, y_filter_lo, 1);
    sum           = vmlal_lane_s16(sum, s2, y_filter_lo, 2);
    sum           = vmlal_lane_s16(sum, s3, y_filter_lo, 3);
    sum           = vmlal_lane_s16(sum, s4, y_filter_hi, 0);
    sum           = vmlal_lane_s16(sum, s5, y_filter_hi, 1);
    sum           = vmlal_lane_s16(sum, s6, y_filter_hi, 2);
    sum           = vmlal_lane_s16(sum, s7, y_filter_hi, 3);

    sum            = vshlq_s32(sum, round_shift);
    uint16x4_t res = vqmovun_s32(sum);
    return vmin_u16(res, max);
}

static INLINE uint16x8_t highbd_convolve8_8_2d_v(const int16x8_t s0, const int16x8_t s1, const int16x8_t s2,
                                                 const int16x8_t s3, const int16x8_t s4, const int16x8_t s5,
                                                 const int16x8_t s6, const int16x8_t s7, const int16x8_t y_filter,
                                                 const int32x4_t round_shift, const int32x4_t offset,
                                                 const uint16x8_t max) {
    const int16x4_t y_filter_lo = vget_low_s16(y_filter);
    const int16x4_t y_filter_hi = vget_high_s16(y_filter);

    int32x4_t sum0 = vmlal_lane_s16(offset, vget_low_s16(s0), y_filter_lo, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s1), y_filter_lo, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s2), y_filter_lo, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s3), y_filter_lo, 3);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s4), y_filter_hi, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s5), y_filter_hi, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s6), y_filter_hi, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s7), y_filter_hi, 3);

    int32x4_t sum1 = vmlal_lane_s16(offset, vget_high_s16(s0), y_filter_lo, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s1), y_filter_lo, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s2), y_filter_lo, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s3), y_filter_lo, 3);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s4), y_filter_hi, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s5), y_filter_hi, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s6), y_filter_hi, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s7), y_filter_hi, 3);

    sum0 = vshlq_s32(sum0, round_shift);
    sum1 = vshlq_s32(sum1, round_shift);

    uint16x8_t res = vcombine_u16(vqmovun_s32(sum0), vqmovun_s32(sum1));
    return vminq_u16(res, max);
}

static INLINE void highbd_convolve_2d_sr_vert_8tap_neon(const uint16_t *src_ptr, int src_stride, uint16_t *dst_ptr,
                                                        int dst_stride, int w, int h, const int16_t *y_filter_ptr,
                                                        ConvolveParams *conv_params, int bd, const int offset) {
    const int16x8_t y_filter         = vld1q_s16(y_filter_ptr);
    const int32x4_t offset_s32       = vdupq_n_s32(offset);
    const int       round1_shift     = conv_params->round_1;
    const int32x4_t round1_shift_s32 = vdupq_n_s32(-round1_shift);

    if (w == 4) {
        const uint16x4_t max = vdup_n_u16((1 << bd) - 1);
        const int16_t   *s   = (const int16_t *)src_ptr;
        uint16_t        *d   = dst_ptr;

        int16x4_t s0, s1, s2, s3, s4, s5, s6;
        load_s16_4x7(s, src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6);
        s += 7 * src_stride;

        do {
            int16x4_t s7, s8, s9, s10;
            load_s16_4x4(s, src_stride, &s7, &s8, &s9, &s10);

            uint16x4_t d0 = highbd_convolve8_4_2d_v(
                s0, s1, s2, s3, s4, s5, s6, s7, y_filter, round1_shift_s32, offset_s32, max);
            uint16x4_t d1 = highbd_convolve8_4_2d_v(
                s1, s2, s3, s4, s5, s6, s7, s8, y_filter, round1_shift_s32, offset_s32, max);
            uint16x4_t d2 = highbd_convolve8_4_2d_v(
                s2, s3, s4, s5, s6, s7, s8, s9, y_filter, round1_shift_s32, offset_s32, max);
            uint16x4_t d3 = highbd_convolve8_4_2d_v(
                s3, s4, s5, s6, s7, s8, s9, s10, y_filter, round1_shift_s32, offset_s32, max);

            store_u16_4x4(d, dst_stride, d0, d1, d2, d3);

            s0 = s4;
            s1 = s5;
            s2 = s6;
            s3 = s7;
            s4 = s8;
            s5 = s9;
            s6 = s10;
            s += 4 * src_stride;
            d += 4 * dst_stride;
            h -= 4;
        } while (h != 0);
    } else {
        const uint16x8_t max = vdupq_n_u16((1 << bd) - 1);

        do {
            int            height = h;
            const int16_t *s      = (const int16_t *)src_ptr;
            uint16_t      *d      = dst_ptr;

            int16x8_t s0, s1, s2, s3, s4, s5, s6;
            load_s16_8x7(s, src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6);
            s += 7 * src_stride;

            do {
                int16x8_t s7, s8, s9, s10;
                load_s16_8x4(s, src_stride, &s7, &s8, &s9, &s10);

                uint16x8_t d0 = highbd_convolve8_8_2d_v(
                    s0, s1, s2, s3, s4, s5, s6, s7, y_filter, round1_shift_s32, offset_s32, max);
                uint16x8_t d1 = highbd_convolve8_8_2d_v(
                    s1, s2, s3, s4, s5, s6, s7, s8, y_filter, round1_shift_s32, offset_s32, max);
                uint16x8_t d2 = highbd_convolve8_8_2d_v(
                    s2, s3, s4, s5, s6, s7, s8, s9, y_filter, round1_shift_s32, offset_s32, max);
                uint16x8_t d3 = highbd_convolve8_8_2d_v(
                    s3, s4, s5, s6, s7, s8, s9, s10, y_filter, round1_shift_s32, offset_s32, max);

                store_u16_8x4(d, dst_stride, d0, d1, d2, d3);

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

static INLINE uint16x8_t highbd_convolve6_8_2d_h(const int16x8_t s[6], const int16x8_t x_filter,
                                                 const int32x4_t shift_s32, const int32x4_t offset) {
    // Values at indices 0 and 7 of y_filter are zero.
    const int16x4_t x_filter_0_3 = vget_low_s16(x_filter);
    const int16x4_t x_filter_4_7 = vget_high_s16(x_filter);

    int32x4_t sum0 = vmlal_lane_s16(offset, vget_low_s16(s[0]), x_filter_0_3, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[1]), x_filter_0_3, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[2]), x_filter_0_3, 3);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[3]), x_filter_4_7, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[4]), x_filter_4_7, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[5]), x_filter_4_7, 2);

    int32x4_t sum1 = vmlal_lane_s16(offset, vget_high_s16(s[0]), x_filter_0_3, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[1]), x_filter_0_3, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[2]), x_filter_0_3, 3);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[3]), x_filter_4_7, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[4]), x_filter_4_7, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[5]), x_filter_4_7, 2);

    sum0 = vqrshlq_s32(sum0, shift_s32);
    sum1 = vqrshlq_s32(sum1, shift_s32);

    return vcombine_u16(vqmovun_s32(sum0), vqmovun_s32(sum1));
}

static INLINE void highbd_convolve_2d_sr_horiz_6tap_neon(const uint16_t *src_ptr, int src_stride, uint16_t *dst_ptr,
                                                         int dst_stride, int w, int h, const int16_t *x_filter_ptr,
                                                         ConvolveParams *conv_params, const int offset) {
    // The smallest block height processed by the SIMD functions is 4, and the
    // horizontal convolution needs to process an extra (filter_taps/2 - 1) lines
    // for the vertical convolution.
    assert(h >= 5);
    const int32x4_t shift_s32  = vdupq_n_s32(-conv_params->round_0);
    const int32x4_t offset_s32 = vdupq_n_s32(offset);

    const int16x8_t x_filter = vld1q_s16(x_filter_ptr);
    int             height   = h;

    do {
        int            width = w;
        const int16_t *s     = (const int16_t *)src_ptr;
        uint16_t      *d     = dst_ptr;

        do {
            int16x8_t s0[6], s1[6], s2[6], s3[6];
            load_s16_8x6(s + 0 * src_stride, 1, &s0[0], &s0[1], &s0[2], &s0[3], &s0[4], &s0[5]);
            load_s16_8x6(s + 1 * src_stride, 1, &s1[0], &s1[1], &s1[2], &s1[3], &s1[4], &s1[5]);
            load_s16_8x6(s + 2 * src_stride, 1, &s2[0], &s2[1], &s2[2], &s2[3], &s2[4], &s2[5]);
            load_s16_8x6(s + 3 * src_stride, 1, &s3[0], &s3[1], &s3[2], &s3[3], &s3[4], &s3[5]);

            uint16x8_t d0 = highbd_convolve6_8_2d_h(s0, x_filter, shift_s32, offset_s32);
            uint16x8_t d1 = highbd_convolve6_8_2d_h(s1, x_filter, shift_s32, offset_s32);
            uint16x8_t d2 = highbd_convolve6_8_2d_h(s2, x_filter, shift_s32, offset_s32);
            uint16x8_t d3 = highbd_convolve6_8_2d_h(s3, x_filter, shift_s32, offset_s32);

            store_u16_8x4(d, dst_stride, d0, d1, d2, d3);

            s += 8;
            d += 8;
            width -= 8;
        } while (width != 0);
        src_ptr += 4 * src_stride;
        dst_ptr += 4 * dst_stride;
        height -= 4;
    } while (height > 4);
    do {
        int            width = w;
        const int16_t *s     = (const int16_t *)src_ptr;
        uint16_t      *d     = dst_ptr;

        do {
            int16x8_t s0[6];
            load_s16_8x6(s, 1, &s0[0], &s0[1], &s0[2], &s0[3], &s0[4], &s0[5]);

            uint16x8_t d0 = highbd_convolve6_8_2d_h(s0, x_filter, shift_s32, offset_s32);
            vst1q_u16(d, d0);

            s += 8;
            d += 8;
            width -= 8;
        } while (width != 0);
        src_ptr += src_stride;
        dst_ptr += dst_stride;
    } while (--height != 0);
}

static INLINE uint16x4_t highbd_convolve4_4_2d_h(const int16x4_t s[4], const int16x4_t x_filter,
                                                 const int32x4_t shift_s32, const int32x4_t offset) {
    int32x4_t sum = vmlal_lane_s16(offset, s[0], x_filter, 0);
    sum           = vmlal_lane_s16(sum, s[1], x_filter, 1);
    sum           = vmlal_lane_s16(sum, s[2], x_filter, 2);
    sum           = vmlal_lane_s16(sum, s[3], x_filter, 3);

    sum = vqrshlq_s32(sum, shift_s32);
    return vqmovun_s32(sum);
}

static INLINE uint16x8_t highbd_convolve8_8_2d_h(const int16x8_t s[8], const int16x8_t x_filter,
                                                 const int32x4_t shift_s32, const int32x4_t offset) {
    const int16x4_t x_filter_0_3 = vget_low_s16(x_filter);
    const int16x4_t x_filter_4_7 = vget_high_s16(x_filter);

    int32x4_t sum0 = vmlal_lane_s16(offset, vget_low_s16(s[0]), x_filter_0_3, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[1]), x_filter_0_3, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[2]), x_filter_0_3, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[3]), x_filter_0_3, 3);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[4]), x_filter_4_7, 0);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[5]), x_filter_4_7, 1);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[6]), x_filter_4_7, 2);
    sum0           = vmlal_lane_s16(sum0, vget_low_s16(s[7]), x_filter_4_7, 3);

    int32x4_t sum1 = vmlal_lane_s16(offset, vget_high_s16(s[0]), x_filter_0_3, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[1]), x_filter_0_3, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[2]), x_filter_0_3, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[3]), x_filter_0_3, 3);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[4]), x_filter_4_7, 0);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[5]), x_filter_4_7, 1);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[6]), x_filter_4_7, 2);
    sum1           = vmlal_lane_s16(sum1, vget_high_s16(s[7]), x_filter_4_7, 3);

    sum0 = vqrshlq_s32(sum0, shift_s32);
    sum1 = vqrshlq_s32(sum1, shift_s32);

    return vcombine_u16(vqmovun_s32(sum0), vqmovun_s32(sum1));
}

static INLINE void highbd_convolve_2d_sr_horiz_neon(const uint16_t *src_ptr, int src_stride, uint16_t *dst_ptr,
                                                    int dst_stride, int w, int h, const int16_t *x_filter_ptr,
                                                    ConvolveParams *conv_params, const int offset) {
    // The smallest block height processed by the SIMD functions is 4, and the
    // horizontal convolution needs to process an extra (filter_taps/2 - 1) lines
    // for the vertical convolution.
    assert(h >= 5);
    const int32x4_t shift_s32  = vdupq_n_s32(-conv_params->round_0);
    const int32x4_t offset_s32 = vdupq_n_s32(offset);

    if (w == 4) {
        // 4-tap filters are used for blocks having width <= 4.
        const int16x4_t x_filter = vld1_s16(x_filter_ptr + 2);
        const int16_t  *s        = (const int16_t *)(src_ptr + 1);
        uint16_t       *d        = dst_ptr;

        do {
            int16x4_t s0[4], s1[4], s2[4], s3[4];
            load_s16_4x4(s + 0 * src_stride, 1, &s0[0], &s0[1], &s0[2], &s0[3]);
            load_s16_4x4(s + 1 * src_stride, 1, &s1[0], &s1[1], &s1[2], &s1[3]);
            load_s16_4x4(s + 2 * src_stride, 1, &s2[0], &s2[1], &s2[2], &s2[3]);
            load_s16_4x4(s + 3 * src_stride, 1, &s3[0], &s3[1], &s3[2], &s3[3]);

            uint16x4_t d0 = highbd_convolve4_4_2d_h(s0, x_filter, shift_s32, offset_s32);
            uint16x4_t d1 = highbd_convolve4_4_2d_h(s1, x_filter, shift_s32, offset_s32);
            uint16x4_t d2 = highbd_convolve4_4_2d_h(s2, x_filter, shift_s32, offset_s32);
            uint16x4_t d3 = highbd_convolve4_4_2d_h(s3, x_filter, shift_s32, offset_s32);

            store_u16_4x4(d, dst_stride, d0, d1, d2, d3);

            s += 4 * src_stride;
            d += 4 * dst_stride;
            h -= 4;
        } while (h > 4);

        do {
            int16x4_t s0[4];
            load_s16_4x4(s, 1, &s0[0], &s0[1], &s0[2], &s0[3]);

            uint16x4_t d0 = highbd_convolve4_4_2d_h(s0, x_filter, shift_s32, offset_s32);

            vst1_u16(d, d0);

            s += src_stride;
            d += dst_stride;
        } while (--h != 0);
    } else {
        const int16x8_t x_filter = vld1q_s16(x_filter_ptr);
        int             height   = h;

        do {
            int            width = w;
            const int16_t *s     = (const int16_t *)src_ptr;
            uint16_t      *d     = dst_ptr;

            do {
                int16x8_t s0[8], s1[8], s2[8], s3[8];
                load_s16_8x8(s + 0 * src_stride, 1, &s0[0], &s0[1], &s0[2], &s0[3], &s0[4], &s0[5], &s0[6], &s0[7]);
                load_s16_8x8(s + 1 * src_stride, 1, &s1[0], &s1[1], &s1[2], &s1[3], &s1[4], &s1[5], &s1[6], &s1[7]);
                load_s16_8x8(s + 2 * src_stride, 1, &s2[0], &s2[1], &s2[2], &s2[3], &s2[4], &s2[5], &s2[6], &s2[7]);
                load_s16_8x8(s + 3 * src_stride, 1, &s3[0], &s3[1], &s3[2], &s3[3], &s3[4], &s3[5], &s3[6], &s3[7]);

                uint16x8_t d0 = highbd_convolve8_8_2d_h(s0, x_filter, shift_s32, offset_s32);
                uint16x8_t d1 = highbd_convolve8_8_2d_h(s1, x_filter, shift_s32, offset_s32);
                uint16x8_t d2 = highbd_convolve8_8_2d_h(s2, x_filter, shift_s32, offset_s32);
                uint16x8_t d3 = highbd_convolve8_8_2d_h(s3, x_filter, shift_s32, offset_s32);

                store_u16_8x4(d, dst_stride, d0, d1, d2, d3);

                s += 8;
                d += 8;
                width -= 8;
            } while (width != 0);
            src_ptr += 4 * src_stride;
            dst_ptr += 4 * dst_stride;
            height -= 4;
        } while (height > 4);

        do {
            int            width = w;
            const int16_t *s     = (const int16_t *)src_ptr;
            uint16_t      *d     = dst_ptr;

            do {
                int16x8_t s0[8];
                load_s16_8x8(s + 0 * src_stride, 1, &s0[0], &s0[1], &s0[2], &s0[3], &s0[4], &s0[5], &s0[6], &s0[7]);

                uint16x8_t d0 = highbd_convolve8_8_2d_h(s0, x_filter, shift_s32, offset_s32);
                vst1q_u16(d, d0);

                s += 8;
                d += 8;
                width -= 8;
            } while (width != 0);
            src_ptr += src_stride;
            dst_ptr += dst_stride;
        } while (--height != 0);
    }
}

void svt_av1_highbd_convolve_2d_sr_neon(const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
                                        int h, const InterpFilterParams *filter_params_x,
                                        const InterpFilterParams *filter_params_y, const int subpel_x_qn,
                                        const int subpel_y_qn, ConvolveParams *conv_params, int bd) {
    if (w == 2 || h == 2) {
        svt_av1_highbd_convolve_2d_sr_c(src,
                                        src_stride,
                                        dst,
                                        dst_stride,
                                        w,
                                        h,
                                        filter_params_x,
                                        filter_params_y,
                                        subpel_x_qn,
                                        subpel_y_qn,
                                        conv_params,
                                        bd);
        return;
    }
    DECLARE_ALIGNED(16, uint16_t, im_block[(MAX_SB_SIZE + MAX_FILTER_TAP) * MAX_SB_SIZE]);
    const int x_filter_taps  = get_filter_tap(filter_params_x, subpel_x_qn);
    const int clamped_x_taps = x_filter_taps < 6 ? 6 : x_filter_taps;

    const int y_filter_taps    = get_filter_tap(filter_params_y, subpel_y_qn);
    const int clamped_y_taps   = y_filter_taps < 6 ? 6 : y_filter_taps;
    const int im_h             = h + clamped_y_taps - 1;
    const int im_stride        = MAX_SB_SIZE;
    const int vert_offset      = clamped_y_taps / 2 - 1;
    const int horiz_offset     = clamped_x_taps / 2 - 1;
    const int x_offset_initial = (1 << (bd + FILTER_BITS - 1));
    const int y_offset_bits    = bd + 2 * FILTER_BITS - conv_params->round_0;
    // The extra shim of (1 << (conv_params->round_1 - 1)) allows us to do a
    // simple shift left instead of a rounding saturating shift left.
    const int y_offset = (1 << (conv_params->round_1 - 1)) - (1 << (y_offset_bits - 1));

    const uint16_t *src_ptr = src - vert_offset * src_stride - horiz_offset;

    const int16_t *x_filter_ptr = av1_get_interp_filter_subpel_kernel(*filter_params_x, subpel_x_qn & SUBPEL_MASK);
    const int16_t *y_filter_ptr = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_qn & SUBPEL_MASK);

    if (x_filter_taps <= 6 && w != 4) {
        highbd_convolve_2d_sr_horiz_6tap_neon(
            src_ptr, src_stride, im_block, im_stride, w, im_h, x_filter_ptr, conv_params, x_offset_initial);
    } else {
        highbd_convolve_2d_sr_horiz_neon(
            src_ptr, src_stride, im_block, im_stride, w, im_h, x_filter_ptr, conv_params, x_offset_initial);
    }

    if (y_filter_taps <= 6) {
        highbd_convolve_2d_sr_vert_6tap_neon(
            im_block, im_stride, dst, dst_stride, w, h, y_filter_ptr, conv_params, bd, y_offset);
    } else {
        highbd_convolve_2d_sr_vert_8tap_neon(
            im_block, im_stride, dst, dst_stride, w, h, y_filter_ptr, conv_params, bd, y_offset);
    }
}
