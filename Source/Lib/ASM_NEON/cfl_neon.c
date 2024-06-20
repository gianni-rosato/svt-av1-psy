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
#include "definitions.h"

/* Store half of a vector. */
static INLINE void vsth_u8(uint8_t *ptr, uint8x8_t val) { vst1_lane_u32((uint32_t *)ptr, vreinterpret_u32_u8(val), 0); }

/* Saturating negate 16-bit integers in a when the corresponding signed 16-bit
integer in b is negative.
Notes:
  * Negating INT16_MIN results in INT16_MIN. However, this cannot occur in
  practice, as scaled_luma is the multiplication of two absolute values.
  * In the Intel equivalent, elements in a are zeroed out when the
  corresponding elements in b are zero. Because vsign is used twice in a
  row, with b in the first call becoming a in the second call, there's no
  impact from not zeroing out. */
static int16x4_t vsign_s16(int16x4_t a, int16x4_t b) {
    const int16x4_t mask = vshr_n_s16(b, 15);
    return veor_s16(vadd_s16(a, mask), mask);
}

/* Saturating negate 16-bit integers in a when the corresponding signed 16-bit
integer in b is negative.
Notes:
  * Negating INT16_MIN results in INT16_MIN. However, this cannot occur in
  practice, as scaled_luma is the multiplication of two absolute values.
  * In the Intel equivalent, elements in a are zeroed out when the
  corresponding elements in b are zero. Because vsignq is used twice in a
  row, with b in the first call becoming a in the second call, there's no
  impact from not zeroing out. */
static int16x8_t vsignq_s16(int16x8_t a, int16x8_t b) {
    const int16x8_t mask = vshrq_n_s16(b, 15);
    return veorq_s16(vaddq_s16(a, mask), mask);
}

static INLINE int16x4_t predict_w4(const int16_t *pred_buf_q3, int16x4_t alpha_sign, int abs_alpha_q12, int16x4_t dc) {
    const int16x4_t ac_q3       = vld1_s16(pred_buf_q3);
    const int16x4_t ac_sign     = veor_s16(alpha_sign, ac_q3);
    int16x4_t       scaled_luma = vqrdmulh_n_s16(vabs_s16(ac_q3), abs_alpha_q12);
    return vadd_s16(vsign_s16(scaled_luma, ac_sign), dc);
}

static INLINE int16x8_t predict_w8(const int16_t *pred_buf_q3, int16x8_t alpha_sign, int abs_alpha_q12, int16x8_t dc) {
    const int16x8_t ac_q3       = vld1q_s16(pred_buf_q3);
    const int16x8_t ac_sign     = veorq_s16(alpha_sign, ac_q3);
    int16x8_t       scaled_luma = vqrdmulhq_n_s16(vabsq_s16(ac_q3), abs_alpha_q12);
    return vaddq_s16(vsignq_s16(scaled_luma, ac_sign), dc);
}

static INLINE int16x8x2_t predict_w16(const int16_t *pred_buf_q3, int16x8_t alpha_sign, int abs_alpha_q12,
                                      int16x8_t dc) {
    /* vld2q_s16 interleaves, which is not useful for prediction. vst1q_s16_x2
    does not interleave, but is not currently available in the compilier used
    by the AOM build system. */
    const int16x8x2_t ac_q3         = vld2q_s16(pred_buf_q3);
    const int16x8_t   ac_sign_0     = veorq_s16(alpha_sign, ac_q3.val[0]);
    const int16x8_t   ac_sign_1     = veorq_s16(alpha_sign, ac_q3.val[1]);
    const int16x8_t   scaled_luma_0 = vqrdmulhq_n_s16(vabsq_s16(ac_q3.val[0]), abs_alpha_q12);
    const int16x8_t   scaled_luma_1 = vqrdmulhq_n_s16(vabsq_s16(ac_q3.val[1]), abs_alpha_q12);
    int16x8x2_t       result;
    result.val[0] = vaddq_s16(vsignq_s16(scaled_luma_0, ac_sign_0), dc);
    result.val[1] = vaddq_s16(vsignq_s16(scaled_luma_1, ac_sign_1), dc);
    return result;
}

static INLINE int16x8x4_t predict_w32(const int16_t *pred_buf_q3, int16x8_t alpha_sign, int abs_alpha_q12,
                                      int16x8_t dc) {
    /* vld4q_s16 interleaves, which is not useful for prediction. vst1q_s16_x4
    does not interleave, but is not currently available in the compilier used
    by the AOM build system. */
    const int16x8x4_t ac_q3         = vld4q_s16(pred_buf_q3);
    const int16x8_t   ac_sign_0     = veorq_s16(alpha_sign, ac_q3.val[0]);
    const int16x8_t   ac_sign_1     = veorq_s16(alpha_sign, ac_q3.val[1]);
    const int16x8_t   ac_sign_2     = veorq_s16(alpha_sign, ac_q3.val[2]);
    const int16x8_t   ac_sign_3     = veorq_s16(alpha_sign, ac_q3.val[3]);
    const int16x8_t   scaled_luma_0 = vqrdmulhq_n_s16(vabsq_s16(ac_q3.val[0]), abs_alpha_q12);
    const int16x8_t   scaled_luma_1 = vqrdmulhq_n_s16(vabsq_s16(ac_q3.val[1]), abs_alpha_q12);
    const int16x8_t   scaled_luma_2 = vqrdmulhq_n_s16(vabsq_s16(ac_q3.val[2]), abs_alpha_q12);
    const int16x8_t   scaled_luma_3 = vqrdmulhq_n_s16(vabsq_s16(ac_q3.val[3]), abs_alpha_q12);
    int16x8x4_t       result;
    result.val[0] = vaddq_s16(vsignq_s16(scaled_luma_0, ac_sign_0), dc);
    result.val[1] = vaddq_s16(vsignq_s16(scaled_luma_1, ac_sign_1), dc);
    result.val[2] = vaddq_s16(vsignq_s16(scaled_luma_2, ac_sign_2), dc);
    result.val[3] = vaddq_s16(vsignq_s16(scaled_luma_3, ac_sign_3), dc);
    return result;
}

void svt_aom_cfl_predict_lbd_neon(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst,
                                  int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width,
                                  int32_t height) {
    (void)bit_depth;
    (void)pred_stride;
    const int16_t        abs_alpha_q12 = abs(alpha_q3) << 9;
    const int16_t *const end           = pred_buf_q3 + height * CFL_BUF_LINE;
    if (width == 4) {
        const int16x4_t alpha_sign = vdup_n_s16(alpha_q3);
        const int16x4_t dc         = vdup_n_s16(*pred);
        do {
            const int16x4_t pred_vector = predict_w4(pred_buf_q3, alpha_sign, abs_alpha_q12, dc);
            vsth_u8(dst, vqmovun_s16(vcombine_s16(pred_vector, pred_vector)));
            dst += dst_stride;
        } while ((pred_buf_q3 += CFL_BUF_LINE) < end);
    } else {
        const int16x8_t alpha_sign = vdupq_n_s16(alpha_q3);
        const int16x8_t dc         = vdupq_n_s16(*pred);
        do {
            if (width == 8) {
                vst1_u8(dst, vqmovun_s16(predict_w8(pred_buf_q3, alpha_sign, abs_alpha_q12, dc)));
            } else if (width == 16) {
                const int16x8x2_t pred_vector = predict_w16(pred_buf_q3, alpha_sign, abs_alpha_q12, dc);
                const uint8x8x2_t predun      = {{vqmovun_s16(pred_vector.val[0]), vqmovun_s16(pred_vector.val[1])}};
                vst2_u8(dst, predun);
            } else {
                const int16x8x4_t pred_vector = predict_w32(pred_buf_q3, alpha_sign, abs_alpha_q12, dc);
                const uint8x8x4_t predun      = {{vqmovun_s16(pred_vector.val[0]),
                                                  vqmovun_s16(pred_vector.val[1]),
                                                  vqmovun_s16(pred_vector.val[2]),
                                                  vqmovun_s16(pred_vector.val[3])}};
                vst4_u8(dst, predun);
            }
            dst += dst_stride;
        } while ((pred_buf_q3 += CFL_BUF_LINE) < end);
    }
}
