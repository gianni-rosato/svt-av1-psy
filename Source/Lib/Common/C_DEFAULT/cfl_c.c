/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"

#define ROUND_POWER_OF_TWO_SIGNED(value, n) \
    (((value) < 0) ? -ROUND_POWER_OF_TWO(-(value), (n)) : ROUND_POWER_OF_TWO((value), (n)))

static INLINE int32_t get_scaled_luma_q0(int32_t alpha_q3, int16_t pred_buf_q3) {
int32_t scaled_luma_q6 = alpha_q3 * pred_buf_q3;
return ROUND_POWER_OF_TWO_SIGNED(scaled_luma_q6, 6);
}

void eb_cfl_predict_lbd_c(
        const int16_t *pred_buf_q3,
        uint8_t *pred,// AMIR ADDED
        int32_t pred_stride,
        uint8_t *dst,// AMIR changed to 8 bit
        int32_t dst_stride,
        int32_t alpha_q3,
        int32_t bit_depth,
        int32_t width,
        int32_t height) {
    for (int32_t j = 0; j < height; j++) {
        for (int32_t i = 0; i < width; i++) {
            dst[i] = (uint8_t)clip_pixel_highbd(
                    get_scaled_luma_q0(alpha_q3, pred_buf_q3[i]) + (int16_t)pred[i], bit_depth);
        }
        dst += dst_stride;
        pred += pred_stride;
        pred_buf_q3 += CFL_BUF_LINE;
    }
}

void eb_cfl_predict_hbd_c(
        const int16_t *pred_buf_q3,
        uint16_t *pred,// AMIR ADDED
        int32_t pred_stride,
        uint16_t *dst,// AMIR changed to 8 bit
        int32_t dst_stride,
        int32_t alpha_q3,
        int32_t bit_depth,
        int32_t width,
        int32_t height) {
    for (int32_t j = 0; j < height; j++) {
        for (int32_t i = 0; i < width; i++) {
            dst[i] = clip_pixel_highbd(
                    get_scaled_luma_q0(alpha_q3, pred_buf_q3[i]) + (int16_t)pred[i], bit_depth);
        }
        dst += dst_stride;
        pred += pred_stride;
        pred_buf_q3 += CFL_BUF_LINE;
    }
}
