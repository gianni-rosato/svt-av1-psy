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
#include "mem_neon.h"

void svt_av1_convolve_2d_copy_sr_neon(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride,
                                      int32_t w, int32_t h, InterpFilterParams *filter_params_x,
                                      InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
                                      const int32_t subpel_y_q4, ConvolveParams *conv_params) {
    (void)filter_params_x;
    (void)filter_params_y;
    (void)subpel_x_q4;
    (void)subpel_y_q4;
    (void)conv_params;

    const uint8_t *src1;
    uint8_t       *dst1;
    int            y;

    if (!(w & 0x0F)) {
        for (y = 0; y < h; ++y) {
            src1 = src;
            dst1 = dst;
            for (int x = 0; x < (w >> 4); ++x) {
                vst1q_u8(dst1, vld1q_u8(src1));
                src1 += 16;
                dst1 += 16;
            }
            src += src_stride;
            dst += dst_stride;
        }
    } else if (!(w & 0x07)) {
        for (y = 0; y < h; ++y) {
            vst1_u8(dst, vld1_u8(src));
            src += src_stride;
            dst += dst_stride;
        }
    } else if (!(w & 0x03)) {
        for (y = 0; y < h; ++y) {
            vst1_lane_u32((uint32_t *)(dst), vreinterpret_u32_u8(vld1_u8(src)), 0);
            src += src_stride;
            dst += dst_stride;
        }
    } else if (!(w & 0x01)) {
        for (y = 0; y < h; ++y) {
            vst1_lane_u16((uint16_t *)(dst), vreinterpret_u16_u8(vld1_u8(src)), 0);
            src += src_stride;
            dst += dst_stride;
        }
    }
}
