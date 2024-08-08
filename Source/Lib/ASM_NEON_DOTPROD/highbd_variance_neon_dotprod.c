/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved.
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>

#include "aom_dsp_rtcd.h"

void svt_aom_highbd_8_mse16x16_neon_dotprod(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr,
                                            int ref_stride, uint32_t *sse) {
    uint16_t *src = CONVERT_TO_SHORTPTR(src_ptr);
    uint16_t *ref = CONVERT_TO_SHORTPTR(ref_ptr);

    uint32x4_t sse_u32 = vdupq_n_u32(0);

    int i = 16;
    do {
        uint16x8_t s0 = vld1q_u16(src);
        uint16x8_t s1 = vld1q_u16(src + 8);
        uint16x8_t r0 = vld1q_u16(ref);
        uint16x8_t r1 = vld1q_u16(ref + 8);

        uint8x16_t s = vcombine_u8(vmovn_u16(s0), vmovn_u16(s1));
        uint8x16_t r = vcombine_u8(vmovn_u16(r0), vmovn_u16(r1));

        uint8x16_t diff = vabdq_u8(s, r);
        sse_u32         = vdotq_u32(sse_u32, diff, diff);

        src += src_stride;
        ref += ref_stride;
    } while (--i != 0);

    *sse = vaddvq_u32(sse_u32);
}
