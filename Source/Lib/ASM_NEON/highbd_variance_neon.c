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

#include "aom_dsp_rtcd.h"

void svt_aom_highbd_8_mse16x16_neon(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride,
                                    uint32_t *sse) {
    uint16_t *src = CONVERT_TO_SHORTPTR(src_ptr);
    uint16_t *ref = CONVERT_TO_SHORTPTR(ref_ptr);

    uint32x4_t sse_u32[2] = {vdupq_n_u32(0), vdupq_n_u32(0)};

    int i = 16;
    do {
        int j = 0;
        do {
            uint16x8_t s = vld1q_u16(src + j);
            uint16x8_t r = vld1q_u16(ref + j);

            uint16x8_t diff = vabdq_u16(s, r);

            sse_u32[0] = vmlal_u16(sse_u32[0], vget_low_u16(diff), vget_low_u16(diff));
            sse_u32[1] = vmlal_u16(sse_u32[1], vget_high_u16(diff), vget_high_u16(diff));

            j += 8;
        } while (j < 16);

        src += src_stride;
        ref += ref_stride;
    } while (--i != 0);

    *sse = vaddvq_u32(vaddq_u32(sse_u32[0], sse_u32[1]));
}
