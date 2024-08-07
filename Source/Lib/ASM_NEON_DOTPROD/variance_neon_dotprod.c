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

unsigned int svt_aom_mse16x16_neon_dotprod(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                           unsigned int *sse) {
    uint32x4_t sse_u32[2] = {vdupq_n_u32(0), vdupq_n_u32(0)};

    int h = 16;
    do {
        uint8x16_t s0 = vld1q_u8(src);
        uint8x16_t s1 = vld1q_u8(src + src_stride);
        uint8x16_t r0 = vld1q_u8(ref);
        uint8x16_t r1 = vld1q_u8(ref + ref_stride);

        uint8x16_t abs_diff0 = vabdq_u8(s0, r0);
        uint8x16_t abs_diff1 = vabdq_u8(s1, r1);

        sse_u32[0] = vdotq_u32(sse_u32[0], abs_diff0, abs_diff0);
        sse_u32[1] = vdotq_u32(sse_u32[1], abs_diff1, abs_diff1);

        src += 2 * src_stride;
        ref += 2 * ref_stride;
        h -= 2;
    } while (h != 0);

    *sse = vaddvq_u32(vaddq_u32(sse_u32[0], sse_u32[1]));
    return *sse;
}
