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
#include "EbDefinitions.h"

static INLINE void residual_kernel4_neon(const uint8_t *input, const uint32_t input_stride, const uint8_t *pred,
                                         const uint32_t pred_stride, int16_t *residual, const uint32_t residual_stride,
                                         const uint32_t area_height) {
    uint32x2_t in, pr;
    int64x2_t  re;

    uint32_t y = area_height;

    do {
        in = vdup_n_u32(*(uint32_t *)(input + 0 * input_stride));
        in = vset_lane_u32(*(uint32_t *)(input + 1 * input_stride), in, 1);
        pr = vdup_n_u32(*(uint32_t *)(pred + 0 * pred_stride));
        pr = vset_lane_u32(*(uint32_t *)(pred + 1 * pred_stride), pr, 1);

        re = vreinterpretq_s64_u16(vsubl_u8(vreinterpret_u8_u32(in), vreinterpret_u8_u32(pr)));
        vst1q_lane_s64((int64_t *)residual, re, 0);
        vst1q_lane_s64((int64_t *)(residual + residual_stride), re, 1);

        input += 2 * input_stride;
        pred += 2 * pred_stride;
        residual += 2 * residual_stride;
        y -= 2;
    } while (y);
}

static INLINE void residual_kernel_neon(const uint8_t *input, const uint8_t *pred, int16_t *residual) {
    const uint8x8_t in = vld1_u8(input);
    const uint8x8_t pr = vld1_u8(pred);
    const int16x8_t re = vreinterpretq_s16_u16(vsubl_u8(in, pr));
    vst1q_s16(residual, re);
}

static INLINE void residual_kernel8_neon(const uint8_t *input, const uint32_t input_stride, const uint8_t *pred,
                                         const uint32_t pred_stride, int16_t *residual, const uint32_t residual_stride,
                                         const uint32_t area_height) {
    uint32_t y = area_height;
    int      i;
    do {
        for (i = 0; i < 8; i += 8) { residual_kernel_neon(input + i, pred + i, residual + i); }

        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
        y -= 1;
    } while (y);
}

static INLINE void residual_kernel16_neon(const uint8_t *input, const uint32_t input_stride, const uint8_t *pred,
                                          const uint32_t pred_stride, int16_t *residual, const uint32_t residual_stride,
                                          const uint32_t area_height) {
    uint32_t y = area_height;
    int      i;
    do {
        for (i = 0; i < 16; i += 8) { residual_kernel_neon(input + i, pred + i, residual + i); }

        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
        y -= 1;
    } while (y);
}

static INLINE void residual_kernel32_neon(const uint8_t *input, const uint32_t input_stride, const uint8_t *pred,
                                          const uint32_t pred_stride, int16_t *residual, const uint32_t residual_stride,
                                          const uint32_t area_height) {
    uint32_t y = area_height;
    int      i;
    do {
        for (i = 0; i < 32; i += 8) { residual_kernel_neon(input + i, pred + i, residual + i); }

        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
        y -= 1;
    } while (y);
}

static INLINE void residual_kernel64_neon(const uint8_t *input, const uint32_t input_stride, const uint8_t *pred,
                                          const uint32_t pred_stride, int16_t *residual, const uint32_t residual_stride,
                                          const uint32_t area_height) {
    uint32_t y = area_height;
    int      i;
    do {
        for (i = 0; i < 64; i += 8) { residual_kernel_neon(input + i, pred + i, residual + i); }

        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
        y -= 1;
    } while (y);
}

static INLINE void residual_kernel128_neon(const uint8_t *input, const uint32_t input_stride, const uint8_t *pred,
                                           const uint32_t pred_stride, int16_t *residual,
                                           const uint32_t residual_stride, const uint32_t area_height) {
    uint32_t y = area_height;
    int      i;
    do {
        for (i = 0; i < 128; i += 8) { residual_kernel_neon(input + i, pred + i, residual + i); }

        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
        y -= 1;
    } while (y);
}

void svt_residual_kernel8bit_neon(uint8_t *input, uint32_t input_stride, uint8_t *pred, uint32_t pred_stride,
                                  int16_t *residual, uint32_t residual_stride, uint32_t area_width,
                                  uint32_t area_height) {
    switch (area_width) {
    case 4: {
        residual_kernel4_neon(input, input_stride, pred, pred_stride, residual, residual_stride, area_height);
        break;
    }

    case 8: {
        residual_kernel8_neon(input, input_stride, pred, pred_stride, residual, residual_stride, area_height);
        break;
    }

    case 16: {
        residual_kernel16_neon(input, input_stride, pred, pred_stride, residual, residual_stride, area_height);
        break;
    }

    case 32: {
        residual_kernel32_neon(input, input_stride, pred, pred_stride, residual, residual_stride, area_height);
        break;
    }

    case 64: {
        residual_kernel64_neon(input, input_stride, pred, pred_stride, residual, residual_stride, area_height);
        break;
    }

    default: // 128
    {
        residual_kernel128_neon(input, input_stride, pred, pred_stride, residual, residual_stride, area_height);
        break;
    }
    }
}
