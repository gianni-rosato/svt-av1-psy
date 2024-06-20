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
#include "mem_neon.h"
#include "transpose_neon.h"

static INLINE void residual_kernel_neon(const uint8_t *restrict input, const uint8_t *restrict pred,
                                        int16x8_t *residual) {
    const uint8x8_t in = vld1_u8(input);
    const uint8x8_t pr = vld1_u8(pred);
    *residual          = vreinterpretq_s16_u16(vsubl_u8(in, pr));
}

static INLINE void hadamard_4x4_one_pass(int16x4_t *a0, int16x4_t *a1, int16x4_t *a2, int16x4_t *a3) {
    const int16x4_t b0 = vhadd_s16(*a0, *a1);
    const int16x4_t b1 = vhsub_s16(*a0, *a1);
    const int16x4_t b2 = vhadd_s16(*a2, *a3);
    const int16x4_t b3 = vhsub_s16(*a2, *a3);

    *a0 = vadd_s16(b0, b2);
    *a1 = vadd_s16(b1, b3);
    *a2 = vsub_s16(b0, b2);
    *a3 = vsub_s16(b1, b3);
}

static INLINE int hadamard_path_4x4_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                         const uint8_t *restrict pred, const uint32_t  pred_stride) {
    uint32x2_t in, pr;
    int16x8_t  residual[2];
    int32x4_t  accum;
    int16x4_t  a0, a1, a2, a3;

    /* Residual Kernel */
    in = vdup_n_u32(*(uint32_t *)(input + 0 * input_stride));
    in = vset_lane_u32(*(uint32_t *)(input + 1 * input_stride), in, 1);
    pr = vdup_n_u32(*(uint32_t *)(pred + 0 * pred_stride));
    pr = vset_lane_u32(*(uint32_t *)(pred + 1 * pred_stride), pr, 1);

    residual[0] = vreinterpretq_s16_u16(vsubl_u8(vreinterpret_u8_u32(in), vreinterpret_u8_u32(pr)));

    input += 2 * input_stride;
    pred += 2 * pred_stride;
    in = vdup_n_u32(*(uint32_t *)(input + 0 * input_stride));
    in = vset_lane_u32(*(uint32_t *)(input + 1 * input_stride), in, 1);
    pr = vdup_n_u32(*(uint32_t *)(pred + 0 * pred_stride));
    pr = vset_lane_u32(*(uint32_t *)(pred + 1 * pred_stride), pr, 1);

    residual[1] = vreinterpretq_s16_u16(vsubl_u8(vreinterpret_u8_u32(in), vreinterpret_u8_u32(pr)));

    /* Hadamard */
    a0 = vget_low_s16(residual[0]);
    a1 = vget_high_s16(residual[0]);
    a2 = vget_low_s16(residual[1]);
    a3 = vget_high_s16(residual[1]);

    hadamard_4x4_one_pass(&a0, &a1, &a2, &a3);

    transpose_elems_inplace_s16_4x4(&a0, &a1, &a2, &a3);
    hadamard_4x4_one_pass(&a0, &a1, &a2, &a3);

    /* SATD */
    accum = vpaddlq_s16(vabsq_s16(vcombine_s16(a0, a1)));
    accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(vcombine_s16(a2, a3))));
    return vaddvq_s32(accum);
}

static INLINE void hadamard8x8_one_pass(int16x8_t *a) {
    const int16x8_t b0 = vaddq_s16(a[0], a[1]);
    const int16x8_t b1 = vsubq_s16(a[0], a[1]);
    const int16x8_t b2 = vaddq_s16(a[2], a[3]);
    const int16x8_t b3 = vsubq_s16(a[2], a[3]);
    const int16x8_t b4 = vaddq_s16(a[4], a[5]);
    const int16x8_t b5 = vsubq_s16(a[4], a[5]);
    const int16x8_t b6 = vaddq_s16(a[6], a[7]);
    const int16x8_t b7 = vsubq_s16(a[6], a[7]);

    const int16x8_t c0 = vaddq_s16(b0, b2);
    const int16x8_t c1 = vaddq_s16(b1, b3);
    const int16x8_t c2 = vsubq_s16(b0, b2);
    const int16x8_t c3 = vsubq_s16(b1, b3);
    const int16x8_t c4 = vaddq_s16(b4, b6);
    const int16x8_t c5 = vaddq_s16(b5, b7);
    const int16x8_t c6 = vsubq_s16(b4, b6);
    const int16x8_t c7 = vsubq_s16(b5, b7);

    a[0] = vaddq_s16(c0, c4);
    a[1] = vsubq_s16(c2, c6);
    a[2] = vsubq_s16(c0, c4);
    a[3] = vaddq_s16(c2, c6);
    a[4] = vaddq_s16(c3, c7);
    a[5] = vsubq_s16(c3, c7);
    a[6] = vsubq_s16(c1, c5);
    a[7] = vaddq_s16(c1, c5);
}

static INLINE int hadamard_path_8x8_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                         const uint8_t *restrict pred, const uint32_t  pred_stride) {
    int16x8_t  residual[8];
    int16x8_t *a = residual;
    int32x4_t  accum;

    /* Residual Kernel */
    int i = 0;
    do {
        residual_kernel_neon(input, pred, residual + i);
        input += input_stride;
        pred += pred_stride;
        i++;
    } while (i != 8);

    /* Hadamard Path */
    hadamard8x8_one_pass(a);
    transpose_elems_inplace_s16_8x8(a + 0, a + 1, a + 2, a + 3, a + 4, a + 5, a + 6, a + 7);
    hadamard8x8_one_pass(a);

    /* SATD */
    accum = vpaddlq_s16(vabsq_s16(a[0]));
    accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(a[1])));
    accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(a[2])));
    accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(a[3])));
    accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(a[4])));
    accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(a[5])));
    accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(a[6])));
    accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(a[7])));
    return vaddvq_s32(accum);
}

static INLINE void residual_hadamard_8x8_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                              const uint8_t *restrict pred, const uint32_t  pred_stride,
                                              int16x8_t *coeff) {
    int16x8_t  residual[8];
    int16x8_t *a = residual;

    /* Residual Kernel */
    int i = 0;
    do {
        residual_kernel_neon(input, pred, residual + i);
        input += input_stride;
        pred += pred_stride;
        i++;
    } while (i != 8);

    /* Hadamard */
    hadamard8x8_one_pass(a);
    transpose_elems_inplace_s16_8x8(a + 0, a + 1, a + 2, a + 3, a + 4, a + 5, a + 6, a + 7);
    hadamard8x8_one_pass(a);

    coeff[0] = a[0];
    coeff[1] = a[1];
    coeff[2] = a[2];
    coeff[3] = a[3];
    coeff[4] = a[4];
    coeff[5] = a[5];
    coeff[6] = a[6];
    coeff[7] = a[7];
}

static INLINE void residual_hadamard_16x16_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                                const uint8_t *restrict pred, const uint32_t  pred_stride,
                                                int16x8_t *coeff) {
    int idx;

    for (idx = 0; idx < 4; ++idx) {
        // src_diff: 9 bit, dynamic range [-255, 255]
        const uint8_t *input_ptr = input + (idx >> 1) * input_stride * 8 + (idx & 0x01) * 8;
        const uint8_t *pred_ptr  = pred + (idx >> 1) * pred_stride * 8 + (idx & 0x01) * 8;
        residual_hadamard_8x8_neon(input_ptr, input_stride, pred_ptr, pred_stride, coeff + idx * 8);
    }

    /* Generate 8x4=32 coefficient values in one iteration. With 8 iterations,
    a total of 256 coefficient values will be generated. */
    for (idx = 0; idx < 64; idx += 8) {
        const int16x8_t a0 = coeff[0];
        const int16x8_t a1 = coeff[8];
        const int16x8_t a2 = coeff[16];
        const int16x8_t a3 = coeff[24];

        const int16x8_t b0 = vshrq_n_s16(vaddq_s16(a0, a1), 1);
        const int16x8_t b1 = vshrq_n_s16(vsubq_s16(a0, a1), 1);
        const int16x8_t b2 = vshrq_n_s16(vaddq_s16(a2, a3), 1);
        const int16x8_t b3 = vshrq_n_s16(vsubq_s16(a2, a3), 1);

        const int16x8_t c0 = vaddq_s16(b0, b2);
        const int16x8_t c1 = vaddq_s16(b1, b3);
        const int16x8_t c2 = vsubq_s16(b0, b2);
        const int16x8_t c3 = vsubq_s16(b1, b3);

        coeff[0]  = c0;
        coeff[8]  = c1;
        coeff[16] = c2;
        coeff[24] = c3;
        coeff += 1;
    }
}

static INLINE int hadamard_path_16x16_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                           const uint8_t *restrict pred, const uint32_t  pred_stride) {
    int        idx;
    int16x8_t  coeff_buf[32];
    int32x4_t  accum = vdupq_n_s32(0);
    int16x8_t *coeff = coeff_buf;

    for (idx = 0; idx < 4; ++idx) {
        // src_diff: 9 bit, dynamic range [-255, 255]
        const uint8_t *input_ptr = input + (idx >> 1) * input_stride * 8 + (idx & 0x01) * 8;
        const uint8_t *pred_ptr  = pred + (idx >> 1) * pred_stride * 8 + (idx & 0x01) * 8;
        residual_hadamard_8x8_neon(input_ptr, input_stride, pred_ptr, pred_stride, coeff + idx * 8);
    }

    /* Generate 8x4=32 coefficient values in one iteration. With 8 iterations,
    a total of 256 coefficient values will be generated. */
    for (idx = 0; idx < 64; idx += 8) {
        const int16x8_t a0 = coeff[0];
        const int16x8_t a1 = coeff[8];
        const int16x8_t a2 = coeff[16];
        const int16x8_t a3 = coeff[24];

        const int16x8_t b0 = vshrq_n_s16(vaddq_s16(a0, a1), 1);
        const int16x8_t b1 = vshrq_n_s16(vsubq_s16(a0, a1), 1);
        const int16x8_t b2 = vshrq_n_s16(vaddq_s16(a2, a3), 1);
        const int16x8_t b3 = vshrq_n_s16(vsubq_s16(a2, a3), 1);

        const int16x8_t c0 = vaddq_s16(b0, b2);
        const int16x8_t c1 = vaddq_s16(b1, b3);
        const int16x8_t c2 = vsubq_s16(b0, b2);
        const int16x8_t c3 = vsubq_s16(b1, b3);

        /* SATD */
        accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(c0)));
        accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(c1)));
        accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(c2)));
        accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(c3)));
        coeff += 1;
    }
    return vaddvq_s32(accum);
}

static INLINE int hadamard_path_32x32_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                           const uint8_t *restrict pred, const uint32_t  pred_stride) {
    int        idx;
    int16x8_t  coeff_buf[128];
    int32x4_t  accum = vdupq_n_s32(0);
    int16x8_t *coeff = coeff_buf;

    /* Hadamard */
    for (idx = 0; idx < 4; ++idx) {
        // src_diff: 9 bit, dynamic range [-255, 255]
        const uint8_t *input_ptr = input + (idx >> 1) * 16 * input_stride + (idx & 0x01) * 16;
        const uint8_t *pred_ptr  = pred + (idx >> 1) * pred_stride * 16 + (idx & 0x01) * 16;
        residual_hadamard_16x16_neon(input_ptr, input_stride, pred_ptr, pred_stride, coeff + idx * 32);
    }

    for (idx = 0; idx < 256; idx += 8) {
        const int16x8_t a0 = coeff[0];
        const int16x8_t a1 = coeff[32];
        const int16x8_t a2 = coeff[64];
        const int16x8_t a3 = coeff[96];

        const int16x8_t b0 = vshrq_n_s16(vaddq_s16(a0, a1), 2);
        const int16x8_t b1 = vshrq_n_s16(vsubq_s16(a0, a1), 2);
        const int16x8_t b2 = vshrq_n_s16(vaddq_s16(a2, a3), 2);
        const int16x8_t b3 = vshrq_n_s16(vsubq_s16(a2, a3), 2);

        const int16x8_t c0 = vaddq_s16(b0, b2);
        const int16x8_t c1 = vaddq_s16(b1, b3);
        const int16x8_t c2 = vsubq_s16(b0, b2);
        const int16x8_t c3 = vsubq_s16(b1, b3);

        /* SATD */
        accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(c0)));
        accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(c1)));
        accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(c2)));
        accum = vaddq_s32(accum, vpaddlq_s16(vabsq_s16(c3)));
        coeff += 1;
    }
    return vaddvq_s32(accum);
}

uint32_t hadamard_path_neon(Buf2D residualBuf, Buf2D coeffBuf, Buf2D inputBuf, Buf2D predBuf, BlockSize bsize) {
    assert(residualBuf.buf != NULL && residualBuf.buf0 == NULL && residualBuf.width == 0 && residualBuf.height == 0 &&
           residualBuf.stride != 0);
    assert(coeffBuf.buf != NULL && coeffBuf.buf0 == NULL && coeffBuf.width == 0 && coeffBuf.height == 0 &&
           coeffBuf.stride == block_size_wide[bsize]);
    assert(inputBuf.buf != NULL && inputBuf.buf0 == NULL && inputBuf.width == 0 && inputBuf.height == 0 &&
           inputBuf.stride != 0);
    assert(predBuf.buf != NULL && predBuf.buf0 == NULL && predBuf.width == 0 && predBuf.height == 0 &&
           predBuf.stride != 0);
    (void)residualBuf;
    (void)coeffBuf;
    uint32_t input_idx, pred_idx;
    uint32_t satd_cost = 0;

    const TxSize tx_size = AOMMIN(TX_32X32, max_txsize_lookup[bsize]);

    const int stepr = tx_size_high_unit[tx_size];
    const int stepc = tx_size_wide_unit[tx_size];

    const int max_blocks_wide = block_size_wide[bsize] >> MI_SIZE_LOG2;
    const int max_blocks_high = block_size_wide[bsize] >> MI_SIZE_LOG2;
    int       row, col;

    for (row = 0; row < max_blocks_high; row += stepr) {
        for (col = 0; col < max_blocks_wide; col += stepc) {
            input_idx = ((row * inputBuf.stride) + col) << 2;
            pred_idx  = ((row * predBuf.stride) + col) << 2;

            switch (tx_size) {
            case TX_4X4: {
                satd_cost += hadamard_path_4x4_neon(
                    inputBuf.buf + input_idx, inputBuf.stride, predBuf.buf + pred_idx, predBuf.stride);
                break;
            }
            case TX_8X8: {
                satd_cost += hadamard_path_8x8_neon(
                    inputBuf.buf + input_idx, inputBuf.stride, predBuf.buf + pred_idx, predBuf.stride);
                break;
            }
            case TX_16X16: {
                satd_cost += hadamard_path_16x16_neon(
                    inputBuf.buf + input_idx, inputBuf.stride, predBuf.buf + pred_idx, predBuf.stride);
                break;
            }
            case TX_32X32: {
                satd_cost += hadamard_path_32x32_neon(
                    inputBuf.buf + input_idx, inputBuf.stride, predBuf.buf + pred_idx, predBuf.stride);
                break;
            }
            default: assert(0);
            }
        }
    }

    return (satd_cost);
}
