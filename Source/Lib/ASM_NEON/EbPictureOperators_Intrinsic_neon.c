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
#include "mem_neon.h"
#include "transpose_neon.h"

static INLINE void residual_kernel4_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                         const uint8_t *restrict pred, const uint32_t pred_stride, int16_t *residual,
                                         const uint32_t residual_stride, const uint32_t area_height) {
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

static INLINE void residual_kernel_neon(const uint8_t *restrict input, const uint8_t *restrict pred,
                                        int16_t *residual) {
    const uint8x8_t in = vld1_u8(input);
    const uint8x8_t pr = vld1_u8(pred);
    const int16x8_t re = vreinterpretq_s16_u16(vsubl_u8(in, pr));
    vst1q_s16(residual, re);
}

static INLINE void residual_kernel8_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                         const uint8_t *restrict pred, const uint32_t pred_stride, int16_t *residual,
                                         const uint32_t residual_stride, const uint32_t area_height) {
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

static INLINE void residual_kernel16_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                          const uint8_t *restrict pred, const uint32_t pred_stride, int16_t *residual,
                                          const uint32_t residual_stride, const uint32_t area_height) {
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

static INLINE void residual_kernel32_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                          const uint8_t *restrict pred, const uint32_t pred_stride, int16_t *residual,
                                          const uint32_t residual_stride, const uint32_t area_height) {
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

static INLINE void residual_kernel64_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                          const uint8_t *restrict pred, const uint32_t pred_stride, int16_t *residual,
                                          const uint32_t residual_stride, const uint32_t area_height) {
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

static INLINE void residual_kernel128_neon(const uint8_t *restrict input, const uint32_t input_stride,
                                           const uint8_t *restrict pred, const uint32_t pred_stride, int16_t *residual,
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

static void hadamard8x8_one_pass(int16x8_t *a) {
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

void svt_aom_hadamard_8x8_neon(const int16_t *src_diff, ptrdiff_t src_stride, int32_t *coeff) {
    int16x8_t a[8];

    a[0] = vld1q_s16(src_diff);
    a[1] = vld1q_s16(src_diff + src_stride);
    a[2] = vld1q_s16(src_diff + 2 * src_stride);
    a[3] = vld1q_s16(src_diff + 3 * src_stride);
    a[4] = vld1q_s16(src_diff + 4 * src_stride);
    a[5] = vld1q_s16(src_diff + 5 * src_stride);
    a[6] = vld1q_s16(src_diff + 6 * src_stride);
    a[7] = vld1q_s16(src_diff + 7 * src_stride);

    hadamard8x8_one_pass(a);
    transpose_elems_inplace_s16_8x8(a + 0, a + 1, a + 2, a + 3, a + 4, a + 5, a + 6, a + 7);
    hadamard8x8_one_pass(a);

    store_s16q_to_tran_low(coeff + 0, a[0]);
    store_s16q_to_tran_low(coeff + 8, a[1]);
    store_s16q_to_tran_low(coeff + 16, a[2]);
    store_s16q_to_tran_low(coeff + 24, a[3]);
    store_s16q_to_tran_low(coeff + 32, a[4]);
    store_s16q_to_tran_low(coeff + 40, a[5]);
    store_s16q_to_tran_low(coeff + 48, a[6]);
    store_s16q_to_tran_low(coeff + 56, a[7]);
}

void svt_aom_hadamard_16x16_neon(const int16_t *src_diff, ptrdiff_t src_stride, int32_t *coeff) {
    int idx;

    for (idx = 0; idx < 4; ++idx) {
        // src_diff: 9 bit, dynamic range [-255, 255]
        const int16_t *src_ptr = src_diff + (idx >> 1) * 8 * src_stride + (idx & 0x01) * 8;
        svt_aom_hadamard_8x8_neon(src_ptr, src_stride, coeff + idx * 64);
    }

    /* Generate 8x4=32 coefficient values in one iteration. With 8 iterations,
    a total of 256 coefficient values will be generated. */
    for (idx = 0; idx < 64; idx += 8) {
        const int16x8_t a0 = load_tran_low_to_s16q(coeff);
        const int16x8_t a1 = load_tran_low_to_s16q(coeff + 64);
        const int16x8_t a2 = load_tran_low_to_s16q(coeff + 128);
        const int16x8_t a3 = load_tran_low_to_s16q(coeff + 192);

        const int16x8_t b0 = vshrq_n_s16(vaddq_s16(a0, a1), 1);
        const int16x8_t b1 = vshrq_n_s16(vsubq_s16(a0, a1), 1);
        const int16x8_t b2 = vshrq_n_s16(vaddq_s16(a2, a3), 1);
        const int16x8_t b3 = vshrq_n_s16(vsubq_s16(a2, a3), 1);

        const int16x8_t c0 = vaddq_s16(b0, b2);
        const int16x8_t c1 = vaddq_s16(b1, b3);
        const int16x8_t c2 = vsubq_s16(b0, b2);
        const int16x8_t c3 = vsubq_s16(b1, b3);

        store_s16q_to_tran_low(coeff + 0, c0);
        store_s16q_to_tran_low(coeff + 64, c1);
        store_s16q_to_tran_low(coeff + 128, c2);
        store_s16q_to_tran_low(coeff + 192, c3);
        coeff += 8;
    }
}

void svt_aom_hadamard_32x32_neon(const int16_t *src_diff, ptrdiff_t src_stride, tran_low_t *coeff) {
    int idx;
    for (idx = 0; idx < 4; ++idx) {
        // src_diff: 9 bit, dynamic range [-255, 255]
        const int16_t *src_ptr = src_diff + (idx >> 1) * 16 * src_stride + (idx & 0x01) * 16;
        svt_aom_hadamard_16x16_neon(src_ptr, src_stride, coeff + idx * 256);
    }

    for (idx = 0; idx < 256; idx += 8) {
        const int16x8_t a0 = load_tran_low_to_s16q(coeff);
        const int16x8_t a1 = load_tran_low_to_s16q(coeff + 256);
        const int16x8_t a2 = load_tran_low_to_s16q(coeff + 512);
        const int16x8_t a3 = load_tran_low_to_s16q(coeff + 768);

        const int16x8_t b0 = vshrq_n_s16(vaddq_s16(a0, a1), 2);
        const int16x8_t b1 = vshrq_n_s16(vsubq_s16(a0, a1), 2);
        const int16x8_t b2 = vshrq_n_s16(vaddq_s16(a2, a3), 2);
        const int16x8_t b3 = vshrq_n_s16(vsubq_s16(a2, a3), 2);

        const int16x8_t c0 = vaddq_s16(b0, b2);
        const int16x8_t c1 = vaddq_s16(b1, b3);
        const int16x8_t c2 = vsubq_s16(b0, b2);
        const int16x8_t c3 = vsubq_s16(b1, b3);

        store_s16_to_tran_low(coeff + 0 + 0, vget_low_s16(c0));
        store_s16_to_tran_low(coeff + 0 + 4, vget_high_s16(c0));
        store_s16_to_tran_low(coeff + 256 + 0, vget_low_s16(c1));
        store_s16_to_tran_low(coeff + 256 + 4, vget_high_s16(c1));
        store_s16_to_tran_low(coeff + 512 + 0, vget_low_s16(c2));
        store_s16_to_tran_low(coeff + 512 + 4, vget_high_s16(c2));
        store_s16_to_tran_low(coeff + 768 + 0, vget_low_s16(c3));
        store_s16_to_tran_low(coeff + 768 + 4, vget_high_s16(c3));

        coeff += 8;
    }
}

void svt_full_distortion_kernel32_bits_neon(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff,
                                            uint32_t recon_coeff_stride, uint64_t distortion_result[DIST_CALC_TOTAL],
                                            uint32_t area_width, uint32_t area_height) {
    int64x2_t sum1      = vdupq_n_s64(0);
    int64x2_t sum2      = vdupq_n_s64(0);
    uint32_t  row_count = area_height;
    do {
        int32_t *coeff_temp       = coeff;
        int32_t *recon_coeff_temp = recon_coeff;

        uint32_t col_count = area_width / 4;
        do {
            int32x4_t x0 = vld1q_s32(coeff_temp);
            int32x4_t y0 = vld1q_s32(recon_coeff_temp);

            int32x2_t x_lo = vget_low_s32(x0);
            int32x2_t x_hi = vget_high_s32(x0);
            int32x2_t y_lo = vget_low_s32(y0);
            int32x2_t y_hi = vget_high_s32(y0);

            sum2 = vmlal_s32(sum2, x_lo, x_lo);
            sum2 = vmlal_s32(sum2, x_hi, x_hi);

            int32x2_t x_lo_sub = vsub_s32(x_lo, y_lo);
            int32x2_t x_hi_sub = vsub_s32(x_hi, y_hi);

            int64x2_t x_lo_wide = vmull_s32(x_lo_sub, x_lo_sub);
            int64x2_t x_hi_wide = vmull_s32(x_hi_sub, x_hi_sub);

            sum1 = vaddq_s64(sum1, x_lo_wide);
            sum1 = vaddq_s64(sum1, x_hi_wide);

            coeff_temp += 4;
            recon_coeff_temp += 4;
        } while (--col_count);

        coeff += coeff_stride;
        recon_coeff += recon_coeff_stride;
        row_count -= 1;
    } while (row_count > 0);

    int64x2_t tmp  = vcombine_s64(vget_high_s64(sum1), vget_low_s64(sum1));
    int64x2_t tmp2 = vaddq_s64(sum1, tmp);
    tmp            = vcombine_s64(vget_high_s64(sum2), vget_low_s64(sum2));
    int64x2_t tmp3 = vaddq_s64(sum2, tmp);
    tmp3           = vcombine_s64(vget_low_s64(tmp2), vget_low_s64(tmp3));

    vst1q_s64((int64_t *)distortion_result, tmp3);
}
