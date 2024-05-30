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
#include "definitions.h"
#include <arm_neon.h>

/* Here partial A is a 16-bit vector of the form: [x8 x7 x6 x5 x4 x3 x2 x1]
 * and partial B has the form:[0  y1 y2 y3 y4 y5 y6 y7].
 * This function computes (x1^2+y1^2)*C1 + (x2^2+y2^2)*C2 + ... + (x7^2+y2^7)*C7 + (x8^2+0^2)*C8
 * where the C1..C8 constants are in const1  and const2.
 */
static INLINE int32x4_t fold_mul_and_sum(int16x8_t partiala, int16x8_t partialb, int32x4_t const1, int32x4_t const2) {
    /* Reverse partial B. */
    partialb = vextq_s16(partialb, partialb, 7);
    partialb = vrev64q_s16(partialb);
    partialb = vextq_s16(partialb, partialb, 4);

    /* Interleave the x and y values of identical indices and pair x8 with 0. */
    const int16x8x2_t tmp = vzipq_s16(partialb, partiala);

    /* Square and add the corresponding x and y values. */
    const int32x4_t partiala_s32 = vpaddq_s32(vmull_s16(vget_low_s16(tmp.val[0]), vget_low_s16(tmp.val[0])),
                                              vmull_high_s16(tmp.val[0], tmp.val[0]));
    const int32x4_t partialb_s32 = vpaddq_s32(vmull_s16(vget_low_s16(tmp.val[1]), vget_low_s16(tmp.val[1])),
                                              vmull_high_s16(tmp.val[1], tmp.val[1]));

    /* Multiply by constant. */
    const int32x4_t scaled_partiala_s32 = vmulq_s32(partiala_s32, const1);
    const int32x4_t scaled_partialb_s32 = vmulq_s32(partialb_s32, const2);

    /* Sum all results. */
    return vaddq_s32(scaled_partiala_s32, scaled_partialb_s32);
}

static INLINE int32x4_t hsum4(int32x4_t x0, int32x4_t x1, int32x4_t x2, int32x4_t x3) {
    const int32x4x2_t t0 = vzipq_s32(x0, x1);
    const int32x4x2_t t1 = vzipq_s32(x2, x3);

    x0 = vcombine_s32(vget_low_s32(t0.val[0]), vget_low_s32(t1.val[0]));
    x1 = vcombine_s32(vget_high_s32(t0.val[0]), vget_high_s32(t1.val[0]));
    x2 = vcombine_s32(vget_low_s32(t0.val[1]), vget_low_s32(t1.val[1]));
    x3 = vcombine_s32(vget_high_s32(t0.val[1]), vget_high_s32(t1.val[1]));

    return vaddq_s32(vaddq_s32(x0, x1), vaddq_s32(x2, x3));
}

static INLINE void compute_directions(int16x8_t lines[8], int32_t tmp_cost1[4]) {
    /* Partial sums for lines 0 and 1. */
    int16x8_t partial4a = vextq_s16(vdupq_n_s16(0), lines[0], 8 - 7);
    int16x8_t partial4b = vextq_s16(lines[0], vdupq_n_s16(0), 1);
    partial4a           = vaddq_s16(partial4a, vextq_s16(vdupq_n_s16(0), lines[1], 8 - 6));
    partial4b           = vaddq_s16(partial4b, vextq_s16(lines[1], vdupq_n_s16(0), 2));
    int16x8_t tmp       = vaddq_s16(lines[0], lines[1]);
    int16x8_t partial5a = vextq_s16(vdupq_n_s16(0), tmp, 8 - 5);
    int16x8_t partial5b = vextq_s16(tmp, vdupq_n_s16(0), 3);
    int16x8_t partial7a = vextq_s16(vdupq_n_s16(0), tmp, 8 - 2);
    int16x8_t partial7b = vextq_s16(tmp, vdupq_n_s16(0), 6);
    int16x8_t partial6  = tmp;

    /* Partial sums for lines 2 and 3. */
    partial4a = vaddq_s16(partial4a, vextq_s16(vdupq_n_s16(0), lines[2], 8 - 5));
    partial4b = vaddq_s16(partial4b, vextq_s16(lines[2], vdupq_n_s16(0), 3));
    partial4a = vaddq_s16(partial4a, vextq_s16(vdupq_n_s16(0), lines[3], 8 - 4));
    partial4b = vaddq_s16(partial4b, vextq_s16(lines[3], vdupq_n_s16(0), 4));
    tmp       = vaddq_s16(lines[2], lines[3]);
    partial5a = vaddq_s16(partial5a, vextq_s16(vdupq_n_s16(0), tmp, 8 - 4));
    partial5b = vaddq_s16(partial5b, vextq_s16(tmp, vdupq_n_s16(0), 4));
    partial7a = vaddq_s16(partial7a, vextq_s16(vdupq_n_s16(0), tmp, 8 - 3));
    partial7b = vaddq_s16(partial7b, vextq_s16(tmp, vdupq_n_s16(0), 5));
    partial6  = vaddq_s16(partial6, tmp);

    /* Partial sums for lines 4 and 5. */
    partial4a = vaddq_s16(partial4a, vextq_s16(vdupq_n_s16(0), lines[4], 8 - 3));
    partial4b = vaddq_s16(partial4b, vextq_s16(lines[4], vdupq_n_s16(0), 5));
    partial4a = vaddq_s16(partial4a, vextq_s16(vdupq_n_s16(0), lines[5], 8 - 2));
    partial4b = vaddq_s16(partial4b, vextq_s16(lines[5], vdupq_n_s16(0), 6));
    tmp       = vaddq_s16(lines[4], lines[5]);
    partial5a = vaddq_s16(partial5a, vextq_s16(vdupq_n_s16(0), tmp, 8 - 3));
    partial5b = vaddq_s16(partial5b, vextq_s16(tmp, vdupq_n_s16(0), 5));
    partial7a = vaddq_s16(partial7a, vextq_s16(vdupq_n_s16(0), tmp, 8 - 4));
    partial7b = vaddq_s16(partial7b, vextq_s16(tmp, vdupq_n_s16(0), 4));
    partial6  = vaddq_s16(partial6, tmp);

    /* Partial sums for lines 6 and 7. */
    partial4a = vaddq_s16(partial4a, vextq_s16(vdupq_n_s16(0), lines[6], 8 - 1));
    partial4b = vaddq_s16(partial4b, vextq_s16(lines[6], vdupq_n_s16(0), 7));
    partial4a = vaddq_s16(partial4a, lines[7]);
    tmp       = vaddq_s16(lines[6], lines[7]);
    partial5a = vaddq_s16(partial5a, vextq_s16(vdupq_n_s16(0), tmp, 8 - 2));
    partial5b = vaddq_s16(partial5b, vextq_s16(tmp, vdupq_n_s16(0), 6));
    partial7a = vaddq_s16(partial7a, vextq_s16(vdupq_n_s16(0), tmp, 8 - 5));
    partial7b = vaddq_s16(partial7b, vextq_s16(tmp, vdupq_n_s16(0), 3));
    partial6  = vaddq_s16(partial6, tmp);

    /* Compute costs in terms of partial sums. */
    const int32x4_t c11           = {840, 420, 280, 210};
    const int32x4_t c12           = {168, 140, 120, 105};
    int32x4_t       partial4a_s32 = fold_mul_and_sum(partial4a, partial4b, c11, c12);

    const int32x4_t c21           = {0, 0, 420, 210};
    const int32x4_t c22           = {140, 105, 105, 105};
    int32x4_t       partial7a_s32 = fold_mul_and_sum(partial7a, partial7b, c21, c22);

    const int32x4_t c31           = {0, 0, 420, 210};
    const int32x4_t c32           = {140, 105, 105, 105};
    int32x4_t       partial5a_s32 = fold_mul_and_sum(partial5a, partial5b, c31, c32);

    int32x4_t partial6_s32 = vpaddq_s32(vmull_s16(vget_low_s16(partial6), vget_low_s16(partial6)),
                                        vmull_high_s16(partial6, partial6));

    partial6_s32 = vmulq_s32(partial6_s32, vdupq_n_s32(105));

    partial4a_s32 = hsum4(partial4a_s32, partial5a_s32, partial6_s32, partial7a_s32);
    vst1q_s32(tmp_cost1, partial4a_s32);
}

/* transpose and reverse the order of the lines -- equivalent to a 90-degree
 * counter-clockwise rotation of the pixels. */
static INLINE void array_reverse_transpose_8x8(int16x8_t *in, int16x8_t *res) {
    const int32x4_t tr0_0 = vreinterpretq_s32_s16(vzip1q_s16(in[0], in[1]));
    const int32x4_t tr0_1 = vreinterpretq_s32_s16(vzip1q_s16(in[2], in[3]));
    const int32x4_t tr0_2 = vreinterpretq_s32_s16(vzip2q_s16(in[0], in[1]));
    const int32x4_t tr0_3 = vreinterpretq_s32_s16(vzip2q_s16(in[2], in[3]));
    const int32x4_t tr0_4 = vreinterpretq_s32_s16(vzip1q_s16(in[4], in[5]));
    const int32x4_t tr0_5 = vreinterpretq_s32_s16(vzip1q_s16(in[6], in[7]));
    const int32x4_t tr0_6 = vreinterpretq_s32_s16(vzip2q_s16(in[4], in[5]));
    const int32x4_t tr0_7 = vreinterpretq_s32_s16(vzip2q_s16(in[6], in[7]));

    const int32x4_t tr1_0 = vzip1q_s32(tr0_0, tr0_1);
    const int32x4_t tr1_1 = vzip1q_s32(tr0_4, tr0_5);
    const int32x4_t tr1_2 = vzip2q_s32(tr0_0, tr0_1);
    const int32x4_t tr1_3 = vzip2q_s32(tr0_4, tr0_5);
    const int32x4_t tr1_4 = vzip1q_s32(tr0_2, tr0_3);
    const int32x4_t tr1_5 = vzip1q_s32(tr0_6, tr0_7);
    const int32x4_t tr1_6 = vzip2q_s32(tr0_2, tr0_3);
    const int32x4_t tr1_7 = vzip2q_s32(tr0_6, tr0_7);

    res[7] = vreinterpretq_s16_s32(vcombine_s32(vget_low_s32(tr1_0), vget_low_s32(tr1_1)));
    res[6] = vreinterpretq_s16_s32(vcombine_s32(vget_high_s32(tr1_0), vget_high_s32(tr1_1)));
    res[5] = vreinterpretq_s16_s32(vcombine_s32(vget_low_s32(tr1_2), vget_low_s32(tr1_3)));
    res[4] = vreinterpretq_s16_s32(vcombine_s32(vget_high_s32(tr1_2), vget_high_s32(tr1_3)));
    res[3] = vreinterpretq_s16_s32(vcombine_s32(vget_low_s32(tr1_4), vget_low_s32(tr1_5)));
    res[2] = vreinterpretq_s16_s32(vcombine_s32(vget_high_s32(tr1_4), vget_high_s32(tr1_5)));
    res[1] = vreinterpretq_s16_s32(vcombine_s32(vget_low_s32(tr1_6), vget_low_s32(tr1_7)));
    res[0] = vreinterpretq_s16_s32(vcombine_s32(vget_high_s32(tr1_6), vget_high_s32(tr1_7)));
}

void svt_aom_copy_rect8_8bit_to_16bit_neon(uint16_t *dst, int32_t dstride, const uint8_t *src, int32_t sstride,
                                           int32_t v, int32_t h) {
    int32_t i, j;
    for (i = 0; i < v; i++) {
        for (j = 0; j < (h & ~0x7); j += 8) {
            const uint8x8_t row = vld1_u8(&src[i * sstride + j]);
            vst1q_u16(&dst[i * dstride + j], vmovl_u8(row));
        }
        for (; j < h; j++) { dst[i * dstride + j] = src[i * sstride + j]; }
    }
}

uint8_t svt_aom_cdef_find_dir_neon(const uint16_t *img, int32_t stride, int32_t *var, int32_t coeff_shift) {
    int16x8_t       lines[8];
    const int16x8_t const_128 = vdupq_n_s16(128);

    for (int i = 0; i < 8; ++i) {
        int16x8_t tmp = vld1q_s16((const int16_t *)img + i * stride);
        lines[i]      = vsubq_s16(vshlq_s16(tmp, vdupq_n_s16(-(int16_t)coeff_shift)), const_128);
    }
    int32_t cost[8];

    /* Compute "mostly vertical" directions. */
    compute_directions(lines, cost + 4);

    array_reverse_transpose_8x8(lines, lines);

    /* Compute "mostly horizontal" directions. */
    compute_directions(lines, cost);

    int     best_dir  = 0;
    int32_t best_cost = 0;
    for (int i = 0; i < 8; i++) {
        if (cost[i] > best_cost) {
            best_cost = cost[i];
            best_dir  = i;
        }
    }

    /* Difference between the optimal variance and the variance along the
     orthogonal direction. Again, the sum(x^2) terms cancel out. */
    *var = best_cost - cost[(best_dir + 4) & 7];

    /* We'd normally divide by 840, but dividing by 1024 is close enough
     for what we're going to do with this. */
    *var >>= 10;
    return best_dir;
}

void svt_aom_cdef_find_dir_dual_neon(const uint16_t *img1, const uint16_t *img2, int stride, int32_t *var_out_1st,
                                     int32_t *var_out_2nd, int32_t coeff_shift, uint8_t *out_dir_1st_8x8,
                                     uint8_t *out_dir_2nd_8x8) {
    // Process first 8x8.
    *out_dir_1st_8x8 = svt_aom_cdef_find_dir_neon(img1, stride, var_out_1st, coeff_shift);

    // Process second 8x8.
    *out_dir_2nd_8x8 = svt_aom_cdef_find_dir_neon(img2, stride, var_out_2nd, coeff_shift);
}
