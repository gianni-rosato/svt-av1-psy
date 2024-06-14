/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#ifndef _HIGHBD_TXFM_UTILITY_NEON_H
#define _HIGHBD_TXFM_UTILITY_NEON_H

#include <arm_neon.h>

#define TRANSPOSE_4X4(x0, x1, x2, x3, y0, y1, y2, y3)                                                 \
    do {                                                                                              \
        int32x4_t u0, u1, u2, u3;                                                                     \
        u0 = vzip1q_s32(x0, x1);                                                                      \
        u1 = vzip2q_s32(x0, x1);                                                                      \
        u2 = vzip1q_s32(x2, x3);                                                                      \
        u3 = vzip2q_s32(x2, x3);                                                                      \
        y0 = vreinterpretq_s32_s64(vzip1q_s64(vreinterpretq_s64_s32(u0), vreinterpretq_s64_s32(u2))); \
        y1 = vreinterpretq_s32_s64(vzip2q_s64(vreinterpretq_s64_s32(u0), vreinterpretq_s64_s32(u2))); \
        y2 = vreinterpretq_s32_s64(vzip1q_s64(vreinterpretq_s64_s32(u1), vreinterpretq_s64_s32(u3))); \
        y3 = vreinterpretq_s32_s64(vzip2q_s64(vreinterpretq_s64_s32(u1), vreinterpretq_s64_s32(u3))); \
    } while (0)

static INLINE void transpose_16x16(const int32x4_t *in, int32x4_t *out) {
    // Upper left 8x8
    TRANSPOSE_4X4(in[0], in[4], in[8], in[12], out[0], out[4], out[8], out[12]);
    TRANSPOSE_4X4(in[1], in[5], in[9], in[13], out[16], out[20], out[24], out[28]);
    TRANSPOSE_4X4(in[16], in[20], in[24], in[28], out[1], out[5], out[9], out[13]);
    TRANSPOSE_4X4(in[17], in[21], in[25], in[29], out[17], out[21], out[25], out[29]);

    // Upper right 8x8
    TRANSPOSE_4X4(in[2], in[6], in[10], in[14], out[32], out[36], out[40], out[44]);
    TRANSPOSE_4X4(in[3], in[7], in[11], in[15], out[48], out[52], out[56], out[60]);
    TRANSPOSE_4X4(in[18], in[22], in[26], in[30], out[33], out[37], out[41], out[45]);
    TRANSPOSE_4X4(in[19], in[23], in[27], in[31], out[49], out[53], out[57], out[61]);

    // Lower left 8x8
    TRANSPOSE_4X4(in[32], in[36], in[40], in[44], out[2], out[6], out[10], out[14]);
    TRANSPOSE_4X4(in[33], in[37], in[41], in[45], out[18], out[22], out[26], out[30]);
    TRANSPOSE_4X4(in[48], in[52], in[56], in[60], out[3], out[7], out[11], out[15]);
    TRANSPOSE_4X4(in[49], in[53], in[57], in[61], out[19], out[23], out[27], out[31]);
    // Lower right 8x8
    TRANSPOSE_4X4(in[34], in[38], in[42], in[46], out[34], out[38], out[42], out[46]);
    TRANSPOSE_4X4(in[35], in[39], in[43], in[47], out[50], out[54], out[58], out[62]);
    TRANSPOSE_4X4(in[50], in[54], in[58], in[62], out[35], out[39], out[43], out[47]);
    TRANSPOSE_4X4(in[51], in[55], in[59], in[63], out[51], out[55], out[59], out[63]);
}

static INLINE int32x4_t half_btf_neon(const int32x4_t *w0, const int32x4_t *n0, const int32x4_t *w1,
                                      const int32x4_t *n1, int32_t bit) {
    int32x4_t x, y;
    x = vmulq_s32(*w0, *n0);
    y = vmulq_s32(*w1, *n1);
    x = vaddq_s32(x, y);
    x = vrshlq_s32(x, vdupq_n_s32(-bit));
    return x;
}

static INLINE int32x4_t half_btf_0_neon(const int32x4_t *w0, const int32x4_t *n0, int32_t bit) {
    int32x4_t x;
    x = vmulq_s32(*w0, *n0);
    x = vrshlq_s32(x, vdupq_n_s32(-bit));
    return x;
}

#endif // _HIGHBD_TXFM_UTILITY_NEON_H
