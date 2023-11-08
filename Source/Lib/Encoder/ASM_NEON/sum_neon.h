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

#ifndef AOM_AOM_DSP_ARM_SUM_NEON_H_
#define AOM_AOM_DSP_ARM_SUM_NEON_H_

#include <arm_neon.h>

static INLINE int horizontal_add_s16x8(const int16x8_t a) { return vaddlvq_s16(a); }

static INLINE int horizontal_add_s32x4(const int32x4_t a) { return vaddvq_s32(a); }

static INLINE int64_t horizontal_add_s64x2(const int64x2_t a) { return vaddvq_s64(a); }

static INLINE uint64_t horizontal_add_u64x2(const uint64x2_t a) { return vaddvq_u64(a); }

static INLINE uint64_t horizontal_long_add_u32x4(const uint32x4_t a) { return vaddlvq_u32(a); }

static INLINE int64_t horizontal_long_add_s32x4(const int32x4_t a) { return vaddlvq_s32(a); }

static INLINE unsigned int horizontal_add_u32x4(const uint32x4_t a) { return vaddvq_u32(a); }

static INLINE uint32x4_t horizontal_add_4d_u32x4(const uint32x4_t sum[4]) {
    uint32x4_t res01 = vpaddq_u32(sum[0], sum[1]);
    uint32x4_t res23 = vpaddq_u32(sum[2], sum[3]);
    return vpaddq_u32(res01, res23);
}

static INLINE int32x4_t horizontal_add_4d_s32x4(const int32x4_t sum[4]) {
    int32x4_t res01 = vpaddq_s32(sum[0], sum[1]);
    int32x4_t res23 = vpaddq_s32(sum[2], sum[3]);
    return vpaddq_s32(res01, res23);
}

static INLINE uint32_t horizontal_long_add_u16x8(const uint16x8_t vec_lo, const uint16x8_t vec_hi) {
    return vaddlvq_u16(vec_lo) + vaddlvq_u16(vec_hi);
}

static INLINE uint32x4_t horizontal_long_add_4d_u16x8(const uint16x8_t sum_lo[4], const uint16x8_t sum_hi[4]) {
    const uint32x4_t a0 = vpaddlq_u16(sum_lo[0]);
    const uint32x4_t a1 = vpaddlq_u16(sum_lo[1]);
    const uint32x4_t a2 = vpaddlq_u16(sum_lo[2]);
    const uint32x4_t a3 = vpaddlq_u16(sum_lo[3]);
    const uint32x4_t b0 = vpadalq_u16(a0, sum_hi[0]);
    const uint32x4_t b1 = vpadalq_u16(a1, sum_hi[1]);
    const uint32x4_t b2 = vpadalq_u16(a2, sum_hi[2]);
    const uint32x4_t b3 = vpadalq_u16(a3, sum_hi[3]);
    const uint32x4_t c0 = vpaddq_u32(b0, b1);
    const uint32x4_t c1 = vpaddq_u32(b2, b3);
    return vpaddq_u32(c0, c1);
}

static INLINE uint32_t horizontal_add_u16x8(const uint16x8_t a) { return vaddlvq_u16(a); }

static INLINE uint32x4_t horizontal_add_4d_u16x8(const uint16x8_t sum[4]) {
    const uint16x8_t a0 = vpaddq_u16(sum[0], sum[1]);
    const uint16x8_t a1 = vpaddq_u16(sum[2], sum[3]);
    const uint16x8_t b0 = vpaddq_u16(a0, a1);
    return vpaddlq_u16(b0);
}

static INLINE uint32_t horizontal_add_u32x2(const uint32x2_t a) { return vaddv_u32(a); }

static INLINE uint64_t horizontal_long_add_u32x2(const uint32x2_t a) { return vaddlv_u32(a); }

static INLINE uint32_t horizontal_add_u16x4(const uint16x4_t a) { return vaddlv_u16(a); }

static INLINE int32x4_t horizontal_add_2d_s32(int32x4_t a, int32x4_t b) { return vpaddq_s32(a, b); }

static INLINE int32x2_t add_pairwise_s32x4(int32x4_t a) { return vget_low_s32(vpaddq_s32(a, a)); }

#endif // AOM_AOM_DSP_ARM_SUM_NEON_H_
