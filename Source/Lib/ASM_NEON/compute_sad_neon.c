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
#include "sum_neon.h"
#include "common_dsp_rtcd.h"
#include "aom_dsp_rtcd.h"
#include "compute_sad_c.h"
#include "utility.h"
#include "mcomp.h"

static INLINE uint16x8_t sad16_neon_init(uint8x16_t src, uint8x16_t ref) {
    const uint8x16_t abs_diff = vabdq_u8(src, ref);
    return vpaddlq_u8(abs_diff);
}

static INLINE void sad16_neon(uint8x16_t src, uint8x16_t ref, uint16x8_t *const sad_sum) {
    const uint8x16_t abs_diff = vabdq_u8(src, ref);
    *sad_sum                  = vpadalq_u8(*sad_sum, abs_diff);
}

/* Find the position of the first occurrence of 'value' in the vector 'x'.
 * Returns the position (index) of the first occurrence of 'value' in the vector 'x'. */
static INLINE uint16_t findposq_u16(uint16x8_t x, uint16_t value) {
    uint16x8_t val_mask;
    uint8x16_t is_one;
    uint64_t   idx;

    val_mask = vdupq_n_u16(value);

    /* Divide each lane in half by using vreinterpret. */
    is_one = vreinterpretq_u8_u16(vceqq_u16(x, val_mask));

    /* Pack the information in the lower 64 bits of the register by considering only alternate
     * 8-bit lanes. */
    is_one = vuzp1q_u8(is_one, is_one);

    /* Get the lower 64 bits from the 128-bit register. */
    idx = vgetq_lane_u64(vreinterpretq_u64_u8(vrev64q_u8(is_one)), 0);

    /* Calculate the position as an index, dividing by 8 to account for 8-bit lanes. */
    return __builtin_clzl(idx) / 8;
}

/* Find the position of the first occurrence of 'value' in the vector 'x'.
 * Returns the position (index) of the first occurrence of 'value' in the vector 'x'. */
static INLINE uint16_t findposq_u32(uint32x4_t x, uint32_t value) {
    uint32x4_t val_mask;
    uint16x8_t is_one;
    uint64_t   idx;

    val_mask = vdupq_n_u32(value);

    /* Divide each lane in half by using vreinterpret. */
    is_one = vreinterpretq_u16_u32(vceqq_u32(x, val_mask));

    /* Pack the information in the lower 64 bits of the register by considering only alternate
     * 16-bit lanes. */
    is_one = vuzp1q_u16(is_one, is_one);

    /* Get the lower 64 bits from the 128-bit register. */
    idx = vgetq_lane_u64(vreinterpretq_u64_u16(vrev64q_u16(is_one)), 0);

    /* Calculate the position as an index, dividing by 16 to account for 16-bit lanes. */
    return __builtin_clzl(idx) / 16;
}

/* Calculate and accumulate into 'res' the absolute differences between 8 values from 'arr' and
 * 'val' */
static INLINE uint16x8_t array_abal_with_value(uint16x8_t res, const uint8_t *arr, uint8_t val) {
    uint8x8_t a, b;
    a = vld1_u8(arr);
    b = vdup_n_u8(val);
    return vabal_u8(res, a, b);
}

static INLINE void update_best_sad_u32(uint32x4_t sad4, uint64_t *best_sad, int16_t *x_search_center,
                                       int16_t *y_search_center, int16_t x_search_index, int16_t y_search_index) {
    uint64_t temp_sad;

    /* Find the minimum SAD value out of the 4 search spaces. */
    temp_sad = vminvq_u32(sad4);

    if (temp_sad < *best_sad) {
        *best_sad        = temp_sad;
        *x_search_center = (int16_t)(x_search_index + findposq_u32(sad4, temp_sad));
        *y_search_center = y_search_index;
    }
}

static INLINE void update_best_sad_u16(uint16x8_t sad8, uint64_t *best_sad, int16_t *x_search_center,
                                       int16_t *y_search_center, int16_t x_search_index, int16_t y_search_index) {
    uint64_t temp_sad;

    /* Find the minimum SAD value out of the 8 search spaces. */
    temp_sad = vminvq_u16(sad8);

    if (temp_sad < *best_sad) {
        *best_sad        = temp_sad;
        *x_search_center = (int16_t)(x_search_index + findposq_u16(sad8, temp_sad));
        *y_search_center = y_search_index;
    }
}

static INLINE void update_best_sad(uint64_t temp_sad, uint64_t *best_sad, int16_t *x_search_center,
                                   int16_t *y_search_center, int16_t x_search_index, int16_t y_search_index) {
    if (temp_sad < *best_sad) {
        *best_sad        = temp_sad;
        *x_search_center = x_search_index;
        *y_search_center = y_search_index;
    }
}

/* Return a uint16x8 vector with 'n' lanes filled with 0 and the others filled with 65535
 * The valid range for 'n' is 0 to 7 */
static INLINE uint16x8_t prepare_maskq_u16(uint16_t n) {
    uint64_t mask    = UINT64_MAX;
    mask             = mask << (8 * n);
    uint8x8_t _mask8 = vcreate_u8(mask);
    return vreinterpretq_u16_u8(vcombine_u8(vzip1_u8(_mask8, _mask8), vzip2_u8(_mask8, _mask8)));
}

/* Return a uint32x4 vector with 'n' lanes filled with 0 and the others filled with 4294967295
 * The valid range for 'n' is 0 to 4 */
static INLINE uint32x4_t prepare_maskq_u32(uint16_t n) {
    uint64_t mask     = UINT64_MAX;
    mask              = n < 4 ? (mask << (16 * n)) : 0;
    uint16x4_t _mask8 = vcreate_u16(mask);
    return vreinterpretq_u32_u16(vcombine_u16(vzip1_u16(_mask8, _mask8), vzip2_u16(_mask8, _mask8)));
}

static INLINE uint32_t sad128xh_neon(const uint8_t *src_ptr, uint32_t src_stride, const uint8_t *ref_ptr,
                                     uint32_t ref_stride, uint32_t h) {
    /* We use 8 accumulators to prevent overflow for large values of 'h', as well
    as enabling optimal UADALP instruction throughput on CPUs that have either
    2 or 4 Neon pipes. */
    uint16x8_t sum[8];
    uint8x16_t s0, s1, s2, s3, s4, s5, s6, s7;
    uint8x16_t r0, r1, r2, r3, r4, r5, r6, r7;
    uint8x16_t diff0, diff1, diff2, diff3, diff4, diff5, diff6, diff7;
    uint32x4_t sum_u32;
    uint32_t   i;

    for (i = 0; i < 8; i++) { sum[i] = vdupq_n_u16(0); }

    i = h;
    do {
        s0     = vld1q_u8(src_ptr);
        r0     = vld1q_u8(ref_ptr);
        diff0  = vabdq_u8(s0, r0);
        sum[0] = vpadalq_u8(sum[0], diff0);

        s1     = vld1q_u8(src_ptr + 16);
        r1     = vld1q_u8(ref_ptr + 16);
        diff1  = vabdq_u8(s1, r1);
        sum[1] = vpadalq_u8(sum[1], diff1);

        s2     = vld1q_u8(src_ptr + 32);
        r2     = vld1q_u8(ref_ptr + 32);
        diff2  = vabdq_u8(s2, r2);
        sum[2] = vpadalq_u8(sum[2], diff2);

        s3     = vld1q_u8(src_ptr + 48);
        r3     = vld1q_u8(ref_ptr + 48);
        diff3  = vabdq_u8(s3, r3);
        sum[3] = vpadalq_u8(sum[3], diff3);

        s4     = vld1q_u8(src_ptr + 64);
        r4     = vld1q_u8(ref_ptr + 64);
        diff4  = vabdq_u8(s4, r4);
        sum[4] = vpadalq_u8(sum[4], diff4);

        s5     = vld1q_u8(src_ptr + 80);
        r5     = vld1q_u8(ref_ptr + 80);
        diff5  = vabdq_u8(s5, r5);
        sum[5] = vpadalq_u8(sum[5], diff5);

        s6     = vld1q_u8(src_ptr + 96);
        r6     = vld1q_u8(ref_ptr + 96);
        diff6  = vabdq_u8(s6, r6);
        sum[6] = vpadalq_u8(sum[6], diff6);

        s7     = vld1q_u8(src_ptr + 112);
        r7     = vld1q_u8(ref_ptr + 112);
        diff7  = vabdq_u8(s7, r7);
        sum[7] = vpadalq_u8(sum[7], diff7);

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    } while (--i != 0);

    sum_u32 = vpaddlq_u16(sum[0]);
    sum_u32 = vpadalq_u16(sum_u32, sum[1]);
    sum_u32 = vpadalq_u16(sum_u32, sum[2]);
    sum_u32 = vpadalq_u16(sum_u32, sum[3]);
    sum_u32 = vpadalq_u16(sum_u32, sum[4]);
    sum_u32 = vpadalq_u16(sum_u32, sum[5]);
    sum_u32 = vpadalq_u16(sum_u32, sum[6]);
    sum_u32 = vpadalq_u16(sum_u32, sum[7]);

    return vaddvq_u32(sum_u32);
}

static INLINE uint32_t sad64xh_neon(const uint8_t *src_ptr, uint32_t src_stride, const uint8_t *ref_ptr,
                                    uint32_t ref_stride, uint32_t h) {
    uint16x8_t sum[4];
    uint8x16_t s0, s1, s2, s3, r0, r1, r2, r3;
    uint8x16_t diff0, diff1, diff2, diff3;
    uint32x4_t sum_u32;
    uint32_t   i;

    sum[0] = vdupq_n_u16(0);
    sum[1] = vdupq_n_u16(0);
    sum[2] = vdupq_n_u16(0);
    sum[3] = vdupq_n_u16(0);

    i = h;
    do {
        s0     = vld1q_u8(src_ptr);
        r0     = vld1q_u8(ref_ptr);
        diff0  = vabdq_u8(s0, r0);
        sum[0] = vpadalq_u8(sum[0], diff0);

        s1     = vld1q_u8(src_ptr + 16);
        r1     = vld1q_u8(ref_ptr + 16);
        diff1  = vabdq_u8(s1, r1);
        sum[1] = vpadalq_u8(sum[1], diff1);

        s2     = vld1q_u8(src_ptr + 32);
        r2     = vld1q_u8(ref_ptr + 32);
        diff2  = vabdq_u8(s2, r2);
        sum[2] = vpadalq_u8(sum[2], diff2);

        s3     = vld1q_u8(src_ptr + 48);
        r3     = vld1q_u8(ref_ptr + 48);
        diff3  = vabdq_u8(s3, r3);
        sum[3] = vpadalq_u8(sum[3], diff3);

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    } while (--i != 0);

    sum_u32 = vpaddlq_u16(sum[0]);
    sum_u32 = vpadalq_u16(sum_u32, sum[1]);
    sum_u32 = vpadalq_u16(sum_u32, sum[2]);
    sum_u32 = vpadalq_u16(sum_u32, sum[3]);

    return vaddvq_u32(sum_u32);
}

static INLINE uint32_t sad32xh_neon(const uint8_t *src_ptr, uint32_t src_stride, const uint8_t *ref_ptr,
                                    uint32_t ref_stride, uint32_t h) {
    uint16x8_t sum[2];
    uint8x16_t s0, r0, s1, r1, diff0, diff1;
    uint32_t   i;

    sum[0] = vdupq_n_u16(0);
    sum[1] = vdupq_n_u16(0);

    i = h;
    do {
        s0     = vld1q_u8(src_ptr);
        r0     = vld1q_u8(ref_ptr);
        diff0  = vabdq_u8(s0, r0);
        sum[0] = vpadalq_u8(sum[0], diff0);

        s1     = vld1q_u8(src_ptr + 16);
        r1     = vld1q_u8(ref_ptr + 16);
        diff1  = vabdq_u8(s1, r1);
        sum[1] = vpadalq_u8(sum[1], diff1);

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    } while (--i != 0);

    return vaddlvq_u16(vaddq_u16(sum[0], sum[1]));
}

static INLINE uint32_t sad16xh_neon(const uint8_t *src_ptr, uint32_t src_stride, const uint8_t *ref_ptr,
                                    uint32_t ref_stride, uint32_t h) {
    uint16x8_t sum;
    uint8x16_t s, r, diff;
    uint32_t   i;

    sum = vdupq_n_u16(0);

    i = h;
    do {
        s = vld1q_u8(src_ptr);
        r = vld1q_u8(ref_ptr);

        diff = vabdq_u8(s, r);
        sum  = vpadalq_u8(sum, diff);

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    } while (--i != 0);

    return vaddlvq_u16(sum);
}

static INLINE uint32_t sad8xh_neon(const uint8_t *src_ptr, uint32_t src_stride, const uint8_t *ref_ptr,
                                   uint32_t ref_stride, uint32_t h) {
    uint16x8_t sum;
    uint8x8_t  s, r;
    uint32_t   i;

    sum = vdupq_n_u16(0);
    i   = h;
    do {
        s = vld1_u8(src_ptr);
        r = vld1_u8(ref_ptr);

        sum = vabal_u8(sum, s, r);

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    } while (--i != 0);

    return vaddlvq_u16(sum);
}

static INLINE uint32_t sad4xh_neon(const uint8_t *src_ptr, uint32_t src_stride, const uint8_t *ref_ptr,
                                   uint32_t ref_stride, uint32_t h) {
    uint16x8_t sum;
    uint8x8_t  s, r;
    uint32_t   i;

    sum = vdupq_n_u16(0);
    i   = h / 2;
    do {
        s = load_unaligned_u8(src_ptr, src_stride);
        r = load_unaligned_u8(ref_ptr, ref_stride);

        sum = vabal_u8(sum, s, r); // add and accumulate

        src_ptr += 2 * src_stride;
        ref_ptr += 2 * ref_stride;
    } while (--i != 0);
    return vaddlvq_u16(sum);
}

static INLINE uint32_t sad24xh_neon(const uint8_t *src_ptr, uint32_t src_stride, const uint8_t *ref_ptr,
                                    uint32_t ref_stride, uint32_t h) {
    uint32_t temp_sad;
    temp_sad = sad16xh_neon(src_ptr + 0, src_stride, ref_ptr + 0, ref_stride, h);
    temp_sad += sad8xh_neon(src_ptr + 16, src_stride, ref_ptr + 16, ref_stride, h);
    return temp_sad;
}

static INLINE uint32_t sad48xh_neon(const uint8_t *src_ptr, uint32_t src_stride, const uint8_t *ref_ptr,
                                    uint32_t ref_stride, uint32_t h) {
    uint32_t temp_sad;
    temp_sad = sad32xh_neon(src_ptr + 0, src_stride, ref_ptr + 0, ref_stride, h);
    temp_sad += sad16xh_neon(src_ptr + 32, src_stride, ref_ptr + 32, ref_stride, h);
    return temp_sad;
}

static INLINE uint32_t sad40xh_neon(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride,
                                    int h) {
    uint32_t res1, res2;
    res1 = sad32xh_neon(src_ptr, src_stride, ref_ptr, ref_stride, h);
    res2 = sad8xh_neon(src_ptr + 32, src_stride, ref_ptr + 32, ref_stride, h);
    return (res1 + res2);
}

static INLINE uint32_t sad56xh_neon(const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride,
                                    int h) {
    uint32_t res1, res2;
    res1 = sad48xh_neon(src_ptr, src_stride, ref_ptr, ref_stride, h);
    res2 = sad8xh_neon(src_ptr + 48, src_stride, ref_ptr + 48, ref_stride, h);
    return (res1 + res2);
}

static INLINE uint32x4_t sadwxhx4d_large_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                              uint32_t ref_stride, uint32_t w, uint32_t h, uint32_t h_overflow) {
    uint32x4_t sum[4];
    uint16x8_t sum_lo[4], sum_hi[4];
    uint32_t   h_limit;

    sum[0] = vdupq_n_u32(0);
    sum[1] = vdupq_n_u32(0);
    sum[2] = vdupq_n_u32(0);
    sum[3] = vdupq_n_u32(0);

    h_limit = h > h_overflow ? h_overflow : h;

    uint32_t i = 0;
    do {
        sum_lo[0] = vdupq_n_u16(0);
        sum_lo[1] = vdupq_n_u16(0);
        sum_lo[2] = vdupq_n_u16(0);
        sum_lo[3] = vdupq_n_u16(0);
        sum_hi[0] = vdupq_n_u16(0);
        sum_hi[1] = vdupq_n_u16(0);
        sum_hi[2] = vdupq_n_u16(0);
        sum_hi[3] = vdupq_n_u16(0);
        do {
            const uint8_t       *loop_src       = src;
            const uint8_t       *loop_ref       = ref;
            const uint8_t *const loop_src_limit = loop_src + w;
            do {
                const uint8x16_t s0 = vld1q_u8(loop_src);
                sad16_neon(s0, vld1q_u8(loop_ref + 0), &sum_lo[0]);
                sad16_neon(s0, vld1q_u8(loop_ref + 1), &sum_lo[1]);
                sad16_neon(s0, vld1q_u8(loop_ref + 2), &sum_lo[2]);
                sad16_neon(s0, vld1q_u8(loop_ref + 3), &sum_lo[3]);

                const uint8x16_t s1 = vld1q_u8(loop_src + 16);
                sad16_neon(s1, vld1q_u8(loop_ref + 0 + 16), &sum_hi[0]);
                sad16_neon(s1, vld1q_u8(loop_ref + 1 + 16), &sum_hi[1]);
                sad16_neon(s1, vld1q_u8(loop_ref + 2 + 16), &sum_hi[2]);
                sad16_neon(s1, vld1q_u8(loop_ref + 3 + 16), &sum_hi[3]);

                loop_src += 32;
                loop_ref += 32;
            } while (loop_src < loop_src_limit);

            src += src_stride;
            ref += ref_stride;
        } while (++i < h_limit);

        sum[0] = vpadalq_u16(sum[0], sum_lo[0]);
        sum[0] = vpadalq_u16(sum[0], sum_hi[0]);
        sum[1] = vpadalq_u16(sum[1], sum_lo[1]);
        sum[1] = vpadalq_u16(sum[1], sum_hi[1]);
        sum[2] = vpadalq_u16(sum[2], sum_lo[2]);
        sum[2] = vpadalq_u16(sum[2], sum_hi[2]);
        sum[3] = vpadalq_u16(sum[3], sum_lo[3]);
        sum[3] = vpadalq_u16(sum[3], sum_hi[3]);

        h_limit += h_overflow;
    } while (i < h);

    return horizontal_add_4d_u32x4(sum);
}

static INLINE uint32x4_t sad128xhx4d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                          uint32_t ref_stride, uint32_t h) {
    return sadwxhx4d_large_neon(src, src_stride, ref, ref_stride, 128, h, 32);
}

static INLINE uint32x4_t sad64xhx4d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                         uint32_t ref_stride, uint32_t h) {
    return sadwxhx4d_large_neon(src, src_stride, ref, ref_stride, 64, h, 64);
}

static INLINE uint32x4_t sad32xhx4d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                         uint32_t ref_stride, uint32_t h) {
    uint16x8_t sum_lo[4], sum_hi[4];

    const uint8x16_t s0_init = vld1q_u8(src);
    sum_lo[0]                = sad16_neon_init(s0_init, vld1q_u8(ref + 0));
    sum_lo[1]                = sad16_neon_init(s0_init, vld1q_u8(ref + 1));
    sum_lo[2]                = sad16_neon_init(s0_init, vld1q_u8(ref + 2));
    sum_lo[3]                = sad16_neon_init(s0_init, vld1q_u8(ref + 3));

    const uint8x16_t s1_init = vld1q_u8(src + 16);
    sum_hi[0]                = sad16_neon_init(s1_init, vld1q_u8(ref + 0 + 16));
    sum_hi[1]                = sad16_neon_init(s1_init, vld1q_u8(ref + 1 + 16));
    sum_hi[2]                = sad16_neon_init(s1_init, vld1q_u8(ref + 2 + 16));
    sum_hi[3]                = sad16_neon_init(s1_init, vld1q_u8(ref + 3 + 16));

    const uint8_t *const src_limit = src + src_stride * h;

    src += src_stride;
    ref += ref_stride;

    while (src < src_limit) {
        const uint8x16_t s0 = vld1q_u8(src);
        sad16_neon(s0, vld1q_u8(ref + 0), &sum_lo[0]);
        sad16_neon(s0, vld1q_u8(ref + 1), &sum_lo[1]);
        sad16_neon(s0, vld1q_u8(ref + 2), &sum_lo[2]);
        sad16_neon(s0, vld1q_u8(ref + 3), &sum_lo[3]);

        const uint8x16_t s1 = vld1q_u8(src + 16);
        sad16_neon(s1, vld1q_u8(ref + 0 + 16), &sum_hi[0]);
        sad16_neon(s1, vld1q_u8(ref + 1 + 16), &sum_hi[1]);
        sad16_neon(s1, vld1q_u8(ref + 2 + 16), &sum_hi[2]);
        sad16_neon(s1, vld1q_u8(ref + 3 + 16), &sum_hi[3]);

        src += src_stride;
        ref += ref_stride;
    }

    return horizontal_long_add_4d_u16x8(sum_lo, sum_hi);
}

static INLINE uint32x4_t sad16xhx4d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                         uint32_t ref_stride, uint32_t h) {
    uint16x8_t sum_u16[4];
    uint32x4_t sum_u32[4];

    const uint8x16_t s_init = vld1q_u8(src);
    sum_u16[0]              = sad16_neon_init(s_init, vld1q_u8(ref + 0));
    sum_u16[1]              = sad16_neon_init(s_init, vld1q_u8(ref + 1));
    sum_u16[2]              = sad16_neon_init(s_init, vld1q_u8(ref + 2));
    sum_u16[3]              = sad16_neon_init(s_init, vld1q_u8(ref + 3));

    const uint8_t *const src_limit = src + src_stride * h;

    src += src_stride;
    ref += ref_stride;

    while (src < src_limit) {
        const uint8x16_t s = vld1q_u8(src);
        sad16_neon(s, vld1q_u8(ref + 0), &sum_u16[0]);
        sad16_neon(s, vld1q_u8(ref + 1), &sum_u16[1]);
        sad16_neon(s, vld1q_u8(ref + 2), &sum_u16[2]);
        sad16_neon(s, vld1q_u8(ref + 3), &sum_u16[3]);

        src += src_stride;
        ref += ref_stride;
    }

    sum_u32[0] = vpaddlq_u16(sum_u16[0]);
    sum_u32[1] = vpaddlq_u16(sum_u16[1]);
    sum_u32[2] = vpaddlq_u16(sum_u16[2]);
    sum_u32[3] = vpaddlq_u16(sum_u16[3]);

    return horizontal_add_4d_u32x4(sum_u32);
}

static INLINE uint32x4_t sad8xhx4d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                        uint32_t ref_stride, uint32_t h) {
    uint16x8_t sum[4];
    uint8x8_t  s;
    uint32_t   i, ref_offset;

    sum[0] = vdupq_n_u16(0);
    sum[1] = vdupq_n_u16(0);
    sum[2] = vdupq_n_u16(0);
    sum[3] = vdupq_n_u16(0);

    ref_offset = 0;
    i          = h;
    do {
        s      = vld1_u8(src);
        sum[0] = vabal_u8(sum[0], s, vld1_u8(ref + 0 + ref_offset));
        sum[1] = vabal_u8(sum[1], s, vld1_u8(ref + 1 + ref_offset));
        sum[2] = vabal_u8(sum[2], s, vld1_u8(ref + 2 + ref_offset));
        sum[3] = vabal_u8(sum[3], s, vld1_u8(ref + 3 + ref_offset));

        src += src_stride;
        ref_offset += ref_stride;
    } while (--i != 0);

    return horizontal_add_4d_u16x8(sum);
}

static INLINE uint32x4_t sad24xhx4d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                         uint32_t ref_stride, uint32_t h) {
    uint32x4_t sad4_0, sad4_1;
    sad4_0 = sad16xhx4d_neon(src + 0, src_stride, ref + 0, ref_stride, h);
    sad4_1 = sad8xhx4d_neon(src + 16, src_stride, ref + 16, ref_stride, h);
    return vaddq_u32(sad4_0, sad4_1);
}

static INLINE uint32x4_t sad48xhx4d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                         uint32_t ref_stride, uint32_t h) {
    uint32x4_t sad4_0, sad4_1;
    sad4_0 = sad32xhx4d_neon(src + 0, src_stride, ref + 0, ref_stride, h);
    sad4_1 = sad16xhx4d_neon(src + 32, src_stride, ref + 32, ref_stride, h);
    return vaddq_u32(sad4_0, sad4_1);
}

static INLINE uint16x8_t sad4xhx8d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                        uint32_t ref_stride, uint32_t h) {
    uint32_t i;
    /* Initialize 'res' to store the sum of absolute differences (SAD) of 8 search spaces. */
    uint16x8_t res;

    res = vdupq_n_u16(0);

    for (i = 0; i < h; i++) {
        /* Take one value from each search space and one value from 'src', then accumulate the
         * absolute differences in 'res'. Repeat this step 4 times to cover a 4-width block. */
        res = array_abal_with_value(res, ref + 0, src[0]);
        res = array_abal_with_value(res, ref + 1, src[1]);
        res = array_abal_with_value(res, ref + 2, src[2]);
        res = array_abal_with_value(res, ref + 3, src[3]);
        src += src_stride;
        ref += ref_stride;
    }
    return res;
}

static INLINE uint16x8_t sad6xhx8d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                        uint32_t ref_stride, uint32_t h) {
    uint32_t i;
    /* Initialize 'res' to store the sum of absolute differences (SAD) of 8 search spaces. */
    uint16x8_t res;

    res = vdupq_n_u16(0);

    for (i = 0; i < h; i++) {
        /* Take one value from each search space and one value from 'src', then accumulate the
         * absolute differences in 'res'. Repeat this step 6 times to cover a 6-width block. */
        res = array_abal_with_value(res, ref + 0, src[0]);
        res = array_abal_with_value(res, ref + 1, src[1]);
        res = array_abal_with_value(res, ref + 2, src[2]);
        res = array_abal_with_value(res, ref + 3, src[3]);
        res = array_abal_with_value(res, ref + 4, src[4]);
        res = array_abal_with_value(res, ref + 5, src[5]);
        src += src_stride;
        ref += ref_stride;
    }
    return res;
}

static INLINE void sad12xhx8d_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref, uint32_t ref_stride,
                                   uint32_t h, uint32x4_t *res) {
    /* 'sad8' will store 8d SAD for block_width = 4 */
    uint16x8_t sad8;
    sad8   = sad4xhx8d_neon(src + 0, src_stride, ref + 0, ref_stride, h);
    res[0] = sad8xhx4d_neon(src + 4, src_stride, ref + 4, ref_stride, h);
    res[1] = sad8xhx4d_neon(src + 4, src_stride, ref + 8, ref_stride, h);
    res[0] = vaddw_u16(res[0], vget_low_u16(sad8));
    res[1] = vaddw_high_u16(res[1], sad8);
}

static INLINE void svt_sad_loop_kernel4xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                               uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                               int16_t *y_search_center, uint32_t src_stride_raw,
                                               int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint16x8_t sad8;
    uint32_t   leftover      = search_area_width & 7;
    uint16x8_t leftover_mask = prepare_maskq_u16(leftover);

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 8; x_search_index += 8) {
            /* Get the SAD of 8 search spaces aligned along the width and store it in 'sad8'. */
            sad8 = sad4xhx8d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);

            /* Update 'best_sad' */
            update_best_sad_u16(sad8, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        if (leftover) {
            /* Get the SAD of 8 search spaces aligned along the width and store it in 'sad8'. */
            sad8 = sad4xhx8d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);

            /* set undesired lanes to maximum value */
            sad8 = vorrq_u16(sad8, leftover_mask);

            /* Update 'best_sad' */
            update_best_sad_u16(sad8, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel6xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                               uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                               int16_t *y_search_center, uint32_t src_stride_raw,
                                               int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint16x8_t sad8;
    uint32_t   leftover      = search_area_width & 7;
    uint16x8_t leftover_mask = prepare_maskq_u16(leftover);

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 8; x_search_index += 8) {
            /* Get the SAD of 8 search spaces aligned along the width and store it in 'sad8'. */
            sad8 = sad6xhx8d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);

            /* Update 'best_sad' */
            update_best_sad_u16(sad8, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        if (leftover) {
            /* Get the SAD of 8 search spaces aligned along the width and store it in 'sad8'. */
            sad8 = sad6xhx8d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);

            /* set undesired lanes to maximum value */
            sad8 = vorrq_u16(sad8, leftover_mask);

            /* Update 'best_sad' */
            update_best_sad_u16(sad8, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel8xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                               uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                               int16_t *y_search_center, uint32_t src_stride_raw,
                                               int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint32x4_t sad4;
    uint64_t   temp_sad;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 4; x_search_index += 4) {
            /* Get the SAD of 4 search spaces aligned along the width and store it in 'sad4'. */
            sad4 = sad8xhx4d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad_u32(sad4, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        for (; x_search_index < search_area_width; x_search_index++) {
            /* Get the SAD of 1 search spaces aligned along the width and store it in 'temp_sad'. */
            temp_sad = sad8xh_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad(temp_sad, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }
        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel12xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                                int16_t *y_search_center, uint32_t src_stride_raw,
                                                int16_t search_area_width, int16_t search_area_height) {
    int16_t x_search_index, y_search_index;

    /* To accumulate the SAD of 8 search spaces */
    uint32x4_t sad8[2];

    uint32_t   leftover = search_area_width & 7;
    uint32x4_t leftover_mask[2];
    leftover_mask[0] = prepare_maskq_u32(MIN(leftover, 4));
    leftover_mask[1] = prepare_maskq_u32(leftover - MIN(leftover, 4));
    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 8; x_search_index += 8) {
            /* Get the SAD of 8 search spaces aligned along the width and store it in 'sad8'. */
            sad12xhx8d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height, sad8);

            /* Update 'best_sad' */
            update_best_sad_u32(sad8[0], best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
            update_best_sad_u32(
                sad8[1], best_sad, x_search_center, y_search_center, x_search_index + 4, y_search_index);
        }

        if (leftover) {
            /* Get the SAD of 8 search spaces aligned along the width and store it in 'sad8'. */
            sad12xhx8d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height, sad8);

            /* set undesired lanes to maximum value */
            sad8[0] = vorrq_u32(sad8[0], leftover_mask[0]);
            sad8[1] = vorrq_u32(sad8[1], leftover_mask[1]);

            /* Update 'best_sad' */
            update_best_sad_u32(sad8[0], best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
            update_best_sad_u32(
                sad8[1], best_sad, x_search_center, y_search_center, x_search_index + 4, y_search_index);
        }

        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel16xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                uint32_t block_height, uint32_t block_width, uint64_t *best_sad,
                                                int16_t *x_search_center, int16_t *y_search_center,
                                                uint32_t src_stride_raw, uint8_t skip_search_line,
                                                int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint32x4_t sad4;
    uint64_t   temp_sad;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        /* Skip search line */
        if (block_width == 16 && block_height <= 16 && skip_search_line) {
            if ((y_search_index & 1) == 0) {
                ref += src_stride_raw;
                continue;
            }
        }

        for (x_search_index = 0; x_search_index <= search_area_width - 4; x_search_index += 4) {
            /* Get the SAD of 4 search spaces aligned along the width and store it in 'sad4'. */
            sad4 = sad16xhx4d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad_u32(sad4, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        for (; x_search_index < search_area_width; x_search_index++) {
            /* Get the SAD of 1 search spaces aligned along the width and store it in 'temp_sad'. */
            temp_sad = sad16xh_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad(temp_sad, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }
        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel24xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                                int16_t *y_search_center, uint32_t src_stride_raw,
                                                int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint32x4_t sad4;
    uint64_t   temp_sad;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 4; x_search_index += 4) {
            /* Get the SAD of 4 search spaces aligned along the width and store it in 'sad4'. */
            sad4 = sad24xhx4d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad_u32(sad4, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        for (; x_search_index < search_area_width; x_search_index++) {
            /* Get the SAD of 1 search spaces aligned along the width and store it in 'temp_sad'. */
            temp_sad = sad24xh_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad(temp_sad, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }
        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel32xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                                int16_t *y_search_center, uint32_t src_stride_raw,
                                                int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint32x4_t sad4;
    uint64_t   temp_sad;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 4; x_search_index += 4) {
            /* Get the SAD of 4 search spaces aligned along the width and store it in 'sad4'. */
            sad4 = sad32xhx4d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad_u32(sad4, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        for (; x_search_index < search_area_width; x_search_index++) {
            /* Get the SAD of 1 search spaces aligned along the width and store it in 'temp_sad'. */
            temp_sad = sad32xh_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad(temp_sad, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }
        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel48xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                                int16_t *y_search_center, uint32_t src_stride_raw,
                                                int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint32x4_t sad4;
    uint64_t   temp_sad;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 4; x_search_index += 4) {
            /* Get the SAD of 4 search spaces aligned along the width and store it in 'sad4'. */
            sad4 = sad48xhx4d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad_u32(sad4, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        for (; x_search_index < search_area_width; x_search_index++) {
            /* Get the SAD of 1 search spaces aligned along the width and store it in 'temp_sad'. */
            temp_sad = sad48xh_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad(temp_sad, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }
        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel64xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                                int16_t *y_search_center, uint32_t src_stride_raw,
                                                int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint32x4_t sad4;
    uint64_t   temp_sad;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 4; x_search_index += 4) {
            /* Get the SAD of 4 search spaces aligned along the width and store it in 'sad4'. */
            sad4 = sad64xhx4d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad_u32(sad4, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        for (; x_search_index < search_area_width; x_search_index++) {
            /* Get the SAD of 1 search spaces aligned along the width and store it in 'temp_sad'. */
            temp_sad = sad64xh_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad(temp_sad, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }
        ref += src_stride_raw;
    }
}

static INLINE void svt_sad_loop_kernel128xh_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                 uint32_t block_height, uint64_t *best_sad, int16_t *x_search_center,
                                                 int16_t *y_search_center, uint32_t src_stride_raw,
                                                 int16_t search_area_width, int16_t search_area_height) {
    int16_t    x_search_index, y_search_index;
    uint32x4_t sad4;
    uint64_t   temp_sad;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index <= search_area_width - 4; x_search_index += 4) {
            /* Get the SAD of 4 search spaces aligned along the width and store it in 'sad4'. */
            sad4 = sad128xhx4d_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad_u32(sad4, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }

        for (; x_search_index < search_area_width; x_search_index++) {
            /* Get the SAD of 1 search spaces aligned along the width and store it in 'temp_sad'. */
            temp_sad = sad128xh_neon(src, src_stride, ref + x_search_index, ref_stride, block_height);
            update_best_sad(temp_sad, best_sad, x_search_center, y_search_center, x_search_index, y_search_index);
        }
        ref += src_stride_raw;
    }
}

void svt_sad_loop_kernel_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                              uint32_t block_height, uint32_t block_width, uint64_t *best_sad, int16_t *x_search_center,
                              int16_t *y_search_center, uint32_t src_stride_raw, uint8_t skip_search_line,
                              int16_t search_area_width, int16_t search_area_height) {
    *best_sad = UINT64_MAX;
    switch (block_width) {
    case 4: {
        svt_sad_loop_kernel4xh_neon(src,
                                    src_stride,
                                    ref,
                                    ref_stride,
                                    block_height,
                                    best_sad,
                                    x_search_center,
                                    y_search_center,
                                    src_stride_raw,
                                    search_area_width,
                                    search_area_height);
        break;
    }
    case 6: {
        svt_sad_loop_kernel6xh_neon(src,
                                    src_stride,
                                    ref,
                                    ref_stride,
                                    block_height,
                                    best_sad,
                                    x_search_center,
                                    y_search_center,
                                    src_stride_raw,
                                    search_area_width,
                                    search_area_height);
        break;
    }
    case 8: {
        svt_sad_loop_kernel8xh_neon(src,
                                    src_stride,
                                    ref,
                                    ref_stride,
                                    block_height,
                                    best_sad,
                                    x_search_center,
                                    y_search_center,
                                    src_stride_raw,
                                    search_area_width,
                                    search_area_height);
        break;
    }
    case 12: {
        svt_sad_loop_kernel12xh_neon(src,
                                     src_stride,
                                     ref,
                                     ref_stride,
                                     block_height,
                                     best_sad,
                                     x_search_center,
                                     y_search_center,
                                     src_stride_raw,
                                     search_area_width,
                                     search_area_height);
        break;
    }
    case 16: {
        svt_sad_loop_kernel16xh_neon(src,
                                     src_stride,
                                     ref,
                                     ref_stride,
                                     block_height,
                                     block_width,
                                     best_sad,
                                     x_search_center,
                                     y_search_center,
                                     src_stride_raw,
                                     skip_search_line,
                                     search_area_width,
                                     search_area_height);
        break;
    }
    case 24: {
        svt_sad_loop_kernel24xh_neon(src,
                                     src_stride,
                                     ref,
                                     ref_stride,
                                     block_height,
                                     best_sad,
                                     x_search_center,
                                     y_search_center,
                                     src_stride_raw,
                                     search_area_width,
                                     search_area_height);
        break;
    }
    case 32: {
        svt_sad_loop_kernel32xh_neon(src,
                                     src_stride,
                                     ref,
                                     ref_stride,
                                     block_height,
                                     best_sad,
                                     x_search_center,
                                     y_search_center,
                                     src_stride_raw,
                                     search_area_width,
                                     search_area_height);
        break;
    }
    case 48: {
        svt_sad_loop_kernel48xh_neon(src,
                                     src_stride,
                                     ref,
                                     ref_stride,
                                     block_height,
                                     best_sad,
                                     x_search_center,
                                     y_search_center,
                                     src_stride_raw,
                                     search_area_width,
                                     search_area_height);
        break;
    }
    case 64: {
        svt_sad_loop_kernel64xh_neon(src,
                                     src_stride,
                                     ref,
                                     ref_stride,
                                     block_height,
                                     best_sad,
                                     x_search_center,
                                     y_search_center,
                                     src_stride_raw,
                                     search_area_width,
                                     search_area_height);
        break;
    }
    case 128: {
        svt_sad_loop_kernel128xh_neon(src,
                                      src_stride,
                                      ref,
                                      ref_stride,
                                      block_height,
                                      best_sad,
                                      x_search_center,
                                      y_search_center,
                                      src_stride_raw,
                                      search_area_width,
                                      search_area_height);
        break;
    }
    default: {
        svt_sad_loop_kernel_c(src,
                              src_stride,
                              ref,
                              ref_stride,
                              block_height,
                              block_width,
                              best_sad,
                              x_search_center,
                              y_search_center,
                              src_stride_raw,
                              skip_search_line,
                              search_area_width,
                              search_area_height);
        break;
    }
    }
}

static INLINE uint32x4_t get_mv_cost_vector(const struct svt_mv_cost_param *mv_cost_params, int16_t row, int16_t col,
                                            int16_t mvx, int16_t mvy, int16_t search_position_start_x,
                                            int16_t search_position_start_y) {
    const MV   baseMv = {(int16_t)(mvy + (search_position_start_y + row) * 8),
                         (int16_t)(mvx + (search_position_start_x + col) * 8)};
    const MV   mvs[4] = {{baseMv.row, baseMv.col + 8 * 0},
                         {baseMv.row, baseMv.col + 8 * 1},
                         {baseMv.row, baseMv.col + 8 * 2},
                         {baseMv.row, baseMv.col + 8 * 3}};
    uint32x4_t cost   = {(uint32_t)svt_aom_fp_mv_err_cost(&mvs[0], mv_cost_params),
                         (uint32_t)svt_aom_fp_mv_err_cost(&mvs[1], mv_cost_params),
                         (uint32_t)svt_aom_fp_mv_err_cost(&mvs[2], mv_cost_params),
                         (uint32_t)svt_aom_fp_mv_err_cost(&mvs[3], mv_cost_params)};
    return cost;
}

static INLINE void update_best_cost_u32(uint32x4_t cost4, uint32_t *best_cost, int16_t *x_search_center,
                                        int16_t *y_search_center, int16_t x_search_index, int16_t y_search_index,
                                        int16_t mvx, int16_t mvy) {
    uint64_t temp_cost;

    /* Find the minimum SAD value out of the 4 search spaces. */
    temp_cost = vminvq_u32(cost4);

    if (temp_cost < *best_cost) {
        *best_cost       = temp_cost;
        *x_search_center = mvx + (int16_t)(x_search_index + findposq_u32(cost4, temp_cost)) * 8;
        *y_search_center = mvy + y_search_index * 8;
    }
}

static INLINE void svt_pme_sad_loop_kernel4xh_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src,
                                                   uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                   uint32_t block_height, uint32_t *best_cost, int16_t *best_mvx,
                                                   int16_t *best_mvy, int16_t search_position_start_x,
                                                   int16_t search_position_start_y, int16_t search_area_width,
                                                   int16_t search_area_height, int16_t search_step, int16_t mvx,
                                                   int16_t mvy) {
    int16_t        i, j;
    const uint8_t *p_ref, *p_src;
    uint16x8_t     res8;
    uint32x4_t     cost4_low, cost4_high, res4_1, res4_2;

    for (i = 0; i < search_area_height; i += search_step) {
        for (j = 0; j <= search_area_width - 8; j += (8 + search_step - 1)) {
            p_src                = src;
            p_ref                = ref + j;
            res8                 = sad4xhx8d_neon(p_src, src_stride, p_ref, ref_stride, block_height);
            uint32x4_t res8_low  = vmovl_u16(vget_low_u16(res8));
            uint32x4_t res8_high = vmovl_u16(vget_high_u16(res8));
            cost4_low            = get_mv_cost_vector(
                mv_cost_params, i, j, mvx, mvy, search_position_start_x, search_position_start_y);
            cost4_high = get_mv_cost_vector(
                mv_cost_params, i, j + 4, mvx, mvy, search_position_start_x, search_position_start_y);
            res4_1 = vaddq_u32(res8_low, cost4_low);
            res4_2 = vaddq_u32(res8_high, cost4_high);
            update_best_cost_u32(res4_1,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
            update_best_cost_u32(res4_2,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + 4 + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
        }

        ref += search_step * ref_stride;
    }
}

static INLINE void svt_pme_sad_loop_kernel8xh_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src,
                                                   uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                   uint32_t block_height, uint32_t *best_cost, int16_t *best_mvx,
                                                   int16_t *best_mvy, int16_t search_position_start_x,
                                                   int16_t search_position_start_y, int16_t search_area_width,
                                                   int16_t search_area_height, int16_t search_step, int16_t mvx,
                                                   int16_t mvy) {
    int16_t        i, j;
    const uint8_t *p_ref, *p_src;
    uint32x4_t     sad4_low, sad4_high, cost4_low, cost4_high, res4_1, res4_2;

    for (i = 0; i < search_area_height; i += search_step) {
        for (j = 0; j <= search_area_width - 8; j += (8 + search_step - 1)) {
            p_src     = src;
            p_ref     = ref + j;
            sad4_low  = sad8xhx4d_neon(p_src, src_stride, p_ref, ref_stride, block_height);
            sad4_high = sad8xhx4d_neon(p_src, src_stride, p_ref + 4, ref_stride, block_height);

            cost4_low = get_mv_cost_vector(
                mv_cost_params, i, j, mvx, mvy, search_position_start_x, search_position_start_y);
            cost4_high = get_mv_cost_vector(
                mv_cost_params, i, j + 4, mvx, mvy, search_position_start_x, search_position_start_y);
            res4_1 = vaddq_u32(sad4_low, cost4_low);
            res4_2 = vaddq_u32(sad4_high, cost4_high);
            update_best_cost_u32(res4_1,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
            update_best_cost_u32(res4_2,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + 4 + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
        }

        ref += search_step * ref_stride;
    }
}

static INLINE void svt_pme_sad_loop_kernel16xh_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src,
                                                    uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                    uint32_t block_height, uint32_t *best_cost, int16_t *best_mvx,
                                                    int16_t *best_mvy, int16_t search_position_start_x,
                                                    int16_t search_position_start_y, int16_t search_area_width,
                                                    int16_t search_area_height, int16_t search_step, int16_t mvx,
                                                    int16_t mvy) {
    int16_t        i, j;
    const uint8_t *p_ref, *p_src;
    uint32x4_t     sad4_low, sad4_high, cost4_low, cost4_high, res4_1, res4_2;

    for (i = 0; i < search_area_height; i += search_step) {
        for (j = 0; j <= search_area_width - 8; j += (8 + search_step - 1)) {
            p_src     = src;
            p_ref     = ref + j;
            sad4_low  = sad16xhx4d_neon(p_src, src_stride, p_ref, ref_stride, block_height);
            sad4_high = sad16xhx4d_neon(p_src, src_stride, p_ref + 4, ref_stride, block_height);

            cost4_low = get_mv_cost_vector(
                mv_cost_params, i, j, mvx, mvy, search_position_start_x, search_position_start_y);
            cost4_high = get_mv_cost_vector(
                mv_cost_params, i, j + 4, mvx, mvy, search_position_start_x, search_position_start_y);
            res4_1 = vaddq_u32(sad4_low, cost4_low);
            res4_2 = vaddq_u32(sad4_high, cost4_high);
            update_best_cost_u32(res4_1,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
            update_best_cost_u32(res4_2,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + 4 + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
        }

        ref += search_step * ref_stride;
    }
}

static INLINE void svt_pme_sad_loop_kernel24xh_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src,
                                                    uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                    uint32_t block_height, uint32_t *best_cost, int16_t *best_mvx,
                                                    int16_t *best_mvy, int16_t search_position_start_x,
                                                    int16_t search_position_start_y, int16_t search_area_width,
                                                    int16_t search_area_height, int16_t search_step, int16_t mvx,
                                                    int16_t mvy) {
    int16_t        i, j;
    const uint8_t *p_ref, *p_src;
    uint32x4_t     sad4_low, sad4_high, cost4_low, cost4_high, res4_1, res4_2;

    for (i = 0; i < search_area_height; i += search_step) {
        for (j = 0; j <= search_area_width - 8; j += (8 + search_step - 1)) {
            p_src     = src;
            p_ref     = ref + j;
            sad4_low  = sad24xhx4d_neon(p_src, src_stride, p_ref, ref_stride, block_height);
            sad4_high = sad24xhx4d_neon(p_src, src_stride, p_ref + 4, ref_stride, block_height);

            cost4_low = get_mv_cost_vector(
                mv_cost_params, i, j, mvx, mvy, search_position_start_x, search_position_start_y);
            cost4_high = get_mv_cost_vector(
                mv_cost_params, i, j + 4, mvx, mvy, search_position_start_x, search_position_start_y);
            res4_1 = vaddq_u32(sad4_low, cost4_low);
            res4_2 = vaddq_u32(sad4_high, cost4_high);
            update_best_cost_u32(res4_1,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
            update_best_cost_u32(res4_2,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + 4 + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
        }

        ref += search_step * ref_stride;
    }
}

static INLINE void svt_pme_sad_loop_kernel32xh_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src,
                                                    uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                    uint32_t block_height, uint32_t *best_cost, int16_t *best_mvx,
                                                    int16_t *best_mvy, int16_t search_position_start_x,
                                                    int16_t search_position_start_y, int16_t search_area_width,
                                                    int16_t search_area_height, int16_t search_step, int16_t mvx,
                                                    int16_t mvy) {
    int16_t        i, j;
    const uint8_t *p_ref, *p_src;
    uint32x4_t     sad4_low, sad4_high, cost4_low, cost4_high, res4_1, res4_2;

    for (i = 0; i < search_area_height; i += search_step) {
        for (j = 0; j <= search_area_width - 8; j += (8 + search_step - 1)) {
            p_src     = src;
            p_ref     = ref + j;
            sad4_low  = sad32xhx4d_neon(p_src, src_stride, p_ref, ref_stride, block_height);
            sad4_high = sad32xhx4d_neon(p_src, src_stride, p_ref + 4, ref_stride, block_height);

            cost4_low = get_mv_cost_vector(
                mv_cost_params, i, j, mvx, mvy, search_position_start_x, search_position_start_y);
            cost4_high = get_mv_cost_vector(
                mv_cost_params, i, j + 4, mvx, mvy, search_position_start_x, search_position_start_y);
            res4_1 = vaddq_u32(sad4_low, cost4_low);
            res4_2 = vaddq_u32(sad4_high, cost4_high);
            update_best_cost_u32(res4_1,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
            update_best_cost_u32(res4_2,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + 4 + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
        }

        ref += search_step * ref_stride;
    }
}

static INLINE void svt_pme_sad_loop_kernel48xh_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src,
                                                    uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                    uint32_t block_height, uint32_t *best_cost, int16_t *best_mvx,
                                                    int16_t *best_mvy, int16_t search_position_start_x,
                                                    int16_t search_position_start_y, int16_t search_area_width,
                                                    int16_t search_area_height, int16_t search_step, int16_t mvx,
                                                    int16_t mvy) {
    int16_t        i, j;
    const uint8_t *p_ref, *p_src;
    uint32x4_t     sad4_low, sad4_high, cost4_low, cost4_high, res4_1, res4_2;

    for (i = 0; i < search_area_height; i += search_step) {
        for (j = 0; j <= search_area_width - 8; j += (8 + search_step - 1)) {
            p_src     = src;
            p_ref     = ref + j;
            sad4_low  = sad48xhx4d_neon(p_src, src_stride, p_ref, ref_stride, block_height);
            sad4_high = sad48xhx4d_neon(p_src, src_stride, p_ref + 4, ref_stride, block_height);

            cost4_low = get_mv_cost_vector(
                mv_cost_params, i, j, mvx, mvy, search_position_start_x, search_position_start_y);
            cost4_high = get_mv_cost_vector(
                mv_cost_params, i, j + 4, mvx, mvy, search_position_start_x, search_position_start_y);
            res4_1 = vaddq_u32(sad4_low, cost4_low);
            res4_2 = vaddq_u32(sad4_high, cost4_high);
            update_best_cost_u32(res4_1,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
            update_best_cost_u32(res4_2,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + 4 + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
        }

        ref += search_step * ref_stride;
    }
}

static INLINE void svt_pme_sad_loop_kernel64xh_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src,
                                                    uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                    uint32_t block_height, uint32_t *best_cost, int16_t *best_mvx,
                                                    int16_t *best_mvy, int16_t search_position_start_x,
                                                    int16_t search_position_start_y, int16_t search_area_width,
                                                    int16_t search_area_height, int16_t search_step, int16_t mvx,
                                                    int16_t mvy) {
    int16_t        i, j;
    const uint8_t *p_ref, *p_src;
    uint32x4_t     sad4_low, sad4_high, cost4_low, cost4_high, res4_1, res4_2;

    for (i = 0; i < search_area_height; i += search_step) {
        for (j = 0; j <= search_area_width - 8; j += (8 + search_step - 1)) {
            p_src     = src;
            p_ref     = ref + j;
            sad4_low  = sad64xhx4d_neon(p_src, src_stride, p_ref, ref_stride, block_height);
            sad4_high = sad64xhx4d_neon(p_src, src_stride, p_ref + 4, ref_stride, block_height);

            cost4_low = get_mv_cost_vector(
                mv_cost_params, i, j, mvx, mvy, search_position_start_x, search_position_start_y);
            cost4_high = get_mv_cost_vector(
                mv_cost_params, i, j + 4, mvx, mvy, search_position_start_x, search_position_start_y);

            res4_1 = vaddq_u32(sad4_low, cost4_low);
            res4_2 = vaddq_u32(sad4_high, cost4_high);
            update_best_cost_u32(res4_1,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
            update_best_cost_u32(res4_2,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + 4 + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
        }

        ref += search_step * ref_stride;
    }
}

static INLINE void svt_pme_sad_loop_kernel128xh_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src,
                                                     uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                     uint32_t block_height, uint32_t *best_cost, int16_t *best_mvx,
                                                     int16_t *best_mvy, int16_t search_position_start_x,
                                                     int16_t search_position_start_y, int16_t search_area_width,
                                                     int16_t search_area_height, int16_t search_step, int16_t mvx,
                                                     int16_t mvy) {
    int16_t        i, j;
    const uint8_t *p_ref, *p_src;
    uint32x4_t     sad4_low, sad4_high, cost4_low, cost4_high, res4_1, res4_2;

    for (i = 0; i < search_area_height; i += search_step) {
        for (j = 0; j <= search_area_width - 8; j += (8 + search_step - 1)) {
            p_src     = src;
            p_ref     = ref + j;
            sad4_low  = sad128xhx4d_neon(p_src, src_stride, p_ref, ref_stride, block_height);
            sad4_high = sad128xhx4d_neon(p_src, src_stride, p_ref + 4, ref_stride, block_height);

            cost4_low = get_mv_cost_vector(
                mv_cost_params, i, j, mvx, mvy, search_position_start_x, search_position_start_y);
            cost4_high = get_mv_cost_vector(
                mv_cost_params, i, j + 4, mvx, mvy, search_position_start_x, search_position_start_y);
            res4_1 = vaddq_u32(sad4_low, cost4_low);
            res4_2 = vaddq_u32(sad4_high, cost4_high);
            update_best_cost_u32(res4_1,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
            update_best_cost_u32(res4_2,
                                 best_cost,
                                 best_mvx,
                                 best_mvy,
                                 j + 4 + search_position_start_x,
                                 i + search_position_start_y,
                                 mvx,
                                 mvy);
        }

        ref += search_step * ref_stride;
    }
}

void svt_pme_sad_loop_kernel_neon(const struct svt_mv_cost_param *mv_cost_params, uint8_t *src, uint32_t src_stride,
                                  uint8_t *ref, uint32_t ref_stride, uint32_t block_height, uint32_t block_width,
                                  uint32_t *best_cost, int16_t *best_mvx, int16_t *best_mvy,
                                  int16_t search_position_start_x, int16_t search_position_start_y,
                                  int16_t search_area_width, int16_t search_area_height, int16_t search_step,
                                  int16_t mvx, int16_t mvy) {
    switch (block_width) {
    case 4: {
        svt_pme_sad_loop_kernel4xh_neon(mv_cost_params,
                                        src,
                                        src_stride,
                                        ref,
                                        ref_stride,
                                        block_height,
                                        best_cost,
                                        best_mvx,
                                        best_mvy,
                                        search_position_start_x,
                                        search_position_start_y,
                                        search_area_width,
                                        search_area_height,
                                        search_step,
                                        mvx,
                                        mvy);
        break;
    }
    case 8: {
        svt_pme_sad_loop_kernel8xh_neon(mv_cost_params,
                                        src,
                                        src_stride,
                                        ref,
                                        ref_stride,
                                        block_height,
                                        best_cost,
                                        best_mvx,
                                        best_mvy,
                                        search_position_start_x,
                                        search_position_start_y,
                                        search_area_width,
                                        search_area_height,
                                        search_step,
                                        mvx,
                                        mvy);
        break;
    }
    case 16: {
        svt_pme_sad_loop_kernel16xh_neon(mv_cost_params,
                                         src,
                                         src_stride,
                                         ref,
                                         ref_stride,
                                         block_height,
                                         best_cost,
                                         best_mvx,
                                         best_mvy,
                                         search_position_start_x,
                                         search_position_start_y,
                                         search_area_width,
                                         search_area_height,
                                         search_step,
                                         mvx,
                                         mvy);
        break;
    }
    case 24: {
        svt_pme_sad_loop_kernel24xh_neon(mv_cost_params,
                                         src,
                                         src_stride,
                                         ref,
                                         ref_stride,
                                         block_height,
                                         best_cost,
                                         best_mvx,
                                         best_mvy,
                                         search_position_start_x,
                                         search_position_start_y,
                                         search_area_width,
                                         search_area_height,
                                         search_step,
                                         mvx,
                                         mvy);
        break;
    }
    case 32: {
        svt_pme_sad_loop_kernel32xh_neon(mv_cost_params,
                                         src,
                                         src_stride,
                                         ref,
                                         ref_stride,
                                         block_height,
                                         best_cost,
                                         best_mvx,
                                         best_mvy,
                                         search_position_start_x,
                                         search_position_start_y,
                                         search_area_width,
                                         search_area_height,
                                         search_step,
                                         mvx,
                                         mvy);
        break;
    }
    case 48: {
        svt_pme_sad_loop_kernel48xh_neon(mv_cost_params,
                                         src,
                                         src_stride,
                                         ref,
                                         ref_stride,
                                         block_height,
                                         best_cost,
                                         best_mvx,
                                         best_mvy,
                                         search_position_start_x,
                                         search_position_start_y,
                                         search_area_width,
                                         search_area_height,
                                         search_step,
                                         mvx,
                                         mvy);
        break;
    }
    case 64: {
        svt_pme_sad_loop_kernel64xh_neon(mv_cost_params,
                                         src,
                                         src_stride,
                                         ref,
                                         ref_stride,
                                         block_height,
                                         best_cost,
                                         best_mvx,
                                         best_mvy,
                                         search_position_start_x,
                                         search_position_start_y,
                                         search_area_width,
                                         search_area_height,
                                         search_step,
                                         mvx,
                                         mvy);
        break;
    }
    case 128: {
        svt_pme_sad_loop_kernel128xh_neon(mv_cost_params,
                                          src,
                                          src_stride,
                                          ref,
                                          ref_stride,
                                          block_height,
                                          best_cost,
                                          best_mvx,
                                          best_mvy,
                                          search_position_start_x,
                                          search_position_start_y,
                                          search_area_width,
                                          search_area_height,
                                          search_step,
                                          mvx,
                                          mvy);
        break;
    }
    default: {
        svt_pme_sad_loop_kernel_c(mv_cost_params,
                                  src,
                                  src_stride,
                                  ref,
                                  ref_stride,
                                  block_height,
                                  block_width,
                                  best_cost,
                                  best_mvx,
                                  best_mvy,
                                  search_position_start_x,
                                  search_position_start_y,
                                  search_area_width,
                                  search_area_height,
                                  search_step,
                                  mvx,
                                  mvy);
        break;
    }
    }
}

uint32_t svt_nxm_sad_kernel_helper_neon(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                        uint32_t ref_stride, uint32_t height, uint32_t width) {
    uint32_t res = 0;
    switch (width) {
    case 4: {
        res = sad4xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 8: {
        res = sad8xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 16: {
        res = sad16xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 24: {
        res = sad24xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 32: {
        res = sad32xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 40: {
        res = sad40xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 48: {
        res = sad48xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 56: {
        res = sad56xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 64: {
        res = sad64xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    case 128: {
        res = sad128xh_neon(src, src_stride, ref, ref_stride, height);
        break;
    }
    default: {
        res = svt_fast_loop_nxm_sad_kernel(src, src_stride, ref, ref_stride, height, width);
        break;
    }
    }
    return res;
}

void svt_ext_eight_sad_calculation_32x32_64x64_neon(uint32_t p_sad16x16[16][8], uint32_t *p_best_sad_32x32,
                                                    uint32_t *p_best_sad_64x64, uint32_t *p_best_mv32x32,
                                                    uint32_t *p_best_mv64x64, uint32_t mv, uint32_t p_sad32x32[4][8]) {
    uint32x4_t       tmp0     = vaddq_u32(vld1q_u32(p_sad16x16[0]), vld1q_u32(p_sad16x16[1]));
    uint32x4_t       tmp1     = vaddq_u32(vld1q_u32(p_sad16x16[2]), vld1q_u32(p_sad16x16[3]));
    const uint32x4_t sad32_a1 = vaddq_u32(tmp0, tmp1);
    vst1q_u32(p_sad32x32[0], sad32_a1);

    tmp0                      = vaddq_u32(vld1q_u32((p_sad16x16[0] + 4)), vld1q_u32((p_sad16x16[1] + 4)));
    tmp1                      = vaddq_u32(vld1q_u32((p_sad16x16[2] + 4)), vld1q_u32((p_sad16x16[3] + 4)));
    const uint32x4_t sad32_a2 = vaddq_u32(tmp0, tmp1);
    vst1q_u32((p_sad32x32[0] + 4), sad32_a2);

    tmp0                      = vaddq_u32(vld1q_u32(p_sad16x16[4]), vld1q_u32(p_sad16x16[5]));
    tmp1                      = vaddq_u32(vld1q_u32(p_sad16x16[6]), vld1q_u32(p_sad16x16[7]));
    const uint32x4_t sad32_b1 = vaddq_u32(tmp0, tmp1);
    vst1q_u32(p_sad32x32[1], sad32_b1);

    tmp0                      = vaddq_u32(vld1q_u32((p_sad16x16[4] + 4)), vld1q_u32((p_sad16x16[5] + 4)));
    tmp1                      = vaddq_u32(vld1q_u32((p_sad16x16[6] + 4)), vld1q_u32((p_sad16x16[7] + 4)));
    const uint32x4_t sad32_b2 = vaddq_u32(tmp0, tmp1);
    vst1q_u32((p_sad32x32[1] + 4), sad32_b2);

    tmp0                      = vaddq_u32(vld1q_u32(p_sad16x16[8]), vld1q_u32(p_sad16x16[9]));
    tmp1                      = vaddq_u32(vld1q_u32(p_sad16x16[10]), vld1q_u32(p_sad16x16[11]));
    const uint32x4_t sad32_c1 = vaddq_u32(tmp0, tmp1);
    vst1q_u32(p_sad32x32[2], sad32_c1);

    tmp0                      = vaddq_u32(vld1q_u32((p_sad16x16[8] + 4)), vld1q_u32((p_sad16x16[9] + 4)));
    tmp1                      = vaddq_u32(vld1q_u32((p_sad16x16[10] + 4)), vld1q_u32((p_sad16x16[11] + 4)));
    const uint32x4_t sad32_c2 = vaddq_u32(tmp0, tmp1);
    vst1q_u32((p_sad32x32[2] + 4), sad32_c2);

    tmp0                      = vaddq_u32(vld1q_u32(p_sad16x16[12]), vld1q_u32(p_sad16x16[13]));
    tmp1                      = vaddq_u32(vld1q_u32(p_sad16x16[14]), vld1q_u32(p_sad16x16[15]));
    const uint32x4_t sad32_d1 = vaddq_u32(tmp0, tmp1);
    vst1q_u32(p_sad32x32[3], sad32_d1);

    tmp0                      = vaddq_u32(vld1q_u32((p_sad16x16[12] + 4)), vld1q_u32((p_sad16x16[13] + 4)));
    tmp1                      = vaddq_u32(vld1q_u32((p_sad16x16[14] + 4)), vld1q_u32((p_sad16x16[15] + 4)));
    const uint32x4_t sad32_d2 = vaddq_u32(tmp0, tmp1);
    vst1q_u32((p_sad32x32[3] + 4), sad32_d2);

    DECLARE_ALIGNED(32, uint32_t, p_sad64x64[8]);
    tmp0 = vaddq_u32(sad32_a1, sad32_b1);
    tmp1 = vaddq_u32(sad32_c1, sad32_d1);
    vst1q_u32(p_sad64x64, vaddq_u32(tmp0, tmp1));
    tmp0 = vaddq_u32(sad32_a2, sad32_b2);
    tmp1 = vaddq_u32(sad32_c2, sad32_d2);
    vst1q_u32((p_sad64x64 + 4), vaddq_u32(tmp0, tmp1));

    DECLARE_ALIGNED(32, uint32_t, computed_idx[8]);
    uint32x4_t       search_idx = vmovl_u16(vcreate_u16(0x0003000200010000));
    const uint32x4_t mv_sse     = vdupq_n_u32(mv);
    uint32x4_t       new_mv_sse = vaddq_u32(search_idx, mv_sse);
    new_mv_sse                  = vandq_u32(new_mv_sse, vdupq_n_u32(0xffff));
    vst1q_u32(computed_idx, vorrq_u32(new_mv_sse, vandq_u32(mv_sse, vdupq_n_u32(0xffff0000))));

    search_idx = vmovl_u16(vcreate_u16(0x0007000600050004));
    new_mv_sse = vaddq_u32(search_idx, mv_sse);
    new_mv_sse = vandq_u32(new_mv_sse, vdupq_n_u32(0xffff));
    vst1q_u32((computed_idx + 4), vorrq_u32(new_mv_sse, vandq_u32(mv_sse, vdupq_n_u32(0xffff0000))));

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            if (p_sad32x32[i][j] < p_best_sad_32x32[i]) {
                p_best_sad_32x32[i] = p_sad32x32[i][j];
                p_best_mv32x32[i]   = computed_idx[j];
            }
        }
    }

    for (int j = 0; j < 8; j++) {
        if (p_sad64x64[j] < p_best_sad_64x64[0]) {
            p_best_sad_64x64[0] = p_sad64x64[j];
            p_best_mv64x64[0]   = computed_idx[j];
        }
    }
}
