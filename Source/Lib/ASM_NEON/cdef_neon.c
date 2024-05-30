/*
* Copyright(c) 2024 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <arm_neon.h>

#include "definitions.h"
#include <math.h>

static INLINE uint32_t sum32(const int32x4_t src) {
    int32x4_t dst;

    dst = vpaddq_s32(src, src);
    dst = vpaddq_s32(dst, dst);

    return (uint32_t)vgetq_lane_s32(dst, 0);
}

static INLINE uint64_t sum64(const int64x2_t src) {
    const int64x2_t dst = vaddq_s64(src, vextq_s64(src, vdupq_n_s64(0), 1));
    return (uint64_t)vgetq_lane_s64(dst, 0);
}

static INLINE void sum_32_to_64(const int32x4_t src, int64x2_t *dst) {
    const int32x4_t src_l = vzip1q_s32(src, vdupq_n_s32(0));
    const int32x4_t src_h = vzip2q_s32(src, vdupq_n_s32(0));
    *dst                  = vaddq_s64(*dst, vreinterpretq_s64_s32(src_l));
    *dst                  = vaddq_s64(*dst, vreinterpretq_s64_s32(src_h));
}

static INLINE uint64_t dist_8xn_8bit_neon(const uint8_t **src, const uint8_t *dst, const int32_t dstride,
                                          const int32_t coeff_shift, uint8_t height, uint8_t subsampling_factor) {
    int16x8_t ss = vdupq_n_s16(0);
    int16x8_t dd = vdupq_n_s16(0);
    int32x4_t s2 = vdupq_n_s32(0);
    int32x4_t sd = vdupq_n_s32(0);
    int32x4_t d2 = vdupq_n_s32(0);

    for (int32_t r = 0; r < height; r += 2 * subsampling_factor) {
        const int16x8_t s_16_0 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(*src + subsampling_factor * 8)));
        const int16x8_t s_16_1 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(*src + 0 * 8)));
        const int16x8_t d_16_0 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(dst + (r + subsampling_factor) * dstride)));
        const int16x8_t d_16_1 = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(dst + r * dstride)));

        ss = vaddq_s16(ss, s_16_0);
        ss = vaddq_s16(ss, s_16_1);
        dd = vaddq_s16(dd, d_16_0);
        dd = vaddq_s16(dd, d_16_1);

        s2 = vaddq_s32(s2,
                       vpaddq_s32(vmull_s16(vget_low_s16(s_16_0), vget_low_s16(s_16_0)),
                                  vmull_s16(vget_high_s16(s_16_0), vget_high_s16(s_16_0))));
        s2 = vaddq_s32(s2,
                       vpaddq_s32(vmull_s16(vget_low_s16(s_16_1), vget_low_s16(s_16_1)),
                                  vmull_s16(vget_high_s16(s_16_1), vget_high_s16(s_16_1))));
        sd = vaddq_s32(sd,
                       vpaddq_s32(vmull_s16(vget_low_s16(s_16_0), vget_low_s16(d_16_0)),
                                  vmull_s16(vget_high_s16(s_16_0), vget_high_s16(d_16_0))));
        sd = vaddq_s32(sd,
                       vpaddq_s32(vmull_s16(vget_low_s16(s_16_1), vget_low_s16(d_16_1)),
                                  vmull_s16(vget_high_s16(s_16_1), vget_high_s16(d_16_1))));
        d2 = vaddq_s32(d2,
                       vpaddq_s32(vmull_s16(vget_low_s16(d_16_0), vget_low_s16(d_16_0)),
                                  vmull_s16(vget_high_s16(d_16_0), vget_high_s16(d_16_0))));
        d2 = vaddq_s32(d2,
                       vpaddq_s32(vmull_s16(vget_low_s16(d_16_1), vget_low_s16(d_16_1)),
                                  vmull_s16(vget_high_s16(d_16_1), vget_high_s16(d_16_1))));

        *src += 8 * 2 * subsampling_factor; // width * 2 lines per iter. * subsampling
    }

    int16x8_t ssdd    = vpaddq_s16(ss, dd);
    ssdd              = vpaddq_s16(ssdd, ssdd);
    int32x4_t ssdd_32 = vreinterpretq_s32_u32(vmovl_u16(vreinterpret_u16_s16(vget_low_s16(ssdd))));

    int32x4_t sum = vpaddq_s32(ssdd_32, ssdd_32);

    uint64_t sum_s = vgetq_lane_s32(sum, 0);
    uint64_t sum_d = vgetq_lane_s32(sum, 1);

    uint64_t sum_s2 = sum32(s2);
    uint64_t sum_d2 = sum32(d2);
    uint64_t sum_sd = sum32(sd);

    /* Compute the variance -- the calculation cannot go negative. */
    uint64_t svar = sum_s2 - ((sum_s * sum_s + 32) >> 6);
    uint64_t dvar = sum_d2 - ((sum_d * sum_d + 32) >> 6);
    return (uint64_t)floor(.5 +
                           (sum_d2 + sum_s2 - 2 * sum_sd) * .5 * (svar + dvar + (400 << 2 * coeff_shift)) /
                               (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));
}

static INLINE void mse_4xn_8bit_neon(const uint8_t **src, const uint8_t *dst, const int32_t dstride, int32x4_t *sum,
                                     uint8_t height, uint8_t subsampling_factor) {
    for (int32_t r = 0; r < height; r += 4 * subsampling_factor) {
        const uint32_t   aa = *(uint32_t *)(*src + (0 * subsampling_factor) * 4);
        const uint32_t   ab = *(uint32_t *)(*src + (1 * subsampling_factor) * 4);
        const uint32_t   ac = *(uint32_t *)(*src + (2 * subsampling_factor) * 4);
        const uint32_t   ad = *(uint32_t *)(*src + (3 * subsampling_factor) * 4);
        const uint8x16_t s  = vcombine_u8( //
            vreinterpret_u8_u32(vzip1_u32(vdup_n_u32(aa), vdup_n_u32(ab))),
            vreinterpret_u8_u32(vzip1_u32(vdup_n_u32(ac), vdup_n_u32(ad))));

        const uint32_t   ba = *(uint32_t *)(dst + (r + (0 * subsampling_factor)) * dstride);
        const uint32_t   bb = *(uint32_t *)(dst + (r + (1 * subsampling_factor)) * dstride);
        const uint32_t   bc = *(uint32_t *)(dst + (r + (2 * subsampling_factor)) * dstride);
        const uint32_t   bd = *(uint32_t *)(dst + (r + (3 * subsampling_factor)) * dstride);
        const uint8x16_t d  = vcombine_u8( //
            vreinterpret_u8_u32(vzip1_u32(vdup_n_u32(ba), vdup_n_u32(bb))),
            vreinterpret_u8_u32(vzip1_u32(vdup_n_u32(bc), vdup_n_u32(bd))));

        const int16x8_t s_16_0 = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(s)));
        const int16x8_t s_16_1 = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(s)));
        const int16x8_t d_16_0 = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(d)));
        const int16x8_t d_16_1 = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(d)));

        const int16x8_t diff_0 = vsubq_s16(d_16_0, s_16_0);
        const int16x8_t diff_1 = vsubq_s16(d_16_1, s_16_1);
        const int32x4_t mse_0  = vpaddq_s32(
            vmulq_s32(vmovl_s16(vget_low_s16(diff_0)), vmovl_s16(vget_low_s16(diff_0))),
            vmulq_s32(vmovl_s16(vget_high_s16(diff_0)), vmovl_s16(vget_high_s16(diff_0))));
        const int32x4_t mse_1 = vpaddq_s32(
            vmulq_s32(vmovl_s16(vget_low_s16(diff_1)), vmovl_s16(vget_low_s16(diff_1))),
            vmulq_s32(vmovl_s16(vget_high_s16(diff_1)), vmovl_s16(vget_high_s16(diff_1))));

        *sum = vaddq_s32(*sum, mse_0);
        *sum = vaddq_s32(*sum, mse_1);

        *src += 4 * 4 * subsampling_factor; // with * 4 rows per iter * subsampling
    }
}

static INLINE void mse_8xn_8bit_neon(const uint8_t **src, const uint8_t *dst, const int32_t dstride, int32x4_t *sum,
                                     uint8_t height, uint8_t subsampling_factor) {
    for (int32_t r = 0; r < height; r += 2 * subsampling_factor) {
        const int16x8_t s_16_0 = vreinterpretq_s16_u16(vmovl_u8((vld1_u8(*src + subsampling_factor * 8))));
        const int16x8_t s_16_1 = vreinterpretq_s16_u16(vmovl_u8((vld1_u8(*src + 0 * 8))));
        const int16x8_t d_16_0 = vreinterpretq_s16_u16(vmovl_u8((vld1_u8(dst + (r + subsampling_factor) * dstride))));
        const int16x8_t d_16_1 = vreinterpretq_s16_u16(vmovl_u8((vld1_u8(dst + r * dstride))));

        const int16x8_t diff_0 = vsubq_s16(d_16_0, s_16_0);
        const int16x8_t diff_1 = vsubq_s16(d_16_1, s_16_1);
        const int32x4_t mse_0  = vpaddq_s32(
            vmulq_s32(vmovl_s16(vget_low_s16(diff_0)), vmovl_s16(vget_low_s16(diff_0))),
            vmulq_s32(vmovl_s16(vget_high_s16(diff_0)), vmovl_s16(vget_high_s16(diff_0))));
        const int32x4_t mse_1 = vpaddq_s32(
            vmulq_s32(vmovl_s16(vget_low_s16(diff_1)), vmovl_s16(vget_low_s16(diff_1))),
            vmulq_s32(vmovl_s16(vget_high_s16(diff_1)), vmovl_s16(vget_high_s16(diff_1))));

        *sum = vaddq_s32(*sum, mse_0);
        *sum = vaddq_s32(*sum, mse_1);

        *src += 8 * 2 * subsampling_factor;
    }
}

static INLINE void mse_4x4_8bit_2x_subsampled_neon(const uint8_t **src, const uint8_t *dst, const int32_t dstride,
                                                   int32x4_t *sum) {
    const uint8x16_t s = vld1q_u8(*src);

    const uint32_t   aa = *(uint32_t *)(dst + 0 * dstride);
    const uint32_t   ab = *(uint32_t *)(*src + 1 * 4);
    const uint32_t   ac = *(uint32_t *)(dst + 2 * dstride);
    const uint32_t   ad = *(uint32_t *)(*src + 3 * 4);
    const uint8x16_t d  = vcombine_u8( // don't add r * dstride b/c add it at end of loop iterations
        vreinterpret_u8_u32(vzip1_u32(vdup_n_u32(aa), vdup_n_u32(ab))),
        vreinterpret_u8_u32(vzip1_u32(vdup_n_u32(ac), vdup_n_u32(ad))));

    const int16x8_t s_16_0 = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(s)));
    const int16x8_t s_16_1 = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(s)));
    const int16x8_t d_16_0 = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(d)));
    const int16x8_t d_16_1 = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(d)));

    const int16x8_t diff_0 = vsubq_s16(d_16_0, s_16_0);
    const int16x8_t diff_1 = vsubq_s16(d_16_1, s_16_1);
    const int32x4_t mse_0  = vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(diff_0)), vmovl_s16(vget_low_s16(diff_0))),
                                       vmulq_s32(vmovl_s16(vget_high_s16(diff_0)), vmovl_s16(vget_high_s16(diff_0))));
    const int32x4_t mse_1  = vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(diff_1)), vmovl_s16(vget_low_s16(diff_1))),
                                       vmulq_s32(vmovl_s16(vget_high_s16(diff_1)), vmovl_s16(vget_high_s16(diff_1))));

    *sum = vaddq_s32(*sum, mse_0);
    *sum = vaddq_s32(*sum, mse_1);

    *src += 16;
}

uint64_t svt_aom_compute_cdef_dist_8bit_neon(const uint8_t *dst8, int32_t dstride, const uint8_t *src8,
                                             const CdefList *dlist, int32_t cdef_count, BlockSize bsize,
                                             int32_t coeff_shift, int32_t pli, uint8_t subsampling_factor) {
    uint64_t sum;
    int32_t  bi, bx, by;

    if ((bsize == BLOCK_8X8) && (pli == 0)) {
        sum = 0;
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += dist_8xn_8bit_neon(
                &src8, dst8 + 8 * by * dstride + 8 * bx, dstride, coeff_shift, 8, subsampling_factor);
        }
    } else {
        int64x2_t mse64 = vdupq_n_s64(0);

        if (bsize == BLOCK_8X8) {
            for (bi = 0; bi < cdef_count; bi++) {
                int32x4_t mse32 = vdupq_n_s32(0);
                by              = dlist[bi].by;
                bx              = dlist[bi].bx;
                mse_8xn_8bit_neon(
                    &src8, dst8 + (8 * by + 0) * dstride + 8 * bx, dstride, &mse32, 8, subsampling_factor);
                sum_32_to_64(mse32, &mse64);
            }
        } else if (bsize == BLOCK_4X8) {
            for (bi = 0; bi < cdef_count; bi++) {
                int32x4_t mse32 = vdupq_n_s32(0);
                by              = dlist[bi].by;
                bx              = dlist[bi].bx;
                mse_4xn_8bit_neon(
                    &src8, dst8 + (8 * by + 0) * dstride + 4 * bx, dstride, &mse32, 8, subsampling_factor);
                sum_32_to_64(mse32, &mse64);
            }
        } else if (bsize == BLOCK_8X4) {
            for (bi = 0; bi < cdef_count; bi++) {
                int32x4_t mse32 = vdupq_n_s32(0);
                by              = dlist[bi].by;
                bx              = dlist[bi].bx;
                mse_8xn_8bit_neon(&src8, dst8 + 4 * by * dstride + 8 * bx, dstride, &mse32, 4, subsampling_factor);
                sum_32_to_64(mse32, &mse64);
            }
        } else {
            assert(bsize == BLOCK_4X4);
            for (bi = 0; bi < cdef_count; bi++) {
                int32x4_t mse32 = vdupq_n_s32(0);
                by              = dlist[bi].by;
                bx              = dlist[bi].bx;
                // For 4x4 blocks, all points can be computed at once.  Subsampling is done in a special function
                // to avoid accessing memory that doesn't belong to the current picture (since subsampling is implemented
                // as a multiplier to the step size).
                if (subsampling_factor == 2)
                    mse_4x4_8bit_2x_subsampled_neon(&src8, dst8 + 4 * by * dstride + 4 * bx, dstride, &mse32);
                else
                    mse_4xn_8bit_neon(&src8, dst8 + 4 * by * dstride + 4 * bx, dstride, &mse32, 4,
                                      1); // no subsampling
                sum_32_to_64(mse32, &mse64);
            }
        }

        sum = sum64(mse64);
    }
    return sum >> 2 * coeff_shift;
}

static INLINE uint64_t dist_8xn_16bit_neon(const uint16_t **src, const uint16_t *dst, const int32_t dstride,
                                           const int32_t coeff_shift, uint8_t height, uint8_t subsampling_factor) {
    int16x8_t ss = vdupq_n_s16(0);
    int16x8_t dd = vdupq_n_s16(0);
    int32x4_t s2 = vdupq_n_s32(0);
    int32x4_t sd = vdupq_n_s32(0);
    int32x4_t d2 = vdupq_n_s32(0);
    int32x4_t ssdd;
    int32x4_t sum;

    for (int32_t r = 0; r < height; r += 2 * subsampling_factor) {
        const int16x8_t s0 = vld1q_s16(
            (const int16_t *)(*src + 0 * 8)); // don't add r * dstride b/c add it at end of loop iterations
        const int16x8_t s1 = vld1q_s16(
            (const int16_t *)(*src +
                              subsampling_factor * 8)); // don't add r * dstride b/c add it at end of loop iterations

        const int16x8_t d0 = vld1q_s16((const int16_t *)(dst + r * dstride));
        const int16x8_t d1 = vld1q_s16((const int16_t *)(dst + (r + subsampling_factor) * dstride));

        ss = vaddq_s16(ss, s0);
        ss = vaddq_s16(ss, s1);
        dd = vaddq_s16(dd, d0);
        dd = vaddq_s16(dd, d1);
        s2 = vaddq_s32(s2,
                       vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(s0)), vmovl_s16(vget_low_s16(s0))),
                                  vmulq_s32(vmovl_s16(vget_high_s16(s0)), vmovl_s16(vget_high_s16(s0)))));
        s2 = vaddq_s32(s2,
                       vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(s1)), vmovl_s16(vget_low_s16(s1))),
                                  vmulq_s32(vmovl_s16(vget_high_s16(s1)), vmovl_s16(vget_high_s16(s1)))));
        sd = vaddq_s32(sd,
                       vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(s0)), vmovl_s16(vget_low_s16(d0))),
                                  vmulq_s32(vmovl_s16(vget_high_s16(s0)), vmovl_s16(vget_high_s16(d0)))));
        sd = vaddq_s32(sd,
                       vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(s1)), vmovl_s16(vget_low_s16(d1))),
                                  vmulq_s32(vmovl_s16(vget_high_s16(s1)), vmovl_s16(vget_high_s16(d1)))));
        d2 = vaddq_s32(d2,
                       vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(d0)), vmovl_s16(vget_low_s16(d0))),
                                  vmulq_s32(vmovl_s16(vget_high_s16(d0)), vmovl_s16(vget_high_s16(d0)))));
        d2 = vaddq_s32(d2,
                       vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(d1)), vmovl_s16(vget_low_s16(d1))),
                                  vmulq_s32(vmovl_s16(vget_high_s16(d1)), vmovl_s16(vget_high_s16(d1)))));

        *src += 8 * 2 * subsampling_factor;
    }

    int32x4_t ss_lo = vreinterpretq_s32_s16(vzip1q_s16(ss, vdupq_n_s16(0)));
    int32x4_t ss_hi = vreinterpretq_s32_s16(vzip2q_s16(ss, vdupq_n_s16(0)));
    int32x4_t dd_lo = vreinterpretq_s32_s16(vzip1q_s16(dd, vdupq_n_s16(0)));
    int32x4_t dd_hi = vreinterpretq_s32_s16(vzip2q_s16(dd, vdupq_n_s16(0)));
    int32x4_t ss_32 = vaddq_s32(ss_lo, ss_hi);
    int32x4_t dd_32 = vaddq_s32(dd_lo, dd_hi);
    ssdd            = vpaddq_s32(ss_32, dd_32);
    sum             = vpaddq_s32(ssdd, ssdd);

    uint64_t sum_s  = vgetq_lane_s32(sum, 0);
    uint64_t sum_d  = vgetq_lane_s32(sum, 1);
    uint64_t sum_s2 = sum32(s2);
    uint64_t sum_d2 = sum32(d2);
    uint64_t sum_sd = sum32(sd);

    /* Compute the variance -- the calculation cannot go negative. */
    uint64_t svar = sum_s2 - ((sum_s * sum_s + 32) >> 6);
    uint64_t dvar = sum_d2 - ((sum_d * sum_d + 32) >> 6);
    return (uint64_t)floor(.5 +
                           (sum_d2 + sum_s2 - 2 * sum_sd) * .5 * (svar + dvar + (400 << 2 * coeff_shift)) /
                               (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));
}

static INLINE void mse_8xn_16bit_neon(const uint16_t **src, const uint16_t *dst, const int32_t dstride, int32x4_t *sum,
                                      uint8_t height, uint8_t subsampling_factor) {
    for (int32_t r = 0; r < height; r += 2 * subsampling_factor) {
        const int16x8_t s0 = vld1q_s16(
            (const int16_t *)(*src + 0 * 8)); // don't add r * dstride b/c add it at end of loop iterations
        const int16x8_t s1 = vld1q_s16((const int16_t *)(*src + subsampling_factor * 8));
        const int16x8_t d0 = vld1q_s16((const int16_t *)(dst + r * dstride));
        const int16x8_t d1 = vld1q_s16((const int16_t *)(dst + (r + subsampling_factor) * dstride));

        const int16x8_t diff_0 = vsubq_s16(d0, s0);
        const int16x8_t diff_1 = vsubq_s16(d1, s1);
        const int32x4_t mse_0  = vpaddq_s32(
            vmulq_s32(vmovl_s16(vget_low_s16(diff_0)), vmovl_s16(vget_low_s16(diff_0))),
            vmulq_s32(vmovl_s16(vget_high_s16(diff_0)), vmovl_s16(vget_high_s16(diff_0))));
        const int32x4_t mse_1 = vpaddq_s32(
            vmulq_s32(vmovl_s16(vget_low_s16(diff_1)), vmovl_s16(vget_low_s16(diff_1))),
            vmulq_s32(vmovl_s16(vget_high_s16(diff_1)), vmovl_s16(vget_high_s16(diff_1))));
        *sum = vaddq_s32(*sum, mse_0);
        *sum = vaddq_s32(*sum, mse_1);

        *src += 8 * 2 * subsampling_factor;
    }
}

static INLINE void mse_4xn_16bit_neon(const uint16_t **src, const uint16_t *dst, const int32_t dstride, int32x4_t *sum,
                                      uint8_t height, uint8_t subsampling_factor) {
    for (int32_t r = 0; r < height; r += 4 * subsampling_factor) {
        const int64_t   s0_val[] = {*(uint64_t *)(*src + (1 * subsampling_factor) * 4), *(uint64_t *)(*src + 0 * 4)};
        const int16x8_t s0       = vreinterpretq_s16_s64(vld1q_s64(s0_val));
        const int64_t   s1_val[] = {
            *(uint64_t *)(*src + (3 * subsampling_factor) * 4),
            *(uint64_t *)(*src +
                          (2 * subsampling_factor) * 4)}; // don't add r * dstride b/c add it at end of loop iterations
        const int16x8_t s1       = vreinterpretq_s16_s64(vld1q_s64(s1_val));
        const int64_t   d0_val[] = {*(uint64_t *)(dst + (r + (1 * subsampling_factor)) * dstride),
                                    *(uint64_t *)(dst + r * dstride)};
        const int16x8_t d0       = vreinterpretq_s16_s64(vld1q_s64(d0_val));
        const int64_t   d1_val[] = {*(uint64_t *)(dst + (r + (3 * subsampling_factor)) * dstride),
                                    *(uint64_t *)(dst + (r + (2 * subsampling_factor)) * dstride)};
        const int16x8_t d1       = vreinterpretq_s16_s64(vld1q_s64(d1_val));

        const int16x8_t diff_0 = vsubq_s16(d0, s0);
        const int16x8_t diff_1 = vsubq_s16(d1, s1);
        const int32x4_t mse_0  = vpaddq_s32(
            vmulq_s32(vmovl_s16(vget_low_s16(diff_0)), vmovl_s16(vget_low_s16(diff_0))),
            vmulq_s32(vmovl_s16(vget_high_s16(diff_0)), vmovl_s16(vget_high_s16(diff_0))));
        const int32x4_t mse_1 = vpaddq_s32(
            vmulq_s32(vmovl_s16(vget_low_s16(diff_1)), vmovl_s16(vget_low_s16(diff_1))),
            vmulq_s32(vmovl_s16(vget_high_s16(diff_1)), vmovl_s16(vget_high_s16(diff_1))));
        *sum = vaddq_s32(*sum, mse_0);
        *sum = vaddq_s32(*sum, mse_1);

        *src += 4 * 4 * subsampling_factor; // with * 4 rows per iter * subsampling
    }
}

static INLINE void mse_4x4_16bit_2x_subsampled_neon(const uint16_t **src, const uint16_t *dst, const int32_t dstride,
                                                    int32x4_t *sum) {
    const int16x8_t s0 = vld1q_s16((const int16_t *)*src);
    const int16x8_t s1 = vld1q_s16((const int16_t *)(*src + 8));

    // set every line to src so distortion will be 0

    const int64_t   d0_val[] = {*(uint64_t *)(dst + 0 * dstride), *(uint64_t *)(*src + 1 * 4)};
    const int16x8_t d0       = vreinterpretq_s16_s64(vld1q_s64(d0_val));
    const int64_t   d1_val[] = {*(uint64_t *)(dst + 2 * dstride), *(uint64_t *)(*src + 3 * 4)};
    const int16x8_t d1       = vreinterpretq_s16_s64(vld1q_s64(d1_val));

    const int16x8_t diff_0 = vsubq_s16(d0, s0);
    const int16x8_t diff_1 = vsubq_s16(d1, s1);
    const int32x4_t mse_0  = vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(diff_0)), vmovl_s16(vget_low_s16(diff_0))),
                                       vmulq_s32(vmovl_s16(vget_high_s16(diff_0)), vmovl_s16(vget_high_s16(diff_0))));
    const int32x4_t mse_1  = vpaddq_s32(vmulq_s32(vmovl_s16(vget_low_s16(diff_1)), vmovl_s16(vget_low_s16(diff_1))),
                                       vmulq_s32(vmovl_s16(vget_high_s16(diff_1)), vmovl_s16(vget_high_s16(diff_1))));
    *sum                   = vaddq_s32(*sum, mse_0);
    *sum                   = vaddq_s32(*sum, mse_1);

    *src += 16;
}

uint64_t svt_aom_compute_cdef_dist_16bit_neon(const uint16_t *dst, int32_t dstride, const uint16_t *src,
                                              const CdefList *dlist, int32_t cdef_count, BlockSize bsize,
                                              int32_t coeff_shift, int32_t pli, uint8_t subsampling_factor) {
    uint64_t sum;
    int32_t  bi, bx, by;

    if ((bsize == BLOCK_8X8) && (pli == 0)) {
        sum = 0;
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += dist_8xn_16bit_neon(
                &src, dst + 8 * by * dstride + 8 * bx, dstride, coeff_shift, 8, subsampling_factor);
        }
    } else {
        int64x2_t mse64 = vdupq_n_s64(0);

        if (bsize == BLOCK_8X8) {
            for (bi = 0; bi < cdef_count; bi++) {
                int32x4_t mse32 = vdupq_n_s32(0);
                by              = dlist[bi].by;
                bx              = dlist[bi].bx;
                mse_8xn_16bit_neon(&src, dst + (8 * by + 0) * dstride + 8 * bx, dstride, &mse32, 8, subsampling_factor);
                sum_32_to_64(mse32, &mse64);
            }
        } else if (bsize == BLOCK_4X8) {
            for (bi = 0; bi < cdef_count; bi++) {
                int32x4_t mse32 = vdupq_n_s32(0);
                by              = dlist[bi].by;
                bx              = dlist[bi].bx;
                mse_4xn_16bit_neon(&src, dst + (8 * by + 0) * dstride + 4 * bx, dstride, &mse32, 8, subsampling_factor);
                sum_32_to_64(mse32, &mse64);
            }
        } else if (bsize == BLOCK_8X4) {
            for (bi = 0; bi < cdef_count; bi++) {
                int32x4_t mse32 = vdupq_n_s32(0);
                by              = dlist[bi].by;
                bx              = dlist[bi].bx;
                mse_8xn_16bit_neon(&src, dst + 4 * by * dstride + 8 * bx, dstride, &mse32, 4, subsampling_factor);
                sum_32_to_64(mse32, &mse64);
            }
        } else {
            assert(bsize == BLOCK_4X4);
            for (bi = 0; bi < cdef_count; bi++) {
                int32x4_t mse32 = vdupq_n_s32(0);
                by              = dlist[bi].by;
                bx              = dlist[bi].bx;
                // For 4x4 blocks, all points can be computed at once.  Subsampling is done in a special function
                // to avoid accessing memory that doesn't belong to the current picture (since subsampling is implemented
                // as a multiplier to the step size).
                if (subsampling_factor == 2)
                    mse_4x4_16bit_2x_subsampled_neon(&src, dst + 4 * by * dstride + 4 * bx, dstride, &mse32);
                else
                    mse_4xn_16bit_neon(&src, dst + 4 * by * dstride + 4 * bx, dstride, &mse32, 4,
                                       1); // no subsampling
                sum_32_to_64(mse32, &mse64);
            }
        }

        sum = sum64(mse64);
    }

    return sum >> 2 * coeff_shift;
}
