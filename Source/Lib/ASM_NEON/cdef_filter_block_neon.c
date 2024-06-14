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
#include "EbCdef.h"
#include "common_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbBitstreamUnit.h"

static INLINE int16x8_t constrain_neon(uint16x8_t a, uint16x8_t b, unsigned int threshold, int adjdamp) {
    uint16x8_t       diff   = vabdq_u16(a, b);
    const uint16x8_t a_gt_b = vcgtq_u16(a, b);
    const uint16x8_t s      = vqsubq_u16(vdupq_n_u16(threshold), vshlq_u16(diff, vdupq_n_s16(-adjdamp)));
    const int16x8_t  clip   = vreinterpretq_s16_u16(vminq_u16(diff, s));
    return vbslq_s16(a_gt_b, clip, vnegq_s16(clip));
}

void svt_av1_cdef_filter_block_8xn_8_neon(uint8_t *dst, int32_t dstride, const uint16_t *in, int32_t pri_strength,
                                          int32_t sec_strength, int32_t dir, int32_t pri_damping, int32_t sec_damping,
                                          int32_t coeff_shift, uint8_t height, uint8_t subsampling_factor) {
    int16x8_t  sum, res;
    uint16x8_t max, min, tap, row;
    uint16x8_t large = vdupq_n_u16(CDEF_VERY_LARGE);
    int16x8_t  p0, p1, p2, p3;
    uint8x8_t  ans;
    uint8_t    i;

    const int *pri_taps = svt_aom_eb_cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int *sec_taps = svt_aom_eb_cdef_sec_taps[(pri_strength >> coeff_shift) & 1];

    if (pri_strength) {
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    }
    if (sec_strength) {
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));
    }

    for (i = 0; i < height; i += subsampling_factor) {
        sum = vdupq_n_s16(0);
        row = vld1q_u16(in + i * CDEF_BSTRIDE);

        max = min = row;
        if (pri_strength) {
            /*primary near taps*/
            tap = vld1q_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir][0]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p0  = constrain_neon(tap, row, pri_strength, pri_damping);

            tap = vld1q_u16(in + (i * CDEF_BSTRIDE) - svt_aom_eb_cdef_directions[dir][0]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p1  = constrain_neon(tap, row, pri_strength, pri_damping);

            /*sum += (pri_taps[0]*(p0+p1))*/
            sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(pri_taps[0]), vaddq_s16(p0, p1)));

            /*primary far taps*/
            tap = vld1q_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir][1]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p0  = constrain_neon(tap, row, pri_strength, pri_damping);

            tap = vld1q_u16(in + (i * CDEF_BSTRIDE) - svt_aom_eb_cdef_directions[dir][1]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p1  = constrain_neon(tap, row, pri_strength, pri_damping);

            /*sum += (pri_taps[1]*(p0+p1))*/
            sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(pri_taps[1]), vaddq_s16(p0, p1)));
        }
        if (sec_strength) {
            /*secondary near taps*/
            tap = vld1q_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir + 2][0]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p0  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vld1q_u16(in + (i * CDEF_BSTRIDE) - svt_aom_eb_cdef_directions[dir + 2][0]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p1  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vld1q_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir - 2][0]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p2  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vld1q_u16(in + (i * CDEF_BSTRIDE) - svt_aom_eb_cdef_directions[dir - 2][0]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p3  = constrain_neon(tap, row, sec_strength, sec_damping);

            /*sum += (sec_taps[0]*(p0+p1+p2+p3))*/
            p0  = vaddq_s16(p0, p1);
            p2  = vaddq_s16(p2, p3);
            sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(sec_taps[0]), vaddq_s16(p0, p2)));

            /*secondary far taps*/
            tap = vld1q_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir + 2][1]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p0  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vld1q_u16(in + (i * CDEF_BSTRIDE) - svt_aom_eb_cdef_directions[dir + 2][1]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p1  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vld1q_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir - 2][1]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p2  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vld1q_u16(in + (i * CDEF_BSTRIDE) - svt_aom_eb_cdef_directions[dir - 2][1]);
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p3  = constrain_neon(tap, row, sec_strength, sec_damping);

            /*sum += (sec_taps[1]*(p0+p1+p2+p3))*/
            p0  = vaddq_s16(p0, p1);
            p2  = vaddq_s16(p2, p3);
            sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(sec_taps[1]), vaddq_s16(p0, p2)));
        }

        /*res =min(max(row+(sum+8+(sum<0))>>4 ,min),max)*/
        sum = vaddq_s16(sum, vreinterpretq_s16_u16(vcltq_s16(sum, vdupq_n_s16(0))));
        res = vaddq_s16(sum, vdupq_n_s16(8));
        res = vshrq_n_s16(res, 4);
        res = vaddq_s16(vreinterpretq_s16_u16(row), res);
        res = vminq_s16(vmaxq_s16(res, vreinterpretq_s16_u16(min)), vreinterpretq_s16_u16(max));

        ans = vqmovun_s16(res);
        vst1_u8(dst + (i * dstride), ans);
    }
}

void svt_av1_cdef_filter_block_4xn_8_neon(uint8_t *dst, int32_t dstride, const uint16_t *in, int32_t pri_strength,
                                          int32_t sec_strength, int32_t dir, int32_t pri_damping, int32_t sec_damping,
                                          int32_t coeff_shift, uint8_t height, uint8_t subsampling_factor) {
    int16x8_t  sum, res;
    uint16x8_t max, min, tap, row, large = vdupq_n_u16(CDEF_VERY_LARGE);
    int16x8_t  p0, p1, p2, p3;
    uint8x8_t  ans;
    uint8_t    i;
    uint32_t  *dst_r1_u32;
    uint32_t  *dst_r2_u32;

    const int *pri_taps = svt_aom_eb_cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int *sec_taps = svt_aom_eb_cdef_sec_taps[(pri_strength >> coeff_shift) & 1];

    if (pri_strength) {
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    }
    if (sec_strength) {
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));
    }

    for (i = 0; i < height; i += (2 * subsampling_factor)) {
        sum = vdupq_n_s16(0);
        row = vcombine_u16(vld1_u16(in + i * CDEF_BSTRIDE), vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE));
        max = min = row;
        if (pri_strength) {
            /*primary near taps*/
            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir][0]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir][0]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p0  = constrain_neon(tap, row, pri_strength, pri_damping);

            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir][0]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir][0]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p1  = constrain_neon(tap, row, pri_strength, pri_damping);

            /*sum += (pri_taps[0]*(p0+p1))*/
            sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(pri_taps[0]), vaddq_s16(p0, p1)));

            /*primary far taps*/
            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir][1]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir][1]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p0  = constrain_neon(tap, row, pri_strength, pri_damping);

            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir][1]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir][1]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p1  = constrain_neon(tap, row, pri_strength, pri_damping);

            /*sum += (pri_taps[1]*(p0+p1))*/
            sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(pri_taps[1]), vaddq_s16(p0, p1)));
        }

        if (sec_strength) {
            /*secondary near taps*/
            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir + 2][0]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir + 2][0]));

            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p0  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir + 2][0]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir + 2][0]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p1  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir - 2][0]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir - 2][0]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p2  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir - 2][0]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir - 2][0]));

            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p3  = constrain_neon(tap, row, sec_strength, sec_damping);

            /*sum += (sec_taps[0]*(p0+p1+p2+p3))*/
            p0  = vaddq_s16(p0, p1);
            p2  = vaddq_s16(p2, p3);
            sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(sec_taps[0]), vaddq_s16(p0, p2)));

            /*secondary far taps*/
            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir + 2][1]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir + 2][1]));

            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p0  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir + 2][1]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir + 2][1]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p1  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir - 2][1]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE + svt_aom_eb_cdef_directions[dir - 2][1]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p2  = constrain_neon(tap, row, sec_strength, sec_damping);

            tap = vcombine_u16(
                vld1_u16(in + i * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir - 2][1]),
                vld1_u16(in + (i + subsampling_factor) * CDEF_BSTRIDE - svt_aom_eb_cdef_directions[dir - 2][1]));
            max = vmaxq_u16(max, vandq_u16(tap, vmvnq_u16(vceqq_u16(tap, large))));
            min = vminq_u16(min, tap);
            p3  = constrain_neon(tap, row, sec_strength, sec_damping);

            /*sum += (sec_taps[1]*(p0+p1+p2+p3))*/
            p0  = vaddq_s16(p0, p1);
            p2  = vaddq_s16(p2, p3);
            sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(sec_taps[1]), vaddq_s16(p0, p2)));
        }
        /*res =min(max(row+(sum+8+(sum<0))>>4 ,min),max)*/
        sum = vaddq_s16(sum, vreinterpretq_s16_u16(vcltq_s16(sum, vdupq_n_s16(0))));
        res = vaddq_s16(sum, vdupq_n_s16(8));
        res = vshrq_n_s16(res, 4);
        res = vaddq_s16(vreinterpretq_s16_u16(row), res);
        res = vminq_s16(vmaxq_s16(res, vreinterpretq_s16_u16(min)), vreinterpretq_s16_u16(max));

        ans = vqmovun_s16(res);
        /*storing 32 bits in the destination buffer of type uint8_t*/
        dst_r1_u32  = (uint32_t *)(dst + (i * dstride));
        dst_r2_u32  = (uint32_t *)(dst + ((i + subsampling_factor) * dstride));
        *dst_r1_u32 = vget_lane_u32(vreinterpret_u32_u8(ans), 0);
        *dst_r2_u32 = vget_lane_u32(vreinterpret_u32_u8(ans), 1);
    }
}

SIMD_INLINE int16x8_t constrain16(int16x8_t a, int16x8_t b, unsigned int threshold, unsigned int adjdamp) {
    int16x8_t       diff = vreinterpretq_s16_u16(vsubq_u16(vreinterpretq_u16_s16(a), vreinterpretq_u16_s16(b)));
    const int16x8_t sign = vshrq_n_s16(diff, 15);
    diff                 = vabsq_s16(diff);
    const int16x8_t s    = vreinterpretq_s16_u16(
        vqsubq_u16(vdupq_n_u16(threshold), vreinterpretq_u16_s16(vshlq_s16(diff, vdupq_n_s16(-adjdamp)))));
    return veorq_s16(vaddq_s16(sign, vminq_s16(diff, s)), sign);
}

void svt_av1_cdef_filter_block_8xn_16_neon(uint16_t *dst, int dstride, const uint16_t *in, int pri_strength,
                                           int sec_strength, int dir, int pri_damping, int sec_damping, int coeff_shift,
                                           uint8_t height, uint8_t subsampling_factor) {
    int           i;
    int16x8_t     sum, p0, p1, p2, p3, row, res;
    int16x8_t     max, min, large = vdupq_n_s16(CDEF_VERY_LARGE);
    const int32_t po1  = svt_aom_eb_cdef_directions[dir][0];
    const int32_t po2  = svt_aom_eb_cdef_directions[dir][1];
    const int32_t s1o1 = svt_aom_eb_cdef_directions[(dir + 2)][0];
    const int32_t s1o2 = svt_aom_eb_cdef_directions[(dir + 2)][1];
    const int32_t s2o1 = svt_aom_eb_cdef_directions[(dir - 2)][0];
    const int32_t s2o2 = svt_aom_eb_cdef_directions[(dir - 2)][1];

    const int *pri_taps = svt_aom_eb_cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int *sec_taps = svt_aom_eb_cdef_sec_taps[(pri_strength >> coeff_shift) & 1];

    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));

    for (i = 0; i < height; i += subsampling_factor) {
        sum = vdupq_n_s16(0);

        const int16_t *ina = (const int16_t *)(in + i * CDEF_BSTRIDE);

        row = vld1q_s16(ina);

        min = max = row;
        // Primary near taps
        p0  = vld1q_s16(&ina[po1]);
        p1  = vld1q_s16(&ina[-po1]);
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p0, vreinterpretq_s16_u16(vceqq_s16(p0, large)))),
                        vbicq_s16(p1, vreinterpretq_s16_u16(vceqq_s16(p1, large))));
        min = vminq_s16(vminq_s16(min, p0), p1);
        p0  = constrain16(p0, row, pri_strength, pri_damping);
        p1  = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(pri_taps[0]), vaddq_s16(p0, p1)));

        // Primary far taps
        p0  = vld1q_s16(&ina[po2]);
        p1  = vld1q_s16(&ina[-po2]);
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p0, vreinterpretq_s16_u16(vceqq_s16(p0, large)))),
                        vbicq_s16(p1, vreinterpretq_s16_u16(vceqq_s16(p1, large))));
        min = vminq_s16(vminq_s16(min, p0), p1);
        p0  = constrain16(p0, row, pri_strength, pri_damping);
        p1  = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[1] * (p0 + p1)
        sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(pri_taps[1]), vaddq_s16(p0, p1)));

        // Secondary near taps
        p0  = vld1q_s16(&ina[s1o1]);
        p1  = vld1q_s16(&ina[-s1o1]);
        p2  = vld1q_s16(&ina[s2o1]);
        p3  = vld1q_s16(&ina[-s2o1]);
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p0, vreinterpretq_s16_u16(vceqq_s16(p0, large)))),
                        vbicq_s16(p1, vreinterpretq_s16_u16(vceqq_s16(p1, large))));
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p2, vreinterpretq_s16_u16(vceqq_s16(p2, large)))),
                        vbicq_s16(p3, vreinterpretq_s16_u16(vceqq_s16(p3, large))));
        min = vminq_s16(vminq_s16(vminq_s16(vminq_s16(min, p0), p1), p2), p3);
        p0  = constrain16(p0, row, sec_strength, sec_damping);
        p1  = constrain16(p1, row, sec_strength, sec_damping);
        p2  = constrain16(p2, row, sec_strength, sec_damping);
        p3  = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(sec_taps[0]), vaddq_s16(vaddq_s16(p0, p1), vaddq_s16(p2, p3))));

        // Secondary far taps
        p0  = vld1q_s16(&ina[s1o2]);
        p1  = vld1q_s16(&ina[-s1o2]);
        p2  = vld1q_s16(&ina[s2o2]);
        p3  = vld1q_s16(&ina[-s2o2]);
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p0, vreinterpretq_s16_u16(vceqq_s16(p0, large)))),
                        vbicq_s16(p1, vreinterpretq_s16_u16(vceqq_s16(p1, large))));
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p2, vreinterpretq_s16_u16(vceqq_s16(p2, large)))),
                        vbicq_s16(p3, vreinterpretq_s16_u16(vceqq_s16(p3, large))));
        min = vminq_s16(vminq_s16(vminq_s16(vminq_s16(min, p0), p1), p2), p3);
        p0  = constrain16(p0, row, sec_strength, sec_damping);
        p1  = constrain16(p1, row, sec_strength, sec_damping);
        p2  = constrain16(p2, row, sec_strength, sec_damping);
        p3  = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(sec_taps[1]), vaddq_s16(vaddq_s16(p0, p1), vaddq_s16(p2, p3))));

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = vaddq_s16(sum, vreinterpretq_s16_u16(vcltq_s16(sum, vdupq_n_s16(0))));
        res = vaddq_s16(sum, vdupq_n_s16(8));
        res = vshrq_n_s16(res, 4);
        res = vaddq_s16(row, res);
        res = vminq_s16(vmaxq_s16(res, min), max);

        vst1q_s16((int16_t *)(dst + i * dstride), res);
    }
}

void svt_av1_cdef_filter_block_4xn_16_neon(uint16_t *dst, int dstride, const uint16_t *in, int pri_strength,
                                           int sec_strength, int dir, int pri_damping, int sec_damping, int coeff_shift,
                                           uint8_t height, uint8_t subsampling_factor) {
    int           i;
    int16x8_t     p0, p1, p2, p3, sum, row, res;
    int16x8_t     max, min, large = vdupq_n_s16(CDEF_VERY_LARGE);
    const int32_t po1  = svt_aom_eb_cdef_directions[dir][0];
    const int32_t po2  = svt_aom_eb_cdef_directions[dir][1];
    const int32_t s1o1 = svt_aom_eb_cdef_directions[(dir + 2)][0];
    const int32_t s1o2 = svt_aom_eb_cdef_directions[(dir + 2)][1];
    const int32_t s2o1 = svt_aom_eb_cdef_directions[(dir - 2)][0];
    const int32_t s2o2 = svt_aom_eb_cdef_directions[(dir - 2)][1];

    const int *pri_taps = svt_aom_eb_cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int *sec_taps = svt_aom_eb_cdef_sec_taps[(pri_strength >> coeff_shift) & 1];

    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));

    for (i = 0; i < height; i += (2 * subsampling_factor)) {
        sum = vdupq_n_s16(0);

        const int16_t *ina = (const int16_t *)(in + i * CDEF_BSTRIDE);
        const int16_t *inb = (const int16_t *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE);

        row = vcombine_s16(vld1_s16(ina), vld1_s16(inb));
        min = max = row;

        // Primary near taps
        p0  = vcombine_s16(vld1_s16(&ina[po1]), vld1_s16(&inb[po1]));
        p1  = vcombine_s16(vld1_s16(&ina[-po1]), vld1_s16(&inb[-po1]));
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p0, vreinterpretq_s16_u16(vceqq_s16(p0, large)))),
                        vbicq_s16(p1, vreinterpretq_s16_u16(vceqq_s16(p1, large))));
        min = vminq_s16(vminq_s16(min, p0), p1);
        p0  = constrain16(p0, row, pri_strength, pri_damping);
        p1  = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(pri_taps[0]), vaddq_s16(p0, p1)));

        // Primary far taps
        p0  = vcombine_s16(vld1_s16(&ina[po2]), vld1_s16(&inb[po2]));
        p1  = vcombine_s16(vld1_s16(&ina[-po2]), vld1_s16(&inb[-po2]));
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p0, vreinterpretq_s16_u16(vceqq_s16(p0, large)))),
                        vbicq_s16(p1, vreinterpretq_s16_u16(vceqq_s16(p1, large))));
        min = vminq_s16(vminq_s16(min, p0), p1);
        p0  = constrain16(p0, row, pri_strength, pri_damping);
        p1  = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[1] * (p0 + p1)
        sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(pri_taps[1]), vaddq_s16(p0, p1)));

        // Secondary near taps
        p0  = vcombine_s16(vld1_s16(&ina[s1o1]), vld1_s16(&inb[s1o1]));
        p1  = vcombine_s16(vld1_s16(&ina[-s1o1]), vld1_s16(&inb[-s1o1]));
        p2  = vcombine_s16(vld1_s16(&ina[s2o1]), vld1_s16(&inb[s2o1]));
        p3  = vcombine_s16(vld1_s16(&ina[-s2o1]), vld1_s16(&inb[-s2o1]));
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p0, vreinterpretq_s16_u16(vceqq_s16(p0, large)))),
                        vbicq_s16(p1, vreinterpretq_s16_u16(vceqq_s16(p1, large))));
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p2, vreinterpretq_s16_u16(vceqq_s16(p2, large)))),
                        vbicq_s16(p3, vreinterpretq_s16_u16(vceqq_s16(p3, large))));
        min = vminq_s16(vminq_s16(vminq_s16(vminq_s16(min, p0), p1), p2), p3);
        p0  = constrain16(p0, row, sec_strength, sec_damping);
        p1  = constrain16(p1, row, sec_strength, sec_damping);
        p2  = constrain16(p2, row, sec_strength, sec_damping);
        p3  = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(sec_taps[0]), vaddq_s16(vaddq_s16(p0, p1), vaddq_s16(p2, p3))));

        // Secondary far taps
        p0  = vcombine_s16(vld1_s16(&ina[s1o2]), vld1_s16(&inb[s1o2]));
        p1  = vcombine_s16(vld1_s16(&ina[-s1o2]), vld1_s16(&inb[-s1o2]));
        p2  = vcombine_s16(vld1_s16(&ina[s2o2]), vld1_s16(&inb[s2o2]));
        p3  = vcombine_s16(vld1_s16(&ina[-s2o2]), vld1_s16(&inb[-s2o2]));
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p0, vreinterpretq_s16_u16(vceqq_s16(p0, large)))),
                        vbicq_s16(p1, vreinterpretq_s16_u16(vceqq_s16(p1, large))));
        max = vmaxq_s16(vmaxq_s16(max, vbicq_s16(p2, vreinterpretq_s16_u16(vceqq_s16(p2, large)))),
                        vbicq_s16(p3, vreinterpretq_s16_u16(vceqq_s16(p3, large))));
        min = vminq_s16(vminq_s16(vminq_s16(vminq_s16(min, p0), p1), p2), p3);
        p0  = constrain16(p0, row, sec_strength, sec_damping);
        p1  = constrain16(p1, row, sec_strength, sec_damping);
        p2  = constrain16(p2, row, sec_strength, sec_damping);
        p3  = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        sum = vaddq_s16(sum, vmulq_s16(vdupq_n_s16(sec_taps[1]), vaddq_s16(vaddq_s16(p0, p1), vaddq_s16(p2, p3))));

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = vaddq_s16(sum, vreinterpretq_s16_u16(vcltq_s16(sum, vdupq_n_s16(0))));
        res = vaddq_s16(sum, vdupq_n_s16(8));
        res = vshrq_n_s16(res, 4);
        res = vaddq_s16(row, res);
        res = vminq_s16(vmaxq_s16(res, min), max);

        vst1_s16((int16_t *)(dst + i * dstride), vget_low_s16(res));
        vst1_s16((int16_t *)(dst + (i + 1 * subsampling_factor) * dstride), vget_high_s16(res));
    }
}

void svt_cdef_filter_block_neon(uint8_t *dst8, uint16_t *dst16, int32_t dstride, const uint16_t *in,
                                int32_t pri_strength, int32_t sec_strength, int32_t dir, int32_t pri_damping,
                                int32_t sec_damping, int32_t bsize, int32_t coeff_shift, uint8_t subsampling_factor) {
    if (dst8) {
        if (bsize == BLOCK_8X8) {
            svt_av1_cdef_filter_block_8xn_8_neon(dst8,
                                                 dstride,
                                                 in,
                                                 pri_strength,
                                                 sec_strength,
                                                 dir,
                                                 pri_damping,
                                                 sec_damping,
                                                 coeff_shift,
                                                 8,
                                                 subsampling_factor);
        } else if (bsize == BLOCK_4X8) {
            svt_av1_cdef_filter_block_4xn_8_neon(dst8,
                                                 dstride,
                                                 in,
                                                 pri_strength,
                                                 sec_strength,
                                                 dir,
                                                 pri_damping,
                                                 sec_damping,
                                                 coeff_shift,
                                                 8,
                                                 subsampling_factor);
        } else if (bsize == BLOCK_8X4) {
            svt_av1_cdef_filter_block_8xn_8_neon(dst8,
                                                 dstride,
                                                 in,
                                                 pri_strength,
                                                 sec_strength,
                                                 dir,
                                                 pri_damping,
                                                 sec_damping,
                                                 coeff_shift,
                                                 4,
                                                 subsampling_factor);
        } else {
            svt_av1_cdef_filter_block_4xn_8_neon(
                dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping, sec_damping, coeff_shift, 4, 1);
        }
    } else {
        if (bsize == BLOCK_8X8) {
            svt_av1_cdef_filter_block_8xn_16_neon(dst16,
                                                  dstride,
                                                  in,
                                                  pri_strength,
                                                  sec_strength,
                                                  dir,
                                                  pri_damping,
                                                  sec_damping,
                                                  coeff_shift,
                                                  8,
                                                  subsampling_factor);
        } else if (bsize == BLOCK_4X8) {
            svt_av1_cdef_filter_block_4xn_16_neon(dst16,
                                                  dstride,
                                                  in,
                                                  pri_strength,
                                                  sec_strength,
                                                  dir,
                                                  pri_damping,
                                                  sec_damping,
                                                  coeff_shift,
                                                  8,
                                                  subsampling_factor);
        } else if (bsize == BLOCK_8X4) {
            svt_av1_cdef_filter_block_8xn_16_neon(dst16,
                                                  dstride,
                                                  in,
                                                  pri_strength,
                                                  sec_strength,
                                                  dir,
                                                  pri_damping,
                                                  sec_damping,
                                                  coeff_shift,
                                                  4,
                                                  subsampling_factor);
        } else {
            assert(bsize == BLOCK_4X4);
            svt_av1_cdef_filter_block_4xn_16_neon(
                dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping, sec_damping, coeff_shift, 4, 1);
        }
    }
}
