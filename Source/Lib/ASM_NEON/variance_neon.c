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
#include "aom_dsp_rtcd.h"
#include "var_filter_neon.h"

static INLINE void variance_4xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int h,
                                     uint32_t *sse, int *sum) {
    int16x8_t sum_s16 = vdupq_n_s16(0);
    int32x4_t sse_s32 = vdupq_n_s32(0);

    // Number of rows we can process before 'sum_s16' overflows:
    // 32767 / 255 ~= 128, but we use an 8-wide accumulator; so 256 4-wide rows.
    assert(h <= 256);

    int i = h;
    do {
        uint8x8_t s    = load_unaligned_u8(src, src_stride);
        uint8x8_t r    = load_unaligned_u8(ref, ref_stride);
        int16x8_t diff = vreinterpretq_s16_u16(vsubl_u8(s, r));

        sum_s16 = vaddq_s16(sum_s16, diff);

        sse_s32 = vmlal_s16(sse_s32, vget_low_s16(diff), vget_low_s16(diff));
        sse_s32 = vmlal_s16(sse_s32, vget_high_s16(diff), vget_high_s16(diff));

        src += 2 * src_stride;
        ref += 2 * ref_stride;
        i -= 2;
    } while (i != 0);

    *sum = vaddlvq_s16(sum_s16);
    *sse = (uint32_t)vaddvq_s32(sse_s32);
}

static INLINE void variance_8xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int h,
                                     uint32_t *sse, int *sum) {
    int16x8_t sum_s16    = vdupq_n_s16(0);
    int32x4_t sse_s32[2] = {vdupq_n_s32(0), vdupq_n_s32(0)};

    // Number of rows we can process before 'sum_s16' overflows:
    // 32767 / 255 ~= 128
    assert(h <= 128);

    int i = h;
    do {
        uint8x8_t s    = vld1_u8(src);
        uint8x8_t r    = vld1_u8(ref);
        int16x8_t diff = vreinterpretq_s16_u16(vsubl_u8(s, r));

        sum_s16 = vaddq_s16(sum_s16, diff);

        sse_s32[0] = vmlal_s16(sse_s32[0], vget_low_s16(diff), vget_low_s16(diff));
        sse_s32[1] = vmlal_s16(sse_s32[1], vget_high_s16(diff), vget_high_s16(diff));

        src += src_stride;
        ref += ref_stride;
    } while (--i != 0);

    *sum = vaddlvq_s16(sum_s16);
    *sse = (uint32_t)vaddvq_s32(vaddq_s32(sse_s32[0], sse_s32[1]));
}

static INLINE void variance_16xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int h,
                                      uint32_t *sse, int *sum) {
    int16x8_t sum_s16[2] = {vdupq_n_s16(0), vdupq_n_s16(0)};
    int32x4_t sse_s32[2] = {vdupq_n_s32(0), vdupq_n_s32(0)};

    // Number of rows we can process before 'sum_s16' accumulators overflow:
    // 32767 / 255 ~= 128, so 128 16-wide rows.
    assert(h <= 128);

    int i = h;
    do {
        uint8x16_t s = vld1q_u8(src);
        uint8x16_t r = vld1q_u8(ref);

        int16x8_t diff_l = vreinterpretq_s16_u16(vsubl_u8(vget_low_u8(s), vget_low_u8(r)));
        int16x8_t diff_h = vreinterpretq_s16_u16(vsubl_u8(vget_high_u8(s), vget_high_u8(r)));

        sum_s16[0] = vaddq_s16(sum_s16[0], diff_l);
        sum_s16[1] = vaddq_s16(sum_s16[1], diff_h);

        sse_s32[0] = vmlal_s16(sse_s32[0], vget_low_s16(diff_l), vget_low_s16(diff_l));
        sse_s32[1] = vmlal_s16(sse_s32[1], vget_high_s16(diff_l), vget_high_s16(diff_l));
        sse_s32[0] = vmlal_s16(sse_s32[0], vget_low_s16(diff_h), vget_low_s16(diff_h));
        sse_s32[1] = vmlal_s16(sse_s32[1], vget_high_s16(diff_h), vget_high_s16(diff_h));

        src += src_stride;
        ref += ref_stride;
    } while (--i != 0);

    *sum = vaddlvq_s16(vaddq_s16(sum_s16[0], sum_s16[1]));
    *sse = (uint32_t)vaddvq_s32(vaddq_s32(sse_s32[0], sse_s32[1]));
}

static INLINE void variance_large_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int w,
                                       int h, int h_limit, uint32_t *sse, int *sum) {
    int32x4_t sum_s32 = vdupq_n_s32(0);

    int32x4_t sse_s32_0 = vdupq_n_s32(0);
    int32x4_t sse_s32_1 = vdupq_n_s32(0);
    int32x4_t sse_s32_2 = vdupq_n_s32(0);
    int32x4_t sse_s32_3 = vdupq_n_s32(0);

    // 'h_limit' is the number of 'w'-width rows we can process before our 16-bit
    // accumulator overflows. After hitting this limit we accumulate into 32-bit
    // elements.
    int h_tmp = h > h_limit ? h_limit : h;

    int i = 0;
    do {
        int16x8_t sum_s16_0 = vdupq_n_s16(0);
        int16x8_t sum_s16_1 = vdupq_n_s16(0);

        do {
            int j = 0;
            do {
                uint8x16_t s = vld1q_u8(src + j);
                uint8x16_t r = vld1q_u8(ref + j);

                const int16x8_t diff_l = vreinterpretq_s16_u16(vsubl_u8(vget_low_u8(s), vget_low_u8(r)));
                const int16x8_t diff_h = vreinterpretq_s16_u16(vsubl_u8(vget_high_u8(s), vget_high_u8(r)));

                sum_s16_0 = vaddq_s16(sum_s16_0, diff_l);
                sum_s16_1 = vaddq_s16(sum_s16_1, diff_h);

                sse_s32_0 = vmlal_s16(sse_s32_0, vget_low_s16(diff_l), vget_low_s16(diff_l));
                sse_s32_1 = vmlal_s16(sse_s32_1, vget_high_s16(diff_l), vget_high_s16(diff_l));
                sse_s32_2 = vmlal_s16(sse_s32_2, vget_low_s16(diff_h), vget_low_s16(diff_h));
                sse_s32_3 = vmlal_s16(sse_s32_3, vget_high_s16(diff_h), vget_high_s16(diff_h));

                j += 16;
            } while (j < w);

            src += src_stride;
            ref += ref_stride;
            i++;
        } while (i < h_tmp);

        sum_s32 = vpadalq_s16(sum_s32, sum_s16_0);
        sum_s32 = vpadalq_s16(sum_s32, sum_s16_1);

        h_tmp += h_limit;
    } while (i < h);

    *sum = vaddvq_s32(sum_s32);
    *sse = (uint32_t)vaddvq_s32(vaddq_s32(vaddq_s32(sse_s32_0, sse_s32_1), vaddq_s32(sse_s32_2, sse_s32_3)));
}

static INLINE void variance_32xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int h,
                                      uint32_t *sse, int *sum) {
    variance_large_neon(src, src_stride, ref, ref_stride, 32, h, 64, sse, sum);
}

static INLINE void variance_64xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int h,
                                      uint32_t *sse, int *sum) {
    variance_large_neon(src, src_stride, ref, ref_stride, 64, h, 32, sse, sum);
}

static INLINE void variance_128xh_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, int h,
                                       uint32_t *sse, int *sum) {
    variance_large_neon(src, src_stride, ref, ref_stride, 128, h, 16, sse, sum);
}

unsigned int svt_aom_variance4x4_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                      unsigned int *sse) {
    int sum;
    variance_4xh_neon(src, src_stride, ref, ref_stride, 4, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 4);
}
unsigned int svt_aom_variance4x8_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                      unsigned int *sse) {
    int sum;
    variance_4xh_neon(src, src_stride, ref, ref_stride, 8, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 5);
}
unsigned int svt_aom_variance4x16_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                       unsigned int *sse) {
    int sum;
    variance_4xh_neon(src, src_stride, ref, ref_stride, 16, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 6);
}

unsigned int svt_aom_variance8x4_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                      unsigned int *sse) {
    int sum;
    variance_8xh_neon(src, src_stride, ref, ref_stride, 4, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 5);
}
unsigned int svt_aom_variance8x8_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                      unsigned int *sse) {
    int sum;
    variance_8xh_neon(src, src_stride, ref, ref_stride, 8, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 6);
}
unsigned int svt_aom_variance8x16_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                       unsigned int *sse) {
    int sum;
    variance_8xh_neon(src, src_stride, ref, ref_stride, 16, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 7);
}
unsigned int svt_aom_variance8x32_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                       unsigned int *sse) {
    int sum;
    variance_8xh_neon(src, src_stride, ref, ref_stride, 32, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 8);
}

unsigned int svt_aom_variance16x4_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                       unsigned int *sse) {
    int sum;
    variance_16xh_neon(src, src_stride, ref, ref_stride, 4, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 6);
}
unsigned int svt_aom_variance16x8_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                       unsigned int *sse) {
    int sum;
    variance_16xh_neon(src, src_stride, ref, ref_stride, 8, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 7);
}
unsigned int svt_aom_variance16x16_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_16xh_neon(src, src_stride, ref, ref_stride, 16, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 8);
}
unsigned int svt_aom_variance16x32_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_16xh_neon(src, src_stride, ref, ref_stride, 32, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 9);
}
unsigned int svt_aom_variance16x64_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_16xh_neon(src, src_stride, ref, ref_stride, 64, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 10);
}

unsigned int svt_aom_variance32x8_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                       unsigned int *sse) {
    int sum;
    variance_32xh_neon(src, src_stride, ref, ref_stride, 8, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 8);
}
unsigned int svt_aom_variance32x16_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_32xh_neon(src, src_stride, ref, ref_stride, 16, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 9);
}
unsigned int svt_aom_variance32x32_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_32xh_neon(src, src_stride, ref, ref_stride, 32, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 10);
}
unsigned int svt_aom_variance32x64_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_32xh_neon(src, src_stride, ref, ref_stride, 64, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 11);
}

unsigned int svt_aom_variance64x16_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_64xh_neon(src, src_stride, ref, ref_stride, 16, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 10);
}
unsigned int svt_aom_variance64x32_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_64xh_neon(src, src_stride, ref, ref_stride, 32, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 11);
}
unsigned int svt_aom_variance64x64_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                        unsigned int *sse) {
    int sum;
    variance_64xh_neon(src, src_stride, ref, ref_stride, 64, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 12);
}
unsigned int svt_aom_variance64x128_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                         unsigned int *sse) {
    int sum;
    variance_64xh_neon(src, src_stride, ref, ref_stride, 128, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 13);
}

unsigned int svt_aom_variance128x64_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                         unsigned int *sse) {
    int sum;
    variance_128xh_neon(src, src_stride, ref, ref_stride, 64, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 13);
}
unsigned int svt_aom_variance128x128_neon(const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride,
                                          unsigned int *sse) {
    int sum;
    variance_128xh_neon(src, src_stride, ref, ref_stride, 128, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) >> 14);
}

unsigned int svt_aom_sub_pixel_variance4x4_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[4 * (4 + 2)];
    uint8_t tmp1[4 * 4];
    var_filter_block2d_bil_w4(src, tmp0, src_stride, 1, (4 + 2), xoffset);
    var_filter_block2d_bil_w4(tmp0, tmp1, 4, 4, 4, yoffset);
    return svt_aom_variance4x4(tmp1, 4, ref, ref_stride, sse);
}
unsigned int svt_aom_sub_pixel_variance4x8_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[4 * (8 + 2)];
    uint8_t tmp1[4 * 8];
    var_filter_block2d_bil_w4(src, tmp0, src_stride, 1, (8 + 2), xoffset);
    var_filter_block2d_bil_w4(tmp0, tmp1, 4, 4, 8, yoffset);
    return svt_aom_variance4x8(tmp1, 4, ref, ref_stride, sse);
}
unsigned int svt_aom_sub_pixel_variance4x16_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                 const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[4 * (16 + 2)];
    uint8_t tmp1[4 * 16];
    var_filter_block2d_bil_w4(src, tmp0, src_stride, 1, (16 + 2), xoffset);
    var_filter_block2d_bil_w4(tmp0, tmp1, 4, 4, 16, yoffset);
    return svt_aom_variance4x16(tmp1, 4, ref, ref_stride, sse);
}

unsigned int svt_aom_sub_pixel_variance8x4_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[8 * (4 + 1)];
    uint8_t tmp1[8 * 4];
    var_filter_block2d_bil_w8(src, tmp0, src_stride, 1, (4 + 1), xoffset);
    var_filter_block2d_bil_w8(tmp0, tmp1, 8, 8, 4, yoffset);
    return svt_aom_variance8x4(tmp1, 8, ref, ref_stride, sse);
}
unsigned int svt_aom_sub_pixel_variance8x8_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[8 * (8 + 1)];
    uint8_t tmp1[8 * 8];
    var_filter_block2d_bil_w8(src, tmp0, src_stride, 1, (8 + 1), xoffset);
    var_filter_block2d_bil_w8(tmp0, tmp1, 8, 8, 8, yoffset);
    return svt_aom_variance8x8(tmp1, 8, ref, ref_stride, sse);
}
unsigned int svt_aom_sub_pixel_variance8x16_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                 const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[8 * (16 + 1)];
    uint8_t tmp1[8 * 16];
    var_filter_block2d_bil_w8(src, tmp0, src_stride, 1, (16 + 1), xoffset);
    var_filter_block2d_bil_w8(tmp0, tmp1, 8, 8, 16, yoffset);
    return svt_aom_variance8x16(tmp1, 8, ref, ref_stride, sse);
}
unsigned int svt_aom_sub_pixel_variance8x32_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                 const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[8 * (32 + 1)];
    uint8_t tmp1[8 * 32];
    var_filter_block2d_bil_w8(src, tmp0, src_stride, 1, (32 + 1), xoffset);
    var_filter_block2d_bil_w8(tmp0, tmp1, 8, 8, 32, yoffset);
    return svt_aom_variance8x32(tmp1, 8, ref, ref_stride, sse);
}

unsigned int svt_aom_sub_pixel_variance16x8_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                 const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[16 * (8 + 1)];
    uint8_t tmp1[16 * 8];
    var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, (8 + 1), xoffset);
    var_filter_block2d_bil_w16(tmp0, tmp1, 16, 16, 8, yoffset);
    return svt_aom_variance16x8(tmp1, 16, ref, ref_stride, sse);
}
unsigned int svt_aom_sub_pixel_variance16x4_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                 const uint8_t *ref, int ref_stride, uint32_t *sse) {
    uint8_t tmp0[16 * (4 + 1)];
    uint8_t tmp1[16 * 4];
    var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, (4 + 1), xoffset);
    var_filter_block2d_bil_w16(tmp0, tmp1, 16, 16, 4, yoffset);
    return svt_aom_variance16x4(tmp1, 16, ref, ref_stride, sse);
}
unsigned int svt_aom_sub_pixel_variance16x16_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance16x16(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[16 * 16];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 16, 16);
            return svt_aom_variance16x16(tmp, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp[16 * 16];
            var_filter_block2d_bil_w16(src, tmp, src_stride, src_stride, 16, yoffset);
            return svt_aom_variance16x16(tmp, 16, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[16 * (16 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, 16);
            return svt_aom_variance16x16(tmp0, 16, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[16 * (16 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, (16 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 16, 16, 16, 16);
            return svt_aom_variance16x16(tmp1, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[16 * (16 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, (16 + 1));
            var_filter_block2d_bil_w16(tmp0, tmp1, 16, 16, 16, yoffset);
            return svt_aom_variance16x16(tmp1, 16, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[16 * (16 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, 16, xoffset);
            return svt_aom_variance16x16(tmp0, 16, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[16 * 16];
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, (16 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 16, 16, 16, 16);
            return svt_aom_variance16x16(tmp1, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[16 * 16];
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, (16 + 1), xoffset);
            var_filter_block2d_bil_w16(tmp0, tmp1, 16, 16, 16, yoffset);
            return svt_aom_variance16x16(tmp1, 16, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance16x32_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance16x32(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[16 * 32];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 16, 32);
            return svt_aom_variance16x32(tmp, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp[16 * 32];
            var_filter_block2d_bil_w16(src, tmp, src_stride, src_stride, 32, yoffset);
            return svt_aom_variance16x32(tmp, 16, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[16 * (32 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, 32);
            return svt_aom_variance16x32(tmp0, 16, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[16 * (32 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, (32 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 16, 16, 16, 32);
            return svt_aom_variance16x32(tmp1, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[16 * (32 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, (32 + 1));
            var_filter_block2d_bil_w16(tmp0, tmp1, 16, 16, 32, yoffset);
            return svt_aom_variance16x32(tmp1, 16, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[16 * (32 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, 32, xoffset);
            return svt_aom_variance16x32(tmp0, 16, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[16 * 32];
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, (32 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 16, 16, 16, 32);
            return svt_aom_variance16x32(tmp1, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[16 * 32];
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, (32 + 1), xoffset);
            var_filter_block2d_bil_w16(tmp0, tmp1, 16, 16, 32, yoffset);
            return svt_aom_variance16x32(tmp1, 16, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance16x64_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance16x64(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[16 * 64];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 16, 64);
            return svt_aom_variance16x64(tmp, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp[16 * 64];
            var_filter_block2d_bil_w16(src, tmp, src_stride, src_stride, 64, yoffset);
            return svt_aom_variance16x64(tmp, 16, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[16 * (64 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, 64);
            return svt_aom_variance16x64(tmp0, 16, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[16 * (64 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, (64 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 16, 16, 16, 64);
            return svt_aom_variance16x64(tmp1, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[16 * (64 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 16, (64 + 1));
            var_filter_block2d_bil_w16(tmp0, tmp1, 16, 16, 64, yoffset);
            return svt_aom_variance16x64(tmp1, 16, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[16 * (64 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, 64, xoffset);
            return svt_aom_variance16x64(tmp0, 16, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[16 * 64];
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, (64 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 16, 16, 16, 64);
            return svt_aom_variance16x64(tmp1, 16, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[16 * 64];
            var_filter_block2d_bil_w16(src, tmp0, src_stride, 1, (64 + 1), xoffset);
            var_filter_block2d_bil_w16(tmp0, tmp1, 16, 16, 64, yoffset);
            return svt_aom_variance16x64(tmp1, 16, ref, ref_stride, sse);
        }
    }
}

unsigned int svt_aom_sub_pixel_variance32x8_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                 const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance32x8(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[32 * 8];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 32, 8);
            return svt_aom_variance32x8(tmp, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp[32 * 8];
            var_filter_block2d_bil_w32(src, tmp, src_stride, src_stride, 8, yoffset);
            return svt_aom_variance32x8(tmp, 32, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[32 * (8 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, 8);
            return svt_aom_variance32x8(tmp0, 32, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[32 * (8 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, (8 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 32, 32, 32, 8);
            return svt_aom_variance32x8(tmp1, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[32 * (8 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, (8 + 1));
            var_filter_block2d_bil_w32(tmp0, tmp1, 32, 32, 8, yoffset);
            return svt_aom_variance32x8(tmp1, 32, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[32 * (8 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, 8, xoffset);
            return svt_aom_variance32x8(tmp0, 32, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[32 * 8];
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, (8 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 32, 32, 32, 8);
            return svt_aom_variance32x8(tmp1, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[32 * 8];
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, (8 + 1), xoffset);
            var_filter_block2d_bil_w32(tmp0, tmp1, 32, 32, 8, yoffset);
            return svt_aom_variance32x8(tmp1, 32, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance32x16_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance32x16(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[32 * 16];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 32, 16);
            return svt_aom_variance32x16(tmp, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp[32 * 16];
            var_filter_block2d_bil_w32(src, tmp, src_stride, src_stride, 16, yoffset);
            return svt_aom_variance32x16(tmp, 32, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[32 * (16 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, 16);
            return svt_aom_variance32x16(tmp0, 32, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[32 * (16 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, (16 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 32, 32, 32, 16);
            return svt_aom_variance32x16(tmp1, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[32 * (16 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, (16 + 1));
            var_filter_block2d_bil_w32(tmp0, tmp1, 32, 32, 16, yoffset);
            return svt_aom_variance32x16(tmp1, 32, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[32 * (16 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, 16, xoffset);
            return svt_aom_variance32x16(tmp0, 32, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[32 * 16];
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, (16 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 32, 32, 32, 16);
            return svt_aom_variance32x16(tmp1, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[32 * 16];
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, (16 + 1), xoffset);
            var_filter_block2d_bil_w32(tmp0, tmp1, 32, 32, 16, yoffset);
            return svt_aom_variance32x16(tmp1, 32, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance32x32_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance32x32(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[32 * 32];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 32, 32);
            return svt_aom_variance32x32(tmp, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp[32 * 32];
            var_filter_block2d_bil_w32(src, tmp, src_stride, src_stride, 32, yoffset);
            return svt_aom_variance32x32(tmp, 32, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[32 * (32 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, 32);
            return svt_aom_variance32x32(tmp0, 32, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[32 * (32 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, (32 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 32, 32, 32, 32);
            return svt_aom_variance32x32(tmp1, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[32 * (32 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, (32 + 1));
            var_filter_block2d_bil_w32(tmp0, tmp1, 32, 32, 32, yoffset);
            return svt_aom_variance32x32(tmp1, 32, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[32 * (32 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, 32, xoffset);
            return svt_aom_variance32x32(tmp0, 32, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[32 * 32];
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, (32 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 32, 32, 32, 32);
            return svt_aom_variance32x32(tmp1, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[32 * 32];
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, (32 + 1), xoffset);
            var_filter_block2d_bil_w32(tmp0, tmp1, 32, 32, 32, yoffset);
            return svt_aom_variance32x32(tmp1, 32, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance32x64_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance32x64(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[32 * 64];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 32, 64);
            return svt_aom_variance32x64(tmp, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp[32 * 64];
            var_filter_block2d_bil_w32(src, tmp, src_stride, src_stride, 64, yoffset);
            return svt_aom_variance32x64(tmp, 32, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[32 * (64 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, 64);
            return svt_aom_variance32x64(tmp0, 32, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[32 * (64 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, (64 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 32, 32, 32, 64);
            return svt_aom_variance32x64(tmp1, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[32 * (64 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 32, (64 + 1));
            var_filter_block2d_bil_w32(tmp0, tmp1, 32, 32, 64, yoffset);
            return svt_aom_variance32x64(tmp1, 32, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[32 * (64 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, 64, xoffset);
            return svt_aom_variance32x64(tmp0, 32, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[32 * 64];
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, (64 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 32, 32, 32, 64);
            return svt_aom_variance32x64(tmp1, 32, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[32 * 64];
            var_filter_block2d_bil_w32(src, tmp0, src_stride, 1, (64 + 1), xoffset);
            var_filter_block2d_bil_w32(tmp0, tmp1, 32, 32, 64, yoffset);
            return svt_aom_variance32x64(tmp1, 32, ref, ref_stride, sse);
        }
    }
}

unsigned int svt_aom_sub_pixel_variance64x16_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance64x16(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[64 * 16];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 64, 16);
            return svt_aom_variance64x16(tmp, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp[64 * 16];
            var_filter_block2d_bil_w64(src, tmp, src_stride, src_stride, 16, yoffset);
            return svt_aom_variance64x16(tmp, 64, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[64 * (16 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, 16);
            return svt_aom_variance64x16(tmp0, 64, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[64 * (16 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, (16 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 64, 64, 64, 16);
            return svt_aom_variance64x16(tmp1, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[64 * (16 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, (16 + 1));
            var_filter_block2d_bil_w64(tmp0, tmp1, 64, 64, 16, yoffset);
            return svt_aom_variance64x16(tmp1, 64, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[64 * (16 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, 16, xoffset);
            return svt_aom_variance64x16(tmp0, 64, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[64 * 16];
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, (16 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 64, 64, 64, 16);
            return svt_aom_variance64x16(tmp1, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[64 * 16];
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, (16 + 1), xoffset);
            var_filter_block2d_bil_w64(tmp0, tmp1, 64, 64, 16, yoffset);
            return svt_aom_variance64x16(tmp1, 64, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance64x32_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance64x32(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[64 * 32];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 64, 32);
            return svt_aom_variance64x32(tmp, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp[64 * 32];
            var_filter_block2d_bil_w64(src, tmp, src_stride, src_stride, 32, yoffset);
            return svt_aom_variance64x32(tmp, 64, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[64 * (32 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, 32);
            return svt_aom_variance64x32(tmp0, 64, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[64 * (32 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, (32 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 64, 64, 64, 32);
            return svt_aom_variance64x32(tmp1, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[64 * (32 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, (32 + 1));
            var_filter_block2d_bil_w64(tmp0, tmp1, 64, 64, 32, yoffset);
            return svt_aom_variance64x32(tmp1, 64, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[64 * (32 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, 32, xoffset);
            return svt_aom_variance64x32(tmp0, 64, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[64 * 32];
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, (32 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 64, 64, 64, 32);
            return svt_aom_variance64x32(tmp1, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[64 * 32];
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, (32 + 1), xoffset);
            var_filter_block2d_bil_w64(tmp0, tmp1, 64, 64, 32, yoffset);
            return svt_aom_variance64x32(tmp1, 64, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance64x64_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                  const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance64x64(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[64 * 64];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 64, 64);
            return svt_aom_variance64x64(tmp, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp[64 * 64];
            var_filter_block2d_bil_w64(src, tmp, src_stride, src_stride, 64, yoffset);
            return svt_aom_variance64x64(tmp, 64, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[64 * (64 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, 64);
            return svt_aom_variance64x64(tmp0, 64, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[64 * (64 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, (64 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 64, 64, 64, 64);
            return svt_aom_variance64x64(tmp1, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[64 * (64 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, (64 + 1));
            var_filter_block2d_bil_w64(tmp0, tmp1, 64, 64, 64, yoffset);
            return svt_aom_variance64x64(tmp1, 64, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[64 * (64 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, 64, xoffset);
            return svt_aom_variance64x64(tmp0, 64, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[64 * 64];
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, (64 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 64, 64, 64, 64);
            return svt_aom_variance64x64(tmp1, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[64 * 64];
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, (64 + 1), xoffset);
            var_filter_block2d_bil_w64(tmp0, tmp1, 64, 64, 64, yoffset);
            return svt_aom_variance64x64(tmp1, 64, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance64x128_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                   const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance64x128(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[64 * 128];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 64, 128);
            return svt_aom_variance64x128(tmp, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp[64 * 128];
            var_filter_block2d_bil_w64(src, tmp, src_stride, src_stride, 128, yoffset);
            return svt_aom_variance64x128(tmp, 64, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[64 * (128 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, 128);
            return svt_aom_variance64x128(tmp0, 64, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[64 * (128 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, (128 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 64, 64, 64, 128);
            return svt_aom_variance64x128(tmp1, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[64 * (128 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 64, (128 + 1));
            var_filter_block2d_bil_w64(tmp0, tmp1, 64, 64, 128, yoffset);
            return svt_aom_variance64x128(tmp1, 64, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[64 * (128 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, 128, xoffset);
            return svt_aom_variance64x128(tmp0, 64, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[64 * 128];
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, (128 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 64, 64, 64, 128);
            return svt_aom_variance64x128(tmp1, 64, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[64 * 128];
            var_filter_block2d_bil_w64(src, tmp0, src_stride, 1, (128 + 1), xoffset);
            var_filter_block2d_bil_w64(tmp0, tmp1, 64, 64, 128, yoffset);
            return svt_aom_variance64x128(tmp1, 64, ref, ref_stride, sse);
        }
    }
}

unsigned int svt_aom_sub_pixel_variance128x64_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                   const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance128x64(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[128 * 64];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 128, 64);
            return svt_aom_variance128x64(tmp, 128, ref, ref_stride, sse);
        } else {
            uint8_t tmp[128 * 64];
            var_filter_block2d_bil_w128(src, tmp, src_stride, src_stride, 64, yoffset);
            return svt_aom_variance128x64(tmp, 128, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[128 * (64 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 128, 64);
            return svt_aom_variance128x64(tmp0, 128, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[128 * (64 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 128, (64 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 128, 128, 128, 64);
            return svt_aom_variance128x64(tmp1, 128, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[128 * (64 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 128, (64 + 1));
            var_filter_block2d_bil_w128(tmp0, tmp1, 128, 128, 64, yoffset);
            return svt_aom_variance128x64(tmp1, 128, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[128 * (64 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w128(src, tmp0, src_stride, 1, 64, xoffset);
            return svt_aom_variance128x64(tmp0, 128, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[128 * 64];
            var_filter_block2d_bil_w128(src, tmp0, src_stride, 1, (64 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 128, 128, 128, 64);
            return svt_aom_variance128x64(tmp1, 128, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[128 * 64];
            var_filter_block2d_bil_w128(src, tmp0, src_stride, 1, (64 + 1), xoffset);
            var_filter_block2d_bil_w128(tmp0, tmp1, 128, 128, 64, yoffset);
            return svt_aom_variance128x64(tmp1, 128, ref, ref_stride, sse);
        }
    }
}
unsigned int svt_aom_sub_pixel_variance128x128_neon(const uint8_t *src, int src_stride, int xoffset, int yoffset,
                                                    const uint8_t *ref, int ref_stride, unsigned int *sse) {
    if (xoffset == 0) {
        if (yoffset == 0) {
            return svt_aom_variance128x128(src, src_stride, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp[128 * 128];
            var_filter_block2d_avg(src, tmp, src_stride, src_stride, 128, 128);
            return svt_aom_variance128x128(tmp, 128, ref, ref_stride, sse);
        } else {
            uint8_t tmp[128 * 128];
            var_filter_block2d_bil_w128(src, tmp, src_stride, src_stride, 128, yoffset);
            return svt_aom_variance128x128(tmp, 128, ref, ref_stride, sse);
        }
    } else if (xoffset == 4) {
        uint8_t tmp0[128 * (128 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 128, 128);
            return svt_aom_variance128x128(tmp0, 128, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[128 * (128 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 128, (128 + 1));
            var_filter_block2d_avg(tmp0, tmp1, 128, 128, 128, 128);
            return svt_aom_variance128x128(tmp1, 128, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[128 * (128 + 1)];
            var_filter_block2d_avg(src, tmp0, src_stride, 1, 128, (128 + 1));
            var_filter_block2d_bil_w128(tmp0, tmp1, 128, 128, 128, yoffset);
            return svt_aom_variance128x128(tmp1, 128, ref, ref_stride, sse);
        }
    } else {
        uint8_t tmp0[128 * (128 + 1)];
        if (yoffset == 0) {
            var_filter_block2d_bil_w128(src, tmp0, src_stride, 1, 128, xoffset);
            return svt_aom_variance128x128(tmp0, 128, ref, ref_stride, sse);
        } else if (yoffset == 4) {
            uint8_t tmp1[128 * 128];
            var_filter_block2d_bil_w128(src, tmp0, src_stride, 1, (128 + 1), xoffset);
            var_filter_block2d_avg(tmp0, tmp1, 128, 128, 128, 128);
            return svt_aom_variance128x128(tmp1, 128, ref, ref_stride, sse);
        } else {
            uint8_t tmp1[128 * 128];
            var_filter_block2d_bil_w128(src, tmp0, src_stride, 1, (128 + 1), xoffset);
            var_filter_block2d_bil_w128(tmp0, tmp1, 128, 128, 128, yoffset);
            return svt_aom_variance128x128(tmp1, 128, ref, ref_stride, sse);
        }
    }
}
