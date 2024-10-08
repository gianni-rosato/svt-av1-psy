/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */
#include "definitions.h"
#include "aom_dsp_rtcd.h"
#include <tmmintrin.h>
#include "synonyms.h"
#include "synonyms_avx2.h"
#include "convolve_avx2.h"

static INLINE int64_t summary_all_avx2(const __m256i *sum_all) {
    int64_t       sum;
    __m256i       zero      = _mm256_setzero_si256();
    const __m256i sum0_4x64 = _mm256_unpacklo_epi32(*sum_all, zero);
    const __m256i sum1_4x64 = _mm256_unpackhi_epi32(*sum_all, zero);
    const __m256i sum_4x64  = _mm256_add_epi64(sum0_4x64, sum1_4x64);
    const __m128i sum_2x64  = _mm_add_epi64(_mm256_castsi256_si128(sum_4x64), _mm256_extracti128_si256(sum_4x64, 1));
    const __m128i sum_1x64  = _mm_add_epi64(sum_2x64, _mm_srli_si128(sum_2x64, 8));
    xx_storel_64(&sum, sum_1x64);
    return sum;
}

static INLINE void summary_32_avx2(const __m256i *sum32, __m256i *sum) {
    const __m256i sum0_4x64 = _mm256_cvtepu32_epi64(_mm256_castsi256_si128(*sum32));
    const __m256i sum1_4x64 = _mm256_cvtepu32_epi64(_mm256_extracti128_si256(*sum32, 1));
    const __m256i sum_4x64  = _mm256_add_epi64(sum0_4x64, sum1_4x64);
    *sum                    = _mm256_add_epi64(*sum, sum_4x64);
}

static INLINE int64_t summary_4x64_avx2(const __m256i sum_4x64) {
    int64_t       sum;
    const __m128i sum_2x64 = _mm_add_epi64(_mm256_castsi256_si128(sum_4x64), _mm256_extracti128_si256(sum_4x64, 1));
    const __m128i sum_1x64 = _mm_add_epi64(sum_2x64, _mm_srli_si128(sum_2x64, 8));

    xx_storel_64(&sum, sum_1x64);
    return sum;
}
static INLINE void highbd_sse_w16_avx2(__m256i *sum, const uint16_t *a, const uint16_t *b) {
    const __m256i v_a_w = yy_loadu_256(a);
    const __m256i v_b_w = yy_loadu_256(b);
    const __m256i v_d_w = _mm256_sub_epi16(v_a_w, v_b_w);
    *sum                = _mm256_add_epi32(*sum, _mm256_madd_epi16(v_d_w, v_d_w));
}

static INLINE void highbd_sse_w4x4_avx2(__m256i *sum, const uint16_t *a, int a_stride, const uint16_t *b,
                                        int b_stride) {
    const __m128i v_a0  = xx_loadl_64(a);
    const __m128i v_a1  = xx_loadl_64(a + a_stride);
    const __m128i v_a2  = xx_loadl_64(a + a_stride * 2);
    const __m128i v_a3  = xx_loadl_64(a + a_stride * 3);
    const __m128i v_b0  = xx_loadl_64(b);
    const __m128i v_b1  = xx_loadl_64(b + b_stride);
    const __m128i v_b2  = xx_loadl_64(b + b_stride * 2);
    const __m128i v_b3  = xx_loadl_64(b + b_stride * 3);
    const __m256i v_a_w = yy_set_m128i(_mm_unpacklo_epi64(v_a0, v_a1), _mm_unpacklo_epi64(v_a2, v_a3));
    const __m256i v_b_w = yy_set_m128i(_mm_unpacklo_epi64(v_b0, v_b1), _mm_unpacklo_epi64(v_b2, v_b3));
    const __m256i v_d_w = _mm256_sub_epi16(v_a_w, v_b_w);
    *sum                = _mm256_add_epi32(*sum, _mm256_madd_epi16(v_d_w, v_d_w));
}

static INLINE void highbd_sse_w8x2_avx2(__m256i *sum, const uint16_t *a, int a_stride, const uint16_t *b,
                                        int b_stride) {
    const __m256i v_a_w = yy_loadu2_128(a + a_stride, a);
    const __m256i v_b_w = yy_loadu2_128(b + b_stride, b);
    const __m256i v_d_w = _mm256_sub_epi16(v_a_w, v_b_w);
    *sum                = _mm256_add_epi32(*sum, _mm256_madd_epi16(v_d_w, v_d_w));
}
int64_t svt_aom_highbd_sse_avx2(const uint8_t *a8, int a_stride, const uint8_t *b8, int b_stride, int width,
                                int height) {
    int32_t   y   = 0;
    int64_t   sse = 0;
    uint16_t *a   = (uint16_t *)(a8);
    uint16_t *b   = (uint16_t *)(b8);
    __m256i   sum = _mm256_setzero_si256();
    switch (width) {
    case 4:
        do {
            highbd_sse_w4x4_avx2(&sum, a, a_stride, b, b_stride);
            a += a_stride << 2;
            b += b_stride << 2;
            y += 4;
        } while (y < height);
        sse = summary_all_avx2(&sum);
        break;
    case 8:
        do {
            highbd_sse_w8x2_avx2(&sum, a, a_stride, b, b_stride);
            a += a_stride << 1;
            b += b_stride << 1;
            y += 2;
        } while (y < height);
        sse = summary_all_avx2(&sum);
        break;
    case 16:
        do {
            highbd_sse_w16_avx2(&sum, a, b);
            a += a_stride;
            b += b_stride;
            y += 1;
        } while (y < height);
        sse = summary_all_avx2(&sum);
        break;
    case 32:
        do {
            int     l     = 0;
            __m256i sum32 = _mm256_setzero_si256();
            do {
                highbd_sse_w16_avx2(&sum32, a, b);
                highbd_sse_w16_avx2(&sum32, a + 16, b + 16);
                a += a_stride;
                b += b_stride;
                l += 1;
            } while (l < 64 && l < (height - y));
            summary_32_avx2(&sum32, &sum);
            y += 64;
        } while (y < height);
        sse = summary_4x64_avx2(sum);
        break;
    case 64:
        do {
            int     l     = 0;
            __m256i sum32 = _mm256_setzero_si256();
            do {
                highbd_sse_w16_avx2(&sum32, a, b);
                highbd_sse_w16_avx2(&sum32, a + 16 * 1, b + 16 * 1);
                highbd_sse_w16_avx2(&sum32, a + 16 * 2, b + 16 * 2);
                highbd_sse_w16_avx2(&sum32, a + 16 * 3, b + 16 * 3);
                a += a_stride;
                b += b_stride;
                l += 1;
            } while (l < 32 && l < (height - y));
            summary_32_avx2(&sum32, &sum);
            y += 32;
        } while (y < height);
        sse = summary_4x64_avx2(sum);
        break;
    case 128:
        do {
            int     l     = 0;
            __m256i sum32 = _mm256_setzero_si256();
            do {
                highbd_sse_w16_avx2(&sum32, a, b);
                highbd_sse_w16_avx2(&sum32, a + 16 * 1, b + 16 * 1);
                highbd_sse_w16_avx2(&sum32, a + 16 * 2, b + 16 * 2);
                highbd_sse_w16_avx2(&sum32, a + 16 * 3, b + 16 * 3);
                highbd_sse_w16_avx2(&sum32, a + 16 * 4, b + 16 * 4);
                highbd_sse_w16_avx2(&sum32, a + 16 * 5, b + 16 * 5);
                highbd_sse_w16_avx2(&sum32, a + 16 * 6, b + 16 * 6);
                highbd_sse_w16_avx2(&sum32, a + 16 * 7, b + 16 * 7);
                a += a_stride;
                b += b_stride;
                l += 1;
            } while (l < 16 && l < (height - y));
            summary_32_avx2(&sum32, &sum);
            y += 16;
        } while (y < height);
        sse = summary_4x64_avx2(sum);
        break;
    default:
        if (width & 0x7) {
            do {
                int     i     = 0;
                __m256i sum32 = _mm256_setzero_si256();
                do {
                    highbd_sse_w8x2_avx2(&sum32, a + i, a_stride, b + i, b_stride);
                    const uint16_t *a2 = a + i + (a_stride << 1);
                    const uint16_t *b2 = b + i + (b_stride << 1);
                    highbd_sse_w8x2_avx2(&sum32, a2, a_stride, b2, b_stride);
                    i += 8;
                } while (i + 4 < width);
                highbd_sse_w4x4_avx2(&sum32, a + i, a_stride, b + i, b_stride);
                summary_32_avx2(&sum32, &sum);
                a += a_stride << 2;
                b += b_stride << 2;
                y += 4;
            } while (y < height);
        } else {
            do {
                int     l     = 0;
                __m256i sum32 = _mm256_setzero_si256();
                do {
                    int i = 0;
                    do {
                        highbd_sse_w8x2_avx2(&sum32, a + i, a_stride, b + i, b_stride);
                        i += 8;
                    } while (i < width);
                    a += a_stride << 1;
                    b += b_stride << 1;
                    l += 2;
                } while (l < 8 && l < (height - y));
                summary_32_avx2(&sum32, &sum);
                y += 8;
            } while (y < height);
        }
        sse = summary_4x64_avx2(sum);
        break;
    }
    return sse;
}
// CONFIG_AV1_HIGHBITDEPTH
