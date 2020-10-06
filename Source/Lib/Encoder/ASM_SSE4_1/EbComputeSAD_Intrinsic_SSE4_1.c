/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <assert.h>

#include "EbDefinitions.h"
#include "EbComputeSAD_C.h"
#include "EbUtility.h"
#include <smmintrin.h>

#define UPDATE_BEST(s, k, offset)      \
    tem_sum = _mm_extract_epi32(s, k); \
    if (tem_sum < low_sum) {           \
        low_sum = tem_sum;             \
        x_best  = j + offset + k;      \
        y_best  = i;                   \
    }

void svt_ext_sad_calculation_32x32_64x64_sse4_intrin(uint32_t *p_sad16x16, uint32_t *p_best_sad_32x32,
                                                     uint32_t *p_best_sad_64x64,
                                                     uint32_t *p_best_mv32x32, uint32_t *p_best_mv64x64,
                                                     uint32_t mv, uint32_t *p_sad32x32) {
    __m128i xmm_n1, sad32x32_greater_than_bitmask, sad32x32_less_than_or_eq_bitmask, best_sad32x32,
        best_mv_32x32, xmm_mv;
    __m128i sad_16x16_0_7_lo, sad_16x16_0_7_hi, sad_16x16_8_15_lo, sad_16x16_8_15_hi, xmm_sad64x64,
        xmm_sad64x64_total, xmm_p_best_sad_32x32, xmm_p_best_mv_32x32;

    sad_16x16_0_7_lo  = _mm_unpacklo_epi32(_mm_loadu_si128((__m128i *)p_sad16x16),
                                          _mm_loadu_si128((__m128i *)(p_sad16x16 + 4)));
    sad_16x16_0_7_hi  = _mm_unpackhi_epi32(_mm_loadu_si128((__m128i *)p_sad16x16),
                                          _mm_loadu_si128((__m128i *)(p_sad16x16 + 4)));
    sad_16x16_8_15_lo = _mm_unpacklo_epi32(_mm_loadu_si128((__m128i *)(p_sad16x16 + 8)),
                                           _mm_loadu_si128((__m128i *)(p_sad16x16 + 12)));
    sad_16x16_8_15_hi = _mm_unpackhi_epi32(_mm_loadu_si128((__m128i *)(p_sad16x16 + 8)),
                                           _mm_loadu_si128((__m128i *)(p_sad16x16 + 12)));

    xmm_sad64x64 =
        _mm_add_epi32(_mm_add_epi32(_mm_unpacklo_epi64(sad_16x16_0_7_lo, sad_16x16_8_15_lo),
                                    _mm_unpackhi_epi64(sad_16x16_0_7_lo, sad_16x16_8_15_lo)),
                      _mm_add_epi32(_mm_unpacklo_epi64(sad_16x16_0_7_hi, sad_16x16_8_15_hi),
                                    _mm_unpackhi_epi64(sad_16x16_0_7_hi, sad_16x16_8_15_hi)));

    p_sad32x32[0] = _mm_extract_epi32(xmm_sad64x64, 0);
    p_sad32x32[1] = _mm_extract_epi32(xmm_sad64x64, 1);
    p_sad32x32[2] = _mm_extract_epi32(xmm_sad64x64, 2);
    p_sad32x32[3] = _mm_extract_epi32(xmm_sad64x64, 3);

    xmm_sad64x64_total = _mm_add_epi32(_mm_srli_si128(xmm_sad64x64, 8), xmm_sad64x64);

    xmm_sad64x64_total = _mm_add_epi32(_mm_srli_si128(xmm_sad64x64_total, 4), xmm_sad64x64_total);

    xmm_mv = _mm_cvtsi32_si128(mv);
    xmm_mv = _mm_unpacklo_epi32(xmm_mv, xmm_mv);
    xmm_mv = _mm_unpacklo_epi64(xmm_mv, xmm_mv);

    xmm_p_best_sad_32x32 = _mm_loadu_si128((__m128i *)p_best_sad_32x32);
    xmm_p_best_mv_32x32  = _mm_loadu_si128((__m128i *)p_best_mv32x32);

    sad32x32_greater_than_bitmask = _mm_cmpgt_epi32(
        xmm_p_best_sad_32x32, xmm_sad64x64); // _mm_cmplt_epi32(xmm_p_best_sad_32x32, xmm_sad64x64);

    xmm_n1 =
        _mm_cmpeq_epi8(xmm_mv, xmm_mv); // anything compared to itself is equal (get 0xFFFFFFFF)
    sad32x32_less_than_or_eq_bitmask = _mm_sub_epi32(xmm_n1, sad32x32_greater_than_bitmask);

    best_sad32x32 =
        _mm_or_si128(_mm_and_si128(xmm_p_best_sad_32x32, sad32x32_less_than_or_eq_bitmask),
                     _mm_and_si128(xmm_sad64x64, sad32x32_greater_than_bitmask));
    best_mv_32x32 =
        _mm_or_si128(_mm_and_si128(xmm_p_best_mv_32x32, sad32x32_less_than_or_eq_bitmask),
                     _mm_and_si128(xmm_mv, sad32x32_greater_than_bitmask));

    _mm_storeu_si128((__m128i *)p_best_sad_32x32, best_sad32x32);
    _mm_storeu_si128((__m128i *)p_best_mv32x32, best_mv_32x32);

    uint32_t sad64x64 = _mm_cvtsi128_si32(xmm_sad64x64_total);
    if (sad64x64 < p_best_sad_64x64[0]) {
        p_best_sad_64x64[0] = sad64x64;
        p_best_mv64x64[0]   = _mm_cvtsi128_si32(xmm_mv);
    }
}

/*******************************************************************************
 * Requirement: width   = 4, 6, 8, 12, 16, 24, 32, 48 or 64 to use SIMD
 * otherwise C version is used
 * Requirement: block_height <= 64
 * Requirement: block_height % 2 = 0 when width = 4, 6 or 8
*******************************************************************************/
void svt_sad_loop_kernel_sse4_1_intrin(
    uint8_t * src, // input parameter, source samples Ptr
    uint32_t  src_stride, // input parameter, source stride
    uint8_t * ref, // input parameter, reference samples Ptr
    uint32_t  ref_stride, // input parameter, reference stride
    uint32_t  block_height, // input parameter, block height (M)
    uint32_t  block_width, // input parameter, block width (N)
    uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center,
    uint32_t src_stride_raw, // input parameter, source stride (no line skipping)
    int16_t search_area_width, int16_t search_area_height) {
    int16_t        x_best = *x_search_center, y_best = *y_search_center;
    uint32_t       low_sum = 0xffffff;
    uint32_t       tem_sum = 0;
    int16_t        i, j;
    uint32_t       k, l;
    uint32_t       leftover = search_area_width & 7;
    const uint8_t *p_ref, *p_src;
    __m128i        s0, s1, s2, s3, s4, s5, s6, s7, s8 = _mm_set1_epi32(-1);

    if (leftover) {
        for (k = 0; k < leftover; k++) s8 = _mm_slli_si128(s8, 2);
    }

    switch (block_width) {
    case 4:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                p_src = src;
                p_ref = ref + j;
                s3    = _mm_setzero_si128();
                for (k = 0; k + 2 <= block_height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                if (k < block_height) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                s3      = _mm_minpos_epu16(s3);
                tem_sum = _mm_extract_epi16(s3, 0);
                if (tem_sum < low_sum) {
                    low_sum = tem_sum;
                    x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    y_best  = i;
                }
            }

            if (leftover) {
                p_src = src;
                p_ref = ref + j;
                s3    = _mm_setzero_si128();
                for (k = 0; k + 2 <= block_height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                if (k < block_height) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                s3      = _mm_or_si128(s3, s8);
                s3      = _mm_minpos_epu16(s3);
                tem_sum = _mm_extract_epi16(s3, 0);
                if (tem_sum < low_sum) {
                    low_sum = tem_sum;
                    x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    y_best  = i;
                }
            }
            ref += src_stride_raw;
        }
        break;

    case 6:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                p_src = src;
                p_ref = ref + j;
                s3    = _mm_setzero_si128();
                for (k = 0; k + 2 <= block_height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                if (k < block_height) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                DECLARE_ALIGNED(16, uint16_t, tsum[8]);
                memset(tsum, 0, 8 * sizeof(uint16_t));
                p_ref = ref + j;
                for (uint32_t search_area = 0; search_area < 8; search_area++) {
                    for (uint32_t y = 0; y < block_height; y++) {
                        tsum[search_area] +=
                            EB_ABS_DIFF(src[y * src_stride + 4],
                                        p_ref[y * ref_stride + 4]) +
                            EB_ABS_DIFF(src[y * src_stride + 5],
                                        p_ref[y * ref_stride + 5]);
                    }
                    p_ref += 1;
                }
                s4 = _mm_loadu_si128((__m128i *)tsum);
                s3 = _mm_adds_epu16(s3, s4);

                s3      = _mm_minpos_epu16(s3);
                tem_sum = _mm_extract_epi16(s3, 0);
                if (tem_sum < low_sum) {
                    low_sum = tem_sum;
                    x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    y_best  = i;
                }
            }

            if (leftover) {
                p_src = src;
                p_ref = ref + j;
                s3    = _mm_setzero_si128();
                for (k = 0; k + 2 <= block_height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(p_src + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                if (k < block_height) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)p_src);
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                s3      = _mm_or_si128(s3, s8);

                DECLARE_ALIGNED(16, uint16_t, tsum[8]);
                memset(tsum, 0, 8 * sizeof(uint16_t));
                p_ref = ref + j;
                for (uint32_t search_area = 0; search_area < leftover; search_area++) {
                    for (uint32_t y = 0; y < block_height; y++) {
                        tsum[search_area] +=
                            EB_ABS_DIFF(src[y * src_stride + 4],
                                        p_ref[y * ref_stride + 4]) +
                            EB_ABS_DIFF(src[y * src_stride + 5],
                                        p_ref[y * ref_stride + 5]);
                    }
                    p_ref += 1;
                }
                s4 = _mm_loadu_si128((__m128i *)tsum);
                s3 = _mm_adds_epu16(s3, s4);

                s3      = _mm_minpos_epu16(s3);
                tem_sum = _mm_extract_epi16(s3, 0);
                if (tem_sum < low_sum) {
                    low_sum = tem_sum;
                    x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    y_best  = i;
                }
            }
            ref += src_stride_raw;
        }
        break;

    case 8:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                p_src = src;
                p_ref = ref + j;
                s3 = s4 = _mm_setzero_si128();
                for (k = 0; k + 2 <= block_height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                    s2 = _mm_loadl_epi64((__m128i *)p_src);
                    s5 = _mm_loadl_epi64((__m128i *)(p_src + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                if (k < block_height) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s2 = _mm_loadl_epi64((__m128i *)p_src);
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                s3      = _mm_adds_epu16(s3, s4);
                s3      = _mm_minpos_epu16(s3);
                tem_sum = _mm_extract_epi16(s3, 0);
                if (tem_sum < low_sum) {
                    low_sum = tem_sum;
                    x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    y_best  = i;
                }
            }

            if (leftover) {
                p_src = src;
                p_ref = ref + j;
                s3 = s4 = _mm_setzero_si128();
                for (k = 0; k + 2 <= block_height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s1 = _mm_loadu_si128((__m128i *)(p_ref + ref_stride));
                    s2 = _mm_loadl_epi64((__m128i *)p_src);
                    s5 = _mm_loadl_epi64((__m128i *)(p_src + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                if (k < block_height) {
                    s0 = _mm_loadu_si128((__m128i *)p_ref);
                    s2 = _mm_loadl_epi64((__m128i *)p_src);
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    p_src += src_stride << 1;
                    p_ref += ref_stride << 1;
                }

                s3      = _mm_adds_epu16(s3, s4);
                s3      = _mm_or_si128(s3, s8);
                s3      = _mm_minpos_epu16(s3);
                tem_sum = _mm_extract_epi16(s3, 0);
                if (tem_sum < low_sum) {
                    low_sum = tem_sum;
                    x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    y_best  = i;
                }
            }
            ref += src_stride_raw;
        }
        break;

    case 12:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum = _mm_extract_epi16(s3, 0);
                    if (tem_sum < low_sum) {
                        low_sum = tem_sum;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s3      = _mm_or_si128(s3, s8);
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum = _mm_extract_epi16(s3, 0);
                    if (tem_sum < low_sum) {
                        low_sum = tem_sum;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), s2);
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), s5);
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), s5);
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), s2);
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), s5);
                            k  = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s0, 0);
                                    s0      = _mm_srli_si128(s0, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 16:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum = _mm_extract_epi16(s3, 0);
                    if (tem_sum < low_sum) {
                        low_sum = tem_sum;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s3      = _mm_or_si128(s3, s8);
                    s3      = _mm_minpos_epu16(s3);
                    tem_sum = _mm_extract_epi16(s3, 0);
                    if (tem_sum < low_sum) {
                        low_sum = tem_sum;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s6);
                    s4      = _mm_minpos_epu16(s3);
                    s6      = _mm_minpos_epu16(s5);
                    s4      = _mm_unpacklo_epi16(s4, s4);
                    s4      = _mm_unpacklo_epi32(s4, s4);
                    s4      = _mm_unpacklo_epi64(s4, s4);
                    s6      = _mm_unpacklo_epi16(s6, s6);
                    s6      = _mm_unpacklo_epi32(s6, s6);
                    s6      = _mm_unpacklo_epi64(s6, s6);
                    s3      = _mm_sub_epi16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s3);
                    s5      = _mm_sub_epi16(s5, s6);
                    s5      = _mm_minpos_epu16(s5);
                    tem_sum = _mm_extract_epi16(s5, 0);
                    tem_sum += _mm_extract_epi16(s4, 0);
                    tem_sum += _mm_extract_epi16(s6, 0);
                    if (tem_sum < low_sum) {
                        low_sum = tem_sum;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s6);
                    s3      = _mm_or_si128(s3, s8);
                    s5      = _mm_or_si128(s5, s8);
                    s4      = _mm_minpos_epu16(s3);
                    s6      = _mm_minpos_epu16(s5);
                    s4      = _mm_unpacklo_epi16(s4, s4);
                    s4      = _mm_unpacklo_epi32(s4, s4);
                    s4      = _mm_unpacklo_epi64(s4, s4);
                    s6      = _mm_unpacklo_epi16(s6, s6);
                    s6      = _mm_unpacklo_epi32(s6, s6);
                    s6      = _mm_unpacklo_epi64(s6, s6);
                    s3      = _mm_sub_epi16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s3);
                    s5      = _mm_sub_epi16(s5, s6);
                    s5      = _mm_minpos_epu16(s5);
                    tem_sum = _mm_extract_epi16(s5, 0);
                    tem_sum += _mm_extract_epi16(s4, 0);
                    tem_sum += _mm_extract_epi16(s6, 0);
                    if (tem_sum < low_sum) {
                        low_sum = tem_sum;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            k  = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s0, 0);
                                    s0      = _mm_srli_si128(s0, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 24:
        if (block_height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s2 = _mm_loadl_epi64((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s6);
                    s4      = _mm_minpos_epu16(s3);
                    s6      = _mm_minpos_epu16(s5);
                    s4      = _mm_unpacklo_epi16(s4, s4);
                    s4      = _mm_unpacklo_epi32(s4, s4);
                    s4      = _mm_unpacklo_epi64(s4, s4);
                    s6      = _mm_unpacklo_epi16(s6, s6);
                    s6      = _mm_unpacklo_epi32(s6, s6);
                    s6      = _mm_unpacklo_epi64(s6, s6);
                    s3      = _mm_sub_epi16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s3);
                    s5      = _mm_sub_epi16(s5, s6);
                    s5      = _mm_minpos_epu16(s5);
                    tem_sum = _mm_extract_epi16(s5, 0);
                    tem_sum += _mm_extract_epi16(s4, 0);
                    tem_sum += _mm_extract_epi16(s6, 0);
                    if (tem_sum < low_sum) {
                        low_sum = tem_sum;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s2 = _mm_loadl_epi64((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s3      = _mm_adds_epu16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s6);
                    s3      = _mm_or_si128(s3, s8);
                    s5      = _mm_or_si128(s5, s8);
                    s4      = _mm_minpos_epu16(s3);
                    s6      = _mm_minpos_epu16(s5);
                    s4      = _mm_unpacklo_epi16(s4, s4);
                    s4      = _mm_unpacklo_epi32(s4, s4);
                    s4      = _mm_unpacklo_epi64(s4, s4);
                    s6      = _mm_unpacklo_epi16(s6, s6);
                    s6      = _mm_unpacklo_epi32(s6, s6);
                    s6      = _mm_unpacklo_epi64(s6, s6);
                    s3      = _mm_sub_epi16(s3, s4);
                    s5      = _mm_adds_epu16(s5, s3);
                    s5      = _mm_sub_epi16(s5, s6);
                    s5      = _mm_minpos_epu16(s5);
                    tem_sum = _mm_extract_epi16(s5, 0);
                    tem_sum += _mm_extract_epi16(s4, 0);
                    tem_sum += _mm_extract_epi16(s6, 0);
                    if (tem_sum < low_sum) {
                        low_sum = tem_sum;
                        x_best  = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        y_best  = i;
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s2 = _mm_loadl_epi64((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s2 = _mm_loadl_epi64((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            k  = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s0, 0);
                                    s0      = _mm_srli_si128(s0, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 32:
        if (block_height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < block_height; k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0      = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0 = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3 = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            k  = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s0, 0);
                                    s0      = _mm_srli_si128(s0, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < (block_height >> 1); k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < block_height; k++) {
                        s0  = _mm_loadu_si128((__m128i *)p_ref);
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2  = _mm_loadu_si128((__m128i *)p_src);
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(
                        s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0  = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3  = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1  = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4  = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5  = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7  = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6  = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0  = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3  = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            s1  = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9  = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4  = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5  = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0  = _mm_add_epi32(
                                s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(
                                s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < (block_height >> 1); k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < block_height; k++) {
                        s0  = _mm_loadu_si128((__m128i *)p_ref);
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2  = _mm_loadu_si128((__m128i *)p_src);
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(
                        s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0  = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3  = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1  = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4  = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5  = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7  = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6  = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0  = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3  = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            s1  = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9  = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4  = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5  = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0  = _mm_add_epi32(
                                s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(
                                s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s0, 0);
                                    s0      = _mm_srli_si128(s0, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 48:
        if (block_height <= 32) {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < (block_height >> 1); k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < block_height; k++) {
                        s0  = _mm_loadu_si128((__m128i *)p_ref);
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2  = _mm_loadu_si128((__m128i *)p_src);
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(
                        s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0  = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3  = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1  = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4  = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5  = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7  = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6  = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0  = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3  = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            s1  = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9  = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4  = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5  = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0  = _mm_add_epi32(
                                s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(
                                s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < (block_height >> 1); k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < block_height; k++) {
                        s0  = _mm_loadu_si128((__m128i *)p_ref);
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2  = _mm_loadu_si128((__m128i *)p_src);
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(
                        s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0  = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3  = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1  = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4  = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5  = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7  = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6  = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0  = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3  = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            s1  = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9  = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4  = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5  = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0  = _mm_add_epi32(
                                s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(
                                s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s0, 0);
                                    s0      = _mm_srli_si128(s0, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            __m128i s9, s10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k        = 0;
                    while (k < block_height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 21 && k < block_height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i *)p_ref);
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                            s2 = _mm_loadu_si128((__m128i *)p_src);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 32));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 40));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            p_src += src_stride;
                            p_ref += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(
                            s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(
                            s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0      = _mm_packus_epi32(s9, s10);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            UPDATE_BEST(s9, 0, 0);
                            UPDATE_BEST(s9, 1, 0);
                            UPDATE_BEST(s9, 2, 0);
                            UPDATE_BEST(s9, 3, 0);
                            UPDATE_BEST(s10, 0, 4);
                            UPDATE_BEST(s10, 1, 4);
                            UPDATE_BEST(s10, 2, 4);
                            UPDATE_BEST(s10, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k        = 0;
                    while (k < block_height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 21 && k < block_height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i *)p_ref);
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                            s2 = _mm_loadu_si128((__m128i *)p_src);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 32));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 40));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            p_src += src_stride;
                            p_ref += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(
                            s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(
                            s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0      = _mm_packus_epi32(s9, s10);
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s9, 0);
                                    s9      = _mm_srli_si128(s9, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s9 = s10;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    case 64:
        if (block_height <= 32) {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < (block_height >> 1); k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 48));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 56));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 48));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < block_height; k++) {
                        s0  = _mm_loadu_si128((__m128i *)p_ref);
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2  = _mm_loadu_si128((__m128i *)p_src);
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 48));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 56));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 48));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(
                        s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0  = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3  = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1  = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4  = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5  = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7  = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6  = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0  = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3  = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            s1  = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9  = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4  = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5  = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0  = _mm_add_epi32(
                                s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(
                                s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            UPDATE_BEST(s0, 0, 0);
                            UPDATE_BEST(s0, 1, 0);
                            UPDATE_BEST(s0, 2, 0);
                            UPDATE_BEST(s0, 3, 0);
                            UPDATE_BEST(s3, 0, 4);
                            UPDATE_BEST(s3, 1, 4);
                            UPDATE_BEST(s3, 2, 4);
                            UPDATE_BEST(s3, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < (block_height >> 1); k++) {
                        s0 = _mm_loadu_si128((__m128i *)p_ref);
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2 = _mm_loadu_si128((__m128i *)p_src);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i *)(p_ref + 48));
                        s1 = _mm_loadu_si128((__m128i *)(p_ref + 56));
                        s2 = _mm_loadu_si128((__m128i *)(p_src + 48));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < block_height; k++) {
                        s0  = _mm_loadu_si128((__m128i *)p_ref);
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 8));
                        s2  = _mm_loadu_si128((__m128i *)p_src);
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 16));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 24));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 16));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 32));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 40));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 32));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0  = _mm_loadu_si128((__m128i *)(p_ref + 48));
                        s1  = _mm_loadu_si128((__m128i *)(p_ref + 56));
                        s2  = _mm_loadu_si128((__m128i *)(p_src + 48));
                        s9  = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        p_src += src_stride;
                        p_ref += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(
                        s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            s0  = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                            s3  = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                            s1  = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                            s4  = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                            s5  = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                            s7  = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                            s6  = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                            s0  = _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7));
                            s3  = _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6));
                            s1  = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9  = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2  = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4  = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5  = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0  = _mm_add_epi32(
                                s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(
                                s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s0, 0);
                                    s0      = _mm_srli_si128(s0, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s0 = s3;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        } else {
            __m128i s9, s10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    p_src = src;
                    p_ref = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k        = 0;
                    while (k < block_height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 16 && k < block_height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i *)p_ref);
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                            s2 = _mm_loadu_si128((__m128i *)p_src);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 32));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 40));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 48));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 56));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 48));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            p_src += src_stride;
                            p_ref += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(
                            s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(
                            s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0      = _mm_packus_epi32(s9, s10);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            UPDATE_BEST(s9, 0, 0);
                            UPDATE_BEST(s9, 1, 0);
                            UPDATE_BEST(s9, 2, 0);
                            UPDATE_BEST(s9, 3, 0);
                            UPDATE_BEST(s10, 0, 4);
                            UPDATE_BEST(s10, 1, 4);
                            UPDATE_BEST(s10, 2, 4);
                            UPDATE_BEST(s10, 3, 4);
                        }
                    }
                }

                if (leftover) {
                    p_src = src;
                    p_ref = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k        = 0;
                    while (k < block_height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 16 && k < block_height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i *)p_ref);
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 8));
                            s2 = _mm_loadu_si128((__m128i *)p_src);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 16));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 24));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 32));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 40));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i *)(p_ref + 48));
                            s1 = _mm_loadu_si128((__m128i *)(p_ref + 56));
                            s2 = _mm_loadu_si128((__m128i *)(p_src + 48));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            p_src += src_stride;
                            p_ref += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(
                            s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(
                            s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0      = _mm_packus_epi32(s9, s10);
                    s0      = _mm_or_si128(s0, s8);
                    s0      = _mm_minpos_epu16(s0);
                    tem_sum = _mm_extract_epi16(s0, 0);
                    tem_sum &= 0x0000FFFF;
                    if (tem_sum < low_sum) {
                        if (tem_sum != 0xFFFF) { // no overflow
                            low_sum = tem_sum;
                            x_best  = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            y_best  = i;
                        } else {
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    tem_sum = _mm_extract_epi32(s9, 0);
                                    s9      = _mm_srli_si128(s9, 4);
                                    if (tem_sum < low_sum) {
                                        low_sum = tem_sum;
                                        x_best  = (int16_t)(j + leftover - k);
                                        y_best  = i;
                                    }
                                }
                                s9 = s10;
                            }
                        }
                    }
                }
                ref += src_stride_raw;
            }
        }
        break;

    default:
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
                              search_area_width,
                              search_area_height);
        return;
    }

    *best_sad        = low_sum;
    *x_search_center = x_best;
    *y_search_center = y_best;
}
