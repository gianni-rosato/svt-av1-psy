/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <assert.h>

#include "EbComputeSAD_SSE4_1.h"
#include "EbDefinitions.h"
#include "smmintrin.h"


#define UPDATE_BEST(s, k, offset) \
  temSum = _mm_extract_epi32(s, k); \
  if (temSum < lowSum) { \
    lowSum = temSum; \
    xBest = j + offset + k; \
    yBest = i; \
  }

void ext_sad_calculation_8x8_16x16_sse4_intrin(
    uint8_t   *src,
    uint32_t   src_stride,
    uint8_t   *ref,
    uint32_t   ref_stride,
    uint32_t  *p_best_sad8x8,
    uint32_t  *p_best_sad16x16,
    uint32_t  *p_best_mv8x8,
    uint32_t  *p_best_mv16x16,
    uint32_t   mv,
    uint32_t  *p_sad16x16,
    uint32_t  *p_sad8x8)
{
    __m128i xmm_sad16x16, xmm_sad8x8_0_1, xmm_sad8x8_2_3, xmm_sad16x16_total, sad8x8_0_3, sad8x8_less_than_bitmask, xmm_N1;
    __m128i sad8x8_greater_or_eq_bitmask, BestMV8x8, BestSad8x8, xmm_pBestSad8x8, xmm_pBestMV8x8, xmm_mv;

    uint32_t sad8x8_0, sad8x8_1, sad8x8_2, sad8x8_3;
    uint32_t sad16x16;

    (void)sad8x8_0;
    (void)sad8x8_1;
    (void)sad8x8_2;
    (void)sad8x8_3;



    src_stride <<= 1;
    ref_stride <<= 1;

    //sad8x8_0, sad8x8_1

    xmm_sad8x8_0_1 = _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i*)src), _mm_loadu_si128((__m128i*)ref)),
        _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src + src_stride)), _mm_loadu_si128((__m128i*)(ref + ref_stride))));
    xmm_sad8x8_0_1 = _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i*)(src + (2 * src_stride))), _mm_loadu_si128((__m128i*)(ref + (2 * ref_stride)))), xmm_sad8x8_0_1);
    xmm_sad8x8_0_1 = _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i*)(src + (3 * src_stride))), _mm_loadu_si128((__m128i*)(ref + (3 * ref_stride)))), xmm_sad8x8_0_1);

    src += src_stride << 2;
    ref += ref_stride << 2;

    p_sad8x8[0] = _mm_extract_epi32(xmm_sad8x8_0_1, 0) << 1;
    p_sad8x8[1] = _mm_extract_epi32(xmm_sad8x8_0_1, 2) << 1;

    //sad8x8_2, sad8x8_3

    xmm_sad8x8_2_3 = _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i*)src), _mm_loadu_si128((__m128i*)ref)),
        _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src + src_stride)), _mm_loadu_si128((__m128i*)(ref + ref_stride))));
    xmm_sad8x8_2_3 = _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i*)(src + (2 * src_stride))), _mm_loadu_si128((__m128i*)(ref + (2 * ref_stride)))), xmm_sad8x8_2_3);
    xmm_sad8x8_2_3 = _mm_add_epi32(_mm_sad_epu8(_mm_loadu_si128((__m128i*)(src + (3 * src_stride))), _mm_loadu_si128((__m128i*)(ref + (3 * ref_stride)))), xmm_sad8x8_2_3);


    p_sad8x8[2] = _mm_extract_epi32(xmm_sad8x8_2_3, 0) << 1;
    p_sad8x8[3] = _mm_extract_epi32(xmm_sad8x8_2_3, 2) << 1;


    sad8x8_0_3 = _mm_slli_epi32(_mm_packs_epi32(xmm_sad8x8_0_1, xmm_sad8x8_2_3), 1);
    xmm_pBestSad8x8 = _mm_loadu_si128((__m128i*)p_best_sad8x8);
    xmm_pBestMV8x8 = _mm_loadu_si128((__m128i*)p_best_mv8x8);


    xmm_mv = _mm_cvtsi64_si128(mv);
    xmm_mv = _mm_unpacklo_epi32(xmm_mv, xmm_mv);
    xmm_mv = _mm_unpacklo_epi64(xmm_mv, xmm_mv);

    // sad8x8_0 < p_best_sad8x8[0] for 0 to 3
    sad8x8_less_than_bitmask = _mm_cmplt_epi32(sad8x8_0_3, xmm_pBestSad8x8);

    xmm_N1 = _mm_cmpeq_epi8(xmm_sad8x8_0_1, xmm_sad8x8_0_1);

    sad8x8_greater_or_eq_bitmask = _mm_sub_epi32(xmm_N1, sad8x8_less_than_bitmask);

    BestSad8x8 = _mm_or_si128(_mm_and_si128(xmm_pBestSad8x8, sad8x8_greater_or_eq_bitmask), _mm_and_si128(sad8x8_less_than_bitmask, sad8x8_0_3));
    BestMV8x8 = _mm_or_si128(_mm_and_si128(xmm_pBestMV8x8, sad8x8_greater_or_eq_bitmask), _mm_and_si128(sad8x8_less_than_bitmask, xmm_mv));

    _mm_storeu_si128((__m128i*)p_best_sad8x8, BestSad8x8);
    _mm_storeu_si128((__m128i*)p_best_mv8x8, BestMV8x8);

    xmm_sad16x16 = _mm_add_epi32(xmm_sad8x8_0_1, xmm_sad8x8_2_3);

    xmm_sad16x16_total = _mm_slli_epi32(_mm_add_epi32(_mm_srli_si128(xmm_sad16x16, 8), xmm_sad16x16), 1);

    *p_sad16x16 = sad16x16 = _mm_cvtsi128_si32(xmm_sad16x16_total);

    if (sad16x16 < p_best_sad16x16[0]) {
        p_best_sad16x16[0] = (uint32_t)sad16x16;
        p_best_mv16x16[0] = mv;
    }

}

void ext_sad_calculation_32x32_64x64_sse4_intrin(
    uint32_t  *p_sad16x16,
    uint32_t  *p_best_sad32x32,
    uint32_t  *p_best_sad64x64,
    uint32_t  *p_best_mv32x32,
    uint32_t  *p_best_mv64x64,
    uint32_t   mv,
    uint32_t  *p_sad32x32)
{
    __m128i xmm_N1, sad32x32_greater_than_bitmask, sad32x32_less_than_or_eq_bitmask, BestSad32x32, BestMV32x32, xmm_mv;
    __m128i Sad16x16_0_7_lo, Sad16x16_0_7_hi, Sad16x16_8_15_lo, Sad16x16_8_15_hi, xmm_sad64x64, xmm_sad64x64_total, xmm_pBestSad32x32, xmm_pBestMV32x32;

    Sad16x16_0_7_lo = _mm_unpacklo_epi32(_mm_loadu_si128((__m128i*)p_sad16x16), _mm_loadu_si128((__m128i*)(p_sad16x16 + 4)));
    Sad16x16_0_7_hi = _mm_unpackhi_epi32(_mm_loadu_si128((__m128i*)p_sad16x16), _mm_loadu_si128((__m128i*)(p_sad16x16 + 4)));
    Sad16x16_8_15_lo = _mm_unpacklo_epi32(_mm_loadu_si128((__m128i*)(p_sad16x16 + 8)), _mm_loadu_si128((__m128i*)(p_sad16x16 + 12)));
    Sad16x16_8_15_hi = _mm_unpackhi_epi32(_mm_loadu_si128((__m128i*)(p_sad16x16 + 8)), _mm_loadu_si128((__m128i*)(p_sad16x16 + 12)));

    xmm_sad64x64 = _mm_add_epi32(_mm_add_epi32(_mm_unpacklo_epi64(Sad16x16_0_7_lo, Sad16x16_8_15_lo), _mm_unpackhi_epi64(Sad16x16_0_7_lo, Sad16x16_8_15_lo)),
        _mm_add_epi32(_mm_unpacklo_epi64(Sad16x16_0_7_hi, Sad16x16_8_15_hi), _mm_unpackhi_epi64(Sad16x16_0_7_hi, Sad16x16_8_15_hi)));


    p_sad32x32[0] = _mm_extract_epi32(xmm_sad64x64, 0);
    p_sad32x32[1] = _mm_extract_epi32(xmm_sad64x64, 1);
    p_sad32x32[2] = _mm_extract_epi32(xmm_sad64x64, 2);
    p_sad32x32[3] = _mm_extract_epi32(xmm_sad64x64, 3);

    xmm_sad64x64_total = _mm_add_epi32(_mm_srli_si128(xmm_sad64x64, 8), xmm_sad64x64);

    xmm_sad64x64_total = _mm_add_epi32(_mm_srli_si128(xmm_sad64x64_total, 4), xmm_sad64x64_total);

    xmm_mv = _mm_cvtsi32_si128(mv);
    xmm_mv = _mm_unpacklo_epi32(xmm_mv, xmm_mv);
    xmm_mv = _mm_unpacklo_epi64(xmm_mv, xmm_mv);

    xmm_pBestSad32x32 = _mm_loadu_si128((__m128i*)p_best_sad32x32);
    xmm_pBestMV32x32 = _mm_loadu_si128((__m128i*)p_best_mv32x32);

    sad32x32_greater_than_bitmask = _mm_cmpgt_epi32(xmm_pBestSad32x32, xmm_sad64x64);// _mm_cmplt_epi32(xmm_pBestSad32x32, xmm_sad64x64);

    xmm_N1 = _mm_cmpeq_epi8(xmm_mv, xmm_mv); // anything compared to itself is equal (get 0xFFFFFFFF)
    sad32x32_less_than_or_eq_bitmask = _mm_sub_epi32(xmm_N1, sad32x32_greater_than_bitmask);

    BestSad32x32 = _mm_or_si128(_mm_and_si128(xmm_pBestSad32x32, sad32x32_less_than_or_eq_bitmask), _mm_and_si128(xmm_sad64x64, sad32x32_greater_than_bitmask));
    BestMV32x32 = _mm_or_si128(_mm_and_si128(xmm_pBestMV32x32, sad32x32_less_than_or_eq_bitmask), _mm_and_si128(xmm_mv, sad32x32_greater_than_bitmask));

    _mm_storeu_si128((__m128i*)p_best_sad32x32, BestSad32x32);
    _mm_storeu_si128((__m128i*)p_best_mv32x32, BestMV32x32);


    uint32_t sad64x64 = _mm_cvtsi128_si32(xmm_sad64x64_total);
    if (sad64x64 < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = sad64x64;
        p_best_mv64x64[0] = _mm_cvtsi128_si32(xmm_mv);
    }
}
/*******************************************************************************
 * Requirement: width   = 4, 8, 16, 24, 32, 48 or 64
 * Requirement: height <= 64
 * Requirement: height % 2 = 0 when width = 4 or 8
*******************************************************************************/
void sad_loop_kernel_sse4_1_intrin(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  ref_stride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width,                          // input parameter, block width (N)
    uint64_t *best_sad,
    int16_t *x_search_center,
    int16_t *y_search_center,
    uint32_t  src_stride_raw,                   // input parameter, source stride (no line skipping)
    int16_t search_area_width,
    int16_t search_area_height)
{
    int16_t xBest = *x_search_center, yBest = *y_search_center;
    uint32_t lowSum = 0xffffff;
    uint32_t temSum = 0;
    int16_t i, j;
    uint32_t k, l;
    uint32_t leftover = search_area_width & 7;
    const uint8_t *pRef, *pSrc;
    __m128i s0, s1, s2, s3, s4, s5, s6, s7, s8 = _mm_set1_epi32(-1);

    if (leftover) {
        for (k = 0; k < leftover; k++) {
            s8 = _mm_slli_si128(s8, 2);
        }
    }

    switch (width) {
    case 4:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                pSrc = src;
                pRef = ref + j;
                s3 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)pSrc);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }

            if (leftover) {
                pSrc = src;
                pRef = ref + j;
                s3 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)pSrc);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_or_si128(s3, s8);
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }
            ref += src_stride_raw;
        }
        break;

    case 8:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                pSrc = src;
                pRef = ref + j;
                s3 = s4 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_loadl_epi64((__m128i*)pSrc);
                    s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_adds_epu16(s3, s4);
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }

            if (leftover) {
                pSrc = src;
                pRef = ref + j;
                s3 = s4 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_loadl_epi64((__m128i*)pSrc);
                    s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_adds_epu16(s3, s4);
                s3 = _mm_or_si128(s3, s8);
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }
            ref += src_stride_raw;
        }
        break;

    case 16:
        if (height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s3 = _mm_minpos_epu16(s3);
                    temSum = _mm_extract_epi16(s3, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        yBest = i;
                    }
                }

                if (leftover) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s3 = _mm_or_si128(s3, s8);
                    s3 = _mm_minpos_epu16(s3);
                    temSum = _mm_extract_epi16(s3, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        yBest = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        else if (height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }

                if (leftover) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s3 = _mm_or_si128(s3, s8);
                    s5 = _mm_or_si128(s5, s8);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        if (height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }

                if (leftover) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s3 = _mm_or_si128(s3, s8);
                    s5 = _mm_or_si128(s5, s8);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        if (height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        else {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        if (height <= 32) {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        else {
            __m128i s9, s10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 21 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 21 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s9, 0);
                                    s9 = _mm_srli_si128(s9, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        if (height <= 32) {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        else {
            __m128i s9, s10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 16 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 16 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s9, 0);
                                    s9 = _mm_srli_si128(s9, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        assert(0);
        break;
    }

    *best_sad = lowSum;
    *x_search_center = xBest;
    *y_search_center = yBest;
}

void sad_loop_kernel_sparse_sse4_1_intrin(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  ref_stride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width,                          // input parameter, block width (N)
    uint64_t *best_sad,
    int16_t *x_search_center,
    int16_t *y_search_center,
    uint32_t  src_stride_raw,                   // input parameter, source stride (no line skipping)
    int16_t search_area_width,
    int16_t search_area_height)
{
    int16_t xBest = *x_search_center, yBest = *y_search_center;
    uint32_t lowSum = 0xffffff;
    uint32_t temSum = 0;
    int16_t i, j;
    uint32_t k, l;
    uint32_t leftover = search_area_width & 7;
    const uint8_t *pRef, *pSrc;
    __m128i s0, s1, s2, s3, s4, s5, s6, s7, s8 = _mm_set1_epi32(-1);

    if (leftover) {
        for (k = 0; k < leftover; k++) {
            s8 = _mm_slli_si128(s8, 2);
        }
    }

    switch (width) {
    case 4:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                pSrc = src;
                pRef = ref + j;
                s3 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)pSrc);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }

            if (leftover) {
                pSrc = src;
                pRef = ref + j;
                s3 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)pSrc);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_or_si128(s3, s8);
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }
            ref += src_stride_raw;
        }
        break;

    case 8:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                pSrc = src;
                pRef = ref + j;
                s3 = s4 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_loadl_epi64((__m128i*)pSrc);
                    s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_adds_epu16(s3, s4);
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }

            if (leftover) {
                pSrc = src;
                pRef = ref + j;
                s3 = s4 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_loadl_epi64((__m128i*)pSrc);
                    s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_adds_epu16(s3, s4);
                s3 = _mm_or_si128(s3, s8);
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }
            ref += src_stride_raw;
        }
        break;

    case 16:
        if (height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                uint32_t startW = (i & 1) << 3;
                for (j = startW; j <= search_area_width - 8; j += 16) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s3 = _mm_minpos_epu16(s3);
                    temSum = _mm_extract_epi16(s3, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        yBest = i;
                    }
                }

                if (leftover) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s3 = _mm_or_si128(s3, s8);
                    s3 = _mm_minpos_epu16(s3);
                    temSum = _mm_extract_epi16(s3, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        yBest = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        else if (height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }

                if (leftover) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s3 = _mm_or_si128(s3, s8);
                    s5 = _mm_or_si128(s5, s8);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        if (height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }

                if (leftover) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s3 = _mm_or_si128(s3, s8);
                    s5 = _mm_or_si128(s5, s8);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        if (height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        else {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        if (height <= 32) {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        else {
            __m128i s9, s10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 21 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 21 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s9, 0);
                                    s9 = _mm_srli_si128(s9, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        if (height <= 32) {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s0, 0);
                                    s0 = _mm_srli_si128(s0, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        else {
            __m128i s9, s10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 16 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 16 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_or_si128(s0, s8);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
                            k = leftover;
                            while (k > 0) {
                                for (l = 0; l < 4 && k; l++, k--) {
                                    temSum = _mm_extract_epi32(s9, 0);
                                    s9 = _mm_srli_si128(s9, 4);
                                    if (temSum < lowSum) {
                                        lowSum = temSum;
                                        xBest = (int16_t)(j + leftover - k);
                                        yBest = i;
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
        assert(0);
        break;
    }

    *best_sad = lowSum;
    *x_search_center = xBest;
    *y_search_center = yBest;
}

/*******************************************************************************
* Requirement: width   = 4, 8, 16, 24, 32, 48 or 64
* Requirement: height <= 64
* Requirement: height % 2 = 0 when width = 4 or 8
*******************************************************************************/
void sad_loop_kernel_sse4_1_hme_l0_intrin(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  ref_stride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width,                          // input parameter, block width (N)
    uint64_t *best_sad,
    int16_t *x_search_center,
    int16_t *y_search_center,
    uint32_t  src_stride_raw,                   // input parameter, source stride (no line skipping)
    int16_t search_area_width,
    int16_t search_area_height)
{
    int16_t xBest = *x_search_center, yBest = *y_search_center;
    uint32_t lowSum = 0xffffff;
    uint32_t temSum = 0;
    int16_t i, j;
    uint32_t k, l;
    const uint8_t *pRef, *pSrc;
    __m128i s0, s1, s2, s3, s4, s5, s6, s7, s9, s10, s11;


    switch (width) {
    case 4:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                pSrc = src;
                pRef = ref + j;
                s3 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_cvtsi32_si128(*(uint32_t *)pSrc);
                    s5 = _mm_cvtsi32_si128(*(uint32_t *)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }
            ref += src_stride_raw;
        }
        break;

    case 8:
        for (i = 0; i < search_area_height; i++) {
            for (j = 0; j <= search_area_width - 8; j += 8) {
                pSrc = src;
                pRef = ref + j;
                s3 = s4 = _mm_setzero_si128();
                for (k = 0; k < height; k += 2) {
                    s0 = _mm_loadu_si128((__m128i*)pRef);
                    s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride));
                    s2 = _mm_loadl_epi64((__m128i*)pSrc);
                    s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                    s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
                    s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));
                    pSrc += src_stride << 1;
                    pRef += ref_stride << 1;
                }
                s3 = _mm_adds_epu16(s3, s4);
                s3 = _mm_minpos_epu16(s3);
                temSum = _mm_extract_epi16(s3, 0);
                if (temSum < lowSum) {
                    lowSum = temSum;
                    xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                    yBest = i;
                }
            }

            ref += src_stride_raw;
        }
        break;

    case 16:
        if (height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 16; j += 16) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    s7 = s9 = s10 = s11 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s7 = _mm_adds_epu16(s7, _mm_mpsadbw_epu8(s1, s2, 0));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 5));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 2));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s3 = _mm_minpos_epu16(s3);
                    temSum = _mm_extract_epi16(s3, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s3, 1));
                        yBest = i;
                    }

                    s7 = _mm_adds_epu16(_mm_adds_epu16(s7, s11), _mm_adds_epu16(s9, s10));
                    s7 = _mm_minpos_epu16(s7);
                    temSum = _mm_extract_epi16(s7, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + 8 + _mm_extract_epi16(s7, 1));
                        yBest = i;
                    }
                }

                ref += src_stride_raw;
            }
        }
        else if (height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }

                ref += src_stride_raw;
            }
        }
        else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                ref += src_stride_raw;
            }
        }
        break;

    case 24:
        if (height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }
                ref += src_stride_raw;
            }
        }
        else {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s2 = _mm_loadl_epi64((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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

                ref += src_stride_raw;
            }
        }
        break;

    case 32:
        if (height <= 16) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s3 = _mm_adds_epu16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s6);
                    s4 = _mm_minpos_epu16(s3);
                    s6 = _mm_minpos_epu16(s5);
                    s4 = _mm_unpacklo_epi16(s4, s4);
                    s4 = _mm_unpacklo_epi32(s4, s4);
                    s4 = _mm_unpacklo_epi64(s4, s4);
                    s6 = _mm_unpacklo_epi16(s6, s6);
                    s6 = _mm_unpacklo_epi32(s6, s6);
                    s6 = _mm_unpacklo_epi64(s6, s6);
                    s3 = _mm_sub_epi16(s3, s4);
                    s5 = _mm_adds_epu16(s5, s3);
                    s5 = _mm_sub_epi16(s5, s6);
                    s5 = _mm_minpos_epu16(s5);
                    temSum = _mm_extract_epi16(s5, 0);
                    temSum += _mm_extract_epi16(s4, 0);
                    temSum += _mm_extract_epi16(s6, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        lowSum = temSum;
                        xBest = (int16_t)(j + _mm_extract_epi16(s5, 1));
                        yBest = i;
                    }
                }

                ref += src_stride_raw;
            }
        }
        else if (height <= 32) {
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                ref += src_stride_raw;
            }
        }
        else {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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

                ref += src_stride_raw;
            }
        }
        break;

    case 48:
        if (height <= 32) {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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

                ref += src_stride_raw;
            }
        }
        else {
            __m128i s9, s10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 21 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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

                ref += src_stride_raw;
            }
        }
        break;

    case 64:
        if (height <= 32) {
            __m128i s9, s10, s11, s12;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s3 = s4 = s5 = s6 = _mm_setzero_si128();
                    for (k = 0; k < height >> 1; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                        s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                        s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s9 = s10 = s11 = s12 = _mm_setzero_si128();
                    for (; k < height; k++) {
                        s0 = _mm_loadu_si128((__m128i*)pRef);
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                        s2 = _mm_loadu_si128((__m128i*)pSrc);
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                        s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                        s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                        s9 = _mm_adds_epu16(s9, _mm_mpsadbw_epu8(s0, s2, 0));
                        s10 = _mm_adds_epu16(s10, _mm_mpsadbw_epu8(s0, s2, 5));
                        s11 = _mm_adds_epu16(s11, _mm_mpsadbw_epu8(s1, s2, 2));
                        s12 = _mm_adds_epu16(s12, _mm_mpsadbw_epu8(s1, s2, 7));
                        pSrc += src_stride;
                        pRef += ref_stride;
                    }
                    s0 = _mm_adds_epu16(_mm_adds_epu16(s3, s4), _mm_adds_epu16(s5, s6));
                    s0 = _mm_adds_epu16(s0, _mm_adds_epu16(_mm_adds_epu16(s9, s10), _mm_adds_epu16(s11, s12)));
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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
                            s1 = _mm_unpacklo_epi16(s9, _mm_setzero_si128());
                            s9 = _mm_unpackhi_epi16(s9, _mm_setzero_si128());
                            s2 = _mm_unpacklo_epi16(s10, _mm_setzero_si128());
                            s10 = _mm_unpackhi_epi16(s10, _mm_setzero_si128());
                            s4 = _mm_unpacklo_epi16(s11, _mm_setzero_si128());
                            s11 = _mm_unpackhi_epi16(s11, _mm_setzero_si128());
                            s5 = _mm_unpacklo_epi16(s12, _mm_setzero_si128());
                            s12 = _mm_unpackhi_epi16(s12, _mm_setzero_si128());
                            s0 = _mm_add_epi32(s0, _mm_add_epi32(_mm_add_epi32(s1, s2), _mm_add_epi32(s4, s5)));
                            s3 = _mm_add_epi32(s3, _mm_add_epi32(_mm_add_epi32(s9, s10), _mm_add_epi32(s11, s12)));
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

                ref += src_stride_raw;
            }
        }
        else {
            __m128i s9, s10;
            for (i = 0; i < search_area_height; i++) {
                for (j = 0; j <= search_area_width - 8; j += 8) {
                    pSrc = src;
                    pRef = ref + j;
                    s9 = s10 = _mm_setzero_si128();
                    k = 0;
                    while (k < height) {
                        s3 = s4 = s5 = s6 = _mm_setzero_si128();
                        for (l = 0; l < 16 && k < height; k++, l++) {
                            s0 = _mm_loadu_si128((__m128i*)pRef);
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 8));
                            s2 = _mm_loadu_si128((__m128i*)pSrc);
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 16));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 24));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 16));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 32));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 40));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 32));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            s0 = _mm_loadu_si128((__m128i*)(pRef + 48));
                            s1 = _mm_loadu_si128((__m128i*)(pRef + 56));
                            s2 = _mm_loadu_si128((__m128i*)(pSrc + 48));
                            s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
                            s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
                            s5 = _mm_adds_epu16(s5, _mm_mpsadbw_epu8(s1, s2, 2));
                            s6 = _mm_adds_epu16(s6, _mm_mpsadbw_epu8(s1, s2, 7));
                            pSrc += src_stride;
                            pRef += ref_stride;
                        }
                        s0 = _mm_unpacklo_epi16(s3, _mm_setzero_si128());
                        s3 = _mm_unpackhi_epi16(s3, _mm_setzero_si128());
                        s1 = _mm_unpacklo_epi16(s4, _mm_setzero_si128());
                        s4 = _mm_unpackhi_epi16(s4, _mm_setzero_si128());
                        s2 = _mm_unpacklo_epi16(s5, _mm_setzero_si128());
                        s5 = _mm_unpackhi_epi16(s5, _mm_setzero_si128());
                        s7 = _mm_unpacklo_epi16(s6, _mm_setzero_si128());
                        s6 = _mm_unpackhi_epi16(s6, _mm_setzero_si128());
                        s9 = _mm_add_epi32(s9, _mm_add_epi32(_mm_add_epi32(s0, s1), _mm_add_epi32(s2, s7)));
                        s10 = _mm_add_epi32(s10, _mm_add_epi32(_mm_add_epi32(s3, s4), _mm_add_epi32(s5, s6)));
                    }
                    s0 = _mm_packus_epi32(s9, s10);
                    s0 = _mm_minpos_epu16(s0);
                    temSum = _mm_extract_epi16(s0, 0);
                    temSum &= 0x0000FFFF;
                    if (temSum < lowSum) {
                        if (temSum != 0xFFFF) { // no overflow
                            lowSum = temSum;
                            xBest = (int16_t)(j + _mm_extract_epi16(s0, 1));
                            yBest = i;
                        }
                        else {
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

                ref += src_stride_raw;
            }
        }
        break;

    default:
        assert(0);
        break;
    }

    *best_sad = lowSum;
    *x_search_center = xBest;
    *y_search_center = yBest;
}

/*******************************************
 * GetEightHorizontalSearchPointResults_8x8_16x16_PU
 *******************************************/
void get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin(
    uint8_t   *src,
    uint32_t   src_stride,
    uint8_t   *ref,
    uint32_t   ref_stride,
    uint32_t  *p_best_sad8x8,
    uint32_t  *p_best_mv8x8,
    uint32_t  *p_best_sad16x16,
    uint32_t  *p_best_mv16x16,
    uint32_t   mv,
    uint16_t  *p_sad16x16)
{

    int16_t x_mv, y_mv;
    const uint8_t *pRef, *pSrc;
    __m128i s0, s1, s2, s3, s4, s5;
    __m128i sad_0, sad_1, sad_2, sad_3;
    uint32_t temSum;


    /*
   -------------------------------------   -----------------------------------
   | 8x8_00 | 8x8_01 | 8x8_04 | 8x8_05 |   8x8_16 | 8x8_17 | 8x8_20 | 8x8_21 |
   -------------------------------------   -----------------------------------
   | 8x8_02 | 8x8_03 | 8x8_06 | 8x8_07 |   8x8_18 | 8x8_19 | 8x8_22 | 8x8_23 |
   -----------------------   -----------   ----------------------   ----------
   | 8x8_08 | 8x8_09 | 8x8_12 | 8x8_13 |   8x8_24 | 8x8_25 | 8x8_29 | 8x8_29 |
   ----------------------    -----------   ---------------------    ----------
   | 8x8_10 | 8x8_11 | 8x8_14 | 8x8_15 |   8x8_26 | 8x8_27 | 8x8_30 | 8x8_31 |
   -------------------------------------   -----------------------------------

   -------------------------------------   -----------------------------------
   | 8x8_32 | 8x8_33 | 8x8_36 | 8x8_37 |   8x8_48 | 8x8_49 | 8x8_52 | 8x8_53 |
   -------------------------------------   -----------------------------------
   | 8x8_34 | 8x8_35 | 8x8_38 | 8x8_39 |   8x8_50 | 8x8_51 | 8x8_54 | 8x8_55 |
   -----------------------   -----------   ----------------------   ----------
   | 8x8_40 | 8x8_41 | 8x8_44 | 8x8_45 |   8x8_56 | 8x8_57 | 8x8_60 | 8x8_61 |
   ----------------------    -----------   ---------------------    ----------
   | 8x8_42 | 8x8_43 | 8x8_46 | 8x8_48 |   8x8_58 | 8x8_59 | 8x8_62 | 8x8_63 |
   -------------------------------------   -----------------------------------
   */

   /*
   ----------------------    ----------------------
   |  16x16_0  |  16x16_1  |  16x16_4  |  16x16_5  |
   ----------------------    ----------------------
   |  16x16_2  |  16x16_3  |  16x16_6  |  16x16_7  |
   -----------------------   -----------------------
   |  16x16_8  |  16x16_9  |  16x16_12 |  16x16_13 |
   ----------------------    ----------------------
   |  16x16_10 |  16x16_11 |  16x16_14 |  16x16_15 |
   -----------------------   -----------------------
   */


   //8x8_0
    {
        pSrc = src;
        pRef = ref;
        s3 = s4 = _mm_setzero_si128();

        s0 = _mm_loadu_si128((__m128i*)pRef);
        s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride * 2));
        s2 = _mm_loadl_epi64((__m128i*)pSrc);
        s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride * 2));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));

        pSrc += src_stride * 4;
        pRef += ref_stride * 4;

        s0 = _mm_loadu_si128((__m128i*)pRef);
        s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride * 2));
        s2 = _mm_loadl_epi64((__m128i*)pSrc);
        s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride * 2));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));

        //final 8x4 SAD
        sad_0 = _mm_adds_epu16(s3, s4);

        //find the best for 8x8_0
        s3 = _mm_minpos_epu16(sad_0);
        temSum = _mm_extract_epi16(s3, 0);
        if (2 * temSum < p_best_sad8x8[0]) {
            p_best_sad8x8[0] = 2 * temSum;
            x_mv = _MVXT(mv) + (int16_t)(_mm_extract_epi16(s3, 1) * 4);
            y_mv = _MVYT(mv);
            p_best_mv8x8[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }

    //8x8_1
    {
        pSrc = src + 8;
        pRef = ref + 8;
        s3 = s4 = _mm_setzero_si128();

        s0 = _mm_loadu_si128((__m128i*)pRef);
        s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride * 2));
        s2 = _mm_loadl_epi64((__m128i*)pSrc);
        s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride * 2));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));

        pSrc += src_stride * 4;
        pRef += ref_stride * 4;

        s0 = _mm_loadu_si128((__m128i*)pRef);
        s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride * 2));
        s2 = _mm_loadl_epi64((__m128i*)pSrc);
        s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride * 2));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));

        //final 8x4 SAD
        sad_1 = _mm_adds_epu16(s3, s4);

        //find the best for 8x8_1
        s3 = _mm_minpos_epu16(sad_1);
        temSum = _mm_extract_epi16(s3, 0);
        if (2 * temSum < p_best_sad8x8[1]) {
            p_best_sad8x8[1] = 2 * temSum;
            x_mv = _MVXT(mv) + (int16_t)(_mm_extract_epi16(s3, 1) * 4);
            y_mv = _MVYT(mv);
            p_best_mv8x8[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }

    //8x8_2
    {
        pSrc = src + 8 * src_stride;
        pRef = ref + 8 * ref_stride;
        s3 = s4 = _mm_setzero_si128();

        s0 = _mm_loadu_si128((__m128i*)pRef);
        s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride * 2));
        s2 = _mm_loadl_epi64((__m128i*)pSrc);
        s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride * 2));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));

        pSrc += src_stride * 4;
        pRef += ref_stride * 4;

        s0 = _mm_loadu_si128((__m128i*)pRef);
        s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride * 2));
        s2 = _mm_loadl_epi64((__m128i*)pSrc);
        s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride * 2));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));

        //final 8x4 SAD
        sad_2 = _mm_adds_epu16(s3, s4);

        //find the best for 8x8_2
        s3 = _mm_minpos_epu16(sad_2);
        temSum = _mm_extract_epi16(s3, 0);
        if (2 * temSum < p_best_sad8x8[2]) {
            p_best_sad8x8[2] = 2 * temSum;
            x_mv = _MVXT(mv) + (int16_t)(_mm_extract_epi16(s3, 1) * 4);
            y_mv = _MVYT(mv);
            p_best_mv8x8[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }

    //8x8_3
    {
        pSrc = src + 8 + 8 * src_stride;
        pRef = ref + 8 + 8 * ref_stride;
        s3 = s4 = _mm_setzero_si128();

        s0 = _mm_loadu_si128((__m128i*)pRef);
        s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride * 2));
        s2 = _mm_loadl_epi64((__m128i*)pSrc);
        s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride * 2));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));

        pSrc += src_stride * 4;
        pRef += ref_stride * 4;

        s0 = _mm_loadu_si128((__m128i*)pRef);
        s1 = _mm_loadu_si128((__m128i*)(pRef + ref_stride * 2));
        s2 = _mm_loadl_epi64((__m128i*)pSrc);
        s5 = _mm_loadl_epi64((__m128i*)(pSrc + src_stride * 2));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s0, s2, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s0, s2, 5));
        s3 = _mm_adds_epu16(s3, _mm_mpsadbw_epu8(s1, s5, 0));
        s4 = _mm_adds_epu16(s4, _mm_mpsadbw_epu8(s1, s5, 5));

        //final 8x4 SAD
        sad_3 = _mm_adds_epu16(s3, s4);

        //find the best for 8x8_3
        s3 = _mm_minpos_epu16(sad_3);
        temSum = _mm_extract_epi16(s3, 0);
        if (2 * temSum < p_best_sad8x8[3]) {
            p_best_sad8x8[3] = 2 * temSum;
            x_mv = _MVXT(mv) + (int16_t)(_mm_extract_epi16(s3, 1) * 4);
            y_mv = _MVYT(mv);
            p_best_mv8x8[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }

    //16x16
    {
        s0 = _mm_adds_epu16(sad_0, sad_1);
        s1 = _mm_adds_epu16(sad_2, sad_3);
        s3 = _mm_adds_epu16(s0, s1);
        //sotore the 8 SADs(16x8 SADs)
        _mm_store_si128((__m128i*)p_sad16x16, s3);
        //find the best for 16x16
        s3 = _mm_minpos_epu16(s3);
        temSum = _mm_extract_epi16(s3, 0);
        if (2 * temSum < p_best_sad16x16[0]) {
            p_best_sad16x16[0] = 2 * temSum;
            x_mv = _MVXT(mv) + (int16_t)(_mm_extract_epi16(s3, 1) * 4);
            y_mv = _MVYT(mv);
            p_best_mv16x16[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }


}

/*******************************************
 Calcualte SAD for 32x32,64x64 from 16x16
 and check if there is improvement, if yes keep
 the best SAD+MV
 *******************************************/
void get_eight_horizontal_search_point_results_32x32_64x64_pu_sse41_intrin(
    uint16_t  *p_sad16x16,
    uint32_t  *p_best_sad32x32,
    uint32_t  *p_best_sad64x64,
    uint32_t  *p_best_mv32x32,
    uint32_t  *p_best_mv64x64,
    uint32_t   mv)
{
    int16_t x_mv, y_mv;

    uint32_t temSum;
    __m128i s0, s1, s2, s3, s4, s5, sad_0, sad_1, s6, s7;
    __m128i sad_00, sad_01, sad_10, sad_11, sad_20, sad_21, sad_30, sad_31;
    __m128i Zero = _mm_setzero_si128();

    /*--------------------
    |  32x32_0  |  32x32_1
    ----------------------
    |  32x32_2  |  32x32_3
    ----------------------*/


    /*  data ordering in p_sad16x16 buffer

                  Search    Search            Search
                  Point 0   Point 1           Point 7
                ---------------------------------------
     16x16_0    |    x    |    x    | ...... |    x    |
                ---------------------------------------
     16x16_1    |    x    |    x    | ...... |    x    |

     16x16_n    |    x    |    x    | ...... |    x    |

                ---------------------------------------
     16x16_15   |    x    |    x    | ...... |    x    |
                ---------------------------------------
    */


    //32x32_0
    s0 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 0 * 8));
    s1 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 1 * 8));
    s2 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 2 * 8));
    s3 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 3 * 8));

    s4 = _mm_unpackhi_epi16(s0, Zero);
    s5 = _mm_unpacklo_epi16(s0, Zero);
    s6 = _mm_unpackhi_epi16(s1, Zero);
    s7 = _mm_unpacklo_epi16(s1, Zero);
    s0 = _mm_add_epi32(s4, s6);
    s1 = _mm_add_epi32(s5, s7);

    s4 = _mm_unpackhi_epi16(s2, Zero);
    s5 = _mm_unpacklo_epi16(s2, Zero);
    s6 = _mm_unpackhi_epi16(s3, Zero);
    s7 = _mm_unpacklo_epi16(s3, Zero);
    s2 = _mm_add_epi32(s4, s6);
    s3 = _mm_add_epi32(s5, s7);

    sad_01 = _mm_add_epi32(s0, s2);
    sad_00 = _mm_add_epi32(s1, s3);

    //sad_00
    temSum = _mm_extract_epi32(sad_00, 0);
    if (2 * temSum < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_00, 1);
    if (2 * temSum < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_00, 2);
    if (2 * temSum < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_00, 3);
    if (2 * temSum < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }

    //sad_01
    temSum = _mm_extract_epi32(sad_01, 0);
    if (2 * temSum < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_01, 1);
    if (2 * temSum < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_01, 2);
    if (2 * temSum < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_01, 3);
    if (2 * temSum < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }

    //32x32_1
    s0 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 4 * 8));
    s1 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 5 * 8));
    s2 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 6 * 8));
    s3 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 7 * 8));

    s4 = _mm_unpackhi_epi16(s0, Zero);
    s5 = _mm_unpacklo_epi16(s0, Zero);
    s6 = _mm_unpackhi_epi16(s1, Zero);
    s7 = _mm_unpacklo_epi16(s1, Zero);
    s0 = _mm_add_epi32(s4, s6);
    s1 = _mm_add_epi32(s5, s7);

    s4 = _mm_unpackhi_epi16(s2, Zero);
    s5 = _mm_unpacklo_epi16(s2, Zero);
    s6 = _mm_unpackhi_epi16(s3, Zero);
    s7 = _mm_unpacklo_epi16(s3, Zero);
    s2 = _mm_add_epi32(s4, s6);
    s3 = _mm_add_epi32(s5, s7);

    sad_11 = _mm_add_epi32(s0, s2);
    sad_10 = _mm_add_epi32(s1, s3);

    //sad_10
    temSum = _mm_extract_epi32(sad_10, 0);
    if (2 * temSum < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_10, 1);
    if (2 * temSum < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_10, 2);
    if (2 * temSum < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_10, 3);
    if (2 * temSum < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }

    //sad_11
    temSum = _mm_extract_epi32(sad_11, 0);
    if (2 * temSum < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_11, 1);
    if (2 * temSum < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_11, 2);
    if (2 * temSum < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_11, 3);
    if (2 * temSum < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }



    //32x32_2
    s0 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 8 * 8));
    s1 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 9 * 8));
    s2 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 10 * 8));
    s3 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 11 * 8));

    s4 = _mm_unpackhi_epi16(s0, Zero);
    s5 = _mm_unpacklo_epi16(s0, Zero);
    s6 = _mm_unpackhi_epi16(s1, Zero);
    s7 = _mm_unpacklo_epi16(s1, Zero);
    s0 = _mm_add_epi32(s4, s6);
    s1 = _mm_add_epi32(s5, s7);

    s4 = _mm_unpackhi_epi16(s2, Zero);
    s5 = _mm_unpacklo_epi16(s2, Zero);
    s6 = _mm_unpackhi_epi16(s3, Zero);
    s7 = _mm_unpacklo_epi16(s3, Zero);
    s2 = _mm_add_epi32(s4, s6);
    s3 = _mm_add_epi32(s5, s7);

    sad_21 = _mm_add_epi32(s0, s2);
    sad_20 = _mm_add_epi32(s1, s3);

    //sad_20
    temSum = _mm_extract_epi32(sad_20, 0);
    if (2 * temSum < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_20, 1);
    if (2 * temSum < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_20, 2);
    if (2 * temSum < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_20, 3);
    if (2 * temSum < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }

    //sad_21
    temSum = _mm_extract_epi32(sad_21, 0);
    if (2 * temSum < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_21, 1);
    if (2 * temSum < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_21, 2);
    if (2 * temSum < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_21, 3);
    if (2 * temSum < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }



    //32x32_3
    s0 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 12 * 8));
    s1 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 13 * 8));
    s2 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 14 * 8));
    s3 = _mm_loadu_si128((__m128i*)(p_sad16x16 + 15 * 8));

    s4 = _mm_unpackhi_epi16(s0, Zero);
    s5 = _mm_unpacklo_epi16(s0, Zero);
    s6 = _mm_unpackhi_epi16(s1, Zero);
    s7 = _mm_unpacklo_epi16(s1, Zero);
    s0 = _mm_add_epi32(s4, s6);
    s1 = _mm_add_epi32(s5, s7);

    s4 = _mm_unpackhi_epi16(s2, Zero);
    s5 = _mm_unpacklo_epi16(s2, Zero);
    s6 = _mm_unpackhi_epi16(s3, Zero);
    s7 = _mm_unpacklo_epi16(s3, Zero);
    s2 = _mm_add_epi32(s4, s6);
    s3 = _mm_add_epi32(s5, s7);

    sad_31 = _mm_add_epi32(s0, s2);
    sad_30 = _mm_add_epi32(s1, s3);

    //sad_30
    temSum = _mm_extract_epi32(sad_30, 0);
    if (2 * temSum < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_30, 1);
    if (2 * temSum < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_30, 2);
    if (2 * temSum < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_30, 3);
    if (2 * temSum < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }

    //sad_31
    temSum = _mm_extract_epi32(sad_31, 0);
    if (2 * temSum < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_31, 1);
    if (2 * temSum < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_31, 2);
    if (2 * temSum < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_31, 3);
    if (2 * temSum < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }


    sad_0 = _mm_add_epi32(_mm_add_epi32(sad_00, sad_10), _mm_add_epi32(sad_20, sad_30));
    sad_1 = _mm_add_epi32(_mm_add_epi32(sad_01, sad_11), _mm_add_epi32(sad_21, sad_31));

    //sad_0
    temSum = _mm_extract_epi32(sad_0, 0);
    if (2 * temSum < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_0, 1);
    if (2 * temSum < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_0, 2);
    if (2 * temSum < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_0, 3);
    if (2 * temSum < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (0 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }

    //sad_1
    temSum = _mm_extract_epi32(sad_1, 0);
    if (2 * temSum < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 0) * 4;   y_mv = _MVYT(mv);
        p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_1, 1);
    if (2 * temSum < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 1) * 4;  y_mv = _MVYT(mv);
        p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_1, 2);
    if (2 * temSum < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 2) * 4;  y_mv = _MVYT(mv);
        p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }
    temSum = _mm_extract_epi32(sad_1, 3);
    if (2 * temSum < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = 2 * temSum;
        x_mv = _MVXT(mv) + (4 + 3) * 4;  y_mv = _MVYT(mv);
        p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
    }

}
