/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbCombinedAveragingSAD_Inline_AVX2_h
#define EbCombinedAveragingSAD_Inline_AVX2_h

#include "immintrin.h"
#include "EbDefinitions.h"
#include "EbMemory_AVX2.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE void ssd8x2_avx2(const uint8_t *const src, const ptrdiff_t src_stride,
                               const uint8_t *const ref1, const ptrdiff_t ref1_stride,
                               const uint8_t *const ref2, const ptrdiff_t ref2_stride,
                               __m256i *const sum) {
    const __m256i zero    = _mm256_setzero_si256();
    const __m256i s       = load_u8_8x2_avx2(src, src_stride);
    const __m256i r1      = load_u8_8x2_avx2(ref1, ref1_stride);
    const __m256i r2      = load_u8_8x2_avx2(ref2, ref2_stride);
    const __m256i avg     = _mm256_avg_epu8(r1, r2);
    const __m256i s_256   = _mm256_unpacklo_epi8(s, zero);
    const __m256i avg_256 = _mm256_unpacklo_epi8(avg, zero);
    const __m256i dif     = _mm256_sub_epi16(s_256, avg_256);
    const __m256i sqr     = _mm256_madd_epi16(dif, dif);
    *sum                  = _mm256_add_epi32(*sum, sqr);
}

static INLINE void ssd32_avx2(const uint8_t *const src, const uint8_t *const ref1,
                              const uint8_t *const ref2, __m256i *const sum) {
    const __m256i zero = _mm256_setzero_si256();
    const __m256i s    = _mm256_loadu_si256((__m256i *)src);
    const __m256i r1   = _mm256_loadu_si256((__m256i *)ref1);
    const __m256i r2   = _mm256_loadu_si256((__m256i *)ref2);
    const __m256i avg  = _mm256_avg_epu8(r1, r2);
    const __m256i s0   = _mm256_unpacklo_epi8(s, zero);
    const __m256i s1   = _mm256_unpackhi_epi8(s, zero);
    const __m256i avg0 = _mm256_unpacklo_epi8(avg, zero);
    const __m256i avg1 = _mm256_unpackhi_epi8(avg, zero);
    const __m256i dif0 = _mm256_sub_epi16(s0, avg0);
    const __m256i dif1 = _mm256_sub_epi16(s1, avg1);
    const __m256i sqr0 = _mm256_madd_epi16(dif0, dif0);
    const __m256i sqr1 = _mm256_madd_epi16(dif1, dif1);
    *sum               = _mm256_add_epi32(*sum, sqr0);
    *sum               = _mm256_add_epi32(*sum, sqr1);
}

#ifdef __cplusplus
}
#endif
#endif // EbCombinedAveragingSAD_Inline_AVX2_h
