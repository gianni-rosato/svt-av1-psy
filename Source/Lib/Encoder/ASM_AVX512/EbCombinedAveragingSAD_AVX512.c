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

#include "EbDefinitions.h"

#ifndef NON_AVX512_SUPPORT

#include <immintrin.h>
#include "EbCombinedAveragingSAD_Inline_AVX2.h"
#include "EbMemory_AVX2.h"
#include "EbMemory_SSE4_1.h"

static INLINE void ssd32_avx512(const uint8_t *const src, const uint8_t *const ref1,
                                const uint8_t *const ref2, __m512i *const sum) {
    const __m256i s    = _mm256_loadu_si256((__m256i *)src);
    const __m256i r1   = _mm256_loadu_si256((__m256i *)ref1);
    const __m256i r2   = _mm256_loadu_si256((__m256i *)ref2);
    const __m256i avg  = _mm256_avg_epu8(r1, r2);
    const __m512i s0   = _mm512_cvtepu8_epi16(s);
    const __m512i avg0 = _mm512_cvtepu8_epi16(avg);
    const __m512i dif0 = _mm512_sub_epi16(s0, avg0);
    const __m512i sqr0 = _mm512_madd_epi16(dif0, dif0);
    *sum               = _mm512_add_epi32(*sum, sqr0);
}

static INLINE void ssd64_avx512(const uint8_t *const src, const uint8_t *const ref1,
                                const uint8_t *const ref2, __m512i *const sum) {
    const __m512i zero = _mm512_setzero_si512();
    const __m512i s    = _mm512_loadu_si512((__m512i *)src);
    const __m512i r1   = _mm512_loadu_si512((__m512i *)ref1);
    const __m512i r2   = _mm512_loadu_si512((__m512i *)ref2);
    const __m512i avg  = _mm512_avg_epu8(r1, r2);
    const __m512i s0   = _mm512_unpacklo_epi8(s, zero);
    const __m512i s1   = _mm512_unpackhi_epi8(s, zero);
    const __m512i avg0 = _mm512_unpacklo_epi8(avg, zero);
    const __m512i avg1 = _mm512_unpackhi_epi8(avg, zero);
    const __m512i dif0 = _mm512_sub_epi16(s0, avg0);
    const __m512i dif1 = _mm512_sub_epi16(s1, avg1);
    const __m512i sqr0 = _mm512_madd_epi16(dif0, dif0);
    const __m512i sqr1 = _mm512_madd_epi16(dif1, dif1);
    *sum               = _mm512_add_epi32(*sum, sqr0);
    *sum               = _mm512_add_epi32(*sum, sqr1);
}

uint32_t combined_averaging_ssd_avx512(uint8_t *src, ptrdiff_t src_stride, uint8_t *ref1,
                                       ptrdiff_t ref1_stride, uint8_t *ref2, ptrdiff_t ref2_stride,
                                       uint32_t height, uint32_t width) {
    uint32_t y = height;
    __m128i  sum_128;

    if (width & 4) {
        const __m128i zero = _mm_setzero_si128();

        sum_128 = _mm_setzero_si128();

        do {
            uint32_t x = 0;
            do {
                const __m128i s     = load_u8_4x2_sse4_1(src + x, src_stride);
                const __m128i r1    = load_u8_4x2_sse4_1(ref1 + x, ref1_stride);
                const __m128i r2    = load_u8_4x2_sse4_1(ref2 + x, ref2_stride);
                const __m128i avg   = _mm_avg_epu8(r1, r2);
                const __m128i s16   = _mm_unpacklo_epi8(s, zero);
                const __m128i avg16 = _mm_unpacklo_epi8(avg, zero);
                const __m128i dif   = _mm_sub_epi16(s16, avg16);
                const __m128i sqr   = _mm_madd_epi16(dif, dif);
                sum_128             = _mm_add_epi32(sum_128, sqr);
                x += 4;
            } while (x < width);

            src += 2 * src_stride;
            ref1 += 2 * ref1_stride;
            ref2 += 2 * ref2_stride;
            y -= 2;
        } while (y);
    } else {
        __m256i sum;

        if (width == 8) {
            sum = _mm256_setzero_si256();

            do {
                ssd8x2_avx2(src, src_stride, ref1, ref1_stride, ref2, ref2_stride, &sum);
                src += 2 * src_stride;
                ref1 += 2 * ref1_stride;
                ref2 += 2 * ref2_stride;
                y -= 2;
            } while (y);
        } else if (width == 16) {
            sum = _mm256_setzero_si256();

            do {
                const __m128i s       = _mm_loadu_si128((__m128i *)src);
                const __m128i r1      = _mm_loadu_si128((__m128i *)ref1);
                const __m128i r2      = _mm_loadu_si128((__m128i *)ref2);
                const __m128i avg     = _mm_avg_epu8(r1, r2);
                const __m256i s_256   = _mm256_cvtepu8_epi16(s);
                const __m256i avg_256 = _mm256_cvtepu8_epi16(avg);
                const __m256i dif     = _mm256_sub_epi16(s_256, avg_256);
                const __m256i sqr     = _mm256_madd_epi16(dif, dif);
                sum                   = _mm256_add_epi32(sum, sqr);

                src += src_stride;
                ref1 += ref1_stride;
                ref2 += ref2_stride;
            } while (--y);
        } else if (width == 32) {
            __m512i sum_512 = _mm512_setzero_si512();

            do {
                ssd32_avx512(src, ref1, ref2, &sum_512);
                src += src_stride;
                ref1 += ref1_stride;
                ref2 += ref2_stride;
            } while (--y);

            const __m256i sum0_256 = _mm512_castsi512_si256(sum_512);
            const __m256i sum1_256 = _mm512_extracti64x4_epi64(sum_512, 1);
            sum                    = _mm256_add_epi32(sum0_256, sum1_256);
        } else if (width == 64) {
            __m512i sum_512 = _mm512_setzero_si512();

            do {
                ssd64_avx512(src, ref1, ref2, &sum_512);
                src += src_stride;
                ref1 += ref1_stride;
                ref2 += ref2_stride;
            } while (--y);

            const __m256i sum0_256 = _mm512_castsi512_si256(sum_512);
            const __m256i sum1_256 = _mm512_extracti64x4_epi64(sum_512, 1);
            sum                    = _mm256_add_epi32(sum0_256, sum1_256);
        } else {
            sum = _mm256_setzero_si256();

            do {
                uint32_t x = 0;
                do {
                    ssd8x2_avx2(
                        src + x, src_stride, ref1 + x, ref1_stride, ref2 + x, ref2_stride, &sum);
                    x += 8;
                } while (x < width);

                src += 2 * src_stride;
                ref1 += 2 * ref1_stride;
                ref2 += 2 * ref2_stride;
                y -= 2;
            } while (y);
        }

        const __m128i sum0_128 = _mm256_castsi256_si128(sum);
        const __m128i sum1_128 = _mm256_extracti128_si256(sum, 1);
        sum_128                = _mm_add_epi32(sum0_128, sum1_128);
    }

    sum_128 = _mm_add_epi32(sum_128, _mm_srli_si128(sum_128, 8));
    sum_128 = _mm_add_epi32(sum_128, _mm_srli_si128(sum_128, 4));
    return _mm_cvtsi128_si32(sum_128);
}

#endif // !NON_AVX512_SUPPORT
