/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <immintrin.h> /* AVX2 */

#include "EbDefinitions.h"
#include "synonyms.h"
#include "synonyms_avx2.h"

static INLINE __m256i txb_init_levels_avx2(const TranLow *const coeff) {
    const __m256i idx   = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);
    const __m256i c0    = yy_loadu_256(coeff + 0 * 8);
    const __m256i c1    = yy_loadu_256(coeff + 1 * 8);
    const __m256i c2    = yy_loadu_256(coeff + 2 * 8);
    const __m256i c3    = yy_loadu_256(coeff + 3 * 8);
    const __m256i c01   = _mm256_packs_epi32(c0, c1);
    const __m256i c23   = _mm256_packs_epi32(c2, c3);
    const __m256i abs01 = _mm256_abs_epi16(c01);
    const __m256i abs23 = _mm256_abs_epi16(c23);
    const __m256i res   = _mm256_packs_epi16(abs01, abs23);
    return _mm256_permutevar8x32_epi32(res, idx);
}

void eb_av1_txb_init_levels_avx2(const TranLow *const coeff, const int32_t width,
                                 const int32_t height, uint8_t *const levels) {
    const TranLow *cf      = coeff;
    const __m128i  x_zeros = _mm_setzero_si128();
    const __m256i  y_zeros = _mm256_setzero_si256();
    uint8_t *      ls      = levels;
    int32_t        i       = height;

    if (width == 4) {
        xx_storeu_128(ls - 16, x_zeros);

        do {
            const __m256i idx   = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
            const __m256i c0    = yy_loadu_256(cf);
            const __m256i c1    = yy_loadu_256(cf + 8);
            const __m256i c01   = _mm256_packs_epi32(c0, c1);
            const __m256i abs01 = _mm256_abs_epi16(c01);
            const __m256i res_  = _mm256_packs_epi16(abs01, y_zeros);
            const __m256i res   = _mm256_permutevar8x32_epi32(res_, idx);
            yy_storeu_256(ls, res);
            cf += 4 * 4;
            ls += 4 * 8;
            i -= 4;
        } while (i);

        yy_storeu_256(ls, y_zeros);
    } else if (width == 8) {
        yy_storeu_256(ls - 24, y_zeros);

        do {
            const __m256i res  = txb_init_levels_avx2(cf);
            const __m128i res0 = _mm256_castsi256_si128(res);
            const __m128i res1 = _mm256_extracti128_si256(res, 1);
            xx_storel_64(ls + 0 * 12 + 0, res0);
            *(int32_t *)(ls + 0 * 12 + 8) = 0;
            _mm_storeh_epi64((__m128i *)(ls + 1 * 12 + 0), res0);
            *(int32_t *)(ls + 1 * 12 + 8) = 0;
            xx_storel_64(ls + 2 * 12 + 0, res1);
            *(int32_t *)(ls + 2 * 12 + 8) = 0;
            _mm_storeh_epi64((__m128i *)(ls + 3 * 12 + 0), res1);
            *(int32_t *)(ls + 3 * 12 + 8) = 0;
            cf += 4 * 8;
            ls += 4 * 12;
            i -= 4;
        } while (i);

        yy_storeu_256(ls + 0 * 32, y_zeros);
        xx_storeu_128(ls + 1 * 32, x_zeros);
    } else if (width == 16) {
        yy_storeu_256(ls - 40, y_zeros);
        xx_storel_64(ls - 8, x_zeros);

        do {
            const __m256i res  = txb_init_levels_avx2(cf);
            const __m128i res0 = _mm256_castsi256_si128(res);
            const __m128i res1 = _mm256_extracti128_si256(res, 1);
            xx_storeu_128(ls, res0);
            *(int32_t *)(ls + 16) = 0;
            xx_storeu_128(ls + 20, res1);
            *(int32_t *)(ls + 20 + 16) = 0;
            cf += 2 * 16;
            ls += 2 * 20;
            i -= 2;
        } while (i);

        yy_storeu_256(ls + 0 * 32, y_zeros);
        yy_storeu_256(ls + 1 * 32, y_zeros);
        xx_storeu_128(ls + 2 * 32, x_zeros);
    } else {
        yy_storeu_256(ls - 72, y_zeros);
        yy_storeu_256(ls - 40, y_zeros);
        xx_storel_64(ls - 8, x_zeros);

        do {
            const __m256i res = txb_init_levels_avx2(cf);
            yy_storeu_256(ls, res);
            *(int32_t *)(ls + 32) = 0;
            cf += 32;
            ls += 36;
        } while (--i);

        yy_storeu_256(ls + 0 * 32, y_zeros);
        yy_storeu_256(ls + 1 * 32, y_zeros);
        yy_storeu_256(ls + 2 * 32, y_zeros);
        yy_storeu_256(ls + 3 * 32, y_zeros);
        xx_storeu_128(ls + 4 * 32, x_zeros);
    }
}
