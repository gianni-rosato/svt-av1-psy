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
#include "synonyms.h"

#include <assert.h>
#include <immintrin.h>

// Loads and stores to do away with the tedium of casting the address
// to the right type.
static INLINE int32_t xx_hsum_epi32_si32(__m128i v_d) {
    v_d = _mm_hadd_epi32(v_d, v_d);
    v_d = _mm_hadd_epi32(v_d, v_d);
    return _mm_cvtsi128_si32(v_d);
}
////////////////////////////////////////////////////////////////////////////////
// 8 bit
////////////////////////////////////////////////////////////////////////////////

static INLINE unsigned int obmc_sad_w4_avx2(const uint8_t *pre, const int pre_stride, const int32_t *wsrc,
                                            const int32_t *mask, const int height) {
    int           n        = 0;
    __m256i       v_sad_d  = _mm256_setzero_si256();
    const __m256i v_bias_d = _mm256_set1_epi32((1 << 12) >> 1);

    do {
        const __m128i v_p_b_0 = xx_loadl_32(pre);
        const __m128i v_p_b_1 = xx_loadl_32(pre + pre_stride);
        const __m128i v_p_b   = _mm_unpacklo_epi32(v_p_b_0, v_p_b_1);
        const __m256i v_m_d   = _mm256_lddqu_si256((__m256i *)(mask + n));
        const __m256i v_w_d   = _mm256_lddqu_si256((__m256i *)(wsrc + n));

        const __m256i v_p_d = _mm256_cvtepu8_epi32(v_p_b);

        // Values in both pre and mask fit in 15 bits, and are packed at 32 bit
        // boundaries. We use pmaddwd, as it has lower latency on Haswell
        // than pmulld but produces the same result with these inputs.
        const __m256i v_pm_d = _mm256_madd_epi16(v_p_d, v_m_d);

        const __m256i v_diff_d    = _mm256_sub_epi32(v_w_d, v_pm_d);
        const __m256i v_absdiff_d = _mm256_abs_epi32(v_diff_d);

        // Rounded absolute difference
        const __m256i v_tmp_d = _mm256_add_epi32(v_absdiff_d, v_bias_d);
        const __m256i v_rad_d = _mm256_srli_epi32(v_tmp_d, 12);

        v_sad_d = _mm256_add_epi32(v_sad_d, v_rad_d);

        n += 8;
        pre += pre_stride << 1;
    } while (n < 8 * (height >> 1));

    __m128i v_sad_d_0 = _mm256_castsi256_si128(v_sad_d);
    __m128i v_sad_d_1 = _mm256_extracti128_si256(v_sad_d, 1);
    v_sad_d_0         = _mm_add_epi32(v_sad_d_0, v_sad_d_1);
    return xx_hsum_epi32_si32(v_sad_d_0);
}

static INLINE unsigned int obmc_sad_w8n_avx2(const uint8_t *pre, const int pre_stride, const int32_t *wsrc,
                                             const int32_t *mask, const int width, const int height) {
    const int     pre_step = pre_stride - width;
    int           n        = 0;
    __m256i       v_sad_d  = _mm256_setzero_si256();
    const __m256i v_bias_d = _mm256_set1_epi32((1 << 12) >> 1);
    assert(width >= 8);
    assert(IS_POWER_OF_TWO(width));

    do {
        const __m128i v_p0_b = xx_loadl_64(pre + n);
        const __m256i v_m0_d = _mm256_lddqu_si256((__m256i *)(mask + n));
        const __m256i v_w0_d = _mm256_lddqu_si256((__m256i *)(wsrc + n));

        const __m256i v_p0_d = _mm256_cvtepu8_epi32(v_p0_b);

        // Values in both pre and mask fit in 15 bits, and are packed at 32 bit
        // boundaries. We use pmaddwd, as it has lower latency on Haswell
        // than pmulld but produces the same result with these inputs.
        const __m256i v_pm0_d = _mm256_madd_epi16(v_p0_d, v_m0_d);

        const __m256i v_diff0_d    = _mm256_sub_epi32(v_w0_d, v_pm0_d);
        const __m256i v_absdiff0_d = _mm256_abs_epi32(v_diff0_d);

        // Rounded absolute difference
        const __m256i v_tmp_d  = _mm256_add_epi32(v_absdiff0_d, v_bias_d);
        const __m256i v_rad0_d = _mm256_srli_epi32(v_tmp_d, 12);

        v_sad_d = _mm256_add_epi32(v_sad_d, v_rad0_d);

        n += 8;

        if ((n & (width - 1)) == 0)
            pre += pre_step;
    } while (n < width * height);

    __m128i v_sad_d_0 = _mm256_castsi256_si128(v_sad_d);
    __m128i v_sad_d_1 = _mm256_extracti128_si256(v_sad_d, 1);
    v_sad_d_0         = _mm_add_epi32(v_sad_d_0, v_sad_d_1);
    return xx_hsum_epi32_si32(v_sad_d_0);
}

#define OBMCSADWXH(w, h)                                                               \
    unsigned int svt_aom_obmc_sad##w##x##h##_avx2(                                     \
        const uint8_t *pre, int pre_stride, const int32_t *wsrc, const int32_t *msk) { \
        if (w == 4) {                                                                  \
            return obmc_sad_w4_avx2(pre, pre_stride, wsrc, msk, h);                    \
        } else {                                                                       \
            return obmc_sad_w8n_avx2(pre, pre_stride, wsrc, msk, w, h);                \
        }                                                                              \
    }

OBMCSADWXH(128, 128)
OBMCSADWXH(128, 64)
OBMCSADWXH(64, 128)
OBMCSADWXH(64, 64)
OBMCSADWXH(64, 32)
OBMCSADWXH(32, 64)
OBMCSADWXH(32, 32)
OBMCSADWXH(32, 16)
OBMCSADWXH(16, 32)
OBMCSADWXH(16, 16)
OBMCSADWXH(16, 8)
OBMCSADWXH(8, 16)
OBMCSADWXH(8, 8)
OBMCSADWXH(8, 4)
OBMCSADWXH(4, 8)
OBMCSADWXH(4, 4)
OBMCSADWXH(4, 16)
OBMCSADWXH(16, 4)
OBMCSADWXH(8, 32)
OBMCSADWXH(32, 8)
OBMCSADWXH(16, 64)
OBMCSADWXH(64, 16)
