/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include "EbDefinitions.h"
#include <assert.h>
#include <emmintrin.h>  // SSE2
#include "aom_dsp_rtcd.h"
#include "EbVariance_SSE2.h"
#include "synonyms.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
}
#endif

uint32_t aom_get_mb_ss_sse2(const int16_t *src) {
    __m128i vsum = _mm_setzero_si128();
    int32_t i;

    for (i = 0; i < 32; ++i) {
        const __m128i v = xx_loadu_128(src);
        vsum = _mm_add_epi32(vsum, _mm_madd_epi16(v, v));
        src += 8;
    }

    vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 4));
    return _mm_cvtsi128_si32(vsum);
}

