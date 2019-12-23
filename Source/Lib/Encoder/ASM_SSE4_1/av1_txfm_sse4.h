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

#ifndef AV1_TXFM_SSE4_H_
#define AV1_TXFM_SSE4_H_

#include <smmintrin.h>
#include "EbTransforms.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE __m128i av1_round_shift_32_sse4_1(__m128i vec, int32_t bit) {
    __m128i tmp, round;
    round = _mm_set1_epi32(1 << (bit - 1));
    tmp   = _mm_add_epi32(vec, round);
    return _mm_srai_epi32(tmp, bit);
}

static INLINE void av1_round_shift_array_32_sse4_1(__m128i *input, __m128i *output,
                                                   const int32_t size, const int32_t bit) {
    if (bit > 0) {
        int32_t i;
        for (i = 0; i < size; i++) output[i] = av1_round_shift_32_sse4_1(input[i], bit);
    } else {
        int32_t i;
        for (i = 0; i < size; i++) output[i] = _mm_slli_epi32(input[i], -bit);
    }
}

static INLINE void av1_round_shift_rect_array_32_sse4_1(__m128i *input, __m128i *output,
                                                        const int32_t size, const int32_t bit) {
    const __m128i sqrt2 = _mm_set1_epi32(new_sqrt2);
    if (bit > 0) {
        int32_t i;
        for (i = 0; i < size; i++) {
            const __m128i r0 = av1_round_shift_32_sse4_1(input[i], bit);
            const __m128i r1 = _mm_mullo_epi32(sqrt2, r0);
            output[i]        = av1_round_shift_32_sse4_1(r1, new_sqrt2_bits);
        }
    } else {
        int32_t i;
        for (i = 0; i < size; i++) {
            const __m128i r0 = _mm_slli_epi32(input[i], -bit);
            const __m128i r1 = _mm_mullo_epi32(sqrt2, r0);
            output[i]        = av1_round_shift_32_sse4_1(r1, new_sqrt2_bits);
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // AV1_TXFM_SSE4_H_
