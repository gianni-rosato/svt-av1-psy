/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_SSE2_h
#define EbPictureOperators_SSE2_h

#include <emmintrin.h>
#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif


static INLINE int32_t hadd32_sse2_intrin(const __m128i src) {
    const __m128i dst0 = _mm_add_epi32(src, _mm_srli_si128(src, 8));
    const __m128i dst1 = _mm_add_epi32(dst0, _mm_srli_si128(dst0, 4));

    return _mm_cvtsi128_si32(dst1);
}


#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_SSE2_h
