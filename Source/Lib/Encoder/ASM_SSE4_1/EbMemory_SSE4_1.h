/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMemory_SSE4_1_h
#define EbMemory_SSE4_1_h

#include "EbDefinitions.h"
#include "smmintrin.h"

static INLINE __m128i load8bit_4x2_sse4_1(const void *const src, const ptrdiff_t strideInByte) {
    const __m128i s = _mm_cvtsi32_si128(*(int32_t *)((uint8_t *)src));
    return _mm_insert_epi32(s, *(int32_t *)((uint8_t *)src + strideInByte), 1);
}

static INLINE __m128i load_u8_4x2_sse4_1(const uint8_t *const src, const ptrdiff_t stride) {
    return load8bit_4x2_sse4_1(src, sizeof(*src) * stride);
}

static INLINE __m128i load_u16_2x2_sse4_1(const uint16_t *const src, const ptrdiff_t stride) {
    return load8bit_4x2_sse4_1(src, sizeof(*src) * stride);
}

#endif // EbMemory_SSE4_1_h
