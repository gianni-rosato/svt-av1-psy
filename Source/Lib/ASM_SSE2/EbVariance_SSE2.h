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

#ifndef EBVARIANCE_SSE2_H
#define EBVARIANCE_SSE2_H

#include "EbDefinitions.h"
#include <assert.h>
#include <emmintrin.h>  // SSE2
#include "aom_dsp_rtcd.h"
#include "synonyms.h"

// Read 4 samples from each of row and row + 1. Interleave the two rows and
// zero-extend them to 16 bit samples stored in the lower half of an SSE
// register.
//static __m128i read64(const uint8_t *p, int32_t stride, int32_t row) {
//  __m128i row1 = xx_loadl_32(p + (row + 1) * stride);
//  return _mm_unpacklo_epi8(_mm_unpacklo_epi8(row0, row1), _mm_setzero_si128());
//}

static INLINE __m128i load4x2_sse2(const uint8_t *const p, const int32_t stride) {
    const __m128i p0 = _mm_cvtsi32_si128(*(const uint32_t *)(p + 0 * stride));
    const __m128i p1 = _mm_cvtsi32_si128(*(const uint32_t *)(p + 1 * stride));
    return _mm_unpacklo_epi8(_mm_unpacklo_epi32(p0, p1), _mm_setzero_si128());
}

static INLINE __m128i load8_8to16_sse2(const uint8_t *const p) {
    const __m128i p0 = _mm_loadl_epi64((const __m128i *)p);
    return _mm_unpacklo_epi8(p0, _mm_setzero_si128());
}

// Accumulate 4 32bit numbers in val to 1 32bit number
static INLINE uint32_t add32x4_sse2(__m128i val) {
    val = _mm_add_epi32(val, _mm_srli_si128(val, 8));
    val = _mm_add_epi32(val, _mm_srli_si128(val, 4));
    return _mm_cvtsi128_si32(val);
}

static INLINE void variance_kernel_sse2(const __m128i src, const __m128i ref,
    __m128i *const sse) {
    const __m128i diff = _mm_sub_epi16(src, ref);
    *sse = _mm_add_epi32(*sse, _mm_madd_epi16(diff, diff));
}

// Can handle 128 pixels' diff sum (such as 8x16 or 16x8)
// Slightly faster than variance_final_256_pel_sse2()
// diff sum of 128 pixels can still fit in 16bit integer
static INLINE void variance_final_128_pel_sse2(__m128i vsse,
    uint32_t *const sse) {
    *sse = add32x4_sse2(vsse);
}

// Can handle 256 pixels' diff sum (such as 16x16)
static INLINE void variance_final_256_pel_sse2(__m128i vsse,
    uint32_t *const sse) {
    *sse = add32x4_sse2(vsse);
}

static INLINE void variance4_sse2(const uint8_t *src, const int32_t src_stride,
    const uint8_t *ref, const int32_t ref_stride,
    const int32_t h, __m128i *const sse) {
    assert(h <= 256);  // May overflow for larger height.

    for (int32_t i = 0; i < h; i += 2) {
        const __m128i s = load4x2_sse2(src, src_stride);
        const __m128i r = load4x2_sse2(ref, ref_stride);

        variance_kernel_sse2(s, r, sse);
        src += 2 * src_stride;
        ref += 2 * ref_stride;
    }
}

static INLINE void variance8_sse2(const uint8_t *src, const int32_t src_stride,
    const uint8_t *ref, const int32_t ref_stride,
    const int32_t h, __m128i *const sse) {
    assert(h <= 128);  // May overflow for larger height.
    for (int32_t i = 0; i < h; i++) {
        const __m128i s = load8_8to16_sse2(src);
        const __m128i r = load8_8to16_sse2(ref);

        variance_kernel_sse2(s, r, sse);
        src += src_stride;
        ref += ref_stride;
    }
}

#endif  // EBVARIANCE_SSE2_H
