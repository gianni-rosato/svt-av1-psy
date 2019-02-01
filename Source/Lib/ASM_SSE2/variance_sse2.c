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

// Accumulate 8 16bit in sum to 4 32bit number
static INLINE __m128i sum_to_32bit_sse2(const __m128i sum) {
    const __m128i sum_lo = _mm_srai_epi32(_mm_unpacklo_epi16(sum, sum), 16);
    const __m128i sum_hi = _mm_srai_epi32(_mm_unpackhi_epi16(sum, sum), 16);
    return _mm_add_epi32(sum_lo, sum_hi);
}

// Can handle 512 pixels' diff sum (such as 16x32 or 32x16)
static INLINE void variance_final_512_pel_sse2(__m128i vsse, __m128i vsum,
    uint32_t *const sse,
    int32_t *const sum) {
    *sse = add32x4_sse2(vsse);

    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_unpacklo_epi16(vsum, vsum);
    vsum = _mm_srai_epi32(vsum, 16);
    *sum = add32x4_sse2(vsum);
}

// Can handle 1024 pixels' diff sum (such as 32x32)
static INLINE void variance_final_1024_pel_sse2(__m128i vsse, __m128i vsum,
    uint32_t *const sse,
    int32_t *const sum) {
    *sse = add32x4_sse2(vsse);

    vsum = sum_to_32bit_sse2(vsum);
    *sum = add32x4_sse2(vsum);
}

static INLINE void variance16_kernel_sse2(const uint8_t *const src,
    const uint8_t *const ref,
    __m128i *const sse,
    __m128i *const sum) {
    const __m128i zero = _mm_setzero_si128();
    const __m128i s = _mm_loadu_si128((const __m128i *)src);
    const __m128i r = _mm_loadu_si128((const __m128i *)ref);
    const __m128i src0 = _mm_unpacklo_epi8(s, zero);
    const __m128i ref0 = _mm_unpacklo_epi8(r, zero);
    const __m128i src1 = _mm_unpackhi_epi8(s, zero);
    const __m128i ref1 = _mm_unpackhi_epi8(r, zero);

    variance_kernel_sse2(src0, ref0, sse, sum);
    variance_kernel_sse2(src1, ref1, sse, sum);
}

static INLINE void variance16_sse2(const uint8_t *src, const int32_t src_stride,
    const uint8_t *ref, const int32_t ref_stride,
    const int32_t h, __m128i *const sse,
    __m128i *const sum) {
    assert(h <= 64);  // May overflow for larger height.
    *sum = _mm_setzero_si128();

    for (int32_t i = 0; i < h; ++i) {
        variance16_kernel_sse2(src, ref, sse, sum);
        src += src_stride;
        ref += ref_stride;
    }
}

static INLINE void variance32_sse2(const uint8_t *src, const int32_t src_stride,
    const uint8_t *ref, const int32_t ref_stride,
    const int32_t h, __m128i *const sse,
    __m128i *const sum) {
    assert(h <= 32);  // May overflow for larger height.
    // Don't initialize sse here since it's an accumulation.
    *sum = _mm_setzero_si128();

    for (int32_t i = 0; i < h; ++i) {
        variance16_kernel_sse2(src + 0, ref + 0, sse, sum);
        variance16_kernel_sse2(src + 16, ref + 16, sse, sum);
        src += src_stride;
        ref += ref_stride;
    }
}

static INLINE void variance64_sse2(const uint8_t *src, const int32_t src_stride,
    const uint8_t *ref, const int32_t ref_stride,
    const int32_t h, __m128i *const sse,
    __m128i *const sum) {
    assert(h <= 16);  // May overflow for larger height.
    *sum = _mm_setzero_si128();

    for (int32_t i = 0; i < h; ++i) {
        variance16_kernel_sse2(src + 0, ref + 0, sse, sum);
        variance16_kernel_sse2(src + 16, ref + 16, sse, sum);
        variance16_kernel_sse2(src + 32, ref + 32, sse, sum);
        variance16_kernel_sse2(src + 48, ref + 48, sse, sum);
        src += src_stride;
        ref += ref_stride;
    }
}

static INLINE void variance128_sse2(const uint8_t *src, const int32_t src_stride,
    const uint8_t *ref, const int32_t ref_stride,
    const int32_t h, __m128i *const sse,
    __m128i *const sum) {
    assert(h <= 8);  // May overflow for larger height.
    *sum = _mm_setzero_si128();

    for (int32_t i = 0; i < h; ++i) {
        for (int32_t j = 0; j < 4; ++j) {
            const int32_t offset0 = j << 5;
            const int32_t offset1 = offset0 + 16;
            variance16_kernel_sse2(src + offset0, ref + offset0, sse, sum);
            variance16_kernel_sse2(src + offset1, ref + offset1, sse, sum);
        }
        src += src_stride;
        ref += ref_stride;
    }
}