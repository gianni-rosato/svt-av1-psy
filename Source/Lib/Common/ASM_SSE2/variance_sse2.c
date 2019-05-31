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

// Can handle 128 pixels' diff sum (such as 8x16 or 16x8)
// Slightly faster than variance_final_256_pel_sse2()
// diff sum of 128 pixels can still fit in 16bit integer
static INLINE void variance_final_128_pel_sse2(__m128i vsse, __m128i vsum,
    unsigned int *const sse,
    int *const sum) {
    *sse = add32x4_sse2(vsse);

    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 4));
    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 2));
    *sum = (int16_t)_mm_extract_epi16(vsum, 0);
}

// Can handle 256 pixels' diff sum (such as 16x16)
static INLINE void variance_final_256_pel_sse2(__m128i vsse, __m128i vsum,
    unsigned int *const sse,
    int *const sum) {
    *sse = add32x4_sse2(vsse);

    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 4));
    *sum = (int16_t)_mm_extract_epi16(vsum, 0);
    *sum += (int16_t)_mm_extract_epi16(vsum, 1);
}

static INLINE void variance_kernel_sse2(const __m128i src, const __m128i ref,
    __m128i *const sse,
    __m128i *const sum) {
    const __m128i diff = _mm_sub_epi16(src, ref);
    *sse = _mm_add_epi32(*sse, _mm_madd_epi16(diff, diff));
    *sum = _mm_add_epi16(*sum, diff);
}

static INLINE void variance4_sse2(const uint8_t *src, const int src_stride,
    const uint8_t *ref, const int ref_stride,
    const int h, __m128i *const sse,
    __m128i *const sum) {
    assert(h <= 256);  // May overflow for larger height.
    *sum = _mm_setzero_si128();

    for (int i = 0; i < h; i += 2) {
        const __m128i s = load4x2_sse2(src, src_stride);
        const __m128i r = load4x2_sse2(ref, ref_stride);

        variance_kernel_sse2(s, r, sse, sum);
        src += 2 * src_stride;
        ref += 2 * ref_stride;
    }
}

static INLINE void variance8_sse2(const uint8_t *src, const int src_stride,
    const uint8_t *ref, const int ref_stride,
    const int h, __m128i *const sse,
    __m128i *const sum) {
    assert(h <= 128);  // May overflow for larger height.
    *sum = _mm_setzero_si128();
    for (int i = 0; i < h; i++) {
        const __m128i s = load8_8to16_sse2(src);
        const __m128i r = load8_8to16_sse2(ref);

        variance_kernel_sse2(s, r, sse, sum);
        src += src_stride;
        ref += ref_stride;
    }
}

#define AOM_VAR_NO_LOOP_SSE2(bw, bh, bits, max_pixels)                        \
  unsigned int aom_variance##bw##x##bh##_sse2(                                \
      const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride, \
      unsigned int *sse) {                                                    \
    __m128i vsse = _mm_setzero_si128();                                       \
    __m128i vsum;                                                             \
    int sum = 0;                                                              \
    variance##bw##_sse2(src, src_stride, ref, ref_stride, bh, &vsse, &vsum);  \
    variance_final_##max_pixels##_pel_sse2(vsse, vsum, sse, &sum);            \
    assert(sum <= 255 * bw * bh);                                             \
    assert(sum >= -255 * bw * bh);                                            \
    return *sse - (uint32_t)(((int64_t)sum * sum) >> bits);                   \
  }

AOM_VAR_NO_LOOP_SSE2(4, 4, 4, 128);
AOM_VAR_NO_LOOP_SSE2(4, 8, 5, 128);
AOM_VAR_NO_LOOP_SSE2(4, 16, 6, 128);

AOM_VAR_NO_LOOP_SSE2(8, 4, 5, 128);
AOM_VAR_NO_LOOP_SSE2(8, 8, 6, 128);
AOM_VAR_NO_LOOP_SSE2(8, 16, 7, 128);
AOM_VAR_NO_LOOP_SSE2(8, 32, 8, 256);
