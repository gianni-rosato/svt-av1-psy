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

#ifndef AOM_DSP_X86_SYNONYMS_H_
#define AOM_DSP_X86_SYNONYMS_H_

#include <immintrin.h>
#include "EbDefinitions.h"

 /**
  * Various reusable shorthands for x86 SIMD intrinsics.
  *
  * Intrinsics prefixed with xx_ operate on or return 128bit XMM registers.
  * Intrinsics prefixed with yy_ operate on or return 256bit YMM registers.
  */

static INLINE __m128i xx_loadl_64(const void *a) {
    return _mm_loadl_epi64((const __m128i *)a);
}

static INLINE __m128i xx_loadu_128(const void *a) {
    return _mm_loadu_si128((const __m128i *)a);
}

static INLINE void xx_storel_32(void *const a, const __m128i v) {
    *(uint32_t *)a = _mm_cvtsi128_si32(v);
}

static INLINE void xx_storel_64(void *const a, const __m128i v) {
    _mm_storel_epi64((__m128i *)a, v);
}

static INLINE void xx_storeu_128(void *const a, const __m128i v) {
    _mm_storeu_si128((__m128i *)a, v);
}

// The _mm_set_epi64x() intrinsic is undefined for some Visual Studio
// compilers. The following function is equivalent to _mm_set_epi64x()
// acting on 32-bit integers.
static INLINE __m128i xx_set_64_from_32i(int32_t e1, int32_t e0) {
#if defined(_MSC_VER) && _MSC_VER < 1900
    return _mm_set_epi32(0, e1, 0, e0);
#else
    return _mm_set_epi64x((uint32_t)e1, (uint32_t)e0);
#endif
}

// The _mm_set1_epi64x() intrinsic is undefined for some Visual Studio
// compilers. The following function is equivalent to _mm_set1_epi64x()
// acting on a 32-bit integer.
static INLINE __m128i xx_set1_64_from_32i(int32_t a) {
#if defined(_MSC_VER) && _MSC_VER < 1900
    return _mm_set_epi32(0, a, 0, a);
#else
    return _mm_set1_epi64x((uint32_t)a);
#endif
}

static INLINE __m128i xx_round_epu16(__m128i v_val_w) {
    return _mm_avg_epu16(v_val_w, _mm_setzero_si128());
}

static INLINE __m128i xx_roundn_epu16(__m128i v_val_w, int32_t bits) {
    const __m128i v_s_w = _mm_srli_epi16(v_val_w, bits - 1);
    return _mm_avg_epu16(v_s_w, _mm_setzero_si128());
}

static INLINE __m128i xx_roundn_epu32(__m128i v_val_d, int32_t bits) {
    const __m128i v_bias_d = _mm_set1_epi32((1 << bits) >> 1);
    const __m128i v_tmp_d = _mm_add_epi32(v_val_d, v_bias_d);
    return _mm_srli_epi32(v_tmp_d, bits);
}

// This is equivalent to ROUND_POWER_OF_TWO(v_val_d, bits)
static INLINE __m128i xx_roundn_epi32_unsigned(__m128i v_val_d, int32_t bits) {
    const __m128i v_bias_d = _mm_set1_epi32((1 << bits) >> 1);
    const __m128i v_tmp_d = _mm_add_epi32(v_val_d, v_bias_d);
    return _mm_srai_epi32(v_tmp_d, bits);
}

// This is equivalent to ROUND_POWER_OF_TWO_SIGNED(v_val_d, bits)
static INLINE __m128i xx_roundn_epi32(__m128i v_val_d, int32_t bits) {
    const __m128i v_bias_d = _mm_set1_epi32((1 << bits) >> 1);
    const __m128i v_sign_d = _mm_srai_epi32(v_val_d, 31);
    const __m128i v_tmp_d =
        _mm_add_epi32(_mm_add_epi32(v_val_d, v_bias_d), v_sign_d);
    return _mm_srai_epi32(v_tmp_d, bits);
}

static INLINE __m128i xx_roundn_epi16(__m128i v_val_d, int32_t bits) {
    const __m128i v_bias_d = _mm_set1_epi16((1 << bits) >> 1);
    const __m128i v_sign_d = _mm_srai_epi16(v_val_d, 15);
    const __m128i v_tmp_d =
        _mm_add_epi16(_mm_add_epi16(v_val_d, v_bias_d), v_sign_d);
    return _mm_srai_epi16(v_tmp_d, bits);
}

// Note:
// _mm256_insert_epi16 intrinsics is available from vs2017.
// We define this macro for vs2015 and earlier. The
// intrinsics used here are in vs2015 document:
// https://msdn.microsoft.com/en-us/library/hh977022.aspx
// Input parameters:
// a: __m256i,
// d: int16_t,
// indx: imm8 (0 - 15)
//#if _MSC_VER <= 1900
#ifdef _WIN32
#if _MSC_VER < 1910
#define _mm256_insert_epi16(a, d, indx)                                      \
  _mm256_insertf128_si256(                                                   \
      a,                                                                     \
      _mm_insert_epi16(_mm256_extractf128_si256(a, indx >> 3), d, indx % 8), \
      indx >> 3)

static INLINE int32_t _mm256_extract_epi32(__m256i a, const int32_t i) {
    return a.m256i_i32[i & 7];
}

static INLINE int32_t _mm256_extract_epi16(__m256i a, const int32_t i) {
    return a.m256i_i16[i & 15];
}

static INLINE __m256i _mm256_insert_epi32(__m256i a, int32_t b, const int32_t i) {
    __m256i c = a;
    c.m256i_i32[i & 7] = b;
    return c;
}
#endif
#endif

#endif  // AOM_DSP_X86_SYNONYMS_H_
