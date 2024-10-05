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
#ifndef AV1_COMMON_X86_AV1_TXFM_COMMON_AVX2_H_
#define AV1_COMMON_X86_AV1_TXFM_COMMON_AVX2_H_

#include <immintrin.h>
#include <stdint.h>
#include "inv_transforms.h"

#ifdef __cplusplus
extern "C" {
#endif

#define pair_set_epi16(a, b) _mm_set1_epi32((int32_t)(((uint16_t)(a)) | (((uint32_t)(b)) << 16)))

// Reverse the 8 16 bit words in __m128i
static INLINE __m128i mm_reverse_epi16(const __m128i x) {
    const __m128i a = _mm_shufflelo_epi16(x, 0x1b);
    const __m128i b = _mm_shufflehi_epi16(a, 0x1b);
    return _mm_shuffle_epi32(b, 0x4e);
}

static INLINE __m256i pair_set_w16_epi16(int16_t a, int16_t b) {
    return _mm256_set1_epi32((int32_t)(((uint16_t)(a)) | (((uint32_t)(b)) << 16)));
}

static INLINE void btf_16_w16_avx2(const __m256i w0, const __m256i w1, __m256i *in0, __m256i *in1, const __m256i _r,
                                   const int32_t cos_bit) {
    __m256i t0 = _mm256_unpacklo_epi16(*in0, *in1);
    __m256i t1 = _mm256_unpackhi_epi16(*in0, *in1);
    __m256i u0 = _mm256_madd_epi16(t0, w0);
    __m256i u1 = _mm256_madd_epi16(t1, w0);
    __m256i v0 = _mm256_madd_epi16(t0, w1);
    __m256i v1 = _mm256_madd_epi16(t1, w1);

    __m256i a0 = _mm256_add_epi32(u0, _r);
    __m256i a1 = _mm256_add_epi32(u1, _r);
    __m256i b0 = _mm256_add_epi32(v0, _r);
    __m256i b1 = _mm256_add_epi32(v1, _r);

    __m256i c0 = _mm256_srai_epi32(a0, cos_bit);
    __m256i c1 = _mm256_srai_epi32(a1, cos_bit);
    __m256i d0 = _mm256_srai_epi32(b0, cos_bit);
    __m256i d1 = _mm256_srai_epi32(b1, cos_bit);

    *in0 = _mm256_packs_epi32(c0, c1);
    *in1 = _mm256_packs_epi32(d0, d1);
}

static INLINE void btf_16_adds_subs_avx2(__m256i *in0, __m256i *in1) {
    const __m256i _in0 = *in0;
    const __m256i _in1 = *in1;
    *in0               = _mm256_adds_epi16(_in0, _in1);
    *in1               = _mm256_subs_epi16(_in0, _in1);
}

static INLINE void btf_16_adds_subs_out_avx2(__m256i *out0, __m256i *out1, __m256i in0, __m256i in1) {
    const __m256i _in0 = in0;
    const __m256i _in1 = in1;
    *out0              = _mm256_adds_epi16(_in0, _in1);
    *out1              = _mm256_subs_epi16(_in0, _in1);
}

static INLINE __m256i load_32bit_to_16bit_w16_avx2(const int32_t *a) {
    const __m256i a_low  = _mm256_lddqu_si256((const __m256i *)a);
    const __m256i a_high = _mm256_lddqu_si256((const __m256i *)(a + 8));
    const __m256i b      = _mm256_packs_epi32(a_low, a_high);
    return _mm256_permute4x64_epi64(b, 0xD8);
}

static INLINE void load_buffer_32bit_to_16bit_w16_avx2(const int32_t *in, int stride, __m256i *out, int out_size) {
    for (int i = 0; i < out_size; ++i) { out[i] = load_32bit_to_16bit_w16_avx2(in + i * stride); }
}

static INLINE void transpose2_8x8_avx2(const __m256i *const in,
                                       __m256i *const out) {
  __m256i t[16], u[16];
  // (1st, 2nd) ==> (lo, hi)
  //   (0, 1)   ==>  (0, 1)
  //   (2, 3)   ==>  (2, 3)
  //   (4, 5)   ==>  (4, 5)
  //   (6, 7)   ==>  (6, 7)
  for (int i = 0; i < 4; i++) {
    t[2 * i] = _mm256_unpacklo_epi16(in[2 * i], in[2 * i + 1]);
    t[2 * i + 1] = _mm256_unpackhi_epi16(in[2 * i], in[2 * i + 1]);
  }
  // (1st, 2nd) ==> (lo, hi)
  //   (0, 2)   ==>  (0, 2)
  //   (1, 3)   ==>  (1, 3)
  //   (4, 6)   ==>  (4, 6)
  //   (5, 7)   ==>  (5, 7)
  for (int i = 0; i < 2; i++) {
    u[i] = _mm256_unpacklo_epi32(t[i], t[i + 2]);
    u[i + 2] = _mm256_unpackhi_epi32(t[i], t[i + 2]);
    u[i + 4] = _mm256_unpacklo_epi32(t[i + 4], t[i + 6]);
    u[i + 6] = _mm256_unpackhi_epi32(t[i + 4], t[i + 6]);
  }
  // (1st, 2nd) ==> (lo, hi)
  //   (0, 4)   ==>  (0, 1)
  //   (1, 5)   ==>  (4, 5)
  //   (2, 6)   ==>  (2, 3)
  //   (3, 7)   ==>  (6, 7)
  for (int i = 0; i < 2; i++) {
    out[2 * i] = _mm256_unpacklo_epi64(u[2 * i], u[2 * i + 4]);
    out[2 * i + 1] = _mm256_unpackhi_epi64(u[2 * i], u[2 * i + 4]);
    out[2 * i + 4] = _mm256_unpacklo_epi64(u[2 * i + 1], u[2 * i + 5]);
    out[2 * i + 5] = _mm256_unpackhi_epi64(u[2 * i + 1], u[2 * i + 5]);
  }
}

static INLINE void transpose_16bit_16x16_avx2(const __m256i *const in, __m256i *const out) {
      __m256i t[16];
#define LOADL(idx)                                                            \
  t[idx] = _mm256_castsi128_si256(_mm_load_si128((__m128i const *)&in[idx])); \
  t[idx] = _mm256_inserti128_si256(                                           \
      t[idx], _mm_load_si128((__m128i const *)&in[idx + 8]), 1);
#define LOADR(idx)                                                           \
  t[8 + idx] =                                                               \
      _mm256_castsi128_si256(_mm_load_si128((__m128i const *)&in[idx] + 1)); \
  t[8 + idx] = _mm256_inserti128_si256(                                      \
      t[8 + idx], _mm_load_si128((__m128i const *)&in[idx + 8] + 1), 1);
  // load left 8x16
  LOADL(0)
  LOADL(1)
  LOADL(2)
  LOADL(3)
  LOADL(4)
  LOADL(5)
  LOADL(6)
  LOADL(7)
  // load right 8x16
  LOADR(0)
  LOADR(1)
  LOADR(2)
  LOADR(3)
  LOADR(4)
  LOADR(5)
  LOADR(6)
  LOADR(7)
  // get the top 16x8 result
  transpose2_8x8_avx2(t, out);
  // get the bottom 16x8 result
  transpose2_8x8_avx2(&t[8], &out[8]);
}

static INLINE void flip_buf_avx2(__m256i *in, __m256i *out, int size) {
    for (int i = 0; i < size; ++i) { out[size - i - 1] = in[i]; }
}

static INLINE void round_shift_16bit_w16_avx2(__m256i *in, int size, int bit) {
    if (bit < 0) {
        bit           = -bit;
        __m256i round = _mm256_set1_epi16(1 << (bit - 1));
        for (int i = 0; i < size; ++i) {
            in[i] = _mm256_adds_epi16(in[i], round);
            in[i] = _mm256_srai_epi16(in[i], bit);
        }
    } else if (bit > 0) {
        for (int i = 0; i < size; ++i) { in[i] = _mm256_slli_epi16(in[i], bit); }
    }
}

static INLINE void av1_round_shift_array_32_avx2(__m256i *input, __m256i *output, const int32_t size,
                                                 const int32_t bit) {
    int32_t i;
    if (bit > 0) {
        const __m256i round = _mm256_set1_epi32(1 << (bit - 1));
        __m256i       r0;
        for (i = 0; i < size; i++) {
            r0        = _mm256_add_epi32(input[i], round);
            output[i] = _mm256_srai_epi32(r0, bit);
        }
    } else {
        for (i = 0; i < size; i++) output[i] = _mm256_slli_epi32(input[i], -bit);
    }
}

static INLINE void av1_round_shift_rect_array_32_avx2(__m256i *input, __m256i *output, const int32_t size,
                                                      const int32_t bit, const int32_t val) {
    const __m256i sqrt2  = _mm256_set1_epi32(val);
    const __m256i round2 = _mm256_set1_epi32(1 << (new_sqrt2_bits - 1));
    int32_t       i;
    if (bit > 0) {
        const __m256i round1 = _mm256_set1_epi32(1 << (bit - 1));
        __m256i       r0, r1, r2, r3;
        for (i = 0; i < size; i++) {
            r0        = _mm256_add_epi32(input[i], round1);
            r1        = _mm256_srai_epi32(r0, bit);
            r2        = _mm256_mullo_epi32(sqrt2, r1);
            r3        = _mm256_add_epi32(r2, round2);
            output[i] = _mm256_srai_epi32(r3, new_sqrt2_bits);
        }
    } else {
        __m256i r0, r1, r2;
        for (i = 0; i < size; i++) {
            r0        = _mm256_slli_epi32(input[i], -bit);
            r1        = _mm256_mullo_epi32(sqrt2, r0);
            r2        = _mm256_add_epi32(r1, round2);
            output[i] = _mm256_srai_epi32(r2, new_sqrt2_bits);
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif // AV1_COMMON_X86_AV1_TXFM_COMMON_AVX2_H_
