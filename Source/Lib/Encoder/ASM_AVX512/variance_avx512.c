/*
* Copyright(c) 2021 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include "EbDefinitions.h"

#if EN_AVX512_SUPPORT

#include <immintrin.h>

#if VARIANCE_AVX512

static INLINE __m128i mm512_add_hi_lo_epi16(const __m512i val) {
    const __m256i val_lo = _mm512_castsi512_si256(val);
    const __m256i val_hi = _mm512_extracti64x4_epi64(val, 1);
    const __m256i val2   = _mm256_add_epi16(val_lo, val_hi);

    return _mm_add_epi16(_mm256_castsi256_si128(val2), _mm256_extractf128_si256(val2, 1));
}

static INLINE __m128i mm512_sum_to_128_epi32(const __m512i val) {
    const __m256i val_lo = _mm512_castsi512_si256(val);
    const __m256i val_hi = _mm512_extracti64x4_epi64(val, 1);
    const __m256i val2   = _mm256_add_epi32(val_lo, val_hi);

    return _mm_add_epi32(_mm256_castsi256_si128(val2), _mm256_extractf128_si256(val2, 1));
}

static INLINE int variance_final_from_32bit_sum_avx512(__m512i vsse, __m128i vsum,
                                                       unsigned int *const sse) {
    // extract the low lane and add it to the high lane
    const __m128i sse_reg_128 = mm512_sum_to_128_epi32(vsse);

    // unpack sse and sum registers and add
    const __m128i sse_sum_lo = _mm_unpacklo_epi32(sse_reg_128, vsum);
    const __m128i sse_sum_hi = _mm_unpackhi_epi32(sse_reg_128, vsum);
    const __m128i sse_sum    = _mm_add_epi32(sse_sum_lo, sse_sum_hi);

    // perform the final summation and extract the results
    const __m128i res = _mm_add_epi32(sse_sum, _mm_srli_si128(sse_sum, 8));
    *((int *)sse)     = _mm_cvtsi128_si32(res);
    return _mm_extract_epi32(res, 1);
}

static INLINE __m512i sum_to_32bit_avx512(const __m512i sum) {
    const __m512i sum_lo = _mm512_cvtepi16_epi32(_mm512_castsi512_si256(sum));
    const __m512i sum_hi = _mm512_cvtepi16_epi32(_mm512_extracti64x4_epi64(sum, 1));
    return _mm512_add_epi32(sum_lo, sum_hi);
}

static INLINE void variance_kernel_avx512(const __m512i src, const __m512i ref, __m512i *const sse,
                                          __m512i *const sum) {
    const __m512i adj_sub = _mm512_set1_epi16(0xff01); // (1,-1)

    // unpack into pairs of source and reference values
    const __m512i src_ref0 = _mm512_unpacklo_epi8(src, ref);
    const __m512i src_ref1 = _mm512_unpackhi_epi8(src, ref);

    // subtract adjacent elements using src*1 + ref*-1
    const __m512i diff0 = _mm512_maddubs_epi16(src_ref0, adj_sub);
    const __m512i diff1 = _mm512_maddubs_epi16(src_ref1, adj_sub);
    const __m512i madd0 = _mm512_madd_epi16(diff0, diff0);
    const __m512i madd1 = _mm512_madd_epi16(diff1, diff1);

    // add to the running totals
    *sum = _mm512_add_epi16(*sum, _mm512_add_epi16(diff0, diff1));
    *sse = _mm512_add_epi32(*sse, _mm512_add_epi32(madd0, madd1));
}

static INLINE int variance_final_512_avx512(__m512i vsse, __m512i vsum, unsigned int *const sse) {
    // extract the low lane and add it to the high lane
    const __m128i vsum_128  = mm512_add_hi_lo_epi16(vsum);
    const __m128i vsum_64   = _mm_add_epi16(vsum_128, _mm_srli_si128(vsum_128, 8));
    const __m128i sum_int32 = _mm_cvtepi16_epi32(vsum_64);
    return variance_final_from_32bit_sum_avx512(vsse, sum_int32, sse);
}

static INLINE int variance_final_1024_avx512(__m512i vsse, __m512i vsum, unsigned int *const sse) {
    // extract the low lane and add it to the high lane
    const __m128i vsum_128 = mm512_add_hi_lo_epi16(vsum);
    const __m128i vsum_64  = _mm_add_epi32(_mm_cvtepi16_epi32(vsum_128),
                                          _mm_cvtepi16_epi32(_mm_srli_si128(vsum_128, 8)));
    return variance_final_from_32bit_sum_avx512(vsse, vsum_64, sse);
}

static INLINE int variance_final_2048_avx512(__m512i vsse, __m512i vsum, unsigned int *const sse) {
    vsum                   = sum_to_32bit_avx512(vsum);
    const __m128i vsum_128 = mm512_sum_to_128_epi32(vsum);
    return variance_final_from_32bit_sum_avx512(vsse, vsum_128, sse);
}

static INLINE void variance32_avx512(const uint8_t *src, const int src_stride, const uint8_t *ref,
                                     const int ref_stride, const int h, __m512i *const vsse,
                                     __m512i *const vsum) {
    *vsum = _mm512_setzero_si512();

    for (int i = 0; i < (h >> 1); i++) {
        const __m512i s = _mm512_inserti64x4(
            _mm512_castsi256_si512(_mm256_loadu_si256((__m256i const *)src)),
            _mm256_loadu_si256((__m256i const *)(src + src_stride)),
            1);
        const __m512i r = _mm512_inserti64x4(
            _mm512_castsi256_si512(_mm256_loadu_si256((__m256i const *)ref)),
            _mm256_loadu_si256((__m256i const *)(ref + ref_stride)),
            1);

        variance_kernel_avx512(s, r, vsse, vsum);

        src += (src_stride << 1);
        ref += (ref_stride << 1);
    }
}

static INLINE void variance64_avx512(const uint8_t *src, const int src_stride, const uint8_t *ref,
                                     const int ref_stride, const int h, __m512i *const vsse,
                                     __m512i *const vsum) {
    *vsum = _mm512_setzero_si512();

    for (int i = 0; i < h; i++) {
        const __m512i s = _mm512_loadu_si512((__m512i const *)(src));
        const __m512i r = _mm512_loadu_si512((__m512i const *)(ref));
        variance_kernel_avx512(s, r, vsse, vsum);

        src += src_stride;
        ref += ref_stride;
    }
}

static INLINE void variance128_avx512(const uint8_t *src, const int src_stride, const uint8_t *ref,
                                      const int ref_stride, const int h, __m512i *const vsse,
                                      __m512i *const vsum) {
    *vsum = _mm512_setzero_si512();

    for (int i = 0; i < h; i++) {
        __m512i s = _mm512_loadu_si512((__m512i const *)(src));
        __m512i r = _mm512_loadu_si512((__m512i const *)(ref));
        variance_kernel_avx512(s, r, vsse, vsum);

        s = _mm512_loadu_si512((__m512i const *)(src + 64));
        r = _mm512_loadu_si512((__m512i const *)(ref + 64));
        variance_kernel_avx512(s, r, vsse, vsum);

        src += src_stride;
        ref += ref_stride;
    }
}

#define AOM_VAR_LOOP_AVX512(bw, bh, bits, uh)                                               \
    unsigned int svt_aom_variance##bw##x##bh##_avx512(const uint8_t *src,                   \
                                                      int            src_stride,            \
                                                      const uint8_t *ref,                   \
                                                      int            ref_stride,            \
                                                      unsigned int * sse) {                  \
        __m512i vsse = _mm512_setzero_si512();                                              \
        __m512i vsum = _mm512_setzero_si512();                                              \
        for (int i = 0; i < (bh / uh); i++) {                                               \
            __m512i vsum16;                                                                 \
            variance##bw##_avx512(src, src_stride, ref, ref_stride, uh, &vsse, &vsum16);    \
            vsum = _mm512_add_epi32(vsum, sum_to_32bit_avx512(vsum16));                     \
            src += uh * src_stride;                                                         \
            ref += uh * ref_stride;                                                         \
        }                                                                                   \
        const __m128i vsum_128 = mm512_sum_to_128_epi32(vsum);                              \
        const int     sum      = variance_final_from_32bit_sum_avx512(vsse, vsum_128, sse); \
        return *sse - (unsigned int)(((int64_t)sum * sum) >> bits);                         \
    }

AOM_VAR_LOOP_AVX512(64, 64, 12, 32); // 64x32 * ( 64/32)
AOM_VAR_LOOP_AVX512(64, 128, 13, 32); // 64x32 * (128/32)
AOM_VAR_LOOP_AVX512(128, 64, 13, 16); // 128x16 * ( 64/16)
AOM_VAR_LOOP_AVX512(128, 128, 14, 16); // 128x16 * (128/16)

#define AOM_VAR_NO_LOOP_AVX512(bw, bh, bits, max_pixel)                            \
    unsigned int svt_aom_variance##bw##x##bh##_avx512(const uint8_t *src,          \
                                                      int            src_stride,   \
                                                      const uint8_t *ref,          \
                                                      int            ref_stride,   \
                                                      unsigned int * sse) {         \
        __m512i vsse = _mm512_setzero_si512();                                     \
        __m512i vsum = _mm512_setzero_si512();                                     \
        variance##bw##_avx512(src, src_stride, ref, ref_stride, bh, &vsse, &vsum); \
        const int sum = variance_final_##max_pixel##_avx512(vsse, vsum, sse);      \
        return *sse - (uint32_t)(((int64_t)sum * sum) >> bits);                    \
    }

AOM_VAR_NO_LOOP_AVX512(32, 8, 8, 512);
AOM_VAR_NO_LOOP_AVX512(32, 16, 9, 512);
AOM_VAR_NO_LOOP_AVX512(32, 32, 10, 1024);
AOM_VAR_NO_LOOP_AVX512(32, 64, 11, 2048);

AOM_VAR_NO_LOOP_AVX512(64, 16, 10, 1024);
AOM_VAR_NO_LOOP_AVX512(64, 32, 11, 2048);

#endif // VARIANCE_AVX512

#endif // !EN_AVX512_SUPPORT
