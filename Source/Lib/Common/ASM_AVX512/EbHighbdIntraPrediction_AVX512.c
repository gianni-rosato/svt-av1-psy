/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#ifndef NON_AVX512_SUPPORT
#include <immintrin.h>
#include "EbHighbdIntraPrediction_SSE2.h"
#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"

// =============================================================================

// DC RELATED PRED

// Handle number of elements: up to 64.
static INLINE __m128i dc_sum_large(const __m256i src) {
    const __m128i s_lo = _mm256_extracti128_si256(src, 0);
    const __m128i s_hi = _mm256_extracti128_si256(src, 1);
    __m128i sum, sum_hi;
    sum = _mm_add_epi16(s_lo, s_hi);
    sum_hi = _mm_srli_si128(sum, 8);
    sum = _mm_add_epi16(sum, sum_hi);
    // Unpack to avoid 12-bit overflow.
    sum = _mm_unpacklo_epi16(sum, _mm_setzero_si128());

    return dc_sum_4x32bit(sum);
}

static INLINE void dc_common_predictor_32xh_kernel_avx512(uint16_t *dst,
    const ptrdiff_t stride, const int32_t h, const __m512i dc) {
    for (int32_t i = 0; i < h; i++) {
        _mm512_storeu_si512((__m512i *)dst, dc);
        dst += stride;
    }
}

static INLINE void dc_common_predictor_32xh(uint16_t *const dst,
    const ptrdiff_t stride, const int32_t h, const __m128i dc) {
    const __m512i expected_dc = _mm512_broadcastw_epi16(dc);
    dc_common_predictor_32xh_kernel_avx512(dst, stride, h, expected_dc);
}

static INLINE void dc_common_predictor_64xh_kernel_avx512(uint16_t *dst,
    const ptrdiff_t stride, const int32_t h, const __m512i dc) {
    for (int32_t i = 0; i < h; i++) {
        _mm512_storeu_si512((__m512i *)(dst + 0x00), dc);
        _mm512_storeu_si512((__m512i *)(dst + 0x20), dc);
        dst += stride;
    }
}

static INLINE void dc_common_predictor_64xh(uint16_t *const dst,
    const ptrdiff_t stride, const int32_t h, const __m128i dc) {
    const __m512i expected_dc = _mm512_broadcastw_epi16(dc);
    dc_common_predictor_64xh_kernel_avx512(dst, stride, h, expected_dc);
}

static INLINE __m128i dc_sum_16(const uint16_t *const src) {
    const __m256i s = _mm256_loadu_si256((const __m256i *) src);
    const __m128i s_lo = _mm256_extracti128_si256(s, 0);
    const __m128i s_hi = _mm256_extracti128_si256(s, 1);
    const __m128i sum = _mm_add_epi16(s_lo, s_hi);
    return dc_sum_8x16bit(sum);
}

static INLINE __m128i dc_sum_32(const uint16_t *const src) {
    const __m512i s32 = _mm512_loadu_si512((const __m512i *) src);
    const __m256i s0 = _mm512_extracti64x4_epi64(s32, 0);
    const __m256i s1 = _mm512_extracti64x4_epi64(s32, 1);
    const __m256i sum = _mm256_add_epi16(s0, s1);
    return dc_sum_large(sum);
}

static INLINE __m128i dc_sum_64(const uint16_t *const src) {
    const __m512i s0 = _mm512_loadu_si512((const __m512i *)(src + 0x00));
    const __m512i s1 = _mm512_loadu_si512((const __m512i *)(src + 0x20));
    const __m512i s01 = _mm512_add_epi16(s0, s1);

    const __m256i s2 = _mm512_extracti64x4_epi64(s01, 0);
    const __m256i s3 = _mm512_extracti64x4_epi64(s01, 1);

    const __m256i sum = _mm256_add_epi16(s2, s3);
    return dc_sum_large(sum);
}

// 32xN

void aom_highbd_dc_left_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(4);
    __m128i sum;
    (void)above;
    (void)bd;

    sum = dc_sum_8(left);
    sum = _mm_add_epi16(sum, round);
    sum = _mm_srli_epi16(sum, 3);
    dc_common_predictor_32xh(dst, stride, 8, sum);
}

void aom_highbd_dc_left_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(8);
    __m128i sum;
    (void)above;
    (void)bd;

    sum = dc_sum_16(left);
    sum = _mm_add_epi16(sum, round);
    sum = _mm_srli_epi16(sum, 4);
    dc_common_predictor_32xh(dst, stride, 16, sum);
}

void aom_highbd_dc_left_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(16);
    __m128i sum;
    (void)above;
    (void)bd;

    sum = dc_sum_32(left);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 5);
    dc_common_predictor_32xh(dst, stride, 32, sum);
}

void aom_highbd_dc_left_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(32);
    __m128i sum;
    (void)above;
    (void)bd;

    sum = dc_sum_64(left);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 6);
    dc_common_predictor_32xh(dst, stride, 64, sum);
}

// 64xN

void aom_highbd_dc_left_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(8);
    __m128i sum;
    (void)above;
    (void)bd;

    sum = dc_sum_16(left);
    sum = _mm_add_epi16(sum, round);
    sum = _mm_srli_epi16(sum, 4);
    dc_common_predictor_64xh(dst, stride, 16, sum);
}

void aom_highbd_dc_left_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(16);
    __m128i sum;
    (void)above;
    (void)bd;

    sum = dc_sum_32(left);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 5);
    dc_common_predictor_64xh(dst, stride, 32, sum);
}

void aom_highbd_dc_left_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    const __m128i round = _mm_cvtsi32_si128(32);
    __m128i sum;
    (void)above;
    (void)bd;

    sum = dc_sum_64(left);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 6);
    dc_common_predictor_64xh(dst, stride, 64, sum);
}

/*highbd dc top predictors */

// 32xN

static INLINE void dc_top_predictor_32xh(uint16_t *const dst,
    const ptrdiff_t stride, const uint16_t *const above,
    const int32_t h, const int32_t bd)
{
    const __m128i round = _mm_cvtsi32_si128(16);
    __m128i sum;
    (void)bd;

    sum = dc_sum_32(above);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 5);
    dc_common_predictor_32xh(dst, stride, h, sum);
}

void aom_highbd_dc_top_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)left;

    dc_top_predictor_32xh(dst, stride, above, 8, bd);
}

void aom_highbd_dc_top_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)left;
    dc_top_predictor_32xh(dst, stride, above, 16, bd);
}

void aom_highbd_dc_top_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)left;
    dc_top_predictor_32xh(dst, stride, above, 32, bd);
}

void aom_highbd_dc_top_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)left;
    dc_top_predictor_32xh(dst, stride, above, 64, bd);
}

// 64xN

static INLINE void dc_top_predictor_64xh(uint16_t *const dst,
    const ptrdiff_t stride, const uint16_t *const above,
    const int32_t h, const int32_t bd)
{
    const __m128i round = _mm_cvtsi32_si128(32);
    __m128i sum;
    (void)bd;

    sum = dc_sum_64(above);
    sum = _mm_add_epi32(sum, round);
    sum = _mm_srli_epi32(sum, 6);
    dc_common_predictor_64xh(dst, stride, h, sum);
}

void aom_highbd_dc_top_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)left;
    dc_top_predictor_64xh(dst, stride, above, 16, bd);
}

void aom_highbd_dc_top_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)left;
    dc_top_predictor_64xh(dst, stride, above, 32, bd);
}

void aom_highbd_dc_top_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)left;
    dc_top_predictor_64xh(dst, stride, above, 64, bd);
}

/* highbd dc predictor */

// 32xN

static INLINE __m128i dc_sum_8_32(const uint16_t *const src_8,
    const uint16_t *const src_32) {
    const __m128i s_8 = _mm_loadu_si128((const __m128i *)src_8);
    const __m512i s32_01 = _mm512_loadu_si512((const __m512i *)(src_32 + 0x00));
    const __m256i s_32_0 = _mm512_extracti64x4_epi64(s32_01,0);
    const __m256i s_32_1 = _mm512_extracti64x4_epi64(s32_01,1);
    const __m256i s_32 = _mm256_add_epi16(s_32_0, s_32_1);
    const __m128i s_lo = _mm256_extracti128_si256(s_32, 0);
    const __m128i s_hi = _mm256_extracti128_si256(s_32, 1);
    const __m128i s_16_sum = _mm_add_epi16(s_lo, s_hi);
    const __m128i sum = _mm_add_epi16(s_8, s_16_sum);
    return dc_sum_8x16bit_large(sum);
}

static INLINE __m128i dc_sum_16_32(const uint16_t *const src_16,
    const uint16_t *const src_32) {
    const __m256i s_16 = _mm256_loadu_si256((const __m256i *)src_16);
    const __m512i s32_01 = _mm512_loadu_si512((const __m512i *)(src_32 + 0x00));
    const __m256i s_32_0 = _mm512_extracti64x4_epi64(s32_01, 0);
    const __m256i s_32_1 = _mm512_extracti64x4_epi64(s32_01, 1);
    const __m256i sum0 = _mm256_add_epi16(s_16, s_32_0);
    const __m256i sum = _mm256_add_epi16(sum0, s_32_1);
    return dc_sum_large(sum);
}

// Handle number of elements: 65 to 128.
static INLINE __m128i dc_sum_larger(const __m256i src) {
    const __m128i s_lo = _mm256_extracti128_si256(src, 0);
    const __m128i s_hi = _mm256_extracti128_si256(src, 1);
    __m128i sum, sum_hi;
    sum = _mm_add_epi16(s_lo, s_hi);
    // Unpack to avoid 12-bit overflow.
    sum_hi = _mm_unpackhi_epi16(sum, _mm_setzero_si128());
    sum = _mm_unpacklo_epi16(sum, _mm_setzero_si128());
    sum = _mm_add_epi32(sum, sum_hi);

    return dc_sum_4x32bit(sum);
}

static INLINE __m128i dc_sum_32_32(const uint16_t *const src0,
    const uint16_t *const src1) {
    const __m512i s_32_0 = _mm512_loadu_si512((const __m512i *)(src0 + 0x00));
    const __m512i s_32_1 = _mm512_loadu_si512((const __m512i *)(src1 + 0x00));
    const __m512i sum_32_01 = _mm512_add_epi16(s_32_0, s_32_1);
    const __m256i sum_16_0 = _mm512_extracti64x4_epi64(sum_32_01,0);
    const __m256i sum_16_1 = _mm512_extracti64x4_epi64(sum_32_01,1);
    const __m256i sum = _mm256_add_epi16(sum_16_0, sum_16_1);
    return dc_sum_large(sum);
}

static INLINE __m128i dc_sum_16_64(const uint16_t *const src_16,
    const uint16_t *const src_64)
{
    const __m256i s_16 = _mm256_loadu_si256((const __m256i *)src_16);
    const __m512i s_64_0 = _mm512_loadu_si512((const __m512i *)(src_64 + 0x00));
    const __m512i s_64_1 = _mm512_loadu_si512((const __m512i *)(src_64 + 0x20));
    const __m512i sum_64_01 = _mm512_add_epi16(s_64_1, s_64_0);
    const __m256i s1 = _mm512_extracti64x4_epi64(sum_64_01, 0);
    const __m256i s2 = _mm512_extracti64x4_epi64(sum_64_01, 1);
    const __m256i s3 = _mm256_add_epi16(s1, s_16);
    const __m256i sum = _mm256_add_epi16(s2, s3);
    return dc_sum_larger(sum);
}


static INLINE __m128i dc_sum_32_64(const uint16_t *const src_32,
    const uint16_t *const src_64) {
    const __m512i s_32_0 = _mm512_loadu_si512((const __m512i *)(src_32 + 0x00));
    const __m512i s_64_0 = _mm512_loadu_si512((const __m256i *)(src_64 + 0x00));
    const __m512i s_64_1 = _mm512_loadu_si512((const __m256i *)(src_64 + 0x20));

    const __m512i sum0 = _mm512_add_epi16(s_32_0, s_64_0);
    const __m512i sum1 = _mm512_add_epi16(sum0, s_64_1);

    const __m256i sum2 = _mm512_extracti64x4_epi64(sum1, 0);
    const __m256i sum3 = _mm512_extracti64x4_epi64(sum1, 1);
    const __m256i sum = _mm256_add_epi16(sum2, sum3);
    return dc_sum_larger(sum);
}

static INLINE __m128i dc_sum_64_64(const uint16_t *const src0,
    const uint16_t *const src1) {
    const __m512i s0 = _mm512_loadu_si512((const __m512i *)(src0 + 0x00));
    const __m512i s1 = _mm512_loadu_si512((const __m256i *)(src0 + 0x20));
    const __m512i s2 = _mm512_loadu_si512((const __m256i *)(src1 + 0x00));
    const __m512i s3 = _mm512_loadu_si512((const __m256i *)(src1 + 0x20));

    const __m512i sum01 = _mm512_add_epi16(s0, s1);
    const __m512i sum23 = _mm512_add_epi16(s2, s3);
    const __m512i sum03 = _mm512_add_epi16(sum01, sum23);

    const __m256i sum03_1 = _mm512_extracti64x4_epi64(sum03, 0);
    const __m256i sum03_2 = _mm512_extracti64x4_epi64(sum03, 1);

    const __m256i sum = _mm256_add_epi16(sum03_1, sum03_2);
    return dc_sum_larger(sum);
}

void aom_highbd_dc_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    (void)bd;
    __m128i sum = dc_sum_8_32(left, above);
    uint32_t sum32 = _mm_cvtsi128_si32(sum);
    sum32 += 20;
    sum32 /= 40;
    const __m512i dc = _mm512_set1_epi16((int16_t)sum32);

    dc_common_predictor_32xh_kernel_avx512(dst, stride, 8, dc);
}

void aom_highbd_dc_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    (void)bd;
    __m128i sum = dc_sum_16_32(left, above);
    uint32_t sum32 = _mm_cvtsi128_si32(sum);
    sum32 += 24;
    sum32 /= 48;
    const __m512i dc = _mm512_set1_epi16((int16_t)sum32);

    dc_common_predictor_32xh_kernel_avx512(dst, stride, 16, dc);
}

void aom_highbd_dc_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    (void)bd;
    __m128i sum = dc_sum_32_32(above, left);
    sum = _mm_add_epi32(sum, _mm_set1_epi32(32));
    sum = _mm_srli_epi32(sum, 6);
    dc_common_predictor_32xh(dst, stride, 32, sum);
}

void aom_highbd_dc_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    (void)bd;
    __m128i sum = dc_sum_32_64(above, left);
    uint32_t sum32 = _mm_cvtsi128_si32(sum);
    sum32 += 48;
    sum32 /= 96;
    const __m512i dc = _mm512_set1_epi16((int16_t)sum32);

    dc_common_predictor_32xh_kernel_avx512(dst, stride, 64, dc);
}

// 64xN

void aom_highbd_dc_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    (void)bd;
    __m128i sum = dc_sum_16_64(left, above);
    uint32_t sum32 = _mm_cvtsi128_si32(sum);
    sum32 += 40;
    sum32 /= 80;
    const __m512i dc = _mm512_set1_epi16((int16_t)sum32);

    dc_common_predictor_64xh_kernel_avx512(dst, stride, 16, dc);
}

void aom_highbd_dc_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    (void)bd;
    __m128i sum = dc_sum_32_64(left, above);
    uint32_t sum32 = _mm_cvtsi128_si32(sum);
    sum32 += 48;
    sum32 /= 96;
    const __m512i dc = _mm512_set1_epi16((int16_t)sum32);

    dc_common_predictor_64xh_kernel_avx512(dst, stride, 32, dc);
}

void aom_highbd_dc_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd) {
    (void)bd;
    __m128i sum = dc_sum_64_64(above, left);
    sum = _mm_add_epi32(sum, _mm_set1_epi32(64));
    sum = _mm_srli_epi32(sum, 7);
    dc_common_predictor_64xh(dst, stride, 64, sum);
}
#endif
