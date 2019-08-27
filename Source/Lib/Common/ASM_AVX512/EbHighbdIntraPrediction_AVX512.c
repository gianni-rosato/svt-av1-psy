/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include <immintrin.h>
#include "EbHighbdIntraPrediction_SSE2.h"
#include "EbDefinitions.h"
#include "aom_dsp_rtcd.h"

#ifndef NON_AVX512_SUPPORT
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

// =============================================================================

// H_PRED

// -----------------------------------------------------------------------------

// 32xN

static INLINE void h_pred_32(uint16_t **const dst, const ptrdiff_t stride,
    const __m128i left)
{
    // Broadcast the 16-bit left pixel to 256-bit register.
    const __m512i row = _mm512_broadcastw_epi16(left);

    _mm512_storeu_si512((__m512i *)(*dst + 0x00), row);
    *dst += stride;
}

// Process 8 rows.
static INLINE void h_pred_32x8(uint16_t **dst, const ptrdiff_t stride,
    const uint16_t *const left)
{
    // dst and it's stride must be 32-byte aligned.
    assert(!((intptr_t)*dst % 32));
    assert(!(stride % 32));

    const __m128i left_u16 = _mm_load_si128((const __m128i *)left);

    h_pred_32(dst, stride, _mm_srli_si128(left_u16, 0));
    h_pred_32(dst, stride, _mm_srli_si128(left_u16, 2));
    h_pred_32(dst, stride, _mm_srli_si128(left_u16, 4));
    h_pred_32(dst, stride, _mm_srli_si128(left_u16, 6));
    h_pred_32(dst, stride, _mm_srli_si128(left_u16, 8));
    h_pred_32(dst, stride, _mm_srli_si128(left_u16, 10));
    h_pred_32(dst, stride, _mm_srli_si128(left_u16, 12));
    h_pred_32(dst, stride, _mm_srli_si128(left_u16, 14));
}

// 32x8

void aom_highbd_h_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)above;
    (void)bd;

    h_pred_32x8(&dst, stride, left);
}

// 32x64

void aom_highbd_h_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)above;
    (void)bd;

    for (int32_t i = 0; i < 8; i++, left += 8) {
        h_pred_32x8(&dst, stride, left);
    }
}

//32x16 32x32
static INLINE void h_store_32_unpacklo(uint16_t **dst, const ptrdiff_t stride,
    const __m128i *row) {
    const __m128i val = _mm_unpacklo_epi64(*row, *row);
    const __m512i val2 = _mm512_broadcast_i32x4(val);
    _mm512_storeu_si512((__m512i *)(*dst), val2);
    *dst += stride;
}

static INLINE void h_store_32_unpackhi(uint16_t **dst, const ptrdiff_t stride,
    const __m128i *row) {
    const __m128i val = _mm_unpackhi_epi64(*row, *row);
    const __m512i val2 = _mm512_broadcast_i32x4(val);
    _mm512_storeu_si512((__m512i *)(*dst), val2);
    *dst += stride;
}

static INLINE void h_predictor_32x8(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left) {
    const __m128i left_u16 = _mm_load_si128((const __m128i *)left);
    const __m128i row0 = _mm_shufflelo_epi16(left_u16, 0x0);
    const __m128i row1 = _mm_shufflelo_epi16(left_u16, 0x55);
    const __m128i row2 = _mm_shufflelo_epi16(left_u16, 0xaa);
    const __m128i row3 = _mm_shufflelo_epi16(left_u16, 0xff);
    const __m128i row4 = _mm_shufflehi_epi16(left_u16, 0x0);
    const __m128i row5 = _mm_shufflehi_epi16(left_u16, 0x55);
    const __m128i row6 = _mm_shufflehi_epi16(left_u16, 0xaa);
    const __m128i row7 = _mm_shufflehi_epi16(left_u16, 0xff);
    h_store_32_unpacklo(&dst, stride, &row0);
    h_store_32_unpacklo(&dst, stride, &row1);
    h_store_32_unpacklo(&dst, stride, &row2);
    h_store_32_unpacklo(&dst, stride, &row3);
    h_store_32_unpackhi(&dst, stride, &row4);
    h_store_32_unpackhi(&dst, stride, &row5);
    h_store_32_unpackhi(&dst, stride, &row6);
    h_store_32_unpackhi(&dst, stride, &row7);
}

void aom_highbd_h_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    int32_t i;
    (void)above;
    (void)bd;

    for (i = 0; i < 2; i++, left += 8) {
        h_predictor_32x8(dst, stride, left);
        dst += stride << 3;
    }
}

void aom_highbd_h_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    int32_t i;
    (void)above;
    (void)bd;

    for (i = 0; i < 4; i++, left += 8) {
        h_predictor_32x8(dst, stride, left);
        dst += stride << 3;
    }
}

// ---------------------------------------------------------------------------- -

// 64xN

static INLINE void h_pred_64(uint16_t **const dst, const ptrdiff_t stride,
    const __m128i left)
{
    // Broadcast the 16-bit left pixel to 256-bit register.
    const __m512i row = _mm512_broadcastw_epi16(left);

    _mm512_store_si512((__m256i *)(*dst + 0x00), row);
    _mm512_store_si512((__m256i *)(*dst + 0x20), row);

    *dst += stride;
}

// Process 8 rows.
static INLINE void h_pred_64x8(uint16_t **dst, const ptrdiff_t stride,
    const uint16_t *const left)
{
    // dst and it's stride must be 32-byte aligned.
    assert(!((intptr_t)*dst % 32));
    assert(!(stride % 32));

    __m128i left_u16 = _mm_load_si128((const __m128i *)left);

    for (int16_t j = 0; j < 8; j++) {
        h_pred_64(dst, stride, left_u16);
        left_u16 = _mm_srli_si128(left_u16, 2);
    }
}

// 64x16

void aom_highbd_h_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)above;
    (void)bd;

    for (int32_t i = 0; i < 2; i++, left += 8) {
        h_pred_64x8(&dst, stride, left);
    }
}

// 64x32

void aom_highbd_h_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)above;
    (void)bd;

    for (int32_t i = 0; i < 4; i++, left += 8) {
        h_pred_64x8(&dst, stride, left);
    }
}

// 64x64

void aom_highbd_h_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    (void)above;
    (void)bd;

    for (int32_t i = 0; i < 8; i++, left += 8) {
        h_pred_64x8(&dst, stride, left);
    }
}

// =============================================================================

// V_PRED

// -----------------------------------------------------------------------------

// 32xN

static INLINE void v_pred_32(uint16_t **const dst, const ptrdiff_t stride,
    const __m512i above01)
{
    _mm512_storeu_si512((__m512i *)(*dst + 0x00), above01);
    *dst += stride;
}

// Process 8 rows.
static INLINE void v_pred_32x8(uint16_t **const dst, const ptrdiff_t stride,
    const __m512i above01)
{
    // dst and it's stride must be 32-byte aligned.
    assert(!((intptr_t)*dst % 32));
    assert(!(stride % 32));

    v_pred_32(dst, stride, above01);
    v_pred_32(dst, stride, above01);
    v_pred_32(dst, stride, above01);
    v_pred_32(dst, stride, above01);
    v_pred_32(dst, stride, above01);
    v_pred_32(dst, stride, above01);
    v_pred_32(dst, stride, above01);
    v_pred_32(dst, stride, above01);
}

// 32x8

void aom_highbd_v_predictor_32x8_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    // Load all 32 pixels in a row into 256-bit registers.
    const __m512i above01 = _mm512_loadu_si512((const __m512i *)(above + 0x00));

    (void)left;
    (void)bd;

    v_pred_32x8(&dst, stride, above01);
}

// 32x16

void aom_highbd_v_predictor_32x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    // Load all 32 pixels in a row into 512-bit registers.
    const __m512i above01 = _mm512_loadu_si512((const __m512i *)(above + 0x00));

    (void)left;
    (void)bd;

    for (int32_t i = 0; i < 2; i++) {
        v_pred_32x8(&dst, stride, above01);
    }
}

// 32x32

void aom_highbd_v_predictor_32x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    // Load all 32 pixels in a row into 512-bit registers.
    const __m512i above01 = _mm512_loadu_si512((const __m512i *)(above + 0x00));

    (void)left;
    (void)bd;

    for (int32_t i = 0; i < 4; i++) {
        v_pred_32x8(&dst, stride, above01);
    }
}

// 32x64

void aom_highbd_v_predictor_32x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    // Load all 32 pixels in a row into 512-bit registers.
    const __m512i above01 = _mm512_loadu_si512((const __m512i *)(above + 0x00));

    (void)left;
    (void)bd;

    for (int32_t i = 0; i < 8; i++) {
        v_pred_32x8(&dst, stride, above01);
    }
}

// -----------------------------------------------------------------------------

// 64xN

static INLINE void v_pred_64(uint16_t **const dst, const ptrdiff_t stride,
    const __m512i above0, const __m512i above1)
{
    _mm512_storeu_si512((__m512i *)(*dst + 0x00), above0);
    _mm512_storeu_si512((__m512i *)(*dst + 0x20), above1);
    *dst += stride;
}

// Process 8 rows.
static INLINE void v_pred_64x8(uint16_t **const dst, const ptrdiff_t stride,
    const __m512i above0, const __m512i above1)
{
    // dst and it's stride must be 32-byte aligned.
    assert(!((intptr_t)*dst % 32));
    assert(!(stride % 32));

    v_pred_64(dst, stride, above0, above1);
    v_pred_64(dst, stride, above0, above1);
    v_pred_64(dst, stride, above0, above1);
    v_pred_64(dst, stride, above0, above1);
    v_pred_64(dst, stride, above0, above1);
    v_pred_64(dst, stride, above0, above1);
    v_pred_64(dst, stride, above0, above1);
    v_pred_64(dst, stride, above0, above1);
}

// 64x16

void aom_highbd_v_predictor_64x16_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    // Load all 64 pixels in a row into 512-bit registers.
    const __m512i above0 = _mm512_loadu_si512((const __m512i *)(above + 0x00));
    const __m512i above1 = _mm512_loadu_si512((const __m512i *)(above + 0x20));

    (void)left;
    (void)bd;

    for (int32_t i = 0; i < 2; i++) {
        v_pred_64x8(&dst, stride, above0, above1);
    }
}

// 64x32

void aom_highbd_v_predictor_64x32_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    // Load all 64 pixels in a row into 512-bit registers.
    const __m512i above0 = _mm512_loadu_si512((const __m512i *)(above + 0x00));
    const __m512i above1 = _mm512_loadu_si512((const __m512i *)(above + 0x20));

    (void)left;
    (void)bd;

    for (int32_t i = 0; i < 4; i++) {
        v_pred_64x8(&dst, stride, above0, above1);
    }
}

// 64x64

void aom_highbd_v_predictor_64x64_avx512(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *above, const uint16_t *left, int32_t bd)
{
    // Load all 64 pixels in a row into 512-bit registers.
    const __m512i above0 = _mm512_loadu_si512((const __m512i *)(above + 0x00));
    const __m512i above1 = _mm512_loadu_si512((const __m512i *)(above + 0x20));

    (void)left;
    (void)bd;

    for (int32_t i = 0; i < 8; i++) {
        v_pred_64x8(&dst, stride, above0, above1);
    }
}
#endif
