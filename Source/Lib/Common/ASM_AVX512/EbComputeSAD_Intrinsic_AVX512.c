/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <assert.h>

#include "EbComputeSAD_AVX2.h"
#include "EbDefinitions.h"
#include "immintrin.h"
#include "EbMemory_AVX2.h"
#include "transpose_avx2.h"

#ifndef NON_AVX512_SUPPORT

static INLINE void sad64_kernel_avx512(const __m512i s,
    const uint8_t *const ref, __m512i *const sum)
{
    const __m512i r = _mm512_loadu_si512((__m512i*)ref);
    *sum = _mm512_add_epi32(*sum, _mm512_sad_epu8(s, r));
}

static INLINE void sad64_avx512(const uint8_t *const src, const uint8_t *ref,
    __m512i *const sum)
{
    const __m512i s = _mm512_loadu_si512((__m512i*)src);
    sad64_kernel_avx512(s, ref, sum);
}

static INLINE uint32_t sad_final_avx512(const __m512i zmm) {
    const __m256i zmm_L = _mm512_castsi512_si256(zmm);
    const __m256i zmm_H = _mm512_extracti64x4_epi64(zmm, 1);
    const __m256i ymm = _mm256_add_epi32(zmm_L, zmm_H);
    const __m128i ymm_L = _mm256_castsi256_si128(ymm);
    const __m128i ymm_H = _mm256_extracti128_si256(ymm, 1);
    const __m128i xmm0 = _mm_add_epi32(ymm_L, ymm_H);
    const __m128i xmm0_H = _mm_srli_si128(xmm0, 8);
    const __m128i xmm1 = _mm_add_epi32(xmm0, xmm0_H);

    return _mm_extract_epi32(xmm1, 0);
}

/*******************************************************************************
* Requirement: height % 2 = 0
*******************************************************************************/
static AOM_FORCE_INLINE uint32_t compute64x_m_sad_avx512_intrin(
    const uint8_t *src, const uint32_t src_stride, const uint8_t *ref,
    const uint32_t ref_stride, const uint32_t height)
{
    uint32_t y = height;
    __m512i zmm = _mm512_setzero_si512();

    do {
        sad64_avx512(src + 0 * src_stride, ref + 0 * ref_stride, &zmm);
        sad64_avx512(src + 1 * src_stride, ref + 1 * ref_stride, &zmm);
        src += src_stride << 1;
        ref += ref_stride << 1;
        y -= 2;
    } while (y);

    return sad_final_avx512(zmm);
}

/*******************************************************************************
* Requirement: height % 2 = 0
*******************************************************************************/
static INLINE uint32_t compute128x_m_sad_avx512_intrin(
    const uint8_t  *src,         // input parameter, source samples Ptr
    const uint32_t  src_stride,  // input parameter, source stride
    const uint8_t  *ref,         // input parameter, reference samples Ptr
    const uint32_t  ref_stride,  // input parameter, reference stride
    const uint32_t  height)      // input parameter, block height (M)
{
    uint32_t y = height;
    __m512i zmm = _mm512_setzero_si512();

    do {
        sad64_avx512(src + 0 * src_stride + 0 * 64, ref + 0 * ref_stride + 0 * 64, &zmm);
        sad64_avx512(src + 0 * src_stride + 1 * 64, ref + 0 * ref_stride + 1 * 64, &zmm);
        sad64_avx512(src + 1 * src_stride + 0 * 64, ref + 1 * ref_stride + 0 * 64, &zmm);
        sad64_avx512(src + 1 * src_stride + 1 * 64, ref + 1 * ref_stride + 1 * 64, &zmm);
        src += src_stride << 1;
        ref += ref_stride << 1;
        y -= 2;
    } while (y);

    return sad_final_avx512(zmm);
}

uint32_t eb_aom_sad64x16_avx512(const uint8_t *src_ptr, int src_stride,
    const uint8_t *ref_ptr, int ref_stride) {
    return compute64x_m_sad_avx512_intrin(src_ptr, src_stride, ref_ptr,
        ref_stride, 16);
}

uint32_t eb_aom_sad64x32_avx512(const uint8_t *src_ptr, int src_stride,
    const uint8_t *ref_ptr, int ref_stride) {
    return compute64x_m_sad_avx512_intrin(src_ptr, src_stride, ref_ptr,
        ref_stride, 32);
}

uint32_t eb_aom_sad64x64_avx512(const uint8_t *src_ptr, int src_stride,
    const uint8_t *ref_ptr, int ref_stride) {
    return compute64x_m_sad_avx512_intrin(src_ptr, src_stride, ref_ptr,
        ref_stride, 64);
}

uint32_t eb_aom_sad64x128_avx512(const uint8_t *src_ptr, int src_stride,
    const uint8_t *ref_ptr, int ref_stride) {
    return compute64x_m_sad_avx512_intrin(src_ptr, src_stride, ref_ptr,
        ref_stride, 128);
}

uint32_t eb_aom_sad128x64_avx512(const uint8_t *src_ptr, int src_stride,
    const uint8_t *ref_ptr, int ref_stride) {
    return compute128x_m_sad_avx512_intrin(src_ptr, src_stride, ref_ptr,
        ref_stride, 64);
}

uint32_t eb_aom_sad128x128_avx512(const uint8_t *src_ptr, int src_stride,
    const uint8_t *ref_ptr, int ref_stride) {
    return compute128x_m_sad_avx512_intrin(src_ptr, src_stride, ref_ptr,
        ref_stride, 128);
}

static INLINE void sad64_4d_avx512(const uint8_t *const src,
    const uint8_t *const ref_array[4], const uint32_t offset,
    __m512i sum[4])
{
    const __m512i s = _mm512_loadu_si512((__m512i*)(src + offset));
    sad64_kernel_avx512(s, ref_array[0] + offset, &sum[0]);
    sad64_kernel_avx512(s, ref_array[1] + offset, &sum[1]);
    sad64_kernel_avx512(s, ref_array[2] + offset, &sum[2]);
    sad64_kernel_avx512(s, ref_array[3] + offset, &sum[3]);
}

static INLINE __m256i add_hi_lo_32_avx512(const __m512i src) {
    const __m256i s0 = _mm512_castsi512_si256(src);
    const __m256i s1 = _mm512_extracti64x4_epi64(src, 1);
    return _mm256_add_epi32(s0, s1);
}

static INLINE __m128i hadd_four_32_avx2(const __m256i src0, const __m256i src1,
    const __m256i src2, const __m256i src3)
{
    const __m256i s01 = _mm256_hadd_epi32(src0, src1);       // 0 0 1 1  0 0 1 1
    const __m256i s23 = _mm256_hadd_epi32(src2, src3);       // 2 2 3 3  2 2 3 3
    const __m256i s0123 = _mm256_hadd_epi32(s01, s23);       // 0 1 2 3  0 1 2 3
    const __m128i sum0 = _mm256_castsi256_si128(s0123);      // 0 1 2 3
    const __m128i sum1 = _mm256_extracti128_si256(s0123, 1); // 0 1 2 3
    return _mm_add_epi32(sum0, sum1);                        // 0 1 2 3
}

static INLINE __m128i hadd_four_32_avx512(const __m512i src0, const __m512i src1,
    const __m512i src2, const __m512i src3)
{
    __m256i s[4];

    s[0] = add_hi_lo_32_avx512(src0);
    s[1] = add_hi_lo_32_avx512(src1);
    s[2] = add_hi_lo_32_avx512(src2);
    s[3] = add_hi_lo_32_avx512(src3);

    return hadd_four_32_avx2(s[0], s[1], s[2], s[3]);
}

#if 0 // For sad64xMx4d, AVX512 is Slower than AVX2 because of worse AVX512 compiler.

static AOM_FORCE_INLINE void compute64x_m_4d_sad_avx512_intrin(
    const uint8_t *src, const uint32_t src_stride,
    const uint8_t *const ref_array[4], const uint32_t ref_stride,
    uint32_t sad_array[4], const uint32_t height)
{
    const uint8_t *ref[4];
    uint32_t y = height;
    __m512i zmm[4] = { 0 };

    ref[0] = ref_array[0];
    ref[1] = ref_array[1];
    ref[2] = ref_array[2];
    ref[3] = ref_array[3];

    do {
        sad64_4d_avx512(src, ref, 0, zmm);
        src += src_stride;
        ref[0] += ref_stride;
        ref[1] += ref_stride;
        ref[2] += ref_stride;
        ref[3] += ref_stride;
    } while (--y);

    const __m128i sum = hadd_four_32_avx512(zmm[0], zmm[1], zmm[2], zmm[3]);
    _mm_storeu_si128((__m128i *)sad_array, sum);
}

void eb_aom_sad64x16x4d_avx512(const uint8_t *src, int src_stride,
    const uint8_t *const ref_array[4], int ref_stride,
    uint32_t sad_array[4]) {
    compute64x_m_4d_sad_avx512_intrin(src, src_stride, ref_array, ref_stride, sad_array, 16);
}

void eb_aom_sad64x32x4d_avx512(const uint8_t *src, int src_stride,
    const uint8_t *const ref_array[4], int ref_stride,
    uint32_t sad_array[4]) {
    compute64x_m_4d_sad_avx512_intrin(src, src_stride, ref_array, ref_stride, sad_array, 32);
}

void eb_aom_sad64x64x4d_avx512(const uint8_t *src, int src_stride,
    const uint8_t *const ref_array[4], int ref_stride,
    uint32_t sad_array[4]) {
    compute64x_m_4d_sad_avx512_intrin(src, src_stride, ref_array, ref_stride, sad_array, 64);
}

void eb_aom_sad64x128x4d_avx512(const uint8_t *src, int src_stride,
    const uint8_t *const ref_array[4], int ref_stride,
    uint32_t sad_array[4]) {
    compute64x_m_4d_sad_avx512_intrin(src, src_stride, ref_array, ref_stride, sad_array, 128);
}

#endif

static INLINE void compute128x_m_4d_sad_avx512_intrin(const uint8_t *src,
    const uint32_t src_stride, const uint8_t *const ref_array[4],
    const uint32_t ref_stride, uint32_t sad_array[4], const uint32_t height)
{
    const uint8_t *ref[4];
    uint32_t y = height;
    __m512i zmm[4] = { 0 };

    ref[0] = ref_array[0];
    ref[1] = ref_array[1];
    ref[2] = ref_array[2];
    ref[3] = ref_array[3];

    do {
        sad64_4d_avx512(src, ref, 0 * 64, zmm);
        sad64_4d_avx512(src, ref, 1 * 64, zmm);
        src += src_stride;
        ref[0] += ref_stride;
        ref[1] += ref_stride;
        ref[2] += ref_stride;
        ref[3] += ref_stride;
    } while (--y);

    const __m128i sum = hadd_four_32_avx512(zmm[0], zmm[1], zmm[2], zmm[3]);
    _mm_storeu_si128((__m128i *)sad_array, sum);
}

void eb_aom_sad128x64x4d_avx512(const uint8_t *src, int src_stride,
    const uint8_t *const ref_array[4], int ref_stride,
    uint32_t sad_array[4]) {
    compute128x_m_4d_sad_avx512_intrin(src, src_stride, ref_array, ref_stride, sad_array, 64);
}

void eb_aom_sad128x128x4d_avx512(const uint8_t *src, int src_stride,
    const uint8_t *const ref_array[4], int ref_stride,
    uint32_t sad_array[4]) {
    compute128x_m_4d_sad_avx512_intrin(src, src_stride, ref_array, ref_stride, sad_array, 128);
}

#endif // !NON_AVX512_SUPPORT
