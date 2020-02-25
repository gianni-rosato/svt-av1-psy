/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"

#include <emmintrin.h>
#include <immintrin.h>

static INLINE void energy_computation_kernel_avx2(const int32_t *const in, __m256i *const sum256) {
    const __m256i zero      = _mm256_setzero_si256();
    const __m256i input     = _mm256_loadu_si256((__m256i *)in);
    const __m256i in_lo     = _mm256_unpacklo_epi32(input, zero);
    const __m256i in_hi     = _mm256_unpackhi_epi32(input, zero);
    const __m256i energy_lo = _mm256_mul_epi32(in_lo, in_lo);
    const __m256i energy_hi = _mm256_mul_epi32(in_hi, in_hi);
    *sum256                 = _mm256_add_epi64(*sum256, energy_lo);
    *sum256                 = _mm256_add_epi64(*sum256, energy_hi);
}

static INLINE uint64_t hadd64_avx2(const __m256i sum256) {
    const __m128i sum256_lo = _mm256_castsi256_si128(sum256);
    const __m128i sum256_hi = _mm256_extracti128_si256(sum256, 1);
    const __m128i sum128    = _mm_add_epi64(sum256_lo, sum256_hi);
    const __m128i sum128_hi = _mm_srli_si128(sum128, 8);
    const __m128i sum       = _mm_add_epi64(sum128, sum128_hi);

    return _mm_extract_epi64(sum, 0);
}

static INLINE uint64_t energy_computation_avx2(const int32_t *const in, const uint32_t size) {
    const __m256i zero = _mm256_setzero_si256();
    uint32_t      i    = 0;
    __m256i       sum  = zero;

    do {
        energy_computation_kernel_avx2(in + i, &sum);
        i += 8;
    } while (i < size);

    return hadd64_avx2(sum);
}

static INLINE uint64_t energy_computation_64_avx2(const int32_t *in, const uint32_t height) {
    const __m256i zero = _mm256_setzero_si256();
    uint32_t      i    = height;
    __m256i       sum  = zero;

    do {
        energy_computation_kernel_avx2(in + 0 * 8, &sum);
        energy_computation_kernel_avx2(in + 1 * 8, &sum);
        energy_computation_kernel_avx2(in + 2 * 8, &sum);
        energy_computation_kernel_avx2(in + 3 * 8, &sum);
        in += 64;
    } while (--i);

    return hadd64_avx2(sum);
}

static INLINE void clean_256_bytes_avx2(int32_t *buf, const uint32_t height) {
    const __m256i zero = _mm256_setzero_si256();
    uint32_t      h    = height;

    do {
        _mm256_storeu_si256((__m256i *)(buf + 0 * 8), zero);
        _mm256_storeu_si256((__m256i *)(buf + 1 * 8), zero);
        _mm256_storeu_si256((__m256i *)(buf + 2 * 8), zero);
        _mm256_storeu_si256((__m256i *)(buf + 3 * 8), zero);
        buf += 64;
    } while (--h);
}

static INLINE void copy_32_bytes_avx2(const int32_t *src, int32_t *dst) {
    const __m256i val = _mm256_loadu_si256((__m256i *)(src + 0 * 8));
    _mm256_storeu_si256((__m256i *)(dst + 0 * 8), val);
}

static INLINE void copy_256x_bytes_avx2(const int32_t *src, int32_t *dst, const uint32_t height) {
    uint32_t h = height;

    do {
        copy_32_bytes_avx2(src + 0 * 8, dst + 0 * 8);
        copy_32_bytes_avx2(src + 1 * 8, dst + 1 * 8);
        copy_32_bytes_avx2(src + 2 * 8, dst + 2 * 8);
        copy_32_bytes_avx2(src + 3 * 8, dst + 3 * 8);
        src += 64;
        dst += 32;
    } while (--h);
}

uint64_t handle_transform16x64_avx2(int32_t *output) {
    //bottom 16x32 area.
    const uint64_t three_quad_energy = energy_computation_avx2(output + 16 * 32, 16 * 32);

    // zero out the bottom 16x32 area.
    memset(output + 16 * 32, 0, 16 * 32 * sizeof(*output));

    return three_quad_energy;
}

uint64_t handle_transform32x64_avx2(int32_t *output) {
    //bottom 32x32 area.
    const uint64_t three_quad_energy = energy_computation_avx2(output + 32 * 32, 32 * 32);

    // zero out the bottom 32x32 area.
    memset(output + 32 * 32, 0, 32 * 32 * sizeof(*output));

    return three_quad_energy;
}

uint64_t handle_transform64x16_avx2(int32_t *output) {
    // top - right 32x16 area.
    const uint64_t three_quad_energy = energy_computation_64_avx2(output + 32, 16);

    // zero out right 32x16 area.
    clean_256_bytes_avx2(output + 32, 16);

    // Re-pack non-zero coeffs in the first 32x16 indices.
    copy_256x_bytes_avx2(output + 64, output + 32, 15);

    return three_quad_energy;
}

uint64_t handle_transform64x32_avx2(int32_t *output) {
    // top - right 32x32 area.
    const uint64_t three_quad_energy = energy_computation_64_avx2(output + 32, 32);

    // zero out right 32x32 area.
    clean_256_bytes_avx2(output + 32, 32);

    // Re-pack non-zero coeffs in the first 32x32 indices.
    copy_256x_bytes_avx2(output + 64, output + 32, 31);

    return three_quad_energy;
}

uint64_t handle_transform64x64_avx2(int32_t *output) {
    uint64_t three_quad_energy;

    // top - right 32x32 area.
    three_quad_energy = energy_computation_64_avx2(output + 32, 32);
    //bottom 64x32 area.
    three_quad_energy += energy_computation_avx2(output + 32 * 64, 64 * 32);

    // zero out top-right 32x32 area.
    clean_256_bytes_avx2(output + 32, 32);

    // zero out the bottom 64x32 area.
    memset(output + 32 * 64, 0, 32 * 64 * sizeof(*output));

    // Re-pack non-zero coeffs in the first 32x32 indices.
    copy_256x_bytes_avx2(output + 64, output + 32, 31);

    return three_quad_energy;
}
