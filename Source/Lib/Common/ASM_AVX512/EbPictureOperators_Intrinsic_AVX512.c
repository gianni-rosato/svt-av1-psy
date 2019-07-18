/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include <immintrin.h>
#include "EbPictureOperators_AVX2.h"
#include "EbPictureOperators_SSE2.h"
#include "EbMemory_AVX2.h"

#ifndef NON_AVX512_SUPPORT

static INLINE int32_t Hadd32_AVX512_INTRIN(const __m512i src) {
    const __m256i src_L = _mm512_extracti64x4_epi64(src, 0);
    const __m256i src_H = _mm512_extracti64x4_epi64(src, 1);
    const __m256i sum = _mm256_add_epi32(src_L, src_H);

    return Hadd32_AVX2_INTRIN(sum);
}

static INLINE void Distortion_AVX512_INTRIN(const __m256i input,
    const __m256i recon, __m512i *const sum) {
    const __m512i in = _mm512_cvtepu8_epi16(input);
    const __m512i re = _mm512_cvtepu8_epi16(recon);
    const __m512i diff = _mm512_sub_epi16(in, re);
    const __m512i dist = _mm512_madd_epi16(diff, diff);
    *sum = _mm512_add_epi32(*sum, dist);
}

#if 0
// Slightly slower than AVX2 version for small area_height. Disabled.
uint64_t SpatialFullDistortionKernel16xN_AVX512_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    int32_t row_count = area_height;
    __m512i sum = _mm512_setzero_si512();

    (void)area_width;

    do
    {
        const __m128i in0 = _mm_loadu_si128((__m128i *)input);
        const __m128i in1 = _mm_loadu_si128((__m128i *)(input +input_stride));
        const __m128i re0 = _mm_loadu_si128((__m128i *)recon);
        const __m128i re1 = _mm_loadu_si128((__m128i *)(recon + recon_stride));
        const __m256i in = _mm256_setr_m128i(in0, in1);
        const __m256i re = _mm256_setr_m128i(re0, re1);
        Distortion_AVX512_INTRIN(in, re, &sum);
        input += 2 * input_stride;
        recon += 2 * recon_stride;
        row_count -= 2;
    } while (row_count);

    return Hadd32_AVX512_INTRIN(sum);
}
#endif

static INLINE void SpatialFullDistortionKernel32_AVX512_INTRIN(
    const uint8_t *const input, const uint8_t *const recon, __m512i *const sum)
{
    const __m256i in = _mm256_loadu_si256((__m256i *)input);
    const __m256i re = _mm256_loadu_si256((__m256i *)recon);
    Distortion_AVX512_INTRIN(in, re, sum);
}

static INLINE void SpatialFullDistortionKernel64_AVX512_INTRIN(
    const uint8_t *const input, const uint8_t *const recon, __m512i *const sum)
{
    const __m512i in = _mm512_loadu_si512((__m512i *)input);
    const __m512i re = _mm512_loadu_si512((__m512i *)recon);
    const __m512i max = _mm512_max_epu8(in, re);
    const __m512i min = _mm512_min_epu8(in, re);
    const __m512i diff = _mm512_sub_epi8(max, min);
    const __m512i diff_L = _mm512_unpacklo_epi8(diff, _mm512_setzero_si512());
    const __m512i diff_H = _mm512_unpackhi_epi8(diff, _mm512_setzero_si512());
    const __m512i dist_L = _mm512_madd_epi16(diff_L, diff_L);
    const __m512i dist_H = _mm512_madd_epi16(diff_H, diff_H);
    const __m512i dist = _mm512_add_epi32(dist_L, dist_H);
    *sum = _mm512_add_epi32(*sum, dist);
}

uint64_t spatial_full_distortion_kernel32x_n_avx512_intrin(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    int32_t row_count = area_height;
    __m512i sum = _mm512_setzero_si512();

    (void)area_width;

    do
    {
        const __m256i in = _mm256_loadu_si256((__m256i *)input);
        const __m256i re = _mm256_loadu_si256((__m256i *)recon);
        Distortion_AVX512_INTRIN(in, re, &sum);
        input += input_stride;
        recon += recon_stride;
    } while (--row_count);

    return Hadd32_AVX512_INTRIN(sum);
}

uint64_t spatial_full_distortion_kernel64x_n_avx512_intrin(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    int32_t row_count = area_height;
    __m512i sum = _mm512_setzero_si512();

    (void)area_width;

    do
    {
        SpatialFullDistortionKernel64_AVX512_INTRIN(input, recon, &sum);
        input += input_stride;
        recon += recon_stride;
    } while (--row_count);

    return Hadd32_AVX512_INTRIN(sum);
}

uint64_t spatial_full_distortion_kernel128x_n_avx512_intrin(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    int32_t row_count = area_height;
    __m512i sum = _mm512_setzero_si512();

    (void)area_width;

    do
    {
        SpatialFullDistortionKernel64_AVX512_INTRIN(input, recon, &sum);
        SpatialFullDistortionKernel64_AVX512_INTRIN(input + 64, recon + 64, &sum);
        input += input_stride;
        recon += recon_stride;
    } while (--row_count);

    return Hadd32_AVX512_INTRIN(sum);
}

uint64_t spatial_full_distortion_kernel_avx512(
    uint8_t   *input,
    uint32_t   input_offset,
    uint32_t   input_stride,
    uint8_t   *recon,
    uint32_t   recon_offset,
    uint32_t   recon_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const uint32_t leftover = area_width & 31;
    int32_t h;
    __m256i sum = _mm256_setzero_si256();
    __m128i sum_L, sum_H, s;
    uint64_t spatialDistortion = 0;
    input += input_offset;
    recon += recon_offset;

    if (leftover) {
        const uint8_t *inp = input + area_width - leftover;
        const uint8_t *rec = recon + area_width - leftover;

        if (leftover == 4) {
            h = area_height;
            do {
                const __m128i in0 = _mm_cvtsi32_si128(*(uint32_t *)inp);
                const __m128i in1 = _mm_cvtsi32_si128(*(uint32_t *)(inp + input_stride));
                const __m128i re0 = _mm_cvtsi32_si128(*(uint32_t *)rec);
                const __m128i re1 = _mm_cvtsi32_si128(*(uint32_t *)(rec + recon_stride));
                const __m256i in = _mm256_setr_m128i(in0, in1);
                const __m256i re = _mm256_setr_m128i(re0, re1);
                Distortion_AVX2_INTRIN(in, re, &sum);
                inp += 2 * input_stride;
                rec += 2 * recon_stride;
                h -= 2;
            } while (h);

            if (area_width == 4) {
                sum_L = _mm256_extracti128_si256(sum, 0);
                sum_H = _mm256_extracti128_si256(sum, 1);
                s = _mm_add_epi32(sum_L, sum_H);
                s = _mm_add_epi32(s, _mm_srli_si128(s, 4));
                spatialDistortion = _mm_cvtsi128_si32(s);
                return spatialDistortion;
            }
        }
        else if (leftover == 8) {
            h = area_height;
            do {
                const __m128i in0 = _mm_loadl_epi64((__m128i *)inp);
                const __m128i in1 = _mm_loadl_epi64((__m128i *)(inp + input_stride));
                const __m128i re0 = _mm_loadl_epi64((__m128i *)rec);
                const __m128i re1 = _mm_loadl_epi64((__m128i *)(rec + recon_stride));
                const __m256i in = _mm256_setr_m128i(in0, in1);
                const __m256i re = _mm256_setr_m128i(re0, re1);
                Distortion_AVX2_INTRIN(in, re, &sum);
                inp += 2 * input_stride;
                rec += 2 * recon_stride;
                h -= 2;
            } while (h);
        }
        else if (leftover <= 16) {
            h = area_height;
            do {
                SpatialFullDistortionKernel16_AVX2_INTRIN(inp, rec, &sum);
                inp += input_stride;
                rec += recon_stride;
            } while (--h);

            if (leftover == 12) {
                const __m256i mask = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, 0, 0);
                sum = _mm256_and_si256(sum, mask);
            }
        }
        else {
            __m256i sum1 = _mm256_setzero_si256();
            h = area_height;
            do {
                SpatialFullDistortionKernel32Leftover_AVX2_INTRIN(inp, rec, &sum, &sum1);
                inp += input_stride;
                rec += recon_stride;
            } while (--h);

            __m256i mask[2];
            if (leftover == 20) {
                mask[0] = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, 0, 0);
                mask[1] = _mm256_setr_epi32(-1, -1, -1, -1, 0, 0, 0, 0);
            }
            else if (leftover == 24) {
                mask[0] = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, -1, -1);
                mask[1] = _mm256_setr_epi32(-1, -1, -1, -1, 0, 0, 0, 0);
            }
            else { // leftover = 28
                mask[0] = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, -1, -1);
                mask[1] = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, 0, 0);
            }

            sum = _mm256_and_si256(sum, mask[0]);
            sum1 = _mm256_and_si256(sum1, mask[1]);
            sum = _mm256_add_epi32(sum, sum1);
        }
    }

    area_width -= leftover;

    if (area_width) {
        const uint8_t *inp = input;
        const uint8_t *rec = recon;
        h = area_height;

        if (area_width == 32) {
            do {
                SpatialFullDistortionKernel32_AVX2_INTRIN(inp, rec, &sum);
                inp += input_stride;
                rec += recon_stride;
            } while (--h);
        }
        else {
            __m512i sum512 = _mm512_setzero_si512();

            if (area_width == 64) {
                do {
                    SpatialFullDistortionKernel64_AVX512_INTRIN(inp, rec, &sum512);
                    inp += input_stride;
                    rec += recon_stride;
                } while (--h);
            }
            else if (area_width == 96) {
                do {
                    SpatialFullDistortionKernel64_AVX512_INTRIN(inp, rec, &sum512);
                    SpatialFullDistortionKernel32_AVX512_INTRIN(inp + 64, rec + 64, &sum512);
                    inp += input_stride;
                    rec += recon_stride;
                } while (--h);
            }
            else { // 128
                do {
                    SpatialFullDistortionKernel64_AVX512_INTRIN(inp, rec, &sum512);
                    SpatialFullDistortionKernel64_AVX512_INTRIN(inp + 64, rec + 64, &sum512);
                    inp += input_stride;
                    rec += recon_stride;
                } while (--h);
            }

            const __m256i sum512_L = _mm512_extracti64x4_epi64(sum512, 0);
            const __m256i sum512_H = _mm512_extracti64x4_epi64(sum512, 1);
            sum = _mm256_add_epi32(sum, sum512_L);
            sum = _mm256_add_epi32(sum, sum512_H);
        }
    }

    return Hadd32_AVX2_INTRIN(sum);
}

#endif // !NON_AVX512_SUPPORT
