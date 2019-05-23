/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/


#include "EbCombinedAveragingSAD_Intrinsic_AVX2.h"
#include "immintrin.h"
#include "EbMemory_AVX2.h"


#define _mm256_set_m128i(/* __m128i */ hi, /* __m128i */ lo) \
    _mm256_insertf128_si256(_mm256_castsi128_si256(lo), (hi), 0x1)


uint32_t combined_averaging8x_msad_avx2_intrin(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1_stride,
    uint8_t  *ref2,
    uint32_t  ref2_stride,
    uint32_t  height,
    uint32_t  width)
{
    __m256i sum = _mm256_setzero_si256();
    __m128i sad;
    uint32_t y;
    (void)width;

    for (y = 0; y < height; y += 4) {
        const __m256i s = load8bit_8x4_avx2(src, src_stride);
        const __m256i r1 = load8bit_8x4_avx2(ref1, ref1_stride);
        const __m256i r2 = load8bit_8x4_avx2(ref2, ref2_stride);
        const __m256i avg = _mm256_avg_epu8(r1, r2);
        const __m256i sad = _mm256_sad_epu8(s, avg);
        sum = _mm256_add_epi32(sum, sad);
        src += src_stride << 2;
        ref1 += ref1_stride << 2;
        ref2 += ref2_stride << 2;
    }

    sad = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm256_extracti128_si256(sum, 1));
    sad = _mm_add_epi32(sad, _mm_srli_si128(sad, 8));

    return _mm_cvtsi128_si32(sad);
}

static INLINE __m256i CombinedAveragingSad16x2_AVX2(const uint8_t *const src,
    const uint32_t src_stride, const uint8_t *const ref1, const uint32_t ref1_stride,
    const uint8_t *const ref2, const uint32_t ref2_stride, const __m256i sum)
{
    const __m256i s = load8bit_16x2_unaligned_avx2(src, src_stride);
    const __m256i r1 = load8bit_16x2_unaligned_avx2(ref1, ref1_stride);
    const __m256i r2 = load8bit_16x2_unaligned_avx2(ref2, ref2_stride);
    const __m256i avg = _mm256_avg_epu8(r1, r2);
    const __m256i sad = _mm256_sad_epu8(s, avg);
    return _mm256_add_epi32(sum, sad);
}

uint32_t combined_averaging16x_msad_avx2_intrin(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1_stride,
    uint8_t  *ref2,
    uint32_t  ref2_stride,
    uint32_t  height,
    uint32_t  width)
{
    __m256i sum = _mm256_setzero_si256();
    __m128i sad;
    uint32_t y;
    (void)width;

    for (y = 0; y < height; y += 2) {
        sum = CombinedAveragingSad16x2_AVX2(src, src_stride, ref1, ref1_stride,
            ref2, ref2_stride, sum);
        src += src_stride << 1;
        ref1 += ref1_stride << 1;
        ref2 += ref2_stride << 1;
    }

    sad = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm256_extracti128_si256(sum, 1));
    sad = _mm_add_epi32(sad, _mm_srli_si128(sad, 8));

    return _mm_cvtsi128_si32(sad);
}

static INLINE __m256i CombinedAveragingSad24_AVX2(const uint8_t *const src,
    const uint8_t *const ref1, const uint8_t *const ref2, const __m256i sum)
{
    const __m256i s = _mm256_loadu_si256((__m256i*)src);
    const __m256i r1 = _mm256_loadu_si256((__m256i*)ref1);
    const __m256i r2 = _mm256_loadu_si256((__m256i*)ref2);
    const __m256i avg = _mm256_avg_epu8(r1, r2);
    const __m256i sad = _mm256_sad_epu8(s, avg);
    return _mm256_add_epi32(sum, sad);
}

uint32_t combined_averaging24x_msad_avx2_intrin(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1_stride,
    uint8_t  *ref2,
    uint32_t  ref2_stride,
    uint32_t  height,
    uint32_t  width)
{
    __m256i sum = _mm256_setzero_si256();
    __m128i sad;
    uint32_t y;
    (void)width;

    for (y = 0; y < height; y += 2) {
        sum = CombinedAveragingSad24_AVX2(src + 0 * src_stride,
            ref1 + 0 * ref1_stride, ref2 + 0 * ref2_stride, sum);
        sum = CombinedAveragingSad24_AVX2(src + 1 * src_stride,
            ref1 + 1 * ref1_stride, ref2 + 1 * ref2_stride, sum);
        src += src_stride << 1;
        ref1 += ref1_stride << 1;
        ref2 += ref2_stride << 1;
    }

    sad = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm_slli_si128(_mm256_extracti128_si256(sum, 1), 8));
    sad = _mm_add_epi32(sad, _mm_srli_si128(sad, 8));

    return _mm_cvtsi128_si32(sad);
}

static INLINE __m256i CombinedAveragingSad32_AVX2(const uint8_t *const src,
    const uint8_t *const ref1, const uint8_t *const ref2, const __m256i sum)
{
    const __m256i s = _mm256_loadu_si256((__m256i*)src);
    const __m256i r1 = _mm256_loadu_si256((__m256i*)ref1);
    const __m256i r2 = _mm256_loadu_si256((__m256i*)ref2);
    const __m256i avg = _mm256_avg_epu8(r1, r2);
    const __m256i sad = _mm256_sad_epu8(s, avg);
    return _mm256_add_epi32(sum, sad);
}

uint32_t combined_averaging32x_msad_avx2_intrin(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1_stride,
    uint8_t  *ref2,
    uint32_t  ref2_stride,
    uint32_t  height,
    uint32_t  width)
{
    __m256i sum = _mm256_setzero_si256();
    __m128i sad;
    uint32_t y;
    (void)width;

    for (y = 0; y < height; y += 2) {
        sum = CombinedAveragingSad32_AVX2(src + 0 * src_stride,
            ref1 + 0 * ref1_stride, ref2 + 0 * ref2_stride, sum);
        sum = CombinedAveragingSad32_AVX2(src + 1 * src_stride,
            ref1 + 1 * ref1_stride, ref2 + 1 * ref2_stride, sum);
        src += src_stride << 1;
        ref1 += ref1_stride << 1;
        ref2 += ref2_stride << 1;
    }

    sad = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm256_extracti128_si256(sum, 1));
    sad = _mm_add_epi32(sad, _mm_srli_si128(sad, 8));

    return _mm_cvtsi128_si32(sad);
}

uint32_t combined_averaging48x_msad_avx2_intrin(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1_stride,
    uint8_t  *ref2,
    uint32_t  ref2_stride,
    uint32_t  height,
    uint32_t  width)
{
    __m256i sum = _mm256_setzero_si256();
    __m128i sad;
    uint32_t y;
    (void)width;

    for (y = 0; y < height; y += 2) {
        sum = CombinedAveragingSad32_AVX2(src + 0 * src_stride,
            ref1 + 0 * ref1_stride, ref2 + 0 * ref2_stride, sum);
        sum = CombinedAveragingSad32_AVX2(src + 1 * src_stride,
            ref1 + 1 * ref1_stride, ref2 + 1 * ref2_stride, sum);
        sum = CombinedAveragingSad16x2_AVX2(src + 32, src_stride, ref1 + 32,
            ref1_stride, ref2 + 32, ref2_stride, sum);

        src += src_stride << 1;
        ref1 += ref1_stride << 1;
        ref2 += ref2_stride << 1;
    }

    sad = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm256_extracti128_si256(sum, 1));
    sad = _mm_add_epi32(sad, _mm_srli_si128(sad, 8));

    return _mm_cvtsi128_si32(sad);
}

uint32_t combined_averaging64x_msad_avx2_intrin(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1_stride,
    uint8_t  *ref2,
    uint32_t  ref2_stride,
    uint32_t  height,
    uint32_t  width)
{
    __m256i sum = _mm256_setzero_si256();
    __m128i sad;
    uint32_t y;
    (void)width;

    for (y = 0; y < height; y++) {
        sum = CombinedAveragingSad32_AVX2(src + 0x00,
            ref1 + 0x00, ref2 + 0x00, sum);
        sum = CombinedAveragingSad32_AVX2(src + 0x20,
            ref1 + 0x20, ref2 + 0x20, sum);
        src += src_stride;
        ref1 += ref1_stride;
        ref2 += ref2_stride;
    }

    sad = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm256_extracti128_si256(sum, 1));
    sad = _mm_add_epi32(sad, _mm_srli_si128(sad, 8));

    return _mm_cvtsi128_si32(sad);
}
uint64_t compute_mean8x8_avx2_intrin(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint32_t   input_stride,       // input parameter, input stride
    uint32_t   input_area_width,    // input parameter, input area width
    uint32_t   input_area_height)   // input parameter, input area height
{
    __m256i sum, sum2, xmm2, xmm1, sum1, xmm0 = _mm256_setzero_si256();
    __m128i  upper, lower, mean = _mm_setzero_si128();
    uint64_t result;
    xmm1 = _mm256_sad_epu8(xmm0, _mm256_set_m128i(_mm_loadl_epi64((__m128i *)(input_samples + input_stride)), _mm_loadl_epi64((__m128i *)(input_samples))));
    xmm2 = _mm256_sad_epu8(xmm0, _mm256_set_m128i(_mm_loadl_epi64((__m128i *)(input_samples + 3 * input_stride)), _mm_loadl_epi64((__m128i *)(input_samples + 2 * input_stride))));
    sum1 = _mm256_add_epi16(xmm1, xmm2);

    input_samples += 4 * input_stride;

    xmm1 = _mm256_sad_epu8(xmm0, _mm256_set_m128i(_mm_loadl_epi64((__m128i *)(input_samples + input_stride)), _mm_loadl_epi64((__m128i *)(input_samples))));
    xmm2 = _mm256_sad_epu8(xmm0, _mm256_set_m128i(_mm_loadl_epi64((__m128i *)(input_samples + 3 * input_stride)), _mm_loadl_epi64((__m128i *)(input_samples + 2 * input_stride))));
    sum2 = _mm256_add_epi16(xmm1, xmm2);

    sum = _mm256_add_epi16(sum1, sum2);
    upper = _mm256_extractf128_si256(sum, 1); //extract upper 128 bit
    upper = _mm_add_epi32(upper, _mm_srli_si128(upper, 8)); // shift 2nd 16 bits to the 1st and sum both

    lower = _mm256_extractf128_si256(sum, 0); //extract lower 128 bit
    lower = _mm_add_epi32(lower, _mm_srli_si128(lower, 8)); // shift 2nd 16 bits to the 1st and sum both

    mean = _mm_add_epi32(lower, upper);

    (void)input_area_width;
    (void)input_area_height;

    result = (uint64_t)_mm_cvtsi128_si32(mean) << 2;
    return result;

}

/********************************************************************************************************************************/
void  compute_interm_var_four8x8_avx2_intrin(
    uint8_t *  input_samples,
    uint16_t   input_stride,
    uint64_t * mean_of8x8_blocks,      // mean of four  8x8
    uint64_t * mean_of_squared8x8_blocks)  // meanSquared
{

    __m256i ymm1, ymm2, ymm3, ymm4, ymm_sum1, ymm_sum2, ymm_FinalSum, ymm_shift,/* ymm_blockMeanSquared*///,
        ymm_in, ymm_in_2S, ymm_in_second, ymm_in_2S_second, ymm_shiftSquared, ymm_permute8,
        ymm_result, ymm_blockMeanSquaredlow, ymm_blockMeanSquaredHi, ymm_inputlo, ymm_inputhi;

    __m128i ymm_blockMeanSquaredlo, ymm_blockMeanSquaredhi, ymm_resultlo, ymm_resulthi;

    __m256i ymm_zero = _mm256_setzero_si256();
    __m128i xmm_zero = _mm_setzero_si128();

    ymm_in = _mm256_loadu_si256((__m256i *) input_samples);
    ymm_in_2S = _mm256_loadu_si256((__m256i *)(input_samples + 2 * input_stride));

    ymm1 = _mm256_sad_epu8(ymm_in, ymm_zero);
    ymm2 = _mm256_sad_epu8(ymm_in_2S, ymm_zero);

    ymm_sum1 = _mm256_add_epi16(ymm1, ymm2);

    input_samples += 4 * input_stride;
    ymm_in_second = _mm256_loadu_si256((__m256i *)input_samples);
    ymm_in_2S_second = _mm256_loadu_si256((__m256i *)(input_samples + 2 * input_stride));

    ymm3 = _mm256_sad_epu8(ymm_in_second, ymm_zero);
    ymm4 = _mm256_sad_epu8(ymm_in_2S_second, ymm_zero);

    ymm_sum2 = _mm256_add_epi16(ymm3, ymm4);

    ymm_FinalSum = _mm256_add_epi16(ymm_sum1, ymm_sum2);

    ymm_shift = _mm256_set_epi64x(3, 3, 3, 3);
    ymm_FinalSum = _mm256_sllv_epi64(ymm_FinalSum, ymm_shift);

    _mm256_storeu_si256((__m256i *)(mean_of8x8_blocks), ymm_FinalSum);

    /*******************************Squared Mean******************************/

    ymm_inputlo = _mm256_unpacklo_epi8(ymm_in, ymm_zero);
    ymm_inputhi = _mm256_unpackhi_epi8(ymm_in, ymm_zero);

    ymm_blockMeanSquaredlow = _mm256_madd_epi16(ymm_inputlo, ymm_inputlo);
    ymm_blockMeanSquaredHi = _mm256_madd_epi16(ymm_inputhi, ymm_inputhi);

    ymm_inputlo = _mm256_unpacklo_epi8(ymm_in_2S, ymm_zero);
    ymm_inputhi = _mm256_unpackhi_epi8(ymm_in_2S, ymm_zero);

    ymm_blockMeanSquaredlow = _mm256_add_epi32(ymm_blockMeanSquaredlow, _mm256_madd_epi16(ymm_inputlo, ymm_inputlo));
    ymm_blockMeanSquaredHi = _mm256_add_epi32(ymm_blockMeanSquaredHi, _mm256_madd_epi16(ymm_inputhi, ymm_inputhi));

    ymm_inputlo = _mm256_unpacklo_epi8(ymm_in_second, ymm_zero);
    ymm_inputhi = _mm256_unpackhi_epi8(ymm_in_second, ymm_zero);

    ymm_blockMeanSquaredlow = _mm256_add_epi32(ymm_blockMeanSquaredlow, _mm256_madd_epi16(ymm_inputlo, ymm_inputlo));
    ymm_blockMeanSquaredHi = _mm256_add_epi32(ymm_blockMeanSquaredHi, _mm256_madd_epi16(ymm_inputhi, ymm_inputhi));

    ymm_inputlo = _mm256_unpacklo_epi8(ymm_in_2S_second, ymm_zero);
    ymm_inputhi = _mm256_unpackhi_epi8(ymm_in_2S_second, ymm_zero);

    ymm_blockMeanSquaredlow = _mm256_add_epi32(ymm_blockMeanSquaredlow, _mm256_madd_epi16(ymm_inputlo, ymm_inputlo));
    ymm_blockMeanSquaredHi = _mm256_add_epi32(ymm_blockMeanSquaredHi, _mm256_madd_epi16(ymm_inputhi, ymm_inputhi));

    ymm_blockMeanSquaredlow = _mm256_add_epi32(ymm_blockMeanSquaredlow, _mm256_srli_si256(ymm_blockMeanSquaredlow, 8));
    ymm_blockMeanSquaredHi = _mm256_add_epi32(ymm_blockMeanSquaredHi, _mm256_srli_si256(ymm_blockMeanSquaredHi, 8));

    ymm_blockMeanSquaredlow = _mm256_add_epi32(ymm_blockMeanSquaredlow, _mm256_srli_si256(ymm_blockMeanSquaredlow, 4));
    ymm_blockMeanSquaredHi = _mm256_add_epi32(ymm_blockMeanSquaredHi, _mm256_srli_si256(ymm_blockMeanSquaredHi, 4));

    ymm_permute8 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 4, 0);
    ymm_blockMeanSquaredlow = _mm256_permutevar8x32_epi32(ymm_blockMeanSquaredlow, ymm_permute8/*8*/);
    ymm_blockMeanSquaredHi = _mm256_permutevar8x32_epi32(ymm_blockMeanSquaredHi, ymm_permute8);

    ymm_blockMeanSquaredlo = _mm256_extracti128_si256(ymm_blockMeanSquaredlow, 0); //lower 128
    ymm_blockMeanSquaredhi = _mm256_extracti128_si256(ymm_blockMeanSquaredHi, 0); //lower 128

    ymm_result = _mm256_unpacklo_epi32(_mm256_castsi128_si256(ymm_blockMeanSquaredlo), _mm256_castsi128_si256(ymm_blockMeanSquaredhi));
    ymm_resultlo = _mm_unpacklo_epi64(_mm256_castsi256_si128(ymm_result), xmm_zero);
    ymm_resulthi = _mm_unpackhi_epi64(_mm256_castsi256_si128(ymm_result), xmm_zero);


    ymm_result = _mm256_set_m128i(ymm_resulthi, ymm_resultlo);

    ymm_permute8 = _mm256_set_epi32(7, 5, 6, 4, 3, 1, 2, 0);
    ymm_result = _mm256_permutevar8x32_epi32(ymm_result, ymm_permute8);

    ymm_shiftSquared = _mm256_set1_epi64x(11);

    ymm_result = _mm256_sllv_epi64(ymm_result, ymm_shiftSquared);


    _mm256_storeu_si256((__m256i *)(mean_of_squared8x8_blocks), ymm_result);


}

uint32_t combined_averaging_ssd_avx2(
    uint8_t   *src,
    ptrdiff_t  src_stride,
    uint8_t   *ref1,
    ptrdiff_t  ref1_stride,
    uint8_t   *ref2,
    ptrdiff_t  ref2_stride,
    uint32_t   height,
    uint32_t   width)
{
    __m256i sum = _mm256_setzero_si256();
    __m256i s, r1, r2, avg, ssd;
    __m128i sum1_128, sum2_128;
    uint32_t row, col;

    for (row = 0; row < height; row++)
    {
        for (col = 0; col < width; col += 4)
        {
            s = _mm256_setzero_si256();
            r1 = _mm256_setzero_si256();
            r2 = _mm256_setzero_si256();
            s = _mm256_insert_epi32(s, *(int32_t *)(src + col), 0);
            r1 = _mm256_insert_epi32(r1, *(int32_t *)(ref1 + col), 0);
            r2 = _mm256_insert_epi32(r2, *(int32_t *)(ref2 + col), 0);

            avg = _mm256_avg_epu8(r1, r2);
            ssd = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(avg));
            s = _mm256_cvtepu8_epi64(_mm256_castsi256_si128(s));
            ssd = _mm256_sub_epi64(s, ssd);
            ssd = _mm256_mul_epi32(ssd, ssd);

            sum = _mm256_add_epi64(sum, ssd);
        }
        src += src_stride;
        ref1 += ref1_stride;
        ref2 += ref2_stride;
    }
    sum1_128 = _mm256_extracti128_si256(sum, 0);
    sum2_128 = _mm256_extracti128_si256(sum, 1);
    sum1_128 = _mm_add_epi64(sum1_128, sum2_128);
    sum2_128 = _mm_shuffle_epi32(sum1_128, 0xc6);
    sum1_128 = _mm_add_epi64(sum1_128, sum2_128);
    return _mm_cvtsi128_si32(sum1_128);
}



