/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "emmintrin.h"
#include "EbComputeMean_SSE2.h"

uint64_t ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint16_t   inputStride)       // input parameter, input stride

{
    __m128i xmm0, xmm_blockMean, xmm_input;

    xmm0 = _mm_setzero_si128();
    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)input_samples), xmm0);
    xmm_blockMean = _mm_madd_epi16(xmm_input, xmm_input);

    /*xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));*/

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 2 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    /*xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples+3*inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));*/

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 4 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    //xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples+5*inputStride)), xmm0);
    //xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 6 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    /*xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples+7*inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));*/

    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_srli_si128(xmm_blockMean, 8));
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_srli_si128(xmm_blockMean, 4));


    return (uint64_t)_mm_cvtsi128_si32(xmm_blockMean) << 11;




}

uint64_t ComputeSubMean8x8_SSE2_INTRIN(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint16_t   inputStride)       // input parameter, input stride

{

    __m128i xmm0 = _mm_setzero_si128(), xmm1, xmm3, xmm_sum1, xmm_sum2;

    xmm1 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples)), xmm0);
    //xmm2 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + inputStride)), xmm0);
    xmm3 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + 2 * inputStride)), xmm0);
    //xmm4 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + 3 * inputStride)), xmm0);
    xmm_sum1 = _mm_add_epi16(xmm1, xmm3);

    input_samples += 4 * inputStride;
    xmm1 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples)), xmm0);
    //xmm2 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + inputStride)), xmm0);
    xmm3 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + 2 * inputStride)), xmm0);
    //xmm4 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + 3 * inputStride)), xmm0);
    xmm_sum2 = _mm_add_epi16(xmm1, xmm3);
    xmm_sum2 = _mm_add_epi16(xmm_sum1, xmm_sum2);

    return (uint64_t)_mm_cvtsi128_si32(xmm_sum2) << 3;

}



uint64_t ComputeMeanOfSquaredValues8x8_SSE2_INTRIN(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint32_t   inputStride,       // input parameter, input stride
    uint32_t   inputAreaWidth,    // input parameter, input area width
    uint32_t   inputAreaHeight)   // input parameter, input area height
{
    __m128i xmm0, xmm_blockMean, xmm_input;
    (void)inputAreaWidth;
    (void)inputAreaHeight;
    xmm0 = _mm_setzero_si128();
    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)input_samples), xmm0);
    xmm_blockMean = _mm_madd_epi16(xmm_input, xmm_input);

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 2 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 3 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 4 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 5 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 6 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    xmm_input = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(input_samples + 7 * inputStride)), xmm0);
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_madd_epi16(xmm_input, xmm_input));

    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_srli_si128(xmm_blockMean, 8));
    xmm_blockMean = _mm_add_epi32(xmm_blockMean, _mm_srli_si128(xmm_blockMean, 4));

    return (uint64_t)_mm_cvtsi128_si32(xmm_blockMean) << 10;
}

uint64_t ComputeMean8x8_SSE2_INTRIN(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint32_t   inputStride,       // input parameter, input stride
    uint32_t   inputAreaWidth,    // input parameter, input area width
    uint32_t   inputAreaHeight)   // input parameter, input area height
{
    __m128i xmm0 = _mm_setzero_si128(), xmm1, xmm2, xmm3, xmm4, xmm_sum1, xmm_sum2;

    xmm1 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples)), xmm0);
    xmm2 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + inputStride)), xmm0);
    xmm3 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + 2 * inputStride)), xmm0);
    xmm4 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + 3 * inputStride)), xmm0);
    xmm_sum1 = _mm_add_epi16(_mm_add_epi16(xmm1, xmm2), _mm_add_epi16(xmm3, xmm4));

    input_samples += 4 * inputStride;
    xmm1 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples)), xmm0);
    xmm2 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + inputStride)), xmm0);
    xmm3 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + 2 * inputStride)), xmm0);
    xmm4 = _mm_sad_epu8(_mm_loadl_epi64((__m128i *)(input_samples + 3 * inputStride)), xmm0);
    xmm_sum2 = _mm_add_epi16(_mm_add_epi16(xmm1, xmm2), _mm_add_epi16(xmm3, xmm4));
    xmm_sum2 = _mm_add_epi16(xmm_sum1, xmm_sum2);

    (void)inputAreaWidth;
    (void)inputAreaHeight;

    return (uint64_t)_mm_cvtsi128_si32(xmm_sum2) << 2;
}
