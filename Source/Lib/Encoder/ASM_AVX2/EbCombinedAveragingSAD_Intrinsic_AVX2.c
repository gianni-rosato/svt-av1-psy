/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <immintrin.h>
#include "EbCombinedAveragingSAD_Inline_AVX2.h"
#include "EbMemory_AVX2.h"
#include "EbMemory_SSE4_1.h"
/********************************************************************************************************************************/
void svt_compute_interm_var_four8x8_avx2_intrin(uint8_t *input_samples, uint16_t input_stride,
                                                uint64_t *mean_of8x8_blocks, // mean of four  8x8
                                                uint64_t *mean_of_squared8x8_blocks) // meanSquared
{
    __m256i ymm1, ymm2, ymm3, ymm4, ymm_sum1, ymm_sum2, ymm_final_sum, ymm_shift,
        /* ymm_blockMeanSquared*/ //,
        ymm_in, ymm_in_2s, ymm_in_second, ymm_in_2s_second, ymm_shift_squared, ymm_permute8,
        ymm_result, ymm_block_mean_squared_low, ymm_block_mean_squared_high, ymm_inputlo,
        ymm_inputhi;

    __m128i ymm_block_mean_squared_lo, ymm_block_mean_squared_hi, ymm_resultlo, ymm_resulthi;

    __m256i ymm_zero = _mm256_setzero_si256();
    __m128i xmm_zero = _mm_setzero_si128();

    ymm_in    = _mm256_loadu_si256((__m256i *)input_samples);
    ymm_in_2s = _mm256_loadu_si256((__m256i *)(input_samples + 2 * input_stride));

    ymm1 = _mm256_sad_epu8(ymm_in, ymm_zero);
    ymm2 = _mm256_sad_epu8(ymm_in_2s, ymm_zero);

    ymm_sum1 = _mm256_add_epi16(ymm1, ymm2);

    input_samples += 4 * input_stride;
    ymm_in_second    = _mm256_loadu_si256((__m256i *)input_samples);
    ymm_in_2s_second = _mm256_loadu_si256((__m256i *)(input_samples + 2 * input_stride));

    ymm3 = _mm256_sad_epu8(ymm_in_second, ymm_zero);
    ymm4 = _mm256_sad_epu8(ymm_in_2s_second, ymm_zero);

    ymm_sum2 = _mm256_add_epi16(ymm3, ymm4);

    ymm_final_sum = _mm256_add_epi16(ymm_sum1, ymm_sum2);

    ymm_shift     = _mm256_set_epi64x(3, 3, 3, 3);
    ymm_final_sum = _mm256_sllv_epi64(ymm_final_sum, ymm_shift);

    _mm256_storeu_si256((__m256i *)(mean_of8x8_blocks), ymm_final_sum);

    /*******************************Squared Mean******************************/

    ymm_inputlo = _mm256_unpacklo_epi8(ymm_in, ymm_zero);
    ymm_inputhi = _mm256_unpackhi_epi8(ymm_in, ymm_zero);

    ymm_block_mean_squared_low  = _mm256_madd_epi16(ymm_inputlo, ymm_inputlo);
    ymm_block_mean_squared_high = _mm256_madd_epi16(ymm_inputhi, ymm_inputhi);

    ymm_inputlo = _mm256_unpacklo_epi8(ymm_in_2s, ymm_zero);
    ymm_inputhi = _mm256_unpackhi_epi8(ymm_in_2s, ymm_zero);

    ymm_block_mean_squared_low =
        _mm256_add_epi32(ymm_block_mean_squared_low, _mm256_madd_epi16(ymm_inputlo, ymm_inputlo));
    ymm_block_mean_squared_high =
        _mm256_add_epi32(ymm_block_mean_squared_high, _mm256_madd_epi16(ymm_inputhi, ymm_inputhi));

    ymm_inputlo = _mm256_unpacklo_epi8(ymm_in_second, ymm_zero);
    ymm_inputhi = _mm256_unpackhi_epi8(ymm_in_second, ymm_zero);

    ymm_block_mean_squared_low =
        _mm256_add_epi32(ymm_block_mean_squared_low, _mm256_madd_epi16(ymm_inputlo, ymm_inputlo));
    ymm_block_mean_squared_high =
        _mm256_add_epi32(ymm_block_mean_squared_high, _mm256_madd_epi16(ymm_inputhi, ymm_inputhi));

    ymm_inputlo = _mm256_unpacklo_epi8(ymm_in_2s_second, ymm_zero);
    ymm_inputhi = _mm256_unpackhi_epi8(ymm_in_2s_second, ymm_zero);

    ymm_block_mean_squared_low =
        _mm256_add_epi32(ymm_block_mean_squared_low, _mm256_madd_epi16(ymm_inputlo, ymm_inputlo));
    ymm_block_mean_squared_high =
        _mm256_add_epi32(ymm_block_mean_squared_high, _mm256_madd_epi16(ymm_inputhi, ymm_inputhi));

    ymm_block_mean_squared_low  = _mm256_add_epi32(ymm_block_mean_squared_low,
                                                  _mm256_srli_si256(ymm_block_mean_squared_low, 8));
    ymm_block_mean_squared_high = _mm256_add_epi32(
        ymm_block_mean_squared_high, _mm256_srli_si256(ymm_block_mean_squared_high, 8));

    ymm_block_mean_squared_low  = _mm256_add_epi32(ymm_block_mean_squared_low,
                                                  _mm256_srli_si256(ymm_block_mean_squared_low, 4));
    ymm_block_mean_squared_high = _mm256_add_epi32(
        ymm_block_mean_squared_high, _mm256_srli_si256(ymm_block_mean_squared_high, 4));

    ymm_permute8 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 4, 0);
    ymm_block_mean_squared_low =
        _mm256_permutevar8x32_epi32(ymm_block_mean_squared_low, ymm_permute8 /*8*/);
    ymm_block_mean_squared_high =
        _mm256_permutevar8x32_epi32(ymm_block_mean_squared_high, ymm_permute8);

    ymm_block_mean_squared_lo = _mm256_castsi256_si128(ymm_block_mean_squared_low); //lower 128
    ymm_block_mean_squared_hi =
        _mm256_extracti128_si256(ymm_block_mean_squared_high, 0); //lower 128

    ymm_result   = _mm256_unpacklo_epi32(_mm256_castsi128_si256(ymm_block_mean_squared_lo),
                                       _mm256_castsi128_si256(ymm_block_mean_squared_hi));
    ymm_resultlo = _mm_unpacklo_epi64(_mm256_castsi256_si128(ymm_result), xmm_zero);
    ymm_resulthi = _mm_unpackhi_epi64(_mm256_castsi256_si128(ymm_result), xmm_zero);

    ymm_result = _mm256_set_m128i(ymm_resulthi, ymm_resultlo);

    ymm_permute8 = _mm256_set_epi32(7, 5, 6, 4, 3, 1, 2, 0);
    ymm_result   = _mm256_permutevar8x32_epi32(ymm_result, ymm_permute8);

    ymm_shift_squared = _mm256_set1_epi64x(11);

    ymm_result = _mm256_sllv_epi64(ymm_result, ymm_shift_squared);

    _mm256_storeu_si256((__m256i *)(mean_of_squared8x8_blocks), ymm_result);
}
