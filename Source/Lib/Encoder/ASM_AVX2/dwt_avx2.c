/*
* Copyright(c) 2020 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include "transpose_avx2.h"

int svt_av1_haar_ac_sad_8x8_uint8_input_avx2(uint8_t *input, int stride, int hbd) {
    DECLARE_ALIGNED(32, int32_t, output[64]);

    int32_t *out_ptr = output;
    int      i;

    if (hbd) {
        uint16_t *x16 = CONVERT_TO_SHORTPTR(input);
        for (i = 0; i < 8; i++) {
            __m256i inp = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)x16));
            _mm256_storeu_si256((__m256i *)out_ptr, _mm256_slli_epi32(inp, 2));

            x16 += stride;
            out_ptr += 8;
        }
    } else {
        for (i = 0; i < 8; i++) {
            __m256i inp = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i *)input));
            _mm256_storeu_si256((__m256i *)out_ptr, _mm256_slli_epi32(inp, 2));

            input += stride;
            out_ptr += 8;
        }
    }

    const __m256i indices = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    const __m256i one     = _mm256_set1_epi32(1);
    const __m256i two     = _mm256_set1_epi32(2);

    for (i = 0; i < 4; i++) {
        __m256i x_2n   = _mm256_i32gather_epi32(output + i * 16, indices, 8);
        __m256i x_2np1 = _mm256_i32gather_epi32(output + i * 16 + 1, indices, 8);
        __m256i x_2np2 = _mm256_i32gather_epi32(output + i * 16 + 2, indices, 8);

        __m256i a = _mm256_slli_epi32(x_2n, 1);
        __m256i b = _mm256_add_epi32(_mm256_add_epi32(x_2n, x_2np2), one);
        b         = _mm256_sub_epi32(x_2np1, _mm256_srli_epi32(b, 1));

        int32_t idx_3 = output[i * 16 + 7] - output[i * 16 + 6];
        int32_t idx_7 = output[i * 16 + 15] - output[i * 16 + 14];
        b             = _mm256_insert_epi32(b, idx_3, 3);
        b             = _mm256_insert_epi32(b, idx_7, 7);

        __m256i r = _mm256_shuffle_epi32(b, (1 << 4) + (2 << 6));

        a = _mm256_add_epi32(a,
                             _mm256_srai_epi32(_mm256_add_epi32(_mm256_add_epi32(r, b), one), 1));

        _mm_storeu_si128((__m128i *)(output + i * 16), _mm256_castsi256_si128(a));
        _mm_storeu_si128((__m128i *)(output + i * 16 + 4), _mm256_castsi256_si128(b));
        _mm_storeu_si128((__m128i *)(output + i * 16 + 8), _mm256_extracti128_si256(a, 0x1));
        _mm_storeu_si128((__m128i *)(output + i * 16 + 12), _mm256_extracti128_si256(b, 0x1));
    }

    transpose_32bit_8x8_avx2((const __m256i *const)output, (__m256i *const)output);

    __m256i sum[8];
    for (i = 0; i < 4; i++) {
        __m256i x_2n   = _mm256_i32gather_epi32(output + i * 16, indices, 8);
        __m256i x_2np1 = _mm256_i32gather_epi32(output + i * 16 + 1, indices, 8);
        __m256i x_2np2 = _mm256_i32gather_epi32(output + i * 16 + 2, indices, 8);

        __m256i a = x_2n;
        __m256i b = _mm256_sub_epi32(_mm256_slli_epi32(x_2np1, 1), x_2n);
        b         = _mm256_add_epi32(_mm256_sub_epi32(b, x_2np2), two);
        b         = _mm256_srai_epi32(b, 2);

        int32_t idx_3 = (output[i * 16 + 7] - output[i * 16 + 6] + 1) >> 1;
        int32_t idx_7 = (output[i * 16 + 15] - output[i * 16 + 14] + 1) >> 1;
        b             = _mm256_insert_epi32(b, idx_3, 3);
        b             = _mm256_insert_epi32(b, idx_7, 7);

        __m256i r = _mm256_shuffle_epi32(b, (1 << 4) + (2 << 6));

        a = _mm256_add_epi32(a,
                             _mm256_srai_epi32(_mm256_add_epi32(_mm256_add_epi32(r, b), one), 1));

        sum[2 * i]     = _mm256_abs_epi32(a);
        sum[2 * i + 1] = _mm256_abs_epi32(b);
    }

    __m256i sum13, sum45, sum67;

    sum13 = _mm256_add_epi32(sum[1], sum[3]);
    sum45 = _mm256_add_epi32(sum[4], sum[5]);
    sum67 = _mm256_add_epi32(sum[6], sum[7]);
    sum13 = _mm256_add_epi32(sum13, _mm256_add_epi32(sum45, sum67));

    __m128i sum_128 = _mm_add_epi32(_mm256_castsi256_si128(sum13),
                                    _mm256_extracti128_si256(sum13, 1));

    sum_128 = _mm_hadd_epi32(sum_128, sum_128);
    sum_128 = _mm_hadd_epi32(sum_128, sum_128);

    return _mm_cvtsi128_si32(sum_128);
}
