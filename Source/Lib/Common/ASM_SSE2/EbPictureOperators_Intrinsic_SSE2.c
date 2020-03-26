/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPictureOperators_SSE2.h"
#include <emmintrin.h>
#include "EbDefinitions.h"

/******************************************************************************************************
                                       residual_kernel16bit_sse2_intrin
******************************************************************************************************/
void residual_kernel16bit_sse2_intrin(uint16_t *input, uint32_t input_stride, uint16_t *pred,
                                      uint32_t pred_stride, int16_t *residual,
                                      uint32_t residual_stride, uint32_t area_width,
                                      uint32_t area_height) {
    uint32_t x, y;
    __m128i  residual0, residual1;

    if (area_width == 4) {
        for (y = 0; y < area_height; y += 2) {
            residual0 =
                _mm_sub_epi16(_mm_loadl_epi64((__m128i *)input), _mm_loadl_epi64((__m128i *)pred));
            residual1 = _mm_sub_epi16(_mm_loadl_epi64((__m128i *)(input + input_stride)),
                                      _mm_loadl_epi64((__m128i *)(pred + pred_stride)));

            _mm_storel_epi64((__m128i *)residual, residual0);
            _mm_storel_epi64((__m128i *)(residual + residual_stride), residual1);

            input += input_stride << 1;
            pred += pred_stride << 1;
            residual += residual_stride << 1;
        }
    } else if (area_width == 8) {
        for (y = 0; y < area_height; y += 2) {
            residual0 =
                _mm_sub_epi16(_mm_loadu_si128((__m128i *)input), _mm_loadu_si128((__m128i *)pred));
            residual1 = _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride)),
                                      _mm_loadu_si128((__m128i *)(pred + pred_stride)));

            _mm_storeu_si128((__m128i *)residual, residual0);
            _mm_storeu_si128((__m128i *)(residual + residual_stride), residual1);

            input += input_stride << 1;
            pred += pred_stride << 1;
            residual += residual_stride << 1;
        }
    } else if (area_width == 16) {
        __m128i residual2, residual3;

        for (y = 0; y < area_height; y += 2) {
            residual0 =
                _mm_sub_epi16(_mm_loadu_si128((__m128i *)input), _mm_loadu_si128((__m128i *)pred));
            residual1 = _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 8)),
                                      _mm_loadu_si128((__m128i *)(pred + 8)));
            residual2 = _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride)),
                                      _mm_loadu_si128((__m128i *)(pred + pred_stride)));
            residual3 = _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 8)),
                                      _mm_loadu_si128((__m128i *)(pred + pred_stride + 8)));

            _mm_storeu_si128((__m128i *)residual, residual0);
            _mm_storeu_si128((__m128i *)(residual + 8), residual1);
            _mm_storeu_si128((__m128i *)(residual + residual_stride), residual2);
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 8), residual3);

            input += input_stride << 1;
            pred += pred_stride << 1;
            residual += residual_stride << 1;
        }
    } else if (area_width == 32) {
        for (y = 0; y < area_height; y += 2) {
            //residual[column_index] = ((int16_t)input[column_index]) - ((int16_t)pred[column_index]);
            _mm_storeu_si128(
                (__m128i *)residual,
                _mm_sub_epi16(_mm_loadu_si128((__m128i *)input), _mm_loadu_si128((__m128i *)pred)));
            _mm_storeu_si128((__m128i *)(residual + 8),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 8)),
                                           _mm_loadu_si128((__m128i *)(pred + 8))));
            _mm_storeu_si128((__m128i *)(residual + 16),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 16)),
                                           _mm_loadu_si128((__m128i *)(pred + 16))));
            _mm_storeu_si128((__m128i *)(residual + 24),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 24)),
                                           _mm_loadu_si128((__m128i *)(pred + 24))));

            _mm_storeu_si128((__m128i *)(residual + residual_stride),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 8),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 8)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 8))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 16),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 16)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 16))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 24),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 24)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 24))));

            input += input_stride << 1;
            pred += pred_stride << 1;
            residual += residual_stride << 1;
        }
    } else if (area_width ==
               64) { // Branch was not tested because the encoder had max txb_size of 32

        for (y = 0; y < area_height; y += 2) {
            //residual[column_index] = ((int16_t)input[column_index]) - ((int16_t)pred[column_index]) 8 indices per _mm_sub_epi16
            _mm_storeu_si128(
                (__m128i *)residual,
                _mm_sub_epi16(_mm_loadu_si128((__m128i *)input), _mm_loadu_si128((__m128i *)pred)));
            _mm_storeu_si128((__m128i *)(residual + 8),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 8)),
                                           _mm_loadu_si128((__m128i *)(pred + 8))));
            _mm_storeu_si128((__m128i *)(residual + 16),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 16)),
                                           _mm_loadu_si128((__m128i *)(pred + 16))));
            _mm_storeu_si128((__m128i *)(residual + 24),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 24)),
                                           _mm_loadu_si128((__m128i *)(pred + 24))));
            _mm_storeu_si128((__m128i *)(residual + 32),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 32)),
                                           _mm_loadu_si128((__m128i *)(pred + 32))));
            _mm_storeu_si128((__m128i *)(residual + 40),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 40)),
                                           _mm_loadu_si128((__m128i *)(pred + 40))));
            _mm_storeu_si128((__m128i *)(residual + 48),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 48)),
                                           _mm_loadu_si128((__m128i *)(pred + 48))));
            _mm_storeu_si128((__m128i *)(residual + 56),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + 56)),
                                           _mm_loadu_si128((__m128i *)(pred + 56))));

            _mm_storeu_si128((__m128i *)(residual + residual_stride),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 8),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 8)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 8))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 16),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 16)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 16))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 24),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 24)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 24))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 32),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 32)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 32))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 40),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 40)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 40))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 48),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 48)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 48))));
            _mm_storeu_si128((__m128i *)(residual + residual_stride + 56),
                             _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride + 56)),
                                           _mm_loadu_si128((__m128i *)(pred + pred_stride + 56))));

            input += input_stride << 1;
            pred += pred_stride << 1;
            residual += residual_stride << 1;
        }
    } else {
        uint32_t input_stride_diff    = 2 * input_stride;
        uint32_t pred_stride_diff     = 2 * pred_stride;
        uint32_t residual_stride_diff = 2 * residual_stride;
        input_stride_diff -= area_width;
        pred_stride_diff -= area_width;
        residual_stride_diff -= area_width;

        if (!(area_width & 7)) {
            for (x = 0; x < area_height; x += 2) {
                for (y = 0; y < area_width; y += 8) {
                    _mm_storeu_si128((__m128i *)residual,
                                     _mm_sub_epi16(_mm_loadu_si128((__m128i *)input),
                                                   _mm_loadu_si128((__m128i *)pred)));
                    _mm_storeu_si128(
                        (__m128i *)(residual + residual_stride),
                        _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride)),
                                      _mm_loadu_si128((__m128i *)(pred + pred_stride))));

                    input += 8;
                    pred += 8;
                    residual += 8;
                }
                input    = input + input_stride_diff;
                pred     = pred + pred_stride_diff;
                residual = residual + residual_stride_diff;
            }
        } else {
            for (x = 0; x < area_height; x += 2) {
                for (y = 0; y < area_width; y += 4) {
                    _mm_storel_epi64((__m128i *)residual,
                                     _mm_sub_epi16(_mm_loadu_si128((__m128i *)input),
                                                   _mm_loadu_si128((__m128i *)pred)));
                    _mm_storel_epi64(
                        (__m128i *)(residual + residual_stride),
                        _mm_sub_epi16(_mm_loadu_si128((__m128i *)(input + input_stride)),
                                      _mm_loadu_si128((__m128i *)(pred + pred_stride))));

                    input += 4;
                    pred += 4;
                    residual += 4;
                }
                input    = input + input_stride_diff;
                pred     = pred + pred_stride_diff;
                residual = residual + residual_stride_diff;
            }
        }
    }
    return;
}




