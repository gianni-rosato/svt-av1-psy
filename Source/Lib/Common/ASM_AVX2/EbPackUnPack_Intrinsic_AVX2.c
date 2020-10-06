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

#include <emmintrin.h>
#include <immintrin.h>
#include <stdint.h>

void svt_enc_un_pack8_bit_data_avx2_intrin(uint16_t *in_16bit_buffer, uint32_t in_stride,
                                           uint8_t *out_8bit_buffer, uint32_t out_stride,
                                           uint32_t width, uint32_t height) {
    __m256i ymm_00ff = _mm256_set1_epi16(0x00FF);
    __m128i xmm_00ff = _mm_set1_epi16(0x00FF);
    switch (width) {
    case 8:
        for (uint32_t y = 0; y < height; y += 2) {
            __m128i in_pixel0 = _mm_loadu_si128((__m128i *)in_16bit_buffer),
                    in_pixel1 = _mm_loadu_si128((__m128i *)(in_16bit_buffer + in_stride));

            __m128i in_pixel0_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff),
                    in_pixel1_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff);

            __m128i in_pixel0_shft_r_2_u8 = _mm_packus_epi16(in_pixel0_shft_r_2,
                                                             in_pixel0_shft_r_2),
                    in_pixel1_shftR_2_u8 = _mm_packus_epi16(in_pixel1_shft_r_2, in_pixel1_shft_r_2);

            _mm_storel_epi64((__m128i *)out_8bit_buffer, in_pixel0_shft_r_2_u8);
            _mm_storel_epi64((__m128i *)(out_8bit_buffer + out_stride), in_pixel1_shftR_2_u8);

            out_8bit_buffer += 2 * out_stride;
            in_16bit_buffer += 2 * in_stride;
        }
        return;
    case 16:
        for (uint32_t y = 0; y < height; y += 2) {
            __m256i in_pixel0 = _mm256_loadu_si256((__m256i *)in_16bit_buffer),
                    in_pixel1 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + in_stride));

            __m256i in_pixel0_shft_r_2_u8 = _mm256_packus_epi16(
                _mm256_and_si256(_mm256_srli_epi16(in_pixel0, 2), ymm_00ff),
                _mm256_and_si256(_mm256_srli_epi16(in_pixel1, 2), ymm_00ff));

            *(uint64_t *)out_8bit_buffer       = _mm256_extract_epi64(in_pixel0_shft_r_2_u8, 0);
            *(uint64_t *)(out_8bit_buffer + 8) = _mm256_extract_epi64(in_pixel0_shft_r_2_u8, 2);
            *(uint64_t *)(out_8bit_buffer + out_stride) = _mm256_extract_epi64(
                in_pixel0_shft_r_2_u8, 1);
            *(uint64_t *)(out_8bit_buffer + out_stride + 8) = _mm256_extract_epi64(
                in_pixel0_shft_r_2_u8, 3);

            out_8bit_buffer += 2 * out_stride;
            in_16bit_buffer += 2 * in_stride;
        }
        return;
    case 32:
        for (uint32_t y = 0; y < height; y += 2) {
            __m256i in_pixel0 = _mm256_loadu_si256((__m256i *)in_16bit_buffer),
                    in_pixel1 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + 16)),
                    in_pixel2 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + in_stride)),
                    in_pixel3 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + in_stride + 16));

            __m256i out8_0_u8 = _mm256_packus_epi16(
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel0, 2), ymm_00ff),
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel1, 2), ymm_00ff)),
                    out8_1_u8 = _mm256_packus_epi16(
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel2, 2), ymm_00ff),
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel3, 2), ymm_00ff));

            *(uint64_t *)out_8bit_buffer        = _mm256_extract_epi64(out8_0_u8, 0);
            *(uint64_t *)(out_8bit_buffer + 8)  = _mm256_extract_epi64(out8_0_u8, 2);
            *(uint64_t *)(out_8bit_buffer + 16) = _mm256_extract_epi64(out8_0_u8, 1);
            *(uint64_t *)(out_8bit_buffer + 24) = _mm256_extract_epi64(out8_0_u8, 3);
            out_8bit_buffer += out_stride;

            *(uint64_t *)out_8bit_buffer        = _mm256_extract_epi64(out8_1_u8, 0);
            *(uint64_t *)(out_8bit_buffer + 8)  = _mm256_extract_epi64(out8_1_u8, 2);
            *(uint64_t *)(out_8bit_buffer + 16) = _mm256_extract_epi64(out8_1_u8, 1);
            *(uint64_t *)(out_8bit_buffer + 24) = _mm256_extract_epi64(out8_1_u8, 3);
            out_8bit_buffer += out_stride;
            in_16bit_buffer += 2 * in_stride;
        }
        return;
    case 64:
        for (uint32_t y = 0; y < height; ++y) {
            __m256i in_pixel0 = _mm256_loadu_si256((__m256i *)in_16bit_buffer),
                    in_pixel1 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + 16)),
                    in_pixel2 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + 32)),
                    in_pixel3 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + 48));

            __m256i out8_0_u8 = _mm256_packus_epi16(
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel0, 2), ymm_00ff),
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel1, 2), ymm_00ff)),
                    out8_1_u8 = _mm256_packus_epi16(
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel2, 2), ymm_00ff),
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel3, 2), ymm_00ff));

            *(uint64_t *)out_8bit_buffer        = _mm256_extract_epi64(out8_0_u8, 0);
            *(uint64_t *)(out_8bit_buffer + 8)  = _mm256_extract_epi64(out8_0_u8, 2);
            *(uint64_t *)(out_8bit_buffer + 16) = _mm256_extract_epi64(out8_0_u8, 1);
            *(uint64_t *)(out_8bit_buffer + 24) = _mm256_extract_epi64(out8_0_u8, 3);
            *(uint64_t *)(out_8bit_buffer + 32) = _mm256_extract_epi64(out8_1_u8, 0);
            *(uint64_t *)(out_8bit_buffer + 40) = _mm256_extract_epi64(out8_1_u8, 2);
            *(uint64_t *)(out_8bit_buffer + 48) = _mm256_extract_epi64(out8_1_u8, 1);
            *(uint64_t *)(out_8bit_buffer + 56) = _mm256_extract_epi64(out8_1_u8, 3);

            out_8bit_buffer += out_stride;
            in_16bit_buffer += in_stride;
        }
        return;
    default:
        if (!(width & 63)) {
            uint32_t in_stride_diff64 = in_stride - width;
            uint32_t out_stride_diff64 = out_stride - width;
            for (uint32_t x = 0; x < height; x += 1) {
                for (uint32_t y = 0; y < width; y += 64) {
                    __m256i in_pixel0 = _mm256_loadu_si256((__m256i *)in_16bit_buffer),
                            in_pixel1 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + 16)),
                            in_pixel2 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + 32)),
                            in_pixel3 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + 48));

                    __m256i out8_0_u8 = _mm256_packus_epi16(
                                _mm256_and_si256(_mm256_srli_epi16(in_pixel0, 2), ymm_00ff),
                                _mm256_and_si256(_mm256_srli_epi16(in_pixel1, 2), ymm_00ff)),
                            out8_1_u8 = _mm256_packus_epi16(
                                _mm256_and_si256(_mm256_srli_epi16(in_pixel2, 2), ymm_00ff),
                                _mm256_and_si256(_mm256_srli_epi16(in_pixel3, 2), ymm_00ff));

                    *(uint64_t *)out_8bit_buffer        = _mm256_extract_epi64(out8_0_u8, 0);
                    *(uint64_t *)(out_8bit_buffer + 8)  = _mm256_extract_epi64(out8_0_u8, 2);
                    *(uint64_t *)(out_8bit_buffer + 16) = _mm256_extract_epi64(out8_0_u8, 1);
                    *(uint64_t *)(out_8bit_buffer + 24) = _mm256_extract_epi64(out8_0_u8, 3);
                    *(uint64_t *)(out_8bit_buffer + 32) = _mm256_extract_epi64(out8_1_u8, 0);
                    *(uint64_t *)(out_8bit_buffer + 40) = _mm256_extract_epi64(out8_1_u8, 2);
                    *(uint64_t *)(out_8bit_buffer + 48) = _mm256_extract_epi64(out8_1_u8, 1);
                    *(uint64_t *)(out_8bit_buffer + 56) = _mm256_extract_epi64(out8_1_u8, 3);

                    out_8bit_buffer += 64;
                    in_16bit_buffer += 64;
                }
                in_16bit_buffer += in_stride_diff64;
                out_8bit_buffer += out_stride_diff64;
            }
        } else if (!(width & 31)) {
            uint32_t in_stride_diff = (2 * in_stride) - width;
            uint32_t out_stride_diff = (2 * out_stride) - width;
            for (uint32_t x = 0; x < height; x += 2) {
                for (uint32_t y = 0; y < width; y += 32) {
                    __m256i in_pixel0 = _mm256_loadu_si256((__m256i *)in_16bit_buffer),
                            in_pixel1 = _mm256_loadu_si256((__m256i *)(in_16bit_buffer + 16)),
                            in_pixel2 = _mm256_loadu_si256(
                                (__m256i *)(in_16bit_buffer + in_stride)),
                            in_pixel3 = _mm256_loadu_si256(
                                (__m256i *)(in_16bit_buffer + in_stride + 16));

                    __m256i out8_0_u8 = _mm256_packus_epi16(
                                _mm256_and_si256(_mm256_srli_epi16(in_pixel0, 2), ymm_00ff),
                                _mm256_and_si256(_mm256_srli_epi16(in_pixel1, 2), ymm_00ff)),
                            out8_1_u8 = _mm256_packus_epi16(
                                _mm256_and_si256(_mm256_srli_epi16(in_pixel2, 2), ymm_00ff),
                                _mm256_and_si256(_mm256_srli_epi16(in_pixel3, 2), ymm_00ff));

                    *(uint64_t *)out_8bit_buffer        = _mm256_extract_epi64(out8_0_u8, 0);
                    *(uint64_t *)(out_8bit_buffer + 8)  = _mm256_extract_epi64(out8_0_u8, 2);
                    *(uint64_t *)(out_8bit_buffer + 16) = _mm256_extract_epi64(out8_0_u8, 1);
                    *(uint64_t *)(out_8bit_buffer + 24) = _mm256_extract_epi64(out8_0_u8, 3);

                    *(uint64_t *)(out_8bit_buffer + out_stride) = _mm256_extract_epi64(out8_1_u8,
                                                                                       0);
                    *(uint64_t *)(out_8bit_buffer + out_stride + 8) = _mm256_extract_epi64(
                        out8_1_u8, 2);
                    *(uint64_t *)(out_8bit_buffer + out_stride + 16) = _mm256_extract_epi64(
                        out8_1_u8, 1);
                    *(uint64_t *)(out_8bit_buffer + out_stride + 24) = _mm256_extract_epi64(
                        out8_1_u8, 3);

                    out_8bit_buffer += 32;
                    in_16bit_buffer += 32;
                }
                in_16bit_buffer += in_stride_diff;
                out_8bit_buffer += out_stride_diff;
            }
        } else if (!(width & 15)) {
            uint32_t in_stride_diff = (2 * in_stride) - width;
            uint32_t out_stride_diff = (2 * out_stride) - width;
            for (uint32_t x = 0; x < height; x += 2) {
                for (uint32_t y = 0; y < width; y += 16) {
                    __m256i in_pixel0 = _mm256_loadu_si256((__m256i *)in_16bit_buffer),
                            in_pixel1 = _mm256_loadu_si256(
                                (__m256i *)(in_16bit_buffer + in_stride));

                    __m256i in_pixel0_shft_r_2_u8 = _mm256_packus_epi16(
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel0, 2), ymm_00ff),
                        _mm256_and_si256(_mm256_srli_epi16(in_pixel1, 2), ymm_00ff));

                    *(uint64_t *)out_8bit_buffer = _mm256_extract_epi64(in_pixel0_shft_r_2_u8, 0);
                    *(uint64_t *)(out_8bit_buffer + 8) = _mm256_extract_epi64(in_pixel0_shft_r_2_u8,
                                                                              2);
                    *(uint64_t *)(out_8bit_buffer + out_stride) = _mm256_extract_epi64(
                        in_pixel0_shft_r_2_u8, 1);
                    *(uint64_t *)(out_8bit_buffer + out_stride + 8) = _mm256_extract_epi64(
                        in_pixel0_shft_r_2_u8, 3);

                    out_8bit_buffer += 16;
                    in_16bit_buffer += 16;
                }
                in_16bit_buffer += in_stride_diff;
                out_8bit_buffer += out_stride_diff;
            }
        } else if (!(width & 7)) {
            uint32_t in_stride_diff = (2 * in_stride) - width;
            uint32_t out_stride_diff = (2 * out_stride) - width;
            for (uint32_t x = 0; x < height; x += 2) {
                for (uint32_t y = 0; y < width; y += 8) {
                    __m128i in_pixel0 = _mm_loadu_si128((__m128i *)in_16bit_buffer),
                            in_pixel1 = _mm_loadu_si128((__m128i *)(in_16bit_buffer + in_stride));

                    __m128i in_pixel0_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel0, 2),
                                                               xmm_00ff),
                            in_pixel1_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel1, 2),
                                                               xmm_00ff);

                    __m128i in_pixel0_shft_r_2_u8 = _mm_packus_epi16(in_pixel0_shft_r_2,
                                                                     in_pixel0_shft_r_2),
                            in_pixel1_shftR_2_u8  = _mm_packus_epi16(in_pixel1_shft_r_2,
                                                                    in_pixel1_shft_r_2);

                    _mm_storel_epi64((__m128i *)out_8bit_buffer, in_pixel0_shft_r_2_u8);
                    _mm_storel_epi64((__m128i *)(out_8bit_buffer + out_stride),
                                     in_pixel1_shftR_2_u8);

                    out_8bit_buffer += 8;
                    in_16bit_buffer += 8;
                }
                in_16bit_buffer += in_stride_diff;
                out_8bit_buffer += out_stride_diff;
            }
        } else {
            uint32_t in_stride_diff = (2 * in_stride) - width;
            uint32_t out_stride_diff = (2 * out_stride) - width;
            uint32_t width_down4 = width & (~0x3);
            for (uint32_t x = 0; x < height; x += 2) {
                uint32_t y = 0;
                for (; y < width_down4; y += 4) {
                    __m128i in_pixel0 = _mm_loadl_epi64((__m128i *)in_16bit_buffer),
                            in_pixel1 = _mm_loadl_epi64((__m128i *)(in_16bit_buffer + in_stride));

                    __m128i in_pixel0_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel0, 2),
                                                               xmm_00ff),
                            in_pixel1_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel1, 2),
                                                               xmm_00ff);

                    __m128i in_pixel0_shft_r_2_u8 = _mm_packus_epi16(in_pixel0_shft_r_2,
                                                                     in_pixel0_shft_r_2),
                            in_pixel1_shftR_2_u8  = _mm_packus_epi16(in_pixel1_shft_r_2,
                                                                    in_pixel1_shft_r_2);

                    *(uint32_t *)out_8bit_buffer = _mm_cvtsi128_si32(in_pixel0_shft_r_2_u8);
                    *(uint32_t *)(out_8bit_buffer +
                                  out_stride)    = _mm_cvtsi128_si32(in_pixel1_shftR_2_u8);

                    out_8bit_buffer += 4;
                    in_16bit_buffer += 4;
                }

                /* Calculate lefts pixels in 2 lines,
                 * when width is not divided by 4.
                 */
                for (; y < width; y++) {
                    uint16_t in_pixel               = *in_16bit_buffer;
                    *out_8bit_buffer                = (uint8_t)(in_pixel >> 2);
                    in_pixel                        = *(in_16bit_buffer + in_stride);
                    *(out_8bit_buffer + out_stride) = (uint8_t)(in_pixel >> 2);
                    ++out_8bit_buffer;
                    ++in_16bit_buffer;
                }

                in_16bit_buffer += in_stride_diff;
                out_8bit_buffer += out_stride_diff;
            }
        }
    }
}
