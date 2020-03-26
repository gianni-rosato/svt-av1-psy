/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <emmintrin.h>
#include <stdint.h>

/****************************************************************************************
eb_enc_msb_un_pack2d_sse2_intrin
******************************************************************************************/

void eb_enc_msb_un_pack2d_sse2_intrin(uint16_t *in16_bit_buffer, uint32_t in_stride,
                                      uint8_t *out8_bit_buffer, uint8_t *outn_bit_buffer,
                                      uint32_t out8_stride, uint32_t outn_stride, uint32_t width,
                                      uint32_t height) {
    uint32_t x, y;

    __m128i xmm_3, xmm_00ff, in_pixel0, in_pixel1, temp_pixel0, temp_pixel1, in_pixel1_shft_r_2_u8,
        in_pixel0_shft_r_2_u8, in_pixel0_shft_r_2, in_pixel1_shft_r_2;
    __m128i temp_pixel0_u8, temp_pixel1_u8;

    xmm_3    = _mm_set1_epi16(0x0003);
    xmm_00ff = _mm_set1_epi16(0x00FF);

    if (width == 4) {
        for (y = 0; y < height; y += 2) {
            in_pixel0 = _mm_loadl_epi64((__m128i *)in16_bit_buffer);
            in_pixel1 = _mm_loadl_epi64((__m128i *)(in16_bit_buffer + in_stride));

            temp_pixel0 = _mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6);
            temp_pixel1 = _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6);

            temp_pixel0_u8 = _mm_packus_epi16(temp_pixel0, temp_pixel0);
            temp_pixel1_u8 = _mm_packus_epi16(temp_pixel1, temp_pixel1);

            in_pixel0_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff);
            in_pixel1_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff);

            in_pixel0_shft_r_2_u8 = _mm_packus_epi16(in_pixel0_shft_r_2, in_pixel0_shft_r_2);
            in_pixel1_shft_r_2_u8 = _mm_packus_epi16(in_pixel1_shft_r_2, in_pixel1_shft_r_2);

            *(uint32_t *)outn_bit_buffer                 = _mm_cvtsi128_si32(temp_pixel0_u8);
            *(uint32_t *)(outn_bit_buffer + outn_stride) = _mm_cvtsi128_si32(temp_pixel1_u8);
            *(uint32_t *)out8_bit_buffer                 = _mm_cvtsi128_si32(in_pixel0_shft_r_2_u8);
            *(uint32_t *)(out8_bit_buffer + out8_stride) = _mm_cvtsi128_si32(in_pixel1_shft_r_2_u8);

            outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    } else if (width == 8) {
        for (y = 0; y < height; y += 2) {
            in_pixel0 = _mm_loadu_si128((__m128i *)in16_bit_buffer);
            in_pixel1 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride));

            temp_pixel0 = _mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6);
            temp_pixel1 = _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6);

            temp_pixel0_u8 = _mm_packus_epi16(temp_pixel0, temp_pixel0);
            temp_pixel1_u8 = _mm_packus_epi16(temp_pixel1, temp_pixel1);

            in_pixel0_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff);
            in_pixel1_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff);

            in_pixel0_shft_r_2_u8 = _mm_packus_epi16(in_pixel0_shft_r_2, in_pixel0_shft_r_2);
            in_pixel1_shft_r_2_u8 = _mm_packus_epi16(in_pixel1_shft_r_2, in_pixel1_shft_r_2);

            _mm_storel_epi64((__m128i *)outn_bit_buffer, temp_pixel0_u8);
            _mm_storel_epi64((__m128i *)(outn_bit_buffer + outn_stride), temp_pixel1_u8);
            _mm_storel_epi64((__m128i *)out8_bit_buffer, in_pixel0_shft_r_2_u8);
            _mm_storel_epi64((__m128i *)(out8_bit_buffer + out8_stride), in_pixel1_shft_r_2_u8);

            outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    } else if (width == 16) {
        __m128i in_pixel2, in_pixel3;

        for (y = 0; y < height; y += 2) {
            in_pixel0 = _mm_loadu_si128((__m128i *)in16_bit_buffer);
            in_pixel1 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 8));
            in_pixel2 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride));
            in_pixel3 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride + 8));

            temp_pixel0_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6),
                                             _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6));
            temp_pixel1_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel2, xmm_3), 6),
                                             _mm_slli_epi16(_mm_and_si128(in_pixel3, xmm_3), 6));

            in_pixel0_shft_r_2_u8 =
                _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff),
                                 _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff));
            in_pixel1_shft_r_2_u8 =
                _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel2, 2), xmm_00ff),
                                 _mm_and_si128(_mm_srli_epi16(in_pixel3, 2), xmm_00ff));

            _mm_storeu_si128((__m128i *)outn_bit_buffer, temp_pixel0_u8);
            _mm_storeu_si128((__m128i *)(outn_bit_buffer + outn_stride), temp_pixel1_u8);
            _mm_storeu_si128((__m128i *)out8_bit_buffer, in_pixel0_shft_r_2_u8);
            _mm_storeu_si128((__m128i *)(out8_bit_buffer + out8_stride), in_pixel1_shft_r_2_u8);

            outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    } else if (width == 32) {
        __m128i in_pixel2, in_pixel3, in_pixel4, in_pixel5, in_pixel6, in_pixel7;
        __m128i outn0_u8, outn1_u8, outn2_u8, outn3_u8, out8_0_u8, out8_1_u8, out8_2_u8, out8_3_u8;

        for (y = 0; y < height; y += 2) {
            in_pixel0 = _mm_loadu_si128((__m128i *)in16_bit_buffer);
            in_pixel1 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 8));
            in_pixel2 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 16));
            in_pixel3 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 24));
            in_pixel4 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride));
            in_pixel5 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride + 8));
            in_pixel6 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride + 16));
            in_pixel7 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride + 24));

            outn0_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6),
                                        _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6));
            outn1_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel2, xmm_3), 6),
                                        _mm_slli_epi16(_mm_and_si128(in_pixel3, xmm_3), 6));
            outn2_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel4, xmm_3), 6),
                                        _mm_slli_epi16(_mm_and_si128(in_pixel5, xmm_3), 6));
            outn3_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel6, xmm_3), 6),
                                        _mm_slli_epi16(_mm_and_si128(in_pixel7, xmm_3), 6));

            out8_0_u8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff));
            out8_1_u8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel2, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel3, 2), xmm_00ff));
            out8_2_u8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel4, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel5, 2), xmm_00ff));
            out8_3_u8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel6, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel7, 2), xmm_00ff));

            _mm_storeu_si128((__m128i *)outn_bit_buffer, outn0_u8);
            _mm_storeu_si128((__m128i *)(outn_bit_buffer + 16), outn1_u8);
            _mm_storeu_si128((__m128i *)(outn_bit_buffer + outn_stride), outn2_u8);
            _mm_storeu_si128((__m128i *)(outn_bit_buffer + outn_stride + 16), outn3_u8);

            _mm_storeu_si128((__m128i *)out8_bit_buffer, out8_0_u8);
            _mm_storeu_si128((__m128i *)(out8_bit_buffer + 16), out8_1_u8);
            _mm_storeu_si128((__m128i *)(out8_bit_buffer + out8_stride), out8_2_u8);
            _mm_storeu_si128((__m128i *)(out8_bit_buffer + out8_stride + 16), out8_3_u8);

            outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    } else if (width == 64) {
        __m128i in_pixel2, in_pixel3, in_pixel4, in_pixel5, in_pixel6, in_pixel7;
        __m128i outn0_u8, outn1_u8, outn2_u8, outn3_u8, out8_0_u8, out8_1_u8, out8_2_u8, out8_3_u8;

        for (y = 0; y < height; ++y) {
            in_pixel0 = _mm_loadu_si128((__m128i *)in16_bit_buffer);
            in_pixel1 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 8));
            in_pixel2 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 16));
            in_pixel3 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 24));
            in_pixel4 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 32));
            in_pixel5 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 40));
            in_pixel6 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 48));
            in_pixel7 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 56));

            outn0_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6),
                                        _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6));
            outn1_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel2, xmm_3), 6),
                                        _mm_slli_epi16(_mm_and_si128(in_pixel3, xmm_3), 6));
            outn2_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel4, xmm_3), 6),
                                        _mm_slli_epi16(_mm_and_si128(in_pixel5, xmm_3), 6));
            outn3_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel6, xmm_3), 6),
                                        _mm_slli_epi16(_mm_and_si128(in_pixel7, xmm_3), 6));

            out8_0_u8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff));
            out8_1_u8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel2, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel3, 2), xmm_00ff));
            out8_2_u8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel4, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel5, 2), xmm_00ff));
            out8_3_u8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel6, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel7, 2), xmm_00ff));

            _mm_storeu_si128((__m128i *)outn_bit_buffer, outn0_u8);
            _mm_storeu_si128((__m128i *)(outn_bit_buffer + 16), outn1_u8);
            _mm_storeu_si128((__m128i *)(outn_bit_buffer + 32), outn2_u8);
            _mm_storeu_si128((__m128i *)(outn_bit_buffer + 48), outn3_u8);

            _mm_storeu_si128((__m128i *)out8_bit_buffer, out8_0_u8);
            _mm_storeu_si128((__m128i *)(out8_bit_buffer + 16), out8_1_u8);
            _mm_storeu_si128((__m128i *)(out8_bit_buffer + 32), out8_2_u8);
            _mm_storeu_si128((__m128i *)(out8_bit_buffer + 48), out8_3_u8);

            outn_bit_buffer += outn_stride;
            out8_bit_buffer += out8_stride;
            in16_bit_buffer += in_stride;
        }
    } else {
        uint32_t in_stride_diff    = (2 * in_stride) - width;
        uint32_t out8_stride_diff  = (2 * out8_stride) - width;
        uint32_t out_n_stride_diff = (2 * outn_stride) - width;

        uint32_t in_stride_diff64    = in_stride - width;
        uint32_t out8_stride_diff64  = out8_stride - width;
        uint32_t out_n_stride_diff64 = outn_stride - width;

        if (!(width & 63)) {
            __m128i in_pixel2, in_pixel3, in_pixel4, in_pixel5, in_pixel6, in_pixel7;
            __m128i outn0_u8, outn1_u8, outn2_u8, outn3_u8, out8_0_u8, out8_1_u8, out8_2_u8,
                out8_3_u8;

            for (x = 0; x < height; x += 1) {
                for (y = 0; y < width; y += 64) {
                    in_pixel0 = _mm_loadu_si128((__m128i *)in16_bit_buffer);
                    in_pixel1 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 8));
                    in_pixel2 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 16));
                    in_pixel3 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 24));
                    in_pixel4 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 32));
                    in_pixel5 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 40));
                    in_pixel6 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 48));
                    in_pixel7 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 56));

                    outn0_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6),
                                                _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6));
                    outn1_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel2, xmm_3), 6),
                                                _mm_slli_epi16(_mm_and_si128(in_pixel3, xmm_3), 6));
                    outn2_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel4, xmm_3), 6),
                                                _mm_slli_epi16(_mm_and_si128(in_pixel5, xmm_3), 6));
                    outn3_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel6, xmm_3), 6),
                                                _mm_slli_epi16(_mm_and_si128(in_pixel7, xmm_3), 6));

                    out8_0_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff));
                    out8_1_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel2, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel3, 2), xmm_00ff));
                    out8_2_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel4, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel5, 2), xmm_00ff));
                    out8_3_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel6, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel7, 2), xmm_00ff));

                    _mm_storeu_si128((__m128i *)outn_bit_buffer, outn0_u8);
                    _mm_storeu_si128((__m128i *)(outn_bit_buffer + 16), outn1_u8);
                    _mm_storeu_si128((__m128i *)(outn_bit_buffer + 32), outn2_u8);
                    _mm_storeu_si128((__m128i *)(outn_bit_buffer + 48), outn3_u8);

                    _mm_storeu_si128((__m128i *)out8_bit_buffer, out8_0_u8);
                    _mm_storeu_si128((__m128i *)(out8_bit_buffer + 16), out8_1_u8);
                    _mm_storeu_si128((__m128i *)(out8_bit_buffer + 32), out8_2_u8);
                    _mm_storeu_si128((__m128i *)(out8_bit_buffer + 48), out8_3_u8);

                    outn_bit_buffer += 64;
                    out8_bit_buffer += 64;
                    in16_bit_buffer += 64;
                }
                in16_bit_buffer += in_stride_diff64;
                outn_bit_buffer += out_n_stride_diff64;
                out8_bit_buffer += out8_stride_diff64;
            }
        } else if (!(width & 31)) {
            __m128i in_pixel2, in_pixel3, in_pixel4, in_pixel5, in_pixel6, in_pixel7;
            __m128i outn0_u8, outn1_u8, outn2_u8, outn3_u8, out8_0_u8, out8_1_u8, out8_2_u8,
                out8_3_u8;

            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 32) {
                    in_pixel0 = _mm_loadu_si128((__m128i *)in16_bit_buffer);
                    in_pixel1 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 8));
                    in_pixel2 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 16));
                    in_pixel3 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 24));
                    in_pixel4 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride));
                    in_pixel5 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride + 8));
                    in_pixel6 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride + 16));
                    in_pixel7 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride + 24));

                    outn0_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6),
                                                _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6));
                    outn1_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel2, xmm_3), 6),
                                                _mm_slli_epi16(_mm_and_si128(in_pixel3, xmm_3), 6));
                    outn2_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel4, xmm_3), 6),
                                                _mm_slli_epi16(_mm_and_si128(in_pixel5, xmm_3), 6));
                    outn3_u8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel6, xmm_3), 6),
                                                _mm_slli_epi16(_mm_and_si128(in_pixel7, xmm_3), 6));

                    out8_0_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff));
                    out8_1_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel2, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel3, 2), xmm_00ff));
                    out8_2_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel4, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel5, 2), xmm_00ff));
                    out8_3_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel6, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel7, 2), xmm_00ff));

                    _mm_storeu_si128((__m128i *)outn_bit_buffer, outn0_u8);
                    _mm_storeu_si128((__m128i *)(outn_bit_buffer + 16), outn1_u8);
                    _mm_storeu_si128((__m128i *)(outn_bit_buffer + outn_stride), outn2_u8);
                    _mm_storeu_si128((__m128i *)(outn_bit_buffer + outn_stride + 16), outn3_u8);

                    _mm_storeu_si128((__m128i *)out8_bit_buffer, out8_0_u8);
                    _mm_storeu_si128((__m128i *)(out8_bit_buffer + 16), out8_1_u8);
                    _mm_storeu_si128((__m128i *)(out8_bit_buffer + out8_stride), out8_2_u8);
                    _mm_storeu_si128((__m128i *)(out8_bit_buffer + out8_stride + 16), out8_3_u8);

                    outn_bit_buffer += 32;
                    out8_bit_buffer += 32;
                    in16_bit_buffer += 32;
                }
                in16_bit_buffer += in_stride_diff;
                outn_bit_buffer += out_n_stride_diff;
                out8_bit_buffer += out8_stride_diff;
            }
        } else if (!(width & 15)) {
            __m128i in_pixel2, in_pixel3;

            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 16) {
                    in_pixel0 = _mm_loadu_si128((__m128i *)in16_bit_buffer);
                    in_pixel1 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + 8));
                    in_pixel2 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride));
                    in_pixel3 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride + 8));

                    temp_pixel0_u8 =
                        _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6),
                                         _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6));
                    temp_pixel1_u8 =
                        _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(in_pixel2, xmm_3), 6),
                                         _mm_slli_epi16(_mm_and_si128(in_pixel3, xmm_3), 6));

                    in_pixel0_shft_r_2_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff));
                    in_pixel1_shft_r_2_u8 =
                        _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(in_pixel2, 2), xmm_00ff),
                                         _mm_and_si128(_mm_srli_epi16(in_pixel3, 2), xmm_00ff));

                    _mm_storeu_si128((__m128i *)outn_bit_buffer, temp_pixel0_u8);
                    _mm_storeu_si128((__m128i *)(outn_bit_buffer + outn_stride), temp_pixel1_u8);
                    _mm_storeu_si128((__m128i *)out8_bit_buffer, in_pixel0_shft_r_2_u8);
                    _mm_storeu_si128((__m128i *)(out8_bit_buffer + out8_stride),
                                     in_pixel1_shft_r_2_u8);

                    outn_bit_buffer += 16;
                    out8_bit_buffer += 16;
                    in16_bit_buffer += 16;
                }
                in16_bit_buffer += in_stride_diff;
                outn_bit_buffer += out_n_stride_diff;
                out8_bit_buffer += out8_stride_diff;
            }
        } else if (!(width & 7)) {
            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 8) {
                    in_pixel0 = _mm_loadu_si128((__m128i *)in16_bit_buffer);
                    in_pixel1 = _mm_loadu_si128((__m128i *)(in16_bit_buffer + in_stride));

                    temp_pixel0 = _mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6);
                    temp_pixel1 = _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6);

                    temp_pixel0_u8 = _mm_packus_epi16(temp_pixel0, temp_pixel0);
                    temp_pixel1_u8 = _mm_packus_epi16(temp_pixel1, temp_pixel1);

                    in_pixel0_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff);
                    in_pixel1_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff);

                    in_pixel0_shft_r_2_u8 =
                        _mm_packus_epi16(in_pixel0_shft_r_2, in_pixel0_shft_r_2);
                    in_pixel1_shft_r_2_u8 =
                        _mm_packus_epi16(in_pixel1_shft_r_2, in_pixel1_shft_r_2);

                    _mm_storel_epi64((__m128i *)outn_bit_buffer, temp_pixel0_u8);
                    _mm_storel_epi64((__m128i *)(outn_bit_buffer + outn_stride), temp_pixel1_u8);
                    _mm_storel_epi64((__m128i *)out8_bit_buffer, in_pixel0_shft_r_2_u8);
                    _mm_storel_epi64((__m128i *)(out8_bit_buffer + out8_stride),
                                     in_pixel1_shft_r_2_u8);

                    outn_bit_buffer += 8;
                    out8_bit_buffer += 8;
                    in16_bit_buffer += 8;
                }
                in16_bit_buffer += in_stride_diff;
                outn_bit_buffer += out_n_stride_diff;
                out8_bit_buffer += out8_stride_diff;
            }
        } else {
            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 4) {
                    in_pixel0 = _mm_loadl_epi64((__m128i *)in16_bit_buffer);
                    in_pixel1 = _mm_loadl_epi64((__m128i *)(in16_bit_buffer + in_stride));

                    temp_pixel0 = _mm_slli_epi16(_mm_and_si128(in_pixel0, xmm_3), 6);
                    temp_pixel1 = _mm_slli_epi16(_mm_and_si128(in_pixel1, xmm_3), 6);

                    temp_pixel0_u8 = _mm_packus_epi16(temp_pixel0, temp_pixel0);
                    temp_pixel1_u8 = _mm_packus_epi16(temp_pixel1, temp_pixel1);

                    in_pixel0_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel0, 2), xmm_00ff);
                    in_pixel1_shft_r_2 = _mm_and_si128(_mm_srli_epi16(in_pixel1, 2), xmm_00ff);

                    in_pixel0_shft_r_2_u8 =
                        _mm_packus_epi16(in_pixel0_shft_r_2, in_pixel0_shft_r_2);
                    in_pixel1_shft_r_2_u8 =
                        _mm_packus_epi16(in_pixel1_shft_r_2, in_pixel1_shft_r_2);

                    *(uint32_t *)outn_bit_buffer                 = _mm_cvtsi128_si32(temp_pixel0_u8);
                    *(uint32_t *)(outn_bit_buffer + outn_stride) = _mm_cvtsi128_si32(temp_pixel1_u8);
                    *(uint32_t *)out8_bit_buffer = _mm_cvtsi128_si32(in_pixel0_shft_r_2_u8);
                    *(uint32_t *)(out8_bit_buffer + out8_stride) =
                        _mm_cvtsi128_si32(in_pixel1_shft_r_2_u8);

                    outn_bit_buffer += 4;
                    out8_bit_buffer += 4;
                    in16_bit_buffer += 4;
                }
                in16_bit_buffer += in_stride_diff;
                outn_bit_buffer += out_n_stride_diff;
                out8_bit_buffer += out8_stride_diff;
            }
        }
    }
    return;
}

void unpack_avg_sse2_intrin(uint16_t *ref16_l0, uint32_t ref_l0_stride, uint16_t *ref16_l1,
                            uint32_t ref_l1_stride, uint8_t *dst_ptr, uint32_t dst_stride,
                            uint32_t width, uint32_t height) {
    uint32_t y;
    __m128i  in_pixel0, in_pixel1;

    if (width == 4) {
        __m128i out8_0_u8_l0, out8_0_u8_l1;
        __m128i avg8_0_u8;

        for (y = 0; y < height; y += 2) {
            //--------
            //Line One
            //--------

            //List0
            in_pixel0    = _mm_loadl_epi64((__m128i *)ref16_l0);
            in_pixel1    = _mm_srli_epi16(in_pixel0, 2);
            out8_0_u8_l0 = _mm_packus_epi16(in_pixel1, in_pixel1);

            //List1
            in_pixel0    = _mm_loadl_epi64((__m128i *)ref16_l1);
            in_pixel1    = _mm_srli_epi16(in_pixel0, 2);
            out8_0_u8_l1 = _mm_packus_epi16(in_pixel1, in_pixel1);

            //AVG
            avg8_0_u8 = _mm_avg_epu8(out8_0_u8_l0, out8_0_u8_l1);

            *(uint32_t *)dst_ptr = _mm_cvtsi128_si32(avg8_0_u8);

            //--------
            //Line Two
            //--------

            //List0
            in_pixel0    = _mm_loadl_epi64((__m128i *)(ref16_l0 + ref_l0_stride));
            in_pixel1    = _mm_srli_epi16(in_pixel0, 2);
            out8_0_u8_l0 = _mm_packus_epi16(in_pixel1, in_pixel1);

            //List1

            in_pixel0    = _mm_loadl_epi64((__m128i *)(ref16_l1 + ref_l1_stride));
            in_pixel1    = _mm_srli_epi16(in_pixel0, 2);
            out8_0_u8_l1 = _mm_packus_epi16(in_pixel1, in_pixel1);

            //AVG
            avg8_0_u8 = _mm_avg_epu8(out8_0_u8_l0, out8_0_u8_l1);

            *(uint32_t *)(dst_ptr + dst_stride) = _mm_cvtsi128_si32(avg8_0_u8);

            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;
        }
    } else if (width == 8) {
        __m128i out8_0_u8_l0, out8_0_u8_l1, out8_2_u8_l0, out8_2_u8_l1;
        __m128i avg8_0_u8, avg8_2_u8;

        for (y = 0; y < height; y += 2) {
            //--------
            //Line One
            //--------

            //List0

            in_pixel0 = _mm_loadu_si128((__m128i *)ref16_l0);

            in_pixel1    = _mm_srli_epi16(in_pixel0, 2);
            out8_0_u8_l0 = _mm_packus_epi16(in_pixel1, in_pixel1);

            //List1

            in_pixel0 = _mm_loadu_si128((__m128i *)ref16_l1);

            in_pixel1    = _mm_srli_epi16(in_pixel0, 2);
            out8_0_u8_l1 = _mm_packus_epi16(in_pixel1, in_pixel1);

            //AVG
            avg8_0_u8 = _mm_avg_epu8(out8_0_u8_l0, out8_0_u8_l1);

            _mm_storel_epi64((__m128i *)dst_ptr, avg8_0_u8);

            //--------
            //Line Two
            //--------

            //List0

            in_pixel0 = _mm_loadu_si128((__m128i *)(ref16_l0 + ref_l0_stride));

            in_pixel1    = _mm_srli_epi16(in_pixel0, 2);
            out8_2_u8_l0 = _mm_packus_epi16(in_pixel1, in_pixel1);

            //List1

            in_pixel0 = _mm_loadu_si128((__m128i *)(ref16_l1 + ref_l1_stride));

            in_pixel1    = _mm_srli_epi16(in_pixel0, 2);
            out8_2_u8_l1 = _mm_packus_epi16(in_pixel1, in_pixel1);

            //AVG
            avg8_2_u8 = _mm_avg_epu8(out8_2_u8_l0, out8_2_u8_l1);

            _mm_storel_epi64((__m128i *)(dst_ptr + dst_stride), avg8_2_u8);

            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;
        }
    } else if (width == 16) {
        __m128i in_pixel4, in_pixel5;
        __m128i out8_0_u8_l0, out8_0_u8_l1, out8_2_u8_l0, out8_2_u8_l1;
        __m128i avg8_0_u8, avg8_2_u8;

        for (y = 0; y < height; y += 2) {
            //--------
            //Line One
            //--------

            //List0

            in_pixel0 = _mm_loadu_si128((__m128i *)ref16_l0);
            in_pixel1 = _mm_loadu_si128((__m128i *)(ref16_l0 + 8));

            out8_0_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel0, 2), _mm_srli_epi16(in_pixel1, 2));

            //List1

            in_pixel0 = _mm_loadu_si128((__m128i *)ref16_l1);
            in_pixel1 = _mm_loadu_si128((__m128i *)(ref16_l1 + 8));

            out8_0_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel0, 2), _mm_srli_epi16(in_pixel1, 2));

            //AVG
            avg8_0_u8 = _mm_avg_epu8(out8_0_u8_l0, out8_0_u8_l1);

            _mm_storeu_si128((__m128i *)dst_ptr, avg8_0_u8);

            //--------
            //Line Two
            //--------

            //List0

            in_pixel4 = _mm_loadu_si128((__m128i *)(ref16_l0 + ref_l0_stride));
            in_pixel5 = _mm_loadu_si128((__m128i *)(ref16_l0 + ref_l0_stride + 8));

            out8_2_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel4, 2), _mm_srli_epi16(in_pixel5, 2));

            //List1

            in_pixel4 = _mm_loadu_si128((__m128i *)(ref16_l1 + ref_l1_stride));
            in_pixel5 = _mm_loadu_si128((__m128i *)(ref16_l1 + ref_l1_stride + 8));

            out8_2_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel4, 2), _mm_srli_epi16(in_pixel5, 2));

            //AVG
            avg8_2_u8 = _mm_avg_epu8(out8_2_u8_l0, out8_2_u8_l1);

            _mm_storeu_si128((__m128i *)(dst_ptr + dst_stride), avg8_2_u8);

            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;
        }
    } else if (width == 32) {
        __m128i in_pixel2, in_pixel3, in_pixel4, in_pixel5, in_pixel6, in_pixel7;
        __m128i out8_0_u8_l0, out8_1_u8_l0, out8_2_u8_l0, out8_3_u8_l0;
        __m128i out8_0_u8_l1, out8_1_u8_l1, out8_2_u8_l1, out8_3_u8_l1;
        __m128i avg8_0_u8, avg8_1_u8, avg8_2_u8, avg8_3_u8;

        for (y = 0; y < height; y += 2) {
            //--------
            //Line One
            //--------

            //List0

            in_pixel0 = _mm_loadu_si128((__m128i *)ref16_l0);
            in_pixel1 = _mm_loadu_si128((__m128i *)(ref16_l0 + 8));
            in_pixel2 = _mm_loadu_si128((__m128i *)(ref16_l0 + 16));
            in_pixel3 = _mm_loadu_si128((__m128i *)(ref16_l0 + 24));

            out8_0_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel0, 2), _mm_srli_epi16(in_pixel1, 2));
            out8_1_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel2, 2), _mm_srli_epi16(in_pixel3, 2));

            //List1

            in_pixel0 = _mm_loadu_si128((__m128i *)ref16_l1);
            in_pixel1 = _mm_loadu_si128((__m128i *)(ref16_l1 + 8));
            in_pixel2 = _mm_loadu_si128((__m128i *)(ref16_l1 + 16));
            in_pixel3 = _mm_loadu_si128((__m128i *)(ref16_l1 + 24));

            out8_0_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel0, 2), _mm_srli_epi16(in_pixel1, 2));
            out8_1_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel2, 2), _mm_srli_epi16(in_pixel3, 2));

            //AVG
            avg8_0_u8 = _mm_avg_epu8(out8_0_u8_l0, out8_0_u8_l1);
            avg8_1_u8 = _mm_avg_epu8(out8_1_u8_l0, out8_1_u8_l1);

            _mm_storeu_si128((__m128i *)dst_ptr, avg8_0_u8);
            _mm_storeu_si128((__m128i *)(dst_ptr + 16), avg8_1_u8);

            //--------
            //Line Two
            //--------

            //List0

            in_pixel4 = _mm_loadu_si128((__m128i *)(ref16_l0 + ref_l0_stride));
            in_pixel5 = _mm_loadu_si128((__m128i *)(ref16_l0 + ref_l0_stride + 8));
            in_pixel6 = _mm_loadu_si128((__m128i *)(ref16_l0 + ref_l0_stride + 16));
            in_pixel7 = _mm_loadu_si128((__m128i *)(ref16_l0 + ref_l0_stride + 24));

            out8_2_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel4, 2), _mm_srli_epi16(in_pixel5, 2));
            out8_3_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel6, 2), _mm_srli_epi16(in_pixel7, 2));

            //List1

            in_pixel4 = _mm_loadu_si128((__m128i *)(ref16_l1 + ref_l1_stride));
            in_pixel5 = _mm_loadu_si128((__m128i *)(ref16_l1 + ref_l1_stride + 8));
            in_pixel6 = _mm_loadu_si128((__m128i *)(ref16_l1 + ref_l1_stride + 16));
            in_pixel7 = _mm_loadu_si128((__m128i *)(ref16_l1 + ref_l1_stride + 24));

            out8_2_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel4, 2), _mm_srli_epi16(in_pixel5, 2));
            out8_3_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel6, 2), _mm_srli_epi16(in_pixel7, 2));

            //AVG
            avg8_2_u8 = _mm_avg_epu8(out8_2_u8_l0, out8_2_u8_l1);
            avg8_3_u8 = _mm_avg_epu8(out8_3_u8_l0, out8_3_u8_l1);

            _mm_storeu_si128((__m128i *)(dst_ptr + dst_stride), avg8_2_u8);
            _mm_storeu_si128((__m128i *)(dst_ptr + dst_stride + 16), avg8_3_u8);

            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;
        }
    } else if (width == 64) {
        __m128i in_pixel2, in_pixel3, in_pixel4, in_pixel5, in_pixel6, in_pixel7;
        __m128i out8_0_u8_l0, out8_1_u8_l0, out8_2_u8_l0, out8_3_u8_l0;
        __m128i out8_0_u8_l1, out8_1_u8_l1, out8_2_u8_l1, out8_3_u8_l1;
        __m128i avg8_0_u8, avg8_1_u8, avg8_2_u8, avg8_3_u8;

        for (y = 0; y < height; ++y) {
            //List0

            in_pixel0 = _mm_loadu_si128((__m128i *)ref16_l0);
            in_pixel1 = _mm_loadu_si128((__m128i *)(ref16_l0 + 8));
            in_pixel2 = _mm_loadu_si128((__m128i *)(ref16_l0 + 16));
            in_pixel3 = _mm_loadu_si128((__m128i *)(ref16_l0 + 24));
            in_pixel4 = _mm_loadu_si128((__m128i *)(ref16_l0 + 32));
            in_pixel5 = _mm_loadu_si128((__m128i *)(ref16_l0 + 40));
            in_pixel6 = _mm_loadu_si128((__m128i *)(ref16_l0 + 48));
            in_pixel7 = _mm_loadu_si128((__m128i *)(ref16_l0 + 56));

            out8_0_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel0, 2), _mm_srli_epi16(in_pixel1, 2));
            out8_1_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel2, 2), _mm_srli_epi16(in_pixel3, 2));
            out8_2_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel4, 2), _mm_srli_epi16(in_pixel5, 2));
            out8_3_u8_l0 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel6, 2), _mm_srli_epi16(in_pixel7, 2));

            //List1

            in_pixel0 = _mm_loadu_si128((__m128i *)ref16_l1);
            in_pixel1 = _mm_loadu_si128((__m128i *)(ref16_l1 + 8));
            in_pixel2 = _mm_loadu_si128((__m128i *)(ref16_l1 + 16));
            in_pixel3 = _mm_loadu_si128((__m128i *)(ref16_l1 + 24));
            in_pixel4 = _mm_loadu_si128((__m128i *)(ref16_l1 + 32));
            in_pixel5 = _mm_loadu_si128((__m128i *)(ref16_l1 + 40));
            in_pixel6 = _mm_loadu_si128((__m128i *)(ref16_l1 + 48));
            in_pixel7 = _mm_loadu_si128((__m128i *)(ref16_l1 + 56));

            //Note: old Version used to use _mm_and_si128 to mask the MSB bits of the pixels
            out8_0_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel0, 2), _mm_srli_epi16(in_pixel1, 2));
            out8_1_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel2, 2), _mm_srli_epi16(in_pixel3, 2));
            out8_2_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel4, 2), _mm_srli_epi16(in_pixel5, 2));
            out8_3_u8_l1 =
                _mm_packus_epi16(_mm_srli_epi16(in_pixel6, 2), _mm_srli_epi16(in_pixel7, 2));

            //AVG
            avg8_0_u8 = _mm_avg_epu8(out8_0_u8_l0, out8_0_u8_l1);
            avg8_1_u8 = _mm_avg_epu8(out8_1_u8_l0, out8_1_u8_l1);
            avg8_2_u8 = _mm_avg_epu8(out8_2_u8_l0, out8_2_u8_l1);
            avg8_3_u8 = _mm_avg_epu8(out8_3_u8_l0, out8_3_u8_l1);

            _mm_storeu_si128((__m128i *)dst_ptr, avg8_0_u8);
            _mm_storeu_si128((__m128i *)(dst_ptr + 16), avg8_1_u8);
            _mm_storeu_si128((__m128i *)(dst_ptr + 32), avg8_2_u8);
            _mm_storeu_si128((__m128i *)(dst_ptr + 48), avg8_3_u8);

            dst_ptr += dst_stride;
            ref16_l0 += ref_l0_stride;
            ref16_l1 += ref_l1_stride;
        }
    }

    return;
}
/********************************************************************************************************************
eb_enc_msb_pack2d_sse2_intrin
*********************************************************************************************************************/
void eb_enc_msb_pack2d_sse2_intrin(uint8_t *in8_bit_buffer, uint32_t in8_stride,
                                   uint8_t *inn_bit_buffer, uint16_t *out16_bit_buffer,
                                   uint32_t inn_stride, uint32_t out_stride, uint32_t width,
                                   uint32_t height) {
    uint32_t count_width, count_height;

    if (width == 4) {
        for (count_height = 0; count_height < height; count_height += 2) {
            _mm_storel_epi64(
                (__m128i *)(out16_bit_buffer),
                _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(inn_bit_buffer)),
                                                 _mm_cvtsi32_si128(*(uint32_t *)(in8_bit_buffer))),
                               6));
            _mm_storel_epi64(
                (__m128i *)(out16_bit_buffer + out_stride),
                _mm_srli_epi16(_mm_unpacklo_epi8(
                                   _mm_cvtsi32_si128(*(uint32_t *)(inn_bit_buffer + inn_stride)),
                                   _mm_cvtsi32_si128(*(uint32_t *)(in8_bit_buffer + in8_stride))),
                               6));
            out16_bit_buffer += (out_stride << 1);
            in8_bit_buffer += (in8_stride << 1);
            inn_bit_buffer += (inn_stride << 1);
        }
    } else if (width == 8) {
        for (count_height = 0; count_height < height; count_height += 2) {
            _mm_storeu_si128(
                (__m128i *)(out16_bit_buffer),
                _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(inn_bit_buffer)),
                                                 _mm_loadl_epi64((__m128i *)(in8_bit_buffer))),
                               6));
            _mm_storeu_si128(
                (__m128i *)(out16_bit_buffer + out_stride),
                _mm_srli_epi16(
                    _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(inn_bit_buffer + inn_stride)),
                                      _mm_loadl_epi64((__m128i *)(in8_bit_buffer + in8_stride))),
                    6));
            out16_bit_buffer += (out_stride << 1);
            in8_bit_buffer += (in8_stride << 1);
            inn_bit_buffer += (inn_stride << 1);
        }
    } else if (width == 16) {
        __m128i out_pixel_1, out_pixel_2, out_pixel_3, out_pixel_4, inn_bit_buffer_lo,
            inn_bit_buffer_hi, in_8bit_buffer_lo, in_8bit_buffer_hi;

        for (count_height = 0; count_height < height; count_height += 2) {
            inn_bit_buffer_lo = _mm_loadu_si128((__m128i *)inn_bit_buffer);
            inn_bit_buffer_hi = _mm_loadu_si128((__m128i *)(inn_bit_buffer + inn_stride));
            in_8bit_buffer_lo = _mm_loadu_si128((__m128i *)in8_bit_buffer);
            in_8bit_buffer_hi = _mm_loadu_si128((__m128i *)(in8_bit_buffer + in8_stride));

            out_pixel_1 =
                _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_lo, in_8bit_buffer_lo), 6);
            out_pixel_2 =
                _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_lo, in_8bit_buffer_lo), 6);
            out_pixel_3 =
                _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_hi, in_8bit_buffer_hi), 6);
            out_pixel_4 =
                _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_hi, in_8bit_buffer_hi), 6);

            _mm_storeu_si128((__m128i *)out16_bit_buffer, out_pixel_1);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 8), out_pixel_2);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + out_stride), out_pixel_3);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + out_stride + 8), out_pixel_4);

            in8_bit_buffer += (in8_stride << 1);
            inn_bit_buffer += (inn_stride << 1);
            out16_bit_buffer += (out_stride << 1);
        }
    } else if (width == 32) {
        __m128i inn_bit_buffer_1, inn_bit_buffer_2, inn_bit_buffer_3, inn_bit_buffer_4,
            in_8bit_buffer1, in_8bit_buffer2, in_8bit_buffer3, in_8bit_buffer4;
        __m128i out_pixel_1, out_pixel_2, out_pixel_3, out_pixel_4, out_pixel_5, out_pixel_6,
            out_pixel_7, out_pixel_8;

        for (count_height = 0; count_height < height; count_height += 2) {
            inn_bit_buffer_1 = _mm_loadu_si128((__m128i *)inn_bit_buffer);
            inn_bit_buffer_2 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + 16));
            inn_bit_buffer_3 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + inn_stride));
            inn_bit_buffer_4 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + inn_stride + 16));

            in_8bit_buffer1 = _mm_loadu_si128((__m128i *)in8_bit_buffer);
            in_8bit_buffer2 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + 16));
            in_8bit_buffer3 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + in8_stride));
            in_8bit_buffer4 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + in8_stride + 16));

            out_pixel_1 = _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_1, in_8bit_buffer1), 6);
            out_pixel_2 = _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_1, in_8bit_buffer1), 6);
            out_pixel_3 = _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_2, in_8bit_buffer2), 6);
            out_pixel_4 = _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_2, in_8bit_buffer2), 6);
            out_pixel_5 = _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_3, in_8bit_buffer3), 6);
            out_pixel_6 = _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_3, in_8bit_buffer3), 6);
            out_pixel_7 = _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_4, in_8bit_buffer4), 6);
            out_pixel_8 = _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_4, in_8bit_buffer4), 6);

            _mm_storeu_si128((__m128i *)out16_bit_buffer, out_pixel_1);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 8), out_pixel_2);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 16), out_pixel_3);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 24), out_pixel_4);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + out_stride), out_pixel_5);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + out_stride + 8), out_pixel_6);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + out_stride + 16), out_pixel_7);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + out_stride + 24), out_pixel_8);

            in8_bit_buffer += (in8_stride << 1);
            inn_bit_buffer += (inn_stride << 1);
            out16_bit_buffer += (out_stride << 1);
        }
    } else if (width == 64) {
        __m128i inn_bit_buffer_1, inn_bit_buffer_2, inn_bit_buffer_3, inn_bit_buffer_4,
            in_8bit_buffer1, in_8bit_buffer2, in_8bit_buffer3, in_8bit_buffer4;
        __m128i out_pixel_1, out_pixel_2, out_pixel_3, out_pixel_4, out_pixel_5, out_pixel_6,
            out_pixel_7, out_pixel_8;

        for (count_height = 0; count_height < height; ++count_height) {
            inn_bit_buffer_1 = _mm_loadu_si128((__m128i *)inn_bit_buffer);
            inn_bit_buffer_2 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + 16));
            inn_bit_buffer_3 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + 32));
            inn_bit_buffer_4 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + 48));

            in_8bit_buffer1 = _mm_loadu_si128((__m128i *)in8_bit_buffer);
            in_8bit_buffer2 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + 16));
            in_8bit_buffer3 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + 32));
            in_8bit_buffer4 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + 48));

            out_pixel_1 = _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_1, in_8bit_buffer1), 6);
            out_pixel_2 = _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_1, in_8bit_buffer1), 6);
            out_pixel_3 = _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_2, in_8bit_buffer2), 6);
            out_pixel_4 = _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_2, in_8bit_buffer2), 6);
            out_pixel_5 = _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_3, in_8bit_buffer3), 6);
            out_pixel_6 = _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_3, in_8bit_buffer3), 6);
            out_pixel_7 = _mm_srli_epi16(_mm_unpacklo_epi8(inn_bit_buffer_4, in_8bit_buffer4), 6);
            out_pixel_8 = _mm_srli_epi16(_mm_unpackhi_epi8(inn_bit_buffer_4, in_8bit_buffer4), 6);

            _mm_storeu_si128((__m128i *)out16_bit_buffer, out_pixel_1);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 8), out_pixel_2);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 16), out_pixel_3);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 24), out_pixel_4);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 32), out_pixel_5);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 40), out_pixel_6);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 48), out_pixel_7);
            _mm_storeu_si128((__m128i *)(out16_bit_buffer + 56), out_pixel_8);

            in8_bit_buffer += in8_stride;
            inn_bit_buffer += inn_stride;
            out16_bit_buffer += out_stride;
        }
    } else {
        uint32_t in_n_stride_diff = (inn_stride << 1) - width;
        uint32_t in_8_stride_diff = (in8_stride << 1) - width;
        uint32_t out_stride_diff  = (out_stride << 1) - width;

        if (!(width & 7)) {
            for (count_height = 0; count_height < height; count_height += 2) {
                for (count_width = 0; count_width < width; count_width += 8) {
                    _mm_storeu_si128(
                        (__m128i *)(out16_bit_buffer),
                        _mm_srli_epi16(
                            _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(inn_bit_buffer)),
                                              _mm_loadl_epi64((__m128i *)(in8_bit_buffer))),
                            6));
                    _mm_storeu_si128(
                        (__m128i *)(out16_bit_buffer + out_stride),
                        _mm_srli_epi16(
                            _mm_unpacklo_epi8(
                                _mm_loadl_epi64((__m128i *)(inn_bit_buffer + inn_stride)),
                                _mm_loadl_epi64((__m128i *)(in8_bit_buffer + in8_stride))),
                            6));
                    out16_bit_buffer += 8;
                    in8_bit_buffer += 8;
                    inn_bit_buffer += 8;
                }
                in8_bit_buffer += in_8_stride_diff;
                inn_bit_buffer += in_n_stride_diff;
                out16_bit_buffer += out_stride_diff;
            }
        } else {
            for (count_height = 0; count_height < height; count_height += 2) {
                for (count_width = 0; count_width < width; count_width += 4) {
                    _mm_storel_epi64(
                        (__m128i *)(out16_bit_buffer),
                        _mm_srli_epi16(
                            _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(inn_bit_buffer)),
                                              _mm_cvtsi32_si128(*(uint32_t *)(in8_bit_buffer))),
                            6));
                    _mm_storel_epi64(
                        (__m128i *)(out16_bit_buffer + out_stride),
                        _mm_srli_epi16(
                            _mm_unpacklo_epi8(
                                _mm_cvtsi32_si128(*(uint32_t *)(inn_bit_buffer + inn_stride)),
                                _mm_cvtsi32_si128(*(uint32_t *)(in8_bit_buffer + in8_stride))),
                            6));
                    out16_bit_buffer += 4;
                    in8_bit_buffer += 4;
                    inn_bit_buffer += 4;
                }
                in8_bit_buffer += in_8_stride_diff;
                inn_bit_buffer += in_n_stride_diff;
                out16_bit_buffer += out_stride_diff;
            }
        }
    }
}
