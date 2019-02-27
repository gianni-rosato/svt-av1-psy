/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPackUnPack_SSE2.h"


#include <emmintrin.h>

#include "EbDefinitions.h"

/****************************************************************************************
eb_enc_msb_un_pack2d_sse2_intrin
******************************************************************************************/

void eb_enc_msb_un_pack2d_sse2_intrin(
    uint16_t      *in16_bit_buffer,
    uint32_t       in_stride,
    uint8_t       *out8_bit_buffer,
    uint8_t       *outn_bit_buffer,
    uint32_t       out8_stride,
    uint32_t       outn_stride,
    uint32_t       width,
    uint32_t       height)
{
    uint32_t x, y;

    __m128i xmm_3, xmm_00FF, inPixel0, inPixel1, tempPixel0, tempPixel1, inPixel1_shftR_2_U8, inPixel0_shftR_2_U8, inPixel0_shftR_2, inPixel1_shftR_2;
    __m128i tempPixel0_U8, tempPixel1_U8;

    xmm_3 = _mm_set1_epi16(0x0003);
    xmm_00FF = _mm_set1_epi16(0x00FF);

    if (width == 4)
    {
        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadl_epi64((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadl_epi64((__m128i*)(in16_bit_buffer + in_stride));

            tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
            tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);

            tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
            tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

            *(uint32_t*)outn_bit_buffer = _mm_cvtsi128_si32(tempPixel0_U8);
            *(uint32_t*)(outn_bit_buffer + outn_stride) = _mm_cvtsi128_si32(tempPixel1_U8);
            *(uint32_t*)out8_bit_buffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
            *(uint32_t*)(out8_bit_buffer + out8_stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

            outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    }
    else if (width == 8)
    {
        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));

            tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
            tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);

            tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
            tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

            _mm_storel_epi64((__m128i*)outn_bit_buffer, tempPixel0_U8);
            _mm_storel_epi64((__m128i*)(outn_bit_buffer + outn_stride), tempPixel1_U8);
            _mm_storel_epi64((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
            _mm_storel_epi64((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

            outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    }
    else if (width == 16)
    {
        __m128i inPixel2, inPixel3;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));
            inPixel3 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 8));

            tempPixel0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            tempPixel1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));

            inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));

            _mm_storeu_si128((__m128i*)outn_bit_buffer, tempPixel0_U8);
            _mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride), tempPixel1_U8);
            _mm_storeu_si128((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

            outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    }
    else if (width == 32)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i outn0_U8, outn1_U8, outn2_U8, outn3_U8, out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 24));

            outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
            outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
            outn3_U8 = _mm_packus_epi16(_mm_and_si128(inPixel6, xmm_3), _mm_and_si128(inPixel7, xmm_3));

            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

            _mm_storeu_si128((__m128i*)outn_bit_buffer, outn0_U8);
            _mm_storeu_si128((__m128i*)(outn_bit_buffer + 16), outn1_U8);
            _mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride), outn2_U8);
            _mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride + 16), outn3_U8);

            _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride + 16), out8_3_U8);

            outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    }
    else if (width == 64)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i outn0_U8, outn1_U8, outn2_U8, outn3_U8, out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; ++y)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 56));

            outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
            outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
            outn3_U8 = _mm_packus_epi16(_mm_and_si128(inPixel6, xmm_3), _mm_and_si128(inPixel7, xmm_3));

            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

            _mm_storeu_si128((__m128i*)outn_bit_buffer, outn0_U8);
            _mm_storeu_si128((__m128i*)(outn_bit_buffer + 16), outn1_U8);
            _mm_storeu_si128((__m128i*)(outn_bit_buffer + 32), outn2_U8);
            _mm_storeu_si128((__m128i*)(outn_bit_buffer + 48), outn3_U8);

            _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 32), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 48), out8_3_U8);

            outn_bit_buffer += outn_stride;
            out8_bit_buffer += out8_stride;
            in16_bit_buffer += in_stride;
        }
    }
    else
    {
        uint32_t inStrideDiff = (2 * in_stride) - width;
        uint32_t out8StrideDiff = (2 * out8_stride) - width;
        uint32_t outnStrideDiff = (2 * outn_stride) - width;

        uint32_t inStrideDiff64 = in_stride - width;
        uint32_t out8StrideDiff64 = out8_stride - width;
        uint32_t outnStrideDiff64 = outn_stride - width;

        if (!(width & 63))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i outn0_U8, outn1_U8, outn2_U8, outn3_U8, out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 1) {
                for (y = 0; y < width; y += 64) {

                    inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 32));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 40));
                    inPixel6 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 48));
                    inPixel7 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 56));

                    outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
                    outn3_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel6, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel7, xmm_3), 6));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    _mm_storeu_si128((__m128i*)outn_bit_buffer, outn0_U8);
                    _mm_storeu_si128((__m128i*)(outn_bit_buffer + 16), outn1_U8);
                    _mm_storeu_si128((__m128i*)(outn_bit_buffer + 32), outn2_U8);
                    _mm_storeu_si128((__m128i*)(outn_bit_buffer + 48), outn3_U8);

                    _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 32), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 48), out8_3_U8);

                    outn_bit_buffer += 64;
                    out8_bit_buffer += 64;
                    in16_bit_buffer += 64;
                }
                in16_bit_buffer += inStrideDiff64;
                outn_bit_buffer += outnStrideDiff64;
                out8_bit_buffer += out8StrideDiff64;
            }
        }
        else if (!(width & 31))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i outn0_U8, outn1_U8, outn2_U8, outn3_U8, out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 32)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride + 8));
                    inPixel6 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 16));
                    inPixel7 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 24));

                    outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
                    outn3_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel6, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel7, xmm_3), 6));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    _mm_storeu_si128((__m128i*)outn_bit_buffer, outn0_U8);
                    _mm_storeu_si128((__m128i*)(outn_bit_buffer + 16), outn1_U8);
                    _mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride), outn2_U8);
                    _mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride + 16), outn3_U8);

                    _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride + 16), out8_3_U8);

                    outn_bit_buffer += 32;
                    out8_bit_buffer += 32;
                    in16_bit_buffer += 32;
                }
                in16_bit_buffer += inStrideDiff;
                outn_bit_buffer += outnStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else if (!(width & 15))
        {
            __m128i inPixel2, inPixel3;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 16)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));
                    inPixel3 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 8));

                    tempPixel0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    tempPixel1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));

                    _mm_storeu_si128((__m128i*)outn_bit_buffer, tempPixel0_U8);
                    _mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride), tempPixel1_U8);
                    _mm_storeu_si128((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

                    outn_bit_buffer += 16;
                    out8_bit_buffer += 16;
                    in16_bit_buffer += 16;
                }
                in16_bit_buffer += inStrideDiff;
                outn_bit_buffer += outnStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else if (!(width & 7))
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 8)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));

                    tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
                    tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);

                    tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
                    tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    _mm_storel_epi64((__m128i*)outn_bit_buffer, tempPixel0_U8);
                    _mm_storel_epi64((__m128i*)(outn_bit_buffer + outn_stride), tempPixel1_U8);
                    _mm_storel_epi64((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
                    _mm_storel_epi64((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

                    outn_bit_buffer += 8;
                    out8_bit_buffer += 8;
                    in16_bit_buffer += 8;
                }
                in16_bit_buffer += inStrideDiff;
                outn_bit_buffer += outnStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 4)
                {
                    inPixel0 = _mm_loadl_epi64((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadl_epi64((__m128i*)(in16_bit_buffer + in_stride));

                    tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
                    tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);

                    tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
                    tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    *(uint32_t*)outn_bit_buffer = _mm_cvtsi128_si32(tempPixel0_U8);
                    *(uint32_t*)(outn_bit_buffer + outn_stride) = _mm_cvtsi128_si32(tempPixel1_U8);
                    *(uint32_t*)out8_bit_buffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
                    *(uint32_t*)(out8_bit_buffer + out8_stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

                    outn_bit_buffer += 4;
                    out8_bit_buffer += 4;
                    in16_bit_buffer += 4;
                }
                in16_bit_buffer += inStrideDiff;
                outn_bit_buffer += outnStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
    }
    return;
}

void eb_enc_un_pack8_bit_data_safe_sub_sse2_intrin(
    uint16_t      *in16_bit_buffer,
    uint32_t       in_stride,
    uint8_t       *out8_bit_buffer,
    uint32_t       out8_stride,
    uint32_t       width,
    uint32_t       height,
    EbBool      sub_pred
)
{

    uint32_t x, y;

    __m128i xmm_00FF, inPixel0, inPixel1, inPixel1_shftR_2_U8, inPixel0_shftR_2_U8, inPixel0_shftR_2, inPixel1_shftR_2;


    xmm_00FF = _mm_set1_epi16(0x00FF);

    if (width == 8)
    {
        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));


            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);


            _mm_storel_epi64((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
            _mm_storel_epi64((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }

        if (sub_pred) {
            in16_bit_buffer -= (in_stride >> 1);
            out8_bit_buffer -= (out8_stride >> 1);
            inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            _mm_storel_epi64((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
        }



    }
    else if (width == 16)
    {
        __m128i inPixel2, inPixel3;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));
            inPixel3 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 8));


            inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));


            _mm_storeu_si128((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);


            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }

        if (sub_pred) {
            in16_bit_buffer -= (in_stride >> 1);
            out8_bit_buffer -= (out8_stride >> 1);
            inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + 8));
            inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            _mm_storeu_si128((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
        }

    }
    else if (width == 32)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i  out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 24));


            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));


            _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride + 16), out8_3_U8);


            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }

        if (sub_pred) {
            in16_bit_buffer -= (in_stride >> 1);
            out8_bit_buffer -= (out8_stride >> 1);
            inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);

        }

    }
    else if (width == 64)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i  out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; ++y)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 56));


            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));


            _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 32), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 48), out8_3_U8);


            out8_bit_buffer += out8_stride;
            in16_bit_buffer += in_stride;
        }

        if (sub_pred) {
            in16_bit_buffer -= (in_stride >> 1);
            out8_bit_buffer -= (out8_stride >> 1);
            inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 56));


            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));


            _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 32), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 48), out8_3_U8);
        }

    }
    else
    {
        uint32_t inStrideDiff = (2 * in_stride) - width;
        uint32_t out8StrideDiff = (2 * out8_stride) - width;

        uint32_t inStrideDiff64 = in_stride - width;
        uint32_t out8StrideDiff64 = out8_stride - width;

        if (!(width & 63))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 1) {
                for (y = 0; y < width; y += 64) {

                    inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 32));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 40));
                    inPixel6 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 48));
                    inPixel7 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 56));


                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));


                    _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 32), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 48), out8_3_U8);

                    out8_bit_buffer += 64;
                    in16_bit_buffer += 64;
                }
                in16_bit_buffer += inStrideDiff64;
                out8_bit_buffer += out8StrideDiff64;
            }
        }
        else if (!(width & 31))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 32)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride + 8));
                    inPixel6 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 16));
                    inPixel7 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 24));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride + 16), out8_3_U8);

                    out8_bit_buffer += 32;
                    in16_bit_buffer += 32;
                }
                in16_bit_buffer += inStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else if (!(width & 15))
        {
            __m128i inPixel2, inPixel3;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 16)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));
                    inPixel3 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 8));

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));


                    _mm_storeu_si128((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

                    out8_bit_buffer += 16;
                    in16_bit_buffer += 16;
                }
                in16_bit_buffer += inStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else if (!(width & 7))
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 8)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    _mm_storel_epi64((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
                    _mm_storel_epi64((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

                    out8_bit_buffer += 8;
                    in16_bit_buffer += 8;
                }
                in16_bit_buffer += inStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 4)
                {
                    inPixel0 = _mm_loadl_epi64((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadl_epi64((__m128i*)(in16_bit_buffer + in_stride));


                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    *(uint32_t*)out8_bit_buffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
                    *(uint32_t*)(out8_bit_buffer + out8_stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

                    out8_bit_buffer += 4;
                    in16_bit_buffer += 4;
                }
                in16_bit_buffer += inStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }

    }


    return;
}

void eb_enc_un_pack8_bit_data_sse2_intrin(
    uint16_t      *in16_bit_buffer,
    uint32_t       in_stride,
    uint8_t       *out8_bit_buffer,
    uint32_t       out8_stride,
    uint32_t       width,
    uint32_t       height)
{
    uint32_t x, y;

    __m128i xmm_00FF, inPixel0, inPixel1, inPixel1_shftR_2_U8, inPixel0_shftR_2_U8, inPixel0_shftR_2, inPixel1_shftR_2;
    //    __m128i tempPixel0_U8, tempPixel1_U8;

    xmm_00FF = _mm_set1_epi16(0x00FF);

    if (width == 4)
    {
        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadl_epi64((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadl_epi64((__m128i*)(in16_bit_buffer + in_stride));

            //tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
            //tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);
            //
            //tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
            //tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

            //*(uint32_t*)outn_bit_buffer = _mm_cvtsi128_si32(tempPixel0_U8);
            //*(uint32_t*)(outn_bit_buffer + outn_stride) = _mm_cvtsi128_si32(tempPixel1_U8);
            *(uint32_t*)out8_bit_buffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
            *(uint32_t*)(out8_bit_buffer + out8_stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

            //outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    }
    else if (width == 8)
    {
        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));

            //tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
            //tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);
            //
            //tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
            //tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

            //_mm_storel_epi64((__m128i*)outn_bit_buffer, tempPixel0_U8);
            //_mm_storel_epi64((__m128i*)(outn_bit_buffer + outn_stride), tempPixel1_U8);
            _mm_storel_epi64((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
            _mm_storel_epi64((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

            //outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    }
    else if (width == 16)
    {
        __m128i inPixel2, inPixel3;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));
            inPixel3 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 8));

            //tempPixel0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            //tempPixel1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));

            inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));

            //_mm_storeu_si128((__m128i*)outn_bit_buffer, tempPixel0_U8);
            //_mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride), tempPixel1_U8);
            _mm_storeu_si128((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

            //outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    }
    else if (width == 32)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i /*outn0_U8, outn1_U8, outn2_U8, outn3_U8,*/ out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 24));

            //outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            //outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
            //outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
            //outn3_U8 = _mm_packus_epi16(_mm_and_si128(inPixel6, xmm_3), _mm_and_si128(inPixel7, xmm_3));

            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

            //_mm_storeu_si128((__m128i*)outn_bit_buffer, outn0_U8);
            //_mm_storeu_si128((__m128i*)(outn_bit_buffer + 16), outn1_U8);
            //_mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride), outn2_U8);
            //_mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride + 16), outn3_U8);

            _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride + 16), out8_3_U8);

            //outn_bit_buffer += 2 * outn_stride;
            out8_bit_buffer += 2 * out8_stride;
            in16_bit_buffer += 2 * in_stride;
        }
    }
    else if (width == 64)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i /*outn0_U8, outn1_U8, outn2_U8, outn3_U8,*/ out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; ++y)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 56));

            //outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            //outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
            //outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
            //outn3_U8 = _mm_packus_epi16(_mm_and_si128(inPixel6, xmm_3), _mm_and_si128(inPixel7, xmm_3));

            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

            //_mm_storeu_si128((__m128i*)outn_bit_buffer, outn0_U8);
            //_mm_storeu_si128((__m128i*)(outn_bit_buffer + 16), outn1_U8);
            //_mm_storeu_si128((__m128i*)(outn_bit_buffer + 32), outn2_U8);
            //_mm_storeu_si128((__m128i*)(outn_bit_buffer + 48), outn3_U8);

            _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 32), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8_bit_buffer + 48), out8_3_U8);

            //outn_bit_buffer += outn_stride;
            out8_bit_buffer += out8_stride;
            in16_bit_buffer += in_stride;
        }
    }
    else
    {
        uint32_t inStrideDiff = (2 * in_stride) - width;
        uint32_t out8StrideDiff = (2 * out8_stride) - width;
        //uint32_t outnStrideDiff = (2 * outn_stride) - width;

        uint32_t inStrideDiff64 = in_stride - width;
        uint32_t out8StrideDiff64 = out8_stride - width;
        //uint32_t outnStrideDiff64 = outn_stride - width;

        if (!(width & 63))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i /*outn0_U8, outn1_U8, outn2_U8, outn3_U8, */out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 1) {
                for (y = 0; y < width; y += 64) {

                    inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 32));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 40));
                    inPixel6 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 48));
                    inPixel7 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 56));

                    //outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    //outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    //outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
                    //outn3_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel6, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel7, xmm_3), 6));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    //_mm_storeu_si128((__m128i*)outn_bit_buffer, outn0_U8);
                    //_mm_storeu_si128((__m128i*)(outn_bit_buffer + 16), outn1_U8);
                    //_mm_storeu_si128((__m128i*)(outn_bit_buffer + 32), outn2_U8);
                    //_mm_storeu_si128((__m128i*)(outn_bit_buffer + 48), outn3_U8);

                    _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 32), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 48), out8_3_U8);

                    //outn_bit_buffer += 64;
                    out8_bit_buffer += 64;
                    in16_bit_buffer += 64;
                }
                in16_bit_buffer += inStrideDiff64;
                //outn_bit_buffer += outnStrideDiff64;
                out8_bit_buffer += out8StrideDiff64;
            }
        }
        else if (!(width & 31))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i /*outn0_U8, outn1_U8, outn2_U8, outn3_U8,*/ out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 32)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16_bit_buffer + in_stride + 8));
                    inPixel6 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 16));
                    inPixel7 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 24));

                    //outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    //outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    //outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
                    //outn3_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel6, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel7, xmm_3), 6));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    //_mm_storeu_si128((__m128i*)outn_bit_buffer, outn0_U8);
                    //_mm_storeu_si128((__m128i*)(outn_bit_buffer + 16), outn1_U8);
                    //_mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride), outn2_U8);
                    //_mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride + 16), outn3_U8);

                    _mm_storeu_si128((__m128i*)out8_bit_buffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride + 16), out8_3_U8);

                    //outn_bit_buffer += 32;
                    out8_bit_buffer += 32;
                    in16_bit_buffer += 32;
                }
                in16_bit_buffer += inStrideDiff;
                //outn_bit_buffer += outnStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else if (!(width & 15))
        {
            __m128i inPixel2, inPixel3;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 16)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));
                    inPixel3 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride + 8));

                    //tempPixel0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    //tempPixel1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    //
                    inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));

                    //_mm_storeu_si128((__m128i*)outn_bit_buffer, tempPixel0_U8);
                    //_mm_storeu_si128((__m128i*)(outn_bit_buffer + outn_stride), tempPixel1_U8);
                    _mm_storeu_si128((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
                    _mm_storeu_si128((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

                    //outn_bit_buffer += 16;
                    out8_bit_buffer += 16;
                    in16_bit_buffer += 16;
                }
                in16_bit_buffer += inStrideDiff;
                //outn_bit_buffer += outnStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else if (!(width & 7))
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 8)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16_bit_buffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16_bit_buffer + in_stride));

                    //tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
                    //tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    //_mm_storel_epi64((__m128i*)outn_bit_buffer, tempPixel0_U8);
                    //_mm_storel_epi64((__m128i*)(outn_bit_buffer + outn_stride), tempPixel1_U8);
                    _mm_storel_epi64((__m128i*)out8_bit_buffer, inPixel0_shftR_2_U8);
                    _mm_storel_epi64((__m128i*)(out8_bit_buffer + out8_stride), inPixel1_shftR_2_U8);

                    //outn_bit_buffer += 8;
                    out8_bit_buffer += 8;
                    in16_bit_buffer += 8;
                }
                in16_bit_buffer += inStrideDiff;
                //outn_bit_buffer += outnStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
        else
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 4)
                {
                    inPixel0 = _mm_loadl_epi64((__m128i*)in16_bit_buffer);
                    inPixel1 = _mm_loadl_epi64((__m128i*)(in16_bit_buffer + in_stride));

                    //tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
                    //tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    //*(uint32_t*)outn_bit_buffer = _mm_cvtsi128_si32(tempPixel0_U8);
                    //*(uint32_t*)(outn_bit_buffer + outn_stride) = _mm_cvtsi128_si32(tempPixel1_U8);
                    *(uint32_t*)out8_bit_buffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
                    *(uint32_t*)(out8_bit_buffer + out8_stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

                    //outn_bit_buffer += 4;
                    out8_bit_buffer += 4;
                    in16_bit_buffer += 4;
                }
                in16_bit_buffer += inStrideDiff;
                //outn_bit_buffer += outnStrideDiff;
                out8_bit_buffer += out8StrideDiff;
            }
        }
    }

    return;
}
void unpack_avg_sse2_intrin(
    uint16_t *ref16_l0,
    uint32_t  ref_l0_stride,
    uint16_t *ref16_l1,
    uint32_t  ref_l1_stride,
    uint8_t  *dst_ptr,
    uint32_t  dst_stride,
    uint32_t  width,
    uint32_t  height)
{

    uint32_t  y;
    __m128i inPixel0, inPixel1;



    if (width == 4)
    {
        __m128i out8_0_U8_L0, out8_0_U8_L1;
        __m128i avg8_0_U8;

        for (y = 0; y < height; y += 2)
        {
            //--------
            //Line One
            //--------

            //List0
            inPixel0 = _mm_loadl_epi64((__m128i*)ref16_l0);
            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L0 = _mm_packus_epi16(inPixel1, inPixel1);

            //List1
            inPixel0 = _mm_loadl_epi64((__m128i*)ref16_l1);
            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L1 = _mm_packus_epi16(inPixel1, inPixel1);

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);

            *(uint32_t*)dst_ptr = _mm_cvtsi128_si32(avg8_0_U8);

            //--------
            //Line Two
            //--------

            //List0
            inPixel0 = _mm_loadl_epi64((__m128i*)(ref16_l0 + ref_l0_stride));
            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L0 = _mm_packus_epi16(inPixel1, inPixel1);

            //List1

            inPixel0 = _mm_loadl_epi64((__m128i*)(ref16_l1 + ref_l1_stride));
            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L1 = _mm_packus_epi16(inPixel1, inPixel1);

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);

            *(uint32_t*)(dst_ptr + dst_stride) = _mm_cvtsi128_si32(avg8_0_U8);

            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;

        }

    }
    else if (width == 8)
    {

        __m128i out8_0_U8_L0, out8_0_U8_L1, out8_2_U8_L0, out8_2_U8_L1;
        __m128i avg8_0_U8, avg8_2_U8;

        for (y = 0; y < height; y += 2)
        {
            //--------
            //Line One
            //--------

            //List0

            inPixel0 = _mm_loadu_si128((__m128i*) ref16_l0);

            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L0 = _mm_packus_epi16(inPixel1, inPixel1);

            //List1

            inPixel0 = _mm_loadu_si128((__m128i*) ref16_l1);

            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L1 = _mm_packus_epi16(inPixel1, inPixel1);

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);

            _mm_storel_epi64((__m128i*) dst_ptr, avg8_0_U8);


            //--------
            //Line Two
            //--------

            //List0

            inPixel0 = _mm_loadu_si128((__m128i*)(ref16_l0 + ref_l0_stride));

            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_2_U8_L0 = _mm_packus_epi16(inPixel1, inPixel1);

            //List1

            inPixel0 = _mm_loadu_si128((__m128i*)(ref16_l1 + ref_l1_stride));

            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_2_U8_L1 = _mm_packus_epi16(inPixel1, inPixel1);

            //AVG
            avg8_2_U8 = _mm_avg_epu8(out8_2_U8_L0, out8_2_U8_L1);

            _mm_storel_epi64((__m128i*)(dst_ptr + dst_stride), avg8_2_U8);


            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;
        }

    }
    else if (width == 16)
    {
        __m128i inPixel4, inPixel5;
        __m128i out8_0_U8_L0, out8_0_U8_L1, out8_2_U8_L0, out8_2_U8_L1;
        __m128i avg8_0_U8, avg8_2_U8;

        for (y = 0; y < height; y += 2)
        {
            //--------
            //Line One
            //--------

            //List0

            inPixel0 = _mm_loadu_si128((__m128i*)  ref16_l0);
            inPixel1 = _mm_loadu_si128((__m128i*) (ref16_l0 + 8));

            out8_0_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));

            //List1

            inPixel0 = _mm_loadu_si128((__m128i*) ref16_l1);
            inPixel1 = _mm_loadu_si128((__m128i*)(ref16_l1 + 8));

            out8_0_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));


            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);

            _mm_storeu_si128((__m128i*) dst_ptr, avg8_0_U8);


            //--------
            //Line Two
            //--------

            //List0

            inPixel4 = _mm_loadu_si128((__m128i*) (ref16_l0 + ref_l0_stride));
            inPixel5 = _mm_loadu_si128((__m128i*) (ref16_l0 + ref_l0_stride + 8));

            out8_2_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));

            //List1

            inPixel4 = _mm_loadu_si128((__m128i*) (ref16_l1 + ref_l1_stride));
            inPixel5 = _mm_loadu_si128((__m128i*) (ref16_l1 + ref_l1_stride + 8));

            out8_2_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));


            //AVG
            avg8_2_U8 = _mm_avg_epu8(out8_2_U8_L0, out8_2_U8_L1);

            _mm_storeu_si128((__m128i*)(dst_ptr + dst_stride), avg8_2_U8);

            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;

        }

    }
    else if (width == 32)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i out8_0_U8_L0, out8_1_U8_L0, out8_2_U8_L0, out8_3_U8_L0;
        __m128i out8_0_U8_L1, out8_1_U8_L1, out8_2_U8_L1, out8_3_U8_L1;
        __m128i avg8_0_U8, avg8_1_U8, avg8_2_U8, avg8_3_U8;


        for (y = 0; y < height; y += 2)
        {
            //--------
            //Line One
            //--------

            //List0

            inPixel0 = _mm_loadu_si128((__m128i*)  ref16_l0);
            inPixel1 = _mm_loadu_si128((__m128i*) (ref16_l0 + 8));
            inPixel2 = _mm_loadu_si128((__m128i*) (ref16_l0 + 16));
            inPixel3 = _mm_loadu_si128((__m128i*) (ref16_l0 + 24));

            out8_0_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));
            out8_1_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel2, 2), _mm_srli_epi16(inPixel3, 2));

            //List1

            inPixel0 = _mm_loadu_si128((__m128i*) ref16_l1);
            inPixel1 = _mm_loadu_si128((__m128i*)(ref16_l1 + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(ref16_l1 + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(ref16_l1 + 24));

            out8_0_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));
            out8_1_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel2, 2), _mm_srli_epi16(inPixel3, 2));

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);
            avg8_1_U8 = _mm_avg_epu8(out8_1_U8_L0, out8_1_U8_L1);

            _mm_storeu_si128((__m128i*) dst_ptr, avg8_0_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + 16), avg8_1_U8);


            //--------
            //Line Two
            //--------

            //List0

            inPixel4 = _mm_loadu_si128((__m128i*) (ref16_l0 + ref_l0_stride));
            inPixel5 = _mm_loadu_si128((__m128i*) (ref16_l0 + ref_l0_stride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (ref16_l0 + ref_l0_stride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (ref16_l0 + ref_l0_stride + 24));

            out8_2_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));
            out8_3_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel6, 2), _mm_srli_epi16(inPixel7, 2));

            //List1

            inPixel4 = _mm_loadu_si128((__m128i*) (ref16_l1 + ref_l1_stride));
            inPixel5 = _mm_loadu_si128((__m128i*) (ref16_l1 + ref_l1_stride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (ref16_l1 + ref_l1_stride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (ref16_l1 + ref_l1_stride + 24));

            out8_2_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));
            out8_3_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel6, 2), _mm_srli_epi16(inPixel7, 2));

            //AVG
            avg8_2_U8 = _mm_avg_epu8(out8_2_U8_L0, out8_2_U8_L1);
            avg8_3_U8 = _mm_avg_epu8(out8_3_U8_L0, out8_3_U8_L1);

            _mm_storeu_si128((__m128i*)(dst_ptr + dst_stride), avg8_2_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + dst_stride + 16), avg8_3_U8);

            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;

        }

    }
    else if (width == 64)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i out8_0_U8_L0, out8_1_U8_L0, out8_2_U8_L0, out8_3_U8_L0;
        __m128i out8_0_U8_L1, out8_1_U8_L1, out8_2_U8_L1, out8_3_U8_L1;
        __m128i avg8_0_U8, avg8_1_U8, avg8_2_U8, avg8_3_U8;

        for (y = 0; y < height; ++y)
        {
            //List0

            inPixel0 = _mm_loadu_si128((__m128i*) ref16_l0);
            inPixel1 = _mm_loadu_si128((__m128i*)(ref16_l0 + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(ref16_l0 + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(ref16_l0 + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(ref16_l0 + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(ref16_l0 + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(ref16_l0 + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(ref16_l0 + 56));



            out8_0_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));
            out8_1_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel2, 2), _mm_srli_epi16(inPixel3, 2));
            out8_2_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));
            out8_3_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel6, 2), _mm_srli_epi16(inPixel7, 2));


            //List1

            inPixel0 = _mm_loadu_si128((__m128i*) ref16_l1);
            inPixel1 = _mm_loadu_si128((__m128i*)(ref16_l1 + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(ref16_l1 + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(ref16_l1 + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(ref16_l1 + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(ref16_l1 + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(ref16_l1 + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(ref16_l1 + 56));


            //Note: old Version used to use _mm_and_si128 to mask the MSB bits of the pixels
            out8_0_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));
            out8_1_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel2, 2), _mm_srli_epi16(inPixel3, 2));
            out8_2_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));
            out8_3_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel6, 2), _mm_srli_epi16(inPixel7, 2));

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);
            avg8_1_U8 = _mm_avg_epu8(out8_1_U8_L0, out8_1_U8_L1);
            avg8_2_U8 = _mm_avg_epu8(out8_2_U8_L0, out8_2_U8_L1);
            avg8_3_U8 = _mm_avg_epu8(out8_3_U8_L0, out8_3_U8_L1);

            _mm_storeu_si128((__m128i*) dst_ptr, avg8_0_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + 16), avg8_1_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + 32), avg8_2_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + 48), avg8_3_U8);


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
void eb_enc_msb_pack2d_sse2_intrin(
    uint8_t     *in8_bit_buffer,
    uint32_t     in8_stride,
    uint8_t     *inn_bit_buffer,
    uint16_t    *out16_bit_buffer,
    uint32_t     inn_stride,
    uint32_t     out_stride,
    uint32_t     width,
    uint32_t     height)
{
    uint32_t count_width, count_height;

    if (width == 4) {

        for (count_height = 0; count_height < height; count_height += 2) {
            _mm_storel_epi64((__m128i*)(out16_bit_buffer), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(inn_bit_buffer)),
                _mm_cvtsi32_si128(*(uint32_t*)(in8_bit_buffer))), 6));
            _mm_storel_epi64((__m128i*)(out16_bit_buffer + out_stride), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(inn_bit_buffer + inn_stride)),
                _mm_cvtsi32_si128(*(uint32_t*)(in8_bit_buffer + in8_stride))), 6));
            out16_bit_buffer += (out_stride << 1);
            in8_bit_buffer += (in8_stride << 1);
            inn_bit_buffer += (inn_stride << 1);
        }
    }
    else if (width == 8) {

        for (count_height = 0; count_height < height; count_height += 2) {

            _mm_storeu_si128((__m128i*)(out16_bit_buffer), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(inn_bit_buffer)),
                _mm_loadl_epi64((__m128i*)(in8_bit_buffer))), 6));
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + out_stride), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(inn_bit_buffer + inn_stride)),
                _mm_loadl_epi64((__m128i*)(in8_bit_buffer + in8_stride))), 6));
            out16_bit_buffer += (out_stride << 1);
            in8_bit_buffer += (in8_stride << 1);
            inn_bit_buffer += (inn_stride << 1);
        }
    }
    else if (width == 16) {
        __m128i outPixel1, outPixel2, outPixel3, outPixel4, innBitBuffer_lo, innBitBuffer_hi, in8BitBuffer_lo, in8BitBuffer_hi;

        for (count_height = 0; count_height < height; count_height += 2) {

            innBitBuffer_lo = _mm_loadu_si128((__m128i *)inn_bit_buffer);
            innBitBuffer_hi = _mm_loadu_si128((__m128i *)(inn_bit_buffer + inn_stride));
            in8BitBuffer_lo = _mm_loadu_si128((__m128i *)in8_bit_buffer);
            in8BitBuffer_hi = _mm_loadu_si128((__m128i *)(in8_bit_buffer + in8_stride));

            outPixel1 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer_lo, in8BitBuffer_lo), 6);
            outPixel2 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer_lo, in8BitBuffer_lo), 6);
            outPixel3 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer_hi, in8BitBuffer_hi), 6);
            outPixel4 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer_hi, in8BitBuffer_hi), 6);

            _mm_storeu_si128((__m128i*)out16_bit_buffer, outPixel1);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 8), outPixel2);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + out_stride), outPixel3);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + out_stride + 8), outPixel4);

            in8_bit_buffer += (in8_stride << 1);
            inn_bit_buffer += (inn_stride << 1);
            out16_bit_buffer += (out_stride << 1);
        }
    }
    else if (width == 32) {
        __m128i innBitBuffer1, innBitBuffer2, innBitBuffer3, innBitBuffer4, in8BitBuffer1, in8BitBuffer2, in8BitBuffer3, in8BitBuffer4;
        __m128i outPixel1, outPixel2, outPixel3, outPixel4, outPixel5, outPixel6, outPixel7, outPixel8;

        for (count_height = 0; count_height < height; count_height += 2)
        {
            innBitBuffer1 = _mm_loadu_si128((__m128i *)inn_bit_buffer);
            innBitBuffer2 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + 16));
            innBitBuffer3 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + inn_stride));
            innBitBuffer4 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + inn_stride + 16));

            in8BitBuffer1 = _mm_loadu_si128((__m128i *)in8_bit_buffer);
            in8BitBuffer2 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + 16));
            in8BitBuffer3 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + in8_stride));
            in8BitBuffer4 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + in8_stride + 16));

            outPixel1 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer1, in8BitBuffer1), 6);
            outPixel2 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer1, in8BitBuffer1), 6);
            outPixel3 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer2, in8BitBuffer2), 6);
            outPixel4 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer2, in8BitBuffer2), 6);
            outPixel5 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer3, in8BitBuffer3), 6);
            outPixel6 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer3, in8BitBuffer3), 6);
            outPixel7 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer4, in8BitBuffer4), 6);
            outPixel8 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer4, in8BitBuffer4), 6);

            _mm_storeu_si128((__m128i*)out16_bit_buffer, outPixel1);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 8), outPixel2);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 16), outPixel3);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 24), outPixel4);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + out_stride), outPixel5);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + out_stride + 8), outPixel6);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + out_stride + 16), outPixel7);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + out_stride + 24), outPixel8);

            in8_bit_buffer += (in8_stride << 1);
            inn_bit_buffer += (inn_stride << 1);
            out16_bit_buffer += (out_stride << 1);
        }
    }
    else if (width == 64) {

        __m128i innBitBuffer1, innBitBuffer2, innBitBuffer3, innBitBuffer4, in8BitBuffer1, in8BitBuffer2, in8BitBuffer3, in8BitBuffer4;
        __m128i outPixel1, outPixel2, outPixel3, outPixel4, outPixel5, outPixel6, outPixel7, outPixel8;

        for (count_height = 0; count_height < height; ++count_height)
        {
            innBitBuffer1 = _mm_loadu_si128((__m128i *)inn_bit_buffer);
            innBitBuffer2 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + 16));
            innBitBuffer3 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + 32));
            innBitBuffer4 = _mm_loadu_si128((__m128i *)(inn_bit_buffer + 48));

            in8BitBuffer1 = _mm_loadu_si128((__m128i *)in8_bit_buffer);
            in8BitBuffer2 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + 16));
            in8BitBuffer3 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + 32));
            in8BitBuffer4 = _mm_loadu_si128((__m128i *)(in8_bit_buffer + 48));

            outPixel1 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer1, in8BitBuffer1), 6);
            outPixel2 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer1, in8BitBuffer1), 6);
            outPixel3 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer2, in8BitBuffer2), 6);
            outPixel4 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer2, in8BitBuffer2), 6);
            outPixel5 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer3, in8BitBuffer3), 6);
            outPixel6 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer3, in8BitBuffer3), 6);
            outPixel7 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer4, in8BitBuffer4), 6);
            outPixel8 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer4, in8BitBuffer4), 6);

            _mm_storeu_si128((__m128i*)out16_bit_buffer, outPixel1);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 8), outPixel2);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 16), outPixel3);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 24), outPixel4);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 32), outPixel5);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 40), outPixel6);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 48), outPixel7);
            _mm_storeu_si128((__m128i*)(out16_bit_buffer + 56), outPixel8);

            in8_bit_buffer += in8_stride;
            inn_bit_buffer += inn_stride;
            out16_bit_buffer += out_stride;
        }
    }
    else {
        uint32_t innStrideDiff = (inn_stride << 1) - width;
        uint32_t in8StrideDiff = (in8_stride << 1) - width;
        uint32_t outStrideDiff = (out_stride << 1) - width;

        if (!(width & 7)) {
            for (count_height = 0; count_height < height; count_height += 2) {
                for (count_width = 0; count_width < width; count_width += 8) {
                    _mm_storeu_si128((__m128i*)(out16_bit_buffer), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(inn_bit_buffer)),
                        _mm_loadl_epi64((__m128i*)(in8_bit_buffer))), 6));
                    _mm_storeu_si128((__m128i*)(out16_bit_buffer + out_stride), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(inn_bit_buffer + inn_stride)),
                        _mm_loadl_epi64((__m128i*)(in8_bit_buffer + in8_stride))), 6));
                    out16_bit_buffer += 8;
                    in8_bit_buffer += 8;
                    inn_bit_buffer += 8;
                }
                in8_bit_buffer += in8StrideDiff;
                inn_bit_buffer += innStrideDiff;
                out16_bit_buffer += outStrideDiff;
            }
        }
        else {
            for (count_height = 0; count_height < height; count_height += 2) {
                for (count_width = 0; count_width < width; count_width += 4) {
                    _mm_storel_epi64((__m128i*)(out16_bit_buffer), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(inn_bit_buffer)),
                        _mm_cvtsi32_si128(*(uint32_t*)(in8_bit_buffer))), 6));
                    _mm_storel_epi64((__m128i*)(out16_bit_buffer + out_stride), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(inn_bit_buffer + inn_stride)),
                        _mm_cvtsi32_si128(*(uint32_t*)(in8_bit_buffer + in8_stride))), 6));
                    out16_bit_buffer += 4;
                    in8_bit_buffer += 4;
                    inn_bit_buffer += 4;
                }
                in8_bit_buffer += in8StrideDiff;
                inn_bit_buffer += innStrideDiff;
                out16_bit_buffer += outStrideDiff;
            }
        }
    }
}