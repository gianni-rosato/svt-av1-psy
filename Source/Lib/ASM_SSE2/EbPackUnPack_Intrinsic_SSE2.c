/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPackUnPack_SSE2.h"


#include <emmintrin.h>

#include "EbDefinitions.h"

/****************************************************************************************
EB_ENC_msbUnPack2D_SSE2_INTRIN
******************************************************************************************/

void EB_ENC_msbUnPack2D_SSE2_INTRIN(
    uint16_t      *in16BitBuffer,
    uint32_t       inStride,
    uint8_t       *out8BitBuffer,
    uint8_t       *outnBitBuffer,
    uint32_t       out8Stride,
    uint32_t       outnStride,
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
            inPixel0 = _mm_loadl_epi64((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadl_epi64((__m128i*)(in16BitBuffer + inStride));

            tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
            tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);

            tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
            tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

            *(uint32_t*)outnBitBuffer = _mm_cvtsi128_si32(tempPixel0_U8);
            *(uint32_t*)(outnBitBuffer + outnStride) = _mm_cvtsi128_si32(tempPixel1_U8);
            *(uint32_t*)out8BitBuffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
            *(uint32_t*)(out8BitBuffer + out8Stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

            outnBitBuffer += 2 * outnStride;
            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }
    }
    else if (width == 8)
    {
        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));

            tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
            tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);

            tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
            tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

            _mm_storel_epi64((__m128i*)outnBitBuffer, tempPixel0_U8);
            _mm_storel_epi64((__m128i*)(outnBitBuffer + outnStride), tempPixel1_U8);
            _mm_storel_epi64((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
            _mm_storel_epi64((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

            outnBitBuffer += 2 * outnStride;
            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }
    }
    else if (width == 16)
    {
        __m128i inPixel2, inPixel3;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));
            inPixel3 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 8));

            tempPixel0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            tempPixel1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));

            inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));

            _mm_storeu_si128((__m128i*)outnBitBuffer, tempPixel0_U8);
            _mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride), tempPixel1_U8);
            _mm_storeu_si128((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

            outnBitBuffer += 2 * outnStride;
            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }
    }
    else if (width == 32)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i outn0_U8, outn1_U8, outn2_U8, outn3_U8, out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 24));

            outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
            outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
            outn3_U8 = _mm_packus_epi16(_mm_and_si128(inPixel6, xmm_3), _mm_and_si128(inPixel7, xmm_3));

            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

            _mm_storeu_si128((__m128i*)outnBitBuffer, outn0_U8);
            _mm_storeu_si128((__m128i*)(outnBitBuffer + 16), outn1_U8);
            _mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride), outn2_U8);
            _mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride + 16), outn3_U8);

            _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride + 16), out8_3_U8);

            outnBitBuffer += 2 * outnStride;
            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }
    }
    else if (width == 64)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i outn0_U8, outn1_U8, outn2_U8, outn3_U8, out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; ++y)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 56));

            outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
            outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
            outn3_U8 = _mm_packus_epi16(_mm_and_si128(inPixel6, xmm_3), _mm_and_si128(inPixel7, xmm_3));

            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

            _mm_storeu_si128((__m128i*)outnBitBuffer, outn0_U8);
            _mm_storeu_si128((__m128i*)(outnBitBuffer + 16), outn1_U8);
            _mm_storeu_si128((__m128i*)(outnBitBuffer + 32), outn2_U8);
            _mm_storeu_si128((__m128i*)(outnBitBuffer + 48), outn3_U8);

            _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 32), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 48), out8_3_U8);

            outnBitBuffer += outnStride;
            out8BitBuffer += out8Stride;
            in16BitBuffer += inStride;
        }
    }
    else
    {
        uint32_t inStrideDiff = (2 * inStride) - width;
        uint32_t out8StrideDiff = (2 * out8Stride) - width;
        uint32_t outnStrideDiff = (2 * outnStride) - width;

        uint32_t inStrideDiff64 = inStride - width;
        uint32_t out8StrideDiff64 = out8Stride - width;
        uint32_t outnStrideDiff64 = outnStride - width;

        if (!(width & 63))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i outn0_U8, outn1_U8, outn2_U8, outn3_U8, out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 1) {
                for (y = 0; y < width; y += 64) {

                    inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 32));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 40));
                    inPixel6 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 48));
                    inPixel7 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 56));

                    outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
                    outn3_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel6, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel7, xmm_3), 6));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    _mm_storeu_si128((__m128i*)outnBitBuffer, outn0_U8);
                    _mm_storeu_si128((__m128i*)(outnBitBuffer + 16), outn1_U8);
                    _mm_storeu_si128((__m128i*)(outnBitBuffer + 32), outn2_U8);
                    _mm_storeu_si128((__m128i*)(outnBitBuffer + 48), outn3_U8);

                    _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 32), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 48), out8_3_U8);

                    outnBitBuffer += 64;
                    out8BitBuffer += 64;
                    in16BitBuffer += 64;
                }
                in16BitBuffer += inStrideDiff64;
                outnBitBuffer += outnStrideDiff64;
                out8BitBuffer += out8StrideDiff64;
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
                    inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride + 8));
                    inPixel6 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 16));
                    inPixel7 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 24));

                    outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
                    outn3_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel6, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel7, xmm_3), 6));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    _mm_storeu_si128((__m128i*)outnBitBuffer, outn0_U8);
                    _mm_storeu_si128((__m128i*)(outnBitBuffer + 16), outn1_U8);
                    _mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride), outn2_U8);
                    _mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride + 16), outn3_U8);

                    _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride + 16), out8_3_U8);

                    outnBitBuffer += 32;
                    out8BitBuffer += 32;
                    in16BitBuffer += 32;
                }
                in16BitBuffer += inStrideDiff;
                outnBitBuffer += outnStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else if (!(width & 15))
        {
            __m128i inPixel2, inPixel3;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 16)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));
                    inPixel3 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 8));

                    tempPixel0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    tempPixel1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));

                    _mm_storeu_si128((__m128i*)outnBitBuffer, tempPixel0_U8);
                    _mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride), tempPixel1_U8);
                    _mm_storeu_si128((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

                    outnBitBuffer += 16;
                    out8BitBuffer += 16;
                    in16BitBuffer += 16;
                }
                in16BitBuffer += inStrideDiff;
                outnBitBuffer += outnStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else if (!(width & 7))
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 8)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));

                    tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
                    tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);

                    tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
                    tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    _mm_storel_epi64((__m128i*)outnBitBuffer, tempPixel0_U8);
                    _mm_storel_epi64((__m128i*)(outnBitBuffer + outnStride), tempPixel1_U8);
                    _mm_storel_epi64((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
                    _mm_storel_epi64((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

                    outnBitBuffer += 8;
                    out8BitBuffer += 8;
                    in16BitBuffer += 8;
                }
                in16BitBuffer += inStrideDiff;
                outnBitBuffer += outnStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 4)
                {
                    inPixel0 = _mm_loadl_epi64((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadl_epi64((__m128i*)(in16BitBuffer + inStride));

                    tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
                    tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);

                    tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
                    tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    *(uint32_t*)outnBitBuffer = _mm_cvtsi128_si32(tempPixel0_U8);
                    *(uint32_t*)(outnBitBuffer + outnStride) = _mm_cvtsi128_si32(tempPixel1_U8);
                    *(uint32_t*)out8BitBuffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
                    *(uint32_t*)(out8BitBuffer + out8Stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

                    outnBitBuffer += 4;
                    out8BitBuffer += 4;
                    in16BitBuffer += 4;
                }
                in16BitBuffer += inStrideDiff;
                outnBitBuffer += outnStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
    }
    return;
}

void EB_ENC_UnPack8BitDataSafeSub_SSE2_INTRIN(
    uint16_t      *in16BitBuffer,
    uint32_t       inStride,
    uint8_t       *out8BitBuffer,
    uint32_t       out8Stride,
    uint32_t       width,
    uint32_t       height,
    EbBool      subPred
)
{

    uint32_t x, y;

    __m128i xmm_00FF, inPixel0, inPixel1, inPixel1_shftR_2_U8, inPixel0_shftR_2_U8, inPixel0_shftR_2, inPixel1_shftR_2;


    xmm_00FF = _mm_set1_epi16(0x00FF);

    if (width == 8)
    {
        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));


            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);


            _mm_storel_epi64((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
            _mm_storel_epi64((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }

        if (subPred) {
            in16BitBuffer -= (inStride >> 1);
            out8BitBuffer -= (out8Stride >> 1);
            inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            _mm_storel_epi64((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
        }



    }
    else if (width == 16)
    {
        __m128i inPixel2, inPixel3;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));
            inPixel3 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 8));


            inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));


            _mm_storeu_si128((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);


            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }

        if (subPred) {
            in16BitBuffer -= (inStride >> 1);
            out8BitBuffer -= (out8Stride >> 1);
            inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + 8));
            inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            _mm_storeu_si128((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
        }

    }
    else if (width == 32)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i  out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 24));


            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));


            _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride + 16), out8_3_U8);


            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }

        if (subPred) {
            in16BitBuffer -= (inStride >> 1);
            out8BitBuffer -= (out8Stride >> 1);
            inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);

        }

    }
    else if (width == 64)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i  out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; ++y)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 56));


            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));


            _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 32), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 48), out8_3_U8);


            out8BitBuffer += out8Stride;
            in16BitBuffer += inStride;
        }

        if (subPred) {
            in16BitBuffer -= (inStride >> 1);
            out8BitBuffer -= (out8Stride >> 1);
            inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 56));


            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));


            _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 32), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 48), out8_3_U8);
        }

    }
    else
    {
        uint32_t inStrideDiff = (2 * inStride) - width;
        uint32_t out8StrideDiff = (2 * out8Stride) - width;

        uint32_t inStrideDiff64 = inStride - width;
        uint32_t out8StrideDiff64 = out8Stride - width;

        if (!(width & 63))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 1) {
                for (y = 0; y < width; y += 64) {

                    inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 32));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 40));
                    inPixel6 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 48));
                    inPixel7 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 56));


                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));


                    _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 32), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 48), out8_3_U8);

                    out8BitBuffer += 64;
                    in16BitBuffer += 64;
                }
                in16BitBuffer += inStrideDiff64;
                out8BitBuffer += out8StrideDiff64;
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
                    inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride + 8));
                    inPixel6 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 16));
                    inPixel7 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 24));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride + 16), out8_3_U8);

                    out8BitBuffer += 32;
                    in16BitBuffer += 32;
                }
                in16BitBuffer += inStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else if (!(width & 15))
        {
            __m128i inPixel2, inPixel3;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 16)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));
                    inPixel3 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 8));

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));


                    _mm_storeu_si128((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

                    out8BitBuffer += 16;
                    in16BitBuffer += 16;
                }
                in16BitBuffer += inStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else if (!(width & 7))
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 8)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    _mm_storel_epi64((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
                    _mm_storel_epi64((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

                    out8BitBuffer += 8;
                    in16BitBuffer += 8;
                }
                in16BitBuffer += inStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 4)
                {
                    inPixel0 = _mm_loadl_epi64((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadl_epi64((__m128i*)(in16BitBuffer + inStride));


                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    *(uint32_t*)out8BitBuffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
                    *(uint32_t*)(out8BitBuffer + out8Stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

                    out8BitBuffer += 4;
                    in16BitBuffer += 4;
                }
                in16BitBuffer += inStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }

    }


    return;
}

void EB_ENC_UnPack8BitData_SSE2_INTRIN(
    uint16_t      *in16BitBuffer,
    uint32_t       inStride,
    uint8_t       *out8BitBuffer,
    uint32_t       out8Stride,
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
            inPixel0 = _mm_loadl_epi64((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadl_epi64((__m128i*)(in16BitBuffer + inStride));

            //tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
            //tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);
            //
            //tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
            //tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

            //*(uint32_t*)outnBitBuffer = _mm_cvtsi128_si32(tempPixel0_U8);
            //*(uint32_t*)(outnBitBuffer + outnStride) = _mm_cvtsi128_si32(tempPixel1_U8);
            *(uint32_t*)out8BitBuffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
            *(uint32_t*)(out8BitBuffer + out8Stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

            //outnBitBuffer += 2 * outnStride;
            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }
    }
    else if (width == 8)
    {
        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));

            //tempPixel0 = _mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6);
            //tempPixel1 = _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6);
            //
            //tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
            //tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

            inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
            inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

            inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
            inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

            //_mm_storel_epi64((__m128i*)outnBitBuffer, tempPixel0_U8);
            //_mm_storel_epi64((__m128i*)(outnBitBuffer + outnStride), tempPixel1_U8);
            _mm_storel_epi64((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
            _mm_storel_epi64((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

            //outnBitBuffer += 2 * outnStride;
            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }
    }
    else if (width == 16)
    {
        __m128i inPixel2, inPixel3;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));
            inPixel3 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 8));

            //tempPixel0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            //tempPixel1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));

            inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));

            //_mm_storeu_si128((__m128i*)outnBitBuffer, tempPixel0_U8);
            //_mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride), tempPixel1_U8);
            _mm_storeu_si128((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

            //outnBitBuffer += 2 * outnStride;
            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }
    }
    else if (width == 32)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i /*outn0_U8, outn1_U8, outn2_U8, outn3_U8,*/ out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; y += 2)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 24));

            //outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            //outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
            //outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
            //outn3_U8 = _mm_packus_epi16(_mm_and_si128(inPixel6, xmm_3), _mm_and_si128(inPixel7, xmm_3));

            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

            //_mm_storeu_si128((__m128i*)outnBitBuffer, outn0_U8);
            //_mm_storeu_si128((__m128i*)(outnBitBuffer + 16), outn1_U8);
            //_mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride), outn2_U8);
            //_mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride + 16), outn3_U8);

            _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride + 16), out8_3_U8);

            //outnBitBuffer += 2 * outnStride;
            out8BitBuffer += 2 * out8Stride;
            in16BitBuffer += 2 * inStride;
        }
    }
    else if (width == 64)
    {
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i /*outn0_U8, outn1_U8, outn2_U8, outn3_U8,*/ out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

        for (y = 0; y < height; ++y)
        {
            inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
            inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 56));

            //outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
            //outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
            //outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
            //outn3_U8 = _mm_packus_epi16(_mm_and_si128(inPixel6, xmm_3), _mm_and_si128(inPixel7, xmm_3));

            out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
            out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
            out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
            out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

            //_mm_storeu_si128((__m128i*)outnBitBuffer, outn0_U8);
            //_mm_storeu_si128((__m128i*)(outnBitBuffer + 16), outn1_U8);
            //_mm_storeu_si128((__m128i*)(outnBitBuffer + 32), outn2_U8);
            //_mm_storeu_si128((__m128i*)(outnBitBuffer + 48), outn3_U8);

            _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 32), out8_2_U8);
            _mm_storeu_si128((__m128i*)(out8BitBuffer + 48), out8_3_U8);

            //outnBitBuffer += outnStride;
            out8BitBuffer += out8Stride;
            in16BitBuffer += inStride;
        }
    }
    else
    {
        uint32_t inStrideDiff = (2 * inStride) - width;
        uint32_t out8StrideDiff = (2 * out8Stride) - width;
        //uint32_t outnStrideDiff = (2 * outnStride) - width;

        uint32_t inStrideDiff64 = inStride - width;
        uint32_t out8StrideDiff64 = out8Stride - width;
        //uint32_t outnStrideDiff64 = outnStride - width;

        if (!(width & 63))
        {
            __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
            __m128i /*outn0_U8, outn1_U8, outn2_U8, outn3_U8, */out8_0_U8, out8_1_U8, out8_2_U8, out8_3_U8;

            for (x = 0; x < height; x += 1) {
                for (y = 0; y < width; y += 64) {

                    inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 32));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 40));
                    inPixel6 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 48));
                    inPixel7 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 56));

                    //outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    //outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    //outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
                    //outn3_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel6, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel7, xmm_3), 6));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    //_mm_storeu_si128((__m128i*)outnBitBuffer, outn0_U8);
                    //_mm_storeu_si128((__m128i*)(outnBitBuffer + 16), outn1_U8);
                    //_mm_storeu_si128((__m128i*)(outnBitBuffer + 32), outn2_U8);
                    //_mm_storeu_si128((__m128i*)(outnBitBuffer + 48), outn3_U8);

                    _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 32), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 48), out8_3_U8);

                    //outnBitBuffer += 64;
                    out8BitBuffer += 64;
                    in16BitBuffer += 64;
                }
                in16BitBuffer += inStrideDiff64;
                //outnBitBuffer += outnStrideDiff64;
                out8BitBuffer += out8StrideDiff64;
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
                    inPixel0 = _mm_loadu_si128((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 16));
                    inPixel3 = _mm_loadu_si128((__m128i*)(in16BitBuffer + 24));
                    inPixel4 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride));
                    inPixel5 = _mm_loadu_si128((__m128i*)(in16BitBuffer + inStride + 8));
                    inPixel6 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 16));
                    inPixel7 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 24));

                    //outn0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    //outn1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    //outn2_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel4, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel5, xmm_3), 6));
                    //outn3_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel6, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel7, xmm_3), 6));

                    out8_0_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    out8_1_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));
                    out8_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel4, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel5, 2), xmm_00FF));
                    out8_3_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel6, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel7, 2), xmm_00FF));

                    //_mm_storeu_si128((__m128i*)outnBitBuffer, outn0_U8);
                    //_mm_storeu_si128((__m128i*)(outnBitBuffer + 16), outn1_U8);
                    //_mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride), outn2_U8);
                    //_mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride + 16), outn3_U8);

                    _mm_storeu_si128((__m128i*)out8BitBuffer, out8_0_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + 16), out8_1_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), out8_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride + 16), out8_3_U8);

                    //outnBitBuffer += 32;
                    out8BitBuffer += 32;
                    in16BitBuffer += 32;
                }
                in16BitBuffer += inStrideDiff;
                //outnBitBuffer += outnStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else if (!(width & 15))
        {
            __m128i inPixel2, inPixel3;

            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 16)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + 8));
                    inPixel2 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));
                    inPixel3 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride + 8));

                    //tempPixel0_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel0, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel1, xmm_3), 6));
                    //tempPixel1_U8 = _mm_packus_epi16(_mm_slli_epi16(_mm_and_si128(inPixel2, xmm_3), 6), _mm_slli_epi16(_mm_and_si128(inPixel3, xmm_3), 6));
                    //
                    inPixel0_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF));
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(_mm_and_si128(_mm_srli_epi16(inPixel2, 2), xmm_00FF), _mm_and_si128(_mm_srli_epi16(inPixel3, 2), xmm_00FF));

                    //_mm_storeu_si128((__m128i*)outnBitBuffer, tempPixel0_U8);
                    //_mm_storeu_si128((__m128i*)(outnBitBuffer + outnStride), tempPixel1_U8);
                    _mm_storeu_si128((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
                    _mm_storeu_si128((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

                    //outnBitBuffer += 16;
                    out8BitBuffer += 16;
                    in16BitBuffer += 16;
                }
                in16BitBuffer += inStrideDiff;
                //outnBitBuffer += outnStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else if (!(width & 7))
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 8)
                {
                    inPixel0 = _mm_loadu_si128((__m128i*) in16BitBuffer);
                    inPixel1 = _mm_loadu_si128((__m128i*) (in16BitBuffer + inStride));

                    //tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
                    //tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    //_mm_storel_epi64((__m128i*)outnBitBuffer, tempPixel0_U8);
                    //_mm_storel_epi64((__m128i*)(outnBitBuffer + outnStride), tempPixel1_U8);
                    _mm_storel_epi64((__m128i*)out8BitBuffer, inPixel0_shftR_2_U8);
                    _mm_storel_epi64((__m128i*)(out8BitBuffer + out8Stride), inPixel1_shftR_2_U8);

                    //outnBitBuffer += 8;
                    out8BitBuffer += 8;
                    in16BitBuffer += 8;
                }
                in16BitBuffer += inStrideDiff;
                //outnBitBuffer += outnStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
        else
        {
            for (x = 0; x < height; x += 2)
            {
                for (y = 0; y < width; y += 4)
                {
                    inPixel0 = _mm_loadl_epi64((__m128i*)in16BitBuffer);
                    inPixel1 = _mm_loadl_epi64((__m128i*)(in16BitBuffer + inStride));

                    //tempPixel0_U8 = _mm_packus_epi16(tempPixel0, tempPixel0);
                    //tempPixel1_U8 = _mm_packus_epi16(tempPixel1, tempPixel1);

                    inPixel0_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel0, 2), xmm_00FF);
                    inPixel1_shftR_2 = _mm_and_si128(_mm_srli_epi16(inPixel1, 2), xmm_00FF);

                    inPixel0_shftR_2_U8 = _mm_packus_epi16(inPixel0_shftR_2, inPixel0_shftR_2);
                    inPixel1_shftR_2_U8 = _mm_packus_epi16(inPixel1_shftR_2, inPixel1_shftR_2);

                    //*(uint32_t*)outnBitBuffer = _mm_cvtsi128_si32(tempPixel0_U8);
                    //*(uint32_t*)(outnBitBuffer + outnStride) = _mm_cvtsi128_si32(tempPixel1_U8);
                    *(uint32_t*)out8BitBuffer = _mm_cvtsi128_si32(inPixel0_shftR_2_U8);
                    *(uint32_t*)(out8BitBuffer + out8Stride) = _mm_cvtsi128_si32(inPixel1_shftR_2_U8);

                    //outnBitBuffer += 4;
                    out8BitBuffer += 4;
                    in16BitBuffer += 4;
                }
                in16BitBuffer += inStrideDiff;
                //outnBitBuffer += outnStrideDiff;
                out8BitBuffer += out8StrideDiff;
            }
        }
    }

    return;
}
void UnpackAvg_SSE2_INTRIN(
    uint16_t *ref16L0,
    uint32_t  refL0Stride,
    uint16_t *ref16L1,
    uint32_t  refL1Stride,
    uint8_t  *dstPtr,
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
            inPixel0 = _mm_loadl_epi64((__m128i*)ref16L0);
            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L0 = _mm_packus_epi16(inPixel1, inPixel1);

            //List1
            inPixel0 = _mm_loadl_epi64((__m128i*)ref16L1);
            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L1 = _mm_packus_epi16(inPixel1, inPixel1);

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);

            *(uint32_t*)dstPtr = _mm_cvtsi128_si32(avg8_0_U8);

            //--------
            //Line Two
            //--------

            //List0
            inPixel0 = _mm_loadl_epi64((__m128i*)(ref16L0 + refL0Stride));
            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L0 = _mm_packus_epi16(inPixel1, inPixel1);

            //List1

            inPixel0 = _mm_loadl_epi64((__m128i*)(ref16L1 + refL1Stride));
            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L1 = _mm_packus_epi16(inPixel1, inPixel1);

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);

            *(uint32_t*)(dstPtr + dst_stride) = _mm_cvtsi128_si32(avg8_0_U8);

            dstPtr += 2 * dst_stride;
            ref16L0 += 2 * refL0Stride;
            ref16L1 += 2 * refL1Stride;

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

            inPixel0 = _mm_loadu_si128((__m128i*) ref16L0);

            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L0 = _mm_packus_epi16(inPixel1, inPixel1);

            //List1

            inPixel0 = _mm_loadu_si128((__m128i*) ref16L1);

            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_0_U8_L1 = _mm_packus_epi16(inPixel1, inPixel1);

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);

            _mm_storel_epi64((__m128i*) dstPtr, avg8_0_U8);


            //--------
            //Line Two
            //--------

            //List0

            inPixel0 = _mm_loadu_si128((__m128i*)(ref16L0 + refL0Stride));

            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_2_U8_L0 = _mm_packus_epi16(inPixel1, inPixel1);

            //List1

            inPixel0 = _mm_loadu_si128((__m128i*)(ref16L1 + refL1Stride));

            inPixel1 = _mm_srli_epi16(inPixel0, 2);
            out8_2_U8_L1 = _mm_packus_epi16(inPixel1, inPixel1);

            //AVG
            avg8_2_U8 = _mm_avg_epu8(out8_2_U8_L0, out8_2_U8_L1);

            _mm_storel_epi64((__m128i*)(dstPtr + dst_stride), avg8_2_U8);


            dstPtr += 2 * dst_stride;
            ref16L0 += 2 * refL0Stride;
            ref16L1 += 2 * refL1Stride;
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

            inPixel0 = _mm_loadu_si128((__m128i*)  ref16L0);
            inPixel1 = _mm_loadu_si128((__m128i*) (ref16L0 + 8));

            out8_0_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));

            //List1

            inPixel0 = _mm_loadu_si128((__m128i*) ref16L1);
            inPixel1 = _mm_loadu_si128((__m128i*)(ref16L1 + 8));

            out8_0_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));


            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);

            _mm_storeu_si128((__m128i*) dstPtr, avg8_0_U8);


            //--------
            //Line Two
            //--------

            //List0

            inPixel4 = _mm_loadu_si128((__m128i*) (ref16L0 + refL0Stride));
            inPixel5 = _mm_loadu_si128((__m128i*) (ref16L0 + refL0Stride + 8));

            out8_2_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));

            //List1

            inPixel4 = _mm_loadu_si128((__m128i*) (ref16L1 + refL1Stride));
            inPixel5 = _mm_loadu_si128((__m128i*) (ref16L1 + refL1Stride + 8));

            out8_2_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));


            //AVG
            avg8_2_U8 = _mm_avg_epu8(out8_2_U8_L0, out8_2_U8_L1);

            _mm_storeu_si128((__m128i*)(dstPtr + dst_stride), avg8_2_U8);

            dstPtr += 2 * dst_stride;
            ref16L0 += 2 * refL0Stride;
            ref16L1 += 2 * refL1Stride;

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

            inPixel0 = _mm_loadu_si128((__m128i*)  ref16L0);
            inPixel1 = _mm_loadu_si128((__m128i*) (ref16L0 + 8));
            inPixel2 = _mm_loadu_si128((__m128i*) (ref16L0 + 16));
            inPixel3 = _mm_loadu_si128((__m128i*) (ref16L0 + 24));

            out8_0_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));
            out8_1_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel2, 2), _mm_srli_epi16(inPixel3, 2));

            //List1

            inPixel0 = _mm_loadu_si128((__m128i*) ref16L1);
            inPixel1 = _mm_loadu_si128((__m128i*)(ref16L1 + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(ref16L1 + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(ref16L1 + 24));

            out8_0_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));
            out8_1_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel2, 2), _mm_srli_epi16(inPixel3, 2));

            //AVG
            avg8_0_U8 = _mm_avg_epu8(out8_0_U8_L0, out8_0_U8_L1);
            avg8_1_U8 = _mm_avg_epu8(out8_1_U8_L0, out8_1_U8_L1);

            _mm_storeu_si128((__m128i*) dstPtr, avg8_0_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + 16), avg8_1_U8);


            //--------
            //Line Two
            //--------

            //List0

            inPixel4 = _mm_loadu_si128((__m128i*) (ref16L0 + refL0Stride));
            inPixel5 = _mm_loadu_si128((__m128i*) (ref16L0 + refL0Stride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (ref16L0 + refL0Stride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (ref16L0 + refL0Stride + 24));

            out8_2_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));
            out8_3_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel6, 2), _mm_srli_epi16(inPixel7, 2));

            //List1

            inPixel4 = _mm_loadu_si128((__m128i*) (ref16L1 + refL1Stride));
            inPixel5 = _mm_loadu_si128((__m128i*) (ref16L1 + refL1Stride + 8));
            inPixel6 = _mm_loadu_si128((__m128i*) (ref16L1 + refL1Stride + 16));
            inPixel7 = _mm_loadu_si128((__m128i*) (ref16L1 + refL1Stride + 24));

            out8_2_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));
            out8_3_U8_L1 = _mm_packus_epi16(_mm_srli_epi16(inPixel6, 2), _mm_srli_epi16(inPixel7, 2));

            //AVG
            avg8_2_U8 = _mm_avg_epu8(out8_2_U8_L0, out8_2_U8_L1);
            avg8_3_U8 = _mm_avg_epu8(out8_3_U8_L0, out8_3_U8_L1);

            _mm_storeu_si128((__m128i*)(dstPtr + dst_stride), avg8_2_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + dst_stride + 16), avg8_3_U8);

            dstPtr += 2 * dst_stride;
            ref16L0 += 2 * refL0Stride;
            ref16L1 += 2 * refL1Stride;

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

            inPixel0 = _mm_loadu_si128((__m128i*) ref16L0);
            inPixel1 = _mm_loadu_si128((__m128i*)(ref16L0 + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(ref16L0 + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(ref16L0 + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(ref16L0 + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(ref16L0 + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(ref16L0 + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(ref16L0 + 56));



            out8_0_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel0, 2), _mm_srli_epi16(inPixel1, 2));
            out8_1_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel2, 2), _mm_srli_epi16(inPixel3, 2));
            out8_2_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel4, 2), _mm_srli_epi16(inPixel5, 2));
            out8_3_U8_L0 = _mm_packus_epi16(_mm_srli_epi16(inPixel6, 2), _mm_srli_epi16(inPixel7, 2));


            //List1

            inPixel0 = _mm_loadu_si128((__m128i*) ref16L1);
            inPixel1 = _mm_loadu_si128((__m128i*)(ref16L1 + 8));
            inPixel2 = _mm_loadu_si128((__m128i*)(ref16L1 + 16));
            inPixel3 = _mm_loadu_si128((__m128i*)(ref16L1 + 24));
            inPixel4 = _mm_loadu_si128((__m128i*)(ref16L1 + 32));
            inPixel5 = _mm_loadu_si128((__m128i*)(ref16L1 + 40));
            inPixel6 = _mm_loadu_si128((__m128i*)(ref16L1 + 48));
            inPixel7 = _mm_loadu_si128((__m128i*)(ref16L1 + 56));


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

            _mm_storeu_si128((__m128i*) dstPtr, avg8_0_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + 16), avg8_1_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + 32), avg8_2_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + 48), avg8_3_U8);


            dstPtr += dst_stride;
            ref16L0 += refL0Stride;
            ref16L1 += refL1Stride;
        }
    }


    return;
}
/********************************************************************************************************************
EB_ENC_msbPack2D_SSE2_INTRIN
*********************************************************************************************************************/
void EB_ENC_msbPack2D_SSE2_INTRIN(
    uint8_t     *in8BitBuffer,
    uint32_t     in8Stride,
    uint8_t     *innBitBuffer,
    uint16_t    *out16BitBuffer,
    uint32_t     innStride,
    uint32_t     outStride,
    uint32_t     width,
    uint32_t     height)
{
    uint32_t count_width, count_height;

    if (width == 4) {

        for (count_height = 0; count_height < height; count_height += 2) {
            _mm_storel_epi64((__m128i*)(out16BitBuffer), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(innBitBuffer)),
                _mm_cvtsi32_si128(*(uint32_t*)(in8BitBuffer))), 6));
            _mm_storel_epi64((__m128i*)(out16BitBuffer + outStride), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(innBitBuffer + innStride)),
                _mm_cvtsi32_si128(*(uint32_t*)(in8BitBuffer + in8Stride))), 6));
            out16BitBuffer += (outStride << 1);
            in8BitBuffer += (in8Stride << 1);
            innBitBuffer += (innStride << 1);
        }
    }
    else if (width == 8) {

        for (count_height = 0; count_height < height; count_height += 2) {

            _mm_storeu_si128((__m128i*)(out16BitBuffer), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(innBitBuffer)),
                _mm_loadl_epi64((__m128i*)(in8BitBuffer))), 6));
            _mm_storeu_si128((__m128i*)(out16BitBuffer + outStride), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(innBitBuffer + innStride)),
                _mm_loadl_epi64((__m128i*)(in8BitBuffer + in8Stride))), 6));
            out16BitBuffer += (outStride << 1);
            in8BitBuffer += (in8Stride << 1);
            innBitBuffer += (innStride << 1);
        }
    }
    else if (width == 16) {
        __m128i outPixel1, outPixel2, outPixel3, outPixel4, innBitBuffer_lo, innBitBuffer_hi, in8BitBuffer_lo, in8BitBuffer_hi;

        for (count_height = 0; count_height < height; count_height += 2) {

            innBitBuffer_lo = _mm_loadu_si128((__m128i *)innBitBuffer);
            innBitBuffer_hi = _mm_loadu_si128((__m128i *)(innBitBuffer + innStride));
            in8BitBuffer_lo = _mm_loadu_si128((__m128i *)in8BitBuffer);
            in8BitBuffer_hi = _mm_loadu_si128((__m128i *)(in8BitBuffer + in8Stride));

            outPixel1 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer_lo, in8BitBuffer_lo), 6);
            outPixel2 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer_lo, in8BitBuffer_lo), 6);
            outPixel3 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer_hi, in8BitBuffer_hi), 6);
            outPixel4 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer_hi, in8BitBuffer_hi), 6);

            _mm_storeu_si128((__m128i*)out16BitBuffer, outPixel1);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 8), outPixel2);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + outStride), outPixel3);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + outStride + 8), outPixel4);

            in8BitBuffer += (in8Stride << 1);
            innBitBuffer += (innStride << 1);
            out16BitBuffer += (outStride << 1);
        }
    }
    else if (width == 32) {
        __m128i innBitBuffer1, innBitBuffer2, innBitBuffer3, innBitBuffer4, in8BitBuffer1, in8BitBuffer2, in8BitBuffer3, in8BitBuffer4;
        __m128i outPixel1, outPixel2, outPixel3, outPixel4, outPixel5, outPixel6, outPixel7, outPixel8;

        for (count_height = 0; count_height < height; count_height += 2)
        {
            innBitBuffer1 = _mm_loadu_si128((__m128i *)innBitBuffer);
            innBitBuffer2 = _mm_loadu_si128((__m128i *)(innBitBuffer + 16));
            innBitBuffer3 = _mm_loadu_si128((__m128i *)(innBitBuffer + innStride));
            innBitBuffer4 = _mm_loadu_si128((__m128i *)(innBitBuffer + innStride + 16));

            in8BitBuffer1 = _mm_loadu_si128((__m128i *)in8BitBuffer);
            in8BitBuffer2 = _mm_loadu_si128((__m128i *)(in8BitBuffer + 16));
            in8BitBuffer3 = _mm_loadu_si128((__m128i *)(in8BitBuffer + in8Stride));
            in8BitBuffer4 = _mm_loadu_si128((__m128i *)(in8BitBuffer + in8Stride + 16));

            outPixel1 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer1, in8BitBuffer1), 6);
            outPixel2 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer1, in8BitBuffer1), 6);
            outPixel3 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer2, in8BitBuffer2), 6);
            outPixel4 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer2, in8BitBuffer2), 6);
            outPixel5 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer3, in8BitBuffer3), 6);
            outPixel6 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer3, in8BitBuffer3), 6);
            outPixel7 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer4, in8BitBuffer4), 6);
            outPixel8 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer4, in8BitBuffer4), 6);

            _mm_storeu_si128((__m128i*)out16BitBuffer, outPixel1);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 8), outPixel2);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 16), outPixel3);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 24), outPixel4);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + outStride), outPixel5);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + outStride + 8), outPixel6);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + outStride + 16), outPixel7);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + outStride + 24), outPixel8);

            in8BitBuffer += (in8Stride << 1);
            innBitBuffer += (innStride << 1);
            out16BitBuffer += (outStride << 1);
        }
    }
    else if (width == 64) {

        __m128i innBitBuffer1, innBitBuffer2, innBitBuffer3, innBitBuffer4, in8BitBuffer1, in8BitBuffer2, in8BitBuffer3, in8BitBuffer4;
        __m128i outPixel1, outPixel2, outPixel3, outPixel4, outPixel5, outPixel6, outPixel7, outPixel8;

        for (count_height = 0; count_height < height; ++count_height)
        {
            innBitBuffer1 = _mm_loadu_si128((__m128i *)innBitBuffer);
            innBitBuffer2 = _mm_loadu_si128((__m128i *)(innBitBuffer + 16));
            innBitBuffer3 = _mm_loadu_si128((__m128i *)(innBitBuffer + 32));
            innBitBuffer4 = _mm_loadu_si128((__m128i *)(innBitBuffer + 48));

            in8BitBuffer1 = _mm_loadu_si128((__m128i *)in8BitBuffer);
            in8BitBuffer2 = _mm_loadu_si128((__m128i *)(in8BitBuffer + 16));
            in8BitBuffer3 = _mm_loadu_si128((__m128i *)(in8BitBuffer + 32));
            in8BitBuffer4 = _mm_loadu_si128((__m128i *)(in8BitBuffer + 48));

            outPixel1 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer1, in8BitBuffer1), 6);
            outPixel2 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer1, in8BitBuffer1), 6);
            outPixel3 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer2, in8BitBuffer2), 6);
            outPixel4 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer2, in8BitBuffer2), 6);
            outPixel5 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer3, in8BitBuffer3), 6);
            outPixel6 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer3, in8BitBuffer3), 6);
            outPixel7 = _mm_srli_epi16(_mm_unpacklo_epi8(innBitBuffer4, in8BitBuffer4), 6);
            outPixel8 = _mm_srli_epi16(_mm_unpackhi_epi8(innBitBuffer4, in8BitBuffer4), 6);

            _mm_storeu_si128((__m128i*)out16BitBuffer, outPixel1);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 8), outPixel2);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 16), outPixel3);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 24), outPixel4);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 32), outPixel5);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 40), outPixel6);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 48), outPixel7);
            _mm_storeu_si128((__m128i*)(out16BitBuffer + 56), outPixel8);

            in8BitBuffer += in8Stride;
            innBitBuffer += innStride;
            out16BitBuffer += outStride;
        }
    }
    else {
        uint32_t innStrideDiff = (innStride << 1) - width;
        uint32_t in8StrideDiff = (in8Stride << 1) - width;
        uint32_t outStrideDiff = (outStride << 1) - width;

        if (!(width & 7)) {
            for (count_height = 0; count_height < height; count_height += 2) {
                for (count_width = 0; count_width < width; count_width += 8) {
                    _mm_storeu_si128((__m128i*)(out16BitBuffer), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(innBitBuffer)),
                        _mm_loadl_epi64((__m128i*)(in8BitBuffer))), 6));
                    _mm_storeu_si128((__m128i*)(out16BitBuffer + outStride), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(innBitBuffer + innStride)),
                        _mm_loadl_epi64((__m128i*)(in8BitBuffer + in8Stride))), 6));
                    out16BitBuffer += 8;
                    in8BitBuffer += 8;
                    innBitBuffer += 8;
                }
                in8BitBuffer += in8StrideDiff;
                innBitBuffer += innStrideDiff;
                out16BitBuffer += outStrideDiff;
            }
        }
        else {
            for (count_height = 0; count_height < height; count_height += 2) {
                for (count_width = 0; count_width < width; count_width += 4) {
                    _mm_storel_epi64((__m128i*)(out16BitBuffer), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(innBitBuffer)),
                        _mm_cvtsi32_si128(*(uint32_t*)(in8BitBuffer))), 6));
                    _mm_storel_epi64((__m128i*)(out16BitBuffer + outStride), _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(innBitBuffer + innStride)),
                        _mm_cvtsi32_si128(*(uint32_t*)(in8BitBuffer + in8Stride))), 6));
                    out16BitBuffer += 4;
                    in8BitBuffer += 4;
                    innBitBuffer += 4;
                }
                in8BitBuffer += in8StrideDiff;
                innBitBuffer += innStrideDiff;
                out16BitBuffer += outStrideDiff;
            }
        }
    }
}