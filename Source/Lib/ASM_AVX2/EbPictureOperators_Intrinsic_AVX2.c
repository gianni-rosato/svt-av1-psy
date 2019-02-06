/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include <stdio.h>
#include <immintrin.h>
#include "EbPictureOperators_AVX2.h"
#include "EbMemory_AVX2.h"

#define _mm256_set_m128i(/* __m128i */ hi, /* __m128i */ lo) \
    _mm256_insertf128_si256(_mm256_castsi128_si256(lo), (hi), 0x1)

void CompressedPackmsb_AVX2_INTRIN(
    uint8_t     *in8BitBuffer,
    uint32_t     in8Stride,
    uint8_t     *innBitBuffer,
    uint16_t    *out16BitBuffer,
    uint32_t     innStride,
    uint32_t     outStride,
    uint32_t     width,
    uint32_t     height)
{

    uint32_t y;

    if (width == 32)
    {
        __m256i inNBit, in8Bit, inNBitStride, in8BitStride, concat0, concat1, concat2, concat3;
        __m256i out0_15, out16_31, out_s0_s15, out_s16_s31;


        __m128i in2Bit, ext0, ext1, ext2, ext3, ext01, ext23, ext01h, ext23h, ext0_15, ext16_31, ext32_47, ext48_63;
        __m128i msk0;

        msk0 = _mm_set1_epi8((int8_t)0xC0);//1100.000

        //processing 2 lines for chroma
        for (y = 0; y < height; y += 2)
        {


            in2Bit = _mm_loadu_si128((__m128i*)innBitBuffer); //2 Lines Stored in 1D format-Could be replaced by 2 _mm_loadl_epi64

            ext0 = _mm_and_si128(in2Bit, msk0);
            ext1 = _mm_and_si128(_mm_slli_epi16(in2Bit, 2), msk0);
            ext2 = _mm_and_si128(_mm_slli_epi16(in2Bit, 4), msk0);
            ext3 = _mm_and_si128(_mm_slli_epi16(in2Bit, 6), msk0);

            ext01 = _mm_unpacklo_epi8(ext0, ext1);
            ext23 = _mm_unpacklo_epi8(ext2, ext3);
            ext0_15 = _mm_unpacklo_epi16(ext01, ext23);
            ext16_31 = _mm_unpackhi_epi16(ext01, ext23);

            ext01h = _mm_unpackhi_epi8(ext0, ext1);
            ext23h = _mm_unpackhi_epi8(ext2, ext3);
            ext32_47 = _mm_unpacklo_epi16(ext01h, ext23h);
            ext48_63 = _mm_unpackhi_epi16(ext01h, ext23h);

            inNBit = _mm256_set_m128i(ext16_31, ext0_15);
            inNBitStride = _mm256_set_m128i(ext48_63, ext32_47);

            in8Bit = _mm256_loadu_si256((__m256i*)in8BitBuffer);
            in8BitStride = _mm256_loadu_si256((__m256i*)(in8BitBuffer + in8Stride));


            //(outPixel | nBitPixel) concatenation is done with unpacklo_epi8 and unpackhi_epi8
            concat0 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBit, in8Bit), 6);
            concat1 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBit, in8Bit), 6);
            concat2 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBitStride, in8BitStride), 6);
            concat3 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBitStride, in8BitStride), 6);

            //Re-organize the packing for writing to the out buffer
            out0_15 = _mm256_inserti128_si256(concat0, _mm256_extracti128_si256(concat1, 0), 1);
            out16_31 = _mm256_inserti128_si256(concat1, _mm256_extracti128_si256(concat0, 1), 0);
            out_s0_s15 = _mm256_inserti128_si256(concat2, _mm256_extracti128_si256(concat3, 0), 1);
            out_s16_s31 = _mm256_inserti128_si256(concat3, _mm256_extracti128_si256(concat2, 1), 0);

            _mm256_store_si256((__m256i*) out16BitBuffer, out0_15);
            _mm256_store_si256((__m256i*) (out16BitBuffer + 16), out16_31);
            _mm256_store_si256((__m256i*) (out16BitBuffer + outStride), out_s0_s15);
            _mm256_store_si256((__m256i*) (out16BitBuffer + outStride + 16), out_s16_s31);

            in8BitBuffer += in8Stride << 1;
            innBitBuffer += innStride << 1;
            out16BitBuffer += outStride << 1;
        }
    }
    else if (width == 64)
    {
        __m256i inNBit, in8Bit, inNBit32, in8Bit32;
        __m256i concat0, concat1, concat2, concat3;
        __m256i out_0_15, out16_31, out32_47, out_48_63;


        __m128i in2Bit, ext0, ext1, ext2, ext3, ext01, ext23, ext01h, ext23h, ext0_15, ext16_31, ext32_47, ext48_63;
        __m128i msk;

        msk = _mm_set1_epi8((int8_t)0xC0);//1100.000

        //One row per iter
        for (y = 0; y < height; y++)
        {

            in2Bit = _mm_loadu_si128((__m128i*)innBitBuffer);

            ext0 = _mm_and_si128(in2Bit, msk);
            ext1 = _mm_and_si128(_mm_slli_epi16(in2Bit, 2), msk);
            ext2 = _mm_and_si128(_mm_slli_epi16(in2Bit, 4), msk);
            ext3 = _mm_and_si128(_mm_slli_epi16(in2Bit, 6), msk);

            ext01 = _mm_unpacklo_epi8(ext0, ext1);
            ext23 = _mm_unpacklo_epi8(ext2, ext3);
            ext0_15 = _mm_unpacklo_epi16(ext01, ext23);
            ext16_31 = _mm_unpackhi_epi16(ext01, ext23);

            ext01h = _mm_unpackhi_epi8(ext0, ext1);
            ext23h = _mm_unpackhi_epi8(ext2, ext3);
            ext32_47 = _mm_unpacklo_epi16(ext01h, ext23h);
            ext48_63 = _mm_unpackhi_epi16(ext01h, ext23h);

            inNBit = _mm256_set_m128i(ext16_31, ext0_15);
            inNBit32 = _mm256_set_m128i(ext48_63, ext32_47);

            in8Bit = _mm256_loadu_si256((__m256i*)in8BitBuffer);
            in8Bit32 = _mm256_loadu_si256((__m256i*)(in8BitBuffer + 32));

            //(outPixel | nBitPixel) concatenation
            concat0 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBit, in8Bit), 6);
            concat1 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBit, in8Bit), 6);
            concat2 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBit32, in8Bit32), 6);
            concat3 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBit32, in8Bit32), 6);

            //Re-organize the packing for writing to the out buffer
            out_0_15 = _mm256_inserti128_si256(concat0, _mm256_extracti128_si256(concat1, 0), 1);
            out16_31 = _mm256_inserti128_si256(concat1, _mm256_extracti128_si256(concat0, 1), 0);
            out32_47 = _mm256_inserti128_si256(concat2, _mm256_extracti128_si256(concat3, 0), 1);
            out_48_63 = _mm256_inserti128_si256(concat3, _mm256_extracti128_si256(concat2, 1), 0);

            _mm256_store_si256((__m256i*) out16BitBuffer, out_0_15);
            _mm256_store_si256((__m256i*) (out16BitBuffer + 16), out16_31);
            _mm256_store_si256((__m256i*) (out16BitBuffer + 32), out32_47);
            _mm256_store_si256((__m256i*) (out16BitBuffer + 48), out_48_63);

            in8BitBuffer += in8Stride;
            innBitBuffer += innStride;
            out16BitBuffer += outStride;

        }

    }

}
#if defined(_MSC_VER)
//nclude <intrin.h>
#endif

void CPack_AVX2_INTRIN(
    const uint8_t     *innBitBuffer,
    uint32_t     innStride,
    uint8_t     *inCompnBitBuffer,
    uint32_t     outStride,
    uint8_t    *localCache,

    uint32_t     width,
    uint32_t     height)
{

    uint32_t y;

    if (width == 32)
    {
        __m256i inNBit;

        __m256i ext0, ext1, ext2, ext3, ext0123, ext0123n, extp;
        __m256i msk0, msk1, msk2, msk3;

        msk0 = _mm256_set1_epi32(0x000000C0);//1100.0000
        msk1 = _mm256_set1_epi32(0x00000030);//0011.0000
        msk2 = _mm256_set1_epi32(0x0000000C);//0000.1100
        msk3 = _mm256_set1_epi32(0x00000003);//0000.0011

        //One row per iter
        for (y = 0; y < height; y++)
        {


            inNBit = _mm256_loadu_si256((__m256i*)innBitBuffer);

            ext0 = _mm256_and_si256(inNBit, msk0);
            ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
            ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
            ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

            ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

            ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));


            extp = _mm256_packus_epi32(ext0123, ext0123n);
            extp = _mm256_packus_epi16(extp, extp);

            _mm_storel_epi64((__m128i*) inCompnBitBuffer, _mm256_castsi256_si128(extp));
            inCompnBitBuffer += 8;
            innBitBuffer += innStride;


        }

    }
    else if (width == 64)
    {
        __m256i inNBit;
        __m256i ext0, ext1, ext2, ext3, ext0123, ext0123n, extp, extp1;
        __m256i msk0, msk1, msk2, msk3;

        msk0 = _mm256_set1_epi32(0x000000C0);//1100.0000
        msk1 = _mm256_set1_epi32(0x00000030);//0011.0000
        msk2 = _mm256_set1_epi32(0x0000000C);//0000.1100
        msk3 = _mm256_set1_epi32(0x00000003);//0000.0011
        if (height == 64)
        {

            uint8_t* localPtr = localCache;


            for (y = 0; y < height; y++)
            {


                inNBit = _mm256_loadu_si256((__m256i*)innBitBuffer);


                ext0 = _mm256_and_si256(inNBit, msk0);
                ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
                ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
                ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

                ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

                ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));

                extp = _mm256_packus_epi32(ext0123, ext0123n);
                extp = _mm256_packus_epi16(extp, extp);


                inNBit = _mm256_loadu_si256((__m256i*)(innBitBuffer + 32));

                ext0 = _mm256_and_si256(inNBit, msk0);
                ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
                ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
                ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

                ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

                ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));


                extp1 = _mm256_packus_epi32(ext0123, ext0123n);
                extp1 = _mm256_packus_epi16(extp1, extp1);

                extp = _mm256_unpacklo_epi64(extp, extp1);

                _mm_storeu_si128((__m128i*)  (localPtr + 16 * (y & 3)), _mm256_castsi256_si128(extp));

                if ((y & 3) == 3)
                {
                    __m256i c0 = _mm256_loadu_si256((__m256i*)(localPtr));
                    __m256i c1 = _mm256_loadu_si256((__m256i*)(localPtr + 32));
                    _mm256_stream_si256((__m256i*)&inCompnBitBuffer[0], c0);
                    _mm256_stream_si256((__m256i*)&inCompnBitBuffer[32], c1);
                    inCompnBitBuffer += 4 * outStride;
                }

                innBitBuffer += innStride;

            }

        }
        else {

            //One row per iter
            for (y = 0; y < height; y++)
            {


                inNBit = _mm256_loadu_si256((__m256i*)innBitBuffer);


                ext0 = _mm256_and_si256(inNBit, msk0);
                ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
                ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
                ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

                ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

                ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));

                extp = _mm256_packus_epi32(ext0123, ext0123n);
                extp = _mm256_packus_epi16(extp, extp);


                inNBit = _mm256_loadu_si256((__m256i*)(innBitBuffer + 32));

                ext0 = _mm256_and_si256(inNBit, msk0);
                ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
                ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
                ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

                ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

                ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));


                extp1 = _mm256_packus_epi32(ext0123, ext0123n);
                extp1 = _mm256_packus_epi16(extp1, extp1);

                extp = _mm256_unpacklo_epi64(extp, extp1);

                _mm_storeu_si128((__m128i*)  inCompnBitBuffer, _mm256_castsi256_si128(extp));

                inCompnBitBuffer += outStride;

                innBitBuffer += innStride;

            }

        }

    }

}


void EB_ENC_msbPack2D_AVX2_INTRIN_AL(
    uint8_t     *in8BitBuffer,
    uint32_t     in8Stride,
    uint8_t     *innBitBuffer,
    uint16_t    *out16BitBuffer,
    uint32_t     innStride,
    uint32_t     outStride,
    uint32_t     width,
    uint32_t     height)
{
    //(outPixel | nBitPixel) concatenation is done with unpacklo_epi8 and unpackhi_epi8

    uint32_t y, x;

    __m128i out0, out1;


    if (width == 4)
    {
        for (y = 0; y < height; y += 2) {

            out0 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)innBitBuffer), _mm_cvtsi32_si128(*(uint32_t *)in8BitBuffer)), 6);
            out1 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(innBitBuffer + innStride)), _mm_cvtsi32_si128(*(uint32_t *)(in8BitBuffer + in8Stride))), 6);

            _mm_storel_epi64((__m128i*) out16BitBuffer, out0);
            _mm_storel_epi64((__m128i*) (out16BitBuffer + outStride), out1);

            in8BitBuffer += in8Stride << 1;
            innBitBuffer += innStride << 1;
            out16BitBuffer += outStride << 1;
        }
    }
    else if (width == 8)
    {
        for (y = 0; y < height; y += 2) {

            out0 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)innBitBuffer), _mm_loadl_epi64((__m128i*)in8BitBuffer)), 6);
            out1 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(innBitBuffer + innStride)), _mm_loadl_epi64((__m128i*)(in8BitBuffer + in8Stride))), 6);

            _mm_storeu_si128((__m128i*) out16BitBuffer, out0);
            _mm_storeu_si128((__m128i*) (out16BitBuffer + outStride), out1);

            in8BitBuffer += in8Stride << 1;
            innBitBuffer += innStride << 1;
            out16BitBuffer += outStride << 1;
        }
    }
    else if (width == 16)
    {
        __m128i inNBit, in8Bit, inNBitStride, in8BitStride, out0, out1, out2, out3;

        for (y = 0; y < height; y += 2) {
            inNBit = _mm_loadu_si128((__m128i*)innBitBuffer);
            in8Bit = _mm_loadu_si128((__m128i*)in8BitBuffer);
            inNBitStride = _mm_loadu_si128((__m128i*)(innBitBuffer + innStride));
            in8BitStride = _mm_loadu_si128((__m128i*)(in8BitBuffer + in8Stride));

            out0 = _mm_srli_epi16(_mm_unpacklo_epi8(inNBit, in8Bit), 6);
            out1 = _mm_srli_epi16(_mm_unpackhi_epi8(inNBit, in8Bit), 6);
            out2 = _mm_srli_epi16(_mm_unpacklo_epi8(inNBitStride, in8BitStride), 6);
            out3 = _mm_srli_epi16(_mm_unpackhi_epi8(inNBitStride, in8BitStride), 6);

            _mm_storeu_si128((__m128i*) out16BitBuffer, out0);
            _mm_storeu_si128((__m128i*) (out16BitBuffer + 8), out1);
            _mm_storeu_si128((__m128i*) (out16BitBuffer + outStride), out2);
            _mm_storeu_si128((__m128i*) (out16BitBuffer + outStride + 8), out3);

            in8BitBuffer += in8Stride << 1;
            innBitBuffer += innStride << 1;
            out16BitBuffer += outStride << 1;
        }
    }
    else if (width == 32)
    {
        __m256i inNBit, in8Bit, inNBitStride, in8BitStride, concat0, concat1, concat2, concat3;
        __m256i out0_15, out16_31, out_s0_s15, out_s16_s31;

        for (y = 0; y < height; y += 2) {

            inNBit = _mm256_loadu_si256((__m256i*)innBitBuffer);
            in8Bit = _mm256_loadu_si256((__m256i*)in8BitBuffer);
            inNBitStride = _mm256_loadu_si256((__m256i*)(innBitBuffer + innStride));
            in8BitStride = _mm256_loadu_si256((__m256i*)(in8BitBuffer + in8Stride));

            //(outPixel | nBitPixel) concatenation is done with unpacklo_epi8 and unpackhi_epi8
            concat0 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBit, in8Bit), 6);
            concat1 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBit, in8Bit), 6);
            concat2 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBitStride, in8BitStride), 6);
            concat3 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBitStride, in8BitStride), 6);

            //Re-organize the packing for writing to the out buffer
            out0_15 = _mm256_inserti128_si256(concat0, _mm256_extracti128_si256(concat1, 0), 1);
            out16_31 = _mm256_inserti128_si256(concat1, _mm256_extracti128_si256(concat0, 1), 0);
            out_s0_s15 = _mm256_inserti128_si256(concat2, _mm256_extracti128_si256(concat3, 0), 1);
            out_s16_s31 = _mm256_inserti128_si256(concat3, _mm256_extracti128_si256(concat2, 1), 0);

            _mm256_store_si256((__m256i*) out16BitBuffer, out0_15);
            _mm256_store_si256((__m256i*) (out16BitBuffer + 16), out16_31);
            _mm256_store_si256((__m256i*) (out16BitBuffer + outStride), out_s0_s15);
            _mm256_store_si256((__m256i*) (out16BitBuffer + outStride + 16), out_s16_s31);

            in8BitBuffer += in8Stride << 1;
            //innBitBuffer += innStride << 1;
            innBitBuffer += innStride * 2;
            out16BitBuffer += outStride << 1;
        }
    }
    else if (width == 64)
    {
        __m256i inNBit, in8Bit, inNBitStride, in8BitStride, inNBit32, in8Bit32, inNBitStride32, in8BitStride32;
        __m256i concat0, concat1, concat2, concat3, concat4, concat5, concat6, concat7;
        __m256i out_0_15, out16_31, out32_47, out_48_63, out_s0_s15, out_s16_s31, out_s32_s47, out_s48_s63;

        for (y = 0; y < height; y += 2) {

            inNBit = _mm256_loadu_si256((__m256i*)innBitBuffer);
            in8Bit = _mm256_loadu_si256((__m256i*)in8BitBuffer);
            inNBit32 = _mm256_loadu_si256((__m256i*)(innBitBuffer + 32));
            in8Bit32 = _mm256_loadu_si256((__m256i*)(in8BitBuffer + 32));
            inNBitStride = _mm256_loadu_si256((__m256i*)(innBitBuffer + innStride));
            in8BitStride = _mm256_loadu_si256((__m256i*)(in8BitBuffer + in8Stride));
            inNBitStride32 = _mm256_loadu_si256((__m256i*)(innBitBuffer + innStride + 32));
            in8BitStride32 = _mm256_loadu_si256((__m256i*)(in8BitBuffer + in8Stride + 32));
            //(outPixel | nBitPixel) concatenation is done with unpacklo_epi8 and unpackhi_epi8
            concat0 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBit, in8Bit), 6);
            concat1 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBit, in8Bit), 6);
            concat2 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBit32, in8Bit32), 6);
            concat3 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBit32, in8Bit32), 6);
            concat4 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBitStride, in8BitStride), 6);
            concat5 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBitStride, in8BitStride), 6);
            concat6 = _mm256_srli_epi16(_mm256_unpacklo_epi8(inNBitStride32, in8BitStride32), 6);
            concat7 = _mm256_srli_epi16(_mm256_unpackhi_epi8(inNBitStride32, in8BitStride32), 6);

            //Re-organize the packing for writing to the out buffer
            out_0_15 = _mm256_inserti128_si256(concat0, _mm256_extracti128_si256(concat1, 0), 1);
            out16_31 = _mm256_inserti128_si256(concat1, _mm256_extracti128_si256(concat0, 1), 0);
            out32_47 = _mm256_inserti128_si256(concat2, _mm256_extracti128_si256(concat3, 0), 1);
            out_48_63 = _mm256_inserti128_si256(concat3, _mm256_extracti128_si256(concat2, 1), 0);
            out_s0_s15 = _mm256_inserti128_si256(concat4, _mm256_extracti128_si256(concat5, 0), 1);
            out_s16_s31 = _mm256_inserti128_si256(concat5, _mm256_extracti128_si256(concat4, 1), 0);
            out_s32_s47 = _mm256_inserti128_si256(concat6, _mm256_extracti128_si256(concat7, 0), 1);
            out_s48_s63 = _mm256_inserti128_si256(concat7, _mm256_extracti128_si256(concat6, 1), 0);

            _mm256_store_si256((__m256i*) out16BitBuffer, out_0_15);
            _mm256_store_si256((__m256i*) (out16BitBuffer + 16), out16_31);
            _mm256_store_si256((__m256i*) (out16BitBuffer + 32), out32_47);
            _mm256_store_si256((__m256i*) (out16BitBuffer + 48), out_48_63);

            _mm256_store_si256((__m256i*) (out16BitBuffer + outStride), out_s0_s15);
            _mm256_store_si256((__m256i*) (out16BitBuffer + outStride + 16), out_s16_s31);
            _mm256_store_si256((__m256i*) (out16BitBuffer + outStride + 32), out_s32_s47);
            _mm256_store_si256((__m256i*) (out16BitBuffer + outStride + 48), out_s48_s63);

            in8BitBuffer += in8Stride << 1;
            //innBitBuffer += innStride << 1;
            innBitBuffer += innStride * 2;
            out16BitBuffer += outStride << 1;
        }
    }
    else
    {
        uint32_t innStrideDiff = 2 * innStride;
        uint32_t in8StrideDiff = 2 * in8Stride;
        uint32_t outStrideDiff = 2 * outStride;
        innStrideDiff -= width;
        in8StrideDiff -= width;
        outStrideDiff -= width;

        if (!(width & 7)) {

            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 8) {

                    out0 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)innBitBuffer), _mm_loadl_epi64((__m128i*)in8BitBuffer)), 6);
                    out1 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(innBitBuffer + innStride)), _mm_loadl_epi64((__m128i*)(in8BitBuffer + in8Stride))), 6);

                    _mm_storeu_si128((__m128i*) out16BitBuffer, out0);
                    _mm_storeu_si128((__m128i*) (out16BitBuffer + outStride), out1);

                    in8BitBuffer += 8;
                    innBitBuffer += 8;
                    out16BitBuffer += 8;
                }
                in8BitBuffer += in8StrideDiff;
                innBitBuffer += innStrideDiff;
                out16BitBuffer += outStrideDiff;
            }
        }
        else {
            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 4) {

                    out0 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)innBitBuffer), _mm_cvtsi32_si128(*(uint32_t *)in8BitBuffer)), 6);
                    out1 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(innBitBuffer + innStride)), _mm_cvtsi32_si128(*(uint32_t *)(in8BitBuffer + in8Stride))), 6);

                    _mm_storel_epi64((__m128i*) out16BitBuffer, out0);
                    _mm_storel_epi64((__m128i*) (out16BitBuffer + outStride), out1);

                    in8BitBuffer += 4;
                    innBitBuffer += 4;
                    out16BitBuffer += 4;
                }
                in8BitBuffer += in8StrideDiff;
                innBitBuffer += innStrideDiff;
                out16BitBuffer += outStrideDiff;
            }
        }
    }
}

#define ALSTORE  1
#define B256     1

void UnpackAvg_AVX2_INTRIN(
    uint16_t *ref16L0,
    uint32_t  refL0Stride,
    uint16_t *ref16L1,
    uint32_t  refL1Stride,
    uint8_t  *dstPtr,
    uint32_t  dst_stride,
    uint32_t  width,
    uint32_t  height)
{

    uint32_t   y;
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
#if ALSTORE
            _mm_store_si128((__m128i*) dstPtr, avg8_0_U8);
#else
            _mm_storeu_si128((__m128i*) dstPtr, avg8_0_U8);
#endif

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
#if ALSTORE
            _mm_store_si128((__m128i*)(dstPtr + dst_stride), avg8_2_U8);
#else
            _mm_storeu_si128((__m128i*)(dstPtr + dst_stride), avg8_2_U8);
#endif
            dstPtr += 2 * dst_stride;
            ref16L0 += 2 * refL0Stride;
            ref16L1 += 2 * refL1Stride;

        }

    }
    else if (width == 32)
    {

#if B256
        __m256i inVal16b_0, inVal16b_1;
        __m256i data8b_32_0_L0, data8b_32_0_L1;
        __m256i avg8b_32_0;
#else
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i out8_0_U8_L0, out8_1_U8_L0, out8_2_U8_L0, out8_3_U8_L0;
        __m128i out8_0_U8_L1, out8_1_U8_L1, out8_2_U8_L1, out8_3_U8_L1;
        __m128i avg8_0_U8, avg8_1_U8, avg8_2_U8, avg8_3_U8;
#endif

        for (y = 0; y < height; y += 2)
        {

#if B256
            //--------
            //Line One
            //--------

            //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L0 + 16));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L1 + 16));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);

            _mm256_storeu_si256((__m256i *)(dstPtr), avg8b_32_0);

            //--------
            //Line Two
            //--------
              //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*)(ref16L0 + refL0Stride));
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L0 + refL0Stride + 16));

            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));

            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*)(ref16L1 + refL1Stride));
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L1 + refL1Stride + 16));

            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);

            _mm256_storeu_si256((__m256i *)(dstPtr + dst_stride), avg8b_32_0);

#else
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
#if ALSTORE
            _mm_store_si128((__m128i*) dstPtr, avg8_0_U8);
            _mm_store_si128((__m128i*)(dstPtr + 16), avg8_1_U8);
#else
            _mm_storeu_si128((__m128i*) dstPtr, avg8_0_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + 16), avg8_1_U8);
#endif

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
#if ALSTORE
            _mm_store_si128((__m128i*)(dstPtr + dst_stride), avg8_2_U8);
            _mm_store_si128((__m128i*)(dstPtr + dst_stride + 16), avg8_3_U8);
#else
            _mm_storeu_si128((__m128i*)(dstPtr + dst_stride), avg8_2_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + dst_stride + 16), avg8_3_U8);
#endif

#endif
            dstPtr += 2 * dst_stride;
            ref16L0 += 2 * refL0Stride;
            ref16L1 += 2 * refL1Stride;

        }

    }
    else if (width == 64)
    {


#if B256
        __m256i inVal16b_0, inVal16b_1, inVal16b_2, inVal16b_3;
        __m256i data8b_32_0_L0, data8b_32_1_L0, data8b_32_0_L1, data8b_32_1_L1;
        __m256i avg8b_32_0, avg8b_32_1;
#else
        __m128i inPixel2, inPixel3, inPixel4, inPixel5, inPixel6, inPixel7;
        __m128i out8_0_U8_L0, out8_1_U8_L0, out8_2_U8_L0, out8_3_U8_L0;
        __m128i out8_0_U8_L1, out8_1_U8_L1, out8_2_U8_L1, out8_3_U8_L1;
        __m128i avg8_0_U8, avg8_1_U8, avg8_2_U8, avg8_3_U8;

#endif

        for (y = 0; y < height; ++y)
        {

#if B256       // _mm256_lddqu_si256

            //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L0 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16L0 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16L0 + 48));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L1 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16L1 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16L1 + 48));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);
            avg8b_32_1 = _mm256_avg_epu8(data8b_32_1_L0, data8b_32_1_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);
            avg8b_32_1 = _mm256_permute4x64_epi64(avg8b_32_1, 216);

            _mm256_storeu_si256((__m256i *)(dstPtr), avg8b_32_0);
            _mm256_storeu_si256((__m256i *)(dstPtr + 32), avg8b_32_1);
#else
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
#if ALSTORE
            _mm_store_si128((__m128i*) dstPtr, avg8_0_U8);
            _mm_store_si128((__m128i*)(dstPtr + 16), avg8_1_U8);
            _mm_store_si128((__m128i*)(dstPtr + 32), avg8_2_U8);
            _mm_store_si128((__m128i*)(dstPtr + 48), avg8_3_U8);
#else
            _mm_storeu_si128((__m128i*) dstPtr, avg8_0_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + 16), avg8_1_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + 32), avg8_2_U8);
            _mm_storeu_si128((__m128i*)(dstPtr + 48), avg8_3_U8);
#endif

#endif
            dstPtr += dst_stride;
            ref16L0 += refL0Stride;
            ref16L1 += refL1Stride;
        }
    }


    return;
}

int32_t  sumResidual8bit_AVX2_INTRIN(
    int16_t * inPtr,
    uint32_t   size,
    uint32_t   strideIn)
{

    int32_t  sumBlock;

    __m128i in0, in1, in01, in2, in3, in23, sum, sumL, sumH;
    __m256i sum0, sum1, sum2, sum3, sum0L, sum0H, sumT, sum01, sumTPerm;
    uint32_t rowIndex;

    //Assumption: 9bit or 11bit residual data . for bigger block sizes or bigger bit depths , re-asses the dynamic range of the internal calculation

    if (size == 4) { //SSSE3

        __m128i zer = _mm_setzero_si128();

        in0 = _mm_loadl_epi64((__m128i*)inPtr);
        in1 = _mm_loadl_epi64((__m128i*)(inPtr + strideIn));
        in1 = _mm_shuffle_epi32(in1, 0x4A); //01.00.10.10
        in01 = _mm_or_si128(in1, in0);

        in2 = _mm_loadl_epi64((__m128i*)(inPtr + 2 * strideIn));
        in3 = _mm_loadl_epi64((__m128i*)(inPtr + 3 * strideIn));
        in3 = _mm_shuffle_epi32(in3, 0x4A); //01.00.10.10
        in23 = _mm_or_si128(in3, in2);

        sum = _mm_add_epi16(in01, in23);
        sum = _mm_hadd_epi16(sum, zer);
        sum = _mm_hadd_epi16(sum, zer);
        sum = _mm_hadd_epi16(sum, zer);

        sum = _mm_cvtepi16_epi32(sum);
        sumBlock = _mm_cvtsi128_si32(sum);

        return sumBlock;

    }
    else if (size == 8) {//SSSE3

        __m128i zer = _mm_setzero_si128();

        sum = _mm_add_epi16(_mm_loadu_si128((__m128i*)(inPtr + 0 * strideIn)), _mm_loadu_si128((__m128i*)(inPtr + 1 * strideIn)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(inPtr + 2 * strideIn)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(inPtr + 3 * strideIn)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(inPtr + 4 * strideIn)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(inPtr + 5 * strideIn)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(inPtr + 6 * strideIn)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(inPtr + 7 * strideIn)));

        sum = _mm_hadd_epi16(sum, zer);
        sum = _mm_hadd_epi16(sum, zer);
        sum = _mm_hadd_epi16(sum, zer);

        sum = _mm_cvtepi16_epi32(sum); //the sum is on 16bit, for negative values, we need to extend the sign to the next 16bit, so that the next extraction to int32_t is fine.
        sumBlock = _mm_cvtsi128_si32(sum);

        return sumBlock;

    }
    else if (size == 16) {//AVX2

        sum0 = _mm256_add_epi16(_mm256_loadu_si256((__m256i *)(inPtr + 0 * strideIn)), _mm256_loadu_si256((__m256i *)(inPtr + 1 * strideIn)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(inPtr + 2 * strideIn)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(inPtr + 3 * strideIn)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(inPtr + 4 * strideIn)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(inPtr + 5 * strideIn)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(inPtr + 6 * strideIn)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(inPtr + 7 * strideIn)));

        inPtr += 8 * strideIn;
        sum1 = _mm256_add_epi16(_mm256_loadu_si256((__m256i *)(inPtr + 0 * strideIn)), _mm256_loadu_si256((__m256i *)(inPtr + 1 * strideIn)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(inPtr + 2 * strideIn)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(inPtr + 3 * strideIn)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(inPtr + 4 * strideIn)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(inPtr + 5 * strideIn)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(inPtr + 6 * strideIn)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(inPtr + 7 * strideIn)));

        sum01 = _mm256_add_epi16(sum0, sum1);

        //go from 16bit to 32bit (to support big values)
        sumL = _mm256_castsi256_si128(sum01);
        sumH = _mm256_extracti128_si256(sum01, 1);
        sum0L = _mm256_cvtepi16_epi32(sumL);
        sum0H = _mm256_cvtepi16_epi32(sumH);

        sumT = _mm256_add_epi32(sum0L, sum0H);

        sumT = _mm256_hadd_epi32(sumT, sumT);
        sumT = _mm256_hadd_epi32(sumT, sumT);
        sumTPerm = _mm256_permute4x64_epi64(sumT, 2); //00.00.00.10
        sumT = _mm256_add_epi32(sumT, sumTPerm);

        sum = _mm256_castsi256_si128(sumT);
        sumBlock = _mm_cvtsi128_si32(sum);

        return sumBlock;

    }
    else if (size == 32) {//AVX2
        int16_t *inPtrTemp = inPtr;

        sum0 = sum1 = sum2 = sum3 = _mm256_setzero_si256();
        for (rowIndex = 0; rowIndex < size; rowIndex += 2) { // Parse every two rows
            sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(inPtrTemp)));
            sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(inPtrTemp + 16)));
            inPtrTemp += strideIn;
            sum2 = _mm256_add_epi16(sum2, _mm256_loadu_si256((__m256i *)(inPtrTemp)));
            sum3 = _mm256_add_epi16(sum3, _mm256_loadu_si256((__m256i *)(inPtrTemp + 16)));
            inPtrTemp += strideIn;

        }
        //go from 16bit to 32bit (to support big values)
        sumL = _mm256_castsi256_si128(sum0);
        sumH = _mm256_extracti128_si256(sum0, 1);
        sum0L = _mm256_cvtepi16_epi32(sumL);
        sum0H = _mm256_cvtepi16_epi32(sumH);
        sumT = _mm256_add_epi32(sum0L, sum0H);

        sumL = _mm256_castsi256_si128(sum1);
        sumH = _mm256_extracti128_si256(sum1, 1);
        sum0L = _mm256_cvtepi16_epi32(sumL);
        sum0H = _mm256_cvtepi16_epi32(sumH);
        sumT = _mm256_add_epi32(sumT, sum0L);
        sumT = _mm256_add_epi32(sumT, sum0H);

        sumL = _mm256_castsi256_si128(sum2);
        sumH = _mm256_extracti128_si256(sum2, 1);
        sum0L = _mm256_cvtepi16_epi32(sumL);
        sum0H = _mm256_cvtepi16_epi32(sumH);
        sumT = _mm256_add_epi32(sumT, sum0L);
        sumT = _mm256_add_epi32(sumT, sum0H);

        sumL = _mm256_castsi256_si128(sum3);
        sumH = _mm256_extracti128_si256(sum3, 1);
        sum0L = _mm256_cvtepi16_epi32(sumL);
        sum0H = _mm256_cvtepi16_epi32(sumH);
        sumT = _mm256_add_epi32(sumT, sum0L);
        sumT = _mm256_add_epi32(sumT, sum0H);

        sumT = _mm256_hadd_epi32(sumT, sumT);
        sumT = _mm256_hadd_epi32(sumT, sumT);
        sumTPerm = _mm256_permute4x64_epi64(sumT, 2); //00.00.00.10
        sumT = _mm256_add_epi32(sumT, sumTPerm);

        sum = _mm256_castsi256_si128(sumT);
        sumBlock = _mm_cvtsi128_si32(sum);

        return sumBlock;
    }

    else {
        printf("\n add the rest \n");
        return 0;
    }


}

void memset16bitBlock_AVX2_INTRIN(
    int16_t * inPtr,
    uint32_t   strideIn,
    uint32_t   size,
    int16_t   value
)
{


    if (size == 4) {

        __m128i line = _mm_set1_epi16(value);

        _mm_storel_epi64((__m128i *)(inPtr + 0 * strideIn), line);
        _mm_storel_epi64((__m128i *)(inPtr + 1 * strideIn), line);
        _mm_storel_epi64((__m128i *)(inPtr + 2 * strideIn), line);
        _mm_storel_epi64((__m128i *)(inPtr + 3 * strideIn), line);

    }
    else if (size == 8) {

        __m128i line = _mm_set1_epi16(value);

        _mm_storeu_si128((__m128i *)(inPtr + 0 * strideIn), line);
        _mm_storeu_si128((__m128i *)(inPtr + 1 * strideIn), line);
        _mm_storeu_si128((__m128i *)(inPtr + 2 * strideIn), line);
        _mm_storeu_si128((__m128i *)(inPtr + 3 * strideIn), line);
        _mm_storeu_si128((__m128i *)(inPtr + 4 * strideIn), line);
        _mm_storeu_si128((__m128i *)(inPtr + 5 * strideIn), line);
        _mm_storeu_si128((__m128i *)(inPtr + 6 * strideIn), line);
        _mm_storeu_si128((__m128i *)(inPtr + 7 * strideIn), line);

    }
    else if (size == 16) {

        __m256i line = _mm256_set1_epi16(value);

        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn), line);

        inPtr += 8 * strideIn;

        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn), line);


    }
    else if (size == 32) {

        __m256i line = _mm256_set1_epi16(value);

        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn + 16), line);

        inPtr += 8 * strideIn;

        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn + 16), line);

        inPtr += 8 * strideIn;

        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn + 16), line);

        inPtr += 8 * strideIn;

        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 0 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 1 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 2 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 3 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 4 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 5 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 6 * strideIn + 16), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn), line);
        _mm256_storeu_si256((__m256i *)(inPtr + 7 * strideIn + 16), line);

    }


    else {
        printf("\n add the rest \n");
    }

}


void UnpackAvgSafeSub_AVX2_INTRIN(
    uint16_t *ref16L0,
    uint32_t  refL0Stride,
    uint16_t *ref16L1,
    uint32_t  refL1Stride,
    uint8_t  *dstPtr,
    uint32_t  dst_stride,
    EbBool  subPred,
    uint32_t  width,
    uint32_t  height)
{

    uint32_t   y;
    __m128i inPixel0, inPixel1;


    if (width == 8)
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

        if (subPred) {
            ref16L0 -= (refL0Stride >> 1);
            ref16L1 -= (refL1Stride >> 1);
            dstPtr -= (dst_stride >> 1);
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

            _mm_store_si128((__m128i*) dstPtr, avg8_0_U8);

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

            _mm_store_si128((__m128i*)(dstPtr + dst_stride), avg8_2_U8);

            dstPtr += 2 * dst_stride;
            ref16L0 += 2 * refL0Stride;
            ref16L1 += 2 * refL1Stride;

        }

        if (subPred) {
            ref16L0 -= (refL0Stride >> 1);
            ref16L1 -= (refL1Stride >> 1);
            dstPtr -= (dst_stride >> 1);
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
            _mm_store_si128((__m128i*) dstPtr, avg8_0_U8);
        }
    }
    else if (width == 32)
    {

        __m256i inVal16b_0, inVal16b_1;
        __m256i data8b_32_0_L0, data8b_32_0_L1;
        __m256i avg8b_32_0;

        for (y = 0; y < height; y += 2)
        {

            //--------
            //Line One
            //--------

            //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L0 + 16));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L1 + 16));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);

            _mm256_storeu_si256((__m256i *)(dstPtr), avg8b_32_0);

            //--------
            //Line Two
            //--------
              //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*)(ref16L0 + refL0Stride));
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L0 + refL0Stride + 16));

            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));

            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*)(ref16L1 + refL1Stride));
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L1 + refL1Stride + 16));

            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);

            _mm256_storeu_si256((__m256i *)(dstPtr + dst_stride), avg8b_32_0);


            dstPtr += 2 * dst_stride;
            ref16L0 += 2 * refL0Stride;
            ref16L1 += 2 * refL1Stride;

        }

        if (subPred) {
            ref16L0 -= (refL0Stride >> 1);
            ref16L1 -= (refL1Stride >> 1);
            dstPtr -= (dst_stride >> 1);
            //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L0 + 16));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L1 + 16));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);
            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);
            _mm256_storeu_si256((__m256i *)(dstPtr), avg8b_32_0);

        }

    }
    else if (width == 64)
    {
        __m256i inVal16b_0, inVal16b_1, inVal16b_2, inVal16b_3;
        __m256i data8b_32_0_L0, data8b_32_1_L0, data8b_32_0_L1, data8b_32_1_L1;
        __m256i avg8b_32_0, avg8b_32_1;


        for (y = 0; y < height; ++y)
        {


            //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L0 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16L0 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16L0 + 48));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L1 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16L1 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16L1 + 48));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);
            avg8b_32_1 = _mm256_avg_epu8(data8b_32_1_L0, data8b_32_1_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);
            avg8b_32_1 = _mm256_permute4x64_epi64(avg8b_32_1, 216);

            _mm256_storeu_si256((__m256i *)(dstPtr), avg8b_32_0);
            _mm256_storeu_si256((__m256i *)(dstPtr + 32), avg8b_32_1);

            dstPtr += dst_stride;
            ref16L0 += refL0Stride;
            ref16L1 += refL1Stride;
        }

        if (subPred) {
            ref16L0 -= (refL0Stride >> 1);
            ref16L1 -= (refL1Stride >> 1);
            dstPtr -= (dst_stride >> 1);
            //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L0 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16L0 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16L0 + 48));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16L1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16L1 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16L1 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16L1 + 48));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);
            avg8b_32_1 = _mm256_avg_epu8(data8b_32_1_L0, data8b_32_1_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);
            avg8b_32_1 = _mm256_permute4x64_epi64(avg8b_32_1, 216);

            _mm256_storeu_si256((__m256i *)(dstPtr), avg8b_32_0);
            _mm256_storeu_si256((__m256i *)(dstPtr + 32), avg8b_32_1);
        }


    }


    return;
}


void PictureAdditionKernel4x4_AV1_SSE2_INTRIN(
    uint8_t  *predPtr,
    uint32_t  predStride,
    int32_t *residual_ptr,
    uint32_t  residualStride,
    uint8_t  *reconPtr,
    uint32_t  reconStride,
    uint32_t  width,
    uint32_t  height,
    int32_t     bd)
{

    __m128i predReg, xmm0, recon_0_7, resReg;
    uint32_t y;
    xmm0 = _mm_setzero_si128();

    for (y = 0; y < 4; ++y) {
        predReg = _mm_cvtsi32_si128(*(uint32_t *)predPtr);
        predReg = _mm_unpacklo_epi8(predReg, xmm0);
        predReg = _mm_unpacklo_epi16(predReg, xmm0);
        resReg = _mm_loadu_si128((__m128i *)residual_ptr);
        resReg = _mm_add_epi32(resReg, predReg);
        recon_0_7 = _mm_packus_epi32(resReg, xmm0);
        recon_0_7 = _mm_packus_epi16(recon_0_7, xmm0);
        *(uint32_t *)reconPtr = _mm_cvtsi128_si32(recon_0_7);
        predPtr += predStride;
        residual_ptr += residualStride;
        reconPtr += reconStride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}
void PictureAdditionKernel8x8_AV1_SSE2_INTRIN(
    uint8_t  *predPtr,
    uint32_t  predStride,
    int32_t *residual_ptr,
    uint32_t  residualStride,
    uint8_t  *reconPtr,
    uint32_t  reconStride,
    uint32_t  width,
    uint32_t  height,
    int32_t     bd)
{

    __m256i predReg, resReg, recon_0_7, xmm0;
    __m128i predReg_128, predReg_128Lo, predReg_128Hi, xmm0_128, recon_0_7_128;
    uint32_t y;
    xmm0_128 = _mm_setzero_si128();
    xmm0 = _mm256_setzero_si256();

    for (y = 0; y < 8; ++y) {
        predReg_128 = _mm_cvtsi64_si128(*(uint64_t *)predPtr);
        predReg_128 = _mm_unpacklo_epi8(predReg_128, xmm0_128);
        predReg_128Lo = _mm_unpacklo_epi16(predReg_128, xmm0_128);
        predReg_128Hi = _mm_unpackhi_epi16(predReg_128, xmm0_128);
        predReg = _mm256_set_m128i(predReg_128Hi, predReg_128Lo);
        resReg = _mm256_loadu_si256((__m256i*)residual_ptr);
        resReg = _mm256_add_epi32(predReg, resReg);
        recon_0_7 = _mm256_packus_epi32(resReg, xmm0);
        recon_0_7_128 = _mm_slli_si128(_mm256_extracti128_si256(recon_0_7, 1), 8);
        recon_0_7_128 = _mm_or_si128(_mm256_castsi256_si128(recon_0_7), recon_0_7_128);
        recon_0_7_128 = _mm_packus_epi16(recon_0_7_128, xmm0_128);
        *(uint64_t *)reconPtr = _mm_cvtsi128_si64(recon_0_7_128);
        predPtr += predStride;
        residual_ptr += residualStride;
        reconPtr += reconStride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}
void PictureAdditionKernel16x16_AV1_SSE2_INTRIN(
    uint8_t  *predPtr,
    uint32_t  predStride,
    int32_t *residual_ptr,
    uint32_t  residualStride,
    uint8_t  *reconPtr,
    uint32_t  reconStride,
    uint32_t  width,
    uint32_t  height,
    int32_t     bd)
{

    __m256i resReg, recon_0_7, xmm0, predRegLo, predRegHi, resRegLo, resRegHi;
    __m128i predReg_128, predReg_128Lo, predReg_128Hi, xmm0_128, predReg_128Lo16Lo, predReg_128Lo16Hi, predReg_128Hi16Lo, predReg_128Hi16Hi;
    uint32_t y;
    xmm0_128 = _mm_setzero_si128();
    xmm0 = _mm256_setzero_si256();

    for (y = 0; y < 16; ++y) {

        predReg_128 = _mm_loadu_si128((__m128i *)predPtr);
        predReg_128Lo = _mm_unpacklo_epi8(predReg_128, xmm0_128);
        predReg_128Hi = _mm_unpackhi_epi8(predReg_128, xmm0_128);
        predReg_128Lo16Lo = _mm_unpacklo_epi16(predReg_128Lo, xmm0_128);
        predReg_128Lo16Hi = _mm_unpackhi_epi16(predReg_128Lo, xmm0_128);
        predReg_128Hi16Lo = _mm_unpacklo_epi16(predReg_128Hi, xmm0_128);
        predReg_128Hi16Hi = _mm_unpackhi_epi16(predReg_128Hi, xmm0_128);
        predRegLo = _mm256_set_m128i(predReg_128Lo16Hi, predReg_128Lo16Lo);
        predRegHi = _mm256_set_m128i(predReg_128Hi16Hi, predReg_128Hi16Lo);
        resRegLo = _mm256_loadu_si256((__m256i*)residual_ptr);
        resRegHi = _mm256_loadu_si256((__m256i*)(residual_ptr + 8));
        predRegLo = _mm256_add_epi32(predRegLo, resRegLo);
        predRegHi = _mm256_add_epi32(predRegHi, resRegHi);
        resReg = _mm256_packus_epi32(predRegLo, predRegHi);
        resReg = _mm256_packus_epi16(resReg, xmm0);
        recon_0_7 = _mm256_shuffle_epi32(resReg, 0xD8);
        predReg_128Hi = _mm256_extracti128_si256(recon_0_7, 1);
        predReg_128Lo = _mm_slli_epi64(predReg_128Hi, 32);
        predReg_128Lo = _mm_or_si128(_mm256_castsi256_si128(recon_0_7), predReg_128Lo);
        _mm_storeu_si128((__m128i*) (reconPtr), predReg_128Lo);
        predPtr += predStride;
        residual_ptr += residualStride;
        reconPtr += reconStride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}

void PictureAdditionKernel32x32_AV1_SSE2_INTRIN(
    uint8_t  *predPtr,
    uint32_t  predStride,
    int32_t *residual_ptr,
    uint32_t  residualStride,
    uint8_t  *reconPtr,
    uint32_t  reconStride,
    uint32_t  width,
    uint32_t  height,
    int32_t     bd)
{

    __m256i predReg, recon_0_7, xmm0, resReg, predReg_Lo, predReg_Hi,
        predReg_Lo16Lo, predReg_Lo16Hi, predReg_Hi16Lo, predReg_Hi16Hi, resReg1, resReg2, resReg3, resReg4;
    __m128i  predReg_128Lo, predReg_128Hi;
    uint32_t y;
    xmm0 = _mm256_setzero_si256();

    for (y = 0; y < 32; ++y) {

        predReg = _mm256_loadu_si256((__m256i*)predPtr);
        predReg_Lo = _mm256_unpacklo_epi8(predReg, xmm0);
        predReg_Hi = _mm256_unpackhi_epi8(predReg, xmm0);
        predReg_Lo16Lo = _mm256_unpacklo_epi16(predReg_Lo, xmm0);
        predReg_Lo16Hi = _mm256_unpackhi_epi16(predReg_Lo, xmm0);
        predReg_Hi16Lo = _mm256_unpacklo_epi16(predReg_Hi, xmm0);
        predReg_Hi16Hi = _mm256_unpackhi_epi16(predReg_Hi, xmm0);

        predReg_Lo = _mm256_set_m128i(_mm256_castsi256_si128(predReg_Lo16Hi), _mm256_castsi256_si128(predReg_Lo16Lo));
        predReg_Hi = _mm256_set_m128i(_mm256_castsi256_si128(predReg_Hi16Hi), _mm256_castsi256_si128(predReg_Hi16Lo));
        predReg_Lo16Lo = _mm256_set_m128i(_mm256_extracti128_si256(predReg_Lo16Hi, 1), _mm256_extracti128_si256(predReg_Lo16Lo, 1));
        predReg_Hi16Hi = _mm256_set_m128i(_mm256_extracti128_si256(predReg_Hi16Hi, 1), _mm256_extracti128_si256(predReg_Hi16Lo, 1));

        resReg1 = _mm256_loadu_si256((__m256i*)residual_ptr);
        resReg2 = _mm256_loadu_si256((__m256i*)(residual_ptr + 8));
        resReg3 = _mm256_loadu_si256((__m256i*)(residual_ptr + 16));
        resReg4 = _mm256_loadu_si256((__m256i*)(residual_ptr + 24));

        resReg1 = _mm256_add_epi32(predReg_Lo, resReg1);
        resReg2 = _mm256_add_epi32(predReg_Hi, resReg2);
        resReg3 = _mm256_add_epi32(predReg_Lo16Lo, resReg3);
        resReg4 = _mm256_add_epi32(predReg_Hi16Hi, resReg4);


        resReg = _mm256_packus_epi32(resReg1, resReg2);
        resReg = _mm256_packus_epi16(resReg, xmm0);
        recon_0_7 = _mm256_shuffle_epi32(resReg, 0xD8);
        predReg_128Hi = _mm256_extracti128_si256(recon_0_7, 1);
        predReg_128Lo = _mm_slli_epi64(predReg_128Hi, 32);
        predReg_128Lo = _mm_or_si128(_mm256_castsi256_si128(recon_0_7), predReg_128Lo);

        resReg = _mm256_packus_epi32(resReg3, resReg4);
        resReg = _mm256_packus_epi16(resReg, xmm0);
        recon_0_7 = _mm256_shuffle_epi32(resReg, 0xD8);
        predReg_128Hi = _mm256_extracti128_si256(recon_0_7, 1);
        predReg_128Hi = _mm_slli_epi64(predReg_128Hi, 32);
        predReg_128Hi = _mm_or_si128(_mm256_castsi256_si128(recon_0_7), predReg_128Hi);
        recon_0_7 = _mm256_set_m128i(predReg_128Hi, predReg_128Lo);
        _mm256_storeu_si256((__m256i*)reconPtr, recon_0_7);

        predPtr += predStride;
        residual_ptr += residualStride;
        reconPtr += reconStride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}
void PictureAdditionKernel64x64_AV1_SSE2_INTRIN(
    uint8_t  *predPtr,
    uint32_t  predStride,
    int32_t *residual_ptr,
    uint32_t  residualStride,
    uint8_t  *reconPtr,
    uint32_t  reconStride,
    uint32_t  width,
    uint32_t  height,
    int32_t     bd)
{

    __m256i predReg, recon_0_7, xmm0, resReg, predReg_Lo, predReg_Hi,
        predReg_Lo16Lo, predReg_Lo16Hi, predReg_Hi16Lo, predReg_Hi16Hi, resReg1, resReg2, resReg3, resReg4;
    __m128i  predReg_128Lo, predReg_128Hi;
    uint32_t y;
    xmm0 = _mm256_setzero_si256();

    for (y = 0; y < 64; ++y) {

        predReg = _mm256_loadu_si256((__m256i*)predPtr);
        predReg_Lo = _mm256_unpacklo_epi8(predReg, xmm0);
        predReg_Hi = _mm256_unpackhi_epi8(predReg, xmm0);
        predReg_Lo16Lo = _mm256_unpacklo_epi16(predReg_Lo, xmm0);
        predReg_Lo16Hi = _mm256_unpackhi_epi16(predReg_Lo, xmm0);
        predReg_Hi16Lo = _mm256_unpacklo_epi16(predReg_Hi, xmm0);
        predReg_Hi16Hi = _mm256_unpackhi_epi16(predReg_Hi, xmm0);

        predReg_Lo = _mm256_set_m128i(_mm256_castsi256_si128(predReg_Lo16Hi), _mm256_castsi256_si128(predReg_Lo16Lo));
        predReg_Hi = _mm256_set_m128i(_mm256_castsi256_si128(predReg_Hi16Hi), _mm256_castsi256_si128(predReg_Hi16Lo));
        predReg_Lo16Lo = _mm256_set_m128i(_mm256_extracti128_si256(predReg_Lo16Hi, 1), _mm256_extracti128_si256(predReg_Lo16Lo, 1));
        predReg_Hi16Hi = _mm256_set_m128i(_mm256_extracti128_si256(predReg_Hi16Hi, 1), _mm256_extracti128_si256(predReg_Hi16Lo, 1));

        resReg1 = _mm256_loadu_si256((__m256i*)residual_ptr);
        resReg2 = _mm256_loadu_si256((__m256i*)(residual_ptr + 8));
        resReg3 = _mm256_loadu_si256((__m256i*)(residual_ptr + 16));
        resReg4 = _mm256_loadu_si256((__m256i*)(residual_ptr + 24));

        resReg1 = _mm256_add_epi32(predReg_Lo, resReg1);
        resReg2 = _mm256_add_epi32(predReg_Hi, resReg2);
        resReg3 = _mm256_add_epi32(predReg_Lo16Lo, resReg3);
        resReg4 = _mm256_add_epi32(predReg_Hi16Hi, resReg4);


        resReg = _mm256_packus_epi32(resReg1, resReg2);
        resReg = _mm256_packus_epi16(resReg, xmm0);
        recon_0_7 = _mm256_shuffle_epi32(resReg, 0xD8);
        predReg_128Hi = _mm256_extracti128_si256(recon_0_7, 1);
        predReg_128Lo = _mm_slli_epi64(predReg_128Hi, 32);
        predReg_128Lo = _mm_or_si128(_mm256_castsi256_si128(recon_0_7), predReg_128Lo);

        resReg = _mm256_packus_epi32(resReg3, resReg4);
        resReg = _mm256_packus_epi16(resReg, xmm0);
        recon_0_7 = _mm256_shuffle_epi32(resReg, 0xD8);
        predReg_128Hi = _mm256_extracti128_si256(recon_0_7, 1);
        predReg_128Hi = _mm_slli_epi64(predReg_128Hi, 32);
        predReg_128Hi = _mm_or_si128(_mm256_castsi256_si128(recon_0_7), predReg_128Hi);
        recon_0_7 = _mm256_set_m128i(predReg_128Hi, predReg_128Lo);
        _mm256_storeu_si256((__m256i*)reconPtr, recon_0_7);

        predReg = _mm256_loadu_si256((__m256i*)(predPtr + 32));
        predReg_Lo = _mm256_unpacklo_epi8(predReg, xmm0);
        predReg_Hi = _mm256_unpackhi_epi8(predReg, xmm0);
        predReg_Lo16Lo = _mm256_unpacklo_epi16(predReg_Lo, xmm0);
        predReg_Lo16Hi = _mm256_unpackhi_epi16(predReg_Lo, xmm0);
        predReg_Hi16Lo = _mm256_unpacklo_epi16(predReg_Hi, xmm0);
        predReg_Hi16Hi = _mm256_unpackhi_epi16(predReg_Hi, xmm0);

        predReg_Lo = _mm256_set_m128i(_mm256_castsi256_si128(predReg_Lo16Hi), _mm256_castsi256_si128(predReg_Lo16Lo));
        predReg_Hi = _mm256_set_m128i(_mm256_castsi256_si128(predReg_Hi16Hi), _mm256_castsi256_si128(predReg_Hi16Lo));
        predReg_Lo16Lo = _mm256_set_m128i(_mm256_extracti128_si256(predReg_Lo16Hi, 1), _mm256_extracti128_si256(predReg_Lo16Lo, 1));
        predReg_Hi16Hi = _mm256_set_m128i(_mm256_extracti128_si256(predReg_Hi16Hi, 1), _mm256_extracti128_si256(predReg_Hi16Lo, 1));

        resReg1 = _mm256_loadu_si256((__m256i*)(residual_ptr + 32));
        resReg2 = _mm256_loadu_si256((__m256i*)(residual_ptr + 40));
        resReg3 = _mm256_loadu_si256((__m256i*)(residual_ptr + 48));
        resReg4 = _mm256_loadu_si256((__m256i*)(residual_ptr + 56));

        resReg1 = _mm256_add_epi32(predReg_Lo, resReg1);
        resReg2 = _mm256_add_epi32(predReg_Hi, resReg2);
        resReg3 = _mm256_add_epi32(predReg_Lo16Lo, resReg3);
        resReg4 = _mm256_add_epi32(predReg_Hi16Hi, resReg4);


        resReg = _mm256_packus_epi32(resReg1, resReg2);
        resReg = _mm256_packus_epi16(resReg, xmm0);
        recon_0_7 = _mm256_shuffle_epi32(resReg, 0xD8);
        predReg_128Hi = _mm256_extracti128_si256(recon_0_7, 1);
        predReg_128Lo = _mm_slli_epi64(predReg_128Hi, 32);
        predReg_128Lo = _mm_or_si128(_mm256_castsi256_si128(recon_0_7), predReg_128Lo);

        resReg = _mm256_packus_epi32(resReg3, resReg4);
        resReg = _mm256_packus_epi16(resReg, xmm0);
        recon_0_7 = _mm256_shuffle_epi32(resReg, 0xD8);
        predReg_128Hi = _mm256_extracti128_si256(recon_0_7, 1);
        predReg_128Hi = _mm_slli_epi64(predReg_128Hi, 32);
        predReg_128Hi = _mm_or_si128(_mm256_castsi256_si128(recon_0_7), predReg_128Hi);
        recon_0_7 = _mm256_set_m128i(predReg_128Hi, predReg_128Lo);
        _mm256_storeu_si256((__m256i*)(reconPtr + 32), recon_0_7);

        predPtr += predStride;
        residual_ptr += residualStride;
        reconPtr += reconStride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}

void FullDistortionKernel32Bits_AVX2(
    int32_t  *coeff,
    uint32_t   coeffStride,
    int32_t  *reconCoeff,
    uint32_t   reconCoeffStride,
    uint64_t   distortionResult[DIST_CALC_TOTAL],
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t rowCount, colCount;
    __m256i sum1 = _mm256_setzero_si256();
    __m256i sum2 = _mm256_setzero_si256();
    __m128i temp1, temp2, temp3;

    rowCount = areaHeight;
    do
    {
        int32_t *coeffTemp = coeff;
        int32_t *reconCoeffTemp = reconCoeff;

        colCount = areaWidth / 4;
        do
        {
            __m128i x0, y0;
            __m256i x, y, z;
            x0 = _mm_loadu_si128((__m128i *)(coeffTemp));
            y0 = _mm_loadu_si128((__m128i *)(reconCoeffTemp));
            x = _mm256_cvtepi32_epi64(x0);
            y = _mm256_cvtepi32_epi64(y0);
            z= _mm256_mul_epi32(x, x);
            sum2 = _mm256_add_epi64(sum2, z);
            x = _mm256_sub_epi64(x, y);
            x = _mm256_mul_epi32(x, x);
            sum1 = _mm256_add_epi32(sum1, x);
            coeffTemp += 4;
            reconCoeffTemp += 4;
        } while (--colCount);

        coeff += coeffStride;
        reconCoeff += reconCoeffStride;
        rowCount -= 1;
    } while (rowCount > 0);

    temp1 = _mm256_extracti128_si256(sum1, 0);
    temp2 = _mm256_extracti128_si256(sum1, 1);
    temp1 = _mm_add_epi64(temp1, temp2);
    temp2 = _mm_shuffle_epi32(temp1, 0x4e);
    temp3 = _mm_add_epi64(temp1, temp2);
    temp1 = _mm256_extracti128_si256(sum2, 0);
    temp2 = _mm256_extracti128_si256(sum2, 1);
    temp1 = _mm_add_epi64(temp1, temp2);
    temp2 = _mm_shuffle_epi32(temp1, 0x4e);
    temp1 = _mm_add_epi64(temp1, temp2);
    temp1 = _mm_unpacklo_epi64(temp3, temp1);

    _mm_storeu_si128((__m128i *)distortionResult, temp1);
}

void FullDistortionKernelCbfZero32Bits_AVX2(
    int32_t  *coeff,
    uint32_t   coeffStride,
    int32_t  *reconCoeff,
    uint32_t   reconCoeffStride,
    uint64_t   distortionResult[DIST_CALC_TOTAL],
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t rowCount, colCount;
    __m256i sum = _mm256_setzero_si256();
    __m128i temp1, temp2;

    rowCount = areaHeight;
    do
    {
        int32_t *coeffTemp = coeff;

        colCount = areaWidth / 4;
        do
        {
            __m128i x0;
            __m256i y0, z0;
            x0 = _mm_loadu_si128((__m128i *)(coeffTemp));
            coeffTemp += 4;
            y0 = _mm256_cvtepi32_epi64(x0);
            z0 = _mm256_mul_epi32(y0, y0);
            sum = _mm256_add_epi64(sum, z0);
        } while (--colCount);

        coeff += coeffStride;
        reconCoeff += coeffStride;
        rowCount -= 1;
    } while (rowCount > 0);

    temp1 = _mm256_extracti128_si256(sum, 0);
    temp2 = _mm256_extracti128_si256(sum, 1);
    temp1 = _mm_add_epi64(temp1, temp2);
    temp2 = _mm_shuffle_epi32(temp1, 0x4e);
    temp1 = _mm_add_epi64(temp1, temp2);
    _mm_storeu_si128((__m128i *)distortionResult, temp1);
    (void)reconCoeffStride;
}

static INLINE void _mm_storeh_epi64(__m128i *const d, const __m128i s) {
    _mm_storeh_pi((__m64 *)d, _mm_castsi128_ps(s));
}

void ResidualKernel4x4_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256i in = load8bit_4x4_avx2(input, inputStride);
    const __m256i pr = load8bit_4x4_avx2(pred, predStride);
    const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
    const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
    const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
    const __m128i re_01 = _mm256_extracti128_si256(re_lo, 0);
    const __m128i re_23 = _mm256_extracti128_si256(re_lo, 1);
    (void)areaWidth;
    (void)areaHeight;

    _mm_storel_epi64((__m128i*)(residual + 0 * residualStride), re_01);
    _mm_storeh_epi64((__m128i*)(residual + 1 * residualStride), re_01);
    _mm_storel_epi64((__m128i*)(residual + 2 * residualStride), re_23);
    _mm_storeh_epi64((__m128i*)(residual + 3 * residualStride), re_23);
}
void ResidualKernel4x8_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 8; y += 4) {
        const __m256i in = load8bit_4x4_avx2(input, inputStride);
        const __m256i pr = load8bit_4x4_avx2(pred, predStride);
        const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m128i re_01 = _mm256_extracti128_si256(re_lo, 0);
        const __m128i re_23 = _mm256_extracti128_si256(re_lo, 1);

        _mm_storel_epi64((__m128i*)(residual + 0 * residualStride), re_01);
        _mm_storeh_epi64((__m128i*)(residual + 1 * residualStride), re_01);
        _mm_storel_epi64((__m128i*)(residual + 2 * residualStride), re_23);
        _mm_storeh_epi64((__m128i*)(residual + 3 * residualStride), re_23);
        input += 4 * inputStride;
        pred += 4 * predStride;
        residual += 4 * residualStride;
    }
}
void ResidualKernel4x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 16; y += 4) {
        const __m256i in = load8bit_4x4_avx2(input, inputStride);
        const __m256i pr = load8bit_4x4_avx2(pred, predStride);
        const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m128i re_01 = _mm256_extracti128_si256(re_lo, 0);
        const __m128i re_23 = _mm256_extracti128_si256(re_lo, 1);

        _mm_storel_epi64((__m128i*)(residual + 0 * residualStride), re_01);
        _mm_storeh_epi64((__m128i*)(residual + 1 * residualStride), re_01);
        _mm_storel_epi64((__m128i*)(residual + 2 * residualStride), re_23);
        _mm_storeh_epi64((__m128i*)(residual + 3 * residualStride), re_23);
        input += 4 * inputStride;
        pred += 4 * predStride;
        residual += 4 * residualStride;
    }
}
void ResidualKernel8x32_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 32; y += 4) {
        const __m256i in = load8bit_8x4_avx2(input, inputStride);
        const __m256i pr = load8bit_8x4_avx2(pred, predStride);
        const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);
        const __m128i re_0 = _mm256_extracti128_si256(re_lo, 0);
        const __m128i re_1 = _mm256_extracti128_si256(re_hi, 0);
        const __m128i re_2 = _mm256_extracti128_si256(re_lo, 1);
        const __m128i re_3 = _mm256_extracti128_si256(re_hi, 1);

        _mm_storeu_si128((__m128i*)(residual + 0 * residualStride), re_0);
        _mm_storeu_si128((__m128i*)(residual + 1 * residualStride), re_1);
        _mm_storeu_si128((__m128i*)(residual + 2 * residualStride), re_2);
        _mm_storeu_si128((__m128i*)(residual + 3 * residualStride), re_3);
        input += 4 * inputStride;
        pred += 4 * predStride;
        residual += 4 * residualStride;
    }
}
void ResidualKernel8x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 16; y += 4) {
        const __m256i in = load8bit_8x4_avx2(input, inputStride);
        const __m256i pr = load8bit_8x4_avx2(pred, predStride);
        const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);
        const __m128i re_0 = _mm256_extracti128_si256(re_lo, 0);
        const __m128i re_1 = _mm256_extracti128_si256(re_hi, 0);
        const __m128i re_2 = _mm256_extracti128_si256(re_lo, 1);
        const __m128i re_3 = _mm256_extracti128_si256(re_hi, 1);

        _mm_storeu_si128((__m128i*)(residual + 0 * residualStride), re_0);
        _mm_storeu_si128((__m128i*)(residual + 1 * residualStride), re_1);
        _mm_storeu_si128((__m128i*)(residual + 2 * residualStride), re_2);
        _mm_storeu_si128((__m128i*)(residual + 3 * residualStride), re_3);
        input += 4 * inputStride;
        pred += 4 * predStride;
        residual += 4 * residualStride;
    }
}

void ResidualKernel8x8_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 8; y += 4) {
        const __m256i in = load8bit_8x4_avx2(input, inputStride);
        const __m256i pr = load8bit_8x4_avx2(pred, predStride);
        const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);
        const __m128i re_0 = _mm256_extracti128_si256(re_lo, 0);
        const __m128i re_1 = _mm256_extracti128_si256(re_hi, 0);
        const __m128i re_2 = _mm256_extracti128_si256(re_lo, 1);
        const __m128i re_3 = _mm256_extracti128_si256(re_hi, 1);

        _mm_storeu_si128((__m128i*)(residual + 0 * residualStride), re_0);
        _mm_storeu_si128((__m128i*)(residual + 1 * residualStride), re_1);
        _mm_storeu_si128((__m128i*)(residual + 2 * residualStride), re_2);
        _mm_storeu_si128((__m128i*)(residual + 3 * residualStride), re_3);
        input += 4 * inputStride;
        pred += 4 * predStride;
        residual += 4 * residualStride;
    }
}

void ResidualKernel8x4_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    //for (y = 0; y < 8; y += 4) {
        const __m256i in = load8bit_8x4_avx2(input, inputStride);
        const __m256i pr = load8bit_8x4_avx2(pred, predStride);
        const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);
        const __m128i re_0 = _mm256_extracti128_si256(re_lo, 0);
        const __m128i re_1 = _mm256_extracti128_si256(re_hi, 0);
        const __m128i re_2 = _mm256_extracti128_si256(re_lo, 1);
        const __m128i re_3 = _mm256_extracti128_si256(re_hi, 1);

        _mm_storeu_si128((__m128i*)(residual + 0 * residualStride), re_0);
        _mm_storeu_si128((__m128i*)(residual + 1 * residualStride), re_1);
        _mm_storeu_si128((__m128i*)(residual + 2 * residualStride), re_2);
        _mm_storeu_si128((__m128i*)(residual + 3 * residualStride), re_3);
        //input += 4 * inputStride;
        //pred += 4 * predStride;
        //residual += 4 * residualStride;
    //}
}

void ResidualKernel16x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 16; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, inputStride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, predStride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residualStride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residualStride), re_hi);
        input += 2 * inputStride;
        pred += 2 * predStride;
        residual += 2 * residualStride;
    }
}

void ResidualKernel16x4_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 4; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, inputStride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, predStride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residualStride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residualStride), re_hi);
        input += 2 * inputStride;
        pred += 2 * predStride;
        residual += 2 * residualStride;
    }
}

void ResidualKernel16x8_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 8; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, inputStride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, predStride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residualStride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residualStride), re_hi);
        input += 2 * inputStride;
        pred += 2 * predStride;
        residual += 2 * residualStride;
    }
}

void ResidualKernel16x32_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 32; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, inputStride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, predStride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residualStride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residualStride), re_hi);
        input += 2 * inputStride;
        pred += 2 * predStride;
        residual += 2 * residualStride;
    }
}

void ResidualKernel16x64_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 64; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, inputStride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, predStride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residualStride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residualStride), re_hi);
        input += 2 * inputStride;
        pred += 2 * predStride;
        residual += 2 * residualStride;
    }
}
static INLINE void ResidualKernel32_AVX2(const uint8_t *const input,
    const uint8_t *const pred, int16_t *const residual)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256i in0 = _mm256_loadu_si256((__m256i *)input);
    const __m256i pr0 = _mm256_loadu_si256((__m256i *)pred);
    const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
    const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
    const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
    const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
    const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
    const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
    const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
    const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);
    _mm256_storeu_si256((__m256i*)(residual + 0x00), re_lo);
    _mm256_storeu_si256((__m256i*)(residual + 0x10), re_hi);
}

void ResidualKernel32x8_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 8; ++y) {
        ResidualKernel32_AVX2(input, pred, residual);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}

void ResidualKernel32x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 16; ++y) {
        ResidualKernel32_AVX2(input, pred, residual);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}

void ResidualKernel32x32_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 32; ++y) {
        ResidualKernel32_AVX2(input, pred, residual);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}


void ResidualKernel32x64_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 64; ++y) {
        ResidualKernel32_AVX2(input, pred, residual);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}

void ResidualKernel64x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 16; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}

void ResidualKernel64x32_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 32; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}

void ResidualKernel64x64_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 64; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}

void ResidualKernel64x128_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 128; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}


void ResidualKernel128x128_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 128; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        ResidualKernel32_AVX2(input + 0x40, pred + 0x40, residual + 0x40);
        ResidualKernel32_AVX2(input + 0x60, pred + 0x60, residual + 0x60);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}

void ResidualKernel128x64_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t y;
    (void)areaWidth;
    (void)areaHeight;

    for (y = 0; y < 64; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        ResidualKernel32_AVX2(input + 0x40, pred + 0x40, residual + 0x40);
        ResidualKernel32_AVX2(input + 0x60, pred + 0x60, residual + 0x60);
        input += inputStride;
        pred += predStride;
        residual += residualStride;
    }
}
void ResidualKernel_avx2(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    if (areaWidth == areaHeight) {
        switch (areaWidth) {
        case 4:
            ResidualKernel4x4_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
            break;
        case 8:
            ResidualKernel8x8_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
            break;
        case 16:
            ResidualKernel16x16_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
            break;
        case 32:
            ResidualKernel32x32_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
            break;
        case 64:
            ResidualKernel64x64_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
            break;
        case 128:
            ResidualKernel128x128_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
            break;
        }
    }
    else {
        if (areaWidth < areaHeight) {
            if (areaWidth + areaWidth == areaHeight) {
                switch (areaWidth) {
                case 4:
                    ResidualKernel4x8_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 8:
                    ResidualKernel8x16_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 16:
                    ResidualKernel16x32_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 32:
                    ResidualKernel32x64_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 64:
                    ResidualKernel64x128_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                }
            }
            else {
                switch (areaWidth) {
                case 4:
                    ResidualKernel4x16_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 8:
                    ResidualKernel8x32_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 16:
                    ResidualKernel16x64_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                }
            }
        }
        else {
            if (areaHeight + areaHeight == areaWidth) {
                switch (areaHeight) {
                case 4:
                    ResidualKernel8x4_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 8:
                    ResidualKernel16x8_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 16:
                    ResidualKernel32x16_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 32:
                    ResidualKernel64x32_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 64:
                    ResidualKernel128x64_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                }
            }
            else {
                switch (areaHeight) {
                case 4:
                    ResidualKernel16x4_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 8:
                    ResidualKernel32x8_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                case 16:
                    ResidualKernel64x16_AVX2_INTRIN(input, inputStride, pred, predStride, residual, residualStride, areaWidth, areaHeight);
                    break;
                }
            }
        }
    }
}
