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

void compressed_packmsb_avx2_intrin(
    uint8_t     *in8_bit_buffer,
    uint32_t     in8_stride,
    uint8_t     *inn_bit_buffer,
    uint16_t    *out16_bit_buffer,
    uint32_t     inn_stride,
    uint32_t     out_stride,
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


            in2Bit = _mm_loadu_si128((__m128i*)inn_bit_buffer); //2 Lines Stored in 1D format-Could be replaced by 2 _mm_loadl_epi64

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

            in8Bit = _mm256_loadu_si256((__m256i*)in8_bit_buffer);
            in8BitStride = _mm256_loadu_si256((__m256i*)(in8_bit_buffer + in8_stride));


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

            _mm256_store_si256((__m256i*) out16_bit_buffer, out0_15);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + 16), out16_31);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + out_stride), out_s0_s15);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + out_stride + 16), out_s16_s31);

            in8_bit_buffer += in8_stride << 1;
            inn_bit_buffer += inn_stride << 1;
            out16_bit_buffer += out_stride << 1;
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

            in2Bit = _mm_loadu_si128((__m128i*)inn_bit_buffer);

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

            in8Bit = _mm256_loadu_si256((__m256i*)in8_bit_buffer);
            in8Bit32 = _mm256_loadu_si256((__m256i*)(in8_bit_buffer + 32));

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

            _mm256_store_si256((__m256i*) out16_bit_buffer, out_0_15);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + 16), out16_31);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + 32), out32_47);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + 48), out_48_63);

            in8_bit_buffer += in8_stride;
            inn_bit_buffer += inn_stride;
            out16_bit_buffer += out_stride;

        }

    }

}
#if defined(_MSC_VER)
//nclude <intrin.h>
#endif

void c_pack_avx2_intrin(
    const uint8_t     *inn_bit_buffer,
    uint32_t     inn_stride,
    uint8_t     *in_compn_bit_buffer,
    uint32_t     out_stride,
    uint8_t    *local_cache,

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


            inNBit = _mm256_loadu_si256((__m256i*)inn_bit_buffer);

            ext0 = _mm256_and_si256(inNBit, msk0);
            ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
            ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
            ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

            ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

            ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));


            extp = _mm256_packus_epi32(ext0123, ext0123n);
            extp = _mm256_packus_epi16(extp, extp);

            _mm_storel_epi64((__m128i*) in_compn_bit_buffer, _mm256_castsi256_si128(extp));
            in_compn_bit_buffer += 8;
            inn_bit_buffer += inn_stride;


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

            uint8_t* localPtr = local_cache;


            for (y = 0; y < height; y++)
            {


                inNBit = _mm256_loadu_si256((__m256i*)inn_bit_buffer);


                ext0 = _mm256_and_si256(inNBit, msk0);
                ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
                ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
                ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

                ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

                ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));

                extp = _mm256_packus_epi32(ext0123, ext0123n);
                extp = _mm256_packus_epi16(extp, extp);


                inNBit = _mm256_loadu_si256((__m256i*)(inn_bit_buffer + 32));

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
                    _mm256_stream_si256((__m256i*)&in_compn_bit_buffer[0], c0);
                    _mm256_stream_si256((__m256i*)&in_compn_bit_buffer[32], c1);
                    in_compn_bit_buffer += 4 * out_stride;
                }

                inn_bit_buffer += inn_stride;

            }

        }
        else {

            //One row per iter
            for (y = 0; y < height; y++)
            {


                inNBit = _mm256_loadu_si256((__m256i*)inn_bit_buffer);


                ext0 = _mm256_and_si256(inNBit, msk0);
                ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
                ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
                ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

                ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

                ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));

                extp = _mm256_packus_epi32(ext0123, ext0123n);
                extp = _mm256_packus_epi16(extp, extp);


                inNBit = _mm256_loadu_si256((__m256i*)(inn_bit_buffer + 32));

                ext0 = _mm256_and_si256(inNBit, msk0);
                ext1 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 1 * 8 + 2), msk1);
                ext2 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 2 * 8 + 4), msk2);
                ext3 = _mm256_and_si256(_mm256_srli_epi32(inNBit, 3 * 8 + 6), msk3);

                ext0123 = _mm256_or_si256(_mm256_or_si256(ext0, ext1), _mm256_or_si256(ext2, ext3));

                ext0123n = _mm256_castsi128_si256(_mm256_extracti128_si256(ext0123, 1));


                extp1 = _mm256_packus_epi32(ext0123, ext0123n);
                extp1 = _mm256_packus_epi16(extp1, extp1);

                extp = _mm256_unpacklo_epi64(extp, extp1);

                _mm_storeu_si128((__m128i*)  in_compn_bit_buffer, _mm256_castsi256_si128(extp));

                in_compn_bit_buffer += out_stride;

                inn_bit_buffer += inn_stride;

            }

        }

    }

}


void eb_enc_msb_pack2d_avx2_intrin_al(
    uint8_t     *in8_bit_buffer,
    uint32_t     in8_stride,
    uint8_t     *inn_bit_buffer,
    uint16_t    *out16_bit_buffer,
    uint32_t     inn_stride,
    uint32_t     out_stride,
    uint32_t     width,
    uint32_t     height)
{
    //(outPixel | nBitPixel) concatenation is done with unpacklo_epi8 and unpackhi_epi8

    uint32_t y, x;

    __m128i out0, out1;


    if (width == 4)
    {
        for (y = 0; y < height; y += 2) {

            out0 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)inn_bit_buffer), _mm_cvtsi32_si128(*(uint32_t *)in8_bit_buffer)), 6);
            out1 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(inn_bit_buffer + inn_stride)), _mm_cvtsi32_si128(*(uint32_t *)(in8_bit_buffer + in8_stride))), 6);

            _mm_storel_epi64((__m128i*) out16_bit_buffer, out0);
            _mm_storel_epi64((__m128i*) (out16_bit_buffer + out_stride), out1);

            in8_bit_buffer += in8_stride << 1;
            inn_bit_buffer += inn_stride << 1;
            out16_bit_buffer += out_stride << 1;
        }
    }
    else if (width == 8)
    {
        for (y = 0; y < height; y += 2) {

            out0 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)inn_bit_buffer), _mm_loadl_epi64((__m128i*)in8_bit_buffer)), 6);
            out1 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(inn_bit_buffer + inn_stride)), _mm_loadl_epi64((__m128i*)(in8_bit_buffer + in8_stride))), 6);

            _mm_storeu_si128((__m128i*) out16_bit_buffer, out0);
            _mm_storeu_si128((__m128i*) (out16_bit_buffer + out_stride), out1);

            in8_bit_buffer += in8_stride << 1;
            inn_bit_buffer += inn_stride << 1;
            out16_bit_buffer += out_stride << 1;
        }
    }
    else if (width == 16)
    {
        __m128i inNBit, in8Bit, inNBitStride, in8BitStride, out0, out1, out2, out3;

        for (y = 0; y < height; y += 2) {
            inNBit = _mm_loadu_si128((__m128i*)inn_bit_buffer);
            in8Bit = _mm_loadu_si128((__m128i*)in8_bit_buffer);
            inNBitStride = _mm_loadu_si128((__m128i*)(inn_bit_buffer + inn_stride));
            in8BitStride = _mm_loadu_si128((__m128i*)(in8_bit_buffer + in8_stride));

            out0 = _mm_srli_epi16(_mm_unpacklo_epi8(inNBit, in8Bit), 6);
            out1 = _mm_srli_epi16(_mm_unpackhi_epi8(inNBit, in8Bit), 6);
            out2 = _mm_srli_epi16(_mm_unpacklo_epi8(inNBitStride, in8BitStride), 6);
            out3 = _mm_srli_epi16(_mm_unpackhi_epi8(inNBitStride, in8BitStride), 6);

            _mm_storeu_si128((__m128i*) out16_bit_buffer, out0);
            _mm_storeu_si128((__m128i*) (out16_bit_buffer + 8), out1);
            _mm_storeu_si128((__m128i*) (out16_bit_buffer + out_stride), out2);
            _mm_storeu_si128((__m128i*) (out16_bit_buffer + out_stride + 8), out3);

            in8_bit_buffer += in8_stride << 1;
            inn_bit_buffer += inn_stride << 1;
            out16_bit_buffer += out_stride << 1;
        }
    }
    else if (width == 32)
    {
        __m256i inNBit, in8Bit, inNBitStride, in8BitStride, concat0, concat1, concat2, concat3;
        __m256i out0_15, out16_31, out_s0_s15, out_s16_s31;

        for (y = 0; y < height; y += 2) {

            inNBit = _mm256_loadu_si256((__m256i*)inn_bit_buffer);
            in8Bit = _mm256_loadu_si256((__m256i*)in8_bit_buffer);
            inNBitStride = _mm256_loadu_si256((__m256i*)(inn_bit_buffer + inn_stride));
            in8BitStride = _mm256_loadu_si256((__m256i*)(in8_bit_buffer + in8_stride));

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

            _mm256_store_si256((__m256i*) out16_bit_buffer, out0_15);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + 16), out16_31);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + out_stride), out_s0_s15);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + out_stride + 16), out_s16_s31);

            in8_bit_buffer += in8_stride << 1;
            //inn_bit_buffer += inn_stride << 1;
            inn_bit_buffer += inn_stride * 2;
            out16_bit_buffer += out_stride << 1;
        }
    }
    else if (width == 64)
    {
        __m256i inNBit, in8Bit, inNBitStride, in8BitStride, inNBit32, in8Bit32, inNBitStride32, in8BitStride32;
        __m256i concat0, concat1, concat2, concat3, concat4, concat5, concat6, concat7;
        __m256i out_0_15, out16_31, out32_47, out_48_63, out_s0_s15, out_s16_s31, out_s32_s47, out_s48_s63;

        for (y = 0; y < height; y += 2) {

            inNBit = _mm256_loadu_si256((__m256i*)inn_bit_buffer);
            in8Bit = _mm256_loadu_si256((__m256i*)in8_bit_buffer);
            inNBit32 = _mm256_loadu_si256((__m256i*)(inn_bit_buffer + 32));
            in8Bit32 = _mm256_loadu_si256((__m256i*)(in8_bit_buffer + 32));
            inNBitStride = _mm256_loadu_si256((__m256i*)(inn_bit_buffer + inn_stride));
            in8BitStride = _mm256_loadu_si256((__m256i*)(in8_bit_buffer + in8_stride));
            inNBitStride32 = _mm256_loadu_si256((__m256i*)(inn_bit_buffer + inn_stride + 32));
            in8BitStride32 = _mm256_loadu_si256((__m256i*)(in8_bit_buffer + in8_stride + 32));
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

            _mm256_store_si256((__m256i*) out16_bit_buffer, out_0_15);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + 16), out16_31);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + 32), out32_47);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + 48), out_48_63);

            _mm256_store_si256((__m256i*) (out16_bit_buffer + out_stride), out_s0_s15);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + out_stride + 16), out_s16_s31);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + out_stride + 32), out_s32_s47);
            _mm256_store_si256((__m256i*) (out16_bit_buffer + out_stride + 48), out_s48_s63);

            in8_bit_buffer += in8_stride << 1;
            //inn_bit_buffer += inn_stride << 1;
            inn_bit_buffer += inn_stride * 2;
            out16_bit_buffer += out_stride << 1;
        }
    }
    else
    {
        uint32_t innStrideDiff = 2 * inn_stride;
        uint32_t in8StrideDiff = 2 * in8_stride;
        uint32_t outStrideDiff = 2 * out_stride;
        innStrideDiff -= width;
        in8StrideDiff -= width;
        outStrideDiff -= width;

        if (!(width & 7)) {

            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 8) {

                    out0 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)inn_bit_buffer), _mm_loadl_epi64((__m128i*)in8_bit_buffer)), 6);
                    out1 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i*)(inn_bit_buffer + inn_stride)), _mm_loadl_epi64((__m128i*)(in8_bit_buffer + in8_stride))), 6);

                    _mm_storeu_si128((__m128i*) out16_bit_buffer, out0);
                    _mm_storeu_si128((__m128i*) (out16_bit_buffer + out_stride), out1);

                    in8_bit_buffer += 8;
                    inn_bit_buffer += 8;
                    out16_bit_buffer += 8;
                }
                in8_bit_buffer += in8StrideDiff;
                inn_bit_buffer += innStrideDiff;
                out16_bit_buffer += outStrideDiff;
            }
        }
        else {
            for (x = 0; x < height; x += 2) {
                for (y = 0; y < width; y += 4) {

                    out0 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)inn_bit_buffer), _mm_cvtsi32_si128(*(uint32_t *)in8_bit_buffer)), 6);
                    out1 = _mm_srli_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(inn_bit_buffer + inn_stride)), _mm_cvtsi32_si128(*(uint32_t *)(in8_bit_buffer + in8_stride))), 6);

                    _mm_storel_epi64((__m128i*) out16_bit_buffer, out0);
                    _mm_storel_epi64((__m128i*) (out16_bit_buffer + out_stride), out1);

                    in8_bit_buffer += 4;
                    inn_bit_buffer += 4;
                    out16_bit_buffer += 4;
                }
                in8_bit_buffer += in8StrideDiff;
                inn_bit_buffer += innStrideDiff;
                out16_bit_buffer += outStrideDiff;
            }
        }
    }
}

#define ALSTORE  1
#define B256     1

void unpack_avg_avx2_intrin(
    uint16_t *ref16_l0,
    uint32_t  ref_l0_stride,
    uint16_t *ref16_l1,
    uint32_t  ref_l1_stride,
    uint8_t  *dst_ptr,
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
#if ALSTORE
            _mm_store_si128((__m128i*) dst_ptr, avg8_0_U8);
#else
            _mm_storeu_si128((__m128i*) dst_ptr, avg8_0_U8);
#endif

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
#if ALSTORE
            _mm_store_si128((__m128i*)(dst_ptr + dst_stride), avg8_2_U8);
#else
            _mm_storeu_si128((__m128i*)(dst_ptr + dst_stride), avg8_2_U8);
#endif
            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;

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
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 16));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 16));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);

            _mm256_storeu_si256((__m256i *)(dst_ptr), avg8b_32_0);

            //--------
            //Line Two
            //--------
              //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*)(ref16_l0 + ref_l0_stride));
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l0 + ref_l0_stride + 16));

            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));

            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*)(ref16_l1 + ref_l1_stride));
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l1 + ref_l1_stride + 16));

            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);

            _mm256_storeu_si256((__m256i *)(dst_ptr + dst_stride), avg8b_32_0);

#else
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
#if ALSTORE
            _mm_store_si128((__m128i*) dst_ptr, avg8_0_U8);
            _mm_store_si128((__m128i*)(dst_ptr + 16), avg8_1_U8);
#else
            _mm_storeu_si128((__m128i*) dst_ptr, avg8_0_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + 16), avg8_1_U8);
#endif

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
#if ALSTORE
            _mm_store_si128((__m128i*)(dst_ptr + dst_stride), avg8_2_U8);
            _mm_store_si128((__m128i*)(dst_ptr + dst_stride + 16), avg8_3_U8);
#else
            _mm_storeu_si128((__m128i*)(dst_ptr + dst_stride), avg8_2_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + dst_stride + 16), avg8_3_U8);
#endif

#endif
            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;

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
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 48));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 48));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);
            avg8b_32_1 = _mm256_avg_epu8(data8b_32_1_L0, data8b_32_1_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);
            avg8b_32_1 = _mm256_permute4x64_epi64(avg8b_32_1, 216);

            _mm256_storeu_si256((__m256i *)(dst_ptr), avg8b_32_0);
            _mm256_storeu_si256((__m256i *)(dst_ptr + 32), avg8b_32_1);
#else
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
#if ALSTORE
            _mm_store_si128((__m128i*) dst_ptr, avg8_0_U8);
            _mm_store_si128((__m128i*)(dst_ptr + 16), avg8_1_U8);
            _mm_store_si128((__m128i*)(dst_ptr + 32), avg8_2_U8);
            _mm_store_si128((__m128i*)(dst_ptr + 48), avg8_3_U8);
#else
            _mm_storeu_si128((__m128i*) dst_ptr, avg8_0_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + 16), avg8_1_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + 32), avg8_2_U8);
            _mm_storeu_si128((__m128i*)(dst_ptr + 48), avg8_3_U8);
#endif

#endif
            dst_ptr += dst_stride;
            ref16_l0 += ref_l0_stride;
            ref16_l1 += ref_l1_stride;
        }
    }


    return;
}

int32_t  sum_residual8bit_avx2_intrin(
    int16_t * in_ptr,
    uint32_t   size,
    uint32_t   stride_in)
{

    int32_t  sumBlock;

    __m128i in0, in1, in01, in2, in3, in23, sum, sumL, sumH;
    __m256i sum0, sum1, sum2, sum3, sum0L, sum0H, sumT, sum01, sumTPerm;
    uint32_t row_index;

    //Assumption: 9bit or 11bit residual data . for bigger block sizes or bigger bit depths , re-asses the dynamic range of the internal calculation

    if (size == 4) { //SSSE3

        __m128i zer = _mm_setzero_si128();

        in0 = _mm_loadl_epi64((__m128i*)in_ptr);
        in1 = _mm_loadl_epi64((__m128i*)(in_ptr + stride_in));
        in1 = _mm_shuffle_epi32(in1, 0x4A); //01.00.10.10
        in01 = _mm_or_si128(in1, in0);

        in2 = _mm_loadl_epi64((__m128i*)(in_ptr + 2 * stride_in));
        in3 = _mm_loadl_epi64((__m128i*)(in_ptr + 3 * stride_in));
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

        sum = _mm_add_epi16(_mm_loadu_si128((__m128i*)(in_ptr + 0 * stride_in)), _mm_loadu_si128((__m128i*)(in_ptr + 1 * stride_in)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(in_ptr + 2 * stride_in)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(in_ptr + 3 * stride_in)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(in_ptr + 4 * stride_in)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(in_ptr + 5 * stride_in)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(in_ptr + 6 * stride_in)));
        sum = _mm_add_epi16(sum, _mm_loadu_si128((__m128i*)(in_ptr + 7 * stride_in)));

        sum = _mm_hadd_epi16(sum, zer);
        sum = _mm_hadd_epi16(sum, zer);
        sum = _mm_hadd_epi16(sum, zer);

        sum = _mm_cvtepi16_epi32(sum); //the sum is on 16bit, for negative values, we need to extend the sign to the next 16bit, so that the next extraction to int32_t is fine.
        sumBlock = _mm_cvtsi128_si32(sum);

        return sumBlock;

    }
    else if (size == 16) {//AVX2

        sum0 = _mm256_add_epi16(_mm256_loadu_si256((__m256i *)(in_ptr + 0 * stride_in)), _mm256_loadu_si256((__m256i *)(in_ptr + 1 * stride_in)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(in_ptr + 2 * stride_in)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(in_ptr + 3 * stride_in)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(in_ptr + 4 * stride_in)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(in_ptr + 5 * stride_in)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(in_ptr + 6 * stride_in)));
        sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(in_ptr + 7 * stride_in)));

        in_ptr += 8 * stride_in;
        sum1 = _mm256_add_epi16(_mm256_loadu_si256((__m256i *)(in_ptr + 0 * stride_in)), _mm256_loadu_si256((__m256i *)(in_ptr + 1 * stride_in)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(in_ptr + 2 * stride_in)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(in_ptr + 3 * stride_in)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(in_ptr + 4 * stride_in)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(in_ptr + 5 * stride_in)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(in_ptr + 6 * stride_in)));
        sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(in_ptr + 7 * stride_in)));

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
        int16_t *inPtrTemp = in_ptr;

        sum0 = sum1 = sum2 = sum3 = _mm256_setzero_si256();
        for (row_index = 0; row_index < size; row_index += 2) { // Parse every two rows
            sum0 = _mm256_add_epi16(sum0, _mm256_loadu_si256((__m256i *)(inPtrTemp)));
            sum1 = _mm256_add_epi16(sum1, _mm256_loadu_si256((__m256i *)(inPtrTemp + 16)));
            inPtrTemp += stride_in;
            sum2 = _mm256_add_epi16(sum2, _mm256_loadu_si256((__m256i *)(inPtrTemp)));
            sum3 = _mm256_add_epi16(sum3, _mm256_loadu_si256((__m256i *)(inPtrTemp + 16)));
            inPtrTemp += stride_in;

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

void memset16bit_block_avx2_intrin(
    int16_t * in_ptr,
    uint32_t   stride_in,
    uint32_t   size,
    int16_t   value
)
{


    if (size == 4) {

        __m128i line = _mm_set1_epi16(value);

        _mm_storel_epi64((__m128i *)(in_ptr + 0 * stride_in), line);
        _mm_storel_epi64((__m128i *)(in_ptr + 1 * stride_in), line);
        _mm_storel_epi64((__m128i *)(in_ptr + 2 * stride_in), line);
        _mm_storel_epi64((__m128i *)(in_ptr + 3 * stride_in), line);

    }
    else if (size == 8) {

        __m128i line = _mm_set1_epi16(value);

        _mm_storeu_si128((__m128i *)(in_ptr + 0 * stride_in), line);
        _mm_storeu_si128((__m128i *)(in_ptr + 1 * stride_in), line);
        _mm_storeu_si128((__m128i *)(in_ptr + 2 * stride_in), line);
        _mm_storeu_si128((__m128i *)(in_ptr + 3 * stride_in), line);
        _mm_storeu_si128((__m128i *)(in_ptr + 4 * stride_in), line);
        _mm_storeu_si128((__m128i *)(in_ptr + 5 * stride_in), line);
        _mm_storeu_si128((__m128i *)(in_ptr + 6 * stride_in), line);
        _mm_storeu_si128((__m128i *)(in_ptr + 7 * stride_in), line);

    }
    else if (size == 16) {

        __m256i line = _mm256_set1_epi16(value);

        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in), line);

        in_ptr += 8 * stride_in;

        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in), line);


    }
    else if (size == 32) {

        __m256i line = _mm256_set1_epi16(value);

        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in + 16), line);

        in_ptr += 8 * stride_in;

        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in + 16), line);

        in_ptr += 8 * stride_in;

        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in + 16), line);

        in_ptr += 8 * stride_in;

        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 0 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 1 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 2 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 3 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 4 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 5 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 6 * stride_in + 16), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in), line);
        _mm256_storeu_si256((__m256i *)(in_ptr + 7 * stride_in + 16), line);

    }


    else {
        printf("\n add the rest \n");
    }

}


void unpack_avg_safe_sub_avx2_intrin(
    uint16_t *ref16_l0,
    uint32_t  ref_l0_stride,
    uint16_t *ref16_l1,
    uint32_t  ref_l1_stride,
    uint8_t  *dst_ptr,
    uint32_t  dst_stride,
    EbBool  sub_pred,
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

        if (sub_pred) {
            ref16_l0 -= (ref_l0_stride >> 1);
            ref16_l1 -= (ref_l1_stride >> 1);
            dst_ptr -= (dst_stride >> 1);
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

            _mm_store_si128((__m128i*) dst_ptr, avg8_0_U8);

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

            _mm_store_si128((__m128i*)(dst_ptr + dst_stride), avg8_2_U8);

            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;

        }

        if (sub_pred) {
            ref16_l0 -= (ref_l0_stride >> 1);
            ref16_l1 -= (ref_l1_stride >> 1);
            dst_ptr -= (dst_stride >> 1);
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
            _mm_store_si128((__m128i*) dst_ptr, avg8_0_U8);
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
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 16));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 16));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);

            _mm256_storeu_si256((__m256i *)(dst_ptr), avg8b_32_0);

            //--------
            //Line Two
            //--------
              //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*)(ref16_l0 + ref_l0_stride));
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l0 + ref_l0_stride + 16));

            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));

            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*)(ref16_l1 + ref_l1_stride));
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l1 + ref_l1_stride + 16));

            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);

            _mm256_storeu_si256((__m256i *)(dst_ptr + dst_stride), avg8b_32_0);


            dst_ptr += 2 * dst_stride;
            ref16_l0 += 2 * ref_l0_stride;
            ref16_l1 += 2 * ref_l1_stride;

        }

        if (sub_pred) {
            ref16_l0 -= (ref_l0_stride >> 1);
            ref16_l1 -= (ref_l1_stride >> 1);
            dst_ptr -= (dst_stride >> 1);
            //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 16));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 16));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);
            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);
            _mm256_storeu_si256((__m256i *)(dst_ptr), avg8b_32_0);

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
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 48));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 48));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);
            avg8b_32_1 = _mm256_avg_epu8(data8b_32_1_L0, data8b_32_1_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);
            avg8b_32_1 = _mm256_permute4x64_epi64(avg8b_32_1, 216);

            _mm256_storeu_si256((__m256i *)(dst_ptr), avg8b_32_0);
            _mm256_storeu_si256((__m256i *)(dst_ptr + 32), avg8b_32_1);

            dst_ptr += dst_stride;
            ref16_l0 += ref_l0_stride;
            ref16_l1 += ref_l1_stride;
        }

        if (sub_pred) {
            ref16_l0 -= (ref_l0_stride >> 1);
            ref16_l1 -= (ref_l1_stride >> 1);
            dst_ptr -= (dst_stride >> 1);
            //List0
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l0);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16_l0 + 48));
            data8b_32_0_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L0 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));
            //List1
            inVal16b_0 = _mm256_loadu_si256((__m256i*) ref16_l1);
            inVal16b_1 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 16));
            inVal16b_2 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 32));
            inVal16b_3 = _mm256_loadu_si256((__m256i*)(ref16_l1 + 48));
            data8b_32_0_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_0, 2), _mm256_srli_epi16(inVal16b_1, 2));
            data8b_32_1_L1 = _mm256_packus_epi16(_mm256_srli_epi16(inVal16b_2, 2), _mm256_srli_epi16(inVal16b_3, 2));


            //Avg
            avg8b_32_0 = _mm256_avg_epu8(data8b_32_0_L0, data8b_32_0_L1);
            avg8b_32_1 = _mm256_avg_epu8(data8b_32_1_L0, data8b_32_1_L1);

            avg8b_32_0 = _mm256_permute4x64_epi64(avg8b_32_0, 216);
            avg8b_32_1 = _mm256_permute4x64_epi64(avg8b_32_1, 216);

            _mm256_storeu_si256((__m256i *)(dst_ptr), avg8b_32_0);
            _mm256_storeu_si256((__m256i *)(dst_ptr + 32), avg8b_32_1);
        }


    }


    return;
}


void picture_addition_kernel4x4_av1_sse2_intrin(
    uint8_t  *pred_ptr,
    uint32_t  pred_stride,
    int32_t *residual_ptr,
    uint32_t  residual_stride,
    uint8_t  *recon_ptr,
    uint32_t  recon_stride,
    uint32_t  width,
    uint32_t  height,
    int32_t     bd)
{

    __m128i predReg, xmm0, recon_0_7, resReg;
    uint32_t y;
    xmm0 = _mm_setzero_si128();

    for (y = 0; y < 4; ++y) {
        predReg = _mm_cvtsi32_si128(*(uint32_t *)pred_ptr);
        predReg = _mm_unpacklo_epi8(predReg, xmm0);
        predReg = _mm_unpacklo_epi16(predReg, xmm0);
        resReg = _mm_loadu_si128((__m128i *)residual_ptr);
        resReg = _mm_add_epi32(resReg, predReg);
        recon_0_7 = _mm_packus_epi32(resReg, xmm0);
        recon_0_7 = _mm_packus_epi16(recon_0_7, xmm0);
        *(uint32_t *)recon_ptr = _mm_cvtsi128_si32(recon_0_7);
        pred_ptr += pred_stride;
        residual_ptr += residual_stride;
        recon_ptr += recon_stride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}
void picture_addition_kernel8x8_av1_sse2_intrin(
    uint8_t  *pred_ptr,
    uint32_t  pred_stride,
    int32_t *residual_ptr,
    uint32_t  residual_stride,
    uint8_t  *recon_ptr,
    uint32_t  recon_stride,
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
        predReg_128 = _mm_cvtsi64_si128(*(uint64_t *)pred_ptr);
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
        *(uint64_t *)recon_ptr = _mm_cvtsi128_si64(recon_0_7_128);
        pred_ptr += pred_stride;
        residual_ptr += residual_stride;
        recon_ptr += recon_stride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}
void picture_addition_kernel16x16_av1_sse2_intrin(
    uint8_t  *pred_ptr,
    uint32_t  pred_stride,
    int32_t *residual_ptr,
    uint32_t  residual_stride,
    uint8_t  *recon_ptr,
    uint32_t  recon_stride,
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

        predReg_128 = _mm_loadu_si128((__m128i *)pred_ptr);
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
        _mm_storeu_si128((__m128i*) (recon_ptr), predReg_128Lo);
        pred_ptr += pred_stride;
        residual_ptr += residual_stride;
        recon_ptr += recon_stride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}

void picture_addition_kernel32x32_av1_sse2_intrin(
    uint8_t  *pred_ptr,
    uint32_t  pred_stride,
    int32_t *residual_ptr,
    uint32_t  residual_stride,
    uint8_t  *recon_ptr,
    uint32_t  recon_stride,
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

        predReg = _mm256_loadu_si256((__m256i*)pred_ptr);
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
        _mm256_storeu_si256((__m256i*)recon_ptr, recon_0_7);

        pred_ptr += pred_stride;
        residual_ptr += residual_stride;
        recon_ptr += recon_stride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}
void picture_addition_kernel64x64_av1_sse2_intrin(
    uint8_t  *pred_ptr,
    uint32_t  pred_stride,
    int32_t *residual_ptr,
    uint32_t  residual_stride,
    uint8_t  *recon_ptr,
    uint32_t  recon_stride,
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

        predReg = _mm256_loadu_si256((__m256i*)pred_ptr);
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
        _mm256_storeu_si256((__m256i*)recon_ptr, recon_0_7);

        predReg = _mm256_loadu_si256((__m256i*)(pred_ptr + 32));
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
        _mm256_storeu_si256((__m256i*)(recon_ptr + 32), recon_0_7);

        pred_ptr += pred_stride;
        residual_ptr += residual_stride;
        recon_ptr += recon_stride;
    }
    (void)width;
    (void)height;
    (void)bd;

    return;
}

void full_distortion_kernel32_bits_avx2(
    int32_t  *coeff,
    uint32_t   coeff_stride,
    int32_t  *recon_coeff,
    uint32_t   recon_coeff_stride,
    uint64_t   distortion_result[DIST_CALC_TOTAL],
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t rowCount, colCount;
    __m256i sum1 = _mm256_setzero_si256();
    __m256i sum2 = _mm256_setzero_si256();
    __m128i temp1, temp2, temp3;

    rowCount = area_height;
    do
    {
        int32_t *coeffTemp = coeff;
        int32_t *reconCoeffTemp = recon_coeff;

        colCount = area_width / 4;
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

        coeff += coeff_stride;
        recon_coeff += recon_coeff_stride;
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

    _mm_storeu_si128((__m128i *)distortion_result, temp1);
}

void full_distortion_kernel_cbf_zero32_bits_avx2(
    int32_t  *coeff,
    uint32_t   coeff_stride,
    int32_t  *recon_coeff,
    uint32_t   recon_coeff_stride,
    uint64_t   distortion_result[DIST_CALC_TOTAL],
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t rowCount, colCount;
    __m256i sum = _mm256_setzero_si256();
    __m128i temp1, temp2;

    rowCount = area_height;
    do
    {
        int32_t *coeffTemp = coeff;

        colCount = area_width / 4;
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

        coeff += coeff_stride;
        recon_coeff += coeff_stride;
        rowCount -= 1;
    } while (rowCount > 0);

    temp1 = _mm256_extracti128_si256(sum, 0);
    temp2 = _mm256_extracti128_si256(sum, 1);
    temp1 = _mm_add_epi64(temp1, temp2);
    temp2 = _mm_shuffle_epi32(temp1, 0x4e);
    temp1 = _mm_add_epi64(temp1, temp2);
    _mm_storeu_si128((__m128i *)distortion_result, temp1);
    (void)recon_coeff_stride;
}

static INLINE void _mm_storeh_epi64(__m128i *const d, const __m128i s) {
    _mm_storeh_pi((__m64 *)d, _mm_castsi128_ps(s));
}

void ResidualKernel4x4_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256i in = load8bit_4x4_avx2(input, input_stride);
    const __m256i pr = load8bit_4x4_avx2(pred, pred_stride);
    const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
    const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
    const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
    const __m128i re_01 = _mm256_extracti128_si256(re_lo, 0);
    const __m128i re_23 = _mm256_extracti128_si256(re_lo, 1);
    (void)area_width;
    (void)area_height;

    _mm_storel_epi64((__m128i*)(residual + 0 * residual_stride), re_01);
    _mm_storeh_epi64((__m128i*)(residual + 1 * residual_stride), re_01);
    _mm_storel_epi64((__m128i*)(residual + 2 * residual_stride), re_23);
    _mm_storeh_epi64((__m128i*)(residual + 3 * residual_stride), re_23);
}
void ResidualKernel4x8_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 8; y += 4) {
        const __m256i in = load8bit_4x4_avx2(input, input_stride);
        const __m256i pr = load8bit_4x4_avx2(pred, pred_stride);
        const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m128i re_01 = _mm256_extracti128_si256(re_lo, 0);
        const __m128i re_23 = _mm256_extracti128_si256(re_lo, 1);

        _mm_storel_epi64((__m128i*)(residual + 0 * residual_stride), re_01);
        _mm_storeh_epi64((__m128i*)(residual + 1 * residual_stride), re_01);
        _mm_storel_epi64((__m128i*)(residual + 2 * residual_stride), re_23);
        _mm_storeh_epi64((__m128i*)(residual + 3 * residual_stride), re_23);
        input += 4 * input_stride;
        pred += 4 * pred_stride;
        residual += 4 * residual_stride;
    }
}
void ResidualKernel4x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 16; y += 4) {
        const __m256i in = load8bit_4x4_avx2(input, input_stride);
        const __m256i pr = load8bit_4x4_avx2(pred, pred_stride);
        const __m256i in_lo = _mm256_unpacklo_epi8(in, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m128i re_01 = _mm256_extracti128_si256(re_lo, 0);
        const __m128i re_23 = _mm256_extracti128_si256(re_lo, 1);

        _mm_storel_epi64((__m128i*)(residual + 0 * residual_stride), re_01);
        _mm_storeh_epi64((__m128i*)(residual + 1 * residual_stride), re_01);
        _mm_storel_epi64((__m128i*)(residual + 2 * residual_stride), re_23);
        _mm_storeh_epi64((__m128i*)(residual + 3 * residual_stride), re_23);
        input += 4 * input_stride;
        pred += 4 * pred_stride;
        residual += 4 * residual_stride;
    }
}
void ResidualKernel8x32_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 32; y += 4) {
        const __m256i in = load8bit_8x4_avx2(input, input_stride);
        const __m256i pr = load8bit_8x4_avx2(pred, pred_stride);
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

        _mm_storeu_si128((__m128i*)(residual + 0 * residual_stride), re_0);
        _mm_storeu_si128((__m128i*)(residual + 1 * residual_stride), re_1);
        _mm_storeu_si128((__m128i*)(residual + 2 * residual_stride), re_2);
        _mm_storeu_si128((__m128i*)(residual + 3 * residual_stride), re_3);
        input += 4 * input_stride;
        pred += 4 * pred_stride;
        residual += 4 * residual_stride;
    }
}
void ResidualKernel8x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 16; y += 4) {
        const __m256i in = load8bit_8x4_avx2(input, input_stride);
        const __m256i pr = load8bit_8x4_avx2(pred, pred_stride);
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

        _mm_storeu_si128((__m128i*)(residual + 0 * residual_stride), re_0);
        _mm_storeu_si128((__m128i*)(residual + 1 * residual_stride), re_1);
        _mm_storeu_si128((__m128i*)(residual + 2 * residual_stride), re_2);
        _mm_storeu_si128((__m128i*)(residual + 3 * residual_stride), re_3);
        input += 4 * input_stride;
        pred += 4 * pred_stride;
        residual += 4 * residual_stride;
    }
}

void ResidualKernel8x8_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 8; y += 4) {
        const __m256i in = load8bit_8x4_avx2(input, input_stride);
        const __m256i pr = load8bit_8x4_avx2(pred, pred_stride);
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

        _mm_storeu_si128((__m128i*)(residual + 0 * residual_stride), re_0);
        _mm_storeu_si128((__m128i*)(residual + 1 * residual_stride), re_1);
        _mm_storeu_si128((__m128i*)(residual + 2 * residual_stride), re_2);
        _mm_storeu_si128((__m128i*)(residual + 3 * residual_stride), re_3);
        input += 4 * input_stride;
        pred += 4 * pred_stride;
        residual += 4 * residual_stride;
    }
}

void ResidualKernel8x4_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    //for (y = 0; y < 8; y += 4) {
        const __m256i in = load8bit_8x4_avx2(input, input_stride);
        const __m256i pr = load8bit_8x4_avx2(pred, pred_stride);
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

        _mm_storeu_si128((__m128i*)(residual + 0 * residual_stride), re_0);
        _mm_storeu_si128((__m128i*)(residual + 1 * residual_stride), re_1);
        _mm_storeu_si128((__m128i*)(residual + 2 * residual_stride), re_2);
        _mm_storeu_si128((__m128i*)(residual + 3 * residual_stride), re_3);
        //input += 4 * input_stride;
        //pred += 4 * pred_stride;
        //residual += 4 * residual_stride;
    //}
}

void ResidualKernel16x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 16; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, input_stride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, pred_stride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residual_stride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residual_stride), re_hi);
        input += 2 * input_stride;
        pred += 2 * pred_stride;
        residual += 2 * residual_stride;
    }
}

void ResidualKernel16x4_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 4; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, input_stride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, pred_stride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residual_stride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residual_stride), re_hi);
        input += 2 * input_stride;
        pred += 2 * pred_stride;
        residual += 2 * residual_stride;
    }
}

void ResidualKernel16x8_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 8; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, input_stride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, pred_stride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residual_stride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residual_stride), re_hi);
        input += 2 * input_stride;
        pred += 2 * pred_stride;
        residual += 2 * residual_stride;
    }
}

void ResidualKernel16x32_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 32; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, input_stride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, pred_stride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residual_stride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residual_stride), re_hi);
        input += 2 * input_stride;
        pred += 2 * pred_stride;
        residual += 2 * residual_stride;
    }
}

void ResidualKernel16x64_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 64; y += 2) {
        const __m256i in0 = load8bit_16x2_unaligned_avx2(input, input_stride);
        const __m256i pr0 = load8bit_16x2_unaligned_avx2(pred, pred_stride);
        const __m256i in1 = _mm256_permute4x64_epi64(in0, 0xD8);
        const __m256i pr1 = _mm256_permute4x64_epi64(pr0, 0xD8);
        const __m256i in_lo = _mm256_unpacklo_epi8(in1, zero);
        const __m256i in_hi = _mm256_unpackhi_epi8(in1, zero);
        const __m256i pr_lo = _mm256_unpacklo_epi8(pr1, zero);
        const __m256i pr_hi = _mm256_unpackhi_epi8(pr1, zero);
        const __m256i re_lo = _mm256_sub_epi16(in_lo, pr_lo);
        const __m256i re_hi = _mm256_sub_epi16(in_hi, pr_hi);

        _mm256_storeu_si256((__m256i*)(residual + 0 * residual_stride), re_lo);
        _mm256_storeu_si256((__m256i*)(residual + 1 * residual_stride), re_hi);
        input += 2 * input_stride;
        pred += 2 * pred_stride;
        residual += 2 * residual_stride;
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
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 8; ++y) {
        ResidualKernel32_AVX2(input, pred, residual);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}

void ResidualKernel32x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 16; ++y) {
        ResidualKernel32_AVX2(input, pred, residual);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}

void ResidualKernel32x32_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 32; ++y) {
        ResidualKernel32_AVX2(input, pred, residual);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}


void ResidualKernel32x64_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 64; ++y) {
        ResidualKernel32_AVX2(input, pred, residual);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}

void ResidualKernel64x16_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 16; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}

void ResidualKernel64x32_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 32; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}

void ResidualKernel64x64_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 64; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}

void ResidualKernel64x128_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 128; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}


void ResidualKernel128x128_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 128; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        ResidualKernel32_AVX2(input + 0x40, pred + 0x40, residual + 0x40);
        ResidualKernel32_AVX2(input + 0x60, pred + 0x60, residual + 0x60);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}

void ResidualKernel128x64_AVX2_INTRIN(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t y;
    (void)area_width;
    (void)area_height;

    for (y = 0; y < 64; ++y) {
        ResidualKernel32_AVX2(input + 0x00, pred + 0x00, residual + 0x00);
        ResidualKernel32_AVX2(input + 0x20, pred + 0x20, residual + 0x20);
        ResidualKernel32_AVX2(input + 0x40, pred + 0x40, residual + 0x40);
        ResidualKernel32_AVX2(input + 0x60, pred + 0x60, residual + 0x60);
        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
    }
}
void ResidualKernel_avx2(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    if (area_width == area_height) {
        switch (area_width) {
        case 4:
            ResidualKernel4x4_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
            break;
        case 8:
            ResidualKernel8x8_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
            break;
        case 16:
            ResidualKernel16x16_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
            break;
        case 32:
            ResidualKernel32x32_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
            break;
        case 64:
            ResidualKernel64x64_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
            break;
        case 128:
            ResidualKernel128x128_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
            break;
        }
    }
    else {
        if (area_width < area_height) {
            if (area_width + area_width == area_height) {
                switch (area_width) {
                case 4:
                    ResidualKernel4x8_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 8:
                    ResidualKernel8x16_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 16:
                    ResidualKernel16x32_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 32:
                    ResidualKernel32x64_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 64:
                    ResidualKernel64x128_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                }
            }
            else {
                switch (area_width) {
                case 4:
                    ResidualKernel4x16_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 8:
                    ResidualKernel8x32_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 16:
                    ResidualKernel16x64_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                }
            }
        }
        else {
            if (area_height + area_height == area_width) {
                switch (area_height) {
                case 4:
                    ResidualKernel8x4_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 8:
                    ResidualKernel16x8_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 16:
                    ResidualKernel32x16_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 32:
                    ResidualKernel64x32_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 64:
                    ResidualKernel128x64_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                }
            }
            else {
                switch (area_height) {
                case 4:
                    ResidualKernel16x4_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 8:
                    ResidualKernel32x8_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                case 16:
                    ResidualKernel64x16_AVX2_INTRIN(input, input_stride, pred, pred_stride, residual, residual_stride, area_width, area_height);
                    break;
                }
            }
        }
    }
}
