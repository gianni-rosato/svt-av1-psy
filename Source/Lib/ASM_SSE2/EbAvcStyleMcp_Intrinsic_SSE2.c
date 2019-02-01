/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/


#include "EbAvcStyleMcp_SSE2.h"
#include "EbMcp_SSE2.h" // THIS SHOULD BE _SSE2 in the future
#include "emmintrin.h"
void AvcStyleCopy_SSE2(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    EbByte               tempBuf,
    EbBool               skip,
    uint32_t                fracPos)
{
    (void)tempBuf;
    (void)fracPos;
    if (skip) {
        //do the last row too.
        EB_MEMCPY(dst + (pu_height - 1)*dst_stride, ref_pic + (pu_height - 1)*src_stride, pu_width);

        src_stride <<= 1;
        dst_stride <<= 1;
        pu_height >>= 1;
    }

    PictureCopyKernel_SSE2(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height);
}

//This function should be removed and replace by AvcStyleCopy_SSE2

void PictureAverageKernel_SSE2_INTRIN(
    EbByte                  src0,
    uint32_t                   src0Stride,
    EbByte                  src1,
    uint32_t                   src1Stride,
    EbByte                  dst,
    uint32_t                   dst_stride,
    uint32_t                   areaWidth,
    uint32_t                   areaHeight)
{
    __m128i xmm_avg1, xmm_avg2, xmm_avg3, xmm_avg4, xmm_avg5, xmm_avg6, xmm_avg7, xmm_avg8;
    uint32_t y;

    if (areaWidth > 16)
    {
        if (areaWidth == 24)
        {
            for (y = 0; y < areaHeight; y += 2) {
                xmm_avg1 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)src0), _mm_loadu_si128((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)(src0 + 16)), _mm_loadl_epi64((__m128i*)(src1 + 16)));
                xmm_avg3 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride)), _mm_loadu_si128((__m128i*)(src1 + src1Stride)));
                xmm_avg4 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)(src0 + src0Stride + 16)), _mm_loadl_epi64((__m128i*)(src1 + src1Stride + 16)));

                _mm_storeu_si128((__m128i*) dst, xmm_avg1);
                _mm_storel_epi64((__m128i*) (dst + 16), xmm_avg2);
                _mm_storeu_si128((__m128i*) (dst + dst_stride), xmm_avg3);
                _mm_storel_epi64((__m128i*) (dst + dst_stride + 16), xmm_avg4);

                src0 += src0Stride << 1;
                src1 += src1Stride << 1;
                dst += dst_stride << 1;
            }
        }
        else if (areaWidth == 32)
        {
            for (y = 0; y < areaHeight; y += 2) {

                xmm_avg1 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)src0), _mm_loadu_si128((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 16)), _mm_loadu_si128((__m128i*)(src1 + 16)));
                xmm_avg3 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride)), _mm_loadu_si128((__m128i*)(src1 + src1Stride)));
                xmm_avg4 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 16)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 16)));

                _mm_storeu_si128((__m128i*) dst, xmm_avg1);
                _mm_storeu_si128((__m128i*) (dst + 16), xmm_avg2);
                _mm_storeu_si128((__m128i*) (dst + dst_stride), xmm_avg3);
                _mm_storeu_si128((__m128i*) (dst + dst_stride + 16), xmm_avg4);

                src0 += src0Stride << 1;
                src1 += src1Stride << 1;
                dst += dst_stride << 1;
            }
        }
        else if (areaWidth == 48)
        {
            for (y = 0; y < areaHeight; y += 2) {
                xmm_avg1 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)src0), _mm_loadu_si128((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 16)), _mm_loadu_si128((__m128i*)(src1 + 16)));
                xmm_avg3 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 32)), _mm_loadu_si128((__m128i*)(src1 + 32)));

                xmm_avg4 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride)), _mm_loadu_si128((__m128i*)(src1 + src1Stride)));
                xmm_avg5 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 16)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 16)));
                xmm_avg6 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 32)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 32)));

                _mm_storeu_si128((__m128i*) dst, xmm_avg1);
                _mm_storeu_si128((__m128i*) (dst + 16), xmm_avg2);
                _mm_storeu_si128((__m128i*) (dst + 32), xmm_avg3);
                _mm_storeu_si128((__m128i*) (dst + dst_stride), xmm_avg4);
                _mm_storeu_si128((__m128i*) (dst + dst_stride + 16), xmm_avg5);
                _mm_storeu_si128((__m128i*) (dst + dst_stride + 32), xmm_avg6);

                src0 += src0Stride << 1;
                src1 += src1Stride << 1;
                dst += dst_stride << 1;

            }
        }
        else
        {
            for (y = 0; y < areaHeight; y += 2) {
                xmm_avg1 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)src0), _mm_loadu_si128((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 16)), _mm_loadu_si128((__m128i*)(src1 + 16)));
                xmm_avg3 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 32)), _mm_loadu_si128((__m128i*)(src1 + 32)));
                xmm_avg4 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 48)), _mm_loadu_si128((__m128i*)(src1 + 48)));

                xmm_avg5 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride)), _mm_loadu_si128((__m128i*)(src1 + src1Stride)));
                xmm_avg6 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 16)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 16)));
                xmm_avg7 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 32)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 32)));
                xmm_avg8 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 48)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 48)));

                _mm_storeu_si128((__m128i*) dst, xmm_avg1);
                _mm_storeu_si128((__m128i*) (dst + 16), xmm_avg2);
                _mm_storeu_si128((__m128i*) (dst + 32), xmm_avg3);
                _mm_storeu_si128((__m128i*) (dst + 48), xmm_avg4);

                _mm_storeu_si128((__m128i*) (dst + dst_stride), xmm_avg5);
                _mm_storeu_si128((__m128i*) (dst + dst_stride + 16), xmm_avg6);
                _mm_storeu_si128((__m128i*) (dst + dst_stride + 32), xmm_avg7);
                _mm_storeu_si128((__m128i*) (dst + dst_stride + 48), xmm_avg8);

                src0 += src0Stride << 1;
                src1 += src1Stride << 1;
                dst += dst_stride << 1;
            }
        }
    }
    else
    {
        if (areaWidth == 16)
        {
            for (y = 0; y < areaHeight; y += 2) {
                xmm_avg1 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)src0), _mm_loadu_si128((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride)), _mm_loadu_si128((__m128i*)(src1 + src1Stride)));

                _mm_storeu_si128((__m128i*) dst, xmm_avg1);
                _mm_storeu_si128((__m128i*) (dst + dst_stride), xmm_avg2);

                src0 += src0Stride << 1;
                src1 += src1Stride << 1;
                dst += dst_stride << 1;
            }
        }
        else if (areaWidth == 4)
        {
            for (y = 0; y < areaHeight; y += 2) {

                xmm_avg1 = _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)src0), _mm_cvtsi32_si128(*(uint32_t *)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)(src0 + src0Stride)), _mm_cvtsi32_si128(*(uint32_t *)(src1 + src1Stride)));

                *(uint32_t *)dst = _mm_cvtsi128_si32(xmm_avg1);
                *(uint32_t *)(dst + dst_stride) = _mm_cvtsi128_si32(xmm_avg2);

                src0 += src0Stride << 1;
                src1 += src1Stride << 1;
                dst += dst_stride << 1;
            }
        }
        else if (areaWidth == 8)
        {
            for (y = 0; y < areaHeight; y += 2) {

                xmm_avg1 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)src0), _mm_loadl_epi64((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)(src0 + src0Stride)), _mm_loadl_epi64((__m128i*)(src1 + src1Stride)));

                _mm_storel_epi64((__m128i*) dst, xmm_avg1);
                _mm_storel_epi64((__m128i*) (dst + dst_stride), xmm_avg2);

                src0 += src0Stride << 1;
                src1 += src1Stride << 1;
                dst += dst_stride << 1;
            }
        }
        else
        {
            for (y = 0; y < areaHeight; y += 2) {

                xmm_avg1 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)src0), _mm_loadl_epi64((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)(src0 + 8)), _mm_cvtsi32_si128(*(uint32_t *)(src1 + 8)));

                xmm_avg3 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)(src0 + src0Stride)), _mm_loadl_epi64((__m128i*)(src1 + src1Stride)));
                xmm_avg4 = _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)(src0 + src0Stride + 8)), _mm_cvtsi32_si128(*(uint32_t *)(src1 + src1Stride + 8)));

                _mm_storel_epi64((__m128i*) dst, xmm_avg1);
                *(uint32_t *)(dst + 8) = _mm_cvtsi128_si32(xmm_avg2);
                _mm_storel_epi64((__m128i*) (dst + dst_stride), xmm_avg3);
                *(uint32_t *)(dst + dst_stride + 8) = _mm_cvtsi128_si32(xmm_avg4);

                src0 += src0Stride << 1;
                src1 += src1Stride << 1;
                dst += dst_stride << 1;
            }
        }
    }
}
void PictureAverageKernel1Line_SSE2_INTRIN(
    EbByte                  src0,
    EbByte                  src1,
    EbByte                  dst,
    uint32_t                   areaWidth)
{
    __m128i xmm_avg1, xmm_avg2, xmm_avg3, xmm_avg4;

    if (areaWidth > 16)
    {
        if (areaWidth == 32)
        {
            //for (y = 0; y < areaHeight; y += 2)
            {

                xmm_avg1 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)src0), _mm_loadu_si128((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 16)), _mm_loadu_si128((__m128i*)(src1 + 16)));
                //xmm_avg3 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride)), _mm_loadu_si128((__m128i*)(src1 + src1Stride)));
                //xmm_avg4 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 16)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 16)));

                _mm_storeu_si128((__m128i*) dst, xmm_avg1);
                _mm_storeu_si128((__m128i*) (dst + 16), xmm_avg2);
                //_mm_storeu_si128((__m128i*) (dst + dst_stride), xmm_avg3);
                //_mm_storeu_si128((__m128i*) (dst + dst_stride + 16), xmm_avg4);

                //src0 += src0Stride << 1;
                //src1 += src1Stride << 1;
                //dst += dst_stride << 1;
            }
        }
        else
        {
            //for (y = 0; y < areaHeight; y += 2)
            {
                xmm_avg1 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)src0), _mm_loadu_si128((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 16)), _mm_loadu_si128((__m128i*)(src1 + 16)));
                xmm_avg3 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 32)), _mm_loadu_si128((__m128i*)(src1 + 32)));
                xmm_avg4 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + 48)), _mm_loadu_si128((__m128i*)(src1 + 48)));

                //xmm_avg5 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride)), _mm_loadu_si128((__m128i*)(src1 + src1Stride)));
                //xmm_avg6 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 16)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 16)));
                //xmm_avg7 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 32)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 32)));
                //xmm_avg8 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride + 48)), _mm_loadu_si128((__m128i*)(src1 + src1Stride + 48)));

                _mm_storeu_si128((__m128i*) dst, xmm_avg1);
                _mm_storeu_si128((__m128i*) (dst + 16), xmm_avg2);
                _mm_storeu_si128((__m128i*) (dst + 32), xmm_avg3);
                _mm_storeu_si128((__m128i*) (dst + 48), xmm_avg4);

                //_mm_storeu_si128((__m128i*) (dst + dst_stride), xmm_avg5);
                //_mm_storeu_si128((__m128i*) (dst + dst_stride + 16), xmm_avg6);
                //_mm_storeu_si128((__m128i*) (dst + dst_stride + 32), xmm_avg7);
                //_mm_storeu_si128((__m128i*) (dst + dst_stride + 48), xmm_avg8);

                //src0 += src0Stride << 1;
                //src1 += src1Stride << 1;
                //dst += dst_stride << 1;
            }
        }
    }
    else
    {
        if (areaWidth == 16)
        {
            //for (y = 0; y < areaHeight; y += 2)
            {
                xmm_avg1 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)src0), _mm_loadu_si128((__m128i*)src1));
                //xmm_avg2 = _mm_avg_epu8(_mm_loadu_si128((__m128i*)(src0 + src0Stride)), _mm_loadu_si128((__m128i*)(src1 + src1Stride)));

                _mm_storeu_si128((__m128i*) dst, xmm_avg1);
                //_mm_storeu_si128((__m128i*) (dst + dst_stride), xmm_avg2);

                //src0 += src0Stride << 1;
                //src1 += src1Stride << 1;
                //dst += dst_stride << 1;
            }
        }
        else if (areaWidth == 4)
        {
            //for (y = 0; y < areaHeight; y += 2)
            {

                xmm_avg1 = _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)src0), _mm_cvtsi32_si128(*(uint32_t *)src1));
                //xmm_avg2 = _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)(src0 + src0Stride)), _mm_cvtsi32_si128(*(uint32_t *)(src1 + src1Stride)));

                *(uint32_t *)dst = _mm_cvtsi128_si32(xmm_avg1);
                //*(uint32_t *)(dst + dst_stride) = _mm_cvtsi128_si32(xmm_avg2);

                //src0 += src0Stride << 1;
                //src1 += src1Stride << 1;
                //dst += dst_stride << 1;
            }
        }
        else if (areaWidth == 8)
        {
            //for (y = 0; y < areaHeight; y += 2)
            {

                xmm_avg1 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)src0), _mm_loadl_epi64((__m128i*)src1));
                //xmm_avg2 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)(src0 + src0Stride)), _mm_loadl_epi64((__m128i*)(src1 + src1Stride)));

                _mm_storel_epi64((__m128i*) dst, xmm_avg1);
                //_mm_storel_epi64((__m128i*) (dst + dst_stride), xmm_avg2);

                //src0 += src0Stride << 1;
                //src1 += src1Stride << 1;
                //dst += dst_stride << 1;
            }
        }
        else
        {
            //for (y = 0; y < areaHeight; y += 2)
            {

                xmm_avg1 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)src0), _mm_loadl_epi64((__m128i*)src1));
                xmm_avg2 = _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)(src0 + 8)), _mm_cvtsi32_si128(*(uint32_t *)(src1 + 8)));

                //xmm_avg3 = _mm_avg_epu8(_mm_loadl_epi64((__m128i*)(src0 + src0Stride)), _mm_loadl_epi64((__m128i*)(src1 + src1Stride)));
                //xmm_avg4 = _mm_avg_epu8(_mm_cvtsi32_si128(*(uint32_t *)(src0 + src0Stride + 8)), _mm_cvtsi32_si128(*(uint32_t *)(src1 + src1Stride + 8)));

                _mm_storel_epi64((__m128i*) dst, xmm_avg1);
                *(uint32_t *)(dst + 8) = _mm_cvtsi128_si32(xmm_avg2);
                //_mm_storel_epi64((__m128i*) (dst + dst_stride), xmm_avg3);
                //*(uint32_t *)(dst + dst_stride + 8) = _mm_cvtsi128_si32(xmm_avg4);

                //src0 += src0Stride << 1;
                //src1 += src1Stride << 1;
                //dst += dst_stride << 1;
            }
        }
    }
}
