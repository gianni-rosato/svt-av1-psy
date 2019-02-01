/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbComputeSAD_C.h"
#include "EbUtility.h"

/*******************************************
* CombinedAveragingSAD
*
*******************************************/
uint32_t CombinedAveragingSAD(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1Stride,
    uint8_t  *ref2,
    uint32_t  ref2Stride,
    uint32_t  height,
    uint32_t  width)
{
    uint32_t x, y;
    uint32_t sad = 0;
    uint8_t avgpel;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            avgpel = (ref1[x] + ref2[x] + 1) >> 1;
            sad += EB_ABS_DIFF(src[x], avgpel);
        }
        src += src_stride;
        ref1 += ref1Stride;
        ref2 += ref2Stride;
    }

    return sad;
}

/*******************************************
*   returns NxM Sum of Absolute Differences
Note: moved from picture operators.
keep this function here for profiling
issues.
*******************************************/
uint32_t FastLoop_NxMSadKernel(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  refStride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width)                          // input parameter, block width (N)
{
    uint32_t x, y;
    uint32_t sad = 0;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            sad += EB_ABS_DIFF(src[x], ref[x]);
        }
        src += src_stride;
        ref += refStride;
    }

    return sad;
}

void SadLoopKernel(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  refStride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width,                          // input parameter, block width (N)
    uint64_t *bestSad,
    int16_t *xSearchCenter,
    int16_t *ySearchCenter,
    uint32_t  srcStrideRaw,                   // input parameter, source stride (no line skipping)
    int16_t search_area_width,
    int16_t search_area_height)
{
    int16_t xSearchIndex;
    int16_t ySearchIndex;

    *bestSad = 0xffffff;

    for (ySearchIndex = 0; ySearchIndex < search_area_height; ySearchIndex++)
    {
        for (xSearchIndex = 0; xSearchIndex < search_area_width; xSearchIndex++)
        {
            uint32_t x, y;
            uint32_t sad = 0;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++)
                {
                    sad += EB_ABS_DIFF(src[y*src_stride + x], ref[xSearchIndex + y * refStride + x]);
                }

            }

            // Update results
            if (sad < *bestSad)
            {
                *bestSad = sad;
                *xSearchCenter = xSearchIndex;
                *ySearchCenter = ySearchIndex;
            }
        }

        ref += srcStrideRaw;
    }

    return;
}