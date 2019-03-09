/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbComputeSAD_C.h"
#include "EbUtility.h"

/*******************************************
* combined_averaging_sad
*
*******************************************/
uint32_t combined_averaging_sad(
    uint8_t  *src,
    uint32_t  src_stride,
    uint8_t  *ref1,
    uint32_t  ref1_stride,
    uint8_t  *ref2,
    uint32_t  ref2_stride,
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
        ref1 += ref1_stride;
        ref2 += ref2_stride;
    }

    return sad;
}

/*******************************************
*   returns NxM Sum of Absolute Differences
Note: moved from picture operators.
keep this function here for profiling
issues.
*******************************************/
uint32_t fast_loop_nx_m_sad_kernel(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  ref_stride,                      // input parameter, reference stride
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
        ref += ref_stride;
    }

    return sad;
}

void sad_loop_kernel(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  ref_stride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width,                          // input parameter, block width (N)
    uint64_t *best_sad,
    int16_t *x_search_center,
    int16_t *y_search_center,
    uint32_t  src_stride_raw,                   // input parameter, source stride (no line skipping)
    int16_t search_area_width,
    int16_t search_area_height)
{
    int16_t xSearchIndex;
    int16_t ySearchIndex;

    *best_sad = 0xffffff;

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
                    sad += EB_ABS_DIFF(src[y*src_stride + x], ref[xSearchIndex + y * ref_stride + x]);
                }

            }

            // Update results
            if (sad < *best_sad)
            {
                *best_sad = sad;
                *x_search_center = xSearchIndex;
                *y_search_center = ySearchIndex;
            }
        }

        ref += src_stride_raw;
    }

    return;
}