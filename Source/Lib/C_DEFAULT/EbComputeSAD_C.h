/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbComputeSAD_C_h
#define EbComputeSAD_C_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    uint32_t FastLoop_NxMSadKernel(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  refStride,                      // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    uint32_t CombinedAveragingSAD(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

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
        int16_t search_area_height);

#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_C_h