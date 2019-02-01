/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbComputeSAD_SSE2_h
#define EbComputeSAD_SSE2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    extern uint32_t Compute4xMSadSub_SSE2_INTRIN(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  refStride,                      // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    extern uint32_t CombinedAveraging4xMSAD_SSE2_INTRIN(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);


#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_SSE2_h