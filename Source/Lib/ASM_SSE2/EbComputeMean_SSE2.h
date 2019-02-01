/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbComputeMean_SS2_h
#define EbComputeMean_SS2_h
#ifdef __cplusplus
extern "C" {
#endif

#include "EbDefinitions.h"

    uint64_t ComputeSubMean8x8_SSE2_INTRIN(
        uint8_t *  input_samples,      // input parameter, input samples Ptr
        uint16_t   inputStride);
    uint64_t ComputeSubdMeanOfSquaredValues8x8_SSE2_INTRIN(
        uint8_t *  input_samples,      // input parameter, input samples Ptr
        uint16_t   inputStride);

    uint64_t ComputeMean8x8_SSE2_INTRIN(
        uint8_t *  input_samples,      // input parameter, input samples Ptr
        uint32_t   inputStride,       // input parameter, input stride
        uint32_t   inputAreaWidth,    // input parameter, input area width
        uint32_t   inputAreaHeight);   // input parameter, input area height

    uint64_t ComputeMeanOfSquaredValues8x8_SSE2_INTRIN(
        uint8_t *  input_samples,      // input parameter, input samples Ptr
        uint32_t   inputStride,       // input parameter, input stride
        uint32_t   inputAreaWidth,    // input parameter, input area width
        uint32_t   inputAreaHeight);   // input parameter, input area height


#ifdef __cplusplus
}
#endif

#endif