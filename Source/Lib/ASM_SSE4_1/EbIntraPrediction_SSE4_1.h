/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbIntraPrediction_SSE4_1_h
#define EbIntraPrediction_SSE4_1_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    extern void IntraModeDCLuma16bit_SSE4_1_INTRIN(
        const uint32_t   size,                       //input parameter, denotes the size of the current PU
        uint16_t         *refSamples,                 //input parameter, pointer to the reference samples
        uint16_t         *prediction_ptr,              //output parameter, pointer to the prediction
        const uint32_t   predictionBufferStride,     //input parameter, denotes the stride for the prediction ptr
        const EbBool  skip);                       //skip half rows

#ifdef __cplusplus
}
#endif
#endif // EbIntraPrediction_SSE4_1_h