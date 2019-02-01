/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbIntraPrediction_SSSE3_h
#define EbIntraPrediction_SSSE3_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void IntraModeAngular_Vertical_Kernel_SSSE3_INTRIN(
    uint32_t            size,
    uint8_t            *refSampMain,
    uint8_t            *prediction_ptr,
    uint32_t            predictionBufferStride,
    const EbBool     skip,
    int32_t            intraPredAngle);

extern void IntraModeAngular_Horizontal_Kernel_SSSE3_INTRIN(
    uint32_t            size,
    uint8_t            *refSampMain,
    uint8_t            *prediction_ptr,
    uint32_t            predictionBufferStride,
    const EbBool     skip,
    int32_t            intraPredAngle);


#ifdef __cplusplus
}
#endif
#endif // EbIntraPrediction_SSSE3_h