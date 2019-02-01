/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbTransforms_SSE4_1_h
#define EbTransforms_SSE4_1_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
    extern EB_ALIGN(16) const int16_t TransformAsmConst_SSE4_1[1632];

    extern void Transform8x8_SSE4_1_INTRIN(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);
    extern void PfreqTransform8x8_SSE4_1_INTRIN(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);
    extern void PfreqN4Transform8x8_SSE4_1_INTRIN(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

#ifdef __cplusplus
}
#endif
#endif // EbTransforms_SSE4_1_h