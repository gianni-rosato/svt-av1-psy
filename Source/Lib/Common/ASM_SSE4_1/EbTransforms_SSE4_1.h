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
    extern EB_ALIGN(16) const int16_t transform_asm_const_sse4_1[1632];

    extern void transform8x8_sse4_1_intrin(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void pfreq_transform8x8_sse4_1_intrin(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void pfreq_n4_transform8x8_sse4_1_intrin(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

#ifdef __cplusplus
}
#endif
#endif // EbTransforms_SSE4_1_h