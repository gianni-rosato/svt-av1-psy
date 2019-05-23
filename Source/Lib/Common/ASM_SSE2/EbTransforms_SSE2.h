/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbTransforms_SSE2_h
#define EbTransforms_SSE2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    extern EB_ALIGN(16) const int16_t dst_transform_asm_const_sse2[88];
    extern EB_ALIGN(16) const int16_t inv_transform_asm_const_sse2[1512];
    extern EB_ALIGN(16) const int16_t inv_dst_transform_asm_const_sse2[72];

    extern void transform4x4_sse2_intrin(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void dst_transform4x4_sse2_intrin(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void transform8x8_sse2_intrin(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void pfreq_transform8x8_sse2_intrin(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void pfreq_n4_transform8x8_sse2_intrin(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void transform16x16_sse2(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void transform32x32_sse2(
        int16_t        *residual,
        const uint32_t  src_stride,
        int16_t        *transform_coefficients,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void inv_transform8x8_sse2_intrin(
        int16_t        *transform_coefficients,
        const uint32_t  src_stride,
        int16_t        *residual,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void inv_transform4x4_sse2_intrin(
        int16_t        *transform_coefficients,
        const uint32_t  src_stride,
        int16_t        *residual,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    extern void inv_dst_transform4x4_sse2_intrin(
        int16_t        *transform_coefficients,
        const uint32_t  src_stride,
        int16_t        *residual,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    void pfreq_transform32x32_sse2(
        int16_t       *src,
        const uint32_t src_stride,
        int16_t       *dst,
        const uint32_t dst_stride,
        int16_t       *intermediate,
        uint32_t       addshift);

    void pfreq_transform16x16_sse2(
        int16_t       *src,
        const uint32_t src_stride,
        int16_t       *dst,
        const uint32_t dst_stride,
        int16_t       *intermediate,
        uint32_t       addshift);

    void pfreq_n4_transform32x32_sse2(
        int16_t       *src,
        const uint32_t src_stride,
        int16_t       *dst,
        const uint32_t dst_stride,
        int16_t       *intermediate,
        uint32_t       addshift);

    void pfreq_n4_transform16x16_sse2(
        int16_t       *src,
        const uint32_t src_stride,
        int16_t       *dst,
        const uint32_t dst_stride,
        int16_t       *intermediate,
        uint32_t       addshift);

    void estimate_inv_transform32x32_sse2(
        int16_t        *transform_coefficients,
        const uint32_t  src_stride,
        int16_t        *residual,
        const uint32_t  dst_stride,
        int16_t        *transform_inner_array_ptr,
        uint32_t        bit_increment);

    void estimate_inv_transform16x16_sse2(
        int16_t  *src,
        uint32_t  src_stride,
        int16_t  *dst,
        uint32_t  dst_stride,
        int16_t  *intermediate,
        uint32_t  addshift);

#ifdef __cplusplus
}
#endif
#endif // EbTransforms_SSE2_h

