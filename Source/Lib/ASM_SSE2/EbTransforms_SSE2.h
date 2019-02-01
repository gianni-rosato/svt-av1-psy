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

    extern EB_ALIGN(16) const int16_t DstTransformAsmConst_SSE2[88];
    extern EB_ALIGN(16) const int16_t InvTransformAsmConst_SSE2[1512];
    extern EB_ALIGN(16) const int16_t InvDstTransformAsmConst_SSE2[72];


    extern void Transform4x4_SSE2_INTRIN(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    extern void DstTransform4x4_SSE2_INTRIN(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    extern void Transform8x8_SSE2_INTRIN(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    extern void PfreqTransform8x8_SSE2_INTRIN(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    extern void PfreqN4Transform8x8_SSE2_INTRIN(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    extern void Transform16x16_SSE2(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);


    extern void Transform32x32_SSE2(
        int16_t                  *residual,
        const uint32_t             src_stride,
        int16_t                  *transformCoefficients,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    extern void InvTransform8x8_SSE2_INTRIN(
        int16_t                  *transformCoefficients,
        const uint32_t             src_stride,
        int16_t                  *residual,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    extern void InvTransform4x4_SSE2_INTRIN(
        int16_t                  *transformCoefficients,
        const uint32_t             src_stride,
        int16_t                  *residual,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    extern void InvDstTransform4x4_SSE2_INTRIN(
        int16_t                  *transformCoefficients,
        const uint32_t             src_stride,
        int16_t                  *residual,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    void PfreqTransform32x32_SSE2(
        int16_t *src,
        const uint32_t src_stride,
        int16_t *dst,
        const uint32_t dst_stride,
        int16_t *intermediate,
        uint32_t addshift);

    void PfreqTransform16x16_SSE2(
        int16_t *src,
        const uint32_t src_stride,
        int16_t *dst,
        const uint32_t dst_stride,
        int16_t *intermediate,
        uint32_t addshift);

    void PfreqN4Transform32x32_SSE2(
        int16_t *src,
        const uint32_t src_stride,
        int16_t *dst,
        const uint32_t dst_stride,
        int16_t *intermediate,
        uint32_t addshift);
    void PfreqN4Transform16x16_SSE2(
        int16_t *src,
        const uint32_t src_stride,
        int16_t *dst,
        const uint32_t dst_stride,
        int16_t *intermediate,
        uint32_t addshift);

    void EstimateInvTransform32x32_SSE2(
        int16_t                  *transformCoefficients,
        const uint32_t             src_stride,
        int16_t                  *residual,
        const uint32_t             dst_stride,
        int16_t                  *transform_inner_array_ptr,
        uint32_t                   bitIncrement);

    void EstimateInvTransform16x16_SSE2(
        int16_t                  *src,
        uint32_t                   src_stride,
        int16_t                  *dst,
        uint32_t                   dst_stride,
        int16_t                  *intermediate,
        uint32_t                   addshift);

#ifdef __cplusplus
}
#endif
#endif // EbTransforms_SSE2_h

