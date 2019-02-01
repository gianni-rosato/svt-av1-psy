/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPackUnPack_asm_h
#define EbPackUnPack_asm_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    void EB_ENC_msbPack2D_SSE2_INTRIN(
        uint8_t     *in8BitBuffer,
        uint32_t     in8Stride,
        uint8_t     *innBitBuffer,
        uint16_t    *out16BitBuffer,
        uint32_t     innStride,
        uint32_t     outStride,
        uint32_t     width,
        uint32_t     height);

    void EB_ENC_msbUnPack2D_SSE2_INTRIN(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint8_t       *outnBitBuffer,
        uint32_t       out8Stride,
        uint32_t       outnStride,
        uint32_t       width,
        uint32_t       height);
    void EB_ENC_UnPack8BitData_SSE2_INTRIN(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint32_t       out8Stride,
        uint32_t       width,
        uint32_t       height);

    void EB_ENC_UnPack8BitDataSafeSub_SSE2_INTRIN(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint32_t       out8Stride,
        uint32_t       width,
        uint32_t       height,
        EbBool      subPred
    );

    void UnpackAvg_SSE2_INTRIN(
        uint16_t *ref16L0,
        uint32_t  refL0Stride,
        uint16_t *ref16L1,
        uint32_t  refL1Stride,
        uint8_t  *dstPtr,
        uint32_t  dst_stride,
        uint32_t  width,
        uint32_t  height);


#ifdef __cplusplus
}
#endif
#endif // EbPackUnPack_asm_h

