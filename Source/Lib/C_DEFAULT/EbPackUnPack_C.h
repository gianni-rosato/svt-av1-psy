/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbPackUnPack_C_h
#define EbPackUnPack_C_h
#ifdef __cplusplus
extern "C" {
#endif

#include "EbDefinitions.h"


    void EB_ENC_msbPack2D(
        uint8_t     *in8BitBuffer,
        uint32_t     in8Stride,
        uint8_t     *innBitBuffer,
        uint16_t    *out16BitBuffer,
        uint32_t     innStride,
        uint32_t     outStride,
        uint32_t     width,
        uint32_t     height);

    void CompressedPackmsb(
        uint8_t     *in8BitBuffer,
        uint32_t     in8Stride,
        uint8_t     *innBitBuffer,
        uint16_t    *out16BitBuffer,
        uint32_t     innStride,
        uint32_t     outStride,
        uint32_t     width,
        uint32_t     height);
    void CPack_C(
        const uint8_t     *innBitBuffer,
        uint32_t     innStride,
        uint8_t     *inCompnBitBuffer,
        uint32_t     outStride,
        uint8_t    *localCache,
        uint32_t     width,
        uint32_t     height);



    void EB_ENC_msbUnPack2D(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint8_t       *outnBitBuffer,
        uint32_t       out8Stride,
        uint32_t       outnStride,
        uint32_t       width,
        uint32_t       height);
    void UnPack8BitData(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint32_t       out8Stride,
        uint32_t       width,
        uint32_t       height);
    void UnpackAvg(
        uint16_t      *ref16L0,
        uint32_t       refL0Stride,
        uint16_t      *ref16L1,
        uint32_t       refL1Stride,
        uint8_t       *dstPtr,
        uint32_t       dst_stride,
        uint32_t       width,
        uint32_t       height);

    void UnpackAvgSafeSub(
        uint16_t      *ref16L0,
        uint32_t       refL0Stride,
        uint16_t      *ref16L1,
        uint32_t       refL1Stride,
        uint8_t       *dstPtr,
        uint32_t       dst_stride,
        EbBool      subPred,
        uint32_t       width,
        uint32_t       height);


#ifdef __cplusplus
}
#endif
#endif // EbPackUnPack_C_h