/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_AVX2
#define EbPictureOperators_AVX2

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    extern void EB_ENC_msbPack2D_AVX2_INTRIN_AL(
        uint8_t     *in8BitBuffer,
        uint32_t     in8Stride,
        uint8_t     *innBitBuffer,
        uint16_t    *out16BitBuffer,
        uint32_t     innStride,
        uint32_t     outStride,
        uint32_t     width,
        uint32_t     height);

    extern void CompressedPackmsb_AVX2_INTRIN(
        uint8_t     *in8BitBuffer,
        uint32_t     in8Stride,
        uint8_t     *innBitBuffer,
        uint16_t    *out16BitBuffer,
        uint32_t     innStride,
        uint32_t     outStride,
        uint32_t     width,
        uint32_t     height);

    void CPack_AVX2_INTRIN(
        const uint8_t     *innBitBuffer,
        uint32_t     innStride,
        uint8_t     *inCompnBitBuffer,
        uint32_t     outStride,
        uint8_t    *localCache,
        uint32_t     width,
        uint32_t     height);


    void UnpackAvg_AVX2_INTRIN(
        uint16_t *ref16L0,
        uint32_t  refL0Stride,
        uint16_t *ref16L1,
        uint32_t  refL1Stride,
        uint8_t  *dstPtr,
        uint32_t  dst_stride,
        uint32_t  width,
        uint32_t  height);

    int32_t  sumResidual8bit_AVX2_INTRIN(
        int16_t * inPtr,
        uint32_t   size,
        uint32_t   strideIn);
    void memset16bitBlock_AVX2_INTRIN(
        int16_t * inPtr,
        uint32_t   strideIn,
        uint32_t   size,
        int16_t   value
    );


    void UnpackAvgSafeSub_AVX2_INTRIN(
        uint16_t *ref16L0,
        uint32_t  refL0Stride,
        uint16_t *ref16L1,
        uint32_t  refL1Stride,
        uint8_t  *dstPtr,
        uint32_t  dst_stride,
        EbBool  subPred,
        uint32_t  width,
        uint32_t  height);

    void PictureAdditionKernel4x4_AV1_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int32_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height,
        int32_t     bd);
    void PictureAdditionKernel8x8_AV1_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int32_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height,
        int32_t     bd);
    void PictureAdditionKernel16x16_AV1_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int32_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height,
        int32_t     bd);
    void PictureAdditionKernel32x32_AV1_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int32_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height,
        int32_t     bd);
    void PictureAdditionKernel64x64_AV1_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int32_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height,
        int32_t     bd);

    void FullDistortionKernelCbfZero32Bits_AVX2(
        int32_t  *coeff,
        uint32_t   coeffStride,
        int32_t  *reconCoeff,
        uint32_t   reconCoeffStride,
        uint64_t   distortionResult[DIST_CALC_TOTAL],
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    void FullDistortionKernel32Bits_AVX2(
        int32_t  *coeff,
        uint32_t   coeffStride,
        int32_t  *reconCoeff,
        uint32_t   reconCoeffStride,
        uint64_t   distortionResult[DIST_CALC_TOTAL],
        uint32_t   areaWidth,
        uint32_t   areaHeight);


#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_AVX2