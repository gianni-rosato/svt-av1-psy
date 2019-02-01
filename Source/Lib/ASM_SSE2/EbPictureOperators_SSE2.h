/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_SSE2_h
#define EbPictureOperators_SSE2_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    //-----
    extern void ZeroOutCoeff4x4_SSE(
        int16_t*                  coeffbuffer,
        uint32_t                   coeffStride,
        uint32_t                   coeffOriginIndex,
        uint32_t                   areaWidth,
        uint32_t                   areaHeight);
    extern void ZeroOutCoeff8x8_SSE2(
        int16_t*                  coeffbuffer,
        uint32_t                   coeffStride,
        uint32_t                   coeffOriginIndex,
        uint32_t                   areaWidth,
        uint32_t                   areaHeight);
    extern void ZeroOutCoeff16x16_SSE2(
        int16_t*                  coeffbuffer,
        uint32_t                   coeffStride,
        uint32_t                   coeffOriginIndex,
        uint32_t                   areaWidth,
        uint32_t                   areaHeight);
    extern void ZeroOutCoeff32x32_SSE2(
        int16_t*                  coeffbuffer,
        uint32_t                   coeffStride,
        uint32_t                   coeffOriginIndex,
        uint32_t                   areaWidth,
        uint32_t                   areaHeight);

    extern void ResidualKernel16bit_SSE2_INTRIN(
        uint16_t   *input,
        uint32_t   inputStride,
        uint16_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    void PictureAdditionKernel4x4_SSE_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height);

    void PictureAdditionKernel8x8_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height);

    void PictureAdditionKernel16x16_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height);

    void PictureAdditionKernel32x32_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height);

    void PictureAdditionKernel64x64_SSE2_INTRIN(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height);

    void ResidualKernelSubSampled4x4_SSE_INTRIN(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight,
        uint8_t    lastLine);

    void ResidualKernelSubSampled8x8_SSE2_INTRIN(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight,
        uint8_t    lastLine);

    void ResidualKernelSubSampled16x16_SSE2_INTRIN(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight,
        uint8_t    lastLine);

    void ResidualKernelSubSampled32x32_SSE2_INTRIN(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight,
        uint8_t    lastLine);

    void ResidualKernelSubSampled64x64_SSE2_INTRIN(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight,
        uint8_t    lastLine);

    void PictureAdditionKernel16bit_SSE2_INTRIN(
        uint16_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint16_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height);



#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_SSE2_h
