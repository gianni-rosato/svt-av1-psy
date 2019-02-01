/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbTransforms_AVX2_h
#define EbTransforms_AVX2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    void QuantizeInvQuantize8x8_AVX2_INTRIN(
        int16_t          *coeff,
        const uint32_t     coeffStride,
        int16_t          *quantCoeff,
        int16_t          *reconCoeff,
        const uint32_t     qFunc,
        const uint32_t     q_offset,
        const int32_t     shiftedQBits,
        const int32_t     shiftedFFunc,
        const int32_t     iq_offset,
        const int32_t     shiftNum,
        const uint32_t     areaSize,
        uint32_t          *nonzerocoeff);

    void QuantizeInvQuantizeNxN_AVX2_INTRIN(
        int16_t          *coeff,
        const uint32_t     coeffStride,
        int16_t          *quantCoeff,
        int16_t          *reconCoeff,
        const uint32_t     qFunc,
        const uint32_t     q_offset,
        const int32_t     shiftedQBits,
        const int32_t     shiftedFFunc,
        const int32_t     iq_offset,
        const int32_t     shiftNum,
        const uint32_t     areaSize,
        uint32_t          *nonzerocoeff);

    void lowPrecisionTransform16x16_AVX2_INTRIN(int16_t *src, uint32_t src_stride, int16_t *dst, uint32_t dst_stride, int16_t *intermediate, uint32_t addshift);
    void lowPrecisionTransform32x32_AVX2_INTRIN(int16_t *src, uint32_t src_stride, int16_t *dst, uint32_t dst_stride, int16_t *intermediate, uint32_t addshift);

    void PfreqTransform32x32_AVX2_INTRIN(
        int16_t *src,
        const uint32_t src_stride,
        int16_t *dst,
        const uint32_t dst_stride,
        int16_t *intermediate,
        uint32_t addshift);

    void PfreqN4Transform32x32_AVX2_INTRIN(
        int16_t *src,
        const uint32_t src_stride,
        int16_t *dst,
        const uint32_t dst_stride,
        int16_t *intermediate,
        uint32_t addshift);

    void MatMult4x4_AVX2_INTRIN(
        int16_t*              coeff,
        const uint32_t         coeffStride,
        const uint16_t        *maskingMatrix,
        const uint32_t         maskingMatrixStride,  //Matrix size
        const uint32_t         computeSize,  //Computation area size
        const int32_t         offset,     //(PMP_MAX >> 1)
        const int32_t         shiftNum, //PMP_PRECISION
        uint32_t*              nonzerocoeff);

    void MatMult4x4_OutBuff_AVX2_INTRIN(
        int16_t*              coeff,
        const uint32_t         coeffStride,
        int16_t*              coeffOut,
        const uint32_t         coeffOutStride,

        const uint16_t        *maskingMatrix,
        const uint32_t         maskingMatrixStride,
        const uint32_t         computeSize,
        const int32_t         offset,
        const int32_t         shiftNum,
        uint32_t*              nonzerocoeff);


    void MatMult8x8_AVX2_INTRIN(
        int16_t*              coeff,
        const uint32_t         coeffStride,
        const uint16_t        *maskingMatrix,
        const uint32_t         maskingMatrixStride,  //Matrix size
        const uint32_t         computeSize,  //Computation area size
        const int32_t         offset,     //(PMP_MAX >> 1)
        const int32_t         shiftNum, //PMP_PRECISION
        uint32_t*              nonzerocoeff);

    void MatMultNxN_AVX2_INTRIN(
        int16_t*              coeff,
        const uint32_t         coeffStride,
        const uint16_t        *maskingMatrix,
        const uint32_t         maskingMatrixStride,  //Matrix size
        const uint32_t         computeSize,  //Computation area size
        const int32_t         offset,     //(PMP_MAX >> 1)
        const int32_t         shiftNum, //PMP_PRECISION
        uint32_t*              nonzerocoeff);



#ifdef __cplusplus
}
#endif
#endif // EbTransforms_AVX2_h


