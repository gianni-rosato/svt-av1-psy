/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbCombinedAveragingSAD_AVX2_h
#define EbCombinedAveragingSAD_AVX2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    uint32_t CombinedAveraging8xMSAD_AVX2_INTRIN(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t CombinedAveraging16xMSAD_AVX2_INTRIN(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t CombinedAveraging24xMSAD_AVX2_INTRIN(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t CombinedAveraging32xMSAD_AVX2_INTRIN(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t CombinedAveraging48xMSAD_AVX2_INTRIN(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t CombinedAveraging64xMSAD_AVX2_INTRIN(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

    uint64_t ComputeMean8x8_AVX2_INTRIN(
        uint8_t *  input_samples,      // input parameter, input samples Ptr
        uint32_t   inputStride,       // input parameter, input stride
        uint32_t   inputAreaWidth,    // input parameter, input area width
        uint32_t   inputAreaHeight);

    void ComputeIntermVarFour8x8_AVX2_INTRIN(
        uint8_t *  input_samples,
        uint16_t   inputStride,
        uint64_t * meanOf8x8Blocks,      // mean of four  8x8
        uint64_t * meanOfSquared8x8Blocks);

    uint32_t combined_averaging_ssd_avx2(
        uint8_t   *src,
        ptrdiff_t  src_stride,
        uint8_t   *ref1,
        ptrdiff_t  ref1_stride,
        uint8_t   *ref2,
        ptrdiff_t  ref2_stride,
        uint32_t   height,
        uint32_t   width);
#ifdef __cplusplus
}
#endif
#endif