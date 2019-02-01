/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMeSadCalculation_asm_h
#define EbMeSadCalculation_asm_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif


    extern void InitializeBuffer_32bits_SSE2_INTRIN(
        uint32_t*        Pointer,
        uint32_t        Count128,
        uint32_t        Count32,
        uint32_t        Value);

    void SadCalculation_8x8_16x16_SSE2_INTRIN(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   refStride,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint32_t  *p_sad16x16);


    void SadCalculation_32x32_64x64_SSE2_INTRIN(
        uint32_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

#ifdef __cplusplus
}
#endif
#endif // EbMeSadCalculation_asm_h

