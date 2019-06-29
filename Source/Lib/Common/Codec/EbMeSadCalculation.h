/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMeSadCalculation_h
#define EbMeSadCalculation_h

#include "EbMeSadCalculation_SSE2.h"
#include "EbComputeSAD_SSE4_1.h"

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    /***************************************
    * Function Types
    ***************************************/
    typedef void(*EbSadCalculation8x8and16x16Type)(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint32_t  *p_sad16x16,
        EbBool     sub_sad);

    typedef void(*EbExtSadCalculation8x8and16x16Type)(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint32_t  *p_sad16x16,
        uint32_t  *p_sad8x8,
        EbBool    sub_sad);

    typedef void(*EbExtSadCalculation32x32and64x64Type)(
        uint32_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv,
        uint32_t  *p_sad32x32);

    typedef void(*EbExtSadCalculationType)(
        uint32_t  *p_sad8x8,
        uint32_t  *p_sad16x16,
        uint32_t  *p_sad32x32,
        uint32_t  *p_best_sad64x32,
        uint32_t  *p_best_mv64x32,
        uint32_t  *p_best_sad32x16,
        uint32_t  *p_best_mv32x16,
        uint32_t  *p_best_sad16x8,
        uint32_t  *p_best_mv16x8,
        uint32_t  *p_best_sad32x64,
        uint32_t  *p_best_mv32x64,
        uint32_t  *p_best_sad16x32,
        uint32_t  *p_best_mv16x32,
        uint32_t  *p_best_sad8x16,
        uint32_t  *p_best_mv8x16,
        uint32_t  *p_best_sad32x8,
        uint32_t  *p_best_mv32x8,
        uint32_t  *p_best_sad8x32,
        uint32_t  *p_best_mv8x32,
        uint32_t  *p_best_sad64x16,
        uint32_t  *p_best_mv64x16,
        uint32_t  *p_best_sad16x64,
        uint32_t  *p_best_mv16x64,
        uint32_t   mv);

    typedef void(*EbSadCalculation32x32and64x64Type)(
        uint32_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

    typedef void(*EbInializeBuffer32Bits)(
        uint32_t*        pointer,
        uint32_t        count128,
        uint32_t        count32,
        uint32_t        value);

    typedef void(*EB_EXT_ALL_SAD_CALCULATION_8x8_16x16_TYPE)(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t   mv,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_mv16x16,
        uint32_t   p_eight_sad16x16[16][8],
        uint32_t   p_eight_sad8x8[64][8]);

    typedef void(*EbEightSadCalculationNsqType)(
        uint32_t   p_sad8x8[64][8],
        uint32_t   p_sad16x16[16][8],
        uint32_t   p_sad32x32[4][8],
        uint32_t  *p_best_sad64x32,
        uint32_t  *p_best_mv64x32,
        uint32_t  *p_best_sad32x16,
        uint32_t  *p_best_mv32x16,
        uint32_t  *p_best_sad16x8,
        uint32_t  *p_best_mv16x8,
        uint32_t  *p_best_sad32x64,
        uint32_t  *p_best_mv32x64,
        uint32_t  *p_best_sad16x32,
        uint32_t  *p_best_mv16x32,
        uint32_t  *p_best_sad8x16,
        uint32_t  *p_best_mv8x16,
        uint32_t  *p_best_sad32x8,
        uint32_t  *p_best_mv32x8,
        uint32_t  *p_best_sad8x32,
        uint32_t  *p_best_mv8x32,
        uint32_t  *p_best_sad64x16,
        uint32_t  *p_best_mv64x16,
        uint32_t  *p_best_sad16x64,
        uint32_t  *p_best_mv16x64,
        uint32_t   mv);

    typedef void(*EbExtEightSadCalculation32x3264x64Type)(
        uint32_t  p_sad16x16[16][8],
        uint32_t *p_best_sad32x32,
        uint32_t *p_best_sad64x64,
        uint32_t *p_best_mv32x32,
        uint32_t *p_best_mv64x64,
        uint32_t  mv,
        uint32_t p_sad32x32[4][8]);
    static EbInializeBuffer32Bits initialize_buffer32bits_func_ptr_array[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        initialize_buffer_32bits_sse2_intrin,
        // AVX2
        initialize_buffer_32bits_sse2_intrin
    };

#ifdef __cplusplus
}
#endif
#endif // EbMeSadCalculation_h
