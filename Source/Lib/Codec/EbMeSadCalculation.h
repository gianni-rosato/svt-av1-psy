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
    typedef void(*EB_SADCALCULATION8X8AND16X16_TYPE)(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint32_t  *p_sad16x16);

    typedef void(*EB_EXTSADCALCULATION8X8AND16X16_TYPE)(
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
        uint32_t  *p_sad8x8);

    typedef void(*EB_EXTSADCALCULATION32X32AND64X64_TYPE)(
        uint32_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv,
        uint32_t  *p_sad32x32);

    typedef void(*EB_EXTSADCALCULATION_TYPE)(
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

    typedef void(*EB_SADCALCULATION32X32AND64X64_TYPE)(
        uint32_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

    typedef void(*EB_INIALIZEBUFFER_32BITS)(
        uint32_t*        pointer,
        uint32_t        count128,
        uint32_t        count32,
        uint32_t        value);

    typedef void(*EB_RECTAMPSADCALCULATION8X8AND16X16_TYPE)(
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
        uint32_t  *p_sad8x8);

    typedef void(*EB_RECTAMPSADCALCULATION32X32AND64X64_TYPE)(
        uint32_t  *p_sad16x16,
        uint32_t  *p_sad32x32,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

    typedef void(*EB_RECTAMPSADCALCULATION_TYPE)(
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
        uint32_t  *p_best_sad64x16,
        uint32_t  *p_best_mv64x16,
        uint32_t  *p_best_sad32x8,
        uint32_t  *p_best_mv32x8,
        uint32_t  *pBestSad64x48,
        uint32_t  *pBestMV64x48,
        uint32_t  *pBestSad32x24,
        uint32_t  *pBestMV32x24,
        uint32_t  *p_best_sad16x64,
        uint32_t  *p_best_mv16x64,
        uint32_t  *p_best_sad8x32,
        uint32_t  *p_best_mv8x32,
        uint32_t  *pBestSad48x64,
        uint32_t  *pBestMV48x64,
        uint32_t  *pBestSad24x32,
        uint32_t  *pBestMV24x32,
        uint32_t   mv);

    static EB_INIALIZEBUFFER_32BITS InitializeBuffer_32bits_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        initialize_buffer_32bits_sse2_intrin,
        // AVX2
        initialize_buffer_32bits_sse2_intrin
    };

#ifdef __cplusplus
}
#endif
#endif // EbMeSadCalculation_h