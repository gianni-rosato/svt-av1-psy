/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbComputeSAD_h
#define EbComputeSAD_h

#include "EbDefinitions.h"

#include "EbCombinedAveragingSAD_Intrinsic_AVX2.h"
#include "EbComputeSAD_C.h"
#include "EbComputeSAD_SSE2.h"
#include "EbComputeSAD_SSE4_1.h"
#include "EbComputeSAD_AVX2.h"
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
#include "EbUtility.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif


    /***************************************
    * Function Ptr Types
    ***************************************/
    typedef uint32_t(*EB_SADKERNELNxM_TYPE)(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref,
        uint32_t  refStride,
        uint32_t  height,
        uint32_t  width);

    static void NxMSadKernelVoidFunc() {}

    typedef void(*EB_SADLOOPKERNELNxM_TYPE)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  refStride,                      // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width,                          // input parameter, block width (N)
        uint64_t *bestSad,
        int16_t *xSearchCenter,
        int16_t *ySearchCenter,
        uint32_t  srcStrideRaw,                   // input parameter, source stride (no line skipping)
        int16_t search_area_width,
        int16_t search_area_height);

    typedef uint32_t(*EB_SADAVGKERNELNxM_TYPE)(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1Stride,
        uint8_t  *ref2,
        uint32_t  ref2Stride,
        uint32_t  height,
        uint32_t  width);

    typedef uint32_t(*EB_COMPUTE8X4SAD_TYPE)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  refStride);                     // input parameter, reference stride

    typedef uint32_t(*EB_COMPUTE8X8SAD_TYPE)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  refStride);                     // input parameter, reference stride

    typedef void(*EB_GETEIGHTSAD8x8)(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   refStride,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint16_t  *p_sad16x16);

    typedef void(*EB_GETEIGHTSAD32x32)(
        uint16_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

    typedef uint32_t(*CombinedAveragingSsd)(
        uint8_t   *src,
        ptrdiff_t  src_stride,
        uint8_t   *ref1,
        ptrdiff_t  ref1_stride,
        uint8_t   *ref2,
        ptrdiff_t  ref2_stride,
        uint32_t   height,
        uint32_t   width
        );

    /***************************************
    * Function Tables
    ***************************************/
    static EB_SADKERNELNxM_TYPE FUNC_TABLE NxMSadKernelSubSampled_funcPtrArray[ASM_TYPE_TOTAL][17] =   // [ASMTYPE][SAD - block height]
    {
        // NON_AVX2
        {
            /*0 4xM  */ Compute4xMSadSub_SSE2_INTRIN,
            /*1 8xM  */ FastLoop_NxMSadKernel,
            /*2 16xM */ FastLoop_NxMSadKernel,
            /*3 24xM */ FastLoop_NxMSadKernel,
            /*4 32xM */ FastLoop_NxMSadKernel,
            /*5      */ 0,
            /*6 48xM */ FastLoop_NxMSadKernel,
            /*7      */ 0,
            /*8 64xM */ FastLoop_NxMSadKernel,
            0,0,0,0,0,0,0,FastLoop_NxMSadKernel
        },
        // AVX2
        {
            /*0 4xM  */ Compute4xMSadSub_SSE2_INTRIN,
            /*1 8xM  */ Compute8xMSad_AVX2_INTRIN,
            /*2 16xM */ Compute16xMSad_AVX2_INTRIN,
            /*3 24xM */ FastLoop_NxMSadKernel,
            /*4 32xM */ Compute32xMSad_AVX2_INTRIN,
            /*5      */ 0,
            /*6 48xM */ FastLoop_NxMSadKernel,
            /*7      */ 0,
            /*8 64xM */ Compute64xMSad_AVX2_INTRIN,
            0,0,0,0,0,0,0,FastLoop_NxMSadKernel
        },
    };
    static EB_SADKERNELNxM_TYPE FUNC_TABLE NxMSadKernel_funcPtrArray[ASM_TYPE_TOTAL][9] =   // [ASMTYPE][SAD - block height]
    {
        // NON_AVX2
        {
            /*0 4xM  */ FastLoop_NxMSadKernel,
            /*1 8xM  */ FastLoop_NxMSadKernel,
            /*2 16xM */ FastLoop_NxMSadKernel,
            /*3 24xM */ FastLoop_NxMSadKernel,
            /*4 32xM */ FastLoop_NxMSadKernel,
            /*5      */ FastLoop_NxMSadKernel,  // size not supported in asm
            /*6 48xM */ FastLoop_NxMSadKernel,
            /*7      */ FastLoop_NxMSadKernel,  // size not supported in asm
            /*8 64xM */ FastLoop_NxMSadKernel
        },
        // AVX2
        {
            /*0 4xM  */ Compute4xMSad_AVX2_INTRIN,
            /*1 8xM  */ Compute8xMSad_AVX2_INTRIN,
            /*2 16xM */ Compute16xMSad_AVX2_INTRIN,//Compute16xMSad_AVX2_INTRIN is slower than the SSE2 version
            /*3 24xM */ Compute24xMSad_AVX2_INTRIN,
            /*4 32xM */ Compute32xMSad_AVX2_INTRIN,
            /*5      */ (EB_SADKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*6 48xM */ Compute48xMSad_AVX2_INTRIN,
            /*7      */ (EB_SADKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*8 64xM */ Compute64xMSad_AVX2_INTRIN,
        },
    };

    static EB_SADAVGKERNELNxM_TYPE FUNC_TABLE NxMSadAveragingKernel_funcPtrArray[ASM_TYPE_TOTAL][9] =   // [ASMTYPE][SAD - block height]
    {
        // NON_AVX2
        {
            /*0 4xM  */     CombinedAveragingSAD,
            /*1 8xM  */     CombinedAveragingSAD,
            /*2 16xM */     CombinedAveragingSAD,
            /*3 24xM */     CombinedAveragingSAD,
            /*4 32xM */     CombinedAveragingSAD,
            /*5      */     (EB_SADAVGKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*6 48xM */     CombinedAveragingSAD,
            /*7      */     (EB_SADAVGKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*8 64xM */     CombinedAveragingSAD
        },
        // AVX2
        {
            /*0 4xM  */     CombinedAveraging4xMSAD_SSE2_INTRIN,
            /*1 8xM  */     CombinedAveraging8xMSAD_AVX2_INTRIN,
            /*2 16xM */     CombinedAveraging16xMSAD_AVX2_INTRIN,
            /*3 24xM */     CombinedAveraging24xMSAD_AVX2_INTRIN,
            /*4 32xM */     CombinedAveraging32xMSAD_AVX2_INTRIN,
            /*5      */     (EB_SADAVGKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*6 48xM */     CombinedAveraging48xMSAD_AVX2_INTRIN,
            /*7      */     (EB_SADAVGKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*8 64xM */     CombinedAveraging64xMSAD_AVX2_INTRIN
        },
    };

    static EB_SADLOOPKERNELNxM_TYPE FUNC_TABLE NxMSadLoopKernelSparse_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        SadLoopKernelSparse_SSE4_1_INTRIN,
        // AVX2
        SadLoopKernelSparse_AVX2_INTRIN,
    };


    static EB_SADLOOPKERNELNxM_TYPE FUNC_TABLE NxMSadLoopKernel_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        SadLoopKernel_SSE4_1_INTRIN,
        // AVX2
        SadLoopKernel_AVX2_INTRIN,
    };

    static EB_GETEIGHTSAD8x8 FUNC_TABLE GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        GetEightHorizontalSearchPointResults_8x8_16x16_PU_SSE41_INTRIN,
        // AVX2
        GetEightHorizontalSearchPointResults_8x8_16x16_PU_AVX2_INTRIN,
    };

    static EB_GETEIGHTSAD32x32 FUNC_TABLE GetEightHorizontalSearchPointResults_32x32_64x64_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        GetEightHorizontalSearchPointResults_32x32_64x64_PU_SSE41_INTRIN,
        // AVX2
        GetEightHorizontalSearchPointResults_32x32_64x64_PU_AVX2_INTRIN,
    };

    uint32_t combined_averaging_ssd_c(
        uint8_t   *src,
        ptrdiff_t  src_stride,
        uint8_t   *ref1,
        ptrdiff_t  ref1_stride,
        uint8_t   *ref2,
        ptrdiff_t  ref2_stride,
        uint32_t   height,
        uint32_t   width);

    static CombinedAveragingSsd FUNC_TABLE combined_averaging_ssd_func_ptr_array[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        combined_averaging_ssd_c,
        // AVX2
        combined_averaging_ssd_avx2,
    };

#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_h
