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
        const uint8_t  *src,
        uint32_t  src_stride,
        const uint8_t  *ref,
        uint32_t  ref_stride,
        uint32_t  height,
        uint32_t  width);

    static void NxMSadKernelVoidFunc() {}

    typedef void(*EB_SADLOOPKERNELNxM_TYPE)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride,                      // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width,                          // input parameter, block width (N)
        uint64_t *best_sad,
        int16_t *x_search_center,
        int16_t *y_search_center,
        uint32_t  src_stride_raw,                   // input parameter, source stride (no line skipping)
        int16_t search_area_width,
        int16_t search_area_height);

    typedef uint32_t(*EB_SADAVGKERNELNxM_TYPE)(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

    typedef uint32_t(*EB_COMPUTE8X4SAD_TYPE)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride);                     // input parameter, reference stride

    typedef uint32_t(*EB_COMPUTE8X8SAD_TYPE)(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                      // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride);                     // input parameter, reference stride

    typedef void(*EB_GETEIGHTSAD8x8)(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
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
            /*0 4xM  */ fast_loop_nx_m_sad_kernel,
            /*1 8xM  */ fast_loop_nx_m_sad_kernel,
            /*2 16xM */ fast_loop_nx_m_sad_kernel,
            /*3 24xM */ fast_loop_nx_m_sad_kernel,
            /*4 32xM */ fast_loop_nx_m_sad_kernel,
            /*5      */ 0,
            /*6 48xM */ fast_loop_nx_m_sad_kernel,
            /*7      */ 0,
            /*8 64xM */ fast_loop_nx_m_sad_kernel,
            0,0,0,0,0,0,0,fast_loop_nx_m_sad_kernel
        },
        // AVX2
        {
            /*0 4xM  */ Compute4xMSadSub_AVX2_INTRIN,
            /*1 8xM  */ compute8x_m_sad_avx2_intrin,
            /*2 16xM */ compute16x_m_sad_avx2_intrin,
            /*3 24xM */ fast_loop_nx_m_sad_kernel,
            /*4 32xM */ compute32x_m_sad_avx2_intrin,
            /*5      */ 0,
            /*6 48xM */ fast_loop_nx_m_sad_kernel,
            /*7      */ 0,
            /*8 64xM */ compute64x_m_sad_avx2_intrin,
            0,0,0,0,0,0,0,fast_loop_nx_m_sad_kernel
        },
    };
    static EB_SADKERNELNxM_TYPE FUNC_TABLE NxMSadKernel_funcPtrArray[ASM_TYPE_TOTAL][9] =   // [ASMTYPE][SAD - block height]
    {
        // NON_AVX2
        {
            /*0 4xM  */ fast_loop_nx_m_sad_kernel,
            /*1 8xM  */ fast_loop_nx_m_sad_kernel,
            /*2 16xM */ fast_loop_nx_m_sad_kernel,
            /*3 24xM */ fast_loop_nx_m_sad_kernel,
            /*4 32xM */ fast_loop_nx_m_sad_kernel,
            /*5      */ fast_loop_nx_m_sad_kernel,  // size not supported in asm
            /*6 48xM */ fast_loop_nx_m_sad_kernel,
            /*7      */ fast_loop_nx_m_sad_kernel,  // size not supported in asm
            /*8 64xM */ fast_loop_nx_m_sad_kernel
        },
        // AVX2
        {
            /*0 4xM  */ compute4x_m_sad_avx2_intrin,
            /*1 8xM  */ compute8x_m_sad_avx2_intrin,
            /*2 16xM */ compute16x_m_sad_avx2_intrin,//compute16x_m_sad_avx2_intrin is slower than the SSE2 version
            /*3 24xM */ compute24x_m_sad_avx2_intrin,
            /*4 32xM */ compute32x_m_sad_avx2_intrin,
            /*5      */ (EB_SADKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*6 48xM */ compute48x_m_sad_avx2_intrin,
            /*7      */ (EB_SADKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*8 64xM */ compute64x_m_sad_avx2_intrin,
        },
    };

    static EB_SADAVGKERNELNxM_TYPE FUNC_TABLE NxMSadAveragingKernel_funcPtrArray[ASM_TYPE_TOTAL][9] =   // [ASMTYPE][SAD - block height]
    {
        // NON_AVX2
        {
            /*0 4xM  */     combined_averaging_sad,
            /*1 8xM  */     combined_averaging_sad,
            /*2 16xM */     combined_averaging_sad,
            /*3 24xM */     combined_averaging_sad,
            /*4 32xM */     combined_averaging_sad,
            /*5      */     (EB_SADAVGKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*6 48xM */     combined_averaging_sad,
            /*7      */     (EB_SADAVGKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*8 64xM */     combined_averaging_sad
        },
        // AVX2
        {
            /*0 4xM  */     CombinedAveraging4xMSAD_SSE2_INTRIN,
            /*1 8xM  */     combined_averaging8x_msad_avx2_intrin,
            /*2 16xM */     combined_averaging16x_msad_avx2_intrin,
            /*3 24xM */     combined_averaging24x_msad_avx2_intrin,
            /*4 32xM */     combined_averaging32x_msad_avx2_intrin,
            /*5      */     (EB_SADAVGKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*6 48xM */     combined_averaging48x_msad_avx2_intrin,
            /*7      */     (EB_SADAVGKERNELNxM_TYPE)NxMSadKernelVoidFunc,
            /*8 64xM */     combined_averaging64x_msad_avx2_intrin
        },
    };

    static EB_SADLOOPKERNELNxM_TYPE FUNC_TABLE NxMSadLoopKernelSparse_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        sad_loop_kernel_sparse_sse4_1_intrin,
        // AVX2
        sad_loop_kernel_sparse_avx2_intrin,
    };


    static EB_SADLOOPKERNELNxM_TYPE FUNC_TABLE NxMSadLoopKernel_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        sad_loop_kernel_sse4_1_intrin,
        // AVX2
        sad_loop_kernel_avx2_intrin,
    };

    static EB_GETEIGHTSAD8x8 FUNC_TABLE GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin,
        // AVX2
        get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin,
    };

    static EB_GETEIGHTSAD32x32 FUNC_TABLE GetEightHorizontalSearchPointResults_32x32_64x64_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        get_eight_horizontal_search_point_results_32x32_64x64_pu_sse41_intrin,
        // AVX2
        get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin,
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
