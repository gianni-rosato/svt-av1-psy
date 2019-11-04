/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#define RTCD_C
#include "aom_dsp_rtcd.h"
#include "EbComputeSAD_C.h"
#include "EbComputeSAD_SSE4_1.h"
#include "EbComputeSAD_AVX2.h"
#include "EbPictureAnalysisProcess.h"
#include "EbTemporalFiltering.h"
#include "EbTemporalFiltering_sse4.h"
#include "EbComputeSAD.h"
#include "EbMotionEstimation.h"
#include "EbPictureOperators.h"
#include "EbPackUnPack_C.h"
#include "EbPackUnPack_SSE2.h"
#include "EbPackUnPack_AVX2.h"
#include "EbMcp_SSE2.h"
#include "EbAvcStyleMcp_SSSE3.h"
#include "EbComputeMean_SSE2.h"
#include "EbCombinedAveragingSAD_Intrinsic_AVX2.h"
#include "EbComputeMean.h"
#include "EbHmCode.h"
#include "EbMeSadCalculation.h"
#include "EbAvcStyleMcp.h"

/*
 * DSP deprecated flags
 */
#define HAS_MMX CPU_FLAGS_MMX
#define HAS_SSE CPU_FLAGS_SSE
#define HAS_SSE2 CPU_FLAGS_SSE2
#define HAS_SSE3 CPU_FLAGS_SSE3
#define HAS_SSSE3 CPU_FLAGS_SSSE3
#define HAS_SSE4_1 CPU_FLAGS_SSE4_1
#define HAS_SSE4_2 CPU_FLAGS_SSE4_2
#define HAS_AVX CPU_FLAGS_AVX
#define HAS_AVX2 CPU_FLAGS_AVX2
#define HAS_AVX512F CPU_FLAGS_AVX512F
#define HAS_AVX512CD CPU_FLAGS_AVX512CD
#define HAS_AVX512DQ CPU_FLAGS_AVX512DQ
#define HAS_AVX512ER CPU_FLAGS_AVX512ER
#define HAS_AVX512PF CPU_FLAGS_AVX512PF
#define HAS_AVX512BW CPU_FLAGS_AVX512BW
#define HAS_AVX512VL CPU_FLAGS_AVX512VL


 /**************************************
 * Instruction Set Support
 **************************************/

 // Helper Functions
static INLINE void RunCpuid(uint32_t eax, uint32_t ecx, int32_t* abcd)
{
#ifdef _WIN32
    __cpuidex(abcd, eax, ecx);
#else
    uint32_t ebx = 0, edx = 0;
# if defined( __i386__ ) && defined ( __PIC__ )
    /* in case of PIC under 32-bit EBX cannot be clobbered */
    __asm__("movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
# else
    __asm__("cpuid" : "+b" (ebx),
# endif
        "+a" (eax), "+c" (ecx), "=d" (edx));
    abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;
#endif
}

static INLINE int32_t CheckXcr0Ymm()
{
    uint32_t xcr0;
#ifdef _WIN32
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}

static int32_t Check4thGenIntelCoreFeatures()
{
    int32_t abcd[4];
    int32_t fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
    int32_t avx2_bmi12_mask = (1 << 5) | (1 << 3) | (1 << 8);

    /* CPUID.(EAX=01H, ECX=0H):ECX.FMA[bit 12]==1   &&
    CPUID.(EAX=01H, ECX=0H):ECX.MOVBE[bit 22]==1 &&
    CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    RunCpuid(1, 0, abcd);
    if ((abcd[2] & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask)
        return 0;

    if (!CheckXcr0Ymm())
        return 0;

    /*  CPUID.(EAX=07H, ECX=0H):EBX.AVX2[bit 5]==1  &&
    CPUID.(EAX=07H, ECX=0H):EBX.BMI1[bit 3]==1  &&
    CPUID.(EAX=07H, ECX=0H):EBX.BMI2[bit 8]==1  */
    RunCpuid(7, 0, abcd);
    if ((abcd[1] & avx2_bmi12_mask) != avx2_bmi12_mask)
        return 0;
    /* CPUID.(EAX=80000001H):ECX.LZCNT[bit 5]==1 */
    RunCpuid(0x80000001, 0, abcd);
    if ((abcd[2] & (1 << 5)) == 0)
        return 0;
    return 1;
}

static INLINE int CheckXcr0Zmm()
{
    uint32_t xcr0;
    uint32_t zmm_ymm_xmm = (7 << 5) | (1 << 2) | (1 << 1);
#ifdef _WIN32
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
#endif
    return ((xcr0 & zmm_ymm_xmm) == zmm_ymm_xmm); /* check if xmm, ymm and zmm state are enabled in XCR0 */
}

static int32_t CanUseIntelAVX512()
{
    int abcd[4];

    /*  CPUID.(EAX=07H, ECX=0):EBX[bit 16]==1 AVX512F
    CPUID.(EAX=07H, ECX=0):EBX[bit 17] AVX512DQ
    CPUID.(EAX=07H, ECX=0):EBX[bit 28] AVX512CD
    CPUID.(EAX=07H, ECX=0):EBX[bit 30] AVX512BW
    CPUID.(EAX=07H, ECX=0):EBX[bit 31] AVX512VL */

    int avx512_ebx_mask =
        (1 << 16)  // AVX-512F
        | (1 << 17)  // AVX-512DQ
        | (1 << 28)  // AVX-512CD
        | (1 << 30)  // AVX-512BW
        | (1 << 31); // AVX-512VL

    if (!Check4thGenIntelCoreFeatures())
        return 0;

    // ensure OS supports ZMM registers (and YMM, and XMM)
    if (!CheckXcr0Zmm())
        return 0;

    RunCpuid(7, 0, abcd);
    if ((abcd[1] & avx512_ebx_mask) != avx512_ebx_mask)
        return 0;

    return 1;
}

CPU_FLAGS get_cpu_flags() {
    CPU_FLAGS flags = 0;

    /* To detail tests CPU features, requires more accurate implementation.
        Documentation help:
        https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex?redirectedfrom=MSDN&view=vs-2019
    */

    if (Check4thGenIntelCoreFeatures()) {
        flags |= CPU_FLAGS_MMX | CPU_FLAGS_SSE | CPU_FLAGS_SSE2
              | CPU_FLAGS_SSE3 | CPU_FLAGS_SSSE3 | CPU_FLAGS_SSE4_1
              | CPU_FLAGS_SSE4_2 | CPU_FLAGS_AVX | CPU_FLAGS_AVX2;
    }

    if (CanUseIntelAVX512()) {
        flags |= CPU_FLAGS_AVX512F | CPU_FLAGS_AVX512DQ | CPU_FLAGS_AVX512CD
              | CPU_FLAGS_AVX512BW | CPU_FLAGS_AVX512VL;
    }

    return flags;
}

CPU_FLAGS get_cpu_flags_to_use() {
    CPU_FLAGS flags = get_cpu_flags();
#ifdef NON_AVX512_SUPPORT
    /* Remove AVX512 flags. */
    flags &= (CPU_FLAGS_AVX512F - 1);
#endif
    return flags;
}

#ifndef NON_AVX512_SUPPORT
#define SET_FUNCTIONS(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512) \
do { ptr = c;\
    if (((uintptr_t)NULL != (uintptr_t)mmx) && (flags & HAS_MMX)) ptr = mmx; \
    if (((uintptr_t)NULL != (uintptr_t)sse) && (flags & HAS_SSE)) ptr = sse; \
    if (((uintptr_t)NULL != (uintptr_t)sse2) && (flags & HAS_SSE2)) ptr = sse2; \
    if (((uintptr_t)NULL != (uintptr_t)sse3) && (flags & HAS_SSE3)) ptr = sse3; \
    if (((uintptr_t)NULL != (uintptr_t)ssse3) && (flags & HAS_SSSE3)) ptr = ssse3; \
    if (((uintptr_t)NULL != (uintptr_t)sse4_1) && (flags & HAS_SSE4_1)) ptr = sse4_1; \
    if (((uintptr_t)NULL != (uintptr_t)sse4_2) && (flags & HAS_SSE4_2)) ptr = sse4_2; \
    if (((uintptr_t)NULL != (uintptr_t)avx) && (flags & HAS_AVX)) ptr = avx; \
    if (((uintptr_t)NULL != (uintptr_t)avx2) && (flags & HAS_AVX2)) ptr = avx2; \
    if (((uintptr_t)NULL != (uintptr_t)avx512) && (flags & HAS_AVX512F)) ptr = avx512; \
} while(0)
#else
#define SET_FUNCTIONS(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512) \
do { ptr = c;\
    if (((uintptr_t)NULL != (uintptr_t)mmx) && (flags & HAS_MMX)) ptr = mmx; \
    if (((uintptr_t)NULL != (uintptr_t)sse) && (flags & HAS_SSE)) ptr = sse; \
    if (((uintptr_t)NULL != (uintptr_t)sse2) && (flags & HAS_SSE2)) ptr = sse2; \
    if (((uintptr_t)NULL != (uintptr_t)sse3) && (flags & HAS_SSE3)) ptr = sse3; \
    if (((uintptr_t)NULL != (uintptr_t)ssse3) && (flags & HAS_SSSE3)) ptr = ssse3; \
    if (((uintptr_t)NULL != (uintptr_t)sse4_1) && (flags & HAS_SSE4_1)) ptr = sse4_1; \
    if (((uintptr_t)NULL != (uintptr_t)sse4_2) && (flags & HAS_SSE4_2)) ptr = sse4_2; \
    if (((uintptr_t)NULL != (uintptr_t)avx) && (flags & HAS_AVX)) ptr = avx; \
    if (((uintptr_t)NULL != (uintptr_t)avx2) && (flags & HAS_AVX2)) ptr = avx2; \
} while(0)
#endif

#define SET_SSE2(ptr, c, sse2) SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, 0, 0)
#define SET_SSE2_AVX2(ptr, c, sse2, avx2) SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, avx2, 0)
#define SET_SSSE3(ptr, c, ssse3) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, ssse3, 0, 0, 0, 0, 0)
#define SET_SSE41(ptr, c, sse4_1) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, 0, 0)
#define SET_SSE41(ptr, c, sse4_1) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, 0, 0)
#define SET_SSE41_AVX2(ptr, c, sse4_1, avx2) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, avx2, 0)
#define SET_SSE41_AVX2_AVX512(ptr, c, sse4_1, avx2, avx512) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, avx2, avx512)
#define SET_AVX2(ptr, c, avx2) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, avx2, 0)
#define SET_AVX2_AVX512(ptr, c, avx2, avx512) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, avx2, avx512)


void setup_rtcd_internal(CPU_FLAGS flags)
{
    /** Should be done during library initialization,
        but for safe limiting cpu flags again. */
    flags &= get_cpu_flags_to_use();

    //to use C: flags=0

    eb_apply_selfguided_restoration = eb_apply_selfguided_restoration_c;
    if (flags & HAS_AVX2) eb_apply_selfguided_restoration = eb_apply_selfguided_restoration_avx2;

    eb_av1_wiener_convolve_add_src = eb_av1_wiener_convolve_add_src_c;
    if (flags & HAS_AVX2) eb_av1_wiener_convolve_add_src = eb_av1_wiener_convolve_add_src_avx2;

    eb_av1_highbd_wiener_convolve_add_src = eb_av1_highbd_wiener_convolve_add_src_c;
    if (flags & HAS_AVX2) eb_av1_highbd_wiener_convolve_add_src = eb_av1_highbd_wiener_convolve_add_src_avx2;

    eb_av1_selfguided_restoration = eb_av1_selfguided_restoration_c;
    if (flags & HAS_AVX2) eb_av1_selfguided_restoration = eb_av1_selfguided_restoration_avx2;
    av1_build_compound_diffwtd_mask = av1_build_compound_diffwtd_mask_c;
    if (flags & HAS_AVX2) av1_build_compound_diffwtd_mask = av1_build_compound_diffwtd_mask_avx2;
#if COMP_HBD
     av1_build_compound_diffwtd_mask_highbd = av1_build_compound_diffwtd_mask_highbd_c;
     if (flags & HAS_AVX2) av1_build_compound_diffwtd_mask_highbd = av1_build_compound_diffwtd_mask_highbd_avx2;
#endif
    av1_wedge_sse_from_residuals = av1_wedge_sse_from_residuals_c;
    if (flags & HAS_AVX2) av1_wedge_sse_from_residuals = av1_wedge_sse_from_residuals_avx2;
    aom_subtract_block = aom_subtract_block_c;
    if (flags & HAS_AVX2) aom_subtract_block = aom_subtract_block_avx2;
    aom_sse = aom_sse_c;
#if COMP_HBD
    aom_highbd_subtract_block = aom_highbd_subtract_block_c;
    if (flags & HAS_AVX2) aom_highbd_subtract_block = aom_highbd_subtract_block_sse2;
#endif
    if (flags & HAS_AVX2) aom_sse = aom_sse_avx2;
    av1_build_compound_diffwtd_mask_d16 = av1_build_compound_diffwtd_mask_d16_c;
#if COMP_HBD
     aom_highbd_sse = aom_highbd_sse_c;
     if (flags & HAS_AVX2) aom_highbd_sse = aom_highbd_sse_avx2;
#endif
    if (flags & HAS_AVX2) av1_build_compound_diffwtd_mask_d16 = av1_build_compound_diffwtd_mask_d16_avx2;
    aom_lowbd_blend_a64_d16_mask = aom_lowbd_blend_a64_d16_mask_c;
    if (flags & HAS_AVX2) aom_lowbd_blend_a64_d16_mask = aom_lowbd_blend_a64_d16_mask_avx2;
    aom_highbd_blend_a64_d16_mask = aom_highbd_blend_a64_d16_mask_c;
    if (flags & HAS_AVX2) aom_highbd_blend_a64_d16_mask = aom_highbd_blend_a64_d16_mask_avx2;

    av1_wedge_compute_delta_squares = av1_wedge_compute_delta_squares_c;
    if (flags & HAS_AVX2) av1_wedge_compute_delta_squares = av1_wedge_compute_delta_squares_avx2;
    av1_wedge_sign_from_residuals = av1_wedge_sign_from_residuals_c;
    if (flags & HAS_AVX2) av1_wedge_sign_from_residuals = av1_wedge_sign_from_residuals_avx2;

    eb_cdef_find_dir = eb_cdef_find_dir_c;
    if (flags & HAS_AVX2) eb_cdef_find_dir = eb_cdef_find_dir_avx2;

    eb_cdef_filter_block = eb_cdef_filter_block_c;
    if (flags & HAS_AVX2) eb_cdef_filter_block = eb_cdef_filter_block_avx2;
    eb_compute_cdef_dist = compute_cdef_dist_c;
    if (flags & HAS_AVX2) eb_compute_cdef_dist = compute_cdef_dist_avx2;
    eb_compute_cdef_dist_8bit = compute_cdef_dist_8bit_c;
    if (flags & HAS_AVX2) eb_compute_cdef_dist_8bit = compute_cdef_dist_8bit_avx2;

#if PREDICT_NSQ_SHAPE
    spatial_full_distortion = spatial_full_distortion_helper;
    if (flags & HAS_AVX2) spatial_full_distortion = spatial_full_distortion_avx2_helper;
#endif

    eb_copy_rect8_8bit_to_16bit = eb_copy_rect8_8bit_to_16bit_c;
    if (flags & HAS_AVX2) eb_copy_rect8_8bit_to_16bit = eb_copy_rect8_8bit_to_16bit_avx2;

    eb_av1_compute_stats = eb_av1_compute_stats_c;
    if (flags & HAS_AVX2) eb_av1_compute_stats = eb_av1_compute_stats_avx2;
    eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_c;
    if (flags & HAS_AVX2) eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_avx2;
    eb_cdef_filter_block_8x8_16 = eb_cdef_filter_block_8x8_16_avx2; // It has no c version, and is only called in parent avx2 function, so it's safe to initialize to avx2 version.
#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_cdef_filter_block_8x8_16 = eb_cdef_filter_block_8x8_16_avx512;
        eb_av1_compute_stats = eb_av1_compute_stats_avx512;
        eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_avx512;
    }
#endif

    eb_av1_lowbd_pixel_proj_error = eb_av1_lowbd_pixel_proj_error_c;
    eb_av1_highbd_pixel_proj_error = eb_av1_highbd_pixel_proj_error_c;
    if (flags & HAS_AVX2) eb_av1_highbd_pixel_proj_error = eb_av1_highbd_pixel_proj_error_avx2;
    eb_av1_filter_intra_edge_high = eb_av1_filter_intra_edge_high_c;
    if (flags & HAS_SSE4_1) eb_av1_filter_intra_edge_high = eb_av1_filter_intra_edge_high_sse4_1;
    eb_av1_calc_frame_error = eb_av1_calc_frame_error_c;
    if (flags & HAS_AVX2) eb_av1_calc_frame_error = eb_av1_calc_frame_error_avx2;
    eb_av1_highbd_convolve_2d_copy_sr = eb_av1_highbd_convolve_2d_copy_sr_c;
    if (flags & HAS_AVX2) eb_av1_highbd_convolve_2d_copy_sr = eb_av1_highbd_convolve_2d_copy_sr_avx2;
    eb_av1_highbd_jnt_convolve_2d_copy = eb_av1_highbd_jnt_convolve_2d_copy_c;
    if (flags & HAS_AVX2) eb_av1_highbd_jnt_convolve_2d_copy = eb_av1_highbd_jnt_convolve_2d_copy_avx2;
    eb_av1_highbd_convolve_y_sr = eb_av1_highbd_convolve_y_sr_c;
    if (flags & HAS_AVX2) eb_av1_highbd_convolve_y_sr = eb_av1_highbd_convolve_y_sr_avx2;
    eb_av1_highbd_convolve_2d_sr = eb_av1_highbd_convolve_2d_sr_c;
    if (flags & HAS_AVX2) eb_av1_highbd_convolve_2d_sr = eb_av1_highbd_convolve_2d_sr_avx2;

    eb_av1_highbd_convolve_2d_scale = eb_av1_highbd_convolve_2d_scale_c;
    //if (flags & HAS_SSE4_1) eb_av1_highbd_convolve_2d_scale = eb_av1_highbd_convolve_2d_scale_sse4_1

    eb_av1_highbd_jnt_convolve_2d = eb_av1_highbd_jnt_convolve_2d_c;
    if (flags & HAS_AVX2) eb_av1_highbd_jnt_convolve_2d = eb_av1_highbd_jnt_convolve_2d_avx2;
    eb_av1_highbd_jnt_convolve_x = eb_av1_highbd_jnt_convolve_x_c;
    if (flags & HAS_AVX2) eb_av1_highbd_jnt_convolve_x = eb_av1_highbd_jnt_convolve_x_avx2;
    eb_av1_highbd_jnt_convolve_y = eb_av1_highbd_jnt_convolve_y_c;
    if (flags & HAS_AVX2) eb_av1_highbd_jnt_convolve_y = eb_av1_highbd_jnt_convolve_y_avx2;
    eb_av1_highbd_convolve_x_sr = eb_av1_highbd_convolve_x_sr_c;
    if (flags & HAS_AVX2) eb_av1_highbd_convolve_x_sr = eb_av1_highbd_convolve_x_sr_avx2;
    eb_subtract_average = eb_subtract_average_c;
    if (flags & HAS_AVX2) eb_subtract_average = eb_subtract_average_avx2;

    eb_av1_filter_intra_edge = eb_av1_filter_intra_edge_high_c_old;

    if (flags & HAS_SSE4_1) eb_av1_filter_intra_edge = eb_av1_filter_intra_edge_sse4_1;

    eb_smooth_v_predictor = smooth_v_predictor_c;
    if (flags & HAS_SSSE3) eb_smooth_v_predictor = eb_smooth_v_predictor_all_ssse3;

    eb_smooth_h_predictor = smooth_h_predictor_c;
    if (flags & HAS_SSSE3) eb_smooth_h_predictor = eb_smooth_h_predictor_all_ssse3;

    get_proj_subspace = get_proj_subspace_c;
    if (flags & HAS_AVX2) get_proj_subspace = get_proj_subspace_avx2;

    search_one_dual = search_one_dual_c;
    if (flags & HAS_AVX2) search_one_dual = search_one_dual_avx2;

    eb_aom_mse16x16 = eb_aom_mse16x16_c;
    if (flags & HAS_AVX2) eb_aom_mse16x16 = eb_aom_mse16x16_avx2;

    eb_av1_convolve_2d_copy_sr = eb_av1_convolve_2d_copy_sr_c;
    if (flags & HAS_AVX2) eb_av1_convolve_2d_copy_sr = eb_av1_convolve_2d_copy_sr_avx2;

    eb_av1_convolve_2d_scale = eb_av1_convolve_2d_scale_c;
    //if (flags & HAS_SSE4_1) eb_av1_convolve_2d_scale = eb_av1_convolve_2d_scale_sse4_1;

    eb_aom_quantize_b = eb_aom_quantize_b_c_II;
    if (flags & HAS_AVX2) eb_aom_quantize_b = eb_aom_quantize_b_avx2;

    eb_aom_quantize_b_32x32 = eb_aom_quantize_b_32x32_c_II;
    if (flags & HAS_AVX2) eb_aom_quantize_b_32x32 = eb_aom_quantize_b_32x32_avx2;

    eb_aom_highbd_quantize_b_32x32 = eb_aom_highbd_quantize_b_32x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_quantize_b_32x32 = eb_aom_highbd_quantize_b_32x32_avx2;

    eb_aom_highbd_quantize_b = eb_aom_highbd_quantize_b_c;
    if (flags & HAS_AVX2) eb_aom_highbd_quantize_b = eb_aom_highbd_quantize_b_avx2;

    eb_av1_inv_txfm2d_add_16x16 = eb_av1_inv_txfm2d_add_16x16_c;
    eb_av1_inv_txfm2d_add_32x32 = eb_av1_inv_txfm2d_add_32x32_c;
    eb_av1_inv_txfm2d_add_4x4 = eb_av1_inv_txfm2d_add_4x4_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_4x4 = eb_av1_inv_txfm2d_add_4x4_avx2;
    eb_av1_inv_txfm2d_add_64x64 = eb_av1_inv_txfm2d_add_64x64_c;
    eb_av1_inv_txfm2d_add_8x8 = eb_av1_inv_txfm2d_add_8x8_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x8 = eb_av1_inv_txfm2d_add_8x8_avx2;

    eb_av1_inv_txfm2d_add_8x16 = eb_av1_inv_txfm2d_add_8x16_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x16 = eb_av1_highbd_inv_txfm_add_avx2;
    eb_av1_inv_txfm2d_add_16x8 = eb_av1_inv_txfm2d_add_16x8_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x8 = eb_av1_highbd_inv_txfm_add_avx2;
    eb_av1_inv_txfm2d_add_16x32 = eb_av1_inv_txfm2d_add_16x32_c;
    eb_av1_inv_txfm2d_add_32x16 = eb_av1_inv_txfm2d_add_32x16_c;
    eb_av1_inv_txfm2d_add_32x8 = eb_av1_inv_txfm2d_add_32x8_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x8 = eb_av1_highbd_inv_txfm_add_avx2;
    eb_av1_inv_txfm2d_add_8x32 = eb_av1_inv_txfm2d_add_8x32_c;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_8x32 = eb_av1_highbd_inv_txfm_add_avx2;
    eb_av1_inv_txfm2d_add_32x64 = eb_av1_inv_txfm2d_add_32x64_c;
    eb_av1_inv_txfm2d_add_64x32 = eb_av1_inv_txfm2d_add_64x32_c;
    eb_av1_inv_txfm2d_add_16x64 = eb_av1_inv_txfm2d_add_16x64_c;
    eb_av1_inv_txfm2d_add_64x16 = eb_av1_inv_txfm2d_add_64x16_c;
    eb_av1_inv_txfm2d_add_4x8 = eb_av1_inv_txfm2d_add_4x8_c;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_4x8 = eb_av1_inv_txfm2d_add_4x8_sse4_1;
    eb_av1_inv_txfm2d_add_8x4 = eb_av1_inv_txfm2d_add_8x4_c;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_8x4 = eb_av1_inv_txfm2d_add_8x4_sse4_1;
    eb_av1_inv_txfm2d_add_4x16 = eb_av1_inv_txfm2d_add_4x16_c;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_4x16 = eb_av1_inv_txfm2d_add_4x16_sse4_1;
    eb_av1_inv_txfm2d_add_16x4 = eb_av1_inv_txfm2d_add_16x4_c;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_16x4 = eb_av1_inv_txfm2d_add_16x4_sse4_1;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_av1_inv_txfm2d_add_16x16 = eb_av1_inv_txfm2d_add_16x16_avx512;
        eb_av1_inv_txfm2d_add_32x32 = eb_av1_inv_txfm2d_add_32x32_avx512;
        eb_av1_inv_txfm2d_add_64x64 = eb_av1_inv_txfm2d_add_64x64_avx512;
        eb_av1_inv_txfm2d_add_16x64 = eb_av1_inv_txfm2d_add_16x64_avx512;
        eb_av1_inv_txfm2d_add_64x16 = eb_av1_inv_txfm2d_add_64x16_avx512;
        eb_av1_inv_txfm2d_add_32x64 = eb_av1_inv_txfm2d_add_32x64_avx512;
        eb_av1_inv_txfm2d_add_64x32 = eb_av1_inv_txfm2d_add_64x32_avx512;
        eb_av1_inv_txfm2d_add_16x32 = eb_av1_inv_txfm2d_add_16x32_avx512;
        eb_av1_inv_txfm2d_add_32x16 = eb_av1_inv_txfm2d_add_32x16_avx512;
        eb_av1_lowbd_pixel_proj_error = eb_av1_lowbd_pixel_proj_error_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x16 = eb_av1_inv_txfm2d_add_16x16_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x32 = eb_av1_inv_txfm2d_add_32x32_avx2;
    if (flags & HAS_SSE4_1) eb_av1_inv_txfm2d_add_64x64 = eb_av1_inv_txfm2d_add_64x64_sse4_1;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x64 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_64x16 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x64 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_64x32 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_16x32 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_inv_txfm2d_add_32x16 = eb_av1_highbd_inv_txfm_add_avx2;
    if (flags & HAS_AVX2) eb_av1_lowbd_pixel_proj_error = eb_av1_lowbd_pixel_proj_error_avx2;
#endif
    eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_c;
    if (flags & HAS_SSSE3) eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_ssse3;
    if (flags & HAS_AVX2) eb_av1_inv_txfm_add = eb_av1_inv_txfm_add_avx2;

    eb_av1_quantize_fp = eb_av1_quantize_fp_c;
    if (flags & HAS_AVX2) eb_av1_quantize_fp = eb_av1_quantize_fp_avx2;

    eb_av1_quantize_fp_32x32 = eb_av1_quantize_fp_32x32_c;
    if (flags & HAS_AVX2) eb_av1_quantize_fp_32x32 = eb_av1_quantize_fp_32x32_avx2;

    eb_av1_quantize_fp_64x64 = eb_av1_quantize_fp_64x64_c;
    if (flags & HAS_AVX2) eb_av1_quantize_fp_64x64 = eb_av1_quantize_fp_64x64_avx2;

    eb_av1_highbd_quantize_fp = eb_av1_highbd_quantize_fp_c;
    if (flags & HAS_AVX2) eb_av1_highbd_quantize_fp = eb_av1_highbd_quantize_fp_avx2;

    highbd_variance64 = highbd_variance64_c;
    if (flags & HAS_AVX2) highbd_variance64 = highbd_variance64_avx2;

    eb_aom_highbd_8_mse16x16 = eb_aom_highbd_8_mse16x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_8_mse16x16 = eb_aom_highbd_8_mse16x16_sse2;

    eb_av1_upsample_intra_edge = eb_av1_upsample_intra_edge_c;
    if (flags & HAS_SSE4_1) eb_av1_upsample_intra_edge = eb_av1_upsample_intra_edge_sse4_1;

    eb_av1_warp_affine = eb_av1_warp_affine_c;
    if (flags & HAS_AVX2) eb_av1_warp_affine = eb_av1_warp_affine_avx2;

    eb_av1_filter_intra_predictor = eb_av1_filter_intra_predictor_c;
#if FILTER_INTRA_FLAG
    if (flags & HAS_SSE4_1) eb_av1_filter_intra_predictor = eb_av1_filter_intra_predictor_sse4_1;
#endif

    eb_aom_highbd_smooth_v_predictor_16x16 = eb_aom_highbd_smooth_v_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x16 = eb_aom_highbd_smooth_v_predictor_16x16_avx2;
    eb_aom_highbd_smooth_v_predictor_16x32 = eb_aom_highbd_smooth_v_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x32 = eb_aom_highbd_smooth_v_predictor_16x32_avx2;
    eb_aom_highbd_smooth_v_predictor_16x4 = eb_aom_highbd_smooth_v_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x4 = eb_aom_highbd_smooth_v_predictor_16x4_avx2;
    eb_aom_highbd_smooth_v_predictor_16x64 = eb_aom_highbd_smooth_v_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x64 = eb_aom_highbd_smooth_v_predictor_16x64_avx2;
    eb_aom_highbd_smooth_v_predictor_16x8 = eb_aom_highbd_smooth_v_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_16x8 = eb_aom_highbd_smooth_v_predictor_16x8_avx2;
    eb_aom_highbd_smooth_v_predictor_2x2 = eb_aom_highbd_smooth_v_predictor_2x2_c;
    eb_aom_highbd_smooth_v_predictor_32x16 = eb_aom_highbd_smooth_v_predictor_32x16_c;
    eb_aom_highbd_smooth_v_predictor_32x32 = eb_aom_highbd_smooth_v_predictor_32x32_c;
    eb_aom_highbd_smooth_v_predictor_32x64 = eb_aom_highbd_smooth_v_predictor_32x64_c;
    eb_aom_highbd_smooth_v_predictor_32x8 = eb_aom_highbd_smooth_v_predictor_32x8_c;
    eb_aom_highbd_smooth_v_predictor_4x16 = eb_aom_highbd_smooth_v_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x16 = eb_aom_highbd_smooth_v_predictor_4x16_ssse3;
    eb_aom_highbd_smooth_v_predictor_4x4 = eb_aom_highbd_smooth_v_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x4 = eb_aom_highbd_smooth_v_predictor_4x4_ssse3;
    eb_aom_highbd_smooth_v_predictor_4x8 = eb_aom_highbd_smooth_v_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x8 = eb_aom_highbd_smooth_v_predictor_4x8_ssse3;
    eb_aom_highbd_smooth_v_predictor_64x16 = eb_aom_highbd_smooth_v_predictor_64x16_c;
    eb_aom_highbd_smooth_v_predictor_64x32 = eb_aom_highbd_smooth_v_predictor_64x32_c;
    eb_aom_highbd_smooth_v_predictor_64x64 = eb_aom_highbd_smooth_v_predictor_64x64_c;
    eb_aom_highbd_smooth_v_predictor_8x16 = eb_aom_highbd_smooth_v_predictor_8x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x16 = eb_aom_highbd_smooth_v_predictor_8x16_avx2;
    eb_aom_highbd_smooth_v_predictor_8x32 = eb_aom_highbd_smooth_v_predictor_8x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x32 = eb_aom_highbd_smooth_v_predictor_8x32_avx2;
    eb_aom_highbd_smooth_v_predictor_8x4 = eb_aom_highbd_smooth_v_predictor_8x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x4 = eb_aom_highbd_smooth_v_predictor_8x4_avx2;
    eb_aom_highbd_smooth_v_predictor_8x8 = eb_aom_highbd_smooth_v_predictor_8x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x8 = eb_aom_highbd_smooth_v_predictor_8x8_avx2;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_highbd_smooth_v_predictor_32x8 = aom_highbd_smooth_v_predictor_32x8_avx512;
        eb_aom_highbd_smooth_v_predictor_32x16 = aom_highbd_smooth_v_predictor_32x16_avx512;
        eb_aom_highbd_smooth_v_predictor_32x32 = aom_highbd_smooth_v_predictor_32x32_avx512;
        eb_aom_highbd_smooth_v_predictor_32x64 = aom_highbd_smooth_v_predictor_32x64_avx512;
        eb_aom_highbd_smooth_v_predictor_64x16 = aom_highbd_smooth_v_predictor_64x16_avx512;
        eb_aom_highbd_smooth_v_predictor_64x32 = aom_highbd_smooth_v_predictor_64x32_avx512;
        eb_aom_highbd_smooth_v_predictor_64x64 = aom_highbd_smooth_v_predictor_64x64_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x8 = eb_aom_highbd_smooth_v_predictor_32x8_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x16 = eb_aom_highbd_smooth_v_predictor_32x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x32 = eb_aom_highbd_smooth_v_predictor_32x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x64 = eb_aom_highbd_smooth_v_predictor_32x64_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x16 = eb_aom_highbd_smooth_v_predictor_64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x32 = eb_aom_highbd_smooth_v_predictor_64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x64 = eb_aom_highbd_smooth_v_predictor_64x64_avx2;
#endif // !NON_AVX512_SUPPORT

    eb_cfl_predict_lbd = eb_cfl_predict_lbd_c;
    if (flags & HAS_AVX2) eb_cfl_predict_lbd = eb_cfl_predict_lbd_avx2;
    eb_cfl_predict_hbd = eb_cfl_predict_hbd_c;
    if (flags & HAS_AVX2) eb_cfl_predict_hbd = eb_cfl_predict_hbd_avx2;

    eb_av1_dr_prediction_z1 = eb_av1_dr_prediction_z1_c;
    if (flags & HAS_AVX2) eb_av1_dr_prediction_z1 = eb_av1_dr_prediction_z1_avx2;
    eb_av1_dr_prediction_z2 = eb_av1_dr_prediction_z2_c;
    if (flags & HAS_AVX2) eb_av1_dr_prediction_z2 = eb_av1_dr_prediction_z2_avx2;
    eb_av1_dr_prediction_z3 = eb_av1_dr_prediction_z3_c;
    if (flags & HAS_AVX2) eb_av1_dr_prediction_z3 = eb_av1_dr_prediction_z3_avx2;
    eb_av1_highbd_dr_prediction_z1 = eb_av1_highbd_dr_prediction_z1_c;
    if (flags & HAS_AVX2) eb_av1_highbd_dr_prediction_z1 = eb_av1_highbd_dr_prediction_z1_avx2;
    eb_av1_highbd_dr_prediction_z2 = eb_av1_highbd_dr_prediction_z2_c;
    if (flags & HAS_AVX2) eb_av1_highbd_dr_prediction_z2 = eb_av1_highbd_dr_prediction_z2_avx2;
    eb_av1_highbd_dr_prediction_z3 = eb_av1_highbd_dr_prediction_z3_c;
    if (flags & HAS_AVX2) eb_av1_highbd_dr_prediction_z3 = eb_av1_highbd_dr_prediction_z3_avx2;
    eb_av1_get_nz_map_contexts = eb_av1_get_nz_map_contexts_c;
    if (flags & HAS_SSE2) eb_av1_get_nz_map_contexts = eb_av1_get_nz_map_contexts_sse2;

#if II_COMP_FLAG
    aom_blend_a64_mask = aom_blend_a64_mask_c;
    if (flags & HAS_SSE4_1) aom_blend_a64_mask = aom_blend_a64_mask_sse4_1;
    if (flags & HAS_AVX2) aom_blend_a64_mask = aom_blend_a64_mask_avx2;
#endif //II_COMP_FLAG
    aom_blend_a64_hmask = aom_blend_a64_hmask_c;
    if (flags & HAS_SSE4_1) aom_blend_a64_hmask = aom_blend_a64_hmask_sse4_1;
    aom_blend_a64_vmask = aom_blend_a64_vmask_c;
    if (flags & HAS_SSE4_1) aom_blend_a64_vmask = aom_blend_a64_vmask_sse4_1;

    aom_highbd_blend_a64_mask = aom_highbd_blend_a64_mask_c;
    if (flags & HAS_SSE4_1) aom_highbd_blend_a64_mask = aom_highbd_blend_a64_mask_sse4_1;
    aom_highbd_blend_a64_hmask = aom_highbd_blend_a64_hmask_c;
    if (flags & HAS_SSE4_1) aom_highbd_blend_a64_hmask = aom_highbd_blend_a64_hmask_sse4_1;
    aom_highbd_blend_a64_vmask = aom_highbd_blend_a64_vmask_c;
    if (flags & HAS_SSE4_1) aom_highbd_blend_a64_vmask = aom_highbd_blend_a64_vmask_sse4_1;

    eb_aom_paeth_predictor_16x16 = eb_aom_paeth_predictor_16x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x16 = eb_aom_paeth_predictor_16x16_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_16x16 = eb_aom_paeth_predictor_16x16_avx2;
    eb_aom_paeth_predictor_16x32 = eb_aom_paeth_predictor_16x32_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x32 = eb_aom_paeth_predictor_16x32_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_16x32 = eb_aom_paeth_predictor_16x32_avx2;
    eb_aom_paeth_predictor_16x4 = eb_aom_paeth_predictor_16x4_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x4 = eb_aom_paeth_predictor_16x4_ssse3;
    eb_aom_paeth_predictor_16x64 = eb_aom_paeth_predictor_16x64_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x64 = eb_aom_paeth_predictor_16x64_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_16x64 = eb_aom_paeth_predictor_16x64_avx2;
    eb_aom_paeth_predictor_16x8 = eb_aom_paeth_predictor_16x8_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_16x8 = eb_aom_paeth_predictor_16x8_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_16x8 = eb_aom_paeth_predictor_16x8_avx2;
    eb_aom_paeth_predictor_32x16 = eb_aom_paeth_predictor_32x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_32x16 = eb_aom_paeth_predictor_32x16_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_32x16 = eb_aom_paeth_predictor_32x16_avx2;
    eb_aom_paeth_predictor_32x32 = eb_aom_paeth_predictor_32x32_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_32x32 = eb_aom_paeth_predictor_32x32_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_32x32 = eb_aom_paeth_predictor_32x32_avx2;
    eb_aom_paeth_predictor_32x64 = eb_aom_paeth_predictor_32x64_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_32x64 = eb_aom_paeth_predictor_32x64_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_32x64 = eb_aom_paeth_predictor_32x64_avx2;
    eb_aom_paeth_predictor_32x8 = eb_aom_paeth_predictor_32x8_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_32x8 = eb_aom_paeth_predictor_32x8_ssse3;
    eb_aom_paeth_predictor_4x16 = eb_aom_paeth_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_4x16 = eb_aom_paeth_predictor_4x16_ssse3;
    eb_aom_paeth_predictor_4x4 = eb_aom_paeth_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_4x4 = eb_aom_paeth_predictor_4x4_ssse3;
    eb_aom_paeth_predictor_4x8 = eb_aom_paeth_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_4x8 = eb_aom_paeth_predictor_4x8_ssse3;
    eb_aom_paeth_predictor_64x16 = eb_aom_paeth_predictor_64x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_64x16 = eb_aom_paeth_predictor_64x16_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_64x16 = eb_aom_paeth_predictor_64x16_avx2;
    eb_aom_paeth_predictor_64x32 = eb_aom_paeth_predictor_64x32_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_64x32 = eb_aom_paeth_predictor_64x32_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_64x32 = eb_aom_paeth_predictor_64x32_avx2;
    eb_aom_paeth_predictor_64x64 = eb_aom_paeth_predictor_64x64_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_64x64 = eb_aom_paeth_predictor_64x64_ssse3;
    if (flags & HAS_AVX2) eb_aom_paeth_predictor_64x64 = eb_aom_paeth_predictor_64x64_avx2;
    eb_aom_paeth_predictor_8x16 = eb_aom_paeth_predictor_8x16_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_8x16 = eb_aom_paeth_predictor_8x16_ssse3;
    eb_aom_paeth_predictor_8x32 = eb_aom_paeth_predictor_8x32_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_8x32 = eb_aom_paeth_predictor_8x32_ssse3;
    eb_aom_paeth_predictor_8x4 = eb_aom_paeth_predictor_8x4_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_8x4 = eb_aom_paeth_predictor_8x4_ssse3;
    eb_aom_paeth_predictor_8x8 = eb_aom_paeth_predictor_8x8_c;
    if (flags & HAS_SSSE3) eb_aom_paeth_predictor_8x8 = eb_aom_paeth_predictor_8x8_ssse3;

    eb_aom_highbd_paeth_predictor_16x16 = eb_aom_highbd_paeth_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x16 = eb_aom_highbd_paeth_predictor_16x16_avx2;
    eb_aom_highbd_paeth_predictor_16x32 = eb_aom_highbd_paeth_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x32 = eb_aom_highbd_paeth_predictor_16x32_avx2;
    eb_aom_highbd_paeth_predictor_16x4 = eb_aom_highbd_paeth_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x4 = eb_aom_highbd_paeth_predictor_16x4_avx2;
    eb_aom_highbd_paeth_predictor_16x64 = eb_aom_highbd_paeth_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x64 = eb_aom_highbd_paeth_predictor_16x64_avx2;
    eb_aom_highbd_paeth_predictor_16x8 = eb_aom_highbd_paeth_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_16x8 = eb_aom_highbd_paeth_predictor_16x8_avx2;
    eb_aom_highbd_paeth_predictor_2x2 = eb_aom_highbd_paeth_predictor_2x2_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_2x2 = eb_aom_highbd_paeth_predictor_2x2_avx2;
    eb_aom_highbd_paeth_predictor_32x16 = eb_aom_highbd_paeth_predictor_32x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_32x16 = eb_aom_highbd_paeth_predictor_32x16_avx2;
    eb_aom_highbd_paeth_predictor_32x32 = eb_aom_highbd_paeth_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_32x32 = eb_aom_highbd_paeth_predictor_32x32_avx2;
    eb_aom_highbd_paeth_predictor_32x64 = eb_aom_highbd_paeth_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_32x64 = eb_aom_highbd_paeth_predictor_32x64_avx2;
    eb_aom_highbd_paeth_predictor_32x8 = eb_aom_highbd_paeth_predictor_32x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_32x8 = eb_aom_highbd_paeth_predictor_32x8_avx2;
    eb_aom_highbd_paeth_predictor_4x16 = eb_aom_highbd_paeth_predictor_4x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_4x16 = eb_aom_highbd_paeth_predictor_4x16_avx2;
    eb_aom_highbd_paeth_predictor_4x4 = eb_aom_highbd_paeth_predictor_4x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_4x4 = eb_aom_highbd_paeth_predictor_4x4_avx2;
    eb_aom_highbd_paeth_predictor_4x8 = eb_aom_highbd_paeth_predictor_4x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_4x8 = eb_aom_highbd_paeth_predictor_4x8_avx2;
    eb_aom_highbd_paeth_predictor_64x16 = eb_aom_highbd_paeth_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_64x16 = eb_aom_highbd_paeth_predictor_64x16_avx2;
    eb_aom_highbd_paeth_predictor_64x32 = eb_aom_highbd_paeth_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_64x32 = eb_aom_highbd_paeth_predictor_64x32_avx2;
    eb_aom_highbd_paeth_predictor_64x64 = eb_aom_highbd_paeth_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_64x64 = eb_aom_highbd_paeth_predictor_64x64_avx2;
    eb_aom_highbd_paeth_predictor_8x16 = eb_aom_highbd_paeth_predictor_8x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_8x16 = eb_aom_highbd_paeth_predictor_8x16_avx2;
    eb_aom_highbd_paeth_predictor_8x32 = eb_aom_highbd_paeth_predictor_8x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_8x32 = eb_aom_highbd_paeth_predictor_8x32_avx2;
    eb_aom_highbd_paeth_predictor_8x4 = eb_aom_highbd_paeth_predictor_8x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_8x4 = eb_aom_highbd_paeth_predictor_8x4_avx2;
    eb_aom_highbd_paeth_predictor_8x8 = eb_aom_highbd_paeth_predictor_8x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_paeth_predictor_8x8 = eb_aom_highbd_paeth_predictor_8x8_avx2;

    eb_aom_dc_predictor_4x4 = eb_aom_dc_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_4x4 = eb_aom_dc_predictor_4x4_sse2;
    eb_aom_dc_predictor_8x8 = eb_aom_dc_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_8x8 = eb_aom_dc_predictor_8x8_sse2;
    eb_aom_dc_predictor_16x16 = eb_aom_dc_predictor_16x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_16x16 = eb_aom_dc_predictor_16x16_sse2;
    eb_aom_dc_predictor_32x32 = eb_aom_dc_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_dc_predictor_32x32 = eb_aom_dc_predictor_32x32_avx2;
    eb_aom_dc_predictor_64x64 = eb_aom_dc_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_dc_predictor_64x64 = eb_aom_dc_predictor_64x64_avx2;
    eb_aom_dc_predictor_32x16 = eb_aom_dc_predictor_32x16_c;
    if (flags & HAS_AVX2) eb_aom_dc_predictor_32x16 = eb_aom_dc_predictor_32x16_avx2;
    eb_aom_dc_predictor_32x64 = eb_aom_dc_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_dc_predictor_32x64 = eb_aom_dc_predictor_32x64_avx2;
    eb_aom_dc_predictor_64x16 = eb_aom_dc_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_dc_predictor_64x16 = eb_aom_dc_predictor_64x16_avx2;
    eb_aom_dc_predictor_8x16 = eb_aom_dc_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_8x16 = eb_aom_dc_predictor_8x16_sse2;
    eb_aom_dc_predictor_8x32 = eb_aom_dc_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_8x32 = eb_aom_dc_predictor_8x32_sse2;
    eb_aom_dc_predictor_8x4 = eb_aom_dc_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_8x4 = eb_aom_dc_predictor_8x4_sse2;
    eb_aom_dc_predictor_64x32 = eb_aom_dc_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_dc_predictor_64x32 = eb_aom_dc_predictor_64x32_avx2;
    eb_aom_dc_predictor_16x32 = eb_aom_dc_predictor_16x32_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_16x32 = eb_aom_dc_predictor_16x32_sse2;
    eb_aom_dc_predictor_16x4 = eb_aom_dc_predictor_16x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_16x4 = eb_aom_dc_predictor_16x4_sse2;
    eb_aom_dc_predictor_16x64 = eb_aom_dc_predictor_16x64_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_16x64 = eb_aom_dc_predictor_16x64_sse2;
    eb_aom_dc_predictor_16x8 = eb_aom_dc_predictor_16x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_16x8 = eb_aom_dc_predictor_16x8_sse2;
    eb_aom_dc_predictor_32x8 = eb_aom_dc_predictor_32x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_32x8 = eb_aom_dc_predictor_32x8_sse2;
    eb_aom_dc_predictor_4x16 = eb_aom_dc_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_4x16 = eb_aom_dc_predictor_4x16_sse2;
    eb_aom_dc_predictor_4x8 = eb_aom_dc_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_predictor_4x8 = eb_aom_dc_predictor_4x8_sse2;

    eb_aom_dc_top_predictor_4x4 = eb_aom_dc_top_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_4x4 = eb_aom_dc_top_predictor_4x4_sse2;
    eb_aom_dc_top_predictor_8x8 = eb_aom_dc_top_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_8x8 = eb_aom_dc_top_predictor_8x8_sse2;
    eb_aom_dc_top_predictor_16x16 = eb_aom_dc_top_predictor_16x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x16 = eb_aom_dc_top_predictor_16x16_sse2;
    eb_aom_dc_top_predictor_32x32 = eb_aom_dc_top_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_dc_top_predictor_32x32 = eb_aom_dc_top_predictor_32x32_avx2;
    eb_aom_dc_top_predictor_64x64 = eb_aom_dc_top_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_dc_top_predictor_64x64 = eb_aom_dc_top_predictor_64x64_avx2;
    eb_aom_dc_top_predictor_16x32 = eb_aom_dc_top_predictor_16x32_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x32 = eb_aom_dc_top_predictor_16x32_sse2;
    eb_aom_dc_top_predictor_16x4 = eb_aom_dc_top_predictor_16x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x4 = eb_aom_dc_top_predictor_16x4_sse2;
    eb_aom_dc_top_predictor_16x64 = eb_aom_dc_top_predictor_16x64_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x64 = eb_aom_dc_top_predictor_16x64_sse2;
    eb_aom_dc_top_predictor_16x8 = eb_aom_dc_top_predictor_16x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_16x8 = eb_aom_dc_top_predictor_16x8_sse2;
    eb_aom_dc_top_predictor_32x16 = eb_aom_dc_top_predictor_32x16_c;
    if (flags & HAS_AVX2) eb_aom_dc_top_predictor_32x16 = eb_aom_dc_top_predictor_32x16_avx2;
    eb_aom_dc_top_predictor_32x64 = eb_aom_dc_top_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_dc_top_predictor_32x64 = eb_aom_dc_top_predictor_32x64_avx2;
    eb_aom_dc_top_predictor_32x8 = eb_aom_dc_top_predictor_32x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_32x8 = eb_aom_dc_top_predictor_32x8_sse2;
    eb_aom_dc_top_predictor_4x16 = eb_aom_dc_top_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_4x16 = eb_aom_dc_top_predictor_4x16_sse2;
    eb_aom_dc_top_predictor_4x8 = eb_aom_dc_top_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_4x8 = eb_aom_dc_top_predictor_4x8_sse2;
    eb_aom_dc_top_predictor_64x16 = eb_aom_dc_top_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_dc_top_predictor_64x16 = eb_aom_dc_top_predictor_64x16_avx2;
    eb_aom_dc_top_predictor_64x32 = eb_aom_dc_top_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_dc_top_predictor_64x32 = eb_aom_dc_top_predictor_64x32_avx2;
    eb_aom_dc_top_predictor_8x16 = eb_aom_dc_top_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_8x16 = eb_aom_dc_top_predictor_8x16_sse2;
    eb_aom_dc_top_predictor_8x32 = eb_aom_dc_top_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_8x32 = eb_aom_dc_top_predictor_8x32_sse2;
    eb_aom_dc_top_predictor_8x4 = eb_aom_dc_top_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_top_predictor_8x4 = eb_aom_dc_top_predictor_8x4_sse2;

    eb_aom_dc_left_predictor_4x4 = eb_aom_dc_left_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_4x4 = eb_aom_dc_left_predictor_4x4_sse2;
    eb_aom_dc_left_predictor_8x8 = eb_aom_dc_left_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_8x8 = eb_aom_dc_left_predictor_8x8_sse2;
    eb_aom_dc_left_predictor_16x16 = eb_aom_dc_left_predictor_16x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_16x16 = eb_aom_dc_left_predictor_16x16_sse2;
    eb_aom_dc_left_predictor_32x32 = eb_aom_dc_left_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_dc_left_predictor_32x32 = eb_aom_dc_left_predictor_32x32_avx2;
    eb_aom_dc_left_predictor_64x64 = eb_aom_dc_left_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_dc_left_predictor_64x64 = eb_aom_dc_left_predictor_64x64_avx2;
    eb_aom_dc_left_predictor_16x32 = eb_aom_dc_left_predictor_16x32_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_16x32 = eb_aom_dc_left_predictor_16x32_sse2;
    eb_aom_dc_left_predictor_16x4 = eb_aom_dc_left_predictor_16x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_16x4 = eb_aom_dc_left_predictor_16x4_sse2;
    eb_aom_dc_left_predictor_16x64 = eb_aom_dc_left_predictor_16x64_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_16x64 = eb_aom_dc_left_predictor_16x64_sse2;
    eb_aom_dc_left_predictor_16x8 = eb_aom_dc_left_predictor_16x8_c;
    if (flags & HAS_SSE2)  eb_aom_dc_left_predictor_16x8 = eb_aom_dc_left_predictor_16x8_sse2;
    eb_aom_dc_left_predictor_32x16 = eb_aom_dc_left_predictor_32x16_c;
    if (flags & HAS_AVX2) eb_aom_dc_left_predictor_32x16 = eb_aom_dc_left_predictor_32x16_avx2;
    eb_aom_dc_left_predictor_32x64 = eb_aom_dc_left_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_dc_left_predictor_32x64 = eb_aom_dc_left_predictor_32x64_avx2;
    eb_aom_dc_left_predictor_64x16 = eb_aom_dc_left_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_dc_left_predictor_64x16 = eb_aom_dc_left_predictor_64x16_avx2;
    eb_aom_dc_left_predictor_64x32 = eb_aom_dc_left_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_dc_left_predictor_64x32 = eb_aom_dc_left_predictor_64x32_avx2;
    eb_aom_dc_left_predictor_32x8 = eb_aom_dc_left_predictor_32x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_32x8 = eb_aom_dc_left_predictor_32x8_sse2;
    eb_aom_dc_left_predictor_4x16 = eb_aom_dc_left_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_4x16 = eb_aom_dc_left_predictor_4x16_sse2;
    eb_aom_dc_left_predictor_4x8 = eb_aom_dc_left_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_4x8 = eb_aom_dc_left_predictor_4x8_sse2;
    eb_aom_dc_left_predictor_8x16 = eb_aom_dc_left_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_8x16 = eb_aom_dc_left_predictor_8x16_sse2;
    eb_aom_dc_left_predictor_8x32 = eb_aom_dc_left_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_8x32 = eb_aom_dc_left_predictor_8x32_sse2;
    eb_aom_dc_left_predictor_8x4 = eb_aom_dc_left_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_left_predictor_8x4 = eb_aom_dc_left_predictor_8x4_sse2;

    eb_aom_dc_128_predictor_4x4 = eb_aom_dc_128_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_4x4 = eb_aom_dc_128_predictor_4x4_sse2;
    eb_aom_dc_128_predictor_8x8 = eb_aom_dc_128_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_8x8 = eb_aom_dc_128_predictor_8x8_sse2;
    eb_aom_dc_128_predictor_16x16 = eb_aom_dc_128_predictor_16x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x16 = eb_aom_dc_128_predictor_16x16_sse2;
    eb_aom_dc_128_predictor_32x32 = eb_aom_dc_128_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_dc_128_predictor_32x32 = eb_aom_dc_128_predictor_32x32_avx2;
    eb_aom_dc_128_predictor_64x64 = eb_aom_dc_128_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_dc_128_predictor_64x64 = eb_aom_dc_128_predictor_64x64_avx2;
    eb_aom_dc_128_predictor_16x32 = eb_aom_dc_128_predictor_16x32_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x32 = eb_aom_dc_128_predictor_16x32_sse2;
    eb_aom_dc_128_predictor_16x4 = eb_aom_dc_128_predictor_16x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x4 = eb_aom_dc_128_predictor_16x4_sse2;
    eb_aom_dc_128_predictor_16x64 = eb_aom_dc_128_predictor_16x64_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x64 = eb_aom_dc_128_predictor_16x64_sse2;
    eb_aom_dc_128_predictor_16x8 = eb_aom_dc_128_predictor_16x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_16x8 = eb_aom_dc_128_predictor_16x8_sse2;
    eb_aom_dc_128_predictor_32x16 = eb_aom_dc_128_predictor_32x16_c;
    if (flags & HAS_AVX2) eb_aom_dc_128_predictor_32x16 = eb_aom_dc_128_predictor_32x16_avx2;
    eb_aom_dc_128_predictor_32x64 = eb_aom_dc_128_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_dc_128_predictor_32x64 = eb_aom_dc_128_predictor_32x64_avx2;
    eb_aom_dc_128_predictor_32x8 = eb_aom_dc_128_predictor_32x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_32x8 = eb_aom_dc_128_predictor_32x8_sse2;
    eb_aom_dc_128_predictor_4x16 = eb_aom_dc_128_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_4x16 = eb_aom_dc_128_predictor_4x16_sse2;
    eb_aom_dc_128_predictor_4x8 = eb_aom_dc_128_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_4x8 = eb_aom_dc_128_predictor_4x8_sse2;
    eb_aom_dc_128_predictor_64x16 = eb_aom_dc_128_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_dc_128_predictor_64x16 = eb_aom_dc_128_predictor_64x16_avx2;
    eb_aom_dc_128_predictor_64x32 = eb_aom_dc_128_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_dc_128_predictor_64x32 = eb_aom_dc_128_predictor_64x32_avx2;
    eb_aom_dc_128_predictor_8x16 = eb_aom_dc_128_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_8x16 = eb_aom_dc_128_predictor_8x16_sse2;
    eb_aom_dc_128_predictor_8x32 = eb_aom_dc_128_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_8x32 = eb_aom_dc_128_predictor_8x32_sse2;
    eb_aom_dc_128_predictor_8x4 = eb_aom_dc_128_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_dc_128_predictor_8x4 = eb_aom_dc_128_predictor_8x4_sse2;

    eb_aom_smooth_h_predictor_16x32 = eb_aom_smooth_h_predictor_16x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x32 = eb_aom_smooth_h_predictor_16x32_ssse3;
    eb_aom_smooth_h_predictor_16x4 = eb_aom_smooth_h_predictor_16x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x4 = eb_aom_smooth_h_predictor_16x4_ssse3;
    eb_aom_smooth_h_predictor_16x64 = eb_aom_smooth_h_predictor_16x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x64 = eb_aom_smooth_h_predictor_16x64_ssse3;
    eb_aom_smooth_h_predictor_16x8 = eb_aom_smooth_h_predictor_16x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x8 = eb_aom_smooth_h_predictor_16x8_ssse3;
    eb_aom_smooth_h_predictor_32x16 = eb_aom_smooth_h_predictor_32x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_32x16 = eb_aom_smooth_h_predictor_32x16_ssse3;
    eb_aom_smooth_h_predictor_32x64 = eb_aom_smooth_h_predictor_32x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_32x64 = eb_aom_smooth_h_predictor_32x64_ssse3;
    eb_aom_smooth_h_predictor_32x8 = eb_aom_smooth_h_predictor_32x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_32x8 = eb_aom_smooth_h_predictor_32x8_ssse3;
    eb_aom_smooth_h_predictor_4x16 = eb_aom_smooth_h_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_4x16 = eb_aom_smooth_h_predictor_4x16_ssse3;
    eb_aom_smooth_h_predictor_4x8 = eb_aom_smooth_h_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_4x8 = eb_aom_smooth_h_predictor_4x8_ssse3;
    eb_aom_smooth_h_predictor_64x16 = eb_aom_smooth_h_predictor_64x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_64x16 = eb_aom_smooth_h_predictor_64x16_ssse3;
    eb_aom_smooth_h_predictor_64x32 = eb_aom_smooth_h_predictor_64x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_64x32 = eb_aom_smooth_h_predictor_64x32_ssse3;
    eb_aom_smooth_h_predictor_8x16 = eb_aom_smooth_h_predictor_8x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_8x16 = eb_aom_smooth_h_predictor_8x16_ssse3;
    eb_aom_smooth_h_predictor_8x32 = eb_aom_smooth_h_predictor_8x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_8x32 = eb_aom_smooth_h_predictor_8x32_ssse3;
    eb_aom_smooth_h_predictor_8x4 = eb_aom_smooth_h_predictor_8x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_8x4 = eb_aom_smooth_h_predictor_8x4_ssse3;
    eb_aom_smooth_h_predictor_64x64 = eb_aom_smooth_h_predictor_64x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_64x64 = eb_aom_smooth_h_predictor_64x64_ssse3;
    eb_aom_smooth_h_predictor_32x32 = eb_aom_smooth_h_predictor_32x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_32x32 = eb_aom_smooth_h_predictor_32x32_ssse3;
    eb_aom_smooth_h_predictor_16x16 = eb_aom_smooth_h_predictor_16x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_16x16 = eb_aom_smooth_h_predictor_16x16_ssse3;
    eb_aom_smooth_h_predictor_8x8 = eb_aom_smooth_h_predictor_8x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_8x8 = eb_aom_smooth_h_predictor_8x8_ssse3;
    eb_aom_smooth_h_predictor_4x4 = eb_aom_smooth_h_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_h_predictor_4x4 = eb_aom_smooth_h_predictor_4x4_ssse3;
    eb_aom_smooth_v_predictor_16x32 = eb_aom_smooth_v_predictor_16x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x32 = eb_aom_smooth_v_predictor_16x32_ssse3;
    eb_aom_smooth_v_predictor_16x4 = eb_aom_smooth_v_predictor_16x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x4 = eb_aom_smooth_v_predictor_16x4_ssse3;
    eb_aom_smooth_v_predictor_16x64 = eb_aom_smooth_v_predictor_16x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x64 = eb_aom_smooth_v_predictor_16x64_ssse3;
    eb_aom_smooth_v_predictor_16x8 = eb_aom_smooth_v_predictor_16x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x8 = eb_aom_smooth_v_predictor_16x8_ssse3;
    eb_aom_smooth_v_predictor_32x16 = eb_aom_smooth_v_predictor_32x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_32x16 = eb_aom_smooth_v_predictor_32x16_ssse3;
    eb_aom_smooth_v_predictor_32x64 = eb_aom_smooth_v_predictor_32x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_32x64 = eb_aom_smooth_v_predictor_32x64_ssse3;
    eb_aom_smooth_v_predictor_32x8 = eb_aom_smooth_v_predictor_32x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_32x8 = eb_aom_smooth_v_predictor_32x8_ssse3;
    eb_aom_smooth_v_predictor_4x16 = eb_aom_smooth_v_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_4x16 = eb_aom_smooth_v_predictor_4x16_ssse3;
    eb_aom_smooth_v_predictor_4x8 = eb_aom_smooth_v_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_4x8 = eb_aom_smooth_v_predictor_4x8_ssse3;
    eb_aom_smooth_v_predictor_64x16 = eb_aom_smooth_v_predictor_64x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_64x16 = eb_aom_smooth_v_predictor_64x16_ssse3;
    eb_aom_smooth_v_predictor_64x32 = eb_aom_smooth_v_predictor_64x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_64x32 = eb_aom_smooth_v_predictor_64x32_ssse3;
    eb_aom_smooth_v_predictor_8x16 = eb_aom_smooth_v_predictor_8x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_8x16 = eb_aom_smooth_v_predictor_8x16_ssse3;
    eb_aom_smooth_v_predictor_8x32 = eb_aom_smooth_v_predictor_8x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_8x32 = eb_aom_smooth_v_predictor_8x32_ssse3;
    eb_aom_smooth_v_predictor_8x4 = eb_aom_smooth_v_predictor_8x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_8x4 = eb_aom_smooth_v_predictor_8x4_ssse3;
    eb_aom_smooth_v_predictor_64x64 = eb_aom_smooth_v_predictor_64x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_64x64 = eb_aom_smooth_v_predictor_64x64_ssse3;
    eb_aom_smooth_v_predictor_32x32 = eb_aom_smooth_v_predictor_32x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_32x32 = eb_aom_smooth_v_predictor_32x32_ssse3;
    eb_aom_smooth_v_predictor_16x16 = eb_aom_smooth_v_predictor_16x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_16x16 = eb_aom_smooth_v_predictor_16x16_ssse3;
    eb_aom_smooth_v_predictor_8x8 = eb_aom_smooth_v_predictor_8x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_8x8 = eb_aom_smooth_v_predictor_8x8_ssse3;
    eb_aom_smooth_v_predictor_4x4 = eb_aom_smooth_v_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_v_predictor_4x4 = eb_aom_smooth_v_predictor_4x4_ssse3;

    eb_aom_smooth_predictor_16x32 = eb_aom_smooth_predictor_16x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x32 = eb_aom_smooth_predictor_16x32_ssse3;
    eb_aom_smooth_predictor_16x4 = eb_aom_smooth_predictor_16x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x4 = eb_aom_smooth_predictor_16x4_ssse3;
    eb_aom_smooth_predictor_16x64 = eb_aom_smooth_predictor_16x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x64 = eb_aom_smooth_predictor_16x64_ssse3;
    eb_aom_smooth_predictor_16x8 = eb_aom_smooth_predictor_16x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x8 = eb_aom_smooth_predictor_16x8_ssse3;
    eb_aom_smooth_predictor_32x16 = eb_aom_smooth_predictor_32x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_32x16 = eb_aom_smooth_predictor_32x16_ssse3;
    eb_aom_smooth_predictor_32x64 = eb_aom_smooth_predictor_32x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_32x64 = eb_aom_smooth_predictor_32x64_ssse3;
    eb_aom_smooth_predictor_32x8 = eb_aom_smooth_predictor_32x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_32x8 = eb_aom_smooth_predictor_32x8_ssse3;
    eb_aom_smooth_predictor_4x16 = eb_aom_smooth_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_4x16 = eb_aom_smooth_predictor_4x16_ssse3;
    eb_aom_smooth_predictor_4x8 = eb_aom_smooth_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_4x8 = eb_aom_smooth_predictor_4x8_ssse3;
    eb_aom_smooth_predictor_64x16 = eb_aom_smooth_predictor_64x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_64x16 = eb_aom_smooth_predictor_64x16_ssse3;
    eb_aom_smooth_predictor_64x32 = eb_aom_smooth_predictor_64x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_64x32 = eb_aom_smooth_predictor_64x32_ssse3;
    eb_aom_smooth_predictor_8x16 = eb_aom_smooth_predictor_8x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_8x16 = eb_aom_smooth_predictor_8x16_ssse3;
    eb_aom_smooth_predictor_8x32 = eb_aom_smooth_predictor_8x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_8x32 = eb_aom_smooth_predictor_8x32_ssse3;
    eb_aom_smooth_predictor_8x4 = eb_aom_smooth_predictor_8x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_8x4 = eb_aom_smooth_predictor_8x4_ssse3;
    eb_aom_smooth_predictor_64x64 = eb_aom_smooth_predictor_64x64_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_64x64 = eb_aom_smooth_predictor_64x64_ssse3;
    eb_aom_smooth_predictor_32x32 = eb_aom_smooth_predictor_32x32_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_32x32 = eb_aom_smooth_predictor_32x32_ssse3;
    eb_aom_smooth_predictor_16x16 = eb_aom_smooth_predictor_16x16_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_16x16 = eb_aom_smooth_predictor_16x16_ssse3;
    eb_aom_smooth_predictor_8x8 = eb_aom_smooth_predictor_8x8_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_8x8 = eb_aom_smooth_predictor_8x8_ssse3;
    eb_aom_smooth_predictor_4x4 = eb_aom_smooth_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_smooth_predictor_4x4 = eb_aom_smooth_predictor_4x4_ssse3;

    eb_aom_v_predictor_4x4 = eb_aom_v_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_4x4 = eb_aom_v_predictor_4x4_sse2;
    eb_aom_v_predictor_8x8 = eb_aom_v_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_8x8 = eb_aom_v_predictor_8x8_sse2;
    eb_aom_v_predictor_16x16 = eb_aom_v_predictor_16x16_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_16x16 = eb_aom_v_predictor_16x16_sse2;
    eb_aom_v_predictor_32x32 = eb_aom_v_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_v_predictor_32x32 = eb_aom_v_predictor_32x32_avx2;
    eb_aom_v_predictor_64x64 = eb_aom_v_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_v_predictor_64x64 = eb_aom_v_predictor_64x64_avx2;
    eb_aom_v_predictor_16x32 = eb_aom_v_predictor_16x32_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_16x32 = eb_aom_v_predictor_16x32_sse2;
    eb_aom_v_predictor_16x4 = eb_aom_v_predictor_16x4_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_16x4 = eb_aom_v_predictor_16x4_sse2;
    eb_aom_v_predictor_16x64 = eb_aom_v_predictor_16x64_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_16x64 = eb_aom_v_predictor_16x64_sse2;
    eb_aom_v_predictor_16x8 = eb_aom_v_predictor_16x8_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_16x8 = eb_aom_v_predictor_16x8_sse2;
    eb_aom_v_predictor_32x16 = eb_aom_v_predictor_32x16_c;
    if (flags & HAS_AVX2) eb_aom_v_predictor_32x16 = eb_aom_v_predictor_32x16_avx2;
    eb_aom_v_predictor_32x64 = eb_aom_v_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_v_predictor_32x64 = eb_aom_v_predictor_32x64_avx2;
    eb_aom_v_predictor_32x8 = eb_aom_v_predictor_32x8_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_32x8 = eb_aom_v_predictor_32x8_sse2;
    eb_aom_v_predictor_4x16 = eb_aom_v_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_4x16 = eb_aom_v_predictor_4x16_sse2;
    eb_aom_v_predictor_4x8 = eb_aom_v_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_4x8 = eb_aom_v_predictor_4x8_sse2;
    eb_aom_v_predictor_64x16 = eb_aom_v_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_v_predictor_64x16 = eb_aom_v_predictor_64x16_avx2;
    eb_aom_v_predictor_64x32 = eb_aom_v_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_v_predictor_64x32 = eb_aom_v_predictor_64x32_avx2;
    eb_aom_v_predictor_8x16 = eb_aom_v_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_8x16 = eb_aom_v_predictor_8x16_sse2;
    eb_aom_v_predictor_8x32 = eb_aom_v_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_8x32 = eb_aom_v_predictor_8x32_sse2;
    eb_aom_v_predictor_8x4 = eb_aom_v_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_v_predictor_8x4 = eb_aom_v_predictor_8x4_sse2;

    eb_aom_h_predictor_4x4 = eb_aom_h_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_4x4 = eb_aom_h_predictor_4x4_sse2;
    eb_aom_h_predictor_8x8 = eb_aom_h_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_8x8 = eb_aom_h_predictor_8x8_sse2;
    eb_aom_h_predictor_16x16 = eb_aom_h_predictor_16x16_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_16x16 = eb_aom_h_predictor_16x16_sse2;
    eb_aom_h_predictor_32x32 = eb_aom_h_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_h_predictor_32x32 = eb_aom_h_predictor_32x32_avx2;
    eb_aom_h_predictor_64x64 = eb_aom_h_predictor_64x64_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_64x64 = eb_aom_h_predictor_64x64_sse2;
    eb_aom_h_predictor_16x32 = eb_aom_h_predictor_16x32_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_16x32 = eb_aom_h_predictor_16x32_sse2;
    eb_aom_h_predictor_16x4 = eb_aom_h_predictor_16x4_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_16x4 = eb_aom_h_predictor_16x4_sse2;
    eb_aom_h_predictor_16x64 = eb_aom_h_predictor_16x64_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_16x64 = eb_aom_h_predictor_16x64_sse2;
    eb_aom_h_predictor_16x8 = eb_aom_h_predictor_16x8_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_16x8 = eb_aom_h_predictor_16x8_sse2;
    eb_aom_h_predictor_32x16 = eb_aom_h_predictor_32x16_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_32x16 = eb_aom_h_predictor_32x16_sse2;
    eb_aom_h_predictor_32x64 = eb_aom_h_predictor_32x64_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_32x64 = eb_aom_h_predictor_32x64_sse2;
    eb_aom_h_predictor_32x8 = eb_aom_h_predictor_32x8_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_32x8 = eb_aom_h_predictor_32x8_sse2;
    eb_aom_h_predictor_4x16 = eb_aom_h_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_4x16 = eb_aom_h_predictor_4x16_sse2;
    eb_aom_h_predictor_4x8 = eb_aom_h_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_4x8 = eb_aom_h_predictor_4x8_sse2;
    eb_aom_h_predictor_64x16 = eb_aom_h_predictor_64x16_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_64x16 = eb_aom_h_predictor_64x16_sse2;
    eb_aom_h_predictor_64x32 = eb_aom_h_predictor_64x32_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_64x32 = eb_aom_h_predictor_64x32_sse2;
    eb_aom_h_predictor_8x16 = eb_aom_h_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_8x16 = eb_aom_h_predictor_8x16_sse2;
    eb_aom_h_predictor_8x32 = eb_aom_h_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_8x32 = eb_aom_h_predictor_8x32_sse2;
    eb_aom_h_predictor_8x4 = eb_aom_h_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_h_predictor_8x4 = eb_aom_h_predictor_8x4_sse2;

    //SAD
    eb_aom_sad4x4 = eb_aom_sad4x4_c;
    if (flags & HAS_AVX2) eb_aom_sad4x4 = eb_aom_sad4x4_avx2;
    eb_aom_sad4x4x4d = eb_aom_sad4x4x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad4x4x4d = eb_aom_sad4x4x4d_avx2;
    eb_aom_sad4x16 = eb_aom_sad4x16_c;
    if (flags & HAS_AVX2) eb_aom_sad4x16 = eb_aom_sad4x16_avx2;
    eb_aom_sad4x16x4d = eb_aom_sad4x16x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad4x16x4d = eb_aom_sad4x16x4d_avx2;
    eb_aom_sad4x8 = eb_aom_sad4x8_c;
    if (flags & HAS_AVX2) eb_aom_sad4x8 = eb_aom_sad4x8_avx2;
    eb_aom_sad4x8x4d = eb_aom_sad4x8x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad4x8x4d = eb_aom_sad4x8x4d_avx2;
    eb_aom_sad64x128x4d = eb_aom_sad64x128x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad64x128x4d = eb_aom_sad64x128x4d_avx2;
    eb_aom_sad64x16x4d = eb_aom_sad64x16x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad64x16x4d = eb_aom_sad64x16x4d_avx2;
    eb_aom_sad64x32x4d = eb_aom_sad64x32x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad64x32x4d = eb_aom_sad64x32x4d_avx2;
    eb_aom_sad64x64x4d = eb_aom_sad64x64x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad64x64x4d = eb_aom_sad64x64x4d_avx2;
    eb_aom_sad8x16 = eb_aom_sad8x16_c;
    if (flags & HAS_AVX2) eb_aom_sad8x16 = eb_aom_sad8x16_avx2;
    eb_aom_sad8x16x4d = eb_aom_sad8x16x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad8x16x4d = eb_aom_sad8x16x4d_avx2;
    eb_aom_sad8x32 = eb_aom_sad8x32_c;
    if (flags & HAS_AVX2) eb_aom_sad8x32 = eb_aom_sad8x32_avx2;
    eb_aom_sad8x32x4d = eb_aom_sad8x32x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad8x32x4d = eb_aom_sad8x32x4d_avx2;
    eb_aom_sad8x8 = eb_aom_sad8x8_c;
    if (flags & HAS_AVX2) eb_aom_sad8x8 = eb_aom_sad8x8_avx2;
    eb_aom_sad8x8x4d = eb_aom_sad8x8x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad8x8x4d = eb_aom_sad8x8x4d_avx2;
    eb_aom_sad16x4 = eb_aom_sad16x4_c;
    if (flags & HAS_AVX2) eb_aom_sad16x4 = eb_aom_sad16x4_avx2;
    eb_aom_sad16x4x4d = eb_aom_sad16x4x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad16x4x4d = eb_aom_sad16x4x4d_avx2;
    eb_aom_sad32x8 = eb_aom_sad32x8_c;
    if (flags & HAS_AVX2) eb_aom_sad32x8 = eb_aom_sad32x8_avx2;
    eb_aom_sad32x8x4d = eb_aom_sad32x8x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad32x8x4d = eb_aom_sad32x8x4d_avx2;
    eb_aom_sad16x64 = eb_aom_sad16x64_c;
    if (flags & HAS_AVX2) eb_aom_sad16x64 = eb_aom_sad16x64_avx2;
    eb_aom_sad16x64x4d = eb_aom_sad16x64x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad16x64x4d = eb_aom_sad16x64x4d_avx2;
    eb_aom_sad32x16 = eb_aom_sad32x16_c;
    if (flags & HAS_AVX2) eb_aom_sad32x16 = eb_aom_sad32x16_avx2;
    eb_aom_sad32x16x4d = eb_aom_sad32x16x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad32x16x4d = eb_aom_sad32x16x4d_avx2;
    eb_aom_sad16x32 = eb_aom_sad16x32_c;
    if (flags & HAS_AVX2) eb_aom_sad16x32 = eb_aom_sad16x32_avx2;
    eb_aom_sad16x32x4d = eb_aom_sad16x32x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad16x32x4d = eb_aom_sad16x32x4d_avx2;
    eb_aom_sad32x64 = eb_aom_sad32x64_c;
    if (flags & HAS_AVX2) eb_aom_sad32x64 = eb_aom_sad32x64_avx2;
    eb_aom_sad32x64x4d = eb_aom_sad32x64x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad32x64x4d = eb_aom_sad32x64x4d_avx2;
    eb_aom_sad32x32 = eb_aom_sad32x32_c;
    if (flags & HAS_AVX2) eb_aom_sad32x32 = eb_aom_sad32x32_avx2;
    eb_aom_sad32x32x4d = eb_aom_sad32x32x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad32x32x4d = eb_aom_sad32x32x4d_avx2;
    eb_aom_sad16x16 = eb_aom_sad16x16_c;
    if (flags & HAS_AVX2) eb_aom_sad16x16 = eb_aom_sad16x16_avx2;
    eb_aom_sad16x16x4d = eb_aom_sad16x16x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad16x16x4d = eb_aom_sad16x16x4d_avx2;
    eb_aom_sad16x8 = eb_aom_sad16x8_c;
    if (flags & HAS_AVX2) eb_aom_sad16x8 = eb_aom_sad16x8_avx2;
    eb_aom_sad16x8x4d = eb_aom_sad16x8x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad16x8x4d = eb_aom_sad16x8x4d_avx2;
    eb_aom_sad8x4 = eb_aom_sad8x4_c;
    if (flags & HAS_AVX2) eb_aom_sad8x4 = eb_aom_sad8x4_avx2;
    eb_aom_sad8x4x4d = eb_aom_sad8x4x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad8x4x4d = eb_aom_sad8x4x4d_avx2;

    eb_aom_sad64x128 = eb_aom_sad64x128_c;
    eb_aom_sad64x16 = eb_aom_sad64x16_c;
    eb_aom_sad64x32 = eb_aom_sad64x32_c;
    eb_aom_sad64x64 = eb_aom_sad64x64_c;
    eb_aom_sad128x128 = eb_aom_sad128x128_c;
    eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_c;
    eb_aom_sad128x64 = eb_aom_sad128x64_c;
    eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_c;
    eb_av1_txb_init_levels = eb_av1_txb_init_levels_c;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_sad64x128 = eb_aom_sad64x128_avx512;
        eb_aom_sad64x16 = eb_aom_sad64x16_avx512;
        eb_aom_sad64x32 = eb_aom_sad64x32_avx512;
        eb_aom_sad64x64 = eb_aom_sad64x64_avx512;
        eb_aom_sad128x128 = eb_aom_sad128x128_avx512;
        eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_avx512;
        eb_aom_sad128x64 = eb_aom_sad128x64_avx512;
        eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_avx512;
        eb_av1_txb_init_levels = eb_av1_txb_init_levels_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_aom_sad64x128 = eb_aom_sad64x128_avx2;
    if (flags & HAS_AVX2) eb_aom_sad64x16 = eb_aom_sad64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_sad64x32 = eb_aom_sad64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_sad64x64 = eb_aom_sad64x64_avx2;
    if (flags & HAS_AVX2) eb_aom_sad128x128 = eb_aom_sad128x128_avx2;
    if (flags & HAS_AVX2) eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_avx2;
    if (flags & HAS_AVX2) eb_aom_sad128x64 = eb_aom_sad128x64_avx2;
    if (flags & HAS_AVX2) eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_avx2;
    if (flags & HAS_AVX2) eb_av1_txb_init_levels = eb_av1_txb_init_levels_avx2;
#endif // !NON_AVX512_SUPPORT

#if OBMC_FLAG
    eb_aom_highbd_blend_a64_vmask = eb_aom_highbd_blend_a64_vmask_c;
    if (flags & HAS_SSE4_1) eb_aom_highbd_blend_a64_vmask = eb_aom_highbd_blend_a64_vmask_sse4_1;
    eb_aom_highbd_blend_a64_hmask = eb_aom_highbd_blend_a64_hmask_c;
    if (flags & HAS_SSE4_1) eb_aom_highbd_blend_a64_hmask = eb_aom_highbd_blend_a64_hmask_sse4_1;

    aom_convolve8_horiz = aom_convolve8_horiz_c;
    if (flags & HAS_AVX2) aom_convolve8_horiz = aom_convolve8_horiz_avx2;
    aom_convolve8_vert = aom_convolve8_vert_c;
    if (flags & HAS_AVX2) aom_convolve8_vert = aom_convolve8_vert_avx2;
    aom_upsampled_pred = aom_upsampled_pred_c;
    if (flags & HAS_AVX2) aom_upsampled_pred = aom_upsampled_pred_sse2;

    aom_obmc_sad128x128 = aom_obmc_sad128x128_c;
    if (flags & HAS_AVX2) aom_obmc_sad128x128 = aom_obmc_sad128x128_avx2;
    aom_obmc_sad128x64 = aom_obmc_sad128x64_c;
    if (flags & HAS_AVX2) aom_obmc_sad128x64 = aom_obmc_sad128x64_avx2;
    aom_obmc_sad16x16 = aom_obmc_sad16x16_c;
    if (flags & HAS_AVX2) aom_obmc_sad16x16 = aom_obmc_sad16x16_avx2;
    aom_obmc_sad16x32 = aom_obmc_sad16x32_c;
    if (flags & HAS_AVX2) aom_obmc_sad16x32 = aom_obmc_sad16x32_avx2;
    aom_obmc_sad16x4 = aom_obmc_sad16x4_c;
    if (flags & HAS_AVX2) aom_obmc_sad16x4 = aom_obmc_sad16x4_avx2;
    aom_obmc_sad16x64 = aom_obmc_sad16x64_c;
    if (flags & HAS_AVX2) aom_obmc_sad16x64 = aom_obmc_sad16x64_avx2;
    aom_obmc_sad16x8 = aom_obmc_sad16x8_c;
    if (flags & HAS_AVX2) aom_obmc_sad16x8 = aom_obmc_sad16x8_avx2;
    aom_obmc_sad32x16 = aom_obmc_sad32x16_c;
    if (flags & HAS_AVX2) aom_obmc_sad32x16 = aom_obmc_sad32x16_avx2;
    aom_obmc_sad32x32 = aom_obmc_sad32x32_c;
    if (flags & HAS_AVX2) aom_obmc_sad32x32 = aom_obmc_sad32x32_avx2;
    aom_obmc_sad32x64 = aom_obmc_sad32x64_c;
    if (flags & HAS_AVX2) aom_obmc_sad32x64 = aom_obmc_sad32x64_avx2;
    aom_obmc_sad32x8 = aom_obmc_sad32x8_c;
    if (flags & HAS_AVX2) aom_obmc_sad32x8 = aom_obmc_sad32x8_avx2;
    aom_obmc_sad4x16 = aom_obmc_sad4x16_c;
    if (flags & HAS_AVX2) aom_obmc_sad4x16 = aom_obmc_sad4x16_avx2;
    aom_obmc_sad4x4 = aom_obmc_sad4x4_c;
    if (flags & HAS_AVX2) aom_obmc_sad4x4 = aom_obmc_sad4x4_avx2;
    aom_obmc_sad4x8 = aom_obmc_sad4x8_c;
    if (flags & HAS_AVX2) aom_obmc_sad4x8 = aom_obmc_sad4x8_avx2;
    aom_obmc_sad64x128 = aom_obmc_sad64x128_c;
    if (flags & HAS_AVX2) aom_obmc_sad64x128 = aom_obmc_sad64x128_avx2;
    aom_obmc_sad64x16 = aom_obmc_sad64x16_c;
    if (flags & HAS_AVX2) aom_obmc_sad64x16 = aom_obmc_sad64x16_avx2;
    aom_obmc_sad64x32 = aom_obmc_sad64x32_c;
    if (flags & HAS_AVX2) aom_obmc_sad64x32 = aom_obmc_sad64x32_avx2;
    aom_obmc_sad64x64 = aom_obmc_sad64x64_c;
    if (flags & HAS_AVX2) aom_obmc_sad64x64 = aom_obmc_sad64x64_avx2;
    aom_obmc_sad8x16 = aom_obmc_sad8x16_c;
    if (flags & HAS_AVX2) aom_obmc_sad8x16 = aom_obmc_sad8x16_avx2;
    aom_obmc_sad8x32 = aom_obmc_sad8x32_c;
    if (flags & HAS_AVX2) aom_obmc_sad8x32 = aom_obmc_sad8x32_avx2;
    aom_obmc_sad8x4 = aom_obmc_sad8x4_c;
    if (flags & HAS_AVX2) aom_obmc_sad8x4 = aom_obmc_sad8x4_avx2;
    aom_obmc_sad8x8 = aom_obmc_sad8x8_c;
    if (flags & HAS_AVX2) aom_obmc_sad8x8 = aom_obmc_sad8x8_avx2;
    aom_obmc_sub_pixel_variance128x128 = aom_obmc_sub_pixel_variance128x128_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance128x128 = aom_obmc_sub_pixel_variance128x128_sse4_1;
    aom_obmc_sub_pixel_variance128x64 = aom_obmc_sub_pixel_variance128x64_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance128x64 = aom_obmc_sub_pixel_variance128x64_sse4_1;
    aom_obmc_sub_pixel_variance16x16 = aom_obmc_sub_pixel_variance16x16_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x16 = aom_obmc_sub_pixel_variance16x16_sse4_1;
    aom_obmc_sub_pixel_variance16x32 = aom_obmc_sub_pixel_variance16x32_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x32 = aom_obmc_sub_pixel_variance16x32_sse4_1;
    aom_obmc_sub_pixel_variance16x4 = aom_obmc_sub_pixel_variance16x4_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x4 = aom_obmc_sub_pixel_variance16x4_sse4_1;
    aom_obmc_sub_pixel_variance16x64 = aom_obmc_sub_pixel_variance16x64_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x64 = aom_obmc_sub_pixel_variance16x64_sse4_1;
    aom_obmc_sub_pixel_variance16x8 = aom_obmc_sub_pixel_variance16x8_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x8 = aom_obmc_sub_pixel_variance16x8_sse4_1;
    aom_obmc_sub_pixel_variance32x16 = aom_obmc_sub_pixel_variance32x16_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance32x16 = aom_obmc_sub_pixel_variance32x16_sse4_1;
    aom_obmc_sub_pixel_variance32x32 = aom_obmc_sub_pixel_variance32x32_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance32x32 = aom_obmc_sub_pixel_variance32x32_sse4_1;
    aom_obmc_sub_pixel_variance32x64 = aom_obmc_sub_pixel_variance32x64_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance32x64 = aom_obmc_sub_pixel_variance32x64_sse4_1;
    aom_obmc_sub_pixel_variance32x8 = aom_obmc_sub_pixel_variance32x8_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance32x8 = aom_obmc_sub_pixel_variance32x8_sse4_1;
    aom_obmc_sub_pixel_variance4x16 = aom_obmc_sub_pixel_variance4x16_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance4x16 = aom_obmc_sub_pixel_variance4x16_sse4_1;
    aom_obmc_sub_pixel_variance4x4 = aom_obmc_sub_pixel_variance4x4_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance4x4 = aom_obmc_sub_pixel_variance4x4_sse4_1;
    aom_obmc_sub_pixel_variance4x8 = aom_obmc_sub_pixel_variance4x8_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance4x8 = aom_obmc_sub_pixel_variance4x8_sse4_1;
    aom_obmc_sub_pixel_variance64x128 = aom_obmc_sub_pixel_variance64x128_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance64x128 = aom_obmc_sub_pixel_variance64x128_sse4_1;
    aom_obmc_sub_pixel_variance64x16 = aom_obmc_sub_pixel_variance64x16_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance64x16 = aom_obmc_sub_pixel_variance64x16_sse4_1;
    aom_obmc_sub_pixel_variance64x32 = aom_obmc_sub_pixel_variance64x32_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance64x32 = aom_obmc_sub_pixel_variance64x32_sse4_1;
    aom_obmc_sub_pixel_variance64x64 = aom_obmc_sub_pixel_variance64x64_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance64x64 = aom_obmc_sub_pixel_variance64x64_sse4_1;
    aom_obmc_sub_pixel_variance8x16 = aom_obmc_sub_pixel_variance8x16_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance8x16 = aom_obmc_sub_pixel_variance8x16_sse4_1;
    aom_obmc_sub_pixel_variance8x32 = aom_obmc_sub_pixel_variance8x32_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance8x32 = aom_obmc_sub_pixel_variance8x32_sse4_1;
    aom_obmc_sub_pixel_variance8x4 = aom_obmc_sub_pixel_variance8x4_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance8x4 = aom_obmc_sub_pixel_variance8x4_sse4_1;
    aom_obmc_sub_pixel_variance8x8 = aom_obmc_sub_pixel_variance8x8_c;
    if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance8x8 = aom_obmc_sub_pixel_variance8x8_sse4_1;
    aom_obmc_variance128x128 = aom_obmc_variance128x128_c;
    if (flags & HAS_AVX2) aom_obmc_variance128x128 = aom_obmc_variance128x128_avx2;
    aom_obmc_variance128x64 = aom_obmc_variance128x64_c;
    if (flags & HAS_AVX2) aom_obmc_variance128x64 = aom_obmc_variance128x64_avx2;
    aom_obmc_variance16x16 = aom_obmc_variance16x16_c;
    if (flags & HAS_AVX2) aom_obmc_variance16x16 = aom_obmc_variance16x16_avx2;
    aom_obmc_variance16x32 = aom_obmc_variance16x32_c;
    if (flags & HAS_AVX2) aom_obmc_variance16x32 = aom_obmc_variance16x32_avx2;
    aom_obmc_variance16x4 = aom_obmc_variance16x4_c;
    if (flags & HAS_AVX2) aom_obmc_variance16x4 = aom_obmc_variance16x4_avx2;
    aom_obmc_variance16x64 = aom_obmc_variance16x64_c;
    if (flags & HAS_AVX2) aom_obmc_variance16x64 = aom_obmc_variance16x64_avx2;
    aom_obmc_variance16x8 = aom_obmc_variance16x8_c;
    if (flags & HAS_AVX2) aom_obmc_variance16x8 = aom_obmc_variance16x8_avx2;
    aom_obmc_variance32x16 = aom_obmc_variance32x16_c;
    if (flags & HAS_AVX2) aom_obmc_variance32x16 = aom_obmc_variance32x16_avx2;
    aom_obmc_variance32x32 = aom_obmc_variance32x32_c;
    if (flags & HAS_AVX2) aom_obmc_variance32x32 = aom_obmc_variance32x32_avx2;
    aom_obmc_variance32x64 = aom_obmc_variance32x64_c;
    if (flags & HAS_AVX2) aom_obmc_variance32x64 = aom_obmc_variance32x64_avx2;
    aom_obmc_variance32x8 = aom_obmc_variance32x8_c;
    if (flags & HAS_AVX2) aom_obmc_variance32x8 = aom_obmc_variance32x8_avx2;
    aom_obmc_variance4x16 = aom_obmc_variance4x16_c;
    if (flags & HAS_AVX2) aom_obmc_variance4x16 = aom_obmc_variance4x16_avx2;
    aom_obmc_variance4x4 = aom_obmc_variance4x4_c;
    if (flags & HAS_AVX2) aom_obmc_variance4x4 = aom_obmc_variance4x4_avx2;
    aom_obmc_variance4x8 = aom_obmc_variance4x8_c;
    if (flags & HAS_AVX2) aom_obmc_variance4x8 = aom_obmc_variance4x8_avx2;
    aom_obmc_variance64x128 = aom_obmc_variance64x128_c;
    if (flags & HAS_AVX2) aom_obmc_variance64x128 = aom_obmc_variance64x128_avx2;
    aom_obmc_variance64x16 = aom_obmc_variance64x16_c;
    if (flags & HAS_AVX2) aom_obmc_variance64x16 = aom_obmc_variance64x16_avx2;
    aom_obmc_variance64x32 = aom_obmc_variance64x32_c;
    if (flags & HAS_AVX2) aom_obmc_variance64x32 = aom_obmc_variance64x32_avx2;
    aom_obmc_variance64x64 = aom_obmc_variance64x64_c;
    if (flags & HAS_AVX2) aom_obmc_variance64x64 = aom_obmc_variance64x64_avx2;
    aom_obmc_variance8x16 = aom_obmc_variance8x16_c;
    if (flags & HAS_AVX2) aom_obmc_variance8x16 = aom_obmc_variance8x16_avx2;
    aom_obmc_variance8x32 = aom_obmc_variance8x32_c;
    if (flags & HAS_AVX2) aom_obmc_variance8x32 = aom_obmc_variance8x32_avx2;
    aom_obmc_variance8x4 = aom_obmc_variance8x4_c;
    if (flags & HAS_AVX2) aom_obmc_variance8x4 = aom_obmc_variance8x4_avx2;
    aom_obmc_variance8x8 = aom_obmc_variance8x8_c;
    if (flags & HAS_AVX2) aom_obmc_variance8x8 = aom_obmc_variance8x8_avx2;
#endif
    //VARIANCE
    eb_aom_variance4x4 = eb_aom_variance4x4_c;
    if (flags & HAS_AVX2) eb_aom_variance4x4 = eb_aom_variance4x4_sse2;
    eb_aom_variance4x8 = eb_aom_variance4x8_c;
    if (flags & HAS_AVX2) eb_aom_variance4x8 = eb_aom_variance4x8_sse2;
    eb_aom_variance4x16 = eb_aom_variance4x16_c;
    if (flags & HAS_AVX2) eb_aom_variance4x16 = eb_aom_variance4x16_sse2;
    eb_aom_variance8x4 = eb_aom_variance8x4_c;
    if (flags & HAS_AVX2) eb_aom_variance8x4 = eb_aom_variance8x4_sse2;
    eb_aom_variance8x8 = eb_aom_variance8x8_c;
    if (flags & HAS_AVX2) eb_aom_variance8x8 = eb_aom_variance8x8_sse2;
    eb_aom_variance8x16 = eb_aom_variance8x16_c;
    if (flags & HAS_AVX2) eb_aom_variance8x16 = eb_aom_variance8x16_sse2;
    eb_aom_variance8x32 = eb_aom_variance8x32_c;
    if (flags & HAS_AVX2) eb_aom_variance8x32 = eb_aom_variance8x32_sse2;
    eb_aom_variance16x4 = eb_aom_variance16x4_c;
    if (flags & HAS_AVX2) eb_aom_variance16x4 = eb_aom_variance16x4_avx2;
    eb_aom_variance16x8 = eb_aom_variance16x8_c;
    if (flags & HAS_AVX2) eb_aom_variance16x8 = eb_aom_variance16x8_avx2;
    eb_aom_variance16x16 = eb_aom_variance16x16_c;
    if (flags & HAS_AVX2) eb_aom_variance16x16 = eb_aom_variance16x16_avx2;
    eb_aom_variance16x32 = eb_aom_variance16x32_c;
    if (flags & HAS_AVX2) eb_aom_variance16x32 = eb_aom_variance16x32_avx2;
    eb_aom_variance16x64 = eb_aom_variance16x64_c;
    if (flags & HAS_AVX2) eb_aom_variance16x64 = eb_aom_variance16x64_avx2;
    eb_aom_variance32x8 = eb_aom_variance32x8_c;
    if (flags & HAS_AVX2) eb_aom_variance32x8 = eb_aom_variance32x8_avx2;
    eb_aom_variance32x16 = eb_aom_variance32x16_c;
    if (flags & HAS_AVX2) eb_aom_variance32x16 = eb_aom_variance32x16_avx2;
    eb_aom_variance32x32 = eb_aom_variance32x32_c;
    if (flags & HAS_AVX2) eb_aom_variance32x32 = eb_aom_variance32x32_avx2;
    eb_aom_variance32x64 = eb_aom_variance32x64_c;
    if (flags & HAS_AVX2) eb_aom_variance32x64 = eb_aom_variance32x64_avx2;
    eb_aom_variance64x16 = eb_aom_variance64x16_c;
    if (flags & HAS_AVX2) eb_aom_variance64x16 = eb_aom_variance64x16_avx2;
    eb_aom_variance64x32 = eb_aom_variance64x32_c;
    if (flags & HAS_AVX2) eb_aom_variance64x32 = eb_aom_variance64x32_avx2;
    eb_aom_variance64x64 = eb_aom_variance64x64_c;
    if (flags & HAS_AVX2) eb_aom_variance64x64 = eb_aom_variance64x64_avx2;
    eb_aom_variance64x128 = eb_aom_variance64x128_c;
    if (flags & HAS_AVX2) eb_aom_variance64x128 = eb_aom_variance64x128_avx2;
    eb_aom_variance128x64 = eb_aom_variance128x64_c;
    if (flags & HAS_AVX2) eb_aom_variance128x64 = eb_aom_variance128x64_avx2;
    eb_aom_variance128x128 = eb_aom_variance128x128_c;
    if (flags & HAS_AVX2) eb_aom_variance128x128 = eb_aom_variance128x128_avx2;

    //QIQ
    eb_aom_quantize_b_64x64 = eb_aom_quantize_b_64x64_c_II;
    if (flags & HAS_AVX2) eb_aom_quantize_b_64x64 = eb_aom_quantize_b_64x64_avx2;

    eb_aom_highbd_quantize_b_64x64 = eb_aom_highbd_quantize_b_64x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_quantize_b_64x64 = eb_aom_highbd_quantize_b_64x64_avx2;
    // transform
    eb_av1_fwd_txfm2d_16x8 = eb_av1_fwd_txfm2d_16x8_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x8 = eb_av1_fwd_txfm2d_16x8_avx2;
    eb_av1_fwd_txfm2d_8x16 = eb_av1_fwd_txfm2d_8x16_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x16 = eb_av1_fwd_txfm2d_8x16_avx2;

    eb_av1_fwd_txfm2d_16x4 = eb_av1_fwd_txfm2d_16x4_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x4 = eb_av1_fwd_txfm2d_16x4_avx2;
    eb_av1_fwd_txfm2d_4x16 = eb_av1_fwd_txfm2d_4x16_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_4x16 = eb_av1_fwd_txfm2d_4x16_avx2;

    eb_av1_fwd_txfm2d_8x4 = eb_av1_fwd_txfm2d_8x4_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x4 = eb_av1_fwd_txfm2d_8x4_avx2;
    eb_av1_fwd_txfm2d_4x8 = eb_av1_fwd_txfm2d_4x8_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_4x8 = eb_av1_fwd_txfm2d_4x8_avx2;

    eb_av1_fwd_txfm2d_32x16 = eb_av1_fwd_txfm2d_32x16_c;
    eb_av1_fwd_txfm2d_32x8 = eb_av1_fwd_txfm2d_32x8_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x8 = eb_av1_fwd_txfm2d_32x8_avx2;
    eb_av1_fwd_txfm2d_8x32 = eb_av1_fwd_txfm2d_8x32_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x32 = eb_av1_fwd_txfm2d_8x32_avx2;
    eb_av1_fwd_txfm2d_16x32 = eb_av1_fwd_txfm2d_16x32_c;
    eb_av1_fwd_txfm2d_32x64 = eb_av1_fwd_txfm2d_32x64_c;
    eb_av1_fwd_txfm2d_64x32 = eb_av1_fwd_txfm2d_64x32_c;
    eb_av1_fwd_txfm2d_16x64 = eb_av1_fwd_txfm2d_16x64_c;
    eb_av1_fwd_txfm2d_64x16 = eb_av1_fwd_txfm2d_64x16_c;
    eb_av1_fwd_txfm2d_64x64 = Av1TransformTwoD_64x64_c;
    eb_av1_fwd_txfm2d_32x32 = Av1TransformTwoD_32x32_c;
    eb_av1_fwd_txfm2d_16x16 = Av1TransformTwoD_16x16_c;
#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_av1_fwd_txfm2d_64x64 = av1_fwd_txfm2d_64x64_avx512;
        eb_av1_fwd_txfm2d_32x32 = av1_fwd_txfm2d_32x32_avx512;
        eb_av1_fwd_txfm2d_16x16 = av1_fwd_txfm2d_16x16_avx512;
        eb_av1_fwd_txfm2d_32x64 = av1_fwd_txfm2d_32x64_avx512;
        eb_av1_fwd_txfm2d_64x32 = av1_fwd_txfm2d_64x32_avx512;
        eb_av1_fwd_txfm2d_16x64 = av1_fwd_txfm2d_16x64_avx512;
        eb_av1_fwd_txfm2d_64x16 = av1_fwd_txfm2d_64x16_avx512;
        eb_av1_fwd_txfm2d_32x16 = av1_fwd_txfm2d_32x16_avx512;
        eb_av1_fwd_txfm2d_16x32 = av1_fwd_txfm2d_16x32_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x64 = eb_av1_fwd_txfm2d_64x64_avx2;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x32 = eb_av1_fwd_txfm2d_32x32_avx2;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x16 = eb_av1_fwd_txfm2d_16x16_avx2;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x64 = eb_av1_fwd_txfm2d_32x64_avx2;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x32 = eb_av1_fwd_txfm2d_64x32_avx2;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x64 = eb_av1_fwd_txfm2d_16x64_avx2;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x16 = eb_av1_fwd_txfm2d_64x16_avx2;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x16 = eb_av1_fwd_txfm2d_32x16_avx2;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x32 = eb_av1_fwd_txfm2d_16x32_avx2;
#endif
    eb_av1_fwd_txfm2d_8x8 = Av1TransformTwoD_8x8_c;
    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x8 = eb_av1_fwd_txfm2d_8x8_avx2;
    eb_av1_fwd_txfm2d_4x4 = Av1TransformTwoD_4x4_c;
    if (flags & HAS_SSE4_1) eb_av1_fwd_txfm2d_4x4 = eb_av1_fwd_txfm2d_4x4_sse4_1;

    HandleTransform16x64 = HandleTransform16x64_c;
    if (flags & HAS_AVX2) HandleTransform16x64 = HandleTransform16x64_avx2;
    HandleTransform32x64 = HandleTransform32x64_c;
    if (flags & HAS_AVX2) HandleTransform32x64 = HandleTransform32x64_avx2;
    HandleTransform64x16 = HandleTransform64x16_c;
    if (flags & HAS_AVX2) HandleTransform64x16 = HandleTransform64x16_avx2;
    HandleTransform64x32 = HandleTransform64x32_c;
    if (flags & HAS_AVX2) HandleTransform64x32 = HandleTransform64x32_avx2;
    HandleTransform64x64 = HandleTransform64x64_c;
    if (flags & HAS_AVX2) HandleTransform64x64 = HandleTransform64x64_avx2;

    // eb_aom_highbd_v_predictor
    eb_aom_highbd_v_predictor_16x16 = eb_aom_highbd_v_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x16 = eb_aom_highbd_v_predictor_16x16_avx2;
    eb_aom_highbd_v_predictor_16x32 = eb_aom_highbd_v_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x32 = eb_aom_highbd_v_predictor_16x32_avx2;
    eb_aom_highbd_v_predictor_16x4 = eb_aom_highbd_v_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x4 = eb_aom_highbd_v_predictor_16x4_avx2;
    eb_aom_highbd_v_predictor_16x64 = eb_aom_highbd_v_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x64 = eb_aom_highbd_v_predictor_16x64_avx2;
    eb_aom_highbd_v_predictor_16x8 = eb_aom_highbd_v_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_16x8 = eb_aom_highbd_v_predictor_16x8_avx2;
    eb_aom_highbd_v_predictor_2x2 = eb_aom_highbd_v_predictor_2x2_c;
    eb_aom_highbd_v_predictor_32x16 = eb_aom_highbd_v_predictor_32x16_c;
    eb_aom_highbd_v_predictor_32x32 = eb_aom_highbd_v_predictor_32x32_c;
    eb_aom_highbd_v_predictor_32x64 = eb_aom_highbd_v_predictor_32x64_c;
    eb_aom_highbd_v_predictor_32x8 = eb_aom_highbd_v_predictor_32x8_c;
    eb_aom_highbd_v_predictor_4x16 = eb_aom_highbd_v_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_4x16 = eb_aom_highbd_v_predictor_4x16_sse2;
    eb_aom_highbd_v_predictor_4x4 = eb_aom_highbd_v_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_4x4 = eb_aom_highbd_v_predictor_4x4_sse2;
    eb_aom_highbd_v_predictor_4x8 = eb_aom_highbd_v_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_4x8 = eb_aom_highbd_v_predictor_4x8_sse2;
    eb_aom_highbd_v_predictor_64x16 = eb_aom_highbd_v_predictor_64x16_c;
    eb_aom_highbd_v_predictor_64x32 = eb_aom_highbd_v_predictor_64x32_c;
    eb_aom_highbd_v_predictor_8x32 = eb_aom_highbd_v_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_8x32 = eb_aom_highbd_v_predictor_8x32_sse2;
    eb_aom_highbd_v_predictor_64x64 = eb_aom_highbd_v_predictor_64x64_c;
    eb_aom_highbd_v_predictor_8x16 = eb_aom_highbd_v_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_8x16 = eb_aom_highbd_v_predictor_8x16_sse2;
    eb_aom_highbd_v_predictor_8x4 = eb_aom_highbd_v_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_8x4 = eb_aom_highbd_v_predictor_8x4_sse2;
    eb_aom_highbd_v_predictor_8x8 = eb_aom_highbd_v_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_v_predictor_8x8 = eb_aom_highbd_v_predictor_8x8_sse2;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_highbd_v_predictor_32x8 = aom_highbd_v_predictor_32x8_avx512;
        eb_aom_highbd_v_predictor_32x16 = aom_highbd_v_predictor_32x16_avx512;
        eb_aom_highbd_v_predictor_32x32 = aom_highbd_v_predictor_32x32_avx512;
        eb_aom_highbd_v_predictor_32x64 = aom_highbd_v_predictor_32x64_avx512;
        eb_aom_highbd_v_predictor_64x16 = aom_highbd_v_predictor_64x16_avx512;
        eb_aom_highbd_v_predictor_64x32 = aom_highbd_v_predictor_64x32_avx512;
        eb_aom_highbd_v_predictor_64x64 = aom_highbd_v_predictor_64x64_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_32x8 = eb_aom_highbd_v_predictor_32x8_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_32x16 = eb_aom_highbd_v_predictor_32x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_32x32 = eb_aom_highbd_v_predictor_32x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_32x64 = eb_aom_highbd_v_predictor_32x64_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_64x16 = eb_aom_highbd_v_predictor_64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_64x32 = eb_aom_highbd_v_predictor_64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_v_predictor_64x64 = eb_aom_highbd_v_predictor_64x64_avx2;
#endif // !NON_AVX512_SUPPORT

    //aom_highbd_smooth_predictor
    eb_aom_highbd_smooth_predictor_16x16 = eb_aom_highbd_smooth_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x16 = eb_aom_highbd_smooth_predictor_16x16_avx2;
    eb_aom_highbd_smooth_predictor_16x32 = eb_aom_highbd_smooth_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x32 = eb_aom_highbd_smooth_predictor_16x32_avx2;
    eb_aom_highbd_smooth_predictor_16x4 = eb_aom_highbd_smooth_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x4 = eb_aom_highbd_smooth_predictor_16x4_avx2;
    eb_aom_highbd_smooth_predictor_16x64 = eb_aom_highbd_smooth_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x64 = eb_aom_highbd_smooth_predictor_16x64_avx2;
    eb_aom_highbd_smooth_predictor_16x8 = eb_aom_highbd_smooth_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_16x8 = eb_aom_highbd_smooth_predictor_16x8_avx2;
    eb_aom_highbd_smooth_predictor_2x2 = eb_aom_highbd_smooth_predictor_2x2_c;
    eb_aom_highbd_smooth_predictor_32x16 = eb_aom_highbd_smooth_predictor_32x16_c;
    eb_aom_highbd_smooth_predictor_32x32 = eb_aom_highbd_smooth_predictor_32x32_c;
    eb_aom_highbd_smooth_predictor_32x64 = eb_aom_highbd_smooth_predictor_32x64_c;
    eb_aom_highbd_smooth_predictor_32x8 = eb_aom_highbd_smooth_predictor_32x8_c;
    eb_aom_highbd_smooth_predictor_4x16 = eb_aom_highbd_smooth_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x16 = eb_aom_highbd_smooth_predictor_4x16_ssse3;
    eb_aom_highbd_smooth_predictor_4x4 = eb_aom_highbd_smooth_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x4 = eb_aom_highbd_smooth_predictor_4x4_ssse3;
    eb_aom_highbd_smooth_predictor_4x8 = eb_aom_highbd_smooth_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x8 = eb_aom_highbd_smooth_predictor_4x8_ssse3;
    eb_aom_highbd_smooth_predictor_64x16 = eb_aom_highbd_smooth_predictor_64x16_c;
    eb_aom_highbd_smooth_predictor_64x32 = eb_aom_highbd_smooth_predictor_64x32_c;
    eb_aom_highbd_smooth_predictor_64x64 = eb_aom_highbd_smooth_predictor_64x64_c;
    eb_aom_highbd_smooth_predictor_8x16 = eb_aom_highbd_smooth_predictor_8x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x16 = eb_aom_highbd_smooth_predictor_8x16_avx2;
    eb_aom_highbd_smooth_predictor_8x32 = eb_aom_highbd_smooth_predictor_8x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x32 = eb_aom_highbd_smooth_predictor_8x32_avx2;
    eb_aom_highbd_smooth_predictor_8x4 = eb_aom_highbd_smooth_predictor_8x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x4 = eb_aom_highbd_smooth_predictor_8x4_avx2;
    eb_aom_highbd_smooth_predictor_8x8 = eb_aom_highbd_smooth_predictor_8x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x8 = eb_aom_highbd_smooth_predictor_8x8_avx2;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_highbd_smooth_predictor_32x8 = aom_highbd_smooth_predictor_32x8_avx512;
        eb_aom_highbd_smooth_predictor_32x16 = aom_highbd_smooth_predictor_32x16_avx512;
        eb_aom_highbd_smooth_predictor_32x32 = aom_highbd_smooth_predictor_32x32_avx512;
        eb_aom_highbd_smooth_predictor_32x64 = aom_highbd_smooth_predictor_32x64_avx512;
        eb_aom_highbd_smooth_predictor_64x16 = aom_highbd_smooth_predictor_64x16_avx512;
        eb_aom_highbd_smooth_predictor_64x32 = aom_highbd_smooth_predictor_64x32_avx512;
        eb_aom_highbd_smooth_predictor_64x64 = aom_highbd_smooth_predictor_64x64_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x8 = eb_aom_highbd_smooth_predictor_32x8_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x16 = eb_aom_highbd_smooth_predictor_32x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x32 = eb_aom_highbd_smooth_predictor_32x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x64 = eb_aom_highbd_smooth_predictor_32x64_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x16 = eb_aom_highbd_smooth_predictor_64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x32 = eb_aom_highbd_smooth_predictor_64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x64 = eb_aom_highbd_smooth_predictor_64x64_avx2;
#endif // !NON_AVX512_SUPPORT

    //aom_highbd_smooth_h_predictor
    eb_aom_highbd_smooth_h_predictor_16x16 = eb_aom_highbd_smooth_h_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x16 = eb_aom_highbd_smooth_h_predictor_16x16_avx2;
    eb_aom_highbd_smooth_h_predictor_16x32 = eb_aom_highbd_smooth_h_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x32 = eb_aom_highbd_smooth_h_predictor_16x32_avx2;
    eb_aom_highbd_smooth_h_predictor_16x4 = eb_aom_highbd_smooth_h_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x4 = eb_aom_highbd_smooth_h_predictor_16x4_avx2;
    eb_aom_highbd_smooth_h_predictor_16x64 = eb_aom_highbd_smooth_h_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x64 = eb_aom_highbd_smooth_h_predictor_16x64_avx2;
    eb_aom_highbd_smooth_h_predictor_16x8 = eb_aom_highbd_smooth_h_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_16x8 = eb_aom_highbd_smooth_h_predictor_16x8_avx2;
    eb_aom_highbd_smooth_h_predictor_2x2 = eb_aom_highbd_smooth_h_predictor_2x2_c;
    eb_aom_highbd_smooth_h_predictor_32x16 = eb_aom_highbd_smooth_h_predictor_32x16_c;
    eb_aom_highbd_smooth_h_predictor_32x32 = eb_aom_highbd_smooth_h_predictor_32x32_c;
    eb_aom_highbd_smooth_h_predictor_32x64 = eb_aom_highbd_smooth_h_predictor_32x64_c;
    eb_aom_highbd_smooth_h_predictor_32x8 = eb_aom_highbd_smooth_h_predictor_32x8_c;
    eb_aom_highbd_smooth_h_predictor_4x16 = eb_aom_highbd_smooth_h_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x16 = eb_aom_highbd_smooth_h_predictor_4x16_ssse3;
    eb_aom_highbd_smooth_h_predictor_4x4 = eb_aom_highbd_smooth_h_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x4 = eb_aom_highbd_smooth_h_predictor_4x4_ssse3;
    eb_aom_highbd_smooth_h_predictor_4x8 = eb_aom_highbd_smooth_h_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x8 = eb_aom_highbd_smooth_h_predictor_4x8_ssse3;
    eb_aom_highbd_smooth_h_predictor_64x16 = eb_aom_highbd_smooth_h_predictor_64x16_c;
    eb_aom_highbd_smooth_h_predictor_64x32 = eb_aom_highbd_smooth_h_predictor_64x32_c;
    eb_aom_highbd_smooth_h_predictor_64x64 = eb_aom_highbd_smooth_h_predictor_64x64_c;
    eb_aom_highbd_smooth_h_predictor_8x16 = eb_aom_highbd_smooth_h_predictor_8x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x16 = eb_aom_highbd_smooth_h_predictor_8x16_avx2;
    eb_aom_highbd_smooth_h_predictor_8x32 = eb_aom_highbd_smooth_h_predictor_8x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x32 = eb_aom_highbd_smooth_h_predictor_8x32_avx2;
    eb_aom_highbd_smooth_h_predictor_8x4 = eb_aom_highbd_smooth_h_predictor_8x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x4 = eb_aom_highbd_smooth_h_predictor_8x4_avx2;
    eb_aom_highbd_smooth_h_predictor_8x8 = eb_aom_highbd_smooth_h_predictor_8x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x8 = eb_aom_highbd_smooth_h_predictor_8x8_avx2;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_highbd_smooth_h_predictor_32x8 = aom_highbd_smooth_h_predictor_32x8_avx512;
        eb_aom_highbd_smooth_h_predictor_32x16 = aom_highbd_smooth_h_predictor_32x16_avx512;
        eb_aom_highbd_smooth_h_predictor_32x32 = aom_highbd_smooth_h_predictor_32x32_avx512;
        eb_aom_highbd_smooth_h_predictor_32x64 = aom_highbd_smooth_h_predictor_32x64_avx512;
        eb_aom_highbd_smooth_h_predictor_64x16 = aom_highbd_smooth_h_predictor_64x16_avx512;
        eb_aom_highbd_smooth_h_predictor_64x32 = aom_highbd_smooth_h_predictor_64x32_avx512;
        eb_aom_highbd_smooth_h_predictor_64x64 = aom_highbd_smooth_h_predictor_64x64_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x8 = eb_aom_highbd_smooth_h_predictor_32x8_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x16 = eb_aom_highbd_smooth_h_predictor_32x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x32 = eb_aom_highbd_smooth_h_predictor_32x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x64 = eb_aom_highbd_smooth_h_predictor_32x64_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x16 = eb_aom_highbd_smooth_h_predictor_64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x32 = eb_aom_highbd_smooth_h_predictor_64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x64 = eb_aom_highbd_smooth_h_predictor_64x64_avx2;
#endif

    //aom_highbd_dc_128_predictor
    eb_aom_highbd_dc_128_predictor_16x16 = eb_aom_highbd_dc_128_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x16 = eb_aom_highbd_dc_128_predictor_16x16_avx2;
    eb_aom_highbd_dc_128_predictor_16x32 = eb_aom_highbd_dc_128_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x32 = eb_aom_highbd_dc_128_predictor_16x32_avx2;
    eb_aom_highbd_dc_128_predictor_16x4 = eb_aom_highbd_dc_128_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x4 = eb_aom_highbd_dc_128_predictor_16x4_avx2;
    eb_aom_highbd_dc_128_predictor_16x64 = eb_aom_highbd_dc_128_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x64 = eb_aom_highbd_dc_128_predictor_16x64_avx2;
    eb_aom_highbd_dc_128_predictor_16x8 = eb_aom_highbd_dc_128_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_16x8 = eb_aom_highbd_dc_128_predictor_16x8_avx2;
    eb_aom_highbd_dc_128_predictor_2x2 = eb_aom_highbd_dc_128_predictor_2x2_c;
    eb_aom_highbd_dc_128_predictor_32x16 = eb_aom_highbd_dc_128_predictor_32x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_32x16 = eb_aom_highbd_dc_128_predictor_32x16_avx2;
    eb_aom_highbd_dc_128_predictor_32x32 = eb_aom_highbd_dc_128_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_32x32 = eb_aom_highbd_dc_128_predictor_32x32_avx2;
    eb_aom_highbd_dc_128_predictor_32x64 = eb_aom_highbd_dc_128_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_32x64 = eb_aom_highbd_dc_128_predictor_32x64_avx2;
    eb_aom_highbd_dc_128_predictor_32x8 = eb_aom_highbd_dc_128_predictor_32x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_32x8 = eb_aom_highbd_dc_128_predictor_32x8_avx2;
    eb_aom_highbd_dc_128_predictor_4x16 = eb_aom_highbd_dc_128_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_4x16 = eb_aom_highbd_dc_128_predictor_4x16_sse2;
    eb_aom_highbd_dc_128_predictor_4x4 = eb_aom_highbd_dc_128_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_4x4 = eb_aom_highbd_dc_128_predictor_4x4_sse2;
    eb_aom_highbd_dc_128_predictor_4x8 = eb_aom_highbd_dc_128_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_4x8 = eb_aom_highbd_dc_128_predictor_4x8_sse2;
    eb_aom_highbd_dc_128_predictor_8x32 = eb_aom_highbd_dc_128_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_8x32 = eb_aom_highbd_dc_128_predictor_8x32_sse2;
    eb_aom_highbd_dc_128_predictor_64x16 = eb_aom_highbd_dc_128_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_64x16 = eb_aom_highbd_dc_128_predictor_64x16_avx2;
    eb_aom_highbd_dc_128_predictor_64x32 = eb_aom_highbd_dc_128_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_64x32 = eb_aom_highbd_dc_128_predictor_64x32_avx2;
    eb_aom_highbd_dc_128_predictor_64x64 = eb_aom_highbd_dc_128_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_128_predictor_64x64 = eb_aom_highbd_dc_128_predictor_64x64_avx2;
    eb_aom_highbd_dc_128_predictor_8x16 = eb_aom_highbd_dc_128_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_8x16 = eb_aom_highbd_dc_128_predictor_8x16_sse2;
    eb_aom_highbd_dc_128_predictor_8x4 = eb_aom_highbd_dc_128_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_8x4 = eb_aom_highbd_dc_128_predictor_8x4_sse2;
    eb_aom_highbd_dc_128_predictor_8x8 = eb_aom_highbd_dc_128_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_128_predictor_8x8 = eb_aom_highbd_dc_128_predictor_8x8_sse2;

    //aom_highbd_dc_left_predictor
    eb_aom_highbd_dc_left_predictor_16x16 = eb_aom_highbd_dc_left_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x16 = eb_aom_highbd_dc_left_predictor_16x16_avx2;
    eb_aom_highbd_dc_left_predictor_16x32 = eb_aom_highbd_dc_left_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x32 = eb_aom_highbd_dc_left_predictor_16x32_avx2;
    eb_aom_highbd_dc_left_predictor_16x4 = eb_aom_highbd_dc_left_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x4 = eb_aom_highbd_dc_left_predictor_16x4_avx2;
    eb_aom_highbd_dc_left_predictor_16x64 = eb_aom_highbd_dc_left_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x64 = eb_aom_highbd_dc_left_predictor_16x64_avx2;
    eb_aom_highbd_dc_left_predictor_16x8 = eb_aom_highbd_dc_left_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_16x8 = eb_aom_highbd_dc_left_predictor_16x8_avx2;
    eb_aom_highbd_dc_left_predictor_2x2 = eb_aom_highbd_dc_left_predictor_2x2_c;
    eb_aom_highbd_dc_left_predictor_32x16 = eb_aom_highbd_dc_left_predictor_32x16_c;
    eb_aom_highbd_dc_left_predictor_32x32 = eb_aom_highbd_dc_left_predictor_32x32_c;
    eb_aom_highbd_dc_left_predictor_32x64 = eb_aom_highbd_dc_left_predictor_32x64_c;
    eb_aom_highbd_dc_left_predictor_32x8 = eb_aom_highbd_dc_left_predictor_32x8_c;
    eb_aom_highbd_dc_left_predictor_4x16 = eb_aom_highbd_dc_left_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_4x16 = eb_aom_highbd_dc_left_predictor_4x16_sse2;
    eb_aom_highbd_dc_left_predictor_4x4 = eb_aom_highbd_dc_left_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_4x4 = eb_aom_highbd_dc_left_predictor_4x4_sse2;
    eb_aom_highbd_dc_left_predictor_4x8 = eb_aom_highbd_dc_left_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_4x8 = eb_aom_highbd_dc_left_predictor_4x8_sse2;
    eb_aom_highbd_dc_left_predictor_8x32 = eb_aom_highbd_dc_left_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_8x32 = eb_aom_highbd_dc_left_predictor_8x32_sse2;
    eb_aom_highbd_dc_left_predictor_64x16 = eb_aom_highbd_dc_left_predictor_64x16_c;
    eb_aom_highbd_dc_left_predictor_64x32 = eb_aom_highbd_dc_left_predictor_64x32_c;
    eb_aom_highbd_dc_left_predictor_64x64 = eb_aom_highbd_dc_left_predictor_64x64_c;
    eb_aom_highbd_dc_left_predictor_8x16 = eb_aom_highbd_dc_left_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_8x16 = eb_aom_highbd_dc_left_predictor_8x16_sse2;
    eb_aom_highbd_dc_left_predictor_8x4 = eb_aom_highbd_dc_left_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_8x4 = eb_aom_highbd_dc_left_predictor_8x4_sse2;
    eb_aom_highbd_dc_left_predictor_8x8 = eb_aom_highbd_dc_left_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_left_predictor_8x8 = eb_aom_highbd_dc_left_predictor_8x8_sse2;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_highbd_dc_left_predictor_32x8 = aom_highbd_dc_left_predictor_32x8_avx512;
        eb_aom_highbd_dc_left_predictor_32x16 = aom_highbd_dc_left_predictor_32x16_avx512;
        eb_aom_highbd_dc_left_predictor_32x32 = aom_highbd_dc_left_predictor_32x32_avx512;
        eb_aom_highbd_dc_left_predictor_32x64 = aom_highbd_dc_left_predictor_32x64_avx512;
        eb_aom_highbd_dc_left_predictor_64x16 = aom_highbd_dc_left_predictor_64x16_avx512;
        eb_aom_highbd_dc_left_predictor_64x32 = aom_highbd_dc_left_predictor_64x32_avx512;
        eb_aom_highbd_dc_left_predictor_64x64 = aom_highbd_dc_left_predictor_64x64_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x8 = eb_aom_highbd_dc_left_predictor_32x8_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x16 = eb_aom_highbd_dc_left_predictor_32x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x32 = eb_aom_highbd_dc_left_predictor_32x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_32x64 = eb_aom_highbd_dc_left_predictor_32x64_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x16 = eb_aom_highbd_dc_left_predictor_64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x32 = eb_aom_highbd_dc_left_predictor_64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_left_predictor_64x64 = eb_aom_highbd_dc_left_predictor_64x64_avx2;
#endif // !NON_AVX512_SUPPORT

    eb_aom_highbd_dc_predictor_16x16 = eb_aom_highbd_dc_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x16 = eb_aom_highbd_dc_predictor_16x16_avx2;
    eb_aom_highbd_dc_predictor_16x32 = eb_aom_highbd_dc_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x32 = eb_aom_highbd_dc_predictor_16x32_avx2;
    eb_aom_highbd_dc_predictor_16x4 = eb_aom_highbd_dc_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x4 = eb_aom_highbd_dc_predictor_16x4_avx2;
    eb_aom_highbd_dc_predictor_16x64 = eb_aom_highbd_dc_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x64 = eb_aom_highbd_dc_predictor_16x64_avx2;
    eb_aom_highbd_dc_predictor_16x8 = eb_aom_highbd_dc_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_16x8 = eb_aom_highbd_dc_predictor_16x8_avx2;
    eb_aom_highbd_dc_predictor_2x2 = eb_aom_highbd_dc_predictor_2x2_c;
    eb_aom_highbd_dc_predictor_32x16 = eb_aom_highbd_dc_predictor_32x16_c;
    eb_aom_highbd_dc_predictor_32x32 = eb_aom_highbd_dc_predictor_32x32_c;
    eb_aom_highbd_dc_predictor_32x64 = eb_aom_highbd_dc_predictor_32x64_c;
    eb_aom_highbd_dc_predictor_32x8 = eb_aom_highbd_dc_predictor_32x8_c;
    eb_aom_highbd_dc_predictor_4x16 = eb_aom_highbd_dc_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_4x16 = eb_aom_highbd_dc_predictor_4x16_sse2;
    eb_aom_highbd_dc_predictor_4x4 = eb_aom_highbd_dc_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_4x4 = eb_aom_highbd_dc_predictor_4x4_sse2;
    eb_aom_highbd_dc_predictor_4x8 = eb_aom_highbd_dc_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_4x8 = eb_aom_highbd_dc_predictor_4x8_sse2;
    eb_aom_highbd_dc_predictor_64x16 = eb_aom_highbd_dc_predictor_64x16_c;
    eb_aom_highbd_dc_predictor_64x32 = eb_aom_highbd_dc_predictor_64x32_c;
    eb_aom_highbd_dc_predictor_64x64 = eb_aom_highbd_dc_predictor_64x64_c;
    eb_aom_highbd_dc_predictor_8x16 = eb_aom_highbd_dc_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_8x16 = eb_aom_highbd_dc_predictor_8x16_sse2;
    eb_aom_highbd_dc_predictor_8x4 = eb_aom_highbd_dc_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_8x4 = eb_aom_highbd_dc_predictor_8x4_sse2;
    eb_aom_highbd_dc_predictor_8x8 = eb_aom_highbd_dc_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_8x8 = eb_aom_highbd_dc_predictor_8x8_sse2;
    eb_aom_highbd_dc_predictor_8x32 = eb_aom_highbd_dc_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_predictor_8x32 = eb_aom_highbd_dc_predictor_8x32_sse2;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_highbd_dc_predictor_32x8 = aom_highbd_dc_predictor_32x8_avx512;
        eb_aom_highbd_dc_predictor_32x16 = aom_highbd_dc_predictor_32x16_avx512;
        eb_aom_highbd_dc_predictor_32x32 = aom_highbd_dc_predictor_32x32_avx512;
        eb_aom_highbd_dc_predictor_32x64 = aom_highbd_dc_predictor_32x64_avx512;
        eb_aom_highbd_dc_predictor_64x16 = aom_highbd_dc_predictor_64x16_avx512;
        eb_aom_highbd_dc_predictor_64x32 = aom_highbd_dc_predictor_64x32_avx512;
        eb_aom_highbd_dc_predictor_64x64 = aom_highbd_dc_predictor_64x64_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x8 = eb_aom_highbd_dc_predictor_32x8_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x16 = eb_aom_highbd_dc_predictor_32x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x32 = eb_aom_highbd_dc_predictor_32x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_32x64 = eb_aom_highbd_dc_predictor_32x64_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x16 = eb_aom_highbd_dc_predictor_64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x32 = eb_aom_highbd_dc_predictor_64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_predictor_64x64 = eb_aom_highbd_dc_predictor_64x64_avx2;
#endif // !NON_AVX512_SUPPORT
    //aom_highbd_dc_top_predictor
    eb_aom_highbd_dc_top_predictor_16x16 = eb_aom_highbd_dc_top_predictor_16x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x16 = eb_aom_highbd_dc_top_predictor_16x16_avx2;
    eb_aom_highbd_dc_top_predictor_16x32 = eb_aom_highbd_dc_top_predictor_16x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x32 = eb_aom_highbd_dc_top_predictor_16x32_avx2;
    eb_aom_highbd_dc_top_predictor_16x4 = eb_aom_highbd_dc_top_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x4 = eb_aom_highbd_dc_top_predictor_16x4_avx2;
    eb_aom_highbd_dc_top_predictor_16x64 = eb_aom_highbd_dc_top_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x64 = eb_aom_highbd_dc_top_predictor_16x64_avx2;
    eb_aom_highbd_dc_top_predictor_16x8 = eb_aom_highbd_dc_top_predictor_16x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_16x8 = eb_aom_highbd_dc_top_predictor_16x8_avx2;
    eb_aom_highbd_dc_top_predictor_2x2 = eb_aom_highbd_dc_top_predictor_2x2_c;
    eb_aom_highbd_dc_top_predictor_32x16 = eb_aom_highbd_dc_top_predictor_32x16_c;
    eb_aom_highbd_dc_top_predictor_32x32 = eb_aom_highbd_dc_top_predictor_32x32_c;
    eb_aom_highbd_dc_top_predictor_32x64 = eb_aom_highbd_dc_top_predictor_32x64_c;
    eb_aom_highbd_dc_top_predictor_32x8 = eb_aom_highbd_dc_top_predictor_32x8_c;
    eb_aom_highbd_dc_top_predictor_4x16 = eb_aom_highbd_dc_top_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_4x16 = eb_aom_highbd_dc_top_predictor_4x16_sse2;
    eb_aom_highbd_dc_top_predictor_4x4 = eb_aom_highbd_dc_top_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_4x4 = eb_aom_highbd_dc_top_predictor_4x4_sse2;
    eb_aom_highbd_dc_top_predictor_4x8 = eb_aom_highbd_dc_top_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_4x8 = eb_aom_highbd_dc_top_predictor_4x8_sse2;
    eb_aom_highbd_dc_top_predictor_64x16 = eb_aom_highbd_dc_top_predictor_64x16_c;
    eb_aom_highbd_dc_top_predictor_64x32 = eb_aom_highbd_dc_top_predictor_64x32_c;
    eb_aom_highbd_dc_top_predictor_64x64 = eb_aom_highbd_dc_top_predictor_64x64_c;
    eb_aom_highbd_dc_top_predictor_8x16 = eb_aom_highbd_dc_top_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_8x16 = eb_aom_highbd_dc_top_predictor_8x16_sse2;
    /*if (flags & HAS_SSE2) */eb_aom_highbd_dc_top_predictor_8x32 = eb_aom_highbd_dc_top_predictor_8x32_c;
    eb_aom_highbd_dc_top_predictor_8x4 = eb_aom_highbd_dc_top_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_8x4 = eb_aom_highbd_dc_top_predictor_8x4_sse2;
    eb_aom_highbd_dc_top_predictor_8x8 = eb_aom_highbd_dc_top_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_dc_top_predictor_8x8 = eb_aom_highbd_dc_top_predictor_8x8_sse2;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_highbd_dc_top_predictor_32x8 = aom_highbd_dc_top_predictor_32x8_avx512;
        eb_aom_highbd_dc_top_predictor_32x16 = aom_highbd_dc_top_predictor_32x16_avx512;
        eb_aom_highbd_dc_top_predictor_32x32 = aom_highbd_dc_top_predictor_32x32_avx512;
        eb_aom_highbd_dc_top_predictor_32x64 = aom_highbd_dc_top_predictor_32x64_avx512;
        eb_aom_highbd_dc_top_predictor_64x16 = aom_highbd_dc_top_predictor_64x16_avx512;
        eb_aom_highbd_dc_top_predictor_64x32 = aom_highbd_dc_top_predictor_64x32_avx512;
        eb_aom_highbd_dc_top_predictor_64x64 = aom_highbd_dc_top_predictor_64x64_avx512;
    }
#else
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x8 = eb_aom_highbd_dc_top_predictor_32x8_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x16 = eb_aom_highbd_dc_top_predictor_32x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x32 = eb_aom_highbd_dc_top_predictor_32x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_32x64 = eb_aom_highbd_dc_top_predictor_32x64_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x16 = eb_aom_highbd_dc_top_predictor_64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x32 = eb_aom_highbd_dc_top_predictor_64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_dc_top_predictor_64x64 = eb_aom_highbd_dc_top_predictor_64x64_avx2;
#endif
    // eb_aom_highbd_h_predictor
    eb_aom_highbd_h_predictor_16x4 = eb_aom_highbd_h_predictor_16x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_16x4 = eb_aom_highbd_h_predictor_16x4_avx2;
    eb_aom_highbd_h_predictor_16x64 = eb_aom_highbd_h_predictor_16x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_16x64 = eb_aom_highbd_h_predictor_16x64_avx2;
    eb_aom_highbd_h_predictor_16x8 = eb_aom_highbd_h_predictor_16x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_16x8 = eb_aom_highbd_h_predictor_16x8_sse2;
    eb_aom_highbd_h_predictor_2x2 = eb_aom_highbd_h_predictor_2x2_c;
    eb_aom_highbd_h_predictor_32x16 = eb_aom_highbd_h_predictor_32x16_c;
    eb_aom_highbd_h_predictor_32x32 = eb_aom_highbd_h_predictor_32x32_c;
    eb_aom_highbd_h_predictor_32x64 = eb_aom_highbd_h_predictor_32x64_c;
    eb_aom_highbd_h_predictor_32x8 = eb_aom_highbd_h_predictor_32x8_c;
    eb_aom_highbd_h_predictor_4x16 = eb_aom_highbd_h_predictor_4x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_4x16 = eb_aom_highbd_h_predictor_4x16_sse2;
    eb_aom_highbd_h_predictor_4x4 = eb_aom_highbd_h_predictor_4x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_4x4 = eb_aom_highbd_h_predictor_4x4_sse2;
    eb_aom_highbd_h_predictor_4x8 = eb_aom_highbd_h_predictor_4x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_4x8 = eb_aom_highbd_h_predictor_4x8_sse2;
    eb_aom_highbd_h_predictor_64x16 = eb_aom_highbd_h_predictor_64x16_c;
    eb_aom_highbd_h_predictor_64x32 = eb_aom_highbd_h_predictor_64x32_c;
    eb_aom_highbd_h_predictor_8x32 = eb_aom_highbd_h_predictor_8x32_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_8x32 = eb_aom_highbd_h_predictor_8x32_sse2;
    eb_aom_highbd_h_predictor_64x64 = eb_aom_highbd_h_predictor_64x64_c;
    eb_aom_highbd_h_predictor_8x16 = eb_aom_highbd_h_predictor_8x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_8x16 = eb_aom_highbd_h_predictor_8x16_sse2;
    eb_aom_highbd_h_predictor_8x4 = eb_aom_highbd_h_predictor_8x4_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_8x4 = eb_aom_highbd_h_predictor_8x4_sse2;
    eb_aom_highbd_h_predictor_8x8 = eb_aom_highbd_h_predictor_8x8_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_8x8 = eb_aom_highbd_h_predictor_8x8_sse2;
    eb_aom_highbd_h_predictor_16x16 = eb_aom_highbd_h_predictor_16x16_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_16x16 = eb_aom_highbd_h_predictor_16x16_sse2;
    eb_aom_highbd_h_predictor_16x32 = eb_aom_highbd_h_predictor_16x32_c;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_16x32 = eb_aom_highbd_h_predictor_16x32_sse2;

#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_aom_highbd_h_predictor_32x16 = aom_highbd_h_predictor_32x16_avx512;
        eb_aom_highbd_h_predictor_32x32 = aom_highbd_h_predictor_32x32_avx512;
        eb_aom_highbd_h_predictor_32x64 = aom_highbd_h_predictor_32x64_avx512;
        eb_aom_highbd_h_predictor_32x8 = aom_highbd_h_predictor_32x8_avx512;
        eb_aom_highbd_h_predictor_64x16 = aom_highbd_h_predictor_64x16_avx512;
        eb_aom_highbd_h_predictor_64x32 = aom_highbd_h_predictor_64x32_avx512;
        eb_aom_highbd_h_predictor_64x64 = aom_highbd_h_predictor_64x64_avx512;
    }
#else
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_32x16 = eb_aom_highbd_h_predictor_32x16_sse2;
    if (flags & HAS_SSE2) eb_aom_highbd_h_predictor_32x32 = eb_aom_highbd_h_predictor_32x32_sse2;
    if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_32x64 = eb_aom_highbd_h_predictor_32x64_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_32x8 = eb_aom_highbd_h_predictor_32x8_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_64x16 = eb_aom_highbd_h_predictor_64x16_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_64x32 = eb_aom_highbd_h_predictor_64x32_avx2;
    if (flags & HAS_AVX2) eb_aom_highbd_h_predictor_64x64 = eb_aom_highbd_h_predictor_64x64_avx2;
#endif

    eb_aom_fft2x2_float = eb_aom_fft2x2_float_c;
    eb_aom_fft4x4_float = eb_aom_fft4x4_float_c;
    if (flags & HAS_SSE2) eb_aom_fft4x4_float = eb_aom_fft4x4_float_sse2;
    eb_aom_fft16x16_float = eb_aom_fft16x16_float_c;
    if (flags & HAS_AVX2) eb_aom_fft16x16_float = eb_aom_fft16x16_float_avx2;
    eb_aom_fft32x32_float = eb_aom_fft32x32_float_c;
    if (flags & HAS_AVX2) eb_aom_fft32x32_float = eb_aom_fft32x32_float_avx2;
    eb_aom_fft8x8_float = eb_aom_fft8x8_float_c;
    if (flags & HAS_AVX2) eb_aom_fft8x8_float = eb_aom_fft8x8_float_avx2;


    eb_aom_ifft16x16_float = eb_aom_ifft16x16_float_c;
    if (flags & HAS_AVX2) eb_aom_ifft16x16_float = eb_aom_ifft16x16_float_avx2;
    eb_aom_ifft32x32_float = eb_aom_ifft32x32_float_c;
    if (flags & HAS_AVX2) eb_aom_ifft32x32_float = eb_aom_ifft32x32_float_avx2;
    eb_aom_ifft8x8_float = eb_aom_ifft8x8_float_c;
    if (flags & HAS_AVX2) eb_aom_ifft8x8_float = eb_aom_ifft8x8_float_avx2;
    eb_aom_ifft2x2_float = eb_aom_ifft2x2_float_c;
    eb_aom_ifft4x4_float = eb_aom_ifft4x4_float_c;
    if (flags & HAS_SSE2) eb_aom_ifft4x4_float = eb_aom_ifft4x4_float_sse2;
    av1_get_gradient_hist = av1_get_gradient_hist_c;
    if (flags & HAS_AVX2) av1_get_gradient_hist = av1_get_gradient_hist_avx2;

    SET_AVX2_AVX512(eb_av1_convolve_2d_sr,
                    eb_av1_convolve_2d_sr_c,
                    eb_av1_convolve_2d_sr_avx2,
                    eb_av1_convolve_2d_sr_avx512);
    SET_AVX2_AVX512(eb_av1_convolve_2d_copy_sr,
                    eb_av1_convolve_2d_copy_sr_c,
                    eb_av1_convolve_2d_copy_sr_avx2,
                    eb_av1_convolve_2d_copy_sr_avx512);
    SET_AVX2_AVX512(eb_av1_convolve_x_sr,
                    eb_av1_convolve_x_sr_c,
                    eb_av1_convolve_x_sr_avx2,
                    eb_av1_convolve_x_sr_avx512);
    SET_AVX2_AVX512(eb_av1_convolve_y_sr,
                    eb_av1_convolve_y_sr_c,
                    eb_av1_convolve_y_sr_avx2,
                    eb_av1_convolve_y_sr_avx512);
    SET_AVX2_AVX512(eb_av1_jnt_convolve_2d,
                    eb_av1_jnt_convolve_2d_c,
                    eb_av1_jnt_convolve_2d_avx2,
                    eb_av1_jnt_convolve_2d_avx512);
    SET_AVX2_AVX512(eb_av1_jnt_convolve_2d_copy,
                    eb_av1_jnt_convolve_2d_copy_c,
                    eb_av1_jnt_convolve_2d_copy_avx2,
                    eb_av1_jnt_convolve_2d_copy_avx512);
    SET_AVX2_AVX512(eb_av1_jnt_convolve_x,
                    eb_av1_jnt_convolve_x_c,
                    eb_av1_jnt_convolve_x_avx2,
                    eb_av1_jnt_convolve_x_avx512);
    SET_AVX2_AVX512(eb_av1_jnt_convolve_y,
                    eb_av1_jnt_convolve_y_c,
                    eb_av1_jnt_convolve_y_avx2,
                    eb_av1_jnt_convolve_y_avx512);
    SET_AVX2_AVX512(search_one_dual,
                    search_one_dual_c,
                    search_one_dual_avx2,
                    search_one_dual_avx512);
    SET_AVX2_AVX512(spatial_full_distortion_kernel,
                    spatial_full_distortion_kernel_c,
                    spatial_full_distortion_kernel_avx2,
                    spatial_full_distortion_kernel_avx512);
    SET_AVX2(full_distortion_kernel16_bits,
             full_distortion_kernel16_bits_c,
             full_distortion_kernel16_bits_avx2);
    SET_AVX2_AVX512(residual_kernel8bit,
                    residual_kernel8bit_c,
                    residual_kernel8bit_avx2,
                    residual_kernel8bit_avx512);
    SET_SSE41_AVX2(sad_loop_kernel_sparse,
                   sad_loop_kernel_sparse_c,
                   sad_loop_kernel_sparse_sse4_1_intrin,
                   sad_loop_kernel_sparse_avx2_intrin);
    SET_SSE41_AVX2_AVX512(sad_loop_kernel,
                          sad_loop_kernel_c,
                          sad_loop_kernel_sse4_1_intrin,
                          sad_loop_kernel_avx2_intrin,
                          sad_loop_kernel_avx512_intrin);
    SET_SSE41_AVX2(sad_loop_kernel_hme_l0,
                   sad_loop_kernel_c,
                   sad_loop_kernel_sse4_1_hme_l0_intrin,
                   sad_loop_kernel_avx2_hme_l0_intrin);
    SET_AVX2(noise_extract_luma_weak,
             noise_extract_luma_weak_c,
             noise_extract_luma_weak_avx2_intrin);
    SET_AVX2(noise_extract_luma_weak_lcu,
             noise_extract_luma_weak_lcu_c,
             noise_extract_luma_weak_lcu_avx2_intrin);
    SET_AVX2(noise_extract_luma_strong,
             noise_extract_luma_strong_c,
             noise_extract_luma_strong_avx2_intrin);
    SET_AVX2(noise_extract_chroma_strong,
             noise_extract_chroma_strong_c,
             noise_extract_chroma_strong_avx2_intrin);
    SET_AVX2(noise_extract_chroma_weak,
             noise_extract_chroma_weak_c,
             noise_extract_chroma_weak_avx2_intrin);
    SET_SSE41(svt_av1_apply_filtering,
              svt_av1_apply_filtering_c,
              svt_av1_apply_temporal_filter_sse4_1);
    SET_SSE41(svt_av1_apply_filtering_highbd,
              svt_av1_apply_filtering_highbd_c,
              svt_av1_highbd_apply_temporal_filter_sse4_1);
    SET_AVX2_AVX512(combined_averaging_ssd,
                    combined_averaging_ssd_c,
                    combined_averaging_ssd_avx2,
                    combined_averaging_ssd_avx512);
    SET_AVX2(ext_sad_calculation_8x8_16x16,
             ext_sad_calculation_8x8_16x16_c,
             ext_sad_calculation_8x8_16x16_avx2_intrin);
    SET_SSE41(ext_sad_calculation_32x32_64x64,
              ext_sad_calculation_32x32_64x64_c,
              ext_sad_calculation_32x32_64x64_sse4_intrin);
    SET_SSE2(sad_calculation_8x8_16x16,
             sad_calculation_8x8_16x16_c,
             sad_calculation_8x8_16x16_sse2_intrin);
    SET_SSE2(sad_calculation_32x32_64x64,
             sad_calculation_32x32_64x64_c,
             sad_calculation_32x32_64x64_sse2_intrin);
    SET_AVX2(ext_all_sad_calculation_8x8_16x16,
             ext_all_sad_calculation_8x8_16x16_c,
             ext_all_sad_calculation_8x8_16x16_avx2);
    SET_AVX2(ext_eigth_sad_calculation_nsq,
             ext_eigth_sad_calculation_nsq_c,
             ext_eigth_sad_calculation_nsq_avx2);
    SET_AVX2(ext_eight_sad_calculation_32x32_64x64,
             ext_eight_sad_calculation_32x32_64x64_c,
             ext_eight_sad_calculation_32x32_64x64_avx2);
    SET_AVX2(eb_sad_kernel4x4,
             fast_loop_nxm_sad_kernel,
             eb_compute4x_m_sad_avx2_intrin);
    SET_AVX2(sum_residual8bit,
             sum_residual_c,
             sum_residual8bit_avx2_intrin);
    SET_AVX2(full_distortion_kernel_cbf_zero32_bits,
             full_distortion_kernel_cbf_zero32_bits_c,
             full_distortion_kernel_cbf_zero32_bits_avx2);
    SET_AVX2(full_distortion_kernel32_bits,
             full_distortion_kernel32_bits_c,
             full_distortion_kernel32_bits_avx2);
    SET_AVX2(compressed_packmsb,
             compressed_packmsb_c,
             compressed_packmsb_avx2_intrin);
    SET_AVX2(c_pack,
             c_pack_c,
             c_pack_avx2_intrin);
    SET_SSE2_AVX2(unpack_avg,
                  unpack_avg_c,
                  unpack_avg_sse2_intrin,
                  unpack_avg_avx2_intrin);
    SET_AVX2(unpack_avg_safe_sub,
             unpack_avg_safe_sub_c,
             unpack_avg_safe_sub_avx2_intrin);
    SET_AVX2(un_pack8_bit_data,
             un_pack8_bit_data_c,
             eb_enc_un_pack8_bit_data_avx2_intrin);
    SET_SSE2(picture_average_kernel,
             picture_average_kernel_c,
             picture_average_kernel_sse2_intrin);
    SET_SSE2(picture_average_kernel1_line,
             picture_average_kernel1_line_c,
             picture_average_kernel1_line_sse2_intrin);
    SET_SSE41_AVX2(get_eight_horizontal_search_point_results_8x8_16x16_pu,
                   get_eight_horizontal_search_point_results_8x8_16x16_pu_c,
                   get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin,
                   get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin);
    SET_SSE41_AVX2(get_eight_horizontal_search_point_results_32x32_64x64_pu,
                   get_eight_horizontal_search_point_results_32x32_64x64_pu_c,
                   get_eight_horizontal_search_point_results_32x32_64x64_pu_sse41_intrin,
                   get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin);
    SET_SSE2(initialize_buffer_32bits,
             initialize_buffer_32bits_c,
             initialize_buffer_32bits_sse2_intrin);
    SET_SSE41(compute8x8_satd_u8,
              compute8x8_satd_u8_c,
              compute8x8_satd_u8_sse4);
    SET_AVX2(nxm_sad_kernel_sub_sampled,
             nxm_sad_kernel_helper_c,
             nxm_sad_kernel_sub_sampled_helper_avx2);
    SET_AVX2(nxm_sad_kernel,
             nxm_sad_kernel_helper_c,
             nxm_sad_kernel_helper_avx2);
    SET_AVX2(nxm_sad_avg_kernel,
             nxm_sad_avg_kernel_helper_c,
             nxm_sad_avg_kernel_helper_avx2);
    SET_SSSE3(avc_style_luma_interpolation_filter,
              avc_style_luma_interpolation_filter_helper_c,
              avc_style_luma_interpolation_filter_helper_ssse3);
    SET_SSE2_AVX2(compute_mean_8x8,
                  compute_mean_c,
                  compute_mean8x8_sse2_intrin,
                  compute_mean8x8_avx2_intrin);
    SET_SSE2(compute_mean_square_values_8x8,
             compute_mean_squared_values_c,
             compute_mean_of_squared_values8x8_sse2_intrin);
    SET_SSE2_AVX2(pack2d_16_bit_src_mul4,
                  eb_enc_msb_pack2_d,
                  eb_enc_msb_pack2d_sse2_intrin,
                  eb_enc_msb_pack2d_avx2_intrin_al);
    SET_SSE2(un_pack2d_16_bit_src_mul4,
             eb_enc_msb_un_pack2_d,
             eb_enc_msb_un_pack2d_sse2_intrin);
    SET_SSE2_AVX2(compute_interm_var_four8x8,
             compute_interm_var_four8x8_c,
             compute_interm_var_four8x8_helper_sse2,
             compute_interm_var_four8x8_avx2_intrin);
    SET_AVX2(sad_16b_kernel,
             sad_16b_kernel_c,
             sad_16bit_kernel_avx2);
    SET_SSE41(eb_av1_highbd_warp_affine,
        eb_av1_highbd_warp_affine_c,
        eb_av1_highbd_warp_affine_sse4_1);
    SET_AVX2(av1_compute_cross_correlation,
             av1_compute_cross_correlation_c,
             av1_compute_cross_correlation_avx2);

#if AUTO_MAX_PARTITION
    av1_nn_predict = av1_nn_predict_c;
    if (flags & HAS_SSE3) av1_nn_predict = av1_nn_predict_sse3;
#endif
}


