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

int32_t Check4thGenIntelCoreFeatures()
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

int32_t CanUseIntelAVX512()
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

void setup_rtcd_internal(EbAsm asm_type)
{
    int32_t flags = HAS_MMX | HAS_SSE | HAS_SSE2 | HAS_SSE3 | HAS_SSSE3 | HAS_SSE4_1 | HAS_SSE4_2 | HAS_AVX;

    if (asm_type == ASM_AVX2)
        flags |= HAS_AVX2;
    //if (asm_type == ASM_NON_AVX2)
    //    flags = ~HAS_AVX2;

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
    av1_wedge_sse_from_residuals = av1_wedge_sse_from_residuals_c;
    if (flags & HAS_AVX2) av1_wedge_sse_from_residuals = av1_wedge_sse_from_residuals_avx2;
    aom_subtract_block = aom_subtract_block_c;
    if (flags & HAS_AVX2) aom_subtract_block = aom_subtract_block_avx2;
    aom_sse = aom_sse_c;
    if (flags & HAS_AVX2) aom_sse = aom_sse_avx2;
    av1_build_compound_diffwtd_mask_d16 = av1_build_compound_diffwtd_mask_d16_c;
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
#ifndef NON_AVX512_SUPPORT
    if (CanUseIntelAVX512()) {
        eb_av1_compute_stats = eb_av1_compute_stats_avx512;
        eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_avx512;
        spatial_full_distortion_kernel_func_ptr_array[ASM_AVX2] = spatial_full_distortion_kernel_avx512;
        nxm_sad_loop_kernel_func_ptr_array[ASM_AVX2] = sad_loop_kernel_avx512_intrin;
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

    eb_av1_convolve_2d_sr = eb_av1_convolve_2d_sr_c;
    if (flags & HAS_AVX2) eb_av1_convolve_2d_sr = eb_av1_convolve_2d_sr_avx2;

    eb_av1_jnt_convolve_2d_copy = eb_av1_jnt_convolve_2d_copy_c;
    if (flags & HAS_AVX2) eb_av1_jnt_convolve_2d_copy = eb_av1_jnt_convolve_2d_copy_avx2;

    eb_av1_convolve_x_sr = eb_av1_convolve_x_sr_c;
    if (flags & HAS_AVX2) eb_av1_convolve_x_sr = eb_av1_convolve_x_sr_avx2;
    eb_av1_convolve_y_sr = eb_av1_convolve_y_sr_c;
    if (flags & HAS_AVX2) eb_av1_convolve_y_sr = eb_av1_convolve_y_sr_avx2;

    eb_av1_convolve_2d_scale = eb_av1_convolve_2d_scale_c;
    //if (flags & HAS_SSE4_1) eb_av1_convolve_2d_scale = eb_av1_convolve_2d_scale_sse4_1;

    eb_av1_jnt_convolve_x = eb_av1_jnt_convolve_x_c;
    if (flags & HAS_AVX2) eb_av1_jnt_convolve_x = eb_av1_jnt_convolve_x_avx2;
    eb_av1_jnt_convolve_y = eb_av1_jnt_convolve_y_c;
    if (flags & HAS_AVX2) eb_av1_jnt_convolve_y = eb_av1_jnt_convolve_y_avx2;

    eb_av1_jnt_convolve_2d = eb_av1_jnt_convolve_2d_c;
    if (flags & HAS_AVX2) eb_av1_jnt_convolve_2d = eb_av1_jnt_convolve_2d_avx2;

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
    if (CanUseIntelAVX512()) {
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
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x16 = eb_aom_highbd_smooth_v_predictor_32x16_avx2;
    eb_aom_highbd_smooth_v_predictor_32x32 = eb_aom_highbd_smooth_v_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x32 = eb_aom_highbd_smooth_v_predictor_32x32_avx2;
    eb_aom_highbd_smooth_v_predictor_32x64 = eb_aom_highbd_smooth_v_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x64 = eb_aom_highbd_smooth_v_predictor_32x64_avx2;
    eb_aom_highbd_smooth_v_predictor_32x8 = eb_aom_highbd_smooth_v_predictor_32x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_32x8 = eb_aom_highbd_smooth_v_predictor_32x8_avx2;
    eb_aom_highbd_smooth_v_predictor_4x16 = eb_aom_highbd_smooth_v_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x16 = eb_aom_highbd_smooth_v_predictor_4x16_ssse3;
    eb_aom_highbd_smooth_v_predictor_4x4 = eb_aom_highbd_smooth_v_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x4 = eb_aom_highbd_smooth_v_predictor_4x4_ssse3;
    eb_aom_highbd_smooth_v_predictor_4x8 = eb_aom_highbd_smooth_v_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_v_predictor_4x8 = eb_aom_highbd_smooth_v_predictor_4x8_ssse3;
    eb_aom_highbd_smooth_v_predictor_64x16 = eb_aom_highbd_smooth_v_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x16 = eb_aom_highbd_smooth_v_predictor_64x16_avx2;
    eb_aom_highbd_smooth_v_predictor_64x32 = eb_aom_highbd_smooth_v_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x32 = eb_aom_highbd_smooth_v_predictor_64x32_avx2;
    eb_aom_highbd_smooth_v_predictor_64x64 = eb_aom_highbd_smooth_v_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_64x64 = eb_aom_highbd_smooth_v_predictor_64x64_avx2;
    eb_aom_highbd_smooth_v_predictor_8x16 = eb_aom_highbd_smooth_v_predictor_8x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x16 = eb_aom_highbd_smooth_v_predictor_8x16_avx2;
    eb_aom_highbd_smooth_v_predictor_8x32 = eb_aom_highbd_smooth_v_predictor_8x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x32 = eb_aom_highbd_smooth_v_predictor_8x32_avx2;
    eb_aom_highbd_smooth_v_predictor_8x4 = eb_aom_highbd_smooth_v_predictor_8x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x4 = eb_aom_highbd_smooth_v_predictor_8x4_avx2;
    eb_aom_highbd_smooth_v_predictor_8x8 = eb_aom_highbd_smooth_v_predictor_8x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_v_predictor_8x8 = eb_aom_highbd_smooth_v_predictor_8x8_avx2;

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
    //av1_get_nz_map_contexts = eb_av1_get_nz_map_contexts_c;
    /*if (flags & HAS_SSE2)*/ eb_av1_get_nz_map_contexts = eb_av1_get_nz_map_contexts_sse2;

    ResidualKernel = residual_kernel_c;
    if (flags & HAS_AVX2) ResidualKernel = ResidualKernel_avx2;
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

    eb_av1_txb_init_levels = eb_av1_txb_init_levels_c;
    if (flags & HAS_AVX2) eb_av1_txb_init_levels = eb_av1_txb_init_levels_avx2;
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

#ifndef NON_AVX512_SUPPORT
    eb_aom_sad64x128 = eb_aom_sad64x128_avx512;
    eb_aom_sad64x16 = eb_aom_sad64x16_avx512;
    eb_aom_sad64x32 = eb_aom_sad64x32_avx512;
    eb_aom_sad64x64 = eb_aom_sad64x64_avx512;
    eb_aom_sad128x128 = eb_aom_sad128x128_avx512;
    eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_avx512;
    eb_aom_sad128x64 = eb_aom_sad128x64_avx512;
    eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_avx512;
#else
    eb_aom_sad64x128 = eb_aom_sad64x128_c;
    if (flags & HAS_AVX2) eb_aom_sad64x128 = eb_aom_sad64x128_avx2;
    eb_aom_sad64x16 = eb_aom_sad64x16_c;
    if (flags & HAS_AVX2) eb_aom_sad64x16 = eb_aom_sad64x16_avx2;
    eb_aom_sad64x32 = eb_aom_sad64x32_c;
    if (flags & HAS_AVX2) eb_aom_sad64x32 = eb_aom_sad64x32_avx2;
    eb_aom_sad64x64 = eb_aom_sad64x64_c;
    if (flags & HAS_AVX2) eb_aom_sad64x64 = eb_aom_sad64x64_avx2;
    eb_aom_sad128x128 = eb_aom_sad128x128_c;
    if (flags & HAS_AVX2) eb_aom_sad128x128 = eb_aom_sad128x128_avx2;
    eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_avx2;
    eb_aom_sad128x64 = eb_aom_sad128x64_c;
    if (flags & HAS_AVX2) eb_aom_sad128x64 = eb_aom_sad128x64_avx2;
    eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_c;
    if (flags & HAS_AVX2) eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_avx2;
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
    if (CanUseIntelAVX512()) {
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
    if (CanUseIntelAVX512()) {
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
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x16 = eb_aom_highbd_smooth_predictor_32x16_avx2;
    eb_aom_highbd_smooth_predictor_32x32 = eb_aom_highbd_smooth_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x32 = eb_aom_highbd_smooth_predictor_32x32_avx2;
    eb_aom_highbd_smooth_predictor_32x64 = eb_aom_highbd_smooth_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x64 = eb_aom_highbd_smooth_predictor_32x64_avx2;
    eb_aom_highbd_smooth_predictor_32x8 = eb_aom_highbd_smooth_predictor_32x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_32x8 = eb_aom_highbd_smooth_predictor_32x8_avx2;
    eb_aom_highbd_smooth_predictor_4x16 = eb_aom_highbd_smooth_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x16 = eb_aom_highbd_smooth_predictor_4x16_ssse3;
    eb_aom_highbd_smooth_predictor_4x4 = eb_aom_highbd_smooth_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x4 = eb_aom_highbd_smooth_predictor_4x4_ssse3;
    eb_aom_highbd_smooth_predictor_4x8 = eb_aom_highbd_smooth_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_predictor_4x8 = eb_aom_highbd_smooth_predictor_4x8_ssse3;
    eb_aom_highbd_smooth_predictor_64x16 = eb_aom_highbd_smooth_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x16 = eb_aom_highbd_smooth_predictor_64x16_avx2;
    eb_aom_highbd_smooth_predictor_64x32 = eb_aom_highbd_smooth_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x32 = eb_aom_highbd_smooth_predictor_64x32_avx2;
    eb_aom_highbd_smooth_predictor_64x64 = eb_aom_highbd_smooth_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_64x64 = eb_aom_highbd_smooth_predictor_64x64_avx2;
    eb_aom_highbd_smooth_predictor_8x16 = eb_aom_highbd_smooth_predictor_8x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x16 = eb_aom_highbd_smooth_predictor_8x16_avx2;
    eb_aom_highbd_smooth_predictor_8x32 = eb_aom_highbd_smooth_predictor_8x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x32 = eb_aom_highbd_smooth_predictor_8x32_avx2;
    eb_aom_highbd_smooth_predictor_8x4 = eb_aom_highbd_smooth_predictor_8x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x4 = eb_aom_highbd_smooth_predictor_8x4_avx2;
    eb_aom_highbd_smooth_predictor_8x8 = eb_aom_highbd_smooth_predictor_8x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_predictor_8x8 = eb_aom_highbd_smooth_predictor_8x8_avx2;

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
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x16 = eb_aom_highbd_smooth_h_predictor_32x16_avx2;
    eb_aom_highbd_smooth_h_predictor_32x32 = eb_aom_highbd_smooth_h_predictor_32x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x32 = eb_aom_highbd_smooth_h_predictor_32x32_avx2;
    eb_aom_highbd_smooth_h_predictor_32x64 = eb_aom_highbd_smooth_h_predictor_32x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x64 = eb_aom_highbd_smooth_h_predictor_32x64_avx2;
    eb_aom_highbd_smooth_h_predictor_32x8 = eb_aom_highbd_smooth_h_predictor_32x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_32x8 = eb_aom_highbd_smooth_h_predictor_32x8_avx2;
    eb_aom_highbd_smooth_h_predictor_4x16 = eb_aom_highbd_smooth_h_predictor_4x16_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x16 = eb_aom_highbd_smooth_h_predictor_4x16_ssse3;
    eb_aom_highbd_smooth_h_predictor_4x4 = eb_aom_highbd_smooth_h_predictor_4x4_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x4 = eb_aom_highbd_smooth_h_predictor_4x4_ssse3;
    eb_aom_highbd_smooth_h_predictor_4x8 = eb_aom_highbd_smooth_h_predictor_4x8_c;
    if (flags & HAS_SSSE3) eb_aom_highbd_smooth_h_predictor_4x8 = eb_aom_highbd_smooth_h_predictor_4x8_ssse3;
    eb_aom_highbd_smooth_h_predictor_64x16 = eb_aom_highbd_smooth_h_predictor_64x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x16 = eb_aom_highbd_smooth_h_predictor_64x16_avx2;
    eb_aom_highbd_smooth_h_predictor_64x32 = eb_aom_highbd_smooth_h_predictor_64x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x32 = eb_aom_highbd_smooth_h_predictor_64x32_avx2;
    eb_aom_highbd_smooth_h_predictor_64x64 = eb_aom_highbd_smooth_h_predictor_64x64_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_64x64 = eb_aom_highbd_smooth_h_predictor_64x64_avx2;
    eb_aom_highbd_smooth_h_predictor_8x16 = eb_aom_highbd_smooth_h_predictor_8x16_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x16 = eb_aom_highbd_smooth_h_predictor_8x16_avx2;
    eb_aom_highbd_smooth_h_predictor_8x32 = eb_aom_highbd_smooth_h_predictor_8x32_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x32 = eb_aom_highbd_smooth_h_predictor_8x32_avx2;
    eb_aom_highbd_smooth_h_predictor_8x4 = eb_aom_highbd_smooth_h_predictor_8x4_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x4 = eb_aom_highbd_smooth_h_predictor_8x4_avx2;
    eb_aom_highbd_smooth_h_predictor_8x8 = eb_aom_highbd_smooth_h_predictor_8x8_c;
    if (flags & HAS_AVX2) eb_aom_highbd_smooth_h_predictor_8x8 = eb_aom_highbd_smooth_h_predictor_8x8_avx2;

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
    if (CanUseIntelAVX512()) {
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
    if (CanUseIntelAVX512()) {
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
    if (CanUseIntelAVX512()) {
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
    if (CanUseIntelAVX512()) {
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

}

