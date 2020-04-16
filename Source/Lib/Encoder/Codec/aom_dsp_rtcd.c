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
#include "EbPictureAnalysisProcess.h"
#include "EbTemporalFiltering.h"
#include "EbComputeSAD.h"
#include "EbMotionEstimation.h"
#include "EbPictureOperators.h"
#include "EbComputeMean.h"
#include "EbMeSadCalculation.h"
#include "EbAvcStyleMcp.h"


/**************************************
 * Instruction Set Support
 **************************************/

#ifndef NON_AVX512_SUPPORT
#define SET_FUNCTIONS(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512) \
    do {                                                                                      \
        ptr = c;                                                                              \
        if (((uintptr_t)NULL != (uintptr_t)mmx) && (flags & HAS_MMX)) ptr = mmx;              \
        if (((uintptr_t)NULL != (uintptr_t)sse) && (flags & HAS_SSE)) ptr = sse;              \
        if (((uintptr_t)NULL != (uintptr_t)sse2) && (flags & HAS_SSE2)) ptr = sse2;           \
        if (((uintptr_t)NULL != (uintptr_t)sse3) && (flags & HAS_SSE3)) ptr = sse3;           \
        if (((uintptr_t)NULL != (uintptr_t)ssse3) && (flags & HAS_SSSE3)) ptr = ssse3;        \
        if (((uintptr_t)NULL != (uintptr_t)sse4_1) && (flags & HAS_SSE4_1)) ptr = sse4_1;     \
        if (((uintptr_t)NULL != (uintptr_t)sse4_2) && (flags & HAS_SSE4_2)) ptr = sse4_2;     \
        if (((uintptr_t)NULL != (uintptr_t)avx) && (flags & HAS_AVX)) ptr = avx;              \
        if (((uintptr_t)NULL != (uintptr_t)avx2) && (flags & HAS_AVX2)) ptr = avx2;           \
        if (((uintptr_t)NULL != (uintptr_t)avx512) && (flags & HAS_AVX512F)) ptr = avx512;    \
    } while (0)
#else
#define SET_FUNCTIONS(ptr, c, mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, avx512) \
    do {                                                                                      \
        ptr = c;                                                                              \
        if (((uintptr_t)NULL != (uintptr_t)mmx) && (flags & HAS_MMX)) ptr = mmx;              \
        if (((uintptr_t)NULL != (uintptr_t)sse) && (flags & HAS_SSE)) ptr = sse;              \
        if (((uintptr_t)NULL != (uintptr_t)sse2) && (flags & HAS_SSE2)) ptr = sse2;           \
        if (((uintptr_t)NULL != (uintptr_t)sse3) && (flags & HAS_SSE3)) ptr = sse3;           \
        if (((uintptr_t)NULL != (uintptr_t)ssse3) && (flags & HAS_SSSE3)) ptr = ssse3;        \
        if (((uintptr_t)NULL != (uintptr_t)sse4_1) && (flags & HAS_SSE4_1)) ptr = sse4_1;     \
        if (((uintptr_t)NULL != (uintptr_t)sse4_2) && (flags & HAS_SSE4_2)) ptr = sse4_2;     \
        if (((uintptr_t)NULL != (uintptr_t)avx) && (flags & HAS_AVX)) ptr = avx;              \
        if (((uintptr_t)NULL != (uintptr_t)avx2) && (flags & HAS_AVX2)) ptr = avx2;           \
    } while (0)
#endif

void setup_rtcd_internal(CPU_FLAGS flags) {
    /** Should be done during library initialization,
        but for safe limiting cpu flags again. */
    (void)flags;
    //to use C: flags=0
    aom_sse = aom_sse_c;

    aom_highbd_sse = aom_highbd_sse_c;

    av1_wedge_compute_delta_squares = av1_wedge_compute_delta_squares_c;
    av1_wedge_sign_from_residuals = av1_wedge_sign_from_residuals_c;

    eb_compute_cdef_dist = compute_cdef_dist_c;
    eb_compute_cdef_dist_8bit = compute_cdef_dist_8bit_c;

    eb_av1_compute_stats = eb_av1_compute_stats_c;
    eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_c;

    eb_av1_lowbd_pixel_proj_error = eb_av1_lowbd_pixel_proj_error_c;
    eb_av1_highbd_pixel_proj_error = eb_av1_highbd_pixel_proj_error_c;
    eb_av1_calc_frame_error = eb_av1_calc_frame_error_c;

    eb_subtract_average = eb_subtract_average_c;

    get_proj_subspace = get_proj_subspace_c;

    search_one_dual = search_one_dual_c;

    eb_aom_mse16x16 = eb_aom_mse16x16_c;

    eb_aom_quantize_b = eb_aom_quantize_b_c_ii;

    eb_aom_quantize_b_32x32 = eb_aom_quantize_b_32x32_c_ii;

    eb_aom_highbd_quantize_b_32x32 = eb_aom_highbd_quantize_b_32x32_c;

    eb_aom_highbd_quantize_b = eb_aom_highbd_quantize_b_c;

    eb_av1_quantize_fp = eb_av1_quantize_fp_c;

    eb_av1_quantize_fp_32x32 = eb_av1_quantize_fp_32x32_c;

    eb_av1_quantize_fp_64x64 = eb_av1_quantize_fp_64x64_c;

    eb_av1_highbd_quantize_fp = eb_av1_highbd_quantize_fp_c;

    eb_aom_highbd_8_mse16x16 = eb_aom_highbd_8_mse16x16_c;


    //SAD
    eb_aom_sad4x4 = eb_aom_sad4x4_c;
    eb_aom_sad4x4x4d = eb_aom_sad4x4x4d_c;
    eb_aom_sad4x16 = eb_aom_sad4x16_c;
    eb_aom_sad4x16x4d = eb_aom_sad4x16x4d_c;
    eb_aom_sad4x8 = eb_aom_sad4x8_c;
    eb_aom_sad4x8x4d = eb_aom_sad4x8x4d_c;
    eb_aom_sad64x128x4d = eb_aom_sad64x128x4d_c;
    eb_aom_sad64x16x4d = eb_aom_sad64x16x4d_c;
    eb_aom_sad64x32x4d = eb_aom_sad64x32x4d_c;
    eb_aom_sad64x64x4d = eb_aom_sad64x64x4d_c;
    eb_aom_sad8x16 = eb_aom_sad8x16_c;
    eb_aom_sad8x16x4d = eb_aom_sad8x16x4d_c;
    eb_aom_sad8x32 = eb_aom_sad8x32_c;
    eb_aom_sad8x32x4d = eb_aom_sad8x32x4d_c;
    eb_aom_sad8x8 = eb_aom_sad8x8_c;
    eb_aom_sad8x8x4d = eb_aom_sad8x8x4d_c;
    eb_aom_sad16x4 = eb_aom_sad16x4_c;
    eb_aom_sad16x4x4d = eb_aom_sad16x4x4d_c;
    eb_aom_sad32x8 = eb_aom_sad32x8_c;
    eb_aom_sad32x8x4d = eb_aom_sad32x8x4d_c;
    eb_aom_sad16x64 = eb_aom_sad16x64_c;
    eb_aom_sad16x64x4d = eb_aom_sad16x64x4d_c;
    eb_aom_sad32x16 = eb_aom_sad32x16_c;
    eb_aom_sad32x16x4d = eb_aom_sad32x16x4d_c;
    eb_aom_sad16x32 = eb_aom_sad16x32_c;
    eb_aom_sad16x32x4d = eb_aom_sad16x32x4d_c;
    eb_aom_sad32x64 = eb_aom_sad32x64_c;
    eb_aom_sad32x64x4d = eb_aom_sad32x64x4d_c;
    eb_aom_sad32x32 = eb_aom_sad32x32_c;
    eb_aom_sad32x32x4d = eb_aom_sad32x32x4d_c;
    eb_aom_sad16x16 = eb_aom_sad16x16_c;
    eb_aom_sad16x16x4d = eb_aom_sad16x16x4d_c;
    eb_aom_sad16x8 = eb_aom_sad16x8_c;
    eb_aom_sad16x8x4d = eb_aom_sad16x8x4d_c;
    eb_aom_sad8x4 = eb_aom_sad8x4_c;
    eb_aom_sad8x4x4d = eb_aom_sad8x4x4d_c;

    eb_aom_sad64x128 = eb_aom_sad64x128_c;
    eb_aom_sad64x16 = eb_aom_sad64x16_c;
    eb_aom_sad64x32 = eb_aom_sad64x32_c;
    eb_aom_sad64x64 = eb_aom_sad64x64_c;
    eb_aom_sad128x128 = eb_aom_sad128x128_c;
    eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_c;
    eb_aom_sad128x64 = eb_aom_sad128x64_c;
    eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_c;
    eb_av1_txb_init_levels = eb_av1_txb_init_levels_c;

    aom_upsampled_pred = aom_upsampled_pred_c;

    aom_obmc_sad128x128 = aom_obmc_sad128x128_c;
    aom_obmc_sad128x64 = aom_obmc_sad128x64_c;
    aom_obmc_sad16x16 = aom_obmc_sad16x16_c;
    aom_obmc_sad16x32 = aom_obmc_sad16x32_c;
    aom_obmc_sad16x4 = aom_obmc_sad16x4_c;
    aom_obmc_sad16x64 = aom_obmc_sad16x64_c;
    aom_obmc_sad16x8 = aom_obmc_sad16x8_c;
    aom_obmc_sad32x16 = aom_obmc_sad32x16_c;
    aom_obmc_sad32x32 = aom_obmc_sad32x32_c;
    aom_obmc_sad32x64 = aom_obmc_sad32x64_c;
    aom_obmc_sad32x8 = aom_obmc_sad32x8_c;
    aom_obmc_sad4x16 = aom_obmc_sad4x16_c;
    aom_obmc_sad4x4 = aom_obmc_sad4x4_c;
    aom_obmc_sad4x8 = aom_obmc_sad4x8_c;
    aom_obmc_sad64x128 = aom_obmc_sad64x128_c;
    aom_obmc_sad64x16 = aom_obmc_sad64x16_c;
    aom_obmc_sad64x32 = aom_obmc_sad64x32_c;
    aom_obmc_sad64x64 = aom_obmc_sad64x64_c;
    aom_obmc_sad8x16 = aom_obmc_sad8x16_c;
    aom_obmc_sad8x32 = aom_obmc_sad8x32_c;
    aom_obmc_sad8x4 = aom_obmc_sad8x4_c;
    aom_obmc_sad8x8 = aom_obmc_sad8x8_c;
    aom_obmc_sub_pixel_variance128x128 = aom_obmc_sub_pixel_variance128x128_c;
    aom_obmc_sub_pixel_variance128x64 = aom_obmc_sub_pixel_variance128x64_c;
    aom_obmc_sub_pixel_variance16x16 = aom_obmc_sub_pixel_variance16x16_c;
    aom_obmc_sub_pixel_variance16x32 = aom_obmc_sub_pixel_variance16x32_c;
    aom_obmc_sub_pixel_variance16x4 = aom_obmc_sub_pixel_variance16x4_c;
    aom_obmc_sub_pixel_variance16x64 = aom_obmc_sub_pixel_variance16x64_c;
    aom_obmc_sub_pixel_variance16x8 = aom_obmc_sub_pixel_variance16x8_c;
    aom_obmc_sub_pixel_variance32x16 = aom_obmc_sub_pixel_variance32x16_c;
    aom_obmc_sub_pixel_variance32x32 = aom_obmc_sub_pixel_variance32x32_c;
    aom_obmc_sub_pixel_variance32x64 = aom_obmc_sub_pixel_variance32x64_c;
    aom_obmc_sub_pixel_variance32x8 = aom_obmc_sub_pixel_variance32x8_c;
    aom_obmc_sub_pixel_variance4x16 = aom_obmc_sub_pixel_variance4x16_c;
    aom_obmc_sub_pixel_variance4x4 = aom_obmc_sub_pixel_variance4x4_c;
    aom_obmc_sub_pixel_variance4x8 = aom_obmc_sub_pixel_variance4x8_c;
    aom_obmc_sub_pixel_variance64x128 = aom_obmc_sub_pixel_variance64x128_c;
    aom_obmc_sub_pixel_variance64x16 = aom_obmc_sub_pixel_variance64x16_c;
    aom_obmc_sub_pixel_variance64x32 = aom_obmc_sub_pixel_variance64x32_c;
    aom_obmc_sub_pixel_variance64x64 = aom_obmc_sub_pixel_variance64x64_c;
    aom_obmc_sub_pixel_variance8x16 = aom_obmc_sub_pixel_variance8x16_c;
    aom_obmc_sub_pixel_variance8x32 = aom_obmc_sub_pixel_variance8x32_c;
    aom_obmc_sub_pixel_variance8x4 = aom_obmc_sub_pixel_variance8x4_c;
    aom_obmc_sub_pixel_variance8x8 = aom_obmc_sub_pixel_variance8x8_c;
    aom_obmc_variance128x128 = aom_obmc_variance128x128_c;
    aom_obmc_variance128x64 = aom_obmc_variance128x64_c;
    aom_obmc_variance16x16 = aom_obmc_variance16x16_c;
    aom_obmc_variance16x32 = aom_obmc_variance16x32_c;
    aom_obmc_variance16x4 = aom_obmc_variance16x4_c;
    aom_obmc_variance16x64 = aom_obmc_variance16x64_c;
    aom_obmc_variance16x8 = aom_obmc_variance16x8_c;
    aom_obmc_variance32x16 = aom_obmc_variance32x16_c;
    aom_obmc_variance32x32 = aom_obmc_variance32x32_c;
    aom_obmc_variance32x64 = aom_obmc_variance32x64_c;
    aom_obmc_variance32x8 = aom_obmc_variance32x8_c;
    aom_obmc_variance4x16 = aom_obmc_variance4x16_c;
    aom_obmc_variance4x4 = aom_obmc_variance4x4_c;
    aom_obmc_variance4x8 = aom_obmc_variance4x8_c;
    aom_obmc_variance64x128 = aom_obmc_variance64x128_c;
    aom_obmc_variance64x16 = aom_obmc_variance64x16_c;
    aom_obmc_variance64x32 = aom_obmc_variance64x32_c;
    aom_obmc_variance64x64 = aom_obmc_variance64x64_c;
    aom_obmc_variance8x16 = aom_obmc_variance8x16_c;
    aom_obmc_variance8x32 = aom_obmc_variance8x32_c;
    aom_obmc_variance8x4 = aom_obmc_variance8x4_c;
    aom_obmc_variance8x8 = aom_obmc_variance8x8_c;

    //VARIANCE
    eb_aom_variance4x4 = eb_aom_variance4x4_c;
    eb_aom_variance4x8 = eb_aom_variance4x8_c;
    eb_aom_variance4x16 = eb_aom_variance4x16_c;
    eb_aom_variance8x4 = eb_aom_variance8x4_c;
    eb_aom_variance8x8 = eb_aom_variance8x8_c;
    eb_aom_variance8x16 = eb_aom_variance8x16_c;
    eb_aom_variance8x32 = eb_aom_variance8x32_c;
    eb_aom_variance16x4 = eb_aom_variance16x4_c;
    eb_aom_variance16x8 = eb_aom_variance16x8_c;
    eb_aom_variance16x16 = eb_aom_variance16x16_c;
    eb_aom_variance16x32 = eb_aom_variance16x32_c;
    eb_aom_variance16x64 = eb_aom_variance16x64_c;
    eb_aom_variance32x8 = eb_aom_variance32x8_c;
    eb_aom_variance32x16 = eb_aom_variance32x16_c;
    eb_aom_variance32x32 = eb_aom_variance32x32_c;
    eb_aom_variance32x64 = eb_aom_variance32x64_c;
    eb_aom_variance64x16 = eb_aom_variance64x16_c;
    eb_aom_variance64x32 = eb_aom_variance64x32_c;
    eb_aom_variance64x64 = eb_aom_variance64x64_c;
    eb_aom_variance64x128 = eb_aom_variance64x128_c;
    eb_aom_variance128x64 = eb_aom_variance128x64_c;
    eb_aom_variance128x128 = eb_aom_variance128x128_c;

    // VARIANCE HBP
    eb_aom_highbd_10_variance4x4 = eb_aom_highbd_10_variance4x4_c;
    eb_aom_highbd_10_variance4x8 = eb_aom_highbd_10_variance4x8_c;
    eb_aom_highbd_10_variance4x16 = eb_aom_highbd_10_variance4x16_c;
    eb_aom_highbd_10_variance8x4 = eb_aom_highbd_10_variance8x4_c;
    eb_aom_highbd_10_variance8x8 = eb_aom_highbd_10_variance8x8_c;
    eb_aom_highbd_10_variance8x16 = eb_aom_highbd_10_variance8x16_c;
    eb_aom_highbd_10_variance8x32 = eb_aom_highbd_10_variance8x32_c;
    eb_aom_highbd_10_variance16x4 = eb_aom_highbd_10_variance16x4_c;
    eb_aom_highbd_10_variance16x8 = eb_aom_highbd_10_variance16x8_c;
    eb_aom_highbd_10_variance16x16 = eb_aom_highbd_10_variance16x16_c;
    eb_aom_highbd_10_variance16x32 = eb_aom_highbd_10_variance16x32_c;
    eb_aom_highbd_10_variance16x64 = eb_aom_highbd_10_variance16x64_c;
    eb_aom_highbd_10_variance32x8 = eb_aom_highbd_10_variance32x8_c;
    eb_aom_highbd_10_variance32x16 = eb_aom_highbd_10_variance32x16_c;
    eb_aom_highbd_10_variance32x32 = eb_aom_highbd_10_variance32x32_c;
    eb_aom_highbd_10_variance32x64 = eb_aom_highbd_10_variance32x64_c;
    eb_aom_highbd_10_variance64x16 = eb_aom_highbd_10_variance64x16_c;
    eb_aom_highbd_10_variance64x32 = eb_aom_highbd_10_variance64x32_c;
    eb_aom_highbd_10_variance64x64 = eb_aom_highbd_10_variance64x64_c;
    eb_aom_highbd_10_variance64x128 = eb_aom_highbd_10_variance64x128_c;
    eb_aom_highbd_10_variance128x64 = eb_aom_highbd_10_variance128x64_c;
    eb_aom_highbd_10_variance128x128 = eb_aom_highbd_10_variance128x128_c;

    //QIQ
    eb_aom_quantize_b_64x64 = eb_aom_quantize_b_64x64_c_ii;

    eb_aom_highbd_quantize_b_64x64 = eb_aom_highbd_quantize_b_64x64_c;
    // transform
    eb_av1_fwd_txfm2d_16x8 = eb_av1_fwd_txfm2d_16x8_c;
    eb_av1_fwd_txfm2d_8x16 = eb_av1_fwd_txfm2d_8x16_c;

    eb_av1_fwd_txfm2d_16x4 = eb_av1_fwd_txfm2d_16x4_c;
    eb_av1_fwd_txfm2d_4x16 = eb_av1_fwd_txfm2d_4x16_c;

    eb_av1_fwd_txfm2d_8x4 = eb_av1_fwd_txfm2d_8x4_c;
    eb_av1_fwd_txfm2d_4x8 = eb_av1_fwd_txfm2d_4x8_c;

    eb_av1_fwd_txfm2d_32x16 = eb_av1_fwd_txfm2d_32x16_c;
    eb_av1_fwd_txfm2d_32x8 = eb_av1_fwd_txfm2d_32x8_c;
    eb_av1_fwd_txfm2d_8x32 = eb_av1_fwd_txfm2d_8x32_c;
    eb_av1_fwd_txfm2d_16x32 = eb_av1_fwd_txfm2d_16x32_c;
    eb_av1_fwd_txfm2d_32x64 = eb_av1_fwd_txfm2d_32x64_c;
    eb_av1_fwd_txfm2d_64x32 = eb_av1_fwd_txfm2d_64x32_c;
    eb_av1_fwd_txfm2d_16x64 = eb_av1_fwd_txfm2d_16x64_c;
    eb_av1_fwd_txfm2d_64x16 = eb_av1_fwd_txfm2d_64x16_c;
    eb_av1_fwd_txfm2d_64x64 = av1_transform_two_d_64x64_c;
    eb_av1_fwd_txfm2d_32x32 = av1_transform_two_d_32x32_c;
    eb_av1_fwd_txfm2d_16x16 = av1_transform_two_d_16x16_c;

    eb_av1_fwd_txfm2d_8x8 = av1_transform_two_d_8x8_c;
    eb_av1_fwd_txfm2d_4x4 = av1_transform_two_d_4x4_c;

    handle_transform16x64 = handle_transform16x64_c;
    handle_transform32x64 = handle_transform32x64_c;
    handle_transform64x16 = handle_transform64x16_c;
    handle_transform64x32 = handle_transform64x32_c;
    handle_transform64x64 = handle_transform64x64_c;

    eb_aom_fft2x2_float = eb_aom_fft2x2_float_c;
    eb_aom_fft4x4_float = eb_aom_fft4x4_float_c;
    eb_aom_fft16x16_float = eb_aom_fft16x16_float_c;
    eb_aom_fft32x32_float = eb_aom_fft32x32_float_c;
    eb_aom_fft8x8_float = eb_aom_fft8x8_float_c;

    eb_aom_ifft16x16_float = eb_aom_ifft16x16_float_c;
    eb_aom_ifft32x32_float = eb_aom_ifft32x32_float_c;
    eb_aom_ifft8x8_float = eb_aom_ifft8x8_float_c;
    eb_aom_ifft2x2_float = eb_aom_ifft2x2_float_c;
    eb_aom_ifft4x4_float = eb_aom_ifft4x4_float_c;
    av1_get_gradient_hist = av1_get_gradient_hist_c;

    search_one_dual = search_one_dual_c;
    sad_loop_kernel_sparse = sad_loop_kernel_sparse_c;
    sad_loop_kernel = sad_loop_kernel_c;
    sad_loop_kernel_hme_l0 = sad_loop_kernel_c;

    svt_av1_apply_filtering = svt_av1_apply_filtering_c;
#if ENHANCED_TF
    svt_av1_apply_temporal_filter_planewise = svt_av1_apply_temporal_filter_planewise_c;
#endif
    svt_av1_apply_filtering_highbd = svt_av1_apply_filtering_highbd_c;
    combined_averaging_ssd = combined_averaging_ssd_c;
    ext_sad_calculation_8x8_16x16 = ext_sad_calculation_8x8_16x16_c;
    ext_sad_calculation_32x32_64x64 = ext_sad_calculation_32x32_64x64_c;
    sad_calculation_8x8_16x16 = sad_calculation_8x8_16x16_c;
    sad_calculation_32x32_64x64 = sad_calculation_32x32_64x64_c;
    ext_all_sad_calculation_8x8_16x16 = ext_all_sad_calculation_8x8_16x16_c;
    ext_eigth_sad_calculation_nsq = ext_eigth_sad_calculation_nsq_c;
    ext_eight_sad_calculation_32x32_64x64 = ext_eight_sad_calculation_32x32_64x64_c;
    eb_sad_kernel4x4 = fast_loop_nxm_sad_kernel;
    get_eight_horizontal_search_point_results_8x8_16x16_pu = get_eight_horizontal_search_point_results_8x8_16x16_pu_c;
    get_eight_horizontal_search_point_results_32x32_64x64_pu = get_eight_horizontal_search_point_results_32x32_64x64_pu_c;

    initialize_buffer_32bits = initialize_buffer_32bits_c;
    nxm_sad_kernel_sub_sampled = nxm_sad_kernel_helper_c;
    nxm_sad_kernel = nxm_sad_kernel_helper_c;
    nxm_sad_avg_kernel = nxm_sad_avg_kernel_helper_c;

    compute_mean_8x8 = compute_mean_c;
    compute_mean_square_values_8x8 = compute_mean_squared_values_c;
    compute_sub_mean_8x8 = compute_sub_mean_8x8_c;
    compute_interm_var_four8x8 = compute_interm_var_four8x8_c;
    sad_16b_kernel = sad_16b_kernel_c;
    av1_compute_cross_correlation = av1_compute_cross_correlation_c;
    av1_k_means_dim1 = av1_k_means_dim1_c;
    av1_k_means_dim2 = av1_k_means_dim2_c;
    av1_calc_indices_dim1 = av1_calc_indices_dim1_c;
    av1_calc_indices_dim2 = av1_calc_indices_dim2_c;

    eb_av1_get_nz_map_contexts = eb_av1_get_nz_map_contexts_c;

#if RESTRUCTURE_SAD
    pme_sad_loop_kernel = pme_sad_loop_kernel_c;
#endif

#ifdef ARCH_X86
    flags &= get_cpu_flags_to_use();
    if (flags & HAS_AVX2) aom_sse = aom_sse_avx2;
    if (flags & HAS_AVX2) aom_highbd_sse = aom_highbd_sse_avx2;
    if (flags & HAS_AVX2) av1_wedge_compute_delta_squares = av1_wedge_compute_delta_squares_avx2;
    if (flags & HAS_AVX2) av1_wedge_sign_from_residuals = av1_wedge_sign_from_residuals_avx2;
    if (flags & HAS_AVX2) eb_compute_cdef_dist = compute_cdef_dist_avx2;
    if (flags & HAS_AVX2) eb_compute_cdef_dist_8bit = compute_cdef_dist_8bit_avx2;
    if (flags & HAS_AVX2) eb_av1_compute_stats = eb_av1_compute_stats_avx2;
    if (flags & HAS_AVX2) eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_avx2;
#ifndef NON_AVX512_SUPPORT
    if (flags & HAS_AVX512F) {
        eb_av1_compute_stats = eb_av1_compute_stats_avx512;
        eb_av1_compute_stats_highbd = eb_av1_compute_stats_highbd_avx512;
    }
#endif
        if (flags & HAS_AVX2) eb_av1_highbd_pixel_proj_error = eb_av1_highbd_pixel_proj_error_avx2;
        if (flags & HAS_AVX2) eb_av1_calc_frame_error = eb_av1_calc_frame_error_avx2;
        if (flags & HAS_AVX2) eb_subtract_average = eb_subtract_average_avx2;
        if (flags & HAS_AVX2) get_proj_subspace = get_proj_subspace_avx2;
        if (flags & HAS_AVX2) search_one_dual = search_one_dual_avx2;
        if (flags & HAS_AVX2) eb_aom_mse16x16 = eb_aom_mse16x16_avx2;
        if (flags & HAS_AVX2) eb_aom_quantize_b = eb_aom_quantize_b_avx2;
        if (flags & HAS_AVX2) eb_aom_quantize_b_32x32 = eb_aom_quantize_b_32x32_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_quantize_b_32x32 = eb_aom_highbd_quantize_b_32x32_avx2;
        if (flags & HAS_AVX2) eb_aom_highbd_quantize_b = eb_aom_highbd_quantize_b_avx2;
        if (flags & HAS_AVX2) eb_av1_lowbd_pixel_proj_error = eb_av1_lowbd_pixel_proj_error_avx2;
#ifndef NON_AVX512_SUPPORT
        if (flags & HAS_AVX512F) {
            eb_av1_lowbd_pixel_proj_error = eb_av1_lowbd_pixel_proj_error_avx512;
        }
#endif
            if (flags & HAS_AVX2) eb_av1_quantize_fp = eb_av1_quantize_fp_avx2;
            if (flags & HAS_AVX2) eb_av1_quantize_fp_32x32 = eb_av1_quantize_fp_32x32_avx2;
            if (flags & HAS_AVX2) eb_av1_quantize_fp_64x64 = eb_av1_quantize_fp_64x64_avx2;
            if (flags & HAS_AVX2) eb_av1_highbd_quantize_fp = eb_av1_highbd_quantize_fp_avx2;
            if (flags & HAS_SSE2) eb_aom_highbd_8_mse16x16 = eb_aom_highbd_8_mse16x16_sse2;
            if (flags & HAS_AVX2) eb_aom_sad4x4 = eb_aom_sad4x4_avx2;
            if (flags & HAS_AVX2) eb_aom_sad4x4x4d = eb_aom_sad4x4x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad4x16 = eb_aom_sad4x16_avx2;
            if (flags & HAS_AVX2) eb_aom_sad4x16x4d = eb_aom_sad4x16x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad4x8 = eb_aom_sad4x8_avx2;
            if (flags & HAS_AVX2) eb_aom_sad4x8x4d = eb_aom_sad4x8x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad64x128x4d = eb_aom_sad64x128x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad64x16x4d = eb_aom_sad64x16x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad64x32x4d = eb_aom_sad64x32x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad64x64x4d = eb_aom_sad64x64x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad8x16 = eb_aom_sad8x16_avx2;
            if (flags & HAS_AVX2) eb_aom_sad8x16x4d = eb_aom_sad8x16x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad8x32 = eb_aom_sad8x32_avx2;
            if (flags & HAS_AVX2) eb_aom_sad8x32x4d = eb_aom_sad8x32x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad8x8 = eb_aom_sad8x8_avx2;
            if (flags & HAS_AVX2) eb_aom_sad8x8x4d = eb_aom_sad8x8x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x4 = eb_aom_sad16x4_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x4x4d = eb_aom_sad16x4x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad32x8 = eb_aom_sad32x8_avx2;
            if (flags & HAS_AVX2) eb_aom_sad32x8x4d = eb_aom_sad32x8x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x64 = eb_aom_sad16x64_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x64x4d = eb_aom_sad16x64x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad32x16 = eb_aom_sad32x16_avx2;
            if (flags & HAS_AVX2) eb_aom_sad32x16x4d = eb_aom_sad32x16x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x32 = eb_aom_sad16x32_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x32x4d = eb_aom_sad16x32x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad32x64 = eb_aom_sad32x64_avx2;
            if (flags & HAS_AVX2) eb_aom_sad32x64x4d = eb_aom_sad32x64x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad32x32 = eb_aom_sad32x32_avx2;
            if (flags & HAS_AVX2) eb_aom_sad32x32x4d = eb_aom_sad32x32x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x16 = eb_aom_sad16x16_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x16x4d = eb_aom_sad16x16x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x8 = eb_aom_sad16x8_avx2;
            if (flags & HAS_AVX2) eb_aom_sad16x8x4d = eb_aom_sad16x8x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad8x4 = eb_aom_sad8x4_avx2;
            if (flags & HAS_AVX2) eb_aom_sad8x4x4d = eb_aom_sad8x4x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad64x128 = eb_aom_sad64x128_avx2;
            if (flags & HAS_AVX2) eb_aom_sad64x16 = eb_aom_sad64x16_avx2;
            if (flags & HAS_AVX2) eb_aom_sad64x32 = eb_aom_sad64x32_avx2;
            if (flags & HAS_AVX2) eb_aom_sad64x64 = eb_aom_sad64x64_avx2;
            if (flags & HAS_AVX2) eb_aom_sad128x128 = eb_aom_sad128x128_avx2;
            if (flags & HAS_AVX2) eb_aom_sad128x128x4d = eb_aom_sad128x128x4d_avx2;
            if (flags & HAS_AVX2) eb_aom_sad128x64 = eb_aom_sad128x64_avx2;
            if (flags & HAS_AVX2) eb_aom_sad128x64x4d = eb_aom_sad128x64x4d_avx2;
            if (flags & HAS_AVX2) eb_av1_txb_init_levels = eb_av1_txb_init_levels_avx2;
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
#endif // !NON_AVX512_SUPPORT
                if (flags & HAS_AVX2) aom_upsampled_pred = aom_upsampled_pred_sse2;
                if (flags & HAS_AVX2) aom_obmc_sad128x128 = aom_obmc_sad128x128_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad128x64 = aom_obmc_sad128x64_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad16x16 = aom_obmc_sad16x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad16x32 = aom_obmc_sad16x32_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad16x4 = aom_obmc_sad16x4_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad16x64 = aom_obmc_sad16x64_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad16x8 = aom_obmc_sad16x8_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad32x16 = aom_obmc_sad32x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad32x32 = aom_obmc_sad32x32_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad32x64 = aom_obmc_sad32x64_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad32x8 = aom_obmc_sad32x8_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad4x16 = aom_obmc_sad4x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad4x4 = aom_obmc_sad4x4_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad4x8 = aom_obmc_sad4x8_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad64x128 = aom_obmc_sad64x128_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad64x16 = aom_obmc_sad64x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad64x32 = aom_obmc_sad64x32_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad64x64 = aom_obmc_sad64x64_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad8x16 = aom_obmc_sad8x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad8x32 = aom_obmc_sad8x32_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad8x4 = aom_obmc_sad8x4_avx2;
                if (flags & HAS_AVX2) aom_obmc_sad8x8 = aom_obmc_sad8x8_avx2;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance128x128 = aom_obmc_sub_pixel_variance128x128_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance128x64 = aom_obmc_sub_pixel_variance128x64_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x16 = aom_obmc_sub_pixel_variance16x16_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x32 = aom_obmc_sub_pixel_variance16x32_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x4 = aom_obmc_sub_pixel_variance16x4_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x64 = aom_obmc_sub_pixel_variance16x64_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance16x8 = aom_obmc_sub_pixel_variance16x8_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance32x16 = aom_obmc_sub_pixel_variance32x16_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance32x32 = aom_obmc_sub_pixel_variance32x32_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance32x64 = aom_obmc_sub_pixel_variance32x64_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance32x8 = aom_obmc_sub_pixel_variance32x8_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance4x16 = aom_obmc_sub_pixel_variance4x16_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance4x4 = aom_obmc_sub_pixel_variance4x4_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance4x8 = aom_obmc_sub_pixel_variance4x8_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance64x128 = aom_obmc_sub_pixel_variance64x128_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance64x16 = aom_obmc_sub_pixel_variance64x16_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance64x32 = aom_obmc_sub_pixel_variance64x32_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance64x64 = aom_obmc_sub_pixel_variance64x64_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance8x16 = aom_obmc_sub_pixel_variance8x16_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance8x32 = aom_obmc_sub_pixel_variance8x32_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance8x4 = aom_obmc_sub_pixel_variance8x4_sse4_1;
                if (flags & HAS_SSE4_1) aom_obmc_sub_pixel_variance8x8 = aom_obmc_sub_pixel_variance8x8_sse4_1;
                if (flags & HAS_AVX2) aom_obmc_variance128x128 = aom_obmc_variance128x128_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance128x64 = aom_obmc_variance128x64_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance16x16 = aom_obmc_variance16x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance16x32 = aom_obmc_variance16x32_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance16x4 = aom_obmc_variance16x4_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance16x64 = aom_obmc_variance16x64_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance16x8 = aom_obmc_variance16x8_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance32x16 = aom_obmc_variance32x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance32x32 = aom_obmc_variance32x32_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance32x64 = aom_obmc_variance32x64_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance32x8 = aom_obmc_variance32x8_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance4x16 = aom_obmc_variance4x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance4x4 = aom_obmc_variance4x4_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance4x8 = aom_obmc_variance4x8_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance64x128 = aom_obmc_variance64x128_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance64x16 = aom_obmc_variance64x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance64x32 = aom_obmc_variance64x32_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance64x64 = aom_obmc_variance64x64_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance8x16 = aom_obmc_variance8x16_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance8x32 = aom_obmc_variance8x32_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance8x4 = aom_obmc_variance8x4_avx2;
                if (flags & HAS_AVX2) aom_obmc_variance8x8 = aom_obmc_variance8x8_avx2;
                if (flags & HAS_AVX2) eb_aom_variance4x4 = eb_aom_variance4x4_sse2;
                if (flags & HAS_AVX2) eb_aom_variance4x8 = eb_aom_variance4x8_sse2;
                if (flags & HAS_AVX2) eb_aom_variance4x16 = eb_aom_variance4x16_sse2;
                if (flags & HAS_AVX2) eb_aom_variance8x4 = eb_aom_variance8x4_sse2;
                if (flags & HAS_AVX2) eb_aom_variance8x8 = eb_aom_variance8x8_sse2;
                if (flags & HAS_AVX2) eb_aom_variance8x16 = eb_aom_variance8x16_sse2;
                if (flags & HAS_AVX2) eb_aom_variance8x32 = eb_aom_variance8x32_sse2;
                if (flags & HAS_AVX2) eb_aom_variance16x4 = eb_aom_variance16x4_avx2;
                if (flags & HAS_AVX2) eb_aom_variance16x8 = eb_aom_variance16x8_avx2;
                if (flags & HAS_AVX2) eb_aom_variance16x16 = eb_aom_variance16x16_avx2;
                if (flags & HAS_AVX2) eb_aom_variance16x32 = eb_aom_variance16x32_avx2;
                if (flags & HAS_AVX2) eb_aom_variance16x64 = eb_aom_variance16x64_avx2;
                if (flags & HAS_AVX2) eb_aom_variance32x8 = eb_aom_variance32x8_avx2;
                if (flags & HAS_AVX2) eb_aom_variance32x16 = eb_aom_variance32x16_avx2;
                if (flags & HAS_AVX2) eb_aom_variance32x32 = eb_aom_variance32x32_avx2;
                if (flags & HAS_AVX2) eb_aom_variance32x64 = eb_aom_variance32x64_avx2;
                if (flags & HAS_AVX2) eb_aom_variance64x16 = eb_aom_variance64x16_avx2;
                if (flags & HAS_AVX2) eb_aom_variance64x32 = eb_aom_variance64x32_avx2;
                if (flags & HAS_AVX2) eb_aom_variance64x64 = eb_aom_variance64x64_avx2;
                if (flags & HAS_AVX2) eb_aom_variance64x128 = eb_aom_variance64x128_avx2;
                if (flags & HAS_AVX2) eb_aom_variance128x64 = eb_aom_variance128x64_avx2;
                if (flags & HAS_AVX2) eb_aom_variance128x128 = eb_aom_variance128x128_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance8x8 = eb_aom_highbd_10_variance8x8_sse2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance8x16 = eb_aom_highbd_10_variance8x16_sse2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance8x32 = eb_aom_highbd_10_variance8x32_sse2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance16x4 = eb_aom_highbd_10_variance16x4_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance16x8 = eb_aom_highbd_10_variance16x8_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance16x16 = eb_aom_highbd_10_variance16x16_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance16x32 = eb_aom_highbd_10_variance16x32_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance16x64 = eb_aom_highbd_10_variance16x64_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance32x8 = eb_aom_highbd_10_variance32x8_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance32x16 = eb_aom_highbd_10_variance32x16_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance32x32 = eb_aom_highbd_10_variance32x32_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance32x64 = eb_aom_highbd_10_variance32x64_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance64x16 = eb_aom_highbd_10_variance64x16_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance64x32 = eb_aom_highbd_10_variance64x32_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance64x64 = eb_aom_highbd_10_variance64x64_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance64x128 = eb_aom_highbd_10_variance64x128_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance128x64 = eb_aom_highbd_10_variance128x64_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_10_variance128x128 = eb_aom_highbd_10_variance128x128_avx2;
                if (flags & HAS_AVX2) eb_aom_quantize_b_64x64 = eb_aom_quantize_b_64x64_avx2;
                if (flags & HAS_AVX2) eb_aom_highbd_quantize_b_64x64 = eb_aom_highbd_quantize_b_64x64_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x8 = eb_av1_fwd_txfm2d_16x8_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x16 = eb_av1_fwd_txfm2d_8x16_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x4 = eb_av1_fwd_txfm2d_16x4_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_4x16 = eb_av1_fwd_txfm2d_4x16_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x4 = eb_av1_fwd_txfm2d_8x4_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_4x8 = eb_av1_fwd_txfm2d_4x8_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x8 = eb_av1_fwd_txfm2d_32x8_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x32 = eb_av1_fwd_txfm2d_8x32_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x64 = eb_av1_fwd_txfm2d_64x64_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x32 = eb_av1_fwd_txfm2d_32x32_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x16 = eb_av1_fwd_txfm2d_16x16_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x64 = eb_av1_fwd_txfm2d_32x64_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x32 = eb_av1_fwd_txfm2d_64x32_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x64 = eb_av1_fwd_txfm2d_16x64_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_64x16 = eb_av1_fwd_txfm2d_64x16_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_32x16 = eb_av1_fwd_txfm2d_32x16_avx2;
                if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_16x32 = eb_av1_fwd_txfm2d_16x32_avx2;
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
#endif
                    if (flags & HAS_AVX2) eb_av1_fwd_txfm2d_8x8 = eb_av1_fwd_txfm2d_8x8_avx2;
                    if (flags & HAS_SSE4_1) eb_av1_fwd_txfm2d_4x4 = eb_av1_fwd_txfm2d_4x4_sse4_1;
                    if (flags & HAS_AVX2) handle_transform16x64 = handle_transform16x64_avx2;
                    if (flags & HAS_AVX2) handle_transform32x64 = handle_transform32x64_avx2;
                    if (flags & HAS_AVX2) handle_transform64x16 = handle_transform64x16_avx2;
                    if (flags & HAS_AVX2) handle_transform64x32 = handle_transform64x32_avx2;
                    if (flags & HAS_AVX2) handle_transform64x64 = handle_transform64x64_avx2;
                    if (flags & HAS_SSE2) eb_aom_fft4x4_float = eb_aom_fft4x4_float_sse2;
                    if (flags & HAS_AVX2) eb_aom_fft16x16_float = eb_aom_fft16x16_float_avx2;
                    if (flags & HAS_AVX2) eb_aom_fft32x32_float = eb_aom_fft32x32_float_avx2;
                    if (flags & HAS_AVX2) eb_aom_fft8x8_float = eb_aom_fft8x8_float_avx2;
                    if (flags & HAS_AVX2) eb_aom_ifft16x16_float = eb_aom_ifft16x16_float_avx2;
                    if (flags & HAS_AVX2) eb_aom_ifft32x32_float = eb_aom_ifft32x32_float_avx2;
                    if (flags & HAS_AVX2) eb_aom_ifft8x8_float = eb_aom_ifft8x8_float_avx2;
                    if (flags & HAS_SSE2) eb_aom_ifft4x4_float = eb_aom_ifft4x4_float_sse2;
                    if (flags & HAS_AVX2) av1_get_gradient_hist = av1_get_gradient_hist_avx2;
                    SET_AVX2_AVX512(
                        search_one_dual, search_one_dual_c, search_one_dual_avx2, search_one_dual_avx512);
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

                    SET_SSE41(
                        svt_av1_apply_filtering, svt_av1_apply_filtering_c, svt_av1_apply_temporal_filter_sse4_1);
#if ENHANCED_TF
                    SET_AVX2(svt_av1_apply_temporal_filter_planewise,
                        svt_av1_apply_temporal_filter_planewise_c,
                        svt_av1_apply_temporal_filter_planewise_c);
#endif
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
                    SET_AVX2(eb_sad_kernel4x4, fast_loop_nxm_sad_kernel, eb_compute4x_m_sad_avx2_intrin);
                    SET_SSE41_AVX2_AVX512(get_eight_horizontal_search_point_results_8x8_16x16_pu,
                        get_eight_horizontal_search_point_results_8x8_16x16_pu_c,
                        get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin,
                        get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin,
                        get_eight_horizontal_search_point_results_8x8_16x16_pu_avx512_intrin);
                    SET_SSE41_AVX2(get_eight_horizontal_search_point_results_32x32_64x64_pu,
                        get_eight_horizontal_search_point_results_32x32_64x64_pu_c,
                        get_eight_horizontal_search_point_results_32x32_64x64_pu_sse41_intrin,
                        get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin);
                    SET_SSE2(
                        initialize_buffer_32bits, initialize_buffer_32bits_c, initialize_buffer_32bits_sse2_intrin);
                    SET_AVX2(nxm_sad_kernel_sub_sampled,
                        nxm_sad_kernel_helper_c,
                        nxm_sad_kernel_sub_sampled_helper_avx2);

                    SET_AVX2(nxm_sad_kernel, nxm_sad_kernel_helper_c, nxm_sad_kernel_helper_avx2);
                    SET_AVX2(nxm_sad_avg_kernel, nxm_sad_avg_kernel_helper_c, nxm_sad_avg_kernel_helper_avx2);
                    SET_SSE2_AVX2(
                        compute_mean_8x8, compute_mean_c, compute_mean8x8_sse2_intrin, compute_mean8x8_avx2_intrin);
                    SET_SSE2(compute_mean_square_values_8x8,
                        compute_mean_squared_values_c,
                        compute_mean_of_squared_values8x8_sse2_intrin);
                    SET_SSE2_AVX2(compute_interm_var_four8x8,
                        compute_interm_var_four8x8_c,
                        compute_interm_var_four8x8_helper_sse2,
                        compute_interm_var_four8x8_avx2_intrin);
                    SET_AVX2(sad_16b_kernel, sad_16b_kernel_c, sad_16bit_kernel_avx2);
                    SET_AVX2(av1_compute_cross_correlation,
                        av1_compute_cross_correlation_c,
                        av1_compute_cross_correlation_avx2);
                    SET_AVX2(av1_k_means_dim1, av1_k_means_dim1_c, av1_k_means_dim1_avx2);
                    SET_AVX2(av1_k_means_dim2, av1_k_means_dim2_c, av1_k_means_dim2_avx2);
                    SET_AVX2(av1_calc_indices_dim1, av1_calc_indices_dim1_c, av1_calc_indices_dim1_avx2);
                    SET_AVX2(av1_calc_indices_dim2, av1_calc_indices_dim2_c, av1_calc_indices_dim2_avx2);
                    if (flags & HAS_SSE2) eb_av1_get_nz_map_contexts = eb_av1_get_nz_map_contexts_sse2;

#if RESTRUCTURE_SAD
                    SET_AVX2(pme_sad_loop_kernel, pme_sad_loop_kernel_c, pme_sad_loop_kernel_avx2);
#endif

#endif

}
