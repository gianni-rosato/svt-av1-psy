// clang-format off
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

// This file is generated. Do not edit.
#ifndef COMMON_DSP_RTCD_H_
#define COMMON_DSP_RTCD_H_

//#include "common_dsp_rtcd.h"
#include "EbDefinitions.h"
//#include "EbCodingUnit.h"
//# include "ml.h"


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

#ifdef RTCD_C
#define RTCD_EXTERN                //CHKN RTCD call in effect. declare the function pointers in  encHandle.
#else
#define RTCD_EXTERN extern         //CHKN run time externing the fucntion pointers.
#endif

#ifdef __GNUC__
#define LIKELY(v) __builtin_expect(v, 1)
#define UNLIKELY(v) __builtin_expect(v, 0)
#else
#define LIKELY(v) (v)
#define UNLIKELY(v) (v)
#endif


 /**************************************
 * Instruction Set Support
 **************************************/

#ifdef _WIN32
# include <intrin.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif

#define SET_SSE2(ptr, c, sse2) SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, 0, 0)
#define SET_SSE2_AVX2(ptr, c, sse2, avx2) SET_FUNCTIONS(ptr, c, 0, 0, sse2, 0, 0, 0, 0, 0, avx2, 0)
#define SET_SSSE3(ptr, c, ssse3) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, ssse3, 0, 0, 0, 0, 0)
#define SET_SSE41(ptr, c, sse4_1) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, 0, 0)
#define SET_SSE41(ptr, c, sse4_1) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, 0, 0)
#define SET_SSE41_AVX2(ptr, c, sse4_1, avx2) \
    SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, avx2, 0)
#define SET_SSE41_AVX2_AVX512(ptr, c, sse4_1, avx2, avx512) \
    SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, sse4_1, 0, 0, avx2, avx512)
#define SET_AVX2(ptr, c, avx2) SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, avx2, 0)
#define SET_AVX2_AVX512(ptr, c, avx2, avx512) \
    SET_FUNCTIONS(ptr, c, 0, 0, 0, 0, 0, 0, 0, 0, avx2, avx512)


    // Helper Functions
    CPU_FLAGS get_cpu_flags();
    CPU_FLAGS get_cpu_flags_to_use();
    void setup_common_rtcd_internal(CPU_FLAGS flags);

    void aom_blend_a64_vmask_c(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h);
    void aom_blend_a64_vmask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h);
    RTCD_EXTERN void(*aom_blend_a64_vmask)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h);

    void aom_highbd_blend_a64_vmask_c(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);
    void aom_highbd_blend_a64_vmask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);
    RTCD_EXTERN void(*aom_highbd_blend_a64_vmask)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);

    void aom_highbd_blend_a64_hmask_c(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);
    void aom_highbd_blend_a64_hmask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);
    RTCD_EXTERN void(*aom_highbd_blend_a64_hmask)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);

    void aom_blend_a64_hmask_c(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h);
    void aom_blend_a64_hmask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h);
    RTCD_EXTERN void(*aom_blend_a64_hmask)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h);

    void aom_blend_a64_mask_c(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int w, int h, int subx, int suby);
    void aom_blend_a64_mask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int w, int h, int subx, int suby);
    void aom_blend_a64_mask_avx2(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int w, int h, int subx, int suby);
    RTCD_EXTERN void(*aom_blend_a64_mask)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int w, int h, int subx, int suby);

    void aom_highbd_blend_a64_mask_c(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int w, int h, int subx, int suby, int bd);
    void aom_highbd_blend_a64_mask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int w, int h, int subx, int suby, int bd);
    RTCD_EXTERN void(*aom_highbd_blend_a64_mask)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int w, int h, int subx, int suby, int bd);

    void eb_aom_highbd_blend_a64_vmask_c(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0, uint32_t src0_stride, const uint16_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);
    void eb_aom_highbd_blend_a64_vmask_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0, uint32_t src0_stride, const uint16_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_blend_a64_vmask)(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0, uint32_t src0_stride, const uint16_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);

    void eb_aom_highbd_blend_a64_hmask_c(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0, uint32_t src0_stride, const uint16_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);
    void eb_aom_highbd_blend_a64_hmask_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0, uint32_t src0_stride, const uint16_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);
    RTCD_EXTERN void(*eb_aom_highbd_blend_a64_hmask)(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0, uint32_t src0_stride, const uint16_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h, int bd);

    void eb_cfl_predict_lbd_c(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    void eb_cfl_predict_lbd_avx2(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    RTCD_EXTERN void(*eb_cfl_predict_lbd)(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);

    void eb_cfl_predict_hbd_c(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    void eb_cfl_predict_hbd_avx2(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    RTCD_EXTERN void(*eb_cfl_predict_hbd)(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);

    void eb_av1_filter_intra_predictor_c(uint8_t *dst, ptrdiff_t stride, TxSize tx_size, const uint8_t *above, const uint8_t *left, int32_t mode);
    void eb_av1_filter_intra_predictor_sse4_1(uint8_t *dst, ptrdiff_t stride, TxSize tx_size, const uint8_t *above, const uint8_t *left, int mode);
    RTCD_EXTERN void (*eb_av1_filter_intra_predictor) (uint8_t *dst, ptrdiff_t stride, TxSize tx_size, const uint8_t *above, const uint8_t *left, int32_t mode);

    void eb_av1_filter_intra_edge_high_c_old(uint8_t *p, int32_t sz, int32_t strength);
    void eb_av1_filter_intra_edge_sse4_1(uint8_t *p, int32_t sz, int32_t strength);
    RTCD_EXTERN void(*eb_av1_filter_intra_edge)(uint8_t *p, int32_t sz, int32_t strength);

    void eb_av1_filter_intra_edge_high_c(uint16_t *p, int32_t sz, int32_t strength);
    void eb_av1_filter_intra_edge_high_sse4_1(uint16_t *p, int32_t sz, int32_t strength);
    RTCD_EXTERN void(*eb_av1_filter_intra_edge_high)(uint16_t *p, int32_t sz, int32_t strength);

    void eb_av1_upsample_intra_edge_c(uint8_t *p, int32_t sz);
    void eb_av1_upsample_intra_edge_sse4_1(uint8_t *p, int32_t sz);
    RTCD_EXTERN void(*eb_av1_upsample_intra_edge)(uint8_t *p, int32_t sz);

    // AMIR
    void eb_av1_upsample_intra_edge_high_c(uint16_t *p, int32_t sz, int32_t bd);
    //void eb_av1_upsample_intra_edge_high_sse4_1(uint16_t *p, int32_t sz, int32_t bd);
    //RTCD_EXTERN void(*eb_av1_upsample_intra_edge_high)(uint16_t *p, int32_t sz, int32_t bd);

    void eb_av1_highbd_dr_prediction_z2_c(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    RTCD_EXTERN void(*eb_av1_highbd_dr_prediction_z2)(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);

    void av1_build_compound_diffwtd_mask_d16_c(uint8_t *mask, DIFFWTD_MASK_TYPE mask_type, const CONV_BUF_TYPE *src0, int src0_stride, const CONV_BUF_TYPE *src1, int src1_stride, int h, int w, ConvolveParams *conv_params, int bd);
    void av1_build_compound_diffwtd_mask_d16_avx2(uint8_t *mask, DIFFWTD_MASK_TYPE mask_type, const CONV_BUF_TYPE *src0, int src0_stride, const CONV_BUF_TYPE *src1, int src1_stride, int h, int w, ConvolveParams *conv_params, int bd);
    RTCD_EXTERN void (*av1_build_compound_diffwtd_mask_d16)(uint8_t *mask, DIFFWTD_MASK_TYPE mask_type, const CONV_BUF_TYPE *src0, int src0_stride, const CONV_BUF_TYPE *src1, int src1_stride, int h, int w, ConvolveParams *conv_params, int bd);

/////////------------------------

    void eb_av1_inv_txfm2d_add_4x4_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_4x4_sse4_1(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_4x4_avx2(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_4x4)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);

    void eb_av1_inv_txfm2d_add_8x8_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_8x8_sse4_1(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_8x8_avx2(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_8x8)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x16_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x16_sse4_1(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x16_avx2(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x16_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x16)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);

    void eb_av1_inv_txfm2d_add_32x32_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_32x32_avx2(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_32x32_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_32x32)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);

    void eb_av1_inv_txfm2d_add_64x64_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_64x64_sse4_1(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    void eb_av1_inv_txfm2d_add_64x64_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_64x64)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, int32_t bd);

    void eb_av1_highbd_inv_txfm_add_avx2(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_8x16_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_8x16)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x8_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x8)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_we, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x32_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x32_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x32)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_32x16_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    void eb_av1_inv_txfm2d_add_32x16_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_32x16)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_32x8_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_32x8)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_8x32_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_8x32)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_32x64_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    void eb_av1_inv_txfm2d_add_32x64_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_32x64)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_64x32_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    void eb_av1_inv_txfm2d_add_64x32_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_64x32)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x64_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x64_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x64)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_64x16_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    void eb_av1_inv_txfm2d_add_64x16_avx512(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_64x16)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void eb_av1_inv_txfm2d_add_4x8_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);
    void eb_av1_inv_txfm2d_add_4x8_sse4_1(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_4x8)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);

    void eb_av1_inv_txfm2d_add_8x4_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);
    void eb_av1_inv_txfm2d_add_8x4_sse4_1(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_8x4)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);

    void eb_av1_inv_txfm2d_add_4x16_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);
    void eb_av1_inv_txfm2d_add_4x16_sse4_1(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_4x16)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);

    void eb_av1_inv_txfm2d_add_16x4_c(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);
    void eb_av1_inv_txfm2d_add_16x4_sse4_1(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*eb_av1_inv_txfm2d_add_16x4)(const int32_t *input, uint16_t *output_r, int32_t stride_r, uint16_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size, int32_t bd);

    void eb_av1_inv_txfm_add_c(const TranLow *dqcoeff, uint8_t *dst_r, int32_t stride_r, uint8_t *dst_w, int32_t stride_w, const TxfmParam *txfm_param);
    void eb_av1_inv_txfm_add_ssse3(const TranLow *dqcoeff, uint8_t *dst_r, int32_t stride_r, uint8_t *dst_w, int32_t stride_w, const TxfmParam *txfm_param);
    void eb_av1_inv_txfm_add_avx2(const TranLow *dqcoeff, uint8_t *dst_r, int32_t stride_r, uint8_t *dst_w, int32_t stride_w, const TxfmParam *txfm_param);
    RTCD_EXTERN void(*eb_av1_inv_txfm_add)(const TranLow *dqcoeff, uint8_t *dst_r, int32_t stride_r, uint8_t *dst_w, int32_t stride_w, const TxfmParam *txfm_param);

    RTCD_EXTERN void(*compressed_packmsb)(uint8_t *in8_bit_buffer, uint32_t in8_stride, uint8_t *inn_bit_buffer, uint16_t *out16_bit_buffer, uint32_t inn_stride, uint32_t out_stride, uint32_t width, uint32_t height);
    RTCD_EXTERN void(*c_pack)(const uint8_t *inn_bit_buffer, uint32_t inn_stride, uint8_t *in_compn_bit_buffer, uint32_t out_stride, uint8_t *local_cache, uint32_t width, uint32_t height);
    RTCD_EXTERN void(*unpack_avg)(uint16_t *ref16_l0, uint32_t ref_l0_stride, uint16_t *ref16_l1, uint32_t ref_l1_stride, uint8_t *dst_ptr, uint32_t dst_stride, uint32_t width, uint32_t height);
    RTCD_EXTERN void(*unpack_avg_safe_sub)(uint16_t *ref16_l0, uint32_t ref_l0_stride, uint16_t *ref16_l1, uint32_t ref_l1_stride, uint8_t *dst_ptr, uint32_t dst_stride, EbBool sub_pred, uint32_t width, uint32_t height);
    RTCD_EXTERN void(*un_pack8_bit_data)(uint16_t *in16_bit_buffer, uint32_t in_stride, uint8_t *out8_bit_buffer, uint32_t out8_stride, uint32_t width, uint32_t height);

    RTCD_EXTERN void(*pack2d_16_bit_src_mul4)(uint8_t *in8_bit_buffer, uint32_t in8_stride, uint8_t *inn_bit_buffer, uint16_t *out16_bit_buffer, uint32_t inn_stride, uint32_t out_stride, uint32_t width, uint32_t height);
    RTCD_EXTERN void(*un_pack2d_16_bit_src_mul4)(uint16_t *in16_bit_buffer, uint32_t in_stride, uint8_t *out8_bit_buffer, uint8_t *outn_bit_buffer, uint32_t out8_stride, uint32_t outn_stride, uint32_t width, uint32_t height);

    void residual_kernel8bit_c(uint8_t *input, uint32_t input_stride, uint8_t *pred, uint32_t pred_stride, int16_t *residual, uint32_t residual_stride, uint32_t area_width, uint32_t area_height);
    void residual_kernel8bit_avx2(uint8_t *input, uint32_t input_stride, uint8_t *pred, uint32_t pred_stride, int16_t *residual, uint32_t residual_stride, uint32_t area_width, uint32_t area_height);
    void residual_kernel8bit_avx512(uint8_t *input, uint32_t input_stride, uint8_t *pred, uint32_t pred_stride, int16_t *residual, uint32_t residual_stride, uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void(*residual_kernel8bit)(uint8_t *input, uint32_t input_stride, uint8_t *pred, uint32_t pred_stride, int16_t *residual, uint32_t residual_stride, uint32_t area_width, uint32_t area_height);

    RTCD_EXTERN uint64_t(*compute8x8_satd_u8)(uint8_t *diff, uint64_t *dc_value, uint32_t src_stride);
    RTCD_EXTERN int32_t(*sum_residual8bit)(int16_t *in_ptr, uint32_t size, uint32_t stride_in);
    RTCD_EXTERN void(*full_distortion_kernel_cbf_zero32_bits)(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff, uint32_t recon_coeff_stride, uint64_t distortion_result[DIST_CALC_TOTAL], uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void(*full_distortion_kernel32_bits)(int32_t *coeff, uint32_t coeff_stride, int32_t *recon_coeff, uint32_t recon_coeff_stride, uint64_t distortion_result[DIST_CALC_TOTAL], uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void(*picture_average_kernel)(EbByte src0, uint32_t src0_stride, EbByte src1, uint32_t src1_stride, EbByte dst, uint32_t dst_stride, uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void(*picture_average_kernel1_line)(EbByte src0, EbByte src1, EbByte dst, uint32_t area_width);

    uint64_t spatial_full_distortion_kernel_c(uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);
    uint64_t spatial_full_distortion_kernel_avx2(uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);
    uint64_t spatial_full_distortion_kernel_avx512(uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN uint64_t(*spatial_full_distortion_kernel)(uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);

    uint64_t full_distortion_kernel16_bits_c(uint8_t* input, uint32_t input_offset, uint32_t input_stride, uint8_t* recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);
    uint64_t full_distortion_kernel16_bits_avx2(uint8_t* input, uint32_t input_offset, uint32_t input_stride, uint8_t* recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN uint64_t(*full_distortion_kernel16_bits)(uint8_t* input, uint32_t input_offset, uint32_t input_stride, uint8_t* recon, uint32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);
    RTCD_EXTERN void(*residual_kernel16bit)(uint16_t *input, uint32_t input_stride, uint16_t *pred, uint32_t pred_stride, int16_t *residual, uint32_t residual_stride, uint32_t area_width, uint32_t area_height);

    void avc_style_luma_interpolation_filter_ssse3_helper(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte temp_buf, EbBool skip, uint32_t frac_pos, uint8_t choice);
    RTCD_EXTERN void(*avc_style_luma_interpolation_filter)(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte temp_buf, EbBool skip, uint32_t frac_pos, uint8_t choice);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
// clang-format on
