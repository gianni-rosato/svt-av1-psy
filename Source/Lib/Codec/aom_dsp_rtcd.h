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
#ifndef AOM_DSP_RTCD_H_
#define AOM_DSP_RTCD_H_

#include "EbDefinitions.h"

#ifdef RTCD_C
#define RTCD_EXTERN                //CHKN RTCD call in effect. declare the function pointers in  encHandle.
#else
#define RTCD_EXTERN extern         //CHKN run time externing the fucntion pointers.
#endif

/*
 * DSP
 */


#define HAS_MMX 0x01
#define HAS_SSE 0x02
#define HAS_SSE2 0x04
#define HAS_SSE3 0x08
#define HAS_SSSE3 0x10
#define HAS_SSE4_1 0x20
#define HAS_AVX 0x40
#define HAS_AVX2 0x80
#define HAS_SSE4_2 0x100


#ifdef __cplusplus
extern "C" {
#endif

    //to not include convolve.h, just forward declare what's needed.
    struct ConvolveParams;
    struct InterpFilterParams;

    void apply_selfguided_restoration_c(const uint8_t *dat, int32_t width, int32_t height, int32_t stride, int32_t eps, const int32_t *xqd, uint8_t *dst, int32_t dst_stride, int32_t *tmpbuf, int32_t bit_depth, int32_t highbd);
    void apply_selfguided_restoration_avx2(const uint8_t *dat, int32_t width, int32_t height, int32_t stride, int32_t eps, const int32_t *xqd, uint8_t *dst, int32_t dst_stride, int32_t *tmpbuf, int32_t bit_depth, int32_t highbd);
    RTCD_EXTERN void(*apply_selfguided_restoration)(const uint8_t *dat, int32_t width, int32_t height, int32_t stride, int32_t eps, const int32_t *xqd, uint8_t *dst, int32_t dst_stride, int32_t *tmpbuf, int32_t bit_depth, int32_t highbd);

    void av1_wiener_convolve_add_src_c(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params);
    void av1_wiener_convolve_add_src_avx2(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_wiener_convolve_add_src)(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params);


    void av1_highbd_wiener_convolve_add_src_c(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params, int32_t bps);
    void av1_highbd_wiener_convolve_add_src_avx2(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params, int32_t bps);
    RTCD_EXTERN void(*av1_highbd_wiener_convolve_add_src)(const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4, const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h, const ConvolveParams *conv_params, int32_t bps);


    void av1_selfguided_restoration_c(const uint8_t *dgd8, int32_t width, int32_t height,
        int32_t dgd_stride, int32_t *flt0, int32_t *flt1, int32_t flt_stride,
        int32_t sgr_params_idx, int32_t bit_depth, int32_t highbd);
    void av1_selfguided_restoration_avx2(const uint8_t *dgd8, int32_t width, int32_t height,
        int32_t dgd_stride, int32_t *flt0, int32_t *flt1, int32_t flt_stride,
        int32_t sgr_params_idx, int32_t bit_depth, int32_t highbd);
    RTCD_EXTERN void(*av1_selfguided_restoration)(const uint8_t *dgd8, int32_t width, int32_t height,
        int32_t dgd_stride, int32_t *flt0, int32_t *flt1, int32_t flt_stride,
        int32_t sgr_params_idx, int32_t bit_depth, int32_t highbd);


    int32_t cdef_find_dir_c(const uint16_t *img, int32_t stride, int32_t *var, int32_t coeff_shift);
    int32_t cdef_find_dir_avx2(const uint16_t *img, int32_t stride, int32_t *var, int32_t coeff_shift);
    RTCD_EXTERN int32_t(*cdef_find_dir)(const uint16_t *img, int32_t stride, int32_t *var, int32_t coeff_shift);

    void cdef_filter_block_c(uint8_t *dst8, uint16_t *dst16, int32_t dstride, const uint16_t *in, int32_t pri_strength, int32_t sec_strength, int32_t dir, int32_t pri_damping, int32_t sec_damping, int32_t bsize, int32_t max, int32_t coeff_shift);
    void cdef_filter_block_avx2(uint8_t *dst8, uint16_t *dst16, int32_t dstride, const uint16_t *in, int32_t pri_strength, int32_t sec_strength, int32_t dir, int32_t pri_damping, int32_t sec_damping, int32_t bsize, int32_t max, int32_t coeff_shift);
    RTCD_EXTERN void(*cdef_filter_block)(uint8_t *dst8, uint16_t *dst16, int32_t dstride, const uint16_t *in, int32_t pri_strength, int32_t sec_strength, int32_t dir, int32_t pri_damping, int32_t sec_damping, int32_t bsize, int32_t max, int32_t coeff_shift);

    void copy_rect8_8bit_to_16bit_c(uint16_t *dst, int32_t dstride, const uint8_t *src, int32_t sstride, int32_t v, int32_t h);
    void copy_rect8_8bit_to_16bit_avx2(uint16_t *dst, int32_t dstride, const uint8_t *src, int32_t sstride, int32_t v, int32_t h);
    RTCD_EXTERN void(*copy_rect8_8bit_to_16bit)(uint16_t *dst, int32_t dstride, const uint8_t *src, int32_t sstride, int32_t v, int32_t h);

    void av1_compute_stats_c(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H);
    void av1_compute_stats_avx2(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H);
    RTCD_EXTERN void(*av1_compute_stats)(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H);

    void av1_compute_stats_highbd_c(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H, aom_bit_depth_t bit_depth);
    void av1_compute_stats_highbd_avx2(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H, aom_bit_depth_t bit_depth);
    RTCD_EXTERN void(*av1_compute_stats_highbd)(int32_t wiener_win, const uint8_t *dgd8, const uint8_t *src8, int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride, int64_t *M, int64_t *H, aom_bit_depth_t bit_depth);

    int64_t av1_lowbd_pixel_proj_error_c(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const sgr_params_type *params);
    int64_t av1_lowbd_pixel_proj_error_avx2(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const sgr_params_type *params);
    RTCD_EXTERN int64_t(*av1_lowbd_pixel_proj_error)(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const sgr_params_type *params);

    int64_t av1_highbd_pixel_proj_error_c(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const sgr_params_type *params);
    int64_t av1_highbd_pixel_proj_error_avx2(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const sgr_params_type *params);
    RTCD_EXTERN int64_t(*av1_highbd_pixel_proj_error)(const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride, const uint8_t *dat8, int32_t dat_stride, int32_t *flt0, int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t xq[2], const sgr_params_type *params);

    void av1_highbd_convolve_2d_copy_sr_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void av1_highbd_convolve_2d_copy_sr_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_convolve_2d_copy_sr)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void av1_highbd_jnt_convolve_2d_copy_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void av1_highbd_jnt_convolve_2d_copy_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_jnt_convolve_2d_copy)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void av1_highbd_convolve_y_sr_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void av1_highbd_convolve_y_sr_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_convolve_y_sr)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void av1_highbd_convolve_2d_sr_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void av1_highbd_convolve_2d_sr_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_convolve_2d_sr)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void av1_highbd_jnt_convolve_2d_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void av1_highbd_jnt_convolve_2d_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_jnt_convolve_2d)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void av1_highbd_jnt_convolve_x_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void av1_highbd_jnt_convolve_x_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_jnt_convolve_x)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void av1_highbd_jnt_convolve_y_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void av1_highbd_jnt_convolve_y_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_jnt_convolve_y)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void av1_highbd_convolve_x_sr_c(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    void av1_highbd_convolve_x_sr_avx2(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_convolve_x_sr)(const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h, const InterpFilterParams *filter_params_x, const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    void subtract_average_c(int16_t *pred_buf_q3, int32_t width, int32_t height, int32_t round_offset, int32_t num_pel_log2);
    void subtract_average_avx2(int16_t *pred_buf_q3, int32_t width, int32_t height, int32_t round_offset, int32_t num_pel_log2);
    RTCD_EXTERN void(*subtract_average)(int16_t *pred_buf_q3, int32_t width, int32_t height, int32_t round_offset, int32_t num_pel_log2);

    void cfl_predict_lbd_c(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    void cfl_predict_lbd_avx2(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    RTCD_EXTERN void(*cfl_predict_lbd)(const int16_t *pred_buf_q3, uint8_t *pred, int32_t pred_stride, uint8_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);

    void cfl_predict_hbd_c(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    void cfl_predict_hbd_avx2(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);
    RTCD_EXTERN void(*cfl_predict_hbd)(const int16_t *pred_buf_q3, uint16_t *pred, int32_t pred_stride, uint16_t *dst, int32_t dst_stride, int32_t alpha_q3, int32_t bit_depth, int32_t width, int32_t height);

#if QT_10BIT_SUPPORT
    void av1_filter_intra_edge_high_c_old(uint8_t *p, int32_t sz, int32_t strength);
#else
    void av1_filter_intra_edge_high_c(uint8_t *p, int32_t sz, int32_t strength);
#endif
    void av1_filter_intra_edge_sse4_1(uint8_t *p, int32_t sz, int32_t strength);
    RTCD_EXTERN void(*av1_filter_intra_edge)(uint8_t *p, int32_t sz, int32_t strength);

#if INTRA_10BIT_SUPPORT
    void av1_filter_intra_edge_high_c(uint16_t *p, int32_t sz, int32_t strength);
    void av1_filter_intra_edge_high_sse4_1(uint16_t *p, int32_t sz, int32_t strength);
    RTCD_EXTERN void(*av1_filter_intra_edge_high)(uint16_t *p, int32_t sz, int32_t strength);
#endif


    void av1_fwd_txfm2d_8x16_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_8x16_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_8x16)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_16x8_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_16x8_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_16x8)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_4x16_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_4x16_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_4x16)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_16x4_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_16x4_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_16x4)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_4x8_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_4x8_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_4x8)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_8x4_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_8x4_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_8x4)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_32x16_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_32x16_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_32x16)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_32x8_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_32x8_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_32x8)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_8x32_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_8x32_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_8x32)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_16x32_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_16x32_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_16x32)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_32x64_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_32x64_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_32x64)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_64x32_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_64x32_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_64x32)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_16x64_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_16x64_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_16x64)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void av1_fwd_txfm2d_64x16_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_64x16_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_64x16)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void Av1TransformTwoD_64x64_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_64x64_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_64x64)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void Av1TransformTwoD_32x32_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_32x32_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_32x32)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);


    void Av1TransformTwoD_16x16_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_16x16_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_16x16)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void Av1TransformTwoD_8x8_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_8x8_avx2(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_8x8)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void Av1TransformTwoD_4x4_c(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    void av1_fwd_txfm2d_4x4_sse4_1(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);
    RTCD_EXTERN void(*av1_fwd_txfm2d_4x4)(int16_t *input, int32_t *output, uint32_t inputStride, TxType transform_type, uint8_t  bit_depth);

    void smooth_v_predictor_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);
    void eb_smooth_v_predictor_all_ssse3(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_smooth_v_predictor)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);


    void smooth_h_predictor_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);
    void eb_smooth_h_predictor_all_ssse3(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*eb_smooth_h_predictor)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left);

    void get_proj_subspace_c(const uint8_t *src8, int width, int height, int src_stride, const uint8_t *dat8, int dat_stride, int use_highbitdepth, int32_t *flt0, int flt0_stride, int32_t *flt1, int flt1_stride, int *xq, const sgr_params_type *params);
    void get_proj_subspace_avx2(const uint8_t *src8, int width, int height, int src_stride, const uint8_t *dat8, int dat_stride, int use_highbitdepth, int32_t *flt0, int flt0_stride, int32_t *flt1, int flt1_stride, int *xq, const sgr_params_type *params);
    RTCD_EXTERN void(*get_proj_subspace)(const uint8_t *src8, int width, int height, int src_stride, const uint8_t *dat8, int dat_stride, int use_highbitdepth, int32_t *flt0, int flt0_stride, int32_t *flt1, int flt1_stride, int *xq, const sgr_params_type *params);

    uint64_t mse_4x4_16bit_c(uint16_t *dst, int dstride, uint16_t *src, int sstride);
    uint64_t mse_4x4_16bit_avx2(uint16_t *dst, int dstride, uint16_t *src, int sstride);
    RTCD_EXTERN uint64_t(*mse_4x4_16bit)(uint16_t *dst, int dstride, uint16_t *src, int sstride);

    uint64_t dist_8x8_16bit_c(uint16_t *dst, int dstride, uint16_t *src, int sstride, int coeff_shift);
    uint64_t dist_8x8_16bit_avx2(uint16_t *dst, int dstride, uint16_t *src, int sstride, int coeff_shift);
    RTCD_EXTERN uint64_t(*dist_8x8_16bit)(uint16_t *dst, int dstride, uint16_t *src, int sstride, int coeff_shift);
#if FAST_CDEF
    uint64_t search_one_dual_c(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast, int start_gi, int end_gi);
    uint64_t search_one_dual_avx2(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast, int start_gi, int end_gi);
    RTCD_EXTERN uint64_t(*search_one_dual)(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast, int start_gi, int end_gi);
#else
    uint64_t search_one_dual_c(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast);
    uint64_t search_one_dual_avx2(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast);
    RTCD_EXTERN uint64_t(*search_one_dual)(int *lev0, int *lev1, int nb_strengths, uint64_t(**mse)[64], int sb_count, int fast);
#endif

    uint32_t aom_mse16x16_c(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
    uint32_t aom_mse16x16_avx2(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
    RTCD_EXTERN uint32_t (*aom_mse16x16)(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);

    void av1_convolve_2d_copy_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void av1_convolve_2d_copy_sr_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_convolve_2d_copy_sr)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);



    void av1_convolve_2d_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void av1_convolve_2d_sr_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_convolve_2d_sr)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void av1_jnt_convolve_2d_copy_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void av1_jnt_convolve_2d_copy_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_jnt_convolve_2d_copy)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);


    void av1_convolve_x_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void av1_convolve_x_sr_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_convolve_x_sr)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void av1_convolve_y_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void av1_convolve_y_sr_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_convolve_y_sr)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void av1_jnt_convolve_x_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void av1_jnt_convolve_x_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_jnt_convolve_x)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void av1_jnt_convolve_y_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void av1_jnt_convolve_y_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_jnt_convolve_y)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);

    void av1_jnt_convolve_2d_c(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    void av1_jnt_convolve_2d_avx2(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);
    RTCD_EXTERN void(*av1_jnt_convolve_2d)(const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h, InterpFilterParams *filter_params_x, InterpFilterParams *filter_params_y, const int32_t subpel_x_q4, const int32_t subpel_y_q4, ConvolveParams *conv_params);


    void aom_quantize_b_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void aom_highbd_quantize_b_avx2(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*aom_quantize_b)(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);


    void aom_quantize_b_32x32_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void aom_highbd_quantize_b_32x32_avx2(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*aom_quantize_b_32x32)(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void aom_quantize_b_64x64_c_II(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    void aom_highbd_quantize_b_64x64_avx2(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*aom_quantize_b_64x64)(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

#if QT_10BIT_SUPPORT

    void aom_highbd_quantize_b_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*aom_highbd_quantize_b)(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void aom_highbd_quantize_b_32x32_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*aom_highbd_quantize_b_32x32)(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

    void aom_highbd_quantize_b_64x64_c(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);
    RTCD_EXTERN void(*aom_highbd_quantize_b_64x64)(const tran_low_t *coeff_ptr, intptr_t n_coeffs, int32_t skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan);

#endif


    void av1_inv_txfm2d_add_4x4_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void av1_inv_txfm2d_add_4x4_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void av1_inv_txfm2d_add_4x4_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_4x4)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void av1_inv_txfm2d_add_8x8_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void av1_inv_txfm2d_add_8x8_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void av1_inv_txfm2d_add_8x8_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_8x8)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void av1_inv_txfm2d_add_16x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void av1_inv_txfm2d_add_16x16_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void av1_inv_txfm2d_add_16x16_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_16x16)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void av1_inv_txfm2d_add_32x32_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void av1_inv_txfm2d_add_32x32_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_32x32)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);

    void av1_inv_txfm2d_add_64x64_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    void av1_inv_txfm2d_add_64x64_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_64x64)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, int32_t bd);


    void av1_highbd_inv_txfm_add_avx2(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_8x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_8x16)(const int32_t *input, uint16_t *output, int32_t stride, const TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_16x8_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_16x8)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_16x32_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_16x32)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_32x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_32x16)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_32x8_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_32x8)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_8x32_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_8x32)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_32x64_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_32x64)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_64x32_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_64x32)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_16x64_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_16x64)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_64x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_64x16)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t eob, int32_t bd);

    void av1_inv_txfm2d_add_4x8_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    void av1_inv_txfm2d_add_4x8_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_4x8)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);

    void av1_inv_txfm2d_add_8x4_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    void av1_inv_txfm2d_add_8x4_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_8x4)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);

    void av1_inv_txfm2d_add_4x16_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    void av1_inv_txfm2d_add_4x16_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_4x16)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);

    void av1_inv_txfm2d_add_16x4_c(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    void av1_inv_txfm2d_add_16x4_sse4_1(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);
    RTCD_EXTERN void(*av1_inv_txfm2d_add_16x4)(const int32_t *input, uint16_t *output, int32_t stride, TxType tx_type, TxSize tx_size, int32_t bd);



    void av1_inv_txfm_add_c(const tran_low_t *dqcoeff, uint8_t *dst, int32_t stride, const TxfmParam *txfm_param);
    void av1_inv_txfm_add_ssse3(const tran_low_t *dqcoeff, uint8_t *dst, int32_t stride, const TxfmParam *txfm_param);
    RTCD_EXTERN void(*av1_inv_txfm_add)(const tran_low_t *dqcoeff, uint8_t *dst, int32_t stride, const TxfmParam *txfm_param);

#if RS_10BIT_FIX
    //uint32_t aom_highbd_8_mse16x16_c(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
    void      aom_highbd_8_mse16x16_sse2(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
    RTCD_EXTERN void(*aom_highbd_8_mse16x16)(const uint8_t *src_ptr, int32_t  source_stride, const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse);
#endif

#if INTRA_ASM
    void av1_upsample_intra_edge_c(uint8_t *p, int32_t sz);
    void av1_upsample_intra_edge_sse4_1(uint8_t *p, int32_t sz);
    RTCD_EXTERN void(*av1_upsample_intra_edge)(uint8_t *p, int32_t sz);


#if INTRA_10BIT_SUPPORT
    // AMIR
    //void av1_upsample_intra_edge_high_c(uint16_t *p, int32_t sz, int32_t bd);
    //void av1_upsample_intra_edge_high_sse4_1(uint16_t *p, int32_t sz, int32_t bd);
    //RTCD_EXTERN void(*av1_upsample_intra_edge_high)(uint16_t *p, int32_t sz, int32_t bd);
#endif
/* DC_PRED top */

    void aom_dc_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);



    /* DC_PRED top */

    void aom_dc_left_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_left_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_left_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_left_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);



    /* DC_PRED top */

    void aom_dc_top_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_top_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_top_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_top_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);




    /* DC_PRED 128 */
    void aom_dc_128_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_dc_128_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_dc_128_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_dc_128_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);




    /* SMOOTH_H_PRED */
    void aom_smooth_h_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_64x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_32x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_16x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_8x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_4x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_16x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_16x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_16x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_16x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_h_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_32x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_h_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_32x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_32x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_4x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_h_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_4x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_64x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_64x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_h_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_8x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_8x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_h_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_h_predictor_8x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_h_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    /* SMOOTH_V_PRED */


    void aom_smooth_v_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_64x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_32x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_16x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_8x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_4x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_16x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_16x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_16x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_16x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_v_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_32x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_v_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_32x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_32x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_4x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_v_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_4x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_64x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_64x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_v_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_8x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_8x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_v_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_v_predictor_8x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_v_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);



    /* SMOOTH_PRED */

    void aom_smooth_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_64x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_32x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_16x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_8x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_4x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_16x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_16x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_16x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_16x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_32x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_32x64_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_32x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_4x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_4x8_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_64x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_64x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_smooth_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_8x16_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_8x32_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_smooth_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_smooth_predictor_8x4_ssse3(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_smooth_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);







    /* V_PRED */
    void aom_v_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_v_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_v_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_v_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_v_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);


    void aom_v_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_v_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_v_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_v_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);



    /* H_PRED */

    void aom_h_predictor_4x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_4x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_4x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_8x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_8x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_8x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_16x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_16x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_16x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_32x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_32x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_64x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_64x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_64x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_16x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_16x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_16x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_16x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_16x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_16x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_16x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_16x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_16x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_16x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_16x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_16x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_32x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_32x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_32x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_32x64_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_32x64_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_32x64)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_32x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_32x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_32x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_4x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_4x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_4x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_4x8_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_4x8_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_4x8)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_64x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_64x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_64x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_64x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_64x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_64x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_8x16_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_8x16_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_8x16)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_8x32_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_8x32_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_8x32)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);

    void aom_h_predictor_8x4_c(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    void aom_h_predictor_8x4_sse2(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);
    RTCD_EXTERN void(*aom_h_predictor_8x4)(uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left);




#endif

    void aom_ifft16x16_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_ifft16x16_float)(const float *input, float *temp, float *output);

    void aom_ifft2x2_float_c(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_ifft2x2_float)(const float *input, float *temp, float *output);

    void aom_ifft32x32_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_ifft32x32_float)(const float *input, float *temp, float *output);

    void aom_ifft4x4_float_sse2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_ifft4x4_float)(const float *input, float *temp, float *output);

    void aom_ifft8x8_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_ifft8x8_float)(const float *input, float *temp, float *output);

    void aom_fft16x16_float_c(const float *input, float *temp, float *output);
    void aom_fft16x16_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_fft16x16_float)(const float *input, float *temp, float *output);

    void aom_fft2x2_float_c(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_fft2x2_float)(const float *input, float *temp, float *output);

    void aom_fft32x32_float_c(const float *input, float *temp, float *output);
    void aom_fft32x32_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_fft32x32_float)(const float *input, float *temp, float *output);

    void aom_fft4x4_float_c(const float *input, float *temp, float *output);
    void aom_fft4x4_float_sse2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_fft4x4_float)(const float *input, float *temp, float *output);

    void aom_fft8x8_float_c(const float *input, float *temp, float *output);
    void aom_fft8x8_float_avx2(const float *input, float *temp, float *output);
    RTCD_EXTERN void(*aom_fft8x8_float)(const float *input, float *temp, float *output);

#if INTRA_10BIT_SUPPORT
    void aom_highbd_dc_128_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_128_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_128_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_128_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_left_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_left_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_left_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_dc_top_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_dc_top_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_dc_top_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);


    void aom_highbd_h_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_16x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_16x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_16x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_32x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_32x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_h_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_h_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_h_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_4x16_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_4x4_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_4x8_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_8x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_8x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_8x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_h_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_h_predictor_8x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_h_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_4x16_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_4x4_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_4x8_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_8x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_8x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_8x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_predictor_8x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);


    void aom_highbd_smooth_v_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_4x16_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_4x4_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_4x8_ssse3(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_8x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_8x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_8x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_smooth_v_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_smooth_v_predictor_8x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_smooth_v_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

 
    void aom_highbd_v_predictor_16x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_16x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_16x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_16x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_16x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_16x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_16x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_16x4_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_16x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);


    void aom_highbd_v_predictor_16x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_16x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_16x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_16x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_16x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_16x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_2x2_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_2x2)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_32x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_32x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_32x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_32x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_32x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_32x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_32x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_32x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_32x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_32x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_32x8_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_32x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_4x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_4x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_4x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_4x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_4x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_4x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_4x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_4x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_4x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_64x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_64x16_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_64x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_64x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_64x32_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_64x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_64x64_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_64x64_avx2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_64x64)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_8x16_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_8x16_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_8x16)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_8x32_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_8x32_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_8x32)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_8x4_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_8x4_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_8x4)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void aom_highbd_v_predictor_8x8_c(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    void aom_highbd_v_predictor_8x8_sse2(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);
    RTCD_EXTERN void(*aom_highbd_v_predictor_8x8)(uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int32_t bd);

    void av1_dr_prediction_z1_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t dx, int32_t dy);
    void av1_dr_prediction_z1_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t dx, int32_t dy);
    RTCD_EXTERN void(*av1_dr_prediction_z1)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t dx, int32_t dy);
    
    void av1_dr_prediction_z2_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy);
    void av1_dr_prediction_z2_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy);
    RTCD_EXTERN void(*av1_dr_prediction_z2)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy);
    
    void av1_dr_prediction_z3_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_left, int32_t dx, int32_t dy);
    void av1_dr_prediction_z3_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_left, int32_t dx, int32_t dy);
    RTCD_EXTERN void(*av1_dr_prediction_z3)(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above, const uint8_t *left, int32_t upsample_left, int32_t dx, int32_t dy);
    
    void av1_highbd_dr_prediction_z1_c(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t dx, int32_t dy, int32_t bd);
    void av1_highbd_dr_prediction_z1_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t dx, int32_t dy, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_dr_prediction_z1)(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t dx, int32_t dy, int32_t bd);
    
    void av1_highbd_dr_prediction_z2_c(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    void av1_highbd_dr_prediction_z2_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_dr_prediction_z2)(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    
    void av1_highbd_dr_prediction_z3_c(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    void av1_highbd_dr_prediction_z3_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    RTCD_EXTERN void(*av1_highbd_dr_prediction_z3)(uint16_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd);
    
    //void av1_get_nz_map_contexts_c(const uint8_t *const levels, const int16_t *const scan, const uint16_t eob, const TxSize tx_size, const TX_CLASS tx_class, int8_t *const coeff_contexts);
    void av1_get_nz_map_contexts_sse2(const uint8_t *const levels, const int16_t *const scan, const uint16_t eob, const TxSize tx_size, const TX_CLASS tx_class, int8_t *const coeff_contexts);
    RTCD_EXTERN void(*av1_get_nz_map_contexts)(const uint8_t *const levels, const int16_t *const scan, const uint16_t eob, const TxSize tx_size, const TX_CLASS tx_class, int8_t *const coeff_contexts);
    
    void highbd_variance64_c(const uint8_t *a8, int32_t a_stride, const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h, uint64_t *sse, int64_t *sum);
    void highbd_variance64_avx2(const uint8_t *a8, int32_t a_stride, const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h, uint64_t *sse, int64_t *sum);
    RTCD_EXTERN void (*highbd_variance64)(const uint8_t *a8, int32_t a_stride, const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h, uint64_t *sse, int64_t *sum);
    
    void ResidualKernel_c(uint8_t *input, uint32_t inputStride, uint8_t *pred, uint32_t predStride, int16_t *residual, uint32_t residualStride, uint32_t areaWidth, uint32_t areaHeight);
    void ResidualKernel_avx2(uint8_t *input, uint32_t inputStride, uint8_t *pred, uint32_t predStride, int16_t *residual, uint32_t residualStride, uint32_t areaWidth, uint32_t areaHeight);
    RTCD_EXTERN void(*ResidualKernel)(uint8_t *input, uint32_t inputStride, uint8_t *pred, uint32_t predStride, int16_t *residual, uint32_t residualStride, uint32_t areaWidth, uint32_t areaHeight);
    
    void av1_txb_init_levels_c(const tran_low_t *const coeff, const int32_t width, const int32_t height, uint8_t *const levels);
    void av1_txb_init_levels_avx2(const tran_low_t *const coeff, const int32_t width, const int32_t height, uint8_t *const levels);
    RTCD_EXTERN void(*av1_txb_init_levels)(const tran_low_t *const coeff, const int32_t width, const int32_t height, uint8_t *const levels);

#endif

    void aom_dsp_rtcd(void);



#ifdef RTCD_C

    static void setup_rtcd_internal(EbAsm asm_type)
    {
        int32_t flags = HAS_MMX | HAS_SSE | HAS_SSE2 | HAS_SSE3 | HAS_SSSE3 | HAS_SSE4_1 | HAS_SSE4_2 | HAS_AVX;

        if (asm_type == ASM_AVX2)
            flags |= HAS_AVX2;
        //if (asm_type == ASM_NON_AVX2)
        //    flags = ~HAS_AVX2;

        //to use C: flags=0

        apply_selfguided_restoration = apply_selfguided_restoration_c;
        if (flags & HAS_AVX2) apply_selfguided_restoration = apply_selfguided_restoration_avx2;

        av1_wiener_convolve_add_src = av1_wiener_convolve_add_src_c;
        if (flags & HAS_AVX2) av1_wiener_convolve_add_src = av1_wiener_convolve_add_src_avx2;

        av1_highbd_wiener_convolve_add_src = av1_highbd_wiener_convolve_add_src_c;
        if (flags & HAS_AVX2) av1_highbd_wiener_convolve_add_src = av1_highbd_wiener_convolve_add_src_avx2;

        av1_selfguided_restoration = av1_selfguided_restoration_c;
        if (flags & HAS_AVX2) av1_selfguided_restoration = av1_selfguided_restoration_avx2;


        cdef_find_dir = cdef_find_dir_c;
        if (flags & HAS_AVX2) cdef_find_dir = cdef_find_dir_avx2;

        cdef_filter_block = cdef_filter_block_c;
        if (flags & HAS_AVX2) cdef_filter_block = cdef_filter_block_avx2;

        copy_rect8_8bit_to_16bit = copy_rect8_8bit_to_16bit_c;
        if (flags & HAS_AVX2) copy_rect8_8bit_to_16bit = copy_rect8_8bit_to_16bit_avx2;

        av1_compute_stats = av1_compute_stats_c;
        if (flags & HAS_AVX2) av1_compute_stats = av1_compute_stats_avx2;
        av1_compute_stats_highbd = av1_compute_stats_highbd_c;
        if (flags & HAS_AVX2) av1_compute_stats_highbd = av1_compute_stats_highbd_avx2;
        av1_lowbd_pixel_proj_error = av1_lowbd_pixel_proj_error_c;
        if (flags & HAS_AVX2) av1_lowbd_pixel_proj_error = av1_lowbd_pixel_proj_error_avx2;
        av1_highbd_pixel_proj_error = av1_highbd_pixel_proj_error_c;
        if (flags & HAS_AVX2) av1_highbd_pixel_proj_error = av1_highbd_pixel_proj_error_avx2;
        av1_filter_intra_edge_high = av1_filter_intra_edge_high_c;
        if (flags & HAS_SSE4_1) av1_filter_intra_edge_high = av1_filter_intra_edge_high_sse4_1;
        av1_highbd_convolve_2d_copy_sr = av1_highbd_convolve_2d_copy_sr_c;
        if (flags & HAS_AVX2) av1_highbd_convolve_2d_copy_sr = av1_highbd_convolve_2d_copy_sr_avx2;
        av1_highbd_jnt_convolve_2d_copy = av1_highbd_jnt_convolve_2d_copy_c;
        if (flags & HAS_AVX2) av1_highbd_jnt_convolve_2d_copy = av1_highbd_jnt_convolve_2d_copy_avx2;
        av1_highbd_convolve_y_sr = av1_highbd_convolve_y_sr_c;
        if (flags & HAS_AVX2) av1_highbd_convolve_y_sr = av1_highbd_convolve_y_sr_avx2;
        av1_highbd_convolve_2d_sr = av1_highbd_convolve_2d_sr_c;
        if (flags & HAS_AVX2) av1_highbd_convolve_2d_sr = av1_highbd_convolve_2d_sr_avx2;
        av1_highbd_jnt_convolve_2d = av1_highbd_jnt_convolve_2d_c;
        if (flags & HAS_AVX2) av1_highbd_jnt_convolve_2d = av1_highbd_jnt_convolve_2d_avx2;
        av1_highbd_jnt_convolve_x = av1_highbd_jnt_convolve_x_c;
        if (flags & HAS_AVX2) av1_highbd_jnt_convolve_x = av1_highbd_jnt_convolve_x_avx2;
        av1_highbd_jnt_convolve_y = av1_highbd_jnt_convolve_y_c;
        if (flags & HAS_AVX2) av1_highbd_jnt_convolve_y = av1_highbd_jnt_convolve_y_avx2;
        av1_highbd_convolve_x_sr = av1_highbd_convolve_x_sr_c;
        if (flags & HAS_AVX2) av1_highbd_convolve_x_sr = av1_highbd_convolve_x_sr_avx2;
        subtract_average = subtract_average_c;
        if (flags & HAS_AVX2) subtract_average = subtract_average_avx2;


#if INTRA_10BIT_SUPPORT
        av1_filter_intra_edge = av1_filter_intra_edge_high_c_old;
#else
        av1_filter_intra_edge = av1_filter_intra_edge_high_c;
#endif
        if (flags & HAS_SSE4_1) av1_filter_intra_edge = av1_filter_intra_edge_sse4_1;

        eb_smooth_v_predictor = smooth_v_predictor_c;
        if (flags & HAS_SSSE3) eb_smooth_v_predictor = eb_smooth_v_predictor_all_ssse3;

        eb_smooth_h_predictor = smooth_h_predictor_c;
        if (flags & HAS_SSSE3) eb_smooth_h_predictor = eb_smooth_h_predictor_all_ssse3;

        get_proj_subspace = get_proj_subspace_c;
        if (flags & HAS_AVX2) get_proj_subspace = get_proj_subspace_avx2;

        mse_4x4_16bit = mse_4x4_16bit_c;
        if (flags & HAS_AVX2) mse_4x4_16bit = mse_4x4_16bit_avx2;

        dist_8x8_16bit = dist_8x8_16bit_c;
        if (flags & HAS_AVX2) dist_8x8_16bit = dist_8x8_16bit_avx2;

        search_one_dual = search_one_dual_c;
        if (flags & HAS_AVX2) search_one_dual = search_one_dual_avx2;

        aom_mse16x16 = aom_mse16x16_c;
        if (flags & HAS_AVX2) aom_mse16x16 = aom_mse16x16_avx2;

        av1_convolve_2d_copy_sr = av1_convolve_2d_copy_sr_c;
        if (flags & HAS_AVX2) av1_convolve_2d_copy_sr = av1_convolve_2d_copy_sr_avx2;


        av1_convolve_2d_sr = av1_convolve_2d_sr_c;
        if (flags & HAS_AVX2) av1_convolve_2d_sr = av1_convolve_2d_sr_avx2;

        av1_jnt_convolve_2d_copy = av1_jnt_convolve_2d_copy_c;
        if (flags & HAS_AVX2) av1_jnt_convolve_2d_copy = av1_jnt_convolve_2d_copy_avx2;

        av1_convolve_x_sr = av1_convolve_x_sr_c;
        if (flags & HAS_AVX2) av1_convolve_x_sr = av1_convolve_x_sr_avx2;
        av1_convolve_y_sr = av1_convolve_y_sr_c;
        if (flags & HAS_AVX2) av1_convolve_y_sr = av1_convolve_y_sr_avx2;

        av1_jnt_convolve_x = av1_jnt_convolve_x_c;
        if (flags & HAS_AVX2) av1_jnt_convolve_x = av1_jnt_convolve_x_avx2;
        av1_jnt_convolve_y = av1_jnt_convolve_y_c;
        if (flags & HAS_AVX2) av1_jnt_convolve_y = av1_jnt_convolve_y_avx2;

        av1_jnt_convolve_2d = av1_jnt_convolve_2d_c;
        if (flags & HAS_AVX2) av1_jnt_convolve_2d = av1_jnt_convolve_2d_avx2;

        aom_quantize_b = aom_quantize_b_c_II;
        if (flags & HAS_AVX2) aom_quantize_b = aom_highbd_quantize_b_avx2;

        aom_quantize_b_32x32 = aom_quantize_b_32x32_c_II;
        if (flags & HAS_AVX2) aom_quantize_b_32x32 = aom_highbd_quantize_b_32x32_avx2;

        aom_highbd_quantize_b_32x32 = aom_highbd_quantize_b_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_quantize_b_32x32 = aom_highbd_quantize_b_32x32_avx2;

        aom_highbd_quantize_b = aom_highbd_quantize_b_c;
        if (flags & HAS_AVX2) aom_highbd_quantize_b = aom_highbd_quantize_b_avx2;

        av1_inv_txfm2d_add_16x16 = av1_inv_txfm2d_add_16x16_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_16x16 = av1_inv_txfm2d_add_16x16_avx2;
        av1_inv_txfm2d_add_32x32 = av1_inv_txfm2d_add_32x32_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_32x32 = av1_inv_txfm2d_add_32x32_avx2;
        av1_inv_txfm2d_add_4x4 = av1_inv_txfm2d_add_4x4_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_4x4 = av1_inv_txfm2d_add_4x4_avx2;
        av1_inv_txfm2d_add_64x64 = av1_inv_txfm2d_add_64x64_c;
        if (flags & HAS_SSE4_1) av1_inv_txfm2d_add_64x64 = av1_inv_txfm2d_add_64x64_sse4_1;
        av1_inv_txfm2d_add_8x8 = av1_inv_txfm2d_add_8x8_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_8x8 = av1_inv_txfm2d_add_8x8_avx2;

        av1_inv_txfm2d_add_8x16 = av1_inv_txfm2d_add_8x16_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_8x16 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_16x8 = av1_inv_txfm2d_add_16x8_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_16x8 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_16x32 = av1_inv_txfm2d_add_16x32_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_16x32 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_32x16 = av1_inv_txfm2d_add_32x16_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_32x16 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_32x8 = av1_inv_txfm2d_add_32x8_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_32x8 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_8x32 = av1_inv_txfm2d_add_8x32_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_8x32 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_32x64 = av1_inv_txfm2d_add_32x64_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_32x64 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_64x32 = av1_inv_txfm2d_add_64x32_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_64x32 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_16x64 = av1_inv_txfm2d_add_16x64_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_16x64 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_64x16 = av1_inv_txfm2d_add_64x16_c;
        if (flags & HAS_AVX2) av1_inv_txfm2d_add_64x16 = av1_highbd_inv_txfm_add_avx2;
        av1_inv_txfm2d_add_4x8 = av1_inv_txfm2d_add_4x8_c;
        if (flags & HAS_SSE4_1) av1_inv_txfm2d_add_4x8 = av1_inv_txfm2d_add_4x8_sse4_1;
        av1_inv_txfm2d_add_8x4 = av1_inv_txfm2d_add_8x4_c;
        if (flags & HAS_SSE4_1) av1_inv_txfm2d_add_8x4 = av1_inv_txfm2d_add_8x4_sse4_1;
        av1_inv_txfm2d_add_4x16 = av1_inv_txfm2d_add_4x16_c;
        if (flags & HAS_SSE4_1) av1_inv_txfm2d_add_4x16 = av1_inv_txfm2d_add_4x16_sse4_1;
        av1_inv_txfm2d_add_16x4 = av1_inv_txfm2d_add_16x4_c;
        if (flags & HAS_SSE4_1) av1_inv_txfm2d_add_16x4 = av1_inv_txfm2d_add_16x4_sse4_1;

        av1_inv_txfm_add = av1_inv_txfm_add_c;
        if (flags & HAS_SSSE3) av1_inv_txfm_add = av1_inv_txfm_add_ssse3;

        highbd_variance64 = highbd_variance64_c;
        if (flags & HAS_AVX2) highbd_variance64 = highbd_variance64_avx2;

        //aom_highbd_8_mse16x16 = aom_highbd_8_mse16x16_c;
        if (flags & HAS_SSE2) aom_highbd_8_mse16x16 = aom_highbd_8_mse16x16_sse2;

#if INTRA_ASM
        av1_upsample_intra_edge = av1_upsample_intra_edge_c;
        if (flags & HAS_SSE4_1) av1_upsample_intra_edge = av1_upsample_intra_edge_sse4_1;

    
#if INTRINSIC_OPT_2

        aom_highbd_smooth_v_predictor_16x16 = aom_highbd_smooth_v_predictor_16x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_16x16 = aom_highbd_smooth_v_predictor_16x16_avx2;
        aom_highbd_smooth_v_predictor_16x32 = aom_highbd_smooth_v_predictor_16x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_16x32 = aom_highbd_smooth_v_predictor_16x32_avx2;
        aom_highbd_smooth_v_predictor_16x4 = aom_highbd_smooth_v_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_16x4 = aom_highbd_smooth_v_predictor_16x4_avx2;
        aom_highbd_smooth_v_predictor_16x64 = aom_highbd_smooth_v_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_16x64 = aom_highbd_smooth_v_predictor_16x64_avx2;
        aom_highbd_smooth_v_predictor_16x8 = aom_highbd_smooth_v_predictor_16x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_16x8 = aom_highbd_smooth_v_predictor_16x8_avx2;
        aom_highbd_smooth_v_predictor_2x2 = aom_highbd_smooth_v_predictor_2x2_c;
        aom_highbd_smooth_v_predictor_32x16 = aom_highbd_smooth_v_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_32x16 = aom_highbd_smooth_v_predictor_32x16_avx2;
        aom_highbd_smooth_v_predictor_32x32 = aom_highbd_smooth_v_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_32x32 = aom_highbd_smooth_v_predictor_32x32_avx2;
        aom_highbd_smooth_v_predictor_32x64 = aom_highbd_smooth_v_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_32x64 = aom_highbd_smooth_v_predictor_32x64_avx2;
        aom_highbd_smooth_v_predictor_32x8 = aom_highbd_smooth_v_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_32x8 = aom_highbd_smooth_v_predictor_32x8_avx2;
        aom_highbd_smooth_v_predictor_4x16 = aom_highbd_smooth_v_predictor_4x16_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_v_predictor_4x16 = aom_highbd_smooth_v_predictor_4x16_ssse3;
        aom_highbd_smooth_v_predictor_4x4 = aom_highbd_smooth_v_predictor_4x4_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_v_predictor_4x4 = aom_highbd_smooth_v_predictor_4x4_ssse3;
        aom_highbd_smooth_v_predictor_4x8 = aom_highbd_smooth_v_predictor_4x8_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_v_predictor_4x8 = aom_highbd_smooth_v_predictor_4x8_ssse3;
        aom_highbd_smooth_v_predictor_64x16 = aom_highbd_smooth_v_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_64x16 = aom_highbd_smooth_v_predictor_64x16_avx2;
        aom_highbd_smooth_v_predictor_64x32 = aom_highbd_smooth_v_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_64x32 = aom_highbd_smooth_v_predictor_64x32_avx2;
        aom_highbd_smooth_v_predictor_64x64 = aom_highbd_smooth_v_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_64x64 = aom_highbd_smooth_v_predictor_64x64_avx2;
        aom_highbd_smooth_v_predictor_8x16 = aom_highbd_smooth_v_predictor_8x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_8x16 = aom_highbd_smooth_v_predictor_8x16_avx2;
        aom_highbd_smooth_v_predictor_8x32 = aom_highbd_smooth_v_predictor_8x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_8x32 = aom_highbd_smooth_v_predictor_8x32_avx2;
        aom_highbd_smooth_v_predictor_8x4 = aom_highbd_smooth_v_predictor_8x4_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_8x4 = aom_highbd_smooth_v_predictor_8x4_avx2;
        aom_highbd_smooth_v_predictor_8x8 = aom_highbd_smooth_v_predictor_8x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_v_predictor_8x8 = aom_highbd_smooth_v_predictor_8x8_avx2;

        cfl_predict_lbd = cfl_predict_lbd_c;
        if (flags & HAS_AVX2) cfl_predict_lbd = cfl_predict_lbd_avx2;
        cfl_predict_hbd = cfl_predict_hbd_c;
        if (flags & HAS_AVX2) cfl_predict_hbd = cfl_predict_hbd_avx2;

        av1_dr_prediction_z1 = av1_dr_prediction_z1_c;
        if (flags & HAS_AVX2) av1_dr_prediction_z1 = av1_dr_prediction_z1_avx2;
        av1_dr_prediction_z2 = av1_dr_prediction_z2_c;
        if (flags & HAS_AVX2) av1_dr_prediction_z2 = av1_dr_prediction_z2_avx2;
        av1_dr_prediction_z3 = av1_dr_prediction_z3_c;
        if (flags & HAS_AVX2) av1_dr_prediction_z3 = av1_dr_prediction_z3_avx2;
        av1_highbd_dr_prediction_z1 = av1_highbd_dr_prediction_z1_c;
        if (flags & HAS_AVX2) av1_highbd_dr_prediction_z1 = av1_highbd_dr_prediction_z1_avx2;
        av1_highbd_dr_prediction_z2 = av1_highbd_dr_prediction_z2_c;
        if (flags & HAS_AVX2) av1_highbd_dr_prediction_z2 = av1_highbd_dr_prediction_z2_avx2;
        av1_highbd_dr_prediction_z3 = av1_highbd_dr_prediction_z3_c;
        if (flags & HAS_AVX2) av1_highbd_dr_prediction_z3 = av1_highbd_dr_prediction_z3_avx2;
        //av1_get_nz_map_contexts = av1_get_nz_map_contexts_c;
        if (flags & HAS_SSE2) av1_get_nz_map_contexts = av1_get_nz_map_contexts_sse2;

        ResidualKernel = ResidualKernel_c;
        if (flags & HAS_AVX2) ResidualKernel = ResidualKernel_avx2;

        av1_txb_init_levels = av1_txb_init_levels_c;
        if (flags & HAS_AVX2) av1_txb_init_levels = av1_txb_init_levels_avx2;

#else
        aom_highbd_dc_128_predictor_16x16 = aom_highbd_dc_128_predictor_16x16_c;
        aom_highbd_dc_128_predictor_16x32 = aom_highbd_dc_128_predictor_16x32_c;
        aom_highbd_dc_128_predictor_16x4 = aom_highbd_dc_128_predictor_16x4_c;
        aom_highbd_dc_128_predictor_16x64 = aom_highbd_dc_128_predictor_16x64_c;
        aom_highbd_dc_128_predictor_16x8 = aom_highbd_dc_128_predictor_16x8_c;
        aom_highbd_dc_128_predictor_32x16 = aom_highbd_dc_128_predictor_32x16_c;
        aom_highbd_dc_128_predictor_32x32 = aom_highbd_dc_128_predictor_32x32_c;
        aom_highbd_dc_128_predictor_32x64 = aom_highbd_dc_128_predictor_32x64_c;
        aom_highbd_dc_128_predictor_32x8 = aom_highbd_dc_128_predictor_32x8_c;
        aom_highbd_dc_128_predictor_4x16 = aom_highbd_dc_128_predictor_4x16_c;
        aom_highbd_dc_128_predictor_8x32 = aom_highbd_dc_128_predictor_8x32_c;


        aom_highbd_dc_left_predictor_16x16 = aom_highbd_dc_left_predictor_16x16_c;
        aom_highbd_dc_left_predictor_16x32 = aom_highbd_dc_left_predictor_16x32_c;
        aom_highbd_dc_left_predictor_16x4 = aom_highbd_dc_left_predictor_16x4_c;
        aom_highbd_dc_left_predictor_16x64 = aom_highbd_dc_left_predictor_16x64_c;
        aom_highbd_dc_left_predictor_16x8 = aom_highbd_dc_left_predictor_16x8_c;
        aom_highbd_dc_left_predictor_32x16 = aom_highbd_dc_left_predictor_32x16_c;
        aom_highbd_dc_left_predictor_32x32 = aom_highbd_dc_left_predictor_32x32_c;
        aom_highbd_dc_left_predictor_32x64 = aom_highbd_dc_left_predictor_32x64_c;
        aom_highbd_dc_left_predictor_32x8 = aom_highbd_dc_left_predictor_32x8_c;
        aom_highbd_dc_left_predictor_4x16 = aom_highbd_dc_left_predictor_4x16_c;
        aom_highbd_dc_left_predictor_8x32 = aom_highbd_dc_left_predictor_8x32_c;


        aom_highbd_dc_predictor_16x16 = aom_highbd_dc_predictor_16x16_c;
        aom_highbd_dc_predictor_16x32 = aom_highbd_dc_predictor_16x32_c;
        aom_highbd_dc_predictor_16x4 = aom_highbd_dc_predictor_16x4_c;
        aom_highbd_dc_predictor_16x64 = aom_highbd_dc_predictor_16x64_c;
        aom_highbd_dc_predictor_16x8 = aom_highbd_dc_predictor_16x8_c;
        aom_highbd_dc_predictor_32x16 = aom_highbd_dc_predictor_32x16_c;
        aom_highbd_dc_predictor_32x32 = aom_highbd_dc_predictor_32x32_c;
        aom_highbd_dc_predictor_32x64 = aom_highbd_dc_predictor_32x64_c;
        aom_highbd_dc_predictor_32x8 = aom_highbd_dc_predictor_32x8_c;
        aom_highbd_dc_predictor_4x16 = aom_highbd_dc_predictor_4x16_c;
        aom_highbd_dc_predictor_64x16 = aom_highbd_dc_predictor_64x16_c;
        aom_highbd_dc_predictor_64x32 = aom_highbd_dc_predictor_64x32_c;
        aom_highbd_dc_predictor_64x64 = aom_highbd_dc_predictor_64x64_c;
        aom_highbd_dc_predictor_8x32 = aom_highbd_dc_predictor_8x32_c;

        aom_highbd_dc_top_predictor_16x16 = aom_highbd_dc_top_predictor_16x16_c;
        aom_highbd_dc_top_predictor_16x32 = aom_highbd_dc_top_predictor_16x32_c;
        aom_highbd_dc_top_predictor_16x4 = aom_highbd_dc_top_predictor_16x4_c;
        aom_highbd_dc_top_predictor_16x64 = aom_highbd_dc_top_predictor_16x64_c;
        aom_highbd_dc_top_predictor_16x8 = aom_highbd_dc_top_predictor_16x8_c;
        aom_highbd_dc_top_predictor_32x16 = aom_highbd_dc_top_predictor_32x16_c;
        aom_highbd_dc_top_predictor_32x32 = aom_highbd_dc_top_predictor_32x32_c;
        aom_highbd_dc_top_predictor_32x64 = aom_highbd_dc_top_predictor_32x64_c;
        aom_highbd_dc_top_predictor_32x8 = aom_highbd_dc_top_predictor_32x8_c;
        aom_highbd_dc_top_predictor_4x16 = aom_highbd_dc_top_predictor_4x16_c;


        aom_highbd_smooth_v_predictor_16x16 = aom_highbd_smooth_v_predictor_16x16_c;
        aom_highbd_smooth_v_predictor_16x32 = aom_highbd_smooth_v_predictor_16x32_c;
        aom_highbd_smooth_v_predictor_16x4 = aom_highbd_smooth_v_predictor_16x4_c;
        aom_highbd_smooth_v_predictor_16x64 = aom_highbd_smooth_v_predictor_16x64_c;
        aom_highbd_smooth_v_predictor_16x8 = aom_highbd_smooth_v_predictor_16x8_c;
        aom_highbd_smooth_v_predictor_32x16 = aom_highbd_smooth_v_predictor_32x16_c;
        aom_highbd_smooth_v_predictor_32x32 = aom_highbd_smooth_v_predictor_32x32_c;
        aom_highbd_smooth_v_predictor_32x64 = aom_highbd_smooth_v_predictor_32x64_c;
        aom_highbd_smooth_v_predictor_32x8 = aom_highbd_smooth_v_predictor_32x8_c;
        aom_highbd_smooth_v_predictor_4x16 = aom_highbd_smooth_v_predictor_4x16_c;
        aom_highbd_smooth_v_predictor_4x4 = aom_highbd_smooth_v_predictor_4x4_c;
        aom_highbd_smooth_v_predictor_4x8 = aom_highbd_smooth_v_predictor_4x8_c;
        aom_highbd_smooth_v_predictor_64x16 = aom_highbd_smooth_v_predictor_64x16_c;
        aom_highbd_smooth_v_predictor_64x32 = aom_highbd_smooth_v_predictor_64x32_c;
        aom_highbd_smooth_v_predictor_64x64 = aom_highbd_smooth_v_predictor_64x64_c;
        aom_highbd_smooth_v_predictor_8x16 = aom_highbd_smooth_v_predictor_8x16_c;
        aom_highbd_smooth_v_predictor_8x32 = aom_highbd_smooth_v_predictor_8x32_c;
        aom_highbd_smooth_v_predictor_8x4 = aom_highbd_smooth_v_predictor_8x4_c;
        aom_highbd_smooth_v_predictor_8x8 = aom_highbd_smooth_v_predictor_8x8_c;

        cfl_predict_lbd = cfl_predict_lbd_c;
        cfl_predict_hbd = cfl_predict_hbd_c;

        av1_dr_prediction_z1 = av1_dr_prediction_z1_c;
        av1_dr_prediction_z2 = av1_dr_prediction_z2_c;
        av1_dr_prediction_z3 = av1_dr_prediction_z3_c;

        av1_highbd_dr_prediction_z1 = av1_highbd_dr_prediction_z1_c;
        av1_highbd_dr_prediction_z2 = av1_highbd_dr_prediction_z2_c;
        av1_highbd_dr_prediction_z3 = av1_highbd_dr_prediction_z3_c;
        ResidualKernel = ResidualKernel_c;

        av1_txb_init_levels = av1_txb_init_levels_c;
#endif
        aom_dc_predictor_4x4 = aom_dc_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_dc_predictor_4x4 = aom_dc_predictor_4x4_sse2;
        aom_dc_predictor_8x8 = aom_dc_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_dc_predictor_8x8 = aom_dc_predictor_8x8_sse2;
        aom_dc_predictor_16x16 = aom_dc_predictor_16x16_c;
        if (flags & HAS_SSE2) aom_dc_predictor_16x16 = aom_dc_predictor_16x16_sse2;
        aom_dc_predictor_32x32 = aom_dc_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_dc_predictor_32x32 = aom_dc_predictor_32x32_avx2;
        aom_dc_predictor_64x64 = aom_dc_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_dc_predictor_64x64 = aom_dc_predictor_64x64_avx2;
        aom_dc_predictor_32x16 = aom_dc_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_dc_predictor_32x16 = aom_dc_predictor_32x16_avx2;
        aom_dc_predictor_32x64 = aom_dc_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_dc_predictor_32x64 = aom_dc_predictor_32x64_avx2;
        aom_dc_predictor_64x16 = aom_dc_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_dc_predictor_64x16 = aom_dc_predictor_64x16_avx2;
        aom_dc_predictor_8x16 = aom_dc_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_dc_predictor_8x16 = aom_dc_predictor_8x16_sse2;
        aom_dc_predictor_8x32 = aom_dc_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_dc_predictor_8x32 = aom_dc_predictor_8x32_sse2;
        aom_dc_predictor_8x4 = aom_dc_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_dc_predictor_8x4 = aom_dc_predictor_8x4_sse2;
        aom_dc_predictor_64x32 = aom_dc_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_dc_predictor_64x32 = aom_dc_predictor_64x32_avx2;
        aom_dc_predictor_16x32 = aom_dc_predictor_16x32_c;
        if (flags & HAS_SSE2) aom_dc_predictor_16x32 = aom_dc_predictor_16x32_sse2;
        aom_dc_predictor_16x4 = aom_dc_predictor_16x4_c;
        if (flags & HAS_SSE2) aom_dc_predictor_16x4 = aom_dc_predictor_16x4_sse2;
        aom_dc_predictor_16x64 = aom_dc_predictor_16x64_c;
        if (flags & HAS_SSE2) aom_dc_predictor_16x64 = aom_dc_predictor_16x64_sse2;
        aom_dc_predictor_16x8 = aom_dc_predictor_16x8_c;
        if (flags & HAS_SSE2) aom_dc_predictor_16x8 = aom_dc_predictor_16x8_sse2;
        aom_dc_predictor_32x8 = aom_dc_predictor_32x8_c;
        if (flags & HAS_SSE2) aom_dc_predictor_32x8 = aom_dc_predictor_32x8_sse2;
        aom_dc_predictor_4x16 = aom_dc_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_dc_predictor_4x16 = aom_dc_predictor_4x16_sse2;
        aom_dc_predictor_4x8 = aom_dc_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_dc_predictor_4x8 = aom_dc_predictor_4x8_sse2;

        aom_dc_top_predictor_4x4 = aom_dc_top_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_4x4 = aom_dc_top_predictor_4x4_sse2;
        aom_dc_top_predictor_8x8 = aom_dc_top_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_8x8 = aom_dc_top_predictor_8x8_sse2;
        aom_dc_top_predictor_16x16 = aom_dc_top_predictor_16x16_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_16x16 = aom_dc_top_predictor_16x16_sse2;
        aom_dc_top_predictor_32x32 = aom_dc_top_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_dc_top_predictor_32x32 = aom_dc_top_predictor_32x32_avx2;
        aom_dc_top_predictor_64x64 = aom_dc_top_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_dc_top_predictor_64x64 = aom_dc_top_predictor_64x64_avx2;
        aom_dc_top_predictor_16x32 = aom_dc_top_predictor_16x32_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_16x32 = aom_dc_top_predictor_16x32_sse2;
        aom_dc_top_predictor_16x4 = aom_dc_top_predictor_16x4_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_16x4 = aom_dc_top_predictor_16x4_sse2;
        aom_dc_top_predictor_16x64 = aom_dc_top_predictor_16x64_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_16x64 = aom_dc_top_predictor_16x64_sse2;
        aom_dc_top_predictor_16x8 = aom_dc_top_predictor_16x8_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_16x8 = aom_dc_top_predictor_16x8_sse2;
        aom_dc_top_predictor_32x16 = aom_dc_top_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_dc_top_predictor_32x16 = aom_dc_top_predictor_32x16_avx2;
        aom_dc_top_predictor_32x64 = aom_dc_top_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_dc_top_predictor_32x64 = aom_dc_top_predictor_32x64_avx2;
        aom_dc_top_predictor_32x8 = aom_dc_top_predictor_32x8_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_32x8 = aom_dc_top_predictor_32x8_sse2;
        aom_dc_top_predictor_4x16 = aom_dc_top_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_4x16 = aom_dc_top_predictor_4x16_sse2;
        aom_dc_top_predictor_4x8 = aom_dc_top_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_4x8 = aom_dc_top_predictor_4x8_sse2;
        aom_dc_top_predictor_64x16 = aom_dc_top_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_dc_top_predictor_64x16 = aom_dc_top_predictor_64x16_avx2;
        aom_dc_top_predictor_64x32 = aom_dc_top_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_dc_top_predictor_64x32 = aom_dc_top_predictor_64x32_avx2;
        aom_dc_top_predictor_8x16 = aom_dc_top_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_8x16 = aom_dc_top_predictor_8x16_sse2;
        aom_dc_top_predictor_8x32 = aom_dc_top_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_8x32 = aom_dc_top_predictor_8x32_sse2;
        aom_dc_top_predictor_8x4 = aom_dc_top_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_dc_top_predictor_8x4 = aom_dc_top_predictor_8x4_sse2;


        aom_dc_left_predictor_4x4 = aom_dc_left_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_4x4 = aom_dc_left_predictor_4x4_sse2;
        aom_dc_left_predictor_8x8 = aom_dc_left_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_8x8 = aom_dc_left_predictor_8x8_sse2;
        aom_dc_left_predictor_16x16 = aom_dc_left_predictor_16x16_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_16x16 = aom_dc_left_predictor_16x16_sse2;
        aom_dc_left_predictor_32x32 = aom_dc_left_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_dc_left_predictor_32x32 = aom_dc_left_predictor_32x32_avx2;
        aom_dc_left_predictor_64x64 = aom_dc_left_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_dc_left_predictor_64x64 = aom_dc_left_predictor_64x64_avx2;
        aom_dc_left_predictor_16x32 = aom_dc_left_predictor_16x32_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_16x32 = aom_dc_left_predictor_16x32_sse2;
        aom_dc_left_predictor_16x4 = aom_dc_left_predictor_16x4_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_16x4 = aom_dc_left_predictor_16x4_sse2;
        aom_dc_left_predictor_16x64 = aom_dc_left_predictor_16x64_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_16x64 = aom_dc_left_predictor_16x64_sse2;
        aom_dc_left_predictor_16x8 = aom_dc_left_predictor_16x8_c;
        if (flags & HAS_SSE2)  aom_dc_left_predictor_16x8 = aom_dc_left_predictor_16x8_sse2;
        aom_dc_left_predictor_32x16 = aom_dc_left_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_dc_left_predictor_32x16 = aom_dc_left_predictor_32x16_avx2;
        aom_dc_left_predictor_32x64 = aom_dc_left_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_dc_left_predictor_32x64 = aom_dc_left_predictor_32x64_avx2;
        aom_dc_left_predictor_64x16 = aom_dc_left_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_dc_left_predictor_64x16 = aom_dc_left_predictor_64x16_avx2;
        aom_dc_left_predictor_64x32 = aom_dc_left_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_dc_left_predictor_64x32 = aom_dc_left_predictor_64x32_avx2;
        aom_dc_left_predictor_32x8 = aom_dc_left_predictor_32x8_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_32x8 = aom_dc_left_predictor_32x8_sse2;
        aom_dc_left_predictor_4x16 = aom_dc_left_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_4x16 = aom_dc_left_predictor_4x16_sse2;
        aom_dc_left_predictor_4x8 = aom_dc_left_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_4x8 = aom_dc_left_predictor_4x8_sse2;
        aom_dc_left_predictor_8x16 = aom_dc_left_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_8x16 = aom_dc_left_predictor_8x16_sse2;
        aom_dc_left_predictor_8x32 = aom_dc_left_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_8x32 = aom_dc_left_predictor_8x32_sse2;
        aom_dc_left_predictor_8x4 = aom_dc_left_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_dc_left_predictor_8x4 = aom_dc_left_predictor_8x4_sse2;


        aom_dc_128_predictor_4x4 = aom_dc_128_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_4x4 = aom_dc_128_predictor_4x4_sse2;
        aom_dc_128_predictor_8x8 = aom_dc_128_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_8x8 = aom_dc_128_predictor_8x8_sse2;
        aom_dc_128_predictor_16x16 = aom_dc_128_predictor_16x16_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_16x16 = aom_dc_128_predictor_16x16_sse2;
        aom_dc_128_predictor_32x32 = aom_dc_128_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_dc_128_predictor_32x32 = aom_dc_128_predictor_32x32_avx2;
        aom_dc_128_predictor_64x64 = aom_dc_128_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_dc_128_predictor_64x64 = aom_dc_128_predictor_64x64_avx2;
        aom_dc_128_predictor_16x32 = aom_dc_128_predictor_16x32_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_16x32 = aom_dc_128_predictor_16x32_sse2;
        aom_dc_128_predictor_16x4 = aom_dc_128_predictor_16x4_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_16x4 = aom_dc_128_predictor_16x4_sse2;
        aom_dc_128_predictor_16x64 = aom_dc_128_predictor_16x64_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_16x64 = aom_dc_128_predictor_16x64_sse2;
        aom_dc_128_predictor_16x8 = aom_dc_128_predictor_16x8_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_16x8 = aom_dc_128_predictor_16x8_sse2;
        aom_dc_128_predictor_32x16 = aom_dc_128_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_dc_128_predictor_32x16 = aom_dc_128_predictor_32x16_avx2;
        aom_dc_128_predictor_32x64 = aom_dc_128_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_dc_128_predictor_32x64 = aom_dc_128_predictor_32x64_avx2;
        aom_dc_128_predictor_32x8 = aom_dc_128_predictor_32x8_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_32x8 = aom_dc_128_predictor_32x8_sse2;
        aom_dc_128_predictor_4x16 = aom_dc_128_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_4x16 = aom_dc_128_predictor_4x16_sse2;
        aom_dc_128_predictor_4x8 = aom_dc_128_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_4x8 = aom_dc_128_predictor_4x8_sse2;
        aom_dc_128_predictor_64x16 = aom_dc_128_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_dc_128_predictor_64x16 = aom_dc_128_predictor_64x16_avx2;
        aom_dc_128_predictor_64x32 = aom_dc_128_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_dc_128_predictor_64x32 = aom_dc_128_predictor_64x32_avx2;
        aom_dc_128_predictor_8x16 = aom_dc_128_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_8x16 = aom_dc_128_predictor_8x16_sse2;
        aom_dc_128_predictor_8x32 = aom_dc_128_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_8x32 = aom_dc_128_predictor_8x32_sse2;
        aom_dc_128_predictor_8x4 = aom_dc_128_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_dc_128_predictor_8x4 = aom_dc_128_predictor_8x4_sse2;

        aom_smooth_h_predictor_16x32 = aom_smooth_h_predictor_16x32_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_16x32 = aom_smooth_h_predictor_16x32_ssse3;
        aom_smooth_h_predictor_16x4 = aom_smooth_h_predictor_16x4_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_16x4 = aom_smooth_h_predictor_16x4_ssse3;
        aom_smooth_h_predictor_16x64 = aom_smooth_h_predictor_16x64_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_16x64 = aom_smooth_h_predictor_16x64_ssse3;
        aom_smooth_h_predictor_16x8 = aom_smooth_h_predictor_16x8_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_16x8 = aom_smooth_h_predictor_16x8_ssse3;
        aom_smooth_h_predictor_32x16 = aom_smooth_h_predictor_32x16_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_32x16 = aom_smooth_h_predictor_32x16_ssse3;
        aom_smooth_h_predictor_32x64 = aom_smooth_h_predictor_32x64_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_32x64 = aom_smooth_h_predictor_32x64_ssse3;
        aom_smooth_h_predictor_32x8 = aom_smooth_h_predictor_32x8_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_32x8 = aom_smooth_h_predictor_32x8_ssse3;
        aom_smooth_h_predictor_4x16 = aom_smooth_h_predictor_4x16_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_4x16 = aom_smooth_h_predictor_4x16_ssse3;
        aom_smooth_h_predictor_4x8 = aom_smooth_h_predictor_4x8_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_4x8 = aom_smooth_h_predictor_4x8_ssse3;
        aom_smooth_h_predictor_64x16 = aom_smooth_h_predictor_64x16_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_64x16 = aom_smooth_h_predictor_64x16_ssse3;
        aom_smooth_h_predictor_64x32 = aom_smooth_h_predictor_64x32_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_64x32 = aom_smooth_h_predictor_64x32_ssse3;
        aom_smooth_h_predictor_8x16 = aom_smooth_h_predictor_8x16_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_8x16 = aom_smooth_h_predictor_8x16_ssse3;
        aom_smooth_h_predictor_8x32 = aom_smooth_h_predictor_8x32_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_8x32 = aom_smooth_h_predictor_8x32_ssse3;
        aom_smooth_h_predictor_8x4 = aom_smooth_h_predictor_8x4_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_8x4 = aom_smooth_h_predictor_8x4_ssse3;
        aom_smooth_h_predictor_64x64 = aom_smooth_h_predictor_64x64_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_64x64 = aom_smooth_h_predictor_64x64_ssse3;
        aom_smooth_h_predictor_32x32 = aom_smooth_h_predictor_32x32_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_32x32 = aom_smooth_h_predictor_32x32_ssse3;
        aom_smooth_h_predictor_16x16 = aom_smooth_h_predictor_16x16_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_16x16 = aom_smooth_h_predictor_16x16_ssse3;
        aom_smooth_h_predictor_8x8 = aom_smooth_h_predictor_8x8_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_8x8 = aom_smooth_h_predictor_8x8_ssse3;
        aom_smooth_h_predictor_4x4 = aom_smooth_h_predictor_4x4_c;
        if (flags & HAS_SSSE3) aom_smooth_h_predictor_4x4 = aom_smooth_h_predictor_4x4_ssse3;


        if (flags & HAS_SSSE3) aom_smooth_v_predictor_16x32 = aom_smooth_v_predictor_16x32_ssse3;
        aom_smooth_v_predictor_16x4 = aom_smooth_v_predictor_16x4_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_16x4 = aom_smooth_v_predictor_16x4_ssse3;
        aom_smooth_v_predictor_16x64 = aom_smooth_v_predictor_16x64_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_16x64 = aom_smooth_v_predictor_16x64_ssse3;
        aom_smooth_v_predictor_16x8 = aom_smooth_v_predictor_16x8_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_16x8 = aom_smooth_v_predictor_16x8_ssse3;
        aom_smooth_v_predictor_32x16 = aom_smooth_v_predictor_32x16_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_32x16 = aom_smooth_v_predictor_32x16_ssse3;
        aom_smooth_v_predictor_32x64 = aom_smooth_v_predictor_32x64_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_32x64 = aom_smooth_v_predictor_32x64_ssse3;
        aom_smooth_v_predictor_32x8 = aom_smooth_v_predictor_32x8_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_32x8 = aom_smooth_v_predictor_32x8_ssse3;
        aom_smooth_v_predictor_4x16 = aom_smooth_v_predictor_4x16_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_4x16 = aom_smooth_v_predictor_4x16_ssse3;
        aom_smooth_v_predictor_4x8 = aom_smooth_v_predictor_4x8_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_4x8 = aom_smooth_v_predictor_4x8_ssse3;
        aom_smooth_v_predictor_64x16 = aom_smooth_v_predictor_64x16_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_64x16 = aom_smooth_v_predictor_64x16_ssse3;
        aom_smooth_v_predictor_64x32 = aom_smooth_v_predictor_64x32_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_64x32 = aom_smooth_v_predictor_64x32_ssse3;
        aom_smooth_v_predictor_8x16 = aom_smooth_v_predictor_8x16_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_8x16 = aom_smooth_v_predictor_8x16_ssse3;
        aom_smooth_v_predictor_8x32 = aom_smooth_v_predictor_8x32_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_8x32 = aom_smooth_v_predictor_8x32_ssse3;
        aom_smooth_v_predictor_8x4 = aom_smooth_v_predictor_8x4_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_8x4 = aom_smooth_v_predictor_8x4_ssse3;
        aom_smooth_v_predictor_64x64 = aom_smooth_v_predictor_64x64_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_64x64 = aom_smooth_v_predictor_64x64_ssse3;
        aom_smooth_v_predictor_32x32 = aom_smooth_v_predictor_32x32_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_32x32 = aom_smooth_v_predictor_32x32_ssse3;
        aom_smooth_v_predictor_16x16 = aom_smooth_v_predictor_16x16_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_16x16 = aom_smooth_v_predictor_16x16_ssse3;
        aom_smooth_v_predictor_8x8 = aom_smooth_v_predictor_8x8_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_8x8 = aom_smooth_v_predictor_8x8_ssse3;
        aom_smooth_v_predictor_4x4 = aom_smooth_v_predictor_4x4_c;
        if (flags & HAS_SSSE3) aom_smooth_v_predictor_4x4 = aom_smooth_v_predictor_4x4_ssse3;

        aom_smooth_predictor_16x32 = aom_smooth_predictor_16x32_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_16x32 = aom_smooth_predictor_16x32_ssse3;
        aom_smooth_predictor_16x4 = aom_smooth_predictor_16x4_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_16x4 = aom_smooth_predictor_16x4_ssse3;
        aom_smooth_predictor_16x64 = aom_smooth_predictor_16x64_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_16x64 = aom_smooth_predictor_16x64_ssse3;
        aom_smooth_predictor_16x8 = aom_smooth_predictor_16x8_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_16x8 = aom_smooth_predictor_16x8_ssse3;
        aom_smooth_predictor_32x16 = aom_smooth_predictor_32x16_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_32x16 = aom_smooth_predictor_32x16_ssse3;
        aom_smooth_predictor_32x64 = aom_smooth_predictor_32x64_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_32x64 = aom_smooth_predictor_32x64_ssse3;
        aom_smooth_predictor_32x8 = aom_smooth_predictor_32x8_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_32x8 = aom_smooth_predictor_32x8_ssse3;
        aom_smooth_predictor_4x16 = aom_smooth_predictor_4x16_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_4x16 = aom_smooth_predictor_4x16_ssse3;
        aom_smooth_predictor_4x8 = aom_smooth_predictor_4x8_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_4x8 = aom_smooth_predictor_4x8_ssse3;
        aom_smooth_predictor_64x16 = aom_smooth_predictor_64x16_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_64x16 = aom_smooth_predictor_64x16_ssse3;
        aom_smooth_predictor_64x32 = aom_smooth_predictor_64x32_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_64x32 = aom_smooth_predictor_64x32_ssse3;
        aom_smooth_predictor_8x16 = aom_smooth_predictor_8x16_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_8x16 = aom_smooth_predictor_8x16_ssse3;
        aom_smooth_predictor_8x32 = aom_smooth_predictor_8x32_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_8x32 = aom_smooth_predictor_8x32_ssse3;
        aom_smooth_predictor_8x4 = aom_smooth_predictor_8x4_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_8x4 = aom_smooth_predictor_8x4_ssse3;
        aom_smooth_predictor_64x64 = aom_smooth_predictor_64x64_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_64x64 = aom_smooth_predictor_64x64_ssse3;
        aom_smooth_predictor_32x32 = aom_smooth_predictor_32x32_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_32x32 = aom_smooth_predictor_32x32_ssse3;
        aom_smooth_predictor_16x16 = aom_smooth_predictor_16x16_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_16x16 = aom_smooth_predictor_16x16_ssse3;
        aom_smooth_predictor_8x8 = aom_smooth_predictor_8x8_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_8x8 = aom_smooth_predictor_8x8_ssse3;
        aom_smooth_predictor_4x4 = aom_smooth_predictor_4x4_c;
        if (flags & HAS_SSSE3) aom_smooth_predictor_4x4 = aom_smooth_predictor_4x4_ssse3;

        aom_v_predictor_4x4 = aom_v_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_v_predictor_4x4 = aom_v_predictor_4x4_sse2;
        aom_v_predictor_8x8 = aom_v_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_v_predictor_8x8 = aom_v_predictor_8x8_sse2;
        aom_v_predictor_16x16 = aom_v_predictor_16x16_c;
        if (flags & HAS_SSE2) aom_v_predictor_16x16 = aom_v_predictor_16x16_sse2;
        aom_v_predictor_32x32 = aom_v_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_v_predictor_32x32 = aom_v_predictor_32x32_avx2;
        aom_v_predictor_64x64 = aom_v_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_v_predictor_64x64 = aom_v_predictor_64x64_avx2;
        aom_v_predictor_16x32 = aom_v_predictor_16x32_c;
        if (flags & HAS_SSE2) aom_v_predictor_16x32 = aom_v_predictor_16x32_sse2;
        aom_v_predictor_16x4 = aom_v_predictor_16x4_c;
        if (flags & HAS_SSE2) aom_v_predictor_16x4 = aom_v_predictor_16x4_sse2;
        aom_v_predictor_16x64 = aom_v_predictor_16x64_c;
        if (flags & HAS_SSE2) aom_v_predictor_16x64 = aom_v_predictor_16x64_sse2;
        aom_v_predictor_16x8 = aom_v_predictor_16x8_c;
        if (flags & HAS_SSE2) aom_v_predictor_16x8 = aom_v_predictor_16x8_sse2;
        aom_v_predictor_32x16 = aom_v_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_v_predictor_32x16 = aom_v_predictor_32x16_avx2;
        aom_v_predictor_32x64 = aom_v_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_v_predictor_32x64 = aom_v_predictor_32x64_avx2;
        aom_v_predictor_32x8 = aom_v_predictor_32x8_c;
        if (flags & HAS_SSE2) aom_v_predictor_32x8 = aom_v_predictor_32x8_sse2;
        aom_v_predictor_4x16 = aom_v_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_v_predictor_4x16 = aom_v_predictor_4x16_sse2;
        aom_v_predictor_4x8 = aom_v_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_v_predictor_4x8 = aom_v_predictor_4x8_sse2;
        aom_v_predictor_64x16 = aom_v_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_v_predictor_64x16 = aom_v_predictor_64x16_avx2;
        aom_v_predictor_64x32 = aom_v_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_v_predictor_64x32 = aom_v_predictor_64x32_avx2;
        aom_v_predictor_8x16 = aom_v_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_v_predictor_8x16 = aom_v_predictor_8x16_sse2;
        aom_v_predictor_8x32 = aom_v_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_v_predictor_8x32 = aom_v_predictor_8x32_sse2;
        aom_v_predictor_8x4 = aom_v_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_v_predictor_8x4 = aom_v_predictor_8x4_sse2;

        aom_h_predictor_4x4 = aom_h_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_h_predictor_4x4 = aom_h_predictor_4x4_sse2;
        aom_h_predictor_8x8 = aom_h_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_h_predictor_8x8 = aom_h_predictor_8x8_sse2;
        aom_h_predictor_16x16 = aom_h_predictor_16x16_c;
        if (flags & HAS_SSE2) aom_h_predictor_16x16 = aom_h_predictor_16x16_sse2;
        aom_h_predictor_32x32 = aom_h_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_h_predictor_32x32 = aom_h_predictor_32x32_avx2;
        aom_h_predictor_64x64 = aom_h_predictor_64x64_c;
        if (flags & HAS_SSE2) aom_h_predictor_64x64 = aom_h_predictor_64x64_sse2;
        aom_h_predictor_16x32 = aom_h_predictor_16x32_c;
        if (flags & HAS_SSE2) aom_h_predictor_16x32 = aom_h_predictor_16x32_sse2;
        aom_h_predictor_16x4 = aom_h_predictor_16x4_c;
        if (flags & HAS_SSE2) aom_h_predictor_16x4 = aom_h_predictor_16x4_sse2;
        aom_h_predictor_16x64 = aom_h_predictor_16x64_c;
        if (flags & HAS_SSE2) aom_h_predictor_16x64 = aom_h_predictor_16x64_sse2;
        aom_h_predictor_16x8 = aom_h_predictor_16x8_c;
        if (flags & HAS_SSE2) aom_h_predictor_16x8 = aom_h_predictor_16x8_sse2;
        aom_h_predictor_32x16 = aom_h_predictor_32x16_c;
        if (flags & HAS_SSE2) aom_h_predictor_32x16 = aom_h_predictor_32x16_sse2;
        aom_h_predictor_32x64 = aom_h_predictor_32x64_c;
        if (flags & HAS_SSE2) aom_h_predictor_32x64 = aom_h_predictor_32x64_sse2;
        aom_h_predictor_32x8 = aom_h_predictor_32x8_c;
        if (flags & HAS_SSE2) aom_h_predictor_32x8 = aom_h_predictor_32x8_sse2;
        aom_h_predictor_4x16 = aom_h_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_h_predictor_4x16 = aom_h_predictor_4x16_sse2;
        aom_h_predictor_4x8 = aom_h_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_h_predictor_4x8 = aom_h_predictor_4x8_sse2;
        aom_h_predictor_64x16 = aom_h_predictor_64x16_c;
        if (flags & HAS_SSE2) aom_h_predictor_64x16 = aom_h_predictor_64x16_sse2;
        aom_h_predictor_64x32 = aom_h_predictor_64x32_c;
        if (flags & HAS_SSE2) aom_h_predictor_64x32 = aom_h_predictor_64x32_sse2;
        aom_h_predictor_8x16 = aom_h_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_h_predictor_8x16 = aom_h_predictor_8x16_sse2;
        aom_h_predictor_8x32 = aom_h_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_h_predictor_8x32 = aom_h_predictor_8x32_sse2;
        aom_h_predictor_8x4 = aom_h_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_h_predictor_8x4 = aom_h_predictor_8x4_sse2;

 
        //QIQ
#if INTRINSIC_OPT_2
        aom_quantize_b_64x64 = aom_quantize_b_64x64_c_II;
        if (flags & HAS_AVX2) aom_quantize_b_64x64 = aom_highbd_quantize_b_64x64_avx2;

        aom_highbd_quantize_b_64x64 = aom_highbd_quantize_b_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_quantize_b_64x64 = aom_highbd_quantize_b_64x64_avx2;
#else
        aom_quantize_b_64x64 = aom_quantize_b_64x64_c_II;
        aom_highbd_quantize_b_64x64 = aom_highbd_quantize_b_64x64_c;

#endif
        // transform
#if INTRINSIC_OPT_2
        av1_fwd_txfm2d_16x8 = av1_fwd_txfm2d_16x8_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_16x8 = av1_fwd_txfm2d_16x8_avx2;
        av1_fwd_txfm2d_8x16 = av1_fwd_txfm2d_8x16_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_8x16 = av1_fwd_txfm2d_8x16_avx2;

        av1_fwd_txfm2d_16x4 = av1_fwd_txfm2d_16x4_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_16x4 = av1_fwd_txfm2d_16x4_avx2;
        av1_fwd_txfm2d_4x16 = av1_fwd_txfm2d_4x16_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_4x16 = av1_fwd_txfm2d_4x16_avx2;

        av1_fwd_txfm2d_8x4 = av1_fwd_txfm2d_8x4_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_8x4 = av1_fwd_txfm2d_8x4_avx2;
        av1_fwd_txfm2d_4x8 = av1_fwd_txfm2d_4x8_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_4x8 = av1_fwd_txfm2d_4x8_avx2;
#else
        av1_fwd_txfm2d_16x8 = av1_fwd_txfm2d_16x8_c;
        av1_fwd_txfm2d_8x16 = av1_fwd_txfm2d_8x16_c;
        av1_fwd_txfm2d_16x4 = av1_fwd_txfm2d_16x4_c;
        av1_fwd_txfm2d_4x16 = av1_fwd_txfm2d_4x16_c;
        av1_fwd_txfm2d_8x4 = av1_fwd_txfm2d_8x4_c;
        av1_fwd_txfm2d_4x8 = av1_fwd_txfm2d_4x8_c;
#endif

        av1_fwd_txfm2d_32x16 = av1_fwd_txfm2d_32x16_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_32x16 = av1_fwd_txfm2d_32x16_avx2;
        av1_fwd_txfm2d_32x8 = av1_fwd_txfm2d_32x8_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_32x8 = av1_fwd_txfm2d_32x8_avx2;
        av1_fwd_txfm2d_8x32 = av1_fwd_txfm2d_8x32_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_8x32 = av1_fwd_txfm2d_8x32_avx2;
        av1_fwd_txfm2d_16x32 = av1_fwd_txfm2d_16x32_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_16x32 = av1_fwd_txfm2d_16x32_avx2;
        av1_fwd_txfm2d_32x64 = av1_fwd_txfm2d_32x64_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_32x64 = av1_fwd_txfm2d_32x64_avx2;
        av1_fwd_txfm2d_64x32 = av1_fwd_txfm2d_64x32_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_64x32 = av1_fwd_txfm2d_64x32_avx2;
        av1_fwd_txfm2d_16x64 = av1_fwd_txfm2d_16x64_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_16x64 = av1_fwd_txfm2d_16x64_avx2;
        av1_fwd_txfm2d_64x16 = av1_fwd_txfm2d_64x16_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_64x16 = av1_fwd_txfm2d_64x16_avx2;
        av1_fwd_txfm2d_64x64 = Av1TransformTwoD_64x64_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_64x64 = av1_fwd_txfm2d_64x64_avx2;
        av1_fwd_txfm2d_32x32 = Av1TransformTwoD_32x32_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_32x32 = av1_fwd_txfm2d_32x32_avx2;
        av1_fwd_txfm2d_16x16 = Av1TransformTwoD_16x16_c;
#if INTRINSIC_OPT_2
        if (flags & HAS_AVX2) av1_fwd_txfm2d_16x16 = av1_fwd_txfm2d_16x16_avx2;
#endif
        av1_fwd_txfm2d_8x8 = Av1TransformTwoD_8x8_c;
        if (flags & HAS_AVX2) av1_fwd_txfm2d_8x8 = av1_fwd_txfm2d_8x8_avx2;
        av1_fwd_txfm2d_4x4 = Av1TransformTwoD_4x4_c;
        if (flags & HAS_SSE4_1) av1_fwd_txfm2d_4x4 = av1_fwd_txfm2d_4x4_sse4_1;


        // aom_highbd_v_predictor
        aom_highbd_v_predictor_16x16 = aom_highbd_v_predictor_16x16_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_16x16 = aom_highbd_v_predictor_16x16_avx2;
        aom_highbd_v_predictor_16x32 = aom_highbd_v_predictor_16x32_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_16x32 = aom_highbd_v_predictor_16x32_avx2;
        aom_highbd_v_predictor_16x4 = aom_highbd_v_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_16x4 = aom_highbd_v_predictor_16x4_avx2;
        aom_highbd_v_predictor_16x64 = aom_highbd_v_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_16x64 = aom_highbd_v_predictor_16x64_avx2;
        aom_highbd_v_predictor_16x8 = aom_highbd_v_predictor_16x8_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_16x8 = aom_highbd_v_predictor_16x8_avx2;
        aom_highbd_v_predictor_2x2 = aom_highbd_v_predictor_2x2_c;
        aom_highbd_v_predictor_32x16 = aom_highbd_v_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_32x16 = aom_highbd_v_predictor_32x16_avx2;
        aom_highbd_v_predictor_32x32 = aom_highbd_v_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_32x32 = aom_highbd_v_predictor_32x32_avx2;
        aom_highbd_v_predictor_32x64 = aom_highbd_v_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_32x64 = aom_highbd_v_predictor_32x64_avx2;
        aom_highbd_v_predictor_32x8 = aom_highbd_v_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_32x8 = aom_highbd_v_predictor_32x8_avx2;
        aom_highbd_v_predictor_4x16 = aom_highbd_v_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_highbd_v_predictor_4x16 = aom_highbd_v_predictor_4x16_sse2;
        aom_highbd_v_predictor_4x4 = aom_highbd_v_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_highbd_v_predictor_4x4 = aom_highbd_v_predictor_4x4_sse2;
        aom_highbd_v_predictor_4x8 = aom_highbd_v_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_highbd_v_predictor_4x8 = aom_highbd_v_predictor_4x8_sse2;
        aom_highbd_v_predictor_64x16 = aom_highbd_v_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_64x16 = aom_highbd_v_predictor_64x16_avx2;
        aom_highbd_v_predictor_64x32 = aom_highbd_v_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_64x32 = aom_highbd_v_predictor_64x32_avx2;
        aom_highbd_v_predictor_8x32 = aom_highbd_v_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_highbd_v_predictor_8x32 = aom_highbd_v_predictor_8x32_sse2;
        aom_highbd_v_predictor_64x64 = aom_highbd_v_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_v_predictor_64x64 = aom_highbd_v_predictor_64x64_avx2;
        aom_highbd_v_predictor_8x16 = aom_highbd_v_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_highbd_v_predictor_8x16 = aom_highbd_v_predictor_8x16_sse2;
        aom_highbd_v_predictor_8x4 = aom_highbd_v_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_highbd_v_predictor_8x4 = aom_highbd_v_predictor_8x4_sse2;
        aom_highbd_v_predictor_8x8 = aom_highbd_v_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_highbd_v_predictor_8x8 = aom_highbd_v_predictor_8x8_sse2;

        //aom_highbd_smooth_predictor
        aom_highbd_smooth_predictor_16x16 = aom_highbd_smooth_predictor_16x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_16x16 = aom_highbd_smooth_predictor_16x16_avx2;
        aom_highbd_smooth_predictor_16x32 = aom_highbd_smooth_predictor_16x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_16x32 = aom_highbd_smooth_predictor_16x32_avx2;
        aom_highbd_smooth_predictor_16x4 = aom_highbd_smooth_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_16x4 = aom_highbd_smooth_predictor_16x4_avx2;
        aom_highbd_smooth_predictor_16x64 = aom_highbd_smooth_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_16x64 = aom_highbd_smooth_predictor_16x64_avx2;
        aom_highbd_smooth_predictor_16x8 = aom_highbd_smooth_predictor_16x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_16x8 = aom_highbd_smooth_predictor_16x8_avx2;
        aom_highbd_smooth_predictor_2x2 = aom_highbd_smooth_predictor_2x2_c;
        aom_highbd_smooth_predictor_32x16 = aom_highbd_smooth_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_32x16 = aom_highbd_smooth_predictor_32x16_avx2;
        aom_highbd_smooth_predictor_32x32 = aom_highbd_smooth_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_32x32 = aom_highbd_smooth_predictor_32x32_avx2;
        aom_highbd_smooth_predictor_32x64 = aom_highbd_smooth_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_32x64 = aom_highbd_smooth_predictor_32x64_avx2;
        aom_highbd_smooth_predictor_32x8 = aom_highbd_smooth_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_32x8 = aom_highbd_smooth_predictor_32x8_avx2;
        aom_highbd_smooth_predictor_4x16 = aom_highbd_smooth_predictor_4x16_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_predictor_4x16 = aom_highbd_smooth_predictor_4x16_ssse3;
        aom_highbd_smooth_predictor_4x4 = aom_highbd_smooth_predictor_4x4_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_predictor_4x4 = aom_highbd_smooth_predictor_4x4_ssse3;
        aom_highbd_smooth_predictor_4x8 = aom_highbd_smooth_predictor_4x8_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_predictor_4x8 = aom_highbd_smooth_predictor_4x8_ssse3;
        aom_highbd_smooth_predictor_64x16 = aom_highbd_smooth_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_64x16 = aom_highbd_smooth_predictor_64x16_avx2;
        aom_highbd_smooth_predictor_64x32 = aom_highbd_smooth_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_64x32 = aom_highbd_smooth_predictor_64x32_avx2;
        aom_highbd_smooth_predictor_64x64 = aom_highbd_smooth_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_64x64 = aom_highbd_smooth_predictor_64x64_avx2;
        aom_highbd_smooth_predictor_8x16 = aom_highbd_smooth_predictor_8x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_8x16 = aom_highbd_smooth_predictor_8x16_avx2;
        aom_highbd_smooth_predictor_8x32 = aom_highbd_smooth_predictor_8x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_8x32 = aom_highbd_smooth_predictor_8x32_avx2;
        aom_highbd_smooth_predictor_8x4 = aom_highbd_smooth_predictor_8x4_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_8x4 = aom_highbd_smooth_predictor_8x4_avx2;
        aom_highbd_smooth_predictor_8x8 = aom_highbd_smooth_predictor_8x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_predictor_8x8 = aom_highbd_smooth_predictor_8x8_avx2;

        //aom_highbd_smooth_h_predictor
        aom_highbd_smooth_h_predictor_16x16 = aom_highbd_smooth_h_predictor_16x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_16x16 = aom_highbd_smooth_h_predictor_16x16_avx2;
        aom_highbd_smooth_h_predictor_16x32 = aom_highbd_smooth_h_predictor_16x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_16x32 = aom_highbd_smooth_h_predictor_16x32_avx2;
        aom_highbd_smooth_h_predictor_16x4 = aom_highbd_smooth_h_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_16x4 = aom_highbd_smooth_h_predictor_16x4_avx2;
        aom_highbd_smooth_h_predictor_16x64 = aom_highbd_smooth_h_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_16x64 = aom_highbd_smooth_h_predictor_16x64_avx2;
        aom_highbd_smooth_h_predictor_16x8 = aom_highbd_smooth_h_predictor_16x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_16x8 = aom_highbd_smooth_h_predictor_16x8_avx2;
        aom_highbd_smooth_h_predictor_2x2 = aom_highbd_smooth_h_predictor_2x2_c;
        aom_highbd_smooth_h_predictor_32x16 = aom_highbd_smooth_h_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_32x16 = aom_highbd_smooth_h_predictor_32x16_avx2;
        aom_highbd_smooth_h_predictor_32x32 = aom_highbd_smooth_h_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_32x32 = aom_highbd_smooth_h_predictor_32x32_avx2;
        aom_highbd_smooth_h_predictor_32x64 = aom_highbd_smooth_h_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_32x64 = aom_highbd_smooth_h_predictor_32x64_avx2;
        aom_highbd_smooth_h_predictor_32x8 = aom_highbd_smooth_h_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_32x8 = aom_highbd_smooth_h_predictor_32x8_avx2;
        aom_highbd_smooth_h_predictor_4x16 = aom_highbd_smooth_h_predictor_4x16_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_h_predictor_4x16 = aom_highbd_smooth_h_predictor_4x16_ssse3;
        aom_highbd_smooth_h_predictor_4x4 = aom_highbd_smooth_h_predictor_4x4_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_h_predictor_4x4 = aom_highbd_smooth_h_predictor_4x4_ssse3;
        aom_highbd_smooth_h_predictor_4x8 = aom_highbd_smooth_h_predictor_4x8_c;
        if (flags & HAS_SSSE3) aom_highbd_smooth_h_predictor_4x8 = aom_highbd_smooth_h_predictor_4x8_ssse3;
        aom_highbd_smooth_h_predictor_64x16 = aom_highbd_smooth_h_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_64x16 = aom_highbd_smooth_h_predictor_64x16_avx2;
        aom_highbd_smooth_h_predictor_64x32 = aom_highbd_smooth_h_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_64x32 = aom_highbd_smooth_h_predictor_64x32_avx2;
        aom_highbd_smooth_h_predictor_64x64 = aom_highbd_smooth_h_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_64x64 = aom_highbd_smooth_h_predictor_64x64_avx2;
        aom_highbd_smooth_h_predictor_8x16 = aom_highbd_smooth_h_predictor_8x16_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_8x16 = aom_highbd_smooth_h_predictor_8x16_avx2;
        aom_highbd_smooth_h_predictor_8x32 = aom_highbd_smooth_h_predictor_8x32_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_8x32 = aom_highbd_smooth_h_predictor_8x32_avx2;
        aom_highbd_smooth_h_predictor_8x4 = aom_highbd_smooth_h_predictor_8x4_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_8x4 = aom_highbd_smooth_h_predictor_8x4_avx2;
        aom_highbd_smooth_h_predictor_8x8 = aom_highbd_smooth_h_predictor_8x8_c;
        if (flags & HAS_AVX2) aom_highbd_smooth_h_predictor_8x8 = aom_highbd_smooth_h_predictor_8x8_avx2;

        //aom_highbd_dc_128_predictor
        aom_highbd_dc_128_predictor_16x16 = aom_highbd_dc_128_predictor_16x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_16x16 = aom_highbd_dc_128_predictor_16x16_avx2;
        aom_highbd_dc_128_predictor_16x32 = aom_highbd_dc_128_predictor_16x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_16x32 = aom_highbd_dc_128_predictor_16x32_avx2;
        aom_highbd_dc_128_predictor_16x4 = aom_highbd_dc_128_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_16x4 = aom_highbd_dc_128_predictor_16x4_avx2;
        aom_highbd_dc_128_predictor_16x64 = aom_highbd_dc_128_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_16x64 = aom_highbd_dc_128_predictor_16x64_avx2;
        aom_highbd_dc_128_predictor_16x8 = aom_highbd_dc_128_predictor_16x8_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_16x8 = aom_highbd_dc_128_predictor_16x8_avx2;
        aom_highbd_dc_128_predictor_2x2 = aom_highbd_dc_128_predictor_2x2_c;
        aom_highbd_dc_128_predictor_32x16 = aom_highbd_dc_128_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_32x16 = aom_highbd_dc_128_predictor_32x16_avx2;
        aom_highbd_dc_128_predictor_32x32 = aom_highbd_dc_128_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_32x32 = aom_highbd_dc_128_predictor_32x32_avx2;
        aom_highbd_dc_128_predictor_32x64 = aom_highbd_dc_128_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_32x64 = aom_highbd_dc_128_predictor_32x64_avx2;
        aom_highbd_dc_128_predictor_32x8 = aom_highbd_dc_128_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_32x8 = aom_highbd_dc_128_predictor_32x8_avx2;
        aom_highbd_dc_128_predictor_4x16 = aom_highbd_dc_128_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_highbd_dc_128_predictor_4x16 = aom_highbd_dc_128_predictor_4x16_sse2;
        aom_highbd_dc_128_predictor_4x4 = aom_highbd_dc_128_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_highbd_dc_128_predictor_4x4 = aom_highbd_dc_128_predictor_4x4_sse2;
        aom_highbd_dc_128_predictor_4x8 = aom_highbd_dc_128_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_highbd_dc_128_predictor_4x8 = aom_highbd_dc_128_predictor_4x8_sse2;
        aom_highbd_dc_128_predictor_8x32 = aom_highbd_dc_128_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_highbd_dc_128_predictor_8x32 = aom_highbd_dc_128_predictor_8x32_sse2;
        aom_highbd_dc_128_predictor_64x16 = aom_highbd_dc_128_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_64x16 = aom_highbd_dc_128_predictor_64x16_avx2;
        aom_highbd_dc_128_predictor_64x32 = aom_highbd_dc_128_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_64x32 = aom_highbd_dc_128_predictor_64x32_avx2;
        aom_highbd_dc_128_predictor_64x64 = aom_highbd_dc_128_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_128_predictor_64x64 = aom_highbd_dc_128_predictor_64x64_avx2;
        aom_highbd_dc_128_predictor_8x16 = aom_highbd_dc_128_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_highbd_dc_128_predictor_8x16 = aom_highbd_dc_128_predictor_8x16_sse2;
        aom_highbd_dc_128_predictor_8x4 = aom_highbd_dc_128_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_highbd_dc_128_predictor_8x4 = aom_highbd_dc_128_predictor_8x4_sse2;
        aom_highbd_dc_128_predictor_8x8 = aom_highbd_dc_128_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_highbd_dc_128_predictor_8x8 = aom_highbd_dc_128_predictor_8x8_sse2;

        //aom_highbd_dc_left_predictor
        aom_highbd_dc_left_predictor_16x16 = aom_highbd_dc_left_predictor_16x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_16x16 = aom_highbd_dc_left_predictor_16x16_avx2;
        aom_highbd_dc_left_predictor_16x32 = aom_highbd_dc_left_predictor_16x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_16x32 = aom_highbd_dc_left_predictor_16x32_avx2;
        aom_highbd_dc_left_predictor_16x4 = aom_highbd_dc_left_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_16x4 = aom_highbd_dc_left_predictor_16x4_avx2;
        aom_highbd_dc_left_predictor_16x64 = aom_highbd_dc_left_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_16x64 = aom_highbd_dc_left_predictor_16x64_avx2;
        aom_highbd_dc_left_predictor_16x8 = aom_highbd_dc_left_predictor_16x8_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_16x8 = aom_highbd_dc_left_predictor_16x8_avx2;
        aom_highbd_dc_left_predictor_2x2 = aom_highbd_dc_left_predictor_2x2_c;
        aom_highbd_dc_left_predictor_32x16 = aom_highbd_dc_left_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_32x16 = aom_highbd_dc_left_predictor_32x16_avx2;
        aom_highbd_dc_left_predictor_32x32 = aom_highbd_dc_left_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_32x32 = aom_highbd_dc_left_predictor_32x32_avx2;
        aom_highbd_dc_left_predictor_32x64 = aom_highbd_dc_left_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_32x64 = aom_highbd_dc_left_predictor_32x64_avx2;
        aom_highbd_dc_left_predictor_32x8 = aom_highbd_dc_left_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_32x8 = aom_highbd_dc_left_predictor_32x8_avx2;
        aom_highbd_dc_left_predictor_4x16 = aom_highbd_dc_left_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_highbd_dc_left_predictor_4x16 = aom_highbd_dc_left_predictor_4x16_sse2;
        aom_highbd_dc_left_predictor_4x4 = aom_highbd_dc_left_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_highbd_dc_left_predictor_4x4 = aom_highbd_dc_left_predictor_4x4_sse2;
        aom_highbd_dc_left_predictor_4x8 = aom_highbd_dc_left_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_highbd_dc_left_predictor_4x8 = aom_highbd_dc_left_predictor_4x8_sse2;
        aom_highbd_dc_left_predictor_8x32 = aom_highbd_dc_left_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_highbd_dc_left_predictor_8x32 = aom_highbd_dc_left_predictor_8x32_sse2;
        aom_highbd_dc_left_predictor_64x16 = aom_highbd_dc_left_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_64x16 = aom_highbd_dc_left_predictor_64x16_avx2;
        aom_highbd_dc_left_predictor_64x32 = aom_highbd_dc_left_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_64x32 = aom_highbd_dc_left_predictor_64x32_avx2;
        aom_highbd_dc_left_predictor_64x64 = aom_highbd_dc_left_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_left_predictor_64x64 = aom_highbd_dc_left_predictor_64x64_avx2;
        aom_highbd_dc_left_predictor_8x16 = aom_highbd_dc_left_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_highbd_dc_left_predictor_8x16 = aom_highbd_dc_left_predictor_8x16_sse2;
        aom_highbd_dc_left_predictor_8x4 = aom_highbd_dc_left_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_highbd_dc_left_predictor_8x4 = aom_highbd_dc_left_predictor_8x4_sse2;
        aom_highbd_dc_left_predictor_8x8 = aom_highbd_dc_left_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_highbd_dc_left_predictor_8x8 = aom_highbd_dc_left_predictor_8x8_sse2;


        aom_highbd_dc_predictor_16x16 = aom_highbd_dc_predictor_16x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_16x16 = aom_highbd_dc_predictor_16x16_avx2;
        aom_highbd_dc_predictor_16x32 = aom_highbd_dc_predictor_16x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_16x32 = aom_highbd_dc_predictor_16x32_avx2;
        aom_highbd_dc_predictor_16x4 = aom_highbd_dc_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_16x4 = aom_highbd_dc_predictor_16x4_avx2;
        aom_highbd_dc_predictor_16x64 = aom_highbd_dc_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_16x64 = aom_highbd_dc_predictor_16x64_avx2;
        aom_highbd_dc_predictor_16x8 = aom_highbd_dc_predictor_16x8_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_16x8 = aom_highbd_dc_predictor_16x8_avx2;
        aom_highbd_dc_predictor_2x2 = aom_highbd_dc_predictor_2x2_c;
        aom_highbd_dc_predictor_32x16 = aom_highbd_dc_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_32x16 = aom_highbd_dc_predictor_32x16_avx2;
        aom_highbd_dc_predictor_32x32 = aom_highbd_dc_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_32x32 = aom_highbd_dc_predictor_32x32_avx2;
        aom_highbd_dc_predictor_32x64 = aom_highbd_dc_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_32x64 = aom_highbd_dc_predictor_32x64_avx2;
        aom_highbd_dc_predictor_32x8 = aom_highbd_dc_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_32x8 = aom_highbd_dc_predictor_32x8_avx2;
        aom_highbd_dc_predictor_4x16 = aom_highbd_dc_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_highbd_dc_predictor_4x16 = aom_highbd_dc_predictor_4x16_sse2;
        aom_highbd_dc_predictor_4x4 = aom_highbd_dc_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_highbd_dc_predictor_4x4 = aom_highbd_dc_predictor_4x4_sse2;
        aom_highbd_dc_predictor_4x8 = aom_highbd_dc_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_highbd_dc_predictor_4x8 = aom_highbd_dc_predictor_4x8_sse2;
        aom_highbd_dc_predictor_64x16 = aom_highbd_dc_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_64x16 = aom_highbd_dc_predictor_64x16_avx2;
        aom_highbd_dc_predictor_64x32 = aom_highbd_dc_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_64x32 = aom_highbd_dc_predictor_64x32_avx2;
        aom_highbd_dc_predictor_64x64 = aom_highbd_dc_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_predictor_64x64 = aom_highbd_dc_predictor_64x64_avx2;
        aom_highbd_dc_predictor_8x16 = aom_highbd_dc_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_highbd_dc_predictor_8x16 = aom_highbd_dc_predictor_8x16_sse2;
        aom_highbd_dc_predictor_8x4 = aom_highbd_dc_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_highbd_dc_predictor_8x4 = aom_highbd_dc_predictor_8x4_sse2;
        aom_highbd_dc_predictor_8x8 = aom_highbd_dc_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_highbd_dc_predictor_8x8 = aom_highbd_dc_predictor_8x8_sse2;
        aom_highbd_dc_predictor_8x32 = aom_highbd_dc_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_highbd_dc_predictor_8x32 = aom_highbd_dc_predictor_8x32_sse2;

        //aom_highbd_dc_top_predictor
        aom_highbd_dc_top_predictor_16x16 = aom_highbd_dc_top_predictor_16x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_16x16 = aom_highbd_dc_top_predictor_16x16_avx2;
        aom_highbd_dc_top_predictor_16x32 = aom_highbd_dc_top_predictor_16x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_16x32 = aom_highbd_dc_top_predictor_16x32_avx2;
        aom_highbd_dc_top_predictor_16x4 = aom_highbd_dc_top_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_16x4 = aom_highbd_dc_top_predictor_16x4_avx2;
        aom_highbd_dc_top_predictor_16x64 = aom_highbd_dc_top_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_16x64 = aom_highbd_dc_top_predictor_16x64_avx2;
        aom_highbd_dc_top_predictor_16x8 = aom_highbd_dc_top_predictor_16x8_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_16x8 = aom_highbd_dc_top_predictor_16x8_avx2;
        aom_highbd_dc_top_predictor_2x2 = aom_highbd_dc_top_predictor_2x2_c;
        aom_highbd_dc_top_predictor_32x16 = aom_highbd_dc_top_predictor_32x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_32x16 = aom_highbd_dc_top_predictor_32x16_avx2;
        aom_highbd_dc_top_predictor_32x32 = aom_highbd_dc_top_predictor_32x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_32x32 = aom_highbd_dc_top_predictor_32x32_avx2;
        aom_highbd_dc_top_predictor_32x64 = aom_highbd_dc_top_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_32x64 = aom_highbd_dc_top_predictor_32x64_avx2;
        aom_highbd_dc_top_predictor_32x8 = aom_highbd_dc_top_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_32x8 = aom_highbd_dc_top_predictor_32x8_avx2;
        aom_highbd_dc_top_predictor_4x16 = aom_highbd_dc_top_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_highbd_dc_top_predictor_4x16 = aom_highbd_dc_top_predictor_4x16_sse2;
        aom_highbd_dc_top_predictor_4x4 = aom_highbd_dc_top_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_highbd_dc_top_predictor_4x4 = aom_highbd_dc_top_predictor_4x4_sse2;
        aom_highbd_dc_top_predictor_4x8 = aom_highbd_dc_top_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_highbd_dc_top_predictor_4x8 = aom_highbd_dc_top_predictor_4x8_sse2;
        aom_highbd_dc_top_predictor_64x16 = aom_highbd_dc_top_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_64x16 = aom_highbd_dc_top_predictor_64x16_avx2;
        aom_highbd_dc_top_predictor_64x32 = aom_highbd_dc_top_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_64x32 = aom_highbd_dc_top_predictor_64x32_avx2;
        aom_highbd_dc_top_predictor_64x64 = aom_highbd_dc_top_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_dc_top_predictor_64x64 = aom_highbd_dc_top_predictor_64x64_avx2;
        aom_highbd_dc_top_predictor_8x16 = aom_highbd_dc_top_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_highbd_dc_top_predictor_8x16 = aom_highbd_dc_top_predictor_8x16_sse2;
        if (flags & HAS_SSE2) aom_highbd_dc_top_predictor_8x32 = aom_highbd_dc_top_predictor_8x32_c;
        aom_highbd_dc_top_predictor_8x4 = aom_highbd_dc_top_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_highbd_dc_top_predictor_8x4 = aom_highbd_dc_top_predictor_8x4_sse2;
        aom_highbd_dc_top_predictor_8x8 = aom_highbd_dc_top_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_highbd_dc_top_predictor_8x8 = aom_highbd_dc_top_predictor_8x8_sse2;


        // aom_highbd_h_predictor
        aom_highbd_h_predictor_16x4 = aom_highbd_h_predictor_16x4_c;
        if (flags & HAS_AVX2) aom_highbd_h_predictor_16x4 = aom_highbd_h_predictor_16x4_avx2;
        aom_highbd_h_predictor_16x64 = aom_highbd_h_predictor_16x64_c;
        if (flags & HAS_AVX2) aom_highbd_h_predictor_16x64 = aom_highbd_h_predictor_16x64_avx2;
        aom_highbd_h_predictor_16x8 = aom_highbd_h_predictor_16x8_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_16x8 = aom_highbd_h_predictor_16x8_sse2;
        aom_highbd_h_predictor_2x2 = aom_highbd_h_predictor_2x2_c;
        aom_highbd_h_predictor_32x16 = aom_highbd_h_predictor_32x16_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_32x16 = aom_highbd_h_predictor_32x16_sse2;
        aom_highbd_h_predictor_32x32 = aom_highbd_h_predictor_32x32_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_32x32 = aom_highbd_h_predictor_32x32_sse2;
        aom_highbd_h_predictor_32x64 = aom_highbd_h_predictor_32x64_c;
        if (flags & HAS_AVX2) aom_highbd_h_predictor_32x64 = aom_highbd_h_predictor_32x64_avx2;
        aom_highbd_h_predictor_32x8 = aom_highbd_h_predictor_32x8_c;
        if (flags & HAS_AVX2) aom_highbd_h_predictor_32x8 = aom_highbd_h_predictor_32x8_avx2;
        aom_highbd_h_predictor_4x16 = aom_highbd_h_predictor_4x16_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_4x16 = aom_highbd_h_predictor_4x16_sse2;
        aom_highbd_h_predictor_4x4 = aom_highbd_h_predictor_4x4_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_4x4 = aom_highbd_h_predictor_4x4_sse2;
        aom_highbd_h_predictor_4x8 = aom_highbd_h_predictor_4x8_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_4x8 = aom_highbd_h_predictor_4x8_sse2;
        aom_highbd_h_predictor_64x16 = aom_highbd_h_predictor_64x16_c;
        if (flags & HAS_AVX2) aom_highbd_h_predictor_64x16 = aom_highbd_h_predictor_64x16_avx2;
        aom_highbd_h_predictor_64x32 = aom_highbd_h_predictor_64x32_c;
        if (flags & HAS_AVX2) aom_highbd_h_predictor_64x32 = aom_highbd_h_predictor_64x32_avx2;
        aom_highbd_h_predictor_8x32 = aom_highbd_h_predictor_8x32_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_8x32 = aom_highbd_h_predictor_8x32_sse2;
        aom_highbd_h_predictor_64x64 = aom_highbd_h_predictor_64x64_c;
        if (flags & HAS_AVX2) aom_highbd_h_predictor_64x64 = aom_highbd_h_predictor_64x64_avx2;
        aom_highbd_h_predictor_8x16 = aom_highbd_h_predictor_8x16_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_8x16 = aom_highbd_h_predictor_8x16_sse2;
        aom_highbd_h_predictor_8x4 = aom_highbd_h_predictor_8x4_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_8x4 = aom_highbd_h_predictor_8x4_sse2;
        aom_highbd_h_predictor_8x8 = aom_highbd_h_predictor_8x8_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_8x8 = aom_highbd_h_predictor_8x8_sse2;
        aom_highbd_h_predictor_16x16 = aom_highbd_h_predictor_16x16_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_16x16 = aom_highbd_h_predictor_16x16_sse2;
        aom_highbd_h_predictor_16x32 = aom_highbd_h_predictor_16x32_c;
        if (flags & HAS_SSE2) aom_highbd_h_predictor_16x32 = aom_highbd_h_predictor_16x32_sse2;

#endif
        aom_fft2x2_float = aom_fft2x2_float_c;
        aom_fft4x4_float = aom_fft4x4_float_c;
        if (flags & HAS_SSE2) aom_fft4x4_float = aom_fft4x4_float_sse2;
        aom_fft16x16_float = aom_fft16x16_float_c;
        if (flags & HAS_AVX2) aom_fft16x16_float = aom_fft16x16_float_avx2;
        aom_fft32x32_float = aom_fft32x32_float_c;
        if (flags & HAS_AVX2) aom_fft32x32_float = aom_fft32x32_float_avx2;
        aom_fft8x8_float = aom_fft8x8_float_c;
        if (flags & HAS_AVX2) aom_fft8x8_float = aom_fft8x8_float_avx2;

        /*if (flags & HAS_AVX2)*/ aom_ifft16x16_float = aom_ifft16x16_float_avx2;
        /*if (flags & HAS_AVX2)*/ aom_ifft32x32_float = aom_ifft32x32_float_avx2;
        /*if (flags & HAS_AVX2)*/ aom_ifft8x8_float = aom_ifft8x8_float_avx2;
        aom_ifft2x2_float = aom_ifft2x2_float_c;
        /*if (flags & HAS_SSE2)*/ aom_ifft4x4_float = aom_ifft4x4_float_sse2;

    }
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
