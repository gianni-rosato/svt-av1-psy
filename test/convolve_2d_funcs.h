/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */
/******************************************************************************
 * @file Convolve2dTest.cc
 *
 * @brief Function declaration for interpolation in inter prediction:
 * - av1_highbd_convolve_2d_copy_sr_avx2
 * - av1_highbd_jnt_convolve_2d_copy_avx2
 * - av1_highbd_convolve_x_sr_avx2
 * - av1_highbd_convolve_y_sr_avx2
 * - av1_highbd_convolve_2d_sr_avx2
 * - av1_highbd_jnt_convolve_x_avx2
 * - av1_highbd_jnt_convolve_y_avx2
 * - av1_highbd_jnt_convolve_2d_avx2
 * @author Cidana-Wenyao
 *
 ******************************************************************************/
#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
InterpFilterParams av1_get_interp_filter_params_with_block_size(
    const InterpFilter interp_filter, const int32_t w);
// higbbd sr version
void av1_highbd_convolve_x_sr_avx2(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
void av1_highbd_convolve_y_sr_avx2(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
void av1_highbd_convolve_2d_sr_avx2(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
void av1_highbd_convolve_2d_copy_sr_avx2(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
void av1_highbd_convolve_2d_sr_c(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
// highbd jnt version
void av1_highbd_jnt_convolve_x_avx2(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
void av1_highbd_jnt_convolve_y_avx2(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
void av1_highbd_jnt_convolve_2d_avx2(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
void av1_highbd_jnt_convolve_2d_copy_avx2(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
void av1_highbd_jnt_convolve_2d_c(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);
// Lowbd sr version
void av1_convolve_2d_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst,
                          int32_t dst_stride, int32_t w, int32_t h,
                          InterpFilterParams *filter_params_x,
                          InterpFilterParams *filter_params_y,
                          const int32_t subpel_x_q4, const int32_t subpel_y_q4,
                          ConvolveParams *conv_params);
void av1_convolve_x_sr_avx2(const uint8_t *src, int32_t src_stride,
                            uint8_t *dst, int32_t dst_stride, int32_t w,
                            int32_t h, InterpFilterParams *filter_params_x,
                            InterpFilterParams *filter_params_y,
                            const int32_t subpel_x_q4,
                            const int32_t subpel_y_q4,
                            ConvolveParams *conv_params);
void av1_convolve_y_sr_avx2(const uint8_t *src, int32_t src_stride,
                            uint8_t *dst, int32_t dst_stride, int32_t w,
                            int32_t h, InterpFilterParams *filter_params_x,
                            InterpFilterParams *filter_params_y,
                            const int32_t subpel_x_q4,
                            const int32_t subpel_y_q4,
                            ConvolveParams *conv_params);
void av1_convolve_2d_copy_sr_avx2(
    const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params);
void av1_convolve_2d_sr_avx2(const uint8_t *src, int32_t src_stride,
                             uint8_t *dst, int32_t dst_stride, int32_t w,
                             int32_t h, InterpFilterParams *filter_params_x,
                             InterpFilterParams *filter_params_y,
                             const int32_t subpel_x_q4,
                             const int32_t subpel_y_q4,
                             ConvolveParams *conv_params);
// Lowbd jnt version
void av1_jnt_convolve_x_avx2(const uint8_t *src, int32_t src_stride,
                             uint8_t *dst, int32_t dst_stride, int32_t w,
                             int32_t h, InterpFilterParams *filter_params_x,
                             InterpFilterParams *filter_params_y,
                             const int32_t subpel_x_q4,
                             const int32_t subpel_y_q4,
                             ConvolveParams *conv_params);
void av1_jnt_convolve_y_avx2(const uint8_t *src, int32_t src_stride,
                             uint8_t *dst, int32_t dst_stride, int32_t w,
                             int32_t h, InterpFilterParams *filter_params_x,
                             InterpFilterParams *filter_params_y,
                             const int32_t subpel_x_q4,
                             const int32_t subpel_y_q4,
                             ConvolveParams *conv_params);
void av1_jnt_convolve_2d_avx2(const uint8_t *src, int32_t src_stride,
                              uint8_t *dst, int32_t dst_stride, int32_t w,
                              int32_t h, InterpFilterParams *filter_params_x,
                              InterpFilterParams *filter_params_y,
                              const int32_t subpel_x_q4,
                              const int32_t subpel_y_q4,
                              ConvolveParams *conv_params);
void av1_jnt_convolve_2d_copy_avx2(
    const uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride,
    int32_t w, int32_t h, InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params);
void av1_jnt_convolve_2d_c(const uint8_t *src, int32_t src_stride, uint8_t *dst,
                           int32_t dst_stride, int32_t w, int32_t h,
                           InterpFilterParams *filter_params_x,
                           InterpFilterParams *filter_params_y,
                           const int32_t subpel_x_q4, const int32_t subpel_y_q4,
                           ConvolveParams *conv_params);
#ifdef __cplusplus
}
#endif
