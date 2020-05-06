/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
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
#ifndef EbTemporalFiltering_h
#define EbTemporalFiltering_h

#include "EbPictureControlSet.h"
#include "EbMotionEstimationProcess.h"

// ALT-REF debug-specific defines
#define DEBUG_TF 0

#define COLOR_CHANNELS 3
#define C_Y 0
#define C_U 1
#define C_V 2

#define EDGE_THRESHOLD 50
#define SQRT_PI_BY_2 1.25331413732
#define SMOOTH_THRESHOLD 16
// Block size used in temporal filtering
#define BW 64
#define BH 64
#define BLK_PELS 4096 // Pixels in the block
#define TF_ENABLE_PLANEWISE_STRATEGY 1
// Window size for plane-wise temporal filtering.
// This is particually used for function `av1_apply_temporal_filter_planewise()`
#define TF_PLANEWISE_FILTER_WINDOW_LENGTH 5
// A scale factor used in plane-wise temporal filtering to raise the filter
// weight from `double` with range [0, 1] to `int` with range [0, 1000].
#define TF_PLANEWISE_FILTER_WEIGHT_SCALE 1000
#define N_16X16_BLOCKS 16
#define N_32X32_BLOCKS 4

#define INT_MAX_TF 2147483647 //max value for an int
#define INT_MIN_TF (-2147483647 - 1) //min value for an int
#define THR_SHIFT 2 // should be 2

#define INIT_WEIGHT 2
#define WEIGHT_MULTIPLIER 16

#define THRES_LOW 10000
#define THRES_HIGH 20000
#define THRES_DIFF_LOW 6000
#define THRES_DIFF_HIGH 12000

#define OD_DIVU_DMAX (1024)
#define AHD_TH_WEIGHT 33
#ifdef __cplusplus
extern "C" {
#endif

int svt_av1_init_temporal_filtering(PictureParentControlSet ** list_picture_control_set_ptr,
                                    PictureParentControlSet *  picture_control_set_ptr_central,
                                    MotionEstimationContext_t *me_context_ptr,
                                    int32_t                    segment_index);

void svt_av1_apply_filtering_c(const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
                               int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src,
                               int uv_src_stride, const uint8_t *u_pre, const uint8_t *v_pre,
                               int uv_pre_stride, unsigned int block_width,
                               unsigned int block_height, int ss_x, int ss_y, int strength,
                               const int *blk_fw, int use_whole_blk, uint32_t *y_accum,
                               uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count,
                               uint32_t *v_accum, uint16_t *v_count);

void svt_av1_apply_filtering_highbd_c(
    const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre, int y_pre_stride,
    const uint16_t *u_src, const uint16_t *v_src, int uv_src_stride, const uint16_t *u_pre,
    const uint16_t *v_pre, int uv_pre_stride, unsigned int block_width, unsigned int block_height,
    int ss_x, int ss_y, int strength, const int *blk_fw, int use_whole_blk, uint32_t *y_accum,
    uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count);
void svt_av1_apply_temporal_filter_planewise_c(
    const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre, int y_pre_stride,
    const uint8_t *u_src, const uint8_t *v_src, int uv_src_stride, const uint8_t *u_pre,
    const uint8_t *v_pre, int uv_pre_stride, unsigned int block_width, unsigned int block_height,
    int ss_x, int ss_y, const double *noise_levels, const int decay_control, uint32_t *y_accum,
    uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count);
#ifdef __cplusplus
}
#endif
#endif //EbTemporalFiltering_h
