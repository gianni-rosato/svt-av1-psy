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

#include <stdlib.h>
#include <stdio.h>

void svt_av1_apply_temporal_filter_sse4_1(
        const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
        int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src,
        int uv_src_stride, const uint8_t *u_pre, const uint8_t *v_pre,
        int uv_pre_stride, unsigned int block_width, unsigned int block_height,
        int ss_x, int ss_y, int strength, const int *blk_fw, int use_whole_blk,
        uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count,
        uint32_t *v_accum, uint16_t *v_count);
void svt_av1_highbd_apply_temporal_filter_sse4_1(
        const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
        int y_pre_stride, const uint16_t *u_src, const uint16_t *v_src,
        int uv_src_stride, const uint16_t *u_pre, const uint16_t *v_pre,
        int uv_pre_stride, unsigned int block_width, unsigned int block_height,
        int ss_x, int ss_y, int strength, const int *blk_fw, int use_whole_blk,
        uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count,
        uint32_t *v_accum, uint16_t *v_count);
