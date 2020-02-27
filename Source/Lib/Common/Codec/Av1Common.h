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

#include "EbDefinitions.h"
#include "EbRestoration.h"

#ifndef EbAV1Common_h
#define EbAV1Common_h

#ifdef __cplusplus
extern "C" {
#endif

// Segment Macros
#define SEGMENT_MAX_COUNT 64
#define SEGMENT_COMPLETION_MASK_SET(mask, index)                \
    MULTI_LINE_MACRO_BEGIN(mask) |= (((uint64_t)1) << (index)); \
    MULTI_LINE_MACRO_END
#define SEGMENT_COMPLETION_MASK_TEST(mask, total_count) \
    (mask) == ((((uint64_t)1) << (total_count)) - 1)
#define SEGMENT_ROW_COMPLETION_TEST(mask, row_index, width) \
    ((((mask) >> ((row_index) * (width))) & ((1ull << (width)) - 1)) == ((1ull << (width)) - 1))
#define SEGMENT_CONVERT_IDX_TO_XY(index, x, y, pic_width_in_sb)    \
    MULTI_LINE_MACRO_BEGIN(y) = (index) / (pic_width_in_sb);       \
    (x)                       = (index) - (y) * (pic_width_in_sb); \
    MULTI_LINE_MACRO_END
#define SEGMENT_START_IDX(index, pic_size_in_sb, num_of_seg) \
    (((index) * (pic_size_in_sb)) / (num_of_seg))
#define SEGMENT_END_IDX(index, pic_size_in_sb, num_of_seg) \
    ((((index) + 1) * (pic_size_in_sb)) / (num_of_seg))


typedef struct Av1Common {
    int32_t      mi_rows;
    int32_t      mi_cols;
    int32_t      ref_frame_sign_bias[REF_FRAMES]; /* Two state 0, 1 */
    uint8_t *    last_frame_seg_map;
    InterpFilter interp_filter;
    int32_t      mi_stride;

    // Marks if we need to use 16bit frame buffers (1: yes, 0: no).
    int32_t         use_highbitdepth;
    int32_t         bit_depth;
    int32_t         color_format;
    int32_t         subsampling_x;
    int32_t         subsampling_y;
    RestorationInfo rst_info[MAX_MB_PLANE];
    // rst_end_stripe[i] is one more than the index of the bottom stripe
    // for tile row i.
    int32_t rst_end_stripe[MAX_TILE_ROWS];
    // Output of loop restoration
    Yv12BufferConfig rst_frame;
    // pointer to a scratch buffer used by self-guided restoration
    int32_t *                       rst_tmpbuf;
    Yv12BufferConfig *              frame_to_show;
    int32_t                         byte_alignment;
    int32_t                         last_tile_cols, last_tile_rows;
    int32_t                         log2_tile_cols; // only valid for uniform tiles
    int32_t                         log2_tile_rows; // only valid for uniform tiles
    int32_t                         tile_width, tile_height; // In MI units

//    SeqHeader                       *seq_header_ptr;
    int8_t                          sg_filter_mode;
    int32_t                         sg_frame_ep_cnt[SGRPROJ_PARAMS];
    int32_t                         sg_frame_ep;
    int8_t                          sg_ref_frame_ep[2];
    int8_t                          wn_filter_mode;


    FrameSize frm_size;
    TilesInfo tiles_info;

} Av1Common;

#ifdef __cplusplus
}
#endif
#endif // EbAV1Common_h
