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


#include "EbDefinitions.h"
#include "EbAv1Structs.h"

#ifndef EbDeblockingCommon_h
#define EbDeblockingCommon_h
#ifdef __cplusplus
extern "C" {
#endif

#define BLK4X4_ADDR_TO_VERTICAL_EDGE_BS_ARRAY_IDX(blk_4x4_addr) \
    (((blk_4x4_addr) & (MAX_SB_SIZE_IN_4X4BLK - 1)) +           \
     (((blk_4x4_addr) / MAX_SB_SIZE_IN_4X4BLK) * MAX_SB_SIZE_IN_4X4BLK))
#define BLK4X4_ADDR_TO_HORIZONTAL_EDGE_BS_ARRAY_IDX(blk_4x4_addr) \
    (((blk_4x4_addr) & (MAX_SB_SIZE_IN_4X4BLK - 1)) +             \
     (((blk_4x4_addr) / MAX_SB_SIZE_IN_4X4BLK) * MAX_SB_SIZE_IN_4X4BLK))
#define GET_LUMA_4X4BLK_ADDR(                                                       \
    luma_sb_wise4x4_blk_pos_x, luma_sb_wise4x4_blk_pos_y, log_max_sb_size_in4x4blk) \
    (((luma_sb_wise4x4_blk_pos_x) >> 2) +                                           \
     (((luma_sb_wise4x4_blk_pos_y) >> 2) << (log_max_sb_size_in4x4blk)))
#define GET_CHROMA_4X4BLK_ADDR(                                                         \
    chroma_sb_wise2x2_blk_pos_x, chroma_sb_wise2x2_blk_pos_y, log_max_sb_size_in4x4blk) \
    (((chroma_sb_wise2x2_blk_pos_x) >> 1) +                                             \
     (((chroma_sb_wise2x2_blk_pos_y) >> 1) << (log_max_sb_size_in4x4blk)))
#define LUMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(pos_x, pos_y, qp_array_stride) \
    (((pos_x) >> LOG_MIN_BLOCK_SIZE) + ((pos_y) >> LOG_MIN_BLOCK_SIZE) * (qp_array_stride))
#define CHROMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(pos_x, pos_y, qp_array_stride) \
    ((2 * (pos_x) >> LOG_MIN_BLOCK_SIZE) + (2 * (pos_y) >> LOG_MIN_BLOCK_SIZE) * (qp_array_stride))
#define CHECK_MV_COMPONENT_EQUAL_OR_GREATER_THAN_4(pu1Ptr, pu2Ptr, pu1RefList, pu2RefList) \
    (EB_ABS_DIFF((pu1Ptr)->mv[(pu1RefList)].x, (pu2Ptr)->mv[(pu2RefList)].x) >= 4 ||       \
     EB_ABS_DIFF((pu1Ptr)->mv[(pu1RefList)].y, (pu2Ptr)->mv[(pu2RefList)].y) >= 4)

// Precision macros used in the mode decision
#define BIT_ESTIMATE_PRECISION 15
#define LAMBDA_PRECISION 16
#define COST_PRECISION 8
#define MD_SHIFT (BIT_ESTIMATE_PRECISION + LAMBDA_PRECISION - COST_PRECISION)
#define MD_OFFSET (1 << (MD_SHIFT - 1))
#define VAR_QP 1
#define MAX_QP_VALUE_PLUS_INTRA_TC_OFFSET 53
#define BETA_OFFSET_VALUE 12 // range -12 to 12
#define TC_OFFSET_VALUE 12 //12 // range -12 to 12

typedef enum EdgeDir { VERT_EDGE = 0, HORZ_EDGE = 1, NUM_EDGE_DIRS } EdgeDir;

typedef struct Av1DeblockingParameters {
    // length of the filter applied to the outer edge
    uint32_t filter_length;
    // deblocking limits
    const uint8_t *lim;
    const uint8_t *mblim;
    const uint8_t *hev_thr;
} Av1DeblockingParameters;

static const int32_t mode_lf_lut[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // INTRA_MODES
        1, 1, 0, 1, // INTER_MODES (GLOBALMV == 0)
        1, 1, 1, 1, 1, 1, 0, 1 // INTER_COMPOUND_MODES (GLOBAL_GLOBALMV == 0)
};

uint8_t get_filter_level_delta_lf(FrameHeader* frm_hdr, const int32_t dir_idx,
                                  int32_t plane, int32_t *sb_delta_lf, uint8_t seg_id,
                                  PredictionMode pred_mode, MvReferenceFrame ref_frame_0);

static INLINE int32_t is_inter_block_no_intrabc(MvReferenceFrame ref_frame_0) {
    return /*is_intrabc_block(mbmi) ||*/ ref_frame_0 > INTRA_FRAME;
}

void eb_av1_loop_filter_frame_init(FrameHeader *frm_hdr, LoopFilterInfoN *lf_info,
                                   int32_t plane_start, int32_t plane_end);



void update_sharpness(LoopFilterInfoN *lfi, int32_t sharpness_lvl);

#ifdef __cplusplus
}
#endif
#endif // EbDeblockingCommon_h
