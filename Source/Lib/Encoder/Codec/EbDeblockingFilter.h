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

#include "EbDeblockingFilter_SSE2.h"

#include "EbPredictionUnit.h"
#include "EbNeighborArrays.h"
#include "EbEncDecProcess.h"
#include "EbDlfProcess.h"

#ifndef EbDeblockingFilter_h
#define EbDeblockingFilter_h
#ifdef __cplusplus
extern "C" {
#endif
#define BLK4X4_ADDR_TO_VERTICAL_EDGE_BS_ARRAY_IDX(blk_4x4_addr)                   (((blk_4x4_addr) & (MAX_LCU_SIZE_IN_4X4BLK - 1)) + (((blk_4x4_addr) / MAX_LCU_SIZE_IN_4X4BLK) * MAX_LCU_SIZE_IN_4X4BLK))
#define BLK4X4_ADDR_TO_HORIZONTAL_EDGE_BS_ARRAY_IDX(blk_4x4_addr)                 (((blk_4x4_addr) & (MAX_LCU_SIZE_IN_4X4BLK - 1)) + (((blk_4x4_addr) / MAX_LCU_SIZE_IN_4X4BLK) * MAX_LCU_SIZE_IN_4X4BLK))
#define GET_LUMA_4X4BLK_ADDR(luma_lcu_wise4x4_blk_pos_x, luma_lcu_wise4x4_blk_pos_y, log_max_lcu_size_in4x4blk)           (((luma_lcu_wise4x4_blk_pos_x)>>2) + (((luma_lcu_wise4x4_blk_pos_y)>>2) << (log_max_lcu_size_in4x4blk)))
#define GET_CHROMA_4X4BLK_ADDR(chroma_lcu_wise2x2_blk_pos_x, chroma_lcu_wise2x2_blk_pos_y, log_max_lcu_size_in4x4blk)     (((chroma_lcu_wise2x2_blk_pos_x)>>1) + (((chroma_lcu_wise2x2_blk_pos_y)>>1) << (log_max_lcu_size_in4x4blk)))
#define LUMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(pos_x, pos_y, qp_array_stride)                            (((pos_x) >> LOG_MIN_BLOCK_SIZE) + ((pos_y) >> LOG_MIN_BLOCK_SIZE) * (qp_array_stride))
#define CHROMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(pos_x, pos_y, qp_array_stride)                          ((2*(pos_x) >> LOG_MIN_BLOCK_SIZE) + (2*(pos_y) >> LOG_MIN_BLOCK_SIZE) * (qp_array_stride))
#define CHECK_MV_COMPONENT_EQUAL_OR_GREATER_THAN_4(pu1Ptr, pu2Ptr, pu1RefList, pu2RefList)                    (                         \
                            EB_ABS_DIFF((pu1Ptr)->mv[(pu1RefList)].x, (pu2Ptr)->mv[(pu2RefList)].x) >= 4 ||         \
                            EB_ABS_DIFF((pu1Ptr)->mv[(pu1RefList)].y, (pu2Ptr)->mv[(pu2RefList)].y) >= 4            \
                            )

    // Precision macros used in the mode decision
#define BIT_ESTIMATE_PRECISION              15
#define LAMBDA_PRECISION                    16
#define COST_PRECISION                      8
#define MD_SHIFT                            (BIT_ESTIMATE_PRECISION + LAMBDA_PRECISION - COST_PRECISION)
#define MD_OFFSET                           (1 << (MD_SHIFT-1))
#define VAR_QP                              1
#define MAX_QP_VALUE_PLUS_INTRA_TC_OFFSET   53
#define BETA_OFFSET_VALUE                   12 // range -12 to 12
#define TC_OFFSET_VALUE                     12//12 // range -12 to 12

#if AV1_LF
    typedef enum LpfPickMethod
    {
        // Try the full image with different values.
        LPF_PICK_FROM_FULL_IMAGE,
        // Try a small portion of the image with different values.
        LPF_PICK_FROM_SUBIMAGE,
        // Estimate the level based on quantizer and frame type
        LPF_PICK_FROM_Q,
        // Pick 0 to disable LPF if LPF was enabled last frame
        LPF_PICK_MINIMAL_LPF
    } LpfPickMethod;
#endif

    typedef enum EDGE_DIR { VERT_EDGE = 0, HORZ_EDGE = 1, NUM_EDGE_DIRS } EDGE_DIR;

    typedef struct AV1_DEBLOCKING_PARAMETERS {
        // length of the filter applied to the outer edge
        uint32_t filter_length;
        // deblocking limits
        const uint8_t *lim;
        const uint8_t *mblim;
        const uint8_t *hev_thr;
    } AV1_DEBLOCKING_PARAMETERS;

    void set_qp_array_based_on_cu(
        PictureControlSet *picture_control_set_ptr,          //input parameter
        uint32_t               cuPos_x,                       //input parameter, sample-based horizontal picture-wise locatin of the CU
        uint32_t               cuPos_y,                       //input parameter, sample-based vertical picture-wise locatin of the CU
        uint32_t               cu_size_in_min_cu_size,             //input parameter
        uint32_t               cu_qp);                         //input parameter, Qp of the CU

    /* assorted LoopFilter functions which get used elsewhere */
    struct AV1Common;
    struct macroblockd;
    struct AV1LfSyncData;

    void eb_av1_loop_filter_init(PictureControlSet *pcs_ptr);
    void eb_av1_loop_filter_frame_init(FrameHeader *frm_hdr,
        LoopFilterInfoN *lf_info, int32_t plane_start, int32_t plane_end);

    void loop_filter_sb(
        EbPictureBufferDesc *frame_buffer,//reconpicture,
        //Yv12BufferConfig *frame_buffer,
        PictureControlSet *pcs_ptr,
        MacroBlockD *xd, int32_t mi_row, int32_t mi_col,
        int32_t plane_start, int32_t plane_end,
        uint8_t LastCol);

    void eb_av1_loop_filter_frame(
        EbPictureBufferDesc *frame_buffer,//reconpicture,
        //Yv12BufferConfig *frame_buffer,
        PictureControlSet *pcs_ptr,
        /*MacroBlockD *xd,*/ int32_t plane_start, int32_t plane_end/*,
        int32_t partial_frame*/);

    void eb_av1_pick_filter_level(
        DlfContext            *context_ptr,
        EbPictureBufferDesc   *srcBuffer, // source input
        PictureControlSet     *pcs_ptr,
        LpfPickMethod          method);

    void eb_av1_filter_block_plane_vert(
        const PictureControlSet *const  pcs_ptr,
        const MacroBlockD *const xd,
        const int32_t plane,
        const MacroblockdPlane *const plane_ptr,
        const uint32_t mi_row, const uint32_t mi_col);

    void eb_av1_filter_block_plane_horz(
        const PictureControlSet *const  pcs_ptr,
        const MacroBlockD *const xd, const int32_t plane,
        const MacroblockdPlane *const plane_ptr,
        const uint32_t mi_row, const uint32_t mi_col);

    typedef struct LoopFilterWorkerData
    {
        EbPictureBufferDesc *frame_buffer;//reconpicture,
        PictureControlSet *pcs_ptr;
        struct MacroblockdPlane planes[MAX_MB_PLANE];
        // TODO(Ranjit): When the filter functions are modified to use xd->lossless
        // add lossless as a member here.
        MacroBlockD *xd;
    } LFWorkerData;

    static INLINE int32_t is_inter_block_no_intrabc(MvReferenceFrame ref_frame_0) {
        return /*is_intrabc_block(mbmi) ||*/ ref_frame_0 > INTRA_FRAME;
    }

    void update_sharpness(LoopFilterInfoN *lfi, int32_t sharpness_lvl);

    uint8_t get_filter_level(FrameHeader* frm_hdr, const LoopFilterInfoN *lfi_n,
        const int32_t dir_idx, int32_t plane, int32_t *sb_delta_lf, uint8_t seg_id,
        PredictionMode pred_mode, MvReferenceFrame ref_frame_0);

    void aom_highbd_lpf_horizontal_14_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_horizontal_14 aom_highbd_lpf_horizontal_14_sse2

    void aom_highbd_lpf_horizontal_4_c(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_horizontal_4_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_horizontal_4 aom_highbd_lpf_horizontal_4_sse2

    void aom_highbd_lpf_horizontal_6_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_horizontal_6 aom_highbd_lpf_horizontal_6_sse2

    void aom_highbd_lpf_horizontal_8_c(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_horizontal_8_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_horizontal_8 aom_highbd_lpf_horizontal_8_sse2

    void aom_highbd_lpf_vertical_14_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_vertical_14 aom_highbd_lpf_vertical_14_sse2

    void aom_highbd_lpf_vertical_4_c(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_vertical_4_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_vertical_4 aom_highbd_lpf_vertical_4_sse2

    void aom_highbd_lpf_vertical_6_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_vertical_6 aom_highbd_lpf_vertical_6_sse2

    void aom_highbd_lpf_vertical_8_c(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
    void aom_highbd_lpf_vertical_8_sse2(uint16_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int32_t bd);
#define aom_highbd_lpf_vertical_8 aom_highbd_lpf_vertical_8_sse2

    void aom_lpf_horizontal_14_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_14 aom_lpf_horizontal_14_sse2

    void aom_lpf_horizontal_4_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_horizontal_4_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_4 aom_lpf_horizontal_4_sse2

    void aom_lpf_horizontal_6_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_horizontal_6_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_6 aom_lpf_horizontal_6_sse2

    void aom_lpf_horizontal_8_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_horizontal_8_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_horizontal_8 aom_lpf_horizontal_8_sse2

    void aom_lpf_vertical_14_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_14 aom_lpf_vertical_14_sse2

    void aom_lpf_vertical_4_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_vertical_4_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_4 aom_lpf_vertical_4_sse2

    void aom_lpf_vertical_6_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_6 aom_lpf_vertical_6_sse2

    void aom_lpf_vertical_8_c(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
    void aom_lpf_vertical_8_sse2(uint8_t *s, int32_t pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh);
#define aom_lpf_vertical_8 aom_lpf_vertical_8_sse2

#ifdef __cplusplus
}
#endif
#endif
// clang-format on
