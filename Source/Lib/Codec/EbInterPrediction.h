/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbInterPrediction_h
#define EbInterPrediction_h

#include "EbDefinitions.h"
#include "EbPictureControlSet.h"
#include "EbCodingUnit.h"
#include "EbPredictionUnit.h"
#include "EbModeDecision.h"
#include "EbMcp.h"
#include "EbMvMerge.h"

#ifdef __cplusplus
extern "C" {
#endif

    struct ModeDecisionContext_s;
    typedef struct InterPredictionContext_s {
        // mcp context
        MotionCompensationPredictionContext_t  *mcp_context;
    } InterPredictionContext_t;

    extern EbErrorType inter_prediction_context_ctor(
        InterPredictionContext_t   **inter_prediction_context,
        uint16_t                     max_cu_width,
        uint16_t                     max_cu_height);

    extern EbErrorType inter2_nx2_n_pu_prediction_avc(
        struct ModeDecisionContext_s           *context_ptr,
        uint32_t                                component_mask,
        PictureControlSet_t                    *picture_control_set_ptr,
        ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
        EbAsm                                   asm_type);

    EbErrorType inter2_nx2_n_pu_prediction_avc_style(
        struct ModeDecisionContext_s           *context_ptr,
        uint32_t                                component_mask,
        PictureControlSet_t                    *picture_control_set_ptr,
        ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
        EbAsm                                   asm_type);

    EbErrorType av1_inter_prediction(
        PictureControlSet_t                    *picture_control_set_ptr,
        uint32_t                                interp_filters,
        CodingUnit_t                           *cu_ptr,
        uint8_t                                 ref_frame_type,
        MvUnit_t                               *mv_unit,
        uint16_t                                pu_origin_x,
        uint16_t                                pu_origin_y,
        uint8_t                                 bwidth,
        uint8_t                                 bheight,
        EbPictureBufferDesc_t                  *ref_pic_list0,
        EbPictureBufferDesc_t                  *ref_pic_list1,
        EbPictureBufferDesc_t                  *prediction_ptr,
        uint16_t                                dst_origin_x,
        uint16_t                                dst_origin_y,
        EbAsm                                   asm_type);

    EbErrorType inter_pu_prediction_av1(
        struct ModeDecisionContext_s           *md_context_ptr,
        uint32_t                                component_mask,
        PictureControlSet_t                    *picture_control_set_ptr,
        ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
        EbAsm                                   asm_type);

    EbErrorType av1_inter_prediction_hbd(
        PictureControlSet_t                    *picture_control_set_ptr,
        uint8_t                                 ref_frame_type,
        CodingUnit_t                           *cu_ptr,
        MvUnit_t                               *mv_unit,
        uint16_t                                pu_origin_x,
        uint16_t                                pu_origin_y,
        uint8_t                                 bwidth,
        uint8_t                                 bheight,
        EbPictureBufferDesc_t                  *ref_pic_list0,
        EbPictureBufferDesc_t                  *ref_pic_list1,
        EbPictureBufferDesc_t                  *prediction_ptr,
        uint16_t                                dst_origin_x,
        uint16_t                                dst_origin_y,
        uint8_t                                 bit_depth,
        EbAsm                                   asm_type);

    EbErrorType choose_mvp_idx_v2(
        ModeDecisionCandidate_t               *candidate_ptr,
        uint32_t                               cu_origin_x,
        uint32_t                               cu_origin_y,
        uint32_t                               pu_index,
        uint32_t                               tbSize,
        int16_t                               *ref0_amvp_cand_array_x,
        int16_t                               *ref0_amvp_cand_array_y,
        uint32_t                               ref0_num_available_amvp_cand,
        int16_t                               *ref1_amvp_cand_array_x,
        int16_t                               *ref1_amvp_cand_array_y,
        uint32_t                               ref1_num_available_amvp_cand,
        PictureControlSet_t                   *picture_control_set_ptr);


    EbErrorType warped_motion_prediction(
        MvUnit_t                               *mv_unit,
        uint16_t                                pu_origin_x,
        uint16_t                                pu_origin_y,
        CodingUnit_t                           *cu_ptr,
        const BlockGeom                        *blk_geom,
        EbPictureBufferDesc_t                  *ref_pic_list0,
        EbPictureBufferDesc_t                  *prediction_ptr,
        uint16_t                                dst_origin_x,
        uint16_t                                dst_origin_y,
        EbWarpedMotionParams                   *wm_params,
        uint8_t                                 bit_depth,
        EbAsm                                   asm_type);

#ifdef __cplusplus
}
#endif
#endif //EbInterPrediction_h