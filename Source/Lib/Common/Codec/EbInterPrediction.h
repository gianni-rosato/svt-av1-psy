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
#include "EbObject.h"

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct SubpelParams {
        int32_t xs;
        int32_t ys;
        int32_t subpel_x;
        int32_t subpel_y;
    } SubpelParams;

    struct ModeDecisionContext;

    typedef struct InterPredictionContext {
        EbDctor                               dctor;
        MotionCompensationPredictionContext  *mcp_context;
    } InterPredictionContext;

    void svt_inter_predictor(const uint8_t *src, int32_t src_stride,
        uint8_t *dst, int32_t dst_stride, const SubpelParams *subpel_params,
        const ScaleFactors *sf, int32_t w, int32_t h, ConvolveParams *conv_params,
        InterpFilters interp_filters, int32_t is_intrabc);

    void svt_highbd_inter_predictor(const uint16_t *src, int32_t src_stride,
        uint16_t *dst, int32_t dst_stride, const SubpelParams *subpel_params,
        const ScaleFactors *sf, int32_t w, int32_t h, ConvolveParams *conv_params,
        InterpFilters interp_filters, int32_t is_intrabc, int32_t bd);

    EbErrorType av1_inter_prediction(
        PictureControlSet                    *picture_control_set_ptr,
        uint32_t                                interp_filters,
        CodingUnit                           *cu_ptr,
        uint8_t                                 ref_frame_type,
        MvUnit                               *mv_unit,
        uint8_t                                  use_intrabc,
        uint16_t                                pu_origin_x,
        uint16_t                                pu_origin_y,
        uint8_t                                 bwidth,
        uint8_t                                 bheight,
        EbPictureBufferDesc                  *ref_pic_list0,
        EbPictureBufferDesc                  *ref_pic_list1,
        EbPictureBufferDesc                  *prediction_ptr,
        uint16_t                                dst_origin_x,
        uint16_t                                dst_origin_y,
        EbBool                                  perform_chroma,
        EbAsm                                   asm_type);

    EbErrorType inter_pu_prediction_av1(
        struct ModeDecisionContext           *md_context_ptr,
        PictureControlSet                    *picture_control_set_ptr,
        ModeDecisionCandidateBuffer          *candidate_buffer_ptr,
        EbAsm                                   asm_type);

    EbErrorType av1_inter_prediction_hbd(
        PictureControlSet                    *picture_control_set_ptr,
        uint8_t                                 ref_frame_type,
        CodingUnit                           *cu_ptr,
        MvUnit                               *mv_unit,
        uint8_t                                  use_intrabc,
        uint16_t                                pu_origin_x,
        uint16_t                                pu_origin_y,
        uint8_t                                 bwidth,
        uint8_t                                 bheight,
        EbPictureBufferDesc                  *ref_pic_list0,
        EbPictureBufferDesc                  *ref_pic_list1,
        EbPictureBufferDesc                  *prediction_ptr,
        uint16_t                                dst_origin_x,
        uint16_t                                dst_origin_y,
        uint8_t                                 bit_depth,
        EbAsm                                   asm_type);

    EbErrorType choose_mvp_idx_v2(
        ModeDecisionCandidate               *candidate_ptr,
        uint32_t                               cu_origin_x,
        uint32_t                               cu_origin_y,
        uint32_t                               pu_index,
        uint32_t                               tb_size,
        int16_t                               *ref0_amvp_cand_array_x,
        int16_t                               *ref0_amvp_cand_array_y,
        uint32_t                               ref0_num_available_amvp_cand,
        int16_t                               *ref1_amvp_cand_array_x,
        int16_t                               *ref1_amvp_cand_array_y,
        uint32_t                               ref1_num_available_amvp_cand,
        PictureControlSet                   *picture_control_set_ptr);

    EbErrorType warped_motion_prediction(
        MvUnit                               *mv_unit,
        uint16_t                                pu_origin_x,
        uint16_t                                pu_origin_y,
        CodingUnit                           *cu_ptr,
        const BlockGeom                        *blk_geom,
        EbPictureBufferDesc                  *ref_pic_list0,
        EbPictureBufferDesc                  *prediction_ptr,
        uint16_t                                dst_origin_x,
        uint16_t                                dst_origin_y,
        EbWarpedMotionParams                   *wm_params,
        uint8_t                                 bit_depth,
        EbBool                                  perform_chroma,
        EbAsm                                   asm_type);

#ifdef __cplusplus
}
#endif
#endif //EbInterPrediction_h
