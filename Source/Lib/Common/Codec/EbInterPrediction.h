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

    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_8[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_4[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_8sharp[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_8smooth[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, bilinear_filters[SUBPEL_SHIFTS]);
    extern DECLARE_ALIGNED(256, const InterpKernel, sub_pel_filters_4smooth[SUBPEL_SHIFTS]);

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
        uint8_t                                compound_idx,
    InterInterCompoundData                     *interinter_comp,
#if II_COMP_FLAG
        TileInfo                                * tile,
        NeighborArrayUnit                       *luma_recon_neighbor_array,
        NeighborArrayUnit                       *cb_recon_neighbor_array ,
        NeighborArrayUnit                       *cr_recon_neighbor_array ,
        uint8_t                                 is_interintra_used ,
        INTERINTRA_MODE                         interintra_mode,
        uint8_t                                 use_wedge_interintra,
        int32_t                                 interintra_wedge_index,
#endif
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
    void search_compound_diff_wedge(
        PictureControlSet                    *picture_control_set_ptr,
        struct ModeDecisionContext                  *context_ptr,
        ModeDecisionCandidate                *candidate_ptr);
    void search_compound_avg_dist(
        PictureControlSet                    *picture_control_set_ptr,
        struct ModeDecisionContext                    *context_ptr,
        ModeDecisionCandidate                *candidate_ptr);

    void av1_dist_wtd_comp_weight_assign(
        SeqHeader *seq_header,
        int cur_frame_index,
        int bck_frame_index,
        int fwd_frame_index,
        int compound_idx,
        int order_idx,
        int *fwd_offset, int *bck_offset,
        int *use_dist_wtd_comp_avg,
        int is_compound);

    void build_masked_compound_no_round(uint8_t *dst, int dst_stride,
        const CONV_BUF_TYPE *src0, int src0_stride,
        const CONV_BUF_TYPE *src1, int src1_stride,
        const InterInterCompoundData *const comp_data, uint8_t *seg_mask,
        BlockSize sb_type, int h, int w, ConvolveParams *conv_params,
        uint8_t bd);

    void av1_get_convolve_filter_params(uint32_t interp_filters,
        InterpFilterParams *params_x, InterpFilterParams *params_y,
        int32_t w, int32_t h);

#if COMP_INTERINTRA
    /* Mapping of interintra to intra mode for use in the intra component */
    static const PredictionMode interintra_to_intra_mode[INTERINTRA_MODES] = {
      DC_PRED, V_PRED, H_PRED, SMOOTH_PRED
    };
    static INLINE int is_interintra_wedge_used(BlockSize sb_type) {
        return wedge_params_lookup[sb_type].bits > 0;
    }

    void combine_interintra(INTERINTRA_MODE mode,
        int8_t use_wedge_interintra, int wedge_index,
        int wedge_sign, BlockSize bsize,
        BlockSize plane_bsize, uint8_t *comppred,
        int compstride, const uint8_t *interpred,
        int interstride, const uint8_t *intrapred,
        int intrastride);

    void combine_interintra_highbd(
        InterIntraMode mode, uint8_t use_wedge_interintra, uint8_t wedge_index,
        uint8_t wedge_sign, BlockSize bsize, BlockSize plane_bsize,
        uint8_t *comppred8, int compstride, const uint8_t *interpred8,
        int interstride, const uint8_t *intrapred8, int intrastride, int bd);

#endif //comp_interintra

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

    extern aom_highbd_convolve_fn_t convolveHbd[/*subX*/2][/*subY*/2][/*bi*/2];
    extern aom_convolve_fn_t convolve[/*subX*/2][/*subY*/2][/*bi*/2];

#ifdef __cplusplus
}
#endif
#endif //EbInterPrediction_h
