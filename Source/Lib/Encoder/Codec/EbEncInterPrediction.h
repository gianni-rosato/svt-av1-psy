/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbEncInterPrediction_h
#define EbEncInterPrediction_h

#include "EbModeDecision.h"
#include "EbMcp.h"
#include "filter.h"
#include "convolve.h"
#include "EbInterPrediction.h"
#include "EbEncIntraPrediction.h"

#ifdef __cplusplus
extern "C" {
#endif

extern aom_highbd_convolve_fn_t convolveHbd[/*subX*/ 2][/*subY*/ 2][/*bi*/ 2];

void svt_av1_init_inter_params(InterPredParams *inter_pred_params, int block_width,
                       int block_height, int pix_row, int pix_col,
                       int subsampling_x, int subsampling_y, int bit_depth,
                       int use_hbd_buf, int is_intrabc,
                       const struct ScaleFactors *sf,
                       const struct Buf2D *ref_buf,
                       uint32_t interp_filters);
void av1_build_inter_predictor(Av1Common *cm, const uint8_t *src, int src_stride, uint8_t *dst,
                               int dst_stride, const MV *src_mv, int pix_col, int pix_row,
                               InterPredParams *inter_pred_params);
EbErrorType av1_inter_prediction(
        PictureControlSet              *pcs_ptr,
        uint32_t                        interp_filters,
        BlkStruct                      *blk_ptr,
        uint8_t                         ref_frame_type,
        MvUnit                         *mv_unit,
        uint8_t                         use_intrabc,
        MotionMode                      motion_mode,
        uint8_t                         use_precomputed_obmc,
        struct ModeDecisionContext     *md_context,
        uint8_t                         compound_idx,
        InterInterCompoundData         *interinter_comp,
        TileInfo                       * tile,
        NeighborArrayUnit              *luma_recon_neighbor_array,
        NeighborArrayUnit              *cb_recon_neighbor_array ,
        NeighborArrayUnit              *cr_recon_neighbor_array ,
        uint8_t                         is_interintra_used ,
        InterIntraMode                 interintra_mode,
        uint8_t                         use_wedge_interintra,
        int32_t                         interintra_wedge_index,
        uint16_t                        pu_origin_x,
        uint16_t                        pu_origin_y,
        uint8_t                         bwidth,
        uint8_t                         bheight,
        EbPictureBufferDesc             *ref_pic_list0,
        EbPictureBufferDesc             *ref_pic_list1,
        EbPictureBufferDesc             *prediction_ptr,
        uint16_t                        dst_origin_x,
        uint16_t                        dst_origin_y,
        EbBool                          perform_chroma,
        uint8_t                         bit_depth);

EbErrorType av1_inter_prediction_16bit_pipeline(
    PictureControlSet              *pcs_ptr,
    uint32_t                        interp_filters,
    BlkStruct                      *blk_ptr,
    uint8_t                         ref_frame_type,
    MvUnit                         *mv_unit,
    uint8_t                         use_intrabc,
    MotionMode                      motion_mode,
    uint8_t                         use_precomputed_obmc,
    struct ModeDecisionContext     *md_context,
    uint8_t                         compound_idx,
    InterInterCompoundData         *interinter_comp,
    TileInfo                       * tile,
    NeighborArrayUnit              *luma_recon_neighbor_array,
    NeighborArrayUnit              *cb_recon_neighbor_array,
    NeighborArrayUnit              *cr_recon_neighbor_array,
    uint8_t                         is_interintra_used,
    InterIntraMode                 interintra_mode,
    uint8_t                         use_wedge_interintra,
    int32_t                         interintra_wedge_index,
    uint16_t                        pu_origin_x,
    uint16_t                        pu_origin_y,
    uint8_t                         bwidth,
    uint8_t                         bheight,
    EbPictureBufferDesc             *ref_pic_list0,
    EbPictureBufferDesc             *ref_pic_list1,
    EbPictureBufferDesc             *prediction_ptr,
    uint16_t                        dst_origin_x,
    uint16_t                        dst_origin_y,
    EbBool                          perform_chroma,
    uint8_t                         bit_depth);

void search_compound_diff_wedge(
        PictureControlSet                    *pcs_ptr,
        struct ModeDecisionContext                  *context_ptr,
        ModeDecisionCandidate                *candidate_ptr);
void calc_pred_masked_compound(
        PictureControlSet                    *pcs_ptr,
        struct ModeDecisionContext           *context_ptr,
        ModeDecisionCandidate                *candidate_ptr);


EbErrorType inter_pu_prediction_av1(
        uint8_t                              hbd_mode_decision,
        struct ModeDecisionContext           *md_context_ptr,
        PictureControlSet                    *pcs_ptr,
        ModeDecisionCandidateBuffer          *candidate_buffer_ptr);

EbErrorType warped_motion_prediction(
        PictureControlSet                    *pcs_ptr,
        MvUnit                               *mv_unit,
        uint8_t                               ref_frame_type,
        uint8_t                               compound_idx,
        InterInterCompoundData               *interinter_comp,
        uint16_t                              pu_origin_x,
        uint16_t                              pu_origin_y,
        BlkStruct                           *blk_ptr,
        const BlockGeom                      *blk_geom,
        EbPictureBufferDesc                  *ref_pic_list0,
        EbPictureBufferDesc                  *ref_pic_list1,
        EbPictureBufferDesc                  *prediction_ptr,
        uint16_t                              dst_origin_x,
        uint16_t                              dst_origin_y,
        EbWarpedMotionParams                 *wm_params_l0,
        EbWarpedMotionParams                 *wm_params_l1,
        uint8_t                               bit_depth,
        EbBool                                perform_chroma,
        EbBool                                is_encode_pass);

const uint8_t *svt_av1_get_obmc_mask(int length);

int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);

void model_rd_from_sse(BlockSize bsize, int16_t quantizer, uint8_t bit_depth, uint64_t sse,
                       uint32_t *rate, uint64_t *dist);

#ifdef __cplusplus
}
#endif
#endif //EbEncInterPrediction_h
