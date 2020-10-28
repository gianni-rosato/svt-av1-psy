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

#ifndef EbEncIntraPrediction_h
#define EbEncIntraPrediction_h

#include "EbMotionEstimationProcess.h"
#include "EbModeDecision.h"
#include "EbIntraPrediction.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ModeDecisionCandidateBuffer;
struct ModeDecisionCandidate;

/////////..............................................//////////////////////////

extern EbErrorType svt_av1_intra_prediction_cl(uint8_t                      hbd_mode_decision,
                                               struct ModeDecisionContext * context_ptr,
                                               PictureControlSet *          pcs_ptr,
                                               ModeDecisionCandidateBuffer *candidate_buffer_ptr);
extern EbErrorType update_neighbor_samples_array_open_loop_mb(uint8_t *above_ref, uint8_t *left_ref,
                                                           EbPictureBufferDesc *input_ptr,
                                                           uint32_t stride, uint32_t srcOriginX,
                                                           uint32_t srcOriginY, uint8_t bwidth,
                                                           uint8_t bheight);
extern EbErrorType update_neighbor_samples_array_open_loop_mb_recon(
    uint8_t *above_ref, uint8_t *left_ref, uint8_t *recon_ptr, uint32_t stride,
    uint32_t src_origin_x, uint32_t src_origin_y, uint8_t bwidth, uint8_t bheight, uint32_t width,
    uint32_t height);



void svt_av1_predict_intra_block(TileInfo *tile, STAGE stage, const BlockGeom *blk_geom,
                                 const Av1Common *cm, int32_t wpx, int32_t hpx, TxSize tx_size,
                                 PredictionMode mode, int32_t angle_delta, int32_t use_palette,
                                 PaletteInfo *palette_info, FilterIntraMode filter_intra_mode,
                                 uint8_t *top_neigh_array, uint8_t *left_neigh_array,
                                 EbPictureBufferDesc *recon_buffer, int32_t col_off, int32_t row_off,
                                 int32_t plane, BlockSize bsize, uint32_t txb_org_x_pict,
                                 uint32_t txb_org_y_pict, uint32_t bl_org_x_pict,
                                 uint32_t bl_org_y_pict, uint32_t bl_org_x_mb, uint32_t bl_org_y_mb,
                                 ModeInfo **mi_grid_base, SeqHeader *seq_header_ptr);

void svt_av1_predict_intra_block_16bit(
        EbBitDepthEnum bit_depth, TileInfo *tile, STAGE stage, const BlockGeom *blk_geom, const Av1Common *cm, int32_t wpx,
        int32_t hpx, TxSize tx_size, PredictionMode mode, int32_t angle_delta, int32_t use_palette,
        PaletteInfo *palette_info, FilterIntraMode filter_intra_mode, uint16_t *top_neigh_array,
        uint16_t *left_neigh_array, EbPictureBufferDesc *recon_buffer, int32_t col_off, int32_t row_off,
        int32_t plane, BlockSize bsize, uint32_t txb_org_x_pict, uint32_t txb_org_y_pict,
        uint32_t bl_org_x_pict, uint32_t bl_org_y_pict, uint32_t bl_org_x_mb, uint32_t bl_org_y_mb,
        ModeInfo **mi_grid_base, SeqHeader *seq_header_ptr);


#ifdef __cplusplus
}
#endif
#endif // EbEncIntraPrediction_h
