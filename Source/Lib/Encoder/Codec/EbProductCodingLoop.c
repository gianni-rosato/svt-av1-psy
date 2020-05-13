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

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbTransformUnit.h"
#include "EbRateDistortionCost.h"
#include "EbFullLoop.h"
#include "EbPictureOperators.h"
#include "EbModeDecisionProcess.h"
#include "EbTransforms.h"
#include "EbMotionEstimation.h"
#include "aom_dsp_rtcd.h"
#include "EbCodingLoop.h"
#include "EbLog.h"
#include "EbCommonUtils.h"
#include "EbResize.h"

#define DIVIDE_AND_ROUND(x, y) (((x) + ((y) >> 1)) / (y))

EbErrorType generate_md_stage_0_cand(SuperBlock *sb_ptr, ModeDecisionContext *context_ptr,
                                     uint32_t *         fast_candidate_total_count,
                                     PictureControlSet *pcs_ptr);

int16_t eb_av1_dc_quant_qtx(int32_t qindex, int32_t delta, AomBitDepth bit_depth);

static INLINE int is_interintra_allowed_bsize(const BlockSize bsize) {
    return (bsize >= BLOCK_8X8) && (bsize <= BLOCK_32X32);
}
void precompute_intra_pred_for_inter_intra(PictureControlSet *  pcs_ptr,
                                           ModeDecisionContext *context_ptr);

int svt_av1_allow_palette(int allow_palette, BlockSize sb_type);

/*******************************************
* set Penalize Skip Flag
*
* Summary: Set the penalize_skipflag to true
* When there is luminance/chrominance change
* or in noisy clip with low motion at meduim
* varince area
*
*******************************************/

const EbPredictionFunc product_prediction_fun_table[3] = {
    NULL, inter_pu_prediction_av1, eb_av1_intra_prediction_cl};

const EbFastCostFunc av1_product_fast_cost_func_table[3] = {
    NULL,
    av1_inter_fast_cost, /*INTER */
    av1_intra_fast_cost /*INTRA */
};

const EbAv1FullCostFunc av1_product_full_cost_func_table[3] = {
    NULL,
    av1_inter_full_cost, /*INTER */
    av1_intra_full_cost /*INTRA */
};

/***************************************************
* Update Recon Samples Neighbor Arrays
***************************************************/
void mode_decision_update_neighbor_arrays(PictureControlSet *  pcs_ptr,
                                          ModeDecisionContext *context_ptr, uint32_t index_mds,
                                          EbBool intra4x4Selected) {
    uint32_t bwdith  = context_ptr->blk_geom->bwidth;
    uint32_t bheight = context_ptr->blk_geom->bheight;

    uint32_t origin_x = context_ptr->blk_origin_x;
    uint32_t origin_y = context_ptr->blk_origin_y;
    (void)intra4x4Selected;

    uint32_t blk_origin_x_uv = context_ptr->round_origin_x >> 1;
    uint32_t blk_origin_y_uv = context_ptr->round_origin_y >> 1;
    uint32_t bwdith_uv       = context_ptr->blk_geom->bwidth_uv;
    uint32_t bwheight_uv     = context_ptr->blk_geom->bheight_uv;

    uint8_t mode_type       = context_ptr->blk_ptr->prediction_mode_flag;
    uint8_t intra_luma_mode = (uint8_t)context_ptr->blk_ptr->pred_mode;
    uint8_t chroma_mode = (uint8_t)context_ptr->blk_ptr->prediction_unit_array->intra_chroma_mode;
    uint8_t skip_flag   = (uint8_t)context_ptr->blk_ptr->skip_flag;

    context_ptr->mv_unit.pred_direction = (uint8_t)(
        context_ptr->md_blk_arr_nsq[index_mds].prediction_unit_array[0].inter_pred_direction_index);
    context_ptr->mv_unit.mv[REF_LIST_0].mv_union =
        context_ptr->md_blk_arr_nsq[index_mds].prediction_unit_array[0].mv[REF_LIST_0].mv_union;
    context_ptr->mv_unit.mv[REF_LIST_1].mv_union =
        context_ptr->md_blk_arr_nsq[index_mds].prediction_unit_array[0].mv[REF_LIST_1].mv_union;
    uint8_t inter_pred_direction_index =
        (uint8_t)context_ptr->blk_ptr->prediction_unit_array->inter_pred_direction_index;
    uint8_t ref_frame_type = (uint8_t)context_ptr->blk_ptr->prediction_unit_array[0].ref_frame_type;
    int32_t is_inter       = (context_ptr->blk_ptr->prediction_mode_flag == INTER_MODE ||
                        context_ptr->blk_ptr->av1xd->use_intrabc)
        ? EB_TRUE
        : EB_FALSE;

    uint16_t tile_idx = context_ptr->tile_index;

    if (context_ptr->interpolation_search_level != IT_SEARCH_OFF)
        neighbor_array_unit_mode_write32(context_ptr->interpolation_type_neighbor_array,
                                         context_ptr->blk_ptr->interp_filters,
                                         origin_x,
                                         origin_y,
                                         bwdith,
                                         bheight,
                                         NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    {
        struct PartitionContext partition;
        partition.above = partition_context_lookup[context_ptr->blk_geom->bsize].above;
        partition.left  = partition_context_lookup[context_ptr->blk_geom->bsize].left;

        neighbor_array_unit_mode_write(context_ptr->leaf_partition_neighbor_array,
                                       (uint8_t *)(&partition), // NaderM
                                       origin_x,
                                       origin_y,
                                       bwdith,
                                       bheight,
                                       NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        // Mode Type Update
        neighbor_array_unit_mode_write(context_ptr->mode_type_neighbor_array,
                                       &mode_type,
                                       origin_x,
                                       origin_y,
                                       bwdith,
                                       bheight,
                                       NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        // Intra Luma Mode Update
        neighbor_array_unit_mode_write(context_ptr->intra_luma_mode_neighbor_array,
                                       &intra_luma_mode, //(uint8_t*)luma_mode,
                                       origin_x,
                                       origin_y,
                                       bwdith,
                                       bheight,
                                       NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        uint16_t txb_count = context_ptr->blk_geom->txb_count[context_ptr->blk_ptr->tx_depth];
        for (uint8_t txb_itr = 0; txb_itr < txb_count; txb_itr++) {
            uint8_t dc_sign_level_coeff = (int32_t)context_ptr->blk_ptr->quantized_dc[0][txb_itr];

            neighbor_array_unit_mode_write(
                context_ptr->luma_dc_sign_level_coeff_neighbor_array,
                (uint8_t *)&dc_sign_level_coeff,
                context_ptr->sb_origin_x +
                    context_ptr->blk_geom
                        ->tx_org_x[is_inter][context_ptr->blk_ptr->tx_depth][txb_itr],
                context_ptr->sb_origin_y +
                    context_ptr->blk_geom
                        ->tx_org_y[is_inter][context_ptr->blk_ptr->tx_depth][txb_itr],
                context_ptr->blk_geom->tx_width[context_ptr->blk_ptr->tx_depth][txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->blk_ptr->tx_depth][txb_itr],
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

            neighbor_array_unit_mode_write(
                pcs_ptr
                    ->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                           [tile_idx],
                (uint8_t *)&dc_sign_level_coeff,
                context_ptr->sb_origin_x +
                    context_ptr->blk_geom
                        ->tx_org_x[is_inter][context_ptr->blk_ptr->tx_depth][txb_itr],
                context_ptr->sb_origin_y +
                    context_ptr->blk_geom
                        ->tx_org_y[is_inter][context_ptr->blk_ptr->tx_depth][txb_itr],
                context_ptr->blk_geom->tx_width[context_ptr->blk_ptr->tx_depth][txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->blk_ptr->tx_depth][txb_itr],
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }
    }

    // Hsan: chroma mode rate estimation is kept even for chroma blind
    if (context_ptr->blk_geom->has_uv) {
        // Intra Chroma Mode Update
        neighbor_array_unit_mode_write(context_ptr->intra_chroma_mode_neighbor_array,
                                       &chroma_mode,
                                       blk_origin_x_uv,
                                       blk_origin_y_uv,
                                       bwdith_uv,
                                       bwheight_uv,
                                       NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    neighbor_array_unit_mode_write(context_ptr->skip_flag_neighbor_array,
                                   &skip_flag,
                                   origin_x,
                                   origin_y,
                                   bwdith,
                                   bheight,
                                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
        //  Update chroma CB cbf and Dc context
        {
            uint8_t dc_sign_level_coeff = (int32_t)context_ptr->blk_ptr->quantized_dc[1][0];
            neighbor_array_unit_mode_write(context_ptr->cb_dc_sign_level_coeff_neighbor_array,
                                           (uint8_t *)&dc_sign_level_coeff,
                                           blk_origin_x_uv,
                                           blk_origin_y_uv,
                                           bwdith_uv,
                                           bwheight_uv,
                                           NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }

        //  Update chroma CR cbf and Dc context
        {
            uint8_t dc_sign_level_coeff = (int32_t)context_ptr->blk_ptr->quantized_dc[2][0];
            neighbor_array_unit_mode_write(context_ptr->cr_dc_sign_level_coeff_neighbor_array,
                                           (uint8_t *)&dc_sign_level_coeff,
                                           blk_origin_x_uv,
                                           blk_origin_y_uv,
                                           bwdith_uv,
                                           bwheight_uv,
                                           NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }
    }
    uint8_t tx_size =
        tx_depth_to_tx_size[context_ptr->blk_ptr->tx_depth][context_ptr->blk_geom->bsize];
    uint8_t bw = tx_size_wide[tx_size];
    uint8_t bh = tx_size_high[tx_size];

    neighbor_array_unit_mode_write(context_ptr->txfm_context_array,
                                   &bw,
                                   origin_x,
                                   origin_y,
                                   bwdith,
                                   bheight,
                                   NEIGHBOR_ARRAY_UNIT_TOP_MASK);

    neighbor_array_unit_mode_write(context_ptr->txfm_context_array,
                                   &bh,
                                   origin_x,
                                   origin_y,
                                   bwdith,
                                   bheight,
                                   NEIGHBOR_ARRAY_UNIT_LEFT_MASK);

    // Update the Inter Pred Type Neighbor Array

    neighbor_array_unit_mode_write(context_ptr->inter_pred_dir_neighbor_array,
                                   &inter_pred_direction_index,
                                   origin_x,
                                   origin_y,
                                   bwdith,
                                   bheight,
                                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    // Update the refFrame Type Neighbor Array
    neighbor_array_unit_mode_write(context_ptr->ref_frame_type_neighbor_array,
                                   &ref_frame_type,
                                   origin_x,
                                   origin_y,
                                   bwdith,
                                   bheight,
                                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (!context_ptr->hbd_mode_decision) {
        update_recon_neighbor_array(
            context_ptr->luma_recon_neighbor_array,
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].neigh_top_recon[0],
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].neigh_left_recon[0],
            origin_x,
            origin_y,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight);
        if (context_ptr->md_tx_size_search_mode) {
            update_recon_neighbor_array(
                pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon[0],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[0],
                origin_x,
                origin_y,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight);
            update_recon_neighbor_array(
                pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon[0],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[0],
                origin_x,
                origin_y,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight);
        }

        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
            update_recon_neighbor_array(
                context_ptr->cb_recon_neighbor_array,
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon[1],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[1],
                blk_origin_x_uv,
                blk_origin_y_uv,
                bwdith_uv,
                bwheight_uv);
            update_recon_neighbor_array(
                context_ptr->cr_recon_neighbor_array,
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon[2],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[2],
                blk_origin_x_uv,
                blk_origin_y_uv,
                bwdith_uv,
                bwheight_uv);
        }
    } else {
        update_recon_neighbor_array16bit(
            context_ptr->luma_recon_neighbor_array16bit,
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                .neigh_top_recon_16bit[0],
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                .neigh_left_recon_16bit[0],
            origin_x,
            origin_y,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight);
        if (context_ptr->md_tx_size_search_mode) {
            update_recon_neighbor_array16bit(
                pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX]
                                                                     [tile_idx],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon_16bit[0],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon_16bit[0],
                origin_x,
                origin_y,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight);
            update_recon_neighbor_array16bit(
                pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX]
                                                                     [tile_idx],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon_16bit[0],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon_16bit[0],
                origin_x,
                origin_y,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight);
        }

        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
            update_recon_neighbor_array16bit(
                context_ptr->cb_recon_neighbor_array16bit,
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon_16bit[1],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon_16bit[1],
                blk_origin_x_uv,
                blk_origin_y_uv,
                bwdith_uv,
                bwheight_uv);
            update_recon_neighbor_array16bit(
                context_ptr->cr_recon_neighbor_array16bit,
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon_16bit[2],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon_16bit[2],
                blk_origin_x_uv,
                blk_origin_y_uv,
                bwdith_uv,
                bwheight_uv);
        }
    }

    return;
}

void copy_neighbour_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                           uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds, uint32_t sb_org_x,
                           uint32_t sb_org_y) {
    uint16_t tile_idx = context_ptr->tile_index;

    const BlockGeom *blk_geom = get_blk_geom_mds(blk_mds);

    uint32_t blk_org_x    = sb_org_x + blk_geom->origin_x;
    uint32_t blk_org_y    = sb_org_y + blk_geom->origin_y;
    uint32_t blk_org_x_uv = (blk_org_x >> 3 << 3) >> 1;
    uint32_t blk_org_y_uv = (blk_org_y >> 3 << 3) >> 1;
    uint32_t bwidth_uv    = blk_geom->bwidth_uv;
    uint32_t bheight_uv   = blk_geom->bheight_uv;

    copy_neigh_arr(pcs_ptr->md_intra_luma_mode_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_intra_luma_mode_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //neighbor_array_unit_reset(pcs_ptr->md_intra_chroma_mode_neighbor_array[depth]);
    copy_neigh_arr(pcs_ptr->md_intra_chroma_mode_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_intra_chroma_mode_neighbor_array[dst_idx][tile_idx],
                   blk_org_x_uv,
                   blk_org_y_uv,
                   bwidth_uv,
                   bheight_uv,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //neighbor_array_unit_reset(pcs_ptr->md_skip_flag_neighbor_array[depth]);
    copy_neigh_arr(pcs_ptr->md_skip_flag_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_skip_flag_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //neighbor_array_unit_reset(pcs_ptr->md_mode_type_neighbor_array[depth]);
    copy_neigh_arr(pcs_ptr->md_mode_type_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_mode_type_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    //neighbor_array_unit_reset(pcs_ptr->md_leaf_depth_neighbor_array[depth]);
    copy_neigh_arr(pcs_ptr->md_leaf_depth_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_leaf_depth_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    copy_neigh_arr(pcs_ptr->mdleaf_partition_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->mdleaf_partition_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (!context_ptr->hbd_mode_decision) {
        copy_neigh_arr(pcs_ptr->md_luma_recon_neighbor_array[src_idx][tile_idx],
                       pcs_ptr->md_luma_recon_neighbor_array[dst_idx][tile_idx],
                       blk_org_x,
                       blk_org_y,
                       blk_geom->bwidth,
                       blk_geom->bheight,
                       NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        if (context_ptr->md_tx_size_search_mode) {
            copy_neigh_arr(pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[src_idx][tile_idx],
                           pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[dst_idx][tile_idx],
                           blk_org_x,
                           blk_org_y,
                           blk_geom->bwidth,
                           blk_geom->bheight,
                           NEIGHBOR_ARRAY_UNIT_FULL_MASK);
            copy_neigh_arr(pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array[src_idx][tile_idx],
                           pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array[dst_idx][tile_idx],
                           blk_org_x,
                           blk_org_y,
                           blk_geom->bwidth,
                           blk_geom->bheight,
                           NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        }
        if (blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
            copy_neigh_arr(pcs_ptr->md_cb_recon_neighbor_array[src_idx][tile_idx],
                           pcs_ptr->md_cb_recon_neighbor_array[dst_idx][tile_idx],
                           blk_org_x_uv,
                           blk_org_y_uv,
                           bwidth_uv,
                           bheight_uv,
                           NEIGHBOR_ARRAY_UNIT_FULL_MASK);

            copy_neigh_arr(pcs_ptr->md_cr_recon_neighbor_array[src_idx][tile_idx],
                           pcs_ptr->md_cr_recon_neighbor_array[dst_idx][tile_idx],
                           blk_org_x_uv,
                           blk_org_y_uv,
                           bwidth_uv,
                           bheight_uv,
                           NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        }
    } else {
        copy_neigh_arr(pcs_ptr->md_luma_recon_neighbor_array16bit[src_idx][tile_idx],
                       pcs_ptr->md_luma_recon_neighbor_array16bit[dst_idx][tile_idx],
                       blk_org_x,
                       blk_org_y,
                       blk_geom->bwidth,
                       blk_geom->bheight,
                       NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        if (context_ptr->md_tx_size_search_mode) {
            copy_neigh_arr(pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[src_idx][tile_idx],
                           pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[dst_idx][tile_idx],
                           blk_org_x,
                           blk_org_y,
                           blk_geom->bwidth,
                           blk_geom->bheight,
                           NEIGHBOR_ARRAY_UNIT_FULL_MASK);
            copy_neigh_arr(pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array16bit[src_idx][tile_idx],
                           pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array16bit[dst_idx][tile_idx],
                           blk_org_x,
                           blk_org_y,
                           blk_geom->bwidth,
                           blk_geom->bheight,
                           NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        }

        if (blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
            copy_neigh_arr(pcs_ptr->md_cb_recon_neighbor_array16bit[src_idx][tile_idx],
                           pcs_ptr->md_cb_recon_neighbor_array16bit[dst_idx][tile_idx],
                           blk_org_x_uv,
                           blk_org_y_uv,
                           bwidth_uv,
                           bheight_uv,
                           NEIGHBOR_ARRAY_UNIT_FULL_MASK);

            copy_neigh_arr(pcs_ptr->md_cr_recon_neighbor_array16bit[src_idx][tile_idx],
                           pcs_ptr->md_cr_recon_neighbor_array16bit[dst_idx][tile_idx],
                           blk_org_x_uv,
                           blk_org_y_uv,
                           bwidth_uv,
                           bheight_uv,
                           NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        }
    }

    //neighbor_array_unit_reset(pcs_ptr->md_skip_coeff_neighbor_array[depth]);
    copy_neigh_arr(pcs_ptr->md_skip_coeff_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_skip_coeff_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    //neighbor_array_unit_reset(pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth]);
    copy_neigh_arr(pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    copy_neigh_arr(
        pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[src_idx][tile_idx],
        pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[dst_idx][tile_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
        copy_neigh_arr(pcs_ptr->md_cb_dc_sign_level_coeff_neighbor_array[src_idx][tile_idx],
                       pcs_ptr->md_cb_dc_sign_level_coeff_neighbor_array[dst_idx][tile_idx],
                       blk_org_x_uv,
                       blk_org_y_uv,
                       bwidth_uv,
                       bheight_uv,
                       NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        //neighbor_array_unit_reset(pcs_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth]);

        copy_neigh_arr(pcs_ptr->md_cr_dc_sign_level_coeff_neighbor_array[src_idx][tile_idx],
                       pcs_ptr->md_cr_dc_sign_level_coeff_neighbor_array[dst_idx][tile_idx],
                       blk_org_x_uv,
                       blk_org_y_uv,
                       bwidth_uv,
                       bheight_uv,
                       NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    //neighbor_array_unit_reset(pcs_ptr->md_txfm_context_array[depth]);
    copy_neigh_arr(pcs_ptr->md_txfm_context_array[src_idx][tile_idx],
                   pcs_ptr->md_txfm_context_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    //neighbor_array_unit_reset(pcs_ptr->md_inter_pred_dir_neighbor_array[depth]);
    copy_neigh_arr(pcs_ptr->md_inter_pred_dir_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_inter_pred_dir_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    //neighbor_array_unit_reset(pcs_ptr->md_ref_frame_type_neighbor_array[depth]);
    copy_neigh_arr(pcs_ptr->md_ref_frame_type_neighbor_array[src_idx][tile_idx],
                   pcs_ptr->md_ref_frame_type_neighbor_array[dst_idx][tile_idx],
                   blk_org_x,
                   blk_org_y,
                   blk_geom->bwidth,
                   blk_geom->bheight,
                   NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    copy_neigh_arr_32(pcs_ptr->md_interpolation_type_neighbor_array[src_idx][tile_idx],
                      pcs_ptr->md_interpolation_type_neighbor_array[dst_idx][tile_idx],
                      blk_org_x,
                      blk_org_y,
                      blk_geom->bwidth,
                      blk_geom->bheight,
                      NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
}

void md_update_all_neighbour_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                    uint32_t last_blk_index_mds, uint32_t sb_origin_x,
                                    uint32_t sb_origin_y) {
    context_ptr->blk_geom       = get_blk_geom_mds(last_blk_index_mds);
    context_ptr->blk_origin_x   = sb_origin_x + context_ptr->blk_geom->origin_x;
    context_ptr->blk_origin_y   = sb_origin_y + context_ptr->blk_geom->origin_y;
    context_ptr->round_origin_x = ((context_ptr->blk_origin_x >> 3) << 3);
    context_ptr->round_origin_y = ((context_ptr->blk_origin_y >> 3) << 3);

    context_ptr->blk_ptr   = &context_ptr->md_blk_arr_nsq[last_blk_index_mds];
    uint8_t avail_blk_flag = context_ptr->md_local_blk_unit[last_blk_index_mds].avail_blk_flag;

    if (avail_blk_flag) {
        mode_decision_update_neighbor_arrays(pcs_ptr, context_ptr, last_blk_index_mds, EB_FALSE);

        update_mi_map(context_ptr,
                      context_ptr->blk_ptr,
                      context_ptr->blk_origin_x,
                      context_ptr->blk_origin_y,
                      context_ptr->blk_geom,
                      avail_blk_flag,
                      pcs_ptr);
    }
}

void md_update_all_neighbour_arrays_multiple(PictureControlSet *  pcs_ptr,
                                             ModeDecisionContext *context_ptr, uint32_t blk_mds,
                                             uint32_t sb_origin_x, uint32_t sb_origin_y) {
    context_ptr->blk_geom = get_blk_geom_mds(blk_mds);

    uint32_t blk_it;
    for (blk_it = 0; blk_it < context_ptr->blk_geom->totns; blk_it++) {
        md_update_all_neighbour_arrays(
            pcs_ptr, context_ptr, blk_mds + blk_it, sb_origin_x, sb_origin_y);
    }
}

#define TOTAL_SQ_BLOCK_COUNT 341
int sq_block_index[TOTAL_SQ_BLOCK_COUNT] = {
    0,    25,   50,   75,   80,   81,   82,   83,   84,   89,   90,   91,   92,   93,   98,   99,
    100,  101,  102,  107,  108,  109,  110,  111,  136,  141,  142,  143,  144,  145,  150,  151,
    152,  153,  154,  159,  160,  161,  162,  163,  168,  169,  170,  171,  172,  197,  202,  203,
    204,  205,  206,  211,  212,  213,  214,  215,  220,  221,  222,  223,  224,  229,  230,  231,
    232,  233,  258,  263,  264,  265,  266,  267,  272,  273,  274,  275,  276,  281,  282,  283,
    284,  285,  290,  291,  292,  293,  294,  319,  344,  349,  350,  351,  352,  353,  358,  359,
    360,  361,  362,  367,  368,  369,  370,  371,  376,  377,  378,  379,  380,  405,  410,  411,
    412,  413,  414,  419,  420,  421,  422,  423,  428,  429,  430,  431,  432,  437,  438,  439,
    440,  441,  466,  471,  472,  473,  474,  475,  480,  481,  482,  483,  484,  489,  490,  491,
    492,  493,  498,  499,  500,  501,  502,  527,  532,  533,  534,  535,  536,  541,  542,  543,
    544,  545,  550,  551,  552,  553,  554,  559,  560,  561,  562,  563,  588,  613,  618,  619,
    620,  621,  622,  627,  628,  629,  630,  631,  636,  637,  638,  639,  640,  645,  646,  647,
    648,  649,  674,  679,  680,  681,  682,  683,  688,  689,  690,  691,  692,  697,  698,  699,
    700,  701,  706,  707,  708,  709,  710,  735,  740,  741,  742,  743,  744,  749,  750,  751,
    752,  753,  758,  759,  760,  761,  762,  767,  768,  769,  770,  771,  796,  801,  802,  803,
    804,  805,  810,  811,  812,  813,  814,  819,  820,  821,  822,  823,  828,  829,  830,  831,
    832,  857,  882,  887,  888,  889,  890,  891,  896,  897,  898,  899,  900,  905,  906,  907,
    908,  909,  914,  915,  916,  917,  918,  943,  948,  949,  950,  951,  952,  957,  958,  959,
    960,  961,  966,  967,  968,  969,  970,  975,  976,  977,  978,  979,  1004, 1009, 1010, 1011,
    1012, 1013, 1018, 1019, 1020, 1021, 1022, 1027, 1028, 1029, 1030, 1031, 1036, 1037, 1038, 1039,
    1040, 1065, 1070, 1071, 1072, 1073, 1074, 1079, 1080, 1081, 1082, 1083, 1088, 1089, 1090, 1091,
    1092, 1097, 1098, 1099, 1100};
void init_sq_nsq_block(SequenceControlSet *scs_ptr, ModeDecisionContext *context_ptr) {
    uint32_t blk_idx = 0;
    do {
        const BlockGeom *blk_geom                                       = get_blk_geom_mds(blk_idx);
        context_ptr->md_local_blk_unit[blk_idx].avail_blk_flag          = EB_FALSE;
        context_ptr->md_local_blk_unit[blk_idx].left_neighbor_partition = INVALID_NEIGHBOR_DATA;
        context_ptr->md_local_blk_unit[blk_idx].above_neighbor_partition = INVALID_NEIGHBOR_DATA;
        if (blk_geom->shape == PART_N) {
            context_ptr->md_blk_arr_nsq[blk_idx].split_flag         = EB_TRUE;
            context_ptr->md_blk_arr_nsq[blk_idx].part               = PARTITION_SPLIT;
            context_ptr->md_local_blk_unit[blk_idx].tested_blk_flag = EB_FALSE;
        }
        context_ptr->md_blk_arr_nsq[blk_idx].do_not_process_block = 0;
        ++blk_idx;
    } while (blk_idx < scs_ptr->max_block_cnt);
}
static INLINE TranHigh check_range(TranHigh input, int32_t bd) {
    // AV1 TX case
    // - 8 bit: signed 16 bit integer
    // - 10 bit: signed 18 bit integer
    // - 12 bit: signed 20 bit integer
    // - max quantization error = 1828 << (bd - 8)
    const int32_t int_max = (1 << (7 + bd)) - 1 + (914 << (bd - 7));
    const int32_t int_min = -int_max - 1;
#if CONFIG_COEFFICIENT_RANGE_CHECKING
    assert(int_min <= input);
    assert(input <= int_max);
#endif // CONFIG_COEFFICIENT_RANGE_CHECKING
    return (TranHigh)clamp64(input, int_min, int_max);
}

#define HIGHBD_WRAPLOW(x, bd) ((int32_t)check_range((x), bd))
static INLINE uint16_t highbd_clip_pixel_add(uint16_t dest, TranHigh trans, int32_t bd) {
    trans = HIGHBD_WRAPLOW(trans, bd);
    return clip_pixel_highbd(dest + (int32_t)trans, bd);
}

/*********************************
* Picture Single Channel Kernel
*********************************/
void picture_addition_kernel(uint8_t *pred_ptr, uint32_t pred_stride, int32_t *residual_ptr,
                             uint32_t residual_stride, uint8_t *recon_ptr, uint32_t recon_stride,
                             uint32_t width, uint32_t height, int32_t bd) {
    uint32_t column_index;
    uint32_t row_index = 0;
    //    const int32_t    maxValue = 0xFF;

    //SVT_LOG("\n");
    //SVT_LOG("Reconstruction---------------------------------------------------\n");

    while (row_index < height) {
        column_index = 0;
        while (column_index < width) {
            //recon_ptr[column_index] = (uint8_t)CLIP3(0, maxValue, ((int32_t)residual_ptr[column_index]) + ((int32_t)pred_ptr[column_index]));
            uint16_t rec            = (uint16_t)pred_ptr[column_index];
            recon_ptr[column_index] = (uint8_t)highbd_clip_pixel_add(
                rec, (TranLow)residual_ptr[column_index], bd);

            //SVT_LOG("%d\t", recon_ptr[column_index]);
            ++column_index;
        }

        //SVT_LOG("\n");
        residual_ptr += residual_stride;
        pred_ptr += pred_stride;
        recon_ptr += recon_stride;
        ++row_index;
    }
    //SVT_LOG("-----------------------------------------------------------------\n");
    //SVT_LOG("\n");
    //SVT_LOG("\n");
    return;
}

void picture_addition_kernel16_bit(uint16_t *pred_ptr, uint32_t pred_stride, int32_t *residual_ptr,
                                   uint32_t residual_stride, uint16_t *recon_ptr,
                                   uint32_t recon_stride, uint32_t width, uint32_t height,
                                   int32_t bd) {
    uint32_t column_index;
    uint32_t row_index = 0;
    //    const int32_t    maxValue = 0xFF;

    //SVT_LOG("\n");
    //SVT_LOG("Reconstruction---------------------------------------------------\n");

    while (row_index < height) {
        column_index = 0;
        while (column_index < width) {
            //recon_ptr[column_index] = (uint8_t)CLIP3(0, maxValue, ((int32_t)residual_ptr[column_index]) + ((int32_t)pred_ptr[column_index]));
            uint16_t rec            = (uint16_t)pred_ptr[column_index];
            recon_ptr[column_index] = highbd_clip_pixel_add(
                rec, (TranLow)residual_ptr[column_index], bd);

            //SVT_LOG("%d\t", recon_ptr[column_index]);
            ++column_index;
        }

        //SVT_LOG("\n");
        residual_ptr += residual_stride;
        pred_ptr += pred_stride;
        recon_ptr += recon_stride;
        ++row_index;
    }
    //    SVT_LOG("-----------------------------------------------------------------\n");
    //    SVT_LOG("\n");
    //    SVT_LOG("\n");
    return;
}

void av1_perform_inverse_transform_recon_luma(ModeDecisionContext *        context_ptr,
                                              ModeDecisionCandidateBuffer *candidate_buffer) {
    uint32_t txb_width;
    uint32_t txb_height;
    uint32_t txb_origin_x;
    uint32_t txb_origin_y;
    uint32_t txb_origin_index;
    uint32_t tu_total_count;
    uint32_t txb_itr;

    uint8_t tx_depth       = candidate_buffer->candidate_ptr->tx_depth;
    tu_total_count         = context_ptr->blk_geom->txb_count[tx_depth];
    txb_itr                = 0;
    uint32_t txb_1d_offset = 0;
    int32_t  is_inter      = (candidate_buffer->candidate_ptr->type == INTER_MODE ||
                        candidate_buffer->candidate_ptr->use_intrabc)
        ? EB_TRUE
        : EB_FALSE;
    do {
        txb_origin_x     = context_ptr->blk_geom->tx_org_x[is_inter][tx_depth][txb_itr];
        txb_origin_y     = context_ptr->blk_geom->tx_org_y[is_inter][tx_depth][txb_itr];
        txb_width        = context_ptr->blk_geom->tx_width[tx_depth][txb_itr];
        txb_height       = context_ptr->blk_geom->tx_height[tx_depth][txb_itr];
        txb_origin_index = txb_origin_x + txb_origin_y * candidate_buffer->prediction_ptr->stride_y;
        uint32_t rec_luma_offset = txb_origin_x +
            txb_origin_y * candidate_buffer->recon_ptr->stride_y;
        uint32_t y_has_coeff = (candidate_buffer->candidate_ptr->y_has_coeff & (1 << txb_itr)) > 0;

        if (y_has_coeff)
            inv_transform_recon_wrapper(candidate_buffer->prediction_ptr->buffer_y,
                                        txb_origin_index,
                                        candidate_buffer->prediction_ptr->stride_y,
                                        context_ptr->hbd_mode_decision
                                            ? (uint8_t *)context_ptr->cfl_temp_luma_recon16bit
                                            : context_ptr->cfl_temp_luma_recon,
                                        rec_luma_offset,
                                        candidate_buffer->recon_ptr->stride_y,
                                        (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_y,
                                        txb_1d_offset,
                                        context_ptr->hbd_mode_decision,
                                        context_ptr->blk_geom->txsize[tx_depth][txb_itr],
                                        candidate_buffer->candidate_ptr->transform_type[txb_itr],
                                        PLANE_TYPE_Y,
                                        (uint32_t)candidate_buffer->candidate_ptr->eob[0][txb_itr]);
        else {
            if (context_ptr->hbd_mode_decision) {
                pic_copy_kernel_16bit(
                    ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y) + txb_origin_index,
                    candidate_buffer->prediction_ptr->stride_y,
                    context_ptr->cfl_temp_luma_recon16bit + rec_luma_offset,
                    candidate_buffer->recon_ptr->stride_y,
                    txb_width,
                    txb_height);
            } else {
                pic_copy_kernel_8bit(
                    &(candidate_buffer->prediction_ptr->buffer_y[txb_origin_index]),
                    candidate_buffer->prediction_ptr->stride_y,
                    &(context_ptr->cfl_temp_luma_recon[rec_luma_offset]),
                    candidate_buffer->recon_ptr->stride_y,
                    txb_width,
                    txb_height);
            }
        }
        txb_1d_offset += context_ptr->blk_geom->tx_width[tx_depth][txb_itr] *
            context_ptr->blk_geom->tx_height[tx_depth][txb_itr];
        ++txb_itr;
    } while (txb_itr < tu_total_count);
}
void av1_perform_inverse_transform_recon(ModeDecisionContext *        context_ptr,
                                         ModeDecisionCandidateBuffer *candidate_buffer,
                                         const BlockGeom *            blk_geom) {
    uint32_t txb_width;
    uint32_t txb_height;
    uint32_t txb_origin_x;
    uint32_t txb_origin_y;
    uint32_t txb_origin_index;
    uint32_t tu_total_count;
    uint32_t txb_index;
    uint32_t txb_itr;

    UNUSED(blk_geom);

    uint8_t tx_depth       = candidate_buffer->candidate_ptr->tx_depth;
    tu_total_count         = context_ptr->blk_geom->txb_count[tx_depth];
    txb_index              = 0;
    txb_itr                = 0;
    uint32_t txb_1d_offset = 0, txb_1d_offset_uv = 0;
    uint32_t rec_luma_offset, rec_cb_offset, rec_cr_offset;
    int32_t  is_inter = (candidate_buffer->candidate_ptr->type == INTER_MODE ||
                        candidate_buffer->candidate_ptr->use_intrabc)
        ? EB_TRUE
        : EB_FALSE;

    do {
        txb_origin_x     = context_ptr->blk_geom->tx_org_x[is_inter][tx_depth][txb_itr];
        txb_origin_y     = context_ptr->blk_geom->tx_org_y[is_inter][tx_depth][txb_itr];
        txb_width        = context_ptr->blk_geom->tx_width[tx_depth][txb_itr];
        txb_height       = context_ptr->blk_geom->tx_height[tx_depth][txb_itr];
        rec_luma_offset  = txb_origin_x + txb_origin_y * candidate_buffer->recon_ptr->stride_y;
        rec_cb_offset    = ((((txb_origin_x >> 3) << 3) +
                          ((txb_origin_y >> 3) << 3) * candidate_buffer->recon_ptr->stride_cb) >>
                         1);
        rec_cr_offset    = ((((txb_origin_x >> 3) << 3) +
                          ((txb_origin_y >> 3) << 3) * candidate_buffer->recon_ptr->stride_cr) >>
                         1);
        txb_origin_index = txb_origin_x + txb_origin_y * candidate_buffer->prediction_ptr->stride_y;
        if (context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].y_has_coeff[txb_itr])
            inv_transform_recon_wrapper(candidate_buffer->prediction_ptr->buffer_y,
                                        txb_origin_index,
                                        candidate_buffer->prediction_ptr->stride_y,
                                        candidate_buffer->recon_ptr->buffer_y,
                                        rec_luma_offset,
                                        candidate_buffer->recon_ptr->stride_y,
                                        (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_y,
                                        txb_1d_offset,
                                        context_ptr->hbd_mode_decision,
                                        context_ptr->blk_geom->txsize[tx_depth][txb_itr],
                                        candidate_buffer->candidate_ptr->transform_type[txb_itr],
                                        PLANE_TYPE_Y,
                                        (uint32_t)candidate_buffer->candidate_ptr->eob[0][txb_itr]);
        else
            picture_copy(candidate_buffer->prediction_ptr,
                         txb_origin_index,
                         0, //txb_chroma_origin_index,
                         candidate_buffer->recon_ptr,
                         rec_luma_offset,
                         0, //txb_chroma_origin_index,
                         txb_width,
                         txb_height,
                         0, //chromaTuSize,
                         0, //chromaTuSize,
                         PICTURE_BUFFER_DESC_Y_FLAG,
                         context_ptr->hbd_mode_decision);

        //CHROMA
        uint8_t tx_depth = candidate_buffer->candidate_ptr->tx_depth;
        if (tx_depth == 0 || txb_itr == 0) {
            if (context_ptr->chroma_level <= CHROMA_MODE_1) {
                uint32_t chroma_txb_width =
                    tx_size_wide[context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr]];
                uint32_t chroma_txb_height =
                    tx_size_high[context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr]];
                uint32_t cb_tu_chroma_origin_index =
                    ((((txb_origin_x >> 3) << 3) +
                      ((txb_origin_y >> 3) << 3) * candidate_buffer->recon_coeff_ptr->stride_cb) >>
                     1);
                uint32_t cr_tu_chroma_origin_index =
                    ((((txb_origin_x >> 3) << 3) +
                      ((txb_origin_y >> 3) << 3) * candidate_buffer->recon_coeff_ptr->stride_cr) >>
                     1);

                if (context_ptr->blk_geom->has_uv &&
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .u_has_coeff[txb_index])
                    inv_transform_recon_wrapper(
                        candidate_buffer->prediction_ptr->buffer_cb,
                        cb_tu_chroma_origin_index,
                        candidate_buffer->prediction_ptr->stride_cb,
                        candidate_buffer->recon_ptr->buffer_cb,
                        rec_cb_offset,
                        candidate_buffer->recon_ptr->stride_cb,
                        (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cb,
                        txb_1d_offset_uv,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                        candidate_buffer->candidate_ptr->transform_type_uv,
                        PLANE_TYPE_UV,
                        (uint32_t)candidate_buffer->candidate_ptr->eob[1][txb_itr]);
                else
                    picture_copy(candidate_buffer->prediction_ptr,
                                 0,
                                 cb_tu_chroma_origin_index,
                                 candidate_buffer->recon_ptr,
                                 0,
                                 rec_cb_offset,
                                 0,
                                 0,
                                 chroma_txb_width,
                                 chroma_txb_height,
                                 PICTURE_BUFFER_DESC_Cb_FLAG,
                                 context_ptr->hbd_mode_decision);

                if (context_ptr->blk_geom->has_uv &&
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .v_has_coeff[txb_index])
                    inv_transform_recon_wrapper(
                        candidate_buffer->prediction_ptr->buffer_cr,
                        cr_tu_chroma_origin_index,
                        candidate_buffer->prediction_ptr->stride_cr,
                        candidate_buffer->recon_ptr->buffer_cr,
                        rec_cr_offset,
                        candidate_buffer->recon_ptr->stride_cr,
                        (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cr,
                        txb_1d_offset_uv,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                        candidate_buffer->candidate_ptr->transform_type_uv,
                        PLANE_TYPE_UV,
                        (uint32_t)candidate_buffer->candidate_ptr->eob[2][txb_itr]);
                else
                    picture_copy(candidate_buffer->prediction_ptr,
                                 0,
                                 cr_tu_chroma_origin_index,
                                 candidate_buffer->recon_ptr,
                                 0,
                                 rec_cr_offset,
                                 0,
                                 0,
                                 chroma_txb_width,
                                 chroma_txb_height,
                                 PICTURE_BUFFER_DESC_Cr_FLAG,
                                 context_ptr->hbd_mode_decision);

                if (context_ptr->blk_geom->has_uv)
                    txb_1d_offset_uv += context_ptr->blk_geom->tx_width_uv[tx_depth][txb_itr] *
                        context_ptr->blk_geom->tx_height_uv[tx_depth][txb_itr];
            }
        }
        txb_1d_offset += context_ptr->blk_geom->tx_width[tx_depth][txb_itr] *
            context_ptr->blk_geom->tx_height[tx_depth][txb_itr];
        ++txb_index;
        ++txb_itr;
    } while (txb_itr < tu_total_count);
}

/*******************************************
* Coding Loop - Fast Loop Initialization
*******************************************/
void product_coding_loop_init_fast_loop(ModeDecisionContext *context_ptr,
                                        NeighborArrayUnit *  skip_coeff_neighbor_array,
                                        NeighborArrayUnit *  inter_pred_dir_neighbor_array,
                                        NeighborArrayUnit *  ref_frame_type_neighbor_array,
                                        NeighborArrayUnit *  intra_luma_mode_neighbor_array,
                                        NeighborArrayUnit *  skip_flag_neighbor_array,
                                        NeighborArrayUnit *  mode_type_neighbor_array,
                                        NeighborArrayUnit *  leaf_depth_neighbor_array,
                                        NeighborArrayUnit *  leaf_partition_neighbor_array) {
    context_ptr->tx_depth = context_ptr->blk_ptr->tx_depth = 0;
    // Generate Split, Skip and intra mode contexts for the rate estimation
    coding_loop_context_generation(context_ptr,
                                   context_ptr->blk_ptr,
                                   context_ptr->blk_origin_x,
                                   context_ptr->blk_origin_y,
                                   BLOCK_SIZE_64,
                                   skip_coeff_neighbor_array,
                                   inter_pred_dir_neighbor_array,
                                   ref_frame_type_neighbor_array,
                                   intra_luma_mode_neighbor_array,
                                   skip_flag_neighbor_array,
                                   mode_type_neighbor_array,
                                   leaf_depth_neighbor_array,
                                   leaf_partition_neighbor_array);
    for (uint32_t index = 0; index < MAX_NFL_BUFF; ++index)
        context_ptr->fast_cost_array[index] = MAX_CU_COST;
    return;
}

void fast_loop_core(ModeDecisionCandidateBuffer *candidate_buffer, PictureControlSet *pcs_ptr,
                    ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                    uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                    uint32_t input_cr_origin_in_index, BlkStruct *blk_ptr, uint32_t cu_origin_index,
                    uint32_t cu_chroma_origin_index, EbBool use_ssd) {
    uint64_t luma_fast_distortion;
    uint64_t chroma_fast_distortion;
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];
    uint32_t fast_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->fast_lambda_md[EB_10_BIT_MD]
        : context_ptr->fast_lambda_md[EB_8_BIT_MD];

    ModeDecisionCandidate *candidate_ptr  = candidate_buffer->candidate_ptr;
    EbPictureBufferDesc *  prediction_ptr = candidate_buffer->prediction_ptr;
    context_ptr->pu_itr                   = 0;
    // Prediction
    // Set default interp_filters
    candidate_buffer->candidate_ptr->interp_filters = (context_ptr->md_staging_use_bilinear)
        ? av1_make_interp_filters(BILINEAR, BILINEAR)
        : 0;
    product_prediction_fun_table[candidate_buffer->candidate_ptr->use_intrabc
                                     ? INTER_MODE
                                     : candidate_ptr->type](
        context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);

    // Distortion
    // Y
    if (use_ssd) {
        EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
            ? full_distortion_kernel16_bits
            : spatial_full_distortion_kernel;

        candidate_buffer->candidate_ptr->luma_fast_distortion = (uint32_t)(
            luma_fast_distortion = spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                              input_origin_index,
                                                              input_picture_ptr->stride_y,
                                                              prediction_ptr->buffer_y,
                                                              (int32_t)cu_origin_index,
                                                              prediction_ptr->stride_y,
                                                              context_ptr->blk_geom->bwidth,
                                                              context_ptr->blk_geom->bheight));
    } else {
        assert((context_ptr->blk_geom->bwidth >> 3) < 17);
        if (!context_ptr->hbd_mode_decision) {
            candidate_buffer->candidate_ptr->luma_fast_distortion = (uint32_t)(
                luma_fast_distortion = nxm_sad_kernel_sub_sampled(
                    input_picture_ptr->buffer_y + input_origin_index,
                    input_picture_ptr->stride_y,
                    prediction_ptr->buffer_y + cu_origin_index,
                    prediction_ptr->stride_y,
                    context_ptr->blk_geom->bheight,
                    context_ptr->blk_geom->bwidth));
        } else {
            candidate_buffer->candidate_ptr->luma_fast_distortion = (uint32_t)(
                luma_fast_distortion = sad_16b_kernel(
                    ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                    input_picture_ptr->stride_y,
                    ((uint16_t *)prediction_ptr->buffer_y) + cu_origin_index,
                    prediction_ptr->stride_y,
                    context_ptr->blk_geom->bheight,
                    context_ptr->blk_geom->bwidth));
        }
    }

    if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1 &&
        context_ptr->md_staging_skip_inter_chroma_pred == EB_FALSE) {
        if (use_ssd) {
            EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
                ? full_distortion_kernel16_bits
                : spatial_full_distortion_kernel;

            chroma_fast_distortion = spatial_full_dist_type_fun(
                input_picture_ptr->buffer_cb,
                input_cb_origin_in_index,
                input_picture_ptr->stride_cb,
                candidate_buffer->prediction_ptr->buffer_cb,
                (int32_t)cu_chroma_origin_index,
                prediction_ptr->stride_cb,
                context_ptr->blk_geom->bwidth_uv,
                context_ptr->blk_geom->bheight_uv);

            chroma_fast_distortion += spatial_full_dist_type_fun(
                input_picture_ptr->buffer_cr,
                input_cr_origin_in_index,
                input_picture_ptr->stride_cb,
                candidate_buffer->prediction_ptr->buffer_cr,
                (int32_t)cu_chroma_origin_index,
                prediction_ptr->stride_cr,
                context_ptr->blk_geom->bwidth_uv,
                context_ptr->blk_geom->bheight_uv);
        } else {
            assert((context_ptr->blk_geom->bwidth_uv >> 3) < 17);

            if (!context_ptr->hbd_mode_decision) {
                chroma_fast_distortion = nxm_sad_kernel_sub_sampled(
                    input_picture_ptr->buffer_cb + input_cb_origin_in_index,
                    input_picture_ptr->stride_cb,
                    candidate_buffer->prediction_ptr->buffer_cb + cu_chroma_origin_index,
                    prediction_ptr->stride_cb,
                    context_ptr->blk_geom->bheight_uv,
                    context_ptr->blk_geom->bwidth_uv);

                chroma_fast_distortion += nxm_sad_kernel_sub_sampled(
                    input_picture_ptr->buffer_cr + input_cr_origin_in_index,
                    input_picture_ptr->stride_cr,
                    candidate_buffer->prediction_ptr->buffer_cr + cu_chroma_origin_index,
                    prediction_ptr->stride_cr,
                    context_ptr->blk_geom->bheight_uv,
                    context_ptr->blk_geom->bwidth_uv);
            } else {
                chroma_fast_distortion = sad_16b_kernel(
                    ((uint16_t *)input_picture_ptr->buffer_cb) + input_cb_origin_in_index,
                    input_picture_ptr->stride_cb,
                    ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cb) +
                        cu_chroma_origin_index,
                    prediction_ptr->stride_cb,
                    context_ptr->blk_geom->bheight_uv,
                    context_ptr->blk_geom->bwidth_uv);

                chroma_fast_distortion += sad_16b_kernel(
                    ((uint16_t *)input_picture_ptr->buffer_cr) + input_cr_origin_in_index,
                    input_picture_ptr->stride_cr,
                    ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cr) +
                        cu_chroma_origin_index,
                    prediction_ptr->stride_cr,
                    context_ptr->blk_geom->bheight_uv,
                    context_ptr->blk_geom->bwidth_uv);
            }
        }
    } else
        chroma_fast_distortion = 0;
    // Fast Cost
    *(candidate_buffer->fast_cost_ptr) = av1_product_fast_cost_func_table[candidate_ptr->type](
        blk_ptr,
        candidate_buffer->candidate_ptr,
        blk_ptr->qp,
        luma_fast_distortion,
        chroma_fast_distortion,
        use_ssd ? full_lambda : fast_lambda,
        use_ssd,
        pcs_ptr,
        &(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
              .ed_ref_mv_stack[candidate_ptr->ref_frame_type][0]),
        context_ptr->blk_geom,
        context_ptr->blk_origin_y >> MI_SIZE_LOG2,
        context_ptr->blk_origin_x >> MI_SIZE_LOG2,
        context_ptr->md_enable_inter_intra,
        context_ptr->full_cost_shut_fast_rate_flag,
        1,
        context_ptr->intra_luma_left_mode,
        context_ptr->intra_luma_top_mode);

    // Init full cost in case we by pass stage1/stage2
    if (context_ptr->md_staging_mode == MD_STAGING_MODE_0)
        *(candidate_buffer->full_cost_ptr) = *(candidate_buffer->fast_cost_ptr);
}
static const int32_t pd0_nic[MD_STAGE_TOTAL - 1][MAX_FRAME_TYPE][CAND_CLASS_TOTAL] = {
    {
        //MD_STAGE_1
        // C0  C1  C2  C3  C4  C5  C6  C7  C8
        {1, 1, 1, 1, 1, 1, 1, 1, 1}, // I_SLICE
        {1, 1, 1, 1, 1, 1, 1, 1, 1}, // REFERENCE_FRAME
        {1, 1, 1, 1, 1, 1, 1, 1, 1} // non-REFERENCE_FRAME
    },
    {//MD_STAGE_2
     {1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1}},
    {//MD_STAGE_3
     {1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1}}};

static const int32_t pd1_nic[MD_STAGE_TOTAL - 1][MAX_FRAME_TYPE][CAND_CLASS_TOTAL] = {
    {
        //MD_STAGE_1
        // C0       C1    C2    C3    C4    C5    C6    C7    C8
        {ALL_S0, 0, 0, 0, 0, 16, 10, 14, 4}, // I_SLICE
        {16, 16, 16, 16, 14, 16, 10, 14, 4}, // REFERENCE_FRAME
        {8, 8, 8, 8, 6, 16, 5, 14, 2} // non-REFERENCE_FRAME
    },
    {//MD_STAGE_2
     {10, 0, 0, 0, 0, 12, 5, 7, 1},
     {10, 6, 6, 6, 12, 12, 5, 7, 1},
     {4, 3, 3, 3, 4, 4, 2, 7, 1}},
    {//MD_STAGE_3
     {10, 0, 0, 0, 0, 12, 5, 7, 1},
     {10, 6, 6, 6, 12, 12, 5, 7, 1},
     {4, 3, 3, 3, 4, 4, 2, 7, 1}}};
static const int32_t pd2_nic[MD_STAGE_TOTAL - 1][MAX_FRAME_TYPE][CAND_CLASS_TOTAL] = {
    {
        //MD_STAGE_1
        // C0       C1    C2    C3    C4    C5     C6    C7   C8
        {ALL_S0, 0, 0, 0, 0, 16, 5, 16, 8}, // I_SLICE
        {48, 12, 12, 8, 8, 8, 5, 8, 6}, // REFERENCE_FRAME
        {8, 6, 6, 4, 4, 4, 3, 4, 4} // non-REFERENCE_FRAME
    },
    {//MD_STAGE_2
     {32, 0, 0, 0, 0, 8, 5, 8, 6},
     {16, 6, 6, 4, 4, 4, 3, 4, 4},
     {4, 3, 3, 2, 2, 2, 2, 2, 2}},
    {//MD_STAGE_3
     {32, 0, 0, 0, 0, 8, 5, 8, 6},
     {16, 6, 6, 4, 4, 4, 3, 4, 4},
     {4, 3, 3, 2, 2, 2, 2, 2, 2}}};
void set_md_stage_counts(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                         uint32_t fastCandidateTotalCount) {
    // Step 1: derive bypass_stage1 flags
    if (context_ptr->md_staging_mode == MD_STAGING_MODE_1 ||
        context_ptr->md_staging_mode == MD_STAGING_MODE_2)
        memset(context_ptr->bypass_md_stage_1, EB_FALSE, CAND_CLASS_TOTAL);
    else
        memset(context_ptr->bypass_md_stage_1, EB_TRUE, CAND_CLASS_TOTAL);

    // Step 1: derive bypass_stage_2 flags
    if (context_ptr->md_staging_mode == MD_STAGING_MODE_2)
        memset(context_ptr->bypass_md_stage_2, EB_FALSE, CAND_CLASS_TOTAL);
    else
        memset(context_ptr->bypass_md_stage_2, EB_TRUE, CAND_CLASS_TOTAL);
    uint8_t cand_it  = 0;
    uint8_t frm_type = pcs_ptr->slice_type == I_SLICE
        ? 0
        : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 2;
    // nic scale factor
    uint8_t scale_num;
    uint8_t scale_denum;
    if (context_ptr->nic_level == 0) {
        scale_num   = 1;
        scale_denum = 1;
    } else if (context_ptr->nic_level == 1) {
        scale_num   = 3;
        scale_denum = 4;
    } else if (context_ptr->nic_level == 2) {
        scale_num   = 2;
        scale_denum = 3;
    } else {
        scale_num   = 1;
        scale_denum = 2;
    }
    if (context_ptr->md_staging_count_level == 0) {
        // Set NICS for PD0
        for (cand_it = CAND_CLASS_0; cand_it < CAND_CLASS_TOTAL; ++cand_it) {
            context_ptr->md_stage_1_count[cand_it] = (pd0_nic[MD_STAGE_1 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd0_nic[MD_STAGE_1 - 1][frm_type][cand_it];
            context_ptr->md_stage_2_count[cand_it] = (pd0_nic[MD_STAGE_2 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd0_nic[MD_STAGE_2 - 1][frm_type][cand_it];
            context_ptr->md_stage_3_count[cand_it] = (pd0_nic[MD_STAGE_3 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd0_nic[MD_STAGE_3 - 1][frm_type][cand_it];
        }
    } else if (context_ptr->md_staging_count_level == 1) {
        // Set NICS for PD1
        for (cand_it = CAND_CLASS_0; cand_it < CAND_CLASS_TOTAL; ++cand_it) {
            context_ptr->md_stage_1_count[cand_it] = (pd1_nic[MD_STAGE_1 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd1_nic[MD_STAGE_1 - 1][frm_type][cand_it];
            context_ptr->md_stage_2_count[cand_it] = (pd1_nic[MD_STAGE_2 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd1_nic[MD_STAGE_2 - 1][frm_type][cand_it];
            context_ptr->md_stage_3_count[cand_it] = (pd1_nic[MD_STAGE_3 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd1_nic[MD_STAGE_3 - 1][frm_type][cand_it];
        }
    } else {
        // Set NICS for PD2
        for (cand_it = CAND_CLASS_0; cand_it < CAND_CLASS_TOTAL; ++cand_it) {
            context_ptr->md_stage_1_count[cand_it] = (pd2_nic[MD_STAGE_1 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd2_nic[MD_STAGE_1 - 1][frm_type][cand_it];
            // apply scale factor and set MIN to 1 candidate
            context_ptr->md_stage_1_count[cand_it] = MAX(
                DIVIDE_AND_ROUND(scale_num * context_ptr->md_stage_1_count[cand_it], scale_denum),
                1);

            context_ptr->md_stage_2_count[cand_it] = (pd2_nic[MD_STAGE_2 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd2_nic[MD_STAGE_2 - 1][frm_type][cand_it];
            // apply scale factor and set MIN to 1 candidate
            context_ptr->md_stage_2_count[cand_it] = MAX(
                DIVIDE_AND_ROUND(scale_num * context_ptr->md_stage_2_count[cand_it],
                                 scale_denum),
                1);

            context_ptr->md_stage_3_count[cand_it] = (pd2_nic[MD_STAGE_3 - 1][frm_type][cand_it] ==
                                                      ALL_S0)
                ? fastCandidateTotalCount
                : (uint32_t)pd2_nic[MD_STAGE_3 - 1][frm_type][cand_it];
            // apply scale factor and set MIN to 1 candidate
            context_ptr->md_stage_3_count[cand_it] = MAX(
                DIVIDE_AND_ROUND(scale_num * context_ptr->md_stage_3_count[cand_it],
                                 scale_denum),
                1);
        }

        if ((pcs_ptr->enc_mode <= ENC_M0 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected)) ||
            ((pcs_ptr->enc_mode <= ENC_M0 ||
              (pcs_ptr->enc_mode <= ENC_M1 && pcs_ptr->parent_pcs_ptr->sc_content_detected)) &&
             context_ptr->blk_geom->shape == PART_N)) {
            uint8_t mult_factor_num   = 5;
            uint8_t mult_factor_denum = 4;
            for (uint8_t i = 0; i < CAND_CLASS_TOTAL; ++i) {
                if (i == CAND_CLASS_0 || i == CAND_CLASS_6 || i == CAND_CLASS_7) {
                    // INTRA scaling
                    if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
                        mult_factor_num   = 5;
                        mult_factor_denum = 4;
                    } else {
                        mult_factor_num   = 1;
                        mult_factor_denum = 1;
                    }
                } else {
                    // INTER scaling
                    if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
                        mult_factor_num   = 1;
                        mult_factor_denum = 1;
                    } else {
                        mult_factor_num   = 5;
                        mult_factor_denum = 4;
                    }
                }
                context_ptr->md_stage_1_count[i] = DIVIDE_AND_ROUND(
                    mult_factor_num * context_ptr->md_stage_1_count[i],
                    mult_factor_denum);
                context_ptr->md_stage_3_count[i] = DIVIDE_AND_ROUND(
                    mult_factor_num * context_ptr->md_stage_3_count[i],
                    mult_factor_denum);
            }
        }

        if (!(pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)) {
            uint8_t mult_factor_num   = 4;
            uint8_t mult_factor_denum = 3;
            for (uint8_t i = 0; i < CAND_CLASS_TOTAL; ++i) {
                context_ptr->md_stage_1_count[i] = DIVIDE_AND_ROUND(
                    mult_factor_num * context_ptr->md_stage_1_count[i],
                    mult_factor_denum);
                context_ptr->md_stage_3_count[i] = DIVIDE_AND_ROUND(
                    mult_factor_num * context_ptr->md_stage_3_count[i],
                    mult_factor_denum);
            }
        }

        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            ////DIVIDE
            uint8_t division_factor_num   = 1;
            uint8_t division_factor_denum = 1;
            if (pcs_ptr->enc_mode <= ENC_M0) {
                division_factor_num   = 1;
                division_factor_denum = 1;
            } else if (pcs_ptr->enc_mode <= ENC_M1) {
                division_factor_num   = 1;
                division_factor_denum = 1;
            } else {
                division_factor_num   = 7;
                division_factor_denum = 8;
            }

            for (uint8_t i = 0; i < CAND_CLASS_TOTAL; ++i) {
                if (i != CAND_CLASS_0 && i != CAND_CLASS_6 && i != CAND_CLASS_7) {
                    context_ptr->md_stage_1_count[i] = DIVIDE_AND_ROUND(
                        division_factor_num * context_ptr->md_stage_1_count[i],
                        division_factor_denum);
                    context_ptr->md_stage_1_count[i] = MAX(context_ptr->md_stage_1_count[i], 1);
                    context_ptr->md_stage_3_count[i] = DIVIDE_AND_ROUND(
                        division_factor_num * context_ptr->md_stage_3_count[i],
                        division_factor_denum);
                    context_ptr->md_stage_3_count[i] = MAX(context_ptr->md_stage_3_count[i], 1);
                }
            }
        } else {
            ////DIVIDE
            uint8_t division_factor_num   = 1;
            uint8_t division_factor_denum = 1;
            if (pcs_ptr->enc_mode <= ENC_M0) {
                division_factor_num   = 1;
                division_factor_denum = 1;
            } else if (pcs_ptr->enc_mode <= ENC_M1) {
                division_factor_num   = 1;
                division_factor_denum = 1;
            } else {
                division_factor_num   = 3;
                division_factor_denum = 4;
            }

            for (uint8_t i = 0; i < CAND_CLASS_TOTAL; ++i) {
                if (i != CAND_CLASS_0 && i != CAND_CLASS_6 && i != CAND_CLASS_7) {
                    context_ptr->md_stage_1_count[i] = DIVIDE_AND_ROUND(
                        division_factor_num * context_ptr->md_stage_1_count[i],
                        division_factor_denum);
                    context_ptr->md_stage_1_count[i] = MAX(context_ptr->md_stage_1_count[i], 1);
                    context_ptr->md_stage_3_count[i] = DIVIDE_AND_ROUND(
                        division_factor_num * context_ptr->md_stage_3_count[i],
                        division_factor_denum);
                    context_ptr->md_stage_3_count[i] = MAX(context_ptr->md_stage_3_count[i], 1);
                }
            }
        }
        if (pcs_ptr->enc_mode > ENC_M0 || pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            uint8_t division_factor_num   = 1;
            uint8_t division_factor_denum = 1;
            if (context_ptr->blk_geom->bheight <= 8 && context_ptr->blk_geom->bwidth <= 8) {
                division_factor_num   = 2;
                division_factor_denum = 3;
            } else if (context_ptr->blk_geom->bheight <= 16 &&
                       context_ptr->blk_geom->bwidth <= 16) {
                division_factor_num   = 3;
                division_factor_denum = 4;
            } else if (context_ptr->blk_geom->bheight <= 32 &&
                       context_ptr->blk_geom->bwidth <= 32) {
                division_factor_num   = 7;
                division_factor_denum = 8;
            } else {
                division_factor_num   = 1;
                division_factor_denum = 1;
            }

            for (uint8_t i = 0; i < CAND_CLASS_TOTAL; ++i) {
                context_ptr->md_stage_1_count[i] = DIVIDE_AND_ROUND(
                    division_factor_num * context_ptr->md_stage_1_count[i],
                    division_factor_denum);
                context_ptr->md_stage_1_count[i] = MAX(context_ptr->md_stage_1_count[i], 1);
                context_ptr->md_stage_3_count[i] = DIVIDE_AND_ROUND(
                    division_factor_num * context_ptr->md_stage_3_count[i],
                    division_factor_denum);
                context_ptr->md_stage_3_count[i] = MAX(context_ptr->md_stage_3_count[i], 1);
            }
        }
    }

    //update NICs of bypassed stage n using NICs from stage n+1
    for (cand_it = CAND_CLASS_0; cand_it < CAND_CLASS_TOTAL; ++cand_it) {
        context_ptr->md_stage_2_count[cand_it] = context_ptr->bypass_md_stage_2[cand_it]
            ? context_ptr->md_stage_3_count[cand_it]
            : context_ptr->md_stage_2_count[cand_it];
        context_ptr->md_stage_1_count[cand_it] = context_ptr->bypass_md_stage_1[cand_it]
            ? context_ptr->md_stage_2_count[cand_it]
            : context_ptr->md_stage_1_count[cand_it];
    }
}
void sort_fast_cost_based_candidates(
    struct ModeDecisionContext *context_ptr, uint32_t input_buffer_start_idx,
    uint32_t
              input_buffer_count, //how many cand buffers to sort. one of the buffer can have max cost.
    uint32_t *cand_buff_indices) {
    ModeDecisionCandidateBuffer **buffer_ptr_array = context_ptr->candidate_buffer_ptr_array;
    uint32_t input_buffer_end_idx = input_buffer_start_idx + input_buffer_count - 1;
    uint32_t buffer_index, i, j;
    uint32_t k = 0;
    for (buffer_index = input_buffer_start_idx; buffer_index <= input_buffer_end_idx;
         buffer_index++, k++) {
        cand_buff_indices[k] = buffer_index;
    }
    for (i = 0; i < input_buffer_count - 1; ++i) {
        for (j = i + 1; j < input_buffer_count; ++j) {
            if (*(buffer_ptr_array[cand_buff_indices[j]]->fast_cost_ptr) <
                *(buffer_ptr_array[cand_buff_indices[i]]->fast_cost_ptr)) {
                buffer_index         = cand_buff_indices[i];
                cand_buff_indices[i] = (uint32_t)cand_buff_indices[j];
                cand_buff_indices[j] = (uint32_t)buffer_index;
            }
        }
    }
}

static INLINE void heap_sort_stage_max_node_fast_cost_ptr(ModeDecisionCandidateBuffer **buffer_ptr,
                                                          uint32_t *sort_index, uint32_t i,
                                                          uint32_t num) {
    uint32_t left, right, max;

    /* Loop for removing recursion. */
    while (1) {
        left  = 2 * i;
        right = 2 * i + 1;
        max   = i;

        if (left <= num &&
            *(buffer_ptr[sort_index[left]]->fast_cost_ptr) >
                *(buffer_ptr[sort_index[i]]->fast_cost_ptr)) {
            max = left;
        }

        if (right <= num &&
            *(buffer_ptr[sort_index[right]]->fast_cost_ptr) >
                *(buffer_ptr[sort_index[max]]->fast_cost_ptr)) {
            max = right;
        }

        if (max == i) {
            break;
        }

        uint32_t swap   = sort_index[i];
        sort_index[i]   = sort_index[max];
        sort_index[max] = swap;
        i               = max;
    }
}

static void qsort_stage_max_node_fast_cost_ptr(ModeDecisionCandidateBuffer **buffer_ptr_array,
                                               uint32_t *dst, uint32_t *a, uint32_t *b, int num) {
    if (num < 4) {
        if (num < 2) {
            if (num) {
                //num = 1
                dst[0] = a[0];
            }
            return;
        }
        if (num > 2) {
            //num = 3
            uint32_t tmp_a = a[0];
            uint32_t tmp_b = a[1];
            uint32_t tmp_c = a[2];
            uint64_t val_a = *(buffer_ptr_array[tmp_a]->fast_cost_ptr);
            uint64_t val_b = *(buffer_ptr_array[tmp_b]->fast_cost_ptr);
            uint64_t val_c = *(buffer_ptr_array[tmp_c]->fast_cost_ptr);

            if (val_a < val_b) {
                if (val_b < val_c) {
                    //Sorted abc
                    dst[0] = tmp_a;
                    dst[1] = tmp_b;
                    dst[2] = tmp_c;
                } else {
                    //xcx
                    if (val_a < val_c) {
                        //Sorted 132
                        dst[0] = tmp_a;
                        dst[1] = tmp_c;
                        dst[2] = tmp_b;
                    } else {
                        //Sorted 231
                        dst[0] = tmp_c;
                        dst[1] = tmp_a;
                        dst[2] = tmp_b;
                    }
                }
            } else {
                //a>b
                if (val_b > val_c) {
                    //Sorted cba
                    dst[0] = tmp_c;
                    dst[1] = tmp_b;
                    dst[2] = tmp_a;
                } else {
                    //bxx
                    if (val_a < val_c) {
                        //Sorted bac
                        dst[0] = tmp_b;
                        dst[1] = tmp_a;
                        dst[2] = tmp_c;
                    } else {
                        //Sorted bca
                        dst[0] = tmp_b;
                        dst[1] = tmp_c;
                        dst[2] = tmp_a;
                    }
                }
            }
            return;
        }

        /* bacuse a and dst can point on this same array, copy temporary values*/
        uint32_t tmp_a = a[0];
        uint32_t tmp_b = a[1];
        if (*(buffer_ptr_array[tmp_a]->fast_cost_ptr) < *(buffer_ptr_array[tmp_b]->fast_cost_ptr)) {
            dst[0] = tmp_a;
            dst[1] = tmp_b;
        } else {
            dst[0] = tmp_b;
            dst[1] = tmp_a;
        }
        return;
    }

    int sorted_down = 0;
    int sorted_up   = num - 1;

    uint64_t pivot_val = *(buffer_ptr_array[a[0]]->fast_cost_ptr);
    for (int i = 1; i < num; ++i) {
        if (pivot_val < *(buffer_ptr_array[a[i]]->fast_cost_ptr)) {
            b[sorted_up] = a[i];
            sorted_up--;
        } else {
            b[sorted_down] = a[i];
            sorted_down++;
        }
    }

    dst[sorted_down] = a[0];

    qsort_stage_max_node_fast_cost_ptr(buffer_ptr_array, dst, b, a, sorted_down);

    qsort_stage_max_node_fast_cost_ptr(buffer_ptr_array,
                                       dst + (sorted_down + 1),
                                       b + (sorted_down + 1),
                                       a + (sorted_down + 1),
                                       num - (sorted_down)-1);
}

static INLINE void sort_array_index_fast_cost_ptr(ModeDecisionCandidateBuffer **buffer_ptr,
                                                  uint32_t *sort_index, uint32_t num) {
    if (num <= 60) {
        //For small array uses 'quick sort', work much faster for small array,
        //but required alloc temporary memory.
        uint32_t sorted_tmp[60];
        qsort_stage_max_node_fast_cost_ptr(buffer_ptr, sort_index, sort_index, sorted_tmp, num);
        return;
    }

    //For big arrays uses 'heap sort', not need allocate memory
    //For small array less that 40 elements heap sort work slower than 'insertion sort'
    uint32_t i;
    for (i = (num - 1) / 2; i > 0; i--) {
        heap_sort_stage_max_node_fast_cost_ptr(buffer_ptr, sort_index, i, num - 1);
    }

    heap_sort_stage_max_node_fast_cost_ptr(buffer_ptr, sort_index, 0, num - 1);

    for (i = num - 1; i > 0; i--) {
        uint32_t swap = sort_index[i];
        sort_index[i] = sort_index[0];
        sort_index[0] = swap;
        heap_sort_stage_max_node_fast_cost_ptr(buffer_ptr, sort_index, 0, i - 1);
    }
}
void sort_full_cost_based_candidates(struct ModeDecisionContext *context_ptr,
                                     uint32_t num_of_cand_to_sort, uint32_t *cand_buff_indices) {
    uint32_t                      i, j, index;
    ModeDecisionCandidateBuffer **buffer_ptr_array = context_ptr->candidate_buffer_ptr_array;
    for (i = 0; i < num_of_cand_to_sort - 1; ++i) {
        for (j = i + 1; j < num_of_cand_to_sort; ++j) {
            if (*(buffer_ptr_array[cand_buff_indices[j]]->full_cost_ptr) <
                *(buffer_ptr_array[cand_buff_indices[i]]->full_cost_ptr)) {
                index                = cand_buff_indices[i];
                cand_buff_indices[i] = (uint32_t)cand_buff_indices[j];
                cand_buff_indices[j] = (uint32_t)index;
            }
        }
    }
}
void construct_best_sorted_arrays_md_stage_1(struct ModeDecisionContext *  context_ptr,
                                             ModeDecisionCandidateBuffer **buffer_ptr_array,
                                             uint32_t *best_candidate_index_array,
                                             uint32_t *sorted_candidate_index_array,
                                             uint64_t *ref_fast_cost) {
    //best = union from all classes
    uint32_t best_candi = 0;
    for (CandClass class_i = CAND_CLASS_0; class_i < CAND_CLASS_TOTAL; class_i++)
        for (uint32_t candi = 0; candi < context_ptr->md_stage_1_count[class_i]; candi++)
            sorted_candidate_index_array[best_candi++] =
                context_ptr->cand_buff_indices[class_i][candi];

    assert(best_candi == context_ptr->md_stage_1_total_count);
    uint32_t full_recon_candidate_count = context_ptr->md_stage_1_total_count;

    //sort best: inter, then intra
    uint32_t i, id;
    uint32_t id_inter = 0;
    uint32_t id_intra = full_recon_candidate_count - 1;
    for (i = 0; i < full_recon_candidate_count; ++i) {
        id = sorted_candidate_index_array[i];
        if (buffer_ptr_array[id]->candidate_ptr->type == INTER_MODE) {
            best_candidate_index_array[id_inter++] = id;
        } else {
            assert(buffer_ptr_array[id]->candidate_ptr->type == INTRA_MODE);
            best_candidate_index_array[id_intra--] = id;
        }
    }

    //sorted best: *(buffer_ptr_array[sorted_candidate_index_array[?]]->fast_cost_ptr)
    sort_array_index_fast_cost_ptr(
        buffer_ptr_array, sorted_candidate_index_array, full_recon_candidate_count);

    // tx search
    *ref_fast_cost = *(buffer_ptr_array[sorted_candidate_index_array[0]]->fast_cost_ptr);
}

void construct_best_sorted_arrays_md_stage_3(struct ModeDecisionContext *  context_ptr,
                                             ModeDecisionCandidateBuffer **buffer_ptr_array,
                                             uint32_t *best_candidate_index_array,
                                             uint32_t *sorted_candidate_index_array) {
    //best = union from all classes
    uint32_t best_candi = 0;
    for (CandClass class_i = CAND_CLASS_0; class_i < CAND_CLASS_TOTAL; class_i++)
        for (uint32_t candi = 0; candi < context_ptr->md_stage_3_count[class_i]; candi++)
            sorted_candidate_index_array[best_candi++] =
                context_ptr->cand_buff_indices[class_i][candi];

    assert(best_candi == context_ptr->md_stage_3_total_count);
    uint32_t fullReconCandidateCount = context_ptr->md_stage_3_total_count;

    //sort best: inter, then intra
    uint32_t i, id;
    uint32_t id_inter = 0;
    uint32_t id_intra = fullReconCandidateCount - 1;
    if (context_ptr->chroma_at_last_md_stage)
        context_ptr->md_stage_3_total_intra_count = 0;
    for (i = 0; i < fullReconCandidateCount; ++i) {
        id = sorted_candidate_index_array[i];
        if (buffer_ptr_array[id]->candidate_ptr->type == INTER_MODE) {
            best_candidate_index_array[id_inter++] = id;
        } else {
            assert(buffer_ptr_array[id]->candidate_ptr->type == INTRA_MODE);
            best_candidate_index_array[id_intra--] = id;
        }
        if (context_ptr->chroma_at_last_md_stage) {
            uint8_t is_inter = (buffer_ptr_array[id]->candidate_ptr->type == INTER_MODE ||
                                buffer_ptr_array[id]->candidate_ptr->use_intrabc)
                ? EB_TRUE
                : EB_FALSE;
            context_ptr->md_stage_3_total_intra_count += !is_inter ? 1 : 0;
        }
    }
    context_ptr->best_intra_cost = MAX_MODE_COST;
    context_ptr->best_inter_cost = MAX_MODE_COST;
    for (i = 0; i < fullReconCandidateCount; ++i) {
        id               = sorted_candidate_index_array[i];
        int32_t is_inter = (buffer_ptr_array[id]->candidate_ptr->type == INTER_MODE ||
                            buffer_ptr_array[id]->candidate_ptr->use_intrabc)
            ? EB_TRUE
            : EB_FALSE;
        if (!is_inter)
            if (*buffer_ptr_array[id]->full_cost_ptr < context_ptr->best_intra_cost)
                context_ptr->best_intra_cost = *buffer_ptr_array[id]->full_cost_ptr;
        if (is_inter)
            if (*buffer_ptr_array[id]->full_cost_ptr < context_ptr->best_inter_cost)
                context_ptr->best_inter_cost = *buffer_ptr_array[id]->full_cost_ptr;
    }

    // prune intra count based on inter/intra cost deviation;
    // check that th != INF to avoid overflow errors
    if ((context_ptr->mds3_intra_prune_th != (uint16_t)~0) &&
        (context_ptr->best_inter_cost * (100 + context_ptr->mds3_intra_prune_th)) <
            (context_ptr->best_intra_cost * 100))
        context_ptr->md_stage_3_total_intra_count = 0;
    sort_array_index_fast_cost_ptr(
        buffer_ptr_array, sorted_candidate_index_array, fullReconCandidateCount);
}

void md_stage_0(

    PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base,
    ModeDecisionCandidate *fast_candidate_array, int32_t fast_candidate_start_index,
    int32_t fast_candidate_end_index, EbPictureBufferDesc *input_picture_ptr,
    uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
    uint32_t input_cr_origin_in_index, BlkStruct *blk_ptr, uint32_t blk_origin_index,
    uint32_t blk_chroma_origin_index, uint32_t candidate_buffer_start_index, uint32_t max_buffers,
    EbBool scratch_buffer_pesent_flag) {
    int32_t  fast_loop_cand_index;
    uint64_t luma_fast_distortion;
    uint32_t highest_cost_index;
    uint64_t highest_cost;
    uint64_t best_first_fast_cost_search_candidate_cost  = MAX_CU_COST;
    int32_t  best_first_fast_cost_search_candidate_index = INVALID_FAST_CANDIDATE_INDEX;
    EbBool   use_ssd                                     = EB_FALSE;
    uint32_t fast_lambda                                 = context_ptr->hbd_mode_decision
        ? context_ptr->fast_lambda_md[EB_10_BIT_MD]
        : context_ptr->fast_lambda_md[EB_8_BIT_MD];
    // Set MD Staging fast_loop_core settings
    context_ptr->md_staging_skip_interpolation_search = (context_ptr->md_staging_mode ==
                                                             MD_STAGING_MODE_1 ||
                                                         context_ptr->md_staging_mode ==
                                                             MD_STAGING_MODE_2)
        ? EB_TRUE
        : context_ptr->interpolation_search_level >= IT_SEARCH_FAST_LOOP_UV_BLIND ? EB_FALSE
                                                                                  : EB_TRUE;

    context_ptr->md_staging_skip_inter_chroma_pred =
        ((context_ptr->md_staging_mode == MD_STAGING_MODE_1 ||
          context_ptr->md_staging_mode == MD_STAGING_MODE_2) &&
         context_ptr->target_class != CAND_CLASS_0 && context_ptr->target_class != CAND_CLASS_6 &&
         context_ptr->target_class != CAND_CLASS_7)
        ? EB_TRUE
        : EB_FALSE;

    context_ptr->md_staging_use_bilinear = (context_ptr->md_staging_mode == MD_STAGING_MODE_1 ||
                                            context_ptr->md_staging_mode == MD_STAGING_MODE_2)
        ? EB_TRUE
        : EB_FALSE;
    // 1st fast loop: src-to-src
    fast_loop_cand_index = fast_candidate_end_index;
    while (fast_loop_cand_index >= fast_candidate_start_index) {
        if (fast_candidate_array[fast_loop_cand_index].cand_class == context_ptr->target_class) {
            // Set the Candidate Buffer
            ModeDecisionCandidateBuffer *candidate_buffer =
                candidate_buffer_ptr_array_base[candidate_buffer_start_index];
            ModeDecisionCandidate *candidate_ptr = candidate_buffer->candidate_ptr =
                &fast_candidate_array[fast_loop_cand_index];
            // Initialize tx_depth
            candidate_buffer->candidate_ptr->tx_depth = 0;
            // Only check (src - src) candidates (Tier0 candidates)
            if (candidate_ptr->distortion_ready) {
                // Distortion
                luma_fast_distortion = candidate_ptr->me_distortion;

                // Fast Cost
                *(candidate_buffer->fast_cost_ptr) =
                    av1_product_fast_cost_func_table[candidate_ptr->type](
                        blk_ptr,
                        candidate_buffer->candidate_ptr,
                        blk_ptr->qp,
                        luma_fast_distortion,
                        0,
                        fast_lambda,
                        0,
                        pcs_ptr,
                        &(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                              .ed_ref_mv_stack[candidate_ptr->ref_frame_type][0]),
                        context_ptr->blk_geom,
                        context_ptr->blk_origin_y >> MI_SIZE_LOG2,
                        context_ptr->blk_origin_x >> MI_SIZE_LOG2,
                        context_ptr->md_enable_inter_intra,
                        context_ptr->full_cost_shut_fast_rate_flag,
                        1,
                        context_ptr->intra_luma_left_mode,
                        context_ptr->intra_luma_top_mode);

                // Keep track of the candidate index of the best  (src - src) candidate
                if (*(candidate_buffer->fast_cost_ptr) <=
                    best_first_fast_cost_search_candidate_cost) {
                    best_first_fast_cost_search_candidate_index = fast_loop_cand_index;
                    best_first_fast_cost_search_candidate_cost = *(candidate_buffer->fast_cost_ptr);
                }

                // Initialize Fast Cost - to do not interact with the second Fast-Cost Search
                *(candidate_buffer->fast_cost_ptr) = MAX_CU_COST;
            }
        }
        --fast_loop_cand_index;
    }

    // 2nd fast loop: src-to-recon
    highest_cost_index   = candidate_buffer_start_index;
    fast_loop_cand_index = fast_candidate_end_index;
    while (fast_loop_cand_index >= fast_candidate_start_index) {
        if (fast_candidate_array[fast_loop_cand_index].cand_class == context_ptr->target_class) {
            ModeDecisionCandidateBuffer *candidate_buffer =
                candidate_buffer_ptr_array_base[highest_cost_index];
            ModeDecisionCandidate *candidate_ptr = candidate_buffer->candidate_ptr =
                &fast_candidate_array[fast_loop_cand_index];
            // Initialize tx_depth
            candidate_buffer->candidate_ptr->tx_depth = 0;
            if (!candidate_ptr->distortion_ready ||
                fast_loop_cand_index == best_first_fast_cost_search_candidate_index) {
                // Prediction
                fast_loop_core(candidate_buffer,
                               pcs_ptr,
                               context_ptr,
                               input_picture_ptr,
                               input_origin_index,
                               input_cb_origin_in_index,
                               input_cr_origin_in_index,
                               blk_ptr,
                               blk_origin_index,
                               blk_chroma_origin_index,
                               use_ssd);
            }

            // Find the buffer with the highest cost
            if (fast_loop_cand_index || scratch_buffer_pesent_flag) {
                // max_cost is volatile to prevent the compiler from loading 0xFFFFFFFFFFFFFF
                //   as a const at the early-out. Loading a large constant on intel x64 processors
                //   clogs the i-cache/intstruction decode. This still reloads the variable from
                //   the stack each pass, so a better solution would be to register the variable,
                //   but this might require asm.
                volatile uint64_t max_cost           = MAX_CU_COST;
                const uint64_t *  fast_cost_array    = context_ptr->fast_cost_array;
                const uint32_t    buffer_index_start = candidate_buffer_start_index;
                const uint32_t    buffer_index_end   = buffer_index_start + max_buffers;
                uint32_t          buffer_index;

                highest_cost_index = buffer_index_start;
                buffer_index       = buffer_index_start + 1;

                do {
                    highest_cost = fast_cost_array[highest_cost_index];
                    if (highest_cost == max_cost)
                        break;

                    if (fast_cost_array[buffer_index] > highest_cost)
                        highest_cost_index = buffer_index;
                } while (++buffer_index < buffer_index_end);
            }
        }
        --fast_loop_cand_index;
    }

    // Set the cost of the scratch canidate to max to get discarded @ the sorting phase
    *(candidate_buffer_ptr_array_base[highest_cost_index]->fast_cost_ptr) =
        (scratch_buffer_pesent_flag)
        ? MAX_CU_COST
        : *(candidate_buffer_ptr_array_base[highest_cost_index]->fast_cost_ptr);
}
void md_full_pel_search(ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                        EbPictureBufferDesc *ref_pic, uint32_t input_origin_index, EbBool use_ssd,
                        int16_t mvx, int16_t mvy, int16_t search_position_start_x,
                        int16_t search_position_end_x, int16_t search_position_start_y,
                        int16_t search_position_end_y, int16_t search_step, int16_t *best_mvx,
                        int16_t *best_mvy, uint32_t *best_distortion) {
    uint8_t hbd_mode_decision = context_ptr->hbd_mode_decision == EB_DUAL_BIT_MD
        ? EB_8_BIT_MD
        : context_ptr->hbd_mode_decision;
    uint32_t                     distortion;
    ModeDecisionCandidateBuffer *candidate_buffer = &(
        context_ptr->candidate_buffer_ptr_array[0][0]);
    candidate_buffer->candidate_ptr = &(context_ptr->fast_candidate_array[0]);

    // Search area adjustment
    if ((context_ptr->blk_origin_x + (mvx >> 3) + search_position_start_x) <
        (-ref_pic->origin_x + 1))
        search_position_start_x = (-ref_pic->origin_x + 1) -
            (context_ptr->blk_origin_x + (mvx >> 3));

    if ((context_ptr->blk_origin_x + context_ptr->blk_geom->bwidth + (mvx >> 3) +
         search_position_end_x) > (ref_pic->origin_x + ref_pic->max_width - 1))
        search_position_end_x = (ref_pic->origin_x + ref_pic->max_width - 1) -
            (context_ptr->blk_origin_x + context_ptr->blk_geom->bwidth + (mvx >> 3));

    if ((context_ptr->blk_origin_y + (mvy >> 3) + search_position_start_y) <
        (-ref_pic->origin_y + 1))
        search_position_start_y = (-ref_pic->origin_y + 1) -
            (context_ptr->blk_origin_y + (mvy >> 3));

    if ((context_ptr->blk_origin_y + context_ptr->blk_geom->bheight + (mvy >> 3) +
         search_position_end_y) > (ref_pic->origin_y + ref_pic->max_height - 1))
        search_position_end_y = (ref_pic->origin_y + ref_pic->max_height - 1) -
            (context_ptr->blk_origin_y + context_ptr->blk_geom->bheight + (mvy >> 3));

#if RESTRUCTURE_SAD
    if (use_ssd) {
#if SWITCH_XY_LOOPS_PME_SAD_SSD
        for (int32_t refinement_pos_y = search_position_start_y;
             refinement_pos_y <= search_position_end_y;
             ++refinement_pos_y) {
            for (int32_t refinement_pos_x = search_position_start_x;
                 refinement_pos_x <= search_position_end_x;
                 ++refinement_pos_x) {
#else
        for (int32_t refinement_pos_x = search_position_start_x;
             refinement_pos_x <= search_position_end_x;
             ++refinement_pos_x) {
            for (int32_t refinement_pos_y = search_position_start_y;
                 refinement_pos_y <= search_position_end_y;
                 ++refinement_pos_y) {
#endif
                int32_t ref_origin_index = ref_pic->origin_x +
                    (context_ptr->blk_origin_x + (mvx >> 3) + refinement_pos_x) +
                    (context_ptr->blk_origin_y + (mvy >> 3) + ref_pic->origin_y +
                     refinement_pos_y) *
                        ref_pic->stride_y;

                EbSpatialFullDistType spatial_full_dist_type_fun = hbd_mode_decision
                    ? full_distortion_kernel16_bits
                    : spatial_full_distortion_kernel;

                distortion = (uint32_t)spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                                  input_origin_index,
                                                                  input_picture_ptr->stride_y,
                                                                  ref_pic->buffer_y,
                                                                  ref_origin_index,
                                                                  ref_pic->stride_y,
                                                                  context_ptr->blk_geom->bwidth,
                                                                  context_ptr->blk_geom->bheight);

                if (distortion < *best_distortion) {
                    *best_mvx        = mvx + (refinement_pos_x * search_step);
                    *best_mvy        = mvy + (refinement_pos_y * search_step);
                    *best_distortion = distortion;
                }
            }
        }
    } else {
        uint32_t ref_origin_index = ref_pic->origin_x +
            (context_ptr->blk_origin_x + (mvx >> 3) + search_position_start_x) +
            (context_ptr->blk_origin_y + (mvy >> 3) + ref_pic->origin_y + search_position_start_y) *
                ref_pic->stride_y;
        assert((context_ptr->blk_geom->bwidth >> 3) < 17);
        uint32_t search_area_width  = search_position_end_x - search_position_start_x + 1;
        uint32_t search_area_height = search_position_end_y - search_position_start_y + 1;
        if (search_area_width & 0xfffffff8) {
            pme_sad_loop_kernel(
                input_picture_ptr->buffer_y + input_origin_index,
                input_picture_ptr->stride_y,
                ref_pic->buffer_y + ref_origin_index,
                ref_pic->stride_y,
                context_ptr->blk_geom->bheight,
                context_ptr->blk_geom->bwidth,
                best_distortion,
                best_mvx,
                best_mvy,
                search_position_start_x,
                search_position_start_y,
                (search_area_width & 0xfffffff8), //pass search_area_width multiple by 8
                search_area_height,
                search_step,
                mvx,
                mvy);
        }
        if (search_area_width & 7) {
#if SWITCH_XY_LOOPS_PME_SAD_SSD
            for (int32_t refinement_pos_y = search_position_start_y;
                 refinement_pos_y <= search_position_end_y;
                 ++refinement_pos_y) {
                int32_t refinement_pos_x = search_position_start_x +
                    (search_area_width & 0xfffffff8);
                for (; refinement_pos_x <= search_position_end_x; ++refinement_pos_x) {
#else
            int32_t refinement_pos_x = search_position_start_x + (search_area_width & 0xfffffff8);
            for (; refinement_pos_x <= search_position_end_x; ++refinement_pos_x) {
                for (int32_t refinement_pos_y = search_position_start_y;
                     refinement_pos_y <= search_position_end_y;
                     ++refinement_pos_y) {
#endif
                    ref_origin_index = ref_pic->origin_x +
                        (context_ptr->blk_origin_x + (mvx >> 3) + refinement_pos_x) +
                        (context_ptr->blk_origin_y + (mvy >> 3) + ref_pic->origin_y +
                         refinement_pos_y) *
                            ref_pic->stride_y;
                    if (hbd_mode_decision) {
                        distortion = sad_16b_kernel(
                            ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                            input_picture_ptr->stride_y,
                            ((uint16_t *)ref_pic->buffer_y) + ref_origin_index,
                            ref_pic->stride_y,
                            context_ptr->blk_geom->bheight,
                            context_ptr->blk_geom->bwidth);
                    } else {
                        distortion = nxm_sad_kernel_sub_sampled(
                            input_picture_ptr->buffer_y + input_origin_index,
                            input_picture_ptr->stride_y,
                            ref_pic->buffer_y + ref_origin_index,
                            ref_pic->stride_y,
                            context_ptr->blk_geom->bheight,
                            context_ptr->blk_geom->bwidth);
                    }

                    if (distortion < *best_distortion) {
                        *best_mvx        = mvx + (refinement_pos_x * search_step);
                        *best_mvy        = mvy + (refinement_pos_y * search_step);
                        *best_distortion = distortion;
                    }
                }
            }
        }
    }
#else
#if SWITCH_XY_LOOPS_PME_SAD_SSD
    for (int32_t refinement_pos_y = search_position_start_y;
         refinement_pos_y <= search_position_end_y;
         ++refinement_pos_y) {
        for (int32_t refinement_pos_x = search_position_start_x;
             refinement_pos_x <= search_position_end_x;
             ++refinement_pos_x) {
#else
    for (int32_t refinement_pos_x = search_position_start_x;
         refinement_pos_x <= search_position_end_x;
         ++refinement_pos_x) {
        for (int32_t refinement_pos_y = search_position_start_y;
             refinement_pos_y <= search_position_end_y;
             ++refinement_pos_y) {
#endif
            int32_t ref_origin_index = ref_pic->origin_x +
                (context_ptr->blk_origin_x + (mvx >> 3) + refinement_pos_x) +
                (context_ptr->blk_origin_y + (mvy >> 3) + ref_pic->origin_y + refinement_pos_y) *
                    ref_pic->stride_y;
            if (use_ssd) {
                EbSpatialFullDistType spatial_full_dist_type_fun = hbd_mode_decision
                    ? full_distortion_kernel16_bits
                    : spatial_full_distortion_kernel;

                distortion = (uint32_t)spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                                  input_origin_index,
                                                                  input_picture_ptr->stride_y,
                                                                  ref_pic->buffer_y,
                                                                  ref_origin_index,
                                                                  ref_pic->stride_y,
                                                                  context_ptr->blk_geom->bwidth,
                                                                  context_ptr->blk_geom->bheight);
            } else {
                assert((context_ptr->blk_geom->bwidth >> 3) < 17);

                if (hbd_mode_decision) {
                    distortion = sad_16b_kernel(
                        ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                        input_picture_ptr->stride_y,
                        ((uint16_t *)ref_pic->buffer_y) + ref_origin_index,
                        ref_pic->stride_y,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                } else {
                    distortion = nxm_sad_kernel_sub_sampled(
                        input_picture_ptr->buffer_y + input_origin_index,
                        input_picture_ptr->stride_y,
                        ref_pic->buffer_y + ref_origin_index,
                        ref_pic->stride_y,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                }
            }

            if (distortion < *best_distortion) {
                *best_mvx        = mvx + (refinement_pos_x * search_step);
                *best_mvy        = mvy + (refinement_pos_y * search_step);
                *best_distortion = distortion;
            }
        }
    }
#endif
}

void md_sub_pel_search(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                       EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index,
                       uint32_t blk_origin_index, EbBool use_ssd, uint8_t list_idx, int8_t ref_idx,
                       int16_t mvx, int16_t mvy, int16_t search_position_start_x,
                       int16_t search_position_end_x, int16_t search_position_start_y,
                       int16_t search_position_end_y, int16_t search_step, int16_t *best_mvx,
                       int16_t *best_mvy, uint32_t *best_distortion,
                       uint8_t search_central_position, uint8_t search_pattern) {
    uint8_t hbd_mode_decision = context_ptr->hbd_mode_decision == EB_DUAL_BIT_MD
        ? EB_8_BIT_MD
        : context_ptr->hbd_mode_decision;
    uint32_t                     distortion;
    ModeDecisionCandidateBuffer *candidate_buffer = &(
        context_ptr->candidate_buffer_ptr_array[0][0]);
    candidate_buffer->candidate_ptr = &(context_ptr->fast_candidate_array[0]);

    for (int32_t refinement_pos_x = search_position_start_x;
         refinement_pos_x <= search_position_end_x;
         ++refinement_pos_x) {
        for (int32_t refinement_pos_y = search_position_start_y;
             refinement_pos_y <= search_position_end_y;
             ++refinement_pos_y) {
            if (refinement_pos_x == 0 && refinement_pos_y == 0 && !search_central_position)
                continue;
            if (search_pattern == 1 && refinement_pos_x != 0 && refinement_pos_y != 0)
                continue;

            if (search_pattern == 2 && refinement_pos_y != 0)
                continue;

            if (search_pattern == 3 && refinement_pos_x != 0)
                continue;

            ModeDecisionCandidate *candidate_ptr  = candidate_buffer->candidate_ptr;
            EbPictureBufferDesc *  prediction_ptr = candidate_buffer->prediction_ptr;

            candidate_ptr->type                         = INTER_MODE;
            candidate_ptr->distortion_ready             = 0;
            candidate_ptr->use_intrabc                  = 0;
            candidate_ptr->merge_flag                   = EB_FALSE;
            candidate_ptr->prediction_direction[0]      = (EbPredDirection)list_idx;
            candidate_ptr->inter_mode                   = NEWMV;
            candidate_ptr->pred_mode                    = NEWMV;
            candidate_ptr->motion_mode                  = SIMPLE_TRANSLATION;
            candidate_ptr->is_interintra_used           = 0;
            candidate_ptr->is_compound                  = 0;
            candidate_ptr->is_new_mv                    = 1;
            candidate_ptr->drl_index                    = 0;
            candidate_ptr->ref_mv_index                 = 0;
            candidate_ptr->ref_frame_type               = svt_get_ref_frame_type(list_idx, ref_idx);
            candidate_ptr->transform_type[PLANE_TYPE_Y] = DCT_DCT;
            candidate_ptr->transform_type[PLANE_TYPE_UV] = DCT_DCT;
            candidate_ptr->motion_vector_xl0             = list_idx == 0
                ? mvx + (refinement_pos_x * search_step)
                : 0;
            candidate_ptr->motion_vector_yl0 = list_idx == 0
                ? mvy + (refinement_pos_y * search_step)
                : 0;
            candidate_ptr->motion_vector_xl1 = list_idx == 1
                ? mvx + (refinement_pos_x * search_step)
                : 0;
            candidate_ptr->motion_vector_yl1 = list_idx == 1
                ? mvy + (refinement_pos_y * search_step)
                : 0;
            candidate_ptr->ref_frame_index_l0 = list_idx == 0 ? ref_idx : -1;
            candidate_ptr->ref_frame_index_l1 = list_idx == 1 ? ref_idx : -1;
            candidate_ptr->interp_filters     = 0;

            // Prediction
            context_ptr->md_staging_skip_interpolation_search = EB_TRUE;
            context_ptr->md_staging_skip_inter_chroma_pred    = EB_TRUE;
            product_prediction_fun_table[INTER_MODE](
                hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);

            // Distortion
            if (use_ssd) {
                EbSpatialFullDistType spatial_full_dist_type_fun = hbd_mode_decision
                    ? full_distortion_kernel16_bits
                    : spatial_full_distortion_kernel;

                distortion = (uint32_t)spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                                  input_origin_index,
                                                                  input_picture_ptr->stride_y,
                                                                  prediction_ptr->buffer_y,
                                                                  (int32_t)blk_origin_index,
                                                                  prediction_ptr->stride_y,
                                                                  context_ptr->blk_geom->bwidth,
                                                                  context_ptr->blk_geom->bheight);
            } else {
                assert((context_ptr->blk_geom->bwidth >> 3) < 17);

                if (hbd_mode_decision) {
                    distortion = sad_16b_kernel(
                        ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                        input_picture_ptr->stride_y,
                        ((uint16_t *)prediction_ptr->buffer_y) + blk_origin_index,
                        prediction_ptr->stride_y,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                } else {
                    distortion = nxm_sad_kernel_sub_sampled(
                        input_picture_ptr->buffer_y + input_origin_index,
                        input_picture_ptr->stride_y,
                        prediction_ptr->buffer_y + blk_origin_index,
                        prediction_ptr->stride_y,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                }
            }
            if (distortion < *best_distortion) {
                *best_mvx        = mvx + (refinement_pos_x * search_step);
                *best_mvy        = mvy + (refinement_pos_y * search_step);
                *best_distortion = distortion;
            }
        }
    }
}

void    av1_set_ref_frame(MvReferenceFrame *rf, int8_t ref_frame_type);
uint8_t get_max_drl_index(uint8_t refmvCnt, PredictionMode mode);
uint8_t is_me_data_present(struct ModeDecisionContext *context_ptr, const MeSbResults *me_results,
                           uint8_t list_idx, uint8_t ref_idx);
// Derive me_sb_addr and me_block_offset used to access ME_MV
void derive_me_offsets(const SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                       ModeDecisionContext *context_ptr) {
    context_ptr->geom_offset_x = 0;
    context_ptr->geom_offset_y = 0;

    if (scs_ptr->seq_header.sb_size == BLOCK_128X128) {
        uint32_t me_sb_size         = scs_ptr->sb_sz;
        uint32_t me_pic_width_in_sb = (pcs_ptr->parent_pcs_ptr->aligned_width + scs_ptr->sb_sz -
                                       1) /
            me_sb_size;
        uint32_t me_sb_x           = (context_ptr->blk_origin_x / me_sb_size);
        uint32_t me_sb_y           = (context_ptr->blk_origin_y / me_sb_size);
        context_ptr->me_sb_addr    = me_sb_x + me_sb_y * me_pic_width_in_sb;
        context_ptr->geom_offset_x = (me_sb_x & 0x1) * me_sb_size;
        context_ptr->geom_offset_y = (me_sb_y & 0x1) * me_sb_size;
    } else
        context_ptr->me_sb_addr = context_ptr->sb_ptr->index;

    // Derive whether if current block would need to have offsets made
    uint32_t bwidth_offset_to_8  = (context_ptr->blk_geom->bwidth == 4) << 2;
    uint32_t bheight_offset_to_8 = (context_ptr->blk_geom->bheight == 4) << 2;

    // if there is an offset needed to set either dimension to 8
    if (bwidth_offset_to_8 || bheight_offset_to_8) {
        // Align parent block has dimensions inherited by current block, if current block has a dimension of 4
        // add 4 so the resulting block follows an 8x8 basis
        uint32_t bwidth_to_search  = context_ptr->blk_geom->bwidth + bwidth_offset_to_8;
        uint32_t bheight_to_search = context_ptr->blk_geom->bheight + bheight_offset_to_8;

        // Align parent block has origin inherited by current block
        uint32_t x_to_search = context_ptr->blk_geom->origin_x -
            (context_ptr->geom_offset_x + ((context_ptr->blk_geom->origin_x & 0x7) ? 4 : 0));
        uint32_t y_to_search = context_ptr->blk_geom->origin_y -
            (context_ptr->geom_offset_y + ((context_ptr->blk_geom->origin_y & 0x7) ? 4 : 0));

        // Search the me_block_offset to the parent block
        for (uint32_t block_index = 0;
             block_index < pcs_ptr->parent_pcs_ptr->max_number_of_pus_per_sb;
             block_index++) {
            if ((bwidth_to_search == partition_width[block_index]) &&
                (bheight_to_search == partition_height[block_index]) &&
                (x_to_search == pu_search_index_map[block_index][0]) &&
                (y_to_search == pu_search_index_map[block_index][1])) {
                context_ptr->me_block_offset = block_index;
                break;
            }
        }
    } else if (context_ptr->blk_geom->bwidth == 128 || context_ptr->blk_geom->bheight == 128) {
        context_ptr->me_block_offset = 0;
    } else {
        context_ptr->me_block_offset = get_me_info_index(
            pcs_ptr->parent_pcs_ptr->max_number_of_pus_per_sb,
            context_ptr->blk_geom,
            context_ptr->geom_offset_x,
            context_ptr->geom_offset_y);
    }
}
// Copy ME_MVs (generated @ PA) from input buffer (pcs_ptr-> .. ->me_results) to local
// MD buffers (context_ptr->sb_me_mv)
void read_refine_me_mvs(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                        EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index,
                        uint32_t blk_origin_index) {
    const SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;

    derive_me_offsets(scs_ptr, pcs_ptr, context_ptr);

    EbBool  use_ssd           = EB_TRUE;
    uint8_t hbd_mode_decision = context_ptr->hbd_mode_decision == EB_DUAL_BIT_MD
        ? EB_8_BIT_MD
        : context_ptr->hbd_mode_decision;
    input_picture_ptr = hbd_mode_decision ? pcs_ptr->input_frame16bit
                                          : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;

    //Update input origin
    input_origin_index = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);

    for (uint32_t ref_it = 0; ref_it < pcs_ptr->parent_pcs_ptr->tot_ref_frame_types; ++ref_it) {
        MvReferenceFrame ref_pair = pcs_ptr->parent_pcs_ptr->ref_frame_type_arr[ref_it];

        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_pair);

        if (rf[1] == NONE_FRAME) {
            uint8_t list_idx = get_list_idx(rf[0]);
            uint8_t ref_idx  = get_ref_frame_idx(rf[0]);

            // Get the ME MV
            const MeSbResults *me_results =
                pcs_ptr->parent_pcs_ptr->me_results[context_ptr->me_sb_addr];
            if (is_me_data_present(context_ptr, me_results, list_idx, ref_idx)) {
                int16_t me_mv_x;
                int16_t me_mv_y;
                if (list_idx == 0) {
                    me_mv_x = (me_results->me_mv_array[context_ptr->me_block_offset][ref_idx].x_mv)
                        << 1;
                    me_mv_y = (me_results->me_mv_array[context_ptr->me_block_offset][ref_idx].y_mv)
                        << 1;
                } else {
                    me_mv_x = (me_results
                                   ->me_mv_array[context_ptr->me_block_offset]
                                                [((scs_ptr->mrp_mode == 0) ? 4 : 2) + ref_idx]
                                   .x_mv)
                        << 1;
                    me_mv_y = (me_results
                                   ->me_mv_array[context_ptr->me_block_offset]
                                                [((scs_ptr->mrp_mode == 0) ? 4 : 2) + ref_idx]
                                   .y_mv)
                        << 1;
                }

                if (context_ptr->perform_me_mv_1_8_pel_ref) {
                    int16_t  best_search_mvx        = (int16_t)~0;
                    int16_t  best_search_mvy        = (int16_t)~0;
                    uint32_t best_search_distortion = (int32_t)~0;
                    uint8_t  search_pattern         = 0;
                    md_sub_pel_search(pcs_ptr,
                                      context_ptr,
                                      input_picture_ptr,
                                      input_origin_index,
                                      blk_origin_index,
                                      use_ssd,
                                      list_idx,
                                      ref_idx,
                                      me_mv_x,
                                      me_mv_y,
                                      -(REFINE_ME_MV_EIGHT_PEL_REF_WINDOW >> 1),
                                      +(REFINE_ME_MV_EIGHT_PEL_REF_WINDOW >> 1),
                                      -(REFINE_ME_MV_EIGHT_PEL_REF_WINDOW >> 1),
                                      +(REFINE_ME_MV_EIGHT_PEL_REF_WINDOW >> 1),
                                      1,
                                      &best_search_mvx,
                                      &best_search_mvy,
                                      &best_search_distortion,
                                      1,
                                      search_pattern);

                    me_mv_x = best_search_mvx;
                    me_mv_y = best_search_mvy;
                }

                context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][0] =
                    me_mv_x;
                context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][1] =
                    me_mv_y;
            }
        }
    }
}
void predictive_me_search(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                          EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index,
                          uint32_t blk_origin_index) {
    EbBool use_ssd = EB_TRUE;
#if !MR_MODE
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        use_ssd = EB_FALSE;
#endif
    uint8_t hbd_mode_decision = context_ptr->hbd_mode_decision == EB_DUAL_BIT_MD
        ? EB_8_BIT_MD
        : context_ptr->hbd_mode_decision;
    input_picture_ptr = hbd_mode_decision ? pcs_ptr->input_frame16bit
                                          : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;

    //Update input origin
    input_origin_index = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);

    // Reset valid_refined_mv
    memset(context_ptr->valid_refined_mv, 0, 8); // [2][4]

    for (uint32_t ref_it = 0; ref_it < pcs_ptr->parent_pcs_ptr->tot_ref_frame_types; ++ref_it) {
        MvReferenceFrame ref_pair = pcs_ptr->parent_pcs_ptr->ref_frame_type_arr[ref_it];

        MacroBlockD *xd = context_ptr->blk_ptr->av1xd;
        uint8_t      drli, max_drl_index;
        IntMv        nearestmv[2], nearmv[2], ref_mv[2];

        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_pair);

        // Reset search variable(s)
        uint32_t best_mvp_distortion = (int32_t)~0;
        uint32_t mvp_distortion;

        int16_t  best_search_mvx        = (int16_t)~0;
        int16_t  best_search_mvy        = (int16_t)~0;
        uint32_t best_search_distortion = (int32_t)~0;

        // Step 0: derive the MVP list; 1 nearest and up to 3 near
        int16_t mvp_x_array[PRED_ME_MAX_MVP_CANIDATES];
        int16_t mvp_y_array[PRED_ME_MAX_MVP_CANIDATES];
        int8_t  mvp_count = 0;
        if (rf[1] == NONE_FRAME) {
            MvReferenceFrame frame_type = rf[0];
            uint8_t          list_idx   = get_list_idx(rf[0]);
            uint8_t          ref_idx    = get_ref_frame_idx(rf[0]);
            if (ref_idx > 1 && context_ptr->predictive_me_level <= 5)
                continue;
            if (ref_idx > context_ptr->md_max_ref_count - 1)
                continue;
            // Get the ME MV
            const MeSbResults *me_results =
                pcs_ptr->parent_pcs_ptr->me_results[context_ptr->me_sb_addr];
            uint32_t pa_me_distortion = ~0; //any non zero value
            if (is_me_data_present(context_ptr, me_results, list_idx, ref_idx)) {
                int16_t me_mv_x;
                int16_t me_mv_y;
                if (list_idx == 0) {
                    me_mv_x = context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][REF_LIST_0]
                                                   [ref_idx][0];
                    me_mv_y = context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][REF_LIST_0]
                                                   [ref_idx][1];
                } else {
                    me_mv_x = context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][REF_LIST_1]
                                                   [ref_idx][0];
                    me_mv_y = context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][REF_LIST_1]
                                                   [ref_idx][1];
                }
                // Round-up to the closest integer the ME MV
                me_mv_x = (me_mv_x + 4) & ~0x07;
                me_mv_y = (me_mv_y + 4) & ~0x07;
                EbReferenceObject *ref_obj =
                    pcs_ptr->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
                EbPictureBufferDesc *ref_pic = hbd_mode_decision ? ref_obj->reference_picture16bit
                                                                 : ref_obj->reference_picture;

                // -------
                // Use scaled references if resolution of the reference is different from that of the input
                // -------
                use_scaled_rec_refs_if_needed(pcs_ptr, input_picture_ptr, ref_obj, &ref_pic);

                int32_t ref_origin_index = ref_pic->origin_x +
                    (context_ptr->blk_origin_x + (me_mv_x >> 3)) +
                    (context_ptr->blk_origin_y + (me_mv_y >> 3) + ref_pic->origin_y) *
                        ref_pic->stride_y;
                if (use_ssd) {
                    EbSpatialFullDistType spatial_full_dist_type_fun = hbd_mode_decision
                        ? full_distortion_kernel16_bits
                        : spatial_full_distortion_kernel;

                    pa_me_distortion = (uint32_t)spatial_full_dist_type_fun(
                        input_picture_ptr->buffer_y,
                        input_origin_index,
                        input_picture_ptr->stride_y,
                        ref_pic->buffer_y,
                        ref_origin_index,
                        ref_pic->stride_y,
                        context_ptr->blk_geom->bwidth,
                        context_ptr->blk_geom->bheight);
                } else {
                    assert((context_ptr->blk_geom->bwidth >> 3) < 17);

                    if (hbd_mode_decision) {
                        pa_me_distortion = sad_16b_kernel(
                            ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                            input_picture_ptr->stride_y,
                            ((uint16_t *)ref_pic->buffer_y) + ref_origin_index,
                            ref_pic->stride_y,
                            context_ptr->blk_geom->bheight,
                            context_ptr->blk_geom->bwidth);
                    } else {
                        pa_me_distortion = nxm_sad_kernel_sub_sampled(
                            input_picture_ptr->buffer_y + input_origin_index,
                            input_picture_ptr->stride_y,
                            ref_pic->buffer_y + ref_origin_index,
                            ref_pic->stride_y,
                            context_ptr->blk_geom->bheight,
                            context_ptr->blk_geom->bwidth);
                    }
                }
            }
            if (pa_me_distortion != 0 || context_ptr->predictive_me_level >= 5) {
                //NEAREST
                mvp_x_array[mvp_count] = (context_ptr
                                              ->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                                              .ref_mvs[frame_type][0]
                                              .as_mv.col +
                                          4) &
                    ~0x07;
                mvp_y_array[mvp_count] = (context_ptr
                                              ->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                                              .ref_mvs[frame_type][0]
                                              .as_mv.row +
                                          4) &
                    ~0x07;

                mvp_count++;

                //NEAR
                max_drl_index = get_max_drl_index(xd->ref_mv_count[frame_type], NEARMV);

                for (drli = 0; drli < max_drl_index; drli++) {
                    get_av1_mv_pred_drl(context_ptr,
                                        context_ptr->blk_ptr,
                                        frame_type,
                                        0,
                                        NEARMV,
                                        drli,
                                        nearestmv,
                                        nearmv,
                                        ref_mv);

                    if (((nearmv[0].as_mv.col + 4) & ~0x07) != mvp_x_array[0] &&
                        ((nearmv[0].as_mv.row + 4) & ~0x07) != mvp_y_array[0]) {
                        mvp_x_array[mvp_count] = (nearmv[0].as_mv.col + 4) & ~0x07;
                        mvp_y_array[mvp_count] = (nearmv[0].as_mv.row + 4) & ~0x07;
                        mvp_count++;
                    }
                }
                // Step 1: derive the best MVP in term of distortion
                int16_t best_mvp_x = 0;
                int16_t best_mvp_y = 0;

                EbReferenceObject *ref_obj =
                    pcs_ptr->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
                EbPictureBufferDesc *ref_pic = hbd_mode_decision ? ref_obj->reference_picture16bit
                                                                 : ref_obj->reference_picture;

                // -------
                // Use scaled references if resolution of the reference is different from that of the input
                // -------
                use_scaled_rec_refs_if_needed(pcs_ptr, input_picture_ptr, ref_obj, &ref_pic);

                for (int8_t mvp_index = 0; mvp_index < mvp_count; mvp_index++) {
                    // MVP Distortion
                    // Skip the pred_me at the boundary
                    if (context_ptr->blk_origin_x + (mvp_x_array[mvp_index] >> 3) +
                                context_ptr->blk_geom->bwidth >
                            ref_pic->max_width + ref_pic->origin_x ||
                        context_ptr->blk_origin_y + (mvp_y_array[mvp_index] >> 3) +
                                context_ptr->blk_geom->bheight >
                            ref_pic->max_height + ref_pic->origin_y ||
                        context_ptr->blk_origin_x + (mvp_x_array[mvp_index] >> 3) <
                            -ref_pic->origin_x ||
                        context_ptr->blk_origin_y + (mvp_y_array[mvp_index] >> 3) <
                            -ref_pic->origin_y)
                        continue;
                    int32_t ref_origin_index = ref_pic->origin_x +
                        (context_ptr->blk_origin_x + (mvp_x_array[mvp_index] >> 3)) +
                        (context_ptr->blk_origin_y + (mvp_y_array[mvp_index] >> 3) +
                         ref_pic->origin_y) *
                            ref_pic->stride_y;
                    if (use_ssd) {
                        EbSpatialFullDistType spatial_full_dist_type_fun = hbd_mode_decision
                            ? full_distortion_kernel16_bits
                            : spatial_full_distortion_kernel;

                        mvp_distortion = (uint32_t)spatial_full_dist_type_fun(
                            input_picture_ptr->buffer_y,
                            input_origin_index,
                            input_picture_ptr->stride_y,
                            ref_pic->buffer_y,
                            ref_origin_index,
                            ref_pic->stride_y,
                            context_ptr->blk_geom->bwidth,
                            context_ptr->blk_geom->bheight);
                    } else {
                        assert((context_ptr->blk_geom->bwidth >> 3) < 17);

                        if (hbd_mode_decision) {
                            mvp_distortion = sad_16b_kernel(
                                ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                                input_picture_ptr->stride_y,
                                ((uint16_t *)ref_pic->buffer_y) + ref_origin_index,
                                ref_pic->stride_y,
                                context_ptr->blk_geom->bheight,
                                context_ptr->blk_geom->bwidth);
                        } else {
                            mvp_distortion = nxm_sad_kernel_sub_sampled(
                                input_picture_ptr->buffer_y + input_origin_index,
                                input_picture_ptr->stride_y,
                                ref_pic->buffer_y + ref_origin_index,
                                ref_pic->stride_y,
                                context_ptr->blk_geom->bheight,
                                context_ptr->blk_geom->bwidth);
                        }
                    }

                    if (mvp_distortion < best_mvp_distortion) {
                        best_mvp_distortion = mvp_distortion;
                        best_mvp_x          = mvp_x_array[mvp_index];
                        best_mvp_y          = mvp_y_array[mvp_index];
                    }
                }

                // Step 2: perform full pel search around the best MVP
                best_mvp_x = (best_mvp_x + 4) & ~0x07;
                best_mvp_y = (best_mvp_y + 4) & ~0x07;

                md_full_pel_search(context_ptr,
                                   input_picture_ptr,
                                   ref_pic,
                                   input_origin_index,
#if ENABLE_PME_SAD
                                   0,
#else
                                   use_ssd,
#endif
                                   best_mvp_x,
                                   best_mvp_y,
                                   -(context_ptr->pred_me_full_pel_search_width >> 1),
                                   +(context_ptr->pred_me_full_pel_search_width >> 1),
                                   -(context_ptr->pred_me_full_pel_search_height >> 1),
                                   +(context_ptr->pred_me_full_pel_search_height >> 1),
                                   8,
                                   &best_search_mvx,
                                   &best_search_mvy,
                                   &best_search_distortion);

                EbBool exit_predictive_me_sub_pel;

                if (pa_me_distortion == 0)
                    exit_predictive_me_sub_pel = EB_TRUE;
                else if (best_search_distortion <= pa_me_distortion)
                    exit_predictive_me_sub_pel = EB_FALSE;
                else {
                    exit_predictive_me_sub_pel = ((((best_search_distortion - pa_me_distortion) *
                                                    100) /
                                                   pa_me_distortion) < PRED_ME_DEVIATION_TH)
                        ? EB_FALSE
                        : EB_TRUE;
                }

                if (exit_predictive_me_sub_pel == EB_FALSE ||
                    context_ptr->predictive_me_level >= 5) {
                    if (context_ptr->predictive_me_level >= 2) {
                        uint8_t search_pattern;
                        // 0: all possible position(s): horizontal, vertical, diagonal
                        // 1: horizontal, vertical
                        // 2: horizontal only
                        // 3: vertical only

                        // Step 3: perform half pel search around the best full pel position
                        search_pattern = (context_ptr->predictive_me_level >= 4) ? 0 : 1;

                        md_sub_pel_search(pcs_ptr,
                                          context_ptr,
                                          input_picture_ptr,
                                          input_origin_index,
                                          blk_origin_index,
                                          use_ssd,
                                          list_idx,
                                          ref_idx,
                                          best_search_mvx,
                                          best_search_mvy,
                                          -(PRED_ME_HALF_PEL_REF_WINDOW >> 1),
                                          +(PRED_ME_HALF_PEL_REF_WINDOW >> 1),
                                          -(PRED_ME_HALF_PEL_REF_WINDOW >> 1),
                                          +(PRED_ME_HALF_PEL_REF_WINDOW >> 1),
                                          4,
                                          &best_search_mvx,
                                          &best_search_mvy,
                                          &best_search_distortion,
                                          0,
                                          search_pattern);

                        if (context_ptr->predictive_me_level == 3) {
                            if ((best_search_mvx & 0x07) != 0 || (best_search_mvy & 0x07) != 0) {
                                if ((best_search_mvx & 0x07) == 0)
                                    search_pattern = 2;
                                else // if(best_search_mvy & 0x07 == 0)
                                    search_pattern = 3;

                                md_sub_pel_search(pcs_ptr,
                                                  context_ptr,
                                                  input_picture_ptr,
                                                  input_origin_index,
                                                  blk_origin_index,
                                                  use_ssd,
                                                  list_idx,
                                                  ref_idx,
                                                  best_search_mvx,
                                                  best_search_mvy,
                                                  -(PRED_ME_HALF_PEL_REF_WINDOW >> 1),
                                                  +(PRED_ME_HALF_PEL_REF_WINDOW >> 1),
                                                  -(PRED_ME_HALF_PEL_REF_WINDOW >> 1),
                                                  +(PRED_ME_HALF_PEL_REF_WINDOW >> 1),
                                                  4,
                                                  &best_search_mvx,
                                                  &best_search_mvy,
                                                  &best_search_distortion,
                                                  0,
                                                  search_pattern);
                            }
                        }

                        // Step 4: perform quarter pel search around the best half pel position
                        search_pattern = (context_ptr->predictive_me_level >= 4) ? 0 : 1;
                        md_sub_pel_search(pcs_ptr,
                                          context_ptr,
                                          input_picture_ptr,
                                          input_origin_index,
                                          blk_origin_index,
                                          use_ssd,
                                          list_idx,
                                          ref_idx,
                                          best_search_mvx,
                                          best_search_mvy,
                                          -(PRED_ME_QUARTER_PEL_REF_WINDOW >> 1),
                                          +(PRED_ME_QUARTER_PEL_REF_WINDOW >> 1),
                                          -(PRED_ME_QUARTER_PEL_REF_WINDOW >> 1),
                                          +(PRED_ME_QUARTER_PEL_REF_WINDOW >> 1),
                                          2,
                                          &best_search_mvx,
                                          &best_search_mvy,
                                          &best_search_distortion,
                                          0,
                                          search_pattern);

                        if (context_ptr->predictive_me_level == 3) {
                            if ((best_search_mvx & 0x03) != 0 || (best_search_mvy & 0x03) != 0) {
                                if ((best_search_mvx & 0x03) == 0)
                                    search_pattern = 2;
                                else // if(best_search_mvy & 0x03 == 0)
                                    search_pattern = 3;

                                md_sub_pel_search(pcs_ptr,
                                                  context_ptr,
                                                  input_picture_ptr,
                                                  input_origin_index,
                                                  blk_origin_index,
                                                  use_ssd,
                                                  list_idx,
                                                  ref_idx,
                                                  best_search_mvx,
                                                  best_search_mvy,
                                                  -(PRED_ME_QUARTER_PEL_REF_WINDOW >> 1),
                                                  +(PRED_ME_QUARTER_PEL_REF_WINDOW >> 1),
                                                  -(PRED_ME_QUARTER_PEL_REF_WINDOW >> 1),
                                                  +(PRED_ME_QUARTER_PEL_REF_WINDOW >> 1),
                                                  2,
                                                  &best_search_mvx,
                                                  &best_search_mvy,
                                                  &best_search_distortion,
                                                  0,
                                                  search_pattern);
                            }
                        }
                    }
                    // Step 5: perform eigh pel search around the best quarter pel position
                    if (pcs_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv) {
                        uint8_t search_pattern = 0;
                        md_sub_pel_search(pcs_ptr,
                                          context_ptr,
                                          input_picture_ptr,
                                          input_origin_index,
                                          blk_origin_index,
                                          use_ssd,
                                          list_idx,
                                          ref_idx,
                                          best_search_mvx,
                                          best_search_mvy,
                                          -(PRED_ME_EIGHT_PEL_REF_WINDOW >> 1),
                                          +(PRED_ME_EIGHT_PEL_REF_WINDOW >> 1),
                                          -(PRED_ME_EIGHT_PEL_REF_WINDOW >> 1),
                                          +(PRED_ME_EIGHT_PEL_REF_WINDOW >> 1),
                                          1,
                                          &best_search_mvx,
                                          &best_search_mvy,
                                          &best_search_distortion,
                                          0,
                                          search_pattern);
                    }
                    context_ptr->best_spatial_pred_mv[list_idx][ref_idx][0] = best_search_mvx;
                    context_ptr->best_spatial_pred_mv[list_idx][ref_idx][1] = best_search_mvy;
                    context_ptr->valid_refined_mv[list_idx][ref_idx]        = 1;
                }
            }
        }
    }
}
void av1_cost_calc_cfl(PictureControlSet *pcs_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
                       SuperBlock *sb_ptr, ModeDecisionContext *context_ptr,
                       uint32_t component_mask, EbPictureBufferDesc *input_picture_ptr,
                       uint32_t input_cb_origin_in_index, uint32_t blk_chroma_origin_index,
                       uint64_t full_distortion[DIST_CALC_TOTAL], uint64_t *coeff_bits,
                       EbBool check_dc) {
    ModeDecisionCandidate *candidate_ptr = candidate_buffer->candidate_ptr;
    uint32_t               count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];
    uint64_t               cb_full_distortion[DIST_CALC_TOTAL];
    uint64_t               cr_full_distortion[DIST_CALC_TOTAL];
    uint64_t               cb_coeff_bits = 0;
    uint64_t               cr_coeff_bits = 0;
    uint32_t               chroma_width  = context_ptr->blk_geom->bwidth_uv;
    uint32_t               chroma_height = context_ptr->blk_geom->bheight_uv;
    // FullLoop and TU search
    int32_t  alpha_q3;
    uint16_t cb_qp = context_ptr->qp;
    uint16_t cr_qp = context_ptr->qp;

    full_distortion[DIST_CALC_RESIDUAL]   = 0;
    full_distortion[DIST_CALC_PREDICTION] = 0;
    *coeff_bits                           = 0;

    // Loop over alphas and find the best
    if (component_mask == COMPONENT_CHROMA_CB || component_mask == COMPONENT_CHROMA ||
        component_mask == COMPONENT_ALL) {
        cb_full_distortion[DIST_CALC_RESIDUAL]   = 0;
        cr_full_distortion[DIST_CALC_RESIDUAL]   = 0;
        cb_full_distortion[DIST_CALC_PREDICTION] = 0;
        cr_full_distortion[DIST_CALC_PREDICTION] = 0;
        cb_coeff_bits                            = 0;
        cr_coeff_bits                            = 0;
        alpha_q3                                 = (check_dc) ? 0
                              : cfl_idx_to_alpha(candidate_ptr->cfl_alpha_idx,
                                                 candidate_ptr->cfl_alpha_signs,
                                                 CFL_PRED_U); // once for U, once for V
        assert(chroma_width * CFL_BUF_LINE + chroma_height <= CFL_BUF_SQUARE);

        if (!context_ptr->hbd_mode_decision) {
            eb_cfl_predict_lbd(
                context_ptr->pred_buf_q3,
                &(candidate_buffer->prediction_ptr->buffer_cb[blk_chroma_origin_index]),
                candidate_buffer->prediction_ptr->stride_cb,
                &(candidate_buffer->cfl_temp_prediction_ptr->buffer_cb[blk_chroma_origin_index]),
                candidate_buffer->cfl_temp_prediction_ptr->stride_cb,
                alpha_q3,
                8,
                chroma_width,
                chroma_height);
        } else {
            eb_cfl_predict_hbd(
                context_ptr->pred_buf_q3,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cb) + blk_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cb,
                ((uint16_t *)candidate_buffer->cfl_temp_prediction_ptr->buffer_cb) +
                    blk_chroma_origin_index,
                candidate_buffer->cfl_temp_prediction_ptr->stride_cb,
                alpha_q3,
                10,
                chroma_width,
                chroma_height);
        }

        // Cb Residual
        residual_kernel(input_picture_ptr->buffer_cb,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cb,
                        candidate_buffer->cfl_temp_prediction_ptr->buffer_cb,
                        blk_chroma_origin_index,
                        candidate_buffer->cfl_temp_prediction_ptr->stride_cb,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cb,
                        blk_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cb,
                        context_ptr->hbd_mode_decision,
                        chroma_width,
                        chroma_height);

        full_loop_r(sb_ptr,
                    candidate_buffer,
                    context_ptr,
                    input_picture_ptr,
                    pcs_ptr,
                    PICTURE_BUFFER_DESC_Cb_FLAG,
                    cb_qp,
                    cr_qp,
                    &(*count_non_zero_coeffs[1]),
                    &(*count_non_zero_coeffs[2]));

        // Create new function
        cu_full_distortion_fast_txb_mode_r(sb_ptr,
                                           candidate_buffer,
                                           context_ptr,
                                           candidate_ptr,
                                           pcs_ptr,
                                           input_picture_ptr,
                                           cb_full_distortion,
                                           cr_full_distortion,
                                           count_non_zero_coeffs,
                                           COMPONENT_CHROMA_CB,
                                           &cb_coeff_bits,
                                           &cr_coeff_bits,
                                           0);

        full_distortion[DIST_CALC_RESIDUAL] += cb_full_distortion[DIST_CALC_RESIDUAL];
        full_distortion[DIST_CALC_PREDICTION] += cb_full_distortion[DIST_CALC_PREDICTION];
        *coeff_bits += cb_coeff_bits;
    }
    if (component_mask == COMPONENT_CHROMA_CR || component_mask == COMPONENT_CHROMA ||
        component_mask == COMPONENT_ALL) {
        cb_full_distortion[DIST_CALC_RESIDUAL]   = 0;
        cr_full_distortion[DIST_CALC_RESIDUAL]   = 0;
        cb_full_distortion[DIST_CALC_PREDICTION] = 0;
        cr_full_distortion[DIST_CALC_PREDICTION] = 0;

        cb_coeff_bits = 0;
        cr_coeff_bits = 0;
        alpha_q3      = (check_dc) ? 0
                              : cfl_idx_to_alpha(candidate_ptr->cfl_alpha_idx,
                                                 candidate_ptr->cfl_alpha_signs,
                                                 CFL_PRED_V); // once for U, once for V
        assert(chroma_width * CFL_BUF_LINE + chroma_height <= CFL_BUF_SQUARE);

        if (!context_ptr->hbd_mode_decision) {
            eb_cfl_predict_lbd(
                context_ptr->pred_buf_q3,
                &(candidate_buffer->prediction_ptr->buffer_cr[blk_chroma_origin_index]),
                candidate_buffer->prediction_ptr->stride_cr,
                &(candidate_buffer->cfl_temp_prediction_ptr->buffer_cr[blk_chroma_origin_index]),
                candidate_buffer->cfl_temp_prediction_ptr->stride_cr,
                alpha_q3,
                8,
                chroma_width,
                chroma_height);
        } else {
            eb_cfl_predict_hbd(
                context_ptr->pred_buf_q3,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cr) + blk_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cr,
                ((uint16_t *)candidate_buffer->cfl_temp_prediction_ptr->buffer_cr) +
                    blk_chroma_origin_index,
                candidate_buffer->cfl_temp_prediction_ptr->stride_cr,
                alpha_q3,
                10,
                chroma_width,
                chroma_height);
        }

        // Cr Residual
        residual_kernel(input_picture_ptr->buffer_cr,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cr,
                        candidate_buffer->cfl_temp_prediction_ptr->buffer_cr,
                        blk_chroma_origin_index,
                        candidate_buffer->cfl_temp_prediction_ptr->stride_cr,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cr,
                        blk_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cr,
                        context_ptr->hbd_mode_decision,
                        chroma_width,
                        chroma_height);

        full_loop_r(sb_ptr,
                    candidate_buffer,
                    context_ptr,
                    input_picture_ptr,
                    pcs_ptr,
                    PICTURE_BUFFER_DESC_Cr_FLAG,
                    cb_qp,
                    cr_qp,
                    &(*count_non_zero_coeffs[1]),
                    &(*count_non_zero_coeffs[2]));
        candidate_ptr->v_has_coeff = *count_non_zero_coeffs[2] ? EB_TRUE : EB_FALSE;

        // Create new function
        cu_full_distortion_fast_txb_mode_r(sb_ptr,
                                           candidate_buffer,
                                           context_ptr,
                                           candidate_ptr,
                                           pcs_ptr,
                                           input_picture_ptr,
                                           cb_full_distortion,
                                           cr_full_distortion,
                                           count_non_zero_coeffs,
                                           COMPONENT_CHROMA_CR,
                                           &cb_coeff_bits,
                                           &cr_coeff_bits,
                                           0);

        full_distortion[DIST_CALC_RESIDUAL] += cr_full_distortion[DIST_CALC_RESIDUAL];
        full_distortion[DIST_CALC_PREDICTION] += cr_full_distortion[DIST_CALC_PREDICTION];
        *coeff_bits += cr_coeff_bits;
    }
}

#define PLANE_SIGN_TO_JOINT_SIGN(plane, a, b) \
    (plane == CFL_PRED_U ? a * CFL_SIGNS + b - 1 : b * CFL_SIGNS + a - 1)
/*************************Pick the best alpha for cfl mode  or Choose DC******************************************************/
void cfl_rd_pick_alpha(PictureControlSet *pcs_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
                       SuperBlock *sb_ptr, ModeDecisionContext *context_ptr,
                       EbPictureBufferDesc *input_picture_ptr, uint32_t input_cb_origin_in_index,
                       uint32_t blk_chroma_origin_index) {
    int64_t  best_rd = INT64_MAX;
    uint64_t full_distortion[DIST_CALC_TOTAL];
    uint64_t coeff_bits;
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];

    const int64_t mode_rd = RDCOST(
        full_lambda,
        (uint64_t)candidate_buffer->candidate_ptr->md_rate_estimation_ptr
            ->intra_uv_mode_fac_bits[CFL_ALLOWED][candidate_buffer->candidate_ptr->intra_luma_mode]
                                    [UV_CFL_PRED],
        0);

    int64_t best_rd_uv[CFL_JOINT_SIGNS][CFL_PRED_PLANES];
    int32_t best_c[CFL_JOINT_SIGNS][CFL_PRED_PLANES];

    for (int32_t plane = 0; plane < CFL_PRED_PLANES; plane++) {
        coeff_bits                          = 0;
        full_distortion[DIST_CALC_RESIDUAL] = 0;
        for (int32_t joint_sign = 0; joint_sign < CFL_JOINT_SIGNS; joint_sign++) {
            best_rd_uv[joint_sign][plane] = INT64_MAX;
            best_c[joint_sign][plane]     = 0;
        }
        // Collect RD stats for an alpha value of zero in this plane.
        // Skip i == CFL_SIGN_ZERO as (0, 0) is invalid.
        for (int32_t i = CFL_SIGN_NEG; i < CFL_SIGNS; i++) {
            const int32_t joint_sign = PLANE_SIGN_TO_JOINT_SIGN(plane, CFL_SIGN_ZERO, i);
            if (i == CFL_SIGN_NEG) {
                candidate_buffer->candidate_ptr->cfl_alpha_idx   = 0;
                candidate_buffer->candidate_ptr->cfl_alpha_signs = joint_sign;

                av1_cost_calc_cfl(pcs_ptr,
                                  candidate_buffer,
                                  sb_ptr,
                                  context_ptr,
                                  (plane == 0) ? COMPONENT_CHROMA_CB : COMPONENT_CHROMA_CR,
                                  input_picture_ptr,
                                  input_cb_origin_in_index,
                                  blk_chroma_origin_index,
                                  full_distortion,
                                  &coeff_bits,
                                  0);

                if (coeff_bits == INT64_MAX)
                    break;
            }

            const int32_t alpha_rate = candidate_buffer->candidate_ptr->md_rate_estimation_ptr
                                           ->cfl_alpha_fac_bits[joint_sign][plane][0];

            best_rd_uv[joint_sign][plane] = RDCOST(
                full_lambda, coeff_bits + alpha_rate, full_distortion[DIST_CALC_RESIDUAL]);
        }
    }

    int32_t best_joint_sign = -1;

    for (int32_t plane = 0; plane < CFL_PRED_PLANES; plane++) {
        for (int32_t pn_sign = CFL_SIGN_NEG; pn_sign < CFL_SIGNS; pn_sign++) {
            int32_t progress = 0;
            for (int32_t c = 0; c < CFL_ALPHABET_SIZE; c++) {
                int32_t flag = 0;
                if (c > 2 && progress < c)
                    break;
                coeff_bits                          = 0;
                full_distortion[DIST_CALC_RESIDUAL] = 0;
                for (int32_t i = 0; i < CFL_SIGNS; i++) {
                    const int32_t joint_sign = PLANE_SIGN_TO_JOINT_SIGN(plane, pn_sign, i);
                    if (i == 0) {
                        candidate_buffer->candidate_ptr->cfl_alpha_idx =
                            (c << CFL_ALPHABET_SIZE_LOG2) + c;
                        candidate_buffer->candidate_ptr->cfl_alpha_signs = joint_sign;

                        av1_cost_calc_cfl(pcs_ptr,
                                          candidate_buffer,
                                          sb_ptr,
                                          context_ptr,
                                          (plane == 0) ? COMPONENT_CHROMA_CB : COMPONENT_CHROMA_CR,
                                          input_picture_ptr,
                                          input_cb_origin_in_index,
                                          blk_chroma_origin_index,
                                          full_distortion,
                                          &coeff_bits,
                                          0);

                        if (coeff_bits == INT64_MAX)
                            break;
                    }

                    const int32_t alpha_rate =
                        candidate_buffer->candidate_ptr->md_rate_estimation_ptr
                            ->cfl_alpha_fac_bits[joint_sign][plane][c];

                    int64_t this_rd = RDCOST(
                        full_lambda, coeff_bits + alpha_rate, full_distortion[DIST_CALC_RESIDUAL]);
                    if (this_rd >= best_rd_uv[joint_sign][plane])
                        continue;
                    best_rd_uv[joint_sign][plane] = this_rd;
                    best_c[joint_sign][plane]     = c;

                    flag = 2;
                    if (best_rd_uv[joint_sign][!plane] == INT64_MAX)
                        continue;
                    this_rd += mode_rd + best_rd_uv[joint_sign][!plane];
                    if (this_rd >= best_rd)
                        continue;
                    best_rd         = this_rd;
                    best_joint_sign = joint_sign;
                }
                progress += flag;
            }
        }
    }

    // Compare with DC Chroma
    coeff_bits                          = 0;
    full_distortion[DIST_CALC_RESIDUAL] = 0;

    candidate_buffer->candidate_ptr->cfl_alpha_idx   = 0;
    candidate_buffer->candidate_ptr->cfl_alpha_signs = 0;

    const int64_t dc_mode_rd = RDCOST(
        full_lambda,
        candidate_buffer->candidate_ptr->md_rate_estimation_ptr
            ->intra_uv_mode_fac_bits[CFL_ALLOWED][candidate_buffer->candidate_ptr->intra_luma_mode]
                                    [UV_DC_PRED],
        0);

    av1_cost_calc_cfl(pcs_ptr,
                      candidate_buffer,
                      sb_ptr,
                      context_ptr,
                      COMPONENT_CHROMA,
                      input_picture_ptr,
                      input_cb_origin_in_index,
                      blk_chroma_origin_index,
                      full_distortion,
                      &coeff_bits,
                      1);

    int64_t dc_rd = RDCOST(full_lambda, coeff_bits, full_distortion[DIST_CALC_RESIDUAL]);

    dc_rd += dc_mode_rd;
    if (dc_rd <= best_rd) {
        candidate_buffer->candidate_ptr->intra_chroma_mode = UV_DC_PRED;
        candidate_buffer->candidate_ptr->cfl_alpha_idx     = 0;
        candidate_buffer->candidate_ptr->cfl_alpha_signs   = 0;
    } else {
        candidate_buffer->candidate_ptr->intra_chroma_mode = UV_CFL_PRED;
        int32_t ind                                        = 0;
        if (best_joint_sign >= 0) {
            const int32_t u = best_c[best_joint_sign][CFL_PRED_U];
            const int32_t v = best_c[best_joint_sign][CFL_PRED_V];
            ind             = (u << CFL_ALPHABET_SIZE_LOG2) + v;
        } else
            best_joint_sign = 0;
        candidate_buffer->candidate_ptr->cfl_alpha_idx   = ind;
        candidate_buffer->candidate_ptr->cfl_alpha_signs = best_joint_sign;
    }
}

// If mode is CFL:
// 1: recon the Luma
// 2: Form the pred_buf_q3
// 3: Loop over alphas and find the best or choose DC
// 4: Recalculate the residual for chroma
static void cfl_prediction(PictureControlSet *          pcs_ptr,
                           ModeDecisionCandidateBuffer *candidate_buffer, SuperBlock *sb_ptr,
                           ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                           uint32_t input_cb_origin_in_index, uint32_t blk_chroma_origin_index) {
    if (context_ptr->blk_geom->has_uv) {
        // 1: recon the Luma
        av1_perform_inverse_transform_recon_luma(context_ptr, candidate_buffer);

        uint32_t rec_luma_offset = ((context_ptr->blk_geom->origin_y >> 3) << 3) *
                candidate_buffer->recon_ptr->stride_y +
            ((context_ptr->blk_geom->origin_x >> 3) << 3);
        // 2: Form the pred_buf_q3
        uint32_t chroma_width  = context_ptr->blk_geom->bwidth_uv;
        uint32_t chroma_height = context_ptr->blk_geom->bheight_uv;

        // Down sample Luma
        if (!context_ptr->hbd_mode_decision) {
            cfl_luma_subsampling_420_lbd(
                &(context_ptr->cfl_temp_luma_recon[rec_luma_offset]),
                candidate_buffer->recon_ptr->stride_y,
                context_ptr->pred_buf_q3,
                context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bwidth
                    ? (context_ptr->blk_geom->bwidth_uv << 1)
                    : context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight_uv == context_ptr->blk_geom->bheight
                    ? (context_ptr->blk_geom->bheight_uv << 1)
                    : context_ptr->blk_geom->bheight);
        } else {
            cfl_luma_subsampling_420_hbd(
                context_ptr->cfl_temp_luma_recon16bit + rec_luma_offset,
                candidate_buffer->recon_ptr->stride_y,
                context_ptr->pred_buf_q3,
                context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bwidth
                    ? (context_ptr->blk_geom->bwidth_uv << 1)
                    : context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight_uv == context_ptr->blk_geom->bheight
                    ? (context_ptr->blk_geom->bheight_uv << 1)
                    : context_ptr->blk_geom->bheight);
        }
        int32_t round_offset = chroma_width * chroma_height / 2;

        eb_subtract_average(context_ptr->pred_buf_q3,
                            chroma_width,
                            chroma_height,
                            round_offset,
                            eb_log2f(chroma_width) + eb_log2f(chroma_height));

        // 3: Loop over alphas and find the best or choose DC
        cfl_rd_pick_alpha(pcs_ptr,
                          candidate_buffer,
                          sb_ptr,
                          context_ptr,
                          input_picture_ptr,
                          input_cb_origin_in_index,
                          blk_chroma_origin_index);

        if (candidate_buffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {
            // 4: Recalculate the prediction and the residual
            int32_t alpha_q3_cb = cfl_idx_to_alpha(candidate_buffer->candidate_ptr->cfl_alpha_idx,
                                                   candidate_buffer->candidate_ptr->cfl_alpha_signs,
                                                   CFL_PRED_U);
            int32_t alpha_q3_cr = cfl_idx_to_alpha(candidate_buffer->candidate_ptr->cfl_alpha_idx,
                                                   candidate_buffer->candidate_ptr->cfl_alpha_signs,
                                                   CFL_PRED_V);

            assert(chroma_height * CFL_BUF_LINE + chroma_width <= CFL_BUF_SQUARE);

            if (!context_ptr->hbd_mode_decision) {
                eb_cfl_predict_lbd(
                    context_ptr->pred_buf_q3,
                    &(candidate_buffer->prediction_ptr->buffer_cb[blk_chroma_origin_index]),
                    candidate_buffer->prediction_ptr->stride_cb,
                    &(candidate_buffer->prediction_ptr->buffer_cb[blk_chroma_origin_index]),
                    candidate_buffer->prediction_ptr->stride_cb,
                    alpha_q3_cb,
                    8,
                    chroma_width,
                    chroma_height);

                eb_cfl_predict_lbd(
                    context_ptr->pred_buf_q3,
                    &(candidate_buffer->prediction_ptr->buffer_cr[blk_chroma_origin_index]),
                    candidate_buffer->prediction_ptr->stride_cr,
                    &(candidate_buffer->prediction_ptr->buffer_cr[blk_chroma_origin_index]),
                    candidate_buffer->prediction_ptr->stride_cr,
                    alpha_q3_cr,
                    8,
                    chroma_width,
                    chroma_height);
            } else {
                eb_cfl_predict_hbd(context_ptr->pred_buf_q3,
                                   ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cb) +
                                       blk_chroma_origin_index,
                                   candidate_buffer->prediction_ptr->stride_cb,
                                   ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cb) +
                                       blk_chroma_origin_index,
                                   candidate_buffer->prediction_ptr->stride_cb,
                                   alpha_q3_cb,
                                   10,
                                   chroma_width,
                                   chroma_height);

                eb_cfl_predict_hbd(context_ptr->pred_buf_q3,
                                   ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cr) +
                                       blk_chroma_origin_index,
                                   candidate_buffer->prediction_ptr->stride_cr,
                                   ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cr) +
                                       blk_chroma_origin_index,
                                   candidate_buffer->prediction_ptr->stride_cr,
                                   alpha_q3_cr,
                                   10,
                                   chroma_width,
                                   chroma_height);
            }

            // Cb Residual
            residual_kernel(input_picture_ptr->buffer_cb,
                            input_cb_origin_in_index,
                            input_picture_ptr->stride_cb,
                            candidate_buffer->prediction_ptr->buffer_cb,
                            blk_chroma_origin_index,
                            candidate_buffer->prediction_ptr->stride_cb,
                            (int16_t *)candidate_buffer->residual_ptr->buffer_cb,
                            blk_chroma_origin_index,
                            candidate_buffer->residual_ptr->stride_cb,
                            context_ptr->hbd_mode_decision,
                            context_ptr->blk_geom->bwidth_uv,
                            context_ptr->blk_geom->bheight_uv);

            // Cr Residual
            residual_kernel(input_picture_ptr->buffer_cr,
                            input_cb_origin_in_index,
                            input_picture_ptr->stride_cr,
                            candidate_buffer->prediction_ptr->buffer_cr,
                            blk_chroma_origin_index,
                            candidate_buffer->prediction_ptr->stride_cr,
                            (int16_t *)candidate_buffer->residual_ptr->buffer_cr,
                            blk_chroma_origin_index,
                            candidate_buffer->residual_ptr->stride_cr,
                            context_ptr->hbd_mode_decision,
                            context_ptr->blk_geom->bwidth_uv,
                            context_ptr->blk_geom->bheight_uv);
        } else {
            // Alphas = 0, Preds are the same as DC. Switch to DC mode
            candidate_buffer->candidate_ptr->intra_chroma_mode = UV_DC_PRED;
        }
    }
}
uint8_t get_skip_tx_search_flag(int32_t sq_size, uint64_t ref_fast_cost, uint64_t cu_cost,
                                uint64_t weight) {
    //NM: Skip tx search when the fast cost of the current mode candidate is substansially
    // Larger than the best fast_cost (
    uint8_t tx_search_skip_flag = cu_cost >= ((ref_fast_cost * weight) / 100) ? 1 : 0;
    tx_search_skip_flag         = sq_size >= 128 ? 1 : tx_search_skip_flag;
    return tx_search_skip_flag;
}

static INLINE TxType av1_get_tx_type(BlockSize sb_type, int32_t is_inter, PredictionMode pred_mode,
                                     UvPredictionMode pred_mode_uv, PlaneType plane_type,
                                     const MacroBlockD *xd, int32_t blk_row, int32_t blk_col,
                                     TxSize tx_size, int32_t reduced_tx_set) {
    UNUSED(sb_type);
    UNUSED(*xd);
    UNUSED(blk_row);
    UNUSED(blk_col);

    // block_size  sb_type = BLOCK_8X8;

    MbModeInfo mbmi;
    mbmi.block_mi.mode    = pred_mode;
    mbmi.block_mi.uv_mode = pred_mode_uv;

    // const MbModeInfo *const mbmi = xd->mi[0];
    // const struct MacroblockdPlane *const pd = &xd->plane[plane_type];
    const TxSetType tx_set_type =
        /*av1_*/ get_ext_tx_set_type(tx_size, is_inter, reduced_tx_set);

    TxType tx_type = DCT_DCT;
    if (/*xd->lossless[mbmi->segment_id] ||*/ txsize_sqr_up_map[tx_size] > TX_32X32)
        tx_type = DCT_DCT;
    else {
        if (plane_type == PLANE_TYPE_Y) {
            //const int32_t txk_type_idx =
            //    av1_get_txk_type_index(/*mbmi->*/sb_type, blk_row, blk_col);
            //tx_type = mbmi->txk_type[txk_type_idx];
        } else if (is_inter /*is_inter_block(mbmi)*/) {
            // scale back to y plane's coordinate
            //blk_row <<= pd->subsampling_y;
            //blk_col <<= pd->subsampling_x;
            //const int32_t txk_type_idx =
            //    av1_get_txk_type_index(mbmi->sb_type, blk_row, blk_col);
            //tx_type = mbmi->txk_type[txk_type_idx];
        } else {
            // In intra mode, uv planes don't share the same prediction mode as y
            // plane, so the tx_type should not be shared
            tx_type = intra_mode_to_tx_type(&mbmi.block_mi, PLANE_TYPE_UV);
        }
    }
    ASSERT(tx_type < TX_TYPES);
    if (!av1_ext_tx_used[tx_set_type][tx_type])
        return DCT_DCT;
    return tx_type;
}

void check_best_indepedant_cfl(PictureControlSet *pcs_ptr, EbPictureBufferDesc *input_picture_ptr,
                               ModeDecisionContext *context_ptr, uint32_t input_cb_origin_in_index,
                               uint32_t                     blk_chroma_origin_index,
                               ModeDecisionCandidateBuffer *candidate_buffer, uint8_t cb_qp,
                               uint8_t cr_qp, uint64_t *cb_full_distortion,
                               uint64_t *cr_full_distortion, uint64_t *cb_coeff_bits,
                               uint64_t *cr_coeff_bits) {
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];
    if (candidate_buffer->candidate_ptr->filter_intra_mode != FILTER_INTRA_MODES)
        assert(candidate_buffer->candidate_ptr->intra_luma_mode == DC_PRED);
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    // cfl cost
    uint64_t chroma_rate = 0;
    if (candidate_buffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {
        chroma_rate +=
            candidate_buffer->candidate_ptr->md_rate_estimation_ptr
                ->cfl_alpha_fac_bits[candidate_buffer->candidate_ptr->cfl_alpha_signs][CFL_PRED_U]
                                    [CFL_IDX_U(candidate_buffer->candidate_ptr->cfl_alpha_idx)] +
            candidate_buffer->candidate_ptr->md_rate_estimation_ptr
                ->cfl_alpha_fac_bits[candidate_buffer->candidate_ptr->cfl_alpha_signs][CFL_PRED_V]
                                    [CFL_IDX_V(candidate_buffer->candidate_ptr->cfl_alpha_idx)];

        chroma_rate +=
            (uint64_t)candidate_buffer->candidate_ptr->md_rate_estimation_ptr
                ->intra_uv_mode_fac_bits[CFL_ALLOWED][candidate_buffer->candidate_ptr
                                                          ->intra_luma_mode][UV_CFL_PRED];
        chroma_rate -= (uint64_t)candidate_buffer->candidate_ptr->md_rate_estimation_ptr
                           ->intra_uv_mode_fac_bits[CFL_ALLOWED][candidate_buffer->candidate_ptr
                                                                     ->intra_luma_mode][UV_DC_PRED];
    } else
        chroma_rate = (uint64_t)candidate_buffer->candidate_ptr->md_rate_estimation_ptr
                          ->intra_uv_mode_fac_bits[CFL_ALLOWED][candidate_buffer->candidate_ptr
                                                                    ->intra_luma_mode][UV_DC_PRED];
    int coeff_rate = (int)(*cb_coeff_bits + *cr_coeff_bits);
    int distortion = (int)(cb_full_distortion[DIST_CALC_RESIDUAL] +
                           cr_full_distortion[DIST_CALC_RESIDUAL]);
    int rate = (int)(coeff_rate + chroma_rate + candidate_buffer->candidate_ptr->fast_luma_rate);
    uint64_t cfl_uv_cost = RDCOST(full_lambda, rate, distortion);

    // cfl vs. best independant
    if (context_ptr->best_uv_cost[candidate_buffer->candidate_ptr->intra_luma_mode]
                                 [3 + candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]] <
        cfl_uv_cost) {
        // Update the current candidate
        candidate_buffer->candidate_ptr->intra_chroma_mode =
            context_ptr->best_uv_mode[candidate_buffer->candidate_ptr->intra_luma_mode]
                                     [MAX_ANGLE_DELTA +
                                      candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]];
        candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_UV] =
            context_ptr->best_uv_angle[candidate_buffer->candidate_ptr->intra_luma_mode]
                                      [MAX_ANGLE_DELTA +
                                       candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]];
        candidate_buffer->candidate_ptr
            ->is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)(
            context_ptr->best_uv_mode[candidate_buffer->candidate_ptr->intra_luma_mode]
                                     [MAX_ANGLE_DELTA +
                                      candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]]));

        // check if candidate_buffer->candidate_ptr->fast_luma_rate = context_ptr->fast_luma_rate[candidate_buffer->candidate_ptr->intra_luma_mode];
        candidate_buffer->candidate_ptr->fast_chroma_rate =
            context_ptr
                ->fast_chroma_rate[candidate_buffer->candidate_ptr->intra_luma_mode]
                                  [MAX_ANGLE_DELTA +
                                   candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]];

        candidate_buffer->candidate_ptr->transform_type_uv = av1_get_tx_type(
            context_ptr->blk_geom->bsize,
            0,
            (PredictionMode)NULL,
            (UvPredictionMode)context_ptr
                ->best_uv_mode[candidate_buffer->candidate_ptr->intra_luma_mode]
                              [3 + candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]],
            PLANE_TYPE_UV,
            0,
            0,
            0,
            context_ptr->blk_geom->txsize_uv[0][0],
            frm_hdr->reduced_tx_set);

        // Start uv search path
        context_ptr->uv_search_path = EB_TRUE;

        memset(candidate_buffer->candidate_ptr->eob[1], 0, sizeof(uint16_t));
        memset(candidate_buffer->candidate_ptr->eob[2], 0, sizeof(uint16_t));
        candidate_buffer->candidate_ptr->u_has_coeff = 0;
        candidate_buffer->candidate_ptr->v_has_coeff = 0;
        cb_full_distortion[DIST_CALC_RESIDUAL]       = 0;
        cr_full_distortion[DIST_CALC_RESIDUAL]       = 0;
        cb_full_distortion[DIST_CALC_PREDICTION]     = 0;
        cr_full_distortion[DIST_CALC_PREDICTION]     = 0;

        *cb_coeff_bits = 0;
        *cr_coeff_bits = 0;

        uint32_t count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];
        context_ptr->md_staging_skip_inter_chroma_pred = EB_FALSE;
        product_prediction_fun_table[candidate_buffer->candidate_ptr->type](
            context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);

        // Cb Residual
        residual_kernel(input_picture_ptr->buffer_cb,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cb,
                        candidate_buffer->prediction_ptr->buffer_cb,
                        blk_chroma_origin_index,
                        candidate_buffer->prediction_ptr->stride_cb,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cb,
                        blk_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cb,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->bwidth_uv,
                        context_ptr->blk_geom->bheight_uv);

        // Cr Residual
        residual_kernel(input_picture_ptr->buffer_cr,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cr,
                        candidate_buffer->prediction_ptr->buffer_cr,
                        blk_chroma_origin_index,
                        candidate_buffer->prediction_ptr->stride_cr,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cr,
                        blk_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cr,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->bwidth_uv,
                        context_ptr->blk_geom->bheight_uv);

        full_loop_r(context_ptr->sb_ptr,
                    candidate_buffer,
                    context_ptr,
                    input_picture_ptr,
                    pcs_ptr,
                    PICTURE_BUFFER_DESC_CHROMA_MASK,
                    cb_qp,
                    cr_qp,
                    &(*count_non_zero_coeffs[1]),
                    &(*count_non_zero_coeffs[2]));

        cu_full_distortion_fast_txb_mode_r(context_ptr->sb_ptr,
                                           candidate_buffer,
                                           context_ptr,
                                           candidate_buffer->candidate_ptr,
                                           pcs_ptr,
                                           input_picture_ptr,
                                           cb_full_distortion,
                                           cr_full_distortion,
                                           count_non_zero_coeffs,
                                           COMPONENT_CHROMA,
                                           cb_coeff_bits,
                                           cr_coeff_bits,
                                           1);

        // End uv search path
        context_ptr->uv_search_path = EB_FALSE;
    }
}

// double check the usage of tx_search_luma_recon_neighbor_array16bit
EbErrorType av1_intra_luma_prediction(ModeDecisionContext *        md_context_ptr,
                                      PictureControlSet *          pcs_ptr,
                                      ModeDecisionCandidateBuffer *candidate_buffer_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    uint8_t     is_inter     = 0; // set to 0 b/c this is an intra path

    uint16_t txb_origin_x = md_context_ptr->blk_origin_x +
        md_context_ptr->blk_geom
            ->tx_org_x[is_inter][md_context_ptr->tx_depth][md_context_ptr->txb_itr] -
        md_context_ptr->blk_geom->origin_x;
    uint16_t txb_origin_y = md_context_ptr->blk_origin_y +
        md_context_ptr->blk_geom
            ->tx_org_y[is_inter][md_context_ptr->tx_depth][md_context_ptr->txb_itr] -
        md_context_ptr->blk_geom->origin_y;
    uint8_t tx_width =
        md_context_ptr->blk_geom->tx_width[md_context_ptr->tx_depth][md_context_ptr->txb_itr];
    uint8_t tx_height =
        md_context_ptr->blk_geom->tx_height[md_context_ptr->tx_depth][md_context_ptr->txb_itr];

    uint32_t mode_type_left_neighbor_index = get_neighbor_array_unit_left_index(
        md_context_ptr->mode_type_neighbor_array, txb_origin_y);
    uint32_t mode_type_top_neighbor_index = get_neighbor_array_unit_top_index(
        md_context_ptr->mode_type_neighbor_array, txb_origin_x);
    uint32_t intra_luma_mode_left_neighbor_index = get_neighbor_array_unit_left_index(
        md_context_ptr->intra_luma_mode_neighbor_array, txb_origin_y);
    uint32_t intra_luma_mode_top_neighbor_index = get_neighbor_array_unit_top_index(
        md_context_ptr->intra_luma_mode_neighbor_array, txb_origin_x);

    md_context_ptr->intra_luma_left_mode = (uint32_t)(
        (md_context_ptr->mode_type_neighbor_array->left_array[mode_type_left_neighbor_index] !=
         INTRA_MODE)
            ? DC_PRED /*EB_INTRA_DC*/
            : (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array
                  ->left_array[intra_luma_mode_left_neighbor_index]);

    md_context_ptr->intra_luma_top_mode = (uint32_t)(
        (md_context_ptr->mode_type_neighbor_array->top_array[mode_type_top_neighbor_index] !=
         INTRA_MODE)
            ? DC_PRED /*EB_INTRA_DC*/
            : (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->top_array
                  [intra_luma_mode_top_neighbor_index]); //   use DC. This seems like we could use a SB-width

    TxSize tx_size =
        md_context_ptr->blk_geom->txsize[md_context_ptr->tx_depth][md_context_ptr->txb_itr];

    PredictionMode mode;
    if (!md_context_ptr->hbd_mode_decision) {
        uint8_t top_neigh_array[64 * 2 + 1];
        uint8_t left_neigh_array[64 * 2 + 1];

        if (txb_origin_y != 0)
            memcpy(top_neigh_array + 1,
                   md_context_ptr->tx_search_luma_recon_neighbor_array->top_array + txb_origin_x,
                   tx_width * 2);
        if (txb_origin_x != 0)
            memcpy(left_neigh_array + 1,
                   md_context_ptr->tx_search_luma_recon_neighbor_array->left_array + txb_origin_y,
                   tx_height * 2);
        if (txb_origin_y != 0 && txb_origin_x != 0)
            top_neigh_array[0] = left_neigh_array[0] =
                md_context_ptr->tx_search_luma_recon_neighbor_array
                    ->top_left_array[MAX_PICTURE_HEIGHT_SIZE + txb_origin_x - txb_origin_y];

        mode = candidate_buffer_ptr->candidate_ptr->pred_mode;
        eb_av1_predict_intra_block(
            &md_context_ptr->sb_ptr->tile_info,
            !ED_STAGE,
            md_context_ptr->blk_geom,
            pcs_ptr->parent_pcs_ptr->av1_cm, //const Av1Common *cm,
            md_context_ptr->blk_geom->bwidth,
            md_context_ptr->blk_geom->bheight,
            tx_size,
            mode, //PredictionMode mode,
            candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
            candidate_buffer_ptr->candidate_ptr->palette_info.pmi.palette_size[0] > 0,
            &candidate_buffer_ptr->candidate_ptr->palette_info, //ATB MD
            candidate_buffer_ptr->candidate_ptr->filter_intra_mode,
            top_neigh_array + 1,
            left_neigh_array + 1,
            candidate_buffer_ptr->prediction_ptr, //uint8_t *dst,
            (md_context_ptr->blk_geom
                 ->tx_org_x[is_inter][md_context_ptr->tx_depth][md_context_ptr->txb_itr] -
             md_context_ptr->blk_geom->origin_x) >>
                2,
            (md_context_ptr->blk_geom
                 ->tx_org_y[is_inter][md_context_ptr->tx_depth][md_context_ptr->txb_itr] -
             md_context_ptr->blk_geom->origin_y) >>
                2,
            PLANE_TYPE_Y, //int32_t plane,
            md_context_ptr->blk_geom->bsize,
            md_context_ptr->blk_origin_x,
            md_context_ptr->blk_origin_y,
            md_context_ptr->blk_origin_x,
            md_context_ptr->blk_origin_y,
            md_context_ptr->blk_geom
                ->tx_org_x[is_inter][md_context_ptr->tx_depth]
                          [md_context_ptr->txb_itr], //uint32_t cuOrgX used only for prediction Ptr
            md_context_ptr->blk_geom
                ->tx_org_y[is_inter][md_context_ptr->tx_depth]
                          [md_context_ptr->txb_itr], //uint32_t cuOrgY used only for prediction Ptr
            pcs_ptr->mi_grid_base,
            &((SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr)->seq_header);
    } else {
        uint16_t top_neigh_array[64 * 2 + 1];
        uint16_t left_neigh_array[64 * 2 + 1];

        if (txb_origin_y != 0)
            memcpy(
                top_neigh_array + 1,
                (uint16_t *)(md_context_ptr->tx_search_luma_recon_neighbor_array16bit->top_array) +
                    txb_origin_x,
                sizeof(uint16_t) * tx_width * 2);
        if (txb_origin_x != 0)
            memcpy(
                left_neigh_array + 1,
                (uint16_t *)(md_context_ptr->tx_search_luma_recon_neighbor_array16bit->left_array) +
                    txb_origin_y,
                sizeof(uint16_t) * tx_height * 2);
        if (txb_origin_y != 0 && txb_origin_x != 0)
            top_neigh_array[0] = left_neigh_array[0] =
                ((uint16_t *)(md_context_ptr->tx_search_luma_recon_neighbor_array16bit
                                  ->top_left_array) +
                 MAX_PICTURE_HEIGHT_SIZE + txb_origin_x - txb_origin_y)[0];

        mode = candidate_buffer_ptr->candidate_ptr->pred_mode;
        eb_av1_predict_intra_block_16bit(
            EB_10BIT,
            &md_context_ptr->sb_ptr->tile_info,
            !ED_STAGE,
            md_context_ptr->blk_geom,
            pcs_ptr->parent_pcs_ptr->av1_cm,
            md_context_ptr->blk_geom->bwidth,
            md_context_ptr->blk_geom->bheight,
            tx_size,
            mode,
            candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
            candidate_buffer_ptr->candidate_ptr->palette_info.pmi.palette_size[0] > 0,
            &candidate_buffer_ptr->candidate_ptr->palette_info, //ATB MD
            candidate_buffer_ptr->candidate_ptr->filter_intra_mode,
            top_neigh_array + 1,
            left_neigh_array + 1,
            candidate_buffer_ptr->prediction_ptr,
            (md_context_ptr->blk_geom
                 ->tx_org_x[is_inter][md_context_ptr->tx_depth][md_context_ptr->txb_itr] -
             md_context_ptr->blk_geom->origin_x) >>
                2,
            (md_context_ptr->blk_geom
                 ->tx_org_y[is_inter][md_context_ptr->tx_depth][md_context_ptr->txb_itr] -
             md_context_ptr->blk_geom->origin_y) >>
                2,
            PLANE_TYPE_Y,
            md_context_ptr->blk_geom->bsize,
            md_context_ptr->blk_origin_x,
            md_context_ptr->blk_origin_y,
            md_context_ptr->blk_origin_x,
            md_context_ptr->blk_origin_y,
            md_context_ptr->blk_geom
                ->tx_org_x[is_inter][md_context_ptr->tx_depth]
                          [md_context_ptr->txb_itr], //uint32_t cuOrgX used only for prediction Ptr
            md_context_ptr->blk_geom
                ->tx_org_y[is_inter][md_context_ptr->tx_depth]
                          [md_context_ptr->txb_itr], //uint32_t cuOrgY used only for prediction Ptr
            pcs_ptr->mi_grid_base,
            &((SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr)->seq_header);
    }

    return return_error;
}

static void tx_search_update_recon_sample_neighbor_array(
    NeighborArrayUnit *lumaReconSampleNeighborArray, EbPictureBufferDesc *recon_buffer,
    uint32_t txb_origin_x, uint32_t txb_origin_y, uint32_t input_origin_x, uint32_t input_origin_y,
    uint32_t width, uint32_t height, EbBool hbd) {
    if (hbd) {
        neighbor_array_unit16bit_sample_write(lumaReconSampleNeighborArray,
                                              (uint16_t *)recon_buffer->buffer_y,
                                              recon_buffer->stride_y,
                                              recon_buffer->origin_x + txb_origin_x,
                                              recon_buffer->origin_y + txb_origin_y,
                                              input_origin_x,
                                              input_origin_y,
                                              width,
                                              height,
                                              NEIGHBOR_ARRAY_UNIT_FULL_MASK);
    } else {
        neighbor_array_unit_sample_write(lumaReconSampleNeighborArray,
                                         recon_buffer->buffer_y,
                                         recon_buffer->stride_y,
                                         recon_buffer->origin_x + txb_origin_x,
                                         recon_buffer->origin_y + txb_origin_y,
                                         input_origin_x,
                                         input_origin_y,
                                         width,
                                         height,
                                         NEIGHBOR_ARRAY_UNIT_FULL_MASK);
    }

    return;
}
uint8_t get_end_tx_depth(BlockSize bsize) {
    uint8_t tx_depth = 0;
    if (bsize == BLOCK_64X64 || bsize == BLOCK_32X32 || bsize == BLOCK_16X16 ||
        bsize == BLOCK_64X32 || bsize == BLOCK_32X64 || bsize == BLOCK_16X32 ||
        bsize == BLOCK_32X16 || bsize == BLOCK_16X8 || bsize == BLOCK_8X16)
        tx_depth = 2;
    else if (bsize == BLOCK_8X8 || bsize == BLOCK_64X16 || bsize == BLOCK_16X64 ||
             bsize == BLOCK_32X8 || bsize == BLOCK_8X32 || bsize == BLOCK_16X4 ||
             bsize == BLOCK_4X16)
        tx_depth = 1;
    return tx_depth;
}

extern uint8_t allowed_tx_set_a[TX_SIZES_ALL][TX_TYPES];

void tx_initialize_neighbor_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                   EbBool is_inter) {
    uint16_t tile_idx = context_ptr->tile_index;
    // Set recon neighbor array to be used @ intra compensation
    if (!is_inter) {
        if (context_ptr->hbd_mode_decision)
            context_ptr->tx_search_luma_recon_neighbor_array16bit = context_ptr->tx_depth == 2
                ? pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX]
                                                                       [tile_idx]
                : context_ptr->tx_depth == 1
                    ? pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX]
                                                                           [tile_idx]
                    : pcs_ptr->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        else
            context_ptr->tx_search_luma_recon_neighbor_array = context_ptr->tx_depth == 2
                ? pcs_ptr
                      ->md_tx_depth_2_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx]
                : context_ptr->tx_depth == 1
                    ? pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                      [tile_idx]
                    : pcs_ptr->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    }
    // Set luma dc sign level coeff
    context_ptr->full_loop_luma_dc_sign_level_coeff_neighbor_array = (context_ptr->tx_depth)
        ? pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                        [tile_idx]
        : pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
}

void tx_update_neighbor_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                               ModeDecisionCandidateBuffer *candidate_buffer, EbBool is_inter) {
    uint16_t tile_idx = context_ptr->tile_index;
    if (context_ptr->tx_depth) {
        if (!is_inter)
            tx_search_update_recon_sample_neighbor_array(
                context_ptr->hbd_mode_decision
                    ? context_ptr->tx_search_luma_recon_neighbor_array16bit
                    : context_ptr->tx_search_luma_recon_neighbor_array,
                candidate_buffer->recon_ptr,
                context_ptr->blk_geom
                    ->tx_org_x[is_inter][context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom
                    ->tx_org_y[is_inter][context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->sb_origin_x +
                    context_ptr->blk_geom
                        ->tx_org_x[is_inter][context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->sb_origin_y +
                    context_ptr->blk_geom
                        ->tx_org_y[is_inter][context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->hbd_mode_decision);

        int8_t dc_sign_level_coeff =
            candidate_buffer->candidate_ptr->quantized_dc[0][context_ptr->txb_itr];
        neighbor_array_unit_mode_write(
            pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                          [tile_idx],
            (uint8_t *)&dc_sign_level_coeff,
            context_ptr->sb_origin_x +
                context_ptr->blk_geom
                    ->tx_org_x[is_inter][context_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->sb_origin_y +
                context_ptr->blk_geom
                    ->tx_org_y[is_inter][context_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr],
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }
}
void tx_reset_neighbor_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                              EbBool is_inter, uint8_t tx_depth) {
    int      sb_size  = pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.super_block_size;
    uint16_t tile_idx = context_ptr->tile_index;
    if (tx_depth) {
        if (!is_inter) {
            if (context_ptr->hbd_mode_decision) {
                if (tx_depth == 2) {
                    copy_neigh_arr(
                        pcs_ptr
                            ->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                        pcs_ptr
                            ->md_tx_depth_2_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX]
                                                                          [tile_idx],
                        context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
                        context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
                        context_ptr->blk_geom->bwidth,
                        context_ptr->blk_geom->bheight,
                        NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK);

                    copy_neigh_arr(
                        pcs_ptr
                            ->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                        pcs_ptr
                            ->md_tx_depth_2_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX]
                                                                          [tile_idx],
                        context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
                        context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
                        context_ptr->blk_geom->bwidth * 2,
                        MIN(context_ptr->blk_geom->bheight * 2,
                            sb_size - context_ptr->blk_geom->origin_y),
                        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
                } else {
                    copy_neigh_arr(
                        pcs_ptr
                            ->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                        pcs_ptr
                            ->md_tx_depth_1_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX]
                                                                          [tile_idx],
                        context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
                        context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
                        context_ptr->blk_geom->bwidth,
                        context_ptr->blk_geom->bheight,
                        NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK);

                    copy_neigh_arr(
                        pcs_ptr
                            ->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                        pcs_ptr
                            ->md_tx_depth_1_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX]
                                                                          [tile_idx],
                        context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
                        context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
                        context_ptr->blk_geom->bwidth * 2,
                        MIN(context_ptr->blk_geom->bheight * 2,
                            sb_size - context_ptr->blk_geom->origin_y),
                        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
                }
            } else {
                if (tx_depth == 2) {
                    copy_neigh_arr(
                        pcs_ptr->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                        pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                        [tile_idx],
                        context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
                        context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
                        context_ptr->blk_geom->bwidth,
                        context_ptr->blk_geom->bheight,
                        NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK);
                    copy_neigh_arr(
                        pcs_ptr->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                        pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                        [tile_idx],
                        context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
                        context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
                        context_ptr->blk_geom->bwidth * 2,
                        MIN(context_ptr->blk_geom->bheight * 2,
                            sb_size - context_ptr->blk_geom->origin_y),
                        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
                } else {
                    copy_neigh_arr(
                        pcs_ptr->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                        pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                        [tile_idx],
                        context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
                        context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
                        context_ptr->blk_geom->bwidth,
                        context_ptr->blk_geom->bheight,
                        NEIGHBOR_ARRAY_UNIT_TOPLEFT_MASK);
                    copy_neigh_arr(
                        pcs_ptr->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
                        pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                        [tile_idx],
                        context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
                        context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
                        context_ptr->blk_geom->bwidth * 2,
                        MIN(context_ptr->blk_geom->bheight * 2,
                            sb_size - context_ptr->blk_geom->origin_y),
                        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
                }
            }
        }
        copy_neigh_arr(
            pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
            pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                          [tile_idx],
            context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x,
            context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }
}

void tx_type_search(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                    ModeDecisionCandidateBuffer *candidate_buffer, uint32_t qp) {
    EbPictureBufferDesc *input_picture_ptr = context_ptr->hbd_mode_decision
        ? pcs_ptr->input_frame16bit
        : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    int32_t seg_qp = pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params.segmentation_enabled
        ? pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params
              .feature_data[context_ptr->blk_ptr->segment_id][SEG_LVL_ALT_Q]
        : 0;

    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];
    TxType   txk_start           = DCT_DCT;
    TxType   txk_end             = TX_TYPES;
    uint64_t best_cost_tx_search = (uint64_t)~0;
    int32_t  tx_type;
    TxSize   tx_size  = context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr];
    int32_t  is_inter = (candidate_buffer->candidate_ptr->type == INTER_MODE ||
                        candidate_buffer->candidate_ptr->use_intrabc)
        ? EB_TRUE
        : EB_FALSE;
    const TxSetType tx_set_type = get_ext_tx_set_type(
        tx_size, is_inter, pcs_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set);
    uint8_t txb_origin_x = (uint8_t)context_ptr->blk_geom
                               ->tx_org_x[is_inter][context_ptr->tx_depth][context_ptr->txb_itr];
    uint8_t txb_origin_y = (uint8_t)context_ptr->blk_geom
                               ->tx_org_y[is_inter][context_ptr->tx_depth][context_ptr->txb_itr];
    uint32_t txb_origin_index = txb_origin_x +
        (txb_origin_y * candidate_buffer->residual_ptr->stride_y);
    uint32_t input_txb_origin_index = (context_ptr->sb_origin_x + txb_origin_x +
                                       input_picture_ptr->origin_x) +
        ((context_ptr->sb_origin_y + txb_origin_y + input_picture_ptr->origin_y) *
         input_picture_ptr->stride_y);

    context_ptr->luma_txb_skip_context = 0;
    context_ptr->luma_dc_sign_context  = 0;
    get_txb_ctx(pcs_ptr,
                COMPONENT_LUMA,
                context_ptr->full_loop_luma_dc_sign_level_coeff_neighbor_array,
                context_ptr->sb_origin_x + txb_origin_x,
                context_ptr->sb_origin_y + txb_origin_y,
                context_ptr->blk_geom->bsize,
                context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr],
                &context_ptr->luma_txb_skip_context,
                &context_ptr->luma_dc_sign_context);
    if (context_ptr->tx_search_reduced_set == 2)
        txk_end = 2;
    TxType best_tx_type = DCT_DCT;
    for (tx_type = txk_start; tx_type < txk_end; ++tx_type) {
        uint64_t txb_full_distortion[3][DIST_CALC_TOTAL];
        uint64_t y_txb_coeff_bits = 0;
        uint32_t y_count_non_zero_coeffs;

        if (context_ptr->tx_search_reduced_set == 2)
            tx_type = (tx_type == 1) ? IDTX : tx_type;

        context_ptr->three_quad_energy = 0;
        if (tx_type != DCT_DCT) {
            if (is_inter) {
                TxSize          max_tx_size       = context_ptr->blk_geom->txsize[0][0];
                const TxSetType tx_set_type_inter = get_ext_tx_set_type(
                    max_tx_size, is_inter, pcs_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set);
                int32_t eset = get_ext_tx_set(
                    max_tx_size, is_inter, pcs_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set);
                // eset == 0 should correspond to a set with only DCT_DCT and there
                // is no need to send the tx_type
                if (eset <= 0)
                    continue;
                else if (av1_ext_tx_used[tx_set_type_inter][tx_type] == 0)
                    continue;
                else if (context_ptr->blk_geom
                                 ->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] > 32 ||
                         context_ptr->blk_geom
                                 ->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] > 32)
                    continue;
            }
            int32_t eset = get_ext_tx_set(
                context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr],
                is_inter,
                context_ptr->tx_search_reduced_set);
            // eset == 0 should correspond to a set with only DCT_DCT and there
            // is no need to send the tx_type
            if (eset <= 0)
                continue;
            else if (av1_ext_tx_used[tx_set_type][tx_type] == 0)
                continue;
            else if (context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] >
                         32 ||
                     context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] >
                         32)
                continue;
        }
        if (context_ptr->tx_search_reduced_set)
            if (!allowed_tx_set_a[context_ptr->blk_geom->txsize[context_ptr->tx_depth]
                                                               [context_ptr->txb_itr]][tx_type])
                continue;

        // For Inter blocks, transform type of chroma follows luma transfrom type
        if (is_inter)
            candidate_buffer->candidate_ptr->transform_type_uv = (context_ptr->txb_itr == 0)
                ? candidate_buffer->candidate_ptr->transform_type[context_ptr->txb_itr]
                : candidate_buffer->candidate_ptr->transform_type_uv;

        // Y: T Q i_q
        av1_estimate_transform(
            &(((int16_t *)candidate_buffer->residual_ptr->buffer_y)[txb_origin_index]),
            candidate_buffer->residual_ptr->stride_y,
            &(((int32_t *)context_ptr->trans_quant_buffers_ptr->txb_trans_coeff2_nx2_n_ptr
                   ->buffer_y)[context_ptr->txb_1d_offset]),
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr],
            &context_ptr->three_quad_energy,
            context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT,
            tx_type,
            PLANE_TYPE_Y,
            DEFAULT_SHAPE);

        av1_quantize_inv_quantize(
            pcs_ptr,
            context_ptr,
            &(((int32_t *)context_ptr->trans_quant_buffers_ptr->txb_trans_coeff2_nx2_n_ptr
                   ->buffer_y)[context_ptr->txb_1d_offset]),
            NOT_USED_VALUE,
            &(((int32_t *)candidate_buffer->residual_quant_coeff_ptr
                   ->buffer_y)[context_ptr->txb_1d_offset]),
            &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_y)[context_ptr->txb_1d_offset]),
            qp,
            seg_qp,
            context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr],
            &candidate_buffer->candidate_ptr->eob[0][context_ptr->txb_itr],
            &y_count_non_zero_coeffs,
            COMPONENT_LUMA,
            context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT,
            tx_type,
            candidate_buffer,
            context_ptr->luma_txb_skip_context,
            context_ptr->luma_dc_sign_context,
            candidate_buffer->candidate_ptr->pred_mode,
            candidate_buffer->candidate_ptr->use_intrabc,
            full_lambda,
            EB_FALSE);

        candidate_buffer->candidate_ptr->quantized_dc[0][context_ptr->txb_itr] = ((
            (int32_t *)
                candidate_buffer->residual_quant_coeff_ptr->buffer_y)[context_ptr->txb_1d_offset]);
        uint32_t y_has_coeff = y_count_non_zero_coeffs > 0;

        // tx_type not equal to DCT_DCT and no coeff is not an acceptable option in AV1.
        if (y_has_coeff == 0 && tx_type != DCT_DCT)
            continue;

        if (context_ptr->md_staging_spatial_sse_full_loop) {
            if (y_has_coeff)
                inv_transform_recon_wrapper(
                    candidate_buffer->prediction_ptr->buffer_y,
                    txb_origin_index,
                    candidate_buffer->prediction_ptr->stride_y,
                    candidate_buffer->recon_ptr->buffer_y,
                    txb_origin_index,
                    candidate_buffer->recon_ptr->stride_y,
                    (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_y,
                    context_ptr->txb_1d_offset,
                    context_ptr->hbd_mode_decision,
                    context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr],
                    tx_type,
                    PLANE_TYPE_Y,
                    (uint16_t)candidate_buffer->candidate_ptr->eob[0][context_ptr->txb_itr]);
            else
                picture_copy(
                    candidate_buffer->prediction_ptr,
                    txb_origin_index,
                    0,
                    candidate_buffer->recon_ptr,
                    txb_origin_index,
                    0,
                    context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr],
                    0,
                    0,
                    PICTURE_BUFFER_DESC_Y_FLAG,
                    context_ptr->hbd_mode_decision);

            EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
                ? full_distortion_kernel16_bits
                : spatial_full_distortion_kernel;

            txb_full_distortion[0][DIST_CALC_PREDICTION] = spatial_full_dist_type_fun(
                input_picture_ptr->buffer_y,
                input_txb_origin_index,
                input_picture_ptr->stride_y,
                candidate_buffer->prediction_ptr->buffer_y,
                (int32_t)txb_origin_index,
                candidate_buffer->prediction_ptr->stride_y,
                context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr]);

            txb_full_distortion[0][DIST_CALC_RESIDUAL] = spatial_full_dist_type_fun(
                input_picture_ptr->buffer_y,
                input_txb_origin_index,
                input_picture_ptr->stride_y,
                candidate_buffer->recon_ptr->buffer_y,
                (int32_t)txb_origin_index,
                candidate_buffer->recon_ptr->stride_y,
                context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr]);

            txb_full_distortion[0][DIST_CALC_PREDICTION] <<= 4;
            txb_full_distortion[0][DIST_CALC_RESIDUAL] <<= 4;
        } else {
            // LUMA DISTORTION
            picture_full_distortion32_bits(
                context_ptr->trans_quant_buffers_ptr->txb_trans_coeff2_nx2_n_ptr,
                context_ptr->txb_1d_offset,
                0,
                candidate_buffer->recon_coeff_ptr,
                context_ptr->txb_1d_offset,
                0,
                context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr],
                NOT_USED_VALUE,
                NOT_USED_VALUE,
                txb_full_distortion[0],
                NOT_USED_VALUE,
                NOT_USED_VALUE,
                y_count_non_zero_coeffs,
                0,
                0,
                COMPONENT_LUMA);

            txb_full_distortion[0][DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
            txb_full_distortion[0][DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;
            //assert(context_ptr->three_quad_energy == 0 && context_ptr->cu_stats->size < 64);
            TxSize tx_size =
                context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr];
            int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(tx_size)) * 2;
            txb_full_distortion[0][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(
                txb_full_distortion[0][DIST_CALC_RESIDUAL], shift);
            txb_full_distortion[0][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(
                txb_full_distortion[0][DIST_CALC_PREDICTION], shift);
        }

        //LUMA-ONLY
        av1_txb_estimate_coeff_bits(
            context_ptr,
            0, //allow_update_cdf,
            NULL, //FRAME_CONTEXT *ec_ctx,
            pcs_ptr,
            candidate_buffer,
            context_ptr->txb_1d_offset,
            0,
            context_ptr->coeff_est_entropy_coder_ptr,
            candidate_buffer->residual_quant_coeff_ptr,
            y_count_non_zero_coeffs,
            0,
            0,
            &y_txb_coeff_bits,
            &y_txb_coeff_bits,
            &y_txb_coeff_bits,
            context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->tx_depth][context_ptr->txb_itr],
            tx_type,
            candidate_buffer->candidate_ptr->transform_type_uv,
            COMPONENT_LUMA);

        uint64_t cost = RDCOST(
            full_lambda, y_txb_coeff_bits, txb_full_distortion[0][DIST_CALC_RESIDUAL]);
        if (cost < best_cost_tx_search) {
            best_cost_tx_search = cost;
            best_tx_type        = tx_type;
        }
    }

    //  Best Tx Type Pass
    candidate_buffer->candidate_ptr->transform_type[context_ptr->txb_itr] = best_tx_type;

    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
        candidate_buffer->candidate_ptr->transform_type_uv = (context_ptr->txb_itr == 0)
            ? candidate_buffer->candidate_ptr->transform_type[context_ptr->txb_itr]
            : candidate_buffer->candidate_ptr->transform_type_uv;
}

static INLINE int block_signals_txsize(BlockSize bsize) { return bsize > BLOCK_4X4; }

static INLINE int get_vartx_max_txsize(/*const MbModeInfo *xd,*/ BlockSize bsize, int plane) {
    /* if (xd->lossless[xd->mi[0]->segment_id]) return TX_4X4;*/
    const TxSize max_txsize = max_txsize_rect_lookup[bsize];
    if (plane == 0)
        return max_txsize; // luma
    return av1_get_adjusted_tx_size(max_txsize); // chroma
}

static INLINE int max_block_wide(const MacroBlockD *xd, BlockSize bsize, int plane) {
    int max_blocks_wide = block_size_wide[bsize];

    if (xd->mb_to_right_edge < 0)
        max_blocks_wide += xd->mb_to_right_edge >> (3 + (plane == 0 ? 0 : 1));

    // Scale the width in the transform block unit.
    return max_blocks_wide >> tx_size_wide_log2[0];
}

static INLINE int max_block_high(const MacroBlockD *xd, BlockSize bsize, int plane) {
    int max_blocks_high = block_size_high[bsize];

    if (xd->mb_to_bottom_edge < 0)
        max_blocks_high += xd->mb_to_bottom_edge >> (3 + (plane == 0 ? 0 : 1));

    // Scale the height in the transform block unit.
    return max_blocks_high >> tx_size_high_log2[0];
}

static INLINE void txfm_partition_update(TXFM_CONTEXT *above_ctx, TXFM_CONTEXT *left_ctx,
                                         TxSize tx_size, TxSize txb_size) {
    BlockSize bsize = txsize_to_bsize[txb_size];
    assert(bsize < BlockSizeS_ALL);
    int     bh  = mi_size_high[bsize];
    int     bw  = mi_size_wide[bsize];
    uint8_t txw = tx_size_wide[tx_size];
    uint8_t txh = tx_size_high[tx_size];
    int     i;
    for (i = 0; i < bh; ++i) left_ctx[i] = txh;
    for (i = 0; i < bw; ++i) above_ctx[i] = txw;
}

static INLINE TxSize get_sqr_tx_size(int tx_dim) {
    switch (tx_dim) {
    case 128:
    case 64: return TX_64X64; break;
    case 32: return TX_32X32; break;
    case 16: return TX_16X16; break;
    case 8: return TX_8X8; break;
    default: return TX_4X4;
    }
}
static INLINE int txfm_partition_context(TXFM_CONTEXT *above_ctx, TXFM_CONTEXT *left_ctx,
                                         BlockSize bsize, TxSize tx_size) {
    const uint8_t txw      = tx_size_wide[tx_size];
    const uint8_t txh      = tx_size_high[tx_size];
    const int     above    = *above_ctx < txw;
    const int     left     = *left_ctx < txh;
    int           category = TXFM_PARTITION_CONTEXTS;

    // dummy return, not used by others.
    if (tx_size <= TX_4X4)
        return 0;

    TxSize max_tx_size = get_sqr_tx_size(AOMMAX(block_size_wide[bsize], block_size_high[bsize]));

    if (max_tx_size >= TX_8X8) {
        category = (txsize_sqr_up_map[tx_size] != max_tx_size && max_tx_size > TX_8X8) +
            (TX_SIZES - 1 - max_tx_size) * 2;
    }
    assert(category != TXFM_PARTITION_CONTEXTS);
    return category * 3 + above + left;
}

static uint64_t cost_tx_size_vartx(MacroBlockD *xd, const MbModeInfo *mbmi, TxSize tx_size,
                                   int depth, int blk_row, int blk_col,
                                   MdRateEstimationContext *md_rate_estimation_ptr) {
    uint64_t  bits            = 0;
    const int max_blocks_high = max_block_high(xd, mbmi->block_mi.sb_type, 0);
    const int max_blocks_wide = max_block_wide(xd, mbmi->block_mi.sb_type, 0);

    if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
        return bits;

    if (depth == MAX_VARTX_DEPTH) {
        txfm_partition_update(
            xd->above_txfm_context + blk_col, xd->left_txfm_context + blk_row, tx_size, tx_size);

        return bits;
    }

    const int ctx = txfm_partition_context(xd->above_txfm_context + blk_col,
                                           xd->left_txfm_context + blk_row,
                                           mbmi->block_mi.sb_type,
                                           tx_size);

    const int write_txfm_partition = (tx_size ==
                                      tx_depth_to_tx_size[mbmi->tx_depth][mbmi->block_mi.sb_type]);

    if (write_txfm_partition) {
        bits += md_rate_estimation_ptr->txfm_partition_fac_bits[ctx][0];

        txfm_partition_update(
            xd->above_txfm_context + blk_col, xd->left_txfm_context + blk_row, tx_size, tx_size);

    } else {
        const TxSize sub_txs = sub_tx_size_map[tx_size];
        const int    bsw     = tx_size_wide_unit[sub_txs];
        const int    bsh     = tx_size_high_unit[sub_txs];

        bits += md_rate_estimation_ptr->txfm_partition_fac_bits[ctx][1];
        if (sub_txs == TX_4X4) {
            txfm_partition_update(xd->above_txfm_context + blk_col,
                                  xd->left_txfm_context + blk_row,
                                  sub_txs,
                                  tx_size);

            return bits;
        }

        assert(bsw > 0 && bsh > 0);
        for (int row = 0; row < tx_size_high_unit[tx_size]; row += bsh)
            for (int col = 0; col < tx_size_wide_unit[tx_size]; col += bsw) {
                int offsetr = blk_row + row;
                int offsetc = blk_col + col;
                bits += cost_tx_size_vartx(
                    xd, mbmi, sub_txs, depth + 1, offsetr, offsetc, md_rate_estimation_ptr);
            }
    }
    return bits;
}

static INLINE void set_txfm_ctx(TXFM_CONTEXT *txfm_ctx, uint8_t txs, int len) {
    int i;
    for (i = 0; i < len; ++i) txfm_ctx[i] = txs;
}

static INLINE void set_txfm_ctxs(TxSize tx_size, int n8_w, int n8_h, int skip,
                                 const MacroBlockD *xd) {
    uint8_t bw = tx_size_wide[tx_size];
    uint8_t bh = tx_size_high[tx_size];

    if (skip) {
        bw = n8_w * MI_SIZE;
        bh = n8_h * MI_SIZE;
    }

    set_txfm_ctx(xd->above_txfm_context, bw, n8_w);
    set_txfm_ctx(xd->left_txfm_context, bh, n8_h);
}

static INLINE int tx_size_to_depth(TxSize tx_size, BlockSize bsize) {
    TxSize ctx_size = max_txsize_rect_lookup[bsize];
    int    depth    = 0;
    while (tx_size != ctx_size) {
        depth++;
        ctx_size = sub_tx_size_map[ctx_size];
        assert(depth <= MAX_TX_DEPTH);
    }
    return depth;
}

#define BLOCK_SIZES_ALL 22

int is_inter_block(const BlockModeInfo *bloc_mi);

// Returns a context number for the given MB prediction signal
// The mode info data structure has a one element border above and to the
// left of the entries corresponding to real blocks.
// The prediction flags in these dummy entries are initialized to 0.
static INLINE int get_tx_size_context(const MacroBlockD *xd) {
    const ModeInfo *        mi          = xd->mi[0];
    const MbModeInfo *      mbmi        = &mi->mbmi;
    const MbModeInfo *const above_mbmi  = xd->above_mbmi;
    const MbModeInfo *const left_mbmi   = xd->left_mbmi;
    const TxSize            max_tx_size = max_txsize_rect_lookup[mbmi->block_mi.sb_type];
    const int               max_tx_wide = tx_size_wide[max_tx_size];
    const int               max_tx_high = tx_size_high[max_tx_size];
    const int               has_above   = xd->up_available;
    const int               has_left    = xd->left_available;

    int above = xd->above_txfm_context[0] >= max_tx_wide;
    int left  = xd->left_txfm_context[0] >= max_tx_high;

    if (has_above)
        if (is_inter_block(&above_mbmi->block_mi))
            above = block_size_wide[above_mbmi->block_mi.sb_type] >= max_tx_wide;

    if (has_left)
        if (is_inter_block(&left_mbmi->block_mi))
            left = block_size_high[left_mbmi->block_mi.sb_type] >= max_tx_high;

    if (has_above && has_left)
        return (above + left);
    else if (has_above)
        return above;
    else if (has_left)
        return left;
    else
        return 0;
}

static uint64_t cost_selected_tx_size(const MacroBlockD *      xd,
                                      MdRateEstimationContext *md_rate_estimation_ptr) {
    const ModeInfo *const   mi    = xd->mi[0];
    const MbModeInfo *const mbmi  = &mi->mbmi;
    const BlockSize         bsize = mbmi->block_mi.sb_type;
    uint64_t                bits  = 0;

    if (block_signals_txsize(bsize)) {
        const TxSize tx_size     = mbmi->tx_size;
        const int    tx_size_ctx = get_tx_size_context(xd);
        assert(bsize < BlockSizeS_ALL);
        const int     depth       = tx_size_to_depth(tx_size, bsize);
        const int32_t tx_size_cat = bsize_to_tx_size_cat(bsize);
        bits += md_rate_estimation_ptr->tx_size_fac_bits[tx_size_cat][tx_size_ctx][depth];
    }

    return bits;
}

static uint64_t tx_size_bits(MdRateEstimationContext *md_rate_estimation_ptr, MacroBlockD *xd,
                             const MbModeInfo *mbmi, TxMode tx_mode, BlockSize bsize,
                             uint8_t skip) {
    uint64_t bits = 0;

    int is_inter_tx = is_inter_block(&mbmi->block_mi) || is_intrabc_block(&mbmi->block_mi);
    if (tx_mode == TX_MODE_SELECT && block_signals_txsize(bsize) &&
        !(is_inter_tx && skip) /*&& !xd->lossless[segment_id]*/) {
        if (is_inter_tx) { // This implies skip flag is 0.
            const TxSize max_tx_size = get_vartx_max_txsize(/*xd,*/ bsize, 0);
            const int    txbh        = tx_size_high_unit[max_tx_size];
            const int    txbw        = tx_size_wide_unit[max_tx_size];
            const int    width       = block_size_wide[bsize] >> tx_size_wide_log2[0];
            const int    height      = block_size_high[bsize] >> tx_size_high_log2[0];
            int          idx, idy;
            for (idy = 0; idy < height; idy += txbh)
                for (idx = 0; idx < width; idx += txbw)
                    bits += cost_tx_size_vartx(
                        xd, mbmi, max_tx_size, 0, idy, idx, md_rate_estimation_ptr);
        } else {
            bits += cost_selected_tx_size(xd, md_rate_estimation_ptr);
            set_txfm_ctxs(mbmi->tx_size, xd->n8_w, xd->n8_h, 0, xd);
        }
    } else {
        set_txfm_ctxs(
            mbmi->tx_size, xd->n8_w, xd->n8_h, skip && is_inter_block(&mbmi->block_mi), xd);
    }
    return bits;
}

void set_mi_row_col(PictureControlSet *pcs_ptr, MacroBlockD *xd, TileInfo *tile, int mi_row, int bh,
                    int mi_col, int bw, uint32_t mi_stride, int mi_rows, int mi_cols);

uint64_t estimate_tx_size_bits(PictureControlSet *pcsPtr, ModeDecisionContext *context_ptr,
                               ModeDecisionCandidate *candidate_ptr, EbBool skip_flag,
                               uint32_t blk_origin_x, uint32_t blk_origin_y, BlkStruct *blk_ptr,
                               const BlockGeom *blk_geom, NeighborArrayUnit *txfm_context_array,
                               uint8_t tx_depth, MdRateEstimationContext *md_rate_estimation_ptr) {
    uint32_t txfm_context_left_index  = get_neighbor_array_unit_left_index(txfm_context_array,
                                                                          blk_origin_y);
    uint32_t txfm_context_above_index = get_neighbor_array_unit_top_index(txfm_context_array,
                                                                          blk_origin_x);

    TxMode        tx_mode   = pcsPtr->parent_pcs_ptr->frm_hdr.tx_mode;
    Av1Common *   cm        = pcsPtr->parent_pcs_ptr->av1_cm;
    MacroBlockD * xd        = blk_ptr->av1xd;
    TileInfo *    tile      = &xd->tile;
    int32_t       mi_row    = blk_origin_y >> MI_SIZE_LOG2;
    int32_t       mi_col    = blk_origin_x >> MI_SIZE_LOG2;
    BlockSize     bsize     = blk_geom->bsize;
    const int32_t bw        = mi_size_wide[bsize];
    const int32_t bh        = mi_size_high[bsize];
    uint32_t      mi_stride = pcsPtr->mi_stride;

    set_mi_row_col(pcsPtr, xd, tile, mi_row, bh, mi_col, bw, mi_stride, cm->mi_rows, cm->mi_cols);

    MbModeInfo *mbmi = &xd->mi[0]->mbmi;

    memcpy(context_ptr->above_txfm_context,
           &(txfm_context_array->top_array[txfm_context_above_index]),
           (blk_geom->bwidth >> MI_SIZE_LOG2) * sizeof(TXFM_CONTEXT));
    memcpy(context_ptr->left_txfm_context,
           &(txfm_context_array->left_array[txfm_context_left_index]),
           (blk_geom->bheight >> MI_SIZE_LOG2) * sizeof(TXFM_CONTEXT));

    xd->above_txfm_context = context_ptr->above_txfm_context;
    xd->left_txfm_context  = context_ptr->left_txfm_context;

    mbmi->tx_size               = blk_geom->txsize[tx_depth][0];
    mbmi->block_mi.sb_type      = blk_geom->bsize;
    mbmi->block_mi.use_intrabc  = candidate_ptr->use_intrabc;
    mbmi->block_mi.ref_frame[0] = candidate_ptr->ref_frame_type;
    mbmi->tx_depth              = tx_depth;

    uint64_t bits = tx_size_bits(md_rate_estimation_ptr, xd, mbmi, tx_mode, bsize, skip_flag);

    return bits;
}

uint64_t get_tx_size_bits(ModeDecisionCandidateBuffer *candidateBuffer,
                          ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                          uint8_t tx_depth, EbBool block_has_coeff) {
    uint64_t tx_size_bits = 0;

    tx_size_bits = estimate_tx_size_bits(pcs_ptr,
                                         context_ptr,
                                         candidateBuffer->candidate_ptr,
                                         block_has_coeff ? 0 : 1,
                                         context_ptr->blk_origin_x,
                                         context_ptr->blk_origin_y,
                                         context_ptr->blk_ptr,
                                         context_ptr->blk_geom,
                                         context_ptr->txfm_context_array,
                                         tx_depth,
                                         context_ptr->md_rate_estimation_ptr);

    return tx_size_bits;
}

void init_tx_candidate_buffer(ModeDecisionCandidateBuffer *candidate_buffer,
                              ModeDecisionContext *context_ptr, uint8_t end_tx_depth) {
    uint32_t block_index = context_ptr->blk_geom->origin_x +
        (context_ptr->blk_geom->origin_y * MAX_SB_SIZE);
    if (end_tx_depth) {
        memcpy(context_ptr->candidate_buffer_tx_depth_1->candidate_ptr,
               candidate_buffer->candidate_ptr,
               sizeof(ModeDecisionCandidate));
        memcpy(context_ptr->candidate_buffer_tx_depth_2->candidate_ptr,
               candidate_buffer->candidate_ptr,
               sizeof(ModeDecisionCandidate));
    }
    if (end_tx_depth) {
        if (context_ptr->hbd_mode_decision) {
            // Copy pred to tx_depth_1 candidate_buffer
            {
                uint16_t *src = &(
                    ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y)[block_index]);
                uint16_t *dst = &(((uint16_t *)context_ptr->candidate_buffer_tx_depth_1
                                       ->prediction_ptr->buffer_y)[block_index]);
                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
                    src += candidate_buffer->prediction_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_1->prediction_ptr->stride_y;
                }
            }
            // Copy residual to tx_depth_1 candidate_buffer
            {
                int16_t *src = &(
                    ((int16_t *)candidate_buffer->residual_ptr->buffer_y)[block_index]);
                int16_t *dst = &(((int16_t *)context_ptr->candidate_buffer_tx_depth_1->residual_ptr
                                      ->buffer_y)[block_index]);

                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    memcpy(dst, src, context_ptr->blk_geom->bwidth << 1);
                    src += candidate_buffer->residual_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_1->residual_ptr->stride_y;
                }
            }
        } else {
            // Copy pred to tx_depth_1 candidate_buffer
            {
                EbByte src = &(candidate_buffer->prediction_ptr->buffer_y[block_index]);
                EbByte dst = &(context_ptr->candidate_buffer_tx_depth_1->prediction_ptr
                                   ->buffer_y[block_index]);
                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    memcpy(dst, src, context_ptr->blk_geom->bwidth);
                    src += candidate_buffer->prediction_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_1->prediction_ptr->stride_y;
                }
            }
            // Copy residual to tx_depth_1 candidate_buffer
            {
                int16_t *src = &(
                    ((int16_t *)candidate_buffer->residual_ptr->buffer_y)[block_index]);
                int16_t *dst = &(((int16_t *)context_ptr->candidate_buffer_tx_depth_1->residual_ptr
                                      ->buffer_y)[block_index]);

                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    memcpy(dst, src, context_ptr->blk_geom->bwidth << 1);
                    src += candidate_buffer->residual_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_1->residual_ptr->stride_y;
                }
            }
        }
    }
    if (end_tx_depth == 2) {
        if (context_ptr->hbd_mode_decision) {
            // Copy pred to tx_depth_1 candidate_buffer
            {
                uint16_t *src = &(
                    ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y)[block_index]);
                uint16_t *dst = &(((uint16_t *)context_ptr->candidate_buffer_tx_depth_2
                                       ->prediction_ptr->buffer_y)[block_index]);

                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
                    src += candidate_buffer->prediction_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_2->prediction_ptr->stride_y;
                }
            }
            // Copy residual to tx_depth_1 candidate_buffer
            {
                int16_t *src = &(
                    ((int16_t *)candidate_buffer->residual_ptr->buffer_y)[block_index]);
                int16_t *dst = &(((int16_t *)context_ptr->candidate_buffer_tx_depth_2->residual_ptr
                                      ->buffer_y)[block_index]);

                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    memcpy(dst, src, context_ptr->blk_geom->bwidth << 1);
                    src += candidate_buffer->residual_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_2->residual_ptr->stride_y;
                }
            }
        } else {
            // Copy pred to tx_depth_2 candidate_buffer
            {
                EbByte src = &(candidate_buffer->prediction_ptr->buffer_y[block_index]);
                EbByte dst = &(context_ptr->candidate_buffer_tx_depth_2->prediction_ptr
                                   ->buffer_y[block_index]);
                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    memcpy(dst, src, context_ptr->blk_geom->bwidth);
                    src += candidate_buffer->prediction_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_2->prediction_ptr->stride_y;
                }
            }
            // Copy residual to tx_depth_2 candidate_buffer
            {
                int16_t *src = &(
                    ((int16_t *)candidate_buffer->residual_ptr->buffer_y)[block_index]);
                int16_t *dst = &(((int16_t *)context_ptr->candidate_buffer_tx_depth_2->residual_ptr
                                      ->buffer_y)[block_index]);

                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    memcpy(dst, src, context_ptr->blk_geom->bwidth << 1);
                    src += candidate_buffer->residual_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_2->residual_ptr->stride_y;
                }
            }
        }
    }
}

void update_tx_candidate_buffer(ModeDecisionCandidateBuffer *candidate_buffer,
                                ModeDecisionContext *context_ptr, uint8_t best_tx_depth) {
    uint32_t block_index = context_ptr->blk_geom->origin_x +
        (context_ptr->blk_geom->origin_y * MAX_SB_SIZE);
    if (best_tx_depth == 1) {
        // Copy depth 1 mode/type/eob ..
        memcpy(candidate_buffer->candidate_ptr,
               context_ptr->candidate_buffer_tx_depth_1->candidate_ptr,
               sizeof(ModeDecisionCandidate));
        // Copy depth 1 pred
        if (context_ptr->hbd_mode_decision) {
            uint16_t *src = &(((uint16_t *)context_ptr->candidate_buffer_tx_depth_1->prediction_ptr
                                   ->buffer_y)[block_index]);
            uint16_t *dst = &(
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y)[block_index]);
            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
                src += context_ptr->candidate_buffer_tx_depth_1->prediction_ptr->stride_y;
                dst += candidate_buffer->prediction_ptr->stride_y;
            }
        } else {
            EbByte src = &(
                context_ptr->candidate_buffer_tx_depth_1->prediction_ptr->buffer_y[block_index]);
            EbByte dst = &(candidate_buffer->prediction_ptr->buffer_y[block_index]);
            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                memcpy(dst, src, context_ptr->blk_geom->bwidth);
                src += context_ptr->candidate_buffer_tx_depth_1->prediction_ptr->stride_y;
                dst += candidate_buffer->prediction_ptr->stride_y;
            }
        }
        // Copy depth 1 recon coeff
        memcpy(candidate_buffer->recon_coeff_ptr->buffer_y,
               context_ptr->candidate_buffer_tx_depth_1->recon_coeff_ptr->buffer_y,
               (context_ptr->blk_geom->bwidth * context_ptr->blk_geom->bheight << 2));
    }
    if (best_tx_depth == 2) {
        // Copy depth 2 mode/type/eob ..
        memcpy(candidate_buffer->candidate_ptr,
               context_ptr->candidate_buffer_tx_depth_2->candidate_ptr,
               sizeof(ModeDecisionCandidate));
        // Copy depth 2 pred
        if (context_ptr->hbd_mode_decision) {
            uint16_t *src = &(((uint16_t *)context_ptr->candidate_buffer_tx_depth_2->prediction_ptr
                                   ->buffer_y)[block_index]);
            uint16_t *dst = &(
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y)[block_index]);
            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
                src += context_ptr->candidate_buffer_tx_depth_2->prediction_ptr->stride_y;
                dst += candidate_buffer->prediction_ptr->stride_y;
            }
        } else {
            EbByte src = &(
                context_ptr->candidate_buffer_tx_depth_2->prediction_ptr->buffer_y[block_index]);
            EbByte dst = &(candidate_buffer->prediction_ptr->buffer_y[block_index]);
            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                memcpy(dst, src, context_ptr->blk_geom->bwidth);
                src += context_ptr->candidate_buffer_tx_depth_2->prediction_ptr->stride_y;
                dst += candidate_buffer->prediction_ptr->stride_y;
            }
        }
        // Copy depth 2 recon coeff
        memcpy(candidate_buffer->recon_coeff_ptr->buffer_y,
               context_ptr->candidate_buffer_tx_depth_2->recon_coeff_ptr->buffer_y,
               (context_ptr->blk_geom->bwidth * context_ptr->blk_geom->bheight << 2));
    }
}
void perform_tx_partitioning(ModeDecisionCandidateBuffer *candidate_buffer,
                             ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                             uint64_t ref_fast_cost, uint8_t start_tx_depth, uint8_t end_tx_depth,
                             uint32_t qp, uint32_t *y_count_non_zero_coeffs, uint64_t *y_coeff_bits,
                             uint64_t *y_full_distortion) {
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];
    EbPictureBufferDesc *input_picture_ptr = context_ptr->hbd_mode_decision
        ? pcs_ptr->input_frame16bit
        : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    int32_t is_inter = (candidate_buffer->candidate_ptr->type == INTER_MODE ||
                        candidate_buffer->candidate_ptr->use_intrabc)
        ? EB_TRUE
        : EB_FALSE;

    uint8_t  best_tx_depth     = 0;
    uint64_t best_cost_search  = (uint64_t)~0;
    uint8_t  is_best_has_coeff = 1;
    init_tx_candidate_buffer(candidate_buffer, context_ptr, end_tx_depth);

    uint8_t tx_search_skip_flag;
    if (context_ptr->md_staging_tx_search == 0)
        tx_search_skip_flag = EB_TRUE;
    else if (context_ptr->md_staging_tx_search == 1) {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected && context_ptr->blk_geom->shape == PART_N)
            tx_search_skip_flag = context_ptr->tx_search_level == TX_SEARCH_FULL_LOOP ? EB_FALSE
                                                                                      : EB_TRUE;
        else
            tx_search_skip_flag = context_ptr->tx_search_level == TX_SEARCH_FULL_LOOP
                ? get_skip_tx_search_flag(context_ptr->blk_geom->sq_size,
                                          ref_fast_cost,
                                          *candidate_buffer->fast_cost_ptr,
                                          context_ptr->tx_weight)
                : EB_TRUE;
    } else
        tx_search_skip_flag = context_ptr->tx_search_level == TX_SEARCH_FULL_LOOP ? EB_FALSE
                                                                                  : EB_TRUE;
    // Transform Depth Loop
    for (context_ptr->tx_depth = start_tx_depth; context_ptr->tx_depth <= end_tx_depth;
         context_ptr->tx_depth++) {
        if (pcs_ptr->parent_pcs_ptr->tx_size_early_exit) {
            if (!is_best_has_coeff)
                continue;
        }
        tx_reset_neighbor_arrays(pcs_ptr, context_ptr, is_inter, context_ptr->tx_depth);
        ModeDecisionCandidateBuffer *tx_candidate_buffer = (context_ptr->tx_depth == 0)
            ? candidate_buffer
            : (context_ptr->tx_depth == 1) ? context_ptr->candidate_buffer_tx_depth_1
                                           : context_ptr->candidate_buffer_tx_depth_2;
        tx_candidate_buffer->candidate_ptr->tx_depth = context_ptr->tx_depth;

        tx_initialize_neighbor_arrays(pcs_ptr, context_ptr, is_inter);

        // Initialize TU Split
        uint32_t tx_y_count_non_zero_coeffs[MAX_NUM_OF_TU_PER_CU];
        uint64_t tx_y_coeff_bits                       = 0;
        uint64_t tx_y_full_distortion[DIST_CALC_TOTAL] = {0};

        // Current tx_cost while performing tx loop
        uint64_t current_tx_cost = 0;

        context_ptr->txb_1d_offset                      = 0;
        context_ptr->three_quad_energy                  = 0;
        tx_candidate_buffer->candidate_ptr->y_has_coeff = 0;

        uint16_t txb_count = context_ptr->blk_geom->txb_count[context_ptr->tx_depth];

        uint32_t block_has_coeff = EB_FALSE;
        for (context_ptr->txb_itr = 0; context_ptr->txb_itr < txb_count; context_ptr->txb_itr++) {
            uint16_t tx_org_x = context_ptr->blk_geom->tx_org_x[is_inter][context_ptr->tx_depth]
                                                               [context_ptr->txb_itr];
            uint16_t tx_org_y = context_ptr->blk_geom->tx_org_y[is_inter][context_ptr->tx_depth]
                                                               [context_ptr->txb_itr];
            uint32_t txb_origin_index = tx_org_x +
                (tx_org_y * tx_candidate_buffer->residual_ptr->stride_y);
            uint32_t input_txb_origin_index = (context_ptr->sb_origin_x + tx_org_x +
                                               input_picture_ptr->origin_x) +
                ((context_ptr->sb_origin_y + tx_org_y + input_picture_ptr->origin_y) *
                 input_picture_ptr->stride_y);

            // Y Prediction

            if (!is_inter) {
                if (context_ptr->tx_depth)
                    av1_intra_luma_prediction(context_ptr, pcs_ptr, tx_candidate_buffer);

                // Y Residual
                residual_kernel(
                    input_picture_ptr->buffer_y,
                    input_txb_origin_index,
                    input_picture_ptr->stride_y,
                    tx_candidate_buffer->prediction_ptr->buffer_y,
                    txb_origin_index,
                    tx_candidate_buffer->prediction_ptr->stride_y,
                    (int16_t *)tx_candidate_buffer->residual_ptr->buffer_y,
                    txb_origin_index,
                    tx_candidate_buffer->residual_ptr->stride_y,
                    context_ptr->hbd_mode_decision,
                    context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr]);
            }

            if (context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] <=
                    32 &&
                context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] <= 32)
                if (!tx_search_skip_flag) {
                    tx_type_search(pcs_ptr, context_ptr, tx_candidate_buffer, qp);
                }

            product_full_loop(tx_candidate_buffer,
                              context_ptr,
                              pcs_ptr,
                              input_picture_ptr,
                              context_ptr->blk_ptr->qp,
                              &(tx_y_count_non_zero_coeffs[0]),
                              &tx_y_coeff_bits,
                              &tx_y_full_distortion[0]);

            uint32_t y_has_coeff = tx_y_count_non_zero_coeffs[context_ptr->txb_itr] > 0;

            tx_update_neighbor_arrays(pcs_ptr, context_ptr, tx_candidate_buffer, is_inter);

            if (y_has_coeff)
                block_has_coeff = EB_TRUE;

            current_tx_cost = RDCOST(
                full_lambda, tx_y_coeff_bits, tx_y_full_distortion[DIST_CALC_RESIDUAL]);
            if (current_tx_cost > best_cost_search)
                break;

        } // Transform Loop

        if (end_tx_depth) {
            uint64_t tx_size_bits = 0;
            if (pcs_ptr->parent_pcs_ptr->frm_hdr.tx_mode == TX_MODE_SELECT)
                tx_size_bits = get_tx_size_bits(tx_candidate_buffer,
                                                context_ptr,
                                                pcs_ptr,
                                                context_ptr->tx_depth,
                                                block_has_coeff);

            uint64_t cost = RDCOST(full_lambda,
                                   (tx_y_coeff_bits + tx_size_bits),
                                   tx_y_full_distortion[DIST_CALC_RESIDUAL]);

            if (cost < best_cost_search) {
                best_cost_search                      = cost;
                best_tx_depth                         = context_ptr->tx_depth;
                is_best_has_coeff                     = block_has_coeff;
                y_full_distortion[DIST_CALC_RESIDUAL] = tx_y_full_distortion[DIST_CALC_RESIDUAL];
                y_full_distortion[DIST_CALC_PREDICTION] =
                    tx_y_full_distortion[DIST_CALC_PREDICTION];
                *y_coeff_bits = tx_y_coeff_bits;
                for (context_ptr->txb_itr = 0; context_ptr->txb_itr < txb_count;
                     context_ptr->txb_itr++) {
                    y_count_non_zero_coeffs[context_ptr->txb_itr] =
                        tx_y_count_non_zero_coeffs[context_ptr->txb_itr];
                }
            }
        } else {
            y_full_distortion[DIST_CALC_RESIDUAL]   = tx_y_full_distortion[DIST_CALC_RESIDUAL];
            y_full_distortion[DIST_CALC_PREDICTION] = tx_y_full_distortion[DIST_CALC_PREDICTION];
            *y_coeff_bits                           = tx_y_coeff_bits;
            for (context_ptr->txb_itr = 0; context_ptr->txb_itr < txb_count;
                 context_ptr->txb_itr++) {
                y_count_non_zero_coeffs[context_ptr->txb_itr] =
                    tx_y_count_non_zero_coeffs[context_ptr->txb_itr];
            }
        }
    } // Transform Depth Loop

    update_tx_candidate_buffer(candidate_buffer, context_ptr, best_tx_depth);
}

void full_loop_core(PictureControlSet *pcs_ptr, SuperBlock *sb_ptr, BlkStruct *blk_ptr,
                    ModeDecisionContext *context_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
                    ModeDecisionCandidate *candidate_ptr, EbPictureBufferDesc *input_picture_ptr,
                    uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                    uint32_t blk_origin_index, uint32_t blk_chroma_origin_index,
                    uint64_t ref_fast_cost) {
    uint64_t y_full_distortion[DIST_CALC_TOTAL];
    uint32_t count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];

    uint64_t cb_full_distortion[DIST_CALC_TOTAL];
    uint64_t cr_full_distortion[DIST_CALC_TOTAL];

    uint64_t y_coeff_bits;
    uint64_t cb_coeff_bits = 0;
    uint64_t cr_coeff_bits = 0;
    uint32_t full_lambda   = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];

    // initialize TU Split
    y_full_distortion[DIST_CALC_RESIDUAL]   = 0;
    y_full_distortion[DIST_CALC_PREDICTION] = 0;
    y_coeff_bits                            = 0;

    candidate_ptr->full_distortion = 0;

    memset(candidate_ptr->eob[0], 0, sizeof(uint16_t));
    memset(candidate_ptr->eob[1], 0, sizeof(uint16_t));
    memset(candidate_ptr->eob[2], 0, sizeof(uint16_t));

    candidate_ptr->chroma_distortion             = 0;
    candidate_ptr->chroma_distortion_inter_depth = 0;
    // Set Skip Flag
    candidate_ptr->skip_flag = EB_FALSE;

    if (candidate_ptr->type != INTRA_MODE) {
        if (context_ptr->md_staging_skip_full_pred == EB_FALSE) {
            product_prediction_fun_table[candidate_ptr->type](
                context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
        }
    }

    // Initialize luma CBF
    candidate_ptr->y_has_coeff = 0;
    candidate_ptr->u_has_coeff = 0;
    candidate_ptr->v_has_coeff = 0;

    // Initialize tx type
    for (int tu_index = 0; tu_index < MAX_TXB_COUNT; tu_index++)
        candidate_ptr->transform_type[tu_index] = DCT_DCT;
    uint8_t start_tx_depth = 0;
    uint8_t end_tx_depth   = 0;
    if (context_ptr->md_tx_size_search_mode == 0) {
        start_tx_depth = end_tx_depth = 0;
    } else if (context_ptr->md_staging_tx_size_mode == 0) {
        start_tx_depth = end_tx_depth = candidate_buffer->candidate_ptr->tx_depth;
    } else {
        // end_tx_depth set to zero for blocks which go beyond the picture boundaries
        if ((context_ptr->sb_origin_x + context_ptr->blk_geom->origin_x +
                     context_ptr->blk_geom->bwidth <
                 pcs_ptr->parent_pcs_ptr->aligned_width &&
             context_ptr->sb_origin_y + context_ptr->blk_geom->origin_y +
                     context_ptr->blk_geom->bheight <
                 pcs_ptr->parent_pcs_ptr->aligned_height))
            end_tx_depth = get_end_tx_depth(context_ptr->blk_geom->bsize);
        else
            end_tx_depth = 0;
    }
    // Transform partitioning path (INTRA Luma)
    int32_t is_inter = (candidate_buffer->candidate_ptr->type == INTER_MODE ||
                        candidate_buffer->candidate_ptr->use_intrabc)
        ? EB_TRUE
        : EB_FALSE;

    //Y Residual: residual for INTRA is computed inside the TU loop
    if (is_inter)
        //Y Residual
        residual_kernel(input_picture_ptr->buffer_y,
                        input_origin_index,
                        input_picture_ptr->stride_y,
                        candidate_buffer->prediction_ptr->buffer_y,
                        blk_origin_index,
                        candidate_buffer->prediction_ptr->stride_y,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_y,
                        blk_origin_index,
                        candidate_buffer->residual_ptr->stride_y,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->bwidth,
                        context_ptr->blk_geom->bheight);

    perform_tx_partitioning(candidate_buffer,
                            context_ptr,
                            pcs_ptr,
                            ref_fast_cost,
                            start_tx_depth,
                            end_tx_depth,
                            context_ptr->blk_ptr->qp,
                            &(*count_non_zero_coeffs[0]),
                            &y_coeff_bits,
                            &y_full_distortion[0]);

    candidate_ptr->chroma_distortion_inter_depth = 0;
    candidate_ptr->chroma_distortion             = 0;

    //CHROMA

    cb_full_distortion[DIST_CALC_RESIDUAL]   = 0;
    cr_full_distortion[DIST_CALC_RESIDUAL]   = 0;
    cb_full_distortion[DIST_CALC_PREDICTION] = 0;
    cr_full_distortion[DIST_CALC_PREDICTION] = 0;

    cb_coeff_bits = 0;
    cr_coeff_bits = 0;

    // FullLoop and TU search
    uint16_t cb_qp = context_ptr->qp;
    uint16_t cr_qp = context_ptr->qp;
    if (context_ptr->md_staging_skip_full_chroma == EB_FALSE) {
        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
            //Cb Residual
            residual_kernel(input_picture_ptr->buffer_cb,
                            input_cb_origin_in_index,
                            input_picture_ptr->stride_cb,
                            candidate_buffer->prediction_ptr->buffer_cb,
                            blk_chroma_origin_index,
                            candidate_buffer->prediction_ptr->stride_cb,
                            (int16_t *)candidate_buffer->residual_ptr->buffer_cb,
                            blk_chroma_origin_index,
                            candidate_buffer->residual_ptr->stride_cb,
                            context_ptr->hbd_mode_decision,
                            context_ptr->blk_geom->bwidth_uv,
                            context_ptr->blk_geom->bheight_uv);

            //Cr Residual
            residual_kernel(input_picture_ptr->buffer_cr,
                            input_cb_origin_in_index,
                            input_picture_ptr->stride_cr,
                            candidate_buffer->prediction_ptr->buffer_cr,
                            blk_chroma_origin_index,
                            candidate_buffer->prediction_ptr->stride_cr,
                            (int16_t *)candidate_buffer->residual_ptr->buffer_cr,
                            blk_chroma_origin_index,
                            candidate_buffer->residual_ptr->stride_cr,
                            context_ptr->hbd_mode_decision,
                            context_ptr->blk_geom->bwidth_uv,
                            context_ptr->blk_geom->bheight_uv);
        }

        if (!is_inter)
            if (candidate_buffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED)
                // If mode is CFL:
                // 1: recon the Luma
                // 2: Form the pred_buf_q3
                // 3: Loop over alphas and find the best or choose DC
                // 4: Recalculate the residual for chroma
                cfl_prediction(pcs_ptr,
                               candidate_buffer,
                               sb_ptr,
                               context_ptr,
                               input_picture_ptr,
                               input_cb_origin_in_index,
                               blk_chroma_origin_index);

        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
            full_loop_r(sb_ptr,
                        candidate_buffer,
                        context_ptr,
                        input_picture_ptr,
                        pcs_ptr,
                        PICTURE_BUFFER_DESC_CHROMA_MASK,
                        cb_qp,
                        cr_qp,
                        &(*count_non_zero_coeffs[1]),
                        &(*count_non_zero_coeffs[2]));

            cu_full_distortion_fast_txb_mode_r(sb_ptr,
                                               candidate_buffer,
                                               context_ptr,
                                               candidate_ptr,
                                               pcs_ptr,
                                               input_picture_ptr,
                                               cb_full_distortion,
                                               cr_full_distortion,
                                               count_non_zero_coeffs,
                                               COMPONENT_CHROMA,
                                               &cb_coeff_bits,
                                               &cr_coeff_bits,
                                               1);
        }

        // Check independant chroma vs. cfl
        if (!is_inter)
            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0)
                if (candidate_buffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED ||
                    candidate_buffer->candidate_ptr->intra_chroma_mode == UV_DC_PRED)
                    check_best_indepedant_cfl(pcs_ptr,
                                              input_picture_ptr,
                                              context_ptr,
                                              input_cb_origin_in_index,
                                              blk_chroma_origin_index,
                                              candidate_buffer,
                                              (uint8_t)cb_qp,
                                              (uint8_t)cr_qp,
                                              cb_full_distortion,
                                              cr_full_distortion,
                                              &cb_coeff_bits,
                                              &cr_coeff_bits);
    }

    candidate_ptr->block_has_coeff = (candidate_ptr->y_has_coeff | candidate_ptr->u_has_coeff |
                                      candidate_ptr->v_has_coeff)
        ? EB_TRUE
        : EB_FALSE;

    //ALL PLANE
    av1_product_full_cost_func_table[candidate_ptr->type](pcs_ptr,
                                                          context_ptr,
                                                          candidate_buffer,
                                                          blk_ptr,
                                                          y_full_distortion,
                                                          cb_full_distortion,
                                                          cr_full_distortion,
                                                          full_lambda,
                                                          &y_coeff_bits,
                                                          &cb_coeff_bits,
                                                          &cr_coeff_bits,
                                                          context_ptr->blk_geom->bsize);
}
void md_stage_1(PictureControlSet *pcs_ptr, SuperBlock *sb_ptr, BlkStruct *blk_ptr,
                ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                uint32_t blk_origin_index, uint32_t blk_chroma_origin_index,
                uint64_t ref_fast_cost) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array = &(
        candidate_buffer_ptr_array_base[0]);
    ModeDecisionCandidateBuffer *candidate_buffer;
    ModeDecisionCandidate *      candidate_ptr;

    uint32_t full_loop_candidate_index;
    uint32_t cand_index;

    // Set MD Staging full_loop_core settings
    context_ptr->md_staging_tx_size_mode          = 0;
    context_ptr->md_staging_tx_search             = 0;
    context_ptr->md_staging_skip_full_chroma      = EB_TRUE;
    context_ptr->md_staging_skip_rdoq             = EB_TRUE;
    context_ptr->md_staging_spatial_sse_full_loop = EB_FALSE;
    for (full_loop_candidate_index = 0;
         full_loop_candidate_index < context_ptr->md_stage_1_count[context_ptr->target_class];
         ++full_loop_candidate_index) {
        cand_index =
            context_ptr->cand_buff_indices[context_ptr->target_class][full_loop_candidate_index];
        candidate_buffer = candidate_buffer_ptr_array[cand_index];
        candidate_ptr    = candidate_buffer->candidate_ptr;

        context_ptr->md_staging_skip_full_pred            = EB_FALSE;
        context_ptr->md_staging_skip_interpolation_search = EB_TRUE;
        context_ptr->md_staging_skip_inter_chroma_pred    = EB_TRUE;
        candidate_buffer->candidate_ptr->interp_filters   = 0;
        full_loop_core(pcs_ptr,
                       sb_ptr,
                       blk_ptr,
                       context_ptr,
                       candidate_buffer,
                       candidate_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index,
                       ref_fast_cost);
    }
}

void md_stage_2(PictureControlSet *pcs_ptr, SuperBlock *sb_ptr, BlkStruct *blk_ptr,
                ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                uint32_t blk_origin_index, uint32_t blk_chroma_origin_index,
                uint64_t ref_fast_cost) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array = &(
        candidate_buffer_ptr_array_base[0]);
    ModeDecisionCandidateBuffer *candidate_buffer;
    ModeDecisionCandidate *      candidate_ptr;

    uint32_t fullLoopCandidateIndex;
    uint32_t candidateIndex;

    // Set MD Staging full_loop_core settings
    for (fullLoopCandidateIndex = 0;
         fullLoopCandidateIndex < context_ptr->md_stage_2_count[context_ptr->target_class];
         ++fullLoopCandidateIndex) {
        candidateIndex =
            context_ptr->cand_buff_indices[context_ptr->target_class][fullLoopCandidateIndex];
        candidate_buffer = candidate_buffer_ptr_array[candidateIndex];
        candidate_ptr    = candidate_buffer->candidate_ptr;

        context_ptr->md_staging_tx_size_mode = 0;

        context_ptr->md_staging_tx_search = (candidate_ptr->cand_class == CAND_CLASS_0 ||
                                             candidate_ptr->cand_class == CAND_CLASS_6 ||
                                             candidate_ptr->cand_class == CAND_CLASS_7)
            ? 2
            : 1;
        context_ptr->md_staging_skip_rdoq                 = EB_FALSE;
        context_ptr->md_staging_skip_full_chroma          = EB_TRUE;
        context_ptr->md_staging_skip_full_pred            = EB_TRUE;
        context_ptr->md_staging_skip_interpolation_search = EB_TRUE;
        context_ptr->md_staging_skip_inter_chroma_pred    = EB_TRUE;
        context_ptr->md_staging_spatial_sse_full_loop     = EB_FALSE;
        full_loop_core(pcs_ptr,
                       sb_ptr,
                       blk_ptr,
                       context_ptr,
                       candidate_buffer,
                       candidate_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index,
                       ref_fast_cost);
    }
}

void md_stage_3(PictureControlSet *pcs_ptr, SuperBlock *sb_ptr, BlkStruct *blk_ptr,
                ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                uint32_t blk_origin_index, uint32_t blk_chroma_origin_index,
                uint32_t fullCandidateTotalCount, uint64_t ref_fast_cost) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array = &(
        candidate_buffer_ptr_array_base[0]);
    ModeDecisionCandidateBuffer *candidate_buffer;
    ModeDecisionCandidate *      candidate_ptr;

    uint32_t best_inter_luma_zero_coeff = 1;
    uint64_t best_full_cost             = 0xFFFFFFFFull;
    uint32_t full_loop_candidate_index;
    uint32_t cand_index;

    for (full_loop_candidate_index = 0; full_loop_candidate_index < fullCandidateTotalCount;
         ++full_loop_candidate_index) {
        cand_index = (context_ptr->full_loop_escape == 2)
            ? context_ptr->sorted_candidate_index_array[full_loop_candidate_index]
            : context_ptr->best_candidate_index_array[full_loop_candidate_index];
        candidate_buffer = candidate_buffer_ptr_array[cand_index];
        candidate_ptr    = candidate_buffer->candidate_ptr;

        // Set MD Staging full_loop_core settings
        context_ptr->md_staging_skip_full_pred = context_ptr->md_staging_mode == MD_STAGING_MODE_0;
        context_ptr->md_staging_skip_interpolation_search = EB_FALSE;
        context_ptr->md_staging_skip_inter_chroma_pred    = EB_FALSE;
        // only perform TXS for intra frames
        context_ptr->md_staging_tx_size_mode = (MR_MODE ||
                                                candidate_ptr->cand_class == CAND_CLASS_0 ||
                                                candidate_ptr->cand_class == CAND_CLASS_6 ||
                                                candidate_ptr->cand_class == CAND_CLASS_7)
            ? 1
            : 0;
        context_ptr->md_staging_tx_search = (candidate_ptr->cand_class == CAND_CLASS_0 ||
                                             candidate_ptr->cand_class == CAND_CLASS_6 ||
                                             candidate_ptr->cand_class == CAND_CLASS_7)
            ? 2
            : 1;
        context_ptr->md_staging_skip_full_chroma = EB_FALSE;

        context_ptr->md_staging_skip_rdoq             = EB_FALSE;
        context_ptr->md_staging_spatial_sse_full_loop = context_ptr->spatial_sse_full_loop;

        if (pcs_ptr->slice_type != I_SLICE) {
            if ((candidate_ptr->type == INTRA_MODE || context_ptr->full_loop_escape == 2) &&
                best_inter_luma_zero_coeff == 0) {
                context_ptr->md_stage_3_total_count = full_loop_candidate_index;
                return;
            }
        }

        if (context_ptr->chroma_at_last_md_stage) {
            if (context_ptr->blk_geom->sq_size < 128) {
                if (context_ptr->blk_geom->has_uv) {
                    if (candidate_ptr->type == INTRA_MODE) {
                        // skip CFL based on inter/intra cost deviation; check that th != INF to avoid overflow errors
                        uint32_t intra_chroma_mode;
                        int32_t  angle_delta;
                        uint8_t  is_directional_chroma_mode_flag;
                        if ((context_ptr->skip_cfl_cost_dev_th != (uint16_t)~0) &&
                            ((context_ptr->best_inter_cost *
                              (100 + context_ptr->skip_cfl_cost_dev_th)) <
                             (context_ptr->best_intra_cost * 100)) &&
                            !(pcs_ptr->parent_pcs_ptr->sc_content_detected)) {
                            intra_chroma_mode =
                                context_ptr->best_uv_mode[candidate_ptr->intra_luma_mode]
                                                         [MAX_ANGLE_DELTA +
                                                          candidate_ptr->angle_delta[PLANE_TYPE_Y]];
                            angle_delta =
                                context_ptr
                                    ->best_uv_angle[candidate_ptr->intra_luma_mode]
                                                   [MAX_ANGLE_DELTA +
                                                    candidate_ptr->angle_delta[PLANE_TYPE_Y]];
                            is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode(
                                (PredictionMode)(
                                    context_ptr
                                        ->best_uv_mode[candidate_ptr->intra_luma_mode]
                                                      [MAX_ANGLE_DELTA +
                                                       candidate_ptr->angle_delta[PLANE_TYPE_Y]]));
                        } else {
                            intra_chroma_mode = candidate_ptr->intra_chroma_mode != UV_CFL_PRED
                                ? context_ptr
                                      ->best_uv_mode[candidate_ptr->intra_luma_mode]
                                                    [MAX_ANGLE_DELTA +
                                                     candidate_ptr->angle_delta[PLANE_TYPE_Y]]
                                : UV_CFL_PRED;
                            angle_delta = candidate_ptr->intra_chroma_mode != UV_CFL_PRED
                                ? context_ptr
                                      ->best_uv_angle[candidate_ptr->intra_luma_mode]
                                                     [MAX_ANGLE_DELTA +
                                                      candidate_ptr->angle_delta[PLANE_TYPE_Y]]
                                : 0;
                            is_directional_chroma_mode_flag = candidate_ptr->intra_chroma_mode !=
                                    UV_CFL_PRED
                                ? (uint8_t)av1_is_directional_mode((PredictionMode)(
                                      context_ptr
                                          ->best_uv_mode[candidate_ptr->intra_luma_mode]
                                                        [MAX_ANGLE_DELTA +
                                                         candidate_ptr->angle_delta[PLANE_TYPE_Y]]))
                                : 0;
                        }
                        candidate_ptr->intra_chroma_mode          = intra_chroma_mode;
                        candidate_ptr->angle_delta[PLANE_TYPE_UV] = angle_delta;
                        candidate_ptr->is_directional_chroma_mode_flag =
                            is_directional_chroma_mode_flag;
                    }
                }
            }
        }
        full_loop_core(pcs_ptr,
                       sb_ptr,
                       blk_ptr,
                       context_ptr,
                       candidate_buffer,
                       candidate_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index,
                       ref_fast_cost);

        if (context_ptr->full_loop_escape) {
            if (pcs_ptr->slice_type != I_SLICE) {
                if (candidate_ptr->type == INTER_MODE) {
                    if (*candidate_buffer->full_cost_ptr < best_full_cost) {
                        best_inter_luma_zero_coeff = candidate_ptr->y_has_coeff;
                        best_full_cost             = *candidate_buffer->full_cost_ptr;
                    }
                }
            }
        }
    }
}

void move_blk_data(PictureControlSet *pcs, EncDecContext *context_ptr, BlkStruct *src_cu,
                   BlkStruct *dst_cu) {
    memcpy(&dst_cu->palette_info.pmi, &src_cu->palette_info.pmi, sizeof(PaletteModeInfo));
    if (svt_av1_allow_palette(pcs->parent_pcs_ptr->palette_mode, context_ptr->blk_geom->bsize)) {
        dst_cu->palette_info.color_idx_map = (uint8_t *)malloc(MAX_PALETTE_SQUARE);
        assert(dst_cu->palette_info.color_idx_map != NULL && "palette:Not-Enough-Memory");
        if (dst_cu->palette_info.color_idx_map != NULL)
            memcpy(dst_cu->palette_info.color_idx_map,
                   src_cu->palette_info.color_idx_map,
                   MAX_PALETTE_SQUARE);
        else
            SVT_LOG("ERROR palette:Not-Enough-Memory\n");
    }
    dst_cu->interp_filters              = src_cu->interp_filters;
    dst_cu->interinter_comp.type        = src_cu->interinter_comp.type;
    dst_cu->interinter_comp.mask_type   = src_cu->interinter_comp.mask_type;
    dst_cu->interinter_comp.wedge_index = src_cu->interinter_comp.wedge_index;
    dst_cu->interinter_comp.wedge_sign  = src_cu->interinter_comp.wedge_sign;
    dst_cu->compound_idx                = src_cu->compound_idx;
    dst_cu->comp_group_idx              = src_cu->comp_group_idx;

    dst_cu->is_interintra_used     = src_cu->is_interintra_used;
    dst_cu->interintra_mode        = src_cu->interintra_mode;
    dst_cu->use_wedge_interintra   = src_cu->use_wedge_interintra;
    dst_cu->interintra_wedge_index = src_cu->interintra_wedge_index; //inter_intra wedge index
    //CHKN TransformUnit             txb_array[TRANSFORM_UNIT_MAX_COUNT]; // 2-bytes * 21 = 42-bytes
    memcpy(dst_cu->txb_array, src_cu->txb_array, TRANSFORM_UNIT_MAX_COUNT * sizeof(TransformUnit));

    //CHKN PredictionUnit            prediction_unit_array[MAX_NUM_OF_PU_PER_CU];    // 35-bytes * 4 = 140 bytes
    memcpy(dst_cu->prediction_unit_array,
           src_cu->prediction_unit_array,
           MAX_NUM_OF_PU_PER_CU * sizeof(PredictionUnit));

    dst_cu->skip_flag_context    = src_cu->skip_flag_context;
    dst_cu->prediction_mode_flag = src_cu->prediction_mode_flag;
    dst_cu->block_has_coeff      = src_cu->block_has_coeff;
    dst_cu->split_flag_context   = src_cu->split_flag_context;
    dst_cu->qp                   = src_cu->qp;
    dst_cu->tx_depth             = src_cu->tx_depth;
    dst_cu->split_flag           = src_cu->split_flag;
    dst_cu->skip_flag            = src_cu->skip_flag;

    //CHKN    MacroBlockD*  av1xd;
    memcpy(dst_cu->av1xd, src_cu->av1xd, sizeof(MacroBlockD));

    // uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];

    //CHKN int16_t inter_mode_ctx[MODE_CTX_REF_FRAMES];
    memcpy(dst_cu->inter_mode_ctx, src_cu->inter_mode_ctx, MODE_CTX_REF_FRAMES * sizeof(int16_t));

    //CHKN uint8_t  drl_index;
    //CHKN PredictionMode               pred_mode;
    dst_cu->drl_index = src_cu->drl_index;
    dst_cu->pred_mode = src_cu->pred_mode;

    //CHKN IntMv  predmv[2];

    memcpy(dst_cu->predmv, src_cu->predmv, 2 * sizeof(IntMv));
    //CHKN uint8_t                         skip_coeff_context;
    //CHKN int16_t                        luma_txb_skip_context;
    //CHKN int16_t                        luma_dc_sign_context;
    //CHKN int16_t                        cb_txb_skip_context;
    //CHKN int16_t                        cb_dc_sign_context;
    //CHKN int16_t                        cr_txb_skip_context;
    //CHKN int16_t                        cr_dc_sign_context;
    //CHKN uint8_t                         reference_mode_context;
    //CHKN uint8_t                         compoud_reference_type_context;
    //CHKN uint32_t                        partitionContext;

    dst_cu->skip_coeff_context             = src_cu->skip_coeff_context;
    dst_cu->reference_mode_context         = src_cu->reference_mode_context;
    dst_cu->compoud_reference_type_context = src_cu->compoud_reference_type_context;
    dst_cu->segment_id                     = src_cu->segment_id;

    memcpy(dst_cu->quantized_dc, src_cu->quantized_dc, 3 * MAX_TXB_COUNT * sizeof(int32_t));
    //CHKN uint32_t   is_inter_ctx;
    //CHKN uint32_t                     interp_filters;

    dst_cu->is_inter_ctx   = src_cu->is_inter_ctx;
    dst_cu->interp_filters = src_cu->interp_filters;

    dst_cu->part              = src_cu->part;
    dst_cu->mds_idx           = src_cu->mds_idx;
    dst_cu->filter_intra_mode = src_cu->filter_intra_mode;
}
void move_blk_data_redund(PictureControlSet *pcs, ModeDecisionContext *context_ptr,
                          BlkStruct *src_cu, BlkStruct *dst_cu) {
    dst_cu->segment_id       = src_cu->segment_id;
    dst_cu->seg_id_predicted = src_cu->seg_id_predicted;
    memcpy(&dst_cu->palette_info.pmi, &src_cu->palette_info.pmi, sizeof(PaletteModeInfo));
    if (svt_av1_allow_palette(pcs->parent_pcs_ptr->palette_mode, context_ptr->blk_geom->bsize))
        memcpy(dst_cu->palette_info.color_idx_map,
               src_cu->palette_info.color_idx_map,
               MAX_PALETTE_SQUARE);
    dst_cu->interp_filters              = src_cu->interp_filters;
    dst_cu->interinter_comp.type        = src_cu->interinter_comp.type;
    dst_cu->interinter_comp.mask_type   = src_cu->interinter_comp.mask_type;
    dst_cu->interinter_comp.wedge_index = src_cu->interinter_comp.wedge_index;
    dst_cu->interinter_comp.wedge_sign  = src_cu->interinter_comp.wedge_sign;
    dst_cu->compound_idx                = src_cu->compound_idx;
    dst_cu->comp_group_idx              = src_cu->comp_group_idx;
    dst_cu->is_interintra_used          = src_cu->is_interintra_used;
    dst_cu->interintra_mode             = src_cu->interintra_mode;
    dst_cu->use_wedge_interintra        = src_cu->use_wedge_interintra;
    dst_cu->interintra_wedge_index      = src_cu->interintra_wedge_index; //inter_intra wedge index
    dst_cu->filter_intra_mode           = src_cu->filter_intra_mode;
    //CHKN TransformUnit_t             txb_array[TRANSFORM_UNIT_MAX_COUNT]; // 2-bytes * 21 = 42-bytes
    memcpy(dst_cu->txb_array, src_cu->txb_array, TRANSFORM_UNIT_MAX_COUNT * sizeof(TransformUnit));

    //CHKN PredictionUnit_t            prediction_unit_array[MAX_NUM_OF_PU_PER_CU];    // 35-bytes * 4 = 140 bytes
    memcpy(dst_cu->prediction_unit_array,
           src_cu->prediction_unit_array,
           MAX_NUM_OF_PU_PER_CU * sizeof(PredictionUnit));
    dst_cu->skip_flag_context    = src_cu->skip_flag_context;
    dst_cu->prediction_mode_flag = src_cu->prediction_mode_flag;
    dst_cu->block_has_coeff      = src_cu->block_has_coeff;
    dst_cu->split_flag_context   = src_cu->split_flag_context;
    dst_cu->qp                   = src_cu->qp;
    dst_cu->skip_flag            = src_cu->skip_flag;
    dst_cu->tx_depth             = src_cu->tx_depth;
    //CHKN    MacroBlockD*  av1xd;
    memcpy(dst_cu->av1xd, src_cu->av1xd, sizeof(MacroBlockD));

    // uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];

    //CHKN int16_t inter_mode_ctx[MODE_CTX_REF_FRAMES];
    memcpy(dst_cu->inter_mode_ctx, src_cu->inter_mode_ctx, MODE_CTX_REF_FRAMES * sizeof(int16_t));

    //CHKN uint8_t  drl_index;
    //CHKN PredictionMode               pred_mode;
    dst_cu->drl_index = src_cu->drl_index;
    dst_cu->pred_mode = src_cu->pred_mode;

    //CHKN IntMv  predmv[2];

    memcpy(dst_cu->predmv, src_cu->predmv, 2 * sizeof(IntMv));

    //CHKN uint8_t                         skip_coeff_context;
    //CHKN int16_t                        luma_txb_skip_context;
    //CHKN int16_t                        luma_dc_sign_context;
    //CHKN int16_t                        cb_txb_skip_context;
    //CHKN int16_t                        cb_dc_sign_context;
    //CHKN int16_t                        cr_txb_skip_context;
    //CHKN int16_t                        cr_dc_sign_context;
    //CHKN uint8_t                         reference_mode_context;
    //CHKN uint8_t                         compoud_reference_type_context;
    //CHKN uint32_t                        partitionContext;

    dst_cu->skip_coeff_context             = src_cu->skip_coeff_context;
    dst_cu->reference_mode_context         = src_cu->reference_mode_context;
    dst_cu->compoud_reference_type_context = src_cu->compoud_reference_type_context;
    memcpy(dst_cu->quantized_dc, src_cu->quantized_dc, 3 * MAX_TXB_COUNT * sizeof(int32_t));
    //CHKN uint32_t   is_inter_ctx;
    //CHKN uint32_t                     interp_filters;

    dst_cu->is_inter_ctx   = src_cu->is_inter_ctx;
    dst_cu->interp_filters = src_cu->interp_filters;

    dst_cu->part = src_cu->part;
}

void check_redundant_block(const BlockGeom *blk_geom, ModeDecisionContext *context_ptr,
                           uint8_t *redundant_blk_avail, uint16_t *redundant_blk_mds) {
    if (blk_geom->redund) {
        for (int it = 0; it < blk_geom->redund_list.list_size; it++) {
            if (context_ptr->md_local_blk_unit[blk_geom->redund_list.blk_mds_table[it]]
                    .avail_blk_flag) {
                *redundant_blk_mds   = blk_geom->redund_list.blk_mds_table[it];
                *redundant_blk_avail = 1;
                break;
            }
        }
    }
}

/*
   search for a valid previously encoded similar
   block (block having the same location and shape as the current block,
   but where neighboring blocks are different from those for the current block)
*/
void check_similar_block(const BlockGeom *blk_geom, ModeDecisionContext *context_ptr,
                         uint8_t *similar_blk_avail, uint16_t *similar_blk_mds) {
    if (blk_geom->similar) {
        for (int it = 0; it < blk_geom->similar_list.list_size; it++) {
            if (context_ptr->md_local_blk_unit[blk_geom->similar_list.blk_mds_table[it]]
                    .avail_blk_flag) {
                *similar_blk_mds   = blk_geom->similar_list.blk_mds_table[it];
                *similar_blk_avail = 1;
                break;
            }
        }
    }
}
/******************************************************
* Derive md Settings(feature signals) that could be
  changed  at the block level
******************************************************/
EbErrorType signal_derivation_block(PictureControlSet *pcs, ModeDecisionContext *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    // set compound_types_to_try
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->compound_types_to_try = MD_COMP_AVG;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->compound_types_to_try = MD_COMP_AVG;
    else {
        if (pcs->parent_pcs_ptr->compound_mode)
            context_ptr->compound_types_to_try = pcs->parent_pcs_ptr->compound_mode == 1
                ? MD_COMP_DIFF0
                : MD_COMP_WEDGE;
        else
            context_ptr->compound_types_to_try = MD_COMP_AVG;
    }

    BlkStruct *similar_cu = &context_ptr->md_blk_arr_nsq[context_ptr->similar_blk_mds];
    if (context_ptr->compound_types_to_try > MD_COMP_AVG && context_ptr->similar_blk_avail) {
        int32_t is_src_compound = similar_cu->pred_mode >= NEAREST_NEARESTMV;
        if (context_ptr->comp_similar_mode == 1) {
            context_ptr->compound_types_to_try = !is_src_compound
                ? MD_COMP_AVG
                : context_ptr->compound_types_to_try;
        } else if (context_ptr->comp_similar_mode == 2) {
            context_ptr->compound_types_to_try = !is_src_compound
                ? MD_COMP_AVG
                : similar_cu->interinter_comp.type;
        }
    }

    context_ptr->inject_inter_candidates = 1;
    if (context_ptr->pd_pass > PD_PASS_1 && context_ptr->similar_blk_avail) {
        int32_t is_src_intra = similar_cu->pred_mode <= PAETH_PRED;
        if (context_ptr->intra_similar_mode)
            context_ptr->inject_inter_candidates = is_src_intra
                ? 0
                : context_ptr->inject_inter_candidates;
    }

    return return_error;
}
/****************************************************
* generate the the size in pixel for partition code
****************************************************/
uint8_t get_part_side(PartitionContextType part) {
    switch (part) {
    case 31: return 4; break;
    case 30: return 8; break;
    case 28: return 16; break;
    case 24: return 32; break;
    case 16: return 64; break;
    case 0: return 128; break;
    default:
        return 255;
        SVT_LOG("error: non supported partition!!\n");
        break;
    }
}
/****************************************************
* Return a predicted Shape based on the above and
* left partitions
****************************************************/
Part get_partition_shape(PartitionContextType above, PartitionContextType left, uint8_t width,
                         uint8_t height) {
    uint8_t above_size = get_part_side(above);
    uint8_t left_size  = get_part_side(left);
    Part    part       = PART_N;

    if (above_size == width && left_size == height)
        part = PART_N;
    else if (above_size > width && left_size > height)
        part = PART_N;
    else if (above_size > width) {
        if (left_size == height)
            part = PART_N;
        else if (left_size < (height / 2))
            part = PART_H4;
        else if (left_size < height)
            part = PART_H;
        else
            SVT_LOG("error: unsupported left_size\n");
    } else if (left_size > height) {
        if (above_size == width)
            part = PART_N;
        else if (above_size < (width / 2))
            part = PART_V4;
        else if (above_size < width)
            part = PART_V;
        else
            SVT_LOG("error: unsupported above_size\n");
    } else if (above_size < width) {
        if (left_size == height)
            part = PART_VA;
        else if (left_size < height)
            part = PART_S;
        else
            SVT_LOG("error: unsupported left_size\n");
    } else if (left_size < height) {
        if (above_size == width)
            part = PART_HA;
        else if (above_size < width)
            part = PART_S;
        else
            SVT_LOG("error: unsupported above_size\n");
    } else if (above_size == width) {
        if (left_size < height)
            part = PART_HB;
        else
            SVT_LOG("error: unsupported left_size\n");
    } else if (left_size == height) {
        if (above_size == width)
            part = PART_HB;
        else
            SVT_LOG("error: unsupported above_size\n");
    } else
        SVT_LOG("error: unsupported above_size && left_size\n");
    return part;
};
void init_chroma_mode(ModeDecisionContext *context_ptr) {
    context_ptr->uv_search_path = EB_TRUE;
    EbBool use_angle_delta      = av1_use_angle_delta(context_ptr->blk_geom->bsize,
                                                 context_ptr->md_intra_angle_delta);
    for (uint8_t intra_mode = DC_PRED; intra_mode <= PAETH_PRED; ++intra_mode) {
        uint8_t angleDeltaCandidateCount = (use_angle_delta &&
                                            av1_is_directional_mode((PredictionMode)intra_mode))
            ? 7
            : 1;
        uint8_t angle_delta_shift = 1;
        for (uint8_t angleDeltaCounter = 0; angleDeltaCounter < angleDeltaCandidateCount;
             ++angleDeltaCounter) {
            int32_t angle_delta = CLIP(angle_delta_shift *
                                           (angleDeltaCandidateCount == 1 ? 0
                                                                          : angleDeltaCounter -
                                                    (angleDeltaCandidateCount >> 1)),
                                       -MAX_ANGLE_DELTA,
                                       MAX_ANGLE_DELTA);
            context_ptr->best_uv_mode[intra_mode][MAX_ANGLE_DELTA + angle_delta]  = intra_mode;
            context_ptr->best_uv_angle[intra_mode][MAX_ANGLE_DELTA + angle_delta] = angle_delta;
            context_ptr->best_uv_cost[intra_mode][MAX_ANGLE_DELTA + angle_delta]  = (uint64_t)~0;
        }
    }
    // End uv search path
    context_ptr->uv_search_path = EB_FALSE;
}
void search_best_independent_uv_mode(PictureControlSet *  pcs_ptr,
                                     EbPictureBufferDesc *input_picture_ptr,
                                     uint32_t             input_cb_origin_in_index,
                                     uint32_t             input_cr_origin_in_index,
                                     uint32_t             cu_chroma_origin_index,
                                     ModeDecisionContext *context_ptr) {
    FrameHeader *frm_hdr     = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    uint32_t     full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];
    // Start uv search path
    context_ptr->uv_search_path = EB_TRUE;

    EbBool use_angle_delta = av1_use_angle_delta(context_ptr->blk_geom->bsize,
                                                 context_ptr->md_intra_angle_delta);

    UvPredictionMode uv_mode;

    int coeff_rate[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
    int distortion[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];

    ModeDecisionCandidate *candidate_array          = context_ptr->fast_candidate_array;
    uint32_t               start_fast_buffer_index  = MODE_DECISION_CANDIDATE_MAX_COUNT_Y;
    uint32_t               start_full_buffer_index  = MAX_NFL_BUFF_Y;
    uint32_t               uv_mode_total_count      = start_fast_buffer_index;
    EbBool                 tem_md_staging_skip_rdoq = context_ptr->md_staging_skip_rdoq;
    if (context_ptr->chroma_at_last_md_stage) {
        context_ptr->md_staging_skip_rdoq = 0;
    }
    UvPredictionMode uv_mode_end = context_ptr->md_enable_paeth
        ? UV_PAETH_PRED
        : context_ptr->md_enable_smooth ? UV_SMOOTH_H_PRED : UV_D67_PRED;

    for (uv_mode = UV_DC_PRED; uv_mode <= uv_mode_end; uv_mode++) {
        uint8_t uv_angle_delta_candidate_count = (use_angle_delta &&
                                                  av1_is_directional_mode((PredictionMode)uv_mode))
            ? 7
            : 1;
        uint8_t uv_angle_delta_shift = 1;

        for (uint8_t uv_angle_delta_counter = 0;
             uv_angle_delta_counter < uv_angle_delta_candidate_count;
             ++uv_angle_delta_counter) {
            int32_t uv_angle_delta = CLIP(
                uv_angle_delta_shift *
                    (uv_angle_delta_candidate_count == 1
                         ? 0
                         : uv_angle_delta_counter - (uv_angle_delta_candidate_count >> 1)),
                -MAX_ANGLE_DELTA,
                MAX_ANGLE_DELTA);

            candidate_array[uv_mode_total_count].type                            = INTRA_MODE;
            candidate_array[uv_mode_total_count].distortion_ready                = 0;
            candidate_array[uv_mode_total_count].use_intrabc                     = 0;
            candidate_array[uv_mode_total_count].angle_delta[PLANE_TYPE_UV]      = 0;
            candidate_array[uv_mode_total_count].pred_mode                       = DC_PRED;
            candidate_array[uv_mode_total_count].intra_chroma_mode               = uv_mode;
            candidate_array[uv_mode_total_count].is_directional_chroma_mode_flag = (uint8_t)
                av1_is_directional_mode((PredictionMode)uv_mode);
            candidate_array[uv_mode_total_count].angle_delta[PLANE_TYPE_UV]       = uv_angle_delta;
            candidate_array[uv_mode_total_count].tx_depth                         = 0;
            candidate_array[uv_mode_total_count].palette_info.pmi.palette_size[0] = 0;
            candidate_array[uv_mode_total_count].palette_info.pmi.palette_size[1] = 0;
            candidate_array[uv_mode_total_count].filter_intra_mode = FILTER_INTRA_MODES;
            candidate_array[uv_mode_total_count].cfl_alpha_signs   = 0;
            candidate_array[uv_mode_total_count].cfl_alpha_idx     = 0;
            candidate_array[uv_mode_total_count].transform_type[0] = DCT_DCT;
            candidate_array[uv_mode_total_count].ref_frame_type    = INTRA_FRAME;
            candidate_array[uv_mode_total_count].motion_mode       = SIMPLE_TRANSLATION;

            candidate_array[uv_mode_total_count].transform_type_uv = av1_get_tx_type(
                context_ptr->blk_geom->bsize,
                0,
                (PredictionMode)NULL,
                (UvPredictionMode)uv_mode,
                PLANE_TYPE_UV,
                0,
                0,
                0,
                context_ptr->blk_geom->txsize_uv[0][0],
                frm_hdr->reduced_tx_set);

            uv_mode_total_count++;
        }
    }
    uv_mode_total_count = uv_mode_total_count - start_fast_buffer_index;
    // Fast-loop search uv_mode
    for (uint8_t uv_mode_count = 0; uv_mode_count < uv_mode_total_count; uv_mode_count++) {
        ModeDecisionCandidateBuffer *candidate_buffer =
            context_ptr->candidate_buffer_ptr_array[uv_mode_count + start_full_buffer_index];
        candidate_buffer->candidate_ptr =
            &context_ptr->fast_candidate_array[uv_mode_count + start_fast_buffer_index];

        context_ptr->md_staging_skip_inter_chroma_pred = EB_FALSE;
        product_prediction_fun_table[candidate_buffer->candidate_ptr->type](
            context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);

        uint32_t chroma_fast_distortion = 0;
        if (!context_ptr->hbd_mode_decision) {
            chroma_fast_distortion = nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_cb + input_cb_origin_in_index,
                input_picture_ptr->stride_cb,
                candidate_buffer->prediction_ptr->buffer_cb + cu_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cb,
                context_ptr->blk_geom->bheight_uv,
                context_ptr->blk_geom->bwidth_uv);

            chroma_fast_distortion += nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_cr + input_cr_origin_in_index,
                input_picture_ptr->stride_cr,
                candidate_buffer->prediction_ptr->buffer_cr + cu_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cr,
                context_ptr->blk_geom->bheight_uv,
                context_ptr->blk_geom->bwidth_uv);
        } else {
            chroma_fast_distortion = sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_cb) + input_cb_origin_in_index,
                input_picture_ptr->stride_cb,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cb) + cu_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cb,
                context_ptr->blk_geom->bheight_uv,
                context_ptr->blk_geom->bwidth_uv);

            chroma_fast_distortion += sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_cr) + input_cr_origin_in_index,
                input_picture_ptr->stride_cr,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cr) + cu_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cr,
                context_ptr->blk_geom->bheight_uv,
                context_ptr->blk_geom->bwidth_uv);
        }
        // Do not consider rate @ this stage
        *(candidate_buffer->fast_cost_ptr) = chroma_fast_distortion;
    }

    // Sort uv_mode (in terms of distortion only)
    uint32_t uv_cand_buff_indices[MAX_NFL_BUFF_Y];
    memset(uv_cand_buff_indices, 0xFFFFFFFF, MAX_NFL_BUFF_Y * sizeof(uint32_t));
    sort_fast_cost_based_candidates(
        context_ptr,
        start_full_buffer_index,
        uv_mode_total_count, //how many cand buffers to sort. one of the buffers can have max cost.
        uv_cand_buff_indices);

    // Reset *(candidate_buffer->fast_cost_ptr)
    for (uint8_t uv_mode_count = 0; uv_mode_count < uv_mode_total_count; uv_mode_count++) {
        ModeDecisionCandidateBuffer *candidate_buffer =
            context_ptr->candidate_buffer_ptr_array[uv_mode_count + start_full_buffer_index];
        *(candidate_buffer->fast_cost_ptr) = MAX_CU_COST;
    }

    // Derive uv_mode_nfl_count
    uint8_t uv_mode_nfl_count;
    if (pcs_ptr->temporal_layer_index == 0)
        uv_mode_nfl_count = uv_mode_total_count;
    else if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
        uv_mode_nfl_count = 16;
    else
        uv_mode_nfl_count = 8;

    // Full-loop search uv_mode
    for (uint8_t uv_mode_count = 0; uv_mode_count < MIN(uv_mode_total_count, uv_mode_nfl_count);
         uv_mode_count++) {
        ModeDecisionCandidateBuffer *candidate_buffer =
            context_ptr->candidate_buffer_ptr_array[uv_cand_buff_indices[uv_mode_count]];
        candidate_buffer->candidate_ptr =
            &context_ptr->fast_candidate_array[uv_cand_buff_indices[uv_mode_count] -
                                               start_full_buffer_index + start_fast_buffer_index];
        uint16_t cb_qp                               = context_ptr->qp;
        uint16_t cr_qp                               = context_ptr->qp;
        uint64_t cb_coeff_bits                       = 0;
        uint64_t cr_coeff_bits                       = 0;
        uint64_t cb_full_distortion[DIST_CALC_TOTAL] = {0, 0};
        uint64_t cr_full_distortion[DIST_CALC_TOTAL] = {0, 0};

        uint32_t count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];

        //Cb Residual
        residual_kernel(input_picture_ptr->buffer_cb,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cb,
                        candidate_buffer->prediction_ptr->buffer_cb,
                        cu_chroma_origin_index,
                        candidate_buffer->prediction_ptr->stride_cb,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cb,
                        cu_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cb,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->bwidth_uv,
                        context_ptr->blk_geom->bheight_uv);

        //Cr Residual
        residual_kernel(input_picture_ptr->buffer_cr,
                        input_cr_origin_in_index,
                        input_picture_ptr->stride_cr,
                        candidate_buffer->prediction_ptr->buffer_cr,
                        cu_chroma_origin_index,
                        candidate_buffer->prediction_ptr->stride_cr,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cr,
                        cu_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cr,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->bwidth_uv,
                        context_ptr->blk_geom->bheight_uv);

        full_loop_r(context_ptr->sb_ptr,
                    candidate_buffer,
                    context_ptr,
                    input_picture_ptr,
                    pcs_ptr,
                    PICTURE_BUFFER_DESC_CHROMA_MASK,
                    cb_qp,
                    cr_qp,
                    &(*count_non_zero_coeffs[1]),
                    &(*count_non_zero_coeffs[2]));

        cu_full_distortion_fast_txb_mode_r(context_ptr->sb_ptr,
                                           candidate_buffer,
                                           context_ptr,
                                           candidate_buffer->candidate_ptr,
                                           pcs_ptr,
                                           input_picture_ptr,
                                           cb_full_distortion,
                                           cr_full_distortion,
                                           count_non_zero_coeffs,
                                           COMPONENT_CHROMA,
                                           &cb_coeff_bits,
                                           &cr_coeff_bits,
                                           1);

        coeff_rate[candidate_buffer->candidate_ptr->intra_chroma_mode]
                  [MAX_ANGLE_DELTA + candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_UV]] =
                      (int)(cb_coeff_bits + cr_coeff_bits);
        distortion[candidate_buffer->candidate_ptr->intra_chroma_mode]
                  [MAX_ANGLE_DELTA + candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_UV]] =
                      (int)(cb_full_distortion[DIST_CALC_RESIDUAL] +
                            cr_full_distortion[DIST_CALC_RESIDUAL]);
    }

    // Loop over all intra mode, then over all uv move to derive the best uv mode for a given intra mode in term of rate

    uint8_t intra_mode_end = context_ptr->md_enable_paeth
        ? PAETH_PRED
        : context_ptr->md_enable_smooth ? SMOOTH_H_PRED : D67_PRED;
    // intra_mode loop (luma mode loop)
    for (uint8_t intra_mode = DC_PRED; intra_mode <= intra_mode_end; ++intra_mode) {
        uint8_t angle_delta_candidate_count = (use_angle_delta &&
                                               av1_is_directional_mode((PredictionMode)intra_mode))
            ? 7
            : 1;
        uint8_t angle_delta_shift = 1;

        for (uint8_t angle_delta_counter = 0; angle_delta_counter < angle_delta_candidate_count;
             ++angle_delta_counter) {
            int32_t angle_delta = CLIP(angle_delta_shift *
                                           (angle_delta_candidate_count == 1 ? 0
                                                                             : angle_delta_counter -
                                                    (angle_delta_candidate_count >> 1)),
                                       -MAX_ANGLE_DELTA,
                                       MAX_ANGLE_DELTA);

            // uv mode loop
            context_ptr->best_uv_cost[intra_mode][MAX_ANGLE_DELTA + angle_delta] = (uint64_t)~0;

            for (uint8_t uv_mode_count = 0;
                 uv_mode_count < MIN(uv_mode_total_count, uv_mode_nfl_count);
                 uv_mode_count++) {
                ModeDecisionCandidate *candidate_ptr = &(
                    context_ptr
                        ->fast_candidate_array[uv_cand_buff_indices[uv_mode_count] -
                                               start_full_buffer_index + start_fast_buffer_index]);

                candidate_ptr->intra_luma_mode          = intra_mode;
                candidate_ptr->is_directional_mode_flag = (uint8_t)av1_is_directional_mode(
                    (PredictionMode)intra_mode);
                candidate_ptr->angle_delta[PLANE_TYPE_Y] = angle_delta;
                candidate_ptr->pred_mode                 = (PredictionMode)intra_mode;

                // Fast Cost
                av1_product_fast_cost_func_table[candidate_ptr->type](
                    context_ptr->blk_ptr,
                    candidate_ptr,
                    context_ptr->qp,
                    0,
                    0,
                    0,
                    0,
                    pcs_ptr,
                    &(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                          .ed_ref_mv_stack[candidate_ptr->ref_frame_type][0]),
                    context_ptr->blk_geom,
                    context_ptr->blk_origin_y >> MI_SIZE_LOG2,
                    context_ptr->blk_origin_x >> MI_SIZE_LOG2,
                    context_ptr->md_enable_inter_intra,
                    context_ptr->full_cost_shut_fast_rate_flag,
                    1,
                    context_ptr->intra_luma_left_mode,
                    context_ptr->intra_luma_top_mode);

                uint64_t rate =
                    coeff_rate[candidate_ptr->intra_chroma_mode]
                              [MAX_ANGLE_DELTA + candidate_ptr->angle_delta[PLANE_TYPE_UV]] +
                    candidate_ptr->fast_luma_rate + candidate_ptr->fast_chroma_rate;
                uint64_t uv_cost = RDCOST(
                    full_lambda,
                    rate,
                    distortion[candidate_ptr->intra_chroma_mode]
                              [MAX_ANGLE_DELTA + candidate_ptr->angle_delta[PLANE_TYPE_UV]]);

                if (uv_cost <
                    context_ptr->best_uv_cost[intra_mode][MAX_ANGLE_DELTA + angle_delta]) {
                    context_ptr->best_uv_mode[intra_mode][MAX_ANGLE_DELTA + angle_delta] =
                        candidate_ptr->intra_chroma_mode;
                    context_ptr->best_uv_angle[intra_mode][MAX_ANGLE_DELTA + angle_delta] =
                        candidate_ptr->angle_delta[PLANE_TYPE_UV];

                    context_ptr->best_uv_cost[intra_mode][MAX_ANGLE_DELTA + angle_delta] = uv_cost;
                    context_ptr->fast_luma_rate[intra_mode][MAX_ANGLE_DELTA + angle_delta] =
                        candidate_ptr->fast_luma_rate;
                    context_ptr->fast_chroma_rate[intra_mode][MAX_ANGLE_DELTA + angle_delta] =
                        candidate_ptr->fast_chroma_rate;
                }
            }
        }
    }

    if (context_ptr->chroma_at_last_md_stage) {
        context_ptr->md_staging_skip_rdoq = tem_md_staging_skip_rdoq;
    }
    // End uv search path
    context_ptr->uv_search_path = EB_FALSE;
}
extern AomVarianceFnPtr mefn_ptr[BlockSizeS_ALL];
unsigned int eb_av1_get_sby_perpixel_variance(const AomVarianceFnPtr *fn_ptr, const uint8_t *src,
                                              int stride, BlockSize bs);

void interintra_class_pruning_1(ModeDecisionContext *context_ptr, uint64_t best_md_stage_cost) {
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL;
         cand_class_it++) {
        if (context_ptr->md_stage_1_cand_prune_th != (uint64_t)~0 ||
            context_ptr->md_stage_1_class_prune_th != (uint64_t)~0)
            if (context_ptr->md_stage_0_count[cand_class_it] > 0 &&
                context_ptr->md_stage_1_count[cand_class_it] > 0) {
                uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
                uint64_t  class_best_cost   = *(
                    context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->fast_cost_ptr);

                // inter class pruning
                if (best_md_stage_cost && class_best_cost &&
                    ((((class_best_cost - best_md_stage_cost) * 100) / best_md_stage_cost) >
                     context_ptr->md_stage_1_class_prune_th)) {
                    context_ptr->md_stage_1_count[cand_class_it] = 0;
                    continue;
                }
                // intra class pruning
                uint32_t cand_count = 1;
                if (class_best_cost)
                    while (
                        cand_count < context_ptr->md_stage_1_count[cand_class_it] &&
                        ((((*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[cand_count]]
                                  ->fast_cost_ptr) -
                            class_best_cost) *
                           100) /
                          class_best_cost) < context_ptr->md_stage_1_cand_prune_th)) {
                        cand_count++;
                    }
                context_ptr->md_stage_1_count[cand_class_it] = cand_count;
            }
        context_ptr->md_stage_1_total_count += context_ptr->md_stage_1_count[cand_class_it];
    }
}

void interintra_class_pruning_2(ModeDecisionContext *context_ptr, uint64_t best_md_stage_cost) {
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL;
         cand_class_it++) {
        if (context_ptr->md_stage_2_3_cand_prune_th != (uint64_t)~0 ||
            context_ptr->md_stage_2_3_class_prune_th != (uint64_t)~0)
            if (context_ptr->md_stage_1_count[cand_class_it] > 0 &&
                context_ptr->md_stage_2_count[cand_class_it] > 0 &&
                context_ptr->bypass_md_stage_1[cand_class_it] == EB_FALSE) {
                uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
                uint64_t  class_best_cost   = *(
                    context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->full_cost_ptr);

                // inter class pruning
                if (best_md_stage_cost && class_best_cost &&
                    ((((class_best_cost - best_md_stage_cost) * 100) / best_md_stage_cost) >
                     context_ptr->md_stage_2_3_class_prune_th)) {
                    context_ptr->md_stage_2_count[cand_class_it] = 0;
                    continue;
                }

                // intra class pruning
                uint32_t cand_count                 = 1;
                uint64_t md_stage_2_3_cand_prune_th = context_ptr->md_stage_2_3_cand_prune_th;
                md_stage_2_3_cand_prune_th          = (cand_class_it == CAND_CLASS_0 ||
                                              cand_class_it == CAND_CLASS_6 ||
                                              cand_class_it == CAND_CLASS_7)
                    ? (uint64_t)~0
                    : (context_ptr->blk_geom->shape == PART_N) ? (uint64_t)~0
                                                               : md_stage_2_3_cand_prune_th;
                if (class_best_cost)
                    while (
                        cand_count < context_ptr->md_stage_2_count[cand_class_it] &&
                        ((((*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[cand_count]]
                                  ->full_cost_ptr) -
                            class_best_cost) *
                           100) /
                          class_best_cost) < md_stage_2_3_cand_prune_th)) {
                        cand_count++;
                    }
                context_ptr->md_stage_2_count[cand_class_it] = cand_count;
            }
        context_ptr->md_stage_2_total_count += context_ptr->md_stage_2_count[cand_class_it];
    }
}

void interintra_class_pruning_3(ModeDecisionContext *context_ptr, uint64_t best_md_stage_cost) {
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL;
         cand_class_it++) {
        if (context_ptr->md_stage_2_3_cand_prune_th != (uint64_t)~0 ||
            context_ptr->md_stage_2_3_class_prune_th != (uint64_t)~0)
            if (context_ptr->md_stage_2_count[cand_class_it] > 0 &&
                context_ptr->md_stage_3_count[cand_class_it] > 0 &&
                context_ptr->bypass_md_stage_2[cand_class_it] == EB_FALSE) {
                uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
                uint64_t  class_best_cost   = *(
                    context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->full_cost_ptr);

                // inter class pruning
                if (best_md_stage_cost && class_best_cost &&
                    ((((class_best_cost - best_md_stage_cost) * 100) / best_md_stage_cost) >
                     context_ptr->md_stage_2_3_class_prune_th)) {
                    context_ptr->md_stage_3_count[cand_class_it] = 0;
                    continue;
                }

                // intra class pruning
                uint32_t cand_count = 1;
                if (class_best_cost)
                    while (
                        cand_count < context_ptr->md_stage_3_count[cand_class_it] &&
                        ((((*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[cand_count]]
                                  ->full_cost_ptr) -
                            class_best_cost) *
                           100) /
                          class_best_cost) < context_ptr->md_stage_2_3_cand_prune_th)) {
                        cand_count++;
                    }
                context_ptr->md_stage_3_count[cand_class_it] = cand_count;
            }
        context_ptr->md_stage_3_total_count += context_ptr->md_stage_3_count[cand_class_it];
    }
}

void md_encode_block(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                     EbPictureBufferDesc *        input_picture_ptr,
                     ModeDecisionCandidateBuffer *bestcandidate_buffers[5]) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array;
    const BlockGeom *             blk_geom = context_ptr->blk_geom;
    ModeDecisionCandidateBuffer * candidate_buffer;
    ModeDecisionCandidate *       fast_candidate_array = context_ptr->fast_candidate_array;
    uint32_t                      candidate_index;
    uint32_t                      fast_candidate_total_count;
    uint32_t                      best_intra_mode = EB_INTRA_MODE_INVALID;
    const uint32_t input_origin_index = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);

    const uint32_t input_cb_origin_in_index = ((context_ptr->round_origin_y >> 1) +
                                               (input_picture_ptr->origin_y >> 1)) *
            input_picture_ptr->stride_cb +
        ((context_ptr->round_origin_x >> 1) + (input_picture_ptr->origin_x >> 1));
    const uint32_t blk_origin_index        = blk_geom->origin_x + blk_geom->origin_y * SB_STRIDE_Y;
    const uint32_t blk_chroma_origin_index = ROUND_UV(blk_geom->origin_x) / 2 +
        ROUND_UV(blk_geom->origin_y) / 2 * SB_STRIDE_UV;
    BlkStruct *blk_ptr         = context_ptr->blk_ptr;
    candidate_buffer_ptr_array = &(candidate_buffer_ptr_array_base[0]);

    signal_derivation_block(pcs_ptr, context_ptr);

    for (uint8_t ref_idx = 0; ref_idx < MAX_REF_TYPE_CAND; ref_idx++)
        context_ptr->ref_best_cost_sq_table[ref_idx] = MAX_CU_COST;

    const AomVarianceFnPtr *fn_ptr = &mefn_ptr[context_ptr->blk_geom->bsize];
    context_ptr->source_variance   = eb_av1_get_sby_perpixel_variance(
        fn_ptr,
        (input_picture_ptr->buffer_y + input_origin_index),
        input_picture_ptr->stride_y,
        context_ptr->blk_geom->bsize);
    blk_ptr->av1xd->tile.mi_col_start = context_ptr->sb_ptr->tile_info.mi_col_start;
    blk_ptr->av1xd->tile.mi_col_end   = context_ptr->sb_ptr->tile_info.mi_col_end;
    blk_ptr->av1xd->tile.mi_row_start = context_ptr->sb_ptr->tile_info.mi_row_start;
    blk_ptr->av1xd->tile.mi_row_end   = context_ptr->sb_ptr->tile_info.mi_row_end;

    product_coding_loop_init_fast_loop(context_ptr,
                                       context_ptr->skip_coeff_neighbor_array,
                                       context_ptr->inter_pred_dir_neighbor_array,
                                       context_ptr->ref_frame_type_neighbor_array,
                                       context_ptr->intra_luma_mode_neighbor_array,
                                       context_ptr->skip_flag_neighbor_array,
                                       context_ptr->mode_type_neighbor_array,
                                       context_ptr->leaf_depth_neighbor_array,
                                       context_ptr->leaf_partition_neighbor_array);

    // Initialize uv_search_path
    context_ptr->uv_search_path = EB_FALSE;
    if (context_ptr->chroma_at_last_md_stage) {
        if (context_ptr->blk_geom->sq_size < 128) {
            if (context_ptr->blk_geom->has_uv) {
                init_chroma_mode(context_ptr);
            }
        }
    } else {
        // Search the best independent intra chroma mode
        if (context_ptr->chroma_level == CHROMA_MODE_0) {
            if (context_ptr->blk_geom->sq_size < 128) {
                if (context_ptr->blk_geom->has_uv) {
                    search_best_independent_uv_mode(pcs_ptr,
                                                    input_picture_ptr,
                                                    input_cb_origin_in_index,
                                                    input_cb_origin_in_index,
                                                    blk_chroma_origin_index,
                                                    context_ptr);
                }
            }
        }
    }
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    // Generate MVP(s)
    if (!context_ptr->md_skip_mvp_generation) {
        if (frm_hdr->allow_intrabc) // pcs_ptr->slice_type == I_SLICE
            generate_av1_mvp_table(&context_ptr->sb_ptr->tile_info,
                                   context_ptr,
                                   context_ptr->blk_ptr,
                                   context_ptr->blk_geom,
                                   context_ptr->blk_origin_x,
                                   context_ptr->blk_origin_y,
                                   pcs_ptr->parent_pcs_ptr->ref_frame_type_arr,
                                   1,
                                   pcs_ptr);
        else if (pcs_ptr->slice_type != I_SLICE)
            generate_av1_mvp_table(&context_ptr->sb_ptr->tile_info,
                                   context_ptr,
                                   context_ptr->blk_ptr,
                                   context_ptr->blk_geom,
                                   context_ptr->blk_origin_x,
                                   context_ptr->blk_origin_y,
                                   pcs_ptr->parent_pcs_ptr->ref_frame_type_arr,
                                   pcs_ptr->parent_pcs_ptr->tot_ref_frame_types,
                                   pcs_ptr);
    } else {
        mvp_bypass_init(pcs_ptr, context_ptr);
    }

    // Read and (if needed) perform 1/8 Pel ME MVs refinement
    read_refine_me_mvs(
        pcs_ptr, context_ptr, input_picture_ptr, input_origin_index, blk_origin_index);
    // Perform ME search around the best MVP
    if (context_ptr->predictive_me_level)
        predictive_me_search(
            pcs_ptr, context_ptr, input_picture_ptr, input_origin_index, blk_origin_index);
    //for every CU, perform Luma DC/V/H/S intra prediction to be used later in inter-intra search

    int allow_ii = is_interintra_allowed_bsize(context_ptr->blk_geom->bsize);
    if (context_ptr->md_enable_inter_intra && allow_ii)
        precompute_intra_pred_for_inter_intra(pcs_ptr, context_ptr);

    generate_md_stage_0_cand(
        context_ptr->sb_ptr, context_ptr, &fast_candidate_total_count, pcs_ptr);

    //MD Stages
    //The first stage(old fast loop) and the last stage(old full loop) should remain at their locations, new stages could be created between those two.
    //a bypass mechanism should be added to skip one or all of the intermediate stages, in a way to to be able to fall back to org design (FastLoop->FullLoop)
    set_md_stage_counts(pcs_ptr, context_ptr, fast_candidate_total_count);

    CandClass cand_class_it;
    uint32_t  buffer_start_idx = 0;
    uint32_t  buffer_count_for_curr_class;
    uint32_t  buffer_total_count        = 0;
    context_ptr->md_stage_1_total_count = 0;
    context_ptr->md_stage_2_total_count = 0;
    context_ptr->md_stage_3_total_count = 0;
    uint64_t best_md_stage_cost         = (uint64_t)~0;
    context_ptr->md_stage               = MD_STAGE_0;

    for (cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++) {
        //number of next level candidates could not exceed number of curr level candidates
        context_ptr->md_stage_1_count[cand_class_it] = MIN(
            context_ptr->md_stage_0_count[cand_class_it],
            context_ptr->md_stage_1_count[cand_class_it]);

        if (context_ptr->md_stage_0_count[cand_class_it] > 0 &&
            context_ptr->md_stage_1_count[cand_class_it] > 0) {
            buffer_count_for_curr_class = context_ptr->md_stage_0_count[cand_class_it] >
                    context_ptr->md_stage_1_count[cand_class_it]
                ? (context_ptr->md_stage_1_count[cand_class_it] + 1)
                : context_ptr->md_stage_1_count[cand_class_it];

            buffer_total_count += buffer_count_for_curr_class;
            assert(buffer_total_count <= MAX_NFL_BUFF && "not enough cand buffers");

            //Input: md_stage_0_count[cand_class_it]  Output:  md_stage_1_count[cand_class_it]
            context_ptr->target_class = cand_class_it;

            md_stage_0(
                pcs_ptr,
                context_ptr,
                candidate_buffer_ptr_array_base,
                fast_candidate_array,
                0,
                fast_candidate_total_count - 1,
                input_picture_ptr,
                input_origin_index,
                input_cb_origin_in_index,
                input_cb_origin_in_index,
                blk_ptr,
                blk_origin_index,
                blk_chroma_origin_index,
                buffer_start_idx,
                buffer_count_for_curr_class,
                context_ptr->md_stage_0_count[cand_class_it] >
                    context_ptr
                        ->md_stage_1_count[cand_class_it]); //is there need to max the temp buffer

            //Sort:  md_stage_1_count[cand_class_it]
            memset(context_ptr->cand_buff_indices[cand_class_it],
                   0xFFFFFFFF,
                   MAX_NFL_BUFF * sizeof(uint32_t));
            sort_fast_cost_based_candidates(
                context_ptr,
                buffer_start_idx,
                buffer_count_for_curr_class, //how many cand buffers to sort. one of the buffers can have max cost.
                context_ptr->cand_buff_indices[cand_class_it]);
            uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
            best_md_stage_cost          = MIN(
                (*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->fast_cost_ptr)),
                best_md_stage_cost);

            buffer_start_idx += buffer_count_for_curr_class; //for next iteration.
        }
    }
    interintra_class_pruning_1(context_ptr, best_md_stage_cost);
    memset(context_ptr->best_candidate_index_array, 0xFFFFFFFF, MAX_NFL_BUFF * sizeof(uint32_t));
    memset(context_ptr->sorted_candidate_index_array, 0xFFFFFFFF, MAX_NFL * sizeof(uint32_t));

    uint64_t ref_fast_cost = MAX_MODE_COST;
    construct_best_sorted_arrays_md_stage_1(context_ptr,
                                            candidate_buffer_ptr_array,
                                            context_ptr->best_candidate_index_array,
                                            context_ptr->sorted_candidate_index_array,
                                            &ref_fast_cost);

    // 1st Full-Loop
    best_md_stage_cost    = (uint64_t)~0;
    context_ptr->md_stage = MD_STAGE_1;
    for (cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++) {
        //number of next level candidates could not exceed number of curr level candidates
        context_ptr->md_stage_2_count[cand_class_it] = MIN(
            context_ptr->md_stage_1_count[cand_class_it],
            context_ptr->md_stage_2_count[cand_class_it]);
        if (context_ptr->bypass_md_stage_1[cand_class_it] == EB_FALSE &&
            context_ptr->md_stage_1_count[cand_class_it] > 0 &&
            context_ptr->md_stage_2_count[cand_class_it] > 0) {
            context_ptr->target_class = cand_class_it;
            md_stage_1(pcs_ptr,
                       context_ptr->sb_ptr,
                       blk_ptr,
                       context_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index,
                       ref_fast_cost);

            // Sort the candidates of the target class based on the 1st full loop cost

            //sort the new set of candidates
            if (context_ptr->md_stage_1_count[cand_class_it])
                sort_full_cost_based_candidates(context_ptr,
                                                context_ptr->md_stage_1_count[cand_class_it],
                                                context_ptr->cand_buff_indices[cand_class_it]);
            uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
            best_md_stage_cost          = MIN(
                (*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->full_cost_ptr)),
                best_md_stage_cost);
        }
    }
    interintra_class_pruning_2(context_ptr, best_md_stage_cost);

    // 2nd Full-Loop
    best_md_stage_cost    = (uint64_t)~0;
    context_ptr->md_stage = MD_STAGE_2;
    for (cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++) {
        //number of next level candidates could not exceed number of curr level candidates
        context_ptr->md_stage_3_count[cand_class_it] = MIN(
            context_ptr->md_stage_2_count[cand_class_it],
            context_ptr->md_stage_3_count[cand_class_it]);

        if (context_ptr->bypass_md_stage_2[cand_class_it] == EB_FALSE &&
            context_ptr->md_stage_2_count[cand_class_it] > 0 &&
            context_ptr->md_stage_3_count[cand_class_it] > 0) {
            context_ptr->target_class = cand_class_it;

            md_stage_2(pcs_ptr,
                       context_ptr->sb_ptr,
                       blk_ptr,
                       context_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index,
                       ref_fast_cost);

            // Sort the candidates of the target class based on the 1st full loop cost

            //sort the new set of candidates
            if (context_ptr->md_stage_2_count[cand_class_it])
                sort_full_cost_based_candidates(context_ptr,
                                                context_ptr->md_stage_2_count[cand_class_it],
                                                context_ptr->cand_buff_indices[cand_class_it]);

            uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
            best_md_stage_cost          = MIN(
                (*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->full_cost_ptr)),
                best_md_stage_cost);
        }
    }

    interintra_class_pruning_3(context_ptr, best_md_stage_cost);

    assert(context_ptr->md_stage_3_total_count <= MAX_NFL);
    assert(context_ptr->md_stage_3_total_count > 0);
    construct_best_sorted_arrays_md_stage_3(context_ptr,
                                            candidate_buffer_ptr_array,
                                            context_ptr->best_candidate_index_array,
                                            context_ptr->sorted_candidate_index_array);
    // Search the best independent intra chroma mode
    if (context_ptr->chroma_at_last_md_stage) {
        // Initialize uv_search_path
        context_ptr->uv_search_path = EB_FALSE;
        if (context_ptr->blk_geom->sq_size < 128) {
            if (context_ptr->blk_geom->has_uv) {
                if (context_ptr->md_stage_3_total_intra_count)
                    search_best_independent_uv_mode(pcs_ptr,
                                                    input_picture_ptr,
                                                    input_cb_origin_in_index,
                                                    input_cb_origin_in_index,
                                                    blk_chroma_origin_index,
                                                    context_ptr);
            }
        }
    }

    // 3rd Full-Loop
    context_ptr->md_stage = MD_STAGE_3;
    md_stage_3(pcs_ptr,
               context_ptr->sb_ptr,
               blk_ptr,
               context_ptr,
               input_picture_ptr,
               input_origin_index,
               input_cb_origin_in_index,
               blk_origin_index,
               blk_chroma_origin_index,
               context_ptr->md_stage_3_total_count,
               ref_fast_cost); // fullCandidateTotalCount to number of buffers to process

    // Full Mode Decision (choose the best mode)
    candidate_index  = product_full_mode_decision(context_ptr,
                                                 blk_ptr,
                                                 candidate_buffer_ptr_array,
                                                 context_ptr->md_stage_3_total_count,
                                                 (context_ptr->full_loop_escape == 2)
                                                     ? context_ptr->sorted_candidate_index_array
                                                     : context_ptr->best_candidate_index_array,
                                                 context_ptr->prune_ref_frame_for_rec_partitions,
                                                 &best_intra_mode);
    candidate_buffer = candidate_buffer_ptr_array[candidate_index];

    bestcandidate_buffers[0] = candidate_buffer;
    uint8_t sq_index         = eb_log2f(context_ptr->blk_geom->sq_size) - 2;
    if (context_ptr->blk_geom->shape == PART_N) {
        context_ptr->parent_sq_type[sq_index] = candidate_buffer->candidate_ptr->type;

        context_ptr->parent_sq_has_coeff[sq_index] =
            (candidate_buffer->candidate_ptr->y_has_coeff ||
             candidate_buffer->candidate_ptr->u_has_coeff ||
             candidate_buffer->candidate_ptr->v_has_coeff)
            ? 1
            : 0;

        context_ptr->parent_sq_pred_mode[sq_index] = candidate_buffer->candidate_ptr->pred_mode;
    }

    av1_perform_inverse_transform_recon(context_ptr, candidate_buffer, context_ptr->blk_geom);

    if (!context_ptr->blk_geom->has_uv) {
        // Store the luma data for 4x* and *x4 blocks to be used for CFL
        EbPictureBufferDesc *recon_ptr       = candidate_buffer->recon_ptr;
        uint32_t             rec_luma_offset = context_ptr->blk_geom->origin_x +
            context_ptr->blk_geom->origin_y * recon_ptr->stride_y;
        if (context_ptr->hbd_mode_decision) {
            for (uint32_t j = 0; j < context_ptr->blk_geom->bheight; ++j)
                memcpy(
                    context_ptr->cfl_temp_luma_recon16bit + rec_luma_offset +
                        j * recon_ptr->stride_y,
                    ((uint16_t *)recon_ptr->buffer_y) + (rec_luma_offset + j * recon_ptr->stride_y),
                    sizeof(uint16_t) * context_ptr->blk_geom->bwidth);
        } else {
            for (uint32_t j = 0; j < context_ptr->blk_geom->bheight; ++j)
                memcpy(&context_ptr->cfl_temp_luma_recon[rec_luma_offset + j * recon_ptr->stride_y],
                       recon_ptr->buffer_y + rec_luma_offset + j * recon_ptr->stride_y,
                       context_ptr->blk_geom->bwidth);
        }
    }
    //copy neigh recon data in blk_ptr
    {
        uint32_t             j;
        EbPictureBufferDesc *recon_ptr       = candidate_buffer->recon_ptr;
        uint32_t             rec_luma_offset = context_ptr->blk_geom->origin_x +
            context_ptr->blk_geom->origin_y * recon_ptr->stride_y;

        uint32_t rec_cb_offset = ((((context_ptr->blk_geom->origin_x >> 3) << 3) +
                                   ((context_ptr->blk_geom->origin_y >> 3) << 3) *
                                       candidate_buffer->recon_ptr->stride_cb) >>
                                  1);
        uint32_t rec_cr_offset = ((((context_ptr->blk_geom->origin_x >> 3) << 3) +
                                   ((context_ptr->blk_geom->origin_y >> 3) << 3) *
                                       candidate_buffer->recon_ptr->stride_cr) >>
                                  1);

        if (!context_ptr->hbd_mode_decision) {
            memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                       .neigh_top_recon[0],
                   recon_ptr->buffer_y + rec_luma_offset +
                       (context_ptr->blk_geom->bheight - 1) * recon_ptr->stride_y,
                   context_ptr->blk_geom->bwidth);
            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
                memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon[1],
                       recon_ptr->buffer_cb + rec_cb_offset +
                           (context_ptr->blk_geom->bheight_uv - 1) * recon_ptr->stride_cb,
                       context_ptr->blk_geom->bwidth_uv);
                memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon[2],
                       recon_ptr->buffer_cr + rec_cr_offset +
                           (context_ptr->blk_geom->bheight_uv - 1) * recon_ptr->stride_cr,
                       context_ptr->blk_geom->bwidth_uv);
            }

            for (j = 0; j < context_ptr->blk_geom->bheight; ++j)
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[0][j] =
                    recon_ptr->buffer_y[rec_luma_offset + context_ptr->blk_geom->bwidth - 1 +
                                        j * recon_ptr->stride_y];

            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
                for (j = 0; j < context_ptr->blk_geom->bheight_uv; ++j) {
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon[1][j] =
                        recon_ptr->buffer_cb[rec_cb_offset + context_ptr->blk_geom->bwidth_uv - 1 +
                                             j * recon_ptr->stride_cb];
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon[2][j] =
                        recon_ptr->buffer_cr[rec_cr_offset + context_ptr->blk_geom->bwidth_uv - 1 +
                                             j * recon_ptr->stride_cr];
                }
            }
        } else {
            uint16_t sz = sizeof(uint16_t);
            memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                       .neigh_top_recon_16bit[0],
                   recon_ptr->buffer_y +
                       sz *
                           (rec_luma_offset +
                            (context_ptr->blk_geom->bheight - 1) * recon_ptr->stride_y),
                   sz * context_ptr->blk_geom->bwidth);
            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
                memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon_16bit[1],
                       recon_ptr->buffer_cb +
                           sz *
                               (rec_cb_offset +
                                (context_ptr->blk_geom->bheight_uv - 1) * recon_ptr->stride_cb),
                       sz * context_ptr->blk_geom->bwidth_uv);
                memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon_16bit[2],
                       recon_ptr->buffer_cr +
                           sz *
                               (rec_cr_offset +
                                (context_ptr->blk_geom->bheight_uv - 1) * recon_ptr->stride_cr),
                       sz * context_ptr->blk_geom->bwidth_uv);
            }

            for (j = 0; j < context_ptr->blk_geom->bheight; ++j)
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon_16bit[0][j] = ((uint16_t *)recon_ptr->buffer_y)
                    [rec_luma_offset + context_ptr->blk_geom->bwidth - 1 + j * recon_ptr->stride_y];

            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1) {
                for (j = 0; j < context_ptr->blk_geom->bheight_uv; ++j) {
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon_16bit[1][j] =
                        ((uint16_t *)recon_ptr
                             ->buffer_cb)[rec_cb_offset + context_ptr->blk_geom->bwidth_uv - 1 +
                                          j * recon_ptr->stride_cb];
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon_16bit[2][j] =
                        ((uint16_t *)recon_ptr
                             ->buffer_cr)[rec_cr_offset + context_ptr->blk_geom->bwidth_uv - 1 +
                                          j * recon_ptr->stride_cr];
                }
            }
        }
    }

#if NO_ENCDEC
    //copy recon
    uint32_t txb_origin_index = context_ptr->blk_geom->origin_x +
        (context_ptr->blk_geom->origin_y * 128);
    uint32_t bwidth  = context_ptr->blk_geom->bwidth;
    uint32_t bheight = context_ptr->blk_geom->bheight;

    if (!context_ptr->hbd_mode_decision) {
        uint8_t *src_ptr = &(((uint8_t *)candidate_buffer->recon_ptr->buffer_y)[txb_origin_index]);
        uint8_t *dst_ptr = &(((uint8_t *)context_ptr->blk_ptr->recon_tmp->buffer_y)[0]);

        uint32_t j;
        for (j = 0; j < bheight; j++)
            memcpy(dst_ptr + j * 128, src_ptr + j * 128, bwidth * sizeof(uint8_t));

        if (context_ptr->blk_geom->has_uv) {
            uint32_t txb_origin_index = ((((context_ptr->blk_geom->origin_x >> 3) << 3) +
                                          ((context_ptr->blk_geom->origin_y >> 3) << 3) *
                                              candidate_buffer->recon_ptr->stride_cb) >>
                                         1);
            bwidth                    = context_ptr->blk_geom->bwidth_uv;
            bheight                   = context_ptr->blk_geom->bheight_uv;

            // Cb
            src_ptr = &(((uint8_t *)candidate_buffer->recon_ptr->buffer_cb)[txb_origin_index]);
            dst_ptr = &(((uint8_t *)context_ptr->blk_ptr->recon_tmp->buffer_cb)[0]);

            for (j = 0; j < bheight; j++)
                memcpy(dst_ptr + j * 64, src_ptr + j * 64, bwidth * sizeof(uint8_t));

            // Cr
            src_ptr = &(((uint8_t *)candidate_buffer->recon_ptr->buffer_cr)[txb_origin_index]);
            dst_ptr = &(((uint8_t *)context_ptr->blk_ptr->recon_tmp->buffer_cr)[0]);

            for (j = 0; j < bheight; j++)
                memcpy(dst_ptr + j * 64, src_ptr + j * 64, bwidth * sizeof(uint8_t));
        }
    } else {
        uint16_t *src_ptr = ((uint16_t *)candidate_buffer->recon_ptr->buffer_y) + txb_origin_index;
        uint16_t *dst_ptr = (uint16_t *)context_ptr->blk_ptr->recon_tmp->buffer_y;
        for (uint32_t j = 0; j < bheight; j++)
            memcpy(dst_ptr + j * 128, src_ptr + j * 128, bwidth * sizeof(uint16_t));

        if (context_ptr->blk_geom->has_uv) {
            txb_origin_index = ((((context_ptr->blk_geom->origin_x >> 3) << 3) +
                                 ((context_ptr->blk_geom->origin_y >> 3) << 3) *
                                     candidate_buffer->recon_ptr->stride_cb) >>
                                1);
            bwidth           = context_ptr->blk_geom->bwidth_uv;
            bheight          = context_ptr->blk_geom->bheight_uv;

            // Cb
            src_ptr = ((uint16_t *)candidate_buffer->recon_ptr->buffer_cb) + txb_origin_index;
            dst_ptr = (uint16_t *)context_ptr->blk_ptr->recon_tmp->buffer_cb;
            for (uint32_t j = 0; j < bheight; j++)
                memcpy(dst_ptr + j * 64, src_ptr + j * 64, bwidth * sizeof(uint16_t));

            // Cr
            src_ptr = ((uint16_t *)candidate_buffer->recon_ptr->buffer_cr) + txb_origin_index;
            dst_ptr = (uint16_t *)context_ptr->blk_ptr->recon_tmp->buffer_cr;
            for (uint32_t j = 0; j < bheight; j++)
                memcpy(dst_ptr + j * 64, src_ptr + j * 64, bwidth * sizeof(uint16_t));
        }
    }
#endif

    context_ptr->md_local_blk_unit[blk_ptr->mds_idx].avail_blk_flag = EB_TRUE;
}

/*
 * Determine if the evaluation of nsq blocks (HA, HB, VA, VB, H4, V4) can be skipped
 * based on the relative cost of the SQ, H, and V blocks.  The scaling factor sq_weight
 * determines how likely it is to skip blocks, and is a function of the qp, block shape,
 * prediction mode, block coeffs, and encode mode.
 *
 * skip HA, HB and H4 if (valid SQ and H) and (H_COST > (SQ_WEIGHT * SQ_COST) / 100)
 * skip VA, VB and V4 if (valid SQ and V) and (V_COST > (SQ_WEIGHT * SQ_COST) / 100)
 *
 * Returns TRUE if the blocks should be skipped; FALSE otherwise.
 */
uint8_t update_skip_nsq_shapes(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                               ModeDecisionContext *context_ptr) {
    uint8_t  skip_nsq  = 0;
    uint32_t sq_weight = context_ptr->sq_weight;

    // return immediately if the skip nsq threshold is infinite
    if (sq_weight == (uint32_t)~0)
        return skip_nsq;

    // use an aggressive threshold for QP 20
    if (scs_ptr->static_config.qp <= QP_20)
        sq_weight += AGGRESSIVE_OFFSET_1;

    // use a conservative threshold for H4, V4 blocks
    if (context_ptr->blk_geom->shape == PART_H4 || context_ptr->blk_geom->shape == PART_V4)
        sq_weight += CONSERVATIVE_OFFSET_0;

    uint32_t     sqi           = context_ptr->blk_geom->sqi_mds;
    MdBlkStruct *local_cu_unit = context_ptr->md_local_blk_unit;

    if (context_ptr->blk_geom->shape == PART_HA || context_ptr->blk_geom->shape == PART_HB ||
        context_ptr->blk_geom->shape == PART_H4) {
        if (context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds].avail_blk_flag &&
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds + 1].avail_blk_flag &&
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds + 2].avail_blk_flag) {
            // Use aggressive thresholds for inter blocks
            if (pcs_ptr->slice_type != I_SLICE) {
                if (context_ptr->blk_geom->shape == PART_HA) {
                    if (context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds + 1]
                            .prediction_mode_flag == INTRA_MODE)
                        sq_weight += CONSERVATIVE_OFFSET_0;
                }
                if (context_ptr->blk_geom->shape == PART_HB) {
                    if (context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds + 2]
                            .prediction_mode_flag == INTRA_MODE)
                        sq_weight += CONSERVATIVE_OFFSET_0;
                }
            }

            // Use aggressive thresholds for blocks without coeffs
            if (context_ptr->blk_geom->shape == PART_HA) {
                if (!context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds + 1]
                         .block_has_coeff)
                    sq_weight += AGGRESSIVE_OFFSET_1;
            }
            if (context_ptr->blk_geom->shape == PART_HB) {
                if (!context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds + 2]
                         .block_has_coeff)
                    sq_weight += AGGRESSIVE_OFFSET_1;
            }

            // compute the cost of the SQ block and H block
            uint64_t sq_cost =
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds].default_cost;
            uint64_t h_cost =
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds + 1].default_cost +
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds + 2].default_cost;

            // Determine if nsq shapes can be skipped based on the relative cost of SQ and H blocks
            skip_nsq = (h_cost > ((sq_cost * sq_weight) / 100));

            if (!skip_nsq && context_ptr->nsq_hv_level > 0) {
                if (local_cu_unit[sqi + 3].avail_blk_flag &&
                    local_cu_unit[sqi + 4].avail_blk_flag) {
                    //compute the cost of V partition
                    uint64_t v_cost = local_cu_unit[sqi + 3].default_cost +
                        local_cu_unit[sqi + 4].default_cost;
                    uint32_t offset = 10;
                    if (context_ptr->nsq_hv_level == 2 && context_ptr->blk_geom->shape == PART_H4)
                        offset = 5;
                    if (offset >= 5 && scs_ptr->static_config.qp <= 20)
                        offset -= 5;
                    uint32_t v_weight = 100 + offset;

                    //if the cost of H partition is bigger than the V partition by a certain percentage, skip HA/HB
                    //use 10% to skip HA/HB, use 5% to skip H4, also for very low QP be more aggressive to skip
                    skip_nsq = (h_cost > ((v_cost * v_weight) / 100));
                }
            }
        }
    }
    if (context_ptr->blk_geom->shape == PART_VA || context_ptr->blk_geom->shape == PART_VB ||
        context_ptr->blk_geom->shape == PART_V4) {
        if (context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds].avail_blk_flag &&
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds + 3].avail_blk_flag &&
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds + 4].avail_blk_flag) {
            // Use aggressive thresholds for inter blocks
            if (pcs_ptr->slice_type != I_SLICE) {
                if (context_ptr->blk_geom->shape == PART_VA) {
                    if (context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds + 3]
                            .prediction_mode_flag == INTRA_MODE)
                        sq_weight += CONSERVATIVE_OFFSET_0;
                }
                if (context_ptr->blk_geom->shape == PART_VB) {
                    if (context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds + 4]
                            .prediction_mode_flag == INTRA_MODE)
                        sq_weight += CONSERVATIVE_OFFSET_0;
                }
            }

            // Use aggressive thresholds for blocks without coeffs
            if (context_ptr->blk_geom->shape == PART_VA) {
                if (!context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds + 3]
                         .block_has_coeff)
                    sq_weight += AGGRESSIVE_OFFSET_1;
            }
            if (context_ptr->blk_geom->shape == PART_VB) {
                if (!context_ptr->md_blk_arr_nsq[context_ptr->blk_geom->sqi_mds + 4]
                         .block_has_coeff)
                    sq_weight += AGGRESSIVE_OFFSET_1;
            }

            // compute the cost of the SQ block and V block
            uint64_t sq_cost =
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds].default_cost;
            uint64_t v_cost =
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds + 3].default_cost +
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds + 4].default_cost;

            // Determine if nsq shapes can be skipped based on the relative cost of SQ and V blocks
            skip_nsq = (v_cost > ((sq_cost * sq_weight) / 100));
            if (!skip_nsq && context_ptr->nsq_hv_level > 0) {
                if (local_cu_unit[sqi + 1].avail_blk_flag &&
                    local_cu_unit[sqi + 2].avail_blk_flag) {
                    uint64_t h_cost = local_cu_unit[sqi + 1].default_cost +
                        local_cu_unit[sqi + 2].default_cost;
                    uint32_t offset = 10;

                    if (context_ptr->nsq_hv_level == 2 && context_ptr->blk_geom->shape == PART_V4)
                        offset = 5;
                    if (offset >= 5 && scs_ptr->static_config.qp <= 20)
                        offset -= 5;

                    uint32_t h_weight = 100 + offset;
                    skip_nsq          = (v_cost > ((h_cost * h_weight) / 100));
                }
            }
        }
    }

    return skip_nsq;
}
/***********************************
get the number of total block in a
branch
***********************************/
uint32_t get_number_of_blocks(uint32_t block_idx) {
    const BlockGeom *blk_geom      = get_blk_geom_mds(block_idx);
    uint32_t         tot_d1_blocks = blk_geom->sq_size == 128
        ? 17
        : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;
    return tot_d1_blocks;
}
/***********************************
Mark the blocks of the lower depth
to be skipped
***********************************/
static void set_child_to_be_skipped(ModeDecisionContext *context_ptr, uint32_t blk_index,
                                    int32_t sb_size, int8_t depth_step) {
    uint32_t         child_block_idx_1, child_block_idx_2, child_block_idx_3, child_block_idx_4;
    uint32_t         block_1d_idx;
    const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

    if (context_ptr->md_blk_arr_nsq[blk_index].split_flag && blk_geom->sq_size > 4) {
        //Set first child to be considered
        child_block_idx_1 = blk_index + d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
        uint32_t child1_tot_d1_blocks = get_number_of_blocks(child_block_idx_1);
        for (block_1d_idx = 0; block_1d_idx < child1_tot_d1_blocks; block_1d_idx++)
            context_ptr->md_blk_arr_nsq[child_block_idx_1 + block_1d_idx].do_not_process_block = 1;
        if (depth_step > 1)
            set_child_to_be_skipped(context_ptr, child_block_idx_1, sb_size, depth_step - 1);
        //Set second child to be considered
        child_block_idx_2 = child_block_idx_1 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        uint32_t child2_tot_d1_blocks = get_number_of_blocks(child_block_idx_2);
        for (block_1d_idx = 0; block_1d_idx < child2_tot_d1_blocks; block_1d_idx++)
            context_ptr->md_blk_arr_nsq[child_block_idx_2 + block_1d_idx].do_not_process_block = 1;
        if (depth_step > 1)
            set_child_to_be_skipped(context_ptr, child_block_idx_2, sb_size, depth_step - 1);
        //Set third child to be considered
        child_block_idx_3 = child_block_idx_2 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        uint32_t child3_tot_d1_blocks = get_number_of_blocks(child_block_idx_3);
        for (block_1d_idx = 0; block_1d_idx < child3_tot_d1_blocks; block_1d_idx++)
            context_ptr->md_blk_arr_nsq[child_block_idx_3 + block_1d_idx].do_not_process_block = 1;
        if (depth_step > 1)
            set_child_to_be_skipped(context_ptr, child_block_idx_3, sb_size, depth_step - 1);
        //Set forth child to be considered
        child_block_idx_4 = child_block_idx_3 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        uint32_t child4_tot_d1_blocks = get_number_of_blocks(child_block_idx_4);
        for (block_1d_idx = 0; block_1d_idx < child4_tot_d1_blocks; block_1d_idx++)
            context_ptr->md_blk_arr_nsq[child_block_idx_4 + block_1d_idx].do_not_process_block = 1;
        if (depth_step > 1)
            set_child_to_be_skipped(context_ptr, child_block_idx_4, sb_size, depth_step - 1);
    }
}
EB_EXTERN EbErrorType mode_decision_sb(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                       const MdcSbData *const mdcResultTbPtr, SuperBlock *sb_ptr,
                                       uint16_t sb_origin_x, uint16_t sb_origin_y, uint32_t sb_addr,
                                       ModeDecisionContext *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    //printf("sb_origin_x = %d, sb_origin_y = %d\n", sb_origin_x, sb_origin_y);

    uint32_t                     blk_index;
    ModeDecisionCandidateBuffer *bestcandidate_buffers[5];
    // Pre Intra Search
    uint32_t                   leaf_count      = mdcResultTbPtr->leaf_count;
    const EbMdcLeafData *const leaf_data_array = mdcResultTbPtr->leaf_data_array;
    const uint16_t             tile_idx        = context_ptr->tile_index;
    context_ptr->sb_ptr                        = sb_ptr;

    EbBool all_blk_init = (pcs_ptr->parent_pcs_ptr->pic_depth_mode <= PIC_SQ_DEPTH_MODE);
    init_sq_nsq_block(scs_ptr, context_ptr);

    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];
    // Mode Decision Neighbor Arrays
    context_ptr->intra_luma_mode_neighbor_array =
        pcs_ptr->md_intra_luma_mode_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->intra_chroma_mode_neighbor_array =
        pcs_ptr->md_intra_chroma_mode_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->mv_neighbor_array =
        pcs_ptr->md_mv_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->skip_flag_neighbor_array =
        pcs_ptr->md_skip_flag_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->mode_type_neighbor_array =
        pcs_ptr->md_mode_type_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->leaf_depth_neighbor_array =
        pcs_ptr->md_leaf_depth_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->leaf_partition_neighbor_array =
        pcs_ptr->mdleaf_partition_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];

    if (!context_ptr->hbd_mode_decision) {
        context_ptr->luma_recon_neighbor_array =
            pcs_ptr->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        context_ptr->cb_recon_neighbor_array =
            pcs_ptr->md_cb_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        context_ptr->cr_recon_neighbor_array =
            pcs_ptr->md_cr_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    } else {
        context_ptr->luma_recon_neighbor_array16bit =
            pcs_ptr->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        context_ptr->cb_recon_neighbor_array16bit =
            pcs_ptr->md_cb_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        context_ptr->cr_recon_neighbor_array16bit =
            pcs_ptr->md_cr_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    }
    context_ptr->skip_coeff_neighbor_array =
        pcs_ptr->md_skip_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->luma_dc_sign_level_coeff_neighbor_array =
        pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->cb_dc_sign_level_coeff_neighbor_array =
        pcs_ptr->md_cb_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->cr_dc_sign_level_coeff_neighbor_array =
        pcs_ptr->md_cr_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->txfm_context_array =
        pcs_ptr->md_txfm_context_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->inter_pred_dir_neighbor_array =
        pcs_ptr->md_inter_pred_dir_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->ref_frame_type_neighbor_array =
        pcs_ptr->md_ref_frame_type_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    context_ptr->interpolation_type_neighbor_array =
        pcs_ptr->md_interpolation_type_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    uint32_t             d1_block_itr      = 0;
    uint32_t             d1_first_block    = 1;
    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    if (context_ptr->hbd_mode_decision) {
        const uint32_t input_luma_offset = ((sb_origin_y + input_picture_ptr->origin_y) *
                                            input_picture_ptr->stride_y) +
            (sb_origin_x + input_picture_ptr->origin_x);
        const uint32_t input_bit_inc_luma_offset = ((sb_origin_y + input_picture_ptr->origin_y) *
                                                    input_picture_ptr->stride_bit_inc_y) +
            (sb_origin_x + input_picture_ptr->origin_x);
        const uint32_t input_cb_offset = (((sb_origin_y + input_picture_ptr->origin_y) >> 1) *
                                          input_picture_ptr->stride_cb) +
            ((sb_origin_x + input_picture_ptr->origin_x) >> 1);
        const uint32_t input_bit_inc_cb_offset = (((sb_origin_y + input_picture_ptr->origin_y) >>
                                                   1) *
                                                  input_picture_ptr->stride_bit_inc_cb) +
            ((sb_origin_x + input_picture_ptr->origin_x) >> 1);
        const uint32_t input_cr_offset = (((sb_origin_y + input_picture_ptr->origin_y) >> 1) *
                                          input_picture_ptr->stride_cr) +
            ((sb_origin_x + input_picture_ptr->origin_x) >> 1);
        const uint32_t input_bit_inc_cr_offset = (((sb_origin_y + input_picture_ptr->origin_y) >>
                                                   1) *
                                                  input_picture_ptr->stride_bit_inc_cr) +
            ((sb_origin_x + input_picture_ptr->origin_x) >> 1);

        uint32_t sb_width  = MIN(scs_ptr->sb_size_pix,
                                pcs_ptr->parent_pcs_ptr->aligned_width - sb_origin_x);
        uint32_t sb_height = MIN(scs_ptr->sb_size_pix,
                                 pcs_ptr->parent_pcs_ptr->aligned_height - sb_origin_y);

        pack2d_src(input_picture_ptr->buffer_y + input_luma_offset,
                   input_picture_ptr->stride_y,
                   input_picture_ptr->buffer_bit_inc_y + input_bit_inc_luma_offset,
                   input_picture_ptr->stride_bit_inc_y,
                   (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_y,
                   context_ptr->input_sample16bit_buffer->stride_y,
                   sb_width,
                   sb_height);

        pack2d_src(input_picture_ptr->buffer_cb + input_cb_offset,
                   input_picture_ptr->stride_cb,
                   input_picture_ptr->buffer_bit_inc_cb + input_bit_inc_cb_offset,
                   input_picture_ptr->stride_bit_inc_cb,
                   (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_cb,
                   context_ptr->input_sample16bit_buffer->stride_cb,
                   sb_width >> 1,
                   sb_height >> 1);

        pack2d_src(input_picture_ptr->buffer_cr + input_cr_offset,
                   input_picture_ptr->stride_cr,
                   input_picture_ptr->buffer_bit_inc_cr + input_bit_inc_cr_offset,
                   input_picture_ptr->stride_bit_inc_cr,
                   (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_cr,
                   context_ptr->input_sample16bit_buffer->stride_cr,
                   sb_width >> 1,
                   sb_height >> 1);

        store16bit_input_src(context_ptr->input_sample16bit_buffer,
                             pcs_ptr,
                             sb_origin_x,
                             sb_origin_y,
                             sb_width,
                             sb_height);
        //input_picture_ptr = context_ptr->input_sample16bit_buffer;
        input_picture_ptr = pcs_ptr->input_frame16bit;
    }

    //CU Loop
    blk_index = 0; //index over mdc array

    uint32_t blk_idx_mds                       = 0;
    uint32_t d1_blocks_accumlated              = 0;
    int      skip_next_nsq                     = 0;
    int      skip_next_sq                      = 0;
    uint32_t next_non_skip_blk_idx_mds         = 0;
    int64_t  depth_cost[NUMBER_OF_DEPTH]       = {-1, -1, -1, -1, -1, -1};
    uint64_t nsq_cost[NUMBER_OF_SHAPES]        = {MAX_CU_COST,
                                           MAX_CU_COST,
                                           MAX_CU_COST,
                                           MAX_CU_COST,
                                           MAX_CU_COST,
                                           MAX_CU_COST,
                                           MAX_CU_COST,
                                           MAX_CU_COST,
                                           MAX_CU_COST,
                                           MAX_CU_COST};
    Part     nsq_shape_table[NUMBER_OF_SHAPES] = {
        PART_N, PART_H, PART_V, PART_HA, PART_HB, PART_VA, PART_VB, PART_H4, PART_V4, PART_S};
    uint8_t skip_next_depth = 0;
    do {
        blk_idx_mds = leaf_data_array[blk_index].mds_idx;

        const BlockGeom *blk_geom = context_ptr->blk_geom = get_blk_geom_mds(blk_idx_mds);
        BlkStruct *      blk_ptr = context_ptr->blk_ptr = &context_ptr->md_blk_arr_nsq[blk_idx_mds];

        context_ptr->cu_size_log2 = blk_geom->bwidth_log2;
        context_ptr->blk_origin_x = sb_origin_x + blk_geom->origin_x;
        context_ptr->blk_origin_y = sb_origin_y + blk_geom->origin_y;

        const EbMdcLeafData *const leaf_data_ptr = &mdcResultTbPtr->leaf_data_array[blk_index];
        context_ptr->sb_sz                       = BLOCK_SIZE_64;
        context_ptr->round_origin_x              = ((context_ptr->blk_origin_x >> 3) << 3);
        context_ptr->round_origin_y              = ((context_ptr->blk_origin_y >> 3) << 3);
        context_ptr->sb_origin_x                 = sb_origin_x;
        context_ptr->sb_origin_y                 = sb_origin_y;
        context_ptr->md_local_blk_unit[blk_idx_mds].tested_blk_flag = EB_TRUE;
        context_ptr->md_ep_pipe_sb[blk_idx_mds].merge_cost          = 0;
        context_ptr->md_ep_pipe_sb[blk_idx_mds].skip_cost           = 0;
        blk_ptr->av1xd->sb_type                                     = blk_geom->bsize;
        blk_ptr->mds_idx                                            = blk_idx_mds;
        context_ptr->md_blk_arr_nsq[blk_idx_mds].mdc_split_flag     = (uint16_t)
                                                                      leaf_data_ptr->split_flag;
        context_ptr->md_blk_arr_nsq[blk_geom->sqi_mds].split_flag = (uint16_t)
                                                                        leaf_data_ptr->split_flag;
        blk_ptr->split_flag =
            (uint16_t)leaf_data_ptr
                ->split_flag; //mdc indicates smallest or non valid CUs with split flag=
        blk_ptr->qp                                             = context_ptr->qp;
        context_ptr->md_local_blk_unit[blk_idx_mds].best_d1_blk = blk_idx_mds;
        if (leaf_data_ptr->tot_d1_blocks != 1) {
            // We need to get the index of the sq_block for each NSQ branch
            if (d1_first_block) {
                copy_neighbour_arrays( //save a clean neigh in [1], encode uses [0], reload the clean in [0] after done last ns block in a partition
                    pcs_ptr,
                    context_ptr,
                    0,
                    1,
                    blk_geom->sqi_mds,
                    sb_origin_x,
                    sb_origin_y);
            }
        }
        int32_t       mi_row           = context_ptr->blk_origin_y >> MI_SIZE_LOG2;
        int32_t       mi_col           = context_ptr->blk_origin_x >> MI_SIZE_LOG2;
        int           mi_stride        = pcs_ptr->parent_pcs_ptr->av1_cm->mi_stride;
        const int32_t offset           = mi_row * mi_stride + mi_col;
        blk_ptr->av1xd->mi             = pcs_ptr->mi_grid_base + offset;
        ModeInfo *mi_ptr               = *blk_ptr->av1xd->mi;
        blk_ptr->av1xd->up_available   = (mi_row > sb_ptr->tile_info.mi_row_start);
        blk_ptr->av1xd->left_available = (mi_col > sb_ptr->tile_info.mi_col_start);
        if (blk_ptr->av1xd->up_available)
            blk_ptr->av1xd->above_mbmi = &mi_ptr[-mi_stride].mbmi;
        else
            blk_ptr->av1xd->above_mbmi = NULL;
        if (blk_ptr->av1xd->left_available)
            blk_ptr->av1xd->left_mbmi = &mi_ptr[-1].mbmi;
        else
            blk_ptr->av1xd->left_mbmi = NULL;

        uint8_t  redundant_blk_avail = 0;
        uint16_t redundant_blk_mds;
        if (all_blk_init)
            check_redundant_block(blk_geom, context_ptr, &redundant_blk_avail, &redundant_blk_mds);

        context_ptr->similar_blk_avail = 0;
        if (all_blk_init)
            check_similar_block(blk_geom,
                                context_ptr,
                                &context_ptr->similar_blk_avail,
                                &context_ptr->similar_blk_mds);
        if (redundant_blk_avail && context_ptr->redundant_blk) {
            // Copy results
            BlkStruct *src_cu = &context_ptr->md_blk_arr_nsq[redundant_blk_mds];
            BlkStruct *dst_cu = blk_ptr;
            move_blk_data_redund(pcs_ptr, context_ptr, src_cu, dst_cu);
            memcpy(&context_ptr->md_local_blk_unit[blk_ptr->mds_idx],
                   &context_ptr->md_local_blk_unit[redundant_blk_mds],
                   sizeof(MdBlkStruct));

            if (!context_ptr->hbd_mode_decision) {
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_left_recon[0],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_left_recon[0],
                       128);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_left_recon[1],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_left_recon[1],
                       128);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_left_recon[2],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_left_recon[2],
                       128);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_top_recon[0],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_top_recon[0],
                       128);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_top_recon[1],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_top_recon[1],
                       128);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_top_recon[2],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_top_recon[2],
                       128);
            } else {
                uint16_t sz = sizeof(uint16_t);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_left_recon_16bit[0],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_left_recon_16bit[0],
                       128 * sz);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_left_recon_16bit[1],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_left_recon_16bit[1],
                       128 * sz);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_left_recon_16bit[2],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_left_recon_16bit[2],
                       128 * sz);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_top_recon_16bit[0],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_top_recon_16bit[0],
                       128 * sz);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_top_recon_16bit[1],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_top_recon_16bit[1],
                       128 * sz);
                memcpy(&context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                            .neigh_top_recon_16bit[2],
                       &context_ptr->md_local_blk_unit[redundant_blk_mds].neigh_top_recon_16bit[2],
                       128 * sz);
            }

            memcpy(&context_ptr->md_ep_pipe_sb[blk_ptr->mds_idx],
                   &context_ptr->md_ep_pipe_sb[redundant_blk_mds],
                   sizeof(MdEncPassCuData));

            if (context_ptr->blk_geom->shape == PART_N) {
                uint8_t sq_index = eb_log2f(context_ptr->blk_geom->sq_size) - 2;
                context_ptr->parent_sq_type[sq_index]      = src_cu->prediction_mode_flag;
                context_ptr->parent_sq_has_coeff[sq_index] = src_cu->block_has_coeff;
                context_ptr->parent_sq_pred_mode[sq_index] = src_cu->pred_mode;
            }
        } else {
            // Initialize tx_depth
            blk_ptr->tx_depth = 0;
            if (blk_geom->quadi > 0 && d1_block_itr == 0 && !skip_next_sq) {
                uint32_t            blk_mds           = context_ptr->blk_geom->sqi_mds;
                uint64_t            parent_depth_cost = 0, current_depth_cost = 0;
                SequenceControlSet *scs_ptr = (SequenceControlSet *)
                                                  pcs_ptr->scs_wrapper_ptr->object_ptr;
                uint32_t parent_depth_idx_mds = blk_mds;

                // from a given child index, derive the index of the parent
                parent_depth_idx_mds =
                    (context_ptr->blk_geom->sqi_mds -
                     (context_ptr->blk_geom->quadi - 3) *
                         ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128]
                                        [context_ptr->blk_geom->depth]) -
                    parent_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128]
                                       [blk_geom->depth];

                if (pcs_ptr->slice_type == I_SLICE && parent_depth_idx_mds == 0 &&
                    scs_ptr->seq_header.sb_size == BLOCK_128X128)
                    parent_depth_cost = MAX_MODE_COST;
                else
                    compute_depth_costs_md_skip(
                        context_ptr,
                        scs_ptr,
                        pcs_ptr->parent_pcs_ptr,
                        parent_depth_idx_mds,
                        ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128]
                                       [context_ptr->blk_geom->depth],
                        &parent_depth_cost,
                        &current_depth_cost);

                if (!pcs_ptr->parent_pcs_ptr->sb_geom[sb_addr]
                         .block_is_allowed[parent_depth_idx_mds])
                    parent_depth_cost = MAX_MODE_COST;

                // compare the cost of the parent to the cost of the already encoded child + an estimated cost for the remaining child @ the current depth
                // if the total child cost is higher than the parent cost then skip the remaining  child @ the current depth
                // when md_exit_th=0 the estimated cost for the remaining child is not taken into account and the action will be lossless compared to no exit
                // MD_EXIT_THSL could be tuned toward a faster encoder but lossy
                if (parent_depth_cost != MAX_MODE_COST &&
                    parent_depth_cost <= current_depth_cost +
                            (current_depth_cost * (4 - context_ptr->blk_geom->quadi) *
                             context_ptr->md_exit_th / context_ptr->blk_geom->quadi / 100)) {
                    skip_next_sq              = 1;
                    next_non_skip_blk_idx_mds = parent_depth_idx_mds +
                        ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128]
                                       [context_ptr->blk_geom->depth - 1];
                } else
                    skip_next_sq = 0;
            }
            // skip until we reach the next block @ the parent block depth
            if (blk_ptr->mds_idx >= next_non_skip_blk_idx_mds && skip_next_sq == 1)
                skip_next_sq = 0;

            uint8_t sq_weight_based_nsq_skip = update_skip_nsq_shapes(
                scs_ptr, pcs_ptr, context_ptr);
            skip_next_depth = context_ptr->blk_ptr->do_not_process_block;
            if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[blk_ptr->mds_idx] &&
                !skip_next_nsq && !skip_next_sq && !sq_weight_based_nsq_skip && !skip_next_depth) {

                if (context_ptr->blk_geom->shape != PART_N) {
                    // Don't do the prune if parent sq block is out of boundary
                    if (!pcs_ptr->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[context_ptr->blk_geom->sqi_mds])
                        context_ptr->prune_ref_frame_for_rec_partitions = 0;
                }
                md_encode_block(pcs_ptr, context_ptr, input_picture_ptr, bestcandidate_buffers);
            } else if (sq_weight_based_nsq_skip || skip_next_depth) {
                if (context_ptr->blk_geom->shape != PART_N)
                    context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].cost =
                        (MAX_MODE_COST >> 4);
                else
                    context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].cost =
                        (MAX_MODE_COST >> 10);
                context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].default_cost =
                    MAX_MODE_COST;
            } else if (skip_next_sq) {
                context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].cost =
                    (MAX_MODE_COST >> 10);
                context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].default_cost =
                    MAX_MODE_COST;
            } else {
                // If the block is out of the boundaries, md is not performed.
                // - For square blocks, since the blocks can be further splitted, they are considered in d2_inter_depth_block_decision with cost of zero.
                // - For non square blocks, since they can not be splitted further the cost is set to a large value (MAX_MODE_COST >> 4) to make sure they are not selected.
                //   The value is set to MAX_MODE_COST >> 4 to make sure there is not overflow when adding costs.
                if (context_ptr->blk_geom->shape != PART_N)
                    context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].cost =
                        (MAX_MODE_COST >> 4);
                else
                    context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].cost = 0;
                if (context_ptr->blk_geom->shape != PART_N)
                    context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].default_cost =
                        MAX_MODE_COST;
                else
                    context_ptr->md_local_blk_unit[context_ptr->blk_ptr->mds_idx].default_cost = 0;
            }
        }
        skip_next_nsq = 0;
        if (blk_geom->nsi + 1 == blk_geom->totns) {
            nsq_cost[context_ptr->blk_geom->shape] = d1_non_square_block_decision(context_ptr,
                                                                                  d1_block_itr);
            d1_block_itr++;
        } else if (d1_block_itr) {
            uint64_t tot_cost      = 0;
            uint32_t first_blk_idx = context_ptr->blk_ptr->mds_idx -
                (blk_geom->nsi); //index of first block in this partition
            for (int blk_it = 0; blk_it < blk_geom->nsi + 1; blk_it++)
                tot_cost += context_ptr->md_local_blk_unit[first_blk_idx + blk_it].cost;
            nsq_cost[context_ptr->blk_geom->shape] = tot_cost;
            if ((tot_cost +
                 tot_cost * (blk_geom->totns - (blk_geom->nsi + 1)) * context_ptr->md_exit_th /
                     (blk_geom->nsi + 1) / 100) >
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->sqi_mds].cost)
                skip_next_nsq = 1;
        }
        if (blk_geom->shape != PART_N) {
            if (blk_geom->nsi + 1 < blk_geom->totns)
                md_update_all_neighbour_arrays(
                    pcs_ptr, context_ptr, blk_idx_mds, sb_origin_x, sb_origin_y);
            else
                copy_neighbour_arrays( //restore [1] in [0] after done last ns block
                    pcs_ptr,
                    context_ptr,
                    1,
                    0,
                    blk_geom->sqi_mds,
                    sb_origin_x,
                    sb_origin_y);
        }
        d1_blocks_accumlated = d1_first_block == 1 ? 1 : d1_blocks_accumlated + 1;
        if (d1_blocks_accumlated == leaf_data_ptr->tot_d1_blocks) {
            //Sorting
            {
                uint32_t i, j, index;
                for (i = 0; i < NUMBER_OF_SHAPES - 1; ++i) {
                    for (j = i + 1; j < NUMBER_OF_SHAPES; ++j) {
                        if (nsq_cost[nsq_shape_table[j]] < nsq_cost[nsq_shape_table[i]]) {
                            index              = nsq_shape_table[i];
                            nsq_shape_table[i] = nsq_shape_table[j];
                            nsq_shape_table[j] = index;
                        }
                    }
                }
                depth_cost[scs_ptr->static_config.super_block_size == 128
                               ? context_ptr->blk_geom->depth
                               : context_ptr->blk_geom->depth + 1] += nsq_cost[nsq_shape_table[0]];
                if (context_ptr->skip_depth && scs_ptr->sb_geom[sb_addr].is_complete_sb) {
                    if (context_ptr->pd_pass > PD_PASS_1) {
                        uint64_t sq_cost       = nsq_cost[0]; // sq cost
                        uint64_t best_nsq_cost = MAX_CU_COST;
                        skip_next_depth        = 0;
                        // Derive best nsq cost
                        for (i = 1; i < NUMBER_OF_SHAPES; ++i)
                            if (nsq_cost[i] < best_nsq_cost)
                                best_nsq_cost = nsq_cost[i];
                        // Compare sq vs best nsq
                        uint64_t th = 30;
                        if (best_nsq_cost != MAX_CU_COST) {
                            if (sq_cost < best_nsq_cost) {
                                if ((best_nsq_cost - sq_cost) * 100 > (sq_cost * th)) {
                                    set_child_to_be_skipped(context_ptr,
                                                            context_ptr->blk_geom->sqi_mds,
                                                            scs_ptr->seq_header.sb_size,
                                                            1);
                                }
                            }
                        }
                    }
                }
            }

            uint32_t last_blk_index_mds = d2_inter_depth_block_decision(
                context_ptr,
                blk_geom->sqi_mds, //input is parent square
                sb_ptr,
                sb_addr,
                sb_origin_x,
                sb_origin_y,
                full_lambda,
                context_ptr->md_rate_estimation_ptr,
                pcs_ptr);
            d1_block_itr   = 0;
            d1_first_block = 1;
            if (context_ptr->md_blk_arr_nsq[last_blk_index_mds].split_flag == EB_FALSE) {
                md_update_all_neighbour_arrays_multiple(
                    pcs_ptr,
                    context_ptr,
                    context_ptr->md_local_blk_unit[last_blk_index_mds].best_d1_blk,
                    sb_origin_x,
                    sb_origin_y);
            }
        } else if (d1_first_block)
            d1_first_block = 0;
        blk_index++;
    } while (blk_index < leaf_count); // End of CU loop

    if (scs_ptr->seq_header.sb_size == BLOCK_64X64)
        depth_cost[0] = MAX_CU_COST;

    for (uint8_t depth_idx = 0; depth_idx < NUMBER_OF_DEPTH; depth_idx++) {
        sb_ptr->depth_cost[depth_idx] = depth_cost[depth_idx] < 0 ? MAX_MODE_COST
                                                                  : depth_cost[depth_idx];
    }

    return return_error;
}
#define MAX_SEARCH_POINT_WIDTH 128
#define MAX_SEARCH_POINT_HEIGHT 128

#define MAX_TATAL_SEARCH_AREA_WIDTH (MAX_SB_SIZE + MAX_SEARCH_POINT_WIDTH + ME_FILTER_TAP)
#define MAX_TATAL_SEARCH_AREA_HEIGHT (MAX_SB_SIZE + MAX_SEARCH_POINT_HEIGHT + ME_FILTER_TAP)

#define MAX_SEARCH_AREA_SIZE MAX_TATAL_SEARCH_AREA_WIDTH *MAX_TATAL_SEARCH_AREA_HEIGHT
