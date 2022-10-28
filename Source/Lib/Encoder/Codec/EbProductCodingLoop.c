/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 3-Clause Clear License and
* the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
#include "mv.h"
#include "mcomp.h"
#include "av1me.h"
#include "limits.h"

#include "EbPackUnPack_C.h"

#include "EbInterPrediction.h"
#define DIVIDE_AND_ROUND(x, y) (((x) + ((y) >> 1)) / (y))
void pack_block(uint8_t *in8_bit_buffer, uint32_t in8_stride, uint8_t *inn_bit_buffer,
                uint32_t inn_stride, uint16_t *out16_bit_buffer, uint32_t out_stride,
                uint32_t width, uint32_t height);
void svt_init_mv_cost_params(MV_COST_PARAMS *mv_cost_params, ModeDecisionContext *context_ptr,
                             const MV *ref_mv, uint8_t base_q_idx, uint32_t rdmult,
                             uint8_t hbd_mode_decision);
int  fp_mv_err_cost(const MV *mv, const MV_COST_PARAMS *mv_cost_params);
EbErrorType generate_md_stage_0_cand(SuperBlock *sb_ptr, ModeDecisionContext *context_ptr,
                                     uint32_t          *fast_candidate_total_count,
                                     PictureControlSet *pcs_ptr);
void        generate_md_stage_0_cand_light_pd1(SuperBlock *sb_ptr, ModeDecisionContext *context_ptr,
                                               uint32_t          *fast_candidate_total_count,
                                               PictureControlSet *pcs_ptr);
EbErrorType generate_md_stage_0_cand_light_pd0(ModeDecisionContext *context_ptr,
                                               uint32_t            *fast_candidate_total_count,
                                               PictureControlSet   *pcs_ptr);

static INLINE int is_interintra_allowed_bsize(const BlockSize bsize) {
    return (bsize >= BLOCK_8X8) && (bsize <= BLOCK_32X32);
}
void precompute_intra_pred_for_inter_intra(PictureControlSet   *pcs_ptr,
                                           ModeDecisionContext *context_ptr);

int svt_av1_allow_palette(int allow_palette, BlockSize sb_type);

void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr, Bool is_highbd);
extern IntraSize       intra_unit[];
EbPictureBufferDesc   *get_ref_pic_buffer(PictureControlSet *pcs_ptr, uint8_t is_highbd,
                                          uint8_t list_idx, uint8_t ref_idx);
const EbPredictionFunc svt_product_prediction_fun_table_light_pd0[2] = {
    svt_av1_intra_prediction_cl, inter_pu_prediction_av1_light_pd0};
const EbPredictionFunc svt_product_prediction_fun_table_light_pd1[2] = {
    svt_av1_intra_prediction_cl, inter_pu_prediction_av1_light_pd1};
const EbPredictionFunc svt_product_prediction_fun_table[2] = {svt_av1_intra_prediction_cl,
                                                              inter_pu_prediction_av1};

const EbFastCostFunc av1_product_fast_cost_func_table[2] = {
    av1_intra_fast_cost, /*INTRA */
    av1_inter_fast_cost /*INTER */
};

const EbAv1FullCostFunc svt_av1_product_full_cost_func_table[2] = {
    av1_intra_full_cost, /*INTRA */
    av1_inter_full_cost /*INTER */
};
/***************************************************
* Update Recon Samples Neighbor Arrays
***************************************************/
void mode_decision_update_neighbor_arrays_light_pd0(ModeDecisionContext *context_ptr) {
    uint32_t bwidth  = context_ptr->blk_geom->bwidth;
    uint32_t bheight = context_ptr->blk_geom->bheight;

    uint32_t origin_x = context_ptr->blk_origin_x;
    uint32_t origin_y = context_ptr->blk_origin_y;

    if (!context_ptr->skip_intra) {
        if (!context_ptr->hbd_mode_decision) {
            update_recon_neighbor_array(
                context_ptr->luma_recon_neighbor_array,
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon[0],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[0],
                origin_x,
                origin_y,
                bwidth,
                bheight);
        } else {
            update_recon_neighbor_array16bit(
                context_ptr->luma_recon_neighbor_array16bit,
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon_16bit[0],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon_16bit[0],
                origin_x,
                origin_y,
                bwidth,
                bheight);
        }
    }

    return;
}
/***************************************************
* Update Recon Samples Neighbor Arrays
***************************************************/
void mode_decision_update_neighbor_arrays(PictureControlSet   *pcs_ptr,
                                          ModeDecisionContext *context_ptr, uint32_t index_mds) {
    uint32_t bwdith  = context_ptr->blk_geom->bwidth;
    uint32_t bheight = context_ptr->blk_geom->bheight;

    uint32_t origin_x        = context_ptr->blk_origin_x;
    uint32_t origin_y        = context_ptr->blk_origin_y;
    uint32_t blk_origin_x_uv = context_ptr->round_origin_x >> 1;
    uint32_t blk_origin_y_uv = context_ptr->round_origin_y >> 1;
    uint32_t bwdith_uv       = context_ptr->blk_geom->bwidth_uv;
    uint32_t bwheight_uv     = context_ptr->blk_geom->bheight_uv;
    (void)index_mds;

    int32_t is_inter = (context_ptr->blk_ptr->prediction_mode_flag == INTER_MODE ||
                        context_ptr->blk_ptr->use_intrabc)
        ? TRUE
        : FALSE;

    uint16_t tile_idx = context_ptr->tile_index;

    {
        if (!(context_ptr->pd_pass == PD_PASS_1 && context_ptr->pic_pred_depth_only)) {
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
        }
        if (context_ptr->rate_est_ctrls.update_skip_coeff_ctx) {
            uint8_t is_skip_coeff = !context_ptr->blk_ptr->block_has_coeff;

            neighbor_array_unit_mode_write(context_ptr->skip_coeff_neighbor_array,
                                           &is_skip_coeff,
                                           origin_x,
                                           origin_y,
                                           bwdith,
                                           bheight,
                                           NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }
        if (context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx) {
            uint16_t txb_count = context_ptr->blk_geom->txb_count[context_ptr->blk_ptr->tx_depth];
            for (uint8_t txb_itr = 0; txb_itr < txb_count; txb_itr++) {
                uint8_t dc_sign_level_coeff =
                    (uint8_t)context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .quantized_dc[0][txb_itr];
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
                    pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array
                        [MD_NEIGHBOR_ARRAY_INDEX][tile_idx],
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
    }
    if (context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx)
        if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            //  Update chroma CB cbf and Dc context
            {
                uint8_t dc_sign_level_coeff =
                    (uint8_t)context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .quantized_dc[1][0];
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
                uint8_t dc_sign_level_coeff =
                    (uint8_t)context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .quantized_dc[2][0];
                neighbor_array_unit_mode_write(context_ptr->cr_dc_sign_level_coeff_neighbor_array,
                                               (uint8_t *)&dc_sign_level_coeff,
                                               blk_origin_x_uv,
                                               blk_origin_y_uv,
                                               bwdith_uv,
                                               bwheight_uv,
                                               NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
            }
        }
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.tx_mode == TX_MODE_SELECT) {
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
    }
    if (!context_ptr->skip_intra) {
        if (context_ptr->encoder_bit_depth > EB_EIGHT_BIT && context_ptr->bypass_encdec &&
            !context_ptr->hbd_mode_decision && context_ptr->pd_pass == PD_PASS_1) {
            // copy HBD
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

            if (context_ptr->txs_ctrls.enabled) {
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

            if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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

            // copy 8 bit
            update_recon_neighbor_array(
                context_ptr->luma_recon_neighbor_array,
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon[0],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[0],
                origin_x,
                origin_y,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight);

            if (context_ptr->txs_ctrls.enabled) {
                update_recon_neighbor_array(
                    pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                    [tile_idx],
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_top_recon[0],
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon[0],
                    origin_x,
                    origin_y,
                    context_ptr->blk_geom->bwidth,
                    context_ptr->blk_geom->bheight);
                update_recon_neighbor_array(
                    pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                    [tile_idx],
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_top_recon[0],
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon[0],
                    origin_x,
                    origin_y,
                    context_ptr->blk_geom->bwidth,
                    context_ptr->blk_geom->bheight);
            }

            if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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
        } else if (!context_ptr->hbd_mode_decision) {
            update_recon_neighbor_array(
                context_ptr->luma_recon_neighbor_array,
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_top_recon[0],
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[0],
                origin_x,
                origin_y,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight);
            if (context_ptr->txs_ctrls.enabled) {
                update_recon_neighbor_array(
                    pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                    [tile_idx],
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_top_recon[0],
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon[0],
                    origin_x,
                    origin_y,
                    context_ptr->blk_geom->bwidth,
                    context_ptr->blk_geom->bheight);
                update_recon_neighbor_array(
                    pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                    [tile_idx],
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_top_recon[0],
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .neigh_left_recon[0],
                    origin_x,
                    origin_y,
                    context_ptr->blk_geom->bwidth,
                    context_ptr->blk_geom->bheight);
            }
            if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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
            if (context_ptr->txs_ctrls.enabled) {
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
            if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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
    }
    return;
}
void copy_neighbour_arrays_light_pd0(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                     uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds,
                                     uint32_t sb_org_x, uint32_t sb_org_y) {
    uint16_t tile_idx = context_ptr->tile_index;

    const BlockGeom *blk_geom = get_blk_geom_mds(blk_mds);

    uint32_t blk_org_x = sb_org_x + blk_geom->origin_x;
    uint32_t blk_org_y = sb_org_y + blk_geom->origin_y;

    if (!context_ptr->hbd_mode_decision) {
        copy_neigh_arr(pcs_ptr->md_luma_recon_neighbor_array[src_idx][tile_idx],
                       pcs_ptr->md_luma_recon_neighbor_array[dst_idx][tile_idx],
                       blk_org_x,
                       blk_org_y,
                       blk_geom->bwidth,
                       blk_geom->bheight,
                       NEIGHBOR_ARRAY_UNIT_FULL_MASK);
    } else {
        copy_neigh_arr(pcs_ptr->md_luma_recon_neighbor_array16bit[src_idx][tile_idx],
                       pcs_ptr->md_luma_recon_neighbor_array16bit[dst_idx][tile_idx],
                       blk_org_x,
                       blk_org_y,
                       blk_geom->bwidth,
                       blk_geom->bheight,
                       NEIGHBOR_ARRAY_UNIT_FULL_MASK);
    }
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
        if (context_ptr->txs_ctrls.enabled) {
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
        if (blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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
        if (context_ptr->txs_ctrls.enabled) {
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
        if (blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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
    if (blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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
    if (context_ptr->rate_est_ctrls.update_skip_coeff_ctx) {
        copy_neigh_arr(pcs_ptr->md_skip_coeff_neighbor_array[src_idx][tile_idx],
                       pcs_ptr->md_skip_coeff_neighbor_array[dst_idx][tile_idx],
                       blk_org_x,
                       blk_org_y,
                       blk_geom->bwidth,
                       blk_geom->bheight,
                       NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }
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
    uint8_t avail_blk_flag = context_ptr->avail_blk_flag[last_blk_index_mds];
    if (avail_blk_flag) {
        mode_decision_update_neighbor_arrays(pcs_ptr, context_ptr, last_blk_index_mds);
        if (!context_ptr->shut_fast_rate ||
            context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx ||
            context_ptr->cand_reduction_ctrls.use_neighbouring_mode_ctrls.enabled)
            update_mi_map(context_ptr->blk_ptr,
                          context_ptr->blk_origin_x,
                          context_ptr->blk_origin_y,
                          context_ptr->blk_geom,
                          pcs_ptr);
    }
}

void md_update_all_neighbour_arrays_multiple(PictureControlSet   *pcs_ptr,
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
/************************************************************************************************
* av1_perform_inverse_transform_recon_luma
* Apply inverse transform for Luma samples
************************************************************************************************/
void av1_perform_inverse_transform_recon_luma(ModeDecisionContext         *context_ptr,
                                              ModeDecisionCandidateBuffer *candidate_buffer) {
    uint32_t tu_total_count;
    uint32_t txb_itr;

    uint8_t tx_depth         = candidate_buffer->candidate_ptr->tx_depth;
    tu_total_count           = context_ptr->blk_geom->txb_count[tx_depth];
    txb_itr                  = 0;
    uint32_t   txb_1d_offset = 0;
    const Bool is_inter      = (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) ||
                           candidate_buffer->candidate_ptr->use_intrabc)
             ? TRUE
             : FALSE;
    do {
        uint32_t txb_origin_x     = context_ptr->blk_geom->tx_org_x[is_inter][tx_depth][txb_itr];
        uint32_t txb_origin_y     = context_ptr->blk_geom->tx_org_y[is_inter][tx_depth][txb_itr];
        uint32_t txb_width        = context_ptr->blk_geom->tx_width[tx_depth][txb_itr];
        uint32_t txb_height       = context_ptr->blk_geom->tx_height[tx_depth][txb_itr];
        uint32_t txb_origin_index = txb_origin_x +
            txb_origin_y * candidate_buffer->prediction_ptr->stride_y;
        uint32_t rec_luma_offset = txb_origin_x +
            txb_origin_y * candidate_buffer->recon_ptr->stride_y;
        uint32_t y_has_coeff = (candidate_buffer->y_has_coeff & (1 << txb_itr)) > 0;
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
                                        (uint32_t)candidate_buffer->eob[0][txb_itr]);
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
void av1_perform_inverse_transform_recon(PictureControlSet *pcs, ModeDecisionContext *context_ptr,
                                         ModeDecisionCandidateBuffer *candidate_buffer,
                                         const BlockGeom             *blk_geom) {
    uint32_t tu_total_count;
    uint32_t txb_index;
    uint32_t txb_itr;

    UNUSED(blk_geom);

    const uint8_t tx_depth   = candidate_buffer->candidate_ptr->tx_depth;
    tu_total_count           = context_ptr->blk_geom->txb_count[tx_depth];
    txb_index                = 0;
    txb_itr                  = 0;
    uint32_t   txb_1d_offset = 0, txb_1d_offset_uv = 0;
    const Bool is_inter = is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) ||
        candidate_buffer->candidate_ptr->use_intrabc;
    do {
        uint32_t txb_origin_x = context_ptr->blk_geom->tx_org_x[is_inter][tx_depth][txb_itr];
        uint32_t txb_origin_y = context_ptr->blk_geom->tx_org_y[is_inter][tx_depth][txb_itr];
        uint32_t txb_width    = context_ptr->blk_geom->tx_width[tx_depth][txb_itr];
        TxSize   tx_size      = context_ptr->blk_geom->txsize[tx_depth][txb_itr];
        if (context_ptr->md_staging_subres_step == 2) {
            if (tx_size == TX_64X64)
                tx_size = TX_64X16;
            else if (tx_size == TX_32X32)
                tx_size = TX_32X8;
            else if (tx_size == TX_16X16)
                tx_size = TX_16X4;
            else
                assert(0);
        } else if (context_ptr->md_staging_subres_step == 1) {
            if (tx_size == TX_64X64)
                tx_size = TX_64X32;
            else if (tx_size == TX_32X32)
                tx_size = TX_32X16;
            else if (tx_size == TX_16X16)
                tx_size = TX_16X8;
            else if (tx_size == TX_8X8)
                tx_size = TX_8X4;
            else
                assert(0);
        }
        uint32_t txb_height = context_ptr->blk_geom->tx_height[tx_depth][txb_itr] >>
            context_ptr->md_staging_subres_step;
        uint32_t rec_luma_offset = txb_origin_x +
            txb_origin_y * candidate_buffer->recon_ptr->stride_y;
        uint32_t rec_cb_offset    = ((((txb_origin_x >> 3) << 3) +
                                   ((txb_origin_y >> 3) << 3) *
                                       candidate_buffer->recon_ptr->stride_cb) >>
                                  1);
        uint32_t rec_cr_offset    = ((((txb_origin_x >> 3) << 3) +
                                   ((txb_origin_y >> 3) << 3) *
                                       candidate_buffer->recon_ptr->stride_cr) >>
                                  1);
        uint32_t txb_origin_index = txb_origin_x +
            txb_origin_y * candidate_buffer->prediction_ptr->stride_y;
        EbPictureBufferDesc *recon_buffer = candidate_buffer->recon_ptr;

        // If bypassing encdec, update the recon pointer to copy the recon directly
        // into the buffer used for EncDec; avoids copy after this function call.
        // candidate_buffer->recon_ptr is only used to update other buffers after this point.
        if (context_ptr->bypass_encdec && context_ptr->pd_pass == PD_PASS_1) {
            if (context_ptr->pred_depth_only && context_ptr->md_disallow_nsq) {
                get_recon_pic(pcs, &recon_buffer, context_ptr->hbd_mode_decision);
                uint16_t org_x = context_ptr->blk_origin_x +
                    (blk_geom->tx_org_x[is_inter][tx_depth][txb_itr] - blk_geom->origin_x);
                uint16_t org_y = context_ptr->blk_origin_y +
                    (blk_geom->tx_org_y[is_inter][tx_depth][txb_itr] - blk_geom->origin_y);

                rec_luma_offset = (recon_buffer->origin_y + org_y) * recon_buffer->stride_y +
                    (recon_buffer->origin_x + org_x);

                uint32_t round_origin_x = (org_x >> 3) << 3; // for Chroma blocks with size of 4
                uint32_t round_origin_y = (org_y >> 3) << 3; // for Chroma blocks with size of 4
                rec_cb_offset = rec_cr_offset = (((recon_buffer->origin_y + round_origin_y) >> 1) *
                                                 recon_buffer->stride_cb) +
                    ((recon_buffer->origin_x + round_origin_x) >> 1);
            } else {
                recon_buffer    = context_ptr->blk_ptr->recon_tmp;
                rec_luma_offset = (txb_origin_x - blk_geom->origin_x) +
                    (txb_origin_y - blk_geom->origin_y) * recon_buffer->stride_y;
                // Chroma only copied for txb_itr == 0; not offset for the recon_tmp buffer b/c it only stores the block size
                rec_cb_offset = rec_cr_offset = 0;
            }
        }
        if (context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                .y_has_coeff[txb_itr]) {
            inv_transform_recon_wrapper(
                candidate_buffer->prediction_ptr->buffer_y,
                txb_origin_index,
                candidate_buffer->prediction_ptr->stride_y << context_ptr->md_staging_subres_step,
                recon_buffer->buffer_y,
                rec_luma_offset,
                recon_buffer->stride_y << context_ptr->md_staging_subres_step,
                (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_y,
                txb_1d_offset,
                context_ptr->hbd_mode_decision,
                tx_size,
                candidate_buffer->candidate_ptr->transform_type[txb_itr],
                PLANE_TYPE_Y,
                (uint32_t)candidate_buffer->eob[0][txb_itr]);
            if (context_ptr->md_staging_subres_step == 2) {
                for (uint32_t i = 0; i < (txb_height * 4); i += 4) {
                    if (context_ptr->hbd_mode_decision) {
                        EB_MEMCPY(((uint16_t *)recon_buffer->buffer_y) + rec_luma_offset +
                                      (i + 1) * recon_buffer->stride_y,
                                  ((uint16_t *)recon_buffer->buffer_y) + rec_luma_offset +
                                      i * recon_buffer->stride_y,
                                  txb_width * sizeof(uint16_t));
                        EB_MEMCPY(((uint16_t *)recon_buffer->buffer_y) + rec_luma_offset +
                                      (i + 2) * recon_buffer->stride_y,
                                  ((uint16_t *)recon_buffer->buffer_y) + rec_luma_offset +
                                      i * recon_buffer->stride_y,
                                  txb_width * sizeof(uint16_t));
                        EB_MEMCPY(((uint16_t *)recon_buffer->buffer_y) + rec_luma_offset +
                                      (i + 3) * recon_buffer->stride_y,
                                  ((uint16_t *)recon_buffer->buffer_y) + rec_luma_offset +
                                      i * recon_buffer->stride_y,
                                  txb_width * sizeof(uint16_t));
                    } else {
                        EB_MEMCPY(
                            recon_buffer->buffer_y + rec_luma_offset +
                                (i + 1) * recon_buffer->stride_y,
                            recon_buffer->buffer_y + rec_luma_offset + i * recon_buffer->stride_y,
                            txb_width);
                        EB_MEMCPY(
                            recon_buffer->buffer_y + rec_luma_offset +
                                (i + 2) * recon_buffer->stride_y,
                            recon_buffer->buffer_y + rec_luma_offset + i * recon_buffer->stride_y,
                            txb_width);
                        EB_MEMCPY(
                            recon_buffer->buffer_y + rec_luma_offset +
                                (i + 3) * recon_buffer->stride_y,
                            recon_buffer->buffer_y + rec_luma_offset + i * recon_buffer->stride_y,
                            txb_width);
                    }
                }
            } else if (context_ptr->md_staging_subres_step) {
                for (uint32_t i = 0; i < (txb_height * 2); i += 2) {
                    if (context_ptr->hbd_mode_decision)
                        EB_MEMCPY(((uint16_t *)recon_buffer->buffer_y) + rec_luma_offset +
                                      (i + 1) * recon_buffer->stride_y,
                                  ((uint16_t *)recon_buffer->buffer_y) + rec_luma_offset +
                                      i * recon_buffer->stride_y,
                                  txb_width * sizeof(uint16_t));
                    else
                        EB_MEMCPY(
                            recon_buffer->buffer_y + rec_luma_offset +
                                (i + 1) * recon_buffer->stride_y,
                            recon_buffer->buffer_y + rec_luma_offset + i * recon_buffer->stride_y,
                            txb_width);
                }
            }
        } else
            svt_av1_picture_copy(candidate_buffer->prediction_ptr,
                                 txb_origin_index,
                                 0, //txb_chroma_origin_index,
                                 recon_buffer,
                                 rec_luma_offset,
                                 0, //txb_chroma_origin_index,
                                 txb_width,
                                 txb_height << context_ptr->md_staging_subres_step,
                                 0, //chromaTuSize,
                                 0, //chromaTuSize,
                                 PICTURE_BUFFER_DESC_Y_FLAG,
                                 context_ptr->hbd_mode_decision);

        //CHROMA
        if (tx_depth == 0 || txb_itr == 0) {
            if (context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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
                        recon_buffer->buffer_cb,
                        rec_cb_offset,
                        recon_buffer->stride_cb,
                        (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cb,
                        txb_1d_offset_uv,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                        candidate_buffer->candidate_ptr->transform_type_uv,
                        PLANE_TYPE_UV,
                        (uint32_t)candidate_buffer->eob[1][txb_itr]);
                else
                    svt_av1_picture_copy(candidate_buffer->prediction_ptr,
                                         0,
                                         cb_tu_chroma_origin_index,
                                         recon_buffer,
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
                        recon_buffer->buffer_cr,
                        rec_cr_offset,
                        recon_buffer->stride_cr,
                        (int32_t *)candidate_buffer->recon_coeff_ptr->buffer_cr,
                        txb_1d_offset_uv,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                        candidate_buffer->candidate_ptr->transform_type_uv,
                        PLANE_TYPE_UV,
                        (uint32_t)candidate_buffer->eob[2][txb_itr]);
                else
                    svt_av1_picture_copy(candidate_buffer->prediction_ptr,
                                         0,
                                         cr_tu_chroma_origin_index,
                                         recon_buffer,
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
void product_coding_loop_init_fast_loop(PictureControlSet   *pcs_ptr,
                                        ModeDecisionContext *context_ptr,
                                        NeighborArrayUnit   *skip_coeff_neighbor_array,
                                        NeighborArrayUnit   *leaf_partition_neighbor_array) {
    context_ptr->tx_depth = context_ptr->blk_ptr->tx_depth = 0;
    // Generate Split, Skip and intra mode contexts for the rate estimation
    coding_loop_context_generation(pcs_ptr,
                                   context_ptr,
                                   context_ptr->blk_ptr,
                                   context_ptr->blk_origin_x,
                                   context_ptr->blk_origin_y,
                                   skip_coeff_neighbor_array,
                                   leaf_partition_neighbor_array);

    return;
}

static uint32_t hadamard_path(ModeDecisionContext *ctx, uint8_t *input, uint32_t input_origin_index,
                              uint32_t input_stride, uint8_t *pred, uint32_t blk_origin_index,
                              uint32_t pred_stride) {
    uint32_t input_idx, pred_idx, res_idx;

    uint32_t satd_cost = 0;

    const TxSize tx_size = AOMMIN(TX_32X32, max_txsize_lookup[ctx->blk_geom->bsize]);

    const int stepr = tx_size_high_unit[tx_size];
    const int stepc = tx_size_wide_unit[tx_size];
    const int txbw  = tx_size_wide[tx_size];
    const int txbh  = tx_size_high[tx_size];

    const int max_blocks_wide = block_size_wide[ctx->blk_geom->bsize] >> MI_SIZE_LOG2;
    const int max_blocks_high = block_size_wide[ctx->blk_geom->bsize] >> MI_SIZE_LOG2;
    int       row, col;

    for (row = 0; row < max_blocks_high; row += stepr) {
        for (col = 0; col < max_blocks_wide; col += stepc) {
            input_idx = input_origin_index + (((row * input_stride) + col) << 2);
            pred_idx  = blk_origin_index + (((row * pred_stride) + col) << 2);
            res_idx   = 0;

            residual_kernel(input,
                            input_idx,
                            input_stride,
                            pred,
                            pred_idx,
                            pred_stride,
                            (int16_t *)ctx->temp_residual_ptr->buffer_y,
                            res_idx,
                            ctx->temp_residual_ptr->stride_y,
                            0, // input and pred 8-bit
                            txbw,
                            txbh);

            switch (tx_size) {
            case TX_4X4:
                svt_aom_hadamard_4x4((int16_t *)ctx->temp_residual_ptr->buffer_y,
                                     ctx->temp_residual_ptr->stride_y,
                                     &(((int32_t *)ctx->tx_coeffs->buffer_y)[0]));
                break;

            case TX_8X8:
                svt_aom_hadamard_8x8((int16_t *)ctx->temp_residual_ptr->buffer_y,
                                     ctx->temp_residual_ptr->stride_y,
                                     &(((int32_t *)ctx->tx_coeffs->buffer_y)[0]));
                break;

            case TX_16X16:
                svt_aom_hadamard_16x16((int16_t *)ctx->temp_residual_ptr->buffer_y,
                                       ctx->temp_residual_ptr->stride_y,
                                       &(((int32_t *)ctx->tx_coeffs->buffer_y)[0]));
                break;

            case TX_32X32:
                svt_aom_hadamard_32x32((int16_t *)ctx->temp_residual_ptr->buffer_y,
                                       ctx->temp_residual_ptr->stride_y,
                                       &(((int32_t *)ctx->tx_coeffs->buffer_y)[0]));
                break;

            default: assert(0);
            }
            satd_cost += svt_aom_satd(&(((int32_t *)ctx->tx_coeffs->buffer_y)[0]),
                                      tx_size_2d[tx_size]);
        }
    }
    return (satd_cost);
}

/*
*/
void fast_loop_core_light_pd0(ModeDecisionCandidateBuffer *candidate_buffer,
                              PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                              EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index,
                              uint32_t cu_origin_index) {
    ModeDecisionCandidate *candidate_ptr  = candidate_buffer->candidate_ptr;
    EbPictureBufferDesc   *prediction_ptr = candidate_buffer->prediction_ptr;

    if (context_ptr->lpd0_ctrls.pd0_level == VERY_LIGHT_PD0) {
        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, candidate_ptr->ref_frame_type);

        assert(rf[1] == NONE_FRAME);

        EbPictureBufferDesc *ref_pic;
        const int8_t         ref_idx_first  = get_ref_frame_idx(rf[0]);
        const int8_t         list_idx_first = get_list_idx(rf[0]);
        int32_t              ref_origin_index;
        int16_t              mv_x, mv_y;

        EbReferenceObject *ref_obj = (EbReferenceObject *)pcs_ptr
                                         ->ref_pic_ptr_array[list_idx_first][ref_idx_first]
                                         ->object_ptr;
        if (list_idx_first == 0) {
            ref_pic = get_ref_pic_buffer(pcs_ptr, context_ptr->hbd_mode_decision, 0, ref_idx_first);
            mv_x    = candidate_ptr->mv[REF_LIST_0].x >> 3;
            mv_y    = candidate_ptr->mv[REF_LIST_0].y >> 3;
        } else {
            ref_pic = get_ref_pic_buffer(pcs_ptr, context_ptr->hbd_mode_decision, 1, ref_idx_first);
            mv_x    = candidate_ptr->mv[REF_LIST_1].x >> 3;
            mv_y    = candidate_ptr->mv[REF_LIST_1].y >> 3;
        }
        // -------
        // Use scaled references if resolution of the reference is different from that of the input
        // -------
        use_scaled_rec_refs_if_needed(
            pcs_ptr, input_picture_ptr, ref_obj, &ref_pic, context_ptr->hbd_mode_decision);
        ref_origin_index = ref_pic->origin_x + (context_ptr->blk_origin_x + mv_x) +
            (context_ptr->blk_origin_y + mv_y + ref_pic->origin_y) * ref_pic->stride_y;

        EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
            ? svt_full_distortion_kernel16_bits
            : svt_spatial_full_distortion_kernel;

        *(candidate_buffer->fast_cost_ptr) = (uint32_t)(spatial_full_dist_type_fun(
                                                 input_picture_ptr->buffer_y,
                                                 input_origin_index,
                                                 input_picture_ptr->stride_y << 1,
                                                 ref_pic->buffer_y,
                                                 ref_origin_index,
                                                 ref_pic->stride_y << 1,
                                                 context_ptr->blk_geom->bwidth,
                                                 context_ptr->blk_geom->bheight >> 1))
            << 1;
    } else {
        // intrabc not allowed in light_pd0
        svt_product_prediction_fun_table_light_pd0[is_inter_mode(candidate_ptr->pred_mode)](
            context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
        if (context_ptr->mds0_ctrls.mds0_dist_type == MDS0_VAR) {
            if (!context_ptr->hbd_mode_decision) {
                const AomVarianceFnPtr *fn_ptr = &mefn_ptr[context_ptr->blk_geom->bsize];
                unsigned int            sse;
                uint8_t                *pred_y = prediction_ptr->buffer_y + cu_origin_index;
                uint8_t                *src_y  = input_picture_ptr->buffer_y + input_origin_index;
                *(candidate_buffer->fast_cost_ptr) = fn_ptr->vf(pred_y,
                                                                prediction_ptr->stride_y,
                                                                src_y,
                                                                input_picture_ptr->stride_y,
                                                                &sse) >>
                    2;
            } else {
                const AomVarianceFnPtr *fn_ptr = &mefn_ptr[context_ptr->blk_geom->bsize];
                unsigned int            sse;
                uint16_t *pred_y = ((uint16_t *)prediction_ptr->buffer_y) + cu_origin_index;
                uint16_t *src_y  = ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index;
                *(candidate_buffer->fast_cost_ptr) = fn_ptr->vf_hbd_10(CONVERT_TO_BYTEPTR(pred_y),
                                                                       prediction_ptr->stride_y,
                                                                       CONVERT_TO_BYTEPTR(src_y),
                                                                       input_picture_ptr->stride_y,
                                                                       &sse) >>
                    1;
            }

        } else {
            assert(context_ptr->mds0_ctrls.mds0_dist_type == MDS0_SAD);
            assert((context_ptr->blk_geom->bwidth >> 3) < 17);
            if (!context_ptr->hbd_mode_decision) {
                *(candidate_buffer->fast_cost_ptr) = svt_nxm_sad_kernel_sub_sampled(
                    input_picture_ptr->buffer_y + input_origin_index,
                    input_picture_ptr->stride_y,
                    prediction_ptr->buffer_y + cu_origin_index,
                    prediction_ptr->stride_y,
                    context_ptr->blk_geom->bheight,
                    context_ptr->blk_geom->bwidth);
            } else {
                *(candidate_buffer->fast_cost_ptr) = sad_16b_kernel(
                    ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                    input_picture_ptr->stride_y,
                    ((uint16_t *)prediction_ptr->buffer_y) + cu_origin_index,
                    prediction_ptr->stride_y,
                    context_ptr->blk_geom->bheight,
                    context_ptr->blk_geom->bwidth);
            }
        }
    }
}
// Light PD1 fast loop core; assumes luma only, 8bit only, and that SSD is not used.
void fast_loop_core_light_pd1(ModeDecisionCandidateBuffer *candidate_buffer,
                              PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                              EbPictureBufferDesc *input_picture_ptr, BlockLocation *loc,
                              BlkStruct *blk_ptr) {
    uint32_t       luma_fast_distortion;
    const uint32_t fast_lambda = context_ptr->fast_lambda_md[EB_8_BIT_MD];

    ModeDecisionCandidate *candidate_ptr  = candidate_buffer->candidate_ptr;
    EbPictureBufferDesc   *prediction_ptr = candidate_buffer->prediction_ptr;
    context_ptr->pu_itr                   = 0;
    // Prediction
    context_ptr->uv_intra_comp_only = FALSE;
    svt_product_prediction_fun_table_light_pd1[is_inter_mode(candidate_ptr->pred_mode)](
        0, context_ptr, pcs_ptr, candidate_buffer);
    // Distortion
    if (context_ptr->mds0_ctrls.mds0_dist_type == MDS0_VAR) {
        const AomVarianceFnPtr *fn_ptr = &mefn_ptr[context_ptr->blk_geom->bsize];
        unsigned int            sse;
        uint8_t                *pred_y = prediction_ptr->buffer_y + loc->blk_origin_index;
        uint8_t                *src_y  = input_picture_ptr->buffer_y + loc->input_origin_index;

        // The variance is shifted because fast_lambda is used, and variance is much larger than SAD (for which
        // fast_lambda was designed), so a scaling is needed to make the values closer.  3 was chosen empirically.
        candidate_buffer->luma_fast_distortion = luma_fast_distortion =
            fn_ptr->vf(
                pred_y, prediction_ptr->stride_y, src_y, input_picture_ptr->stride_y, &sse) >>
            3;
    } else {
        assert(context_ptr->mds0_ctrls.mds0_dist_type == MDS0_SAD);
        assert((context_ptr->blk_geom->bwidth >> 3) < 17);
        candidate_buffer->luma_fast_distortion =
            (uint32_t)(luma_fast_distortion = svt_nxm_sad_kernel_sub_sampled(
                           input_picture_ptr->buffer_y + loc->input_origin_index,
                           input_picture_ptr->stride_y,
                           prediction_ptr->buffer_y + loc->blk_origin_index,
                           prediction_ptr->stride_y,
                           context_ptr->blk_geom->bheight,
                           context_ptr->blk_geom->bwidth));
    }
    // If distortion cost is greater than the best cost, exit early. This candidate will never be
    // selected b/c only one candidate is sent to MDS3
    if (context_ptr->mds0_best_cost != (uint32_t)~0) {
        const uint64_t distortion_cost = RDCOST(fast_lambda, 0, luma_fast_distortion);
        if (distortion_cost > context_ptr->mds0_best_cost) {
            *(candidate_buffer->fast_cost_ptr) = MAX_MODE_COST;
            return;
        }
    }
    // Fast Cost
    if (context_ptr->shut_fast_rate) {
        *(candidate_buffer->fast_cost_ptr) = luma_fast_distortion;
        candidate_buffer->fast_luma_rate   = 0;
        candidate_buffer->fast_chroma_rate = 0;
    } else {
        *(candidate_buffer->fast_cost_ptr) =
            av1_product_fast_cost_func_table[is_inter_mode(candidate_ptr->pred_mode)](
                context_ptr,
                blk_ptr,
                candidate_buffer,
                NOT_USED_VALUE,
                luma_fast_distortion,
                0, //chroma_fast_distortion,
                fast_lambda,
                pcs_ptr,
                &(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                      .ed_ref_mv_stack[candidate_ptr->ref_frame_type][0]),
                context_ptr->blk_geom,
                context_ptr->blk_origin_y >> MI_SIZE_LOG2,
                context_ptr->blk_origin_x >> MI_SIZE_LOG2,
                0, // inter-intra enabled
                context_ptr->intra_luma_left_mode,
                context_ptr->intra_luma_top_mode);
    }
}
void fast_loop_core(ModeDecisionCandidateBuffer *candidate_buffer, PictureControlSet *pcs_ptr,
                    ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                    uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                    uint32_t input_cr_origin_in_index, BlkStruct *blk_ptr, uint32_t cu_origin_index,
                    uint32_t cu_chroma_origin_index) {
    uint32_t luma_fast_distortion;
    uint32_t chroma_fast_distortion = 0;
    uint32_t full_lambda            = context_ptr->hbd_mode_decision
                   ? context_ptr->full_lambda_md[EB_10_BIT_MD]
                   : context_ptr->full_lambda_md[EB_8_BIT_MD];
    uint32_t fast_lambda            = context_ptr->hbd_mode_decision
                   ? context_ptr->fast_lambda_md[EB_10_BIT_MD]
                   : context_ptr->fast_lambda_md[EB_8_BIT_MD];

    ModeDecisionCandidate *candidate_ptr  = candidate_buffer->candidate_ptr;
    EbPictureBufferDesc   *prediction_ptr = candidate_buffer->prediction_ptr;
    context_ptr->pu_itr                   = 0;
    // Prediction
    context_ptr->uv_intra_comp_only = FALSE;
    svt_product_prediction_fun_table[is_inter_mode(candidate_ptr->pred_mode) ||
                                     candidate_ptr->use_intrabc](
        context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
    // Distortion
    // Y
    if (context_ptr->mds0_ctrls.mds0_dist_type == MDS0_SSD) {
        EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
            ? svt_full_distortion_kernel16_bits
            : svt_spatial_full_distortion_kernel;
        candidate_buffer->luma_fast_distortion           = luma_fast_distortion =
            (uint32_t)(spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                                                  input_origin_index,
                                                  input_picture_ptr->stride_y,
                                                  prediction_ptr->buffer_y,
                                                  (int32_t)cu_origin_index,
                                                  prediction_ptr->stride_y,
                                                  context_ptr->blk_geom->bwidth,
                                                  context_ptr->blk_geom->bheight));
    } else if (context_ptr->mds0_ctrls.mds0_dist_type == MDS0_VAR) {
        if (!context_ptr->hbd_mode_decision) {
            const AomVarianceFnPtr *fn_ptr = &mefn_ptr[context_ptr->blk_geom->bsize];
            unsigned int            sse;
            uint8_t                *pred_y = prediction_ptr->buffer_y + cu_origin_index;
            uint8_t                *src_y  = input_picture_ptr->buffer_y + input_origin_index;
            candidate_buffer->luma_fast_distortion = luma_fast_distortion =
                fn_ptr->vf(
                    pred_y, prediction_ptr->stride_y, src_y, input_picture_ptr->stride_y, &sse) >>
                2;
        } else {
            const AomVarianceFnPtr *fn_ptr = &mefn_ptr[context_ptr->blk_geom->bsize];
            unsigned int            sse;
            uint16_t *pred_y = ((uint16_t *)prediction_ptr->buffer_y) + cu_origin_index;
            uint16_t *src_y  = ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index;
            candidate_buffer->luma_fast_distortion = luma_fast_distortion =
                fn_ptr->vf_hbd_10(CONVERT_TO_BYTEPTR(pred_y),
                                  prediction_ptr->stride_y,
                                  CONVERT_TO_BYTEPTR(src_y),
                                  input_picture_ptr->stride_y,
                                  &sse) >>
                1;
        }
    } else {
        assert(context_ptr->mds0_ctrls.mds0_dist_type == MDS0_SAD);
        assert((context_ptr->blk_geom->bwidth >> 3) < 17);
        if (!context_ptr->hbd_mode_decision) {
            candidate_buffer->luma_fast_distortion =
                (uint32_t)(luma_fast_distortion = svt_nxm_sad_kernel_sub_sampled(
                               input_picture_ptr->buffer_y + input_origin_index,
                               input_picture_ptr->stride_y,
                               prediction_ptr->buffer_y + cu_origin_index,
                               prediction_ptr->stride_y,
                               context_ptr->blk_geom->bheight,
                               context_ptr->blk_geom->bwidth));
        } else {
            candidate_buffer->luma_fast_distortion =
                (uint32_t)(luma_fast_distortion = sad_16b_kernel(
                               ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                               input_picture_ptr->stride_y,
                               ((uint16_t *)prediction_ptr->buffer_y) + cu_origin_index,
                               prediction_ptr->stride_y,
                               context_ptr->blk_geom->bheight,
                               context_ptr->blk_geom->bwidth));
        }
    }
    if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1 &&
        context_ptr->md_staging_skip_chroma_pred == FALSE) {
        if (context_ptr->mds0_ctrls.mds0_dist_type == MDS0_SSD) {
            EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
                ? svt_full_distortion_kernel16_bits
                : svt_spatial_full_distortion_kernel;
            chroma_fast_distortion                           = (uint32_t)spatial_full_dist_type_fun(
                input_picture_ptr->buffer_cb,
                input_cb_origin_in_index,
                input_picture_ptr->stride_cb,
                candidate_buffer->prediction_ptr->buffer_cb,
                (int32_t)cu_chroma_origin_index,
                prediction_ptr->stride_cb,
                context_ptr->blk_geom->bwidth_uv,
                context_ptr->blk_geom->bheight_uv);
            chroma_fast_distortion += (uint32_t)spatial_full_dist_type_fun(
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
                chroma_fast_distortion = svt_nxm_sad_kernel_sub_sampled(
                    input_picture_ptr->buffer_cb + input_cb_origin_in_index,
                    input_picture_ptr->stride_cb,
                    candidate_buffer->prediction_ptr->buffer_cb + cu_chroma_origin_index,
                    prediction_ptr->stride_cb,
                    context_ptr->blk_geom->bheight_uv,
                    context_ptr->blk_geom->bwidth_uv);

                chroma_fast_distortion += svt_nxm_sad_kernel_sub_sampled(
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
    }
    if (context_ptr->mds0_ctrls.enable_cost_based_early_exit &&
        context_ptr->mds0_best_cost != (uint32_t)~0) {
        const uint64_t distortion_cost = RDCOST(
            (context_ptr->mds0_ctrls.mds0_dist_type == MDS0_SSD) ? full_lambda : fast_lambda,
            0,
            luma_fast_distortion + chroma_fast_distortion);
        if (distortion_cost > context_ptr->mds0_best_cost &&
            (100 * (distortion_cost - context_ptr->mds0_best_cost)) >
                (context_ptr->mds0_best_cost * context_ptr->mds0_ctrls.mds0_distortion_th)) {
            *(candidate_buffer->fast_cost_ptr) = MAX_MODE_COST;
            return;
        }
    }
    // Fast Cost
    if (context_ptr->shut_fast_rate) {
        *(candidate_buffer->fast_cost_ptr) = luma_fast_distortion + chroma_fast_distortion;
        candidate_buffer->fast_luma_rate   = 0;
        candidate_buffer->fast_chroma_rate = 0;
    } else {
        *(candidate_buffer->fast_cost_ptr) =
            av1_product_fast_cost_func_table[is_inter_mode(candidate_ptr->pred_mode)](
                context_ptr,
                blk_ptr,
                candidate_buffer,
                NOT_USED_VALUE,
                luma_fast_distortion,
                chroma_fast_distortion,
                (context_ptr->mds0_ctrls.mds0_dist_type == MDS0_SSD) ? full_lambda : fast_lambda,
                pcs_ptr,
                &(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                      .ed_ref_mv_stack[candidate_ptr->ref_frame_type][0]),
                context_ptr->blk_geom,
                context_ptr->blk_origin_y >> MI_SIZE_LOG2,
                context_ptr->blk_origin_x >> MI_SIZE_LOG2,
                context_ptr->inter_intra_comp_ctrls.enabled,
                context_ptr->intra_luma_left_mode,
                context_ptr->intra_luma_top_mode);
    }
    // Init full cost in case we by pass stage1/stage2
    if (context_ptr->nic_ctrls.md_staging_mode == MD_STAGING_MODE_0)
        *(candidate_buffer->full_cost_ptr) = *(candidate_buffer->fast_cost_ptr);
}
void set_inter_comp_controls(ModeDecisionContext *ctx, uint8_t inter_comp_mode) {
    InterCompCtrls *inter_comp_ctrls = &ctx->inter_comp_ctrls;

    switch (inter_comp_mode) {
    case 0: //OFF
        inter_comp_ctrls->tot_comp_types = 1;
        break;
    case 1: //FULL
        inter_comp_ctrls->tot_comp_types      = (get_wedge_params_bits(ctx->blk_geom->bsize) == 0)
                 ? MD_COMP_WEDGE
                 : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 1;
        inter_comp_ctrls->do_pme              = 1;
        inter_comp_ctrls->do_nearest_near_new = 1;
        inter_comp_ctrls->do_3x3_bi           = 1;

        inter_comp_ctrls->skip_mvp_on_ref_info = 0;
        inter_comp_ctrls->use_rate             = 1;
        inter_comp_ctrls->pred0_to_pred1_mult  = 0;

        break;
    case 2:
        inter_comp_ctrls->tot_comp_types      = (get_wedge_params_bits(ctx->blk_geom->bsize) == 0)
                 ? MD_COMP_WEDGE
                 : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 1;
        inter_comp_ctrls->do_pme              = 1;
        inter_comp_ctrls->do_nearest_near_new = 1;
        inter_comp_ctrls->do_3x3_bi           = 0;

        inter_comp_ctrls->skip_mvp_on_ref_info = 0;
        inter_comp_ctrls->use_rate             = 0;
        inter_comp_ctrls->pred0_to_pred1_mult  = 0;

        break;

    case 3:
        inter_comp_ctrls->tot_comp_types      = (get_wedge_params_bits(ctx->blk_geom->bsize) == 0)
                 ? MD_COMP_WEDGE
                 : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 1;
        inter_comp_ctrls->do_pme              = 1;
        inter_comp_ctrls->do_nearest_near_new = 0;
        inter_comp_ctrls->do_3x3_bi           = 0;

        inter_comp_ctrls->skip_mvp_on_ref_info = 0;
        inter_comp_ctrls->use_rate             = 0;
        inter_comp_ctrls->pred0_to_pred1_mult  = 1;

        break;
    case 4:
        inter_comp_ctrls->tot_comp_types      = (get_wedge_params_bits(ctx->blk_geom->bsize) == 0)
                 ? MD_COMP_WEDGE
                 : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 0;
        inter_comp_ctrls->do_pme              = 0;
        inter_comp_ctrls->do_nearest_near_new = 0;
        inter_comp_ctrls->do_3x3_bi           = 0;

        inter_comp_ctrls->skip_mvp_on_ref_info = 1;
        inter_comp_ctrls->use_rate             = 0;
        inter_comp_ctrls->pred0_to_pred1_mult  = 4;

        break;
    case 5:
        inter_comp_ctrls->tot_comp_types      = (get_wedge_params_bits(ctx->blk_geom->bsize) == 0)
                 ? MD_COMP_WEDGE
                 : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 0;
        inter_comp_ctrls->do_me               = 0;
        inter_comp_ctrls->do_pme              = 0;
        inter_comp_ctrls->do_nearest_near_new = 0;
        inter_comp_ctrls->do_3x3_bi           = 0;

        inter_comp_ctrls->skip_mvp_on_ref_info = 1;
        inter_comp_ctrls->use_rate             = 0;
        inter_comp_ctrls->pred0_to_pred1_mult  = 2;

        break;
    default: assert(0); break;
    }
}
/* Set the max number of NICs for each MD stage, based on the picture type and scaling settings.

   pic_type = I_SLICE ? 0 : REF ? 1 : 2;
*/
void set_nics(NicScalingCtrls *scaling_ctrls, uint32_t mds1_count[CAND_CLASS_TOTAL],
              uint32_t mds2_count[CAND_CLASS_TOTAL], uint32_t mds3_count[CAND_CLASS_TOTAL],
              uint8_t pic_type) {
    for (CandClass cidx = CAND_CLASS_0; cidx < CAND_CLASS_TOTAL; cidx++) {
        mds1_count[cidx] = MD_STAGE_NICS[pic_type][cidx];
        mds2_count[cidx] = MD_STAGE_NICS[pic_type][cidx] >> 1;
        mds3_count[cidx] = MD_STAGE_NICS[pic_type][cidx] >> 2;
    }

    // minimum nics allowed
    uint32_t min_mds1_nics = (pic_type < 2 && scaling_ctrls->stage1_scaling_num) ? 2 : 1;
    uint32_t min_mds2_nics = (pic_type < 2 && scaling_ctrls->stage2_scaling_num) ? 2 : 1;
    uint32_t min_mds3_nics = (pic_type < 2 && scaling_ctrls->stage3_scaling_num) ? 2 : 1;

    // Set the scaling numerators
    uint32_t stage1_num = scaling_ctrls->stage1_scaling_num;
    uint32_t stage2_num = scaling_ctrls->stage2_scaling_num;
    uint32_t stage3_num = scaling_ctrls->stage3_scaling_num;
    // The scaling denominator is 16 for all stages
    uint32_t scale_denum = MD_STAGE_NICS_SCAL_DENUM;
    // no NIC setting should be done beyond this point
    for (CandClass cidx = 0; cidx < CAND_CLASS_TOTAL; ++cidx) {
        mds1_count[cidx] = MAX(min_mds1_nics,
                               DIVIDE_AND_ROUND(mds1_count[cidx] * stage1_num, scale_denum));
        mds2_count[cidx] = MAX(min_mds2_nics,
                               DIVIDE_AND_ROUND(mds2_count[cidx] * stage2_num, scale_denum));
        mds3_count[cidx] = MAX(min_mds3_nics,
                               DIVIDE_AND_ROUND(mds3_count[cidx] * stage3_num, scale_denum));
    }
}

void set_md_stage_counts(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    // Step 1: Set the number of NICs for each stage
    // no NIC setting should be done beyond this point
    // Set md_stage count
    uint8_t pic_type = pcs->slice_type == I_SLICE        ? 0
        : pcs->parent_pcs_ptr->is_used_as_reference_flag ? 1
                                                         : 2;
    set_nics(&ctx->nic_ctrls.scaling_ctrls,
             ctx->md_stage_1_count,
             ctx->md_stage_2_count,
             ctx->md_stage_3_count,
             pic_type);

    // Step 2: derive bypass_stage1 flags
    ctx->bypass_md_stage_1 = (ctx->nic_ctrls.md_staging_mode == MD_STAGING_MODE_1 ||
                              ctx->nic_ctrls.md_staging_mode == MD_STAGING_MODE_2)
        ? FALSE
        : TRUE;
    ctx->bypass_md_stage_2 = (ctx->nic_ctrls.md_staging_mode == MD_STAGING_MODE_2) ? FALSE : TRUE;

    // Step 3: update count for md_stage_1 and d_stage_2 if bypassed (no NIC
    // setting should be done beyond this point)
    if (ctx->bypass_md_stage_1) {
        for (CandClass cidx = CAND_CLASS_0; cidx < CAND_CLASS_TOTAL; cidx++)
            ctx->md_stage_2_count[cidx] = ctx->md_stage_1_count[cidx];
    }
    if (ctx->bypass_md_stage_2) {
        for (CandClass cidx = CAND_CLASS_0; cidx < CAND_CLASS_TOTAL; cidx++)
            ctx->md_stage_3_count[cidx] = ctx->md_stage_2_count[cidx];
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
void construct_best_sorted_arrays_md_stage_3(
    struct ModeDecisionContext *context_ptr,
    uint32_t                   *best_candidate_index_array) { //best = union from all classes

    uint32_t best_candi = 0;
    for (CandClass class_i = CAND_CLASS_0; class_i < CAND_CLASS_TOTAL; class_i++)
        for (uint32_t candi = 0; candi < context_ptr->md_stage_3_count[class_i]; candi++)
            best_candidate_index_array[best_candi++] =
                context_ptr->cand_buff_indices[class_i][candi];

    assert(best_candi == context_ptr->md_stage_3_total_count);
}
void get_mds3_intra_count_for_chroma(struct ModeDecisionContext   *context_ptr,
                                     ModeDecisionCandidateBuffer **buffer_ptr_array,
                                     uint32_t                     *best_candidate_index_array) {
    uint32_t fullReconCandidateCount = context_ptr->md_stage_3_total_count;

    // Only if chroma_at_last_md_stage
    uint32_t i, id;
    context_ptr->md_stage_3_total_intra_count = 0;
    for (i = 0; i < fullReconCandidateCount; ++i) {
        id                  = best_candidate_index_array[i];
        const Bool is_inter = (is_inter_mode(buffer_ptr_array[id]->candidate_ptr->pred_mode) ||
                               buffer_ptr_array[id]->candidate_ptr->use_intrabc)
            ? TRUE
            : FALSE;
        context_ptr->md_stage_3_total_intra_count += !is_inter ? 1 : 0;
    }

    // Derive best_intra_cost and best_inter_cost
    context_ptr->best_intra_cost = MAX_MODE_COST;
    context_ptr->best_inter_cost = MAX_MODE_COST;
    for (i = 0; i < fullReconCandidateCount; ++i) {
        id                  = best_candidate_index_array[i];
        const Bool is_inter = (is_inter_mode(buffer_ptr_array[id]->candidate_ptr->pred_mode) ||
                               buffer_ptr_array[id]->candidate_ptr->use_intrabc)
            ? TRUE
            : FALSE;
        if (!is_inter)
            if (*buffer_ptr_array[id]->full_cost_ptr < context_ptr->best_intra_cost)
                context_ptr->best_intra_cost = *buffer_ptr_array[id]->full_cost_ptr;
        if (is_inter)
            if (*buffer_ptr_array[id]->full_cost_ptr < context_ptr->best_inter_cost)
                context_ptr->best_inter_cost = *buffer_ptr_array[id]->full_cost_ptr;
    }

    // Update md_stage_3_total_intra_count based based on inter/intra cost deviation
    if ((context_ptr->best_inter_cost * context_ptr->uv_ctrls.uv_intra_th) <
        (context_ptr->best_intra_cost * 100))
        context_ptr->md_stage_3_total_intra_count = 0;
}
/*
*/
void md_stage_0_light_pd0(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                          int32_t fast_candidate_end_index, EbPictureBufferDesc *input_picture_ptr,
                          uint32_t input_origin_index, uint32_t blk_origin_index) {
    int32_t fast_loop_cand_index = fast_candidate_end_index;
    while (fast_loop_cand_index >= 0) {
        ModeDecisionCandidateBuffer *candidate_buffer =
            context_ptr->candidate_buffer_ptr_array[fast_loop_cand_index];
        candidate_buffer->candidate_ptr = &context_ptr->fast_candidate_array[fast_loop_cand_index];

        // Initialize tx_depth
        candidate_buffer->candidate_ptr->tx_depth = 0;
        fast_loop_core_light_pd0(candidate_buffer,
                                 pcs_ptr,
                                 context_ptr,
                                 input_picture_ptr,
                                 input_origin_index,
                                 blk_origin_index);

        if (*candidate_buffer->fast_cost_ptr < context_ptr->mds0_best_cost) {
            context_ptr->mds0_best_cost = *candidate_buffer->fast_cost_ptr;
            context_ptr->mds0_best_idx  = fast_loop_cand_index;
        }
        --fast_loop_cand_index;
    }
}
void md_stage_0_light_pd1(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                          ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base,
                          ModeDecisionCandidate *fast_candidate_array, int32_t fast_cand_idx,
                          EbPictureBufferDesc *input_picture_ptr, BlockLocation *loc,
                          BlkStruct *blk_ptr) {
    // Set MD Staging fast_loop_core settings
    context_ptr->md_staging_skip_chroma_pred = TRUE;
    context_ptr->end_plane                   = 1;

    // 2nd fast loop: src-to-recon
    uint32_t cand_buff_idx = 0;

    while (fast_cand_idx >= 0) {
        ModeDecisionCandidateBuffer *candidate_buffer =
            candidate_buffer_ptr_array_base[cand_buff_idx];
        candidate_buffer->candidate_ptr = &fast_candidate_array[fast_cand_idx];

        // Initialize tx_depth
        candidate_buffer->candidate_ptr->tx_depth = 0;

        // Prediction
        fast_loop_core_light_pd1(
            candidate_buffer, pcs_ptr, context_ptr, input_picture_ptr, loc, blk_ptr);

        if (*candidate_buffer->fast_cost_ptr < context_ptr->mds0_best_cost) {
            context_ptr->mds0_best_cost = *candidate_buffer->fast_cost_ptr;
            context_ptr->mds0_best_idx  = cand_buff_idx;
            cand_buff_idx               = cand_buff_idx ? 0 : 1;
        }

        --fast_cand_idx;
    }
}
void md_stage_0(

    PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base,
    ModeDecisionCandidate *fast_candidate_array, int32_t fast_candidate_start_index,
    int32_t fast_candidate_end_index, EbPictureBufferDesc *input_picture_ptr,
    uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
    uint32_t input_cr_origin_in_index, BlkStruct *blk_ptr, uint32_t blk_origin_index,
    uint32_t blk_chroma_origin_index, uint32_t candidate_buffer_start_index, uint32_t max_buffers,
    Bool scratch_buffer_pesent_flag) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    int32_t             fast_loop_cand_index;
    uint32_t            highest_cost_index;
    // Set MD Staging fast_loop_core settings
    context_ptr->md_staging_skip_interpolation_search = (context_ptr->ifs_ctrls.level == IFS_MDS0)
        ? FALSE
        : TRUE;
    context_ptr->md_staging_skip_chroma_pred          = TRUE;
    // 2nd fast loop: src-to-recon
    highest_cost_index     = candidate_buffer_start_index;
    fast_loop_cand_index   = fast_candidate_end_index;
    context_ptr->end_plane = (context_ptr->blk_geom->has_uv &&
                              context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1 &&
                              !context_ptr->md_staging_skip_chroma_pred)
        ? (int)MAX_MB_PLANE
        : 1;
    /* Variables for sorting the scratch buffer and getting the next index. */
    // max_cost is volatile to prevent the compiler from loading 0xFFFFFFFFFFFFFF
    //   as a const at the early-out. Loading a large constant on intel x64 processors
    //   clogs the i-cache/intstruction decode. This still reloads the variable from
    //   the stack each pass, so a better solution would be to register the variable,
    //   but this might require asm.
    volatile uint64_t max_cost = MAX_CU_COST;

    // Process CLASS_0 candidates through iterations at mds0;
    // 1st iteration : evaluate the regular mode(s).
    //
    // 2nd iteration : evaluate the angular (pAngle!=0) mode(s)
    //     skip pAngle !=0 of a given mode, if pAngle==0 of the same mode is not the best after the 1st iteration
    //     and if the cost to the best is significant
    //
    // 3rd iteration : evaluate the filter - intra mode(s)
    //     if the winner is DC/H/V/D157/Paeth, then only test filter-DC/filter-H/filter-V/filter-D157/filter-Paeth (i.e.at most 1 out of the 5 possible filter-intra)
    //     else test all filter - intra.

    uint8_t tot_itr = (context_ptr->target_class != CAND_CLASS_0 ||
                       !context_ptr->cand_reduction_ctrls.mds0_reduce_intra)
        ? 1
        : (context_ptr->md_filter_intra_level) ? 3
                                               : 2;

    uint64_t       best_reg_intra_cost                = MAX_CU_COST; // Derived at the 1st itr
    PredictionMode best_reg_intra_mode                = INTRA_INVALID; // Derived at the 1st itr
    uint64_t       regular_intra_cost[PAETH_PRED + 1] = {MAX_CU_COST};

    uint32_t tot_processed_cand = 0;

    for (uint8_t itr = 0; itr < tot_itr; itr++) {
        for (fast_loop_cand_index = fast_candidate_end_index;
             fast_loop_cand_index >= fast_candidate_start_index;
             fast_loop_cand_index--) {
            if (fast_candidate_array[fast_loop_cand_index].cand_class ==
                context_ptr->target_class) {
                ModeDecisionCandidateBuffer *candidate_buffer =
                    candidate_buffer_ptr_array_base[highest_cost_index];
                candidate_buffer->candidate_ptr = &fast_candidate_array[fast_loop_cand_index];
                // Initialize tx_depth
                candidate_buffer->candidate_ptr->tx_depth       = 0;
                candidate_buffer->candidate_ptr->interp_filters = 0;
                if (tot_itr > 1) {
                    if (itr == 0) {
                        // Eval regular only if itr=0 (i.e. skip angular and skip filter)
                        if (fast_candidate_array[fast_loop_cand_index].angle_delta[PLANE_TYPE_Y] !=
                                0 ||
                            fast_candidate_array[fast_loop_cand_index].filter_intra_mode !=
                                FILTER_INTRA_MODES)
                            continue;
                    } else if (itr == 1) {
                        // Eval angular only if itr=1 (i.e. skip regular and skip filter)
                        if (fast_candidate_array[fast_loop_cand_index].angle_delta[PLANE_TYPE_Y] ==
                                0 ||
                            fast_candidate_array[fast_loop_cand_index].filter_intra_mode !=
                                FILTER_INTRA_MODES)
                            continue;
                        // Use regular info to reduce angular processing
                        // Eval the child-angular if the parent-angular is the best out of the regular modes (i.e. itr=0 eval)
                        else {
                            if (best_reg_intra_mode !=
                                fast_candidate_array[fast_loop_cand_index].pred_mode) {
                                if ((((regular_intra_cost[fast_candidate_array[fast_loop_cand_index]
                                                              .pred_mode] -
                                       MAX(best_reg_intra_cost, 1)) *
                                      100) /
                                     MAX(best_reg_intra_cost, 1)) > MDS0_REDUCE_ANGULAR_INTRA_TH)
                                    continue;
                            }
                        }
                    } else {
                        // Eval filter only if itr=1 (i.e. skip regular and skip angular)
                        // Eval regular only if itr=0 (i.e. skip angular and skip filter)
                        if (fast_candidate_array[fast_loop_cand_index].filter_intra_mode ==
                            FILTER_INTRA_MODES)
                            continue;

                        if (best_reg_intra_mode == PAETH_PRED) {
                            if (fast_candidate_array[fast_loop_cand_index].filter_intra_mode !=
                                FILTER_PAETH_PRED)
                                continue;
                        } else if (best_reg_intra_mode == D157_PRED) {
                            if (fast_candidate_array[fast_loop_cand_index].filter_intra_mode !=
                                FILTER_D157_PRED)
                                continue;
                        } else if (best_reg_intra_mode == H_PRED) {
                            if (fast_candidate_array[fast_loop_cand_index].filter_intra_mode !=
                                FILTER_H_PRED)
                                continue;
                        } else if (best_reg_intra_mode == V_PRED) {
                            if (fast_candidate_array[fast_loop_cand_index].filter_intra_mode !=
                                FILTER_V_PRED)
                                continue;
                        } else if (best_reg_intra_mode == DC_PRED) {
                            if (fast_candidate_array[fast_loop_cand_index].filter_intra_mode !=
                                FILTER_DC_PRED)
                                continue;
                        }
                    }
                }
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
                               blk_chroma_origin_index);

                tot_processed_cand++;

                if (scs_ptr->vq_ctrls.sharpness_ctrls.unipred_bias &&
                    pcs_ptr->parent_pcs_ptr->is_noise_level &&
                    is_inter_singleref_mode(candidate_buffer->candidate_ptr->pred_mode)) {
                    *candidate_buffer->fast_cost_ptr = (*candidate_buffer->fast_cost_ptr *
                                                        uni_psy_bias[pcs_ptr->picture_qp]) /
                        100;
                }
                if (*candidate_buffer->fast_cost_ptr < context_ptr->mds0_best_cost) {
                    context_ptr->mds0_best_cost = *candidate_buffer->fast_cost_ptr;
                    context_ptr->mds0_best_class =
                        fast_candidate_array[fast_loop_cand_index].cand_class;
                }
                if (tot_itr > 1) {
                    if (itr == 0) {
                        regular_intra_cost[fast_candidate_array[fast_loop_cand_index].pred_mode] =
                            *candidate_buffer->fast_cost_ptr;

                        if (*candidate_buffer->fast_cost_ptr < best_reg_intra_cost) {
                            best_reg_intra_cost = *candidate_buffer->fast_cost_ptr;
                            best_reg_intra_mode =
                                fast_candidate_array[fast_loop_cand_index].pred_mode;
                        }
                    }
                }
                // Find the buffer with the highest cost
                if (fast_loop_cand_index || scratch_buffer_pesent_flag) {
                    const uint64_t *fast_cost_array    = context_ptr->fast_cost_array;
                    const uint32_t  buffer_index_start = candidate_buffer_start_index;
                    const uint32_t  buffer_index_end   = buffer_index_start + max_buffers;
                    if (buffer_index_end == 2) {
                        highest_cost_index = fast_cost_array[0] < fast_cost_array[1] ? 1 : 0;
                    } else {
                        highest_cost_index    = buffer_index_start;
                        uint32_t buffer_index = buffer_index_start + 1;
                        do {
                            uint64_t highest_cost = fast_cost_array[highest_cost_index];
                            if (highest_cost == max_cost)
                                break;

                            if (fast_cost_array[buffer_index] > highest_cost)
                                highest_cost_index = buffer_index;
                        } while (++buffer_index < buffer_index_end);
                    }
                }
            }
        }
    }

    //if pruning happened, update MDS1 count accordingly to not process invalid candidates in subsequent MD stages
    context_ptr->md_stage_1_count[context_ptr->target_class] = MIN(
        context_ptr->md_stage_1_count[context_ptr->target_class], tot_processed_cand);

    // Set the cost of the scratch canidate to max to get discarded @ the sorting phase
    *(candidate_buffer_ptr_array_base[highest_cost_index]->fast_cost_ptr) =
        (scratch_buffer_pesent_flag)
        ? MAX_CU_COST
        : *(candidate_buffer_ptr_array_base[highest_cost_index]->fast_cost_ptr);
}
void svt_pme_sad_loop_kernel_c(const struct svt_mv_cost_param *mv_cost_params,
                               uint8_t  *src, // input parameter, source samples Ptr
                               uint32_t  src_stride, // input parameter, source stride
                               uint8_t  *ref, // input parameter, reference samples Ptr
                               uint32_t  ref_stride, // input parameter, reference stride
                               uint32_t  block_height, // input parameter, block height (M)
                               uint32_t  block_width, // input parameter, block width (N)
                               uint32_t *best_cost, int16_t *best_mvx, int16_t *best_mvy,
                               int16_t search_position_start_x, int16_t search_position_start_y,
                               int16_t search_area_width, int16_t search_area_height,
                               int16_t search_step, int16_t mvx, int16_t mvy) {
    int16_t xSearchIndex;
    int16_t ySearchIndex;
    int16_t col_num       = 0;
    int16_t search_step_x = 1;
    for (ySearchIndex = 0; ySearchIndex < search_area_height; ySearchIndex += search_step) {
        for (xSearchIndex = 0; xSearchIndex < search_area_width; xSearchIndex += search_step_x) {
            if (((search_area_width - xSearchIndex) < 8) && (col_num == 0))
                continue;
            if (col_num == 7) {
                col_num       = 0;
                search_step_x = search_step;
            } else {
                col_num++;
                search_step_x = 1;
            }
            uint32_t x, y;
            uint32_t cost = 0;

            for (y = 0; y < block_height; y++) {
                for (x = 0; x < block_width; x++)
                    cost += EB_ABS_DIFF(src[y * src_stride + x],
                                        ref[xSearchIndex + y * ref_stride + x]);
            }

            MV       best_mv;
            uint32_t refinement_pos_x = search_position_start_x + xSearchIndex;
            uint32_t refinement_pos_y = search_position_start_y + ySearchIndex;
            best_mv.col               = mvx + (refinement_pos_x * 8);
            best_mv.row               = mvy + (refinement_pos_y * 8);
            cost += fp_mv_err_cost(&best_mv, mv_cost_params);
            if (cost < *best_cost) {
                *best_mvx  = mvx + (refinement_pos_x * 8);
                *best_mvy  = mvy + (refinement_pos_y * 8);
                *best_cost = cost;
            }
        }

        ref += search_step * ref_stride;
    }

    return;
}

void md_full_pel_search_large_lbd(MV_COST_PARAMS *mv_cost_params, ModeDecisionContext *context_ptr,
                                  EbPictureBufferDesc *input_picture_ptr,
                                  EbPictureBufferDesc *ref_pic, uint32_t input_origin_index,
                                  int16_t mvx, int16_t mvy, int16_t search_position_start_x,
                                  int16_t search_position_end_x, int16_t search_position_start_y,
                                  int16_t search_position_end_y, int16_t sparse_search_step,
                                  int16_t *best_mvx, int16_t *best_mvy, uint32_t *best_cost) {
    //We cannot use sparse_search_step with mpsad for search_position_start_x/search_position_end_x,
    //So for x dimension we assume sparse_search_step is always 1

    int32_t ref_origin_index = ref_pic->origin_x +
        (context_ptr->blk_origin_x + (mvx >> 3) + search_position_start_x) +
        (context_ptr->blk_origin_y + (mvy >> 3) + ref_pic->origin_y + search_position_start_y) *
            ref_pic->stride_y;

    int16_t remain_search_area = 8 - ((search_position_end_x - search_position_start_x) % 8);
    remain_search_area         = remain_search_area == 8 ? 0 : remain_search_area;
    search_position_end_x = MAX(search_position_end_x, search_position_end_x + remain_search_area);
    uint32_t search_area_width  = search_position_end_x - search_position_start_x;
    uint32_t search_area_height = search_position_end_y - search_position_start_y + 1;
    assert(!(search_area_width & 7));
    if (search_area_width & 0xfffffff8) {
        svt_pme_sad_loop_kernel(
            mv_cost_params,
            input_picture_ptr->buffer_y + input_origin_index,
            input_picture_ptr->stride_y,
            ref_pic->buffer_y + ref_origin_index,
            ref_pic->stride_y,
            context_ptr->blk_geom->bheight,
            context_ptr->blk_geom->bwidth,
            best_cost,
            best_mvx,
            best_mvy,
            search_position_start_x,
            search_position_start_y,
            (search_area_width & 0xfffffff8), //pass search_area_width multiple by 8
            search_area_height,
            sparse_search_step,
            mvx,
            mvy);
    }

    if (search_area_width & 7) {
        uint32_t cost;
        printf("Error 1\n");
        for (int32_t refinement_pos_y = search_position_start_y;
             refinement_pos_y <= search_position_end_y;
             refinement_pos_y = refinement_pos_y + sparse_search_step) {
            int32_t refinement_pos_x = search_position_start_x + (search_area_width & 0xfffffff8);
            for (; refinement_pos_x <= search_position_end_x; refinement_pos_x++) {
                ref_origin_index = ref_pic->origin_x +
                    (context_ptr->blk_origin_x + (mvx >> 3) + refinement_pos_x) +
                    (context_ptr->blk_origin_y + (mvy >> 3) + ref_pic->origin_y +
                     refinement_pos_y) *
                        ref_pic->stride_y;

                assert((context_ptr->blk_geom->bwidth >> 3) < 17);

                cost = svt_nxm_sad_kernel_sub_sampled(
                    input_picture_ptr->buffer_y + input_origin_index,
                    input_picture_ptr->stride_y,
                    ref_pic->buffer_y + ref_origin_index,
                    ref_pic->stride_y,
                    context_ptr->blk_geom->bheight,
                    context_ptr->blk_geom->bwidth);

                MV best_mv;
                best_mv.col = mvx + (refinement_pos_x * 8);
                best_mv.row = mvy + (refinement_pos_y * 8);
                cost += fp_mv_err_cost(&best_mv, mv_cost_params);
                if (cost < *best_cost) {
                    *best_mvx  = mvx + (refinement_pos_x * 8);
                    *best_mvy  = mvy + (refinement_pos_y * 8);
                    *best_cost = cost;
                }
            }
        }
    }
}

void md_full_pel_search(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                        EbPictureBufferDesc *input_picture_ptr, EbPictureBufferDesc *ref_pic,
                        uint32_t input_origin_index, Bool use_ssd, int16_t mvx, int16_t mvy,
                        int16_t search_position_start_x, int16_t search_position_end_x,
                        int16_t search_position_start_y, int16_t search_position_end_y,
                        int16_t sparse_search_step, uint8_t is_sprs_lev0_performed,
                        int16_t *best_mvx, int16_t *best_mvy, uint32_t *best_cost,
                        uint8_t hbd_mode_decision) {
    // Mvcost params
    MV_COST_PARAMS mv_cost_params;
    FrameHeader   *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;

    uint32_t rdmult = use_ssd
        ? context_ptr->full_lambda_md[hbd_mode_decision ? EB_10_BIT_MD : EB_8_BIT_MD]
        : context_ptr->fast_lambda_md[hbd_mode_decision ? EB_10_BIT_MD : EB_8_BIT_MD];
    svt_init_mv_cost_params(&mv_cost_params,
                            context_ptr,
                            &context_ptr->ref_mv,
                            frm_hdr->quantization_params.base_q_idx,
                            rdmult,
                            hbd_mode_decision);
    uint32_t cost;
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

    if (context_ptr->enable_psad) {
        if (!use_ssd && !hbd_mode_decision &&
            (search_position_end_x - search_position_start_x) >= 7) {
            md_full_pel_search_large_lbd(&mv_cost_params,
                                         context_ptr,
                                         input_picture_ptr,
                                         ref_pic,
                                         input_origin_index,
                                         mvx,
                                         mvy,
                                         search_position_start_x,
                                         search_position_end_x,
                                         search_position_start_y,
                                         search_position_end_y,
                                         sparse_search_step,
                                         best_mvx,
                                         best_mvy,
                                         best_cost);
            return;
        }
    }
    for (int32_t refinement_pos_x = search_position_start_x;
         refinement_pos_x <= search_position_end_x;
         refinement_pos_x = refinement_pos_x + sparse_search_step) {
        for (int32_t refinement_pos_y = search_position_start_y;
             refinement_pos_y <= search_position_end_y;
             refinement_pos_y = refinement_pos_y + sparse_search_step) {
            // If sparse search level_1, and if search level_0 previously performed
            if (sparse_search_step == 2 && is_sprs_lev0_performed) {
                // If level_0 range
                if ((refinement_pos_x + (mvx >> 3)) >= context_ptr->sprs_lev0_start_x &&
                    (refinement_pos_x + (mvx >> 3)) <= context_ptr->sprs_lev0_end_x &&
                    (refinement_pos_y + (mvy >> 3)) >= context_ptr->sprs_lev0_start_y &&
                    (refinement_pos_y + (mvy >> 3)) <= context_ptr->sprs_lev0_end_y)
                    // If level_0 position
                    if (refinement_pos_x % 4 == 0 && refinement_pos_y % 4 == 0)
                        continue;
            }
            int32_t ref_origin_index = ref_pic->origin_x +
                (context_ptr->blk_origin_x + (mvx >> 3) + refinement_pos_x) +
                (context_ptr->blk_origin_y + (mvy >> 3) + ref_pic->origin_y + refinement_pos_y) *
                    ref_pic->stride_y;
            if (use_ssd) {
                EbSpatialFullDistType spatial_full_dist_type_fun = hbd_mode_decision
                    ? svt_full_distortion_kernel16_bits
                    : svt_spatial_full_distortion_kernel;

                cost = (uint32_t)spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
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
                    cost = sad_16b_kernel(
                        ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                        input_picture_ptr->stride_y,
                        ((uint16_t *)ref_pic->buffer_y) + ref_origin_index,
                        ref_pic->stride_y,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                } else {
                    cost = svt_nxm_sad_kernel_sub_sampled(
                        input_picture_ptr->buffer_y + input_origin_index,
                        input_picture_ptr->stride_y,
                        ref_pic->buffer_y + ref_origin_index,
                        ref_pic->stride_y,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                }
            }
            MV best_mv;
            best_mv.col = mvx + (refinement_pos_x * 8);
            best_mv.row = mvy + (refinement_pos_y * 8);
            cost += fp_mv_err_cost(&best_mv, &mv_cost_params);
            if (cost < *best_cost) {
                *best_mvx  = mvx + (refinement_pos_x * 8);
                *best_mvy  = mvy + (refinement_pos_y * 8);
                *best_cost = cost;
            }
        }
    }
}
uint8_t get_max_drl_index(uint8_t refmvCnt, PredictionMode mode);
uint8_t is_me_data_present(struct ModeDecisionContext *context_ptr, const MeSbResults *me_results,
                           uint8_t list_idx, uint8_t ref_idx);
// Derive me_sb_addr and me_block_offset used to access ME_MV
void derive_me_offsets(const SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                       ModeDecisionContext *context_ptr) {
    // @ this stage NSQ block(s) are inheriting SQ block(s) ME results; MV(s), pruning PA_ME results

    // Get parent_depth_idx_mds
    uint16_t         parent_depth_idx_mds = context_ptr->blk_geom->parent_depth_idx_mds;
    const BlockGeom *sq_blk_geom = (context_ptr->blk_geom->bwidth != context_ptr->blk_geom->bheight)
        ? get_blk_geom_mds(context_ptr->blk_geom
                               ->sqi_mds) // Use parent block SQ info as ME not performed for NSQ
        : (context_ptr->blk_geom->bwidth == 4 ||
           context_ptr->blk_geom->bheight ==
               4) // Use parent_depth SQ block info as ME not performed for 4x4
        ? get_blk_geom_mds(parent_depth_idx_mds)
        : context_ptr->blk_geom;

    context_ptr->geom_offset_x = 0;
    context_ptr->geom_offset_y = 0;

    if (scs_ptr->seq_header.sb_size == BLOCK_128X128) {
        uint32_t me_sb_size         = scs_ptr->b64_size;
        uint32_t me_pic_width_in_sb = (pcs_ptr->parent_pcs_ptr->aligned_width + scs_ptr->b64_size -
                                       1) /
            me_sb_size;
        uint32_t me_sb_x             = (context_ptr->blk_origin_x / me_sb_size);
        uint32_t me_sb_y             = (context_ptr->blk_origin_y / me_sb_size);
        context_ptr->me_sb_addr      = me_sb_x + me_sb_y * me_pic_width_in_sb;
        context_ptr->geom_offset_x   = (me_sb_x & 0x1) * me_sb_size;
        context_ptr->geom_offset_y   = (me_sb_y & 0x1) * me_sb_size;
        context_ptr->me_block_offset = (uint32_t)
            me_idx_128x128[((context_ptr->geom_offset_y / me_sb_size) * 2) +
                           (context_ptr->geom_offset_x / me_sb_size)]
                          [context_ptr->blk_geom->blkidx_mds];
    } else {
        context_ptr->me_sb_addr = context_ptr->sb_ptr->index;

        if (context_ptr->blk_geom->geom_idx == GEOM_0) {
            context_ptr->me_block_offset = me_idx_85[context_ptr->blk_geom->blkidx_mds];
            if (!context_ptr->sb_ptr->pcs_ptr->parent_pcs_ptr->enable_me_8x8) {
                if (context_ptr->me_block_offset >= MAX_SB64_PU_COUNT_NO_8X8)
                    context_ptr->me_block_offset =
                        me_idx_85_8x8_to_16x16_conversion[context_ptr->me_block_offset -
                                                          MAX_SB64_PU_COUNT_NO_8X8];
                if (!context_ptr->sb_ptr->pcs_ptr->parent_pcs_ptr->enable_me_16x16)
                    if (context_ptr->me_block_offset >= MAX_SB64_PU_COUNT_WO_16X16) {
                        assert(context_ptr->me_block_offset < 21);
                        context_ptr->me_block_offset =
                            me_idx_16x16_to_parent_32x32_conversion[context_ptr->me_block_offset -
                                                                    MAX_SB64_PU_COUNT_WO_16X16];
                    }
            }
        } else
            context_ptr->me_block_offset = me_idx[context_ptr->blk_geom->blkidx_mds];
    }

    if (sq_blk_geom->bwidth == 128 || sq_blk_geom->bheight == 128) {
        context_ptr->me_block_offset = 0;
    }
    assert(context_ptr->me_block_offset != (uint32_t)(-1));
    context_ptr->me_cand_offset = context_ptr->me_block_offset *
        pcs_ptr->parent_pcs_ptr->pa_me_data->max_cand;
}
#define MAX_MD_NSQ_SARCH_MVC_CNT 5
void md_nsq_motion_search(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                          EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index,
                          uint8_t list_idx, uint8_t ref_idx, const MeSbResults *me_results,
                          int16_t *me_mv_x, int16_t *me_mv_y) {
    // Step 0: derive the MVC list for the NSQ search; 1 SQ MV (default MV for NSQ) and up to 4 sub-block MV(s) (e.g. if 16x8 then 2 8x8, if 32x8 then 4 8x8)
    int16_t       mvc_x_array[MAX_MD_NSQ_SARCH_MVC_CNT];
    int16_t       mvc_y_array[MAX_MD_NSQ_SARCH_MVC_CNT];
    int8_t        mvc_count  = 0;
    const uint8_t max_refs   = pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs;
    const uint8_t max_l0     = pcs_ptr->parent_pcs_ptr->pa_me_data->max_l0;
    context_ptr->enable_psad = context_ptr->md_nsq_motion_search_ctrls.enable_psad;
    // SQ MV (default MVC for NSQ)
    mvc_x_array[mvc_count] = *me_mv_x;
    mvc_y_array[mvc_count] = *me_mv_y;
    mvc_count++;
    if ((context_ptr->blk_geom->bwidth != 4 && context_ptr->blk_geom->bheight != 4) &&
        context_ptr->blk_geom->sq_size >= 16) {
        uint8_t min_size = MIN(context_ptr->blk_geom->bwidth, context_ptr->blk_geom->bheight);
        // Derive the sub-block(s) MVs (additional MVC for NSQ)
        for (uint32_t block_index = 0;
             block_index < pcs_ptr->parent_pcs_ptr->max_number_of_pus_per_sb;
             block_index++) {
            if ((min_size == partition_width[block_index] ||
                 min_size == partition_height[block_index]) &&
                ((pu_search_index_map[block_index][0] >=
                  (context_ptr->blk_geom->origin_x - context_ptr->geom_offset_x)) &&
                 (pu_search_index_map[block_index][0] < context_ptr->blk_geom->bwidth +
                      (context_ptr->blk_geom->origin_x - context_ptr->geom_offset_x))) &&
                ((pu_search_index_map[block_index][1] >=
                  (context_ptr->blk_geom->origin_y - context_ptr->geom_offset_y)) &&
                 (pu_search_index_map[block_index][1] < context_ptr->blk_geom->bheight +
                      (context_ptr->blk_geom->origin_y - context_ptr->geom_offset_y)))) {
                if (list_idx == 0) {
                    mvc_x_array[mvc_count] =
                        (me_results->me_mv_array[block_index * max_refs + ref_idx].x_mv) << 1;
                    mvc_y_array[mvc_count] =
                        (me_results->me_mv_array[block_index * max_refs + ref_idx].y_mv) << 1;
                } else {
                    mvc_x_array[mvc_count] =
                        (me_results->me_mv_array[block_index * max_refs + max_l0 + ref_idx].x_mv)
                        << 1;
                    mvc_y_array[mvc_count] =
                        (me_results->me_mv_array[block_index * max_refs + max_l0 + ref_idx].y_mv)
                        << 1;
                }

                mvc_count++;
            }
        }
    }

    // Search Center
    int16_t              search_center_mvx  = mvc_x_array[0];
    int16_t              search_center_mvy  = mvc_y_array[0];
    uint32_t             search_center_cost = (uint32_t)~0;
    uint8_t              hbd_mode_decision  = EB_8_BIT_MD;
    EbReferenceObject   *ref_obj = pcs_ptr->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
    EbPictureBufferDesc *ref_pic = get_ref_pic_buffer(
        pcs_ptr, hbd_mode_decision, list_idx, ref_idx);
    // -------
    // Use scaled references if resolution of the reference is different from that of the input
    // -------
    use_scaled_rec_refs_if_needed(pcs_ptr, input_picture_ptr, ref_obj, &ref_pic, hbd_mode_decision);

    for (int16_t mvc_index = 0; mvc_index < mvc_count; mvc_index++) {
        // Round-up the search center to the closest integer
        mvc_x_array[mvc_index] = (mvc_x_array[mvc_index] + 4) & ~0x07;
        mvc_y_array[mvc_index] = (mvc_y_array[mvc_index] + 4) & ~0x07;

        md_full_pel_search(pcs_ptr,
                           context_ptr,
                           input_picture_ptr,
                           ref_pic,
                           input_origin_index,
                           context_ptr->md_nsq_motion_search_ctrls.use_ssd,
                           mvc_x_array[mvc_index],
                           mvc_y_array[mvc_index],
                           0,
                           0,
                           0,
                           0,
                           1,
                           0,
                           &search_center_mvx,
                           &search_center_mvy,
                           &search_center_cost,
                           hbd_mode_decision);
    }

    *me_mv_x                  = search_center_mvx;
    *me_mv_y                  = search_center_mvy;
    int16_t  best_search_mvx  = (int16_t)~0;
    int16_t  best_search_mvy  = (int16_t)~0;
    uint32_t best_search_cost = (uint32_t)~0;

    md_full_pel_search(pcs_ptr,
                       context_ptr,
                       input_picture_ptr,
                       ref_pic,
                       input_origin_index,
                       context_ptr->md_nsq_motion_search_ctrls.use_ssd,
                       search_center_mvx,
                       search_center_mvy,
                       -(context_ptr->md_nsq_motion_search_ctrls.full_pel_search_width >> 1),
                       +(context_ptr->md_nsq_motion_search_ctrls.full_pel_search_width >> 1),
                       -(context_ptr->md_nsq_motion_search_ctrls.full_pel_search_height >> 1),
                       +(context_ptr->md_nsq_motion_search_ctrls.full_pel_search_height >> 1),
                       1,
                       0,
                       &best_search_mvx,
                       &best_search_mvy,
                       &best_search_cost,
                       hbd_mode_decision);
    if (best_search_cost < search_center_cost) {
        *me_mv_x = best_search_mvx;
        *me_mv_y = best_search_mvy;
    }
}
/*
   clips input MV (in 1/8 precision) to stay within boundaries of a given ref pic
*/
void clip_mv_on_pic_boundary(int32_t blk_origin_x, int32_t blk_origin_y, int32_t bwidth,
                             int32_t bheight, EbPictureBufferDesc *ref_pic, int16_t *mvx,
                             int16_t *mvy) {
    if (blk_origin_x + (*mvx >> 3) + bwidth > ref_pic->max_width + ref_pic->origin_x)
        *mvx = (ref_pic->max_width - blk_origin_x) << 3;

    if (blk_origin_y + (*mvy >> 3) + bheight > ref_pic->max_height + ref_pic->origin_y)
        *mvy = (ref_pic->max_height - blk_origin_y) << 3;

    if (blk_origin_x + (*mvx >> 3) < -ref_pic->origin_x)
        *mvx = (-blk_origin_x - bwidth) << 3;

    if (blk_origin_y + (*mvy >> 3) < -ref_pic->origin_y)
        *mvy = (-blk_origin_y - bheight) << 3;
}
/*
 * Check the size of the spatial MVs and MVPs of the given block
 *
 * Return a motion category, based on the MV size.
 */
uint8_t check_spatial_mv_size(ModeDecisionContext *ctx, uint8_t list_idx, uint8_t ref_idx,
                              int16_t *me_mv_x, int16_t *me_mv_y) {
    uint8_t search_area_multiplier = 0;

    // Iterate over all MVPs; if large, set high search_area_multiplier
    for (int8_t mvp_index = 0; mvp_index < ctx->mvp_count[list_idx][ref_idx]; mvp_index++) {
        if (ctx->mvp_array[list_idx][ref_idx][mvp_index].col > HIGH_SPATIAL_MV_TH ||
            ctx->mvp_array[list_idx][ref_idx][mvp_index].row > HIGH_SPATIAL_MV_TH ||
            *me_mv_x > HIGH_SPATIAL_MV_TH || *me_mv_y > HIGH_SPATIAL_MV_TH) {
            search_area_multiplier = MAX(3, search_area_multiplier);
            return search_area_multiplier; // reached MAX value already
        } else if (ctx->mvp_array[list_idx][ref_idx][mvp_index].col > MEDIUM_SPATIAL_MV_TH ||
                   ctx->mvp_array[list_idx][ref_idx][mvp_index].row > MEDIUM_SPATIAL_MV_TH ||
                   *me_mv_x > MEDIUM_SPATIAL_MV_TH || *me_mv_y > MEDIUM_SPATIAL_MV_TH) {
            search_area_multiplier = MAX(2, search_area_multiplier);
        } else if (ctx->mvp_array[list_idx][ref_idx][mvp_index].col > LOW_SPATIAL_MV_TH ||
                   ctx->mvp_array[list_idx][ref_idx][mvp_index].row > LOW_SPATIAL_MV_TH ||
                   *me_mv_x > LOW_SPATIAL_MV_TH || *me_mv_y > LOW_SPATIAL_MV_TH) {
            search_area_multiplier = MAX(1, search_area_multiplier);
        }
    }
    return search_area_multiplier;
}

/*
 * Check the size of the temporal MVs
 *
 * Return a motion category, based on the MV size.
 */
uint8_t check_temporal_mv_size(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    uint8_t search_area_multiplier = 0;

    Av1Common  *cm             = pcs->parent_pcs_ptr->av1_cm;
    int32_t     mi_row         = ctx->blk_origin_y >> MI_SIZE_LOG2;
    int32_t     mi_col         = ctx->blk_origin_x >> MI_SIZE_LOG2;
    TPL_MV_REF *prev_frame_mvs = pcs->tpl_mvs + (mi_row >> 1) * (cm->mi_stride >> 1) +
        (mi_col >> 1);
    TPL_MV_REF *mv = prev_frame_mvs;
    if (prev_frame_mvs->mfmv0.as_int != INVALID_MV) {
        if (ABS(mv->mfmv0.as_mv.row) > MEDIUM_TEMPORAL_MV_TH ||
            ABS(mv->mfmv0.as_mv.col) > MEDIUM_TEMPORAL_MV_TH) {
            search_area_multiplier = MAX(2, search_area_multiplier);
        } else if (ABS(mv->mfmv0.as_mv.row) > LOW_TEMPORAL_MV_TH ||
                   ABS(mv->mfmv0.as_mv.col) > LOW_TEMPORAL_MV_TH) {
            search_area_multiplier = MAX(1, search_area_multiplier);
        }
    }

    return search_area_multiplier;
}
/*
 * Detect if block has high motion, and if so, perform an expanded ME search.
 */
void md_sq_motion_search(PictureControlSet *pcs, ModeDecisionContext *ctx,
                         EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index,
                         uint8_t list_idx, uint8_t ref_idx, int16_t *me_mv_x, int16_t *me_mv_y) {
    uint8_t              hbd_mode_decision = EB_8_BIT_MD;
    EbReferenceObject   *ref_obj           = pcs->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
    EbPictureBufferDesc *ref_pic = get_ref_pic_buffer(pcs, hbd_mode_decision, list_idx, ref_idx);
    // -------
    // Use scaled references if resolution of the reference is different from that of the input
    // -------
    use_scaled_rec_refs_if_needed(pcs, input_picture_ptr, ref_obj, &ref_pic, hbd_mode_decision);

    MdSqMotionSearchCtrls *md_sq_me_ctrls = &ctx->md_sq_me_ctrls;
    uint16_t               dist           = ABS(
        (int16_t)((int64_t)pcs->picture_number -
                  (int64_t)pcs->parent_pcs_ptr->ref_pic_poc_array[list_idx][ref_idx]));
    uint8_t search_area_multiplier = 0;
    ctx->enable_psad               = md_sq_me_ctrls->enable_psad;
    // Get pa_me distortion and MVs
    int16_t  pa_me_mvx  = (int16_t)~0;
    int16_t  pa_me_mvy  = (int16_t)~0;
    uint32_t pa_me_cost = (uint32_t)~0;
    md_full_pel_search(pcs,
                       ctx,
                       input_picture_ptr,
                       ref_pic,
                       input_origin_index,
                       md_sq_me_ctrls->use_ssd,
                       *me_mv_x,
                       *me_mv_y,
                       0,
                       0,
                       0,
                       0,
                       1,
                       0,
                       &pa_me_mvx,
                       &pa_me_mvy,
                       &pa_me_cost,
                       hbd_mode_decision);

    // Identify potential high active block(s) and ME failure using 2 checks : (1) high ME_MV distortion, (2) active co - located block for non - intra ref(Temporal - MV(s)) or active surrounding block(s) for intra ref(Spatial - MV(s))
    if (ctx->blk_geom->sq_size <= 64) {
        uint32_t fast_lambda = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD]
                                                      : ctx->fast_lambda_md[EB_8_BIT_MD];

        // Check if pa_me distortion is above the per-pixel threshold.  Rate is set to 16.
        if (RDCOST(fast_lambda, 16, pa_me_cost) >
            RDCOST(fast_lambda,
                   16,
                   md_sq_me_ctrls->pame_distortion_th * ctx->blk_geom->bwidth *
                       ctx->blk_geom->bheight)) {
            ref_obj = (EbReferenceObject *)pcs->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;

            search_area_multiplier = !(ref_obj == NULL || ref_obj->frame_type == KEY_FRAME ||
                                       ref_obj->frame_type == INTRA_ONLY_FRAME)
                ? check_temporal_mv_size(pcs, ctx)
                : check_spatial_mv_size(ctx, list_idx, ref_idx, me_mv_x, me_mv_y);
        }
    }

    // If high motion was detected, perform an expanded ME search
    if (search_area_multiplier) {
        int16_t  best_search_mvx  = (int16_t)~0;
        int16_t  best_search_mvy  = (int16_t)~0;
        uint32_t best_search_cost = (uint32_t)~0;

        dist = get_scaled_picture_distance(dist);

        // Sparse-search Level_0
        if (md_sq_me_ctrls->sprs_lev0_enabled) {
            uint16_t sprs_lev0_w = (md_sq_me_ctrls->sprs_lev0_multiplier *
                                    MIN((md_sq_me_ctrls->sprs_lev0_w * search_area_multiplier *
                                         dist),
                                        md_sq_me_ctrls->max_sprs_lev0_w)) /
                100;
            uint16_t sprs_lev0_h = (md_sq_me_ctrls->sprs_lev0_multiplier *
                                    MIN((md_sq_me_ctrls->sprs_lev0_h * search_area_multiplier *
                                         dist),
                                        md_sq_me_ctrls->max_sprs_lev0_h)) /
                100;
            uint8_t sprs_lev0_step = md_sq_me_ctrls->sprs_lev0_step;

            // Derive start/end position of sparse search (must be a multiple of the step size)
            int16_t search_position_start_x = -(((sprs_lev0_w >> 1) / sprs_lev0_step) *
                                                sprs_lev0_step);
            int16_t search_position_end_x   = +(((sprs_lev0_w >> 1) / sprs_lev0_step) *
                                              sprs_lev0_step);
            int16_t search_position_start_y = -(((sprs_lev0_h >> 1) / sprs_lev0_step) *
                                                sprs_lev0_step);
            int16_t search_position_end_y   = +(((sprs_lev0_h >> 1) / sprs_lev0_step) *
                                              sprs_lev0_step);

            ctx->sprs_lev0_start_x = (*me_mv_x >> 3) + search_position_start_x;
            ctx->sprs_lev0_end_x   = (*me_mv_x >> 3) + search_position_end_x;
            ctx->sprs_lev0_start_y = (*me_mv_y >> 3) + search_position_start_y;
            ctx->sprs_lev0_end_y   = (*me_mv_y >> 3) + search_position_end_y;

            md_full_pel_search(pcs,
                               ctx,
                               input_picture_ptr,
                               ref_pic,
                               input_origin_index,
                               md_sq_me_ctrls->use_ssd,
                               *me_mv_x,
                               *me_mv_y,
                               search_position_start_x,
                               search_position_end_x,
                               search_position_start_y,
                               search_position_end_y,
                               sprs_lev0_step,
                               0,
                               &best_search_mvx,
                               &best_search_mvy,
                               &best_search_cost,
                               hbd_mode_decision);

            *me_mv_x = best_search_mvx;
            *me_mv_y = best_search_mvy;
        }

        // Sparse-search Level_1
        if (md_sq_me_ctrls->sprs_lev1_enabled) {
            uint16_t sprs_lev1_w = (md_sq_me_ctrls->sprs_lev1_multiplier *
                                    MIN((md_sq_me_ctrls->sprs_lev1_w * search_area_multiplier *
                                         dist),
                                        md_sq_me_ctrls->max_sprs_lev1_w)) /
                100;
            uint16_t sprs_lev1_h = (md_sq_me_ctrls->sprs_lev1_multiplier *
                                    MIN((md_sq_me_ctrls->sprs_lev1_h * search_area_multiplier *
                                         dist),
                                        md_sq_me_ctrls->max_sprs_lev1_h)) /
                100;
            uint8_t sprs_lev1_step = md_sq_me_ctrls->sprs_lev1_step;

            // Derive start/end position of sparse search (must be a multiple of the step size)
            int16_t search_position_start_x = -(((sprs_lev1_w >> 1) / sprs_lev1_step) *
                                                sprs_lev1_step);
            int16_t search_position_end_x   = +(((sprs_lev1_w >> 1) / sprs_lev1_step) *
                                              sprs_lev1_step);
            int16_t search_position_start_y = -(((sprs_lev1_h >> 1) / sprs_lev1_step) *
                                                sprs_lev1_step);
            int16_t search_position_end_y   = +(((sprs_lev1_h >> 1) / sprs_lev1_step) *
                                              sprs_lev1_step);

            search_position_start_x = (search_position_start_x % 4 == 0)
                ? search_position_start_x - 2
                : search_position_start_x;
            search_position_end_x   = (search_position_end_x % 4 == 0) ? search_position_end_x + 2
                                                                       : search_position_end_x;
            search_position_start_y = (search_position_start_y % 4 == 0)
                ? search_position_start_y - 2
                : search_position_start_y;
            search_position_end_y   = (search_position_end_y % 4 == 0) ? search_position_end_y + 2
                                                                       : search_position_end_y;

            md_full_pel_search(
                pcs,
                ctx,
                input_picture_ptr,
                ref_pic,
                input_origin_index,
                md_sq_me_ctrls->use_ssd,
                *me_mv_x,
                *me_mv_y,
                search_position_start_x,
                search_position_end_x,
                search_position_start_y,
                search_position_end_y,
                sprs_lev1_step,
                (ctx->md_sq_me_ctrls.sprs_lev0_enabled && ctx->md_sq_me_ctrls.sprs_lev0_step == 4)
                    ? 1
                    : 0,
                &best_search_mvx,
                &best_search_mvy,
                &best_search_cost,
                hbd_mode_decision);

            *me_mv_x = best_search_mvx;
            *me_mv_y = best_search_mvy;
        }

        // Sparse-search Level_2
        if (md_sq_me_ctrls->sprs_lev2_enabled) {
            md_full_pel_search(
                pcs,
                ctx,
                input_picture_ptr,
                ref_pic,
                input_origin_index,
                md_sq_me_ctrls->use_ssd,
                *me_mv_x,
                *me_mv_y,
                -(((md_sq_me_ctrls->sprs_lev2_w >> 1) / md_sq_me_ctrls->sprs_lev2_step) *
                  md_sq_me_ctrls->sprs_lev2_step),
                +(((md_sq_me_ctrls->sprs_lev2_w >> 1) / md_sq_me_ctrls->sprs_lev2_step) *
                  md_sq_me_ctrls->sprs_lev2_step),
                -(((md_sq_me_ctrls->sprs_lev2_h >> 1) / md_sq_me_ctrls->sprs_lev2_step) *
                  md_sq_me_ctrls->sprs_lev2_step),
                +(((md_sq_me_ctrls->sprs_lev2_h >> 1) / md_sq_me_ctrls->sprs_lev2_step) *
                  md_sq_me_ctrls->sprs_lev2_step),
                md_sq_me_ctrls->sprs_lev2_step,
                0,
                &best_search_mvx,
                &best_search_mvy,
                &best_search_cost,
                hbd_mode_decision);

            *me_mv_x = best_search_mvx;
            *me_mv_y = best_search_mvy;
        }
    }
}
/*
 * Perform 1/2-Pel, 1/4-Pel, and 1/8-Pel search around the best Full-Pel position
 */
int md_subpel_search(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                     MdSubPelSearchCtrls md_subpel_ctrls, EbPictureBufferDesc *input_picture_ptr,
                     uint8_t list_idx, uint8_t ref_idx, int16_t *me_mv_x, int16_t *me_mv_y) {
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;

    const Av1Common *const cm = pcs_ptr->parent_pcs_ptr->av1_cm;
    MacroBlockD           *xd = context_ptr->blk_ptr->av1xd;

    // ref_mv is used to calculate the cost of the motion vector
    MV ref_mv;
    ref_mv.col = context_ptr->ref_mv.col;
    ref_mv.row = context_ptr->ref_mv.row;
    // High level params
    SUBPEL_MOTION_SEARCH_PARAMS  ms_params_struct;
    SUBPEL_MOTION_SEARCH_PARAMS *ms_params = &ms_params_struct;

    ms_params->allow_hp = md_subpel_ctrls.max_precision == EIGHTH_PEL &&
        pcs_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv;
    ms_params->forced_stop = md_subpel_ctrls.max_precision;
    ms_params->iters_per_step =
        md_subpel_ctrls
            .subpel_iters_per_step; // Maximum number of steps in logarithmic subpel search before giving up.
    // Derive mv_limits (TODO Hsan_Subpel should be derived under md_context @ eack block)
    // Set up limit values for MV components.
    // Mv beyond the range do not produce new/different prediction block.
    MvLimits mv_limits;
    int      mi_row    = xd->mi_row;
    int      mi_col    = xd->mi_col;
    int      mi_width  = mi_size_wide[context_ptr->blk_geom->bsize];
    int      mi_height = mi_size_high[context_ptr->blk_geom->bsize];
    mv_limits.row_min  = -(((mi_row + mi_height) * MI_SIZE) + AOM_INTERP_EXTEND);
    mv_limits.col_min  = -(((mi_col + mi_width) * MI_SIZE) + AOM_INTERP_EXTEND);
    mv_limits.row_max  = (cm->mi_rows - mi_row) * MI_SIZE + AOM_INTERP_EXTEND;
    mv_limits.col_max  = (cm->mi_cols - mi_col) * MI_SIZE + AOM_INTERP_EXTEND;
    svt_av1_set_mv_search_range(&mv_limits, &ref_mv);
    svt_av1_set_subpel_mv_search_range(&ms_params->mv_limits, (FullMvLimits *)&mv_limits, &ref_mv);
    // Mvcost params
    svt_init_mv_cost_params(&ms_params->mv_cost_params,
                            context_ptr,
                            &ref_mv,
                            frm_hdr->quantization_params.base_q_idx,
                            context_ptr->full_lambda_md[EB_8_BIT_MD],
                            0); // 10BIT not supported
    // Subpel variance params
    ms_params->var_params.vfp                = &mefn_ptr[context_ptr->blk_geom->bsize];
    ms_params->var_params.subpel_search_type = md_subpel_ctrls.subpel_search_type;
    ms_params->var_params.w                  = block_size_wide[context_ptr->blk_geom->bsize];
    ms_params->var_params.h                  = block_size_high[context_ptr->blk_geom->bsize];

    // Ref and src buffers
    MSBuffers *ms_buffers = &ms_params->var_params.ms_buffers;

    // Ref buffer
    EbReferenceObject   *ref_obj = pcs_ptr->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
    EbPictureBufferDesc *ref_pic = get_ref_pic_buffer(
        pcs_ptr, 0 /* 10BIT not supported*/, list_idx, ref_idx);
    // -------
    // Use scaled references if resolution of the reference is different from that of the input
    // -------
    use_scaled_rec_refs_if_needed(pcs_ptr, input_picture_ptr, ref_obj, &ref_pic, EB_8_BIT_MD);

    int32_t ref_origin_index = ref_pic->origin_x + context_ptr->blk_origin_x +
        (context_ptr->blk_origin_y + ref_pic->origin_y) * ref_pic->stride_y;

    // Ref buffer
    struct svt_buf_2d ref_struct;
    ref_struct.buf    = ref_pic->buffer_y + ref_origin_index;
    ref_struct.width  = ref_pic->width;
    ref_struct.height = ref_pic->height;
    ref_struct.stride = ref_pic->stride_y;
    ms_buffers->ref   = &ref_struct;
    // Src buffer
    uint32_t input_origin_index = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);
    struct svt_buf_2d src_struct;
    src_struct.buf    = input_picture_ptr->buffer_y + input_origin_index;
    src_struct.width  = input_picture_ptr->width;
    src_struct.height = input_picture_ptr->height;
    src_struct.stride = input_picture_ptr->stride_y;
    ms_buffers->src   = &src_struct;
    int_mv best_mv;
    best_mv.as_mv.col = *me_mv_x >> 3;
    best_mv.as_mv.row = *me_mv_y >> 3;

    int          not_used        = 0;
    MV           subpel_start_mv = get_mv_from_fullmv(&best_mv.as_fullmv);
    unsigned int pred_sse        = 0; // not used
    // Assign which subpel search method to use
    fractional_mv_step_fp *subpel_search_method = md_subpel_ctrls.subpel_search_method ==
            SUBPEL_TREE
        ? svt_av1_find_best_sub_pixel_tree
        : svt_av1_find_best_sub_pixel_tree_pruned;
    ms_params->pred_variance_th                 = md_subpel_ctrls.pred_variance_th;
    ms_params->abs_th_mult                      = md_subpel_ctrls.abs_th_mult;
    ms_params->round_dev_th                     = md_subpel_ctrls.round_dev_th;
    ms_params->skip_diag_refinement             = md_subpel_ctrls.skip_diag_refinement;
    uint8_t early_exit                          = (context_ptr->is_intra_bordered &&
                          context_ptr->cand_reduction_ctrls.use_neighbouring_mode_ctrls.enabled) ||
        (md_subpel_ctrls.skip_zz_mv && best_mv.as_mv.col == 0 && best_mv.as_mv.row == 0);

    int besterr = subpel_search_method(xd,
                                       (const struct AV1Common *const)cm,
                                       ms_params,
                                       subpel_start_mv,
                                       &best_mv.as_mv,
                                       &not_used,
                                       &pred_sse,
                                       pcs_ptr->picture_qp,
                                       context_ptr->blk_geom->bsize,
                                       early_exit);
    *me_mv_x    = best_mv.as_mv.col;
    *me_mv_y    = best_mv.as_mv.row;

    return besterr;
}
// Copy ME_MVs (generated @ PA) from input buffer (pcs_ptr-> .. ->me_results) to local
// MD buffers (context_ptr->sb_me_mv) - simplified for LPD1
void read_refine_me_mvs_light_pd1(PictureControlSet   *pcs_ptr,
                                  EbPictureBufferDesc *input_picture_ptr,
                                  ModeDecisionContext *context_ptr) {
    // derive me offsets
    context_ptr->geom_offset_x = 0;
    context_ptr->geom_offset_y = 0;

    context_ptr->me_cand_offset = context_ptr->me_block_offset *
        pcs_ptr->parent_pcs_ptr->pa_me_data->max_cand;

    // init best ME cost to MAX
    context_ptr->md_me_dist = (uint32_t)~0;

    for (uint32_t ref_it = 0; ref_it < pcs_ptr->parent_pcs_ptr->tot_ref_frame_types; ++ref_it) {
        MvReferenceFrame ref_pair = pcs_ptr->parent_pcs_ptr->ref_frame_type_arr[ref_it];
        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_pair);

        if (rf[1] == NONE_FRAME) {
            uint8_t list_idx = get_list_idx(rf[0]);
            uint8_t ref_idx  = get_ref_frame_idx(rf[0]);

            // Get the ME MV
            const MeSbResults *me_results =
                pcs_ptr->parent_pcs_ptr->pa_me_data->me_results[context_ptr->me_sb_addr];

            if (is_me_data_present(context_ptr, me_results, list_idx, ref_idx)) {
                EbPictureBufferDesc *ref_pic = get_ref_pic_buffer(pcs_ptr, 0, list_idx, ref_idx);
                EbReferenceObject   *ref_obj =
                    (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
                // -------
                // Use scaled references if resolution of the reference is different from that of the input
                // -------
                use_scaled_rec_refs_if_needed(
                    pcs_ptr, input_picture_ptr, ref_obj, &ref_pic, context_ptr->hbd_mode_decision);

                int16_t me_mv_x;
                int16_t me_mv_y;
                if (list_idx == 0) {
                    me_mv_x = (me_results
                                   ->me_mv_array[context_ptr->me_block_offset *
                                                     pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs +
                                                 ref_idx]
                                   .x_mv)
                        << 1;
                    me_mv_y = (me_results
                                   ->me_mv_array[context_ptr->me_block_offset *
                                                     pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs +
                                                 ref_idx]
                                   .y_mv)
                        << 1;
                } else {
                    me_mv_x =
                        (me_results
                             ->me_mv_array[context_ptr->me_block_offset *
                                               pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs +
                                           pcs_ptr->parent_pcs_ptr->pa_me_data->max_l0 + ref_idx]
                             .x_mv)
                        << 1;
                    me_mv_y =
                        (me_results
                             ->me_mv_array[context_ptr->me_block_offset *
                                               pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs +
                                           pcs_ptr->parent_pcs_ptr->pa_me_data->max_l0 + ref_idx]
                             .y_mv)
                        << 1;
                }
                uint8_t skip_subpel =
                    (context_ptr->is_intra_bordered &&
                     context_ptr->cand_reduction_ctrls.use_neighbouring_mode_ctrls.enabled) ||
                    (context_ptr->md_subpel_me_ctrls.skip_zz_mv && me_mv_x == 0 && me_mv_y == 0);
                // can only skip if using dc only b/c otherwise need cost at candidate generation
                skip_subpel &= (!context_ptr->intra_ctrls.enable_intra) ||
                    (context_ptr->intra_ctrls.intra_mode_end == DC_PRED);

                if (context_ptr->md_subpel_me_ctrls.enabled && !skip_subpel) {
                    // Could use context_ptr->mvp_array[list_idx][ref_idx][0], but that requires the single ref MVP array to be init'd, but it is not in light-PD1 path
                    context_ptr->ref_mv.col = context_ptr->shut_fast_rate ||
                            context_ptr->cand_reduction_ctrls.reduce_unipred_candidates >= 3
                        ? 0
                        : (context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                               .ed_ref_mv_stack[rf[0]][0]
                               .this_mv.as_mv.col +
                           4) &
                            ~0x07;
                    context_ptr->ref_mv.row = context_ptr->shut_fast_rate ||
                            context_ptr->cand_reduction_ctrls.reduce_unipred_candidates >= 3
                        ? 0
                        : (context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                               .ed_ref_mv_stack[rf[0]][0]
                               .this_mv.as_mv.row +
                           4) &
                            ~0x07;
                    context_ptr->post_subpel_me_mv_cost[list_idx][ref_idx] = (uint32_t)
                        md_subpel_search(
                            pcs_ptr,
                            context_ptr,
                            context_ptr->md_subpel_me_ctrls,
                            pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr, // 10BIT not supported
                            list_idx,
                            ref_idx,
                            &me_mv_x,
                            &me_mv_y);

                    if (context_ptr->post_subpel_me_mv_cost[list_idx][ref_idx] <
                        context_ptr->md_me_dist)
                        context_ptr->md_me_dist =
                            context_ptr->post_subpel_me_mv_cost[list_idx][ref_idx];
                }
                context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][0] =
                    me_mv_x;
                context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][1] =
                    me_mv_y;
                clip_mv_on_pic_boundary(
                    context_ptr->blk_origin_x,
                    context_ptr->blk_origin_y,
                    context_ptr->blk_geom->bwidth,
                    context_ptr->blk_geom->bheight,
                    ref_pic,
                    &context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][0],
                    &context_ptr
                         ->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][1]);
            }
        }
    }
}
// Copy ME_MVs (generated @ PA) from input buffer (pcs_ptr-> .. ->me_results) to local
// MD buffers (context_ptr->sb_me_mv)
void read_refine_me_mvs(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                        EbPictureBufferDesc *input_picture_ptr) {
    const SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    derive_me_offsets(scs_ptr, pcs_ptr, context_ptr);
    uint8_t hbd_mode_decision = EB_8_BIT_MD;
    input_picture_ptr         = hbd_mode_decision ? pcs_ptr->input_frame16bit
                                                  : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;

    //Update input origin
    uint32_t input_origin_index = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);
    // Get parent_depth_idx_mds
    uint16_t parent_depth_idx_mds = context_ptr->blk_geom->parent_depth_idx_mds;
    context_ptr->md_me_dist       = (uint32_t)~0;
    for (uint32_t ref_it = 0; ref_it < pcs_ptr->parent_pcs_ptr->tot_ref_frame_types; ++ref_it) {
        MvReferenceFrame ref_pair = pcs_ptr->parent_pcs_ptr->ref_frame_type_arr[ref_it];

        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_pair);

        if (rf[1] == NONE_FRAME) {
            uint8_t            list_idx = get_list_idx(rf[0]);
            uint8_t            ref_idx  = get_ref_frame_idx(rf[0]);
            EbReferenceObject *ref_obj  = pcs_ptr->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
            EbPictureBufferDesc *ref_pic = get_ref_pic_buffer(
                pcs_ptr, hbd_mode_decision, list_idx, ref_idx);
            // -------
            // Use scaled references if resolution of the reference is different from that of the input
            // -------
            use_scaled_rec_refs_if_needed(
                pcs_ptr, input_picture_ptr, ref_obj, &ref_pic, hbd_mode_decision);

            // Get the ME MV
            const MeSbResults *me_results =
                pcs_ptr->parent_pcs_ptr->pa_me_data->me_results[context_ptr->me_sb_addr];
            if (is_me_data_present(context_ptr, me_results, list_idx, ref_idx)) {
                int16_t me_mv_x;
                int16_t me_mv_y;
                if (context_ptr->avail_blk_flag[context_ptr->blk_geom->sqi_mds] &&
                    // If NSQ then use the MV of SQ as default MV center
                    (context_ptr->blk_geom->bwidth != context_ptr->blk_geom->bheight) &&
                    // Not applicable for BLOCK_128X64 and BLOCK_64X128 as the 2nd part of each and BLOCK_128X128 do not share the same me_results
                    context_ptr->blk_geom->bsize != BLOCK_64X128 &&
                    context_ptr->blk_geom->bsize != BLOCK_128X64) {
                    me_mv_x = (context_ptr->sb_me_mv[context_ptr->blk_geom->sqi_mds][list_idx]
                                                    [ref_idx][0] +
                               4) &
                        ~0x07;
                    me_mv_y = (context_ptr->sb_me_mv[context_ptr->blk_geom->sqi_mds][list_idx]
                                                    [ref_idx][1] +
                               4) &
                        ~0x07;

                    clip_mv_on_pic_boundary(context_ptr->blk_origin_x,
                                            context_ptr->blk_origin_y,
                                            context_ptr->blk_geom->bwidth,
                                            context_ptr->blk_geom->bheight,
                                            ref_pic,
                                            &me_mv_x,
                                            &me_mv_y);

                } else if (context_ptr->blk_geom->bsize == BLOCK_4X4 &&
                           context_ptr->avail_blk_flag[parent_depth_idx_mds]) {
                    me_mv_x = (context_ptr->sb_me_mv[parent_depth_idx_mds][list_idx][ref_idx][0] +
                               4) &
                        ~0x07;
                    me_mv_y = (context_ptr->sb_me_mv[parent_depth_idx_mds][list_idx][ref_idx][1] +
                               4) &
                        ~0x07;

                    clip_mv_on_pic_boundary(context_ptr->blk_origin_x,
                                            context_ptr->blk_origin_y,
                                            context_ptr->blk_geom->bwidth,
                                            context_ptr->blk_geom->bheight,
                                            ref_pic,
                                            &me_mv_x,
                                            &me_mv_y);

                } else {
                    if (list_idx == 0) {
                        me_mv_x =
                            (me_results
                                 ->me_mv_array[context_ptr->me_block_offset *
                                                   pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs +
                                               ref_idx]
                                 .x_mv)
                            << 1;
                        me_mv_y =
                            (me_results
                                 ->me_mv_array[context_ptr->me_block_offset *
                                                   pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs +
                                               ref_idx]
                                 .y_mv)
                            << 1;
                    } else {
                        me_mv_x =
                            (me_results
                                 ->me_mv_array[context_ptr->me_block_offset *
                                                   pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs +
                                               pcs_ptr->parent_pcs_ptr->pa_me_data->max_l0 +
                                               ref_idx]
                                 .x_mv)
                            << 1;
                        me_mv_y =
                            (me_results
                                 ->me_mv_array[context_ptr->me_block_offset *
                                                   pcs_ptr->parent_pcs_ptr->pa_me_data->max_refs +
                                               pcs_ptr->parent_pcs_ptr->pa_me_data->max_l0 +
                                               ref_idx]
                                 .y_mv)
                            << 1;
                    }
                    clip_mv_on_pic_boundary(context_ptr->blk_origin_x,
                                            context_ptr->blk_origin_y,
                                            context_ptr->blk_geom->bwidth,
                                            context_ptr->blk_geom->bheight,
                                            ref_pic,
                                            &me_mv_x,
                                            &me_mv_y);
                }
                // Set ref MV
                // Could use context_ptr->mvp_array[list_idx][ref_idx][0], but that requires the single ref MVP array to be init'd, but it is not in light-PD1 path
                context_ptr->ref_mv.col = context_ptr->shut_fast_rate
                    ? 0
                    : (context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .ed_ref_mv_stack[rf[0]][0]
                           .this_mv.as_mv.col +
                       4) &
                        ~0x07;
                context_ptr->ref_mv.row = context_ptr->shut_fast_rate
                    ? 0
                    : (context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .ed_ref_mv_stack[rf[0]][0]
                           .this_mv.as_mv.row +
                       4) &
                        ~0x07;
                if ((context_ptr->blk_geom->bwidth != context_ptr->blk_geom->bheight) &&
                    context_ptr->md_nsq_motion_search_ctrls.enabled) {
                    md_nsq_motion_search(pcs_ptr,
                                         context_ptr,
                                         input_picture_ptr,
                                         input_origin_index,
                                         list_idx,
                                         ref_idx,
                                         me_results,
                                         &me_mv_x,
                                         &me_mv_y);
                } else if (context_ptr->md_sq_me_ctrls.enabled) {
                    md_sq_motion_search(pcs_ptr,
                                        context_ptr,
                                        input_picture_ptr,
                                        input_origin_index,
                                        list_idx,
                                        ref_idx,
                                        &me_mv_x,
                                        &me_mv_y);
                }
                context_ptr->post_subpel_me_mv_cost[list_idx][ref_idx] = (int32_t)~0;
                context_ptr->fp_me_mv[list_idx][ref_idx].col           = me_mv_x;
                context_ptr->fp_me_mv[list_idx][ref_idx].row           = me_mv_y;
                if (context_ptr->md_subpel_me_ctrls.enabled) {
                    context_ptr->post_subpel_me_mv_cost[list_idx][ref_idx] = (uint32_t)
                        md_subpel_search(
                            pcs_ptr,
                            context_ptr,
                            context_ptr->md_subpel_me_ctrls,
                            pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr, // 10BIT not supported
                            list_idx,
                            ref_idx,
                            &me_mv_x,
                            &me_mv_y);
                    if (context_ptr->post_subpel_me_mv_cost[list_idx][ref_idx] <
                        context_ptr->md_me_dist)
                        context_ptr->md_me_dist =
                            context_ptr->post_subpel_me_mv_cost[list_idx][ref_idx];
                }
                // Copy ME MV after subpel
                context_ptr->sub_me_mv[list_idx][ref_idx].col = me_mv_x;
                context_ptr->sub_me_mv[list_idx][ref_idx].row = me_mv_y;
                context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][0] =
                    me_mv_x;
                context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][1] =
                    me_mv_y;
                clip_mv_on_pic_boundary(
                    context_ptr->blk_origin_x,
                    context_ptr->blk_origin_y,
                    context_ptr->blk_geom->bwidth,
                    context_ptr->blk_geom->bheight,
                    ref_pic,
                    &context_ptr->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][0],
                    &context_ptr
                         ->sb_me_mv[context_ptr->blk_geom->blkidx_mds][list_idx][ref_idx][1]);
            }
        }
    }
}

/*
Loop over TPL blocks in the SB to update inter information.  Return 1 if the stats for the SB are valid; else return 0.

sb_max_rf_idx: The maximum rf_idx selected by any TPL block in the SB
*/

static Bool get_sb_tpl_inter_stats(PictureControlSet *pcs, ModeDecisionContext *ctx,
                                   uint8_t *sb_inter_selection, uint8_t *sb_max_list0_ref_idx,
                                   uint8_t *sb_max_list1_ref_idx) {
    PictureParentControlSet *ppcs = pcs->parent_pcs_ptr;

    // Check that TPL data is available and that INTRA was tested in TPL.
    // Note that not all INTRA modes may be tested in TPL.
    if (ppcs->tpl_ctrls.enable && ppcs->tpl_src_data_ready &&
        (ppcs->is_used_as_reference_flag || !ppcs->tpl_ctrls.disable_intra_pred_nref)) {
        const int      aligned16_width = (ppcs->aligned_width + 15) >> 4;
        const uint32_t mb_origin_x     = ctx->sb_origin_x;
        const uint32_t mb_origin_y     = ctx->sb_origin_y;
        const int      tpl_blk_size    = ppcs->tpl_ctrls.dispenser_search_level == 0 ? 16
                    : ppcs->tpl_ctrls.dispenser_search_level == 1                    ? 32
                                                                                     : 64;

        // Get actual SB width (for cases of incomplete SBs)
        SbGeom   *sb_geom = &ppcs->sb_geom[ctx->sb_index];
        const int sb_cols = MAX(1, sb_geom->width / tpl_blk_size);
        const int sb_rows = MAX(1, sb_geom->height / tpl_blk_size);

        uint8_t tot_cnt           = 0;
        uint8_t inter_cnt         = 0;
        uint8_t max_list0_ref_idx = 0;
        uint8_t max_list1_ref_idx = 0;

        // Loop over all blocks in the SB
        for (int i = 0; i < sb_rows; i++) {
            TplSrcStats *tpl_src_stats_buffer =
                &ppcs->pa_me_data->tpl_src_stats_buffer[((mb_origin_y >> 4) + i) * aligned16_width +
                                                        (mb_origin_x >> 4)];
            for (int j = 0; j < sb_cols; j++) {
                tot_cnt++;

                if (!is_intra_mode(tpl_src_stats_buffer->best_mode)) {
                    uint8_t list_index    = tpl_src_stats_buffer->best_rf_idx < 4 ? 0 : 1;
                    uint8_t ref_pic_index = tpl_src_stats_buffer->best_rf_idx >= 4
                        ? (tpl_src_stats_buffer->best_rf_idx - 4)
                        : tpl_src_stats_buffer->best_rf_idx;

                    if (list_index)
                        max_list1_ref_idx = MAX(max_list1_ref_idx, ref_pic_index);
                    else
                        max_list0_ref_idx = MAX(max_list0_ref_idx, ref_pic_index);

                    inter_cnt++;
                }

                tpl_src_stats_buffer++;
            }
        }

        *sb_inter_selection   = (inter_cnt * 100) / tot_cnt;
        *sb_max_list0_ref_idx = max_list0_ref_idx;
        *sb_max_list1_ref_idx = max_list1_ref_idx;
        return 1;
    }
    return 0;
}
void perform_md_reference_pruning(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                  EbPictureBufferDesc *input_picture_ptr) {
    uint32_t early_inter_distortion_array[MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH];
    memset(early_inter_distortion_array,
           0xFE,
           sizeof(early_inter_distortion_array[0]) * MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);
    uint32_t offset_tab[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH] = {{0}};
    memset(context_ptr->ref_filtering_res,
           0,
           sizeof(context_ptr->ref_filtering_res[0][0][0]) * TOT_INTER_GROUP *
               MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);
    uint32_t min_dist          = (uint32_t)~0;
    uint8_t  hbd_mode_decision = EB_8_BIT_MD;

    input_picture_ptr = hbd_mode_decision ? pcs_ptr->input_frame16bit
                                          : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;

    // Update input origin
    uint32_t input_origin_index = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);
    int     use_tpl_info = 0;
    uint8_t sb_max_list0_ref_idx;
    uint8_t sb_max_list1_ref_idx;
    uint8_t sb_inter_selection;
    if (context_ptr->ref_pruning_ctrls.use_tpl_info_offset &&
        pcs_ptr->parent_pcs_ptr->tpl_ctrls.enable) {
        if (get_sb_tpl_inter_stats(pcs_ptr,
                                   context_ptr,
                                   &sb_inter_selection,
                                   &sb_max_list0_ref_idx,
                                   &sb_max_list1_ref_idx)) {
            use_tpl_info = 1;
        }
    }
    for (uint32_t ref_it = 0; ref_it < pcs_ptr->parent_pcs_ptr->tot_ref_frame_types; ++ref_it) {
        MvReferenceFrame ref_pair = pcs_ptr->parent_pcs_ptr->ref_frame_type_arr[ref_it];
        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_pair);

        if (rf[1] == NONE_FRAME) {
            uint32_t best_mvp_distortion = (int32_t)~0;
            uint8_t  list_idx            = get_list_idx(rf[0]);
            uint8_t  ref_idx             = get_ref_frame_idx(rf[0]);
            // Only use TPL info if all references are tested
            if (use_tpl_info) {
                if ((list_idx == 0 && ref_idx > sb_max_list0_ref_idx) ||
                    (list_idx == 1 && ref_idx > sb_max_list1_ref_idx))

                    offset_tab[list_idx][ref_idx] +=
                        context_ptr->ref_pruning_ctrls.use_tpl_info_offset;
            }
            // Step 1: derive the best MVP in term of distortion
            for (int8_t mvp_index = 0; mvp_index < context_ptr->mvp_count[list_idx][ref_idx];
                 mvp_index++) {
                // MVP Distortion
                EbReferenceObject *ref_obj =
                    pcs_ptr->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;

                EbPictureBufferDesc *ref_pic = get_ref_pic_buffer(
                    pcs_ptr, hbd_mode_decision, list_idx, ref_idx);
                // -------
                // Use scaled references if resolution of the reference is different from that of the input
                // -------
                use_scaled_rec_refs_if_needed(
                    pcs_ptr, input_picture_ptr, ref_obj, &ref_pic, hbd_mode_decision);

                clip_mv_on_pic_boundary(context_ptr->blk_origin_x,
                                        context_ptr->blk_origin_y,
                                        context_ptr->blk_geom->bwidth,
                                        context_ptr->blk_geom->bheight,
                                        ref_pic,
                                        &context_ptr->mvp_array[list_idx][ref_idx][mvp_index].col,
                                        &context_ptr->mvp_array[list_idx][ref_idx][mvp_index].row);
                // Never be negative here
                int32_t ref_origin_index = ref_pic->origin_x +
                    (context_ptr->blk_origin_x +
                     (context_ptr->mvp_array[list_idx][ref_idx][mvp_index].col >> 3)) +
                    (context_ptr->blk_origin_y +
                     (context_ptr->mvp_array[list_idx][ref_idx][mvp_index].row >> 3) +
                     ref_pic->origin_y) *
                        ref_pic->stride_y;
                assert((context_ptr->blk_geom->bwidth >> 3) < 17);
                uint32_t mvp_distortion;
                if (hbd_mode_decision) {
                    uint16_t *src_10b;
                    DECLARE_ALIGNED(16, uint16_t, packed_buf[PACKED_BUFFER_SIZE]);
                    // pack the reference into temp 16bit buffer
                    uint8_t offset = INTERPOLATION_OFFSET;
                    int32_t stride = STRIDE_PACK;

                    pack_block(ref_pic->buffer_y + ref_origin_index - offset -
                                   (offset * ref_pic->stride_y),
                               ref_pic->stride_y,
                               ref_pic->buffer_bit_inc_y + ref_origin_index - offset -
                                   (offset * ref_pic->stride_bit_inc_y),
                               ref_pic->stride_bit_inc_y,
                               (uint16_t *)packed_buf,
                               stride,
                               context_ptr->blk_geom->bwidth + (offset << 1),
                               context_ptr->blk_geom->bheight + (offset << 1));

                    src_10b = (uint16_t *)packed_buf + offset + (offset * stride);

                    mvp_distortion = sad_16b_kernel(
                        ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                        input_picture_ptr->stride_y,
                        src_10b,
                        stride,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                } else
                    mvp_distortion = svt_nxm_sad_kernel_sub_sampled(
                        input_picture_ptr->buffer_y + input_origin_index,
                        input_picture_ptr->stride_y,
                        ref_pic->buffer_y + ref_origin_index,
                        ref_pic->stride_y,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                if (mvp_distortion < best_mvp_distortion)
                    best_mvp_distortion = mvp_distortion;
            }

            // Evaluate the PA_ME MVs (if available)
            const MeSbResults *me_results =
                pcs_ptr->parent_pcs_ptr->pa_me_data->me_results[context_ptr->me_sb_addr];
            uint32_t pa_me_distortion = (uint32_t)~0; //any non zero value
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
                EbPictureBufferDesc *ref_pic = get_ref_pic_buffer(
                    pcs_ptr, hbd_mode_decision, list_idx, ref_idx);
                // -------
                // Use scaled references if resolution of the reference is different from that of the input
                // -------
                use_scaled_rec_refs_if_needed(
                    pcs_ptr, input_picture_ptr, ref_obj, &ref_pic, hbd_mode_decision);

                clip_mv_on_pic_boundary(context_ptr->blk_origin_x,
                                        context_ptr->blk_origin_y,
                                        context_ptr->blk_geom->bwidth,
                                        context_ptr->blk_geom->bheight,
                                        ref_pic,
                                        &me_mv_x,
                                        &me_mv_y);
                // Never be negative here
                int32_t ref_origin_index = ref_pic->origin_x +
                    (context_ptr->blk_origin_x + (me_mv_x >> 3)) +
                    (context_ptr->blk_origin_y + (me_mv_y >> 3) + ref_pic->origin_y) *
                        ref_pic->stride_y;
                assert((context_ptr->blk_geom->bwidth >> 3) < 17);
                if (hbd_mode_decision) {
                    uint16_t *src_10b;
                    DECLARE_ALIGNED(16, uint16_t, packed_buf[PACKED_BUFFER_SIZE]);
                    // pack the reference into temp 16bit buffer
                    uint8_t offset = INTERPOLATION_OFFSET;
                    int32_t stride = STRIDE_PACK;

                    pack_block(ref_pic->buffer_y + ref_origin_index - offset -
                                   (offset * ref_pic->stride_y),
                               ref_pic->stride_y,
                               ref_pic->buffer_bit_inc_y + ref_origin_index - offset -
                                   (offset * ref_pic->stride_bit_inc_y),
                               ref_pic->stride_bit_inc_y,
                               (uint16_t *)packed_buf,
                               stride,
                               context_ptr->blk_geom->bwidth + (offset << 1),
                               context_ptr->blk_geom->bheight + (offset << 1));

                    src_10b = (uint16_t *)packed_buf + offset + (offset * stride);

                    pa_me_distortion = sad_16b_kernel(
                        ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                        input_picture_ptr->stride_y,
                        src_10b,
                        stride,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
                } else
                    pa_me_distortion = svt_nxm_sad_kernel_sub_sampled(
                        input_picture_ptr->buffer_y + input_origin_index,
                        input_picture_ptr->stride_y,
                        ref_pic->buffer_y + ref_origin_index,
                        ref_pic->stride_y,
                        context_ptr->blk_geom->bheight,
                        context_ptr->blk_geom->bwidth);
            }

            // early_inter_distortion_array
            early_inter_distortion_array[list_idx * REF_LIST_MAX_DEPTH + ref_idx] = MIN(
                pa_me_distortion, best_mvp_distortion);
            if (early_inter_distortion_array[list_idx * REF_LIST_MAX_DEPTH + ref_idx] < min_dist)
                min_dist = early_inter_distortion_array[list_idx * REF_LIST_MAX_DEPTH + ref_idx];
        }
    }
    uint32_t th = (context_ptr->ref_pruning_ctrls.check_closest_multiplier *
                   (context_ptr->blk_geom->bheight * context_ptr->blk_geom->bwidth) *
                   pcs_ptr->picture_qp) /
        24;
    if (context_ptr->ref_pruning_ctrls.check_closest_multiplier &&
        early_inter_distortion_array[0] < th &&
        early_inter_distortion_array[REF_LIST_MAX_DEPTH] < th) {
        for (unsigned li = 0; li < MAX_NUM_OF_REF_PIC_LIST; li++) {
            for (unsigned ri = 0; ri < REF_LIST_MAX_DEPTH; ri++) {
                for (unsigned gi = 0; gi < TOT_INTER_GROUP; gi++) {
                    if (ri == 0 ||
                        context_ptr->ref_pruning_ctrls.max_dev_to_best[gi] == (uint32_t)~0) {
                        context_ptr->ref_filtering_res[gi][li][ri].do_ref = 1;
                    }
                }
            }
        }
    } else {
        // Sort early_inter_distortion_array
        unsigned num_of_cand_to_sort = MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH;
        uint32_t dev_to_the_best[MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH] = {0};
        for (unsigned i = 0; i < num_of_cand_to_sort - 1; ++i)
            dev_to_the_best[i] =
                (((int64_t)MAX(early_inter_distortion_array[i], 1) - MAX(min_dist, 1)) * 100) /
                MAX(min_dist, 1);
        for (unsigned li = 0; li < MAX_NUM_OF_REF_PIC_LIST; li++) {
            for (unsigned ri = 0; ri < REF_LIST_MAX_DEPTH; ri++) {
                for (unsigned gi = 0; gi < TOT_INTER_GROUP; gi++) {
                    uint32_t offset     = offset_tab[li][ri];
                    uint32_t pruning_th = (offset == (uint32_t)~0 ||
                                           context_ptr->ref_pruning_ctrls.max_dev_to_best[gi] == 0)
                        ? 0
                        : (context_ptr->ref_pruning_ctrls.max_dev_to_best[gi] == (uint32_t)~0)
                        ? (uint32_t)~0
                        : MAX(0,
                              ((int64_t)context_ptr->ref_pruning_ctrls.max_dev_to_best[gi] -
                               (int64_t)offset));

                    if (dev_to_the_best[li * REF_LIST_MAX_DEPTH + ri] < pruning_th) {
                        context_ptr->ref_filtering_res[gi][li][ri].do_ref = 1;
                    }
                }
            }
        }
    }
}
/*
 * Read/store all nearest/near MVs for a block for single ref case, and save the best distortion for each ref.
 */
void build_single_ref_mvp_array(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    for (uint32_t ref_it = 0; ref_it < pcs->parent_pcs_ptr->tot_ref_frame_types; ++ref_it) {
        MvReferenceFrame ref_pair = pcs->parent_pcs_ptr->ref_frame_type_arr[ref_it];

        MacroBlockD     *xd = ctx->blk_ptr->av1xd;
        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_pair);
        // Single ref
        if (rf[1] == NONE_FRAME) {
            MvReferenceFrame frame_type = rf[0];
            uint8_t          list_idx   = get_list_idx(rf[0]);
            uint8_t          ref_idx    = get_ref_frame_idx(rf[0]);
            if (ctx->shut_fast_rate) {
                ctx->mvp_array[list_idx][ref_idx][0].col = 0;
                ctx->mvp_array[list_idx][ref_idx][0].row = 0;
                ctx->mvp_count[list_idx][ref_idx]        = 1;
                continue;
            }
            uint8_t drli, max_drl_index;
            int8_t  mvp_count = 0;

            //NEAREST
            ctx->mvp_array[list_idx][ref_idx][mvp_count].col =
                (ctx->md_local_blk_unit[ctx->blk_geom->blkidx_mds]
                     .ed_ref_mv_stack[frame_type][0]
                     .this_mv.as_mv.col +
                 4) &
                ~0x07;
            ctx->mvp_array[list_idx][ref_idx][mvp_count].row =
                (ctx->md_local_blk_unit[ctx->blk_geom->blkidx_mds]
                     .ed_ref_mv_stack[frame_type][0]
                     .this_mv.as_mv.row +
                 4) &
                ~0x07;
            mvp_count++;

            //NEAR
            max_drl_index = get_max_drl_index(xd->ref_mv_count[frame_type], NEARMV);

            for (drli = 0; drli < max_drl_index; drli++) {
                IntMv nearmv = ctx->md_local_blk_unit[ctx->blk_geom->blkidx_mds]
                                   .ed_ref_mv_stack[frame_type][1 + drli]
                                   .this_mv;
                if (((nearmv.as_mv.col + 4) & ~0x07) != ctx->mvp_array[list_idx][ref_idx][0].col &&
                    ((nearmv.as_mv.row + 4) & ~0x07) != ctx->mvp_array[list_idx][ref_idx][0].row) {
                    ctx->mvp_array[list_idx][ref_idx][mvp_count].col = (nearmv.as_mv.col + 4) &
                        ~0x07;
                    ctx->mvp_array[list_idx][ref_idx][mvp_count].row = (nearmv.as_mv.row + 4) &
                        ~0x07;
                    mvp_count++;
                }
            }
            ctx->mvp_count[list_idx][ref_idx] = mvp_count;
        }
    }
}
Bool is_valid_unipred_ref(struct ModeDecisionContext *context_ptr, uint8_t inter_cand_group,
                          uint8_t list_idx, uint8_t ref_idx);
/*
* Performs an ME search around MVP(s)
* For a given (block, list_idx, ref_idx), if PME search is skipped then set ME_MV=ME_MV to preserve PME candidate = (ME_MV, PME_MV)
*/
void pme_search(PictureControlSet *pcs, ModeDecisionContext *ctx,
                EbPictureBufferDesc *input_picture_ptr) {
    uint8_t hbd_mode_decision = EB_8_BIT_MD;

    input_picture_ptr = hbd_mode_decision ? pcs->input_frame16bit
                                          : pcs->parent_pcs_ptr->enhanced_picture_ptr;

    uint32_t input_origin_index = (ctx->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (ctx->blk_origin_x + input_picture_ptr->origin_x);

    ctx->enable_psad = 0;
    for (uint32_t ref_it = 0; ref_it < pcs->parent_pcs_ptr->tot_ref_frame_types; ++ref_it) {
        MvReferenceFrame ref_pair = pcs->parent_pcs_ptr->ref_frame_type_arr[ref_it];
        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_pair);

        // Reset search variable(s)
        uint32_t best_mvp_cost           = (int32_t)~0;
        int16_t  best_search_mvx         = (int16_t)~0;
        int16_t  best_search_mvy         = (int16_t)~0;
        uint32_t pme_mv_cost             = (int32_t)~0;
        uint32_t me_mv_cost              = (int32_t)~0;
        uint32_t post_subpel_pme_mv_cost = (int32_t)~0;

        if (rf[1] == NONE_FRAME) {
            uint8_t list_idx                     = get_list_idx(rf[0]);
            uint8_t ref_idx                      = get_ref_frame_idx(rf[0]);
            ctx->valid_pme_mv[list_idx][ref_idx] = 0;
            EbReferenceObject   *ref_obj = pcs->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
            EbPictureBufferDesc *ref_pic = get_ref_pic_buffer(
                pcs, hbd_mode_decision, list_idx, ref_idx);
            // -------
            // Use scaled references if resolution of the reference is different from that of the input
            // -------
            use_scaled_rec_refs_if_needed(
                pcs, input_picture_ptr, ref_obj, &ref_pic, hbd_mode_decision);
            if (!is_valid_unipred_ref(ctx, PRED_ME_GROUP, list_idx, ref_idx))
                continue;
            // Get the ME MV
            const MeSbResults *me_results =
                pcs->parent_pcs_ptr->pa_me_data->me_results[ctx->me_sb_addr];

            uint8_t me_data_present = is_me_data_present(ctx, me_results, list_idx, ref_idx);

            if (me_data_present) {
                // Early MVP vs. ME_MV check; do not perform PME search for blocks that have a valid ME_MV unless the ME_MV has a different direction than all MVP(s) and the ME_MV mag is higher than MV_TH (not around(0,0))
                if (ctx->md_pme_ctrls.early_check_mv_th_multiplier != MIN_SIGNED_VALUE) {
                    uint8_t is_me_mv_diffrent_than_mvp = 0;
                    for (int8_t mvp_index = 0; mvp_index < ctx->mvp_count[list_idx][ref_idx];
                         mvp_index++) {
                        int16_t mvp_x = ctx->mvp_array[list_idx][ref_idx][mvp_index].col;
                        int16_t mvp_y = ctx->mvp_array[list_idx][ref_idx][mvp_index].row;

                        int mv_th = (((pcs->parent_pcs_ptr->enhanced_picture_ptr->width *
                                       pcs->parent_pcs_ptr->enhanced_picture_ptr->height) >>
                                      17) *
                                     ctx->md_pme_ctrls.early_check_mv_th_multiplier) /
                            10;

                        // Check x direction
                        if (ABS(mvp_x) > mv_th) {
                            if (ctx->fp_me_mv[list_idx][ref_idx].col * mvp_x < 0) {
                                is_me_mv_diffrent_than_mvp = 1;
                                break;
                            }
                        }

                        // Check y direction
                        if (ABS(mvp_y) > mv_th) {
                            if (ctx->fp_me_mv[list_idx][ref_idx].row * mvp_y < 0) {
                                is_me_mv_diffrent_than_mvp = 1;
                                break;
                            }
                        }
                    }

                    if (is_me_mv_diffrent_than_mvp == 0) {
                        ctx->valid_pme_mv[list_idx][ref_idx] = 1;
                        ctx->pme_res[list_idx][ref_idx].dist =
                            ctx->post_subpel_me_mv_cost[list_idx][ref_idx];
                        ctx->best_pme_mv[list_idx][ref_idx][0] =
                            ctx->sub_me_mv[list_idx][ref_idx].col;
                        ctx->best_pme_mv[list_idx][ref_idx][1] =
                            ctx->sub_me_mv[list_idx][ref_idx].row;
                        continue;
                    }
                }
                int16_t me_mv_x;
                int16_t me_mv_y;
                if (list_idx == 0) {
                    me_mv_x = ctx->sb_me_mv[ctx->blk_geom->blkidx_mds][REF_LIST_0][ref_idx][0];
                    me_mv_y = ctx->sb_me_mv[ctx->blk_geom->blkidx_mds][REF_LIST_0][ref_idx][1];
                } else {
                    me_mv_x = ctx->sb_me_mv[ctx->blk_geom->blkidx_mds][REF_LIST_1][ref_idx][0];
                    me_mv_y = ctx->sb_me_mv[ctx->blk_geom->blkidx_mds][REF_LIST_1][ref_idx][1];
                }
                // Round-up to the closest integer the ME MV
                me_mv_x = (me_mv_x + 4) & ~0x07;
                me_mv_y = (me_mv_y + 4) & ~0x07;

                // Set a ref MV (nearest) for the ME MV
                ctx->ref_mv.col = ctx->mvp_array[list_idx][ref_idx][0].col;
                ctx->ref_mv.row = ctx->mvp_array[list_idx][ref_idx][0].row;
                md_full_pel_search(pcs,
                                   ctx,
                                   input_picture_ptr,
                                   ref_pic,
                                   input_origin_index,
                                   ctx->md_pme_ctrls.use_ssd,
                                   me_mv_x,
                                   me_mv_y,
                                   0,
                                   0,
                                   0,
                                   0,
                                   1,
                                   0,
                                   &me_mv_x,
                                   &me_mv_y,
                                   &me_mv_cost,
                                   hbd_mode_decision);
            }

            // Step 1: derive the best MVP in term of distortion
            int16_t best_mvp_x = 0;
            int16_t best_mvp_y = 0;

            for (int8_t mvp_index = 0; mvp_index < ctx->mvp_count[list_idx][ref_idx]; mvp_index++) {
                // Set a ref MV (MVP under eval) for the MVP under eval
                ctx->ref_mv.col = ctx->mvp_array[list_idx][ref_idx][mvp_index].col;
                ctx->ref_mv.row = ctx->mvp_array[list_idx][ref_idx][mvp_index].row;

                clip_mv_on_pic_boundary(ctx->blk_origin_x,
                                        ctx->blk_origin_y,
                                        ctx->blk_geom->bwidth,
                                        ctx->blk_geom->bheight,
                                        ref_pic,
                                        &ctx->mvp_array[list_idx][ref_idx][mvp_index].col,
                                        &ctx->mvp_array[list_idx][ref_idx][mvp_index].row);

                md_full_pel_search(pcs,
                                   ctx,
                                   input_picture_ptr,
                                   ref_pic,
                                   input_origin_index,
                                   ctx->md_pme_ctrls.use_ssd,
                                   ctx->mvp_array[list_idx][ref_idx][mvp_index].col,
                                   ctx->mvp_array[list_idx][ref_idx][mvp_index].row,
                                   0,
                                   0,
                                   0,
                                   0,
                                   1,
                                   0,
                                   &best_mvp_x,
                                   &best_mvp_y,
                                   &best_mvp_cost,
                                   hbd_mode_decision);
            }
            if (me_data_present) {
                int64_t pme_to_me_cost_dev =
                    (((int64_t)MAX(best_mvp_cost, 1) - (int64_t)MAX(me_mv_cost, 1)) * 100) /
                    (int64_t)MAX(me_mv_cost, 1);

                if ((ABS(ctx->fp_me_mv[list_idx][ref_idx].col - best_mvp_x) <=
                         ctx->md_pme_ctrls.pre_fp_pme_to_me_mv_th &&
                     ABS(ctx->fp_me_mv[list_idx][ref_idx].row - best_mvp_y) <=
                         ctx->md_pme_ctrls.pre_fp_pme_to_me_mv_th) ||
                    pme_to_me_cost_dev >= ctx->md_pme_ctrls.pre_fp_pme_to_me_cost_th) {
                    ctx->valid_pme_mv[list_idx][ref_idx] = 1;
                    ctx->pme_res[list_idx][ref_idx].dist =
                        ctx->post_subpel_me_mv_cost[list_idx][ref_idx];
                    ctx->best_pme_mv[list_idx][ref_idx][0] = ctx->sub_me_mv[list_idx][ref_idx].col;
                    ctx->best_pme_mv[list_idx][ref_idx][1] = ctx->sub_me_mv[list_idx][ref_idx].row;
                    continue;
                }
            }
            // Set ref MV
            ctx->ref_mv.col  = best_mvp_x;
            ctx->ref_mv.row  = best_mvp_y;
            ctx->enable_psad = ctx->md_pme_ctrls.enable_psad;
            md_full_pel_search(pcs,
                               ctx,
                               input_picture_ptr,
                               ref_pic,
                               input_origin_index,
                               ctx->md_pme_ctrls.use_ssd,
                               best_mvp_x,
                               best_mvp_y,
                               -(ctx->md_pme_ctrls.full_pel_search_width >> 1),
                               +(ctx->md_pme_ctrls.full_pel_search_width >> 1),
                               -(ctx->md_pme_ctrls.full_pel_search_height >> 1),
                               +(ctx->md_pme_ctrls.full_pel_search_height >> 1),
                               1,
                               0,
                               &best_search_mvx,
                               &best_search_mvy,
                               &pme_mv_cost,
                               hbd_mode_decision);
            if (me_data_present) {
                int64_t pme_to_me_cost_dev =
                    (((int64_t)MAX(pme_mv_cost, 1) - (int64_t)MAX(me_mv_cost, 1)) * 100) /
                    (int64_t)MAX(me_mv_cost, 1);

                if ((ABS(ctx->fp_me_mv[list_idx][ref_idx].col - best_search_mvx) <=
                         ctx->md_pme_ctrls.post_fp_pme_to_me_mv_th &&
                     ABS(ctx->fp_me_mv[list_idx][ref_idx].row - best_search_mvy) <=
                         ctx->md_pme_ctrls.post_fp_pme_to_me_mv_th) ||
                    pme_to_me_cost_dev >= ctx->md_pme_ctrls.post_fp_pme_to_me_cost_th) {
                    ctx->valid_pme_mv[list_idx][ref_idx] = 1;
                    ctx->pme_res[list_idx][ref_idx].dist =
                        ctx->post_subpel_me_mv_cost[list_idx][ref_idx];
                    ctx->best_pme_mv[list_idx][ref_idx][0] = ctx->sub_me_mv[list_idx][ref_idx].col;
                    ctx->best_pme_mv[list_idx][ref_idx][1] = ctx->sub_me_mv[list_idx][ref_idx].row;
                    continue;
                }
            }
            if (ctx->md_subpel_pme_ctrls.enabled) {
                post_subpel_pme_mv_cost = (uint32_t)md_subpel_search(
                    pcs,
                    ctx,
                    ctx->md_subpel_pme_ctrls,
                    pcs->parent_pcs_ptr->enhanced_picture_ptr, // 10BIT not supported
                    list_idx,
                    ref_idx,
                    &best_search_mvx,
                    &best_search_mvy);
            }

            ctx->best_pme_mv[list_idx][ref_idx][0] = best_search_mvx;
            ctx->best_pme_mv[list_idx][ref_idx][1] = best_search_mvy;
            ctx->valid_pme_mv[list_idx][ref_idx]   = 1;
            ctx->pme_res[list_idx][ref_idx].dist   = post_subpel_pme_mv_cost;
        }
    }
}
void av1_cost_calc_cfl(PictureControlSet *pcs_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
                       ModeDecisionContext *context_ptr, uint32_t component_mask,
                       EbPictureBufferDesc *input_picture_ptr, uint32_t input_cb_origin_in_index,
                       uint32_t blk_chroma_origin_index, uint64_t full_distortion[DIST_CALC_TOTAL],
                       uint64_t *coeff_bits, Bool check_dc) {
    ModeDecisionCandidate *candidate_ptr = candidate_buffer->candidate_ptr;
    uint32_t               count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];
    uint64_t               cb_full_distortion[DIST_CALC_TOTAL];
    uint64_t               cr_full_distortion[DIST_CALC_TOTAL];
    uint32_t               chroma_width  = context_ptr->blk_geom->bwidth_uv;
    uint32_t               chroma_height = context_ptr->blk_geom->bheight_uv;
    // FullLoop and TU search
    uint16_t cb_qindex = context_ptr->qp_index;

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
        uint64_t cb_coeff_bits                   = 0;
        uint64_t cr_coeff_bits                   = 0;
        int32_t  alpha_q3                        = (check_dc) ? 0
                                                              : cfl_idx_to_alpha(candidate_ptr->cfl_alpha_idx,
                                                         candidate_ptr->cfl_alpha_signs,
                                                         CFL_PRED_U); // once for U, once for V
        assert(chroma_width * CFL_BUF_LINE + chroma_height <= CFL_BUF_SQUARE);

        if (!context_ptr->hbd_mode_decision) {
            svt_cfl_predict_lbd(
                context_ptr->pred_buf_q3,
                &(candidate_buffer->prediction_ptr->buffer_cb[blk_chroma_origin_index]),
                candidate_buffer->prediction_ptr->stride_cb,
                &(context_ptr->scratch_prediction_ptr->buffer_cb[blk_chroma_origin_index]),
                context_ptr->scratch_prediction_ptr->stride_cb,
                alpha_q3,
                8,
                chroma_width,
                chroma_height);
        } else {
            svt_cfl_predict_hbd(
                context_ptr->pred_buf_q3,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cb) + blk_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cb,
                ((uint16_t *)context_ptr->scratch_prediction_ptr->buffer_cb) +
                    blk_chroma_origin_index,
                context_ptr->scratch_prediction_ptr->stride_cb,
                alpha_q3,
                10,
                chroma_width,
                chroma_height);
        }

        // Cb Residual
        residual_kernel(input_picture_ptr->buffer_cb,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cb,
                        context_ptr->scratch_prediction_ptr->buffer_cb,
                        blk_chroma_origin_index,
                        context_ptr->scratch_prediction_ptr->stride_cb,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cb,
                        blk_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cb,
                        context_ptr->hbd_mode_decision,
                        chroma_width,
                        chroma_height);
        svt_aom_full_loop_uv(pcs_ptr,
                             context_ptr,
                             candidate_buffer,
                             input_picture_ptr,
                             COMPONENT_CHROMA_CB,
                             cb_qindex,
                             count_non_zero_coeffs,
                             cb_full_distortion,
                             cr_full_distortion,
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

        uint64_t cb_coeff_bits = 0;
        uint64_t cr_coeff_bits = 0;
        int32_t  alpha_q3      = check_dc ? 0
                                          : cfl_idx_to_alpha(candidate_ptr->cfl_alpha_idx,
                                                       candidate_ptr->cfl_alpha_signs,
                                                       CFL_PRED_V); // once for U, once for V
        assert(chroma_width * CFL_BUF_LINE + chroma_height <= CFL_BUF_SQUARE);

        if (!context_ptr->hbd_mode_decision) {
            svt_cfl_predict_lbd(
                context_ptr->pred_buf_q3,
                &(candidate_buffer->prediction_ptr->buffer_cr[blk_chroma_origin_index]),
                candidate_buffer->prediction_ptr->stride_cr,
                &(context_ptr->scratch_prediction_ptr->buffer_cr[blk_chroma_origin_index]),
                context_ptr->scratch_prediction_ptr->stride_cr,
                alpha_q3,
                8,
                chroma_width,
                chroma_height);
        } else {
            svt_cfl_predict_hbd(
                context_ptr->pred_buf_q3,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cr) + blk_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cr,
                ((uint16_t *)context_ptr->scratch_prediction_ptr->buffer_cr) +
                    blk_chroma_origin_index,
                context_ptr->scratch_prediction_ptr->stride_cr,
                alpha_q3,
                10,
                chroma_width,
                chroma_height);
        }

        // Cr Residual
        residual_kernel(input_picture_ptr->buffer_cr,
                        input_cb_origin_in_index,
                        input_picture_ptr->stride_cr,
                        context_ptr->scratch_prediction_ptr->buffer_cr,
                        blk_chroma_origin_index,
                        context_ptr->scratch_prediction_ptr->stride_cr,
                        (int16_t *)candidate_buffer->residual_ptr->buffer_cr,
                        blk_chroma_origin_index,
                        candidate_buffer->residual_ptr->stride_cr,
                        context_ptr->hbd_mode_decision,
                        chroma_width,
                        chroma_height);
        svt_aom_full_loop_uv(pcs_ptr,
                             context_ptr,
                             candidate_buffer,
                             input_picture_ptr,
                             COMPONENT_CHROMA_CR,
                             cb_qindex,
                             count_non_zero_coeffs,
                             cb_full_distortion,
                             cr_full_distortion,
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
/************************************************************************************************
* md_cfl_rd_pick_alpha
* Pick the best alpha for cfl mode or Choose DC
************************************************************************************************/
void md_cfl_rd_pick_alpha(PictureControlSet *pcs_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
                          ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                          uint32_t input_cb_origin_in_index, uint32_t blk_chroma_origin_index) {
    int64_t  best_rd = INT64_MAX;
    uint64_t full_distortion[DIST_CALC_TOTAL];
    uint64_t coeff_bits;

    uint32_t      full_lambda = context_ptr->hbd_mode_decision
             ? context_ptr->full_lambda_md[EB_10_BIT_MD]
             : context_ptr->full_lambda_md[EB_8_BIT_MD];
    const int64_t mode_rd     = RDCOST(
        full_lambda,
        (uint64_t)context_ptr->md_rate_estimation_ptr
            ->intra_uv_mode_fac_bits[CFL_ALLOWED][candidate_buffer->candidate_ptr->pred_mode]
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
                context_ptr->md_rate_estimation_ptr->cfl_alpha_fac_bits[joint_sign][plane][0];
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
                if (c > context_ptr->cfl_ctrls.itr_th && progress < c)
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

                    const int32_t alpha_rate = context_ptr->md_rate_estimation_ptr
                                                   ->cfl_alpha_fac_bits[joint_sign][plane][c];
                    int64_t this_rd = RDCOST(
                        full_lambda, coeff_bits + alpha_rate, full_distortion[DIST_CALC_RESIDUAL]);
                    if (this_rd >= best_rd_uv[joint_sign][plane])
                        continue;
                    best_rd_uv[joint_sign][plane] = this_rd;
                    best_c[joint_sign][plane]     = c;
                    flag                          = context_ptr->cfl_ctrls.itr_th;
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
    const int64_t dc_mode_rd                         = RDCOST(
        full_lambda,
        context_ptr->md_rate_estimation_ptr
            ->intra_uv_mode_fac_bits[CFL_ALLOWED][candidate_buffer->candidate_ptr->pred_mode]
                                    [UV_DC_PRED],
        0);
    av1_cost_calc_cfl(pcs_ptr,
                      candidate_buffer,
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
    if (dc_rd <= best_rd || best_rd == INT64_MAX) {
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
/************************************************************************************************
* cfl_prediction
* Performed cfl prediction in the following steps:
* 1: recon the Luma
* 2: Form the pred_buf_q3
* 3: Loop over alphas and find the best or choose DC
* 4: Recalculate the residual for chroma
************************************************************************************************/
static void cfl_prediction(PictureControlSet           *pcs_ptr,
                           ModeDecisionCandidateBuffer *candidate_buffer,
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
            svt_cfl_luma_subsampling_420_lbd(
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
            svt_cfl_luma_subsampling_420_hbd(
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

        svt_subtract_average(context_ptr->pred_buf_q3,
                             chroma_width,
                             chroma_height,
                             round_offset,
                             svt_log2f(chroma_width) + svt_log2f(chroma_height));

        // 3: Loop over alphas and find the best or choose DC
        md_cfl_rd_pick_alpha(pcs_ptr,
                             candidate_buffer,
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
                svt_cfl_predict_lbd(
                    context_ptr->pred_buf_q3,
                    &(candidate_buffer->prediction_ptr->buffer_cb[blk_chroma_origin_index]),
                    candidate_buffer->prediction_ptr->stride_cb,
                    &(candidate_buffer->prediction_ptr->buffer_cb[blk_chroma_origin_index]),
                    candidate_buffer->prediction_ptr->stride_cb,
                    alpha_q3_cb,
                    8,
                    chroma_width,
                    chroma_height);

                svt_cfl_predict_lbd(
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
                svt_cfl_predict_hbd(context_ptr->pred_buf_q3,
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

                svt_cfl_predict_hbd(context_ptr->pred_buf_q3,
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
static INLINE TxType av1_get_tx_type(int32_t is_inter, PredictionMode pred_mode,
                                     UvPredictionMode pred_mode_uv, PlaneType plane_type,
                                     TxSize tx_size, int32_t reduced_tx_set) {
    if (txsize_sqr_up_map[tx_size] > TX_32X32 || plane_type == PLANE_TYPE_Y || is_inter) {
        return DCT_DCT;
    }
    // In intra mode, uv planes don't share the same prediction mode as y
    // plane, so the tx_type should not be shared
    TxType tx_type = intra_mode_to_tx_type(pred_mode, pred_mode_uv, PLANE_TYPE_UV);
    assert(tx_type < TX_TYPES);
    const TxSetType tx_set_type = get_ext_tx_set_type(tx_size, is_inter, reduced_tx_set);
    return !av1_ext_tx_used[tx_set_type][tx_type] ? DCT_DCT : tx_type;
}
void check_best_indepedant_cfl(PictureControlSet *pcs_ptr, EbPictureBufferDesc *input_picture_ptr,
                               ModeDecisionContext *context_ptr, uint32_t input_cb_origin_in_index,
                               uint32_t                     blk_chroma_origin_index,
                               ModeDecisionCandidateBuffer *candidate_buffer, uint8_t cb_qindex,
                               uint64_t *cb_full_distortion, uint64_t *cr_full_distortion,
                               uint64_t *cb_coeff_bits, uint64_t *cr_coeff_bits) {
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_lambda_md[EB_8_BIT_MD];
    if (candidate_buffer->candidate_ptr->filter_intra_mode != FILTER_INTRA_MODES)
        assert(candidate_buffer->candidate_ptr->pred_mode == DC_PRED);
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    // cfl cost
    uint64_t chroma_rate = 0;
    if (candidate_buffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {
        chroma_rate +=
            context_ptr->md_rate_estimation_ptr
                ->cfl_alpha_fac_bits[candidate_buffer->candidate_ptr->cfl_alpha_signs][CFL_PRED_U]
                                    [CFL_IDX_U(candidate_buffer->candidate_ptr->cfl_alpha_idx)] +
            context_ptr->md_rate_estimation_ptr
                ->cfl_alpha_fac_bits[candidate_buffer->candidate_ptr->cfl_alpha_signs][CFL_PRED_V]
                                    [CFL_IDX_V(candidate_buffer->candidate_ptr->cfl_alpha_idx)];
        chroma_rate += (uint64_t)context_ptr->md_rate_estimation_ptr->intra_uv_mode_fac_bits
                           [CFL_ALLOWED][candidate_buffer->candidate_ptr->pred_mode][UV_CFL_PRED];
    } else
        chroma_rate = (uint64_t)context_ptr->md_rate_estimation_ptr->intra_uv_mode_fac_bits
                          [CFL_ALLOWED][candidate_buffer->candidate_ptr->pred_mode][UV_DC_PRED];
    int      coeff_rate  = (int)(*cb_coeff_bits + *cr_coeff_bits);
    int      distortion  = (int)(cb_full_distortion[DIST_CALC_RESIDUAL] +
                           cr_full_distortion[DIST_CALC_RESIDUAL]);
    int      rate        = (int)(coeff_rate + chroma_rate + candidate_buffer->fast_luma_rate);
    uint64_t cfl_uv_cost = RDCOST(full_lambda, rate, distortion);
    // cfl vs. best independant
    if (context_ptr->best_uv_cost[candidate_buffer->candidate_ptr->pred_mode]
                                 [3 + candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]] <
        cfl_uv_cost) {
        // Update the current candidate
        candidate_buffer->candidate_ptr->intra_chroma_mode =
            context_ptr->best_uv_mode[candidate_buffer->candidate_ptr->pred_mode]
                                     [MAX_ANGLE_DELTA +
                                      candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]];
        candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_UV] =
            context_ptr->best_uv_angle[candidate_buffer->candidate_ptr->pred_mode]
                                      [MAX_ANGLE_DELTA +
                                       candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]];
        // check if candidate_buffer->candidate_ptr->fast_luma_rate = context_ptr->fast_luma_rate[candidate_buffer->candidate_ptr->intra_luma_mode];
        candidate_buffer->fast_chroma_rate =
            context_ptr
                ->fast_chroma_rate[candidate_buffer->candidate_ptr->pred_mode]
                                  [MAX_ANGLE_DELTA +
                                   candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]];
        candidate_buffer->candidate_ptr->transform_type_uv = av1_get_tx_type(
            0, // is_inter
            (PredictionMode)0,
            (UvPredictionMode)context_ptr
                ->best_uv_mode[candidate_buffer->candidate_ptr->pred_mode]
                              [3 + candidate_buffer->candidate_ptr->angle_delta[PLANE_TYPE_Y]],
            PLANE_TYPE_UV,
            context_ptr->blk_geom->txsize_uv[0][0],
            frm_hdr->reduced_tx_set);
        context_ptr->uv_intra_comp_only = TRUE;
        memset(candidate_buffer->eob[1], 0, sizeof(uint16_t));
        memset(candidate_buffer->eob[2], 0, sizeof(uint16_t));
        candidate_buffer->u_has_coeff            = 0;
        candidate_buffer->v_has_coeff            = 0;
        cb_full_distortion[DIST_CALC_RESIDUAL]   = 0;
        cr_full_distortion[DIST_CALC_RESIDUAL]   = 0;
        cb_full_distortion[DIST_CALC_PREDICTION] = 0;
        cr_full_distortion[DIST_CALC_PREDICTION] = 0;

        *cb_coeff_bits = 0;
        *cr_coeff_bits = 0;

        uint32_t count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];
        context_ptr->md_staging_skip_chroma_pred = FALSE;
        svt_product_prediction_fun_table[is_inter_mode(candidate_buffer->candidate_ptr->pred_mode)](
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
        svt_aom_full_loop_uv(pcs_ptr,
                             context_ptr,
                             candidate_buffer,
                             input_picture_ptr,
                             COMPONENT_CHROMA,
                             cb_qindex,
                             count_non_zero_coeffs,
                             cb_full_distortion,
                             cr_full_distortion,
                             cb_coeff_bits,
                             cr_coeff_bits,
                             1);
    }
}

// double check the usage of tx_search_luma_recon_neighbor_array16bit
EbErrorType av1_intra_luma_prediction(ModeDecisionContext         *md_context_ptr,
                                      PictureControlSet           *pcs_ptr,
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

    MacroBlockD *xd                      = md_context_ptr->blk_ptr->av1xd;
    md_context_ptr->intra_luma_left_mode = DC_PRED;
    md_context_ptr->intra_luma_top_mode  = DC_PRED;
    if (xd->left_available)
        md_context_ptr->intra_luma_left_mode = xd->mi[-1]->mbmi.block_mi.mode >= NEARESTMV
            ? DC_PRED
            : xd->mi[-1]->mbmi.block_mi.mode;
    if (xd->up_available)
        md_context_ptr->intra_luma_top_mode = xd->mi[-xd->mi_stride]->mbmi.block_mi.mode >=
                NEARESTMV
            ? DC_PRED
            : xd->mi[-xd->mi_stride]->mbmi.block_mi.mode;
    TxSize tx_size =
        md_context_ptr->blk_geom->txsize[md_context_ptr->tx_depth][md_context_ptr->txb_itr];
    uint32_t sb_size_luma = pcs_ptr->parent_pcs_ptr->scs_ptr->sb_size;

    PredictionMode mode;
    if (!md_context_ptr->hbd_mode_decision) {
        uint8_t top_neigh_array[64 * 2 + 1];
        uint8_t left_neigh_array[64 * 2 + 1];

        mode = candidate_buffer_ptr->candidate_ptr->pred_mode;
        if (candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y] == 0) {
            IntraSize intra_size = intra_unit[mode];
            if (txb_origin_y != 0 && intra_size.top)
                svt_memcpy(
                    top_neigh_array + 1,
                    md_context_ptr->tx_search_luma_recon_neighbor_array->top_array + txb_origin_x,
                    tx_width * intra_size.top);
            if (txb_origin_x != 0 && intra_size.left) {
                uint16_t multipler = (txb_origin_y % sb_size_luma + tx_height * intra_size.left) >
                        sb_size_luma
                    ? 1
                    : intra_size.left;
                svt_memcpy(
                    left_neigh_array + 1,
                    md_context_ptr->tx_search_luma_recon_neighbor_array->left_array + txb_origin_y,
                    tx_height * multipler);
            }

        } else {
            if (txb_origin_y != 0)
                svt_memcpy(
                    top_neigh_array + 1,
                    md_context_ptr->tx_search_luma_recon_neighbor_array->top_array + txb_origin_x,
                    tx_width * 2);
            if (txb_origin_x != 0) {
                uint16_t multipler = (txb_origin_y % sb_size_luma + tx_height * 2) > sb_size_luma
                    ? 1
                    : 2;
                svt_memcpy(
                    left_neigh_array + 1,
                    md_context_ptr->tx_search_luma_recon_neighbor_array->left_array + txb_origin_y,
                    tx_height * multipler);
            }
        }

        if (txb_origin_y != 0 && txb_origin_x != 0)
            top_neigh_array[0] = left_neigh_array[0] =
                md_context_ptr->tx_search_luma_recon_neighbor_array->top_left_array
                    [md_context_ptr->tx_search_luma_recon_neighbor_array->max_pic_h + txb_origin_x -
                     txb_origin_y];
        svt_av1_predict_intra_block(
            !ED_STAGE,
            md_context_ptr->blk_geom,
            md_context_ptr->blk_ptr->av1xd,
            md_context_ptr->blk_geom->bwidth,
            md_context_ptr->blk_geom->bheight,
            tx_size,
            mode, //PredictionMode mode,
            candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
            candidate_buffer_ptr->candidate_ptr->palette_info
                ? (candidate_buffer_ptr->candidate_ptr->palette_size[0] > 0)
                : 0,
            candidate_buffer_ptr->candidate_ptr->palette_info, //ATB MD
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
            &pcs_ptr->scs_ptr->seq_header);
    } else {
        uint16_t top_neigh_array[64 * 2 + 1];
        uint16_t left_neigh_array[64 * 2 + 1];

        mode = candidate_buffer_ptr->candidate_ptr->pred_mode;

        if (candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y] == 0) {
            IntraSize intra_size = intra_unit[mode];

            if (txb_origin_y != 0 && intra_size.top)
                svt_memcpy(top_neigh_array + 1,
                           (uint16_t *)(md_context_ptr->tx_search_luma_recon_neighbor_array16bit
                                            ->top_array) +
                               txb_origin_x,
                           sizeof(uint16_t) * tx_width * intra_size.top);
            if (txb_origin_x != 0 && intra_size.left) {
                uint16_t multipler = (txb_origin_y % sb_size_luma + tx_height * intra_size.left) >
                        sb_size_luma
                    ? 1
                    : intra_size.left;
                svt_memcpy(left_neigh_array + 1,
                           (uint16_t *)(md_context_ptr->tx_search_luma_recon_neighbor_array16bit
                                            ->left_array) +
                               txb_origin_y,
                           sizeof(uint16_t) * tx_height * multipler);
            }

        } else {
            if (txb_origin_y != 0)
                svt_memcpy(top_neigh_array + 1,
                           (uint16_t *)(md_context_ptr->tx_search_luma_recon_neighbor_array16bit
                                            ->top_array) +
                               txb_origin_x,
                           sizeof(uint16_t) * tx_width * 2);
            if (txb_origin_x != 0) {
                uint16_t multipler = (txb_origin_y % sb_size_luma + tx_height * 2) > sb_size_luma
                    ? 1
                    : 2;
                svt_memcpy(left_neigh_array + 1,
                           (uint16_t *)(md_context_ptr->tx_search_luma_recon_neighbor_array16bit
                                            ->left_array) +
                               txb_origin_y,
                           sizeof(uint16_t) * tx_height * multipler);
            }
        }

        if (txb_origin_y != 0 && txb_origin_x != 0)
            top_neigh_array[0] = left_neigh_array[0] =
                ((uint16_t *)(md_context_ptr->tx_search_luma_recon_neighbor_array16bit
                                  ->top_left_array) +
                 md_context_ptr->tx_search_luma_recon_neighbor_array16bit->max_pic_h +
                 txb_origin_x - txb_origin_y)[0];

        svt_av1_predict_intra_block_16bit(
            EB_TEN_BIT,
            !ED_STAGE,
            md_context_ptr->blk_geom,
            md_context_ptr->blk_ptr->av1xd,
            md_context_ptr->blk_geom->bwidth,
            md_context_ptr->blk_geom->bheight,
            tx_size,
            mode,
            candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
            candidate_buffer_ptr->candidate_ptr->palette_info
                ? (candidate_buffer_ptr->candidate_ptr->palette_size[0] > 0)
                : 0,
            candidate_buffer_ptr->candidate_ptr->palette_info, //ATB MD
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
            &pcs_ptr->scs_ptr->seq_header);
    }

    return return_error;
}

static void tx_search_update_recon_sample_neighbor_array(
    NeighborArrayUnit *lumaReconSampleNeighborArray, EbPictureBufferDesc *recon_buffer,
    uint32_t txb_origin_x, uint32_t txb_origin_y, uint32_t input_origin_x, uint32_t input_origin_y,
    uint32_t width, uint32_t height, Bool hbd) {
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
        bsize == BLOCK_32X16 || bsize == BLOCK_16X8 || bsize == BLOCK_8X16 ||
        bsize == BLOCK_64X16 || bsize == BLOCK_16X64 || bsize == BLOCK_32X8 ||
        bsize == BLOCK_8X32 || bsize == BLOCK_16X4 || bsize == BLOCK_4X16)
        tx_depth = 2;
    else if (bsize == BLOCK_8X8)
        tx_depth = 1;
    // tx_depth=0 if BLOCK_8X4, BLOCK_4X8, BLOCK_4X4, BLOCK_128X128, BLOCK_128X64, BLOCK_64X128
    return tx_depth;
}

uint8_t allowed_txt[6][TX_SIZES_ALL][TX_TYPES] = {
    {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
    //txt_th2
    {{1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    // th4
    {{1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0},
     {1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    //th_35d
    {{1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0},
     {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1},
     {1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    // th5d
    {{1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    // dct_dct and IDXT for SC
    {{1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}}};

void tx_initialize_neighbor_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                   Bool is_inter) {
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
                ? pcs_ptr
                      ->md_tx_depth_1_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx]
                : pcs_ptr->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    }
    // Set luma dc sign level coeff
    context_ptr->full_loop_luma_dc_sign_level_coeff_neighbor_array = (context_ptr->tx_depth)
        ? pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX]
                                                                        [tile_idx]
        : pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
}

void tx_update_neighbor_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                               ModeDecisionCandidateBuffer *candidate_buffer, Bool is_inter) {
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
        int8_t dc_sign_level_coeff = candidate_buffer->quantized_dc[0][context_ptr->txb_itr];
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
                              Bool is_inter, uint8_t tx_depth) {
    int      sb_size  = pcs_ptr->parent_pcs_ptr->scs_ptr->super_block_size;
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
void copy_txt_data(ModeDecisionCandidateBuffer *candidate_buffer, ModeDecisionContext *context_ptr,
                   uint32_t txb_origin_index, TxType best_tx_type) {
    uint8_t  tx_depth      = context_ptr->tx_depth;
    uint8_t  txb_itr       = context_ptr->txb_itr;
    uint32_t txb_1d_offset = context_ptr->txb_1d_offset;
    uint8_t  tx_width      = context_ptr->blk_geom->tx_width[tx_depth][txb_itr];
    uint8_t  tx_height     = context_ptr->blk_geom->tx_height[tx_depth][txb_itr];
    // copy recon_coeff_ptr
    memcpy(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_y) + txb_1d_offset,
           ((int32_t *)context_ptr->recon_coeff_ptr[best_tx_type]->buffer_y) + txb_1d_offset,
           (tx_width * tx_height * sizeof(uint32_t)));
    // copy quant_coeff_ptr
    memcpy(((int32_t *)candidate_buffer->quant_coeff_ptr->buffer_y) + txb_1d_offset,
           ((int32_t *)context_ptr->quant_coeff_ptr[best_tx_type]->buffer_y) + txb_1d_offset,
           (tx_width * tx_height * sizeof(uint32_t)));
    // copy recon_ptr
    EbPictureBufferDesc *recon_ptr = candidate_buffer->recon_ptr;
    if (context_ptr->hbd_mode_decision) {
        for (uint32_t j = 0; j < tx_height; ++j)
            memcpy(((uint16_t *)recon_ptr->buffer_y) + txb_origin_index + j * recon_ptr->stride_y,
                   ((uint16_t *)context_ptr->recon_ptr[best_tx_type]->buffer_y) + txb_origin_index +
                       j * recon_ptr->stride_y,
                   tx_width * sizeof(uint16_t));
    } else {
        for (uint32_t j = 0; j < tx_height; ++j)
            memcpy(recon_ptr->buffer_y + txb_origin_index + j * recon_ptr->stride_y,
                   context_ptr->recon_ptr[best_tx_type]->buffer_y + txb_origin_index +
                       j * recon_ptr->stride_y,
                   context_ptr->blk_geom->tx_width[tx_depth][txb_itr]);
    }
}
uint8_t get_tx_type_group(ModeDecisionContext         *context_ptr,
                          ModeDecisionCandidateBuffer *candidate_buffer, Bool only_dct_dct) {
    int tx_group = 1;
    if (!only_dct_dct) {
        if (is_intra_mode(candidate_buffer->candidate_ptr->pred_mode)) {
            tx_group =
                (context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] <
                     16 ||
                 context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] < 16)
                ? context_ptr->txt_ctrls.txt_group_intra_lt_16x16
                : context_ptr->txt_ctrls.txt_group_intra_gt_eq_16x16;
        } else {
            tx_group =
                (context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] <
                     16 ||
                 context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] < 16)
                ? context_ptr->txt_ctrls.txt_group_inter_lt_16x16
                : context_ptr->txt_ctrls.txt_group_inter_gt_eq_16x16;
        }
    }
    if (context_ptr->tx_depth == 1)
        tx_group = MAX(tx_group - context_ptr->txs_ctrls.depth1_txt_group_offset, 1);
    else if (context_ptr->tx_depth == 2)
        tx_group = MAX(tx_group - context_ptr->txs_ctrls.depth2_txt_group_offset, 1);
    return tx_group;
}
/*
 **************
*/
void perform_tx_light_pd0(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                          ModeDecisionCandidateBuffer *candidate_buffer, uint32_t qindex,
                          uint32_t *y_count_non_zero_coeffs, uint64_t *y_coeff_bits,
                          uint64_t *y_full_distortion) {
    context_ptr->three_quad_energy = 0;

    TxSize tx_size = context_ptr->blk_geom->txsize[0][0];

    if (context_ptr->md_staging_subres_step == 2) {
        if (tx_size == TX_64X64)
            tx_size = TX_64X16;
        else if (tx_size == TX_32X32)
            tx_size = TX_32X8;
        else if (tx_size == TX_16X16)
            tx_size = TX_16X4;
        else
            assert(0);
    } else if (context_ptr->md_staging_subres_step == 1) {
        if (tx_size == TX_64X64)
            tx_size = TX_64X32;
        else if (tx_size == TX_32X32)
            tx_size = TX_32X16;
        else if (tx_size == TX_16X16)
            tx_size = TX_16X8;
        else if (tx_size == TX_8X8)
            tx_size = TX_8X4;
        else
            assert(0);
    }
    assert(tx_size < TX_SIZES_ALL);
    const int32_t  tx_type          = DCT_DCT;
    const uint32_t txb_origin_index = context_ptr->blk_geom->origin_x +
        (context_ptr->blk_geom->origin_y * candidate_buffer->residual_ptr->stride_y);

    int32_t *const transf_coeff = &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[0]);
    int32_t *const recon_coeff  = &(((int32_t *)candidate_buffer->recon_coeff_ptr->buffer_y)[0]);

    EB_TRANS_COEFF_SHAPE pf_shape = context_ptr->pf_ctrls.pf_shape;

    // Y: T Q i_q
    av1_estimate_transform(
        &(((int16_t *)candidate_buffer->residual_ptr->buffer_y)[txb_origin_index]),
        candidate_buffer->residual_ptr->stride_y,
        transf_coeff,
        NOT_USED_VALUE,
        tx_size,
        &context_ptr->three_quad_energy,
        context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
        tx_type,
        PLANE_TYPE_Y,
        pf_shape);

    av1_quantize_inv_quantize_light(pcs_ptr,
                                    transf_coeff,
                                    &(((int32_t *)candidate_buffer->quant_coeff_ptr->buffer_y)[0]),
                                    recon_coeff,
                                    MIN(255, qindex + context_ptr->rate_est_ctrls.lpd0_qp_offset),
                                    tx_size,
                                    &candidate_buffer->eob[0][0],
                                    y_count_non_zero_coeffs,
                                    context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
                                    DCT_DCT);

    // LUMA DISTORTION
    uint32_t txbwidth  = context_ptr->blk_geom->tx_width[0][0];
    uint32_t txbheight = (context_ptr->blk_geom->tx_height[0][0] >>
                          context_ptr->md_staging_subres_step);
    uint32_t bwidth, bheight;
    if (pf_shape) {
        bwidth  = MAX((txbwidth >> pf_shape), 4);
        bheight = (txbheight >> pf_shape);
    } else {
        bwidth  = txbwidth < 64 ? txbwidth : 32;
        bheight = txbheight < 64 ? txbheight : 32;
    }
    picture_full_distortion32_bits_single(transf_coeff,
                                          recon_coeff,
                                          txbwidth < 64 ? txbwidth : 32,
                                          bwidth, // bwidth
                                          bheight, // bheight
                                          y_full_distortion,
                                          *y_count_non_zero_coeffs);
    y_full_distortion[DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
    const int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale_tab[tx_size]) * 2;
    y_full_distortion[DIST_CALC_RESIDUAL] =
        RIGHT_SIGNED_SHIFT(y_full_distortion[DIST_CALC_RESIDUAL], shift)
        << context_ptr->md_staging_subres_step;
    //LUMA-ONLY

    const uint32_t th = ((bwidth * bheight) >> 5);
    if (context_ptr->rate_est_ctrls.coeff_rate_est_lvl == 0) {
        uint8_t input_resolution_factor[INPUT_SIZE_COUNT] = {0, 1, 2, 3, 4, 4, 4};
        *y_coeff_bits                                     = 5000 +
            (input_resolution_factor[pcs_ptr->parent_pcs_ptr->input_resolution] * 1600) +
            (candidate_buffer->eob[0][0] * 100);
    } else if (context_ptr->rate_est_ctrls.coeff_rate_est_lvl >= 2 &&
               (candidate_buffer->eob[0][0] < th))
        *y_coeff_bits = 6000 + candidate_buffer->eob[0][0] * 500;
    else
        av1_txb_estimate_coeff_bits_light_pd0(context_ptr,
                                              candidate_buffer,
                                              context_ptr->txb_1d_offset,
                                              candidate_buffer->quant_coeff_ptr,
                                              *y_count_non_zero_coeffs,
                                              y_coeff_bits,
                                              tx_size);
    // Needed for generating recon
    candidate_buffer->y_has_coeff = (y_count_non_zero_coeffs[0] > 0);
}
// Return true if DCT_DCT is the only TX type to search
static INLINE Bool search_dct_dct_only(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                       ModeDecisionCandidateBuffer *candidate_buffer,
                                       uint8_t tx_depth, uint8_t txb_itr, uint8_t is_inter) {
    if (context_ptr->md_staging_txt_level == 0)
        return 1;

    // If previous MD stages have 0 coeffs, use DCT_DCT only
    if (context_ptr->md_stage == MD_STAGE_3 && context_ptr->use_tx_shortcuts_mds3) {
        return 1;
    } else if (context_ptr->tx_shortcut_ctrls.bypass_tx_when_zcoeff &&
               context_ptr->md_stage == MD_STAGE_3 && context_ptr->perform_mds1 &&
               !candidate_buffer->block_has_coeff) {
        return 1;
    }

    // Turn OFF TXT search for disallowed cases
    const TxSize max_tx_size = context_ptr->blk_geom->txsize[0][0];
    const TxSize tx_size     = context_ptr->blk_geom->txsize[tx_depth][txb_itr];

    // get_ext_tx_set() == 0 should correspond to a set with only DCT_DCT and there is no need to send the tx_type
    if (context_ptr->blk_geom->tx_height[tx_depth][txb_itr] > 32 ||
        context_ptr->blk_geom->tx_width[tx_depth][txb_itr] > 32 ||
        get_ext_tx_types(tx_size, is_inter, pcs_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set) == 1 ||
        get_ext_tx_set(tx_size, is_inter, pcs_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set) <= 0 ||
        (is_inter &&
         get_ext_tx_set(max_tx_size, is_inter, pcs_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set) <=
             0))
        return 1;
    return 0;
}
void tx_type_search(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                    ModeDecisionCandidateBuffer *candidate_buffer, uint32_t qindex,
                    uint8_t tx_search_skip_flag, uint32_t *y_count_non_zero_coeffs,
                    uint64_t *y_coeff_bits, uint64_t *y_full_distortion) {
    EbPictureBufferDesc *input_picture_ptr = context_ptr->hbd_mode_decision
        ? pcs_ptr->input_frame16bit
        : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    int32_t seg_qp = pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params.segmentation_enabled
        ? pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params
              .feature_data[context_ptr->blk_ptr->segment_id][SEG_LVL_ALT_Q]
        : 0;

    uint32_t   full_lambda = context_ptr->hbd_mode_decision
          ? context_ptr->full_lambda_md[EB_10_BIT_MD]
          : context_ptr->full_lambda_md[EB_8_BIT_MD];
    TxSize     tx_size = context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr];
    const Bool is_inter = (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) ||
                           candidate_buffer->candidate_ptr->use_intrabc)
        ? TRUE
        : FALSE;
    // Do not turn ON TXT search beyond this point
    const uint8_t only_dct_dct = search_dct_dct_only(pcs_ptr,
                                                     context_ptr,
                                                     candidate_buffer,
                                                     context_ptr->tx_depth,
                                                     context_ptr->txb_itr,
                                                     is_inter) ||
        tx_search_skip_flag;

    const TxSize    max_tx_size       = context_ptr->blk_geom->txsize[0][0];
    const TxSetType tx_set_type_inter = get_ext_tx_set_type(
        max_tx_size, is_inter, pcs_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set);
    const TxSetType tx_set_type = get_ext_tx_set_type(
        tx_size, is_inter, pcs_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set);

    // resize after checks on allowable TX types
    if (context_ptr->md_staging_subres_step == 2) {
        if (tx_size == TX_64X64)
            tx_size = TX_64X16;
        else if (tx_size == TX_32X32)
            tx_size = TX_32X8;
        else if (tx_size == TX_16X16)
            tx_size = TX_16X4;
        else
            assert(0);
    } else if (context_ptr->md_staging_subres_step == 1) {
        if (tx_size == TX_64X64)
            tx_size = TX_64X32;
        else if (tx_size == TX_32X32)
            tx_size = TX_32X16;
        else if (tx_size == TX_16X16)
            tx_size = TX_16X8;
        else if (tx_size == TX_8X8)
            tx_size = TX_8X4;
        else
            assert(0);
    }
    EB_TRANS_COEFF_SHAPE pf_shape = context_ptr->pf_ctrls.pf_shape;
    if (context_ptr->md_stage == MD_STAGE_3 && context_ptr->use_tx_shortcuts_mds3) {
        pf_shape = N4_SHAPE;
    }
    // only have prev. stage coeff info if mds1/2 were performed
    else if (context_ptr->tx_shortcut_ctrls.apply_pf_on_coeffs &&
             context_ptr->md_stage == MD_STAGE_3 && context_ptr->perform_mds1) {
        uint8_t use_pfn4_cond = 0;

        const uint16_t th =
            ((context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] >> 4) *
             (context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] >> 4));
        use_pfn4_cond = (candidate_buffer->count_non_zero_coeffs < th) ||
                !candidate_buffer->block_has_coeff
            ? 1
            : 0;
        if (use_pfn4_cond)
            pf_shape = N4_SHAPE;
    } else if (context_ptr->md_stage == MD_STAGE_3 && !context_ptr->perform_mds1 &&
               context_ptr->tx_shortcut_ctrls.use_neighbour_info) {
        MacroBlockD *xd = context_ptr->blk_ptr->av1xd;
        if (xd->left_available && xd->up_available && context_ptr->blk_geom->sq_size > 16 &&
            context_ptr->is_subres_safe == 1) {
            const BlockModeInfoEnc *const left_mi  = &xd->left_mbmi->block_mi;
            const BlockModeInfoEnc *const above_mi = &xd->above_mbmi->block_mi;
            if (left_mi->skip && above_mi->skip) {
                context_ptr->use_tx_shortcuts_mds3 = 1;
                pf_shape                           = N2_SHAPE;
                // PF for MVP cands
                if (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) &&
                    (candidate_buffer->candidate_ptr->pred_mode != NEWMV &&
                     candidate_buffer->candidate_ptr->pred_mode != NEW_NEWMV))
                    pf_shape = N4_SHAPE;
                // PF for NRST_NRST
                else if (candidate_buffer->candidate_ptr->pred_mode == NEAREST_NEARESTMV)
                    pf_shape = N4_SHAPE;
            }
        }
    }
    uint64_t       best_cost_tx_search = (uint64_t)~0;
    int            best_satd_tx_search = INT_MAX;
    const uint16_t satd_early_exit_th  = only_dct_dct ? 0
         : is_inter ? context_ptr->txt_ctrls.satd_early_exit_th_inter
                    : context_ptr->txt_ctrls
                         .satd_early_exit_th_intra; // only compute satd when using TXT search
    int32_t        tx_type;
    uint16_t       txb_origin_x =
        context_ptr->blk_geom->tx_org_x[is_inter][context_ptr->tx_depth][context_ptr->txb_itr];
    uint16_t txb_origin_y =
        context_ptr->blk_geom->tx_org_y[is_inter][context_ptr->tx_depth][context_ptr->txb_itr];
    uint32_t txb_origin_index = txb_origin_x +
        (txb_origin_y * candidate_buffer->residual_ptr->stride_y);
    uint32_t input_txb_origin_index = (context_ptr->sb_origin_x + txb_origin_x +
                                       input_picture_ptr->origin_x) +
        ((context_ptr->sb_origin_y + txb_origin_y + input_picture_ptr->origin_y) *
         input_picture_ptr->stride_y);
    int32_t cropped_tx_width = MIN(
        context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
        pcs_ptr->parent_pcs_ptr->aligned_width - (context_ptr->sb_origin_x + txb_origin_x));
    int32_t cropped_tx_height = MIN(
        (uint8_t)(context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] >>
                  context_ptr->md_staging_subres_step),
        pcs_ptr->parent_pcs_ptr->aligned_height - (context_ptr->sb_origin_y + txb_origin_y));
    const uint32_t txbwidth =
        context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr];
    const uint32_t txbheight =
        (context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] >>
         context_ptr->md_staging_subres_step);
    context_ptr->luma_txb_skip_context = 0;
    context_ptr->luma_dc_sign_context  = 0;
    if (context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx)
        get_txb_ctx(pcs_ptr,
                    COMPONENT_LUMA,
                    context_ptr->full_loop_luma_dc_sign_level_coeff_neighbor_array,
                    context_ptr->sb_origin_x + txb_origin_x,
                    context_ptr->sb_origin_y + txb_origin_y,
                    context_ptr->blk_geom->bsize,
                    tx_size,
                    &context_ptr->luma_txb_skip_context,
                    &context_ptr->luma_dc_sign_context);
    TxType best_tx_type = DCT_DCT;
    // local variables for all TX types
    uint16_t eob_txt[TX_TYPES]                                  = {0};
    int32_t  quantized_dc_txt[TX_TYPES]                         = {0};
    uint32_t y_count_non_zero_coeffs_txt[TX_TYPES]              = {0};
    uint64_t y_txb_coeff_bits_txt[TX_TYPES]                     = {0};
    uint64_t txb_full_distortion_txt[TX_TYPES][DIST_CALC_TOTAL] = {{0}};
    int      tx_type_tot_group = get_tx_type_group(context_ptr, candidate_buffer, only_dct_dct);
    for (int tx_type_group_idx = 0; tx_type_group_idx < tx_type_tot_group; ++tx_type_group_idx) {
        uint32_t best_tx_non_coeff = 64 * 64;
        for (int tx_type_idx = 0; tx_type_idx < TX_TYPES; ++tx_type_idx) {
            if (pcs_ptr->parent_pcs_ptr->sc_class1)
                tx_type = tx_type_group_sc[tx_type_group_idx][tx_type_idx];
            else
                tx_type = tx_type_group[tx_type_group_idx][tx_type_idx];

            if (tx_type == INVALID_TX_TYPE)
                break;

            if (only_dct_dct && tx_type != DCT_DCT)
                continue;
            if (tx_type != DCT_DCT) {
                if (is_inter) {
                    if (av1_ext_tx_used[tx_set_type_inter][tx_type] == 0)
                        continue;
                }

                if (av1_ext_tx_used[tx_set_type][tx_type] == 0)
                    continue;
            }
            // Do not use temporary buffers when TXT is OFF
            EbPictureBufferDesc *recon_coeff_ptr = (tx_type == DCT_DCT)
                ? candidate_buffer->recon_coeff_ptr
                : context_ptr->recon_coeff_ptr[tx_type];

            EbPictureBufferDesc *recon_ptr = (tx_type == DCT_DCT) ? candidate_buffer->recon_ptr
                                                                  : context_ptr->recon_ptr[tx_type];
            EbPictureBufferDesc *quant_coeff_ptr = (tx_type == DCT_DCT)
                ? candidate_buffer->quant_coeff_ptr
                : context_ptr->quant_coeff_ptr[tx_type];
            context_ptr->three_quad_energy       = 0;
            if (!tx_search_skip_flag) {
                // Y: T Q i_q
                av1_estimate_transform(
                    &(((int16_t *)candidate_buffer->residual_ptr->buffer_y)[txb_origin_index]),
                    candidate_buffer->residual_ptr->stride_y,
                    &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[context_ptr->txb_1d_offset]),
                    NOT_USED_VALUE,
                    tx_size,
                    &context_ptr->three_quad_energy,
                    context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
                    tx_type,
                    PLANE_TYPE_Y,
                    pf_shape);
                if (satd_early_exit_th) {
                    int satd = svt_aom_satd(&(((int32_t *)context_ptr->tx_coeffs
                                                   ->buffer_y)[context_ptr->txb_1d_offset]),
                                            (txbwidth * txbheight))
                        << context_ptr->md_staging_subres_step;

                    // If SATD of current type is better than the prevous best, update best, and continue evaluating tx_type
                    if (satd < best_satd_tx_search) {
                        best_satd_tx_search = satd;
                    } else {
                        // If SATD of current type is much worse than the best then stop evaluating current tx_type
                        if ((satd - best_satd_tx_search) * 100 >
                            best_satd_tx_search * satd_early_exit_th)
                            continue;
                    }
                }

                quantized_dc_txt[tx_type] = av1_quantize_inv_quantize(
                    pcs_ptr,
                    context_ptr,
                    &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[context_ptr->txb_1d_offset]),
                    NOT_USED_VALUE,
                    &(((int32_t *)quant_coeff_ptr->buffer_y)[context_ptr->txb_1d_offset]),
                    &(((int32_t *)recon_coeff_ptr->buffer_y)[context_ptr->txb_1d_offset]),
                    qindex,
                    seg_qp,
                    context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] >>
                        context_ptr->md_staging_subres_step,
                    tx_size,
                    &eob_txt[tx_type],
                    &(y_count_non_zero_coeffs_txt[tx_type]),
                    COMPONENT_LUMA,
                    context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
                    tx_type,
                    candidate_buffer,
                    context_ptr->luma_txb_skip_context,
                    context_ptr->luma_dc_sign_context,
                    candidate_buffer->candidate_ptr->pred_mode,
                    candidate_buffer->candidate_ptr->use_intrabc,
                    full_lambda,
                    FALSE);
            }
            uint32_t y_has_coeff = y_count_non_zero_coeffs_txt[tx_type] > 0;

            // tx_type not equal to DCT_DCT and no coeff is not an acceptable option in AV1.
            if (y_has_coeff == 0 && tx_type != DCT_DCT)
                continue;

            // Perform T-1 if md_staging_spatial_sse_full_loop_level or  INTRA and tx_depth > 0 or
            if (context_ptr->md_staging_spatial_sse_full_loop_level ||
                (!is_inter && candidate_buffer->candidate_ptr->tx_depth)) {
                if (y_has_coeff)
                    inv_transform_recon_wrapper(
                        candidate_buffer->prediction_ptr->buffer_y,
                        txb_origin_index,
                        candidate_buffer->prediction_ptr->stride_y,
                        recon_ptr->buffer_y,
                        txb_origin_index,
                        candidate_buffer->recon_ptr->stride_y,
                        (int32_t *)recon_coeff_ptr->buffer_y,
                        context_ptr->txb_1d_offset,
                        context_ptr->hbd_mode_decision,
                        context_ptr->blk_geom->txsize[context_ptr->tx_depth][context_ptr->txb_itr],
                        tx_type,
                        PLANE_TYPE_Y,
                        (uint32_t)eob_txt[tx_type]);
                else
                    svt_av1_picture_copy(
                        candidate_buffer->prediction_ptr,
                        txb_origin_index,
                        0,
                        recon_ptr,
                        txb_origin_index,
                        0,
                        context_ptr->blk_geom
                            ->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                        context_ptr->blk_geom
                            ->tx_height[context_ptr->tx_depth][context_ptr->txb_itr],
                        0,
                        0,
                        PICTURE_BUFFER_DESC_Y_FLAG,
                        context_ptr->hbd_mode_decision);

                EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
                    ? svt_full_distortion_kernel16_bits
                    : svt_spatial_full_distortion_kernel;
                txb_full_distortion_txt[tx_type][DIST_CALC_PREDICTION] = spatial_full_dist_type_fun(
                    input_picture_ptr->buffer_y,
                    input_txb_origin_index,
                    input_picture_ptr->stride_y,
                    candidate_buffer->prediction_ptr->buffer_y,
                    (int32_t)txb_origin_index,
                    candidate_buffer->prediction_ptr->stride_y,
                    cropped_tx_width,
                    cropped_tx_height);
                txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL] = spatial_full_dist_type_fun(
                    input_picture_ptr->buffer_y,
                    input_txb_origin_index,
                    input_picture_ptr->stride_y,
                    recon_ptr->buffer_y,
                    (int32_t)txb_origin_index,
                    candidate_buffer->recon_ptr->stride_y,
                    cropped_tx_width,
                    cropped_tx_height);
                txb_full_distortion_txt[tx_type][DIST_CALC_PREDICTION] <<= 4;
                txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL] <<= 4;
            } else {
                // LUMA DISTORTION
                uint32_t bwidth, bheight;
                if (pf_shape && !tx_search_skip_flag) {
                    bwidth  = MAX((txbwidth >> pf_shape), 4);
                    bheight = (txbheight >> pf_shape);
                } else {
                    bwidth  = txbwidth < 64 ? txbwidth : 32;
                    bheight = txbheight < 64 ? txbheight : 32;
                }
                picture_full_distortion32_bits_single(
                    &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[context_ptr->txb_1d_offset]),
                    &(((int32_t *)recon_coeff_ptr->buffer_y)[context_ptr->txb_1d_offset]),
                    txbwidth < 64 ? txbwidth : 32,
                    bwidth,
                    bheight,
                    txb_full_distortion_txt[tx_type],
                    y_count_non_zero_coeffs_txt[tx_type]);
                txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL] +=
                    context_ptr->three_quad_energy;
                txb_full_distortion_txt[tx_type][DIST_CALC_PREDICTION] +=
                    context_ptr->three_quad_energy;
                //assert(context_ptr->three_quad_energy == 0 && context_ptr->cu_stats->size < 64);
                const int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale_tab[tx_size]) * 2;
                txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(
                    txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL], shift);
                txb_full_distortion_txt[tx_type][DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(
                    txb_full_distortion_txt[tx_type][DIST_CALC_PREDICTION], shift);
            }
            txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL] =
                txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL]
                << context_ptr->md_staging_subres_step;
            txb_full_distortion_txt[tx_type][DIST_CALC_PREDICTION] =
                txb_full_distortion_txt[tx_type][DIST_CALC_PREDICTION]
                << context_ptr->md_staging_subres_step;
            // Do not perform rate estimation @ tx_type search if current tx_type dist is higher than best_cost
            uint64_t early_cost = RDCOST(
                full_lambda, 0, txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL]);
            if (early_cost > best_cost_tx_search) {
                continue;
            }
            //LUMA-ONLY
            uint64_t th =
                ((context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] *
                  context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr]) >>
                 6);
            if ((context_ptr->rate_est_ctrls.coeff_rate_est_lvl >= 2 ||
                 context_ptr->rate_est_ctrls.coeff_rate_est_lvl == 0) &&
                (y_count_non_zero_coeffs_txt[tx_type] < (th)))
                y_txb_coeff_bits_txt[tx_type] = 6000 + y_count_non_zero_coeffs_txt[tx_type] * 1000;
            else if (context_ptr->rate_est_ctrls.coeff_rate_est_lvl == 0)
                y_txb_coeff_bits_txt[tx_type] = 3000 + y_count_non_zero_coeffs_txt[tx_type] * 100;
            else
                av1_txb_estimate_coeff_bits(
                    context_ptr,
                    0, //allow_update_cdf,
                    NULL, //FRAME_CONTEXT *ec_ctx,
                    pcs_ptr,
                    candidate_buffer,
                    context_ptr->txb_1d_offset,
                    0,
                    quant_coeff_ptr,
                    y_count_non_zero_coeffs_txt[tx_type],
                    0,
                    0,
                    &(y_txb_coeff_bits_txt[tx_type]),
                    &(y_txb_coeff_bits_txt[tx_type]),
                    &(y_txb_coeff_bits_txt[tx_type]),
                    tx_size,
                    context_ptr->blk_geom->txsize_uv[context_ptr->tx_depth][context_ptr->txb_itr],
                    tx_type,
                    NOT_USED_VALUE,
                    COMPONENT_LUMA);
            uint64_t cost = RDCOST(full_lambda,
                                   y_txb_coeff_bits_txt[tx_type],
                                   txb_full_distortion_txt[tx_type][DIST_CALC_RESIDUAL]);
            if (cost < best_cost_tx_search) {
                best_cost_tx_search = cost;
                best_tx_type        = tx_type;
                best_tx_non_coeff   = y_count_non_zero_coeffs_txt[tx_type];
            }

            // Skip remaining TX types based on absolute cost TH and absolute # of coeffs TH
            if (!only_dct_dct) {
                uint32_t coeff_th      = context_ptr->txt_ctrls.early_exit_coeff_th;
                uint32_t dist_err_unit = context_ptr->txt_ctrls.early_exit_dist_th;
                uint32_t dist_err =
                    context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] *
                    context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] *
                    dist_err_unit;

                uint64_t cost_th = dist_err_unit ? RDCOST(full_lambda, 1, dist_err)
                                                 : 0; // if distortion th=0, set cost_th to 0

                if (best_tx_non_coeff < coeff_th || best_cost_tx_search < cost_th) {
                    tx_type_idx       = TX_TYPES;
                    tx_type_group_idx = tx_type_tot_group;
                }
            }
        }
    }
    //  Best Tx Type Pass
    candidate_buffer->candidate_ptr->transform_type[context_ptr->txb_itr] = best_tx_type;
    // update with best_tx_type data
    (*y_coeff_bits) += y_txb_coeff_bits_txt[best_tx_type];
    y_full_distortion[DIST_CALC_RESIDUAL] +=
        txb_full_distortion_txt[best_tx_type][DIST_CALC_RESIDUAL];
    y_full_distortion[DIST_CALC_PREDICTION] +=
        txb_full_distortion_txt[best_tx_type][DIST_CALC_PREDICTION];

    y_count_non_zero_coeffs[context_ptr->txb_itr] = y_count_non_zero_coeffs_txt[best_tx_type];
    candidate_buffer->y_has_coeff |= ((y_count_non_zero_coeffs_txt[best_tx_type] > 0)
                                      << context_ptr->txb_itr);
    candidate_buffer->quantized_dc[0][context_ptr->txb_itr] = quantized_dc_txt[best_tx_type];
    candidate_buffer->eob[0][context_ptr->txb_itr]          = eob_txt[best_tx_type];
    // Do not copy when TXT is OFF
    // Data is already in candidate_buffer
    if (best_tx_type != DCT_DCT) {
        // copy best_tx_type data
        copy_txt_data(candidate_buffer, context_ptr, txb_origin_index, best_tx_type);
    }
    context_ptr->txb_1d_offset +=
        context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr] *
        (context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] >>
         context_ptr->md_staging_subres_step);
    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter && context_ptr->txb_itr == 0)
        candidate_buffer->candidate_ptr->transform_type_uv =
            candidate_buffer->candidate_ptr->transform_type[context_ptr->txb_itr];
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
        max_blocks_wide += gcc_right_shift(xd->mb_to_right_edge, 3 + !!plane);

    // Scale the width in the transform block unit.
    return max_blocks_wide >> tx_size_wide_log2[0];
}

static INLINE int max_block_high(const MacroBlockD *xd, BlockSize bsize, int plane) {
    int max_blocks_high = block_size_high[bsize];

    if (xd->mb_to_bottom_edge < 0)
        max_blocks_high += gcc_right_shift(xd->mb_to_bottom_edge, 3 + !!plane);

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
    if (tx_size == TX_4X4)
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
    const int write_txfm_partition =
        (tx_size == tx_depth_to_tx_size[mbmi->block_mi.tx_depth][mbmi->block_mi.sb_type]);
    if (write_txfm_partition) {
        bits += md_rate_estimation_ptr->txfm_partition_fac_bits[ctx][0];

        txfm_partition_update(
            xd->above_txfm_context + blk_col, xd->left_txfm_context + blk_row, tx_size, tx_size);

    } else {
        assert(tx_size < TX_SIZES_ALL);
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

int is_inter_block(const BlockModeInfoEnc *bloc_mi);

// Returns a context number for the given MB prediction signal
// The mode info data structure has a one element border above and to the
// left of the entries corresponding to real blocks.
// The prediction flags in these dummy entries are initialized to 0.
static INLINE int get_tx_size_context(const MacroBlockD *xd) {
    const ModeInfo         *mi          = xd->mi[0];
    const MbModeInfo       *mbmi        = &mi->mbmi;
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
static uint64_t cost_selected_tx_size(const MacroBlockD       *xd,
                                      MdRateEstimationContext *md_rate_estimation_ptr,
                                      TxSize                   tx_size) {
    const ModeInfo *const   mi    = xd->mi[0];
    const MbModeInfo *const mbmi  = &mi->mbmi;
    const BlockSize         bsize = mbmi->block_mi.sb_type;
    uint64_t                bits  = 0;

    if (block_signals_txsize(bsize)) {
        const int tx_size_ctx = get_tx_size_context(xd);
        assert(bsize < BlockSizeS_ALL);
        const int     depth       = tx_size_to_depth(tx_size, bsize);
        const int32_t tx_size_cat = bsize_to_tx_size_cat(bsize);
        bits += md_rate_estimation_ptr->tx_size_fac_bits[tx_size_cat][tx_size_ctx][depth];
    }

    return bits;
}
static uint64_t tx_size_bits(MdRateEstimationContext *md_rate_estimation_ptr, MacroBlockD *xd,
                             const MbModeInfo *mbmi, TxSize tx_size, TxMode tx_mode,
                             BlockSize bsize, uint8_t skip) {
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
            bits += cost_selected_tx_size(xd, md_rate_estimation_ptr, tx_size);
            set_txfm_ctxs(tx_size, xd->n8_w, xd->n8_h, 0, xd);
        }
    } else {
        set_txfm_ctxs(tx_size, xd->n8_w, xd->n8_h, skip && is_inter_block(&mbmi->block_mi), xd);
    }
    return bits;
}

void set_mi_row_col(PictureControlSet *pcs_ptr, MacroBlockD *xd, TileInfo *tile, int mi_row, int bh,
                    int mi_col, int bw, uint32_t mi_stride, int mi_rows, int mi_cols);

uint64_t estimate_tx_size_bits(PictureControlSet *pcsPtr, ModeDecisionContext *context_ptr,
                               ModeDecisionCandidate *candidate_ptr, Bool skip_flag,
                               uint32_t blk_origin_x, uint32_t blk_origin_y, BlkStruct *blk_ptr,
                               const BlockGeom *blk_geom, NeighborArrayUnit *txfm_context_array,
                               uint8_t tx_depth, MdRateEstimationContext *md_rate_estimation_ptr) {
    uint32_t txfm_context_left_index  = get_neighbor_array_unit_left_index(txfm_context_array,
                                                                          blk_origin_y);
    uint32_t txfm_context_above_index = get_neighbor_array_unit_top_index(txfm_context_array,
                                                                          blk_origin_x);

    TxMode       tx_mode = pcsPtr->parent_pcs_ptr->frm_hdr.tx_mode;
    MacroBlockD *xd      = blk_ptr->av1xd;
    BlockSize    bsize   = blk_geom->bsize;
    MbModeInfo  *mbmi    = &xd->mi[0]->mbmi;

    svt_memcpy(context_ptr->above_txfm_context,
               &(txfm_context_array->top_array[txfm_context_above_index]),
               (blk_geom->bwidth >> MI_SIZE_LOG2) * sizeof(TXFM_CONTEXT));
    svt_memcpy(context_ptr->left_txfm_context,
               &(txfm_context_array->left_array[txfm_context_left_index]),
               (blk_geom->bheight >> MI_SIZE_LOG2) * sizeof(TXFM_CONTEXT));

    xd->above_txfm_context      = context_ptr->above_txfm_context;
    xd->left_txfm_context       = context_ptr->left_txfm_context;
    mbmi->block_mi.sb_type      = blk_geom->bsize;
    mbmi->block_mi.use_intrabc  = candidate_ptr->use_intrabc;
    mbmi->block_mi.ref_frame[0] = candidate_ptr->ref_frame_type;
    mbmi->block_mi.tx_depth     = tx_depth;
    uint64_t bits               = tx_size_bits(
        md_rate_estimation_ptr, xd, mbmi, blk_geom->txsize[tx_depth][0], tx_mode, bsize, skip_flag);
    return bits;
}

uint64_t get_tx_size_bits(ModeDecisionCandidateBuffer *candidateBuffer,
                          ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                          uint8_t tx_depth, Bool block_has_coeff) {
    return estimate_tx_size_bits(pcs_ptr,
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
}

void init_tx_candidate_buffer(ModeDecisionCandidateBuffer *candidate_buffer,
                              ModeDecisionContext *context_ptr, uint8_t end_tx_depth) {
    uint32_t block_index = context_ptr->blk_geom->origin_x +
        (context_ptr->blk_geom->origin_y * context_ptr->sb_size);
    if (end_tx_depth) {
        svt_memcpy(context_ptr->candidate_buffer_tx_depth_1->candidate_ptr,
                   candidate_buffer->candidate_ptr,
                   sizeof(ModeDecisionCandidate));
        context_ptr->candidate_buffer_tx_depth_1->block_has_coeff =
            candidate_buffer->block_has_coeff;
        if (context_ptr->hbd_mode_decision) {
            // Copy pred to tx_depth_1 candidate_buffer
            {
                uint16_t *src = &(
                    ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y)[block_index]);
                uint16_t *dst = &(((uint16_t *)context_ptr->candidate_buffer_tx_depth_1
                                       ->prediction_ptr->buffer_y)[block_index]);
                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    svt_memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
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
                    svt_memcpy(dst, src, context_ptr->blk_geom->bwidth << 1);
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
                    svt_memcpy(dst, src, context_ptr->blk_geom->bwidth);
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
                    svt_memcpy(dst, src, context_ptr->blk_geom->bwidth << 1);
                    src += candidate_buffer->residual_ptr->stride_y;
                    dst += context_ptr->candidate_buffer_tx_depth_1->residual_ptr->stride_y;
                }
            }
        }
    }
    if (end_tx_depth == 2) {
        svt_memcpy(context_ptr->candidate_buffer_tx_depth_2->candidate_ptr,
                   candidate_buffer->candidate_ptr,
                   sizeof(ModeDecisionCandidate));

        context_ptr->candidate_buffer_tx_depth_2->block_has_coeff =
            candidate_buffer->block_has_coeff;
        if (context_ptr->hbd_mode_decision) {
            // Copy pred to tx_depth_1 candidate_buffer
            {
                uint16_t *src = &(
                    ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y)[block_index]);
                uint16_t *dst = &(((uint16_t *)context_ptr->candidate_buffer_tx_depth_2
                                       ->prediction_ptr->buffer_y)[block_index]);

                for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                    svt_memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
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
                    svt_memcpy(dst, src, context_ptr->blk_geom->bwidth << 1);
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
                    svt_memcpy(dst, src, context_ptr->blk_geom->bwidth);
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
                    svt_memcpy(dst, src, context_ptr->blk_geom->bwidth << 1);
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
        (context_ptr->blk_geom->origin_y * context_ptr->sb_size);
    if (best_tx_depth == 1) {
        // Copy depth 1 mode/type/eob ..
        svt_memcpy(candidate_buffer->candidate_ptr,
                   context_ptr->candidate_buffer_tx_depth_1->candidate_ptr,
                   sizeof(ModeDecisionCandidate));
        svt_memcpy(candidate_buffer->eob,
                   context_ptr->candidate_buffer_tx_depth_1->eob,
                   sizeof(uint16_t) * 1 /*copy luma*/ * MAX_TXB_COUNT);
        svt_memcpy(candidate_buffer->quantized_dc,
                   context_ptr->candidate_buffer_tx_depth_1->quantized_dc,
                   sizeof(int32_t) * 1 /*copy luma*/ * MAX_TXB_COUNT);

        candidate_buffer->y_has_coeff = context_ptr->candidate_buffer_tx_depth_1->y_has_coeff;
        // Copy depth 1 pred
        if (context_ptr->hbd_mode_decision) {
            uint16_t *src = &(((uint16_t *)context_ptr->candidate_buffer_tx_depth_1->prediction_ptr
                                   ->buffer_y)[block_index]);
            uint16_t *dst = &(
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y)[block_index]);
            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                svt_memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
                src += context_ptr->candidate_buffer_tx_depth_1->prediction_ptr->stride_y;
                dst += candidate_buffer->prediction_ptr->stride_y;
            }
        } else {
            EbByte src = &(
                context_ptr->candidate_buffer_tx_depth_1->prediction_ptr->buffer_y[block_index]);
            EbByte dst = &(candidate_buffer->prediction_ptr->buffer_y[block_index]);
            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                svt_memcpy(dst, src, context_ptr->blk_geom->bwidth);
                src += context_ptr->candidate_buffer_tx_depth_1->prediction_ptr->stride_y;
                dst += candidate_buffer->prediction_ptr->stride_y;
            }
        }
        // Copy depth 1 recon coeff
        svt_memcpy(candidate_buffer->recon_coeff_ptr->buffer_y,
                   context_ptr->candidate_buffer_tx_depth_1->recon_coeff_ptr->buffer_y,
                   (context_ptr->blk_geom->bwidth * context_ptr->blk_geom->bheight << 2));
        svt_memcpy(candidate_buffer->quant_coeff_ptr->buffer_y,
                   context_ptr->candidate_buffer_tx_depth_1->quant_coeff_ptr->buffer_y,
                   (context_ptr->blk_geom->bwidth * context_ptr->blk_geom->bheight << 2));
    }
    if (best_tx_depth == 2) {
        // Copy depth 2 mode/type/eob ..
        svt_memcpy(candidate_buffer->candidate_ptr,
                   context_ptr->candidate_buffer_tx_depth_2->candidate_ptr,
                   sizeof(ModeDecisionCandidate));
        svt_memcpy(candidate_buffer->eob,
                   context_ptr->candidate_buffer_tx_depth_2->eob,
                   sizeof(uint16_t) * 1 /*copy luma*/ * MAX_TXB_COUNT);
        svt_memcpy(candidate_buffer->quantized_dc,
                   context_ptr->candidate_buffer_tx_depth_2->quantized_dc,
                   sizeof(int32_t) * 1 /*copy luma*/ * MAX_TXB_COUNT);

        candidate_buffer->y_has_coeff = context_ptr->candidate_buffer_tx_depth_2->y_has_coeff;
        // Copy depth 2 pred
        if (context_ptr->hbd_mode_decision) {
            uint16_t *src = &(((uint16_t *)context_ptr->candidate_buffer_tx_depth_2->prediction_ptr
                                   ->buffer_y)[block_index]);
            uint16_t *dst = &(
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y)[block_index]);
            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                svt_memcpy(dst, src, context_ptr->blk_geom->bwidth * sizeof(uint16_t));
                src += context_ptr->candidate_buffer_tx_depth_2->prediction_ptr->stride_y;
                dst += candidate_buffer->prediction_ptr->stride_y;
            }
        } else {
            EbByte src = &(
                context_ptr->candidate_buffer_tx_depth_2->prediction_ptr->buffer_y[block_index]);
            EbByte dst = &(candidate_buffer->prediction_ptr->buffer_y[block_index]);
            for (int i = 0; i < context_ptr->blk_geom->bheight; i++) {
                svt_memcpy(dst, src, context_ptr->blk_geom->bwidth);
                src += context_ptr->candidate_buffer_tx_depth_2->prediction_ptr->stride_y;
                dst += candidate_buffer->prediction_ptr->stride_y;
            }
        }
        // Copy depth 2 recon coeff
        svt_memcpy(candidate_buffer->recon_coeff_ptr->buffer_y,
                   context_ptr->candidate_buffer_tx_depth_2->recon_coeff_ptr->buffer_y,
                   (context_ptr->blk_geom->bwidth * context_ptr->blk_geom->bheight << 2));
        svt_memcpy(candidate_buffer->quant_coeff_ptr->buffer_y,
                   context_ptr->candidate_buffer_tx_depth_2->quant_coeff_ptr->buffer_y,
                   (context_ptr->blk_geom->bwidth * context_ptr->blk_geom->bheight << 2));
    }
}
void perform_tx_partitioning(ModeDecisionCandidateBuffer *candidate_buffer,
                             ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                             uint8_t start_tx_depth, uint8_t end_tx_depth, uint32_t qindex,
                             uint32_t *y_count_non_zero_coeffs, uint64_t *y_coeff_bits,
                             uint64_t *y_full_distortion) {
    uint32_t             full_lambda       = context_ptr->hbd_mode_decision
                          ? context_ptr->full_lambda_md[EB_10_BIT_MD]
                          : context_ptr->full_lambda_md[EB_8_BIT_MD];
    EbPictureBufferDesc *input_picture_ptr = context_ptr->hbd_mode_decision
        ? pcs_ptr->input_frame16bit
        : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    const Bool           is_inter = (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) ||
                           candidate_buffer->candidate_ptr->use_intrabc)
                  ? TRUE
                  : FALSE;
    uint8_t              best_tx_depth     = 0;
    uint64_t             best_cost_search  = (uint64_t)~0;
    uint8_t              is_best_has_coeff = 1;
    if (end_tx_depth)
        init_tx_candidate_buffer(candidate_buffer, context_ptr, end_tx_depth);
    // Transform Depth Loop
    for (context_ptr->tx_depth = start_tx_depth; context_ptr->tx_depth <= end_tx_depth;
         context_ptr->tx_depth++) {
        if (context_ptr->txs_ctrls.prev_depth_coeff_exit) {
            if (!is_best_has_coeff)
                continue;
        }
        if (context_ptr->tx_depth)
            tx_reset_neighbor_arrays(pcs_ptr, context_ptr, is_inter, context_ptr->tx_depth);
        ModeDecisionCandidateBuffer *tx_candidate_buffer = (context_ptr->tx_depth == 0)
            ? candidate_buffer
            : (context_ptr->tx_depth == 1) ? context_ptr->candidate_buffer_tx_depth_1
                                           : context_ptr->candidate_buffer_tx_depth_2;
        tx_candidate_buffer->candidate_ptr->tx_depth     = context_ptr->tx_depth;
        if (context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx || !is_inter)
            tx_initialize_neighbor_arrays(pcs_ptr, context_ptr, is_inter);

        // Initialize TU Split
        uint32_t tx_y_count_non_zero_coeffs[MAX_NUM_OF_TU_PER_CU];
        uint64_t tx_y_coeff_bits                       = 0;
        uint64_t tx_y_full_distortion[DIST_CALC_TOTAL] = {0};

        context_ptr->txb_1d_offset       = 0;
        tx_candidate_buffer->y_has_coeff = 0;

        uint16_t txb_count = context_ptr->blk_geom->txb_count[context_ptr->tx_depth];

        uint32_t block_has_coeff = FALSE;
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
                // This check assumes no txs search @ a previous md_stage()
                if (context_ptr->tx_depth)
                    av1_intra_luma_prediction(context_ptr, pcs_ptr, tx_candidate_buffer);

                // Y Residual
                residual_kernel(
                    input_picture_ptr->buffer_y,
                    input_txb_origin_index,
                    input_picture_ptr->stride_y << context_ptr->md_staging_subres_step,
                    tx_candidate_buffer->prediction_ptr->buffer_y,
                    txb_origin_index,
                    tx_candidate_buffer->prediction_ptr->stride_y
                        << context_ptr->md_staging_subres_step,
                    (int16_t *)tx_candidate_buffer->residual_ptr->buffer_y,
                    txb_origin_index,
                    tx_candidate_buffer->residual_ptr->stride_y,
                    context_ptr->hbd_mode_decision,
                    context_ptr->blk_geom->tx_width[context_ptr->tx_depth][context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height[context_ptr->tx_depth][context_ptr->txb_itr] >>
                        context_ptr->md_staging_subres_step);
            }
            uint8_t tx_search_skip_flag = 0;
            // only have prev. stage coeff info if mds1/2 were performed
            if (context_ptr->tx_shortcut_ctrls.bypass_tx_when_zcoeff &&
                context_ptr->md_stage == MD_STAGE_3 && context_ptr->perform_mds1 &&
                !candidate_buffer->block_has_coeff)
                tx_search_skip_flag = 1;
            tx_type_search(pcs_ptr,
                           context_ptr,
                           tx_candidate_buffer,
                           qindex,
                           tx_search_skip_flag,
                           &(tx_y_count_non_zero_coeffs[0]),
                           &tx_y_coeff_bits,
                           &tx_y_full_distortion[0]);

            uint32_t y_has_coeff = tx_y_count_non_zero_coeffs[context_ptr->txb_itr] > 0;

            if (context_ptr->tx_depth)
                tx_update_neighbor_arrays(pcs_ptr, context_ptr, tx_candidate_buffer, is_inter);

            if (y_has_coeff)
                block_has_coeff = TRUE;

            uint64_t current_tx_cost = RDCOST(
                full_lambda, tx_y_coeff_bits, tx_y_full_distortion[DIST_CALC_RESIDUAL]);
            if (current_tx_cost > best_cost_search)
                break;

        } // Transform Loop

        if (end_tx_depth) {
            const uint64_t tx_size_bit = pcs_ptr->parent_pcs_ptr->frm_hdr.tx_mode == TX_MODE_SELECT
                ? get_tx_size_bits(tx_candidate_buffer,
                                   context_ptr,
                                   pcs_ptr,
                                   context_ptr->tx_depth,
                                   block_has_coeff)
                : 0;

            const uint64_t cost = RDCOST(full_lambda,
                                         tx_y_coeff_bits + tx_size_bit,
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

    if (best_tx_depth)
        update_tx_candidate_buffer(candidate_buffer, context_ptr, best_tx_depth);
}
/*
   DCT_DCT path for light PD1
*/
void perform_dct_dct_tx_light_pd1(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                  ModeDecisionCandidateBuffer *candidate_buffer, BlockLocation *loc,
                                  uint64_t *y_coeff_bits, uint64_t *y_full_distortion) {
    uint32_t             full_lambda       = context_ptr->hbd_mode_decision
                          ? context_ptr->full_lambda_md[EB_10_BIT_MD]
                          : context_ptr->full_lambda_md[EB_8_BIT_MD];
    EbPictureBufferDesc *input_picture_ptr = context_ptr->hbd_mode_decision
        ? pcs_ptr->input_frame16bit
        : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    const Bool is_inter = is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) ? TRUE : FALSE;
    context_ptr->three_quad_energy = 0;
    residual_kernel(input_picture_ptr->buffer_y,
                    loc->input_origin_index,
                    input_picture_ptr->stride_y,
                    candidate_buffer->prediction_ptr->buffer_y,
                    loc->blk_origin_index,
                    candidate_buffer->prediction_ptr->stride_y,
                    (int16_t *)candidate_buffer->residual_ptr->buffer_y,
                    loc->blk_origin_index,
                    candidate_buffer->residual_ptr->stride_y,
                    context_ptr->hbd_mode_decision,
                    context_ptr->blk_geom->bwidth,
                    context_ptr->blk_geom->bheight);
    TxSize tx_size = context_ptr->blk_geom->txsize[0][0];
    assert(tx_size < TX_SIZES_ALL);
    EB_TRANS_COEFF_SHAPE pf_shape = context_ptr->pf_ctrls.pf_shape;
    if (context_ptr->use_tx_shortcuts_mds3) {
        pf_shape = N4_SHAPE;
    } else if (context_ptr->lpd1_tx_ctrls.use_neighbour_info) {
        MacroBlockD *xd = context_ptr->blk_ptr->av1xd;
        if (xd->left_available && xd->up_available && context_ptr->blk_geom->sq_size > 16 &&
            context_ptr->is_subres_safe == 1) {
            const BlockModeInfoEnc *const left_mi  = &xd->left_mbmi->block_mi;
            const BlockModeInfoEnc *const above_mi = &xd->above_mbmi->block_mi;
            if (left_mi->skip && above_mi->skip) {
                context_ptr->use_tx_shortcuts_mds3 = 1;
                pf_shape                           = N2_SHAPE;
                // PF for MVP cands
                if (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) &&
                    (candidate_buffer->candidate_ptr->pred_mode != NEWMV &&
                     candidate_buffer->candidate_ptr->pred_mode != NEW_NEWMV))
                    pf_shape = N4_SHAPE;
                // PF for NRST_NRST
                else if (candidate_buffer->candidate_ptr->pred_mode == NEAREST_NEARESTMV)
                    pf_shape = N4_SHAPE;
            }
        }
    }
    // local variables for performing TX
    uint32_t y_count_non_zero_coeffs_txt;

    EbPictureBufferDesc *const recon_coeff_ptr = candidate_buffer->recon_coeff_ptr;
    EbPictureBufferDesc *const quant_coeff_ptr = candidate_buffer->quant_coeff_ptr;

    // Y: T Q i_q
    av1_estimate_transform(
        &(((int16_t *)candidate_buffer->residual_ptr->buffer_y)[loc->blk_origin_index]),
        candidate_buffer->residual_ptr->stride_y,
        &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[0]),
        NOT_USED_VALUE,
        tx_size,
        &context_ptr->three_quad_energy,
        context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
        DCT_DCT,
        PLANE_TYPE_Y,
        pf_shape);
    candidate_buffer->quantized_dc[0][0] = av1_quantize_inv_quantize(
        pcs_ptr,
        context_ptr,
        &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[0]),
        NOT_USED_VALUE,
        &(((int32_t *)quant_coeff_ptr->buffer_y)[0]),
        &(((int32_t *)recon_coeff_ptr->buffer_y)[0]),
        context_ptr->blk_ptr->qindex,
        0,
        context_ptr->blk_geom->tx_width[0][0],
        context_ptr->blk_geom->tx_height[0][0],
        tx_size,
        &(candidate_buffer->eob[0][0]),
        &(y_count_non_zero_coeffs_txt),
        COMPONENT_LUMA,
        context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
        DCT_DCT,
        candidate_buffer,
        0,
        0,
        candidate_buffer->candidate_ptr->pred_mode,
        0,
        full_lambda,
        FALSE);
    // LUMA DISTORTION
    const uint32_t txbwidth  = context_ptr->blk_geom->tx_width[0][0];
    const uint32_t txbheight = context_ptr->blk_geom->tx_height[0][0];
    uint32_t       bwidth, bheight;
    if (pf_shape) {
        bwidth  = MAX((txbwidth >> pf_shape), 4);
        bheight = (txbheight >> pf_shape);
    } else {
        bwidth  = txbwidth < 64 ? txbwidth : 32;
        bheight = txbheight < 64 ? txbheight : 32;
    }

    picture_full_distortion32_bits_single(&(((int32_t *)context_ptr->tx_coeffs->buffer_y)[0]),
                                          &(((int32_t *)recon_coeff_ptr->buffer_y)[0]),
                                          txbwidth < 64 ? txbwidth : 32,
                                          bwidth,
                                          bheight,
                                          y_full_distortion,
                                          candidate_buffer->eob[0][0]);
    const int32_t shift                   = (MAX_TX_SCALE - av1_get_tx_scale_tab[tx_size]) * 2;
    y_full_distortion[DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(
        y_full_distortion[DIST_CALC_RESIDUAL] + context_ptr->three_quad_energy, shift);
    y_full_distortion[DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(
        y_full_distortion[DIST_CALC_PREDICTION] + context_ptr->three_quad_energy, shift);
    //LUMA-ONLY
    const uint32_t th = ((txbwidth * txbheight) >> 6);
    if ((context_ptr->rate_est_ctrls.coeff_rate_est_lvl >= 2 ||
         context_ptr->rate_est_ctrls.coeff_rate_est_lvl == 0) &&
        candidate_buffer->eob[0][0] < th)
        *y_coeff_bits = 6000 + candidate_buffer->eob[0][0] * 1000;
    else if (context_ptr->rate_est_ctrls.coeff_rate_est_lvl == 0)
        *y_coeff_bits = 6000 + candidate_buffer->eob[0][0] * 400;
    else
        av1_txb_estimate_coeff_bits(context_ptr,
                                    0,
                                    NULL,
                                    pcs_ptr,
                                    candidate_buffer,
                                    0,
                                    0,
                                    quant_coeff_ptr,
                                    candidate_buffer->eob[0][0],
                                    0,
                                    0,
                                    y_coeff_bits,
                                    NULL,
                                    NULL,
                                    tx_size,
                                    context_ptr->blk_geom->txsize_uv[0][0],
                                    DCT_DCT,
                                    NOT_USED_VALUE,
                                    COMPONENT_LUMA);
    //Update with best_tx_type data
    candidate_buffer->candidate_ptr->transform_type[0] = DCT_DCT;
    candidate_buffer->y_has_coeff                      = (candidate_buffer->eob[0][0] > 0);
    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
        candidate_buffer->candidate_ptr->transform_type_uv = DCT_DCT;
}
// TX path when TXT and TXS are off
void perform_dct_dct_tx(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                        ModeDecisionCandidateBuffer *candidate_buffer, uint8_t tx_search_skip_flag,
                        uint32_t qindex, uint32_t *y_count_non_zero_coeffs, uint64_t *y_coeff_bits,
                        uint64_t *y_full_distortion) {
    const uint32_t             full_lambda       = context_ptr->hbd_mode_decision
                          ? context_ptr->full_lambda_md[EB_10_BIT_MD]
                          : context_ptr->full_lambda_md[EB_8_BIT_MD];
    EbPictureBufferDesc *const input_picture_ptr = context_ptr->hbd_mode_decision
        ? pcs_ptr->input_frame16bit
        : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    const Bool is_inter            = (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) ||
                           candidate_buffer->candidate_ptr->use_intrabc)
                   ? TRUE
                   : FALSE;
    context_ptr->tx_depth          = 0;
    context_ptr->txb_itr           = 0;
    context_ptr->txb_1d_offset     = 0;
    context_ptr->three_quad_energy = 0;
    const int tx_depth             = 0;
    const int txb_itr              = 0;
    const int txb_1d_offset        = 0;
    const int tx_type              = DCT_DCT;

    const uint16_t tx_org_x         = context_ptr->blk_geom->tx_org_x[is_inter][tx_depth][txb_itr];
    const uint16_t tx_org_y         = context_ptr->blk_geom->tx_org_y[is_inter][tx_depth][txb_itr];
    const uint32_t txb_origin_index = tx_org_x +
        (tx_org_y * candidate_buffer->residual_ptr->stride_y);
    const uint32_t input_txb_origin_index = (context_ptr->sb_origin_x + tx_org_x +
                                             input_picture_ptr->origin_x) +
        ((context_ptr->sb_origin_y + tx_org_y + input_picture_ptr->origin_y) *
         input_picture_ptr->stride_y);

    // Y Residual
    if (!is_inter) {
        residual_kernel(
            input_picture_ptr->buffer_y,
            input_txb_origin_index,
            input_picture_ptr->stride_y << context_ptr->md_staging_subres_step,
            candidate_buffer->prediction_ptr->buffer_y,
            txb_origin_index,
            candidate_buffer->prediction_ptr->stride_y << context_ptr->md_staging_subres_step,
            (int16_t *)candidate_buffer->residual_ptr->buffer_y,
            txb_origin_index,
            candidate_buffer->residual_ptr->stride_y,
            context_ptr->hbd_mode_decision,
            context_ptr->blk_geom->tx_width[tx_depth][txb_itr],
            context_ptr->blk_geom->tx_height[tx_depth][txb_itr] >>
                context_ptr->md_staging_subres_step);
    }

    // TX search

    const int seg_qp  = pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params.segmentation_enabled
         ? pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params
              .feature_data[context_ptr->blk_ptr->segment_id][SEG_LVL_ALT_Q]
         : 0;
    TxSize    tx_size = context_ptr->blk_geom->txsize[tx_depth][txb_itr];

    if (context_ptr->md_staging_subres_step == 2) {
        if (tx_size == TX_64X64)
            tx_size = TX_64X16;
        else if (tx_size == TX_32X32)
            tx_size = TX_32X8;
        else if (tx_size == TX_16X16)
            tx_size = TX_16X4;
        else
            assert(0);
    } else if (context_ptr->md_staging_subres_step == 1) {
        if (tx_size == TX_64X64)
            tx_size = TX_64X32;
        else if (tx_size == TX_32X32)
            tx_size = TX_32X16;
        else if (tx_size == TX_16X16)
            tx_size = TX_16X8;
        else if (tx_size == TX_8X8)
            tx_size = TX_8X4;
        else
            assert(0);
    }
    assert(tx_size < TX_SIZES_ALL);
    EB_TRANS_COEFF_SHAPE pf_shape = context_ptr->pf_ctrls.pf_shape;
    if (context_ptr->md_stage == MD_STAGE_3 && context_ptr->use_tx_shortcuts_mds3) {
        pf_shape = N4_SHAPE;
    }
    // only have prev. stage coeff info if mds1/2 were performed
    else if (context_ptr->tx_shortcut_ctrls.apply_pf_on_coeffs &&
             context_ptr->md_stage == MD_STAGE_3 && context_ptr->perform_mds1) {
        uint8_t use_pfn4_cond = 0;

        const uint16_t th = ((context_ptr->blk_geom->tx_width[tx_depth][txb_itr] >> 4) *
                             (context_ptr->blk_geom->tx_height[tx_depth][txb_itr] >> 4));
        use_pfn4_cond     = (candidate_buffer->count_non_zero_coeffs < th) ||
                !candidate_buffer->block_has_coeff
                ? 1
                : 0;
        if (use_pfn4_cond)
            pf_shape = N4_SHAPE;
    }
    context_ptr->luma_txb_skip_context = 0;
    context_ptr->luma_dc_sign_context  = 0;
    if (context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx) {
        // Arrays updated here only necessary if DC sign context update is enabled,
        // or for intra_luma_preiction when TXS is on. TXS is assumed off in this path
        // so only update if DC sign array is needed
        tx_initialize_neighbor_arrays(pcs_ptr, context_ptr, is_inter);
        get_txb_ctx(pcs_ptr,
                    COMPONENT_LUMA,
                    context_ptr->full_loop_luma_dc_sign_level_coeff_neighbor_array,
                    context_ptr->sb_origin_x + tx_org_x,
                    context_ptr->sb_origin_y + tx_org_y,
                    context_ptr->blk_geom->bsize,
                    tx_size,
                    &context_ptr->luma_txb_skip_context,
                    &context_ptr->luma_dc_sign_context);
    }

    // local variables for performing TX
    uint32_t y_count_non_zero_coeffs_txt = 0;

    EbPictureBufferDesc *const recon_coeff_ptr = candidate_buffer->recon_coeff_ptr;
    EbPictureBufferDesc *const recon_ptr       = candidate_buffer->recon_ptr;
    EbPictureBufferDesc *const quant_coeff_ptr = candidate_buffer->quant_coeff_ptr;

    if (!tx_search_skip_flag) {
        // Y: T Q i_q
        av1_estimate_transform(
            &(((int16_t *)candidate_buffer->residual_ptr->buffer_y)[txb_origin_index]),
            candidate_buffer->residual_ptr->stride_y,
            &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[txb_1d_offset]),
            NOT_USED_VALUE,
            tx_size,
            &context_ptr->three_quad_energy,
            context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
            tx_type,
            PLANE_TYPE_Y,
            pf_shape);
        candidate_buffer->quantized_dc[0][txb_itr] = av1_quantize_inv_quantize(
            pcs_ptr,
            context_ptr,
            &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[txb_1d_offset]),
            NOT_USED_VALUE,
            &(((int32_t *)quant_coeff_ptr->buffer_y)[txb_1d_offset]),
            &(((int32_t *)recon_coeff_ptr->buffer_y)[txb_1d_offset]),
            qindex,
            seg_qp,
            context_ptr->blk_geom->tx_width[tx_depth][txb_itr],
            context_ptr->blk_geom->tx_height[tx_depth][txb_itr] >>
                context_ptr->md_staging_subres_step,
            tx_size,
            &(candidate_buffer->eob[0][txb_itr]),
            &(y_count_non_zero_coeffs_txt),
            COMPONENT_LUMA,
            context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
            tx_type,
            candidate_buffer,
            context_ptr->luma_txb_skip_context,
            context_ptr->luma_dc_sign_context,
            candidate_buffer->candidate_ptr->pred_mode,
            candidate_buffer->candidate_ptr->use_intrabc,
            full_lambda,
            FALSE);
    } else {
        // Init params
        candidate_buffer->quantized_dc[0][txb_itr] = 0;
        candidate_buffer->eob[0][txb_itr]          = 0;
    }

    // Perform T-1 if md_staging_spatial_sse_full_loop_level (assumes TX depth is 0 b/c TXS is assumed off)
    if (context_ptr->md_staging_spatial_sse_full_loop_level) {
        if (y_count_non_zero_coeffs_txt)
            inv_transform_recon_wrapper(candidate_buffer->prediction_ptr->buffer_y,
                                        txb_origin_index,
                                        candidate_buffer->prediction_ptr->stride_y,
                                        recon_ptr->buffer_y,
                                        txb_origin_index,
                                        candidate_buffer->recon_ptr->stride_y,
                                        (int32_t *)recon_coeff_ptr->buffer_y,
                                        txb_1d_offset,
                                        context_ptr->hbd_mode_decision,
                                        context_ptr->blk_geom->txsize[tx_depth][txb_itr],
                                        tx_type,
                                        PLANE_TYPE_Y,
                                        (uint32_t)candidate_buffer->eob[0][txb_itr]);
        else
            svt_av1_picture_copy(candidate_buffer->prediction_ptr,
                                 txb_origin_index,
                                 0,
                                 recon_ptr,
                                 txb_origin_index,
                                 0,
                                 context_ptr->blk_geom->tx_width[tx_depth][txb_itr],
                                 context_ptr->blk_geom->tx_height[tx_depth][txb_itr],
                                 0,
                                 0,
                                 PICTURE_BUFFER_DESC_Y_FLAG,
                                 context_ptr->hbd_mode_decision);

        const int32_t cropped_tx_width = MIN(
            context_ptr->blk_geom->tx_width[tx_depth][txb_itr],
            pcs_ptr->parent_pcs_ptr->aligned_width - (context_ptr->sb_origin_x + tx_org_x));
        const int32_t cropped_tx_height = MIN(
            (uint8_t)(context_ptr->blk_geom->tx_height[tx_depth][txb_itr] >>
                      context_ptr->md_staging_subres_step),
            pcs_ptr->parent_pcs_ptr->aligned_height - (context_ptr->sb_origin_y + tx_org_y));
        EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
            ? svt_full_distortion_kernel16_bits
            : svt_spatial_full_distortion_kernel;
        y_full_distortion[DIST_CALC_PREDICTION]          = spatial_full_dist_type_fun(
            input_picture_ptr->buffer_y,
            input_txb_origin_index,
            input_picture_ptr->stride_y,
            candidate_buffer->prediction_ptr->buffer_y,
            (int32_t)txb_origin_index,
            candidate_buffer->prediction_ptr->stride_y,
            cropped_tx_width,
            cropped_tx_height);
        y_full_distortion[DIST_CALC_RESIDUAL] = spatial_full_dist_type_fun(
            input_picture_ptr->buffer_y,
            input_txb_origin_index,
            input_picture_ptr->stride_y,
            recon_ptr->buffer_y,
            (int32_t)txb_origin_index,
            candidate_buffer->recon_ptr->stride_y,
            cropped_tx_width,
            cropped_tx_height);
        y_full_distortion[DIST_CALC_PREDICTION] <<= 4;
        y_full_distortion[DIST_CALC_RESIDUAL] <<= 4;
    } else {
        // LUMA DISTORTION
        uint32_t txbwidth  = context_ptr->blk_geom->tx_width[tx_depth][txb_itr];
        uint32_t txbheight = (context_ptr->blk_geom->tx_height[tx_depth][txb_itr] >>
                              context_ptr->md_staging_subres_step);

        uint32_t bwidth, bheight;
        if (pf_shape && !tx_search_skip_flag) {
            bwidth  = MAX((txbwidth >> pf_shape), 4);
            bheight = (txbheight >> pf_shape);
        } else {
            bwidth  = txbwidth < 64 ? txbwidth : 32;
            bheight = txbheight < 64 ? txbheight : 32;
        }

        picture_full_distortion32_bits_single(
            &(((int32_t *)context_ptr->tx_coeffs->buffer_y)[txb_1d_offset]),
            &(((int32_t *)recon_coeff_ptr->buffer_y)[txb_1d_offset]),
            txbwidth < 64 ? txbwidth : 32,
            bwidth,
            bheight,
            y_full_distortion,
            y_count_non_zero_coeffs_txt);

        y_full_distortion[DIST_CALC_RESIDUAL] += context_ptr->three_quad_energy;
        y_full_distortion[DIST_CALC_PREDICTION] += context_ptr->three_quad_energy;

        const int32_t shift                   = (MAX_TX_SCALE - av1_get_tx_scale_tab[tx_size]) * 2;
        y_full_distortion[DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(
            y_full_distortion[DIST_CALC_RESIDUAL], shift);
        y_full_distortion[DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(
            y_full_distortion[DIST_CALC_PREDICTION], shift);
    }
    y_full_distortion[DIST_CALC_RESIDUAL] <<= context_ptr->md_staging_subres_step;
    y_full_distortion[DIST_CALC_PREDICTION] <<= context_ptr->md_staging_subres_step;

    //LUMA-ONLY
    const uint32_t th = ((context_ptr->blk_geom->tx_width[tx_depth][txb_itr] *
                          context_ptr->blk_geom->tx_height[tx_depth][txb_itr]) >>
                         6);
    if ((context_ptr->rate_est_ctrls.coeff_rate_est_lvl >= 2 ||
         context_ptr->rate_est_ctrls.coeff_rate_est_lvl == 0) &&
        (y_count_non_zero_coeffs_txt < th))
        *y_coeff_bits = 6000 + y_count_non_zero_coeffs_txt * 1000;
    else if (context_ptr->rate_est_ctrls.coeff_rate_est_lvl == 0)
        *y_coeff_bits = 6000 + y_count_non_zero_coeffs_txt * 400;
    else
        av1_txb_estimate_coeff_bits(context_ptr,
                                    0, //allow_update_cdf,
                                    NULL, //FRAME_CONTEXT *ec_ctx,
                                    pcs_ptr,
                                    candidate_buffer,
                                    txb_1d_offset,
                                    0,
                                    quant_coeff_ptr,
                                    y_count_non_zero_coeffs_txt,
                                    0,
                                    0,
                                    y_coeff_bits,
                                    NULL,
                                    NULL,
                                    tx_size,
                                    context_ptr->blk_geom->txsize_uv[tx_depth][txb_itr],
                                    tx_type,
                                    NOT_USED_VALUE,
                                    COMPONENT_LUMA);

    // Update with best_tx_type data
    candidate_buffer->candidate_ptr->transform_type[txb_itr] = tx_type;
    y_count_non_zero_coeffs[txb_itr]                         = y_count_non_zero_coeffs_txt;
    candidate_buffer->y_has_coeff = ((y_count_non_zero_coeffs_txt > 0) << txb_itr);
    // For Inter blocks, transform type of chroma follows luma transfrom type
    if (is_inter)
        candidate_buffer->candidate_ptr->transform_type_uv =
            candidate_buffer->candidate_ptr->transform_type[txb_itr];
}
void full_loop_core_light_pd0(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                              ModeDecisionCandidateBuffer *candidate_buffer,
                              EbPictureBufferDesc *input_picture_ptr, uint32_t input_origin_index,
                              uint32_t blk_origin_index) {
    uint64_t y_full_distortion[DIST_CALC_TOTAL];
    uint32_t count_non_zero_coeffs;
    uint64_t y_coeff_bits;
    uint32_t full_lambda = context_ptr->hbd_mode_decision
        ? context_ptr->full_sb_lambda_md[EB_10_BIT_MD]
        : context_ptr->full_sb_lambda_md[EB_8_BIT_MD];
    if (context_ptr->subres_ctrls.odd_to_even_deviation_th && context_ptr->pd_pass == PD_PASS_0 &&
        context_ptr->md_stage == MD_STAGE_3 &&
        context_ptr->is_subres_safe == (uint8_t)~0 /* only if invalid*/ &&
        context_ptr->blk_geom->bheight == 64 && context_ptr->blk_geom->bwidth == 64) {
        uint32_t sad_even, sad_odd;
        if (!context_ptr->hbd_mode_decision) {
            sad_even = svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_y + input_origin_index,
                input_picture_ptr->stride_y << 1,
                candidate_buffer->prediction_ptr->buffer_y + blk_origin_index,
                candidate_buffer->prediction_ptr->stride_y << 1,
                context_ptr->blk_geom->bheight >> 1,
                context_ptr->blk_geom->bwidth);

            sad_odd = svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_y + input_origin_index + input_picture_ptr->stride_y,
                input_picture_ptr->stride_y << 1,
                candidate_buffer->prediction_ptr->buffer_y + blk_origin_index +
                    candidate_buffer->prediction_ptr->stride_y,
                candidate_buffer->prediction_ptr->stride_y << 1,
                context_ptr->blk_geom->bheight >> 1,
                context_ptr->blk_geom->bwidth);

        } else {
            sad_even = sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                input_picture_ptr->stride_y << 1,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y) + blk_origin_index,
                candidate_buffer->prediction_ptr->stride_y << 1,
                context_ptr->blk_geom->bheight >> 1,
                context_ptr->blk_geom->bwidth);

            sad_odd = sad_16b_kernel(((uint16_t *)input_picture_ptr->buffer_y) +
                                         input_origin_index + input_picture_ptr->stride_y,
                                     input_picture_ptr->stride_y << 1,
                                     ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y) +
                                         blk_origin_index +
                                         candidate_buffer->prediction_ptr->stride_y,
                                     candidate_buffer->prediction_ptr->stride_y << 1,
                                     context_ptr->blk_geom->bheight >> 1,
                                     context_ptr->blk_geom->bwidth);
        }

        int deviation = (int)(((int)MAX(sad_even, 1) - (int)MAX(sad_odd, 1)) * 100) /
            (int)MAX(sad_odd, 1);
        if (ABS(deviation) <= context_ptr->subres_ctrls.odd_to_even_deviation_th) {
            context_ptr->is_subres_safe = 1;
        } else {
            context_ptr->is_subres_safe = 0;
        }
    }
    if (context_ptr->is_subres_safe != 1)
        context_ptr->md_staging_subres_step = 0;

    // If using 4x subsampling, can't have 8x8 b/c no 8x2 transform
    // If using 2x subsampling, can't have 4x4 b/c no 4x2 transform
    // subres tx assumes NSQ is off
    assert(IMPLIES(context_ptr->md_staging_subres_step == 2, context_ptr->blk_geom->sq_size >= 16));
    assert(IMPLIES(context_ptr->md_staging_subres_step == 1, context_ptr->blk_geom->sq_size >= 8));

    //Y Residual
    residual_kernel(
        input_picture_ptr->buffer_y,
        input_origin_index,
        input_picture_ptr->stride_y << context_ptr->md_staging_subres_step,
        candidate_buffer->prediction_ptr->buffer_y,
        blk_origin_index,
        candidate_buffer->prediction_ptr->stride_y << context_ptr->md_staging_subres_step,
        (int16_t *)candidate_buffer->residual_ptr->buffer_y,
        blk_origin_index,
        candidate_buffer->residual_ptr->stride_y,
        context_ptr->hbd_mode_decision,
        context_ptr->blk_geom->bwidth,
        context_ptr->blk_geom->bheight >> context_ptr->md_staging_subres_step);

    perform_tx_light_pd0(pcs_ptr,
                         context_ptr,
                         candidate_buffer,
                         context_ptr->blk_ptr->qindex,
                         &count_non_zero_coeffs,
                         &y_coeff_bits,
                         &y_full_distortion[0]);
    candidate_buffer->count_non_zero_coeffs = count_non_zero_coeffs;
    av1_full_cost_light_pd0(
        context_ptr, candidate_buffer, y_full_distortion, full_lambda, &y_coeff_bits);
}
/*
  check if we need to do inverse transform and recon
*/
uint8_t do_md_recon(PictureParentControlSet *pcs, ModeDecisionContext *ctxt) {
    uint8_t encdec_bypass = ctxt->bypass_encdec &&
        (ctxt->pd_pass == PD_PASS_1); //if enc dec is bypassed MD has to produce the final recon
    uint8_t need_md_rec_for_intra_pred = !ctxt->skip_intra; //for intra prediction of current frame
    uint8_t need_md_rec_for_ref        = (pcs->is_used_as_reference_flag ||
                                   pcs->scs_ptr->static_config.recon_enabled) &&
        encdec_bypass; //for inter prediction of future frame or if recon is being output
    uint8_t need_md_rec_for_dlf_search  = pcs->dlf_ctrls.enabled; //for DLF levels
    uint8_t need_md_rec_for_cdef_search = pcs->cdef_ctrls.enabled &&
        !pcs->cdef_ctrls.use_reference_cdef_fs; // CDEF search levels needing the recon samples
    uint8_t need_md_rec_for_restoration_search =
        pcs->enable_restoration; // any resoration search level

    uint8_t do_recon;
    if (need_md_rec_for_intra_pred || need_md_rec_for_ref || need_md_rec_for_dlf_search ||
        need_md_rec_for_cdef_search || need_md_rec_for_restoration_search)
        do_recon = 1;
    else
        do_recon = 0;

    return do_recon;
}
static const uint8_t eb_av1_var_offs[MAX_SB_SIZE] = {
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128};
static const uint16_t eb_av1_var_offs_hbd[MAX_SB_SIZE] = {
    512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,
    512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,
    512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,
    512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,
    512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,
    512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,
    512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512};

// Detect blocks that whose chroma component is important (used as a detector for chroma TX shortcuts in reg. PD1 and LPD1)
// Update ctx->chroma_complexity accordingly
void chroma_complexity_check_pred(ModeDecisionContext         *ctx,
                                  ModeDecisionCandidateBuffer *cand_buffer,
                                  EbPictureBufferDesc *input_pic, BlockLocation *loc,
                                  uint8_t use_var) {
    if (ctx->chroma_complexity == COMPONENT_CHROMA)
        return;

    uint32_t y_dist = 0, cb_dist = 0, cr_dist = 0;
    uint8_t  shift = 0;
    shift          = ctx->blk_geom->bheight_uv > 8 ? 2
                 : ctx->blk_geom->bheight_uv > 4   ? 1
                                                   : 0; // no shift for 4x4

    if (!ctx->hbd_mode_decision) {
        y_dist = svt_nxm_sad_kernel_sub_sampled(
            input_pic->buffer_y + loc->input_origin_index,
            input_pic->stride_y << shift,
            cand_buffer->prediction_ptr->buffer_y + loc->blk_origin_index,
            cand_buffer->prediction_ptr->stride_y << shift,
            ctx->blk_geom->bheight_uv >> shift,
            ctx->blk_geom->bwidth_uv);
        // Only need to check Cb component if not already identified as complex
        if (ctx->chroma_complexity == COMPONENT_LUMA ||
            ctx->chroma_complexity == COMPONENT_CHROMA_CR)
            cb_dist = svt_nxm_sad_kernel_sub_sampled(
                input_pic->buffer_cb + loc->input_cb_origin_in_index,
                input_pic->stride_cb << shift,
                cand_buffer->prediction_ptr->buffer_cb + loc->blk_chroma_origin_index,
                cand_buffer->prediction_ptr->stride_cb << shift,
                ctx->blk_geom->bheight_uv >> shift,
                ctx->blk_geom->bwidth_uv);
        // Only need to check Cr component if not already identified as complex
        if (ctx->chroma_complexity == COMPONENT_LUMA ||
            ctx->chroma_complexity == COMPONENT_CHROMA_CB)
            cr_dist = svt_nxm_sad_kernel_sub_sampled(
                input_pic->buffer_cr + loc->input_cb_origin_in_index,
                input_pic->stride_cr << shift,
                cand_buffer->prediction_ptr->buffer_cr + loc->blk_chroma_origin_index,
                cand_buffer->prediction_ptr->stride_cr << shift,
                ctx->blk_geom->bheight_uv >> shift,
                ctx->blk_geom->bwidth_uv);

    } else {
        y_dist = sad_16b_kernel(
            ((uint16_t *)input_pic->buffer_y) + loc->input_origin_index,
            input_pic->stride_y << shift,
            ((uint16_t *)cand_buffer->prediction_ptr->buffer_y) + loc->blk_origin_index,
            cand_buffer->prediction_ptr->stride_y << shift,
            ctx->blk_geom->bheight_uv >> shift,
            ctx->blk_geom->bwidth_uv);
        // Only need to check Cb component if not already identified as complex
        if (ctx->chroma_complexity == COMPONENT_LUMA ||
            ctx->chroma_complexity == COMPONENT_CHROMA_CR)
            cb_dist = sad_16b_kernel(
                ((uint16_t *)input_pic->buffer_cb) + loc->input_cb_origin_in_index,
                input_pic->stride_cb << shift,
                ((uint16_t *)cand_buffer->prediction_ptr->buffer_cb) + loc->blk_chroma_origin_index,
                cand_buffer->prediction_ptr->stride_cb << shift,
                ctx->blk_geom->bheight_uv >> shift,
                ctx->blk_geom->bwidth_uv);
        // Only need to check Cr component if not already identified as complex
        if (ctx->chroma_complexity == COMPONENT_LUMA ||
            ctx->chroma_complexity == COMPONENT_CHROMA_CB)
            cr_dist = sad_16b_kernel(
                ((uint16_t *)input_pic->buffer_cr) + loc->input_cb_origin_in_index,
                input_pic->stride_cr << shift,
                ((uint16_t *)cand_buffer->prediction_ptr->buffer_cr) + loc->blk_chroma_origin_index,
                cand_buffer->prediction_ptr->stride_cr << shift,
                ctx->blk_geom->bheight_uv >> shift,
                ctx->blk_geom->bwidth_uv);
    }
    y_dist <<= 1;

    if (cb_dist > y_dist && cr_dist > y_dist) {
        ctx->chroma_complexity = COMPONENT_CHROMA;
    } else if (cb_dist > y_dist) {
        ctx->chroma_complexity = (ctx->chroma_complexity == COMPONENT_CHROMA_CR)
            ? COMPONENT_CHROMA
            : COMPONENT_CHROMA_CB;
    } else if (cr_dist > y_dist) {
        ctx->chroma_complexity = (ctx->chroma_complexity == COMPONENT_CHROMA_CB)
            ? COMPONENT_CHROMA
            : COMPONENT_CHROMA_CR;
    }

    if (use_var) {
        const AomVarianceFnPtr *fn_ptr = &mefn_ptr[ctx->blk_geom->bsize_uv];
        unsigned int            sse;
        unsigned int            var_cb;
        unsigned int            var_cr;
        if (ctx->hbd_mode_decision) {
            var_cb = fn_ptr->vf_hbd_10(CONVERT_TO_BYTEPTR(((uint16_t *)input_pic->buffer_cb) +
                                                          loc->input_cb_origin_in_index),
                                       input_pic->stride_cb,
                                       CONVERT_TO_BYTEPTR(eb_av1_var_offs_hbd),
                                       0,
                                       &sse);
            var_cr = fn_ptr->vf_hbd_10(CONVERT_TO_BYTEPTR(((uint16_t *)input_pic->buffer_cr) +
                                                          loc->input_cb_origin_in_index),
                                       input_pic->stride_cr,
                                       CONVERT_TO_BYTEPTR(eb_av1_var_offs_hbd),
                                       0,
                                       &sse);
        } else {
            var_cb = fn_ptr->vf(input_pic->buffer_cb + loc->input_cb_origin_in_index,
                                input_pic->stride_cb,
                                eb_av1_var_offs,
                                0,
                                &sse);
            var_cr = fn_ptr->vf(input_pic->buffer_cr + loc->input_cb_origin_in_index,
                                input_pic->stride_cr,
                                eb_av1_var_offs,
                                0,
                                &sse);
        }

        int block_var_cb = ROUND_POWER_OF_TWO(var_cb,
                                              num_pels_log2_lookup[ctx->blk_geom->bsize_uv]);
        int block_var_cr = ROUND_POWER_OF_TWO(var_cr,
                                              num_pels_log2_lookup[ctx->blk_geom->bsize_uv]);

        // th controls how safe the detector is (can be changed in the future, or made a parameter)
        uint16_t th = 150;
        if (block_var_cb > th && block_var_cr > th) {
            ctx->chroma_complexity = COMPONENT_CHROMA;
        } else if (block_var_cb > th) {
            ctx->chroma_complexity = (ctx->chroma_complexity == COMPONENT_CHROMA_CR)
                ? COMPONENT_CHROMA
                : COMPONENT_CHROMA_CB;
        } else if (block_var_cr > th) {
            ctx->chroma_complexity = (ctx->chroma_complexity == COMPONENT_CHROMA_CB)
                ? COMPONENT_CHROMA
                : COMPONENT_CHROMA_CR;
        }
    }
}

// Detect blocks that whose chroma component is important (used as a detector for skipping the chroma TX path in LPD1)
COMPONENT_TYPE chroma_complexity_check(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                       ModeDecisionCandidate *candidate_ptr,
                                       EbPictureBufferDesc *input_picture_ptr, BlockLocation *loc) {
    /* For INTER blocks, compute the luma/chroma full-pel distortions; if chroma distortion is much higher, then block is complex
    in chroma, and chroma should be performed. */
    if (is_inter_mode(candidate_ptr->pred_mode)) {
        EbPictureBufferDesc *ref_pic;
        EbReferenceObject   *ref_obj;
        int16_t              mv_x, mv_y;
        MvReferenceFrame     rf[2];
        av1_set_ref_frame(rf, candidate_ptr->ref_frame_type);
        const int8_t ref_idx_first = get_ref_frame_idx(rf[0]);
        if (rf[1] != NONE_FRAME || get_list_idx(rf[0]) == UNI_PRED_LIST_0) {
            ref_obj = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[0][ref_idx_first]->object_ptr;
            ref_pic = get_ref_pic_buffer(pcs_ptr, context_ptr->hbd_mode_decision, 0, ref_idx_first);
            mv_x    = candidate_ptr->mv[REF_LIST_0].x >> 3;
            mv_y    = candidate_ptr->mv[REF_LIST_0].y >> 3;
        } else {
            ref_obj = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[1][ref_idx_first]->object_ptr;
            ref_pic = get_ref_pic_buffer(pcs_ptr, context_ptr->hbd_mode_decision, 1, ref_idx_first);
            mv_x    = candidate_ptr->mv[REF_LIST_1].x >> 3;
            mv_y    = candidate_ptr->mv[REF_LIST_1].y >> 3;
        }
        // -------
        // Use scaled references if resolution of the reference is different from that of the input
        // -------
        use_scaled_rec_refs_if_needed(
            pcs_ptr, input_picture_ptr, ref_obj, &ref_pic, context_ptr->hbd_mode_decision);

        uint32_t src_y_offset = ref_pic->origin_x + context_ptr->blk_origin_x + mv_x +
            (ref_pic->origin_y + context_ptr->blk_origin_y + mv_y) * ref_pic->stride_y;
        uint32_t src_cb_offset = ((ref_pic->origin_x + context_ptr->blk_origin_x + mv_x) >> 1) +
            (((ref_pic->origin_y + context_ptr->blk_origin_y + mv_y) >> 1)) * ref_pic->stride_cb;
        uint32_t src_cr_offset = ((ref_pic->origin_x + context_ptr->blk_origin_x + mv_x) >> 1) +
            (((ref_pic->origin_y + context_ptr->blk_origin_y + mv_y) >> 1)) * ref_pic->stride_cr;
        uint8_t shift = 0;
        if (context_ptr->lpd1_tx_ctrls.chroma_detector_level >= 2)
            shift = context_ptr->blk_geom->bheight_uv > 8 ? 2
                : context_ptr->blk_geom->bheight_uv > 4   ? 1
                                                          : 0; // no shift for 4x4
        else
            shift = context_ptr->blk_geom->bheight_uv > 4 ? 1 : 0; // no shift for 4x4

        uint32_t y_dist, cb_dist, cr_dist;

        if (context_ptr->hbd_mode_decision) {
            uint16_t *src_10b;
            DECLARE_ALIGNED(16, uint16_t, packed_buf[PACKED_BUFFER_SIZE]);
            // pack the reference into temp 16bit buffer
            uint8_t offset = 0;
            int32_t stride;

            pack_block(
                ref_pic->buffer_y + src_y_offset - offset - (offset * (ref_pic->stride_y << shift)),
                ref_pic->stride_y << shift,
                ref_pic->buffer_bit_inc_y + src_y_offset - offset -
                    (offset * (ref_pic->stride_bit_inc_y << shift)),
                ref_pic->stride_bit_inc_y << shift,
                (uint16_t *)packed_buf,
                MAX_SB_SIZE,
                context_ptr->blk_geom->bwidth_uv + (offset << 1),
                (context_ptr->blk_geom->bheight_uv >> shift) + (offset << 1));

            src_10b = (uint16_t *)packed_buf + offset + (offset * MAX_SB_SIZE);
            stride  = MAX_SB_SIZE;

            // Y dist only computed over UV size so SADs are comparable
            y_dist = sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_y) + loc->input_origin_index,
                input_picture_ptr->stride_y << shift,
                src_10b,
                stride,
                context_ptr->blk_geom->bheight_uv >> shift,
                context_ptr->blk_geom->bwidth_uv);

            // pack the reference into temp 16bit buffer

            pack_block(ref_pic->buffer_cb + src_cb_offset - offset -
                           (offset * (ref_pic->stride_cb << shift)),
                       ref_pic->stride_cb << shift,
                       ref_pic->buffer_bit_inc_cb + src_cb_offset - offset -
                           (offset * (ref_pic->stride_bit_inc_cb << shift)),
                       ref_pic->stride_bit_inc_cb << shift,
                       (uint16_t *)packed_buf,
                       MAX_SB_SIZE,
                       context_ptr->blk_geom->bwidth_uv + (offset << 1),
                       (context_ptr->blk_geom->bheight_uv >> shift) + (offset << 1));

            src_10b = (uint16_t *)packed_buf + offset + (offset * MAX_SB_SIZE);
            stride  = MAX_SB_SIZE;

            cb_dist = sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_cb) + loc->input_cb_origin_in_index,
                input_picture_ptr->stride_cb << shift,
                src_10b,
                stride,
                context_ptr->blk_geom->bheight_uv >> shift,
                context_ptr->blk_geom->bwidth_uv);

            // pack the reference into temp 16bit buffer
            pack_block(ref_pic->buffer_cr + src_cr_offset - offset -
                           (offset * (ref_pic->stride_cr << shift)),
                       ref_pic->stride_cr << shift,
                       ref_pic->buffer_bit_inc_cr + src_cr_offset - offset -
                           (offset * (ref_pic->stride_bit_inc_cr << shift)),
                       ref_pic->stride_bit_inc_cr << shift,
                       (uint16_t *)packed_buf,
                       MAX_SB_SIZE,
                       context_ptr->blk_geom->bwidth_uv + (offset << 1),
                       (context_ptr->blk_geom->bheight_uv >> shift) + (offset << 1));

            src_10b = (uint16_t *)packed_buf + offset + (offset * MAX_SB_SIZE);
            stride  = MAX_SB_SIZE;

            cr_dist = sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_cr) + loc->input_cb_origin_in_index,
                input_picture_ptr->stride_cr << shift,
                src_10b,
                stride,
                context_ptr->blk_geom->bheight_uv >> shift,
                context_ptr->blk_geom->bwidth_uv);
        } else {
            // Y dist only computed over UV size so SADs are comparable
            y_dist = svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_y + loc->input_origin_index,
                input_picture_ptr->stride_y << shift,
                ref_pic->buffer_y + src_y_offset,
                ref_pic->stride_y << shift,
                context_ptr->blk_geom->bheight_uv >> shift,
                context_ptr->blk_geom->bwidth_uv);

            cb_dist = svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_cb + loc->input_cb_origin_in_index,
                input_picture_ptr->stride_cb << shift,
                ref_pic->buffer_cb + src_cb_offset,
                ref_pic->stride_cb << shift,
                context_ptr->blk_geom->bheight_uv >> shift,
                context_ptr->blk_geom->bwidth_uv);

            cr_dist = svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_cr + loc->input_cb_origin_in_index,
                input_picture_ptr->stride_cr << shift,
                ref_pic->buffer_cr + src_cr_offset,
                ref_pic->stride_cr << shift,
                context_ptr->blk_geom->bheight_uv >> shift,
                context_ptr->blk_geom->bwidth_uv);
        }
        // shift y_dist by to ensure chroma is much higher than luma
        if (context_ptr->lpd1_tx_ctrls.chroma_detector_level >= 2)
            y_dist <<= 2;
        else
            y_dist <<= 1;

        if (cb_dist > y_dist && cr_dist > y_dist) {
            return COMPONENT_CHROMA;
        } else if (cb_dist > y_dist) {
            return COMPONENT_CHROMA_CB;
        } else if (cr_dist > y_dist) {
            return COMPONENT_CHROMA_CR;
        }
    }

    /* For INTRA blocks, if the chroma variance of the block is high, perform chroma. Can also use variance check as an additional
    check for INTER blocks. */
    if (is_intra_mode(candidate_ptr->pred_mode) ||
        context_ptr->lpd1_tx_ctrls.chroma_detector_level <= 2) {
        const AomVarianceFnPtr *fn_ptr = &mefn_ptr[context_ptr->blk_geom->bsize_uv];
        unsigned int            sse;
        unsigned int            var_cb;
        unsigned int            var_cr;
        if (context_ptr->hbd_mode_decision) {
            var_cb = fn_ptr->vf_hbd_10(
                CONVERT_TO_BYTEPTR(((uint16_t *)input_picture_ptr->buffer_cb) +
                                   loc->input_cb_origin_in_index),
                input_picture_ptr->stride_cb,
                CONVERT_TO_BYTEPTR(eb_av1_var_offs_hbd),
                0,
                &sse);
            var_cr = fn_ptr->vf_hbd_10(
                CONVERT_TO_BYTEPTR(((uint16_t *)input_picture_ptr->buffer_cr) +
                                   loc->input_cb_origin_in_index),
                input_picture_ptr->stride_cr,
                CONVERT_TO_BYTEPTR(eb_av1_var_offs_hbd),
                0,
                &sse);
        } else {
            var_cb = fn_ptr->vf(input_picture_ptr->buffer_cb + loc->input_cb_origin_in_index,
                                input_picture_ptr->stride_cb,
                                eb_av1_var_offs,
                                0,
                                &sse);
            var_cr = fn_ptr->vf(input_picture_ptr->buffer_cr + loc->input_cb_origin_in_index,
                                input_picture_ptr->stride_cr,
                                eb_av1_var_offs,
                                0,
                                &sse);
        }
        int block_var_cb = ROUND_POWER_OF_TWO(
            var_cb, num_pels_log2_lookup[context_ptr->blk_geom->bsize_uv]);
        int block_var_cr = ROUND_POWER_OF_TWO(
            var_cr, num_pels_log2_lookup[context_ptr->blk_geom->bsize_uv]);

        // th controls how safe the detector is
        uint16_t th = context_ptr->lpd1_tx_ctrls.chroma_detector_level <= 1 ? 75 : 150;
        if (block_var_cb > th && block_var_cr > th) {
            return COMPONENT_CHROMA;
        } else if (block_var_cb > th) {
            return COMPONENT_CHROMA_CB;
        } else if (block_var_cr > th) {
            return COMPONENT_CHROMA_CR;
        }
    }

    // At end, complex chroma was not detected, so only chroma path can be skipped
    return COMPONENT_LUMA;
}
/*
   full loop core for light PD1 path
*/
void full_loop_core_light_pd1(PictureControlSet *pcs_ptr, BlkStruct *blk_ptr,
                              ModeDecisionContext         *context_ptr,
                              ModeDecisionCandidateBuffer *candidate_buffer,
                              ModeDecisionCandidate       *candidate_ptr,
                              EbPictureBufferDesc *input_picture_ptr, BlockLocation *loc) {
    uint64_t y_full_distortion[DIST_CALC_TOTAL];
    uint64_t cb_full_distortion[DIST_CALC_TOTAL];
    uint64_t cr_full_distortion[DIST_CALC_TOTAL];
    uint64_t y_coeff_bits;
    uint64_t cb_coeff_bits;
    uint64_t cr_coeff_bits;
    candidate_ptr->skip_mode = FALSE;
    Bool perform_tx          = 1;
    if (context_ptr->lpd1_allow_skipping_tx) {
        if (context_ptr->lpd1_skip_inter_tx_level == 2 && is_inter_mode(candidate_ptr->pred_mode))
            perform_tx = 0;

        MacroBlockD *xd = context_ptr->blk_ptr->av1xd;
        if (xd->left_available && xd->up_available) {
            const BlockModeInfoEnc *const left_mi  = &xd->left_mbmi->block_mi;
            const BlockModeInfoEnc *const above_mi = &xd->above_mbmi->block_mi;
            if (left_mi->skip && above_mi->skip &&
                ((left_mi->mode == NEAREST_NEARESTMV && above_mi->mode == NEAREST_NEARESTMV) ||
                 context_ptr->lpd1_skip_inter_tx_level)) {
                /* For M12 and below, do not skip TX for candidates other than NRST_NRST and do not remove the check on neighbouring
                    coeffs, as that may introduce blocking artifacts in certain clips. */

                // Skip TX for NRST_NRST
                if (context_ptr->lpd1_tx_ctrls.skip_nrst_nrst_luma_tx &&
                    candidate_ptr->pred_mode == NEAREST_NEARESTMV)
                    perform_tx = 0;
                // Skip TX for INTER - should only be true for M13
                else if (context_ptr->lpd1_skip_inter_tx_level == 1 &&
                         is_inter_mode(candidate_ptr->pred_mode))
                    perform_tx = 0;
            }
        }
    }
    const uint8_t recon_needed = do_md_recon(pcs_ptr->parent_pcs_ptr, context_ptr);

    // If need 10bit prediction, perform luma compensation before TX
    if ((perform_tx || recon_needed) && context_ptr->hbd_mode_decision) {
        context_ptr->md_stage                    = MD_STAGE_0;
        context_ptr->end_plane                   = 1;
        context_ptr->md_staging_skip_chroma_pred = 1;
        context_ptr->uv_intra_comp_only          = FALSE;
        svt_product_prediction_fun_table_light_pd1[is_inter_mode(candidate_ptr->pred_mode)](
            context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
        context_ptr->md_stage                    = MD_STAGE_3;
        context_ptr->end_plane                   = MAX_MB_PLANE;
        context_ptr->uv_intra_comp_only          = TRUE;
        context_ptr->md_staging_skip_chroma_pred = 0;
    }
    if (perform_tx) {
        perform_dct_dct_tx_light_pd1(
            pcs_ptr, context_ptr, candidate_buffer, loc, &y_coeff_bits, &y_full_distortion[0]);
    } else {
        candidate_buffer->eob[0][0]                        = 0;
        candidate_buffer->quantized_dc[0][0]               = 0;
        candidate_buffer->y_has_coeff                      = 0;
        y_full_distortion[DIST_CALC_RESIDUAL]              = 0;
        y_full_distortion[DIST_CALC_PREDICTION]            = 0;
        y_coeff_bits                                       = 6000;
        candidate_buffer->candidate_ptr->transform_type[0] = DCT_DCT;
        // For Inter blocks, transform type of chroma follows luma transfrom type
        if (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode))
            candidate_buffer->candidate_ptr->transform_type_uv = DCT_DCT;
    }
    // Update coeff info based on luma TX so that chroma can take advantage of most accurate info
    candidate_buffer->block_has_coeff       = (candidate_buffer->y_has_coeff) ? 1 : 0;
    candidate_buffer->count_non_zero_coeffs = candidate_buffer->eob[0][0];
    uint8_t perform_chroma                  = candidate_buffer->block_has_coeff ||
        !(context_ptr->lpd1_tx_ctrls.zero_y_coeff_exit);
    COMPONENT_TYPE chroma_component = COMPONENT_CHROMA;
    context_ptr->chroma_complexity  = COMPONENT_LUMA;

    // If going to skip chroma TX, detect if block is complex in chroma, and if so, force chroma to be performed.
    if (!perform_chroma) {
        if (context_ptr->lpd1_tx_ctrls.chroma_detector_level) {
            chroma_component = chroma_complexity_check(
                pcs_ptr, context_ptr, candidate_ptr, input_picture_ptr, loc);

            if (context_ptr->lpd1_tx_ctrls.chroma_detector_level <= 3)
                context_ptr->chroma_complexity = chroma_component;
        } else {
            chroma_component = COMPONENT_LUMA;
        }

        perform_chroma = chroma_component > COMPONENT_LUMA;

        if (chroma_component == COMPONENT_CHROMA_CB) {
            cr_full_distortion[DIST_CALC_RESIDUAL]   = 0;
            cr_full_distortion[DIST_CALC_PREDICTION] = 0;
            cr_coeff_bits                            = 0;
            candidate_buffer->v_has_coeff            = 0;
            candidate_buffer->eob[2][0]              = 0;
        } else if (chroma_component == COMPONENT_CHROMA_CR) {
            cb_full_distortion[DIST_CALC_RESIDUAL]   = 0;
            cb_full_distortion[DIST_CALC_PREDICTION] = 0;
            cb_coeff_bits                            = 0;
            candidate_buffer->u_has_coeff            = 0;
            candidate_buffer->eob[1][0]              = 0;
        }
    }
    context_ptr->lpd1_chroma_comp = recon_needed ? COMPONENT_CHROMA : chroma_component;
    // If no luma coeffs, may skip chroma TX and full cost calc (skip chroma compensation
    // if not needed for recon)
    if (perform_chroma) {
        // If using chroma pred samples in the next chroma complexity detector, need to generate pred samples for all components
        if (!recon_needed)
            context_ptr->lpd1_chroma_comp = context_ptr->lpd1_tx_ctrls.chroma_detector_level <= 3
                ? COMPONENT_CHROMA
                : chroma_component;
        //Chroma Prediction
        svt_product_prediction_fun_table_light_pd1[is_inter_mode(candidate_ptr->pred_mode)](
            context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
        // Perform additional check to detect complex chroma blocks
        if (context_ptr->lpd1_tx_ctrls.chroma_detector_level &&
            context_ptr->lpd1_tx_ctrls.chroma_detector_level <= 3 &&
            context_ptr->chroma_complexity != COMPONENT_CHROMA &&
            (context_ptr->use_tx_shortcuts_mds3 ||
             context_ptr->lpd1_tx_ctrls.use_uv_shortcuts_on_y_coeffs)) {
            chroma_complexity_check_pred(
                context_ptr, candidate_buffer, input_picture_ptr, loc, 0 /*use_var*/);
        }

        //CHROMA
        full_loop_chroma_light_pd1(pcs_ptr,
                                   context_ptr,
                                   candidate_buffer,
                                   input_picture_ptr,
                                   loc->input_cb_origin_in_index,
                                   loc->blk_chroma_origin_index,
                                   chroma_component,
                                   context_ptr->qp_index,
                                   cb_full_distortion,
                                   cr_full_distortion,
                                   &cb_coeff_bits,
                                   &cr_coeff_bits);
        candidate_buffer->block_has_coeff = (candidate_buffer->y_has_coeff ||
                                             candidate_buffer->u_has_coeff ||
                                             candidate_buffer->v_has_coeff)
            ? TRUE
            : FALSE;
        //ALL PLANE
        svt_av1_product_full_cost_func_table[is_inter_mode(candidate_ptr->pred_mode)](
            pcs_ptr,
            context_ptr,
            candidate_buffer,
            blk_ptr,
            y_full_distortion,
            cb_full_distortion,
            cr_full_distortion,
            context_ptr->hbd_mode_decision ? context_ptr->full_lambda_md[EB_10_BIT_MD]
                                           : context_ptr->full_lambda_md[EB_8_BIT_MD],
            &y_coeff_bits,
            &cb_coeff_bits,
            &cr_coeff_bits,
            context_ptr->blk_geom->bsize);
    } else {
        // Only need chroma pred if generating recon
        if (context_ptr->lpd1_chroma_comp > COMPONENT_LUMA) {
            //Chroma Prediction
            svt_product_prediction_fun_table_light_pd1[is_inter_mode(candidate_ptr->pred_mode)](
                context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
        }
        candidate_buffer->u_has_coeff = candidate_buffer->v_has_coeff = 0;
        if (candidate_ptr->skip_mode_allowed)
            candidate_ptr->skip_mode = TRUE;
    }
}

// Derive the start and end TX depths based on block characteristics
static INLINE void get_start_end_tx_depth(PictureParentControlSet *ppcs, ModeDecisionContext *ctx,
                                          ModeDecisionCandidate *candidate_ptr,
                                          uint8_t *start_tx_depth, uint8_t *end_tx_depth) {
    TxsControls           *txs_ctrls = &ctx->txs_ctrls;
    const BlockGeom *const blk_geom  = ctx->blk_geom;

    if (txs_ctrls->enabled == 0) {
        *start_tx_depth = *end_tx_depth = 0;
    } else if (ctx->md_staging_tx_size_mode == 0) {
        *start_tx_depth = *end_tx_depth = candidate_ptr->tx_depth;
    } else {
        *start_tx_depth = 0;
        // end_tx_depth set to zero for blocks which go beyond the picture boundaries
        if ((ctx->sb_origin_x + blk_geom->origin_x + ctx->blk_geom->bwidth <= ppcs->aligned_width &&
             ctx->sb_origin_y + blk_geom->origin_y + ctx->blk_geom->bheight <=
                 ppcs->aligned_height))
            *end_tx_depth = get_end_tx_depth(blk_geom->bsize);
        else
            *end_tx_depth = 0;
    }

    *end_tx_depth = MIN(*end_tx_depth,
                        (is_intra_mode(candidate_ptr->pred_mode)
                             ? (blk_geom->shape == PART_N ? txs_ctrls->intra_class_max_depth_sq
                                                          : txs_ctrls->intra_class_max_depth_nsq)
                             : (blk_geom->shape == PART_N ? txs_ctrls->inter_class_max_depth_sq
                                                          : txs_ctrls->inter_class_max_depth_nsq)));
}

void full_loop_core(PictureControlSet *pcs_ptr, BlkStruct *blk_ptr,
                    ModeDecisionContext *context_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
                    ModeDecisionCandidate *candidate_ptr, EbPictureBufferDesc *input_picture_ptr,
                    uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                    uint32_t blk_origin_index, uint32_t blk_chroma_origin_index) {
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
    int32_t  is_inter      = (is_inter_mode(candidate_buffer->candidate_ptr->pred_mode) ||
                        candidate_buffer->candidate_ptr->use_intrabc)
              ? TRUE
              : FALSE;
    // initialize TU Split
    y_full_distortion[DIST_CALC_RESIDUAL]   = 0;
    y_full_distortion[DIST_CALC_PREDICTION] = 0;
    y_coeff_bits                            = 0;
    candidate_buffer->full_distortion       = 0;
    // Set Skip Flag
    candidate_ptr->skip_mode = FALSE;
    if (is_inter_mode(candidate_ptr->pred_mode)) {
        if (context_ptr->md_staging_perform_inter_pred) {
            // Perform INTER prediction
            svt_product_prediction_fun_table[1](
                context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
        }
    } else if (context_ptr->md_staging_skip_full_chroma == FALSE) {
        if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            // Cb/Cr Prediction
            if (context_ptr->md_staging_perform_intra_chroma_pred) {
                context_ptr->uv_intra_comp_only = context_ptr->need_hbd_comp_mds3 ? FALSE : TRUE;
                // Here, the mode is INTRA, but if intra_bc is used, must use inter prediction function
                svt_product_prediction_fun_table[candidate_buffer->candidate_ptr->use_intrabc](
                    context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
            }
        }
    }
    // Initialize luma CBF
    candidate_buffer->y_has_coeff = 0;
    candidate_buffer->u_has_coeff = 0;
    candidate_buffer->v_has_coeff = 0;
    uint8_t start_tx_depth        = 0;
    uint8_t end_tx_depth          = 0;
    get_start_end_tx_depth(
        pcs_ptr->parent_pcs_ptr, context_ptr, candidate_ptr, &start_tx_depth, &end_tx_depth);
    if (context_ptr->subres_ctrls.odd_to_even_deviation_th && context_ptr->pd_pass == PD_PASS_0 &&
        context_ptr->md_stage == MD_STAGE_3 &&
        context_ptr->is_subres_safe == (uint8_t)~0 /* only if invalid*/ &&
        context_ptr->blk_geom->bheight == 64 && context_ptr->blk_geom->bwidth == 64) {
        uint32_t sad_even, sad_odd;
        if (!context_ptr->hbd_mode_decision) {
            sad_even = svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_y + input_origin_index,
                input_picture_ptr->stride_y << 1,
                candidate_buffer->prediction_ptr->buffer_y + blk_origin_index,
                candidate_buffer->prediction_ptr->stride_y << 1,
                context_ptr->blk_geom->bheight >> 1,
                context_ptr->blk_geom->bwidth);

            sad_odd = svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_y + input_origin_index + input_picture_ptr->stride_y,
                input_picture_ptr->stride_y << 1,
                candidate_buffer->prediction_ptr->buffer_y + blk_origin_index +
                    candidate_buffer->prediction_ptr->stride_y,
                candidate_buffer->prediction_ptr->stride_y << 1,
                context_ptr->blk_geom->bheight >> 1,
                context_ptr->blk_geom->bwidth);

        } else {
            sad_even = sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_y) + input_origin_index,
                input_picture_ptr->stride_y << 1,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y) + blk_origin_index,
                candidate_buffer->prediction_ptr->stride_y << 1,
                context_ptr->blk_geom->bheight >> 1,
                context_ptr->blk_geom->bwidth);

            sad_odd = sad_16b_kernel(((uint16_t *)input_picture_ptr->buffer_y) +
                                         input_origin_index + input_picture_ptr->stride_y,
                                     input_picture_ptr->stride_y << 1,
                                     ((uint16_t *)candidate_buffer->prediction_ptr->buffer_y) +
                                         blk_origin_index +
                                         candidate_buffer->prediction_ptr->stride_y,
                                     candidate_buffer->prediction_ptr->stride_y << 1,
                                     context_ptr->blk_geom->bheight >> 1,
                                     context_ptr->blk_geom->bwidth);
        }

        int deviation = (int)(((int)MAX(sad_even, 1) - (int)MAX(sad_odd, 1)) * 100) /
            (int)MAX(sad_odd, 1);
        if (ABS(deviation) <= context_ptr->subres_ctrls.odd_to_even_deviation_th) {
            context_ptr->is_subres_safe = 1;
        } else {
            context_ptr->is_subres_safe = 0;
        }
    }

    context_ptr->md_staging_subres_step = (context_ptr->is_subres_safe == 1)
        ? context_ptr->md_staging_subres_step
        : 0;
    //Y Residual: residual for INTRA is computed inside the TU loop
    if (is_inter)
        //Y Residual
        residual_kernel(
            input_picture_ptr->buffer_y,
            input_origin_index,
            input_picture_ptr->stride_y << context_ptr->md_staging_subres_step,
            candidate_buffer->prediction_ptr->buffer_y,
            blk_origin_index,
            candidate_buffer->prediction_ptr->stride_y << context_ptr->md_staging_subres_step,
            (int16_t *)candidate_buffer->residual_ptr->buffer_y,
            blk_origin_index,
            candidate_buffer->residual_ptr->stride_y,
            context_ptr->hbd_mode_decision,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight >> context_ptr->md_staging_subres_step);
    if (context_ptr->perform_mds1 && context_ptr->md_stage == MD_STAGE_3 &&
        context_ptr->tx_shortcut_ctrls.bypass_tx_when_zcoeff &&
        !candidate_buffer->block_has_coeff) {
        start_tx_depth = 0;
        end_tx_depth   = 0;
    }
    // Check if should perform TX type search
    if (pcs_ptr->scs_ptr->super_block_size == 64 && start_tx_depth == 0 &&
        end_tx_depth == 0 && // TXS off
        !pcs_ptr->parent_pcs_ptr
             ->sc_class1 && // Can't be SC b/c SC tries DCT_DCT and IDTX when only_dct_dct is 1
        search_dct_dct_only(pcs_ptr,
                            context_ptr,
                            candidate_buffer,
                            0 /*tx_depth*/,
                            0 /*txb_itr*/,
                            is_inter)) { // TXT off

        uint8_t tx_search_skip_flag = 0;
        if (context_ptr->perform_mds1 && context_ptr->md_stage == MD_STAGE_3 &&
            context_ptr->tx_shortcut_ctrls.bypass_tx_when_zcoeff &&
            !candidate_buffer->block_has_coeff)
            tx_search_skip_flag = 1;
        perform_dct_dct_tx(pcs_ptr,
                           context_ptr,
                           candidate_buffer,
                           tx_search_skip_flag,
                           context_ptr->blk_ptr->qindex,
                           &(*count_non_zero_coeffs[0]),
                           &y_coeff_bits,
                           &y_full_distortion[0]);
    } else
        perform_tx_partitioning(candidate_buffer,
                                context_ptr,
                                pcs_ptr,
                                start_tx_depth,
                                end_tx_depth,
                                context_ptr->blk_ptr->qindex,
                                &(*count_non_zero_coeffs[0]),
                                &y_coeff_bits,
                                &y_full_distortion[0]);
    // Update coeff info based on luma TX so that chroma can take advantage of most accurate info
    candidate_buffer->block_has_coeff = (candidate_buffer->y_has_coeff) ? 1 : 0;

    uint16_t txb_count =
        context_ptr->blk_geom->txb_count[candidate_buffer->candidate_ptr->tx_depth];
    candidate_buffer->count_non_zero_coeffs = 0;
    for (uint8_t txb_itr = 0; txb_itr < txb_count; txb_itr++)
        candidate_buffer->count_non_zero_coeffs += count_non_zero_coeffs[0][txb_itr];
    //CHROMA

    cb_full_distortion[DIST_CALC_RESIDUAL]   = 0;
    cr_full_distortion[DIST_CALC_RESIDUAL]   = 0;
    cb_full_distortion[DIST_CALC_PREDICTION] = 0;
    cr_full_distortion[DIST_CALC_PREDICTION] = 0;

    cb_coeff_bits = 0;
    cr_coeff_bits = 0;

    context_ptr->chroma_complexity = COMPONENT_LUMA;
    if (context_ptr->tx_shortcut_ctrls.chroma_detector_level &&
        context_ptr->md_stage == MD_STAGE_3 &&
        (context_ptr->tx_shortcut_ctrls.apply_pf_on_coeffs || context_ptr->use_tx_shortcuts_mds3)) {
        BlockLocation loc;
        loc.input_origin_index       = input_origin_index;
        loc.input_cb_origin_in_index = input_cb_origin_in_index;
        loc.blk_origin_index         = blk_origin_index;
        loc.blk_chroma_origin_index  = blk_chroma_origin_index;

        chroma_complexity_check_pred(
            context_ptr, candidate_buffer, input_picture_ptr, &loc, 1 /*use_var*/);
    }

    // FullLoop and TU search
    uint16_t cb_qindex = context_ptr->qp_index;
    if (context_ptr->md_staging_skip_full_chroma == FALSE) {
        if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
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
        Bool cfl_performed = FALSE;
        if (!is_inter)
            if (candidate_buffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {
                cfl_performed = TRUE;
                // If mode is CFL:
                // 1: recon the Luma
                // 2: Form the pred_buf_q3
                // 3: Loop over alphas and find the best or choose DC
                // 4: Recalculate the residual for chroma
                cfl_prediction(pcs_ptr,
                               candidate_buffer,
                               context_ptr,
                               input_picture_ptr,
                               input_cb_origin_in_index,
                               blk_chroma_origin_index);
            }
        if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            svt_aom_full_loop_uv(pcs_ptr,
                                 context_ptr,
                                 candidate_buffer,
                                 input_picture_ptr,
                                 COMPONENT_CHROMA,
                                 cb_qindex,
                                 count_non_zero_coeffs,
                                 cb_full_distortion,
                                 cr_full_distortion,
                                 &cb_coeff_bits,
                                 &cr_coeff_bits,
                                 1);
        }

        // Check independant chroma vs. cfl
        if (!is_inter)
            if (candidate_ptr->palette_info == NULL || candidate_ptr->palette_size[0] == 0)
                if (context_ptr->blk_geom->has_uv && context_ptr->uv_ctrls.uv_mode == CHROMA_MODE_0)
                    if (cfl_performed)
                        check_best_indepedant_cfl(pcs_ptr,
                                                  input_picture_ptr,
                                                  context_ptr,
                                                  input_cb_origin_in_index,
                                                  blk_chroma_origin_index,
                                                  candidate_buffer,
                                                  (uint8_t)cb_qindex,
                                                  cb_full_distortion,
                                                  cr_full_distortion,
                                                  &cb_coeff_bits,
                                                  &cr_coeff_bits);
    }
    candidate_buffer->block_has_coeff = (candidate_buffer->y_has_coeff ||
                                         candidate_buffer->u_has_coeff ||
                                         candidate_buffer->v_has_coeff)
        ? TRUE
        : FALSE;
    //ALL PLANE
    svt_av1_product_full_cost_func_table[is_inter_mode(candidate_ptr->pred_mode)](
        pcs_ptr,
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
static void md_stage_1(PictureControlSet *pcs_ptr, BlkStruct *blk_ptr,
                       ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                       uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                       uint32_t blk_origin_index, uint32_t blk_chroma_origin_index) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array = &(
        candidate_buffer_ptr_array_base[0]);

    // Set MD Staging full_loop_core settings
    context_ptr->md_staging_tx_size_mode     = 0;
    context_ptr->md_staging_txt_level        = 0;
    context_ptr->md_staging_skip_full_chroma = TRUE;

    context_ptr->md_staging_spatial_sse_full_loop_level = FALSE;
    context_ptr->md_staging_fast_coeff_est_level        = (context_ptr->pd_pass == PD_PASS_1)
               ? 1
               : context_ptr->rate_est_ctrls.pd0_fast_coeff_est_level;
    context_ptr->md_staging_subres_step                 = context_ptr->subres_ctrls.step;
    context_ptr->md_staging_skip_chroma_pred            = TRUE;
    context_ptr->md_staging_perform_intra_chroma_pred   = FALSE;
    context_ptr->md_staging_perform_inter_pred = (context_ptr->ifs_ctrls.level == IFS_MDS1) ? TRUE
                                                                                            : FALSE;
    context_ptr->md_staging_skip_interpolation_search = (context_ptr->ifs_ctrls.level == IFS_MDS1)
        ? FALSE
        : TRUE;
    context_ptr->end_plane                            = (context_ptr->blk_geom->has_uv &&
                              context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1 &&
                              !context_ptr->md_staging_skip_chroma_pred)
                                   ? (int)MAX_MB_PLANE
                                   : 1;
    for (uint32_t full_loop_candidate_index = 0;
         full_loop_candidate_index < context_ptr->md_stage_1_count[context_ptr->target_class];
         ++full_loop_candidate_index) {
        uint32_t cand_index =
            context_ptr->cand_buff_indices[context_ptr->target_class][full_loop_candidate_index];
        ModeDecisionCandidateBuffer *candidate_buffer = candidate_buffer_ptr_array[cand_index];
        ModeDecisionCandidate       *candidate_ptr    = candidate_buffer->candidate_ptr;
        // Use RDOQ in MDS1 for intra candidates if MDS2 is bypassed (otherwise RDOQ is already used for intra candidates in MDS2)
        context_ptr->md_staging_skip_rdoq = (context_ptr->bypass_md_stage_2 &&
                                             is_intra_mode(candidate_ptr->pred_mode))
            ? FALSE
            : TRUE;
        full_loop_core(pcs_ptr,
                       blk_ptr,
                       context_ptr,
                       candidate_buffer,
                       candidate_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index);
    }
}
static void md_stage_2(PictureControlSet *pcs_ptr, BlkStruct *blk_ptr,
                       ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                       uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                       uint32_t blk_origin_index, uint32_t blk_chroma_origin_index) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array = &(
        candidate_buffer_ptr_array_base[0]);

    context_ptr->md_staging_perform_inter_pred = (context_ptr->ifs_ctrls.level == IFS_MDS2) ? TRUE
                                                                                            : FALSE;
    context_ptr->md_staging_skip_interpolation_search = (context_ptr->ifs_ctrls.level == IFS_MDS2)
        ? FALSE
        : TRUE;
    context_ptr->md_staging_skip_chroma_pred          = TRUE;
    context_ptr->md_staging_perform_intra_chroma_pred = FALSE;
    context_ptr->end_plane                            = (context_ptr->blk_geom->has_uv &&
                              context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1 &&
                              !context_ptr->md_staging_skip_chroma_pred)
                                   ? (int)MAX_MB_PLANE
                                   : 1;
    // Set MD Staging full_loop_core settings
    for (uint32_t fullLoopCandidateIndex = 0;
         fullLoopCandidateIndex < context_ptr->md_stage_2_count[context_ptr->target_class];
         ++fullLoopCandidateIndex) {
        uint32_t candidateIndex =
            context_ptr->cand_buff_indices[context_ptr->target_class][fullLoopCandidateIndex];
        ModeDecisionCandidateBuffer *candidate_buffer = candidate_buffer_ptr_array[candidateIndex];
        ModeDecisionCandidate       *candidate_ptr    = candidate_buffer->candidate_ptr;
        context_ptr->md_staging_tx_size_mode          = 0;
        context_ptr->md_staging_txt_level             = is_intra_mode(candidate_ptr->pred_mode)
                        ? context_ptr->txt_ctrls.enabled
                        : 0;
        context_ptr->md_staging_skip_rdoq = is_intra_mode(candidate_ptr->pred_mode) ? TRUE : FALSE;
        context_ptr->md_staging_skip_full_chroma = TRUE;

        context_ptr->md_staging_spatial_sse_full_loop_level =
            context_ptr->spatial_sse_full_loop_level;
        context_ptr->md_staging_fast_coeff_est_level = (context_ptr->pd_pass == PD_PASS_1)
            ? 1
            : context_ptr->rate_est_ctrls.pd0_fast_coeff_est_level;
        context_ptr->md_staging_subres_step          = (context_ptr->pd_pass == PD_PASS_1)
                     ? 0
                     : context_ptr->subres_ctrls.step;

        full_loop_core(pcs_ptr,
                       blk_ptr,
                       context_ptr,
                       candidate_buffer,
                       candidate_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index);
    }
}
void update_intra_chroma_mode(ModeDecisionContext         *context_ptr,
                              ModeDecisionCandidateBuffer *candidate_buffer,
                              PictureControlSet           *pcs_ptr) {
    ModeDecisionCandidate *candidate_ptr = candidate_buffer->candidate_ptr;
    int32_t is_inter = (is_inter_mode(candidate_ptr->pred_mode) || candidate_ptr->use_intrabc)
        ? TRUE
        : FALSE;
    if (context_ptr->blk_geom->sq_size < 128) {
        if (context_ptr->blk_geom->has_uv) {
            if (!is_inter) {
                if (candidate_ptr->palette_info == NULL || candidate_ptr->palette_size[0] == 0) {
                    uint32_t intra_chroma_mode;
                    int32_t  angle_delta;
                    if (((context_ptr->best_inter_cost * context_ptr->uv_ctrls.uv_cfl_th) <
                         (context_ptr->best_intra_cost * 100))) {
                        intra_chroma_mode =
                            context_ptr->best_uv_mode[candidate_ptr->pred_mode]
                                                     [MAX_ANGLE_DELTA +
                                                      candidate_ptr->angle_delta[PLANE_TYPE_Y]];
                        angle_delta =
                            context_ptr->best_uv_angle[candidate_ptr->pred_mode]
                                                      [MAX_ANGLE_DELTA +
                                                       candidate_ptr->angle_delta[PLANE_TYPE_Y]];
                    } else {
                        intra_chroma_mode = candidate_ptr->intra_chroma_mode != UV_CFL_PRED
                            ? context_ptr->best_uv_mode[candidate_ptr->pred_mode]
                                                       [MAX_ANGLE_DELTA +
                                                        candidate_ptr->angle_delta[PLANE_TYPE_Y]]
                            : UV_CFL_PRED;
                        angle_delta       = candidate_ptr->intra_chroma_mode != UV_CFL_PRED
                                  ? context_ptr->best_uv_angle[candidate_ptr->pred_mode]
                                                        [MAX_ANGLE_DELTA +
                                                         candidate_ptr->angle_delta[PLANE_TYPE_Y]]
                                  : 0;
                    }
                    // If CFL OFF or not applicable, and intra_chroma_mode used @ md_stage_0() (first stage intra_mode)
                    // and the best independant intra mode (final stage intra_mode) are not matching then the chroma pred
                    // should be re-performed using best independant chroma pred
                    if (candidate_ptr->intra_chroma_mode != UV_CFL_PRED)
                        if (candidate_ptr->intra_chroma_mode != intra_chroma_mode ||
                            candidate_ptr->angle_delta[PLANE_TYPE_UV] != angle_delta) {
                            // Set to TRUE to redo INTRA CHROMA compensation
                            context_ptr->md_staging_perform_intra_chroma_pred = TRUE;
                            // Update fast_chroma_rate
                            candidate_buffer->fast_chroma_rate =
                                context_ptr
                                    ->fast_chroma_rate[candidate_ptr->pred_mode]
                                                      [MAX_ANGLE_DELTA +
                                                       candidate_ptr->angle_delta[PLANE_TYPE_Y]];
                            // Update intra_chroma_mode
                            candidate_ptr->intra_chroma_mode          = intra_chroma_mode;
                            candidate_ptr->angle_delta[PLANE_TYPE_UV] = angle_delta;
                            // Update transform_type_uv
                            FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
                            if (candidate_ptr->intra_chroma_mode == UV_CFL_PRED)
                                candidate_ptr->transform_type_uv = DCT_DCT;
                            else
                                candidate_ptr->transform_type_uv = av1_get_tx_type(
                                    0, // is_inter
                                    candidate_ptr->pred_mode,
                                    (UvPredictionMode)candidate_ptr->intra_chroma_mode,
                                    PLANE_TYPE_UV,
                                    context_ptr->blk_geom->txsize_uv[0][0],
                                    frm_hdr->reduced_tx_set);
                        }
                }
            }
        }
    }
}

static void md_stage_3_light_pd0(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                 EbPictureBufferDesc *input_picture_ptr,
                                 uint32_t input_origin_index, uint32_t blk_origin_index) {
    ModeDecisionCandidateBuffer *candidate_buffer =
        context_ptr->candidate_buffer_ptr_array[context_ptr->mds0_best_idx];
    //TO MOVE ALL OF THESE at SB LEVEL
    context_ptr->md_staging_perform_inter_pred          = 0;
    context_ptr->md_staging_skip_interpolation_search   = 1;
    context_ptr->md_staging_skip_chroma_pred            = FALSE;
    context_ptr->md_staging_tx_size_mode                = 0;
    context_ptr->md_staging_txt_level                   = 0;
    context_ptr->md_staging_skip_full_chroma            = FALSE;
    context_ptr->md_staging_skip_rdoq                   = 0;
    context_ptr->md_staging_spatial_sse_full_loop_level = 0;
    context_ptr->md_staging_fast_coeff_est_level =
        context_ptr->rate_est_ctrls.pd0_fast_coeff_est_level; //ON  !!!

    // For 8x8 blocks, can't use 4x subsampling b/c no 8x2 transform
    context_ptr->md_staging_subres_step = context_ptr->blk_geom->sq_size >= 16
        ? context_ptr->subres_ctrls.step
        : MIN(1, context_ptr->subres_ctrls.step); //ON  !!!

    assert(IMPLIES(context_ptr->md_staging_subres_step == 2, context_ptr->blk_geom->sq_size >= 16));
    assert(IMPLIES(context_ptr->md_staging_subres_step == 1, context_ptr->blk_geom->sq_size >= 8));
    assert_err(IMPLIES(!context_ptr->disallow_4x4, context_ptr->md_staging_subres_step == 0),
               "residual subsampling cannot be used with 4x4 blocks");

    context_ptr->md_staging_perform_intra_chroma_pred = TRUE;

    full_loop_core_light_pd0(pcs_ptr,
                             context_ptr,
                             candidate_buffer,
                             input_picture_ptr,
                             input_origin_index,
                             blk_origin_index);
}
/*
   md stage 3 for light PD1 path
*/
static void md_stage_3_light_pd1(PictureControlSet *pcs_ptr, BlkStruct *blk_ptr,
                                 ModeDecisionContext *context_ptr,
                                 EbPictureBufferDesc *input_picture_ptr, BlockLocation *loc) {
    ModeDecisionCandidateBuffer *candidate_buffer =
        context_ptr->candidate_buffer_ptr_array[context_ptr->mds0_best_idx];
    ModeDecisionCandidate *candidate_ptr = candidate_buffer->candidate_ptr;

    context_ptr->end_plane          = MAX_MB_PLANE;
    context_ptr->uv_intra_comp_only = TRUE;
    // If EncDec is bypassed, disable features affecting the TX that are usually disabled in EncDec
    if (context_ptr->bypass_encdec) {
        context_ptr->rdoq_ctrls.skip_uv      = 0;
        context_ptr->rdoq_ctrls.dct_dct_only = 0;
    }
    context_ptr->md_staging_skip_rdoq            = FALSE;
    context_ptr->md_staging_fast_coeff_est_level = 1;
    context_ptr->md_staging_subres_step          = 0;

    full_loop_core_light_pd1(
        pcs_ptr, blk_ptr, context_ptr, candidate_buffer, candidate_ptr, input_picture_ptr, loc);
}
static void md_stage_3(PictureControlSet *pcs_ptr, BlkStruct *blk_ptr,
                       ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                       uint32_t input_origin_index, uint32_t input_cb_origin_in_index,
                       uint32_t blk_origin_index, uint32_t blk_chroma_origin_index,
                       uint32_t fullCandidateTotalCount) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array = &(
        candidate_buffer_ptr_array_base[0]);
    context_ptr->md_staging_perform_intra_chroma_pred = TRUE;
    context_ptr->md_staging_skip_chroma_pred          = FALSE;
    context_ptr->end_plane                            = (context_ptr->blk_geom->has_uv &&
                              context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1 &&
                              !context_ptr->md_staging_skip_chroma_pred)
                                   ? (int)MAX_MB_PLANE
                                   : 1;
    for (uint32_t full_loop_candidate_index = 0;
         full_loop_candidate_index < fullCandidateTotalCount;
         ++full_loop_candidate_index) {
        uint32_t cand_index = context_ptr->best_candidate_index_array[full_loop_candidate_index];
        ModeDecisionCandidateBuffer *candidate_buffer = candidate_buffer_ptr_array[cand_index];
        ModeDecisionCandidate       *candidate_ptr    = candidate_buffer->candidate_ptr;
        if (candidate_buffer->candidate_ptr->pred_mode == DC_PRED)
            if (context_ptr->scale_palette)
                if (candidate_ptr->palette_info != NULL && candidate_ptr->palette_size[0] > 0)
                    // if MD is done on 8bit( when HBD is 0 and bypass encdec is ON)
                    // Scale  palette colors to 10bit
                    for (uint8_t col = 0; col < candidate_ptr->palette_size[0]; col++)
                        candidate_ptr->palette_info->pmi.palette_colors[col] *= 4;
        // If EncDec is bypassed, disable features affecting the TX that are usually disabled in EncDec
        if (context_ptr->bypass_encdec && context_ptr->pd_pass == PD_PASS_1) {
            context_ptr->pf_ctrls.pf_shape       = DEFAULT_SHAPE;
            context_ptr->rdoq_ctrls.skip_uv      = 0;
            context_ptr->rdoq_ctrls.dct_dct_only = 0;
        }
        // Set MD Staging full_loop_core settings
        context_ptr->md_staging_perform_inter_pred = context_ptr->nic_ctrls.md_staging_mode !=
                MD_STAGING_MODE_0 ||
            (context_ptr->pd_pass == PD_PASS_1 && context_ptr->need_hbd_comp_mds3);
        context_ptr->md_staging_skip_interpolation_search = (context_ptr->ifs_ctrls.level ==
                                                             IFS_MDS3)
            ? FALSE
            : TRUE;
        context_ptr->md_staging_tx_size_mode              = context_ptr->txs_ctrls.enabled &&
            (context_ptr->blk_geom->sq_size >= context_ptr->txs_ctrls.min_sq_size);
        context_ptr->md_staging_txt_level        = context_ptr->txt_ctrls.enabled;
        context_ptr->md_staging_skip_full_chroma = FALSE;
        context_ptr->md_staging_skip_rdoq        = FALSE;
        context_ptr->md_staging_spatial_sse_full_loop_level =
            context_ptr->spatial_sse_full_loop_level;
        context_ptr->md_staging_fast_coeff_est_level = (context_ptr->pd_pass == PD_PASS_1)
            ? 1
            : context_ptr->rate_est_ctrls.pd0_fast_coeff_est_level;
        context_ptr->md_staging_subres_step          = (context_ptr->pd_pass == PD_PASS_1)
                     ? 0
                     : context_ptr->subres_ctrls.step;
        if (context_ptr->uv_ctrls.nd_uv_serach_mode)
            update_intra_chroma_mode(context_ptr, candidate_buffer, pcs_ptr);
        full_loop_core(pcs_ptr,
                       blk_ptr,
                       context_ptr,
                       candidate_buffer,
                       candidate_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index);
    }
}

void move_blk_data(PictureControlSet *pcs, EncDecContext *context_ptr, BlkStruct *src_cu,
                   BlkStruct *dst_cu) {
    dst_cu->palette_size[0] = src_cu->palette_size[0];
    dst_cu->palette_size[1] = src_cu->palette_size[1];
    if (svt_av1_allow_palette(pcs->parent_pcs_ptr->palette_level, context_ptr->blk_geom->bsize)) {
        svt_memcpy(&dst_cu->palette_info->pmi, &src_cu->palette_info->pmi, sizeof(PaletteModeInfo));
        assert(dst_cu->palette_info->color_idx_map != NULL && "palette: Not-Enough-Memory");
        if (dst_cu->palette_info->color_idx_map != NULL)
            svt_memcpy(dst_cu->palette_info->color_idx_map,
                       src_cu->palette_info->color_idx_map,
                       MAX_PALETTE_SQUARE);
        else
            SVT_ERROR("palette: Not-Enough-Memory\n");
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
    svt_memcpy(
        dst_cu->txb_array, src_cu->txb_array, TRANSFORM_UNIT_MAX_COUNT * sizeof(TransformUnit));

    //CHKN PredictionUnit            prediction_unit_array[MAX_NUM_OF_PU_PER_CU];    // 35-bytes * 4 = 140 bytes
    svt_memcpy(dst_cu->prediction_unit_array,
               src_cu->prediction_unit_array,
               MAX_NUM_OF_PU_PER_CU * sizeof(PredictionUnit));

    dst_cu->skip_flag_context    = src_cu->skip_flag_context;
    dst_cu->prediction_mode_flag = src_cu->prediction_mode_flag;
    dst_cu->block_has_coeff      = src_cu->block_has_coeff;
    dst_cu->split_flag_context   = src_cu->split_flag_context;
    dst_cu->qindex               = src_cu->qindex;
    dst_cu->tx_depth             = src_cu->tx_depth;
    dst_cu->split_flag           = src_cu->split_flag;
    dst_cu->skip_mode            = src_cu->skip_mode;

    //CHKN    MacroBlockD*  av1xd;
    // Don't copy if dest. is NULL
    if (dst_cu->av1xd != NULL)
        svt_memcpy(dst_cu->av1xd, src_cu->av1xd, sizeof(MacroBlockD));

    // uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];

    //CHKN int16_t inter_mode_ctx[MODE_CTX_REF_FRAMES];
    svt_memcpy(
        dst_cu->inter_mode_ctx, src_cu->inter_mode_ctx, MODE_CTX_REF_FRAMES * sizeof(int16_t));

    //CHKN uint8_t  drl_index;
    //CHKN PredictionMode               pred_mode;
    dst_cu->drl_index = src_cu->drl_index;
    dst_cu->pred_mode = src_cu->pred_mode;

    //CHKN IntMv  predmv[2];

    svt_memcpy(dst_cu->predmv, src_cu->predmv, 2 * sizeof(IntMv));
    dst_cu->skip_coeff_context = src_cu->skip_coeff_context;
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
    dst_cu->segment_id = src_cu->segment_id;

    //CHKN uint32_t   is_inter_ctx;
    //CHKN uint32_t                     interp_filters;

    dst_cu->is_inter_ctx   = src_cu->is_inter_ctx;
    dst_cu->interp_filters = src_cu->interp_filters;

    dst_cu->part              = src_cu->part;
    dst_cu->mds_idx           = src_cu->mds_idx;
    dst_cu->filter_intra_mode = src_cu->filter_intra_mode;
    dst_cu->use_intrabc       = src_cu->use_intrabc;
    dst_cu->drl_ctx[0]        = src_cu->drl_ctx[0];
    dst_cu->drl_ctx[1]        = src_cu->drl_ctx[1];
    dst_cu->drl_ctx_near[0]   = src_cu->drl_ctx_near[0];
    dst_cu->drl_ctx_near[1]   = src_cu->drl_ctx_near[1];
}
void move_blk_data_redund(PictureControlSet *pcs, ModeDecisionContext *context_ptr,
                          BlkStruct *src_cu, BlkStruct *dst_cu) {
    dst_cu->segment_id = src_cu->segment_id;
    if (svt_av1_allow_palette(pcs->parent_pcs_ptr->palette_level, context_ptr->blk_geom->bsize)) {
        svt_memcpy(&dst_cu->palette_info->pmi, &src_cu->palette_info->pmi, sizeof(PaletteModeInfo));
        svt_memcpy(dst_cu->palette_info->color_idx_map,
                   src_cu->palette_info->color_idx_map,
                   MAX_PALETTE_SQUARE);
    }

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
    svt_memcpy(
        dst_cu->txb_array, src_cu->txb_array, TRANSFORM_UNIT_MAX_COUNT * sizeof(TransformUnit));

    //CHKN PredictionUnit_t            prediction_unit_array[MAX_NUM_OF_PU_PER_CU];    // 35-bytes * 4 = 140 bytes
    svt_memcpy(dst_cu->prediction_unit_array,
               src_cu->prediction_unit_array,
               MAX_NUM_OF_PU_PER_CU * sizeof(PredictionUnit));
    dst_cu->skip_flag_context    = src_cu->skip_flag_context;
    dst_cu->prediction_mode_flag = src_cu->prediction_mode_flag;
    dst_cu->block_has_coeff      = src_cu->block_has_coeff;
    dst_cu->split_flag_context   = src_cu->split_flag_context;
    dst_cu->qindex               = src_cu->qindex;
    dst_cu->skip_mode            = src_cu->skip_mode;
    dst_cu->tx_depth             = src_cu->tx_depth;
    //CHKN    MacroBlockD*  av1xd;
    svt_memcpy(dst_cu->av1xd, src_cu->av1xd, sizeof(MacroBlockD));

    // uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];

    //CHKN int16_t inter_mode_ctx[MODE_CTX_REF_FRAMES];
    svt_memcpy(
        dst_cu->inter_mode_ctx, src_cu->inter_mode_ctx, MODE_CTX_REF_FRAMES * sizeof(int16_t));

    //CHKN uint8_t  drl_index;
    //CHKN PredictionMode               pred_mode;
    dst_cu->drl_index = src_cu->drl_index;
    dst_cu->pred_mode = src_cu->pred_mode;

    //CHKN IntMv  predmv[2];

    svt_memcpy(dst_cu->predmv, src_cu->predmv, 2 * sizeof(IntMv));
    dst_cu->skip_coeff_context = src_cu->skip_coeff_context;
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

    //CHKN uint32_t   is_inter_ctx;
    //CHKN uint32_t                     interp_filters;

    dst_cu->is_inter_ctx   = src_cu->is_inter_ctx;
    dst_cu->interp_filters = src_cu->interp_filters;

    dst_cu->part            = src_cu->part;
    dst_cu->use_intrabc     = src_cu->use_intrabc;
    dst_cu->drl_ctx[0]      = src_cu->drl_ctx[0];
    dst_cu->drl_ctx[1]      = src_cu->drl_ctx[1];
    dst_cu->drl_ctx_near[0] = src_cu->drl_ctx_near[0];
    dst_cu->drl_ctx_near[1] = src_cu->drl_ctx_near[1];
    for (int list_idx = 0; list_idx < MAX_NUM_OF_REF_PIC_LIST; list_idx++) {
        for (int ref_idx = 0; ref_idx < MAX_REF_IDX; ref_idx++) {
            context_ptr->sb_me_mv[dst_cu->mds_idx][list_idx][ref_idx][0] =
                context_ptr->sb_me_mv[src_cu->mds_idx][list_idx][ref_idx][0];
            context_ptr->sb_me_mv[dst_cu->mds_idx][list_idx][ref_idx][1] =
                context_ptr->sb_me_mv[src_cu->mds_idx][list_idx][ref_idx][1];
        }
    }
}

void check_redundant_block(const BlockGeom *blk_geom, ModeDecisionContext *context_ptr,
                           uint8_t *redundant_blk_avail, uint16_t *redundant_blk_mds) {
    if (blk_geom->redund) {
        for (int it = 0; it < blk_geom->redund_list.list_size; it++) {
            if (context_ptr->avail_blk_flag[blk_geom->redund_list.blk_mds_table[it]]) {
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
            if (context_ptr->avail_blk_flag[blk_geom->similar_list.blk_mds_table[it]]) {
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
    context_ptr->spatial_sse_full_loop_level =
        context_ptr->spatial_sse_ctrls.spatial_sse_full_loop_level;
    UNUSED(pcs);
    // set compound_types_to_try
    set_inter_comp_controls(context_ptr, context_ptr->inter_compound_mode);
    return return_error;
}

void init_chroma_mode(ModeDecisionContext *context_ptr) {
    Bool use_angle_delta = av1_use_angle_delta(context_ptr->blk_geom->bsize,
                                               context_ptr->intra_ctrls.angular_pred_level);
    for (uint8_t intra_mode = DC_PRED; intra_mode <= PAETH_PRED; ++intra_mode) {
        uint8_t angleDeltaCandidateCount = (use_angle_delta &&
                                            av1_is_directional_mode((PredictionMode)intra_mode))
            ? 7
            : 1;
        uint8_t angle_delta_shift        = 1;
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
}

static INLINE void rtime_alloc_uv_cand_buff_indices(uint32_t **uv_cand_buff_indices,
                                                    uint32_t   max_nics_uv) {
    (*uv_cand_buff_indices) = (uint32_t *)malloc(max_nics_uv * sizeof(*uv_cand_buff_indices));
}
/*
Perform search for the best chroma mode (intra modes only).  The search involves
the following main parts:

1. Prepare all the chroma candidates to be tested
2. Perform compensation and compute the distortion for each candidate
3. Sort the candidates, and for the best n candidates, perform the full loop (TX, quant, inv. quant)
4. Compute the full cost for the remaining candidates and select the best chroma mode (to be combined
   with luma modes in future MD stages).

*/
static void search_best_independent_uv_mode(PictureControlSet   *pcs,
                                            EbPictureBufferDesc *input_picture_ptr,
                                            uint32_t             input_cb_origin_in_index,
                                            uint32_t             input_cr_origin_in_index,
                                            uint32_t             cu_chroma_origin_index,
                                            ModeDecisionContext *ctx) {
    PictureParentControlSet *ppcs    = pcs->parent_pcs_ptr;
    FrameHeader             *frm_hdr = &ppcs->frm_hdr;
    uint32_t full_lambda = ctx->full_lambda_md[ctx->hbd_mode_decision ? EB_10_BIT_MD : EB_8_BIT_MD];

    uint64_t coeff_rate[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
    uint64_t distortion[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];

    ModeDecisionCandidate *candidate_array         = ctx->fast_candidate_array;
    uint32_t               start_fast_buffer_index = ppcs->max_can_count;
    uint32_t               start_full_buffer_index = ctx->max_nics;
    unsigned int           uv_mode_total_count     = start_fast_buffer_index;
    UvPredictionMode       uv_mode_end     = (UvPredictionMode)ctx->intra_ctrls.intra_mode_end;
    uint8_t                uv_mode_start   = UV_DC_PRED;
    Bool                   use_angle_delta = av1_use_angle_delta(ctx->blk_geom->bsize,
                                               ctx->intra_ctrls.angular_pred_level);
    uint8_t                disable_angle_prediction = (ctx->intra_ctrls.angular_pred_level == 0);
    const int              uv_angle_delta_shift     = 1;
    uint8_t                directional_mode_skip_mask[INTRA_MODES] = {0};
    // For aggressive angular levels, don't test angular candidate for certain modes
    if (ctx->intra_ctrls.angular_pred_level >= 4) {
        for (uint8_t i = D45_PRED; i < INTRA_MODE_END; i++) directional_mode_skip_mask[i] = 1;
    }

    // Prepare chroma candidates to be tested
    for (UvPredictionMode uv_mode = uv_mode_start; uv_mode <= uv_mode_end; ++uv_mode) {
        // If mode is not directional, or is enabled directional mode, proceed with injection
        if (!av1_is_directional_mode((PredictionMode)uv_mode) ||
            (!disable_angle_prediction &&
             directional_mode_skip_mask[(PredictionMode)uv_mode] == 0)) {
            int uv_angle_delta_candidate_count =
                (use_angle_delta && av1_is_directional_mode((PredictionMode)uv_mode) &&
                 ctx->intra_ctrls.angular_pred_level <= 2)
                ? 7
                : 1;

            for (int uv_angle_delta_counter = 0;
                 uv_angle_delta_counter < uv_angle_delta_candidate_count;
                 ++uv_angle_delta_counter) {
                int32_t uv_angle_delta = CLIP(
                    uv_angle_delta_shift *
                        (uv_angle_delta_candidate_count == 1
                             ? 0
                             : uv_angle_delta_counter - (uv_angle_delta_candidate_count >> 1)),
                    -MAX_ANGLE_DELTA,
                    MAX_ANGLE_DELTA);
                if (ctx->intra_ctrls.angular_pred_level >= 2 &&
                    (uv_angle_delta == -1 || uv_angle_delta == 1 || uv_angle_delta == -2 ||
                     uv_angle_delta == 2))
                    continue;
                candidate_array[uv_mode_total_count].use_intrabc                = 0;
                candidate_array[uv_mode_total_count].angle_delta[PLANE_TYPE_UV] = 0;
                candidate_array[uv_mode_total_count].pred_mode                  = DC_PRED;
                candidate_array[uv_mode_total_count].intra_chroma_mode          = uv_mode;
                candidate_array[uv_mode_total_count].angle_delta[PLANE_TYPE_UV] = uv_angle_delta;
                candidate_array[uv_mode_total_count].tx_depth                   = 0;
                candidate_array[uv_mode_total_count].palette_info               = NULL;
                candidate_array[uv_mode_total_count].filter_intra_mode = FILTER_INTRA_MODES;
                candidate_array[uv_mode_total_count].cfl_alpha_signs   = 0;
                candidate_array[uv_mode_total_count].cfl_alpha_idx     = 0;
                candidate_array[uv_mode_total_count].transform_type[0] = DCT_DCT;
                candidate_array[uv_mode_total_count].ref_frame_type    = INTRA_FRAME;
                candidate_array[uv_mode_total_count].motion_mode       = SIMPLE_TRANSLATION;
                candidate_array[uv_mode_total_count].transform_type_uv = av1_get_tx_type(
                    0, // is_inter
                    (PredictionMode)0,
                    (UvPredictionMode)uv_mode,
                    PLANE_TYPE_UV,
                    ctx->blk_geom->txsize_uv[0][0],
                    frm_hdr->reduced_tx_set);
                uv_mode_total_count++;
            }
        }
    }
    uv_mode_total_count = uv_mode_total_count - start_fast_buffer_index;

    // Prepare fast-loop search settings
    ctx->md_staging_skip_rdoq        = 0;
    ctx->uv_intra_comp_only          = TRUE;
    ctx->md_staging_skip_chroma_pred = FALSE;
    ctx->end_plane = (ctx->blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1)
        ? (int)MAX_MB_PLANE
        : 1;
    assert(ctx->md_staging_skip_chroma_pred == FALSE);

    // Perform fast-loop search for all candidates
    for (unsigned int uv_mode_count = 0; uv_mode_count < uv_mode_total_count; uv_mode_count++) {
        ModeDecisionCandidateBuffer *candidate_buffer =
            ctx->candidate_buffer_ptr_array[uv_mode_count + start_full_buffer_index];
        candidate_buffer->candidate_ptr =
            &ctx->fast_candidate_array[uv_mode_count + start_fast_buffer_index];
        svt_product_prediction_fun_table[is_inter_mode(candidate_buffer->candidate_ptr->pred_mode)](
            ctx->hbd_mode_decision, ctx, pcs, candidate_buffer);
        uint32_t chroma_fast_distortion;
        if (ctx->mds0_ctrls.mds0_dist_type == MDS0_VAR) {
            if (!ctx->hbd_mode_decision) {
                const AomVarianceFnPtr *fn_ptr = &mefn_ptr[ctx->blk_geom->bsize_uv];
                unsigned int            sse;
                uint8_t                *pred_cb = candidate_buffer->prediction_ptr->buffer_cb +
                    cu_chroma_origin_index;
                uint8_t *src_cb        = input_picture_ptr->buffer_cb + input_cb_origin_in_index;
                chroma_fast_distortion = fn_ptr->vf(pred_cb,
                                                    candidate_buffer->prediction_ptr->stride_cb,
                                                    src_cb,
                                                    input_picture_ptr->stride_cb,
                                                    &sse) >>
                    2;
                uint8_t *pred_cr = candidate_buffer->prediction_ptr->buffer_cr +
                    cu_chroma_origin_index;
                uint8_t *src_cr = input_picture_ptr->buffer_cr + input_cr_origin_in_index;
                chroma_fast_distortion += (fn_ptr->vf(pred_cr,
                                                      candidate_buffer->prediction_ptr->stride_cr,
                                                      src_cr,
                                                      input_picture_ptr->stride_cr,
                                                      &sse) >>
                                           2);
            } else {
                const AomVarianceFnPtr *fn_ptr = &mefn_ptr[ctx->blk_geom->bsize_uv];
                unsigned int            sse;
                uint16_t *pred_cb = ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cb) +
                    cu_chroma_origin_index;
                uint16_t *src_cb = ((uint16_t *)input_picture_ptr->buffer_cb) +
                    input_cb_origin_in_index;
                chroma_fast_distortion = fn_ptr->vf_hbd_10(
                                             CONVERT_TO_BYTEPTR(pred_cb),
                                             candidate_buffer->prediction_ptr->stride_cb,
                                             CONVERT_TO_BYTEPTR(src_cb),
                                             input_picture_ptr->stride_cb,
                                             &sse) >>
                    1;
                uint16_t *pred_cr = ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cr) +
                    cu_chroma_origin_index;
                uint16_t *src_cr = ((uint16_t *)input_picture_ptr->buffer_cr) +
                    input_cr_origin_in_index;
                chroma_fast_distortion += (fn_ptr->vf_hbd_10(
                                               CONVERT_TO_BYTEPTR(pred_cr),
                                               candidate_buffer->prediction_ptr->stride_cr,
                                               CONVERT_TO_BYTEPTR(src_cr),
                                               input_picture_ptr->stride_cr,
                                               &sse) >>
                                           1);
            }
        } else if (!ctx->hbd_mode_decision) {
            chroma_fast_distortion = svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_cb + input_cb_origin_in_index,
                input_picture_ptr->stride_cb,
                candidate_buffer->prediction_ptr->buffer_cb + cu_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cb,
                ctx->blk_geom->bheight_uv,
                ctx->blk_geom->bwidth_uv);

            chroma_fast_distortion += svt_nxm_sad_kernel_sub_sampled(
                input_picture_ptr->buffer_cr + input_cr_origin_in_index,
                input_picture_ptr->stride_cr,
                candidate_buffer->prediction_ptr->buffer_cr + cu_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cr,
                ctx->blk_geom->bheight_uv,
                ctx->blk_geom->bwidth_uv);
        } else {
            chroma_fast_distortion = sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_cb) + input_cb_origin_in_index,
                input_picture_ptr->stride_cb,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cb) + cu_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cb,
                ctx->blk_geom->bheight_uv,
                ctx->blk_geom->bwidth_uv);

            chroma_fast_distortion += sad_16b_kernel(
                ((uint16_t *)input_picture_ptr->buffer_cr) + input_cr_origin_in_index,
                input_picture_ptr->stride_cr,
                ((uint16_t *)candidate_buffer->prediction_ptr->buffer_cr) + cu_chroma_origin_index,
                candidate_buffer->prediction_ptr->stride_cr,
                ctx->blk_geom->bheight_uv,
                ctx->blk_geom->bwidth_uv);
        }
        // Do not consider rate @ this stage
        *(candidate_buffer->fast_cost_ptr) = chroma_fast_distortion;
    }

    // Sort uv_mode candidates (in terms of distortion only)
    uint32_t *uv_cand_buff_indices;
    rtime_alloc_uv_cand_buff_indices(&uv_cand_buff_indices, ctx->max_nics_uv);
    memset(uv_cand_buff_indices, 0xFF, ctx->max_nics_uv * sizeof(*uv_cand_buff_indices));

    sort_fast_cost_based_candidates(
        ctx,
        start_full_buffer_index,
        uv_mode_total_count, //how many cand buffers to sort. one of the buffers can have max cost.
        uv_cand_buff_indices);

    // Reset *(candidate_buffer->fast_cost_ptr)
    for (unsigned int uv_mode_count = 0; uv_mode_count < uv_mode_total_count; uv_mode_count++) {
        ModeDecisionCandidateBuffer *candidate_buffer =
            ctx->candidate_buffer_ptr_array[uv_mode_count + start_full_buffer_index];
        *(candidate_buffer->fast_cost_ptr) = MAX_CU_COST;
    }

    // Set number of UV candidates to be tested in the full loop
    unsigned int uv_mode_nfl_count = pcs->slice_type == I_SLICE ? 64
        : ppcs->is_used_as_reference_flag                       ? 32
                                                                : 16;
    uv_mode_nfl_count              = MAX(
        1, DIVIDE_AND_ROUND(uv_mode_nfl_count * ctx->uv_ctrls.uv_nic_scaling_num, 16));
    uv_mode_nfl_count = MIN(uv_mode_nfl_count, uv_mode_total_count);
    uv_mode_nfl_count = MAX(uv_mode_nfl_count, 1);

    // Full-loop search uv_mode
    for (unsigned int uv_mode_count = 0;
         uv_mode_count < MIN(uv_mode_total_count, uv_mode_nfl_count);
         uv_mode_count++) {
        ModeDecisionCandidateBuffer *candidate_buffer =
            ctx->candidate_buffer_ptr_array[uv_cand_buff_indices[uv_mode_count]];
        ModeDecisionCandidate *cand = candidate_buffer->candidate_ptr =
            &ctx->fast_candidate_array[uv_cand_buff_indices[uv_mode_count] -
                                       start_full_buffer_index + start_fast_buffer_index];
        uint16_t cb_qindex                           = ctx->qp_index;
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
                        ctx->hbd_mode_decision,
                        ctx->blk_geom->bwidth_uv,
                        ctx->blk_geom->bheight_uv);

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
                        ctx->hbd_mode_decision,
                        ctx->blk_geom->bwidth_uv,
                        ctx->blk_geom->bheight_uv);
        svt_aom_full_loop_uv(pcs,
                             ctx,
                             candidate_buffer,
                             input_picture_ptr,
                             COMPONENT_CHROMA,
                             cb_qindex,
                             count_non_zero_coeffs,
                             cb_full_distortion,
                             cr_full_distortion,
                             &cb_coeff_bits,
                             &cr_coeff_bits,
                             1);

        coeff_rate[cand->intra_chroma_mode][MAX_ANGLE_DELTA + cand->angle_delta[PLANE_TYPE_UV]] =
            cb_coeff_bits + cr_coeff_bits;
        distortion[cand->intra_chroma_mode][MAX_ANGLE_DELTA + cand->angle_delta[PLANE_TYPE_UV]] =
            cb_full_distortion[DIST_CALC_RESIDUAL] + cr_full_distortion[DIST_CALC_RESIDUAL];
    }

    // Loop over all intra modes, then loop over all uv_modes to derive the best uv_mode for a given intra mode (in term of rate)
    uint8_t   intra_mode_end    = ctx->intra_ctrls.intra_mode_end;
    const int angle_delta_shift = 1;
    // intra_mode loop (luma mode loop)
    for (PredictionMode intra_mode = DC_PRED; intra_mode <= intra_mode_end; ++intra_mode) {
        // If mode is not directional, or is enabled directional mode, proceed with injection
        if (!av1_is_directional_mode(intra_mode) ||
            (!disable_angle_prediction && directional_mode_skip_mask[intra_mode] == 0)) {
            int angle_delta_candidate_count = (use_angle_delta &&
                                               av1_is_directional_mode(intra_mode) &&
                                               ctx->intra_ctrls.angular_pred_level <= 2)
                ? 7
                : 1;

            for (int angle_delta_counter = 0; angle_delta_counter < angle_delta_candidate_count;
                 ++angle_delta_counter) {
                int32_t angle_delta = CLIP(
                    angle_delta_shift *
                        (angle_delta_candidate_count == 1
                             ? 0
                             : angle_delta_counter - (angle_delta_candidate_count >> 1)),
                    -MAX_ANGLE_DELTA,
                    MAX_ANGLE_DELTA);

                if (ctx->intra_ctrls.angular_pred_level >= 2 &&
                    (angle_delta == -1 || angle_delta == 1 || angle_delta == -2 ||
                     angle_delta == 2))
                    continue;

                // uv mode loop
                ctx->best_uv_cost[intra_mode][MAX_ANGLE_DELTA + angle_delta] = (uint64_t)~0;
                for (unsigned int uv_mode_count = 0;
                     uv_mode_count < MIN(uv_mode_total_count, uv_mode_nfl_count);
                     uv_mode_count++) {
                    ModeDecisionCandidateBuffer *candidate_buffer =
                        ctx->candidate_buffer_ptr_array[uv_cand_buff_indices[uv_mode_count]];
                    ModeDecisionCandidate *candidate_ptr = &(
                        ctx->fast_candidate_array[uv_cand_buff_indices[uv_mode_count] -
                                                  start_full_buffer_index +
                                                  start_fast_buffer_index]);
                    candidate_ptr->angle_delta[PLANE_TYPE_Y] = angle_delta;
                    candidate_ptr->pred_mode                 = intra_mode;

                    // Fast Cost
                    av1_product_fast_cost_func_table[is_inter_mode(candidate_ptr->pred_mode)](
                        ctx,
                        ctx->blk_ptr,
                        candidate_buffer,
                        NOT_USED_VALUE,
                        0,
                        0,
                        0,
                        pcs,
                        &(ctx->md_local_blk_unit[ctx->blk_geom->blkidx_mds]
                              .ed_ref_mv_stack[candidate_ptr->ref_frame_type][0]),
                        ctx->blk_geom,
                        ctx->blk_origin_y >> MI_SIZE_LOG2,
                        ctx->blk_origin_x >> MI_SIZE_LOG2,
                        ctx->inter_intra_comp_ctrls.enabled,
                        ctx->intra_luma_left_mode,
                        ctx->intra_luma_top_mode);
                    uint64_t rate =
                        coeff_rate[candidate_ptr->intra_chroma_mode]
                                  [MAX_ANGLE_DELTA + candidate_ptr->angle_delta[PLANE_TYPE_UV]] +
                        candidate_buffer->fast_luma_rate + candidate_buffer->fast_chroma_rate;
                    uint64_t uv_cost = RDCOST(
                        full_lambda,
                        rate,
                        distortion[candidate_ptr->intra_chroma_mode]
                                  [MAX_ANGLE_DELTA + candidate_ptr->angle_delta[PLANE_TYPE_UV]]);

                    if (uv_cost < ctx->best_uv_cost[intra_mode][MAX_ANGLE_DELTA + angle_delta]) {
                        ctx->best_uv_mode[intra_mode][MAX_ANGLE_DELTA + angle_delta] =
                            candidate_ptr->intra_chroma_mode;
                        ctx->best_uv_angle[intra_mode][MAX_ANGLE_DELTA + angle_delta] =
                            candidate_ptr->angle_delta[PLANE_TYPE_UV];

                        ctx->best_uv_cost[intra_mode][MAX_ANGLE_DELTA + angle_delta] = uv_cost;
                        ctx->fast_luma_rate[intra_mode][MAX_ANGLE_DELTA + angle_delta] =
                            candidate_buffer->fast_luma_rate;
                        ctx->fast_chroma_rate[intra_mode][MAX_ANGLE_DELTA + angle_delta] =
                            candidate_buffer->fast_chroma_rate;
                    }
                }
            }
        }
    }

    free(uv_cand_buff_indices);
}

// Perform the NIC class pruning and candidate pruning after MSD0
void post_mds0_nic_pruning(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                           uint64_t best_md_stage_cost, uint64_t best_md_stage_dist) {
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL;
         cand_class_it++) {
        uint64_t mds1_cand_th = is_intra_class(cand_class_it)
            ? context_ptr->nic_ctrls.pruning_ctrls.mds1_cand_base_th_intra
            : context_ptr->nic_ctrls.pruning_ctrls.mds1_cand_base_th_inter;
        if (mds1_cand_th != (uint64_t)~0 ||
            context_ptr->nic_ctrls.pruning_ctrls.mds1_class_th != (uint64_t)~0)
            if (context_ptr->md_stage_0_count[cand_class_it] > 0 &&
                context_ptr->md_stage_1_count[cand_class_it] > 0) {
                uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
                uint64_t  class_best_cost   = *(
                    context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->fast_cost_ptr);
                // inter class pruning
                if (class_best_cost && best_md_stage_cost &&
                    (class_best_cost != best_md_stage_cost)) {
                    if (context_ptr->nic_ctrls.pruning_ctrls.mds1_class_th == 0) {
                        context_ptr->md_stage_1_count[cand_class_it] = 0;
                        continue;
                    }
                    uint64_t dev = ((class_best_cost - best_md_stage_cost) * 100) /
                        best_md_stage_cost;
                    if (dev) {
                        if (dev >= context_ptr->nic_ctrls.pruning_ctrls.mds1_class_th) {
                            context_ptr->md_stage_1_count[cand_class_it] = 0;
                            continue;
                        } else if (context_ptr->nic_ctrls.pruning_ctrls.mds1_band_cnt >= 3 &&
                                   context_ptr->md_stage_1_count[cand_class_it] > 1) {
                            uint8_t band_idx =
                                (uint8_t)(dev *
                                          (context_ptr->nic_ctrls.pruning_ctrls.mds1_band_cnt - 1) /
                                          context_ptr->nic_ctrls.pruning_ctrls.mds1_class_th);
                            context_ptr->md_stage_1_count[cand_class_it] = DIVIDE_AND_ROUND(
                                context_ptr->md_stage_1_count[cand_class_it], band_idx + 1);
                        }
                    }
                }
                // intra class pruning
                uint32_t cand_count = 1;
                if (class_best_cost) {
                    uint16_t mds1_cand_th_rank_factor =
                        context_ptr->nic_ctrls.pruning_ctrls.mds1_cand_th_rank_factor;
                    while (
                        cand_count < context_ptr->md_stage_1_count[cand_class_it] &&
                        ((((*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[cand_count]]
                                  ->fast_cost_ptr) -
                            class_best_cost) *
                           100) /
                          class_best_cost) <
                         (mds1_cand_th /
                          (mds1_cand_th_rank_factor ? (mds1_cand_th_rank_factor * cand_count)
                                                    : 1)))) {
                        cand_count++;
                    }
                }
                context_ptr->md_stage_1_count[cand_class_it] = cand_count;
            }
        context_ptr->md_stage_1_total_count += context_ptr->md_stage_1_count[cand_class_it];
    }

    // skip MDS1 ONLY when we have only one candidate OR that ((candidate cost) / QP) > TH
    if (context_ptr->nic_ctrls.pruning_ctrls.enable_skipping_mds1) {
        uint64_t th_normalizer = context_ptr->blk_geom->bheight * context_ptr->blk_geom->bwidth *
            (pcs_ptr->picture_qp >> 1);

        if (context_ptr->md_stage_1_total_count > 1) {
            if ((100 * best_md_stage_dist) <
                (context_ptr->nic_ctrls.pruning_ctrls.force_1_cand_th * th_normalizer)) {
                context_ptr->perform_mds1 = 0;

                for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL;
                     cand_class_it++) {
                    context_ptr->md_stage_1_count[cand_class_it] = (cand_class_it ==
                                                                    context_ptr->mds0_best_class_it)
                        ? 1
                        : 0;
                }
                context_ptr->md_stage_1_total_count = 1;
            }
        } else {
            context_ptr->perform_mds1 = 0;
        }
    }
}
// Perform the NIC class pruning and candidate pruning after MSD1
void post_mds1_nic_pruning(ModeDecisionContext *context_ptr, uint64_t best_md_stage_cost) {
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL;
         cand_class_it++) {
        uint64_t mds2_cand_th = context_ptr->nic_ctrls.pruning_ctrls.mds2_cand_base_th;
        if (mds2_cand_th != (uint64_t)~0 ||
            context_ptr->nic_ctrls.pruning_ctrls.mds2_class_th != (uint64_t)~0)
            if (context_ptr->md_stage_1_count[cand_class_it] > 0 &&
                context_ptr->md_stage_2_count[cand_class_it] > 0 &&
                context_ptr->bypass_md_stage_1 == FALSE) {
                uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
                uint64_t  class_best_cost   = *(
                    context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->full_cost_ptr);

                // class pruning
                if (class_best_cost && best_md_stage_cost &&
                    (class_best_cost != best_md_stage_cost)) {
                    if (context_ptr->nic_ctrls.pruning_ctrls.mds2_class_th == 0) {
                        context_ptr->md_stage_2_count[cand_class_it] = 0;
                        continue;
                    }
                    uint64_t dev = ((class_best_cost - best_md_stage_cost) * 100) /
                        best_md_stage_cost;
                    if (dev) {
                        if (dev >= context_ptr->nic_ctrls.pruning_ctrls.mds2_class_th) {
                            context_ptr->md_stage_2_count[cand_class_it] = 0;
                            continue;
                        } else if (context_ptr->nic_ctrls.pruning_ctrls.mds2_band_cnt >= 3 &&
                                   context_ptr->md_stage_2_count[cand_class_it] > 1) {
                            uint8_t band_idx =
                                (uint8_t)(dev *
                                          (context_ptr->nic_ctrls.pruning_ctrls.mds2_band_cnt - 1) /
                                          context_ptr->nic_ctrls.pruning_ctrls.mds2_class_th);
                            context_ptr->md_stage_2_count[cand_class_it] = DIVIDE_AND_ROUND(
                                context_ptr->md_stage_2_count[cand_class_it], band_idx + 1);
                        }
                    }
                }
                // intra class pruning
                // candidate pruning
                if (context_ptr->md_stage_2_count[cand_class_it] > 0) {
                    uint32_t cand_count = 1;
                    if (class_best_cost &&
                        cand_count < context_ptr->md_stage_2_count[cand_class_it]) {
                        uint16_t mds2_cand_th_rank_factor =
                            context_ptr->nic_ctrls.pruning_ctrls.mds2_cand_th_rank_factor;
                        // When enabled, modify the rank factor based on info from previous MD stages
                        if (mds2_cand_th_rank_factor) {
                            if (cand_class_it != context_ptr->mds1_best_class_it)
                                mds2_cand_th_rank_factor += 3;
                            else if (context_ptr->mds0_best_idx == context_ptr->mds1_best_idx)
                                mds2_cand_th_rank_factor += 2;
                        }
                        uint16_t mds2_relative_dev_th =
                            context_ptr->nic_ctrls.pruning_ctrls.mds2_relative_dev_th;
                        uint64_t dev =
                            (((*(context_ptr
                                     ->candidate_buffer_ptr_array[cand_buff_indices[cand_count]]
                                     ->full_cost_ptr) -
                               class_best_cost) *
                              100) /
                             class_best_cost);
                        uint64_t prev_dev = dev;
                        while (
                            (!mds2_relative_dev_th || (dev <= (prev_dev + mds2_relative_dev_th))) &&
                            (dev <
                             (mds2_cand_th /
                              (mds2_cand_th_rank_factor ? (mds2_cand_th_rank_factor * cand_count)
                                                        : 1)))) {
                            cand_count++;
                            // Break out of loop if reached max cand_count to avoid accessing unallocated candidate buffer
                            if (cand_count >= context_ptr->md_stage_2_count[cand_class_it])
                                break;
                            prev_dev = dev;
                            dev =
                                (((*(context_ptr
                                         ->candidate_buffer_ptr_array[cand_buff_indices[cand_count]]
                                         ->full_cost_ptr) -
                                   class_best_cost) *
                                  100) /
                                 class_best_cost);
                        }
                    }
                    context_ptr->md_stage_2_count[cand_class_it] = cand_count;
                }
            }
        context_ptr->md_stage_2_total_count += context_ptr->md_stage_2_count[cand_class_it];
    }
}
// Perform the NIC class pruning and candidate pruning after MSD2
void post_mds2_nic_pruning(ModeDecisionContext *context_ptr, uint64_t best_md_stage_cost) {
    context_ptr->md_stage_3_total_count = 0;
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL;
         cand_class_it++) {
        uint64_t mds3_cand_th = context_ptr->nic_ctrls.pruning_ctrls.mds3_cand_base_th;
        if (mds3_cand_th != (uint64_t)~0 ||
            context_ptr->nic_ctrls.pruning_ctrls.mds3_class_th != (uint64_t)~0)
            if (context_ptr->md_stage_2_count[cand_class_it] > 0 &&
                context_ptr->md_stage_3_count[cand_class_it] > 0 &&
                context_ptr->bypass_md_stage_2 == FALSE) {
                uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
                uint64_t  class_best_cost   = *(
                    context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->full_cost_ptr);

                // inter class pruning
                if (class_best_cost && best_md_stage_cost &&
                    (class_best_cost != best_md_stage_cost)) {
                    if (context_ptr->nic_ctrls.pruning_ctrls.mds3_class_th == 0) {
                        context_ptr->md_stage_3_count[cand_class_it] = 0;
                        continue;
                    }
                    uint64_t dev = ((class_best_cost - best_md_stage_cost) * 100) /
                        best_md_stage_cost;
                    if (dev) {
                        if (dev >= context_ptr->nic_ctrls.pruning_ctrls.mds3_class_th) {
                            context_ptr->md_stage_3_count[cand_class_it] = 0;
                            continue;
                        } else if (context_ptr->nic_ctrls.pruning_ctrls.mds3_band_cnt >= 3 &&
                                   context_ptr->md_stage_3_count[cand_class_it] > 1) {
                            uint8_t band_idx =
                                (uint8_t)(dev *
                                          (context_ptr->nic_ctrls.pruning_ctrls.mds3_band_cnt - 1) /
                                          context_ptr->nic_ctrls.pruning_ctrls.mds3_class_th);
                            context_ptr->md_stage_3_count[cand_class_it] = DIVIDE_AND_ROUND(
                                context_ptr->md_stage_3_count[cand_class_it], band_idx + 1);
                        }
                    }
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
                          class_best_cost) < mds3_cand_th)) {
                        cand_count++;
                    }
                context_ptr->md_stage_3_count[cand_class_it] = cand_count;
            }
        context_ptr->md_stage_3_total_count += context_ptr->md_stage_3_count[cand_class_it];
    }
}
int      av1_get_reference_mode_context_new(const MacroBlockD *xd);
uint64_t estimate_ref_frame_type_bits(ModeDecisionContext *ctx, BlkStruct *blk_ptr,
                                      uint8_t ref_frame_type, Bool is_compound);
/*
 * Estimate the rate of signaling all available ref_frame_type
 */
void estimate_ref_frames_num_bits(struct ModeDecisionContext *context_ptr,
                                  PictureControlSet          *pcs_ptr) {
    uint64_t     comp_inter_fac_bits_uni = 0;
    uint64_t     comp_inter_fac_bits_bi  = 0;
    FrameHeader *frm_hdr                 = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    // does the feature use compound prediction or not
    // (if not specified at the frame/segment level)
    if (frm_hdr->reference_mode == REFERENCE_MODE_SELECT) {
        if (MIN(context_ptr->blk_geom->bwidth, context_ptr->blk_geom->bheight) >= 8) {
            int32_t reference_mode_context;
            //aom_write_symbol(w, is_compound, av1_get_reference_mode_cdf(blk_ptr->av1xd), 2);
            reference_mode_context = av1_get_reference_mode_context_new(
                context_ptr->blk_ptr->av1xd);
            comp_inter_fac_bits_uni =
                context_ptr->md_rate_estimation_ptr->comp_inter_fac_bits[reference_mode_context][0];
            comp_inter_fac_bits_bi =
                context_ptr->md_rate_estimation_ptr->comp_inter_fac_bits[reference_mode_context][1];
        }
    }
    for (uint32_t ref_it = 0; ref_it < pcs_ptr->parent_pcs_ptr->tot_ref_frame_types; ++ref_it) {
        MvReferenceFrame ref_pair = pcs_ptr->parent_pcs_ptr->ref_frame_type_arr[ref_it];
        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_pair);

        //single ref/list
        if (rf[1] == NONE_FRAME) {
            MvReferenceFrame ref_frame_type = rf[0];
            context_ptr->estimate_ref_frames_num_bits[ref_frame_type] =
                estimate_ref_frame_type_bits(context_ptr, context_ptr->blk_ptr, ref_frame_type, 0) +
                comp_inter_fac_bits_uni;
        } else {
            context_ptr->estimate_ref_frames_num_bits[ref_pair] =
                estimate_ref_frame_type_bits(context_ptr, context_ptr->blk_ptr, ref_pair, 1) +
                comp_inter_fac_bits_bi;
        }
    }
}

void calc_scr_to_recon_dist_per_quadrant(ModeDecisionContext         *context_ptr,
                                         EbPictureBufferDesc         *input_picture_ptr,
                                         const uint32_t               input_origin_index,
                                         const uint32_t               input_cb_origin_in_index,
                                         ModeDecisionCandidateBuffer *candidate_buffer,
                                         const uint32_t               blk_origin_index,
                                         const uint32_t               blk_chroma_origin_index) {
    // if a non-4x4 SQ
    if ((context_ptr->blk_geom->bwidth == context_ptr->blk_geom->bheight) &&
        context_ptr->blk_geom->sq_size > 4) {
        EbPictureBufferDesc *recon_ptr = candidate_buffer->recon_ptr;

        EbSpatialFullDistType spatial_full_dist_type_fun = context_ptr->hbd_mode_decision
            ? svt_full_distortion_kernel16_bits
            : svt_spatial_full_distortion_kernel;

        uint8_t r, c;
        int32_t quadrant_size = context_ptr->blk_geom->sq_size >> 1;

        for (r = 0; r < 2; r++) {
            for (c = 0; c < 2; c++) {
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .rec_dist_per_quadrant[c + (r << 1)] = spatial_full_dist_type_fun(
                    input_picture_ptr->buffer_y,
                    input_origin_index + c * quadrant_size +
                        (r * quadrant_size) * input_picture_ptr->stride_y,
                    input_picture_ptr->stride_y,
                    recon_ptr->buffer_y,
                    blk_origin_index + c * quadrant_size +
                        (r * quadrant_size) * recon_ptr->stride_y,
                    recon_ptr->stride_y,
                    (uint32_t)quadrant_size,
                    (uint32_t)quadrant_size);
                // If quadrant_size == 4 then rec_dist_per_quadrant will have luma only because spatial_full_dist_type_fun does not support smaller than 4x4
                if (context_ptr->blk_geom->has_uv &&
                    context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1 && quadrant_size > 4) {
                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .rec_dist_per_quadrant[c + (r << 1)] += spatial_full_dist_type_fun(
                        input_picture_ptr->buffer_cb,
                        input_cb_origin_in_index + c * (quadrant_size >> 1) +
                            (r * (quadrant_size >> 1)) * input_picture_ptr->stride_cb,
                        input_picture_ptr->stride_cb,
                        recon_ptr->buffer_cb,
                        blk_chroma_origin_index + c * (quadrant_size >> 1) +
                            (r * (quadrant_size >> 1)) * recon_ptr->stride_cb,
                        recon_ptr->stride_cb,
                        (uint32_t)(quadrant_size >> 1),
                        (uint32_t)(quadrant_size >> 1));

                    context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                        .rec_dist_per_quadrant[c + (r << 1)] += spatial_full_dist_type_fun(
                        input_picture_ptr->buffer_cr,
                        input_cb_origin_in_index + c * (quadrant_size >> 1) +
                            (r * (quadrant_size >> 1)) * input_picture_ptr->stride_cr,
                        input_picture_ptr->stride_cr,
                        recon_ptr->buffer_cr,
                        blk_chroma_origin_index + c * (quadrant_size >> 1) +
                            (r * (quadrant_size >> 1)) * recon_ptr->stride_cr,
                        recon_ptr->stride_cr,
                        (uint32_t)(quadrant_size >> 1),
                        (uint32_t)(quadrant_size >> 1));
                }
            }
        }
    }
}
uint8_t is_intra_bordered(const ModeDecisionContext *ctx) {
    MacroBlockD *xd = ctx->blk_ptr->av1xd;

    const MbModeInfo *const above_mbmi = xd->above_mbmi;
    const MbModeInfo *const left_mbmi  = xd->left_mbmi;
    const int               has_above  = xd->up_available;
    const int               has_left   = xd->left_available;

    if (has_above && has_left)
        if ((!is_inter_block(&above_mbmi->block_mi)) && (!is_inter_block(&left_mbmi->block_mi)))
            return 1;
        else
            return 0;
    else
        return 0;
}
/*

*/
void md_encode_block_light_pd0(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                               EbPictureBufferDesc *input_picture_ptr) {
    const BlockGeom *blk_geom = context_ptr->blk_geom;
    uint32_t         fast_candidate_total_count;
    const uint32_t input_origin_index = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);
    const uint32_t blk_origin_index = blk_geom->origin_x +
        blk_geom->origin_y * context_ptr->sb_size;

    BlkStruct *blk_ptr = context_ptr->blk_ptr;
    if (!context_ptr->skip_intra) {
        init_xd(pcs_ptr, context_ptr);
        context_ptr->end_plane          = 1;
        context_ptr->uv_intra_comp_only = FALSE;
    }
    if (pcs_ptr->slice_type != I_SLICE) {
        context_ptr->me_sb_addr = context_ptr->sb_ptr->index;
        if (context_ptr->blk_geom->geom_idx == GEOM_0) {
            context_ptr->me_block_offset = me_idx_85[context_ptr->blk_geom->blkidx_mds];
            if (!pcs_ptr->parent_pcs_ptr->enable_me_8x8) {
                if (context_ptr->me_block_offset >= MAX_SB64_PU_COUNT_NO_8X8)
                    context_ptr->me_block_offset =
                        me_idx_85_8x8_to_16x16_conversion[context_ptr->me_block_offset -
                                                          MAX_SB64_PU_COUNT_NO_8X8];
                if (!pcs_ptr->parent_pcs_ptr->enable_me_16x16)
                    if (context_ptr->me_block_offset >= MAX_SB64_PU_COUNT_WO_16X16) {
                        assert(context_ptr->me_block_offset < 21);
                        context_ptr->me_block_offset =
                            me_idx_16x16_to_parent_32x32_conversion[context_ptr->me_block_offset -
                                                                    MAX_SB64_PU_COUNT_WO_16X16];
                    }
            }
        } else
            context_ptr->me_block_offset = me_idx[context_ptr->blk_geom->blkidx_mds];

        context_ptr->me_cand_offset = context_ptr->me_block_offset *
            pcs_ptr->parent_pcs_ptr->pa_me_data->max_cand;
    }

    generate_md_stage_0_cand_light_pd0(context_ptr, &fast_candidate_total_count, pcs_ptr);

    context_ptr->md_stage       = MD_STAGE_0;
    context_ptr->mds0_best_idx  = 0;
    context_ptr->mds0_best_cost = (uint64_t)~0;
    assert(fast_candidate_total_count <= context_ptr->max_nics && "not enough cand buffers");
    // If only one candidate, only need to perform compensation, not distortion calc
    // unless if VLPD0 where mds0 will become the last stage and SSD is needed
    if (fast_candidate_total_count == 1 && context_ptr->lpd0_ctrls.pd0_level != VERY_LIGHT_PD0) {
        ModeDecisionCandidateBuffer *candidate_buffer = context_ptr->candidate_buffer_ptr_array[0];
        candidate_buffer->candidate_ptr               = &context_ptr->fast_candidate_array[0];
        candidate_buffer->candidate_ptr->tx_depth     = 0;
        svt_product_prediction_fun_table_light_pd0[is_inter_mode(
            candidate_buffer->candidate_ptr->pred_mode)](
            context_ptr->hbd_mode_decision, context_ptr, pcs_ptr, candidate_buffer);
    } else
        md_stage_0_light_pd0(pcs_ptr,
                             context_ptr,
                             fast_candidate_total_count - 1,
                             input_picture_ptr,
                             input_origin_index,
                             blk_origin_index);

    if (context_ptr->lpd0_ctrls.pd0_level == VERY_LIGHT_PD0) {
        uint32_t rate = 0;
        uint64_t dist = context_ptr->mds0_best_cost;
        context_ptr->md_local_blk_unit[blk_ptr->mds_idx].cost =
            context_ptr->md_local_blk_unit[blk_ptr->mds_idx].default_cost = RDCOST(
                context_ptr->full_lambda_md[EB_8_BIT_MD], rate, dist);
    } else {
        context_ptr->md_stage = MD_STAGE_3;
        md_stage_3_light_pd0(
            pcs_ptr, context_ptr, input_picture_ptr, input_origin_index, blk_origin_index);
        // Update the cost
        context_ptr->md_local_blk_unit[blk_ptr->mds_idx].cost =
            context_ptr->md_local_blk_unit[blk_ptr->mds_idx].default_cost = *(
                context_ptr->candidate_buffer_ptr_array[context_ptr->mds0_best_idx]->full_cost_ptr);
    }
    assert(context_ptr->lpd1_ctrls.pd1_level < LPD1_LEVELS);
    // Save info used by depth refinemetn and the light-PD1 detector (detector uses 64x64 block info only)
    if (context_ptr->lpd0_ctrls.pd0_level != VERY_LIGHT_PD0) {
        // Save info needed only for LPD1 detector
        if (context_ptr->lpd1_ctrls.pd1_level > REGULAR_PD1 &&
            context_ptr->lpd1_ctrls.use_lpd1_detector[context_ptr->lpd1_ctrls.pd1_level] &&
            blk_geom->sq_size == 64) {
            ModeDecisionCandidate *candidate_ptr =
                context_ptr->candidate_buffer_ptr_array[context_ptr->mds0_best_idx]->candidate_ptr;
            blk_ptr->prediction_mode_flag = is_inter_mode(candidate_ptr->pred_mode) ? INTER_MODE
                                                                                    : INTRA_MODE;
            PredictionUnit *pu_ptr        = blk_ptr->prediction_unit_array;
            if (is_inter_mode(candidate_ptr->pred_mode)) {
                pu_ptr->inter_pred_direction_index = av1_get_pred_dir(
                    candidate_ptr->ref_frame_type);
                // Set MVs
                if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0) {
                    pu_ptr->mv[REF_LIST_0].as_int = candidate_ptr->mv[REF_LIST_0].as_int;
                } else if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1) {
                    pu_ptr->mv[REF_LIST_1].as_int = candidate_ptr->mv[REF_LIST_1].as_int;
                } else //if (pu_ptr->inter_pred_direction_index == BI_PRED)
                {
                    assert(pu_ptr->inter_pred_direction_index == BI_PRED);
                    pu_ptr->mv[REF_LIST_0].as_int = candidate_ptr->mv[REF_LIST_0].as_int;
                    pu_ptr->mv[REF_LIST_1].as_int = candidate_ptr->mv[REF_LIST_1].as_int;
                }
            }
        }

        // Save info needed for depth refinement
        context_ptr->md_local_blk_unit[blk_ptr->mds_idx].count_non_zero_coeffs =
            context_ptr->candidate_buffer_ptr_array[context_ptr->mds0_best_idx]
                ->count_non_zero_coeffs;
    }
    // If intra is used, generate recon and copy to necessary buffers
    if (!context_ptr->skip_intra) {
        uint32_t                     candidate_index = context_ptr->mds0_best_idx;
        ModeDecisionCandidateBuffer *candidate_buffer =
            context_ptr->candidate_buffer_ptr_array[candidate_index];

        // Update the variables needed for recon
        candidate_buffer->candidate_ptr->transform_type[0] = DCT_DCT;
        context_ptr->md_local_blk_unit[blk_ptr->mds_idx].y_has_coeff[0] =
            (uint8_t)candidate_buffer->y_has_coeff;
        // generate recon
        av1_perform_inverse_transform_recon(
            pcs_ptr, context_ptr, candidate_buffer, context_ptr->blk_geom);

        //copy neigh recon data in blk_ptr
        uint32_t             j;
        EbPictureBufferDesc *recon_ptr       = candidate_buffer->recon_ptr;
        uint32_t             rec_luma_offset = context_ptr->blk_geom->origin_x +
            context_ptr->blk_geom->origin_y * recon_ptr->stride_y;

        if (!context_ptr->hbd_mode_decision) {
            svt_memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon[0],
                       recon_ptr->buffer_y + rec_luma_offset +
                           (context_ptr->blk_geom->bheight - 1) * recon_ptr->stride_y,
                       context_ptr->blk_geom->bwidth);

            for (j = 0; j < context_ptr->blk_geom->bheight; ++j)
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon[0][j] =
                    recon_ptr->buffer_y[rec_luma_offset + context_ptr->blk_geom->bwidth - 1 +
                                        j * recon_ptr->stride_y];
        } else {
            uint16_t sz = sizeof(uint16_t);
            svt_memcpy(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                           .neigh_top_recon_16bit[0],
                       recon_ptr->buffer_y +
                           sz *
                               (rec_luma_offset +
                                (context_ptr->blk_geom->bheight - 1) * recon_ptr->stride_y),
                       sz * context_ptr->blk_geom->bwidth);

            for (j = 0; j < context_ptr->blk_geom->bheight; ++j)
                context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                    .neigh_left_recon_16bit[0][j] = ((uint16_t *)recon_ptr->buffer_y)
                    [rec_luma_offset + context_ptr->blk_geom->bwidth - 1 + j * recon_ptr->stride_y];
        }
    }

    context_ptr->avail_blk_flag[blk_ptr->mds_idx] = TRUE;
}

int get_comp_group_idx_context_enc(const MacroBlockD *xd);
int get_comp_index_context_enc(PictureParentControlSet *pcs_ptr, int cur_frame_index,
                               int bck_frame_index, int fwd_frame_index, const MacroBlockD *xd);
// Copy the recon samples to update the neighbour arrays
void copy_recon_md(PictureControlSet *pcs, ModeDecisionContext *ctx,
                   ModeDecisionCandidateBuffer *candidate_buffer) {
    const BlockGeom *blk_geom = ctx->blk_geom;
    if (!ctx->blk_geom->has_uv && ctx->cfl_ctrls.enabled) {
        // Store the luma data for 4x* and *x4 blocks to be used for CFL
        uint32_t dst_offset = blk_geom->origin_x +
            blk_geom->origin_y * candidate_buffer->recon_ptr->stride_y;
        uint32_t             dst_stride = candidate_buffer->recon_ptr->stride_y;
        EbPictureBufferDesc *recon_ptr  = candidate_buffer->recon_ptr;
        uint32_t rec_luma_offset = blk_geom->origin_x + blk_geom->origin_y * recon_ptr->stride_y;
        if (ctx->bypass_encdec && ctx->pd_pass == PD_PASS_1) {
            // If using only pred depth and no NSQ, can copy directly to final buffer b/c no d1 or d2 decision
            if (ctx->pred_depth_only && ctx->md_disallow_nsq) {
                get_recon_pic(pcs, &recon_ptr, ctx->hbd_mode_decision);
                rec_luma_offset = (recon_ptr->origin_y + ctx->blk_origin_y) * recon_ptr->stride_y +
                    (recon_ptr->origin_x + ctx->blk_origin_x);
            } else {
                recon_ptr       = ctx->blk_ptr->recon_tmp;
                rec_luma_offset = 0;
            }
        }
        // if using 8bit MD and bypassing encdec, need to save 8bit and 10bit recon
        if (ctx->encoder_bit_depth > EB_EIGHT_BIT && ctx->bypass_encdec &&
            !ctx->hbd_mode_decision && ctx->pd_pass == PD_PASS_1) {
            // copy 10bit
            if (ctx->pred_depth_only && ctx->md_disallow_nsq)
                get_recon_pic(pcs, &recon_ptr, 1);
            else
                recon_ptr = ctx->blk_ptr->recon_tmp;

            for (uint32_t j = 0; j < blk_geom->bheight; ++j)
                svt_memcpy(
                    ctx->cfl_temp_luma_recon16bit + dst_offset + j * dst_stride,
                    ((uint16_t *)recon_ptr->buffer_y) + (rec_luma_offset + j * recon_ptr->stride_y),
                    sizeof(uint16_t) * blk_geom->bwidth);

            // Copy 8bit
            if (ctx->pred_depth_only && ctx->md_disallow_nsq)
                recon_ptr = (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    ? ((EbReferenceObject *)
                           pcs->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                          ->reference_picture
                    : pcs->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr;
            else
                recon_ptr = ctx->blk_ptr->recon_tmp;

            for (uint32_t j = 0; j < blk_geom->bheight; ++j)
                svt_memcpy(&ctx->cfl_temp_luma_recon[dst_offset + j * dst_stride],
                           recon_ptr->buffer_y + rec_luma_offset + j * recon_ptr->stride_y,
                           blk_geom->bwidth);
        } else if (ctx->hbd_mode_decision) {
            for (uint32_t j = 0; j < blk_geom->bheight; ++j) {
                svt_memcpy(
                    ctx->cfl_temp_luma_recon16bit + dst_offset + j * dst_stride,
                    ((uint16_t *)recon_ptr->buffer_y) + (rec_luma_offset + j * recon_ptr->stride_y),
                    sizeof(uint16_t) * blk_geom->bwidth);
            }
        } else {
            for (uint32_t j = 0; j < blk_geom->bheight; ++j) {
                svt_memcpy(&ctx->cfl_temp_luma_recon[dst_offset + j * dst_stride],
                           recon_ptr->buffer_y + rec_luma_offset + j * recon_ptr->stride_y,
                           blk_geom->bwidth);
            }
        }
    } // END CFL COPIES
    //copy neigh recon data in blk_ptr
    EbPictureBufferDesc *recon_ptr = candidate_buffer->recon_ptr;
    uint32_t rec_luma_offset       = blk_geom->origin_x + blk_geom->origin_y * recon_ptr->stride_y;

    uint32_t rec_cb_offset = ((((blk_geom->origin_x >> 3) << 3) +
                               ((blk_geom->origin_y >> 3) << 3) *
                                   candidate_buffer->recon_ptr->stride_cb) >>
                              1);
    uint32_t rec_cr_offset = ((((blk_geom->origin_x >> 3) << 3) +
                               ((blk_geom->origin_y >> 3) << 3) *
                                   candidate_buffer->recon_ptr->stride_cr) >>
                              1);

    // If bypassing MD, recon is stored in different buffer; need to update the buffer to copy from
    if (ctx->bypass_encdec && ctx->pd_pass == PD_PASS_1) {
        // If using only pred depth and no NSQ, can copy directly to final buffer b/c no d1 or d2 decision
        if (ctx->pred_depth_only && ctx->md_disallow_nsq) {
            get_recon_pic(pcs, &recon_ptr, ctx->hbd_mode_decision);
            rec_luma_offset = (recon_ptr->origin_y + ctx->blk_origin_y) * recon_ptr->stride_y +
                (recon_ptr->origin_x + ctx->blk_origin_x);

            uint32_t round_origin_x = (ctx->blk_origin_x >> 3)
                << 3; // for Chroma blocks with size of 4
            uint32_t round_origin_y = (ctx->blk_origin_y >> 3)
                << 3; // for Chroma blocks with size of 4
            rec_cb_offset = rec_cr_offset = ((round_origin_x + recon_ptr->origin_x +
                                              (round_origin_y + recon_ptr->origin_y) *
                                                  recon_ptr->stride_cb) >>
                                             1);
        } else {
            recon_ptr       = ctx->blk_ptr->recon_tmp;
            rec_luma_offset = rec_cb_offset = rec_cr_offset = 0;
        }
    }
    // if using 8bit MD and bypassing encdec, need to save 8bit and 10bit recon
    if (ctx->encoder_bit_depth > EB_EIGHT_BIT && ctx->bypass_encdec && !ctx->hbd_mode_decision &&
        ctx->pd_pass == PD_PASS_1) {
        // copy 16bit recon
        if (ctx->pred_depth_only && ctx->md_disallow_nsq)
            get_recon_pic(pcs, &recon_ptr, 1);
        else
            recon_ptr = ctx->blk_ptr->recon_tmp;
        uint16_t sz = sizeof(uint16_t);
        // Copy bottom row (used for intra pred of the below block)
        svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[0],
                   recon_ptr->buffer_y +
                       sz * (rec_luma_offset + (blk_geom->bheight - 1) * recon_ptr->stride_y),
                   sz * blk_geom->bwidth);

        if (blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[1],
                       recon_ptr->buffer_cb +
                           sz * (rec_cb_offset + (blk_geom->bheight_uv - 1) * recon_ptr->stride_cb),
                       sz * blk_geom->bwidth_uv);

            svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[2],
                       recon_ptr->buffer_cr +
                           sz * (rec_cr_offset + (blk_geom->bheight_uv - 1) * recon_ptr->stride_cr),
                       sz * blk_geom->bwidth_uv);
        }

        // Copy right column (used for intra pred of the right block)
        for (uint32_t j = 0; j < blk_geom->bheight; ++j) {
            ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[0][j] =
                ((uint16_t *)recon_ptr
                     ->buffer_y)[rec_luma_offset + blk_geom->bwidth - 1 + j * recon_ptr->stride_y];
        }

        if (blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            for (uint32_t j = 0; j < blk_geom->bheight_uv; ++j) {
                ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[1][j] =
                    ((uint16_t *)recon_ptr->buffer_cb)[rec_cb_offset + blk_geom->bwidth_uv - 1 +
                                                       j * recon_ptr->stride_cb];
                ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[2][j] =
                    ((uint16_t *)recon_ptr->buffer_cr)[rec_cr_offset + blk_geom->bwidth_uv - 1 +
                                                       j * recon_ptr->stride_cr];
            }
        }

        // Copy 8bit recon
        if (ctx->pred_depth_only && ctx->md_disallow_nsq)
            recon_ptr = (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                ? ((EbReferenceObject *)
                       pcs->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                      ->reference_picture
                : pcs->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr;
        else
            recon_ptr = ctx->blk_ptr->recon_tmp;

        // Copy bottom row (used for intra pred of the below block)
        svt_memcpy(
            ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[0],
            recon_ptr->buffer_y + rec_luma_offset + (blk_geom->bheight - 1) * recon_ptr->stride_y,
            blk_geom->bwidth);

        if (blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[1],
                       recon_ptr->buffer_cb + rec_cb_offset +
                           (blk_geom->bheight_uv - 1) * recon_ptr->stride_cb,
                       blk_geom->bwidth_uv);
            svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[2],
                       recon_ptr->buffer_cr + rec_cr_offset +
                           (blk_geom->bheight_uv - 1) * recon_ptr->stride_cr,
                       blk_geom->bwidth_uv);
        }

        // Copy right column (used for intra pred of the right block)
        for (uint32_t j = 0; j < blk_geom->bheight; ++j)
            ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[0][j] =
                recon_ptr
                    ->buffer_y[rec_luma_offset + blk_geom->bwidth - 1 + j * recon_ptr->stride_y];

        if (blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            for (uint32_t j = 0; j < blk_geom->bheight_uv; ++j) {
                ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[1][j] =
                    recon_ptr->buffer_cb[rec_cb_offset + blk_geom->bwidth_uv - 1 +
                                         j * recon_ptr->stride_cb];
                ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[2][j] =
                    recon_ptr->buffer_cr[rec_cr_offset + blk_geom->bwidth_uv - 1 +
                                         j * recon_ptr->stride_cr];
            }
        }
    } else if (!ctx->hbd_mode_decision) {
        // Copy 8bit recon
        // Copy bottom row (used for intra pred of the below block)
        svt_memcpy(
            ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[0],
            recon_ptr->buffer_y + rec_luma_offset + (blk_geom->bheight - 1) * recon_ptr->stride_y,
            blk_geom->bwidth);

        if (blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[1],
                       recon_ptr->buffer_cb + rec_cb_offset +
                           (blk_geom->bheight_uv - 1) * recon_ptr->stride_cb,
                       blk_geom->bwidth_uv);
            svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[2],
                       recon_ptr->buffer_cr + rec_cr_offset +
                           (blk_geom->bheight_uv - 1) * recon_ptr->stride_cr,
                       blk_geom->bwidth_uv);
        }

        // Copy right column (used for intra pred of the right block)
        for (uint32_t j = 0; j < blk_geom->bheight; ++j) {
            ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[0][j] =
                recon_ptr
                    ->buffer_y[rec_luma_offset + blk_geom->bwidth - 1 + j * recon_ptr->stride_y];
        }

        if (blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            for (uint32_t j = 0; j < blk_geom->bheight_uv; ++j) {
                ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[1][j] =
                    recon_ptr->buffer_cb[rec_cb_offset + blk_geom->bwidth_uv - 1 +
                                         j * recon_ptr->stride_cb];
                ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[2][j] =
                    recon_ptr->buffer_cr[rec_cr_offset + blk_geom->bwidth_uv - 1 +
                                         j * recon_ptr->stride_cr];
            }
        }
    } else {
        // Copy 16bit recon
        uint16_t sz = sizeof(uint16_t);

        // Copy bottom row (used for intra pred of the below block)
        svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[0],
                   recon_ptr->buffer_y +
                       sz * (rec_luma_offset + (blk_geom->bheight - 1) * recon_ptr->stride_y),
                   sz * blk_geom->bwidth);

        if (blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[1],
                       recon_ptr->buffer_cb +
                           sz * (rec_cb_offset + (blk_geom->bheight_uv - 1) * recon_ptr->stride_cb),
                       sz * blk_geom->bwidth_uv);
            svt_memcpy(ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[2],
                       recon_ptr->buffer_cr +
                           sz * (rec_cr_offset + (blk_geom->bheight_uv - 1) * recon_ptr->stride_cr),
                       sz * blk_geom->bwidth_uv);
        }

        // Copy right column (used for intra pred of the right block)
        for (uint32_t j = 0; j < blk_geom->bheight; ++j) {
            ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[0][j] =
                ((uint16_t *)recon_ptr
                     ->buffer_y)[rec_luma_offset + blk_geom->bwidth - 1 + j * recon_ptr->stride_y];
        }

        if (blk_geom->has_uv && ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1) {
            for (uint32_t j = 0; j < blk_geom->bheight_uv; ++j) {
                ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[1][j] =
                    ((uint16_t *)recon_ptr->buffer_cb)[rec_cb_offset + blk_geom->bwidth_uv - 1 +
                                                       j * recon_ptr->stride_cb];
                ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[2][j] =
                    ((uint16_t *)recon_ptr->buffer_cr)[rec_cr_offset + blk_geom->bwidth_uv - 1 +
                                                       j * recon_ptr->stride_cr];
            }
        }
    } // END RECON COPIES
}
// Since light-PD1 uses pred_depth_only, the recon pixels can be copied directly to the recon buffer (no need
// to copy to a temp buffer and copy after d2 decision)
void copy_recon_light_pd1(PictureControlSet *pcs, ModeDecisionContext *ctx,
                          ModeDecisionCandidateBuffer *candidate_buffer) {
    const uint32_t blk_origin_x    = ctx->blk_origin_x;
    const uint32_t blk_origin_y    = ctx->blk_origin_y;
    const uint32_t bwidth          = ctx->blk_geom->bwidth;
    const uint32_t bheight         = ctx->blk_geom->bheight;
    const uint32_t blk_origin_x_uv = ctx->round_origin_x >> 1;
    const uint32_t blk_origin_y_uv = ctx->round_origin_y >> 1;
    const uint32_t bwidth_uv       = ctx->blk_geom->bwidth_uv;
    const uint32_t bheight_uv      = ctx->blk_geom->bheight_uv;

    //copy neigh recon data in blk_ptr
    uint32_t             j;
    EbPictureBufferDesc *recon_ptr;
    uint32_t             rec_luma_offset;
    uint32_t             rec_cb_offset;
    uint32_t             rec_cr_offset;

    // If bypassing MD, recon is stored in different buffer; need to update the buffer to copy from
    if (ctx->bypass_encdec) {
        // If using only pred depth and no NSQ, can copy directly to final buffer b/c no d1 or d2 decision
        // Assume non-16bit
        recon_ptr       = (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                  ? ((EbReferenceObject *)pcs->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                  ->reference_picture
                  : pcs->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr;
        rec_luma_offset = (recon_ptr->origin_y + blk_origin_y) * recon_ptr->stride_y +
            (recon_ptr->origin_x + blk_origin_x);
        rec_cb_offset = rec_cr_offset = ((blk_origin_x + recon_ptr->origin_x +
                                          (blk_origin_y + recon_ptr->origin_y) *
                                              recon_ptr->stride_cb) >>
                                         1);

    } else {
        recon_ptr       = candidate_buffer->recon_ptr;
        rec_luma_offset = ctx->blk_geom->origin_x + ctx->blk_geom->origin_y * recon_ptr->stride_y;
        rec_cb_offset   = ((ctx->blk_geom->origin_x +
                          ctx->blk_geom->origin_y * candidate_buffer->recon_ptr->stride_cb) >>
                         1);
        rec_cr_offset   = ((ctx->blk_geom->origin_x +
                          ctx->blk_geom->origin_y * candidate_buffer->recon_ptr->stride_cr) >>
                         1);
    }

    // Y
    // Copy top and bottom rows
    uint8_t *dst_ptr_top_left = ctx->luma_recon_neighbor_array->top_array +
        get_neighbor_array_unit_top_index(ctx->luma_recon_neighbor_array, blk_origin_x) *
            ctx->luma_recon_neighbor_array->unit_size;
    uint8_t *dst_ptr_bot_right = ctx->luma_recon_neighbor_array->top_left_array +
        get_neighbor_array_unit_top_left_index(
            ctx->luma_recon_neighbor_array, blk_origin_x, blk_origin_y + (bheight - 1)) *
            ctx->luma_recon_neighbor_array->unit_size;
    uint8_t *src_ptr = recon_ptr->buffer_y + rec_luma_offset + (bheight - 1) * recon_ptr->stride_y;
    svt_memcpy(dst_ptr_top_left, src_ptr, bwidth);
    svt_memcpy(dst_ptr_bot_right, src_ptr, bwidth);

    // Copy right and left columns
    dst_ptr_top_left = ctx->luma_recon_neighbor_array->left_array +
        get_neighbor_array_unit_left_index(ctx->luma_recon_neighbor_array, blk_origin_y) *
            ctx->luma_recon_neighbor_array->unit_size;
    dst_ptr_bot_right = ctx->luma_recon_neighbor_array->top_left_array +
        get_neighbor_array_unit_top_left_index(
            ctx->luma_recon_neighbor_array, blk_origin_x + (bwidth - 1), blk_origin_y) *
            ctx->luma_recon_neighbor_array->unit_size;
    src_ptr = recon_ptr->buffer_y + rec_luma_offset + bwidth - 1;
    for (j = 0; j < bheight; ++j) {
        *dst_ptr_bot_right = dst_ptr_top_left[j] = src_ptr[j * recon_ptr->stride_y];
        dst_ptr_bot_right -= 1;
    }

    // Cb
    // Copy top and bottom rows
    dst_ptr_top_left = ctx->cb_recon_neighbor_array->top_array +
        get_neighbor_array_unit_top_index(ctx->cb_recon_neighbor_array, blk_origin_x_uv) *
            ctx->cb_recon_neighbor_array->unit_size;
    dst_ptr_bot_right = ctx->cb_recon_neighbor_array->top_left_array +
        get_neighbor_array_unit_top_left_index(
            ctx->cb_recon_neighbor_array, blk_origin_x_uv, blk_origin_y_uv + (bheight_uv - 1)) *
            ctx->cb_recon_neighbor_array->unit_size;
    src_ptr = recon_ptr->buffer_cb + rec_cb_offset + (bheight_uv - 1) * recon_ptr->stride_cb;
    svt_memcpy(dst_ptr_top_left, src_ptr, bwidth_uv);
    svt_memcpy(dst_ptr_bot_right, src_ptr, bwidth_uv);

    // Copy right and left columns
    dst_ptr_top_left = ctx->cb_recon_neighbor_array->left_array +
        get_neighbor_array_unit_left_index(ctx->cb_recon_neighbor_array, blk_origin_y_uv) *
            ctx->cb_recon_neighbor_array->unit_size;
    dst_ptr_bot_right = ctx->cb_recon_neighbor_array->top_left_array +
        get_neighbor_array_unit_top_left_index(
            ctx->cb_recon_neighbor_array, blk_origin_x_uv + (bwidth_uv - 1), blk_origin_y_uv) *
            ctx->cb_recon_neighbor_array->unit_size;
    src_ptr = recon_ptr->buffer_cb + rec_cb_offset + bwidth_uv - 1;
    for (j = 0; j < bheight_uv; ++j) {
        *dst_ptr_bot_right = dst_ptr_top_left[j] = src_ptr[j * recon_ptr->stride_cb];
        dst_ptr_bot_right -= 1;
    }

    // Cr
    // Copy top and bottom rows
    dst_ptr_top_left = ctx->cr_recon_neighbor_array->top_array +
        get_neighbor_array_unit_top_index(ctx->cr_recon_neighbor_array, blk_origin_x_uv) *
            ctx->cr_recon_neighbor_array->unit_size;
    dst_ptr_bot_right = ctx->cr_recon_neighbor_array->top_left_array +
        get_neighbor_array_unit_top_left_index(
            ctx->cr_recon_neighbor_array, blk_origin_x_uv, blk_origin_y_uv + (bheight_uv - 1)) *
            ctx->cr_recon_neighbor_array->unit_size;
    src_ptr = recon_ptr->buffer_cr + rec_cr_offset + (bheight_uv - 1) * recon_ptr->stride_cr;
    svt_memcpy(dst_ptr_top_left, src_ptr, bwidth_uv);
    svt_memcpy(dst_ptr_bot_right, src_ptr, bwidth_uv);

    // Copy right and left columns
    dst_ptr_top_left = ctx->cr_recon_neighbor_array->left_array +
        get_neighbor_array_unit_left_index(ctx->cr_recon_neighbor_array, blk_origin_y_uv) *
            ctx->cr_recon_neighbor_array->unit_size;
    dst_ptr_bot_right = ctx->cr_recon_neighbor_array->top_left_array +
        get_neighbor_array_unit_top_left_index(
            ctx->cr_recon_neighbor_array, blk_origin_x_uv + (bwidth_uv - 1), blk_origin_y_uv) *
            ctx->cr_recon_neighbor_array->unit_size;
    src_ptr = recon_ptr->buffer_cr + rec_cr_offset + bwidth_uv - 1;
    for (j = 0; j < bheight_uv; ++j) {
        *dst_ptr_bot_right = dst_ptr_top_left[j] = src_ptr[j * recon_ptr->stride_cr];
        dst_ptr_bot_right -= 1;
    }

    // If bypassing EncDec for 10bit, need to save 8bit and 10bit recon
    if (ctx->encoder_bit_depth > EB_EIGHT_BIT && ctx->bypass_encdec) {
        get_recon_pic(pcs, &recon_ptr, 1);
        // Y
        // Copy top and bottom rows
        uint16_t *dst_ptr_top_left_16bit =
            (uint16_t *)(ctx->luma_recon_neighbor_array16bit->top_array +
                         get_neighbor_array_unit_top_index(ctx->luma_recon_neighbor_array16bit,
                                                           blk_origin_x) *
                             ctx->luma_recon_neighbor_array16bit->unit_size);
        uint16_t *dst_ptr_bot_right_16bit =
            (uint16_t *)(ctx->luma_recon_neighbor_array16bit->top_left_array +
                         get_neighbor_array_unit_top_left_index(ctx->luma_recon_neighbor_array16bit,
                                                                blk_origin_x,
                                                                blk_origin_y + (bheight - 1)) *
                             ctx->luma_recon_neighbor_array16bit->unit_size);
        uint16_t *src_ptr_16bit = ((uint16_t *)recon_ptr->buffer_y) + rec_luma_offset +
            (bheight - 1) * recon_ptr->stride_y;
        svt_memcpy(dst_ptr_top_left_16bit, src_ptr_16bit, bwidth * sizeof(uint16_t));
        svt_memcpy(dst_ptr_bot_right_16bit, src_ptr_16bit, bwidth * sizeof(uint16_t));

        // Copy right and left columns
        dst_ptr_top_left_16bit  = (uint16_t *)(ctx->luma_recon_neighbor_array16bit->left_array +
                                              get_neighbor_array_unit_left_index(
                                                  ctx->luma_recon_neighbor_array16bit,
                                                  blk_origin_y) *
                                                  ctx->luma_recon_neighbor_array16bit->unit_size);
        dst_ptr_bot_right_16bit = (uint16_t *)(ctx->luma_recon_neighbor_array16bit->top_left_array +
                                               get_neighbor_array_unit_top_left_index(
                                                   ctx->luma_recon_neighbor_array16bit,
                                                   blk_origin_x + (bwidth - 1),
                                                   blk_origin_y) *
                                                   ctx->luma_recon_neighbor_array16bit->unit_size);
        src_ptr_16bit           = ((uint16_t *)recon_ptr->buffer_y) + rec_luma_offset + bwidth - 1;
        for (j = 0; j < bheight; ++j) {
            *dst_ptr_bot_right_16bit = dst_ptr_top_left_16bit[j] =
                src_ptr_16bit[j * recon_ptr->stride_y];
            dst_ptr_bot_right_16bit -= 1;
        }

        // Cb
        // Copy top and bottom rows
        dst_ptr_top_left_16bit  = (uint16_t *)(ctx->cb_recon_neighbor_array16bit->top_array +
                                              get_neighbor_array_unit_top_index(
                                                  ctx->cb_recon_neighbor_array16bit,
                                                  blk_origin_x_uv) *
                                                  ctx->cb_recon_neighbor_array16bit->unit_size);
        dst_ptr_bot_right_16bit = (uint16_t *)(ctx->cb_recon_neighbor_array16bit->top_left_array +
                                               get_neighbor_array_unit_top_left_index(
                                                   ctx->cb_recon_neighbor_array16bit,
                                                   blk_origin_x_uv,
                                                   blk_origin_y_uv + (bheight_uv - 1)) *
                                                   ctx->cb_recon_neighbor_array16bit->unit_size);
        src_ptr_16bit           = ((uint16_t *)recon_ptr->buffer_cb) + rec_cb_offset +
            (bheight_uv - 1) * recon_ptr->stride_cb;
        svt_memcpy(dst_ptr_top_left_16bit, src_ptr_16bit, bwidth_uv * sizeof(uint16_t));
        svt_memcpy(dst_ptr_bot_right_16bit, src_ptr_16bit, bwidth_uv * sizeof(uint16_t));

        // Copy right and left columns
        dst_ptr_top_left_16bit  = (uint16_t *)(ctx->cb_recon_neighbor_array16bit->left_array +
                                              get_neighbor_array_unit_left_index(
                                                  ctx->cb_recon_neighbor_array16bit,
                                                  blk_origin_y_uv) *
                                                  ctx->cb_recon_neighbor_array16bit->unit_size);
        dst_ptr_bot_right_16bit = (uint16_t *)(ctx->cb_recon_neighbor_array16bit->top_left_array +
                                               get_neighbor_array_unit_top_left_index(
                                                   ctx->cb_recon_neighbor_array16bit,
                                                   blk_origin_x_uv + (bwidth_uv - 1),
                                                   blk_origin_y_uv) *
                                                   ctx->cb_recon_neighbor_array16bit->unit_size);
        src_ptr_16bit = ((uint16_t *)recon_ptr->buffer_cb) + rec_cb_offset + bwidth_uv - 1;
        for (j = 0; j < bheight_uv; ++j) {
            *dst_ptr_bot_right_16bit = dst_ptr_top_left_16bit[j] =
                src_ptr_16bit[j * recon_ptr->stride_cb];
            dst_ptr_bot_right_16bit -= 1;
        }

        // Cr
        // Copy top and bottom rows
        dst_ptr_top_left_16bit  = (uint16_t *)(ctx->cr_recon_neighbor_array16bit->top_array +
                                              get_neighbor_array_unit_top_index(
                                                  ctx->cr_recon_neighbor_array16bit,
                                                  blk_origin_x_uv) *
                                                  ctx->cr_recon_neighbor_array16bit->unit_size);
        dst_ptr_bot_right_16bit = (uint16_t *)(ctx->cr_recon_neighbor_array16bit->top_left_array +
                                               get_neighbor_array_unit_top_left_index(
                                                   ctx->cr_recon_neighbor_array16bit,
                                                   blk_origin_x_uv,
                                                   blk_origin_y_uv + (bheight_uv - 1)) *
                                                   ctx->cr_recon_neighbor_array16bit->unit_size);
        src_ptr_16bit           = ((uint16_t *)recon_ptr->buffer_cr) + rec_cr_offset +
            (bheight_uv - 1) * recon_ptr->stride_cr;
        svt_memcpy(dst_ptr_top_left_16bit, src_ptr_16bit, bwidth_uv * sizeof(uint16_t));
        svt_memcpy(dst_ptr_bot_right_16bit, src_ptr_16bit, bwidth_uv * sizeof(uint16_t));

        // Copy right and left columns
        dst_ptr_top_left_16bit  = (uint16_t *)(ctx->cr_recon_neighbor_array16bit->left_array +
                                              get_neighbor_array_unit_left_index(
                                                  ctx->cr_recon_neighbor_array16bit,
                                                  blk_origin_y_uv) *
                                                  ctx->cr_recon_neighbor_array16bit->unit_size);
        dst_ptr_bot_right_16bit = (uint16_t *)(ctx->cr_recon_neighbor_array16bit->top_left_array +
                                               get_neighbor_array_unit_top_left_index(
                                                   ctx->cr_recon_neighbor_array16bit,
                                                   blk_origin_x_uv + (bwidth_uv - 1),
                                                   blk_origin_y_uv) *
                                                   ctx->cr_recon_neighbor_array16bit->unit_size);
        src_ptr_16bit = ((uint16_t *)recon_ptr->buffer_cr) + rec_cr_offset + bwidth_uv - 1;
        for (j = 0; j < bheight_uv; ++j) {
            *dst_ptr_bot_right_16bit = dst_ptr_top_left_16bit[j] =
                src_ptr_16bit[j * recon_ptr->stride_cr];
            dst_ptr_bot_right_16bit -= 1;
        }
    }
}
/*
 * Convert the recon picture from 16bit to 8bit.  Recon pic is passed through the pcs.
 */
static void convert_recon_16bit_to_8bit(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    EbPictureBufferDesc *recon_buffer_16bit;
    EbPictureBufferDesc *recon_buffer_8bit;
    get_recon_pic(pcs, &recon_buffer_16bit, 1);
    get_recon_pic(pcs, &recon_buffer_8bit, 0);

    uint32_t pred_buf_x_offest = ctx->blk_origin_x;
    uint32_t pred_buf_y_offest = ctx->blk_origin_y;

    uint16_t *dst_16bit = (uint16_t *)(recon_buffer_16bit->buffer_y) + pred_buf_x_offest +
        recon_buffer_16bit->origin_x +
        (pred_buf_y_offest + recon_buffer_16bit->origin_y) * recon_buffer_16bit->stride_y;

    int32_t dst_stride_16bit = recon_buffer_16bit->stride_y;

    uint8_t *dst;
    int32_t  dst_stride;

    dst = recon_buffer_8bit->buffer_y + pred_buf_x_offest + recon_buffer_8bit->origin_x +
        (pred_buf_y_offest + recon_buffer_8bit->origin_y) * recon_buffer_8bit->stride_y;
    dst_stride = recon_buffer_8bit->stride_y;

    uint8_t *dst_nbit        = recon_buffer_8bit->buffer_bit_inc_y
               ? recon_buffer_8bit->buffer_bit_inc_y + pred_buf_x_offest + recon_buffer_8bit->origin_x +
            (pred_buf_y_offest + recon_buffer_8bit->origin_y) * recon_buffer_8bit->stride_y
               : recon_buffer_8bit->buffer_bit_inc_y;
    int32_t  dst_nbit_stride = recon_buffer_8bit->stride_bit_inc_y;

    un_pack2d(dst_16bit,
              dst_stride_16bit,
              dst,
              dst_stride,
              dst_nbit,
              dst_nbit_stride,
              ctx->blk_geom->bwidth,
              ctx->blk_geom->bheight);

    //copy recon from 16bit to 8bit
    pred_buf_x_offest = ((ctx->blk_origin_x >> 3) << 3) >> 1;
    pred_buf_y_offest = ((ctx->blk_origin_y >> 3) << 3) >> 1;

    dst_16bit = (uint16_t *)(recon_buffer_16bit->buffer_cb) + pred_buf_x_offest +
        recon_buffer_16bit->origin_x / 2 +
        (pred_buf_y_offest + recon_buffer_16bit->origin_y / 2) * recon_buffer_16bit->stride_cb;
    dst_stride_16bit = recon_buffer_16bit->stride_cb;

    dst = recon_buffer_8bit->buffer_cb + pred_buf_x_offest + recon_buffer_8bit->origin_x / 2 +
        (pred_buf_y_offest + recon_buffer_8bit->origin_y / 2) * recon_buffer_8bit->stride_cb;
    dst_stride = recon_buffer_8bit->stride_cb;

    dst_nbit        = recon_buffer_8bit->buffer_bit_inc_cb ? recon_buffer_8bit->buffer_bit_inc_cb +
            pred_buf_x_offest + recon_buffer_8bit->origin_x / 2 +
            (pred_buf_y_offest + recon_buffer_8bit->origin_y / 2) * recon_buffer_8bit->stride_cb
                                                           : recon_buffer_8bit->buffer_bit_inc_cb;
    dst_nbit_stride = recon_buffer_8bit->stride_bit_inc_cb;

    un_pack2d(dst_16bit,
              dst_stride_16bit,
              dst,
              dst_stride,
              dst_nbit,
              dst_nbit_stride,
              ctx->blk_geom->bwidth_uv,
              ctx->blk_geom->bheight_uv);

    dst_16bit = (uint16_t *)(recon_buffer_16bit->buffer_cr) +
        (pred_buf_x_offest + recon_buffer_16bit->origin_x / 2 +
         (pred_buf_y_offest + recon_buffer_16bit->origin_y / 2) * recon_buffer_16bit->stride_cr);
    dst_stride_16bit = recon_buffer_16bit->stride_cr;
    dst = recon_buffer_8bit->buffer_cr + pred_buf_x_offest + recon_buffer_8bit->origin_x / 2 +
        (pred_buf_y_offest + recon_buffer_8bit->origin_y / 2) * recon_buffer_8bit->stride_cr;
    dst_stride = recon_buffer_8bit->stride_cr;

    dst_nbit        = recon_buffer_8bit->buffer_bit_inc_cr ? recon_buffer_8bit->buffer_bit_inc_cr +
            pred_buf_x_offest + recon_buffer_8bit->origin_x / 2 +
            (pred_buf_y_offest + recon_buffer_8bit->origin_y / 2) * recon_buffer_8bit->stride_cr
                                                           : recon_buffer_8bit->buffer_bit_inc_cr;
    dst_nbit_stride = recon_buffer_8bit->stride_bit_inc_cr;

    un_pack2d(dst_16bit,
              dst_stride_16bit,
              dst,
              dst_stride,
              dst_nbit,
              dst_nbit_stride,
              ctx->blk_geom->bwidth_uv,
              ctx->blk_geom->bheight_uv);
}
// Check if reference frame pair of the current block matches with the given
// block.
static INLINE int match_ref_frame_pair(const BlockModeInfoEnc *mbmi,
                                       const MvReferenceFrame *ref_frames) {
    return ((ref_frames[0] == mbmi->ref_frame[0]) && (ref_frames[1] == mbmi->ref_frame[1]));
}
void lpd1_tx_shortcut_detector(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                               ModeDecisionCandidateBuffer **candidate_buffer_ptr_array) {
    const BlockGeom       *blk_geom = context_ptr->blk_geom;
    ModeDecisionCandidate *candidate_ptr =
        candidate_buffer_ptr_array[context_ptr->mds0_best_idx]->candidate_ptr;
    const uint32_t best_md_stage_dist =
        candidate_buffer_ptr_array[context_ptr->mds0_best_idx]->luma_fast_distortion;
    const uint32_t th_normalizer = blk_geom->bheight * blk_geom->bwidth *
        (pcs_ptr->picture_qp >> 1);
    context_ptr->use_tx_shortcuts_mds3 = (100 * best_md_stage_dist) <
        (context_ptr->lpd1_tx_ctrls.use_mds3_shortcuts_th * th_normalizer);
    context_ptr->lpd1_allow_skipping_tx = (100 * best_md_stage_dist) <
        (context_ptr->lpd1_tx_ctrls.skip_tx_th * th_normalizer);

    if (context_ptr->lpd1_tx_ctrls.use_neighbour_info && context_ptr->depth_removal_ctrls.enabled &&
        context_ptr->depth_removal_ctrls.disallow_below_64x64) {
        context_ptr->use_tx_shortcuts_mds3 = 1;
    }

    if ((!context_ptr->use_tx_shortcuts_mds3 || !context_ptr->lpd1_allow_skipping_tx) &&
        context_ptr->lpd1_tx_ctrls.use_neighbour_info && context_ptr->is_subres_safe) {
        MacroBlockD *xd = context_ptr->blk_ptr->av1xd;
        if (xd->left_available && xd->up_available) {
            const BlockModeInfoEnc *const left_mi  = &xd->left_mbmi->block_mi;
            const BlockModeInfoEnc *const above_mi = &xd->above_mbmi->block_mi;
            if (left_mi->skip && above_mi->skip) {
                MvReferenceFrame rf[2];
                av1_set_ref_frame(rf, candidate_ptr->ref_frame_type);
                int num_ref_frame_pair_match = match_ref_frame_pair(left_mi, rf);
                num_ref_frame_pair_match += match_ref_frame_pair(above_mi, rf);

                uint16_t use_tx_shortcuts_mds3_mult = 2 *
                    context_ptr->lpd1_tx_ctrls.use_neighbour_info; // is halved below
                uint16_t lpd1_allow_skipping_tx_mult = 2 *
                    context_ptr->lpd1_tx_ctrls.use_neighbour_info; // is halved below

                if (num_ref_frame_pair_match == 2) {
                    if (left_mi->mode == candidate_ptr->pred_mode &&
                        above_mi->mode == candidate_ptr->pred_mode) {
                        use_tx_shortcuts_mds3_mult  = 4;
                        lpd1_allow_skipping_tx_mult = 3;
                        if (is_inter_mode(candidate_ptr->pred_mode)) {
                            if (rf[1] != NONE_FRAME) { // BI_PRED
                                if (left_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_0].x &&
                                    left_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_0].y &&
                                    left_mi->mv[1].as_mv.col == candidate_ptr->mv[REF_LIST_1].x &&
                                    left_mi->mv[1].as_mv.row == candidate_ptr->mv[REF_LIST_1].y &&
                                    above_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_0].x &&
                                    above_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_0].y &&
                                    above_mi->mv[1].as_mv.col == candidate_ptr->mv[REF_LIST_1].x &&
                                    above_mi->mv[1].as_mv.row == candidate_ptr->mv[REF_LIST_1].y) {
                                    use_tx_shortcuts_mds3_mult  = 6;
                                    lpd1_allow_skipping_tx_mult = 4;
                                }
                            } else if (get_list_idx(rf[0]) == UNI_PRED_LIST_0) { // List 0 unipred
                                if (left_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_0].x &&
                                    left_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_0].y &&
                                    above_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_0].x &&
                                    above_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_0].y) {
                                    use_tx_shortcuts_mds3_mult  = 6;
                                    lpd1_allow_skipping_tx_mult = 4;
                                }
                            } else { // List 1 unipred
                                if (left_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_1].x &&
                                    left_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_1].y &&
                                    above_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_1].x &&
                                    above_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_1].y) {
                                    use_tx_shortcuts_mds3_mult  = 6;
                                    lpd1_allow_skipping_tx_mult = 4;
                                }
                            }
                        }
                    }
                }
                context_ptr->use_tx_shortcuts_mds3  = (100 * best_md_stage_dist) <
                        (((use_tx_shortcuts_mds3_mult *
                           context_ptr->lpd1_tx_ctrls.use_mds3_shortcuts_th) >>
                          1) *
                         th_normalizer)
                     ? 1
                     : context_ptr->use_tx_shortcuts_mds3;
                context_ptr->lpd1_allow_skipping_tx = (100 * best_md_stage_dist) <
                        (((lpd1_allow_skipping_tx_mult * context_ptr->lpd1_tx_ctrls.skip_tx_th) >>
                          1) *
                         th_normalizer)
                    ? 1
                    : context_ptr->lpd1_allow_skipping_tx;
            }
        }
    }
}
void md_encode_block_light_pd1(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                               uint32_t sb_addr, EbPictureBufferDesc *input_picture_ptr) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array;
    const BlockGeom              *blk_geom = context_ptr->blk_geom;
    ModeDecisionCandidateBuffer  *candidate_buffer;
    ModeDecisionCandidate        *fast_candidate_array = context_ptr->fast_candidate_array;
    uint32_t                      fast_candidate_total_count;

    BlockLocation loc;
    loc.input_origin_index = context_ptr->blk_origin_x + input_picture_ptr->origin_x +
        (context_ptr->blk_origin_y + input_picture_ptr->origin_y) * input_picture_ptr->stride_y;
    loc.input_cb_origin_in_index = ((context_ptr->blk_origin_x + input_picture_ptr->origin_x) >>
                                    1) +
        ((context_ptr->blk_origin_y + input_picture_ptr->origin_y) >> 1) *
            input_picture_ptr->stride_cb;
    loc.blk_origin_index        = blk_geom->origin_x + blk_geom->origin_y * context_ptr->sb_size;
    loc.blk_chroma_origin_index = (blk_geom->origin_x >> 1) +
        (blk_geom->origin_y >> 1) * (context_ptr->sb_size >> 1);

    BlkStruct *blk_ptr             = context_ptr->blk_ptr;
    candidate_buffer_ptr_array     = &(candidate_buffer_ptr_array_base[0]);
    context_ptr->blk_lambda_tuning = pcs_ptr->parent_pcs_ptr->blk_lambda_tuning;
    //Get the new lambda for current block
    if (pcs_ptr->parent_pcs_ptr->blk_lambda_tuning) {
        set_tuned_blk_lambda(context_ptr, pcs_ptr);
    }

    // need to init xd before product_coding_loop_init_fast_loop()
    init_xd(pcs_ptr, context_ptr);
    context_ptr->me_sb_addr = context_ptr->sb_ptr->index;
    if (context_ptr->blk_geom->geom_idx == GEOM_0) {
        context_ptr->me_block_offset = me_idx_85[context_ptr->blk_geom->blkidx_mds];
        if (!pcs_ptr->parent_pcs_ptr->enable_me_8x8) {
            if (context_ptr->me_block_offset >= MAX_SB64_PU_COUNT_NO_8X8)
                context_ptr->me_block_offset =
                    me_idx_85_8x8_to_16x16_conversion[context_ptr->me_block_offset -
                                                      MAX_SB64_PU_COUNT_NO_8X8];
            if (!pcs_ptr->parent_pcs_ptr->enable_me_16x16)
                if (context_ptr->me_block_offset >= MAX_SB64_PU_COUNT_WO_16X16) {
                    assert(context_ptr->me_block_offset < 21);
                    context_ptr->me_block_offset =
                        me_idx_16x16_to_parent_32x32_conversion[context_ptr->me_block_offset -
                                                                MAX_SB64_PU_COUNT_WO_16X16];
                }
        }
    } else
        context_ptr->me_block_offset = me_idx[context_ptr->blk_geom->blkidx_mds];
    if (!context_ptr->shut_fast_rate && pcs_ptr->slice_type != I_SLICE) {
        generate_av1_mvp_table(context_ptr,
                               context_ptr->blk_ptr,
                               context_ptr->blk_geom,
                               context_ptr->blk_origin_x,
                               context_ptr->blk_origin_y,
                               pcs_ptr->parent_pcs_ptr->ref_frame_type_arr,
                               pcs_ptr->parent_pcs_ptr->tot_ref_frame_types,
                               pcs_ptr);
    }

    product_coding_loop_init_fast_loop(pcs_ptr,
                                       context_ptr,
                                       context_ptr->skip_coeff_neighbor_array,
                                       context_ptr->leaf_partition_neighbor_array);
    context_ptr->is_intra_bordered =
        context_ptr->cand_reduction_ctrls.use_neighbouring_mode_ctrls.enabled
        ? is_intra_bordered(context_ptr)
        : 0;
    // Read and (if needed) perform 1/8 Pel ME MVs refinement
    if (pcs_ptr->slice_type != I_SLICE)
        read_refine_me_mvs_light_pd1(pcs_ptr, input_picture_ptr, context_ptr);
    generate_md_stage_0_cand_light_pd1(
        context_ptr->sb_ptr, context_ptr, &fast_candidate_total_count, pcs_ptr);

    if (pcs_ptr->slice_type != I_SLICE) {
        if (!context_ptr->shut_fast_rate) {
            estimate_ref_frames_num_bits(context_ptr, pcs_ptr);
        }
    }

    context_ptr->md_stage       = MD_STAGE_0;
    context_ptr->mds0_best_idx  = 0;
    context_ptr->mds0_best_cost = (uint64_t)~0;
    uint8_t perform_md_recon    = do_md_recon(pcs_ptr->parent_pcs_ptr, context_ptr);

    // If there is only a single candidate, skip compensation if transform will be skipped (unless compensation is needed for recon)
    if (fast_candidate_total_count > 1 || perform_md_recon ||
        context_ptr->lpd1_skip_inter_tx_level < 2 ||
        is_intra_mode(fast_candidate_array[0].pred_mode)) {
        md_stage_0_light_pd1(pcs_ptr,
                             context_ptr,
                             candidate_buffer_ptr_array_base,
                             fast_candidate_array,
                             fast_candidate_total_count - 1,
                             input_picture_ptr,
                             &loc,
                             blk_ptr);

        context_ptr->perform_mds1 = 0;
        lpd1_tx_shortcut_detector(pcs_ptr, context_ptr, candidate_buffer_ptr_array);
    } else {
        candidate_buffer                          = candidate_buffer_ptr_array_base[0];
        candidate_buffer->candidate_ptr           = &fast_candidate_array[0];
        *(candidate_buffer->fast_cost_ptr)        = 0;
        candidate_buffer->fast_luma_rate          = 0;
        candidate_buffer->fast_chroma_rate        = 0;
        candidate_buffer->candidate_ptr->tx_depth = 0;
        context_ptr->use_tx_shortcuts_mds3        = 1;
        context_ptr->lpd1_allow_skipping_tx       = 1;
    }
    // For 10bit content, when recon is not needed, hbd_mode_decision can stay =0,
    // and the 8bit prediction is used to produce the residual (with 8bit source).
    // When recon is needed, the prediction must be re-done in 10bit,
    // and the residual will be generated with the 10bit pred and 10bit source

    // Using the 8bit residual for the TX will cause different streams compared to using the 10bit residual.
    // To generate the same streams, compute the 10bit prediction before computing the recon

    if (context_ptr->encoder_bit_depth > EB_EIGHT_BIT && context_ptr->bypass_encdec &&
        perform_md_recon) {
        context_ptr->hbd_mode_decision = 2;

        // Update input pic and offsets
        input_picture_ptr      = pcs_ptr->input_frame16bit;
        loc.input_origin_index = context_ptr->blk_origin_x + input_picture_ptr->origin_x +
            (context_ptr->blk_origin_y + input_picture_ptr->origin_y) * input_picture_ptr->stride_y;
        loc.input_cb_origin_in_index = ((context_ptr->blk_origin_x + input_picture_ptr->origin_x) >>
                                        1) +
            ((context_ptr->blk_origin_y + input_picture_ptr->origin_y) >> 1) *
                input_picture_ptr->stride_cb;
        loc.blk_origin_index = blk_geom->origin_x + blk_geom->origin_y * context_ptr->sb_size;
        loc.blk_chroma_origin_index = (blk_geom->origin_x >> 1) +
            (blk_geom->origin_y >> 1) * (context_ptr->sb_size >> 1);
    }
    context_ptr->md_stage = MD_STAGE_3;
    md_stage_3_light_pd1(pcs_ptr, blk_ptr, context_ptr, input_picture_ptr, &loc);
    candidate_buffer = candidate_buffer_ptr_array[context_ptr->mds0_best_idx];

    // Full Mode Decision (choose the best mode)
    product_full_mode_decision_light_pd1(context_ptr, blk_ptr, pcs_ptr, sb_addr, candidate_buffer);
    // Perform inverse transform recon when needed
    if (perform_md_recon)
        av1_perform_inverse_transform_recon(
            pcs_ptr, context_ptr, candidate_buffer, context_ptr->blk_geom);

    // Convert 10bit recon (used as final EncDec recon) to 8bit recon (used for MD intra pred)
    if (context_ptr->encoder_bit_depth > EB_EIGHT_BIT && context_ptr->bypass_encdec &&
        context_ptr->hbd_mode_decision) {
        if (!context_ptr->skip_intra)
            convert_recon_16bit_to_8bit(pcs_ptr, context_ptr);
        context_ptr->hbd_mode_decision = 0;
    }
    if (!context_ptr->skip_intra) {
        copy_recon_light_pd1(pcs_ptr, context_ptr, candidate_buffer);
    }

    context_ptr->avail_blk_flag[blk_ptr->mds_idx] = TRUE;
}
void tx_shortcut_detector(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                          ModeDecisionCandidateBuffer **candidate_buffer_ptr_array) {
    const BlockGeom       *blk_geom = context_ptr->blk_geom;
    ModeDecisionCandidate *candidate_ptr =
        candidate_buffer_ptr_array[context_ptr->mds0_best_idx]->candidate_ptr;
    const uint32_t best_md_stage_dist =
        candidate_buffer_ptr_array[context_ptr->mds0_best_idx]->luma_fast_distortion;
    const uint32_t th_normalizer = blk_geom->bheight * blk_geom->bwidth *
        (pcs_ptr->picture_qp >> 1);
    context_ptr->use_tx_shortcuts_mds3 = (100 * best_md_stage_dist) <
        (context_ptr->tx_shortcut_ctrls.use_mds3_shortcuts_th * th_normalizer);

    if (!context_ptr->use_tx_shortcuts_mds3 && context_ptr->tx_shortcut_ctrls.use_neighbour_info &&
        context_ptr->is_subres_safe) {
        MacroBlockD *xd = context_ptr->blk_ptr->av1xd;
        if (xd->left_available && xd->up_available) {
            const BlockModeInfoEnc *const left_mi  = &xd->left_mbmi->block_mi;
            const BlockModeInfoEnc *const above_mi = &xd->above_mbmi->block_mi;
            if (left_mi->skip && above_mi->skip) {
                MvReferenceFrame rf[2];
                av1_set_ref_frame(rf, candidate_ptr->ref_frame_type);
                int num_ref_frame_pair_match = match_ref_frame_pair(left_mi, rf);
                num_ref_frame_pair_match += match_ref_frame_pair(above_mi, rf);

                uint16_t use_tx_shortcuts_mds3_mult = 2 *
                    context_ptr->tx_shortcut_ctrls.use_neighbour_info; // is halved below
                if (num_ref_frame_pair_match == 2) {
                    if (left_mi->mode == candidate_ptr->pred_mode &&
                        above_mi->mode == candidate_ptr->pred_mode) {
                        use_tx_shortcuts_mds3_mult = 4;
                        if (is_inter_mode(candidate_ptr->pred_mode)) {
                            if (rf[1] != NONE_FRAME) { // BI_PRED
                                if (left_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_0].x &&
                                    left_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_0].y &&
                                    left_mi->mv[1].as_mv.col == candidate_ptr->mv[REF_LIST_1].x &&
                                    left_mi->mv[1].as_mv.row == candidate_ptr->mv[REF_LIST_1].y &&
                                    above_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_0].x &&
                                    above_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_0].y &&
                                    above_mi->mv[1].as_mv.col == candidate_ptr->mv[REF_LIST_1].x &&
                                    above_mi->mv[1].as_mv.row == candidate_ptr->mv[REF_LIST_1].y) {
                                    use_tx_shortcuts_mds3_mult = 6;
                                }
                            } else if (get_list_idx(rf[0]) == UNI_PRED_LIST_0) { // List 0 unipred
                                if (left_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_0].x &&
                                    left_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_0].y &&
                                    above_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_0].x &&
                                    above_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_0].y) {
                                    use_tx_shortcuts_mds3_mult = 6;
                                }
                            } else { // List 1 unipred
                                if (left_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_1].x &&
                                    left_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_1].y &&
                                    above_mi->mv[0].as_mv.col == candidate_ptr->mv[REF_LIST_1].x &&
                                    above_mi->mv[0].as_mv.row == candidate_ptr->mv[REF_LIST_1].y) {
                                    use_tx_shortcuts_mds3_mult = 6;
                                }
                            }
                        }
                    }
                }
                context_ptr->use_tx_shortcuts_mds3 = (100 * best_md_stage_dist) <
                        (((use_tx_shortcuts_mds3_mult *
                           context_ptr->tx_shortcut_ctrls.use_mds3_shortcuts_th) >>
                          1) *
                         th_normalizer)
                    ? 1
                    : context_ptr->use_tx_shortcuts_mds3;
            }
        }
    }
}
void md_encode_block(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr, uint32_t sb_addr,
                     EbPictureBufferDesc *input_picture_ptr) {
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array_base =
        context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array;
    const BlockGeom              *blk_geom                 = context_ptr->blk_geom;
    uint32_t                      input_cb_origin_in_index = ((context_ptr->round_origin_y >> 1) +
                                         (input_picture_ptr->origin_y >> 1)) *
            input_picture_ptr->stride_cb +
        ((context_ptr->round_origin_x >> 1) + (input_picture_ptr->origin_x >> 1));
    uint32_t blk_chroma_origin_index = ROUND_UV(blk_geom->origin_x) / 2 +
        ROUND_UV(blk_geom->origin_y) / 2 * (context_ptr->sb_size >> 1);
    BlkStruct *blk_ptr             = context_ptr->blk_ptr;
    candidate_buffer_ptr_array     = &(candidate_buffer_ptr_array_base[0]);
    context_ptr->blk_lambda_tuning = pcs_ptr->parent_pcs_ptr->blk_lambda_tuning;
    //Get the new lambda for current block
    if (pcs_ptr->parent_pcs_ptr->blk_lambda_tuning) {
        set_tuned_blk_lambda(context_ptr, pcs_ptr);
    }

    init_xd(pcs_ptr, context_ptr);
    if (!context_ptr->shut_fast_rate) {
        FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
        // Generate MVP(s)
        if (frm_hdr->allow_intrabc) { // pcs_ptr->slice_type == I_SLICE
            MvReferenceFrame ref_frame = INTRA_FRAME;
            generate_av1_mvp_table(context_ptr,
                                   context_ptr->blk_ptr,
                                   context_ptr->blk_geom,
                                   context_ptr->blk_origin_x,
                                   context_ptr->blk_origin_y,
                                   &ref_frame,
                                   1,
                                   pcs_ptr);
        } else if (pcs_ptr->slice_type != I_SLICE) {
            generate_av1_mvp_table(context_ptr,
                                   context_ptr->blk_ptr,
                                   context_ptr->blk_geom,
                                   context_ptr->blk_origin_x,
                                   context_ptr->blk_origin_y,
                                   pcs_ptr->parent_pcs_ptr->ref_frame_type_arr,
                                   pcs_ptr->parent_pcs_ptr->tot_ref_frame_types,
                                   pcs_ptr);
        }
    }
    product_coding_loop_init_fast_loop(pcs_ptr,
                                       context_ptr,
                                       context_ptr->skip_coeff_neighbor_array,
                                       context_ptr->leaf_partition_neighbor_array);

    // Initialize uv_search_path
    if (context_ptr->uv_ctrls.nd_uv_serach_mode) {
        if (context_ptr->blk_geom->sq_size < 128) {
            if (context_ptr->blk_geom->has_uv) {
                init_chroma_mode(context_ptr);
            }
        }
    } else {
        // Search the best independent intra chroma mode
        if (context_ptr->uv_ctrls.uv_mode == CHROMA_MODE_0) {
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
    context_ptr->is_intra_bordered =
        context_ptr->cand_reduction_ctrls.use_neighbouring_mode_ctrls.enabled
        ? is_intra_bordered(context_ptr)
        : 0;
    if (context_ptr->md_pme_ctrls.modulate_pme_for_blk_size_res)
        context_ptr->updated_enable_pme = context_ptr->md_pme_ctrls.enabled &&
            context_ptr->blk_geom->sq_size >= 32;
    else
        context_ptr->updated_enable_pme = context_ptr->md_pme_ctrls.enabled;
    context_ptr->updated_enable_pme = context_ptr->is_intra_bordered &&
            context_ptr->cand_reduction_ctrls.use_neighbouring_mode_ctrls.enabled
        ? 0
        : context_ptr->updated_enable_pme;
    // Read MVPs (rounded-up to the closest integer) for use in md_sq_motion_search() and/or predictive_me_search() and/or perform_md_reference_pruning()
    if (pcs_ptr->slice_type != I_SLICE &&
        (context_ptr->md_sq_me_ctrls.enabled || context_ptr->updated_enable_pme ||
         context_ptr->ref_pruning_ctrls.enabled))
        build_single_ref_mvp_array(pcs_ptr, context_ptr);
    if (pcs_ptr->slice_type != I_SLICE)
        // Read and (if needed) perform 1/8 Pel ME MVs refinement
        read_refine_me_mvs(pcs_ptr, context_ptr, input_picture_ptr);
    // Initialized for eliminate_candidate_based_on_pme_me_results()
    context_ptr->pme_res[0][0].dist = context_ptr->pme_res[1][0].dist = 0xFFFFFFFF;
    // Perform md reference pruning
    if (context_ptr->ref_pruning_ctrls.enabled)
        perform_md_reference_pruning(pcs_ptr, context_ptr, input_picture_ptr);
    // Perform ME search around the best MVP
    if (context_ptr->updated_enable_pme) {
        pme_search(pcs_ptr, context_ptr, input_picture_ptr);
    }

    if (context_ptr->md_inter_intra_level) {
        int allow_ii = is_interintra_allowed_bsize(context_ptr->blk_geom->bsize);
        if (allow_ii)
            precompute_intra_pred_for_inter_intra(pcs_ptr, context_ptr);
    }

    uint32_t fast_candidate_total_count;
    generate_md_stage_0_cand(
        context_ptr->sb_ptr, context_ptr, &fast_candidate_total_count, pcs_ptr);
    if (pcs_ptr->slice_type != I_SLICE) {
        if (!context_ptr->shut_fast_rate) {
            estimate_ref_frames_num_bits(context_ptr, pcs_ptr);
        }
    }
    CandClass cand_class_it;
    uint32_t  buffer_start_idx = 0;
    uint32_t  buffer_count_for_curr_class;
    uint32_t  buffer_total_count        = 0;
    context_ptr->md_stage_1_total_count = 0;
    context_ptr->md_stage_2_total_count = 0;
    context_ptr->md_stage_3_total_count = 0;
    // Derive NIC(s)
    set_md_stage_counts(pcs_ptr, context_ptr);
    uint64_t best_md_stage_cost        = (uint64_t)~0;
    uint64_t best_md_stage_dist        = (uint64_t)~0;
    context_ptr->md_stage              = MD_STAGE_0;
    context_ptr->mds0_best_idx         = 0;
    context_ptr->mds0_best_class_it    = 0;
    context_ptr->mds1_best_idx         = 0;
    context_ptr->mds1_best_class_it    = 0;
    context_ptr->perform_mds1          = 1;
    context_ptr->use_tx_shortcuts_mds3 = 0;
    context_ptr->mds0_best_cost        = (uint64_t)~0;
    context_ptr->mds0_best_class       = 0;
    uint32_t input_origin_index        = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
            input_picture_ptr->stride_y +
        (context_ptr->blk_origin_x + input_picture_ptr->origin_x);
    uint32_t blk_origin_index = blk_geom->origin_x + blk_geom->origin_y * context_ptr->sb_size;
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
            assert_err(buffer_total_count <= context_ptr->max_nics, "not enough cand buffers");
            //Input: md_stage_0_count[cand_class_it]  Output:  md_stage_1_count[cand_class_it]
            context_ptr->target_class = cand_class_it;
            md_stage_0(
                pcs_ptr,
                context_ptr,
                candidate_buffer_ptr_array_base,
                context_ptr->fast_candidate_array,
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
            uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
            if (context_ptr->md_stage_1_count[cand_class_it] == 1) {
                cand_buff_indices[0] =
                    *(candidate_buffer_ptr_array[buffer_start_idx]->fast_cost_ptr) <
                        *(candidate_buffer_ptr_array[buffer_start_idx + 1]->fast_cost_ptr)
                    ? buffer_start_idx
                    : buffer_start_idx + 1;
            } else {
                sort_fast_cost_based_candidates(
                    context_ptr,
                    buffer_start_idx,
                    buffer_count_for_curr_class, //how many cand buffers to sort. one of the buffers can have max cost.
                    context_ptr->cand_buff_indices[cand_class_it]);
            }
            if (*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->fast_cost_ptr) <
                best_md_stage_cost) {
                best_md_stage_cost = *(
                    context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]->fast_cost_ptr);
                best_md_stage_dist = context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]
                                         ->luma_fast_distortion;
                context_ptr->mds0_best_idx      = cand_buff_indices[0];
                context_ptr->mds0_best_class_it = cand_class_it;
            }

            buffer_start_idx += buffer_count_for_curr_class; //for next iteration.
        }
    }
    post_mds0_nic_pruning(pcs_ptr, context_ptr, best_md_stage_cost, best_md_stage_dist);
    // Use detector for applying TX shortcuts at MDS3; if MDS1 is performed, use that info to apply
    // shortcuts instead of MDS0 info
    if (context_ptr->perform_mds1 == 0 && context_ptr->tx_shortcut_ctrls.use_mds3_shortcuts_th > 0)
        tx_shortcut_detector(pcs_ptr, context_ptr, candidate_buffer_ptr_array);
    // 1st Full-Loop
    best_md_stage_cost    = (uint64_t)~0;
    context_ptr->md_stage = MD_STAGE_1;
    for (cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++) {
        //number of next level candidates could not exceed number of curr level candidates
        context_ptr->md_stage_2_count[cand_class_it] = MIN(
            context_ptr->md_stage_1_count[cand_class_it],
            context_ptr->md_stage_2_count[cand_class_it]);
        if (context_ptr->perform_mds1) {
            if (context_ptr->bypass_md_stage_1 == FALSE &&
                context_ptr->md_stage_1_count[cand_class_it] > 0 &&
                context_ptr->md_stage_2_count[cand_class_it] > 0) {
                context_ptr->target_class = cand_class_it;
                md_stage_1(pcs_ptr,
                           blk_ptr,
                           context_ptr,
                           input_picture_ptr,
                           input_origin_index,
                           input_cb_origin_in_index,
                           blk_origin_index,
                           blk_chroma_origin_index);

                // Sort the candidates of the target class based on the 1st full loop cost

                //sort the new set of candidates
                if (context_ptr->md_stage_1_count[cand_class_it])
                    sort_full_cost_based_candidates(context_ptr,
                                                    context_ptr->md_stage_1_count[cand_class_it],
                                                    context_ptr->cand_buff_indices[cand_class_it]);
                uint32_t *cand_buff_indices = context_ptr->cand_buff_indices[cand_class_it];
                if (*(context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]
                          ->full_cost_ptr) < best_md_stage_cost) {
                    best_md_stage_cost = *(
                        context_ptr->candidate_buffer_ptr_array[cand_buff_indices[0]]
                            ->full_cost_ptr);
                    context_ptr->mds1_best_idx      = cand_buff_indices[0];
                    context_ptr->mds1_best_class_it = cand_class_it;
                }
            }
        } else {
            context_ptr->mds1_best_idx      = context_ptr->mds0_best_idx;
            context_ptr->mds1_best_class_it = context_ptr->mds0_best_class_it;
        }
    }
    if (context_ptr->perform_mds1)
        post_mds1_nic_pruning(context_ptr, best_md_stage_cost);
    // 2nd Full-Loop
    best_md_stage_cost    = (uint64_t)~0;
    context_ptr->md_stage = MD_STAGE_2;
    for (cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++) {
        //number of next level candidates could not exceed number of curr level candidates
        context_ptr->md_stage_3_count[cand_class_it] = MIN(
            context_ptr->md_stage_2_count[cand_class_it],
            context_ptr->md_stage_3_count[cand_class_it]);
        if (context_ptr->bypass_md_stage_2 == FALSE &&
            context_ptr->md_stage_2_count[cand_class_it] > 0 &&
            context_ptr->md_stage_3_count[cand_class_it] > 0) {
            context_ptr->target_class = cand_class_it;

            md_stage_2(pcs_ptr,
                       blk_ptr,
                       context_ptr,
                       input_picture_ptr,
                       input_origin_index,
                       input_cb_origin_in_index,
                       blk_origin_index,
                       blk_chroma_origin_index);

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

    if (context_ptr->perform_mds1) {
        post_mds2_nic_pruning(context_ptr, best_md_stage_cost);
        construct_best_sorted_arrays_md_stage_3(context_ptr,
                                                context_ptr->best_candidate_index_array);
    } else {
        context_ptr->md_stage_3_total_count = 1;
        context_ptr->best_candidate_index_array[0] =
            context_ptr->cand_buff_indices[context_ptr->mds1_best_class_it][0];
    }
    assert(context_ptr->md_stage_3_total_count > 0);
    // Search the best independent intra chroma mode
    if (context_ptr->uv_ctrls.nd_uv_serach_mode) {
        get_mds3_intra_count_for_chroma(
            context_ptr, candidate_buffer_ptr_array, context_ptr->best_candidate_index_array);
        // Initialize uv_search_path
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

    uint8_t org_hbd          = context_ptr->hbd_mode_decision;
    uint8_t perform_md_recon = do_md_recon(pcs_ptr->parent_pcs_ptr, context_ptr);

    // For 10bit content, when recon is not needed, hbd_mode_decision can stay =0,
    // and the 8bit prediction is used to produce the residual (with 8bit source).
    // When recon is needed, the prediction must be re-done in 10bit,
    // and the residual will be generated with the 10bit pred and 10bit source

    // Using the 8bit residual for the TX will cause different streams compared to using the 10bit residual.
    // To generate the same streams, compute the 10bit prediction before computing the recon

    if (context_ptr->encoder_bit_depth > EB_EIGHT_BIT && context_ptr->bypass_encdec &&
        !context_ptr->hbd_mode_decision && context_ptr->pd_pass == PD_PASS_1 && perform_md_recon) {
        context_ptr->hbd_mode_decision  = 2;
        context_ptr->need_hbd_comp_mds3 = 1;
        context_ptr->scale_palette      = 1;
        // Set the new input picture and offsets
        input_picture_ptr        = pcs_ptr->input_frame16bit;
        input_cb_origin_in_index = ((context_ptr->round_origin_y >> 1) +
                                    (input_picture_ptr->origin_y >> 1)) *
                input_picture_ptr->stride_cb +
            ((context_ptr->round_origin_x >> 1) + (input_picture_ptr->origin_x >> 1));
        blk_chroma_origin_index = ROUND_UV(blk_geom->origin_x) / 2 +
            ROUND_UV(blk_geom->origin_y) / 2 * (context_ptr->sb_size >> 1);
        input_origin_index = (context_ptr->blk_origin_y + input_picture_ptr->origin_y) *
                input_picture_ptr->stride_y +
            (context_ptr->blk_origin_x + input_picture_ptr->origin_x);
        blk_origin_index = blk_geom->origin_x + blk_geom->origin_y * context_ptr->sb_size;
    }
    // 3rd Full-Loop
    context_ptr->md_stage = MD_STAGE_3;
    md_stage_3(pcs_ptr,
               blk_ptr,
               context_ptr,
               input_picture_ptr,
               input_origin_index,
               input_cb_origin_in_index,
               blk_origin_index,
               blk_chroma_origin_index,
               context_ptr->md_stage_3_total_count);

    // Full Mode Decision (choose the best mode)
    uint32_t                     candidate_index  = product_full_mode_decision(context_ptr,
                                                          blk_ptr,
                                                          pcs_ptr,
                                                          sb_addr,
                                                          candidate_buffer_ptr_array,
                                                          context_ptr->md_stage_3_total_count,
                                                          context_ptr->best_candidate_index_array);
    ModeDecisionCandidateBuffer *candidate_buffer = candidate_buffer_ptr_array[candidate_index];
    //perform inverse transform recon when needed
    if (perform_md_recon)
        av1_perform_inverse_transform_recon(
            pcs_ptr, context_ptr, candidate_buffer, context_ptr->blk_geom);
    if ((!context_ptr->md_disallow_nsq && context_ptr->max_part0_to_part1_dev) ||
        (!context_ptr->disallow_4x4 && context_ptr->blk_geom->bsize == BLOCK_8X8))
        calc_scr_to_recon_dist_per_quadrant(context_ptr,
                                            input_picture_ptr,
                                            input_origin_index,
                                            input_cb_origin_in_index,
                                            candidate_buffer,
                                            blk_origin_index,
                                            blk_chroma_origin_index);
    if (context_ptr->encoder_bit_depth > EB_EIGHT_BIT && context_ptr->bypass_encdec && !org_hbd &&
        context_ptr->pd_pass == PD_PASS_1 && context_ptr->hbd_mode_decision) {
        if (!context_ptr->skip_intra)
            convert_recon_16bit_to_8bit(pcs_ptr, context_ptr);
        context_ptr->hbd_mode_decision  = 0;
        context_ptr->need_hbd_comp_mds3 = 0;
    }

    if (!context_ptr->skip_intra) {
        copy_recon_md(pcs_ptr, context_ptr, candidate_buffer);
    }
    context_ptr->avail_blk_flag[blk_ptr->mds_idx] = TRUE;
}
Bool update_skip_nsq_based_on_sq_recon_dist(ModeDecisionContext *ctx) {
    Bool             skip_nsq               = 0;
    uint32_t         max_part0_to_part1_dev = ctx->max_part0_to_part1_dev;
    const BlockGeom *blk_geom               = ctx->blk_geom;

    // return immediately if SQ, or NSQ but Parent not available, or max_part0_to_part1_dev is off
    if (blk_geom->shape == PART_N || ctx->avail_blk_flag[blk_geom->sqi_mds] == FALSE ||
        max_part0_to_part1_dev == 0)
        return skip_nsq;

    BlkStruct   *sq_blk_ptr           = &ctx->md_blk_arr_nsq[blk_geom->sqi_mds];
    MdBlkStruct *sq_md_local_blk_unit = &ctx->md_local_blk_unit[blk_geom->sqi_mds];

    // Derive the distortion/cost ratio
    const uint32_t full_lambda     = ctx->hbd_mode_decision ? ctx->full_lambda_md[EB_10_BIT_MD]
                                                            : ctx->full_lambda_md[EB_8_BIT_MD];
    const uint64_t dist            = RDCOST(full_lambda, 0, sq_md_local_blk_unit->full_distortion);
    const uint64_t dist_cost_ratio = (dist * 100) / sq_md_local_blk_unit->default_cost;
    const uint64_t min_ratio       = 0;
    const uint64_t max_ratio       = 100;
    const uint64_t modulated_th = (100 * (dist_cost_ratio - min_ratio)) / (max_ratio - min_ratio);

    // Modulate TH based on parent SQ pred_mode
    switch (sq_blk_ptr->pred_mode) {
    case NEWMV:
    case NEW_NEWMV:
        max_part0_to_part1_dev = max_part0_to_part1_dev - ((max_part0_to_part1_dev * 75) / 100);
        break;
    case NEAREST_NEARESTMV:
    case NEAR_NEARMV: max_part0_to_part1_dev *= 3; break;
    case GLOBALMV:
    case GLOBAL_GLOBALMV: max_part0_to_part1_dev <<= 2; break;
    default: break;
    }

    if (blk_geom->shape == PART_H || blk_geom->shape == PART_HA || blk_geom->shape == PART_HB ||
        blk_geom->shape == PART_H4) {
        // multiply the TH by 4 when Parent is D45 or D135 (diagonal) or when Parent is D67 / V / D113 (H_path)
        if (sq_blk_ptr->pred_mode == V_PRED || sq_blk_ptr->pred_mode == D67_PRED ||
            sq_blk_ptr->pred_mode == D113_PRED || sq_blk_ptr->pred_mode == D45_PRED ||
            sq_blk_ptr->pred_mode == D135_PRED)
            max_part0_to_part1_dev <<= 2;
        else if (sq_blk_ptr->pred_mode == H_PRED)
            max_part0_to_part1_dev = 0;

        const uint64_t dist_q0 = MAX(1, sq_md_local_blk_unit->rec_dist_per_quadrant[0]);
        const uint64_t dist_q1 = MAX(1, sq_md_local_blk_unit->rec_dist_per_quadrant[1]);
        const uint64_t dist_q2 = MAX(1, sq_md_local_blk_unit->rec_dist_per_quadrant[2]);
        const uint64_t dist_q3 = MAX(1, sq_md_local_blk_unit->rec_dist_per_quadrant[3]);

        const uint64_t dist_h0 = dist_q0 + dist_q1;
        const uint64_t dist_h1 = dist_q2 + dist_q3;

        const uint32_t dev = (uint32_t)((ABS((int64_t)dist_h0 - (int64_t)dist_h1) * 100) /
                                        MIN(dist_h0, dist_h1));
        // TH = TH + TH * Min(dev_0,dev_1); dev_0 is q0 - to - q1 deviation, and dev_1 is q2 - to - q3 deviation
        const uint32_t quad_dev_t = (uint32_t)((ABS((int64_t)dist_q0 - (int64_t)dist_q1) * 100) /
                                               MIN(dist_q0, dist_q1));
        const uint32_t quad_dev_b = (uint32_t)((ABS((int64_t)dist_q2 - (int64_t)dist_q3) * 100) /
                                               MIN(dist_q2, dist_q3));
        max_part0_to_part1_dev    = max_part0_to_part1_dev +
            (((uint64_t)max_part0_to_part1_dev * MIN(quad_dev_t, quad_dev_b)) / 100);

        max_part0_to_part1_dev = (uint32_t)((dist_cost_ratio <= min_ratio) ? 0
                                                : (dist_cost_ratio <= max_ratio)
                                                ? (max_part0_to_part1_dev * modulated_th) / 100
                                                : dist_cost_ratio);
        if (dev < max_part0_to_part1_dev)
            return TRUE;
    }

    if (blk_geom->shape == PART_V || blk_geom->shape == PART_VA || blk_geom->shape == PART_VB ||
        blk_geom->shape == PART_V4) {
        // multiply the TH by 4 when Parent is D45 or D135 (diagonal) or when Parent is D157 / H / D203 (V_path)
        if (sq_blk_ptr->pred_mode == H_PRED || sq_blk_ptr->pred_mode == D157_PRED ||
            sq_blk_ptr->pred_mode == D203_PRED || sq_blk_ptr->pred_mode == D45_PRED ||
            sq_blk_ptr->pred_mode == D135_PRED)
            max_part0_to_part1_dev <<= 2;
        else if (sq_blk_ptr->pred_mode == V_PRED)
            max_part0_to_part1_dev = 0;

        const uint64_t dist_q0 = MAX(1, sq_md_local_blk_unit->rec_dist_per_quadrant[0]);
        const uint64_t dist_q1 = MAX(1, sq_md_local_blk_unit->rec_dist_per_quadrant[1]);
        const uint64_t dist_q2 = MAX(1, sq_md_local_blk_unit->rec_dist_per_quadrant[2]);
        const uint64_t dist_q3 = MAX(1, sq_md_local_blk_unit->rec_dist_per_quadrant[3]);

        const uint64_t dist_v0 = dist_q0 + dist_q2;
        const uint64_t dist_v1 = dist_q1 + dist_q3;

        const uint32_t dev = (uint32_t)((ABS((int64_t)dist_v0 - (int64_t)dist_v1) * 100) /
                                        MIN(dist_v0, dist_v1));

        // TH = TH + TH * Min(dev_0,dev_1); dev_0 is q0-to-q2 deviation, and dev_1 is q1-to-q3 deviation
        const uint32_t quad_dev_l = (uint32_t)((ABS((int64_t)dist_q0 - (int64_t)dist_q2) * 100) /
                                               MIN(dist_q0, dist_q2));
        const uint32_t quad_dev_r = (uint32_t)((ABS((int64_t)dist_q1 - (int64_t)dist_q3) * 100) /
                                               MIN(dist_q1, dist_q3));
        max_part0_to_part1_dev    = max_part0_to_part1_dev +
            (((uint64_t)max_part0_to_part1_dev * MIN(quad_dev_l, quad_dev_r)) / 100);

        max_part0_to_part1_dev = (uint32_t)((dist_cost_ratio <= min_ratio) ? 0
                                                : (dist_cost_ratio <= max_ratio)
                                                ? (max_part0_to_part1_dev * modulated_th) / 100
                                                : dist_cost_ratio);
        if (dev < max_part0_to_part1_dev)
            return TRUE;
    }

    return skip_nsq;
}

/*
 * Determine if the evaluation of nsq blocks (HA, HB, VA, VB, H4, V4) can be skipped
 * based on the relative cost of the SQ, H, and V blocks.  The scaling factor sq_weight
 * determines how likely it is to skip blocks, and is a function of the qp, block shape,
 * prediction mode, block coeffs, and encode mode.  If skip_hv4_on_best_part is enabled
 * H4/V4 blocks will be skipped if the best partition so far is not H/V.
 *
 * skip HA, HB and H4 if (valid SQ and H) and (H_COST > (SQ_WEIGHT * SQ_COST) / 100)
 * skip VA, VB and V4 if (valid SQ and V) and (V_COST > (SQ_WEIGHT * SQ_COST) / 100)
 *
 * Returns TRUE if the blocks should be skipped; FALSE otherwise.
 */
static uint8_t update_skip_nsq_shapes(ModeDecisionContext *ctx) {
    const BlockGeom *blk_geom = ctx->blk_geom;
    const uint16_t   sqi      = blk_geom->sqi_mds;
    const Part       shape    = blk_geom->shape;

    if (ctx->skip_hv4_on_best_part) {
        // Skip H4/V4 shapes when best partition so far is not H/V
        if (shape == PART_H4) {
            if (ctx->md_blk_arr_nsq[sqi].part != PARTITION_HORZ)
                return 1;
        } else if (blk_geom->shape == PART_V4) {
            if (ctx->md_blk_arr_nsq[sqi].part != PARTITION_VERT)
                return 1;
        }
    }

    uint8_t  skip_nsq  = 0;
    uint32_t sq_weight = ctx->sq_weight;

    // return immediately if the skip nsq threshold is infinite
    if (sq_weight == (uint32_t)~0)
        return skip_nsq;

    // use a conservative threshold for H4, V4 blocks
    if (shape == PART_H4 || shape == PART_V4)
        sq_weight += CONSERVATIVE_OFFSET_0;

    MdBlkStruct *local_cu_unit = ctx->md_local_blk_unit;

    if (shape == PART_HA || shape == PART_HB || shape == PART_H4) {
        if (ctx->avail_blk_flag[sqi] && ctx->avail_blk_flag[sqi + 1] &&
            ctx->avail_blk_flag[sqi + 2]) {
            // Use aggressive thresholds for blocks without coeffs
            if (shape == PART_HA) {
                if (!ctx->md_blk_arr_nsq[sqi + 1].block_has_coeff)
                    sq_weight = (int32_t)sq_weight + AGGRESSIVE_OFFSET_1;
            }
            if (shape == PART_HB) {
                if (!ctx->md_blk_arr_nsq[sqi + 2].block_has_coeff)
                    sq_weight = (int32_t)sq_weight + AGGRESSIVE_OFFSET_1;
            }

            // compute the cost of the SQ block and H block
            const uint64_t sq_cost = local_cu_unit[sqi].default_cost;
            const uint64_t h_cost  = local_cu_unit[sqi + 1].default_cost +
                local_cu_unit[sqi + 2].default_cost;

            // Determine if nsq shapes can be skipped based on the relative cost of SQ and H blocks
            skip_nsq = (h_cost > ((sq_cost * sq_weight) / 100));
            // If not skipping, perform a check on the relative H/V costs
            if (!skip_nsq && ctx->avail_blk_flag[sqi + 3] && ctx->avail_blk_flag[sqi + 4]) {
                //compute the cost of V partition
                const uint64_t v_cost = local_cu_unit[sqi + 3].default_cost +
                    local_cu_unit[sqi + 4].default_cost;
                const uint32_t v_weight = 110;

                //if the cost of H partition is bigger than the V partition by a certain percentage
                skip_nsq = (h_cost > ((v_cost * v_weight) / 100));
            }
        }
    }

    if (shape == PART_VA || shape == PART_VB || shape == PART_V4) {
        if (ctx->avail_blk_flag[sqi] && ctx->avail_blk_flag[sqi + 3] &&
            ctx->avail_blk_flag[sqi + 4]) {
            // Use aggressive thresholds for blocks without coeffs
            if (shape == PART_VA) {
                if (!ctx->md_blk_arr_nsq[sqi + 3].block_has_coeff)
                    sq_weight = (int32_t)sq_weight + AGGRESSIVE_OFFSET_1;
            }
            if (shape == PART_VB) {
                if (!ctx->md_blk_arr_nsq[sqi + 4].block_has_coeff)
                    sq_weight = (int32_t)sq_weight + AGGRESSIVE_OFFSET_1;
            }

            // compute the cost of the SQ block and V block
            const uint64_t sq_cost = local_cu_unit[sqi].default_cost;
            const uint64_t v_cost  = local_cu_unit[sqi + 3].default_cost +
                local_cu_unit[sqi + 4].default_cost;

            // Determine if nsq shapes can be skipped based on the relative cost of SQ and V blocks
            skip_nsq = (v_cost > ((sq_cost * sq_weight) / 100));

            // If not skipping, perform a check on the relative H/V costs
            if (!skip_nsq && ctx->avail_blk_flag[sqi + 1] && ctx->avail_blk_flag[sqi + 2]) {
                const uint64_t h_cost = local_cu_unit[sqi + 1].default_cost +
                    local_cu_unit[sqi + 2].default_cost;
                const uint32_t h_weight = 110;

                //if the cost of V partition is bigger than the H partition by a certain percentage
                skip_nsq = (v_cost > ((h_cost * h_weight) / 100));
            }
        }
    }

    return skip_nsq;
}

void md_pme_search_controls(ModeDecisionContext *ctx, uint8_t md_pme_level);
void set_inter_intra_ctrls(ModeDecisionContext *ctx, uint8_t inter_intra_level);
void set_txt_controls(ModeDecisionContext *mdctxt, uint8_t txt_level);

// Set the levels used by features which apply aggressive settings for certain blocks (e.g. NSQ stats)
// Return 1 to skip, else 0
//
// Level 0 - skip
// Level 1-4 - adjust certain settings
Bool update_md_settings(ModeDecisionContext *ctx, uint8_t level) {
    // Level 0 is skip
    if (level == 0)
        return 1;

    if (level >= 1) {
        // Don't make NICs more conservative
        ctx->nic_ctrls.scaling_ctrls.stage1_scaling_num = MIN(
            ctx->nic_ctrls.scaling_ctrls.stage1_scaling_num, 5);
        ctx->nic_ctrls.scaling_ctrls.stage2_scaling_num = MIN(
            ctx->nic_ctrls.scaling_ctrls.stage2_scaling_num, 3);
        ctx->nic_ctrls.scaling_ctrls.stage3_scaling_num = MIN(
            ctx->nic_ctrls.scaling_ctrls.stage3_scaling_num, 3);
        ctx->txs_ctrls.enabled = 0;
    }
    if (level >= 2) {
        ctx->inter_comp_ctrls.tot_comp_types = 1;
        ctx->md_inter_intra_level            = 0;
        set_inter_intra_ctrls(ctx, ctx->md_inter_intra_level);
        md_pme_search_controls(ctx, 4);
    }
    if (level >= 3) {
        ctx->dist_based_ref_pruning = 8;
        svt_aom_set_dist_based_ref_pruning_controls(ctx, ctx->dist_based_ref_pruning);
        ctx->nic_ctrls.scaling_ctrls.stage1_scaling_num = MIN(
            ctx->nic_ctrls.scaling_ctrls.stage1_scaling_num, 2);
        ctx->nic_ctrls.scaling_ctrls.stage2_scaling_num = MIN(
            ctx->nic_ctrls.scaling_ctrls.stage2_scaling_num, 1);
        ctx->nic_ctrls.scaling_ctrls.stage3_scaling_num = MIN(
            ctx->nic_ctrls.scaling_ctrls.stage3_scaling_num, 1);
    }
    if (level >= 4) {
        set_txt_controls(ctx, 7);
        ctx->uv_ctrls.uv_mode = CHROMA_MODE_1;
    }
    return 0;
}

// Update MD settings or skip NSQ block based on the coeff-area of the parent SQ block
// Returns 1 to skip the NSQ block; 0 otherwise
uint8_t update_md_settings_based_on_sq_coeff_area(ModeDecisionContext *ctx) {
    uint8_t                                     skip_nsq = 0;
    ParentSqCoeffAreaBasedCyclesReductionCtrls *cycles_red_ctrls =
        &ctx->parent_sq_coeff_area_based_cycles_reduction_ctrls;
    if (cycles_red_ctrls->enabled) {
        if (ctx->blk_geom->shape != PART_N) {
            if (ctx->avail_blk_flag[ctx->blk_geom->sqi_mds]) {
                uint32_t count_non_zero_coeffs =
                    ctx->md_local_blk_unit[ctx->blk_geom->sqi_mds].count_non_zero_coeffs;
                uint32_t total_samples = (ctx->blk_geom->sq_size * ctx->blk_geom->sq_size);
                // High frequency band actions
                if (count_non_zero_coeffs >=
                    ((total_samples * cycles_red_ctrls->high_freq_band1_th) / 100))
                    skip_nsq = update_md_settings(ctx, cycles_red_ctrls->high_freq_band1_level);
                else if (count_non_zero_coeffs >=
                         ((total_samples * cycles_red_ctrls->high_freq_band2_th) / 100))
                    skip_nsq = update_md_settings(ctx, cycles_red_ctrls->high_freq_band2_level);
                else if (count_non_zero_coeffs >=
                         ((total_samples * cycles_red_ctrls->high_freq_band3_th) / 100))
                    skip_nsq = update_md_settings(ctx, cycles_red_ctrls->high_freq_band3_level);
                // Low frequency band actions
                else if (cycles_red_ctrls->enable_zero_coeff_action && count_non_zero_coeffs == 0) {
                    skip_nsq = update_md_settings(ctx, cycles_red_ctrls->zero_coeff_action);
                    set_txt_controls(ctx, 0);
                } else if (cycles_red_ctrls->enable_one_coeff_action && count_non_zero_coeffs == 1)
                    skip_nsq = update_md_settings(ctx, cycles_red_ctrls->one_coeff_action);
                else if (count_non_zero_coeffs <
                         ((total_samples * cycles_red_ctrls->low_freq_band1_th) / 100))
                    skip_nsq = update_md_settings(ctx, cycles_red_ctrls->low_freq_band1_level);
                else if (count_non_zero_coeffs <
                         ((total_samples * cycles_red_ctrls->low_freq_band2_th) / 100))
                    skip_nsq = update_md_settings(ctx, cycles_red_ctrls->low_freq_band2_level);
            }
        }
    }
    return skip_nsq;
}
/*
 * Pad high bit depth pictures.
 *
 * Returns pointer to padded data.
 */
EbPictureBufferDesc *pad_hbd_pictures(SequenceControlSet *scs, PictureControlSet *pcs,
                                      ModeDecisionContext *ctx, EbPictureBufferDesc *in_pic,
                                      uint16_t sb_org_x, uint16_t sb_org_y) {
    //perform the packing of 10bit if not done in previous PD passes
    if (!ctx->hbd_pack_done) {
        const uint32_t input_luma_offset = ((sb_org_y + in_pic->origin_y) * in_pic->stride_y) +
            (sb_org_x + in_pic->origin_x);
        const uint32_t input_cb_offset = (((sb_org_y + in_pic->origin_y) >> 1) *
                                          in_pic->stride_cb) +
            ((sb_org_x + in_pic->origin_x) >> 1);
        const uint32_t input_cr_offset = (((sb_org_y + in_pic->origin_y) >> 1) *
                                          in_pic->stride_cr) +
            ((sb_org_x + in_pic->origin_x) >> 1);

        uint32_t sb_width  = MIN(scs->sb_size, pcs->parent_pcs_ptr->aligned_width - sb_org_x);
        uint32_t sb_height = MIN(scs->sb_size, pcs->parent_pcs_ptr->aligned_height - sb_org_y);

        //sb_width is n*8 so the 2bit-decompression kernel works properly
        uint32_t comp_stride_y           = in_pic->stride_y / 4;
        uint32_t comp_luma_buffer_offset = comp_stride_y * in_pic->origin_y + in_pic->origin_x / 4;
        comp_luma_buffer_offset += sb_org_x / 4 + sb_org_y * comp_stride_y;

        compressed_pack_sb(in_pic->buffer_y + input_luma_offset,
                           in_pic->stride_y,
                           in_pic->buffer_bit_inc_y + comp_luma_buffer_offset,
                           comp_stride_y,
                           (uint16_t *)ctx->input_sample16bit_buffer->buffer_y,
                           ctx->input_sample16bit_buffer->stride_y,
                           sb_width,
                           sb_height);

        uint32_t comp_stride_uv            = in_pic->stride_cb / 4;
        uint32_t comp_chroma_buffer_offset = comp_stride_uv * (in_pic->origin_y / 2) +
            in_pic->origin_x / 2 / 4;
        comp_chroma_buffer_offset += sb_org_x / 4 / 2 + sb_org_y / 2 * comp_stride_uv;

        compressed_pack_sb(in_pic->buffer_cb + input_cb_offset,
                           in_pic->stride_cb,
                           in_pic->buffer_bit_inc_cb + comp_chroma_buffer_offset,
                           comp_stride_uv,
                           (uint16_t *)ctx->input_sample16bit_buffer->buffer_cb,
                           ctx->input_sample16bit_buffer->stride_cb,
                           sb_width / 2,
                           sb_height / 2);

        compressed_pack_sb(in_pic->buffer_cr + input_cr_offset,
                           in_pic->stride_cr,
                           in_pic->buffer_bit_inc_cr + comp_chroma_buffer_offset,
                           comp_stride_uv,
                           (uint16_t *)ctx->input_sample16bit_buffer->buffer_cr,
                           ctx->input_sample16bit_buffer->stride_cr,
                           sb_width / 2,
                           sb_height / 2);

        // PAD the packed source in incomplete sb up to max SB size
        pad_input_picture_16bit((uint16_t *)ctx->input_sample16bit_buffer->buffer_y,
                                ctx->input_sample16bit_buffer->stride_y,
                                sb_width,
                                sb_height,
                                scs->sb_size - sb_width,
                                scs->sb_size - sb_height);

        pad_input_picture_16bit((uint16_t *)ctx->input_sample16bit_buffer->buffer_cb,
                                ctx->input_sample16bit_buffer->stride_cb,
                                sb_width >> 1,
                                sb_height >> 1,
                                (scs->sb_size - sb_width) >> 1,
                                (scs->sb_size - sb_height) >> 1);

        pad_input_picture_16bit((uint16_t *)ctx->input_sample16bit_buffer->buffer_cr,
                                ctx->input_sample16bit_buffer->stride_cr,
                                sb_width >> 1,
                                sb_height >> 1,
                                (scs->sb_size - sb_width) >> 1,
                                (scs->sb_size - sb_height) >> 1);
        store16bit_input_src(
            ctx->input_sample16bit_buffer, pcs, sb_org_x, sb_org_y, scs->sb_size, scs->sb_size);

        ctx->hbd_pack_done = 1;
    }
    return ((scs->static_config.pass != ENC_FIRST_PASS) ? pcs->input_frame16bit : in_pic);
}

/*
 * Update the neighbour arrays before starting block processing.
 */
void update_neighbour_arrays_light_pd0(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    const uint16_t tile_idx = ctx->tile_index;

    // only need recon if using luma
    if (!ctx->skip_intra) {
        if (!ctx->hbd_mode_decision) {
            ctx->luma_recon_neighbor_array =
                pcs->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
            ctx->cb_recon_neighbor_array =
                pcs->md_cb_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
            ctx->cr_recon_neighbor_array =
                pcs->md_cr_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        } else {
            ctx->luma_recon_neighbor_array16bit =
                pcs->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
            ctx->cb_recon_neighbor_array16bit =
                pcs->md_cb_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
            ctx->cr_recon_neighbor_array16bit =
                pcs->md_cr_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        }
    }
}
/*
 * Update the neighbour arrays before starting block processing.
 */
void update_neighbour_arrays(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    const uint16_t tile_idx = ctx->tile_index;
    ctx->leaf_partition_neighbor_array =
        pcs->mdleaf_partition_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    if (ctx->encoder_bit_depth > EB_EIGHT_BIT && ctx->bypass_encdec && !ctx->hbd_mode_decision &&
        ctx->pd_pass == PD_PASS_1) {
        ctx->luma_recon_neighbor_array =
            pcs->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        ctx->cb_recon_neighbor_array =
            pcs->md_cb_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        ctx->cr_recon_neighbor_array =
            pcs->md_cr_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];

        ctx->luma_recon_neighbor_array16bit =
            pcs->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        ctx->cb_recon_neighbor_array16bit =
            pcs->md_cb_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        ctx->cr_recon_neighbor_array16bit =
            pcs->md_cr_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    } else if (!ctx->hbd_mode_decision) {
        ctx->luma_recon_neighbor_array =
            pcs->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        ctx->cb_recon_neighbor_array =
            pcs->md_cb_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        ctx->cr_recon_neighbor_array =
            pcs->md_cr_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    } else {
        ctx->luma_recon_neighbor_array16bit =
            pcs->md_luma_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        ctx->cb_recon_neighbor_array16bit =
            pcs->md_cb_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
        ctx->cr_recon_neighbor_array16bit =
            pcs->md_cr_recon_neighbor_array16bit[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    }
    ctx->skip_coeff_neighbor_array =
        pcs->md_skip_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    ctx->luma_dc_sign_level_coeff_neighbor_array =
        pcs->md_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    ctx->cb_dc_sign_level_coeff_neighbor_array =
        pcs->md_cb_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    ctx->cr_dc_sign_level_coeff_neighbor_array =
        pcs->md_cr_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
    ctx->txfm_context_array = pcs->md_txfm_context_array[MD_NEIGHBOR_ARRAY_INDEX][tile_idx];
}

/*
 * Initialize data needed for processing each block.  Update neighbour array if the block
 * is the first d1 block.  Called before process each block.
 */
void init_block_data(PictureControlSet *pcs, ModeDecisionContext *ctx,
                     const EbMdcLeafData *const leaf_data_ptr, const uint8_t blk_split_flag,
                     uint16_t sb_org_x, uint16_t sb_org_y, uint32_t blk_idx_mds,
                     uint8_t first_d1_blk) {
    const BlockGeom *blk_geom         = ctx->blk_geom;
    BlkStruct       *blk_ptr          = ctx->blk_ptr;
    ctx->scale_palette                = 0;
    ctx->blk_origin_x                 = sb_org_x + blk_geom->origin_x;
    ctx->blk_origin_y                 = sb_org_y + blk_geom->origin_y;
    ctx->round_origin_x               = ((ctx->blk_origin_x >> 3) << 3);
    ctx->round_origin_y               = ((ctx->blk_origin_y >> 3) << 3);
    ctx->tested_blk_flag[blk_idx_mds] = TRUE;
    blk_ptr->mds_idx                  = blk_idx_mds;
    blk_ptr->split_flag = blk_split_flag; //mdc indicates smallest or non valid CUs with split flag=
    blk_ptr->qindex     = ctx->qp_index;
    ctx->md_local_blk_unit[blk_idx_mds].left_neighbor_partition  = INVALID_NEIGHBOR_DATA;
    ctx->md_local_blk_unit[blk_idx_mds].above_neighbor_partition = INVALID_NEIGHBOR_DATA;
    ctx->sb64_sq_no4xn_geom                                      = 0;
    if (pcs->parent_pcs_ptr->scs_ptr->super_block_size == 64 &&
        blk_geom->bwidth == blk_geom->bheight && blk_geom->bsize > BLOCK_8X4)
        ctx->sb64_sq_no4xn_geom = 1;

    if (leaf_data_ptr->tot_d1_blocks != 1) {
        // We need to get the index of the sq_block for each NSQ branch
        if (first_d1_blk) {
            copy_neighbour_arrays( //save a clean neigh in [1], encode uses [0], reload the clean in [0] after done last ns block in a partition
                pcs,
                ctx,
                0,
                1,
                blk_geom->sqi_mds,
                sb_org_x,
                sb_org_y);
        }
    }
}
void check_curr_to_parent_cost_light_pd0(SequenceControlSet *scs, PictureControlSet *pcs,
                                         ModeDecisionContext *ctx, uint32_t sb_addr,
                                         uint32_t *next_non_skip_blk_idx_mds,
                                         Bool     *md_early_exit_sq) {
    const BlockGeom *blk_geom = ctx->blk_geom;
    BlkStruct       *blk_ptr  = ctx->blk_ptr;

    if (blk_geom->quadi > 0 && !(*md_early_exit_sq)) {
        uint64_t parent_depth_cost = 0, current_depth_cost = 0;

        // from a given child index, derive the index of the parent
        uint32_t parent_depth_idx_mds = blk_geom->parent_depth_idx_mds;

        if ((pcs->slice_type == I_SLICE && parent_depth_idx_mds == 0 &&
             scs->seq_header.sb_size == BLOCK_128X128) ||
            !pcs->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[parent_depth_idx_mds] ||
            ctx->tested_blk_flag[parent_depth_idx_mds] == 0) {
            *md_early_exit_sq = 0;
            return;
        } else {
            compute_depth_costs_md_skip_light_pd0(ctx,
                                                  parent_depth_idx_mds,
                                                  blk_geom->ns_depth_offset,
                                                  &parent_depth_cost,
                                                  &current_depth_cost);
            const uint32_t th = ctx->md_depth_early_exit_th == 0 ? 100
                                                                 : ctx->md_depth_early_exit_th;
            if (parent_depth_cost != MAX_MODE_COST &&
                (parent_depth_cost * th) <= (current_depth_cost * 100)) {
                *md_early_exit_sq          = 1;
                *next_non_skip_blk_idx_mds = parent_depth_idx_mds +
                    ns_depth_offset[blk_geom->geom_idx][blk_geom->depth - 1];
            } else {
                *md_early_exit_sq = 0;
            }
        }
    }
    // skip until we reach the next block @ the parent block depth
    if (blk_ptr->mds_idx >= *next_non_skip_blk_idx_mds && *md_early_exit_sq == 1)
        *md_early_exit_sq = 0;
}
/*
 * Check cost of current depth to parent depth, and if current cost is larger, signal to exit
 * processing depth early.
 */
void check_curr_to_parent_cost(SequenceControlSet *scs, PictureControlSet *pcs,
                               ModeDecisionContext *ctx, uint32_t sb_addr,
                               uint32_t *next_non_skip_blk_idx_mds, Bool *md_early_exit_sq,
                               uint8_t d1_blk_count) {
    const BlockGeom *blk_geom = ctx->blk_geom;
    BlkStruct       *blk_ptr  = ctx->blk_ptr;

    if (blk_geom->quadi > 0 && d1_blk_count == 0 && !(*md_early_exit_sq)) {
        uint64_t parent_depth_cost = 0, current_depth_cost = 0;

        // from a given child index, derive the index of the parent
        uint32_t parent_depth_idx_mds = blk_geom->parent_depth_idx_mds;
        assert(parent_depth_idx_mds == blk_geom->parent_depth_idx_mds);
        if ((pcs->slice_type == I_SLICE && parent_depth_idx_mds == 0 &&
             scs->seq_header.sb_size == BLOCK_128X128) ||
            !pcs->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[parent_depth_idx_mds])
            parent_depth_cost = MAX_MODE_COST;
        else
            compute_depth_costs_md_skip(ctx,
                                        pcs->parent_pcs_ptr,
                                        parent_depth_idx_mds,
                                        blk_geom->ns_depth_offset,
                                        &parent_depth_cost,
                                        &current_depth_cost);
        const uint32_t th = ctx->md_depth_early_exit_th == 0 ? 100 : ctx->md_depth_early_exit_th;
        if (parent_depth_cost != MAX_MODE_COST &&
            (parent_depth_cost * th) <= (current_depth_cost * 100)) {
            *md_early_exit_sq          = 1;
            *next_non_skip_blk_idx_mds = parent_depth_idx_mds +
                ns_depth_offset[blk_geom->geom_idx][blk_geom->depth - 1];
        } else {
            *md_early_exit_sq = 0;
        }
    }
    // skip until we reach the next block @ the parent block depth
    if (blk_ptr->mds_idx >= *next_non_skip_blk_idx_mds && *md_early_exit_sq == 1)
        *md_early_exit_sq = 0;
}

/*
 * Check if a block is redundant, and if so, copy the data from the original block
 * return 1 if block is redundant and updated, 0 otherwise
 */
Bool update_redundant(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    uint8_t          redundant_blk_avail = 0;
    uint16_t         redundant_blk_mds;
    const BlockGeom *blk_geom   = ctx->blk_geom;
    BlkStruct       *blk_ptr    = ctx->blk_ptr;
    uint8_t          bwidth     = blk_geom->bwidth;
    uint8_t          bheight    = blk_geom->bheight;
    uint8_t          bwidth_uv  = blk_geom->bwidth_uv;
    uint8_t          bheight_uv = blk_geom->bheight_uv;

    if (!ctx->md_disallow_nsq)
        check_redundant_block(blk_geom, ctx, &redundant_blk_avail, &redundant_blk_mds);

    ctx->similar_blk_avail = 0;

    if (!ctx->md_disallow_nsq)
        check_similar_block(blk_geom, ctx, &ctx->similar_blk_avail, &ctx->similar_blk_mds);

    if (redundant_blk_avail && ctx->redundant_blk) {
        // Copy results
        BlkStruct *src_cu = &ctx->md_blk_arr_nsq[redundant_blk_mds];
        BlkStruct *dst_cu = blk_ptr;
        move_blk_data_redund(pcs, ctx, src_cu, dst_cu);
        memcpy(&ctx->md_local_blk_unit[blk_ptr->mds_idx],
               &ctx->md_local_blk_unit[redundant_blk_mds],
               sizeof(MdBlkStruct));
        ctx->avail_blk_flag[dst_cu->mds_idx] = ctx->avail_blk_flag[redundant_blk_mds];

        if (!ctx->hbd_mode_decision) {
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[0],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_left_recon[0],
                   bheight);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[1],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_left_recon[1],
                   bheight_uv);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon[2],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_left_recon[2],
                   bheight_uv);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[0],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_top_recon[0],
                   bwidth);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[1],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_top_recon[1],
                   bwidth_uv);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon[2],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_top_recon[2],
                   bwidth_uv);
        } else {
            uint16_t sz = sizeof(uint16_t);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[0],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_left_recon_16bit[0],
                   bheight * sz);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[1],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_left_recon_16bit[1],
                   bheight_uv * sz);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_left_recon_16bit[2],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_left_recon_16bit[2],
                   bheight_uv * sz);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[0],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_top_recon_16bit[0],
                   bwidth * sz);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[1],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_top_recon_16bit[1],
                   bwidth_uv * sz);
            memcpy(&ctx->md_local_blk_unit[blk_geom->blkidx_mds].neigh_top_recon_16bit[2],
                   &ctx->md_local_blk_unit[redundant_blk_mds].neigh_top_recon_16bit[2],
                   bwidth_uv * sz);
        }
        return 1;
    }
    return 0;
}
/*

*/
void process_block_light_pd0(SequenceControlSet *scs, PictureControlSet *pcs,
                             ModeDecisionContext *ctx, const EbMdcLeafData *const leaf_data_ptr,
                             const uint8_t blk_split_flag, EbPictureBufferDesc *in_pic,
                             uint32_t sb_addr, uint16_t sb_org_x, uint16_t sb_org_y,
                             uint32_t blk_idx_mds, uint32_t *next_non_skip_blk_idx_mds,
                             Bool *md_early_exit_sq, uint8_t first_d1_blk) {
    const BlockGeom *blk_geom = ctx->blk_geom = get_blk_geom_mds(blk_idx_mds);
    BlkStruct       *blk_ptr = ctx->blk_ptr = &ctx->md_blk_arr_nsq[blk_idx_mds];
    init_block_data(
        pcs, ctx, leaf_data_ptr, blk_split_flag, sb_org_x, sb_org_y, blk_idx_mds, first_d1_blk);
    // Check current depth cost; if larger than parent, exit early
    check_curr_to_parent_cost_light_pd0(
        scs, pcs, ctx, sb_addr, next_non_skip_blk_idx_mds, md_early_exit_sq);

    // encode the current block only if it's not redundant
    {
        Bool skip_processing_block = *md_early_exit_sq;
        if (!skip_processing_block &&
            pcs->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[blk_ptr->mds_idx]) {
            // Encode the block
            md_encode_block_light_pd0(pcs, ctx, in_pic);
        } else if (!pcs->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[blk_ptr->mds_idx]) {
            // If the block is out of the boundaries, MD is not performed.
            // - For square blocks, since the blocks can be split further, they are considered in d2_inter_depth_block_decision() with cost of zero.
            // - For non-square blocks, since they can not be further split, the cost is set to the MAX value (MAX_MODE_COST) to ensure they are not selected.
            ctx->md_local_blk_unit[blk_ptr->mds_idx].cost         = (blk_geom->shape != PART_N)
                        ? MAX_MODE_COST >> 4
                        : 0;
            ctx->md_local_blk_unit[blk_ptr->mds_idx].default_cost = (blk_geom->shape != PART_N)
                ? MAX_MODE_COST >> 4
                : 0;
        } else {
            ctx->md_local_blk_unit[blk_ptr->mds_idx].cost         = MAX_MODE_COST >> 4;
            ctx->md_local_blk_unit[blk_ptr->mds_idx].default_cost = MAX_MODE_COST >> 4;
        }
    }
}

/*
 * Determined if a block should be processed, and if so, perform the MD pass on the block.
 */
void process_block_light_pd1(PictureControlSet *pcs, ModeDecisionContext *ctx,
                             const EbMdcLeafData *const leaf_data_ptr, EbPictureBufferDesc *in_pic,
                             uint32_t sb_addr, uint16_t sb_org_x, uint16_t sb_org_y,
                             uint32_t blk_idx_mds) {
    ctx->blk_geom = get_blk_geom_mds(blk_idx_mds);
    ctx->blk_ptr  = &ctx->md_blk_arr_nsq[blk_idx_mds];

    init_block_data(pcs,
                    ctx,
                    leaf_data_ptr,
                    FALSE, // blk_split_flag, - pred depth only; NSQ off
                    sb_org_x,
                    sb_org_y,
                    blk_idx_mds,
                    1); // first_d1_blk - NSQ off

    // Encode the block
    if (pcs->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[blk_idx_mds]) {
        md_encode_block_light_pd1(pcs, ctx, sb_addr, in_pic);
    }
}
/*
 * Determined if a block should be processed, and if so, perform the MD pass on the block.
 */
void process_block(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx,
                   const EbMdcLeafData *const leaf_data_ptr, const uint8_t blk_split_flag,
                   EbPictureBufferDesc *in_pic, uint32_t sb_addr, uint16_t sb_org_x,
                   uint16_t sb_org_y, uint32_t blk_idx_mds, uint32_t *next_non_skip_blk_idx_mds,
                   Bool *md_early_exit_sq, Bool *md_early_exit_nsq, uint8_t first_d1_blk,
                   uint8_t d1_blk_count) {
    const BlockGeom *blk_geom = ctx->blk_geom = get_blk_geom_mds(blk_idx_mds);
    BlkStruct       *blk_ptr = ctx->blk_ptr = &ctx->md_blk_arr_nsq[blk_idx_mds];
    init_block_data(
        pcs, ctx, leaf_data_ptr, blk_split_flag, sb_org_x, sb_org_y, blk_idx_mds, first_d1_blk);

    // Reset settings, in case they were over-written by previous block
    // Only reset settings when features that change settings are used.
    if (!ctx->md_disallow_nsq)
        signal_derivation_enc_dec_kernel_oq(scs, pcs, ctx);
    signal_derivation_block(pcs, ctx);
    // Check current depth cost; if larger than parent, exit early
    // if using pred depth only, you won't skip, so no need to check
    if (!(ctx->pd_pass == PD_PASS_1 && ctx->pred_depth_only))
        check_curr_to_parent_cost(
            scs, pcs, ctx, sb_addr, next_non_skip_blk_idx_mds, md_early_exit_sq, d1_blk_count);

    // encode the current block only if it's not redundant
    if (!ctx->redundant_blk || !update_redundant(pcs, ctx)) {
        Bool skip_processing_block = (*md_early_exit_nsq) || (*md_early_exit_sq) ||
            ctx->do_not_process_blk[blk_idx_mds];
        // call nsq-reduction func if NSQ is on
        if (!ctx->md_disallow_nsq) {
            skip_processing_block |= update_skip_nsq_based_on_sq_recon_dist(ctx);
            skip_processing_block |= update_skip_nsq_shapes(ctx);
            skip_processing_block |= update_md_settings_based_on_sq_coeff_area(ctx);
        }
        if (!skip_processing_block &&
            pcs->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[blk_ptr->mds_idx]) {
            // Encode the block
            md_encode_block(pcs, ctx, sb_addr, in_pic);
        } else if (!pcs->parent_pcs_ptr->sb_geom[sb_addr].block_is_allowed[blk_ptr->mds_idx]) {
            // If the block is out of the boundaries, MD is not performed.
            // - For square blocks, since the blocks can be split further, they are considered in d2_inter_depth_block_decision() with cost of zero.
            // - For non-square blocks, since they can not be further split, the cost is set to the MAX value (MAX_MODE_COST) to ensure they are not selected.
            ctx->md_local_blk_unit[blk_ptr->mds_idx].cost         = (blk_geom->shape != PART_N)
                        ? MAX_MODE_COST >> 4
                        : 0;
            ctx->md_local_blk_unit[blk_ptr->mds_idx].default_cost = (blk_geom->shape != PART_N)
                ? MAX_MODE_COST >> 4
                : 0;
        } else {
            ctx->md_local_blk_unit[blk_ptr->mds_idx].cost         = MAX_MODE_COST >> 4;
            ctx->md_local_blk_unit[blk_ptr->mds_idx].default_cost = MAX_MODE_COST >> 4;
        }
    }
}
/*
 * Update d1 data (including d1 decision) after each processed block, determine if should use early exit.
 */
void update_d1_data(PictureControlSet *pcs, ModeDecisionContext *ctx, uint16_t sb_org_x,
                    uint16_t sb_org_y, uint32_t blk_idx_mds, Bool *skip_next_nsq,
                    uint8_t *d1_blk_count) {
    const BlockGeom *blk_geom = ctx->blk_geom;
    BlkStruct       *blk_ptr  = ctx->blk_ptr;

    *skip_next_nsq = 0;
    if (blk_geom->nsi + 1 == blk_geom->totns) {
        d1_non_square_block_decision(ctx, *d1_blk_count);
        (*d1_blk_count)++;
    } else if (*d1_blk_count) {
        uint64_t tot_cost      = 0;
        uint32_t first_blk_idx = blk_ptr->mds_idx -
            (blk_geom->nsi); //index of first block in this partition
        for (int blk_it = 0; blk_it < blk_geom->nsi + 1; blk_it++)
            tot_cost += ctx->md_local_blk_unit[first_blk_idx + blk_it].cost;
        uint32_t full_lambda = ctx->hbd_mode_decision ? ctx->full_sb_lambda_md[EB_10_BIT_MD]
                                                      : ctx->full_sb_lambda_md[EB_8_BIT_MD];
        uint64_t part_cost   = svt_aom_partition_rate_cost(pcs->parent_pcs_ptr,
                                                         ctx,
                                                         &ctx->md_blk_arr_nsq[blk_geom->sqi_mds],
                                                         from_shape_to_part[blk_geom->shape],
                                                         full_lambda,
                                                         ctx->md_rate_estimation_ptr);

        tot_cost += part_cost;
        if (tot_cost > ctx->md_local_blk_unit[blk_geom->sqi_mds].cost)
            *skip_next_nsq = 1;
    }

    if (blk_geom->shape != PART_N) {
        if (blk_geom->nsi + 1 < blk_geom->totns)
            md_update_all_neighbour_arrays(pcs, ctx, blk_idx_mds, sb_org_x, sb_org_y);
        else
            copy_neighbour_arrays( //restore [1] in [0] after done last ns block
                pcs,
                ctx,
                1,
                0,
                blk_geom->sqi_mds,
                sb_org_x,
                sb_org_y);
    }
}

/*
 * Update d2 data (including d2 decision) after processing the last d1 block of a given square.
 */
void update_d2_decision_light_pd0(SequenceControlSet *scs, PictureControlSet *pcs,
                                  ModeDecisionContext *ctx, uint32_t sb_addr, uint16_t sb_org_x,
                                  uint16_t sb_org_y) {
    uint32_t last_blk_index_mds = d2_inter_depth_block_decision(
        scs,
        pcs,
        ctx,
        ctx->blk_geom->sqi_mds, //input is parent square
        sb_addr);

    // only needed to update recon
    if (!ctx->skip_intra && ctx->md_blk_arr_nsq[last_blk_index_mds].split_flag == FALSE) {
        ctx->blk_geom     = get_blk_geom_mds(last_blk_index_mds);
        ctx->blk_origin_x = sb_org_x + ctx->blk_geom->origin_x;
        ctx->blk_origin_y = sb_org_y + ctx->blk_geom->origin_y;
        ctx->blk_ptr      = &ctx->md_blk_arr_nsq[last_blk_index_mds];

        if (ctx->avail_blk_flag[last_blk_index_mds]) {
            mode_decision_update_neighbor_arrays_light_pd0(ctx);
        }
    }
}
/*
* Get the number of total block in a branch
*/
static uint32_t get_number_of_blocks(uint32_t block_idx) {
    const BlockGeom *blk_geom      = get_blk_geom_mds(block_idx);
    uint32_t         tot_d1_blocks = blk_geom->sq_size == 128 ? 17
                : blk_geom->sq_size > 8                       ? 25
                : blk_geom->sq_size == 8                      ? 5
                                                              : 1;
    return tot_d1_blocks;
}
/*
* Mark the blocks of the lower depth to be skipped
*/
static void set_child_to_be_skipped(ModeDecisionContext *context_ptr, uint32_t blk_index,
                                    int32_t sb_size, int8_t depth_step) {
    const BlockGeom *const blk_geom = get_blk_geom_mds(blk_index);

    if (context_ptr->md_blk_arr_nsq[blk_index].split_flag && blk_geom->sq_size > 4) {
        //Set first child to be considered
        uint32_t child_block_idx_1 = blk_index +
            d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
        uint32_t child1_tot_d1_blocks = get_number_of_blocks(child_block_idx_1);
        for (uint32_t block_1d_idx = 0; block_1d_idx < child1_tot_d1_blocks; block_1d_idx++)
            context_ptr->do_not_process_blk[child_block_idx_1 + block_1d_idx] = 1;
        if (depth_step > 1)
            set_child_to_be_skipped(context_ptr, child_block_idx_1, sb_size, depth_step - 1);
        //Set second child to be considered
        uint32_t child_block_idx_2 = child_block_idx_1 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        uint32_t child2_tot_d1_blocks = get_number_of_blocks(child_block_idx_2);
        for (uint32_t block_1d_idx = 0; block_1d_idx < child2_tot_d1_blocks; block_1d_idx++)
            context_ptr->do_not_process_blk[child_block_idx_2 + block_1d_idx] = 1;
        if (depth_step > 1)
            set_child_to_be_skipped(context_ptr, child_block_idx_2, sb_size, depth_step - 1);
        //Set third child to be considered
        uint32_t child_block_idx_3 = child_block_idx_2 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        uint32_t child3_tot_d1_blocks = get_number_of_blocks(child_block_idx_3);
        for (uint32_t block_1d_idx = 0; block_1d_idx < child3_tot_d1_blocks; block_1d_idx++)
            context_ptr->do_not_process_blk[child_block_idx_3 + block_1d_idx] = 1;
        if (depth_step > 1)
            set_child_to_be_skipped(context_ptr, child_block_idx_3, sb_size, depth_step - 1);
        //Set forth child to be considered
        uint32_t child_block_idx_4 = child_block_idx_3 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        uint32_t child4_tot_d1_blocks = get_number_of_blocks(child_block_idx_4);
        for (uint32_t block_1d_idx = 0; block_1d_idx < child4_tot_d1_blocks; block_1d_idx++)
            context_ptr->do_not_process_blk[child_block_idx_4 + block_1d_idx] = 1;
        if (depth_step > 1)
            set_child_to_be_skipped(context_ptr, child_block_idx_4, sb_size, depth_step - 1);
    }
}

static void eval_sub_depth_skip(SequenceControlSet *scs, ModeDecisionContext *ctx) {
    uint8_t n = 4;
    float   average, variance, std_deviation, sum = 0, sum1 = 0;

    // Compute the sum of all dist
    for (uint8_t q_idx = 0; q_idx < n; q_idx++) {
        sum = sum + ctx->md_local_blk_unit[ctx->blk_geom->sqi_mds].rec_dist_per_quadrant[q_idx];
    }
    average = sum / (float)n;

    // Compute variance and standard deviation
    for (uint8_t q_idx = 0; q_idx < n; q_idx++) {
        sum1 = sum1 +
            (float)pow(
                   (ctx->md_local_blk_unit[ctx->blk_geom->sqi_mds].rec_dist_per_quadrant[q_idx] -
                    average),
                   2);
    }
    variance      = sum1 / n;
    std_deviation = sqrtf(variance);

    // Derive the distortion/cost ratio
    uint32_t full_lambda = ctx->hbd_mode_decision ? ctx->full_lambda_md[EB_10_BIT_MD]
                                                  : ctx->full_lambda_md[EB_8_BIT_MD];
    uint64_t dist        = RDCOST(
        full_lambda, 0, ctx->md_local_blk_unit[ctx->blk_geom->sqi_mds].full_distortion);
    uint64_t dist_cost_ratio = (dist * 100) / ctx->md_local_blk_unit[ctx->blk_geom->sqi_mds].cost;
    float    min_ratio       = ctx->skip_4x4_depth_ctrls.min_distortion_cost_ratio;
    float    max_ratio       = 100;
    float    modulated_th    = (100 * (dist_cost_ratio - min_ratio)) / (max_ratio - min_ratio);
    float    quand_deviation_th = (dist_cost_ratio <= min_ratio) ? 0
           : (dist_cost_ratio <= max_ratio)
           ? (ctx->skip_4x4_depth_ctrls.quad_deviation_th * modulated_th) / 100
           : dist_cost_ratio;
    if (std_deviation < quand_deviation_th) {
        set_child_to_be_skipped(ctx,
                                ctx->blk_geom->sqi_mds,
                                scs->seq_header.sb_size,
                                ctx->skip_4x4_depth_ctrls.skip_all ? 6 : 1);
    }
}
/*
 * Update d2 data (including d2 decision) after processing the last d1 block of a given square.
 */
void update_d2_decision(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx,
                        uint32_t sb_addr, uint16_t sb_org_x, uint16_t sb_org_y) {
    uint32_t last_blk_index_mds;
    if (ctx->pd_pass == PD_PASS_1 && ctx->pred_depth_only)
        last_blk_index_mds = ctx->blk_geom->sqi_mds;
    else
        last_blk_index_mds = d2_inter_depth_block_decision(
            scs,
            pcs,
            ctx,
            ctx->blk_geom->sqi_mds, //input is parent square
            sb_addr);
    if (ctx->md_blk_arr_nsq[last_blk_index_mds].split_flag == FALSE) {
        md_update_all_neighbour_arrays_multiple(
            pcs, ctx, ctx->md_blk_arr_nsq[last_blk_index_mds].best_d1_blk, sb_org_x, sb_org_y);
    }
    // Here d1 is already performed but not d2
    if (ctx->skip_4x4_depth_ctrls.enabled && ctx->blk_geom->bsize == BLOCK_8X8 &&
        ctx->md_blk_arr_nsq[ctx->blk_geom->sqi_mds].split_flag && // could be further splitted
        ctx->avail_blk_flag[ctx->blk_geom->sqi_mds]) { // valid block

        eval_sub_depth_skip(scs, ctx);
    }
}

EB_EXTERN EbErrorType svt_aom_check_high_freq(PictureControlSet *pcs, SuperBlock *sb_ptr,
                                              ModeDecisionContext *ctx) {
    EbErrorType         return_error   = EB_ErrorNone;
    SequenceControlSet *scs_ptr        = pcs->scs_ptr;
    uint8_t             ref_frame_type = NONE_FRAME;

    ctx->sb_ptr = sb_ptr;

    if (pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] == 0 ||
        pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] <
            ctx->detect_high_freq_ctrls.me_8x8_sad_var_th)
        return return_error;

    EbPictureBufferDesc *input_pic = pcs->parent_pcs_ptr->enhanced_picture_ptr;

    MotionEstimationData *pa_me_data = pcs->parent_pcs_ptr->pa_me_data;

    uint32_t blk32_idx_tab[2][4] = {{1, 22, 43, 64}, {25, 294, 563, 832}};
    uint32_t sum_b32_satd        = 0;
    uint8_t  is_high_satd        = 0; // the b64 is detected if at least one b32 is detected

    for (uint32_t blk_idx = 0; blk_idx < 4; blk_idx++) {
        ctx->b32_satd[blk_idx] = (uint32_t)~0;

        uint32_t blk_idx_mds = blk32_idx_tab[scs_ptr->geom_idx][blk_idx];

        // block position should be calculated from the values in MD context,
        // because sb params are different since frames might be downscaled
        // if super-res or resize is enabled
        const BlockGeom *blk_geom = ctx->blk_geom = get_blk_geom_mds(blk_idx_mds);
        ctx->blk_origin_x                 = (uint16_t)(ctx->sb_origin_x + blk_geom->origin_x);
        ctx->blk_origin_y                 = (uint16_t)(ctx->sb_origin_y + blk_geom->origin_y);
        const uint32_t input_origin_index = (ctx->blk_origin_y + input_pic->origin_y) *
                input_pic->stride_y +
            (ctx->blk_origin_x + input_pic->origin_x);

        ctx->me_sb_addr = ctx->sb_ptr->index;

        ctx->me_block_offset = (ctx->blk_geom->geom_idx == GEOM_0)
            ? me_idx_85[ctx->blk_geom->blkidx_mds]
            : me_idx[ctx->blk_geom->blkidx_mds];

        ctx->me_cand_offset = ctx->me_block_offset * pa_me_data->max_cand;

        // ME offset(s)
        FrameHeader *frm_hdr             = &pcs->parent_pcs_ptr->frm_hdr;
        Bool         is_compound_enabled = (frm_hdr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;

        const MeSbResults *me_results = pa_me_data->me_results[ctx->me_sb_addr];
        uint8_t total_me_cnt          = me_results->total_me_candidate_index[ctx->me_block_offset];

        const MeCandidate *me_block_results = &me_results->me_candidate_array[ctx->me_cand_offset];

        const uint8_t max_refs = pa_me_data->max_refs;
        const uint8_t max_l0   = pa_me_data->max_l0;

        for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt;
             ++me_candidate_index) {
            const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
            const uint8_t      inter_direction      = me_block_results_ptr->direction;

            if (inter_direction == 2)
                break;

            const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
            const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;
            Mv            mv[MAX_NUM_OF_REF_PIC_LIST];
            /**************
                NEWMV L0
            ************* */
            mv[REF_LIST_0] = (Mv){{0, 0}};
            if (inter_direction == 0) {
                const int16_t mv_x =
                    (me_results->me_mv_array[ctx->me_block_offset * max_refs + list0_ref_index]
                         .x_mv)
                    << 1;
                const int16_t mv_y =
                    (me_results->me_mv_array[ctx->me_block_offset * max_refs + list0_ref_index]
                         .y_mv)
                    << 1;
                mv[REF_LIST_0] = (Mv){{mv_x, mv_y}};
                ref_frame_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
            }

            /**************
                NEWMV L1
            ************* */
            mv[REF_LIST_1] = (Mv){{0, 0}};
            if (is_compound_enabled) {
                if (inter_direction == 1) {
                    const int16_t mv_x = (me_results
                                              ->me_mv_array[ctx->me_block_offset * max_refs +
                                                            max_l0 + list0_ref_index]
                                              .x_mv)
                        << 1;
                    const int16_t mv_y = (me_results
                                              ->me_mv_array[ctx->me_block_offset * max_refs +
                                                            max_l0 + list0_ref_index]
                                              .y_mv)
                        << 1;
                    mv[REF_LIST_1] = (Mv){{mv_x, mv_y}};
                    ref_frame_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                }
            }

            MvReferenceFrame rf[2];
            av1_set_ref_frame(rf, ref_frame_type);

            assert(rf[1] == NONE_FRAME);

            EbPictureBufferDesc *ref_pic;
            const int8_t         ref_idx_first  = get_ref_frame_idx(rf[0]);
            const int8_t         list_idx_first = get_list_idx(rf[0]);
            int32_t              ref_origin_index;
            int16_t              mv_x, mv_y;

            EbReferenceObject *ref_obj = (EbReferenceObject *)pcs
                                             ->ref_pic_ptr_array[list_idx_first][ref_idx_first]
                                             ->object_ptr;

            if (list_idx_first == 0) {
                ref_pic = get_ref_pic_buffer(pcs, 0, 0, ref_idx_first);
                mv_x    = mv[REF_LIST_0].x >> 3;
                mv_y    = mv[REF_LIST_0].y >> 3;
            } else {
                ref_pic = get_ref_pic_buffer(pcs, 0, 1, ref_idx_first);
                mv_x    = mv[REF_LIST_1].x >> 3;
                mv_y    = mv[REF_LIST_1].y >> 3;
            }
            // -------
            // Use scaled references if resolution of the reference is different from that of the input
            // -------
            use_scaled_rec_refs_if_needed(pcs, input_pic, ref_obj, &ref_pic, 0);
            ref_origin_index = ref_pic->origin_x + (ctx->blk_origin_x + mv_x) +
                (ctx->blk_origin_y + mv_y + ref_pic->origin_y) * ref_pic->stride_y;

            uint32_t satd = hadamard_path(ctx,
                                          input_pic->buffer_y,
                                          input_origin_index,
                                          input_pic->stride_y,
                                          ref_pic->buffer_y,
                                          ref_origin_index,
                                          ref_pic->stride_y);

            ctx->b32_satd[blk_idx] = MIN(ctx->b32_satd[blk_idx], satd);
        }

        if (ctx->b32_satd[blk_idx] >= ctx->detect_high_freq_ctrls.high_satd_th) {
            is_high_satd = 1;
        }

        sum_b32_satd += ctx->b32_satd[blk_idx];

        if (is_high_satd &&
            sum_b32_satd > pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index]) {
            int dev = ((sum_b32_satd - pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index]) *
                       100) /
                pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index];
            if (dev >= ctx->detect_high_freq_ctrls.satd_to_sad_dev_th) {
                ctx->high_freq_present = 1;
                return return_error;
            }
        }
    }

    return return_error;
}

EB_EXTERN EbErrorType mode_decision_sb_light_pd0(SequenceControlSet *scs, PictureControlSet *pcs,
                                                 const MdcSbData *const mdc_sb_data,
                                                 SuperBlock *sb_ptr, uint16_t sb_org_x,
                                                 uint16_t sb_org_y, uint32_t sb_addr,
                                                 ModeDecisionContext *ctx) {
    EbErrorType return_error = EB_ErrorNone;
    ctx->sb_ptr              = sb_ptr;
    ctx->sb_origin_x         = sb_org_x;
    ctx->sb_origin_y         = sb_org_y;
    // Set SB-level variables here
    ctx->tx_depth              = 0;
    ctx->txb_1d_offset         = 0;
    ctx->txb_itr               = 0;
    ctx->luma_txb_skip_context = 0;
    ctx->luma_dc_sign_context  = 0;
    // Update neighbour arrays for the SB
    if (!ctx->skip_intra)
        update_neighbour_arrays_light_pd0(pcs, ctx);

    // get the input picture; if high bit-depth, pad the input pic
    EbPictureBufferDesc *input_pic = pcs->parent_pcs_ptr->enhanced_picture_ptr;
    if (ctx->hbd_mode_decision) {
        input_pic = pad_hbd_pictures(scs, pcs, ctx, input_pic, sb_org_x, sb_org_y);
    }

    // Initialize variables used to track blocks
    uint32_t                   leaf_count      = mdc_sb_data->leaf_count;
    const EbMdcLeafData *const leaf_data_array = mdc_sb_data->leaf_data_array;

    Bool     md_early_exit_sq          = 0;
    uint32_t next_non_skip_blk_idx_mds = 0;

    // Iterate over all blocks which are flagged to be considered
    for (uint32_t blk_idx = 0; blk_idx < leaf_count; blk_idx++) {
        uint32_t                   blk_idx_mds    = leaf_data_array[blk_idx].mds_idx;
        const EbMdcLeafData *const leaf_data_ptr  = &leaf_data_array[blk_idx];
        const uint8_t              blk_split_flag = mdc_sb_data->split_flag[blk_idx];

        process_block_light_pd0(scs,
                                pcs,
                                ctx,
                                leaf_data_ptr,
                                blk_split_flag,
                                input_pic,
                                sb_addr,
                                sb_org_x,
                                sb_org_y,
                                blk_idx_mds,
                                &next_non_skip_blk_idx_mds,
                                &md_early_exit_sq,
                                1); // first_d1_blk

        // Only using SQ, so always at tot_d1_blocks
        update_d2_decision_light_pd0(scs, pcs, ctx, sb_addr, sb_org_x, sb_org_y);
    }

    return return_error;
}
/*
 * Loop over all passed blocks in an SB and perform mode decision for each block,
 * then output the optimal mode distribution/partitioning for the given SB.
 *
 * For each block, selects the best mode through multiple MD stages (accuracy increases
 * while the number of mode candidates decreases as you move from one stage to another).
 * Based on the block costs, selects the best partition for a parent block (if NSQ
 * shapes are present). Finally, performs inter-depth decision towards a final partitiioning.
 */
EB_EXTERN void mode_decision_sb_light_pd1(SequenceControlSet *scs, PictureControlSet *pcs,
                                          const MdcSbData *const mdc_sb_data, SuperBlock *sb_ptr,
                                          uint16_t sb_org_x, uint16_t sb_org_y, uint32_t sb_addr,
                                          ModeDecisionContext *ctx) {
    ctx->sb_ptr      = sb_ptr;
    ctx->sb_origin_x = sb_org_x;
    ctx->sb_origin_y = sb_org_y;

    // Update neighbour arrays for the SB
    update_neighbour_arrays(pcs, ctx);

    // get the input picture; if high bit-depth, pad the input pic
    EbPictureBufferDesc *input_pic = pcs->parent_pcs_ptr->enhanced_picture_ptr;
    // If will need the 16bit picture, pad the input pic.  Done once for SB.
    if (ctx->encoder_bit_depth > EB_EIGHT_BIT && ctx->bypass_encdec) {
        // If using 8bit MD but bypassing EncDec, will need th 16bit pic later, but don't change input_pic
        pad_hbd_pictures(scs, pcs, ctx, input_pic, sb_org_x, sb_org_y);
    }

    // Initialize variables used to track blocks
    uint32_t                   leaf_count      = mdc_sb_data->leaf_count;
    const EbMdcLeafData *const leaf_data_array = mdc_sb_data->leaf_data_array;

    ctx->coded_area_sb         = 0;
    ctx->coded_area_sb_uv      = 0;
    ctx->tx_depth              = 0;
    ctx->txb_itr               = 0;
    ctx->txb_1d_offset         = 0;
    ctx->luma_txb_skip_context = 0;
    ctx->luma_dc_sign_context  = 0;
    ctx->cb_txb_skip_context   = 0;
    ctx->cb_dc_sign_context    = 0;
    ctx->cr_txb_skip_context   = 0;
    ctx->cr_dc_sign_context    = 0;

    // Iterate over all blocks which are flagged to be considered
    for (uint32_t blk_idx = 0; blk_idx < leaf_count; blk_idx++) {
        uint32_t                   blk_idx_mds   = leaf_data_array[blk_idx].mds_idx;
        const EbMdcLeafData *const leaf_data_ptr = &leaf_data_array[blk_idx];

        process_block_light_pd1(
            pcs, ctx, leaf_data_ptr, input_pic, sb_addr, sb_org_x, sb_org_y, blk_idx_mds);

        // If using pred_depth only and NSQ off, no need to update cost with skip flag
        ctx->md_blk_arr_nsq[blk_idx_mds].part        = PARTITION_NONE;
        ctx->md_blk_arr_nsq[blk_idx_mds].best_d1_blk = blk_idx_mds;

        // The current block is the last at a given d1 level; if so update d2 info
        if (ctx->avail_blk_flag[blk_idx_mds]) {
            // If TXS enabled at picture level, there are necessary context updates
            if (pcs->parent_pcs_ptr->frm_hdr.tx_mode == TX_MODE_SELECT) {
                uint8_t tx_size = tx_depth_to_tx_size[ctx->blk_ptr->tx_depth][ctx->blk_geom->bsize];
                uint8_t bw      = tx_size_wide[tx_size];
                uint8_t bh      = tx_size_high[tx_size];

                neighbor_array_unit_mode_write(ctx->txfm_context_array,
                                               &bw,
                                               ctx->blk_origin_x,
                                               ctx->blk_origin_y,
                                               ctx->blk_geom->bwidth,
                                               ctx->blk_geom->bheight,
                                               NEIGHBOR_ARRAY_UNIT_TOP_MASK);

                neighbor_array_unit_mode_write(ctx->txfm_context_array,
                                               &bh,
                                               ctx->blk_origin_x,
                                               ctx->blk_origin_y,
                                               ctx->blk_geom->bwidth,
                                               ctx->blk_geom->bheight,
                                               NEIGHBOR_ARRAY_UNIT_LEFT_MASK);
            }
            if (!ctx->shut_fast_rate ||
                ctx->cand_reduction_ctrls.use_neighbouring_mode_ctrls.enabled)
                update_mi_map(
                    ctx->blk_ptr, ctx->blk_origin_x, ctx->blk_origin_y, ctx->blk_geom, pcs);
        }
    }
}
/*
 * Loop over all passed blocks in an SB and perform mode decision for each block,
 * then output the optimal mode distribution/partitioning for the given SB.
 *
 * For each block, selects the best mode through multiple MD stages (accuracy increases
 * while the number of mode candidates decreases as you move from one stage to another).
 * Based on the block costs, selects the best partition for a parent block (if NSQ
 * shapes are present). Finally, performs inter-depth decision towards a final partitiioning.
 */
EB_EXTERN EbErrorType mode_decision_sb(SequenceControlSet *scs, PictureControlSet *pcs,
                                       const MdcSbData *const mdc_sb_data, SuperBlock *sb_ptr,
                                       uint16_t sb_org_x, uint16_t sb_org_y, uint32_t sb_addr,
                                       ModeDecisionContext *ctx) {
    EbErrorType return_error = EB_ErrorNone;
    ctx->sb_ptr              = sb_ptr;
    ctx->sb_origin_x         = sb_org_x;
    ctx->sb_origin_y         = sb_org_y;

    // Update neighbour arrays for the SB
    update_neighbour_arrays(pcs, ctx);

    // get the input picture; if high bit-depth, pad the input pic
    EbPictureBufferDesc *input_pic = pcs->parent_pcs_ptr->enhanced_picture_ptr;
    // If will need the 16bit picture, pad the input pic.  Done once for SB.
    if (ctx->hbd_mode_decision) {
        input_pic = pad_hbd_pictures(scs, pcs, ctx, input_pic, sb_org_x, sb_org_y);
    } else if (ctx->encoder_bit_depth > EB_EIGHT_BIT && ctx->bypass_encdec &&
               ctx->pd_pass == PD_PASS_1) {
        // If using 8bit MD but bypassing EncDec, will need th 16bit pic later, but don't change input_pic
        pad_hbd_pictures(scs, pcs, ctx, input_pic, sb_org_x, sb_org_y);
    }
    // Initialize variables used to track blocks
    uint32_t                   leaf_count      = mdc_sb_data->leaf_count;
    const EbMdcLeafData *const leaf_data_array = mdc_sb_data->leaf_data_array;

    Bool     md_early_exit_sq          = 0;
    Bool     md_early_exit_nsq         = 0;
    uint32_t next_non_skip_blk_idx_mds = 0;

    uint8_t  first_d1_blk         = 1;
    uint8_t  d1_blk_count         = 0;
    uint32_t d1_blocks_accumlated = 0;
    ctx->coded_area_sb            = 0;
    ctx->coded_area_sb_uv         = 0;
    // Iterate over all blocks which are flagged to be considered
    for (uint32_t blk_idx = 0; blk_idx < leaf_count; blk_idx++) {
        uint32_t                   blk_idx_mds    = leaf_data_array[blk_idx].mds_idx;
        const EbMdcLeafData *const leaf_data_ptr  = &leaf_data_array[blk_idx];
        const uint8_t              blk_split_flag = mdc_sb_data->split_flag[blk_idx];
        process_block(scs,
                      pcs,
                      ctx,
                      leaf_data_ptr,
                      blk_split_flag,
                      input_pic,
                      sb_addr,
                      sb_org_x,
                      sb_org_y,
                      blk_idx_mds,
                      &next_non_skip_blk_idx_mds,
                      &md_early_exit_sq,
                      &md_early_exit_nsq,
                      first_d1_blk,
                      d1_blk_count);

        // If using pred_depth only and NSQ off, no need to update cost with skip flag
        if (ctx->pd_pass == PD_PASS_1 && ctx->md_disallow_nsq && ctx->pred_depth_only) {
            ctx->md_blk_arr_nsq[blk_idx_mds].part        = PARTITION_NONE;
            ctx->md_blk_arr_nsq[blk_idx_mds].best_d1_blk = blk_idx_mds;
        } else
            update_d1_data(
                pcs, ctx, sb_org_x, sb_org_y, blk_idx_mds, &md_early_exit_nsq, &d1_blk_count);

        // Check if the current block is the last at a given d1 level; if so update d2 info
        d1_blocks_accumlated = (first_d1_blk == 1) ? 1 : d1_blocks_accumlated + 1;
        if (d1_blocks_accumlated == leaf_data_ptr->tot_d1_blocks) {
            // Perform d2 inter-depth decision after final d1 block
            update_d2_decision(scs, pcs, ctx, sb_addr, sb_org_x, sb_org_y);

            first_d1_blk = 1;
            d1_blk_count = 0;
        } else if (first_d1_blk) {
            first_d1_blk = 0;
        }
    }

    return return_error;
}
#define MAX_SEARCH_POINT_WIDTH 128
#define MAX_SEARCH_POINT_HEIGHT 128

#define MAX_TATAL_SEARCH_AREA_WIDTH (MAX_SB_SIZE + MAX_SEARCH_POINT_WIDTH + ME_FILTER_TAP)
#define MAX_TATAL_SEARCH_AREA_HEIGHT (MAX_SB_SIZE + MAX_SEARCH_POINT_HEIGHT + ME_FILTER_TAP)

#define MAX_SEARCH_AREA_SIZE MAX_TATAL_SEARCH_AREA_WIDTH *MAX_TATAL_SEARCH_AREA_HEIGHT
