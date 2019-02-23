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
#include "EbComputeSAD.h"
#include "EbTransforms.h"
#include "EbMeSadCalculation.h"
#include "EbMotionEstimation.h"
#include "EbAvcStyleMcp.h"
#include "aom_dsp_rtcd.h"
#if TX_SEARCH_LEVELS
#include "EbCodingLoop.h"
#endif

#define TH_NFL_BIAS             7
extern void av1_predict_intra_block_md(
    ModeDecisionContext_t       *cu_ptr,
    const Av1Common *cm,

    int32_t wpx,
    int32_t hpx,
    TxSize tx_size,
    PredictionMode mode,
    int32_t angle_delta,
    int32_t use_palette,
    FILTER_INTRA_MODE filter_intra_mode,
    uint8_t* topNeighArray,
    uint8_t* leftNeighArray,
    EbPictureBufferDesc_t  *reconBuffer,
    int32_t col_off,
    int32_t row_off,
    int32_t plane,
    BlockSize bsize,
    uint32_t cuOrgX,
    uint32_t cuOrgY,
    uint32_t OrgX,
    uint32_t OrgY
);
EbErrorType ProductGenerateMdCandidatesCu(
    LargestCodingUnit_t             *sb_ptr,
    ModeDecisionContext_t           *context_ptr,
    SsMeContext_t                  *ss_mecontext,
    const uint32_t                    leaf_index,

    const uint32_t                    lcuAddr,
    uint32_t                         *buffer_total_count,
    uint32_t                         *fastCandidateTotalCount,
    EbPtr                           interPredContextPtr,
    PictureControlSet_t            *picture_control_set_ptr);



void PfZeroOutUselessQuadrants(
    int16_t* transformCoeffBuffer,
    uint32_t  transformCoeffStride,
    uint32_t  quadrantSize,
    EbAsm  asm_type);



/*******************************************
* set Penalize Skip Flag
*
* Summary: Set the PenalizeSkipFlag to true
* When there is luminance/chrominance change
* or in noisy clip with low motion at meduim
* varince area
*
*******************************************/

#if CHROMA_BLIND
const EB_PREDICTION_FUNC  ProductPredictionFunTable[3] = { NULL, inter_pu_prediction_av1, AV1IntraPredictionCL};
#else
const EB_PREDICTION_FUNC  ProductPredictionFunTableCL[3][3] = {
    { NULL, inter2_nx2_n_pu_prediction_avc, AV1IntraPredictionCL },
    { NULL, inter2_nx2_n_pu_prediction_avc_style, AV1IntraPredictionCL },
    { NULL, inter_pu_prediction_av1, AV1IntraPredictionCL }
};
#endif

const EB_FAST_COST_FUNC   Av1ProductFastCostFuncTable[3] =
{
    NULL,

    Av1InterFastCost, /*INTER */
    Av1IntraFastCost /*INTRA */
};

const EB_AV1_FULL_COST_FUNC   Av1ProductFullCostFuncTable[3] =
{
    NULL,
    Av1InterFullCost, /*INTER */
    Av1IntraFullCost/*INTRA */

};


/***************************************************
* Update Recon Samples Neighbor Arrays
***************************************************/
void mode_decision_update_neighbor_arrays(
    ModeDecisionContext_t   *context_ptr,
    uint32_t                   index_mds,
    EbBool                  intraMdOpenLoop,
    EbBool                  intra4x4Selected
)

{
    uint32_t  bwdith = context_ptr->blk_geom->bwidth;
    uint32_t  bheight = context_ptr->blk_geom->bheight;

    uint32_t                   origin_x = context_ptr->cu_origin_x;
    uint32_t                   origin_y = context_ptr->cu_origin_y;
    (void)intra4x4Selected;

    uint32_t  cu_origin_x_uv = context_ptr->round_origin_x >> 1;
    uint32_t  cu_origin_y_uv = context_ptr->round_origin_y >> 1;
    uint32_t  bwdith_uv = context_ptr->blk_geom->bwidth_uv;
    uint32_t  bwheight_uv = context_ptr->blk_geom->bheight_uv;

    uint8_t modeType = context_ptr->cu_ptr->prediction_mode_flag;
    uint8_t intra_luma_mode = (uint8_t)context_ptr->cu_ptr->pred_mode;
    uint8_t chroma_mode = (uint8_t)context_ptr->cu_ptr->prediction_unit_array->intra_chroma_mode;
    uint8_t skip_flag = (uint8_t)context_ptr->cu_ptr->skip_flag;


    EbBool availableCoeff =
        (context_ptr->cu_ptr->transform_unit_array[0].y_has_coeff ||
            context_ptr->cu_ptr->transform_unit_array[0].v_has_coeff ||
            context_ptr->cu_ptr->transform_unit_array[0].u_has_coeff) ? EB_TRUE : EB_FALSE;

    uint8_t skipCoeff = !availableCoeff;


    context_ptr->mv_unit.predDirection = (uint8_t)(context_ptr->md_cu_arr_nsq[index_mds].prediction_unit_array[0].inter_pred_direction_index);
    context_ptr->mv_unit.mv[REF_LIST_0].mvUnion = context_ptr->md_cu_arr_nsq[index_mds].prediction_unit_array[0].mv[REF_LIST_0].mvUnion;
    context_ptr->mv_unit.mv[REF_LIST_1].mvUnion = context_ptr->md_cu_arr_nsq[index_mds].prediction_unit_array[0].mv[REF_LIST_1].mvUnion;
    MvUnit_t                *mv_unit = &context_ptr->mv_unit;
    (void)mv_unit;



    uint8_t                    y_has_coeff = context_ptr->cu_ptr->transform_unit_array[0].y_has_coeff;
    int32_t                   lumaDcCoeff = (int32_t)context_ptr->cu_ptr->quantized_dc[0];
    uint8_t                    u_has_coeff = context_ptr->cu_ptr->transform_unit_array[0].u_has_coeff;
    int32_t                   cbDcCoeff = (int32_t)context_ptr->cu_ptr->quantized_dc[1];
    uint8_t                    v_has_coeff = context_ptr->cu_ptr->transform_unit_array[0].v_has_coeff;
    int32_t                   crDcCoeff = (int32_t)context_ptr->cu_ptr->quantized_dc[2];
    uint8_t                    inter_pred_direction_index = (uint8_t)context_ptr->cu_ptr->prediction_unit_array->inter_pred_direction_index;
    uint8_t                    ref_frame_type = (uint8_t)context_ptr->cu_ptr->prediction_unit_array[0].ref_frame_type;


    NeighborArrayUnitModeWrite32(
        context_ptr->interpolation_type_neighbor_array,
        context_ptr->cu_ptr->interp_filters,
        origin_x,
        origin_y,
        bwdith,
        bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    {
        struct PartitionContext partition;
        partition.above = partition_context_lookup[context_ptr->blk_geom->bsize].above;
        partition.left = partition_context_lookup[context_ptr->blk_geom->bsize].left;

        NeighborArrayUnitModeWrite(
            context_ptr->leaf_partition_neighbor_array,
            (uint8_t*)(&partition), // NaderM
            origin_x,
            origin_y,
            bwdith,
            bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        // Mode Type Update
        NeighborArrayUnitModeWrite(
            context_ptr->mode_type_neighbor_array,
            &modeType,
            origin_x,
            origin_y,
            bwdith,
            bheight,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);

        // Intra Luma Mode Update
        NeighborArrayUnitModeWrite(
            context_ptr->intra_luma_mode_neighbor_array,
            &intra_luma_mode,//(uint8_t*)lumaMode,
            origin_x,
            origin_y,
            bwdith,
            bheight,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        {
            uint8_t dcSignCtx = 0;
            if (lumaDcCoeff > 0)
                dcSignCtx = 2;
            else if (lumaDcCoeff < 0)
                dcSignCtx = 1;
            else
                dcSignCtx = 0;
            uint8_t dcSignLevelCoeff = (uint8_t)((dcSignCtx << COEFF_CONTEXT_BITS) | y_has_coeff);
            if (!y_has_coeff)
                dcSignLevelCoeff = 0;

            NeighborArrayUnitModeWrite(
                context_ptr->luma_dc_sign_level_coeff_neighbor_array,
                (uint8_t*)&dcSignLevelCoeff,
                origin_x,
                origin_y,
                bwdith,
                bheight,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }

    }

    // Hsan: chroma mode rate estimation is kept even for chroma blind 
    if (context_ptr->blk_geom->has_uv) {

        // Intra Chroma Mode Update
        NeighborArrayUnitModeWrite(
            context_ptr->intra_chroma_mode_neighbor_array,
            &chroma_mode,
            cu_origin_x_uv,
            cu_origin_y_uv,
            bwdith_uv,
            bwheight_uv,

            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    NeighborArrayUnitModeWrite(
        context_ptr->skip_flag_neighbor_array,
        &skip_flag,
        origin_x,
        origin_y,
        bwdith,
        bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //  Update skip_coeff_neighbor_array,
    NeighborArrayUnitModeWrite(
        context_ptr->skip_coeff_neighbor_array,
        &skipCoeff,
        origin_x,
        origin_y,
        bwdith,
        bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

#if CHROMA_BLIND
    if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
    if (context_ptr->blk_geom->has_uv) {
#endif
        //  Update chroma CB cbf and Dc context
        {
            uint8_t dcSignCtx = 0;
            if (cbDcCoeff > 0)
                dcSignCtx = 2;
            else if (cbDcCoeff < 0)
                dcSignCtx = 1;
            else
                dcSignCtx = 0;
            uint8_t dcSignLevelCoeff = (uint8_t)((dcSignCtx << COEFF_CONTEXT_BITS) | u_has_coeff);
            if (!u_has_coeff)
                dcSignLevelCoeff = 0;

            NeighborArrayUnitModeWrite(
                context_ptr->cb_dc_sign_level_coeff_neighbor_array,
                (uint8_t*)&dcSignLevelCoeff,
                cu_origin_x_uv,
                cu_origin_y_uv,
                bwdith_uv,
                bwheight_uv,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }

        //  Update chroma CR cbf and Dc context
        {
            uint8_t dcSignCtx = 0;
            if (crDcCoeff > 0)
                dcSignCtx = 2;
            else if (crDcCoeff < 0)
                dcSignCtx = 1;
            else
                dcSignCtx = 0;
            uint8_t dcSignLevelCoeff = (uint8_t)((dcSignCtx << COEFF_CONTEXT_BITS) | v_has_coeff);
            if (!v_has_coeff)
                dcSignLevelCoeff = 0;

            NeighborArrayUnitModeWrite(
                context_ptr->cr_dc_sign_level_coeff_neighbor_array,
                (uint8_t*)&dcSignLevelCoeff,
                cu_origin_x_uv,
                cu_origin_y_uv,
                bwdith_uv,
                bwheight_uv,
                NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        }

    }

    // Update the Inter Pred Type Neighbor Array

    NeighborArrayUnitModeWrite(
        context_ptr->inter_pred_dir_neighbor_array,
        &inter_pred_direction_index,
        origin_x,
        origin_y,
        bwdith,
        bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    // Update the refFrame Type Neighbor Array
    NeighborArrayUnitModeWrite(
        context_ptr->ref_frame_type_neighbor_array,
        &ref_frame_type,
        origin_x,
        origin_y,
        bwdith,
        bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (intraMdOpenLoop == EB_FALSE)
    {
        update_recon_neighbor_array(
            context_ptr->luma_recon_neighbor_array,
            context_ptr->cu_ptr->neigh_top_recon[0],
            context_ptr->cu_ptr->neigh_left_recon[0],
            origin_x,
            origin_y,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight);

    }


    if (intraMdOpenLoop == EB_FALSE) {
#if CHROMA_BLIND
        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
        if (context_ptr->blk_geom->has_uv) {
#endif
            update_recon_neighbor_array(
                context_ptr->cb_recon_neighbor_array,
                context_ptr->cu_ptr->neigh_top_recon[1],
                context_ptr->cu_ptr->neigh_left_recon[1],
                cu_origin_x_uv,
                cu_origin_y_uv,
                bwdith_uv,
                bwheight_uv);
            update_recon_neighbor_array(
                context_ptr->cr_recon_neighbor_array,
                context_ptr->cu_ptr->neigh_top_recon[2],
                context_ptr->cu_ptr->neigh_left_recon[2],
                cu_origin_x_uv,
                cu_origin_y_uv,
                bwdith_uv,
                bwheight_uv);
        }
    }


    return;
}

void copy_neighbour_arrays(
    PictureControlSet_t                *picture_control_set_ptr,
    ModeDecisionContext_t               *context_ptr,
    uint32_t                            src_idx,
    uint32_t                            dst_idx,
    uint32_t                            blk_mds,
    uint32_t                            sb_org_x,
    uint32_t                            sb_org_y)
{
    (void)*context_ptr;

    const BlockGeom * blk_geom = Get_blk_geom_mds(blk_mds);

    uint32_t                            blk_org_x = sb_org_x + blk_geom->origin_x;
    uint32_t                            blk_org_y = sb_org_y + blk_geom->origin_y;
    uint32_t                            blk_org_x_uv = (blk_org_x >> 3 << 3) >> 1;
    uint32_t                            blk_org_y_uv = (blk_org_y >> 3 << 3) >> 1;
    uint32_t                            bwidth_uv = blk_geom->bwidth_uv;
    uint32_t                            bheight_uv = blk_geom->bheight_uv;

    copy_neigh_arr(
        picture_control_set_ptr->md_intra_luma_mode_neighbor_array[src_idx],
        picture_control_set_ptr->md_intra_luma_mode_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //NeighborArrayUnitReset(picture_control_set_ptr->md_intra_chroma_mode_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_intra_chroma_mode_neighbor_array[src_idx],
        picture_control_set_ptr->md_intra_chroma_mode_neighbor_array[dst_idx],
        blk_org_x_uv,
        blk_org_y_uv,
        bwidth_uv,
        bheight_uv,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //NeighborArrayUnitReset(picture_control_set_ptr->md_skip_flag_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_skip_flag_neighbor_array[src_idx],
        picture_control_set_ptr->md_skip_flag_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //NeighborArrayUnitReset(picture_control_set_ptr->md_mode_type_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_mode_type_neighbor_array[src_idx],
        picture_control_set_ptr->md_mode_type_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    //NeighborArrayUnitReset(picture_control_set_ptr->md_leaf_depth_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_leaf_depth_neighbor_array[src_idx],
        picture_control_set_ptr->md_leaf_depth_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    copy_neigh_arr(
        picture_control_set_ptr->mdleaf_partition_neighbor_array[src_idx],
        picture_control_set_ptr->mdleaf_partition_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    //NeighborArrayUnitReset(picture_control_set_ptr->md_luma_recon_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_luma_recon_neighbor_array[src_idx],
        picture_control_set_ptr->md_luma_recon_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

#if CHROMA_BLIND
    if (blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
    if (blk_geom->has_uv) {
#endif
        //NeighborArrayUnitReset(picture_control_set_ptr->md_cb_recon_neighbor_array[depth]);
        copy_neigh_arr(
            picture_control_set_ptr->md_cb_recon_neighbor_array[src_idx],
            picture_control_set_ptr->md_cb_recon_neighbor_array[dst_idx],
            blk_org_x_uv,
            blk_org_y_uv,
            bwidth_uv,
            bheight_uv,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);

        //NeighborArrayUnitReset(picture_control_set_ptr->md_cr_recon_neighbor_array[depth]);
        copy_neigh_arr(
            picture_control_set_ptr->md_cr_recon_neighbor_array[src_idx],
            picture_control_set_ptr->md_cr_recon_neighbor_array[dst_idx],
            blk_org_x_uv,
            blk_org_y_uv,
            bwidth_uv,
            bheight_uv,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);


    }

    //NeighborArrayUnitReset(picture_control_set_ptr->md_skip_coeff_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_skip_coeff_neighbor_array[src_idx],
        picture_control_set_ptr->md_skip_coeff_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    //NeighborArrayUnitReset(picture_control_set_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_luma_dc_sign_level_coeff_neighbor_array[src_idx],
        picture_control_set_ptr->md_luma_dc_sign_level_coeff_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

#if CHROMA_BLIND
    if (blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
    if (blk_geom->has_uv) {
#endif
        copy_neigh_arr(
            picture_control_set_ptr->md_cb_dc_sign_level_coeff_neighbor_array[src_idx],
            picture_control_set_ptr->md_cb_dc_sign_level_coeff_neighbor_array[dst_idx],
            blk_org_x_uv,
            blk_org_y_uv,
            bwidth_uv,
            bheight_uv,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
        //NeighborArrayUnitReset(picture_control_set_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth]);

        copy_neigh_arr(
            picture_control_set_ptr->md_cr_dc_sign_level_coeff_neighbor_array[src_idx],
            picture_control_set_ptr->md_cr_dc_sign_level_coeff_neighbor_array[dst_idx],
            blk_org_x_uv,
            blk_org_y_uv,
            bwidth_uv,
            bheight_uv,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }
    //NeighborArrayUnitReset(picture_control_set_ptr->md_inter_pred_dir_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_inter_pred_dir_neighbor_array[src_idx],
        picture_control_set_ptr->md_inter_pred_dir_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    //NeighborArrayUnitReset(picture_control_set_ptr->md_ref_frame_type_neighbor_array[depth]);
    copy_neigh_arr(
        picture_control_set_ptr->md_ref_frame_type_neighbor_array[src_idx],
        picture_control_set_ptr->md_ref_frame_type_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    copy_neigh_arr_32(
        picture_control_set_ptr->md_interpolation_type_neighbor_array[src_idx],
        picture_control_set_ptr->md_interpolation_type_neighbor_array[dst_idx],
        blk_org_x,
        blk_org_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
}

void md_update_all_neighbour_arrays(
    PictureControlSet_t                *picture_control_set_ptr,
    ModeDecisionContext_t               *context_ptr,
    uint32_t                             lastCuIndex_mds,
    uint32_t                            sb_origin_x,
    uint32_t                            sb_origin_y)
{

    context_ptr->blk_geom = Get_blk_geom_mds(lastCuIndex_mds);
    context_ptr->cu_origin_x = sb_origin_x + context_ptr->blk_geom->origin_x;
    context_ptr->cu_origin_y = sb_origin_y + context_ptr->blk_geom->origin_y;
    context_ptr->round_origin_x = ((context_ptr->cu_origin_x >> 3) << 3);
    context_ptr->round_origin_y = ((context_ptr->cu_origin_y >> 3) << 3);


    context_ptr->cu_ptr = &context_ptr->md_cu_arr_nsq[lastCuIndex_mds];

    mode_decision_update_neighbor_arrays(
        context_ptr,
        lastCuIndex_mds,
        picture_control_set_ptr->intra_md_open_loop_flag,
        EB_FALSE);

    update_mi_map(
#if CHROMA_BLIND
        context_ptr,
#endif
        context_ptr->cu_ptr,
        context_ptr->cu_origin_x,
        context_ptr->cu_origin_y,
        context_ptr->blk_geom,
        0,
        picture_control_set_ptr);



}

void md_update_all_neighbour_arrays_multiple(
    PictureControlSet_t                *picture_control_set_ptr,
    ModeDecisionContext_t               *context_ptr,
    uint32_t                            blk_mds,
    uint32_t                            sb_origin_x,
    uint32_t                            sb_origin_y){

    context_ptr->blk_geom = Get_blk_geom_mds(blk_mds);

    uint32_t blk_it;
    for (blk_it = 0; blk_it < context_ptr->blk_geom->totns; blk_it++)
    {

        md_update_all_neighbour_arrays(
            picture_control_set_ptr,
            context_ptr,
            blk_mds + blk_it,
            sb_origin_x,
            sb_origin_y);
    }
}

//*************************//
// set_nfl
// Based on the MDStage and the encodeMode
// the NFL candidates numbers are set
//*************************//
void set_nfl(
    ModeDecisionContext_t     *context_ptr,
    PictureControlSet_t       *picture_control_set_ptr){

    // Set NFL Candidates
    // NFL Level MD         Settings
    // 0                    MAX_NFL 12
    // 1                    8
    // 2                    6
    // 3                    4
    // 4                    4/3/2

    if (context_ptr->nfl_level == 0)
        context_ptr->full_recon_search_count = MAX_NFL;
    else if (context_ptr->nfl_level == 1)
#if TUNED_SETTINGS_FOR_M1
        context_ptr->full_recon_search_count = 10;
#else
        context_ptr->full_recon_search_count = 8;
#endif
    else if (context_ptr->nfl_level == 2)
        context_ptr->full_recon_search_count = 6;
    else if (context_ptr->nfl_level == 3)
        context_ptr->full_recon_search_count = 4;
    else
        context_ptr->full_recon_search_count =
            (picture_control_set_ptr->slice_type == I_SLICE) ? 4 :
            (context_ptr->blk_geom->bwidth == 32 && context_ptr->blk_geom->bheight == 32 && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 3 : 2;

        //if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_ptr->index] == LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE)
        //    context_ptr->full_recon_search_count = 1;
    ASSERT(context_ptr->full_recon_search_count <= MAX_NFL);
}

//*************************//
// SetNmm
// Based on the MDStage and the encodeMode
// the NMM candidates numbers are set
//*************************//
void Initialize_cu_data_structure(
    ModeDecisionContext_t   *context_ptr,
    SequenceControlSet_t    *sequence_control_set_ptr,
    LargestCodingUnit_t        *sb_ptr,
    const MdcLcuData_t        * const mdcResultTbPtr)
{
    UNUSED(*sequence_control_set_ptr);
    UNUSED(*sb_ptr);
    UNUSED(*mdcResultTbPtr);
    uint32_t blk_idx = 0;

    blk_idx = 0;
    do
    {
        const BlockGeom * blk_geom = Get_blk_geom_mds(blk_idx);

        if (blk_geom->shape == PART_N)
        {
            context_ptr->md_cu_arr_nsq[blk_idx].split_flag = EB_TRUE;  //this means that all the CUs at init time are not finals. only idd makes them final.
                                                                     //MDC would give us the split flag info for all the CUs, and we store in cu_ptr at the atrt of MD.
                                                                     //splitFalg=1 : to be tested CU.
                                                                     //split_flag=0 : to be tested CU + the CU could a final depth(smallest CU) or an invalid CU(out of pic bound)

            context_ptr->md_cu_arr_nsq[blk_idx].part = PARTITION_SPLIT;

            context_ptr->md_local_cu_unit[blk_idx].tested_cu_flag = EB_FALSE;
#if FIX_47
            //TODO: try to move this whole function to init
            context_ptr->md_cu_arr_nsq[blk_idx].mds_idx = blk_geom->blkidx_mds;
#endif
        }

        ++blk_idx;
    } while (blk_idx < sequence_control_set_ptr->max_block_cnt);
}

static INLINE tran_high_t check_range(tran_high_t input, int32_t bd) {
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
#endif  // CONFIG_COEFFICIENT_RANGE_CHECKING
    return (tran_high_t)clamp64(input, int_min, int_max);
}

#define HIGHBD_WRAPLOW(x, bd) ((int32_t)check_range((x), bd))
static INLINE uint16_t highbd_clip_pixel_add(uint16_t dest, tran_high_t trans,
    int32_t bd) {
    trans = HIGHBD_WRAPLOW(trans, bd);
    return clip_pixel_highbd(dest + (int32_t)trans, bd);
}

/*********************************
* Picture Single Channel Kernel
*********************************/
void PictureAdditionKernel(
    uint8_t  *predPtr,
    uint32_t  predStride,
    int32_t *residual_ptr,
    uint32_t  residualStride,
    uint8_t  *reconPtr,
    uint32_t  reconStride,
    uint32_t  width,
    uint32_t  height,
    int32_t     bd)
{
    uint32_t          columnIndex;
    uint32_t          rowIndex = 0;
    //    const int32_t    maxValue = 0xFF;

        //printf("\n");
        //printf("Reconstruction---------------------------------------------------\n");

    while (rowIndex < height) {

        columnIndex = 0;
        while (columnIndex < width) {
            //reconPtr[columnIndex] = (uint8_t)CLIP3(0, maxValue, ((int32_t)residual_ptr[columnIndex]) + ((int32_t)predPtr[columnIndex]));
            uint16_t rec = (uint16_t)predPtr[columnIndex];
            reconPtr[columnIndex] = (uint8_t)highbd_clip_pixel_add(rec, (tran_low_t)residual_ptr[columnIndex], bd);

            //printf("%d\t", reconPtr[columnIndex]);
            ++columnIndex;
        }

        //printf("\n");
        residual_ptr += residualStride;
        predPtr += predStride;
        reconPtr += reconStride;
        ++rowIndex;
    }
    //printf("-----------------------------------------------------------------\n");
    //printf("\n");
    //printf("\n");
    return;
}

void PictureAdditionKernel16Bit(
    uint16_t  *predPtr,
    uint32_t  predStride,
    int32_t *residual_ptr,
    uint32_t  residualStride,
    uint16_t  *reconPtr,
    uint32_t  reconStride,
    uint32_t  width,
    uint32_t  height,
    int32_t     bd)
{
    uint32_t          columnIndex;
    uint32_t          rowIndex = 0;
    //    const int32_t    maxValue = 0xFF;

        //printf("\n");
        //printf("Reconstruction---------------------------------------------------\n");

    while (rowIndex < height) {

        columnIndex = 0;
        while (columnIndex < width) {
            //reconPtr[columnIndex] = (uint8_t)CLIP3(0, maxValue, ((int32_t)residual_ptr[columnIndex]) + ((int32_t)predPtr[columnIndex]));
            uint16_t rec = (uint16_t)predPtr[columnIndex];
            reconPtr[columnIndex] = highbd_clip_pixel_add(rec, (tran_low_t)residual_ptr[columnIndex], bd);

            //printf("%d\t", reconPtr[columnIndex]);
            ++columnIndex;
        }

        //printf("\n");
        residual_ptr += residualStride;
        predPtr += predStride;
        reconPtr += reconStride;
        ++rowIndex;
    }
    //    printf("-----------------------------------------------------------------\n");
    //    printf("\n");
    //    printf("\n");
    return;
}

void AV1PerformInverseTransformReconLuma(
    PictureControlSet_t               *picture_control_set_ptr,
    ModeDecisionContext_t             *context_ptr,
    ModeDecisionCandidateBuffer_t     *candidateBuffer,
    CodingUnit_t                      *cu_ptr,
    const BlockGeom                   *blk_geom,
    EbAsm                              asm_type) {
    (void)cu_ptr;
    uint32_t                              tu_width;
    uint32_t                              tu_height;
    uint32_t                              txb_origin_x;
    uint32_t                              txb_origin_y;
    uint32_t                              tuOriginIndex;
    uint32_t                              tuTotalCount;

    uint32_t                              txb_itr;

    if (picture_control_set_ptr->intra_md_open_loop_flag == EB_FALSE) {
        tuTotalCount = blk_geom->txb_count;
        txb_itr = 0;
        uint32_t txb_1d_offset = 0;
        uint32_t recLumaOffset = (blk_geom->origin_y) * candidateBuffer->reconPtr->strideY +
            (blk_geom->origin_x);
        do {
            txb_origin_x = context_ptr->blk_geom->tx_org_x[txb_itr];
            txb_origin_y = context_ptr->blk_geom->tx_org_y[txb_itr];
            tu_width = context_ptr->blk_geom->tx_width[txb_itr];
            tu_height = context_ptr->blk_geom->tx_height[txb_itr];

            tuOriginIndex = txb_origin_x + txb_origin_y * candidateBuffer->prediction_ptr->strideY;

            uint32_t y_has_coeff = (candidateBuffer->candidate_ptr->y_has_coeff & (1 << txb_itr)) > 0;

            if (y_has_coeff) {
                (void)context_ptr;
                uint8_t     *predBuffer = &(candidateBuffer->prediction_ptr->bufferY[tuOriginIndex]);
                uint8_t     *recBuffer = &(candidateBuffer->reconPtr->bufferY[recLumaOffset]);

                uint32_t j;

                for (j = 0; j < tu_height; j++)
                    memcpy(recBuffer + j * candidateBuffer->reconPtr->strideY, predBuffer + j * candidateBuffer->prediction_ptr->strideY, tu_width);

                Av1InvTransformRecon8bit(

                    &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferY)[txb_1d_offset]),
                    recBuffer,
                    candidateBuffer->reconPtr->strideY,
                    context_ptr->blk_geom->txsize[txb_itr],
                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    (uint16_t)candidateBuffer->candidate_ptr->eob[0][txb_itr]);

            }
            else {

                PictureCopy8Bit(
                    candidateBuffer->prediction_ptr,
                    tuOriginIndex,
                    0,//tuChromaOriginIndex,
                    candidateBuffer->reconPtr,
                    recLumaOffset,
                    0,//tuChromaOriginIndex,
                    tu_width,
                    tu_height,
                    0,//chromaTuSize,
                    0,//chromaTuSize,
                    PICTURE_BUFFER_DESC_Y_FLAG,
                    asm_type);
            }

            txb_1d_offset += context_ptr->blk_geom->tx_width[txb_itr] * context_ptr->blk_geom->tx_height[txb_itr];
            ++txb_itr;

        } while (txb_itr < tuTotalCount);
    }
}
void AV1PerformInverseTransformRecon(
    PictureControlSet_t               *picture_control_set_ptr,
    ModeDecisionContext_t             *context_ptr,
    ModeDecisionCandidateBuffer_t     *candidateBuffer,
    CodingUnit_t                      *cu_ptr,
    const BlockGeom                   *blk_geom,
    EbAsm                              asm_type) {

    uint32_t                           tu_width;
    uint32_t                           tu_height;
    uint32_t                           txb_origin_x;
    uint32_t                           txb_origin_y;
    uint32_t                           tuOriginIndex;
    uint32_t                           tuTotalCount;
    uint32_t                           tu_index;
    uint32_t                           txb_itr;
    TransformUnit_t                   *txb_ptr;
    
    UNUSED(blk_geom);

    if (picture_control_set_ptr->intra_md_open_loop_flag == EB_FALSE) {
        tuTotalCount = context_ptr->blk_geom->txb_count;
        tu_index = 0;
        txb_itr = 0;
        uint32_t txb_1d_offset = 0, txb_1d_offset_uv = 0;
        uint32_t recLumaOffset, recCbOffset, recCrOffset;

        do {
            txb_origin_x = context_ptr->blk_geom->tx_org_x[txb_itr];
            txb_origin_y = context_ptr->blk_geom->tx_org_y[txb_itr];
            tu_width = context_ptr->blk_geom->tx_width[txb_itr];
            tu_height = context_ptr->blk_geom->tx_height[txb_itr];
            txb_ptr = &cu_ptr->transform_unit_array[tu_index];
            recLumaOffset = context_ptr->blk_geom->tx_org_x[txb_itr] + context_ptr->blk_geom->tx_org_y[txb_itr] * candidateBuffer->reconPtr->strideY;
            recCbOffset = ((((context_ptr->blk_geom->tx_org_x[txb_itr] >> 3) << 3) + ((context_ptr->blk_geom->tx_org_y[txb_itr] >> 3) << 3) * candidateBuffer->reconPtr->strideCb) >> 1);
            recCrOffset = ((((context_ptr->blk_geom->tx_org_x[txb_itr] >> 3) << 3) + ((context_ptr->blk_geom->tx_org_y[txb_itr] >> 3) << 3) * candidateBuffer->reconPtr->strideCr) >> 1);
            tuOriginIndex = txb_origin_x + txb_origin_y * candidateBuffer->prediction_ptr->strideY;
            if (txb_ptr->y_has_coeff) {
                uint8_t     *predBuffer = &(candidateBuffer->prediction_ptr->bufferY[tuOriginIndex]);
                uint8_t     *recBuffer = &(candidateBuffer->reconPtr->bufferY[recLumaOffset]);
                uint32_t     j;

                for (j = 0; j < tu_height; j++)
                    memcpy(recBuffer + j * candidateBuffer->reconPtr->strideY, predBuffer + j * candidateBuffer->prediction_ptr->strideY, tu_width);

                Av1InvTransformRecon8bit(
                    &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferY)[txb_1d_offset]),
                    recBuffer,
                    candidateBuffer->reconPtr->strideY,
                    context_ptr->blk_geom->txsize[txb_itr],
                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    (uint16_t)candidateBuffer->candidate_ptr->eob[0][txb_itr]);
            }
            else {
                PictureCopy8Bit(
                    candidateBuffer->prediction_ptr,
                    tuOriginIndex,
                    0,//tuChromaOriginIndex,
                    candidateBuffer->reconPtr,
                    recLumaOffset,
                    0,//tuChromaOriginIndex,
                    tu_width,
                    tu_height,
                    0,//chromaTuSize,
                    0,//chromaTuSize,
                    PICTURE_BUFFER_DESC_Y_FLAG,
                    asm_type);
            }
#if CHROMA_BLIND
            if (context_ptr->chroma_level == CHROMA_MODE_0) 
            {
#endif
                //CHROMA
                uint32_t chroma_tu_width = tx_size_wide[context_ptr->blk_geom->txsize_uv[txb_itr]];
                uint32_t chroma_tu_height = tx_size_high[context_ptr->blk_geom->txsize_uv[txb_itr]];
                uint32_t cbTuChromaOriginIndex = ((((txb_origin_x >> 3) << 3) + ((txb_origin_y >> 3) << 3) * candidateBuffer->reconCoeffPtr->strideCb) >> 1);
                uint32_t crTuChromaOriginIndex = ((((txb_origin_x >> 3) << 3) + ((txb_origin_y >> 3) << 3) * candidateBuffer->reconCoeffPtr->strideCr) >> 1);

                if (context_ptr->blk_geom->has_uv && txb_ptr->u_has_coeff) {

                    uint8_t     *predBuffer = &(candidateBuffer->prediction_ptr->bufferCb[cbTuChromaOriginIndex]);
                    uint8_t     *recBuffer = &(candidateBuffer->reconPtr->bufferCb[recCbOffset]);
                    uint32_t j;
                    for (j = 0; j < chroma_tu_height; j++)
                        memcpy(recBuffer + j * candidateBuffer->reconPtr->strideCb, predBuffer + j * candidateBuffer->prediction_ptr->strideCb, chroma_tu_width);

                    Av1InvTransformRecon8bit(
                        &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferCb)[txb_1d_offset_uv]),
                        recBuffer,
                        candidateBuffer->reconPtr->strideCb,
                        context_ptr->blk_geom->txsize_uv[txb_itr],
                        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                        PLANE_TYPE_UV,
                        (uint16_t)candidateBuffer->candidate_ptr->eob[1][txb_itr]);
                }
                else {

                    PictureCopy8Bit(
                        candidateBuffer->prediction_ptr,
                        0,
                        cbTuChromaOriginIndex,
                        candidateBuffer->reconPtr,
                        0,
                        recCbOffset,
                        0,
                        0,
                        chroma_tu_width,
                        chroma_tu_height,
                        PICTURE_BUFFER_DESC_Cb_FLAG,
                        asm_type);
                }

                if (context_ptr->blk_geom->has_uv && txb_ptr->v_has_coeff) {

                    uint8_t     *predBuffer = &(candidateBuffer->prediction_ptr->bufferCr[crTuChromaOriginIndex]);
                    uint8_t     *recBuffer = &(candidateBuffer->reconPtr->bufferCr[recCrOffset]);
                    uint32_t j;
                    for (j = 0; j < chroma_tu_height; j++)
                        memcpy(recBuffer + j * candidateBuffer->reconPtr->strideCr, predBuffer + j * candidateBuffer->prediction_ptr->strideCr, chroma_tu_width);

                    Av1InvTransformRecon8bit(
                        &(((int32_t*)candidateBuffer->reconCoeffPtr->bufferCr)[txb_1d_offset_uv]),

                        recBuffer,
                        candidateBuffer->reconPtr->strideCr,
                        context_ptr->blk_geom->txsize_uv[txb_itr],
                        candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV],
                        PLANE_TYPE_UV,
                        (uint16_t)candidateBuffer->candidate_ptr->eob[2][txb_itr]);
                }
                else {

                    PictureCopy8Bit(
                        candidateBuffer->prediction_ptr,
                        0,
                        crTuChromaOriginIndex,
                        candidateBuffer->reconPtr,
                        0,
                        recCrOffset,
                        0,
                        0,
                        chroma_tu_width,
                        chroma_tu_height,
                        PICTURE_BUFFER_DESC_Cr_FLAG,
                        asm_type);

                }
                //CHROMA END
#if CHROMA_BLIND
                if (context_ptr->blk_geom->has_uv)
                    txb_1d_offset_uv += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];
            }
#endif
            txb_1d_offset += context_ptr->blk_geom->tx_width[txb_itr] * context_ptr->blk_geom->tx_height[txb_itr];
#if !CHROMA_BLIND
            if (context_ptr->blk_geom->has_uv)
                txb_1d_offset_uv += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];
#endif
            ++tu_index;
            ++txb_itr;

        } while (txb_itr < tuTotalCount);
    }
}

/*******************************************
* Coding Loop - Fast Loop Initialization
*******************************************/
void ProductCodingLoopInitFastLoop(
    ModeDecisionContext_t      *context_ptr,
    NeighborArrayUnit_t        *skip_coeff_neighbor_array,
    NeighborArrayUnit_t        *luma_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *cb_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *cr_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *inter_pred_dir_neighbor_array,
    NeighborArrayUnit_t        *ref_frame_type_neighbor_array,
    NeighborArrayUnit_t        *intraLumaNeighborArray,
    NeighborArrayUnit_t        *skip_flag_neighbor_array,
    NeighborArrayUnit_t        *mode_type_neighbor_array,
    NeighborArrayUnit_t        *leaf_depth_neighbor_array,
    NeighborArrayUnit_t        *leaf_partition_neighbor_array
)
{
    // Keep track of the SB Ptr
    context_ptr->luma_intra_ref_samples_gen_done = EB_FALSE;
    context_ptr->chroma_intra_ref_samples_gen_done = EB_FALSE;

    // Generate Split, Skip and intra mode contexts for the rate estimation
    CodingLoopContextGeneration(
        context_ptr,
        context_ptr->cu_ptr,
        context_ptr->cu_origin_x,
        context_ptr->cu_origin_y,
        BLOCK_SIZE_64,
        skip_coeff_neighbor_array,
        luma_dc_sign_level_coeff_neighbor_array,
        cb_dc_sign_level_coeff_neighbor_array,
        cr_dc_sign_level_coeff_neighbor_array,
        inter_pred_dir_neighbor_array,
        ref_frame_type_neighbor_array,
        intraLumaNeighborArray,
        skip_flag_neighbor_array,
        mode_type_neighbor_array,
        leaf_depth_neighbor_array,
        leaf_partition_neighbor_array);

    // *Notes
    // -Instead of creating function pointers in a static array, put the func pointers in a queue
    // -This function should also do the SAD calc for each mode
    // -The best independent intra chroma mode should be determined here
    // -Modify PictureBufferDesc to be able to create Luma, Cb, and/or Cr only buffers via flags
    // -Have one PictureBufferDesc that points to the intra chroma buffer to be used.
    // -Properly signal the DM mode at this point

    // Initialize the candidate buffer costs
    {
        uint32_t buffer_depth_index_start = context_ptr->buffer_depth_index_start[0];
        uint32_t buffer_depth_index_width = context_ptr->buffer_depth_index_width[0];
        uint32_t index = 0;

        for (index = 0; index < buffer_depth_index_width; ++index) {
            context_ptr->fast_cost_array[buffer_depth_index_start + index] = 0xFFFFFFFFFFFFFFFFull;
            context_ptr->full_cost_array[buffer_depth_index_start + index] = 0xFFFFFFFFFFFFFFFFull;
        }
    }
    return;
}

uint64_t ProductGenerateChromaWeight(
    PictureControlSet_t                 *picture_control_set_ptr,
    uint32_t                               qp)
{
    uint64_t weight;

    if (picture_control_set_ptr->slice_type == I_SLICE) {
        weight = ChromaWeightFactorLd[qp];
    }
    else {
        // Random Access
        if (picture_control_set_ptr->temporal_layer_index == 0) {
            weight = ChromaWeightFactorRa[qp];
        }
        else if (picture_control_set_ptr->temporal_layer_index < 3) {
            weight = ChromaWeightFactorRaQpScalingL1[qp];
        }
        else {
            weight = ChromaWeightFactorRaQpScalingL3[qp];
        }
    }
    return (weight << 1);
}

uint64_t SpatialFullDistortionKernel(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *recon,
    uint32_t   reconStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight);

uint64_t SpatialFullDistortionKernel8x8_SSSE3_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *recon,
    uint32_t   reconStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight);

uint64_t SpatialFullDistortionKernel16MxN_SSSE3_INTRIN(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *recon,
    uint32_t   reconStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight);

void ProductMdFastPuPrediction(
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t       *candidateBuffer,
    ModeDecisionContext_t               *context_ptr,
#if !CHROMA_BLIND
    uint32_t                             use_chroma_information_in_fast_loop,
#endif
    uint32_t                             modeType,
    ModeDecisionCandidate_t             *const candidate_ptr,
    uint32_t                             fastLoopCandidateIndex,
    uint32_t                             bestFirstFastCostSearchCandidateIndex,
    EbAsm                                asm_type)
{
#if !CHROMA_BLIND
    UNUSED(use_chroma_information_in_fast_loop);
#endif
    UNUSED(candidate_ptr);
    UNUSED(fastLoopCandidateIndex);
    UNUSED(bestFirstFastCostSearchCandidateIndex);
    context_ptr->pu_itr = 0;
#if !CHROMA_BLIND
    EbBool enableSubPelFlag = picture_control_set_ptr->parent_pcs_ptr->use_subpel_flag;
    enableSubPelFlag = 2;
#endif
    // Prediction
#if INTERPOLATION_SEARCH_LEVELS
    context_ptr->skip_interpolation_search = picture_control_set_ptr->parent_pcs_ptr->interpolation_search_level == IT_SEARCH_FAST_LOOP ? 0 : 1;
#endif
    candidateBuffer->candidate_ptr->prediction_is_ready_luma = EB_TRUE;
    candidateBuffer->candidate_ptr->interp_filters = 0;

#if CHROMA_BLIND
    ProductPredictionFunTable[modeType](
        context_ptr,
        picture_control_set_ptr,
        candidateBuffer,
        asm_type);
#else
    ProductPredictionFunTableCL[enableSubPelFlag][modeType](
        context_ptr,
        context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
        picture_control_set_ptr,
        candidateBuffer,
        asm_type);
#endif
}
void generate_intra_reference_samples(
    const Av1Common         *cm,
    ModeDecisionContext_t   *md_context_ptr);

void ProductPerformFastLoop(
    PictureControlSet_t                 *picture_control_set_ptr,
    LargestCodingUnit_t                 *sb_ptr,
    ModeDecisionContext_t               *context_ptr,
    ModeDecisionCandidateBuffer_t      **candidateBufferPtrArrayBase,
    ModeDecisionCandidate_t             *fast_candidate_array,
    uint32_t                             fastCandidateTotalCount,
    EbPictureBufferDesc_t               *inputPicturePtr,
    uint32_t                             inputOriginIndex,
    uint32_t                             inputCbOriginIndex,
    uint32_t                             inputCrOriginIndex,
    CodingUnit_t                        *cu_ptr,
    uint32_t                             cuOriginIndex,
    uint32_t                             cuChromaOriginIndex,
    uint32_t                             maxBuffers,
    uint32_t                            *secondFastCostSearchCandidateTotalCount,
    EbAsm                                asm_type) {

    int32_t                          fastLoopCandidateIndex;
    uint64_t                          lumaFastDistortion;
    uint64_t                          chromaFastDistortion;
    ModeDecisionCandidateBuffer_t  *candidateBuffer;
    //    const EB_SLICE                  slice_type = picture_control_set_ptr->slice_type;
    uint32_t                          highestCostIndex;
    uint64_t                          highestCost;
    uint32_t                            isCandzz = 0;
    const uint8_t bwidth = context_ptr->blk_geom->bwidth;
    const uint8_t bheight = context_ptr->blk_geom->bheight;
    const BlockSize bsize = context_ptr->blk_geom->bsize;
    const uint8_t bwidth_uv = context_ptr->blk_geom->bwidth_uv;
    const uint8_t bheight_uv = context_ptr->blk_geom->bheight_uv;

    uint32_t firstFastCandidateTotalCount;
    // Initialize first fast cost loop variables
    uint64_t bestFirstFastCostSearchCandidateCost = 0xFFFFFFFFFFFFFFFFull;
    int32_t bestFirstFastCostSearchCandidateIndex = INVALID_FAST_CANDIDATE_INDEX;
    //    SequenceControlSet_t           *sequence_control_set_ptr = ((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr);

#if INTRA_CORE_OPT
    if (context_ptr->blk_geom->sq_size < 128) {
        generate_intra_reference_samples(
            picture_control_set_ptr->parent_pcs_ptr->av1_cm,
            context_ptr);
    }
#endif

    {

        firstFastCandidateTotalCount = 0;
        // First Fast-Cost Search Candidate Loop
        fastLoopCandidateIndex = fastCandidateTotalCount - 1;
        do
        {
            lumaFastDistortion = 0;

            // Set the Candidate Buffer
            candidateBuffer = candidateBufferPtrArrayBase[0];
            ModeDecisionCandidate_t *const candidate_ptr = candidateBuffer->candidate_ptr = &fast_candidate_array[fastLoopCandidateIndex];
            const unsigned distortion_ready = candidate_ptr->distortion_ready;


            // Only check (src - src) candidates (Tier0 candidates)
            if (!!distortion_ready)
            {
                const uint32_t type = candidate_ptr->type;

                lumaFastDistortion = candidate_ptr->me_distortion;
                firstFastCandidateTotalCount++;
                // Visual favor for DC/Planar
                //if( candidateBuffer->candidate_ptr->type == INTRA_MODE && (candidateBuffer->candidate_ptr->intra_luma_mode == EB_INTRA_PLANAR||candidateBuffer->candidate_ptr->intra_luma_mode == EB_INTRA_DC)){
                //    lumaFastDistortion = lumaFastDistortion * 90 / 100;
                //}

                {
                    // Fast Cost Calc
                    Av1ProductFastCostFuncTable[type](
                        context_ptr,
                        cu_ptr,
                        candidateBuffer,
                        cu_ptr->qp,
                        lumaFastDistortion,
                        0,
                        context_ptr->fast_lambda,
                        picture_control_set_ptr);

                    // Keep track of the candidate index of the best  (src - src) candidate
                    if (*(candidateBuffer->fast_cost_ptr) <= bestFirstFastCostSearchCandidateCost) {
                        bestFirstFastCostSearchCandidateIndex = fastLoopCandidateIndex;
                        bestFirstFastCostSearchCandidateCost = *(candidateBuffer->fast_cost_ptr);
                    }
                    // Initialize Fast Cost - to do not interact with the second Fast-Cost Search
                    *(candidateBuffer->fast_cost_ptr) = 0xFFFFFFFFFFFFFFFFull;
                }
            }
        } while (--fastLoopCandidateIndex >= 0);
    }

    // Second Fast-Cost Search Candidate Loop
    *secondFastCostSearchCandidateTotalCount = 0;
    highestCostIndex = context_ptr->buffer_depth_index_start[0];
    fastLoopCandidateIndex = fastCandidateTotalCount - 1;

    uint16_t                         lcuAddr = sb_ptr->index;
    do
    {
        candidateBuffer = candidateBufferPtrArrayBase[highestCostIndex];
        ModeDecisionCandidate_t *const  candidate_ptr = candidateBuffer->candidate_ptr = &fast_candidate_array[fastLoopCandidateIndex];
        const unsigned                  distortion_ready = candidate_ptr->distortion_ready;
        EbPictureBufferDesc_t * const   prediction_ptr = candidateBuffer->prediction_ptr;

        {
            candidateBuffer->sub_sampled_pred = EB_FALSE;
            candidateBuffer->sub_sampled_pred_chroma = EB_FALSE;

        }

        candidate_ptr->prediction_is_ready_luma = EB_FALSE;


        if ((!distortion_ready) || fastLoopCandidateIndex == bestFirstFastCostSearchCandidateIndex) {
#if !CHROMA_BLIND
            context_ptr->round_mv_to_integer = (candidate_ptr->merge_flag == EB_TRUE) ?

                EB_TRUE :
                EB_FALSE;

            context_ptr->round_mv_to_integer = picture_control_set_ptr->parent_pcs_ptr->use_subpel_flag ? EB_FALSE : context_ptr->round_mv_to_integer;
#endif
            lumaFastDistortion = 0;
            chromaFastDistortion = 0;
            // Set the Candidate Buffer

            ProductMdFastPuPrediction(
                picture_control_set_ptr,
                candidateBuffer,
                context_ptr, 
#if !CHROMA_BLIND
                EB_TRUE/*use_chroma_information_in_fast_loop*/,
#endif
                candidate_ptr->type,
                candidate_ptr,
                fastLoopCandidateIndex,
                bestFirstFastCostSearchCandidateIndex,
                asm_type);

            //Distortion
            uint8_t * const inputBufferY = inputPicturePtr->bufferY + inputOriginIndex;
            const unsigned inputStrideY = inputPicturePtr->strideY;
            uint8_t * const predBufferY = prediction_ptr->bufferY + cuOriginIndex;
            // Skip distortion computation if the candidate is MPM
            if (candidateBuffer->candidate_ptr->mpm_flag == EB_FALSE) {
                if (fastLoopCandidateIndex == bestFirstFastCostSearchCandidateIndex && candidate_ptr->type == INTRA_MODE)
                    lumaFastDistortion = candidate_ptr->me_distortion;
                else {
                    // Y
                    lumaFastDistortion += (NxMSadKernelSubSampled_funcPtrArray[asm_type][bwidth >> 3](
                        inputBufferY,
                        inputStrideY << candidateBuffer->sub_sampled_pred,
                        predBufferY,
                        prediction_ptr->strideY,
                        bheight >> candidateBuffer->sub_sampled_pred,
                        bwidth)) << candidateBuffer->sub_sampled_pred;
                }
#if CHROMA_BLIND
                if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
                // Cb
                if (context_ptr->blk_geom->has_uv) {
#endif

                    uint8_t * const inputBufferCb = inputPicturePtr->bufferCb + inputCbOriginIndex;
                    uint8_t *  const predBufferCb = candidateBuffer->prediction_ptr->bufferCb + cuChromaOriginIndex;

                    chromaFastDistortion += NxMSadKernelSubSampled_funcPtrArray[asm_type][bwidth >> 4](
                        inputBufferCb,
                        inputPicturePtr->strideCb << candidateBuffer->sub_sampled_pred_chroma,
                        predBufferCb,
                        prediction_ptr->strideCb,
                        bheight_uv >> candidateBuffer->sub_sampled_pred_chroma,
                        bwidth_uv) << candidateBuffer->sub_sampled_pred_chroma;


                    uint8_t * const inputBufferCr = inputPicturePtr->bufferCr + inputCrOriginIndex;
                    uint8_t * const predBufferCr = candidateBuffer->prediction_ptr->bufferCr + cuChromaOriginIndex;

                    chromaFastDistortion += NxMSadKernelSubSampled_funcPtrArray[asm_type][bwidth >> 4](
                        inputBufferCr,
                        inputPicturePtr->strideCb << candidateBuffer->sub_sampled_pred_chroma,
                        predBufferCr,
                        prediction_ptr->strideCr,
                        bheight_uv >> candidateBuffer->sub_sampled_pred_chroma,
                        bwidth_uv) << candidateBuffer->sub_sampled_pred_chroma;

                }
            }

            if (picture_control_set_ptr->parent_pcs_ptr->cmplx_status_sb[lcuAddr] == CMPLX_NOISE) {

                if (bsize == BLOCK_64X64 && candidate_ptr->type == INTER_MODE) { // Nader - to be reviewed for 128x128 sb

                    uint32_t  predDirection = (uint32_t)candidate_ptr->prediction_direction[0];
                    EbBool list0ZZ = (predDirection & 1) ? EB_TRUE : (EbBool)(candidate_ptr->motionVector_x_L0 == 0 && candidate_ptr->motionVector_y_L0 == 0);
                    EbBool list1ZZ = (predDirection > 0) ? (EbBool)(candidate_ptr->motionVector_x_L1 == 0 && candidate_ptr->motionVector_y_L1 == 0) : EB_TRUE;

                    isCandzz = (list0ZZ && list1ZZ) ? 1 : 0;
                    chromaFastDistortion = isCandzz ? chromaFastDistortion >> 2 : chromaFastDistortion;
                }

            }

            // Fast Cost Calc
            Av1ProductFastCostFuncTable[candidate_ptr->type](
                context_ptr,
                cu_ptr,
                candidateBuffer,
                cu_ptr->qp,
                lumaFastDistortion,
                chromaFastDistortion,
                context_ptr->fast_lambda,
                picture_control_set_ptr);
            (*secondFastCostSearchCandidateTotalCount)++;
        }

        // Find the buffer with the highest cost
        if (fastLoopCandidateIndex)
        {
            // maxCost is volatile to prevent the compiler from loading 0xFFFFFFFFFFFFFF
            //   as a const at the early-out. Loading a large constant on intel x64 processors
            //   clogs the i-cache/intstruction decode. This still reloads the variable from
            //   the stack each pass, so a better solution would be to register the variable,
            //   but this might require asm.

            volatile uint64_t maxCost = ~0ull;
            const uint64_t *fast_cost_array = context_ptr->fast_cost_array;
            const uint32_t bufferIndexStart = context_ptr->buffer_depth_index_start[0];
            const uint32_t bufferIndexEnd = bufferIndexStart + maxBuffers;
            uint32_t bufferIndex;

            highestCostIndex = bufferIndexStart;
            bufferIndex = bufferIndexStart + 1;

            do {
                highestCost = fast_cost_array[highestCostIndex];
                if (highestCost == maxCost)
                    break;

                if (fast_cost_array[bufferIndex] > highestCost)
                    highestCostIndex = bufferIndex;

            } while (++bufferIndex < bufferIndexEnd);
        }
    } while (--fastLoopCandidateIndex >= 0);// End Second FastLoop

}

void ProductConfigureChroma(
    PictureControlSet_t                 *picture_control_set_ptr,
    ModeDecisionContext_t               *context_ptr,
    LargestCodingUnit_t                 *sb_ptr) {

    uint32_t  lcuAddr = sb_ptr->index;
    uint32_t  lcuEdgeNum = picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[lcuAddr].edge_block_num;
    uint64_t  chroma_weight = 1;
    UNUSED(lcuEdgeNum);

    chroma_weight = ProductGenerateChromaWeight(
        picture_control_set_ptr,
        sb_ptr->qp);

    context_ptr->chroma_weight = (picture_control_set_ptr->parent_pcs_ptr->failing_motion_sb_flag[lcuAddr]) ? chroma_weight << 1 : chroma_weight;
}

void ProductDerivePartialFrequencyN2Flag(
    SequenceControlSet_t               *sequence_control_set_ptr,
    PictureControlSet_t                *picture_control_set_ptr,
    ModeDecisionContext_t              *context_ptr)
{
#if ENCODER_MODE_CLEANUP
    context_ptr->pf_md_mode = PF_OFF;
    UNUSED(sequence_control_set_ptr);
    UNUSED(picture_control_set_ptr);
#else
    if (sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE)
    {
        if (picture_control_set_ptr->parent_pcs_ptr->enc_mode == ENC_M2)
        {
            context_ptr->pf_md_mode = (picture_control_set_ptr->temporal_layer_index > 0) ? PF_N2 : PF_OFF;
        }

        else if (picture_control_set_ptr->parent_pcs_ptr->enc_mode >= ENC_M3 && picture_control_set_ptr->parent_pcs_ptr->enc_mode < ENC_M6)
        {
            context_ptr->pf_md_mode = (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE) ? PF_N2 : (context_ptr->blk_geom->bwidth <= 8 && context_ptr->blk_geom->bheight <= 8) ? PF_N2 : PF_N4;
        }
        else if (picture_control_set_ptr->parent_pcs_ptr->enc_mode >= ENC_M6)

        {

            if ((picture_control_set_ptr->slice_type == I_SLICE) || (picture_control_set_ptr->parent_pcs_ptr->uncovered_area_sb_flag[context_ptr->sb_ptr->index])) {
                context_ptr->pf_md_mode = PF_OFF;
            }

            else if (((picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE) && (picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[context_ptr->sb_ptr->index].edge_block_num > 0)) || (context_ptr->blk_geom->bwidth <= 8 && context_ptr->blk_geom->bheight <= 8)) {
                context_ptr->pf_md_mode = PF_N2;
            }

            else
                context_ptr->pf_md_mode = PF_N4;

        }
        else
            context_ptr->pf_md_mode = PF_OFF;

    }
    else
    {
        if (picture_control_set_ptr->parent_pcs_ptr->enc_mode <= ENC_M3 && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_FALSE)
            context_ptr->pf_md_mode = PF_N2;
        else if (picture_control_set_ptr->parent_pcs_ptr->enc_mode < ENC_M6 && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_FALSE)
            context_ptr->pf_md_mode = (context_ptr->blk_geom->bwidth <= 16 && context_ptr->blk_geom->bheight <= 16) ? PF_N2 : PF_N4;
        else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_FALSE)
            context_ptr->pf_md_mode = (context_ptr->blk_geom->bwidth <= 8 && context_ptr->blk_geom->bheight <= 8) ? PF_N2 : PF_N4;
        else
            context_ptr->pf_md_mode = PF_OFF;

    }
#endif
    context_ptr->pf_md_mode = PF_OFF;


}

void AV1CostCalcCfl(
    PictureControlSet_t                *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t      *candidateBuffer,
    LargestCodingUnit_t                *sb_ptr,
    ModeDecisionContext_t              *context_ptr,
    uint32_t                            component_mask,
    EbPictureBufferDesc_t              *inputPicturePtr,
    uint32_t                            inputCbOriginIndex,
    uint32_t                            cuChromaOriginIndex,
    uint64_t                            full_distortion[DIST_CALC_TOTAL],
    uint64_t                           *coeffBits,
    EbAsm                               asm_type) {

    ModeDecisionCandidate_t            *candidate_ptr = candidateBuffer->candidate_ptr;
    uint32_t                            count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];
    uint64_t                            cbFullDistortion[DIST_CALC_TOTAL];
    uint64_t                            crFullDistortion[DIST_CALC_TOTAL];
    uint64_t                            cb_coeff_bits = 0;
    uint64_t                            cr_coeff_bits = 0;
    uint32_t                            chroma_width = context_ptr->blk_geom->bwidth_uv;
    uint32_t                            chroma_height = context_ptr->blk_geom->bheight_uv;
    // FullLoop and TU search
    int32_t                             alpha_q3;
    uint8_t                             cbQp = context_ptr->qp;
    uint8_t                             crQp = context_ptr->qp;

    full_distortion[DIST_CALC_RESIDUAL] = 0;
    full_distortion[DIST_CALC_PREDICTION] = 0;
    *coeffBits = 0;
    
    // Loop over alphas and find the best
    if (component_mask == COMPONENT_CHROMA_CB || component_mask == COMPONENT_CHROMA || component_mask == COMPONENT_ALL) {
        cbFullDistortion[DIST_CALC_RESIDUAL] = 0;
        crFullDistortion[DIST_CALC_RESIDUAL] = 0;
        cbFullDistortion[DIST_CALC_PREDICTION] = 0;
        crFullDistortion[DIST_CALC_PREDICTION] = 0;
        cb_coeff_bits = 0;
        cr_coeff_bits = 0;

        alpha_q3 =
            cfl_idx_to_alpha(candidate_ptr->cfl_alpha_idx, candidate_ptr->cfl_alpha_signs, CFL_PRED_U); // once for U, once for V
        if (candidate_ptr->cfl_alpha_idx == 0 && candidate_ptr->cfl_alpha_signs == 0)// To check DC
            alpha_q3 = 0;

        assert(chroma_width * CFL_BUF_LINE + chroma_height <=
            CFL_BUF_SQUARE);

        cfl_predict_lbd(
            context_ptr->pred_buf_q3,
            &(candidateBuffer->prediction_ptr->bufferCb[cuChromaOriginIndex]),
            candidateBuffer->prediction_ptr->strideCb,
            //dst_16,
            &(candidateBuffer->cflTempPredictionPtr->bufferCb[cuChromaOriginIndex]),
            candidateBuffer->cflTempPredictionPtr->strideCb,
            alpha_q3,
            8,
            chroma_width,
            chroma_height);
        //Cb Residual

        ResidualKernel(
            &(inputPicturePtr->bufferCb[inputCbOriginIndex]),
            inputPicturePtr->strideCb,
            &(candidateBuffer->cflTempPredictionPtr->bufferCb[cuChromaOriginIndex]),
            candidateBuffer->cflTempPredictionPtr->strideCb,
            &(((int16_t*)candidateBuffer->residual_ptr->bufferCb)[cuChromaOriginIndex]),
            candidateBuffer->residual_ptr->strideCb,
            chroma_width,
            chroma_height);

        FullLoop_R(
            sb_ptr,
            candidateBuffer,
            context_ptr,
            inputPicturePtr,
            picture_control_set_ptr,
            PICTURE_BUFFER_DESC_Cb_FLAG,
            cbQp,
            crQp,
            &(*count_non_zero_coeffs[1]),
            &(*count_non_zero_coeffs[2]));


        // Create new function
        CuFullDistortionFastTuMode_R(
            sb_ptr,
            candidateBuffer,
            context_ptr,
            candidate_ptr,
            picture_control_set_ptr,
            cbFullDistortion,
            crFullDistortion,
            count_non_zero_coeffs,
            COMPONENT_CHROMA_CB,
            &cb_coeff_bits,
            &cr_coeff_bits,
            asm_type);

        full_distortion[DIST_CALC_RESIDUAL] += cbFullDistortion[DIST_CALC_RESIDUAL];
        full_distortion[DIST_CALC_PREDICTION] += cbFullDistortion[DIST_CALC_PREDICTION];
        *coeffBits += cb_coeff_bits;

    }
    if (component_mask == COMPONENT_CHROMA_CR || component_mask == COMPONENT_CHROMA || component_mask == COMPONENT_ALL) {

        cbFullDistortion[DIST_CALC_RESIDUAL] = 0;
        crFullDistortion[DIST_CALC_RESIDUAL] = 0;
        cbFullDistortion[DIST_CALC_PREDICTION] = 0;
        crFullDistortion[DIST_CALC_PREDICTION] = 0;

        cb_coeff_bits = 0;
        cr_coeff_bits = 0;

        alpha_q3 =
            cfl_idx_to_alpha(candidate_ptr->cfl_alpha_idx, candidate_ptr->cfl_alpha_signs, CFL_PRED_V); // once for U, once for V

        if (candidate_ptr->cfl_alpha_idx == 0 && candidate_ptr->cfl_alpha_signs == 0) // To check DC
            alpha_q3 = 0;

        assert(chroma_width * CFL_BUF_LINE + chroma_height <=
            CFL_BUF_SQUARE);

        cfl_predict_lbd(
            context_ptr->pred_buf_q3,
            &(candidateBuffer->prediction_ptr->bufferCr[cuChromaOriginIndex]),
            candidateBuffer->prediction_ptr->strideCr,
            //dst_16,
            &(candidateBuffer->cflTempPredictionPtr->bufferCr[cuChromaOriginIndex]),
            candidateBuffer->cflTempPredictionPtr->strideCr,
            alpha_q3,
            8,
            chroma_width,
            chroma_height);

        //Cr Residual
        ResidualKernel(
            &(inputPicturePtr->bufferCr[inputCbOriginIndex]),
            inputPicturePtr->strideCr,
            &(candidateBuffer->cflTempPredictionPtr->bufferCr[cuChromaOriginIndex]),
            candidateBuffer->cflTempPredictionPtr->strideCr,
            &(((int16_t*)candidateBuffer->residual_ptr->bufferCr)[cuChromaOriginIndex]),
            candidateBuffer->residual_ptr->strideCr,
            chroma_width,
            chroma_height);

        FullLoop_R(
            sb_ptr,
            candidateBuffer,
            context_ptr,
            inputPicturePtr,
            picture_control_set_ptr,
            PICTURE_BUFFER_DESC_Cr_FLAG,
            cbQp,
            crQp,
            &(*count_non_zero_coeffs[1]),
            &(*count_non_zero_coeffs[2]));
        candidate_ptr->v_has_coeff = *count_non_zero_coeffs[2] ? EB_TRUE : EB_FALSE;

        // Create new function
        CuFullDistortionFastTuMode_R(
            sb_ptr,
            candidateBuffer,
            context_ptr,
            candidate_ptr,
            picture_control_set_ptr,
            cbFullDistortion,
            crFullDistortion,
            count_non_zero_coeffs,
            COMPONENT_CHROMA_CR,
            &cb_coeff_bits,
            &cr_coeff_bits,
            asm_type);

        full_distortion[DIST_CALC_RESIDUAL] += crFullDistortion[DIST_CALC_RESIDUAL];
        full_distortion[DIST_CALC_PREDICTION] += crFullDistortion[DIST_CALC_PREDICTION];
        *coeffBits += cr_coeff_bits;
    }
}

#define PLANE_SIGN_TO_JOINT_SIGN(plane, a, b) \
  (plane == CFL_PRED_U ? a * CFL_SIGNS + b - 1 : b * CFL_SIGNS + a - 1)
/*************************Pick the best alpha for cfl mode  or Choose DC******************************************************/
#if CHROMA_BLIND 
void cfl_rd_pick_alpha(
#else
static void cfl_rd_pick_alpha(
#endif
    PictureControlSet_t     *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    LargestCodingUnit_t     *sb_ptr,
    ModeDecisionContext_t   *context_ptr,
    EbPictureBufferDesc_t   *inputPicturePtr,
    uint32_t                   inputCbOriginIndex,
    uint32_t                     cuChromaOriginIndex,
    EbAsm                    asm_type) {

    int64_t                  best_rd = INT64_MAX;
    uint64_t                  full_distortion[DIST_CALC_TOTAL];
    uint64_t                  coeffBits;

    const int64_t mode_rd =
        RDCOST(context_ptr->full_lambda,
        (uint64_t)candidateBuffer->candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[CFL_ALLOWED][candidateBuffer->candidate_ptr->intra_luma_mode][UV_CFL_PRED], 0);

    int64_t best_rd_uv[CFL_JOINT_SIGNS][CFL_PRED_PLANES];
    int32_t best_c[CFL_JOINT_SIGNS][CFL_PRED_PLANES];

    for (int32_t plane = 0; plane < CFL_PRED_PLANES; plane++) {
        coeffBits = 0;
        full_distortion[DIST_CALC_RESIDUAL] = 0;
        for (int32_t joint_sign = 0; joint_sign < CFL_JOINT_SIGNS; joint_sign++) {
            best_rd_uv[joint_sign][plane] = INT64_MAX;
            best_c[joint_sign][plane] = 0;
        }
        // Collect RD stats for an alpha value of zero in this plane.
        // Skip i == CFL_SIGN_ZERO as (0, 0) is invalid.
        for (int32_t i = CFL_SIGN_NEG; i < CFL_SIGNS; i++) {
            const int32_t joint_sign = PLANE_SIGN_TO_JOINT_SIGN(plane, CFL_SIGN_ZERO, i);
            if (i == CFL_SIGN_NEG) {
                candidateBuffer->candidate_ptr->cfl_alpha_idx = 0;
                candidateBuffer->candidate_ptr->cfl_alpha_signs = joint_sign;

                AV1CostCalcCfl(
                    picture_control_set_ptr,
                    candidateBuffer,
                    sb_ptr,
                    context_ptr,
                    (plane == 0) ? COMPONENT_CHROMA_CB : COMPONENT_CHROMA_CR,
                    inputPicturePtr,
                    inputCbOriginIndex,
                    cuChromaOriginIndex,
                    full_distortion,
                    &coeffBits,
                    asm_type);

                if (coeffBits == INT64_MAX) break;
            }

            const int32_t alpha_rate = candidateBuffer->candidate_ptr->md_rate_estimation_ptr->cflAlphaFacBits[joint_sign][plane][0];

            best_rd_uv[joint_sign][plane] =
                RDCOST(context_ptr->full_lambda, coeffBits + alpha_rate, full_distortion[DIST_CALC_RESIDUAL]);
        }
    }

    int32_t best_joint_sign = -1;

    for (int32_t plane = 0; plane < CFL_PRED_PLANES; plane++) {
        for (int32_t pn_sign = CFL_SIGN_NEG; pn_sign < CFL_SIGNS; pn_sign++) {
            int32_t progress = 0;
            for (int32_t c = 0; c < CFL_ALPHABET_SIZE; c++) {
                int32_t flag = 0;
                if (c > 2 && progress < c) break;
                coeffBits = 0;
                full_distortion[DIST_CALC_RESIDUAL] = 0;
                for (int32_t i = 0; i < CFL_SIGNS; i++) {
                    const int32_t joint_sign = PLANE_SIGN_TO_JOINT_SIGN(plane, pn_sign, i);
                    if (i == 0) {
                        candidateBuffer->candidate_ptr->cfl_alpha_idx = (c << CFL_ALPHABET_SIZE_LOG2) + c;
                        candidateBuffer->candidate_ptr->cfl_alpha_signs = joint_sign;

                        AV1CostCalcCfl(
                            picture_control_set_ptr,
                            candidateBuffer,
                            sb_ptr,
                            context_ptr,
                            (plane == 0) ? COMPONENT_CHROMA_CB : COMPONENT_CHROMA_CR,
                            inputPicturePtr,
                            inputCbOriginIndex,
                            cuChromaOriginIndex,
                            full_distortion,
                            &coeffBits,
                            asm_type);

                        if (coeffBits == INT64_MAX) break;
                    }

                    const int32_t alpha_rate = candidateBuffer->candidate_ptr->md_rate_estimation_ptr->cflAlphaFacBits[joint_sign][plane][c];

                    int64_t this_rd =
                        RDCOST(context_ptr->full_lambda, coeffBits + alpha_rate, full_distortion[DIST_CALC_RESIDUAL]);
                    if (this_rd >= best_rd_uv[joint_sign][plane]) continue;
                    best_rd_uv[joint_sign][plane] = this_rd;
                    best_c[joint_sign][plane] = c;

                    flag = 2;
                    if (best_rd_uv[joint_sign][!plane] == INT64_MAX) continue;
                    this_rd += mode_rd + best_rd_uv[joint_sign][!plane];
                    if (this_rd >= best_rd) continue;
                    best_rd = this_rd;
                    best_joint_sign = joint_sign;
                }
                progress += flag;
            }
        }
    }

    // Compare with DC Chroma
    coeffBits = 0;
    full_distortion[DIST_CALC_RESIDUAL] = 0;

    candidateBuffer->candidate_ptr->cfl_alpha_idx = 0;
    candidateBuffer->candidate_ptr->cfl_alpha_signs = 0;

    const int64_t dc_mode_rd =
        RDCOST(context_ptr->full_lambda,
            candidateBuffer->candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[CFL_ALLOWED][candidateBuffer->candidate_ptr->intra_luma_mode][UV_DC_PRED], 0);

    AV1CostCalcCfl(
        picture_control_set_ptr,
        candidateBuffer,
        sb_ptr,
        context_ptr,
        COMPONENT_CHROMA,
        inputPicturePtr,
        inputCbOriginIndex,
        cuChromaOriginIndex,
        full_distortion,
        &coeffBits,
        asm_type);

    int64_t dc_rd =
        RDCOST(context_ptr->full_lambda, coeffBits, full_distortion[DIST_CALC_RESIDUAL]);


    dc_rd += dc_mode_rd;
    if (dc_rd <= best_rd) {
        candidateBuffer->candidate_ptr->intra_chroma_mode = UV_DC_PRED;
        candidateBuffer->candidate_ptr->cfl_alpha_idx = 0;
        candidateBuffer->candidate_ptr->cfl_alpha_signs = 0;
    }
    else {
#if CHROMA_BLIND
        candidateBuffer->candidate_ptr->intra_chroma_mode = UV_CFL_PRED;
#endif
        int32_t ind = 0;
        if (best_joint_sign >= 0) {
            const int32_t u = best_c[best_joint_sign][CFL_PRED_U];
            const int32_t v = best_c[best_joint_sign][CFL_PRED_V];
            ind = (u << CFL_ALPHABET_SIZE_LOG2) + v;
        }
        else {
            best_joint_sign = 0;
        }
        candidateBuffer->candidate_ptr->cfl_alpha_idx = ind;
        candidateBuffer->candidate_ptr->cfl_alpha_signs = best_joint_sign;
    }

}



// If mode is CFL:
// 1: recon the Luma
// 2: Form the pred_buf_q3
// 3: Loop over alphas and find the best or choose DC
// 4: Recalculate the residual for chroma
static void CflPrediction(
    PictureControlSet_t     *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t  *candidateBuffer,
    LargestCodingUnit_t     *sb_ptr,
    ModeDecisionContext_t   *context_ptr,
    EbPictureBufferDesc_t   *inputPicturePtr,
    uint32_t                   inputCbOriginIndex,
    uint32_t                     cuChromaOriginIndex,
    EbAsm                    asm_type)
{
    // 1: recon the Luma
    AV1PerformInverseTransformReconLuma(
        picture_control_set_ptr,
        context_ptr,
        candidateBuffer,
        context_ptr->cu_ptr,
        context_ptr->blk_geom,
        asm_type);


    // 2: Form the pred_buf_q3
    uint32_t chroma_width = context_ptr->blk_geom->bwidth_uv;
    uint32_t chroma_height = context_ptr->blk_geom->bheight_uv;

    uint32_t recLumaOffset = (context_ptr->blk_geom->origin_y) * candidateBuffer->reconPtr->strideY +
        (context_ptr->blk_geom->origin_x);

    // Down sample Luma
    cfl_luma_subsampling_420_lbd_c(
        &(candidateBuffer->reconPtr->bufferY[recLumaOffset]),
        candidateBuffer->reconPtr->strideY,
        context_ptr->pred_buf_q3,
        context_ptr->blk_geom->bwidth,
        context_ptr->blk_geom->bheight);


    int32_t round_offset = chroma_width * chroma_height / 2;

    subtract_average(
        context_ptr->pred_buf_q3,
        chroma_width,
        chroma_height,
        round_offset,
        LOG2F(chroma_width) + LOG2F(chroma_height));


    // 3: Loop over alphas and find the best or choose DC
    cfl_rd_pick_alpha(
        picture_control_set_ptr,
        candidateBuffer,
        sb_ptr,
        context_ptr,
        inputPicturePtr,
        inputCbOriginIndex,
        cuChromaOriginIndex,
        asm_type);


    if (candidateBuffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {

        // 4: Recalculate the prediction and the residual
        int32_t alpha_q3_cb =
            cfl_idx_to_alpha(candidateBuffer->candidate_ptr->cfl_alpha_idx, candidateBuffer->candidate_ptr->cfl_alpha_signs, CFL_PRED_U);
        int32_t alpha_q3_cr =
            cfl_idx_to_alpha(candidateBuffer->candidate_ptr->cfl_alpha_idx, candidateBuffer->candidate_ptr->cfl_alpha_signs, CFL_PRED_V);

        assert(chroma_height * CFL_BUF_LINE + chroma_width <=
            CFL_BUF_SQUARE);

        cfl_predict_lbd(
            context_ptr->pred_buf_q3,
            &(candidateBuffer->prediction_ptr->bufferCb[cuChromaOriginIndex]),
            candidateBuffer->prediction_ptr->strideCb,
            &(candidateBuffer->prediction_ptr->bufferCb[cuChromaOriginIndex]),
            candidateBuffer->prediction_ptr->strideCb,
            alpha_q3_cb,
            8,
            chroma_width,
            chroma_height);

        cfl_predict_lbd(
            context_ptr->pred_buf_q3,
            &(candidateBuffer->prediction_ptr->bufferCr[cuChromaOriginIndex]),
            candidateBuffer->prediction_ptr->strideCr,
            &(candidateBuffer->prediction_ptr->bufferCr[cuChromaOriginIndex]),
            candidateBuffer->prediction_ptr->strideCr,
            alpha_q3_cr,
            8,
            chroma_width,
            chroma_height);


        //Cb Residual
        ResidualKernel(
            &(inputPicturePtr->bufferCb[inputCbOriginIndex]),
            inputPicturePtr->strideCb,
            &(candidateBuffer->prediction_ptr->bufferCb[cuChromaOriginIndex]),
            candidateBuffer->prediction_ptr->strideCb,
            &(((int16_t*)candidateBuffer->residual_ptr->bufferCb)[cuChromaOriginIndex]),
            candidateBuffer->residual_ptr->strideCb,
            context_ptr->blk_geom->bwidth_uv,
            context_ptr->blk_geom->bheight_uv);


        //Cr Residual
        ResidualKernel(
            &(inputPicturePtr->bufferCr[inputCbOriginIndex]),
            inputPicturePtr->strideCr,
            &(candidateBuffer->prediction_ptr->bufferCr[cuChromaOriginIndex]),
            candidateBuffer->prediction_ptr->strideCr,
            &(((int16_t*)candidateBuffer->residual_ptr->bufferCr)[cuChromaOriginIndex]),
            candidateBuffer->residual_ptr->strideCr,
            context_ptr->blk_geom->bwidth_uv,
            context_ptr->blk_geom->bheight_uv);
    }
    else {
        // Alphas = 0, Preds are the same as DC. Switch to DC mode
        candidateBuffer->candidate_ptr->intra_chroma_mode = UV_DC_PRED;
    }
}
#if TX_SEARCH_LEVELS
uint8_t get_skip_tx_search_flag(
    int32_t                  sq_size,
    uint64_t                 ref_fast_cost,
    uint64_t                 cu_cost,
    uint64_t                 weight)
{
    //NM: Skip tx search when the fast cost of the current mode candidate is substansially 
    // Larger than the best fast_cost (
    uint8_t  tx_search_skip_fag = cu_cost >= ((ref_fast_cost * weight) / 100) ? 1 : 0;
    tx_search_skip_fag = sq_size >= 128 ? 1 : tx_search_skip_fag;

    return tx_search_skip_fag;
}
#endif

void AV1PerformFullLoop(
    PictureControlSet_t     *picture_control_set_ptr,
    LargestCodingUnit_t     *sb_ptr,
    CodingUnit_t            *cu_ptr,
    ModeDecisionContext_t   *context_ptr,
    EbPictureBufferDesc_t   *inputPicturePtr,
    uint32_t                 inputOriginIndex,
    uint32_t                 inputCbOriginIndex,
    uint32_t                 cuOriginIndex,
    uint32_t                 cuChromaOriginIndex,
    uint32_t                 fullCandidateTotalCount,
#if TX_SEARCH_LEVELS
    uint64_t                 ref_fast_cost,
#endif
    EbAsm                    asm_type)
{

    //uint32_t      prevRootCbf;
    uint64_t      bestfullCost;
    uint32_t      fullLoopCandidateIndex;
    uint8_t       candidateIndex;

    uint64_t      y_full_distortion[DIST_CALC_TOTAL];
    uint32_t      count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];

    uint64_t      cbFullDistortion[DIST_CALC_TOTAL];
    uint64_t      crFullDistortion[DIST_CALC_TOTAL];

    uint64_t      y_coeff_bits;
    uint64_t        cb_coeff_bits = 0;
    uint64_t        cr_coeff_bits = 0;

    bestfullCost = 0xFFFFFFFFull;

    ModeDecisionCandidateBuffer_t         **candidateBufferPtrArrayBase = context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer_t         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[context_ptr->buffer_depth_index_start[0]]);
    ModeDecisionCandidateBuffer_t          *candidateBuffer;
    ModeDecisionCandidate_t                *candidate_ptr;

    //  if (context_ptr->blk_geom->origin_x == 16 && context_ptr->blk_geom->origin_y == 0 && context_ptr->blk_geom->bsize == BLOCK_8X8)
    //      printf("NOPPPP");

    for (fullLoopCandidateIndex = 0; fullLoopCandidateIndex < fullCandidateTotalCount; ++fullLoopCandidateIndex) {

        candidateIndex = context_ptr->best_candidate_index_array[fullLoopCandidateIndex];

        // initialize TU Split
        y_full_distortion[DIST_CALC_RESIDUAL] = 0;
        y_full_distortion[DIST_CALC_PREDICTION] = 0;
        y_coeff_bits = 0;

        // Set the Candidate Buffer
        candidateBuffer = candidate_buffer_ptr_array[candidateIndex];
        candidate_ptr = candidateBuffer->candidate_ptr;//this is the FastCandidateStruct

        candidate_ptr->full_distortion = 0;

        memset(candidate_ptr->eob[0], 0, sizeof(uint16_t));
        memset(candidate_ptr->eob[1], 0, sizeof(uint16_t));
        memset(candidate_ptr->eob[2], 0, sizeof(uint16_t));



        candidate_ptr->chroma_distortion = 0;
        candidate_ptr->chroma_distortion_inter_depth = 0;
#if !CHROMA_BLIND
        context_ptr->round_mv_to_integer = (candidate_ptr->type == INTER_MODE && candidate_ptr->merge_flag == EB_TRUE) ?
            EB_TRUE :
            EB_FALSE;

        context_ptr->round_mv_to_integer = picture_control_set_ptr->parent_pcs_ptr->use_subpel_flag ? EB_FALSE : context_ptr->round_mv_to_integer;
#endif
        // Set Skip Flag
        candidate_ptr->skip_flag = EB_FALSE;
#if INTERPOLATION_SEARCH_LEVELS
        if (picture_control_set_ptr->parent_pcs_ptr->interpolation_search_level == IT_SEARCH_FULL_LOOP) {
            context_ptr->skip_interpolation_search = 0;

            if (candidate_ptr->type != INTRA_MODE) {
#else
        if (candidate_ptr->prediction_is_ready_luma == EB_FALSE) {
#endif

#if CHROMA_BLIND
            ProductPredictionFunTable[candidate_ptr->type](
                context_ptr,
                picture_control_set_ptr,
                candidateBuffer,
                asm_type);
#else
            ProductPredictionFunTableCL[2][candidate_ptr->type](
                context_ptr,
                context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                picture_control_set_ptr,
                candidateBuffer,
                asm_type);
#endif                
            }
#if INTERPOLATION_SEARCH_LEVELS
        }
#endif

        //Y Residual
        ResidualKernel(
            &(inputPicturePtr->bufferY[inputOriginIndex]),
            inputPicturePtr->strideY,
            &(candidateBuffer->prediction_ptr->bufferY[cuOriginIndex]),
            candidateBuffer->prediction_ptr->strideY/* 64*/,
            &(((int16_t*)candidateBuffer->residual_ptr->bufferY)[cuOriginIndex]),
            candidateBuffer->residual_ptr->strideY,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight);

        //TOADD
        //Cb Residual
#if CHROMA_BLIND
        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
        if (context_ptr->blk_geom->has_uv) {
#endif

            ResidualKernel(
                &(inputPicturePtr->bufferCb[inputCbOriginIndex]),
                inputPicturePtr->strideCb,
                &(candidateBuffer->prediction_ptr->bufferCb[cuChromaOriginIndex]),
                candidateBuffer->prediction_ptr->strideCb,
                &(((int16_t*)candidateBuffer->residual_ptr->bufferCb)[cuChromaOriginIndex]),
                candidateBuffer->residual_ptr->strideCb,
                context_ptr->blk_geom->bwidth_uv,
                context_ptr->blk_geom->bheight_uv);

            //Cr Residual
            ResidualKernel(
                &(inputPicturePtr->bufferCr[inputCbOriginIndex]),
                inputPicturePtr->strideCr,
                &(candidateBuffer->prediction_ptr->bufferCr[cuChromaOriginIndex]),
                candidateBuffer->prediction_ptr->strideCr,
                &(((int16_t*)candidateBuffer->residual_ptr->bufferCr)[cuChromaOriginIndex]),
                candidateBuffer->residual_ptr->strideCr,
                context_ptr->blk_geom->bwidth_uv,
                context_ptr->blk_geom->bheight_uv);

        }
        // initialize luma CBF

        candidate_ptr->y_has_coeff = 0;
        candidate_ptr->u_has_coeff = 0;
        candidate_ptr->v_has_coeff = 0;

#if TX_SEARCH_LEVELS

        uint8_t  tx_search_skip_fag = picture_control_set_ptr->parent_pcs_ptr->tx_search_level == TX_SEARCH_FULL_LOOP ? get_skip_tx_search_flag(
            context_ptr->blk_geom->sq_size,
            ref_fast_cost,
            *candidateBuffer->fast_cost_ptr,
            picture_control_set_ptr->parent_pcs_ptr->tx_weight) : 1;

        if (!tx_search_skip_fag){
#else

#if TURN_OFF_TX_TYPE_SEARCH
#if ENCODER_MODE_CLEANUP
        if (picture_control_set_ptr->enc_mode <= ENC_M1) {
#endif
            if (context_ptr->blk_geom->sq_size < 128) //no tx search for 128x128 for now
#endif
#endif
                ProductFullLoopTxSearch(
                    candidateBuffer,
                    context_ptr,
                    picture_control_set_ptr);

            candidate_ptr->full_distortion = 0;


            memset(candidate_ptr->eob[0], 0, sizeof(uint16_t));


            //re-init
            candidate_ptr->y_has_coeff = 0;
#if ENCODER_MODE_CLEANUP
        }
#endif

        ProductFullLoop(
            candidateBuffer,
            context_ptr,
            picture_control_set_ptr,
            context_ptr->cu_ptr->qp,
            &(*count_non_zero_coeffs[0]),
            &y_coeff_bits,
            &y_full_distortion[0]);


        if (candidate_ptr->type == INTRA_MODE && candidateBuffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {

            // If mode is CFL:
            // 1: recon the Luma
            // 2: Form the pred_buf_q3
            // 3: Loop over alphas and find the best or choose DC
            // 4: Recalculate the residual for chroma
            CflPrediction(
                picture_control_set_ptr,
                candidateBuffer,
                sb_ptr,
                context_ptr,
                inputPicturePtr,
                inputCbOriginIndex,
                cuChromaOriginIndex,
                asm_type);

        }

        candidate_ptr->chroma_distortion_inter_depth = 0;
        candidate_ptr->chroma_distortion = 0;

        //CHROMA


        cbFullDistortion[DIST_CALC_RESIDUAL] = 0;
        crFullDistortion[DIST_CALC_RESIDUAL] = 0;
        cbFullDistortion[DIST_CALC_PREDICTION] = 0;
        crFullDistortion[DIST_CALC_PREDICTION] = 0;

        cb_coeff_bits = 0;
        cr_coeff_bits = 0;

        // FullLoop and TU search
        uint8_t cbQp = context_ptr->qp;
        uint8_t crQp = context_ptr->qp;

#if CHROMA_BLIND
        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
        if (context_ptr->blk_geom->has_uv) {
#endif
            FullLoop_R(
                sb_ptr,
                candidateBuffer,
                context_ptr,
                inputPicturePtr,
                picture_control_set_ptr,
                PICTURE_BUFFER_DESC_CHROMA_MASK,
                cbQp,
                crQp,
                &(*count_non_zero_coeffs[1]),
                &(*count_non_zero_coeffs[2]));





            CuFullDistortionFastTuMode_R(
                sb_ptr,
                candidateBuffer,
                context_ptr,
                candidate_ptr,
                picture_control_set_ptr,
                cbFullDistortion,
                crFullDistortion,
                count_non_zero_coeffs,
                COMPONENT_CHROMA,
                &cb_coeff_bits,
                &cr_coeff_bits,
                asm_type);
#if !CHROMA_BLIND
            candidate_ptr->block_has_coeff = (candidate_ptr->y_has_coeff | candidate_ptr->u_has_coeff | candidate_ptr->v_has_coeff) ? EB_TRUE : EB_FALSE;
#endif
        }

#if CHROMA_BLIND
        candidate_ptr->block_has_coeff = (candidate_ptr->y_has_coeff | candidate_ptr->u_has_coeff | candidate_ptr->v_has_coeff) ? EB_TRUE : EB_FALSE;
#endif


        //ALL PLANE
        Av1ProductFullCostFuncTable[candidate_ptr->type](
            picture_control_set_ptr,
            context_ptr,
            candidateBuffer,
            cu_ptr,
            y_full_distortion,
            cbFullDistortion,
            crFullDistortion,
            context_ptr->full_lambda,
            &y_coeff_bits,
            &cb_coeff_bits,
            &cr_coeff_bits,
            context_ptr->blk_geom->bsize);




        candidateBuffer->cb_distortion[DIST_CALC_RESIDUAL] = cbFullDistortion[DIST_CALC_RESIDUAL];
        candidateBuffer->cb_distortion[DIST_CALC_PREDICTION] = cbFullDistortion[DIST_CALC_PREDICTION];
        candidateBuffer->cb_coeff_bits = cb_coeff_bits;

        candidateBuffer->cr_distortion[DIST_CALC_RESIDUAL] = crFullDistortion[DIST_CALC_RESIDUAL];
        candidateBuffer->cr_distortion[DIST_CALC_PREDICTION] = crFullDistortion[DIST_CALC_PREDICTION];
        candidateBuffer->cr_coeff_bits = cr_coeff_bits;
        candidateBuffer->candidate_ptr->full_distortion = (uint32_t)(y_full_distortion[0]);
        candidateBuffer->candidate_ptr->luma_distortion = (uint32_t)(y_full_distortion[0]);


        candidateBuffer->y_coeff_bits = y_coeff_bits;
#if !CHROMA_BLIND
        if (picture_control_set_ptr->parent_pcs_ptr->chroma_mode == CHROMA_MODE_BEST)
        {
            candidateBuffer->y_coeff_bits = y_coeff_bits;
            candidateBuffer->y_full_distortion[DIST_CALC_RESIDUAL] = y_full_distortion[DIST_CALC_RESIDUAL];
            candidateBuffer->y_full_distortion[DIST_CALC_PREDICTION] = y_full_distortion[DIST_CALC_PREDICTION];
        }
#endif
#if 0 //AMIR_DEBUG
        //if (picture_control_set_ptr->parent_pcs_ptr->picture_number == 0 && /*context_ptr->cu_size > 16 &&*/ (cuStatsPtr)->origin_x >= 0 && (cuStatsPtr)->origin_x < 64 && (cuStatsPtr)->origin_y >= 0 && (cuStatsPtr)->origin_y < 64){
        printf("POC:%d\t(%d,%d)\t%d\t%d\t%d\t%d\t%lld\t%lld\t%lld\t%lld\t%lld\n",
            picture_control_set_ptr->parent_pcs_ptr->picture_number,
            sb_ptr->origin_x + (cuStatsPtr)->origin_x,
            sb_ptr->origin_y + (cuStatsPtr)->origin_y,
            (context_ptr)->cu_size,
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y],
            fullLoopCandidateIndex,
            candidate_ptr->pred_mode,
            *candidateBuffer->full_cost_ptr,
            candidateBuffer->y_coeff_bits,
            cb_coeff_bits + cr_coeff_bits,
            (uint64_t)candidateBuffer->candidate_ptr->full_distortion,
            candidateBuffer->cb_distortion[DIST_CALC_RESIDUAL] + candidateBuffer->cr_distortion[DIST_CALC_RESIDUAL]
        );

        //}
#endif

        candidate_ptr->full_distortion = (uint32_t)(y_full_distortion[0]);


#if SHUT_CBF_FL_SKIP
#if ENCODER_MODE_CLEANUP
        if(0)
#else
        if (/*sequence_control_set_ptr->static_config.tune == TUNE_VQ ||*/ picture_control_set_ptr->enc_mode > ENC_M1)
#endif
#endif
            if (picture_control_set_ptr->slice_type != I_SLICE) {
                if (candidate_ptr->type == INTER_MODE) {
                    if (*candidateBuffer->full_cost_ptr < bestfullCost) {
                        //prevRootCbf = candidate_ptr->yCbf;
                        bestfullCost = *candidateBuffer->full_cost_ptr;
                    }
                }
            }

    }//end for( full loop)
}

void move_cu_data(
    CodingUnit_t *src_cu,
    CodingUnit_t *dst_cu)
{


    //CHKN TransformUnit_t             transform_unit_array[TRANSFORM_UNIT_MAX_COUNT]; // 2-bytes * 21 = 42-bytes
    memcpy(dst_cu->transform_unit_array, src_cu->transform_unit_array, TRANSFORM_UNIT_MAX_COUNT * sizeof(TransformUnit_t));


    //CHKN PredictionUnit_t            prediction_unit_array[MAX_NUM_OF_PU_PER_CU];    // 35-bytes * 4 = 140 bytes
    memcpy(dst_cu->prediction_unit_array, src_cu->prediction_unit_array, MAX_NUM_OF_PU_PER_CU * sizeof(PredictionUnit_t));

    //CHKN     unsigned                    skip_flag_context : 2;
    //CHKN     unsigned                    prediction_mode_flag : 2;
    //CHKN     unsigned                    rootCbf : 1;
    //CHKN     unsigned                    split_flag_context : 2;
    //CHKN #if !ADD_DELTA_QP_SUPPORT
    //CHKN     unsigned                    qp : 6;
    //CHKN     unsigned                    ref_qp : 6;
    //CHKN
    //CHKN     signed                         delta_qp : 8; // can be signed 8bits
    //CHKN     signed                         org_delta_qp : 8;
    //CHKN #endif
    //CHKN
    //CHKN #if ADD_DELTA_QP_SUPPORT
    //CHKN     uint16_t                       qp;
    //CHKN     uint16_t                       ref_qp;
    //CHKN
    //CHKN     int16_t                          delta_qp; // can be signed 8bits
    //CHKN     int16_t                          org_delta_qp;
    //CHKN #endif

    dst_cu->skip_flag_context = src_cu->skip_flag_context;
    dst_cu->prediction_mode_flag = src_cu->prediction_mode_flag;
    dst_cu->block_has_coeff = src_cu->block_has_coeff;
    dst_cu->split_flag_context = src_cu->split_flag_context;
    dst_cu->qp = src_cu->qp;
    dst_cu->ref_qp = src_cu->ref_qp;
    dst_cu->delta_qp = src_cu->delta_qp;
    dst_cu->org_delta_qp = src_cu->org_delta_qp;




    //CHKN    // Coded Tree
    //CHKN    struct {
    //CHKN        unsigned                   leaf_index : 8;
    //CHKN        unsigned                   split_flag : 1;
    //CHKN        unsigned                   skip_flag : 1;
    //CHKN
    //CHKN    };

    dst_cu->leaf_index = src_cu->leaf_index;
    dst_cu->split_flag = src_cu->split_flag;
    dst_cu->skip_flag = src_cu->skip_flag;


    //CHKN    MacroBlockD*  av1xd;
    memcpy(dst_cu->av1xd, src_cu->av1xd, sizeof(MacroBlockD));


    // uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];

    //CHKN int16_t inter_mode_ctx[MODE_CTX_REF_FRAMES];
    memcpy(dst_cu->inter_mode_ctx, src_cu->inter_mode_ctx, MODE_CTX_REF_FRAMES * sizeof(int16_t));


    //CHKN IntMv ref_mvs[MODE_CTX_REF_FRAMES][MAX_MV_REF_CANDIDATES]; //used only for nonCompound modes.
    memcpy(dst_cu->ref_mvs, src_cu->ref_mvs, MODE_CTX_REF_FRAMES*MAX_MV_REF_CANDIDATES * sizeof(IntMv));

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

    dst_cu->skip_coeff_context = src_cu->skip_coeff_context;
    dst_cu->luma_txb_skip_context = src_cu->luma_txb_skip_context;
    dst_cu->luma_dc_sign_context = src_cu->luma_dc_sign_context;
    dst_cu->cb_txb_skip_context = src_cu->cb_txb_skip_context;
    dst_cu->cb_dc_sign_context = src_cu->cb_dc_sign_context;
    dst_cu->cr_txb_skip_context = src_cu->cr_txb_skip_context;
    dst_cu->cr_dc_sign_context = src_cu->cr_dc_sign_context;
    dst_cu->reference_mode_context = src_cu->reference_mode_context;
    dst_cu->compoud_reference_type_context = src_cu->compoud_reference_type_context;

    //CHKN int32_t                        quantized_dc[3];
    memcpy(dst_cu->quantized_dc, src_cu->quantized_dc, 3 * sizeof(int32_t));

    //CHKN uint32_t   is_inter_ctx;
    //CHKN uint32_t                     interp_filters;

    dst_cu->is_inter_ctx = src_cu->is_inter_ctx;
    dst_cu->interp_filters = src_cu->interp_filters;

    dst_cu->part = src_cu->part;
    dst_cu->shape = src_cu->shape;
    dst_cu->mds_idx = src_cu->mds_idx;
}



/*******************************************
* ModeDecision LCU
*   performs CL (LCU)
*******************************************/
EbBool allowed_ns_cu(
    ModeDecisionContext_t             *context_ptr,
    uint8_t                            is_complete_sb){
  
    EbBool  ret = 1;
    // Disable NSQ for non-complete LCU
    if (!is_complete_sb) {
        if (context_ptr->blk_geom->shape != PART_N) {
            ret = 0;
        }
    }
    return ret;
}

#if TX_SEARCH_LEVELS
void init_candidate_buffer(
    ModeDecisionCandidate_t        *candidate_ptr,
    uint32_t                        count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU])
{
    candidate_ptr->y_has_coeff = 0;
    candidate_ptr->u_has_coeff = 0;
    candidate_ptr->v_has_coeff = 0;

    candidate_ptr->full_distortion = 0;

    memset(candidate_ptr->eob[0], 0, sizeof(uint16_t)*MAX_TXB_COUNT);
    memset(count_non_zero_coeffs[0], 0, sizeof(uint32_t)*MAX_NUM_OF_TU_PER_CU);

    candidate_ptr->chroma_distortion = 0;
    candidate_ptr->chroma_distortion_inter_depth = 0;
    memset(candidate_ptr->eob[1], 0, sizeof(uint16_t)*MAX_TXB_COUNT);
    memset(count_non_zero_coeffs[1], 0, sizeof(uint32_t)*MAX_NUM_OF_TU_PER_CU);
    memset(candidate_ptr->eob[2], 0, sizeof(uint16_t)*MAX_TXB_COUNT);
    memset(count_non_zero_coeffs[2], 0, sizeof(uint32_t)*MAX_NUM_OF_TU_PER_CU);
}
void inter_depth_tx_search(
    PictureControlSet_t                      *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t            *candidateBuffer,
    CodingUnit_t                             *cu_ptr,
    ModeDecisionContext_t                    *context_ptr,
    EbPictureBufferDesc_t                    *inputPicturePtr,
    uint64_t                                  ref_fast_cost,
    EbAsm                                     asm_type)
{

    uint8_t  tx_search_skip_fag = picture_control_set_ptr->parent_pcs_ptr->tx_search_level == TX_SEARCH_INTER_DEPTH ? get_skip_tx_search_flag(
        context_ptr->blk_geom->sq_size,
        ref_fast_cost,
        *candidateBuffer->fast_cost_ptr,
        picture_control_set_ptr->parent_pcs_ptr->tx_weight) : 1;

    if (!tx_search_skip_fag) {

        uint64_t      y_full_distortion[DIST_CALC_TOTAL] = { 0 };
        uint32_t      count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU];

        uint64_t      cbFullDistortion[DIST_CALC_TOTAL];
        uint64_t      crFullDistortion[DIST_CALC_TOTAL];

        uint64_t      y_coeff_bits = 0;
        uint64_t      cb_coeff_bits = 0;
        uint64_t      cr_coeff_bits = 0;

        ModeDecisionCandidate_t                *candidate_ptr = candidateBuffer->candidate_ptr;

        init_candidate_buffer(
            candidate_ptr,
            count_non_zero_coeffs);


        ProductFullLoopTxSearch(
            candidateBuffer,
            context_ptr,
            picture_control_set_ptr
        );


        candidate_ptr->full_distortion = 0;


        memset(candidate_ptr->eob[0], 0, sizeof(uint16_t)*MAX_TXB_COUNT);

        //re-init
        candidate_ptr->y_has_coeff = 0;

        ProductFullLoop(
            candidateBuffer,
            context_ptr,
            picture_control_set_ptr,
            context_ptr->cu_ptr->qp,
            &(*count_non_zero_coeffs[0]),
            &y_coeff_bits,
            &y_full_distortion[0]);


        candidate_ptr->chroma_distortion_inter_depth = 0;
        candidate_ptr->chroma_distortion = 0;

        //CHROMA
        cbFullDistortion[DIST_CALC_RESIDUAL] = 0;
        crFullDistortion[DIST_CALC_RESIDUAL] = 0;
        cbFullDistortion[DIST_CALC_PREDICTION] = 0;
        crFullDistortion[DIST_CALC_PREDICTION] = 0;

        cb_coeff_bits = 0;
        cr_coeff_bits = 0;

        // FullLoop and TU search
        uint8_t cbQp = context_ptr->qp;
        uint8_t crQp = context_ptr->qp;
#if CHROMA_BLIND
        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
        if (context_ptr->blk_geom->has_uv) {
#endif
            FullLoop_R(
                context_ptr->sb_ptr,
                candidateBuffer,
                context_ptr,
                inputPicturePtr,
                picture_control_set_ptr,
                PICTURE_BUFFER_DESC_CHROMA_MASK,
                cbQp,
                crQp,
                &(*count_non_zero_coeffs[1]),
                &(*count_non_zero_coeffs[2]));

            CuFullDistortionFastTuMode_R(
                context_ptr->sb_ptr,
                candidateBuffer,
                context_ptr,
                candidate_ptr,
                picture_control_set_ptr,
                cbFullDistortion,
                crFullDistortion,
                count_non_zero_coeffs,
                COMPONENT_CHROMA,
                &cb_coeff_bits,
                &cr_coeff_bits,
                asm_type);

            candidate_ptr->block_has_coeff = (candidate_ptr->y_has_coeff | candidate_ptr->u_has_coeff | candidate_ptr->v_has_coeff) ? EB_TRUE : EB_FALSE;

        }

        Av1ProductFullCostFuncTable[candidate_ptr->type](
            picture_control_set_ptr,
            context_ptr,
            candidateBuffer,
            cu_ptr,
            y_full_distortion,
            cbFullDistortion,
            crFullDistortion,
            context_ptr->full_lambda,
            &y_coeff_bits,
            &cb_coeff_bits,
            &cr_coeff_bits,
            context_ptr->blk_geom->bsize);

        candidateBuffer->cb_distortion[DIST_CALC_RESIDUAL] = cbFullDistortion[DIST_CALC_RESIDUAL];
        candidateBuffer->cb_distortion[DIST_CALC_PREDICTION] = cbFullDistortion[DIST_CALC_PREDICTION];
        candidateBuffer->cb_coeff_bits = cb_coeff_bits;

        candidateBuffer->cr_distortion[DIST_CALC_RESIDUAL] = crFullDistortion[DIST_CALC_RESIDUAL];
        candidateBuffer->cr_distortion[DIST_CALC_PREDICTION] = crFullDistortion[DIST_CALC_PREDICTION];
        candidateBuffer->cr_coeff_bits = cr_coeff_bits;

        candidateBuffer->candidate_ptr->full_distortion = (uint32_t)(y_full_distortion[0]);
        candidateBuffer->candidate_ptr->luma_distortion = (uint32_t)(y_full_distortion[0]);


        candidateBuffer->y_coeff_bits = y_coeff_bits;
#if !CHROMA_BLIND
        if (picture_control_set_ptr->parent_pcs_ptr->chroma_mode == CHROMA_MODE_BEST)
        {
            candidateBuffer->y_coeff_bits = y_coeff_bits;
            candidateBuffer->y_full_distortion[DIST_CALC_RESIDUAL] = y_full_distortion[DIST_CALC_RESIDUAL];
            candidateBuffer->y_full_distortion[DIST_CALC_PREDICTION] = y_full_distortion[DIST_CALC_PREDICTION];
        }
#endif
        candidate_ptr->full_distortion = (uint32_t)(y_full_distortion[0]);
        //Update tx
        context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost = *(candidateBuffer->full_cost_ptr);
        context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost = (context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost - candidateBuffer->candidate_ptr->chroma_distortion) + candidateBuffer->candidate_ptr->chroma_distortion_inter_depth;

        if (candidate_ptr->type == INTRA_MODE)
            context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost_luma = candidateBuffer->full_cost_luma;


        context_ptr->md_ep_pipe_sb[cu_ptr->mds_idx].merge_cost = *candidateBuffer->full_cost_merge_ptr;
        context_ptr->md_ep_pipe_sb[cu_ptr->mds_idx].skip_cost = *candidateBuffer->full_cost_skip_ptr;


        if (candidate_ptr->type == INTER_MODE && candidate_ptr->merge_flag == EB_TRUE) {
            context_ptr->md_ep_pipe_sb[cu_ptr->leaf_index].chroma_distortion = candidateBuffer->candidate_ptr->chroma_distortion;
        }


        context_ptr->md_local_cu_unit[cu_ptr->mds_idx].full_distortion = candidateBuffer->candidate_ptr->full_distortion;

        context_ptr->md_local_cu_unit[cu_ptr->mds_idx].chroma_distortion = (uint32_t)candidateBuffer->candidate_ptr->chroma_distortion;
        context_ptr->md_local_cu_unit[cu_ptr->mds_idx].chroma_distortion_inter_depth = (uint32_t)candidateBuffer->candidate_ptr->chroma_distortion_inter_depth;

        //cu_ptr->prediction_mode_flag = candidate_ptr->type;
        cu_ptr->skip_flag = candidate_ptr->skip_flag; // note, the skip flag is re-checked in the ENCDEC process
        cu_ptr->block_has_coeff = ((candidate_ptr->block_has_coeff) > 0) ? EB_TRUE : EB_FALSE;
        cu_ptr->quantized_dc[0] = candidateBuffer->candidate_ptr->quantized_dc[0];
        cu_ptr->quantized_dc[1] = candidateBuffer->candidate_ptr->quantized_dc[1];
        cu_ptr->quantized_dc[2] = candidateBuffer->candidate_ptr->quantized_dc[2];

        context_ptr->md_local_cu_unit[cu_ptr->mds_idx].count_non_zero_coeffs = candidate_ptr->count_non_zero_coeffs;

        TransformUnit_t        *txb_ptr;
        uint32_t                txb_itr;
        uint32_t                tu_index;
        uint32_t                tuTotalCount;

        tuTotalCount = context_ptr->blk_geom->txb_count;
        tu_index = 0;
        txb_itr = 0;
        

#if NO_ENCDEC
        int32_t txb_1d_offset = 0, txb_1d_offset_uv = 0;

        cu_ptr->block_has_coeff = 0;
#endif


        // Set TU
        do {

            txb_ptr = &cu_ptr->transform_unit_array[tu_index];

            txb_ptr->split_flag = EB_FALSE;
            txb_ptr->y_has_coeff = (EbBool)(((candidate_ptr->y_has_coeff)  & (1 << tu_index)) > 0);
            txb_ptr->u_has_coeff = (EbBool)(((candidate_ptr->u_has_coeff) & (1 << (tu_index))) > 0);
            txb_ptr->v_has_coeff = (EbBool)(((candidate_ptr->v_has_coeff) & (1 << (tu_index))) > 0);
            txb_ptr->transform_type[PLANE_TYPE_Y] = candidate_ptr->transform_type[PLANE_TYPE_Y];
            txb_ptr->transform_type[PLANE_TYPE_UV] = candidate_ptr->transform_type[PLANE_TYPE_UV];

#if NO_ENCDEC

            if (context_ptr->blk_geom->has_uv) {
                cu_ptr->block_has_coeff |= txb_ptr->y_has_coeff;
                cu_ptr->block_has_coeff |= txb_ptr->u_has_coeff;
                cu_ptr->block_has_coeff |= txb_ptr->v_has_coeff;
            }
            else {
                cu_ptr->block_has_coeff |= txb_ptr->y_has_coeff;
            }


            cu_ptr->cand_buff_index = lowestCostIndex;

            cu_ptr->skip_flag = 0;   //SKIP is turned OFF for this case!!
            txb_ptr->nz_coef_count[0] = candidate_ptr->eob[0][tu_index];
            txb_ptr->nz_coef_count[1] = candidate_ptr->eob[1][tu_index];
            txb_ptr->nz_coef_count[2] = candidate_ptr->eob[2][tu_index];

            if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0) {
                cu_ptr->predmv[0].as_mv.col = candidate_ptr->motion_vector_pred_x[REF_LIST_0];
                cu_ptr->predmv[0].as_mv.row = candidate_ptr->motion_vector_pred_y[REF_LIST_0];
            }
            else if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1) {
                cu_ptr->predmv[0].as_mv.col = candidate_ptr->motion_vector_pred_x[REF_LIST_1];
                cu_ptr->predmv[0].as_mv.row = candidate_ptr->motion_vector_pred_y[REF_LIST_1];
            }
            else if (pu_ptr->inter_pred_direction_index == BI_PRED) {
                cu_ptr->predmv[0].as_mv.col = candidate_ptr->motion_vector_pred_x[REF_LIST_0];
                cu_ptr->predmv[0].as_mv.row = candidate_ptr->motion_vector_pred_y[REF_LIST_0];
                cu_ptr->predmv[1].as_mv.col = candidate_ptr->motion_vector_pred_x[REF_LIST_1];
                cu_ptr->predmv[1].as_mv.row = candidate_ptr->motion_vector_pred_y[REF_LIST_1];
            }
#endif
#if NO_ENCDEC
            //copy coeff
            {

                uint32_t  bwidth = context_ptr->blk_geom->tx_width[txb_itr] < 64 ? context_ptr->blk_geom->tx_width[txb_itr] : 32;
                uint32_t  bheight = context_ptr->blk_geom->tx_height[txb_itr] < 64 ? context_ptr->blk_geom->tx_height[txb_itr] : 32;

                int32_t* srcPtr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residualQuantCoeffPtr->bufferY)[txb_1d_offset]);
                int32_t* dstPtr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->bufferY)[txb_1d_offset]);

                uint32_t j;

                for (j = 0; j < bheight; j++)
                {
                    memcpy(dstPtr + j * bwidth, srcPtr + j * bwidth, bwidth * sizeof(int32_t));
                }

                if (context_ptr->blk_geom->has_uv)
                {
                    // Cb
                    bwidth = context_ptr->blk_geom->tx_width_uv[txb_itr];
                    bheight = context_ptr->blk_geom->tx_height_uv[txb_itr];

                    srcPtr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residualQuantCoeffPtr->bufferCb)[txb_1d_offset_uv]);
                    dstPtr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->bufferCb)[txb_1d_offset_uv]);

                    for (j = 0; j < bheight; j++)
                    {
                        memcpy(dstPtr + j * bwidth, srcPtr + j * bwidth, bwidth * sizeof(int32_t));
                    }

                    srcPtr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residualQuantCoeffPtr->bufferCr)[txb_1d_offset_uv]);
                    dstPtr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->bufferCr)[txb_1d_offset_uv]);

                    for (j = 0; j < bheight; j++)
                    {
                        memcpy(dstPtr + j * bwidth, srcPtr + j * bwidth, bwidth * sizeof(int32_t));
                    }
                }

                txb_1d_offset += context_ptr->blk_geom->tx_width[txb_itr] * context_ptr->blk_geom->tx_height[txb_itr];
                if (context_ptr->blk_geom->has_uv)
                    txb_1d_offset_uv += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];

            }

#endif

            ++tu_index;
            ++txb_itr;

        } while (txb_itr < tuTotalCount);
    }
}
#endif
void md_encode_block(
    SequenceControlSet_t             *sequence_control_set_ptr,
    PictureControlSet_t              *picture_control_set_ptr,
    ModeDecisionContext_t            *context_ptr,
    SsMeContext_t                    *ss_mecontext,
    uint32_t                          leaf_index,
    uint32_t                          lcuAddr,
    ModeDecisionCandidateBuffer_t    *bestCandidateBuffers[5])
{

    ModeDecisionCandidateBuffer_t         **candidateBufferPtrArrayBase = context_ptr->candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer_t         **candidate_buffer_ptr_array;
    const BlockGeom                          *blk_geom = context_ptr->blk_geom;

    uint32_t                                  buffer_total_count;
    ModeDecisionCandidateBuffer_t            *candidateBuffer;
    ModeDecisionCandidate_t                  *fast_candidate_array = context_ptr->fast_candidate_array;
    uint8_t                                   candidateIndex;
    uint32_t                                  fastCandidateTotalCount;
    uint32_t                                  fullCandidateTotalCount;
    uint32_t                                  maxBuffers;
    uint32_t                                  secondFastCostSearchCandidateTotalCount;
    EbAsm                                     asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
    uint32_t                                  best_intra_mode = EB_INTRA_MODE_INVALID;

    EbPictureBufferDesc_t                    *inputPicturePtr = picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    const uint32_t                            inputOriginIndex = (context_ptr->cu_origin_y + inputPicturePtr->origin_y) * inputPicturePtr->strideY + (context_ptr->cu_origin_x + inputPicturePtr->origin_x);

    const uint32_t inputCbOriginIndex = ((context_ptr->round_origin_y >> 1) + (inputPicturePtr->origin_y >> 1)) * inputPicturePtr->strideCb + ((context_ptr->round_origin_x >> 1) + (inputPicturePtr->origin_x >> 1));
    const uint32_t cuOriginIndex = blk_geom->origin_x + blk_geom->origin_y * SB_STRIDE_Y;
    const uint32_t cuChromaOriginIndex = ROUND_UV(blk_geom->origin_x) / 2 + ROUND_UV(blk_geom->origin_y) / 2 * SB_STRIDE_UV;
    CodingUnit_t *  cu_ptr = context_ptr->cu_ptr;
    candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[context_ptr->buffer_depth_index_start[0]]);

    if (allowed_ns_cu(
#if DISABLE_NSQ_FOR_NON_REF || DISABLE_NSQ
#if ENCODER_MODE_CLEANUP
        context_ptr, sequence_control_set_ptr->sb_geom[lcuAddr].is_complete_sb))
#else
        context_ptr, sequence_control_set_ptr->sb_geom[lcuAddr].is_complete_sb))
#endif
#else
        context_ptr, sequence_control_set_ptr->sb_geom[lcuAddr].is_complete_sb))
#endif
    {
        // Set PF Mode - should be done per TU (and not per CU) to avoid the correction
        ProductDerivePartialFrequencyN2Flag(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            context_ptr);

        ProductCodingLoopInitFastLoop(
            context_ptr,
            context_ptr->skip_coeff_neighbor_array,
            context_ptr->luma_dc_sign_level_coeff_neighbor_array,
            context_ptr->cb_dc_sign_level_coeff_neighbor_array,
            context_ptr->cr_dc_sign_level_coeff_neighbor_array,
            context_ptr->inter_pred_dir_neighbor_array,
            context_ptr->ref_frame_type_neighbor_array,

            context_ptr->intra_luma_mode_neighbor_array,
            context_ptr->skip_flag_neighbor_array,
            context_ptr->mode_type_neighbor_array,
            context_ptr->leaf_depth_neighbor_array,
            context_ptr->leaf_partition_neighbor_array);

        set_nfl(
            context_ptr,
            picture_control_set_ptr);

        ProductGenerateMdCandidatesCu(
            context_ptr->sb_ptr,
            context_ptr,
            ss_mecontext,
            leaf_index,
            lcuAddr,
            &buffer_total_count,
            &fastCandidateTotalCount,
            (void*)context_ptr->inter_prediction_context,
            picture_control_set_ptr);

        //if we want to recon N candidate, we would need N+1 buffers
        maxBuffers = MIN((buffer_total_count + 1), context_ptr->buffer_depth_index_width[0]);

        ProductPerformFastLoop(
            picture_control_set_ptr,
            context_ptr->sb_ptr,
            context_ptr,
            candidateBufferPtrArrayBase,
            fast_candidate_array,
            fastCandidateTotalCount,
            inputPicturePtr,
            inputOriginIndex,
            inputCbOriginIndex,
            inputCbOriginIndex,
            cu_ptr,
            cuOriginIndex,
            cuChromaOriginIndex,
            maxBuffers,
            &secondFastCostSearchCandidateTotalCount,
            asm_type);

        // Make sure buffer_total_count is not larger than the number of fast modes
        buffer_total_count = MIN(secondFastCostSearchCandidateTotalCount, buffer_total_count);

        // PreModeDecision
        // -Input is the buffers
        // -Output is list of buffers for full reconstruction
        uint8_t  disable_merge_index = 0;

#if TX_SEARCH_LEVELS
        uint64_t ref_fast_cost = MAX_MODE_COST;
#endif
        PreModeDecision(
            cu_ptr,
            (secondFastCostSearchCandidateTotalCount == buffer_total_count) ? buffer_total_count : maxBuffers,
            candidate_buffer_ptr_array,
            &fullCandidateTotalCount,
            context_ptr->best_candidate_index_array,
            &disable_merge_index,
#if TX_SEARCH_LEVELS
            &ref_fast_cost,
#endif
            (EbBool)(secondFastCostSearchCandidateTotalCount == buffer_total_count)); // The fast loop bug fix is now added to 4K only


        AV1PerformFullLoop(
            picture_control_set_ptr,
            context_ptr->sb_ptr,
            cu_ptr,
            context_ptr,
            inputPicturePtr,
            inputOriginIndex,
            inputCbOriginIndex,
            cuOriginIndex,
            cuChromaOriginIndex,
            MIN(fullCandidateTotalCount, buffer_total_count),
#if TX_SEARCH_LEVELS
            ref_fast_cost,
#endif
            asm_type); // fullCandidateTotalCount to number of buffers to process

        // Full Mode Decision (choose the best mode)
        candidateIndex = product_full_mode_decision(
            context_ptr,
            cu_ptr,
            context_ptr->blk_geom->bwidth,
            context_ptr->blk_geom->bheight,
            candidate_buffer_ptr_array,
            fullCandidateTotalCount,
            context_ptr->best_candidate_index_array,
            &best_intra_mode);

        candidateBuffer = candidate_buffer_ptr_array[candidateIndex];

        bestCandidateBuffers[0] = candidateBuffer;

#if INTERPOLATION_SEARCH_LEVELS
        if (picture_control_set_ptr->parent_pcs_ptr->interpolation_search_level == IT_SEARCH_INTER_DEPTH) {

            if (candidateBuffer->candidate_ptr->type != INTRA_MODE && candidateBuffer->candidate_ptr->motion_mode == SIMPLE_TRANSLATION) {

                context_ptr->skip_interpolation_search = 0;
#if CHROMA_BLIND
                ProductPredictionFunTable[candidateBuffer->candidate_ptr->type](
                    context_ptr,
                    picture_control_set_ptr,
                    candidateBuffer,
                    asm_type);
#else
                ProductPredictionFunTableCL[2][candidateBuffer->candidate_ptr->type](
                    context_ptr,
                    context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                    picture_control_set_ptr,
                    candidateBuffer,
                    asm_type);
#endif
                cu_ptr->interp_filters = candidateBuffer->candidate_ptr->interp_filters;
            }
        }
#endif

#if TX_SEARCH_LEVELS
        inter_depth_tx_search(
            picture_control_set_ptr,
            candidateBuffer,
            cu_ptr,
            context_ptr,
            inputPicturePtr,
            ref_fast_cost,
            asm_type);
#endif

#if NSQ_SEARCH_LEVELS
        uint8_t sq_index = LOG2F(context_ptr->blk_geom->sq_size) - 2;
        if (context_ptr->blk_geom->shape == PART_N) {

            context_ptr->parent_sq_type[sq_index] = candidateBuffer->candidate_ptr->type;

            context_ptr->parent_sq_has_coeff[sq_index] = (candidateBuffer->candidate_ptr->y_has_coeff ||
                candidateBuffer->candidate_ptr->u_has_coeff ||
                candidateBuffer->candidate_ptr->v_has_coeff) ? 1 : 0;

            context_ptr->parent_sq_pred_mode[sq_index] = candidateBuffer->candidate_ptr->pred_mode;
        }
#endif

        AV1PerformInverseTransformRecon(
            picture_control_set_ptr,
            context_ptr,
            candidateBuffer,
            cu_ptr,
            context_ptr->blk_geom,
            asm_type);

        //copy neigh recon data in cu_ptr
        {
            uint32_t j;
            EbPictureBufferDesc_t  *reconPtr = candidateBuffer->reconPtr;
            uint32_t recLumaOffset = context_ptr->blk_geom->origin_x + context_ptr->blk_geom->origin_y * reconPtr->strideY;

            uint32_t recCbOffset = ((((context_ptr->blk_geom->origin_x >> 3) << 3) + ((context_ptr->blk_geom->origin_y >> 3) << 3) * candidateBuffer->reconPtr->strideCb) >> 1);
            uint32_t recCrOffset = ((((context_ptr->blk_geom->origin_x >> 3) << 3) + ((context_ptr->blk_geom->origin_y >> 3) << 3) * candidateBuffer->reconPtr->strideCr) >> 1);

            memcpy(cu_ptr->neigh_top_recon[0], reconPtr->bufferY + recLumaOffset + (context_ptr->blk_geom->bheight - 1)*reconPtr->strideY, context_ptr->blk_geom->bwidth);
#if CHROMA_BLIND
            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0)
#else
            if (context_ptr->blk_geom->has_uv)
#endif
            {
                memcpy(cu_ptr->neigh_top_recon[1], reconPtr->bufferCb + recCbOffset + (context_ptr->blk_geom->bheight_uv - 1)*reconPtr->strideCb, context_ptr->blk_geom->bwidth_uv);
                memcpy(cu_ptr->neigh_top_recon[2], reconPtr->bufferCr + recCrOffset + (context_ptr->blk_geom->bheight_uv - 1)*reconPtr->strideCr, context_ptr->blk_geom->bwidth_uv);
            }

            for (j = 0; j < context_ptr->blk_geom->bheight; ++j)
                cu_ptr->neigh_left_recon[0][j] = reconPtr->bufferY[recLumaOffset + context_ptr->blk_geom->bwidth - 1 + j * reconPtr->strideY];
#if CHROMA_BLIND
            if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
            if (context_ptr->blk_geom->has_uv)
            {
#endif
                for (j = 0; j < context_ptr->blk_geom->bheight_uv; ++j) {
                    cu_ptr->neigh_left_recon[1][j] = reconPtr->bufferCb[recCbOffset + context_ptr->blk_geom->bwidth_uv - 1 + j * reconPtr->strideCb];
                    cu_ptr->neigh_left_recon[2][j] = reconPtr->bufferCr[recCrOffset + context_ptr->blk_geom->bwidth_uv - 1 + j * reconPtr->strideCr];
                }
            }
        }


#if NO_ENCDEC
        //copy recon
        {

            uint32_t  tuOriginIndex = context_ptr->blk_geom->origin_x + (context_ptr->blk_geom->origin_y * 128);
            uint32_t  bwidth = context_ptr->blk_geom->bwidth;
            uint32_t  bheight = context_ptr->blk_geom->bheight;

            uint8_t* srcPtr = &(((uint8_t*)candidateBuffer->reconPtr->bufferY)[tuOriginIndex]);
            uint8_t* dstPtr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->bufferY)[0]);

            uint32_t j;
            for (j = 0; j < bheight; j++)
            {
                memcpy(dstPtr + j * 128, srcPtr + j * 128, bwidth * sizeof(uint8_t));
            }

            // Cb
            if (context_ptr->blk_geom->has_uv)
            {

                uint32_t tuOriginIndex = ((((context_ptr->blk_geom->origin_x >> 3) << 3) + ((context_ptr->blk_geom->origin_y >> 3) << 3) * candidateBuffer->reconPtr->strideCb) >> 1);

                bwidth = context_ptr->blk_geom->bwidth_uv;
                bheight = context_ptr->blk_geom->bheight_uv;

                srcPtr = &(((uint8_t*)candidateBuffer->reconPtr->bufferCb)[tuOriginIndex]);
                dstPtr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->bufferCb)[0]);

                for (j = 0; j < bheight; j++)
                {
                    memcpy(dstPtr + j * 64, srcPtr + j * 64, bwidth * sizeof(uint8_t));
                }

                // Cr

                srcPtr = &(((uint8_t*)candidateBuffer->reconPtr->bufferCr)[tuOriginIndex]);
                dstPtr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->bufferCr)[0]);

                for (j = 0; j < bheight; j++)
                {
                    memcpy(dstPtr + j * 64, srcPtr + j * 64, bwidth * sizeof(uint8_t));
                }

            }

        }
#endif


    }
    else
    {
        context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost = MAX_MODE_COST;
        cu_ptr->prediction_unit_array->ref_frame_type = 0;
    }
}


EB_EXTERN EbErrorType mode_decision_sb(
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureControlSet_t                 *picture_control_set_ptr,
    const MdcLcuData_t * const           mdcResultTbPtr,
    LargestCodingUnit_t                 *sb_ptr,
    uint16_t                             sb_origin_x,
    uint16_t                             sb_origin_y,
    uint32_t                             lcuAddr,
    SsMeContext_t                       *ss_mecontext,
    ModeDecisionContext_t               *context_ptr)
{
    EbErrorType                          return_error = EB_ErrorNone;

    uint32_t                             cuIdx;
    ModeDecisionCandidateBuffer_t       *bestCandidateBuffers[5];

    // CTB merge
    uint32_t                               lastCuIndex;

    // Pre Intra Search
    EbAsm                                  asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
    const uint32_t                         sb_height = MIN(BLOCK_SIZE_64, (uint32_t)(sequence_control_set_ptr->luma_height - sb_origin_y));

    uint32_t                               leaf_count = mdcResultTbPtr->leaf_count;
    const EbMdcLeafData_t *const           leaf_data_array = mdcResultTbPtr->leaf_data_array;
    UNUSED(sb_height);
    UNUSED(asm_type);
    UNUSED(lastCuIndex);

    context_ptr->sb_ptr = sb_ptr;
    context_ptr->group_of8x8_blocks_count = 0;
    context_ptr->group_of16x16_blocks_count = 0;

    ProductConfigureChroma(
        picture_control_set_ptr,
        context_ptr,
        sb_ptr);
    Initialize_cu_data_structure(
        context_ptr,
        sequence_control_set_ptr,
        sb_ptr,
        mdcResultTbPtr);

    // Mode Decision Neighbor Arrays
    context_ptr->intra_luma_mode_neighbor_array = picture_control_set_ptr->md_intra_luma_mode_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->intra_chroma_mode_neighbor_array = picture_control_set_ptr->md_intra_chroma_mode_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->mv_neighbor_array = picture_control_set_ptr->md_mv_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->skip_flag_neighbor_array = picture_control_set_ptr->md_skip_flag_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->mode_type_neighbor_array = picture_control_set_ptr->md_mode_type_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->leaf_depth_neighbor_array = picture_control_set_ptr->md_leaf_depth_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->leaf_partition_neighbor_array = picture_control_set_ptr->mdleaf_partition_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->luma_recon_neighbor_array = picture_control_set_ptr->md_luma_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->cb_recon_neighbor_array = picture_control_set_ptr->md_cb_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->cr_recon_neighbor_array = picture_control_set_ptr->md_cr_recon_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];

    context_ptr->skip_coeff_neighbor_array = picture_control_set_ptr->md_skip_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->luma_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->md_luma_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->cb_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->md_cb_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->cr_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->md_cr_dc_sign_level_coeff_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->inter_pred_dir_neighbor_array = picture_control_set_ptr->md_inter_pred_dir_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->ref_frame_type_neighbor_array = picture_control_set_ptr->md_ref_frame_type_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];
    context_ptr->interpolation_type_neighbor_array = picture_control_set_ptr->md_interpolation_type_neighbor_array[MD_NEIGHBOR_ARRAY_INDEX];

    //CU Loop
    cuIdx = 0;  //index over mdc array

    uint32_t blk_idx_mds = 0;

    EbBool all_d1_blocks_done = 0;
    uint32_t  d1_blocks_accumlated = 0;
    UNUSED(all_d1_blocks_done);

    do
    {

        blk_idx_mds = leaf_data_array[cuIdx].mds_idx;

        const BlockGeom * blk_geom = context_ptr->blk_geom = Get_blk_geom_mds(blk_idx_mds);
        CodingUnit_t *  cu_ptr = context_ptr->cu_ptr = &context_ptr->md_cu_arr_nsq[blk_idx_mds];

        context_ptr->cu_size_log2 = blk_geom->bwidth_log2;
        context_ptr->cu_origin_x = sb_origin_x + blk_geom->origin_x;
        context_ptr->cu_origin_y = sb_origin_y + blk_geom->origin_y;

        const EbMdcLeafData_t * const leafDataPtr = &mdcResultTbPtr->leaf_data_array[cuIdx];
        context_ptr->sb_sz = BLOCK_SIZE_64;
        context_ptr->round_origin_x = ((context_ptr->cu_origin_x >> 3) << 3);
        context_ptr->round_origin_y = ((context_ptr->cu_origin_y >> 3) << 3);
        context_ptr->sb_origin_x = sb_origin_x;
        context_ptr->sb_origin_y = sb_origin_y;
        context_ptr->md_local_cu_unit[blk_idx_mds].tested_cu_flag = EB_TRUE;

        cu_ptr->mds_idx = blk_idx_mds;
#if FIX_INTER_DEPTH
        context_ptr->md_cu_arr_nsq[blk_idx_mds].mdc_split_flag = (uint16_t)leafDataPtr->split_flag;

#endif
        cu_ptr->split_flag = (uint16_t)leafDataPtr->split_flag; //mdc indicates smallest or non valid CUs with split flag=
        cu_ptr->qp = context_ptr->qp;
        cu_ptr->best_d1_blk = blk_idx_mds;
#if INJECT_ONLY_SQ

            if (leafDataPtr->tot_d1_blocks != 1)
            {
                if (blk_geom->shape == PART_N)
                    copy_neighbour_arrays(      //save a clean neigh in [1], encode uses [0], reload the clean in [0] after done last ns block in a partition
                        picture_control_set_ptr,
                        context_ptr,
                        0, 1,
                        blk_idx_mds,
                        sb_origin_x,
                        sb_origin_y);
            }
#else
        if (blk_geom->shape == PART_N)
            copy_neighbour_arrays(      //save a clean neigh in [1], encode uses [0], reload the clean in [0] after done last ns block in a partition
                picture_control_set_ptr,
                context_ptr,
                0, 1,
                blk_idx_mds,
                sb_origin_x,
                sb_origin_y);
#endif

        md_encode_block(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            context_ptr,
            ss_mecontext,
            0xFFFFFFFF,
            lcuAddr,
            bestCandidateBuffers);

        if (blk_geom->nsi + 1 == blk_geom->totns)
            d1_non_square_block_decision(context_ptr);

        if (blk_geom->shape != PART_N) {
            if (blk_geom->nsi + 1 < blk_geom->totns)
                md_update_all_neighbour_arrays(
                    picture_control_set_ptr,
                    context_ptr,
                    blk_idx_mds,
                    sb_origin_x,
                    sb_origin_y);
            else
                copy_neighbour_arrays(      //restore [1] in [0] after done last ns block
                    picture_control_set_ptr,
                    context_ptr,
                    1, 0,
                    blk_geom->sqi_mds,
                    sb_origin_x,
                    sb_origin_y);
        }



        d1_blocks_accumlated = blk_geom->shape == PART_N ? 1 : d1_blocks_accumlated + 1;

        if (d1_blocks_accumlated == leafDataPtr->tot_d1_blocks)
        {

            uint32_t  lastCuIndex_mds = d2_inter_depth_block_decision(
                context_ptr,
                blk_geom->sqi_mds,//input is parent square
                sb_ptr,
                lcuAddr,
                sb_origin_x,
                sb_origin_y,
                context_ptr->full_lambda,
                context_ptr->md_rate_estimation_ptr,
                picture_control_set_ptr);


            if (context_ptr->md_cu_arr_nsq[lastCuIndex_mds].split_flag == EB_FALSE)
            {

                md_update_all_neighbour_arrays_multiple(
                    picture_control_set_ptr,
                    context_ptr,
                    context_ptr->md_cu_arr_nsq[lastCuIndex_mds].best_d1_blk,
                    sb_origin_x,
                    sb_origin_y);

            }


        }

        cuIdx++;

    } while (cuIdx < leaf_count);// End of CU loop




    return return_error;
}



/*******************************************
* Compute4x4SAD_Default
*   Unoptimized 4x4 SAD
*******************************************/
uint32_t Compute4x4SAD_Kernel(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  refStride,                      // input parameter, reference stride
    uint32_t  height,                         // input parameter, block height (M)
    uint32_t  width)                          // input parameter, block width (N)
{
    uint32_t rowNumberInBlock4x4;
    uint32_t sadBlock4x4 = 0;

    for (rowNumberInBlock4x4 = 0; rowNumberInBlock4x4 < 4; ++rowNumberInBlock4x4) {
        sadBlock4x4 += EB_ABS_DIFF(src[0x00], ref[0x00]);
        sadBlock4x4 += EB_ABS_DIFF(src[0x01], ref[0x01]);
        sadBlock4x4 += EB_ABS_DIFF(src[0x02], ref[0x02]);
        sadBlock4x4 += EB_ABS_DIFF(src[0x03], ref[0x03]);

        src += src_stride;
        ref += refStride;
    }
    (void)height;
    (void)width;
    return sadBlock4x4;
}

static EB_SADKERNELNxM_TYPE FUNC_TABLE compute4x4SAD_funcPtrArray[ASM_TYPE_TOTAL] =// [C_DEFAULT/ASM]
{
    // C_DEFAULT
    Compute4x4SAD_Kernel,
#if INTRINSIC_OPT_2
    // SSE2
    Compute4xMSad_AVX2_INTRIN,
#else
    // SSE2
    Compute4x4SAD_Kernel,
#endif


};


static uint32_t tab4x4[256] = {
    0, 1, 4, 5, 16, 17, 20, 21, 64, 65, 68, 69, 80, 81, 84, 85,
    2, 3, 6, 7, 18, 19, 22, 23, 66, 67, 70, 71, 82, 83, 86, 87,
    8, 9, 12, 13, 24, 25, 28, 29, 72, 73, 76, 77, 88, 89, 92, 93,
    10, 11, 14, 15, 26, 27, 30, 31, 74, 75, 78, 79, 90, 91, 94, 95,
    32, 33, 36, 37, 48, 49, 52, 53, 96, 97, 100, 101, 112, 113, 116, 117,
    34, 35, 38, 39, 50, 51, 54, 55, 98, 99, 102, 103, 114, 115, 118, 119,
    40, 41, 44, 45, 56, 57, 60, 61, 104, 105, 108, 109, 120, 121, 124, 125,
    42, 43, 46, 47, 58, 59, 62, 63, 106, 107, 110, 111, 122, 123, 126, 127,
    128, 129, 132, 133, 144, 145, 148, 149, 192, 193, 196, 197, 208, 209, 212, 213,
    130, 131, 134, 135, 146, 147, 150, 151, 194, 195, 198, 199, 210, 211, 214, 215,
    136, 137, 140, 141, 152, 153, 156, 157, 200, 201, 204, 205, 216, 217, 220, 221,
    138, 139, 142, 143, 154, 155, 158, 159, 202, 203, 206, 207, 218, 219, 222, 223,
    160, 161, 164, 165, 176, 177, 180, 181, 224, 225, 228, 229, 240, 241, 244, 245,
    162, 163, 166, 167, 178, 179, 182, 183, 226, 227, 230, 231, 242, 243, 246, 247,
    168, 169, 172, 173, 184, 185, 188, 189, 232, 233, 236, 237, 248, 249, 252, 253,
    170, 171, 174, 175, 186, 187, 190, 191, 234, 235, 238, 239, 250, 251, 254, 255,

};

static uint32_t tab8x4[128] = {
    0, 2, 8, 10, 32, 34, 40, 42,
    1, 3, 9, 11, 33, 35, 41, 43,
    4, 6, 12, 14, 36, 38, 44, 46,
    5, 7, 13, 15, 37, 39, 45, 47,
    16, 18, 24, 26, 48, 50, 56, 58,
    17, 19, 25, 27, 49, 51, 57, 59,
    20, 22, 28, 30, 52, 54, 60, 62,
    21, 23, 29, 31, 53, 55, 61, 63,
    64, 66, 72, 74, 96, 98, 104, 106,
    65, 67, 73, 75, 97, 99, 105, 107,
    68, 70, 76, 78, 100, 102, 108, 110,
    69, 71, 77, 79, 101, 103, 109, 111,
    80, 82, 88, 90, 112, 114, 120, 122,
    81, 83, 89, 91, 113, 115, 121, 123,
    84, 86, 92, 94, 116, 118, 124, 126,
    85, 87, 93, 95, 117, 119, 125, 127
};

static uint32_t tab4x8[128] = {
    0, 1, 2, 3, 8, 9, 10, 11, 32, 33, 34, 35, 40, 41, 42, 43,
    4, 5, 6, 7, 12, 13, 14, 15, 36, 37, 38, 39, 44, 45, 46, 47,
    16, 17, 18, 19, 24, 25, 26, 27, 48, 49, 50, 51, 56, 57, 58, 59,
    20, 21, 22, 23, 28, 29, 30, 31, 52, 53, 54, 55, 60, 61, 62, 63,
    64, 65, 66, 67, 72, 73, 74, 75, 96, 97, 98, 99, 104, 105, 106, 107,
    68, 69, 70, 71, 76, 77, 78, 79, 100, 101, 102, 103, 108, 109, 110, 111,
    80, 81, 82, 83, 88, 89, 90, 91, 112, 113, 114, 115, 120, 121, 122, 123,
    84, 85, 86, 87, 92, 93, 94, 95, 116, 117, 118, 119, 124, 125, 126, 127
};


static uint32_t tab16x4[64] = {
     0    ,        4,        16,            20,
     1    ,        5,        17,            21,
     2    ,        6,        18,            22,
     3    ,        7,        19,            23,
     8    ,        12,        24,            28,
     9    ,        13,        25,            29,
     10,        14,        26,            30,
     11,        15,        27,            31,
     32,        36,        48,            52,
     33,        37,        49,            53,
     34,        38,        50,            54,
     35,        39,        51,            55,
     40,        44,        56,            60,
     41,        45,        57,            61,
     42,        46,        58,            62,
     43,        47,        59,            63

};
static uint32_t tab4x16[64] = {
    0,    1,    2,    3,    8,    9,    10,    11,    16,    17,    18,    19,    24,    25,    26,    27,
    4,    5,    6,    7,    12,    13,    14,    15,    20,    21,    22,    23,    28,    29,    30,    31,
    32,    33,    34,    35,    40,    41,    42,    43,    48,    49,    50,    51,    56,    57,    58,    59,
    36,    37,    38,    39,    44,    45,    46,    47,    52,    53,    54,    55,    60,    61,    62,    63
};

static uint32_t tab64x16[4] = {
    0,    1,    2,    3,
};

static uint32_t tab16x64[4] = {
    0,    1,    2,    3,
};

/***************************************************************
* in_loop_me_8xN_Nx8_distortion_update
*  Compute the distortion at a given position and update
*  the best for the supported 8xN and Nx8 blocks
***************************************************************/
static void in_loop_me_8xN_Nx8_distortion_update(
    //Inputs
    uint32_t  curr_mv,
    uint32_t    block_4x4_index,
    uint32_t    *dist_4x4,
    //Outputs
    uint32_t    *best_mv_8x4,
    uint32_t    *best_dist_8x4,
    uint32_t    *dist_8x4,
    uint32_t    *best_mv_4x8,
    uint32_t    *best_dist_4x8,
    uint32_t    *dist_4x8,
    uint32_t    *best_mv_8x8,
    uint32_t    *best_dist_8x8,
    uint32_t    *dist_8x8)
{
    uint32_t square_block_index;
    uint32_t first_rec_block_index;
    uint32_t second_rec_block_index;


    //8x4
    first_rec_block_index = (block_4x4_index - 3) / 2;
    second_rec_block_index = first_rec_block_index + 1;

    dist_8x4[first_rec_block_index] = dist_4x4[block_4x4_index - 3] + dist_4x4[block_4x4_index - 2];

    if (dist_8x4[first_rec_block_index] < best_dist_8x4[first_rec_block_index]) {
        best_mv_8x4[first_rec_block_index] = curr_mv;
        best_dist_8x4[first_rec_block_index] = dist_8x4[first_rec_block_index];
    }

    dist_8x4[second_rec_block_index] = dist_4x4[block_4x4_index - 1] + dist_4x4[block_4x4_index];

    if (dist_8x4[second_rec_block_index] < best_dist_8x4[second_rec_block_index]) {
        best_mv_8x4[second_rec_block_index] = curr_mv;
        best_dist_8x4[second_rec_block_index] = dist_8x4[second_rec_block_index];
    }

    //4x8
    dist_4x8[first_rec_block_index] = dist_4x4[block_4x4_index - 3] + dist_4x4[block_4x4_index - 1];

    if (dist_4x8[first_rec_block_index] < best_dist_4x8[first_rec_block_index]) {
        best_mv_4x8[first_rec_block_index] = curr_mv;
        best_dist_4x8[first_rec_block_index] = dist_4x8[first_rec_block_index];
    }

    dist_4x8[second_rec_block_index] = dist_4x4[block_4x4_index - 2] + dist_4x4[block_4x4_index];

    if (dist_4x8[second_rec_block_index] < best_dist_4x8[second_rec_block_index]) {
        best_mv_4x8[second_rec_block_index] = curr_mv;
        best_dist_4x8[second_rec_block_index] = dist_4x8[second_rec_block_index];
    }

    //8x8
    square_block_index = (block_4x4_index - 3) / 4;

    dist_8x8[square_block_index] = dist_4x8[first_rec_block_index] + dist_4x8[second_rec_block_index];

    if (dist_8x8[square_block_index] < best_dist_8x8[square_block_index]) {
        best_mv_8x8[square_block_index] = curr_mv;
        best_dist_8x8[square_block_index] = dist_8x8[square_block_index];
    }
}
/***************************************************************
* in_loop_me_16xN_Nx16_distortion_update
*  Compute the distortion at a given position and update
*  the best for the supported 16xN and Nx16 blocks
***************************************************************/
static void in_loop_me_16xN_Nx16_distortion_update(
    //Inputs
    uint32_t  curr_mv,
    uint32_t  block_8x8_index,
    uint32_t    block_4x4_index,
    uint32_t    *dist_8x4,
    uint32_t    *dist_4x8,
    uint32_t    *dist_8x8,
    //Outputs
    uint32_t    *best_mv_16x4,
    uint32_t    *best_dist_16x4,
    uint32_t    *dist_16x4,
    uint32_t    *best_mv_16x8,
    uint32_t    *best_dist_16x8,
    uint32_t    *dist_16x8,
    uint32_t    *best_mv_4x16,
    uint32_t    *best_dist_4x16,
    uint32_t    *dist_4x16,
    uint32_t    *best_mv_8x16,
    uint32_t    *best_dist_8x16,
    uint32_t    *dist_8x16,
    uint32_t    *best_mv_16x16,
    uint32_t    *best_dist_16x16,
    uint32_t    *dist_16x16
)
{
    uint32_t square_block_index;
    uint32_t first_rec_block_index;
    uint32_t second_rec_block_index;
    uint32_t third_rec_block_index;
    uint32_t fourth_rec_block_index;
    uint32_t start_index;
    //16x4
    first_rec_block_index = (block_8x8_index - 3);
    second_rec_block_index = first_rec_block_index + 1;
    third_rec_block_index = second_rec_block_index + 1;
    fourth_rec_block_index = third_rec_block_index + 1;

    start_index = (block_4x4_index - 15) >> 1;

    dist_16x4[first_rec_block_index] = dist_8x4[start_index] + dist_8x4[start_index + 2];

    if (dist_16x4[first_rec_block_index] < best_dist_16x4[first_rec_block_index]) {
        best_mv_16x4[first_rec_block_index] = curr_mv;
        best_dist_16x4[first_rec_block_index] = dist_16x4[first_rec_block_index];
    }

    dist_16x4[second_rec_block_index] = dist_8x4[start_index + 1] + dist_8x4[start_index + 3];

    if (dist_16x4[second_rec_block_index] < best_dist_16x4[second_rec_block_index]) {
        best_mv_16x4[second_rec_block_index] = curr_mv;
        best_dist_16x4[second_rec_block_index] = dist_16x4[second_rec_block_index];
    }

    dist_16x4[third_rec_block_index] = dist_8x4[start_index + 4] + dist_8x4[start_index + 6];

    if (dist_16x4[third_rec_block_index] < best_dist_16x4[third_rec_block_index]) {
        best_mv_16x4[third_rec_block_index] = curr_mv;
        best_dist_16x4[third_rec_block_index] = dist_16x4[third_rec_block_index];
    }

    dist_16x4[fourth_rec_block_index] = dist_8x4[start_index + 5] + dist_8x4[start_index + 7];

    if (dist_16x4[fourth_rec_block_index] < best_dist_16x4[fourth_rec_block_index]) {
        best_mv_16x4[fourth_rec_block_index] = curr_mv;
        best_dist_16x4[fourth_rec_block_index] = dist_16x4[fourth_rec_block_index];
    }

    //4x16

    dist_4x16[first_rec_block_index] = dist_4x8[start_index] + dist_4x8[start_index + 4];

    if (dist_4x16[first_rec_block_index] < best_dist_4x16[first_rec_block_index]) {
        best_mv_4x16[first_rec_block_index] = curr_mv;
        best_dist_4x16[first_rec_block_index] = dist_4x16[first_rec_block_index];
    }

    dist_4x16[second_rec_block_index] = dist_4x8[start_index + 1] + dist_4x8[start_index + 5];

    if (dist_4x16[second_rec_block_index] < best_dist_4x16[second_rec_block_index]) {
        best_mv_4x16[second_rec_block_index] = curr_mv;
        best_dist_4x16[second_rec_block_index] = dist_4x16[second_rec_block_index];
    }

    dist_4x16[third_rec_block_index] = dist_4x8[start_index + 2] + dist_4x8[start_index + 6];

    if (dist_4x16[third_rec_block_index] < best_dist_4x16[third_rec_block_index]) {
        best_mv_4x16[third_rec_block_index] = curr_mv;
        best_dist_4x16[third_rec_block_index] = dist_4x16[third_rec_block_index];
    }

    dist_4x16[fourth_rec_block_index] = dist_4x8[start_index + 3] + dist_4x8[start_index + 7];

    if (dist_4x16[fourth_rec_block_index] < best_dist_4x16[fourth_rec_block_index]) {
        best_mv_4x16[fourth_rec_block_index] = curr_mv;
        best_dist_4x16[fourth_rec_block_index] = dist_4x16[fourth_rec_block_index];
    }

    //16x8
    first_rec_block_index = (block_8x8_index - 3) / 2;
    second_rec_block_index = first_rec_block_index + 1;

    dist_16x8[first_rec_block_index] = dist_8x8[block_8x8_index - 3] + dist_8x8[block_8x8_index - 2];

    if (dist_16x8[first_rec_block_index] < best_dist_16x8[first_rec_block_index]) {
        best_mv_16x8[first_rec_block_index] = curr_mv;
        best_dist_16x8[first_rec_block_index] = dist_16x8[first_rec_block_index];
    }

    dist_16x8[second_rec_block_index] = dist_8x8[block_8x8_index - 1] + dist_8x8[block_8x8_index];

    if (dist_16x8[second_rec_block_index] < best_dist_16x8[second_rec_block_index]) {
        best_mv_16x8[second_rec_block_index] = curr_mv;
        best_dist_16x8[second_rec_block_index] = dist_16x8[second_rec_block_index];
    }

    //8x16
    dist_8x16[first_rec_block_index] = dist_8x8[block_8x8_index - 3] + dist_8x8[block_8x8_index - 1];

    if (dist_8x16[first_rec_block_index] < best_dist_8x16[first_rec_block_index]) {
        best_mv_8x16[first_rec_block_index] = curr_mv;
        best_dist_8x16[first_rec_block_index] = dist_8x16[first_rec_block_index];
    }

    dist_8x16[second_rec_block_index] = dist_8x8[block_8x8_index - 2] + dist_8x8[block_8x8_index];

    if (dist_8x16[second_rec_block_index] < best_dist_8x16[second_rec_block_index]) {
        best_mv_8x16[second_rec_block_index] = curr_mv;
        best_dist_8x16[second_rec_block_index] = dist_8x16[second_rec_block_index];
    }

    //16x16
    square_block_index = (block_8x8_index - 3) / 4;

    dist_16x16[square_block_index] = dist_16x8[first_rec_block_index] + dist_16x8[second_rec_block_index];

    if (dist_16x16[square_block_index] < best_dist_16x16[square_block_index]) {
        best_mv_16x16[square_block_index] = curr_mv;
        best_dist_16x16[square_block_index] = dist_16x16[square_block_index];
    }
}
/***************************************************************
* in_loop_me_32xN_Nx32_distortion_update
*  Compute the distortion at a given position and update
*  the best for the supported 32xN and Nx32 blocks
***************************************************************/
static void in_loop_me_32xN_Nx32_distortion_update(
    //Inputs
    uint32_t  curr_mv,
    uint32_t  block_16x16_index,
    uint32_t    block_8x8_index,
    uint32_t    *dist_16x8,
    uint32_t    *dist_8x16,
    uint32_t    *dist_16x16,
    //Outputs
    uint32_t    *best_mv_32x8,
    uint32_t    *best_dist_32x8,
    uint32_t    *dist_32x8,
    uint32_t    *best_mv_32x16,
    uint32_t    *best_dist_32x16,
    uint32_t    *dist_32x16,
    uint32_t    *best_mv_8x32,
    uint32_t    *best_dist_8x32,
    uint32_t    *dist_8x32,
    uint32_t    *best_mv_16x32,
    uint32_t    *best_dist_16x32,
    uint32_t    *dist_16x32,
    uint32_t    *best_mv_32x32,
    uint32_t    *best_dist_32x32,
    uint32_t    *dist_32x32
)
{
    uint32_t square_block_index;
    uint32_t first_rec_block_index;
    uint32_t second_rec_block_index;
    uint32_t third_rec_block_index;
    uint32_t fourth_rec_block_index;
    uint32_t start_index;

    //32x8
    first_rec_block_index = (block_16x16_index - 3);
    second_rec_block_index = first_rec_block_index + 1;
    third_rec_block_index = second_rec_block_index + 1;
    fourth_rec_block_index = third_rec_block_index + 1;

    start_index = (block_8x8_index - 15) >> 1;

    dist_32x8[first_rec_block_index] = dist_16x8[start_index] + dist_16x8[start_index + 2];

    if (dist_32x8[first_rec_block_index] < best_dist_32x8[first_rec_block_index]) {
        best_mv_32x8[first_rec_block_index] = curr_mv;
        best_dist_32x8[first_rec_block_index] = dist_32x8[first_rec_block_index];
    }

    dist_32x8[second_rec_block_index] = dist_16x8[start_index + 1] + dist_16x8[start_index + 3];

    if (dist_32x8[second_rec_block_index] < best_dist_32x8[second_rec_block_index]) {
        best_mv_32x8[second_rec_block_index] = curr_mv;
        best_dist_32x8[second_rec_block_index] = dist_32x8[second_rec_block_index];
    }

    dist_32x8[third_rec_block_index] = dist_16x8[start_index + 4] + dist_16x8[start_index + 6];

    if (dist_32x8[third_rec_block_index] < best_dist_32x8[third_rec_block_index]) {
        best_mv_32x8[third_rec_block_index] = curr_mv;
        best_dist_32x8[third_rec_block_index] = dist_32x8[third_rec_block_index];
    }

    dist_32x8[fourth_rec_block_index] = dist_16x8[start_index + 5] + dist_16x8[start_index + 7];

    if (dist_32x8[fourth_rec_block_index] < best_dist_32x8[fourth_rec_block_index]) {
        best_mv_32x8[fourth_rec_block_index] = curr_mv;
        best_dist_32x8[fourth_rec_block_index] = dist_32x8[fourth_rec_block_index];
    }

    //8x32

    dist_8x32[first_rec_block_index] = dist_8x16[start_index] + dist_8x16[start_index + 4];

    if (dist_8x32[first_rec_block_index] < best_dist_8x32[first_rec_block_index]) {
        best_mv_8x32[first_rec_block_index] = curr_mv;
        best_dist_8x32[first_rec_block_index] = dist_8x32[first_rec_block_index];
    }

    dist_8x32[second_rec_block_index] = dist_8x16[start_index + 1] + dist_8x16[start_index + 5];

    if (dist_8x32[second_rec_block_index] < best_dist_8x32[second_rec_block_index]) {
        best_mv_8x32[second_rec_block_index] = curr_mv;
        best_dist_8x32[second_rec_block_index] = dist_8x32[second_rec_block_index];
    }

    dist_8x32[third_rec_block_index] = dist_8x16[start_index + 2] + dist_8x16[start_index + 6];

    if (dist_8x32[third_rec_block_index] < best_dist_8x32[third_rec_block_index]) {
        best_mv_8x32[third_rec_block_index] = curr_mv;
        best_dist_8x32[third_rec_block_index] = dist_8x32[third_rec_block_index];
    }

    dist_8x32[fourth_rec_block_index] = dist_8x16[start_index + 3] + dist_8x16[start_index + 7];

    if (dist_8x32[fourth_rec_block_index] < best_dist_8x32[fourth_rec_block_index]) {
        best_mv_8x32[fourth_rec_block_index] = curr_mv;
        best_dist_8x32[fourth_rec_block_index] = dist_8x32[fourth_rec_block_index];
    }

    //32x16
    first_rec_block_index = (block_16x16_index - 3) / 2;
    second_rec_block_index = first_rec_block_index + 1;

    dist_32x16[first_rec_block_index] = dist_16x16[block_16x16_index - 3] + dist_16x16[block_16x16_index - 2];

    if (dist_32x16[first_rec_block_index] < best_dist_32x16[first_rec_block_index]) {
        best_mv_32x16[first_rec_block_index] = curr_mv;
        best_dist_32x16[first_rec_block_index] = dist_32x16[first_rec_block_index];
    }

    dist_32x16[second_rec_block_index] = dist_16x16[block_16x16_index - 1] + dist_16x16[block_16x16_index];

    if (dist_32x16[second_rec_block_index] < best_dist_32x16[second_rec_block_index]) {
        best_mv_32x16[second_rec_block_index] = curr_mv;
        best_dist_32x16[second_rec_block_index] = dist_32x16[second_rec_block_index];
    }

    //16x32
    dist_16x32[first_rec_block_index] = dist_16x16[block_16x16_index - 3] + dist_16x16[block_16x16_index - 1];

    if (dist_16x32[first_rec_block_index] < best_dist_16x32[first_rec_block_index]) {
        best_mv_16x32[first_rec_block_index] = curr_mv;
        best_dist_16x32[first_rec_block_index] = dist_16x32[first_rec_block_index];
    }

    dist_16x32[second_rec_block_index] = dist_16x16[block_16x16_index - 2] + dist_16x16[block_16x16_index];

    if (dist_16x32[second_rec_block_index] < best_dist_16x32[second_rec_block_index]) {
        best_mv_16x32[second_rec_block_index] = curr_mv;
        best_dist_16x32[second_rec_block_index] = dist_16x32[second_rec_block_index];
    }

    //32x32
    square_block_index = (block_16x16_index - 3) / 4;

    dist_32x32[square_block_index] = dist_32x16[first_rec_block_index] + dist_32x16[second_rec_block_index];

    if (dist_32x32[square_block_index] < best_dist_32x32[square_block_index]) {
        best_mv_32x32[square_block_index] = curr_mv;
        best_dist_32x32[square_block_index] = dist_32x32[square_block_index];
    }
}
/***************************************************************
* in_loop_me_64xN_Nx64_distortion_update
*  Compute the distortion at a given position and update
*  the best for the supported 64xN and Nx64 blocks
***************************************************************/
static void in_loop_me_64xN_Nx64_distortion_update(
    uint32_t     curr_mv,
    uint32_t     block_32x32_index,
    uint32_t     block_16x16_index,
    uint32_t    *dist_32x16,
    uint32_t    *dist_16x32,
    uint32_t    *dist_32x32,
    uint32_t    *best_mv_64x16,
    uint32_t    *best_dist_64x16,
    uint32_t    *dist_64x16,
    uint32_t    *best_mv_64x32,
    uint32_t    *best_dist_64x32,
    uint32_t    *dist_64x32,
    uint32_t    *best_mv_16x64,
    uint32_t    *best_dist_16x64,
    uint32_t    *dist_16x64,
    uint32_t    *best_mv_32x64,
    uint32_t    *best_dist_32x64,
    uint32_t    *dist_32x64,
    uint32_t    *best_mv_64x64,
    uint32_t    *best_dist_64x64,
    uint32_t    *dist_64x64)
{
    uint32_t square_block_index;
    uint32_t first_rec_block_index;
    uint32_t second_rec_block_index;
    uint32_t third_rec_block_index;
    uint32_t fourth_rec_block_index;
    uint32_t start_index;
    UNUSED(dist_64x32);
    UNUSED(dist_32x64);
    //64x16
    first_rec_block_index = (block_32x32_index - 3);
    second_rec_block_index = first_rec_block_index + 1;
    third_rec_block_index = second_rec_block_index + 1;
    fourth_rec_block_index = third_rec_block_index + 1;

    start_index = (block_16x16_index - 15) >> 1;

    dist_64x16[first_rec_block_index] = dist_32x16[start_index] + dist_32x16[start_index + 2];

    if (dist_64x16[first_rec_block_index] < best_dist_64x16[first_rec_block_index]) {
        best_mv_64x16[first_rec_block_index] = curr_mv;
        best_dist_64x16[first_rec_block_index] = dist_64x16[first_rec_block_index];
    }

    dist_64x16[second_rec_block_index] = dist_32x16[start_index + 1] + dist_32x16[start_index + 3];

    if (dist_64x16[second_rec_block_index] < best_dist_64x16[second_rec_block_index]) {
        best_mv_64x16[second_rec_block_index] = curr_mv;
        best_dist_64x16[second_rec_block_index] = dist_64x16[second_rec_block_index];
    }

    dist_64x16[third_rec_block_index] = dist_32x16[start_index + 4] + dist_32x16[start_index + 6];

    if (dist_64x16[third_rec_block_index] < best_dist_64x16[third_rec_block_index]) {
        best_mv_64x16[third_rec_block_index] = curr_mv;
        best_dist_64x16[third_rec_block_index] = dist_64x16[third_rec_block_index];
    }

    dist_64x16[fourth_rec_block_index] = dist_32x16[start_index + 5] + dist_32x16[start_index + 7];

    if (dist_64x16[fourth_rec_block_index] < best_dist_64x16[fourth_rec_block_index]) {
        best_mv_64x16[fourth_rec_block_index] = curr_mv;
        best_dist_64x16[fourth_rec_block_index] = dist_64x16[fourth_rec_block_index];
    }

    //16x64

    dist_16x64[first_rec_block_index] = dist_16x32[start_index] + dist_16x32[start_index + 4];

    if (dist_16x64[first_rec_block_index] < best_dist_16x64[first_rec_block_index]) {
        best_mv_16x64[first_rec_block_index] = curr_mv;
        best_dist_16x64[first_rec_block_index] = dist_16x64[first_rec_block_index];
    }

    dist_16x64[second_rec_block_index] = dist_16x32[start_index + 1] + dist_16x32[start_index + 5];

    if (dist_16x64[second_rec_block_index] < best_dist_16x64[second_rec_block_index]) {
        best_mv_16x64[second_rec_block_index] = curr_mv;
        best_dist_16x64[second_rec_block_index] = dist_16x64[second_rec_block_index];
    }

    dist_16x64[third_rec_block_index] = dist_16x32[start_index + 2] + dist_16x32[start_index + 6];

    if (dist_16x64[third_rec_block_index] < best_dist_16x64[third_rec_block_index]) {
        best_mv_16x64[third_rec_block_index] = curr_mv;
        best_dist_16x64[third_rec_block_index] = dist_16x64[third_rec_block_index];
    }

    dist_16x64[fourth_rec_block_index] = dist_16x32[start_index + 3] + dist_16x32[start_index + 7];

    if (dist_16x64[fourth_rec_block_index] < best_dist_16x64[fourth_rec_block_index]) {
        best_mv_16x64[fourth_rec_block_index] = curr_mv;
        best_dist_16x64[fourth_rec_block_index] = dist_16x64[fourth_rec_block_index];
    }

    //64x32
    first_rec_block_index = (block_32x32_index - 3) / 2;
    second_rec_block_index = first_rec_block_index + 1;

    dist_64x32[first_rec_block_index] = dist_32x32[block_32x32_index - 3] + dist_32x32[block_32x32_index - 2];

    if (dist_64x32[first_rec_block_index] < best_dist_64x32[first_rec_block_index]) {
        best_mv_64x32[first_rec_block_index] = curr_mv;
        best_dist_64x32[first_rec_block_index] = dist_64x32[first_rec_block_index];
    }

    dist_64x32[second_rec_block_index] = dist_32x32[block_32x32_index - 1] + dist_32x32[block_32x32_index];

    if (dist_64x32[second_rec_block_index] < best_dist_64x32[second_rec_block_index]) {
        best_mv_64x32[second_rec_block_index] = curr_mv;
        best_dist_64x32[second_rec_block_index] = dist_64x32[second_rec_block_index];
    }

    //32x64
    dist_32x64[first_rec_block_index] = dist_32x32[block_32x32_index - 3] + dist_32x32[block_32x32_index - 1];

    if (dist_32x64[first_rec_block_index] < best_dist_32x64[first_rec_block_index]) {
        best_mv_32x64[first_rec_block_index] = curr_mv;
        best_dist_32x64[first_rec_block_index] = dist_32x64[first_rec_block_index];
    }

    dist_32x64[second_rec_block_index] = dist_32x32[block_32x32_index - 2] + dist_32x32[block_32x32_index];

    if (dist_32x64[second_rec_block_index] < best_dist_32x64[second_rec_block_index]) {
        best_mv_32x64[second_rec_block_index] = curr_mv;
        best_dist_32x64[second_rec_block_index] = dist_32x64[second_rec_block_index];
    }

    //64x64
    square_block_index = (block_32x32_index - 3) / 4;

    dist_64x64[square_block_index] = dist_64x32[first_rec_block_index] + dist_64x32[second_rec_block_index];

    if (dist_64x64[square_block_index] < best_dist_64x64[square_block_index]) {
        best_mv_64x64[square_block_index] = curr_mv;
        best_dist_64x64[square_block_index] = dist_64x64[square_block_index];
    }
}

/***************************************************************
* in_loop_me_128xN_Nx128_distortion_update
*  Compute the distortion at a given position and update
*  the best for the supported 128xN and Nx128 blocks
***************************************************************/
static void in_loop_me_128xN_Nx128_distortion_update(
    uint32_t     curr_mv,
    uint32_t     block_64x64_index,
    uint32_t     block_32x32_index,
    uint32_t    *dist_64x32,
    uint32_t    *dist_32x64,
    uint32_t    *dist_64x64,
    uint32_t    *best_mv_128x64,
    uint32_t    *best_dist_128x64,
    uint32_t    *dist_128x64,
    uint32_t    *best_mv_64x128,
    uint32_t    *best_dist_64x128,
    uint32_t    *dist_64x128,
    uint32_t    *best_mv_128x128,
    uint32_t    *best_dist_128x128,
    uint32_t    *dist_128x128
)
{
    uint32_t square_block_index;
    uint32_t first_rec_block_index;
    uint32_t second_rec_block_index;
    UNUSED(block_32x32_index);
    UNUSED(dist_64x32);
    UNUSED(dist_32x64);
    //128x64
    first_rec_block_index = (block_64x64_index - 3) / 4;
    second_rec_block_index = first_rec_block_index + 1;

    dist_128x64[first_rec_block_index] = dist_64x64[block_64x64_index - 3] + dist_64x64[block_64x64_index - 2];

    if (dist_128x64[first_rec_block_index] < best_dist_128x64[first_rec_block_index]) {
        best_mv_128x64[first_rec_block_index] = curr_mv;
        best_dist_128x64[first_rec_block_index] = dist_128x64[first_rec_block_index];
    }

    dist_128x64[second_rec_block_index] = dist_64x64[block_64x64_index - 1] + dist_64x64[block_64x64_index];

    if (dist_128x64[second_rec_block_index] < best_dist_128x64[second_rec_block_index]) {
        best_mv_128x64[second_rec_block_index] = curr_mv;
        best_dist_128x64[second_rec_block_index] = dist_128x64[second_rec_block_index];
    }

    //64x128
    dist_64x128[first_rec_block_index] = dist_64x64[block_64x64_index - 3] + dist_64x64[block_64x64_index - 1];

    if (dist_64x128[first_rec_block_index] < best_dist_64x128[first_rec_block_index]) {
        best_mv_64x128[first_rec_block_index] = curr_mv;
        best_dist_64x128[first_rec_block_index] = dist_64x128[first_rec_block_index];
    }

    dist_64x128[second_rec_block_index] = dist_64x64[block_64x64_index - 2] + dist_64x64[block_64x64_index];

    if (dist_64x128[second_rec_block_index] < best_dist_64x128[second_rec_block_index]) {
        best_mv_64x128[second_rec_block_index] = curr_mv;
        best_dist_64x128[second_rec_block_index] = dist_64x128[second_rec_block_index];
    }

    //128x128
    square_block_index = (block_64x64_index - 3) / 4;

    *dist_128x128 = dist_128x64[first_rec_block_index] + dist_128x64[second_rec_block_index];

    if (*dist_128x128 < best_dist_128x128[square_block_index]) {
        best_mv_128x128[square_block_index] = curr_mv;
        best_dist_128x128[square_block_index] = *dist_128x128;
    }
}
/***************************************************************
* in_loop_me_get_search_point_results_block
*  Compute the distortion at a given position
***************************************************************/

static void in_loop_me_get_search_point_results_block(
    SsMeContext_t            *context_ptr,                    // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t                   list_index,                      // input parameter, reference list index
    uint32_t                   ref_index,
    int32_t                   x_search_index,                  // input parameter, search region position in the horizontal direction, used to derive xMV
    int32_t                   y_search_index,                  // input parameter, search region position in the vertical direction, used to derive yMV
    uint32_t                   number_of_sb_quad,
    EbAsm                   asm_type)
{
    uint8_t  *src_ptr = context_ptr->sb_buffer;

    // NADER
    uint8_t   *ref_ptr = context_ptr->integer_buffer_ptr[list_index][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[list_index][0]);
    uint32_t   ref_luma_stride = context_ptr->interpolated_full_stride[list_index][0];
    uint32_t   curr_mv_1 = (((uint16_t)y_search_index) << 18);
    uint16_t   curr_mv_2 = (((uint16_t)x_search_index << 2));
    uint32_t   curr_mv = curr_mv_1 | curr_mv_2;
    uint32_t  *best_dist_4x4 = context_ptr->p_best_sad4x4;
    uint32_t  *best_mv_4x4 = context_ptr->p_best_mv4x4;
    uint32_t  *dist_4x4 = context_ptr->p_sad4x4;
    uint32_t  *best_dist_8x4 = context_ptr->p_best_sad8x4;
    uint32_t  *best_mv_8x4 = context_ptr->p_best_mv8x4;
    uint32_t  *dist_8x4 = context_ptr->p_sad8x4;
    uint32_t  *best_dist_4x8 = context_ptr->p_best_sad4x8;
    uint32_t  *best_mv_4x8 = context_ptr->p_best_mv4x8;
    uint32_t  *dist_4x8 = context_ptr->p_sad4x8;
    uint32_t  *best_dist_8x8 = context_ptr->p_best_sad8x8;
    uint32_t  *best_mv_8x8 = context_ptr->p_best_mv8x8;
    uint32_t  *dist_8x8 = context_ptr->p_sad8x8;
    uint32_t  *best_dist_16x16 = context_ptr->p_best_sad16x16;
    uint32_t  *best_mv_16x16 = context_ptr->p_best_mv16x16;
    uint32_t  *dist_16x16 = context_ptr->p_sad16x16;
    uint32_t  *best_dist_16x8 = context_ptr->p_best_sad16x8;
    uint32_t  *best_mv_16x8 = context_ptr->p_best_mv16x8;
    uint32_t  *dist_16x8 = context_ptr->p_sad16x8;
    uint32_t  *best_dist_16x4 = context_ptr->p_best_sad16x4;
    uint32_t  *best_mv_16x4 = context_ptr->p_best_mv16x4;
    uint32_t  *dist_16x4 = context_ptr->p_sad16x4;
    uint32_t  *best_dist_8x16 = context_ptr->p_best_sad8x16;
    uint32_t  *best_mv_8x16 = context_ptr->p_best_mv8x16;
    uint32_t  *dist_8x16 = context_ptr->p_sad8x16;
    uint32_t  *best_dist_4x16 = context_ptr->p_best_sad4x16;
    uint32_t  *best_mv_4x16 = context_ptr->p_best_mv4x16;
    uint32_t  *dist_4x16 = context_ptr->p_sad4x16;
    uint32_t  *best_dist_32x8 = context_ptr->p_best_sad32x8;
    uint32_t  *best_mv_32x8 = context_ptr->p_best_mv32x8;
    uint32_t  *dist_32x8 = context_ptr->p_sad32x8;
    uint32_t  *best_dist_32x16 = context_ptr->p_best_sad32x16;
    uint32_t  *best_mv_32x16 = context_ptr->p_best_mv32x16;
    uint32_t  *dist_32x16 = context_ptr->p_sad32x16;
    uint32_t  *best_dist_16x32 = context_ptr->p_best_sad16x32;
    uint32_t  *best_mv_16x32 = context_ptr->p_best_mv16x32;
    uint32_t  *dist_16x32 = context_ptr->p_sad16x32;
    uint32_t  *best_dist_8x32 = context_ptr->p_best_sad8x32;
    uint32_t  *best_mv_8x32 = context_ptr->p_best_mv8x32;
    uint32_t  *dist_8x32 = context_ptr->p_sad8x32;
    uint32_t  *best_dist_32x32 = context_ptr->p_best_sad32x32;
    uint32_t  *best_mv_32x32 = context_ptr->p_best_mv32x32;
    uint32_t  *dist_32x32 = context_ptr->p_sad32x32;
    uint32_t  *best_dist_64x16 = context_ptr->p_best_sad64x16;
    uint32_t  *best_mv_64x16 = context_ptr->p_best_mv64x16;
    uint32_t  *dist_64x16 = context_ptr->p_sad64x16;
    uint32_t  *best_dist_64x32 = context_ptr->p_best_sad64x32;
    uint32_t  *best_mv_64x32 = context_ptr->p_best_mv64x32;
    uint32_t  *dist_64x32 = context_ptr->p_sad64x32;
    uint32_t  *best_dist_32x64 = context_ptr->p_best_sad32x64;
    uint32_t  *best_mv_32x64 = context_ptr->p_best_mv32x64;
    uint32_t  *dist_32x64 = context_ptr->p_sad32x64;
    uint32_t  *best_dist_16x64 = context_ptr->p_best_sad16x64;
    uint32_t  *best_mv_16x64 = context_ptr->p_best_mv16x64;
    uint32_t  *dist_16x64 = context_ptr->p_sad16x64;
    uint32_t  *best_dist_64x64 = context_ptr->p_best_sad64x64;
    uint32_t  *best_mv_64x64 = context_ptr->p_best_mv64x64;
    uint32_t  *dist_64x64 = context_ptr->p_sad64x64;
    uint32_t  *best_dist_128x64 = context_ptr->p_best_sad128x64;
    uint32_t  *best_mv_128x64 = context_ptr->p_best_mv128x64;
    uint32_t  *dist_128x64 = context_ptr->p_sad128x64;
    uint32_t  *best_dist_64x128 = context_ptr->p_best_sad64x128;
    uint32_t  *best_mv_64x128 = context_ptr->p_best_mv64x128;
    uint32_t  *dist_64x128 = context_ptr->p_sad64x128;
    uint32_t  *best_dist_128x128 = context_ptr->p_best_sad128x128;
    uint32_t  *best_mv_128x128 = context_ptr->p_best_mv128x128;
    uint32_t  dist_128x128 = context_ptr->p_sad128x128;
    const uint32_t  src_stride = context_ptr->sb_buffer_stride;
    uint32_t block_64x64_index;
    uint32_t block_32x32_index;
    uint32_t block_16x16_index;
    uint32_t block_8x8_index;
    uint32_t block_4x4_index;
    uint32_t block_64x64_x;
    uint32_t block_32x32_x;
    uint32_t block_16x16_x;
    uint32_t block_8x8_x;
    uint32_t block_4x4_x;
    uint32_t block_64x64_y;
    uint32_t block_32x32_y;
    uint32_t block_16x16_y;
    uint32_t block_8x8_y;
    uint32_t block_4x4_y;
    uint32_t quad_offset = number_of_sb_quad > 1 ? 2 : 1;

    for (block_64x64_y = 0; block_64x64_y < quad_offset; block_64x64_y++) {
        for (block_64x64_x = 0; block_64x64_x < quad_offset; block_64x64_x++) {

            block_64x64_index = block_64x64_x + (block_64x64_y * 2);

            for (block_32x32_y = 0; block_32x32_y < 2; block_32x32_y++) {
                for (block_32x32_x = 0; block_32x32_x < 2; block_32x32_x++) {

                    block_32x32_index = (block_64x64_index * 4) + block_32x32_x + (block_32x32_y * 2);

                    for (block_16x16_y = 0; block_16x16_y < 2; block_16x16_y++) {
                        for (block_16x16_x = 0; block_16x16_x < 2; block_16x16_x++) {

                            block_16x16_index = (block_32x32_index * 4) + block_16x16_x + (block_16x16_y * 2);

                            for (block_8x8_y = 0; block_8x8_y < 2; block_8x8_y++) {
                                for (block_8x8_x = 0; block_8x8_x < 2; block_8x8_x++) {

                                    block_8x8_index = (block_16x16_index * 4) + block_8x8_x + (block_8x8_y * 2);

                                    for (block_4x4_y = 0; block_4x4_y < 2; block_4x4_y++) {
                                        for (block_4x4_x = 0; block_4x4_x < 2; block_4x4_x++) {

                                            block_4x4_index = (block_8x8_index * 4) + block_4x4_x + (block_4x4_y * 2);

                                            uint32_t block_4x4_addr_y = (block_64x64_y * 64) + (block_32x32_y * 32) + (block_16x16_y * 16) + (block_8x8_y * 8) + (block_4x4_y * 4);
                                            uint32_t block_4x4_addr_x = (block_64x64_x * 64) + (block_32x32_x * 32) + (block_16x16_x * 16) + (block_8x8_x * 8) + (block_4x4_x * 4);
                                            uint32_t block_4x4_addr_src = (block_4x4_addr_y * src_stride) + block_4x4_addr_x;
                                            uint32_t block_4x4_addr_ref = ref_index + ((block_4x4_addr_y * ref_luma_stride) + block_4x4_addr_x);

                                            //4x4
                                            dist_4x4[block_4x4_index] = compute4x4SAD_funcPtrArray[asm_type](
                                                src_ptr + block_4x4_addr_src,
                                                src_stride,
                                                ref_ptr + block_4x4_addr_ref,
                                                ref_luma_stride,
                                                4,
                                                4);

                                            if (dist_4x4[block_4x4_index] < best_dist_4x4[block_4x4_index]) {
                                                best_mv_4x4[block_4x4_index] = curr_mv;
                                                best_dist_4x4[block_4x4_index] = dist_4x4[block_4x4_index];
                                            }
                                        }
                                    }

                                    // Nader - Full-pel search for depth 4 blocks
                                    in_loop_me_8xN_Nx8_distortion_update(
                                        //Inputs
                                        curr_mv,
                                        block_4x4_index,
                                        dist_4x4,
                                        //Outputs
                                        best_mv_8x4,
                                        best_dist_8x4,
                                        dist_8x4,
                                        best_mv_4x8,
                                        best_dist_4x8,
                                        dist_4x8,
                                        best_mv_8x8,
                                        best_dist_8x8,
                                        dist_8x8);
                                }
                            }

                            // Nader - Full-pel search for depth 3 blocks
                            in_loop_me_16xN_Nx16_distortion_update(
                                //Inputs
                                curr_mv,
                                block_8x8_index,
                                block_4x4_index,
                                dist_8x4,
                                dist_4x8,
                                dist_8x8,
                                //Outputs
                                best_mv_16x4,
                                best_dist_16x4,
                                dist_16x4,
                                best_mv_16x8,
                                best_dist_16x8,
                                dist_16x8,
                                best_mv_4x16,
                                best_dist_4x16,
                                dist_4x16,
                                best_mv_8x16,
                                best_dist_8x16,
                                dist_8x16,
                                best_mv_16x16,
                                best_dist_16x16,
                                dist_16x16);
                        }
                    }

                    // Nader - Full-pel search for depth 2 blocks
                    in_loop_me_32xN_Nx32_distortion_update(
                        //Inputs
                        curr_mv,
                        block_16x16_index,
                        block_8x8_index,
                        dist_16x8,
                        dist_8x16,
                        dist_16x16,
                        //Outputs
                        best_mv_32x8,
                        best_dist_32x8,
                        dist_32x8,
                        best_mv_32x16,
                        best_dist_32x16,
                        dist_32x16,
                        best_mv_8x32,
                        best_dist_8x32,
                        dist_8x32,
                        best_mv_16x32,
                        best_dist_16x32,
                        dist_16x32,
                        best_mv_32x32,
                        best_dist_32x32,
                        dist_32x32);
                }
            }

            // Nader - Full-pel search for depth 1 blocks
            in_loop_me_64xN_Nx64_distortion_update(
                //Inputs
                curr_mv,
                block_32x32_index,
                block_16x16_index,
                dist_32x16,
                dist_16x32,
                dist_32x32,
                //Outputs
                best_mv_64x16,
                best_dist_64x16,
                dist_64x16,
                best_mv_64x32,
                best_dist_64x32,
                dist_64x32,
                best_mv_16x64,
                best_dist_16x64,
                dist_16x64,
                best_mv_32x64,
                best_dist_32x64,
                dist_32x64,
                best_mv_64x64,
                best_dist_64x64,
                dist_64x64);
        }
    }


    if (number_of_sb_quad > 1) {

        // Nader - Full-pel search for depth 0 blocks
        in_loop_me_128xN_Nx128_distortion_update(
            //Inputs
            curr_mv,
            block_64x64_index,
            block_32x32_index,
            dist_64x32,
            dist_32x64,
            dist_64x64,
            //Outputs
            best_mv_128x64,
            best_dist_128x64,
            dist_128x64,
            best_mv_64x128,
            best_dist_64x128,
            dist_64x128,
            best_mv_128x128,
            best_dist_128x128,
            &dist_128x128);
    }
}

/***************************************************************
* in_loop_me_fullpel_search_sblock
*  perform the full pel search for the whole super-block
***************************************************************/
static void in_loop_me_fullpel_search_sblock(
    SsMeContext_t            *context_ptr,
    uint32_t                   list_index,
    int16_t                   x_search_area_origin,
    int16_t                     y_search_area_origin,
    uint32_t                   search_area_width,
    uint32_t                   search_area_height,
    uint32_t                   number_of_sb_quad,
    EbAsm                   asm_type)
{
    uint32_t x_search_index, y_search_index;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {

        for (x_search_index = 0; x_search_index < search_area_width; x_search_index++) {

            in_loop_me_get_search_point_results_block(
                context_ptr,
                list_index,
                x_search_index + y_search_index * context_ptr->interpolated_full_stride[list_index][0],
                (int32_t)x_search_index + x_search_area_origin,
                (int32_t)y_search_index + y_search_area_origin,
                number_of_sb_quad,
                asm_type);
        }

    }
}

/***************************************************************
* in_loop_me_context_ctor
*  in-loop motion estimation construtor
***************************************************************/
EbErrorType in_loop_me_context_ctor(
    SsMeContext_t                          **object_dbl_ptr)
{

    uint32_t                   listIndex;
    uint32_t                   refPicIndex;

    EB_MALLOC(SsMeContext_t*, *object_dbl_ptr, sizeof(SsMeContext_t), EB_N_PTR);

    // Intermediate LCU-sized buffer to retain the input samples
    (*object_dbl_ptr)->sb_buffer_stride = MAX_SB_SIZE;


    EB_ALLIGN_MALLOC(uint8_t *, (*object_dbl_ptr)->sb_buffer, sizeof(uint8_t) * MAX_SB_SIZE * (*object_dbl_ptr)->sb_buffer_stride, EB_A_PTR);

    EB_MEMSET((*object_dbl_ptr)->sb_buffer, 0, sizeof(uint8_t) * MAX_SB_SIZE * (*object_dbl_ptr)->sb_buffer_stride);



    (*object_dbl_ptr)->interpolated_stride = MAX_SEARCH_AREA_WIDTH;

    // EB_MALLOC(EB_BitFraction *, (*object_dbl_ptr)->mvd_bits_array, sizeof(EB_BitFraction) * NUMBER_OF_MVD_CASES, EB_N_PTR);
    // 15 intermediate buffers to retain the interpolated reference samples

    //      0    1    2    3
    // 0    A    a    b    c
    // 1    d    e    f    g
    // 2    h    i    j    k
    // 3    n    p    q    r

    //                  _____________
    //                 |             |
    // --I samples --> |Interpolation|-- O samples -->
    //                 | ____________|

    // Before Interpolation: 2 x 3
    //   I   I
    //   I   I
    //   I   I

    // After 1-D Horizontal Interpolation: (2 + 1) x 3 - a, b, and c
    // O I O I O
    // O I O I O
    // O I O I O

    // After 1-D Vertical Interpolation: 2 x (3 + 1) - d, h, and n
    //   O   O
    //   I   I
    //   O   O
    //   I   I
    //   O   O
    //   I   I
    //   O   O

    // After 2-D (Horizontal/Vertical) Interpolation: (2 + 1) x (3 + 1) - e, f, g, i, j, k, n, p, q, and r
    // O   O   O
    //   I   I
    // O   O   O
    //   I   I
    // O   O   O
    //   I   I
    // O   O   O

    for (listIndex = 0; listIndex < MAX_NUM_OF_REF_PIC_LIST; listIndex++) {

        for (refPicIndex = 0; refPicIndex < MAX_REF_IDX; refPicIndex++) {

            EB_MALLOC(uint8_t *, (*object_dbl_ptr)->integer_buffer[listIndex][refPicIndex], sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

            EB_MALLOC(uint8_t *, (*object_dbl_ptr)->pos_b_buffer[listIndex][refPicIndex], sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

            EB_MALLOC(uint8_t *, (*object_dbl_ptr)->pos_h_buffer[listIndex][refPicIndex], sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

            EB_MALLOC(uint8_t *, (*object_dbl_ptr)->pos_j_buffer[listIndex][refPicIndex], sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

        }

    }

    EB_MALLOC(uint8_t *, (*object_dbl_ptr)->avctemp_buffer, sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

    return EB_ErrorNone;
}



/***************************************************************
* in_loop_me_interpolate_search_region_avc_style
*  performs AVC-style interpolation for the whole Search Region
***************************************************************/
static void in_loop_me_interpolate_search_region_avc_style(
    SsMeContext_t           *context_ptr,                       // input/output parameter, ME context ptr, used to get/set interpolated search area Ptr
    uint32_t                   listIndex,                        // Refrence picture list index
    uint8_t                   *searchRegionBuffer,               // input parameter, search region index, used to point to reference samples
    uint32_t                   lumaStride,                       // input parameter, reference Picture stride
    uint32_t                   search_area_width,                  // input parameter, search area width
    uint32_t                   search_area_height,                 // input parameter, search area height
    uint32_t                   inputBitDepth,                    // input parameter, input sample bit depth
    EbAsm                   asm_type)
{

    //      0    1    2    3
    // 0    A    a    b    c
    // 1    d    e    f    g
    // 2    h    i    j    k
    // 3    n    p    q    r

    // Position  Frac-pos Y  Frac-pos X  Horizontal filter  Vertical filter
    // A         0           0           -                  -
    // a         0           1           F0                 -
    // b         0           2           F1                 -
    // c         0           3           F2                 -
    // d         1           0           -                  F0
    // e         1           1           F0                 F0
    // f         1           2           F1                 F0
    // g         1           3           F2                 F0
    // h         2           0           -                  F1
    // i         2           1           F0                 F1
    // j         2           2           F1                 F1
    // k         2           3           F2                 F1
    // n         3           0           -                  F2
    // p         3           1           F0                 F2
    // q         3           2           F1                 F2
    // r         3           3           F2                 F2

    // Start a b c

    // The Search area needs to be a multiple of 8 to align with the ASM kernel
    // Also the search area must be oversized by 2 to account for edge conditions
    uint32_t searchAreaWidthForAsm = ROUND_UP_MUL_8(search_area_width + 2);


    (void)inputBitDepth;
    // Half pel interpolation of the search region using f1 -> pos_b_buffer
    if (searchAreaWidthForAsm) {

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][2](
            searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride - (ME_FILTER_TAP >> 1) + 1,
            lumaStride,
            context_ptr->pos_b_buffer[listIndex][0],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + ME_FILTER_TAP,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2);
    }

    // Half pel interpolation of the search region using f1 -> pos_h_buffer
    if (searchAreaWidthForAsm) {
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][8](
            searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride - 1 + lumaStride,
            lumaStride,
            context_ptr->pos_h_buffer[listIndex][0],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2);
    }

    if (searchAreaWidthForAsm) {
        // Half pel interpolation of the search region using f1 -> pos_j_buffer
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][8](
            context_ptr->pos_b_buffer[listIndex][0] + context_ptr->interpolated_stride,
            context_ptr->interpolated_stride,
            context_ptr->pos_j_buffer[listIndex][0],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2);
    }

    return;
}


/***************************************************************
* in_loop_me_halfpel_refinement_block
*   performs Half Pel refinement for one block
***************************************************************/
static void in_loop_me_halfpel_refinement_block(
    SequenceControlSet_t    *sequence_control_set_ptr,             // input parameter, Sequence control set Ptr
    SsMeContext_t           *context_ptr,                        // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t                   block_index_in_sb_buffer,                  // input parameter, PU origin, used to point to source samples
    uint8_t                   *pos_b_buffer,                        // input parameter, position "b" interpolated search area Ptr
    uint8_t                   *pos_h_buffer,                        // input parameter, position "h" interpolated search area Ptr
    uint8_t                   *pos_j_buffer,                        // input parameter, position "j" interpolated search area Ptr
    uint32_t                   pu_width,                           // input parameter, PU width
    uint32_t                   pu_height,                          // input parameter, PU height
    int16_t                   x_search_area_origin,                 // input parameter, search area origin in the horizontal direction, used to point to reference samples
    int16_t                   y_search_area_origin,                 // input parameter, search area origin in the vertical direction, used to point to reference samples
    EbAsm                   asm_type,
    uint32_t                  *pBestSad,
    uint32_t                  *pBestMV,
    uint8_t                   *psubPelDirection
)
{

    EncodeContext_t         *encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;

    int32_t searchRegionIndex;
    uint64_t bestHalfSad = 0;
    uint64_t distortionLeftPosition = 0;
    uint64_t distortionRightPosition = 0;
    uint64_t distortionTopPosition = 0;
    uint64_t distortionBottomPosition = 0;
    uint64_t distortionTopLeftPosition = 0;
    uint64_t distortionTopRightPosition = 0;
    uint64_t distortionBottomLeftPosition = 0;
    uint64_t distortionBottomRightPosition = 0;

    int16_t xMvHalf[8];
    int16_t yMvHalf[8];

    int16_t x_mv = _MVXT(*pBestMV);
    int16_t y_mv = _MVYT(*pBestMV);
    int16_t xSearchIndex = (x_mv >> 2) - x_search_area_origin;
    int16_t ySearchIndex = (y_mv >> 2) - y_search_area_origin;

    (void)sequence_control_set_ptr;
    (void)encode_context_ptr;

    //TODO : remove these, and update the MV by just shifts

    xMvHalf[0] = x_mv - 2; // L  position
    xMvHalf[1] = x_mv + 2; // R  position
    xMvHalf[2] = x_mv;     // T  position
    xMvHalf[3] = x_mv;     // B  position
    xMvHalf[4] = x_mv - 2; // TL position
    xMvHalf[5] = x_mv + 2; // TR position
    xMvHalf[6] = x_mv + 2; // BR position
    xMvHalf[7] = x_mv - 2; // BL position

    yMvHalf[0] = y_mv;     // L  position
    yMvHalf[1] = y_mv;     // R  position
    yMvHalf[2] = y_mv - 2; // T  position
    yMvHalf[3] = y_mv + 2; // B  position
    yMvHalf[4] = y_mv - 2; // TL position
    yMvHalf[5] = y_mv - 2; // TR position
    yMvHalf[6] = y_mv + 2; // BR position
    yMvHalf[7] = y_mv + 2; // BL position



    // L position
    searchRegionIndex = xSearchIndex + (int16_t)context_ptr->interpolated_stride * ySearchIndex;
#if USE_INLOOP_ME_FULL_SAD
    distortionLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
    distortionLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
    if (distortionLeftPosition < *pBestSad) {
        *pBestSad = (uint32_t)distortionLeftPosition;
        *pBestMV = ((uint16_t)yMvHalf[0] << 16) | ((uint16_t)xMvHalf[0]);
    }

    // R position
    searchRegionIndex++;
#if USE_INLOOP_ME_FULL_SAD
    distortionRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
    distortionRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif

    if (distortionRightPosition < *pBestSad) {
        *pBestSad = (uint32_t)distortionRightPosition;
        *pBestMV = ((uint16_t)yMvHalf[1] << 16) | ((uint16_t)xMvHalf[1]);
    }

    // T position
    searchRegionIndex = xSearchIndex + (int16_t)context_ptr->interpolated_stride * ySearchIndex;
#if USE_INLOOP_ME_FULL_SAD
    distortionTopPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
    distortionTopPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
    if (distortionTopPosition < *pBestSad) {
        *pBestSad = (uint32_t)distortionTopPosition;
        *pBestMV = ((uint16_t)yMvHalf[2] << 16) | ((uint16_t)xMvHalf[2]);
    }


    // B position
    searchRegionIndex += (int16_t)context_ptr->interpolated_stride;
#if USE_INLOOP_ME_FULL_SAD
    distortionBottomPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
    distortionBottomPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
    if (distortionBottomPosition < *pBestSad) {
        *pBestSad = (uint32_t)distortionBottomPosition;
        *pBestMV = ((uint16_t)yMvHalf[3] << 16) | ((uint16_t)xMvHalf[3]);
    }


    //TL position
    searchRegionIndex = xSearchIndex + (int16_t)context_ptr->interpolated_stride * ySearchIndex;
#if USE_INLOOP_ME_FULL_SAD
    distortionTopLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
    distortionTopLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
    if (distortionTopLeftPosition < *pBestSad) {
        *pBestSad = (uint32_t)distortionTopLeftPosition;
        *pBestMV = ((uint16_t)yMvHalf[4] << 16) | ((uint16_t)xMvHalf[4]);
    }


    //TR position
    searchRegionIndex++;
#if USE_INLOOP_ME_FULL_SAD
    distortionTopRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
    distortionTopRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
    if (distortionTopRightPosition < *pBestSad) {
        *pBestSad = (uint32_t)distortionTopRightPosition;
        *pBestMV = ((uint16_t)yMvHalf[5] << 16) | ((uint16_t)xMvHalf[5]);
    }


    //BR position
    searchRegionIndex += (int16_t)context_ptr->interpolated_stride;
#if USE_INLOOP_ME_FULL_SAD
    distortionBottomRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
    distortionBottomRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
    if (distortionBottomRightPosition < *pBestSad) {
        *pBestSad = (uint32_t)distortionBottomRightPosition;
        *pBestMV = ((uint16_t)yMvHalf[6] << 16) | ((uint16_t)xMvHalf[6]);
    }


    //BL position
    searchRegionIndex--;
#if USE_INLOOP_ME_FULL_SAD
    distortionBottomLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
    distortionBottomLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
    if (distortionBottomLeftPosition < *pBestSad) {
        *pBestSad = (uint32_t)distortionBottomLeftPosition;
        *pBestMV = ((uint16_t)yMvHalf[7] << 16) | ((uint16_t)xMvHalf[7]);
    }



    bestHalfSad = MIN(distortionLeftPosition, MIN(distortionRightPosition, MIN(distortionTopPosition, MIN(distortionBottomPosition, MIN(distortionTopLeftPosition, MIN(distortionTopRightPosition, MIN(distortionBottomLeftPosition, distortionBottomRightPosition)))))));


    if (bestHalfSad == distortionLeftPosition) {
        *psubPelDirection = LEFT_POSITION;
    }
    else if (bestHalfSad == distortionRightPosition) {
        *psubPelDirection = RIGHT_POSITION;
    }
    else if (bestHalfSad == distortionTopPosition) {
        *psubPelDirection = TOP_POSITION;
    }
    else if (bestHalfSad == distortionBottomPosition) {
        *psubPelDirection = BOTTOM_POSITION;
    }
    else if (bestHalfSad == distortionTopLeftPosition) {
        *psubPelDirection = TOP_LEFT_POSITION;
    }
    else if (bestHalfSad == distortionTopRightPosition) {
        *psubPelDirection = TOP_RIGHT_POSITION;
    }
    else if (bestHalfSad == distortionBottomLeftPosition) {
        *psubPelDirection = BOTTOM_LEFT_POSITION;
    }
    else if (bestHalfSad == distortionBottomRightPosition) {
        *psubPelDirection = BOTTOM_RIGHT_POSITION;
    }
    return;
}


/***************************************************************
* in_loop_me_halfpel_search_sblock
*   performs Half Pel refinement
***************************************************************/
void in_loop_me_halfpel_search_sblock(
    SequenceControlSet_t    *sequence_control_set_ptr,             // input parameter, Sequence control set Ptr
    SsMeContext_t           *context_ptr,                        // input/output parameter, ME context Ptr, used to get/update ME results
    uint8_t                   *pos_b_buffer,                        // input parameter, position "b" interpolated search area Ptr
    uint8_t                   *pos_h_buffer,                        // input parameter, position "h" interpolated search area Ptr
    uint8_t                   *pos_j_buffer,                        // input parameter, position "j" interpolated search area Ptr
    int16_t                   x_search_area_origin,                 // input parameter, search area origin in the horizontal direction, used to point to reference samples
    int16_t                   y_search_area_origin,                 // input parameter, search area origin in the vertical direction, used to point to reference samples
    EbAsm                   asm_type)
{

    uint32_t idx;
    uint32_t block_index;
    uint32_t block_shift_x;
    uint32_t block_shift_y;
    uint32_t block_index_in_sb_buffer;
    uint32_t posb_buffer_index;
    uint32_t posh_buffer_index;
    uint32_t posj_buffer_index;

    uint32_t block_offset = 0;
    uint32_t x_offset = 0;
    uint32_t y_offset = 0;
    uint32_t quad_index = 0;
    uint32_t number_of_sb_quad = context_ptr->sb_size == BLOCK_128X128 ? 4 : 1;

    // 4x4   [256 4x4 blocks]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 256; ++block_index) {

            block_offset = (quad_index * 256);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab4x4[block_index] + block_offset;
            block_shift_x = ((block_index & 0xf) << 2) + x_offset;
            block_shift_y = ((block_index >> 4) << 2) + y_offset;


            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                4,
                4,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad4x4[idx],
                &context_ptr->p_best_mv4x4[idx],
                &context_ptr->psub_pel_direction4x4[idx]);
        }
    }


    // 8x4   [128 8x4 blocks]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 128; ++block_index) {

            block_offset = (quad_index * 128);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab8x4[block_index] + block_offset;
            block_shift_x = ((block_index & 0x07) << 3) + x_offset;;
            block_shift_y = ((block_index >> 3) << 2) + y_offset;;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                8,
                4,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x4[idx],
                &context_ptr->p_best_mv8x4[idx],
                &context_ptr->psub_pel_direction8x4[idx]);

        }
    }

    // 4x8   [128 4x8 blocks]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 128; ++block_index) {

            block_offset = (quad_index * 128);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab4x8[block_index] + block_offset;
            block_shift_x = ((block_index & 0xf) << 2) + x_offset;
            block_shift_y = ((block_index >> 4) << 3) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                4,
                8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad4x8[idx],
                &context_ptr->p_best_mv4x8[idx],
                &context_ptr->psub_pel_direction4x8[idx]);
        }
    }


    // 8x8   [64 8x8 blocks]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 64; ++block_index) {

            block_offset = (quad_index * 64);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab8x8[block_index] + block_offset;
            block_shift_x = ((block_index & 0x07) << 3) + x_offset;
            block_shift_y = ((block_index >> 3) << 3) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                8,
                8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x8[idx],
                &context_ptr->p_best_mv8x8[idx],
                &context_ptr->psub_pel_direction8x8[idx]);

        }
    }

    // 16x8 [32 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 32; ++block_index) {

            block_offset = (quad_index * 32);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab16x8[block_index] + block_offset;
            block_shift_x = ((block_index & 0x03) << 4) + x_offset;
            block_shift_y = ((block_index >> 2) << 3) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                16,
                8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x8[idx],
                &context_ptr->p_best_mv16x8[idx],
                &context_ptr->psub_pel_direction16x8[idx]);
        }
    }

    // 8x16 [32 partitions]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 32; ++block_index) {

            block_offset = (quad_index * 32);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab8x16[block_index] + block_offset;
            block_shift_x = ((block_index & 0x07) << 3) + x_offset;
            block_shift_y = ((block_index >> 3) << 4) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                8,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x16[idx],
                &context_ptr->p_best_mv8x16[idx],
                &context_ptr->psub_pel_direction8x16[idx]);
        }
    }

    // 32x8 [16 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 16; ++block_index) {

            block_offset = (quad_index * 16);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab32x8[block_index] + block_offset;
            block_shift_x = ((block_index & 0x01) << 5) + x_offset;
            block_shift_y = ((block_index >> 1) << 3) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                32,
                8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x8[idx],
                &context_ptr->p_best_mv32x8[idx],
                &context_ptr->psub_pel_direction32x8[idx]);

        }

    }

    // 8x32 [16 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 16; ++block_index) {

            block_offset = (quad_index * 16);
            idx = tab8x32[block_index] + block_offset;
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            block_shift_x = ((block_index & 0x07) << 3) + x_offset;
            block_shift_y = ((block_index >> 3) << 5) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                8,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x32[idx],
                &context_ptr->p_best_mv8x32[idx],
                &context_ptr->psub_pel_direction8x32[idx]);

        }
    }

    // 16x16 [16 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 16; ++block_index) {



            block_offset = (quad_index * 16);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab16x16[block_index] + block_offset;
            block_shift_x = ((block_index & 0x03) << 4) + x_offset;
            block_shift_y = ((block_index >> 2) << 4) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                16,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x16[idx],
                &context_ptr->p_best_mv16x16[idx],
                &context_ptr->psub_pel_direction16x16[idx]);
        }
    }

    // 32x16 [8 partitions]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 8; ++block_index) {

            block_offset = (quad_index * 8);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab32x16[block_index] + block_offset;
            block_shift_x = ((block_index & 0x01) << 5) + x_offset;
            block_shift_y = ((block_index >> 1) << 4) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                32,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x16[idx],
                &context_ptr->p_best_mv32x16[idx],
                &context_ptr->psub_pel_direction32x16[idx]);

        }
    }

    // 16x32 [8 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 8; ++block_index) {

            block_offset = (quad_index * 8);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab16x32[block_index] + block_offset;
            block_shift_x = ((block_index & 0x03) << 4) + x_offset;
            block_shift_y = ((block_index >> 2) << 5) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                16,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x32[idx],
                &context_ptr->p_best_mv16x32[idx],
                &context_ptr->psub_pel_direction16x32[idx]);

        }
    }

    // 32x32 [4 partitions]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 4; ++block_index) {

            block_offset = (quad_index * 4);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab32x32[block_index] + block_offset;
            block_shift_x = ((block_index & 0x01) << 5) + x_offset;
            block_shift_y = ((block_index >> 1) << 5) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                32,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x32[idx],
                &context_ptr->p_best_mv32x32[idx],
                &context_ptr->psub_pel_direction32x32[idx]);
        }

    }

    // 64x32 [2 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 2; ++block_index) {

            block_offset = (quad_index * 2);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab64x32[block_index] + block_offset;
            block_shift_x = x_offset;
            block_shift_y = (block_index << 5) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                64,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad64x32[idx],
                &context_ptr->p_best_mv64x32[idx],
                &context_ptr->psub_pel_direction64x32[idx]);

        }
    }

    // 32x64 [2 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 2; ++block_index) {

            block_offset = (quad_index * 2);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            idx = tab32x64[block_index] + block_offset;
            block_shift_x = (block_index << 5) + x_offset;
            block_shift_y = y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                32,
                64,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,

                &context_ptr->p_best_sad32x64[idx],
                &context_ptr->p_best_mv32x64[idx],
                &context_ptr->psub_pel_direction32x64[idx]);

        }

    }

    // 64x64 [1 partition]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {
        idx = quad_index;
        x_offset = (quad_index & 0x01) << 6;
        y_offset = (quad_index >> 1) << 6;
        block_shift_x = x_offset;
        block_shift_y = y_offset;

        block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

        posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
        posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
        posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

        in_loop_me_halfpel_refinement_block(
            sequence_control_set_ptr,
            context_ptr,
            block_index_in_sb_buffer,
            &(pos_b_buffer[posb_buffer_index]),
            &(pos_h_buffer[posh_buffer_index]),
            &(pos_j_buffer[posj_buffer_index]),
            64,
            64,
            x_search_area_origin,
            y_search_area_origin,
            asm_type,
            &context_ptr->p_best_sad64x64[idx],
            &context_ptr->p_best_mv64x64[idx],
            &context_ptr->psub_pel_direction64x64[idx]);
    }


#if NO_SUBPEL_FOR_128X128
    if (0) {
#else
    if (context_ptr->sb_size == BLOCK_128X128) {
#endif
        // 128x64 [2 partitions]
        for (block_index = 0; block_index < 2; ++block_index) {

            block_shift_x = 0;
            block_shift_y = block_index << 6;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                128,
                64,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad128x64[block_index],
                &context_ptr->p_best_mv128x64[block_index],
                &context_ptr->psub_pel_direction128x64[block_index]);
        }

        // 64x128 [2 partitions]
        for (block_index = 0; block_index < 2; ++block_index) {


            block_shift_x = block_index << 6;
            block_shift_y = 0;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                64,
                128,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad64x128[block_index],
                &context_ptr->p_best_mv64x128[block_index],
                &context_ptr->psub_pel_direction64x128[block_index]);
        }

        // 128x128 [1 partition]
        {
            block_index = 0;
            block_shift_x = 0;
            block_shift_y = 0;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            posb_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index = block_shift_x + block_shift_y * context_ptr->interpolated_stride;

            in_loop_me_halfpel_refinement_block(
                sequence_control_set_ptr,
                context_ptr,
                block_index_in_sb_buffer,
                &(pos_b_buffer[posb_buffer_index]),
                &(pos_h_buffer[posh_buffer_index]),
                &(pos_j_buffer[posj_buffer_index]),
                128,
                128,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad128x128[block_index],
                &context_ptr->p_best_mv128x128[block_index],
                &context_ptr->psub_pel_direction128x128);
        }
    }
    return;
    }

/***************************************************************
* in_loop_me_quarterpel_refinement_on_the_fly_block
*   performs Quarter Pel refinement for each block
***************************************************************/
static void in_loop_me_quarterpel_refinement_on_the_fly_block(
    SsMeContext_t         *context_ptr,                      // [IN] ME context Ptr, used to get SB Ptr
    uint32_t                 block_index_in_sb_buffer,                // [IN] PU origin, used to point to source samples
    uint8_t                **buf1,                            // [IN]
    uint32_t                *buf1Stride,
    uint8_t                **buf2,                            // [IN]
    uint32_t                *buf2Stride,
    uint32_t                 pu_width,                         // [IN]  PU width
    uint32_t                 pu_height,                        // [IN]  PU height
    int16_t                 x_search_area_origin,               // [IN] search area origin in the horizontal direction, used to point to reference samples
    int16_t                 y_search_area_origin,               // [IN] search area origin in the vertical direction, used to point to reference samples
    EbAsm                 asm_type,
    uint32_t                *pBestSad,
    uint32_t                *pBestMV,
    uint8_t                  sub_pel_direction)
{

    int16_t x_mv = _MVXT(*pBestMV);
    int16_t y_mv = _MVYT(*pBestMV);

    int16_t xSearchIndex = ((x_mv + 2) >> 2) - x_search_area_origin;
    int16_t ySearchIndex = ((y_mv + 2) >> 2) - y_search_area_origin;

    uint64_t dist;

    EbBool validTL, validT, validTR, validR, validBR, validB, validBL, validL;

    int16_t xMvQuarter[8];
    int16_t yMvQuarter[8];
    int32_t searchRegionIndex1 = 0;
    int32_t searchRegionIndex2 = 0;

    if ((y_mv & 2) + ((x_mv & 2) >> 1)) {

        validTL = (EbBool)(sub_pel_direction == RIGHT_POSITION || sub_pel_direction == BOTTOM_RIGHT_POSITION || sub_pel_direction == BOTTOM_POSITION);
        validT = (EbBool)(sub_pel_direction == BOTTOM_RIGHT_POSITION || sub_pel_direction == BOTTOM_POSITION || sub_pel_direction == BOTTOM_LEFT_POSITION);
        validTR = (EbBool)(sub_pel_direction == BOTTOM_POSITION || sub_pel_direction == BOTTOM_LEFT_POSITION || sub_pel_direction == LEFT_POSITION);
        validR = (EbBool)(sub_pel_direction == BOTTOM_LEFT_POSITION || sub_pel_direction == LEFT_POSITION || sub_pel_direction == TOP_LEFT_POSITION);
        validBR = (EbBool)(sub_pel_direction == LEFT_POSITION || sub_pel_direction == TOP_LEFT_POSITION || sub_pel_direction == TOP_POSITION);
        validB = (EbBool)(sub_pel_direction == TOP_LEFT_POSITION || sub_pel_direction == TOP_POSITION || sub_pel_direction == TOP_RIGHT_POSITION);
        validBL = (EbBool)(sub_pel_direction == TOP_POSITION || sub_pel_direction == TOP_RIGHT_POSITION || sub_pel_direction == RIGHT_POSITION);
        validL = (EbBool)(sub_pel_direction == TOP_RIGHT_POSITION || sub_pel_direction == RIGHT_POSITION || sub_pel_direction == BOTTOM_RIGHT_POSITION);

    }
    else {

        validTL = (EbBool)(sub_pel_direction == LEFT_POSITION || sub_pel_direction == TOP_LEFT_POSITION || sub_pel_direction == TOP_POSITION);
        validT = (EbBool)(sub_pel_direction == TOP_LEFT_POSITION || sub_pel_direction == TOP_POSITION || sub_pel_direction == TOP_RIGHT_POSITION);
        validTR = (EbBool)(sub_pel_direction == TOP_POSITION || sub_pel_direction == TOP_RIGHT_POSITION || sub_pel_direction == RIGHT_POSITION);
        validR = (EbBool)(sub_pel_direction == TOP_RIGHT_POSITION || sub_pel_direction == RIGHT_POSITION || sub_pel_direction == BOTTOM_RIGHT_POSITION);
        validBR = (EbBool)(sub_pel_direction == RIGHT_POSITION || sub_pel_direction == BOTTOM_RIGHT_POSITION || sub_pel_direction == BOTTOM_POSITION);
        validB = (EbBool)(sub_pel_direction == BOTTOM_RIGHT_POSITION || sub_pel_direction == BOTTOM_POSITION || sub_pel_direction == BOTTOM_LEFT_POSITION);
        validBL = (EbBool)(sub_pel_direction == BOTTOM_POSITION || sub_pel_direction == BOTTOM_LEFT_POSITION || sub_pel_direction == LEFT_POSITION);
        validL = (EbBool)(sub_pel_direction == BOTTOM_LEFT_POSITION || sub_pel_direction == LEFT_POSITION || sub_pel_direction == TOP_LEFT_POSITION);
    }

    xMvQuarter[0] = x_mv - 1; // L  position
    xMvQuarter[1] = x_mv + 1; // R  position
    xMvQuarter[2] = x_mv;     // T  position
    xMvQuarter[3] = x_mv;     // B  position
    xMvQuarter[4] = x_mv - 1; // TL position
    xMvQuarter[5] = x_mv + 1; // TR position
    xMvQuarter[6] = x_mv + 1; // BR position
    xMvQuarter[7] = x_mv - 1; // BL position

    yMvQuarter[0] = y_mv;     // L  position
    yMvQuarter[1] = y_mv;     // R  position
    yMvQuarter[2] = y_mv - 1; // T  position
    yMvQuarter[3] = y_mv + 1; // B  position
    yMvQuarter[4] = y_mv - 1; // TL position
    yMvQuarter[5] = y_mv - 1; // TR position
    yMvQuarter[6] = y_mv + 1; // BR position
    yMvQuarter[7] = y_mv + 1; // BL position


    // L position
    if (validL) {

        searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[0] * (int32_t)ySearchIndex;
        searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[0] * (int32_t)ySearchIndex;

#if USE_INLOOP_ME_FULL_SAD
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride, buf1[0] + searchRegionIndex1, buf1Stride[0], buf2[0] + searchRegionIndex2, buf2Stride[0], pu_height, pu_width);
#else
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, buf1[0] + searchRegionIndex1, buf1Stride[0] << 1, buf2[0] + searchRegionIndex2, buf2Stride[0] << 1, pu_height >> 1, pu_width);

        dist = dist << 1;
#endif

        if (dist < *pBestSad) {
            *pBestSad = (uint32_t)dist;
            *pBestMV = ((uint16_t)yMvQuarter[0] << 16) | ((uint16_t)xMvQuarter[0]);
        }
    }

    // R positions
    if (validR) {

        searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[1] * (int32_t)ySearchIndex;
        searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[1] * (int32_t)ySearchIndex;
#if USE_INLOOP_ME_FULL_SAD
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride, buf1[1] + searchRegionIndex1, buf1Stride[1], buf2[1] + searchRegionIndex2, buf2Stride[1], pu_height, pu_width);
#else
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, buf1[1] + searchRegionIndex1, buf1Stride[1] << 1, buf2[1] + searchRegionIndex2, buf2Stride[1] << 1, pu_height >> 1, pu_width);
        dist = dist << 1;
#endif


        if (dist < *pBestSad) {
            *pBestSad = (uint32_t)dist;
            *pBestMV = ((uint16_t)yMvQuarter[1] << 16) | ((uint16_t)xMvQuarter[1]);
        }
    }

    // T position
    if (validT) {

        searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[2] * (int32_t)ySearchIndex;
        searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[2] * (int32_t)ySearchIndex;

#if USE_INLOOP_ME_FULL_SAD
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride, buf1[2] + searchRegionIndex1, buf1Stride[2], buf2[2] + searchRegionIndex2, buf2Stride[2], pu_height, pu_width);
#else
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, buf1[2] + searchRegionIndex1, buf1Stride[2] << 1, buf2[2] + searchRegionIndex2, buf2Stride[2] << 1, pu_height >> 1, pu_width);
        dist = dist << 1;
#endif


        if (dist < *pBestSad) {
            *pBestSad = (uint32_t)dist;
            *pBestMV = ((uint16_t)yMvQuarter[2] << 16) | ((uint16_t)xMvQuarter[2]);
        }
    }

    // B position
    if (validB) {

        searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[3] * (int32_t)ySearchIndex;
        searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[3] * (int32_t)ySearchIndex;

#if USE_INLOOP_ME_FULL_SAD
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride, buf1[3] + searchRegionIndex1, buf1Stride[3], buf2[3] + searchRegionIndex2, buf2Stride[3], pu_height, pu_width);
#else
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, buf1[3] + searchRegionIndex1, buf1Stride[3] << 1, buf2[3] + searchRegionIndex2, buf2Stride[3] << 1, pu_height >> 1, pu_width);
        dist = dist << 1;
#endif


        if (dist < *pBestSad) {
            *pBestSad = (uint32_t)dist;
            *pBestMV = ((uint16_t)yMvQuarter[3] << 16) | ((uint16_t)xMvQuarter[3]);
        }
    }

    //TL position
    if (validTL) {

        searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[4] * (int32_t)ySearchIndex;
        searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[4] * (int32_t)ySearchIndex;

#if USE_INLOOP_ME_FULL_SAD
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride, buf1[4] + searchRegionIndex1, buf1Stride[4], buf2[4] + searchRegionIndex2, buf2Stride[4], pu_height, pu_width);
#else
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, buf1[4] + searchRegionIndex1, buf1Stride[4] << 1, buf2[4] + searchRegionIndex2, buf2Stride[4] << 1, pu_height >> 1, pu_width);
        dist = dist << 1;
#endif

        if (dist < *pBestSad) {
            *pBestSad = (uint32_t)dist;
            *pBestMV = ((uint16_t)yMvQuarter[4] << 16) | ((uint16_t)xMvQuarter[4]);
        }
    }

    //TR position
    if (validTR) {

        searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[5] * (int32_t)ySearchIndex;
        searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[5] * (int32_t)ySearchIndex;

#if USE_INLOOP_ME_FULL_SAD
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride, buf1[5] + searchRegionIndex1, buf1Stride[5], buf2[5] + searchRegionIndex2, buf2Stride[5], pu_height, pu_width);
#else
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, buf1[5] + searchRegionIndex1, buf1Stride[5] << 1, buf2[5] + searchRegionIndex2, buf2Stride[5] << 1, pu_height >> 1, pu_width);
        dist = dist << 1;
#endif

        if (dist < *pBestSad) {
            *pBestSad = (uint32_t)dist;
            *pBestMV = ((uint16_t)yMvQuarter[5] << 16) | ((uint16_t)xMvQuarter[5]);
        }
    }

    //BR position
    if (validBR) {

        searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[6] * (int32_t)ySearchIndex;
        searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[6] * (int32_t)ySearchIndex;

#if USE_INLOOP_ME_FULL_SAD
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride, buf1[6] + searchRegionIndex1, buf1Stride[6], buf2[6] + searchRegionIndex2, buf2Stride[6], pu_height, pu_width);
#else
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, buf1[6] + searchRegionIndex1, buf1Stride[6] << 1, buf2[6] + searchRegionIndex2, buf2Stride[6] << 1, pu_height >> 1, pu_width);
        dist = dist << 1;
#endif

        if (dist < *pBestSad) {
            *pBestSad = (uint32_t)dist;
            *pBestMV = ((uint16_t)yMvQuarter[6] << 16) | ((uint16_t)xMvQuarter[6]);
        }
    }

    //BL position
    if (validBL) {

        searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[7] * (int32_t)ySearchIndex;
        searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[7] * (int32_t)ySearchIndex;

#if USE_INLOOP_ME_FULL_SAD
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride, buf1[7] + searchRegionIndex1, buf1Stride[7], buf2[7] + searchRegionIndex2, buf2Stride[7], pu_height, pu_width);
#else
        dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[block_index_in_sb_buffer]), context_ptr->sb_src_stride << 1, buf1[7] + searchRegionIndex1, buf1Stride[7] << 1, buf2[7] + searchRegionIndex2, buf2Stride[7] << 1, pu_height >> 1, pu_width);
        dist = dist << 1;
#endif

        if (dist < *pBestSad) {
            *pBestSad = (uint32_t)dist;
            *pBestMV = ((uint16_t)yMvQuarter[7] << 16) | ((uint16_t)xMvQuarter[7]);
        }
    }


    return;
}


/***************************************************************
* set_quarterpel_refinement_inputs_on_the_fly_block
*   determine the 2 half pel buffers to perform the averaging
*   for Quarter Pel Refinement
***************************************************************/
static void set_quarterpel_refinement_inputs_on_the_fly_block(
    uint8_t   *pos_Full,   //[IN] points to A
    uint32_t   FullStride, //[IN]
    uint8_t   *pos_b,     //[IN] points to b
    uint8_t   *pos_h,     //[IN] points to h
    uint8_t   *pos_j,     //[IN] points to j
    uint32_t   Stride,    //[IN]
    int16_t   x_mv,        //[IN]
    int16_t   y_mv,        //[IN]
    uint8_t   **buf1,       //[OUT]
    uint32_t  *buf1Stride, //[OUT]
    uint8_t   **buf2,       //[OUT]
    uint32_t  *buf2Stride  //[OUT]
)
{

    uint32_t  quarterPelRefinementMethod = (y_mv & 2) + ((x_mv & 2) >> 1);

    //for each one of the 8 postions, we need to determine the 2 half pel buffers to  do averaging

    //     A    a    b    c
    //     d    e    f    g
    //     h    i    j    k
    //     n    p    q    r

    switch (quarterPelRefinementMethod) {

    case EB_QUARTER_IN_FULL:

        /*c=b+A*/ buf1[0] = pos_b;                     buf1Stride[0] = Stride;        buf2[0] = pos_Full;             buf2Stride[0] = FullStride;
        /*a=A+b*/ buf1[1] = pos_Full;                  buf1Stride[1] = FullStride;    buf2[1] = pos_b + 1;             buf2Stride[1] = Stride;
        /*n=h+A*/ buf1[2] = pos_h;                      buf1Stride[2] = Stride;        buf2[2] = pos_Full;              buf2Stride[2] = FullStride;
        /*d=A+h*/ buf1[3] = pos_Full;                   buf1Stride[3] = FullStride;    buf2[3] = pos_h + Stride;        buf2Stride[3] = Stride;
        /*r=b+h*/ buf1[4] = pos_b;                      buf1Stride[4] = Stride;        buf2[4] = pos_h;                 buf2Stride[4] = Stride;
        /*p=h+b*/ buf1[5] = pos_h;                      buf1Stride[5] = Stride;        buf2[5] = pos_b + 1;             buf2Stride[5] = Stride;
        /*e=h+b*/ buf1[6] = pos_h + Stride;             buf1Stride[6] = Stride;        buf2[6] = pos_b + 1;             buf2Stride[6] = Stride;
        /*g=b+h*/ buf1[7] = pos_b;                      buf1Stride[7] = Stride;        buf2[7] = pos_h + Stride;        buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_HORIZONTAL:

        /*a=A+b*/ buf1[0] = pos_Full - 1;               buf1Stride[0] = FullStride;    buf2[0] = pos_b;                buf2Stride[0] = Stride;
        /*c=b+A*/ buf1[1] = pos_b;                     buf1Stride[1] = Stride;        buf2[1] = pos_Full;             buf2Stride[1] = FullStride;
        /*q=j+b*/ buf1[2] = pos_j;                     buf1Stride[2] = Stride;        buf2[2] = pos_b;                buf2Stride[2] = Stride;
        /*f=b+j*/ buf1[3] = pos_b;                     buf1Stride[3] = Stride;        buf2[3] = pos_j + Stride;        buf2Stride[3] = Stride;
        /*p=h+b*/ buf1[4] = pos_h - 1;                  buf1Stride[4] = Stride;        buf2[4] = pos_b;                buf2Stride[4] = Stride;
        /*r=b+h*/ buf1[5] = pos_b;                     buf1Stride[5] = Stride;        buf2[5] = pos_h;                buf2Stride[5] = Stride;
        /*g=b+h*/ buf1[6] = pos_b;                     buf1Stride[6] = Stride;        buf2[6] = pos_h + Stride;        buf2Stride[6] = Stride;
        /*e=h+b*/ buf1[7] = pos_h - 1 + Stride;         buf1Stride[7] = Stride;        buf2[7] = pos_b;                buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_VERTICAL:

        /*k=j+h*/buf1[0] = pos_j;                      buf1Stride[0] = Stride;        buf2[0] = pos_h;                 buf2Stride[0] = Stride;
        /*i=h+j*/buf1[1] = pos_h;                      buf1Stride[1] = Stride;        buf2[1] = pos_j + 1;              buf2Stride[1] = Stride;
        /*d=A+h*/buf1[2] = pos_Full - FullStride;      buf1Stride[2] = FullStride;    buf2[2] = pos_h;                  buf2Stride[2] = Stride;
        /*n=h+A*/buf1[3] = pos_h;                       buf1Stride[3] = Stride;        buf2[3] = pos_Full;               buf2Stride[3] = FullStride;
        /*g=b+h*/buf1[4] = pos_b - Stride;              buf1Stride[4] = Stride;        buf2[4] = pos_h;                  buf2Stride[4] = Stride;
        /*e=h+b*/buf1[5] = pos_h;                      buf1Stride[5] = Stride;        buf2[5] = pos_b + 1 - Stride;     buf2Stride[5] = Stride;
        /*p=h+b*/buf1[6] = pos_h;                      buf1Stride[6] = Stride;        buf2[6] = pos_b + 1;              buf2Stride[6] = Stride;
        /*r=b+h*/buf1[7] = pos_b;                      buf1Stride[7] = Stride;        buf2[7] = pos_h;                 buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_DIAGONAL:

        /*i=h+j*/buf1[0] = pos_h - 1;                   buf1Stride[0] = Stride;        buf2[0] = pos_j;                  buf2Stride[0] = Stride;
        /*k=j+h*/buf1[1] = pos_j;                       buf1Stride[1] = Stride;        buf2[1] = pos_h;                  buf2Stride[1] = Stride;
        /*f=b+j*/buf1[2] = pos_b - Stride;              buf1Stride[2] = Stride;        buf2[2] = pos_j;                  buf2Stride[2] = Stride;
        /*q=j+b*/buf1[3] = pos_j;                       buf1Stride[3] = Stride;        buf2[3] = pos_b;                  buf2Stride[3] = Stride;
        /*e=h+b*/buf1[4] = pos_h - 1;                   buf1Stride[4] = Stride;        buf2[4] = pos_b - Stride;         buf2Stride[4] = Stride;
        /*g=b+h*/buf1[5] = pos_b - Stride;              buf1Stride[5] = Stride;        buf2[5] = pos_h;                  buf2Stride[5] = Stride;
        /*r=b+h*/buf1[6] = pos_b;                      buf1Stride[6] = Stride;        buf2[6] = pos_h;                  buf2Stride[6] = Stride;
        /*p=h+b*/buf1[7] = pos_h - 1;                   buf1Stride[7] = Stride;        buf2[7] = pos_b;                  buf2Stride[7] = Stride;

        break;

    default:
        break;

    }

    return;
}

/***************************************************************
* in_loop_me_quarterpel_search_sblock
*   perform the quarter-pel refinement for the whole super-block
***************************************************************/
static void in_loop_me_quarterpel_search_sblock(
    SsMeContext_t                *context_ptr,                     //[IN/OUT]  ME context Ptr, used to get/update ME results
    uint8_t                        *pos_Full,                       //[IN]
    uint32_t                        full_stride,                      //[IN]
    uint8_t                        *pos_b,                          //[IN]
    uint8_t                        *pos_h,                          //[IN]
    uint8_t                        *pos_j,                          //[IN]
    int16_t                        x_search_area_origin,            //[IN] search area origin in the horizontal direction, used to point to reference samples
    int16_t                        y_search_area_origin,               //[IN] search area origin in the vertical direction, used to point to reference samples
    EbAsm                        asm_type)
{
    uint32_t  block_index;

    uint32_t  block_shift_x;
    uint32_t  block_shift_y;

    uint32_t  block_index_in_sb_buffer;

    //for each one of the 8 positions, we need to determine the 2 buffers to  do averaging
    uint8_t  *buf1[8];
    uint8_t  *buf2[8];

    uint32_t  buf1Stride[8];
    uint32_t  buf2Stride[8];

    int16_t  x_mv, y_mv;
    uint32_t  nidx;

    uint32_t quad_index = 0;
    uint32_t block_offset = 0;
    uint32_t x_offset = 0;
    uint32_t y_offset = 0;
    uint32_t number_of_sb_quad = context_ptr->sb_size == BLOCK_128X128 ? 4 : 1;

    // 4x4   [256 partitions]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 256; ++block_index) {

            block_offset = (quad_index * 256);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab4x4[block_index] + block_offset;
            block_shift_x = ((block_index & 0xf) << 2) + x_offset;
            block_shift_y = ((block_index >> 4) << 2) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv4x4[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv4x4[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                4, 4,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad4x4[nidx],
                &context_ptr->p_best_mv4x4[nidx],
                context_ptr->psub_pel_direction4x4[nidx]);
        }
    }

    // 8x4   [128 8x4 blocks]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 128; ++block_index) {


            block_offset = (quad_index * 128);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab8x4[block_index] + block_offset;
            block_shift_x = ((block_index & 0x07) << 3) + x_offset;
            block_shift_y = ((block_index >> 3) << 2) + y_offset;


            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;


            x_mv = _MVXT(context_ptr->p_best_mv8x4[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x4[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                8, 4,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x4[nidx],
                &context_ptr->p_best_mv8x4[nidx],
                context_ptr->psub_pel_direction8x4[nidx]);

        }

    }

    // 4x8   [128 4x8 blocks]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 128; ++block_index) {

            block_offset = (quad_index * 128);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab4x8[block_index] + block_offset;
            block_shift_x = ((block_index & 0xf) << 2) + x_offset;
            block_shift_y = ((block_index >> 4) << 3) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;


            x_mv = _MVXT(context_ptr->p_best_mv4x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv4x8[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                4, 8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad4x8[nidx],
                &context_ptr->p_best_mv4x8[nidx],
                context_ptr->psub_pel_direction4x8[nidx]);

        }

    }

    // 8x8   [64 8x8 blocks]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 64; ++block_index) {


            block_offset = (quad_index * 64);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab8x8[block_index] + block_offset;
            block_shift_x = ((block_index & 0x07) << 3) + x_offset;
            block_shift_y = ((block_index >> 3) << 3) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv8x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x8[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                8, 8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x8[nidx],
                &context_ptr->p_best_mv8x8[nidx],
                context_ptr->psub_pel_direction8x8[nidx]);
        }
    }

    // 16x8 [32 partitions]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 32; ++block_index) {


            block_offset = (quad_index * 32);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab16x8[block_index] + block_offset;
            block_shift_x = ((block_index & 0x03) << 4) + x_offset;
            block_shift_y = ((block_index >> 2) << 3) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv16x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x8[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                16, 8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x8[nidx],
                &context_ptr->p_best_mv16x8[nidx],
                context_ptr->psub_pel_direction16x8[nidx]);
        }

    }

    // 8x16 [32 partitions]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 32; ++block_index) {

            block_offset = (quad_index * 32);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab8x16[block_index] + block_offset;
            block_shift_x = ((block_index & 0x07) << 3) + x_offset;
            block_shift_y = ((block_index >> 3) << 4) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv8x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x16[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                8, 16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x16[nidx],
                &context_ptr->p_best_mv8x16[nidx],
                context_ptr->psub_pel_direction8x16[nidx]);

        }

    }

    // 32x8 [16 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {
        for (block_index = 0; block_index < 16; ++block_index) {

            block_offset = (quad_index * 16);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab32x8[block_index] + block_offset;
            block_shift_x = ((block_index & 0x01) << 5) + x_offset;
            block_shift_y = ((block_index >> 1) << 3) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv32x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv32x8[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                32, 8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x8[nidx],
                &context_ptr->p_best_mv32x8[nidx],
                context_ptr->psub_pel_direction32x8[nidx]);
        }
    }

    // 8x32 [16 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {
        for (block_index = 0; block_index < 16; ++block_index) {

            block_offset = (quad_index * 16);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab8x32[block_index] + block_offset;
            block_shift_x = ((block_index & 0x07) << 3) + x_offset;
            block_shift_y = ((block_index >> 3) << 5) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv8x32[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x32[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                8, 32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x32[nidx],
                &context_ptr->p_best_mv8x32[nidx],
                context_ptr->psub_pel_direction8x32[nidx]);
        }
    }

    // 16x16 [16 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 16; ++block_index) {

            block_offset = (quad_index * 16);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab16x16[block_index] + block_offset;
            block_shift_x = ((block_index & 0x03) << 4) + x_offset;
            block_shift_y = ((block_index >> 2) << 4) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv16x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x16[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                16, 16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x16[nidx],
                &context_ptr->p_best_mv16x16[nidx],
                context_ptr->psub_pel_direction16x16[nidx]);
        }
    }

    // 32x16 [8 partitions]

    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {
        for (block_index = 0; block_index < 8; ++block_index) {


            block_offset = (quad_index * 8);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab32x16[block_index] + block_offset;
            block_shift_x = ((block_index & 0x01) << 5) + x_offset;
            block_shift_y = ((block_index >> 1) << 4) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv32x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv32x16[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                32, 16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x16[nidx],
                &context_ptr->p_best_mv32x16[nidx],
                context_ptr->psub_pel_direction32x16[nidx]);
        }
    }

    // 16x32 [8 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 8; ++block_index) {

            block_offset = (quad_index * 8);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab16x32[block_index] + block_offset;
            block_shift_x = ((block_index & 0x03) << 4) + x_offset;
            block_shift_y = ((block_index >> 2) << 5) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv16x32[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x32[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                16, 32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x32[nidx],
                &context_ptr->p_best_mv16x32[nidx],
                context_ptr->psub_pel_direction16x32[nidx]);

        }
    }

    // 32x32 [4 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 4; ++block_index) {

            block_offset = (quad_index * 4);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab32x32[block_index] + block_offset;
            block_shift_x = ((block_index & 0x01) << 5) + x_offset;
            block_shift_y = ((block_index >> 1)) + y_offset;
            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv32x32[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv32x32[nidx]);


            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                32, 32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,

                &context_ptr->p_best_sad32x32[nidx],
                &context_ptr->p_best_mv32x32[nidx],
                context_ptr->psub_pel_direction32x32[nidx]);

        }
    }

    // 64x32 [2 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {
        for (block_index = 0; block_index < 2; ++block_index) {

            block_offset = (quad_index * 2);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab64x32[block_index] + block_offset;
            block_shift_x = x_offset;
            block_shift_y = (block_index << 5) + y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv64x32[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv64x32[nidx]);


            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                64, 32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,

                &context_ptr->p_best_sad64x32[nidx],
                &context_ptr->p_best_mv64x32[nidx],
                context_ptr->psub_pel_direction64x32[nidx]);
        }
    }

    // 32x64 [2 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {

        for (block_index = 0; block_index < 2; ++block_index) {

            block_offset = (quad_index * 2);
            x_offset = (quad_index & 0x01) << 6;
            y_offset = (quad_index >> 1) << 6;
            nidx = tab32x64[block_index] + block_offset;
            block_shift_x = (block_index << 5) + x_offset;
            block_shift_y = y_offset;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv32x64[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv32x64[nidx]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            buf1[0] = buf1[0] + block_shift_x + block_shift_y * buf1Stride[0];              buf2[0] = buf2[0] + block_shift_x + block_shift_y * buf2Stride[0];
            buf1[1] = buf1[1] + block_shift_x + block_shift_y * buf1Stride[1];              buf2[1] = buf2[1] + block_shift_x + block_shift_y * buf2Stride[1];
            buf1[2] = buf1[2] + block_shift_x + block_shift_y * buf1Stride[2];              buf2[2] = buf2[2] + block_shift_x + block_shift_y * buf2Stride[2];
            buf1[3] = buf1[3] + block_shift_x + block_shift_y * buf1Stride[3];              buf2[3] = buf2[3] + block_shift_x + block_shift_y * buf2Stride[3];
            buf1[4] = buf1[4] + block_shift_x + block_shift_y * buf1Stride[4];              buf2[4] = buf2[4] + block_shift_x + block_shift_y * buf2Stride[4];
            buf1[5] = buf1[5] + block_shift_x + block_shift_y * buf1Stride[5];              buf2[5] = buf2[5] + block_shift_x + block_shift_y * buf2Stride[5];
            buf1[6] = buf1[6] + block_shift_x + block_shift_y * buf1Stride[6];              buf2[6] = buf2[6] + block_shift_x + block_shift_y * buf2Stride[6];
            buf1[7] = buf1[7] + block_shift_x + block_shift_y * buf1Stride[7];              buf2[7] = buf2[7] + block_shift_x + block_shift_y * buf2Stride[7];

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                32, 64,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x64[nidx],
                &context_ptr->p_best_mv32x64[nidx],
                context_ptr->psub_pel_direction32x64[nidx]);
        }

    }

    // 64x64 [1 partitions]
    for (quad_index = 0; quad_index < number_of_sb_quad; quad_index++) {


        block_index = 0;

        block_offset = quad_index;
        x_offset = (quad_index & 0x01) << 6;
        y_offset = (quad_index >> 1) << 6;
        nidx = block_offset;
        block_shift_x = x_offset;
        block_shift_y = y_offset;

        block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

        x_mv = _MVXT(context_ptr->p_best_mv64x64[nidx]);
        y_mv = _MVYT(context_ptr->p_best_mv64x64[nidx]);


        set_quarterpel_refinement_inputs_on_the_fly_block(
            pos_Full,
            full_stride,
            pos_b,
            pos_h,
            pos_j,
            context_ptr->interpolated_stride,
            x_mv,
            y_mv,
            buf1, buf1Stride,
            buf2, buf2Stride);

        in_loop_me_quarterpel_refinement_on_the_fly_block(
            context_ptr,
            block_index_in_sb_buffer,
            buf1, buf1Stride,
            buf2, buf2Stride,
            64, 64,
            x_search_area_origin,
            y_search_area_origin,
            asm_type,
            &context_ptr->p_best_sad64x64[nidx],
            &context_ptr->p_best_mv64x64[nidx],
            context_ptr->psub_pel_direction64x64[nidx]);

    }


#if NO_SUBPEL_FOR_128X128
    if (0) {
#else
    if (context_ptr->sb_size == BLOCK_128X128) {
#endif
        // 128x64 [2 partitions]
        for (block_index = 0; block_index < 2; ++block_index) {

            block_index = 0;

            block_shift_x = 0;
            block_shift_y = block_index << 6;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv128x64[block_index]);
            y_mv = _MVYT(context_ptr->p_best_mv128x64[block_index]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                128, 64,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad128x64[block_index],
                &context_ptr->p_best_mv128x64[block_index],
                context_ptr->psub_pel_direction128x64[block_index]);

        }
        // 64x128 [2 partitions]
        for (block_index = 0; block_index < 2; ++block_index) {

            block_index = 0;

            block_shift_x = block_index << 6;
            block_shift_y = 0;

            block_index_in_sb_buffer = block_shift_x + block_shift_y * context_ptr->sb_src_stride;

            x_mv = _MVXT(context_ptr->p_best_mv64x128[block_index]);
            y_mv = _MVYT(context_ptr->p_best_mv64x128[block_index]);

            set_quarterpel_refinement_inputs_on_the_fly_block(
                pos_Full,
                full_stride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1, buf1Stride,
                buf2, buf2Stride);

            in_loop_me_quarterpel_refinement_on_the_fly_block(
                context_ptr,
                block_index_in_sb_buffer,
                buf1, buf1Stride,
                buf2, buf2Stride,
                64, 128,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad64x128[block_index],
                &context_ptr->p_best_mv64x128[block_index],
                context_ptr->psub_pel_direction64x128[block_index]);

        }
        // 128x128 [1 partitions]
        block_index = 0;

        block_shift_x = 0;
        block_shift_y = 0;

        block_index_in_sb_buffer = 0;

        x_mv = _MVXT(context_ptr->p_best_mv128x128[block_index]);
        y_mv = _MVYT(context_ptr->p_best_mv128x128[block_index]);

        set_quarterpel_refinement_inputs_on_the_fly_block(
            pos_Full,
            full_stride,
            pos_b,
            pos_h,
            pos_j,
            context_ptr->interpolated_stride,
            x_mv,
            y_mv,
            buf1, buf1Stride,
            buf2, buf2Stride);

        in_loop_me_quarterpel_refinement_on_the_fly_block(
            context_ptr,
            block_index_in_sb_buffer,
            buf1, buf1Stride,
            buf2, buf2Stride,
            128, 128,
            x_search_area_origin,
            y_search_area_origin,
            asm_type,
            &context_ptr->p_best_sad128x128[block_index],
            &context_ptr->p_best_mv128x128[block_index],
            context_ptr->psub_pel_direction128x128);
    }
    return;
    }

#define MAX_SEARCH_POINT_WIDTH  128
#define MAX_SEARCH_POINT_HEIGHT 128

#define MAX_TATAL_SEARCH_AREA_WIDTH        (MAX_SB_SIZE + MAX_SEARCH_POINT_WIDTH  + ME_FILTER_TAP)
#define MAX_TATAL_SEARCH_AREA_HEIGHT       (MAX_SB_SIZE + MAX_SEARCH_POINT_HEIGHT  + ME_FILTER_TAP)


#define MAX_SEARCH_AREA_SIZE     MAX_TATAL_SEARCH_AREA_WIDTH * MAX_TATAL_SEARCH_AREA_HEIGHT
/***************************************************************
* in_loop_motion_estimation_sblock
*  perform the full-pel serach for the whole super-block
*  on the reference reconstructed pictures
***************************************************************/
EB_EXTERN EbErrorType in_loop_motion_estimation_sblock(
    PictureControlSet_t         *picture_control_set_ptr,  // input parameter, Picture Control Set Ptr
    uint32_t                       sb_origin_x,            // input parameter, SB Origin X
    uint32_t                       sb_origin_y,            // input parameter, SB Origin X
    int16_t                       xMvL0,
    int16_t                       yMvL0,
    int16_t                       xMvL1,
    int16_t                       yMvL1,
    SsMeContext_t                 *context_ptr)           // input parameter, ME Context Ptr, used to store decimated/interpolated LCU/SR

{
    EbErrorType return_error = EB_ErrorNone;

    SequenceControlSet_t    *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;

    int16_t                  xTopLeftSearchRegion;
    int16_t                  yTopLeftSearchRegion;
    uint32_t                  searchRegionIndex;
    int16_t                  picture_width = (int16_t)((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr)->luma_width;
    int16_t                  picture_height = (int16_t)((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr)->luma_height;

    int16_t                  padWidth = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t                  padHeight = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t                  search_area_width;
    int16_t                  search_area_height;
    int16_t                  x_search_area_origin;
    int16_t                  y_search_area_origin;
    int16_t                  origin_x = (int16_t)sb_origin_x;
    int16_t                  origin_y = (int16_t)sb_origin_y;

    uint8_t                   refPicIndex = 0;
    // Final ME Search Center
    int16_t                  xSearchCenter = 0;
    int16_t                  ySearchCenter = 0;

    uint32_t                  numOfListToSearch;
    uint32_t                  listIndex;
    EbPictureBufferDesc_t  *refPicPtr;
    EbReferenceObject_t    *referenceObject;

    EbAsm                  asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;

    uint32_t                  number_of_sb_quad = sequence_control_set_ptr->sb_size == BLOCK_128X128 ? 4 : 1;
    context_ptr->sb_size = sequence_control_set_ptr->sb_size;
    context_ptr->sb_side = sequence_control_set_ptr->sb_size == BLOCK_128X128 ? 128 : 64;

    const uint32_t start_idx_8x8 = 256 * number_of_sb_quad;
    const uint32_t start_idx_16x16 = 320 * number_of_sb_quad;
    const uint32_t start_idx_32x32 = 336 * number_of_sb_quad;
    const uint32_t start_idx_64x64 = 340 * number_of_sb_quad;
    const uint32_t start_idx_8x4 = 341 * number_of_sb_quad;
    const uint32_t start_idx_4x8 = 469 * number_of_sb_quad;
    const uint32_t start_idx_4x16 = 597 * number_of_sb_quad;
    const uint32_t start_idx_16x4 = 661 * number_of_sb_quad;
    const uint32_t start_idx_16x8 = 725 * number_of_sb_quad;
    const uint32_t start_idx_8x16 = 757 * number_of_sb_quad;
    const uint32_t start_idx_32x8 = 789 * number_of_sb_quad;
    const uint32_t start_idx_8x32 = 805 * number_of_sb_quad;
    const uint32_t start_idx_32x16 = 821 * number_of_sb_quad;
    const uint32_t start_idx_16x32 = 829 * number_of_sb_quad;
    const uint32_t start_idx_64x16 = 837 * number_of_sb_quad;
    const uint32_t start_idx_16x64 = 841 * number_of_sb_quad;
    const uint32_t start_idx_64x32 = 845 * number_of_sb_quad;
    const uint32_t start_idx_32x64 = 847 * number_of_sb_quad;
    const uint32_t start_idx_128x64 = 849 * number_of_sb_quad;

#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
    context_ptr->fractionalSearchMethod = (picture_control_set_ptr->enc_mode >= ENC_M3) ? FULL_SAD_SEARCH : SSD_SEARCH;
#else
    context_ptr->fractionalSearchMethod = SUB_SAD_SEARCH;
#endif

#if M0_ME_SEARCH_BASE
#if ENCODER_MODE_CLEANUP
    numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
#else
    numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE || (picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->enc_mode > ENC_M1)) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
#endif
#else
    numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE) || (picture_control_set_ptr->temporal_layer_index == 0) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
#endif

    // Uni-Prediction motion estimation loop
    // List Loop
    for (listIndex = REF_LIST_0; listIndex <= numOfListToSearch; ++listIndex) {

        EbBool  is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
        referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[listIndex]->objectPtr;
        refPicPtr = is16bit ? (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit : (EbPictureBufferDesc_t*)referenceObject->referencePicture;
        search_area_width = (int16_t)MIN(context_ptr->search_area_width, 127);
        search_area_height = (int16_t)MIN(context_ptr->search_area_height, 127);
        xSearchCenter = listIndex == REF_LIST_0 ? xMvL0 : xMvL1;
        ySearchCenter = listIndex == REF_LIST_0 ? yMvL0 : yMvL1;

        x_search_area_origin = xSearchCenter - (search_area_width >> 1);
        y_search_area_origin = ySearchCenter - (search_area_height >> 1);

        // Correct the left edge of the Search Area if it is not on the reference Picture
        x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth) ?
            -padWidth - origin_x :
            x_search_area_origin;

        search_area_width = ((origin_x + x_search_area_origin) < -padWidth) ?
            search_area_width - (-padWidth - (origin_x + x_search_area_origin)) :
            search_area_width;

        // Correct the right edge of the Search Area if its not on the reference Picture
        x_search_area_origin = ((origin_x + x_search_area_origin) > picture_width - 1) ?
            x_search_area_origin - ((origin_x + x_search_area_origin) - (picture_width - 1)) :
            x_search_area_origin;

        // //check whether the needed search area is coverd by the reference picture and adjust its origin to satisfy the condition if not.
        if (sequence_control_set_ptr->sb_size == BLOCK_128X128) {

            int32_t righ_sa_pos_x = refPicPtr->origin_x + origin_x + x_search_area_origin + search_area_width + (context_ptr->sb_side - 1) + (ME_FILTER_TAP >> 1);
            int32_t righ_ref_pos_x = picture_width - 1 + (2 * refPicPtr->origin_x);

            x_search_area_origin = righ_sa_pos_x > righ_ref_pos_x ? x_search_area_origin - (righ_sa_pos_x - righ_ref_pos_x) : x_search_area_origin;

            int32_t bottom_sa_pos_x = refPicPtr->origin_y + origin_y + y_search_area_origin + search_area_height + (context_ptr->sb_side - 1) + (ME_FILTER_TAP >> 1);
            int32_t bottom_ref_pos_x = picture_height - 1 + (2 * refPicPtr->origin_y);

            y_search_area_origin = bottom_sa_pos_x > bottom_ref_pos_x ? y_search_area_origin - (bottom_sa_pos_x - bottom_ref_pos_x) : y_search_area_origin;
        }


        search_area_width = ((origin_x + x_search_area_origin + search_area_width) > picture_width) ?
            MAX(1, search_area_width - ((origin_x + x_search_area_origin + search_area_width) - picture_width)) :
            search_area_width;

        // Correct the top edge of the Search Area if it is not on the reference Picture
        y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight) ?
            -padHeight - origin_y :
            y_search_area_origin;

        search_area_height = ((origin_y + y_search_area_origin) < -padHeight) ?
            search_area_height - (-padHeight - (origin_y + y_search_area_origin)) :
            search_area_height;

        // Correct the bottom edge of the Search Area if its not on the reference Picture
        y_search_area_origin = ((origin_y + y_search_area_origin) > picture_height - 1) ?
            y_search_area_origin - ((origin_y + y_search_area_origin) - (picture_height - 1)) :
            y_search_area_origin;

        search_area_height = (origin_y + y_search_area_origin + search_area_height > picture_height) ?
            MAX(1, search_area_height - ((origin_y + y_search_area_origin + search_area_height) - picture_height)) :
            search_area_height;

        context_ptr->x_search_area_origin[listIndex][0] = x_search_area_origin;
        context_ptr->y_search_area_origin[listIndex][0] = y_search_area_origin;

        xTopLeftSearchRegion = (int16_t)(refPicPtr->origin_x + sb_origin_x) - (ME_FILTER_TAP >> 1) + x_search_area_origin;
        yTopLeftSearchRegion = (int16_t)(refPicPtr->origin_y + sb_origin_y) - (ME_FILTER_TAP >> 1) + y_search_area_origin;
        searchRegionIndex = (xTopLeftSearchRegion)+(yTopLeftSearchRegion)* refPicPtr->strideY;

        // Umpack the reference for 16bit reference picture.
        if (is16bit) {

            uint16_t *ptr16 = (uint16_t *)refPicPtr->bufferY + searchRegionIndex;

            uint8_t searchAreaBuffer[MAX_SEARCH_AREA_SIZE];

            extract8_bitdata_safe_sub(
                ptr16,
                refPicPtr->strideY,
                searchAreaBuffer,
                MAX_TATAL_SEARCH_AREA_WIDTH,
#if FIX_ME_SR_10BIT
                search_area_width + context_ptr->sb_side + ME_FILTER_TAP,
                search_area_height + context_ptr->sb_side + ME_FILTER_TAP,
#else
                MAX_TATAL_SEARCH_AREA_WIDTH,
                MAX_TATAL_SEARCH_AREA_HEIGHT,
#endif
                EB_FALSE,
                asm_type);

            context_ptr->integer_buffer_ptr[listIndex][0] = &(searchAreaBuffer[0]);
            context_ptr->interpolated_full_stride[listIndex][0] = MAX_TATAL_SEARCH_AREA_WIDTH;

        }
        else {
            context_ptr->integer_buffer_ptr[listIndex][0] = &(refPicPtr->bufferY[searchRegionIndex]);
            context_ptr->interpolated_full_stride[listIndex][0] = refPicPtr->strideY;
        }

        // Move to the top left of the search region
        xTopLeftSearchRegion = (int16_t)(refPicPtr->origin_x + sb_origin_x) + x_search_area_origin;
        yTopLeftSearchRegion = (int16_t)(refPicPtr->origin_y + sb_origin_y) + y_search_area_origin;
        searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * refPicPtr->strideY;

        //849 * 4 + 5 block are supported
        InitializeBuffer_32bits_funcPtrArray[(uint32_t)asm_type](context_ptr->p_sb_best_sad[listIndex][refPicIndex], (MAX_SS_ME_PU_COUNT / 4), 1, MAX_SAD_VALUE);

        context_ptr->p_best_sad4x4 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][0]);
        context_ptr->p_best_mv4x4 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][0]);

        context_ptr->p_best_sad8x8 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][256 * number_of_sb_quad]);
        context_ptr->p_best_mv8x8 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][256 * number_of_sb_quad]);

        context_ptr->p_best_sad16x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][320 * number_of_sb_quad]);
        context_ptr->p_best_mv16x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][320 * number_of_sb_quad]);

        context_ptr->p_best_sad32x32 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][336 * number_of_sb_quad]);
        context_ptr->p_best_mv32x32 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][336 * number_of_sb_quad]);

        context_ptr->p_best_sad64x64 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][340 * number_of_sb_quad]);
        context_ptr->p_best_mv64x64 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][340 * number_of_sb_quad]);

        context_ptr->p_best_sad8x4 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][341 * number_of_sb_quad]);
        context_ptr->p_best_mv8x4 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][341 * number_of_sb_quad]);

        context_ptr->p_best_sad4x8 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][469 * number_of_sb_quad]);
        context_ptr->p_best_mv4x8 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][469 * number_of_sb_quad]);

        context_ptr->p_best_sad4x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][597 * number_of_sb_quad]);
        context_ptr->p_best_mv4x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][597 * number_of_sb_quad]);

        context_ptr->p_best_sad16x4 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][661 * number_of_sb_quad]);
        context_ptr->p_best_mv16x4 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][661 * number_of_sb_quad]);

        context_ptr->p_best_sad16x8 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][725 * number_of_sb_quad]);
        context_ptr->p_best_mv16x8 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][725 * number_of_sb_quad]);

        context_ptr->p_best_sad8x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][757 * number_of_sb_quad]);
        context_ptr->p_best_mv8x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][757 * number_of_sb_quad]);

        context_ptr->p_best_sad32x8 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][789 * number_of_sb_quad]);
        context_ptr->p_best_mv32x8 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][789 * number_of_sb_quad]);

        context_ptr->p_best_sad8x32 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][805 * number_of_sb_quad]);
        context_ptr->p_best_mv8x32 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][805 * number_of_sb_quad]);

        context_ptr->p_best_sad32x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][821 * number_of_sb_quad]);
        context_ptr->p_best_mv32x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][821 * number_of_sb_quad]);

        context_ptr->p_best_sad16x32 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][829 * number_of_sb_quad]);
        context_ptr->p_best_mv16x32 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][829 * number_of_sb_quad]);

        context_ptr->p_best_sad64x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][837 * number_of_sb_quad]);
        context_ptr->p_best_mv64x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][837 * number_of_sb_quad]);

        context_ptr->p_best_sad16x64 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][841 * number_of_sb_quad]);
        context_ptr->p_best_mv16x64 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][841 * number_of_sb_quad]);

        context_ptr->p_best_sad64x32 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][845 * number_of_sb_quad]);
        context_ptr->p_best_mv64x32 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][845 * number_of_sb_quad]);

        context_ptr->p_best_sad32x64 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][847 * number_of_sb_quad]);
        context_ptr->p_best_mv32x64 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][847 * number_of_sb_quad]);

        if (sequence_control_set_ptr->sb_size == BLOCK_128X128) {

            context_ptr->p_best_sad128x64 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][849 * number_of_sb_quad]);
            context_ptr->p_best_mv128x64 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][849 * number_of_sb_quad]);

            context_ptr->p_best_sad64x128 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][(849 * number_of_sb_quad) + 2]);
            context_ptr->p_best_mv64x128 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][(849 * number_of_sb_quad) + 2]);

            context_ptr->p_best_sad128x128 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][(849 * number_of_sb_quad) + 4]);
            context_ptr->p_best_mv128x128 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][(849 * number_of_sb_quad) + 4]);

        }


        in_loop_me_fullpel_search_sblock(
            context_ptr,
            listIndex,
            x_search_area_origin,
            y_search_area_origin,
            search_area_width,
            search_area_height,
            number_of_sb_quad,
            asm_type);

        if (picture_control_set_ptr->parent_pcs_ptr->use_subpel_flag == 1) {


            // Move to the top left of the search region
            xTopLeftSearchRegion = (int16_t)(refPicPtr->origin_x + sb_origin_x) + x_search_area_origin;
            yTopLeftSearchRegion = (int16_t)(refPicPtr->origin_y + sb_origin_y) + y_search_area_origin;
            searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * refPicPtr->strideY;

            // Interpolate the search region for Half-Pel Refinements
            // H - AVC Style

            in_loop_me_interpolate_search_region_avc_style(
                context_ptr,
                listIndex,
                context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[listIndex][0]),
                context_ptr->interpolated_full_stride[listIndex][0],
                (uint32_t)search_area_width + (context_ptr->sb_side - 1),
                (uint32_t)search_area_height + (context_ptr->sb_side - 1),
                8,
                asm_type);

            // Half-Pel Refinement [8 search positions]
            in_loop_me_halfpel_search_sblock(
                sequence_control_set_ptr,
                context_ptr,
                &(context_ptr->pos_b_buffer[listIndex][0][(ME_FILTER_TAP >> 1) * context_ptr->interpolated_stride]),
                &(context_ptr->pos_h_buffer[listIndex][0][1]),
                &(context_ptr->pos_j_buffer[listIndex][0][0]),
                x_search_area_origin,
                y_search_area_origin,
                asm_type);

            // Quarter-Pel Refinement [8 search positions]
            in_loop_me_quarterpel_search_sblock(
                context_ptr,
                context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[listIndex][0]),
                context_ptr->interpolated_full_stride[listIndex][0],
                &(context_ptr->pos_b_buffer[listIndex][0][(ME_FILTER_TAP >> 1) * context_ptr->interpolated_stride]),  //points to b position of the figure above
                &(context_ptr->pos_h_buffer[listIndex][0][1]),                                                      //points to h position of the figure above
                &(context_ptr->pos_j_buffer[listIndex][0][0]),                                                      //points to j position of the figure above
                x_search_area_origin,
                y_search_area_origin,
                asm_type);

        }

    }

    // Nader - Bipred candidate can be generated here if needed.


    uint32_t max_number_of_block_in_sb = sequence_control_set_ptr->sb_size == BLOCK_128X128 ? MAX_SS_ME_PU_COUNT : 849;

    for (listIndex = REF_LIST_0; listIndex <= numOfListToSearch; ++listIndex) {

        uint32_t block_index;
        uint32_t block_offset;
        uint32_t nidx;
        uint32_t candidate_cnt = 0;

        for (block_index = 0; block_index < max_number_of_block_in_sb; ++block_index) {

            //4x4
            if (block_index < start_idx_8x8) {
                block_offset = (block_index / 256) * 256;
                nidx = tab4x4[block_index - block_offset] + block_offset;
            } //8x8
            else if (block_index < start_idx_16x16) {
                block_offset = ((block_index - start_idx_8x8) / 64) * 64;
                nidx = tab8x8[block_index - start_idx_8x8 - block_offset] + block_offset + start_idx_8x8;
            }//16x16
            else if (block_index < start_idx_32x32) {
                block_offset = ((block_index - start_idx_16x16) / 16) * 16;
                nidx = tab16x16[block_index - start_idx_16x16 - block_offset] + block_offset + start_idx_16x16;
            }//32x32
            else if (block_index < start_idx_64x64) {
                block_offset = ((block_index - start_idx_32x32) / 4) * 4;
                nidx = tab32x32[block_index - start_idx_32x32 - block_offset] + block_offset + start_idx_32x32;
            } //64x64
            else if (block_index < start_idx_8x4) {
                block_offset = (block_index - start_idx_64x64);
                nidx = block_offset + start_idx_64x64;
            } //8x4
            else if (block_index < start_idx_4x8) {
                block_offset = ((block_index - start_idx_8x4) / 128) * 128;
                nidx = tab8x4[block_index - start_idx_8x4 - block_offset] + block_offset + start_idx_8x4;
            }//4x8
            else if (block_index < start_idx_4x16) {
                block_offset = ((block_index - start_idx_4x8) / 128) * 128;
                nidx = tab4x8[block_index - start_idx_4x8 - block_offset] + block_offset + start_idx_4x8;
            }//4x16
            else if (block_index < start_idx_16x4) {
                block_offset = ((block_index - start_idx_4x16) / 64) * 64;
                nidx = tab4x16[block_index - start_idx_4x16 - block_offset] + block_offset + start_idx_4x16;
            }//16x4
            else if (block_index < start_idx_16x8) {
                block_offset = ((block_index - start_idx_16x4) / 64) * 64;
                nidx = tab16x4[block_index - start_idx_16x4 - block_offset] + block_offset + start_idx_16x4;
            }//16x8
            else if (block_index < start_idx_8x16) {
                block_offset = ((block_index - start_idx_16x8) / 32) * 32;
                nidx = tab16x8[block_index - start_idx_16x8 - block_offset] + block_offset + start_idx_16x8;
            }//8x16
            else if (block_index < start_idx_32x8) {
                block_offset = ((block_index - start_idx_8x16) / 32) * 32;
                nidx = tab8x16[block_index - start_idx_8x16 - block_offset] + block_offset + start_idx_8x16;

            }//32x8
            else if (block_index < start_idx_8x32) {
                block_offset = ((block_index - start_idx_32x8) / 16) * 16;
                nidx = tab32x8[block_index - start_idx_32x8 - block_offset] + block_offset + start_idx_32x8;
            }//8x32
            else if (block_index < start_idx_32x16) {
                block_offset = ((block_index - start_idx_8x32) / 16) * 16;
                nidx = tab8x32[block_index - start_idx_8x32 - block_offset] + block_offset + start_idx_8x32;

            }//32x16
            else if (block_index < start_idx_16x32) {
                block_offset = ((block_index - start_idx_32x16) / 8) * 8;
                nidx = tab32x16[block_index - start_idx_32x16 - block_offset] + block_offset + start_idx_32x16;
            }//16x32
            else if (block_index < start_idx_64x16) {
                block_offset = ((block_index - start_idx_16x32) / 8) * 8;
                nidx = tab16x32[block_index - start_idx_16x32 - block_offset] + block_offset + start_idx_16x32;
            }//64x16
            else if (block_index < start_idx_16x64) {
                block_offset = ((block_index - start_idx_64x16) / 4) * 4;
                nidx = tab64x16[block_index - start_idx_64x16 - block_offset] + block_offset + start_idx_64x16;

            }//16x64
            else if (block_index < start_idx_64x32) {
                block_offset = ((block_index - start_idx_16x64) / 4) * 4;
                nidx = tab16x64[block_index - start_idx_16x64 - block_offset] + block_offset + start_idx_16x64;
            }//64x32
            else if (block_index < start_idx_32x64) {
                block_offset = ((block_index - start_idx_64x32) / 2) * 2;
                nidx = tab64x32[block_index - start_idx_64x32 - block_offset] + block_offset + start_idx_64x32;

            }//32x64
            else if (block_index < start_idx_128x64) {
                block_offset = ((block_index - start_idx_32x64) / 2) * 2;
                nidx = tab32x64[block_index - start_idx_32x64 - block_offset] + block_offset + start_idx_32x64;

            }//128x64, //64x128 and 128x128
            else {
                nidx = block_index;
            }

            context_ptr->inloop_me_mv[0][0][candidate_cnt][0] = _MVXT(context_ptr->p_sb_best_mv[0][0][nidx]);
            context_ptr->inloop_me_mv[0][0][candidate_cnt][1] = _MVYT(context_ptr->p_sb_best_mv[0][0][nidx]);
            context_ptr->inloop_me_mv[1][0][candidate_cnt][0] = _MVXT(context_ptr->p_sb_best_mv[1][0][nidx]);
            context_ptr->inloop_me_mv[1][0][candidate_cnt][1] = _MVYT(context_ptr->p_sb_best_mv[1][0][nidx]);
            candidate_cnt++;
        }

    }

    return return_error;
}








