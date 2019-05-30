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
#include <string.h>

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbTransformUnit.h"
#include "EbRateDistortionCost.h"
#include "EbDeblockingFilter.h"
#include "EbPictureOperators.h"

#include "EbModeDecisionProcess.h"
#include "EbEncDecProcess.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbTransforms.h"
#include "EbModeDecisionConfigurationProcess.h"
#include "EbIntraPrediction.h"
#include "aom_dsp_rtcd.h"
#include "EbCodingLoop.h"
void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);
extern void av1_predict_intra_block(
    TileInfo                    *tile,
    STAGE                       stage,
    const BlockGeom            *blk_geom,
    const Av1Common *cm,
    int32_t wpx,
    int32_t hpx,
    TxSize tx_size,
    PredictionMode mode,
    int32_t angle_delta,
    int32_t use_palette,
    FilterIntraMode filter_intra_mode,
    uint8_t* topNeighArray,
    uint8_t* leftNeighArray,
    EbPictureBufferDesc  *recon_buffer,
    int32_t col_off,
    int32_t row_off,
    int32_t plane,
    BlockSize bsize,
#if ATB_EP
    uint32_t tu_org_x_pict,
    uint32_t tu_org_y_pict,
#endif
    uint32_t bl_org_x_pict,
    uint32_t bl_org_y_pict,
    uint32_t bl_org_x_mb,
    uint32_t bl_org_y_mb);

void av1_predict_intra_block_16bit(
    TileInfo               *tile,
    EncDecContext         *context_ptr,
    const Av1Common *cm,
    int32_t wpx,
    int32_t hpx,
    TxSize tx_size,
    PredictionMode mode,
    int32_t angle_delta,
    int32_t use_palette,
    FilterIntraMode filter_intra_mode,
    uint16_t* topNeighArray,
    uint16_t* leftNeighArray,
    EbPictureBufferDesc  *recon_buffer,
    int32_t col_off,
    int32_t row_off,
    int32_t plane,
    BlockSize bsize,
    uint32_t bl_org_x_pict,
    uint32_t bl_org_y_pict);

/*******************************************
* set Penalize Skip Flag
*
* Summary: Set the penalize_skipflag to true
* When there is luminance/chrominance change
* or in noisy clip with low motion at meduim
* varince area
*
*******************************************/

#define S32 32*32
#define S16 16*16
#define S8  8*8
#define S4  4*4

typedef void(*EB_AV1_ENCODE_LOOP_FUNC_PTR)(
    PictureControlSet    *picture_control_set_ptr,
    EncDecContext       *context_ptr,
    LargestCodingUnit   *sb_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    uint32_t                 cb_qp,
    EbPictureBufferDesc *predSamples,             // no basis/offset
    EbPictureBufferDesc *coeffSamplesTB,          // lcu based
    EbPictureBufferDesc *residual16bit,           // no basis/offset
    EbPictureBufferDesc *transform16bit,          // no basis/offset
    EbPictureBufferDesc *inverse_quant_buffer,
    int16_t                *transformScratchBuffer,
    EbAsm                 asm_type,
    uint32_t                  *count_non_zero_coeffs,
    uint32_t                 component_mask,
    uint32_t                   use_delta_qp,
    uint32_t                 dZoffset,
    uint16_t                 *eob,
    MacroblockPlane       *candidate_plane);

typedef void(*EB_AV1_GENERATE_RECON_FUNC_PTR)(
    EncDecContext       *context_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    EbPictureBufferDesc *predSamples,     // no basis/offset
    EbPictureBufferDesc *residual16bit,    // no basis/offset
    int16_t                *transformScratchBuffer,
    uint32_t                 component_mask,
    uint16_t                *eob,
    EbAsm                 asm_type);

/***************************************************
* Update Intra Mode Neighbor Arrays
***************************************************/
static void EncodePassUpdateIntraModeNeighborArrays(
#if DC_SIGN_CONTEXT_EP
    EncDecContext     *context_ptr,
    NeighborArrayUnit *luma_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit *cb_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit *cr_dc_sign_level_coeff_neighbor_array,
#endif
    NeighborArrayUnit *mode_type_neighbor_array,
    NeighborArrayUnit *intra_luma_mode_neighbor_array,
    NeighborArrayUnit *intra_chroma_mode_neighbor_array,
    uint8_t            luma_mode,
    uint8_t            chroma_mode,
    uint32_t           origin_x,
    uint32_t           origin_y,
    uint32_t           width,
    uint32_t           height,
    uint32_t           width_uv,
    uint32_t           height_uv,
    uint32_t           component_mask)
{
    uint8_t modeType = INTRA_MODE;

    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {
        // Mode Type Update
        neighbor_array_unit_mode_write(
            mode_type_neighbor_array,
            &modeType,
            origin_x,
            origin_y,
            width,
            height,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);

        // Intra Luma Mode Update
        neighbor_array_unit_mode_write(
            intra_luma_mode_neighbor_array,
            &luma_mode,
            origin_x,
            origin_y,
            width,
            height,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }
    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        // Intra Luma Mode Update
        neighbor_array_unit_mode_write(
            intra_chroma_mode_neighbor_array,
            &chroma_mode,
            ((origin_x >> 3) << 3) / 2,
            ((origin_y >> 3) << 3) / 2,
            width_uv,
            height_uv,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

#if DC_SIGN_CONTEXT_EP
    uint8_t dc_sign_level_coeff;
    uint16_t txb_count = context_ptr->blk_geom->txb_count[context_ptr->cu_ptr->tx_depth];

    for (uint8_t txb_itr = 0; txb_itr < txb_count; txb_itr++) {
        dc_sign_level_coeff = (int32_t)context_ptr->cu_ptr->quantized_dc[0][txb_itr];

        neighbor_array_unit_mode_write(
            luma_dc_sign_level_coeff_neighbor_array,
            (uint8_t*)&dc_sign_level_coeff,
            context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[context_ptr->cu_ptr->tx_depth][txb_itr],
            context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[context_ptr->cu_ptr->tx_depth][txb_itr],
            context_ptr->blk_geom->tx_width[context_ptr->cu_ptr->tx_depth][txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->cu_ptr->tx_depth][txb_itr],
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
        dc_sign_level_coeff = (int32_t)context_ptr->cu_ptr->quantized_dc[1][0];

        neighbor_array_unit_mode_write(
            cb_dc_sign_level_coeff_neighbor_array,
            (uint8_t*)&dc_sign_level_coeff,
            ((origin_x >> 3) << 3) / 2,
            ((origin_y >> 3) << 3) / 2,
            width_uv,
            height_uv,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        dc_sign_level_coeff = (int32_t)context_ptr->cu_ptr->quantized_dc[2][0];

        neighbor_array_unit_mode_write(
            cr_dc_sign_level_coeff_neighbor_array,
            (uint8_t*)&dc_sign_level_coeff,
            ((origin_x >> 3) << 3) / 2,
            ((origin_y >> 3) << 3) / 2,
            width_uv,
            height_uv,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }
#endif
    return;
}

/***************************************************
* Update Inter Mode Neighbor Arrays
***************************************************/
static void EncodePassUpdateInterModeNeighborArrays(
#if DC_SIGN_CONTEXT_EP
    EncDecContext     *context_ptr,
    NeighborArrayUnit *luma_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit *cb_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit *cr_dc_sign_level_coeff_neighbor_array,
#endif
    NeighborArrayUnit *mode_type_neighbor_array,
    NeighborArrayUnit *mv_neighbor_array,
    NeighborArrayUnit *skipNeighborArray,
    MvUnit            *mv_unit,
    uint8_t           *skip_flag,
    uint32_t           origin_x,
    uint32_t           origin_y,
    uint32_t           bwidth,
    uint32_t           bheight)
{
    uint8_t modeType = INTER_MODE;

    // Mode Type Update
    neighbor_array_unit_mode_write(
        mode_type_neighbor_array,
        &modeType,
        origin_x,
        origin_y,
        bwidth,
        bheight,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    // Motion Vector Unit
    neighbor_array_unit_mode_write(
        mv_neighbor_array,
        (uint8_t*)mv_unit,
        origin_x,
        origin_y,
        bwidth,
        bheight,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    // Skip Flag
    neighbor_array_unit_mode_write(
        skipNeighborArray,
        skip_flag,
        origin_x,
        origin_y,
        bwidth,
        bheight,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

#if DC_SIGN_CONTEXT_EP
    uint8_t dc_sign_level_coeff;
    uint16_t txb_count = context_ptr->blk_geom->txb_count[context_ptr->cu_ptr->tx_depth];

    for (uint8_t txb_itr = 0; txb_itr < txb_count; txb_itr++) {
        dc_sign_level_coeff = (int32_t)context_ptr->cu_ptr->quantized_dc[0][txb_itr];

        neighbor_array_unit_mode_write(
            luma_dc_sign_level_coeff_neighbor_array,
            (uint8_t*)&dc_sign_level_coeff,
            context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[context_ptr->cu_ptr->tx_depth][txb_itr],
            context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[context_ptr->cu_ptr->tx_depth][txb_itr],
            context_ptr->blk_geom->tx_width[context_ptr->cu_ptr->tx_depth][txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->cu_ptr->tx_depth][txb_itr],
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    }

    dc_sign_level_coeff = (int32_t)context_ptr->cu_ptr->quantized_dc[1][0];

    neighbor_array_unit_mode_write(
        cb_dc_sign_level_coeff_neighbor_array,
        (uint8_t*)&dc_sign_level_coeff,
        ((origin_x >> 3) << 3) / 2,
        ((origin_y >> 3) << 3) / 2,
        context_ptr->blk_geom->bwidth_uv,
        context_ptr->blk_geom->bheight_uv,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    dc_sign_level_coeff = (int32_t)context_ptr->cu_ptr->quantized_dc[2][0];

    neighbor_array_unit_mode_write(
        cr_dc_sign_level_coeff_neighbor_array,
        (uint8_t*)&dc_sign_level_coeff,
        ((origin_x >> 3) << 3) / 2,
        ((origin_y >> 3) << 3) / 2,
        context_ptr->blk_geom->bwidth_uv,
        context_ptr->blk_geom->bheight_uv,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

#endif
    return;
}

/***************************************************
* Update Recon Samples Neighbor Arrays
***************************************************/
static void EncodePassUpdateReconSampleNeighborArrays(
    NeighborArrayUnit     *lumaReconSampleNeighborArray,
    NeighborArrayUnit     *cbReconSampleNeighborArray,
    NeighborArrayUnit     *crReconSampleNeighborArray,
    EbPictureBufferDesc   *recon_buffer,
    uint32_t                   origin_x,
    uint32_t                   origin_y,
    uint32_t                   width,
    uint32_t                   height,
    uint32_t                   bwidth_uv,
    uint32_t                   bheight_uv,
    uint32_t                   component_mask,
    EbBool                  is16bit)
{
    uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
    uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4

    if (is16bit == EB_TRUE) {
        if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK)
        {
            // Recon Samples - Luma
            neighbor_array_unit16bit_sample_write(
                lumaReconSampleNeighborArray,
                (uint16_t*)(recon_buffer->buffer_y),
                recon_buffer->stride_y,
                recon_buffer->origin_x + origin_x,
                recon_buffer->origin_y + origin_y,
                origin_x,
                origin_y,
                width,
                height,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        }

        if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
            // Recon Samples - Cb
            neighbor_array_unit16bit_sample_write(
                cbReconSampleNeighborArray,
                (uint16_t*)(recon_buffer->buffer_cb),
                recon_buffer->stride_cb,
                (recon_buffer->origin_x + round_origin_x) >> 1,
                (recon_buffer->origin_y + round_origin_y) >> 1,
                round_origin_x >> 1,
                round_origin_y >> 1,
                bwidth_uv,
                bheight_uv,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK);

            // Recon Samples - Cr
            neighbor_array_unit16bit_sample_write(
                crReconSampleNeighborArray,
                (uint16_t*)(recon_buffer->buffer_cr),
                recon_buffer->stride_cr,
                (recon_buffer->origin_x + round_origin_x) >> 1,
                (recon_buffer->origin_y + round_origin_y) >> 1,
                round_origin_x >> 1,
                round_origin_y >> 1,
                bwidth_uv,
                bheight_uv,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        }
    }
    else {
        if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK)
        {
            // Recon Samples - Luma
            neighbor_array_unit_sample_write(
                lumaReconSampleNeighborArray,
                recon_buffer->buffer_y,
                recon_buffer->stride_y,
                recon_buffer->origin_x + origin_x,
                recon_buffer->origin_y + origin_y,
                origin_x,
                origin_y,
                width,
                height,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        }

        if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
            // Recon Samples - Cb
            neighbor_array_unit_sample_write(
                cbReconSampleNeighborArray,
                recon_buffer->buffer_cb,
                recon_buffer->stride_cb,
                (recon_buffer->origin_x + round_origin_x) >> 1,
                (recon_buffer->origin_y + round_origin_y) >> 1,
                round_origin_x >> 1,
                round_origin_y >> 1,
                bwidth_uv,
                bheight_uv,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK);

            // Recon Samples - Cr
            neighbor_array_unit_sample_write(
                crReconSampleNeighborArray,
                recon_buffer->buffer_cr,
                recon_buffer->stride_cr,
                (recon_buffer->origin_x + round_origin_x) >> 1,
                (recon_buffer->origin_y + round_origin_y) >> 1,
                round_origin_x >> 1,
                round_origin_y >> 1,
                bwidth_uv,
                bheight_uv,
                NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        }
    }

    return;
}

/************************************************************
* Update Intra Luma Neighbor Modes
************************************************************/
void GeneratePuIntraLumaNeighborModes(
    CodingUnit            *cu_ptr,
    uint32_t                   pu_origin_x,
    uint32_t                   pu_origin_y,
    uint32_t                   sb_sz,
    NeighborArrayUnit     *intraLumaNeighborArray,
    NeighborArrayUnit     *intraChromaNeighborArray,
    NeighborArrayUnit     *mode_type_neighbor_array)
{
    (void)sb_sz;

    uint32_t modeTypeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        mode_type_neighbor_array,
        pu_origin_y);
    uint32_t modeTypeTopNeighborIndex = get_neighbor_array_unit_top_index(
        mode_type_neighbor_array,
        pu_origin_x);
    uint32_t intraLumaModeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        intraLumaNeighborArray,
        pu_origin_y);
    uint32_t intraLumaModeTopNeighborIndex = get_neighbor_array_unit_top_index(
        intraLumaNeighborArray,
        pu_origin_x);

    uint32_t puOriginX_round = (pu_origin_x >> 3) << 3;
    uint32_t puOriginY_round = (pu_origin_y >> 3) << 3;

    uint32_t intraChromaModeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        intraChromaNeighborArray,
        puOriginY_round >> 1);
    uint32_t intraChromaModeTopNeighborIndex = get_neighbor_array_unit_top_index(
        intraChromaNeighborArray,
        puOriginX_round >> 1);

    (&cu_ptr->prediction_unit_array[0])->intra_luma_left_mode = (uint32_t)(
        (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] != INTRA_MODE) ? DC_PRED/*EB_INTRA_DC*/ :
        (uint32_t)intraLumaNeighborArray->left_array[intraLumaModeLeftNeighborIndex]);

    (&cu_ptr->prediction_unit_array[0])->intra_luma_top_mode = (uint32_t)(
        (mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] != INTRA_MODE) ? DC_PRED/*EB_INTRA_DC*/ :
        (uint32_t)intraLumaNeighborArray->top_array[intraLumaModeTopNeighborIndex]);       //   use DC. This seems like we could use a LCU-width

    uint32_t modeTypeLeftNeighborIndex_round = get_neighbor_array_unit_left_index(
        mode_type_neighbor_array,
        puOriginY_round);
    uint32_t modeTypeTopNeighborIndex_round = get_neighbor_array_unit_top_index(
        mode_type_neighbor_array,
        puOriginX_round);

    (&cu_ptr->prediction_unit_array[0])->intra_chroma_left_mode = (uint32_t)(
        (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex_round] != INTRA_MODE) ? UV_DC_PRED :
        (uint32_t)intraChromaNeighborArray->left_array[intraChromaModeLeftNeighborIndex]);

    (&cu_ptr->prediction_unit_array[0])->intra_chroma_top_mode = (uint32_t)(
        (mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex_round] != INTRA_MODE) ? UV_DC_PRED :
        (uint32_t)intraChromaNeighborArray->top_array[intraChromaModeTopNeighborIndex]);       //   use DC. This seems like we could use a LCU-width

    return;
}

void encode_pass_tx_search(
    PictureControlSet            *picture_control_set_ptr,
    EncDecContext                *context_ptr,
    LargestCodingUnit            *sb_ptr,
    uint32_t                       cb_qp,
    EbPictureBufferDesc          *coeffSamplesTB,
    EbPictureBufferDesc          *residual16bit,
    EbPictureBufferDesc          *transform16bit,
    EbPictureBufferDesc          *inverse_quant_buffer,
    int16_t                        *transformScratchBuffer,
    EbAsm                          asm_type,
    uint32_t                       *count_non_zero_coeffs,
    uint32_t                       component_mask,
    uint32_t                       use_delta_qp,
    uint32_t                       dZoffset,
    uint16_t                       *eob,
    MacroblockPlane                *candidate_plane);

/**********************************************************
* Encode Loop
*
* Summary: Performs a H.265 conformant
*   Transform, Quantization  and Inverse Quantization of a TU.
*
* Inputs:
*   origin_x
*   origin_y
*   txb_size
*   sb_sz
*   input - input samples (position sensitive)
*   pred - prediction samples (position independent)
*
* Outputs:
*   Inverse quantized coeff - quantization indices (position sensitive)
*
**********************************************************/
static void Av1EncodeLoop(
    PictureControlSet    *picture_control_set_ptr,
    EncDecContext       *context_ptr,
    LargestCodingUnit   *sb_ptr,
    uint32_t                 origin_x,   //pic based tx org x
    uint32_t                 origin_y,   //pic based tx org y
    uint32_t                 cb_qp,
    EbPictureBufferDesc *predSamples,             // no basis/offset
    EbPictureBufferDesc *coeffSamplesTB,          // lcu based
    EbPictureBufferDesc *residual16bit,           // no basis/offset
    EbPictureBufferDesc *transform16bit,          // no basis/offset
    EbPictureBufferDesc *inverse_quant_buffer,
    int16_t                *transformScratchBuffer,
    EbAsm                 asm_type,
    uint32_t                  *count_non_zero_coeffs,
    uint32_t                 component_mask,
    uint32_t                   use_delta_qp,
    uint32_t                 dZoffset,
    uint16_t                 *eob,
    MacroblockPlane       *candidate_plane){
    (void)dZoffset;
    (void)use_delta_qp;
    (void)cb_qp;

    //    uint32_t                 chroma_qp = cb_qp;
    CodingUnit          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    //    EB_SLICE               slice_type = sb_ptr->picture_control_set_ptr->slice_type;
    //    uint32_t                 temporal_layer_index = sb_ptr->picture_control_set_ptr->temporal_layer_index;
    uint32_t                 qp = cu_ptr->qp;
    EbPictureBufferDesc  *input_samples = context_ptr->input_samples;

    uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
    uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4

    const uint32_t input_luma_offset = ((origin_y + input_samples->origin_y)          * input_samples->stride_y) + (origin_x + input_samples->origin_x);
    const uint32_t input_cb_offset = (((round_origin_y + input_samples->origin_y) >> 1)    * input_samples->stride_cb) + ((round_origin_x + input_samples->origin_x) >> 1);
    const uint32_t input_cr_offset = (((round_origin_y + input_samples->origin_y) >> 1)    * input_samples->stride_cr) + ((round_origin_x + input_samples->origin_x) >> 1);
    const uint32_t pred_luma_offset = ((predSamples->origin_y + origin_y)        * predSamples->stride_y) + (predSamples->origin_x + origin_x);
    const uint32_t pred_cb_offset = (((predSamples->origin_y + round_origin_y) >> 1)  * predSamples->stride_cb) + ((predSamples->origin_x + round_origin_x) >> 1);
    const uint32_t pred_cr_offset = (((predSamples->origin_y + round_origin_y) >> 1)  * predSamples->stride_cr) + ((predSamples->origin_x + round_origin_x) >> 1);

#if ATB_SUPPORT
    const uint32_t scratch_luma_offset = context_ptr->blk_geom->tx_org_x[cu_ptr->tx_depth][context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[cu_ptr->tx_depth][context_ptr->txb_itr] * SB_STRIDE_Y;
    const uint32_t scratch_cb_offset = ROUND_UV(context_ptr->blk_geom->tx_org_x[cu_ptr->tx_depth][context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[cu_ptr->tx_depth][context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t scratch_cr_offset = ROUND_UV(context_ptr->blk_geom->tx_org_x[cu_ptr->tx_depth][context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[cu_ptr->tx_depth][context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
#else
    const uint32_t scratch_luma_offset = context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr] * SB_STRIDE_Y;
    const uint32_t scratch_cb_offset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t scratch_cr_offset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
#endif

    const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

    const uint32_t coeff1dOffsetChroma = context_ptr->coded_area_sb_uv;
    UNUSED(coeff1dOffsetChroma);

#if !OPT_LOSSLESS_0
    //uint8_t enable_contouring_qc_update_flag;
    //enable_contouring_qc_update_flag = DeriveContouringClass(
    //    sb_ptr->picture_control_set_ptr->parent_pcs_ptr,
    //    sb_ptr->index,
    //    cu_ptr->leaf_index) && (cu_ptr->qp < sb_ptr->picture_control_set_ptr->picture_qp);
#endif

    context_ptr->three_quad_energy = 0;
    //**********************************
    // Luma
    //**********************************
    if (component_mask == PICTURE_BUFFER_DESC_FULL_MASK || component_mask == PICTURE_BUFFER_DESC_LUMA_MASK)
    {
        ResidualKernel(
            input_samples->buffer_y + input_luma_offset,
            input_samples->stride_y,
            predSamples->buffer_y + pred_luma_offset,
            predSamples->stride_y,
            ((int16_t*)residual16bit->buffer_y) + scratch_luma_offset,
            residual16bit->stride_y,
#if ATB_SUPPORT
            context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
            context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);
#endif
#if ATB_MD
        uint8_t  tx_search_skip_fag = (picture_control_set_ptr->parent_pcs_ptr->tx_search_level == TX_SEARCH_ENC_DEC && (picture_control_set_ptr->parent_pcs_ptr->atb_mode == 0 || cu_ptr->prediction_mode_flag == INTER_MODE)) ? get_skip_tx_search_flag(
#else
        uint8_t  tx_search_skip_fag = picture_control_set_ptr->parent_pcs_ptr->tx_search_level == TX_SEARCH_ENC_DEC ? get_skip_tx_search_flag(
#endif
#if BYPASS_USELESS_TX_SEARCH
            context_ptr->blk_geom,
#else
            context_ptr->blk_geom->sq_size,
#endif
            MAX_MODE_COST,
            0,
            1) : 1;

        if (!tx_search_skip_fag) {
                encode_pass_tx_search(
                    picture_control_set_ptr,
                    context_ptr,
                    sb_ptr,
                    cb_qp,
                    coeffSamplesTB,
                    residual16bit,
                    transform16bit,
                    inverse_quant_buffer,
                    transformScratchBuffer,
                    asm_type,
                    count_non_zero_coeffs,
                    component_mask,
                    use_delta_qp,
                    dZoffset,
                    eob,
                    candidate_plane);
        }

        av1_estimate_transform(
            ((int16_t*)residual16bit->buffer_y) + scratch_luma_offset,
            residual16bit->stride_y,
            ((TranLow*)transform16bit->buffer_y) + coeff1dOffset,
            NOT_USED_VALUE,
#if ATB_SUPPORT
            context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
#endif
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_Y],
            asm_type,
            PLANE_TYPE_Y,
#if PF_N2_SUPPORT
            DEFAULT_SHAPE);
#else
            context_ptr->trans_coeff_shape_luma);
#endif

#if DC_SIGN_CONTEXT_EP
        cu_ptr->quantized_dc[0][context_ptr->txb_itr] = av1_quantize_inv_quantize(
#else
        av1_quantize_inv_quantize(
#endif
            sb_ptr->picture_control_set_ptr,
            context_ptr->md_context,
            ((TranLow*)transform16bit->buffer_y) + coeff1dOffset,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->buffer_y) + coeff1dOffset,
            ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
            qp,
#if ATB_SUPPORT
            context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
            context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
#endif
            &eob[0],
            asm_type,
            &(count_non_zero_coeffs[0]),
#if !PF_N2_SUPPORT
            0,
#endif
            COMPONENT_LUMA,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_Y],
            &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]),
            cu_ptr->luma_txb_skip_context,
#if ATB_DC_CONTEXT_SUPPORT_0
            cu_ptr->luma_dc_sign_context[context_ptr->txb_itr],
#else
            cu_ptr->luma_dc_sign_context,
#endif
            cu_ptr->pred_mode,
            EB_TRUE);

#if BLK_SKIP_DECISION
        if (context_ptr->md_skip_blk) {
            count_non_zero_coeffs[0] = 0;
            eob[0] = 0;
        }
#endif
        txb_ptr->y_has_coeff = count_non_zero_coeffs[0] ? EB_TRUE : EB_FALSE;

        if (count_non_zero_coeffs[0] == 0) {
            // INTER. Chroma follows Luma in transform type
            if (cu_ptr->prediction_mode_flag == INTER_MODE) {
                txb_ptr->transform_type[PLANE_TYPE_Y] = DCT_DCT;
                txb_ptr->transform_type[PLANE_TYPE_UV] = DCT_DCT;
            }
            else { // INTRA
                txb_ptr->transform_type[PLANE_TYPE_Y] = DCT_DCT;
            }
        }
#if !ATB_EP || ATB_EC_NO_CFL
        if (cu_ptr->prediction_mode_flag == INTRA_MODE && (context_ptr->evaluate_cfl_ep || cu_ptr->prediction_unit_array->intra_chroma_mode == UV_CFL_PRED)) {
            EbPictureBufferDesc *reconSamples = predSamples;
            uint32_t reconLumaOffset = (reconSamples->origin_y + origin_y)            * reconSamples->stride_y + (reconSamples->origin_x + origin_x);

            if (txb_ptr->y_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
                uint8_t     *predBuffer = predSamples->buffer_y + pred_luma_offset;

                av1_inv_transform_recon8bit(
                    ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
                    predBuffer,
                    predSamples->stride_y,
#if ATB_SUPPORT
                    context_ptr->blk_geom->txsize[tx_depth][context_ptr->txb_itr],
#else
                    context_ptr->blk_geom->txsize[context_ptr->txb_itr],
#endif
                    txb_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    eob[0]);
            }
#if CFL_FIX
            if (context_ptr->blk_geom->has_uv) {
                reconLumaOffset = (reconSamples->origin_y + round_origin_y)            * reconSamples->stride_y + (reconSamples->origin_x + round_origin_x);
#endif
                // Down sample Luma
                cfl_luma_subsampling_420_lbd_c(
                    reconSamples->buffer_y + reconLumaOffset,
                    reconSamples->stride_y,
                    context_ptr->md_context->pred_buf_q3,
#if CFL_FIX
                    context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bwidth ? (context_ptr->blk_geom->bwidth_uv << 1) : context_ptr->blk_geom->bwidth,
                    context_ptr->blk_geom->bheight_uv == context_ptr->blk_geom->bheight ? (context_ptr->blk_geom->bheight_uv << 1) : context_ptr->blk_geom->bheight);
#else
                    context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);
#endif
#if ATB_SUPPORT
                int32_t round_offset = ((context_ptr->blk_geom->tx_width_uv[tx_depth][context_ptr->txb_itr])*(context_ptr->blk_geom->tx_height_uv[tx_depth][context_ptr->txb_itr])) / 2;
#else
                int32_t round_offset = ((context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr])*(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr])) / 2;
#endif

                subtract_average(
                    context_ptr->md_context->pred_buf_q3,
#if ATB_SUPPORT
                    context_ptr->blk_geom->tx_width_uv[tx_depth][context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[tx_depth][context_ptr->txb_itr],
                    round_offset,
                    LOG2F(context_ptr->blk_geom->tx_width_uv[tx_depth][context_ptr->txb_itr]) + LOG2F(context_ptr->blk_geom->tx_height_uv[tx_depth][context_ptr->txb_itr]));

#else
                    context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                    round_offset,
                    LOG2F(context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr]) + LOG2F(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]));
#endif

                if (context_ptr->evaluate_cfl_ep)
                {
                    // 3: Loop over alphas and find the best or choose DC
                    // Use the 1st spot of the candidate buffer to hold cfl settings: (1) to use same kernel as MD for CFL evaluation: cfl_rd_pick_alpha() (toward unification), (2) to avoid dedicated buffers for CFL evaluation @ EP (toward less memory)
                    ModeDecisionCandidateBuffer  *candidateBuffer = &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]);

                    // Input(s)
                    candidateBuffer->candidate_ptr->type = INTRA_MODE;
                    candidateBuffer->candidate_ptr->intra_luma_mode = cu_ptr->pred_mode;
                    candidateBuffer->candidate_ptr->cfl_alpha_signs = 0;
                    candidateBuffer->candidate_ptr->cfl_alpha_idx = 0;
                    context_ptr->md_context->blk_geom = context_ptr->blk_geom;

                    EbByte src_pred_ptr;
                    EbByte dst_pred_ptr;

                    // Copy Cb pred samples from ep buffer to md buffer
                    src_pred_ptr = predSamples->buffer_cb + pred_cb_offset;
                    dst_pred_ptr = &(candidateBuffer->prediction_ptr->buffer_cb[scratch_cb_offset]);
                    for (int i = 0; i < context_ptr->blk_geom->bheight_uv; i++) {
                        memcpy(dst_pred_ptr, src_pred_ptr, context_ptr->blk_geom->bwidth_uv);
                        src_pred_ptr += predSamples->stride_cb;
                        dst_pred_ptr += candidateBuffer->prediction_ptr->stride_cb;
                    }

                    // Copy Cr pred samples from ep buffer to md buffer
                    src_pred_ptr = predSamples->buffer_cr + pred_cr_offset;
                    dst_pred_ptr = &(candidateBuffer->prediction_ptr->buffer_cr[scratch_cr_offset]);
                    for (int i = 0; i < context_ptr->blk_geom->bheight_uv; i++) {
                        memcpy(dst_pred_ptr, src_pred_ptr, context_ptr->blk_geom->bwidth_uv);
                        src_pred_ptr += predSamples->stride_cr;
                        dst_pred_ptr += candidateBuffer->prediction_ptr->stride_cr;
                    }

                    cfl_rd_pick_alpha(
                        picture_control_set_ptr,
                        candidateBuffer,
                        sb_ptr,
                        context_ptr->md_context,
                        input_samples,
                        input_cb_offset,
                        scratch_cb_offset,
                        asm_type);

                    // Output(s)
                    if (candidateBuffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {
                        cu_ptr->prediction_unit_array->intra_chroma_mode = UV_CFL_PRED;
                        cu_ptr->prediction_unit_array->cfl_alpha_idx = candidateBuffer->candidate_ptr->cfl_alpha_idx;
                        cu_ptr->prediction_unit_array->cfl_alpha_signs = candidateBuffer->candidate_ptr->cfl_alpha_signs;
                        cu_ptr->prediction_unit_array->is_directional_chroma_mode_flag = EB_FALSE;
                    }
                }

                if (cu_ptr->prediction_unit_array->intra_chroma_mode == UV_CFL_PRED) {
                    int32_t alpha_q3 =
                        cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_U); // once for U, once for V

                    //TOCHANGE
                    //assert(chroma_size * CFL_BUF_LINE + chroma_size <= CFL_BUF_SQUARE);

                    cfl_predict_lbd(
                        context_ptr->md_context->pred_buf_q3,
                        predSamples->buffer_cb + pred_cb_offset,
                        predSamples->stride_cb,
                        predSamples->buffer_cb + pred_cb_offset,
                        predSamples->stride_cb,
                        alpha_q3,
                        8,
#if ATB_SUPPORT
                        context_ptr->blk_geom->tx_width_uv[tx_depth][context_ptr->txb_itr],
                        context_ptr->blk_geom->tx_height_uv[tx_depth][context_ptr->txb_itr]);
#else
                        context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                        context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif
                    alpha_q3 =
                        cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_V); // once for U, once for V

                    //TOCHANGE
                    //assert(chroma_size * CFL_BUF_LINE + chroma_size <= CFL_BUF_SQUARE);

                    cfl_predict_lbd(
                        context_ptr->md_context->pred_buf_q3,
                        predSamples->buffer_cr + pred_cr_offset,
                        predSamples->stride_cr,
                        predSamples->buffer_cr + pred_cr_offset,
                        predSamples->stride_cr,
                        alpha_q3,
                        8,
#if ATB_SUPPORT
                        context_ptr->blk_geom->tx_width_uv[tx_depth][context_ptr->txb_itr],
                        context_ptr->blk_geom->tx_height_uv[tx_depth][context_ptr->txb_itr]);
#else
                        context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                        context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif
                }
#if CFL_FIX
            }
#endif
        }
#endif
#if ATB_EP
        txb_ptr->nz_coef_count[0] = (uint16_t)count_non_zero_coeffs[0];
#endif
    }

    if (component_mask == PICTURE_BUFFER_DESC_FULL_MASK || component_mask == PICTURE_BUFFER_DESC_CHROMA_MASK) {
#if ATB_EP
        if (cu_ptr->prediction_mode_flag == INTRA_MODE && (context_ptr->evaluate_cfl_ep || cu_ptr->prediction_unit_array->intra_chroma_mode == UV_CFL_PRED)) {
            EbPictureBufferDesc *reconSamples = predSamples;
            uint32_t reconLumaOffset = (reconSamples->origin_y + round_origin_y) * reconSamples->stride_y + (reconSamples->origin_x + round_origin_x);

            // Down sample Luma
            cfl_luma_subsampling_420_lbd_c(
                reconSamples->buffer_y + reconLumaOffset,
                reconSamples->stride_y,
                context_ptr->md_context->pred_buf_q3,
#if CFL_FIX
                context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bwidth ? (context_ptr->blk_geom->bwidth_uv << 1) : context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight_uv == context_ptr->blk_geom->bheight ? (context_ptr->blk_geom->bheight_uv << 1) : context_ptr->blk_geom->bheight);
#else
                context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);
#endif
#if ATB_SUPPORT
            int32_t round_offset = ((context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr])*(context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr])) / 2;
#else
            int32_t round_offset = ((context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr])*(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr])) / 2;
#endif

            subtract_average(
                context_ptr->md_context->pred_buf_q3,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                round_offset,
                LOG2F(context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr]) + LOG2F(context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]));

#else
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                round_offset,
                LOG2F(context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr]) + LOG2F(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]));
#endif

            if (context_ptr->evaluate_cfl_ep)
            {
                // 3: Loop over alphas and find the best or choose DC
                // Use the 1st spot of the candidate buffer to hold cfl settings: (1) to use same kernel as MD for CFL evaluation: cfl_rd_pick_alpha() (toward unification), (2) to avoid dedicated buffers for CFL evaluation @ EP (toward less memory)
                ModeDecisionCandidateBuffer  *candidateBuffer = &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]);

                // Input(s)
                candidateBuffer->candidate_ptr->type = INTRA_MODE;
                candidateBuffer->candidate_ptr->intra_luma_mode = cu_ptr->pred_mode;
                candidateBuffer->candidate_ptr->cfl_alpha_signs = 0;
                candidateBuffer->candidate_ptr->cfl_alpha_idx = 0;
                context_ptr->md_context->blk_geom = context_ptr->blk_geom;

                EbByte src_pred_ptr;
                EbByte dst_pred_ptr;

                // Copy Cb pred samples from ep buffer to md buffer
                src_pred_ptr = predSamples->buffer_cb + pred_cb_offset;
                dst_pred_ptr = &(candidateBuffer->prediction_ptr->buffer_cb[scratch_cb_offset]);
                for (int i = 0; i < context_ptr->blk_geom->bheight_uv; i++) {
                    memcpy(dst_pred_ptr, src_pred_ptr, context_ptr->blk_geom->bwidth_uv);
                    src_pred_ptr += predSamples->stride_cb;
                    dst_pred_ptr += candidateBuffer->prediction_ptr->stride_cb;
                }

                // Copy Cr pred samples from ep buffer to md buffer
                src_pred_ptr = predSamples->buffer_cr + pred_cr_offset;
                dst_pred_ptr = &(candidateBuffer->prediction_ptr->buffer_cr[scratch_cr_offset]);
                for (int i = 0; i < context_ptr->blk_geom->bheight_uv; i++) {
                    memcpy(dst_pred_ptr, src_pred_ptr, context_ptr->blk_geom->bwidth_uv);
                    src_pred_ptr += predSamples->stride_cr;
                    dst_pred_ptr += candidateBuffer->prediction_ptr->stride_cr;
                }

                cfl_rd_pick_alpha(
                    picture_control_set_ptr,
                    candidateBuffer,
                    sb_ptr,
                    context_ptr->md_context,
                    input_samples,
                    input_cb_offset,
                    scratch_cb_offset,
                    asm_type);

                // Output(s)
                if (candidateBuffer->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {
                    cu_ptr->prediction_unit_array->intra_chroma_mode = UV_CFL_PRED;
                    cu_ptr->prediction_unit_array->cfl_alpha_idx = candidateBuffer->candidate_ptr->cfl_alpha_idx;
                    cu_ptr->prediction_unit_array->cfl_alpha_signs = candidateBuffer->candidate_ptr->cfl_alpha_signs;
                    cu_ptr->prediction_unit_array->is_directional_chroma_mode_flag = EB_FALSE;
                }
            }

            if (cu_ptr->prediction_unit_array->intra_chroma_mode == UV_CFL_PRED) {
                int32_t alpha_q3 =
                    cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_U); // once for U, once for V

                //TOCHANGE
                //assert(chroma_size * CFL_BUF_LINE + chroma_size <= CFL_BUF_SQUARE);

                cfl_predict_lbd(
                    context_ptr->md_context->pred_buf_q3,
                    predSamples->buffer_cb + pred_cb_offset,
                    predSamples->stride_cb,
                    predSamples->buffer_cb + pred_cb_offset,
                    predSamples->stride_cb,
                    alpha_q3,
                    8,
#if ATB_SUPPORT
                    context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
                    context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif
                alpha_q3 =
                    cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_V); // once for U, once for V

                //TOCHANGE
                //assert(chroma_size * CFL_BUF_LINE + chroma_size <= CFL_BUF_SQUARE);

                cfl_predict_lbd(
                    context_ptr->md_context->pred_buf_q3,
                    predSamples->buffer_cr + pred_cr_offset,
                    predSamples->stride_cr,
                    predSamples->buffer_cr + pred_cr_offset,
                    predSamples->stride_cr,
                    alpha_q3,
                    8,
#if ATB_SUPPORT
                    context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
                    context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif
            }
        }
#endif

        //**********************************
        // Cb
        //**********************************

        ResidualKernel(
            input_samples->buffer_cb + input_cb_offset,
            input_samples->stride_cb,
            predSamples->buffer_cb + pred_cb_offset,
            predSamples->stride_cb,
            ((int16_t*)residual16bit->buffer_cb) + scratch_cb_offset,
            residual16bit->stride_cb,
#if ATB_SUPPORT
            context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
            context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif

        ResidualKernel(
            input_samples->buffer_cr + input_cr_offset,
            input_samples->stride_cr,
            predSamples->buffer_cr + pred_cr_offset,
            predSamples->stride_cr,
            ((int16_t*)residual16bit->buffer_cr) + scratch_cr_offset,
            residual16bit->stride_cr,
#if ATB_SUPPORT
            context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
            context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif

        av1_estimate_transform(
            ((int16_t*)residual16bit->buffer_cb) + scratch_cb_offset,
            residual16bit->stride_cb,
            ((TranLow*)transform16bit->buffer_cb) + context_ptr->coded_area_sb_uv,
            NOT_USED_VALUE,
#if ATB_SUPPORT
            context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_UV],
            asm_type,
            PLANE_TYPE_UV,
#if PF_N2_SUPPORT
            DEFAULT_SHAPE);
#else
            context_ptr->trans_coeff_shape_chroma);
#endif
#if DC_SIGN_CONTEXT_EP
        cu_ptr->quantized_dc[1][context_ptr->txb_itr] = av1_quantize_inv_quantize(
#else
        av1_quantize_inv_quantize(
#endif
            sb_ptr->picture_control_set_ptr,
            context_ptr->md_context,
            ((TranLow*)transform16bit->buffer_cb) + context_ptr->coded_area_sb_uv,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->buffer_cb) + context_ptr->coded_area_sb_uv,
            ((int32_t*)inverse_quant_buffer->buffer_cb) + context_ptr->coded_area_sb_uv,
            qp,
#if ATB_SUPPORT
            context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
            context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
            &eob[1],
            asm_type,
            &(count_non_zero_coeffs[1]),
#if !PF_N2_SUPPORT
            0,
#endif
            COMPONENT_CHROMA_CB,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_UV],
            &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]),
            cu_ptr->cb_txb_skip_context,
            cu_ptr->cb_dc_sign_context,
            cu_ptr->pred_mode,
            EB_TRUE);

#if BLK_SKIP_DECISION
        if (context_ptr->md_skip_blk) {
            count_non_zero_coeffs[1] = 0;
            eob[1] = 0;
        }
#endif
        txb_ptr->u_has_coeff = count_non_zero_coeffs[1] ? EB_TRUE : EB_FALSE;

        //**********************************
        // Cr
        //**********************************

        av1_estimate_transform(
            ((int16_t*)residual16bit->buffer_cr) + scratch_cb_offset,
            residual16bit->stride_cr,
            ((TranLow*)transform16bit->buffer_cr) + context_ptr->coded_area_sb_uv,
            NOT_USED_VALUE,
#if ATB_SUPPORT
            context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_UV],
            asm_type,
            PLANE_TYPE_UV,
#if PF_N2_SUPPORT
            DEFAULT_SHAPE);
#else
            context_ptr->trans_coeff_shape_chroma);
#endif
#if DC_SIGN_CONTEXT_EP
        cu_ptr->quantized_dc[2][context_ptr->txb_itr] = av1_quantize_inv_quantize(
#else
        av1_quantize_inv_quantize(
#endif
            sb_ptr->picture_control_set_ptr,
            context_ptr->md_context,
            ((TranLow*)transform16bit->buffer_cr) + context_ptr->coded_area_sb_uv,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->buffer_cr) + context_ptr->coded_area_sb_uv,
            ((TranLow*)inverse_quant_buffer->buffer_cr) + context_ptr->coded_area_sb_uv,
            qp,
#if ATB_SUPPORT
            context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
            context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
            &eob[2],
            asm_type,
            &(count_non_zero_coeffs[2]),
#if !PF_N2_SUPPORT
            0,
#endif
            COMPONENT_CHROMA_CR,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_UV],
            &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]),
            cu_ptr->cr_txb_skip_context,
            cu_ptr->cr_dc_sign_context,
            cu_ptr->pred_mode,
            EB_TRUE);
#if BLK_SKIP_DECISION
        if (context_ptr->md_skip_blk) {
            count_non_zero_coeffs[2] = 0;
            eob[2] = 0;
        }
#endif
        txb_ptr->v_has_coeff = count_non_zero_coeffs[2] ? EB_TRUE : EB_FALSE;

#if ATB_EP
        txb_ptr->nz_coef_count[1] = (uint16_t)count_non_zero_coeffs[1];
        txb_ptr->nz_coef_count[2] = (uint16_t)count_non_zero_coeffs[2];
#endif
    }
#if !PF_N2_SUPPORT
    txb_ptr->trans_coeff_shape_luma = context_ptr->trans_coeff_shape_luma;
    txb_ptr->trans_coeff_shape_chroma = context_ptr->trans_coeff_shape_chroma;
#endif
#if !ATB_EP
    txb_ptr->nz_coef_count[0] = (uint16_t)count_non_zero_coeffs[0];
    txb_ptr->nz_coef_count[1] = (uint16_t)count_non_zero_coeffs[1];
    txb_ptr->nz_coef_count[2] = (uint16_t)count_non_zero_coeffs[2];
#endif
    return;
}

void encode_pass_tx_search_hbd(
    PictureControlSet            *picture_control_set_ptr,
    EncDecContext                *context_ptr,
    LargestCodingUnit            *sb_ptr,
    uint32_t                       cb_qp,
    EbPictureBufferDesc          *coeffSamplesTB,
    EbPictureBufferDesc          *residual16bit,
    EbPictureBufferDesc          *transform16bit,
    EbPictureBufferDesc          *inverse_quant_buffer,
    int16_t                        *transformScratchBuffer,
    EbAsm                          asm_type,
    uint32_t                       *count_non_zero_coeffs,
    uint32_t                       component_mask,
    uint32_t                       use_delta_qp,
    uint32_t                       dZoffset,
    uint16_t                       *eob,
    MacroblockPlane                *candidate_plane);

/**********************************************************
* Encode Loop
*
* Summary: Performs a H.265 conformant
*   Transform, Quantization  and Inverse Quantization of a TU.
*
* Inputs:
*   origin_x
*   origin_y
*   txb_size
*   sb_sz
*   input - input samples (position sensitive)
*   pred - prediction samples (position independent)
*
* Outputs:
*   Inverse quantized coeff - quantization indices (position sensitive)
*
**********************************************************/
static void Av1EncodeLoop16bit(
    PictureControlSet    *picture_control_set_ptr,
    EncDecContext       *context_ptr,
    LargestCodingUnit   *sb_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    uint32_t                 cb_qp,
    EbPictureBufferDesc *predSamples,         // no basis/offset
    EbPictureBufferDesc *coeffSamplesTB,      // lcu based
    EbPictureBufferDesc *residual16bit,       // no basis/offset
    EbPictureBufferDesc *transform16bit,      // no basis/offset
    EbPictureBufferDesc *inverse_quant_buffer,
    int16_t                *transformScratchBuffer,
    EbAsm                 asm_type,
    uint32_t                  *count_non_zero_coeffs,
    uint32_t                 component_mask,
    uint32_t                   use_delta_qp,
    uint32_t                 dZoffset,
    uint16_t                 *eob,
    MacroblockPlane       *candidate_plane)

{
    (void)use_delta_qp;
    (void)dZoffset;
    (void)cb_qp;

    CodingUnit          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    //    EB_SLICE               slice_type = sb_ptr->picture_control_set_ptr->slice_type;
    //    uint32_t                 temporal_layer_index = sb_ptr->picture_control_set_ptr->temporal_layer_index;
    uint32_t                 qp = cu_ptr->qp;

    EbPictureBufferDesc *inputSamples16bit = context_ptr->input_sample16bit_buffer;
    EbPictureBufferDesc *predSamples16bit = predSamples;
#if ATB_SUPPORT
    uint32_t round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
    uint32_t round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4
    const uint32_t input_luma_offset = context_ptr->blk_geom->tx_org_x[cu_ptr->tx_depth][context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[cu_ptr->tx_depth][context_ptr->txb_itr] * SB_STRIDE_Y;
    const uint32_t input_cb_offset = ROUND_UV(context_ptr->blk_geom->tx_org_x[cu_ptr->tx_depth][context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[cu_ptr->tx_depth][context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t input_cr_offset = ROUND_UV(context_ptr->blk_geom->tx_org_x[cu_ptr->tx_depth][context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[cu_ptr->tx_depth][context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t pred_luma_offset = ((predSamples16bit->origin_y + origin_y)        * predSamples16bit->stride_y) + (predSamples16bit->origin_x + origin_x);
    const uint32_t pred_cb_offset = (((predSamples16bit->origin_y + round_origin_y) >> 1)  * predSamples16bit->stride_cb) + ((predSamples16bit->origin_x + round_origin_x) >> 1);
    const uint32_t pred_cr_offset = (((predSamples16bit->origin_y + round_origin_y) >> 1)  * predSamples16bit->stride_cr) + ((predSamples16bit->origin_x + round_origin_x) >> 1);
    const uint32_t scratch_luma_offset = context_ptr->blk_geom->origin_x + context_ptr->blk_geom->origin_y * SB_STRIDE_Y;
    const uint32_t scratch_cb_offset = ROUND_UV(context_ptr->blk_geom->origin_x) / 2 + ROUND_UV(context_ptr->blk_geom->origin_y) / 2 * SB_STRIDE_UV;
    const uint32_t scratch_cr_offset = ROUND_UV(context_ptr->blk_geom->origin_x) / 2 + ROUND_UV(context_ptr->blk_geom->origin_y) / 2 * SB_STRIDE_UV;
#else
    uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
    uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4
    const uint32_t           input_luma_offset = context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr] * SB_STRIDE_Y;
    const uint32_t           input_cb_offset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t           input_cr_offset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t           pred_luma_offset = ((predSamples16bit->origin_y + origin_y)        * predSamples16bit->stride_y) + (predSamples16bit->origin_x + origin_x);
    const uint32_t           pred_cb_offset = (((predSamples16bit->origin_y + round_origin_y) >> 1)  * predSamples16bit->stride_cb) + ((predSamples16bit->origin_x + round_origin_x) >> 1);
    const uint32_t           pred_cr_offset = (((predSamples16bit->origin_y + round_origin_y) >> 1)  * predSamples16bit->stride_cr) + ((predSamples16bit->origin_x + round_origin_x) >> 1);
    const uint32_t scratch_luma_offset = context_ptr->blk_geom->origin_x + context_ptr->blk_geom->origin_y * SB_STRIDE_Y;
    const uint32_t scratch_cb_offset = ROUND_UV(context_ptr->blk_geom->origin_x) / 2 + ROUND_UV(context_ptr->blk_geom->origin_y) / 2 * SB_STRIDE_UV;
    const uint32_t scratch_cr_offset = ROUND_UV(context_ptr->blk_geom->origin_x) / 2 + ROUND_UV(context_ptr->blk_geom->origin_y) / 2 * SB_STRIDE_UV;
#endif
    const uint32_t coeff1dOffset = context_ptr->coded_area_sb;
    const uint32_t coeff1dOffsetChroma = context_ptr->coded_area_sb_uv;
    UNUSED(coeff1dOffsetChroma);

    //Update QP for Quant
    qp += QP_BD_OFFSET;

    {
        //**********************************
        // Luma
        //**********************************
        if (component_mask == PICTURE_BUFFER_DESC_FULL_MASK || component_mask == PICTURE_BUFFER_DESC_LUMA_MASK) {
            residual_kernel16bit(
                ((uint16_t*)inputSamples16bit->buffer_y) + input_luma_offset,
                inputSamples16bit->stride_y,
                ((uint16_t*)predSamples16bit->buffer_y) + pred_luma_offset,
                predSamples16bit->stride_y,
                ((int16_t*)residual16bit->buffer_y) + scratch_luma_offset,
                residual16bit->stride_y,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
                context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);
#endif
#if ATB_MD
            uint8_t  tx_search_skip_fag = (picture_control_set_ptr->parent_pcs_ptr->tx_search_level == TX_SEARCH_ENC_DEC && (picture_control_set_ptr->parent_pcs_ptr->atb_mode == 0 || cu_ptr->prediction_mode_flag == INTER_MODE)) ? get_skip_tx_search_flag(
#else
            uint8_t  tx_search_skip_fag = picture_control_set_ptr->parent_pcs_ptr->tx_search_level == TX_SEARCH_ENC_DEC ? get_skip_tx_search_flag(
#endif
#if BYPASS_USELESS_TX_SEARCH
                context_ptr->blk_geom,
#else
                context_ptr->blk_geom->sq_size,
#endif
                MAX_MODE_COST,
                0,
                1) : 1;

            if (!tx_search_skip_fag) {
                    encode_pass_tx_search_hbd(
                        picture_control_set_ptr,
                        context_ptr,
                        sb_ptr,
                        cb_qp,
                        coeffSamplesTB,
                        residual16bit,
                        transform16bit,
                        inverse_quant_buffer,
                        transformScratchBuffer,
                        asm_type,
                        count_non_zero_coeffs,
                        component_mask,
                        use_delta_qp,
                        dZoffset,
                        eob,
                        candidate_plane);
            }

            av1_estimate_transform(
                ((int16_t*)residual16bit->buffer_y) + scratch_luma_offset,
                residual16bit->stride_y,
                ((TranLow*)transform16bit->buffer_y) + coeff1dOffset,
                NOT_USED_VALUE,
#if ATB_SUPPORT
                context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->txsize[context_ptr->txb_itr],
#endif
                &context_ptr->three_quad_energy,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_Y],
                asm_type,
                PLANE_TYPE_Y,
#if PF_N2_SUPPORT
                DEFAULT_SHAPE);
#else
                context_ptr->trans_coeff_shape_luma);
#endif
#if DC_SIGN_CONTEXT_EP
            cu_ptr->quantized_dc[0][context_ptr->txb_itr] = av1_quantize_inv_quantize(
#else
            av1_quantize_inv_quantize(
#endif
                sb_ptr->picture_control_set_ptr,
                context_ptr->md_context,
                ((int32_t*)transform16bit->buffer_y) + coeff1dOffset,
                NOT_USED_VALUE,
                ((int32_t*)coeffSamplesTB->buffer_y) + coeff1dOffset,
                ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
                qp,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->txb_itr],
                context_ptr->blk_geom->txsize[context_ptr->txb_itr],
#endif
                &eob[0],
                asm_type,
                &(count_non_zero_coeffs[0]),
#if !PF_N2_SUPPORT
                0,
#endif
                COMPONENT_LUMA,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_Y],
                &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]),
                cu_ptr->luma_txb_skip_context,
#if ATB_DC_CONTEXT_SUPPORT_0
                cu_ptr->luma_dc_sign_context[context_ptr->txb_itr],
#else
                cu_ptr->luma_dc_sign_context,
#endif
                cu_ptr->pred_mode,
                EB_TRUE);
#if BLK_SKIP_DECISION
            if (context_ptr->md_skip_blk) {
                count_non_zero_coeffs[0] = 0;
                eob[0] = 0;
            }
#endif
            txb_ptr->y_has_coeff = count_non_zero_coeffs[0] ? EB_TRUE : EB_FALSE;
            if (count_non_zero_coeffs[0] == 0) {
                // INTER. Chroma follows Luma in transform type
                if (cu_ptr->prediction_mode_flag == INTER_MODE) {
                    txb_ptr->transform_type[PLANE_TYPE_Y] = DCT_DCT;
                    txb_ptr->transform_type[PLANE_TYPE_UV] = DCT_DCT;
                }
                else { // INTRA
                    txb_ptr->transform_type[PLANE_TYPE_Y] = DCT_DCT;
                }
            }

#if ATB_EP
            txb_ptr->nz_coef_count[0] = (uint16_t)count_non_zero_coeffs[0];
#endif
        }

        if (cu_ptr->prediction_mode_flag == INTRA_MODE && cu_ptr->prediction_unit_array->intra_chroma_mode == UV_CFL_PRED) {
            EbPictureBufferDesc *reconSamples = predSamples16bit;
            uint32_t reconLumaOffset = (reconSamples->origin_y + origin_y)            * reconSamples->stride_y + (reconSamples->origin_x + origin_x);
            if (txb_ptr->y_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
                uint16_t     *predBuffer = ((uint16_t*)predSamples16bit->buffer_y) + pred_luma_offset;
                av1_inv_transform_recon(
                    ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
                    CONVERT_TO_BYTEPTR(predBuffer),
                    predSamples->stride_y,
#if ATB_SUPPORT
                    context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                    context_ptr->blk_geom->txsize[context_ptr->txb_itr],
#endif
                    BIT_INCREMENT_10BIT,
                    txb_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    eob[0]);
            }
#if CFL_FIX
            if (context_ptr->blk_geom->has_uv) {
                reconLumaOffset = (reconSamples->origin_y + round_origin_y)            * reconSamples->stride_y + (reconSamples->origin_x + round_origin_x);
#endif
            // Down sample Luma
            cfl_luma_subsampling_420_hbd_c(
                ((uint16_t*)reconSamples->buffer_y) + reconLumaOffset,
                reconSamples->stride_y,
                context_ptr->md_context->pred_buf_q3,
#if CFL_FIX
                context_ptr->blk_geom->bwidth_uv == context_ptr->blk_geom->bwidth ? (context_ptr->blk_geom->bwidth_uv << 1) : context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight_uv == context_ptr->blk_geom->bheight ? (context_ptr->blk_geom->bheight_uv << 1) : context_ptr->blk_geom->bheight);
#else
                context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);
#endif
#if ATB_SUPPORT
            int32_t round_offset = ((context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr])*(context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr])) / 2;
#else
            int32_t round_offset = ((context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr])*(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr])) / 2;
#endif

            subtract_average(
                context_ptr->md_context->pred_buf_q3,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                round_offset,
                LOG2F(context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr]) + LOG2F(context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]));

#else
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                round_offset,
                LOG2F(context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr]) + LOG2F(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]));
#endif

            int32_t alpha_q3 =
                cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_U); // once for U, once for V
            // TOCHANGE
            // assert(chroma_size * CFL_BUF_LINE + chroma_size <=                CFL_BUF_SQUARE);

            cfl_predict_hbd(
                context_ptr->md_context->pred_buf_q3,
                ((uint16_t*)predSamples16bit->buffer_cb) + pred_cb_offset,
                predSamples16bit->stride_cb,
                ((uint16_t*)predSamples16bit->buffer_cb) + pred_cb_offset,
                predSamples16bit->stride_cb,
                alpha_q3,
                10,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif

            alpha_q3 =
                cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_V); // once for U, once for V
            // TOCHANGE
            //assert(chroma_size * CFL_BUF_LINE + chroma_size <=                CFL_BUF_SQUARE);

            cfl_predict_hbd(
                context_ptr->md_context->pred_buf_q3,
                ((uint16_t*)predSamples16bit->buffer_cr) + pred_cr_offset,
                predSamples16bit->stride_cr,
                ((uint16_t*)predSamples16bit->buffer_cr) + pred_cr_offset,
                predSamples16bit->stride_cr,
                alpha_q3,
                10,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif
        }
#if CFL_FIX
        }
#endif

        if (component_mask == PICTURE_BUFFER_DESC_FULL_MASK || component_mask == PICTURE_BUFFER_DESC_CHROMA_MASK) {
            //**********************************
            // Cb
            //**********************************
            residual_kernel16bit(
                ((uint16_t*)inputSamples16bit->buffer_cb) + input_cb_offset,
                inputSamples16bit->stride_cb,
                ((uint16_t*)predSamples16bit->buffer_cb) + pred_cb_offset,
                predSamples16bit->stride_cb,
                ((int16_t*)residual16bit->buffer_cb) + scratch_cb_offset,
                residual16bit->stride_cb,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif

            residual_kernel16bit(
                ((uint16_t*)inputSamples16bit->buffer_cr) + input_cr_offset,
                inputSamples16bit->stride_cr,
                ((uint16_t*)predSamples16bit->buffer_cr) + pred_cr_offset,
                predSamples16bit->stride_cr,
                ((int16_t*)residual16bit->buffer_cr) + scratch_cr_offset,
                residual16bit->stride_cr,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr]);
#else
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
#endif

            av1_estimate_transform(
                ((int16_t*)residual16bit->buffer_cb) + scratch_cb_offset,
                residual16bit->stride_cb,
                ((TranLow*)transform16bit->buffer_cb) + context_ptr->coded_area_sb_uv,
                NOT_USED_VALUE,
#if ATB_SUPPORT
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
                &context_ptr->three_quad_energy,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                asm_type,
                PLANE_TYPE_UV,
#if PF_N2_SUPPORT
                DEFAULT_SHAPE);
#else
                context_ptr->trans_coeff_shape_chroma);
#endif
#if DC_SIGN_CONTEXT_EP
            cu_ptr->quantized_dc[1][context_ptr->txb_itr] = av1_quantize_inv_quantize(
#else
            av1_quantize_inv_quantize(
#endif
                sb_ptr->picture_control_set_ptr,
                context_ptr->md_context,
                ((int32_t*)transform16bit->buffer_cb) + context_ptr->coded_area_sb_uv,
                NOT_USED_VALUE,
                ((int32_t*)coeffSamplesTB->buffer_cb) + context_ptr->coded_area_sb_uv,
                ((int32_t*)inverse_quant_buffer->buffer_cb) + context_ptr->coded_area_sb_uv,
                qp,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
                &eob[1],
                asm_type,
                &(count_non_zero_coeffs[1]),
#if !PF_N2_SUPPORT
                0,
#endif
                COMPONENT_CHROMA_CB,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]),
                cu_ptr->cb_txb_skip_context,
                cu_ptr->cb_dc_sign_context,
                cu_ptr->pred_mode,
                EB_TRUE);

#if BLK_SKIP_DECISION
            if (context_ptr->md_skip_blk) {
                count_non_zero_coeffs[1] = 0;
                eob[1] = 0;
            }
#endif
            txb_ptr->u_has_coeff = count_non_zero_coeffs[1] ? EB_TRUE : EB_FALSE;

            //**********************************
            // Cr
            //**********************************

            av1_estimate_transform(
                ((int16_t*)residual16bit->buffer_cr) + scratch_cb_offset,
                residual16bit->stride_cr,
                ((TranLow*)transform16bit->buffer_cr) + context_ptr->coded_area_sb_uv,
                NOT_USED_VALUE,
#if ATB_SUPPORT
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
                &context_ptr->three_quad_energy,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                asm_type,
                PLANE_TYPE_UV,
#if PF_N2_SUPPORT
                DEFAULT_SHAPE);
#else
                context_ptr->trans_coeff_shape_chroma);
#endif

#if DC_SIGN_CONTEXT_EP
            cu_ptr->quantized_dc[2][context_ptr->txb_itr] = av1_quantize_inv_quantize(
#else
            av1_quantize_inv_quantize(
#endif
                sb_ptr->picture_control_set_ptr,
                context_ptr->md_context,
                ((int32_t*)transform16bit->buffer_cr) + context_ptr->coded_area_sb_uv,
                NOT_USED_VALUE,
                ((int32_t*)coeffSamplesTB->buffer_cr) + context_ptr->coded_area_sb_uv,
                ((int32_t*)inverse_quant_buffer->buffer_cr) + context_ptr->coded_area_sb_uv,
                qp,
#if ATB_SUPPORT
                context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
                &eob[2],
                asm_type,
                &(count_non_zero_coeffs[2]),
#if !PF_N2_SUPPORT
                0,
#endif
                COMPONENT_CHROMA_CR,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]),
                cu_ptr->cr_txb_skip_context,
                cu_ptr->cr_dc_sign_context,
                cu_ptr->pred_mode,
                EB_TRUE);
#if BLK_SKIP_DECISION
            if (context_ptr->md_skip_blk) {
                count_non_zero_coeffs[2] = 0;
                eob[2] = 0;
            }
#endif
            txb_ptr->v_has_coeff = count_non_zero_coeffs[2] ? EB_TRUE : EB_FALSE;

#if ATB_EP
            txb_ptr->nz_coef_count[1] = (uint16_t)count_non_zero_coeffs[1];
            txb_ptr->nz_coef_count[2] = (uint16_t)count_non_zero_coeffs[2];
#endif
        }
    }

#if !PF_N2_SUPPORT
    txb_ptr->trans_coeff_shape_luma = context_ptr->trans_coeff_shape_luma;
    txb_ptr->trans_coeff_shape_chroma = context_ptr->trans_coeff_shape_chroma;
#endif

#if !ATB_EP
    txb_ptr->nz_coef_count[0] = (uint16_t)count_non_zero_coeffs[0];
    txb_ptr->nz_coef_count[1] = (uint16_t)count_non_zero_coeffs[1];
    txb_ptr->nz_coef_count[2] = (uint16_t)count_non_zero_coeffs[2];
#endif
    return;
}

/**********************************************************
* Encode Generate Recon
*
* Summary: Performs a H.265 conformant
*   Inverse Transform and generate
*   the reconstructed samples of a TU.
*
* Inputs:
*   origin_x
*   origin_y
*   txb_size
*   sb_sz
*   input - Inverse Qunatized Coeff (position sensitive)
*   pred - prediction samples (position independent)
*
* Outputs:
*   Recon  (position independent)
*
**********************************************************/
static void Av1EncodeGenerateRecon(
    EncDecContext       *context_ptr,
    uint32_t               origin_x,
    uint32_t               origin_y,
    EbPictureBufferDesc *predSamples,     // no basis/offset
    EbPictureBufferDesc *residual16bit,    // no basis/offset
    int16_t               *transformScratchBuffer,
    uint32_t               component_mask,
    uint16_t              *eob,
    EbAsm                  asm_type)
{
    uint32_t               pred_luma_offset;
    uint32_t               predChromaOffset;
    CodingUnit          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];

    // *Note - The prediction is built in-place in the Recon buffer. It is overwritten with Reconstructed
    //   samples if the CBF==1 && SKIP==False

    //**********************************
    // Luma
    //**********************************
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {
#if ATB_EP
        {
#else
        if (cu_ptr->prediction_mode_flag != INTRA_MODE || (cu_ptr->prediction_unit_array->intra_chroma_mode != UV_CFL_PRED && context_ptr->evaluate_cfl_ep == EB_FALSE))
        {
#endif
            pred_luma_offset = (predSamples->origin_y + origin_y)             * predSamples->stride_y + (predSamples->origin_x + origin_x);
            if (txb_ptr->y_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
                (void)asm_type;
                (void)transformScratchBuffer;
                uint8_t     *predBuffer = predSamples->buffer_y + pred_luma_offset;
                av1_inv_transform_recon8bit(
                    ((int32_t*)residual16bit->buffer_y) + context_ptr->coded_area_sb,
                    predBuffer,
                    predSamples->stride_y,
#if ATB_SUPPORT
                    context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                    context_ptr->blk_geom->txsize[context_ptr->txb_itr],
#endif
                    txb_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    eob[0]
                );
            }
        }
    }

    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK)    {
        //**********************************
        // Chroma
        //**********************************

        uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
        uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4
        predChromaOffset = (((predSamples->origin_y + round_origin_y) >> 1)           * predSamples->stride_cb) + ((predSamples->origin_x + round_origin_x) >> 1);

        //**********************************
        // Cb
        //**********************************
        if (txb_ptr->u_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
            uint8_t     *predBuffer = predSamples->buffer_cb + predChromaOffset;

            av1_inv_transform_recon8bit(
                ((int32_t*)residual16bit->buffer_cb) + context_ptr->coded_area_sb_uv,
                predBuffer,
                predSamples->stride_cb,
#if ATB_SUPPORT
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
                txb_ptr->transform_type[PLANE_TYPE_UV],
                PLANE_TYPE_UV,
                eob[1]);
        }

        //**********************************
        // Cr
        //**********************************
        predChromaOffset = (((predSamples->origin_y + round_origin_y) >> 1)           * predSamples->stride_cr) + ((predSamples->origin_x + round_origin_x) >> 1);

        if (txb_ptr->v_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
            uint8_t     *predBuffer = predSamples->buffer_cr + predChromaOffset;

            av1_inv_transform_recon8bit(
                ((int32_t*)residual16bit->buffer_cr) + context_ptr->coded_area_sb_uv,
                predBuffer,
                predSamples->stride_cr,
#if ATB_SUPPORT
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
                txb_ptr->transform_type[PLANE_TYPE_UV],
                PLANE_TYPE_UV,
                eob[2]);
        }
    }

    return;
}

/**********************************************************
* Encode Generate Recon
*
* Summary: Performs a H.265 conformant
*   Inverse Transform and generate
*   the reconstructed samples of a TU.
*
* Inputs:
*   origin_x
*   origin_y
*   txb_size
*   sb_sz
*   input - Inverse Qunatized Coeff (position sensitive)
*   pred - prediction samples (position independent)
*
* Outputs:
*   Recon  (position independent)
*
**********************************************************/
static void Av1EncodeGenerateRecon16bit(
    EncDecContext         *context_ptr,
    uint32_t               origin_x,
    uint32_t               origin_y,
    EbPictureBufferDesc   *predSamples,     // no basis/offset
    EbPictureBufferDesc   *residual16bit,    // no basis/offset
    int16_t               *transformScratchBuffer,
    uint32_t               component_mask,
    uint16_t              *eob,
    EbAsm                  asm_type)
{
    uint32_t pred_luma_offset;
    uint32_t predChromaOffset;

    CodingUnit          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];

    (void)asm_type;
    (void)transformScratchBuffer;
    //**********************************
    // Luma
    //**********************************
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {
        if (cu_ptr->prediction_mode_flag != INTRA_MODE || cu_ptr->prediction_unit_array->intra_chroma_mode != UV_CFL_PRED)

        {
            pred_luma_offset = (predSamples->origin_y + origin_y)* predSamples->stride_y + (predSamples->origin_x + origin_x);
            if (txb_ptr->y_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
                uint16_t     *predBuffer = ((uint16_t*)predSamples->buffer_y) + pred_luma_offset;
                av1_inv_transform_recon(
                    ((int32_t*)residual16bit->buffer_y) + context_ptr->coded_area_sb,
                    CONVERT_TO_BYTEPTR(predBuffer),
                    predSamples->stride_y,
#if ATB_SUPPORT
                    context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                    context_ptr->blk_geom->txsize[context_ptr->txb_itr],
#endif
                    BIT_INCREMENT_10BIT,
                    txb_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    eob[0]
                );
            }
        }
    }

    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK)    {
        //**********************************
        // Chroma
        //**********************************

        //**********************************
        // Cb
        //**********************************

        uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
        uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4

        predChromaOffset = (((predSamples->origin_y + round_origin_y) >> 1)           * predSamples->stride_cb) + ((predSamples->origin_x + round_origin_x) >> 1);

        if (txb_ptr->u_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
            uint16_t     *predBuffer = ((uint16_t*)predSamples->buffer_cb) + predChromaOffset;
            av1_inv_transform_recon(
                ((int32_t*)residual16bit->buffer_cb) + context_ptr->coded_area_sb_uv,
                CONVERT_TO_BYTEPTR(predBuffer),
                predSamples->stride_cb,
#if ATB_SUPPORT
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                PLANE_TYPE_UV,
                eob[1]);
        }

        //**********************************
        // Cr
        //**********************************
        predChromaOffset = (((predSamples->origin_y + round_origin_y) >> 1)           * predSamples->stride_cr) + ((predSamples->origin_x + round_origin_x) >> 1);
        if (txb_ptr->v_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
            uint16_t     *predBuffer = ((uint16_t*)predSamples->buffer_cr) + predChromaOffset;
            av1_inv_transform_recon(
                ((int32_t*)residual16bit->buffer_cr) + context_ptr->coded_area_sb_uv,
                CONVERT_TO_BYTEPTR(predBuffer),
                predSamples->stride_cr,
#if ATB_SUPPORT
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                PLANE_TYPE_UV,
                eob[2]);
        }
    }

    return;
}
static EB_AV1_ENCODE_LOOP_FUNC_PTR   Av1EncodeLoopFunctionTable[2] =
{
    Av1EncodeLoop,
    Av1EncodeLoop16bit
};

EB_AV1_GENERATE_RECON_FUNC_PTR   Av1EncodeGenerateReconFunctionPtr[2] =
{
    Av1EncodeGenerateRecon,
    Av1EncodeGenerateRecon16bit
};
#if !MEMORY_FOOTPRINT_OPT
/*******************************************
* Encode Pass - Assign Delta Qp
*******************************************/
static void EncodePassUpdateQp(
    PictureControlSet     *picture_control_set_ptr,
    EncDecContext         *context_ptr,
    EbBool                  availableCoeff,
    EbBool                  isDeltaQpEnable,
    EbBool                 *isDeltaQpNotCoded,
    uint32_t                   dif_cu_delta_qp_depth,
    uint8_t                   *prev_coded_qp,
    uint8_t                   *prev_quant_group_coded_qp,
    uint32_t                   sb_qp
)
{
    uint32_t ref_qp;
    uint8_t qp;

    uint32_t  log2MinCuQpDeltaSize = LOG2F_MAX_LCU_SIZE - dif_cu_delta_qp_depth;
    int32_t  qpTopNeighbor = 0;
    int32_t  qpLeftNeighbor = 0;
    EbBool newQuantGroup;
    uint32_t  quantGroupX = context_ptr->cu_origin_x - (context_ptr->cu_origin_x & ((1 << log2MinCuQpDeltaSize) - 1));
    uint32_t  quantGroupY = context_ptr->cu_origin_y - (context_ptr->cu_origin_y & ((1 << log2MinCuQpDeltaSize) - 1));
    EbBool sameLcuCheckTop = (((quantGroupY - 1) >> LOG2F_MAX_LCU_SIZE) == ((quantGroupY) >> LOG2F_MAX_LCU_SIZE)) ? EB_TRUE : EB_FALSE;
    EbBool sameLcuCheckLeft = (((quantGroupX - 1) >> LOG2F_MAX_LCU_SIZE) == ((quantGroupX) >> LOG2F_MAX_LCU_SIZE)) ? EB_TRUE : EB_FALSE;
    // Neighbor Array
    uint32_t qpLeftNeighborIndex = 0;
    uint32_t qpTopNeighborIndex = 0;

    // CU larger than the quantization group

    if (Log2f(context_ptr->blk_geom->bwidth) >= log2MinCuQpDeltaSize)
        *isDeltaQpNotCoded = EB_TRUE;
    // At the beginning of a new quantization group
    if (((context_ptr->cu_origin_x & ((1 << log2MinCuQpDeltaSize) - 1)) == 0) &&
        ((context_ptr->cu_origin_y & ((1 << log2MinCuQpDeltaSize) - 1)) == 0))
    {
        *isDeltaQpNotCoded = EB_TRUE;
        newQuantGroup = EB_TRUE;
    }
    else
        newQuantGroup = EB_FALSE;
    // setting the previous Quantization Group QP
    if (newQuantGroup == EB_TRUE)
        *prev_coded_qp = *prev_quant_group_coded_qp;
    if (sameLcuCheckTop) {
        qpTopNeighborIndex =
            LUMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(
                quantGroupX,
                quantGroupY - 1,
                picture_control_set_ptr->qp_array_stride);
        qpTopNeighbor = picture_control_set_ptr->qp_array[qpTopNeighborIndex];
    }
    else
        qpTopNeighbor = *prev_coded_qp;
    if (sameLcuCheckLeft) {
        qpLeftNeighborIndex =
            LUMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(
                quantGroupX - 1,
                quantGroupY,
                picture_control_set_ptr->qp_array_stride);

        qpLeftNeighbor = picture_control_set_ptr->qp_array[qpLeftNeighborIndex];
    }
    else
        qpLeftNeighbor = *prev_coded_qp;
    ref_qp = (qpLeftNeighbor + qpTopNeighbor + 1) >> 1;

    qp = (uint8_t)context_ptr->cu_ptr->qp;
    // Update the State info
    if (isDeltaQpEnable) {
        if (*isDeltaQpNotCoded) {
            if (availableCoeff) {
                qp = (uint8_t)context_ptr->cu_ptr->qp;
                *prev_coded_qp = qp;
                *prev_quant_group_coded_qp = qp;
                *isDeltaQpNotCoded = EB_FALSE;
            }
            else {
                qp = (uint8_t)ref_qp;
                *prev_quant_group_coded_qp = qp;
            }
        }
    }
    else
        qp = (uint8_t)sb_qp;
    context_ptr->cu_ptr->qp = qp;
    return;
}
#endif

#if !MEMORY_FOOTPRINT_OPT
EbErrorType QpmDeriveBeaAndSkipQpmFlagLcu(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    LargestCodingUnit                    *sb_ptr,
    uint32_t                                 sb_index,
    EncDecContext                        *context_ptr)
{
    EbErrorType                    return_error = EB_ErrorNone;
#if ADD_DELTA_QP_SUPPORT
    uint16_t                           picture_qp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    uint16_t                           min_qp_allowed = 0;
    uint16_t                           max_qp_allowed = 255;
    uint16_t                           deltaQpRes = (uint16_t)picture_control_set_ptr->parent_pcs_ptr->delta_q_res;
#else
    uint8_t                           picture_qp = picture_control_set_ptr->picture_qp;

    uint8_t                           min_qp_allowed = (uint8_t)sequence_control_set_ptr->static_config.min_qp_allowed;
    uint8_t                           max_qp_allowed = (uint8_t)sequence_control_set_ptr->static_config.max_qp_allowed;
#endif

    context_ptr->qpm_qp = picture_qp;

    SbStat *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);

    context_ptr->non_moving_delta_qp = 0;

    context_ptr->grass_enhancement_flag = ((picture_control_set_ptr->scene_caracteristic_id == EB_FRAME_CARAC_1) && (sb_stat_ptr->cu_stat_array[0].grass_area)
        && (sb_ptr->picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index].edge_block_num > 0))

        ? EB_TRUE : EB_FALSE;

    context_ptr->backgorund_enhancement = EB_FALSE;

    context_ptr->skip_qpm_flag = sequence_control_set_ptr->static_config.improve_sharpness ? EB_FALSE : EB_TRUE;
#if !DISABLE_OIS_USE
    if ((picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag == EB_FALSE) && ((picture_control_set_ptr->parent_pcs_ptr->pic_noise_class >= PIC_NOISE_CLASS_3_1) || (picture_control_set_ptr->parent_pcs_ptr->high_dark_low_light_area_density_flag) || (picture_control_set_ptr->parent_pcs_ptr->intra_coded_block_probability > 90)))
        context_ptr->skip_qpm_flag = EB_TRUE;
#endif
    if (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE)
        context_ptr->skip_qpm_flag = EB_TRUE;
#if ADD_DELTA_QP_SUPPORT
    context_ptr->skip_qpm_flag = EB_FALSE;
#endif

    if (context_ptr->skip_qpm_flag == EB_FALSE) {
        if (picture_control_set_ptr->parent_pcs_ptr->pic_homogenous_over_time_sb_percentage > 30 && picture_control_set_ptr->slice_type != I_SLICE) {
#if ADD_DELTA_QP_SUPPORT
            context_ptr->qpm_qp = CLIP3(min_qp_allowed, max_qp_allowed, picture_qp + deltaQpRes);
#else
            context_ptr->qpm_qp = CLIP3(min_qp_allowed, max_qp_allowed, picture_qp + 1);
#endif
        }
    }

    return return_error;
}
#endif
#if ADD_DELTA_QP_SUPPORT
/*****************************************************************************
* NM - Note: Clean up
* AV1 QPM is SB based and all sub-Lcu buffers needs to be removed
******************************************************************************/
EbErrorType Av1QpModulationLcu(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    LargestCodingUnit                    *sb_ptr,
    uint32_t                                  sb_index,
    uint8_t                                   type,
    EncDecContext                        *context_ptr)
{
    EbErrorType                            return_error = EB_ErrorNone;

    int64_t                                  complexityDistance;
    int8_t                                   delta_qp = 0;
    uint16_t                                  qpm_qp = context_ptr->qpm_qp;
    uint16_t                                  min_qp_allowed = 0;
    uint16_t                                  max_qp_allowed = 255;
    uint16_t                                  cu_qp;
    EbBool                                 acEnergyBasedAntiContouring = picture_control_set_ptr->slice_type == I_SLICE ? EB_TRUE : EB_FALSE;
    uint8_t                                   lowerQPClass;

    int8_t    non_moving_delta_qp = context_ptr->non_moving_delta_qp;
    int8_t    bea64x64DeltaQp;

    uint8_t   deltaQpRes = picture_control_set_ptr->parent_pcs_ptr->delta_q_res;

    cu_qp = qpm_qp;
    sb_ptr->qp = qpm_qp;

    uint32_t  distortion = 0;

    if (!context_ptr->skip_qpm_flag) {
        // INTRA MODE
        if (type == INTRA_MODE) {
            OisSbResults        *ois_sb_results_ptr = picture_control_set_ptr->parent_pcs_ptr->ois_sb_results[sb_index];
            OisCandidate *OisCuPtr = ois_sb_results_ptr->sorted_ois_candidate[from_1101_to_85[cu_index]];
            distortion = OisCuPtr[ois_sb_results_ptr->best_distortion_index[from_1101_to_85[cu_index]]].distortion;

            distortion = (uint32_t)CLIP3(picture_control_set_ptr->parent_pcs_ptr->intra_complexity_min[0], picture_control_set_ptr->parent_pcs_ptr->intra_complexity_max[0], distortion);
            complexityDistance = ((int32_t)distortion - (int32_t)picture_control_set_ptr->parent_pcs_ptr->intra_complexity_avg[0]);

            if (complexityDistance < 0)
                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->intra_min_distance[0] != 0) ? (int8_t)((context_ptr->min_delta_qp_weight * context_ptr->min_delta_qp[0] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->intra_min_distance[0])) : 0;
            else
                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->intra_max_distance[0] != 0) ? (int8_t)((context_ptr->max_delta_qp_weight * context_ptr->max_delta_qp[0] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->intra_max_distance[0])) : 0;
        }
        // INTER MODE
        else {
            distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][0].distortion_direction[0].distortion;

            distortion = (uint32_t)CLIP3(picture_control_set_ptr->parent_pcs_ptr->inter_complexity_min[0], picture_control_set_ptr->parent_pcs_ptr->inter_complexity_max[0], distortion);
            complexityDistance = ((int32_t)distortion - (int32_t)picture_control_set_ptr->parent_pcs_ptr->inter_complexity_avg[0]);

            if (complexityDistance < 0)
                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->inter_min_distance[0] != 0) ? (int8_t)((context_ptr->min_delta_qp_weight * context_ptr->min_delta_qp[0] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->inter_min_distance[0])) : 0;
            else
                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->inter_max_distance[0] != 0) ? (int8_t)((context_ptr->max_delta_qp_weight * context_ptr->max_delta_qp[0] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->inter_max_distance[0])) : 0;
        }

        if (context_ptr->backgorund_enhancement) {
            // Use the 8x8 background enhancement only for the Intra slice, otherwise, use the existing SB based BEA results
            bea64x64DeltaQp = non_moving_delta_qp;

            if ((picture_control_set_ptr->parent_pcs_ptr->y_mean[sb_index][0] > ANTI_CONTOURING_LUMA_T2) || (picture_control_set_ptr->parent_pcs_ptr->y_mean[sb_index][0] < ANTI_CONTOURING_LUMA_T1)) {
                if (bea64x64DeltaQp < 0)
                    bea64x64DeltaQp = 0;
            }

            delta_qp += bea64x64DeltaQp;
        }

        if ((picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag))
            delta_qp = (delta_qp < context_ptr->min_delta_qp[0]) ? delta_qp : context_ptr->min_delta_qp[0];
        SbStat *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);
        if (sb_stat_ptr->stationary_edge_over_time_flag && delta_qp > 0)
            delta_qp = 0;
        if (acEnergyBasedAntiContouring) {
            lowerQPClass = derive_contouring_class(
                sb_ptr->picture_control_set_ptr->parent_pcs_ptr,
                sb_ptr->index,
                (uint8_t)1/*cu_index*/);

            if (lowerQPClass) {
                if (lowerQPClass == 3)
                    delta_qp = ANTI_CONTOURING_DELTA_QP_0;
                else if (lowerQPClass == 2)
                    delta_qp = ANTI_CONTOURING_DELTA_QP_1;
                else if (lowerQPClass == 1)
                    delta_qp = ANTI_CONTOURING_DELTA_QP_2;
            }
        }

        delta_qp -= context_ptr->grass_enhancement_flag ? 3 : 0;

        delta_qp *= deltaQpRes;

        if (sequence_control_set_ptr->static_config.rate_control_mode == 1 || sequence_control_set_ptr->static_config.rate_control_mode == 2) {
            if (qpm_qp > (RC_QPMOD_MAXQP * deltaQpRes))
                delta_qp = MIN(0, delta_qp);
            cu_qp = (uint32_t)(qpm_qp + delta_qp);

            if ((qpm_qp <= (RC_QPMOD_MAXQP *deltaQpRes))) {
                cu_qp = (uint8_t)CLIP3(
                    min_qp_allowed,
                    RC_QPMOD_MAXQP*deltaQpRes,
                    cu_qp);
            }
        }
        else
            cu_qp = (uint8_t)(qpm_qp + delta_qp);
        cu_qp = (uint8_t)CLIP3(
            min_qp_allowed,
            max_qp_allowed,
            cu_qp);
    }

    sb_ptr->qp = sequence_control_set_ptr->static_config.improve_sharpness ? cu_qp : qpm_qp;

    sb_ptr->delta_qp = (int16_t)sb_ptr->qp - (int16_t)qpm_qp;

    sb_ptr->org_delta_qp = sb_ptr->delta_qp;

    if (sb_ptr->delta_qp % deltaQpRes != 0)
        printf("Qpm_error: delta_qp must be multiplier of deltaQpRes\n");
    if (sb_ptr->qp == 0)
        printf("Qpm_error: qp must be greater than 0 when use_delta_q is ON\n");

    return return_error;
}

#endif
#if !MEMORY_FOOTPRINT_OPT
EbErrorType EncQpmDeriveDeltaQPForEachLeafLcu(
    SequenceControlSet                   *sequence_control_set_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    LargestCodingUnit                    *sb_ptr,
    uint32_t                                  sb_index,
    CodingUnit                           *cu_ptr,
    uint32_t                                  cu_depth,
    uint32_t                                  cu_index,
    uint32_t                                  cu_size,
    uint8_t                                   type,
    uint8_t                                   parent32x32_index,
    EncDecContext                        *context_ptr)
{
    EbErrorType                    return_error = EB_ErrorNone;

    //SbParams                        sb_params;
    int64_t                          complexityDistance;
    int8_t                           delta_qp = 0;
    uint8_t                           qpm_qp = (uint8_t)context_ptr->qpm_qp;
    uint8_t                           min_qp_allowed = (uint8_t)sequence_control_set_ptr->static_config.min_qp_allowed;
    uint8_t                           max_qp_allowed = (uint8_t)sequence_control_set_ptr->static_config.max_qp_allowed;
    uint8_t                           cu_qp;

    EbBool  use16x16Stat = EB_FALSE;

    uint32_t usedDepth = cu_depth;
    if (use16x16Stat)
        usedDepth = 2;

    uint32_t cuIndexInRaterScan = md_scan_to_raster_scan[cu_index];
#if !OPT_LOSSLESS_0
    EbBool                         acEnergyBasedAntiContouring = picture_control_set_ptr->slice_type == I_SLICE ? EB_TRUE : EB_FALSE;
    uint8_t                           lowerQPClass;
#endif

    int8_t    non_moving_delta_qp = context_ptr->non_moving_delta_qp;
    int8_t    bea64x64DeltaQp;

    cu_qp = qpm_qp;
    cu_ptr->qp = qpm_qp;

    uint32_t  distortion = 0;

    if (!context_ptr->skip_qpm_flag) {
        // INTRA MODE
        if (type == INTRA_MODE) {
            OisSbResults        *ois_sb_results_ptr = picture_control_set_ptr->parent_pcs_ptr->ois_sb_results[sb_index];
            OisCandidate *OisCuPtr = ois_sb_results_ptr->ois_candidate_array[ep_to_pa_block_index[cu_index]];
            distortion = OisCuPtr[ois_sb_results_ptr->best_distortion_index[ep_to_pa_block_index[cu_index]]].distortion;
            distortion = (uint32_t)CLIP3(picture_control_set_ptr->parent_pcs_ptr->intra_complexity_min[usedDepth], picture_control_set_ptr->parent_pcs_ptr->intra_complexity_max[usedDepth], distortion);
            complexityDistance = ((int32_t)distortion - (int32_t)picture_control_set_ptr->parent_pcs_ptr->intra_complexity_avg[usedDepth]);

            if (complexityDistance < 0)
                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->intra_min_distance[usedDepth] != 0) ? (int8_t)((context_ptr->min_delta_qp_weight * context_ptr->min_delta_qp[usedDepth] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->intra_min_distance[usedDepth])) : 0;
            else
                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->intra_max_distance[usedDepth] != 0) ? (int8_t)((context_ptr->max_delta_qp_weight * context_ptr->max_delta_qp[usedDepth] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->intra_max_distance[usedDepth])) : 0;
        }
        // INTER MODE
        else {
#if MRP_CONNECTION
            distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index]->me_candidate[cuIndexInRaterScan][0].distortion;
#else
            distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][cuIndexInRaterScan].distortion_direction[0].distortion;
#endif

            if (use16x16Stat) {
                uint32_t cuIndexRScan = md_scan_to_raster_scan[ParentBlockIndex[cu_index]];
#if MRP_CONNECTION
                distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index]->me_candidate[cuIndexRScan][0].distortion;
#else
                distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][cuIndexRScan].distortion_direction[0].distortion;
#endif
            }
            distortion = (uint32_t)CLIP3(picture_control_set_ptr->parent_pcs_ptr->inter_complexity_min[usedDepth], picture_control_set_ptr->parent_pcs_ptr->inter_complexity_max[usedDepth], distortion);
            complexityDistance = ((int32_t)distortion - (int32_t)picture_control_set_ptr->parent_pcs_ptr->inter_complexity_avg[usedDepth]);

            if (complexityDistance < 0)
                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->inter_min_distance[usedDepth] != 0) ? (int8_t)((context_ptr->min_delta_qp_weight * context_ptr->min_delta_qp[usedDepth] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->inter_min_distance[usedDepth])) : 0;
            else
                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->inter_max_distance[usedDepth] != 0) ? (int8_t)((context_ptr->max_delta_qp_weight * context_ptr->max_delta_qp[usedDepth] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->inter_max_distance[usedDepth])) : 0;
        }

        if (context_ptr->backgorund_enhancement) {
            // Use the 8x8 background enhancement only for the Intra slice, otherwise, use the existing SB based BEA results
            bea64x64DeltaQp = non_moving_delta_qp;

            if (((cu_index > 0) && ((picture_control_set_ptr->parent_pcs_ptr->y_mean[sb_index][parent32x32_index]) > ANTI_CONTOURING_LUMA_T2 || (picture_control_set_ptr->parent_pcs_ptr->y_mean[sb_index][parent32x32_index]) < ANTI_CONTOURING_LUMA_T1)) ||
                ((cu_index == 0) && ((picture_control_set_ptr->parent_pcs_ptr->y_mean[sb_index][0]) > ANTI_CONTOURING_LUMA_T2 || (picture_control_set_ptr->parent_pcs_ptr->y_mean[sb_index][0]) < ANTI_CONTOURING_LUMA_T1))) {
                if (bea64x64DeltaQp < 0)
                    bea64x64DeltaQp = 0;
            }

            delta_qp += bea64x64DeltaQp;
        }

        if ((picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag))
            delta_qp = (delta_qp < context_ptr->min_delta_qp[0]) ? delta_qp : context_ptr->min_delta_qp[0];
        SbStat *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);
        if (sb_stat_ptr->stationary_edge_over_time_flag && delta_qp > 0)
            delta_qp = 0;
#if !OPT_LOSSLESS_0
        if (acEnergyBasedAntiContouring) {
            lowerQPClass = derive_contouring_class(
                sb_ptr->picture_control_set_ptr->parent_pcs_ptr,
                sb_ptr->index,
                (uint8_t)cu_index);

            if (lowerQPClass) {
                if (lowerQPClass == 3)
                    delta_qp = ANTI_CONTOURING_DELTA_QP_0;
                else if (lowerQPClass == 2)
                    delta_qp = ANTI_CONTOURING_DELTA_QP_1;
                else if (lowerQPClass == 1)
                    delta_qp = ANTI_CONTOURING_DELTA_QP_2;
            }
        }
#endif

        delta_qp -= context_ptr->grass_enhancement_flag ? 3 : 0;

        if (sequence_control_set_ptr->static_config.rate_control_mode == 1 || sequence_control_set_ptr->static_config.rate_control_mode == 2) {
            if (qpm_qp > RC_QPMOD_MAXQP)
                delta_qp = MIN(0, delta_qp);
            cu_qp = (uint32_t)(qpm_qp + delta_qp);

            if ((qpm_qp <= RC_QPMOD_MAXQP)) {
                cu_qp = (uint8_t)CLIP3(
                    min_qp_allowed,
                    RC_QPMOD_MAXQP,
                    cu_qp);
            }
        }
        else
            cu_qp = (uint8_t)(qpm_qp + delta_qp);
        cu_qp = (uint8_t)CLIP3(
            min_qp_allowed,
            max_qp_allowed,
            cu_qp);
    }

    cu_ptr->qp = sequence_control_set_ptr->static_config.improve_sharpness ? cu_qp : qpm_qp;

    sb_ptr->qp = (cu_size == 64) ? (uint8_t)cu_ptr->qp : sb_ptr->qp;

    cu_ptr->delta_qp = (int16_t)cu_ptr->qp - (int16_t)qpm_qp;

    cu_ptr->org_delta_qp = cu_ptr->delta_qp;

    return return_error;
}
#endif
void Store16bitInputSrc(
    EncDecContext         *context_ptr,
    PictureControlSet     *picture_control_set_ptr,
    uint32_t                 lcuX,
    uint32_t                 lcuY,
    uint32_t                 lcuW,
    uint32_t                 lcuH ){
    uint32_t rowIt;
    uint16_t* fromPtr;
    uint16_t* toPtr;

    fromPtr = (uint16_t*)context_ptr->input_sample16bit_buffer->buffer_y;
    toPtr = (uint16_t*)picture_control_set_ptr->input_frame16bit->buffer_y + (lcuX + picture_control_set_ptr->input_frame16bit->origin_x) + (lcuY + picture_control_set_ptr->input_frame16bit->origin_y)*picture_control_set_ptr->input_frame16bit->stride_y;

    for (rowIt = 0; rowIt < lcuH; rowIt++)
        memcpy(toPtr + rowIt * picture_control_set_ptr->input_frame16bit->stride_y, fromPtr + rowIt * context_ptr->input_sample16bit_buffer->stride_y, lcuW * 2);
    lcuX = lcuX / 2;
    lcuY = lcuY / 2;
    lcuW = lcuW / 2;
    lcuH = lcuH / 2;

    fromPtr = (uint16_t*)context_ptr->input_sample16bit_buffer->buffer_cb;
    toPtr = (uint16_t*)picture_control_set_ptr->input_frame16bit->buffer_cb + (lcuX + picture_control_set_ptr->input_frame16bit->origin_x / 2) + (lcuY + picture_control_set_ptr->input_frame16bit->origin_y / 2)*picture_control_set_ptr->input_frame16bit->stride_cb;

    for (rowIt = 0; rowIt < lcuH; rowIt++)
        memcpy(toPtr + rowIt * picture_control_set_ptr->input_frame16bit->stride_cb, fromPtr + rowIt * context_ptr->input_sample16bit_buffer->stride_cb, lcuW * 2);
    fromPtr = (uint16_t*)context_ptr->input_sample16bit_buffer->buffer_cr;
    toPtr = (uint16_t*)picture_control_set_ptr->input_frame16bit->buffer_cr + (lcuX + picture_control_set_ptr->input_frame16bit->origin_x / 2) + (lcuY + picture_control_set_ptr->input_frame16bit->origin_y / 2)*picture_control_set_ptr->input_frame16bit->stride_cb;

    for (rowIt = 0; rowIt < lcuH; rowIt++)
        memcpy(toPtr + rowIt * picture_control_set_ptr->input_frame16bit->stride_cr, fromPtr + rowIt * context_ptr->input_sample16bit_buffer->stride_cr, lcuW * 2);
}

void update_av1_mi_map(
    CodingUnit        *cu_ptr,
    uint32_t           cu_origin_x,
    uint32_t           cu_origin_y,
    const BlockGeom   *blk_geom,
    PictureControlSet *picture_control_set_ptr);

void move_cu_data(
    CodingUnit *src_cu,
    CodingUnit *dst_cu);

#if ATB_EP
void perform_intra_coding_loop(
    SequenceControlSet *sequence_control_set_ptr,
    PictureControlSet  *picture_control_set_ptr,
    LargestCodingUnit  *sb_ptr,
    uint32_t            tbAddr,
    CodingUnit         *cu_ptr,
    PredictionUnit     *pu_ptr,
    EncDecContext      *context_ptr,
    uint32_t            dZoffset) {
    EbAsm                   asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
    EbBool                  is16bit = context_ptr->is16bit;

    EbPictureBufferDesc    *recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;
    EbPictureBufferDesc    *coeff_buffer_sb = sb_ptr->quantized_coeff;

    NeighborArrayUnit      *ep_luma_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_luma_recon_neighbor_array16bit : picture_control_set_ptr->ep_luma_recon_neighbor_array;
    NeighborArrayUnit      *ep_cb_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_cb_recon_neighbor_array16bit : picture_control_set_ptr->ep_cb_recon_neighbor_array;
    NeighborArrayUnit      *ep_cr_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_cr_recon_neighbor_array16bit : picture_control_set_ptr->ep_cr_recon_neighbor_array;

    EbPictureBufferDesc    *residual_buffer = context_ptr->residual_buffer;
    EbPictureBufferDesc    *transform_buffer = context_ptr->transform_buffer;
    EbPictureBufferDesc    *inverse_quant_buffer = context_ptr->inverse_quant_buffer;
    int16_t                *transform_inner_array_ptr = context_ptr->transform_inner_array_ptr;

    uint32_t                count_non_zero_coeffs[3];
    MacroblockPlane         cuPlane[3];
    uint16_t                eobs[MAX_TXB_COUNT][3];
    uint64_t                y_tu_coeff_bits;
    uint64_t                cb_tu_coeff_bits;
    uint64_t                cr_tu_coeff_bits;
    EntropyCoder           *coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;

    if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
        //get the 16bit form of the input LCU
        if (is16bit)
            recon_buffer = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            recon_buffer = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
    else  // non ref pictures
        recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;

    uint32_t totTu = context_ptr->blk_geom->txb_count[cu_ptr->tx_depth];

    // Luma path
    for (context_ptr->txb_itr = 0; context_ptr->txb_itr < totTu; context_ptr->txb_itr++) {
        uint16_t txb_origin_x = context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[cu_ptr->tx_depth][context_ptr->txb_itr];
        uint16_t txb_origin_y = context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[cu_ptr->tx_depth][context_ptr->txb_itr];

#if DC_SIGN_CONTEXT_EP
        context_ptr->cu_ptr->luma_txb_skip_context = 0;
        context_ptr->cu_ptr->luma_dc_sign_context[context_ptr->txb_itr] = 0;
        get_txb_ctx(
            COMPONENT_LUMA,
            picture_control_set_ptr->ep_luma_dc_sign_level_coeff_neighbor_array,
            txb_origin_x,
            txb_origin_y,
            context_ptr->blk_geom->bsize,
            context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
            &context_ptr->cu_ptr->luma_txb_skip_context,
            &context_ptr->cu_ptr->luma_dc_sign_context[context_ptr->txb_itr]);
#endif
        if (is16bit) {
            uint16_t    topNeighArray[64 * 2 + 1];
            uint16_t    leftNeighArray[64 * 2 + 1];
            PredictionMode mode;

            TxSize  tx_size = context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr];

            if (txb_origin_y != 0)
                memcpy(topNeighArray + 1, (uint16_t*)(ep_luma_recon_neighbor_array->top_array) + txb_origin_x, context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr] * 2 * sizeof(uint16_t));
            if (txb_origin_x != 0)
                memcpy(leftNeighArray + 1, (uint16_t*)(ep_luma_recon_neighbor_array->left_array) + txb_origin_y, context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr] * 2 * sizeof(uint16_t));
            if (txb_origin_y != 0 && txb_origin_x != 0)
                topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_luma_recon_neighbor_array->top_left_array) + MAX_PICTURE_HEIGHT_SIZE + txb_origin_x - txb_origin_y)[0];

            mode = cu_ptr->pred_mode;

            av1_predict_intra_block_16bit(
                &sb_ptr->tile_info,
                context_ptr,
                picture_control_set_ptr->parent_pcs_ptr->av1_cm,
                context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr],
                tx_size,
                mode,
                pu_ptr->angle_delta[PLANE_TYPE_Y],
                0,
                FILTER_INTRA_MODES,
                topNeighArray + 1,
                leftNeighArray + 1,
                recon_buffer,
                0,
                0,
                0,
                context_ptr->blk_geom->bsize,
                context_ptr->cu_origin_x,
                context_ptr->cu_origin_y);
        }
        else {
            uint8_t    topNeighArray[64 * 2 + 1];
            uint8_t    leftNeighArray[64 * 2 + 1];
            PredictionMode mode;

            TxSize  tx_size = context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr];

            if (txb_origin_y != 0)
                memcpy(topNeighArray + 1, ep_luma_recon_neighbor_array->top_array + txb_origin_x, context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr] * 2);

            if (txb_origin_x != 0)
                memcpy(leftNeighArray + 1, ep_luma_recon_neighbor_array->left_array + txb_origin_y, context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr] * 2);

            if (txb_origin_y != 0 && txb_origin_x != 0)
                topNeighArray[0] = leftNeighArray[0] = ep_luma_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE + txb_origin_x - txb_origin_y];

            mode = cu_ptr->pred_mode;

            // Hsan: if CHROMA_MODE_2, then CFL will be evaluated @ EP as no CHROMA @ MD
            // If that's the case then you should ensure than the 1st chroma prediction uses UV_DC_PRED (that's the default configuration for CHROMA_MODE_2 if CFL applicable (set @ fast loop candidates injection) then MD assumes chroma mode always UV_DC_PRED)
            av1_predict_intra_block(
                &sb_ptr->tile_info,
                ED_STAGE,
                context_ptr->blk_geom,
                picture_control_set_ptr->parent_pcs_ptr->av1_cm,
                context_ptr->blk_geom->bwidth,
                context_ptr->blk_geom->bheight,
                tx_size,
                mode,
                pu_ptr->angle_delta[PLANE_TYPE_Y],
                0,
                FILTER_INTRA_MODES,
                topNeighArray + 1,
                leftNeighArray + 1,
                recon_buffer,
                context_ptr->blk_geom->tx_boff_x[cu_ptr->tx_depth][context_ptr->txb_itr] >> 2,
                context_ptr->blk_geom->tx_boff_y[cu_ptr->tx_depth][context_ptr->txb_itr] >> 2,
                0,
                context_ptr->blk_geom->bsize,
                txb_origin_x,
                txb_origin_y,
                context_ptr->cu_origin_x,
                context_ptr->cu_origin_y,
                0,
                0);
        }
        // Encode Transform Unit -INTRA-

        uint8_t cb_qp = cu_ptr->qp;
        Av1EncodeLoopFunctionTable[is16bit](
            picture_control_set_ptr,
            context_ptr,
            sb_ptr,
            txb_origin_x,
            txb_origin_y,
            cb_qp,
            recon_buffer,
            coeff_buffer_sb,
            residual_buffer,
            transform_buffer,
            inverse_quant_buffer,
            transform_inner_array_ptr,
            asm_type,
            count_non_zero_coeffs,
            PICTURE_BUFFER_DESC_LUMA_MASK,
            0,
            cu_ptr->delta_qp > 0 ? 0 : dZoffset,
            eobs[context_ptr->txb_itr],
            cuPlane);

#if  CABAC_UP
        if (picture_control_set_ptr->update_cdf)
        {
            ModeDecisionCandidateBuffer         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
            ModeDecisionCandidateBuffer         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[0]);
            ModeDecisionCandidateBuffer          *candidateBuffer;

            // Set the Candidate Buffer
            candidateBuffer = candidate_buffer_ptr_array[0];
            // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
#if ATB_TX_TYPE_SUPPORT_PER_TU
            candidateBuffer->candidate_ptr->transform_type[context_ptr->txb_itr] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
            candidateBuffer->candidate_ptr->transform_type_uv = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
#else
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
#endif
            candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
            candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

            const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

            av1_tu_estimate_coeff_bits(
                1,//allow_update_cdf,
                &picture_control_set_ptr->ec_ctx_array[tbAddr],
                picture_control_set_ptr,
#if ATB_DC_CONTEXT_SUPPORT_0
                context_ptr->txb_itr,
#endif
                candidateBuffer,
                cu_ptr,
                coeff1dOffset,
                context_ptr->coded_area_sb_uv,
                coeff_est_entropy_coder_ptr,
                coeff_buffer_sb,
                eobs[context_ptr->txb_itr][0],
                eobs[context_ptr->txb_itr][1],
                eobs[context_ptr->txb_itr][2],
                &y_tu_coeff_bits,
                &cb_tu_coeff_bits,
                &cr_tu_coeff_bits,
                context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#if ATB_TX_TYPE_SUPPORT_PER_TU
                candidateBuffer->candidate_ptr->transform_type[context_ptr->txb_itr],
                candidateBuffer->candidate_ptr->transform_type_uv,
#endif
                COMPONENT_LUMA,
                asm_type);
        }
#endif

        Av1EncodeGenerateReconFunctionPtr[is16bit](
            context_ptr,
            txb_origin_x,
            txb_origin_y,
            recon_buffer,
            inverse_quant_buffer,
            transform_inner_array_ptr,
            PICTURE_BUFFER_DESC_LUMA_MASK,
            eobs[context_ptr->txb_itr],
            asm_type);

        // Update Recon Samples-INTRA-
        EncodePassUpdateReconSampleNeighborArrays(
            ep_luma_recon_neighbor_array,
            ep_cb_recon_neighbor_array,
            ep_cr_recon_neighbor_array,
            recon_buffer,
            txb_origin_x,
            txb_origin_y,
            context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            PICTURE_BUFFER_DESC_LUMA_MASK,
            is16bit);

        context_ptr->coded_area_sb += context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr] * context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr];
    } // Transform Loop

    // Chroma path

    if(context_ptr->blk_geom->has_uv)
    {
        context_ptr->txb_itr = 0;
        uint16_t txb_origin_x = context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[cu_ptr->tx_depth][context_ptr->txb_itr];
        uint16_t txb_origin_y = context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[cu_ptr->tx_depth][context_ptr->txb_itr];

        uint32_t cu_originx_uv = (context_ptr->cu_origin_x >> 3 << 3) >> 1;
        uint32_t cu_originy_uv = (context_ptr->cu_origin_y >> 3 << 3) >> 1;

#if DC_SIGN_CONTEXT_EP
        cu_ptr->cb_txb_skip_context = 0;
        cu_ptr->cb_dc_sign_context = 0;
        get_txb_ctx(
            COMPONENT_CHROMA,
            picture_control_set_ptr->ep_cb_dc_sign_level_coeff_neighbor_array,
            cu_originx_uv,
            cu_originy_uv,
            context_ptr->blk_geom->bsize_uv,
            context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            &cu_ptr->cb_txb_skip_context,
            &cu_ptr->cb_dc_sign_context);

        cu_ptr->cr_txb_skip_context = 0;
        cu_ptr->cr_dc_sign_context = 0;
        get_txb_ctx(
            COMPONENT_CHROMA,
            picture_control_set_ptr->ep_cr_dc_sign_level_coeff_neighbor_array,
            cu_originx_uv,
            cu_originy_uv,
            context_ptr->blk_geom->bsize_uv,
            context_ptr->blk_geom->txsize_uv[context_ptr->cu_ptr->tx_depth][context_ptr->txb_itr],
            &cu_ptr->cr_txb_skip_context,
            &cu_ptr->cr_dc_sign_context);
#endif

        if (is16bit) {
            uint16_t    topNeighArray[64 * 2 + 1];
            uint16_t    leftNeighArray[64 * 2 + 1];
            PredictionMode mode;

            int32_t plane_end = 2;

            for (int32_t plane = 1; plane <= plane_end; ++plane) {
                TxSize  tx_size = plane ? context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr] : context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr];

                if (plane == 1) {
                    if (cu_originy_uv != 0)
                        memcpy(topNeighArray + 1, (uint16_t*)(ep_cb_recon_neighbor_array->top_array) + cu_originx_uv, context_ptr->blk_geom->bwidth_uv * 2 * sizeof(uint16_t));
                    if (cu_originx_uv != 0)
                        memcpy(leftNeighArray + 1, (uint16_t*)(ep_cb_recon_neighbor_array->left_array) + cu_originy_uv, context_ptr->blk_geom->bheight_uv * 2 * sizeof(uint16_t));
                    if (cu_originy_uv != 0 && cu_originx_uv != 0)
                        topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_cb_recon_neighbor_array->top_left_array) + MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv)[0];
                }
                else if (plane == 2) {
                    if (cu_originy_uv != 0)
                        memcpy(topNeighArray + 1, (uint16_t*)(ep_cr_recon_neighbor_array->top_array) + cu_originx_uv, context_ptr->blk_geom->bwidth_uv * 2 * sizeof(uint16_t));
                    if (cu_originx_uv != 0)
                        memcpy(leftNeighArray + 1, (uint16_t*)(ep_cr_recon_neighbor_array->left_array) + cu_originy_uv, context_ptr->blk_geom->bheight_uv * 2 * sizeof(uint16_t));
                    if (cu_originy_uv != 0 && cu_originx_uv != 0)
                        topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_cr_recon_neighbor_array->top_left_array) + MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv)[0];
                }

                mode = (pu_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode)UV_DC_PRED : (PredictionMode)pu_ptr->intra_chroma_mode;

                av1_predict_intra_block_16bit(
                    &sb_ptr->tile_info,
                    context_ptr,
                    picture_control_set_ptr->parent_pcs_ptr->av1_cm,
                    plane ? context_ptr->blk_geom->bwidth_uv : context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr],
                    plane ? context_ptr->blk_geom->bheight_uv : context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr],
                    tx_size,
                    mode,
#if SEARCH_UV_MODE // conformance
                    plane ? pu_ptr->angle_delta[PLANE_TYPE_UV] : pu_ptr->angle_delta[PLANE_TYPE_Y],
#else
                    plane ? 0 : pu_ptr->angle_delta[PLANE_TYPE_Y],
#endif
                    0,
                    FILTER_INTRA_MODES,
                    topNeighArray + 1,
                    leftNeighArray + 1,
                    recon_buffer,
                    //int32_t dst_stride,
                    0,
                    0,
                    plane,
                    context_ptr->blk_geom->bsize,
                    plane ? context_ptr->cu_origin_x : context_ptr->cu_origin_x,
                    plane ? context_ptr->cu_origin_y : context_ptr->cu_origin_y);
            }
        }
        else {
            uint8_t    topNeighArray[64 * 2 + 1];
            uint8_t    leftNeighArray[64 * 2 + 1];
            PredictionMode mode;

            // Partition Loop
            int32_t plane_end = 2;

            for (int32_t plane = 1; plane <= plane_end; ++plane) {
                TxSize  tx_size = plane ? context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr] : context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr];

                if (plane == 1) {
                    if (cu_originy_uv != 0)
                        memcpy(topNeighArray + 1, ep_cb_recon_neighbor_array->top_array + cu_originx_uv, context_ptr->blk_geom->bwidth_uv * 2);

                    if (cu_originx_uv != 0)
                        memcpy(leftNeighArray + 1, ep_cb_recon_neighbor_array->left_array + cu_originy_uv, context_ptr->blk_geom->bheight_uv * 2);

                    if (cu_originy_uv != 0 && cu_originx_uv != 0)
                        topNeighArray[0] = leftNeighArray[0] = ep_cb_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv];
                }
                else {
                    if (cu_originy_uv != 0)
                        memcpy(topNeighArray + 1, ep_cr_recon_neighbor_array->top_array + cu_originx_uv, context_ptr->blk_geom->bwidth_uv * 2);

                    if (cu_originx_uv != 0)
                        memcpy(leftNeighArray + 1, ep_cr_recon_neighbor_array->left_array + cu_originy_uv, context_ptr->blk_geom->bheight_uv * 2);

                    if (cu_originy_uv != 0 && cu_originx_uv != 0)
                        topNeighArray[0] = leftNeighArray[0] = ep_cr_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv];
                }

                mode = (pu_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode)UV_DC_PRED : (PredictionMode)pu_ptr->intra_chroma_mode;

                // Hsan: if CHROMA_MODE_2, then CFL will be evaluated @ EP as no CHROMA @ MD
                // If that's the case then you should ensure than the 1st chroma prediction uses UV_DC_PRED (that's the default configuration for CHROMA_MODE_2 if CFL applicable (set @ fast loop candidates injection) then MD assumes chroma mode always UV_DC_PRED)
                av1_predict_intra_block(
                    &sb_ptr->tile_info,
                    ED_STAGE,
                    context_ptr->blk_geom,
                    picture_control_set_ptr->parent_pcs_ptr->av1_cm,
                    plane ? context_ptr->blk_geom->bwidth_uv : context_ptr->blk_geom->bwidth,
                    plane ? context_ptr->blk_geom->bheight_uv : context_ptr->blk_geom->bheight,
                    tx_size,
                    mode,
#if SEARCH_UV_MODE // conformance
                    plane ? pu_ptr->angle_delta[PLANE_TYPE_UV] : pu_ptr->angle_delta[PLANE_TYPE_Y],
#else
                    plane ? 0 : pu_ptr->angle_delta[PLANE_TYPE_Y],
#endif
                    0,
                    FILTER_INTRA_MODES,
                    topNeighArray + 1,
                    leftNeighArray + 1,
                    recon_buffer,
                    plane ? 0 : context_ptr->blk_geom->tx_boff_x[cu_ptr->tx_depth][context_ptr->txb_itr] >> 2,
                    plane ? 0 : context_ptr->blk_geom->tx_boff_y[cu_ptr->tx_depth][context_ptr->txb_itr] >> 2,
                    plane,
                    context_ptr->blk_geom->bsize,
                    txb_origin_x,
                    txb_origin_y,
                    plane ? context_ptr->cu_origin_x : context_ptr->cu_origin_x,
                    plane ? context_ptr->cu_origin_y : context_ptr->cu_origin_y,
                    0,
                    0);
            }
        }

        // Encode Transform Unit -INTRA-
        uint8_t cb_qp = cu_ptr->qp;

        Av1EncodeLoopFunctionTable[is16bit](
            picture_control_set_ptr,
            context_ptr,
            sb_ptr,
            txb_origin_x,
            txb_origin_y,
            cb_qp,
            recon_buffer,
            coeff_buffer_sb,
            residual_buffer,
            transform_buffer,
            inverse_quant_buffer,
            transform_inner_array_ptr,
            asm_type,
            count_non_zero_coeffs,
            PICTURE_BUFFER_DESC_CHROMA_MASK,
            0,
            cu_ptr->delta_qp > 0 ? 0 : dZoffset,
            eobs[context_ptr->txb_itr],
            cuPlane);

#if  CABAC_UP
        if (picture_control_set_ptr->update_cdf)
        {
            ModeDecisionCandidateBuffer         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
            ModeDecisionCandidateBuffer         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[0]);
            ModeDecisionCandidateBuffer          *candidateBuffer;

            // Set the Candidate Buffer
            candidateBuffer = candidate_buffer_ptr_array[0];
            // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
#if ATB_TX_TYPE_SUPPORT_PER_TU
            candidateBuffer->candidate_ptr->transform_type[context_ptr->txb_itr] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
            candidateBuffer->candidate_ptr->transform_type_uv = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
#else
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
            candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
#endif
            candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
            candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

            const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

            av1_tu_estimate_coeff_bits(
                1,//allow_update_cdf,
                &picture_control_set_ptr->ec_ctx_array[tbAddr],
                picture_control_set_ptr,
#if ATB_DC_CONTEXT_SUPPORT_0
                context_ptr->txb_itr,
#endif
                candidateBuffer,
                cu_ptr,
                coeff1dOffset,
                context_ptr->coded_area_sb_uv,
                coeff_est_entropy_coder_ptr,
                coeff_buffer_sb,
                eobs[context_ptr->txb_itr][0],
                eobs[context_ptr->txb_itr][1],
                eobs[context_ptr->txb_itr][2],
                &y_tu_coeff_bits,
                &cb_tu_coeff_bits,
                &cr_tu_coeff_bits,
                context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#if ATB_TX_TYPE_SUPPORT_PER_TU
                candidateBuffer->candidate_ptr->transform_type[context_ptr->txb_itr],
                candidateBuffer->candidate_ptr->transform_type_uv,
#endif
                COMPONENT_CHROMA,
                asm_type);
        }
#endif

        Av1EncodeGenerateReconFunctionPtr[is16bit](
            context_ptr,
            txb_origin_x,
            txb_origin_y,
            recon_buffer,
            inverse_quant_buffer,
            transform_inner_array_ptr,
            PICTURE_BUFFER_DESC_CHROMA_MASK,
            eobs[context_ptr->txb_itr],
            asm_type);

        // Update Recon Samples-INTRA-
        EncodePassUpdateReconSampleNeighborArrays(
            ep_luma_recon_neighbor_array,
            ep_cb_recon_neighbor_array,
            ep_cr_recon_neighbor_array,
            recon_buffer,
            txb_origin_x,
            txb_origin_y,
            context_ptr->blk_geom->tx_width[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
            PICTURE_BUFFER_DESC_CHROMA_MASK,
            is16bit);

        context_ptr->coded_area_sb_uv += context_ptr->blk_geom->tx_width_uv[cu_ptr->tx_depth][context_ptr->txb_itr] * context_ptr->blk_geom->tx_height_uv[cu_ptr->tx_depth][context_ptr->txb_itr];
    } // Transform Loop

    for (context_ptr->txb_itr = 0; context_ptr->txb_itr < totTu; context_ptr->txb_itr++) {
        uint8_t uv_pass = cu_ptr->tx_depth && context_ptr->txb_itr ? 0 : 1;

        if (context_ptr->blk_geom->has_uv && uv_pass) {
            cu_ptr->block_has_coeff = cu_ptr->block_has_coeff |
                cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff |
                cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff |
                cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff;

            if (cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff)
                cu_ptr->transform_unit_array[0].u_has_coeff = EB_TRUE;
            if (cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff)
                cu_ptr->transform_unit_array[0].v_has_coeff = EB_TRUE;
        }
        else {
            cu_ptr->block_has_coeff = cu_ptr->block_has_coeff |
                cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff;
        }
    } // Transform Loop
}

#endif
/*******************************************
* Encode Pass
*
* Summary: Performs a H.265 conformant
*   reconstruction based on the LCU
*   mode decision.
*
* Inputs:
*   SourcePic
*   Coding Results
*   SB Location
*   Sequence Control Set
*   Picture Control Set
*
* Outputs:
*   Reconstructed Samples
*   Coefficient Samples
*
*******************************************/
EB_EXTERN void av1_encode_pass(
    SequenceControlSet      *sequence_control_set_ptr,
    PictureControlSet       *picture_control_set_ptr,
    LargestCodingUnit       *sb_ptr,
    uint32_t                 tbAddr,
    uint32_t                 sb_origin_x,
    uint32_t                 sb_origin_y,
#if !MEMORY_FOOTPRINT_OPT
    uint32_t                 sb_qp,
#endif
    EncDecContext           *context_ptr)
{
    EbBool                    is16bit = context_ptr->is16bit;
    EbPictureBufferDesc    *recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;
    EbPictureBufferDesc    *coeff_buffer_sb = sb_ptr->quantized_coeff;
    EbPictureBufferDesc    *inputPicture;
    ModeDecisionContext    *mdcontextPtr;
    mdcontextPtr = context_ptr->md_context;
    inputPicture = context_ptr->input_samples = (EbPictureBufferDesc*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;

    SbStat                *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[tbAddr]);
#if !MEMORY_FOOTPRINT_OPT
    EbBool                    availableCoeff;
    // QP Neighbor Arrays
    EbBool                    isDeltaQpNotCoded = EB_TRUE;
#endif
    // SB Stats
    uint32_t                  sb_width = MIN(sequence_control_set_ptr->sb_size_pix, sequence_control_set_ptr->luma_width - sb_origin_x);
    uint32_t                  sb_height = MIN(sequence_control_set_ptr->sb_size_pix, sequence_control_set_ptr->luma_height - sb_origin_y);
    // MV merge mode
    uint32_t                  y_has_coeff;
    uint32_t                  u_has_coeff;
    uint32_t                  v_has_coeff;
    EbAsm                     asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
    uint64_t                  y_coeff_bits;
    uint64_t                  cb_coeff_bits;
    uint64_t                  cr_coeff_bits;
    uint64_t                  y_full_distortion[DIST_CALC_TOTAL];
    uint64_t                  yTuFullDistortion[DIST_CALC_TOTAL];
    uint32_t                  count_non_zero_coeffs[3];
    MacroblockPlane           cuPlane[3];
    uint16_t                  eobs[MAX_TXB_COUNT][3];
    uint64_t                  y_tu_coeff_bits;
    uint64_t                  cb_tu_coeff_bits;
    uint64_t                  cr_tu_coeff_bits;
    EncodeContext          *encode_context_ptr;
#if !MEMORY_FOOTPRINT_OPT
    uint32_t                  lcuRowIndex = sb_origin_y / BLOCK_SIZE_64;
#endif
    // Dereferencing early
    NeighborArrayUnit      *ep_mode_type_neighbor_array = picture_control_set_ptr->ep_mode_type_neighbor_array;
    NeighborArrayUnit      *ep_intra_luma_mode_neighbor_array = picture_control_set_ptr->ep_intra_luma_mode_neighbor_array;
    NeighborArrayUnit      *ep_intra_chroma_mode_neighbor_array = picture_control_set_ptr->ep_intra_chroma_mode_neighbor_array;
    NeighborArrayUnit      *ep_mv_neighbor_array = picture_control_set_ptr->ep_mv_neighbor_array;
    NeighborArrayUnit      *ep_luma_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_luma_recon_neighbor_array16bit : picture_control_set_ptr->ep_luma_recon_neighbor_array;
    NeighborArrayUnit      *ep_cb_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_cb_recon_neighbor_array16bit : picture_control_set_ptr->ep_cb_recon_neighbor_array;
    NeighborArrayUnit      *ep_cr_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_cr_recon_neighbor_array16bit : picture_control_set_ptr->ep_cr_recon_neighbor_array;
    NeighborArrayUnit      *ep_skip_flag_neighbor_array = picture_control_set_ptr->ep_skip_flag_neighbor_array;
#if DC_SIGN_CONTEXT_EP
    NeighborArrayUnit      *ep_luma_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->ep_luma_dc_sign_level_coeff_neighbor_array;
    NeighborArrayUnit      *ep_cb_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->ep_cb_dc_sign_level_coeff_neighbor_array;
    NeighborArrayUnit      *ep_cr_dc_sign_level_coeff_neighbor_array = picture_control_set_ptr->ep_cr_dc_sign_level_coeff_neighbor_array;
#endif

    EbBool                 constrained_intra_flag = picture_control_set_ptr->constrained_intra_flag;

    EbBool dlfEnableFlag = (EbBool)(picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode &&
        (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag ||
            sequence_control_set_ptr->static_config.recon_enabled ||
            sequence_control_set_ptr->static_config.stat_report));

    const EbBool isIntraLCU = picture_control_set_ptr->limit_intra ? EB_FALSE : EB_TRUE;

    EbBool doRecon = (EbBool)(
        (picture_control_set_ptr->limit_intra == 0 || isIntraLCU == 1) ||
        picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag ||
        sequence_control_set_ptr->static_config.recon_enabled ||
        sequence_control_set_ptr->static_config.stat_report);

    EntropyCoder  *coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;

    uint32_t           dZoffset = 0;
#if !MEMORY_FOOTPRINT_OPT
    if (!sb_stat_ptr->stationary_edge_over_time_flag && sequence_control_set_ptr->static_config.improve_sharpness && picture_control_set_ptr->parent_pcs_ptr->pic_noise_class < PIC_NOISE_CLASS_3_1) {
        int16_t cuDeltaQp = (int16_t)(sb_ptr->qp - picture_control_set_ptr->parent_pcs_ptr->average_qp);
        uint32_t dzCondition = cuDeltaQp > 0 ? 0 : 1;

        if (sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE) {
            if (!(picture_control_set_ptr->parent_pcs_ptr->is_pan ||
                (picture_control_set_ptr->parent_pcs_ptr->non_moving_index_average < 10 && sb_ptr->aura_status_iii) ||
                (sb_stat_ptr->cu_stat_array[0].skin_area) ||
#if !DISABLE_OIS_USE
                (picture_control_set_ptr->parent_pcs_ptr->intra_coded_block_probability > 90) ||
#endif
                (picture_control_set_ptr->parent_pcs_ptr->high_dark_area_density_flag))) {
                if (picture_control_set_ptr->slice_type != I_SLICE &&
                    picture_control_set_ptr->temporal_layer_index == 0 &&
#if !DISABLE_OIS_USE
                    picture_control_set_ptr->parent_pcs_ptr->intra_coded_block_probability > 60 &&
#endif
                    !picture_control_set_ptr->parent_pcs_ptr->is_tilt &&
                    picture_control_set_ptr->parent_pcs_ptr->pic_homogenous_over_time_sb_percentage > 40)
                {
                    dZoffset = 10;
                }

                if (dzCondition) {
                    if (picture_control_set_ptr->scene_caracteristic_id == EB_FRAME_CARAC_1) {
                        if (picture_control_set_ptr->slice_type == I_SLICE)
                            dZoffset = sb_stat_ptr->cu_stat_array[0].grass_area ? 10 : dZoffset;
                        else if (picture_control_set_ptr->temporal_layer_index == 0)
                            dZoffset = sb_stat_ptr->cu_stat_array[0].grass_area ? 9 : dZoffset;
                        else if (picture_control_set_ptr->temporal_layer_index == 1)
                            dZoffset = sb_stat_ptr->cu_stat_array[0].grass_area ? 5 : dZoffset;
                    }
                }
            }
        }
    }
#endif
#if MEMORY_FOOTPRINT_OPT
    context_ptr->skip_qpm_flag = EB_TRUE;
#else
    if (sequence_control_set_ptr->static_config.improve_sharpness) {
        QpmDeriveBeaAndSkipQpmFlagLcu(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_ptr,
            tbAddr,
            context_ptr);
    }
    else
        context_ptr->skip_qpm_flag = EB_TRUE;
#endif

    encode_context_ptr = ((SequenceControlSet*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr))->encode_context_ptr;

    if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
        //get the 16bit form of the input LCU
        if (is16bit)
            recon_buffer = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            recon_buffer = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
    else  // non ref pictures
        recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;

    EbBool use_delta_qp = (EbBool)sequence_control_set_ptr->static_config.improve_sharpness;
    EbBool oneSegment = (sequence_control_set_ptr->enc_dec_segment_col_count_array[picture_control_set_ptr->temporal_layer_index] == 1) && (sequence_control_set_ptr->enc_dec_segment_row_count_array[picture_control_set_ptr->temporal_layer_index] == 1);
    EbBool useDeltaQpSegments = oneSegment ? 0 : (EbBool)sequence_control_set_ptr->static_config.improve_sharpness;

    // DeriveZeroLumaCbf
    EbBool  highIntraRef = EB_FALSE;
    EbBool  checkZeroLumaCbf = EB_FALSE;

    if (is16bit) {
        //SB128_TODO change 10bit SB creation

        if ((sequence_control_set_ptr->static_config.ten_bit_format == 1) || (sequence_control_set_ptr->static_config.compressed_ten_bit_format == 1))
        {
            const uint32_t input_luma_offset = ((sb_origin_y + inputPicture->origin_y)         * inputPicture->stride_y) + (sb_origin_x + inputPicture->origin_x);
            const uint32_t input_cb_offset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->stride_cb) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint32_t input_cr_offset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->stride_cr) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint16_t luma2BitWidth = inputPicture->width / 4;
            const uint16_t chroma2BitWidth = inputPicture->width / 8;

            compressed_pack_lcu(
                inputPicture->buffer_y + input_luma_offset,
                inputPicture->stride_y,
                inputPicture->buffer_bit_inc_y + sb_origin_y * luma2BitWidth + (sb_origin_x / 4)*sb_height,
                sb_width / 4,
                (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_y,
                context_ptr->input_sample16bit_buffer->stride_y,
                sb_width,
                sb_height,
                asm_type);

            compressed_pack_lcu(
                inputPicture->buffer_cb + input_cb_offset,
                inputPicture->stride_cb,
                inputPicture->buffer_bit_inc_cb + sb_origin_y / 2 * chroma2BitWidth + (sb_origin_x / 8)*(sb_height / 2),
                sb_width / 8,
                (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_cb,
                context_ptr->input_sample16bit_buffer->stride_cb,
                sb_width >> 1,
                sb_height >> 1,
                asm_type);

            compressed_pack_lcu(
                inputPicture->buffer_cr + input_cr_offset,
                inputPicture->stride_cr,
                inputPicture->buffer_bit_inc_cr + sb_origin_y / 2 * chroma2BitWidth + (sb_origin_x / 8)*(sb_height / 2),
                sb_width / 8,
                (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_cr,
                context_ptr->input_sample16bit_buffer->stride_cr,
                sb_width >> 1,
                sb_height >> 1,
                asm_type);
        }
        else {
            const uint32_t input_luma_offset = ((sb_origin_y + inputPicture->origin_y)         * inputPicture->stride_y) + (sb_origin_x + inputPicture->origin_x);
            const uint32_t inputBitIncLumaOffset = ((sb_origin_y + inputPicture->origin_y)         * inputPicture->stride_bit_inc_y) + (sb_origin_x + inputPicture->origin_x);
            const uint32_t input_cb_offset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->stride_cb) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint32_t inputBitIncCbOffset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->stride_bit_inc_cb) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint32_t input_cr_offset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->stride_cr) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint32_t inputBitIncCrOffset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->stride_bit_inc_cr) + ((sb_origin_x + inputPicture->origin_x) >> 1);

            pack2d_src(
                inputPicture->buffer_y + input_luma_offset,
                inputPicture->stride_y,
                inputPicture->buffer_bit_inc_y + inputBitIncLumaOffset,
                inputPicture->stride_bit_inc_y,
                (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_y,
                context_ptr->input_sample16bit_buffer->stride_y,
                sb_width,
                sb_height,
                asm_type);

            pack2d_src(
                inputPicture->buffer_cb + input_cb_offset,
                inputPicture->stride_cr,
                inputPicture->buffer_bit_inc_cb + inputBitIncCbOffset,
                inputPicture->stride_bit_inc_cr,
                (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_cb,
                context_ptr->input_sample16bit_buffer->stride_cb,
                sb_width >> 1,
                sb_height >> 1,
                asm_type);

            pack2d_src(
                inputPicture->buffer_cr + input_cr_offset,
                inputPicture->stride_cr,
                inputPicture->buffer_bit_inc_cr + inputBitIncCrOffset,
                inputPicture->stride_bit_inc_cr,
                (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_cr,
                context_ptr->input_sample16bit_buffer->stride_cr,
                sb_width >> 1,
                sb_height >> 1,
                asm_type);
        }

        Store16bitInputSrc(context_ptr, picture_control_set_ptr, sb_origin_x, sb_origin_y, sb_width, sb_height);
    }

    if ((sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE) && !picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag) {
        if (!((sb_stat_ptr->stationary_edge_over_time_flag) || (picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag)))
        {
            if (picture_control_set_ptr->slice_type == B_SLICE) {
#if MRP_MD
                EbReferenceObject  *refObjL0 = (EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                EbReferenceObject  *refObjL1 = (EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
#else
                EbReferenceObject  *refObjL0 = (EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;
                EbReferenceObject  *refObjL1 = (EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr;
#endif
                uint32_t const TH = (sequence_control_set_ptr->static_config.frame_rate >> 16) < 50 ? 25 : 30;

                if ((refObjL0->tmp_layer_idx == 2 && refObjL0->intra_coded_area > TH) || (refObjL1->tmp_layer_idx == 2 && refObjL1->intra_coded_area > TH))
                    highIntraRef = EB_TRUE;
            }
            if (highIntraRef == EB_FALSE)
                checkZeroLumaCbf = EB_TRUE;
        }
    }
    context_ptr->intra_coded_area_sb[tbAddr] = 0;
#if !PF_N2_SUPPORT
    context_ptr->trans_coeff_shape_luma = 0;
    context_ptr->trans_coeff_shape_chroma = 0;
#endif
    context_ptr->coded_area_sb = 0;
    context_ptr->coded_area_sb_uv = 0;

#if AV1_LF
    if (dlfEnableFlag && picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode == 1){
        if (tbAddr == 0) {
            av1_loop_filter_init(picture_control_set_ptr);

            av1_pick_filter_level(
                0,
                (EbPictureBufferDesc*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                picture_control_set_ptr,
                LPF_PICK_FROM_Q);

            av1_loop_filter_frame_init(picture_control_set_ptr, 0, 3);
        }
    }
#endif
#if ADD_DELTA_QP_SUPPORT
    if (context_ptr->skip_qpm_flag == EB_FALSE && sequence_control_set_ptr->static_config.improve_sharpness) {
        Av1QpModulationLcu(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_ptr,
            tbAddr,
            picture_control_set_ptr->slice_type == I_SLICE ? INTRA_MODE : INTER_MODE,
            context_ptr);
    }
#endif

#if CABAC_UP
    uint8_t allow_update_cdf = picture_control_set_ptr->update_cdf;
#endif

    uint32_t final_cu_itr = 0;

    // CU Loop

    uint32_t    blk_it = 0;

    while (blk_it < sequence_control_set_ptr->max_block_cnt) {
        CodingUnit  *cu_ptr = context_ptr->cu_ptr = &context_ptr->md_context->md_cu_arr_nsq[blk_it];
        PartitionType part = cu_ptr->part;

        const BlockGeom * blk_geom = context_ptr->blk_geom = get_blk_geom_mds(blk_it);
        UNUSED(blk_geom);

        sb_ptr->cu_partition_array[blk_it] = context_ptr->md_context->md_cu_arr_nsq[blk_it].part;

        if (part != PARTITION_SPLIT) {
            int32_t offset_d1 = ns_blk_offset[(int32_t)part]; //cu_ptr->best_d1_blk; // TOCKECK
            int32_t num_d1_block = ns_blk_num[(int32_t)part]; // context_ptr->blk_geom->totns; // TOCKECK

           // for (int32_t d1_itr = blk_it; d1_itr < blk_it + num_d1_block; d1_itr++) {
            for (int32_t d1_itr = (int32_t)blk_it + offset_d1; d1_itr < (int32_t)blk_it + offset_d1 + num_d1_block; d1_itr++) {
                const BlockGeom * blk_geom = context_ptr->blk_geom = get_blk_geom_mds(d1_itr);

                // PU Stack variables
                PredictionUnit        *pu_ptr = (PredictionUnit *)EB_NULL; //  done
                EbPictureBufferDesc   *residual_buffer = context_ptr->residual_buffer;
                EbPictureBufferDesc   *transform_buffer = context_ptr->transform_buffer;

                EbPictureBufferDesc   *inverse_quant_buffer = context_ptr->inverse_quant_buffer;

                int16_t                  *transform_inner_array_ptr = context_ptr->transform_inner_array_ptr;

                CodingUnit            *cu_ptr = context_ptr->cu_ptr = &context_ptr->md_context->md_cu_arr_nsq[d1_itr];

                context_ptr->cu_origin_x = (uint16_t)(sb_origin_x + blk_geom->origin_x);
                context_ptr->cu_origin_y = (uint16_t)(sb_origin_y + blk_geom->origin_y);
                cu_ptr->delta_qp = 0;
#if  BLK_SKIP_DECISION
                context_ptr->md_skip_blk = context_ptr->md_context->blk_skip_decision ? ((cu_ptr->prediction_mode_flag == INTRA_MODE || cu_ptr->block_has_coeff) ? 0 : 1) : 0;
#endif
                cu_ptr->block_has_coeff = 0;

                // if(picture_control_set_ptr->picture_number==4 && context_ptr->cu_origin_x==0 && context_ptr->cu_origin_y==0)
                //     printf("CHEDD");
                uint32_t  coded_area_org = context_ptr->coded_area_sb;
                uint32_t  coded_area_org_uv = context_ptr->coded_area_sb_uv;

                // Derive disable_cfl_flag as evaluate_cfl_ep = f(disable_cfl_flag)
                EbBool disable_cfl_flag = (context_ptr->blk_geom->sq_size > 32 ||
                    context_ptr->blk_geom->bwidth == 4 ||
                    context_ptr->blk_geom->bheight == 4) ? EB_TRUE : EB_FALSE;
                // Evaluate cfl @ EP if applicable, and not done @ MD
                context_ptr->evaluate_cfl_ep = (disable_cfl_flag == EB_FALSE && context_ptr->md_context->chroma_level == CHROMA_MODE_2);

#if ADD_DELTA_QP_SUPPORT
                if (context_ptr->skip_qpm_flag == EB_FALSE && sequence_control_set_ptr->static_config.improve_sharpness) {
                    cu_ptr->qp = sb_ptr->qp;
                    cu_ptr->delta_qp = sb_ptr->delta_qp;
                    cu_ptr->org_delta_qp = sb_ptr->org_delta_qp;
                }
                else {
                    uint16_t                           picture_qp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
                    sb_ptr->qp = picture_qp;
                    cu_ptr->qp = sb_ptr->qp;
                    cu_ptr->delta_qp = 0;
                    cu_ptr->org_delta_qp = 0;
                }

#else
                cu_ptr->qp = (sequence_control_set_ptr->static_config.improve_sharpness) ? context_ptr->qpm_qp : picture_control_set_ptr->picture_qp;
                sb_ptr->qp = (sequence_control_set_ptr->static_config.improve_sharpness) ? context_ptr->qpm_qp : picture_control_set_ptr->picture_qp;
#if !MEMORY_FOOTPRINT_OPT
                cu_ptr->org_delta_qp = cu_ptr->delta_qp;
#endif
#endif
#if !MEMORY_FOOTPRINT_OPT
#if !ADD_DELTA_QP_SUPPORT
                //CHKN remove usage of depth
                if (!context_ptr->skip_qpm_flag && (sequence_control_set_ptr->static_config.improve_sharpness) && (0xFF <= picture_control_set_ptr->dif_cu_delta_qp_depth)) {
                    EncQpmDeriveDeltaQPForEachLeafLcu(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        sb_ptr,
                        tbAddr,
                        cu_ptr,
                        0xFF, //This is obviously not ok
                        d1_itr, // TOCHECK OMK
                        context_ptr->blk_geom->bwidth, // TOCHECK
                        cu_ptr->prediction_mode_flag,
                        context_ptr->cu_stats->parent32x32_index, // TOCHECK not valid
                        context_ptr);
                }

#endif
#endif
                if (cu_ptr->prediction_mode_flag == INTRA_MODE) {
                    context_ptr->is_inter = cu_ptr->av1xd->use_intrabc;
                    context_ptr->tot_intra_coded_area += blk_geom->bwidth* blk_geom->bheight;
                    if (picture_control_set_ptr->slice_type != I_SLICE)
                        context_ptr->intra_coded_area_sb[tbAddr] += blk_geom->bwidth* blk_geom->bheight;
                    // *Note - Transforms are the same size as predictions
                    // Partition Loop
                    context_ptr->txb_itr = 0;
#if ATB_EP
                    // Transform partitioning path (INTRA Luma/Chroma)
                    if (cu_ptr->av1xd->use_intrabc == 0) {
                        // Set the PU Loop Variables
                        pu_ptr = cu_ptr->prediction_unit_array;
                        // Generate Intra Luma Neighbor Modes
                        GeneratePuIntraLumaNeighborModes(
                            cu_ptr,
                            context_ptr->cu_origin_x,
                            context_ptr->cu_origin_y,
                            BLOCK_SIZE_64,
                            ep_intra_luma_mode_neighbor_array,
                            ep_intra_chroma_mode_neighbor_array,
                            ep_mode_type_neighbor_array);

                        perform_intra_coding_loop(
                            sequence_control_set_ptr,
                            picture_control_set_ptr,
                            sb_ptr,
                            tbAddr,
                            cu_ptr,
                            pu_ptr,
                            context_ptr,
                            dZoffset);

                        // Update the Intra-specific Neighbor Arrays
                        EncodePassUpdateIntraModeNeighborArrays(
#if DC_SIGN_CONTEXT_EP
                            context_ptr,
                            ep_luma_dc_sign_level_coeff_neighbor_array,
                            ep_cb_dc_sign_level_coeff_neighbor_array,
                            ep_cr_dc_sign_level_coeff_neighbor_array,
#endif
                            ep_mode_type_neighbor_array,
                            ep_intra_luma_mode_neighbor_array,
                            ep_intra_chroma_mode_neighbor_array,
                            (uint8_t)cu_ptr->pred_mode,
                            (uint8_t)pu_ptr->intra_chroma_mode,
                            context_ptr->cu_origin_x,
                            context_ptr->cu_origin_y,
                            context_ptr->blk_geom->bwidth,
                            context_ptr->blk_geom->bheight,
                            context_ptr->blk_geom->bwidth_uv,
                            context_ptr->blk_geom->bheight_uv,
                            blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK);
                    }
                    // Transform partitioning free patch (except the 128x128 case)
                    else
#endif
                    {
                        // Set the PU Loop Variables
                        pu_ptr = cu_ptr->prediction_unit_array;
                        // Generate Intra Luma Neighbor Modes
                        GeneratePuIntraLumaNeighborModes( // HT done
                            cu_ptr,
                            context_ptr->cu_origin_x,
                            context_ptr->cu_origin_y,
                            BLOCK_SIZE_64,
                            ep_intra_luma_mode_neighbor_array,
                            ep_intra_chroma_mode_neighbor_array,
                            ep_mode_type_neighbor_array);

                        {
                           uint32_t cu_originy_uv = (context_ptr->cu_origin_y >> 3 << 3) >> 1;
                           uint32_t cu_originx_uv = (context_ptr->cu_origin_x >> 3 << 3) >> 1;

#if DC_SIGN_CONTEXT_EP
                            context_ptr->cu_ptr->luma_txb_skip_context = 0;
                            context_ptr->cu_ptr->luma_dc_sign_context[context_ptr->txb_itr] = 0;
                            get_txb_ctx(
                                COMPONENT_LUMA,
                                picture_control_set_ptr->ep_luma_dc_sign_level_coeff_neighbor_array,
                                context_ptr->cu_origin_x,
                                context_ptr->cu_origin_y,
                                context_ptr->blk_geom->bsize,
                                context_ptr->blk_geom->txsize[0][0],
                                &context_ptr->cu_ptr->luma_txb_skip_context,
                                &context_ptr->cu_ptr->luma_dc_sign_context[0]);

                            cu_ptr->cb_txb_skip_context = 0;
                            cu_ptr->cb_dc_sign_context = 0;
                            get_txb_ctx(
                                COMPONENT_CHROMA,
                                picture_control_set_ptr->ep_cb_dc_sign_level_coeff_neighbor_array,
                                cu_originx_uv,
                                cu_originy_uv,
                                context_ptr->blk_geom->bsize_uv,
                                context_ptr->blk_geom->txsize_uv[0][0],
                                &cu_ptr->cb_txb_skip_context,
                                &cu_ptr->cb_dc_sign_context);

                            cu_ptr->cr_txb_skip_context = 0;
                            cu_ptr->cr_dc_sign_context = 0;
                            get_txb_ctx(
                                COMPONENT_CHROMA,
                                picture_control_set_ptr->ep_cr_dc_sign_level_coeff_neighbor_array,
                                cu_originx_uv,
                                cu_originy_uv,
                                context_ptr->blk_geom->bsize_uv,
                                context_ptr->blk_geom->txsize_uv[0][0],
                                &cu_ptr->cr_txb_skip_context,
                                &cu_ptr->cr_dc_sign_context);
#endif

                            if (cu_ptr->av1xd->use_intrabc)
                            {
                                MvReferenceFrame ref_frame = INTRA_FRAME;
                                generate_av1_mvp_table(
                                    &sb_ptr->tile_info,
                                    context_ptr->md_context,
                                    cu_ptr,
                                    context_ptr->blk_geom,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    &ref_frame,
                                    1,
                                    picture_control_set_ptr);

                                IntMv nearestmv, nearmv;
                                av1_find_best_ref_mvs_from_stack(0, context_ptr->md_context->md_local_cu_unit[blk_geom->blkidx_mds].ed_ref_mv_stack, cu_ptr->av1xd, ref_frame, &nearestmv, &nearmv,
                                    0);

                                if (nearestmv.as_int == INVALID_MV)
                                    nearestmv.as_int = 0;
                                if (nearmv.as_int == INVALID_MV)
                                    nearmv.as_int = 0;
                                IntMv dv_ref = nearestmv.as_int == 0 ? nearmv : nearestmv;
                                if (dv_ref.as_int == 0)
                                    av1_find_ref_dv(&dv_ref, &cu_ptr->av1xd->tile, sequence_control_set_ptr->mib_size, context_ptr->cu_origin_y >> MI_SIZE_LOG2, context_ptr->cu_origin_x >> MI_SIZE_LOG2);
                                // Ref DV should not have sub-pel.
                                assert((dv_ref.as_mv.col & 7) == 0);
                                assert((dv_ref.as_mv.row & 7) == 0);
                                context_ptr->md_context->md_local_cu_unit[blk_geom->blkidx_mds].ed_ref_mv_stack[INTRA_FRAME][0].this_mv = dv_ref;
                                cu_ptr->predmv[0] = dv_ref;

                                //keep final usefull mvp for entropy
                                memcpy(cu_ptr->av1xd->final_ref_mv_stack,
                                    context_ptr->md_context->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[cu_ptr->prediction_unit_array[0].ref_frame_type],
                                    sizeof(CandidateMv)*MAX_REF_MV_STACK_SIZE);

                                pu_ptr = cu_ptr->prediction_unit_array;
                                // Set MvUnit
                                context_ptr->mv_unit.pred_direction = (uint8_t)pu_ptr->inter_pred_direction_index;
                                context_ptr->mv_unit.mv[REF_LIST_0].mv_union = pu_ptr->mv[REF_LIST_0].mv_union;
                                context_ptr->mv_unit.mv[REF_LIST_1].mv_union = pu_ptr->mv[REF_LIST_1].mv_union;

                                EbPictureBufferDesc * ref_pic_list0 = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;

                                if (is16bit)
                                    ref_pic_list0 = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;

                                if (is16bit)
                                    av1_inter_prediction_hbd(
                                        picture_control_set_ptr,
                                        cu_ptr->prediction_unit_array->ref_frame_type,
                                        cu_ptr,
                                        &context_ptr->mv_unit,
                                        1,//intrabc
                                        context_ptr->cu_origin_x,
                                        context_ptr->cu_origin_y,
                                        blk_geom->bwidth,
                                        blk_geom->bheight,
                                        ref_pic_list0,
                                        0,
                                        recon_buffer,
                                        context_ptr->cu_origin_x,
                                        context_ptr->cu_origin_y,
                                        (uint8_t)sequence_control_set_ptr->static_config.encoder_bit_depth,
                                        asm_type);
                                else
                                av1_inter_prediction(
                                    picture_control_set_ptr,
                                    cu_ptr->interp_filters,
                                    cu_ptr,
                                    cu_ptr->prediction_unit_array->ref_frame_type,
                                    &context_ptr->mv_unit,
                                    1,// use_intrabc,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    blk_geom->bwidth,
                                    blk_geom->bheight,
                                    ref_pic_list0,
                                    0,
                                    recon_buffer,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    EB_TRUE,
                                    asm_type);
                            }
                            else
                            {
                                if (is16bit) {
                                    uint16_t    topNeighArray[64 * 2 + 1];
                                    uint16_t    leftNeighArray[64 * 2 + 1];
                                    PredictionMode mode;

                                int32_t plane_end = blk_geom->has_uv ? 2 : 0;

                                for (int32_t plane = 0; plane <= plane_end; ++plane) {
#if ATB_SUPPORT
                                    TxSize  tx_size = plane ? blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr] : blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr];
#else
                                    TxSize  tx_size = plane ? blk_geom->txsize_uv[context_ptr->txb_itr] : blk_geom->txsize[context_ptr->txb_itr];
#endif
                                    if (plane == 0) {
                                        if (context_ptr->cu_origin_y != 0)
                                            memcpy(topNeighArray + 1, (uint16_t*)(ep_luma_recon_neighbor_array->top_array) + context_ptr->cu_origin_x, blk_geom->bwidth * 2 * sizeof(uint16_t));
                                        if (context_ptr->cu_origin_x != 0)
                                            memcpy(leftNeighArray + 1, (uint16_t*)(ep_luma_recon_neighbor_array->left_array) + context_ptr->cu_origin_y, blk_geom->bheight * 2 * sizeof(uint16_t));
                                        if (context_ptr->cu_origin_y != 0 && context_ptr->cu_origin_x != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_luma_recon_neighbor_array->top_left_array) + MAX_PICTURE_HEIGHT_SIZE + context_ptr->cu_origin_x - context_ptr->cu_origin_y)[0];
                                    }

                                    else if (plane == 1) {
                                        if (cu_originy_uv != 0)
                                            memcpy(topNeighArray + 1, (uint16_t*)(ep_cb_recon_neighbor_array->top_array) + cu_originx_uv, blk_geom->bwidth_uv * 2 * sizeof(uint16_t));
                                        if (cu_originx_uv != 0)
                                            memcpy(leftNeighArray + 1, (uint16_t*)(ep_cb_recon_neighbor_array->left_array) + cu_originy_uv, blk_geom->bheight_uv * 2 * sizeof(uint16_t));
                                        if (cu_originy_uv != 0 && cu_originx_uv != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_cb_recon_neighbor_array->top_left_array) + MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv)[0];
                                    }
                                    else {
                                        if (cu_originy_uv != 0)
                                            memcpy(topNeighArray + 1, (uint16_t*)(ep_cr_recon_neighbor_array->top_array) + cu_originx_uv, blk_geom->bwidth_uv * 2 * sizeof(uint16_t));
                                        if (cu_originx_uv != 0)
                                            memcpy(leftNeighArray + 1, (uint16_t*)(ep_cr_recon_neighbor_array->left_array) + cu_originy_uv, blk_geom->bheight_uv * 2 * sizeof(uint16_t));
                                        if (cu_originy_uv != 0 && cu_originx_uv != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_cr_recon_neighbor_array->top_left_array) + MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv)[0];
                                    }

                                    if (plane)
                                        mode = (pu_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode)UV_DC_PRED : (PredictionMode)pu_ptr->intra_chroma_mode;
                                    else
                                        mode = cu_ptr->pred_mode; //PredictionMode mode,
                                    av1_predict_intra_block_16bit(
                                        &sb_ptr->tile_info,
                                        context_ptr,
                                        picture_control_set_ptr->parent_pcs_ptr->av1_cm,                  //const Av1Common *cm,
                                        plane ? blk_geom->bwidth_uv : blk_geom->bwidth,                  //int32_t wpx,
                                        plane ? blk_geom->bheight_uv : blk_geom->bheight,                  //int32_t hpx,
                                        tx_size,
                                        mode,                                                       //PredictionMode mode,
#if SEARCH_UV_MODE // conformance
                                        plane ? pu_ptr->angle_delta[PLANE_TYPE_UV] : pu_ptr->angle_delta[PLANE_TYPE_Y],
#else
                                        plane ? 0 : pu_ptr->angle_delta[PLANE_TYPE_Y],                //int32_t angle_delta,
#endif
                                        0,                                                          //int32_t use_palette,
                                        FILTER_INTRA_MODES,                                         //CHKN FilterIntraMode filter_intra_mode,
                                        topNeighArray + 1,
                                        leftNeighArray + 1,
                                        recon_buffer,                                                //uint8_t *dst,
                                        //int32_t dst_stride,
                                        0,                                                          //int32_t col_off,
                                        0,                                                          //int32_t row_off,
                                        plane,                                                      //int32_t plane,
                                        blk_geom->bsize,                  //uint32_t puSize,
                                        context_ptr->cu_origin_x,  //uint32_t cuOrgX,
                                        context_ptr->cu_origin_y);   //uint32_t cuOrgY
                                }
                            }
                            else {
                                uint8_t    topNeighArray[64 * 2 + 1];
                                uint8_t    leftNeighArray[64 * 2 + 1];
                                PredictionMode mode;
                                // Partition Loop
                                int32_t plane_end = blk_geom->has_uv ? 2 : 0;

                                for (int32_t plane = 0; plane <= plane_end; ++plane) {
#if ATB_SUPPORT
                                    TxSize  tx_size = plane ? blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr] : blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr];
#else
                                    TxSize  tx_size = plane ? blk_geom->txsize_uv[context_ptr->txb_itr] : blk_geom->txsize[context_ptr->txb_itr];
#endif
                                    if (plane == 0) {
                                        if (context_ptr->cu_origin_y != 0)
                                            memcpy(topNeighArray + 1, ep_luma_recon_neighbor_array->top_array + context_ptr->cu_origin_x, blk_geom->bwidth * 2);

                                        if (context_ptr->cu_origin_x != 0)
                                            memcpy(leftNeighArray + 1, ep_luma_recon_neighbor_array->left_array + context_ptr->cu_origin_y, blk_geom->bheight * 2);

                                        if (context_ptr->cu_origin_y != 0 && context_ptr->cu_origin_x != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ep_luma_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE + context_ptr->cu_origin_x - context_ptr->cu_origin_y];
                                    }

                                    else if (plane == 1) {
                                        if (cu_originy_uv != 0)
                                            memcpy(topNeighArray + 1, ep_cb_recon_neighbor_array->top_array + cu_originx_uv, blk_geom->bwidth_uv * 2);

                                        if (cu_originx_uv != 0)
                                            memcpy(leftNeighArray + 1, ep_cb_recon_neighbor_array->left_array + cu_originy_uv, blk_geom->bheight_uv * 2);

                                        if (cu_originy_uv != 0 && cu_originx_uv != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ep_cb_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv];
                                    }
                                    else {
                                        if (cu_originy_uv != 0)
                                            memcpy(topNeighArray + 1, ep_cr_recon_neighbor_array->top_array + cu_originx_uv, blk_geom->bwidth_uv * 2);

                                        if (cu_originx_uv != 0)
                                            memcpy(leftNeighArray + 1, ep_cr_recon_neighbor_array->left_array + cu_originy_uv, blk_geom->bheight_uv * 2);

                                        if (cu_originy_uv != 0 && cu_originx_uv != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ep_cr_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv];
                                    }

                                    if (plane)
                                        mode = (pu_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode)UV_DC_PRED : (PredictionMode)pu_ptr->intra_chroma_mode;
                                    else
                                        mode = cu_ptr->pred_mode; //PredictionMode mode,
                                    // Hsan: if CHROMA_MODE_2, then CFL will be evaluated @ EP as no CHROMA @ MD
                                    // If that's the case then you should ensure than the 1st chroma prediction uses UV_DC_PRED (that's the default configuration for CHROMA_MODE_2 if CFL applicable (set @ fast loop candidates injection) then MD assumes chroma mode always UV_DC_PRED)
                                    av1_predict_intra_block(
                                        &sb_ptr->tile_info,
                                        ED_STAGE,
                                        context_ptr->blk_geom,
                                        picture_control_set_ptr->parent_pcs_ptr->av1_cm,                  //const Av1Common *cm,
                                        plane ? blk_geom->bwidth_uv : blk_geom->bwidth,                   //int32_t wpx,
                                        plane ? blk_geom->bheight_uv : blk_geom->bheight,                  //int32_t hpx,
                                        tx_size,
                                        mode,                                                       //PredictionMode mode,
#if SEARCH_UV_MODE // conformance
                                        plane ? pu_ptr->angle_delta[PLANE_TYPE_UV] : pu_ptr->angle_delta[PLANE_TYPE_Y],
#else
                                        plane ? 0 : pu_ptr->angle_delta[PLANE_TYPE_Y],                //int32_t angle_delta,
#endif
                                        0,                                                          //int32_t use_palette,
                                        FILTER_INTRA_MODES,                                         //CHKN FilterIntraMode filter_intra_mode,
                                        topNeighArray + 1,
                                        leftNeighArray + 1,
                                        recon_buffer,                                                //uint8_t *dst,
                                        //int32_t dst_stride,
                                        0,                                                          //int32_t col_off,
                                        0,                                                          //int32_t row_off,
                                        plane,                                                      //int32_t plane,
                                        blk_geom->bsize,                  //uint32_t puSize,
#if ATB_EP
                                        context_ptr->cu_origin_x,
                                        context_ptr->cu_origin_y,
#endif
                                        context_ptr->cu_origin_x,
                                        context_ptr->cu_origin_y,
                                        0,  // MD ONLY - NOT USED BY ENCDEC
                                        0);
                                }
                                }
                            }

                            // Encode Transform Unit -INTRA-
                            {
                                uint8_t             cb_qp = cu_ptr->qp;

                                Av1EncodeLoopFunctionTable[is16bit](
                                    picture_control_set_ptr,
                                    context_ptr,
                                    sb_ptr,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    cb_qp,
                                    recon_buffer,
                                    coeff_buffer_sb,
                                    residual_buffer,
                                    transform_buffer,
                                    inverse_quant_buffer,
                                    transform_inner_array_ptr,
                                    asm_type,
                                    count_non_zero_coeffs,
                                    blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                                    useDeltaQpSegments,
                                    cu_ptr->delta_qp > 0 ? 0 : dZoffset,
                                    eobs[context_ptr->txb_itr],
                                    cuPlane);

#if  CABAC_UP
                                if(allow_update_cdf)
                                {
                                    ModeDecisionCandidateBuffer         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
                                    ModeDecisionCandidateBuffer         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[0]);
                                    ModeDecisionCandidateBuffer          *candidateBuffer;

                                    // Set the Candidate Buffer
                                    candidateBuffer = candidate_buffer_ptr_array[0];
                                    // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
#if !ATB_TX_TYPE_SUPPORT_PER_TU
                                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y];
                                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV];
#endif
                                    candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
                                    candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

                                    const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

                                    av1_tu_estimate_coeff_bits(
                                        1,//allow_update_cdf,
                                        &picture_control_set_ptr->ec_ctx_array[tbAddr],
                                        picture_control_set_ptr,
 #if ATB_DC_CONTEXT_SUPPORT_0
                                        context_ptr->txb_itr,
#endif
                                        candidateBuffer,
                                        cu_ptr,
                                        coeff1dOffset,
                                        context_ptr->coded_area_sb_uv,
                                        coeff_est_entropy_coder_ptr,
                                        coeff_buffer_sb,
                                        eobs[context_ptr->txb_itr][0],
                                        eobs[context_ptr->txb_itr][1],
                                        eobs[context_ptr->txb_itr][2],
                                        &y_tu_coeff_bits,
                                        &cb_tu_coeff_bits,
                                        &cr_tu_coeff_bits,
#if ATB_SUPPORT
                                        context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
                                        context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                                        context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                                        context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
#if ATB_TX_TYPE_SUPPORT_PER_TU
                                        cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y],
                                        cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV],
#endif
                                        context_ptr->blk_geom->has_uv ? COMPONENT_ALL : COMPONENT_LUMA,
                                        asm_type);
                                }
#endif
                                //intra mode
                                Av1EncodeGenerateReconFunctionPtr[is16bit](
                                    context_ptr,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    recon_buffer,
                                    inverse_quant_buffer,
                                    transform_inner_array_ptr,
                                    blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                                    eobs[context_ptr->txb_itr],
                                    asm_type);
                            }

                            // Update the Intra-specific Neighbor Arrays
                            EncodePassUpdateIntraModeNeighborArrays(
#if DC_SIGN_CONTEXT_EP
                                context_ptr,
                                ep_luma_dc_sign_level_coeff_neighbor_array,
                                ep_cb_dc_sign_level_coeff_neighbor_array,
                                ep_cr_dc_sign_level_coeff_neighbor_array,
#endif
                                ep_mode_type_neighbor_array,
                                ep_intra_luma_mode_neighbor_array,
                                ep_intra_chroma_mode_neighbor_array,
                                (uint8_t)cu_ptr->pred_mode,
                                (uint8_t)pu_ptr->intra_chroma_mode,
                                context_ptr->cu_origin_x,
                                context_ptr->cu_origin_y,
                                context_ptr->blk_geom->bwidth,
                                context_ptr->blk_geom->bheight,
                                context_ptr->blk_geom->bwidth_uv,
                                context_ptr->blk_geom->bheight_uv,
                                blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK);

                            // Update Recon Samples-INTRA-
                            EncodePassUpdateReconSampleNeighborArrays(
                                ep_luma_recon_neighbor_array,
                                ep_cb_recon_neighbor_array,
                                ep_cr_recon_neighbor_array,
                                recon_buffer,
                                context_ptr->cu_origin_x,
                                context_ptr->cu_origin_y,
                                context_ptr->blk_geom->bwidth,
                                context_ptr->blk_geom->bheight,
                                context_ptr->blk_geom->bwidth_uv,
                                context_ptr->blk_geom->bheight_uv,
                                blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                                is16bit);

                            if (context_ptr->blk_geom->has_uv) {
                                cu_ptr->block_has_coeff = cu_ptr->block_has_coeff |
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff |
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff |
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff;

                                if (cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff)
                                    cu_ptr->transform_unit_array[0].u_has_coeff = EB_TRUE;
                                if (cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff)
                                    cu_ptr->transform_unit_array[0].v_has_coeff = EB_TRUE;
                            }
                            else {
                                cu_ptr->block_has_coeff = cu_ptr->block_has_coeff |
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff;
                            }
                        } // Transform Loop
#if !ATB_EP
                    } // Partition Loop
#endif
                    context_ptr->coded_area_sb += blk_geom->bwidth * blk_geom->bheight;
                    if (blk_geom->has_uv)
                        context_ptr->coded_area_sb_uv += blk_geom->bwidth_uv * blk_geom->bheight_uv;
#if ATB_EP
                    }
#endif
                }

                // Inter
                else if (cu_ptr->prediction_mode_flag == INTER_MODE) {
                    context_ptr->is_inter = 1;
#if MRP_MD
                    int8_t ref_idx_l0 = (&cu_ptr->prediction_unit_array[0])->ref_frame_index_l0;
                    int8_t ref_idx_l1 = (&cu_ptr->prediction_unit_array[0])->ref_frame_index_l1;
#if MRP_MD_UNI_DIR_BIPRED
                    MvReferenceFrame rf[2];
                    av1_set_ref_frame(rf, (&cu_ptr->prediction_unit_array[0])->ref_frame_type);
                    uint8_t list_idx0, list_idx1;
                    list_idx0 = get_list_idx(rf[0]);
                    if (rf[1] == NONE_FRAME)
                        list_idx1 = get_list_idx(rf[0]);
                    else
                        list_idx1 = get_list_idx(rf[1]);
                    EbReferenceObject* refObj0 = ref_idx_l0 >= 0 ? (EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[list_idx0][ref_idx_l0]->object_ptr : (EbReferenceObject*)EB_NULL;
                    EbReferenceObject* refObj1 = ref_idx_l1 >= 0 ? (EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[list_idx1][ref_idx_l1]->object_ptr : (EbReferenceObject*)EB_NULL;
#else
                    EbReferenceObject_t* refObj0 = ref_idx_l0 >= 0 ? (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0][ref_idx_l0]->object_ptr : (EbReferenceObject_t*)EB_NULL;
                    EbReferenceObject_t* refObj1 = ref_idx_l1 >= 0 ? (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx_l1]->object_ptr : (EbReferenceObject_t*)EB_NULL;
#endif
#else
                    EbReferenceObject* refObj0 = (EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;
                    EbReferenceObject* refObj1 = picture_control_set_ptr->slice_type == B_SLICE ?
                        (EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr : 0;
#endif
                    uint16_t  txb_origin_x;
                    uint16_t  txb_origin_y;
                    EbBool isCuSkip = EB_FALSE;

                    //********************************
                    //        INTER
                    //********************************

                    EbBool  zeroLumaCbfMD = EB_FALSE;
                    //EbBool doLumaMC = EB_TRUE;
                    EbBool doMVpred = EB_TRUE;
                    //if QPM and Segments are used, First Cu in SB row should have at least one coeff.
                    EbBool isFirstCUinRow = (use_delta_qp == 1) &&
                        !oneSegment &&
                        (context_ptr->cu_origin_x == 0 && context_ptr->cu_origin_y == sb_origin_y) ? EB_TRUE : EB_FALSE;
                    zeroLumaCbfMD = (EbBool)(checkZeroLumaCbf && ((&cu_ptr->prediction_unit_array[0])->merge_flag == EB_FALSE && cu_ptr->block_has_coeff == 0 && isFirstCUinRow == EB_FALSE));
                    zeroLumaCbfMD = EB_FALSE;

                    //Motion Compensation could be avoided in the case below
                    EbBool doMC = EB_TRUE;

                    // Perform Merge/Skip Decision if the mode coming from MD is merge. for the First CU in Row merge will remain as is.
                    if (cu_ptr->prediction_unit_array[0].merge_flag == EB_TRUE)
                    {
                        if (isFirstCUinRow == EB_FALSE)
                            isCuSkip = mdcontextPtr->md_ep_pipe_sb[cu_ptr->mds_idx].skip_cost <= mdcontextPtr->md_ep_pipe_sb[cu_ptr->mds_idx].merge_cost ? 1 : 0;
                    }

                    //MC could be avoided in some cases below
                    if (isFirstCUinRow == EB_FALSE) {
                        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_FALSE && constrained_intra_flag == EB_TRUE &&
                            cu_ptr->prediction_unit_array[0].merge_flag == EB_TRUE)
                        {
                            if (isCuSkip)
                            {
                                //here merge is decided to be skip in nonRef frame.
                                doMC = EB_FALSE;
                                doMVpred = EB_FALSE;
                            }
                        }
                        else if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_FALSE && constrained_intra_flag == EB_TRUE &&
                            zeroLumaCbfMD == EB_TRUE)
                        {
                            //MV mode with no Coeff  in nonRef frame.
                            doMC = EB_FALSE;
                        }

                        else if (picture_control_set_ptr->limit_intra && isIntraLCU == EB_FALSE)
                        {
                            if (isCuSkip)
                            {
                                doMC = EB_FALSE;
                                doMVpred = EB_FALSE;
                            }
                        }
                    }

                    doMC = (EbBool)(doRecon | doMC);

                    doMVpred = (EbBool)(doRecon | doMVpred);

                    //IntMv  predmv[2];
                    enc_pass_av1_mv_pred(
                        &sb_ptr->tile_info,
                         context_ptr->md_context,
                        cu_ptr,
                        blk_geom,
                        context_ptr->cu_origin_x,
                        context_ptr->cu_origin_y,
                        picture_control_set_ptr,
                        cu_ptr->prediction_unit_array[0].ref_frame_type,
                        cu_ptr->prediction_unit_array[0].is_compound,
                        cu_ptr->pred_mode,
                        cu_ptr->predmv);
                    //out1:  predmv
                    //out2:   cu_ptr->inter_mode_ctx[ cu_ptr->prediction_unit_array[0].ref_frame_type ]

                    //keep final usefull mvp for entropy
                    memcpy(cu_ptr->av1xd->final_ref_mv_stack,
                       context_ptr->md_context->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[cu_ptr->prediction_unit_array[0].ref_frame_type],
                        sizeof(CandidateMv)*MAX_REF_MV_STACK_SIZE);

                    {
                        // 1st Partition Loop
                        pu_ptr = cu_ptr->prediction_unit_array;

                        // Set MvUnit
                        context_ptr->mv_unit.pred_direction = (uint8_t)pu_ptr->inter_pred_direction_index;
                        context_ptr->mv_unit.mv[REF_LIST_0].mv_union = pu_ptr->mv[REF_LIST_0].mv_union;
                        context_ptr->mv_unit.mv[REF_LIST_1].mv_union = pu_ptr->mv[REF_LIST_1].mv_union;

                        // Inter Prediction
                        if (doMC &&
                            pu_ptr->motion_mode == WARPED_CAUSAL)
                        {
                            warped_motion_prediction(
                                &context_ptr->mv_unit,
                                context_ptr->cu_origin_x,
                                context_ptr->cu_origin_y,
                                cu_ptr,
                                blk_geom,
                                is16bit ? refObj0->reference_picture16bit : refObj0->reference_picture,
                                recon_buffer,
                                context_ptr->cu_origin_x,
                                context_ptr->cu_origin_y,
                                &cu_ptr->prediction_unit_array[0].wm_params,
                                (uint8_t) sequence_control_set_ptr->static_config.encoder_bit_depth,
                                EB_TRUE,
                                asm_type);
                        }

                        if (doMC &&
                            pu_ptr->motion_mode != WARPED_CAUSAL)
                        {
                            if (is16bit) {
                                av1_inter_prediction_hbd(
                                    picture_control_set_ptr,
                                    cu_ptr->prediction_unit_array->ref_frame_type,
                                    cu_ptr,
                                    &context_ptr->mv_unit,
                                    0,// use_intrabc,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    blk_geom->bwidth,
                                    blk_geom->bheight,
#if FIXED_MRP_10BIT
                                    cu_ptr->prediction_unit_array->ref_frame_index_l0 >= 0 ? refObj0->reference_picture16bit : (EbPictureBufferDesc*)EB_NULL,
                                    cu_ptr->prediction_unit_array->ref_frame_index_l1 >= 0 ? refObj1->reference_picture16bit : (EbPictureBufferDesc*)EB_NULL,
#else
                                    refObj0->reference_picture16bit,
                                    picture_control_set_ptr->slice_type == B_SLICE ? refObj1->reference_picture16bit : 0,
#endif
                                    recon_buffer,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    (uint8_t)sequence_control_set_ptr->static_config.encoder_bit_depth,
                                    asm_type);
                            } else {
                                av1_inter_prediction(
                                    picture_control_set_ptr,
                                    cu_ptr->interp_filters,
                                    cu_ptr,
                                    cu_ptr->prediction_unit_array->ref_frame_type,
                                    &context_ptr->mv_unit,
                                    0,//use_intrabc,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    blk_geom->bwidth,
                                    blk_geom->bheight,
#if MRP_MD
                                    cu_ptr->prediction_unit_array->ref_frame_index_l0 >= 0 ? refObj0->reference_picture : (EbPictureBufferDesc*)EB_NULL,
                                    cu_ptr->prediction_unit_array->ref_frame_index_l1 >= 0 ? refObj1->reference_picture : (EbPictureBufferDesc*)EB_NULL,
#else
                                    refObj0->reference_picture,
                                    picture_control_set_ptr->slice_type == B_SLICE ? refObj1->reference_picture : 0,
#endif
                                    recon_buffer,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    EB_TRUE,
                                    asm_type);
                            }
                        }
                    }

                    context_ptr->txb_itr = 0;
                    // Transform Loop
                    cu_ptr->transform_unit_array[0].y_has_coeff = EB_FALSE;
                    cu_ptr->transform_unit_array[0].u_has_coeff = EB_FALSE;
                    cu_ptr->transform_unit_array[0].v_has_coeff = EB_FALSE;

                    // initialize TU Split
                    y_full_distortion[DIST_CALC_RESIDUAL] = 0;
                    y_full_distortion[DIST_CALC_PREDICTION] = 0;

                    y_coeff_bits = 0;
                    cb_coeff_bits = 0;
                    cr_coeff_bits = 0;

#if ATB_SUPPORT
                    uint32_t totTu = context_ptr->blk_geom->txb_count[cu_ptr->tx_depth];
#else
                    uint32_t totTu = context_ptr->blk_geom->txb_count;
#endif
                    uint8_t   tuIt;
                    uint8_t   cb_qp = cu_ptr->qp;
                    uint32_t  component_mask = context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK;

                    if (cu_ptr->prediction_unit_array[0].merge_flag == EB_FALSE) {
                        for (uint8_t tuIt = 0; tuIt < totTu; tuIt++) {
                            context_ptr->txb_itr = tuIt;
#if ATB_SUPPORT
                            uint8_t uv_pass = cu_ptr->tx_depth && tuIt ? 0 : 1; //NM: 128x128 exeption
#endif
#if ATB_SUPPORT
                            txb_origin_x = context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[cu_ptr->tx_depth][tuIt];
                            txb_origin_y = context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[cu_ptr->tx_depth][tuIt];
#else
                            txb_origin_x = context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[tuIt];
                            txb_origin_y = context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[tuIt];
#endif

#if DC_SIGN_CONTEXT_EP
                            uint32_t cu_originy_uv = (context_ptr->cu_origin_y >> 3 << 3) >> 1;
                            uint32_t cu_originx_uv = (context_ptr->cu_origin_x >> 3 << 3) >> 1;

                            context_ptr->cu_ptr->luma_txb_skip_context = 0;
                            context_ptr->cu_ptr->luma_dc_sign_context[context_ptr->txb_itr] = 0;
                            get_txb_ctx(
                                COMPONENT_LUMA,
                                picture_control_set_ptr->ep_luma_dc_sign_level_coeff_neighbor_array,
                                txb_origin_x,
                                txb_origin_y,
                                context_ptr->blk_geom->bsize,
                                context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
                                &context_ptr->cu_ptr->luma_txb_skip_context,
                                &context_ptr->cu_ptr->luma_dc_sign_context[context_ptr->txb_itr]);

                            cu_ptr->cb_txb_skip_context = 0;
                            cu_ptr->cb_dc_sign_context = 0;
                            get_txb_ctx(
                                COMPONENT_CHROMA,
                                picture_control_set_ptr->ep_cb_dc_sign_level_coeff_neighbor_array,
                                cu_originx_uv,
                                cu_originy_uv,
                                context_ptr->blk_geom->bsize_uv,
                                context_ptr->blk_geom->txsize_uv[context_ptr->cu_ptr->tx_depth][context_ptr->txb_itr],
                                &cu_ptr->cb_txb_skip_context,
                                &cu_ptr->cb_dc_sign_context);

                            cu_ptr->cr_txb_skip_context = 0;
                            cu_ptr->cr_dc_sign_context = 0;
                            get_txb_ctx(
                                COMPONENT_CHROMA,
                                picture_control_set_ptr->ep_cr_dc_sign_level_coeff_neighbor_array,
                                cu_originx_uv,
                                cu_originy_uv,
                                context_ptr->blk_geom->bsize_uv,
                                context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
                                &cu_ptr->cr_txb_skip_context,
                                &cu_ptr->cr_dc_sign_context);
#endif
                            if (!zeroLumaCbfMD)
                                //inter mode  1
                                Av1EncodeLoopFunctionTable[is16bit](
                                    picture_control_set_ptr,
                                    context_ptr,
                                    sb_ptr,
                                    txb_origin_x,   //pic org
                                    txb_origin_y,
                                    cb_qp,
                                    recon_buffer,
                                    coeff_buffer_sb,
                                    residual_buffer,
                                    transform_buffer,
                                    inverse_quant_buffer,
                                    transform_inner_array_ptr,
                                    asm_type,
                                    count_non_zero_coeffs,
#if ATB_SUPPORT
                                    context_ptr->blk_geom->has_uv && uv_pass ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
#else
                                    context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
#endif
                                    useDeltaQpSegments,
                                    cu_ptr->delta_qp > 0 ? 0 : dZoffset,
                                    eobs[context_ptr->txb_itr],
                                    cuPlane);

                            // SKIP the CBF zero mode for DC path. There are problems with cost calculations
#if PF_N2_SUPPORT
                            {
#else
                            if (context_ptr->trans_coeff_shape_luma != ONLY_DC_SHAPE) {
#endif
                                // Compute Tu distortion
                                if (!zeroLumaCbfMD)

                                    // LUMA DISTORTION
                                    picture_full_distortion32_bits(
                                        transform_buffer,
                                        context_ptr->coded_area_sb,
                                        0,
                                        inverse_quant_buffer,
                                        context_ptr->coded_area_sb,
                                        0,
#if ATB_SUPPORT
                                        blk_geom->tx_width[cu_ptr->tx_depth][tuIt],
                                        blk_geom->tx_height[cu_ptr->tx_depth][tuIt],
#else
                                        blk_geom->tx_width[tuIt],
                                        blk_geom->tx_height[tuIt],
#endif
                                        context_ptr->blk_geom->bwidth_uv,
                                        context_ptr->blk_geom->bheight_uv,
                                        yTuFullDistortion,
                                        yTuFullDistortion,
                                        yTuFullDistortion,
                                        eobs[context_ptr->txb_itr][0],
                                        0,
                                        0,
                                        COMPONENT_LUMA,
                                        asm_type);
#if ATB_SUPPORT
                                TxSize  txSize = blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr];
#else
                                TxSize  txSize = blk_geom->txsize[context_ptr->txb_itr];
#endif
                                int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
                                yTuFullDistortion[DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(yTuFullDistortion[DIST_CALC_RESIDUAL], shift);
                                yTuFullDistortion[DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(yTuFullDistortion[DIST_CALC_PREDICTION], shift);

                                y_tu_coeff_bits = 0;
                                cb_tu_coeff_bits = 0;
                                cr_tu_coeff_bits = 0;

                                if (!zeroLumaCbfMD) {
                                    ModeDecisionCandidateBuffer         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
                                    ModeDecisionCandidateBuffer         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[0]);
                                    ModeDecisionCandidateBuffer          *candidateBuffer;

                                    // Set the Candidate Buffer
                                    candidateBuffer = candidate_buffer_ptr_array[0];
                                    // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
#if !ATB_TX_TYPE_SUPPORT_PER_TU
                                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[tuIt].transform_type[PLANE_TYPE_Y];
                                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[tuIt].transform_type[PLANE_TYPE_UV];
#endif
                                    candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;

                                    const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

                                    av1_tu_estimate_coeff_bits(
#if CABAC_UP
                                        0,//allow_update_cdf,
                                        NULL,
#endif
                                        picture_control_set_ptr,
#if ATB_DC_CONTEXT_SUPPORT_0
                                        tuIt,
#endif
                                        candidateBuffer,
                                        cu_ptr,
                                        coeff1dOffset,
                                        context_ptr->coded_area_sb_uv,
                                        coeff_est_entropy_coder_ptr,
                                        coeff_buffer_sb,
                                        eobs[context_ptr->txb_itr][0],
                                        eobs[context_ptr->txb_itr][1],
                                        eobs[context_ptr->txb_itr][2],
                                        &y_tu_coeff_bits,
                                        &cb_tu_coeff_bits,
                                        &cr_tu_coeff_bits,
#if ATB_SUPPORT
                                        context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
                                        context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                                        context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                                        context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
#if ATB_TX_TYPE_SUPPORT_PER_TU
                                        cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y],
                                        cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV],
#endif
#if ATB_SUPPORT
                                        context_ptr->blk_geom->has_uv && uv_pass ? COMPONENT_ALL : COMPONENT_LUMA,
#else
                                        context_ptr->blk_geom->has_uv ? COMPONENT_ALL : COMPONENT_LUMA,
#endif
                                        asm_type);
                                }

                                // CBF Tu decision
                                if (zeroLumaCbfMD == EB_FALSE)

                                    av1_encode_tu_calc_cost(
                                        context_ptr,
                                        count_non_zero_coeffs,
                                        yTuFullDistortion,
                                        &y_tu_coeff_bits,
                                        component_mask);

                                else {
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff = 0;
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff = 0;
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff = 0;
                                }
                                // Update count_non_zero_coeffs after CBF decision
                                if (cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff == EB_FALSE)
                                    count_non_zero_coeffs[0] = 0;
#if ATB_SUPPORT
                                if (context_ptr->blk_geom->has_uv && uv_pass) {
#else
                                if (context_ptr->blk_geom->has_uv) {
#endif
                                    if (cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff == EB_FALSE)
                                        count_non_zero_coeffs[1] = 0;
                                    if (cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff == EB_FALSE)
                                        count_non_zero_coeffs[2] = 0;
                                }

                                // Update TU count_non_zero_coeffs
                                cu_ptr->transform_unit_array[context_ptr->txb_itr].nz_coef_count[0] = (uint16_t)count_non_zero_coeffs[0];
                                cu_ptr->transform_unit_array[context_ptr->txb_itr].nz_coef_count[1] = (uint16_t)count_non_zero_coeffs[1];
                                cu_ptr->transform_unit_array[context_ptr->txb_itr].nz_coef_count[2] = (uint16_t)count_non_zero_coeffs[2];

                                y_coeff_bits += y_tu_coeff_bits;
#if ATB_SUPPORT
                                if (context_ptr->blk_geom->has_uv && uv_pass) {
#else
                                if (context_ptr->blk_geom->has_uv) {
#endif
                                    cb_coeff_bits += cb_tu_coeff_bits;
                                    cr_coeff_bits += cr_tu_coeff_bits;
                                }

                                y_full_distortion[DIST_CALC_RESIDUAL] += yTuFullDistortion[DIST_CALC_RESIDUAL];
                                y_full_distortion[DIST_CALC_PREDICTION] += yTuFullDistortion[DIST_CALC_PREDICTION];

#if CABAC_UP
                                if (allow_update_cdf) {
                                    ModeDecisionCandidateBuffer         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
                                    ModeDecisionCandidateBuffer         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[0]);
                                    ModeDecisionCandidateBuffer          *candidateBuffer;

                                    // Set the Candidate Buffer
                                    candidateBuffer = candidate_buffer_ptr_array[0];
                                    // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
#if !ATB_TX_TYPE_SUPPORT_PER_TU
                                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[tuIt].transform_type[PLANE_TYPE_Y];
                                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[tuIt].transform_type[PLANE_TYPE_UV];
#endif
                                    candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
                                    candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

                                    const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

                                    //CHKN add updating eobs[] after CBF decision
                                    if (cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff == EB_FALSE)
                                        eobs[context_ptr->txb_itr][0] = 0;
#if ATB_SUPPORT
                                    if (context_ptr->blk_geom->has_uv && uv_pass) {
#else
                                    if (context_ptr->blk_geom->has_uv) {
#endif
                                        if (cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff == EB_FALSE)
                                            eobs[context_ptr->txb_itr][1] = 0;
                                        if (cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff == EB_FALSE)
                                            eobs[context_ptr->txb_itr][2] = 0;
                                    }

                                    av1_tu_estimate_coeff_bits(
                                        1,//allow_update_cdf,
                                        &picture_control_set_ptr->ec_ctx_array[tbAddr],
                                        picture_control_set_ptr,
#if ATB_DC_CONTEXT_SUPPORT_0
                                        tuIt,
#endif
                                        candidateBuffer,
                                        cu_ptr,
                                        coeff1dOffset,
                                        context_ptr->coded_area_sb_uv,
                                        coeff_est_entropy_coder_ptr,
                                        coeff_buffer_sb,
                                        eobs[context_ptr->txb_itr][0],
                                        eobs[context_ptr->txb_itr][1],
                                        eobs[context_ptr->txb_itr][2],
                                        &y_tu_coeff_bits,
                                        &cb_tu_coeff_bits,
                                        &cr_tu_coeff_bits,
#if ATB_SUPPORT
                                        context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
                                        context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                                        context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                                        context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
#if ATB_TX_TYPE_SUPPORT_PER_TU
                                        cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y],
                                        cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV],
#endif
#if ATB_SUPPORT
                                        context_ptr->blk_geom->has_uv && uv_pass ? COMPONENT_ALL : COMPONENT_LUMA,
#else
                                        context_ptr->blk_geom->has_uv ? COMPONENT_ALL : COMPONENT_LUMA,
#endif
                                        asm_type);
                                }
#endif
                            }
#if ATB_SUPPORT
                            context_ptr->coded_area_sb += blk_geom->tx_width[cu_ptr->tx_depth][tuIt] * blk_geom->tx_height[cu_ptr->tx_depth][tuIt];
#if ATB_SUPPORT
                            if (context_ptr->blk_geom->has_uv && uv_pass)
#else
                            if (blk_geom->has_uv)
#endif
                                context_ptr->coded_area_sb_uv += blk_geom->tx_width_uv[cu_ptr->tx_depth][tuIt] * blk_geom->tx_height_uv[cu_ptr->tx_depth][tuIt];
#else
                            context_ptr->coded_area_sb += blk_geom->tx_width[tuIt] * blk_geom->tx_height[tuIt];
                            if (blk_geom->has_uv)
                                context_ptr->coded_area_sb_uv += blk_geom->tx_width_uv[tuIt] * blk_geom->tx_height_uv[tuIt];
#endif
                        } // Transform Loop
                    }

                    //Set Final CU data flags after skip/Merge decision.
                    if (isFirstCUinRow == EB_FALSE) {
                        if (cu_ptr->prediction_unit_array[0].merge_flag == EB_TRUE) {
                            cu_ptr->skip_flag = (isCuSkip) ? EB_TRUE : EB_FALSE;
                            cu_ptr->prediction_unit_array[0].merge_flag = (isCuSkip) ? EB_FALSE : EB_TRUE;
                        }
                    }

                    // Initialize the Transform Loop

                    context_ptr->txb_itr = 0;
                    y_has_coeff = 0;
                    u_has_coeff = 0;
                    v_has_coeff = 0;
#if ATB_SUPPORT
                    totTu = context_ptr->blk_geom->txb_count[cu_ptr->tx_depth];
#else
                    totTu = context_ptr->blk_geom->txb_count;
#endif

                    //reset coeff buffer offsets at the start of a new Tx loop
                    context_ptr->coded_area_sb = coded_area_org;
                    context_ptr->coded_area_sb_uv = coded_area_org_uv;
                    for (tuIt = 0; tuIt < totTu; tuIt++)
                    {
#if ATB_SUPPORT
                        uint8_t uv_pass = cu_ptr->tx_depth && tuIt ? 0 : 1; //NM: 128x128 exeption
#endif
                        context_ptr->txb_itr = tuIt;
#if ATB_SUPPORT
                        txb_origin_x = context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[cu_ptr->tx_depth][tuIt];
                        txb_origin_y = context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[cu_ptr->tx_depth][tuIt];
#else
                        txb_origin_x = context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[tuIt];
                        txb_origin_y = context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[tuIt];
#endif
                        if (cu_ptr->skip_flag == EB_TRUE) {
                            cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff = EB_FALSE;
                            cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff = EB_FALSE;
                            cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff = EB_FALSE;
                        }
                        else if ((&cu_ptr->prediction_unit_array[0])->merge_flag == EB_TRUE) {
                            //inter mode  2

                            Av1EncodeLoopFunctionTable[is16bit](
                                picture_control_set_ptr,
                                context_ptr,
                                sb_ptr,
                                txb_origin_x, //pic offset
                                txb_origin_y,
                                cb_qp,
                                recon_buffer,
                                coeff_buffer_sb,
                                residual_buffer,
                                transform_buffer,
                                inverse_quant_buffer,
                                transform_inner_array_ptr,
                                asm_type,
                                count_non_zero_coeffs,
#if ATB_SUPPORT
                                context_ptr->blk_geom->has_uv && uv_pass ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
#else
                                context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
#endif
                                useDeltaQpSegments,
                                cu_ptr->delta_qp > 0 ? 0 : dZoffset,
                                eobs[context_ptr->txb_itr],
                                cuPlane);

#if CABAC_UP
                            if (allow_update_cdf) {
                                ModeDecisionCandidateBuffer         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
                                ModeDecisionCandidateBuffer         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[0]);
                                ModeDecisionCandidateBuffer          *candidateBuffer;

                                // Set the Candidate Buffer
                                candidateBuffer = candidate_buffer_ptr_array[0];
                                // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
#if !ATB_TX_TYPE_SUPPORT_PER_TU
                                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[tuIt].transform_type[PLANE_TYPE_Y];
                                candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[tuIt].transform_type[PLANE_TYPE_UV];
#endif
                                candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;
                                candidateBuffer->candidate_ptr->pred_mode = cu_ptr->pred_mode;

                                const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

                                av1_tu_estimate_coeff_bits(
                                    1,//allow_update_cdf,
                                    &picture_control_set_ptr->ec_ctx_array[tbAddr],
                                    picture_control_set_ptr,
#if ATB_DC_CONTEXT_SUPPORT_0
                                    tuIt,
#endif
                                    candidateBuffer,
                                    cu_ptr,
                                    coeff1dOffset,
                                    context_ptr->coded_area_sb_uv,
                                    coeff_est_entropy_coder_ptr,
                                    coeff_buffer_sb,
                                    eobs[context_ptr->txb_itr][0],
                                    eobs[context_ptr->txb_itr][1],
                                    eobs[context_ptr->txb_itr][2],
                                    &y_tu_coeff_bits,
                                    &cb_tu_coeff_bits,
                                    &cr_tu_coeff_bits,
#if ATB_SUPPORT
                                    context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr],
                                    context_ptr->blk_geom->txsize_uv[cu_ptr->tx_depth][context_ptr->txb_itr],
#else
                                    context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                                    context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
#endif
#if ATB_TX_TYPE_SUPPORT_PER_TU
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_Y],
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].transform_type[PLANE_TYPE_UV],
#endif
#if ATB_SUPPORT
                                    context_ptr->blk_geom->has_uv && uv_pass ? COMPONENT_ALL : COMPONENT_LUMA,
#else
                                    context_ptr->blk_geom->has_uv ? COMPONENT_ALL : COMPONENT_LUMA,
#endif
                                    asm_type);
                            }
#endif
                        }
#if ATB_SUPPORT
                        if (context_ptr->blk_geom->has_uv && uv_pass) {
#else
                        if (context_ptr->blk_geom->has_uv) {
#endif
                            cu_ptr->block_has_coeff = cu_ptr->block_has_coeff |
                                cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff |
                                cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff |
                                cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff;
                        }
                        else {
                            cu_ptr->block_has_coeff = cu_ptr->block_has_coeff |
                                cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff;
                        }

                        //inter mode
                        if (doRecon)

                            Av1EncodeGenerateReconFunctionPtr[is16bit](
                                context_ptr,
                                txb_origin_x,  //pic offset
                                txb_origin_y,
                                recon_buffer,
                                inverse_quant_buffer,
                                transform_inner_array_ptr,
#if ATB_SUPPORT
                                context_ptr->blk_geom->has_uv && uv_pass ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
#else
                                context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
#endif
                                eobs[context_ptr->txb_itr],
                                asm_type);
#if ATB_SUPPORT
                        if (context_ptr->blk_geom->has_uv && uv_pass) {
#else
                        if (context_ptr->blk_geom->has_uv) {
#endif
                            y_has_coeff |= cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff;
                            u_has_coeff |= cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff;
                            v_has_coeff |= cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff;
                        }
                        else
                            y_has_coeff |= cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff;
#if ATB_SUPPORT
                        context_ptr->coded_area_sb += blk_geom->tx_width[cu_ptr->tx_depth][tuIt] * blk_geom->tx_height[cu_ptr->tx_depth][tuIt];

#if ATB_SUPPORT
                        if (context_ptr->blk_geom->has_uv && uv_pass)
#else
                        if (blk_geom->has_uv)
#endif
                            context_ptr->coded_area_sb_uv += blk_geom->tx_width_uv[cu_ptr->tx_depth][tuIt] * blk_geom->tx_height_uv[cu_ptr->tx_depth][tuIt];
#else
                        context_ptr->coded_area_sb += blk_geom->tx_width[tuIt] * blk_geom->tx_height[tuIt];
                        if (blk_geom->has_uv)
                            context_ptr->coded_area_sb_uv += blk_geom->tx_width_uv[tuIt] * blk_geom->tx_height_uv[tuIt];
#endif
                    } // Transform Loop

                    // Calculate Root CBF
                    if (context_ptr->blk_geom->has_uv)
                        cu_ptr->block_has_coeff = (y_has_coeff | u_has_coeff | v_has_coeff) ? EB_TRUE : EB_FALSE;
                    else
                        cu_ptr->block_has_coeff = (y_has_coeff) ? EB_TRUE : EB_FALSE;

                    // Force Skip if MergeFlag == TRUE && RootCbf == 0

                    if (cu_ptr->skip_flag == EB_FALSE &&
                        cu_ptr->prediction_unit_array[0].merge_flag == EB_TRUE && cu_ptr->block_has_coeff == EB_FALSE)
                    {
                        cu_ptr->skip_flag = EB_TRUE;
                    }

                    {
                        // Set the PU Loop Variables
                        pu_ptr = cu_ptr->prediction_unit_array;

                        // Set MvUnit
                        context_ptr->mv_unit.pred_direction = (uint8_t)pu_ptr->inter_pred_direction_index;
                        context_ptr->mv_unit.mv[REF_LIST_0].mv_union = pu_ptr->mv[REF_LIST_0].mv_union;
                        context_ptr->mv_unit.mv[REF_LIST_1].mv_union = pu_ptr->mv[REF_LIST_1].mv_union;

                        // Update Neighbor Arrays (Mode Type, mvs, SKIP)
                        {
                            uint8_t skip_flag = (uint8_t)cu_ptr->skip_flag;
                            EncodePassUpdateInterModeNeighborArrays(
#if DC_SIGN_CONTEXT_EP
                                context_ptr,
                                ep_luma_dc_sign_level_coeff_neighbor_array,
                                ep_cb_dc_sign_level_coeff_neighbor_array,
                                ep_cr_dc_sign_level_coeff_neighbor_array,
#endif
                                ep_mode_type_neighbor_array,
                                ep_mv_neighbor_array,
                                ep_skip_flag_neighbor_array,
                                &context_ptr->mv_unit,
                                &skip_flag,
                                context_ptr->cu_origin_x,
                                context_ptr->cu_origin_y,
                                blk_geom->bwidth,
                                blk_geom->bheight);
                        }
                    } // 2nd Partition Loop

                    // Update Recon Samples Neighbor Arrays -INTER-

                    if (doRecon)
                        EncodePassUpdateReconSampleNeighborArrays(
                            ep_luma_recon_neighbor_array,
                            ep_cb_recon_neighbor_array,
                            ep_cr_recon_neighbor_array,
                            recon_buffer,
                            context_ptr->cu_origin_x,
                            context_ptr->cu_origin_y,
                            context_ptr->blk_geom->bwidth,
                            context_ptr->blk_geom->bheight,
                            context_ptr->blk_geom->bwidth_uv,
                            context_ptr->blk_geom->bheight_uv,
                            context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                            is16bit);
                }
                else {
                    CHECK_REPORT_ERROR_NC(
                        encode_context_ptr->app_callback_ptr,
                        EB_ENC_CL_ERROR2);
                }

                update_av1_mi_map(
                    cu_ptr,
                    context_ptr->cu_origin_x,
                    context_ptr->cu_origin_y,
                    blk_geom,
                    picture_control_set_ptr);

                if (dlfEnableFlag)
                {
#if !MEMORY_FOOTPRINT_OPT
                    if (blk_geom->has_uv) {
                        availableCoeff = (cu_ptr->prediction_mode_flag == INTER_MODE) ? (EbBool)cu_ptr->block_has_coeff :
                            (cu_ptr->transform_unit_array[0].y_has_coeff ||
                                cu_ptr->transform_unit_array[0].v_has_coeff ||
                                cu_ptr->transform_unit_array[0].u_has_coeff) ? EB_TRUE : EB_FALSE;
                    }
                    else
                        availableCoeff = (cu_ptr->transform_unit_array[0].y_has_coeff) ? EB_TRUE : EB_FALSE;
                    // Assign the LCU-level QP
                    //NM - To be revisited
                    EncodePassUpdateQp(
                        picture_control_set_ptr,
                        context_ptr,
                        availableCoeff,
                        use_delta_qp,
                        &isDeltaQpNotCoded,
                        picture_control_set_ptr->dif_cu_delta_qp_depth,
                        &(picture_control_set_ptr->enc_prev_coded_qp[oneSegment ? 0 : lcuRowIndex]),
                        &(picture_control_set_ptr->enc_prev_quant_group_coded_qp[oneSegment ? 0 : lcuRowIndex]),
                        sb_qp);
#endif
                }

                {
                    {
                        // Set the PU Loop Variables
                        pu_ptr = cu_ptr->prediction_unit_array;
                        // Set MvUnit
                        context_ptr->mv_unit.pred_direction = (uint8_t)pu_ptr->inter_pred_direction_index;
                        context_ptr->mv_unit.mv[REF_LIST_0].mv_union = pu_ptr->mv[REF_LIST_0].mv_union;
                        context_ptr->mv_unit.mv[REF_LIST_1].mv_union = pu_ptr->mv[REF_LIST_1].mv_union;
                    }
                }

                {
                    CodingUnit *src_cu = &context_ptr->md_context->md_cu_arr_nsq[d1_itr];

                    CodingUnit *dst_cu = &sb_ptr->final_cu_arr[final_cu_itr++];

                    move_cu_data(src_cu, dst_cu);
                }
            }
            blk_it += ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][context_ptr->blk_geom->depth];
        }
        else
            blk_it += d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][context_ptr->blk_geom->depth];
    } // CU Loop
#if !MEMORY_FOOTPRINT_OPT
    sb_ptr->tot_final_cu = final_cu_itr;
#endif
#if AV1_LF
    // First Pass Deblocking
    if (dlfEnableFlag && picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode == 1) {
        if (picture_control_set_ptr->parent_pcs_ptr->lf.filter_level[0] || picture_control_set_ptr->parent_pcs_ptr->lf.filter_level[1]) {
            uint8_t LastCol = ((sb_origin_x)+sb_width == sequence_control_set_ptr->luma_width) ? 1 : 0;
            loop_filter_sb(
                recon_buffer,
                picture_control_set_ptr,
                NULL,
                sb_origin_y >> 2,
                sb_origin_x >> 2,
                0,
                3,
                LastCol);
        }
    }
#endif

    return;
}

#if NO_ENCDEC
EB_EXTERN void no_enc_dec_pass(
    SequenceControlSet    *sequence_control_set_ptr,
    PictureControlSet     *picture_control_set_ptr,
    LargestCodingUnit     *sb_ptr,
    uint32_t                   tbAddr,
    uint32_t                   sb_origin_x,
    uint32_t                   sb_origin_y,
    uint32_t                   sb_qp,
    EncDecContext         *context_ptr)
{
    context_ptr->coded_area_sb = 0;
    context_ptr->coded_area_sb_uv = 0;

    uint32_t      final_cu_itr = 0;

    uint32_t    blk_it = 0;

    while (blk_it < sequence_control_set_ptr->max_block_cnt) {
        CodingUnit  *cu_ptr = context_ptr->cu_ptr = &context_ptr->md_context->md_cu_arr_nsq[blk_it];
        PartitionType part = cu_ptr->part;
        const BlockGeom * blk_geom = context_ptr->blk_geom = get_blk_geom_mds(blk_it);

        sb_ptr->cu_partition_array[blk_it] = context_ptr->md_context->md_cu_arr_nsq[blk_it].part;

        if (part != PARTITION_SPLIT) {
            int32_t offset_d1 = ns_blk_offset[(int32_t)part]; //cu_ptr->best_d1_blk; // TOCKECK
            int32_t num_d1_block = ns_blk_num[(int32_t)part]; // context_ptr->blk_geom->totns; // TOCKECK

            for (int32_t d1_itr = blk_it + offset_d1; d1_itr < blk_it + offset_d1 + num_d1_block; d1_itr++) {
                const BlockGeom * blk_geom = context_ptr->blk_geom = get_blk_geom_mds(d1_itr);
                CodingUnit            *cu_ptr = context_ptr->cu_ptr = &context_ptr->md_context->md_cu_arr_nsq[d1_itr];

                cu_ptr->delta_qp = 0;
                cu_ptr->qp = (sequence_control_set_ptr->static_config.improve_sharpness) ? context_ptr->qpm_qp : picture_control_set_ptr->picture_qp;
                sb_ptr->qp = (sequence_control_set_ptr->static_config.improve_sharpness) ? context_ptr->qpm_qp : picture_control_set_ptr->picture_qp;
                cu_ptr->org_delta_qp = cu_ptr->delta_qp;

                {
                    CodingUnit *src_cu = &context_ptr->md_context->md_cu_arr_nsq[d1_itr];
                    CodingUnit *dst_cu = &sb_ptr->final_cu_arr[final_cu_itr++];

                    move_cu_data(src_cu, dst_cu);
                }

                //copy coeff
                int32_t txb_1d_offset = 0, txb_1d_offset_uv = 0;

                int32_t txb_itr = 0;
                do
                {
                    uint32_t  bwidth = context_ptr->blk_geom->tx_width[txb_itr] < 64 ? context_ptr->blk_geom->tx_width[txb_itr] : 32;
                    uint32_t  bheight = context_ptr->blk_geom->tx_height[txb_itr] < 64 ? context_ptr->blk_geom->tx_height[txb_itr] : 32;

                    int32_t* src_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->buffer_y)[txb_1d_offset]);
                    int32_t* dst_ptr = &(((int32_t*)sb_ptr->quantized_coeff->buffer_y)[context_ptr->coded_area_sb]);

                    uint32_t j;
                    for (j = 0; j < bheight; j++)
                        memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                    if (context_ptr->blk_geom->has_uv)
                    {
                        // Cb
                        bwidth = context_ptr->blk_geom->tx_width_uv[txb_itr];
                        bheight = context_ptr->blk_geom->tx_height_uv[txb_itr];

                        src_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->buffer_cb)[txb_1d_offset_uv]);
                        dst_ptr = &(((int32_t*)sb_ptr->quantized_coeff->buffer_cb)[context_ptr->coded_area_sb_uv]);

                        for (j = 0; j < bheight; j++)
                            memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                        //Cr
                        src_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->buffer_cr)[txb_1d_offset_uv]);
                        dst_ptr = &(((int32_t*)sb_ptr->quantized_coeff->buffer_cr)[context_ptr->coded_area_sb_uv]);

                        for (j = 0; j < bheight; j++)
                            memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                    }

                    context_ptr->coded_area_sb += context_ptr->blk_geom->tx_width[txb_itr] * context_ptr->blk_geom->tx_height[txb_itr];
                    if (context_ptr->blk_geom->has_uv)
                        context_ptr->coded_area_sb_uv += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];

                    txb_1d_offset += context_ptr->blk_geom->tx_width[txb_itr] * context_ptr->blk_geom->tx_height[txb_itr];
                    if (context_ptr->blk_geom->has_uv)
                        txb_1d_offset_uv += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];

                    txb_itr++;
                } while (txb_itr < context_ptr->blk_geom->txb_count);

                //copy recon
                {
                    EbPictureBufferDesc          *ref_pic;
                    if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                    {
                        EbReferenceObject* refObj = (EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
                        ref_pic = refObj->reference_picture;
                    }
                    else
                        ref_pic = picture_control_set_ptr->recon_picture_ptr;
                    context_ptr->cu_origin_x = sb_origin_x + context_ptr->blk_geom->origin_x;
                    context_ptr->cu_origin_y = sb_origin_y + context_ptr->blk_geom->origin_y;

                    uint32_t  bwidth = context_ptr->blk_geom->bwidth;
                    uint32_t  bheight = context_ptr->blk_geom->bheight;

                    uint8_t* src_ptr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->buffer_y)[0]);
                    uint8_t* dst_ptr = ref_pic->buffer_y + ref_pic->origin_x + context_ptr->cu_origin_x + (ref_pic->origin_y + context_ptr->cu_origin_y)*ref_pic->stride_y;

                    uint32_t j;
                    for (j = 0; j < bheight; j++)
                        memcpy(dst_ptr + j * ref_pic->stride_y, src_ptr + j * 128, bwidth * sizeof(uint8_t));
                    if (context_ptr->blk_geom->has_uv)
                    {
                        bwidth = context_ptr->blk_geom->bwidth_uv;
                        bheight = context_ptr->blk_geom->bheight_uv;

                        src_ptr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->buffer_cb)[0]);

                        dst_ptr = ref_pic->buffer_cb + ref_pic->origin_x / 2 + ((context_ptr->cu_origin_x >> 3) << 3) / 2 + (ref_pic->origin_y / 2 + ((context_ptr->cu_origin_y >> 3) << 3) / 2)*ref_pic->stride_cb;

                        for (j = 0; j < bheight; j++)
                            memcpy(dst_ptr + j * ref_pic->stride_cb, src_ptr + j * 64, bwidth * sizeof(uint8_t));
                        src_ptr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->buffer_cr)[0]);

                        dst_ptr = ref_pic->buffer_cr + ref_pic->origin_x / 2 + ((context_ptr->cu_origin_x >> 3) << 3) / 2 + (ref_pic->origin_y / 2 + ((context_ptr->cu_origin_y >> 3) << 3) / 2)*ref_pic->stride_cr;

                        for (j = 0; j < bheight; j++)
                            memcpy(dst_ptr + j * ref_pic->stride_cr, src_ptr + j * 64, bwidth * sizeof(uint8_t));
                    }
                }
            }
            blk_it += ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][context_ptr->blk_geom->depth];
        }
        else
            blk_it += d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][context_ptr->blk_geom->depth];
    } // CU Loop

    return;
}
#endif
