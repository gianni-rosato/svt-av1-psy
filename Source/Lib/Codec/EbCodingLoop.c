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
#include "EbModeDecisionConfiguration.h"
#include "EbIntraPrediction.h"
#include "aom_dsp_rtcd.h"
#include "EbCodingLoop.h"

static const uint32_t me2Nx2NOffset[4] = { 0, 1, 5, 21 };
extern void av1_predict_intra_block(
#if TILES   
    TileInfo                    *tile,
#endif
#if INTRA_CORE_OPT
    ModeDecisionContext_t                  *md_context_ptr,
#endif
    STAGE                       stage,
    uint8_t                     intra_luma_left_mode,
    uint8_t                     intra_luma_top_mode,
    uint8_t                     intra_chroma_left_mode,
    uint8_t                     intra_chroma_top_mode,
    const BlockGeom            *blk_geom,
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
    EbPictureBufferDesc_t  *recon_buffer,
#if !INTRA_CORE_OPT
    int32_t col_off,
    int32_t row_off,
#endif
    int32_t plane,
    block_size bsize,
    uint32_t bl_org_x_pict,
    uint32_t bl_org_y_pict,
    uint32_t bl_org_x_mb,
    uint32_t bl_org_y_mb);

#if INTRA_10BIT_SUPPORT
void av1_predict_intra_block_16bit(
#if TILES   
    TileInfo               *tile,
#endif
    EncDecContext_t         *context_ptr,
    CodingUnit_t *cu_ptr,
    const Av1Common *cm,
    int32_t wpx,
    int32_t hpx,
    TxSize tx_size,
    PredictionMode mode,
    int32_t angle_delta,
    int32_t use_palette,
    FILTER_INTRA_MODE filter_intra_mode,
    uint16_t* topNeighArray,
    uint16_t* leftNeighArray,
    EbPictureBufferDesc_t  *recon_buffer,
    int32_t col_off,
    int32_t row_off,
    int32_t plane,
    block_size bsize,
    uint32_t bl_org_x_pict,
    uint32_t bl_org_y_pict);
#endif

/*******************************************
* set Penalize Skip Flag
*
* Summary: Set the PenalizeSkipFlag to true
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
#if ENCDEC_TX_SEARCH
    PictureControlSet_t    *picture_control_set_ptr,
#endif
    EncDecContext_t       *context_ptr,
    LargestCodingUnit_t   *sb_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    uint32_t                 cbQp,
    EbPictureBufferDesc_t *predSamples,             // no basis/offset
    EbPictureBufferDesc_t *coeffSamplesTB,          // lcu based
    EbPictureBufferDesc_t *residual16bit,           // no basis/offset
    EbPictureBufferDesc_t *transform16bit,          // no basis/offset
    EbPictureBufferDesc_t *inverse_quant_buffer,
    int16_t                *transformScratchBuffer,
    EbAsm                 asm_type,
    uint32_t                  *count_non_zero_coeffs,
    uint32_t                 component_mask,
    uint32_t                   use_delta_qp,
    uint32_t                 dZoffset,
    uint16_t                 *eob,
    MacroblockPlane       *candidate_plane);


typedef void(*EB_AV1_GENERATE_RECON_FUNC_PTR)(
    EncDecContext_t       *context_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    EbPictureBufferDesc_t *predSamples,     // no basis/offset
    EbPictureBufferDesc_t *residual16bit,    // no basis/offset
    int16_t                *transformScratchBuffer,
    uint32_t                 component_mask,
    uint16_t                *eob,
    EbAsm                 asm_type);


typedef void(*EB_GENERATE_RECON_FUNC_PTR)(
    EncDecContext_t       *context_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    EbPictureBufferDesc_t *predSamples,     // no basis/offset
    EbPictureBufferDesc_t *residual16bit,    // no basis/offset
    int16_t                *transformScratchBuffer,
    EbAsm                 asm_type);

typedef void(*EB_GENERATE_RECON_INTRA_4x4_FUNC_PTR)(
    EncDecContext_t       *context_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    EbPictureBufferDesc_t *predSamples,     // no basis/offset
    EbPictureBufferDesc_t *residual16bit,    // no basis/offset
    int16_t                *transformScratchBuffer,
    uint32_t                 component_mask,
    EbAsm                 asm_type);

typedef EbErrorType(*EB_GENERATE_INTRA_SAMPLES_FUNC_PTR)(
    EbBool                         *is_left_availble,
    EbBool                         *is_above_availble,
    EbBool                     constrained_intra_flag,   //input parameter, indicates if constrained intra is switched on/off
    EbBool                     strongIntraSmoothingFlag,
    uint32_t                      origin_x,
    uint32_t                      origin_y,
    uint32_t                      size,
    uint32_t                      cu_depth,
    NeighborArrayUnit_t        *mode_type_neighbor_array,
    NeighborArrayUnit_t        *luma_recon_neighbor_array,
    NeighborArrayUnit_t        *cb_recon_neighbor_array,
    NeighborArrayUnit_t        *cr_recon_neighbor_array,
    void                       *refWrapperPtr,
    EbBool                     pictureLeftBoundary,
    EbBool                     pictureTopBoundary,
    EbBool                     pictureRightBoundary);
typedef EbErrorType(*EB_ENC_PASS_INTRA_FUNC_PTR)(
    uint8_t                          upsample_left,
    uint8_t                          upsample_above,
    uint8_t                          upsample_left_chroma,
    uint8_t                          upsample_above_chroma,
    EbBool                         is_left_availble,
    EbBool                         is_above_availble,
    void                       *referenceSamples,
    uint32_t                      origin_x,
    uint32_t                      origin_y,
    uint32_t                      puSize,
    EbPictureBufferDesc_t      *prediction_ptr,
    uint32_t                      luma_mode,
    uint32_t                      chroma_mode,
    int32_t                      angle_delta,
    uint16_t                      bitdepth,
    EbAsm                      asm_type);


/***************************************************
* Update Intra Mode Neighbor Arrays
***************************************************/
static void EncodePassUpdateIntraModeNeighborArrays(
    NeighborArrayUnit_t     *mode_type_neighbor_array,
    NeighborArrayUnit_t     *intra_luma_mode_neighbor_array,
    NeighborArrayUnit_t     *intra_chroma_mode_neighbor_array,
    uint8_t                    luma_mode,
    uint8_t                    chroma_mode,
    uint32_t                   origin_x,
    uint32_t                   origin_y,
    uint32_t                   width,
    uint32_t                   height,
    uint32_t                   width_uv,
    uint32_t                   height_uv,
    uint32_t                   component_mask)
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

    return;
}

/***************************************************
* Update Inter Mode Neighbor Arrays
***************************************************/
static void EncodePassUpdateInterModeNeighborArrays(
    NeighborArrayUnit_t     *mode_type_neighbor_array,
    NeighborArrayUnit_t     *mv_neighbor_array,
    NeighborArrayUnit_t     *skipNeighborArray,
    MvUnit_t                *mv_unit,
    uint8_t                   *skip_flag,
    uint32_t                   origin_x,
    uint32_t                   origin_y,
    uint32_t                   bwidth,
    uint32_t                   bheight)
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

    return;
}

/***************************************************
* Update Recon Samples Neighbor Arrays
***************************************************/
static void EncodePassUpdateReconSampleNeighborArrays(
    NeighborArrayUnit_t     *lumaReconSampleNeighborArray,
    NeighborArrayUnit_t     *cbReconSampleNeighborArray,
    NeighborArrayUnit_t     *crReconSampleNeighborArray,
    EbPictureBufferDesc_t   *recon_buffer,
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

        if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK)
        {
            // Recon Samples - Cb
            neighbor_array_unit16bit_sample_write(
                cbReconSampleNeighborArray,
                (uint16_t*)(recon_buffer->bufferCb),
                recon_buffer->strideCb,
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
                (uint16_t*)(recon_buffer->bufferCr),
                recon_buffer->strideCr,
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

        if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK)
        {
            // Recon Samples - Cb
            neighbor_array_unit_sample_write(
                cbReconSampleNeighborArray,
                recon_buffer->bufferCb,
                recon_buffer->strideCb,
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
                recon_buffer->bufferCr,
                recon_buffer->strideCr,
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
    CodingUnit_t            *cu_ptr,
    uint32_t                   pu_origin_x,
    uint32_t                   pu_origin_y,
    uint32_t                   sb_sz,
    NeighborArrayUnit_t     *intraLumaNeighborArray,
    NeighborArrayUnit_t     *intraChromaNeighborArray,
    NeighborArrayUnit_t     *mode_type_neighbor_array)
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
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != INTRA_MODE) ? DC_PRED/*EB_INTRA_DC*/ :
        (uint32_t)intraLumaNeighborArray->leftArray[intraLumaModeLeftNeighborIndex]);

    (&cu_ptr->prediction_unit_array[0])->intra_luma_top_mode = (uint32_t)(
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != INTRA_MODE) ? DC_PRED/*EB_INTRA_DC*/ :
        (uint32_t)intraLumaNeighborArray->topArray[intraLumaModeTopNeighborIndex]);       //   use DC. This seems like we could use a LCU-width

    uint32_t modeTypeLeftNeighborIndex_round = get_neighbor_array_unit_left_index(
        mode_type_neighbor_array,
        puOriginY_round);
    uint32_t modeTypeTopNeighborIndex_round = get_neighbor_array_unit_top_index(
        mode_type_neighbor_array,
        puOriginX_round);

    (&cu_ptr->prediction_unit_array[0])->intra_chroma_left_mode = (uint32_t)(
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex_round] != INTRA_MODE) ? UV_DC_PRED :
        (uint32_t)intraChromaNeighborArray->leftArray[intraChromaModeLeftNeighborIndex]);

    (&cu_ptr->prediction_unit_array[0])->intra_chroma_top_mode = (uint32_t)(
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex_round] != INTRA_MODE) ? UV_DC_PRED :
        (uint32_t)intraChromaNeighborArray->topArray[intraChromaModeTopNeighborIndex]);       //   use DC. This seems like we could use a LCU-width


    return;
}


void PfZeroOutUselessQuadrants(
    int16_t* transformCoeffBuffer,
    uint32_t  transformCoeffStride,
    uint32_t  quadrantSize,
    EbAsm  asm_type) {

    pic_zero_out_coef_func_ptr_array[asm_type][quadrantSize >> 3](
        transformCoeffBuffer,
        transformCoeffStride,
        quadrantSize,
        quadrantSize,
        quadrantSize);

    pic_zero_out_coef_func_ptr_array[asm_type][quadrantSize >> 3](
        transformCoeffBuffer,
        transformCoeffStride,
        quadrantSize * transformCoeffStride,
        quadrantSize,
        quadrantSize);

    pic_zero_out_coef_func_ptr_array[asm_type][quadrantSize >> 3](
        transformCoeffBuffer,
        transformCoeffStride,
        quadrantSize * transformCoeffStride + quadrantSize,
        quadrantSize,
        quadrantSize);

}

void encode_pass_tx_search(
    PictureControlSet_t            *picture_control_set_ptr,
    EncDecContext_t                *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                       cbQp,
    EbPictureBufferDesc_t          *coeffSamplesTB,
    EbPictureBufferDesc_t          *residual16bit,
    EbPictureBufferDesc_t          *transform16bit,
    EbPictureBufferDesc_t          *inverse_quant_buffer,
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
#if ENCDEC_TX_SEARCH
    PictureControlSet_t    *picture_control_set_ptr,
#endif
    EncDecContext_t       *context_ptr,
    LargestCodingUnit_t   *sb_ptr,
    uint32_t                 origin_x,   //pic based tx org x
    uint32_t                 origin_y,   //pic based tx org y
    uint32_t                 cbQp,
    EbPictureBufferDesc_t *predSamples,             // no basis/offset
    EbPictureBufferDesc_t *coeffSamplesTB,          // lcu based
    EbPictureBufferDesc_t *residual16bit,           // no basis/offset
    EbPictureBufferDesc_t *transform16bit,          // no basis/offset
    EbPictureBufferDesc_t *inverse_quant_buffer,
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
    (void)cbQp;

    //    uint32_t                 chroma_qp = cbQp;
    CodingUnit_t          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    //    EB_SLICE               slice_type = sb_ptr->picture_control_set_ptr->slice_type;
    //    uint32_t                 temporal_layer_index = sb_ptr->picture_control_set_ptr->temporal_layer_index;
    uint32_t                 qp = cu_ptr->qp;
    EbPictureBufferDesc_t  *input_samples = context_ptr->input_samples;

    uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
    uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4

    const uint32_t inputLumaOffset = ((origin_y + input_samples->origin_y)          * input_samples->stride_y) + (origin_x + input_samples->origin_x);
    const uint32_t inputCbOffset = (((round_origin_y + input_samples->origin_y) >> 1)    * input_samples->strideCb) + ((round_origin_x + input_samples->origin_x) >> 1);
    const uint32_t inputCrOffset = (((round_origin_y + input_samples->origin_y) >> 1)    * input_samples->strideCr) + ((round_origin_x + input_samples->origin_x) >> 1);
    const uint32_t predLumaOffset = ((predSamples->origin_y + origin_y)        * predSamples->stride_y) + (predSamples->origin_x + origin_x);
    const uint32_t predCbOffset = (((predSamples->origin_y + round_origin_y) >> 1)  * predSamples->strideCb) + ((predSamples->origin_x + round_origin_x) >> 1);
    const uint32_t predCrOffset = (((predSamples->origin_y + round_origin_y) >> 1)  * predSamples->strideCr) + ((predSamples->origin_x + round_origin_x) >> 1);


    const uint32_t scratchLumaOffset = context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr] * SB_STRIDE_Y;
    const uint32_t scratchCbOffset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t scratchCrOffset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;


    const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

    const uint32_t coeff1dOffsetChroma = context_ptr->coded_area_sb_uv;
    UNUSED(coeff1dOffsetChroma);


    //uint8_t enable_contouring_qc_update_flag;
    //enable_contouring_qc_update_flag = DeriveContouringClass(
    //    sb_ptr->picture_control_set_ptr->parent_pcs_ptr,
    //    sb_ptr->index,
    //    cu_ptr->leaf_index) && (cu_ptr->qp < sb_ptr->picture_control_set_ptr->picture_qp);

    EbBool clean_sparse_coeff_flag = EB_FALSE;

    context_ptr->three_quad_energy = 0;
    //**********************************
    // Luma
    //**********************************
    if (component_mask == PICTURE_BUFFER_DESC_FULL_MASK || component_mask == PICTURE_BUFFER_DESC_LUMA_MASK)
    {

        ResidualKernel(
            input_samples->buffer_y + inputLumaOffset,
            input_samples->stride_y,
            predSamples->buffer_y + predLumaOffset,
            predSamples->stride_y,
            ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
            residual16bit->stride_y,
            context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);
        
        uint8_t  tx_search_skip_fag = picture_control_set_ptr->parent_pcs_ptr->tx_search_level == TX_SEARCH_ENC_DEC ? get_skip_tx_search_flag(
            context_ptr->blk_geom->sq_size,
            MAX_MODE_COST,
            0,
            1) : 1;

        if (!tx_search_skip_fag) {

                encode_pass_tx_search(
                    picture_control_set_ptr,
                    context_ptr,
                    sb_ptr,
                    cbQp,
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
            ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
            residual16bit->stride_y,
            ((tran_low_t*)transform16bit->buffer_y) + coeff1dOffset,
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_Y],
            asm_type,
            PLANE_TYPE_Y,
            context_ptr->trans_coeff_shape_luma);

        av1_quantize_inv_quantize(
            sb_ptr->picture_control_set_ptr,
            ((tran_low_t*)transform16bit->buffer_y) + coeff1dOffset,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->buffer_y) + coeff1dOffset,
            ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
            qp,
            context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize[context_ptr->txb_itr],
            &eob[0],
            candidate_plane[0],
            asm_type,
            &(count_non_zero_coeffs[0]),
            0,
            0,
            COMPONENT_LUMA,
#if QT_10BIT_SUPPORT
            BIT_INCREMENT_8BIT,
#endif
            txb_ptr->transform_type[PLANE_TYPE_Y],
            clean_sparse_coeff_flag);

        txb_ptr->y_has_coeff = count_non_zero_coeffs[0] ? EB_TRUE : EB_FALSE;


#if TX_TYPE_FIX
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
#endif

#if CHROMA_BLIND
        if (cu_ptr->prediction_mode_flag == INTRA_MODE && (context_ptr->evaluate_cfl_ep || cu_ptr->prediction_unit_array->intra_chroma_mode == UV_CFL_PRED)) {
#else
        if (cu_ptr->prediction_mode_flag == INTRA_MODE && cu_ptr->prediction_unit_array->intra_chroma_mode == UV_CFL_PRED) {
#endif


            EbPictureBufferDesc_t *reconSamples = predSamples;
            uint32_t reconLumaOffset = (reconSamples->origin_y + origin_y)            * reconSamples->stride_y + (reconSamples->origin_x + origin_x);

            if (txb_ptr->y_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

                uint8_t     *predBuffer = predSamples->buffer_y + predLumaOffset;


                av1_inv_transform_recon8bit(
                    ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
                    predBuffer,
                    predSamples->stride_y,
                    context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                    txb_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    eob[0]);

            }

            // Down sample Luma
            cfl_luma_subsampling_420_lbd_c(
                reconSamples->buffer_y + reconLumaOffset,
                reconSamples->stride_y,
#if CHROMA_BLIND
                context_ptr->md_context->pred_buf_q3,
#else
                context_ptr->pred_buf_q3,
#endif
                context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);


            int32_t round_offset = ((context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr])*(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr])) / 2;



            subtract_average(
#if CHROMA_BLIND
                context_ptr->md_context->pred_buf_q3,
#else
                context_ptr->pred_buf_q3,
#endif
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                round_offset,
                LOG2F(context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr]) + LOG2F(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]));

#if CHROMA_BLIND
            if (context_ptr->evaluate_cfl_ep)
            {
                // 3: Loop over alphas and find the best or choose DC
                // Use the 1st spot of the candidate buffer to hold cfl settings: (1) to use same kernel as MD for CFL evaluation: cfl_rd_pick_alpha() (toward unification), (2) to avoid dedicated buffers for CFL evaluation @ EP (toward less memory)
                ModeDecisionCandidateBuffer_t  *candidateBuffer = &(context_ptr->md_context->candidate_buffer_ptr_array[0][0]);

                // Input(s)
                candidateBuffer->candidate_ptr->type = INTRA_MODE;
                candidateBuffer->candidate_ptr->intra_luma_mode = cu_ptr->pred_mode;
                candidateBuffer->candidate_ptr->cfl_alpha_signs = 0;
                candidateBuffer->candidate_ptr->cfl_alpha_idx = 0;
                context_ptr->md_context->blk_geom = context_ptr->blk_geom;

                EbByte src_pred_ptr;
                EbByte dst_pred_ptr;

                // Copy Cb pred samples from ep buffer to md buffer
                src_pred_ptr = predSamples->bufferCb + predCbOffset;
                dst_pred_ptr = &(candidateBuffer->prediction_ptr->bufferCb[scratchCbOffset]);
                for (int i = 0; i < context_ptr->blk_geom->bheight_uv; i++) {
                    memcpy(dst_pred_ptr, src_pred_ptr, context_ptr->blk_geom->bwidth_uv);
                    src_pred_ptr += predSamples->strideCb;
                    dst_pred_ptr += candidateBuffer->prediction_ptr->strideCb;
                }

                // Copy Cr pred samples from ep buffer to md buffer
                src_pred_ptr = predSamples->bufferCr + predCrOffset;
                dst_pred_ptr = &(candidateBuffer->prediction_ptr->bufferCr[scratchCrOffset]);
                for (int i = 0; i < context_ptr->blk_geom->bheight_uv; i++) {
                    memcpy(dst_pred_ptr, src_pred_ptr, context_ptr->blk_geom->bwidth_uv);
                    src_pred_ptr += predSamples->strideCr;
                    dst_pred_ptr += candidateBuffer->prediction_ptr->strideCr;
                }

                cfl_rd_pick_alpha(
                    picture_control_set_ptr,
                    candidateBuffer,
                    sb_ptr,
                    context_ptr->md_context,
                    input_samples,
                    inputCbOffset,
                    scratchCbOffset,
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
#endif
                int32_t alpha_q3 =
                    cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_U); // once for U, once for V

                //TOCHANGE
                //assert(chromaSize * CFL_BUF_LINE + chromaSize <= CFL_BUF_SQUARE);

                cfl_predict_lbd(
#if CHROMA_BLIND
                    context_ptr->md_context->pred_buf_q3,
#else
                    context_ptr->pred_buf_q3,
#endif
                    predSamples->bufferCb + predCbOffset,
                    predSamples->strideCb,
                    predSamples->bufferCb + predCbOffset,
                    predSamples->strideCb,
                    alpha_q3,
                    8,
                    context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
                alpha_q3 =
                    cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_V); // once for U, once for V

                //TOCHANGE
                //assert(chromaSize * CFL_BUF_LINE + chromaSize <= CFL_BUF_SQUARE);

                cfl_predict_lbd(
#if CHROMA_BLIND
                    context_ptr->md_context->pred_buf_q3,
#else
                    context_ptr->pred_buf_q3,
#endif
                    predSamples->bufferCr + predCrOffset,
                    predSamples->strideCr,
                    predSamples->bufferCr + predCrOffset,
                    predSamples->strideCr,
                    alpha_q3,
                    8,
                    context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);

#if CHROMA_BLIND
            }
#endif
        }

    }

    if (component_mask == PICTURE_BUFFER_DESC_FULL_MASK || component_mask == PICTURE_BUFFER_DESC_CHROMA_MASK)
    {
        //**********************************
        // Cb
        //**********************************

        ResidualKernel(
            input_samples->bufferCb + inputCbOffset,
            input_samples->strideCb,
            predSamples->bufferCb + predCbOffset,
            predSamples->strideCb,
            ((int16_t*)residual16bit->bufferCb) + scratchCbOffset,
            residual16bit->strideCb,
            context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);


        ResidualKernel(
            input_samples->bufferCr + inputCrOffset,
            input_samples->strideCr,
            predSamples->bufferCr + predCrOffset,
            predSamples->strideCr,
            ((int16_t*)residual16bit->bufferCr) + scratchCrOffset,
            residual16bit->strideCr,

            context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);

        av1_estimate_transform(
            ((int16_t*)residual16bit->bufferCb) + scratchCbOffset,
            residual16bit->strideCb,
            ((tran_low_t*)transform16bit->bufferCb) + context_ptr->coded_area_sb_uv,
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_UV],
            asm_type,
            PLANE_TYPE_UV,
            context_ptr->trans_coeff_shape_chroma);


        av1_quantize_inv_quantize(
            sb_ptr->picture_control_set_ptr,

            ((tran_low_t*)transform16bit->bufferCb) + context_ptr->coded_area_sb_uv,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->bufferCb) + context_ptr->coded_area_sb_uv,
            ((int32_t*)inverse_quant_buffer->bufferCb) + context_ptr->coded_area_sb_uv,
            qp,
            context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
            &eob[1],
            candidate_plane[1],
            asm_type,
            &(count_non_zero_coeffs[1]),
            0,
            0,
            COMPONENT_CHROMA_CB,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_UV],
            clean_sparse_coeff_flag);


        txb_ptr->u_has_coeff = count_non_zero_coeffs[1] ? EB_TRUE : EB_FALSE;

        //**********************************
        // Cr
        //**********************************

        av1_estimate_transform(
            ((int16_t*)residual16bit->bufferCr) + scratchCbOffset,
            residual16bit->strideCr,
            ((tran_low_t*)transform16bit->bufferCr) + context_ptr->coded_area_sb_uv,
            NOT_USED_VALUE,
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
            &context_ptr->three_quad_energy,
            transformScratchBuffer,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_UV],
            asm_type,
            PLANE_TYPE_UV,
            context_ptr->trans_coeff_shape_chroma);


        av1_quantize_inv_quantize(
            sb_ptr->picture_control_set_ptr,
            ((tran_low_t*)transform16bit->bufferCr) + context_ptr->coded_area_sb_uv,
            NOT_USED_VALUE,
            ((int32_t*)coeffSamplesTB->bufferCr) + context_ptr->coded_area_sb_uv,
            ((tran_low_t*)inverse_quant_buffer->bufferCr) + context_ptr->coded_area_sb_uv,
            qp,
            context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
            context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
            &eob[2],
            candidate_plane[2],
            asm_type,
            &(count_non_zero_coeffs[2]),
            0,
            0,
            COMPONENT_CHROMA_CR,
            BIT_INCREMENT_8BIT,
            txb_ptr->transform_type[PLANE_TYPE_UV],
            clean_sparse_coeff_flag);

        txb_ptr->v_has_coeff = count_non_zero_coeffs[2] ? EB_TRUE : EB_FALSE;
    }

    txb_ptr->trans_coeff_shape_luma = context_ptr->trans_coeff_shape_luma;
    txb_ptr->trans_coeff_shape_chroma = context_ptr->trans_coeff_shape_chroma;
    txb_ptr->nz_coef_count[0] = (uint16_t)count_non_zero_coeffs[0];
    txb_ptr->nz_coef_count[1] = (uint16_t)count_non_zero_coeffs[1];
    txb_ptr->nz_coef_count[2] = (uint16_t)count_non_zero_coeffs[2];
    return;
}

void encode_pass_tx_search_hbd(
    PictureControlSet_t            *picture_control_set_ptr,
    EncDecContext_t                *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                       cbQp,
    EbPictureBufferDesc_t          *coeffSamplesTB,
    EbPictureBufferDesc_t          *residual16bit,
    EbPictureBufferDesc_t          *transform16bit,
    EbPictureBufferDesc_t          *inverse_quant_buffer,
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
#if ENCDEC_TX_SEARCH
    PictureControlSet_t    *picture_control_set_ptr,
#endif
    EncDecContext_t       *context_ptr,
    LargestCodingUnit_t   *sb_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    uint32_t                 cbQp,
    EbPictureBufferDesc_t *predSamples,         // no basis/offset
    EbPictureBufferDesc_t *coeffSamplesTB,      // lcu based
    EbPictureBufferDesc_t *residual16bit,       // no basis/offset
    EbPictureBufferDesc_t *transform16bit,      // no basis/offset
    EbPictureBufferDesc_t *inverse_quant_buffer,
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
#if QT_10BIT_SUPPORT
    (void)cbQp;
#else
    uint32_t                 chroma_qp = cbQp;
#endif

    CodingUnit_t          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    //    EB_SLICE               slice_type = sb_ptr->picture_control_set_ptr->slice_type;
    //    uint32_t                 temporal_layer_index = sb_ptr->picture_control_set_ptr->temporal_layer_index;
    uint32_t                 qp = cu_ptr->qp;

    EbPictureBufferDesc_t *inputSamples16bit = context_ptr->input_sample16bit_buffer;
    EbPictureBufferDesc_t *predSamples16bit = predSamples;
    uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
    uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4
    const uint32_t           inputLumaOffset = context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr] * SB_STRIDE_Y;
    const uint32_t           inputCbOffset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t           inputCrOffset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
    const uint32_t           predLumaOffset = ((predSamples16bit->origin_y + origin_y)        * predSamples16bit->stride_y) + (predSamples16bit->origin_x + origin_x);
    const uint32_t           predCbOffset = (((predSamples16bit->origin_y + round_origin_y) >> 1)  * predSamples16bit->strideCb) + ((predSamples16bit->origin_x + round_origin_x) >> 1);
    const uint32_t           predCrOffset = (((predSamples16bit->origin_y + round_origin_y) >> 1)  * predSamples16bit->strideCr) + ((predSamples16bit->origin_x + round_origin_x) >> 1);
    const uint32_t scratchLumaOffset = context_ptr->blk_geom->origin_x + context_ptr->blk_geom->origin_y * SB_STRIDE_Y;
    const uint32_t scratchCbOffset = ROUND_UV(context_ptr->blk_geom->origin_x) / 2 + ROUND_UV(context_ptr->blk_geom->origin_y) / 2 * SB_STRIDE_UV;
    const uint32_t scratchCrOffset = ROUND_UV(context_ptr->blk_geom->origin_x) / 2 + ROUND_UV(context_ptr->blk_geom->origin_y) / 2 * SB_STRIDE_UV;

#if QT_10BIT_SUPPORT
    const uint32_t coeff1dOffset = context_ptr->coded_area_sb;
    const uint32_t coeff1dOffsetChroma = context_ptr->coded_area_sb_uv;
    UNUSED(coeff1dOffsetChroma);
#endif

    EbBool clean_sparse_coeff_flag = EB_FALSE;

    //Update QP for Quant
    qp += QP_BD_OFFSET;
#if !QT_10BIT_SUPPORT
    chroma_qp += QP_BD_OFFSET;
#endif

    {

        //**********************************
        // Luma
        //**********************************
        if (component_mask == PICTURE_BUFFER_DESC_FULL_MASK || component_mask == PICTURE_BUFFER_DESC_LUMA_MASK) {

            residual_kernel16bit(
                ((uint16_t*)inputSamples16bit->buffer_y) + inputLumaOffset,
                inputSamples16bit->stride_y,
                ((uint16_t*)predSamples16bit->buffer_y) + predLumaOffset,
                predSamples16bit->stride_y,
                ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
                residual16bit->stride_y,
                context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);

            uint8_t  tx_search_skip_fag = picture_control_set_ptr->parent_pcs_ptr->tx_search_level == TX_SEARCH_ENC_DEC ? get_skip_tx_search_flag(
                context_ptr->blk_geom->sq_size,
                MAX_MODE_COST,
                0,
                1) : 1;

            if (!tx_search_skip_fag) {

                    encode_pass_tx_search_hbd(
                        picture_control_set_ptr,
                        context_ptr,
                        sb_ptr,
                        cbQp,
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
                ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
                residual16bit->stride_y,
                ((tran_low_t*)transform16bit->buffer_y) + coeff1dOffset,
                NOT_USED_VALUE,
                context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                &context_ptr->three_quad_energy,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_Y],
                asm_type,
                PLANE_TYPE_Y,
                context_ptr->trans_coeff_shape_luma);

            av1_quantize_inv_quantize(
                sb_ptr->picture_control_set_ptr,
                ((int32_t*)transform16bit->buffer_y) + coeff1dOffset,
                NOT_USED_VALUE,
#if QT_10BIT_SUPPORT
                ((int32_t*)coeffSamplesTB->buffer_y) + coeff1dOffset,
#else
                ((int32_t*)coeffSamplesTB->buffer_y) + scratchLumaOffset,
#endif
                ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
                qp,
                context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height[context_ptr->txb_itr],
                context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                &eob[0],
                candidate_plane[0],
                asm_type,
                &(count_non_zero_coeffs[0]),
                0,
                0,
                COMPONENT_LUMA,
#if QT_10BIT_SUPPORT
                BIT_INCREMENT_10BIT,
#endif
                txb_ptr->transform_type[PLANE_TYPE_Y],
                clean_sparse_coeff_flag);
            txb_ptr->y_has_coeff = count_non_zero_coeffs[0] ? EB_TRUE : EB_FALSE;
#if TX_TYPE_FIX
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
#endif

        }

        if (component_mask == PICTURE_BUFFER_DESC_FULL_MASK || component_mask == PICTURE_BUFFER_DESC_CHROMA_MASK) {

            if (cu_ptr->prediction_mode_flag == INTRA_MODE && cu_ptr->prediction_unit_array->intra_chroma_mode == UV_CFL_PRED) {
                EbPictureBufferDesc_t *reconSamples = predSamples16bit;
                uint32_t reconLumaOffset = (reconSamples->origin_y + origin_y)            * reconSamples->stride_y + (reconSamples->origin_x + origin_x);
                if (txb_ptr->y_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
#if QT_10BIT_SUPPORT
                    uint16_t     *predBuffer = ((uint16_t*)predSamples16bit->buffer_y) + predLumaOffset;
                    av1_inv_transform_recon(
                        ((int32_t*)inverse_quant_buffer->buffer_y) + coeff1dOffset,
                        CONVERT_TO_BYTEPTR(predBuffer),
                        predSamples->stride_y,
                        context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                        BIT_INCREMENT_10BIT,
                        txb_ptr->transform_type[PLANE_TYPE_Y],
                        PLANE_TYPE_Y,
                        eob[0]);
#else
                    av1_estimate_inv_transform(
                        ((int32_t*)inverse_quant_buffer->buffer_y) + scratchLumaOffset,
                        64,
                        ((int32_t*)inverse_quant_buffer->buffer_y) + scratchLumaOffset,
                        64,
                        txb_size,
                        transformScratchBuffer,
                        BIT_INCREMENT_10BIT,
                        txb_ptr->transform_type[PLANE_TYPE_Y],
                        eob[0],
                        asm_type,
                        0);

                    picture_addition_kernel16_bit(
                        ((uint16_t*)predSamples16bit->buffer_y) + predLumaOffset,
                        predSamples16bit->stride_y,
                        ((int32_t*)inverse_quant_buffer->buffer_y) + scratchLumaOffset,
                        64,
                        ((uint16_t*)reconSamples->buffer_y) + reconLumaOffset,
                        reconSamples->stride_y,
                        txb_size,
                        txb_size,
                        10);
#endif
                }

                // Down sample Luma
                cfl_luma_subsampling_420_hbd_c(
                    ((uint16_t*)reconSamples->buffer_y) + reconLumaOffset,
                    reconSamples->stride_y,
#if CHROMA_BLIND
                    context_ptr->md_context->pred_buf_q3,
#else
                    context_ptr->pred_buf_q3,
#endif
                    context_ptr->blk_geom->tx_width[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height[context_ptr->txb_itr]);

                int32_t round_offset = ((context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr])*(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr])) / 2;


                subtract_average(
#if CHROMA_BLIND
                    context_ptr->md_context->pred_buf_q3,
#else
                    context_ptr->pred_buf_q3,
#endif
                    context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                    round_offset,
                    LOG2F(context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr]) + LOG2F(context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]));



                int32_t alpha_q3 =
                    cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_U); // once for U, once for V
                // TOCHANGE
                // assert(chromaSize * CFL_BUF_LINE + chromaSize <=                CFL_BUF_SQUARE);

                cfl_predict_hbd(
#if CHROMA_BLIND
                    context_ptr->md_context->pred_buf_q3,
#else
                    context_ptr->pred_buf_q3,
#endif
                    ((uint16_t*)predSamples16bit->bufferCb) + predCbOffset,
                    predSamples16bit->strideCb,
                    ((uint16_t*)predSamples16bit->bufferCb) + predCbOffset,
                    predSamples16bit->strideCb,
                    alpha_q3,
                    10,
                    context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);

                alpha_q3 =
                    cfl_idx_to_alpha(cu_ptr->prediction_unit_array->cfl_alpha_idx, cu_ptr->prediction_unit_array->cfl_alpha_signs, CFL_PRED_V); // once for U, once for V
                // TOCHANGE
                //assert(chromaSize * CFL_BUF_LINE + chromaSize <=                CFL_BUF_SQUARE);

                cfl_predict_hbd(
#if CHROMA_BLIND
                    context_ptr->md_context->pred_buf_q3,
#else
                    context_ptr->pred_buf_q3,
#endif
                    ((uint16_t*)predSamples16bit->bufferCr) + predCrOffset,
                    predSamples16bit->strideCr,
                    ((uint16_t*)predSamples16bit->bufferCr) + predCrOffset,
                    predSamples16bit->strideCr,
                    alpha_q3,
                    10,
                    context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                    context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);
            }

        }

        {
            //**********************************
            // Cb
            //**********************************
            residual_kernel16bit(
                ((uint16_t*)inputSamples16bit->bufferCb) + inputCbOffset,
                inputSamples16bit->strideCb,
                ((uint16_t*)predSamples16bit->bufferCb) + predCbOffset,
                predSamples16bit->strideCb,
                ((int16_t*)residual16bit->bufferCb) + scratchCbOffset,
                residual16bit->strideCb,
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);


            residual_kernel16bit(

                ((uint16_t*)inputSamples16bit->bufferCr) + inputCrOffset,
                inputSamples16bit->strideCr,
                ((uint16_t*)predSamples16bit->bufferCr) + predCrOffset,
                predSamples16bit->strideCr,
                ((int16_t*)residual16bit->bufferCr) + scratchCrOffset,
                residual16bit->strideCr,
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr]);


            av1_estimate_transform(
                ((int16_t*)residual16bit->bufferCb) + scratchCbOffset,
                residual16bit->strideCb,

                ((tran_low_t*)transform16bit->bufferCb) + context_ptr->coded_area_sb_uv,
                NOT_USED_VALUE,
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
                &context_ptr->three_quad_energy,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                asm_type,
                PLANE_TYPE_UV,
                context_ptr->trans_coeff_shape_chroma);


            av1_quantize_inv_quantize(
                sb_ptr->picture_control_set_ptr,
                ((int32_t*)transform16bit->bufferCb) + context_ptr->coded_area_sb_uv,
                NOT_USED_VALUE,

#if QT_10BIT_SUPPORT
                ((int32_t*)coeffSamplesTB->bufferCb) + context_ptr->coded_area_sb_uv,
#else
                ((int32_t*)coeffSamplesTB->bufferCb) + scratchCbOffset,
#endif
                ((int32_t*)inverse_quant_buffer->bufferCb) + context_ptr->coded_area_sb_uv,
                qp,
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
                &eob[1],
                candidate_plane[1],
                asm_type,
                &(count_non_zero_coeffs[1]),
                0,
                0,
                COMPONENT_CHROMA_CB,
#if QT_10BIT_SUPPORT
                BIT_INCREMENT_10BIT,
#endif
                txb_ptr->transform_type[PLANE_TYPE_UV],
                clean_sparse_coeff_flag);

            txb_ptr->u_has_coeff = count_non_zero_coeffs[1] ? EB_TRUE : EB_FALSE;

            //**********************************
            // Cr
            //**********************************
#if !QT_10BIT_SUPPORT
            encode_transform(
                ((int16_t*)residual16bit->bufferCr) + scratchCrOffset,
                32,
                ((int16_t*)transform16bit->bufferCr) + scratchCrOffset,
                32,
                txb_size >> 1,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                EB_FALSE,
                context_ptr->trans_coeff_shape_chroma,
                asm_type);
#endif

            av1_estimate_transform(
                ((int16_t*)residual16bit->bufferCr) + scratchCbOffset,

                residual16bit->strideCr,

                ((tran_low_t*)transform16bit->bufferCr) + context_ptr->coded_area_sb_uv,
                NOT_USED_VALUE,


                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
                &context_ptr->three_quad_energy,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                asm_type,
                PLANE_TYPE_UV,
                context_ptr->trans_coeff_shape_chroma);


            av1_quantize_inv_quantize(
                sb_ptr->picture_control_set_ptr,
                ((int32_t*)transform16bit->bufferCr) + context_ptr->coded_area_sb_uv,
                NOT_USED_VALUE,
#if QT_10BIT_SUPPORT
                ((int32_t*)coeffSamplesTB->bufferCr) + context_ptr->coded_area_sb_uv,
#else
                ((int32_t*)coeffSamplesTB->bufferCr) + scratchCbOffset,
#endif
                ((int32_t*)inverse_quant_buffer->bufferCr) + context_ptr->coded_area_sb_uv,
                qp,
                context_ptr->blk_geom->tx_width_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->tx_height_uv[context_ptr->txb_itr],
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
                &eob[2],
                candidate_plane[2],
                asm_type,
                &(count_non_zero_coeffs[2]),
                0,
                0,
                COMPONENT_CHROMA_CR,
#if QT_10BIT_SUPPORT
                BIT_INCREMENT_10BIT,
#endif
                txb_ptr->transform_type[PLANE_TYPE_UV],
                clean_sparse_coeff_flag);
            txb_ptr->v_has_coeff = count_non_zero_coeffs[2] ? EB_TRUE : EB_FALSE;

        }
    }


    txb_ptr->trans_coeff_shape_luma = context_ptr->trans_coeff_shape_luma;
    txb_ptr->trans_coeff_shape_chroma = context_ptr->trans_coeff_shape_chroma;
    txb_ptr->nz_coef_count[0] = (uint16_t)count_non_zero_coeffs[0];
    txb_ptr->nz_coef_count[1] = (uint16_t)count_non_zero_coeffs[1];
    txb_ptr->nz_coef_count[2] = (uint16_t)count_non_zero_coeffs[2];
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
    EncDecContext_t       *context_ptr,
    uint32_t               origin_x,
    uint32_t               origin_y,
    EbPictureBufferDesc_t *predSamples,     // no basis/offset
    EbPictureBufferDesc_t *residual16bit,    // no basis/offset
    int16_t               *transformScratchBuffer,
    uint32_t               component_mask,
    uint16_t              *eob,
    EbAsm                  asm_type)
{
    uint32_t               predLumaOffset;
    uint32_t               predChromaOffset;
    CodingUnit_t          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];

    // *Note - The prediction is built in-place in the Recon buffer. It is overwritten with Reconstructed
    //   samples if the CBF==1 && SKIP==False

    //**********************************
    // Luma
    //**********************************
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {

#if CHROMA_BLIND
        if (cu_ptr->prediction_mode_flag != INTRA_MODE || (cu_ptr->prediction_unit_array->intra_chroma_mode != UV_CFL_PRED && context_ptr->evaluate_cfl_ep == EB_FALSE))
#else
        if (cu_ptr->prediction_mode_flag != INTRA_MODE || cu_ptr->prediction_unit_array->intra_chroma_mode != UV_CFL_PRED)
#endif
        {
            predLumaOffset = (predSamples->origin_y + origin_y)             * predSamples->stride_y + (predSamples->origin_x + origin_x);
            if (txb_ptr->y_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {
                (void)asm_type;
                (void)transformScratchBuffer;
                uint8_t     *predBuffer = predSamples->buffer_y + predLumaOffset;
                av1_inv_transform_recon8bit(
                    ((int32_t*)residual16bit->buffer_y) + context_ptr->coded_area_sb,
                    predBuffer,
                    predSamples->stride_y,
                    context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                    txb_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    eob[0]
                );
            }
        }
    }

    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK)

        //**********************************
        // Chroma
        //**********************************

    {
        uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
        uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4
        predChromaOffset = (((predSamples->origin_y + round_origin_y) >> 1)           * predSamples->strideCb) + ((predSamples->origin_x + round_origin_x) >> 1);

        //**********************************
        // Cb
        //**********************************
        if (txb_ptr->u_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {


            uint8_t     *predBuffer = predSamples->bufferCb + predChromaOffset;

            av1_inv_transform_recon8bit(
                ((int32_t*)residual16bit->bufferCb) + context_ptr->coded_area_sb_uv,

                predBuffer,
                predSamples->strideCb,
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
                txb_ptr->transform_type[PLANE_TYPE_UV],
                PLANE_TYPE_UV,
                eob[1]);

        }

        //**********************************
        // Cr
        //**********************************
        predChromaOffset = (((predSamples->origin_y + round_origin_y) >> 1)           * predSamples->strideCr) + ((predSamples->origin_x + round_origin_x) >> 1);

        if (txb_ptr->v_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

            uint8_t     *predBuffer = predSamples->bufferCr + predChromaOffset;

            av1_inv_transform_recon8bit(
                ((int32_t*)residual16bit->bufferCr) + context_ptr->coded_area_sb_uv,
                predBuffer,
                predSamples->strideCr,
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
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
    EncDecContext_t         *context_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    EbPictureBufferDesc_t   *predSamples,     // no basis/offset
    EbPictureBufferDesc_t   *residual16bit,    // no basis/offset
    int16_t                 *transformScratchBuffer,
    uint32_t                 component_mask,
    uint16_t                *eob,
    EbAsm                    asm_type)
{

    uint32_t predLumaOffset;
    uint32_t predChromaOffset;
#if !QT_10BIT_SUPPORT
    uint32_t scratchLumaOffset;
    uint32_t scratchChromaOffset;
    uint32_t reconLumaOffset;
    uint32_t reconChromaOffset;
#endif

    CodingUnit_t          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];

#if QT_10BIT_SUPPORT
    (void)asm_type;
    (void)transformScratchBuffer;
#endif
    //**********************************
    // Luma
    //**********************************
    if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {
        if (cu_ptr->prediction_mode_flag != INTRA_MODE || cu_ptr->prediction_unit_array->intra_chroma_mode != UV_CFL_PRED)

        {
            predLumaOffset = (predSamples->origin_y + origin_y)* predSamples->stride_y + (predSamples->origin_x + origin_x);
#if !QT_10BIT_SUPPORT
            scratchLumaOffset = context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr] + context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr] * SB_STRIDE_Y;
            reconLumaOffset = (predSamples->origin_y + origin_y)* predSamples->stride_y + (predSamples->origin_x + origin_x);
#endif
            if (txb_ptr->y_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

#if QT_10BIT_SUPPORT
                uint16_t     *predBuffer = ((uint16_t*)predSamples->buffer_y) + predLumaOffset;
                av1_inv_transform_recon(
                    ((int32_t*)residual16bit->buffer_y) + context_ptr->coded_area_sb,
                    CONVERT_TO_BYTEPTR(predBuffer),
                    predSamples->stride_y,
                    context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                    BIT_INCREMENT_10BIT,
                    txb_ptr->transform_type[PLANE_TYPE_Y],
                    PLANE_TYPE_Y,
                    eob[0]
                );
#else
                av1_estimate_inv_transform(
                    ((int32_t*)residual16bit->buffer_y) + scratchLumaOffset,
                    64,
                    ((int32_t*)residual16bit->buffer_y) + scratchLumaOffset,
                    64,
                    txb_size,
                    transformScratchBuffer,
                    BIT_INCREMENT_10BIT,
                    txb_ptr->transform_type[PLANE_TYPE_Y],
                    eob[0],
                    asm_type,
                    0);

                picture_addition_kernel16_bit(
                    (uint16_t*)predSamples->buffer_y + predLumaOffset,
                    predSamples->stride_y,
                    ((int32_t*)residual16bit->buffer_y) + scratchLumaOffset,
                    64,
                    (uint16_t*)predSamples->buffer_y + reconLumaOffset,
                    predSamples->stride_y,
                    txb_size,
                    txb_size,
                    10);
#endif
            }

        }
    }

    if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK)

        //**********************************
        // Chroma
        //**********************************

    {

        //**********************************
        // Cb
        //**********************************

        uint32_t                 round_origin_x = (origin_x >> 3) << 3;// for Chroma blocks with size of 4
        uint32_t                 round_origin_y = (origin_y >> 3) << 3;// for Chroma blocks with size of 4

        predChromaOffset = (((predSamples->origin_y + round_origin_y) >> 1)           * predSamples->strideCb) + ((predSamples->origin_x + round_origin_x) >> 1);
#if !QT_10BIT_SUPPORT
        scratchChromaOffset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
        reconChromaOffset = (((predSamples->origin_y + origin_y) >> 1) * predSamples->strideCb) + ((predSamples->origin_x + origin_x) >> 1);
#endif

        if (txb_ptr->u_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {


#if QT_10BIT_SUPPORT
            uint16_t     *predBuffer = ((uint16_t*)predSamples->bufferCb) + predChromaOffset;
            av1_inv_transform_recon(
                ((int32_t*)residual16bit->bufferCb) + context_ptr->coded_area_sb_uv,
                CONVERT_TO_BYTEPTR(predBuffer),
                predSamples->strideCb,
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                PLANE_TYPE_UV,
                eob[1]);
#else
            av1_estimate_inv_transform(
                ((int32_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                ((int32_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                txb_size >> 1,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                eob[1],
                asm_type,
                0);

            picture_addition_kernel16_bit(
                (uint16_t*)predSamples->bufferCb + predChromaOffset,
                predSamples->strideCb,
                ((int32_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                (uint16_t*)predSamples->bufferCb + reconChromaOffset,
                predSamples->strideCb,
                txb_size >> 1,
                txb_size >> 1,
                10);
#endif

        }

        //**********************************
        // Cr
        //**********************************
        predChromaOffset = (((predSamples->origin_y + round_origin_y) >> 1)           * predSamples->strideCr) + ((predSamples->origin_x + round_origin_x) >> 1);
#if !QT_10BIT_SUPPORT
        scratchChromaOffset = ROUND_UV(context_ptr->blk_geom->tx_org_x[context_ptr->txb_itr]) / 2 + ROUND_UV(context_ptr->blk_geom->tx_org_y[context_ptr->txb_itr]) / 2 * SB_STRIDE_UV;
        reconChromaOffset = (((predSamples->origin_y + origin_y) >> 1) * predSamples->strideCr) + ((predSamples->origin_x + origin_x) >> 1);
#endif
        if (txb_ptr->v_has_coeff == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

#if QT_10BIT_SUPPORT
            uint16_t     *predBuffer = ((uint16_t*)predSamples->bufferCr) + predChromaOffset;
            av1_inv_transform_recon(
                ((int32_t*)residual16bit->bufferCr) + context_ptr->coded_area_sb_uv,
                CONVERT_TO_BYTEPTR(predBuffer),
                predSamples->strideCr,
                context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                PLANE_TYPE_UV,
                eob[2]);
#else
            av1_estimate_inv_transform(
                ((int32_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                ((int32_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                txb_size >> 1,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                txb_ptr->transform_type[PLANE_TYPE_UV],
                eob[2],
                asm_type,
                0);

            picture_addition_kernel16_bit(
                (uint16_t*)predSamples->bufferCr + predChromaOffset,
                predSamples->strideCr,
                ((int32_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                (uint16_t*)predSamples->bufferCr + reconChromaOffset,
                predSamples->strideCr,
                txb_size >> 1,
                txb_size >> 1,
                10);
#endif


        }
    }

    return;
}
#if !QT_10BIT_SUPPORT
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
static void EncodeGenerateRecon(
    EncDecContext_t       *context_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    EbPictureBufferDesc_t *predSamples,     // no basis/offset
    EbPictureBufferDesc_t *residual16bit,    // no basis/offset
    int16_t                *transformScratchBuffer,
    EbAsm                 asm_type)
{
    uint32_t predLumaOffset;
    uint32_t predChromaOffset;
    uint32_t scratchLumaOffset;
    uint32_t scratchChromaOffset;
    uint32_t reconLumaOffset;
    uint32_t reconChromaOffset;

    CodingUnit_t          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    uint32_t                 txb_size = context_ptr->cu_stats->size;

    EbPictureBufferDesc_t *reconSamples = predSamples;
    // *Note - The prediction is built in-place in the Recon buffer. It is overwritten with Reconstructed
    //   samples if the CBF==1 && SKIP==False

    //**********************************
    // Luma
    //**********************************

    {
        predLumaOffset = (predSamples->origin_y + origin_y)             * predSamples->stride_y + (predSamples->origin_x + origin_x);
        scratchLumaOffset = ((origin_y & (63)) * 64) + (origin_x & (63));
        reconLumaOffset = (reconSamples->origin_y + origin_y)            * reconSamples->stride_y + (reconSamples->origin_x + origin_x);
        if (txb_ptr->lumaCbf == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

            encode_inv_transform(
                txb_ptr->trans_coeff_shape_luma == ONLY_DC_SHAPE || txb_ptr->is_only_dc[0],
                ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
                64,
                ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
                64,
                txb_size,
                transformScratchBuffer,
                BIT_INCREMENT_8BIT,
                (EbBool)(txb_size == MIN_PU_SIZE),
                asm_type);

            addition_kernel_func_ptr_array[asm_type][txb_size >> 3](
                predSamples->buffer_y + predLumaOffset,
                predSamples->stride_y,
                ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
                64,
                reconSamples->buffer_y + reconLumaOffset,
                reconSamples->stride_y,
                txb_size,
                txb_size);
        }
    }

    //**********************************
    // Chroma
    //**********************************

    {
        predChromaOffset = (((predSamples->origin_y + origin_y) >> 1)           * predSamples->strideCb) + ((predSamples->origin_x + origin_x) >> 1);
        scratchChromaOffset = (((origin_y & (63)) >> 1) * 32) + ((origin_x & (63)) >> 1);
        reconChromaOffset = (((reconSamples->origin_y + origin_y) >> 1)          * reconSamples->strideCb) + ((reconSamples->origin_x + origin_x) >> 1);
        //**********************************
        // Cb
        //**********************************
        if (txb_ptr->cbCbf == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

            encode_inv_transform(
                txb_ptr->trans_coeff_shape_chroma == ONLY_DC_SHAPE || txb_ptr->is_only_dc[1],
                ((int16_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                ((int16_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                txb_size >> 1,
                transformScratchBuffer,
                BIT_INCREMENT_8BIT,
                EB_FALSE,
                asm_type);

            addition_kernel_func_ptr_array[asm_type][txb_size >> 4](
                predSamples->bufferCb + predChromaOffset,
                predSamples->strideCb,
                ((int16_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                reconSamples->bufferCb + reconChromaOffset,
                reconSamples->strideCb,
                txb_size >> 1,
                txb_size >> 1);
        }

        //**********************************
        // Cr
        //**********************************
        predChromaOffset = (((predSamples->origin_y + origin_y) >> 1)           * predSamples->strideCr) + ((predSamples->origin_x + origin_x) >> 1);
        scratchChromaOffset = (((origin_y & (63)) >> 1) * 32) + ((origin_x & (63)) >> 1);
        reconChromaOffset = (((reconSamples->origin_y + origin_y) >> 1)          * reconSamples->strideCr) + ((reconSamples->origin_x + origin_x) >> 1);
        if (txb_ptr->crCbf == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

            encode_inv_transform(
                txb_ptr->trans_coeff_shape_chroma == ONLY_DC_SHAPE || txb_ptr->is_only_dc[2],
                ((int16_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                ((int16_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                txb_size >> 1,
                transformScratchBuffer,
                BIT_INCREMENT_8BIT,
                EB_FALSE,
                asm_type);

            addition_kernel_func_ptr_array[asm_type][txb_size >> 4](
                predSamples->bufferCr + predChromaOffset,
                predSamples->strideCr,
                ((int16_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                reconSamples->bufferCr + reconChromaOffset,
                reconSamples->strideCr,
                txb_size >> 1,
                txb_size >> 1);
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
static void EncodeGenerateRecon16bit(
    EncDecContext_t       *context_ptr,
    uint32_t                 origin_x,
    uint32_t                 origin_y,
    EbPictureBufferDesc_t *predSamples,     // no basis/offset
    EbPictureBufferDesc_t *residual16bit,    // no basis/offset
    int16_t                *transformScratchBuffer,
    EbAsm                 asm_type)
{

    uint32_t predLumaOffset;
    uint32_t predChromaOffset;
    uint32_t scratchLumaOffset;
    uint32_t scratchChromaOffset;
    uint32_t reconLumaOffset;
    uint32_t reconChromaOffset;

    CodingUnit_t          *cu_ptr = context_ptr->cu_ptr;
    TransformUnit_t       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
    uint32_t                 txb_size = context_ptr->cu_stats->size;

    //**********************************
    // Luma
    //**********************************

    {
        predLumaOffset = (predSamples->origin_y + origin_y)* predSamples->stride_y + (predSamples->origin_x + origin_x);
        scratchLumaOffset = ((origin_y & (63)) * 64) + (origin_x & (63));
        reconLumaOffset = (predSamples->origin_y + origin_y)* predSamples->stride_y + (predSamples->origin_x + origin_x);
        if (txb_ptr->lumaCbf == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

            encode_inv_transform(
                txb_ptr->trans_coeff_shape_luma == ONLY_DC_SHAPE || txb_ptr->is_only_dc[0],
                ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
                64,
                ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
                64,
                txb_size,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                (EbBool)(txb_size == MIN_PU_SIZE),
                asm_type);

            addition_kernel_func_ptr_array16bit[asm_type](
                (uint16_t*)predSamples->buffer_y + predLumaOffset,
                predSamples->stride_y,
                ((int16_t*)residual16bit->buffer_y) + scratchLumaOffset,
                64,
                (uint16_t*)predSamples->buffer_y + reconLumaOffset,
                predSamples->stride_y,
                txb_size,
                txb_size);

        }

    }

    //**********************************
    // Chroma
    //**********************************

    {

        //**********************************
        // Cb
        //**********************************
        predChromaOffset = (((predSamples->origin_y + origin_y) >> 1)  * predSamples->strideCb) + ((predSamples->origin_x + origin_x) >> 1);
        scratchChromaOffset = (((origin_y & (63)) >> 1) * 32) + ((origin_x & (63)) >> 1);
        reconChromaOffset = (((predSamples->origin_y + origin_y) >> 1) * predSamples->strideCb) + ((predSamples->origin_x + origin_x) >> 1);
        if (txb_ptr->cbCbf == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

            encode_inv_transform(
                txb_ptr->trans_coeff_shape_chroma == ONLY_DC_SHAPE || txb_ptr->is_only_dc[1],
                ((int16_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                ((int16_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                txb_size >> 1,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                EB_FALSE,
                asm_type);

            addition_kernel_func_ptr_array16bit[asm_type](
                (uint16_t*)predSamples->bufferCb + predChromaOffset,
                predSamples->strideCb,
                ((int16_t*)residual16bit->bufferCb) + scratchChromaOffset,
                32,
                (uint16_t*)predSamples->bufferCb + reconChromaOffset,
                predSamples->strideCb,
                txb_size >> 1,
                txb_size >> 1);

        }

        //**********************************
        // Cr
        //**********************************
        predChromaOffset = (((predSamples->origin_y + origin_y) >> 1)  * predSamples->strideCr) + ((predSamples->origin_x + origin_x) >> 1);
        scratchChromaOffset = (((origin_y & (63)) >> 1) * 32) + ((origin_x & (63)) >> 1);
        reconChromaOffset = (((predSamples->origin_y + origin_y) >> 1) * predSamples->strideCr) + ((predSamples->origin_x + origin_x) >> 1);
        if (txb_ptr->crCbf == EB_TRUE && cu_ptr->skip_flag == EB_FALSE) {

            encode_inv_transform(
                txb_ptr->trans_coeff_shape_chroma == ONLY_DC_SHAPE || txb_ptr->is_only_dc[2],
                ((int16_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                ((int16_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                txb_size >> 1,
                transformScratchBuffer,
                BIT_INCREMENT_10BIT,
                EB_FALSE,
                asm_type);

            addition_kernel_func_ptr_array16bit[asm_type](
                (uint16_t*)predSamples->bufferCr + predChromaOffset,
                predSamples->strideCr,
                ((int16_t*)residual16bit->bufferCr) + scratchChromaOffset,
                32,
                (uint16_t*)predSamples->bufferCr + reconChromaOffset,
                predSamples->strideCr,
                txb_size >> 1,
                txb_size >> 1);

        }
    }

    return;
}

#endif
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
#if !QT_10BIT_SUPPORT

EB_GENERATE_RECON_FUNC_PTR   EncodeGenerateReconFunctionPtr[2] =
{
    EncodeGenerateRecon,
    EncodeGenerateRecon16bit
};
#endif

#if !QT_10BIT_SUPPORT
EB_GENERATE_RECON_INTRA_4x4_FUNC_PTR   EncodeGenerateReconIntra4x4FunctionPtr[2] =
{
    EncodeGenerateReconIntra4x4,
    EncodeGenerateReconIntra4x416bit
};

EB_GENERATE_INTRA_SAMPLES_FUNC_PTR GenerateIntraReferenceSamplesFuncTable[2] =
{
    GenerateIntraReferenceSamplesEncodePass,
    GenerateIntraReference16bitSamplesEncodePass
};

EB_ENC_PASS_INTRA_FUNC_PTR EncodePassIntraPredictionFuncTable[2] =
{
    EncodePassIntraPrediction,
    EncodePassIntraPrediction16bit
};
#endif

/*******************************************
* Encode Pass - Assign Delta Qp
*******************************************/
static void EncodePassUpdateQp(
    PictureControlSet_t     *picture_control_set_ptr,
    EncDecContext_t         *context_ptr,
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

    if (Log2f(context_ptr->blk_geom->bwidth) >= log2MinCuQpDeltaSize) {

        *isDeltaQpNotCoded = EB_TRUE;
    }

    // At the beginning of a new quantization group
    if (((context_ptr->cu_origin_x & ((1 << log2MinCuQpDeltaSize) - 1)) == 0) &&
        ((context_ptr->cu_origin_y & ((1 << log2MinCuQpDeltaSize) - 1)) == 0))
    {
        *isDeltaQpNotCoded = EB_TRUE;
        newQuantGroup = EB_TRUE;
    }
    else {
        newQuantGroup = EB_FALSE;
    }

    // setting the previous Quantization Group QP
    if (newQuantGroup == EB_TRUE) {
        *prev_coded_qp = *prev_quant_group_coded_qp;
    }

    if (sameLcuCheckTop) {
        qpTopNeighborIndex =
            LUMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(
                quantGroupX,
                quantGroupY - 1,
                picture_control_set_ptr->qp_array_stride);
        qpTopNeighbor = picture_control_set_ptr->qp_array[qpTopNeighborIndex];
    }
    else {
        qpTopNeighbor = *prev_coded_qp;
    }

    if (sameLcuCheckLeft) {
        qpLeftNeighborIndex =
            LUMA_SAMPLE_PIC_WISE_LOCATION_TO_QP_ARRAY_IDX(
                quantGroupX - 1,
                quantGroupY,
                picture_control_set_ptr->qp_array_stride);

        qpLeftNeighbor = picture_control_set_ptr->qp_array[qpLeftNeighborIndex];
    }
    else {
        qpLeftNeighbor = *prev_coded_qp;
    }

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
    else {
        qp = (uint8_t)sb_qp;
    }
    context_ptr->cu_ptr->qp = qp;
    return;
}



EbErrorType QpmDeriveBeaAndSkipQpmFlagLcu(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    LargestCodingUnit_t                    *sb_ptr,
    uint32_t                                 sb_index,
    EncDecContext_t                        *context_ptr)
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


    context_ptr->qpmQp = picture_qp;

    SbStat_t *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);


    context_ptr->non_moving_delta_qp = 0;

    context_ptr->grass_enhancement_flag = ((picture_control_set_ptr->scene_caracteristic_id == EB_FRAME_CARAC_1) && (sb_stat_ptr->cu_stat_array[0].grass_area)
        && (sb_ptr->picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index].edge_block_num > 0))

        ? EB_TRUE : EB_FALSE;


    context_ptr->backgorund_enhancement = EB_FALSE;


    context_ptr->skip_qpm_flag = sequence_control_set_ptr->static_config.improve_sharpness ? EB_FALSE : EB_TRUE;

    if ((picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag == EB_FALSE) && ((picture_control_set_ptr->parent_pcs_ptr->pic_noise_class >= PIC_NOISE_CLASS_3_1) || (picture_control_set_ptr->parent_pcs_ptr->high_dark_low_light_area_density_flag) || (picture_control_set_ptr->parent_pcs_ptr->intra_coded_block_probability > 90))) {
        context_ptr->skip_qpm_flag = EB_TRUE;
    }

    if (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE) {
        context_ptr->skip_qpm_flag = EB_TRUE;
    }

#if ADD_DELTA_QP_SUPPORT
    context_ptr->skip_qpm_flag = EB_FALSE;
#endif

    if (context_ptr->skip_qpm_flag == EB_FALSE) {
        if (picture_control_set_ptr->parent_pcs_ptr->pic_homogenous_over_time_sb_percentage > 30 && picture_control_set_ptr->slice_type != I_SLICE) {
#if ADD_DELTA_QP_SUPPORT
            context_ptr->qpmQp = CLIP3(min_qp_allowed, max_qp_allowed, picture_qp + deltaQpRes);
#else
            context_ptr->qpmQp = CLIP3(min_qp_allowed, max_qp_allowed, picture_qp + 1);
#endif
        }
    }


    return return_error;
}
#if ADD_DELTA_QP_SUPPORT
/*****************************************************************************
* NM - Note: Clean up
* AV1 QPM is SB based and all sub-Lcu buffers needs to be removed
******************************************************************************/
EbErrorType Av1QpModulationLcu(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    LargestCodingUnit_t                    *sb_ptr,
    uint32_t                                  sb_index,
    uint8_t                                   type,
    EncDecContext_t                        *context_ptr)
{
    EbErrorType                            return_error = EB_ErrorNone;

    int64_t                                  complexityDistance;
    int8_t                                   delta_qp = 0;
    uint16_t                                  qpmQp = context_ptr->qpmQp;
    uint16_t                                  min_qp_allowed = 0;
    uint16_t                                  max_qp_allowed = 255;
    uint16_t                                  cu_qp;
    EbBool                                 acEnergyBasedAntiContouring = picture_control_set_ptr->slice_type == I_SLICE ? EB_TRUE : EB_FALSE;
    uint8_t                                   lowerQPClass;

    int8_t    non_moving_delta_qp = context_ptr->non_moving_delta_qp;
    int8_t    bea64x64DeltaQp;

    uint8_t   deltaQpRes = picture_control_set_ptr->parent_pcs_ptr->delta_q_res;

    cu_qp = qpmQp;
    sb_ptr->qp = qpmQp;

    uint32_t  distortion = 0;

    if (!context_ptr->skip_qpm_flag) {

        // INTRA MODE
        if (type == INTRA_MODE) {




            OisCu32Cu16Results_t  *oisCu32Cu16ResultsPtr = picture_control_set_ptr->parent_pcs_ptr->ois_cu32_cu16_results[sb_index];
            //OisCu8Results_t         *oisCu8ResultsPtr = picture_control_set_ptr->parent_pcs_ptr->ois_cu8_results[sb_index];

            distortion =
                oisCu32Cu16ResultsPtr->sorted_ois_candidate[1][0].distortion +
                oisCu32Cu16ResultsPtr->sorted_ois_candidate[2][0].distortion +
                oisCu32Cu16ResultsPtr->sorted_ois_candidate[3][0].distortion +
                oisCu32Cu16ResultsPtr->sorted_ois_candidate[4][0].distortion;



            distortion = (uint32_t)CLIP3(picture_control_set_ptr->parent_pcs_ptr->intra_complexity_min[0], picture_control_set_ptr->parent_pcs_ptr->intra_complexity_max[0], distortion);
            complexityDistance = ((int32_t)distortion - (int32_t)picture_control_set_ptr->parent_pcs_ptr->intra_complexity_avg[0]);

            if (complexityDistance < 0) {

                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->intra_min_distance[0] != 0) ? (int8_t)((context_ptr->min_delta_qp_weight * context_ptr->min_delta_qp[0] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->intra_min_distance[0])) : 0;
            }
            else {

                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->intra_max_distance[0] != 0) ? (int8_t)((context_ptr->max_delta_qp_weight * context_ptr->max_delta_qp[0] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->intra_max_distance[0])) : 0;
            }



        }
        // INTER MODE
        else {


            distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][0].distortionDirection[0].distortion;


            distortion = (uint32_t)CLIP3(picture_control_set_ptr->parent_pcs_ptr->inter_complexity_min[0], picture_control_set_ptr->parent_pcs_ptr->inter_complexity_max[0], distortion);
            complexityDistance = ((int32_t)distortion - (int32_t)picture_control_set_ptr->parent_pcs_ptr->inter_complexity_avg[0]);

            if (complexityDistance < 0) {

                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->inter_min_distance[0] != 0) ? (int8_t)((context_ptr->min_delta_qp_weight * context_ptr->min_delta_qp[0] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->inter_min_distance[0])) : 0;
            }
            else {

                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->inter_max_distance[0] != 0) ? (int8_t)((context_ptr->max_delta_qp_weight * context_ptr->max_delta_qp[0] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->inter_max_distance[0])) : 0;
            }

        }

        if (context_ptr->backgorund_enhancement) {
            // Use the 8x8 background enhancement only for the Intra slice, otherwise, use the existing SB based BEA results
            bea64x64DeltaQp = non_moving_delta_qp;

            if ((picture_control_set_ptr->parent_pcs_ptr->yMean[sb_index][0] > ANTI_CONTOURING_LUMA_T2) || (picture_control_set_ptr->parent_pcs_ptr->yMean[sb_index][0] < ANTI_CONTOURING_LUMA_T1)) {

                if (bea64x64DeltaQp < 0) {
                    bea64x64DeltaQp = 0;
                }

            }

            delta_qp += bea64x64DeltaQp;
        }

        if ((picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag)) {
            delta_qp = (delta_qp < context_ptr->min_delta_qp[0]) ? delta_qp : context_ptr->min_delta_qp[0];
        }

        SbStat_t *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);
        if (sb_stat_ptr->stationary_edge_over_time_flag && delta_qp > 0) {
            delta_qp = 0;
        }

        if (acEnergyBasedAntiContouring) {

            lowerQPClass = DeriveContouringClass(
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

            if (qpmQp > (RC_QPMOD_MAXQP * deltaQpRes)) {
                delta_qp = MIN(0, delta_qp);
            }


            cu_qp = (uint32_t)(qpmQp + delta_qp);


            if ((qpmQp <= (RC_QPMOD_MAXQP *deltaQpRes))) {
                cu_qp = (uint8_t)CLIP3(
                    min_qp_allowed,
                    RC_QPMOD_MAXQP*deltaQpRes,
                    cu_qp);
            }
        }
        else {
            cu_qp = (uint8_t)(qpmQp + delta_qp);
        }

        cu_qp = (uint8_t)CLIP3(
            min_qp_allowed,
            max_qp_allowed,
            cu_qp);


    }

    sb_ptr->qp = sequence_control_set_ptr->static_config.improve_sharpness ? cu_qp : qpmQp;


    sb_ptr->delta_qp = (int16_t)sb_ptr->qp - (int16_t)qpmQp;

    sb_ptr->org_delta_qp = sb_ptr->delta_qp;

    if (sb_ptr->delta_qp % deltaQpRes != 0)
        printf("Qpm_error: delta_qp must be multiplier of deltaQpRes\n");
    if (sb_ptr->qp == 0)
        printf("Qpm_error: qp must be greater than 0 when use_delta_q is ON\n");

    return return_error;
}

#endif
EbErrorType EncQpmDeriveDeltaQPForEachLeafLcu(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    LargestCodingUnit_t                    *sb_ptr,
    uint32_t                                  sb_index,
    CodingUnit_t                           *cu_ptr,
    uint32_t                                  cu_depth,
    uint32_t                                  cu_index,
    uint32_t                                  cu_size,
    uint8_t                                   type,
    uint8_t                                   parent32x32Index,
    EncDecContext_t                        *context_ptr)
{
    EbErrorType                    return_error = EB_ErrorNone;


    //SbParams_t                        sb_params;
    int64_t                          complexityDistance;
    int8_t                           delta_qp = 0;
    uint8_t                           qpmQp = (uint8_t)context_ptr->qpmQp;
    uint8_t                           min_qp_allowed = (uint8_t)sequence_control_set_ptr->static_config.min_qp_allowed;
    uint8_t                           max_qp_allowed = (uint8_t)sequence_control_set_ptr->static_config.max_qp_allowed;
    uint8_t                           cu_qp;

    EbBool  use16x16Stat = EB_FALSE;

    uint32_t usedDepth = cu_depth;
    if (use16x16Stat)
        usedDepth = 2;

    uint32_t cuIndexInRaterScan = MD_SCAN_TO_RASTER_SCAN[cu_index];

    EbBool                         acEnergyBasedAntiContouring = picture_control_set_ptr->slice_type == I_SLICE ? EB_TRUE : EB_FALSE;
    uint8_t                           lowerQPClass;

    int8_t    non_moving_delta_qp = context_ptr->non_moving_delta_qp;
    int8_t    bea64x64DeltaQp;

    cu_qp = qpmQp;
    cu_ptr->qp = qpmQp;

    uint32_t  distortion = 0;

    if (!context_ptr->skip_qpm_flag) {

        // INTRA MODE
        if (type == INTRA_MODE) {




            OisCu32Cu16Results_t  *oisCu32Cu16ResultsPtr = picture_control_set_ptr->parent_pcs_ptr->ois_cu32_cu16_results[sb_index];
            OisCu8Results_t         *oisCu8ResultsPtr = picture_control_set_ptr->parent_pcs_ptr->ois_cu8_results[sb_index];

            if (cu_size > 32) {
                distortion =
                    oisCu32Cu16ResultsPtr->sorted_ois_candidate[1][0].distortion +
                    oisCu32Cu16ResultsPtr->sorted_ois_candidate[2][0].distortion +
                    oisCu32Cu16ResultsPtr->sorted_ois_candidate[3][0].distortion +
                    oisCu32Cu16ResultsPtr->sorted_ois_candidate[4][0].distortion;
            }
            else if (cu_size == 32) {
                const uint32_t me2Nx2NTableOffset = context_ptr->cu_stats->cuNumInDepth + me2Nx2NOffset[context_ptr->cu_stats->depth];
                distortion = oisCu32Cu16ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].distortion;
            }
            else {
                if (cu_size > 8) {
                    const uint32_t me2Nx2NTableOffset = context_ptr->cu_stats->cuNumInDepth + me2Nx2NOffset[context_ptr->cu_stats->depth];
                    distortion = oisCu32Cu16ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].distortion;
                }
                else {


                    if (use16x16Stat) {

                        const CodedUnitStats_t  *cu_stats = GetCodedUnitStats(ParentBlockIndex[cu_index]);
                        const uint32_t me2Nx2NTableOffset = cu_stats->cuNumInDepth + me2Nx2NOffset[cu_stats->depth];

                        distortion = oisCu32Cu16ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].distortion;
                    }
                    else {



                        const uint32_t me2Nx2NTableOffset = context_ptr->cu_stats->cuNumInDepth;

                        if (oisCu8ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].valid_distortion) {
                            distortion = oisCu8ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].distortion;
                        }
                        else {

                            const CodedUnitStats_t  *cu_stats = GetCodedUnitStats(ParentBlockIndex[cu_index]);
                            const uint32_t me2Nx2NTableOffset = cu_stats->cuNumInDepth + me2Nx2NOffset[cu_stats->depth];

                            if (oisCu32Cu16ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].valid_distortion) {
                                distortion = oisCu32Cu16ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].distortion;
                            }
                            else {
                                distortion = 0;
                            }
                        }

                    }


                }
            }






            distortion = (uint32_t)CLIP3(picture_control_set_ptr->parent_pcs_ptr->intra_complexity_min[usedDepth], picture_control_set_ptr->parent_pcs_ptr->intra_complexity_max[usedDepth], distortion);
            complexityDistance = ((int32_t)distortion - (int32_t)picture_control_set_ptr->parent_pcs_ptr->intra_complexity_avg[usedDepth]);

            if (complexityDistance < 0) {

                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->intra_min_distance[usedDepth] != 0) ? (int8_t)((context_ptr->min_delta_qp_weight * context_ptr->min_delta_qp[usedDepth] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->intra_min_distance[usedDepth])) : 0;
            }
            else {

                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->intra_max_distance[usedDepth] != 0) ? (int8_t)((context_ptr->max_delta_qp_weight * context_ptr->max_delta_qp[usedDepth] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->intra_max_distance[usedDepth])) : 0;
            }



        }
        // INTER MODE
        else {


            distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][cuIndexInRaterScan].distortionDirection[0].distortion;



            if (use16x16Stat) {
                uint32_t cuIndexRScan = MD_SCAN_TO_RASTER_SCAN[ParentBlockIndex[cu_index]];

                distortion = picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][cuIndexRScan].distortionDirection[0].distortion;

            }
            distortion = (uint32_t)CLIP3(picture_control_set_ptr->parent_pcs_ptr->inter_complexity_min[usedDepth], picture_control_set_ptr->parent_pcs_ptr->inter_complexity_max[usedDepth], distortion);
            complexityDistance = ((int32_t)distortion - (int32_t)picture_control_set_ptr->parent_pcs_ptr->inter_complexity_avg[usedDepth]);

            if (complexityDistance < 0) {

                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->inter_min_distance[usedDepth] != 0) ? (int8_t)((context_ptr->min_delta_qp_weight * context_ptr->min_delta_qp[usedDepth] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->inter_min_distance[usedDepth])) : 0;
            }
            else {

                delta_qp = (picture_control_set_ptr->parent_pcs_ptr->inter_max_distance[usedDepth] != 0) ? (int8_t)((context_ptr->max_delta_qp_weight * context_ptr->max_delta_qp[usedDepth] * complexityDistance) / (100 * picture_control_set_ptr->parent_pcs_ptr->inter_max_distance[usedDepth])) : 0;
            }



        }

        if (context_ptr->backgorund_enhancement) {
            // Use the 8x8 background enhancement only for the Intra slice, otherwise, use the existing SB based BEA results
            bea64x64DeltaQp = non_moving_delta_qp;

            if (((cu_index > 0) && ((picture_control_set_ptr->parent_pcs_ptr->yMean[sb_index][parent32x32Index]) > ANTI_CONTOURING_LUMA_T2 || (picture_control_set_ptr->parent_pcs_ptr->yMean[sb_index][parent32x32Index]) < ANTI_CONTOURING_LUMA_T1)) ||
                ((cu_index == 0) && ((picture_control_set_ptr->parent_pcs_ptr->yMean[sb_index][0]) > ANTI_CONTOURING_LUMA_T2 || (picture_control_set_ptr->parent_pcs_ptr->yMean[sb_index][0]) < ANTI_CONTOURING_LUMA_T1))) {

                if (bea64x64DeltaQp < 0) {
                    bea64x64DeltaQp = 0;
                }

            }

            delta_qp += bea64x64DeltaQp;
        }

        if ((picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag)) {
            delta_qp = (delta_qp < context_ptr->min_delta_qp[0]) ? delta_qp : context_ptr->min_delta_qp[0];
        }

        SbStat_t *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);
        if (sb_stat_ptr->stationary_edge_over_time_flag && delta_qp > 0) {
            delta_qp = 0;
        }

        if (acEnergyBasedAntiContouring) {

            lowerQPClass = DeriveContouringClass(
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


        delta_qp -= context_ptr->grass_enhancement_flag ? 3 : 0;

        if (sequence_control_set_ptr->static_config.rate_control_mode == 1 || sequence_control_set_ptr->static_config.rate_control_mode == 2) {

            if (qpmQp > RC_QPMOD_MAXQP) {
                delta_qp = MIN(0, delta_qp);
            }

            cu_qp = (uint32_t)(qpmQp + delta_qp);


            if ((qpmQp <= RC_QPMOD_MAXQP)) {
                cu_qp = (uint8_t)CLIP3(
                    min_qp_allowed,
                    RC_QPMOD_MAXQP,
                    cu_qp);
            }
        }
        else {
            cu_qp = (uint8_t)(qpmQp + delta_qp);
        }

        cu_qp = (uint8_t)CLIP3(
            min_qp_allowed,
            max_qp_allowed,
            cu_qp);


    }

    cu_ptr->qp = sequence_control_set_ptr->static_config.improve_sharpness ? cu_qp : qpmQp;

    sb_ptr->qp = (cu_size == 64) ? (uint8_t)cu_ptr->qp : sb_ptr->qp;


    cu_ptr->delta_qp = (int16_t)cu_ptr->qp - (int16_t)qpmQp;

    cu_ptr->org_delta_qp = cu_ptr->delta_qp;


    return return_error;
}

void Store16bitInputSrc(
    EncDecContext_t         *context_ptr,
    PictureControlSet_t     *picture_control_set_ptr,
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
    {
        memcpy(toPtr + rowIt * picture_control_set_ptr->input_frame16bit->stride_y, fromPtr + rowIt * context_ptr->input_sample16bit_buffer->stride_y, lcuW * 2);
    }

    lcuX = lcuX / 2;
    lcuY = lcuY / 2;
    lcuW = lcuW / 2;
    lcuH = lcuH / 2;

    fromPtr = (uint16_t*)context_ptr->input_sample16bit_buffer->bufferCb;
    toPtr = (uint16_t*)picture_control_set_ptr->input_frame16bit->bufferCb + (lcuX + picture_control_set_ptr->input_frame16bit->origin_x / 2) + (lcuY + picture_control_set_ptr->input_frame16bit->origin_y / 2)*picture_control_set_ptr->input_frame16bit->strideCb;

    for (rowIt = 0; rowIt < lcuH; rowIt++)
    {
        memcpy(toPtr + rowIt * picture_control_set_ptr->input_frame16bit->strideCb, fromPtr + rowIt * context_ptr->input_sample16bit_buffer->strideCb, lcuW * 2);
    }

    fromPtr = (uint16_t*)context_ptr->input_sample16bit_buffer->bufferCr;
    toPtr = (uint16_t*)picture_control_set_ptr->input_frame16bit->bufferCr + (lcuX + picture_control_set_ptr->input_frame16bit->origin_x / 2) + (lcuY + picture_control_set_ptr->input_frame16bit->origin_y / 2)*picture_control_set_ptr->input_frame16bit->strideCb;

    for (rowIt = 0; rowIt < lcuH; rowIt++)
    {
        memcpy(toPtr + rowIt * picture_control_set_ptr->input_frame16bit->strideCr, fromPtr + rowIt * context_ptr->input_sample16bit_buffer->strideCr, lcuW * 2);
    }


}

void update_av1_mi_map(
    CodingUnit_t                   *cu_ptr,
    uint32_t                          cu_origin_x,
    uint32_t                          cu_origin_y,
    const BlockGeom                 *blk_geom,
    PictureControlSet_t            *picture_control_set_ptr);

void move_cu_data(
    CodingUnit_t *src_cu,
    CodingUnit_t *dst_cu);

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
EB_EXTERN void AV1EncodePass(
    SequenceControlSet_t      *sequence_control_set_ptr,
    PictureControlSet_t       *picture_control_set_ptr,
    LargestCodingUnit_t       *sb_ptr,
    uint32_t                   tbAddr,
    uint32_t                   sb_origin_x,
    uint32_t                   sb_origin_y,
    uint32_t                   sb_qp,
    EncDecContext_t           *context_ptr)
{

    EbBool                    is16bit = context_ptr->is16bit;
    EbPictureBufferDesc_t    *recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;
    EbPictureBufferDesc_t    *coeff_buffer_sb = sb_ptr->quantized_coeff;
    EbPictureBufferDesc_t    *inputPicture;
    ModeDecisionContext_t    *mdcontextPtr;
    mdcontextPtr = context_ptr->md_context;
    inputPicture = context_ptr->input_samples = (EbPictureBufferDesc_t*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;

    SbStat_t                *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[tbAddr]);

    EbBool                    availableCoeff;

    // QP Neighbor Arrays
    EbBool                    isDeltaQpNotCoded = EB_TRUE;

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
    EncodeContext_t          *encode_context_ptr;
    uint32_t                  lcuRowIndex = sb_origin_y / BLOCK_SIZE_64;

    // Dereferencing early
    NeighborArrayUnit_t      *ep_mode_type_neighbor_array = picture_control_set_ptr->ep_mode_type_neighbor_array;
    NeighborArrayUnit_t      *ep_intra_luma_mode_neighbor_array = picture_control_set_ptr->ep_intra_luma_mode_neighbor_array;
    NeighborArrayUnit_t      *ep_intra_chroma_mode_neighbor_array = picture_control_set_ptr->ep_intra_chroma_mode_neighbor_array;
    NeighborArrayUnit_t      *ep_mv_neighbor_array = picture_control_set_ptr->ep_mv_neighbor_array;
    NeighborArrayUnit_t      *ep_luma_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_luma_recon_neighbor_array16bit : picture_control_set_ptr->ep_luma_recon_neighbor_array;
    NeighborArrayUnit_t      *ep_cb_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_cb_recon_neighbor_array16bit : picture_control_set_ptr->ep_cb_recon_neighbor_array;
    NeighborArrayUnit_t      *ep_cr_recon_neighbor_array = is16bit ? picture_control_set_ptr->ep_cr_recon_neighbor_array16bit : picture_control_set_ptr->ep_cr_recon_neighbor_array;
    NeighborArrayUnit_t      *ep_skip_flag_neighbor_array = picture_control_set_ptr->ep_skip_flag_neighbor_array;

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

    EntropyCoder_t  *coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;

    uint32_t           dZoffset = 0;

    if (!sb_stat_ptr->stationary_edge_over_time_flag && sequence_control_set_ptr->static_config.improve_sharpness && picture_control_set_ptr->parent_pcs_ptr->pic_noise_class < PIC_NOISE_CLASS_3_1) {

        int16_t cuDeltaQp = (int16_t)(sb_ptr->qp - picture_control_set_ptr->parent_pcs_ptr->average_qp);
        uint32_t dzCondition = cuDeltaQp > 0 ? 0 : 1;

        if (sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE) {

            if (!(picture_control_set_ptr->parent_pcs_ptr->is_pan ||
                (picture_control_set_ptr->parent_pcs_ptr->non_moving_index_average < 10 && sb_ptr->aura_status_iii) ||
                (sb_stat_ptr->cu_stat_array[0].skin_area) ||
                (picture_control_set_ptr->parent_pcs_ptr->intra_coded_block_probability > 90) ||
                (picture_control_set_ptr->parent_pcs_ptr->high_dark_area_density_flag))) {

                if (picture_control_set_ptr->slice_type != I_SLICE &&
                    picture_control_set_ptr->temporal_layer_index == 0 &&
                    picture_control_set_ptr->parent_pcs_ptr->intra_coded_block_probability > 60 &&
                    !picture_control_set_ptr->parent_pcs_ptr->is_tilt &&
                    picture_control_set_ptr->parent_pcs_ptr->pic_homogenous_over_time_sb_percentage > 40)
                {
                    dZoffset = 10;
                }

                if (dzCondition) {
                    if (picture_control_set_ptr->scene_caracteristic_id == EB_FRAME_CARAC_1) {
                        if (picture_control_set_ptr->slice_type == I_SLICE) {
                            dZoffset = sb_stat_ptr->cu_stat_array[0].grass_area ? 10 : dZoffset;
                        }
                        else if (picture_control_set_ptr->temporal_layer_index == 0) {
                            dZoffset = sb_stat_ptr->cu_stat_array[0].grass_area ? 9 : dZoffset;
                        }
                        else if (picture_control_set_ptr->temporal_layer_index == 1) {
                            dZoffset = sb_stat_ptr->cu_stat_array[0].grass_area ? 5 : dZoffset;
                        }
                    }
                }
            }
        }
    }
    if (sequence_control_set_ptr->static_config.improve_sharpness) {

        QpmDeriveBeaAndSkipQpmFlagLcu(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_ptr,
            tbAddr,
            context_ptr);
    }
    else {
        context_ptr->skip_qpm_flag = EB_TRUE;
    }


    encode_context_ptr = ((SequenceControlSet_t*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr))->encode_context_ptr;

    if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE) {

        //get the 16bit form of the input LCU
        if (is16bit) {

            recon_buffer = ((EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->referencePicture16bit;

        }

        else {
            recon_buffer = ((EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->referencePicture;
        }
    }
    else { // non ref pictures
        recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;
    }


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

            const uint32_t inputLumaOffset = ((sb_origin_y + inputPicture->origin_y)         * inputPicture->stride_y) + (sb_origin_x + inputPicture->origin_x);
            const uint32_t inputCbOffset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->strideCb) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint32_t inputCrOffset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->strideCr) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint16_t luma2BitWidth = inputPicture->width / 4;
            const uint16_t chroma2BitWidth = inputPicture->width / 8;


            compressed_pack_lcu(
                inputPicture->buffer_y + inputLumaOffset,
                inputPicture->stride_y,
                inputPicture->bufferBitIncY + sb_origin_y * luma2BitWidth + (sb_origin_x / 4)*sb_height,
                sb_width / 4,
                (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_y,
                context_ptr->input_sample16bit_buffer->stride_y,
                sb_width,
                sb_height,
                asm_type);

            compressed_pack_lcu(
                inputPicture->bufferCb + inputCbOffset,
                inputPicture->strideCb,
                inputPicture->bufferBitIncCb + sb_origin_y / 2 * chroma2BitWidth + (sb_origin_x / 8)*(sb_height / 2),
                sb_width / 8,
                (uint16_t *)context_ptr->input_sample16bit_buffer->bufferCb,
                context_ptr->input_sample16bit_buffer->strideCb,
                sb_width >> 1,
                sb_height >> 1,
                asm_type);

            compressed_pack_lcu(
                inputPicture->bufferCr + inputCrOffset,
                inputPicture->strideCr,
                inputPicture->bufferBitIncCr + sb_origin_y / 2 * chroma2BitWidth + (sb_origin_x / 8)*(sb_height / 2),
                sb_width / 8,
                (uint16_t *)context_ptr->input_sample16bit_buffer->bufferCr,
                context_ptr->input_sample16bit_buffer->strideCr,
                sb_width >> 1,
                sb_height >> 1,
                asm_type);

        }
        else {

            const uint32_t inputLumaOffset = ((sb_origin_y + inputPicture->origin_y)         * inputPicture->stride_y) + (sb_origin_x + inputPicture->origin_x);
            const uint32_t inputBitIncLumaOffset = ((sb_origin_y + inputPicture->origin_y)         * inputPicture->strideBitIncY) + (sb_origin_x + inputPicture->origin_x);
            const uint32_t inputCbOffset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->strideCb) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint32_t inputBitIncCbOffset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->strideBitIncCb) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint32_t inputCrOffset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->strideCr) + ((sb_origin_x + inputPicture->origin_x) >> 1);
            const uint32_t inputBitIncCrOffset = (((sb_origin_y + inputPicture->origin_y) >> 1)  * inputPicture->strideBitIncCr) + ((sb_origin_x + inputPicture->origin_x) >> 1);

            pack2d_src(
                inputPicture->buffer_y + inputLumaOffset,
                inputPicture->stride_y,
                inputPicture->bufferBitIncY + inputBitIncLumaOffset,
                inputPicture->strideBitIncY,
                (uint16_t *)context_ptr->input_sample16bit_buffer->buffer_y,
                context_ptr->input_sample16bit_buffer->stride_y,
                sb_width,
                sb_height,
                asm_type);


            pack2d_src(
                inputPicture->bufferCb + inputCbOffset,
                inputPicture->strideCr,
                inputPicture->bufferBitIncCb + inputBitIncCbOffset,
                inputPicture->strideBitIncCr,
                (uint16_t *)context_ptr->input_sample16bit_buffer->bufferCb,
                context_ptr->input_sample16bit_buffer->strideCb,
                sb_width >> 1,
                sb_height >> 1,
                asm_type);


            pack2d_src(
                inputPicture->bufferCr + inputCrOffset,
                inputPicture->strideCr,
                inputPicture->bufferBitIncCr + inputBitIncCrOffset,
                inputPicture->strideBitIncCr,
                (uint16_t *)context_ptr->input_sample16bit_buffer->bufferCr,
                context_ptr->input_sample16bit_buffer->strideCr,
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

                EbReferenceObject_t  *refObjL0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;
                EbReferenceObject_t  *refObjL1 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr;
                uint32_t const TH = (sequence_control_set_ptr->static_config.frame_rate >> 16) < 50 ? 25 : 30;

                if ((refObjL0->tmpLayerIdx == 2 && refObjL0->intra_coded_area > TH) || (refObjL1->tmpLayerIdx == 2 && refObjL1->intra_coded_area > TH))
                    highIntraRef = EB_TRUE;

            }
            if (highIntraRef == EB_FALSE) {

                checkZeroLumaCbf = EB_TRUE;
            }
        }
    }
    context_ptr->intra_coded_area_sb[tbAddr] = 0;

    context_ptr->trans_coeff_shape_luma = 0;
    context_ptr->trans_coeff_shape_chroma = 0;
    context_ptr->coded_area_sb = 0;
    context_ptr->coded_area_sb_uv = 0;

#if AV1_LF 
    if (dlfEnableFlag && picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode == 1){        
        if (tbAddr == 0) {
            av1_loop_filter_init(picture_control_set_ptr);

            av1_pick_filter_level(
#if FILT_PROC
                0,
#else
                context_ptr,
#endif
                (EbPictureBufferDesc_t*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
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





    uint32_t final_cu_itr = 0;

    // CU Loop

    uint32_t    blk_it = 0;

    while (blk_it < sequence_control_set_ptr->max_block_cnt) {

        CodingUnit_t  *cu_ptr = context_ptr->cu_ptr = &context_ptr->md_context->md_cu_arr_nsq[blk_it];
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
                PredictionUnit_t        *pu_ptr = (PredictionUnit_t *)EB_NULL; //  done
                EbPictureBufferDesc_t   *residual_buffer = context_ptr->residual_buffer;
                EbPictureBufferDesc_t   *transform_buffer = context_ptr->transform_buffer;

                EbPictureBufferDesc_t   *inverse_quant_buffer = context_ptr->inverse_quant_buffer;

                int16_t                  *transform_inner_array_ptr = context_ptr->transform_inner_array_ptr;

                CodingUnit_t            *cu_ptr = context_ptr->cu_ptr = &context_ptr->md_context->md_cu_arr_nsq[d1_itr];

                context_ptr->cu_origin_x = (uint16_t)(sb_origin_x + blk_geom->origin_x);
                context_ptr->cu_origin_y = (uint16_t)(sb_origin_y + blk_geom->origin_y);
                cu_ptr->delta_qp = 0;
                cu_ptr->block_has_coeff = 0;


                // if(picture_control_set_ptr->picture_number==4 && context_ptr->cu_origin_x==0 && context_ptr->cu_origin_y==0)
                //     printf("CHEDD");
                uint32_t  coded_area_org = context_ptr->coded_area_sb;
                uint32_t  coded_area_org_uv = context_ptr->coded_area_sb_uv;

#if CHROMA_BLIND
                // Derive disable_cfl_flag as evaluate_cfl_ep = f(disable_cfl_flag)
                EbBool disable_cfl_flag = (context_ptr->blk_geom->sq_size > 32 ||
                    context_ptr->blk_geom->bwidth == 4 ||
                    context_ptr->blk_geom->bheight == 4) ? EB_TRUE : EB_FALSE;
                // Evaluate cfl @ EP if applicable, and not done @ MD 
                context_ptr->evaluate_cfl_ep = (disable_cfl_flag == EB_FALSE && context_ptr->md_context->chroma_level == CHROMA_MODE_1);
#endif

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
                cu_ptr->qp = (sequence_control_set_ptr->static_config.improve_sharpness) ? context_ptr->qpmQp : picture_control_set_ptr->picture_qp;
                sb_ptr->qp = (sequence_control_set_ptr->static_config.improve_sharpness) ? context_ptr->qpmQp : picture_control_set_ptr->picture_qp;
                cu_ptr->org_delta_qp = cu_ptr->delta_qp;
#endif
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
                        context_ptr->cu_stats->parent32x32Index, // TOCHECK not valid
                        context_ptr);
                }

#endif

                if (cu_ptr->prediction_mode_flag == INTRA_MODE) {

                    context_ptr->is_inter = 0;
                    context_ptr->tot_intra_coded_area += blk_geom->bwidth* blk_geom->bheight;
                    if (picture_control_set_ptr->slice_type != I_SLICE) {
                        context_ptr->intra_coded_area_sb[tbAddr] += blk_geom->bwidth* blk_geom->bheight;
                    }

                    // *Note - Transforms are the same size as predictions
                    // Partition Loop
                    context_ptr->txb_itr = 0;
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

                            //   if (picture_control_set_ptr->picture_number == 0 && context_ptr->cu_origin_x == 384 && context_ptr->cu_origin_y == 160)
                             //      printf("CHEDD");


#if INTRA_10BIT_SUPPORT

                            uint32_t cu_originy_uv = (context_ptr->cu_origin_y >> 3 << 3) >> 1;
                            uint32_t cu_originx_uv = (context_ptr->cu_origin_x >> 3 << 3) >> 1;
                            if (is16bit) {
                                uint16_t    topNeighArray[64 * 2 + 1];
                                uint16_t    leftNeighArray[64 * 2 + 1];
                                PredictionMode mode;

                                int32_t plane_end = blk_geom->has_uv ? 2 : 0;

                                for (int32_t plane = 0; plane <= plane_end; ++plane) {
                                    TxSize  tx_size = plane ? blk_geom->txsize_uv[context_ptr->txb_itr] : blk_geom->txsize[context_ptr->txb_itr];
                                    if (plane == 0) {
                                        if (context_ptr->cu_origin_y != 0)
                                            memcpy(topNeighArray + 1, (uint16_t*)(ep_luma_recon_neighbor_array->topArray) + context_ptr->cu_origin_x, blk_geom->bwidth * 2 * sizeof(uint16_t));
                                        if (context_ptr->cu_origin_x != 0)
                                            memcpy(leftNeighArray + 1, (uint16_t*)(ep_luma_recon_neighbor_array->leftArray) + context_ptr->cu_origin_y, blk_geom->bheight * 2 * sizeof(uint16_t));
                                        if (context_ptr->cu_origin_y != 0 && context_ptr->cu_origin_x != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_luma_recon_neighbor_array->topLeftArray) + MAX_PICTURE_HEIGHT_SIZE + context_ptr->cu_origin_x - context_ptr->cu_origin_y)[0];
                                    }

                                    else if (plane == 1) {
                                        if (cu_originy_uv != 0)
                                            memcpy(topNeighArray + 1, (uint16_t*)(ep_cb_recon_neighbor_array->topArray) + cu_originx_uv, blk_geom->bwidth_uv * 2 * sizeof(uint16_t));
                                        if (cu_originx_uv != 0)
                                            memcpy(leftNeighArray + 1, (uint16_t*)(ep_cb_recon_neighbor_array->leftArray) + cu_originy_uv, blk_geom->bheight_uv * 2 * sizeof(uint16_t));
                                        if (cu_originy_uv != 0 && cu_originx_uv != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_cb_recon_neighbor_array->topLeftArray) + MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv)[0];
                                    }
                                    else {
                                        if (cu_originy_uv != 0)
                                            memcpy(topNeighArray + 1, (uint16_t*)(ep_cr_recon_neighbor_array->topArray) + cu_originx_uv, blk_geom->bwidth_uv * 2 * sizeof(uint16_t));
                                        if (cu_originx_uv != 0)
                                            memcpy(leftNeighArray + 1, (uint16_t*)(ep_cr_recon_neighbor_array->leftArray) + cu_originy_uv, blk_geom->bheight_uv * 2 * sizeof(uint16_t));
                                        if (cu_originy_uv != 0 && cu_originx_uv != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(ep_cr_recon_neighbor_array->topLeftArray) + MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv)[0];

                                    }


                                    if (plane)
                                        mode = (pu_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode)UV_DC_PRED : (PredictionMode)pu_ptr->intra_chroma_mode;
                                    else
                                        mode = cu_ptr->pred_mode; //PredictionMode mode,

                                    av1_predict_intra_block_16bit(
#if TILES   
                                        &sb_ptr->tile_info,
#endif
                                        context_ptr,
                                        cu_ptr,
                                        picture_control_set_ptr->parent_pcs_ptr->av1_cm,                  //const Av1Common *cm,
                                        plane ? blk_geom->bwidth_uv : blk_geom->bwidth,                  //int32_t wpx,
                                        plane ? blk_geom->bheight_uv : blk_geom->bheight,                  //int32_t hpx,
                                        tx_size,
                                        mode,                                                       //PredictionMode mode,
                                        plane ? 0 : pu_ptr->angle_delta[PLANE_TYPE_Y],                //int32_t angle_delta,
                                        0,                                                          //int32_t use_palette,
                                        FILTER_INTRA_MODES,                                         //CHKN FILTER_INTRA_MODE filter_intra_mode,
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
                                    TxSize  tx_size = plane ? blk_geom->txsize_uv[context_ptr->txb_itr] : blk_geom->txsize[context_ptr->txb_itr];

                                    if (plane == 0) {
                                        if (context_ptr->cu_origin_y != 0)
                                            memcpy(topNeighArray + 1, ep_luma_recon_neighbor_array->topArray + context_ptr->cu_origin_x, blk_geom->bwidth * 2);

                                        if (context_ptr->cu_origin_x != 0)
                                            memcpy(leftNeighArray + 1, ep_luma_recon_neighbor_array->leftArray + context_ptr->cu_origin_y, blk_geom->bheight * 2);

                                        if (context_ptr->cu_origin_y != 0 && context_ptr->cu_origin_x != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ep_luma_recon_neighbor_array->topLeftArray[MAX_PICTURE_HEIGHT_SIZE + context_ptr->cu_origin_x - context_ptr->cu_origin_y];
                                    }

                                    else if (plane == 1) {
                                        if (cu_originy_uv != 0)
                                            memcpy(topNeighArray + 1, ep_cb_recon_neighbor_array->topArray + cu_originx_uv, blk_geom->bwidth_uv * 2);

                                        if (cu_originx_uv != 0)
                                            memcpy(leftNeighArray + 1, ep_cb_recon_neighbor_array->leftArray + cu_originy_uv, blk_geom->bheight_uv * 2);

                                        if (cu_originy_uv != 0 && cu_originx_uv != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ep_cb_recon_neighbor_array->topLeftArray[MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv];
                                    }
                                    else {
                                        if (cu_originy_uv != 0)
                                            memcpy(topNeighArray + 1, ep_cr_recon_neighbor_array->topArray + cu_originx_uv, blk_geom->bwidth_uv * 2);

                                        if (cu_originx_uv != 0)
                                            memcpy(leftNeighArray + 1, ep_cr_recon_neighbor_array->leftArray + cu_originy_uv, blk_geom->bheight_uv * 2);

                                        if (cu_originy_uv != 0 && cu_originx_uv != 0)
                                            topNeighArray[0] = leftNeighArray[0] = ep_cr_recon_neighbor_array->topLeftArray[MAX_PICTURE_HEIGHT_SIZE / 2 + cu_originx_uv - cu_originy_uv];
                                    }

                                    if (plane)
                                        mode = (pu_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode)UV_DC_PRED : (PredictionMode)pu_ptr->intra_chroma_mode;
                                    else
                                        mode = cu_ptr->pred_mode; //PredictionMode mode,

                                    // Hsan: if CHROMA_MODE_1, then CFL will be evaluated @ EP as no CHROMA @ MD 
                                    // If that's the case then you should ensure than the 1st chroma prediction uses UV_DC_PRED (that's the default configuration for CHROMA_MODE_1 if CFL applicable (set @ fast loop candidates injection) then MD assumes chroma mode always UV_DC_PRED)
                                    av1_predict_intra_block(
#if TILES
                                        &sb_ptr->tile_info,
#endif
#if INTRA_CORE_OPT
                                        NULL,
#endif
                                        ED_STAGE,
                                        cu_ptr->prediction_unit_array[0].intra_luma_left_mode,
                                        cu_ptr->prediction_unit_array[0].intra_luma_top_mode,
                                        cu_ptr->prediction_unit_array[0].intra_chroma_left_mode,
                                        cu_ptr->prediction_unit_array[0].intra_chroma_top_mode,
                                        context_ptr->blk_geom,
                                        picture_control_set_ptr->parent_pcs_ptr->av1_cm,                  //const Av1Common *cm,
                                        plane ? blk_geom->bwidth_uv : blk_geom->bwidth,                   //int32_t wpx,
                                        plane ? blk_geom->bheight_uv : blk_geom->bheight,                  //int32_t hpx,
                                        tx_size,
                                        mode,                                                       //PredictionMode mode,
                                        plane ? 0 : pu_ptr->angle_delta[PLANE_TYPE_Y],                //int32_t angle_delta,
                                        0,                                                          //int32_t use_palette,
                                        FILTER_INTRA_MODES,                                         //CHKN FILTER_INTRA_MODE filter_intra_mode,
                                        topNeighArray + 1,
                                        leftNeighArray + 1,
                                        recon_buffer,                                                //uint8_t *dst,
                                        //int32_t dst_stride,
#if !INTRA_CORE_OPT
                                        0,                                                          //int32_t col_off,
                                        0,                                                          //int32_t row_off,
#endif
                                        plane,                                                      //int32_t plane,
                                        blk_geom->bsize,                  //uint32_t puSize,
                                        context_ptr->cu_origin_x,
                                        context_ptr->cu_origin_y,
                                        0,  // MD ONLY - NOT USED BY ENCDEC
                                        0);
                                }
                            }

#else
                            uint8_t    topNeighArray[64 * 2 + 1];
                            uint8_t    leftNeighArray[64 * 2 + 1];
                            PredictionMode mode;

                            int32_t size = cu_stats->size * 2;

                            for (int32_t plane = 0; plane <= 2; ++plane) {
                                if (plane == 0) {

                                    if (context_ptr->cu_origin_y != 0)
                                        memcpy(topNeighArray + 1, ep_luma_recon_neighbor_array->topArray + context_ptr->cu_origin_x, size);
                                    if (context_ptr->cu_origin_x != 0)
                                        memcpy(leftNeighArray + 1, ep_luma_recon_neighbor_array->leftArray + context_ptr->cu_origin_y, size);
                                    if (context_ptr->cu_origin_y != 0 && context_ptr->cu_origin_x != 0)
                                        topNeighArray[0] = leftNeighArray[0] = ep_luma_recon_neighbor_array->topLeftArray[MAX_PICTURE_HEIGHT_SIZE + context_ptr->cu_origin_x - context_ptr->cu_origin_y];
                                }
                                else if (plane == 1) {
                                    if (context_ptr->cu_origin_y != 0)
                                        memcpy(topNeighArray + 1, ep_cb_recon_neighbor_array->topArray + context_ptr->cu_origin_x / 2, size / 2);
                                    if (context_ptr->cu_origin_x != 0)
                                        memcpy(leftNeighArray + 1, ep_cb_recon_neighbor_array->leftArray + context_ptr->cu_origin_y / 2, size / 2);
                                    if (context_ptr->cu_origin_y != 0 && context_ptr->cu_origin_x != 0)
                                        topNeighArray[0] = leftNeighArray[0] = ep_cb_recon_neighbor_array->topLeftArray[MAX_PICTURE_HEIGHT_SIZE / 2 + context_ptr->cu_origin_x / 2 - context_ptr->cu_origin_y / 2];
                                }
                                else {
                                    if (context_ptr->cu_origin_y != 0)
                                        memcpy(topNeighArray + 1, ep_cr_recon_neighbor_array->topArray + context_ptr->cu_origin_x / 2, size / 2);
                                    if (context_ptr->cu_origin_x != 0)
                                        memcpy(leftNeighArray + 1, ep_cr_recon_neighbor_array->leftArray + context_ptr->cu_origin_y / 2, size / 2);
                                    if (context_ptr->cu_origin_y != 0 && context_ptr->cu_origin_x != 0)
                                        topNeighArray[0] = leftNeighArray[0] = ep_cr_recon_neighbor_array->topLeftArray[MAX_PICTURE_HEIGHT_SIZE / 2 + context_ptr->cu_origin_x / 2 - context_ptr->cu_origin_y / 2];

                                }
                                if (plane)
                                    mode = (pu_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode)UV_DC_PRED : (PredictionMode)pu_ptr->intra_chroma_mode;
                                else
                                    mode = cu_ptr->pred_mode; //PredictionMode mode,

                                av1_predict_intra_block(
                                    context_ptr,
                                    cu_ptr,
                                    picture_control_set_ptr->parent_pcs_ptr->av1_cm,                  //const Av1Common *cm,
                                    plane ? cu_stats->size / 2 : cu_stats->size,                  //int32_t wpx,
                                    plane ? cu_stats->size / 2 : cu_stats->size,                  //int32_t hpx,
                                    plane ? tx_size_Chroma : tx_size,                           //TxSize tx_size,
                                    mode,                                                       //PredictionMode mode,
                                    plane ? 0 : pu_ptr->angle_delta[PLANE_TYPE_Y],                //int32_t angle_delta,
                                    0,                                                          //int32_t use_palette,
                                    FILTER_INTRA_MODES,                                         //CHKN FILTER_INTRA_MODE filter_intra_mode,
                                    topNeighArray + 1,
                                    leftNeighArray + 1,
                                    recon_buffer,                                                //uint8_t *dst,
                                                                                                //int32_t dst_stride,
                                    0,                                                          //int32_t col_off,
                                    0,                                                          //int32_t row_off,
                                    plane,                                                      //int32_t plane,
                                    plane ? cu_stats->size / 2 : cu_stats->size,                  //uint32_t puSize,
                                    plane ? context_ptr->cu_origin_x / 2 : context_ptr->cu_origin_x,  //uint32_t cuOrgX,
                                    plane ? context_ptr->cu_origin_y / 2 : context_ptr->cu_origin_y   //uint32_t cuOrgY
                                );
                            }
#endif

                            // Encode Transform Unit -INTRA-
                            {

                                uint8_t             cbQp = cu_ptr->qp;


                                Av1EncodeLoopFunctionTable[is16bit](
#if ENCDEC_TX_SEARCH
                                    picture_control_set_ptr,
#endif
                                    context_ptr,
                                    sb_ptr,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    cbQp,
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

                                if (cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff) {
                                    cu_ptr->transform_unit_array[0].u_has_coeff = EB_TRUE;
                                }

                                if (cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff) {
                                    cu_ptr->transform_unit_array[0].v_has_coeff = EB_TRUE;
                                }

                            }
                            else {
                                cu_ptr->block_has_coeff = cu_ptr->block_has_coeff |
                                    cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff;
                            }


                        } // Transform Loop

                    } // Partition Loop

                    context_ptr->coded_area_sb += blk_geom->bwidth * blk_geom->bheight;
                    if (blk_geom->has_uv)
                        context_ptr->coded_area_sb_uv += blk_geom->bwidth_uv * blk_geom->bheight_uv;
                }

                // Inter
                else if (cu_ptr->prediction_mode_flag == INTER_MODE) {

#if ENCDEC_TX_SEARCH
                    context_ptr->is_inter = 1;
#endif

                    EbReferenceObject_t* refObj0 = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;
                    EbReferenceObject_t* refObj1 = picture_control_set_ptr->slice_type == B_SLICE ?
                        (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr : 0;

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
                        {
                            isCuSkip = mdcontextPtr->md_ep_pipe_sb[cu_ptr->mds_idx].skip_cost <= mdcontextPtr->md_ep_pipe_sb[cu_ptr->mds_idx].merge_cost ? 1 : 0;
                        }
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
#if TILES
                        &sb_ptr->tile_info,
#endif
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
                        context_ptr->mv_unit.predDirection = (uint8_t)pu_ptr->inter_pred_direction_index;
                        context_ptr->mv_unit.mv[REF_LIST_0].mvUnion = pu_ptr->mv[REF_LIST_0].mvUnion;
                        context_ptr->mv_unit.mv[REF_LIST_1].mvUnion = pu_ptr->mv[REF_LIST_1].mvUnion;

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
                                is16bit ? refObj0->referencePicture16bit : refObj0->referencePicture,
                                recon_buffer,
                                context_ptr->cu_origin_x,
                                context_ptr->cu_origin_y,
                                &cu_ptr->prediction_unit_array[0].wm_params,
                                (uint8_t) sequence_control_set_ptr->static_config.encoder_bit_depth,
#if CHROMA_BLIND
                                EB_TRUE,
#endif
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
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    blk_geom->bwidth,
                                    blk_geom->bheight,
                                    refObj0->referencePicture16bit,
                                    picture_control_set_ptr->slice_type == B_SLICE ? refObj1->referencePicture16bit : 0,
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
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
                                    blk_geom->bwidth,
                                    blk_geom->bheight,
                                    refObj0->referencePicture,
                                    picture_control_set_ptr->slice_type == B_SLICE ? refObj1->referencePicture : 0,
                                    recon_buffer,
                                    context_ptr->cu_origin_x,
                                    context_ptr->cu_origin_y,
#if CHROMA_BLIND
                                    EB_TRUE,
#endif
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


                    uint32_t totTu = context_ptr->blk_geom->txb_count;
                    uint8_t   tuIt;
                    uint8_t    cbQp = cu_ptr->qp;
                    uint32_t  component_mask = context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK;

                    if (cu_ptr->prediction_unit_array[0].merge_flag == EB_FALSE) {

                        for (tuIt = 0; tuIt < totTu; tuIt++) {
                            context_ptr->txb_itr = tuIt;
                            txb_origin_x = context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[tuIt];
                            txb_origin_y = context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[tuIt];
                            if (!zeroLumaCbfMD)
                                //inter mode  1
                                Av1EncodeLoopFunctionTable[is16bit](
#if ENCDEC_TX_SEARCH
                                    picture_control_set_ptr,
#endif
                                    context_ptr,
                                    sb_ptr,
                                    txb_origin_x,   //pic org
                                    txb_origin_y,
                                    cbQp,
                                    recon_buffer,
                                    coeff_buffer_sb,
                                    residual_buffer,
                                    transform_buffer,
                                    inverse_quant_buffer,
                                    transform_inner_array_ptr,
                                    asm_type,
                                    count_non_zero_coeffs,
                                    context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                                    useDeltaQpSegments,
                                    cu_ptr->delta_qp > 0 ? 0 : dZoffset,
                                    eobs[context_ptr->txb_itr],
                                    cuPlane);

                            // SKIP the CBF zero mode for DC path. There are problems with cost calculations
                            if (context_ptr->trans_coeff_shape_luma != ONLY_DC_SHAPE) {
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
                                        blk_geom->tx_width[tuIt],
                                        blk_geom->tx_height[tuIt],
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
                                TxSize  txSize = blk_geom->txsize[context_ptr->txb_itr];
                                int32_t shift = (MAX_TX_SCALE - av1_get_tx_scale(txSize)) * 2;
                                yTuFullDistortion[DIST_CALC_RESIDUAL] = RIGHT_SIGNED_SHIFT(yTuFullDistortion[DIST_CALC_RESIDUAL], shift);
                                yTuFullDistortion[DIST_CALC_PREDICTION] = RIGHT_SIGNED_SHIFT(yTuFullDistortion[DIST_CALC_PREDICTION], shift);

                                y_tu_coeff_bits = 0;
                                cb_tu_coeff_bits = 0;
                                cr_tu_coeff_bits = 0;

                                if (!zeroLumaCbfMD) {

                                    ModeDecisionCandidateBuffer_t         **candidateBufferPtrArrayBase = context_ptr->md_context->candidate_buffer_ptr_array;
                                    ModeDecisionCandidateBuffer_t         **candidate_buffer_ptr_array = &(candidateBufferPtrArrayBase[context_ptr->md_context->buffer_depth_index_start[0]]);
                                    ModeDecisionCandidateBuffer_t          *candidateBuffer;

                                    // Set the Candidate Buffer
                                    candidateBuffer = candidate_buffer_ptr_array[0];
                                    // Rate estimation function uses the values from CandidatePtr. The right values are copied from cu_ptr to CandidatePtr
                                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_Y] = cu_ptr->transform_unit_array[tuIt].transform_type[PLANE_TYPE_Y];
                                    candidateBuffer->candidate_ptr->transform_type[PLANE_TYPE_UV] = cu_ptr->transform_unit_array[tuIt].transform_type[PLANE_TYPE_UV];
                                    candidateBuffer->candidate_ptr->type = cu_ptr->prediction_mode_flag;

                                    const uint32_t coeff1dOffset = context_ptr->coded_area_sb;

                                    Av1TuEstimateCoeffBits(
                                        picture_control_set_ptr,
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
                                        context_ptr->blk_geom->txsize[context_ptr->txb_itr],
                                        context_ptr->blk_geom->txsize_uv[context_ptr->txb_itr],
                                        context_ptr->blk_geom->has_uv ? COMPONENT_ALL : COMPONENT_LUMA,
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
                                if (context_ptr->blk_geom->has_uv) {
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

                                if (context_ptr->blk_geom->has_uv) {
                                    cb_coeff_bits += cb_tu_coeff_bits;
                                    cr_coeff_bits += cr_tu_coeff_bits;
                                }

                                y_full_distortion[DIST_CALC_RESIDUAL] += yTuFullDistortion[DIST_CALC_RESIDUAL];
                                y_full_distortion[DIST_CALC_PREDICTION] += yTuFullDistortion[DIST_CALC_PREDICTION];

                            }
                            context_ptr->coded_area_sb += blk_geom->tx_width[tuIt] * blk_geom->tx_height[tuIt];
                            if (blk_geom->has_uv)
                                context_ptr->coded_area_sb_uv += blk_geom->tx_width_uv[tuIt] * blk_geom->tx_height_uv[tuIt];

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
                    totTu = context_ptr->blk_geom->txb_count;


                    //reset coeff buffer offsets at the start of a new Tx loop
                    context_ptr->coded_area_sb = coded_area_org;
                    context_ptr->coded_area_sb_uv = coded_area_org_uv;
                    for (tuIt = 0; tuIt < totTu; tuIt++)
                    {
                        context_ptr->txb_itr = tuIt;
                        txb_origin_x = context_ptr->cu_origin_x + context_ptr->blk_geom->tx_boff_x[tuIt];
                        txb_origin_y = context_ptr->cu_origin_y + context_ptr->blk_geom->tx_boff_y[tuIt];
                        if (cu_ptr->skip_flag == EB_TRUE) {
                            cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff = EB_FALSE;
                            cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff = EB_FALSE;
                            cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff = EB_FALSE;


                        }
                        else if ((&cu_ptr->prediction_unit_array[0])->merge_flag == EB_TRUE) {

                            //inter mode  2

                            Av1EncodeLoopFunctionTable[is16bit](
#if ENCDEC_TX_SEARCH
                                picture_control_set_ptr,
#endif
                                context_ptr,
                                sb_ptr,
                                txb_origin_x, //pic offset
                                txb_origin_y,
                                cbQp,
                                recon_buffer,
                                coeff_buffer_sb,
                                residual_buffer,
                                transform_buffer,
                                inverse_quant_buffer,
                                transform_inner_array_ptr,
                                asm_type,
                                count_non_zero_coeffs,
                                context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                                useDeltaQpSegments,
                                cu_ptr->delta_qp > 0 ? 0 : dZoffset,
                                eobs[context_ptr->txb_itr],
                                cuPlane);



                        }

                        if (context_ptr->blk_geom->has_uv) {
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
                                context_ptr->blk_geom->has_uv ? PICTURE_BUFFER_DESC_FULL_MASK : PICTURE_BUFFER_DESC_LUMA_MASK,
                                eobs[context_ptr->txb_itr],
                                asm_type);
                        if (context_ptr->blk_geom->has_uv) {
                            y_has_coeff |= cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff;
                            u_has_coeff |= cu_ptr->transform_unit_array[context_ptr->txb_itr].u_has_coeff;
                            v_has_coeff |= cu_ptr->transform_unit_array[context_ptr->txb_itr].v_has_coeff;
                        }
                        else {
                            y_has_coeff |= cu_ptr->transform_unit_array[context_ptr->txb_itr].y_has_coeff;
                        }


                        context_ptr->coded_area_sb += blk_geom->tx_width[tuIt] * blk_geom->tx_height[tuIt];
                        if (blk_geom->has_uv)
                            context_ptr->coded_area_sb_uv += blk_geom->tx_width_uv[tuIt] * blk_geom->tx_height_uv[tuIt];


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
                        context_ptr->mv_unit.predDirection = (uint8_t)pu_ptr->inter_pred_direction_index;
                        context_ptr->mv_unit.mv[REF_LIST_0].mvUnion = pu_ptr->mv[REF_LIST_0].mvUnion;
                        context_ptr->mv_unit.mv[REF_LIST_1].mvUnion = pu_ptr->mv[REF_LIST_1].mvUnion;

                        // Update Neighbor Arrays (Mode Type, MVs, SKIP)
                        {
                            uint8_t skip_flag = (uint8_t)cu_ptr->skip_flag;
                            EncodePassUpdateInterModeNeighborArrays(
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

                    if (blk_geom->has_uv) {
                        availableCoeff = (cu_ptr->prediction_mode_flag == INTER_MODE) ? (EbBool)cu_ptr->block_has_coeff :
                            (cu_ptr->transform_unit_array[0].y_has_coeff ||
                                cu_ptr->transform_unit_array[0].v_has_coeff ||
                                cu_ptr->transform_unit_array[0].u_has_coeff) ? EB_TRUE : EB_FALSE;
                    }
                    else {
                        availableCoeff = (cu_ptr->transform_unit_array[0].y_has_coeff) ? EB_TRUE : EB_FALSE;
                    }


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

                }

                {
                    {
                        // Set the PU Loop Variables
                        pu_ptr = cu_ptr->prediction_unit_array;
                        // Set MvUnit
                        context_ptr->mv_unit.predDirection = (uint8_t)pu_ptr->inter_pred_direction_index;
                        context_ptr->mv_unit.mv[REF_LIST_0].mvUnion = pu_ptr->mv[REF_LIST_0].mvUnion;
                        context_ptr->mv_unit.mv[REF_LIST_1].mvUnion = pu_ptr->mv[REF_LIST_1].mvUnion;
                    }

                }


                {

                    CodingUnit_t *src_cu = &context_ptr->md_context->md_cu_arr_nsq[d1_itr];

                    CodingUnit_t *dst_cu = &sb_ptr->final_cu_arr[final_cu_itr++];

                    move_cu_data(src_cu, dst_cu);
                }

            }
            blk_it += ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][context_ptr->blk_geom->depth];
        }
        else {
            blk_it += d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][context_ptr->blk_geom->depth];

        }


    } // CU Loop

    sb_ptr->tot_final_cu = final_cu_itr;
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
    SequenceControlSet_t    *sequence_control_set_ptr,
    PictureControlSet_t     *picture_control_set_ptr,
    LargestCodingUnit_t     *sb_ptr,
    uint32_t                   tbAddr,
    uint32_t                   sb_origin_x,
    uint32_t                   sb_origin_y,
    uint32_t                   sb_qp,
    EncDecContext_t         *context_ptr)
{

    context_ptr->coded_area_sb = 0;
    context_ptr->coded_area_sb_uv = 0;

    uint32_t      final_cu_itr = 0;


    uint32_t    blk_it = 0;

    while (blk_it < sequence_control_set_ptr->max_block_cnt) {


        CodingUnit_t  *cu_ptr = context_ptr->cu_ptr = &context_ptr->md_context->md_cu_arr_nsq[blk_it];
        PartitionType part = cu_ptr->part;
        const BlockGeom * blk_geom = context_ptr->blk_geom = get_blk_geom_mds(blk_it);


        sb_ptr->cu_partition_array[blk_it] = context_ptr->md_context->md_cu_arr_nsq[blk_it].part;

        if (part != PARTITION_SPLIT) {



            int32_t offset_d1 = ns_blk_offset[(int32_t)part]; //cu_ptr->best_d1_blk; // TOCKECK
            int32_t num_d1_block = ns_blk_num[(int32_t)part]; // context_ptr->blk_geom->totns; // TOCKECK

            for (int32_t d1_itr = blk_it + offset_d1; d1_itr < blk_it + offset_d1 + num_d1_block; d1_itr++) {

                const BlockGeom * blk_geom = context_ptr->blk_geom = get_blk_geom_mds(d1_itr);
                CodingUnit_t            *cu_ptr = context_ptr->cu_ptr = &context_ptr->md_context->md_cu_arr_nsq[d1_itr];


                cu_ptr->delta_qp = 0;
                cu_ptr->qp = (sequence_control_set_ptr->static_config.improve_sharpness) ? context_ptr->qpmQp : picture_control_set_ptr->picture_qp;
                sb_ptr->qp = (sequence_control_set_ptr->static_config.improve_sharpness) ? context_ptr->qpmQp : picture_control_set_ptr->picture_qp;
                cu_ptr->org_delta_qp = cu_ptr->delta_qp;


                {
                    CodingUnit_t *src_cu = &context_ptr->md_context->md_cu_arr_nsq[d1_itr];
                    CodingUnit_t *dst_cu = &sb_ptr->final_cu_arr[final_cu_itr++];

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
                    {
                        memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                    }

                    if (context_ptr->blk_geom->has_uv)
                    {
                        // Cb
                        bwidth = context_ptr->blk_geom->tx_width_uv[txb_itr];
                        bheight = context_ptr->blk_geom->tx_height_uv[txb_itr];

                        src_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->bufferCb)[txb_1d_offset_uv]);
                        dst_ptr = &(((int32_t*)sb_ptr->quantized_coeff->bufferCb)[context_ptr->coded_area_sb_uv]);

                        for (j = 0; j < bheight; j++)
                        {
                            memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                        }

                        //Cr
                        src_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->bufferCr)[txb_1d_offset_uv]);
                        dst_ptr = &(((int32_t*)sb_ptr->quantized_coeff->bufferCr)[context_ptr->coded_area_sb_uv]);

                        for (j = 0; j < bheight; j++)
                        {
                            memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                        }

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
                    EbPictureBufferDesc_t          *ref_pic;
                    if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                    {
                        EbReferenceObject_t* refObj = (EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
                        ref_pic = refObj->referencePicture;
                    }
                    else
                    {
                        ref_pic = picture_control_set_ptr->recon_picture_ptr;
                    }

                    context_ptr->cu_origin_x = sb_origin_x + context_ptr->blk_geom->origin_x;
                    context_ptr->cu_origin_y = sb_origin_y + context_ptr->blk_geom->origin_y;

                    uint32_t  bwidth = context_ptr->blk_geom->bwidth;
                    uint32_t  bheight = context_ptr->blk_geom->bheight;

                    uint8_t* src_ptr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->buffer_y)[0]);
                    uint8_t* dst_ptr = ref_pic->buffer_y + ref_pic->origin_x + context_ptr->cu_origin_x + (ref_pic->origin_y + context_ptr->cu_origin_y)*ref_pic->stride_y;

                    uint32_t j;
                    for (j = 0; j < bheight; j++)
                    {
                        memcpy(dst_ptr + j * ref_pic->stride_y, src_ptr + j * 128, bwidth * sizeof(uint8_t));
                    }

                    if (context_ptr->blk_geom->has_uv)
                    {

                        bwidth = context_ptr->blk_geom->bwidth_uv;
                        bheight = context_ptr->blk_geom->bheight_uv;

                        src_ptr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->bufferCb)[0]);

                        dst_ptr = ref_pic->bufferCb + ref_pic->origin_x / 2 + ((context_ptr->cu_origin_x >> 3) << 3) / 2 + (ref_pic->origin_y / 2 + ((context_ptr->cu_origin_y >> 3) << 3) / 2)*ref_pic->strideCb;

                        for (j = 0; j < bheight; j++)
                        {
                            memcpy(dst_ptr + j * ref_pic->strideCb, src_ptr + j * 64, bwidth * sizeof(uint8_t));
                        }

                        src_ptr = &(((uint8_t*)context_ptr->cu_ptr->recon_tmp->bufferCr)[0]);

                        dst_ptr = ref_pic->bufferCr + ref_pic->origin_x / 2 + ((context_ptr->cu_origin_x >> 3) << 3) / 2 + (ref_pic->origin_y / 2 + ((context_ptr->cu_origin_y >> 3) << 3) / 2)*ref_pic->strideCr;


                        for (j = 0; j < bheight; j++)
                        {
                            memcpy(dst_ptr + j * ref_pic->strideCr, src_ptr + j * 64, bwidth * sizeof(uint8_t));
                        }

                    }

                }



            }
            blk_it += ns_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][context_ptr->blk_geom->depth];
        }
        else
        {
            blk_it += d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][context_ptr->blk_geom->depth];
        }

    } // CU Loop



    return;
}
#endif
