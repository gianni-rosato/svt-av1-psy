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

#include "EbEntropyCodingProcess.h"
#include "EbEncDecResults.h"
#include "EbEntropyCodingResults.h"
#include "EbRateControlTasks.h"

#if TILES
#define  AV1_MIN_TILE_SIZE_BYTES 1
void av1_reset_loop_restoration(PictureControlSet_t     *piCSetPtr);
void av1_tile_set_col(TileInfo *tile, PictureParentControlSet_t * pcsPtr, int col);
void av1_tile_set_row(TileInfo *tile, PictureParentControlSet_t * pcsPtr, int row);
#endif

/******************************************************
 * Enc Dec Context Constructor
 ******************************************************/
EbErrorType entropy_coding_context_ctor(
    EntropyCodingContext_t **context_dbl_ptr,
    EbFifo_t                *enc_dec_input_fifo_ptr,
    EbFifo_t                *packetization_output_fifo_ptr,
    EbFifo_t                *rate_control_output_fifo_ptr,
    EbBool                  is16bit)
{
    EntropyCodingContext_t *context_ptr;
    EB_MALLOC(EntropyCodingContext_t*, context_ptr, sizeof(EntropyCodingContext_t), EB_N_PTR);
    *context_dbl_ptr = context_ptr;

    context_ptr->is16bit = is16bit;

    // Input/Output System Resource Manager FIFOs
    context_ptr->enc_dec_input_fifo_ptr = enc_dec_input_fifo_ptr;
    context_ptr->entropy_coding_output_fifo_ptr = packetization_output_fifo_ptr;
    context_ptr->rate_control_output_fifo_ptr = rate_control_output_fifo_ptr;

    return EB_ErrorNone;
}

/***********************************************
 * Entropy Coding Reset Neighbor Arrays
 ***********************************************/
static void EntropyCodingResetNeighborArrays(PictureControlSet_t *picture_control_set_ptr)
{
    neighbor_array_unit_reset(picture_control_set_ptr->mode_type_neighbor_array);

    neighbor_array_unit_reset(picture_control_set_ptr->partition_context_neighbor_array);

    neighbor_array_unit_reset(picture_control_set_ptr->skip_flag_neighbor_array);

    neighbor_array_unit_reset(picture_control_set_ptr->skip_coeff_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->luma_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->cb_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->cr_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->inter_pred_dir_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ref_frame_type_neighbor_array);

    neighbor_array_unit_reset(picture_control_set_ptr->intra_luma_mode_neighbor_array);
    neighbor_array_unit_reset32(picture_control_set_ptr->interpolation_type_neighbor_array);
    return;
}

void av1_get_syntax_rate_from_cdf(
    int32_t                      *costs,
    const aom_cdf_prob       *cdf,
    const int32_t                *inv_map);

void av1_cost_tokens_from_cdf(int32_t *costs, const aom_cdf_prob *cdf,
    const int32_t *inv_map) {
    // int32_t i;
    // aom_cdf_prob prev_cdf = 0;
    // for (i = 0;; ++i) {
    //     aom_cdf_prob p15 = AOM_ICDF(cdf[i]) - prev_cdf;
    //     p15 = (p15 < EC_MIN_PROB) ? EC_MIN_PROB : p15;
    //     prev_cdf = AOM_ICDF(cdf[i]);
    //
    //     if (inv_map)
    //         costs[inv_map[i]] = av1_cost_symbol(p15);
    //     else
    //         costs[i] = av1_cost_symbol(p15);
    //
    //     // Stop once we reach the end of the CDF
    //     if (cdf[i] == AOM_ICDF(CDF_PROB_TOP)) break;
    // }

    av1_get_syntax_rate_from_cdf(costs, cdf, inv_map);

}

static void build_nmv_component_cost_table(int32_t *mvcost,
    const nmv_component *const mvcomp,
    MvSubpelPrecision precision) {
    int32_t i, v;
    int32_t sign_cost[2], class_cost[MV_CLASSES], class0_cost[CLASS0_SIZE];
    int32_t bits_cost[MV_OFFSET_BITS][2];
    int32_t class0_fp_cost[CLASS0_SIZE][MV_FP_SIZE], fp_cost[MV_FP_SIZE];
    int32_t class0_hp_cost[2], hp_cost[2];

    av1_cost_tokens_from_cdf(sign_cost, mvcomp->sign_cdf, NULL);
    av1_cost_tokens_from_cdf(class_cost, mvcomp->classes_cdf, NULL);
    av1_cost_tokens_from_cdf(class0_cost, mvcomp->class0_cdf, NULL);
    for (i = 0; i < MV_OFFSET_BITS; ++i) {
        av1_cost_tokens_from_cdf(bits_cost[i], mvcomp->bits_cdf[i], NULL);
    }

    for (i = 0; i < CLASS0_SIZE; ++i)
        av1_cost_tokens_from_cdf(class0_fp_cost[i], mvcomp->class0_fp_cdf[i], NULL);
    av1_cost_tokens_from_cdf(fp_cost, mvcomp->fp_cdf, NULL);

    if (precision > MV_SUBPEL_LOW_PRECISION) {
        av1_cost_tokens_from_cdf(class0_hp_cost, mvcomp->class0_hp_cdf, NULL);
        av1_cost_tokens_from_cdf(hp_cost, mvcomp->hp_cdf, NULL);
    }
    mvcost[0] = 0;
    for (v = 1; v <= MV_MAX; ++v) {
        int32_t z, c, o, d, e, f, cost = 0;
        z = v - 1;
        c = av1_get_mv_class(z, &o);
        cost += class_cost[c];
        d = (o >> 3);     /* int32_t mv data */
        f = (o >> 1) & 3; /* fractional pel mv data */
        e = (o & 1);      /* high precision mv data */
        if (c == MV_CLASS_0) {
            cost += class0_cost[d];
        }
        else {
            const int32_t b = c + CLASS0_BITS - 1; /* number of bits */
            for (i = 0; i < b; ++i) cost += bits_cost[i][((d >> i) & 1)];
        }
        if (precision > MV_SUBPEL_NONE) {
            if (c == MV_CLASS_0) {
                cost += class0_fp_cost[d][f];
            }
            else {
                cost += fp_cost[f];
            }
            if (precision > MV_SUBPEL_LOW_PRECISION) {
                if (c == MV_CLASS_0) {
                    cost += class0_hp_cost[e];
                }
                else {
                    cost += hp_cost[e];
                }
            }
        }
        mvcost[v] = cost + sign_cost[0];
        mvcost[-v] = cost + sign_cost[1];
    }
}
void av1_build_nmv_cost_table(int32_t *mvjoint, int32_t *mvcost[2],
    const nmv_context *ctx,
    MvSubpelPrecision precision) {
    av1_cost_tokens_from_cdf(mvjoint, ctx->joints_cdf, NULL);
    build_nmv_component_cost_table(mvcost[0], &ctx->comps[0], precision);
    build_nmv_component_cost_table(mvcost[1], &ctx->comps[1], precision);
}


/**************************************************
 * Reset Entropy Coding Picture
 **************************************************/
static void ResetEntropyCodingPicture(
    EntropyCodingContext_t  *context_ptr,
    PictureControlSet_t     *picture_control_set_ptr,
    SequenceControlSet_t    *sequence_control_set_ptr)
{
    ResetBitstream(EntropyCoderGetBitstreamPtr(picture_control_set_ptr->entropy_coder_ptr));

    uint32_t                       entropyCodingQp;

    context_ptr->is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    // QP
#if ADD_DELTA_QP_SUPPORT
    uint16_t picture_qp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    context_ptr->qp = picture_qp;
#else
    context_ptr->qp = picture_control_set_ptr->picture_qp;
#endif
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation

    context_ptr->chroma_qp = context_ptr->qp;
    if (picture_control_set_ptr->use_delta_qp) {
        entropyCodingQp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    }
    else {
        entropyCodingQp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    }


    // Reset CABAC Contexts
    // Reset QP Assignement
    picture_control_set_ptr->prev_coded_qp = picture_control_set_ptr->picture_qp;
    picture_control_set_ptr->prev_quant_group_coded_qp = picture_control_set_ptr->picture_qp;

#if ADD_DELTA_QP_SUPPORT //PART 0
    picture_control_set_ptr->parent_pcs_ptr->prev_qindex = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    if (picture_control_set_ptr->parent_pcs_ptr->allow_intrabc)
        assert(picture_control_set_ptr->parent_pcs_ptr->delta_lf_present_flag == 0);
    /*else
        aom_wb_write_bit(wb, pcsPtr->delta_lf_present_flag);*/
    if (picture_control_set_ptr->parent_pcs_ptr->delta_lf_present_flag) {
        //aom_wb_write_literal(wb, OD_ILOG_NZ(pcsPtr->delta_lf_res) - 1, 2);
        picture_control_set_ptr->parent_pcs_ptr->prev_delta_lf_from_base = 0;
        //aom_wb_write_bit(wb, pcsPtr->delta_lf_multi);
        const int32_t frame_lf_count =
            picture_control_set_ptr->parent_pcs_ptr->monochrome == 0 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
        for (int32_t lf_id = 0; lf_id < frame_lf_count; ++lf_id)
            picture_control_set_ptr->parent_pcs_ptr->prev_delta_lf[lf_id] = 0;
    }
#endif

    // pass the ent
    OutputBitstreamUnit_t *outputBitstreamPtr = (OutputBitstreamUnit_t*)(picture_control_set_ptr->entropy_coder_ptr->ecOutputBitstreamPtr);
    //****************************************************************//

    uint8_t *data = outputBitstreamPtr->bufferAv1;
    picture_control_set_ptr->entropy_coder_ptr->ecWriter.allow_update_cdf = !picture_control_set_ptr->parent_pcs_ptr->large_scale_tile;
    picture_control_set_ptr->entropy_coder_ptr->ecWriter.allow_update_cdf =
        picture_control_set_ptr->entropy_coder_ptr->ecWriter.allow_update_cdf && !picture_control_set_ptr->parent_pcs_ptr->disable_cdf_update;
    aom_start_encode(&picture_control_set_ptr->entropy_coder_ptr->ecWriter, data);

    // ADD Reset here

    ResetEntropyCoder(
        sequence_control_set_ptr->encode_context_ptr,
        picture_control_set_ptr->entropy_coder_ptr,
        entropyCodingQp,
        picture_control_set_ptr->slice_type);

    EntropyCodingResetNeighborArrays(picture_control_set_ptr);


    return;
}


#if TILES
static void reset_ec_tile(
    uint32_t  total_size,
    uint32_t  is_last_tile_in_tg,
    EntropyCodingContext_t  *context_ptr,
    PictureControlSet_t     *picture_control_set_ptr,
    SequenceControlSet_t    *sequence_control_set_ptr)
{
    ResetBitstream(EntropyCoderGetBitstreamPtr(picture_control_set_ptr->entropy_coder_ptr));

    uint32_t                       entropy_coding_qp;

    context_ptr->is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    // QP
#if ADD_DELTA_QP_SUPPORT
    uint16_t picture_qp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    context_ptr->qp = picture_qp;
#else
    context_ptr->qp = picture_control_set_ptr->picture_qp;
#endif
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation

    context_ptr->chroma_qp = context_ptr->qp;
    if (picture_control_set_ptr->use_delta_qp) {
        entropy_coding_qp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    }
    else {
        entropy_coding_qp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    }
    
    // Reset CABAC Contexts
    // Reset QP Assignement
    picture_control_set_ptr->prev_coded_qp = picture_control_set_ptr->picture_qp;
    picture_control_set_ptr->prev_quant_group_coded_qp = picture_control_set_ptr->picture_qp;

#if ADD_DELTA_QP_SUPPORT //PART 0
    picture_control_set_ptr->parent_pcs_ptr->prev_qindex = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    if (picture_control_set_ptr->parent_pcs_ptr->allow_intrabc)
        assert(picture_control_set_ptr->parent_pcs_ptr->delta_lf_present_flag == 0);
    /*else
        aom_wb_write_bit(wb, pcsPtr->delta_lf_present_flag);*/
    if (picture_control_set_ptr->parent_pcs_ptr->delta_lf_present_flag) {
        //aom_wb_write_literal(wb, OD_ILOG_NZ(pcsPtr->delta_lf_res) - 1, 2);
        picture_control_set_ptr->parent_pcs_ptr->prev_delta_lf_from_base = 0;
        //aom_wb_write_bit(wb, pcsPtr->delta_lf_multi);
        const int32_t frame_lf_count =
            picture_control_set_ptr->parent_pcs_ptr->monochrome == 0 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
        for (int32_t lf_id = 0; lf_id < frame_lf_count; ++lf_id)
            picture_control_set_ptr->parent_pcs_ptr->prev_delta_lf[lf_id] = 0;
    }
#endif

    // pass the ent
    OutputBitstreamUnit_t *outputBitstreamPtr = (OutputBitstreamUnit_t*)(picture_control_set_ptr->entropy_coder_ptr->ecOutputBitstreamPtr);
    //****************************************************************//

    uint8_t *data = outputBitstreamPtr->bufferAv1 + total_size;
    picture_control_set_ptr->entropy_coder_ptr->ecWriter.allow_update_cdf = !picture_control_set_ptr->parent_pcs_ptr->large_scale_tile;
    picture_control_set_ptr->entropy_coder_ptr->ecWriter.allow_update_cdf =
        picture_control_set_ptr->entropy_coder_ptr->ecWriter.allow_update_cdf && !picture_control_set_ptr->parent_pcs_ptr->disable_cdf_update;


    //if not last tile, advance buffer by 4B to leave space for tile Size
    if (is_last_tile_in_tg == 0)
        data += 4;

    aom_start_encode(&picture_control_set_ptr->entropy_coder_ptr->ecWriter, data);

    //reset probabilities
    ResetEntropyCoder(
        sequence_control_set_ptr->encode_context_ptr,
        picture_control_set_ptr->entropy_coder_ptr,
        entropy_coding_qp,
        picture_control_set_ptr->slice_type);

    EntropyCodingResetNeighborArrays(picture_control_set_ptr);


    return;
}
#endif
/******************************************************
 * EncDec Configure LCU
 ******************************************************/
static void EntropyCodingConfigureLcu(
    EntropyCodingContext_t  *context_ptr,
    LargestCodingUnit_t     *sb_ptr,
    PictureControlSet_t     *picture_control_set_ptr)
{
#if ADD_DELTA_QP_SUPPORT
    context_ptr->qp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
#else
    context_ptr->qp = picture_control_set_ptr->picture_qp;
#endif

    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation

    context_ptr->chroma_qp = context_ptr->qp;

    sb_ptr->qp = context_ptr->qp;

    return;
}

/******************************************************
 * Entropy Coding Lcu
 ******************************************************/
static void EntropyCodingLcu(
    EntropyCodingContext_t              *context_ptr,
    LargestCodingUnit_t               *sb_ptr,
    PictureControlSet_t               *picture_control_set_ptr,
    SequenceControlSet_t              *sequence_control_set_ptr,
    uint32_t                             sb_origin_x,
    uint32_t                             sb_origin_y,
    EbBool                            terminateSliceFlag,
    uint32_t                             pictureOriginX,
    uint32_t                             pictureOriginY)
{
    UNUSED(sb_origin_x);
    UNUSED(sb_origin_y);
    (void)terminateSliceFlag;
    (void)sequence_control_set_ptr;
    EbPictureBufferDesc_t *coeffPicturePtr = sb_ptr->quantized_coeff;

    //rate Control
    uint32_t                       writtenBitsBeforeQuantizedCoeff;
    uint32_t                       writtenBitsAfterQuantizedCoeff;

    //store the number of written bits before coding quantized coeffs (flush is not called yet):
    // The total number of bits is
    // number of written bits
    // + 32  - bits remaining in interval Low value
    // + number of buffered byte * 8
    // This should be only for coeffs not any flag
    writtenBitsBeforeQuantizedCoeff = ((OutputBitstreamUnit_t*)EntropyCoderGetBitstreamPtr(picture_control_set_ptr->entropy_coder_ptr))->writtenBitsCount;

    (void)pictureOriginX;
    (void)pictureOriginY;

    write_sb(
        context_ptr,
        sb_ptr,
        picture_control_set_ptr,
        picture_control_set_ptr->entropy_coder_ptr,
        coeffPicturePtr);

    //store the number of written bits after coding quantized coeffs (flush is not called yet):
    // The total number of bits is
    // number of written bits
    // + 32  - bits remaining in interval Low value
    // + number of buffered byte * 8
    writtenBitsAfterQuantizedCoeff = ((OutputBitstreamUnit_t*)EntropyCoderGetBitstreamPtr(picture_control_set_ptr->entropy_coder_ptr))->writtenBitsCount;

    sb_ptr->total_bits = writtenBitsAfterQuantizedCoeff - writtenBitsBeforeQuantizedCoeff;

    picture_control_set_ptr->parent_pcs_ptr->quantized_coeff_num_bits += sb_ptr->quantized_coeffs_bits;

    return;
}

/******************************************************
 * Update Entropy Coding Rows
 *
 * This function is responsible for synchronizing the
 *   processing of Entropy Coding LCU-rows and starts
 *   processing of LCU-rows as soon as their inputs are
 *   available and the previous LCU-row has completed.
 *   At any given time, only one segment row per picture
 *   is being processed.
 *
 * The function has two parts:
 *
 * (1) Update the available row index which tracks
 *   which SB Row-inputs are available.
 *
 * (2) Increment the lcu-row counter as the segment-rows
 *   are completed.
 *
 * Since there is the potentential for thread collusion,
 *   a MUTEX a used to protect the sensitive data and
 *   the execution flow is separated into two paths
 *
 * (A) Initial update.
 *  -Update the Completion Mask [see (1) above]
 *  -If the picture is not currently being processed,
 *     check to see if the next segment-row is available
 *     and start processing.
 * (B) Continued processing
 *  -Upon the completion of a segment-row, check
 *     to see if the next segment-row's inputs have
 *     become available and begin processing if so.
 *
 * On last important point is that the thread-safe
 *   code section is kept minimally short. The MUTEX
 *   should NOT be locked for the entire processing
 *   of the segment-row (B) as this would block other
 *   threads from performing an update (A).
 ******************************************************/
static EbBool UpdateEntropyCodingRows(
    PictureControlSet_t *picture_control_set_ptr,
    uint32_t              *row_index,
    uint32_t               row_count,
    EbBool             *initialProcessCall)
{
    EbBool processNextRow = EB_FALSE;

    // Note, any writes & reads to status variables (e.g. in_progress) in MD-CTRL must be thread-safe
    eb_block_on_mutex(picture_control_set_ptr->entropy_coding_mutex);

    // Update availability mask
    if (*initialProcessCall == EB_TRUE) {
        unsigned i;

        for (i = *row_index; i < *row_index + row_count; ++i) {
            picture_control_set_ptr->entropy_coding_row_array[i] = EB_TRUE;
        }

        while (picture_control_set_ptr->entropy_coding_row_array[picture_control_set_ptr->entropy_coding_current_available_row] == EB_TRUE &&
            picture_control_set_ptr->entropy_coding_current_available_row < picture_control_set_ptr->entropy_coding_row_count)
        {
            ++picture_control_set_ptr->entropy_coding_current_available_row;
        }
    }

    // Release in_progress token
    if (*initialProcessCall == EB_FALSE && picture_control_set_ptr->entropy_coding_in_progress == EB_TRUE) {
        picture_control_set_ptr->entropy_coding_in_progress = EB_FALSE;
    }

    // Test if the picture is not already complete AND not currently being worked on by another ENCDEC process
    if (picture_control_set_ptr->entropy_coding_current_row < picture_control_set_ptr->entropy_coding_row_count &&
        picture_control_set_ptr->entropy_coding_row_array[picture_control_set_ptr->entropy_coding_current_row] == EB_TRUE &&
        picture_control_set_ptr->entropy_coding_in_progress == EB_FALSE)
    {
        // Test if the next LCU-row is ready to go
        if (picture_control_set_ptr->entropy_coding_current_row <= picture_control_set_ptr->entropy_coding_current_available_row)
        {
            picture_control_set_ptr->entropy_coding_in_progress = EB_TRUE;
            *row_index = picture_control_set_ptr->entropy_coding_current_row++;
            processNextRow = EB_TRUE;
        }
    }

    *initialProcessCall = EB_FALSE;

    eb_release_mutex(picture_control_set_ptr->entropy_coding_mutex);

    return processNextRow;
}

/******************************************************
 * Entropy Coding Kernel
 ******************************************************/
void* EntropyCodingKernel(void *input_ptr)
{
    // Context & SCS & PCS
    EntropyCodingContext_t                  *context_ptr = (EntropyCodingContext_t*)input_ptr;
    PictureControlSet_t                     *picture_control_set_ptr;
    SequenceControlSet_t                    *sequence_control_set_ptr;

    // Input
    EbObjectWrapper_t                       *encDecResultsWrapperPtr;
    EncDecResults_t                         *encDecResultsPtr;

    // Output
    EbObjectWrapper_t                       *entropyCodingResultsWrapperPtr;
    EntropyCodingResults_t                  *entropyCodingResultsPtr;

    // SB Loop variables
    LargestCodingUnit_t                     *sb_ptr;
    uint16_t                                   sb_index;
    uint8_t                                    sb_sz;
    uint8_t                                    lcuSizeLog2;
    uint32_t                                   xLcuIndex;
    uint32_t                                   yLcuIndex;
    uint32_t                                   sb_origin_x;
    uint32_t                                   sb_origin_y;
    EbBool                                  lastLcuFlag;
    uint32_t                                   picture_width_in_sb;
    // Variables
    EbBool                                  initialProcessCall;
    for (;;) {

        // Get Mode Decision Results
        eb_get_full_object(
            context_ptr->enc_dec_input_fifo_ptr,
            &encDecResultsWrapperPtr);
        encDecResultsPtr = (EncDecResults_t*)encDecResultsWrapperPtr->object_ptr;
        picture_control_set_ptr = (PictureControlSet_t*)encDecResultsPtr->pictureControlSetWrapperPtr->object_ptr;
        sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
        lastLcuFlag = EB_FALSE;

        // SB Constants

        sb_sz = (uint8_t)sequence_control_set_ptr->sb_size_pix;

        lcuSizeLog2 = (uint8_t)Log2f(sb_sz);
        context_ptr->sb_sz = sb_sz;
        picture_width_in_sb = (sequence_control_set_ptr->luma_width + sb_sz - 1) >> lcuSizeLog2;
#if TILES
        if(picture_control_set_ptr->parent_pcs_ptr->av1_cm->tile_cols * picture_control_set_ptr->parent_pcs_ptr->av1_cm->tile_rows == 1)
#endif
        {
            initialProcessCall = EB_TRUE;
            yLcuIndex = encDecResultsPtr->completedLcuRowIndexStart;

            // LCU-loops
            while (UpdateEntropyCodingRows(picture_control_set_ptr, &yLcuIndex, encDecResultsPtr->completedLcuRowCount, &initialProcessCall) == EB_TRUE)
            {
                uint32_t rowTotalBits = 0;

                if (yLcuIndex == 0) {
                    ResetEntropyCodingPicture(
                        context_ptr,
                        picture_control_set_ptr,
                        sequence_control_set_ptr);
                    picture_control_set_ptr->entropy_coding_pic_done = EB_FALSE;
                }

                for (xLcuIndex = 0; xLcuIndex < picture_width_in_sb; ++xLcuIndex)
                {


                    sb_index = (uint16_t)(xLcuIndex + yLcuIndex * picture_width_in_sb);
                    sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];

                    sb_origin_x = xLcuIndex << lcuSizeLog2;
                    sb_origin_y = yLcuIndex << lcuSizeLog2;
                    context_ptr->sb_origin_x = sb_origin_x;
                    context_ptr->sb_origin_y = sb_origin_y;
                    lastLcuFlag = (sb_index == sequence_control_set_ptr->sb_tot_cnt - 1) ? EB_TRUE : EB_FALSE;

#if TILES 
                    if (sb_index == 0)
                        av1_reset_loop_restoration(picture_control_set_ptr);
#endif
                    // Configure the LCU
                    EntropyCodingConfigureLcu(
                        context_ptr,
                        sb_ptr,
                        picture_control_set_ptr);

                    // Entropy Coding
                    EntropyCodingLcu(
                        context_ptr,
                        sb_ptr,
                        picture_control_set_ptr,
                        sequence_control_set_ptr,
                        sb_origin_x,
                        sb_origin_y,
                        lastLcuFlag,
                        0,
                        0);

                    rowTotalBits += sb_ptr->total_bits;
                }

                // At the end of each LCU-row, send the updated bit-count to Entropy Coding
                {
                    EbObjectWrapper_t *rateControlTaskWrapperPtr;
                    RateControlTasks_t *rateControlTaskPtr;

                    // Get Empty EncDec Results
                    eb_get_empty_object(
                        context_ptr->rate_control_output_fifo_ptr,
                        &rateControlTaskWrapperPtr);
                    rateControlTaskPtr = (RateControlTasks_t*)rateControlTaskWrapperPtr->object_ptr;
                    rateControlTaskPtr->taskType = RC_ENTROPY_CODING_ROW_FEEDBACK_RESULT;
                    rateControlTaskPtr->picture_number = picture_control_set_ptr->picture_number;
                    rateControlTaskPtr->rowNumber = yLcuIndex;
                    rateControlTaskPtr->bitCount = rowTotalBits;

                    rateControlTaskPtr->pictureControlSetWrapperPtr = 0;
                    rateControlTaskPtr->segment_index = ~0u;

                    // Post EncDec Results
                    eb_post_full_object(rateControlTaskWrapperPtr);
                }

                eb_block_on_mutex(picture_control_set_ptr->entropy_coding_mutex);
                if (picture_control_set_ptr->entropy_coding_pic_done == EB_FALSE) {

                    // If the picture is complete, terminate the slice
                    if (picture_control_set_ptr->entropy_coding_current_row == picture_control_set_ptr->entropy_coding_row_count)
                    {
                        uint32_t refIdx;

                        picture_control_set_ptr->entropy_coding_pic_done = EB_TRUE;

                        EncodeSliceFinish(picture_control_set_ptr->entropy_coder_ptr);

                        // Release the List 0 Reference Pictures
                        for (refIdx = 0; refIdx < picture_control_set_ptr->parent_pcs_ptr->ref_list0_count; ++refIdx) {
                            if (picture_control_set_ptr->ref_pic_ptr_array[0] != EB_NULL) {

                                eb_release_object(picture_control_set_ptr->ref_pic_ptr_array[0]);
                            }
                        }

                        // Release the List 1 Reference Pictures
                        for (refIdx = 0; refIdx < picture_control_set_ptr->parent_pcs_ptr->ref_list1_count; ++refIdx) {
                            if (picture_control_set_ptr->ref_pic_ptr_array[1] != EB_NULL) {

                                eb_release_object(picture_control_set_ptr->ref_pic_ptr_array[1]);
                            }
                        }

                        // Get Empty Entropy Coding Results
                        eb_get_empty_object(
                            context_ptr->entropy_coding_output_fifo_ptr,
                            &entropyCodingResultsWrapperPtr);
                        entropyCodingResultsPtr = (EntropyCodingResults_t*)entropyCodingResultsWrapperPtr->object_ptr;
                        entropyCodingResultsPtr->pictureControlSetWrapperPtr = encDecResultsPtr->pictureControlSetWrapperPtr;

                        // Post EntropyCoding Results
                        eb_post_full_object(entropyCodingResultsWrapperPtr);

                    } // End if(PictureCompleteFlag)
                }
                eb_release_mutex(picture_control_set_ptr->entropy_coding_mutex);


            }

        }
#if TILES
        else
        {

             struct PictureParentControlSet_s     *ppcs_ptr = picture_control_set_ptr->parent_pcs_ptr;
             Av1Common *const cm = ppcs_ptr->av1_cm;           
             uint32_t total_size = 0;
             int tile_row, tile_col;
             const int tile_cols = ppcs_ptr->av1_cm->tile_cols;
             const int tile_rows = ppcs_ptr->av1_cm->tile_rows;

             //Entropy Tile Loop
             for (tile_row = 0; tile_row < tile_rows; tile_row++)
             {
                 
                 TileInfo tile_info;
                 av1_tile_set_row(&tile_info, ppcs_ptr, tile_row);

                 for (tile_col = 0; tile_col < tile_cols; tile_col++)
                 { 
                     const int tile_idx = tile_row * tile_cols + tile_col;
                     uint32_t is_last_tile_in_tg = 0;

                     if ( tile_idx == (tile_cols * tile_rows - 1)) {
                         is_last_tile_in_tg = 1;                       
                     }
                     else {
                         is_last_tile_in_tg = 0;                        
                     }

                     reset_ec_tile(
                         total_size,
                         is_last_tile_in_tg,
                         context_ptr,
                         picture_control_set_ptr,
                         sequence_control_set_ptr);

                     av1_tile_set_col(&tile_info, ppcs_ptr, tile_col);
   
                     av1_reset_loop_restoration(picture_control_set_ptr);
                   
                     for (yLcuIndex = cm->tile_row_start_sb[tile_row]; yLcuIndex < (uint32_t)cm->tile_row_start_sb[tile_row + 1]; ++yLcuIndex)
                     {
                         for (xLcuIndex = cm->tile_col_start_sb[tile_col]; xLcuIndex < (uint32_t)cm->tile_col_start_sb[tile_col + 1]; ++xLcuIndex)
                         {
                             
                             int sb_index = (uint16_t)(xLcuIndex + yLcuIndex * picture_width_in_sb);
                             sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];
                             sb_origin_x = xLcuIndex << lcuSizeLog2;
                             sb_origin_y = yLcuIndex << lcuSizeLog2;
                             context_ptr->sb_origin_x = sb_origin_x;
                             context_ptr->sb_origin_y = sb_origin_y;
                             lastLcuFlag = (sb_index == sequence_control_set_ptr->sb_tot_cnt - 1) ? EB_TRUE : EB_FALSE;
                            
                             // Configure the LCU
                             EntropyCodingConfigureLcu(
                                 context_ptr,
                                 sb_ptr,
                                 picture_control_set_ptr);                           

                             // Entropy Coding
                             EntropyCodingLcu(
                                 context_ptr,
                                 sb_ptr,
                                 picture_control_set_ptr,
                                 sequence_control_set_ptr,
                                 sb_origin_x,
                                 sb_origin_y,
                                 lastLcuFlag,
                                 0,
                                 0);
                         }
                     }
                                         
                     EncodeSliceFinish(picture_control_set_ptr->entropy_coder_ptr);
                    
                     int tile_size = picture_control_set_ptr->entropy_coder_ptr->ecWriter.pos;
                     assert(tile_size >= AV1_MIN_TILE_SIZE_BYTES);
                    
                     if (!is_last_tile_in_tg) {
                         
                         OutputBitstreamUnit_t *outputBitstreamPtr = (OutputBitstreamUnit_t*)(picture_control_set_ptr->entropy_coder_ptr->ecOutputBitstreamPtr);
                         uint8_t *buf_data = outputBitstreamPtr->bufferAv1 + total_size;
                         mem_put_le32(buf_data, tile_size - AV1_MIN_TILE_SIZE_BYTES);
                     }                   

                     if (is_last_tile_in_tg==0)                     
                         total_size += 4;                     

                     total_size += tile_size;

                 }

             }

             //the picture is complete, terminate the slice            
             {
                 uint32_t refIdx;         
                 picture_control_set_ptr->entropy_coder_ptr->ec_frame_size = total_size;

                 // Release the List 0 Reference Pictures
                 for (refIdx = 0; refIdx < picture_control_set_ptr->parent_pcs_ptr->ref_list0_count; ++refIdx) {
                     if (picture_control_set_ptr->ref_pic_ptr_array[0] != EB_NULL) {
                         eb_release_object(picture_control_set_ptr->ref_pic_ptr_array[0]);
                     }
                 }

                 // Release the List 1 Reference Pictures
                 for (refIdx = 0; refIdx < picture_control_set_ptr->parent_pcs_ptr->ref_list1_count; ++refIdx) {
                     if (picture_control_set_ptr->ref_pic_ptr_array[1] != EB_NULL) {
                         eb_release_object(picture_control_set_ptr->ref_pic_ptr_array[1]);
                     }
                 }

                 // Get Empty Entropy Coding Results
                 eb_get_empty_object(
                     context_ptr->entropy_coding_output_fifo_ptr,
                     &entropyCodingResultsWrapperPtr);
                 entropyCodingResultsPtr = (EntropyCodingResults_t*)entropyCodingResultsWrapperPtr->object_ptr;
                 entropyCodingResultsPtr->pictureControlSetWrapperPtr = encDecResultsPtr->pictureControlSetWrapperPtr;

                 // Post EntropyCoding Results
                 eb_post_full_object(entropyCodingResultsWrapperPtr);

             } 

        }
#endif
        // Release Mode Decision Results
        eb_release_object(encDecResultsWrapperPtr);

    }

    return EB_NULL;
}
