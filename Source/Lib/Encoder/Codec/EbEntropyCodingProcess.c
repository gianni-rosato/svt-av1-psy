// clang-format off
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
#if TWO_PASS
#include <stdio.h>
#endif
#include "EbEncHandle.h"
#include "EbEntropyCodingProcess.h"
#include "EbEncDecResults.h"
#include "EbEntropyCodingResults.h"
#include "EbRateControlTasks.h"
#include "EbCabacContextModel.h"
#include "EbLog.h"
#define  AV1_MIN_TILE_SIZE_BYTES 1
void eb_av1_reset_loop_restoration(PictureControlSet     *piCSetPtr);


static void rest_context_dctor(EbPtr p)
{
    EbThreadContext* thread_context_ptr = (EbThreadContext*)p;
    EntropyCodingContext *obj = (EntropyCodingContext*)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}


/******************************************************
 * Enc Dec Context Constructor
 ******************************************************/
EbErrorType entropy_coding_context_ctor(
    EbThreadContext     *thread_context_ptr,
    const EbEncHandle   *enc_handle_ptr,
    int                 index,
    int                 rate_control_index)
{
    EntropyCodingContext  *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv = context_ptr;
    thread_context_ptr->dctor = rest_context_dctor;


    context_ptr->is16bit = (EbBool)(enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);;

    // Input/Output System Resource Manager FIFOs
    context_ptr->enc_dec_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->rest_results_resource_ptr, index);
    context_ptr->entropy_coding_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->entropy_coding_results_resource_ptr, index);
    context_ptr->rate_control_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->rate_control_tasks_resource_ptr, rate_control_index);

    return EB_ErrorNone;
}

/***********************************************
 * Entropy Coding Reset Neighbor Arrays
 ***********************************************/
static void EntropyCodingResetNeighborArrays(PictureControlSet *picture_control_set_ptr)
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
    neighbor_array_unit_reset(picture_control_set_ptr->txfm_context_array);
    neighbor_array_unit_reset(picture_control_set_ptr->segmentation_id_pred_array);
    return;
}

void av1_get_syntax_rate_from_cdf(
    int32_t                      *costs,
    const AomCdfProb       *cdf,
    const int32_t                *inv_map);

void eb_av1_cost_tokens_from_cdf(int32_t *costs, const AomCdfProb *cdf,
    const int32_t *inv_map) {
    // int32_t i;
    // AomCdfProb prev_cdf = 0;
    // for (i = 0;; ++i) {
    //     AomCdfProb p15 = AOM_ICDF(cdf[i]) - prev_cdf;
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
    const NmvComponent *const mvcomp,
    MvSubpelPrecision precision) {
    int32_t i, v;
    int32_t sign_cost[2], class_cost[MV_CLASSES], class0_cost[CLASS0_SIZE];
    int32_t bits_cost[MV_OFFSET_BITS][2];
    int32_t class0_fp_cost[CLASS0_SIZE][MV_FP_SIZE], fp_cost[MV_FP_SIZE];
    int32_t class0_hp_cost[2], hp_cost[2];

    eb_av1_cost_tokens_from_cdf(sign_cost, mvcomp->sign_cdf, NULL);
    eb_av1_cost_tokens_from_cdf(class_cost, mvcomp->classes_cdf, NULL);
    eb_av1_cost_tokens_from_cdf(class0_cost, mvcomp->class0_cdf, NULL);
    for (i = 0; i < MV_OFFSET_BITS; ++i)
        eb_av1_cost_tokens_from_cdf(bits_cost[i], mvcomp->bits_cdf[i], NULL);
    for (i = 0; i < CLASS0_SIZE; ++i)
        eb_av1_cost_tokens_from_cdf(class0_fp_cost[i], mvcomp->class0_fp_cdf[i], NULL);
    eb_av1_cost_tokens_from_cdf(fp_cost, mvcomp->fp_cdf, NULL);

    if (precision > MV_SUBPEL_LOW_PRECISION) {
        eb_av1_cost_tokens_from_cdf(class0_hp_cost, mvcomp->class0_hp_cdf, NULL);
        eb_av1_cost_tokens_from_cdf(hp_cost, mvcomp->hp_cdf, NULL);
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
        if (c == MV_CLASS_0)
            cost += class0_cost[d];
        else {
            const int32_t b = c + CLASS0_BITS - 1; /* number of bits */
            for (i = 0; i < b; ++i) cost += bits_cost[i][((d >> i) & 1)];
        }
        if (precision > MV_SUBPEL_NONE) {
            if (c == MV_CLASS_0)
                cost += class0_fp_cost[d][f];
            else
                cost += fp_cost[f];
            if (precision > MV_SUBPEL_LOW_PRECISION) {
                if (c == MV_CLASS_0)
                    cost += class0_hp_cost[e];
                else
                    cost += hp_cost[e];
            }
        }
        mvcost[v] = cost + sign_cost[0];
        mvcost[-v] = cost + sign_cost[1];
    }
}
void eb_av1_build_nmv_cost_table(int32_t *mvjoint, int32_t *mvcost[2],
    const NmvContext *ctx,
    MvSubpelPrecision precision) {
    eb_av1_cost_tokens_from_cdf(mvjoint, ctx->joints_cdf, NULL);
    build_nmv_component_cost_table(mvcost[0], &ctx->comps[0], precision);
    build_nmv_component_cost_table(mvcost[1], &ctx->comps[1], precision);
}

/**************************************************
 * Reset Entropy Coding Picture
 **************************************************/
static void reset_entropy_coding_picture(
    EntropyCodingContext  *context_ptr,
    PictureControlSet     *picture_control_set_ptr,
    SequenceControlSet    *sequence_control_set_ptr)
{
    reset_bitstream(entropy_coder_get_bitstream_ptr(picture_control_set_ptr->entropy_coder_ptr));

    uint32_t                       entropyCodingQp;

    context_ptr->is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    // QP
#if ADD_DELTA_QP_SUPPORT
#else
    context_ptr->qp = picture_control_set_ptr->picture_qp;
#endif
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation
    entropyCodingQp = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;

#if ADD_DELTA_QP_SUPPORT
    picture_control_set_ptr->parent_pcs_ptr->prev_qindex = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc)
        assert(picture_control_set_ptr->parent_pcs_ptr->frm_hdr.delta_lf_params.delta_lf_present == 0);
    if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.delta_lf_params.delta_lf_present) {
        picture_control_set_ptr->parent_pcs_ptr->prev_delta_lf_from_base = 0;
        const int32_t frame_lf_count =
            picture_control_set_ptr->parent_pcs_ptr->monochrome == 0 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
        for (int32_t lf_id = 0; lf_id < frame_lf_count; ++lf_id)
            picture_control_set_ptr->parent_pcs_ptr->prev_delta_lf[lf_id] = 0;
    }
#endif

    // pass the ent
    OutputBitstreamUnit *output_bitstream_ptr = (OutputBitstreamUnit*)(picture_control_set_ptr->entropy_coder_ptr->ec_output_bitstream_ptr);
    //****************************************************************//

    uint8_t *data = output_bitstream_ptr->buffer_av1;
    picture_control_set_ptr->entropy_coder_ptr->ec_writer.allow_update_cdf = !picture_control_set_ptr->parent_pcs_ptr->large_scale_tile;
    picture_control_set_ptr->entropy_coder_ptr->ec_writer.allow_update_cdf =
        picture_control_set_ptr->entropy_coder_ptr->ec_writer.allow_update_cdf && !frm_hdr->disable_cdf_update;
    aom_start_encode(&picture_control_set_ptr->entropy_coder_ptr->ec_writer, data);

    // ADD Reset here
    if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame != PRIMARY_REF_NONE)
        memcpy(picture_control_set_ptr->entropy_coder_ptr->fc, &picture_control_set_ptr->ref_frame_context[picture_control_set_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame], sizeof(FRAME_CONTEXT));
    else
        reset_entropy_coder(
            sequence_control_set_ptr->encode_context_ptr,
            picture_control_set_ptr->entropy_coder_ptr,
            entropyCodingQp,
            picture_control_set_ptr->slice_type);
    EntropyCodingResetNeighborArrays(picture_control_set_ptr);

    return;
}

static void reset_ec_tile(
    uint32_t  total_size,
    uint32_t  is_last_tile_in_tg,
    EntropyCodingContext  *context_ptr,
    PictureControlSet     *picture_control_set_ptr,
    SequenceControlSet    *sequence_control_set_ptr)
{
    reset_bitstream(entropy_coder_get_bitstream_ptr(picture_control_set_ptr->entropy_coder_ptr));

    uint32_t                       entropy_coding_qp;

    context_ptr->is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    // QP
#if ADD_DELTA_QP_SUPPORT
#else
    context_ptr->qp = picture_control_set_ptr->picture_qp;
#endif
    entropy_coding_qp = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
#if ADD_DELTA_QP_SUPPORT
    picture_control_set_ptr->parent_pcs_ptr->prev_qindex = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc)
        assert(picture_control_set_ptr->parent_pcs_ptr->frm_hdr.delta_lf_params.delta_lf_present == 0);
    if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.delta_lf_params.delta_lf_present) {
        picture_control_set_ptr->parent_pcs_ptr->prev_delta_lf_from_base = 0;
        const int32_t frame_lf_count =
            picture_control_set_ptr->parent_pcs_ptr->monochrome == 0 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
        for (int32_t lf_id = 0; lf_id < frame_lf_count; ++lf_id)
            picture_control_set_ptr->parent_pcs_ptr->prev_delta_lf[lf_id] = 0;
    }
#endif

    // pass the ent
    OutputBitstreamUnit *output_bitstream_ptr = (OutputBitstreamUnit*)(picture_control_set_ptr->entropy_coder_ptr->ec_output_bitstream_ptr);
    //****************************************************************//

    uint8_t *data = output_bitstream_ptr->buffer_av1 + total_size;
    picture_control_set_ptr->entropy_coder_ptr->ec_writer.allow_update_cdf = !picture_control_set_ptr->parent_pcs_ptr->large_scale_tile;
    picture_control_set_ptr->entropy_coder_ptr->ec_writer.allow_update_cdf =
        picture_control_set_ptr->entropy_coder_ptr->ec_writer.allow_update_cdf && !frm_hdr->disable_cdf_update;

    //if not last tile, advance buffer by 4B to leave space for tile Size
    if (is_last_tile_in_tg == 0)
        data += 4;

    aom_start_encode(&picture_control_set_ptr->entropy_coder_ptr->ec_writer, data);
    if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame != PRIMARY_REF_NONE)
        memcpy(picture_control_set_ptr->entropy_coder_ptr->fc, &picture_control_set_ptr->ref_frame_context[picture_control_set_ptr->parent_pcs_ptr->frm_hdr.primary_ref_frame], sizeof(FRAME_CONTEXT));
    else
        //reset probabilities
        reset_entropy_coder(
            sequence_control_set_ptr->encode_context_ptr,
            picture_control_set_ptr->entropy_coder_ptr,
            entropy_coding_qp,
            picture_control_set_ptr->slice_type);
    EntropyCodingResetNeighborArrays(picture_control_set_ptr);

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
    PictureControlSet *picture_control_set_ptr,
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

        for (i = *row_index; i < *row_index + row_count; ++i)
            picture_control_set_ptr->entropy_coding_row_array[i] = EB_TRUE;
        while (picture_control_set_ptr->entropy_coding_row_array[picture_control_set_ptr->entropy_coding_current_available_row] == EB_TRUE &&
            picture_control_set_ptr->entropy_coding_current_available_row < picture_control_set_ptr->entropy_coding_row_count)
        {
            ++picture_control_set_ptr->entropy_coding_current_available_row;
        }
    }

    // Release in_progress token
    if (*initialProcessCall == EB_FALSE && picture_control_set_ptr->entropy_coding_in_progress == EB_TRUE)
        picture_control_set_ptr->entropy_coding_in_progress = EB_FALSE;
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
#if TWO_PASS
/******************************************************
 * Write Stat to File
 * write stat_struct per frame in the first pass
 ******************************************************/
void write_stat_to_file(
    SequenceControlSet    *sequence_control_set_ptr,
    stat_struct_t          stat_struct,
    uint64_t               ref_poc)
{
    eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->stat_file_mutex);
    int32_t fseek_return_value = fseek(sequence_control_set_ptr->static_config.output_stat_file, (long)ref_poc * sizeof(stat_struct_t), SEEK_SET);
    if (fseek_return_value != 0)
        SVT_LOG("Error in fseek  returnVal %i\n", fseek_return_value);
    fwrite(&stat_struct,
        sizeof(stat_struct_t),
        (size_t)1,
        sequence_control_set_ptr->static_config.output_stat_file);
    eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->stat_file_mutex);
}
#endif
/******************************************************
 * Entropy Coding Kernel
 ******************************************************/
void* entropy_coding_kernel(void *input_ptr)
{
    // Context & SCS & PCS
    EbThreadContext                       *thread_context_ptr = (EbThreadContext*)input_ptr;
    EntropyCodingContext                  *context_ptr = (EntropyCodingContext*)thread_context_ptr->priv;
    PictureControlSet                     *picture_control_set_ptr;
    SequenceControlSet                    *sequence_control_set_ptr;

    // Input
    EbObjectWrapper                       *encDecResultsWrapperPtr;
    EncDecResults                         *encDecResultsPtr;

    // Output
    EbObjectWrapper                       *entropyCodingResultsWrapperPtr;
    EntropyCodingResults                  *entropyCodingResultsPtr;

    // SB Loop variables
    SuperBlock                                 *sb_ptr;
    uint16_t                                   sb_index;
    uint8_t                                    sb_sz;
    uint8_t                                    lcuSizeLog2;
    uint32_t                                   x_lcu_index;
    uint32_t                                   y_lcu_index;
    uint32_t                                   sb_origin_x;
    uint32_t                                   sb_origin_y;
    uint32_t                                   picture_width_in_sb;
    // Variables
    EbBool                                  initialProcessCall;
    for (;;) {
        // Get Mode Decision Results
        eb_get_full_object(
            context_ptr->enc_dec_input_fifo_ptr,
            &encDecResultsWrapperPtr);
        encDecResultsPtr = (EncDecResults*)encDecResultsWrapperPtr->object_ptr;
        picture_control_set_ptr = (PictureControlSet*)encDecResultsPtr->picture_control_set_wrapper_ptr->object_ptr;
        sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
        // SB Constants

        sb_sz = (uint8_t)sequence_control_set_ptr->sb_size_pix;

        lcuSizeLog2 = (uint8_t)Log2f(sb_sz);
        context_ptr->sb_sz = sb_sz;
        picture_width_in_sb = (sequence_control_set_ptr->seq_header.max_frame_width + sb_sz - 1) >> lcuSizeLog2;
        if(picture_control_set_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_cols * picture_control_set_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_rows == 1)

        {
            initialProcessCall = EB_TRUE;
            y_lcu_index = encDecResultsPtr->completed_lcu_row_index_start;

            // LCU-loops
            while (UpdateEntropyCodingRows(picture_control_set_ptr, &y_lcu_index, encDecResultsPtr->completed_lcu_row_count, &initialProcessCall) == EB_TRUE)
            {
                uint32_t rowTotalBits = 0;

                if (y_lcu_index == 0) {
                    reset_entropy_coding_picture(
                        context_ptr,
                        picture_control_set_ptr,
                        sequence_control_set_ptr);
                    picture_control_set_ptr->entropy_coding_pic_done = EB_FALSE;
                }

                for (x_lcu_index = 0; x_lcu_index < picture_width_in_sb; ++x_lcu_index)
                {
                    sb_index = (uint16_t)(x_lcu_index + y_lcu_index * picture_width_in_sb);
                    sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];

                    sb_origin_x = x_lcu_index << lcuSizeLog2;
                    sb_origin_y = y_lcu_index << lcuSizeLog2;
                    context_ptr->sb_origin_x = sb_origin_x;
                    context_ptr->sb_origin_y = sb_origin_y;
                    if (sb_index == 0)
                        eb_av1_reset_loop_restoration(picture_control_set_ptr);
#if PAL_SUP
                    if (sb_index == 0)
                        context_ptr->tok = picture_control_set_ptr->tile_tok[0][0];
#endif
                    sb_ptr->total_bits = 0;
                    uint32_t prev_pos = sb_index ? picture_control_set_ptr->entropy_coder_ptr->ec_writer.ec.offs : 0;//residual_bc.pos
                    EbPictureBufferDesc *coeff_picture_ptr = sb_ptr->quantized_coeff;
                    write_sb(
                        context_ptr,
                        sb_ptr,
                        picture_control_set_ptr,
                        picture_control_set_ptr->entropy_coder_ptr,
                        coeff_picture_ptr);
                    sb_ptr->total_bits = (picture_control_set_ptr->entropy_coder_ptr->ec_writer.ec.offs - prev_pos) << 3;
                    picture_control_set_ptr->parent_pcs_ptr->quantized_coeff_num_bits += sb_ptr->total_bits;
                    rowTotalBits += sb_ptr->total_bits;
                }

                // At the end of each LCU-row, send the updated bit-count to Entropy Coding
                {
                    EbObjectWrapper *rateControlTaskWrapperPtr;
                    RateControlTasks *rateControlTaskPtr;

                    // Get Empty EncDec Results
                    eb_get_empty_object(
                        context_ptr->rate_control_output_fifo_ptr,
                        &rateControlTaskWrapperPtr);
                    rateControlTaskPtr = (RateControlTasks*)rateControlTaskWrapperPtr->object_ptr;
                    rateControlTaskPtr->task_type = RC_ENTROPY_CODING_ROW_FEEDBACK_RESULT;
                    rateControlTaskPtr->picture_number = picture_control_set_ptr->picture_number;
                    rateControlTaskPtr->row_number = y_lcu_index;
                    rateControlTaskPtr->bit_count = rowTotalBits;

                    rateControlTaskPtr->picture_control_set_wrapper_ptr = 0;
                    rateControlTaskPtr->segment_index = ~0u;

                    // Post EncDec Results
                    eb_post_full_object(rateControlTaskWrapperPtr);
                }

                eb_block_on_mutex(picture_control_set_ptr->entropy_coding_mutex);
                if (picture_control_set_ptr->entropy_coding_pic_done == EB_FALSE) {
                    // If the picture is complete, terminate the slice
                    if (picture_control_set_ptr->entropy_coding_current_row == picture_control_set_ptr->entropy_coding_row_count)
                    {
                        uint32_t ref_idx;

                        picture_control_set_ptr->entropy_coding_pic_done = EB_TRUE;

                        encode_slice_finish(picture_control_set_ptr->entropy_coder_ptr);
#if TWO_PASS
                        // for Non Reference frames
                        if (sequence_control_set_ptr->use_output_stat_file &&
                            !picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                            write_stat_to_file(
                                sequence_control_set_ptr,
                                *picture_control_set_ptr->parent_pcs_ptr->stat_struct_first_pass_ptr,
                                picture_control_set_ptr->parent_pcs_ptr->picture_number);
#endif
                        // Release the List 0 Reference Pictures
                        for (ref_idx = 0; ref_idx < picture_control_set_ptr->parent_pcs_ptr->ref_list0_count; ++ref_idx) {
#if TWO_PASS
                            if (sequence_control_set_ptr->use_output_stat_file &&
                                picture_control_set_ptr->ref_pic_ptr_array[0][ref_idx] != EB_NULL && picture_control_set_ptr->ref_pic_ptr_array[0][ref_idx]->live_count == 1)
                                write_stat_to_file(
                                    sequence_control_set_ptr,
                                    ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[0][ref_idx]->object_ptr)->stat_struct,
                                    ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[0][ref_idx]->object_ptr)->ref_poc);
#endif
                            if (picture_control_set_ptr->ref_pic_ptr_array[0][ref_idx] != EB_NULL) {

                                eb_release_object(picture_control_set_ptr->ref_pic_ptr_array[0][ref_idx]);
                            }
                        }

                        // Release the List 1 Reference Pictures
                        for (ref_idx = 0; ref_idx < picture_control_set_ptr->parent_pcs_ptr->ref_list1_count; ++ref_idx) {
#if TWO_PASS
                            if (sequence_control_set_ptr->use_output_stat_file &&
                                picture_control_set_ptr->ref_pic_ptr_array[1][ref_idx] != EB_NULL && picture_control_set_ptr->ref_pic_ptr_array[1][ref_idx]->live_count == 1)
                                write_stat_to_file(
                                    sequence_control_set_ptr,
                                    ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[1][ref_idx]->object_ptr)->stat_struct,
                                    ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[1][ref_idx]->object_ptr)->ref_poc);
#endif

                            if (picture_control_set_ptr->ref_pic_ptr_array[1][ref_idx] != EB_NULL)
                                eb_release_object(picture_control_set_ptr->ref_pic_ptr_array[1][ref_idx]);
                        }

                        // Get Empty Entropy Coding Results
                        eb_get_empty_object(
                            context_ptr->entropy_coding_output_fifo_ptr,
                            &entropyCodingResultsWrapperPtr);
                        entropyCodingResultsPtr = (EntropyCodingResults*)entropyCodingResultsWrapperPtr->object_ptr;
                        entropyCodingResultsPtr->picture_control_set_wrapper_ptr = encDecResultsPtr->picture_control_set_wrapper_ptr;

                        // Post EntropyCoding Results
                        eb_post_full_object(entropyCodingResultsWrapperPtr);
                    } // End if(PictureCompleteFlag)
                }
                eb_release_mutex(picture_control_set_ptr->entropy_coding_mutex);
            }
        }
        else
        {
             struct PictureParentControlSet     *ppcs_ptr = picture_control_set_ptr->parent_pcs_ptr;
             Av1Common *const cm = ppcs_ptr->av1_cm;
             uint32_t total_size = 0;
             int tile_row, tile_col;
             const int tile_cols = ppcs_ptr->av1_cm->tiles_info.tile_cols;
             const int tile_rows = ppcs_ptr->av1_cm->tiles_info.tile_rows;

             //Entropy Tile Loop
             for (tile_row = 0; tile_row < tile_rows; tile_row++)
             {
                 TileInfo tile_info;
                 eb_av1_tile_set_row(&tile_info, &cm->tiles_info, cm->mi_rows, tile_row);

                 for (tile_col = 0; tile_col < tile_cols; tile_col++)
                 {
                     const int tile_idx = tile_row * tile_cols + tile_col;
                     uint32_t is_last_tile_in_tg = 0;

                     if ( tile_idx == (tile_cols * tile_rows - 1))
                         is_last_tile_in_tg = 1;
                     else
                         is_last_tile_in_tg = 0;
                     reset_ec_tile(
                         total_size,
                         is_last_tile_in_tg,
                         context_ptr,
                         picture_control_set_ptr,
                         sequence_control_set_ptr);
#if PAL_SUP
                     context_ptr->tok = picture_control_set_ptr->tile_tok[0][0];
#endif

                     eb_av1_tile_set_col(&tile_info, &cm->tiles_info, cm->mi_cols, tile_col);

                     eb_av1_reset_loop_restoration(picture_control_set_ptr);
                     int sb_size_log2 = sequence_control_set_ptr->seq_header.sb_size_log2;

                     for ((y_lcu_index = cm->tiles_info.tile_row_start_mi[tile_row] >> sb_size_log2);
                          (y_lcu_index < (uint32_t)cm->tiles_info.tile_row_start_mi[tile_row + 1] >> sb_size_log2);
                          y_lcu_index++)
                     {
                         for (x_lcu_index = (cm->tiles_info.tile_col_start_mi[tile_col] >> sb_size_log2);
                              x_lcu_index < ((uint32_t)cm->tiles_info.tile_col_start_mi[tile_col + 1] >> sb_size_log2);
                              x_lcu_index++)
                         {
                             int sb_index = (uint16_t)(x_lcu_index + y_lcu_index * picture_width_in_sb);
                             sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];
                             sb_origin_x = x_lcu_index << lcuSizeLog2;
                             sb_origin_y = y_lcu_index << lcuSizeLog2;
                             context_ptr->sb_origin_x = sb_origin_x;
                             context_ptr->sb_origin_y = sb_origin_y;
                             sb_ptr->total_bits = 0;
                             uint32_t prev_pos = sb_index ? picture_control_set_ptr->entropy_coder_ptr->ec_writer.ec.offs : 0;//residual_bc.pos
                             EbPictureBufferDesc *coeff_picture_ptr = sb_ptr->quantized_coeff;
                             write_sb(
                                 context_ptr,
                                 sb_ptr,
                                 picture_control_set_ptr,
                                 picture_control_set_ptr->entropy_coder_ptr,
                                 coeff_picture_ptr);
                             sb_ptr->total_bits = (picture_control_set_ptr->entropy_coder_ptr->ec_writer.ec.offs - prev_pos) << 3;
                             picture_control_set_ptr->parent_pcs_ptr->quantized_coeff_num_bits += sb_ptr->total_bits;
                         }
                     }

                     encode_slice_finish(picture_control_set_ptr->entropy_coder_ptr);

                     int tile_size = picture_control_set_ptr->entropy_coder_ptr->ec_writer.pos;
                     assert(tile_size >= AV1_MIN_TILE_SIZE_BYTES);

                     if (!is_last_tile_in_tg) {
                         OutputBitstreamUnit *output_bitstream_ptr = (OutputBitstreamUnit*)(picture_control_set_ptr->entropy_coder_ptr->ec_output_bitstream_ptr);
                         uint8_t *buf_data = output_bitstream_ptr->buffer_av1 + total_size;
                         mem_put_le32(buf_data, tile_size - AV1_MIN_TILE_SIZE_BYTES);
                     }

                     if (is_last_tile_in_tg==0)
                         total_size += 4;

                     total_size += tile_size;
                 }
             }

             //the picture is complete, terminate the slice
             {
                 uint32_t ref_idx;
                 picture_control_set_ptr->entropy_coder_ptr->ec_frame_size = total_size;

                 // Release the List 0 Reference Pictures
                 for (ref_idx = 0; ref_idx < picture_control_set_ptr->parent_pcs_ptr->ref_list0_count; ++ref_idx) {
                     if (picture_control_set_ptr->ref_pic_ptr_array[0][ref_idx] != EB_NULL)
                         eb_release_object(picture_control_set_ptr->ref_pic_ptr_array[0][ref_idx]);
                 }

                 // Release the List 1 Reference Pictures
                 for (ref_idx = 0; ref_idx < picture_control_set_ptr->parent_pcs_ptr->ref_list1_count; ++ref_idx) {
                     if (picture_control_set_ptr->ref_pic_ptr_array[1][ref_idx] != EB_NULL)
                         eb_release_object(picture_control_set_ptr->ref_pic_ptr_array[1][ref_idx]);
                 }

                 // Get Empty Entropy Coding Results
                 eb_get_empty_object(
                     context_ptr->entropy_coding_output_fifo_ptr,
                     &entropyCodingResultsWrapperPtr);
                 entropyCodingResultsPtr = (EntropyCodingResults*)entropyCodingResultsWrapperPtr->object_ptr;
                 entropyCodingResultsPtr->picture_control_set_wrapper_ptr = encDecResultsPtr->picture_control_set_wrapper_ptr;

                 // Post EntropyCoding Results
                 eb_post_full_object(entropyCodingResultsWrapperPtr);
             }
        }

        // Release Mode Decision Results
        eb_release_object(encDecResultsWrapperPtr);
    }

    return EB_NULL;
}
// clang-format on
