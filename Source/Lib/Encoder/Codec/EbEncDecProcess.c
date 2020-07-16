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

#include "EbEncHandle.h"
#include "EbEncDecTasks.h"
#include "EbEncDecResults.h"
#include "EbCodingLoop.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbUtility.h"
#include "grainSynthesis.h"
//To fix warning C4013: 'convert_16bit_to_8bit' undefined; assuming extern returning int
#include "common_dsp_rtcd.h"
#include "EbRateDistortionCost.h"
#include "EbPictureDecisionProcess.h"

#if !REMOVE_MR_MACRO
#if MR_MODE
#define MR_MODE_MULTI_PASS_PD 1
#if IMPROVE_SUB_PEL
#define MR_MODE_SUB_PEL 1
#endif
#else
#define MR_MODE_MULTI_PASS_PD 0
#if IMPROVE_SUB_PEL
#define MR_MODE_SUB_PEL 0
#endif
#endif
#endif
#define FC_SKIP_TX_SR_TH025 125 // Fast cost skip tx search threshold.
#define FC_SKIP_TX_SR_TH010 110 // Fast cost skip tx search threshold.
void eb_av1_cdef_search(EncDecContext *context_ptr, SequenceControlSet *scs_ptr,
                        PictureControlSet *pcs_ptr);

void av1_cdef_frame16bit(uint8_t is_16bit, SequenceControlSet *scs_ptr, PictureControlSet *pCs);

void eb_av1_add_film_grain(EbPictureBufferDesc *src, EbPictureBufferDesc *dst,
                           AomFilmGrain *film_grain_ptr);

void eb_av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm,
                                                 int32_t after_cdef);
void eb_av1_loop_restoration_filter_frame(Yv12BufferConfig *frame, Av1Common *cm,
                                          int32_t optimized_lr);

static void enc_dec_context_dctor(EbPtr p) {
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    EncDecContext *  obj                = (EncDecContext *)thread_context_ptr->priv;
    EB_DELETE(obj->md_context);
    EB_DELETE(obj->residual_buffer);
    EB_DELETE(obj->transform_buffer);
    EB_DELETE(obj->inverse_quant_buffer);
    EB_DELETE(obj->input_sample16bit_buffer);
    if (obj->is_md_rate_estimation_ptr_owner) EB_FREE(obj->md_rate_estimation_ptr);
    EB_FREE_ARRAY(obj);
}

/******************************************************
 * Enc Dec Context Constructor
 ******************************************************/
EbErrorType enc_dec_context_ctor(EbThreadContext *  thread_context_ptr,
                                 const EbEncHandle *enc_handle_ptr, int index, int tasks_index,
                                 int demux_index)

{
    const EbSvtAv1EncConfiguration *static_config =
        &enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config;
    EbBool        is_16bit                 = (EbBool)(static_config->encoder_bit_depth > EB_8BIT);
    EbColorFormat color_format             = static_config->encoder_color_format;
    uint8_t       enable_hbd_mode_decision = static_config->enable_hbd_mode_decision;

    EncDecContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = enc_dec_context_dctor;

    context_ptr->is_16bit     = is_16bit;
    context_ptr->color_format = color_format;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->enc_dec_tasks_resource_ptr, index);
    context_ptr->enc_dec_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->enc_dec_results_resource_ptr, index);
    context_ptr->enc_dec_feedback_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_tasks_resource_ptr, tasks_index);
    context_ptr->picture_demux_output_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);

    // MD rate Estimation tables
    EB_MALLOC(context_ptr->md_rate_estimation_ptr, sizeof(MdRateEstimationContext));
    context_ptr->is_md_rate_estimation_ptr_owner = EB_TRUE;

    // Prediction Buffer
    {
        EbPictureBufferDescInitData init_data;

        init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        init_data.max_width          = SB_STRIDE_Y;
        init_data.max_height         = SB_STRIDE_Y;
        init_data.bit_depth          = EB_8BIT;
        init_data.left_padding       = 0;
        init_data.right_padding      = 0;
        init_data.top_padding        = 0;
        init_data.bot_padding        = 0;
        init_data.split_mode         = EB_FALSE;
        init_data.color_format       = color_format;

        context_ptr->input_sample16bit_buffer = (EbPictureBufferDesc *)NULL;
        if (is_16bit || static_config->is_16bit_pipeline) {
            init_data.bit_depth = EB_16BIT;

            EB_NEW(context_ptr->input_sample16bit_buffer,
                   eb_picture_buffer_desc_ctor,
                   (EbPtr)&init_data);
            init_data.bit_depth = static_config->is_16bit_pipeline ? static_config->encoder_bit_depth : init_data.bit_depth;
        }
    }

    // Scratch Coeff Buffer
    {
        EbPictureBufferDescInitData init_data;

        init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        init_data.max_width          = SB_STRIDE_Y;
        init_data.max_height         = SB_STRIDE_Y;
        init_data.bit_depth          = EB_16BIT;
        init_data.color_format       = color_format;
        init_data.left_padding       = 0;
        init_data.right_padding      = 0;
        init_data.top_padding        = 0;
        init_data.bot_padding        = 0;
        init_data.split_mode         = EB_FALSE;

        EbPictureBufferDescInitData init_32bit_data;

        init_32bit_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        init_32bit_data.max_width          = SB_STRIDE_Y;
        init_32bit_data.max_height         = SB_STRIDE_Y;
        init_32bit_data.bit_depth          = EB_32BIT;
        init_32bit_data.color_format       = color_format;
        init_32bit_data.left_padding       = 0;
        init_32bit_data.right_padding      = 0;
        init_32bit_data.top_padding        = 0;
        init_32bit_data.bot_padding        = 0;
        init_32bit_data.split_mode         = EB_FALSE;
        EB_NEW(context_ptr->inverse_quant_buffer,
               eb_picture_buffer_desc_ctor,
               (EbPtr)&init_32bit_data);
        EB_NEW(context_ptr->transform_buffer, eb_picture_buffer_desc_ctor, (EbPtr)&init_32bit_data);
        EB_NEW(context_ptr->residual_buffer, eb_picture_buffer_desc_ctor, (EbPtr)&init_data);
    }

    // Mode Decision Context
#if SB64_MEM_OPT
#if CHANGE_HBD_MODE
    EB_NEW(context_ptr->md_context,
           mode_decision_context_ctor,
           color_format,
           static_config->super_block_size,
           0,
           0,
           enable_hbd_mode_decision == DEFAULT ? 2 : enable_hbd_mode_decision ,
           static_config->screen_content_mode);
#else
    EB_NEW(context_ptr->md_context,
           mode_decision_context_ctor,
           color_format,
           static_config->super_block_size,
           0,
           0,
           enable_hbd_mode_decision,
           static_config->screen_content_mode);
#endif
#else
    EB_NEW(context_ptr->md_context,
           mode_decision_context_ctor,
           color_format,
           0,
           0,
           enable_hbd_mode_decision,
           static_config->screen_content_mode);
#endif
    if (enable_hbd_mode_decision)
        context_ptr->md_context->input_sample16bit_buffer = context_ptr->input_sample16bit_buffer;

    context_ptr->md_context->enc_dec_context_ptr = context_ptr;

    return EB_ErrorNone;
}

/**************************************************
 * Reset Segmentation Map
 *************************************************/
static void reset_segmentation_map(SegmentationNeighborMap *segmentation_map) {
    if (segmentation_map->data != NULL)
        EB_MEMSET(segmentation_map->data, ~0, segmentation_map->map_size);
}

/**************************************************
 * Reset Mode Decision Neighbor Arrays
 *************************************************/
static void reset_encode_pass_neighbor_arrays(PictureControlSet *pcs_ptr, uint16_t tile_idx) {
    neighbor_array_unit_reset(pcs_ptr->ep_intra_luma_mode_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_intra_chroma_mode_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_mv_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_skip_flag_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_mode_type_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_leaf_depth_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_luma_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_luma_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_partition_context_neighbor_array[tile_idx]);
    // TODO(Joel): 8-bit ep_luma_recon_neighbor_array (Cb,Cr) when is_16bit==0?
    EbBool is_16bit =
        (EbBool)(pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    if (is_16bit || pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline) {
        neighbor_array_unit_reset(pcs_ptr->ep_luma_recon_neighbor_array16bit[tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->ep_cb_recon_neighbor_array16bit[tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->ep_cr_recon_neighbor_array16bit[tile_idx]);
    }
    return;
}

/**************************************************
 * Reset Coding Loop
 **************************************************/
static void reset_enc_dec(EncDecContext *context_ptr, PictureControlSet *pcs_ptr,
                          SequenceControlSet *scs_ptr, uint32_t segment_index) {
    context_ptr->is_16bit = (EbBool)(scs_ptr->static_config.encoder_bit_depth > EB_8BIT) || (EbBool)(scs_ptr->static_config.is_16bit_pipeline);
    context_ptr->bit_depth = scs_ptr->static_config.encoder_bit_depth;
    uint16_t picture_qp   = pcs_ptr->picture_qp;
    uint16_t tile_group_idx = context_ptr->tile_group_index;
#if !QP2QINDEX
    context_ptr->qp = picture_qp;
    context_ptr->qp_index =
        pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present
            ? (uint8_t)quantizer_to_qindex[context_ptr->qp]
            : (uint8_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation

    context_ptr->chroma_qp = (uint8_t)context_ptr->qp;
    // Lambda Assignement
    context_ptr->qp_index =
        (uint8_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
#endif
#if TPL_LA_LAMBDA_SCALING
    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
#if QP2QINDEX
        &context_ptr->pic_fast_lambda[EB_8_BIT_MD],
        &context_ptr->pic_full_lambda[EB_8_BIT_MD],
#else
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
#endif
        8,
#if QP2QINDEX
        pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx,
#else
        context_ptr->qp_index,
#endif
        EB_TRUE);

    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
#if QP2QINDEX
        &context_ptr->pic_fast_lambda[EB_10_BIT_MD],
        &context_ptr->pic_full_lambda[EB_10_BIT_MD],
#else
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
#endif
        10,
#if QP2QINDEX
        pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx,
#else
        context_ptr->qp_index,
#endif
        EB_TRUE);
#else
    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
#if QP2QINDEX
        &context_ptr->pic_fast_lambda,
        &context_ptr->pic_full_lambda,
#else
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
#endif
        (uint8_t)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
#if QP2QINDEX
        pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx,
#else
        context_ptr->qp_index,
#endif
        EB_TRUE);
#endif
    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    if (context_ptr->is_md_rate_estimation_ptr_owner) {
        EB_FREE(context_ptr->md_rate_estimation_ptr);
        context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
    }
    context_ptr->md_rate_estimation_ptr = pcs_ptr->md_rate_estimation_array;
    if (segment_index == 0) {
        if (context_ptr->tile_group_index == 0) {
            reset_segmentation_map(pcs_ptr->segmentation_neighbor_map);
        }

        for (uint16_t r =
                 pcs_ptr->parent_pcs_ptr->tile_group_info[tile_group_idx].tile_group_tile_start_y;
             r < pcs_ptr->parent_pcs_ptr->tile_group_info[tile_group_idx].tile_group_tile_end_y;
             r++) {
            for (uint16_t c = pcs_ptr->parent_pcs_ptr->tile_group_info[tile_group_idx]
                                  .tile_group_tile_start_x;
                 c < pcs_ptr->parent_pcs_ptr->tile_group_info[tile_group_idx].tile_group_tile_end_x;
                 c++) {
                uint16_t tile_idx = c + r * pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_cols;
                reset_encode_pass_neighbor_arrays(pcs_ptr, tile_idx);
            }
        }
    }

    return;
}

#if !QP2QINDEX
/******************************************************
 * EncDec Configure SB
 ******************************************************/
static void enc_dec_configure_sb(EncDecContext *context_ptr, SuperBlock *sb_ptr,
                                 PictureControlSet *pcs_ptr, uint8_t sb_qp) {
    context_ptr->qp = sb_qp;

    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation
    context_ptr->chroma_qp = (uint8_t)context_ptr->qp;
    /* Note(CHKN) : when Qp modulation varies QP on a sub-SB(CU) basis,  Lamda has to change based on Cu->QP , and then this code has to move inside the CU loop in MD */
    (void)sb_ptr;
    context_ptr->qp_index =
        (uint8_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
        (uint8_t)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index,
        EB_TRUE);

    return;
}
#endif

/******************************************************
 * Update MD Segments
 *
 * This function is responsible for synchronizing the
 *   processing of MD Segment-rows.
 *   In short, the function starts processing
 *   of MD segment-rows as soon as their inputs are available
 *   and the previous segment-row has completed.  At
 *   any given time, only one segment row per picture
 *   is being processed.
 *
 * The function has two functions:
 *
 * (1) Update the Segment Completion Mask which tracks
 *   which MD Segment inputs are available.
 *
 * (2) Increment the segment-row counter (current_row_idx)
 *   as the segment-rows are completed.
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
 * (b) Continued processing
 *  -Upon the completion of a segment-row, check
 *     to see if the next segment-row's inputs have
 *     become available and begin processing if so.
 *
 * On last important point is that the thread-safe
 *   code section is kept minimally short. The MUTEX
 *   should NOT be locked for the entire processing
 *   of the segment-row (b) as this would block other
 *   threads from performing an update (A).
 ******************************************************/
EbBool assign_enc_dec_segments(EncDecSegments *segmentPtr, uint16_t *segmentInOutIndex,
                               EncDecTasks *taskPtr, EbFifo *srmFifoPtr) {
    EbBool           continue_processing_flag = EB_FALSE;
    uint32_t row_segment_index = 0;
    uint32_t segment_index;
    uint32_t right_segment_index;
    uint32_t bottom_left_segment_index;

    int16_t feedback_row_index = -1;

    uint32_t self_assigned = EB_FALSE;

    //static FILE *trace = 0;
    //
    //if(trace == 0) {
    //    trace = fopen("seg-trace.txt","w");
    //}

    switch (taskPtr->input_type) {
    case ENCDEC_TASKS_MDC_INPUT:

        // The entire picture is provided by the MDC process, so
        //   no logic is necessary to clear input dependencies.

        // Start on Segment 0 immediately
        *segmentInOutIndex  = segmentPtr->row_array[0].current_seg_index;
        taskPtr->input_type = ENCDEC_TASKS_CONTINUE;
        ++segmentPtr->row_array[0].current_seg_index;
        continue_processing_flag = EB_TRUE;

        //fprintf(trace, "Start  Pic: %u Seg: %u\n",
        //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
        //    *segmentInOutIndex);

        break;

    case ENCDEC_TASKS_ENCDEC_INPUT:

        // Setup row_segment_index to release the in_progress token
        //row_segment_index = taskPtr->encDecSegmentRowArray[0];

        // Start on the assigned row immediately
        *segmentInOutIndex  = segmentPtr->row_array[taskPtr->enc_dec_segment_row].current_seg_index;
        taskPtr->input_type = ENCDEC_TASKS_CONTINUE;
        ++segmentPtr->row_array[taskPtr->enc_dec_segment_row].current_seg_index;
        continue_processing_flag = EB_TRUE;

        //fprintf(trace, "Start  Pic: %u Seg: %u\n",
        //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
        //    *segmentInOutIndex);

        break;

    case ENCDEC_TASKS_CONTINUE:

        // Update the Dependency List for Right and Bottom Neighbors
        segment_index     = *segmentInOutIndex;
        row_segment_index = segment_index / segmentPtr->segment_band_count;

        right_segment_index       = segment_index + 1;
        bottom_left_segment_index = segment_index + segmentPtr->segment_band_count;

        // Right Neighbor
        if (segment_index < segmentPtr->row_array[row_segment_index].ending_seg_index) {
            eb_block_on_mutex(segmentPtr->row_array[row_segment_index].assignment_mutex);

            --segmentPtr->dep_map.dependency_map[right_segment_index];

            if (segmentPtr->dep_map.dependency_map[right_segment_index] == 0) {
                *segmentInOutIndex = segmentPtr->row_array[row_segment_index].current_seg_index;
                ++segmentPtr->row_array[row_segment_index].current_seg_index;
                self_assigned            = EB_TRUE;
                continue_processing_flag = EB_TRUE;

                //fprintf(trace, "Start  Pic: %u Seg: %u\n",
                //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
                //    *segmentInOutIndex);
            }

            eb_release_mutex(segmentPtr->row_array[row_segment_index].assignment_mutex);
        }

        // Bottom-left Neighbor
        if (row_segment_index < segmentPtr->segment_row_count - 1 &&
            bottom_left_segment_index >=
                segmentPtr->row_array[row_segment_index + 1].starting_seg_index) {
            eb_block_on_mutex(segmentPtr->row_array[row_segment_index + 1].assignment_mutex);

            --segmentPtr->dep_map.dependency_map[bottom_left_segment_index];

            if (segmentPtr->dep_map.dependency_map[bottom_left_segment_index] == 0) {
                if (self_assigned == EB_TRUE)
                    feedback_row_index = (int16_t)row_segment_index + 1;
                else {
                    *segmentInOutIndex =
                        segmentPtr->row_array[row_segment_index + 1].current_seg_index;
                    ++segmentPtr->row_array[row_segment_index + 1].current_seg_index;
                    continue_processing_flag = EB_TRUE;

                    //fprintf(trace, "Start  Pic: %u Seg: %u\n",
                    //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
                    //    *segmentInOutIndex);
                }
            }
            eb_release_mutex(segmentPtr->row_array[row_segment_index + 1].assignment_mutex);
        }

        if (feedback_row_index > 0) {
            EbObjectWrapper *wrapper_ptr;
            eb_get_empty_object(srmFifoPtr, &wrapper_ptr);
            EncDecTasks *    feedback_task_ptr     = (EncDecTasks *)wrapper_ptr->object_ptr;
            feedback_task_ptr->input_type          = ENCDEC_TASKS_ENCDEC_INPUT;
            feedback_task_ptr->enc_dec_segment_row = feedback_row_index;
            feedback_task_ptr->pcs_wrapper_ptr     = taskPtr->pcs_wrapper_ptr;
            feedback_task_ptr->tile_group_index = taskPtr->tile_group_index;
            eb_post_full_object(wrapper_ptr);
        }

        break;

    default: break;
    }

    return continue_processing_flag;
}
void recon_output(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    EncodeContext *     encode_context_ptr = scs_ptr->encode_context_ptr;
    // The totalNumberOfReconFrames counter has to be write/read protected as
    //   it is used to determine the end of the stream.  If it is not protected
    //   the encoder might not properly terminate.
    eb_block_on_mutex(encode_context_ptr->total_number_of_recon_frame_mutex);

    if (!pcs_ptr->parent_pcs_ptr->is_alt_ref) {
        EbBool           is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
        EbObjectWrapper *output_recon_wrapper_ptr;
        // Get Recon Buffer
        eb_get_empty_object(scs_ptr->encode_context_ptr->recon_output_fifo_ptr,
                            &output_recon_wrapper_ptr);
        EbBufferHeaderType *output_recon_ptr = (EbBufferHeaderType *)output_recon_wrapper_ptr->object_ptr;
        output_recon_ptr->flags = 0;

        // START READ/WRITE PROTECTED SECTION
        if (encode_context_ptr->total_number_of_recon_frames ==
            encode_context_ptr->terminating_picture_number)
            output_recon_ptr->flags = EB_BUFFERFLAG_EOS;

        encode_context_ptr->total_number_of_recon_frames++;

        //eb_release_mutex(encode_context_ptr->terminating_conditions_mutex);

        // STOP READ/WRITE PROTECTED SECTION
        output_recon_ptr->n_filled_len = 0;

        // Copy the Reconstructed Picture to the Output Recon Buffer
        {
            uint32_t sample_total_count;
            uint8_t *recon_read_ptr;
            uint8_t *recon_write_ptr;

            EbPictureBufferDesc *recon_ptr;
            {
                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    recon_ptr = is_16bit ? ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr
                                                ->reference_picture_wrapper_ptr->object_ptr)
                                               ->reference_picture16bit
                                         : ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr
                                                ->reference_picture_wrapper_ptr->object_ptr)
                                               ->reference_picture;
                else {
                    if (is_16bit)
                        recon_ptr = pcs_ptr->recon_picture16bit_ptr;
                    else
                        recon_ptr = pcs_ptr->recon_picture_ptr;
                }
            }

            // FGN: Create a buffer if needed, copy the reconstructed picture and run the film grain synthesis algorithm

            if (scs_ptr->seq_header.film_grain_params_present) {
                EbPictureBufferDesc *intermediate_buffer_ptr;
                {
                    if (is_16bit)
                        intermediate_buffer_ptr = pcs_ptr->film_grain_picture16bit_ptr;
                    else
                        intermediate_buffer_ptr = pcs_ptr->film_grain_picture_ptr;
                }

                AomFilmGrain *film_grain_ptr;

                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    film_grain_ptr =
                        &((EbReferenceObject *)
                              pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                             ->film_grain_params;
                else
                    film_grain_ptr = &pcs_ptr->parent_pcs_ptr->frm_hdr.film_grain_params;

                eb_av1_add_film_grain(recon_ptr, intermediate_buffer_ptr, film_grain_ptr);
                recon_ptr = intermediate_buffer_ptr;
            }

            // End running the film grain
            // Y Recon Samples
            sample_total_count = ((recon_ptr->max_width - scs_ptr->max_input_pad_right) *
                                  (recon_ptr->max_height - scs_ptr->max_input_pad_bottom))
                                 << is_16bit;
            recon_read_ptr = recon_ptr->buffer_y +
                             (recon_ptr->origin_y << is_16bit) * recon_ptr->stride_y +
                             (recon_ptr->origin_x << is_16bit);
            recon_write_ptr = &(output_recon_ptr->p_buffer[output_recon_ptr->n_filled_len]);

            CHECK_REPORT_ERROR((output_recon_ptr->n_filled_len + sample_total_count <=
                                output_recon_ptr->n_alloc_len),
                               encode_context_ptr->app_callback_ptr,
                               EB_ENC_ROB_OF_ERROR);

            // Initialize Y recon buffer
            picture_copy_kernel(recon_read_ptr,
                                recon_ptr->stride_y,
                                recon_write_ptr,
                                recon_ptr->max_width - scs_ptr->max_input_pad_right,
                                recon_ptr->width - scs_ptr->pad_right,
                                recon_ptr->height - scs_ptr->pad_bottom,
                                1 << is_16bit);

            output_recon_ptr->n_filled_len += sample_total_count;

            // U Recon Samples
            sample_total_count = ((recon_ptr->max_width - scs_ptr->max_input_pad_right) *
                                      (recon_ptr->max_height - scs_ptr->max_input_pad_bottom) >>
                                  2)
                                 << is_16bit;
            recon_read_ptr = recon_ptr->buffer_cb +
                             ((recon_ptr->origin_y << is_16bit) >> 1) * recon_ptr->stride_cb +
                             ((recon_ptr->origin_x << is_16bit) >> 1);
            recon_write_ptr = &(output_recon_ptr->p_buffer[output_recon_ptr->n_filled_len]);

            CHECK_REPORT_ERROR((output_recon_ptr->n_filled_len + sample_total_count <=
                                output_recon_ptr->n_alloc_len),
                               encode_context_ptr->app_callback_ptr,
                               EB_ENC_ROB_OF_ERROR);

            // Initialize U recon buffer
            picture_copy_kernel(recon_read_ptr,
                                recon_ptr->stride_cb,
                                recon_write_ptr,
                                (recon_ptr->max_width - scs_ptr->max_input_pad_right) >> 1,
                                (recon_ptr->width - scs_ptr->pad_right) >> 1,
                                (recon_ptr->height - scs_ptr->pad_bottom) >> 1,
                                1 << is_16bit);
            output_recon_ptr->n_filled_len += sample_total_count;

            // V Recon Samples
            sample_total_count = ((recon_ptr->max_width - scs_ptr->max_input_pad_right) *
                                      (recon_ptr->max_height - scs_ptr->max_input_pad_bottom) >>
                                  2)
                                 << is_16bit;
            recon_read_ptr = recon_ptr->buffer_cr +
                             ((recon_ptr->origin_y << is_16bit) >> 1) * recon_ptr->stride_cr +
                             ((recon_ptr->origin_x << is_16bit) >> 1);
            recon_write_ptr = &(output_recon_ptr->p_buffer[output_recon_ptr->n_filled_len]);

            CHECK_REPORT_ERROR((output_recon_ptr->n_filled_len + sample_total_count <=
                                output_recon_ptr->n_alloc_len),
                               encode_context_ptr->app_callback_ptr,
                               EB_ENC_ROB_OF_ERROR);

            // Initialize V recon buffer

            picture_copy_kernel(recon_read_ptr,
                                recon_ptr->stride_cr,
                                recon_write_ptr,
                                (recon_ptr->max_width - scs_ptr->max_input_pad_right) >> 1,
                                (recon_ptr->width - scs_ptr->pad_right) >> 1,
                                (recon_ptr->height - scs_ptr->pad_bottom) >> 1,
                                1 << is_16bit);
            output_recon_ptr->n_filled_len += sample_total_count;
            output_recon_ptr->pts = pcs_ptr->picture_number;
        }

        // Post the Recon object
        eb_post_full_object(output_recon_wrapper_ptr);
    } else {
        // Overlay and altref have 1 recon only, which is from overlay pictures. So the recon of the alt_ref is not sent to the application.
        // However, to hanlde the end of sequence properly, total_number_of_recon_frames is increamented
        encode_context_ptr->total_number_of_recon_frames++;
    }
    eb_release_mutex(encode_context_ptr->total_number_of_recon_frame_mutex);
}

//************************************/
// Calculate Frame SSIM
/************************************/

void aom_ssim_parms_8x8_c(const uint8_t *s, int sp, const uint8_t *r, int rp,
                          uint32_t *sum_s, uint32_t *sum_r, uint32_t *sum_sq_s,
                          uint32_t *sum_sq_r, uint32_t *sum_sxr) {
  int i, j;
  for (i = 0; i < 8; i++, s += sp, r += rp) {
    for (j = 0; j < 8; j++) {
      *sum_s += s[j];
      *sum_r += r[j];
      *sum_sq_s += s[j] * s[j];
      *sum_sq_r += r[j] * r[j];
      *sum_sxr += s[j] * r[j];
    }
  }
}

void aom_highbd_ssim_parms_8x8_c(const uint8_t *s, int sp, const uint8_t *sinc, int spinc, const uint16_t *r,
                                 int rp, uint32_t *sum_s, uint32_t *sum_r,
                                 uint32_t *sum_sq_s, uint32_t *sum_sq_r,
                                 uint32_t *sum_sxr) {
  int i, j;
  uint32_t ss;
  for (i = 0; i < 8; i++, s += sp, sinc += spinc, r += rp) {
    for (j = 0; j < 8; j++) {
      ss = (int64_t)(s[j] << 2) + ((sinc[j]>>6)&0x3);
      *sum_s += ss;
      *sum_r += r[j];
      *sum_sq_s += ss * ss;
      *sum_sq_r += r[j] * r[j];
      *sum_sxr += ss * r[j];
    }
  }
}

static const int64_t cc1 = 26634;        // (64^2*(.01*255)^2
static const int64_t cc2 = 239708;       // (64^2*(.03*255)^2
static const int64_t cc1_10 = 428658;    // (64^2*(.01*1023)^2
static const int64_t cc2_10 = 3857925;   // (64^2*(.03*1023)^2
static const int64_t cc1_12 = 6868593;   // (64^2*(.01*4095)^2
static const int64_t cc2_12 = 61817334;  // (64^2*(.03*4095)^2

static double similarity(uint32_t sum_s, uint32_t sum_r, uint32_t sum_sq_s,
                         uint32_t sum_sq_r, uint32_t sum_sxr, int count,
                         uint32_t bd) {
  int64_t ssim_n, ssim_d;
  int64_t c1, c2;

  if (bd == 8) {
    // scale the constants by number of pixels
    c1 = (cc1 * count * count) >> 12;
    c2 = (cc2 * count * count) >> 12;
  } else if (bd == 10) {
    c1 = (cc1_10 * count * count) >> 12;
    c2 = (cc2_10 * count * count) >> 12;
  } else if (bd == 12) {
    c1 = (cc1_12 * count * count) >> 12;
    c2 = (cc2_12 * count * count) >> 12;
  } else {
    c1 = c2 = 0;
    assert(0);
  }

  ssim_n = ((int64_t)2 * sum_s * sum_r + c1) *
           ((int64_t)2 * count * sum_sxr - (int64_t)2 * sum_s * sum_r + c2);

  ssim_d = ((int64_t)sum_s * sum_s + (int64_t)sum_r * sum_r + c1) *
           ((int64_t)count * sum_sq_s - (int64_t)sum_s * sum_s +
            (int64_t)count * sum_sq_r - (int64_t)sum_r * sum_r + c2);

  return ssim_n * 1.0 / ssim_d;
}

static double ssim_8x8(const uint8_t *s, int sp, const uint8_t *r, int rp) {
  uint32_t sum_s = 0, sum_r = 0, sum_sq_s = 0, sum_sq_r = 0, sum_sxr = 0;
  aom_ssim_parms_8x8_c(s, sp, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
  return similarity(sum_s, sum_r, sum_sq_s, sum_sq_r, sum_sxr, 64, 8);
}

static double highbd_ssim_8x8(const uint8_t *s, int sp, const uint8_t *sinc, int spinc, const uint16_t *r,
                              int rp, uint32_t bd, uint32_t shift) {
  uint32_t sum_s = 0, sum_r = 0, sum_sq_s = 0, sum_sq_r = 0, sum_sxr = 0;
  aom_highbd_ssim_parms_8x8_c(s, sp, sinc, spinc, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
  return similarity(sum_s >> shift, sum_r >> shift, sum_sq_s >> (2 * shift),
                    sum_sq_r >> (2 * shift), sum_sxr >> (2 * shift), 64, bd);
}

// We are using a 8x8 moving window with starting location of each 8x8 window
// on the 4x4 pixel grid. Such arrangement allows the windows to overlap
// block boundaries to penalize blocking artifacts.
static double aom_ssim2(const uint8_t *img1, int stride_img1,
                        const uint8_t *img2, int stride_img2,
                        int width, int height) {
    int i, j;
    int samples = 0;
    double ssim_total = 0;

    // sample point start with each 4x4 location
    for (i = 0; i <= height - 8;
        i += 4, img1 += stride_img1 * 4, img2 += stride_img2 * 4) {
        for (j = 0; j <= width - 8; j += 4) {
            double v = ssim_8x8(img1 + j, stride_img1, img2 + j, stride_img2);
            ssim_total += v;
            samples++;
        }
    }
    ssim_total /= samples;
    return ssim_total;
}

static double aom_highbd_ssim2(const uint8_t *img1, int stride_img1,
                               const uint8_t *img1inc, int stride_img1inc,
                               const uint16_t *img2, int stride_img2,
                               int width, int height, uint32_t bd, uint32_t shift) {
  int i, j;
  int samples = 0;
  double ssim_total = 0;

  // sample point start with each 4x4 location
  for (i = 0; i <= height - 8;
       i += 4, img1 += stride_img1 * 4, img1inc += stride_img1inc * 4, img2 += stride_img2 * 4) {
    for (j = 0; j <= width - 8; j += 4) {
      double v = highbd_ssim_8x8((img1 + j), stride_img1,
                                 (img1inc + j), stride_img1inc,
                                 (img2 + j), stride_img2, bd,
                                 shift);
      ssim_total += v;
      samples++;
    }
  }
  ssim_total /= samples;
  return ssim_total;
}

void ssim_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory) {
    EbBool is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);

    const uint32_t ss_x = scs_ptr->subsampling_x;
    const uint32_t ss_y = scs_ptr->subsampling_y;

    if (!is_16bit) {
        EbPictureBufferDesc *recon_ptr;

        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_ptr = ((EbReferenceObject*)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
        else
            recon_ptr = pcs_ptr->recon_picture_ptr;

        EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc*)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;

        EbByte  input_buffer;
        EbByte  recon_coeff_buffer;

        EbByte buffer_y;
        EbByte buffer_cb;
        EbByte buffer_cr;

        double luma_ssim = 0.0;
        double cb_ssim = 0.0;
        double cr_ssim = 0.0;

        // if current source picture was temporally filtered, use an alternative buffer which stores
        // the original source picture
        if(pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE){
            buffer_y = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[0];
            buffer_cb = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[1];
            buffer_cr = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[2];
        }
        else {
            buffer_y = input_picture_ptr->buffer_y;
            buffer_cb = input_picture_ptr->buffer_cb;
            buffer_cr = input_picture_ptr->buffer_cr;
        }

        recon_coeff_buffer = &((recon_ptr->buffer_y)[recon_ptr->origin_x + recon_ptr->origin_y * recon_ptr->stride_y]);
        input_buffer = &(buffer_y[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
        luma_ssim = aom_ssim2(input_buffer, input_picture_ptr->stride_y, recon_coeff_buffer, recon_ptr->stride_y,
                              scs_ptr->seq_header.max_frame_width, scs_ptr->seq_header.max_frame_height);

        recon_coeff_buffer = &((recon_ptr->buffer_cb)[recon_ptr->origin_x / 2 + recon_ptr->origin_y / 2 * recon_ptr->stride_cb]);
        input_buffer = &(buffer_cb[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
        cb_ssim = aom_ssim2(input_buffer, input_picture_ptr->stride_cb, recon_coeff_buffer, recon_ptr->stride_cb,
                            scs_ptr->chroma_width, scs_ptr->chroma_height);

        recon_coeff_buffer = &((recon_ptr->buffer_cr)[recon_ptr->origin_x / 2 + recon_ptr->origin_y / 2 * recon_ptr->stride_cr]);
        input_buffer = &(buffer_cr[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
        cr_ssim = aom_ssim2(input_buffer, input_picture_ptr->stride_cr, recon_coeff_buffer, recon_ptr->stride_cr,
                            scs_ptr->chroma_width, scs_ptr->chroma_height);

        pcs_ptr->parent_pcs_ptr->luma_ssim = luma_ssim;
        pcs_ptr->parent_pcs_ptr->cb_ssim = cb_ssim;
        pcs_ptr->parent_pcs_ptr->cr_ssim = cr_ssim;

        if (free_memory && pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            EB_FREE_ARRAY(buffer_y);
            EB_FREE_ARRAY(buffer_cb);
            EB_FREE_ARRAY(buffer_cr);
        }
    }
    else {
      EbPictureBufferDesc *recon_ptr;

        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_ptr = ((EbReferenceObject*)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            recon_ptr = pcs_ptr->recon_picture16bit_ptr;
        EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc*)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;

        EbByte    input_buffer;
        uint16_t *recon_coeff_buffer;

        double luma_ssim = 0.0;
        double cb_ssim = 0.0;
        double cr_ssim = 0.0;

        if (scs_ptr->static_config.ten_bit_format == 1) {

            /* SSIM calculation for compressed 10-bit format has not been verified and debugged,
               since this format is not supported elsewhere in this version. See verify_settings(),
               which exits with an error if compressed 10-bit format is enabled. To avoid
               extra complexity of unpacking into a temporary buffer, or having to write
               new core SSIM functions, we ignore the two least signifcant bits in this
               case, and set these to zero. One test shows a difference in SSIM
               of 0.00085 setting the two least significant bits to zero. */

            const uint32_t luma_width        = input_picture_ptr->width - scs_ptr->max_input_pad_right;
            const uint32_t luma_height       = input_picture_ptr->height - scs_ptr->max_input_pad_bottom;
            const uint32_t chroma_width      = luma_width >> ss_x;
            const uint32_t pic_width_in_sb   = (luma_width + 64 - 1) / 64;
            const uint32_t pic_height_in_sb  = (luma_height + 64 - 1) / 64;
            const uint32_t chroma_height     = luma_height >> ss_y;
            uint32_t       sb_num_in_height, sb_num_in_width, bd, shift;
            uint8_t        zero_buffer[64*64];

            bd = 10;
            shift = 0 ; // both input and output are 10 bit (bitdepth - input_bd)
            memset(&zero_buffer[0], 0, sizeof(uint8_t)*64*64);

            EbByte input_buffer_org =
                &((input_picture_ptr
                       ->buffer_y)[input_picture_ptr->origin_x +
                                   input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
            uint16_t *recon_buffer_org = (uint16_t *)(&(
                (recon_ptr->buffer_y)[(recon_ptr->origin_x << is_16bit) +
                                      (recon_ptr->origin_y << is_16bit) * recon_ptr->stride_y]));
            ;

            EbByte input_buffer_org_u = &(
                (input_picture_ptr
                     ->buffer_cb)[input_picture_ptr->origin_x / 2 +
                                  input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
            ;
            uint16_t *recon_buffer_org_u = recon_coeff_buffer =
                (uint16_t *)(&((recon_ptr->buffer_cb)[(recon_ptr->origin_x << is_16bit) / 2 +
                                                      (recon_ptr->origin_y << is_16bit) / 2 *
                                                          recon_ptr->stride_cb]));
            ;

            EbByte input_buffer_org_v = &(
                (input_picture_ptr
                     ->buffer_cr)[input_picture_ptr->origin_x / 2 +
                                  input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            ;
            uint16_t *recon_buffer_org_v = recon_coeff_buffer =
                (uint16_t *)(&((recon_ptr->buffer_cr)[(recon_ptr->origin_x << is_16bit) / 2 +
                                                      (recon_ptr->origin_y << is_16bit) / 2 *
                                                          recon_ptr->stride_cr]));
            ;

           for (sb_num_in_height = 0; sb_num_in_height < pic_height_in_sb; ++sb_num_in_height) {
                for (sb_num_in_width = 0; sb_num_in_width < pic_width_in_sb; ++sb_num_in_width) {
                    uint32_t tb_origin_x = sb_num_in_width * 64;
                    uint32_t tb_origin_y = sb_num_in_height * 64;
                    uint32_t sb_width =
                        (luma_width - tb_origin_x) < 64 ? (luma_width - tb_origin_x) : 64;
                    uint32_t sb_height =
                        (luma_height - tb_origin_y) < 64 ? (luma_height - tb_origin_y) : 64;

                    input_buffer =
                        input_buffer_org + tb_origin_y * input_picture_ptr->stride_y + tb_origin_x;
                    recon_coeff_buffer =
                        recon_buffer_org + tb_origin_y * recon_ptr->stride_y + tb_origin_x;

                    luma_ssim += aom_highbd_ssim2(input_buffer, input_picture_ptr->stride_y, &zero_buffer[0], 64,
                                                  recon_coeff_buffer, recon_ptr->stride_y, sb_width, sb_height, bd, shift);

                    //U+V
                    tb_origin_x = sb_num_in_width * 32;
                    tb_origin_y = sb_num_in_height * 32;
                    sb_width =
                        (chroma_width - tb_origin_x) < 32 ? (chroma_width - tb_origin_x) : 32;
                    sb_height =
                        (chroma_height - tb_origin_y) < 32 ? (chroma_height - tb_origin_y) : 32;

                    input_buffer = input_buffer_org_u + tb_origin_y * input_picture_ptr->stride_cb +
                                   tb_origin_x;
                    recon_coeff_buffer =
                        recon_buffer_org_u + tb_origin_y * recon_ptr->stride_cb + tb_origin_x;

                    cb_ssim += aom_highbd_ssim2(input_buffer, input_picture_ptr->stride_cb, &zero_buffer[0], 64,
                                                recon_coeff_buffer, recon_ptr->stride_cb, sb_width, sb_height, bd, shift);

                    input_buffer = input_buffer_org_v + tb_origin_y * input_picture_ptr->stride_cr +
                                   tb_origin_x;
                    recon_coeff_buffer =
                        recon_buffer_org_v + tb_origin_y * recon_ptr->stride_cr + tb_origin_x;

                    cr_ssim += aom_highbd_ssim2(input_buffer, input_picture_ptr->stride_cr, &zero_buffer[0], 64,
                                                recon_coeff_buffer, recon_ptr->stride_cr, sb_width, sb_height, bd, shift);
                }
            }

            luma_ssim /= pic_height_in_sb * pic_width_in_sb;
            cb_ssim   /= pic_height_in_sb * pic_width_in_sb;
            cr_ssim   /= pic_height_in_sb * pic_width_in_sb;

            pcs_ptr->parent_pcs_ptr->luma_ssim = luma_ssim;
            pcs_ptr->parent_pcs_ptr->cb_ssim = cb_ssim;
            pcs_ptr->parent_pcs_ptr->cr_ssim = cr_ssim;
        }
        else {
            recon_coeff_buffer = (uint16_t*)(&((recon_ptr->buffer_y)[(recon_ptr->origin_x << is_16bit) + (recon_ptr->origin_y << is_16bit) * recon_ptr->stride_y]));

            // if current source picture was temporally filtered, use an alternative buffer which stores
            // the original source picture
            EbByte buffer_y, buffer_bit_inc_y;
            EbByte buffer_cb, buffer_bit_inc_cb;
            EbByte buffer_cr, buffer_bit_inc_cr;
            int bd, shift;

            if(pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE){
                buffer_y = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[0];
                buffer_bit_inc_y = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[0];
                buffer_cb = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[1];
                buffer_bit_inc_cb = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[1];
                buffer_cr = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[2];
                buffer_bit_inc_cr = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[2];
            }else{
                buffer_y = input_picture_ptr->buffer_y;
                buffer_bit_inc_y = input_picture_ptr->buffer_bit_inc_y;
                buffer_cb = input_picture_ptr->buffer_cb;
                buffer_bit_inc_cb = input_picture_ptr->buffer_bit_inc_cb;
                buffer_cr = input_picture_ptr->buffer_cr;
                buffer_bit_inc_cr = input_picture_ptr->buffer_bit_inc_cr;
            }

            bd = 10;
            shift = 0 ; // both input and output are 10 bit (bitdepth - input_bd)

            input_buffer = &((buffer_y)[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
            EbByte input_buffer_bit_inc = &(
                (buffer_bit_inc_y)[input_picture_ptr->origin_x +
                                   input_picture_ptr->origin_y *
                                       input_picture_ptr->stride_bit_inc_y]);
            luma_ssim = aom_highbd_ssim2(input_buffer, input_picture_ptr->stride_y, input_buffer_bit_inc, input_picture_ptr->stride_bit_inc_y,
                                         recon_coeff_buffer, recon_ptr->stride_y, scs_ptr->seq_header.max_frame_width, scs_ptr->seq_header.max_frame_height, bd, shift);

            recon_coeff_buffer = (uint16_t*)(&((recon_ptr->buffer_cb)[(recon_ptr->origin_x << is_16bit) / 2 + (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cb]));
            input_buffer = &((buffer_cb)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
            input_buffer_bit_inc = &((buffer_bit_inc_cb)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_bit_inc_cb]);
            cb_ssim = aom_highbd_ssim2(input_buffer, input_picture_ptr->stride_cb, input_buffer_bit_inc, input_picture_ptr->stride_bit_inc_cb,
                                       recon_coeff_buffer, recon_ptr->stride_cb, scs_ptr->chroma_width, scs_ptr->chroma_height, bd, shift);

            recon_coeff_buffer = (uint16_t*)(&((recon_ptr->buffer_cr)[(recon_ptr->origin_x << is_16bit) / 2 + (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cr]));
            input_buffer = &((buffer_cr)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            input_buffer_bit_inc = &((buffer_bit_inc_cr)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_bit_inc_cr]);
            cr_ssim = aom_highbd_ssim2(input_buffer, input_picture_ptr->stride_cr, input_buffer_bit_inc, input_picture_ptr->stride_bit_inc_cr,
                                    recon_coeff_buffer, recon_ptr->stride_cr, scs_ptr->chroma_width, scs_ptr->chroma_height, bd, shift);

            pcs_ptr->parent_pcs_ptr->luma_ssim = luma_ssim;
            pcs_ptr->parent_pcs_ptr->cb_ssim = cb_ssim;
            pcs_ptr->parent_pcs_ptr->cr_ssim = cr_ssim;

            if (free_memory && pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
                EB_FREE_ARRAY(buffer_y);
                EB_FREE_ARRAY(buffer_cb);
                EB_FREE_ARRAY(buffer_cr);
                EB_FREE_ARRAY(buffer_bit_inc_y);
                EB_FREE_ARRAY(buffer_bit_inc_cb);
                EB_FREE_ARRAY(buffer_bit_inc_cr);
            }
        }
    }

}

void psnr_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory) {
    EbBool is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);

    const uint32_t ss_x = scs_ptr->subsampling_x;
    const uint32_t ss_y = scs_ptr->subsampling_y;

    if (!is_16bit) {
        EbPictureBufferDesc *recon_ptr;

        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_ptr = ((EbReferenceObject *)
                             pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                            ->reference_picture;
        else
            recon_ptr = pcs_ptr->recon_picture_ptr;

        EbPictureBufferDesc *input_picture_ptr =
            (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;

        uint64_t sse_total[3] = {0};
        uint64_t residual_distortion = 0;
        EbByte   input_buffer;
        EbByte   recon_coeff_buffer;

        EbByte buffer_y;
        EbByte buffer_cb;
        EbByte buffer_cr;

        // if current source picture was temporally filtered, use an alternative buffer which stores
        // the original source picture
        if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            buffer_y  = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[0];
            buffer_cb = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[1];
            buffer_cr = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[2];
        } else {
            buffer_y  = input_picture_ptr->buffer_y;
            buffer_cb = input_picture_ptr->buffer_cb;
            buffer_cr = input_picture_ptr->buffer_cr;
        }

        recon_coeff_buffer = &(
            (recon_ptr->buffer_y)[recon_ptr->origin_x + recon_ptr->origin_y * recon_ptr->stride_y]);
        input_buffer = &(buffer_y[input_picture_ptr->origin_x +
                                  input_picture_ptr->origin_y * input_picture_ptr->stride_y]);

        residual_distortion = 0;

        for (int row_index = 0;
             row_index < input_picture_ptr->height - scs_ptr->max_input_pad_bottom;
             ++row_index) {
            for (int column_index = 0;
                 column_index < input_picture_ptr->width - scs_ptr->max_input_pad_right;
                 ++column_index) {
                residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                    (recon_coeff_buffer[column_index]));
            }

            input_buffer += input_picture_ptr->stride_y;
            recon_coeff_buffer += recon_ptr->stride_y;
        }

        sse_total[0] = residual_distortion;

        recon_coeff_buffer =
            &((recon_ptr->buffer_cb)[recon_ptr->origin_x / 2 +
                                     recon_ptr->origin_y / 2 * recon_ptr->stride_cb]);
        input_buffer = &(buffer_cb[input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);

        residual_distortion = 0;
        for (int row_index = 0; row_index <
             (input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y;
             ++row_index) {
            for (int column_index = 0; column_index <
                 (input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x;
                 ++column_index) {
                residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                    (recon_coeff_buffer[column_index]));
            }

            input_buffer += input_picture_ptr->stride_cb;
            recon_coeff_buffer += recon_ptr->stride_cb;
        }

        sse_total[1] = residual_distortion;

        recon_coeff_buffer =
            &((recon_ptr->buffer_cr)[recon_ptr->origin_x / 2 +
                                     recon_ptr->origin_y / 2 * recon_ptr->stride_cr]);
        input_buffer        = &(buffer_cr[input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
        residual_distortion = 0;

        for (int row_index = 0; row_index <
             (input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y;
             ++row_index) {
            for (int column_index = 0; column_index <
                 (input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x;
                 ++column_index) {
                residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                    (recon_coeff_buffer[column_index]));
            }

            input_buffer += input_picture_ptr->stride_cr;
            recon_coeff_buffer += recon_ptr->stride_cr;
        }

        sse_total[2]                      = residual_distortion;
        pcs_ptr->parent_pcs_ptr->luma_sse = (uint32_t)sse_total[0];
        pcs_ptr->parent_pcs_ptr->cb_sse   = (uint32_t)sse_total[1];
        pcs_ptr->parent_pcs_ptr->cr_sse   = (uint32_t)sse_total[2];

        if(free_memory && pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            EB_FREE_ARRAY(buffer_y);
            EB_FREE_ARRAY(buffer_cb);
            EB_FREE_ARRAY(buffer_cr);
        }
    }
    else {
        EbPictureBufferDesc *recon_ptr;

        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_ptr = ((EbReferenceObject *)
                             pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                            ->reference_picture16bit;
        else
            recon_ptr = pcs_ptr->recon_picture16bit_ptr;
        EbPictureBufferDesc *input_picture_ptr =
            (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;

        uint64_t  sse_total[3] = {0};
        uint64_t  residual_distortion = 0;
        EbByte    input_buffer;
        EbByte    input_buffer_bit_inc;
        uint16_t *recon_coeff_buffer;

        if (scs_ptr->static_config.ten_bit_format == 1) {
            const uint32_t luma_width        = input_picture_ptr->width - scs_ptr->max_input_pad_right;
            const uint32_t luma_height       = input_picture_ptr->height - scs_ptr->max_input_pad_bottom;
            const uint32_t chroma_width      = luma_width >> ss_x;
            const uint32_t pic_width_in_sb   = (luma_width + 64 - 1) / 64;
            const uint32_t pic_height_in_sb  = (luma_height + 64 - 1) / 64;
            const uint32_t luma_2bit_width   = luma_width / 4;
            const uint32_t chroma_height     = luma_height >> ss_y;
            const uint32_t chroma_2bit_width = chroma_width / 4;
            uint32_t       sb_num_in_height, sb_num_in_width;

            EbByte input_buffer_org =
                &((input_picture_ptr
                       ->buffer_y)[input_picture_ptr->origin_x +
                                   input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
            uint16_t *recon_buffer_org = (uint16_t *)(&(
                (recon_ptr->buffer_y)[(recon_ptr->origin_x << is_16bit) +
                                      (recon_ptr->origin_y << is_16bit) * recon_ptr->stride_y]));
            ;

            EbByte input_buffer_org_u = &(
                (input_picture_ptr
                     ->buffer_cb)[input_picture_ptr->origin_x / 2 +
                                  input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
            ;
            uint16_t *recon_buffer_org_u = recon_coeff_buffer =
                (uint16_t *)(&((recon_ptr->buffer_cb)[(recon_ptr->origin_x << is_16bit) / 2 +
                                                      (recon_ptr->origin_y << is_16bit) / 2 *
                                                          recon_ptr->stride_cb]));
            ;

            EbByte input_buffer_org_v = &(
                (input_picture_ptr
                     ->buffer_cr)[input_picture_ptr->origin_x / 2 +
                                  input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            ;
            uint16_t *recon_buffer_org_v = recon_coeff_buffer =
                (uint16_t *)(&((recon_ptr->buffer_cr)[(recon_ptr->origin_x << is_16bit) / 2 +
                                                      (recon_ptr->origin_y << is_16bit) / 2 *
                                                          recon_ptr->stride_cr]));
            ;

            residual_distortion            = 0;
            uint64_t residual_distortion_u = 0;
            uint64_t residual_distortion_v = 0;

            for (sb_num_in_height = 0; sb_num_in_height < pic_height_in_sb; ++sb_num_in_height) {
                for (sb_num_in_width = 0; sb_num_in_width < pic_width_in_sb; ++sb_num_in_width) {
                    uint32_t tb_origin_x = sb_num_in_width * 64;
                    uint32_t tb_origin_y = sb_num_in_height * 64;
                    uint32_t sb_width =
                        (luma_width - tb_origin_x) < 64 ? (luma_width - tb_origin_x) : 64;
                    uint32_t sb_height =
                        (luma_height - tb_origin_y) < 64 ? (luma_height - tb_origin_y) : 64;

                    input_buffer =
                        input_buffer_org + tb_origin_y * input_picture_ptr->stride_y + tb_origin_x;
                    input_buffer_bit_inc = input_picture_ptr->buffer_bit_inc_y +
                                           tb_origin_y * luma_2bit_width +
                                           (tb_origin_x / 4) * sb_height;
                    recon_coeff_buffer =
                        recon_buffer_org + tb_origin_y * recon_ptr->stride_y + tb_origin_x;

                    uint64_t j, k;
                    uint16_t out_pixel;
                    uint8_t  n_bit_pixel;
                    uint8_t  four_2bit_pels;
                    uint32_t inn_stride = sb_width / 4;

                    for (j = 0; j < sb_height; j++) {
                        for (k = 0; k < sb_width / 4; k++) {
                            four_2bit_pels = input_buffer_bit_inc[k + j * inn_stride];

                            n_bit_pixel = (four_2bit_pels >> 6) & 3;
                            out_pixel   = input_buffer[k * 4 + 0 + j * input_picture_ptr->stride_y]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 0 + j * recon_ptr->stride_y]);

                            n_bit_pixel = (four_2bit_pels >> 4) & 3;
                            out_pixel   = input_buffer[k * 4 + 1 + j * input_picture_ptr->stride_y]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 1 + j * recon_ptr->stride_y]);

                            n_bit_pixel = (four_2bit_pels >> 2) & 3;
                            out_pixel   = input_buffer[k * 4 + 2 + j * input_picture_ptr->stride_y]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 2 + j * recon_ptr->stride_y]);

                            n_bit_pixel = (four_2bit_pels >> 0) & 3;
                            out_pixel   = input_buffer[k * 4 + 3 + j * input_picture_ptr->stride_y]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 3 + j * recon_ptr->stride_y]);
                        }
                    }

                    //U+V

                    tb_origin_x = sb_num_in_width * 32;
                    tb_origin_y = sb_num_in_height * 32;
                    sb_width =
                        (chroma_width - tb_origin_x) < 32 ? (chroma_width - tb_origin_x) : 32;
                    sb_height =
                        (chroma_height - tb_origin_y) < 32 ? (chroma_height - tb_origin_y) : 32;

                    inn_stride = sb_width / 4;

                    input_buffer = input_buffer_org_u + tb_origin_y * input_picture_ptr->stride_cb +
                                   tb_origin_x;

                    input_buffer_bit_inc = input_picture_ptr->buffer_bit_inc_cb +
                                           tb_origin_y * chroma_2bit_width +
                                           (tb_origin_x / 4) * sb_height;

                    recon_coeff_buffer =
                        recon_buffer_org_u + tb_origin_y * recon_ptr->stride_cb + tb_origin_x;

                    for (j = 0; j < sb_height; j++) {
                        for (k = 0; k < sb_width / 4; k++) {
                            four_2bit_pels = input_buffer_bit_inc[k + j * inn_stride];

                            n_bit_pixel = (four_2bit_pels >> 6) & 3;
                            out_pixel   = input_buffer[k * 4 + 0 + j * input_picture_ptr->stride_cb]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion_u += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 0 + j * recon_ptr->stride_cb]);

                            n_bit_pixel = (four_2bit_pels >> 4) & 3;
                            out_pixel   = input_buffer[k * 4 + 1 + j * input_picture_ptr->stride_cb]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion_u += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 1 + j * recon_ptr->stride_cb]);

                            n_bit_pixel = (four_2bit_pels >> 2) & 3;
                            out_pixel   = input_buffer[k * 4 + 2 + j * input_picture_ptr->stride_cb]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion_u += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 2 + j * recon_ptr->stride_cb]);

                            n_bit_pixel = (four_2bit_pels >> 0) & 3;
                            out_pixel   = input_buffer[k * 4 + 3 + j * input_picture_ptr->stride_cb]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion_u += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 3 + j * recon_ptr->stride_cb]);
                        }
                    }

                    input_buffer = input_buffer_org_v + tb_origin_y * input_picture_ptr->stride_cr +
                                   tb_origin_x;
                    input_buffer_bit_inc = input_picture_ptr->buffer_bit_inc_cr +
                                           tb_origin_y * chroma_2bit_width +
                                           (tb_origin_x / 4) * sb_height;
                    recon_coeff_buffer =
                        recon_buffer_org_v + tb_origin_y * recon_ptr->stride_cr + tb_origin_x;

                    for (j = 0; j < sb_height; j++) {
                        for (k = 0; k < sb_width / 4; k++) {
                            four_2bit_pels = input_buffer_bit_inc[k + j * inn_stride];

                            n_bit_pixel = (four_2bit_pels >> 6) & 3;
                            out_pixel   = input_buffer[k * 4 + 0 + j * input_picture_ptr->stride_cr]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion_v += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 0 + j * recon_ptr->stride_cr]);

                            n_bit_pixel = (four_2bit_pels >> 4) & 3;
                            out_pixel   = input_buffer[k * 4 + 1 + j * input_picture_ptr->stride_cr]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion_v += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 1 + j * recon_ptr->stride_cr]);

                            n_bit_pixel = (four_2bit_pels >> 2) & 3;
                            out_pixel   = input_buffer[k * 4 + 2 + j * input_picture_ptr->stride_cr]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion_v += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 2 + j * recon_ptr->stride_cr]);

                            n_bit_pixel = (four_2bit_pels >> 0) & 3;
                            out_pixel   = input_buffer[k * 4 + 3 + j * input_picture_ptr->stride_cr]
                                        << 2;
                            out_pixel = out_pixel | n_bit_pixel;
                            residual_distortion_v += (int64_t)SQR(
                                (int64_t)out_pixel -
                                (int64_t)recon_coeff_buffer[k * 4 + 3 + j * recon_ptr->stride_cr]);
                        }
                    }
                }
            }

            sse_total[0] = residual_distortion;
            sse_total[1] = residual_distortion_u;
            sse_total[2] = residual_distortion_v;
        } else {
            recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr->buffer_y)[(recon_ptr->origin_x << is_16bit) +
                                      (recon_ptr->origin_y << is_16bit) * recon_ptr->stride_y]));

            // if current source picture was temporally filtered, use an alternative buffer which stores
            // the original source picture
            EbByte buffer_y, buffer_bit_inc_y;
            EbByte buffer_cb, buffer_bit_inc_cb;
            EbByte buffer_cr, buffer_bit_inc_cr;

            if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
                buffer_y          = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[0];
                buffer_bit_inc_y  = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[0];
                buffer_cb         = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[1];
                buffer_bit_inc_cb = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[1];
                buffer_cr         = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[2];
                buffer_bit_inc_cr = pcs_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[2];
            } else {
                buffer_y          = input_picture_ptr->buffer_y;
                buffer_bit_inc_y  = input_picture_ptr->buffer_bit_inc_y;
                buffer_cb         = input_picture_ptr->buffer_cb;
                buffer_bit_inc_cb = input_picture_ptr->buffer_bit_inc_cb;
                buffer_cr         = input_picture_ptr->buffer_cr;
                buffer_bit_inc_cr = input_picture_ptr->buffer_bit_inc_cr;
            }

            input_buffer         = &((buffer_y)[input_picture_ptr->origin_x +
                                        input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
            input_buffer_bit_inc = &((buffer_bit_inc_y)[input_picture_ptr->origin_x +
                                                        input_picture_ptr->origin_y *
                                                            input_picture_ptr->stride_bit_inc_y]);

            residual_distortion = 0;

            for (int row_index = 0;
                 row_index < input_picture_ptr->height - scs_ptr->max_input_pad_bottom;
                 ++row_index) {
                for (int column_index = 0; column_index <
                     input_picture_ptr->width - scs_ptr->max_input_pad_right;
                     ++column_index) {
                    residual_distortion +=
                        (int64_t)SQR((int64_t)((((input_buffer[column_index]) << 2) |
                                                ((input_buffer_bit_inc[column_index] >> 6) & 3))) -
                                     (recon_coeff_buffer[column_index]));
                }

                input_buffer += input_picture_ptr->stride_y;
                input_buffer_bit_inc += input_picture_ptr->stride_bit_inc_y;
                recon_coeff_buffer += recon_ptr->stride_y;
            }

            sse_total[0] = residual_distortion;

            recon_coeff_buffer =
                (uint16_t *)(&((recon_ptr->buffer_cb)[(recon_ptr->origin_x << is_16bit) / 2 +
                                                      (recon_ptr->origin_y << is_16bit) / 2 *
                                                          recon_ptr->stride_cb]));
            input_buffer =
                &((buffer_cb)[input_picture_ptr->origin_x / 2 +
                              input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
            input_buffer_bit_inc = &((buffer_bit_inc_cb)[input_picture_ptr->origin_x / 2 +
                                                         input_picture_ptr->origin_y / 2 *
                                                             input_picture_ptr->stride_bit_inc_cb]);

            residual_distortion = 0;
            for (int row_index = 0; row_index <
                 (input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y;
                 ++row_index) {
                for (int column_index = 0; column_index <
                     (input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x;
                     ++column_index) {
                    residual_distortion +=
                        (int64_t)SQR((int64_t)((((input_buffer[column_index]) << 2) |
                                                ((input_buffer_bit_inc[column_index] >> 6) & 3))) -
                                     (recon_coeff_buffer[column_index]));
                }

                input_buffer += input_picture_ptr->stride_cb;
                input_buffer_bit_inc += input_picture_ptr->stride_bit_inc_cb;
                recon_coeff_buffer += recon_ptr->stride_cb;
            }

            sse_total[1] = residual_distortion;

            recon_coeff_buffer =
                (uint16_t *)(&((recon_ptr->buffer_cr)[(recon_ptr->origin_x << is_16bit) / 2 +
                                                      (recon_ptr->origin_y << is_16bit) / 2 *
                                                          recon_ptr->stride_cr]));
            input_buffer =
                &((buffer_cr)[input_picture_ptr->origin_x / 2 +
                              input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            input_buffer_bit_inc = &((buffer_bit_inc_cr)[input_picture_ptr->origin_x / 2 +
                                                         input_picture_ptr->origin_y / 2 *
                                                             input_picture_ptr->stride_bit_inc_cr]);

            residual_distortion = 0;

            for (int row_index = 0; row_index <
                 (input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y;
                 ++row_index) {
                for (int column_index = 0; column_index <
                     (input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x;
                     ++column_index) {
                    residual_distortion +=
                        (int64_t)SQR((int64_t)((((input_buffer[column_index]) << 2) |
                                                ((input_buffer_bit_inc[column_index] >> 6) & 3))) -
                                     (recon_coeff_buffer[column_index]));
                }

                input_buffer += input_picture_ptr->stride_cr;
                input_buffer_bit_inc += input_picture_ptr->stride_bit_inc_cr;
                recon_coeff_buffer += recon_ptr->stride_cr;
            }

            sse_total[2] = residual_distortion;

            if (free_memory && pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
                EB_FREE_ARRAY(buffer_y);
                EB_FREE_ARRAY(buffer_cb);
                EB_FREE_ARRAY(buffer_cr);
                EB_FREE_ARRAY(buffer_bit_inc_y);
                EB_FREE_ARRAY(buffer_bit_inc_cb);
                EB_FREE_ARRAY(buffer_bit_inc_cr);
           }
        }

        pcs_ptr->parent_pcs_ptr->luma_sse = (uint32_t)sse_total[0];
        pcs_ptr->parent_pcs_ptr->cb_sse   = (uint32_t)sse_total[1];
        pcs_ptr->parent_pcs_ptr->cr_sse   = (uint32_t)sse_total[2];
    }
}

void pad_ref_and_set_flags(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    EbReferenceObject *reference_object =
        (EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
    EbPictureBufferDesc *ref_pic_ptr = (EbPictureBufferDesc *)reference_object->reference_picture;
    EbPictureBufferDesc *ref_pic_16bit_ptr =
        (EbPictureBufferDesc *)reference_object->reference_picture16bit;
    EbBool is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);

    if (!is_16bit) {
        pad_picture_to_multiple_of_min_blk_size_dimensions(scs_ptr, ref_pic_ptr);
        // Y samples
        generate_padding(ref_pic_ptr->buffer_y,
                         ref_pic_ptr->stride_y,
                         ref_pic_ptr->width,
                         ref_pic_ptr->height,
                         ref_pic_ptr->origin_x,
                         ref_pic_ptr->origin_y);

        // Cb samples
        generate_padding(ref_pic_ptr->buffer_cb,
                         ref_pic_ptr->stride_cb,
                         ref_pic_ptr->width >> 1,
                         ref_pic_ptr->height >> 1,
                         ref_pic_ptr->origin_x >> 1,
                         ref_pic_ptr->origin_y >> 1);

        // Cr samples
        generate_padding(ref_pic_ptr->buffer_cr,
                         ref_pic_ptr->stride_cr,
                         ref_pic_ptr->width >> 1,
                         ref_pic_ptr->height >> 1,
                         ref_pic_ptr->origin_x >> 1,
                         ref_pic_ptr->origin_y >> 1);
    }

    //We need this for MCP
    if (is_16bit) {
        pad_picture_to_multiple_of_min_blk_size_dimensions_16bit(scs_ptr, ref_pic_16bit_ptr);
        // Y samples
        generate_padding16_bit(ref_pic_16bit_ptr->buffer_y,
                               ref_pic_16bit_ptr->stride_y << 1,
                               ref_pic_16bit_ptr->width << 1,
                               ref_pic_16bit_ptr->height,
                               ref_pic_16bit_ptr->origin_x << 1,
                               ref_pic_16bit_ptr->origin_y);

        // Cb samples
        generate_padding16_bit(ref_pic_16bit_ptr->buffer_cb,
                               ref_pic_16bit_ptr->stride_cb << 1,
                               ref_pic_16bit_ptr->width,
                               ref_pic_16bit_ptr->height >> 1,
                               ref_pic_16bit_ptr->origin_x,
                               ref_pic_16bit_ptr->origin_y >> 1);

        // Cr samples
        generate_padding16_bit(ref_pic_16bit_ptr->buffer_cr,
                               ref_pic_16bit_ptr->stride_cr << 1,
                               ref_pic_16bit_ptr->width,
                               ref_pic_16bit_ptr->height >> 1,
                               ref_pic_16bit_ptr->origin_x,
                               ref_pic_16bit_ptr->origin_y >> 1);

        // Hsan: unpack ref samples (to be used @ MD)
        un_pack2d((uint16_t *)ref_pic_16bit_ptr->buffer_y,
                  ref_pic_16bit_ptr->stride_y,
                  ref_pic_ptr->buffer_y,
                  ref_pic_ptr->stride_y,
                  ref_pic_ptr->buffer_bit_inc_y,
                  ref_pic_ptr->stride_bit_inc_y,
                  ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1),
                  ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1));

        un_pack2d((uint16_t *)ref_pic_16bit_ptr->buffer_cb,
                  ref_pic_16bit_ptr->stride_cb,
                  ref_pic_ptr->buffer_cb,
                  ref_pic_ptr->stride_cb,
                  ref_pic_ptr->buffer_bit_inc_cb,
                  ref_pic_ptr->stride_bit_inc_cb,
                  (ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1)) >> 1,
                  (ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1)) >> 1);

        un_pack2d((uint16_t *)ref_pic_16bit_ptr->buffer_cr,
                  ref_pic_16bit_ptr->stride_cr,
                  ref_pic_ptr->buffer_cr,
                  ref_pic_ptr->stride_cr,
                  ref_pic_ptr->buffer_bit_inc_cr,
                  ref_pic_ptr->stride_bit_inc_cr,
                  (ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1)) >> 1,
                  (ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1)) >> 1);
    }
    if ((scs_ptr->static_config.is_16bit_pipeline) && (!is_16bit)) {
        // Y samples
        generate_padding16_bit(ref_pic_16bit_ptr->buffer_y,
                               ref_pic_16bit_ptr->stride_y << 1,
                               (ref_pic_16bit_ptr->width - scs_ptr->max_input_pad_right) << 1,
                               ref_pic_16bit_ptr->height - scs_ptr->max_input_pad_bottom,
                               ref_pic_16bit_ptr->origin_x << 1,
                               ref_pic_16bit_ptr->origin_y);

        // Cb samples
        generate_padding16_bit(ref_pic_16bit_ptr->buffer_cb,
                               ref_pic_16bit_ptr->stride_cb << 1,
                               (ref_pic_16bit_ptr->width - scs_ptr->max_input_pad_right),
                               (ref_pic_16bit_ptr->height - scs_ptr->max_input_pad_bottom) >> 1,
                               ref_pic_16bit_ptr->origin_x,
                               ref_pic_16bit_ptr->origin_y >> 1);

        // Cr samples
        generate_padding16_bit(ref_pic_16bit_ptr->buffer_cr,
                               ref_pic_16bit_ptr->stride_cr << 1,
                               (ref_pic_16bit_ptr->width - scs_ptr->max_input_pad_right),
                               (ref_pic_16bit_ptr->height - scs_ptr->max_input_pad_bottom) >> 1,
                               ref_pic_16bit_ptr->origin_x,
                               ref_pic_16bit_ptr->origin_y >> 1);

        // Hsan: unpack ref samples (to be used @ MD)

        //Y
        uint16_t *buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_y);
        uint8_t * buf_8bit = ref_pic_ptr->buffer_y;
        convert_16bit_to_8bit(buf_16bit,
            ref_pic_16bit_ptr->stride_y,
            buf_8bit,
            ref_pic_ptr->stride_y,
            ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1),
            ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1));

        //CB
        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cb);
        buf_8bit = ref_pic_ptr->buffer_cb;
        convert_16bit_to_8bit(buf_16bit,
            ref_pic_16bit_ptr->stride_cb,
            buf_8bit,
            ref_pic_ptr->stride_cb,
            (ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1)) >> 1,
            (ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1)) >> 1);

        //CR
        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cr);
        buf_8bit = ref_pic_ptr->buffer_cr;
        convert_16bit_to_8bit(buf_16bit,
            ref_pic_16bit_ptr->stride_cr,
            buf_8bit,
            ref_pic_ptr->stride_cr,
            (ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1)) >> 1,
            (ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1)) >> 1);
    }
    // set up the ref POC
    reference_object->ref_poc = pcs_ptr->parent_pcs_ptr->picture_number;

    // set up the QP
    reference_object->qp = (uint8_t)pcs_ptr->parent_pcs_ptr->picture_qp;

    // set up the Slice Type
    reference_object->slice_type          = pcs_ptr->parent_pcs_ptr->slice_type;
    reference_object->referenced_area_avg = pcs_ptr->parent_pcs_ptr->referenced_area_avg;
}

void copy_statistics_to_ref_obj_ect(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    pcs_ptr->intra_coded_area =
        (100 * pcs_ptr->intra_coded_area) /
        (pcs_ptr->parent_pcs_ptr->aligned_width * pcs_ptr->parent_pcs_ptr->aligned_height);
    if (pcs_ptr->slice_type == I_SLICE) pcs_ptr->intra_coded_area = 0;

    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->intra_coded_area = (uint8_t)(pcs_ptr->intra_coded_area);

    uint32_t sb_index;
    for (sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index)
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->non_moving_index_array[sb_index] =
            pcs_ptr->parent_pcs_ptr->non_moving_index_array[sb_index];
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->tmp_layer_idx = (uint8_t)pcs_ptr->temporal_layer_index;
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->is_scene_change = pcs_ptr->parent_pcs_ptr->scene_change_flag;

    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->cdef_frame_strength = pcs_ptr->parent_pcs_ptr->cdef_frame_strength;

    Av1Common *cm = pcs_ptr->parent_pcs_ptr->av1_cm;
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->sg_frame_ep = cm->sg_frame_ep;
    if (scs_ptr->mfmv_enabled) {
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->frame_type = pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type;
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->order_hint = pcs_ptr->parent_pcs_ptr->cur_order_hint;
        eb_memcpy(((EbReferenceObject *)
                    pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                   ->ref_order_hint,
               pcs_ptr->parent_pcs_ptr->ref_order_hint,
               7 * sizeof(uint32_t));
    }
}

#if OBMC_FAST
void set_obmc_controls(ModeDecisionContext *mdctxt, uint8_t obmc_mode) {

    ObmcControls*obmc_ctrls = &mdctxt->obmc_ctrls;

    switch (obmc_mode)
    {
    case 0:
        obmc_ctrls->enabled = 0;
#if ON_OFF_FEATURE_MRP
        obmc_ctrls->pme_best_ref = override_feature_level(mdctxt->mrp_level,0,0,0);
#else
        obmc_ctrls->pme_best_ref = 0;
#endif
        obmc_ctrls->me_count = 0;
#if ON_OFF_FEATURE_MRP
        obmc_ctrls->mvp_ref_count = override_feature_level(mdctxt->mrp_level,0,4,1);
#else
        obmc_ctrls->mvp_ref_count = 0;
#endif
        obmc_ctrls->near_count = 0;
        break;
    case 1:
        obmc_ctrls->enabled = 1;
#if ON_OFF_FEATURE_MRP
        obmc_ctrls->pme_best_ref = override_feature_level(mdctxt->mrp_level,0,0,0);
#else
        obmc_ctrls->pme_best_ref = 0;
#endif
        obmc_ctrls->me_count = ~0;
#if ON_OFF_FEATURE_MRP
        obmc_ctrls->mvp_ref_count = override_feature_level(mdctxt->mrp_level,4,4,1);
#else
        obmc_ctrls->mvp_ref_count = 4;
#endif
        obmc_ctrls->near_count = 3;
        break;
    case 2:
        obmc_ctrls->enabled = 1;
#if ON_OFF_FEATURE_MRP
        obmc_ctrls->pme_best_ref = override_feature_level(mdctxt->mrp_level,0,0,0);
#else
        obmc_ctrls->pme_best_ref = 0;
#endif
        obmc_ctrls->me_count = ~0;
#if ON_OFF_FEATURE_MRP
        obmc_ctrls->mvp_ref_count = override_feature_level(mdctxt->mrp_level,4,4,1);
#else
        obmc_ctrls->mvp_ref_count = 4;
#endif
        obmc_ctrls->near_count = 3;
        break;
    case 3:
        obmc_ctrls->enabled = 1;
#if ON_OFF_FEATURE_MRP
        obmc_ctrls->pme_best_ref = override_feature_level(mdctxt->mrp_level,1,0,0);
#else
        obmc_ctrls->pme_best_ref = 1;
#endif
        obmc_ctrls->me_count = 1;
#if ON_OFF_FEATURE_MRP
        obmc_ctrls->mvp_ref_count = override_feature_level(mdctxt->mrp_level,1,4,1);
#else
        obmc_ctrls->mvp_ref_count = 1;
#endif
        obmc_ctrls->near_count = 1;
        break;
    default:
        assert(0);
        break;
    }


}
#endif
#if MD_REFERENCE_MASKING
#if !SOFT_CYCLES_REDUCTION
#if PRUNING_PER_INTER_TYPE
void set_inter_inter_distortion_based_reference_pruning_controls(
    ModeDecisionContext *mdctxt, uint8_t inter_inter_distortion_based_reference_pruning_mode) {
    RefPruningControls *ref_pruning_ctrls = &mdctxt->ref_pruning_ctrls;

    switch (inter_inter_distortion_based_reference_pruning_mode) {
    case 0: ref_pruning_ctrls->inter_to_inter_pruning_enabled = 0; break;
    case 1:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;

        ref_pruning_ctrls->best_refs[PA_ME_GROUP]         = 7;
        ref_pruning_ctrls->best_refs[UNI_3x3_GROUP]       = 2;
        ref_pruning_ctrls->best_refs[BI_3x3_GROUP]        = 2;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = MR_MODE ? 7 : 0;
#else
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 7;
#endif
        ref_pruning_ctrls->best_refs[WARP_GROUP]          = 7;
        ref_pruning_ctrls->best_refs[NRST_NEAR_GROUP]     = 7;
        ref_pruning_ctrls->best_refs[PRED_ME_GROUP]       = 7;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 0;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
#else
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 0;
#endif
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 0;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 0;

        break;
    case 2:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;

        ref_pruning_ctrls->best_refs[PA_ME_GROUP]         = 6;
        ref_pruning_ctrls->best_refs[UNI_3x3_GROUP]       = 2;
        ref_pruning_ctrls->best_refs[BI_3x3_GROUP]        = 2;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 0;
#else
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 6;
#endif
        ref_pruning_ctrls->best_refs[WARP_GROUP]          = 6;
        ref_pruning_ctrls->best_refs[NRST_NEAR_GROUP]     = 6;
        ref_pruning_ctrls->best_refs[PRED_ME_GROUP]       = 6;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 0;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
#else
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 0;
#endif
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 0;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 0;
        break;
    case 3:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;

        ref_pruning_ctrls->best_refs[PA_ME_GROUP]         = 5;
        ref_pruning_ctrls->best_refs[UNI_3x3_GROUP]       = 2;
        ref_pruning_ctrls->best_refs[BI_3x3_GROUP]        = 2;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 0;
#else
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 5;
#endif
        ref_pruning_ctrls->best_refs[WARP_GROUP]          = 5;
        ref_pruning_ctrls->best_refs[NRST_NEAR_GROUP]     = 5;
        ref_pruning_ctrls->best_refs[PRED_ME_GROUP]       = 5;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 0;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
#else
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 0;
#endif
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 0;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 0;
        break;
    case 4:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;

        ref_pruning_ctrls->best_refs[PA_ME_GROUP]         = 4;
        ref_pruning_ctrls->best_refs[UNI_3x3_GROUP]       = 2;
        ref_pruning_ctrls->best_refs[BI_3x3_GROUP]        = 2;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 0;
#else
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 4;
#endif
        ref_pruning_ctrls->best_refs[WARP_GROUP]          = 4;
        ref_pruning_ctrls->best_refs[NRST_NEAR_GROUP]     = 4;
        ref_pruning_ctrls->best_refs[PRED_ME_GROUP]       = 4;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 0;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
#else
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 0;
#endif
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 0;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 0;
        break;
    case 5:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;

        ref_pruning_ctrls->best_refs[PA_ME_GROUP]         = 3;
        ref_pruning_ctrls->best_refs[UNI_3x3_GROUP]       = 2;
        ref_pruning_ctrls->best_refs[BI_3x3_GROUP]        = 2;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 0;
#else
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 3;
#endif
        ref_pruning_ctrls->best_refs[WARP_GROUP]          = 3;
        ref_pruning_ctrls->best_refs[NRST_NEAR_GROUP]     = 3;
        ref_pruning_ctrls->best_refs[PRED_ME_GROUP]       = 3;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 0;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
#else
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 0;
#endif
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 0;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 0;

        break;
    case 6:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;

        ref_pruning_ctrls->best_refs[PA_ME_GROUP]         = 2;
        ref_pruning_ctrls->best_refs[UNI_3x3_GROUP]       = 2;
        ref_pruning_ctrls->best_refs[BI_3x3_GROUP]        = 2;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 0;
#else
        ref_pruning_ctrls->best_refs[NRST_NEW_NEAR_GROUP] = 2;
#endif
        ref_pruning_ctrls->best_refs[WARP_GROUP]          = 2;
        ref_pruning_ctrls->best_refs[NRST_NEAR_GROUP]     = 2;
        ref_pruning_ctrls->best_refs[PRED_ME_GROUP]       = 2;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 0;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
#if OPTIMIZE_NEAREST_NEW_NEAR
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
#else
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 0;
#endif
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 0;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 0;

        break;
    default: assert(0); break;
    }
}
#else
void set_inter_inter_distortion_based_reference_pruning_controls(ModeDecisionContext *mdctxt, uint8_t inter_inter_distortion_based_reference_pruning_mode) {

    RefPruningControls *ref_pruning_ctrls = &mdctxt->ref_pruning_ctrls;

    switch (inter_inter_distortion_based_reference_pruning_mode)
    {
    case 0:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 0;
        ref_pruning_ctrls->best_refs = 7;
        break;
    case 1:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;
        ref_pruning_ctrls->best_refs = 6;
        break;
    case 2:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;
        ref_pruning_ctrls->best_refs = 5;
        break;
    case 3:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;
        ref_pruning_ctrls->best_refs = 4;


        break;
    case 4:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;
        ref_pruning_ctrls->best_refs = 3;


        break;
    case 5:
        ref_pruning_ctrls->inter_to_inter_pruning_enabled = 1;
        ref_pruning_ctrls->best_refs = 2;


        break;
    default:
        assert(0);
        break;
    }
}
#endif
#endif
void set_inter_intra_distortion_based_reference_pruning_controls(ModeDecisionContext *mdctxt, uint8_t inter_intra_distortion_based_reference_pruning_mode) {

    RefPruningControls *ref_pruning_ctrls = &mdctxt->ref_pruning_ctrls;

    switch (inter_intra_distortion_based_reference_pruning_mode)
    {
    case 0:
        ref_pruning_ctrls->intra_to_inter_pruning_enabled = 0;
        break;
    case 1:
        ref_pruning_ctrls->intra_to_inter_pruning_enabled = 1;
        break;
    case 2:
        ref_pruning_ctrls->intra_to_inter_pruning_enabled = 1;
        break;
    case 3:
        ref_pruning_ctrls->intra_to_inter_pruning_enabled = 1;
        break;
    default:
        assert(0);
        break;
    }
}
#endif


#if BLOCK_REDUCTION_ALGORITHM_1 || BLOCK_REDUCTION_ALGORITHM_2
void set_block_based_depth_reduction_controls(ModeDecisionContext *mdctxt, uint8_t block_based_depth_reduction_level) {

    DepthReductionCtrls *depth_reduction_ctrls = &mdctxt->depth_reduction_ctrls;

    switch (block_based_depth_reduction_level)
    {
    case 0:

        depth_reduction_ctrls->enabled = 0;
        break;

    case 1:

        depth_reduction_ctrls->enabled = 1;

        depth_reduction_ctrls->cost_sq_vs_nsq_energy_based_depth_reduction_enabled = 1;
        depth_reduction_ctrls->current_to_parent_deviation_th = 0;
        depth_reduction_ctrls->sq_to_best_nsq_deviation_th = 0;
#if !M8_CLEAN_UP
        depth_reduction_ctrls->quant_coeff_energy_th = 0;
#endif
        depth_reduction_ctrls->nsq_data_based_depth_reduction_enabled = 0;
        depth_reduction_ctrls->sq_to_4_sq_children_th = 0;
        depth_reduction_ctrls->h_v_to_h4_v4_th = 0;

        break;

    case 2:

        depth_reduction_ctrls->enabled = 1;

        depth_reduction_ctrls->cost_sq_vs_nsq_energy_based_depth_reduction_enabled = 1;
        depth_reduction_ctrls->current_to_parent_deviation_th = 0;
        depth_reduction_ctrls->sq_to_best_nsq_deviation_th = 0;
#if !M8_CLEAN_UP
        depth_reduction_ctrls->quant_coeff_energy_th = 0;
#endif
        depth_reduction_ctrls->nsq_data_based_depth_reduction_enabled = 1;
        depth_reduction_ctrls->sq_to_4_sq_children_th = 0;
        depth_reduction_ctrls->h_v_to_h4_v4_th = 0;

        break;

    default:
        assert(0);
        break;
    }
}
#endif
#if ADD_MD_NSQ_SEARCH
void md_nsq_motion_search_controls(ModeDecisionContext *mdctxt, uint8_t md_nsq_mv_search_level) {

    MdNsqMotionSearchCtrls *md_nsq_motion_search_ctrls = &mdctxt->md_nsq_motion_search_ctrls;

    switch (md_nsq_mv_search_level)
    {
    case 0:
        md_nsq_motion_search_ctrls->enabled = 0;
        break;
    case 1:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 31;
        md_nsq_motion_search_ctrls->full_pel_search_height = 31;
#if !PERFORM_SUB_PEL_MD
        md_nsq_motion_search_ctrls->perform_sub_pel = 1;
        md_nsq_motion_search_ctrls->half_pel_search_width = 3;
        md_nsq_motion_search_ctrls->half_pel_search_height = 3;
        md_nsq_motion_search_ctrls->quarter_pel_search_width = 3;
        md_nsq_motion_search_ctrls->quarter_pel_search_height = 3;
#endif
        break;

    case 2:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 15;
        md_nsq_motion_search_ctrls->full_pel_search_height = 15;
#if !PERFORM_SUB_PEL_MD
        md_nsq_motion_search_ctrls->perform_sub_pel = 1;
        md_nsq_motion_search_ctrls->half_pel_search_width = 3;
        md_nsq_motion_search_ctrls->half_pel_search_height = 3;
        md_nsq_motion_search_ctrls->quarter_pel_search_width = 3;
        md_nsq_motion_search_ctrls->quarter_pel_search_height = 3;
#endif
        break;
    case 3:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 11;
        md_nsq_motion_search_ctrls->full_pel_search_height = 11;
#if !PERFORM_SUB_PEL_MD
        md_nsq_motion_search_ctrls->perform_sub_pel = 1;
        md_nsq_motion_search_ctrls->half_pel_search_width = 3;
        md_nsq_motion_search_ctrls->half_pel_search_height = 3;
        md_nsq_motion_search_ctrls->quarter_pel_search_width = 3;
        md_nsq_motion_search_ctrls->quarter_pel_search_height = 3;
#endif
        break;
    case 4:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 7;
        md_nsq_motion_search_ctrls->full_pel_search_height = 7;
#if !PERFORM_SUB_PEL_MD
        md_nsq_motion_search_ctrls->perform_sub_pel = 1;
        md_nsq_motion_search_ctrls->half_pel_search_width = 3;
        md_nsq_motion_search_ctrls->half_pel_search_height = 3;
        md_nsq_motion_search_ctrls->quarter_pel_search_width = 3;
        md_nsq_motion_search_ctrls->quarter_pel_search_height = 3;
#endif
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if PERFORM_SUB_PEL_MD
#if REMOVE_MR_MACRO
void md_subpel_search_controls(ModeDecisionContext *mdctxt, uint8_t md_subpel_search_level, EbEncMode enc_mode) {
#else
void md_subpel_search_controls(ModeDecisionContext *mdctxt, uint8_t md_subpel_search_level) {
#endif
    MdSubPelSearchCtrls *md_subpel_search_ctrls = &mdctxt->md_subpel_search_ctrls;

    switch (md_subpel_search_level) {
    case 0: md_subpel_search_ctrls->enabled = 0; break;
    case 1:
        md_subpel_search_ctrls->enabled = 1;
        md_subpel_search_ctrls->use_ssd = 0;

        md_subpel_search_ctrls->do_4x4 = 1;
        md_subpel_search_ctrls->do_nsq = 1;

        md_subpel_search_ctrls->half_pel_search_enabled = 1;
        md_subpel_search_ctrls->half_pel_search_scan    = 0;
#if SEARCH_TOP_N
#if IMPROVE_HALF_PEL
#if REMOVE_MR_MACRO
        md_subpel_search_ctrls->half_pel_search_pos_cnt = enc_mode <= ENC_MR ? 8 : 5;
#else
        md_subpel_search_ctrls->half_pel_search_pos_cnt = MR_MODE_SUB_PEL ? 8 : 5;
#endif
#else
        md_subpel_search_ctrls->half_pel_search_pos_cnt = 5;
#endif
#endif
        md_subpel_search_ctrls->quarter_pel_search_enabled = 1;
        md_subpel_search_ctrls->quarter_pel_search_scan    = 0;
#if IMPROVE_QUARTER_PEL
#if REMOVE_MR_MACRO
        md_subpel_search_ctrls->quarter_pel_search_pos_cnt = enc_mode <= ENC_MR ? 8 : 1;
#else
        md_subpel_search_ctrls->quarter_pel_search_pos_cnt = MR_MODE_SUB_PEL ? 8 : 1;
#endif
#endif
        md_subpel_search_ctrls->eight_pel_search_enabled = 1;
        md_subpel_search_ctrls->eight_pel_search_scan    = 0;
#if IMPROVE_EIGHT_PEL
#if REMOVE_MR_MACRO
        md_subpel_search_ctrls->eight_pel_search_pos_cnt = enc_mode <= ENC_MR ? 8 : 1;
#else
        md_subpel_search_ctrls->eight_pel_search_pos_cnt = MR_MODE_SUB_PEL ? 8 : 1;
#endif
#endif
        break;
    case 2:
        md_subpel_search_ctrls->enabled = 1;
        md_subpel_search_ctrls->use_ssd = 0;

        md_subpel_search_ctrls->do_4x4 = 1;
        md_subpel_search_ctrls->do_nsq = 1;

        md_subpel_search_ctrls->half_pel_search_enabled = 1;
        md_subpel_search_ctrls->half_pel_search_scan = 0;
#if SEARCH_TOP_N
        md_subpel_search_ctrls->half_pel_search_pos_cnt = 3;
#endif

        md_subpel_search_ctrls->quarter_pel_search_enabled = 1;
        md_subpel_search_ctrls->quarter_pel_search_scan = 0;
#if IMPROVE_QUARTER_PEL
        md_subpel_search_ctrls->quarter_pel_search_pos_cnt = 1;
#endif
        md_subpel_search_ctrls->eight_pel_search_enabled = 1;
        md_subpel_search_ctrls->eight_pel_search_scan = 0;
#if IMPROVE_EIGHT_PEL
        md_subpel_search_ctrls->eight_pel_search_pos_cnt = 1;
#endif
        break;
    case 3:
        md_subpel_search_ctrls->enabled = 1;
        md_subpel_search_ctrls->use_ssd = 0;

        md_subpel_search_ctrls->do_4x4 = 1;
        md_subpel_search_ctrls->do_nsq = 1;

        md_subpel_search_ctrls->half_pel_search_enabled = 1;
        md_subpel_search_ctrls->half_pel_search_scan    = 0;
#if SEARCH_TOP_N
        md_subpel_search_ctrls->half_pel_search_pos_cnt = 1;
#endif

        md_subpel_search_ctrls->quarter_pel_search_enabled = 1;
        md_subpel_search_ctrls->quarter_pel_search_scan    = 0;
#if IMPROVE_QUARTER_PEL
        md_subpel_search_ctrls->quarter_pel_search_pos_cnt = 1;
#endif
        md_subpel_search_ctrls->eight_pel_search_enabled = 1;
        md_subpel_search_ctrls->eight_pel_search_scan    = 0;
#if IMPROVE_EIGHT_PEL
        md_subpel_search_ctrls->eight_pel_search_pos_cnt = 1;
#endif
        break;
    case 4:
        md_subpel_search_ctrls->enabled = 1;
        md_subpel_search_ctrls->use_ssd = 0;

        md_subpel_search_ctrls->do_4x4 = 1;
        md_subpel_search_ctrls->do_nsq = 1;

        md_subpel_search_ctrls->half_pel_search_enabled = 1;
        md_subpel_search_ctrls->half_pel_search_scan    = 1;
#if SEARCH_TOP_N
        md_subpel_search_ctrls->half_pel_search_pos_cnt = 1;
#endif
        md_subpel_search_ctrls->quarter_pel_search_enabled = 1;
        md_subpel_search_ctrls->quarter_pel_search_scan    = 1;
#if IMPROVE_QUARTER_PEL
        md_subpel_search_ctrls->quarter_pel_search_pos_cnt = 1;
#endif
        md_subpel_search_ctrls->eight_pel_search_enabled = 0;

        break;

    default: assert(0); break;
    }

    // Common setting(s)
    md_subpel_search_ctrls->half_pel_search_width = 3;
    md_subpel_search_ctrls->half_pel_search_height = 3;
    md_subpel_search_ctrls->half_pel_interpolation = 0;

    md_subpel_search_ctrls->quarter_pel_search_width = 3;
    md_subpel_search_ctrls->quarter_pel_search_height = 3;
    md_subpel_search_ctrls->quarter_pel_interpolation = 0;

    md_subpel_search_ctrls->eight_pel_search_width = 3;
    md_subpel_search_ctrls->eight_pel_search_height = 3;
    md_subpel_search_ctrls->eight_pel_interpolation = 0;

}
#endif
#if SB_CLASSIFIER
/******************************************************
* Derive SB classifier thresholds
******************************************************/
#if !CLEANUP_CYCLE_ALLOCATION
#if NEW_CYCLE_ALLOCATION
#if MULTI_BAND_ACTIONS
void set_sb_class_controls(ModeDecisionContext *context_ptr) {
    SbClassControls *sb_class_ctrls = &context_ptr->sb_class_ctrls;
    for (uint8_t sb_class_idx = 0; sb_class_idx < NUMBER_OF_SB_CLASS; sb_class_idx++)
        sb_class_ctrls->sb_class_th[sb_class_idx] = 20;
    switch (context_ptr->coeffcients_area_based_cycles_allocation_level) {
    case 0:
#if NON_UNIFORM_NSQ_BANDING
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 100;
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 100;
#else
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 20;  // 1.    [95%;100%]
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 20;  // 2.    [90%;95%]
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 20;  // 3.    [85%;90%]
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 20;  // 4     [80%;85%]
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 20;  // 5     [75%;80%]
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 20;  // 6.    [70%;75%]
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 20;  // 7.    [65%;70%]
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 20;  // 8.    [60%;65%]
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 20;  // 9.    [55%;60%]
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 20;  //10.   [50%;55%]
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 20;  //11.   [45%;50%]
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 20;  //12.   [40%;45%]
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 20;  //13.   [35%;40%]
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 20;  //14.   [30%;35%]
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 20;  //15.   [25%;30%]
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 20;  //16.   [20%;25%]
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 20;  //17    [15%;20%]
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 20;  //18.   [10%;15%]
        sb_class_ctrls->sb_class_th[SB_CLASS_19] = 20;  //19    [5%;10%]
        sb_class_ctrls->sb_class_th[SB_CLASS_20] = 20;  //20    [0%;5%]
#endif
        break;
    case 1: // TH 80%
#if NON_UNIFORM_NSQ_BANDING
    case 2: // TH 70%
    case 3: // TH 60%
    case 4: // TH 50%
    case 5: // TH 40%
    case 6: // TH 30%
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 85;
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 75;
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 65;
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 60;
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 55;
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 50;
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 45;
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 40;
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 35;
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 30;
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 25;
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 20;
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 17;
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 14;
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 10;
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 6;
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 3;
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 0;
#else
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 19;  // 1.    [95%;100%]
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 18;  // 2.    [90%;95%]
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 17;  // 3.    [85%;90%]
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 16;  // 4.    [80%;85%]
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 15;  // 5.    [75%;80%]
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 14;  // 6.    [70%;75%]
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 13;  // 7.    [65%;70%]
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 12;  // 8.    [60%;65%]
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 11;  // 9.    [55%;60%]
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 10;  //10    .[50%;55%]
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 9;  //11.    [45%;50%]
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 8;  //12.    [40%;45%]
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 7;  //13.    [35%;40%]
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 6;  //14.    [30%;35%]
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 5;  //15.    [25%;30%]
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 4;  //16.    [20%;25%]
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 3;  //17.    [15%;20%]
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 2;  //18.    [10%;15%]
        sb_class_ctrls->sb_class_th[SB_CLASS_19] = 1;  //19.    [5%;10%]
        sb_class_ctrls->sb_class_th[SB_CLASS_20] = 0;  //20.    [0%;5%]
        break;
    case 2: // TH 70%
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 19;  // 1.    [95%;100%]
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 18;  // 2.    [90%;95%]
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 17;  // 3.    [85%;90%]
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 16;  // 4.    [80%;85%]
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 15;  // 5.    [75%;80%]
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 14;  // 6.    [70%;75%]
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 13;  // 7.    [65%;70%]
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 12;  // 8.    [60%;65%]
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 11;  // 9.    [55%;60%]
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 10;  //10    .[50%;55%]
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 9;  //11.    [45%;50%]
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 8;  //12.    [40%;45%]
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 7;  //13.    [35%;40%]
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 6;  //14.    [30%;35%]
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 5;  //15.    [25%;30%]
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 4;  //16.    [20%;25%]
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 3;  //17.    [15%;20%]
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 2;  //18.    [10%;15%]
        sb_class_ctrls->sb_class_th[SB_CLASS_19] = 1;  //19.    [5%;10%]
        sb_class_ctrls->sb_class_th[SB_CLASS_20] = 0;  //20.    [0%;5%]
        break;
    case 3: // TH 60%
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 19;  // 1.    [95%;100%]
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 18;  // 2.    [90%;95%]
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 17;  // 3.    [85%;90%]
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 16;  // 4.    [80%;85%]
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 15;  // 5.    [75%;80%]
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 14;  // 6.    [70%;75%]
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 13;  // 7.    [65%;70%]
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 12;  // 8.    [60%;65%]
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 11;  // 9.    [55%;60%]
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 10;  //10    .[50%;55%]
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 9;  //11.    [45%;50%]
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 8;  //12.    [40%;45%]
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 7;  //13.    [35%;40%]
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 6;  //14.    [30%;35%]
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 5;  //15.    [25%;30%]
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 4;  //16.    [20%;25%]
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 3;  //17.    [15%;20%]
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 2;  //18.    [10%;15%]
        sb_class_ctrls->sb_class_th[SB_CLASS_19] = 1;  //19.    [5%;10%]
        sb_class_ctrls->sb_class_th[SB_CLASS_20] = 0;  //20.    [0%;5%]
        break;
    case 4: // TH 50%
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 19;  // 1.     [95%;100%]
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 18;  // 2.     [90%;95%]
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 17;  // 3.     [85%;90%]
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 16;  // 4.     [80%;85%]
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 15;  // 5.     [75%;80%]
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 14;  // 6.     [70%;75%]
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 13;  // 7.     [65%;70%]
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 12;  // 8.     [60%;65%]
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 11;  // 9.     [55%;60%]
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 10;  //10     .[50%;55%]
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 9;  //11.     [45%;50%]
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 8;  //12.     [40%;45%]
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 7;  //13.     [35%;40%]
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 6;  //14.     [30%;35%]
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 5;  //15.     [25%;30%]
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 4;  //16.     [20%;25%]
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 3;  //17.     [15%;20%]
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 2;  //18.     [10%;15%]
        sb_class_ctrls->sb_class_th[SB_CLASS_19] = 1;  //19.     [5%;10%]
        sb_class_ctrls->sb_class_th[SB_CLASS_20] = 0;  //20.     [0%;5%]
        break;
    case 5: // TH 40%
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 19;  // 1.    [95%;100%]
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 18;  // 2.    [90%;95%]
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 17;  // 3.    [85%;90%]
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 16;  // 4.    [80%;85%]
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 15;  // 5.    [75%;80%]
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 14;  // 6.    [70%;75%]
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 13;  // 7.    [65%;70%]
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 12;  // 8.    [60%;65%]
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 11;  // 9.    [55%;60%]
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 10;  //10    .[50%;55%]
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 9;  //11.    [45%;50%]
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 8;  //12.    [40%;45%]
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 7;  //13.    [35%;40%]
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 6;  //14.    [30%;35%]
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 5;  //15.    [25%;30%]
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 4;  //16.    [20%;25%]
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 3;  //17.    [15%;20%]
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 2;  //18.    [10%;15%]
        sb_class_ctrls->sb_class_th[SB_CLASS_19] = 1;  //19.    [5%;10%]
        sb_class_ctrls->sb_class_th[SB_CLASS_20] = 0;  //20.    [0%;5%]
        break;
    case 6:
        sb_class_ctrls->sb_class_th[SB_CLASS_1] = 19;  // 1.    [95%;100%]
        sb_class_ctrls->sb_class_th[SB_CLASS_2] = 18;  // 2.    [90%;95%]
        sb_class_ctrls->sb_class_th[SB_CLASS_3] = 17;  // 3.    [85%;90%]
        sb_class_ctrls->sb_class_th[SB_CLASS_4] = 16;  // 4.    [80%;85%]
        sb_class_ctrls->sb_class_th[SB_CLASS_5] = 15;  // 5.    [75%;80%]
        sb_class_ctrls->sb_class_th[SB_CLASS_6] = 14;  // 6.    [70%;75%]
        sb_class_ctrls->sb_class_th[SB_CLASS_7] = 13;  // 7.    [65%;70%]
        sb_class_ctrls->sb_class_th[SB_CLASS_8] = 12;  // 8.    [60%;65%]
        sb_class_ctrls->sb_class_th[SB_CLASS_9] = 11;  // 9.    [55%;60%]
        sb_class_ctrls->sb_class_th[SB_CLASS_10] = 10;  //10    .[50%;55%]
        sb_class_ctrls->sb_class_th[SB_CLASS_11] = 9;  //11.    [45%;50%]
        sb_class_ctrls->sb_class_th[SB_CLASS_12] = 8;  //12.    [40%;45%]
        sb_class_ctrls->sb_class_th[SB_CLASS_13] = 7;  //13.    [35%;40%]
        sb_class_ctrls->sb_class_th[SB_CLASS_14] = 6;  //14.    [30%;35%]
        sb_class_ctrls->sb_class_th[SB_CLASS_15] = 5;  //15.    [25%;30%]
        sb_class_ctrls->sb_class_th[SB_CLASS_16] = 4;  //16.    [20%;25%]
        sb_class_ctrls->sb_class_th[SB_CLASS_17] = 3;  //17.    [15%;20%]
        sb_class_ctrls->sb_class_th[SB_CLASS_18] = 2;  //18.    [10%;15%]
        sb_class_ctrls->sb_class_th[SB_CLASS_19] = 1;  //19.    [5%;10%]
        sb_class_ctrls->sb_class_th[SB_CLASS_20] = 0;  //20.    [0%;5%]
#endif
        break;
    default:
        printf("Error - Invalid sb_class_level");
        break;
    }
}
#else
void set_sb_class_controls(ModeDecisionContext *context_ptr) {
    SbClassControls *sb_class_ctrls = &context_ptr->sb_class_ctrls;
    for (uint8_t sb_class_idx = 0; sb_class_idx < NUMBER_OF_SB_CLASS; sb_class_idx++)
        sb_class_ctrls->sb_class_th[sb_class_idx] = 20;
    switch (context_ptr->coeffcients_area_based_cycles_allocation_level) {
    case 0:
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 20;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 20;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 20;
        sb_class_ctrls->sb_class_th[VERY_LOW_COMPLEX_CLASS] = 20;
        break;
    case 1: // TH 80%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 18;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 16;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 14;
        sb_class_ctrls->sb_class_th[VERY_LOW_COMPLEX_CLASS] = 0;
        break;
    case 2: // TH 70%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 16;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 14;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 10;
        sb_class_ctrls->sb_class_th[VERY_LOW_COMPLEX_CLASS] = 0;
        break;
    case 3: // TH 60%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 14;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 12;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 8;
        sb_class_ctrls->sb_class_th[VERY_LOW_COMPLEX_CLASS] = 0;
        break;
    case 4: // TH 50%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 12;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 10;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 6;
        sb_class_ctrls->sb_class_th[VERY_LOW_COMPLEX_CLASS] = 0;
        break;
    case 5: // TH 40%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 10;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 8;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 4;
        sb_class_ctrls->sb_class_th[VERY_LOW_COMPLEX_CLASS] = 0;
        break;
    case 6:
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 7;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 6;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 2;
        sb_class_ctrls->sb_class_th[VERY_LOW_COMPLEX_CLASS] = 0;
        break;
    default:
        printf("Error - Invalid sb_class_level");
        break;
    }
}
#endif
#else
void set_sb_class_controls(ModeDecisionContext *context_ptr) {
    SbClassControls *sb_class_ctrls = &context_ptr->sb_class_ctrls;
    for (uint8_t sb_class_idx = 0; sb_class_idx < NUMBER_OF_SB_CLASS; sb_class_idx++)
        sb_class_ctrls->sb_class_th[sb_class_idx] = 20;
    switch (context_ptr->coeffcients_area_based_cycles_allocation_level) {
    case 0:
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 20;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 20;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 20;
        break;
    case 1: // TH 80%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 16;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 14;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 12;
        break;
    case 2: // TH 70%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 14;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 12;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 10;
        break;
    case 3: // TH 60%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 12;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 10;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 8;
        break;
    case 4: // TH 50%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 10;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 8;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 6;
        break;
    case 5: // TH 40%
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 8;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 6;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 4;
        break;
#if UPGRADE_M8
    case 6:
        sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS] = 6;
        sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS] = 4;
        sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS] = 2;
        break;
#endif
    default:
        printf("Error - Invalid sb_class_level");
        break;
    }
}
#endif
#endif
#endif

#if MULTI_BAND_ACTIONS
#if NON_UNIFORM_NSQ_BANDING
uint8_t m0_nsq_cycles_reduction_th[19] = {
 0, // NONE
 17, //[85%;100%]
 15,//[75%;85%]
 14,//[65%;75%]
 13,//[60%;65%]
 12,//[55%;60%]
 11,//[50%;65%]
 10,//[45%;50%]
 9,//[40%;45%]
 8,//[35%;40%]
 7,//[30%;35%]
 6,//[25%;30%]
 6,//[20%;25%]
 5,//[17%;20%]
 5,//[14%;17%]
 4,//[10%;14%]
 3,//[6%;10%]
 2,//[3%;6%]
 1 //[0%;3%]
};
uint8_t m1_nsq_cycles_reduction_th[19] = {
 0, // NONE
 17, //[85%;100%]
 15,//[75%;85%]
 14,//[65%;75%]
 13,//[60%;65%]
 12,//[55%;60%]
 11,//[50%;65%]
 10,//[45%;50%]
 9,//[40%;45%]
 8,//[35%;40%]
 7,//[30%;35%]
 6,//[25%;30%]
 6,//[20%;25%]
 5,//[17%;20%]
 5,//[14%;17%]
 4,//[10%;14%]
 3,//[6%;10%]
 2,//[3%;6%]
 1 //[0%;3%]
};
#else
uint8_t m0_nsq_cycles_reduction_th[21] = {
0, // NONE
50,//[95%;100%]
 35,//[90%;95%]
 25,//[85%;90%]
 20,//[80%;85%]
 17,//[75%;80%]
 15,//[70%;75%]
 14,//[65%;70%]
 13,//[60%;65%]
 12,//[55%;60%]
 11,//[50%;55%]
 10,//[45%;50%]
 9,//[40%;45%]
 8,//[35%;40%]
 7,//[30%;35%]
 6,//[25%;30%]
 5,//[20%;25%]
 4,//[15%;20%]
 3,//[10%;15%]
 2,//[5%;10%]
1 //[0%;5%]
};
uint8_t m1_nsq_cycles_reduction_th[21] = {
0, // NONE
50,//[95%;100%]
 50,//[90%;95%]
 35,//[85%;90%]
 25,//[80%;85%]
 20,//[75%;80%]
 17,//[70%;75%]
 15,//[65%;70%]
 14,//[60%;65%]
 13,//[55%;60%]
 12,//[50%;55%]
 11,//[45%;50%]
 10,//[40%;45%]
 9,//[35%;40%]
 8,//[30%;35%]
 7,//[25%;30%]
 6,//[20%;25%]
 5,//[15%;20%]
 4,//[10%;15%]
 3,//[5%;10%]
2 //[0%;5%]
};
#endif
#endif
#if NSQ_CYCLES_REDUCTION
#if ADAPTIVE_NSQ_CR
void set_nsq_cycle_redcution_controls(ModeDecisionContext *mdctxt, uint16_t nsq_cycles_red_mode) {
    NsqCycleRControls*nsq_cycle_red_ctrls = &mdctxt->nsq_cycles_red_ctrls;

    switch (nsq_cycles_red_mode)
    {
    case 0: // nsq_cycles_reduction Off
        nsq_cycle_red_ctrls->enabled = 0;
        nsq_cycle_red_ctrls->th = 0;
        break;
    case 1:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 1;
        break;
    case 2:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 2;
        break;
    case 3:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 3;
        break;
    case 4:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 4;
        break;
    case 5:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 5;
        break;
    case 6:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 6;
        break;
    case 7:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 7;
        break;
    case 8:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 8;
        break;
    case 9:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 9;
        break;
    case 10:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 10;
        break;
    case 11:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 11;
        break;
    case 12:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 12;
        break;
    case 13:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 13;
        break;
#if JUNE26_ADOPTIONS
    case 14:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 15;
        break;
    case 15:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 50;
        break;
#else
    case 14:
        nsq_cycle_red_ctrls->enabled = 1;
#if NEW_NSQ_RED_LEVEL
        nsq_cycle_red_ctrls->th = 50;
#else
        nsq_cycle_red_ctrls->th = 14;
#endif
        break;
#endif
    default:
        assert(0);
        break;
    }
}
#else
void set_nsq_cycle_redcution_controls(ModeDecisionContext *mdctxt, uint8_t nsq_cycles_red_mode) {

    NsqCycleRControls*nsq_cycle_red_ctrls = &mdctxt->nsq_cycles_red_ctrls;

    switch (nsq_cycles_red_mode)
    {
    case 0: // nsq_cycles_reduction Off
        nsq_cycle_red_ctrls->enabled = 0;
        nsq_cycle_red_ctrls->th = 0;
        break;
    case 1:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 1;
        break;
    case 2:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 3;
        break;
    case 3:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 5;
        break;
    case 4:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 8;
        break;
    case 5:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 10;
        break;
    case 6:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 15;
        break;
    case 7:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 20;
        break;
    case 8:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 30;
        break;
    case 9:
        nsq_cycle_red_ctrls->enabled = 1;
        nsq_cycle_red_ctrls->th = 50;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#endif

#if SOFT_CYCLES_REDUCTION
void adaptive_md_cycles_redcution_controls(ModeDecisionContext *mdctxt, uint8_t adaptive_md_cycles_red_mode) {
    AMdCycleRControls*adaptive_md_cycles_red_ctrls = &mdctxt->admd_cycles_red_ctrls;
    switch (adaptive_md_cycles_red_mode)
    {
    case 0: // soft_cycles_reduction Off
        adaptive_md_cycles_red_ctrls->enabled = 0;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 0;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
    case 1:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 50;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
    case 2:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 150;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
    case 3:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 200;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
    case 4:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 300;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
#if TUNE_ADAPTIVE_MD_CR_TH
    case 5:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 700;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
     case 6:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 1000;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
    case 7:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 1500;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
    case 8:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 2000;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
     case 9:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 5000;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
#else
    case 5:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 1000;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
    case 6:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 1500;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
    case 7:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 2000;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
     case 8:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->sq_weight_th = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 5000;
        adaptive_md_cycles_red_ctrls->nics_th = 0;
        adaptive_md_cycles_red_ctrls->mrp_th = 0;
        adaptive_md_cycles_red_ctrls->compound_th = 0;
        break;
#endif
    default:
        assert(0);
        break;
    }
}
#endif
#if DEPTH_CYCLES_REDUCTION
void set_depth_cycle_redcution_controls(ModeDecisionContext *mdctxt, uint8_t depth_cycles_red_mode) {

    DepthCycleRControls*depth_cycle_red_ctrls = &mdctxt->depth_cycles_red_ctrls;

    switch (depth_cycles_red_mode)
    {
    case 0: // depth_cycles_reduction Off
        depth_cycle_red_ctrls->enabled = 0;
        depth_cycle_red_ctrls->th = 0;
        break;
    case 1:
        depth_cycle_red_ctrls->enabled = 1;
        depth_cycle_red_ctrls->th = 1;
        break;
    case 2:
        depth_cycle_red_ctrls->enabled = 1;
        depth_cycle_red_ctrls->th = 5;
        break;
    case 3:
        depth_cycle_red_ctrls->enabled = 1;
        depth_cycle_red_ctrls->th = 10;
        break;
    case 4:
        depth_cycle_red_ctrls->enabled = 1;
        depth_cycle_red_ctrls->th = 15;
        break;
    case 5:
        depth_cycle_red_ctrls->enabled = 1;
        depth_cycle_red_ctrls->th = 20;
        break;
     case 6:
        depth_cycle_red_ctrls->enabled = 1;
        depth_cycle_red_ctrls->th = 25;
        break;
     case 7:
        depth_cycle_red_ctrls->enabled = 1;
        depth_cycle_red_ctrls->th = 30;
        break;
     case 8:
        depth_cycle_red_ctrls->enabled = 1;
        depth_cycle_red_ctrls->th = 40;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if COEFF_BASED_TXT_BYPASS
void set_txt_cycle_reduction_controls(ModeDecisionContext *mdctxt, uint8_t txt_cycles_red_mode) {

    TxtCycleRControls* txt_cycle_red_ctrls = &mdctxt->txt_cycles_red_ctrls;

    switch (txt_cycles_red_mode)
    {
    case 0: // txt_cycles_reduction Off
        txt_cycle_red_ctrls->enabled = 0;
        txt_cycle_red_ctrls->intra_th = 0;
        txt_cycle_red_ctrls->inter_th = 0;
        break;
    case 1:
        txt_cycle_red_ctrls->enabled = 1;
        txt_cycle_red_ctrls->intra_th = 0;
        txt_cycle_red_ctrls->inter_th = 1;
        break;
    case 2:
        txt_cycle_red_ctrls->enabled = 1;
        txt_cycle_red_ctrls->intra_th = 1;
        txt_cycle_red_ctrls->inter_th = 3;
        break;
    case 3:
        txt_cycle_red_ctrls->enabled = 1;
        txt_cycle_red_ctrls->intra_th = 1;
        txt_cycle_red_ctrls->inter_th = 5;
        break;
    case 4:
        txt_cycle_red_ctrls->enabled = 1;
        txt_cycle_red_ctrls->intra_th = 3;
        txt_cycle_red_ctrls->inter_th = 7;
        break;
    case 5:
        txt_cycle_red_ctrls->enabled = 1;
        txt_cycle_red_ctrls->intra_th = 5;
        txt_cycle_red_ctrls->inter_th = 8;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if COEFF_BASED_TXS_BYPASS
void set_txs_cycle_reduction_controls(ModeDecisionContext *mdctxt, uint8_t txs_cycles_red_mode) {

    TxsCycleRControls* txs_cycle_red_ctrls = &mdctxt->txs_cycles_red_ctrls;

    switch (txs_cycles_red_mode)
    {
    case 0: // txt_cycles_reduction Off
        txs_cycle_red_ctrls->enabled = 0;
        txs_cycle_red_ctrls->intra_th = 0;
        txs_cycle_red_ctrls->inter_th = 0;
        break;
    case 1:
        txs_cycle_red_ctrls->enabled = 1;
        txs_cycle_red_ctrls->intra_th = 0;
        txs_cycle_red_ctrls->inter_th = 30;
        break;
    case 2:
        txs_cycle_red_ctrls->enabled = 1;
        txs_cycle_red_ctrls->intra_th = 15;
        txs_cycle_red_ctrls->inter_th = 50;
        break;
    case 3:
        txs_cycle_red_ctrls->enabled = 1;
        txs_cycle_red_ctrls->intra_th = 25;
        txs_cycle_red_ctrls->inter_th = 75;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
/******************************************************
* Derive EncDec Settings for OQ
Input   : encoder mode and pd pass
Output  : EncDec Kernel signal(s)
******************************************************/
#if MD_CONFIG_SB
EbErrorType signal_derivation_enc_dec_kernel_oq(
    SequenceControlSet *sequence_control_set_ptr,
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
#if REMOVE_MR_MACRO
    EbEncMode enc_mode = pcs_ptr->enc_mode;
#else
    uint8_t enc_mode = pcs_ptr->enc_mode;
#endif
    uint8_t pd_pass = context_ptr->pd_pass;

#if ON_OFF_FEATURE_MRP
    // mrp level
    context_ptr->mrp_level = pcs_ptr->parent_pcs_ptr->mrp_level ;
#endif

#if OPT_BLOCK_INDICES_GEN_2
#if SB_CLASSIFIER
    // sb_classifier levels
    // Level                Settings
    // 0                    Off
    // 1                    TH 80%
    // 2                    TH 70%
    // 3                    TH 60%
    // 4                    TH 50%
    // 5                    TH 40%
    if (pd_pass == PD_PASS_0)
        context_ptr->enable_area_based_cycles_allocation = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->enable_area_based_cycles_allocation = 0;
    else {
        if (pcs_ptr->slice_type == I_SLICE)
            context_ptr->enable_area_based_cycles_allocation = 0;
#if COEFF_BASED_BYPASS_OFF_480P
        // Do not use cycles reduction algorithms in 480p and below
        else if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
            context_ptr->enable_area_based_cycles_allocation = 0;
        else
            context_ptr->enable_area_based_cycles_allocation = 1;
#else
        else
            context_ptr->enable_area_based_cycles_allocation = 1;
#endif
    }
 #if MULTI_BAND_ACTIONS
#if !CLEANUP_CYCLE_ALLOCATION
    context_ptr->coeffcients_area_based_cycles_allocation_level = 1;
#endif
#else
    if (MR_MODE) {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 1;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 0;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 0;
    }
    else if (enc_mode == ENC_M0) {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 3;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
#if MAY17_ADOPTIONS
            context_ptr->coeffcients_area_based_cycles_allocation_level =
            pcs_ptr->parent_pcs_ptr->sc_content_detected ? 3 : 2;
#else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
#endif
#if MAY12_ADOPTIONS
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_360p_RANGE)
#else
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
#endif
            context_ptr->coeffcients_area_based_cycles_allocation_level = 1;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 0;
    }
#if !NEW_M1_CAND
    else if (enc_mode == ENC_M1) {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 4;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 3;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 1;
    }
#endif
#if MAY03_4K_10BIT_ADOPTS
    else if (enc_mode <= ENC_M1) {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 4;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 4;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 3;
#if MAY12_ADOPTIONS
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_480p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_360p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 1;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 0;
#else
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
#endif
    }
#endif
#if APR24_ADOPTIONS_M6_M7
    else if (enc_mode <= ENC_M6) {
#else
#if UPGRADE_M8
    else if (enc_mode <= ENC_M7) {
#else
    else {
#endif
#endif
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 5;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 4;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 3;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
    }
#if UPGRADE_M8
    else {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 6;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 6;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 5;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 4;
    }
#endif
#endif
#endif
#endif
    // Tx_search Level                                Settings
    // 0                                              OFF
    // 1                                              Tx search at encdec
    // 2                                              Tx search at inter-depth
    // 3                                              Tx search at full loop
    if (pd_pass == PD_PASS_0)
        context_ptr->tx_search_level = TX_SEARCH_OFF;
    else if (pd_pass == PD_PASS_1)
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
#if REMOVE_UNUSED_CODE_PH2
    else
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
#else
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR2_M8_ADOPTIONS
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
#else
        if (enc_mode <= ENC_M6)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
#endif
#if MAR25_ADOPTIONS
    else if (enc_mode <= ENC_M8)

#else
#if MAR17_ADOPTIONS
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M4)
#endif
#endif
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;

#if MAR2_M8_ADOPTIONS
    else {
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    }
#else
    else if (enc_mode <= ENC_M7) {
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    }
    else
        context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
#endif
#endif
#if TXT_CONTROL
    // Set MD tx_level
    // md_txt_search_level                            Settings
    // 0                                              FULL
    // 1                                              Tx_weight 1
    // 2                                              Tx_weight 2
    // 3                                              Tx_weight 1 + disabling rdoq and sssse
    // 4                                              Tx_weight 1 + disabling rdoq and sssse + reduced set
    if (pd_pass == PD_PASS_0)
        context_ptr->md_txt_search_level = 1;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_txt_search_level = 1;
#if !MAY19_ADOPTIONS
    else if (MR_MODE)
        context_ptr->md_txt_search_level = 0;
#endif
    else {
#if UNIFY_SC_NSC
#if NEW_M8
        context_ptr->md_txt_search_level = 0;
#else
        if (enc_mode <= ENC_M7)
            context_ptr->md_txt_search_level = 0;
        else
            context_ptr->md_txt_search_level = 3;
#endif
#else
#if APR23_ADOPTIONS_2
        // New adoption levels after UPDATE_TXT_LEVEL
        if (enc_mode <= ENC_M0)
#if UNIFY_SC_NSC
            context_ptr->md_txt_search_level = pcs_ptr->parent_pcs_ptr->sc_content_detected ? 2 : 1;
#else
#if MAY17_ADOPTIONS
            context_ptr->md_txt_search_level = pcs_ptr->parent_pcs_ptr->sc_content_detected ? 2 : 1;
#else
            context_ptr->md_txt_search_level = 1;
#endif
#endif
#if JUNE17_ADOPTIONS
        else if (enc_mode <= ENC_M7)
#else
#if JUNE11_ADOPTIONS
        else if (enc_mode <= ENC_M6)
#else
#if PRESET_SHIFITNG
        else if (enc_mode <= ENC_M5)
#else
        else if (enc_mode <= ENC_M7)
#endif
#endif
#endif
            context_ptr->md_txt_search_level = 2;
        else if (enc_mode <= ENC_M8)
            context_ptr->md_txt_search_level = 3;
#else
#if NEW_M1_CAND
        if(pcs_ptr->parent_pcs_ptr->sc_content_detected)
            if (enc_mode <= ENC_M0)
                context_ptr->md_txt_search_level = 1;
            else
                context_ptr->md_txt_search_level = 2;
#if APR23_ADOPTIONS_2
        else if (enc_mode <= ENC_M2)
#else
        else if (enc_mode <= ENC_M3)
#endif
            context_ptr->md_txt_search_level = 1;
        else
            context_ptr->md_txt_search_level = 2;
#else
#if APR23_ADOPTIONS
        if (enc_mode <= ENC_M3)
#else
        if (enc_mode <= ENC_M0)
#endif
            context_ptr->md_txt_search_level = 1;
#if APR23_ADOPTIONS
        else
            context_ptr->md_txt_search_level = 2;
#else
        else if (enc_mode <= ENC_M2)
            context_ptr->md_txt_search_level = 2;
        else if (enc_mode <= ENC_M8)
            context_ptr->md_txt_search_level = 3;
        else
            context_ptr->md_txt_search_level = 4;
#endif
#endif
#endif
#endif
    }
#else
    // Set tx search skip weights (MAX_MODE_COST: no skipping; 0: always skipping)
    if (pd_pass == PD_PASS_0)
        context_ptr->tx_weight = MAX_MODE_COST;
    else if (pd_pass == PD_PASS_1)
        context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
    else if (MR_MODE) // tx weight
        context_ptr->tx_weight = MAX_MODE_COST;
    else {
        if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
            context_ptr->tx_weight = MAX_MODE_COST;
#if MAR12_ADOPTIONS
        else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR18_ADOPTIONS
            context_ptr->tx_weight = MAX_MODE_COST;
#else
            if (enc_mode <= ENC_M3)
                context_ptr->tx_weight = MAX_MODE_COST;
            else
                context_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
#endif
#endif
#if APR02_ADOPTIONS
        else if (enc_mode <= ENC_M4)
#else
#if MAR10_ADOPTIONS
        else if (enc_mode <= ENC_M1)
#else
        else if (pcs_ptr->parent_pcs_ptr->sc_content_detected && enc_mode <= ENC_M0)
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else if (enc_mode <= ENC_M0)
#endif
#endif
            context_ptr->tx_weight = MAX_MODE_COST;
#if MAR30_ADOPTIONS
#if !APR02_ADOPTIONS
        else if (enc_mode <= ENC_M4)
            context_ptr->tx_weight = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? MAX_MODE_COST : FC_SKIP_TX_SR_TH025;
#endif
        else
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
#else
#if MAR3_M2_ADOPTIONS
#if MAR4_M3_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (enc_mode <= ENC_M8 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
        else if (enc_mode <= ENC_M3 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#else
        else if (enc_mode <= ENC_M2 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#else
        else if (enc_mode <= ENC_M1 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
#endif
    }

    // Set tx search reduced set falg (0: full tx set; 1: reduced tx set; 1: two
    // tx))
    if (pd_pass == PD_PASS_0)
        context_ptr->tx_search_reduced_set = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->tx_search_reduced_set = 0;
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR17_ADOPTIONS
        if (enc_mode <= ENC_M8)
            context_ptr->tx_search_reduced_set = 0;
#else
        if (enc_mode <= ENC_M5)
            context_ptr->tx_search_reduced_set = 0;
        else if (enc_mode <= ENC_M6)
            if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
                context_ptr->tx_search_reduced_set = 0;
            else
                context_ptr->tx_search_reduced_set = 1;
#endif
#if MAR2_M8_ADOPTIONS
        else
            context_ptr->tx_search_reduced_set = 1;
#else
        else if (enc_mode <= ENC_M7)
            context_ptr->tx_search_reduced_set = 1;
        else
            context_ptr->tx_search_reduced_set = 2;
#endif
    else if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
        context_ptr->tx_search_reduced_set = 0;
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
    else if (enc_mode <= ENC_M8)
#else
    else if (enc_mode <= ENC_M5)
#endif
#else
    else if (enc_mode <= ENC_M3)
#endif
        context_ptr->tx_search_reduced_set = 0;
    else
        context_ptr->tx_search_reduced_set = 1;
#endif
#if COEFF_BASED_TXT_BYPASS
    uint8_t txt_cycles_reduction_level = 0;
#if SEPARATE_ADAPTIVE_TXT_INTER_INTRA
    if (pcs_ptr->parent_pcs_ptr->slice_type == I_SLICE) {
        txt_cycles_reduction_level = 0;
    }
    else {
        if (pd_pass == PD_PASS_0)
            txt_cycles_reduction_level = 0;
        else if (pd_pass == PD_PASS_1)
            txt_cycles_reduction_level = 0;
#if UNIFY_SC_NSC
#if JUNE26_ADOPTIONS
        else if (enc_mode <= ENC_M5)
#else
        else if (enc_mode <= ENC_M4)
#endif
            txt_cycles_reduction_level = 0;
        else
            txt_cycles_reduction_level = 5;
#else
        else if (enc_mode <= ENC_M2)
            txt_cycles_reduction_level = 0;
#if JUNE17_ADOPTIONS
        else if (enc_mode <= ENC_M3)
            txt_cycles_reduction_level = 1;
        else
            txt_cycles_reduction_level = 5;
#else
        else
            txt_cycles_reduction_level = 1;
#endif
#endif
    }
#endif
    set_txt_cycle_reduction_controls(context_ptr, txt_cycles_reduction_level);
#endif
    // Interpolation search Level                     Settings
    // 0                                              OFF
    // 1                                              Interpolation search at
    // inter-depth 2                                              Interpolation
    // search at full loop 3                                              Chroma
    // blind interpolation search at fast loop 4 Interpolation search at fast loop
    if (pd_pass == PD_PASS_0)
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
    else if (pd_pass == PD_PASS_1) {
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
    }
#if !UNIFY_SC_NSC
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#endif
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
    else if (enc_mode <= ENC_M8)
#else
    else if (enc_mode <= ENC_M5)
#endif
#else
    else if (enc_mode <= ENC_M3)
#endif
        context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
#if !MAR10_ADOPTIONS
    else if (enc_mode <= ENC_M7)
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
        else
            context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#endif
    else
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;

    // Set Chroma Mode
    // Level                Settings
    // CHROMA_MODE_0  0     Full chroma search @ MD
    // CHROMA_MODE_1  1     Fast chroma search @ MD
    // CHROMA_MODE_2  2     Chroma blind @ MD + CFL @ EP
    // CHROMA_MODE_3  3     Chroma blind @ MD + no CFL @ EP
    if (pd_pass == PD_PASS_0)
        context_ptr->chroma_level = CHROMA_MODE_2; // or CHROMA_MODE_3
    else if (pd_pass == PD_PASS_1) {
        context_ptr->chroma_level = CHROMA_MODE_1;
    }
    else if (sequence_control_set_ptr->static_config.set_chroma_mode ==
        DEFAULT) {
#if UNIFY_SC_NSC
        if (enc_mode <= ENC_M5)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else
            context_ptr->chroma_level = CHROMA_MODE_1;
#else
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR2_M7_ADOPTIONS
#if MAR10_ADOPTIONS
#if FIX_CHROMA_PALETTE_INTERACTION
#if MAY12_ADOPTIONS
#if JUNE17_ADOPTIONS
            if (enc_mode <= ENC_M6)
#else
#if JUNE11_ADOPTIONS
            if (enc_mode <= ENC_M3)
#else
#if PRESET_SHIFITNG
            if (enc_mode <= ENC_M2)
#else
            if (enc_mode <= ENC_M4)
#endif
#endif
#endif
#else
            if (enc_mode <= ENC_M0)
#endif
                context_ptr->chroma_level = CHROMA_MODE_0;
#if REMOVE_UNUSED_CODE_PH2
            else
#else
            else if (enc_mode <= ENC_M8)
#endif
#else
            if (enc_mode <= ENC_M8)
#endif
#else
            if (enc_mode <= ENC_M7)
#endif
#else
            if (enc_mode <= ENC_M6)
#endif
                context_ptr->chroma_level = CHROMA_MODE_1;
#if !REMOVE_UNUSED_CODE_PH2
            else if (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                context_ptr->chroma_level = CHROMA_MODE_1;
            else
                context_ptr->chroma_level =
                (sequence_control_set_ptr->encoder_bit_depth == EB_8BIT)
                ? CHROMA_MODE_2
                : CHROMA_MODE_3;
#endif
#if PRESETS_SHIFT
#if JUNE17_ADOPTIONS
        else if (enc_mode <= ENC_M4)
#else
#if JUNE11_ADOPTIONS
        else if (enc_mode <= ENC_M3)
#else
#if PRESET_SHIFITNG
        else if (enc_mode <= ENC_M2)
#else
        else if (enc_mode <= ENC_M4)
#endif
#endif
#endif
            context_ptr->chroma_level = CHROMA_MODE_0;
#else
#if MAR17_ADOPTIONS
        else if (enc_mode <= ENC_M7)
            context_ptr->chroma_level = CHROMA_MODE_0;
#else
#if MAR10_ADOPTIONS
        else if (enc_mode <= ENC_M4)
#else
        else if (enc_mode <= ENC_M3)
#endif
            context_ptr->chroma_level = CHROMA_MODE_0;
        else if (enc_mode <= ENC_M5 &&
            pcs_ptr->temporal_layer_index == 0)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else if (enc_mode <= ENC_M5)
            context_ptr->chroma_level = CHROMA_MODE_1;
#endif
#endif
        else
#if MAY03_4K_10BIT_ADOPTS
            context_ptr->chroma_level = CHROMA_MODE_1;
#else
            context_ptr->chroma_level =
            (sequence_control_set_ptr->encoder_bit_depth == EB_8BIT)
#if ADOPT_CHROMA_MODE1_CFL_OFF
            ? CHROMA_MODE_1
#else
            ? CHROMA_MODE_2
#endif
            : CHROMA_MODE_3;
#endif
#endif
    }
    else // use specified level
        context_ptr->chroma_level =
        sequence_control_set_ptr->static_config.set_chroma_mode;

    // Chroma independent modes search
    // Level                Settings
    // 0                    post first md_stage
    // 1                    post last md_stage
#if APR22_ADOPTIONS
#if FIXED_LAST_STAGE_SC
#if !UNIFY_SC_NSC
#if JUNE8_ADOPTIONS
    if (MR_MODE) {
        context_ptr->chroma_at_last_md_stage = 0;
        context_ptr->chroma_at_last_md_stage_intra_th = (uint64_t)~0;
        context_ptr->chroma_at_last_md_stage_cfl_th = (uint64_t)~0;
    }
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
#else
#if MAY19_ADOPTIONS
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
#else
    if (MR_MODE) {
        context_ptr->chroma_at_last_md_stage = 0;
    }
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
#endif
#endif
        if (enc_mode <= ENC_M0) {
            context_ptr->chroma_at_last_md_stage = 0;
            context_ptr->chroma_at_last_md_stage_intra_th = (uint64_t) ~0;
            context_ptr->chroma_at_last_md_stage_cfl_th = (uint64_t)~0;
        }
#if MAY12_ADOPTIONS
#if JUNE17_ADOPTIONS
        else if (enc_mode <= ENC_M3) {
#else
#if JUNE8_ADOPTIONS
        else if (enc_mode <= ENC_M2) {
#else
#if PRESET_SHIFITNG
        else if(enc_mode <= ENC_M1) {
#else
        else if(enc_mode <= ENC_M2) {
#endif
#endif
#endif
            context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
            context_ptr->chroma_at_last_md_stage_intra_th = 150;
            context_ptr->chroma_at_last_md_stage_cfl_th = 150;
        }
        else {
            context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
            context_ptr->chroma_at_last_md_stage_intra_th = 130;
            context_ptr->chroma_at_last_md_stage_cfl_th = 130;
    }
#else
        else {
            context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
            context_ptr->chroma_at_last_md_stage_intra_th = (uint64_t)~0;
            context_ptr->chroma_at_last_md_stage_cfl_th = (uint64_t)~0;
        }
#endif
    }
#else
#if REMOVE_MR_MACRO
    if (enc_mode <= ENC_MR) {
#else
    if (MR_MODE) {
#endif
        context_ptr->chroma_at_last_md_stage = 0;
        context_ptr->chroma_at_last_md_stage_intra_th = (uint64_t)~0;
        context_ptr->chroma_at_last_md_stage_cfl_th = (uint64_t)~0;
    }
#endif
#if MAY12_ADOPTIONS
#if JUNE23_ADOPTIONS
    else if (enc_mode <= ENC_M3) {
#else
#if JUNE8_ADOPTIONS
    else if (enc_mode <= ENC_M2) {
#else
#if PRESET_SHIFITNG
    else if (enc_mode <= ENC_M1) {
#else
    else if (enc_mode <= ENC_M2) {
#endif
#endif
#endif
        context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
        context_ptr->chroma_at_last_md_stage_intra_th = 130;
        context_ptr->chroma_at_last_md_stage_cfl_th = 130;
    }
    else {
        context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
        context_ptr->chroma_at_last_md_stage_intra_th = 100;
        context_ptr->chroma_at_last_md_stage_cfl_th = 100;
    }
#else
    else {
        context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
        context_ptr->chroma_at_last_md_stage_intra_th = 130;
        context_ptr->chroma_at_last_md_stage_cfl_th = 130;
    }
#endif
#else
    context_ptr->chroma_at_last_md_stage =
        MR_MODE ? 0 : (context_ptr->chroma_level == CHROMA_MODE_0 && !pcs_ptr->parent_pcs_ptr->sc_content_detected) ? 1 : 0;
#endif
#else
#if MAR30_ADOPTIONS
    context_ptr->chroma_at_last_md_stage =
        MR_MODE ? 0 : (context_ptr->chroma_level == CHROMA_MODE_0 && !pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 1 : 0;
#else
    context_ptr->chroma_at_last_md_stage =
        MR_MODE ? 0 : (context_ptr->chroma_level == CHROMA_MODE_0 && !pcs_ptr->parent_pcs_ptr->sc_content_detected) ? 1 : 0;
#endif
#endif
#if M5_CHROMA_NICS
    // Chroma independent modes nics
    // Level                Settings
    // 0                    All supported modes.
    // 1                    All supported modes in  Intra picture and 4 in inter picture
#if MAY12_ADOPTIONS
    context_ptr->independent_chroma_nics = 0;
#else
#if SHIFT_M5_SC_TO_M3
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        if (enc_mode <= ENC_M2)
            context_ptr->independent_chroma_nics = 0;
        else
            context_ptr->independent_chroma_nics = 1;
    else if (enc_mode <= ENC_M3)
        context_ptr->independent_chroma_nics = 0;
    else
        context_ptr->independent_chroma_nics = 1;
#else
#if PRESETS_SHIFT
    context_ptr->independent_chroma_nics = enc_mode <= ENC_M3 ? 0 : 1;
#else
#if MAR23_ADOPTIONS
    context_ptr->independent_chroma_nics = enc_mode <= ENC_M4 ? 0 : 1;
#else
    context_ptr->independent_chroma_nics = enc_mode == ENC_M5 ? 1 : 0;
#endif
#endif
#endif
#endif
#endif
#if ADDED_CFL_OFF
    // Cfl level
    // Level                Settings
    // 0                    Allow cfl
    // 1                    Disable cfl
#if PRESETS_SHIFT
#if ALLOW_CFL_M8
    context_ptr->md_disable_cfl = EB_FALSE;
#else
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->md_disable_cfl = EB_FALSE;
#if APR23_ADOPTIONS_2
    else if (enc_mode <= ENC_M5)
#else
    else if (enc_mode <= ENC_M4)
#endif
        context_ptr->md_disable_cfl = EB_FALSE;
    else
        context_ptr->md_disable_cfl = EB_TRUE;
#endif
#else
#if ADOPT_CHROMA_MODE1_CFL_OFF
#if MAR17_ADOPTIONS
    if (enc_mode <= ENC_M7)
#else
    if(enc_mode <= ENC_M5)
#endif
        context_ptr->md_disable_cfl = EB_FALSE;
    else if(!pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
        context_ptr->md_disable_cfl = EB_TRUE;
#endif
#if !REFACTOR_SIGNALS
    if (sequence_control_set_ptr->static_config.disable_cfl_flag == 1 && context_ptr->md_disable_cfl == EB_TRUE)
        context_ptr->chroma_at_last_md_stage = 0; // Indeprndent chroma search at last MD stage is not supported when CFL is off
#endif
#endif
#if CFL_REDUCED_ALPHA
    // libaom_short_cuts_ths
    // 1                    faster than libaom
    // 2                    libaom - default
#if MAY12_ADOPTIONS
    context_ptr->libaom_short_cuts_ths = 2;
#else
#if SHIFT_M5_SC_TO_M3
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        if (enc_mode <= ENC_M2)
            context_ptr->libaom_short_cuts_ths = 2;
        else
            context_ptr->libaom_short_cuts_ths = 1;
    else if (enc_mode <= ENC_M3)
        context_ptr->libaom_short_cuts_ths = 2;
    else
        context_ptr->libaom_short_cuts_ths = 1;
#else
#if PRESETS_SHIFT
    if (enc_mode <= ENC_M3)
        context_ptr->libaom_short_cuts_ths = 2;
    else
        context_ptr->libaom_short_cuts_ths = 1;
#else
#if MAR23_ADOPTIONS
    if (enc_mode <= ENC_M4)
        context_ptr->libaom_short_cuts_ths = 2;
    else
        context_ptr->libaom_short_cuts_ths = 1;
#else
    if (enc_mode == ENC_M5)
        context_ptr->libaom_short_cuts_ths = 1;
    else
        context_ptr->libaom_short_cuts_ths = 2;
#endif
#endif
#endif
#endif
#endif
#if UV_SEARCH_MODE_INJCECTION
    // 0                    inject all supprted chroma mode
    // 1                    follow the luma injection
#if MAY12_ADOPTIONS
    context_ptr->intra_chroma_search_follows_intra_luma_injection = 0;
#else
#if SHIFT_M5_SC_TO_M3
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        if (enc_mode <= ENC_M2)
            context_ptr->intra_chroma_search_follows_intra_luma_injection = 0;
        else
            context_ptr->intra_chroma_search_follows_intra_luma_injection = 1;
    else if (enc_mode <= ENC_M3)
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 0;
    else
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 1;
#else
#if PRESETS_SHIFT
    if (enc_mode <= ENC_M3)
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 0;
    else
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 1;
#else
#if MAR23_ADOPTIONS
    if (enc_mode <= ENC_M4)
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 0;
    else
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 1;
#else
    context_ptr->intra_chroma_search_follows_intra_luma_injection = 1;
     if (enc_mode == ENC_M5)
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 1;
    else
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 0;
#endif
#endif
#endif
#endif
#endif
#if M8_4x4
#if APR24_M3_ADOPTIONS
     // Set disallow_4x4
#if !UNIFY_SC_NSC
     if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if APR25_12AM_ADOPTIONS
#if JUNE17_ADOPTIONS
         if (enc_mode <= ENC_M4)
#else
#if PRESET_SHIFITNG
         if (enc_mode <= ENC_M5)
#else
         if (enc_mode <= ENC_M7)
#endif
#endif
#else
#if APR24_ADOPTIONS_M6_M7
         if (enc_mode <= ENC_M6)
#else
         if (enc_mode <= ENC_M5)
#endif
#endif
             context_ptr->disallow_4x4 = EB_FALSE;
#if REVERT_WHITE // disallow_4x4
#if JUNE17_ADOPTIONS
         else if (enc_mode <= ENC_M6)
#else
#if PRESET_SHIFITNG
         else if (enc_mode <= ENC_M5)
#else
         else if (enc_mode <= ENC_M7)
#endif
#endif
#else
         else if (enc_mode <= ENC_M8)
#endif
            context_ptr->disallow_4x4 = (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
        else
            context_ptr->disallow_4x4 = EB_TRUE;
#if JUNE11_ADOPTIONS
     else if (enc_mode <= ENC_M2)
#else
#if M1_C2_ADOPTIONS
     else if (enc_mode <= ENC_M0)
#else
     else if (enc_mode <= ENC_M2)
#endif
#endif
#else
#if JUNE26_ADOPTIONS
     if (enc_mode <= ENC_M1)
#else
     if (enc_mode <= ENC_M2)
#endif
#endif
         context_ptr->disallow_4x4 = EB_FALSE;
#if JUNE25_ADOPTIONS
     else
         context_ptr->disallow_4x4 = (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
#else
#if UNIFY_SC_NSC
#if JUNE23_ADOPTIONS
     else if (enc_mode <= ENC_M4)
#else
     else if (enc_mode <= ENC_M3)
#endif
         context_ptr->disallow_4x4 = (pcs_ptr->temporal_layer_index == 0) ? EB_FALSE : EB_TRUE;
#endif
#if REVERT_WHITE // disallow_4x4
#if PRESET_SHIFITNG
     else if (enc_mode <= ENC_M5)
#else
     else if (enc_mode <= ENC_M7)
#endif
#else
     else if (enc_mode <= ENC_M8)
#endif
         context_ptr->disallow_4x4 = (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
     else
         context_ptr->disallow_4x4 = EB_TRUE;
#endif
#else
#if UPGRADE_M6_M7_M8
     if (enc_mode <= ENC_M5)
         context_ptr->disallow_4x4 = EB_FALSE;
#if M5_I_4x4
     else if (enc_mode <= ENC_M8)
#else
     else if (enc_mode <= ENC_M7)
#endif
         context_ptr->disallow_4x4 = (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
     else
         context_ptr->disallow_4x4 = EB_TRUE;
#else
     context_ptr->disallow_4x4 = pcs_ptr->enc_mode <= ENC_M5 ? EB_FALSE : EB_TRUE;
#endif
#endif
     // If SB non-multiple of 4, then disallow_4x4 could not be used
     // SB Stats
     uint32_t sb_width =
         MIN(sequence_control_set_ptr->sb_size_pix, pcs_ptr->parent_pcs_ptr->aligned_width - context_ptr->sb_ptr->origin_x);
     uint32_t sb_height =
         MIN(sequence_control_set_ptr->sb_size_pix, pcs_ptr->parent_pcs_ptr->aligned_height - context_ptr->sb_ptr->origin_y);
     if (sb_width % 8 != 0 || sb_height % 8 != 0) {
         context_ptr->disallow_4x4 = EB_FALSE;
     }
#endif
#if SB_CLASSIFIER
#if OPT_BLOCK_INDICES_GEN_2

#if PD0_PD1_NSQ_BLIND
#if MAY21_NSQ_OFF_FIX
#if JUNE23_ADOPTIONS
     if (pd_pass == PD_PASS_0)
#if REMOVE_MR_MACRO
         context_ptr->md_disallow_nsq = enc_mode <= ENC_MR ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
#else
         context_ptr->md_disallow_nsq = MR_MODE ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
#endif
     else if (pd_pass == PD_PASS_1)
#if REMOVE_MR_MACRO
         context_ptr->md_disallow_nsq = enc_mode <= ENC_MR ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
#else
         context_ptr->md_disallow_nsq = MR_MODE ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
#endif
#else
     if (pd_pass == PD_PASS_0)
         context_ptr->md_disallow_nsq = (enc_mode <= ENC_M0) ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
     else if (pd_pass == PD_PASS_1)
         context_ptr->md_disallow_nsq = (enc_mode <= ENC_M0) ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
#endif
#else
     if (pd_pass == PD_PASS_0)
         context_ptr->md_disallow_nsq = (enc_mode <= ENC_M0) ? 0 : 1;
     else if (pd_pass == PD_PASS_1)
         context_ptr->md_disallow_nsq = (enc_mode <= ENC_M0) ? 0 : 1;
#endif
     else
         // Update nsq settings based on the sb_class
#if NEW_CYCLE_ALLOCATION
#if DISALLOW_ALL_ACTIONS
         context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;
#else
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
            context_ptr->md_disallow_nsq = (context_ptr->enable_area_based_cycles_allocation &&  ((context_ptr->sb_class == HIGH_COMPLEX_CLASS) || (context_ptr->sb_class == MEDIUM_COMPLEX_CLASS))) ? 1 : pcs_ptr->parent_pcs_ptr->disallow_nsq;
        else
            context_ptr->md_disallow_nsq = (context_ptr->enable_area_based_cycles_allocation &&  context_ptr->sb_class == HIGH_COMPLEX_CLASS) ? 1 : pcs_ptr->parent_pcs_ptr->disallow_nsq;
#endif
#else
         context_ptr->md_disallow_nsq = (context_ptr->enable_area_based_cycles_allocation &&  context_ptr->sb_class == HIGH_COMPLEX_CLASS) ? 1 : pcs_ptr->parent_pcs_ptr->disallow_nsq;
#endif
#else
     // Update nsq settings based on the sb_class
     context_ptr->md_disallow_nsq = (context_ptr->enable_area_based_cycles_allocation &&  context_ptr->sb_class == HIGH_COMPLEX_CLASS ) ? 1 : pcs_ptr->parent_pcs_ptr->disallow_nsq;
#endif
#else
     context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;
#endif
#endif
#if !OPT_BLOCK_INDICES_GEN_2
#if SB_CLASSIFIER
    // sb_classifier levels
    // Level                Settings
    // 0                    Off
    // 1                    TH 80%
    // 2                    TH 70%
    // 3                    TH 60%
    // 4                    TH 50%
    // 5                    TH 40%
    if (pd_pass == PD_PASS_0)
        context_ptr->enable_area_based_cycles_allocation = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->enable_area_based_cycles_allocation = 0;
    else {
        if (pcs_ptr->slice_type == I_SLICE)
            context_ptr->enable_area_based_cycles_allocation = 0;
        else
            context_ptr->enable_area_based_cycles_allocation = 1;
    }

    if (MR_MODE) {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 1;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 0;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 0;
    }
    else if (enc_mode == ENC_M0) {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 3;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 1;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 0;
    }
    else if (enc_mode == ENC_M1) {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 4;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 3;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 1;
    }
    else {
        if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_4K_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 5;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 4;
        else if (pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_720p_RANGE)
            context_ptr->coeffcients_area_based_cycles_allocation_level = 3;
        else
            context_ptr->coeffcients_area_based_cycles_allocation_level = 2;
    }
#endif
#endif
#if !M8_CLEAN_UP
    // Set the full loop escape level
    // Level                Settings
    // 0                    Off
    // 1                    On but only INTRA
    // 2                    On both INTRA and INTER

    if (pd_pass == PD_PASS_0)
        context_ptr->full_loop_escape = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->full_loop_escape = 0;
    else

        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR17_ADOPTIONS
            context_ptr->full_loop_escape = 0;
#else
            if (enc_mode <= ENC_M5)
                context_ptr->full_loop_escape = 0;
            else
                context_ptr->full_loop_escape = 2;
#endif
#if MAR2_M7_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (enc_mode <= ENC_M8)
#else
        else if (enc_mode <= ENC_M7)
#endif
#else
        else if (enc_mode <= ENC_M5)
#endif
            context_ptr->full_loop_escape = 0;
        else
            context_ptr->full_loop_escape = 2;
#endif


    // Set global MV injection
    // Level                Settings
    // 0                    Injection off
    // 1                    On
    if (sequence_control_set_ptr->static_config.enable_global_motion == EB_TRUE) {
        if (pd_pass == PD_PASS_0)
            context_ptr->global_mv_injection = 0;
        else if (pd_pass == PD_PASS_1)
            context_ptr->global_mv_injection = 0;
        else
#if UNIFY_SC_NSC
#if JUNE26_ADOPTIONS
            if (enc_mode <= ENC_M6)
#else
            if (enc_mode <= ENC_M5)
#endif
                context_ptr->global_mv_injection = 1;
            else
                context_ptr->global_mv_injection = 0;
#else
#if MAR4_M6_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAY19_ADOPTIONS
#if JUNE17_ADOPTIONS
                if (enc_mode <= ENC_M5)
#else
#if PRESET_SHIFITNG
                if (enc_mode <= ENC_M4)
#else
                if (enc_mode <= ENC_M6)
#endif
#endif
#else
#if MAY12_ADOPTIONS
                if (enc_mode <= ENC_M4)
#else
#if SHIFT_M5_SC_TO_M3
                if (enc_mode <= ENC_M2)
#else
#if PRESETS_SHIFT
                if (enc_mode <= ENC_M4)
#else
#if MAR12_ADOPTIONS
#if MAR17_ADOPTIONS
                if (enc_mode <= ENC_M7)
#else
                if (enc_mode <= ENC_M3)
#endif
#else
#if MAR10_ADOPTIONS
#if MAR11_ADOPTIONS
                if (enc_mode <= ENC_M1)
#else
                if (enc_mode <= ENC_M2)
#endif
#else
                if (enc_mode <= ENC_M3)
#endif
#endif
#endif
#endif
#endif
#endif
                    context_ptr->global_mv_injection = 1;
                else
                    context_ptr->global_mv_injection = 0;
#if MAY19_ADOPTIONS
#if JUNE17_ADOPTIONS
        else if (enc_mode <= ENC_M5)
#else
#if PRESET_SHIFITNG
        else if (enc_mode <= ENC_M4)
#else
        else if (enc_mode <= ENC_M6)
#endif
#endif
#else
#if PRESETS_SHIFT
            else if (enc_mode <= ENC_M4)
#else
#if MAR17_ADOPTIONS
            else if (enc_mode <= ENC_M7)
#else
            else if (enc_mode <= ENC_M5)
#endif
#endif
#endif
#else
            if (enc_mode <= ENC_M3)
#endif
                context_ptr->global_mv_injection = 1;
            else
                context_ptr->global_mv_injection = 0;
#endif
    }
    else
        context_ptr->global_mv_injection = 0;

    if (pd_pass == PD_PASS_0)
        context_ptr->new_nearest_injection = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->new_nearest_injection = 0;
    else
        context_ptr->new_nearest_injection = 1;

    if (pd_pass == PD_PASS_0)
        context_ptr->new_nearest_near_comb_injection = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->new_nearest_near_comb_injection = 0;
    else

        if (sequence_control_set_ptr->static_config.new_nearest_comb_inject ==
            DEFAULT)
#if OPTIMIZE_NEAREST_NEW_NEAR
#if !UNIFY_SC_NSC
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if JUNE17_ADOPTIONS
                if (enc_mode <= ENC_M5)
#else
                if (enc_mode <= ENC_M4)
#endif
                    context_ptr->new_nearest_near_comb_injection = 1;
                else
                    context_ptr->new_nearest_near_comb_injection = 0;
            else
#endif
#if JUNE23_ADOPTIONS
                if (enc_mode <= ENC_M1)
#else
                if (enc_mode <= ENC_M0)
#endif
                    context_ptr->new_nearest_near_comb_injection = 1;
                else
                    context_ptr->new_nearest_near_comb_injection = 0;
#else
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR30_ADOPTIONS
#if MAY19_ADOPTIONS
#if PRESET_SHIFITNG
                if (enc_mode <= ENC_M4)
#else
                if (enc_mode <= ENC_M6)
#endif
#else
#if MAY12_ADOPTIONS
                if (enc_mode <= ENC_M4)
#else
#if SHIFT_M5_SC_TO_M3
                if (enc_mode <= ENC_M2)
#else
#if APR23_ADOPTIONS_2
                if (enc_mode <= ENC_M4)
#else
#if APR22_ADOPTIONS
                if (enc_mode <= ENC_M2)
#else
                if (enc_mode <= ENC_M0)
#endif
#endif
#endif
#endif
#endif
                    context_ptr->new_nearest_near_comb_injection = 1;
                else
                    context_ptr->new_nearest_near_comb_injection = 0;
#else
                context_ptr->new_nearest_near_comb_injection = 0;
#if MAR10_ADOPTIONS
            else if (enc_mode <= ENC_M1)
#else
            else if (enc_mode <= ENC_M0)
#endif
                context_ptr->new_nearest_near_comb_injection = 1;
#endif
#if APR22_ADOPTIONS
            else if(MR_MODE)
                context_ptr->new_nearest_near_comb_injection = 1;
#endif
            else
                context_ptr->new_nearest_near_comb_injection = 0;
#endif
        else
            context_ptr->new_nearest_near_comb_injection =
            sequence_control_set_ptr->static_config.new_nearest_comb_inject;

    // Set warped motion injection
    // Level                Settings
    // 0                    OFF
    // 1                    On

    if (pd_pass == PD_PASS_0) {
        context_ptr->warped_motion_injection = 0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->warped_motion_injection = 1;
    }
    else
        context_ptr->warped_motion_injection = 1;

    // Set unipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
    if (pd_pass == PD_PASS_0) {
        context_ptr->unipred3x3_injection = 0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->unipred3x3_injection = 2;
    }
    else
#if JUNE11_ADOPTIONS
    {
#if JUNE23_ADOPTIONS
        if (enc_mode <= ENC_M2)
#else
        if (enc_mode <= ENC_M3)
#endif
            context_ptr->unipred3x3_injection = 1;
        else
            context_ptr->unipred3x3_injection = 0;
    }
#else
#if SHIFT_M5_SC_TO_M3
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if JUNE8_ADOPTIONS
            if (enc_mode <= ENC_M2)
#else
#if PRESET_SHIFITNG
            if (enc_mode <= ENC_M1)
#else
            if (enc_mode <= ENC_M2)
#endif
#endif
                context_ptr->unipred3x3_injection = 1;
            else
                context_ptr->unipred3x3_injection = 0;
        else
#endif
#if PRESETS_SHIFT
#if PRESET_SHIFITNG
        if (enc_mode <= ENC_M2)
#else
        if (enc_mode <= ENC_M4)
#endif
#else
        if (enc_mode <= ENC_M7)
#endif
            context_ptr->unipred3x3_injection = 1;

        else
            context_ptr->unipred3x3_injection = 0;
#endif

    // Set bipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
    if (pd_pass == PD_PASS_0) {
        context_ptr->bipred3x3_injection = 0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->bipred3x3_injection = 2;
    }
    else if (sequence_control_set_ptr->static_config.bipred_3x3_inject ==
        DEFAULT)
#if UNIFY_SC_NSC
#if JUNE23_ADOPTIONS
        if (enc_mode <= ENC_M2)
#else
        if (enc_mode <= ENC_M1)
#endif
            context_ptr->bipred3x3_injection = 1;
        else
            context_ptr->bipred3x3_injection = 2;
#else
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAY12_ADOPTIONS
            if (enc_mode <= ENC_M8)
#else
#if SHIFT_M5_SC_TO_M3
            if (enc_mode <= ENC_M2)
#else
#if PRESETS_SHIFT
            if (enc_mode <= ENC_M4)
#else
#if MAR30_ADOPTIONS
            if (enc_mode <= ENC_M7)
#else
#if MAR20_M4_ADOPTIONS
            if (enc_mode <= ENC_M3)
#else
            if (enc_mode <= ENC_M4)
#endif
#endif
#endif
#endif
#endif
                context_ptr->bipred3x3_injection = 1;
            else
#if MAR18_ADOPTIONS
#if M8_BIPRED_3x3 && !UPGRADE_M8
                if (enc_mode <= ENC_M5)
                    context_ptr->bipred3x3_injection = 2;
                else
                    context_ptr->bipred3x3_injection = 0;
#else
                context_ptr->bipred3x3_injection = 2;
#endif
#else
                context_ptr->bipred3x3_injection = 0;
#endif
#if MAR18_ADOPTIONS
#if PRESETS_SHIFT
#if JUNE11_ADOPTIONS
        else if (enc_mode <= ENC_M1)
#else
#if M1_COMBO_1 || NEW_M1_CAND
        else if (enc_mode <= ENC_M0)
#else
        else if (enc_mode <= ENC_M2)
#endif
#endif
#else
#if MAR20_M4_ADOPTIONS
        else if (enc_mode <= ENC_M3)
#else
        else if (enc_mode <= ENC_M4)
#endif
#endif
            context_ptr->bipred3x3_injection = 1;
        else
#if M8_BIPRED_3x3 && !UPGRADE_M8
    if (enc_mode <= ENC_M5)
        context_ptr->bipred3x3_injection = 2;
    else
        context_ptr->bipred3x3_injection = 0;
#else
            context_ptr->bipred3x3_injection = 2;
#endif
#else
#if MAR17_ADOPTIONS
        else if (enc_mode <= ENC_M7)
            context_ptr->bipred3x3_injection = 1;
#else
#if MAR12_ADOPTIONS
        else if (enc_mode <= ENC_M3)
#else
#if MAR3_M2_ADOPTIONS
#if MAR10_ADOPTIONS
#if MAR11_ADOPTIONS
        else if (enc_mode <= ENC_M2)
#else
        else if (enc_mode <= ENC_M3)
#endif
#else
        else if (enc_mode <= ENC_M2)
#endif
#else
        else if (enc_mode <= ENC_M1)
#endif
#endif
            context_ptr->bipred3x3_injection = 1;
#if MAR3_M6_ADOPTIONS
        else if (enc_mode <= ENC_M6)
#else
        else if (enc_mode <= ENC_M4)
#endif
            context_ptr->bipred3x3_injection = 2;
#endif
        else
            context_ptr->bipred3x3_injection = 0;
#endif
#endif
    else
        context_ptr->bipred3x3_injection =
        sequence_control_set_ptr->static_config.bipred_3x3_inject;

    // Level                Settings
    // 0                    Level 0: OFF
    // 1                    Level 1: sub-pel refinement off
    // 2                    Level 2: (H + V) 1/2 & 1/4 refinement only = 4 half-pel + 4 quarter-pel = 8 positions + pred_me_distortion to pa_me_distortion deviation on
    // 3                    Level 3: (H + V + D only ~ the best) 1/2 & 1/4 refinement = up to 6 half-pel + up to 6  quarter-pel = up to 12 positions + pred_me_distortion to pa_me_distortion deviation on
    // 4                    Level 4: (H + V + D) 1/2 & 1/4 refinement = 8 half-pel + 8 quarter-pel = 16 positions + pred_me_distortion to pa_me_distortion deviation on
    // 5                    Level 5: (H + V + D) 1/2 & 1/4 refinement = 8 half-pel + 8 quarter-pel = 16 positions + pred_me_distortion to pa_me_distortion deviation off
    // 6                    Level 6: (H + V + D) 1/2 & 1/4 refinement = 8 half-pel + 8 quarter-pel = 16 positions + pred_me_distortion to pa_me_distortion deviation off

    // NB: if PME_UP_TO_4_REF is ON, levels 1-5 are restricted to using max 4 ref frames, and 1/8 Pel refinement is always performed for the 8 positions for levels 1-6
    if (pcs_ptr->slice_type != I_SLICE) {

        if (pd_pass == PD_PASS_0)
            context_ptr->predictive_me_level = 0;
        else if (pd_pass == PD_PASS_1)

            context_ptr->predictive_me_level = 2;

        else

            if (sequence_control_set_ptr->static_config.pred_me == DEFAULT) {
#if !UNIFY_SC_NSC
                if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAY23_M0_ADOPTIONS
#if JUNE15_ADOPTIONS
                    if (MR_MODE)
#else
                    if (enc_mode <= ENC_M0)
#endif
                        context_ptr->predictive_me_level = 6;
                    else
#endif
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
#if M8_PRED_ME && !UPGRADE_M8
                    if (enc_mode <= ENC_M5)
#else
#if REVERT_WHITE // Pred_ME
#if JUNE17_ADOPTIONS
                    if (enc_mode <= ENC_M6)
#else
#if PRESET_SHIFITNG
                    if (enc_mode <= ENC_M5)
#else
                    if (enc_mode <= ENC_M7)
#endif
#endif
#else
                    if (enc_mode <= ENC_M8)
#endif
#endif
#else
                    if (enc_mode <= ENC_M5)
#endif
#else
                    if (enc_mode <= ENC_M3)
#endif

                        context_ptr->predictive_me_level = 1;

                    else
                        context_ptr->predictive_me_level = 0;
                else
#endif
#if JUNE23_ADOPTIONS
                    if (enc_mode <= ENC_M1)
#else
#if MAY23_M0_ADOPTIONS
                    if (enc_mode <= ENC_M0)
#else
#if MAY16_7PM_ADOPTIONS
                    if (MR_MODE)
#else
#if MAR10_ADOPTIONS
#if M1_COMBO_1 || NEW_M1_CAND
                    if (enc_mode <= ENC_M0)
#else
                    if (enc_mode <= ENC_M1)
#endif
#else
                    if (enc_mode <= ENC_M0)
#endif
#endif
#endif
#endif
                        context_ptr->predictive_me_level = 6;
#if MAR12_M8_ADOPTIONS
#if REVERT_WHITE // Pred_ME
#if JUNE26_ADOPTIONS
                    else if (enc_mode <= ENC_M5)
#else
#if JUNE25_ADOPTIONS
                    else if (enc_mode <= ENC_M6)
#else
#if JUNE17_ADOPTIONS
                    else if (enc_mode <= ENC_M4)
#else
#if PRESET_SHIFITNG
                    else if (enc_mode <= ENC_M5)
#else
                    else if (enc_mode <= ENC_M7)
#endif
#endif
#endif
#endif
#else
                    else
#endif
#if M8_PRED_ME && !UPGRADE_M8
                        if (enc_mode <= ENC_M5)
                            context_ptr->predictive_me_level = 5;
                        else
                            context_ptr->predictive_me_level = 0;
#else
#if M1_C3_ADOPTIONS
                        context_ptr->predictive_me_level = 4;
#else
                        context_ptr->predictive_me_level = 5;
#endif
#endif
#if JUNE17_ADOPTIONS
#if NEW_M8
                else if (enc_mode <= ENC_M8)
#else
#if M7_PRED_ME
                else if (enc_mode <= ENC_M7)
#else
                else if (enc_mode <= ENC_M6)
#endif
#endif
                    context_ptr->predictive_me_level = 2;
#endif
#if REVERT_WHITE // Pred_ME
            else
                context_ptr->predictive_me_level = 0;
#endif
#else
#if MAR4_M6_ADOPTIONS
                    else if (enc_mode <= ENC_M5)
#else
                    else if (enc_mode <= ENC_M3)
#endif
                        context_ptr->predictive_me_level = 5;
#if MAR3_M6_ADOPTIONS
#if MAR10_ADOPTIONS
                    else if (enc_mode <= ENC_M8)
#else
                    else if (enc_mode <= ENC_M6)
#endif
#else
                    else if (enc_mode <= ENC_M4)
#endif
                        context_ptr->predictive_me_level = 2;
                    else
                        context_ptr->predictive_me_level = 0;
#endif

            }
            else
                context_ptr->predictive_me_level =
                sequence_control_set_ptr->static_config.pred_me;
    }
    else
        context_ptr->predictive_me_level = 0;

#if ADD_SAD_AT_PME_SIGNAL
    // Level                    Settings
    // FALSE                    Use SSD at PME
    // TRUE                     Use SAD at PME
    if (pd_pass == PD_PASS_0)
        context_ptr->use_sad_at_pme = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->use_sad_at_pme = EB_FALSE;
#if !MAR18_MR_TESTS_ADOPTIONS // It was found that SAD is better than SSD for SC content
    else if (MR_MODE)
        context_ptr->use_sad_at_pme = EB_FALSE;
#endif
#if !UNIFY_SC_NSC
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->use_sad_at_pme = EB_TRUE;
#endif
    else
        context_ptr->use_sad_at_pme = EB_FALSE;
#endif

    // Derive md_staging_mode
    //
    // MD_STAGING_MODE_1
    //  ____________________________________________________________________________________________________________________________________________________________
    // |        | md_stage_0                  | md_stage_1                     | md_stage_2                              | md_stage_3                              |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    // |CLASS_0 |Prediction for Luma & Chroma |Res, T, Q, Q-1 for Luma Only    |Bypassed                                 |Res, T, Q, Q-1, T-1 or Luma & Chroma     |
    // |CLASS_6 |SAD                          |No RDOQ                         |                                         |RDOQ (f(RDOQ Level))                     |
    // |CLASS_7 |                             |No Tx Type Search               |                                         |Tx Type Search (f(Tx Type Search Level)) |
    // |        |                             |No Tx Size Search               |                                         |Tx Size Search (f(Tx Size Search Level))|
    // |        |                             |SSD @ Frequency Domain          |                                         |CFL vs. Independent                      |
    // |        |                             |                                |                                         |SSD @ Spatial Domain                     |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    // |CLASS_1 |Prediction for Luma Only     |IFS (f(IFS))                    |Bypassed                                 |Prediction for Luma & Chroma  (Best IF)  |
    // |CLASS_2 |Bilinear Only (IFS OFF)      |Res, T, Q, Q-1 for Luma Only    |                                         |Res, T, Q, Q-1, T-1 or Luma & Chroma     |
    // |CLASS_3 |SAD                          |No RDOQ                         |                                         |RDOQ (f(RDOQ Level))                     |
    // |CLASS_4 |                             |No Tx Type Search               |                                         |Tx Type Search (f(Tx Type Search Level)) |
    // |CLASS_5 |                             |No Tx Size Search               |                                         |Tx Size Search  (f(Tx Size Search Level))|
    // |CLASS_8 |                             |SSD @ Frequency Domain          |                                         |SSD @ Spatial Domain                     |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    //
    // MD_STAGING_MODE_2
    //  ____________________________________________________________________________________________________________________________________________________________
    // |        | md_stage_0                  | md_stage_1                     | md_stage_2                              | md_stage_3                              |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    // |CLASS_0 |Prediction for Luma & Chroma |Res, T, Q, Q-1 for Luma Only    |Res, T, Q, Q-1 for Luma Only             |Res, T, Q, Q-1, T-1 or Luma & Chroma     |
    // |CLASS_6 |SAD                          |No RDOQ                         |RDOQ (f(RDOQ Level))                     |RDOQ (f(RDOQ Level))                     |
    // |CLASS_7 |                             |No Tx Type Search               |Tx Type Search (f(Tx Type Search Level)) |Tx Type Search (f(Tx Type Search Level)) |
    // |        |                             |No Tx Size Search               |No Tx Size Search                        |Tx Size Search (f(Tx Size Search Level))|
    // |        |                             |SSD @ Frequency Domain          |SSD @ Frequency Domain                   |CFL vs. Independent                      |
    // |        |                             |                                |                                         |SSD @ Spatial Domain                     |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    // |CLASS_1 |Prediction for Luma Only     |IFS (f(IFS))                    |Res, T, Q, Q-1  for Luma Only            |Prediction for Luma & Chroma  (Best IF)  |
    // |CLASS_2 |Bilinear Only (IFS OFF)      |Res, T, Q, Q-1 for Luma Only    |RDOQ (f(RDOQ Level))                     |Res, T, Q, Q-1, T-1 or Luma & Chroma     |
    // |CLASS_3 |SAD                          |No RDOQ                         |Tx Type Search (f(Tx Type Search Level)) |RDOQ (f(RDOQ Level))                     |
    // |CLASS_4 |                             |No Tx Type Search               |No Tx Size Search                        |Tx Type Search (f(Tx Type Search Level)) |
    // |CLASS_5 |                             |No Tx Size Search               |SSD @ Frequency Domain                   |Tx Size Search  (f(Tx Size Search Level))|
    // |CLASS_8 |                             |SSD @ Frequency Domain          |                                         |SSD @ Spatial Domain                     |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|

    if (pd_pass == PD_PASS_0) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    }
    else
#if MAR17_ADOPTIONS
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
#else
        if (enc_mode <= ENC_M4)
            context_ptr->md_staging_mode = MD_STAGING_MODE_1;
        else
            context_ptr->md_staging_mode = MD_STAGING_MODE_0;
#endif

    // Set md staging count level
    // Level 0              minimum count = 1
    // Level 1              set towards the best possible partitioning (to further optimize)
    // Level 2              HG: breack down or look up-table(s) are required !
    if (pd_pass == PD_PASS_0) {
        context_ptr->md_staging_count_level = 0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->md_staging_count_level = 1;
    }
    else {
        context_ptr->md_staging_count_level = 2;
    }

#if !REMOVE_COMBINE_CLASS12
    // Combine MD Class1&2
    // 0                    OFF
    // 1                    ON
    if (pd_pass == PD_PASS_0) {
        context_ptr->combine_class12 = 0;
    } else if (pd_pass == PD_PASS_1) {
        context_ptr->combine_class12 = 0;
    } else
        if (sequence_control_set_ptr->static_config.combine_class_12 == DEFAULT)
        context_ptr->combine_class12 = (enc_mode <= ENC_M4) ? 0 : 1;
    else
        context_ptr->combine_class12 = sequence_control_set_ptr->static_config.combine_class_12;
#endif
    // Set interpolation filter search blk size
    // Level                Settings
    // 0                    ON for 8x8 and above
    // 1                    ON for 16x16 and above
    // 2                    ON for 32x32 and above
    if (pd_pass == PD_PASS_0) {
        context_ptr->interpolation_filter_search_blk_size = 0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->interpolation_filter_search_blk_size = 0;
    }
    else
#if MAR17_ADOPTIONS
        context_ptr->interpolation_filter_search_blk_size = 0;
#else
#if MAR4_M6_ADOPTIONS
        if (enc_mode <= ENC_M5)
#else
        if (enc_mode <= ENC_M4)
#endif
            context_ptr->interpolation_filter_search_blk_size = 0;
        else
            context_ptr->interpolation_filter_search_blk_size = 1;
#endif

    // Derive Spatial SSE Flag
    if (pd_pass == PD_PASS_0)
        context_ptr->spatial_sse_full_loop = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->spatial_sse_full_loop = EB_FALSE;
    else if (sequence_control_set_ptr->static_config.spatial_sse_fl == DEFAULT)
#if !UNIFY_SC_NSC
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR10_ADOPTIONS
            if (enc_mode <= ENC_M8)
#else
            if (enc_mode <= ENC_M6)
#endif
                context_ptr->spatial_sse_full_loop = EB_TRUE;
            else
                context_ptr->spatial_sse_full_loop = EB_FALSE;
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (enc_mode <= ENC_M8)
#else
        else if (enc_mode <= ENC_M5)
#endif
#else
        else if (enc_mode <= ENC_M4)
#endif
#else
        if (enc_mode <= ENC_M8)
#endif
            context_ptr->spatial_sse_full_loop = EB_TRUE;
        else
            context_ptr->spatial_sse_full_loop = EB_FALSE;
    else
        context_ptr->spatial_sse_full_loop =
        sequence_control_set_ptr->static_config.spatial_sse_fl;

    if (context_ptr->chroma_level <= CHROMA_MODE_1)
        context_ptr->blk_skip_decision = EB_TRUE;
    else
        context_ptr->blk_skip_decision = EB_FALSE;

    // Derive Trellis Quant Coeff Optimization Flag
    if (pd_pass == PD_PASS_0)
        context_ptr->enable_rdoq = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->enable_rdoq = EB_FALSE;
    else
        if (sequence_control_set_ptr->static_config.enable_rdoq == DEFAULT)
#if !UNIFY_SC_NSC
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR17_ADOPTIONS
#if M8_RDOQ
#if UPGRADE_M6_M7_M8
#if JUNE17_ADOPTIONS
                if (enc_mode <= ENC_M6)
#else
#if PRESET_SHIFITNG
                if (enc_mode <= ENC_M4)
#else
                if (enc_mode <= ENC_M6)
#endif
#endif
#else
                if (enc_mode <= ENC_M5)
#endif
#else
                if (enc_mode <= ENC_M8)
#endif
#else
#if MAR4_M6_ADOPTIONS
                if (enc_mode <= ENC_M5)
#else
                if (enc_mode <= ENC_M3)
#endif
#endif
                    context_ptr->enable_rdoq = EB_TRUE;
                else
#if M5_I_RDOQ
                    context_ptr->enable_rdoq = pcs_ptr->slice_type == I_SLICE ? EB_TRUE : EB_FALSE;
#else
                    context_ptr->enable_rdoq = EB_FALSE;
#endif
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
            else if (enc_mode <= ENC_M8)
#else
            else if (enc_mode <= ENC_M5)
#endif
#else
            else if (enc_mode <= ENC_M3)
#endif
#else
            if (enc_mode <= ENC_M8)
#endif
                context_ptr->enable_rdoq = EB_TRUE;
            else
                context_ptr->enable_rdoq = EB_FALSE;
        else
            context_ptr->enable_rdoq =
            sequence_control_set_ptr->static_config.enable_rdoq;

    // Derive redundant block
    if (pd_pass == PD_PASS_0)
        context_ptr->redundant_blk = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->redundant_blk = EB_TRUE;
    else
        if (sequence_control_set_ptr->static_config.enable_redundant_blk ==
            DEFAULT)
#if !UNIFY_SC_NSC
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
                if (enc_mode <= ENC_M8)
                    context_ptr->redundant_blk = EB_TRUE;
                else
                    context_ptr->redundant_blk = EB_FALSE;
#if MAR4_M8_ADOPTIONS
            else if (enc_mode <= ENC_M8)
#else
            else if (enc_mode <= ENC_M5)
#endif
#else
            if (enc_mode <= ENC_M8)
#endif
                context_ptr->redundant_blk = EB_TRUE;
            else
                context_ptr->redundant_blk = EB_FALSE;
        else
            context_ptr->redundant_blk =
            sequence_control_set_ptr->static_config.enable_redundant_blk;

    // Set edge_skp_angle_intra
    if (sequence_control_set_ptr->static_config.encoder_bit_depth == EB_8BIT)
        if (pd_pass == PD_PASS_0)
            context_ptr->edge_based_skip_angle_intra = 0;
        else if (pd_pass == PD_PASS_1)
            context_ptr->edge_based_skip_angle_intra = 1;
        else if (sequence_control_set_ptr->static_config.edge_skp_angle_intra == DEFAULT) {
#if !UNIFY_SC_NSC
#if MAR12_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAY19_ADOPTIONS
#if JUNE17_ADOPTIONS
                if (enc_mode <= ENC_M5)
#else
#if PRESET_SHIFITNG
                if (enc_mode <= ENC_M4)
#else
                if (enc_mode <= ENC_M6)
#endif
#endif
#else
#if MAY12_ADOPTIONS
                if (enc_mode <= ENC_M4)
#else
#if SHIFT_M5_SC_TO_M3
                if (enc_mode <= ENC_M2)
#else
#if PRESETS_SHIFT
                if (enc_mode <= ENC_M4)
#else
#if MAR30_ADOPTIONS
                if (enc_mode <= ENC_M7)
#else
                if (enc_mode <= ENC_M3)
#endif
#endif
#endif
#endif
#endif
                    context_ptr->edge_based_skip_angle_intra = 0;
                else
                    context_ptr->edge_based_skip_angle_intra = 1;
            else
#endif
#endif
#if APR22_ADOPTIONS
#if JUNE23_ADOPTIONS
#if JUNE25_ADOPTIONS
            if (enc_mode <= ENC_M6)
#else
            if (enc_mode <= ENC_M4)
#endif
#else
#if JUNE11_ADOPTIONS
            if (enc_mode <= ENC_M3)
#else
#if JUNE8_ADOPTIONS
            if (enc_mode <= ENC_M2)
#else
#if PRESET_SHIFITNG
            if (enc_mode <= ENC_M1)
#else
            if (enc_mode <= ENC_M2)
#endif
#endif
#endif
#endif
#else
#if APR02_ADOPTIONS
            if (MR_MODE)
#else
#if MAR10_ADOPTIONS
            if (enc_mode <= ENC_M1)
#else
            if (MR_MODE)
#endif
#endif
#endif
                context_ptr->edge_based_skip_angle_intra = 0;
            else
                context_ptr->edge_based_skip_angle_intra = 1;
        } else
            context_ptr->edge_based_skip_angle_intra =
            sequence_control_set_ptr->static_config.edge_skp_angle_intra;
    else
        context_ptr->edge_based_skip_angle_intra = 0;

    // Set prune_ref_frame_for_rec_partitions
    if (pd_pass == PD_PASS_0)
#if ON_OFF_FEATURE_MRP
        context_ptr->prune_ref_frame_for_rec_partitions = override_feature_level(context_ptr->mrp_level,0,0,0);
#else
        context_ptr->prune_ref_frame_for_rec_partitions = 0;
#endif
    else if (pd_pass == PD_PASS_1)
#if ON_OFF_FEATURE_MRP
        context_ptr->prune_ref_frame_for_rec_partitions = override_feature_level(context_ptr->mrp_level,1,0,0);
#else
        context_ptr->prune_ref_frame_for_rec_partitions = 1;
#endif
    else
        if (sequence_control_set_ptr->static_config.prune_ref_rec_part == DEFAULT)
#if UNIFY_SC_NSC
#if JUNE25_ADOPTIONS
            if (enc_mode <= ENC_M8)
#else
            if (enc_mode <= ENC_M5)
#endif
#else
#if PRESETS_SHIFT
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected ||
#if JUNE11_ADOPTIONS
                enc_mode <= ENC_M5)
#else
#if JUNE8_ADOPTIONS
                enc_mode <= ENC_M2)
#else
#if PRESET_SHIFITNG
                enc_mode <= ENC_M1)
#else
                enc_mode <= ENC_M2)
#endif
#endif
#endif
#else
#if MAR4_M3_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected ||
                enc_mode <= ENC_M3)
#else
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected ||
                enc_mode <= ENC_M1)
#endif
#endif
#endif
#if ON_OFF_FEATURE_MRP
                context_ptr->prune_ref_frame_for_rec_partitions = override_feature_level(context_ptr->mrp_level,0,0,0);
#else
                context_ptr->prune_ref_frame_for_rec_partitions = 0;
#endif
            else
#if ON_OFF_FEATURE_MRP
                context_ptr->prune_ref_frame_for_rec_partitions = override_feature_level(context_ptr->mrp_level,1,0,0);
#else
                context_ptr->prune_ref_frame_for_rec_partitions = 1;
#endif
        else
            context_ptr->prune_ref_frame_for_rec_partitions =
            sequence_control_set_ptr->static_config.prune_ref_rec_part;

#if !INTER_COMP_REDESIGN
    // Derive INTER/INTER WEDGE variance TH
    // Phoenix: Active only when inter/inter compound is on
#if MAR10_ADOPTIONS && !MAR18_MR_TESTS_ADOPTIONS
    if (MR_MODE)
        context_ptr->inter_inter_wedge_variance_th = 0;
    else
#endif
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->inter_inter_wedge_variance_th = 0;
#if MAR18_MR_TESTS_ADOPTIONS
    else if (enc_mode <= ENC_M1)
        context_ptr->inter_inter_wedge_variance_th = 0;
#endif
    else
        context_ptr->inter_inter_wedge_variance_th = 100;
#endif
#if !REMOVE_MD_EXIT
    // Derive MD Exit TH
    if (pd_pass == PD_PASS_0)
        context_ptr->md_exit_th = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_exit_th = 18;
    else
        if (MR_MODE)
            context_ptr->md_exit_th = 0;
#if MAR3_M2_ADOPTIONS
#if MAR4_M3_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (enc_mode <= ENC_M8)
#else
        else if (enc_mode <= ENC_M3)
#endif
#else
        else if (enc_mode <= ENC_M0 || (enc_mode <= ENC_M2 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#else
        else if (enc_mode <= ENC_M0 || (enc_mode <= ENC_M1 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif

            context_ptr->md_exit_th = 0;
        else
            context_ptr->md_exit_th = 18;
#endif
    // md_stage_1_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_1_cand_prune_th
    if (pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_stage_1_cand_prune_th = 75;
    else
#if UNIFY_SC_NSC
        if (enc_mode <= ENC_M1)
            context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
        else if (enc_mode <= ENC_M2)
            context_ptr->md_stage_1_cand_prune_th = 75;
        else
            context_ptr->md_stage_1_cand_prune_th = 45;
#else
#if JUNE8_ADOPTIONS
#if JUNE9_ADOPTIONS
        if (enc_mode <= ENC_M1 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && enc_mode <= ENC_M2))
#else
        if (enc_mode <= ENC_M0 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && enc_mode <= ENC_M2))
#endif
            context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
#else
#if MAY16_7PM_ADOPTIONS
        if (MR_MODE)
            context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
#if PRESET_SHIFITNG
        else if (enc_mode <= ENC_M0 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && enc_mode <= ENC_M1))
#else
        else if (enc_mode <= ENC_M0 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && enc_mode <= ENC_M2))
#endif
            context_ptr->md_stage_1_cand_prune_th =
            (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? (uint64_t)~0 : 100;
#else
#if APR22_ADOPTIONS
#if MAY12_ADOPTIONS
        if (enc_mode <= ENC_M0 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && enc_mode <= ENC_M2))
#else
        if (enc_mode <= ENC_M2 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#else
#if APR08_ADOPTIONS
        if (MR_MODE || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
#if MAR30_ADOPTIONS
        if (enc_mode <= ENC_M0 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        if (enc_mode <= ENC_M1 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#endif
#endif
            context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
#endif
#endif
#if PRESETS_SHIFT
#if JUNE17_ADOPTIONS
        else if (enc_mode <= ENC_M2 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && enc_mode <= ENC_M3))
#else
#if PRESET_SHIFITNG
        else if (enc_mode <= ENC_M2)
#else
        else if (enc_mode <= ENC_M4)
#endif
#endif
            context_ptr->md_stage_1_cand_prune_th =
            sequence_control_set_ptr->static_config.md_stage_1_cand_prune_th;
        else
            context_ptr->md_stage_1_cand_prune_th = 45;
#else
#if MAR16_M8_ADOPTIONS
        else if (enc_mode <= ENC_M7)
            context_ptr->md_stage_1_cand_prune_th =
            sequence_control_set_ptr->static_config.md_stage_1_cand_prune_th;
        else
            context_ptr->md_stage_1_cand_prune_th = 45;
#else
    else
        context_ptr->md_stage_1_cand_prune_th =
        sequence_control_set_ptr->static_config.md_stage_1_cand_prune_th;
#endif
#endif
#endif

    // md_stage_1_class_prune_th (for class removal)
    // Remove class if deviation to the best higher than TH_C
    if (pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_stage_1_class_prune_th = 100;
    else
#if PRESETS_SHIFT
#if APR25_10AM_ADOPTIONS
#if UNIFY_SC_NSC
        if (enc_mode <= ENC_M2)
#else
#if PRESET_SHIFITNG
        if (enc_mode <= ENC_M2 ||
#else
        if (enc_mode <= ENC_M4 ||
#endif
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#else
#if APR23_ADOPTIONS_2
        if (enc_mode <= ENC_M5 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        if (enc_mode <= ENC_M2 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#endif
#else
#if MAR12_ADOPTIONS
        if (enc_mode <= ENC_M3 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
#if MAR11_ADOPTIONS
        if (enc_mode <= ENC_M2 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        if (enc_mode <= ENC_M3 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#endif
#endif
            context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
        else
            context_ptr->md_stage_1_class_prune_th =
            sequence_control_set_ptr->static_config.md_stage_1_class_prune_th;

    // md_stage_2_3_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than
    // md_stage_2_3_cand_prune_th
    if (pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_3_cand_prune_th = 5;
    else
#if ADD_MRS_MODE
#if REMOVE_MR_MACRO
        if (enc_mode <= ENC_MRS)
#else
        if (MRS_MODE)
#endif
            context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
        else
#endif
#if MAY23_M0_ADOPTIONS
#if !UNIFY_SC_NSC
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if JUNE17_ADOPTIONS
            if (enc_mode <= ENC_M2)
#else
#if JUNE11_ADOPTIONS
            if (enc_mode <= ENC_M1)
#else
            if (enc_mode <= ENC_M0)
#endif
#endif
                context_ptr->md_stage_2_3_cand_prune_th = 45;
            else
                context_ptr->md_stage_2_3_cand_prune_th = 15;
#if JUNE15_ADOPTIONS
        else if (MR_MODE)
#else
        else if (enc_mode <= ENC_M0)
#endif
#else
#if REMOVE_MR_MACRO
        if (enc_mode <= ENC_MR)
#else
        if (MR_MODE)
#endif
#endif
            context_ptr->md_stage_2_3_cand_prune_th = 45;
#if JUNE11_ADOPTIONS
#if JUNE25_ADOPTIONS
        else if (enc_mode <= ENC_M8)
#else
#if JUNE23_ADOPTIONS
        else if (enc_mode <= ENC_M4)
#else
        else if (enc_mode <= ENC_M1)
#endif
#endif
            context_ptr->md_stage_2_3_cand_prune_th = 15;
#endif
        else
            context_ptr->md_stage_2_3_cand_prune_th = 5;
#else
#if MAR10_ADOPTIONS
        if (MR_MODE)
#if MAY19_ADOPTIONS
            context_ptr->md_stage_2_3_cand_prune_th = 45;
#else
            context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
#endif
#if M1_C2_ADOPTIONS
        else if (enc_mode <= ENC_M0 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
#if PRESETS_SHIFT
#if M1_COMBO_2 || M2_COMBO_3 || NEW_M1_CAND
#if MAY12_ADOPTIONS
        else if (enc_mode <= ENC_M2 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
#if APR23_ADOPTIONS_2
        else if (enc_mode <= ENC_M0 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        else if (enc_mode <= ENC_M0)
#endif
#endif
#else
        else if (enc_mode <= ENC_M2 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#else
#if MAR20_M4_ADOPTIONS
        else if (enc_mode <= ENC_M3 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        else if (enc_mode <= ENC_M4 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#endif
#endif
#else
        if (enc_mode <= ENC_M3 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
            context_ptr->md_stage_2_3_cand_prune_th = 15;
        else
            context_ptr->md_stage_2_3_cand_prune_th = 5;
#endif
    // md_stage_2_3_class_prune_th (for class removal)
    // Remove class if deviation to the best is higher than
    // md_stage_2_3_class_prune_th
    if (pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_3_class_prune_th = 25;
    else
#if UNIFY_SC_NSC
        context_ptr->md_stage_2_3_class_prune_th =
        sequence_control_set_ptr->static_config.md_stage_2_3_class_prune_th;
#else
#if MAR10_ADOPTIONS && !MAR18_MR_TESTS_ADOPTIONS
        if (MR_MODE)
            context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
        else
#endif
#if MAY12_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
                if (enc_mode <= ENC_M0)
                    context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
#if PRESET_SHIFITNG
                else if (enc_mode <= ENC_M1)
#else
                else if (enc_mode <= ENC_M2)
#endif
                    context_ptr->md_stage_2_3_class_prune_th = 50;
                else
                    context_ptr->md_stage_2_3_class_prune_th = 25;
#else
#if SHIFT_M5_SC_TO_M3
            if ((enc_mode <= ENC_M2 &&
                pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if APR23_ADOPTIONS_2
            if ((enc_mode <= ENC_M4 &&
                pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if APR23_ADOPTIONS
            if ((enc_mode <= ENC_M2 &&
                pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if MAR30_ADOPTIONS
            if ((enc_mode <= ENC_M1 &&
                pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if MAR12_ADOPTIONS
            if ((enc_mode <= ENC_M3 &&
                pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if MAR11_ADOPTIONS
            if ((enc_mode <= ENC_M2 &&
                pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
        if ((enc_mode <= ENC_M3 &&
            pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#endif
#endif
#endif
#endif
#endif
            context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
#endif
        else
            context_ptr->md_stage_2_3_class_prune_th =
            sequence_control_set_ptr->static_config.md_stage_2_3_class_prune_th;
#endif
#if COEFF_BASED_BYPASS_NSQ && !MULTI_BAND_ACTIONS
    if (pd_pass == PD_PASS_0)
        context_ptr->coeff_area_based_bypass_nsq_th = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->coeff_area_based_bypass_nsq_th = 0;
    else
#if COEFF_BASED_BYPASS_NSQ_FIX
        if (MR_MODE || pcs_ptr->slice_type == I_SLICE)
#else
        if (MR_MODE || pcs_ptr->slice_type != I_SLICE)
#endif
            context_ptr->coeff_area_based_bypass_nsq_th = 0;
        else if (enc_mode <= ENC_M0)
            context_ptr->coeff_area_based_bypass_nsq_th = 4;
        else if (enc_mode <= ENC_M1)
#if M1_TH4
            context_ptr->coeff_area_based_bypass_nsq_th = 4;
#else
            context_ptr->coeff_area_based_bypass_nsq_th = 5;
#endif
        else
            context_ptr->coeff_area_based_bypass_nsq_th = 0; // TH to be identified for M2-M8
#endif

#if MULTI_BAND_ACTIONS
    if (pd_pass == PD_PASS_0)
        context_ptr->coeff_area_based_bypass_nsq_th = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->coeff_area_based_bypass_nsq_th = 0;
    else if (pd_pass == PD_PASS_2) {
        if (enc_mode == ENC_M0)
            context_ptr->coeff_area_based_bypass_nsq_th = context_ptr->enable_area_based_cycles_allocation ? m0_nsq_cycles_reduction_th[context_ptr->sb_class] : 0;
        else if (enc_mode <= ENC_M1)
            context_ptr->coeff_area_based_bypass_nsq_th = context_ptr->enable_area_based_cycles_allocation ? m1_nsq_cycles_reduction_th[context_ptr->sb_class] : 0;
        else
            context_ptr->coeff_area_based_bypass_nsq_th = context_ptr->enable_area_based_cycles_allocation ? m1_nsq_cycles_reduction_th[context_ptr->sb_class] : 0;
    }
#endif


#if NSQ_CYCLES_REDUCTION
    // NSQ cycles reduction level: TBD
    uint8_t nsq_cycles_red_mode = 0;
    if (pd_pass == PD_PASS_0)
        nsq_cycles_red_mode = 0;
    else if (pd_pass == PD_PASS_1)
        nsq_cycles_red_mode = 0;
    else
#if JUNE26_ADOPTIONS
        if (enc_mode <= ENC_M3)
            nsq_cycles_red_mode = 0;
        else if (enc_mode <= ENC_M4)
            nsq_cycles_red_mode = (pcs_ptr->slice_type == I_SLICE) ? 0 : 14;
        else
            nsq_cycles_red_mode = 15;
#else
#if JUNE11_ADOPTIONS
        if (pcs_ptr->slice_type == I_SLICE) {
#if DISALLOW_CYCLES_REDUCTION_REF
            nsq_cycles_red_mode = 0;
#else
            if (enc_mode <= ENC_M3)
                nsq_cycles_red_mode = 0;
            else
                nsq_cycles_red_mode = 1;
#endif
        }
        else
        {
#if UNIFY_SC_NSC
            if (enc_mode <= ENC_M5)
                nsq_cycles_red_mode = 0;
#else
            if (enc_mode <= ENC_M3)
                nsq_cycles_red_mode = 0;
            else if (enc_mode <= ENC_M4)
                nsq_cycles_red_mode = 2;
#endif
            else
#if FIX_NSQ_CYCLE_RED_LEVEL
                nsq_cycles_red_mode = 14;
#else
                nsq_cycles_red_mode = 15;
#endif
        }
#else
#if ADAPTIVE_NSQ_CR
        nsq_cycles_red_mode = (pcs_ptr->slice_type == I_SLICE ) ?  0 :  0;
#else
         nsq_cycles_red_mode = 0;
#endif
#endif
#endif
#if !ENABLE_ADAPTIVE_NSQ_ALL_FRAMES
#if DISALLOW_CYCLES_REDUCTION_REF
    if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
        nsq_cycles_red_mode = 0;
#endif
#endif
    set_nsq_cycle_redcution_controls(context_ptr, nsq_cycles_red_mode);

    NsqCycleRControls*nsq_cycle_red_ctrls = &context_ptr->nsq_cycles_red_ctrls;
    // Overwrite allcation action when nsq_cycles_reduction th is higher.
#if DECOUPLE_FROM_ALLOCATION

    if (nsq_cycle_red_ctrls->enabled)
        context_ptr->nsq_cycles_reduction_th = nsq_cycle_red_ctrls->th;
    else
        context_ptr->nsq_cycles_reduction_th = 0;

#else
    if(nsq_cycle_red_ctrls->enabled)
        context_ptr->coeff_area_based_bypass_nsq_th = MAX(nsq_cycle_red_ctrls->th,context_ptr->coeff_area_based_bypass_nsq_th);
#endif
#endif
#if DEPTH_CYCLES_REDUCTION
    // Depth cycles reduction level: TBD
    uint8_t depth_cycles_red_mode = 0;
#if ADAPTIVE_DEPTH_CR
#if JUNE11_ADOPTIONS
    if (pcs_ptr->slice_type == I_SLICE) {
#if DISALLOW_CYCLES_REDUCTION_REF
        depth_cycles_red_mode = 0;
#else
        if (enc_mode <= ENC_M4)
            depth_cycles_red_mode = 0;
        else
            depth_cycles_red_mode = 6;
#endif
    }
    else {
#if !UNIFY_SC_NSC
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if JUNE15_ADOPTIONS
            if (enc_mode <= ENC_M0)
#else
            if (MR_MODE)
#endif
                depth_cycles_red_mode = 0;
#if JUNE17_ADOPTIONS
            else if (enc_mode <= ENC_M3)
                depth_cycles_red_mode = 2;
#else
            else if (enc_mode <= ENC_M1)
                depth_cycles_red_mode = 2;
            else if (enc_mode <= ENC_M2)
                depth_cycles_red_mode = 3;
#endif
            else
                depth_cycles_red_mode = 5;
        else if (enc_mode <= ENC_M0)
#if JUNE15_ADOPTIONS
            depth_cycles_red_mode = 5;
#else
            depth_cycles_red_mode = 0;
        else if (enc_mode <= ENC_M1)
            depth_cycles_red_mode = 2;
        else if (enc_mode <= ENC_M2)
            depth_cycles_red_mode = 3;
        else if (enc_mode <= ENC_M4)
            depth_cycles_red_mode = 5;
#endif
        else
            depth_cycles_red_mode = 6;
#else
#if JUNE26_ADOPTIONS
        if (enc_mode <= ENC_M0)
            depth_cycles_red_mode = 0;
        else
            depth_cycles_red_mode = 6;
#else
        depth_cycles_red_mode = 0;
#endif
#endif
    }
#else
    depth_cycles_red_mode = pcs_ptr->slice_type != I_SLICE ? 0 : 0;
#endif
#else
    depth_cycles_red_mode = pcs_ptr->slice_type != I_SLICE ? 0 : 0;
#endif
#if DISALLOW_CYCLES_REDUCTION_REF
    if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
        depth_cycles_red_mode = 0;
#endif
    set_depth_cycle_redcution_controls(context_ptr, depth_cycles_red_mode);
#endif
#if SOFT_CYCLES_REDUCTION
    uint8_t adaptive_md_cycles_level = 0;
#if !JUNE26_ADOPTIONS
    if (pd_pass == PD_PASS_2) {

        if (pcs_ptr->slice_type == I_SLICE) {
            if (enc_mode <= ENC_M4)
                adaptive_md_cycles_level = 0;
#if SOFT_CYCLES_M6M7
#if !JUNE25_ADOPTIONS
            else if (enc_mode <= ENC_M5)
#if ADMDTM5_TUNE
                adaptive_md_cycles_level = 0;
#else
                adaptive_md_cycles_level = 3;
#endif
#endif
            else
                adaptive_md_cycles_level = 4;
#else
            else
                adaptive_md_cycles_level = 3;
#endif
        }
        else {
            if (enc_mode <= ENC_M0)
                adaptive_md_cycles_level = 0;
            else if (enc_mode <= ENC_M1)
                adaptive_md_cycles_level = 1;
#if ADMDTM2_TUNE
            else if (enc_mode <= ENC_M2)
                adaptive_md_cycles_level = pcs_ptr->temporal_layer_index == 0 ? 2 : 4;
#endif
            else if (enc_mode <= ENC_M3)
#if ADMDTM3_TUNE
                adaptive_md_cycles_level = pcs_ptr->temporal_layer_index == 0 ? 2 : 4;
#else
                adaptive_md_cycles_level = 2;
#endif
            else if (enc_mode <= ENC_M4)
#if ADMDTM4_TUNE
                adaptive_md_cycles_level = pcs_ptr->temporal_layer_index == 0 ? 4 : 5;
#else
                adaptive_md_cycles_level = 4;
#endif
#if !JUNE25_ADOPTIONS
            else if (enc_mode <= ENC_M5)
#if ADMDTM5_TUNE
                adaptive_md_cycles_level = 5;
#else
                adaptive_md_cycles_level = 6;
#endif
#endif
#if SOFT_CYCLES_M6M7
#if TUNE_ADAPTIVE_MD_CR_TH
            else if (enc_mode <= ENC_M6)
                adaptive_md_cycles_level = 8;
            else
                adaptive_md_cycles_level = 9;
#else
            else if (enc_mode <= ENC_M6)
                adaptive_md_cycles_level = 7;
            else
                adaptive_md_cycles_level = 8;
#endif
#endif
        }
    }
#endif
    adaptive_md_cycles_redcution_controls(context_ptr, adaptive_md_cycles_level);
#endif
    // Weighting (expressed as a percentage) applied to
    // square shape costs for determining if a and b
    // shapes should be skipped. Namely:
    // skip HA, HB, and H4 if h_cost > (weighted sq_cost)
    // skip VA, VB, and V4 if v_cost > (weighted sq_cost)
    if (pd_pass == PD_PASS_0)
        context_ptr->sq_weight = (uint32_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->sq_weight = 100;

    else
#if ADD_MRS_MODE
#if REMOVE_MR_MACRO
        if (enc_mode <= ENC_MRS)
#else
        if (MRS_MODE)
#endif
            context_ptr->sq_weight = (uint32_t)~0;
        else
#endif
#if REMOVE_MR_MACRO
        if (enc_mode <= ENC_MR)
#else
        if (MR_MODE)
#endif
#if MAY19_ADOPTIONS
            context_ptr->sq_weight = 115;
#else
#if APR22_ADOPTIONS
            context_ptr->sq_weight = (uint32_t)~0;
#else
            context_ptr->sq_weight =
            sequence_control_set_ptr->static_config.sq_weight + 15;
#endif
#endif
        else
#if UNIFY_SC_NSC
            if (enc_mode <= ENC_M0)
                context_ptr->sq_weight = 105;
            else if (enc_mode <= ENC_M1)
                context_ptr->sq_weight = 100;
#if JUNE23_ADOPTIONS
#if JUNE25_ADOPTIONS
            else if (enc_mode <= ENC_M3)
                context_ptr->sq_weight = 95;
            else if (enc_mode <= ENC_M4)
                context_ptr->sq_weight = 85;
#else
            else if (enc_mode <= ENC_M4)
                context_ptr->sq_weight = 95;
#endif
#else
            else if (enc_mode <= ENC_M3)
                context_ptr->sq_weight = 95;
            else if (enc_mode <= ENC_M4)
                context_ptr->sq_weight = 85;
#endif
            else
                context_ptr->sq_weight = 80;
#else
#if MAR12_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if !MAR30_ADOPTIONS
#if MAR18_MR_TESTS_ADOPTIONS
                if (enc_mode <= ENC_M1)
                    context_ptr->sq_weight =
                    sequence_control_set_ptr->static_config.sq_weight + 15;
                else
#endif
#endif
#if PRESETS_SHIFT
#if NEW_M1_CAND
#if JUNE8_ADOPTIONS
#if !JUNE15_ADOPTIONS
#if M0_SQ_WEIGHT_ADOPTION
                    if (enc_mode <= ENC_M0)
                        context_ptr->sq_weight = 110;
                    else
#endif
#endif
                    if (enc_mode <= ENC_M1)
#else
#if M1_C2_ADOPTIONS
                    if (enc_mode <= ENC_M0)
#else
#if COEFF_BASED_BYPASS_NSQ && !REMOVE_SQ_WEIGHT_TOGGLING
                    if (enc_mode <= ENC_M1)
#else
                    if (enc_mode <= ENC_M0)
#endif
#endif
#endif
#else
                    if (enc_mode <= ENC_M3)
#endif
#else
#if MAR30_ADOPTIONS
                    if (enc_mode <= ENC_M4)
#else
                    if (enc_mode <= ENC_M3)
#endif
#endif
#if MAY12_ADOPTIONS
#if COEFF_BASED_BYPASS_NSQ && !REMOVE_SQ_WEIGHT_TOGGLING
                        context_ptr->sq_weight = (uint32_t)~0;
#else
#if MAY15_M0_ADOPTIONS
                        context_ptr->sq_weight = 105;
#else
                    context_ptr->sq_weight =
                    sequence_control_set_ptr->static_config.sq_weight + 10;
#endif
#endif
#else
                    context_ptr->sq_weight =
                    sequence_control_set_ptr->static_config.sq_weight + 5;
#endif
#if !MAY12_ADOPTIONS
#if NEW_M1_CAND
#if SHIFT_M5_SC_TO_M3
                else if (enc_mode <= ENC_M2)
#else
                else if (enc_mode <= ENC_M3)
#endif
                    context_ptr->sq_weight =
                    sequence_control_set_ptr->static_config.sq_weight;
#endif
#endif
#if !JUNE8_ADOPTIONS
#if M1_C2_ADOPTIONS
                else if(enc_mode <= ENC_M1)
                    context_ptr->sq_weight = 100;
#endif
#endif
#if JUNE17_ADOPTIONS
                else if (enc_mode <= ENC_M3)
                        context_ptr->sq_weight = 95;
                else if (enc_mode <= ENC_M4)
                        context_ptr->sq_weight = 85;
                else
                        context_ptr->sq_weight = 80;
#else
                else
                    context_ptr->sq_weight =
                    sequence_control_set_ptr->static_config.sq_weight - 5;
#endif
            else

#endif
#if MAR10_ADOPTIONS
#if MAR30_ADOPTIONS
#if APR02_ADOPTIONS
#if M1_COMBO_1 || NEW_M1_CAND
                if (enc_mode <= ENC_M0)
#else
                if (enc_mode <= ENC_M1)
#endif
#else
                if (enc_mode <= ENC_M1 && pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
#endif
#else
            if (enc_mode <= ENC_M1)
#endif
#else
            if (enc_mode <= ENC_M0 || (enc_mode <= ENC_M1 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected)))
#endif
#if MAY12_ADOPTIONS
#if COEFF_BASED_BYPASS_NSQ && !REMOVE_SQ_WEIGHT_TOGGLING
                context_ptr->sq_weight =(uint32_t)~0;
#else
#if M0_SQ_WEIGHT_ADOPTION && !JUNE15_ADOPTIONS
                context_ptr->sq_weight = 110;
#else
#if MAY15_M0_ADOPTIONS
                context_ptr->sq_weight = 105;
#else
                context_ptr->sq_weight =
                sequence_control_set_ptr->static_config.sq_weight + 10;
#endif
#endif
#endif
#else
                context_ptr->sq_weight =
                sequence_control_set_ptr->static_config.sq_weight + 5;
#endif
#if M1_COMBO_1 || NEW_M1_CAND
            else if (enc_mode <= ENC_M1)
                context_ptr->sq_weight =
#if M1_C2_ADOPTIONS
                100;
#else
#if M1_COMBO_3 || NEW_M1_CAND
#if COEFF_BASED_BYPASS_NSQ && !REMOVE_SQ_WEIGHT_TOGGLING
                (uint32_t)~0;
#else
                pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? sequence_control_set_ptr->static_config.sq_weight : sequence_control_set_ptr->static_config.sq_weight - 5;
#endif
#else
                sequence_control_set_ptr->static_config.sq_weight;
#endif
#endif
#endif
#if JUNE17_ADOPTIONS
            else if (enc_mode <= ENC_M3)
                context_ptr->sq_weight = 95;
            else if (enc_mode <= ENC_M4)
                context_ptr->sq_weight = 85;
            else
                context_ptr->sq_weight = 80;
#else
            else
                context_ptr->sq_weight =
                sequence_control_set_ptr->static_config.sq_weight - 5;
#endif
#endif

#if NEW_CYCLE_ALLOCATION && !DISALLOW_ALL_ACTIONS
    if (context_ptr->enable_area_based_cycles_allocation) {
        if (context_ptr->sb_class == LOW_COMPLEX_CLASS)
            context_ptr->sq_weight = pcs_ptr->parent_pcs_ptr->sc_content_detected && enc_mode <= ENC_M6 ? 70 : 50;
        else if (context_ptr->sb_class == VERY_LOW_COMPLEX_CLASS)
            context_ptr->sq_weight = 100 - (10 * context_ptr->coeffcients_area_based_cycles_allocation_level);
    }
#endif

    // nsq_hv_level  needs sq_weight to be ON
    // 0: OFF
    // 1: ON 10% + skip HA/HB/H4  or skip VA/VB/V4
    // 2: ON 10% + skip HA/HB  or skip VA/VB   ,  5% + skip H4  or skip V4
#if JUNE17_ADOPTIONS
    if (pd_pass < PD_PASS_2 || enc_mode <= ENC_M2)
#else
#if JUNE11_ADOPTIONS
    if (enc_mode <= ENC_M1 || pd_pass < PD_PASS_2 || (enc_mode <= ENC_M2 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if MAY12_ADOPTIONS
#if PRESET_SHIFITNG
    if (enc_mode <= ENC_M0 || pd_pass < PD_PASS_2 || (enc_mode <= ENC_M2 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
    if (enc_mode <= ENC_M0 || pd_pass < PD_PASS_2 || (enc_mode <= ENC_M4 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#else
#if SHIFT_M5_SC_TO_M3
    if (enc_mode <= ENC_M0 || pd_pass < PD_PASS_2 || (enc_mode <= ENC_M2 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if APR22_ADOPTIONS
    if (enc_mode <= ENC_M0 || pd_pass < PD_PASS_2 || (enc_mode <= ENC_M4 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if PRESETS_SHIFT
    if (MR_MODE || pd_pass < PD_PASS_2 || (enc_mode <= ENC_M4 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if MAR30_ADOPTIONS
    if (MR_MODE || pd_pass < PD_PASS_2 || (enc_mode <= ENC_M7 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if MAR18_MR_TESTS_ADOPTIONS
    if (MR_MODE || pd_pass < PD_PASS_2 || (enc_mode <= ENC_M3 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
    if (MR_MODE || pd_pass < PD_PASS_2)
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
        context_ptr->nsq_hv_level = 0;
#if MAR25_ADOPTIONS
    else if (enc_mode <= ENC_M8) {
#else
#if MAR17_ADOPTIONS
    else if (enc_mode <= ENC_M7) {
#else
#if MAR4_M6_ADOPTIONS
    else if (enc_mode <= ENC_M5) {
#else
    else if (enc_mode <= ENC_M3) {
#endif
#endif
#endif
        context_ptr->nsq_hv_level = 1;
        assert(context_ptr->sq_weight != (uint32_t)~0);
    }
    else {
        context_ptr->nsq_hv_level = 2;
        assert(context_ptr->sq_weight != (uint32_t)~0);
    }
    // Set pred ME full search area
#if UNIFY_SC_NSC
    if (pd_pass == PD_PASS_0) {
        context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15;
        context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15;
        context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15;
    }
    else {
#if JUNE23_ADOPTIONS
        context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M2 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
        context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M2 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
        context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M1 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
        context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M1 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#endif
    }
#else
    if (pd_pass == PD_PASS_0) {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15;
        }
    }
    else if (pd_pass == PD_PASS_1) {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15;
        }
    }
    else {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
#if JUNE9_ADOPTIONS
            context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M1 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M1 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
#if JUNE8_ADOPTIONS
            context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M0 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M0 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
#if MAY16_7PM_ADOPTIONS
            context_ptr->pred_me_full_pel_search_width =
                (enc_mode <= ENC_M0 && pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
                ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15
                : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height =
                (enc_mode <= ENC_M0 && pcs_ptr->parent_pcs_ptr->input_resolution >= INPUT_SIZE_1080p_RANGE)
                ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15
                : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
#if M1_C2_ADOPTIONS
            context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M0 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M0 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
#if PRESETS_SHIFT
            context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M2 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M2 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
#if MAR30_ADOPTIONS
            context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M3 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M3 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
#if MAR10_ADOPTIONS
            context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M1 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M1 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
            context_ptr->pred_me_full_pel_search_width = enc_mode <= ENC_M0 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = enc_mode <= ENC_M0 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#endif
#endif
#endif
#endif
#endif
#endif
#endif
        }
    }
#endif

#if !INTER_COMP_REDESIGN
    // comp_similar_mode
    // 0: OFF
    // 1: If previous similar block is not compound, do not inject compound
    // 2: If previous similar block is not compound, do not inject compound else
    // consider the compound modes up the similar s one
#if MAR17_ADOPTIONS
    if (enc_mode <= ENC_M7)
#else
#if MAR11_ADOPTIONS
    if (enc_mode <= ENC_M4)
#else
    if (enc_mode <= ENC_M3)
#endif
#endif
        context_ptr->comp_similar_mode = 1;
    else
        context_ptr->comp_similar_mode = 2;

#endif
    // Set coeff_based_nsq_cand_reduction
    if (pd_pass == PD_PASS_0)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
#if !MAY19_ADOPTIONS
#if MAR18_MR_TESTS_ADOPTIONS
    else if (MR_MODE &&
        pcs_ptr->parent_pcs_ptr->sc_content_detected == 0)
#else
#if MAR10_ADOPTIONS
    else if (MR_MODE)
#else
    else if (MR_MODE &&
        pcs_ptr->parent_pcs_ptr->sc_content_detected == 0)
#endif
#endif
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
#endif
    else
#if JUNE8_ADOPTIONS
#if REMOVE_MR_MACRO
        if (enc_mode <= ENC_MR)
#else
        if (MR_MODE)
#endif
            context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
        else
#endif
        context_ptr->coeff_based_nsq_cand_reduction = EB_TRUE;

#if OBMC_CLI
    // Set pic_obmc_level @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_pic_obmc_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_pic_obmc_level = 0;
    else
        context_ptr->md_pic_obmc_level =
        pcs_ptr->parent_pcs_ptr->pic_obmc_level;

#if OBMC_FAST
    set_obmc_controls(context_ptr, context_ptr->md_pic_obmc_level);
#endif
#else
    // Set pic_obmc_mode @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_pic_obmc_mode = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_pic_obmc_mode = 0;
    else
        context_ptr->md_pic_obmc_mode =
        pcs_ptr->parent_pcs_ptr->pic_obmc_mode;

#if OBMC_FAST
    set_obmc_controls(context_ptr,context_ptr->md_pic_obmc_mode);
#endif
#endif

    // Set enable_inter_intra @ MD
#if  CLEANUP_INTER_INTRA
    //Block level switch, has to follow the picture level
#endif
    if (pd_pass == PD_PASS_0)
        context_ptr->md_enable_inter_intra = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_enable_inter_intra = 0;
    else
        context_ptr->md_enable_inter_intra =
        pcs_ptr->parent_pcs_ptr->enable_inter_intra;

    // Set enable_paeth @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_enable_paeth = 1;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_enable_paeth = 1;
    else if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth == DEFAULT)
        context_ptr->md_enable_paeth = 1;
    else
        context_ptr->md_enable_paeth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth;

    // Set enable_smooth @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_enable_smooth = 1;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_enable_smooth = 1;
    else if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth == DEFAULT)
        context_ptr->md_enable_smooth = 1;
    else
        context_ptr->md_enable_smooth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth;

    // Set md_tx_size_search_mode @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_tx_size_search_mode = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_tx_size_search_mode = 0;
    else
        context_ptr->md_tx_size_search_mode = pcs_ptr->parent_pcs_ptr->tx_size_search_mode;

#if COEFF_BASED_TXS_BYPASS
    uint8_t txs_cycles_reduction_level = 0;
    set_txs_cycle_reduction_controls(context_ptr, txs_cycles_reduction_level);
#endif
#if OPT_BLOCK_INDICES_GEN_2 && !NEW_CYCLE_ALLOCATION
    // Update txs settings based on the sb_class
    context_ptr->md_tx_size_search_mode = (context_ptr->enable_area_based_cycles_allocation && context_ptr->sb_class == MEDIUM_COMPLEX_CLASS) ? 0 : context_ptr->md_tx_size_search_mode;
#endif
    // Set md_filter_intra_mode @ MD
#if FILTER_INTRA_CLI
    // md_filter_intra_level specifies whether filter intra would be active
    // for a given prediction candidate in mode decision.

    // md_filter_intra_level | Settings
    // 0                      | OFF
    // 1                      | ON
    if (pd_pass == PD_PASS_0)
        context_ptr->md_filter_intra_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_filter_intra_level = 0;
    else
        context_ptr->md_filter_intra_level =
        pcs_ptr->pic_filter_intra_level;
#else
    if (pd_pass == PD_PASS_0)
        context_ptr->md_filter_intra_mode = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_filter_intra_mode = 0;
    else
        context_ptr->md_filter_intra_mode =
        pcs_ptr->pic_filter_intra_mode;
#endif
#if SHUT_PALETTE_BC_PD_PASS_0_1
    // Set md_allow_intrabc @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_allow_intrabc = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_allow_intrabc = 0;
    else
        context_ptr->md_allow_intrabc = pcs_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc;

    // Set md_palette_mode @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_palette_mode = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_palette_mode = 0;
    else
        context_ptr->md_palette_mode = pcs_ptr->parent_pcs_ptr->palette_mode;
#endif
    // intra_similar_mode
    // 0: OFF
    // 1: If previous similar block is intra, do not inject any inter
    context_ptr->intra_similar_mode = 1;

#if MD_REFERENCE_MASKING
#if !SOFT_CYCLES_REDUCTION
    // Set inter_inter_distortion_based_reference_pruning
    if (pcs_ptr->slice_type != I_SLICE) {
        if (pd_pass == PD_PASS_0)
#if ON_OFF_FEATURE_MRP
            context_ptr->inter_inter_distortion_based_reference_pruning = override_feature_level(context_ptr->mrp_level,0,0,0);
#else
            context_ptr->inter_inter_distortion_based_reference_pruning = 0;
#endif
        else if (pd_pass == PD_PASS_1)
#if ON_OFF_FEATURE_MRP
            context_ptr->inter_inter_distortion_based_reference_pruning = override_feature_level(context_ptr->mrp_level,0,0,0);
#else
            context_ptr->inter_inter_distortion_based_reference_pruning = 0;
#endif
#if MAY23_M0_ADOPTIONS
#if PRUNING_PER_INTER_TYPE
#if !REDUCE_MR_COMP_CANDS
        else if (MR_MODE)
            context_ptr->inter_inter_distortion_based_reference_pruning = 0;
#endif
#if !JUNE15_ADOPTIONS
#if NEW_MRP_SETTINGS
        else if (enc_mode <= ENC_M0 && pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if ON_OFF_FEATURE_MRP
            context_ptr->inter_inter_distortion_based_reference_pruning = override_feature_level(context_ptr->mrp_level,0,0,0);
#else
            context_ptr->inter_inter_distortion_based_reference_pruning = 0;
#endif
#endif
        else if (enc_mode <= ENC_M0)
#else
        else if (MR_MODE || (enc_mode <= ENC_M0 && !pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#if ON_OFF_FEATURE_MRP
            context_ptr->inter_inter_distortion_based_reference_pruning = override_feature_level(context_ptr->mrp_level,1,0,0);
#else
            context_ptr->inter_inter_distortion_based_reference_pruning = 1;
#endif
        else
#if ON_OFF_FEATURE_MRP
            context_ptr->inter_inter_distortion_based_reference_pruning = override_feature_level(context_ptr->mrp_level,4,0,0);
#else
            context_ptr->inter_inter_distortion_based_reference_pruning = 4;
#endif
#else
       else if (enc_mode <= ENC_M0)
            context_ptr->inter_inter_distortion_based_reference_pruning = 0;
        else
            context_ptr->inter_inter_distortion_based_reference_pruning = 3;
#endif
#else
#if MAY19_ADOPTIONS
        else if (MR_MODE)
            context_ptr->inter_inter_distortion_based_reference_pruning = 0;
#endif
#if MAR25_ADOPTIONS
#if MAY16_7PM_ADOPTIONS
        else if (enc_mode <= ENC_M0 && pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
#if MAY12_ADOPTIONS
        else if (enc_mode <= ENC_M0)
#else
#if PRESETS_SHIFT
#if M1_COMBO_1 || NEW_M1_CAND
        else if (enc_mode <= ENC_M0 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        else if (enc_mode <= ENC_M2 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#else
        else if (enc_mode <= ENC_M3 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#endif
#endif
            context_ptr->inter_inter_distortion_based_reference_pruning = 0;
#if !MAY19_ADOPTIONS
#if M1_COMBO_1 || NEW_M1_CAND
        else if (enc_mode <= ENC_M1)
            context_ptr->inter_inter_distortion_based_reference_pruning = 3;
#endif
#endif
        else
            context_ptr->inter_inter_distortion_based_reference_pruning = 3;
#else
        else
            context_ptr->inter_inter_distortion_based_reference_pruning = 0; // 3 as default mode
#endif
#endif
    }
    else {
#if ON_OFF_FEATURE_MRP
        context_ptr->inter_inter_distortion_based_reference_pruning = override_feature_level(context_ptr->mrp_level,0,0,0);
#else
        context_ptr->inter_inter_distortion_based_reference_pruning = 0;
#endif
    }
    set_inter_inter_distortion_based_reference_pruning_controls(context_ptr, context_ptr->inter_inter_distortion_based_reference_pruning);
#endif
    // Set inter_intra_distortion_based_reference_pruning
    if (pcs_ptr->slice_type != I_SLICE) {
        if (pd_pass == PD_PASS_0)
            context_ptr->inter_intra_distortion_based_reference_pruning = 0;
        else if (pd_pass == PD_PASS_1)
            context_ptr->inter_intra_distortion_based_reference_pruning = 0;
        else
            context_ptr->inter_intra_distortion_based_reference_pruning = 0;  // 1 as default mode
    }
    else {
        context_ptr->inter_intra_distortion_based_reference_pruning = 0;
    }
    set_inter_intra_distortion_based_reference_pruning_controls(context_ptr, context_ptr->inter_intra_distortion_based_reference_pruning);
#endif
#if BLOCK_REDUCTION_ALGORITHM_1 || BLOCK_REDUCTION_ALGORITHM_2
    if (pd_pass == PD_PASS_0)
        context_ptr->block_based_depth_reduction_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->block_based_depth_reduction_level = 0;
    else
#if MAY23_M0_ADOPTIONS
#if JUNE11_ADOPTIONS
#if JUNE25_ADOPTIONS
        if (enc_mode <= ENC_M8)
#else
        if (enc_mode <= ENC_M5)
#endif
#else
#if JUNE8_ADOPTIONS
        if (enc_mode <= ENC_M0 || (enc_mode <= ENC_M1 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
        if (enc_mode <= ENC_M0)
#endif
#endif
            context_ptr->block_based_depth_reduction_level = 0;
        else
            context_ptr->block_based_depth_reduction_level = 2;
#else
#if MAY19_ADOPTIONS
        if (MR_MODE)
            context_ptr->block_based_depth_reduction_level = 0;
        else
#endif
#if NEW_M1_CAND
#if MAY15_M0_ADOPTIONS
        if (enc_mode <= ENC_M0 && pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        if (enc_mode <= ENC_M0)
#endif
            context_ptr->block_based_depth_reduction_level = 0;
        else
            context_ptr->block_based_depth_reduction_level = 2;
#else
#if APR02_ADOPTIONS
        if (enc_mode <= ENC_M0 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
            context_ptr->block_based_depth_reduction_level = 0;
#if PRESETS_SHIFT
        else if (enc_mode <= ENC_M2)
#else
        else if (enc_mode <= ENC_M3)
#endif
#else
        if (MR_MODE || pcs_ptr->parent_pcs_ptr->sc_content_detected)
            context_ptr->block_based_depth_reduction_level = 0;
        else if (enc_mode <= ENC_M1)
#endif
            context_ptr->block_based_depth_reduction_level = 1;
        else
            context_ptr->block_based_depth_reduction_level = 2;
#endif
#endif
    set_block_based_depth_reduction_controls(context_ptr, context_ptr->block_based_depth_reduction_level);
#endif
#if ADD_MD_NSQ_SEARCH
    if (pd_pass == PD_PASS_0)
        context_ptr->md_nsq_mv_search_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_nsq_mv_search_level = 0;
    else
#if JUNE15_ADOPTIONS
#if REMOVE_MR_MACRO
        if (enc_mode <= ENC_MRS)
#else
        if (MRS_MODE)
#endif
#else
#if ADD_MRS_MODE
        if (MRS_MODE || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
#if JUNE8_ADOPTIONS
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        if (MR_MODE || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#endif
#endif
            context_ptr->md_nsq_mv_search_level = 1;
#if JUNE23_ADOPTIONS
        else if (enc_mode <= ENC_M3)
#else
#if UNIFY_SC_NSC
        else if (enc_mode <= ENC_M1)
#else
#if JUNE15_ADOPTIONS
        else if (enc_mode <= ENC_M1 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
#if JUNE9_ADOPTIONS
        else if (enc_mode <= ENC_M1)
#else
        else if (enc_mode <= ENC_M0)
#endif
#endif
#endif
#endif
            context_ptr->md_nsq_mv_search_level = 2;
#if !JUNE23_ADOPTIONS
#if PRESET_SHIFITNG
        else if (enc_mode <= ENC_M3)
#else
        else if (enc_mode <= ENC_M5)
#endif
            context_ptr->md_nsq_mv_search_level = 3;
#endif
        else
            context_ptr->md_nsq_mv_search_level = 4;

    md_nsq_motion_search_controls(context_ptr, context_ptr->md_nsq_mv_search_level);
#endif
#if PERFORM_SUB_PEL_MD
    if (pd_pass == PD_PASS_0)
#if ADD_SKIP_INTRA_SIGNAL
#if JUNE26_ADOPTIONS
        context_ptr->md_subpel_search_level = enc_mode <= ENC_M5 ? 4 : 0;
#else
        context_ptr->md_subpel_search_level = enc_mode <= ENC_M6 ? 4 : 0;
#endif
#else
        context_ptr->md_subpel_search_level = 4;
#endif
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_subpel_search_level = 0;
#if !UNIFY_SC_NSC
    else
#if IMPROVE_SUB_PEL
        if(MR_MODE_SUB_PEL)
            context_ptr->md_subpel_search_level = 1;
        else
#endif
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
            context_ptr->md_subpel_search_level = 3;
#endif
#if JUNE9_ADOPTIONS
#if JUNE17_ADOPTIONS
        else if (enc_mode <= ENC_M0)
#else
        else if (enc_mode <= ENC_M1)
#endif
            context_ptr->md_subpel_search_level = 1;
#else
        else if (enc_mode <= ENC_M0)
            context_ptr->md_subpel_search_level = 1;
        else if (enc_mode <= ENC_M1)
            context_ptr->md_subpel_search_level = 2;
#endif
        else
            context_ptr->md_subpel_search_level = 3;
#if REMOVE_MR_MACRO
    md_subpel_search_controls(context_ptr, context_ptr->md_subpel_search_level,enc_mode);
#else
    md_subpel_search_controls(context_ptr, context_ptr->md_subpel_search_level);
#endif
#endif
    // Set max_ref_count @ MD
    if (pd_pass == PD_PASS_0)
#if ON_OFF_FEATURE_MRP
        context_ptr->md_max_ref_count = override_feature_level(context_ptr->mrp_level,4,4,1);
#else
        context_ptr->md_max_ref_count = 4;
#endif
    else if (pd_pass == PD_PASS_1)
#if ON_OFF_FEATURE_MRP
        context_ptr->md_max_ref_count = override_feature_level(context_ptr->mrp_level,1,4,1);
#else
        context_ptr->md_max_ref_count = 1;
#endif
    else
#if ON_OFF_FEATURE_MRP
        context_ptr->md_max_ref_count = override_feature_level(context_ptr->mrp_level,4,4,1);
#else
        context_ptr->md_max_ref_count = 4;
#endif
#if !PRUNING_PER_INTER_TYPE
#if ADD_BEST_CAND_COUNT_SIGNAL
    if (pd_pass == PD_PASS_0)
        context_ptr->bipred3x3_number_input_mv = 4;
    else if (pd_pass == PD_PASS_1)
        context_ptr->bipred3x3_number_input_mv = 4;
#if M1_C2_ADOPTIONS
    else if (enc_mode <= ENC_M0)
#else
    else if (enc_mode <= ENC_M2)
#endif
        context_ptr->bipred3x3_number_input_mv = 4;
    else
        context_ptr->bipred3x3_number_input_mv = 1;
#endif
#endif
    // Set md_skip_mvp_generation (and use (0,0) as MVP instead)
    if (pd_pass == PD_PASS_0)
        context_ptr->md_skip_mvp_generation = EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_skip_mvp_generation = EB_FALSE;
    else
        context_ptr->md_skip_mvp_generation = EB_FALSE;

    // Set dc_cand_only_flag
    if (pd_pass == PD_PASS_0)
        context_ptr->dc_cand_only_flag = EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->dc_cand_only_flag =
        (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
    else
        context_ptr->dc_cand_only_flag = EB_FALSE;

    // Set intra_angle_delta @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_intra_angle_delta = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_intra_angle_delta = 0;
    else if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta == DEFAULT)
        context_ptr->md_intra_angle_delta = 1;
    else
        context_ptr->md_intra_angle_delta = pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta;

    // Set disable_angle_z2_prediction_flag
    if (pd_pass == PD_PASS_0)
        context_ptr->disable_angle_z2_intra_flag = EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->disable_angle_z2_intra_flag = EB_TRUE;
    else
        context_ptr->disable_angle_z2_intra_flag = EB_FALSE;

    // Set full_cost_derivation_fast_rate_blind_flag
    if (pd_pass == PD_PASS_0)
        context_ptr->full_cost_shut_fast_rate_flag = EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->full_cost_shut_fast_rate_flag = EB_FALSE;
    else
        context_ptr->full_cost_shut_fast_rate_flag = EB_FALSE;
#if !PD0_INTER_CAND
    // Set best_me_cand_only_flag
    if (pd_pass == PD_PASS_0)
        context_ptr->best_me_cand_only_flag = EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->best_me_cand_only_flag = EB_FALSE;
    else
        context_ptr->best_me_cand_only_flag = EB_FALSE;
#endif
#if !CLEANUP_CYCLE_ALLOCATION
    // Set skip_depth
    if (pd_pass == PD_PASS_0)
        context_ptr->skip_depth = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->skip_depth = 0;
#if MAR18_MR_TESTS_ADOPTIONS
#if !APR08_ADOPTIONS
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR30_ADOPTIONS
        if (enc_mode <= ENC_M7)
#else
        if (enc_mode <= ENC_M3)
#endif
            context_ptr->skip_depth = 0;
        else
            context_ptr->skip_depth = 1;
#endif
    else
        context_ptr->skip_depth = 0;
#else
    else if (MR_MODE)
        context_ptr->skip_depth = 0;
    else
        context_ptr->skip_depth =
        pcs_ptr->parent_pcs_ptr->sc_content_detected ? 1 : 0;
#endif
#endif
#if !PERFORM_SUB_PEL_MD
    // Set perform_me_mv_1_8_pel_ref
    if (pd_pass == PD_PASS_0)
        context_ptr->perform_me_mv_1_8_pel_ref = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->perform_me_mv_1_8_pel_ref = EB_FALSE;
    else
        context_ptr->perform_me_mv_1_8_pel_ref =
        (pcs_ptr->parent_pcs_ptr->frm_hdr
            .allow_high_precision_mv);
#endif
#if ADD_SKIP_INTRA_SIGNAL
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->skip_intra = 0;
    else if (pd_pass == PD_PASS_0)
#if JUNE26_ADOPTIONS
        if (enc_mode <= ENC_M5)
#else
        if (enc_mode <= ENC_M6)
#endif
            context_ptr->skip_intra = 0;
        else
            context_ptr->skip_intra = 1;
    else
        context_ptr->skip_intra = 0;
#endif
    // skip cfl based on inter/intra cost deviation (skip if intra_cost is
    // skip_cfl_cost_dev_th % greater than inter_cost)
    if (MR_MODE)
        context_ptr->skip_cfl_cost_dev_th = (uint16_t)~0;
    else
        context_ptr->skip_cfl_cost_dev_th = 30;

    // set intra count to zero for md stage 3 if intra_cost is
    // mds3_intra_prune_th % greater than inter_cost
    if (MR_MODE)
        context_ptr->mds3_intra_prune_th = (uint16_t)~0;
    else
        context_ptr->mds3_intra_prune_th = 30;

    return return_error;
}
#else
EbErrorType signal_derivation_enc_dec_kernel_oq(SequenceControlSet * scs_ptr,
                                                PictureControlSet *  pcs_ptr,
                                                ModeDecisionContext *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    // Tx_search Level                                Settings
    // 0                                              OFF
    // 1                                              Tx search at encdec
    // 2                                              Tx search at inter-depth
    // 3                                              Tx search at full loop
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->tx_search_level = TX_SEARCH_OFF;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR2_M8_ADOPTIONS
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
#else
        if (pcs_ptr->enc_mode <= ENC_M6)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
#endif
#if MAR17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7)
#else
    else if (pcs_ptr->enc_mode <= ENC_M4)
#endif
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
#if MAR2_M8_ADOPTIONS
    else {
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    }
#else
    else if (pcs_ptr->enc_mode <= ENC_M7) {
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    } else
        context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
#endif
    // Set tx search skip weights (MAX_MODE_COST: no skipping; 0: always skipping)
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->tx_weight = MAX_MODE_COST;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
    else if (MR_MODE) // tx weight
        context_ptr->tx_weight = MAX_MODE_COST;
    else {
        if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
            context_ptr->tx_weight = MAX_MODE_COST;
#if MAR12_ADOPTIONS
        else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR18_ADOPTIONS
            context_ptr->tx_weight = MAX_MODE_COST;
#else
            if (pcs_ptr->enc_mode <= ENC_M3)
                context_ptr->tx_weight = MAX_MODE_COST;
            else
                context_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
#endif
#endif
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M1)
#else
        else if (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M0)
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else if (pcs_ptr->enc_mode <= ENC_M0)
#endif
            context_ptr->tx_weight = MAX_MODE_COST;
#if MAR3_M2_ADOPTIONS
#if MAR4_M3_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M8 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
        else if (pcs_ptr->enc_mode <= ENC_M3 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M2 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M1 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
    }

    // Set tx search reduced set falg (0: full tx set; 1: reduced tx set; 1: two tx))
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->tx_search_reduced_set = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->tx_search_reduced_set = 0;
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR17_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M8)
            context_ptr->tx_search_reduced_set = 0;
#else
        if (pcs_ptr->enc_mode <= ENC_M5)
            context_ptr->tx_search_reduced_set = 0;
        else if (pcs_ptr->enc_mode <= ENC_M6)
            if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
                context_ptr->tx_search_reduced_set = 0;
            else
                context_ptr->tx_search_reduced_set = 1;
#endif
#if MAR2_M8_ADOPTIONS
        else
            context_ptr->tx_search_reduced_set = 1;
#else
        else if (pcs_ptr->enc_mode <= ENC_M7)
            context_ptr->tx_search_reduced_set = 1;
        else
            context_ptr->tx_search_reduced_set = 2;
#endif
    else if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
        context_ptr->tx_search_reduced_set = 0;
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M8)
#else
    else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
    else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
        context_ptr->tx_search_reduced_set = 0;
    else
        context_ptr->tx_search_reduced_set = 1;
    // Interpolation search Level                     Settings
    // 0                                              OFF
    // 1                                              Interpolation search at inter-depth
    // 2                                              Interpolation search at full loop
    // 3                                              Chroma blind interpolation search at fast loop
    // 4                                              Interpolation search at fast loop
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
    } else if (MR_MODE)
        context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP;
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M8)
#else
    else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
    else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
        context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
#if !MAR10_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7)
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
        else
            context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#endif
    else
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;

    // Set Chroma Mode
    // Level                Settings
    // CHROMA_MODE_0  0     Full chroma search @ MD
    // CHROMA_MODE_1  1     Fast chroma search @ MD
    // CHROMA_MODE_2  2     Chroma blind @ MD + CFL @ EP
    // CHROMA_MODE_3  3     Chroma blind @ MD + no CFL @ EP
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->chroma_level = CHROMA_MODE_2; // or CHROMA_MODE_3
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->chroma_level = CHROMA_MODE_1;
    } else if (scs_ptr->static_config.set_chroma_mode == DEFAULT) {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR2_M7_ADOPTIONS
#if MAR10_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M8)
#else
            if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M6)
#endif
                context_ptr->chroma_level = CHROMA_MODE_1;
            else if (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                context_ptr->chroma_level = CHROMA_MODE_1;
            else
                context_ptr->chroma_level =
                    (scs_ptr->encoder_bit_depth == EB_8BIT) ? CHROMA_MODE_2 : CHROMA_MODE_3;
#if MAR17_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M7)
            context_ptr->chroma_level = CHROMA_MODE_0;
#else
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M4)
#else
        else if (pcs_ptr->enc_mode <= ENC_M1)
#endif
            context_ptr->chroma_level = CHROMA_MODE_0;
        else if (pcs_ptr->enc_mode <= ENC_M5 && pcs_ptr->temporal_layer_index == 0)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else if (pcs_ptr->enc_mode <= ENC_M5)
            context_ptr->chroma_level = CHROMA_MODE_1;
#endif
        else
            context_ptr->chroma_level =
#if ADOPT_CHROMA_MODE1_CFL_OFF
                (scs_ptr->encoder_bit_depth == EB_8BIT) ? CHROMA_MODE_1 : CHROMA_MODE_3;
#else
                (scs_ptr->encoder_bit_depth == EB_8BIT) ? CHROMA_MODE_2 : CHROMA_MODE_3;
#endif

    } else // use specified level
        context_ptr->chroma_level = scs_ptr->static_config.set_chroma_mode;
    // Chroma independent modes search
    // Level                Settings
    // 0                    post first md_stage
    // 1                    post last md_stage
    context_ptr->chroma_at_last_md_stage =
        MR_MODE ? 0 : (context_ptr->chroma_level == CHROMA_MODE_0 && !pcs_ptr->parent_pcs_ptr->sc_content_detected) ? 1 : 0;
#if M5_CHROMA_NICS
    // Chroma independent modes nics
    // Level                Settings
    // 0                    All supported modes.
    // 1                    All supported modes in  Intra picture and 4 in inter picture
    context_ptr->independent_chroma_nics = pcs_ptr->enc_mode == ENC_M5 ? 1 : 0;
#endif
#if ADDED_CFL_OFF
    // Cfl level
    // Level                Settings
    // 0                    Allow cfl
    // 1                    Disable cfl
#if ADOPT_CHROMA_MODE1_CFL_OFF
#if MAR17_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M7)
#else
    if(pcs_ptr->enc_mode <= ENC_M5)
#endif
        context_ptr->md_disable_cfl = EB_FALSE;
    else if(!pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
        context_ptr->md_disable_cfl = EB_TRUE;

    if (scs_ptr->static_config.disable_cfl_flag == 1 && context_ptr->md_disable_cfl == EB_TRUE)
        context_ptr->chroma_at_last_md_stage = 0; // Indeprndent chroma search at last MD stage is not supported when CFL is off
#endif
#if CFL_REDUCED_ALPHA
    // libaom_short_cuts_ths
    // 1                    faster than libaom
    // 2                    libaom - default
    if (pcs_ptr->enc_mode == ENC_M5)
        context_ptr->libaom_short_cuts_ths = 1;
    else
        context_ptr->libaom_short_cuts_ths = 2;
#endif
#if UV_SEARCH_MODE_INJCECTION
    // 0                    inject all supprted chroma mode
    // 1                    follow the luma injection
    context_ptr->intra_chroma_search_follows_intra_luma_injection = 1;
     if (pcs_ptr->enc_mode == ENC_M5)
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 1;
    else
        context_ptr->intra_chroma_search_follows_intra_luma_injection = 0;
#endif
    // Set the full loop escape level
    // Level                Settings
    // 0                    Off
    // 1                    On but only INTRA
    // 2                    On both INTRA and INTER
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->full_loop_escape = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->full_loop_escape = 0;
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR17_ADOPTIONS
            context_ptr->full_loop_escape = 0;
#else
        if (pcs_ptr->enc_mode <= ENC_M5)
            context_ptr->full_loop_escape = 0;
        else
            context_ptr->full_loop_escape = 2;
#endif
#if MAR2_M7_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M8)
#else
        else if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#else
    else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
        context_ptr->full_loop_escape = 0;
    else
        context_ptr->full_loop_escape = 2;

        // Set global MV injection
        // Level                Settings
        // 0                    Injection off (Hsan: but not derivation as used by MV ref derivation)
        // 1                    On
        // Note: global motion is off for frames with super-res enabled
    if (scs_ptr->static_config.enable_global_motion == EB_TRUE &&
        pcs_ptr->parent_pcs_ptr->frame_superres_enabled == EB_FALSE) {
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->global_mv_injection = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->global_mv_injection = 0;
        else
#if MAR4_M6_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR12_ADOPTIONS
#if MAR17_ADOPTIONS
                if (pcs_ptr->enc_mode <= ENC_M7)
#else
                if (pcs_ptr->enc_mode <= ENC_M3)
#endif
#else
#if MAR10_ADOPTIONS
#if MAR11_ADOPTIONS
                if (pcs_ptr->enc_mode <= ENC_M1)
#else
                if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#else
                if (pcs_ptr->enc_mode <= ENC_M3)
#endif
#endif
                    context_ptr->global_mv_injection = 1;
                else
                    context_ptr->global_mv_injection = 0;
#if MAR17_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M7)
#else
            else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M1)
#endif
            context_ptr->global_mv_injection = 1;
        else
            context_ptr->global_mv_injection = 0;
    } else
        context_ptr->global_mv_injection = 0;
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->new_nearest_injection = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->new_nearest_injection = 0;
    else
        context_ptr->new_nearest_injection = 1;

    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->new_nearest_near_comb_injection = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->new_nearest_near_comb_injection = 0;
    else if (scs_ptr->static_config.new_nearest_comb_inject == DEFAULT)
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
                context_ptr->new_nearest_near_comb_injection = 0;
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M1)
#else
        else if (pcs_ptr->enc_mode <= ENC_M0)
#endif
            context_ptr->new_nearest_near_comb_injection = 1;
        else
            context_ptr->new_nearest_near_comb_injection = 0;
    else
        context_ptr->new_nearest_near_comb_injection =
            scs_ptr->static_config.new_nearest_comb_inject;
    // Set warped motion injection
    // Level                Settings
    // 0                    OFF
    // 1                    On
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    if(frm_hdr->allow_warped_motion)
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->warped_motion_injection = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->warped_motion_injection = 1;
        else
            context_ptr->warped_motion_injection = 1;
    else
        context_ptr->warped_motion_injection = 0;
    // Set unipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->unipred3x3_injection = 0;
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->unipred3x3_injection = 2;
    } else
    if (pcs_ptr->enc_mode <= ENC_M7)
        context_ptr->unipred3x3_injection = 1;
    else
        context_ptr->unipred3x3_injection = 0;
    // Set bipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->bipred3x3_injection = 0;
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->bipred3x3_injection = 2;
    } else if (scs_ptr->static_config.bipred_3x3_inject == DEFAULT)
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR20_M4_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M3)
#else
            if (pcs_ptr->enc_mode <= ENC_M4)
#endif
                context_ptr->bipred3x3_injection = 1;
            else
#if MAR18_ADOPTIONS
                context_ptr->bipred3x3_injection = 2;
#else
                context_ptr->bipred3x3_injection = 0;
#endif
#if MAR18_ADOPTIONS
#if MAR20_M4_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M3)
#else
        else if (pcs_ptr->enc_mode <= ENC_M4)
#endif
            context_ptr->bipred3x3_injection = 1;
        else
            context_ptr->bipred3x3_injection = 2;
#else
#if MAR17_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M7)
            context_ptr->bipred3x3_injection = 1;
#else
#if MAR12_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M3)
#else
#if MAR3_M2_ADOPTIONS
#if MAR10_ADOPTIONS
#if MAR11_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M2)
#else
        else if (pcs_ptr->enc_mode <= ENC_M3)
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#endif
            context_ptr->bipred3x3_injection = 1;
#if MAR3_M6_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M6)
#else
        else if (pcs_ptr->enc_mode <= ENC_M4)
#endif
            context_ptr->bipred3x3_injection = 2;
#endif
        else
            context_ptr->bipred3x3_injection = 0;
#endif
    else
        context_ptr->bipred3x3_injection = scs_ptr->static_config.bipred_3x3_inject;

    // Level                Settings
    // 0                    Level 0: OFF
    // 1                    Level 1: 7x5 full-pel search + sub-pel refinement off
    // 2                    Level 2: 7x5 full-pel search +  (H + V) sub-pel refinement only = 4 half-pel + 4 quarter-pel = 8 positions + pred_me_distortion to pa_me_distortion deviation on
    // 3                    Level 3: 7x5 full-pel search +  (H + V + D only ~ the best) sub-pel refinement = up to 6 half-pel + up to 6  quarter-pel = up to 12 positions + pred_me_distortion to pa_me_distortion deviation on
    // 4                    Level 4: 7x5 full-pel search +  (H + V + D) sub-pel refinement = 8 half-pel + 8 quarter-pel = 16 positions + pred_me_distortion to pa_me_distortion deviation on
    // 5                    Level 5: 7x5 full-pel search +  (H + V + D) sub-pel refinement = 8 half-pel + 8 quarter-pel = 16 positions + pred_me_distortion to pa_me_distortion deviation off
    // 6                    Level 6: (H + V + D) 1/2 & 1/4 refinement = 8 half-pel + 8 quarter-pel = 16 positions + pred_me_distortion to pa_me_distortion deviation off

    // NB: levels 1-5 are restricted to using max 4 ref frames, and 1/8 Pel refinement is always performed for the 8 positions for levels 1-6
    if (pcs_ptr->slice_type != I_SLICE) {
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->predictive_me_level = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->predictive_me_level = 2;
        else if (scs_ptr->static_config.pred_me == DEFAULT) {
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
                    if (pcs_ptr->enc_mode <= ENC_M8)
#else
                    if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
                if (pcs_ptr->enc_mode <= ENC_M4)
#endif
                    context_ptr->predictive_me_level = 1;
                else
                    context_ptr->predictive_me_level = 0;
            else
#if MAR10_ADOPTIONS
                if (pcs_ptr->enc_mode <= ENC_M1)
#else
                if (pcs_ptr->enc_mode <= ENC_M0)
#endif
                context_ptr->predictive_me_level = 6;
#if MAR12_M8_ADOPTIONS
                    else
                        context_ptr->predictive_me_level = 5;
#else
#if MAR4_M6_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M5)
#else
            else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
                context_ptr->predictive_me_level = 5;
#if MAR3_M6_ADOPTIONS
#if MAR10_ADOPTIONS
                    else if (pcs_ptr->enc_mode <= ENC_M8)
#else
                    else if (pcs_ptr->enc_mode <= ENC_M6)
#endif
#else
            else if (pcs_ptr->enc_mode <= ENC_M4)
#endif
                context_ptr->predictive_me_level = 2;
            else
                context_ptr->predictive_me_level = 0;
#endif

        } else
            context_ptr->predictive_me_level = scs_ptr->static_config.pred_me;
    } else
        context_ptr->predictive_me_level = 0;

#if ADD_SAD_AT_PME_SIGNAL
    // Level                    Settings
    // FALSE                    Use SSD at PME
    // TRUE                     Use SAD at PME
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->use_sad_at_pme = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->use_sad_at_pme = EB_FALSE;
#if !MAR18_MR_TESTS_ADOPTIONS // It was found that SAD is better than SSD for SC content
    else if (MR_MODE)
        context_ptr->use_sad_at_pme = EB_FALSE;
#endif
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->use_sad_at_pme = EB_TRUE;
    else
        context_ptr->use_sad_at_pme = EB_FALSE;
#endif
    // Derive md_staging_mode
    //
    // MD_STAGING_MODE_1
    //  ____________________________________________________________________________________________________________________________________________________________
    // |        | md_stage_0                  | md_stage_1                     | md_stage_2                              | md_stage_3                              |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    // |CLASS_0 |Prediction for Luma & Chroma |Res, T, Q, Q-1 for Luma Only    |Bypassed                                 |Res, T, Q, Q-1, T-1 or Luma & Chroma     |
    // |CLASS_6 |SAD                          |No RDOQ                         |                                         |RDOQ (f(RDOQ Level))                     |
    // |CLASS_7 |                             |No Tx Type Search               |                                         |Tx Type Search (f(Tx Type Search Level)) |
    // |        |                             |No Tx Size Search               |                                         |Tx Size Search (f(Tx Size Search Level))|
    // |        |                             |SSD @ Frequency Domain          |                                         |CFL vs. Independent                      |
    // |        |                             |                                |                                         |SSD @ Spatial Domain                     |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    // |CLASS_1 |Prediction for Luma Only     |IFS (f(IFS))                    |Bypassed                                 |Prediction for Luma & Chroma  (Best IF)  |
    // |CLASS_2 |Bilinear Only (IFS OFF)      |Res, T, Q, Q-1 for Luma Only    |                                         |Res, T, Q, Q-1, T-1 or Luma & Chroma     |
    // |CLASS_3 |SAD                          |No RDOQ                         |                                         |RDOQ (f(RDOQ Level))                     |
    // |CLASS_4 |                             |No Tx Type Search               |                                         |Tx Type Search (f(Tx Type Search Level)) |
    // |CLASS_5 |                             |No Tx Size Search               |                                         |Tx Size Search  (f(Tx Size Search Level))|
    // |CLASS_8 |                             |SSD @ Frequency Domain          |                                         |SSD @ Spatial Domain                     |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    //
    // MD_STAGING_MODE_2
    //  ____________________________________________________________________________________________________________________________________________________________
    // |        | md_stage_0                  | md_stage_1                     | md_stage_2                              | md_stage_3                              |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    // |CLASS_0 |Prediction for Luma & Chroma |Res, T, Q, Q-1 for Luma Only    |Res, T, Q, Q-1 for Luma Only             |Res, T, Q, Q-1, T-1 or Luma & Chroma     |
    // |CLASS_6 |SAD                          |No RDOQ                         |RDOQ (f(RDOQ Level))                     |RDOQ (f(RDOQ Level))                     |
    // |CLASS_7 |                             |No Tx Type Search               |Tx Type Search (f(Tx Type Search Level)) |Tx Type Search (f(Tx Type Search Level)) |
    // |        |                             |No Tx Size Search               |No Tx Size Search                        |Tx Size Search (f(Tx Size Search Level))|
    // |        |                             |SSD @ Frequency Domain          |SSD @ Frequency Domain                   |CFL vs. Independent                      |
    // |        |                             |                                |                                         |SSD @ Spatial Domain                     |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    // |CLASS_1 |Prediction for Luma Only     |IFS (f(IFS))                    |Res, T, Q, Q-1  for Luma Only            |Prediction for Luma & Chroma  (Best IF)  |
    // |CLASS_2 |Bilinear Only (IFS OFF)      |Res, T, Q, Q-1 for Luma Only    |RDOQ (f(RDOQ Level))                     |Res, T, Q, Q-1, T-1 or Luma & Chroma     |
    // |CLASS_3 |SAD                          |No RDOQ                         |Tx Type Search (f(Tx Type Search Level)) |RDOQ (f(RDOQ Level))                     |
    // |CLASS_4 |                             |No Tx Type Search               |No Tx Size Search                        |Tx Type Search (f(Tx Type Search Level)) |
    // |CLASS_5 |                             |No Tx Size Search               |SSD @ Frequency Domain                   |Tx Size Search  (f(Tx Size Search Level))|
    // |CLASS_8 |                             |SSD @ Frequency Domain          |                                         |SSD @ Spatial Domain                     |
    // |________|_____________________________|________________________________|_________________________________________|_________________________________________|
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    }
    else
#if MAR17_ADOPTIONS
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
#else
        if (pcs_ptr->enc_mode <= ENC_M5)
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    else
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;
#endif
    // Set md staging count level
    // Level 0              minimum count = 1
    // Level 1              set towards the best possible partitioning (to further optimize)
    // Level 2              HG: breack down or look up-table(s) are required !
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->md_staging_count_level = 0;
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->md_staging_count_level = 1;
    } else {
        context_ptr->md_staging_count_level = 2;
    }

#if !REMOVE_COMBINE_CLASS12
    // Combine MD Class1&2
    // 0                    OFF
    // 1                    ON
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->combine_class12 = 0;
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->combine_class12 = 0;
    } else

        if (scs_ptr->static_config.combine_class_12 == DEFAULT)
        context_ptr->combine_class12 = (pcs_ptr->enc_mode <= ENC_M4) ? 0 : 1;
    else
        context_ptr->combine_class12 = scs_ptr->static_config.combine_class_12;
#endif
    // Set interpolation filter search blk size
    // Level                Settings
    // 0                    ON for 8x8 and above
    // 1                    ON for 16x16 and above
    // 2                    ON for 32x32 and above
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->interpolation_filter_search_blk_size = 0;
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->interpolation_filter_search_blk_size = 0;
    } else
#if MAR17_ADOPTIONS
        context_ptr->interpolation_filter_search_blk_size = 0;
#else
#if MAR4_M6_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M5)
#else
        if (pcs_ptr->enc_mode <= ENC_M4)
#endif
        context_ptr->interpolation_filter_search_blk_size = 0;
    else
        context_ptr->interpolation_filter_search_blk_size = 1;
#endif

    // Set PF MD
    context_ptr->pf_md_mode = PF_OFF;
    // Derive Spatial SSE Flag
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->spatial_sse_full_loop = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->spatial_sse_full_loop = EB_FALSE;
    else if (scs_ptr->static_config.spatial_sse_fl == DEFAULT)
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR10_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M8)
#else
            if (pcs_ptr->enc_mode <= ENC_M6)
#endif
                context_ptr->spatial_sse_full_loop = EB_TRUE;
            else
                context_ptr->spatial_sse_full_loop = EB_FALSE;
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M8)
#else
        else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M4)
#endif
            context_ptr->spatial_sse_full_loop = EB_TRUE;
        else
            context_ptr->spatial_sse_full_loop = EB_FALSE;
    else
        context_ptr->spatial_sse_full_loop = scs_ptr->static_config.spatial_sse_fl;

    if (context_ptr->chroma_level <= CHROMA_MODE_1)
        context_ptr->blk_skip_decision = EB_TRUE;
    else
        context_ptr->blk_skip_decision = EB_FALSE;

    // Derive enable_rdoqFlag
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->enable_rdoq = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->enable_rdoq = EB_FALSE;
    else if (scs_ptr->static_config.enable_rdoq == DEFAULT)
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR17_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M8)
#else
#if MAR4_M6_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M5)
#else
            if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#endif
                context_ptr->enable_rdoq = EB_TRUE;
            else
                context_ptr->enable_rdoq = EB_FALSE;
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M8)
#else
        else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
            context_ptr->enable_rdoq = EB_TRUE;
        else
            context_ptr->enable_rdoq = EB_FALSE;

    else
        context_ptr->enable_rdoq = scs_ptr->static_config.enable_rdoq;

    // Derive redundant block
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->redundant_blk = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->redundant_blk = EB_TRUE;
    else if (scs_ptr->static_config.enable_redundant_blk == DEFAULT)
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
            if (pcs_ptr->enc_mode <= ENC_M8)
                context_ptr->redundant_blk = EB_TRUE;
            else
                context_ptr->redundant_blk = EB_FALSE;
#if MAR4_M8_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M8)
#else
        else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
            context_ptr->redundant_blk = EB_TRUE;
        else
            context_ptr->redundant_blk = EB_FALSE;

    else
        context_ptr->redundant_blk = scs_ptr->static_config.enable_redundant_blk;

    if (scs_ptr->static_config.encoder_bit_depth == EB_8BIT)
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->edge_based_skip_angle_intra = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->edge_based_skip_angle_intra = 1;
        else if (scs_ptr->static_config.edge_skp_angle_intra == DEFAULT) {
#if MAR12_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
                if (pcs_ptr->enc_mode <= ENC_M3)
                    context_ptr->edge_based_skip_angle_intra = 0;
                else
                    context_ptr->edge_based_skip_angle_intra = 1;
            else
#endif
#if MAR10_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M1)
#else
            if (MR_MODE)
#endif
                context_ptr->edge_based_skip_angle_intra = 0;
            else
                context_ptr->edge_based_skip_angle_intra = 1;
        } else
            context_ptr->edge_based_skip_angle_intra = scs_ptr->static_config.edge_skp_angle_intra;
    else
        context_ptr->edge_based_skip_angle_intra = 0;
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->prune_ref_frame_for_rec_partitions = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->prune_ref_frame_for_rec_partitions = 1;
    else if (scs_ptr->static_config.prune_ref_rec_part == DEFAULT)
#if MAR4_M3_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected ||
                pcs_ptr->enc_mode <= ENC_M3)
#else
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected || pcs_ptr->enc_mode <= ENC_M1)
#endif
            context_ptr->prune_ref_frame_for_rec_partitions = 0;
        else
            context_ptr->prune_ref_frame_for_rec_partitions = 1;
    else
        context_ptr->prune_ref_frame_for_rec_partitions = scs_ptr->static_config.prune_ref_rec_part;

#if !INTER_COMP_REDESIGN
    // Derive INTER/INTER WEDGE variance TH
    // Phoenix: Active only when inter/inter compound is on
#if MAR10_ADOPTIONS && !MAR18_MR_TESTS_ADOPTIONS
    if (MR_MODE)
        context_ptr->inter_inter_wedge_variance_th = 0;
    else
#endif
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->inter_inter_wedge_variance_th = 0;
#if MAR18_MR_TESTS_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M1)
        context_ptr->inter_inter_wedge_variance_th = 0;
#endif
    else
        context_ptr->inter_inter_wedge_variance_th = 100;
#endif
#if !REMOVE_MD_EXIT
    // Derive MD Exit TH
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_exit_th = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_exit_th = 18;
    else if (MR_MODE)
        context_ptr->md_exit_th = 0;
#if MAR3_M2_ADOPTIONS
#if MAR4_M3_ADOPTIONS
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M8)
#else
        else if (pcs_ptr->enc_mode <= ENC_M3)
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->enc_mode <= ENC_M2 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#else
     else if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->enc_mode <= ENC_M1 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
        context_ptr->md_exit_th = 0;
    else
        context_ptr->md_exit_th = 18;
#endif
    // md_stage_1_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_1_cand_prune_th
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_1_cand_prune_th = 75;
    else if (pcs_ptr->enc_mode <= ENC_M1 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
#if MAR16_M8_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7)
        context_ptr->md_stage_1_cand_prune_th =
        scs_ptr->static_config.md_stage_1_cand_prune_th;
    else
        context_ptr->md_stage_1_cand_prune_th = 45;
#else
    else if (pcs_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_1_cand_prune_th = scs_ptr->static_config.md_stage_1_cand_prune_th;
    else
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
#endif
    // md_stage_1_class_prune_th (for class removal)
    // Remove class if deviation to the best higher than TH_C
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_1_class_prune_th = 100;

#if MAR12_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M3 ||
        pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
#if MAR11_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M2 ||
        pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
    if (pcs_ptr->enc_mode <= ENC_M1 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#endif
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
    else
        context_ptr->md_stage_1_class_prune_th = scs_ptr->static_config.md_stage_1_class_prune_th;
    // md_stage_2_3_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_2_3_cand_prune_th
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_3_cand_prune_th = 5;
#if MAR10_ADOPTIONS
        if (MR_MODE)
            context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
#if MAR20_M4_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M3 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
        else if (pcs_ptr->enc_mode <= ENC_M4 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
#else
    else if (MR_MODE)
        context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
    else if (pcs_ptr->enc_mode <= ENC_M1 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#endif
        context_ptr->md_stage_2_3_cand_prune_th = 15;
    else if (pcs_ptr->enc_mode <= ENC_M2)
        context_ptr->md_stage_2_3_cand_prune_th =
            scs_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE ? 15 : 12;
    else if (pcs_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_2_3_cand_prune_th =
            scs_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE ? 5 : 3;
    else
        context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;

    // md_stage_2_3_class_prune_th (for class removal)
    // Remove class if deviation to the best is higher than md_stage_2_3_class_prune_th
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_3_class_prune_th = 25;
    else
#if MAR10_ADOPTIONS && !MAR18_MR_TESTS_ADOPTIONS
        if (MR_MODE)
            context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
        else
#endif
#if MAR12_ADOPTIONS
            if ((pcs_ptr->enc_mode <= ENC_M3 &&
                pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
#if MAR11_ADOPTIONS
            if ((pcs_ptr->enc_mode <= ENC_M2 &&
                pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
        if ((pcs_ptr->enc_mode <= ENC_M3 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#endif
#endif
        context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
    else
        context_ptr->md_stage_2_3_class_prune_th = scs_ptr->static_config.md_stage_2_3_class_prune_th;
    // Weighting (expressed as a percentage) applied to
    // square shape costs for determining if a and b
    // shapes should be skipped. Namely:
    // skip HA and HB if h_cost > (weighted sq_cost)
    // skip VA and VB if v_cost > (weighted sq_cost)

    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->sq_weight = (uint32_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->sq_weight = 100;
    if (MR_MODE)
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight + 15;
    else
#if MAR12_ADOPTIONS
            if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#if MAR18_MR_TESTS_ADOPTIONS
                if (pcs_ptr->enc_mode <= ENC_M1)
                    context_ptr->sq_weight =
                    scs_ptr->static_config.sq_weight + 15;
                else
#endif
                if (pcs_ptr->enc_mode <= ENC_M3)
                    context_ptr->sq_weight =
                    scs_ptr->static_config.sq_weight + 5;
                else
                    context_ptr->sq_weight =
                    scs_ptr->static_config.sq_weight - 5;
            else

#endif
#if MAR10_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M1)
#else
        if (pcs_ptr->enc_mode <= ENC_M0 ||
        (pcs_ptr->enc_mode <= ENC_M1 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected)))
#endif
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight + 5;
    else
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight - 5;
    // nsq_hv_level  needs sq_weight to be ON
    // 0: OFF
    // 1: ON
    //        if H > V + TH1% then skip HA / HB / H4
    //        if V > H + TH1% then skip VA / VB / V4
    // 2: ON  for  H4 / V4 use more agressive TH2% for faster mode
#if MAR18_MR_TESTS_ADOPTIONS
    if (MR_MODE || context_ptr->pd_pass < PD_PASS_2 || (pcs_ptr->enc_mode <= ENC_M3 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
#else
    if (MR_MODE || context_ptr->pd_pass < PD_PASS_2)
#endif
        context_ptr->nsq_hv_level = 0;
#if MAR17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7) {
#else
#if MAR4_M6_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M3) {
#endif
#endif
        context_ptr->nsq_hv_level = 1;
        assert(context_ptr->sq_weight != (uint32_t)~0);
    }
    else {
        context_ptr->nsq_hv_level = 2;
        assert(context_ptr->sq_weight != (uint32_t)~0);
    }
    // Set pred ME full search area
    if (context_ptr->pd_pass == PD_PASS_0) {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15;
        }
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15;
        }
    } else {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
#if MAR10_ADOPTIONS
            context_ptr->pred_me_full_pel_search_width = pcs_ptr->enc_mode <= ENC_M1 ? PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = pcs_ptr->enc_mode <= ENC_M1 ? PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#else
            context_ptr->pred_me_full_pel_search_width = pcs_ptr->enc_mode <= ENC_M0 ?
                PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = pcs_ptr->enc_mode <= ENC_M0 ?
                PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
#endif
        }
    }
#if !INTER_COMP_REDESIGN
    //comp_similar_mode
    //0: OFF
    //1: If previous similar block is not compound, do not inject compound
    //2: If previous similar block is not compound, do not inject compound
    //   else consider the compound modes up the mode for the similar block
#if MAR17_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M7)
#else
#if MAR11_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M4)
#else
    if (pcs_ptr->enc_mode <= ENC_M3)
#endif
#endif
        context_ptr->comp_similar_mode = 1;
    else
        context_ptr->comp_similar_mode = 2;
#endif
    //intra_similar_mode
    //0: OFF
    //1: If previous similar block is intra, do not inject any inter
    context_ptr->intra_similar_mode = 1;

    // Set coeff_based_nsq_cand_reduction
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
#if MAR18_MR_TESTS_ADOPTIONS
    else if (MR_MODE &&
        pcs_ptr->parent_pcs_ptr->sc_content_detected == 0)
#else
#if MAR10_ADOPTIONS
    else if (MR_MODE)
#else
    else if (MR_MODE &&
        pcs_ptr->parent_pcs_ptr->sc_content_detected == 0)
#endif
#endif
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
    else
        context_ptr->coeff_based_nsq_cand_reduction = EB_TRUE;

    // Set pic_obmc_mode @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_pic_obmc_mode = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_pic_obmc_mode = 0;
    else
        context_ptr->md_pic_obmc_mode = pcs_ptr->parent_pcs_ptr->pic_obmc_mode;

#if OBMC_FAST
    set_obmc_controls(context_ptr, context_ptr->md_pic_obmc_mode);
#endif
    // Set enable_inter_intra @ MD
    //Block level switch, has to follow the picture level
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_enable_inter_intra = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_enable_inter_intra = 0;
    else
        context_ptr->md_enable_inter_intra = pcs_ptr->parent_pcs_ptr->enable_inter_intra;

    // Set intra_angle_delta @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_intra_angle_delta = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_intra_angle_delta = 0;
    else if (scs_ptr->static_config.intra_angle_delta == DEFAULT)
        context_ptr->md_intra_angle_delta = 1;
    else
        context_ptr->md_intra_angle_delta = scs_ptr->static_config.intra_angle_delta;

    // Set enable_paeth @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_enable_paeth = 1;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_enable_paeth = 1;
    else if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth == DEFAULT)
            context_ptr->md_enable_paeth = 1;
    else
        context_ptr->md_enable_paeth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth;

    // Set enable_smooth @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_enable_smooth = 1;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_enable_smooth = 1;
    else if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth == DEFAULT)
            context_ptr->md_enable_smooth = 1;
    else
        context_ptr->md_enable_smooth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth;

    // Set md_atb_mode @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_tx_size_search_mode = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_tx_size_search_mode = 0;
    else
        context_ptr->md_tx_size_search_mode = pcs_ptr->parent_pcs_ptr->tx_size_search_mode;

    // Set md_filter_intra_mode @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_filter_intra_mode = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_filter_intra_mode = 0;
    else
        context_ptr->md_filter_intra_mode = pcs_ptr->pic_filter_intra_mode;

    // Set max_ref_count @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_max_ref_count = 4;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_max_ref_count = 1;
    else
        context_ptr->md_max_ref_count = 4;

    // Set md_skip_mvp_generation (and use (0,0) as MVP instead)
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_skip_mvp_generation = EB_TRUE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_skip_mvp_generation = EB_FALSE;
    else
        context_ptr->md_skip_mvp_generation = EB_FALSE;

    // Set dc_cand_only_flag
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->dc_cand_only_flag = EB_TRUE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->dc_cand_only_flag = (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
    else
        context_ptr->dc_cand_only_flag = EB_FALSE;

    // Set disable_angle_z2_prediction_flag
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->disable_angle_z2_intra_flag = EB_TRUE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->disable_angle_z2_intra_flag = EB_TRUE;
    else
        context_ptr->disable_angle_z2_intra_flag = EB_FALSE;

    // Set full_cost_derivation_fast_rate_blind_flag
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->full_cost_shut_fast_rate_flag = EB_TRUE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->full_cost_shut_fast_rate_flag = EB_FALSE;
    else
        context_ptr->full_cost_shut_fast_rate_flag = EB_FALSE;

    // Set best_me_cand_only_flag
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->best_me_cand_only_flag = EB_TRUE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->best_me_cand_only_flag = EB_FALSE;
    else
        context_ptr->best_me_cand_only_flag = EB_FALSE;

    // Set skip_depth
#if SKIP_DEPTH_SYNCH
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->skip_depth = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->skip_depth = 0;
#if MAR18_MR_TESTS_ADOPTIONS
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        if (pcs_ptr->enc_mode <= ENC_M3)
            context_ptr->skip_depth = 0;
        else
            context_ptr->skip_depth = 1;
    else
        context_ptr->skip_depth = 0;
#else
    else if (MR_MODE)
        context_ptr->skip_depth = 0;
    else
        context_ptr->skip_depth =
        pcs_ptr->parent_pcs_ptr->sc_content_detected ? 1 : 0;
#endif
#else
    if(MR_MODE || context_ptr->pd_pass <= PD_PASS_1)
        context_ptr->skip_depth = 0;
    else
        context_ptr->skip_depth =
        pcs_ptr->parent_pcs_ptr->sc_content_detected ? 1 : 0;
#endif
    // Set perform_me_mv_1_8_pel_ref
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->perform_me_mv_1_8_pel_ref = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->perform_me_mv_1_8_pel_ref = EB_FALSE;
    else
        context_ptr->perform_me_mv_1_8_pel_ref = (pcs_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv);
    // Set nic_level for PD2 only
    // nic_level        nic scale factor
    // 0                1
    // 1                3/4
    // 2                2/3
    // 3                1/2
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        if (pcs_ptr->enc_mode <= ENC_M1)
            context_ptr->nic_level = 0;
        else
            context_ptr->nic_level = 2;
    else
        if (pcs_ptr->enc_mode <= ENC_M1)
            context_ptr->nic_level = 0;
        else if (pcs_ptr->enc_mode <= ENC_M1)
            context_ptr->nic_level = 2;
        else
            context_ptr->nic_level = 3;

    // skip cfl based on inter/intra cost deviation (skip if intra_cost is
    // skip_cfl_cost_dev_th % greater than inter_cost)
    if (MR_MODE)
        context_ptr->skip_cfl_cost_dev_th = (uint16_t)~0;
    else
        context_ptr->skip_cfl_cost_dev_th = 30;

    // set intra count to zero for md stage 3 if intra_cost is
    // mds3_intra_prune_th % greater than inter_cost
    if (MR_MODE)
        context_ptr->mds3_intra_prune_th = (uint16_t)~0;
    else
        context_ptr->mds3_intra_prune_th = 30;
    return return_error;
}

#endif
void copy_neighbour_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                           uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds, uint32_t sb_org_x,
                           uint32_t sb_org_y);

static void set_parent_to_be_considered(MdcSbData *results_ptr, uint32_t blk_index, int32_t sb_size,
                                        int8_t depth_step) {
    const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
    if (blk_geom->sq_size < ((sb_size == BLOCK_128X128) ? 128 : 64)) {
        //Set parent to be considered
        uint32_t parent_depth_idx_mds =
            (blk_geom->sqi_mds -
             (blk_geom->quadi - 3) * ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth]) -
            parent_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
        const BlockGeom *parent_blk_geom = get_blk_geom_mds(parent_depth_idx_mds);
        uint32_t         parent_tot_d1_blocks =
            parent_blk_geom->sq_size == 128
                ? 17
                : parent_blk_geom->sq_size > 8 ? 25 : parent_blk_geom->sq_size == 8 ? 5 : 1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < parent_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[parent_depth_idx_mds + block_1d_idx].consider_block = 1;
        }

        if (depth_step < -1)
            set_parent_to_be_considered(results_ptr, parent_depth_idx_mds, sb_size, depth_step + 1);
    }
}
#if MULTI_PASS_PD_FOR_INCOMPLETE
static void set_child_to_be_considered(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr, MdcSbData *results_ptr, uint32_t blk_index, uint32_t sb_index, int32_t sb_size,
#if TRACK_PER_DEPTH_DELTA
    int8_t pred_depth,
#endif
#if ADAPTIVE_DEPTH_CR
    uint8_t pred_sq_idx,
#endif
    int8_t depth_step) {
#else
static void set_child_to_be_considered(MdcSbData *results_ptr, uint32_t blk_index, int32_t sb_size,
                                       int8_t depth_step) {
#endif


    const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
    unsigned         tot_d1_blocks = blk_geom->sq_size == 128
        ? 17
        : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;
    if (blk_geom->sq_size > 4) {
        for (uint32_t block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[blk_index + block_1d_idx].consider_block     = 1;
            results_ptr->leaf_data_array[blk_index + block_1d_idx].refined_split_flag = EB_TRUE;
        }
        //Set first child to be considered
        uint32_t child_block_idx_1 = blk_index +
            d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
        const BlockGeom *child1_blk_geom = get_blk_geom_mds(child_block_idx_1);
        uint32_t         child1_tot_d1_blocks =
            child1_blk_geom->sq_size == 128
                ? 17
                : child1_blk_geom->sq_size > 8 ? 25 : child1_blk_geom->sq_size == 8 ? 5 : 1;

        for (uint32_t block_1d_idx = 0; block_1d_idx < child1_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].consider_block = 1;
#if TRACK_PER_DEPTH_DELTA
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].pred_depth_refinement = child1_blk_geom->depth - pred_depth;
#endif
#if ADAPTIVE_DEPTH_CR
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].pred_depth = pred_sq_idx;
#endif
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
#if MULTI_PASS_PD_FOR_INCOMPLETE
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_1])
#if TRACK_PER_DEPTH_DELTA
#if ADAPTIVE_DEPTH_CR
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_1, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
#else
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_1, sb_index, sb_size, pred_depth, depth_step > 1 ? depth_step - 1 : 1);
#endif
#else
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_1, sb_index, sb_size, depth_step > 1 ? depth_step - 1 : 1);
#endif
#else
        if (depth_step > 1)
            set_child_to_be_considered(results_ptr, child_block_idx_1, sb_size, depth_step - 1);
#endif
        //Set second child to be considered
        uint32_t child_block_idx_2 = child_block_idx_1 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom *child2_blk_geom = get_blk_geom_mds(child_block_idx_2);
        uint32_t         child2_tot_d1_blocks =
            child2_blk_geom->sq_size == 128
                ? 17
                : child2_blk_geom->sq_size > 8 ? 25 : child2_blk_geom->sq_size == 8 ? 5 : 1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < child2_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].consider_block = 1;
#if TRACK_PER_DEPTH_DELTA
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].pred_depth_refinement = child2_blk_geom->depth - pred_depth;
#endif
#if ADAPTIVE_DEPTH_CR
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].pred_depth = pred_sq_idx;
#endif
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
#if MULTI_PASS_PD_FOR_INCOMPLETE
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_2])
#if TRACK_PER_DEPTH_DELTA
#if ADAPTIVE_DEPTH_CR
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_2, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
#else
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_2, sb_index, sb_size, pred_depth, depth_step > 1 ? depth_step - 1 : 1);
#endif
#else
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_2, sb_index, sb_size, depth_step > 1 ? depth_step - 1 : 1);
#endif
#else
        if (depth_step > 1)
            set_child_to_be_considered(results_ptr, child_block_idx_2, sb_size, depth_step - 1);
#endif
        //Set third child to be considered
        uint32_t child_block_idx_3 = child_block_idx_2 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom *child3_blk_geom = get_blk_geom_mds(child_block_idx_3);
        uint32_t         child3_tot_d1_blocks =
            child3_blk_geom->sq_size == 128
                ? 17
                : child3_blk_geom->sq_size > 8 ? 25 : child3_blk_geom->sq_size == 8 ? 5 : 1;

        for (uint32_t block_1d_idx = 0; block_1d_idx < child3_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].consider_block = 1;
#if TRACK_PER_DEPTH_DELTA
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].pred_depth_refinement = child3_blk_geom->depth - pred_depth;
#endif
#if ADAPTIVE_DEPTH_CR
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].pred_depth = pred_sq_idx;
#endif
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }

#if MULTI_PASS_PD_FOR_INCOMPLETE
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_3])
#if TRACK_PER_DEPTH_DELTA
#if ADAPTIVE_DEPTH_CR
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_3, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
#else
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_3, sb_index, sb_size, pred_depth, depth_step > 1 ? depth_step - 1 : 1);
#endif
#else
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_3, sb_index, sb_size, depth_step > 1 ? depth_step - 1 : 1);
#endif
#else
        if (depth_step > 1)
            set_child_to_be_considered(results_ptr, child_block_idx_3, sb_size, depth_step - 1);
#endif
        //Set forth child to be considered
        uint32_t child_block_idx_4 = child_block_idx_3 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom *child4_blk_geom = get_blk_geom_mds(child_block_idx_4);
        uint32_t         child4_tot_d1_blocks =
            child4_blk_geom->sq_size == 128
                ? 17
                : child4_blk_geom->sq_size > 8 ? 25 : child4_blk_geom->sq_size == 8 ? 5 : 1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < child4_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].consider_block = 1;
#if TRACK_PER_DEPTH_DELTA
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].pred_depth_refinement = child4_blk_geom->depth - pred_depth;
#endif
#if ADAPTIVE_DEPTH_CR
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].pred_depth = pred_sq_idx;
#endif
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
#if MULTI_PASS_PD_FOR_INCOMPLETE
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_4])
#if TRACK_PER_DEPTH_DELTA
#if ADAPTIVE_DEPTH_CR
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_4, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
#else
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_4, sb_index, sb_size, pred_depth, depth_step > 1 ? depth_step - 1 : 1);
#endif
#else
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_4, sb_index, sb_size, depth_step > 1 ? depth_step - 1 : 1);
#endif
#else
        if (depth_step > 1)
            set_child_to_be_considered(results_ptr, child_block_idx_4, sb_size, depth_step - 1);
#endif
    }
}
#if OPT_BLOCK_INDICES_GEN_0
static void build_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
    uint32_t sb_index) {

    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
    results_ptr->leaf_count = 0;
    uint32_t blk_index = 0;
    uint32_t d1_blocks_accumlated, tot_d1_blocks = 0, d1_block_idx;

    while (blk_index < scs_ptr->max_block_cnt) {
        tot_d1_blocks = 0;
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

        // SQ/NSQ block(s) filter based on the SQ size
        uint8_t is_block_tagged =
#if REMOVE_UNUSED_CODE_PH2
            (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE)
#else
            (blk_geom->sq_size == 128 && (pcs_ptr->slice_type == I_SLICE || pcs_ptr->parent_pcs_ptr->sb_64x64_simulated))
#endif
            ? 0
            : 1;

        // split_flag is f(min_sq_size)
        int32_t min_sq_size = (context_ptr->disallow_4x4) ? 8 : 4;

        // SQ/NSQ block(s) filter based on the block validity
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_block_tagged) {
#if OPT_BLOCK_INDICES_GEN_2
#if OPT_BLOCK_INDICES_GEN_3
            tot_d1_blocks = (context_ptr->md_disallow_nsq) ||
#if NO_NSQ_ABOVE
                (blk_geom->sq_size >= 64 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_64x64) ||
                (blk_geom->sq_size >= 32 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_32x32) ||
                (blk_geom->sq_size >= 16 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_16x16) ||
#endif
#if NO_NSQ_B32
                (blk_geom->sq_size <= 64 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_64x64) ||
                (blk_geom->sq_size <= 32 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_32x32) ||
#endif
                (blk_geom->sq_size <= 8 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_8x8) ||
                (blk_geom->sq_size <= 16 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_16x16) ? 1 :
#else
            tot_d1_blocks = (context_ptr->md_disallow_nsq) ? 1 :
#endif
#else
            tot_d1_blocks = (pcs_ptr->parent_pcs_ptr->disallow_nsq) ? 1 :
#endif
#if OPT_BLOCK_INDICES_GEN_3
                (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_all_non_hv_nsq_blocks_below_16x16) ? 5 :
                (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_all_h4_v4_blocks_below_16x16) ? 17 :
#endif
                blk_geom->sq_size == 128
                ? 17
                : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;

#if NO_AB_HV4
            if (pcs_ptr->parent_pcs_ptr->disallow_HVA_HVB_HV4)
                tot_d1_blocks = MIN(5, tot_d1_blocks);

            if (pcs_ptr->parent_pcs_ptr->disallow_HV4)
                tot_d1_blocks = MIN(17, tot_d1_blocks);
#endif
            d1_blocks_accumlated = 0;
#if SOFT_CYCLES_REDUCTION
            for (d1_block_idx = 0; d1_block_idx < tot_d1_blocks; d1_block_idx++) {
                if (results_ptr->leaf_data_array[blk_index + d1_block_idx].consider_block) {
                    AMdCycleRControls*adaptive_md_cycles_red_ctrls = &context_ptr->admd_cycles_red_ctrls;
                    if (adaptive_md_cycles_red_ctrls->enabled) {
                        if (adaptive_md_cycles_red_ctrls->skip_nsq_th) {
                            const BlockGeom *nsq_blk_geom = get_blk_geom_mds(blk_index + d1_block_idx);
                            if (nsq_blk_geom->shape != PART_N) {
                                int8_t pred_depth_refinement = results_ptr->leaf_data_array[blk_index + d1_block_idx].pred_depth_refinement;
                                pred_depth_refinement = MIN(pred_depth_refinement, 2);
                                pred_depth_refinement = MAX(pred_depth_refinement, -2);
                                pred_depth_refinement += 2;
                                if (context_ptr->ad_md_prob[pred_depth_refinement][nsq_blk_geom->shape] < adaptive_md_cycles_red_ctrls->skip_nsq_th)
                                    results_ptr->leaf_data_array[blk_index + d1_block_idx].consider_block = 0;
                            }
                        }
                    }
                }
            }
            d1_blocks_accumlated = 0;
#endif

            for (d1_block_idx = 0; d1_block_idx < tot_d1_blocks; d1_block_idx++)
                d1_blocks_accumlated +=
                results_ptr->leaf_data_array[blk_index + d1_block_idx].consider_block ? 1 : 0;

            for (uint32_t idx = 0; idx < tot_d1_blocks; ++idx) {
                if (results_ptr->leaf_data_array[blk_index].consider_block) {

                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = blk_index;
#if SOFT_CYCLES_REDUCTION
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = d1_blocks_accumlated;
#else
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = tot_d1_blocks;
#endif
#if TRACK_PER_DEPTH_DELTA
                    results_ptr->leaf_data_array[results_ptr->leaf_count].final_pred_depth_refinement = results_ptr->leaf_data_array[blk_index].pred_depth_refinement;
                    if (results_ptr->leaf_data_array[results_ptr->leaf_count].final_pred_depth_refinement == -8)
                        printf("final_pred_depth_refinement error\n");
#endif
#if ADAPTIVE_DEPTH_CR
                    results_ptr->leaf_data_array[results_ptr->leaf_count].final_pred_depth = results_ptr->leaf_data_array[blk_index].pred_depth;
                    if (results_ptr->leaf_data_array[results_ptr->leaf_count].final_pred_depth == -8)
                        printf("final_pred_depth error\n");

#endif
                    results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = results_ptr->leaf_data_array[blk_index].refined_split_flag;

                }
                blk_index++;
            }
            blk_index +=
                (d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] -
                    tot_d1_blocks);
        }
        else {
            blk_index +=
                (blk_geom->sq_size > min_sq_size)
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
        }
    }
}
#else
#if DEPTH_PART_CLEAN_UP
static void build_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
    uint32_t sb_index) {
    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
#else
static void build_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                   uint32_t sb_index) {
    MdcSbData *results_ptr  = &pcs_ptr->mdc_sb_array[sb_index];
#endif
    results_ptr->leaf_count = 0;
    uint32_t blk_index      = 0;
    uint32_t d1_blocks_accumlated, tot_d1_blocks = 0, d1_block_idx;

    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
        EbBool           split_flag = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        //if the parent sq is inside inject this block
        uint8_t is_blk_allowed =
            pcs_ptr->slice_type != I_SLICE ? 1 : (blk_geom->sq_size < 128) ? 1 : 0;
        //init consider block flag
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed) {
            tot_d1_blocks = blk_geom->sq_size == 128
                                ? 17
                                : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;
            d1_blocks_accumlated = 0;
            for (d1_block_idx = 0; d1_block_idx < tot_d1_blocks; d1_block_idx++)
                d1_blocks_accumlated +=
                    results_ptr->leaf_data_array[blk_index + d1_block_idx].consider_block ? 1 : 0;

            for (uint32_t idx = 0; idx < tot_d1_blocks; ++idx) {
                if (results_ptr->leaf_data_array[blk_index].consider_block) {
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks =
                        d1_blocks_accumlated;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].leaf_index =
                        0; //valid only for square 85 world. will be removed.
                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = blk_index;
                    split_flag = results_ptr->leaf_data_array[blk_index].refined_split_flag;
                    results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = split_flag;
                }
                blk_index++;
            }
        }
#if MULTI_PASS_PD_FOR_INCOMPLETE
        else {
            blk_index +=
                split_flag
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
        }
#else
        blk_index +=
            split_flag
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] -
                      tot_d1_blocks
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] -
                      tot_d1_blocks;
#endif
    }
}
#endif
uint64_t  pd_level_tab[2][9][2][3] =
{
    {
        // Thresholds to use if block is screen content or an I-slice
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}}
    } ,
    {
        // Thresholds to use if block is not screen content or an I-slice
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}},
        {{5,0,0},{5,0,0}}
    }
};
#if FIX_MR_PD1
uint64_t  mr_pd_level_tab[2][9][2][3] =
{
    {
        // Thresholds to use if block is screen content or an I-slice
        {{200,200,200},{200,200,200}},
        {{200,200,200},{200,200,200}},
        {{200,200,200},{200,200,200}},
        {{200,200,200},{200,200,200}},
        {{200,200,200},{200,200,200}},
        {{200,200,200},{200,200,200}},
        {{200,200,200},{200,200,200}},
        {{200,200,200},{200,200,200}},
        {{200,200,200},{200,200,200}},
    } ,
    {
        // Thresholds to use if block is not screen content or an I-slice
        {{100,10,10},{100,10,10}},
        {{100,10,10},{100,10,10}},
        {{100,10,10},{100,10,10}},
        {{100,10,10},{100,10,10}},
        {{100,10,10},{100,10,10}},
        {{100,10,10},{100,10,10}},
        {{100,10,10},{100,10,10}},
        {{100,10,10},{100,10,10}},
        {{100,10,10},{100,10,10}},
    }
};
#endif
void derive_start_end_depth(PictureControlSet *pcs_ptr, SuperBlock *sb_ptr, uint32_t sb_size,
                            int8_t *s_depth, int8_t *e_depth, const BlockGeom *blk_geom) {
    uint8_t encode_mode = pcs_ptr->parent_pcs_ptr->enc_mode;

    int8_t start_depth = sb_size == BLOCK_128X128 ? 0 : 1;
    int8_t end_depth   = 5;
    int8_t depth       = blk_geom->depth + start_depth;

    int8_t depthp1 = MIN(depth + 1, end_depth);
    int8_t depthp2 = MIN(depth + 2, end_depth);
    int8_t depthp3 = MIN(depth + 3, end_depth);

    uint8_t depthm1 = MAX(depth - 1, start_depth);
    uint8_t depthm2 = MAX(depth - 2, start_depth);
    uint8_t depthm3 = MAX(depth - 3, start_depth);

    uint64_t max_distance = 0xFFFFFFFFFFFFFFFF;

#if UNIFY_SC_NSC
    mth01 = pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][0][0];
    mth02 = pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][0][1];
    mth03 = pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][0][2];
    pth01 = pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][1][0];
    pth02 = pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][1][1];
    pth03 = pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][1][2];
#else
    uint64_t mth01 = pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][0][0];
    uint64_t mth02 = pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][0][1];
    uint64_t mth03 = pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][0][2];
    uint64_t pth01 = pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][1][0];
    uint64_t pth02 = pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][1][1];
    uint64_t pth03 = pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][1][2];
#endif
#if FIX_MR_PD1
#if MR_MODE_FOR_PIC_MULTI_PASS_PD_MODE_1
#if REMOVE_MR_MACRO
    if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0 || pcs_ptr->parent_pcs_ptr->enc_mode <= ENC_MR) {
#else
    if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0 || MR_MODE_MULTI_PASS_PD) {
#endif
#else
    if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0) {
#endif
#if UNIFY_SC_NSC
        mth01 = mr_pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][0][0];
        mth02 = mr_pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][0][1];
        mth03 = mr_pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][0][2];
        pth01 = mr_pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][1][0];
        pth02 = mr_pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][1][1];
        pth03 = mr_pd_level_tab[pcs_ptr->slice_type != I_SLICE][encode_mode][1][2];
#else
        mth01 = mr_pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][0][0];
        mth02 = mr_pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][0][1];
        mth03 = mr_pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][0][2];
        pth01 = mr_pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][1][0];
        pth02 = mr_pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][1][1];
        pth03 = mr_pd_level_tab[!pcs_ptr->parent_pcs_ptr->sc_content_detected &&
                                  pcs_ptr->slice_type != I_SLICE][encode_mode][1][2];
#endif
    }
#endif
    uint64_t dist_001 =
        sb_ptr->depth_cost[depth] == 0
            ? max_distance
            : sb_ptr->depth_cost[depthp1] <= sb_ptr->depth_cost[depth]
                  ? 0
                  : (((int64_t)sb_ptr->depth_cost[depthp1] - (int64_t)sb_ptr->depth_cost[depth]) *
                     100) /
                        sb_ptr->depth_cost[depth];

    uint64_t dist_100 =
        sb_ptr->depth_cost[depth] == 0
            ? max_distance
            : sb_ptr->depth_cost[depthm1] <= sb_ptr->depth_cost[depth]
                  ? 0
                  : (((int64_t)sb_ptr->depth_cost[depthm1] - (int64_t)sb_ptr->depth_cost[depth]) *
                     100) /
                        sb_ptr->depth_cost[depth];

    uint64_t dist_002 =
        sb_ptr->depth_cost[depth] == 0
            ? max_distance
            : sb_ptr->depth_cost[depthp2] <= sb_ptr->depth_cost[depth]
                  ? 0
                  : (((int64_t)sb_ptr->depth_cost[depthp2] - (int64_t)sb_ptr->depth_cost[depth]) *
                     100) /
                        sb_ptr->depth_cost[depth];

    uint64_t dist_200 =
        sb_ptr->depth_cost[depth] == 0
            ? max_distance
            : sb_ptr->depth_cost[depthm2] <= sb_ptr->depth_cost[depth]
                  ? 0
                  : (((int64_t)sb_ptr->depth_cost[depthm2] - (int64_t)sb_ptr->depth_cost[depth]) *
                     100) /
                        sb_ptr->depth_cost[depth];

    uint64_t dist_003 =
        sb_ptr->depth_cost[depth] == 0
            ? max_distance
            : sb_ptr->depth_cost[depthp3] <= sb_ptr->depth_cost[depth]
                  ? 0
                  : (((int64_t)sb_ptr->depth_cost[depthp3] - (int64_t)sb_ptr->depth_cost[depth]) *
                     100) /
                        sb_ptr->depth_cost[depth];

    uint64_t dist_300 =
        sb_ptr->depth_cost[depth] == 0
            ? max_distance
            : sb_ptr->depth_cost[depthm3] <= sb_ptr->depth_cost[depth]
                  ? 0
                  : (((int64_t)sb_ptr->depth_cost[depthm3] - (int64_t)sb_ptr->depth_cost[depth]) *
                     100) /
                        sb_ptr->depth_cost[depth];

    if (dist_300 < mth03)
        *s_depth = -3;
    else if (dist_200 < mth02)
        *s_depth = -2;
    else if (dist_100 < mth01)
        *s_depth = -1;
    else
        *s_depth = 0;

    if (dist_003 < pth03)
        *e_depth = 3;
    else if (dist_002 < pth02)
        *e_depth = 2;
    else if (dist_001 < pth01)
        *e_depth = 1;
    else
        *e_depth = 0;
#if DEPTH_PART_CLEAN_UP && !OPT_BLOCK_INDICES_GEN_1  // refinement rest
    *s_depth =  MAX((sb_size == BLOCK_128X128 && pcs_ptr->parent_pcs_ptr->sb_64x64_simulated ? 1 : 0) - blk_geom->depth, *s_depth);
    *e_depth = MIN((pcs_ptr->parent_pcs_ptr->disallow_4x4 ? 4 :5) - blk_geom->depth, *e_depth);
#endif
}

static uint64_t generate_best_part_cost(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t  blk_index = 0;
    uint64_t best_part_cost = 0;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        // if the parent square is inside inject this block
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        // derive split_flag
        EbBool split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
        if (scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE)
                    best_part_cost += context_ptr->md_local_blk_unit[blk_index].cost;
            }
        }
        blk_index += split_flag ?
            d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] :
            ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
    return best_part_cost;
}
#if ADAPTIVE_TXT_CR
void generate_statistics_txt(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t blk_index = 0;
    uint32_t part_cnt[TXT_DEPTH_DELTA_NUM][TX_TYPES] = { {0},{0} };
    EbBool split_flag;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
        if (scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    if (context_ptr->md_local_blk_unit[blk_index].avail_blk_flag) {
                        uint8_t part_idx = context_ptr->md_blk_arr_nsq[blk_index].part;
                        int8_t pred_depth_refinement = context_ptr->md_local_blk_unit[blk_geom->sqi_mds].pred_depth_refinement;
                        // Set the bounds of pred_depth_refinement for array indexing
                        pred_depth_refinement = MIN(pred_depth_refinement, 1);
                        pred_depth_refinement = MAX(pred_depth_refinement, -1);
                        // Add one b/c starts at -1 (need proper offset for array)
                        // Only track whether the refinement is positive or negative
                        pred_depth_refinement++;
                        if (pred_depth_refinement < 0 || pred_depth_refinement >(TXT_DEPTH_DELTA_NUM - 1))
                            printf("pred_depth_refinement array idx error\t%d\n", pred_depth_refinement);
                        // Select the best partition, blk_index refers to the SQ block
                        uint32_t best_idx = 0;
                        uint32_t blks_in_best = 0;
                        switch (part_idx) {
                        case PARTITION_NONE:
                            best_idx = blk_index;
                            blks_in_best = 1;
                            break;
                        case PARTITION_HORZ:
                            best_idx = blk_index + 1;
                            blks_in_best = 2;
                            break;
                        case PARTITION_VERT:
                            best_idx = blk_index + 3;
                            blks_in_best = 2;
                            break;
                        case PARTITION_HORZ_A:
                            best_idx = blk_index + 5;
                            blks_in_best = 3;
                            break;
                        case PARTITION_HORZ_B:
                            best_idx = blk_index + 8;
                            blks_in_best = 3;
                            break;
                        case PARTITION_VERT_A:
                            best_idx = blk_index + 11;
                            blks_in_best = 3;
                            break;
                        case PARTITION_VERT_B:
                            best_idx = blk_index + 14;
                            blks_in_best = 3;
                            break;
                        case PARTITION_HORZ_4:
                            best_idx = blk_index + 17;
                            blks_in_best = 4;
                            break;
                        case PARTITION_VERT_4:
                            best_idx = blk_index + 21;
                            blks_in_best = 4;
                            break;
                        default:
                            assert(0);
                            break;
                        }
                        // Loop over the blocks in the best partition
                        for (uint32_t curr_idx = best_idx; curr_idx < best_idx + blks_in_best; curr_idx++) {
                            // Use the info of the best partition, not the square (only partition, cost,
                            // and split_flag are updated in the SQ block as the best
                            const BlockGeom * best_blk_geom = get_blk_geom_mds(curr_idx);
                            uint8_t tx_depth = context_ptr->md_blk_arr_nsq[curr_idx].tx_depth;
#if STATS_TX_TYPES_FIX
                            for (uint8_t txb_itr = 0; txb_itr < best_blk_geom->txb_count[tx_depth]; txb_itr++) {
#else
                            for (uint8_t txb_itr = 0; txb_itr < best_blk_geom[curr_idx].txb_count[tx_depth]; txb_itr++) {
#endif
                                uint8_t tx_type = context_ptr->md_blk_arr_nsq[curr_idx].txb_array[txb_itr].transform_type[PLANE_TYPE_Y];
                                uint32_t count_unit = (best_blk_geom->tx_width[tx_depth][txb_itr] * best_blk_geom->tx_height[tx_depth][txb_itr]); // count the area, not just the occurence
                                part_cnt[pred_depth_refinement][tx_type] += count_unit;
                            }
                        }
                    }
                }
            }
        }
        blk_index += split_flag ?
            d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] :
            ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
    for (uint8_t depth_delta = 0; depth_delta < TXT_DEPTH_DELTA_NUM; depth_delta++)
        for (uint8_t txs_idx = 0; txs_idx < TX_TYPES; txs_idx++)
            context_ptr->txt_cnt[depth_delta][txs_idx] += part_cnt[depth_delta][txs_idx];
}
#endif
#if ADAPTIVE_NSQ_CR || ADAPTIVE_DEPTH_CR
Part part_to_shape[NUMBER_OF_SHAPES] = {
    PART_N,
    PART_H,
    PART_V,
    PART_S,
    PART_HA,
    PART_HB,
    PART_VA,
    PART_VB,
    PART_H4,
    PART_V4
};
#endif
#if ADAPTIVE_DEPTH_CR
void generate_statistics_depth(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t blk_index = 0;
    // init stat
    EbBool split_flag;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
        if (scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    if (context_ptr->md_local_blk_unit[blk_index].avail_blk_flag) {
                        int8_t pred_depth_refinement = context_ptr->md_local_blk_unit[blk_geom->sqi_mds].pred_depth_refinement;
                        pred_depth_refinement = MIN(pred_depth_refinement, 1);
                        pred_depth_refinement = MAX(pred_depth_refinement, -1);
#if SOFT_CYCLES_REDUCTION
                        uint8_t part_idx = part_to_shape[context_ptr->md_blk_arr_nsq[blk_index].part];
                        context_ptr->pred_depth_count[pred_depth_refinement + 2][part_idx]+= (blk_geom->bwidth*blk_geom->bheight);
#else
                        context_ptr->pred_depth_count[pred_depth_refinement + 2] += (blk_geom->bwidth*blk_geom->bheight);
#endif
                    }
                }
            }
        }
        blk_index += split_flag ?
            d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] :
            ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
}
#endif
#if ADAPTIVE_DEPTH_CR
/******************************************************
* Generate probabilities for the depth_cycles_reduction
******************************************************/
void generate_depth_prob(PictureControlSet * pcs_ptr, ModeDecisionContext *context_ptr)
{
    if (pcs_ptr->parent_pcs_ptr->slice_type != I_SLICE) {
#if SOFT_CYCLES_REDUCTION
        uint32_t pred_depth_count[DEPTH_DELTA_NUM][NUMBER_OF_SHAPES - 1] = { {0},{0},{0},{0},{0} };
#else
        uint32_t pred_depth_count[DEPTH_DELTA_NUM] = { 0 };
#endif
        uint32_t samples_num = 0;
        // Sum statistics from reference list0
#if ON_OFF_FEATURE_MRP
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->mrp_ctrls.ref_list0_count_try; ref_idx++) {
#else
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list0_count_try; ref_idx++) {
#endif
            EbReferenceObject *ref_obj_l0 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][ref_idx]->object_ptr;
#if SOFT_CYCLES_REDUCTION
            for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
                for (uint8_t part_idx = 0; part_idx < (NUMBER_OF_SHAPES - 1); part_idx++) {
                    pred_depth_count[pred_depth][part_idx] += ref_obj_l0->ref_pred_depth_count[pred_depth][part_idx];
                    samples_num += ref_obj_l0->ref_pred_depth_count[pred_depth][part_idx];
                }
            }
#else
            for (uint8_t pred_depth = 0; pred_depth < 5; pred_depth++) {
                pred_depth_count[pred_depth] += ref_obj_l0->ref_pred_depth_count[pred_depth];
                samples_num += ref_obj_l0->ref_pred_depth_count[pred_depth];
            }
#endif
        }
        // Sum statistics from reference list1
#if ON_OFF_FEATURE_MRP
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->mrp_ctrls.ref_list1_count_try; ref_idx++) {
#else
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list1_count_try; ref_idx++) {
#endif
            EbReferenceObject *ref_obj_l1 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx]->object_ptr;
#if SOFT_CYCLES_REDUCTION
            for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
                for (uint8_t part_idx = 0; part_idx < (NUMBER_OF_SHAPES - 1); part_idx++) {
                    pred_depth_count[pred_depth][part_idx] += ref_obj_l1->ref_pred_depth_count[pred_depth][part_idx];
                    samples_num += ref_obj_l1->ref_pred_depth_count[pred_depth][part_idx];
                }
            }
#else
            for (uint8_t pred_depth = 0; pred_depth < 5; pred_depth++) {
                pred_depth_count[pred_depth] += ref_obj_l1->ref_pred_depth_count[pred_depth];
                samples_num += ref_obj_l1->ref_pred_depth_count[pred_depth];
            }
#endif
        }
        // Generate the selection %
#if SOFT_CYCLES_REDUCTION
        uint32_t sum = 0;
#if !REMOVE_PRINT_STATEMENTS
        printf("\nstart \n");
#endif
        for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
            for (uint8_t part_idx = 0; part_idx < (NUMBER_OF_SHAPES - 1); part_idx++) {
                context_ptr->ad_md_prob[pred_depth][part_idx] = (uint32_t)((pred_depth_count[pred_depth][part_idx] * (uint32_t)DEPTH_PROB_PRECISION) / (uint32_t)samples_num);
                sum += context_ptr->ad_md_prob[pred_depth][part_idx];
#if !REMOVE_PRINT_STATEMENTS
                printf("%d\t", context_ptr->ad_md_prob[pred_depth][part_idx]);
#endif
            }
#if !REMOVE_PRINT_STATEMENTS
            printf("\n");
#endif
        }
#if !REMOVE_PRINT_STATEMENTS
        printf("\n");
        printf("\nsum = %d\n",sum/100);
#endif
        sum = 0;
        //Generate depth prob
        for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
            for (uint8_t part_idx = 1; part_idx < (NUMBER_OF_SHAPES - 1); part_idx++) {
                pred_depth_count[pred_depth][0] += pred_depth_count[pred_depth][part_idx];
            }
        }
        for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
            context_ptr->depth_prob[pred_depth] = (uint32_t)((pred_depth_count[pred_depth][0] * (uint32_t)100) / (uint32_t)samples_num);
            sum += context_ptr->depth_prob[pred_depth];
        }
#if !REMOVE_PRINT_STATEMENTS
        printf("\n");
        printf("\nsum2 = %d\n",sum);
#endif
#else
        for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++)
            context_ptr->depth_prob[pred_depth] = (uint32_t)((pred_depth_count[pred_depth] * (uint32_t)100) / (uint32_t)samples_num);
#endif

    }
#if SOFT_CYCLES_REDUCTION
    else {
        memcpy(context_ptr->ad_md_prob, intra_adaptive_md_cycles_reduction_th, sizeof(uint32_t) * DEPTH_DELTA_NUM * (NUMBER_OF_SHAPES - 1));
    }
#endif
}
#endif
#if ADAPTIVE_NSQ_CR
/******************************************************
* Generate probabilities for the nsq_cycles_reduction
******************************************************/
void generate_nsq_prob(PictureControlSet * pcs_ptr,ModeDecisionContext *context_ptr)
{
    if (pcs_ptr->parent_pcs_ptr->slice_type != I_SLICE) {
        uint32_t part_cnt[NUMBER_OF_SHAPES-1][FB_NUM][SSEG_NUM];
        uint8_t band, partidx, sse_idx;
        uint64_t samples_num = 0;
        // init stat
        memset(part_cnt, 0, sizeof(uint32_t) * (NUMBER_OF_SHAPES-1) * FB_NUM * SSEG_NUM);
        // Sum statistics from reference list0
#if ON_OFF_FEATURE_MRP
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->mrp_ctrls.ref_list0_count_try; ref_idx++) {
#else
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list0_count_try; ref_idx++) {
#endif
            EbReferenceObject *ref_obj_l0 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][ref_idx]->object_ptr;
            for (partidx = 0; partidx < NUMBER_OF_SHAPES-1; partidx++) {
                for (band = 0; band < FB_NUM; band++) {
                    for (sse_idx = 0; sse_idx < SSEG_NUM; sse_idx++) {
                        part_cnt[partidx][band][sse_idx] += ref_obj_l0->ref_part_cnt[partidx][band][sse_idx];
                        samples_num += ref_obj_l0->ref_part_cnt[partidx][band][sse_idx];
                    }
                }
            }
        }
        // Sum statistics from reference list1
#if ON_OFF_FEATURE_MRP
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->mrp_ctrls.ref_list1_count_try; ref_idx++) {
#else
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list1_count_try; ref_idx++) {
#endif
            EbReferenceObject *ref_obj_l1 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx]->object_ptr;
            for (partidx = 0; partidx < NUMBER_OF_SHAPES-1; partidx++) {
                for (band = 0; band < FB_NUM; band++) {
                    for (sse_idx = 0; sse_idx < SSEG_NUM; sse_idx++) {
                        part_cnt[partidx][band][sse_idx] += ref_obj_l1->ref_part_cnt[partidx][band][sse_idx];
                        samples_num += ref_obj_l1->ref_part_cnt[partidx][band][sse_idx];
                    }
                }
            }
        }
        for (partidx = 0; partidx < NUMBER_OF_SHAPES-1; partidx++) {
            for (band = 0; band < FB_NUM; band++) {
                for (sse_idx = 0; sse_idx < SSEG_NUM; sse_idx++) {
                    context_ptr->part_prob[partidx][band][sse_idx] = samples_num ? (uint16_t)((part_cnt[partidx][band][sse_idx] * 1000) / samples_num) :
                        block_prob_tab[0][partidx][band][sse_idx];
                }
            }
        }
    }
}
#endif
#if ADAPTIVE_NSQ_CR
void generate_statistics_nsq(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t blk_index = 0;
    uint32_t total_samples = 0;
    uint32_t count_non_zero_coeffs = 0;
    uint32_t part_cnt[NUMBER_OF_SHAPES-1][FB_NUM][SSEG_NUM];
    uint8_t band,partidx,sse_idx;
    for (partidx = 0; partidx < NUMBER_OF_SHAPES-1; partidx++) {
        for (band = 0; band < FB_NUM; band++) {
            memset(part_cnt[partidx][band], 0, sizeof(uint32_t) * SSEG_NUM);
        }
    }
    EbBool split_flag;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
        if (scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    if (context_ptr->md_local_blk_unit[blk_index].avail_blk_flag) {
                        uint8_t band_idx = 0;
                        uint8_t sq_size_idx = 7 - (uint8_t)Log2f((uint8_t)blk_geom->sq_size);
                        uint64_t band_width = (sq_size_idx == 0) ? 100 : (sq_size_idx == 1) ? 50 : 20;
                        uint8_t part_idx = part_to_shape[context_ptr->md_blk_arr_nsq[blk_index].part];
                        uint8_t sse_g_band = context_ptr->md_local_blk_unit[blk_geom->sqi_mds].avail_blk_flag ?
                            context_ptr->md_local_blk_unit[blk_geom->sqi_mds].sse_gradian_band[part_idx] : 1;
                        count_non_zero_coeffs = context_ptr->md_local_blk_unit[blk_index].count_non_zero_coeffs;
                        total_samples = (blk_geom->bwidth*blk_geom->bheight);
                        if (count_non_zero_coeffs >= ((total_samples * 18) / band_width)) {
                            band_idx = 9;
                        }
                        else if (count_non_zero_coeffs >= ((total_samples * 16) / band_width)) {
                            band_idx = 8;
                        }
                        else if (count_non_zero_coeffs >= ((total_samples * 14) / band_width)) {
                            band_idx = 7;
                        }
                        else if (count_non_zero_coeffs >= ((total_samples * 12) / band_width)) {
                            band_idx = 6;
                        }
                        else if (count_non_zero_coeffs >= ((total_samples * 10) / band_width)) {
                            band_idx = 5;
                        }
                        else if (count_non_zero_coeffs >= ((total_samples * 8) / band_width)) {
                            band_idx = 4;
                        }
                        else if (count_non_zero_coeffs >= ((total_samples * 6) / band_width)) {
                            band_idx = 3;
                        }
                        else if (count_non_zero_coeffs >= ((total_samples * 4) / band_width)) {
                            band_idx = 2;
                        }
                        else if (count_non_zero_coeffs >= ((total_samples * 2) / band_width)) {
                            band_idx = 1;
                        }
                        else {
                            band_idx = 0;
                        }
                        if (sq_size_idx == 0)
                            band_idx = band_idx == 0 ? 0 : band_idx <= 2 ? 1 : 2;
                        else if (sq_size_idx == 1)
                            band_idx = band_idx == 0 ? 0 : band_idx <= 3 ? 1 : 2;
                        else
                            band_idx = band_idx == 0 ? 0 : band_idx <= 8 ? 1 : 2;

                        part_cnt[part_to_shape[context_ptr->md_blk_arr_nsq[blk_index].part]][band_idx][sse_g_band] += (blk_geom->bwidth*blk_geom->bheight);

                    }
                }
            }
        }
        blk_index += split_flag ?
            d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] :
            ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
    for (partidx = 0; partidx < NUMBER_OF_SHAPES-1; partidx++) {
        for (band = 0; band < FB_NUM; band++) {
            for (sse_idx = 0; sse_idx < SSEG_NUM; sse_idx++) {
                context_ptr->part_cnt[partidx][band][sse_idx] += part_cnt[partidx][band][sse_idx];
            }
        }
    }
}
#endif
#if ADAPTIVE_TXT_CR

/******************************************************
* Generate probabilities for the txt_cycles_reduction
******************************************************/
void generate_txt_prob(PictureControlSet * pcs_ptr,ModeDecisionContext *context_ptr)
{
    if (pcs_ptr->parent_pcs_ptr->slice_type != I_SLICE) {
        uint32_t txt_cnt[TXT_DEPTH_DELTA_NUM][TX_TYPES] = { {0},{0} };
        uint32_t samples_num = 0;
        // Sum statistics from reference list0
#if ON_OFF_FEATURE_MRP
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->mrp_ctrls.ref_list0_count_try; ref_idx++) {
#else
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list0_count_try; ref_idx++) {
#endif
            EbReferenceObject *ref_obj_l0 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][ref_idx]->object_ptr;
            for (uint8_t depth_delta = 0; depth_delta < TXT_DEPTH_DELTA_NUM; depth_delta++) {
                for (uint8_t txs_idx = 0; txs_idx < TX_TYPES; txs_idx++) {
                    txt_cnt[depth_delta][txs_idx] += ref_obj_l0->ref_txt_cnt[depth_delta][txs_idx];
                    samples_num += ref_obj_l0->ref_txt_cnt[depth_delta][txs_idx];
                }
            }
        }
        // Sum statistics from reference list1
#if ON_OFF_FEATURE_MRP
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->mrp_ctrls.ref_list1_count_try; ref_idx++) {
#else
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list1_count_try; ref_idx++) {
#endif
            EbReferenceObject *ref_obj_l1 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx]->object_ptr;
            for (uint8_t depth_delta = 0; depth_delta < TXT_DEPTH_DELTA_NUM; depth_delta++) {
                for (uint8_t txs_idx = 1; txs_idx < TX_TYPES; txs_idx++) {
                    txt_cnt[depth_delta][txs_idx] += ref_obj_l1->ref_txt_cnt[depth_delta][txs_idx];
                    samples_num += ref_obj_l1->ref_txt_cnt[depth_delta][txs_idx];
                }
            }
        }
        for (uint8_t depth_delta = 0; depth_delta < TXT_DEPTH_DELTA_NUM; depth_delta++) {
            for (uint8_t txs_idx = 1; txs_idx < TX_TYPES; txs_idx++) {
                context_ptr->txt_prob[depth_delta][txs_idx] = (uint32_t)((txt_cnt[depth_delta][txs_idx] * (uint32_t)10000) / (uint32_t)samples_num);
            }
        }
    }
}
#endif
#if SB_CLASSIFIER
#if CLEANUP_CYCLE_ALLOCATION
const uint32_t sb_class_th[NUMBER_OF_SB_CLASS] = { 0,85,75,65,60,55,50,45,40,
                                                   35,30,25,20,17,14,10,6,3,0 };
#endif
static uint8_t determine_sb_class(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t blk_index = 0;
    uint64_t total_samples = 0;
    uint64_t count_non_zero_coeffs = 0;
    uint8_t sb_class = NONE_CLASS;
#if !CLEANUP_CYCLE_ALLOCATION
    SbClassControls *sb_class_ctrls = &context_ptr->sb_class_ctrls;
#endif
    EbBool split_flag;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
        if (scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
#if SB_CLASSIFIER_R2R_FIX
            context_ptr->md_local_blk_unit[blk_index].avail_blk_flag &&
#endif
            is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    count_non_zero_coeffs += context_ptr->md_local_blk_unit[blk_index].count_non_zero_coeffs;
                    total_samples += (blk_geom->bwidth*blk_geom->bheight);
                }
            }
        }
        blk_index += split_flag ?
            d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] :
            ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
#if CLEANUP_CYCLE_ALLOCATION
    if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_1]) / 100))
        sb_class = SB_CLASS_1;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_2]) / 100))
        sb_class = SB_CLASS_2;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_3]) / 100))
        sb_class = SB_CLASS_3;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_4]) / 100))
        sb_class = SB_CLASS_4;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_5]) / 100))
        sb_class = SB_CLASS_5;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_6]) / 100))
        sb_class = SB_CLASS_6;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_7]) / 100))
        sb_class = SB_CLASS_7;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_8]) / 100))
        sb_class = SB_CLASS_8;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_9]) / 100))
        sb_class = SB_CLASS_9;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_10]) / 100))
        sb_class = SB_CLASS_10;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_11]) / 100))
        sb_class = SB_CLASS_11;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_12]) / 100))
        sb_class = SB_CLASS_12;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_13]) / 100))
        sb_class = SB_CLASS_13;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_14]) / 100))
        sb_class = SB_CLASS_14;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_15]) / 100))
        sb_class = SB_CLASS_15;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_16]) / 100))
        sb_class = SB_CLASS_16;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_17]) / 100))
        sb_class = SB_CLASS_17;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_th[SB_CLASS_18]) / 100))
        sb_class = SB_CLASS_18;
#else
#if MULTI_BAND_ACTIONS
#if NON_UNIFORM_NSQ_BANDING
    if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_1]) / 100))
        sb_class = SB_CLASS_1;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_2]) / 100))
        sb_class = SB_CLASS_2;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_3]) / 100))
        sb_class = SB_CLASS_3;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_4]) / 100))
        sb_class = SB_CLASS_4;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_5]) / 100))
        sb_class = SB_CLASS_5;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_6]) / 100))
        sb_class = SB_CLASS_6;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_7]) / 100))
        sb_class = SB_CLASS_7;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_8]) / 100))
        sb_class = SB_CLASS_8;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_9]) / 100))
        sb_class = SB_CLASS_9;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_10]) / 100))
        sb_class = SB_CLASS_10;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_11]) / 100))
        sb_class = SB_CLASS_11;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_12]) / 100))
        sb_class = SB_CLASS_12;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_13]) / 100))
        sb_class = SB_CLASS_13;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_14]) / 100))
        sb_class = SB_CLASS_14;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_15]) / 100))
        sb_class = SB_CLASS_15;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_16]) / 100))
        sb_class = SB_CLASS_16;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_17]) / 100))
        sb_class = SB_CLASS_17;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_18]) / 100))
        sb_class = SB_CLASS_18;
#else
    if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_1]) / 20))
        sb_class = SB_CLASS_1;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_2]) / 20))
        sb_class = SB_CLASS_2;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_3]) / 20))
        sb_class = SB_CLASS_3;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_4]) / 20))
        sb_class = SB_CLASS_4;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_5]) / 20))
        sb_class = SB_CLASS_5;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_6]) / 20))
        sb_class = SB_CLASS_6;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_7]) / 20))
        sb_class = SB_CLASS_7;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_8]) / 20))
        sb_class = SB_CLASS_8;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_9]) / 20))
        sb_class = SB_CLASS_9;
     else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_10]) / 20))
        sb_class = SB_CLASS_10;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_11]) / 20))
        sb_class = SB_CLASS_11;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_12]) / 20))
        sb_class = SB_CLASS_12;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_13]) / 20))
        sb_class = SB_CLASS_13;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_14]) / 20))
        sb_class = SB_CLASS_14;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_15]) / 20))
        sb_class = SB_CLASS_15;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_16]) / 20))
        sb_class = SB_CLASS_16;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_17]) / 20))
        sb_class = SB_CLASS_17;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_18]) / 20))
        sb_class = SB_CLASS_18;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_19]) / 20))
        sb_class = SB_CLASS_19;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[SB_CLASS_20]) / 20))
        sb_class = SB_CLASS_20;
#endif
#else
    if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[HIGH_COMPLEX_CLASS]) / 20))
        sb_class = HIGH_COMPLEX_CLASS;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[MEDIUM_COMPLEX_CLASS]) / 20))
        sb_class = MEDIUM_COMPLEX_CLASS;
    else if (count_non_zero_coeffs >= ((total_samples * sb_class_ctrls->sb_class_th[LOW_COMPLEX_CLASS]) / 20))
        sb_class = LOW_COMPLEX_CLASS;
#if NEW_CYCLE_ALLOCATION
    else if (count_non_zero_coeffs == ((total_samples * sb_class_ctrls->sb_class_th[VERY_LOW_COMPLEX_CLASS]) / 20))
        sb_class = VERY_LOW_COMPLEX_CLASS;
#endif
#endif
#endif
    return sb_class;
}
#endif
#if DEPTH_CYCLES_REDUCTION
#define DEPTH_MAX_PROB 300 // max probabilty value for depth 100 -> 10%
// Depth probabilies per sq_size, pedicted depth and frequency band
// for sc content
#if !ADAPTIVE_DEPTH_CR
const uint16_t depth_cycles_reduction_sc_th[6][5][4] = {
{
{ 0,0,0,0},
{ 0,0,0,0},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 0,0,0,0},
{ 0,0,0,0}
},
{
{ 0,0,0,0},
{ 0,4,14,0},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 0,2,6,0},
{ 1,4,6,0}
},
{
{ 0,0,0,0},
{ 0,2,13,0},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 17,13,20,0},
{ 26,5,3,0}
},
{
{ 0,0,0,0},
{ 10,17,62,0},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 79,44,44,0},
{ 45,12,11,0}
},
{
{ 0,0,0,0},
{ 72,63,112,0},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 79,29,22,0},
{ 0,0,0,0}
},
{
{ 0,0,0,0},
{ 78,47,38,0},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 0,0,0,0},
{ 0,0,0,0}
}
};
// Depth probabilies per sq_size, pedicted depth and frequency band
// for non-sc content
uint16_t depth_cycles_reduction_th[6][5][4] = {
{
{ 0,0,0,0},
{ 0,0,0,0},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 0,3,6,30},
{ 0,1,1,1},
},
{
{ 0,0,0,0},
{ 1,11,26,234},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 7,19,10,12},
{ 6,3,1,0},
},
{
{ 0,0,0,0},
{ 11,44,82,108},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 29,14,5,2},
{ 7,0,0,0},
},
{
{ 0,0,0,0},
{ 24,90,73,32},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 28,6,1,0},
{ 2,0,0,0},
},
{
{ 0,0,0,0},
{ 26,24,8,1},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 4,1,0,0},
{ 0,0,0,0},
},
{
{ 0,0,0,0},
{ 4,1,0,0},
{DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB,DEPTH_MAX_PROB},
{ 0,0,0,0},
{ 0,0,0,0}
}
};
#endif
#endif
static void perform_pred_depth_refinement(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                          ModeDecisionContext *context_ptr, uint32_t sb_index) {
#if DEPTH_PART_CLEAN_UP
    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
#else
    MdcSbData *results_ptr = &pcs_ptr->mdc_sb_array[sb_index];
#endif
    uint32_t   blk_index   = 0;

    // Reset mdc_sb_array data to defaults; it will be updated based on the predicted blocks (stored in md_blk_arr_nsq)
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom                              = get_blk_geom_mds(blk_index);
        results_ptr->leaf_data_array[blk_index].consider_block = 0;
        results_ptr->leaf_data_array[blk_index].split_flag =
            blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        results_ptr->leaf_data_array[blk_index].refined_split_flag =
            blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
#if TRACK_PER_DEPTH_DELTA
        results_ptr->leaf_data_array[blk_index].pred_depth_refinement = -8;
#endif
#if ADAPTIVE_DEPTH_CR
        results_ptr->leaf_data_array[blk_index].pred_depth = -8;
#endif
        blk_index++;
    }

    results_ptr->leaf_count = 0;
    blk_index               = 0;

    SuperBlock *sb_ptr = pcs_ptr->sb_ptr_array[sb_index];

    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
        const unsigned   tot_d1_blocks = blk_geom->sq_size == 128
            ? 17
            : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;

        // if the parent square is inside inject this block
        uint8_t is_blk_allowed =
            pcs_ptr->slice_type != I_SLICE ? 1 : (blk_geom->sq_size < 128) ? 1 : 0;

        // derive split_flag
        EbBool split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;

        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    int8_t s_depth = 0;
                    int8_t e_depth = 0;

                    if (context_ptr->pd_pass == PD_PASS_0) {
                        uint32_t full_lambda =  context_ptr->hbd_mode_decision ?
                            context_ptr->full_lambda_md[EB_10_BIT_MD] :
                            context_ptr->full_lambda_md[EB_8_BIT_MD];

                        uint32_t sb_width = scs_ptr->seq_header.sb_size == BLOCK_128X128 ?
                            128 : 64;
                        uint32_t sb_height = scs_ptr->seq_header.sb_size == BLOCK_128X128 ?
                            128 : 64;
                        uint64_t dist_sum = (sb_width * sb_height * 100);

                        uint64_t early_exit_th = RDCOST(full_lambda, 16, dist_sum);
                        uint64_t best_part_cost = generate_best_part_cost(
                            scs_ptr,
                            pcs_ptr,
                            context_ptr,
                            sb_index);

#if FIX_MR_PD1
#if MR_MODE_FOR_PIC_MULTI_PASS_PD_MODE_1
#if MAR19_ADOPTIONS
                        // Shut thresholds in MR_MODE
#if APR22_ADOPTIONS
#if MAY12_ADOPTIONS
#if JUNE15_ADOPTIONS
#if REMOVE_MR_MACRO
                        if (pcs_ptr->enc_mode <= ENC_MRS) {
#else
                        if (MRS_MODE) {
#endif
                            s_depth = -2;
                            e_depth = 2;
                        }
#if REMOVE_MR_MACRO
                        else if (pcs_ptr->enc_mode <= ENC_MR) {
#else
                        else if (MR_MODE) {
#endif
#else
                        if (MR_MODE_MULTI_PASS_PD) {
#endif
#else
                        if (MR_MODE_MULTI_PASS_PD || (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M0)) {
#endif
#if MAY19_ADOPTIONS
#if MR_DEPTH_REFINEMENT
#if JUNE8_ADOPTIONS
#if UNIFY_SC_NSC
                            if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_240p_RANGE) {
#else
                            if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
                                s_depth = -2;
                                e_depth = 2;
                            }
                            else if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_240p_RANGE) {
#endif
#else
                            if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_240p_RANGE ||
                                pcs_ptr->parent_pcs_ptr->sc_content_detected) {
#endif
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : -1;
                                e_depth = 2;
                            }
                            else if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) {
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : -1;
                                e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : 1;
                            }
                            else {
                                s_depth = -2;
                                e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : 1;
                            }
#else
                            if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
                                s_depth = -2;
                                e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : 1;
                            }
                            else {
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : -1;
                                e_depth = 2;
                            }
#endif
#else
                            s_depth = -2;
                            e_depth = 2;
#endif
                        }
#else
                        if (MR_MODE_MULTI_PASS_PD) {
                            s_depth = -3;
                            e_depth = 3;
                        }
#endif
#if ADOPT_SKIPPING_PD1
                        else if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0) {
#if M8_MPPD
#if !MAY17_ADOPTIONS
#if MAY12_ADOPTIONS
#if MAY16_7PM_ADOPTIONS
                            if (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M0) {
#else
                            if (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M2) {
#endif
                                s_depth = -2;
                                e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : 1;
                            }
                            else
#endif
#endif
#if MAY16_7PM_ADOPTIONS
#if JUNE26_ADOPTIONS
                            if (pcs_ptr->enc_mode <= ENC_M6) {
#else
#if JUNE17_ADOPTIONS
                            if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if JUNE11_ADOPTIONS
                            if (pcs_ptr->enc_mode <= ENC_M3) {
#else
#if JUNE8_ADOPTIONS
                            if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M2)) {
#else
#if PRESET_SHIFITNG
                            if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M1)) {
#else
                            if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M2)) {
#endif
#endif
#endif
#endif
#endif
#else
#if M1_C3_ADOPTIONS
                            if (pcs_ptr->enc_mode <= ENC_M0) {
#else
#if MAY12_ADOPTIONS
                            if (pcs_ptr->enc_mode <= ENC_M2 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M4)) {
#else
#if APR24_ADOPTIONS_M6_M7
                            if (pcs_ptr->enc_mode <= ENC_M6) {
#else
                            if(pcs_ptr->enc_mode <= ENC_M5) {
#endif
#endif
#endif
#endif
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : -1;
                                e_depth = pcs_ptr->slice_type == I_SLICE ?  2 :  1;
                            }
#if !JUNE11_ADOPTIONS
#if M1_C3_ADOPTIONS
#if PRESET_SHIFITNG
#if PRESET_SHIFITNG
                            else if (pcs_ptr->enc_mode <= ENC_M1 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M2)) {
#else
                            else if (pcs_ptr->enc_mode <= ENC_M1 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M4)) {
#endif
#else
                            else if (pcs_ptr->enc_mode <= ENC_M2 || (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M4)) {
#endif
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : -1;
#if MAY16_7PM_ADOPTIONS
                                e_depth = pcs_ptr->slice_type == I_SLICE ? 2 :
                                    (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag && pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                                    ? 1
                                    : 0;
#else
                                e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 0;
#endif
                            }
#endif
#endif
#if !JUNE17_ADOPTIONS
#if MAY12_ADOPTIONS
#if PRESET_SHIFITNG
                            else if (pcs_ptr->enc_mode <= ENC_M4) {
#else
                            else if (pcs_ptr->enc_mode <= ENC_M6) {
#endif
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : -1;
                                e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : 0;
                            }
#endif
#endif
                            else {
#if M5_I_PD
#if UPGRADE_M6_M7_M8
#if JUNE26_ADOPTIONS
                                if (pcs_ptr->enc_mode <= ENC_M8) {
                                    s_depth = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? -1 : 0;
                                    e_depth = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 0;
                                }
#else
#if PRESET_SHIFITNG
                                if (pcs_ptr->enc_mode <= ENC_M5) {
#else
                                if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
                                    s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? -1 : 0;
                                    e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 0;
                                }
#if JUNE17_ADOPTIONS
                                else if (pcs_ptr->enc_mode <= ENC_M6) {
                                    s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? -1: 0;
                                    e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? 1 : 0;
                                }
#endif
#endif
                                else {
#if REVERT_WHITE // MPPD
                                    s_depth = pcs_ptr->slice_type == I_SLICE ? -1 : 0;
                                    e_depth = pcs_ptr->slice_type == I_SLICE ? 1 : 0;
#else
                                    s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : 0;
                                    e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : 0;
#endif
                                }
#else
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : 0;
                                e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : 0;
#endif
#else
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -1 : 0;
                                e_depth = pcs_ptr->slice_type == I_SLICE ?  1 : 0;
#endif
                            }
#else
                            s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : -1;
                            e_depth = pcs_ptr->slice_type == I_SLICE ? 2 : 1;
#endif
                        }
#endif
#if !UNIFY_SC_NSC
#if MAR30_ADOPTIONS
                        else if ((pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M1)) {

#if OPT_BLOCK_INDICES_GEN_1
                            s_depth = -3;
                            e_depth =  3;
#else
                            s_depth = (blk_geom->sq_size == 64 && pcs_ptr->parent_pcs_ptr->sb_64x64_simulated) ? 0
                                    : (blk_geom->sq_size == 32 && pcs_ptr->parent_pcs_ptr->sb_64x64_simulated) ? -1
                                    : (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->sb_64x64_simulated) ? -2
                                    : -3;
                            e_depth = (blk_geom->sq_size == 8 && pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 0
                                    : (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 1
                                    : (blk_geom->sq_size == 32 && pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 2
                                    : 3;
#endif
                        }
#endif
#endif
                        else if (best_part_cost < early_exit_th && pcs_ptr->parent_pcs_ptr->multi_pass_pd_level != MULTI_PASS_PD_LEVEL_0) {
#else
                        if (best_part_cost < early_exit_th && pcs_ptr->parent_pcs_ptr->multi_pass_pd_level != MULTI_PASS_PD_LEVEL_0 && !MR_MODE_MULTI_PASS_PD) {
#endif
#else
                        if (best_part_cost < early_exit_th && pcs_ptr->parent_pcs_ptr->multi_pass_pd_level != MULTI_PASS_PD_LEVEL_0) {
#endif
#else
                        if (best_part_cost < early_exit_th) {
#endif
                            s_depth = 0;
                            e_depth = 0;
                        }
                        else {
                        derive_start_end_depth(pcs_ptr,
                                               sb_ptr,
                                               scs_ptr->seq_header.sb_size,
                                               &s_depth,
                                               &e_depth,
                                               blk_geom);
                        }
#if DEPTH_CYCLES_REDUCTION
 #if ADAPTIVE_DEPTH_CR
                        DepthCycleRControls*depth_cycle_red_ctrls = &context_ptr->depth_cycles_red_ctrls;
                        if (depth_cycle_red_ctrls->enabled) {
                            int8_t addj_s_depth = 0;
                            int8_t addj_e_depth = 0;
                            if (context_ptr->sb_class) {

                                if (depth_cycle_red_ctrls->th) {
                                    addj_s_depth = context_ptr->depth_prob[0] < depth_cycle_red_ctrls->th ? 0 : -2;
                                    if (addj_s_depth == 0)
                                        addj_s_depth = context_ptr->depth_prob[1] < depth_cycle_red_ctrls->th ? 0 : -1;
                                    addj_e_depth = context_ptr->depth_prob[4] < depth_cycle_red_ctrls->th ? 0 : 2;
                                    if (addj_e_depth == 0)
                                        addj_e_depth = context_ptr->depth_prob[3] < depth_cycle_red_ctrls->th ? 0 : 1;
                                }
                                s_depth = MAX(s_depth, addj_s_depth);
                                e_depth = MIN(e_depth, addj_e_depth);
                            }
                        }
#else
                        DepthCycleRControls*depth_cycle_red_ctrls = &context_ptr->depth_cycles_red_ctrls;
                        uint8_t sq_size_idx = 7 - (uint8_t)Log2f((uint8_t)context_ptr->blk_geom->sq_size);
                        if (depth_cycle_red_ctrls->enabled) {
                            int8_t addj_s_depth = 0;
                            int8_t addj_e_depth = 0;
                            if (context_ptr->sb_class) {
                                uint8_t frequency_band = context_ptr->sb_class <= 11 ? 0 : context_ptr->sb_class <= 18 ? 1 : context_ptr->sb_class <= 23 ? 2 : 3;
                                if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
                                    if (depth_cycle_red_ctrls->th) {
                                        addj_s_depth = depth_cycles_reduction_sc_th[sq_size_idx][0][frequency_band] < depth_cycle_red_ctrls->th ? 0 : -2;
                                        if (addj_s_depth == 0)
                                            addj_s_depth = depth_cycles_reduction_sc_th[sq_size_idx][1][frequency_band] < depth_cycle_red_ctrls->th ? 0 : -1;
                                        addj_e_depth = depth_cycles_reduction_sc_th[sq_size_idx][4][frequency_band] < depth_cycle_red_ctrls->th ? 0 : 2;
                                        if (addj_e_depth == 0)
                                            addj_e_depth = depth_cycles_reduction_sc_th[sq_size_idx][3][frequency_band] < depth_cycle_red_ctrls->th ? 0 : 1;
                                    }
                                }else{
                                    if (depth_cycle_red_ctrls->th) {
                                        addj_s_depth = depth_cycles_reduction_th[sq_size_idx][0][frequency_band] < depth_cycle_red_ctrls->th ? 0 : -2;
                                        if (addj_s_depth == 0)
                                            addj_s_depth = depth_cycles_reduction_th[sq_size_idx][1][frequency_band] < depth_cycle_red_ctrls->th ? 0 : -1;
                                        addj_e_depth = depth_cycles_reduction_th[sq_size_idx][4][frequency_band] < depth_cycle_red_ctrls->th ? 0 : 2;
                                        if (addj_e_depth == 0)
                                            addj_e_depth = depth_cycles_reduction_th[sq_size_idx][3][frequency_band] < depth_cycle_red_ctrls->th ? 0 : 1;
                                    }
                                }
                            }
                            s_depth = MAX(s_depth, addj_s_depth);
                            e_depth = MIN(e_depth, addj_e_depth);
                        }
#endif
#endif
                    } else if (context_ptr->pd_pass == PD_PASS_1) {
#if DEPTH_PART_CLEAN_UP
                        EbBool zero_coeff_present_flag =
                            context_ptr->md_blk_arr_nsq[blk_index].block_has_coeff == 0;
#if ADD_NEW_MPPD_LEVEL
#if MAR23_ADOPTIONS
                        if (pcs_ptr->slice_type == I_SLICE) {
#if MAR25_ADOPTIONS
#if OPT_BLOCK_INDICES_GEN_1
                            s_depth =  -1;
#if UNIFY_SC_NSC
#if REMOVE_MR_MACRO
                            e_depth = (pcs_ptr->parent_pcs_ptr->enc_mode <= ENC_MR) ? 3 : 2;
#else
                            e_depth = (MR_MODE_MULTI_PASS_PD) ? 3 : 2;
#endif
#else
                            e_depth =  (MR_MODE_MULTI_PASS_PD || pcs_ptr->parent_pcs_ptr->sc_content_detected) ? 3 : 2;
#endif
#else
                            s_depth = (blk_geom->sq_size == 64 && pcs_ptr->parent_pcs_ptr->sb_64x64_simulated) ? 0 : -1;
                            e_depth = (blk_geom->sq_size == 8  && pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 0
                                    : (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 1
                                    : (blk_geom->sq_size == 32 && pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 2
                                    : (MR_MODE_MULTI_PASS_PD || pcs_ptr->parent_pcs_ptr->sc_content_detected) ? 3 : 2;
#endif
#else
                            s_depth = -1;
                            e_depth = 2;
#endif
                        }
                        else if (zero_coeff_present_flag && (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 || pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4)) {
#else
                        if (zero_coeff_present_flag && (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 || pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4)) {
#endif
#else
                        if (zero_coeff_present_flag && (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 || pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3)) {
#endif
#else
                        EbBool zero_coeff_present_flag = EB_FALSE;

                        if (pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2)
                            zero_coeff_present_flag =
                                context_ptr->md_blk_arr_nsq[blk_index].block_has_coeff == 0;

                        else if (pcs_ptr->parent_pcs_ptr->pic_depth_mode ==
                                 PIC_MULTI_PASS_PD_MODE_3) {
                            switch (blk_geom->bsize) {
                            case BLOCK_128X128:
                                zero_coeff_present_flag =
                                    (context_ptr->md_local_blk_unit[blk_index].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index].block_has_coeff ==
                                         0) || // SQ
                                    (context_ptr->md_local_blk_unit[blk_index + 1].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 1].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 2].block_has_coeff ==
                                         0) || // H
                                    (context_ptr->md_local_blk_unit[blk_index + 3].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 3].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 4].block_has_coeff ==
                                         0) || // V
                                    (context_ptr->md_local_blk_unit[blk_index + 5].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 5].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 6].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 7].block_has_coeff ==
                                         0) ||
                                    (context_ptr->md_local_blk_unit[blk_index + 8].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 8].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 9].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 10].block_has_coeff ==
                                         0) ||
                                    (context_ptr->md_local_blk_unit[blk_index + 11]
                                         .avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 11].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 12].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 13].block_has_coeff ==
                                         0) ||
                                    (context_ptr->md_local_blk_unit[blk_index + 14]
                                         .avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 14].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 15].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 16].block_has_coeff ==
                                         0);
                                break;

                            case BLOCK_64X64:
                            case BLOCK_32X32:
                            case BLOCK_16X16:
                                zero_coeff_present_flag =
                                    (context_ptr->md_local_blk_unit[blk_index].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index].block_has_coeff ==
                                         0) || // SQ
                                    (context_ptr->md_local_blk_unit[blk_index + 1].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 1].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 2].block_has_coeff ==
                                         0) || // H
                                    (context_ptr->md_local_blk_unit[blk_index + 3].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 3].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 4].block_has_coeff ==
                                         0) || // V
                                    (context_ptr->md_local_blk_unit[blk_index + 5].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 5].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 6].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 7].block_has_coeff ==
                                         0) ||
                                    (context_ptr->md_local_blk_unit[blk_index + 8].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 8].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 9].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 10].block_has_coeff ==
                                         0) ||
                                    (context_ptr->md_local_blk_unit[blk_index + 11]
                                         .avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 11].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 12].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 13].block_has_coeff ==
                                         0) ||
                                    (context_ptr->md_local_blk_unit[blk_index + 14]
                                         .avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 14].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 15].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 16].block_has_coeff ==
                                         0) ||
                                    (context_ptr->md_local_blk_unit[blk_index + 17]
                                         .avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 17].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 18].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 19].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 20].block_has_coeff ==
                                         0) ||
                                    (context_ptr->md_local_blk_unit[blk_index + 21]
                                         .avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 21].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 22].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 23].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 24].block_has_coeff ==
                                         0);
                                break;

                            case BLOCK_8X8:
                                zero_coeff_present_flag =
                                    (context_ptr->md_local_blk_unit[blk_index].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index].block_has_coeff ==
                                         0) || // SQ
                                    (context_ptr->md_local_blk_unit[blk_index + 1].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 1].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 2].block_has_coeff ==
                                         0) || // H
                                    (context_ptr->md_local_blk_unit[blk_index + 3].avail_blk_flag &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 3].block_has_coeff ==
                                         0 &&
                                     context_ptr->md_blk_arr_nsq[blk_index + 4].block_has_coeff ==
                                         0); // V
                                break;

                            case BLOCK_4X4:
                                zero_coeff_present_flag =
                                    (context_ptr->md_blk_arr_nsq[blk_index].block_has_coeff ==
                                     0); // SQ
                                break;

                            default: assert(0); break;
                            }
                        }

                        if (zero_coeff_present_flag) {
#endif
                            s_depth = 0;
                            e_depth = 0;
                        } else

#if MR_MODE_FOR_PIC_MULTI_PASS_PD_MODE_1 || ADD_NEW_MPPD_LEVEL
#if ADD_NEW_MPPD_LEVEL
                            // This removes the SQ-versus-NSQ decision for the new MULTI_PASS_PD_LEVEL_1
#if REMOVE_MR_MACRO
                            if (pcs_ptr->parent_pcs_ptr->enc_mode <= ENC_MR || pcs_ptr->parent_pcs_ptr->multi_pass_pd_level <= MULTI_PASS_PD_LEVEL_1) {
#else
                            if (MR_MODE_MULTI_PASS_PD || pcs_ptr->parent_pcs_ptr->multi_pass_pd_level <= MULTI_PASS_PD_LEVEL_1) { // Active when multi_pass_pd_level = PIC_MULTI_PASS_PD_MODE_2 or PIC_MULTI_PASS_PD_MODE_3 or PIC_MULTI_PASS_PD_MODE_4
#endif
#else
                            if (MR_MODE_MULTI_PASS_PD) { // Active when multi_pass_pd_level = PIC_MULTI_PASS_PD_MODE_1 or PIC_MULTI_PASS_PD_MODE_2 or PIC_MULTI_PASS_PD_MODE_3
#endif
#if OPT_BLOCK_INDICES_GEN_1
                                s_depth = -1;
                                e_depth =  1;
#else
                                s_depth = (blk_geom->sq_size == 64 && pcs_ptr->parent_pcs_ptr->sb_64x64_simulated) ? 0 : -1;
                                e_depth = (blk_geom->sq_size == 8 && pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 0 : 1;
#endif
                            }
                            else
#endif
                            if (context_ptr->md_local_blk_unit[blk_index].best_d1_blk == blk_index) {

#if DEPTH_PART_CLEAN_UP // refinement rest
#if REMOVE_UNUSED_CODE_PH2
                            s_depth = -1;
#else
                            s_depth = (blk_geom->sq_size == 64 && pcs_ptr->parent_pcs_ptr->sb_64x64_simulated) ? 0 : -1;
#endif
#else
                                s_depth = -1;
#endif
                                e_depth = 0;
                            } else {
                                s_depth = 0;
#if DEPTH_PART_CLEAN_UP // refinement rest
#if OPT_BLOCK_INDICES_GEN_1
                            e_depth = 1;
#else
                            e_depth = (blk_geom->sq_size == 8 && pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 0 : 1;
#endif
#else
                                e_depth = 1;
#endif
                            }
                    }

#if NEW_CYCLE_ALLOCATION && !DISALLOW_ALL_ACTIONS
                    if (!pcs_ptr->parent_pcs_ptr->sc_content_detected) {
                        s_depth = (context_ptr->sb_class == HIGH_COMPLEX_CLASS || context_ptr->sb_class == MEDIUM_COMPLEX_CLASS) ? 0 : s_depth;
                        e_depth = (context_ptr->sb_class == HIGH_COMPLEX_CLASS || context_ptr->sb_class == MEDIUM_COMPLEX_CLASS) ? 0 : e_depth;
                    }
#endif

#if ADOPT_SKIPPING_PD1
                    // Check that the start and end depth are in allowed range, given other features
                    // which restrict allowable depths
#if !REMOVE_UNUSED_CODE_PH2
                    if (pcs_ptr->parent_pcs_ptr->sb_64x64_simulated) {
                        s_depth = (blk_geom->sq_size == 64) ? 0
                                : (blk_geom->sq_size == 32) ? MAX(-1, s_depth)
                                : (blk_geom->sq_size == 16) ? MAX(-2, s_depth)
                                : s_depth;
                    }
#endif
#if OPT_BLOCK_INDICES_GEN_1
                    if (context_ptr->disallow_4x4) {
#else
                    if (pcs_ptr->parent_pcs_ptr->disallow_4x4) {
#endif
                        e_depth = (blk_geom->sq_size == 8) ? 0
                                : (blk_geom->sq_size == 16) ? MIN(1, e_depth)
                                : (blk_geom->sq_size == 32) ? MIN(2, e_depth)
                                : e_depth;
                    }
#endif
                    // Add current pred depth block(s)
                    for (unsigned block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].consider_block = 1;
#if TRACK_PER_DEPTH_DELTA
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].pred_depth_refinement = 0;
#endif
#if ADAPTIVE_DEPTH_CR
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].pred_depth = (int8_t)blk_geom->depth;
#endif
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].refined_split_flag =
                            EB_FALSE;
                    }

#if ADAPTIVE_DEPTH_CR
                    uint8_t sq_size_idx = 7 - (uint8_t)Log2f((uint8_t)blk_geom->sq_size);
#endif
                    // Add block indices of upper depth(s)
                    if (s_depth != 0)
#if TRACK_PER_DEPTH_DELTA
#if ADAPTIVE_DEPTH_CR
                        set_parent_to_be_considered(
                            results_ptr, blk_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth,sq_size_idx,  s_depth);
#else
                        set_parent_to_be_considered(
                            results_ptr, blk_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth, s_depth);
#endif
#else
                        set_parent_to_be_considered(
                            results_ptr, blk_index, scs_ptr->seq_header.sb_size, s_depth);
#endif
                    // Add block indices of lower depth(s)
                    if (e_depth != 0)
#if TRACK_PER_DEPTH_DELTA
#if ADAPTIVE_DEPTH_CR
                        set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, blk_index, sb_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth,sq_size_idx, e_depth);
#else
                        set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, blk_index, sb_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth, e_depth);
#endif
#else
#if MULTI_PASS_PD_FOR_INCOMPLETE
                        set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, blk_index, sb_index, scs_ptr->seq_header.sb_size, e_depth);
#else
                        set_child_to_be_considered(
                            results_ptr, blk_index, scs_ptr->seq_header.sb_size, e_depth);
#endif
#endif
                }
            }
        }
        blk_index +=
            split_flag
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
}
#if OPT_BLOCK_INDICES_GEN_0
static void build_starting_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr, uint32_t sb_index) {

    MdcSbData *results_ptr = context_ptr->mdc_sb_array;

    results_ptr->leaf_count = 0;
    uint32_t blk_index = 0;
    uint32_t tot_d1_blocks;
    while (blk_index < scs_ptr->max_block_cnt) {
        tot_d1_blocks = 0;
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

        // SQ/NSQ block(s) filter based on the SQ size
        uint8_t is_block_tagged =
#if REMOVE_UNUSED_CODE_PH2
            (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE) ||
#else
            (blk_geom->sq_size == 128 && (pcs_ptr->slice_type == I_SLICE || pcs_ptr->parent_pcs_ptr->sb_64x64_simulated)) ||
#endif
            (blk_geom->sq_size == 4 && context_ptr->disallow_4x4)
            ? 0
            : 1;

        // split_flag is f(min_sq_size)
        int32_t min_sq_size = (context_ptr->disallow_4x4) ? 8 : 4;

        // SQ/NSQ block(s) filter based on the block validity
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_block_tagged) {
#if OPT_BLOCK_INDICES_GEN_2
#if OPT_BLOCK_INDICES_GEN_3
            tot_d1_blocks = (context_ptr->md_disallow_nsq) ||
#if NO_NSQ_ABOVE
                (blk_geom->sq_size >= 64 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_64x64) ||
                (blk_geom->sq_size >= 32 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_32x32) ||
                (blk_geom->sq_size >= 16 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_16x16) ||
#endif
#if NO_NSQ_B32
                (blk_geom->sq_size <= 64 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_64x64) ||
                (blk_geom->sq_size <= 32 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_32x32) ||
#endif
                (blk_geom->sq_size <= 8 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_8x8) ||
                (blk_geom->sq_size <= 16 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_16x16) ? 1 :
#else
            tot_d1_blocks = (context_ptr->md_disallow_nsq) ? 1 :
#endif
#else
            tot_d1_blocks = (pcs_ptr->parent_pcs_ptr->disallow_nsq) ? 1 :
#endif
#if OPT_BLOCK_INDICES_GEN_3
                (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_all_non_hv_nsq_blocks_below_16x16) ? 5 :
                (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_all_h4_v4_blocks_below_16x16) ? 17 :
#endif
                blk_geom->sq_size == 128
                ? 17
                : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;

#if NO_AB_HV4
            if (pcs_ptr->parent_pcs_ptr->disallow_HVA_HVB_HV4)
                tot_d1_blocks = MIN(5, tot_d1_blocks);

            if (pcs_ptr->parent_pcs_ptr->disallow_HV4)
                tot_d1_blocks = MIN(17, tot_d1_blocks);
#endif
            for (uint32_t idx = 0; idx < tot_d1_blocks; ++idx) {
                blk_geom = get_blk_geom_mds(blk_index);

                if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index]) {

                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = blk_index;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = tot_d1_blocks;

                    if (blk_geom->sq_size > min_sq_size)
                        results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag =
                        EB_TRUE;
                    else
                        results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag =
                        EB_FALSE;
                }
                blk_index++;
            }
            blk_index +=
                (d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] -
                    tot_d1_blocks);
        }
        else {
            blk_index +=
                (blk_geom->sq_size > min_sq_size)
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
        }
    }
}
#else
#if DEPTH_PART_CLEAN_UP
// Build the t=0 cand_block_array
void build_starting_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, EncDecContext *context_ptr, MdcSbData *results_ptr) {

    results_ptr->leaf_count = 0;
    uint32_t blk_index = 0;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

        uint8_t is_blk_allowed =
            (blk_geom->sq_size == 128 && (pcs_ptr->slice_type == I_SLICE || pcs_ptr->parent_pcs_ptr->sb_64x64_simulated)) ||
            (blk_geom->sq_size == 4   && pcs_ptr->parent_pcs_ptr->disallow_4x4)
                ? 0
                : 1;

        int32_t min_sq_size = (pcs_ptr->parent_pcs_ptr->disallow_4x4) ? 8 : 4;

        if (pcs_ptr->parent_pcs_ptr->sb_geom[context_ptr->sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed) {
            results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks =
                blk_geom->sq_size == 128
                ? 17
                : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;

            results_ptr->leaf_data_array[results_ptr->leaf_count].leaf_index =
                0; //valid only for square 85 world. will be removed.
            results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = blk_index;
            if (blk_geom->sq_size > min_sq_size)
                results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = EB_TRUE;
            else
                results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = EB_FALSE;
        }

        blk_index++;
    }


    pcs_ptr->parent_pcs_ptr->average_qp = (uint8_t)pcs_ptr->parent_pcs_ptr->picture_qp;
}
#endif
#endif

/* EncDec (Encode Decode) Kernel */
/*********************************************************************************
*
* @brief
*  The EncDec process contains both the mode decision and the encode pass engines
*  of the encoder. The mode decision encapsulates multiple partitioning decision (PD) stages
*  and multiple mode decision (MD) stages. At the end of the last mode decision stage,
*  the winning partition and modes combinations per block get reconstructed in the encode pass
*  operation which is part of the common section between the encoder and the decoder
*  Common encoder and decoder tasks such as Intra Prediction, Motion Compensated Prediction,
*  Transform, Quantization are performed in this process.
*
* @par Description:
*  The EncDec process operates on an SB basis.
*  The EncDec process takes as input the Motion Vector XY pairs candidates
*  and corresponding distortion estimates from the Motion Estimation process,
*  and the picture-level QP from the Rate Control process. All inputs are passed
*  through the picture structures: PictureControlSet and SequenceControlSet.
*  local structures of type EncDecContext and ModeDecisionContext contain all parameters
*  and results corresponding to the SuperBlock being processed.
*  each of the context structures is local to on thread and thus there's no risk of
*  affecting (changing) other SBs data in the process.
*
* @param[in] Vector
*  Motion Vector XY pairs from Motion Estimation process
*
* @param[in] Distortion Estimates
*  Distortion estimates from Motion Estimation process
*
* @param[in] Picture QP
*  Picture Quantization Parameter from Rate Control process
*
* @param[out] Blocks
*  The encode pass takes the selected partitioning and coding modes as input from mode decision for each
*  superblock and produces quantized transfrom coefficients for the residuals and the appropriate syntax
*  elements to be sent to the entropy coding engine
*
********************************************************************************/
void *enc_dec_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext *   thread_context_ptr = (EbThreadContext *)input_ptr;
    EncDecContext *     context_ptr        = (EncDecContext *)thread_context_ptr->priv;
    PictureControlSet * pcs_ptr;
    SequenceControlSet *scs_ptr;

    // Input
    EbObjectWrapper *enc_dec_tasks_wrapper_ptr;
    EncDecTasks *    enc_dec_tasks_ptr;

    // Output
    EbObjectWrapper *enc_dec_results_wrapper_ptr;
    EncDecResults *  enc_dec_results_ptr;
    // SB Loop variables
    SuperBlock *sb_ptr;
    uint16_t    sb_index;
    uint32_t    x_sb_index;
    uint32_t    y_sb_index;
    uint32_t    sb_origin_x;
    uint32_t    sb_origin_y;
    MdcSbData * mdc_ptr;

    // Segments
    uint16_t        segment_index;
    uint32_t        x_sb_start_index;
    uint32_t        y_sb_start_index;
    uint32_t        sb_start_index;
    uint32_t        sb_segment_count;
    uint32_t        sb_segment_index;
    uint32_t        segment_row_index;
    uint32_t        segment_band_index;
    uint32_t        segment_band_size;
    EncDecSegments *segments_ptr;

    segment_index = 0;

    for (;;) {
        // Get Mode Decision Results
        EB_GET_FULL_OBJECT(context_ptr->mode_decision_input_fifo_ptr, &enc_dec_tasks_wrapper_ptr);

        enc_dec_tasks_ptr = (EncDecTasks *)enc_dec_tasks_wrapper_ptr->object_ptr;
        pcs_ptr           = (PictureControlSet *)enc_dec_tasks_ptr->pcs_wrapper_ptr->object_ptr;
        scs_ptr           = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
        context_ptr->tile_group_index = enc_dec_tasks_ptr->tile_group_index;
        context_ptr->coded_sb_count   = 0;
        segments_ptr = pcs_ptr->enc_dec_segment_ctrl[context_ptr->tile_group_index];
        EbBool last_sb_flag           = EB_FALSE;
        // SB Constants
        uint8_t sb_sz      = (uint8_t)scs_ptr->sb_size_pix;
        uint8_t sb_size_log2 = (uint8_t)eb_log2f(sb_sz);
        context_ptr->sb_sz = sb_sz;
        uint32_t pic_width_in_sb = (pcs_ptr->parent_pcs_ptr->aligned_width + sb_sz - 1) >>
            sb_size_log2;
        uint16_t tile_group_width_in_sb = pcs_ptr->parent_pcs_ptr
                                              ->tile_group_info[context_ptr->tile_group_index]
                                              .tile_group_width_in_sb;
        uint32_t sb_row_index_start = 0, sb_row_index_count = 0;
        context_ptr->tot_intra_coded_area       = 0;

#if ADAPTIVE_NSQ_CR
        memset(context_ptr->md_context->part_cnt, 0, sizeof(uint32_t) * SSEG_NUM * (NUMBER_OF_SHAPES-1) * FB_NUM);
        generate_nsq_prob(pcs_ptr, context_ptr->md_context);
#endif
#if ADAPTIVE_DEPTH_CR
#if SOFT_CYCLES_REDUCTION
        memset(context_ptr->md_context->pred_depth_count, 0, sizeof(uint32_t) * DEPTH_DELTA_NUM * (NUMBER_OF_SHAPES-1));
#else
        memset(context_ptr->md_context->pred_depth_count, 0, sizeof(uint32_t) * DEPTH_DELTA_NUM);
#endif
        generate_depth_prob(pcs_ptr, context_ptr->md_context);
#endif
#if ADAPTIVE_TXT_CR
        memset( context_ptr->md_context->txt_cnt, 0, sizeof(uint32_t) * TXT_DEPTH_DELTA_NUM * TX_TYPES);
        generate_txt_prob(pcs_ptr, context_ptr->md_context);
#endif

        // Segment-loop
        while (assign_enc_dec_segments(segments_ptr,
                                       &segment_index,
                                       enc_dec_tasks_ptr,
                                       context_ptr->enc_dec_feedback_fifo_ptr) == EB_TRUE) {
            x_sb_start_index = segments_ptr->x_start_array[segment_index];
            y_sb_start_index = segments_ptr->y_start_array[segment_index];
            sb_start_index = y_sb_start_index * tile_group_width_in_sb + x_sb_start_index;
            sb_segment_count = segments_ptr->valid_sb_count_array[segment_index];

            segment_row_index = segment_index / segments_ptr->segment_band_count;
            segment_band_index =
                segment_index - segment_row_index * segments_ptr->segment_band_count;
            segment_band_size = (segments_ptr->sb_band_count * (segment_band_index + 1) +
                                 segments_ptr->segment_band_count - 1) /
                                segments_ptr->segment_band_count;

            // Reset Coding Loop State
            reset_mode_decision(scs_ptr,
                                context_ptr->md_context,
                                pcs_ptr,
                                context_ptr->tile_group_index,
                                segment_index);

            // Reset EncDec Coding State
            reset_enc_dec( // HT done
                context_ptr,
                pcs_ptr,
                scs_ptr,
                segment_index);

            if (pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL)
                ((EbReferenceObject *)
                     pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                    ->average_intensity = pcs_ptr->parent_pcs_ptr->average_intensity[0];
            for (y_sb_index = y_sb_start_index, sb_segment_index = sb_start_index;
                 sb_segment_index < sb_start_index + sb_segment_count;
                 ++y_sb_index) {
                for (x_sb_index = x_sb_start_index;
                     x_sb_index < tile_group_width_in_sb &&
                     (x_sb_index + y_sb_index < segment_band_size) &&
                     sb_segment_index < sb_start_index + sb_segment_count;
                     ++x_sb_index, ++sb_segment_index) {
                    uint16_t tile_group_y_sb_start =
                        pcs_ptr->parent_pcs_ptr->tile_group_info[context_ptr->tile_group_index]
                            .tile_group_sb_start_y;
                    uint16_t tile_group_x_sb_start =
                        pcs_ptr->parent_pcs_ptr->tile_group_info[context_ptr->tile_group_index]
                            .tile_group_sb_start_x;
                    sb_index = (uint16_t)((y_sb_index + tile_group_y_sb_start) * pic_width_in_sb +
                                          x_sb_index + tile_group_x_sb_start);
#if M8_4x4
                    sb_ptr = context_ptr->md_context->sb_ptr = pcs_ptr->sb_ptr_array[sb_index];
#else
                    sb_ptr   = pcs_ptr->sb_ptr_array[sb_index];
#endif
                    sb_origin_x = (x_sb_index + tile_group_x_sb_start) << sb_size_log2;
                    sb_origin_y = (y_sb_index + tile_group_y_sb_start) << sb_size_log2;
                    //printf("[%ld]:ED sb index %d, (%d, %d), encoded total sb count %d, ctx coded sb count %d\n",
                    //        pcs_ptr->picture_number,
                    //        sb_index, sb_origin_x, sb_origin_y,
                    //        pcs_ptr->enc_dec_coded_sb_count,
                    //        context_ptr->coded_sb_count);
                    context_ptr->tile_index             = sb_ptr->tile_info.tile_rs_index;
                    context_ptr->md_context->tile_index = sb_ptr->tile_info.tile_rs_index;

                    sb_row_index_start =
                        (x_sb_index + 1 == tile_group_width_in_sb && sb_row_index_count == 0)
                            ? y_sb_index
                            : sb_row_index_start;
                    sb_row_index_count = (x_sb_index + 1 == tile_group_width_in_sb)
                                             ? sb_row_index_count + 1
                                             : sb_row_index_count;
#if DEPTH_PART_CLEAN_UP
                    mdc_ptr = context_ptr->md_context->mdc_sb_array;
#else
                    mdc_ptr               = &pcs_ptr->mdc_sb_array[sb_index];
#endif
                    context_ptr->sb_index = sb_index;
#if SB_CLASSIFIER
                    context_ptr->md_context->sb_class = NONE_CLASS;
#endif

                    if (pcs_ptr->update_cdf) {
                        if (scs_ptr->seq_header.pic_based_rate_est &&
                            scs_ptr->enc_dec_segment_row_count_array[pcs_ptr->temporal_layer_index] == 1 &&
                            scs_ptr->enc_dec_segment_col_count_array[pcs_ptr->temporal_layer_index] == 1) {
                            if (sb_index == 0)
#if MD_FRAME_CONTEXT_MEM_OPT
                                pcs_ptr->ec_ctx_array[sb_index] =  pcs_ptr->md_frame_context;
#else
                                pcs_ptr->ec_ctx_array[sb_index] = *pcs_ptr->coeff_est_entropy_coder_ptr->fc;
#endif
                            else
                                pcs_ptr->ec_ctx_array[sb_index] = pcs_ptr->ec_ctx_array[sb_index - 1];
                        }
                        else {
#if REU_UPDATE
                            // Use the latest available CDF for the current SB
                            // Use the weighted average of left (3x) and top right (1x) if available.
                            int8_t top_right_available =
                                ((int32_t)(sb_origin_y >> MI_SIZE_LOG2) >
                                 sb_ptr->tile_info.mi_row_start) &&
                                ((int32_t)((sb_origin_x + (1 << sb_size_log2)) >> MI_SIZE_LOG2) <
                                 sb_ptr->tile_info.mi_col_end);

                            int8_t left_available = ((int32_t)(sb_origin_x >> MI_SIZE_LOG2) >
                                                     sb_ptr->tile_info.mi_col_start);

                            if (!left_available && !top_right_available)
                                pcs_ptr->ec_ctx_array[sb_index] =
#if MD_FRAME_CONTEXT_MEM_OPT
                                  pcs_ptr->md_frame_context;
#else
                                    *pcs_ptr->coeff_est_entropy_coder_ptr->fc;
#endif
                            else if (!left_available)
                                pcs_ptr->ec_ctx_array[sb_index] =
                                    pcs_ptr->ec_ctx_array[sb_index - pic_width_in_sb + 1];
                            else if (!top_right_available)
                                pcs_ptr->ec_ctx_array[sb_index] =
                                    pcs_ptr->ec_ctx_array[sb_index - 1];
                            else {
                                pcs_ptr->ec_ctx_array[sb_index] =
                                    pcs_ptr->ec_ctx_array[sb_index - 1];
                                avg_cdf_symbols(
                                    &pcs_ptr->ec_ctx_array[sb_index],
                                    &pcs_ptr->ec_ctx_array[sb_index - pic_width_in_sb + 1],
                                    AVG_CDF_WEIGHT_LEFT,
                                    AVG_CDF_WEIGHT_TOP);
#else
                            // Use the latest available CDF for the current SB
                            // Use the weighted average of left (3x) and top (1x) if available.
                            int8_t up_available = ((int32_t)(sb_origin_y >> MI_SIZE_LOG2) >
                                sb_ptr->tile_info.mi_row_start);
                            int8_t left_available = ((int32_t)(sb_origin_x >> MI_SIZE_LOG2) >
                                sb_ptr->tile_info.mi_col_start);
                            if (!left_available && !up_available)
                                pcs_ptr->ec_ctx_array[sb_index] =
                                *pcs_ptr->coeff_est_entropy_coder_ptr->fc;
                            else if (!left_available)
                                pcs_ptr->ec_ctx_array[sb_index] =
                                pcs_ptr->ec_ctx_array[sb_index - pic_width_in_sb];
                            else if (!up_available)
                                pcs_ptr->ec_ctx_array[sb_index] = pcs_ptr->ec_ctx_array[sb_index - 1];
                            else {
                                pcs_ptr->ec_ctx_array[sb_index] = pcs_ptr->ec_ctx_array[sb_index - 1];
                                avg_cdf_symbols(&pcs_ptr->ec_ctx_array[sb_index],
                                    &pcs_ptr->ec_ctx_array[sb_index - pic_width_in_sb],
                                    AVG_CDF_WEIGHT_LEFT,
                                    AVG_CDF_WEIGHT_TOP);
#endif
                            }
                        }

                        //in case of using 1 enc-dec segment, point to first SB data
                        uint32_t real_sb_idx = scs_ptr->seq_header.pic_based_rate_est &&
                            scs_ptr->enc_dec_segment_row_count_array[pcs_ptr->temporal_layer_index] == 1 &&
                            scs_ptr->enc_dec_segment_col_count_array[pcs_ptr->temporal_layer_index] == 1 ?
                            0 : sb_index;

                        // Copy all fileds from picture
                        pcs_ptr->rate_est_array[real_sb_idx] = *pcs_ptr->md_rate_estimation_array;

                        // Compute rate using latest CDFs
                        av1_estimate_syntax_rate(&pcs_ptr->rate_est_array[real_sb_idx],
                            pcs_ptr->slice_type == I_SLICE,
                            &pcs_ptr->ec_ctx_array[sb_index]);
                        av1_estimate_mv_rate(pcs_ptr,
                            &pcs_ptr->rate_est_array[real_sb_idx],
                            &pcs_ptr->ec_ctx_array[sb_index]);
                        av1_estimate_coefficients_rate(&pcs_ptr->rate_est_array[real_sb_idx],
                            &pcs_ptr->ec_ctx_array[sb_index]);

                        //let the candidate point to the new rate table.
                        uint32_t cand_index;
                        for (cand_index = 0; cand_index < MODE_DECISION_CANDIDATE_MAX_COUNT;
                            ++cand_index)
                            context_ptr->md_context->fast_candidate_ptr_array[cand_index]
                            ->md_rate_estimation_ptr = &pcs_ptr->rate_est_array[real_sb_idx];
                        context_ptr->md_context->md_rate_estimation_ptr =
                            &pcs_ptr->rate_est_array[real_sb_idx];
                    }
                    // Configure the SB
                    mode_decision_configure_sb(
#if QP2QINDEX
                        context_ptr->md_context, pcs_ptr, (uint8_t)sb_ptr->qindex);
#else
                        context_ptr->md_context, pcs_ptr, (uint8_t)sb_ptr->qp);
#endif
#if DEPTH_PART_CLEAN_UP && !OPT_BLOCK_INDICES_GEN_0
                    // Build the t=0 cand_block_array
                    build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr, mdc_ptr);
#endif
                    // Multi-Pass PD Path
                    // For each SB, all blocks are tested in PD0 (4421 blocks if 128x128 SB, and 1101 blocks if 64x64 SB).
                    // Then the PD0 predicted Partitioning Structure is refined by considering up to three refinements depths away from the predicted depth, both in the direction of smaller block sizes and in the direction of larger block sizes (up to Pred - 3 / Pred + 3 refinement). The selection of the refinement depth is performed using the cost
                    // deviation between the current depth cost and candidate depth cost. The generated blocks are used as input candidates to PD1.
                    // The PD1 predicted Partitioning Structure is also refined (up to Pred - 1 / Pred + 1 refinement) using the square (SQ) vs. non-square (NSQ) decision(s)
                    // inside the predicted depth and using coefficient information. The final set of blocks is evaluated in PD2 to output the final Partitioning Structure
#if DEPTH_PART_CLEAN_UP
#if ADD_NEW_MPPD_LEVEL
                    if ((pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_1 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4)
#if MULTI_PASS_PD_FOR_INCOMPLETE
                        ) {
#else
                        &&
#endif
#else
                    if ((pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_1 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3) &&
#endif
#if !MULTI_PASS_PD_FOR_INCOMPLETE
                        pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].is_complete_sb) {
#endif
#else
                    if ((pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_0 ||
                         pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_1 ||
                         pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2 ||
                         pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3) &&
                        pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].is_complete_sb) {
#endif
                        // Save a clean copy of the neighbor arrays
                        copy_neighbour_arrays(pcs_ptr,
                                              context_ptr->md_context,
                                              MD_NEIGHBOR_ARRAY_INDEX,
                                              MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                              0,
                                              sb_origin_x,
                                              sb_origin_y);

                        // [PD_PASS_0] Signal(s) derivation
                        context_ptr->md_context->pd_pass = PD_PASS_0;
                        signal_derivation_enc_dec_kernel_oq(
                            scs_ptr, pcs_ptr, context_ptr->md_context);

                        // [PD_PASS_0] Mode Decision - Reduce the total number of partitions to be tested in later stages.
                        // Input : mdc_blk_ptr built @ mdc process (up to 4421)
                        // Output: md_blk_arr_nsq reduced set of block(s)

#if OPT_BLOCK_INDICES_GEN_0
                        // Build the t=0 cand_block_array
                        build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif

                        // PD0 MD Tool(s) : Best ME candidate only as INTER candidate(s), DC only as INTRA candidate(s), Chroma blind, Spatial SSE,
                        // no MVP table generation, no fast rate @ full cost derivation, Md-Stage 0 and Md-Stage 2 using count=1 (i.e. only best md-stage-0 candidate)
                        mode_decision_sb(scs_ptr,
                                         pcs_ptr,
                                         mdc_ptr,
                                         sb_ptr,
                                         sb_origin_x,
                                         sb_origin_y,
                                         sb_index,
                                         context_ptr->md_context);
#if SB_CLASSIFIER
#if ADAPTIVE_DEPTH_CR
                        if (1) {
#else
                        if (pcs_ptr->slice_type != I_SLICE) {
#endif
#if !CLEANUP_CYCLE_ALLOCATION
                            set_sb_class_controls(context_ptr->md_context);
#endif
                            context_ptr->md_context->sb_class = determine_sb_class(
                                scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
                        }
#endif

                        // Perform Pred_0 depth refinement - Add blocks to be considered in the next stage(s) of PD based on depth cost.
                        perform_pred_depth_refinement(
                            scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

                        // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
#if !OPT_BLOCK_INDICES_GEN_0
#if DEPTH_PART_CLEAN_UP
                        build_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#else
                        build_cand_block_array(scs_ptr, pcs_ptr, sb_index);
#endif
#endif
                        // Reset neighnor information to current SB @ position (0,0)
                        copy_neighbour_arrays(pcs_ptr,
                                              context_ptr->md_context,
                                              MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                              MD_NEIGHBOR_ARRAY_INDEX,
                                              0,
                                              sb_origin_x,
                                              sb_origin_y);

#if DEPTH_PART_CLEAN_UP
#if ADD_NEW_MPPD_LEVEL
                        if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_1 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4) {
#else
                        if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_1 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3) {
#endif
#else
                        if (pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_1 ||
                            pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2 ||
                            pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3) {
#endif
                            // [PD_PASS_1] Signal(s) derivation
                            context_ptr->md_context->pd_pass = PD_PASS_1;
                            signal_derivation_enc_dec_kernel_oq(
                                scs_ptr, pcs_ptr, context_ptr->md_context);
#if OPT_BLOCK_INDICES_GEN_0
                            // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                            build_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif

                            // [PD_PASS_1] Mode Decision - Further reduce the number of
                            // partitions to be considered in later PD stages. This pass uses more accurate
                            // info than PD0 to give a better PD estimate.
                            // Input : mdc_blk_ptr built @ PD0 refinement
                            // Output: md_blk_arr_nsq reduced set of block(s)

                            // PD1 MD Tool(s) : ME and Predictive ME only as INTER candidate(s) but MRP blind (only reference index 0 for motion compensation),
                            // DC only as INTRA candidate(s)
                            mode_decision_sb(scs_ptr,
                                             pcs_ptr,
                                             mdc_ptr,
                                             sb_ptr,
                                             sb_origin_x,
                                             sb_origin_y,
                                             sb_index,
                                             context_ptr->md_context);

                            // Perform Pred_1 depth refinement - Add blocks to be considered in the next stage(s) of PD based on depth cost.
                            perform_pred_depth_refinement(
                                scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

                            // Re-build mdc_blk_ptr for the 3rd PD Pass [PD_PASS_2]
#if DEPTH_PART_CLEAN_UP
                            build_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#else
                            build_cand_block_array(scs_ptr, pcs_ptr, sb_index);
#endif
                            // Reset neighnor information to current SB @ position (0,0)
                            copy_neighbour_arrays(pcs_ptr,
                                                  context_ptr->md_context,
                                                  MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                                  MD_NEIGHBOR_ARRAY_INDEX,
                                                  0,
                                                  sb_origin_x,
                                                  sb_origin_y);
                        }
                    }
#if OPT_BLOCK_INDICES_GEN_0 && !OPT_BLOCK_INDICES_GEN_4
                    else
                        // Build the t=0 cand_block_array
                        build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif
                    // [PD_PASS_2] Signal(s) derivation
                    context_ptr->md_context->pd_pass = PD_PASS_2;
                    signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);
#if OPT_BLOCK_INDICES_GEN_0
                    // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                    if(pcs_ptr->parent_pcs_ptr->multi_pass_pd_level != MULTI_PASS_PD_OFF)
                    build_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#if OPT_BLOCK_INDICES_GEN_4
                    else
                        // Build the t=0 cand_block_array
                        build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif
#endif

                    // [PD_PASS_2] Mode Decision - Obtain the final partitioning decision using more accurate info
                    // than previous stages.  Reduce the total number of partitions to 1.
                    // Input : mdc_blk_ptr built @ PD1 refinement
                    // Output: md_blk_arr_nsq reduced set of block(s)

                    // PD2 MD Tool(s): default MD Tool(s)

                    mode_decision_sb(scs_ptr,
                                     pcs_ptr,
                                     mdc_ptr,
                                     sb_ptr,
                                     sb_origin_x,
                                     sb_origin_y,
                                     sb_index,
                                     context_ptr->md_context);
#if ADAPTIVE_NSQ_CR
                    generate_statistics_nsq(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif
#if ADAPTIVE_DEPTH_CR
                    generate_statistics_depth(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif
#if ADAPTIVE_TXT_CR
                    generate_statistics_txt(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif
#if !QP2QINDEX
                    // Configure the SB
                    enc_dec_configure_sb(context_ptr, sb_ptr, pcs_ptr, (uint8_t)sb_ptr->qp);
#endif

#if NO_ENCDEC
                    no_enc_dec_pass(scs_ptr,
                                    pcs_ptr,
                                    sb_ptr,
                                    sb_index,
                                    sb_origin_x,
                                    sb_origin_y,
                                    sb_ptr->qp,
                                    context_ptr);
#else
                    // Encode Pass
                    av1_encode_pass(
                        scs_ptr, pcs_ptr, sb_ptr, sb_index, sb_origin_x, sb_origin_y, context_ptr);
#endif

                    context_ptr->coded_sb_count++;
                    if (pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL)
                        ((EbReferenceObject *)
                             pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                            ->intra_coded_area_sb[sb_index] = (uint8_t)(
                            (100 * context_ptr->intra_coded_area_sb[sb_index]) / (64 * 64));
                }
                x_sb_start_index = (x_sb_start_index > 0) ? x_sb_start_index - 1 : 0;
            }
        }

        eb_block_on_mutex(pcs_ptr->intra_mutex);
        pcs_ptr->intra_coded_area += (uint32_t)context_ptr->tot_intra_coded_area;
        pcs_ptr->enc_dec_coded_sb_count += (uint32_t)context_ptr->coded_sb_count;
        last_sb_flag = (pcs_ptr->sb_total_count_pix == pcs_ptr->enc_dec_coded_sb_count);
        eb_release_mutex(pcs_ptr->intra_mutex);

        if (last_sb_flag) {
            // Copy film grain data from parent picture set to the reference object for further reference
            if (scs_ptr->seq_header.film_grain_params_present) {
                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
                    pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr) {
                    ((EbReferenceObject *)
                         pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                        ->film_grain_params = pcs_ptr->parent_pcs_ptr->frm_hdr.film_grain_params;
                }
            }
            if (pcs_ptr->parent_pcs_ptr->frame_end_cdf_update_mode &&
                pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
                pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr)
                for (int frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame)
                    ((EbReferenceObject *)
                         pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                        ->global_motion[frame] = pcs_ptr->parent_pcs_ptr->global_motion[frame];
            eb_memcpy(pcs_ptr->parent_pcs_ptr->av1x->sgrproj_restore_cost,
                      context_ptr->md_rate_estimation_ptr->sgrproj_restore_fac_bits,
                      2 * sizeof(int32_t));
            eb_memcpy(pcs_ptr->parent_pcs_ptr->av1x->switchable_restore_cost,
                      context_ptr->md_rate_estimation_ptr->switchable_restore_fac_bits,
                      3 * sizeof(int32_t));
            eb_memcpy(pcs_ptr->parent_pcs_ptr->av1x->wiener_restore_cost,
                      context_ptr->md_rate_estimation_ptr->wiener_restore_fac_bits,
                      2 * sizeof(int32_t));
#if QP2QINDEX
#if TPL_LA_LAMBDA_SCALING
            pcs_ptr->parent_pcs_ptr->av1x->rdmult =
                context_ptr->pic_full_lambda[(context_ptr->bit_depth == EB_10BIT) ? EB_10_BIT_MD
                                                                                  : EB_8_BIT_MD];
#else
            pcs_ptr->parent_pcs_ptr->av1x->rdmult = context_ptr->pic_full_lambda;
#endif
#else
            pcs_ptr->parent_pcs_ptr->av1x->rdmult = context_ptr->full_lambda;
#endif
#if DECOUPLE_ME_RES
            eb_release_object(pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr);
            pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr = (EbObjectWrapper *)EB_NULL;
#endif
        }

        if (last_sb_flag) {
            // Get Empty EncDec Results
            eb_get_empty_object(context_ptr->enc_dec_output_fifo_ptr, &enc_dec_results_wrapper_ptr);
            enc_dec_results_ptr = (EncDecResults *)enc_dec_results_wrapper_ptr->object_ptr;
            enc_dec_results_ptr->pcs_wrapper_ptr = enc_dec_tasks_ptr->pcs_wrapper_ptr;
            //CHKN these are not needed for DLF
            enc_dec_results_ptr->completed_sb_row_index_start = 0;
            enc_dec_results_ptr->completed_sb_row_count =
                ((pcs_ptr->parent_pcs_ptr->aligned_height + scs_ptr->sb_size_pix - 1) >> sb_size_log2);
            // Post EncDec Results
            eb_post_full_object(enc_dec_results_wrapper_ptr);
        }
        // Release Mode Decision Results
        eb_release_object(enc_dec_tasks_wrapper_ptr);
    }
    return NULL;
}

void eb_av1_add_film_grain(EbPictureBufferDesc *src, EbPictureBufferDesc *dst,
                           AomFilmGrain *film_grain_ptr) {
    uint8_t *luma, *cb, *cr;
    int32_t  height, width, luma_stride, chroma_stride;
    int32_t  use_high_bit_depth = 0;
    int32_t  chroma_subsamp_x   = 0;
    int32_t  chroma_subsamp_y   = 0;

    AomFilmGrain params = *film_grain_ptr;

    switch (src->bit_depth) {
    case EB_8BIT:
        params.bit_depth   = 8;
        use_high_bit_depth = 0;
        chroma_subsamp_x   = 1;
        chroma_subsamp_y   = 1;
        break;
    case EB_10BIT:
        params.bit_depth   = 10;
        use_high_bit_depth = 1;
        chroma_subsamp_x   = 1;
        chroma_subsamp_y   = 1;
        break;
    default: //todo: Throw an error if unknown format?
        params.bit_depth   = 10;
        use_high_bit_depth = 1;
        chroma_subsamp_x   = 1;
        chroma_subsamp_y   = 1;
    }

    dst->max_width  = src->max_width;
    dst->max_height = src->max_height;

    fgn_copy_rect(
        src->buffer_y + ((src->origin_y * src->stride_y + src->origin_x) << use_high_bit_depth),
        src->stride_y,
        dst->buffer_y + ((dst->origin_y * dst->stride_y + dst->origin_x) << use_high_bit_depth),
        dst->stride_y,
        dst->width,
        dst->height,
        use_high_bit_depth);

    fgn_copy_rect(src->buffer_cb + ((src->stride_cb * (src->origin_y >> chroma_subsamp_y) +
                                     (src->origin_x >> chroma_subsamp_x))
                                    << use_high_bit_depth),
                  src->stride_cb,
                  dst->buffer_cb + ((dst->stride_cb * (dst->origin_y >> chroma_subsamp_y) +
                                     (dst->origin_x >> chroma_subsamp_x))
                                    << use_high_bit_depth),
                  dst->stride_cb,
                  dst->width >> chroma_subsamp_x,
                  dst->height >> chroma_subsamp_y,
                  use_high_bit_depth);

    fgn_copy_rect(src->buffer_cr + ((src->stride_cr * (src->origin_y >> chroma_subsamp_y) +
                                     (src->origin_x >> chroma_subsamp_x))
                                    << use_high_bit_depth),
                  src->stride_cr,
                  dst->buffer_cr + ((dst->stride_cr * (dst->origin_y >> chroma_subsamp_y) +
                                     (dst->origin_x >> chroma_subsamp_x))
                                    << use_high_bit_depth),
                  dst->stride_cr,
                  dst->width >> chroma_subsamp_x,
                  dst->height >> chroma_subsamp_y,
                  use_high_bit_depth);

    luma = dst->buffer_y + ((dst->origin_y * dst->stride_y + dst->origin_x) << use_high_bit_depth);
    cb   = dst->buffer_cb + ((dst->stride_cb * (dst->origin_y >> chroma_subsamp_y) +
                            (dst->origin_x >> chroma_subsamp_x))
                           << use_high_bit_depth);
    cr   = dst->buffer_cr + ((dst->stride_cr * (dst->origin_y >> chroma_subsamp_y) +
                            (dst->origin_x >> chroma_subsamp_x))
                           << use_high_bit_depth);

    luma_stride   = dst->stride_y;
    chroma_stride = dst->stride_cb;

    width  = dst->width;
    height = dst->height;

    eb_av1_add_film_grain_run(&params,
                              luma,
                              cb,
                              cr,
                              height,
                              width,
                              luma_stride,
                              chroma_stride,
                              use_high_bit_depth,
                              chroma_subsamp_y,
                              chroma_subsamp_x);
    return;
}
