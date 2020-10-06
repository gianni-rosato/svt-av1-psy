/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>

#include "EbEncHandle.h"
#include "EbEncDecTasks.h"
#include "EbEncDecResults.h"
#include "EbCodingLoop.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbUtility.h"
#include "grainSynthesis.h"
//To fix warning C4013: 'svt_convert_16bit_to_8bit' undefined; assuming extern returning int
#include "common_dsp_rtcd.h"
#include "EbRateDistortionCost.h"
#include "EbPictureDecisionProcess.h"
#include "firstpass.h"

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
    int8_t       enable_hbd_mode_decision = static_config->enable_hbd_mode_decision;

    EncDecContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = enc_dec_context_dctor;

    context_ptr->is_16bit     = is_16bit;
    context_ptr->color_format = color_format;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_input_fifo_ptr = eb_system_resource_get_consumer_fifo(
        enc_handle_ptr->enc_dec_tasks_resource_ptr, index);
    context_ptr->enc_dec_output_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_results_resource_ptr, index);
    context_ptr->enc_dec_feedback_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_tasks_resource_ptr, tasks_index);
    context_ptr->picture_demux_output_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);

    // MD rate Estimation tables
    EB_MALLOC(context_ptr->md_rate_estimation_ptr, sizeof(MdRateEstimationContext));
    context_ptr->is_md_rate_estimation_ptr_owner = EB_TRUE;

    // Prediction Buffer
    context_ptr->input_sample16bit_buffer = NULL;
    if (is_16bit || static_config->is_16bit_pipeline)
        EB_NEW(context_ptr->input_sample16bit_buffer,
               eb_picture_buffer_desc_ctor,
               &(EbPictureBufferDescInitData){
                   .buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK,
                   .max_width          = SB_STRIDE_Y,
                   .max_height         = SB_STRIDE_Y,
                   .bit_depth          = EB_16BIT,
                   .left_padding       = 0,
                   .right_padding      = 0,
                   .top_padding        = 0,
                   .bot_padding        = 0,
                   .split_mode         = EB_FALSE,
                   .color_format       = color_format,
               });

    // Scratch Coeff Buffer
    EbPictureBufferDescInitData init_32bit_data = {
        .buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK,
        .max_width          = SB_STRIDE_Y,
        .max_height         = SB_STRIDE_Y,
        .bit_depth          = EB_32BIT,
        .color_format       = color_format,
        .left_padding       = 0,
        .right_padding      = 0,
        .top_padding        = 0,
        .bot_padding        = 0,
        .split_mode         = EB_FALSE,
    };

    EB_NEW(context_ptr->inverse_quant_buffer, eb_picture_buffer_desc_ctor, &init_32bit_data);
    EB_NEW(context_ptr->transform_buffer, eb_picture_buffer_desc_ctor, &init_32bit_data);
    EB_NEW(context_ptr->residual_buffer,
           eb_picture_buffer_desc_ctor,
           &(EbPictureBufferDescInitData){
               .buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK,
               .max_width          = SB_STRIDE_Y,
               .max_height         = SB_STRIDE_Y,
               .bit_depth          = EB_16BIT,
               .color_format       = color_format,
               .left_padding       = 0,
               .right_padding      = 0,
               .top_padding        = 0,
               .bot_padding        = 0,
               .split_mode         = EB_FALSE,
           });

    // Mode Decision Context
    EB_NEW(context_ptr->md_context,
           mode_decision_context_ctor,
           color_format,
           static_config->super_block_size,
           0,
           0,
           enable_hbd_mode_decision == DEFAULT ? 2 : enable_hbd_mode_decision ,
           static_config->screen_content_mode);
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
    uint16_t tile_group_idx = context_ptr->tile_group_index;
    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
        pcs_ptr,
        &context_ptr->pic_fast_lambda[EB_8_BIT_MD],
        &context_ptr->pic_full_lambda[EB_8_BIT_MD],
        8,
        pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx,
        EB_TRUE);

    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
        pcs_ptr,
        &context_ptr->pic_fast_lambda[EB_10_BIT_MD],
        &context_ptr->pic_full_lambda[EB_10_BIT_MD],
        10,
        pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx,
        EB_TRUE);
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
  double ssim_n, ssim_d;
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

  ssim_n = (2.0 * sum_s * sum_r + c1) *
           (2.0 * count * sum_sxr - 2.0 * sum_s * sum_r + c2);

  ssim_d = ((double)sum_s * sum_s + (double)sum_r * sum_r + c1) *
           ((double)count * sum_sq_s - (double)sum_s * sum_s +
            (double)count * sum_sq_r - (double)sum_r * sum_r + c2);

  return ssim_n / ssim_d;
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
    assert(samples > 0);
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
  assert(samples > 0);
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
        // Non visible Reference samples should be overwritten by the last visible line of pixels
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
        if (pcs_ptr->hbd_mode_decision != EB_10_BIT_MD) {

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
        svt_convert_16bit_to_8bit(buf_16bit,
            ref_pic_16bit_ptr->stride_y,
            buf_8bit,
            ref_pic_ptr->stride_y,
            ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1),
            ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1));

        //CB
        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cb);
        buf_8bit = ref_pic_ptr->buffer_cb;
        svt_convert_16bit_to_8bit(buf_16bit,
            ref_pic_16bit_ptr->stride_cb,
            buf_8bit,
            ref_pic_ptr->stride_cb,
            (ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1)) >> 1,
            (ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1)) >> 1);

        //CR
        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cr);
        buf_8bit = ref_pic_ptr->buffer_cr;
        svt_convert_16bit_to_8bit(buf_16bit,
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
    reference_object->r0 = pcs_ptr->parent_pcs_ptr->r0;
}

void copy_statistics_to_ref_obj_ect(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    pcs_ptr->intra_coded_area =
        (100 * pcs_ptr->intra_coded_area) /
        (pcs_ptr->parent_pcs_ptr->aligned_width * pcs_ptr->parent_pcs_ptr->aligned_height);
    memcpy(((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                    ->ref_part_cnt, pcs_ptr->part_cnt, sizeof(uint32_t) * (NUMBER_OF_SHAPES-1) * FB_NUM *SSEG_NUM);
    memcpy(((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->ref_pred_depth_count, pcs_ptr->pred_depth_count, sizeof(uint32_t) * DEPTH_DELTA_NUM * (NUMBER_OF_SHAPES-1));
    memcpy(((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->ref_txt_cnt, pcs_ptr->txt_cnt, sizeof(uint32_t) * TXT_DEPTH_DELTA_NUM *TX_TYPES);
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

void set_obmc_controls(ModeDecisionContext *mdctxt, uint8_t obmc_mode) {

    ObmcControls*obmc_ctrls = &mdctxt->obmc_ctrls;

    switch (obmc_mode)
    {
    case 0:
        obmc_ctrls->enabled = 0;
        obmc_ctrls->pme_best_ref = 0;
        obmc_ctrls->me_count = 0;
        obmc_ctrls->mvp_ref_count = 0;
        obmc_ctrls->near_count = 0;
        break;
    case 1:
        obmc_ctrls->enabled = 1;
        obmc_ctrls->pme_best_ref = 0;
        obmc_ctrls->me_count = ~0;
        obmc_ctrls->mvp_ref_count = 4;
        obmc_ctrls->near_count = 3;
        break;
    case 2:
        obmc_ctrls->enabled = 1;
        obmc_ctrls->pme_best_ref = 0;
        obmc_ctrls->me_count = ~0;
        obmc_ctrls->mvp_ref_count = 4;
        obmc_ctrls->near_count = 3;
        break;
    case 3:
        obmc_ctrls->enabled = 1;
        obmc_ctrls->pme_best_ref = 1;
        obmc_ctrls->me_count = 1;
        obmc_ctrls->mvp_ref_count = 1;
        obmc_ctrls->near_count = 1;
        break;
    default:
        assert(0);
        break;
    }


}
void set_block_based_depth_refinement_controls(ModeDecisionContext *mdctxt, uint8_t block_based_depth_refinement_level) {

    DepthRefinementCtrls *depth_refinement_ctrls = &mdctxt->depth_refinement_ctrls;

    switch (block_based_depth_refinement_level)
    {
    case 0:
        depth_refinement_ctrls->enabled = 0;
        break;
    case 1:
        depth_refinement_ctrls->enabled = 1;
        depth_refinement_ctrls->parent_to_current_th = -10;
        depth_refinement_ctrls->sub_to_current_th = 5;
        break;
    default:
        assert(0);
        break;
    }
}
/*
 * Control NSQ search
 */
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
        break;

    case 2:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 15;
        md_nsq_motion_search_ctrls->full_pel_search_height = 15;
        break;
    case 3:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 11;
        md_nsq_motion_search_ctrls->full_pel_search_height = 11;
        break;
    case 4:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 7;
        md_nsq_motion_search_ctrls->full_pel_search_height = 7;
        break;
    default:
        assert(0);
        break;
    }
}
void md_pme_search_controls(ModeDecisionContext *mdctxt, uint8_t md_pme_level) {

    MdPmeCtrls *md_pme_ctrls = &mdctxt->md_pme_ctrls;

    switch (md_pme_level)
    {
    case 0:
        md_pme_ctrls->enabled = 0;
        break;

    case 1:
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 1;
        md_pme_ctrls->full_pel_search_width = 15;
        md_pme_ctrls->full_pel_search_height = 15;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
        break;

    case 2:
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 1;
        md_pme_ctrls->full_pel_search_width = 7;
        md_pme_ctrls->full_pel_search_height = 5;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
        break;

    case 3:
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 1;
        md_pme_ctrls->full_pel_search_width = 7;
        md_pme_ctrls->full_pel_search_height = 5;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 100;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
        break;

    default:
        assert(0);
        break;
    }
}
/*
 * Control Adaptive ME search
 */
void md_sq_motion_search_controls(ModeDecisionContext *mdctxt, uint8_t md_sq_mv_search_level) {

    MdSqMotionSearchCtrls *md_sq_me_ctrls = &mdctxt->md_sq_me_ctrls;

    switch (md_sq_mv_search_level) {
    case 0:
        md_sq_me_ctrls->enabled = 0;
        break;
    case 1:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled = 1;
        md_sq_me_ctrls->sprs_lev0_step = 4;
        md_sq_me_ctrls->sprs_lev0_w = 15;
        md_sq_me_ctrls->sprs_lev0_h = 15;
        md_sq_me_ctrls->max_sprs_lev0_w = 150;
        md_sq_me_ctrls->max_sprs_lev0_h = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 500;

        md_sq_me_ctrls->sprs_lev1_enabled = 1;
        md_sq_me_ctrls->sprs_lev1_step = 2;
        md_sq_me_ctrls->sprs_lev1_w = 4;
        md_sq_me_ctrls->sprs_lev1_h = 4;
        md_sq_me_ctrls->max_sprs_lev1_w = 50;
        md_sq_me_ctrls->max_sprs_lev1_h = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 500;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step = 1;
        md_sq_me_ctrls->sprs_lev2_w = 3;
        md_sq_me_ctrls->sprs_lev2_h = 3;
        break;
    case 2:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled = 1;
        md_sq_me_ctrls->sprs_lev0_step = 4;
        md_sq_me_ctrls->sprs_lev0_w = 15;
        md_sq_me_ctrls->sprs_lev0_h = 15;
        md_sq_me_ctrls->max_sprs_lev0_w = 150;
        md_sq_me_ctrls->max_sprs_lev0_h = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 400;

        md_sq_me_ctrls->sprs_lev1_enabled = 1;
        md_sq_me_ctrls->sprs_lev1_step = 2;
        md_sq_me_ctrls->sprs_lev1_w = 4;
        md_sq_me_ctrls->sprs_lev1_h = 4;
        md_sq_me_ctrls->max_sprs_lev1_w = 50;
        md_sq_me_ctrls->max_sprs_lev1_h = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 400;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step = 1;
        md_sq_me_ctrls->sprs_lev2_w = 3;
        md_sq_me_ctrls->sprs_lev2_h = 3;
        break;
    case 3:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled = 1;
        md_sq_me_ctrls->sprs_lev0_step = 4;
        md_sq_me_ctrls->sprs_lev0_w = 15;
        md_sq_me_ctrls->sprs_lev0_h = 15;
        md_sq_me_ctrls->max_sprs_lev0_w = 150;
        md_sq_me_ctrls->max_sprs_lev0_h = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 300;

        md_sq_me_ctrls->sprs_lev1_enabled = 1;
        md_sq_me_ctrls->sprs_lev1_step = 2;
        md_sq_me_ctrls->sprs_lev1_w = 4;
        md_sq_me_ctrls->sprs_lev1_h = 4;
        md_sq_me_ctrls->max_sprs_lev1_w = 50;
        md_sq_me_ctrls->max_sprs_lev1_h = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 300;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step = 1;
        md_sq_me_ctrls->sprs_lev2_w = 3;
        md_sq_me_ctrls->sprs_lev2_h = 3;
        break;
    case 4:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;
        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled = 1;
        md_sq_me_ctrls->sprs_lev0_step = 4;
        md_sq_me_ctrls->sprs_lev0_w = 15;
        md_sq_me_ctrls->sprs_lev0_h = 15;
        md_sq_me_ctrls->max_sprs_lev0_w = 150;
        md_sq_me_ctrls->max_sprs_lev0_h = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 100;

        md_sq_me_ctrls->sprs_lev1_enabled = 1;
        md_sq_me_ctrls->sprs_lev1_step = 2;
        md_sq_me_ctrls->sprs_lev1_w = 4;
        md_sq_me_ctrls->sprs_lev1_h = 4;
        md_sq_me_ctrls->max_sprs_lev1_w = 50;
        md_sq_me_ctrls->max_sprs_lev1_h = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 100;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step = 1;
        md_sq_me_ctrls->sprs_lev2_w = 3;
        md_sq_me_ctrls->sprs_lev2_h = 3;
        break;
    default:
        assert(0);
        break;
    }
}
/*
 * Control Subpel search of ME MV(s)
 */
void md_subpel_me_controls(ModeDecisionContext *mdctxt, uint8_t md_subpel_me_level) {
    MdSubPelSearchCtrls *md_subpel_me_ctrls = &mdctxt->md_subpel_me_ctrls;

    switch (md_subpel_me_level) {
    case 0:
        md_subpel_me_ctrls->enabled = 0;
        break;
    case 1:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->eight_pel_search_enabled = 1;
        break;
    case 2:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        break;
    case 3:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        break;
    default: assert(0); break;
    }
}

/*
 * Control Subpel search of PME MV(s)
 */
void md_subpel_pme_controls(ModeDecisionContext *mdctxt, uint8_t md_subpel_pme_level) {
    MdSubPelSearchCtrls *md_subpel_pme_ctrls = &mdctxt->md_subpel_pme_ctrls;

    switch (md_subpel_pme_level) {
    case 0:
        md_subpel_pme_ctrls->enabled = 0;
        break;
    case 1:
        md_subpel_pme_ctrls->enabled = 1;
        md_subpel_pme_ctrls->subpel_search_type = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->eight_pel_search_enabled = 1;
        break;
    case 2:
        md_subpel_pme_ctrls->enabled = 1;
        md_subpel_pme_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->eight_pel_search_enabled = 0;
        break;
    case 3:
        md_subpel_pme_ctrls->enabled = 1;
        md_subpel_pme_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 1;
        md_subpel_pme_ctrls->eight_pel_search_enabled = 0;
        break;
    default: assert(0); break;
    }
}
void coeff_based_switch_md_controls(ModeDecisionContext *mdctxt, uint8_t switch_md_mode_based_on_sq_coeff_level) {
    CoeffBSwMdCtrls *coeffb_sw_md_ctrls = &mdctxt->cb_sw_md_ctrls;

    switch (switch_md_mode_based_on_sq_coeff_level) {
    case 0: coeffb_sw_md_ctrls->enabled = 0; break;
    case 1:
        coeffb_sw_md_ctrls->enabled = 1;
        coeffb_sw_md_ctrls->mode_offset = 3;
        coeffb_sw_md_ctrls->skip_block = 0;
        break;
    case 2:
        coeffb_sw_md_ctrls->enabled = 1;
        coeffb_sw_md_ctrls->mode_offset = 4;
        coeffb_sw_md_ctrls->skip_block = 0;
        break;
    case 3:
        coeffb_sw_md_ctrls->enabled = 1;
        coeffb_sw_md_ctrls->mode_offset = 4;
        coeffb_sw_md_ctrls->skip_block = 1;
        break;
    default: assert(0); break;
    }
}
/******************************************************
* Derive SB classifier thresholds
******************************************************/

uint8_t nsq_cycles_reduction_th[19] = {
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
void adaptive_md_cycles_redcution_controls(ModeDecisionContext *mdctxt, uint8_t adaptive_md_cycles_red_mode) {
    AMdCycleRControls* adaptive_md_cycles_red_ctrls = &mdctxt->admd_cycles_red_ctrls;
    switch (adaptive_md_cycles_red_mode)
    {
    case 0:
        adaptive_md_cycles_red_ctrls->enabled = 0;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 0;
        adaptive_md_cycles_red_ctrls->switch_mode_th = 0;
        adaptive_md_cycles_red_ctrls->mode_offset = 0;
        break;
    case 1:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 0;
        adaptive_md_cycles_red_ctrls->switch_mode_th = 300;
        adaptive_md_cycles_red_ctrls->mode_offset = 1;
        break;
    case 2:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 100;
        adaptive_md_cycles_red_ctrls->switch_mode_th = 700;
        adaptive_md_cycles_red_ctrls->mode_offset = 2;
        break;
    case 3:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 200;
        adaptive_md_cycles_red_ctrls->switch_mode_th = 1000;
        adaptive_md_cycles_red_ctrls->mode_offset = 2;
        break;
    case 4:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 300;
        adaptive_md_cycles_red_ctrls->switch_mode_th = 300;
        adaptive_md_cycles_red_ctrls->mode_offset = 1;
        break;
    case 5:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 500;
        adaptive_md_cycles_red_ctrls->switch_mode_th = 1500;
        adaptive_md_cycles_red_ctrls->mode_offset = 1;
        break;
    case 6:
        adaptive_md_cycles_red_ctrls->enabled = 1;
        adaptive_md_cycles_red_ctrls->skip_nsq_th = 750;
        adaptive_md_cycles_red_ctrls->switch_mode_th = 1500;
        adaptive_md_cycles_red_ctrls->mode_offset = 1;
        break;
    default:
        assert(0);
        break;
    }
}
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
/******************************************************
* Derive EncDec Settings for OQ
Input   : encoder mode and pd pass
Output  : EncDec Kernel signal(s)
******************************************************/
EbErrorType signal_derivation_enc_dec_kernel_oq(
    SequenceControlSet *sequence_control_set_ptr,
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr,
    EbEncMode mode_offset) {
    EbErrorType return_error = EB_ErrorNone;
    EbEncMode enc_mode;
    if (mode_offset)
        enc_mode = MIN(ENC_M8, pcs_ptr->enc_mode + mode_offset);
    else
        enc_mode = pcs_ptr->enc_mode;
    uint8_t pd_pass = context_ptr->pd_pass;

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
        // Do not use cycles reduction algorithms in 480p and below
        else if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
            context_ptr->enable_area_based_cycles_allocation = 0;
        else
            context_ptr->enable_area_based_cycles_allocation = 1;
    }
    // Tx_search Level for Luma                       Settings
    // TX_SEARCH_DCT_DCT_ONLY                         DCT_DCT only
    // TX_SEARCH_DCT_TX_TYPES                         Tx search DCT type(s): DCT_DCT, V_DCT, H_DCT
    // TX_SEARCH_ALL_TX_TYPES                         Tx search all type(s)
    if (pd_pass == PD_PASS_0)
        context_ptr->tx_search_level = TX_SEARCH_DCT_DCT_ONLY;
    else if (pd_pass == PD_PASS_1)
        context_ptr->tx_search_level = TX_SEARCH_DCT_DCT_ONLY;
    else
        if (enc_mode <= ENC_M6)
            context_ptr->tx_search_level = TX_SEARCH_ALL_TX_TYPES;
        else
            if (pcs_ptr->parent_pcs_ptr->slice_type == I_SLICE)
                context_ptr->tx_search_level = TX_SEARCH_ALL_TX_TYPES;
            else
                context_ptr->tx_search_level = TX_SEARCH_DCT_TX_TYPES;
    uint8_t txt_cycles_reduction_level = 0;
    if (pcs_ptr->parent_pcs_ptr->slice_type == I_SLICE) {
        txt_cycles_reduction_level = 0;
    }
    else {
        if (pd_pass == PD_PASS_0)
            txt_cycles_reduction_level = 0;
        else if (pd_pass == PD_PASS_1)
            txt_cycles_reduction_level = 0;
        else if (enc_mode <= ENC_M4)
            txt_cycles_reduction_level = 0;
        else
            txt_cycles_reduction_level = 5;
    }
    set_txt_cycle_reduction_controls(context_ptr, txt_cycles_reduction_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->interpolation_search_level = IFS_OFF;
    else if (pd_pass == PD_PASS_1)
        context_ptr->interpolation_search_level = IFS_OFF;
    else
        if (enc_mode <= ENC_M6)
            context_ptr->interpolation_search_level = IFS_MDS1;
        else
            context_ptr->interpolation_search_level = IFS_MDS3;
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
        if (enc_mode <= ENC_M5)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else
            context_ptr->chroma_level = CHROMA_MODE_1;
    }
    else // use specified level
        context_ptr->chroma_level =
        sequence_control_set_ptr->static_config.set_chroma_mode;

    // Chroma independent modes search
    // Level                Settings
    // 0                    post first md_stage
    // 1                    post last md_stage
    if (enc_mode <= ENC_MR) {
        context_ptr->chroma_at_last_md_stage = 0;
        context_ptr->chroma_at_last_md_stage_intra_th = (uint64_t)~0;
        context_ptr->chroma_at_last_md_stage_cfl_th = (uint64_t)~0;
    }
    else if (enc_mode <= ENC_M3) {
        context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
        context_ptr->chroma_at_last_md_stage_intra_th = 130;
        context_ptr->chroma_at_last_md_stage_cfl_th = 130;
    }
    else {
        context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
        context_ptr->chroma_at_last_md_stage_intra_th = 100;
        context_ptr->chroma_at_last_md_stage_cfl_th = 100;
    }
    // Cfl level
    // Level                Settings
    // 0                    Allow cfl
    // 1                    Disable cfl
    context_ptr->md_disable_cfl = EB_FALSE;
     // Set disallow_4x4
     if (enc_mode <= ENC_M1)
         context_ptr->disallow_4x4 = EB_FALSE;
     else
         context_ptr->disallow_4x4 = (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
     // If SB non-multiple of 4, then disallow_4x4 could not be used
     // SB Stats
     uint32_t sb_width =
         MIN(sequence_control_set_ptr->sb_size_pix, pcs_ptr->parent_pcs_ptr->aligned_width - context_ptr->sb_ptr->origin_x);
     uint32_t sb_height =
         MIN(sequence_control_set_ptr->sb_size_pix, pcs_ptr->parent_pcs_ptr->aligned_height - context_ptr->sb_ptr->origin_y);
     if (sb_width % 8 != 0 || sb_height % 8 != 0) {
         context_ptr->disallow_4x4 = EB_FALSE;
     }

     if (pd_pass == PD_PASS_0)
         context_ptr->md_disallow_nsq = enc_mode <= ENC_MR ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
     else if (pd_pass == PD_PASS_1)
         context_ptr->md_disallow_nsq = enc_mode <= ENC_MR ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
     else
         // Update nsq settings based on the sb_class
         context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;


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
            if (enc_mode <= ENC_M6)
                context_ptr->global_mv_injection = 1;
            else
                context_ptr->global_mv_injection = 0;
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
                if (enc_mode <= ENC_M1)
                    context_ptr->new_nearest_near_comb_injection = 1;
                else
                    context_ptr->new_nearest_near_comb_injection = 0;
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
    {
        if (enc_mode <= ENC_M2)
            context_ptr->unipred3x3_injection = 1;
        else
            context_ptr->unipred3x3_injection = 0;
    }

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
    else if (sequence_control_set_ptr->static_config.bipred_3x3_inject == DEFAULT) {
        if (enc_mode <= ENC_M2)
            context_ptr->bipred3x3_injection = 1;
        else if (enc_mode <= ENC_M6)
            context_ptr->bipred3x3_injection = 2;
        else
            context_ptr->bipred3x3_injection = 0;
        }
    else{
        context_ptr->bipred3x3_injection =
        sequence_control_set_ptr->static_config.bipred_3x3_inject;
        }
        // Level   Settings
        // 0       OFF: No compound mode search : AVG only
        // 1       ON: Full - AVG/DIST/DIFF/WEDGE
        // 2       ON: Fast - Use AVG only for non-closest ref frames or ref frames with high distortion
    if (sequence_control_set_ptr->compound_mode) {
            if (sequence_control_set_ptr->static_config.compound_level == DEFAULT) {
                if (enc_mode <= ENC_M0)
                    context_ptr->inter_compound_mode = 1;
                else if (enc_mode <= ENC_M3)
                    context_ptr->inter_compound_mode = 2;
                else
                    context_ptr->inter_compound_mode = 0;
            }
            else {
                context_ptr->inter_compound_mode = sequence_control_set_ptr->static_config.compound_level;
            }
        }
    else
            context_ptr->inter_compound_mode = 0;
    if (pd_pass == PD_PASS_0) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    }
    else
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;

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

    // spatial_sse_full_loop_level | Default Encoder Settings            | Command Line Settings
    //             0               | OFF subject to possible constraints | OFF in PD_PASS_2
    //             1               | ON subject to possible constraints  | ON in PD_PASS_2
    if (pd_pass == PD_PASS_0)
        context_ptr->spatial_sse_full_loop_level = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->spatial_sse_full_loop_level = EB_FALSE;
    else if (sequence_control_set_ptr->static_config.spatial_sse_full_loop_level == DEFAULT)
        if (enc_mode <= ENC_M9)
            context_ptr->spatial_sse_full_loop_level = EB_TRUE;
        else
            context_ptr->spatial_sse_full_loop_level = EB_FALSE;
    else
        context_ptr->spatial_sse_full_loop_level =
        sequence_control_set_ptr->static_config.spatial_sse_full_loop_level;

    if (context_ptr->chroma_level <= CHROMA_MODE_1)
        context_ptr->blk_skip_decision = EB_TRUE;
    else
        context_ptr->blk_skip_decision = EB_FALSE;

    if (pd_pass == PD_PASS_0)
        context_ptr->rdoq_level = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->rdoq_level = EB_FALSE;
    else
        if (sequence_control_set_ptr->static_config.rdoq_level == DEFAULT)
            if (enc_mode <= ENC_M9)
                context_ptr->rdoq_level = EB_TRUE;
            else
                context_ptr->rdoq_level = EB_FALSE;
        else
            context_ptr->rdoq_level =
            sequence_control_set_ptr->static_config.rdoq_level;

    // Derive redundant block
    if (pd_pass == PD_PASS_0)
        context_ptr->redundant_blk = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->redundant_blk = EB_TRUE;
    else
        if (sequence_control_set_ptr->static_config.enable_redundant_blk ==
            DEFAULT)
            if (enc_mode <= ENC_M9)
                context_ptr->redundant_blk = EB_TRUE;
            else
                context_ptr->redundant_blk = EB_FALSE;
        else
            context_ptr->redundant_blk =
            sequence_control_set_ptr->static_config.enable_redundant_blk;

    // md_stage_1_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_1_cand_prune_th
    if (pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_stage_1_cand_prune_th = 75;
    else
        if (enc_mode <= ENC_M5)
            context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
        else
            context_ptr->md_stage_1_cand_prune_th = 45;

    // md_stage_1_class_prune_th (for class removal)
    // Remove class if deviation to the best higher than TH_C
    if (pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_stage_1_class_prune_th = 100;
    else
        if (enc_mode <= ENC_M5)
            context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
        else
            context_ptr->md_stage_1_class_prune_th = 100;

    // md_stage_2_3_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than
    // md_stage_2_3_cand_prune_th
    if (pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_3_cand_prune_th = 5;
    else
        if (enc_mode <= ENC_MRS)
            context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
        else
        if (enc_mode <= ENC_MR)
            context_ptr->md_stage_2_3_cand_prune_th = 45;
        else if (enc_mode <= ENC_M9)
            context_ptr->md_stage_2_3_cand_prune_th = 15;
        else
            context_ptr->md_stage_2_3_cand_prune_th = 5;
    // md_stage_2_3_class_prune_th (for class removal)
    // Remove class if deviation to the best is higher than
    // md_stage_2_3_class_prune_th
    if (pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_3_class_prune_th = 25;
    else
        context_ptr->md_stage_2_3_class_prune_th = 25;
    // If using a mode offset, do not modify the NSQ-targeting features
    if (!mode_offset) {
        if (pd_pass == PD_PASS_0)
            context_ptr->coeff_area_based_bypass_nsq_th = 0;
        else if (pd_pass == PD_PASS_1)
            context_ptr->coeff_area_based_bypass_nsq_th = 0;
        else
            context_ptr->coeff_area_based_bypass_nsq_th = context_ptr->enable_area_based_cycles_allocation ? nsq_cycles_reduction_th[context_ptr->sb_class] : 0;
        uint8_t adaptive_md_cycles_level = 0;
        if (pd_pass == PD_PASS_2) {
            if (enc_mode <= ENC_MR)
                adaptive_md_cycles_level = 0;
            else if (enc_mode <= ENC_M0)
                adaptive_md_cycles_level = pcs_ptr->slice_type == I_SLICE ? 0 : 1;
            else if (enc_mode <= ENC_M1)
                adaptive_md_cycles_level = pcs_ptr->slice_type == I_SLICE ? 0 : 2;
            else if (enc_mode <= ENC_M2)
                adaptive_md_cycles_level = pcs_ptr->slice_type == I_SLICE ? 0 : 3;
            else if (enc_mode <= ENC_M3)
                adaptive_md_cycles_level = pcs_ptr->slice_type == I_SLICE ? 0 : 5;
            else
                adaptive_md_cycles_level = pcs_ptr->slice_type == I_SLICE ? 4 : 6;
        }
        adaptive_md_cycles_redcution_controls(context_ptr, adaptive_md_cycles_level);
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
            if (enc_mode <= ENC_MRS)
                context_ptr->sq_weight = (uint32_t)~0;
            else
                if (enc_mode <= ENC_MR)
                    context_ptr->sq_weight = 115;
                else
                    if (enc_mode <= ENC_M0)
                        context_ptr->sq_weight = 105;
                    else if (enc_mode <= ENC_M1)
                        context_ptr->sq_weight = 100;
                    else if (enc_mode <= ENC_M2)
                        context_ptr->sq_weight = 95;
                    else
                        context_ptr->sq_weight = 90;

    }
    // If using a mode offset, do not modify the NSQ-targeting features
    if (!mode_offset) {
        if (pd_pass < PD_PASS_2)
            context_ptr->switch_md_mode_based_on_sq_coeff = 0;
        else if (pcs_ptr->slice_type == I_SLICE)
            context_ptr->switch_md_mode_based_on_sq_coeff = 0;
        else if (enc_mode <= ENC_MR)
            context_ptr->switch_md_mode_based_on_sq_coeff = 0;
        else if (enc_mode <= ENC_M2)
            context_ptr->switch_md_mode_based_on_sq_coeff = 1;
        else if (enc_mode <= ENC_M3)
            context_ptr->switch_md_mode_based_on_sq_coeff = 2;
        else
            context_ptr->switch_md_mode_based_on_sq_coeff = 3;


        coeff_based_switch_md_controls(context_ptr, context_ptr->switch_md_mode_based_on_sq_coeff);

    }
    // Set pic_obmc_level @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_pic_obmc_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_pic_obmc_level = 0;
    else
        context_ptr->md_pic_obmc_level =
        pcs_ptr->parent_pcs_ptr->pic_obmc_level;

    set_obmc_controls(context_ptr, context_ptr->md_pic_obmc_level);

    // Set enable_inter_intra @ MD
#if  CLEANUP_INTER_INTRA
    //Block level switch, has to follow the picture level
#endif
    // inter intra pred                      Settings
    // 0                                     OFF
    // 1                                     FULL
    // 2                                     FAST 1 : Do not inject for unipred3x3 or PME inter candidates
    // 3                                     FAST 2 : Level 1 + do not inject for non-closest ref frames or ref frames with high distortion
    if (pcs_ptr->parent_pcs_ptr->slice_type != I_SLICE && sequence_control_set_ptr->seq_header.enable_interintra_compound) {
        if (pd_pass == PD_PASS_0)
            context_ptr->md_inter_intra_level = 0;
        else if (pd_pass == PD_PASS_1)
            context_ptr->md_inter_intra_level = 0;
        else if (enc_mode <= ENC_M2)
            context_ptr->md_inter_intra_level = 2;
        else
            context_ptr->md_inter_intra_level = 0;
    }
    else
        context_ptr->md_inter_intra_level = 0;

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

    // Assign whether to use TXS in inter classes (if TXS is ON)
    // 0 OFF - Use TXS for intra candidates only
    // 1 ON  - Use TXS for all candidates
    // 2 ON  - INTER TXS restricted to max 1 depth
    if (enc_mode <= ENC_MRS)
        context_ptr->txs_in_inter_classes = 1;
    else if (enc_mode <= ENC_M0)
        context_ptr->txs_in_inter_classes = 2;
    else
        context_ptr->txs_in_inter_classes = 0;

    // Each NIC scaling level corresponds to a scaling factor, given by the below {x,y}
    // combinations, where x is the numerator, and y is the denominator.  e.g. {1,8} corresponds
    // to 1/8x scaling of the base NICs, which are set in set_md_stage_counts().
    //{10,8 },    // level0
    //{ 8,8 },    // level1
    //{ 7,8 },    // level2
    //{ 6,8 },    // level3
    //{ 5,8 },    // level4
    //{ 4,8 },    // level5
    //{ 3,8 },    // level6
    //{ 2,8 },    // level7
    //{ 3,16},    // level8
    //{ 1,8 },    // level9
    //{ 1,16}     // level10
    // If using a mode offset, do not modify the NSQ-targeting features or NICS
    if (!mode_offset) {
        if (enc_mode <= ENC_MR)
            context_ptr->nic_scaling_level = 0;
        else if (enc_mode <= ENC_M0)
            context_ptr->nic_scaling_level = 1;
        else if (enc_mode <= ENC_M1)
            context_ptr->nic_scaling_level = 3;
        else if (enc_mode <= ENC_M2)
            context_ptr->nic_scaling_level = 6;
        else if (enc_mode <= ENC_M4)
            context_ptr->nic_scaling_level = 8;
        else
            context_ptr->nic_scaling_level = 9;
    }
    uint8_t txs_cycles_reduction_level = 0;
    set_txs_cycle_reduction_controls(context_ptr, txs_cycles_reduction_level);
    // Set md_filter_intra_mode @ MD
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
    // Set md_allow_intrabc @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_allow_intrabc = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_allow_intrabc = 0;
    else
        context_ptr->md_allow_intrabc = pcs_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc;

    // Set md_palette_level @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_palette_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_palette_level = 0;
    else
        context_ptr->md_palette_level = pcs_ptr->parent_pcs_ptr->palette_level;

    // Set block_based_depth_refinement_level
    if (enc_mode <= ENC_M6)
        context_ptr->block_based_depth_refinement_level = 0;
    else {
        if (pcs_ptr->slice_type == I_SLICE) {
            context_ptr->block_based_depth_refinement_level = 0;
        }
        else {
            context_ptr->block_based_depth_refinement_level = 1;
        }
    }
    set_block_based_depth_refinement_controls(context_ptr, context_ptr->block_based_depth_refinement_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->md_sq_mv_search_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_sq_mv_search_level = 0;
    else
        if (enc_mode <= ENC_M3)
            context_ptr->md_sq_mv_search_level = 1;
        else if (enc_mode <= ENC_M4)
            context_ptr->md_sq_mv_search_level = 2;
        else if (enc_mode <= ENC_M5)
            context_ptr->md_sq_mv_search_level = 3;
        else
            context_ptr->md_sq_mv_search_level = 4;
    md_sq_motion_search_controls(context_ptr, context_ptr->md_sq_mv_search_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->md_nsq_mv_search_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_nsq_mv_search_level = 0;
    else
        if (enc_mode <= ENC_MRS)
            context_ptr->md_nsq_mv_search_level = 1;
        else if (enc_mode <= ENC_M3)
            context_ptr->md_nsq_mv_search_level = 2;
        else
            context_ptr->md_nsq_mv_search_level = 4;

    md_nsq_motion_search_controls(context_ptr, context_ptr->md_nsq_mv_search_level);
    // Set PME level
    if (pd_pass == PD_PASS_0)
        context_ptr->md_pme_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_pme_level = 3;
    else
        if (enc_mode <= ENC_M2)
            context_ptr->md_pme_level = 1;
        else if (enc_mode <= ENC_M6)
            context_ptr->md_pme_level = 2;
        else
            context_ptr->md_pme_level = 3;
    md_pme_search_controls(context_ptr, context_ptr->md_pme_level);

    if (pd_pass == PD_PASS_0)
        context_ptr->md_subpel_me_level = enc_mode <= ENC_M4 ? 3 : 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_subpel_me_level = 3;
    else
        if (enc_mode <= ENC_M4)
            context_ptr->md_subpel_me_level = 1;
        else
            context_ptr->md_subpel_me_level = 2;

    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);

    if (pd_pass == PD_PASS_0)
        context_ptr->md_subpel_pme_level = enc_mode <= ENC_M4 ? 3 : 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_subpel_pme_level = 3;
    else
        if (enc_mode <= ENC_M4)
            context_ptr->md_subpel_pme_level = 1;
        else
            context_ptr->md_subpel_pme_level = 2;

    md_subpel_pme_controls(context_ptr, context_ptr->md_subpel_pme_level);
    // Set max_ref_count @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_max_ref_count = 4;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_max_ref_count = 1;
    else
        context_ptr->md_max_ref_count = 4;
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
    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    if (pd_pass == PD_PASS_0)
        context_ptr->shut_fast_rate = EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->shut_fast_rate = EB_FALSE;
    else
        context_ptr->shut_fast_rate = EB_FALSE;
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->skip_intra = 0;
    else if (pd_pass == PD_PASS_0)
        if (enc_mode <= ENC_M4)
            context_ptr->skip_intra = 0;
        else
            context_ptr->skip_intra = 1;
    else
        context_ptr->skip_intra = 0;
    // skip cfl based on inter/intra cost deviation (skip if intra_cost is
    // skip_cfl_cost_dev_th % greater than inter_cost)
    if (enc_mode <= ENC_MR)
        context_ptr->skip_cfl_cost_dev_th = (uint16_t)~0;
    else
        context_ptr->skip_cfl_cost_dev_th = 30;

    // set intra count to zero for md stage 3 if intra_cost is
    // mds3_intra_prune_th % greater than inter_cost
    if (enc_mode <= ENC_MR)
        context_ptr->mds3_intra_prune_th = (uint16_t)~0;
    else
        context_ptr->mds3_intra_prune_th = 30;

    return return_error;
}
/******************************************************
* Derive EncDec Settings for first pass
Input   : encoder mode and pd pass
Output  : EncDec Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_enc_dec_kernel(
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr);
void copy_neighbour_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                           uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds, uint32_t sb_org_x,
                           uint32_t sb_org_y);

static void set_parent_to_be_considered(MdcSbData *results_ptr, uint32_t blk_index, int32_t sb_size,
                                        int8_t pred_depth,
                                        uint8_t pred_sq_idx,
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
            results_ptr->leaf_data_array[parent_depth_idx_mds + block_1d_idx].pred_depth_refinement = parent_blk_geom->depth - pred_depth;
            results_ptr->leaf_data_array[parent_depth_idx_mds + block_1d_idx].pred_depth = pred_sq_idx;
            results_ptr->leaf_data_array[parent_depth_idx_mds + block_1d_idx].consider_block = 1;
        }

        if (depth_step < -1)
            set_parent_to_be_considered(results_ptr, parent_depth_idx_mds, sb_size, pred_depth, pred_sq_idx, depth_step + 1);
    }
}
static void set_child_to_be_considered(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr, MdcSbData *results_ptr, uint32_t blk_index, uint32_t sb_index, int32_t sb_size,
    int8_t pred_depth,
    uint8_t pred_sq_idx,
    int8_t depth_step) {


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
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].pred_depth_refinement = child1_blk_geom->depth - pred_depth;
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].pred_depth = pred_sq_idx;
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_1])
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_1, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
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
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].pred_depth_refinement = child2_blk_geom->depth - pred_depth;
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].pred_depth = pred_sq_idx;
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_2])
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_2, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
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
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].pred_depth_refinement = child3_blk_geom->depth - pred_depth;
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].pred_depth = pred_sq_idx;
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }

        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_3])
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_3, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
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
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].pred_depth_refinement = child4_blk_geom->depth - pred_depth;
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].pred_depth = pred_sq_idx;
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_4])
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_4, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
    }
}
static void build_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
    uint32_t sb_index) {

    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
    results_ptr->leaf_count = 0;
    uint32_t blk_index = 0;
    uint32_t d1_blocks_accumlated, d1_block_idx;

    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

        // SQ/NSQ block(s) filter based on the SQ size
        uint8_t is_block_tagged =
            (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE)
            ? 0
            : 1;

        // split_flag is f(min_sq_size)
        int32_t min_sq_size = (context_ptr->disallow_4x4) ? 8 : 4;

        // SQ/NSQ block(s) filter based on the block validity
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_block_tagged) {
            uint32_t tot_d1_blocks = (context_ptr->md_disallow_nsq) ||
                (blk_geom->sq_size >= 64 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_64x64) ||
                (blk_geom->sq_size >= 32 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_32x32) ||
                (blk_geom->sq_size >= 16 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_16x16) ||
                (blk_geom->sq_size <= 64 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_64x64) ||
                (blk_geom->sq_size <= 32 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_32x32) ||
                (blk_geom->sq_size <= 8 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_8x8) ||
                (blk_geom->sq_size <= 16 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_16x16) ? 1 :
                (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_all_non_hv_nsq_blocks_below_16x16) ? 5 :
                (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_all_h4_v4_blocks_below_16x16) ? 17 :
                blk_geom->sq_size == 128
                ? 17
                : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;

            if (pcs_ptr->parent_pcs_ptr->disallow_HVA_HVB_HV4)
                tot_d1_blocks = MIN(5, tot_d1_blocks);

            if (pcs_ptr->parent_pcs_ptr->disallow_HV4)
                tot_d1_blocks = MIN(17, tot_d1_blocks);
            d1_blocks_accumlated = 0;
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

            for (d1_block_idx = 0; d1_block_idx < tot_d1_blocks; d1_block_idx++)
                d1_blocks_accumlated +=
                results_ptr->leaf_data_array[blk_index + d1_block_idx].consider_block ? 1 : 0;

            for (uint32_t idx = 0; idx < tot_d1_blocks; ++idx) {
                if (results_ptr->leaf_data_array[blk_index].consider_block) {

                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = blk_index;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = d1_blocks_accumlated;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].final_pred_depth_refinement = results_ptr->leaf_data_array[blk_index].pred_depth_refinement;
                    if (results_ptr->leaf_data_array[results_ptr->leaf_count].final_pred_depth_refinement == -8)
                        printf("final_pred_depth_refinement error\n");
                    results_ptr->leaf_data_array[results_ptr->leaf_count].final_pred_depth = results_ptr->leaf_data_array[blk_index].pred_depth;
                    if (results_ptr->leaf_data_array[results_ptr->leaf_count].final_pred_depth == -8)
                        printf("final_pred_depth error\n");

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
void generate_statistics_txt(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t blk_index = 0;
    uint32_t part_cnt[TXT_DEPTH_DELTA_NUM][TX_TYPES] = { {0},{0} };
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        EbBool split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
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
                            for (uint8_t txb_itr = 0; txb_itr < best_blk_geom->txb_count[tx_depth]; txb_itr++) {
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
void generate_statistics_depth(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t blk_index = 0;
    // init stat
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        EbBool split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
        if (scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    if (context_ptr->md_local_blk_unit[blk_index].avail_blk_flag) {
                        int8_t pred_depth_refinement = context_ptr->md_local_blk_unit[blk_geom->sqi_mds].pred_depth_refinement;
                        pred_depth_refinement = MIN(pred_depth_refinement, 1);
                        pred_depth_refinement = MAX(pred_depth_refinement, -1);
                        uint8_t part_idx = part_to_shape[context_ptr->md_blk_arr_nsq[blk_index].part];
                        context_ptr->pred_depth_count[pred_depth_refinement + 2][part_idx]+= (blk_geom->bwidth*blk_geom->bheight);
                    }
                }
            }
        }
        blk_index += split_flag ?
            d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] :
            ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
}
/******************************************************
* Generate probabilities for the depth_cycles_reduction
******************************************************/
void generate_depth_prob(PictureControlSet * pcs_ptr, ModeDecisionContext *context_ptr)
{
    if (pcs_ptr->parent_pcs_ptr->slice_type != I_SLICE) {
        uint32_t pred_depth_count[DEPTH_DELTA_NUM][NUMBER_OF_SHAPES - 1] = { {0},{0},{0},{0},{0} };
        uint32_t samples_num = 0;
        // Sum statistics from reference list0
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list0_count_try; ref_idx++) {
            EbReferenceObject *ref_obj_l0 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][ref_idx]->object_ptr;
            for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
                for (uint8_t part_idx = 0; part_idx < (NUMBER_OF_SHAPES - 1); part_idx++) {
                    pred_depth_count[pred_depth][part_idx] += ref_obj_l0->ref_pred_depth_count[pred_depth][part_idx];
                    samples_num += ref_obj_l0->ref_pred_depth_count[pred_depth][part_idx];
                }
            }
        }
        // Sum statistics from reference list1
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list1_count_try; ref_idx++) {
            EbReferenceObject *ref_obj_l1 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx]->object_ptr;
            for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
                for (uint8_t part_idx = 0; part_idx < (NUMBER_OF_SHAPES - 1); part_idx++) {
                    pred_depth_count[pred_depth][part_idx] += ref_obj_l1->ref_pred_depth_count[pred_depth][part_idx];
                    samples_num += ref_obj_l1->ref_pred_depth_count[pred_depth][part_idx];
                }
            }
        }
        // Generate the selection %
        assert(samples_num > 0);
        for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
            for (uint8_t part_idx = 0; part_idx < (NUMBER_OF_SHAPES - 1); part_idx++) {
                context_ptr->ad_md_prob[pred_depth][part_idx] = (uint32_t)((pred_depth_count[pred_depth][part_idx] * (uint32_t)DEPTH_PROB_PRECISION) / (uint32_t)samples_num);
            }
        }
        //Generate depth prob
        for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
            for (uint8_t part_idx = 1; part_idx < (NUMBER_OF_SHAPES - 1); part_idx++) {
                pred_depth_count[pred_depth][0] += pred_depth_count[pred_depth][part_idx];
            }
        }
        for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++) {
            context_ptr->depth_prob[pred_depth] = (uint32_t)((pred_depth_count[pred_depth][0] * (uint32_t)100) / (uint32_t)samples_num);
        }

    }
    else {
        memcpy(context_ptr->ad_md_prob, intra_adaptive_md_cycles_reduction_th, sizeof(uint32_t) * DEPTH_DELTA_NUM * (NUMBER_OF_SHAPES - 1));
    }
}
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
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list0_count_try; ref_idx++) {
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
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list1_count_try; ref_idx++) {
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
void generate_statistics_nsq(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t blk_index = 0;
    uint32_t part_cnt[NUMBER_OF_SHAPES-1][FB_NUM][SSEG_NUM];
    for (uint8_t partidx = 0; partidx < NUMBER_OF_SHAPES-1; partidx++) {
        for (uint8_t band = 0; band < FB_NUM; band++) {
            memset(part_cnt[partidx][band], 0, sizeof(uint32_t) * SSEG_NUM);
        }
    }
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        EbBool split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
        if (scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    if (context_ptr->md_local_blk_unit[blk_index].avail_blk_flag) {
                        uint8_t band_idx = 0;
                        uint8_t sq_size_idx = 7 - (uint8_t)eb_log2f((uint8_t)blk_geom->sq_size);
                        uint64_t band_width = (sq_size_idx == 0) ? 100 : (sq_size_idx == 1) ? 50 : 20;
                        uint8_t part_idx = part_to_shape[context_ptr->md_blk_arr_nsq[blk_index].part];
                        uint8_t sse_g_band = (!context_ptr->md_disallow_nsq && context_ptr->md_local_blk_unit[blk_geom->sqi_mds].avail_blk_flag) ?
                            context_ptr->md_local_blk_unit[blk_geom->sqi_mds].sse_gradian_band[part_idx] : 1;
                        const uint32_t count_non_zero_coeffs = context_ptr->md_local_blk_unit[blk_index].count_non_zero_coeffs;
                        const uint32_t total_samples = (blk_geom->bwidth*blk_geom->bheight);
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
    for (uint8_t partidx = 0; partidx < NUMBER_OF_SHAPES-1; partidx++) {
        for (uint8_t band = 0; band < FB_NUM; band++) {
            for (uint8_t sse_idx = 0; sse_idx < SSEG_NUM; sse_idx++) {
                context_ptr->part_cnt[partidx][band][sse_idx] += part_cnt[partidx][band][sse_idx];
            }
        }
    }
}

/******************************************************
* Generate probabilities for the txt_cycles_reduction
******************************************************/
void generate_txt_prob(PictureControlSet * pcs_ptr,ModeDecisionContext *context_ptr)
{
    if (pcs_ptr->parent_pcs_ptr->slice_type != I_SLICE) {
        uint32_t txt_cnt[TXT_DEPTH_DELTA_NUM][TX_TYPES] = { {0},{0} };
        uint32_t samples_num = 0;
        // Sum statistics from reference list0
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list0_count_try; ref_idx++) {
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
        for (uint8_t ref_idx = 0; ref_idx < pcs_ptr->parent_pcs_ptr->ref_list1_count_try; ref_idx++) {
            EbReferenceObject *ref_obj_l1 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][ref_idx]->object_ptr;
            for (uint8_t depth_delta = 0; depth_delta < TXT_DEPTH_DELTA_NUM; depth_delta++) {
                for (uint8_t txs_idx = 1; txs_idx < TX_TYPES; txs_idx++) {
                    txt_cnt[depth_delta][txs_idx] += ref_obj_l1->ref_txt_cnt[depth_delta][txs_idx];
                    samples_num += ref_obj_l1->ref_txt_cnt[depth_delta][txs_idx];
                }
            }
        }
        assert(samples_num > 0);
        for (uint8_t depth_delta = 0; depth_delta < TXT_DEPTH_DELTA_NUM; depth_delta++) {
            for (uint8_t txs_idx = 1; txs_idx < TX_TYPES; txs_idx++) {
                context_ptr->txt_prob[depth_delta][txs_idx] = (uint32_t)((txt_cnt[depth_delta][txs_idx] * (uint32_t)10000) / (uint32_t)samples_num);
            }
        }
    }
}
const uint32_t sb_class_th[NUMBER_OF_SB_CLASS] = { 0,85,75,65,60,55,50,45,40,
                                                   35,30,25,20,17,14,10,6,3,0 };
static uint8_t determine_sb_class(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t blk_index = 0;
    uint64_t total_samples = 0;
    uint64_t count_non_zero_coeffs = 0;
    uint8_t sb_class = NONE_CLASS;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        EbBool split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
        if (scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            context_ptr->md_local_blk_unit[blk_index].avail_blk_flag &&
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
    return sb_class;
}
#define DEPTH_MAX_PROB 300 // max probabilty value for depth 100 -> 10%
// Depth probabilies per sq_size, pedicted depth and frequency band
// for sc content
uint8_t is_parent_to_current_deviation_small(SequenceControlSet *scs_ptr,
    ModeDecisionContext *mdctxt, const BlockGeom *blk_geom) {

    int64_t parent_to_current_deviation = MIN_SIGNED_VALUE;
    // block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
        // Get the parent of the current block
    uint32_t parent_depth_idx_mds =
        (blk_geom->sqi_mds -
        (blk_geom->quadi - 3) * ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]) -
        parent_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];

    if (mdctxt->md_local_blk_unit[parent_depth_idx_mds].avail_blk_flag) {
        parent_to_current_deviation =
            (int64_t)(((int64_t)MAX(mdctxt->md_local_blk_unit[parent_depth_idx_mds].default_cost, 1) - (int64_t)MAX((mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost * 4), 1)) * 100) /
            (int64_t)MAX((mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost * 4), 1);
    }

    if (parent_to_current_deviation <= mdctxt->depth_refinement_ctrls.parent_to_current_th)
        return EB_TRUE;

    return EB_FALSE;
}

uint8_t is_child_to_current_deviation_small(SequenceControlSet *scs_ptr,
    ModeDecisionContext *mdctxt, const BlockGeom *blk_geom, uint32_t blk_index) {

    int64_t child_to_current_deviation = MIN_SIGNED_VALUE;

    uint32_t child_block_idx_1, child_block_idx_2, child_block_idx_3, child_block_idx_4;
    child_block_idx_1 = blk_index + d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    child_block_idx_2 = child_block_idx_1 + ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth + 1];
    child_block_idx_3 = child_block_idx_2 + ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth + 1];
    child_block_idx_4 = child_block_idx_3 + ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth + 1];

    uint64_t child_cost = 0;
    uint8_t child_cnt = 0;
    if (mdctxt->md_local_blk_unit[child_block_idx_1].avail_blk_flag) {
        child_cost += mdctxt->md_local_blk_unit[child_block_idx_1].default_cost;
        child_cnt++;
    }
    if (mdctxt->md_local_blk_unit[child_block_idx_2].avail_blk_flag) {
        child_cost += mdctxt->md_local_blk_unit[child_block_idx_2].default_cost;
        child_cnt++;
    }
    if (mdctxt->md_local_blk_unit[child_block_idx_3].avail_blk_flag) {
        child_cost += mdctxt->md_local_blk_unit[child_block_idx_3].default_cost;
        child_cnt++;
    }
    if (mdctxt->md_local_blk_unit[child_block_idx_4].avail_blk_flag) {
        child_cost += mdctxt->md_local_blk_unit[child_block_idx_4].default_cost;
        child_cnt++;
    }

    if (child_cnt) {
        child_cost = (child_cost / child_cnt) * 4;
        child_to_current_deviation =
            (int64_t)(((int64_t)MAX(child_cost, 1) - (int64_t)MAX(mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost, 1)) * 100) /
            (int64_t)(MAX(mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost, 1));
    }


    if (child_to_current_deviation <= mdctxt->depth_refinement_ctrls.sub_to_current_th)
        return EB_TRUE;

    return EB_FALSE;
}
static void perform_pred_depth_refinement(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                          ModeDecisionContext *context_ptr, uint32_t sb_index) {
    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
    uint32_t   blk_index   = 0;

    // Reset mdc_sb_array data to defaults; it will be updated based on the predicted blocks (stored in md_blk_arr_nsq)
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom                              = get_blk_geom_mds(blk_index);
        results_ptr->leaf_data_array[blk_index].consider_block = 0;
        results_ptr->leaf_data_array[blk_index].split_flag =
            blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        results_ptr->leaf_data_array[blk_index].refined_split_flag =
            blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        results_ptr->leaf_data_array[blk_index].pred_depth_refinement = -8;
        results_ptr->leaf_data_array[blk_index].pred_depth = -8;
        blk_index++;
    }

    results_ptr->leaf_count = 0;
    blk_index               = 0;
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
                        // Shut thresholds in MR_MODE
                        if (pcs_ptr->enc_mode <= ENC_MRS) {
                            s_depth = -2;
                            e_depth = 2;
                        }
                        else if (pcs_ptr->enc_mode <= ENC_MR) {
                            if (pcs_ptr->parent_pcs_ptr->input_resolution == INPUT_SIZE_240p_RANGE) {
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
                        }
                        else if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0) {
                                if (pcs_ptr->enc_mode <= ENC_M5) {
                                s_depth = pcs_ptr->slice_type == I_SLICE ? -2 : -1;
                                e_depth = pcs_ptr->slice_type == I_SLICE ?  2 :  1;
                            }
                            else {
                                if (pcs_ptr->enc_mode <= ENC_M9) {
                                    s_depth = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? -1 : 0;
                                    e_depth = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 0;
                            }
                                else {
                                    s_depth = pcs_ptr->slice_type == I_SLICE ? -1 : 0;
                                    e_depth = pcs_ptr->slice_type == I_SLICE ? 1 : 0;
                                }
                            }
                        }
                        else {
                            s_depth = 0;
                            e_depth = 0;
                        }
                    } else if (context_ptr->pd_pass == PD_PASS_1) {
                        EbBool zero_coeff_present_flag =
                            context_ptr->md_blk_arr_nsq[blk_index].block_has_coeff == 0;
                        if (pcs_ptr->slice_type == I_SLICE) {
                            s_depth =  -1;
                            e_depth = (pcs_ptr->parent_pcs_ptr->enc_mode <= ENC_MR) ? 3 : 2;
                        }
                        else if (zero_coeff_present_flag && (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 || pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4)) {
                            s_depth = 0;
                            e_depth = 0;
                        } else

                            // This removes the SQ-versus-NSQ decision for the new MULTI_PASS_PD_LEVEL_1
                            if (pcs_ptr->parent_pcs_ptr->enc_mode <= ENC_MR || pcs_ptr->parent_pcs_ptr->multi_pass_pd_level <= MULTI_PASS_PD_LEVEL_1) {
                                s_depth = -1;
                                e_depth =  1;
                            }
                            else
                            if (context_ptr->md_local_blk_unit[blk_index].best_d1_blk == blk_index) {

                            s_depth = -1;
                                e_depth = 0;
                            } else {
                                s_depth = 0;
                            e_depth = 1;
                            }
                    }

                    // Check that the start and end depth are in allowed range, given other features
                    // which restrict allowable depths
                    if (context_ptr->disallow_4x4) {
                        e_depth = (blk_geom->sq_size == 8) ? 0
                                : (blk_geom->sq_size == 16) ? MIN(1, e_depth)
                                : (blk_geom->sq_size == 32) ? MIN(2, e_depth)
                                : e_depth;
                    }
                    // Add current pred depth block(s)
                    for (unsigned block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].consider_block = 1;
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].pred_depth_refinement = 0;
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].pred_depth = (int8_t)blk_geom->depth;
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].refined_split_flag =
                            EB_FALSE;
                    }

                    uint8_t sq_size_idx = 7 - (uint8_t)eb_log2f((uint8_t)blk_geom->sq_size);
                    // Add block indices of upper depth(s)
                    // Block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
                    uint8_t add_parent_depth = 1;
                    if (context_ptr->depth_refinement_ctrls.enabled && s_depth == -1 && pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[blk_index] && blk_geom->sq_size < ((scs_ptr->seq_header.sb_size == BLOCK_128X128) ? 128 : 64)) {
                        add_parent_depth = is_parent_to_current_deviation_small(
                            scs_ptr, context_ptr, blk_geom);
                    }
                    if (add_parent_depth)
                    if (s_depth != 0)
                        set_parent_to_be_considered(
                            results_ptr, blk_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth,sq_size_idx,  s_depth);
                    // Add block indices of lower depth(s)
                    // Block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
                    uint8_t add_sub_depth = 1;
                    if (context_ptr->depth_refinement_ctrls.enabled && e_depth == 1 && pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[blk_index]) {
                        add_sub_depth = is_child_to_current_deviation_small(
                            scs_ptr, context_ptr, blk_geom, blk_index);
                    }
                    if (add_sub_depth)
                    if (e_depth != 0)
                        set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, blk_index, sb_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth,sq_size_idx, e_depth);
                }
            }
        }
        blk_index +=
            split_flag
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
}
static void build_starting_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr, uint32_t sb_index) {

    MdcSbData *results_ptr = context_ptr->mdc_sb_array;

    results_ptr->leaf_count = 0;
    uint32_t blk_index = 0;
    int32_t force_blk_size = FORCED_BLK_SIZE;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
        if (use_output_stat(scs_ptr) && blk_geom->bheight >= FORCED_BLK_SIZE && blk_geom->bwidth >= FORCED_BLK_SIZE) {
            force_blk_size = FORCED_BLK_SIZE;
            if (blk_geom->bheight == FORCED_BLK_SIZE && blk_geom->bwidth == FORCED_BLK_SIZE &&
                !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index]) {
                int32_t cropped_width =
                    MIN(blk_geom->bwidth,
                        pcs_ptr->parent_pcs_ptr->aligned_width - (context_ptr->sb_origin_x + blk_geom->origin_x));
                int32_t cropped_height =
                    MIN(blk_geom->bheight,
                        pcs_ptr->parent_pcs_ptr->aligned_height - (context_ptr->sb_origin_y + blk_geom->origin_y));
                force_blk_size = (cropped_width != blk_geom->bwidth || cropped_height != blk_geom->bheight) ?
                    MAX(4, MIN(FORCED_BLK_SIZE, MIN(cropped_width, cropped_height))) :
                    FORCED_BLK_SIZE;
            }
        }
        // SQ/NSQ block(s) filter based on the SQ size
        uint8_t is_block_tagged =
            (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE) ||
            (blk_geom->sq_size == 4 && context_ptr->disallow_4x4)
            ? 0
            : 1;

        // split_flag is f(min_sq_size)
        int32_t min_sq_size = (context_ptr->disallow_4x4) ? 8 : 4;

        // SQ/NSQ block(s) filter based on the block validity
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_block_tagged) {
            uint32_t tot_d1_blocks = (context_ptr->md_disallow_nsq) ||
                (blk_geom->sq_size >= 64 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_64x64) ||
                (blk_geom->sq_size >= 32 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_32x32) ||
                (blk_geom->sq_size >= 16 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_above_16x16) ||
                (blk_geom->sq_size <= 64 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_64x64) ||
                (blk_geom->sq_size <= 32 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_32x32) ||
                (blk_geom->sq_size <= 8 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_8x8) ||
                (blk_geom->sq_size <= 16 && pcs_ptr->parent_pcs_ptr->disallow_all_nsq_blocks_below_16x16) ? 1 :
                (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_all_non_hv_nsq_blocks_below_16x16) ? 5 :
                (blk_geom->sq_size == 16 && pcs_ptr->parent_pcs_ptr->disallow_all_h4_v4_blocks_below_16x16) ? 17 :
                blk_geom->sq_size == 128
                ? 17
                : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;

            if (pcs_ptr->parent_pcs_ptr->disallow_HVA_HVB_HV4)
                tot_d1_blocks = MIN(5, tot_d1_blocks);

            if (pcs_ptr->parent_pcs_ptr->disallow_HV4)
                tot_d1_blocks = MIN(17, tot_d1_blocks);
            for (uint32_t idx = 0; idx < tot_d1_blocks; ++idx) {
                blk_geom = get_blk_geom_mds(blk_index);

                if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index]) {

                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = blk_index;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = tot_d1_blocks;

                    if (use_output_stat(scs_ptr)) {
                        if (blk_geom->sq_size == force_blk_size)
                            results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag = EB_FALSE;
                    }
                    else {
                    if (blk_geom->sq_size > min_sq_size)
                        results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag =
                        EB_TRUE;
                    else
                        results_ptr->leaf_data_array[results_ptr->leaf_count++].split_flag =
                        EB_FALSE;
                    }
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
void *mode_decision_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext *   thread_context_ptr = (EbThreadContext *)input_ptr;
    EncDecContext *     context_ptr        = (EncDecContext *)thread_context_ptr->priv;

    // Input
    EbObjectWrapper *enc_dec_tasks_wrapper_ptr;

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

        EncDecTasks *    enc_dec_tasks_ptr    = (EncDecTasks *)enc_dec_tasks_wrapper_ptr->object_ptr;
        PictureControlSet * pcs_ptr           = (PictureControlSet *)enc_dec_tasks_ptr->pcs_wrapper_ptr->object_ptr;
        SequenceControlSet *scs_ptr           = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
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

        memset(context_ptr->md_context->part_cnt, 0, sizeof(uint32_t) * SSEG_NUM * (NUMBER_OF_SHAPES-1) * FB_NUM);
        generate_nsq_prob(pcs_ptr, context_ptr->md_context);
        memset(context_ptr->md_context->pred_depth_count, 0, sizeof(uint32_t) * DEPTH_DELTA_NUM * (NUMBER_OF_SHAPES-1));
        generate_depth_prob(pcs_ptr, context_ptr->md_context);
        memset( context_ptr->md_context->txt_cnt, 0, sizeof(uint32_t) * TXT_DEPTH_DELTA_NUM * TX_TYPES);
        generate_txt_prob(pcs_ptr, context_ptr->md_context);

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
                    if (use_output_stat(scs_ptr) && sb_index == 0)
                        setup_firstpass_data(pcs_ptr->parent_pcs_ptr);
                    sb_ptr = context_ptr->md_context->sb_ptr = pcs_ptr->sb_ptr_array[sb_index];
                    sb_origin_x = (x_sb_index + tile_group_x_sb_start) << sb_size_log2;
                    sb_origin_y = (y_sb_index + tile_group_y_sb_start) << sb_size_log2;
                    //printf("[%ld]:ED sb index %d, (%d, %d), encoded total sb count %d, ctx coded sb count %d\n",
                    //        pcs_ptr->picture_number,
                    //        sb_index, sb_origin_x, sb_origin_y,
                    //        pcs_ptr->enc_dec_coded_sb_count,
                    //        context_ptr->coded_sb_count);
                    context_ptr->tile_index             = sb_ptr->tile_info.tile_rs_index;
                    context_ptr->md_context->tile_index = sb_ptr->tile_info.tile_rs_index;
                    context_ptr->md_context->sb_origin_x = sb_origin_x;
                    context_ptr->md_context->sb_origin_y = sb_origin_y;

                    sb_row_index_start =
                        (x_sb_index + 1 == tile_group_width_in_sb && sb_row_index_count == 0)
                            ? y_sb_index
                            : sb_row_index_start;
                    sb_row_index_count = (x_sb_index + 1 == tile_group_width_in_sb)
                                             ? sb_row_index_count + 1
                                             : sb_row_index_count;
                    mdc_ptr = context_ptr->md_context->mdc_sb_array;
                    context_ptr->sb_index = sb_index;
                    context_ptr->md_context->sb_class = NONE_CLASS;

                    if (pcs_ptr->update_cdf) {
                        if (scs_ptr->seq_header.pic_based_rate_est &&
                            scs_ptr->enc_dec_segment_row_count_array[pcs_ptr->temporal_layer_index] == 1 &&
                            scs_ptr->enc_dec_segment_col_count_array[pcs_ptr->temporal_layer_index] == 1) {
                            if (sb_index == 0)
                                pcs_ptr->ec_ctx_array[sb_index] =  pcs_ptr->md_frame_context;
                            else
                                pcs_ptr->ec_ctx_array[sb_index] = pcs_ptr->ec_ctx_array[sb_index - 1];
                        }
                        else {
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
                                  pcs_ptr->md_frame_context;
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
                            }
                        }
                        // Initial Rate Estimation of the syntax elements
                        av1_estimate_syntax_rate(&context_ptr->md_context->rate_est_table,
                            pcs_ptr->slice_type == I_SLICE,
                            &pcs_ptr->ec_ctx_array[sb_index]);
                        // Initial Rate Estimation of the Motion vectors
                        av1_estimate_mv_rate(pcs_ptr,
                            &context_ptr->md_context->rate_est_table,
                            &pcs_ptr->ec_ctx_array[sb_index]);

                        av1_estimate_coefficients_rate(&context_ptr->md_context->rate_est_table,
                            &pcs_ptr->ec_ctx_array[sb_index]);

                        //let the candidate point to the new rate table.
                        uint32_t cand_index;
                        for (cand_index = 0; cand_index < MODE_DECISION_CANDIDATE_MAX_COUNT;
                            ++cand_index)
                            context_ptr->md_context->fast_candidate_ptr_array[cand_index]
                            ->md_rate_estimation_ptr = &context_ptr->md_context->rate_est_table;
                        context_ptr->md_context->md_rate_estimation_ptr =
                            &context_ptr->md_context->rate_est_table;
                    }
                    // Configure the SB
                    mode_decision_configure_sb(
                        context_ptr->md_context, pcs_ptr, (uint8_t)sb_ptr->qindex);
                    // Multi-Pass PD
                    if ((pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_1 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 ||
                         pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4)
                        ) {
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
                            scs_ptr, pcs_ptr, context_ptr->md_context, 0);

                        // [PD_PASS_0]
                        // Input : mdc_blk_ptr built @ mdc process (up to 4421)
                        // Output: md_blk_arr_nsq reduced set of block(s)

                        // Build the t=0 cand_block_array
                        build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

                        // PD0 MD Tool(s) : ME_MV(s) as INTER candidate(s), DC as INTRA candidate, luma only, Frequency domain SSE,
                        // no fast rate (no MVP table generation), MDS0 then MDS3, reduced NIC(s), 1 ref per list,..
                        mode_decision_sb(scs_ptr,
                                         pcs_ptr,
                                         mdc_ptr,
                                         sb_ptr,
                                         sb_origin_x,
                                         sb_origin_y,
                                         sb_index,
                                         context_ptr->md_context);
                            context_ptr->md_context->sb_class = determine_sb_class(
                                scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

                        // Perform Pred_0 depth refinement - add depth(s) to be considered in the next stage(s)
                        perform_pred_depth_refinement(
                            scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

                        // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                        // Reset neighnor information to current SB @ position (0,0)
                        copy_neighbour_arrays(pcs_ptr,
                                              context_ptr->md_context,
                                              MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                              MD_NEIGHBOR_ARRAY_INDEX,
                                              0,
                                              sb_origin_x,
                                              sb_origin_y);

                        if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_1 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4) {
                            // [PD_PASS_1] Signal(s) derivation
                            context_ptr->md_context->pd_pass = PD_PASS_1;
                            signal_derivation_enc_dec_kernel_oq(
                                scs_ptr, pcs_ptr, context_ptr->md_context,0);
                            // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                            build_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

                            // [PD_PASS_1] Mode Decision - Further reduce the number of
                            // depth(s) to be considered in later PD stages. This pass uses more accurate
                            // info than PD0 to give a better PD estimate.
                            // Input : mdc_blk_ptr built @ PD0 refinement
                            // Output: md_blk_arr_nsq reduced set of block(s)

                            // PD1 MD Tool(s): PME,..
                            mode_decision_sb(scs_ptr,
                                             pcs_ptr,
                                             mdc_ptr,
                                             sb_ptr,
                                             sb_origin_x,
                                             sb_origin_y,
                                             sb_index,
                                             context_ptr->md_context);

                            // Perform Pred_1 depth refinement - add depth(s) to be considered in the next stage(s)
                            perform_pred_depth_refinement(
                                scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
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
                    // [PD_PASS_2] Signal(s) derivation
                    context_ptr->md_context->pd_pass = PD_PASS_2;
                    if (use_output_stat(scs_ptr))
                        first_pass_signal_derivation_enc_dec_kernel(pcs_ptr, context_ptr->md_context);
                    else
                        signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context, 0);
                    // Re-build mdc_blk_ptr for the 3rd PD Pass [PD_PASS_2]
                    if(pcs_ptr->parent_pcs_ptr->multi_pass_pd_level != MULTI_PASS_PD_OFF)
                    build_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
                    else
                        // Build the t=0 cand_block_array
                        build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

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
                    generate_statistics_nsq(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
                    generate_statistics_depth(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
                    generate_statistics_txt(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

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
                    if(!use_output_stat(scs_ptr))
                    av1_encode_decode(
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
        // Accumulate block selection
        for (uint8_t partidx = 0; partidx < NUMBER_OF_SHAPES-1; partidx++)
            for (uint8_t band = 0; band < FB_NUM; band++)
                for (uint8_t sse_idx = 0; sse_idx < SSEG_NUM; sse_idx++)
                    pcs_ptr->part_cnt[partidx][band][sse_idx] += context_ptr->md_context->part_cnt[partidx][band][sse_idx];

        // Accumulate pred depth selection
        for (uint8_t pred_depth = 0; pred_depth < DEPTH_DELTA_NUM; pred_depth++)
            for (uint8_t part_idx = 0; part_idx < (NUMBER_OF_SHAPES-1); part_idx++)
                pcs_ptr->pred_depth_count[pred_depth][part_idx] += context_ptr->md_context->pred_depth_count[pred_depth][part_idx];
        // Accumulate tx_type selection
        for (uint8_t depth_delta = 0; depth_delta < TXT_DEPTH_DELTA_NUM; depth_delta++)
            for (uint8_t txs_idx = 0; txs_idx < TX_TYPES; txs_idx++)
                pcs_ptr->txt_cnt[depth_delta][txs_idx] += context_ptr->md_context->txt_cnt[depth_delta][txs_idx];

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
            pcs_ptr->parent_pcs_ptr->av1x->rdmult =
                context_ptr->pic_full_lambda[(context_ptr->bit_depth == EB_10BIT) ? EB_10_BIT_MD
                                                                                  : EB_8_BIT_MD];
            if (use_output_stat(scs_ptr)) {
                first_pass_frame_end(pcs_ptr->parent_pcs_ptr, pcs_ptr->parent_pcs_ptr->ts_duration);
                if(pcs_ptr->parent_pcs_ptr->end_of_sequence_flag)
                    svt_av1_end_first_pass(pcs_ptr->parent_pcs_ptr);
            }
            eb_release_object(pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr);
            pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr = (EbObjectWrapper *)NULL;
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
