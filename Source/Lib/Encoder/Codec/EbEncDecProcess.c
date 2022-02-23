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
#include "EbPictureAnalysisProcess.h"
void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr, EbBool is_highbd);
int  svt_av1_allow_palette(int allow_palette, BlockSize sb_type);
#define FC_SKIP_TX_SR_TH025 125 // Fast cost skip tx search threshold.
#define FC_SKIP_TX_SR_TH010 110 // Fast cost skip tx search threshold.
void copy_mv_rate(PictureControlSet *pcs, MdRateEstimationContext *dst_rate);
void svt_av1_cdef_search(EncDecContext *context_ptr, SequenceControlSet *scs_ptr,
                         PictureControlSet *pcs_ptr);

void av1_cdef_frame16bit(uint8_t is_16bit, SequenceControlSet *scs_ptr, PictureControlSet *pCs);

void svt_av1_add_film_grain(EbPictureBufferDesc *src, EbPictureBufferDesc *dst,
                            AomFilmGrain *film_grain_ptr);

void svt_av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm,
                                                  int32_t after_cdef);
void svt_av1_loop_restoration_filter_frame(Yv12BufferConfig *frame, Av1Common *cm,
                                           int32_t optimized_lr);
extern void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr,
                          EbBool is_highbd);
void        svt_c_unpack_compressed_10bit(const uint8_t *inn_bit_buffer, uint32_t inn_stride,
                                          uint8_t *in_compn_bit_buffer, uint32_t out_stride,
                                          uint32_t height);

/*
* return by-pass encdec
*/
uint8_t get_bypass_encdec(EbEncMode enc_mode, uint8_t hbd_mode_decision,
                          uint8_t encoder_bit_depth) {
    UNUSED(hbd_mode_decision);
    uint8_t bypass_encdec = 1;
    if (encoder_bit_depth == EB_8BIT) {
        // 8bit settings
        if (enc_mode <= ENC_M5)
            bypass_encdec = 0;
        else
            bypass_encdec = 1;
    } else {
        // 10bit settings
        if (enc_mode <= ENC_M9)
            bypass_encdec = 0;
        else
            bypass_encdec = 1;
    }
    return bypass_encdec;
}

static void enc_dec_context_dctor(EbPtr p) {
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    EncDecContext   *obj                = (EncDecContext *)thread_context_ptr->priv;
    EB_DELETE(obj->md_context);
    EB_DELETE(obj->residual_buffer);
    EB_DELETE(obj->transform_buffer);
    EB_DELETE(obj->inverse_quant_buffer);
    EB_DELETE(obj->input_sample16bit_buffer);
    EB_FREE_ARRAY(obj);
}

/******************************************************
 * Enc Dec Context Constructor
 ******************************************************/
EbErrorType enc_dec_context_ctor(EbThreadContext   *thread_context_ptr,
                                 const EbEncHandle *enc_handle_ptr, int index, int tasks_index)

{
    const EbSvtAv1EncConfiguration *static_config =
        &enc_handle_ptr->scs_instance_array[0]->scs_ptr->static_config;
    EbColorFormat color_format = static_config->encoder_color_format;
    int8_t        enable_hbd_mode_decision =
        enc_handle_ptr->scs_instance_array[0]->scs_ptr->enable_hbd_mode_decision;

    EncDecContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = enc_dec_context_dctor;

    context_ptr->is_16bit     = enc_handle_ptr->scs_instance_array[0]->scs_ptr->is_16bit_pipeline;
    context_ptr->color_format = color_format;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->enc_dec_tasks_resource_ptr, index);
    context_ptr->enc_dec_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_results_resource_ptr, index);
    context_ptr->enc_dec_feedback_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_tasks_resource_ptr, tasks_index);

    // Prediction Buffer
    context_ptr->input_sample16bit_buffer = NULL;
    if (context_ptr->is_16bit)
        EB_NEW(context_ptr->input_sample16bit_buffer,
               svt_picture_buffer_desc_ctor,
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

    EB_NEW(context_ptr->inverse_quant_buffer, svt_picture_buffer_desc_ctor, &init_32bit_data);
    EB_NEW(context_ptr->transform_buffer, svt_picture_buffer_desc_ctor, &init_32bit_data);
    EB_NEW(context_ptr->residual_buffer,
           svt_picture_buffer_desc_ctor,
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
           enc_handle_ptr->scs_instance_array[0]->scs_ptr->super_block_size,
           static_config->enc_mode,
           enc_handle_ptr->scs_instance_array[0]->scs_ptr->max_block_cnt,
           static_config->encoder_bit_depth,
           0,
           0,
           enable_hbd_mode_decision == DEFAULT ? 2 : enable_hbd_mode_decision,
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
    neighbor_array_unit_reset(pcs_ptr->ep_luma_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_luma_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_luma_dc_sign_level_coeff_neighbor_array_update[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_dc_sign_level_coeff_neighbor_array_update[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_dc_sign_level_coeff_neighbor_array_update[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_partition_context_neighbor_array[tile_idx]);
    // TODO(Joel): 8-bit ep_luma_recon_neighbor_array (Cb,Cr) when is_16bit==0?
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->is_16bit_pipeline) {
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
    context_ptr->is_16bit   = scs_ptr->is_16bit_pipeline;
    context_ptr->bit_depth  = scs_ptr->static_config.encoder_bit_depth;
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
    EbBool   continue_processing_flag = EB_FALSE;
    uint32_t row_segment_index        = 0;
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
        // Reset enc_dec segments
        for (uint32_t row_index = 0; row_index < segmentPtr->segment_row_count; ++row_index) {
            segmentPtr->row_array[row_index].current_seg_index =
                segmentPtr->row_array[row_index].starting_seg_index;
        }

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
            svt_block_on_mutex(segmentPtr->row_array[row_segment_index].assignment_mutex);

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

            svt_release_mutex(segmentPtr->row_array[row_segment_index].assignment_mutex);
        }

        // Bottom-left Neighbor
        if (row_segment_index < segmentPtr->segment_row_count - 1 &&
            bottom_left_segment_index >=
                segmentPtr->row_array[row_segment_index + 1].starting_seg_index) {
            svt_block_on_mutex(segmentPtr->row_array[row_segment_index + 1].assignment_mutex);

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
            svt_release_mutex(segmentPtr->row_array[row_segment_index + 1].assignment_mutex);
        }

        if (feedback_row_index > 0) {
            EbObjectWrapper *wrapper_ptr;
            svt_get_empty_object(srmFifoPtr, &wrapper_ptr);
            EncDecTasks *feedback_task_ptr         = (EncDecTasks *)wrapper_ptr->object_ptr;
            feedback_task_ptr->input_type          = ENCDEC_TASKS_ENCDEC_INPUT;
            feedback_task_ptr->enc_dec_segment_row = feedback_row_index;
            feedback_task_ptr->pcs_wrapper_ptr     = taskPtr->pcs_wrapper_ptr;
            feedback_task_ptr->tile_group_index    = taskPtr->tile_group_index;
            svt_post_full_object(wrapper_ptr);
        }

        break;

    default: break;
    }

    return continue_processing_flag;
}
void recon_output(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
    // The totalNumberOfReconFrames counter has to be write/read protected as
    //   it is used to determine the end of the stream.  If it is not protected
    //   the encoder might not properly terminate.
    svt_block_on_mutex(encode_context_ptr->total_number_of_recon_frame_mutex);

    if (!pcs_ptr->parent_pcs_ptr->is_alt_ref) {
        EbBool           is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
        EbObjectWrapper *output_recon_wrapper_ptr;
        // Get Recon Buffer
        svt_get_empty_object(scs_ptr->encode_context_ptr->recon_output_fifo_ptr,
                             &output_recon_wrapper_ptr);
        EbBufferHeaderType *output_recon_ptr = (EbBufferHeaderType *)
                                                   output_recon_wrapper_ptr->object_ptr;
        output_recon_ptr->flags = 0;

        // START READ/WRITE PROTECTED SECTION
        if (encode_context_ptr->total_number_of_recon_frames ==
            encode_context_ptr->terminating_picture_number)
            output_recon_ptr->flags = EB_BUFFERFLAG_EOS;

        encode_context_ptr->total_number_of_recon_frames++;

        // STOP READ/WRITE PROTECTED SECTION
        output_recon_ptr->n_filled_len = 0;

        // Copy the Reconstructed Picture to the Output Recon Buffer
        {
            uint32_t sample_total_count;
            uint8_t *recon_read_ptr;
            uint8_t *recon_write_ptr;

            EbPictureBufferDesc *recon_ptr;
            EbPictureBufferDesc *intermediate_buffer_ptr = NULL;
            { get_recon_pic(pcs_ptr, &recon_ptr, is_16bit); }

            // FGN: Create a buffer if needed, copy the reconstructed picture and run the film grain synthesis algorithm
            if ((scs_ptr->static_config.pass != ENC_FIRST_PASS) &&
                scs_ptr->seq_header.film_grain_params_present &&
                pcs_ptr->parent_pcs_ptr->frm_hdr.film_grain_params.apply_grain) {
                AomFilmGrain *film_grain_ptr;

                uint16_t                    padding = scs_ptr->super_block_size + 32;
                EbPictureBufferDescInitData temp_recon_desc_init_data;
                temp_recon_desc_init_data.max_width  = (uint16_t)scs_ptr->max_input_luma_width;
                temp_recon_desc_init_data.max_height = (uint16_t)scs_ptr->max_input_luma_height;
                temp_recon_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

                temp_recon_desc_init_data.left_padding  = padding;
                temp_recon_desc_init_data.right_padding = padding;
                temp_recon_desc_init_data.top_padding   = padding;
                temp_recon_desc_init_data.bot_padding   = padding;
                temp_recon_desc_init_data.split_mode    = EB_FALSE;
                temp_recon_desc_init_data.color_format =
                    scs_ptr->static_config.encoder_color_format;

                if (is_16bit) {
                    temp_recon_desc_init_data.bit_depth = EB_16BIT;
                } else {
                    temp_recon_desc_init_data.bit_depth = EB_8BIT;
                }

                EB_NO_THROW_NEW(intermediate_buffer_ptr,
                                svt_recon_picture_buffer_desc_ctor,
                                (EbPtr)&temp_recon_desc_init_data);

                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    film_grain_ptr = &((EbReferenceObject *)pcs_ptr->parent_pcs_ptr
                                           ->reference_picture_wrapper_ptr->object_ptr)
                                          ->film_grain_params;
                else
                    film_grain_ptr = &pcs_ptr->parent_pcs_ptr->frm_hdr.film_grain_params;

                if (intermediate_buffer_ptr) {
                    svt_av1_add_film_grain(recon_ptr, intermediate_buffer_ptr, film_grain_ptr);
                    recon_ptr = intermediate_buffer_ptr;
                }
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

            if (intermediate_buffer_ptr) {
                EB_DELETE(intermediate_buffer_ptr);
            }
        }

        // Post the Recon object
        svt_post_full_object(output_recon_wrapper_ptr);
    } else {
        // Overlay and altref have 1 recon only, which is from overlay pictures. So the recon of the alt_ref is not sent to the application.
        // However, to hanlde the end of sequence properly, total_number_of_recon_frames is increamented
        encode_context_ptr->total_number_of_recon_frames++;
    }
    svt_release_mutex(encode_context_ptr->total_number_of_recon_frame_mutex);
}

//************************************/
// Calculate Frame SSIM
/************************************/

static void svt_aom_ssim_parms_8x8_c(const uint8_t *s, int sp, const uint8_t *r, int rp,
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

static void svt_aom_highbd_ssim_parms_8x8_c(const uint8_t *s, int sp, const uint8_t *sinc,
                                            int spinc, const uint16_t *r, int rp, uint32_t *sum_s,
                                            uint32_t *sum_r, uint32_t *sum_sq_s, uint32_t *sum_sq_r,
                                            uint32_t *sum_sxr) {
    int      i, j;
    uint32_t ss;
    for (i = 0; i < 8; i++, s += sp, sinc += spinc, r += rp) {
        for (j = 0; j < 8; j++) {
            ss = (int64_t)(s[j] << 2) + ((sinc[j] >> 6) & 0x3);
            *sum_s += ss;
            *sum_r += r[j];
            *sum_sq_s += ss * ss;
            *sum_sq_r += r[j] * r[j];
            *sum_sxr += ss * r[j];
        }
    }
}

static const int64_t cc1    = 26634; // (64^2*(.01*255)^2
static const int64_t cc2    = 239708; // (64^2*(.03*255)^2
static const int64_t cc1_10 = 428658; // (64^2*(.01*1023)^2
static const int64_t cc2_10 = 3857925; // (64^2*(.03*1023)^2
static const int64_t cc1_12 = 6868593; // (64^2*(.01*4095)^2
static const int64_t cc2_12 = 61817334; // (64^2*(.03*4095)^2

static double similarity(uint32_t sum_s, uint32_t sum_r, uint32_t sum_sq_s, uint32_t sum_sq_r,
                         uint32_t sum_sxr, int count, uint32_t bd) {
    double  ssim_n, ssim_d;
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

    ssim_n = (2.0 * sum_s * sum_r + c1) * (2.0 * count * sum_sxr - 2.0 * sum_s * sum_r + c2);

    ssim_d = ((double)sum_s * sum_s + (double)sum_r * sum_r + c1) *
        ((double)count * sum_sq_s - (double)sum_s * sum_s + (double)count * sum_sq_r -
         (double)sum_r * sum_r + c2);

    return ssim_n / ssim_d;
}

static double ssim_8x8(const uint8_t *s, int sp, const uint8_t *r, int rp) {
    uint32_t sum_s = 0, sum_r = 0, sum_sq_s = 0, sum_sq_r = 0, sum_sxr = 0;
    svt_aom_ssim_parms_8x8_c(s, sp, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
    return similarity(sum_s, sum_r, sum_sq_s, sum_sq_r, sum_sxr, 64, 8);
}

static double highbd_ssim_8x8(const uint8_t *s, int sp, const uint8_t *sinc, int spinc,
                              const uint16_t *r, int rp, uint32_t bd, uint32_t shift) {
    uint32_t sum_s = 0, sum_r = 0, sum_sq_s = 0, sum_sq_r = 0, sum_sxr = 0;
    svt_aom_highbd_ssim_parms_8x8_c(
        s, sp, sinc, spinc, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
    return similarity(sum_s >> shift,
                      sum_r >> shift,
                      sum_sq_s >> (2 * shift),
                      sum_sq_r >> (2 * shift),
                      sum_sxr >> (2 * shift),
                      64,
                      bd);
}

// We are using a 8x8 moving window with starting location of each 8x8 window
// on the 4x4 pixel grid. Such arrangement allows the windows to overlap
// block boundaries to penalize blocking artifacts.
static double aom_ssim2(const uint8_t *img1, int stride_img1, const uint8_t *img2, int stride_img2,
                        int width, int height) {
    int    i, j;
    int    samples    = 0;
    double ssim_total = 0;

    // sample point start with each 4x4 location
    for (i = 0; i <= height - 8; i += 4, img1 += stride_img1 * 4, img2 += stride_img2 * 4) {
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

static double aom_highbd_ssim2(const uint8_t *img1, int stride_img1, const uint8_t *img1inc,
                               int stride_img1inc, const uint16_t *img2, int stride_img2, int width,
                               int height, uint32_t bd, uint32_t shift) {
    int    i, j;
    int    samples    = 0;
    double ssim_total = 0;

    // sample point start with each 4x4 location
    for (i = 0; i <= height - 8;
         i += 4, img1 += stride_img1 * 4, img1inc += stride_img1inc * 4, img2 += stride_img2 * 4) {
        for (j = 0; j <= width - 8; j += 4) {
            double v = highbd_ssim_8x8((img1 + j),
                                       stride_img1,
                                       (img1inc + j),
                                       stride_img1inc,
                                       (img2 + j),
                                       stride_img2,
                                       bd,
                                       shift);
            ssim_total += v;
            samples++;
        }
    }
    assert(samples > 0);
    ssim_total /= samples;
    return ssim_total;
}

void free_temporal_filtering_buffer(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    // save_source_picture_ptr will be allocated only if temporal_filtering_on is true in svt_av1_init_temporal_filtering().
    if (!pcs_ptr->parent_pcs_ptr->temporal_filtering_on) {
        return;
    }

    EB_FREE_ARRAY(pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[0]);
    EB_FREE_ARRAY(pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[1]);
    EB_FREE_ARRAY(pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[2]);

    EbBool is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    if (is_16bit) {
        EB_FREE_ARRAY(pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[0]);
        EB_FREE_ARRAY(pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[1]);
        EB_FREE_ARRAY(pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[2]);
    }
}

EbErrorType ssim_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr,
                              EbBool free_memory) {
    EbBool is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);

    const uint32_t ss_x = scs_ptr->subsampling_x;
    const uint32_t ss_y = scs_ptr->subsampling_y;

    EbPictureBufferDesc *recon_ptr;
    get_recon_pic(pcs_ptr, &recon_ptr, is_16bit);
    if (!is_16bit) {
        EbPictureBufferDesc *input_picture_ptr =
            (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;

        EbByte input_buffer;
        EbByte recon_coeff_buffer;

        EbByte buffer_y;
        EbByte buffer_cb;
        EbByte buffer_cr;

        double luma_ssim = 0.0;
        double cb_ssim   = 0.0;
        double cr_ssim   = 0.0;

        // if current source picture was temporally filtered, use an alternative buffer which stores
        // the original source picture
        if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            assert(pcs_ptr->parent_pcs_ptr->save_source_picture_width == input_picture_ptr->width &&
                   pcs_ptr->parent_pcs_ptr->save_source_picture_height ==
                       input_picture_ptr->height);
            buffer_y  = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[0];
            buffer_cb = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[1];
            buffer_cr = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[2];
        } else {
            buffer_y  = input_picture_ptr->buffer_y;
            buffer_cb = input_picture_ptr->buffer_cb;
            buffer_cr = input_picture_ptr->buffer_cr;
        }

        recon_coeff_buffer = &(
            (recon_ptr->buffer_y)[recon_ptr->origin_x + recon_ptr->origin_y * recon_ptr->stride_y]);
        input_buffer = &(buffer_y[input_picture_ptr->origin_x +
                                  input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
        luma_ssim    = aom_ssim2(input_buffer,
                              input_picture_ptr->stride_y,
                              recon_coeff_buffer,
                              recon_ptr->stride_y,
                              scs_ptr->seq_header.max_frame_width,
                              scs_ptr->seq_header.max_frame_height);

        recon_coeff_buffer = &(
            (recon_ptr->buffer_cb)[recon_ptr->origin_x / 2 +
                                   recon_ptr->origin_y / 2 * recon_ptr->stride_cb]);
        input_buffer = &(buffer_cb[input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
        cb_ssim      = aom_ssim2(input_buffer,
                            input_picture_ptr->stride_cb,
                            recon_coeff_buffer,
                            recon_ptr->stride_cb,
                            scs_ptr->chroma_width,
                            scs_ptr->chroma_height);

        recon_coeff_buffer = &(
            (recon_ptr->buffer_cr)[recon_ptr->origin_x / 2 +
                                   recon_ptr->origin_y / 2 * recon_ptr->stride_cr]);
        input_buffer = &(buffer_cr[input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
        cr_ssim      = aom_ssim2(input_buffer,
                            input_picture_ptr->stride_cr,
                            recon_coeff_buffer,
                            recon_ptr->stride_cr,
                            scs_ptr->chroma_width,
                            scs_ptr->chroma_height);

        pcs_ptr->parent_pcs_ptr->luma_ssim = luma_ssim;
        pcs_ptr->parent_pcs_ptr->cb_ssim   = cb_ssim;
        pcs_ptr->parent_pcs_ptr->cr_ssim   = cr_ssim;

        if (free_memory && pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            EB_FREE_ARRAY(buffer_y);
            EB_FREE_ARRAY(buffer_cb);
            EB_FREE_ARRAY(buffer_cr);
        }
    } else {
        get_recon_pic(pcs_ptr, &recon_ptr, is_16bit);
        EbPictureBufferDesc *input_picture_ptr =
            (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;

        EbByte    input_buffer;
        uint16_t *recon_coeff_buffer;

        double luma_ssim = 0.0;
        double cb_ssim   = 0.0;
        double cr_ssim   = 0.0;

        if (scs_ptr->ten_bit_format == 1) {
            /* SSIM calculation for compressed 10-bit format has not been verified and debugged,
               since this format is not supported elsewhere in this version. See verify_settings(),
               which exits with an error if compressed 10-bit format is enabled. To avoid
               extra complexity of unpacking into a temporary buffer, or having to write
               new core SSIM functions, we ignore the two least signifcant bits in this
               case, and set these to zero. One test shows a difference in SSIM
               of 0.00085 setting the two least significant bits to zero. */

            const uint32_t luma_width   = input_picture_ptr->width - scs_ptr->max_input_pad_right;
            const uint32_t luma_height  = input_picture_ptr->height - scs_ptr->max_input_pad_bottom;
            const uint32_t chroma_width = luma_width >> ss_x;
            const uint32_t pic_width_in_sb  = (luma_width + 64 - 1) / 64;
            const uint32_t pic_height_in_sb = (luma_height + 64 - 1) / 64;
            const uint32_t chroma_height    = luma_height >> ss_y;
            uint32_t       sb_num_in_height, sb_num_in_width, bd, shift;
            uint8_t        zero_buffer[64 * 64];

            bd    = 10;
            shift = 0; // both input and output are 10 bit (bitdepth - input_bd)
            memset(&zero_buffer[0], 0, sizeof(uint8_t) * 64 * 64);

            EbByte input_buffer_org = &(
                (input_picture_ptr
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
            uint16_t *recon_buffer_org_u = recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr
                     ->buffer_cb)[(recon_ptr->origin_x << is_16bit) / 2 +
                                  (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cb]));
            ;

            EbByte input_buffer_org_v = &(
                (input_picture_ptr
                     ->buffer_cr)[input_picture_ptr->origin_x / 2 +
                                  input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            ;
            uint16_t *recon_buffer_org_v = recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr
                     ->buffer_cr)[(recon_ptr->origin_x << is_16bit) / 2 +
                                  (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cr]));
            ;

            for (sb_num_in_height = 0; sb_num_in_height < pic_height_in_sb; ++sb_num_in_height) {
                for (sb_num_in_width = 0; sb_num_in_width < pic_width_in_sb; ++sb_num_in_width) {
                    uint32_t tb_origin_x = sb_num_in_width * 64;
                    uint32_t tb_origin_y = sb_num_in_height * 64;
                    uint32_t sb_width = (luma_width - tb_origin_x) < 64 ? (luma_width - tb_origin_x)
                                                                        : 64;
                    uint32_t sb_height = (luma_height - tb_origin_y) < 64
                        ? (luma_height - tb_origin_y)
                        : 64;

                    input_buffer = input_buffer_org + tb_origin_y * input_picture_ptr->stride_y +
                        tb_origin_x;
                    recon_coeff_buffer = recon_buffer_org + tb_origin_y * recon_ptr->stride_y +
                        tb_origin_x;

                    luma_ssim += aom_highbd_ssim2(input_buffer,
                                                  input_picture_ptr->stride_y,
                                                  &zero_buffer[0],
                                                  64,
                                                  recon_coeff_buffer,
                                                  recon_ptr->stride_y,
                                                  sb_width,
                                                  sb_height,
                                                  bd,
                                                  shift);

                    //U+V
                    tb_origin_x = sb_num_in_width * 32;
                    tb_origin_y = sb_num_in_height * 32;
                    sb_width    = (chroma_width - tb_origin_x) < 32 ? (chroma_width - tb_origin_x)
                                                                    : 32;
                    sb_height   = (chroma_height - tb_origin_y) < 32 ? (chroma_height - tb_origin_y)
                                                                     : 32;

                    input_buffer = input_buffer_org_u + tb_origin_y * input_picture_ptr->stride_cb +
                        tb_origin_x;
                    recon_coeff_buffer = recon_buffer_org_u + tb_origin_y * recon_ptr->stride_cb +
                        tb_origin_x;

                    cb_ssim += aom_highbd_ssim2(input_buffer,
                                                input_picture_ptr->stride_cb,
                                                &zero_buffer[0],
                                                64,
                                                recon_coeff_buffer,
                                                recon_ptr->stride_cb,
                                                sb_width,
                                                sb_height,
                                                bd,
                                                shift);

                    input_buffer = input_buffer_org_v + tb_origin_y * input_picture_ptr->stride_cr +
                        tb_origin_x;
                    recon_coeff_buffer = recon_buffer_org_v + tb_origin_y * recon_ptr->stride_cr +
                        tb_origin_x;

                    cr_ssim += aom_highbd_ssim2(input_buffer,
                                                input_picture_ptr->stride_cr,
                                                &zero_buffer[0],
                                                64,
                                                recon_coeff_buffer,
                                                recon_ptr->stride_cr,
                                                sb_width,
                                                sb_height,
                                                bd,
                                                shift);
                }
            }

            luma_ssim /= pic_height_in_sb * pic_width_in_sb;
            cb_ssim /= pic_height_in_sb * pic_width_in_sb;
            cr_ssim /= pic_height_in_sb * pic_width_in_sb;

            pcs_ptr->parent_pcs_ptr->luma_ssim = luma_ssim;
            pcs_ptr->parent_pcs_ptr->cb_ssim   = cb_ssim;
            pcs_ptr->parent_pcs_ptr->cr_ssim   = cr_ssim;
        } else {
            recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr->buffer_y)[(recon_ptr->origin_x << is_16bit) +
                                      (recon_ptr->origin_y << is_16bit) * recon_ptr->stride_y]));

            // if current source picture was temporally filtered, use an alternative buffer which stores
            // the original source picture
            EbByte buffer_y, buffer_bit_inc_y;
            EbByte buffer_cb, buffer_bit_inc_cb;
            EbByte buffer_cr, buffer_bit_inc_cr;
            int    bd, shift;

            if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
                assert(pcs_ptr->parent_pcs_ptr->save_source_picture_width ==
                           input_picture_ptr->width &&
                       pcs_ptr->parent_pcs_ptr->save_source_picture_height ==
                           input_picture_ptr->height);
                buffer_y          = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[0];
                buffer_bit_inc_y  = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[0];
                buffer_cb         = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[1];
                buffer_bit_inc_cb = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[1];
                buffer_cr         = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[2];
                buffer_bit_inc_cr = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[2];
            } else {
                uint32_t height_y  = (uint32_t)(input_picture_ptr->height +
                                               input_picture_ptr->origin_y +
                                               input_picture_ptr->origin_bot_y);
                uint32_t height_uv = (uint32_t)((input_picture_ptr->height +
                                                 input_picture_ptr->origin_y +
                                                 input_picture_ptr->origin_bot_y) >>
                                                ss_y);

                uint8_t *uncompressed_pics[3];
                EB_MALLOC_ARRAY(uncompressed_pics[0],
                                pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr->luma_size);
                EB_MALLOC_ARRAY(
                    uncompressed_pics[1],
                    pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr->chroma_size);
                EB_MALLOC_ARRAY(
                    uncompressed_pics[2],
                    pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr->chroma_size);

                svt_c_unpack_compressed_10bit(input_picture_ptr->buffer_bit_inc_y,
                                              input_picture_ptr->stride_bit_inc_y / 4,
                                              uncompressed_pics[0],
                                              input_picture_ptr->stride_bit_inc_y,
                                              height_y);
                // U
                svt_c_unpack_compressed_10bit(input_picture_ptr->buffer_bit_inc_cb,
                                              input_picture_ptr->stride_bit_inc_cb / 4,
                                              uncompressed_pics[1],
                                              input_picture_ptr->stride_bit_inc_cb,
                                              height_uv);
                // V
                svt_c_unpack_compressed_10bit(input_picture_ptr->buffer_bit_inc_cr,
                                              input_picture_ptr->stride_bit_inc_cr / 4,
                                              uncompressed_pics[2],
                                              input_picture_ptr->stride_bit_inc_cr,
                                              height_uv);

                buffer_y          = input_picture_ptr->buffer_y;
                buffer_bit_inc_y  = uncompressed_pics[0];
                buffer_cb         = input_picture_ptr->buffer_cb;
                buffer_bit_inc_cb = uncompressed_pics[1];
                buffer_cr         = input_picture_ptr->buffer_cr;
                buffer_bit_inc_cr = uncompressed_pics[2];
            }

            bd    = 10;
            shift = 0; // both input and output are 10 bit (bitdepth - input_bd)

            input_buffer                = &((buffer_y)[input_picture_ptr->origin_x +
                                        input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
            EbByte input_buffer_bit_inc = &(
                (buffer_bit_inc_y)[input_picture_ptr->origin_x +
                                   input_picture_ptr->origin_y *
                                       input_picture_ptr->stride_bit_inc_y]);
            luma_ssim = aom_highbd_ssim2(input_buffer,
                                         input_picture_ptr->stride_y,
                                         input_buffer_bit_inc,
                                         input_picture_ptr->stride_bit_inc_y,
                                         recon_coeff_buffer,
                                         recon_ptr->stride_y,
                                         scs_ptr->seq_header.max_frame_width,
                                         scs_ptr->seq_header.max_frame_height,
                                         bd,
                                         shift);

            recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr
                     ->buffer_cb)[(recon_ptr->origin_x << is_16bit) / 2 +
                                  (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cb]));
            input_buffer       = &(
                (buffer_cb)[input_picture_ptr->origin_x / 2 +
                            input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
            input_buffer_bit_inc = &((buffer_bit_inc_cb)[input_picture_ptr->origin_x / 2 +
                                                         input_picture_ptr->origin_y / 2 *
                                                             input_picture_ptr->stride_bit_inc_cb]);
            cb_ssim              = aom_highbd_ssim2(input_buffer,
                                       input_picture_ptr->stride_cb,
                                       input_buffer_bit_inc,
                                       input_picture_ptr->stride_bit_inc_cb,
                                       recon_coeff_buffer,
                                       recon_ptr->stride_cb,
                                       scs_ptr->chroma_width,
                                       scs_ptr->chroma_height,
                                       bd,
                                       shift);

            recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr
                     ->buffer_cr)[(recon_ptr->origin_x << is_16bit) / 2 +
                                  (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cr]));
            input_buffer       = &(
                (buffer_cr)[input_picture_ptr->origin_x / 2 +
                            input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            input_buffer_bit_inc = &((buffer_bit_inc_cr)[input_picture_ptr->origin_x / 2 +
                                                         input_picture_ptr->origin_y / 2 *
                                                             input_picture_ptr->stride_bit_inc_cr]);
            cr_ssim              = aom_highbd_ssim2(input_buffer,
                                       input_picture_ptr->stride_cr,
                                       input_buffer_bit_inc,
                                       input_picture_ptr->stride_bit_inc_cr,
                                       recon_coeff_buffer,
                                       recon_ptr->stride_cr,
                                       scs_ptr->chroma_width,
                                       scs_ptr->chroma_height,
                                       bd,
                                       shift);

            pcs_ptr->parent_pcs_ptr->luma_ssim = luma_ssim;
            pcs_ptr->parent_pcs_ptr->cb_ssim   = cb_ssim;
            pcs_ptr->parent_pcs_ptr->cr_ssim   = cr_ssim;

            if (free_memory && pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
                EB_FREE_ARRAY(buffer_y);
                EB_FREE_ARRAY(buffer_cb);
                EB_FREE_ARRAY(buffer_cr);
                EB_FREE_ARRAY(buffer_bit_inc_y);
                EB_FREE_ARRAY(buffer_bit_inc_cb);
                EB_FREE_ARRAY(buffer_bit_inc_cr);
            }
            if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_FALSE) {
                EB_FREE_ARRAY(buffer_bit_inc_y);
                EB_FREE_ARRAY(buffer_bit_inc_cb);
                EB_FREE_ARRAY(buffer_bit_inc_cr);
            }
        }
    }
    return EB_ErrorNone;
}

EbErrorType psnr_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr,
                              EbBool free_memory) {
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
            recon_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr;

        EbPictureBufferDesc *input_picture_ptr =
            (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;

        uint64_t sse_total[3]        = {0};
        uint64_t residual_distortion = 0;
        EbByte   input_buffer;
        EbByte   recon_coeff_buffer;

        EbByte buffer_y;
        EbByte buffer_cb;
        EbByte buffer_cr;

        // if current source picture was temporally filtered, use an alternative buffer which stores
        // the original source picture
        if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            assert(pcs_ptr->parent_pcs_ptr->save_source_picture_width == input_picture_ptr->width &&
                   pcs_ptr->parent_pcs_ptr->save_source_picture_height ==
                       input_picture_ptr->height);
            buffer_y  = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[0];
            buffer_cb = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[1];
            buffer_cr = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[2];
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

        recon_coeff_buffer = &(
            (recon_ptr->buffer_cb)[recon_ptr->origin_x / 2 +
                                   recon_ptr->origin_y / 2 * recon_ptr->stride_cb]);
        input_buffer = &(buffer_cb[input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);

        residual_distortion = 0;
        for (int row_index = 0;
             row_index < (input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y;
             ++row_index) {
            for (int column_index = 0;
                 column_index < (input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x;
                 ++column_index) {
                residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                    (recon_coeff_buffer[column_index]));
            }

            input_buffer += input_picture_ptr->stride_cb;
            recon_coeff_buffer += recon_ptr->stride_cb;
        }

        sse_total[1] = residual_distortion;

        recon_coeff_buffer = &(
            (recon_ptr->buffer_cr)[recon_ptr->origin_x / 2 +
                                   recon_ptr->origin_y / 2 * recon_ptr->stride_cr]);
        input_buffer        = &(buffer_cr[input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
        residual_distortion = 0;

        for (int row_index = 0;
             row_index < (input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y;
             ++row_index) {
            for (int column_index = 0;
                 column_index < (input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x;
                 ++column_index) {
                residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                    (recon_coeff_buffer[column_index]));
            }

            input_buffer += input_picture_ptr->stride_cr;
            recon_coeff_buffer += recon_ptr->stride_cr;
        }

        sse_total[2]                      = residual_distortion;
        pcs_ptr->parent_pcs_ptr->luma_sse = sse_total[0];
        pcs_ptr->parent_pcs_ptr->cb_sse   = sse_total[1];
        pcs_ptr->parent_pcs_ptr->cr_sse   = sse_total[2];

        if (free_memory && pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            EB_FREE_ARRAY(buffer_y);
            EB_FREE_ARRAY(buffer_cb);
            EB_FREE_ARRAY(buffer_cr);
        }
    } else {
        EbPictureBufferDesc *recon_ptr;

        get_recon_pic(pcs_ptr, &recon_ptr, is_16bit);
        EbPictureBufferDesc *input_picture_ptr =
            (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;

        uint64_t  sse_total[3]        = {0};
        uint64_t  residual_distortion = 0;
        EbByte    input_buffer;
        EbByte    input_buffer_bit_inc;
        uint16_t *recon_coeff_buffer;

        if (scs_ptr->ten_bit_format == 1) {
            const uint32_t luma_width   = input_picture_ptr->width - scs_ptr->max_input_pad_right;
            const uint32_t luma_height  = input_picture_ptr->height - scs_ptr->max_input_pad_bottom;
            const uint32_t chroma_width = luma_width >> ss_x;
            const uint32_t pic_width_in_sb   = (luma_width + 64 - 1) / 64;
            const uint32_t pic_height_in_sb  = (luma_height + 64 - 1) / 64;
            const uint32_t luma_2bit_width   = luma_width / 4;
            const uint32_t chroma_height     = luma_height >> ss_y;
            const uint32_t chroma_2bit_width = chroma_width / 4;
            uint32_t       sb_num_in_height, sb_num_in_width;

            EbByte input_buffer_org = &(
                (input_picture_ptr
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
            uint16_t *recon_buffer_org_u = recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr
                     ->buffer_cb)[(recon_ptr->origin_x << is_16bit) / 2 +
                                  (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cb]));
            ;

            EbByte input_buffer_org_v = &(
                (input_picture_ptr
                     ->buffer_cr)[input_picture_ptr->origin_x / 2 +
                                  input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            ;
            uint16_t *recon_buffer_org_v = recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr
                     ->buffer_cr)[(recon_ptr->origin_x << is_16bit) / 2 +
                                  (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cr]));
            ;

            residual_distortion            = 0;
            uint64_t residual_distortion_u = 0;
            uint64_t residual_distortion_v = 0;

            for (sb_num_in_height = 0; sb_num_in_height < pic_height_in_sb; ++sb_num_in_height) {
                for (sb_num_in_width = 0; sb_num_in_width < pic_width_in_sb; ++sb_num_in_width) {
                    uint32_t tb_origin_x = sb_num_in_width * 64;
                    uint32_t tb_origin_y = sb_num_in_height * 64;
                    uint32_t sb_width = (luma_width - tb_origin_x) < 64 ? (luma_width - tb_origin_x)
                                                                        : 64;
                    uint32_t sb_height = (luma_height - tb_origin_y) < 64
                        ? (luma_height - tb_origin_y)
                        : 64;

                    input_buffer = input_buffer_org + tb_origin_y * input_picture_ptr->stride_y +
                        tb_origin_x;
                    input_buffer_bit_inc = input_picture_ptr->buffer_bit_inc_y +
                        tb_origin_y * luma_2bit_width + (tb_origin_x / 4) * sb_height;
                    recon_coeff_buffer = recon_buffer_org + tb_origin_y * recon_ptr->stride_y +
                        tb_origin_x;

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
                    sb_width    = (chroma_width - tb_origin_x) < 32 ? (chroma_width - tb_origin_x)
                                                                    : 32;
                    sb_height   = (chroma_height - tb_origin_y) < 32 ? (chroma_height - tb_origin_y)
                                                                     : 32;

                    inn_stride = sb_width / 4;

                    input_buffer = input_buffer_org_u + tb_origin_y * input_picture_ptr->stride_cb +
                        tb_origin_x;

                    input_buffer_bit_inc = input_picture_ptr->buffer_bit_inc_cb +
                        tb_origin_y * chroma_2bit_width + (tb_origin_x / 4) * sb_height;

                    recon_coeff_buffer = recon_buffer_org_u + tb_origin_y * recon_ptr->stride_cb +
                        tb_origin_x;

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
                        tb_origin_y * chroma_2bit_width + (tb_origin_x / 4) * sb_height;
                    recon_coeff_buffer = recon_buffer_org_v + tb_origin_y * recon_ptr->stride_cr +
                        tb_origin_x;

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
                assert(pcs_ptr->parent_pcs_ptr->save_source_picture_width ==
                           input_picture_ptr->width &&
                       pcs_ptr->parent_pcs_ptr->save_source_picture_height ==
                           input_picture_ptr->height);
                buffer_y          = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[0];
                buffer_bit_inc_y  = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[0];
                buffer_cb         = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[1];
                buffer_bit_inc_cb = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[1];
                buffer_cr         = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[2];
                buffer_bit_inc_cr = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[2];
            } else {
                uint32_t height_y  = (uint32_t)(input_picture_ptr->height +
                                               input_picture_ptr->origin_y +
                                               input_picture_ptr->origin_bot_y);
                uint32_t height_uv = (uint32_t)((input_picture_ptr->height +
                                                 input_picture_ptr->origin_y +
                                                 input_picture_ptr->origin_bot_y) >>
                                                ss_y);

                uint8_t *uncompressed_pics[3];
                EB_MALLOC_ARRAY(uncompressed_pics[0],
                                pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr->luma_size);
                EB_MALLOC_ARRAY(
                    uncompressed_pics[1],
                    pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr->chroma_size);
                EB_MALLOC_ARRAY(
                    uncompressed_pics[2],
                    pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr->chroma_size);

                svt_c_unpack_compressed_10bit(input_picture_ptr->buffer_bit_inc_y,
                                              input_picture_ptr->stride_bit_inc_y / 4,
                                              uncompressed_pics[0],
                                              input_picture_ptr->stride_bit_inc_y,
                                              height_y);
                // U
                svt_c_unpack_compressed_10bit(input_picture_ptr->buffer_bit_inc_cb,
                                              input_picture_ptr->stride_bit_inc_cb / 4,
                                              uncompressed_pics[1],
                                              input_picture_ptr->stride_bit_inc_cb,
                                              height_uv);
                // V
                svt_c_unpack_compressed_10bit(input_picture_ptr->buffer_bit_inc_cr,
                                              input_picture_ptr->stride_bit_inc_cr / 4,
                                              uncompressed_pics[2],
                                              input_picture_ptr->stride_bit_inc_cr,
                                              height_uv);

                buffer_y          = input_picture_ptr->buffer_y;
                buffer_bit_inc_y  = uncompressed_pics[0];
                buffer_cb         = input_picture_ptr->buffer_cb;
                buffer_bit_inc_cb = uncompressed_pics[1];
                buffer_cr         = input_picture_ptr->buffer_cr;
                buffer_bit_inc_cr = uncompressed_pics[2];
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
                for (int column_index = 0;
                     column_index < input_picture_ptr->width - scs_ptr->max_input_pad_right;
                     ++column_index) {
                    residual_distortion += (int64_t)SQR(
                        (int64_t)((((input_buffer[column_index]) << 2) |
                                   ((input_buffer_bit_inc[column_index] >> 6) & 3))) -
                        (recon_coeff_buffer[column_index]));
                }

                input_buffer += input_picture_ptr->stride_y;
                input_buffer_bit_inc += input_picture_ptr->stride_bit_inc_y;
                recon_coeff_buffer += recon_ptr->stride_y;
            }

            sse_total[0] = residual_distortion;

            recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr
                     ->buffer_cb)[(recon_ptr->origin_x << is_16bit) / 2 +
                                  (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cb]));
            input_buffer       = &(
                (buffer_cb)[input_picture_ptr->origin_x / 2 +
                            input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
            input_buffer_bit_inc = &((buffer_bit_inc_cb)[input_picture_ptr->origin_x / 2 +
                                                         input_picture_ptr->origin_y / 2 *
                                                             input_picture_ptr->stride_bit_inc_cb]);

            residual_distortion = 0;
            for (int row_index = 0;
                 row_index < (input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y;
                 ++row_index) {
                for (int column_index = 0; column_index <
                     (input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x;
                     ++column_index) {
                    residual_distortion += (int64_t)SQR(
                        (int64_t)((((input_buffer[column_index]) << 2) |
                                   ((input_buffer_bit_inc[column_index] >> 6) & 3))) -
                        (recon_coeff_buffer[column_index]));
                }

                input_buffer += input_picture_ptr->stride_cb;
                input_buffer_bit_inc += input_picture_ptr->stride_bit_inc_cb;
                recon_coeff_buffer += recon_ptr->stride_cb;
            }

            sse_total[1] = residual_distortion;

            recon_coeff_buffer = (uint16_t *)(&(
                (recon_ptr
                     ->buffer_cr)[(recon_ptr->origin_x << is_16bit) / 2 +
                                  (recon_ptr->origin_y << is_16bit) / 2 * recon_ptr->stride_cr]));
            input_buffer       = &(
                (buffer_cr)[input_picture_ptr->origin_x / 2 +
                            input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            input_buffer_bit_inc = &((buffer_bit_inc_cr)[input_picture_ptr->origin_x / 2 +
                                                         input_picture_ptr->origin_y / 2 *
                                                             input_picture_ptr->stride_bit_inc_cr]);

            residual_distortion = 0;

            for (int row_index = 0;
                 row_index < (input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y;
                 ++row_index) {
                for (int column_index = 0; column_index <
                     (input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x;
                     ++column_index) {
                    residual_distortion += (int64_t)SQR(
                        (int64_t)((((input_buffer[column_index]) << 2) |
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
            if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_FALSE) {
                EB_FREE_ARRAY(buffer_bit_inc_y);
                EB_FREE_ARRAY(buffer_bit_inc_cb);
                EB_FREE_ARRAY(buffer_bit_inc_cr);
            }
        }
        pcs_ptr->parent_pcs_ptr->luma_sse = sse_total[0];
        pcs_ptr->parent_pcs_ptr->cb_sse   = sse_total[1];
        pcs_ptr->parent_pcs_ptr->cr_sse   = sse_total[2];
    }
    return EB_ErrorNone;
}

void pad_ref_and_set_flags(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    EbReferenceObject *reference_object =
        (EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;

    EbPictureBufferDesc
        *ref_pic_ptr; //= (EbPictureBufferDesc *)reference_object->reference_picture;
    EbPictureBufferDesc
        *ref_pic_16bit_ptr; // =   (EbPictureBufferDesc *)reference_object->reference_picture16bit;

    {
        get_recon_pic(pcs_ptr, &ref_pic_ptr, 0);
        get_recon_pic(pcs_ptr, &ref_pic_16bit_ptr, 1);
    }
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
        generate_padding16_bit((uint16_t *)ref_pic_16bit_ptr->buffer_y,
                               ref_pic_16bit_ptr->stride_y,
                               ref_pic_16bit_ptr->width,
                               ref_pic_16bit_ptr->height,
                               ref_pic_16bit_ptr->origin_x,
                               ref_pic_16bit_ptr->origin_y);

        // Cb samples
        generate_padding16_bit((uint16_t *)ref_pic_16bit_ptr->buffer_cb,
                               ref_pic_16bit_ptr->stride_cb,
                               ref_pic_16bit_ptr->width >> 1,
                               ref_pic_16bit_ptr->height >> 1,
                               ref_pic_16bit_ptr->origin_x >> 1,
                               ref_pic_16bit_ptr->origin_y >> 1);

        // Cr samples
        generate_padding16_bit((uint16_t *)ref_pic_16bit_ptr->buffer_cr,
                               ref_pic_16bit_ptr->stride_cr,
                               ref_pic_16bit_ptr->width >> 1,
                               ref_pic_16bit_ptr->height >> 1,
                               ref_pic_16bit_ptr->origin_x >> 1,
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
    if ((scs_ptr->is_16bit_pipeline) && (!is_16bit)) {
        // Y samples
        generate_padding16_bit((uint16_t *)ref_pic_16bit_ptr->buffer_y,
                               ref_pic_16bit_ptr->stride_y,
                               ref_pic_16bit_ptr->width - scs_ptr->max_input_pad_right,
                               ref_pic_16bit_ptr->height - scs_ptr->max_input_pad_bottom,
                               ref_pic_16bit_ptr->origin_x,
                               ref_pic_16bit_ptr->origin_y);

        // Cb samples
        generate_padding16_bit((uint16_t *)ref_pic_16bit_ptr->buffer_cb,
                               ref_pic_16bit_ptr->stride_cb,
                               (ref_pic_16bit_ptr->width - scs_ptr->max_input_pad_right) >> 1,
                               (ref_pic_16bit_ptr->height - scs_ptr->max_input_pad_bottom) >> 1,
                               ref_pic_16bit_ptr->origin_x >> 1,
                               ref_pic_16bit_ptr->origin_y >> 1);

        // Cr samples
        generate_padding16_bit((uint16_t *)ref_pic_16bit_ptr->buffer_cr,
                               ref_pic_16bit_ptr->stride_cr,
                               (ref_pic_16bit_ptr->width - scs_ptr->max_input_pad_right) >> 1,
                               (ref_pic_16bit_ptr->height - scs_ptr->max_input_pad_bottom) >> 1,
                               ref_pic_16bit_ptr->origin_x >> 1,
                               ref_pic_16bit_ptr->origin_y >> 1);

        // Hsan: unpack ref samples (to be used @ MD)

        //Y
        uint16_t *buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_y);
        uint8_t  *buf_8bit  = ref_pic_ptr->buffer_y;
        svt_convert_16bit_to_8bit(buf_16bit,
                                  ref_pic_16bit_ptr->stride_y,
                                  buf_8bit,
                                  ref_pic_ptr->stride_y,
                                  ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1),
                                  ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1));

        //CB
        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cb);
        buf_8bit  = ref_pic_ptr->buffer_cb;
        svt_convert_16bit_to_8bit(buf_16bit,
                                  ref_pic_16bit_ptr->stride_cb,
                                  buf_8bit,
                                  ref_pic_ptr->stride_cb,
                                  (ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1)) >> 1,
                                  (ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1)) >> 1);

        //CR
        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cr);
        buf_8bit  = ref_pic_ptr->buffer_cr;
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
    reference_object->slice_type = pcs_ptr->parent_pcs_ptr->slice_type;
    reference_object->r0         = pcs_ptr->parent_pcs_ptr->r0;
}

void copy_statistics_to_ref_obj_ect(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    pcs_ptr->intra_coded_area = (100 * pcs_ptr->intra_coded_area) /
        (pcs_ptr->parent_pcs_ptr->aligned_width * pcs_ptr->parent_pcs_ptr->aligned_height);
    pcs_ptr->skip_coded_area = (100 * pcs_ptr->skip_coded_area) /
        (pcs_ptr->parent_pcs_ptr->aligned_width * pcs_ptr->parent_pcs_ptr->aligned_height);
    if (pcs_ptr->slice_type == I_SLICE)
        pcs_ptr->intra_coded_area = 0;
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->intra_coded_area = (uint8_t)(pcs_ptr->intra_coded_area);
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->skip_coded_area                   = (uint8_t)(pcs_ptr->skip_coded_area);
    struct PictureParentControlSet *ppcs    = pcs_ptr->parent_pcs_ptr;
    FrameHeader                    *frm_hdr = &ppcs->frm_hdr;
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->ref_cdef_strengths_num = ppcs->nb_cdef_strengths;
    for (int i = 0; i < ppcs->nb_cdef_strengths; i++) {
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->ref_cdef_strengths[0][i] = frm_hdr->cdef_params.cdef_y_strength[i];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->ref_cdef_strengths[1][i] = frm_hdr->cdef_params.cdef_uv_strength[i];
    }
    uint32_t sb_index;
    for (sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index) {
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->sb_intra[sb_index] = pcs_ptr->sb_intra[sb_index];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->sb_skip[sb_index] = pcs_ptr->sb_skip[sb_index];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->sb_64x64_mvp[sb_index] = pcs_ptr->sb_64x64_mvp[sb_index];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->sb_me_64x64_dist[sb_index] = pcs_ptr->parent_pcs_ptr->me_64x64_distortion[sb_index];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->sb_me_8x8_cost_var[sb_index] =
            pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[sb_index];
    }
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->tmp_layer_idx = (uint8_t)pcs_ptr->temporal_layer_index;
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->is_scene_change = pcs_ptr->parent_pcs_ptr->scene_change_flag;

    Av1Common *cm = pcs_ptr->parent_pcs_ptr->av1_cm;
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->sg_frame_ep = cm->sg_frame_ep;
    if (scs_ptr->mfmv_enabled) {
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->frame_type = pcs_ptr->parent_pcs_ptr->frm_hdr.frame_type;
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->order_hint = pcs_ptr->parent_pcs_ptr->cur_order_hint;
        svt_memcpy(((EbReferenceObject *)
                        pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                       ->ref_order_hint,
                   pcs_ptr->parent_pcs_ptr->ref_order_hint,
                   7 * sizeof(uint32_t));
    }
    // Copy the prev frame wn filter coeffs
    EbReferenceObject *obj = (EbReferenceObject *)
                                 pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
    if (cm->wn_filter_ctrls.enabled && cm->wn_filter_ctrls.use_prev_frame_coeffs) {
        for (int32_t plane = 0; plane < MAX_MB_PLANE; ++plane) {
            int32_t ntiles = pcs_ptr->rst_info[plane].units_per_tile;
            for (int32_t u = 0; u < ntiles; ++u) {
                obj->unit_info[plane][u].restoration_type =
                    pcs_ptr->rst_info[plane].unit_info[u].restoration_type;
                if (pcs_ptr->rst_info[plane].unit_info[u].restoration_type == RESTORE_WIENER)
                    obj->unit_info[plane][u].wiener_info =
                        pcs_ptr->rst_info[plane].unit_info[u].wiener_info;
            }
        }
    }
}

void set_obmc_controls(ModeDecisionContext *mdctxt, uint8_t obmc_mode) {
    ObmcControls *obmc_ctrls = &mdctxt->obmc_ctrls;
    switch (obmc_mode) {
    case 0: obmc_ctrls->enabled = 0; break;
    case 1:
        obmc_ctrls->enabled            = 1;
        obmc_ctrls->max_blk_size_16x16 = 0;
        break;
    case 2:
        obmc_ctrls->enabled            = 1;
        obmc_ctrls->max_blk_size_16x16 = 1;
        break;
    default: assert(0); break;
    }
}
/*
 * Generate depth removal settings
 */

#define LOW_8x8_DIST_VAR_TH 25000
#define HIGH_8x8_DIST_VAR_TH 50000

void set_depth_removal_level_controls(PictureControlSet *pcs_ptr, ModeDecisionContext *mdctxt,
                                      uint8_t depth_removal_level) {
    DepthRemovalCtrls *depth_removal_ctrls = &mdctxt->depth_removal_ctrls;

    if (pcs_ptr->slice_type == I_SLICE) {
        SbParams *sb_params = &pcs_ptr->parent_pcs_ptr->sb_params_array[mdctxt->sb_index];

        uint16_t disallow_below_16x16_variance_th = 0;
        uint16_t disallow_below_32x32_variance_th = 0;
        uint16_t disallow_below_64x64_variance_th = 0;

        switch (depth_removal_level) {
        case 0: depth_removal_ctrls->enabled = 0; break;

        case 1:
            depth_removal_ctrls->enabled     = 1;
            disallow_below_16x16_variance_th = 150;
            disallow_below_32x32_variance_th = 50;
            disallow_below_64x64_variance_th = 25;
            break;
        }

        if (depth_removal_ctrls->enabled) {
            // If variance is available, use in depth removal decision
            if (pcs_ptr->parent_pcs_ptr->variance) {
                depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 &&
                                                             sb_params->height % 16 == 0)
                    ? (depth_removal_ctrls->disallow_below_16x16 ||
                       pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <
                           disallow_below_16x16_variance_th)
                    : 0;

                depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 &&
                                                             sb_params->height % 32 == 0)
                    ? (depth_removal_ctrls->disallow_below_32x32 ||
                       pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <
                           disallow_below_32x32_variance_th)
                    : 0;

                depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 &&
                                                             sb_params->height % 64 == 0)
                    ? (depth_removal_ctrls->disallow_below_64x64 ||
                       pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <
                           disallow_below_64x64_variance_th)
                    : 0;
            } else {
                depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 &&
                                                             sb_params->height % 16 == 0)
                    ? depth_removal_ctrls->disallow_below_16x16
                    : 0;

                depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 &&
                                                             sb_params->height % 32 == 0)
                    ? depth_removal_ctrls->disallow_below_32x32
                    : 0;

                depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 &&
                                                             sb_params->height % 64 == 0)
                    ? depth_removal_ctrls->disallow_below_64x64
                    : 0;
            }
        }
    } else {
        uint32_t me_8x8_cost_variance =
            pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[mdctxt->sb_index];

        SbParams *sb_params = &pcs_ptr->parent_pcs_ptr->sb_params_array[mdctxt->sb_index];

        // me_distortion => EB_8_BIT_MD
        uint32_t fast_lambda = mdctxt->fast_lambda_md[EB_8_BIT_MD];

        uint32_t sb_size = 64 * 64;

        uint64_t cost_th_rate = 1 << 13;

        uint64_t disallow_below_16x16_cost_th_multiplier = 0;
        uint64_t disallow_below_32x32_cost_th_multiplier = 0;
        uint64_t disallow_below_64x64_cost_th_multiplier = 0;

        int64_t dev_16x16_to_8x8_th   = MAX_SIGNED_VALUE;
        int64_t dev_32x32_to_16x16_th = MAX_SIGNED_VALUE;
        int64_t dev_32x32_to_8x8_th   = MAX_SIGNED_VALUE;

        int8_t qp_scale_factor = 0;

        // Modulate depth_removal level for Layer0 frames based on the qp_offset band
        if (pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present) {
            int diff = mdctxt->sb_ptr->qindex -
                quantizer_to_qindex[pcs_ptr->parent_pcs_ptr->picture_qp];
            if (diff <= -12)
                depth_removal_level = MAX(0, (int)depth_removal_level - 4);
            else if (diff <= -6)
                depth_removal_level = MAX(0, (int)depth_removal_level - 3);
            else if (diff <= -3)
                depth_removal_level = MAX(0, (int)depth_removal_level - 2);
            else if (diff < 0)
                depth_removal_level = MAX(0, (int)depth_removal_level - 1);
        }

        switch (depth_removal_level) {
        case 0: depth_removal_ctrls->enabled = 0; break;

        case 1:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 0;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 15;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 1;

            break;

        case 2:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 20;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 1;

            break;
        case 3:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 1;
            break;
        case 4:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 2;
            break;
        case 5:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 2;
            break;
        case 6:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 2;
            break;
        case 7:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 2;
            break;
        case 8:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 100;
            dev_32x32_to_16x16_th                   = 50;
            qp_scale_factor                         = 3;
            break;
        case 9:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 8;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 100;
            dev_32x32_to_16x16_th                   = 50;
            qp_scale_factor                         = 3;

            break;
        case 10:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 32;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 200;
            dev_32x32_to_16x16_th                   = 75;
            qp_scale_factor                         = 3;
            break;
        case 11:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 32;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 250;
            dev_32x32_to_16x16_th                   = 125;
            qp_scale_factor                         = 3;
            break;
        case 12:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 32;
            disallow_below_32x32_cost_th_multiplier = 4;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 250;
            dev_32x32_to_16x16_th                   = 150;
            qp_scale_factor                         = 4;
            break;
        case 13:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 64;
            disallow_below_32x32_cost_th_multiplier = 4;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 250;
            dev_32x32_to_16x16_th                   = 150;
            qp_scale_factor                         = 4;
            break;
        case 14:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 64;
            disallow_below_32x32_cost_th_multiplier = 4;
            disallow_below_64x64_cost_th_multiplier = 4;
            dev_16x16_to_8x8_th                     = 250;
            dev_32x32_to_16x16_th                   = 150;
            qp_scale_factor                         = 4;
            break;
        case 15:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 96;
            disallow_below_32x32_cost_th_multiplier = 6;
            disallow_below_64x64_cost_th_multiplier = 6;
            dev_16x16_to_8x8_th                     = 300;
            dev_32x32_to_16x16_th                   = 200;
            qp_scale_factor                         = 4;
            break;
        }
        if (depth_removal_ctrls->enabled) {
            //dev_16x16_to_8x8_th , dev_32x32_to_16x16_th = f(me_8x8_cost_variance)
            me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);
            if (me_8x8_cost_variance < LOW_8x8_DIST_VAR_TH) {
                dev_16x16_to_8x8_th = dev_16x16_to_8x8_th << 2;
            } else if (me_8x8_cost_variance < HIGH_8x8_DIST_VAR_TH) {
                dev_16x16_to_8x8_th   = dev_16x16_to_8x8_th << 1;
                dev_32x32_to_16x16_th = dev_32x32_to_16x16_th >> 1;
            } else {
                dev_16x16_to_8x8_th   = 0;
                dev_32x32_to_16x16_th = 0;
            }

            //dev_16x16_to_8x8_th , dev_32x32_to_16x16_th = f(QP)
            dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) *
                qp_scale_factor;
            dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) *
                qp_scale_factor;

            // dev_32x32_to_8x8_th = f(dev_32x32_to_16x16_th); a bit higher
            dev_32x32_to_8x8_th = (dev_32x32_to_16x16_th * ((1 << 2) + 1)) >> 2;

            uint64_t disallow_below_16x16_cost_th = disallow_below_16x16_cost_th_multiplier
                ? RDCOST(fast_lambda,
                         cost_th_rate,
                         (sb_size >> 1) * disallow_below_16x16_cost_th_multiplier)
                : 0;
            uint64_t disallow_below_32x32_cost_th = disallow_below_32x32_cost_th_multiplier
                ? RDCOST(fast_lambda,
                         cost_th_rate,
                         (sb_size >> 1) * disallow_below_32x32_cost_th_multiplier)
                : 0;
            uint64_t disallow_below_64x64_cost_th = disallow_below_64x64_cost_th_multiplier
                ? RDCOST(fast_lambda,
                         cost_th_rate,
                         (sb_size >> 1) * disallow_below_64x64_cost_th_multiplier)
                : 0;

            uint64_t cost_64x64 = RDCOST(
                fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[mdctxt->sb_index]);
            uint64_t cost_32x32 = RDCOST(
                fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_32x32_distortion[mdctxt->sb_index]);
            uint64_t cost_16x16 = RDCOST(
                fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_16x16_distortion[mdctxt->sb_index]);
            uint64_t cost_8x8 = RDCOST(
                fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_8x8_distortion[mdctxt->sb_index]);

            int64_t dev_32x32_to_16x16 =
                (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_16x16, 1)) * 1000) /
                (int64_t)MAX(cost_16x16, 1);

            int64_t dev_32x32_to_8x8 =
                (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_8x8, 1)) * 1000) /
                (int64_t)MAX(cost_8x8, 1);

            int64_t dev_16x16_to_8x8 =
                (int64_t)(((int64_t)MAX(cost_16x16, 1) - (int64_t)MAX(cost_8x8, 1)) * 1000) /
                (int64_t)MAX(cost_8x8, 1);
            depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 &&
                                                         sb_params->height % 64 == 0)
                ? (depth_removal_ctrls->disallow_below_64x64 ||
                   cost_64x64 < disallow_below_64x64_cost_th)
                : 0;

            depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 &&
                                                         sb_params->height % 32 == 0)
                ? (depth_removal_ctrls->disallow_below_32x32 ||
                   cost_32x32 < disallow_below_32x32_cost_th ||
                   (dev_32x32_to_16x16 < dev_32x32_to_16x16_th &&
                    dev_32x32_to_8x8 < dev_32x32_to_8x8_th))
                : 0;

            depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 &&
                                                         sb_params->height % 16 == 0)
                ? (depth_removal_ctrls->disallow_below_16x16 ||
                   cost_16x16 < disallow_below_16x16_cost_th ||
                   dev_16x16_to_8x8 < dev_16x16_to_8x8_th)
                : 0;
        }
    }
}
/*
 * Control NSQ search
 */
void md_nsq_motion_search_controls(ModeDecisionContext *mdctxt, uint8_t md_nsq_mv_search_level) {
    MdNsqMotionSearchCtrls *md_nsq_motion_search_ctrls = &mdctxt->md_nsq_motion_search_ctrls;

    switch (md_nsq_mv_search_level) {
    case 0: md_nsq_motion_search_ctrls->enabled = 0; break;
    case 1:
        md_nsq_motion_search_ctrls->enabled                = 1;
        md_nsq_motion_search_ctrls->use_ssd                = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width  = 31;
        md_nsq_motion_search_ctrls->full_pel_search_height = 31;
        md_nsq_motion_search_ctrls->enable_psad            = 1;
        break;

    case 2:
        md_nsq_motion_search_ctrls->enabled                = 1;
        md_nsq_motion_search_ctrls->use_ssd                = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width  = 15;
        md_nsq_motion_search_ctrls->full_pel_search_height = 15;
        md_nsq_motion_search_ctrls->enable_psad            = 1;
        break;
    case 3:
        md_nsq_motion_search_ctrls->enabled                = 1;
        md_nsq_motion_search_ctrls->use_ssd                = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width  = 11;
        md_nsq_motion_search_ctrls->full_pel_search_height = 11;
        md_nsq_motion_search_ctrls->enable_psad            = 1;
        break;
    case 4:
        md_nsq_motion_search_ctrls->enabled                = 1;
        md_nsq_motion_search_ctrls->use_ssd                = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width  = 8;
        md_nsq_motion_search_ctrls->full_pel_search_height = 7;
        md_nsq_motion_search_ctrls->enable_psad            = 1;
        break;
    default: assert(0); break;
    }
}
void md_pme_search_controls(ModeDecisionContext *mdctxt, uint8_t md_pme_level) {
    MdPmeCtrls *md_pme_ctrls = &mdctxt->md_pme_ctrls;

    switch (md_pme_level) {
    case 0: md_pme_ctrls->enabled = 0; break;
    case 1:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 1;
        md_pme_ctrls->full_pel_search_width         = 31;
        md_pme_ctrls->full_pel_search_height        = 31;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = MIN_SIGNED_VALUE;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 2:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 1;
        md_pme_ctrls->full_pel_search_width         = 20;
        md_pme_ctrls->full_pel_search_height        = 15;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = MIN_SIGNED_VALUE;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 3:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 1;
        md_pme_ctrls->full_pel_search_width         = 7;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = MIN_SIGNED_VALUE;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 4:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 1;
        md_pme_ctrls->full_pel_search_width         = 7;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 100;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 25;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 5:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 8;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 100;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 25;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 1;
        break;
    case 6:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 8;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 1;
        break;
    case 7:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 8;
        md_pme_ctrls->full_pel_search_height        = 3;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 1;
        break;
    case 8:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 3;
        md_pme_ctrls->full_pel_search_height        = 3;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 9:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 3;
        md_pme_ctrls->full_pel_search_height        = 3;
        md_pme_ctrls->early_check_mv_th_multiplier  = 50;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 1;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 10:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 2;
        md_pme_ctrls->full_pel_search_height        = 2;
        md_pme_ctrls->early_check_mv_th_multiplier  = 80;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 30;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 20;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 10;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 40;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 1;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    default: assert(0); break;
    }
}
void set_subres_controls(ModeDecisionContext *mdctxt, uint8_t subres_level) {
    SubresCtrls *subres_ctrls = &mdctxt->subres_ctrls;

    switch (subres_level) {
    case 0: subres_ctrls->step = 0; break;
    case 1: subres_ctrls->step = 1; break;
    case 2: subres_ctrls->step = 2; break;
    default: assert(0); break;
    }
    // Set the TH used to determine if subres is safe to use (based on ODD vs. EVEN rows' distortion)
    if (subres_ctrls->step == 0)
        subres_ctrls->odd_to_even_deviation_th = 0;
    else
        subres_ctrls->odd_to_even_deviation_th = 5;
}
void set_pf_controls(ModeDecisionContext *mdctxt, uint8_t pf_level) {
    PfCtrls *pf_ctrls = &mdctxt->pf_ctrls;

    switch (pf_level) {
    case 0: pf_ctrls->pf_shape = ONLY_DC_SHAPE; break;
    case 1: pf_ctrls->pf_shape = DEFAULT_SHAPE; break;
    case 2: pf_ctrls->pf_shape = N2_SHAPE; break;
    case 3: pf_ctrls->pf_shape = N4_SHAPE; break;
    default: assert(0); break;
    }
}
void set_block_based_depth_refinement_controls(ModeDecisionContext *mdctxt,
                                               uint8_t block_based_depth_refinement_level) {
    DepthRefinementCtrls *depth_refinement_ctrls = &mdctxt->depth_refinement_ctrls;

    switch (block_based_depth_refinement_level) {
    case 0: depth_refinement_ctrls->enabled = 0; break;

    case 1:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 25;
        depth_refinement_ctrls->sub_to_current_th          = 25;
        depth_refinement_ctrls->cost_band_based_modulation = 0;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;

    case 2:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 25;
        depth_refinement_ctrls->sub_to_current_th          = 25;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = 15;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;

    case 3:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 25;
        depth_refinement_ctrls->sub_to_current_th          = 25;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;

    case 4:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 20;
        depth_refinement_ctrls->sub_to_current_th          = 20;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;

    case 5:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 15;
        depth_refinement_ctrls->sub_to_current_th          = 15;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;

    case 6:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 10;
        depth_refinement_ctrls->sub_to_current_th          = 10;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;

    case 7:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 5;
        depth_refinement_ctrls->sub_to_current_th          = 5;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;

    case 8:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 5;
        depth_refinement_ctrls->sub_to_current_th          = 5;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 800;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;
    case 9:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = 5;
        depth_refinement_ctrls->sub_to_current_th          = -50;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 800;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 1;
        break;
    case 10:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = -25;
        depth_refinement_ctrls->sub_to_current_th          = -50;
        depth_refinement_ctrls->cost_band_based_modulation = 1;
        depth_refinement_ctrls->max_cost_multiplier        = 800;
        depth_refinement_ctrls->max_band_cnt               = 4;
        depth_refinement_ctrls->decrement_per_band[0]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      = MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      = 10;
        depth_refinement_ctrls->decrement_per_band[3]      = 5;
        depth_refinement_ctrls->up_to_2_depth              = 1;
        break;

    // Pred_Only
    case 11:
        depth_refinement_ctrls->enabled                    = 1;
        depth_refinement_ctrls->parent_to_current_th       = MIN_SIGNED_VALUE;
        depth_refinement_ctrls->sub_to_current_th          = MIN_SIGNED_VALUE;
        depth_refinement_ctrls->cost_band_based_modulation = 0;
        depth_refinement_ctrls->up_to_2_depth              = 0;
        break;
    }
}
/*
 * Control Adaptive ME search
 */
void md_sq_motion_search_controls(ModeDecisionContext *mdctxt, uint8_t md_sq_mv_search_level) {
    MdSqMotionSearchCtrls *md_sq_me_ctrls = &mdctxt->md_sq_me_ctrls;

    switch (md_sq_mv_search_level) {
    case 0: md_sq_me_ctrls->enabled = 0; break;
    case 1:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled    = 1;
        md_sq_me_ctrls->sprs_lev0_step       = 4;
        md_sq_me_ctrls->sprs_lev0_w          = 15;
        md_sq_me_ctrls->sprs_lev0_h          = 15;
        md_sq_me_ctrls->max_sprs_lev0_w      = 150;
        md_sq_me_ctrls->max_sprs_lev0_h      = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 500;

        md_sq_me_ctrls->sprs_lev1_enabled    = 1;
        md_sq_me_ctrls->sprs_lev1_step       = 2;
        md_sq_me_ctrls->sprs_lev1_w          = 4;
        md_sq_me_ctrls->sprs_lev1_h          = 4;
        md_sq_me_ctrls->max_sprs_lev1_w      = 50;
        md_sq_me_ctrls->max_sprs_lev1_h      = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 500;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step    = 1;
        md_sq_me_ctrls->sprs_lev2_w       = 3;
        md_sq_me_ctrls->sprs_lev2_h       = 3;
        md_sq_me_ctrls->enable_psad       = 1;
        break;
    case 2:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled    = 1;
        md_sq_me_ctrls->sprs_lev0_step       = 4;
        md_sq_me_ctrls->sprs_lev0_w          = 15;
        md_sq_me_ctrls->sprs_lev0_h          = 15;
        md_sq_me_ctrls->max_sprs_lev0_w      = 150;
        md_sq_me_ctrls->max_sprs_lev0_h      = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 400;

        md_sq_me_ctrls->sprs_lev1_enabled    = 1;
        md_sq_me_ctrls->sprs_lev1_step       = 2;
        md_sq_me_ctrls->sprs_lev1_w          = 4;
        md_sq_me_ctrls->sprs_lev1_h          = 4;
        md_sq_me_ctrls->max_sprs_lev1_w      = 50;
        md_sq_me_ctrls->max_sprs_lev1_h      = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 400;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step    = 1;
        md_sq_me_ctrls->sprs_lev2_w       = 3;
        md_sq_me_ctrls->sprs_lev2_h       = 3;
        md_sq_me_ctrls->enable_psad       = 1;
        break;
    case 3:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled    = 1;
        md_sq_me_ctrls->sprs_lev0_step       = 4;
        md_sq_me_ctrls->sprs_lev0_w          = 15;
        md_sq_me_ctrls->sprs_lev0_h          = 15;
        md_sq_me_ctrls->max_sprs_lev0_w      = 150;
        md_sq_me_ctrls->max_sprs_lev0_h      = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 300;

        md_sq_me_ctrls->sprs_lev1_enabled    = 1;
        md_sq_me_ctrls->sprs_lev1_step       = 2;
        md_sq_me_ctrls->sprs_lev1_w          = 4;
        md_sq_me_ctrls->sprs_lev1_h          = 4;
        md_sq_me_ctrls->max_sprs_lev1_w      = 50;
        md_sq_me_ctrls->max_sprs_lev1_h      = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 300;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step    = 1;
        md_sq_me_ctrls->sprs_lev2_w       = 3;
        md_sq_me_ctrls->sprs_lev2_h       = 3;
        md_sq_me_ctrls->enable_psad       = 1;
        break;
    case 4:
        md_sq_me_ctrls->enabled            = 1;
        md_sq_me_ctrls->use_ssd            = 0;
        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled    = 1;
        md_sq_me_ctrls->sprs_lev0_step       = 4;
        md_sq_me_ctrls->sprs_lev0_w          = 15;
        md_sq_me_ctrls->sprs_lev0_h          = 15;
        md_sq_me_ctrls->max_sprs_lev0_w      = 150;
        md_sq_me_ctrls->max_sprs_lev0_h      = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 100;

        md_sq_me_ctrls->sprs_lev1_enabled    = 1;
        md_sq_me_ctrls->sprs_lev1_step       = 2;
        md_sq_me_ctrls->sprs_lev1_w          = 4;
        md_sq_me_ctrls->sprs_lev1_h          = 4;
        md_sq_me_ctrls->max_sprs_lev1_w      = 50;
        md_sq_me_ctrls->max_sprs_lev1_h      = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 100;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step    = 1;
        md_sq_me_ctrls->sprs_lev2_w       = 3;
        md_sq_me_ctrls->sprs_lev2_h       = 3;
        md_sq_me_ctrls->enable_psad       = 1;
        break;
    default: assert(0); break;
    }
}
/*
 * Control Subpel search of ME MV(s)
 */
void md_subpel_me_controls(ModeDecisionContext *mdctxt, uint8_t md_subpel_me_level) {
    MdSubPelSearchCtrls *md_subpel_me_ctrls = &mdctxt->md_subpel_me_ctrls;

    switch (md_subpel_me_level) {
    case 0: md_subpel_me_ctrls->enabled = 0; break;
    case 1:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 0;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 2:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 0;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 3:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 0;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 4:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 1;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 5:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 2;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 6:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 2;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 7:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 1;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 8:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 2;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 9:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 3;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 10:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 2;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 11:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 3;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 12:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 50;
        md_subpel_me_ctrls->abs_th_mult           = 2;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 3;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 13:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 100;
        md_subpel_me_ctrls->abs_th_mult           = 2;
        md_subpel_me_ctrls->round_dev_th          = -25;
        md_subpel_me_ctrls->skip_diag_refinement  = 4;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        break;
    case 14:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 100;
        md_subpel_me_ctrls->abs_th_mult           = 2;
        md_subpel_me_ctrls->round_dev_th          = -25;
        md_subpel_me_ctrls->skip_diag_refinement  = 4;
        md_subpel_me_ctrls->skip_zz_mv            = 1;
        break;
    case 15:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = HALF_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 100;
        md_subpel_me_ctrls->abs_th_mult           = 2;
        md_subpel_me_ctrls->round_dev_th          = -25;
        md_subpel_me_ctrls->skip_diag_refinement  = 4;
        md_subpel_me_ctrls->skip_zz_mv            = 1;
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
    case 0: md_subpel_pme_ctrls->enabled = 0; break;
    case 1:
        md_subpel_pme_ctrls->enabled               = 1;
        md_subpel_pme_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_pme_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_pme_ctrls->pred_variance_th      = 0;
        md_subpel_pme_ctrls->abs_th_mult           = 0;
        md_subpel_pme_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv            = 0;
        break;
    case 2:
        md_subpel_pme_ctrls->enabled               = 1;
        md_subpel_pme_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_pme_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_pme_ctrls->pred_variance_th      = 0;
        md_subpel_pme_ctrls->abs_th_mult           = 0;
        md_subpel_pme_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv            = 0;
        break;
    case 3:
        md_subpel_pme_ctrls->enabled               = 1;
        md_subpel_pme_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_pme_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_pme_ctrls->pred_variance_th      = 0;
        md_subpel_pme_ctrls->abs_th_mult           = 0;
        md_subpel_pme_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv            = 0;
        break;
    case 4:
        md_subpel_pme_ctrls->enabled               = 1;
        md_subpel_pme_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->max_precision         = HALF_PEL;
        md_subpel_pme_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_pme_ctrls->pred_variance_th      = 0;
        md_subpel_pme_ctrls->abs_th_mult           = 0;
        md_subpel_pme_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv            = 0;
        break;
    default: assert(0); break;
    }
}
/*
 * Control RDOQ
 */
void set_rdoq_controls(ModeDecisionContext *mdctxt, uint8_t rdoq_level) {
    RdoqCtrls *rdoq_ctrls = &mdctxt->rdoq_ctrls;

    switch (rdoq_level) {
    case 0: rdoq_ctrls->enabled = 0; break;
    case 1:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 0;
        rdoq_ctrls->dct_dct_only      = 0;
        rdoq_ctrls->eob_th            = (uint8_t)~0;
        rdoq_ctrls->eob_fast_th       = (uint8_t)~0;
        break;
    case 2:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 1;
        rdoq_ctrls->dct_dct_only      = 0;
        rdoq_ctrls->eob_th            = (uint8_t)~0;
        rdoq_ctrls->eob_fast_th       = (uint8_t)~0;
        break;
    case 3:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 1;
        rdoq_ctrls->dct_dct_only      = 1;
        rdoq_ctrls->eob_th            = (uint8_t)~0;
        rdoq_ctrls->eob_fast_th       = (uint8_t)~0;
        break;
    case 4:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 1;
        rdoq_ctrls->dct_dct_only      = 1;
        rdoq_ctrls->eob_th            = (uint8_t)~0;
        rdoq_ctrls->eob_fast_th       = 30;
        break;
    case 5:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 1;
        rdoq_ctrls->dct_dct_only      = 1;
        rdoq_ctrls->eob_th            = 85;
        rdoq_ctrls->eob_fast_th       = 0;
        break;
    default: assert(0); break;
    }
}
/*
 * Settings for the parent SQ coeff-area based cycles reduction algorithm.
 */
void set_parent_sq_coeff_area_based_cycles_reduction_ctrls(ModeDecisionContext *ctx,
                                                           uint8_t              resolution,
                                                           uint8_t              cycles_alloc_lvl) {
    ParentSqCoeffAreaBasedCyclesReductionCtrls *cycle_red_ctrls =
        &ctx->parent_sq_coeff_area_based_cycles_reduction_ctrls;
    switch (cycles_alloc_lvl) {
    case 0: cycle_red_ctrls->enabled = 0; break;
    case 1:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = resolution <= INPUT_SIZE_360p_RANGE
               ? UNUSED_HIGH_FREQ_BAND_TH
               : 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th    = resolution <= INPUT_SIZE_360p_RANGE
               ? UNUSED_HIGH_FREQ_BAND_TH
               : 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th    = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 1;
        cycle_red_ctrls->enable_one_coeff_action  = 0;
        cycle_red_ctrls->one_coeff_action         = 0;

        cycle_red_ctrls->low_freq_band1_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 2:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = resolution <= INPUT_SIZE_360p_RANGE
               ? UNUSED_HIGH_FREQ_BAND_TH
               : 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th    = resolution <= INPUT_SIZE_360p_RANGE
               ? UNUSED_HIGH_FREQ_BAND_TH
               : 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th    = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 2;
        cycle_red_ctrls->enable_one_coeff_action  = 0;
        cycle_red_ctrls->one_coeff_action         = 0;

        cycle_red_ctrls->low_freq_band1_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 3:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = 90;
        cycle_red_ctrls->high_freq_band1_level = resolution <= INPUT_SIZE_360p_RANGE ? 2 : 3;
        cycle_red_ctrls->high_freq_band2_th    = 70;
        cycle_red_ctrls->high_freq_band2_level = resolution <= INPUT_SIZE_360p_RANGE ? 1 : 3;
        cycle_red_ctrls->high_freq_band3_th    = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 2;
        cycle_red_ctrls->enable_one_coeff_action  = 0;
        cycle_red_ctrls->one_coeff_action         = 0;

        cycle_red_ctrls->low_freq_band1_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 4:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = 90;
        cycle_red_ctrls->high_freq_band1_level = resolution <= INPUT_SIZE_360p_RANGE ? 2 : 3;
        cycle_red_ctrls->high_freq_band2_th    = 70;
        cycle_red_ctrls->high_freq_band2_level = resolution <= INPUT_SIZE_360p_RANGE ? 1 : 3;
        cycle_red_ctrls->high_freq_band3_th    = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 3;
        cycle_red_ctrls->enable_one_coeff_action  = 1;
        cycle_red_ctrls->one_coeff_action         = 1;

        cycle_red_ctrls->low_freq_band1_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 5:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th    = 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th    = 50;
        cycle_red_ctrls->high_freq_band3_level = 1;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 2;
        cycle_red_ctrls->enable_one_coeff_action  = 0;
        cycle_red_ctrls->one_coeff_action         = 0;

        cycle_red_ctrls->low_freq_band1_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 6:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = 90;
        cycle_red_ctrls->high_freq_band1_level = resolution <= INPUT_SIZE_360p_RANGE ? 2 : 3;
        cycle_red_ctrls->high_freq_band2_th    = 70;
        cycle_red_ctrls->high_freq_band2_level = resolution <= INPUT_SIZE_360p_RANGE ? 1 : 3;
        cycle_red_ctrls->high_freq_band3_th    = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 0;
        cycle_red_ctrls->enable_one_coeff_action  = 1;
        cycle_red_ctrls->one_coeff_action         = 1;

        cycle_red_ctrls->low_freq_band1_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;

    case 7:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = 90;
        cycle_red_ctrls->high_freq_band1_level = 0;
        cycle_red_ctrls->high_freq_band2_th    = 70;
        cycle_red_ctrls->high_freq_band2_level = 3;
        cycle_red_ctrls->high_freq_band3_th    = 50;
        cycle_red_ctrls->high_freq_band3_level = 2;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 0;
        cycle_red_ctrls->enable_one_coeff_action  = 1;
        cycle_red_ctrls->one_coeff_action         = 1;

        cycle_red_ctrls->low_freq_band1_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th    = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    default: assert(0); break;
    }
}
void set_txt_controls(ModeDecisionContext *mdctxt, uint8_t txt_level) {
    TxtControls *txt_ctrls = &mdctxt->txt_ctrls;

    switch (txt_level) {
    case 0:
        txt_ctrls->enabled = 0;

        txt_ctrls->txt_group_inter_lt_16x16    = 1;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;

        txt_ctrls->txt_group_intra_lt_16x16    = 1;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 1;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        break;
    case 1:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = MAX_TX_TYPE_GROUP;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        break;
    case 2:
        txt_ctrls->enabled                     = 1;
        txt_ctrls->txt_group_inter_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 5;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        break;
    case 3:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 5;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 5;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        break;
    case 4:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 5;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 3;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        break;
    case 5:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 4;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        break;
    case 6:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 4;

        txt_ctrls->early_exit_dist_th  = 100;
        txt_ctrls->early_exit_coeff_th = 4;
        break;
    case 7:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 3;

        txt_ctrls->early_exit_dist_th  = 400;
        txt_ctrls->early_exit_coeff_th = 8;
        break;
    case 8:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16    = 4;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
        txt_ctrls->early_exit_dist_th          = 400;
        txt_ctrls->early_exit_coeff_th         = 8;
        break;
    case 9:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;

        txt_ctrls->txt_group_intra_lt_16x16    = 4;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
        txt_ctrls->early_exit_dist_th          = 400;
        txt_ctrls->early_exit_coeff_th         = 8;
        break;
    case 10:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 2;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;
        txt_ctrls->txt_group_intra_lt_16x16    = 3;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
        txt_ctrls->early_exit_dist_th          = 400;
        txt_ctrls->early_exit_coeff_th         = 8;
        break;
    default: assert(0); break;
    }
}
void set_interpolation_search_level_ctrls(ModeDecisionContext *context_ptr,
                                          uint8_t              interpolation_search_level) {
    InterpolationSearchCtrls *ifs_ctrls = &context_ptr->ifs_ctrls;

    switch (interpolation_search_level) {
    case 0:
        ifs_ctrls->level                 = IFS_OFF;
        ifs_ctrls->quarter_pel_only      = 0;
        ifs_ctrls->early_skip            = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 1:
        ifs_ctrls->level                 = IFS_MDS0;
        ifs_ctrls->quarter_pel_only      = 0;
        ifs_ctrls->early_skip            = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 2:
        ifs_ctrls->level                 = IFS_MDS1;
        ifs_ctrls->quarter_pel_only      = 0;
        ifs_ctrls->early_skip            = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 3:
        ifs_ctrls->level                 = IFS_MDS2;
        ifs_ctrls->quarter_pel_only      = 0;
        ifs_ctrls->early_skip            = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 4:
        ifs_ctrls->level                 = IFS_MDS3;
        ifs_ctrls->quarter_pel_only      = 0;
        ifs_ctrls->early_skip            = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 5:
        ifs_ctrls->level                 = IFS_MDS3;
        ifs_ctrls->quarter_pel_only      = 0;
        ifs_ctrls->early_skip            = 1;
        ifs_ctrls->subsampled_distortion = 1;
        ifs_ctrls->skip_sse_rd_model     = 1;
        break;
    case 6:
        ifs_ctrls->level                 = IFS_MDS3;
        ifs_ctrls->quarter_pel_only      = 1;
        ifs_ctrls->early_skip            = 1;
        ifs_ctrls->subsampled_distortion = 1;
        ifs_ctrls->skip_sse_rd_model     = 1;
        break;
    default: assert(0); break;
    }
}
void set_cand_reduction_ctrls(PictureControlSet *pcs_ptr, ModeDecisionContext *mdctxt,
                              uint8_t cand_reduction_level, const uint32_t picture_qp,
                              uint32_t me_8x8_cost_variance, uint32_t me_64x64_distortion,
                              uint8_t l0_was_skip, uint8_t l1_was_skip, uint8_t ref_skip_perc) {
    CandReductionCtrls *cand_reduction_ctrls = &mdctxt->cand_reduction_ctrls;

    switch (cand_reduction_level) {
    case 0:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 0;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 0;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 3;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled = 0;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 0;

        break;

    case 1:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 3;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 0;

        break;

    case 2:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 1;

        break;

    case 3:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th   = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 2;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 1;

        break;

    case 4:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th   = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 1;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 2;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates =
            (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ||
             ((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
              me_8x8_cost_variance < (500 * picture_qp) &&
              me_64x64_distortion < (500 * picture_qp)))
            ? 3
            : 1;

        break;

    case 5:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th   = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 0;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 1;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 2;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates =
            (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ||
             ((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
              me_8x8_cost_variance < (500 * picture_qp) &&
              me_64x64_distortion < (500 * picture_qp)))
            ? 3
            : 1;

        break;

    default: assert(0); break;
    }

    // lpd1_mvp_best_me_list can only use this feature when a single unipred ME candidate is selected,
    if (!(pcs_ptr->parent_pcs_ptr->ref_list0_count_try == 1 &&
          pcs_ptr->parent_pcs_ptr->ref_list1_count_try == 1 &&
          pcs_ptr->parent_pcs_ptr->use_best_me_unipred_cand_only))
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;
}
void set_chroma_controls(ModeDecisionContext *mdctxt, uint8_t uv_level) {
    UvCtrls *uv_ctrls = &mdctxt->uv_ctrls;

    switch (uv_level) {
    case 0:
        uv_ctrls->enabled           = 0;
        uv_ctrls->uv_mode           = CHROMA_MODE_2;
        uv_ctrls->nd_uv_serach_mode = 0;
        break;
    case 1:
        uv_ctrls->enabled           = 1;
        uv_ctrls->uv_mode           = CHROMA_MODE_0;
        uv_ctrls->nd_uv_serach_mode = 0;
        uv_ctrls->uv_intra_th       = (uint64_t)~0;
        uv_ctrls->uv_cfl_th         = (uint64_t)~0;
        break;
    case 2:
        uv_ctrls->enabled           = 1;
        uv_ctrls->uv_mode           = CHROMA_MODE_0;
        uv_ctrls->nd_uv_serach_mode = 1;
        uv_ctrls->uv_intra_th       = 130;
        uv_ctrls->uv_cfl_th         = 130;
        break;
    case 3:
        uv_ctrls->enabled           = 1;
        uv_ctrls->uv_mode           = CHROMA_MODE_0;
        uv_ctrls->nd_uv_serach_mode = 1;
        uv_ctrls->uv_intra_th       = 100;
        uv_ctrls->uv_cfl_th         = 100;
        break;
    case 4:
        uv_ctrls->enabled           = 1;
        uv_ctrls->uv_mode           = CHROMA_MODE_1;
        uv_ctrls->nd_uv_serach_mode = 0;
        break;
    default: assert(0); break;
    }
}
void set_wm_controls(ModeDecisionContext *mdctxt, uint8_t wm_level) {
    WmCtrls *wm_ctrls = &mdctxt->wm_ctrls;

    switch (wm_level) {
    case 0: wm_ctrls->enabled = 0; break;
    case 1:
        wm_ctrls->enabled               = 1;
        wm_ctrls->use_wm_for_mvp        = 1;
        wm_ctrls->num_new_mv_refinement = 12;
        break;
    case 2:
        wm_ctrls->enabled               = 1;
        wm_ctrls->use_wm_for_mvp        = 0;
        wm_ctrls->num_new_mv_refinement = 0;
        break;
    default: assert(0); break;
    }
}
// Get the nic_level used for each preset (to be passed to setting function: set_nic_controls())
uint8_t get_nic_level(EbEncMode enc_mode, uint8_t temporal_layer_index) {
    uint8_t nic_level;

    if (enc_mode <= ENC_MRS)
        nic_level = 0;
    else if (enc_mode <= ENC_MR)
        nic_level = 1;
    else if (enc_mode <= ENC_M0)
        nic_level = (temporal_layer_index == 0) ? 2 : 4;
    else if (enc_mode <= ENC_M1)
        nic_level = 5;
    else if (enc_mode <= ENC_M2)
        nic_level = 9;
    else if (enc_mode <= ENC_M3)
        nic_level = 10;
    else if (enc_mode <= ENC_M4)
        nic_level = 11;
    else if (enc_mode <= ENC_M5)
        nic_level = 12;
    else if (enc_mode <= ENC_M7)
        nic_level = 14;
    else if (enc_mode <= ENC_M11)
        nic_level = 15;
    else
        nic_level = 16;

    return nic_level;
}
/*
* Set the NIC scaling and pruning controls.
*
* This function is used in MD to set the NIC controls and is also used at memory allocation
* to allocate the candidate buffers.  Therefore, the function returns the nic_scaling_level
* (index into MD_STAGE_NICS_SCAL_NUM array).
*
* When called at memory allocation, there is no context (it is passed as NULL) so the signals
* are not set.
*/
uint8_t set_nic_controls(ModeDecisionContext *ctx, uint8_t nic_level) {
    NicPruningCtrls *nic_pruning_ctrls = ctx ? &ctx->nic_ctrls.pruning_ctrls : NULL;
    uint8_t          nic_scaling_level = 0;
    uint8_t          md_staging_mode   = MD_STAGING_MODE_0;
    switch (nic_level) {
    case 0: // MAX NIC scaling; no pruning
        // NIC scaling level
        nic_scaling_level = 0;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_class_th = (uint64_t)~0;
            nic_pruning_ctrls->mds3_class_th = (uint64_t)~0;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_cand_base_th = (uint64_t)~0;
            nic_pruning_ctrls->mds3_cand_base_th = (uint64_t)~0;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 1:
        // NIC scaling level
        nic_scaling_level = 0;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 4;

            nic_pruning_ctrls->mds3_class_th = 25;
            nic_pruning_ctrls->mds3_band_cnt = 4;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_cand_base_th = 50;
            nic_pruning_ctrls->mds3_cand_base_th = 50;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 2:
        // NIC scaling level
        nic_scaling_level = 1;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 4;

            nic_pruning_ctrls->mds3_class_th = 25;
            nic_pruning_ctrls->mds3_band_cnt = 4;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_cand_base_th = 50;
            nic_pruning_ctrls->mds3_cand_base_th = 50;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 3:
        // NIC scaling level
        nic_scaling_level = 1;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 4;

            nic_pruning_ctrls->mds3_class_th = 25;
            nic_pruning_ctrls->mds3_band_cnt = 8;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 500;
            nic_pruning_ctrls->mds2_cand_base_th = 30;
            nic_pruning_ctrls->mds3_cand_base_th = 30;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 4:
        // NIC scaling level
        nic_scaling_level = 2;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;
            nic_pruning_ctrls->mds1_band_cnt = 2;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 4;

            nic_pruning_ctrls->mds3_class_th = 25;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 20;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 5:
        // NIC scaling level
        nic_scaling_level = 2;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 6:
        // NIC scaling level
        nic_scaling_level = 3;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 7:
        // NIC scaling level
        nic_scaling_level = 4;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 8:
        // NIC scaling level
        nic_scaling_level = 5;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 9:
        // NIC scaling level
        nic_scaling_level = 6;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 10:
        // NIC scaling level
        nic_scaling_level = 6;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 300;
            nic_pruning_ctrls->mds1_band_cnt = 4;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 11:
        // NIC scaling level
        nic_scaling_level = 8;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 300;
            nic_pruning_ctrls->mds1_band_cnt = 4;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 12:
        // NIC scaling level
        nic_scaling_level = 11;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 300;
            nic_pruning_ctrls->mds1_band_cnt = 4;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 13:
        // NIC scaling level
        nic_scaling_level = 11;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 200;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 200;
            nic_pruning_ctrls->mds2_cand_base_th = 15;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 14:
        // NIC scaling level
        nic_scaling_level = 12;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 200;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 200;
            nic_pruning_ctrls->mds2_cand_base_th = 15;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 15:
        // NIC scaling level
        nic_scaling_level = 14;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 200;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 10;
            nic_pruning_ctrls->mds2_band_cnt = 2;

            nic_pruning_ctrls->mds3_class_th = 10;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 50;
            nic_pruning_ctrls->mds2_cand_base_th = 5;
            nic_pruning_ctrls->mds3_cand_base_th = 5;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 16:
        // NIC scaling level
        nic_scaling_level = 15;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 75;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 0;
            nic_pruning_ctrls->mds2_band_cnt = 2;

            nic_pruning_ctrls->mds3_class_th = 0;
            nic_pruning_ctrls->mds3_band_cnt = 2;

            nic_pruning_ctrls->enable_skipping_mds1 = 1;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 1;
            nic_pruning_ctrls->mds2_cand_base_th = 1;
            nic_pruning_ctrls->mds3_cand_base_th = 1;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 17:
        // NIC scaling level
        nic_scaling_level = 15;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 75;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 0;
            nic_pruning_ctrls->mds2_band_cnt = 2;

            nic_pruning_ctrls->mds3_class_th = 0;
            nic_pruning_ctrls->mds3_band_cnt = 2;

            nic_pruning_ctrls->enable_skipping_mds1 = 1;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 1;
            nic_pruning_ctrls->mds2_cand_base_th = 1;
            nic_pruning_ctrls->mds3_cand_base_th = 1;
        }
        md_staging_mode = MD_STAGING_MODE_0;
        break;
    default: assert(0); break;
    }

    if (ctx) {
        NicScalingCtrls *nic_scaling_ctrls = &ctx->nic_ctrls.scaling_ctrls;
        nic_scaling_ctrls->stage1_scaling_num =
            MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_1];
        nic_scaling_ctrls->stage2_scaling_num =
            MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_2];
        nic_scaling_ctrls->stage3_scaling_num =
            MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_3];
        ctx->nic_ctrls.md_staging_mode = md_staging_mode;
    }

    return nic_scaling_level;
}
void set_inter_intra_ctrls(ModeDecisionContext *mdctxt, uint8_t inter_intra_level) {
    InterIntraCompCtrls *ii_ctrls = &mdctxt->inter_intra_comp_ctrls;

    switch (inter_intra_level) {
    case 0: ii_ctrls->enabled = 0; break;
    case 1: ii_ctrls->enabled = 1; break;
    default: assert(0); break;
    }
}
void set_depth_ctrls(ModeDecisionContext *ctx, uint8_t depth_level) {
    DepthCtrls *depth_ctrls = &ctx->depth_ctrls;

    switch (depth_level) {
    case 0:
        depth_ctrls->s_depth = 0;
        depth_ctrls->e_depth = 0;
        break;
    case 1:
        depth_ctrls->s_depth = -2;
        depth_ctrls->e_depth = 2;
        break;
    case 2:
        depth_ctrls->s_depth = -1;
        depth_ctrls->e_depth = 1;
        break;
    default: assert(0); break;
    }
}
/*
* return the 4x4 level
Used by signal_derivation_enc_dec_kernel_oq and memory allocation
*/
uint8_t get_disallow_4x4(EbEncMode enc_mode, EB_SLICE slice_type) {
    uint8_t disallow_4x4;
    if (enc_mode <= ENC_M0)
        disallow_4x4 = EB_FALSE;
    else if (enc_mode <= ENC_M5)
        disallow_4x4 = (slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
    else
        disallow_4x4 = EB_TRUE;

    return disallow_4x4;
}
void set_lpd1_ctrls(ModeDecisionContext *ctx, uint8_t lpd1_lvl) {
    Lpd1Ctrls *ctrls = &ctx->lpd1_ctrls;
    switch (lpd1_lvl) {
    case 0:
        ctrls->pd1_level = REGULAR_PD1; // Light-PD1 path not used
        break;
    case 1:
        ctrls->pd1_level = LPD1_LVL_0;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 128;
        ctrls->coeff_th[LPD1_LVL_0]                = 50;
        ctrls->max_mv_length[LPD1_LVL_0]           = 300;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = 250000;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 1024;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 1;
        break;
    case 2:
        ctrls->pd1_level = LPD1_LVL_0;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256;
        ctrls->coeff_th[LPD1_LVL_0]                = 100;
        ctrls->max_mv_length[LPD1_LVL_0]           = 500;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 1024;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 1;
        break;
    case 3:
        ctrls->pd1_level = LPD1_LVL_1;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256 << 9;
        ctrls->coeff_th[LPD1_LVL_0]                = 8192;
        ctrls->max_mv_length[LPD1_LVL_0]           = 2048;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 16384;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 3;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1]       = 1;
        ctrls->use_ref_info[LPD1_LVL_1]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_1]            = 256 << 6;
        ctrls->coeff_th[LPD1_LVL_1]                = 2000;
        ctrls->max_mv_length[LPD1_LVL_1]           = 1600;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = 750000;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1]   = 16384;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1]       = 3;
        break;
    case 4:
        ctrls->pd1_level = LPD1_LVL_3;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_0]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_0]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 5;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1]       = 1;
        ctrls->use_ref_info[LPD1_LVL_1]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_1]            = 256 << 13;
        ctrls->coeff_th[LPD1_LVL_1]                = 8192 * 8;
        ctrls->max_mv_length[LPD1_LVL_1]           = 2048 * 8;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1]   = 16384 * 6;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1]       = 5;

        // Set LPD1 level 2 controls
        ctrls->use_lpd1_detector[LPD1_LVL_2]       = 1;
        ctrls->use_ref_info[LPD1_LVL_2]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_2]            = 256 << 13;
        ctrls->coeff_th[LPD1_LVL_2]                = 8192 * 8;
        ctrls->max_mv_length[LPD1_LVL_2]           = 2048 * 8;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_2] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_2]   = 16384 * 6;
        ctrls->skip_pd0_me_shift[LPD1_LVL_2]       = 5;

        // Set LPD1 level 3 controls
        ctrls->use_lpd1_detector[LPD1_LVL_3]       = 1;
        ctrls->use_ref_info[LPD1_LVL_3]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_3]            = 256 << 13;
        ctrls->coeff_th[LPD1_LVL_3]                = 8192 * 8;
        ctrls->max_mv_length[LPD1_LVL_3]           = 2048 * 8;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_3] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_3]   = 16384 * 6;
        ctrls->skip_pd0_me_shift[LPD1_LVL_3]       = 5;
        break;
    case 5:
        ctrls->pd1_level = LPD1_LVL_4;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_0]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_0]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 5;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1]       = 1;
        ctrls->use_ref_info[LPD1_LVL_1]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_1]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_1]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_1]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1]       = 5;

        // Set LPD1 level 2 controls
        ctrls->use_lpd1_detector[LPD1_LVL_2]       = 1;
        ctrls->use_ref_info[LPD1_LVL_2]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_2]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_2]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_2]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_2] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_2]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_2]       = 5;

        // Set LPD1 level 3 controls
        ctrls->use_lpd1_detector[LPD1_LVL_3]       = 1;
        ctrls->use_ref_info[LPD1_LVL_3]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_3]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_3]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_3]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_3] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_3]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_3]       = 5;

        // Set LPD1 level 4 controls
        ctrls->use_lpd1_detector[LPD1_LVL_4]       = 1;
        ctrls->use_ref_info[LPD1_LVL_4]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_4]            = 256 << 9;
        ctrls->coeff_th[LPD1_LVL_4]                = 8192;
        ctrls->max_mv_length[LPD1_LVL_4]           = 2048;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_4] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_4]   = 16384 * 2;
        ctrls->skip_pd0_me_shift[LPD1_LVL_4]       = 3;
        break;
    case 6:
        ctrls->pd1_level = LPD1_LVL_5;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_0]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_0]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = (uint16_t)~0;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1]       = 1;
        ctrls->use_ref_info[LPD1_LVL_1]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_1]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_1]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_1]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1]       = (uint16_t)~0;

        // Set LPD1 level 2 controls
        ctrls->use_lpd1_detector[LPD1_LVL_2]       = 1;
        ctrls->use_ref_info[LPD1_LVL_2]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_2]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_2]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_2]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_2] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_2]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_2]       = (uint16_t)~0;

        // Set LPD1 level 3 controls
        ctrls->use_lpd1_detector[LPD1_LVL_3]       = 1;
        ctrls->use_ref_info[LPD1_LVL_3]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_3]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_3]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_3]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_3] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_3]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_3]       = (uint16_t)~0;

        // Set LPD1 level 4 controls
        ctrls->use_lpd1_detector[LPD1_LVL_4]       = 1;
        ctrls->use_ref_info[LPD1_LVL_4]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_4]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_4]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_4]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_4] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_4]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_4]       = (uint16_t)~0;

        // Set LPD1 level 5 controls
        ctrls->use_lpd1_detector[LPD1_LVL_5]       = 1;
        ctrls->use_ref_info[LPD1_LVL_5]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_5]            = 256 << 15;
        ctrls->coeff_th[LPD1_LVL_5]                = 8192 * 16;
        ctrls->max_mv_length[LPD1_LVL_5]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_5] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_5]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_5]       = (uint16_t)~0;
        break;
    default: assert(0); break;
    }
}
// use this function to set the disallow_below_16x16 level and to set the accompanying enable_me_8x8 level
uint8_t get_disallow_below_16x16_picture_level(EbEncMode enc_mode, EbInputResolution resolution,
                                               EB_SLICE slice_type, uint8_t sc_class1,
                                               uint8_t is_used_as_reference_flag,
                                               uint8_t temporal_layer_index) {
    uint8_t disallow_below_16x16 = 0;

    if (sc_class1)
        disallow_below_16x16 = 0;
    else if (enc_mode <= ENC_M8)
        disallow_below_16x16 = 0;
    else if (enc_mode <= ENC_M9)
        if (resolution <= INPUT_SIZE_1080p_RANGE)
            disallow_below_16x16 = is_used_as_reference_flag ? 0 : 1;
        else {
            disallow_below_16x16 = (slice_type == I_SLICE) ? 0 : 1;
        }
    else if (enc_mode <= ENC_M11)
        disallow_below_16x16 = (resolution <= INPUT_SIZE_480p_RANGE)
            ? (temporal_layer_index == 0 ? 0 : 1)
            : ((slice_type == I_SLICE) ? 0 : 1);

    else
        disallow_below_16x16 = (slice_type == I_SLICE) ? 0 : 1;

    return disallow_below_16x16;
}
/*
 * Generate per-SB MD settings (do not change per-PD)
 */
EbErrorType signal_derivation_enc_dec_kernel_common(SequenceControlSet  *scs_ptr,
                                                    PictureControlSet   *pcs_ptr,
                                                    ModeDecisionContext *ctx) {
    EbErrorType return_error = EB_ErrorNone;

    EbEncMode enc_mode = pcs_ptr->enc_mode;

    // Level 0: pred depth only
    // Level 1: [-2, +2] depth refinement
    // Level 2: [-1, +1] depth refinement
    uint8_t depth_level = 0;
    if (enc_mode <= ENC_MRS)
        depth_level = 1;
    else if (pcs_ptr->parent_pcs_ptr->sc_class1) {
        if (enc_mode <= ENC_M3)
            depth_level = pcs_ptr->slice_type == I_SLICE ? 1 : 2;
        else if (enc_mode <= ENC_M8)
            depth_level = 2;
        else
            depth_level = pcs_ptr->slice_type == I_SLICE ? 2 : 0;
    } else if (enc_mode <= ENC_M2)
        depth_level = pcs_ptr->slice_type == I_SLICE ? 1 : 2;
    else if (enc_mode <= ENC_M8)
        depth_level = 2;
    else
        depth_level = 0;
    set_depth_ctrls(ctx, depth_level);
    ctx->pred_depth_only = (depth_level == 0);
    ctx->pd0_level       = pcs_ptr->pic_pd0_level;
    SbParams *sb_params  = &pcs_ptr->parent_pcs_ptr->sb_params_array[ctx->sb_index];
    ctx->depth_removal_ctrls.disallow_below_64x64 = 0;
    ctx->depth_removal_ctrls.disallow_below_32x32 = 0;
    /*
if disallow_below_16x16 is turned ON then enable_me_8x8 should be turned OFF for the same preset in order to save memory and cycles as that feature optimizes the me_candidate_array,
me_mv_array and the total_me_candidate_index arrays when 8x8 blocks are not used

if any check other than an I-SLICE check is used on disallow_below_16x16 then the enable_me_8x8 should be turned ON for the entire preset because without the 8x8 me data the non I-SLICE pictures
that use 8x8 blocks will lose significant BD-Rate as the parent 16x16 me data will be used for the 8x8 blocks
*/
    ctx->depth_removal_ctrls.disallow_below_16x16 = pcs_ptr->pic_disallow_below_16x16;

    if (sb_params->width % 32 != 0 || sb_params->height % 32 != 0)
        ctx->depth_removal_ctrls.disallow_below_64x64 = EB_FALSE;
    if (sb_params->width % 16 != 0 || sb_params->height % 16 != 0)
        ctx->depth_removal_ctrls.disallow_below_32x32 = EB_FALSE;
    if (sb_params->width % 8 != 0 || sb_params->height % 8 != 0)
        ctx->depth_removal_ctrls.disallow_below_16x16 = EB_FALSE;

    // me_distortion/variance generated for 64x64 blocks only
    if (scs_ptr->super_block_size == 64) {
        set_depth_removal_level_controls(pcs_ptr, ctx, pcs_ptr->pic_depth_removal_level);
    }
    if (/*scs_ptr->rc_stat_gen_pass_mode || */ ctx->skip_pd0) {
        // ctx->depth_removal_ctrls.disallow_below_64x64 = 1;
        ctx->depth_removal_ctrls.disallow_below_32x32 = 1;
        //ctx->depth_removal_ctrls.disallow_below_16x16 = 1;
        if (sb_params->width % 32 != 0 || sb_params->height % 32 != 0) {
            ctx->depth_removal_ctrls.disallow_below_64x64 = EB_FALSE;
            ctx->depth_removal_ctrls.disallow_below_32x32 = EB_FALSE;
        }
        if (sb_params->width % 16 != 0 || sb_params->height % 16 != 0)
            ctx->depth_removal_ctrls.disallow_below_32x32 = EB_FALSE;
        else
            ctx->depth_removal_ctrls.enabled = 1;
        if (sb_params->width % 8 != 0 || sb_params->height % 8 != 0)
            ctx->depth_removal_ctrls.disallow_below_16x16 = EB_FALSE;
    }

    set_lpd1_ctrls(ctx, pcs_ptr->pic_lpd1_lvl);
    return return_error;
}
/*
 * Generate per-SB/per-PD MD settings
 */
void set_dist_based_ref_pruning_controls(ModeDecisionContext *mdctxt,
                                         uint8_t              dist_based_ref_pruning_level) {
    RefPruningControls *ref_pruning_ctrls = &mdctxt->ref_pruning_ctrls;

    switch (dist_based_ref_pruning_level) {
    case 0: ref_pruning_ctrls->enabled = 0; break;
    case 1:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = (uint32_t)~0;
        ref_pruning_ctrls->ref_idx_2_offset                     = 0;
        ref_pruning_ctrls->ref_idx_3_offset                     = 0;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 2:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 90;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 90;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 60;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 60;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 60;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 60;
        ref_pruning_ctrls->ref_idx_2_offset                     = 10;
        ref_pruning_ctrls->ref_idx_3_offset                     = 20;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 3:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 60;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 60;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 30;
        ref_pruning_ctrls->ref_idx_2_offset                     = 10;
        ref_pruning_ctrls->ref_idx_3_offset                     = 20;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 4:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 60;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 60;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 60;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = 30;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 10;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 10;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 10;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 10;
        ref_pruning_ctrls->ref_idx_2_offset                     = 10;
        ref_pruning_ctrls->ref_idx_3_offset                     = 20;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 5:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 30;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 30;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 10;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = 10;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 0;
        ref_pruning_ctrls->ref_idx_2_offset                     = 10;
        ref_pruning_ctrls->ref_idx_3_offset                     = 20;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;
        break;
    case 6:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 0;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 0;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 0;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = 0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = 0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 0;
        ref_pruning_ctrls->ref_idx_2_offset                     = 0;
        ref_pruning_ctrls->ref_idx_3_offset                     = 0;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;
        break;

    default: assert(0); break;
    }
}

void set_txs_controls(ModeDecisionContext *ctx, uint8_t txs_level) {
    TxsControls *txs_ctrls = &ctx->txs_ctrls;

    switch (txs_level) {
    case 0: txs_ctrls->enabled = 0; break;

    case 1:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 2;
        txs_ctrls->inter_class_max_depth   = 2;
        txs_ctrls->depth1_txt_group_offset = 0;
        txs_ctrls->depth2_txt_group_offset = 0;
        txs_ctrls->min_sq_size             = 0;
        break;

    case 2:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 2;
        txs_ctrls->inter_class_max_depth   = 1;
        txs_ctrls->depth1_txt_group_offset = 0;
        txs_ctrls->depth2_txt_group_offset = 0;
        txs_ctrls->min_sq_size             = 0;
        break;

    case 3:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 2;
        txs_ctrls->inter_class_max_depth   = 0;
        txs_ctrls->depth1_txt_group_offset = 0;
        txs_ctrls->depth2_txt_group_offset = 0;
        txs_ctrls->min_sq_size             = 0;
        break;

    case 4:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 1;
        txs_ctrls->inter_class_max_depth   = 0;
        txs_ctrls->depth1_txt_group_offset = 4;
        txs_ctrls->depth2_txt_group_offset = 4;
        txs_ctrls->min_sq_size             = 0;
        break;
    case 5:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 1;
        txs_ctrls->inter_class_max_depth   = 0;
        txs_ctrls->depth1_txt_group_offset = 4;
        txs_ctrls->depth2_txt_group_offset = 4;
        txs_ctrls->min_sq_size             = 32;
        break;
    default: assert(0); break;
    }
}
void set_spatial_sse_full_loop_level(ModeDecisionContext *ctx,
                                     uint8_t              spatial_sse_full_loop_level) {
    SpatialSSECtrls *spatial_sse_ctrls = &ctx->spatial_sse_ctrls;

    switch (spatial_sse_full_loop_level) {
    case 0: spatial_sse_ctrls->spatial_sse_full_loop_level = EB_FALSE; break;
    case 1: spatial_sse_ctrls->spatial_sse_full_loop_level = EB_TRUE; break;
    default: assert(0); break;
    }
}
// Compute a qp-aware threshold based on the variance of the SB, used to apply selectively subres
uint64_t compute_subres_th(SequenceControlSet *scs, PictureControlSet *pcs,
                           ModeDecisionContext *ctx) {
    uint32_t fast_lambda   = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD]
                                                    : ctx->fast_lambda_md[EB_8_BIT_MD];
    uint32_t sb_size       = scs->super_block_size * scs->super_block_size;
    uint64_t cost_th_rate  = 1 << 13;
    uint64_t use_subres_th = 0;

    if (scs->calculate_variance) {
        if (pcs->parent_pcs_ptr->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 400)
            use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 8);
        else if (pcs->parent_pcs_ptr->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 800)
            use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 7);
        else
            use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 6);
    } else {
        use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 6);
    }
    return use_subres_th;
}
// Compute a qp-aware threshold based on the variance of the SB, used to apply selectively apply PF
uint64_t compute_pf_th(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx) {
    uint32_t fast_lambda  = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD]
                                                   : ctx->fast_lambda_md[EB_8_BIT_MD];
    uint32_t sb_size      = scs->super_block_size * scs->super_block_size;
    uint64_t cost_th_rate = 1 << 13;
    uint64_t use_pf_th    = 0;

    if (scs->calculate_variance) {
        if (pcs->parent_pcs_ptr->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 400)
            use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        else if (pcs->parent_pcs_ptr->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 800)
            use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        else
            use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
    } else {
        use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
    }

    return use_pf_th;
}

void set_lpd1_tx_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t lpd1_tx_level) {
    PictureParentControlSet *ppcs  = pcs->parent_pcs_ptr;
    Lpd1TxCtrls             *ctrls = &ctx->lpd1_tx_ctrls;

    switch (lpd1_tx_level) {
    case 0:
        ctrls->zero_y_coeff_exit            = 0;
        ctrls->skip_nrst_nrst_luma_tx       = 0;
        ctrls->skip_tx_th                   = 0;
        ctrls->use_uv_shortcuts_on_y_coeffs = 0;

        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->use_neighbour_info    = 0;
        break;
    case 1:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = 1;
        ctrls->skip_nrst_nrst_luma_tx       = 0;
        ctrls->skip_tx_th                   = 0;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info    = 0;
        break;
    case 2:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = 1;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 25;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info    = 0;
        break;
    case 3:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = ppcs->is_used_as_reference_flag ? 2 : 3;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 25;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info    = 0;
        break;
    case 4:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = ppcs->is_used_as_reference_flag ? 2 : 3;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 25;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info    = 1;
        break;
    case 5:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = ppcs->is_used_as_reference_flag ? 2 : 3;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 50;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info    = 2;
        break;
    case 6:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = ppcs->is_used_as_reference_flag ? 2 : 3;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 70;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 100;
        ctrls->use_neighbour_info    = 2;
        break;
    default: assert(0); break;
    }
}
void set_cfl_ctrls(ModeDecisionContext *ctx, uint8_t cfl_level) {
    CflCtrls *ctrls = &ctx->cfl_ctrls;

    switch (cfl_level) {
    case 0: ctrls->enabled = 0; break;
    case 1:
        ctrls->enabled = 1;
        ctrls->itr_th  = 2;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->itr_th  = 1;
        break;
    default: assert(0); break;
    }
}
void set_rate_est_ctrls(ModeDecisionContext *ctx, uint8_t rate_est_level) {
    MdRateEstCtrls *ctrls = &ctx->rate_est_ctrls;

    switch (rate_est_level) {
    case 0:
        ctrls->update_skip_ctx_dc_sign_ctx = 0;
        ctrls->update_skip_coeff_ctx       = 0;
        ctrls->coeff_rate_est_lvl          = 0;
        ctrls->lpd0_qp_offset              = 8;
        ctrls->pd0_fast_coeff_est_level    = 2;
        break;
    case 1:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx       = 1;
        ctrls->coeff_rate_est_lvl          = 1;
        ctrls->lpd0_qp_offset              = 0;
        ctrls->pd0_fast_coeff_est_level    = 1;
        break;
    case 2:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx       = 0;
        ctrls->coeff_rate_est_lvl          = 1;
        ctrls->lpd0_qp_offset              = 0;
        ctrls->pd0_fast_coeff_est_level    = 2;
        break;
    case 3:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx       = 0;
        ctrls->coeff_rate_est_lvl          = 2;
        ctrls->lpd0_qp_offset              = 0;
        ctrls->pd0_fast_coeff_est_level    = 2;
        break;
    case 4:
        ctrls->update_skip_ctx_dc_sign_ctx = 0;
        ctrls->update_skip_coeff_ctx       = 0;
        ctrls->coeff_rate_est_lvl          = 2;
        ctrls->lpd0_qp_offset              = 0;
        ctrls->pd0_fast_coeff_est_level    = 2;
        break;
    default: assert(0); break;
    }
}
void set_intra_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t intra_level) {
    IntraCtrls *ctrls = &ctx->intra_ctrls;

    // If intra is disallowed at the pic level, must disallow at SB level
    if (pcs->skip_intra)
        intra_level = 0;

    assert(IMPLIES(pcs->slice_type == I_SLICE, intra_level > 0));

    switch (intra_level) {
    case 0:
        ctrls->enable_intra       = 0;
        ctrls->intra_mode_end     = DC_PRED;
        ctrls->angular_pred_level = 0;
        break;
    case 1:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = PAETH_PRED;
        ctrls->angular_pred_level = 1;
        break;
    case 2:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = PAETH_PRED;
        ctrls->angular_pred_level = 2;
        break;
    case 3:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = SMOOTH_H_PRED;
        ctrls->angular_pred_level = 3;
        break;
    case 4:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = SMOOTH_PRED;
        ctrls->angular_pred_level = 4;
        break;
    case 5:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = DC_PRED;
        ctrls->angular_pred_level = 0;
        break;
    default: assert(0); break;
    }

    /* For PD1, the ability to skip intra must be set at the pic level to ensure all SBs
    perform inverse TX and generate the recon. */
    if (ctx->pd_pass == PD_PASS_1) {
        ctx->skip_intra = pcs->skip_intra;

        // Check user-defined settings
        if (pcs->parent_pcs_ptr->scs_ptr->enable_paeth == 0)
            ctrls->intra_mode_end = MIN(ctrls->intra_mode_end, SMOOTH_H_PRED);

        if (pcs->parent_pcs_ptr->scs_ptr->enable_smooth == 0)
            ctrls->intra_mode_end = MIN(ctrls->intra_mode_end, D67_PRED);

        if (pcs->parent_pcs_ptr->scs_ptr->intra_angle_delta == 0)
            ctrls->angular_pred_level = 0;
    } else {
        ctx->skip_intra = !(ctrls->enable_intra) || pcs->skip_intra;
    }
}

void set_tx_shortcut_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx,
                           uint8_t tx_shortcut_level) {
    PictureParentControlSet *ppcs  = pcs->parent_pcs_ptr;
    TxShortcutCtrls         *ctrls = &ctx->tx_shortcut_ctrls;

    switch (tx_shortcut_level) {
    case 0:
        ctrls->bypass_tx_when_zcoeff = 0;
        ctrls->apply_pf_on_coeffs    = 0;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->use_neighbour_info    = 0;
        ctrls->chroma_detector_level = 0;
        break;
    case 1:
        ctrls->bypass_tx_when_zcoeff = 1;
        ctrls->apply_pf_on_coeffs    = 0;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->chroma_detector_level = 1;
        ctrls->use_neighbour_info    = 0;
        break;
    case 2:
        ctrls->bypass_tx_when_zcoeff = 2;
        ctrls->apply_pf_on_coeffs    = 1;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->chroma_detector_level = 1;
        ctrls->use_neighbour_info    = 0;
        break;
    case 3:
        ctrls->bypass_tx_when_zcoeff = 2;
        ctrls->apply_pf_on_coeffs    = 1;
        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->chroma_detector_level = 1;
        ctrls->use_neighbour_info    = 0;
        break;
    case 4:
        ctrls->bypass_tx_when_zcoeff = 2;
        ctrls->apply_pf_on_coeffs    = 1;
        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->chroma_detector_level = ppcs->is_used_as_reference_flag ? 1 : 0;
        ctrls->use_neighbour_info    = 1;
        break;
    case 5:
        ctrls->bypass_tx_when_zcoeff = 2;
        ctrls->apply_pf_on_coeffs    = 1;
        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->chroma_detector_level = 0;
        ctrls->use_neighbour_info    = 1;
        break;
    default: assert(0); break;
    }

    // Chroma detector should be used in M12 and below (at least in REF frames) to prevent blurring artifacts in some clips
    if (tx_shortcut_level && ppcs->is_used_as_reference_flag && pcs->enc_mode <= ENC_M12)
        assert(ctrls->chroma_detector_level &&
               "Chroma detector should be used for ref frames in low presets to prevent blurring "
               "artifacts.");
}
void set_mds0_controls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t mds0_level) {
    Mds0Ctrls *ctrls = &ctx->mds0_ctrls;

    switch (mds0_level) {
    case 0:
        ctrls->mds0_dist_type               = MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th           = 0;
        break;
    case 1:
        ctrls->mds0_dist_type               = MDS0_SSD;
        ctrls->enable_cost_based_early_exit = 0;
        ctrls->mds0_distortion_th           = 0;
        break;
    case 2:
        ctrls->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR
                                                                               : MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 0;
        ctrls->mds0_distortion_th           = 0;
        break;
    case 3:
        ctrls->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR
                                                                               : MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th           = 50;
        break;
    case 4:
        ctrls->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR
                                                                               : MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th           = 0;
        break;
    default: assert(0); break;
    }
}

// Set signals used for light-pd0 path; only PD0 should call this function
// assumes NSQ OFF, no 4x4, no chroma, no TXT/TXS/RDOQ/SSSE, SB_64x64
EbErrorType signal_derivation_enc_dec_kernel_oq_light_pd0(SequenceControlSet  *scs,
                                                          PictureControlSet   *pcs,
                                                          ModeDecisionContext *ctx) {
    EbErrorType return_error = EB_ErrorNone;

    ctx->md_disallow_nsq         = 1;
    ctx->inject_inter_candidates = 1;

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    ctx->shut_fast_rate = EB_TRUE;

    uint8_t intra_level = 0;
    if (pcs->slice_type == I_SLICE || pcs->parent_pcs_ptr->transition_present)
        intra_level = 1;
    else if (ctx->pd0_level <= LIGHT_PD0_LVL1)
        intra_level = (pcs->temporal_layer_index == 0) ? 1 : 0;
    else
        intra_level = 0;
    set_intra_ctrls(pcs, ctx, intra_level);
    if (ctx->pd0_level == VERY_LIGHT_PD0) {
        // Modulate the inter-depth bias based on the QP and the temporal complexity of the SB
        // towards more split for low QPs or/and complex SBs,
        // and less split for high QPs or/and easy SBs (to compensate for the absence of the coeff rate)
        ctx->inter_depth_bias = 950 + pcs->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
        if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 1000)
            ctx->inter_depth_bias = ctx->inter_depth_bias - 150;
        else if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 500)
            ctx->inter_depth_bias = ctx->inter_depth_bias - 50;
        else if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 250)
            ctx->inter_depth_bias = ctx->inter_depth_bias + 50;
        else
            ctx->inter_depth_bias = ctx->inter_depth_bias + 150;
    } else {
        ctx->inter_depth_bias = 0;
    }
    if (ctx->pd0_level == VERY_LIGHT_PD0)
        return return_error;

    uint8_t mds0_level = 0;
    if (ctx->pd0_level <= LIGHT_PD0_LVL2)
        mds0_level = 4;
    else
        mds0_level = 0;

    set_mds0_controls(pcs, ctx, mds0_level);
    set_chroma_controls(ctx, 0 /*chroma off*/);

    uint8_t pf_level = 1;
    if (pcs->slice_type != I_SLICE && !pcs->parent_pcs_ptr->transition_present) {
        if (ctx->pd0_level <= LIGHT_PD0_LVL2) {
            pf_level = 1;
        } else {
            // Use ME distortion and variance detector to enable PF
            uint64_t use_pf_th   = compute_pf_th(scs, pcs, ctx);
            uint32_t fast_lambda = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD]
                                                          : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64  = RDCOST(
                fast_lambda, 0, pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index]);

            if (ctx->pd0_level <= LIGHT_PD0_LVL3) {
                pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
            } else {
                if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
                else {
                    pf_level = (cost_64x64 < use_pf_th)  ? 3
                        : (cost_64x64 < (2 * use_pf_th)) ? 2
                                                         : 1;

                    if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 16000 &&
                        ctx->depth_removal_ctrls.enabled &&
                        (ctx->depth_removal_ctrls.disallow_below_32x32 ||
                         ctx->depth_removal_ctrls.disallow_below_64x64)) {
                        if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] <
                                4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] <
                                6 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] <
                                4 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index])
                            pf_level = 3;
                        else if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] <
                                 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index])
                            pf_level = (cost_64x64 < (12 * use_pf_th)) ? 3 : 2;

                        pf_level = (cost_64x64 < (8 * use_pf_th)) ? 3 : MAX(2, pf_level);
                    } else if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 16000 &&
                               ctx->depth_removal_ctrls.enabled &&
                               ctx->depth_removal_ctrls.disallow_below_16x16) {
                        if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] <
                                4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] <
                                6 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] <
                                4 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index])
                            pf_level = (cost_64x64 < (4 * use_pf_th)) ? 3 : 2;

                        pf_level = (cost_64x64 < (2 * use_pf_th)) ? 3
                            : (cost_64x64 < (8 * use_pf_th))      ? MAX(2, pf_level)
                                                                  : pf_level;
                    }
                }
            }
        }
    } else {
        pf_level = 1;
    }
    set_pf_controls(ctx, pf_level);

    uint8_t subres_level;
    if (ctx->pd0_level <= LIGHT_PD0_LVL1) {
        subres_level = 0;
    } else {
        subres_level = 0;

        SbParams *sb_params_ptr = &pcs->parent_pcs_ptr->sb_params_array[ctx->sb_index];

        // The controls checks the deviation between: (1) the pred-to-src SAD of even rows and (2) the pred-to-src SAD of odd rows for each 64x64 to decide whether to use subres or not
        // then applies the result to the 64x64 block and to all children, therefore if incomplete 64x64 then shut subres
        if (sb_params_ptr->is_complete_sb) {
            // Use ME distortion and variance detector to enable subres
            uint64_t use_subres_th = compute_subres_th(scs, pcs, ctx);
            uint32_t fast_lambda   = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD]
                                                            : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64    = RDCOST(
                fast_lambda, 0, pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index]);

            if (ctx->pd0_level <= LIGHT_PD0_LVL2) {
                if (pcs->slice_type == I_SLICE || pcs->parent_pcs_ptr->transition_present)
                    subres_level = 1;
                else
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
            } else if (ctx->pd0_level <= LIGHT_PD0_LVL3) {
                if (pcs->slice_type == I_SLICE || pcs->parent_pcs_ptr->transition_present)
                    subres_level = 1;
                else if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
                else
                    subres_level = 2;
            } else {
                if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    subres_level = (ctx->depth_removal_ctrls.enabled &&
                                    (ctx->depth_removal_ctrls.disallow_below_16x16 ||
                                     ctx->depth_removal_ctrls.disallow_below_32x32 ||
                                     ctx->depth_removal_ctrls.disallow_below_64x64))
                        ? 2
                        : 1;
                else
                    subres_level = 2;
            }
        } else {
            if (ctx->pd0_level <= LIGHT_PD0_LVL2)
                subres_level = 0;
            else
                subres_level = pcs->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 2;
        }
    }
    set_subres_controls(ctx, subres_level);

    if (ctx->pd0_level <= LIGHT_PD0_LVL2)
        set_rate_est_ctrls(ctx, 2);
    else
        set_rate_est_ctrls(ctx, 0);

    // set at pic-level b/c feature depends on some pic-level initializations
    ctx->approx_inter_rate = 1;

    return return_error;
}

void signal_derivation_enc_dec_kernel_oq_light_pd1(PictureControlSet   *pcs_ptr,
                                                   ModeDecisionContext *context_ptr) {
    EbEncMode lpd1_level = context_ptr->lpd1_ctrls.pd1_level;

    // Get ref info, used to set some feature levels
    const uint32_t picture_qp           = pcs_ptr->picture_qp;
    uint32_t       me_8x8_cost_variance = (uint32_t)~0;
    uint32_t       me_64x64_distortion  = (uint32_t)~0;
    uint8_t        l0_was_skip = 0, l1_was_skip = 0;
    uint8_t        l0_was_64x64_mvp = 0, l1_was_64x64_mvp = 0;

    if (pcs_ptr->slice_type != I_SLICE) {
        me_8x8_cost_variance = pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[context_ptr->sb_index];
        me_64x64_distortion  = pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index];

        EbReferenceObject *ref_obj_l0 =
            (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
        l0_was_skip = ref_obj_l0->sb_skip[context_ptr->sb_index], l1_was_skip = 1;
        l0_was_64x64_mvp = ref_obj_l0->sb_64x64_mvp[context_ptr->sb_index], l1_was_64x64_mvp = 1;
        if (pcs_ptr->slice_type == B_SLICE) {
            EbReferenceObject *ref_obj_l1 =
                (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
            l1_was_skip      = ref_obj_l1->sb_skip[context_ptr->sb_index];
            l1_was_64x64_mvp = ref_obj_l1->sb_64x64_mvp[context_ptr->sb_index];
        }
    }
    uint8_t ref_skip_perc = pcs_ptr->ref_skip_percentage;

    // Set candidate reduction levels
    uint8_t cand_reduction_level = 0;
    if (pcs_ptr->slice_type == I_SLICE)
        cand_reduction_level = 0;
    else if (lpd1_level <= LPD1_LVL_0)
        cand_reduction_level = 2;
    else if (lpd1_level <= LPD1_LVL_2)
        cand_reduction_level = 3;
    else if (lpd1_level <= LPD1_LVL_3)
        cand_reduction_level = 4;
    else
        cand_reduction_level = 5;

    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
        cand_reduction_level = 5;

    set_cand_reduction_ctrls(pcs_ptr,
                             context_ptr,
                             cand_reduction_level,
                             picture_qp,
                             me_8x8_cost_variance,
                             me_64x64_distortion,
                             l0_was_skip,
                             l1_was_skip,
                             ref_skip_perc);

    if (lpd1_level <= LPD1_LVL_0)
        context_ptr->rdoq_level = 4;
    else if (lpd1_level <= LPD1_LVL_4)
        context_ptr->rdoq_level = 5;
    else
        context_ptr->rdoq_level = 0;
    set_rdoq_controls(context_ptr, context_ptr->rdoq_level);

    if (lpd1_level <= LPD1_LVL_0)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <=
                INPUT_SIZE_1080p_RANGE
            ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 9 : 11)
            : 15;
    else if (lpd1_level <= LPD1_LVL_1)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <=
                INPUT_SIZE_1080p_RANGE
            ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13)
            : 15;
    else if (lpd1_level <= LPD1_LVL_2) {
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <=
                INPUT_SIZE_1080p_RANGE
            ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 13 : 14)
            : 15;
    } else if (lpd1_level <= LPD1_LVL_3) {
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <=
                INPUT_SIZE_1080p_RANGE
            ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 13 : 14)
            : 15;
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 50) ||
             (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
            me_8x8_cost_variance < (200 * picture_qp) && me_64x64_distortion < (200 * picture_qp))
            context_ptr->md_subpel_me_level = 0;
    } else if (lpd1_level <= LPD1_LVL_4) {
        context_ptr->md_subpel_me_level = 15;

        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 30) ||
             (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
            me_8x8_cost_variance < (500 * picture_qp) && me_64x64_distortion < (500 * picture_qp))
            context_ptr->md_subpel_me_level = 0;
    } else {
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 15
                                                                                             : 0;

        if ((((l0_was_skip || l1_was_skip) && ref_skip_perc > 20) ||
             (l0_was_64x64_mvp || l1_was_64x64_mvp)) &&
            me_8x8_cost_variance < (800 * picture_qp) && me_64x64_distortion < (800 * picture_qp))
            context_ptr->md_subpel_me_level = 0;
    }
    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);

    uint8_t mds0_level = 0;
    if (lpd1_level <= LPD1_LVL_2)
        mds0_level = 4;
    else {
        mds0_level = 4;
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 40) ||
             (me_8x8_cost_variance < (250 * picture_qp) &&
              me_64x64_distortion < (250 * picture_qp))))
            mds0_level = 0;
    }

    set_mds0_controls(pcs_ptr, context_ptr, mds0_level);

    uint8_t lpd1_tx_level = 0;
    if (lpd1_level <= LPD1_LVL_2)
        lpd1_tx_level = 3;
    else if (lpd1_level <= LPD1_LVL_3) {
        lpd1_tx_level = 4;
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
             me_8x8_cost_variance < (800 * picture_qp) &&
             me_64x64_distortion < (800 * picture_qp)) ||
            (me_8x8_cost_variance < (100 * picture_qp) && me_64x64_distortion < (100 * picture_qp)))
            lpd1_tx_level = 6;
    } else {
        lpd1_tx_level = 5;
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
             me_8x8_cost_variance < (800 * picture_qp) &&
             me_64x64_distortion < (800 * picture_qp)) ||
            (me_8x8_cost_variance < (100 * picture_qp) && me_64x64_distortion < (100 * picture_qp)))
            lpd1_tx_level = 6;
    }

    set_lpd1_tx_ctrls(pcs_ptr, context_ptr, lpd1_tx_level);

    /* In modes below M13, only use level 1-3 for chroma detector, as more aggressive levels will cause
    blurring artifacts in certain clips.

    Do not test this signal in M12 and below during preset tuning.  This signal should be kept as an enc_mode check
    instead of and LPD1_LEVEL check to ensure that M12 and below do not use it.
    */
    if (pcs_ptr->enc_mode >= ENC_M13) {
        if (lpd1_level <= LPD1_LVL_4)
            context_ptr->lpd1_tx_ctrls.chroma_detector_level = 4;
        else
            context_ptr->lpd1_tx_ctrls.chroma_detector_level = 0;
    }

    /* In modes below M13, only skip non-NEAREST_NEAREST TX b/c skipping all inter TX will cause blocking artifacts
    in certain clips.  This signal is separated from the general lpd1_tx_ctrls (above) to avoid
    accidentally turning this on for modes below M13.

    Do not test this signal in M12 and below during preset tuning.  This signal should be kept as an enc_mode check
    instead of and LPD1_LEVEL check to ensure that M12 and below do not use it.
    */
    if (pcs_ptr->enc_mode <= ENC_M12)
        context_ptr->lpd1_skip_inter_tx_level = 0;
    else {
        assert(pcs_ptr->enc_mode >= ENC_M13 && "Only enable this feature for M13+");
        if (lpd1_level <= LPD1_LVL_4) {
            context_ptr->lpd1_skip_inter_tx_level = 1;
            if (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
                 me_8x8_cost_variance < (800 * picture_qp) &&
                 me_64x64_distortion < (800 * picture_qp)) ||
                (me_8x8_cost_variance < (100 * picture_qp) &&
                 me_64x64_distortion < (100 * picture_qp))) {
                context_ptr->lpd1_skip_inter_tx_level = 2;
            }
        } else {
            context_ptr->lpd1_skip_inter_tx_level =
                pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 2;
            if (((((l0_was_skip || l1_was_skip) && ref_skip_perc > 20) ||
                  (l0_was_64x64_mvp || l1_was_64x64_mvp)) &&
                 me_8x8_cost_variance < (2000 * picture_qp) &&
                 me_64x64_distortion < (2000 * picture_qp)) ||
                (me_8x8_cost_variance < (400 * picture_qp) &&
                 me_64x64_distortion < (400 * picture_qp))) {
                context_ptr->lpd1_skip_inter_tx_level = 2;
            }
        }
    }
    uint8_t rate_est_level = 0;
    if (lpd1_level <= LPD1_LVL_0)
        rate_est_level = 4;
    else
        rate_est_level = 0;
    set_rate_est_ctrls(context_ptr, rate_est_level);

    // If want to turn off approximating inter rate, must ensure that the approximation is also disabled
    // at the pic level (pcs_ptr->approx_inter_rate)
    context_ptr->approx_inter_rate = 1;

    uint8_t pf_level = 1;
    if (lpd1_level <= LPD1_LVL_4)
        pf_level = 1;
    else
        pf_level = 2;
    set_pf_controls(context_ptr, pf_level);

    uint8_t intra_level = 0;
    if (lpd1_level <= LPD1_LVL_2)
        intra_level = 4;
    else
        intra_level = 5;

    set_intra_ctrls(pcs_ptr, context_ptr, intra_level);

    /* Set signals that have assumed values in the light-PD1 path (but need to be initialized as they may be checked) */

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    context_ptr->shut_fast_rate             = EB_FALSE;
    context_ptr->uv_ctrls.enabled           = 1;
    context_ptr->uv_ctrls.uv_mode           = CHROMA_MODE_1;
    context_ptr->uv_ctrls.nd_uv_serach_mode = 0;
    set_cfl_ctrls(context_ptr, 0);
    context_ptr->md_disallow_nsq                            = pcs_ptr->parent_pcs_ptr->disallow_nsq;
    context_ptr->new_nearest_injection                      = 1;
    context_ptr->inject_inter_candidates                    = 1;
    context_ptr->blk_skip_decision                          = EB_TRUE;
    context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx = 0;
    context_ptr->rate_est_ctrls.update_skip_coeff_ctx       = 0;
    context_ptr->subres_ctrls.odd_to_even_deviation_th      = 0;
}
EbErrorType signal_derivation_enc_dec_kernel_oq(SequenceControlSet *scs, PictureControlSet *pcs_ptr,
                                                ModeDecisionContext *context_ptr) {
    EbErrorType              return_error         = EB_ErrorNone;
    EbEncMode                enc_mode             = pcs_ptr->enc_mode;
    uint8_t                  pd_pass              = context_ptr->pd_pass;
    PictureParentControlSet *ppcs                 = pcs_ptr->parent_pcs_ptr;
    const uint8_t            is_ref               = ppcs->is_used_as_reference_flag;
    const uint8_t            is_base              = ppcs->temporal_layer_index == 0;
    const EbInputResolution  input_resolution     = ppcs->input_resolution;
    const EB_SLICE           slice_type           = pcs_ptr->slice_type;
    const uint8_t            fast_decode          = scs->static_config.fast_decode;
    const uint32_t           picture_qp           = pcs_ptr->picture_qp;
    uint32_t                 me_8x8_cost_variance = (uint32_t)~0;
    uint32_t                 me_64x64_distortion  = (uint32_t)~0;
    uint8_t                  l0_was_skip = 0, l1_was_skip = 0;
    uint8_t                  ref_skip_perc = pcs_ptr->ref_skip_percentage;
    set_cand_reduction_ctrls(pcs_ptr,
                             context_ptr,
                             pd_pass == PD_PASS_0 ? 0 : pcs_ptr->cand_reduction_level,
                             picture_qp,
                             me_8x8_cost_variance,
                             me_64x64_distortion,
                             l0_was_skip,
                             l1_was_skip,
                             ref_skip_perc);
    set_txt_controls(context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->txt_level);

    set_tx_shortcut_ctrls(
        pcs_ptr, context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->tx_shortcut_level);

    set_interpolation_search_level_ctrls(
        context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->interpolation_search_level);

    set_chroma_controls(context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->chroma_level);

    set_cfl_ctrls(context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->cfl_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->md_disallow_nsq = enc_mode <= ENC_M0 ? pcs_ptr->parent_pcs_ptr->disallow_nsq
                                                          : 1;
    else
        // Update nsq settings based on the sb_class
        context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;

    if (pd_pass == PD_PASS_0)
        context_ptr->global_mv_injection = 0;
    else
        context_ptr->global_mv_injection = pcs_ptr->parent_pcs_ptr->gm_ctrls.enabled;

    if (pd_pass == PD_PASS_0)
        context_ptr->new_nearest_injection = 0;
    else
        context_ptr->new_nearest_injection = 1;
    context_ptr->new_nearest_near_comb_injection = pd_pass == PD_PASS_0
        ? 0
        : pcs_ptr->new_nearest_near_comb_injection;

    //set Warped-Motion controls from Picture level.
    set_wm_controls(context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->wm_level);

    context_ptr->unipred3x3_injection    = pd_pass == PD_PASS_0 ? 0 : pcs_ptr->unipred3x3_injection;
    context_ptr->bipred3x3_injection     = pd_pass == PD_PASS_0 ? 0 : pcs_ptr->bipred3x3_injection;
    context_ptr->inject_inter_candidates = 1;
    context_ptr->inter_compound_mode     = pd_pass == PD_PASS_0 ? 0 : pcs_ptr->inter_compound_mode;

    set_dist_based_ref_pruning_controls(context_ptr,
                                        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->dist_based_ref_pruning);

    set_spatial_sse_full_loop_level(
        context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->spatial_sse_full_loop_level);
    if (context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1)
        context_ptr->blk_skip_decision = EB_TRUE;
    else
        context_ptr->blk_skip_decision = EB_FALSE;
    if (pd_pass == PD_PASS_0)
        if (enc_mode <= ENC_M1)
            context_ptr->rdoq_level = 1;
        else
            context_ptr->rdoq_level = 0;
    else if (enc_mode <= ENC_M11)
        context_ptr->rdoq_level = 1;
    else
        context_ptr->rdoq_level = 5;
    set_rdoq_controls(context_ptr, context_ptr->rdoq_level);

    // Derive redundant block
    if (pd_pass == PD_PASS_0 || context_ptr->md_disallow_nsq)
        context_ptr->redundant_blk = EB_FALSE;
    else
        context_ptr->redundant_blk = EB_TRUE;

    set_parent_sq_coeff_area_based_cycles_reduction_ctrls(
        context_ptr,
        pcs_ptr->parent_pcs_ptr->input_resolution,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->parent_sq_coeff_area_based_cycles_reduction_level);
    context_ptr->sq_weight = pd_pass == PD_PASS_0 ? (uint32_t)~0 : pcs_ptr->sq_weight;

    context_ptr->max_part0_to_part1_dev = pd_pass == PD_PASS_0 ? 0
                                                               : pcs_ptr->max_part0_to_part1_dev;

    // Set pic_obmc_level @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_pic_obmc_level = 0;
    else
        context_ptr->md_pic_obmc_level = pcs_ptr->parent_pcs_ptr->pic_obmc_level;

    set_obmc_controls(context_ptr, context_ptr->md_pic_obmc_level);

    context_ptr->md_inter_intra_level = pcs_ptr->md_inter_intra_level;
    set_inter_intra_ctrls(context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->md_inter_intra_level);
    set_txs_controls(context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->txs_level);

    set_nic_controls(context_ptr, pd_pass == PD_PASS_0 ? 17 : pcs_ptr->nic_level);
    // Set md_filter_intra_mode @ MD
    // md_filter_intra_level specifies whether filter intra would be active
    // for a given prediction candidate in mode decision.

    // md_filter_intra_level | Settings
    // 0                      | OFF
    // 1                      | ON
    if (pd_pass == PD_PASS_0)
        context_ptr->md_filter_intra_level = 0;
    else
        context_ptr->md_filter_intra_level = pcs_ptr->pic_filter_intra_level;
    // Set md_allow_intrabc @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_allow_intrabc = 0;
    else
        context_ptr->md_allow_intrabc = pcs_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc;

    // Set md_palette_level @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_palette_level = 0;
    else
        context_ptr->md_palette_level = pcs_ptr->parent_pcs_ptr->palette_level;

    uint8_t pf_level = 1;
    set_pf_controls(context_ptr, pf_level);
    md_sq_motion_search_controls(context_ptr,
                                 pd_pass == PD_PASS_0 ? 0 : pcs_ptr->md_sq_mv_search_level);

    md_nsq_motion_search_controls(context_ptr,
                                  pd_pass == PD_PASS_0 ? 0 : pcs_ptr->md_nsq_mv_search_level);

    md_pme_search_controls(context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->md_pme_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->md_subpel_me_level = enc_mode <= ENC_M5 ? 3 : 0;
    else if (enc_mode <= ENC_M0)
        context_ptr->md_subpel_me_level = 1;
    else if (enc_mode <= ENC_M6)
        context_ptr->md_subpel_me_level = input_resolution <= INPUT_SIZE_480p_RANGE ? 1 : 2;
    else if (enc_mode <= ENC_M9)
        context_ptr->md_subpel_me_level = is_base ? 2 : (is_ref ? 4 : 7);
    else if (enc_mode <= ENC_M11)
        context_ptr->md_subpel_me_level = is_ref ? 4 : 7;
    else if (enc_mode <= ENC_M12)
        context_ptr->md_subpel_me_level = is_ref ? 9 : 11;
    else
        context_ptr->md_subpel_me_level = is_base ? 9 : 0;
    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->md_subpel_pme_level = enc_mode <= ENC_M0 ? 3 : 0;
    else if (enc_mode <= ENC_M0)
        context_ptr->md_subpel_pme_level = 1;
    else
        context_ptr->md_subpel_pme_level = 2;

    md_subpel_pme_controls(context_ptr, context_ptr->md_subpel_pme_level);
    uint8_t rate_est_level = 0;
    if (pd_pass == PD_PASS_0) {
        if (enc_mode <= ENC_M2)
            rate_est_level = 1;
        else
            rate_est_level = 2;
    } else if (fast_decode == 0) {
        if (enc_mode <= ENC_M3)
            rate_est_level = 1;
        else if (enc_mode <= ENC_M10)
            rate_est_level = 2;
        else if (enc_mode <= ENC_M12)
            rate_est_level = (slice_type == I_SLICE) ? 3 : 4;
        else
            rate_est_level = (slice_type == I_SLICE) ? 3 : 0;
    } else if (fast_decode <= 2) {
        if (enc_mode <= ENC_M3)
            rate_est_level = 1;
        else if (enc_mode <= ENC_M4)
            rate_est_level = 2;
        else if (enc_mode <= ENC_M10)
            rate_est_level = input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : 3;
        else if (enc_mode <= ENC_M12)
            rate_est_level = (slice_type == I_SLICE) ? 3 : 4;
        else
            rate_est_level = (slice_type == I_SLICE) ? 3 : 0;
    } else {
        if (enc_mode <= ENC_M3)
            rate_est_level = 1;
        else if (enc_mode <= ENC_M4)
            rate_est_level = 2;
        else if (enc_mode <= ENC_M10)
            rate_est_level = input_resolution <= INPUT_SIZE_360p_RANGE ? 2 : 3;
        else if (enc_mode <= ENC_M12)
            rate_est_level = (slice_type == I_SLICE) ? 3 : 4;
        else
            rate_est_level = (slice_type == I_SLICE) ? 3 : 0;
    }
    set_rate_est_ctrls(context_ptr, rate_est_level);

    // set at pic-level b/c feature depends on some pic-level initializations
    context_ptr->approx_inter_rate = pcs_ptr->approx_inter_rate;
    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    if (pd_pass == PD_PASS_0)
        context_ptr->shut_fast_rate = EB_TRUE;
    else
        context_ptr->shut_fast_rate = EB_FALSE;

    // intra_level must be greater than 0 for I_SLICE
    uint8_t intra_level = 0;
    if (pd_pass == PD_PASS_0) {
        if (pcs_ptr->slice_type == I_SLICE || pcs_ptr->parent_pcs_ptr->transition_present)
            intra_level = 5;
        else if (enc_mode <= ENC_M1)
            intra_level = 5;
        else
            intra_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 0;
    } else if (enc_mode <= ENC_M0)
        intra_level = 1;
    else if (enc_mode <= ENC_M2)
        intra_level = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 4;
    else if (enc_mode <= ENC_M5)
        intra_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 4;
    else if (enc_mode <= ENC_M7)
        intra_level = (pcs_ptr->slice_type == I_SLICE ||
                       pcs_ptr->parent_pcs_ptr->transition_present)
            ? 1
            : (pcs_ptr->temporal_layer_index == 0) ? 2
                                                   : 4;
    else if (enc_mode <= ENC_M8)
        intra_level = pcs_ptr->slice_type == I_SLICE ? 1 : 4;
    else
        intra_level = (pcs_ptr->slice_type == I_SLICE ||
                       pcs_ptr->parent_pcs_ptr->transition_present)
            ? 3
            : 4;

    set_intra_ctrls(pcs_ptr, context_ptr, intra_level);
    set_mds0_controls(pcs_ptr, context_ptr, pd_pass == PD_PASS_0 ? 2 : pcs_ptr->mds0_level);
    set_subres_controls(context_ptr, 0);
    context_ptr->inter_depth_bias = 0;
    return return_error;
}
void copy_neighbour_arrays_light_pd0(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                     uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds,
                                     uint32_t sb_org_x, uint32_t sb_org_y);
void copy_neighbour_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                           uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds, uint32_t sb_org_x,
                           uint32_t sb_org_y);

static void set_parent_to_be_considered(MdcSbData *results_ptr, uint32_t blk_index, int32_t sb_size,
                                        int8_t pred_depth, uint8_t pred_sq_idx,
                                        const uint8_t disallow_nsq, int8_t depth_step) {
    const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
    if (blk_geom->sq_size < ((sb_size == BLOCK_128X128) ? 128 : 64)) {
        //Set parent to be considered
        uint32_t         parent_depth_idx_mds = blk_geom->parent_depth_idx_mds;
        const BlockGeom *parent_blk_geom      = get_blk_geom_mds(parent_depth_idx_mds);
        const uint32_t   parent_tot_d1_blocks = disallow_nsq ? 1
              : parent_blk_geom->sq_size == 128              ? 17
              : parent_blk_geom->sq_size > 8                 ? 25
              : parent_blk_geom->sq_size == 8                ? 5
                                                             : 1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < parent_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[parent_depth_idx_mds + block_1d_idx] = 1;
        }

        if (depth_step < -1)
            set_parent_to_be_considered(results_ptr,
                                        parent_depth_idx_mds,
                                        sb_size,
                                        pred_depth,
                                        pred_sq_idx,
                                        disallow_nsq,
                                        depth_step + 1);
    }
}
static void set_child_to_be_considered(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                       MdcSbData *results_ptr, uint32_t blk_index,
                                       uint32_t sb_index, int32_t sb_size, int8_t pred_depth,
                                       uint8_t pred_sq_idx, int8_t depth_step) {
    const BlockGeom *blk_geom      = get_blk_geom_mds(blk_index);
    unsigned         tot_d1_blocks = blk_geom->sq_size == 128 ? 17
                : blk_geom->sq_size > 8                       ? 25
                : blk_geom->sq_size == 8                      ? 5
                                                              : 1;
    if (blk_geom->geom_idx == GEOM_0)
        tot_d1_blocks = 1;

    if (blk_geom->sq_size == 8 && context_ptr->disallow_4x4)
        return;
    if (blk_geom->sq_size > 4) {
        for (uint32_t block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[blk_index + block_1d_idx]     = 1;
            results_ptr->refined_split_flag[blk_index + block_1d_idx] = EB_TRUE;
        }
        //Set first child to be considered
        uint32_t         child_block_idx_1    = blk_index + blk_geom->d1_depth_offset;
        const BlockGeom *child1_blk_geom      = get_blk_geom_mds(child_block_idx_1);
        const uint32_t   child1_tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1
              : child1_blk_geom->sq_size == 128                                       ? 17
              : child1_blk_geom->sq_size > 8                                          ? 25
              : child1_blk_geom->sq_size == 8                                         ? 5
                                                                                      : 1;

        for (uint32_t block_1d_idx = 0; block_1d_idx < child1_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[child_block_idx_1 + block_1d_idx]     = 1;
            results_ptr->refined_split_flag[child_block_idx_1 + block_1d_idx] = EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 ||
            !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_1])
            set_child_to_be_considered(pcs_ptr,
                                       context_ptr,
                                       results_ptr,
                                       child_block_idx_1,
                                       sb_index,
                                       sb_size,
                                       pred_depth,
                                       pred_sq_idx,
                                       depth_step > 1 ? depth_step - 1 : 1);
        //Set second child to be considered
        uint32_t child_block_idx_2 = child_block_idx_1 +
            ns_depth_offset[blk_geom->geom_idx][blk_geom->depth + 1];
        const BlockGeom *child2_blk_geom      = get_blk_geom_mds(child_block_idx_2);
        const uint32_t   child2_tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1
              : child2_blk_geom->sq_size == 128                                       ? 17
              : child2_blk_geom->sq_size > 8                                          ? 25
              : child2_blk_geom->sq_size == 8                                         ? 5
                                                                                      : 1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < child2_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[child_block_idx_2 + block_1d_idx]     = 1;
            results_ptr->refined_split_flag[child_block_idx_2 + block_1d_idx] = EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 ||
            !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_2])
            set_child_to_be_considered(pcs_ptr,
                                       context_ptr,
                                       results_ptr,
                                       child_block_idx_2,
                                       sb_index,
                                       sb_size,
                                       pred_depth,
                                       pred_sq_idx,
                                       depth_step > 1 ? depth_step - 1 : 1);
        //Set third child to be considered
        uint32_t child_block_idx_3 = child_block_idx_2 +
            ns_depth_offset[blk_geom->geom_idx][blk_geom->depth + 1];
        const BlockGeom *child3_blk_geom      = get_blk_geom_mds(child_block_idx_3);
        const uint32_t   child3_tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1
              : child3_blk_geom->sq_size == 128                                       ? 17
              : child3_blk_geom->sq_size > 8                                          ? 25
              : child3_blk_geom->sq_size == 8                                         ? 5
                                                                                      : 1;

        for (uint32_t block_1d_idx = 0; block_1d_idx < child3_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[child_block_idx_3 + block_1d_idx]     = 1;
            results_ptr->refined_split_flag[child_block_idx_3 + block_1d_idx] = EB_FALSE;
        }

        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 ||
            !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_3])
            set_child_to_be_considered(pcs_ptr,
                                       context_ptr,
                                       results_ptr,
                                       child_block_idx_3,
                                       sb_index,
                                       sb_size,
                                       pred_depth,
                                       pred_sq_idx,
                                       depth_step > 1 ? depth_step - 1 : 1);
        //Set forth child to be considered
        uint32_t child_block_idx_4 = child_block_idx_3 +
            ns_depth_offset[blk_geom->geom_idx][blk_geom->depth + 1];
        const BlockGeom *child4_blk_geom      = get_blk_geom_mds(child_block_idx_4);
        const uint32_t   child4_tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1
              : child4_blk_geom->sq_size == 128                                       ? 17
              : child4_blk_geom->sq_size > 8                                          ? 25
              : child4_blk_geom->sq_size == 8                                         ? 5
                                                                                      : 1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < child4_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[child_block_idx_4 + block_1d_idx]     = 1;
            results_ptr->refined_split_flag[child_block_idx_4 + block_1d_idx] = EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 ||
            !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_4])
            set_child_to_be_considered(pcs_ptr,
                                       context_ptr,
                                       results_ptr,
                                       child_block_idx_4,
                                       sb_index,
                                       sb_size,
                                       pred_depth,
                                       pred_sq_idx,
                                       depth_step > 1 ? depth_step - 1 : 1);
    }
}
uint32_t get_tot_1d_blks(struct PictureParentControlSet *ppcs, const int32_t sq_size,
                         const uint8_t disallow_nsq) {
    uint32_t tot_d1_blocks;

    tot_d1_blocks = (disallow_nsq) ||
            (sq_size >= 64 && ppcs->disallow_all_nsq_blocks_above_64x64) ||
            (sq_size >= 32 && ppcs->disallow_all_nsq_blocks_above_32x32) ||
            (sq_size >= 16 && ppcs->disallow_all_nsq_blocks_above_16x16) ||
            (sq_size <= 64 && ppcs->disallow_all_nsq_blocks_below_64x64) ||
            (sq_size <= 32 && ppcs->disallow_all_nsq_blocks_below_32x32) ||
            (sq_size <= 8 && ppcs->disallow_all_nsq_blocks_below_8x8) ||
            (sq_size <= 16 && ppcs->disallow_all_nsq_blocks_below_16x16)
        ? 1
        : (sq_size == 16 && ppcs->disallow_all_non_hv_nsq_blocks_below_16x16) ? 5
        : (sq_size == 16 && ppcs->disallow_all_h4_v4_blocks_below_16x16)      ? 17
        : sq_size == 128                                                      ? 17
        : sq_size > 8                                                         ? 25
        : sq_size == 8                                                        ? 5
                                                                              : 1;

    if (ppcs->disallow_HVA_HVB_HV4)
        tot_d1_blocks = MIN(5, tot_d1_blocks);

    if (ppcs->disallow_HV4)
        tot_d1_blocks = MIN(17, tot_d1_blocks);

    return tot_d1_blocks;
}

EbErrorType rtime_alloc_palette_info(BlkStruct *md_blk_arr_nsq) {
    EB_MALLOC_ARRAY(md_blk_arr_nsq->palette_info, 1);
    EB_MALLOC_ARRAY(md_blk_arr_nsq->palette_info->color_idx_map, MAX_PALETTE_SQUARE);

    return EB_ErrorNone;
}
// Initialize structures used to indicate which blocks will be tested at MD.
// MD data structures should be updated in init_block_data(), not here.
static void build_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                   ModeDecisionContext *context_ptr, uint32_t sb_index,
                                   EbBool is_complete_sb) {
    memset(context_ptr->tested_blk_flag, 0, sizeof(uint8_t) * scs_ptr->max_block_cnt);
    memset(context_ptr->avail_blk_flag, EB_FALSE, sizeof(uint8_t) * scs_ptr->max_block_cnt);

    MdcSbData *results_ptr       = context_ptr->mdc_sb_array;
    results_ptr->leaf_count      = 0;
    uint32_t       blk_index     = 0;
    const uint16_t max_block_cnt = scs_ptr->max_block_cnt;
    int32_t        min_sq_size;
    if (context_ptr->pred_depth_only)
        min_sq_size = (context_ptr->depth_removal_ctrls.enabled &&
                       context_ptr->depth_removal_ctrls.disallow_below_64x64)
            ? 64
            : (context_ptr->depth_removal_ctrls.enabled &&
               context_ptr->depth_removal_ctrls.disallow_below_32x32)
            ? 32
            : (context_ptr->depth_removal_ctrls.enabled &&
               context_ptr->depth_removal_ctrls.disallow_below_16x16)
            ? 16
            : context_ptr->disallow_4x4 ? 8
                                        : 4;
    else
        min_sq_size = context_ptr->disallow_4x4 ? 8 : 4;

    while (blk_index < max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

        // Initialize here because may not be updated at inter-depth decision for incomplete SBs
        if (!is_complete_sb)
            context_ptr->md_blk_arr_nsq[blk_index].part = PARTITION_SPLIT;

        // SQ/NSQ block(s) filter based on the SQ size
        uint8_t is_block_tagged = (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE) ||
                (context_ptr->pred_depth_only && (blk_geom->sq_size < min_sq_size))
            ? 0
            : 1;
        if (context_ptr->skip_pd0)
            is_block_tagged = !context_ptr->depth_removal_ctrls.disallow_below_64x64 &&
                    context_ptr->depth_removal_ctrls.disallow_below_32x32 &&
                    (blk_geom->sq_size != 32)
                ? 0
                : is_block_tagged;

        // SQ/NSQ block(s) filter based on the block validity
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            is_block_tagged) {
            const uint32_t tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq
                ? 1
                : get_tot_1d_blks(
                      pcs_ptr->parent_pcs_ptr, blk_geom->sq_size, context_ptr->md_disallow_nsq);

            for (uint32_t idx = blk_index; idx < (tot_d1_blocks + blk_index); ++idx) {
                //  MD palette info buffer
                if (pcs_ptr->parent_pcs_ptr->palette_level) {
                    if (context_ptr->md_blk_arr_nsq[idx].palette_mem == 0) {
                        rtime_alloc_palette_info(&context_ptr->md_blk_arr_nsq[idx]);
                        context_ptr->md_blk_arr_nsq[idx].palette_mem = 1;
                    }
                }

                context_ptr->md_blk_arr_nsq[idx].palette_size[0] = 0;
                context_ptr->md_blk_arr_nsq[idx].palette_size[1] = 0;

                if (results_ptr->consider_block[idx]) {
                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = idx;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks =
                        tot_d1_blocks;
                    results_ptr->split_flag[results_ptr->leaf_count++] =
                        results_ptr->refined_split_flag[idx];
                }
            }
            blk_index += blk_geom->d1_depth_offset;
        } else {
            blk_index += (blk_geom->sq_size > min_sq_size) ? blk_geom->d1_depth_offset
                                                           : blk_geom->ns_depth_offset;
        }
    }
}
void update_pred_th_offset(ModeDecisionContext *ctx, const BlockGeom *blk_geom, int8_t *s_depth,
                           int8_t *e_depth, int64_t *th_offset) {
    uint32_t full_lambda = ctx->hbd_mode_decision ? ctx->full_lambda_md[EB_10_BIT_MD]
                                                  : ctx->full_lambda_md[EB_8_BIT_MD];

    // cost-band-based modulation
    uint64_t max_cost = RDCOST(
        full_lambda,
        16,
        ctx->depth_refinement_ctrls.max_cost_multiplier * blk_geom->bwidth * blk_geom->bheight);

    if (ctx->md_local_blk_unit[blk_geom->sqi_mds].default_cost <= max_cost) {
        uint64_t band_size = max_cost / ctx->depth_refinement_ctrls.max_band_cnt;
        uint64_t band_idx  = ctx->md_local_blk_unit[blk_geom->sqi_mds].default_cost / band_size;
        if (ctx->depth_refinement_ctrls.decrement_per_band[band_idx] == MAX_SIGNED_VALUE) {
            *s_depth = 0;
            *e_depth = 0;
        } else {
            *th_offset = -ctx->depth_refinement_ctrls.decrement_per_band[band_idx];
        }
    } else {
        *th_offset = 0;
    }
}
uint8_t is_parent_to_current_deviation_small(ModeDecisionContext *mdctxt, const BlockGeom *blk_geom,
                                             int64_t th_offset) {
    if (mdctxt->depth_refinement_ctrls.parent_to_current_th == MIN_SIGNED_VALUE)
        return EB_FALSE;
    // block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
    // Get the parent of the current block
    uint32_t parent_depth_idx_mds = blk_geom->parent_depth_idx_mds;
    if (mdctxt->avail_blk_flag[parent_depth_idx_mds]) {
        mdctxt->parent_to_current_deviation =
            (int64_t)(((int64_t)MAX(mdctxt->md_local_blk_unit[parent_depth_idx_mds].default_cost,
                                    1) -
                       (int64_t)MAX((mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost * 4),
                                    1)) *
                      100) /
            (int64_t)MAX((mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost * 4), 1);
    }
    if (mdctxt->parent_to_current_deviation <=
        (mdctxt->depth_refinement_ctrls.parent_to_current_th + th_offset))
        return EB_TRUE;

    return EB_FALSE;
}

uint8_t is_child_to_current_deviation_small(SequenceControlSet  *scs_ptr,
                                            ModeDecisionContext *mdctxt, const BlockGeom *blk_geom,
                                            uint32_t blk_index, int64_t th_offset) {
    if (mdctxt->depth_refinement_ctrls.sub_to_current_th == MIN_SIGNED_VALUE)
        return EB_FALSE;
    const uint32_t ns_d1_offset = blk_geom->d1_depth_offset;

    (void)scs_ptr;
    assert(blk_geom->depth < 6);
    const uint32_t ns_depth_plus1_offset = ns_depth_offset[blk_geom->geom_idx][blk_geom->depth + 1];
    const uint32_t child_block_idx_1     = blk_index + ns_d1_offset;
    const uint32_t child_block_idx_2     = child_block_idx_1 + ns_depth_plus1_offset;
    const uint32_t child_block_idx_3     = child_block_idx_2 + ns_depth_plus1_offset;
    const uint32_t child_block_idx_4     = child_block_idx_3 + ns_depth_plus1_offset;

    uint64_t child_cost = 0;
    uint8_t  child_cnt  = 0;
    if (mdctxt->avail_blk_flag[child_block_idx_1]) {
        child_cost += mdctxt->md_local_blk_unit[child_block_idx_1].default_cost;
        child_cnt++;
    }
    if (mdctxt->avail_blk_flag[child_block_idx_2]) {
        child_cost += mdctxt->md_local_blk_unit[child_block_idx_2].default_cost;
        child_cnt++;
    }
    if (mdctxt->avail_blk_flag[child_block_idx_3]) {
        child_cost += mdctxt->md_local_blk_unit[child_block_idx_3].default_cost;
        child_cnt++;
    }
    if (mdctxt->avail_blk_flag[child_block_idx_4]) {
        child_cost += mdctxt->md_local_blk_unit[child_block_idx_4].default_cost;
        child_cnt++;
    }

    if (child_cnt) {
        child_cost = (child_cost / child_cnt) * 4;
        mdctxt->child_to_current_deviation =
            (int64_t)(((int64_t)MAX(child_cost, 1) -
                       (int64_t)MAX(mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost, 1)) *
                      100) /
            (int64_t)(MAX(mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost, 1));
    }

    if (mdctxt->child_to_current_deviation <=
        (mdctxt->depth_refinement_ctrls.sub_to_current_th + th_offset))
        return EB_TRUE;

    return EB_FALSE;
}
static void perform_pred_depth_refinement(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                          ModeDecisionContext *context_ptr, uint32_t sb_index) {
    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
    uint32_t   blk_index   = 0;
    uint8_t    use_cost_band_based_modulation =
        (!scs_ptr->vq_ctrls.stability_ctrls.depth_refinement ||
         (pcs_ptr->slice_type != I_SLICE &&
          pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[context_ptr->sb_index] <
              VQ_STABILITY_ME_VAR_TH));

    if (pcs_ptr->parent_pcs_ptr->disallow_nsq) {
        if (context_ptr->disallow_4x4) {
            memset(results_ptr->consider_block, 0, sizeof(uint8_t) * scs_ptr->max_block_cnt);
            memset(results_ptr->split_flag, 1, sizeof(uint8_t) * scs_ptr->max_block_cnt);
            memset(results_ptr->refined_split_flag, 1, sizeof(uint8_t) * scs_ptr->max_block_cnt);
        } else {
            while (blk_index < scs_ptr->max_block_cnt) {
                const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

                EbBool split_flag                      = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
                results_ptr->consider_block[blk_index] = 0;
                results_ptr->split_flag[blk_index]     = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
                results_ptr->refined_split_flag[blk_index] = blk_geom->sq_size > 4 ? EB_TRUE
                                                                                   : EB_FALSE;
                blk_index += split_flag ? blk_geom->d1_depth_offset : blk_geom->ns_depth_offset;
            }
        }
    } else {
        // Reset mdc_sb_array data to defaults; it will be updated based on the predicted blocks (stored in md_blk_arr_nsq)
        while (blk_index < scs_ptr->max_block_cnt) {
            const BlockGeom *blk_geom                  = get_blk_geom_mds(blk_index);
            results_ptr->consider_block[blk_index]     = 0;
            results_ptr->split_flag[blk_index]         = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
            results_ptr->refined_split_flag[blk_index] = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
            blk_index++;
        }
    }
    results_ptr->leaf_count = 0;
    blk_index               = 0;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom      = get_blk_geom_mds(blk_index);
        const unsigned   tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1
              : blk_geom->sq_size == 128                                       ? 17
              : blk_geom->sq_size > 8                                          ? 25
              : blk_geom->sq_size == 8                                         ? 5
                                                                               : 1;

        // if the parent square is inside inject this block
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1
            : (blk_geom->sq_size < 128)                         ? 1
                                                                : 0;

        // derive split_flag
        EbBool split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;

        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    // Add current pred depth block(s)
                    for (unsigned block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
                        results_ptr->consider_block[blk_index + block_1d_idx]     = 1;
                        results_ptr->refined_split_flag[blk_index + block_1d_idx] = EB_FALSE;
                    }

                    int8_t s_depth = context_ptr->depth_ctrls.s_depth;
                    int8_t e_depth = context_ptr->depth_ctrls.e_depth;

                    if (context_ptr->skip_pd0) {
                        SbParams *sb_params =
                            &pcs_ptr->parent_pcs_ptr->sb_params_array[context_ptr->sb_index];
                        if ((sb_params->width % 32 == 0) && (sb_params->height % 32 == 0)) {
                            s_depth = 0;
                            e_depth = 0;
                        }
                    }
                    // If multiple depths are selected, perform refinement
                    if (s_depth != 0 || e_depth != 0) {
                        // Check that the start and end depth are in allowed range, given other features
                        // which restrict allowable depths
                        if (context_ptr->disallow_4x4) {
                            e_depth = (blk_geom->sq_size == 8) ? 0
                                : (blk_geom->sq_size == 16)    ? MIN(1, e_depth)
                                : (blk_geom->sq_size == 32)    ? MIN(2, e_depth)
                                                               : e_depth;
                        }
                        if (context_ptr->depth_removal_ctrls.enabled) {
                            if (context_ptr->depth_removal_ctrls.disallow_below_64x64) {
                                e_depth = (blk_geom->sq_size <= 64) ? 0
                                    : (blk_geom->sq_size == 128)    ? MIN(1, e_depth)
                                                                    : e_depth;
                            } else if (context_ptr->depth_removal_ctrls.disallow_below_32x32) {
                                e_depth = (blk_geom->sq_size <= 32) ? 0
                                    : (blk_geom->sq_size == 64)     ? MIN(1, e_depth)
                                    : (blk_geom->sq_size == 128)    ? MIN(2, e_depth)
                                                                    : e_depth;
                            } else if (context_ptr->depth_removal_ctrls.disallow_below_16x16) {
                                e_depth = (blk_geom->sq_size <= 16) ? 0
                                    : (blk_geom->sq_size == 32)     ? MIN(1, e_depth)
                                    : (blk_geom->sq_size == 64)     ? MIN(2, e_depth)
                                    : (blk_geom->sq_size == 128)    ? MIN(3, e_depth)
                                                                    : e_depth;
                            }
                        }

                        uint8_t sq_size_idx = 7 - (uint8_t)svt_log2f((uint8_t)blk_geom->sq_size);
                        int64_t th_offset   = 0;

                        if (context_ptr->depth_refinement_ctrls.enabled &&
                            context_ptr->depth_refinement_ctrls.cost_band_based_modulation &&
                            use_cost_band_based_modulation && (s_depth != 0 || e_depth != 0)) {
                            update_pred_th_offset(
                                context_ptr, blk_geom, &s_depth, &e_depth, &th_offset);
                        }

                        // Add block indices of upper depth(s)
                        // Block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
                        uint8_t add_parent_depth                 = 1;
                        context_ptr->parent_to_current_deviation = MIN_SIGNED_VALUE;
                        if (context_ptr->depth_refinement_ctrls.enabled && s_depth == -1 &&
                            pcs_ptr->parent_pcs_ptr->sb_geom[sb_index]
                                .block_is_allowed[blk_index] &&
                            blk_geom->sq_size <
                                ((scs_ptr->seq_header.sb_size == BLOCK_128X128) ? 128 : 64)) {
                            add_parent_depth = is_parent_to_current_deviation_small(
                                context_ptr, blk_geom, th_offset);
                        }

                        // Add block indices of lower depth(s)
                        // Block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
                        uint8_t add_sub_depth                   = 1;
                        context_ptr->child_to_current_deviation = MIN_SIGNED_VALUE;
                        if (context_ptr->depth_refinement_ctrls.enabled && e_depth == 1 &&
                            pcs_ptr->parent_pcs_ptr->sb_geom[sb_index]
                                .block_is_allowed[blk_index] &&
                            blk_geom->sq_size > 4) {
                            add_sub_depth = is_child_to_current_deviation_small(
                                scs_ptr, context_ptr, blk_geom, blk_index, th_offset);
                        }

                        // Use a maximum of 2 depth per block (PRED+Parent or PRED+Sub)
                        if (context_ptr->depth_refinement_ctrls.enabled &&
                            context_ptr->depth_refinement_ctrls.up_to_2_depth) {
                            if ((s_depth == -1) && add_parent_depth && (e_depth == 1) &&
                                add_sub_depth) {
                                if (context_ptr->parent_to_current_deviation != MIN_SIGNED_VALUE &&
                                    context_ptr->child_to_current_deviation != MIN_SIGNED_VALUE) {
                                    if (context_ptr->parent_to_current_deviation <=
                                        context_ptr->child_to_current_deviation) {
                                        add_sub_depth = 0;
                                    } else {
                                        add_parent_depth = 0;
                                    }
                                }
                            }
                        }

                        if (s_depth != 0 && add_parent_depth)
                            set_parent_to_be_considered(results_ptr,
                                                        blk_index,
                                                        scs_ptr->seq_header.sb_size,
                                                        (int8_t)blk_geom->depth,
                                                        sq_size_idx,
                                                        pcs_ptr->parent_pcs_ptr->disallow_nsq,
                                                        s_depth);

                        if (e_depth != 0 && add_sub_depth)
                            set_child_to_be_considered(pcs_ptr,
                                                       context_ptr,
                                                       results_ptr,
                                                       blk_index,
                                                       sb_index,
                                                       scs_ptr->seq_header.sb_size,
                                                       (int8_t)blk_geom->depth,
                                                       sq_size_idx,
                                                       e_depth);
                    }
                }
            }
        }
        blk_index += split_flag ? blk_geom->d1_depth_offset : blk_geom->ns_depth_offset;
    }
}
// Initialize structures used to indicate which blocks will be tested at MD.
// MD data structures should be updated in init_block_data(), not here.
EbErrorType build_starting_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                            ModeDecisionContext *context_ptr, uint32_t sb_index) {
    memset(context_ptr->tested_blk_flag, 0, sizeof(uint8_t) * scs_ptr->max_block_cnt);
    memset(context_ptr->avail_blk_flag, EB_FALSE, sizeof(uint8_t) * scs_ptr->max_block_cnt);
    MdcSbData *results_ptr       = context_ptr->mdc_sb_array;
    results_ptr->leaf_count      = 0;
    uint32_t       blk_index     = 0;
    const uint16_t max_block_cnt = scs_ptr->max_block_cnt;
    const int32_t  min_sq_size   = (context_ptr->depth_removal_ctrls.enabled &&
                                 context_ptr->depth_removal_ctrls.disallow_below_64x64)
           ? 64
           : (context_ptr->depth_removal_ctrls.enabled &&
           context_ptr->depth_removal_ctrls.disallow_below_32x32)
           ? 32
           : (context_ptr->depth_removal_ctrls.enabled &&
           context_ptr->depth_removal_ctrls.disallow_below_16x16)
           ? 16
           : context_ptr->disallow_4x4 ? 8
                                       : 4;

    // Loop over all blocks to initialize data for partitions to be tested
    while (blk_index < max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
        // SQ/NSQ block(s) filter based on the SQ size
        uint8_t is_block_tagged = (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE) ||
                (blk_geom->sq_size < min_sq_size)
            ? 0
            : 1;
        if (context_ptr->skip_pd0)
            is_block_tagged = !context_ptr->depth_removal_ctrls.disallow_below_64x64 &&
                    context_ptr->depth_removal_ctrls.disallow_below_32x32 &&
                    (blk_geom->sq_size != 32)
                ? 0
                : is_block_tagged;

        // SQ/NSQ block(s) filter based on the block validity
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] &&
            is_block_tagged) {
            const uint32_t tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq
                ? 1
                : get_tot_1d_blks(
                      pcs_ptr->parent_pcs_ptr, blk_geom->sq_size, context_ptr->md_disallow_nsq);

            for (uint32_t idx = blk_index; idx < (tot_d1_blocks + blk_index); ++idx) {
                if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[idx]) {
                    //  MD palette info buffer
                    if (pcs_ptr->parent_pcs_ptr->palette_level) {
                        if (context_ptr->md_blk_arr_nsq[idx].palette_mem == 0) {
                            rtime_alloc_palette_info(&context_ptr->md_blk_arr_nsq[idx]);
                            context_ptr->md_blk_arr_nsq[idx].palette_mem = 1;
                        }
                    }

                    context_ptr->md_blk_arr_nsq[idx].palette_size[0]              = 0;
                    context_ptr->md_blk_arr_nsq[idx].palette_size[1]              = 0;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = idx;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks =
                        tot_d1_blocks;
                    results_ptr->split_flag[results_ptr->leaf_count++] = (blk_geom->sq_size >
                                                                          min_sq_size)
                        ? EB_TRUE
                        : EB_FALSE;
                }
            }
            blk_index += blk_geom->d1_depth_offset;
        } else {
            if (context_ptr->skip_pd0)
                context_ptr->md_blk_arr_nsq[blk_index].part = (blk_geom->sq_size > min_sq_size)
                    ? PARTITION_SPLIT
                    : PARTITION_NONE;
            blk_index += (blk_geom->sq_size > min_sq_size) ? blk_geom->d1_depth_offset
                                                           : blk_geom->ns_depth_offset;
        }
    }

    return EB_ErrorNone;
}
void recode_loop_update_q(PictureParentControlSet *ppcs_ptr, int *const loop, int *const q,
                          int *const q_low, int *const q_high, const int top_index,
                          const int bottom_index, int *const undershoot_seen,
                          int *const overshoot_seen, int *const low_cr_seen, const int loop_count);
void sb_qp_derivation_tpl_la(PictureControlSet *pcs_ptr);
void mode_decision_configuration_init_qp_update(PictureControlSet *pcs_ptr);
void init_enc_dec_segement(PictureParentControlSet *parentpicture_control_set_ptr);

static void recode_loop_decision_maker(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr,
                                       EbBool *do_recode) {
    PictureParentControlSet *ppcs_ptr           = pcs_ptr->parent_pcs_ptr;
    EncodeContext *const     encode_context_ptr = ppcs_ptr->scs_ptr->encode_context_ptr;
    RATE_CONTROL *const      rc                 = &(encode_context_ptr->rc);
    int32_t                  loop               = 0;
    FrameHeader             *frm_hdr            = &ppcs_ptr->frm_hdr;
    int32_t                  q                  = frm_hdr->quantization_params.base_q_idx;
    if (ppcs_ptr->loop_count == 0) {
        ppcs_ptr->q_low  = ppcs_ptr->bottom_index;
        ppcs_ptr->q_high = ppcs_ptr->top_index;
    }

    // Update q and decide whether to do a recode loop
    recode_loop_update_q(ppcs_ptr,
                         &loop,
                         &q,
                         &ppcs_ptr->q_low,
                         &ppcs_ptr->q_high,
                         ppcs_ptr->top_index,
                         ppcs_ptr->bottom_index,
                         &ppcs_ptr->undershoot_seen,
                         &ppcs_ptr->overshoot_seen,
                         &ppcs_ptr->low_cr_seen,
                         ppcs_ptr->loop_count);

    // Special case for overlay frame.
    if (loop && ppcs_ptr->is_src_frame_alt_ref &&
        ppcs_ptr->projected_frame_size < rc->max_frame_bandwidth) {
        loop = 0;
    }
    *do_recode = loop == 1;

    if (*do_recode) {
        ppcs_ptr->loop_count++;

        frm_hdr->quantization_params.base_q_idx = (uint8_t)CLIP3(
            (int32_t)quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed],
            (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed],
            q);

        ppcs_ptr->picture_qp = (uint8_t)CLIP3((int32_t)scs_ptr->static_config.min_qp_allowed,
                                              (int32_t)scs_ptr->static_config.max_qp_allowed,
                                              (frm_hdr->quantization_params.base_q_idx + 2) >> 2);
        pcs_ptr->picture_qp  = ppcs_ptr->picture_qp;

        // 2pass QPM with tpl_la
        if (scs_ptr->static_config.enable_adaptive_quantization == 2 &&
            ppcs_ptr->tpl_ctrls.enable && ppcs_ptr->r0 != 0)
            sb_qp_derivation_tpl_la(pcs_ptr);
        else {
            ppcs_ptr->frm_hdr.delta_q_params.delta_q_present = 0;
            for (int sb_addr = 0; sb_addr < pcs_ptr->sb_total_count_pix; ++sb_addr) {
                SuperBlock *sb_ptr = pcs_ptr->sb_ptr_array[sb_addr];
                sb_ptr->qindex     = quantizer_to_qindex[pcs_ptr->picture_qp];
            }
        }
    } else {
        ppcs_ptr->loop_count = 0;
    }
}

/* for debug/documentation purposes: list all features assumed off for light pd1*/
void exaustive_light_pd1_features(ModeDecisionContext *md_ctx, PictureParentControlSet *ppcs,
                                  uint8_t use_light_pd1, uint8_t debug_lpd1_features) {
    if (debug_lpd1_features) {
        uint8_t light_pd1;

        // Use light-PD1 path if the assumed features are off
        if (md_ctx->obmc_ctrls.enabled == 0 && md_ctx->md_allow_intrabc == 0 &&
            md_ctx->hbd_mode_decision == 0 && md_ctx->ifs_ctrls.level == IFS_OFF &&
            ppcs->frm_hdr.allow_warped_motion == 0 && md_ctx->inter_intra_comp_ctrls.enabled == 0 &&
            md_ctx->rate_est_ctrls.update_skip_ctx_dc_sign_ctx == 0 &&
            md_ctx->spatial_sse_ctrls.spatial_sse_full_loop_level == 0 &&
            md_ctx->md_sq_me_ctrls.enabled == 0 && md_ctx->md_pme_ctrls.enabled == 0 &&
            md_ctx->txt_ctrls.enabled == 0 && md_ctx->mds0_ctrls.mds0_dist_type != MDS0_SSD &&
            md_ctx->unipred3x3_injection == 0 && md_ctx->bipred3x3_injection == 0 &&
            md_ctx->inter_compound_mode == 0 && md_ctx->md_pic_obmc_level == 0 &&
            md_ctx->md_filter_intra_level == 0 && md_ctx->new_nearest_near_comb_injection == 0 &&
            md_ctx->md_palette_level == 0 && md_ctx->cand_reduction_ctrls.merge_inter_classes &&
            ppcs->gm_ctrls.enabled == 0 &&
            // If TXS enabled at picture level, there are necessary context updates that must be added to LPD1
            ppcs->frm_hdr.tx_mode != TX_MODE_SELECT && md_ctx->txs_ctrls.enabled == 0 &&
            md_ctx->pred_depth_only && md_ctx->md_disallow_nsq == EB_TRUE &&
            md_ctx->disallow_4x4 == EB_TRUE && ppcs->scs_ptr->super_block_size == 64 &&
            ppcs->ref_list0_count_try == 1 && ppcs->ref_list1_count_try == 1 &&
            md_ctx->cfl_ctrls.enabled == 0 && md_ctx->uv_ctrls.nd_uv_serach_mode == 0 &&
            md_ctx->uv_ctrls.uv_mode == CHROMA_MODE_1) {
            light_pd1 = 1;
        } else {
            light_pd1 = 0;
        }

        assert_err(light_pd1 == use_light_pd1,
                   "Warning: light PD1 feature assumption is broken \n");
    }
}
/* Light-PD1 classifier used when cost/coeff info is available.  If PD0 is skipped, or the trasnsform is
not performed, a separate detector (lpd1_detector_skip_pd0) is used. */
void lpd1_detector_post_pd0(PictureControlSet *pcs, ModeDecisionContext *md_ctx) {
    for (int pd1_lvl = LPD1_LEVELS - 1; pd1_lvl > REGULAR_PD1; pd1_lvl--) {
        if (md_ctx->lpd1_ctrls.pd1_level == pd1_lvl) {
            if (md_ctx->lpd1_ctrls.use_lpd1_detector[pd1_lvl]) {
                // Use info from ref frames (if available)
                if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl] && pcs->slice_type != I_SLICE) {
                    EbReferenceObject *ref_obj_l0 =
                        (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                    uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
                    if (pcs->slice_type == B_SLICE) {
                        EbReferenceObject *ref_obj_l1 =
                            (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                        l1_was_intra = ref_obj_l1->sb_intra[md_ctx->sb_index];
                    }
                    if (l0_was_intra && l1_was_intra) {
                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        continue;
                    } else if (l0_was_intra || l1_was_intra) {
                        md_ctx->lpd1_ctrls.cost_th_dist[pd1_lvl] >>= 7;
                        md_ctx->lpd1_ctrls.coeff_th[pd1_lvl] >>= 6;
                        md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] >>= 4;
                    }
                }

                /* Use the cost and coeffs of the 64x64 block to avoid looping over all tested blocks to find
                the selected partitioning. */
                const uint64_t pd0_cost = md_ctx->md_local_blk_unit[0].cost;
                // If block was not tested in PD0, won't have coeff info, so set to max and base detection on cost only (which is set
                // even if 64x64 block is not tested)
                const uint32_t nz_coeffs = md_ctx->avail_blk_flag[0] == EB_TRUE
                    ? md_ctx->md_local_blk_unit[0].count_non_zero_coeffs
                    : (uint32_t)~0;

                const uint32_t lambda =
                    md_ctx->full_sb_lambda_md[EB_8_BIT_MD]; // light-PD1 assumes 8-bit MD
                const uint32_t rate = md_ctx->lpd1_ctrls.coeff_th[pd1_lvl];
                const uint32_t dist = md_ctx->lpd1_ctrls.cost_th_dist[pd1_lvl];
                /* dist << 14 is equivalent to 64 * 64 * 4 * dist (64 * 64 so the distortion is the per-pixel SSD) and 4 because
                the distortion of the 64x64 block is shifted by 2 (same as multiplying by 4) in perform_tx_light_pd0. */
                const uint64_t low_th   = RDCOST(lambda, 6000 + rate * 500, (uint64_t)dist << 14);
                const uint16_t coeff_th = md_ctx->lpd1_ctrls.coeff_th[pd1_lvl];

                // If the PD0 cost is very high and the number of non-zero coeffs is high, the block is difficult, so should use regular PD1
                if (pd0_cost > low_th && nz_coeffs > coeff_th) {
                    md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                }

                // If the best PD0 mode was INTER, check the MV length
                if (md_ctx->avail_blk_flag[0] == EB_TRUE &&
                    md_ctx->md_blk_arr_nsq[0].prediction_mode_flag == INTER_MODE &&
                    md_ctx->lpd1_ctrls.max_mv_length[pd1_lvl] != (uint16_t)~0) {
                    PredictionUnit *pu_ptr        = md_ctx->md_blk_arr_nsq[0].prediction_unit_array;
                    const uint16_t  max_mv_length = md_ctx->lpd1_ctrls.max_mv_length[pd1_lvl];

                    if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0) {
                        if (pu_ptr->mv[REF_LIST_0].x > max_mv_length ||
                            pu_ptr->mv[REF_LIST_0].y > max_mv_length)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                    } else if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1) {
                        if (pu_ptr->mv[REF_LIST_1].x > max_mv_length ||
                            pu_ptr->mv[REF_LIST_1].y > max_mv_length)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                    } else //if (pu_ptr->inter_pred_direction_index == BI_PRED)
                    {
                        assert(pu_ptr->inter_pred_direction_index == BI_PRED);
                        if (pu_ptr->mv[REF_LIST_0].x > max_mv_length ||
                            pu_ptr->mv[REF_LIST_0].y > max_mv_length)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        if (pu_ptr->mv[REF_LIST_1].x > max_mv_length ||
                            pu_ptr->mv[REF_LIST_1].y > max_mv_length)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                    }
                }

                if (pcs->slice_type != I_SLICE) {
                    /* me_8x8_cost_variance_th is shifted by 5 then mulitplied by the pic QP (max 63).  Therefore, the TH must be less than
                       (((uint32_t)~0) >> 1) to avoid overflow issues from the multiplication. */
                    if (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] <
                            (((uint32_t)~0) >> 1) &&
                        pcs->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] >
                            (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] >> 5) *
                                pcs->picture_qp)
                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                }
            }
        }
    }
}

/* Light-PD1 classifier used when cost/coeff info is unavailable.  If PD0 is skipped, or the trasnsform is
not performed, this detector is used (else lpd1_detector_post_pd0() is used). */
void lpd1_detector_skip_pd0(PictureControlSet *pcs, ModeDecisionContext *md_ctx,
                            uint32_t pic_width_in_sb) {
    const uint16_t left_sb_index = md_ctx->sb_index - 1;
    const uint16_t top_sb_index  = md_ctx->sb_index - (uint16_t)pic_width_in_sb;

    for (int pd1_lvl = LPD1_LEVELS - 1; pd1_lvl > REGULAR_PD1; pd1_lvl--) {
        if (md_ctx->lpd1_ctrls.pd1_level == pd1_lvl) {
            if (md_ctx->lpd1_ctrls.use_lpd1_detector[pd1_lvl]) {
                // Use info from ref. frames (if available)
                if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl] && pcs->slice_type != I_SLICE) {
                    EbReferenceObject *ref_obj_l0 =
                        (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                    uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
                    uint8_t l0_was_skip = ref_obj_l0->sb_skip[md_ctx->sb_index], l1_was_skip = 1;
                    uint32_t l0_me_64x64_dist   = ref_obj_l0->slice_type != I_SLICE
                          ? ref_obj_l0->sb_me_64x64_dist[md_ctx->sb_index]
                          : 0,
                             l1_me_64x64_dist   = 0;
                    uint32_t l0_me_8x8_cost_var = ref_obj_l0->slice_type != I_SLICE
                        ? ref_obj_l0->sb_me_8x8_cost_var[md_ctx->sb_index]
                        : 0,
                             l1_me_8x8_cost_var = 0;
                    if (pcs->slice_type == B_SLICE) {
                        EbReferenceObject *ref_obj_l1 =
                            (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                        l1_was_intra       = ref_obj_l1->sb_intra[md_ctx->sb_index];
                        l1_was_skip        = ref_obj_l1->sb_skip[md_ctx->sb_index];
                        l1_me_64x64_dist   = ref_obj_l1->slice_type != I_SLICE
                              ? ref_obj_l1->sb_me_64x64_dist[md_ctx->sb_index]
                              : 0;
                        l1_me_8x8_cost_var = ref_obj_l1->slice_type != I_SLICE
                            ? ref_obj_l1->sb_me_8x8_cost_var[md_ctx->sb_index]
                            : 0;
                    }

                    // Keep a complexity score for the SB, based on available information.
                    // If the score is high, then reduce the lpd1_level to be used
                    int16_t score = 0;

                    if (l0_was_intra)
                        score += 5;
                    if (l1_was_intra)
                        score += 5;

                    if (!l0_was_skip)
                        score += 5;
                    if (!l1_was_skip)
                        score += 5;

                    if (pcs->slice_type == B_SLICE) {
                        if (pcs->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] >
                                (l0_me_64x64_dist * 3) ||
                            pcs->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] >
                                (l1_me_64x64_dist * 3))
                            score += 10;

                        if (pcs->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] >
                                (l0_me_8x8_cost_var * 3) ||
                            pcs->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] >
                                (l1_me_8x8_cost_var * 3))
                            score += 10;
                    } else {
                        score += 20;
                    }

                    if (score >= 20) {
                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        continue;
                    }
                }

                // I_SLICE doesn't have ME info
                if (pcs->slice_type != I_SLICE) {
                    // If the SB origin of one dimension is zero, then this SB is the first block in a row/column, so won't have neighbours
                    if (md_ctx->sb_origin_x == 0 || md_ctx->sb_origin_y == 0) {
                        if (pcs->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] >
                            md_ctx->lpd1_ctrls.skip_pd0_edge_dist_th[pd1_lvl])
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        /* me_8x8_cost_variance_th is shifted by 5 then mulitplied by the pic QP (max 63).  Therefore, the TH must be less than
                           (((uint32_t)~0) >> 1) to avoid overflow issues from the multiplication. */
                        else if (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] <
                                     (((uint32_t)~0) >> 1) &&
                                 pcs->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] >
                                     (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] >> 5) *
                                         pcs->picture_qp)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                    } else {
                        if (md_ctx->lpd1_ctrls.skip_pd0_me_shift[pd1_lvl] != (uint16_t)~0 &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] >
                                ((pcs->parent_pcs_ptr->me_64x64_distortion[left_sb_index] +
                                  pcs->parent_pcs_ptr->me_64x64_distortion[top_sb_index])
                                 << md_ctx->lpd1_ctrls.skip_pd0_me_shift[pd1_lvl]))
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        else if (md_ctx->lpd1_ctrls.skip_pd0_me_shift[pd1_lvl] != (uint16_t)~0 &&
                                 pcs->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] >
                                     ((pcs->parent_pcs_ptr->me_8x8_cost_variance[left_sb_index] +
                                       pcs->parent_pcs_ptr->me_8x8_cost_variance[top_sb_index])
                                      << md_ctx->lpd1_ctrls.skip_pd0_me_shift[pd1_lvl])) {
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        } else if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl]) {
                            // Use info from neighbouring SBs
                            if (pcs->sb_intra[left_sb_index] && pcs->sb_intra[top_sb_index]) {
                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                            } else if (!pcs->sb_skip[left_sb_index] &&
                                       !pcs->sb_skip[top_sb_index] &&
                                       (pcs->sb_intra[left_sb_index] ||
                                        pcs->sb_intra[top_sb_index])) {
                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

/*
* Check whether vlpd0 is safe or not
*/
uint8_t is_vlpd0_safe(PictureControlSet *pcs_ptr, ModeDecisionContext *md_ctx) {
    uint8_t is_vlpd0_safe = EB_TRUE;

    EbReferenceObject *ref_obj_l0 =
        (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
    uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
    if (pcs_ptr->slice_type == B_SLICE) {
        EbReferenceObject *ref_obj_l1 =
            (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
        l1_was_intra = ref_obj_l1->sb_intra[md_ctx->sb_index];
    }
    if (l0_was_intra || l1_was_intra) {
        return EB_FALSE;
    }

    uint32_t me_8x8_cost_variance_th = 250000;
    if (pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] >
        (me_8x8_cost_variance_th >> 5) * pcs_ptr->picture_qp)
        return EB_FALSE;

    return is_vlpd0_safe;
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
    EbThreadContext *thread_context_ptr = (EbThreadContext *)input_ptr;
    EncDecContext   *context_ptr        = (EncDecContext *)thread_context_ptr->priv;

    // Input
    EbObjectWrapper *enc_dec_tasks_wrapper_ptr;

    // Output
    EbObjectWrapper *enc_dec_results_wrapper_ptr;
    EncDecResults   *enc_dec_results_ptr;
    // SB Loop variables
    SuperBlock *sb_ptr;
    uint16_t    sb_index;
    uint32_t    x_sb_index;
    uint32_t    y_sb_index;
    uint32_t    sb_origin_x;
    uint32_t    sb_origin_y;
    MdcSbData  *mdc_ptr;

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

        EncDecTasks       *enc_dec_tasks_ptr = (EncDecTasks *)enc_dec_tasks_wrapper_ptr->object_ptr;
        PictureControlSet *pcs_ptr           = (PictureControlSet *)
                                         enc_dec_tasks_ptr->pcs_wrapper_ptr->object_ptr;
        SequenceControlSet  *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
        ModeDecisionContext *md_ctx  = context_ptr->md_context;
        struct PictureParentControlSet *ppcs = pcs_ptr->parent_pcs_ptr;
        md_ctx->encoder_bit_depth            = (uint8_t)scs_ptr->static_config.encoder_bit_depth;
        md_ctx->corrupted_mv_check           = (pcs_ptr->parent_pcs_ptr->aligned_width >=
                                      (1 << (MV_IN_USE_BITS - 3))) ||
            (pcs_ptr->parent_pcs_ptr->aligned_height >= (1 << (MV_IN_USE_BITS - 3)));
        context_ptr->tile_group_index = enc_dec_tasks_ptr->tile_group_index;
        context_ptr->coded_sb_count   = 0;
        segments_ptr = pcs_ptr->enc_dec_segment_ctrl[context_ptr->tile_group_index];
        // SB Constants
        uint8_t sb_sz            = (uint8_t)scs_ptr->sb_size_pix;
        uint8_t sb_size_log2     = (uint8_t)svt_log2f(sb_sz);
        context_ptr->sb_sz       = sb_sz;
        uint32_t pic_width_in_sb = (pcs_ptr->parent_pcs_ptr->aligned_width + sb_sz - 1) >>
            sb_size_log2;
        uint16_t tile_group_width_in_sb = pcs_ptr->parent_pcs_ptr
                                              ->tile_group_info[context_ptr->tile_group_index]
                                              .tile_group_width_in_sb;
        context_ptr->tot_intra_coded_area = 0;
        context_ptr->tot_skip_coded_area  = 0;
        // Bypass encdec for the first pass
        if (scs_ptr->static_config.pass == ENC_FIRST_PASS ||
            (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag &&
             scs_ptr->rc_stat_gen_pass_mode && !pcs_ptr->parent_pcs_ptr->first_frame_in_minigop)) {
            svt_release_object(pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr);
            pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr = (EbObjectWrapper *)NULL;
            pcs_ptr->parent_pcs_ptr->pa_me_data          = NULL;
            // Get Empty EncDec Results
            svt_get_empty_object(context_ptr->enc_dec_output_fifo_ptr,
                                 &enc_dec_results_wrapper_ptr);
            enc_dec_results_ptr = (EncDecResults *)enc_dec_results_wrapper_ptr->object_ptr;
            enc_dec_results_ptr->pcs_wrapper_ptr = enc_dec_tasks_ptr->pcs_wrapper_ptr;

            // Post EncDec Results
            svt_post_full_object(enc_dec_results_wrapper_ptr);
        } else {
            if (enc_dec_tasks_ptr->input_type == ENCDEC_TASKS_SUPERRES_INPUT) {
                // do as dorecode do
                pcs_ptr->enc_dec_coded_sb_count = 0;
                // re-init mode decision configuration for qp update for re-encode frame
                mode_decision_configuration_init_qp_update(pcs_ptr);
                // init segment for re-encode frame
                init_enc_dec_segement(pcs_ptr->parent_pcs_ptr);

                // post tile based encdec task
                EbObjectWrapper *enc_dec_re_encode_tasks_wrapper_ptr;
                uint16_t         tg_count = pcs_ptr->parent_pcs_ptr->tile_group_cols *
                    pcs_ptr->parent_pcs_ptr->tile_group_rows;
                for (uint16_t tile_group_idx = 0; tile_group_idx < tg_count; tile_group_idx++) {
                    svt_get_empty_object(context_ptr->enc_dec_feedback_fifo_ptr,
                                         &enc_dec_re_encode_tasks_wrapper_ptr);

                    EncDecTasks *enc_dec_re_encode_tasks_ptr =
                        (EncDecTasks *)enc_dec_re_encode_tasks_wrapper_ptr->object_ptr;
                    enc_dec_re_encode_tasks_ptr->pcs_wrapper_ptr =
                        enc_dec_tasks_ptr->pcs_wrapper_ptr;
                    enc_dec_re_encode_tasks_ptr->input_type       = ENCDEC_TASKS_MDC_INPUT;
                    enc_dec_re_encode_tasks_ptr->tile_group_index = tile_group_idx;

                    // Post the Full Results Object
                    svt_post_full_object(enc_dec_re_encode_tasks_wrapper_ptr);
                }

                svt_release_object(enc_dec_tasks_wrapper_ptr);
                continue;
            }

            if (pcs_ptr->cdf_ctrl.enabled) {
                if (!pcs_ptr->cdf_ctrl.update_mv)
                    copy_mv_rate(pcs_ptr, &context_ptr->md_context->rate_est_table);
                if (!pcs_ptr->cdf_ctrl.update_se)

                    av1_estimate_syntax_rate(
                        &context_ptr->md_context->rate_est_table,
                        pcs_ptr->slice_type == I_SLICE ? EB_TRUE : EB_FALSE,
                        pcs_ptr->pic_filter_intra_level,
                        pcs_ptr->parent_pcs_ptr->frm_hdr.allow_screen_content_tools,
                        scs_ptr->seq_header.enable_restoration,
                        pcs_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc,
                        pcs_ptr->parent_pcs_ptr->partition_contexts,
                        &pcs_ptr->md_frame_context);
                if (!pcs_ptr->cdf_ctrl.update_coef)
                    av1_estimate_coefficients_rate(&context_ptr->md_context->rate_est_table,
                                                   &pcs_ptr->md_frame_context);
            }
            // Segment-loop
            while (assign_enc_dec_segments(segments_ptr,
                                           &segment_index,
                                           enc_dec_tasks_ptr,
                                           context_ptr->enc_dec_feedback_fifo_ptr) == EB_TRUE) {
                x_sb_start_index = segments_ptr->x_start_array[segment_index];
                y_sb_start_index = segments_ptr->y_start_array[segment_index];
                sb_start_index   = y_sb_start_index * tile_group_width_in_sb + x_sb_start_index;
                sb_segment_count = segments_ptr->valid_sb_count_array[segment_index];

                segment_row_index  = segment_index / segments_ptr->segment_band_count;
                segment_band_index = segment_index -
                    segment_row_index * segments_ptr->segment_band_count;
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
                for (y_sb_index = y_sb_start_index, sb_segment_index = sb_start_index;
                     sb_segment_index < sb_start_index + sb_segment_count;
                     ++y_sb_index) {
                    for (x_sb_index = x_sb_start_index; x_sb_index < tile_group_width_in_sb &&
                         (x_sb_index + y_sb_index < segment_band_size) &&
                         sb_segment_index < sb_start_index + sb_segment_count;
                         ++x_sb_index, ++sb_segment_index) {
                        uint16_t tile_group_y_sb_start =
                            pcs_ptr->parent_pcs_ptr->tile_group_info[context_ptr->tile_group_index]
                                .tile_group_sb_start_y;
                        uint16_t tile_group_x_sb_start =
                            pcs_ptr->parent_pcs_ptr->tile_group_info[context_ptr->tile_group_index]
                                .tile_group_sb_start_x;
                        sb_index = context_ptr->md_context->sb_index =
                            (uint16_t)((y_sb_index + tile_group_y_sb_start) * pic_width_in_sb +
                                       x_sb_index + tile_group_x_sb_start);
                        sb_ptr = context_ptr->md_context->sb_ptr = pcs_ptr->sb_ptr_array[sb_index];
                        sb_origin_x = (x_sb_index + tile_group_x_sb_start) << sb_size_log2;
                        sb_origin_y = (y_sb_index + tile_group_y_sb_start) << sb_size_log2;
                        //printf("[%ld]:ED sb index %d, (%d, %d), encoded total sb count %d, ctx coded sb count %d\n",
                        //        pcs_ptr->picture_number,
                        //        sb_index, sb_origin_x, sb_origin_y,
                        //        pcs_ptr->enc_dec_coded_sb_count,
                        //        context_ptr->coded_sb_count);
                        context_ptr->tile_index              = sb_ptr->tile_info.tile_rs_index;
                        context_ptr->md_context->tile_index  = sb_ptr->tile_info.tile_rs_index;
                        context_ptr->md_context->sb_origin_x = sb_origin_x;
                        context_ptr->md_context->sb_origin_y = sb_origin_y;
                        mdc_ptr               = context_ptr->md_context->mdc_sb_array;
                        context_ptr->sb_index = sb_index;
                        if (pcs_ptr->cdf_ctrl.enabled) {
                            if (scs_ptr->seq_header.pic_based_rate_est &&
                                scs_ptr->enc_dec_segment_row_count_array
                                        [pcs_ptr->temporal_layer_index] == 1 &&
                                scs_ptr->enc_dec_segment_col_count_array
                                        [pcs_ptr->temporal_layer_index] == 1) {
                                if (sb_index == 0)
                                    pcs_ptr->ec_ctx_array[sb_index] = pcs_ptr->md_frame_context;
                                else
                                    pcs_ptr->ec_ctx_array[sb_index] =
                                        pcs_ptr->ec_ctx_array[sb_index - 1];
                            } else {
                                // Use the latest available CDF for the current SB
                                // Use the weighted average of left (3x) and top right (1x) if available.
                                int8_t top_right_available = ((int32_t)(sb_origin_y >>
                                                                        MI_SIZE_LOG2) >
                                                              sb_ptr->tile_info.mi_row_start) &&
                                    ((int32_t)((sb_origin_x + (1 << sb_size_log2)) >>
                                               MI_SIZE_LOG2) < sb_ptr->tile_info.mi_col_end);

                                int8_t left_available = ((int32_t)(sb_origin_x >> MI_SIZE_LOG2) >
                                                         sb_ptr->tile_info.mi_col_start);

                                if (!left_available && !top_right_available)
                                    pcs_ptr->ec_ctx_array[sb_index] = pcs_ptr->md_frame_context;
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
                            if (pcs_ptr->cdf_ctrl.update_se)
                                av1_estimate_syntax_rate(
                                    &context_ptr->md_context->rate_est_table,
                                    pcs_ptr->slice_type == I_SLICE,
                                    pcs_ptr->pic_filter_intra_level,
                                    pcs_ptr->parent_pcs_ptr->frm_hdr.allow_screen_content_tools,
                                    scs_ptr->seq_header.enable_restoration,
                                    pcs_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc,
                                    pcs_ptr->parent_pcs_ptr->partition_contexts,
                                    &pcs_ptr->ec_ctx_array[sb_index]);
                            // Initial Rate Estimation of the Motion vectors
                            if (pcs_ptr->cdf_ctrl.update_mv)
                                av1_estimate_mv_rate(pcs_ptr,
                                                     &context_ptr->md_context->rate_est_table,
                                                     &pcs_ptr->ec_ctx_array[sb_index]);

                            if (pcs_ptr->cdf_ctrl.update_coef)
                                av1_estimate_coefficients_rate(
                                    &context_ptr->md_context->rate_est_table,
                                    &pcs_ptr->ec_ctx_array[sb_index]);
                            context_ptr->md_context->md_rate_estimation_ptr =
                                &context_ptr->md_context->rate_est_table;
                        }
                        // Configure the SB
                        mode_decision_configure_sb(
                            context_ptr->md_context, pcs_ptr, (uint8_t)sb_ptr->qindex);
                        // signals set once per SB (i.e. not per PD)
                        signal_derivation_enc_dec_kernel_common(
                            scs_ptr, pcs_ptr, context_ptr->md_context);

                        if (pcs_ptr->parent_pcs_ptr->palette_level)
                            // Status of palette info alloc
                            for (int i = 0; i < scs_ptr->max_block_cnt; ++i)
                                context_ptr->md_context->md_blk_arr_nsq[i].palette_mem = 0;

                        // Initialize is_subres_safe
                        context_ptr->md_context->is_subres_safe = (uint8_t)~0;
                        // Signal initialized here; if needed, will be set in md_encode_block before MDS3
                        md_ctx->need_hbd_comp_mds3 = 0;
                        uint8_t skip_pd_pass_0 =
                            (scs_ptr->super_block_size == 64 &&
                             context_ptr->md_context->depth_removal_ctrls.disallow_below_64x64)
                            ? 1
                            : 0;
                        if (context_ptr->md_context->skip_pd0)
                            if (context_ptr->md_context->depth_removal_ctrls.disallow_below_32x32)
                                skip_pd_pass_0 = 1;
                        if (context_ptr->md_context->pd0_level == VERY_LIGHT_PD0) {
                            // Use the next conservative level if not safe to use VLPD0
                            if (!is_vlpd0_safe(pcs_ptr, md_ctx))
                                context_ptr->md_context->pd0_level = VERY_LIGHT_PD0 - 1;
                        }
                        // PD0 is only skipped if there is a single depth to test
                        if (skip_pd_pass_0)
                            md_ctx->pred_depth_only = 1;
                        // Multi-Pass PD
                        if (!skip_pd_pass_0 &&
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_ON) {
                            // [PD_PASS_0]
                            // Input : mdc_blk_ptr built @ mdc process (up to 4421)
                            // Output: md_blk_arr_nsq reduced set of block(s)
                            context_ptr->md_context->pd_pass = PD_PASS_0;
                            // skip_intra much be TRUE for non-I_SLICE pictures to use light_pd0 path
                            if (context_ptr->md_context->pd0_level != REGULAR_PD0) {
                                // [PD_PASS_0] Signal(s) derivation
                                signal_derivation_enc_dec_kernel_oq_light_pd0(
                                    scs_ptr, pcs_ptr, context_ptr->md_context);

                                // Save a clean copy of the neighbor arrays
                                if (!context_ptr->md_context->skip_intra)
                                    copy_neighbour_arrays_light_pd0(
                                        pcs_ptr,
                                        context_ptr->md_context,
                                        MD_NEIGHBOR_ARRAY_INDEX,
                                        MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                        0,
                                        sb_origin_x,
                                        sb_origin_y);

                                // Build the t=0 cand_block_array
                                build_starting_cand_block_array(
                                    scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
                                mode_decision_sb_light_pd0(scs_ptr,
                                                           pcs_ptr,
                                                           mdc_ptr,
                                                           sb_ptr,
                                                           sb_origin_x,
                                                           sb_origin_y,
                                                           sb_index,
                                                           context_ptr->md_context);
                                // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                                // Reset neighnor information to current SB @ position (0,0)
                                if (!context_ptr->md_context->skip_intra)
                                    copy_neighbour_arrays_light_pd0(
                                        pcs_ptr,
                                        context_ptr->md_context,
                                        MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                        MD_NEIGHBOR_ARRAY_INDEX,
                                        0,
                                        sb_origin_x,
                                        sb_origin_y);
                            } else {
                                // [PD_PASS_0] Signal(s) derivation
                                signal_derivation_enc_dec_kernel_oq(
                                    scs_ptr, pcs_ptr, context_ptr->md_context);

                                // Save a clean copy of the neighbor arrays
                                copy_neighbour_arrays(pcs_ptr,
                                                      context_ptr->md_context,
                                                      MD_NEIGHBOR_ARRAY_INDEX,
                                                      MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                                      0,
                                                      sb_origin_x,
                                                      sb_origin_y);

                                // Build the t=0 cand_block_array
                                build_starting_cand_block_array(
                                    scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
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
                                // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                                // Reset neighnor information to current SB @ position (0,0)
                                copy_neighbour_arrays(pcs_ptr,
                                                      context_ptr->md_context,
                                                      MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                                      MD_NEIGHBOR_ARRAY_INDEX,
                                                      0,
                                                      sb_origin_x,
                                                      sb_origin_y);
                            }
                            // This classifier is used for only pd0_level 0 and pd0_level 1
                            // where the count_non_zero_coeffs is derived @ PD0
                            if (context_ptr->md_context->pd0_level != VERY_LIGHT_PD0)
                                lpd1_detector_post_pd0(pcs_ptr, md_ctx);
                            // Force pred depth only for modes where that is not the default
                            if (md_ctx->lpd1_ctrls.pd1_level > REGULAR_PD1) {
                                set_depth_ctrls(md_ctx, 0);
                                md_ctx->pred_depth_only = 1;
                            }
                            // Perform Pred_0 depth refinement - add depth(s) to be considered in the next stage(s)
                            perform_pred_depth_refinement(
                                scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
                        }
                        // [PD_PASS_1] Signal(s) derivation
                        context_ptr->md_context->pd_pass = PD_PASS_1;
                        // This classifier is used for the case PD0 is bypassed and for pd0_level 2
                        // where the count_non_zero_coeffs is not derived @ PD0
                        if (skip_pd_pass_0 ||
                            context_ptr->md_context->pd0_level == VERY_LIGHT_PD0) {
                            lpd1_detector_skip_pd0(pcs_ptr, md_ctx, pic_width_in_sb);
                        }

                        // Can only use light-PD1 under the following conditions
                        if (!(md_ctx->hbd_mode_decision == 0 && md_ctx->pred_depth_only &&
                              ppcs->disallow_nsq == EB_TRUE && md_ctx->disallow_4x4 == EB_TRUE &&
                              scs_ptr->super_block_size == 64)) {
                            md_ctx->lpd1_ctrls.pd1_level = REGULAR_PD1;
                        }
                        exaustive_light_pd1_features(
                            md_ctx, ppcs, md_ctx->lpd1_ctrls.pd1_level > REGULAR_PD1, 0);
                        if (md_ctx->lpd1_ctrls.pd1_level > REGULAR_PD1)
                            signal_derivation_enc_dec_kernel_oq_light_pd1(pcs_ptr,
                                                                          context_ptr->md_context);
                        else
                            signal_derivation_enc_dec_kernel_oq(
                                scs_ptr, pcs_ptr, context_ptr->md_context);
                        if (!skip_pd_pass_0 &&
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level != MULTI_PASS_PD_OFF)
                            build_cand_block_array(
                                scs_ptr,
                                pcs_ptr,
                                context_ptr->md_context,
                                sb_index,
                                pcs_ptr->parent_pcs_ptr->sb_params_array[sb_index].is_complete_sb);
                        else
                            // Build the t=0 cand_block_array
                            build_starting_cand_block_array(
                                scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
                        // [PD_PASS_1] Mode Decision - Obtain the final partitioning decision using more accurate info
                        // than previous stages.  Reduce the total number of partitions to 1.
                        // Input : mdc_blk_ptr built @ PD0 refinement
                        // Output: md_blk_arr_nsq reduced set of block(s)

                        // PD1 MD Tool(s): default MD Tool(s)
                        if (md_ctx->lpd1_ctrls.pd1_level > REGULAR_PD1)
                            mode_decision_sb_light_pd1(scs_ptr,
                                                       pcs_ptr,
                                                       mdc_ptr,
                                                       sb_ptr,
                                                       sb_origin_x,
                                                       sb_origin_y,
                                                       sb_index,
                                                       context_ptr->md_context);
                        else
                            mode_decision_sb(scs_ptr,
                                             pcs_ptr,
                                             mdc_ptr,
                                             sb_ptr,
                                             sb_origin_x,
                                             sb_origin_y,
                                             sb_index,
                                             context_ptr->md_context);
                        //if (/*ppcs->is_used_as_reference_flag &&*/ md_ctx->hbd_mode_decision == 0 && scs_ptr->static_config.encoder_bit_depth > EB_8BIT)
                        //    md_ctx->bypass_encdec = 0;
                        // Encode Pass
                        if (!context_ptr->md_context->bypass_encdec) {
                            av1_encode_decode(scs_ptr,
                                              pcs_ptr,
                                              sb_ptr,
                                              sb_index,
                                              sb_origin_x,
                                              sb_origin_y,
                                              context_ptr);
                        }
                        av1_encdec_update(scs_ptr,
                                          pcs_ptr,
                                          sb_ptr,
                                          sb_index,
                                          sb_origin_x,
                                          sb_origin_y,
                                          context_ptr);

                        context_ptr->coded_sb_count++;
                    }
                    x_sb_start_index = (x_sb_start_index > 0) ? x_sb_start_index - 1 : 0;
                }
            }

            svt_block_on_mutex(pcs_ptr->intra_mutex);
            pcs_ptr->intra_coded_area += (uint32_t)context_ptr->tot_intra_coded_area;
            pcs_ptr->skip_coded_area += (uint32_t)context_ptr->tot_skip_coded_area;
            // Accumulate block selection
            pcs_ptr->enc_dec_coded_sb_count += (uint32_t)context_ptr->coded_sb_count;
            EbBool last_sb_flag = (pcs_ptr->sb_total_count_pix == pcs_ptr->enc_dec_coded_sb_count);
            svt_release_mutex(pcs_ptr->intra_mutex);

            if (last_sb_flag) {
                EbBool do_recode = EB_FALSE;
                if ((scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
                     scs_ptr->static_config.pass == ENC_LAST_PASS || scs_ptr->lap_enabled ||
                     scs_ptr->static_config.max_bit_rate != 0) &&
                    scs_ptr->encode_context_ptr->recode_loop != DISALLOW_RECODE) {
                    recode_loop_decision_maker(pcs_ptr, scs_ptr, &do_recode);
                }

                if (do_recode) {
                    pcs_ptr->enc_dec_coded_sb_count = 0;
                    // re-init mode decision configuration for qp update for re-encode frame
                    mode_decision_configuration_init_qp_update(pcs_ptr);
                    // init segment for re-encode frame
                    init_enc_dec_segement(pcs_ptr->parent_pcs_ptr);
                    EbObjectWrapper *enc_dec_re_encode_tasks_wrapper_ptr;
                    uint16_t         tg_count = pcs_ptr->parent_pcs_ptr->tile_group_cols *
                        pcs_ptr->parent_pcs_ptr->tile_group_rows;
                    for (uint16_t tile_group_idx = 0; tile_group_idx < tg_count; tile_group_idx++) {
                        svt_get_empty_object(context_ptr->enc_dec_feedback_fifo_ptr,
                                             &enc_dec_re_encode_tasks_wrapper_ptr);

                        EncDecTasks *enc_dec_re_encode_tasks_ptr =
                            (EncDecTasks *)enc_dec_re_encode_tasks_wrapper_ptr->object_ptr;
                        enc_dec_re_encode_tasks_ptr->pcs_wrapper_ptr =
                            enc_dec_tasks_ptr->pcs_wrapper_ptr;
                        enc_dec_re_encode_tasks_ptr->input_type       = ENCDEC_TASKS_MDC_INPUT;
                        enc_dec_re_encode_tasks_ptr->tile_group_index = tile_group_idx;

                        // Post the Full Results Object
                        svt_post_full_object(enc_dec_re_encode_tasks_wrapper_ptr);
                    }

                } else {
                    EB_FREE_ARRAY(pcs_ptr->ec_ctx_array);
                    // Copy film grain data from parent picture set to the reference object for further reference
                    if (scs_ptr->seq_header.film_grain_params_present) {
                        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
                            pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr) {
                            ((EbReferenceObject *)
                                 pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                                ->film_grain_params =
                                pcs_ptr->parent_pcs_ptr->frm_hdr.film_grain_params;
                        }
                    }
                    // Force each frame to update their data so future frames can use it,
                    // even if the current frame did not use it.  This enables REF frames to
                    // have the feature off, while NREF frames can have it on.  Used for multi-threading.
                    if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
                        pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr)
                        for (int frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame)
                            ((EbReferenceObject *)
                                 pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                                ->global_motion[frame] =
                                pcs_ptr->parent_pcs_ptr->global_motion[frame];
                    svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->sgrproj_restore_cost,
                               pcs_ptr->md_rate_estimation_array->sgrproj_restore_fac_bits,
                               2 * sizeof(int32_t));
                    svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->switchable_restore_cost,
                               pcs_ptr->md_rate_estimation_array->switchable_restore_fac_bits,
                               3 * sizeof(int32_t));
                    svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->wiener_restore_cost,
                               pcs_ptr->md_rate_estimation_array->wiener_restore_fac_bits,
                               2 * sizeof(int32_t));
                    pcs_ptr->parent_pcs_ptr->av1x->rdmult =
                        context_ptr
                            ->pic_full_lambda[(context_ptr->bit_depth == EB_10BIT) ? EB_10_BIT_MD
                                                                                   : EB_8_BIT_MD];
                    if (pcs_ptr->parent_pcs_ptr->superres_total_recode_loop == 0) {
                        svt_release_object(pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr);
                        pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr = (EbObjectWrapper *)NULL;
                        pcs_ptr->parent_pcs_ptr->pa_me_data          = NULL;
                    }
                    // Get Empty EncDec Results
                    svt_get_empty_object(context_ptr->enc_dec_output_fifo_ptr,
                                         &enc_dec_results_wrapper_ptr);
                    enc_dec_results_ptr = (EncDecResults *)enc_dec_results_wrapper_ptr->object_ptr;
                    enc_dec_results_ptr->pcs_wrapper_ptr = enc_dec_tasks_ptr->pcs_wrapper_ptr;

                    // Post EncDec Results
                    svt_post_full_object(enc_dec_results_wrapper_ptr);
                }
            }
        }
        // Release Mode Decision Results
        svt_release_object(enc_dec_tasks_wrapper_ptr);
    }
    return NULL;
}

void svt_av1_add_film_grain(EbPictureBufferDesc *src, EbPictureBufferDesc *dst,
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

    fgn_copy_rect(src->buffer_cb +
                      ((src->stride_cb * (src->origin_y >> chroma_subsamp_y) +
                        (src->origin_x >> chroma_subsamp_x))
                       << use_high_bit_depth),
                  src->stride_cb,
                  dst->buffer_cb +
                      ((dst->stride_cb * (dst->origin_y >> chroma_subsamp_y) +
                        (dst->origin_x >> chroma_subsamp_x))
                       << use_high_bit_depth),
                  dst->stride_cb,
                  dst->width >> chroma_subsamp_x,
                  dst->height >> chroma_subsamp_y,
                  use_high_bit_depth);

    fgn_copy_rect(src->buffer_cr +
                      ((src->stride_cr * (src->origin_y >> chroma_subsamp_y) +
                        (src->origin_x >> chroma_subsamp_x))
                       << use_high_bit_depth),
                  src->stride_cr,
                  dst->buffer_cr +
                      ((dst->stride_cr * (dst->origin_y >> chroma_subsamp_y) +
                        (dst->origin_x >> chroma_subsamp_x))
                       << use_high_bit_depth),
                  dst->stride_cr,
                  dst->width >> chroma_subsamp_x,
                  dst->height >> chroma_subsamp_y,
                  use_high_bit_depth);

    luma = dst->buffer_y + ((dst->origin_y * dst->stride_y + dst->origin_x) << use_high_bit_depth);
    cb   = dst->buffer_cb +
        ((dst->stride_cb * (dst->origin_y >> chroma_subsamp_y) +
          (dst->origin_x >> chroma_subsamp_x))
         << use_high_bit_depth);
    cr = dst->buffer_cr +
        ((dst->stride_cr * (dst->origin_y >> chroma_subsamp_y) +
          (dst->origin_x >> chroma_subsamp_x))
         << use_high_bit_depth);

    luma_stride   = dst->stride_y;
    chroma_stride = dst->stride_cb;

    width  = dst->width;
    height = dst->height;

    svt_av1_add_film_grain_run(&params,
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
