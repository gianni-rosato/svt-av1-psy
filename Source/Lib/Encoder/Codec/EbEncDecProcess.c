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
#if CS2_ADOPTIONS_1
#include "EbRateDistortionCost.h"
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
    EB_NEW(context_ptr->md_context,
           mode_decision_context_ctor,
           color_format,
           0,
           0,
           enable_hbd_mode_decision,
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
#if TILES_PARALLEL
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
#else
static void reset_encode_pass_neighbor_arrays(PictureControlSet *pcs_ptr) {
    neighbor_array_unit_reset(pcs_ptr->ep_intra_luma_mode_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_intra_chroma_mode_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_mv_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_skip_flag_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_mode_type_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_leaf_depth_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_luma_recon_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_recon_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_recon_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_luma_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(pcs_ptr->ep_partition_context_neighbor_array);
    // TODO(Joel): 8-bit ep_luma_recon_neighbor_array (Cb,Cr) when is_16bit==0?
    EbBool is_16bit =
        (EbBool)(pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    if (is_16bit || pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline) {
        neighbor_array_unit_reset(pcs_ptr->ep_luma_recon_neighbor_array16bit);
        neighbor_array_unit_reset(pcs_ptr->ep_cb_recon_neighbor_array16bit);
        neighbor_array_unit_reset(pcs_ptr->ep_cr_recon_neighbor_array16bit);
    }
    return;
}
#endif

/**************************************************
 * Reset Coding Loop
 **************************************************/
static void reset_enc_dec(EncDecContext *context_ptr, PictureControlSet *pcs_ptr,
                          SequenceControlSet *scs_ptr, uint32_t segment_index) {
    context_ptr->is_16bit = (EbBool)(scs_ptr->static_config.encoder_bit_depth > EB_8BIT) || (EbBool)(scs_ptr->static_config.is_16bit_pipeline);
    context_ptr->bit_depth = scs_ptr->static_config.encoder_bit_depth;
    uint16_t picture_qp   = pcs_ptr->picture_qp;
#if TILES_PARALLEL
    uint16_t tile_group_idx = context_ptr->tile_group_index;
#endif
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
    (*av1_lambda_assignment_function_table[pcs_ptr->parent_pcs_ptr->pred_structure])(
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
#if !NEW_MD_LAMBDA
        &context_ptr->fast_chroma_lambda,
        &context_ptr->full_chroma_lambda,
#endif
        (uint8_t)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index,
#if OMARK_HBD0_ED && OMARK_LAMBDA
        EB_TRUE);
#else
        context_ptr->md_context->hbd_mode_decision);
#endif
    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    if (context_ptr->is_md_rate_estimation_ptr_owner) {
        EB_FREE(context_ptr->md_rate_estimation_ptr);
        context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
    }
    context_ptr->md_rate_estimation_ptr = pcs_ptr->md_rate_estimation_array;
    if (segment_index == 0) {
#if TILES_PARALLEL
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
#else
        reset_encode_pass_neighbor_arrays(pcs_ptr);
        reset_segmentation_map(pcs_ptr->segmentation_neighbor_map);
#endif
    }

    return;
}

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
#if !NEW_MD_LAMBDA
        &context_ptr->fast_chroma_lambda,
        &context_ptr->full_chroma_lambda,
#endif
        (uint8_t)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index,
#if OMARK_HBD0_ED && OMARK_LAMBDA
        EB_TRUE);
#else
        context_ptr->md_context->hbd_mode_decision);
#endif

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
    EbObjectWrapper *wrapper_ptr;
    EncDecTasks *    feedback_task_ptr;

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
                    self_assigned            = EB_TRUE;
                    continue_processing_flag = EB_TRUE;

                    //fprintf(trace, "Start  Pic: %u Seg: %u\n",
                    //    (unsigned) ((PictureControlSet*) taskPtr->pcs_wrapper_ptr->object_ptr)->picture_number,
                    //    *segmentInOutIndex);
                }
            }
            eb_release_mutex(segmentPtr->row_array[row_segment_index + 1].assignment_mutex);
        }

        if (feedback_row_index > 0) {
            eb_get_empty_object(srmFifoPtr, &wrapper_ptr);
            feedback_task_ptr                      = (EncDecTasks *)wrapper_ptr->object_ptr;
            feedback_task_ptr->input_type          = ENCDEC_TASKS_ENCDEC_INPUT;
            feedback_task_ptr->enc_dec_segment_row = feedback_row_index;
            feedback_task_ptr->pcs_wrapper_ptr     = taskPtr->pcs_wrapper_ptr;
#if TILES_PARALLEL
            feedback_task_ptr->tile_group_index = taskPtr->tile_group_index;
#endif
            eb_post_full_object(wrapper_ptr);
        }

        break;

    default: break;
    }

    return continue_processing_flag;
}
void recon_output(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    EbObjectWrapper *   output_recon_wrapper_ptr;
    EbBufferHeaderType *output_recon_ptr;
    EncodeContext *     encode_context_ptr = scs_ptr->encode_context_ptr;
    EbBool              is_16bit           = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    // The totalNumberOfReconFrames counter has to be write/read protected as
    //   it is used to determine the end of the stream.  If it is not protected
    //   the encoder might not properly terminate.
    eb_block_on_mutex(encode_context_ptr->total_number_of_recon_frame_mutex);

    if (!pcs_ptr->parent_pcs_ptr->is_alt_ref) {
        // Get Recon Buffer
        eb_get_empty_object(scs_ptr->encode_context_ptr->recon_output_fifo_ptr,
                            &output_recon_wrapper_ptr);
        output_recon_ptr        = (EbBufferHeaderType *)output_recon_wrapper_ptr->object_ptr;
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

void psnr_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
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
        uint32_t column_index;
        uint32_t row_index           = 0;
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

        while (row_index < (uint32_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom)) {
            column_index = 0;
            while (column_index < (uint32_t)(input_picture_ptr->width - scs_ptr->max_input_pad_right)) {
                residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                    (recon_coeff_buffer[column_index]));
                ++column_index;
            }

            input_buffer += input_picture_ptr->stride_y;
            recon_coeff_buffer += recon_ptr->stride_y;
            ++row_index;
        }

        sse_total[0] = residual_distortion;

        recon_coeff_buffer =
            &((recon_ptr->buffer_cb)[recon_ptr->origin_x / 2 +
                                     recon_ptr->origin_y / 2 * recon_ptr->stride_cb]);
        input_buffer = &(buffer_cb[input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);

        residual_distortion = 0;
        row_index           = 0;
        while (row_index < (uint32_t)((input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y)) {
            column_index = 0;
            while (column_index < (uint32_t)((input_picture_ptr->width - scs_ptr->max_input_pad_bottom) >> ss_x)) {
                residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                    (recon_coeff_buffer[column_index]));
                ++column_index;
            }

            input_buffer += input_picture_ptr->stride_cb;
            recon_coeff_buffer += recon_ptr->stride_cb;
            ++row_index;
        }

        sse_total[1] = residual_distortion;

        recon_coeff_buffer =
            &((recon_ptr->buffer_cr)[recon_ptr->origin_x / 2 +
                                     recon_ptr->origin_y / 2 * recon_ptr->stride_cr]);
        input_buffer        = &(buffer_cr[input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
        residual_distortion = 0;
        row_index           = 0;

        while (row_index < (uint32_t)((input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y)) {
            column_index = 0;
            while (column_index < (uint32_t)((input_picture_ptr->width - scs_ptr->max_input_pad_bottom) >> ss_x)) {
                residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                    (recon_coeff_buffer[column_index]));
                ++column_index;
            }

            input_buffer += input_picture_ptr->stride_cr;
            recon_coeff_buffer += recon_ptr->stride_cr;
            ++row_index;
        }

        sse_total[2]                      = residual_distortion;
        pcs_ptr->parent_pcs_ptr->luma_sse = (uint32_t)sse_total[0];
        pcs_ptr->parent_pcs_ptr->cb_sse   = (uint32_t)sse_total[1];
        pcs_ptr->parent_pcs_ptr->cr_sse   = (uint32_t)sse_total[2];

        if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            EB_FREE_ARRAY(buffer_y);
            EB_FREE_ARRAY(buffer_cb);
            EB_FREE_ARRAY(buffer_cr);
        }
    } else {
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
        uint32_t  column_index;
        uint32_t  row_index           = 0;
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

            while (row_index < (uint32_t)(input_picture_ptr->height - scs_ptr->max_input_pad_bottom)) {
                column_index = 0;
                while (column_index < (uint32_t)(input_picture_ptr->width - scs_ptr->max_input_pad_right)) {
                    residual_distortion +=
                        (int64_t)SQR((int64_t)((((input_buffer[column_index]) << 2) |
                                                ((input_buffer_bit_inc[column_index] >> 6) & 3))) -
                                     (recon_coeff_buffer[column_index]));

                    ++column_index;
                }

                input_buffer += input_picture_ptr->stride_y;
                input_buffer_bit_inc += input_picture_ptr->stride_bit_inc_y;
                recon_coeff_buffer += recon_ptr->stride_y;
                ++row_index;
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
            row_index           = 0;
            while (row_index < (uint32_t)((input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y)) {
                column_index = 0;
                while (column_index < (uint32_t)((input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x)) {
                    residual_distortion +=
                        (int64_t)SQR((int64_t)((((input_buffer[column_index]) << 2) |
                                                ((input_buffer_bit_inc[column_index] >> 6) & 3))) -
                                     (recon_coeff_buffer[column_index]));
                    ++column_index;
                }

                input_buffer += input_picture_ptr->stride_cb;
                input_buffer_bit_inc += input_picture_ptr->stride_bit_inc_cb;
                recon_coeff_buffer += recon_ptr->stride_cb;
                ++row_index;
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
            row_index           = 0;

            while (row_index < (uint32_t)((input_picture_ptr->height - scs_ptr->max_input_pad_bottom) >> ss_y)) {
                column_index = 0;
                while (column_index < (uint32_t)((input_picture_ptr->width - scs_ptr->max_input_pad_right) >> ss_x)) {
                    residual_distortion +=
                        (int64_t)SQR((int64_t)((((input_buffer[column_index]) << 2) |
                                                ((input_buffer_bit_inc[column_index] >> 6) & 3))) -
                                     (recon_coeff_buffer[column_index]));
                    ++column_index;
                }

                input_buffer += input_picture_ptr->stride_cr;
                input_buffer_bit_inc += input_picture_ptr->stride_bit_inc_cr;
                recon_coeff_buffer += recon_ptr->stride_cr;
                ++row_index;
            }

            sse_total[2] = residual_distortion;

            if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
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
        // Y samples
        generate_padding(ref_pic_ptr->buffer_y,
                         ref_pic_ptr->stride_y,
                         ref_pic_ptr->width - scs_ptr->max_input_pad_right,
                         ref_pic_ptr->height - scs_ptr->max_input_pad_bottom,
                         ref_pic_ptr->origin_x,
                         ref_pic_ptr->origin_y);

        // Cb samples
        generate_padding(ref_pic_ptr->buffer_cb,
                         ref_pic_ptr->stride_cb,
                         (ref_pic_ptr->width - scs_ptr->max_input_pad_right) >> 1,
                         (ref_pic_ptr->height - scs_ptr->max_input_pad_bottom) >> 1,
                         ref_pic_ptr->origin_x >> 1,
                         ref_pic_ptr->origin_y >> 1);

        // Cr samples
        generate_padding(ref_pic_ptr->buffer_cr,
                         ref_pic_ptr->stride_cr,
                         (ref_pic_ptr->width - scs_ptr->max_input_pad_right) >> 1,
                         (ref_pic_ptr->height - scs_ptr->max_input_pad_bottom) >> 1,
                         ref_pic_ptr->origin_x >> 1,
                         ref_pic_ptr->origin_y >> 1);
    }

    //We need this for MCP
    if (is_16bit) {
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
        memcpy(((EbReferenceObject *)
                    pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                   ->ref_order_hint,
               pcs_ptr->parent_pcs_ptr->ref_order_hint,
               7 * sizeof(uint32_t));
    }
}

/******************************************************
* Derive EncDec Settings for OQ
Input   : encoder mode and pd pass
Output  : EncDec Kernel signal(s)
******************************************************/
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
        if (pcs_ptr->enc_mode <= ENC_M6)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    else if (pcs_ptr->enc_mode <= ENC_M4)
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
    else if (pcs_ptr->enc_mode <= ENC_M7) {
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    } else
        context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;

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
#if CS2_ADOPTIONS_1
        else if (pcs_ptr->parent_pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M0)
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else if (pcs_ptr->enc_mode <= ENC_M0)
            context_ptr->tx_weight = MAX_MODE_COST;
        else if (pcs_ptr->enc_mode <= ENC_M1 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected))
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
#else
        else if (!MR_MODE && pcs_ptr->enc_mode <= ENC_M5)
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else if (!MR_MODE) {
            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
            else
                context_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
        }
#endif
    }

    // Set tx search reduced set falg (0: full tx set; 1: reduced tx set; 1: two tx))
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->tx_search_reduced_set = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->tx_search_reduced_set = 0;
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        if (pcs_ptr->enc_mode <= ENC_M5)
            context_ptr->tx_search_reduced_set = 0;
        else if (pcs_ptr->enc_mode <= ENC_M6)
            if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
                context_ptr->tx_search_reduced_set = 0;
            else
                context_ptr->tx_search_reduced_set = 1;
        else if (pcs_ptr->enc_mode <= ENC_M7)
            context_ptr->tx_search_reduced_set = 1;
        else
            context_ptr->tx_search_reduced_set = 2;

    else if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
        context_ptr->tx_search_reduced_set = 0;
    else if (pcs_ptr->enc_mode <= ENC_M2)
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
#if ENHANCED_MULTI_PASS_PD_MD_STAGING_SETTINGS
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#else
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
        else
            context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#endif
    } else if (MR_MODE)
        context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP;
    else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
    else if (pcs_ptr->enc_mode <= ENC_M2)

        context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
    else if (pcs_ptr->enc_mode <= ENC_M7)
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
        else
            context_ptr->interpolation_search_level = IT_SEARCH_OFF;
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
#if ENHANCED_MULTI_PASS_PD_MD_STAGING_SETTINGS
        context_ptr->chroma_level = CHROMA_MODE_1;
#else
        if (pcs_ptr->temporal_layer_index == 0)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else
            context_ptr->chroma_level = CHROMA_MODE_1;
#endif
    } else if (scs_ptr->static_config.set_chroma_mode == DEFAULT) {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
            if (pcs_ptr->enc_mode <= ENC_M6)
                context_ptr->chroma_level = CHROMA_MODE_1;
            else if (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                context_ptr->chroma_level = CHROMA_MODE_1;
            else
                context_ptr->chroma_level =
                    (scs_ptr->encoder_bit_depth == EB_8BIT) ? CHROMA_MODE_2 : CHROMA_MODE_3;
#if CS2_ADOPTIONS_1
        else if (pcs_ptr->enc_mode <= ENC_M1)
#else
        else if (pcs_ptr->enc_mode == ENC_M0)
#endif
            context_ptr->chroma_level = CHROMA_MODE_0;
        else if (pcs_ptr->enc_mode <= ENC_M5 && pcs_ptr->temporal_layer_index == 0)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else if (pcs_ptr->enc_mode <= ENC_M5)
            context_ptr->chroma_level = CHROMA_MODE_1;
        else
            context_ptr->chroma_level =
                (scs_ptr->encoder_bit_depth == EB_8BIT) ? CHROMA_MODE_2 : CHROMA_MODE_3;

    } else // use specified level
        context_ptr->chroma_level = scs_ptr->static_config.set_chroma_mode;
#if MOVE_OPT
    // Chroma independent modes search
    // Level                Settings
    // 0                    post first md_stage
    // 1                    post last md_stage
#if CS2_ADOPTIONS_1
    context_ptr->chroma_at_last_md_stage =
        MR_MODE ? 0 : (context_ptr->chroma_level == CHROMA_MODE_0 && !pcs_ptr->parent_pcs_ptr->sc_content_detected) ? 1 : 0;
#else
    context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) /*&& (pcs_ptr->enc_mode == ENC_M1)*/ ? 1 : 0;
#endif
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
        if (pcs_ptr->enc_mode <= ENC_M5)
            context_ptr->full_loop_escape = 0;
        else
            context_ptr->full_loop_escape = 2;
    else if (pcs_ptr->enc_mode <= ENC_M5)
        context_ptr->full_loop_escape = 0;
    else
        context_ptr->full_loop_escape = 2;

        // Set global MV injection
        // Level                Settings
        // 0                    Injection off (Hsan: but not derivation as used by MV ref derivation)
        // 1                    On
#if GLOBAL_WARPED_MOTION
    if (scs_ptr->static_config.enable_global_motion == EB_TRUE) {
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->global_mv_injection = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->global_mv_injection = 0;
        else if (pcs_ptr->enc_mode <= ENC_M1)
            context_ptr->global_mv_injection = 1;
        else
            context_ptr->global_mv_injection = 0;
    } else
        context_ptr->global_mv_injection = 0;
#else
    if (scs_ptr->static_config.enable_global_warped_motion == EB_TRUE) {
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            if (pcs_ptr->enc_mode <= ENC_M1)
                context_ptr->global_mv_injection = 1;
            else
                context_ptr->global_mv_injection = 0;
        } else {
            if (pcs_ptr->enc_mode <= ENC_M7)
                context_ptr->global_mv_injection = 1;
            else
                context_ptr->global_mv_injection = 0;
        }
    } else
        context_ptr->global_mv_injection = 0;
#endif
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
#if CS2_ADOPTIONS_1
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
                context_ptr->new_nearest_near_comb_injection = 0;
        else if (pcs_ptr->enc_mode <= ENC_M0)
            context_ptr->new_nearest_near_comb_injection = 1;
        else
            context_ptr->new_nearest_near_comb_injection = 0;
#else
#if !MR_MODE
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
            context_ptr->new_nearest_near_comb_injection = 0;
        else
#endif
            if (pcs_ptr->enc_mode <= ENC_M0 && pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
            context_ptr->new_nearest_near_comb_injection = 1;
        else
            context_ptr->new_nearest_near_comb_injection = 0;
 #endif
    else
        context_ptr->new_nearest_near_comb_injection =
            scs_ptr->static_config.new_nearest_comb_inject;
#if !ENHANCED_ME_MV
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->nx4_4xn_parent_mv_injection = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->nx4_4xn_parent_mv_injection = 0;
    else if (scs_ptr->static_config.nx4_4xn_parent_mv_inject == DEFAULT)
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
            if (pcs_ptr->enc_mode <= ENC_M2)
                context_ptr->nx4_4xn_parent_mv_injection = 1;
            else
                context_ptr->nx4_4xn_parent_mv_injection = 0;
        else if (pcs_ptr->enc_mode <= ENC_M7)
            context_ptr->nx4_4xn_parent_mv_injection = 1;
        else
            context_ptr->nx4_4xn_parent_mv_injection = 0;

    else
        context_ptr->nx4_4xn_parent_mv_injection = scs_ptr->static_config.nx4_4xn_parent_mv_inject;
#endif
    // Set warped motion injection
    // Level                Settings
    // 0                    OFF
    // 1                    On
#if WARP_IMPROVEMENT
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
#else
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->warped_motion_injection = 0;
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->warped_motion_injection = 1;
    } else if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->warped_motion_injection = 0;
    else
        context_ptr->warped_motion_injection = 1;
#endif
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
#if CS2_ADOPTIONS_1
    if (pcs_ptr->enc_mode <= ENC_M7)
        context_ptr->unipred3x3_injection = 1;
    else
        context_ptr->unipred3x3_injection = 0;
#else
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        if (MR_MODE)
            context_ptr->unipred3x3_injection = 1;
        else if (pcs_ptr->enc_mode <= ENC_M2)
            context_ptr->unipred3x3_injection = 2;
        else
            context_ptr->unipred3x3_injection = 0;
    else if (pcs_ptr->enc_mode <= ENC_M7)
        context_ptr->unipred3x3_injection = 1;
    else
        context_ptr->unipred3x3_injection = 0;
#endif
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
            if (pcs_ptr->enc_mode <= ENC_M4)
                context_ptr->bipred3x3_injection = 1;
            else
                context_ptr->bipred3x3_injection = 0;
        else if (pcs_ptr->enc_mode <= ENC_M2)
            context_ptr->bipred3x3_injection = 1;
        else if (pcs_ptr->enc_mode <= ENC_M4)
            context_ptr->bipred3x3_injection = 2;
        else
            context_ptr->bipred3x3_injection = 0;

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
#if CS2_ADOPTIONS_1
                if (pcs_ptr->enc_mode <= ENC_M4)
                    context_ptr->predictive_me_level = 1;
#else
                if (pcs_ptr->enc_mode <= ENC_M1)
                    if (!MR_MODE)
                        context_ptr->predictive_me_level = (pcs_ptr->enc_mode <= ENC_M0) ? 2 : 4;
                    else
                        context_ptr->predictive_me_level = 4;

                else if (pcs_ptr->enc_mode <= ENC_M4)
                    context_ptr->predictive_me_level = 2;
#endif
                else
                    context_ptr->predictive_me_level = 0;
#if CS2_ADOPTIONS_1
            else if (pcs_ptr->enc_mode <= ENC_M0)
                context_ptr->predictive_me_level = 6;
#endif
            else if (pcs_ptr->enc_mode <= ENC_M2)
                context_ptr->predictive_me_level = 5;
            else if (pcs_ptr->enc_mode <= ENC_M4)
                context_ptr->predictive_me_level = 2;
            else
                context_ptr->predictive_me_level = 0;

        } else
            context_ptr->predictive_me_level = scs_ptr->static_config.pred_me;
    } else
        context_ptr->predictive_me_level = 0;

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
#if CS2_ADOPTIONS_1
    else if (pcs_ptr->enc_mode <= ENC_M5)
#else
    else if (pcs_ptr->enc_mode <= ENC_M4)
#endif
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    else
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;

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

        if (pcs_ptr->enc_mode <= ENC_M4)
        context_ptr->interpolation_filter_search_blk_size = 0;
    else
        context_ptr->interpolation_filter_search_blk_size = 1;

    // Set PF MD
    context_ptr->pf_md_mode = PF_OFF;
    // Derive Spatial SSE Flag
#if CS2_ADOPTIONS_1
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->spatial_sse_full_loop = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->spatial_sse_full_loop = EB_FALSE;
#else
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->spatial_sse_full_loop = EB_TRUE;
    else if (context_ptr->pd_pass == PD_PASS_1)
#endif
    else if (scs_ptr->static_config.spatial_sse_fl == DEFAULT)
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
            if (pcs_ptr->enc_mode <= ENC_M6)
                context_ptr->spatial_sse_full_loop = EB_TRUE;
            else
                context_ptr->spatial_sse_full_loop = EB_FALSE;
        else if (pcs_ptr->enc_mode <= ENC_M4)
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
            if (pcs_ptr->enc_mode <= ENC_M2)
                context_ptr->enable_rdoq = EB_TRUE;
            else
                context_ptr->enable_rdoq = EB_FALSE;
        else if (pcs_ptr->enc_mode <= ENC_M2)
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
        else if (pcs_ptr->enc_mode <= ENC_M5)
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
#if CS2_ADOPTIONS_1
            if (MR_MODE)
                context_ptr->edge_based_skip_angle_intra = 0;
            else
                context_ptr->edge_based_skip_angle_intra = 1;
#else
            if (MR_MODE)
                context_ptr->edge_based_skip_angle_intra = 0;
            else if (pcs_ptr->enc_mode <= ENC_M7 && !pcs_ptr->parent_pcs_ptr->sc_content_detected)
                context_ptr->edge_based_skip_angle_intra = 0;
            else
                context_ptr->edge_based_skip_angle_intra = 1;
#endif
        } else
            context_ptr->edge_based_skip_angle_intra = scs_ptr->static_config.edge_skp_angle_intra;
    else
        context_ptr->edge_based_skip_angle_intra = 0;
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->prune_ref_frame_for_rec_partitions = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->prune_ref_frame_for_rec_partitions = 1;
    else if (scs_ptr->static_config.prune_ref_rec_part == DEFAULT)
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected || pcs_ptr->enc_mode <= ENC_M1)
            context_ptr->prune_ref_frame_for_rec_partitions = 0;
        else
            context_ptr->prune_ref_frame_for_rec_partitions = 1;
    else
        context_ptr->prune_ref_frame_for_rec_partitions = scs_ptr->static_config.prune_ref_rec_part;

    // Derive INTER/INTER WEDGE variance TH
    // Phoenix: Active only when inter/inter compound is on
#if CS2_ADOPTIONS_1
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
    if (pcs_ptr->enc_mode <= ENC_M7)
#endif
        context_ptr->inter_inter_wedge_variance_th = 0;
    else
        context_ptr->inter_inter_wedge_variance_th = 100;

    // Derive MD Exit TH
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_exit_th = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
#if CS2_ADOPTIONS_1
        context_ptr->md_exit_th = 18;
#else
        context_ptr->md_exit_th = (pcs_ptr->parent_pcs_ptr->sc_content_detected) ? 10 : 18;
#endif
#if CS2_ADOPTIONS_1
    else if (MR_MODE)
        context_ptr->md_exit_th = 0;
     else if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->enc_mode <= ENC_M1 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
        context_ptr->md_exit_th = 0;
    else
        context_ptr->md_exit_th = 18;
#else
    else if (MR_MODE ||
             (pcs_ptr->enc_mode == ENC_M0 && pcs_ptr->parent_pcs_ptr->sc_content_detected == 0))
        context_ptr->md_exit_th = 0;
    else
        context_ptr->md_exit_th = (pcs_ptr->parent_pcs_ptr->sc_content_detected) ? 10 : 18;

#endif
    // md_stage_1_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_1_cand_prune_th
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_1_cand_prune_th = 75;
#if CS2_ADOPTIONS_1
    else if (pcs_ptr->enc_mode <= ENC_M1 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
    else if (MR_MODE ||
             (pcs_ptr->enc_mode == ENC_M0 && (pcs_ptr->parent_pcs_ptr->sc_content_detected == 0)) ||
             scs_ptr->input_resolution == INPUT_SIZE_576p_RANGE_OR_LOWER)
#endif
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
    else if (pcs_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_1_cand_prune_th = scs_ptr->static_config.md_stage_1_cand_prune_th;
    else
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;

    // md_stage_1_class_prune_th (for class removal)
    // Remove class if deviation to the best higher than TH_C
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_1_class_prune_th = 100;
#if CS2_ADOPTIONS_1
        if (pcs_ptr->enc_mode <= ENC_M1 ||
            pcs_ptr->parent_pcs_ptr->sc_content_detected)
#else
    else if (MR_MODE ||
             (pcs_ptr->enc_mode == ENC_M0 && (pcs_ptr->parent_pcs_ptr->sc_content_detected == 0)) ||
             scs_ptr->input_resolution == INPUT_SIZE_576p_RANGE_OR_LOWER)
#endif
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;

#if CS2_ADOPTIONS_1
    else
        context_ptr->md_stage_1_class_prune_th = scs_ptr->static_config.md_stage_1_class_prune_th;
#else
    else if (pcs_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_1_class_prune_th = scs_ptr->static_config.md_stage_1_class_prune_th;
    else
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
#endif
    // md_stage_2_3_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_2_3_cand_prune_th
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
#if CS2_ADOPTIONS_1
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_3_cand_prune_th = 5;
    else if (MR_MODE)
        context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
    else if (pcs_ptr->enc_mode <= ENC_M1 || pcs_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->md_stage_2_3_cand_prune_th = 15;
#else
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_3_cand_prune_th =
        scs_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE ? 5 : 3;
    else if (MR_MODE || pcs_ptr->parent_pcs_ptr->sc_content_detected || pcs_ptr->enc_mode <= ENC_M0)
        context_ptr->md_stage_2_3_cand_prune_th = (uint64_t)~0;
#endif
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
#if CS2_ADOPTIONS_1
    else if ((pcs_ptr->enc_mode <= ENC_M3 && pcs_ptr->parent_pcs_ptr->sc_content_detected))
        context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
    else
        context_ptr->md_stage_2_3_class_prune_th = scs_ptr->static_config.md_stage_2_3_class_prune_th;
#else
    else if (MR_MODE)
        context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
    else if (pcs_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_2_3_class_prune_th = scs_ptr->static_config.md_stage_2_3_class_prune_th;
    else // to be tested for m5-m8
        context_ptr->md_stage_2_3_class_prune_th = (uint64_t)~0;
#endif
    // Weighting (expressed as a percentage) applied to
    // square shape costs for determining if a and b
    // shapes should be skipped. Namely:
    // skip HA and HB if h_cost > (weighted sq_cost)
    // skip VA and VB if v_cost > (weighted sq_cost)

    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->sq_weight = (uint32_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->sq_weight = 100;
#if CS2_ADOPTIONS_1
    if (MR_MODE)
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight + 15;
    else if (pcs_ptr->enc_mode <= ENC_M0 ||
        (pcs_ptr->enc_mode <= ENC_M1 && !(pcs_ptr->parent_pcs_ptr->sc_content_detected)))
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight + 5;
    else
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight - 5;
#else
    else if (MR_MODE)
        context_ptr->sq_weight = (uint32_t)~0;
#if ENHANCED_SQ_WEIGHT
    else if (pcs_ptr->enc_mode <= ENC_M1)
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight + 5;
    else if (pcs_ptr->enc_mode <= ENC_M2)
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight;
    else
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight - 5;
#else
    else
        context_ptr->sq_weight = scs_ptr->static_config.sq_weight;
#endif
#endif
#if NSQ_HV
    // nsq_hv_level  needs sq_weight to be ON
    // 0: OFF
    // 1: ON
    //        if H > V + TH1% then skip HA / HB / H4
    //        if V > H + TH1% then skip VA / VB / V4
    // 2: ON  for  H4 / V4 use more agressive TH2% for faster mode
    if (MR_MODE || context_ptr->pd_pass < PD_PASS_2)
        context_ptr->nsq_hv_level = 0;
    else if (pcs_ptr->enc_mode <= ENC_M3) {
        context_ptr->nsq_hv_level = 1;
        assert(context_ptr->sq_weight != (uint32_t)~0);
    }
    else {
        context_ptr->nsq_hv_level = 2;
        assert(context_ptr->sq_weight != (uint32_t)~0);
    }
#endif
    // Set pred ME full search area
    if (context_ptr->pd_pass == PD_PASS_0) {
#if CS2_ADOPTIONS_1
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15;
        }
#else
        context_ptr->pred_me_full_pel_search_width  = PRED_ME_FULL_PEL_SEARCH_WIDTH;
        context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_SEARCH_HEIGHT;
#endif
    } else if (context_ptr->pd_pass == PD_PASS_1) {
#if CS2_ADOPTIONS_1
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15;
        }
#else
        context_ptr->pred_me_full_pel_search_width  = PRED_ME_FULL_PEL_SEARCH_WIDTH;
        context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_SEARCH_HEIGHT;
#endif
    } else {
#if CS2_ADOPTIONS_1
        if (pcs_ptr->parent_pcs_ptr->sc_content_detected) {
            context_ptr->pred_me_full_pel_search_width = PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_7;
        }
        else {
            context_ptr->pred_me_full_pel_search_width = pcs_ptr->enc_mode <= ENC_M0 ?
                PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_15 : PRED_ME_FULL_PEL_REF_WINDOW_WIDTH_7;
            context_ptr->pred_me_full_pel_search_height = pcs_ptr->enc_mode <= ENC_M0 ?
                PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_15 : PRED_ME_FULL_PEL_REF_WINDOW_HEIGHT_5;
        }
#else
        context_ptr->pred_me_full_pel_search_width =
            (pcs_ptr->parent_pcs_ptr->sc_content_detected == 0 && pcs_ptr->enc_mode == ENC_M0)
                ? PRED_ME_FULL_PEL_SEARCH_WIDTH_EXTENDED
                : PRED_ME_FULL_PEL_SEARCH_WIDTH;
        context_ptr->pred_me_full_pel_search_height =
            (pcs_ptr->parent_pcs_ptr->sc_content_detected == 0 && pcs_ptr->enc_mode == ENC_M0)
                ? PRED_ME_FULL_PEL_SEARCH_HEIGHT_EXTENDED
                : PRED_ME_FULL_PEL_SEARCH_HEIGHT;
#endif
    }
#if COMP_SIMILAR
    //comp_similar_mode
    //0: OFF
    //1: If previous similar block is not compound, do not inject compound
    //2: If previous similar block is not compound, do not inject compound
    //   else consider the compound modes up the mode for the similar block
    if (pcs_ptr->enc_mode <= ENC_M3)
        context_ptr->comp_similar_mode = 1;
    else
        context_ptr->comp_similar_mode = 2;
#else

    // set compound_types_to_try
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->compound_types_to_try = MD_COMP_AVG;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->compound_types_to_try = MD_COMP_AVG;
    else {
        if (pcs_ptr->parent_pcs_ptr->compound_mode)
            context_ptr->compound_types_to_try =
                pcs_ptr->parent_pcs_ptr->compound_mode == 1 ? MD_COMP_DIFF0 : MD_COMP_WEDGE;
        else
            context_ptr->compound_types_to_try = MD_COMP_AVG;
    }
#endif

#if  INTRA_SIMILAR
    //intra_similar_mode
    //0: OFF
    //1: If previous similar block is intra, do not inject any inter
    context_ptr->intra_similar_mode = 1;
#endif

    // Set coeff_based_nsq_cand_reduction
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
#if CS2_ADOPTIONS_1
    else if (MR_MODE &&
        pcs_ptr->parent_pcs_ptr->sc_content_detected == 0)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
#endif
    else
        context_ptr->coeff_based_nsq_cand_reduction = EB_TRUE;

    // Set pic_obmc_mode @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_pic_obmc_mode = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_pic_obmc_mode = 0;
    else
        context_ptr->md_pic_obmc_mode = pcs_ptr->parent_pcs_ptr->pic_obmc_mode;

    // Set enable_inter_intra @ MD
#if  CLEANUP_INTER_INTRA
    //Block level switch, has to follow the picture level
#endif
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
    if(MR_MODE || context_ptr->pd_pass <= PD_PASS_1)
        context_ptr->skip_depth = 0;
    else
        context_ptr->skip_depth =
        pcs_ptr->parent_pcs_ptr->sc_content_detected ? 1 : 0;

#if ENHANCED_ME_MV
    // Set perform_me_mv_1_8_pel_ref
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->perform_me_mv_1_8_pel_ref = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->perform_me_mv_1_8_pel_ref = EB_FALSE;
    else
        context_ptr->perform_me_mv_1_8_pel_ref = (pcs_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv);
#endif
#if NICS_CLEANUP
    // Set nic_level for PD2 only
    // nic_level        nic scale factor
    // 0                1
    // 1                3/4
    // 2                2/3
    // 3                1/2
    if (pcs_ptr->parent_pcs_ptr->sc_content_detected)
        if (pcs_ptr->enc_mode <= ENC_M1)
#if CS2_ADOPTIONS_1
            context_ptr->nic_level = 0;
#else
            context_ptr->nic_level = 1;
#endif
        else
            context_ptr->nic_level = 2;
    else
#if CS2_ADOPTIONS_1
        if (pcs_ptr->enc_mode <= ENC_M1)
            context_ptr->nic_level = 0;
#else
        if (pcs_ptr->enc_mode <= ENC_M0)
            context_ptr->nic_level = 1;
#endif
        else if (pcs_ptr->enc_mode <= ENC_M1)
            context_ptr->nic_level = 2;
        else
            context_ptr->nic_level = 3;
#endif

#if CS2_ADOPTIONS_1
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
#endif
    return return_error;
}
void copy_neighbour_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                           uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds, uint32_t sb_org_x,
                           uint32_t sb_org_y);

static void set_parent_to_be_considered(MdcSbData *results_ptr, uint32_t blk_index, int32_t sb_size,
                                        int8_t depth_step) {
    uint32_t         parent_depth_idx_mds, block_1d_idx;
    const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
    if (blk_geom->sq_size < ((sb_size == BLOCK_128X128) ? 128 : 64)) {
        //Set parent to be considered
        parent_depth_idx_mds =
            (blk_geom->sqi_mds -
             (blk_geom->quadi - 3) * ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth]) -
            parent_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
        const BlockGeom *parent_blk_geom = get_blk_geom_mds(parent_depth_idx_mds);
        uint32_t         parent_tot_d1_blocks =
            parent_blk_geom->sq_size == 128
                ? 17
                : parent_blk_geom->sq_size > 8 ? 25 : parent_blk_geom->sq_size == 8 ? 5 : 1;
        for (block_1d_idx = 0; block_1d_idx < parent_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[parent_depth_idx_mds + block_1d_idx].consider_block = 1;
        }

        if (depth_step < -1)
            set_parent_to_be_considered(results_ptr, parent_depth_idx_mds, sb_size, depth_step + 1);
    }
}

static void set_child_to_be_considered(MdcSbData *results_ptr, uint32_t blk_index, int32_t sb_size,
                                       int8_t depth_step) {
    uint32_t         child_block_idx_1, child_block_idx_2, child_block_idx_3, child_block_idx_4;
    uint32_t         tot_d1_blocks, block_1d_idx;
    const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
    tot_d1_blocks =
        blk_geom->sq_size == 128 ? 17 : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;
    if (blk_geom->sq_size > 4) {
        for (block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[blk_index + block_1d_idx].consider_block     = 1;
            results_ptr->leaf_data_array[blk_index + block_1d_idx].refined_split_flag = EB_TRUE;
        }
        //Set first child to be considered
        child_block_idx_1 = blk_index + d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
        const BlockGeom *child1_blk_geom = get_blk_geom_mds(child_block_idx_1);
        uint32_t         child1_tot_d1_blocks =
            child1_blk_geom->sq_size == 128
                ? 17
                : child1_blk_geom->sq_size > 8 ? 25 : child1_blk_geom->sq_size == 8 ? 5 : 1;

        for (block_1d_idx = 0; block_1d_idx < child1_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].consider_block = 1;
            results_ptr->leaf_data_array[child_block_idx_1 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
        if (depth_step > 1)
            set_child_to_be_considered(results_ptr, child_block_idx_1, sb_size, depth_step - 1);
        //Set second child to be considered
        child_block_idx_2 =
            child_block_idx_1 + ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom *child2_blk_geom = get_blk_geom_mds(child_block_idx_2);
        uint32_t         child2_tot_d1_blocks =
            child2_blk_geom->sq_size == 128
                ? 17
                : child2_blk_geom->sq_size > 8 ? 25 : child2_blk_geom->sq_size == 8 ? 5 : 1;
        for (block_1d_idx = 0; block_1d_idx < child2_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].consider_block = 1;
            results_ptr->leaf_data_array[child_block_idx_2 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
        if (depth_step > 1)
            set_child_to_be_considered(results_ptr, child_block_idx_2, sb_size, depth_step - 1);
        //Set third child to be considered
        child_block_idx_3 =
            child_block_idx_2 + ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom *child3_blk_geom = get_blk_geom_mds(child_block_idx_3);
        uint32_t         child3_tot_d1_blocks =
            child3_blk_geom->sq_size == 128
                ? 17
                : child3_blk_geom->sq_size > 8 ? 25 : child3_blk_geom->sq_size == 8 ? 5 : 1;

        for (block_1d_idx = 0; block_1d_idx < child3_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].consider_block = 1;
            results_ptr->leaf_data_array[child_block_idx_3 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
        if (depth_step > 1)
            set_child_to_be_considered(results_ptr, child_block_idx_3, sb_size, depth_step - 1);
        //Set forth child to be considered
        child_block_idx_4 =
            child_block_idx_3 + ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom *child4_blk_geom = get_blk_geom_mds(child_block_idx_4);
        uint32_t         child4_tot_d1_blocks =
            child4_blk_geom->sq_size == 128
                ? 17
                : child4_blk_geom->sq_size > 8 ? 25 : child4_blk_geom->sq_size == 8 ? 5 : 1;
        for (block_1d_idx = 0; block_1d_idx < child4_tot_d1_blocks; block_1d_idx++) {
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].consider_block = 1;
            results_ptr->leaf_data_array[child_block_idx_4 + block_1d_idx].refined_split_flag =
                EB_FALSE;
        }
        if (depth_step > 1)
            set_child_to_be_considered(results_ptr, child_block_idx_4, sb_size, depth_step - 1);
    }
}

static void build_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                   uint32_t sb_index) {
    MdcSbData *results_ptr  = &pcs_ptr->mdc_sb_array[sb_index];
    results_ptr->leaf_count = 0;
    uint32_t blk_index      = 0;
    uint32_t d1_blocks_accumlated, tot_d1_blocks = 0, d1_block_idx;

    EbBool split_flag;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
        split_flag                = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
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
        blk_index +=
            split_flag
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] -
                      tot_d1_blocks
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] -
                      tot_d1_blocks;
    }
}

#if CS2_ADOPTIONS_1
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
#else
uint64_t pd_level_tab[2][9][2][3] = {
    {
        // Thresholds to use if block is screen content or an I-slice
        {{200, 200, 200}, {200, 200, 200}},
        {{200, 200, 200}, {200, 200, 200}},
        {{200, 200, 200}, {200, 200, 200}},
        {{200, 200, 200}, {200, 200, 200}},
        {{200, 200, 200}, {200, 200, 200}},
        {{200, 200, 200}, {200, 200, 200}},
        {{200, 200, 200}, {200, 200, 200}},
        {{200, 200, 200}, {200, 200, 200}},
        {{200, 200, 200}, {200, 200, 200}},
    },
    {
        // Thresholds to use if block is not screen content or an I-slice
        {{100, 10, 10}, {100, 10, 10}},
        {{100, 10, 10}, {100, 10, 10}},
        {{100, 10, 10}, {100, 10, 10}},
        {{100, 10, 10}, {100, 10, 10}},
        {{100, 10, 10}, {100, 10, 10}},
        {{100, 10, 10}, {100, 10, 10}},
        {{100, 10, 10}, {100, 10, 10}},
        {{100, 10, 10}, {100, 10, 10}},
        {{100, 10, 10}, {100, 10, 10}},
    }};
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
}

#if CS2_ADOPTIONS_1
static uint64_t generate_best_part_cost(
    SequenceControlSet  *scs_ptr,
    PictureControlSet   *pcs_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {
    uint32_t  blk_index = 0;
    uint64_t best_part_cost = 0;
    EbBool split_flag;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        // if the parent square is inside inject this block
        uint8_t is_blk_allowed = pcs_ptr->slice_type != I_SLICE ? 1 :
            (blk_geom->sq_size < 128) ? 1 : 0;
        // derive split_flag
        split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;
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
#endif
static void perform_pred_depth_refinement(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                          ModeDecisionContext *context_ptr, uint32_t sb_index) {
    MdcSbData *results_ptr = &pcs_ptr->mdc_sb_array[sb_index];
    uint32_t   blk_index   = 0;

    // Reset mdc_sb_array data to defaults; it will be updated based on the predicted blocks (stored in md_blk_arr_nsq)
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom                              = get_blk_geom_mds(blk_index);
        results_ptr->leaf_data_array[blk_index].consider_block = 0;
        results_ptr->leaf_data_array[blk_index].split_flag =
            blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        results_ptr->leaf_data_array[blk_index].refined_split_flag =
            blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        blk_index++;
    }

    results_ptr->leaf_count = 0;
    blk_index               = 0;

    SuperBlock *sb_ptr = pcs_ptr->sb_ptr_array[sb_index];

    uint32_t tot_d1_blocks, block_1d_idx;
    EbBool   split_flag;

    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
        tot_d1_blocks             = blk_geom->sq_size == 128
                            ? 17
                            : blk_geom->sq_size > 8 ? 25 : blk_geom->sq_size == 8 ? 5 : 1;

        // if the parent square is inside inject this block
        uint8_t is_blk_allowed =
            pcs_ptr->slice_type != I_SLICE ? 1 : (blk_geom->sq_size < 128) ? 1 : 0;

        // derive split_flag
        split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;

        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
                    int8_t s_depth = 0;
                    int8_t e_depth = 0;

                    if (context_ptr->pd_pass == PD_PASS_0) {
#if CS2_ADOPTIONS_1
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

                        if (best_part_cost < early_exit_th) {
                            s_depth = 0;
                            e_depth = 0;
                        }
                        else {
#endif
                        derive_start_end_depth(pcs_ptr,
                                               sb_ptr,
                                               scs_ptr->seq_header.sb_size,
                                               &s_depth,
                                               &e_depth,
                                               blk_geom);
#if CS2_ADOPTIONS_1
                        }
#endif
                    } else if (context_ptr->pd_pass == PD_PASS_1) {
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
                            s_depth = 0;
                            e_depth = 0;
                        } else

                            if (context_ptr->md_local_blk_unit[blk_index].best_d1_blk == blk_index) {
                                s_depth = -1;
                                e_depth = 0;
                            } else {
                                s_depth = 0;
                                e_depth = 1;
                            }
                    }

                    // Add current pred depth block(s)
                    for (block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].consider_block = 1;
                        results_ptr->leaf_data_array[blk_index + block_1d_idx].refined_split_flag =
                            EB_FALSE;
                    }

                    // Add block indices of upper depth(s)
                    if (s_depth != 0)
                        set_parent_to_be_considered(
                            results_ptr, blk_index, scs_ptr->seq_header.sb_size, s_depth);

                    // Add block indices of lower depth(s)
                    if (e_depth != 0)
                        set_child_to_be_considered(
                            results_ptr, blk_index, scs_ptr->seq_header.sb_size, e_depth);
                }
            }
        }
        blk_index +=
            split_flag
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
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
    uint8_t     sb_sz;
    uint8_t     sb_size_log2;
    uint32_t    x_sb_index;
    uint32_t    y_sb_index;
    uint32_t    sb_origin_x;
    uint32_t    sb_origin_y;
    EbBool      last_sb_flag;
    EbBool      end_of_row_flag;
    uint32_t    sb_row_index_start;
    uint32_t    sb_row_index_count;
    uint32_t    pic_width_in_sb;
    MdcSbData * mdc_ptr;

    // Variables
    EbBool is_16bit;

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

#if TILES_PARALLEL
    uint16_t tile_group_width_in_sb;
#endif

    segment_index = 0;

    for (;;) {
        // Get Mode Decision Results
        EB_GET_FULL_OBJECT(context_ptr->mode_decision_input_fifo_ptr, &enc_dec_tasks_wrapper_ptr);

        enc_dec_tasks_ptr = (EncDecTasks *)enc_dec_tasks_wrapper_ptr->object_ptr;
        pcs_ptr           = (PictureControlSet *)enc_dec_tasks_ptr->pcs_wrapper_ptr->object_ptr;
        scs_ptr           = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
#if TILES_PARALLEL
        context_ptr->tile_group_index = enc_dec_tasks_ptr->tile_group_index;
        context_ptr->coded_sb_count   = 0;
        segments_ptr = pcs_ptr->enc_dec_segment_ctrl[context_ptr->tile_group_index];
#else
        segments_ptr                     = pcs_ptr->enc_dec_segment_ctrl;
#endif
        last_sb_flag = EB_FALSE;
        is_16bit     = (EbBool)(scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
        (void)is_16bit;
        (void)end_of_row_flag;
        // SB Constants
        sb_sz              = (uint8_t)scs_ptr->sb_size_pix;
        sb_size_log2       = (uint8_t)Log2f(sb_sz);
        context_ptr->sb_sz = sb_sz;
        pic_width_in_sb    = (pcs_ptr->parent_pcs_ptr->aligned_width + sb_sz - 1) >> sb_size_log2;
#if TILES_PARALLEL
        tile_group_width_in_sb =
            pcs_ptr->parent_pcs_ptr->tile_group_info[context_ptr->tile_group_index]
                .tile_group_width_in_sb;
#endif
        end_of_row_flag    = EB_FALSE;
        sb_row_index_start = sb_row_index_count = 0;
        context_ptr->tot_intra_coded_area       = 0;

        // Segment-loop
        while (assign_enc_dec_segments(segments_ptr,
                                       &segment_index,
                                       enc_dec_tasks_ptr,
                                       context_ptr->enc_dec_feedback_fifo_ptr) == EB_TRUE) {
            x_sb_start_index = segments_ptr->x_start_array[segment_index];
            y_sb_start_index = segments_ptr->y_start_array[segment_index];
#if TILES_PARALLEL
            sb_start_index = y_sb_start_index * tile_group_width_in_sb + x_sb_start_index;
#else
            sb_start_index = y_sb_start_index * pic_width_in_sb + x_sb_start_index;
#endif
            sb_segment_count = segments_ptr->valid_sb_count_array[segment_index];

            segment_row_index = segment_index / segments_ptr->segment_band_count;
            segment_band_index =
                segment_index - segment_row_index * segments_ptr->segment_band_count;
            segment_band_size = (segments_ptr->sb_band_count * (segment_band_index + 1) +
                                 segments_ptr->segment_band_count - 1) /
                                segments_ptr->segment_band_count;

            // Reset Coding Loop State
#if TILES_PARALLEL
            reset_mode_decision(scs_ptr,
                                context_ptr->md_context,
                                pcs_ptr,
                                context_ptr->tile_group_index,
                                segment_index);
#else
            reset_mode_decision(scs_ptr, context_ptr->md_context, pcs_ptr, segment_index);
#endif

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
#if TILES_PARALLEL
                     x_sb_index < tile_group_width_in_sb &&
#else
                     x_sb_index < pic_width_in_sb &&
#endif
                     (x_sb_index + y_sb_index < segment_band_size) &&
                     sb_segment_index < sb_start_index + sb_segment_count;
                     ++x_sb_index, ++sb_segment_index) {
#if TILES_PARALLEL
                    uint16_t tile_group_y_sb_start =
                        pcs_ptr->parent_pcs_ptr->tile_group_info[context_ptr->tile_group_index]
                            .tile_group_sb_start_y;
                    uint16_t tile_group_x_sb_start =
                        pcs_ptr->parent_pcs_ptr->tile_group_info[context_ptr->tile_group_index]
                            .tile_group_sb_start_x;
                    sb_index = (uint16_t)((y_sb_index + tile_group_y_sb_start) * pic_width_in_sb +
                                          x_sb_index + tile_group_x_sb_start);
                    sb_ptr   = pcs_ptr->sb_ptr_array[sb_index];
                    sb_origin_x = (x_sb_index + tile_group_x_sb_start) << sb_size_log2;
                    sb_origin_y = (y_sb_index + tile_group_y_sb_start) << sb_size_log2;
                    //printf("[%ld]:ED sb index %d, (%d, %d), encoded total sb count %d, ctx coded sb count %d\n",
                    //        pcs_ptr->picture_number,
                    //        sb_index, sb_origin_x, sb_origin_y,
                    //        pcs_ptr->enc_dec_coded_sb_count,
                    //        context_ptr->coded_sb_count);
                    context_ptr->tile_index             = sb_ptr->tile_info.tile_rs_index;
                    context_ptr->md_context->tile_index = sb_ptr->tile_info.tile_rs_index;

                    end_of_row_flag =
                        (x_sb_index + 1 == tile_group_width_in_sb) ? EB_TRUE : EB_FALSE;
                    sb_row_index_start =
                        (x_sb_index + 1 == tile_group_width_in_sb && sb_row_index_count == 0)
                            ? y_sb_index
                            : sb_row_index_start;
                    sb_row_index_count = (x_sb_index + 1 == tile_group_width_in_sb)
                                             ? sb_row_index_count + 1
                                             : sb_row_index_count;
#else
                    sb_index        = (uint16_t)(y_sb_index * pic_width_in_sb + x_sb_index);
                    sb_ptr          = pcs_ptr->sb_ptr_array[sb_index];
                    sb_origin_x     = x_sb_index << sb_size_log2;
                    sb_origin_y     = y_sb_index << sb_size_log2;
                    last_sb_flag    = (sb_index == pcs_ptr->sb_total_count_pix - 1) ? EB_TRUE : EB_FALSE;
                    end_of_row_flag = (x_sb_index == pic_width_in_sb - 1) ? EB_TRUE : EB_FALSE;
                    sb_row_index_start =
                        (x_sb_index == pic_width_in_sb - 1 && sb_row_index_count == 0)
                            ? y_sb_index
                            : sb_row_index_start;
                    sb_row_index_count = (x_sb_index == pic_width_in_sb - 1)
                                             ? sb_row_index_count + 1
                                             : sb_row_index_count;
#endif
                    mdc_ptr               = &pcs_ptr->mdc_sb_array[sb_index];
                    context_ptr->sb_index = sb_index;

                    if (pcs_ptr->update_cdf) {
#if MD_RATE_EST_ENH
                        if (scs_ptr->seq_header.pic_based_rate_est &&
                            scs_ptr->enc_dec_segment_row_count_array[pcs_ptr->temporal_layer_index] == 1 &&
                            scs_ptr->enc_dec_segment_col_count_array[pcs_ptr->temporal_layer_index] == 1) {
                            if (sb_index == 0)
                                pcs_ptr->ec_ctx_array[sb_index] = *pcs_ptr->coeff_est_entropy_coder_ptr->fc;
                            else
                                pcs_ptr->ec_ctx_array[sb_index] = pcs_ptr->ec_ctx_array[sb_index - 1];
                        }
                        else {
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
                            }
                        }
#else
                        // Use the latest available CDF for the current SB
                        // Use the weighted average of left (3x) and top (1x) if available.
                        int8_t up_available   = ((int32_t)(sb_origin_y >> MI_SIZE_LOG2) >
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
                        }
#endif

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
                        context_ptr->md_context, pcs_ptr, (uint8_t)sb_ptr->qp);
                    // Multi-Pass PD Path
                    // For each SB, all blocks are tested in PD0 (4421 blocks if 128x128 SB, and 1101 blocks if 64x64 SB).
                    // Then the PD0 predicted Partitioning Structure is refined by considering up to three refinements depths away from the predicted depth, both in the direction of smaller block sizes and in the direction of larger block sizes (up to Pred - 3 / Pred + 3 refinement). The selection of the refinement depth is performed using the cost
                    // deviation between the current depth cost and candidate depth cost. The generated blocks are used as input candidates to PD1.
                    // The PD1 predicted Partitioning Structure is also refined (up to Pred - 1 / Pred + 1 refinement) using the square (SQ) vs. non-square (NSQ) decision(s)
                    // inside the predicted depth and using coefficient information. The final set of blocks is evaluated in PD2 to output the final Partitioning Structure

                    if ((pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_0 ||
                         pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_1 ||
                         pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2 ||
                         pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3) &&
                        pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].is_complete_sb) {
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

                        // Perform Pred_0 depth refinement - Add blocks to be considered in the next stage(s) of PD based on depth cost.
                        perform_pred_depth_refinement(
                            scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

                        // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                        build_cand_block_array(scs_ptr, pcs_ptr, sb_index);

                        // Reset neighnor information to current SB @ position (0,0)
                        copy_neighbour_arrays(pcs_ptr,
                                              context_ptr->md_context,
                                              MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                              MD_NEIGHBOR_ARRAY_INDEX,
                                              0,
                                              sb_origin_x,
                                              sb_origin_y);

                        if (pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_1 ||
                            pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2 ||
                            pcs_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3) {
                            // [PD_PASS_1] Signal(s) derivation
                            context_ptr->md_context->pd_pass = PD_PASS_1;
                            signal_derivation_enc_dec_kernel_oq(
                                scs_ptr, pcs_ptr, context_ptr->md_context);

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
                            build_cand_block_array(scs_ptr, pcs_ptr, sb_index);

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
                    signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);

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

                    // Configure the SB
                    enc_dec_configure_sb(context_ptr, sb_ptr, pcs_ptr, (uint8_t)sb_ptr->qp);

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

#if TILES_PARALLEL
                    context_ptr->coded_sb_count++;
#endif
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
#if TILES_PARALLEL
        pcs_ptr->enc_dec_coded_sb_count += (uint32_t)context_ptr->coded_sb_count;
        last_sb_flag = (pcs_ptr->sb_total_count_pix == pcs_ptr->enc_dec_coded_sb_count);
#endif
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
            pcs_ptr->parent_pcs_ptr->av1x->rdmult = context_ptr->full_lambda;
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
