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
#include "EbPictureAnalysisProcess.h"
#if FTR_MEM_OPT
 void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr, EbBool is_highbd);
#endif
#if OPT_MEM_PALETTE
 int svt_av1_allow_palette(int allow_palette, BlockSize sb_type);
#endif
#define FC_SKIP_TX_SR_TH025 125 // Fast cost skip tx search threshold.
#define FC_SKIP_TX_SR_TH010 110 // Fast cost skip tx search threshold.
void copy_mv_rate(PictureControlSet *pcs, MdRateEstimationContext * dst_rate);
void svt_av1_cdef_search(EncDecContext *context_ptr, SequenceControlSet *scs_ptr,
                         PictureControlSet *pcs_ptr);

void av1_cdef_frame16bit(uint8_t is_16bit, SequenceControlSet *scs_ptr, PictureControlSet *pCs);

void svt_av1_add_film_grain(EbPictureBufferDesc *src, EbPictureBufferDesc *dst,
                            AomFilmGrain *film_grain_ptr);

void svt_av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm,
                                                  int32_t after_cdef);
void svt_av1_loop_restoration_filter_frame(Yv12BufferConfig *frame, Av1Common *cm,
                                           int32_t optimized_lr);
extern void get_recon_pic(PictureControlSet* pcs_ptr, EbPictureBufferDesc** recon_ptr, EbBool is_highbd);
#if SS_2B_COMPRESS
void svt_c_unpack_compressed_10bit(const uint8_t *inn_bit_buffer, uint32_t inn_stride, uint8_t *in_compn_bit_buffer, uint32_t out_stride, uint32_t height);
#endif

#if TUNE_BYPASS_MEM

/*
* return by-pass encdec
*/
uint8_t  get_bypass_encdec( EbEncMode enc_mode ,uint8_t hbd_mode_decision,uint8_t encoder_bit_depth ) {
#if CLN_FEAT_LEVEL
    UNUSED(hbd_mode_decision);
#endif
    uint8_t  bypass_encdec  = 1 ;
#if FTR_10BIT_MDS3_LPD1
        if (encoder_bit_depth == EB_8BIT) {
            // 8bit settings
#if TUNE_ED_BYPASS_8BIT
            if (enc_mode <= ENC_M6)
#else
            if (enc_mode <= ENC_M4)
#endif
                bypass_encdec = 0;
            else
                bypass_encdec = 1;
        }
        else {
            // 10bit settings
#if TUNE_M7_M8_3
#if CLN_FEAT_LEVEL
            if (enc_mode <= ENC_M9)
#else
            if (enc_mode <= ENC_M6)
#endif
#else
            if (enc_mode <= ENC_M7)
#endif
                bypass_encdec = 0;
#if FIX_10BIT_R2R
#if TUNE_M8_M10_4K_SUPER
#if TUNE_M11_SLOWDOWN
#if CLN_REMOVE_UNUSED_FEATS
#if !CLN_FEAT_LEVEL
            else if (enc_mode <= ENC_M10)
#endif
#else
            else if (enc_mode <= ENC_M11)
#endif
#else
            // To fix max BDR issues with clips like cosmos_aom_sdr_12149-12330_588x250
            else if (enc_mode <= ENC_M10)
#endif
#else
            else if (enc_mode <= ENC_M9)
#endif
#else
            if (enc_mode <= ENC_M9)
#endif
#if !CLN_FEAT_LEVEL
                bypass_encdec = hbd_mode_decision ? 1 : 0;
#endif
            else
                bypass_encdec = 1;
        }
#else
#if TUNE_MEGA_M9_M4
#if TUNE_M10_M3_1
        if (enc_mode <= ENC_M4)
#else
        if (enc_mode <= ENC_M5)
#endif
#else
        if (enc_mode <= ENC_M8)
#endif
            bypass_encdec = 0;
        else
            bypass_encdec = 1;
#endif
   return bypass_encdec ;
}



#endif
static void enc_dec_context_dctor(EbPtr p) {
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    EncDecContext *  obj                = (EncDecContext *)thread_context_ptr->priv;
    EB_DELETE(obj->md_context);
    EB_DELETE(obj->residual_buffer);
    EB_DELETE(obj->transform_buffer);
    EB_DELETE(obj->inverse_quant_buffer);
    EB_DELETE(obj->input_sample16bit_buffer);
#if !CLN_MDCONTEXT
    if (obj->is_md_rate_estimation_ptr_owner) EB_FREE(obj->md_rate_estimation_ptr);
#endif
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
    EbColorFormat color_format             = static_config->encoder_color_format;
    int8_t       enable_hbd_mode_decision = static_config->enable_hbd_mode_decision;

    EncDecContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = enc_dec_context_dctor;

    context_ptr->is_16bit     = static_config->is_16bit_pipeline;
    context_ptr->color_format = color_format;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->enc_dec_tasks_resource_ptr, index);
    context_ptr->enc_dec_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_results_resource_ptr, index);
    context_ptr->enc_dec_feedback_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->enc_dec_tasks_resource_ptr, tasks_index);
    context_ptr->picture_demux_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);

#if !CLN_MDCONTEXT
    // MD rate Estimation tables
    EB_MALLOC(context_ptr->md_rate_estimation_ptr, sizeof(MdRateEstimationContext));
    context_ptr->is_md_rate_estimation_ptr_owner = EB_TRUE;
#endif

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
#if TUNE_BYPASS_MEM

    EB_NEW(context_ptr->md_context,
        mode_decision_context_ctor,
        color_format,
        static_config->super_block_size,
        static_config->enc_mode,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr->max_block_cnt,
        static_config->encoder_bit_depth,
        0,
        0,
        enable_hbd_mode_decision == DEFAULT ? 2 : enable_hbd_mode_decision,
        static_config->screen_content_mode);
#elif CLN_GEOM

    EB_NEW(context_ptr->md_context,
        mode_decision_context_ctor,
        color_format,
        static_config->super_block_size,
        static_config->enc_mode,
        enc_handle_ptr->scs_instance_array[0]->scs_ptr->max_block_cnt,
        0,
        0,
        enable_hbd_mode_decision == DEFAULT ? 2 : enable_hbd_mode_decision,
        static_config->screen_content_mode);
#else
    EB_NEW(context_ptr->md_context,
           mode_decision_context_ctor,
           color_format,
           static_config->super_block_size,
           static_config->enc_mode,
           0,
           0,
           enable_hbd_mode_decision == DEFAULT ? 2 : enable_hbd_mode_decision ,
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
    neighbor_array_unit_reset(pcs_ptr->ep_luma_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_recon_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_luma_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_dc_sign_level_coeff_neighbor_array[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_dc_sign_level_coeff_neighbor_array[tile_idx]);
#if REFCTR_SEP_ENCDEC
    neighbor_array_unit_reset(pcs_ptr->ep_luma_dc_sign_level_coeff_neighbor_array_update[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cb_dc_sign_level_coeff_neighbor_array_update[tile_idx]);
    neighbor_array_unit_reset(pcs_ptr->ep_cr_dc_sign_level_coeff_neighbor_array_update[tile_idx]);
#endif
    neighbor_array_unit_reset(pcs_ptr->ep_partition_context_neighbor_array[tile_idx]);
    // TODO(Joel): 8-bit ep_luma_recon_neighbor_array (Cb,Cr) when is_16bit==0?
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline) {
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
    context_ptr->is_16bit = scs_ptr->static_config.is_16bit_pipeline;
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
#if !CLN_MDCONTEXT
    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    if (context_ptr->is_md_rate_estimation_ptr_owner) {
        EB_FREE(context_ptr->md_rate_estimation_ptr);
        context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
    }
    context_ptr->md_rate_estimation_ptr = pcs_ptr->md_rate_estimation_array;
#endif
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
            EncDecTasks *    feedback_task_ptr     = (EncDecTasks *)wrapper_ptr->object_ptr;
            feedback_task_ptr->input_type          = ENCDEC_TASKS_ENCDEC_INPUT;
            feedback_task_ptr->enc_dec_segment_row = feedback_row_index;
            feedback_task_ptr->pcs_wrapper_ptr     = taskPtr->pcs_wrapper_ptr;
            feedback_task_ptr->tile_group_index = taskPtr->tile_group_index;
            svt_post_full_object(wrapper_ptr);
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
    svt_block_on_mutex(encode_context_ptr->total_number_of_recon_frame_mutex);

    if (!pcs_ptr->parent_pcs_ptr->is_alt_ref) {
        EbBool           is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
        EbObjectWrapper *output_recon_wrapper_ptr;
        // Get Recon Buffer
        svt_get_empty_object(scs_ptr->encode_context_ptr->recon_output_fifo_ptr,
                             &output_recon_wrapper_ptr);
        EbBufferHeaderType *output_recon_ptr = (EbBufferHeaderType *)output_recon_wrapper_ptr->object_ptr;
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
            {
#if FTR_MEM_OPT

                get_recon_pic(pcs_ptr, &recon_ptr, is_16bit);
#else
                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    recon_ptr = is_16bit ? ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr
                                                ->reference_picture_wrapper_ptr->object_ptr)
                                               ->reference_picture16bit
                                         : ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr
                                                ->reference_picture_wrapper_ptr->object_ptr)
                                               ->reference_picture;
                else {
                    if (is_16bit)
                        recon_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr;
                    else
                        recon_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr;
                }
#endif
            }

            // FGN: Create a buffer if needed, copy the reconstructed picture and run the film grain synthesis algorithm

            if (scs_ptr->seq_header.film_grain_params_present
                && pcs_ptr->parent_pcs_ptr->frm_hdr.film_grain_params.apply_grain) {
                EbPictureBufferDesc *intermediate_buffer_ptr;
                {   // temp_lf_recon_picture_ptr is finished in use
                    // and reuse it to copy recon and add film grain noise
                    if (is_16bit)
                        intermediate_buffer_ptr = pcs_ptr->temp_lf_recon_picture16bit_ptr;
                    else
                        intermediate_buffer_ptr = pcs_ptr->temp_lf_recon_picture_ptr;
                }

                AomFilmGrain *film_grain_ptr;

                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    film_grain_ptr =
                        &((EbReferenceObject *)
                              pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                             ->film_grain_params;
                else
                    film_grain_ptr = &pcs_ptr->parent_pcs_ptr->frm_hdr.film_grain_params;

                svt_av1_add_film_grain(recon_ptr, intermediate_buffer_ptr, film_grain_ptr);
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

static void svt_aom_highbd_ssim_parms_8x8_c(const uint8_t *s, int sp, const uint8_t *sinc, int spinc, const uint16_t *r,
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
  svt_aom_ssim_parms_8x8_c(s, sp, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
  return similarity(sum_s, sum_r, sum_sq_s, sum_sq_r, sum_sxr, 64, 8);
}

static double highbd_ssim_8x8(const uint8_t *s, int sp, const uint8_t *sinc, int spinc, const uint16_t *r,
                              int rp, uint32_t bd, uint32_t shift) {
  uint32_t sum_s = 0, sum_r = 0, sum_sq_s = 0, sum_sq_r = 0, sum_sxr = 0;
  svt_aom_highbd_ssim_parms_8x8_c(s, sp, sinc, spinc, r, rp, &sum_s, &sum_r, &sum_sq_s, &sum_sq_r, &sum_sxr);
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

void free_temporal_filtering_buffer(PictureControlSet* pcs_ptr, SequenceControlSet* scs_ptr) {
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

#if SS_2B_COMPRESS
EbErrorType ssim_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory) {
#else
void ssim_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory) {
#endif
    EbBool is_16bit = (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);

    const uint32_t ss_x = scs_ptr->subsampling_x;
    const uint32_t ss_y = scs_ptr->subsampling_y;


    EbPictureBufferDesc* recon_ptr;
    get_recon_pic(pcs_ptr, &recon_ptr, is_16bit);
    if (!is_16bit) {
        EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc*)pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;

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
            assert(pcs_ptr->parent_pcs_ptr->save_source_picture_width == input_picture_ptr->width
                && pcs_ptr->parent_pcs_ptr->save_source_picture_height == input_picture_ptr->height);
            buffer_y = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[0];
            buffer_cb = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[1];
            buffer_cr = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[2];
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

#if FTR_MEM_OPT

        get_recon_pic(pcs_ptr, &recon_ptr, is_16bit);
#else
        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_ptr = ((EbReferenceObject*)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            recon_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr;
#endif
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
                assert(pcs_ptr->parent_pcs_ptr->save_source_picture_width == input_picture_ptr->width
                    && pcs_ptr->parent_pcs_ptr->save_source_picture_height == input_picture_ptr->height);
                buffer_y = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[0];
                buffer_bit_inc_y = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[0];
                buffer_cb = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[1];
                buffer_bit_inc_cb = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[1];
                buffer_cr = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[2];
                buffer_bit_inc_cr = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[2];
            }else{
#if SS_2B_COMPRESS
                uint32_t height_y = (uint32_t)(
                    input_picture_ptr->height +
                    input_picture_ptr->origin_y +
                    input_picture_ptr->origin_bot_y);
                uint32_t height_uv = (uint32_t)(
                    (input_picture_ptr->height +
                        input_picture_ptr->origin_y +
                        input_picture_ptr->origin_bot_y) >>
                    ss_y);

                uint8_t *uncompressed_pics[3];
                EB_MALLOC_ARRAY(uncompressed_pics[0], pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->luma_size);
                EB_MALLOC_ARRAY(uncompressed_pics[1], pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->chroma_size);
                EB_MALLOC_ARRAY(uncompressed_pics[2], pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->chroma_size);

                svt_c_unpack_compressed_10bit(
                    input_picture_ptr->buffer_bit_inc_y,
                    input_picture_ptr->stride_bit_inc_y / 4,
                    uncompressed_pics[0],
                    input_picture_ptr->stride_bit_inc_y,
                    height_y);
                // U
                svt_c_unpack_compressed_10bit(
                    input_picture_ptr->buffer_bit_inc_cb,
                    input_picture_ptr->stride_bit_inc_cb / 4,
                    uncompressed_pics[1],
                    input_picture_ptr->stride_bit_inc_cb,
                    height_uv);
                // V
                svt_c_unpack_compressed_10bit(
                    input_picture_ptr->buffer_bit_inc_cr,
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
#else
                buffer_y = input_picture_ptr->buffer_y;
                buffer_bit_inc_y = input_picture_ptr->buffer_bit_inc_y;
                buffer_cb = input_picture_ptr->buffer_cb;
                buffer_bit_inc_cb = input_picture_ptr->buffer_bit_inc_cb;
                buffer_cr = input_picture_ptr->buffer_cr;
                buffer_bit_inc_cr = input_picture_ptr->buffer_bit_inc_cr;
#endif
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
#if SS_2B_COMPRESS
            if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_FALSE) {
                EB_FREE_ARRAY(buffer_bit_inc_y);
                EB_FREE_ARRAY(buffer_bit_inc_cb);
                EB_FREE_ARRAY(buffer_bit_inc_cr);
            }
        }
    }
    return EB_ErrorNone;
#else
        }
    }
#endif
}

#if SS_2B_COMPRESS
EbErrorType psnr_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory) {
#else
void psnr_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory) {
#endif
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
            assert(pcs_ptr->parent_pcs_ptr->save_source_picture_width == input_picture_ptr->width
                && pcs_ptr->parent_pcs_ptr->save_source_picture_height == input_picture_ptr->height);
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
        pcs_ptr->parent_pcs_ptr->luma_sse = sse_total[0];
        pcs_ptr->parent_pcs_ptr->cb_sse   = sse_total[1];
        pcs_ptr->parent_pcs_ptr->cr_sse   = sse_total[2];

        if(free_memory && pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            EB_FREE_ARRAY(buffer_y);
            EB_FREE_ARRAY(buffer_cb);
            EB_FREE_ARRAY(buffer_cr);
        }
    }
    else {
        EbPictureBufferDesc *recon_ptr;


#if FTR_MEM_OPT

        get_recon_pic(pcs_ptr, &recon_ptr, is_16bit);
#else
        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_ptr = ((EbReferenceObject *)
                             pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                            ->reference_picture16bit;
        else
            recon_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr;
#endif
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
                assert(pcs_ptr->parent_pcs_ptr->save_source_picture_width == input_picture_ptr->width
                    && pcs_ptr->parent_pcs_ptr->save_source_picture_height == input_picture_ptr->height);
                buffer_y          = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[0];
                buffer_bit_inc_y  = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[0];
                buffer_cb         = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[1];
                buffer_bit_inc_cb = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[1];
                buffer_cr         = pcs_ptr->parent_pcs_ptr->save_source_picture_ptr[2];
                buffer_bit_inc_cr = pcs_ptr->parent_pcs_ptr->save_source_picture_bit_inc_ptr[2];
            } else {
#if SS_2B_COMPRESS
                uint32_t height_y = (uint32_t)(
                    input_picture_ptr->height +
                    input_picture_ptr->origin_y +
                    input_picture_ptr->origin_bot_y);
                uint32_t height_uv = (uint32_t)(
                    (input_picture_ptr->height +
                        input_picture_ptr->origin_y +
                        input_picture_ptr->origin_bot_y) >>
                    ss_y);

                uint8_t *uncompressed_pics[3];
                EB_MALLOC_ARRAY(uncompressed_pics[0], pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->luma_size);
                EB_MALLOC_ARRAY(uncompressed_pics[1], pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->chroma_size);
                EB_MALLOC_ARRAY(uncompressed_pics[2], pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr->chroma_size);

                svt_c_unpack_compressed_10bit(
                    input_picture_ptr->buffer_bit_inc_y,
                    input_picture_ptr->stride_bit_inc_y / 4,
                    uncompressed_pics[0],
                    input_picture_ptr->stride_bit_inc_y,
                    height_y);
                // U
                svt_c_unpack_compressed_10bit(
                    input_picture_ptr->buffer_bit_inc_cb,
                    input_picture_ptr->stride_bit_inc_cb / 4,
                    uncompressed_pics[1],
                    input_picture_ptr->stride_bit_inc_cb,
                    height_uv);
                // V
                svt_c_unpack_compressed_10bit(
                    input_picture_ptr->buffer_bit_inc_cr,
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
#else
                buffer_y          = input_picture_ptr->buffer_y;
                buffer_bit_inc_y  = input_picture_ptr->buffer_bit_inc_y;
                buffer_cb         = input_picture_ptr->buffer_cb;
                buffer_bit_inc_cb = input_picture_ptr->buffer_bit_inc_cb;
                buffer_cr         = input_picture_ptr->buffer_cr;
                buffer_bit_inc_cr = input_picture_ptr->buffer_bit_inc_cr;
#endif
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
#if SS_2B_COMPRESS
            if (pcs_ptr->parent_pcs_ptr->temporal_filtering_on == EB_FALSE) {
                EB_FREE_ARRAY(buffer_bit_inc_y);
                EB_FREE_ARRAY(buffer_bit_inc_cb);
                EB_FREE_ARRAY(buffer_bit_inc_cr);
            }
        }
        pcs_ptr->parent_pcs_ptr->luma_sse = (uint32_t)sse_total[0];
        pcs_ptr->parent_pcs_ptr->cb_sse   = (uint32_t)sse_total[1];
        pcs_ptr->parent_pcs_ptr->cr_sse   = (uint32_t)sse_total[2];
    }
    return EB_ErrorNone;
#else
        }

        pcs_ptr->parent_pcs_ptr->luma_sse = sse_total[0];
        pcs_ptr->parent_pcs_ptr->cb_sse   = sse_total[1];
        pcs_ptr->parent_pcs_ptr->cr_sse   = sse_total[2];
    }
#endif
}

void pad_ref_and_set_flags(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    EbReferenceObject *reference_object =
        (EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;

#if FTR_MEM_OPT

#if FTR_MEM_OPT
    EbPictureBufferDesc *ref_pic_ptr ;//= (EbPictureBufferDesc *)reference_object->reference_picture;
    EbPictureBufferDesc *ref_pic_16bit_ptr;// =   (EbPictureBufferDesc *)reference_object->reference_picture16bit;

    {
        get_recon_pic(pcs_ptr, &ref_pic_ptr, 0);
        get_recon_pic(pcs_ptr, &ref_pic_16bit_ptr, 1);

    }
#else
    EbPictureBufferDesc *ref_pic_ptr  = (EbPictureBufferDesc *)reference_object->reference_picture;
    EbPictureBufferDesc *ref_pic_16bit_ptr =  (EbPictureBufferDesc *)reference_object->reference_picture16bit;

        svt_memcpy(
            (uint16_t *)(ref_pic_16bit_ptr->buffer_y),  (uint16_t *)(pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr->buffer_y), pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr->luma_size *2);

        svt_memcpy(
            (uint16_t *)(ref_pic_16bit_ptr->buffer_cb), (uint16_t *)(pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr->buffer_cb), pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr->chroma_size*2);

        svt_memcpy(
            (uint16_t *)(ref_pic_16bit_ptr->buffer_cr),  (uint16_t *)(pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr->buffer_cr), pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr->chroma_size*2);


#endif
#else
    EbPictureBufferDesc *ref_pic_ptr = (EbPictureBufferDesc *)reference_object->reference_picture;
    EbPictureBufferDesc *ref_pic_16bit_ptr =
        (EbPictureBufferDesc *)reference_object->reference_picture16bit;
#endif
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
#if !FTR_MEM_OPT
        if (pcs_ptr->hbd_mode_decision != EB_10_BIT_MD) {
#endif
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
#if !FTR_MEM_OPT
        }
#endif
    }
    if ((scs_ptr->static_config.is_16bit_pipeline) && (!is_16bit)) {
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
#if  FTR_INTRA_DETECTOR
    pcs_ptr->intra_coded_area =
        (100 * pcs_ptr->intra_coded_area) /
        (pcs_ptr->parent_pcs_ptr->aligned_width * pcs_ptr->parent_pcs_ptr->aligned_height);
#endif
#if FTR_COEFF_DETECTOR
    pcs_ptr->skip_coded_area =
        (100 * pcs_ptr->skip_coded_area) /
        (pcs_ptr->parent_pcs_ptr->aligned_width * pcs_ptr->parent_pcs_ptr->aligned_height);
#endif
#if  FTR_INTRA_DETECTOR
    if (pcs_ptr->slice_type == I_SLICE) pcs_ptr->intra_coded_area = 0;
#endif
#if  FTR_INTRA_DETECTOR
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->intra_coded_area = (uint8_t)(pcs_ptr->intra_coded_area);
#endif
#if FTR_COEFF_DETECTOR
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->skip_coded_area = (uint8_t)(pcs_ptr->skip_coded_area);
#endif
#if MS_CDEF_OPT3
    struct PictureParentControlSet *ppcs    = pcs_ptr->parent_pcs_ptr;
    FrameHeader *                   frm_hdr = &ppcs->frm_hdr;
    ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
        ->ref_cdef_strengths_num = ppcs->nb_cdef_strengths;
    for (int i = 0; i < ppcs->nb_cdef_strengths; i++) {
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->ref_cdef_strengths[0][i] = frm_hdr->cdef_params.cdef_y_strength[i];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->ref_cdef_strengths[1][i] = frm_hdr->cdef_params.cdef_uv_strength[i];
    }
#endif
    uint32_t sb_index;
#if FTR_VLPD1
    for (sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index) {
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->non_moving_index_array[sb_index] =
            pcs_ptr->parent_pcs_ptr->non_moving_index_array[sb_index];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->sb_intra[sb_index] =
            pcs_ptr->sb_intra[sb_index];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->sb_skip[sb_index] =
            pcs_ptr->sb_skip[sb_index];
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->sb_64x64_mvp[sb_index] =
            pcs_ptr->sb_64x64_mvp[sb_index];
    }
#else
    for (sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index)
        ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
            ->non_moving_index_array[sb_index] =
            pcs_ptr->parent_pcs_ptr->non_moving_index_array[sb_index];
#endif
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
#if FTR_NEW_WN_LVLS
    // Copy the prev frame wn filter coeffs
    EbReferenceObject* obj = (EbReferenceObject*)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
    if (cm->wn_filter_ctrls.enabled && cm->wn_filter_ctrls.use_prev_frame_coeffs) {
        for (int32_t plane = 0; plane < MAX_MB_PLANE; ++plane) {
            int32_t ntiles = pcs_ptr->rst_info[plane].units_per_tile;
            for (int32_t u = 0; u < ntiles; ++u) {
                obj->unit_info[plane][u].restoration_type = pcs_ptr->rst_info[plane].unit_info[u].restoration_type;
                if (pcs_ptr->rst_info[plane].unit_info[u].restoration_type == RESTORE_WIENER)
                    obj->unit_info[plane][u].wiener_info = pcs_ptr->rst_info[plane].unit_info[u].wiener_info;
            }
        }
    }
#endif
}

void set_obmc_controls(ModeDecisionContext *mdctxt, uint8_t obmc_mode) {
    ObmcControls*obmc_ctrls = &mdctxt->obmc_ctrls;
    switch (obmc_mode)
    {
    case 0:
        obmc_ctrls->enabled = 0;
        break;
    case 1:
        obmc_ctrls->enabled = 1;
        obmc_ctrls->max_blk_size_16x16 = 0;
        break;
    case 2:
        obmc_ctrls->enabled = 1;
        obmc_ctrls->max_blk_size_16x16 = 1;
        break;
    default:
        assert(0);
        break;
    }
}
#if !OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS
void set_block_based_depth_refinement_controls(ModeDecisionContext* mdctxt, uint8_t block_based_depth_refinement_level) {

    DepthRefinementCtrls* depth_refinement_ctrls = &mdctxt->depth_refinement_ctrls;

    switch (block_based_depth_refinement_level)
    {
    case 0:
        depth_refinement_ctrls->enabled = 0;
        break;

    case 1:
        depth_refinement_ctrls->enabled = 1;
        depth_refinement_ctrls->parent_to_current_th = 25;
        depth_refinement_ctrls->sub_to_current_th = 25;
        depth_refinement_ctrls->use_pred_block_cost = 0;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;

    case 2:
        depth_refinement_ctrls->enabled = 1;
        depth_refinement_ctrls->parent_to_current_th = 25;
        depth_refinement_ctrls->sub_to_current_th = 25;
        depth_refinement_ctrls->use_pred_block_cost = 1;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;

    case 3:
        depth_refinement_ctrls->enabled = 1;
        depth_refinement_ctrls->parent_to_current_th = 20;
        depth_refinement_ctrls->sub_to_current_th = 20;
        depth_refinement_ctrls->use_pred_block_cost = 1;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;
    case 4:
        depth_refinement_ctrls->enabled = 1;
        depth_refinement_ctrls->parent_to_current_th = 15;
        depth_refinement_ctrls->sub_to_current_th = 15;
        depth_refinement_ctrls->use_pred_block_cost = 1;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;
    case 5:
        depth_refinement_ctrls->enabled = 1;
        depth_refinement_ctrls->parent_to_current_th = 10;
        depth_refinement_ctrls->sub_to_current_th = 10;
        depth_refinement_ctrls->use_pred_block_cost = 1;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;
    case 6:
        depth_refinement_ctrls->enabled = 1;
        depth_refinement_ctrls->parent_to_current_th = 10;
        depth_refinement_ctrls->sub_to_current_th = 10;
        depth_refinement_ctrls->use_pred_block_cost = 2;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;
    case 7:
        depth_refinement_ctrls->enabled = 1;
        depth_refinement_ctrls->parent_to_current_th = 5;
        depth_refinement_ctrls->sub_to_current_th = 5;
        depth_refinement_ctrls->use_pred_block_cost = 2;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;
    }
}
#endif
/*
 * Generate depth removal settings
 */
#if FTR_DEPTH_REMOVAL_QP
#if FTR_SIMPLIFIED_DEPTH_REMOVAL

#define LOW_8x8_DIST_VAR_TH  25000
#define HIGH_8x8_DIST_VAR_TH 50000

void set_depth_removal_level_controls(PictureControlSet *pcs_ptr, ModeDecisionContext *mdctxt, uint8_t depth_removal_level) {

    DepthRemovalCtrls *depth_removal_ctrls = &mdctxt->depth_removal_ctrls;

#if OPT_DEPTH_REMOVAL_I_SLICE
    if (pcs_ptr->slice_type == I_SLICE) {
        SbParams *sb_params = &pcs_ptr->parent_pcs_ptr->sb_params_array[mdctxt->sb_index];

        uint16_t disallow_below_16x16_variance_th = 0;
        uint16_t disallow_below_32x32_variance_th = 0;
        uint16_t disallow_below_64x64_variance_th = 0;

        switch (depth_removal_level) {

        case 0:
            depth_removal_ctrls->enabled = 0;
            break;

        case 1:
            depth_removal_ctrls->enabled = 1;
            disallow_below_16x16_variance_th = 150;
            disallow_below_32x32_variance_th =  50;
            disallow_below_64x64_variance_th =  25;
            break;

        }

#if FTR_M13
        depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 && sb_params->height % 16 == 0)
            ? (depth_removal_ctrls->disallow_below_16x16 || pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_16x16_variance_th)
            : 0;

        depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 && sb_params->height % 32 == 0)
            ? (depth_removal_ctrls->disallow_below_32x32 || pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_32x32_variance_th)
            : 0;

        depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 && sb_params->height % 64 == 0)
            ? (depth_removal_ctrls->disallow_below_64x64 || pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_64x64_variance_th)
            : 0;
#else
        depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 && sb_params->height % 16 == 0)
            ? (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_16x16_variance_th)
            : 0;

        depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 && sb_params->height % 32 == 0)
            ? (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_32x32_variance_th)
            : 0;

        depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 && sb_params->height % 64 == 0)
            ? (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_64x64_variance_th)
            : 0;
#endif
    }
    else {
#endif
    uint32_t me_8x8_cost_variance = pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[mdctxt->sb_index];

    SbParams *sb_params = &pcs_ptr->parent_pcs_ptr->sb_params_array[mdctxt->sb_index];

    // me_distortion => EB_8_BIT_MD
    uint32_t fast_lambda = mdctxt->fast_lambda_md[EB_8_BIT_MD];

    uint32_t sb_size = 64 * 64;

    uint64_t cost_th_rate = 1 << 13;

    uint64_t disallow_below_16x16_cost_th_multiplier = 0;
    uint64_t disallow_below_32x32_cost_th_multiplier = 0;
    uint64_t disallow_below_64x64_cost_th_multiplier = 0;

    int64_t dev_16x16_to_8x8_th = MAX_SIGNED_VALUE;
    int64_t dev_32x32_to_16x16_th = MAX_SIGNED_VALUE;
    int64_t dev_32x32_to_8x8_th = MAX_SIGNED_VALUE;

    int8_t qp_scale_factor = 0;

#if FTR_MOD_DEPTH_REMOVAL_LVL
    // Modulate depth_removal level for Layer0 frames based on the qp_offset band
#else
    // Modulate depth_removal level @ BASE based on the qp_offset band
#endif
#if CLN_TPL_LEVEL_7
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present) {
#else
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present && pcs_ptr->parent_pcs_ptr->tpl_ctrls.modulate_depth_removal_level) {
#endif
        int diff = mdctxt->sb_ptr->qindex - quantizer_to_qindex[pcs_ptr->parent_pcs_ptr->picture_qp];
        if (diff <= -12)
            depth_removal_level = MAX(0, (int)depth_removal_level - 4);
        else if (diff <= -6)
            depth_removal_level = MAX(0, (int)depth_removal_level - 3);
        else if (diff <= -3)
            depth_removal_level = MAX(0, (int)depth_removal_level - 2);
        else if (diff < 0)
            depth_removal_level = MAX(0, (int)depth_removal_level - 1);
#if FTR_MOD_DEPTH_REMOVAL_LVL && !CLN_TPL_LEVEL_7
        if (pcs_ptr->parent_pcs_ptr->tpl_ctrls.modulate_depth_removal_level > 1) {
            if (diff >= 12)
                depth_removal_level = depth_removal_level + 4;
            else if (diff >= 6)
                depth_removal_level = depth_removal_level + 3;
            else if (diff >= 3)
                depth_removal_level = depth_removal_level + 2;
            else if (diff > 0)
                depth_removal_level = depth_removal_level + 1;
         }
#endif
    }

    switch (depth_removal_level) {

    case 0:
        depth_removal_ctrls->enabled = 0;
        break;

    case 1:
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 0;
        disallow_below_32x32_cost_th_multiplier = 0;
        disallow_below_64x64_cost_th_multiplier = 0;
        dev_16x16_to_8x8_th = 15;
        dev_32x32_to_16x16_th = 0;
        qp_scale_factor = 1;

        break;

    case 2:
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 4;
        disallow_below_32x32_cost_th_multiplier = 0;
        disallow_below_64x64_cost_th_multiplier = 0;
        dev_16x16_to_8x8_th = 20;
        dev_32x32_to_16x16_th = 0;
        qp_scale_factor = 1;

        break;
#if TUNE_M7_M8_3
    case 3:
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 4;
        disallow_below_32x32_cost_th_multiplier = 0;
        disallow_below_64x64_cost_th_multiplier = 0;
        dev_16x16_to_8x8_th = 50;
        dev_32x32_to_16x16_th = 0;
        qp_scale_factor = 1;
        break;
#endif
#if TUNE_M7_M8_3
    case 4:
#else
    case 3:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 4;
        disallow_below_32x32_cost_th_multiplier = 0;
        disallow_below_64x64_cost_th_multiplier = 0;
        dev_16x16_to_8x8_th = 50;
        dev_32x32_to_16x16_th = 0;
        qp_scale_factor = 2;
        break;
#if TUNE_M7_M8_3
    case 5:
#else
    case 4:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 4;
        disallow_below_32x32_cost_th_multiplier = 0;
        disallow_below_64x64_cost_th_multiplier = 0;
        dev_16x16_to_8x8_th = 50;
        dev_32x32_to_16x16_th = 0;
        qp_scale_factor = 2;
        break;
#if TUNE_M7_M8_3
    case 6:
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 4;
        disallow_below_32x32_cost_th_multiplier = 2;
        disallow_below_64x64_cost_th_multiplier = 0;
        dev_16x16_to_8x8_th = 50;
        dev_32x32_to_16x16_th = 0;
        qp_scale_factor = 2;
        break;
#endif
#if TUNE_M7_M8_3
    case 7:
#else
    case 5:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 4;
        disallow_below_32x32_cost_th_multiplier = 2;
        disallow_below_64x64_cost_th_multiplier = 2;
        dev_16x16_to_8x8_th = 50;
        dev_32x32_to_16x16_th = 0;
        qp_scale_factor = 2;
        break;
#if TUNE_M7_M8_3
    case 8:
#else
    case 6:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 4;
        disallow_below_32x32_cost_th_multiplier = 2;
        disallow_below_64x64_cost_th_multiplier = 2;
        dev_16x16_to_8x8_th = 100;
        dev_32x32_to_16x16_th = 50;
        qp_scale_factor = 3;
        break;
#if TUNE_M7_M8_3
    case 9:
#else
    case 7:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 8;
        disallow_below_32x32_cost_th_multiplier = 2;
        disallow_below_64x64_cost_th_multiplier = 2;
        dev_16x16_to_8x8_th = 100;
        dev_32x32_to_16x16_th = 50;
        qp_scale_factor = 3;

        break;
#if TUNE_M7_M8_3
    case 10:
#else
    case 8:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 32;
        disallow_below_32x32_cost_th_multiplier = 2;
        disallow_below_64x64_cost_th_multiplier = 2;
        dev_16x16_to_8x8_th = 200;
        dev_32x32_to_16x16_th = 75;
        qp_scale_factor = 3;
        break;
#if TUNE_M7_M8_3
    case 11:
#else
    case 9:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 32;
        disallow_below_32x32_cost_th_multiplier = 2;
        disallow_below_64x64_cost_th_multiplier = 2;
        dev_16x16_to_8x8_th = 250;
        dev_32x32_to_16x16_th = 125;
        qp_scale_factor = 3;
        break;
#if TUNE_M7_M8_3
    case 12:
#else
    case 10:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 32;
        disallow_below_32x32_cost_th_multiplier = 4;
        disallow_below_64x64_cost_th_multiplier = 2;
        dev_16x16_to_8x8_th = 250;
        dev_32x32_to_16x16_th = 150;
        qp_scale_factor = 4;
        break;
#if TUNE_M7_M8_3
    case 13:
#else
    case 11:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 64;
        disallow_below_32x32_cost_th_multiplier = 4;
        disallow_below_64x64_cost_th_multiplier = 2;
        dev_16x16_to_8x8_th = 250;
        dev_32x32_to_16x16_th = 150;
        qp_scale_factor = 4;
        break;
#if TUNE_M11_2
#if TUNE_M7_M8_3
    case 14:
#else
    case 12:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 64;
        disallow_below_32x32_cost_th_multiplier = 4;
        disallow_below_64x64_cost_th_multiplier = 4;
        dev_16x16_to_8x8_th = 250;
        dev_32x32_to_16x16_th = 150;
        qp_scale_factor = 4;
        break;
#if TUNE_4K_M11
#if TUNE_M7_M8_3
    case 15:
#else
    case 13:
#endif
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 96;
        disallow_below_32x32_cost_th_multiplier = 6;
        disallow_below_64x64_cost_th_multiplier = 6;
        dev_16x16_to_8x8_th = 300;
        dev_32x32_to_16x16_th = 200;
        qp_scale_factor = 4;
        break;
#endif
#else
    case 12:
        depth_removal_ctrls->enabled = 1;
        disallow_below_16x16_cost_th_multiplier = 64;
        disallow_below_32x32_cost_th_multiplier = 4;
        disallow_below_64x64_cost_th_multiplier = 2;
        dev_16x16_to_8x8_th = 250;
        dev_32x32_to_16x16_th = 150;
        qp_scale_factor = 5;
        break;
#endif
    }
#if !SHUT_8x8_IF_NON_ISLICE
    depth_removal_ctrls->disallow_below_64x64 = 0;
    depth_removal_ctrls->disallow_below_32x32 = 0;
    depth_removal_ctrls->disallow_below_16x16 = 0;
#endif
    if (depth_removal_ctrls->enabled) {

        //dev_16x16_to_8x8_th , dev_32x32_to_16x16_th = f(me_8x8_cost_variance)
        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);
        if (me_8x8_cost_variance < LOW_8x8_DIST_VAR_TH) {
            dev_16x16_to_8x8_th = dev_16x16_to_8x8_th << 2;
        }
        else if (me_8x8_cost_variance < HIGH_8x8_DIST_VAR_TH) {
            dev_16x16_to_8x8_th = dev_16x16_to_8x8_th << 1;
            dev_32x32_to_16x16_th = dev_32x32_to_16x16_th >> 1;
        }
        else {
            dev_16x16_to_8x8_th = 0;
            dev_32x32_to_16x16_th = 0;
        }

        //dev_16x16_to_8x8_th , dev_32x32_to_16x16_th = f(QP)
        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * qp_scale_factor;
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * qp_scale_factor;

        // dev_32x32_to_8x8_th = f(dev_32x32_to_16x16_th); a bit higher
        dev_32x32_to_8x8_th = (dev_32x32_to_16x16_th * ((1 << 2) + 1)) >> 2;


        uint64_t disallow_below_16x16_cost_th = disallow_below_16x16_cost_th_multiplier ? RDCOST(fast_lambda, cost_th_rate, (sb_size >> 1) *  disallow_below_16x16_cost_th_multiplier) : 0;
        uint64_t disallow_below_32x32_cost_th = disallow_below_32x32_cost_th_multiplier ? RDCOST(fast_lambda, cost_th_rate, (sb_size >> 1) *  disallow_below_32x32_cost_th_multiplier) : 0;
        uint64_t disallow_below_64x64_cost_th = disallow_below_64x64_cost_th_multiplier ? RDCOST(fast_lambda, cost_th_rate, (sb_size >> 1) *  disallow_below_64x64_cost_th_multiplier) : 0;

        uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[mdctxt->sb_index]);
        uint64_t cost_32x32 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_32x32_distortion[mdctxt->sb_index]);
        uint64_t cost_16x16 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_16x16_distortion[mdctxt->sb_index]);
        uint64_t cost_8x8 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_8x8_distortion[mdctxt->sb_index]);

        int64_t dev_32x32_to_16x16 =
            (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_16x16, 1)) * 1000) /
            (int64_t)MAX(cost_16x16, 1);

        int64_t dev_32x32_to_8x8 =
            (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_8x8, 1)) * 1000) /
            (int64_t)MAX(cost_8x8, 1);

        int64_t dev_16x16_to_8x8 =
            (int64_t)(((int64_t)MAX(cost_16x16, 1) - (int64_t)MAX(cost_8x8, 1)) * 1000) /
            (int64_t)MAX(cost_8x8, 1);
#if SHUT_8x8_IF_NON_ISLICE
        depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 && sb_params->height % 64 == 0)
            ? (depth_removal_ctrls->disallow_below_64x64 || cost_64x64 < disallow_below_64x64_cost_th)
            : 0;

        depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 && sb_params->height % 32 == 0)
            ? (depth_removal_ctrls->disallow_below_32x32 || cost_32x32 < disallow_below_32x32_cost_th || (dev_32x32_to_16x16 < dev_32x32_to_16x16_th && dev_32x32_to_8x8 < dev_32x32_to_8x8_th))
            : 0;

        depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 && sb_params->height % 16 == 0)
            ? (depth_removal_ctrls->disallow_below_16x16  || cost_16x16 < disallow_below_16x16_cost_th || dev_16x16_to_8x8 < dev_16x16_to_8x8_th)
            : 0;
#else
        depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 && sb_params->height % 64 == 0)
            ? (cost_64x64 < disallow_below_64x64_cost_th)
            : 0;

        depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 && sb_params->height % 32 == 0)
            ? (cost_32x32 < disallow_below_32x32_cost_th || (dev_32x32_to_16x16 < dev_32x32_to_16x16_th && dev_32x32_to_8x8 < dev_32x32_to_8x8_th))
            : 0;

        depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 && sb_params->height % 16 == 0)
            ? (cost_16x16 < disallow_below_16x16_cost_th || dev_16x16_to_8x8 < dev_16x16_to_8x8_th)
            : 0;
#endif
    }
#if OPT_DEPTH_REMOVAL_I_SLICE
        }
#endif
}
#else
void set_depth_removal_level_controls(PictureControlSet *pcs_ptr, ModeDecisionContext *mdctxt, uint8_t block_based_depth_refinement_level) {
    DepthRemovalCtrls *depth_removal_ctrls = &mdctxt->depth_removal_ctrls;

    uint32_t me_8x8_cost_variance = pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[mdctxt->sb_index];

    SbParams *sb_params = &pcs_ptr->parent_pcs_ptr->sb_params_array[mdctxt->sb_index];
    uint32_t fast_lambda = mdctxt->hbd_mode_decision ?
        mdctxt->fast_lambda_md[EB_10_BIT_MD] :
        mdctxt->fast_lambda_md[EB_8_BIT_MD];
    uint32_t sb_size = 64 * 64;
    uint64_t cost_th_rate = 1 << 13;
    uint64_t disallow_below_64x64_th = 0;
    uint64_t disallow_below_32x32_th = 0;
    uint64_t disallow_below_16x16_th = 0;
    int64_t dev_16x16_to_8x8_th = MAX_SIGNED_VALUE;
    int64_t dev_32x32_to_16x16_th = 0;

    switch (block_based_depth_refinement_level) {
    case 0:
        depth_removal_ctrls->enabled = 0;
        break;
    case 1:
        depth_removal_ctrls->enabled = 1;
        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = 0;
        }

        dev_16x16_to_8x8_th = 2;

        if (me_8x8_cost_variance < 2000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if (me_8x8_cost_variance < 7000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if (me_8x8_cost_variance < 15000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if (me_8x8_cost_variance < 30000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 1) / 5;

        break;
    case 2:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = 0;
        }

        dev_16x16_to_8x8_th = 2;

        if (me_8x8_cost_variance < 2000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if (me_8x8_cost_variance < 7000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if (me_8x8_cost_variance < 15000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if (me_8x8_cost_variance < 30000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;

        break;
    case 3:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = 0;
        }

        dev_16x16_to_8x8_th = 2;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 2), 1);

        if (me_8x8_cost_variance < 2000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if (me_8x8_cost_variance < 7000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if (me_8x8_cost_variance < 15000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if (me_8x8_cost_variance < 30000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;

        break;
    case 4:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 16);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 12);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);

        }

        dev_16x16_to_8x8_th = 2;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 5000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if (me_8x8_cost_variance < 10000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if (me_8x8_cost_variance < 20000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if (me_8x8_cost_variance < 40000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
        break;

    case 5:
        depth_removal_ctrls->enabled = 1;
        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, (sb_size * 3) >> 1);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 16);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 4);
        }

        dev_16x16_to_8x8_th = 2;
        dev_32x32_to_16x16_th = 2;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 25000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 40000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 60000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else if (me_8x8_cost_variance < 80000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 1) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 1) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 1) / 5;
        }

        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1);
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1);

        break;
    case 6:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, (sb_size * 3) >> 1);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);
        }
        else {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, (sb_size * 3) >> 1);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 16);
        }

        dev_16x16_to_8x8_th = 2;
        dev_32x32_to_16x16_th = 2;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 25000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 70000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else if (me_8x8_cost_variance < 90000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }

        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 3), 1);
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 3), 1);

        break;
    case 7:
        depth_removal_ctrls->enabled = 1;

        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);

        dev_16x16_to_8x8_th = 3;
        dev_32x32_to_16x16_th = 2;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 35000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 70000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else if (me_8x8_cost_variance < 90000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th = 0;
        }
        else {
            dev_16x16_to_8x8_th = 0;
            dev_32x32_to_16x16_th = 0;
        }
        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;

        break;
    case 8:
        depth_removal_ctrls->enabled = 1;

        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);

        dev_16x16_to_8x8_th = 3;
        dev_32x32_to_16x16_th = 2;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 100000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 150000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 1) / 5;
        }

        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;

        break;
    case 9:
        depth_removal_ctrls->enabled = 1;

        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);

        dev_16x16_to_8x8_th = 3;
        dev_32x32_to_16x16_th = 3;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 100000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 150000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 1) / 5;
        }

        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;

        break;
    case 10:
        depth_removal_ctrls->enabled = 1;
        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);

        dev_16x16_to_8x8_th = 4;
        dev_32x32_to_16x16_th = 3;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 100000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 150000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }

        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;

        break;
    case 11:
        depth_removal_ctrls->enabled = 1;

        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);

        dev_16x16_to_8x8_th = 4;
        dev_32x32_to_16x16_th = 4;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= dev_16x16_to_8x8_th * 4;
            dev_32x32_to_16x16_th *= dev_32x32_to_16x16_th * 2;
        }
        else {
            dev_16x16_to_8x8_th = dev_16x16_to_8x8_th * 2;
        }

        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 4), 1) * 3;

        break;
    case 12:
        depth_removal_ctrls->enabled = 1;

        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 128);

        dev_16x16_to_8x8_th = 4;
        dev_32x32_to_16x16_th = 4;

        me_8x8_cost_variance /= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1)), 1);

        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= dev_16x16_to_8x8_th * 8;
            dev_32x32_to_16x16_th *= dev_32x32_to_16x16_th * 2;
        }
        else {
            dev_16x16_to_8x8_th = dev_16x16_to_8x8_th * 2;
        }

        dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 3), 1);
        dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs_ptr->picture_qp + 10), 1) >> 3), 1);

        break;
    }

    depth_removal_ctrls->disallow_below_64x64 = 0;
    depth_removal_ctrls->disallow_below_32x32 = 0;
    depth_removal_ctrls->disallow_below_16x16 = 0;

    if (depth_removal_ctrls->enabled) {
        uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[mdctxt->sb_index]);
        uint64_t cost_32x32 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_32x32_distortion[mdctxt->sb_index]);
        uint64_t cost_16x16 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_16x16_distortion[mdctxt->sb_index]);
        uint64_t cost_8x8 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_8x8_distortion[mdctxt->sb_index]);

        int64_t dev_32x32_to_16x16 =
            (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_16x16, 1)) * 100) /
            (int64_t)MAX(cost_16x16, 1);

        int64_t dev_32x32_to_8x8 =
            (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_8x8, 1)) * 100) /
            (int64_t)MAX(cost_8x8, 1);

        int64_t dev_32x32_to_8x8_th = (dev_32x32_to_16x16_th * 5) / 4;

        int64_t dev_16x16_to_8x8 =
            (int64_t)(((int64_t)MAX(cost_16x16, 1) - (int64_t)MAX(cost_8x8, 1)) * 100) /
            (int64_t)MAX(cost_8x8, 1);

        depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 && sb_params->height % 64 == 0)
            ? (cost_64x64 < disallow_below_64x64_th)
            : 0;

        depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 && sb_params->height % 32 == 0)
            ? (cost_32x32 < disallow_below_32x32_th || (dev_32x32_to_16x16 < dev_32x32_to_16x16_th && dev_32x32_to_8x8 < dev_32x32_to_8x8_th))
            : 0;

        depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 && sb_params->height % 16 == 0)
            ? (cost_16x16 < disallow_below_16x16_th || dev_16x16_to_8x8 < dev_16x16_to_8x8_th)
            : 0;
    }
}
#endif
#else
void set_depth_removal_level_controls(PictureControlSet *pcs_ptr, ModeDecisionContext *mdctxt, uint8_t block_based_depth_refinement_level) {
    DepthRemovalCtrls *depth_removal_ctrls = &mdctxt->depth_removal_ctrls;
    const uint32_t me_8x8_cost_variance = pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[mdctxt->sb_index];
    SbParams *sb_params = &pcs_ptr->parent_pcs_ptr->sb_params_array[mdctxt->sb_index];
    uint32_t fast_lambda = mdctxt->hbd_mode_decision ?
        mdctxt->fast_lambda_md[EB_10_BIT_MD] :
        mdctxt->fast_lambda_md[EB_8_BIT_MD];
    uint32_t sb_size = 64 * 64;
    uint64_t cost_th_rate = 1 << 13;
    uint64_t disallow_below_64x64_th = 0;
    uint64_t disallow_below_32x32_th = 0;
    uint64_t disallow_below_16x16_th = 0;
    int64_t dev_16x16_to_8x8_th = MAX_SIGNED_VALUE;
    int64_t dev_32x32_to_16x16_th = 0;

    switch (block_based_depth_refinement_level) {
    case 0:
        depth_removal_ctrls->enabled = 0;
        break;
    case 1:
        depth_removal_ctrls->enabled = 1;
        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = 0;
        }

        dev_16x16_to_8x8_th = 2;
        if(me_8x8_cost_variance < 2000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if(me_8x8_cost_variance <7000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if(me_8x8_cost_variance <15000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if(me_8x8_cost_variance <30000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 1) / 5;

        break;
    case 2:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = 0;
        }

        dev_16x16_to_8x8_th = 2;
        if(me_8x8_cost_variance < 2000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if(me_8x8_cost_variance <7000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if(me_8x8_cost_variance <15000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if(me_8x8_cost_variance <30000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;

        break;
    case 3:

        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 16);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 12);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        }

        dev_16x16_to_8x8_th = 2;
        if(me_8x8_cost_variance <2000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if(me_8x8_cost_variance <8000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if(me_8x8_cost_variance <13000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if(me_8x8_cost_variance <25000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
        break;
    case 4:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 16);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 12);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);

        }

        dev_16x16_to_8x8_th = 2;
        if(me_8x8_cost_variance <5000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if(me_8x8_cost_variance <10000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if(me_8x8_cost_variance <20000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if(me_8x8_cost_variance <40000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
        break;

    case 5:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = 0;
        }

        dev_16x16_to_8x8_th = 2;
        if(me_8x8_cost_variance < 3500)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
        else if(me_8x8_cost_variance <13000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        else if(me_8x8_cost_variance <30000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        else if(me_8x8_cost_variance <50000)
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
        else
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 1) / 5;
        break;
    case 6:
        depth_removal_ctrls->enabled = 1;
        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = 0;
        }
        dev_16x16_to_8x8_th = 2;
        dev_32x32_to_16x16_th = 2;
        if (me_8x8_cost_variance < 9000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 20000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 4) / 5;
        }
        else if (me_8x8_cost_variance < 70000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 3) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        break;
    case 7:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 16);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 12);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        }

        dev_16x16_to_8x8_th = 2;
        dev_32x32_to_16x16_th = 2;
        if (me_8x8_cost_variance < 9000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 20000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 70000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 3) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        break;
    case 8:
        depth_removal_ctrls->enabled = 1;

        if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 200) {
            disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
            disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 16);
        }
        else if (pcs_ptr->parent_pcs_ptr->variance[mdctxt->sb_index][ME_TIER_ZERO_PU_64x64] <= 400) {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 12);
        }
        else {
            disallow_below_64x64_th = 0;
            disallow_below_32x32_th = 0;
            disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);

        }

        dev_16x16_to_8x8_th = 2;
        dev_32x32_to_16x16_th = 2;
        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 100000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 150000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 3) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        break;
    case 9:
        depth_removal_ctrls->enabled = 1;

        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, (sb_size * 3) >> 1);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);

        dev_16x16_to_8x8_th = 2;
        dev_32x32_to_16x16_th = 2;
        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 100000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 150000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        break;
    case 10:
        depth_removal_ctrls->enabled = 1;
        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, (sb_size * 3) >> 1);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);
        dev_16x16_to_8x8_th = 4;
        dev_32x32_to_16x16_th = 3;
        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 100000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 20) / 5;
        }
        else if (me_8x8_cost_variance < 150000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 3) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        else {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 2) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 2) / 5;
        }
        break;
    case 11:
        depth_removal_ctrls->enabled = 1;

        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 32);

        dev_16x16_to_8x8_th = 4;
        dev_32x32_to_16x16_th = 4;

        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 50) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= dev_16x16_to_8x8_th * 4;
            dev_32x32_to_16x16_th *= dev_32x32_to_16x16_th * 2;
        }
        else {
            dev_16x16_to_8x8_th = dev_16x16_to_8x8_th * 2;
        }
        break;
#if TUNE_M10_DEPTH_ME
    case 12:
        depth_removal_ctrls->enabled = 1;

        disallow_below_64x64_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);
        disallow_below_32x32_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
        disallow_below_16x16_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 128);

        dev_16x16_to_8x8_th = 4;
        dev_32x32_to_16x16_th = 4;

        if (me_8x8_cost_variance < 50000) {
            dev_16x16_to_8x8_th *= (dev_16x16_to_8x8_th * 10) / 5;
            dev_32x32_to_16x16_th *= (dev_32x32_to_16x16_th * 10) / 5;
        }
        else if (me_8x8_cost_variance < 200000) {
            dev_16x16_to_8x8_th *= dev_16x16_to_8x8_th * 8;
            dev_32x32_to_16x16_th *= dev_32x32_to_16x16_th * 2;
        }
        else {
            dev_16x16_to_8x8_th = dev_16x16_to_8x8_th * 2;
        }
        break;
#endif
    }

    depth_removal_ctrls->disallow_below_64x64 = 0;
    depth_removal_ctrls->disallow_below_32x32 = 0;
    depth_removal_ctrls->disallow_below_16x16 = 0;

    if (depth_removal_ctrls->enabled) {

        uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[mdctxt->sb_index]);
        uint64_t cost_32x32 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_32x32_distortion[mdctxt->sb_index]);
        uint64_t cost_16x16 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_16x16_distortion[mdctxt->sb_index]);
        uint64_t cost_8x8 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_8x8_distortion[mdctxt->sb_index]);
        int64_t dev_32x32_to_16x16 =
            (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_16x16, 1)) * 100) /
            (int64_t)MAX(cost_16x16, 1);
        int64_t dev_32x32_to_8x8 =
            (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_8x8, 1)) * 100) /
            (int64_t)MAX(cost_8x8, 1);
        int64_t dev_32x32_to_8x8_th = (dev_32x32_to_16x16_th * 5) / 4;
        int64_t dev_16x16_to_8x8 =
            (int64_t)(((int64_t)MAX(cost_16x16, 1) - (int64_t)MAX(cost_8x8, 1)) * 100) /
            (int64_t)MAX(cost_8x8, 1);

        depth_removal_ctrls->disallow_below_64x64 = (sb_params->width % 64 == 0 && sb_params->height % 64 == 0)
            ? (cost_64x64 < disallow_below_64x64_th)
            : 0;

        depth_removal_ctrls->disallow_below_32x32 = (sb_params->width % 32 == 0 && sb_params->height % 32 == 0)
            ? (cost_32x32 < disallow_below_32x32_th || (dev_32x32_to_16x16 < dev_32x32_to_16x16_th && dev_32x32_to_8x8 < dev_32x32_to_8x8_th))
            : 0;

        depth_removal_ctrls->disallow_below_16x16 = (sb_params->width % 16 == 0 && sb_params->height % 16 == 0)
            ? (cost_16x16 < disallow_below_16x16_th || dev_16x16_to_8x8 < dev_16x16_to_8x8_th)
            : 0;
    }
}
#endif
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
#if FTR_USE_PSAD
        md_nsq_motion_search_ctrls->enable_psad = 1;
#endif
        break;

    case 2:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 15;
        md_nsq_motion_search_ctrls->full_pel_search_height = 15;
#if FTR_USE_PSAD
        md_nsq_motion_search_ctrls->enable_psad = 1;
#endif
        break;
    case 3:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
        md_nsq_motion_search_ctrls->full_pel_search_width = 11;
        md_nsq_motion_search_ctrls->full_pel_search_height = 11;
#if FTR_USE_PSAD
        md_nsq_motion_search_ctrls->enable_psad = 1;
#endif
        break;
    case 4:
        md_nsq_motion_search_ctrls->enabled = 1;
        md_nsq_motion_search_ctrls->use_ssd = 0;
#if FTR_USE_PSAD
        md_nsq_motion_search_ctrls->full_pel_search_width = 8;
        md_nsq_motion_search_ctrls->full_pel_search_height = 7;
#else
        md_nsq_motion_search_ctrls->full_pel_search_width = 7;
        md_nsq_motion_search_ctrls->full_pel_search_height = 7;
#endif
#if FTR_USE_PSAD
        md_nsq_motion_search_ctrls->enable_psad = 1;
#endif
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
#if TUNE_PME_M0
    case 1:
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 1;
        md_pme_ctrls->full_pel_search_width = 31;
        md_pme_ctrls->full_pel_search_height = 31;
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad = 0;
        break;
    case 2:
#else
    case 1:
#endif
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 1;
#if TUNE_PME_M0
        md_pme_ctrls->full_pel_search_width = 20;
        md_pme_ctrls->full_pel_search_height = 15;
#else
        md_pme_ctrls->full_pel_search_width = 15;
        md_pme_ctrls->full_pel_search_height = 15;
#endif
#if FTR_IMPROVE_PME
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
#endif
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
#if TUNE_BLOCK_SIZE
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
#endif
#if FTR_USE_PSAD
        md_pme_ctrls->enable_psad = 0;
#endif
        break;
#if TUNE_PME_M0
    case 3:
#else
    case 2:
#endif
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 1;
        md_pme_ctrls->full_pel_search_width = 7;
        md_pme_ctrls->full_pel_search_height = 5;
#if FTR_IMPROVE_PME
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
#endif
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = MIN_SIGNED_VALUE;
#if TUNE_BLOCK_SIZE
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
#endif
#if FTR_USE_PSAD
        md_pme_ctrls->enable_psad = 0;
#endif
        break;
#if TUNE_PME_M0
    case 4:
#else
    case 3:
#endif
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 1;
        md_pme_ctrls->full_pel_search_width = 7;
        md_pme_ctrls->full_pel_search_height = 5;
#if FTR_IMPROVE_PME
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
#endif
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 100;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
#if TUNE_BLOCK_SIZE
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
#endif
#if FTR_USE_PSAD
        md_pme_ctrls->enable_psad = 0;
#endif
        break;
#if FTR_USE_PSAD
#if TUNE_PME_M0
    case 5:
#else
    case 4:
#endif
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 8;
        md_pme_ctrls->full_pel_search_height = 5;
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 100;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad = 1;
        break;
#if TUNE_PME_M0
    case 6:
#else
    case 5:
#endif
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 8;
        md_pme_ctrls->full_pel_search_height = 5;
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad = 1;
        break;
#if TUNE_PME_M0
    case 7:
#else
    case 6:
#endif
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 8;
        md_pme_ctrls->full_pel_search_height = 3;
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad = 1;
        break;
#if TUNE_PME_M0
    case 8:
#else
    case 7:
#endif
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 3;
        md_pme_ctrls->full_pel_search_height = 3;
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad = 0;
        break;
#if TUNE_PME_M0
    case 9:
#else
    case 8:
#endif
#if TUNE_REG_PD1
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 3;
        md_pme_ctrls->full_pel_search_height = 3;
        md_pme_ctrls->early_check_mv_th_multiplier = 50;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 1;
        md_pme_ctrls->enable_psad = 0;
#else
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 3;
        md_pme_ctrls->full_pel_search_height = 3;
        md_pme_ctrls->early_check_mv_th_multiplier = 10;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 1;
        md_pme_ctrls->enable_psad = 0;
#endif
        break;
#if OPT_M11_PME
#if TUNE_PME_M0
    case 10:
#else
    case 9:
#endif
#if TUNE_REG_PD1
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 2;
        md_pme_ctrls->full_pel_search_height = 2;
        md_pme_ctrls->early_check_mv_th_multiplier = 80;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 30;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 20;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 10;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 40;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 1;
        md_pme_ctrls->enable_psad = 0;
        break;
#else
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 3;
        md_pme_ctrls->full_pel_search_height = 3;
        md_pme_ctrls->early_check_mv_th_multiplier = 50;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 1;
        md_pme_ctrls->enable_psad = 0;
        break;
#endif
#endif
#else
    case 4:
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 3;
        md_pme_ctrls->full_pel_search_height = 3;
#if FTR_IMPROVE_PME
        md_pme_ctrls->early_check_mv_th_multiplier = MIN_SIGNED_VALUE;
#endif
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
#if TUNE_BLOCK_SIZE
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
#endif
        break;
#if FTR_IMPROVE_PME
    case 5:
        md_pme_ctrls->enabled = 1;
        md_pme_ctrls->use_ssd = 0;
        md_pme_ctrls->full_pel_search_width = 3;
        md_pme_ctrls->full_pel_search_height = 3;
        md_pme_ctrls->early_check_mv_th_multiplier = 10;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th = 32;
#if TUNE_BLOCK_SIZE
        md_pme_ctrls->modulate_pme_for_blk_size_res = 1;
#endif
        break;
#endif
#endif
    default:
        assert(0);
        break;
    }
}
#if FTR_SUBSAMPLE_RESIDUAL
void set_subres_controls(ModeDecisionContext *mdctxt, uint8_t subres_level) {

    SubresCtrls *subres_ctrls = &mdctxt->subres_ctrls;

    switch (subres_level) {
    case 0:
        subres_ctrls->step = 0;
        break;
    case 1:
        subres_ctrls->step = 1;
        break;
#if TUNE_ADD_SUBRESS_FACTOR4
    case 2:
        subres_ctrls->step = 2;
        break;
#endif
    default:
        assert(0);
        break;
    }
#if FIX_SUBRES_R2R
    // Set the TH used to determine if subres is safe to use (based on ODD vs. EVEN rows' distortion)
    if (subres_ctrls->step == 0)
        subres_ctrls->odd_to_even_deviation_th = 0;
    else
        subres_ctrls->odd_to_even_deviation_th = 5;
#endif
}
#endif
void set_pf_controls(ModeDecisionContext *mdctxt, uint8_t pf_level) {

   PfCtrls *pf_ctrls = &mdctxt->pf_ctrls;

    switch (pf_level) {
    case 0:
        pf_ctrls->pf_shape = ONLY_DC_SHAPE;
        break;
    case 1:
        pf_ctrls->pf_shape = DEFAULT_SHAPE;
        break;
    case 2:
        pf_ctrls->pf_shape = N2_SHAPE;
        break;
    case 3:
        pf_ctrls->pf_shape = N4_SHAPE;
        break;
    default:
        assert(0);
        break;
    }
}
#if !CLN_INDEPTH
/*
 * Control in-depth block skip
 */
void set_in_depth_block_skip_ctrls(ModeDecisionContext *mdctxt, uint8_t in_depth_block_skip_level) {

    InDepthBlockSkipCtrls *in_depth_block_skip_ctrls = &mdctxt->in_depth_block_skip_ctrls;

    switch (in_depth_block_skip_level)
    {
    case 0:
        in_depth_block_skip_ctrls->base_weight                 =   0;
        break;

    case 1:
        in_depth_block_skip_ctrls->base_weight                 = 150;

#if TUNE_IN_DEPTH_LIST0_M9_LEVEL
        in_depth_block_skip_ctrls->cost_band_based_modulation  =   1;
        in_depth_block_skip_ctrls->max_cost_multiplier         =  25;
        in_depth_block_skip_ctrls->max_band_cnt                =   1;
        in_depth_block_skip_ctrls->weight_per_band[0]          = 100;
#else
        in_depth_block_skip_ctrls->cost_band_based_modulation  =   1;
        in_depth_block_skip_ctrls->max_cost_multiplier         = 400;
        in_depth_block_skip_ctrls->max_band_cnt                =   5;
        in_depth_block_skip_ctrls->weight_per_band[0]          = 175;
        in_depth_block_skip_ctrls->weight_per_band[1]          = 150;
        in_depth_block_skip_ctrls->weight_per_band[2]          = 125;
        in_depth_block_skip_ctrls->weight_per_band[3]          = 100;
        in_depth_block_skip_ctrls->weight_per_band[4]          =  75;
#endif

        in_depth_block_skip_ctrls->child_cnt_based_modulation  =   0;
        in_depth_block_skip_ctrls->cnt_based_weight[0]         = 150;
        in_depth_block_skip_ctrls->cnt_based_weight[1]         = 125;
        in_depth_block_skip_ctrls->cnt_based_weight[2]         = 100;
        break;

    case 2:
        in_depth_block_skip_ctrls->base_weight                 = 150;

        in_depth_block_skip_ctrls->cost_band_based_modulation  =   0;

        in_depth_block_skip_ctrls->child_cnt_based_modulation  =   0;
        break;

    default:
        assert(0);
        break;
    }
}
#endif
#if !CLN_REMOVE_UNUSED_FEATS
/*
 * Control lower-depth block skip
 */
void set_lower_depth_block_skip_ctrls(ModeDecisionContext *mdctxt, uint8_t lower_depth_block_skip_level) {

    LowerDepthBlockSkipCtrls *lower_depth_skip_ctrls = &mdctxt->lower_depth_block_skip_ctrls;

    switch (lower_depth_block_skip_level)
    {
    case 0:
        lower_depth_skip_ctrls->enabled             = 0;
        break;

    case 1:
        lower_depth_skip_ctrls->enabled                   =   1;
        lower_depth_skip_ctrls->quad_deviation_th         = 500;
        lower_depth_skip_ctrls->min_distortion_cost_ratio =  50;
        lower_depth_skip_ctrls->skip_all                  =   0;
        break;

    case 2:
        lower_depth_skip_ctrls->enabled                   =   1;
        lower_depth_skip_ctrls->quad_deviation_th         = 500;
        lower_depth_skip_ctrls->min_distortion_cost_ratio =  50;
        lower_depth_skip_ctrls->skip_all                  =   1;
        break;

    default:
        assert(0);
        break;
    }
}
#endif
#if OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS
void set_block_based_depth_refinement_controls(ModeDecisionContext* mdctxt, uint8_t block_based_depth_refinement_level) {

    DepthRefinementCtrls* depth_refinement_ctrls = &mdctxt->depth_refinement_ctrls;

    switch (block_based_depth_refinement_level)
    {
    case 0:
        depth_refinement_ctrls->enabled = 0;
        break;

    case 1:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =  25;
        depth_refinement_ctrls->sub_to_current_th          =  25;
        depth_refinement_ctrls->cost_band_based_modulation =   0;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;

    case 2:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =  25;
        depth_refinement_ctrls->sub_to_current_th          =  25;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  15;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;

    case 3:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =  25;
        depth_refinement_ctrls->sub_to_current_th          =  25;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;

    case 4:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =  20;
        depth_refinement_ctrls->sub_to_current_th          =  20;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;

    case 5:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =  15;
        depth_refinement_ctrls->sub_to_current_th          =  15;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        depth_refinement_ctrls->up_to_2_depth              =  0;
        break;

    case 6:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =  10;
        depth_refinement_ctrls->sub_to_current_th          =  10;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        depth_refinement_ctrls->up_to_2_depth              =  0;
        break;

    case 7:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =   5;
        depth_refinement_ctrls->sub_to_current_th          =   5;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 400;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;

    case 8:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =   5;
        depth_refinement_ctrls->sub_to_current_th          =   5;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 800;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;
#if TUNE_TXS_IFS_MFMV_DEPTH_M9
    case 9:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =   5;
        depth_refinement_ctrls->sub_to_current_th          = -50;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 800;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        depth_refinement_ctrls->up_to_2_depth              =   1;
        break;
#else
    case 9:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =   5;
        depth_refinement_ctrls->sub_to_current_th          =   5;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 1200;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]        = 5;
        break;
#endif
    case 10:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       = -25;
        depth_refinement_ctrls->sub_to_current_th          = -50;
        depth_refinement_ctrls->cost_band_based_modulation =   1;
        depth_refinement_ctrls->max_cost_multiplier        = 800;
        depth_refinement_ctrls->max_band_cnt               =   4;
        depth_refinement_ctrls->decrement_per_band[0]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[1]      =  MAX_SIGNED_VALUE;
        depth_refinement_ctrls->decrement_per_band[2]      =  10;
        depth_refinement_ctrls->decrement_per_band[3]      =   5;
        depth_refinement_ctrls->up_to_2_depth              =   1;
        break;

    // Pred_Only
    case 11:
        depth_refinement_ctrls->enabled                    =   1;
        depth_refinement_ctrls->parent_to_current_th       =   MIN_SIGNED_VALUE;
        depth_refinement_ctrls->sub_to_current_th          =   MIN_SIGNED_VALUE;
        depth_refinement_ctrls->cost_band_based_modulation =   0;
        depth_refinement_ctrls->up_to_2_depth              =   0;
        break;
    }
}
#endif
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
#if FTR_USE_PSAD
        md_sq_me_ctrls->enable_psad = 1;
#endif
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
#if FTR_USE_PSAD
        md_sq_me_ctrls->enable_psad = 1;
#endif
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
#if FTR_USE_PSAD
        md_sq_me_ctrls->enable_psad = 1;
#endif
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
#if FTR_USE_PSAD
        md_sq_me_ctrls->enable_psad = 1;
#endif
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
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = EIGHTH_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 1;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
#if OPT11_SUBPEL
        md_subpel_me_ctrls->skip_diag_refinement = 0;
#endif
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
    case 2:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
#if OPT11_SUBPEL
        md_subpel_me_ctrls->skip_diag_refinement = 0;
#endif
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
    case 3:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
#if OPT11_SUBPEL
        md_subpel_me_ctrls->skip_diag_refinement = 0;
#endif
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
#if FTR_LOW_AC_SUBPEL
    case 4:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = EIGHTH_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 1;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
        md_subpel_me_ctrls->skip_diag_refinement = 1;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
    case 5:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = EIGHTH_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 1;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
        md_subpel_me_ctrls->skip_diag_refinement = 2;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
    case 6:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
        md_subpel_me_ctrls->skip_diag_refinement = 2;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
    case 7:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
        md_subpel_me_ctrls->skip_diag_refinement = 1;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
    case 8:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
        md_subpel_me_ctrls->skip_diag_refinement = 2;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
#if OPT_SUBPEL
    case 9:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
        md_subpel_me_ctrls->skip_diag_refinement = 3;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
    case 10:
#else
    case 9:
#endif
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = HALF_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
        md_subpel_me_ctrls->skip_diag_refinement = 2;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
#if OPT_SUBPEL
    case 11:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = HALF_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_me_ctrls->pred_variance_th = 0;
        md_subpel_me_ctrls->abs_th_mult = 0;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
        md_subpel_me_ctrls->skip_diag_refinement = 3;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
#endif
#if OPT_M11_SUBPEL
    case 12:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th = 50;
        md_subpel_me_ctrls->abs_th_mult = 2;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement = 3;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
#if TUNE_M11_SUBPEL
    case 13:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th = 100;
        md_subpel_me_ctrls->abs_th_mult = 2;
        md_subpel_me_ctrls->round_dev_th = -25;
        md_subpel_me_ctrls->skip_diag_refinement = 4;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_me_ctrls->skip_zz_mv = 0;
#endif
        break;
#if OPT_SKIP_SUBPEL_ZZ
    case 14:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = QUARTER_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th = 100;
        md_subpel_me_ctrls->abs_th_mult = 2;
        md_subpel_me_ctrls->round_dev_th = -25;
        md_subpel_me_ctrls->skip_diag_refinement = 4;
        md_subpel_me_ctrls->skip_zz_mv = 1;
        break;
#endif
#if TUNE_4K_M11
    case 15:
#if FTR_VLPD1 // Fix levels to obey "onion ring" principle
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
#if CLN_SUBPEL_CTRLS
        md_subpel_me_ctrls->max_precision = HALF_PEL;
#else
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = HALF_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = HALF_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th = 100;
        md_subpel_me_ctrls->abs_th_mult = 2;
        md_subpel_me_ctrls->round_dev_th = -25;
        md_subpel_me_ctrls->skip_diag_refinement = 4;
        md_subpel_me_ctrls->skip_zz_mv = 1;
        break;
#else
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = HALF_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = HALF_PEL;
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th = 50;
        md_subpel_me_ctrls->abs_th_mult = 2;
        md_subpel_me_ctrls->round_dev_th = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement = 3;
        md_subpel_me_ctrls->skip_zz_mv = 0;
        break;
    case 16:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = HALF_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = HALF_PEL;
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th = 100;
        md_subpel_me_ctrls->abs_th_mult = 2;
        md_subpel_me_ctrls->round_dev_th = -25;
        md_subpel_me_ctrls->skip_diag_refinement = 4;
        md_subpel_me_ctrls->skip_zz_mv = 0;
        break;
#endif
#endif
#else
    case 13:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = HALF_PEL;
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th = 100;
        md_subpel_me_ctrls->abs_th_mult = 2;
        md_subpel_me_ctrls->round_dev_th = -25;
        md_subpel_me_ctrls->skip_diag_refinement = 3;
        break;
#endif
#endif
#else
    case 4:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->eight_pel_search_enabled = 1;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT11_SUBPEL
        md_subpel_me_ctrls->skip_diag_refinement = 1;
#endif
        break;
    case 5:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT11_SUBPEL
        md_subpel_me_ctrls->skip_diag_refinement = 1;
#endif
        break;
    case 6:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT11_SUBPEL
        md_subpel_me_ctrls->skip_diag_refinement = 1;
#endif
        break;
#if TUNE_CTR_QUARTER_PEL
    case 7:
        md_subpel_me_ctrls->enabled = 1;
        md_subpel_me_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->eight_pel_search_enabled = 0;
        md_subpel_me_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_me_ctrls->forced_stop_gt_32 = HALF_PEL;
        md_subpel_me_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT11_SUBPEL
        md_subpel_me_ctrls->skip_diag_refinement = 1;
#endif
        break;
#endif
#endif
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
#if CLN_SUBPEL_CTRLS
        md_subpel_pme_ctrls->max_precision = EIGHTH_PEL;
#else
        md_subpel_pme_ctrls->eight_pel_search_enabled = 1;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_pme_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_pme_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
#endif
        md_subpel_pme_ctrls->subpel_search_method = SUBPEL_TREE;
#if OPT_M11_SUBPEL
        md_subpel_pme_ctrls->pred_variance_th = 0;
        md_subpel_pme_ctrls->abs_th_mult = 0;
        md_subpel_pme_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_pme_ctrls->skip_zz_mv = 0;
#endif
        break;
    case 2:
        md_subpel_pme_ctrls->enabled = 1;
        md_subpel_pme_ctrls->subpel_search_type = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
#if CLN_SUBPEL_CTRLS
        md_subpel_pme_ctrls->max_precision = EIGHTH_PEL;
#else
        md_subpel_pme_ctrls->eight_pel_search_enabled = 1;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_pme_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_pme_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
#endif
        md_subpel_pme_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
#if OPT_M11_SUBPEL
        md_subpel_pme_ctrls->pred_variance_th = 0;
        md_subpel_pme_ctrls->abs_th_mult = 0;
        md_subpel_pme_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#endif
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_pme_ctrls->skip_zz_mv = 0;
#endif
        break;
#if OPT_M11_SUBPEL
    case 3:
        md_subpel_pme_ctrls->enabled = 1;
        md_subpel_pme_ctrls->subpel_search_type = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
#if CLN_SUBPEL_CTRLS
        md_subpel_pme_ctrls->max_precision = EIGHTH_PEL;
#else
        md_subpel_pme_ctrls->eight_pel_search_enabled = 1;
        md_subpel_pme_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_pme_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_pme_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_pme_ctrls->pred_variance_th = 0;
        md_subpel_pme_ctrls->abs_th_mult = 0;
        md_subpel_pme_ctrls->round_dev_th = MAX_SIGNED_VALUE;
#if OPT_SKIP_SUBPEL_ZZ
        md_subpel_pme_ctrls->skip_zz_mv = 0;
#endif
        break;
#if TUNE_4K_M11
    case 4:
        md_subpel_pme_ctrls->enabled = 1;
        md_subpel_pme_ctrls->subpel_search_type = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
#if CLN_SUBPEL_CTRLS
        md_subpel_pme_ctrls->max_precision = HALF_PEL;
#else
        md_subpel_pme_ctrls->eight_pel_search_enabled = 1;
        md_subpel_pme_ctrls->forced_stop_lt_eq_32 = HALF_PEL;
        md_subpel_pme_ctrls->forced_stop_gt_32 = HALF_PEL;
#endif
        md_subpel_pme_ctrls->subpel_search_method = SUBPEL_TREE_PRUNED;
        md_subpel_pme_ctrls->pred_variance_th = 0;
        md_subpel_pme_ctrls->abs_th_mult = 0;
        md_subpel_pme_ctrls->round_dev_th = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv = 0;
        break;
#endif
#else
    case 3:
        md_subpel_pme_ctrls->enabled = 1;
        md_subpel_pme_ctrls->subpel_search_type = USE_4_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 1;
        md_subpel_pme_ctrls->eight_pel_search_enabled = 0;
#if TUNE_CTR_QUARTER_PEL
        md_subpel_pme_ctrls->forced_stop_lt_eq_32 = EIGHTH_PEL;
        md_subpel_pme_ctrls->forced_stop_gt_32 = EIGHTH_PEL;
#endif
        md_subpel_pme_ctrls->subpel_search_method = SUBPEL_TREE;
        break;
#endif
    default: assert(0); break;
    }
}
#if !CLN_CAND_REDUCTION_CTRLS
void set_cand_elimination_controls(ModeDecisionContext *mdctxt, uint8_t eliminate_candidate_based_on_pme_me_results){
    CandEliminationCtlrs *cand_elimination_ctrls = &mdctxt->cand_elimination_ctrs;
    switch (eliminate_candidate_based_on_pme_me_results) {
    case 0:
        cand_elimination_ctrls->enabled = 0;
        break;
    case 1:
        cand_elimination_ctrls->enabled = 1;
        cand_elimination_ctrls->dc_only = 1;
        cand_elimination_ctrls->inject_new_me = 1;
        cand_elimination_ctrls->inject_new_pme = 1;
        cand_elimination_ctrls->inject_new_warp = 1;
#if OPT_EARLY_ELIM_TH
        cand_elimination_ctrls->th_multiplier = 1;
#endif
        break;
    case 2:
        cand_elimination_ctrls->enabled = 1;
        cand_elimination_ctrls->dc_only = 1;
        cand_elimination_ctrls->inject_new_me = 1;
        cand_elimination_ctrls->inject_new_pme = 1;
        cand_elimination_ctrls->inject_new_warp = 2;
#if OPT_EARLY_ELIM_TH
        cand_elimination_ctrls->th_multiplier = 1;
#endif
        break;
#if OPT_EARLY_ELIM_TH
    case 3:
        cand_elimination_ctrls->enabled = 1;
        cand_elimination_ctrls->dc_only = 1;
        cand_elimination_ctrls->inject_new_me = 1;
        cand_elimination_ctrls->inject_new_pme = 1;
        cand_elimination_ctrls->inject_new_warp = 2;
        cand_elimination_ctrls->th_multiplier = 3;
        break;
#endif
    default: assert(0); break;
    }
}
#endif
#if !CLN_MDS0_CTRLS
#if OPT_EARLY_CAND_ELIM
void set_early_cand_elimination_controls(ModeDecisionContext *mdctxt, uint8_t level) {
    EarlyCandElimCtrls *cand_elimination_ctrls = &mdctxt->early_cand_elimination_ctrls;
    switch (level) {
    case 0:
        cand_elimination_ctrls->enabled = 0;
        break;
    case 1:
        cand_elimination_ctrls->enabled = 1;
        cand_elimination_ctrls->mds0_distortion_th = 50;
        break;
    case 2:
        cand_elimination_ctrls->enabled = 1;
        cand_elimination_ctrls->mds0_distortion_th = 20;
        break;
    case 3:
        cand_elimination_ctrls->enabled = 1;
        cand_elimination_ctrls->mds0_distortion_th = 0;
        break;
    default: assert(0); break;
    }
}
#endif
#endif
/*
 * Control RDOQ
 */
void set_rdoq_controls(ModeDecisionContext *mdctxt, uint8_t rdoq_level) {
    RdoqCtrls *rdoq_ctrls = &mdctxt->rdoq_ctrls;

    switch (rdoq_level) {
    case 0:
        rdoq_ctrls->enabled = 0;
        break;
    case 1:
        rdoq_ctrls->enabled = 1;
#if FTR_RDOQ_ONLY_DCT_DCT
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y = 1;
        rdoq_ctrls->fp_q_uv = 1;
#else
        rdoq_ctrls->eob_fast_l_inter = 0;
        rdoq_ctrls->eob_fast_l_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_c_inter = 0;
#else
        rdoq_ctrls->eob_fast_c_inter = 1;
#endif
        rdoq_ctrls->eob_fast_c_intra = 0;
        rdoq_ctrls->fp_q_l = 1;
        rdoq_ctrls->fp_q_c = 1;
#endif
        rdoq_ctrls->satd_factor = (uint8_t)~0;
        rdoq_ctrls->early_exit_th = 0;
#if FTR_RDOQ_ONLY_DCT_DCT
        rdoq_ctrls->skip_uv          = 0;
        rdoq_ctrls->dct_dct_only     = 0;
#else
        rdoq_ctrls->disallow_md_rdoq_uv = 0;
        rdoq_ctrls->md_satd_factor = (uint8_t)~0;
#endif
#if OPT_RDOQ_EOB
        rdoq_ctrls->eob_th           = (uint8_t)~0;
#endif
#if FTR_RDO_OPT
        rdoq_ctrls->eob_fast_th = (uint8_t)~0;
#endif
        break;
#if TUNE_RDOQ_LEVEL_M9
    case 2:
        rdoq_ctrls->enabled = 1;
#if FTR_RDOQ_ONLY_DCT_DCT
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y = 1;
        rdoq_ctrls->fp_q_uv = 1;
#else
        rdoq_ctrls->eob_fast_l_inter = 0;
        rdoq_ctrls->eob_fast_l_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_c_inter = 0;
#else
        rdoq_ctrls->eob_fast_c_inter = 1;
#endif
        rdoq_ctrls->eob_fast_c_intra = 0;
        rdoq_ctrls->fp_q_l = 1;
        rdoq_ctrls->fp_q_c = 1;
#endif
        rdoq_ctrls->satd_factor = (uint8_t)~0;
        rdoq_ctrls->early_exit_th = 0;
#if FTR_RDOQ_ONLY_DCT_DCT
        rdoq_ctrls->skip_uv             = 1;
        rdoq_ctrls->dct_dct_only        = 0;
#else
        rdoq_ctrls->disallow_md_rdoq_uv = 1;
        rdoq_ctrls->md_satd_factor = (uint8_t)~0;
#endif
#if OPT_RDOQ_EOB
        rdoq_ctrls->eob_th              = (uint8_t)~0;
#endif
#if FTR_RDO_OPT
        rdoq_ctrls->eob_fast_th = (uint8_t)~0;
#endif
        break;
#if FTR_RDOQ_ONLY_DCT_DCT
    case 3:
        rdoq_ctrls->enabled = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y = 1;
        rdoq_ctrls->fp_q_uv = 1;
        rdoq_ctrls->satd_factor = (uint8_t)~0;
        rdoq_ctrls->early_exit_th = 0;
#if FTR_RDOQ_ONLY_DCT_DCT
        rdoq_ctrls->skip_uv             = 1;
        rdoq_ctrls->dct_dct_only        = 1;
#else
        rdoq_ctrls->disallow_md_rdoq_uv = 1;
        rdoq_ctrls->md_satd_factor = (uint8_t)~0;
#endif
#if OPT_RDOQ_EOB
        rdoq_ctrls->eob_th              = (uint8_t)~0;
#endif
#if FTR_RDO_OPT
        rdoq_ctrls->eob_fast_th = (uint8_t)~0;
#endif
        break;
#endif
#if FTR_RDO_OPT
    case 4:
        rdoq_ctrls->enabled = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y = 1;
        rdoq_ctrls->fp_q_uv = 1;
        rdoq_ctrls->satd_factor = (uint8_t)~0;
        rdoq_ctrls->early_exit_th = 0;
#if FTR_RDOQ_ONLY_DCT_DCT
        rdoq_ctrls->skip_uv = 1;
        rdoq_ctrls->dct_dct_only = 1;
#else
        rdoq_ctrls->disallow_md_rdoq_uv = 1;
        rdoq_ctrls->md_satd_factor = (uint8_t)~0;
#endif
#if OPT_RDOQ_EOB
        rdoq_ctrls->eob_th = (uint8_t)~0;
#endif
#if FTR_RDO_OPT
        rdoq_ctrls->eob_fast_th = 30;
#endif
        break;
#endif
#if OPT_RDOQ_EOB
#if FTR_RDO_OPT
    case 5:
#else
    case 4:
#endif
        rdoq_ctrls->enabled = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y = 1;
        rdoq_ctrls->fp_q_uv = 1;
        rdoq_ctrls->satd_factor = (uint8_t)~0;
        rdoq_ctrls->early_exit_th = 0;
        rdoq_ctrls->skip_uv = 1;
        rdoq_ctrls->dct_dct_only = 1;
        rdoq_ctrls->eob_th = 85;
#if CLN_RDOQ_CTRLS
        rdoq_ctrls->eob_fast_th = 0;
#endif
        break;
#endif
#else
    case 2:
        rdoq_ctrls->enabled = 1;
        rdoq_ctrls->eob_fast_l_inter = 0;
        rdoq_ctrls->eob_fast_l_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_c_inter = 0;
#else
        rdoq_ctrls->eob_fast_c_inter = 1;
#endif
        rdoq_ctrls->eob_fast_c_intra = 0;
        rdoq_ctrls->fp_q_l = 1;
        rdoq_ctrls->fp_q_c = 0;
        rdoq_ctrls->satd_factor = 128;
        rdoq_ctrls->early_exit_th = 5;
        rdoq_ctrls->disallow_md_rdoq_uv = 1;
        rdoq_ctrls->md_satd_factor = 64;
        break;
    case 3:
        rdoq_ctrls->enabled = 1;
        rdoq_ctrls->eob_fast_l_inter = 0;
        rdoq_ctrls->eob_fast_l_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_c_inter = 0;
#else
        rdoq_ctrls->eob_fast_c_inter = 1;
#endif
        rdoq_ctrls->eob_fast_c_intra = 0;
        rdoq_ctrls->fp_q_l = 1;
        rdoq_ctrls->fp_q_c = 0;
        rdoq_ctrls->satd_factor = 128;
        rdoq_ctrls->early_exit_th = 5;
        rdoq_ctrls->disallow_md_rdoq_uv = 1;
        rdoq_ctrls->md_satd_factor = 32;
        break;
#endif
    default: assert(0); break;
    }
}
/*
 * Settings for the parent SQ coeff-area based cycles reduction algorithm.
 */
void set_parent_sq_coeff_area_based_cycles_reduction_ctrls(ModeDecisionContext* ctx, uint8_t resolution, uint8_t cycles_alloc_lvl) {
    ParentSqCoeffAreaBasedCyclesReductionCtrls* cycle_red_ctrls = &ctx->parent_sq_coeff_area_based_cycles_reduction_ctrls;
    switch (cycles_alloc_lvl) {
    case 0:
        cycle_red_ctrls->enabled = 0;
        break;
    case 1:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th = resolution <= INPUT_SIZE_360p_RANGE ? UNUSED_HIGH_FREQ_BAND_TH : 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th = resolution <= INPUT_SIZE_360p_RANGE ? UNUSED_HIGH_FREQ_BAND_TH : 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action = 1;
        cycle_red_ctrls->enable_one_coeff_action = 0;
        cycle_red_ctrls->one_coeff_action = 0;

        cycle_red_ctrls->low_freq_band1_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 2:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th = resolution <= INPUT_SIZE_360p_RANGE ? UNUSED_HIGH_FREQ_BAND_TH : 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th = resolution <= INPUT_SIZE_360p_RANGE ? UNUSED_HIGH_FREQ_BAND_TH : 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action = 2;
        cycle_red_ctrls->enable_one_coeff_action = 0;
        cycle_red_ctrls->one_coeff_action = 0;

        cycle_red_ctrls->low_freq_band1_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 3:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th = 90;
        cycle_red_ctrls->high_freq_band1_level = resolution <= INPUT_SIZE_360p_RANGE ? 2 : 3;
        cycle_red_ctrls->high_freq_band2_th = 70;
        cycle_red_ctrls->high_freq_band2_level = resolution <= INPUT_SIZE_360p_RANGE ? 1 : 3;
        cycle_red_ctrls->high_freq_band3_th = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action = 2;
        cycle_red_ctrls->enable_one_coeff_action = 0;
        cycle_red_ctrls->one_coeff_action = 0;

        cycle_red_ctrls->low_freq_band1_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 4:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th = 90;
        cycle_red_ctrls->high_freq_band1_level = resolution <= INPUT_SIZE_360p_RANGE ? 2 : 3;
        cycle_red_ctrls->high_freq_band2_th = 70;
        cycle_red_ctrls->high_freq_band2_level = resolution <= INPUT_SIZE_360p_RANGE ? 1 : 3;
        cycle_red_ctrls->high_freq_band3_th = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action = 3;
        cycle_red_ctrls->enable_one_coeff_action = 1;
        cycle_red_ctrls->one_coeff_action = 1;

        cycle_red_ctrls->low_freq_band1_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 5:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th = 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th = 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th = 50;
        cycle_red_ctrls->high_freq_band3_level = 1;


        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action = 2;
        cycle_red_ctrls->enable_one_coeff_action = 0;
        cycle_red_ctrls->one_coeff_action = 0;

        cycle_red_ctrls->low_freq_band1_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    case 6:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th = 90;
        cycle_red_ctrls->high_freq_band1_level = resolution <= INPUT_SIZE_360p_RANGE ? 2 : 3;
        cycle_red_ctrls->high_freq_band2_th = 70;
        cycle_red_ctrls->high_freq_band2_level = resolution <= INPUT_SIZE_360p_RANGE ? 1 : 3;
        cycle_red_ctrls->high_freq_band3_th = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action = 0;
        cycle_red_ctrls->enable_one_coeff_action = 1;
        cycle_red_ctrls->one_coeff_action = 1;

        cycle_red_ctrls->low_freq_band1_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;

    case 7:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th = 90;
        cycle_red_ctrls->high_freq_band1_level = 0;
        cycle_red_ctrls->high_freq_band2_th = 70;
        cycle_red_ctrls->high_freq_band2_level = 3;
        cycle_red_ctrls->high_freq_band3_th = 50;
        cycle_red_ctrls->high_freq_band3_level = 2;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action = 0;
        cycle_red_ctrls->enable_one_coeff_action = 1;
        cycle_red_ctrls->one_coeff_action = 1;

        cycle_red_ctrls->low_freq_band1_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band1_level = 0;
        cycle_red_ctrls->low_freq_band2_th = UNUSED_LOW_FREQ_BAND_TH;
        cycle_red_ctrls->low_freq_band2_level = 0;
        break;
    default:
        assert(0);
        break;
    }
}
#if TUNE_NEW_TXT_LVLS
void set_txt_controls(ModeDecisionContext *mdctxt, uint8_t txt_level, uint8_t resolution) {
#else
void set_txt_controls(ModeDecisionContext *mdctxt, uint8_t txt_level) {
#endif

    TxtControls * txt_ctrls = &mdctxt->txt_ctrls;

    switch (txt_level)
    {
    case 0:
        txt_ctrls->enabled = 0;

        txt_ctrls->txt_group_inter_lt_16x16 = 1;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;

        txt_ctrls->txt_group_intra_lt_16x16 = 1;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 1;
#if TUNE_NEW_TXT_LVLS
        txt_ctrls->early_exit_dist_th = 0;
        txt_ctrls->early_exit_coeff_th = 0;
#endif
        break;
    case 1:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = MAX_TX_TYPE_GROUP;

        txt_ctrls->txt_group_intra_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
#if TUNE_NEW_TXT_LVLS
        txt_ctrls->early_exit_dist_th = 0;
        txt_ctrls->early_exit_coeff_th = 0;
#endif
        break;
    case 2:
        txt_ctrls->enabled = 1;
        txt_ctrls->txt_group_inter_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 5;

        txt_ctrls->txt_group_intra_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
#if TUNE_NEW_TXT_LVLS
        txt_ctrls->early_exit_dist_th = 0;
        txt_ctrls->early_exit_coeff_th = 0;
#endif
        break;
    case 3:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = 5;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 5;

        txt_ctrls->txt_group_intra_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
#if TUNE_NEW_TXT_LVLS
        txt_ctrls->early_exit_dist_th = 0;
        txt_ctrls->early_exit_coeff_th = 0;
#endif
        break;
    case 4:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = 5;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 3;

        txt_ctrls->txt_group_intra_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
#if TUNE_NEW_TXT_LVLS
        txt_ctrls->early_exit_dist_th = 0;
        txt_ctrls->early_exit_coeff_th = 0;
#endif
        break;
    case 5:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 4;
#if TUNE_NEW_TXT_LVLS
        txt_ctrls->early_exit_dist_th = 0;
        txt_ctrls->early_exit_coeff_th = 0;
#endif
        break;
#if TUNE_NEW_TXT_LVLS
    case 6:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 4;

        txt_ctrls->early_exit_dist_th = 100;
        txt_ctrls->early_exit_coeff_th = (resolution <= INPUT_SIZE_720p_RANGE) ? 4 : 16;
        break;
#if TUNE_REG_PD1
    case 7:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 3;

        txt_ctrls->early_exit_dist_th = 400;
        txt_ctrls->early_exit_coeff_th = (resolution <= INPUT_SIZE_720p_RANGE) ? 8 : 32;
        break;
    case 8:
#else
    case 7:
#endif
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16 = 4;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
#if TUNE_REG_PD1
        txt_ctrls->early_exit_dist_th = 400;
        txt_ctrls->early_exit_coeff_th = (resolution <= INPUT_SIZE_720p_RANGE) ? 8 : 32;
#else
        txt_ctrls->early_exit_dist_th = 300;
        txt_ctrls->early_exit_coeff_th = (resolution <= INPUT_SIZE_720p_RANGE) ? 4 : 16;
#endif
        break;
#if TUNE_TXT_LEVEL
#if TUNE_REG_PD1
    case 9:
#else
    case 8:
#endif
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;

        txt_ctrls->txt_group_intra_lt_16x16 = 4;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
#if TUNE_REG_PD1
        txt_ctrls->early_exit_dist_th = 400;
        txt_ctrls->early_exit_coeff_th = (resolution <= INPUT_SIZE_720p_RANGE) ? 8 : 32;
#else
        txt_ctrls->early_exit_dist_th = 300;
        txt_ctrls->early_exit_coeff_th = (resolution <= INPUT_SIZE_720p_RANGE) ? 4 : 16;
#endif
        break;
#endif
#if TUNE_REG_PD1
    case 10:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16 = 2;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;
        txt_ctrls->txt_group_intra_lt_16x16 = 3;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
        txt_ctrls->early_exit_dist_th = 400;
        txt_ctrls->early_exit_coeff_th = (resolution <= INPUT_SIZE_720p_RANGE) ? 8 : 32;
        break;
#endif
#else
    case 6:
        txt_ctrls->enabled = 1;

#if TUNE_TXT_LEVEL_M9
        txt_ctrls->txt_group_inter_lt_16x16 = 2;
#else
        txt_ctrls->txt_group_inter_lt_16x16 = 3;
#endif
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;
        txt_ctrls->txt_group_intra_lt_16x16 = 3;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 1;

        break;
#endif
    default:
        assert(0);
        break;
    }
}
#if TUNE_BLOCK_SIZE
void set_interpolation_search_level_ctrls(ModeDecisionContext* context_ptr, uint8_t interpolation_search_level) {

    InterpolationSearchCtrls* ifs_ctrls = &context_ptr->ifs_ctrls;

    switch (interpolation_search_level) {
    case 0:
        ifs_ctrls->interpolation_search_level = IFS_OFF;
        ifs_ctrls->quarter_pel_only = 0;
        ifs_ctrls->modulate_filter_per_resolution = 0;
#if OPT_IFS
        ifs_ctrls->early_skip = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model = 0;
#endif
        break;
    case 1:
        ifs_ctrls->interpolation_search_level = IFS_MDS0;
        ifs_ctrls->quarter_pel_only = 0;
        ifs_ctrls->modulate_filter_per_resolution = 0;
#if OPT_IFS
        ifs_ctrls->early_skip = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model = 0;
#endif
        break;
    case 2:
        ifs_ctrls->interpolation_search_level = IFS_MDS1;
        ifs_ctrls->quarter_pel_only = 0;
        ifs_ctrls->modulate_filter_per_resolution = 0;
#if OPT_IFS
        ifs_ctrls->early_skip = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model = 0;
#endif
        break;
    case 3:
        ifs_ctrls->interpolation_search_level = IFS_MDS2;
        ifs_ctrls->quarter_pel_only = 0;
        ifs_ctrls->modulate_filter_per_resolution = 0;
#if OPT_IFS
        ifs_ctrls->early_skip = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model = 0;
#endif
        break;
    case 4:
        ifs_ctrls->interpolation_search_level = IFS_MDS3;
        ifs_ctrls->quarter_pel_only = 0;
        ifs_ctrls->modulate_filter_per_resolution = 0;
#if OPT_IFS
        ifs_ctrls->early_skip = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model = 0;
#endif
        break;
    case 5:
        ifs_ctrls->interpolation_search_level = IFS_MDS3;
        ifs_ctrls->quarter_pel_only = 1;
        ifs_ctrls->modulate_filter_per_resolution = 1;
#if OPT_IFS
        ifs_ctrls->early_skip = 0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model = 0;
#endif
        break;
#if OPT_IFS
    case 6:
        ifs_ctrls->interpolation_search_level = IFS_MDS3;
        ifs_ctrls->quarter_pel_only = 0;
        ifs_ctrls->modulate_filter_per_resolution = 0;
        ifs_ctrls->early_skip = 1;
        ifs_ctrls->subsampled_distortion = 1;
        ifs_ctrls->skip_sse_rd_model = 1;
        break;
    case 7:
        ifs_ctrls->interpolation_search_level = IFS_MDS3;
        ifs_ctrls->quarter_pel_only = 1;
        ifs_ctrls->modulate_filter_per_resolution = 1;
        ifs_ctrls->early_skip = 1;
        ifs_ctrls->subsampled_distortion = 1;
        ifs_ctrls->skip_sse_rd_model = 1;
        break;
#endif
    default:
        assert(0);
        break;
    }
}
#endif
#if CLN_CAND_REDUCTION_CTRLS
void set_cand_reduction_ctrls(PictureControlSet* pcs_ptr, ModeDecisionContext* mdctxt, uint8_t cand_reduction_level,
    const uint32_t picture_qp,
    uint32_t me_8x8_cost_variance,
    uint32_t me_64x64_distortion,
    uint8_t l0_was_skip, uint8_t l1_was_skip, uint8_t ref_skip_perc) {

    CandReductionCtrls* cand_reduction_ctrls = &mdctxt->cand_reduction_ctrls;

    switch (cand_reduction_level)
    {
    case 0:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 0;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 0;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count = 3;
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
        cand_reduction_ctrls->near_count_ctrls.enabled = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count = 3;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier = 1;

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
        cand_reduction_ctrls->near_count_ctrls.enabled = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 1;

        break;

    case 3:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 2;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 1;

        break;

    case 4:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 1;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 2;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates =
            (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag || ((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
                me_8x8_cost_variance < (500 * picture_qp) &&
                me_64x64_distortion < (500 * picture_qp))) ?
            3 : 1;

        break;

    case 5:
        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count = 0;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 1;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 2;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates =
            (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag || ((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
                me_8x8_cost_variance < (500 * picture_qp) &&
                me_64x64_distortion < (500 * picture_qp))) ?
            3 : 1;

        break;

    default:
        assert(0);
        break;
    }

    // lpd1_mvp_best_me_list can only use this feature when a single unipred ME candidate is selected,
    if (!(pcs_ptr->parent_pcs_ptr->ref_list0_count_try == 1 && pcs_ptr->parent_pcs_ptr->ref_list1_count_try == 1 && pcs_ptr->parent_pcs_ptr->use_best_me_unipred_cand_only))
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;
}
#endif
#if !CLN_CAND_REDUCTION_CTRLS
void set_near_count_ctrls(ModeDecisionContext* mdctxt, uint8_t near_count_level) {

    NearCountCtrls* near_count_ctrls = &mdctxt->near_count_ctrls;

    switch (near_count_level)
    {
    case 0:
        near_count_ctrls->enabled = 0;

        near_count_ctrls->near_count = 0;
        near_count_ctrls->near_near_count = 0;

        break;
    case 1:
        near_count_ctrls->enabled = 1;

        near_count_ctrls->near_count = 3;
        near_count_ctrls->near_near_count = 3;

        break;
    case 2:
        near_count_ctrls->enabled = 1;
        near_count_ctrls->near_count = 1;
        near_count_ctrls->near_near_count = 3;
        break;

    case 3:
        near_count_ctrls->enabled = 1;
        near_count_ctrls->near_count = 1;
        near_count_ctrls->near_near_count = 1;
        break;

    case 4:
        near_count_ctrls->enabled = 1;
        near_count_ctrls->near_count = 0;
#if FTR_VLPD1
        near_count_ctrls->near_near_count = 1;
#else
        near_count_ctrls->near_near_count = 0;
#endif
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if CHROMA_CLEANUP
void set_chroma_controls(ModeDecisionContext *mdctxt, uint8_t uv_level) {
    UvCtrls* uv_ctrls = &mdctxt->uv_ctrls;

    switch (uv_level)
    {
    case 0:
        uv_ctrls->enabled = 0;
        uv_ctrls->uv_mode = CHROMA_MODE_2;
        uv_ctrls->nd_uv_serach_mode = 0;
        break;
    case 1:
        uv_ctrls->enabled = 1;
        uv_ctrls->uv_mode = CHROMA_MODE_0;
        uv_ctrls->nd_uv_serach_mode = 0;
        uv_ctrls->uv_intra_th = (uint64_t)~0;
        uv_ctrls->uv_cfl_th = (uint64_t)~0;
        break;
    case 2:
        uv_ctrls->enabled = 1;
        uv_ctrls->uv_mode = CHROMA_MODE_0;
        uv_ctrls->nd_uv_serach_mode = 1;
        uv_ctrls->uv_intra_th = 130;
        uv_ctrls->uv_cfl_th = 130;
        break;
     case 3:
        uv_ctrls->enabled = 1;
        uv_ctrls->uv_mode = CHROMA_MODE_0;
        uv_ctrls->nd_uv_serach_mode = 1;
        uv_ctrls->uv_intra_th = 100;
        uv_ctrls->uv_cfl_th = 100;
        break;
    case 4:
        uv_ctrls->enabled = 1;
        uv_ctrls->uv_mode = CHROMA_MODE_1;
        uv_ctrls->nd_uv_serach_mode = 0;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if FTR_NEW_WM_LVL
void set_wm_controls(ModeDecisionContext *mdctxt, uint8_t wm_level) {
    WmCtrls* wm_ctrls = &mdctxt->wm_ctrls;

    switch (wm_level)
    {
    case 0:
        wm_ctrls->enabled = 0;
        break;
    case 1:
        wm_ctrls->enabled = 1;
        wm_ctrls->use_wm_for_mvp = 1;
        wm_ctrls->num_new_mv_refinement = 12;
        break;
    case 2:
        wm_ctrls->enabled = 1;
        wm_ctrls->use_wm_for_mvp = 0;
        wm_ctrls->num_new_mv_refinement = 0;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if CLN_NIC_SIGS
// Get the nic_level used for each preset (to be passed to setting function: set_nic_controls())
#if CLN_REG_PD_SIG_SET_2
uint8_t get_nic_level(EbEncMode enc_mode, uint8_t temporal_layer_index) {

    uint8_t nic_level;

    if (enc_mode <= ENC_MRS)
#else
uint8_t get_nic_level(PdPass pd_pass, EbEncMode enc_mode, uint8_t temporal_layer_index) {

    uint8_t nic_level;

    if (pd_pass == PD_PASS_0)
#if CLN_MD_STAGING_CTRLS
        nic_level = 17;
#else
        nic_level = 16;
#endif
    else if (enc_mode <= ENC_MRS)
#endif
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
    else if (enc_mode <= ENC_M6)
        nic_level = (temporal_layer_index == 0) ? 12 : 14;
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

    NicPruningCtrls* nic_pruning_ctrls = ctx ? &ctx->nic_ctrls.pruning_ctrls : NULL;
    uint8_t nic_scaling_level = 0;
#if CLN_MD_STAGING_CTRLS
    uint8_t md_staging_mode = MD_STAGING_MODE_0;
#endif
    switch (nic_level)
    {
    case 0: // MAX NIC scaling; no pruning
        // NIC scaling level
        nic_scaling_level = 0;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_class_th = (uint64_t)~0;
            nic_pruning_ctrls->mds3_class_th = (uint64_t)~0;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_cand_base_th = (uint64_t)~0;
            nic_pruning_ctrls->mds3_cand_base_th = (uint64_t)~0;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_1;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_cand_base_th = 50;
            nic_pruning_ctrls->mds3_cand_base_th = 50;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_cand_base_th = 50;
            nic_pruning_ctrls->mds3_cand_base_th = 50;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 500;
            nic_pruning_ctrls->mds2_cand_base_th = 30;
            nic_pruning_ctrls->mds3_cand_base_th = 30;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 20;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_2;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 300;
            nic_pruning_ctrls->mds2_cand_base_th = 20;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_1;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 200;
            nic_pruning_ctrls->mds2_cand_base_th = 15;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_1;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 200;
            nic_pruning_ctrls->mds2_cand_base_th = 15;
            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_1;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 50;
            nic_pruning_ctrls->mds2_cand_base_th = 5;
            nic_pruning_ctrls->mds3_cand_base_th = 5;
        }
#if CLN_MD_STAGING_CTRLS
        md_staging_mode = MD_STAGING_MODE_1;
#endif
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 1;
            nic_pruning_ctrls->mds2_cand_base_th = 1;
            nic_pruning_ctrls->mds3_cand_base_th = 1;
        }
#if CLN_MD_STAGING_CTRLS
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
            nic_pruning_ctrls->force_1_cand_th = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th = 1;
            nic_pruning_ctrls->mds2_cand_base_th = 1;
            nic_pruning_ctrls->mds3_cand_base_th = 1;
        }
        md_staging_mode = MD_STAGING_MODE_0;
#endif
        break;
    default:
        assert(0);
        break;
    }

    if (ctx) {
        NicScalingCtrls* nic_scaling_ctrls = &ctx->nic_ctrls.scaling_ctrls;
        nic_scaling_ctrls->stage1_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_1];
        nic_scaling_ctrls->stage2_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_2];
        nic_scaling_ctrls->stage3_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_3];
#if CLN_MD_STAGING_CTRLS
        ctx->nic_ctrls.md_staging_mode = md_staging_mode;
#endif
    }

    return nic_scaling_level;
}
#else
void set_nic_controls(ModeDecisionContext *mdctxt, uint8_t nic_scaling_level) {

    NicCtrls* nic_ctrls = &mdctxt->nic_ctrls;
        nic_ctrls->stage1_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_1];
        nic_ctrls->stage2_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_2];
        nic_ctrls->stage3_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_3];

}
void set_nic_pruning_controls(ModeDecisionContext *mdctxt, uint8_t nic_pruning_level) {

    NicPruningCtrls* nic_pruning_ctrls = &mdctxt->nic_pruning_ctrls;

    switch (nic_pruning_level)
    {
    case 0:
        nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;
        nic_pruning_ctrls->mds2_class_th = (uint64_t)~0;
        nic_pruning_ctrls->mds3_class_th = (uint64_t)~0;

        nic_pruning_ctrls->mds1_cand_base_th = (uint64_t)~0;
        nic_pruning_ctrls->mds2_cand_base_th = (uint64_t)~0;
        nic_pruning_ctrls->mds3_cand_base_th = (uint64_t)~0;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;

    case 1:
        nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;

        nic_pruning_ctrls->mds2_class_th = 25;
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds2_band_cnt = 4;
#else
        nic_pruning_ctrls->mds2_band_cnt = 3;
#endif

        nic_pruning_ctrls->mds3_class_th = 25;
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds3_band_cnt = 4;
#else
        nic_pruning_ctrls->mds3_band_cnt = 3;
#endif

        nic_pruning_ctrls->mds1_cand_base_th = (uint64_t)~0;

#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds2_cand_base_th = 50;
        nic_pruning_ctrls->mds3_cand_base_th = 50;
#else
        nic_pruning_ctrls->mds2_cand_base_th = 45;

        nic_pruning_ctrls->mds3_cand_base_th = 45;
#endif
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;

    case 2:
        nic_pruning_ctrls->mds1_class_th = 300;
        nic_pruning_ctrls->mds1_band_cnt = 2;

        nic_pruning_ctrls->mds2_class_th = 25;
        nic_pruning_ctrls->mds2_band_cnt = 3;

        nic_pruning_ctrls->mds3_class_th = 15;
        nic_pruning_ctrls->mds3_band_cnt = 3;

        nic_pruning_ctrls->mds1_cand_base_th = 300;

        nic_pruning_ctrls->mds2_cand_base_th = 15;

        nic_pruning_ctrls->mds3_cand_base_th = 15;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;

    case 3:
        nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;
        nic_pruning_ctrls->mds1_band_cnt = 2;


        nic_pruning_ctrls->mds2_class_th = 25;
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds2_band_cnt = 4;
#else
        nic_pruning_ctrls->mds2_band_cnt = 3;
#endif

        nic_pruning_ctrls->mds3_class_th = 25;
        nic_pruning_ctrls->mds3_band_cnt = 12;

        nic_pruning_ctrls->mds1_cand_base_th = 300;

        nic_pruning_ctrls->mds2_cand_base_th = 20;

        nic_pruning_ctrls->mds3_cand_base_th = 20;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;

    case 4:
        nic_pruning_ctrls->mds1_class_th = 300;
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds1_band_cnt = 4;
#else
        nic_pruning_ctrls->mds1_band_cnt = 6;
#endif

        nic_pruning_ctrls->mds2_class_th = 25;
        nic_pruning_ctrls->mds2_band_cnt = 10;

        nic_pruning_ctrls->mds3_class_th = 15;
        nic_pruning_ctrls->mds3_band_cnt = 16;

        nic_pruning_ctrls->mds1_cand_base_th = 300;

        nic_pruning_ctrls->mds2_cand_base_th = 20;

        nic_pruning_ctrls->mds3_cand_base_th = 15;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;

    case 5:
        nic_pruning_ctrls->mds1_class_th = 200;
        nic_pruning_ctrls->mds1_band_cnt = 16;

        nic_pruning_ctrls->mds2_class_th = 25;
        nic_pruning_ctrls->mds2_band_cnt = 10;

        nic_pruning_ctrls->mds3_class_th = 15;
        nic_pruning_ctrls->mds3_band_cnt = 16;

        nic_pruning_ctrls->mds1_cand_base_th = 200;

        nic_pruning_ctrls->mds2_cand_base_th = 15;

        nic_pruning_ctrls->mds3_cand_base_th = 15;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;

    case 6:
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds1_class_th = 200;
        nic_pruning_ctrls->mds1_band_cnt = 16;
#else
        nic_pruning_ctrls->mds1_class_th = 100;
        nic_pruning_ctrls->mds1_band_cnt = 2;
#endif

        nic_pruning_ctrls->mds2_class_th = 25;
        nic_pruning_ctrls->mds2_band_cnt = 3;

        nic_pruning_ctrls->mds3_class_th = 15;
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds3_band_cnt = 16;
#else
        nic_pruning_ctrls->mds3_band_cnt = 3;
#endif
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds1_cand_base_th = 50;
#else
        nic_pruning_ctrls->mds1_cand_base_th = 45;
#endif

        nic_pruning_ctrls->mds2_cand_base_th = 15;

        nic_pruning_ctrls->mds3_cand_base_th = 15;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;

    case 7:
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds1_class_th = 200;
        nic_pruning_ctrls->mds1_band_cnt = 16;
#else
        nic_pruning_ctrls->mds1_class_th = 100;
        nic_pruning_ctrls->mds1_band_cnt = 2;
#endif

        nic_pruning_ctrls->mds2_class_th = 10;
        nic_pruning_ctrls->mds2_band_cnt = 2;

        nic_pruning_ctrls->mds3_class_th = 10;
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds3_band_cnt = 16;
#else
        nic_pruning_ctrls->mds3_band_cnt = 2;
#endif
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds1_cand_base_th = 50;
#else
        nic_pruning_ctrls->mds1_cand_base_th = 45;
#endif

        nic_pruning_ctrls->mds2_cand_base_th = 5;

        nic_pruning_ctrls->mds3_cand_base_th = 5;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;
    case 8:
        nic_pruning_ctrls->mds1_class_th = 100;
        nic_pruning_ctrls->mds1_band_cnt = 2;

        nic_pruning_ctrls->mds2_class_th = 10;
        nic_pruning_ctrls->mds2_band_cnt = 2;

        nic_pruning_ctrls->mds3_class_th = 10;
        nic_pruning_ctrls->mds3_band_cnt = 2;

        nic_pruning_ctrls->mds1_cand_base_th = 45;

        nic_pruning_ctrls->mds2_cand_base_th = 1;

        nic_pruning_ctrls->mds3_cand_base_th = 1;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;
    case 9:
        nic_pruning_ctrls->mds1_class_th = 100;
        nic_pruning_ctrls->mds1_band_cnt = 16;
        nic_pruning_ctrls->mds2_class_th = 5;
        nic_pruning_ctrls->mds2_band_cnt = 10;

        nic_pruning_ctrls->mds3_class_th = 10;
#if CLN_NIC_PRUNING
        nic_pruning_ctrls->mds3_band_cnt = 16;
#else
        nic_pruning_ctrls->mds3_band_cnt = 2;
#endif

        nic_pruning_ctrls->mds1_cand_base_th = 50;
        nic_pruning_ctrls->mds2_cand_base_th = 5;
        nic_pruning_ctrls->mds3_cand_base_th = 1;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif

        break;
    case 10:
        nic_pruning_ctrls->mds1_class_th = 100;
        nic_pruning_ctrls->mds1_band_cnt = 16;
        nic_pruning_ctrls->mds2_class_th = 2;
        nic_pruning_ctrls->mds2_band_cnt = 10;

        nic_pruning_ctrls->mds3_class_th = 10;
        nic_pruning_ctrls->mds3_band_cnt = 2;

        nic_pruning_ctrls->mds1_cand_base_th = 20;
        nic_pruning_ctrls->mds2_cand_base_th = 1;
        nic_pruning_ctrls->mds3_cand_base_th = 1;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;
    case 11:
        nic_pruning_ctrls->mds1_class_th = 75;
        nic_pruning_ctrls->mds1_band_cnt = 16;

        nic_pruning_ctrls->mds2_class_th = 1;
        nic_pruning_ctrls->mds2_band_cnt = 2;

        nic_pruning_ctrls->mds3_class_th = 1;
        nic_pruning_ctrls->mds3_band_cnt = 2;

        nic_pruning_ctrls->mds1_cand_base_th = 1;
        nic_pruning_ctrls->mds2_cand_base_th = 1;
        nic_pruning_ctrls->mds3_cand_base_th = 1;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 0;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;
#if OPT_NIC_PRUNING
    case 12:
        nic_pruning_ctrls->mds1_class_th = 75;
        nic_pruning_ctrls->mds1_band_cnt = 16;

        nic_pruning_ctrls->mds2_class_th = 0;
        nic_pruning_ctrls->mds2_band_cnt = 2;

        nic_pruning_ctrls->mds3_class_th = 0;
        nic_pruning_ctrls->mds3_band_cnt = 2;

        nic_pruning_ctrls->mds1_cand_base_th = 1;
        nic_pruning_ctrls->mds2_cand_base_th = 1;
        nic_pruning_ctrls->mds3_cand_base_th = 1;
#if CLN_NIC_PRUNE_CTRLS
        nic_pruning_ctrls->enable_skipping_mds1 = 1;
        nic_pruning_ctrls->force_1_cand_th = 0;
#endif
        break;
#endif
    default:
        assert(0);
        break;
    }
}
#endif
void set_inter_intra_ctrls(ModeDecisionContext* mdctxt, uint8_t inter_intra_level) {

    InterIntraCompCtrls* ii_ctrls = &mdctxt->inter_intra_comp_ctrls;

    switch (inter_intra_level) {
    case 0:
        ii_ctrls->enabled = 0;
        break;
    case 1:
        ii_ctrls->enabled = 1;
        break;
    default:
        assert(0);
        break;
    }
}
void set_depth_ctrls(ModeDecisionContext* ctx, uint8_t depth_level) {
    DepthCtrls* depth_ctrls = &ctx->depth_ctrls;

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
    default:
        assert(0);
        break;
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
#if TUNE_M5_M6
    else if (enc_mode <= ENC_M4)
#else
    else if (enc_mode <= ENC_M6)
#endif
        disallow_4x4 = (slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
    else
        disallow_4x4 = EB_TRUE;

    return disallow_4x4;
}
#if FTR_LPD1_DETECTOR
#if FTR_VLPD1
#if CLN_LPD1_LVLS
void set_lpd1_ctrls(ModeDecisionContext *ctx, uint8_t lpd1_lvl) {
    Lpd1Ctrls* ctrls = &ctx->lpd1_ctrls;
    switch (lpd1_lvl) {
    case 0:
        ctrls->pd1_level = REGULAR_PD1; // Light-PD1 path not used
        break;
    case 1:
        ctrls->pd1_level = LPD1_LVL_0;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0] = 1;
        ctrls->use_ref_info[LPD1_LVL_0] = 0;
        ctrls->cost_th_dist[LPD1_LVL_0] = 128;
        ctrls->coeff_th[LPD1_LVL_0] = 50;
        ctrls->max_mv_length[LPD1_LVL_0] = 300;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = 250000;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0] = 1024;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0] = 1;
        break;
    case 2:
        ctrls->pd1_level = LPD1_LVL_0;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0] = 1;
        ctrls->use_ref_info[LPD1_LVL_0] = 0;
        ctrls->cost_th_dist[LPD1_LVL_0] = 256;
        ctrls->coeff_th[LPD1_LVL_0] = 100;
        ctrls->max_mv_length[LPD1_LVL_0] = 500;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0] = 1024;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0] = 1;
        break;
    case 3:
        ctrls->pd1_level = LPD1_LVL_2;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0] = 1;
        ctrls->use_ref_info[LPD1_LVL_0] = 0;
        ctrls->cost_th_dist[LPD1_LVL_0] = 256 << 9;
        ctrls->coeff_th[LPD1_LVL_0] = 8192;
        ctrls->max_mv_length[LPD1_LVL_0] = 2048;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0] = 16384;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0] = 3;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1] = 1;
        ctrls->use_ref_info[LPD1_LVL_1] = 0;
        ctrls->cost_th_dist[LPD1_LVL_1] = 256 << 9;
        ctrls->coeff_th[LPD1_LVL_1] = 8192;
        ctrls->max_mv_length[LPD1_LVL_1] = 2048;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1] = 16384;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1] = 3;

        // Set LPD1 level 2 controls
        ctrls->use_lpd1_detector[LPD1_LVL_2] = 1;
        ctrls->use_ref_info[LPD1_LVL_2] = 1;
        ctrls->cost_th_dist[LPD1_LVL_2] = 256 << 6;
        ctrls->coeff_th[LPD1_LVL_2] = 2000;
        ctrls->max_mv_length[LPD1_LVL_2] = 1600;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_2] = 500000;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_2] = 16384;
        ctrls->skip_pd0_me_shift[LPD1_LVL_2] = 2;
        break;
    case 4:
        ctrls->pd1_level = LPD1_LVL_4;

        // LPD1 level 3 doesn't use the detector (will be used for all SBs)
        ctrls->use_lpd1_detector[LPD1_LVL_3] = 0;

        // Set LPD1 level 4 controls
        ctrls->use_lpd1_detector[LPD1_LVL_4] = 1;
        ctrls->use_ref_info[LPD1_LVL_4] = 1;
        ctrls->cost_th_dist[LPD1_LVL_4] = 256 << 9;
        ctrls->coeff_th[LPD1_LVL_4] = 8192;
        ctrls->max_mv_length[LPD1_LVL_4] = 2048;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_4] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_4] = 16384 * 2;
        ctrls->skip_pd0_me_shift[LPD1_LVL_4] = 3;
        break;
    default:
        assert(0);
        break;
    }
}
#else
void set_lpd1_ctrls(ModeDecisionContext *ctx, uint8_t lpd1_lvl) {
    Lpd1Ctrls* ctrls = &ctx->lpd1_ctrls;
    switch (lpd1_lvl) {
    case 0:
        ctrls->pd1_level = REGULAR_PD1; // Light-PD1 paths not used
        break;
    case 1:
        ctrls->pd1_level = LIGHT_PD1;

        // Set light-pd1 controls
        ctrls->use_lpd1_detector[LIGHT_PD1] = 1;
        ctrls->use_ref_info[LIGHT_PD1] = 0;
        ctrls->cost_th_dist[LIGHT_PD1] = 128;
        ctrls->coeff_th[LIGHT_PD1] = 50;
        ctrls->max_mv_length[LIGHT_PD1] = 300;
        ctrls->me_8x8_cost_variance_th[LIGHT_PD1] = 250000;
        ctrls->skip_pd0_edge_dist_th[LIGHT_PD1] = 1024;
        ctrls->skip_pd0_me_shift[LIGHT_PD1] = 1;
        break;
    case 2:
        ctrls->pd1_level = LIGHT_PD1;

        // Set light-pd1 controls
        ctrls->use_lpd1_detector[LIGHT_PD1] = 1;
        ctrls->use_ref_info[LIGHT_PD1] = 0;
        ctrls->cost_th_dist[LIGHT_PD1] = 256;
        ctrls->coeff_th[LIGHT_PD1] = 100;
        ctrls->max_mv_length[LIGHT_PD1] = 500;
        ctrls->me_8x8_cost_variance_th[LIGHT_PD1] = (uint32_t) ~0;
        ctrls->skip_pd0_edge_dist_th[LIGHT_PD1] = 1024;
        ctrls->skip_pd0_me_shift[LIGHT_PD1] = 1;
        break;
    case 3:
        ctrls->pd1_level = LIGHT_PD1;

        // Set light-pd1 controls
        ctrls->use_lpd1_detector[LIGHT_PD1] = 1;
        ctrls->use_ref_info[LIGHT_PD1] = 0;
        ctrls->cost_th_dist[LIGHT_PD1] = 256 * 4;
        ctrls->coeff_th[LIGHT_PD1] = 100 * 8;
        ctrls->max_mv_length[LIGHT_PD1] = 800;
        ctrls->me_8x8_cost_variance_th[LIGHT_PD1] = 250000;
        ctrls->skip_pd0_edge_dist_th[LIGHT_PD1] = 1024;
        ctrls->skip_pd0_me_shift[LIGHT_PD1] = 1;
        break;
    case 4:
        ctrls->pd1_level = LIGHT_PD1;

        // Set light-pd1 controls
        ctrls->use_lpd1_detector[LIGHT_PD1] = 1;
        ctrls->use_ref_info[LIGHT_PD1] = 0;
        ctrls->cost_th_dist[LIGHT_PD1] = 256 << 9;
        ctrls->coeff_th[LIGHT_PD1] = 8192;
        ctrls->max_mv_length[LIGHT_PD1] = 2048;
        ctrls->me_8x8_cost_variance_th[LIGHT_PD1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LIGHT_PD1] = 16384;
        ctrls->skip_pd0_me_shift[LIGHT_PD1] = 3;
        break;
#if FTR_M13
    case 5:
        ctrls->pd1_level = LIGHT_PD1;

        // Set light-pd1 controls
        ctrls->use_lpd1_detector[LIGHT_PD1] = 0;
        ctrls->use_ref_info[LIGHT_PD1] = 0;
        ctrls->cost_th_dist[LIGHT_PD1] = 256 << 9;
        ctrls->coeff_th[LIGHT_PD1] = 8192;
        ctrls->max_mv_length[LIGHT_PD1] = 2048;
        ctrls->me_8x8_cost_variance_th[LIGHT_PD1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LIGHT_PD1] = 16384;
        ctrls->skip_pd0_me_shift[LIGHT_PD1] = 3;
        break;
    case 6:
#else
    case 5:
#endif
        ctrls->pd1_level = VERY_LIGHT_PD1;

        // Set light-pd1 controls
        ctrls->use_lpd1_detector[LIGHT_PD1] = 1;
        ctrls->use_ref_info[LIGHT_PD1] = 0;
        ctrls->cost_th_dist[LIGHT_PD1] = 256 << 9;
        ctrls->coeff_th[LIGHT_PD1] = 8192;
        ctrls->max_mv_length[LIGHT_PD1] = 2048;
        ctrls->me_8x8_cost_variance_th[LIGHT_PD1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LIGHT_PD1] = 16384;
        ctrls->skip_pd0_me_shift[LIGHT_PD1] = 3;

        // set very-light-pd1 controls
        ctrls->use_lpd1_detector[VERY_LIGHT_PD1] = 1;
        ctrls->use_ref_info[VERY_LIGHT_PD1] = 1;
        ctrls->cost_th_dist[VERY_LIGHT_PD1] = 256 << 6;
        ctrls->coeff_th[VERY_LIGHT_PD1] = 2000;
        ctrls->max_mv_length[VERY_LIGHT_PD1] = 1600;
        ctrls->me_8x8_cost_variance_th[VERY_LIGHT_PD1] = 500000;
        ctrls->skip_pd0_edge_dist_th[VERY_LIGHT_PD1] = 16384;
        ctrls->skip_pd0_me_shift[VERY_LIGHT_PD1] = 2;
        break;
#if FTR_M13
    case 7:
        ctrls->pd1_level = VERY_LIGHT_PD1;

        // Set light-pd1 controls
        ctrls->use_lpd1_detector[LIGHT_PD1] = 0;
        ctrls->use_ref_info[LIGHT_PD1] = 0;
        ctrls->cost_th_dist[LIGHT_PD1] = 256 << 9;
        ctrls->coeff_th[LIGHT_PD1] = 8192;
        ctrls->max_mv_length[LIGHT_PD1] = 2048;
        ctrls->me_8x8_cost_variance_th[LIGHT_PD1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LIGHT_PD1] = 16384;
        ctrls->skip_pd0_me_shift[LIGHT_PD1] = 3;

        // set very-light-pd1 controls
        ctrls->use_lpd1_detector[VERY_LIGHT_PD1] = 1;
        ctrls->use_ref_info[VERY_LIGHT_PD1] = 1;
        ctrls->cost_th_dist[VERY_LIGHT_PD1] = 256 << 6;
        ctrls->coeff_th[VERY_LIGHT_PD1] = 2000;
        ctrls->max_mv_length[VERY_LIGHT_PD1] = 1600;
        ctrls->me_8x8_cost_variance_th[VERY_LIGHT_PD1] = 500000;
        ctrls->skip_pd0_edge_dist_th[VERY_LIGHT_PD1] = 16384;
        ctrls->skip_pd0_me_shift[VERY_LIGHT_PD1] = 2;
        break;

#endif
    default:
        assert(0);
        break;
    }
}
#endif
#else
void set_lpd1_ctrls(ModeDecisionContext *ctx, uint8_t lpd1_lvl) {
    Lpd1Ctrls* ctrls = &ctx->lpd1_ctrls;
    switch (lpd1_lvl) {
    case 0:
        ctrls->enabled = 0;
        ctrls->use_light_pd1 = 0;
        break;
    case 1:
        ctrls->enabled = 1;
        ctrls->use_light_pd1 = 1;
        ctrls->use_lpd1_detector = 1;
        ctrls->cost_th_dist = 128;
        ctrls->coeff_th = 50;
#if OPT_USE_MVS_LPD1_DETECTOR
        ctrls->max_mv_length = 300;
#endif
#if OPT_LIGHT_PD1_USE_ME_DIST_VAR
#if TUNE_M7_M8_3
        ctrls->me_8x8_cost_variance_th = 250000;
#else
        ctrls->me_8x8_cost_variance_th = (uint32_t) ~0;
#endif
#endif
#if TUNE_LPD1_DETECTOR
        ctrls->skip_pd0_edge_dist_th = 1024;
        ctrls->skip_pd0_me_dist_shift = 1;
#endif
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->use_light_pd1 = 1;
        ctrls->use_lpd1_detector = 1;
        ctrls->cost_th_dist = 256;
        ctrls->coeff_th = 100;
#if OPT_USE_MVS_LPD1_DETECTOR
        ctrls->max_mv_length = 500;
#endif
#if OPT_LIGHT_PD1_USE_ME_DIST_VAR
        ctrls->me_8x8_cost_variance_th = (uint32_t) ~0;
#endif
#if TUNE_LPD1_DETECTOR
        ctrls->skip_pd0_edge_dist_th = 1024;
        ctrls->skip_pd0_me_dist_shift = 1;
#endif
        break;
    case 3:
#if OPT_LIGHT_PD1
        ctrls->enabled = 1;
        ctrls->use_light_pd1 = 1;
        ctrls->use_lpd1_detector = 1;
        ctrls->cost_th_dist =256 * 4;
        ctrls->coeff_th = 100 * 8;
#if OPT_USE_MVS_LPD1_DETECTOR
#if TUNE_4K_M11
        ctrls->max_mv_length = 800;
#else
        ctrls->max_mv_length = 500;
#endif
#endif
#if OPT_LIGHT_PD1_USE_ME_DIST_VAR
#if TUNE_4K_M11
        ctrls->me_8x8_cost_variance_th = 250000;
#else
        ctrls->me_8x8_cost_variance_th = (uint32_t) ~0;
#endif
#endif
#else
        ctrls->enabled = 1;
        ctrls->use_light_pd1 = 1;
        ctrls->use_lpd1_detector = 0;
        ctrls->cost_th_dist = 256;
        ctrls->coeff_th = 100;
#if OPT_USE_MVS_LPD1_DETECTOR
        ctrls->max_mv_length = 500;
#endif
#endif
#if TUNE_LPD1_DETECTOR
        ctrls->skip_pd0_edge_dist_th = 1024;
        ctrls->skip_pd0_me_dist_shift = 1;
#endif
        break;
#if TUNE_LPD1_DETECTOR
    case 4:
        ctrls->enabled = 1;
        ctrls->use_light_pd1 = 1;
        ctrls->use_lpd1_detector = 1;
#if TUNE_LPD1_DETECTOR_LVL
        ctrls->cost_th_dist = 256 << 9;
        ctrls->coeff_th = 8192;
        ctrls->max_mv_length = 2048;
        ctrls->me_8x8_cost_variance_th = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th = 16384;
        ctrls->skip_pd0_me_dist_shift = 3;
#else
        ctrls->cost_th_dist = 256 * 16;
        ctrls->coeff_th = 1025;
        ctrls->max_mv_length = 1600;
        ctrls->me_8x8_cost_variance_th = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th = 4096;
        ctrls->skip_pd0_me_dist_shift = 2;
#endif
        break;
#endif
    default:
        assert(0);
        break;
    }
}
#endif
#endif
#if RFCT_ME8X8
// use this function to set the disallow_below_16x16 level and to set the accompanying enable_me_8x8 level
uint8_t get_disallow_below_16x16_picture_level(EbEncMode enc_mode, EbInputResolution resolution, EB_SLICE slice_type, uint8_t sc_class1) {

    uint8_t disallow_below_16x16 = 0;

    if (sc_class1)
        disallow_below_16x16 = 0;
    else if (enc_mode <= ENC_M7)
        disallow_below_16x16 = 0;
    else if (enc_mode <= ENC_M11)
        disallow_below_16x16 = (resolution <= INPUT_SIZE_720p_RANGE) ? 0 : ((slice_type == I_SLICE) ? 0 : 1);
    else
        disallow_below_16x16 = (slice_type == I_SLICE) ? 0 : 1;

    return disallow_below_16x16;
}
#endif
/*
 * Generate per-SB MD settings (do not change per-PD)
 */
EbErrorType signal_derivation_enc_dec_kernel_common(
    SequenceControlSet *scs_ptr,
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *ctx) {

    EbErrorType return_error = EB_ErrorNone;

    EbEncMode enc_mode = pcs_ptr->enc_mode;

#if !CLN_4X4_SIG
#if OPT_PD0_PATH
    // Set disallow_4x4
    ctx->disallow_4x4 = get_disallow_4x4(enc_mode, pcs_ptr->slice_type);
#if !CLN_REMOVE_UNUSED_FEATS
#if TUNE_M10_M0
#if TUNE_M9_M10
    if (enc_mode <= ENC_M10)
#else
    if (enc_mode <= ENC_M9)
#endif
#else
    if (enc_mode <= ENC_M8)
#endif
        ctx->ep_use_md_skip_decision = 0;
    else
        ctx->ep_use_md_skip_decision = 1;
#endif
#endif
#endif
    // Level 0: pred depth only
    // Level 1: [-2, +2] depth refinement
    // Level 2: [-1, +1] depth refinement
    uint8_t depth_level = 0;
    if (enc_mode <= ENC_MRS)
        depth_level = 1;
    else if (pcs_ptr->parent_pcs_ptr->sc_class1) {
        if (enc_mode <= ENC_M2)
            depth_level = pcs_ptr->slice_type == I_SLICE ? 1 : 2;
        else
            depth_level = 2;
    }
    else if(enc_mode <= ENC_M2)
        depth_level = pcs_ptr->slice_type == I_SLICE ? 1 : 2;
#if FTR_BYPASS_ENCDEC
#if TUNE_M9_SLOW
#if TUNE_M8_M10 && !TUNE_M9_SLOWDOWN
    else if (enc_mode <= ENC_M8)
#else
    else if (enc_mode <= ENC_M9)
#endif
#else
    else if (enc_mode <= ENC_M8)
#endif
        depth_level = 2;
#if TUNE_M8_M10 && !TUNE_M9_M10
    else if (enc_mode <= ENC_M9)
        depth_level = pcs_ptr->slice_type == I_SLICE ? 2 : 0;
#endif
    else
        depth_level = 0;
#else
    else
        depth_level = 2;

#endif
    set_depth_ctrls(ctx, depth_level);
#if FTR_BYPASS_ENCDEC
    ctx->pred_depth_only = (depth_level == 0);
#endif
#if !CLN_BYP_ED_SIG
#if FTR_BYPASS_ENCDEC
#if FIX_10BIT_R2R
    // TODO: Bypassing EncDec doesn't work if HVA_HVB_HV4 are enabled (for all bit depths; causes non-conformant bitstreams),
    // or if NSQ is enabled for 10bit content (causes r2r).
    // TODO: This signal can only be modified per picture right now, not per SB.  Per SB requires
    // neighbour array updates at EncDec for all SBs, that are currently skipped if EncDec is bypassed.
#else
    // TODO: for now, bypassing doesn't work if HVA_HVB_HV4 are enabled.  Only supported for 8bit
    // TODO: This signal can only be modified per picture right now, not per SB.  Per SB requires
    // neighbour array updates at EncDec for all SBs, that are currently skipped if EncDec is bypassed.
#endif
#if FTR_BYPASS_ENCDEC
#if FTR_10BIT_MDS3_LPD1
#if FIX_10BIT_R2R
    if (pcs_ptr->parent_pcs_ptr->disallow_HVA_HVB_HV4 &&
        (scs_ptr->static_config.encoder_bit_depth == EB_8BIT || pcs_ptr->parent_pcs_ptr->disallow_nsq) &&
#else
    if (pcs_ptr->parent_pcs_ptr->disallow_HVA_HVB_HV4 &&
#endif
#else
    if (scs_ptr->static_config.encoder_bit_depth == EB_8BIT && pcs_ptr->parent_pcs_ptr->disallow_HVA_HVB_HV4 &&
#endif
        !pcs_ptr->parent_pcs_ptr->frm_hdr.segmentation_params.segmentation_enabled) {
#if FTR_10BIT_MDS3_LPD1
#if TUNE_BYPASS_MEM
        ctx->bypass_encdec =  get_bypass_encdec( enc_mode , ctx->hbd_mode_decision , scs_ptr->static_config.encoder_bit_depth ) ;
#else
        if (scs_ptr->static_config.encoder_bit_depth == EB_8BIT) {
            // 8bit settings
            if (enc_mode <= ENC_M4)
                ctx->bypass_encdec = 0;
            else
                ctx->bypass_encdec = 1;
        }
        else {
            // 10bit settings
#if TUNE_M7_M8_3
            if (enc_mode <= ENC_M6)
#else
            if (enc_mode <= ENC_M7)
#endif
                ctx->bypass_encdec = 0;
#if FIX_10BIT_R2R
#if TUNE_M8_M10_4K_SUPER
            // To fix max BDR issues with clips like cosmos_aom_sdr_12149-12330_588x250
            else if (enc_mode <= ENC_M10)
#else
            else if (enc_mode <= ENC_M9)
#endif
#else
            if (enc_mode <= ENC_M9)
#endif
                ctx->bypass_encdec = ctx->hbd_mode_decision ? 1 : 0;
            else
                ctx->bypass_encdec = 1;
        }
#endif
#else
#if TUNE_MEGA_M9_M4
#if TUNE_M10_M3_1
        if (enc_mode <= ENC_M4)
#else
        if (enc_mode <= ENC_M5)
#endif
#else
        if (enc_mode <= ENC_M8)
#endif
            ctx->bypass_encdec = 0;
        else
            ctx->bypass_encdec = 1;
#endif
    }
    else
        ctx->bypass_encdec = 0;
#else
    ctx->bypass_encdec = 0;// (scs_ptr->static_config.encoder_bit_depth == EB_8BIT && pcs_ptr->parent_pcs_ptr->disallow_HVA_HVB_HV4) ? 1 : 0;
#endif
#endif
#endif
#if CLN_LPD0_CTRL
#if CLN_PD0_SIG
    ctx->pd0_level = pcs_ptr->pic_pd0_level;
#else
#if TUNE_M1_M8
    if (enc_mode <= ENC_M4)
#else
    if (enc_mode <= ENC_M5)
#endif
        ctx->pd0_level = REGULAR_PD0;
#if TUNE_M1_M8
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M8)
#endif
        ctx->pd0_level = LIGHT_PD0_LVL1;
    else if (enc_mode <= ENC_M9)
        ctx->pd0_level = LIGHT_PD0_LVL2;
#if !TUNE_IMPROVE_M11_M10
    else if (enc_mode <= ENC_M10)
#if TUNE_4K_STABILITY
        ctx->pd0_level = scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? LIGHT_PD0_LVL3 : LIGHT_PD0_LVL2;
#else
        ctx->pd0_level = LIGHT_PD0_LVL3;
#endif
#endif
    else
#if TUNE_4K_STABILITY
        ctx->pd0_level = scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE? (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? LIGHT_PD0_LVL4 : VERY_LIGHT_PD0) : LIGHT_PD0_LVL2;
#else
        ctx->pd0_level = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? LIGHT_PD0_LVL4 : VERY_LIGHT_PD0;
#endif
#endif
#else
#if LIGHT_PD0
#if FTR_VLPD0
#if OPT_M12_SPEED
    if (enc_mode <= ENC_M5)
        ctx->pd0_level = 0;
#if TUNE_M10_M11
    else if (enc_mode <= ENC_M10)
#else
    else if (enc_mode <= ENC_M11)
#endif
        ctx->pd0_level = 1;
#if FTR_M13
    else
        ctx->pd0_level = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 1 : 2;
#else
    else if (enc_mode <= ENC_M12)
        ctx->pd0_level = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 1 : 2;
    else
        ctx->pd0_level = 2;
#endif
#else
    if (enc_mode <= ENC_M5)
        ctx->pd0_level = 0;
    else if (enc_mode <= ENC_M10)
        ctx->pd0_level = 1;
    else
        ctx->pd0_level = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 1 : 2;
#endif
#else
    // The light PD0 path assumes many features are OFF
#if TUNE_MEGA_M9_M4
    if (enc_mode <= ENC_M5)
#else
    if (enc_mode <= ENC_M7)
#endif
        ctx->use_light_pd0 = 0;
    else
        ctx->use_light_pd0 = 1;
#endif
#if 0//LIGHT_PD0_2 // TODO: PHX - should we force features for light-PD0 path?
    // If using light-PD0 path, must force some features ON/OFF
    if (ctx->use_light_pd0) {
        // Use rate est. opt to avoid cost mismatches between check_curr_to_parent_cost_light_pd0() and d2_inter_depth_block_decision()
        pcs_ptr->parent_pcs_ptr->partition_contexts = 4;
    }
#endif
#endif
#endif
#if !CLN_SKIP_PD0_SIG
#if FTR_M13
    if (enc_mode <= ENC_M12)
        ctx->skip_pd0 = 0;
    else
        ctx->skip_pd0 = 1;
#endif
#endif
#if SHUT_8x8_IF_NON_ISLICE
    SbParams *sb_params = &pcs_ptr->parent_pcs_ptr->sb_params_array[ctx->sb_index];
#endif
    ctx->depth_removal_ctrls.disallow_below_64x64 = 0;
    ctx->depth_removal_ctrls.disallow_below_32x32 = 0;
#if ME_8X8
/*
if disallow_below_16x16 is turned ON then enable_me_8x8 should be turned OFF for the same preset in order to save memory and cycles as that feature optimizes the me_candidate_array,
me_mv_array and the total_me_candidate_index arrays when 8x8 blocks are not used

if any check other than an I-SLICE check is used on disallow_below_16x16 then the enable_me_8x8 should be turned ON for the entire preset because without the 8x8 me data the non I-SLICE pictures
that use 8x8 blocks will lose significant BD-Rate as the parent 16x16 me data will be used for the 8x8 blocks
*/
#endif
#if SHUT_8x8_IF_NON_ISLICE
#if RFCT_ME8X8
#if CLN_DISALLOW_BELOW_16X16_SIG
    ctx->depth_removal_ctrls.disallow_below_16x16 = pcs_ptr->pic_disallow_below_16x16;
#else
    ctx->depth_removal_ctrls.disallow_below_16x16 = get_disallow_below_16x16_picture_level(enc_mode, scs_ptr->input_resolution, pcs_ptr->slice_type, pcs_ptr->parent_pcs_ptr->sc_class1);
#endif
#else
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        ctx->depth_removal_ctrls.disallow_below_16x16 = 0;
#if TUNE_M8_M10
#if TUNE_M7_M8
#if TUNE_M7_MT
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M6)
#endif
#else
    else if (enc_mode <= ENC_M8)
#endif
#else
    else if (enc_mode <= ENC_M10)
#endif
        ctx->depth_removal_ctrls.disallow_below_16x16 = 0;
#if TUNE_M8_M10
#if TUNE_M11_SLOWDOWN
    else if (enc_mode <= ENC_M11)
#else
    else if (enc_mode <= ENC_M10)
#endif
        ctx->depth_removal_ctrls.disallow_below_16x16 = (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 0 : ((pcs_ptr->slice_type == I_SLICE) ? 0 : 1);
#endif
    else
        ctx->depth_removal_ctrls.disallow_below_16x16 = (pcs_ptr->slice_type == I_SLICE) ? 0 : 1;
#endif

    if (sb_params->width % 32 != 0 || sb_params->height % 32 != 0)
        ctx->depth_removal_ctrls.disallow_below_64x64 = EB_FALSE;
    if (sb_params->width % 16 != 0 || sb_params->height % 16 != 0)
        ctx->depth_removal_ctrls.disallow_below_32x32 = EB_FALSE;
    if (sb_params->width % 8 != 0 || sb_params->height % 8 != 0)
        ctx->depth_removal_ctrls.disallow_below_16x16 = EB_FALSE;
#else
    ctx->depth_removal_ctrls.disallow_below_16x16 = 0;
#endif

    // me_distortion/variance generated for 64x64 blocks only
#if CLN_DEPTH_REMOVAL_SIG
    if (scs_ptr->static_config.super_block_size == 64) {
        set_depth_removal_level_controls(pcs_ptr, ctx, pcs_ptr->pic_depth_removal_level);
    }
#else
#if OPT_DEPTH_REMOVAL_I_SLICE
    if (scs_ptr->static_config.super_block_size == 64) {
        uint8_t depth_removal_level;
        if (pcs_ptr->slice_type == I_SLICE) {
            if (pcs_ptr->parent_pcs_ptr->sc_class1)
                depth_removal_level = 0;
#if TUNE_M8_M10
#if TUNE_4K_M8_M11
#if CLN_M6_M12_FEATURES
#if TUNE_IMPROVE_M12
            else if (enc_mode <= ENC_M12)
#else
            else if (enc_mode <= ENC_M11)
#endif
#else
            else if (enc_mode <= ENC_M7)
#endif
#else
            else if (enc_mode <= ENC_M8)
#endif
#else
            else if (enc_mode <= ENC_M10)
#endif
                depth_removal_level = 0;
#if TUNE_4K_M8_M11 && !TUNE_M8_M10_4K_SUPER
            else if (enc_mode <= ENC_M8) {
                if (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                    depth_removal_level = 0;
                else
                    depth_removal_level = 1;
            }
#endif
            else
                depth_removal_level = 1;
        }
        else {
#else
    if (pcs_ptr->slice_type != I_SLICE && scs_ptr->static_config.super_block_size == 64) {
#endif
#if FTR_DEPTH_REMOVAL_QP
        // Set depth_removal_level_controls
#if !OPT_DEPTH_REMOVAL_I_SLICE
        uint8_t depth_removal_level;
#endif
        if (pcs_ptr->parent_pcs_ptr->sc_class1)
            depth_removal_level = 0;
        else if (enc_mode <= ENC_M2)
            depth_removal_level = 0;
#if TUNE_M10_M0
#if TUNE_M1_M8
        else if (enc_mode <= ENC_M5) {
#else
        else if (enc_mode <= ENC_M4) {
#endif
#else
        else if (enc_mode <= ENC_M3) {
#endif
            if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = 0;
            else
                depth_removal_level = 2;
        }
#if TUNE_M7_M10_MT && !TUNE_M7_M8 || TUNE_M7_MT
#if TUNE_M10_M3_1 && !TUNE_M8_M10
        else if (enc_mode <= ENC_M8) {
#else
#if TUNE_M7_M8_3
        else if (enc_mode <= ENC_M6) {
#else
        else if (enc_mode <= ENC_M7) {
#endif
#endif
#else
        else if (enc_mode <= ENC_M6) {
#endif
#if FTR_SIMPLIFIED_DEPTH_REMOVAL
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = 1;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = 1;
            else
                depth_removal_level = 2;
#endif
        }
#if TUNE_M7_M8_3
#if !TUNE_M1_M8
        else if (enc_mode <= ENC_M7) {
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = 1;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 5;
        }
#endif
#endif
#if !TUNE_M10_M3_1
        else if (enc_mode <= ENC_M8) {
#if FTR_SIMPLIFIED_DEPTH_REMOVAL
            if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 1 : 2;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = 1;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = 2;
            else
                depth_removal_level = 3;
#endif
        }
#endif
#if TUNE_4K_M8_M11 && !TUNE_M9_SLOWDOWN
        else if (enc_mode <= ENC_M8) {
#else
#if CLN_DEPTH_REMOVAL
        else if (enc_mode <= ENC_M11) {
#else
        else if (enc_mode <= ENC_M9) {
#endif
#endif
#if FTR_SIMPLIFIED_DEPTH_REMOVAL
#if TUNE_M7_M8_3
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
#if TUNE_M9_11_3
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 3;
#else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
#endif
#if !CLN_DEPTH_REMOVAL
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 6;
#endif
            else if (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 6;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 11;
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 3;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
#endif
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_240p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 3;
            else
                depth_removal_level = 3;
#endif
        }
#if !CLN_DEPTH_REMOVAL
#if TUNE_4K_M8_M11
#if TUNE_M10_SLOWDOWN
#if TUNE_M11_SLOWDOWN
        else if (enc_mode <= ENC_M11) {
#else
        else if (enc_mode <= ENC_M10) {
#endif
#else
        else if (enc_mode <= ENC_M9) {
#endif
#if TUNE_M7_M8_3
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
#if TUNE_M9_11_3
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 5;
#else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
#endif
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 6;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 6;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 11;
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 3;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 9;
#endif
        }
#endif
        else if (enc_mode <= ENC_M10) {
#if FTR_SIMPLIFIED_DEPTH_REMOVAL
#if TUNE_M7_M8_3
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 7;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 9;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 10;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 11;
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 5;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 7;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 8;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 9;
#endif
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 5;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 6;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 7;
#endif
        }
#endif
#if TUNE_NEW_M11_2
        else {
#else
        else if (enc_mode <= ENC_M11) {
#endif
#if TUNE_M11
#if FTR_SIMPLIFIED_DEPTH_REMOVAL
#if TUNE_M7_M8_3
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 8;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 11;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 14;
            else
#if TUNE_4K_STABILITY
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 11;
#else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 11 : 15;
#endif
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 6;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 9;
#if TUNE_4K_M11
            else if (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 12;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 9 : 13;
#else
            else
#if TUNE_M11_2
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 12;
#else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 10;
#endif
#endif
#endif
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 7;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 6;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 7;
            else
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 12;
#endif
#else
            if (scs_ptr->input_resolution <= INPUT_SIZE_240p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 9;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 7 : 10;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 9 : 10;
            else if (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 9 : 12;
            else
                depth_removal_level = 12;
#endif
        }
#if !TUNE_NEW_M11_2
        else {
            depth_removal_level = 12;
        }
#endif
#else
#if CLN_FTR_EARLY_DEPTH_REMOVAL
        // Set depth_removal_level_controls
        uint8_t depth_removal_level;
#endif
    // Set depth_removal_level_controls
        uint8_t depth_removal_level;
        if (pcs_ptr->parent_pcs_ptr->sc_class1)
            depth_removal_level = 0;
        else
#if TUNE_M3_M6_MEM_OPT
            if (enc_mode <= ENC_M2)
                depth_removal_level = 0;
            else if (enc_mode <= ENC_M3) {
#else
            if (enc_mode <= ENC_M4)
#endif
#if TUNE_M3_M6_MEM_OPT
                if (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                    depth_removal_level = 0;
                else
                    depth_removal_level = 2;
            }
#else
            depth_removal_level = 0;
#endif
        else
#if TUNE_MEGA_M9_M4
            if (enc_mode <= ENC_M6) {
                    if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                        depth_removal_level = 1;
                    else
                        depth_removal_level = 2;
            }
            else
#endif
                if (enc_mode <= ENC_M7) {
                if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                    depth_removal_level = 1;
                else if (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE)
                    depth_removal_level = 2;
                else if (scs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE)
#if TUNE_MEGA_M9_M4
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
#else
                        depth_removal_level = 2;
#endif
                else
#if TUNE_M0_M7_MEGA_FEB
#if TUNE_MEGA_M9_M4
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 6;
#else
                        depth_removal_level = 2;
#endif
#else
                    depth_removal_level = 2;
#endif
            }
            else if (enc_mode <= ENC_M8) {
#if TUNE_MEGA_M9_M4
#if TUNE_M10_M7
                    if (scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE)
                        depth_removal_level = 1;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                        depth_removal_level = 2;
#else
                    if (scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 3;
#endif
                    else if (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
                    else
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 6;
#else
                    if (scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 3;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 4;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 7;
                else if (scs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE)
                    depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 7;
                else
                    depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 7 : 8;
#endif
            }
#if TUNE_M10_DEPTH_ME
                else if (enc_mode <= ENC_M9) {
#else
                else {
#endif
#if TUNE_MEGA_M9_M4
#if TUNE_M10_M7
                if (scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE)
                    depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
                else if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                    depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
                else if (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE)
                    depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 6;
                else if (scs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE)
                    depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 6;
                else
                    depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 6;
#else
                    if (scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 4;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 6;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 7 : 8;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 7 : 10;
                    else
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 8 : 11;
#endif
#else
                    if (scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 4;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 6;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 7 : 8;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 9 : 10;
                    else
                        depth_removal_level = 11;
#endif
                }
#if TUNE_M10_DEPTH_ME
#if TUNE_NEW_M10_M11
                else if (enc_mode <= ENC_M11) {
#else
                else {
#endif
                    if (scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 9;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 7 : 10;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 9 : 10;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE)
                        depth_removal_level = (pcs_ptr->temporal_layer_index == 0) ? 9 : 12;
                    else
                        depth_removal_level = 12;
                }
#endif
#if TUNE_NEW_M10_M11
                else {
                    if (scs_ptr->input_resolution < INPUT_SIZE_360p_RANGE)
                        depth_removal_level = 12;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_480p_RANGE)
                        depth_removal_level = 12;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_720p_RANGE)
                        depth_removal_level = 12;
                    else if (scs_ptr->input_resolution < INPUT_SIZE_1080p_RANGE)
                        depth_removal_level = 12;
                    else
                        depth_removal_level = 12;
                }
#endif

#endif
#if OPT_DEPTH_REMOVAL_I_SLICE
        }
#endif
        set_depth_removal_level_controls(pcs_ptr, ctx, depth_removal_level);
    }
#endif
#if FTR_M13
    if (/*scs_ptr->rc_stat_gen_pass_mode || */ctx->skip_pd0) {
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
#endif
#if !CLN_BLOCK_BASED_DEPTH_SIG
#if OPTIMIZE_L6
    // Set block_based_depth_refinement_level
    int8_t block_based_depth_refinement_level;
#else
    // Set block_based_depth_refinement_level
    uint8_t block_based_depth_refinement_level;
#endif
#if OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS
    // do not use feature for SC
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        block_based_depth_refinement_level = 0;
    else if (enc_mode <= ENC_M2)
        block_based_depth_refinement_level = 0;
#if OPT_DEPTH_REFINEMENT_M3_M4
#if TUNE_M0_M7_MEGA_FEB
    else if (enc_mode <= ENC_M4)
#else
    else if (enc_mode <= ENC_M3)
#endif
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 0 : 2;
    else if (enc_mode <= ENC_M5)
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
#if TUNE_M10_M3_1 && !TUNE_M7_M8
    else if (enc_mode <= ENC_M7)
#else
#if TUNE_M7_M8_2
#if TUNE_M8_SLOWDOWN
    else if (enc_mode <= ENC_M8)
#else
    else if (enc_mode <= ENC_M7)
#endif
#else
    else if (enc_mode <= ENC_M6)
#endif
#endif
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
#else
    else if (enc_mode <= ENC_M3)
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 0 : 3;
    else if (enc_mode <= ENC_M4)
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 3;
#endif
#if !TUNE_M10_M3_1
    else if (enc_mode <= ENC_M7)
        block_based_depth_refinement_level = 4;
#endif
#if TUNE_TXS_IFS_MFMV_DEPTH_M9
#if !TUNE_M7_M8_2
#if TUNE_M7_M8
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M8)
#endif
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 8;
#endif
#if !TUNE_M8_SLOWDOWN
#if TUNE_M7_M8
    else if (enc_mode <= ENC_M8)
 #if OPTIMIZE_L6
        block_based_depth_refinement_level = (scs_ptr->static_config.hierarchical_levels == (EB_MAX_TEMPORAL_LAYERS - 1)) ?
        ((pcs_ptr->temporal_layer_index == 0) ? 6 : 8):
        ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? (pcs_ptr->temporal_layer_index == 0) ? 6 : 10 : (pcs_ptr->temporal_layer_index == 0) ? 6 : 11);

#else
        block_based_depth_refinement_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? (pcs_ptr->temporal_layer_index == 0) ? 6 : 10 : (pcs_ptr->temporal_layer_index == 0) ? 6 : 11;
#endif
#endif
#endif
#if TUNE_M10_M0
#if CLN_M10_M12_DIFFS
    else
#else
    else if (enc_mode <= ENC_M10)
#endif
#else
    else if (enc_mode <= ENC_M9)
#endif
#if TUNE_M8_M9_FEB24
        block_based_depth_refinement_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? (pcs_ptr->temporal_layer_index == 0) ? 6 : 10 : (pcs_ptr->slice_type == I_SLICE) ? 6 : 11;
#else
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 9;
#endif
#if !TUNE_M10_M0
    else if (enc_mode <= ENC_M10)
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 9 : 10;
#endif
#if !CLN_M10_M12_DIFFS
    else
        block_based_depth_refinement_level = 11;
#endif
#else
#if TUNE_NEW_M9_LEVEL
    else if (enc_mode <= ENC_M9)
#else
    else if (enc_mode <= ENC_M8)
#endif
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 8;
    else
        block_based_depth_refinement_level = (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ?
        (pcs_ptr->temporal_layer_index == 0) ? 6 : 9
        : 9;
#endif
#else
    // do not use feature for SC
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        block_based_depth_refinement_level = 0;
    else
    if (enc_mode <= ENC_M2)
        block_based_depth_refinement_level = 0;
    else if (enc_mode <= ENC_M4)
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 0 : 2;
    else if (enc_mode <= ENC_M6)
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
    else if (enc_mode <= ENC_M7)
        block_based_depth_refinement_level = 3;
    else
        block_based_depth_refinement_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 7;
#endif
#if OPTIMIZE_L6
#if FIX_DG
    if (scs_ptr->static_config.max_heirachical_level == (EB_MAX_TEMPORAL_LAYERS - 1))
#else
    if (scs_ptr->static_config.hierarchical_levels == (EB_MAX_TEMPORAL_LAYERS - 1))
#endif
        block_based_depth_refinement_level = MAX(0, block_based_depth_refinement_level - 1);
#endif
    set_block_based_depth_refinement_controls(ctx, block_based_depth_refinement_level);
#endif
#if !CLN_RATE_EST_CTRLS
#if  FTR_SIMPLIFIED_MV_COST
    ctx->use_low_precision_cost_estimation = pcs_ptr->use_low_precision_cost_estimation;
#else
#if FTR_LOW_AC_COST_EST
#if OPT_PD0_TXT
    if(enc_mode <= ENC_M8)
#else
    if(enc_mode <= ENC_M9)
#endif
        ctx->use_low_precision_cost_estimation = 0;
    else
        ctx->use_low_precision_cost_estimation = 1;
#endif
#endif
#endif

#if CLN_LPD1_LVL_SIG
    set_lpd1_ctrls(ctx, pcs_ptr->pic_lpd1_lvl);
#else
#if FTR_LPD1_DETECTOR
    uint8_t lpd1_lvl;
#if CLN_LPD1_LVLS
    if (enc_mode <= ENC_M7)
        lpd1_lvl = 0;
#if !TUNE_M1_M8
    else if (enc_mode <= ENC_M8)
        lpd1_lvl = (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 0 : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1);
#endif
    else if (enc_mode <= ENC_M9)
        lpd1_lvl = (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 0 : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1); // used to be 2
    else if (enc_mode <= ENC_M10)
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 2;
    // Possible intermediate level for M11: lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 2;
    else if (enc_mode <= ENC_M12)
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 3;
    else
#if TUNE_M13_LPD1_LVL
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 4;
#else
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 3 : 4;
#endif
#else
#if TUNE_M9_M10
#if TUNE_M7_M8_3
    if (enc_mode <= ENC_M7)
#else
    if (enc_mode <= ENC_M8)
#endif
        lpd1_lvl = 0;
#if TUNE_M7_M8_3
    else if (enc_mode <= ENC_M8)
        lpd1_lvl = (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 0 : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1);
#endif
    else if (enc_mode <= ENC_M9)
        lpd1_lvl = (scs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 0 : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 2);
    else if (enc_mode <= ENC_M10)
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 2;
#else
    if (enc_mode <= ENC_M10)
        lpd1_lvl = 0;
#endif
#if FTR_M13
    else if (enc_mode <= ENC_M12)
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 6;
    else
#if CLN_LPD1_LVLS
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 6 : 7;
#else
        lpd1_lvl = pcs_ptr->slice_type == I_SLICE ? 5 : 7;
#endif
#else
    else
#if TUNE_4K_M11
#if TUNE_LPD1_DETECTOR
#if FTR_VLPD1
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 5;
#else
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 4;
#endif
#else
        lpd1_lvl = (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 2) : (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 3);
#endif
#else
        lpd1_lvl = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : 2;
#endif
#endif
#endif
    // Can only use light-PD1 under the following conditions
    // There is another check before PD1 is called; pred_depth_only is not checked here, because some modes
    // may force pred_depth_only at the light-pd1 detector
    if (lpd1_lvl &&
        !(ctx->hbd_mode_decision == 0 &&
#if !OPT_LPD1_MRP // NB reference pruning not done in light-PD1
            pcs_ptr->parent_pcs_ptr->frm_hdr.allow_warped_motion == 0 &&
            pcs_ptr->parent_pcs_ptr->ref_list0_count_try == 1 &&
            pcs_ptr->parent_pcs_ptr->ref_list1_count_try == 1 &&
#endif
            //ctx->pred_depth_only &&
            pcs_ptr->parent_pcs_ptr->disallow_nsq == EB_TRUE &&
            ctx->disallow_4x4 == EB_TRUE &&
            scs_ptr->static_config.super_block_size == 64)) {

        lpd1_lvl = 0;
    }

    set_lpd1_ctrls(ctx, lpd1_lvl);
#endif
#endif
    return return_error;
}
/*
 * Generate per-SB/per-PD MD settings
 */
#if !CLN_NIC_SIGS
/*
* return the nic scalling level
  Used by nics control and memory allocation
*/
uint8_t  get_nic_scaling_level(PdPass pd_pass, EbEncMode enc_mode ,uint8_t temporal_layer_index ) {
    uint8_t  nic_scaling_level  = 1 ;
    if (pd_pass == PD_PASS_0)
        nic_scaling_level = 15;
#if !FIX_PRESET_TUNING
    else if (pd_pass == PD_PASS_1)
        nic_scaling_level = 12;
#endif
    else
        if (enc_mode <= ENC_MR)
            nic_scaling_level = 0;
        else if (enc_mode <= ENC_M1)
            nic_scaling_level = (temporal_layer_index == 0) ? 1 : 2;
        else if (enc_mode <= ENC_M3)
            nic_scaling_level = 6;
        else if (enc_mode <= ENC_M4)
            nic_scaling_level = 8;
#if !TUNE_MEGA_M9_M4
        else if (enc_mode <= ENC_M5)
            nic_scaling_level = 10;
#else
        else if (enc_mode <= ENC_M6)
            nic_scaling_level = 11;
#endif
        else if (enc_mode <= ENC_M7)
            nic_scaling_level = 12;
#if TUNE_M10_INITIAL_PRESET && !TUNE_M9_M10_MAR
        else if (enc_mode <= ENC_M10)
#else
#if (TUNE_M10_M0 && !TUNE_M8_M10) || TUNE_M10_SLOWDOWN
#if CLN_M10_M12_DIFFS
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M10)
#endif
#else
        else if (enc_mode <= ENC_M9)
#endif
#endif
            nic_scaling_level = 14;
        else
            nic_scaling_level = 15;
   return nic_scaling_level ;
}
#endif
void set_dist_based_ref_pruning_controls(
    ModeDecisionContext *mdctxt, uint8_t dist_based_ref_pruning_level) {
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

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]            = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP]    = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]      = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]             = 1;

        break;
    case 2:
        ref_pruning_ctrls->enabled = 1;


        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP] = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP] = 90;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP] = 90;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP] = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP] = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP] = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP] = 60;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST] = 60;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF] = 60;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE] = 60;
        ref_pruning_ctrls->ref_idx_2_offset = 10;
        ref_pruning_ctrls->ref_idx_3_offset = 20;


        ref_pruning_ctrls->closest_refs[PA_ME_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST] = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF] = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE] = 1;

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

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]            = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP]    = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]      = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]             = 1;

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

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]            = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP]    = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]      = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]             = 1;

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

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]            = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP]    = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]      = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]             = 1;
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

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]            = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP]    = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]           = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]             = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]      = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]              = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]             = 1;
        break;

    default: assert(0); break;
    }
}

#if FIX_REMOVE_PD1
#if OPT_TXS_SEARCH
void set_txs_controls(ModeDecisionContext *ctx, uint8_t txs_level) {

    TxsControls * txs_ctrls = &ctx->txs_ctrls;

    switch (txs_level)
    {
    case 0:
        txs_ctrls->enabled                 = 0;
        break;

    case 1:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 2;
        txs_ctrls->inter_class_max_depth   = 2;
        txs_ctrls->depth1_txt_group_offset = 0;
        txs_ctrls->depth2_txt_group_offset = 0;
#if OPT_TXS_WM
        txs_ctrls->min_sq_size = 0;
#endif
        break;

    case 2:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 2;
        txs_ctrls->inter_class_max_depth   = 1;
        txs_ctrls->depth1_txt_group_offset = 0;
        txs_ctrls->depth2_txt_group_offset = 0;
#if OPT_TXS_WM
        txs_ctrls->min_sq_size = 0;
#endif
        break;

    case 3:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 2;
        txs_ctrls->inter_class_max_depth   = 0;
        txs_ctrls->depth1_txt_group_offset = 0;
        txs_ctrls->depth2_txt_group_offset = 0;
#if OPT_TXS_WM
        txs_ctrls->min_sq_size = 0;
#endif
        break;

    case 4:
        txs_ctrls->enabled                 = 1;
        txs_ctrls->prev_depth_coeff_exit   = 1;
        txs_ctrls->intra_class_max_depth   = 1;
        txs_ctrls->inter_class_max_depth   = 0;
#if TUNE_TXS_IFS_MFMV_DEPTH_M9
        txs_ctrls->depth1_txt_group_offset = 4;
        txs_ctrls->depth2_txt_group_offset = 4;
#else
        txs_ctrls->depth1_txt_group_offset = 1;
        txs_ctrls->depth2_txt_group_offset = 2;
#endif
#if OPT_TXS_WM
        txs_ctrls->min_sq_size = 0;
#endif
        break;
#if OPT_TXS_WM
    case 5:
        txs_ctrls->enabled = 1;
        txs_ctrls->prev_depth_coeff_exit = 1;
        txs_ctrls->intra_class_max_depth = 1;
        txs_ctrls->inter_class_max_depth = 0;
        txs_ctrls->depth1_txt_group_offset = 4;
        txs_ctrls->depth2_txt_group_offset = 4;
        txs_ctrls->min_sq_size = 32;
        break;
#endif
    default:
        assert(0);
        break;
    }
}
#endif
#if TUNE_BLOCK_SIZE
void set_spatial_sse_full_loop_level(ModeDecisionContext* ctx, uint8_t spatial_sse_full_loop_level) {

    SpatialSSECtrls* spatial_sse_ctrls = &ctx->spatial_sse_ctrls;

    switch (spatial_sse_full_loop_level) {
    case 0:
        spatial_sse_ctrls->spatial_sse_full_loop_level = EB_FALSE;
#if !CLN_FEAT_LEVEL
        spatial_sse_ctrls->blk_size_res_based_modulation = 0;
#endif
        break;
    case 1:
        spatial_sse_ctrls->spatial_sse_full_loop_level = EB_TRUE;
#if CLN_FEAT_LEVEL
        break;
#else
        spatial_sse_ctrls->blk_size_res_based_modulation = 0;
        break;
    case 2:
        spatial_sse_ctrls->spatial_sse_full_loop_level = EB_TRUE;
        spatial_sse_ctrls->blk_size_res_based_modulation = 1;
        break;
#endif
    default:
        assert(0);
        break;
    }
}
#endif
#if FTR_SUBSAMPLE_RESIDUAL
// Compute a qp-aware threshold based on the variance of the SB, used to apply selectively subres
uint64_t compute_subres_th(SequenceControlSet *scs,
    PictureControlSet *pcs,
    ModeDecisionContext *ctx) {

    uint32_t fast_lambda = ctx->hbd_mode_decision ?
        ctx->fast_lambda_md[EB_10_BIT_MD] :
        ctx->fast_lambda_md[EB_8_BIT_MD];
    uint32_t sb_size = scs->static_config.super_block_size * scs->static_config.super_block_size;
    uint64_t cost_th_rate = 1 << 13;
    uint64_t use_subres_th = 0;

    if (pcs->parent_pcs_ptr->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 400)
        use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 8);
    else if (pcs->parent_pcs_ptr->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 800)
        use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 7);
    else
        use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 6);

    return use_subres_th;
}
#endif
#if OPT_REDUCE_TX
// Compute a qp-aware threshold based on the variance of the SB, used to apply selectively apply PF
uint64_t compute_pf_th(SequenceControlSet *scs,
                   PictureControlSet *pcs,
                   ModeDecisionContext *ctx) {

    uint32_t fast_lambda = ctx->hbd_mode_decision ?
        ctx->fast_lambda_md[EB_10_BIT_MD] :
        ctx->fast_lambda_md[EB_8_BIT_MD];
    uint32_t sb_size = scs->static_config.super_block_size * scs->static_config.super_block_size;
    uint64_t cost_th_rate = 1 << 13;
    uint64_t use_pf_th = 0;

    if (pcs->parent_pcs_ptr->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 400)
        use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
    else if (pcs->parent_pcs_ptr->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 800)
        use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
    else
        use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);

    return use_pf_th;
}
#endif
#if !OPT_PD0_PATH
/*
* return the 4x4 level
  Used by signal_derivation_enc_dec_kernel_oq and memory allocation
*/
uint8_t get_disallow_4x4(EbEncMode enc_mode, EB_SLICE slice_type) {
    uint8_t disallow_4x4;
    if (enc_mode <= ENC_M0)
        disallow_4x4 = EB_FALSE;
#if NEW_PRESETS && !TUNE_M10_M7
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M6)
#endif
        disallow_4x4 = (slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
    else
        disallow_4x4 = EB_TRUE;
    return disallow_4x4;
}
#endif
#if FTR_SKIP_MDS1 && !CLN_NIC_PRUNE_CTRLS
void set_mds1_skipping_controls(ModeDecisionContext* ctx, uint8_t mds1_skip_level) {

    SkipMDS1Ctrls* ctrls = &ctx->skip_mds1_ctrls;

    switch (mds1_skip_level) {
    case 0:
        ctrls->enabled = 0;
        break;
#if CLN_REG_PD1_TX_CTRLS
    case 1:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 0;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 20;
        break;
#else
#if CLN_DECPL_TX_FEATS
    case 1:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 0;
        ctrls->use_mds3_shortcuts_th = 0;
#if FTR_TX_NEIGH_INFO
        ctrls->use_neighbour_info = 0;
#endif
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 0;
        ctrls->use_mds3_shortcuts_th = 25;
#if FTR_TX_NEIGH_INFO
        ctrls->use_neighbour_info = 0;
#endif
        break;
#if FTR_TX_NEIGH_INFO
    case 3:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 0;
        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info = 1;
        break;
#if !CLN_LPD1_LVLS
    case 4:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 20;
        ctrls->use_mds3_shortcuts_th = 30;
        ctrls->use_neighbour_info = 1;
        break;
#if FTR_VLPD1
    case 5:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 20;
        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info = 1;
        break;
#endif
#endif
#if FTR_M13
#if CLN_LPD1_LVLS
    case 4:
#else
    case 6:
#endif
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 0;
        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info = 2;
        break;
#if CLN_LPD1_LVLS
    case 5:
#else
    case 7:
#endif
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 20;
        ctrls->use_mds3_shortcuts_th = 100;
        ctrls->use_neighbour_info = 2;
        break;
#endif
#else
    case 3:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 20;
        ctrls->use_mds3_shortcuts_th = 30;
        break;
#endif
#else
    case 1:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 0;
        ctrls->use_mds3_shortcuts_th = 25;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 20;
        ctrls->use_mds3_shortcuts_th = 25;
        break;
    case 3:
        ctrls->enabled = 1;
        ctrls->force_1_cand_th = 20;
        ctrls->use_mds3_shortcuts_th = 30;
        break;
#endif
#endif
    default:
        assert(0);
        break;
    }
}
#endif
#if !CLN_CAND_REDUCTION_CTRLS
#if REMOVE_CLOSE_MVS
void set_redundant_cand_controls(ModeDecisionContext* ctx, uint8_t redundant_cand_level) {

    RedundantCandCtrls *ctrls = &ctx->redundant_cand_ctrls;

    switch (redundant_cand_level) {
    case 0:
        ctrls->score_th = 0;
        break;
    case 1:
        ctrls->score_th = 8;
        ctrls->mag_th = 64;
        break;
#if !CLN_LPD1_LVLS
#if FTR_VLPD1
    case 2:
        ctrls->score_th = 16;
        ctrls->mag_th = 128;
        break;
#endif
#if FTR_LESS_BI
    case 3:
        ctrls->score_th = 32;
        ctrls->mag_th = 32;
        break;
#endif
#endif
    default:
        assert(0);
        break;
    }
}
#endif
#endif
#if !CLN_CAND_REDUCTION_CTRLS
#if OPT_USE_INTRA_NEIGHBORING
void set_use_neighbouring_mode_ctrls(ModeDecisionContext* ctx, uint8_t use_neighbouring_mode) {

    UseNeighbouringModeCtrls* ctrls = &ctx->use_neighbouring_mode_ctrls;

    switch (use_neighbouring_mode) {
    case 0:
        ctrls->enabled = 0;
        break;
    case 1:
        ctrls->enabled = 1;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#endif
#if FTR_SKIP_TX_LPD1
#if CLN_LPD1_TX_CTRLS
void set_lpd1_tx_ctrls(ModeDecisionContext* ctx, uint8_t lpd1_tx_level) {

    Lpd1TxCtrls* ctrls = &ctx->lpd1_tx_ctrls;

    switch (lpd1_tx_level) {
    case 0:
        ctrls->zero_y_coeff_exit = 0;
        ctrls->skip_luma_tx_lvl = SKIP_TX_OFF;
        ctrls->skip_tx_th = 0;
        ctrls->use_uv_shortcuts_on_y_coeffs = 0;

        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->use_neighbour_info = 0;
        break;
    case 1:
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_level = 1;
        ctrls->skip_luma_tx_lvl = SKIP_TX_OFF;
        ctrls->skip_tx_th = 0;
        ctrls->use_skip_tx_neigh_coeff_detector = 1;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info = 0;
        break;
    case 2:
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_level = 1;
        ctrls->skip_luma_tx_lvl = SKIP_NRST_NRST_TX;
        ctrls->skip_tx_th = 25;
        ctrls->use_skip_tx_neigh_coeff_detector = 1;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info = 0;
        break;
    case 3:
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_level = 2;
        ctrls->skip_luma_tx_lvl = SKIP_MVP_TX;
        ctrls->skip_tx_th = 25;
        ctrls->use_skip_tx_neigh_coeff_detector = 1;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info = 0;
        break;
    case 4:
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_level = 2;
        ctrls->skip_luma_tx_lvl = SKIP_INTER_TX;
        ctrls->skip_tx_th = 25;
        ctrls->use_skip_tx_neigh_coeff_detector = 1;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info = 1;
        break;
    case 5:
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_level = 2;
        ctrls->skip_luma_tx_lvl = SKIP_INTER_TX;
        ctrls->skip_tx_th = 50;
        ctrls->use_skip_tx_neigh_coeff_detector = 1;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info = 2;
        break;
    case 6:
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_level = 2;
        ctrls->skip_luma_tx_lvl = SKIP_INTER_TX;
        ctrls->skip_tx_th = 70;
        ctrls->use_skip_tx_neigh_coeff_detector = 0;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 100;
        ctrls->use_neighbour_info = 2;
        break;
    default:
        assert(0);
        break;
    }
}
#else
void set_skip_tx_ctrls(ModeDecisionContext* ctx, uint8_t skip_tx_level) {

    SkipTxCtrls* ctrls = &ctx->skip_tx_ctrls;

    switch (skip_tx_level) {
    case 0:
        ctrls->zero_y_coeff_exit = 0;
        ctrls->skip_nrst_nrst_tx = 0;
#if OPT_TX_SKIP
        ctrls->skip_mvp_tx = 0;
#endif
#if FTR_TX_NEIGH_INFO
        ctrls->skip_tx_th = 0;
#endif
        break;
    case 1:
        ctrls->zero_y_coeff_exit = 1;
#if FTR_VLPD1
        ctrls->chroma_detector_dist_shift = 1;
        ctrls->chroma_detector_downsampling = 1;
#endif
        ctrls->skip_nrst_nrst_tx = 0;
#if OPT_TX_SKIP
        ctrls->skip_mvp_tx = 0;
#endif
#if FTR_TX_NEIGH_INFO
        ctrls->skip_tx_th = 0;
#endif
#if FTR_VLPD1
        ctrls->skip_inter_tx = 0;
        ctrls->use_neigh_coeff_info = 1;
#endif
        break;
    case 2:
        ctrls->zero_y_coeff_exit = 1;
#if FTR_VLPD1
        ctrls->chroma_detector_dist_shift = 1;
        ctrls->chroma_detector_downsampling = 1;
#endif
        ctrls->skip_nrst_nrst_tx = 1;
#if OPT_TX_SKIP
        ctrls->skip_mvp_tx = 0;
#endif
#if FTR_TX_NEIGH_INFO
        ctrls->skip_tx_th = 25;
#endif
#if FTR_VLPD1
        ctrls->skip_inter_tx = 0;
        ctrls->use_neigh_coeff_info = 1;
#endif
        break;
#if OPT_TX_SKIP
    case 3:
        ctrls->zero_y_coeff_exit = 1;
#if FTR_VLPD1
        ctrls->chroma_detector_dist_shift = 2;
        ctrls->chroma_detector_downsampling = 2;
#endif
        ctrls->skip_nrst_nrst_tx = 1;
        ctrls->skip_mvp_tx = 1;
#if FTR_TX_NEIGH_INFO
        ctrls->skip_tx_th = 25;
#endif
#if FTR_VLPD1
        ctrls->skip_inter_tx = 0;
        ctrls->use_neigh_coeff_info = 1;
#endif
        break;
#if FTR_VLPD1
#if !CLN_LPD1_LVLS
    case 4:
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_dist_shift = 2;
        ctrls->chroma_detector_downsampling = 2;
        ctrls->skip_nrst_nrst_tx = 1;
        ctrls->skip_mvp_tx = 1;
        ctrls->skip_tx_th = 50;
        ctrls->skip_inter_tx = 0;
        ctrls->use_neigh_coeff_info = 1;
        break;
#endif
#if FTR_M13
#if CLN_LPD1_LVLS
    case 4:
#else
    case 5:
#endif
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_dist_shift = 2;
        ctrls->chroma_detector_downsampling = 2;
        ctrls->skip_nrst_nrst_tx = 1;
        ctrls->skip_mvp_tx = 1;
        ctrls->skip_tx_th = 25;
        ctrls->skip_inter_tx = 1;
        ctrls->use_neigh_coeff_info = 1;
        break;
#if CLN_LPD1_LVLS
    case 5:
#else
    case 6:
#endif
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_dist_shift = 2;
        ctrls->chroma_detector_downsampling = 2;
        ctrls->skip_nrst_nrst_tx = 1;
        ctrls->skip_mvp_tx = 1;
        ctrls->skip_tx_th = 50;
        ctrls->skip_inter_tx = 1;
        ctrls->use_neigh_coeff_info = 1;
        break;
#if CLN_LPD1_LVLS
    case 6:
#else
    case 7:
#endif
#else
    case 5:
#endif
        ctrls->zero_y_coeff_exit = 1;
        ctrls->chroma_detector_dist_shift = 2;
        ctrls->chroma_detector_downsampling = 2;
        ctrls->skip_nrst_nrst_tx = 1;
        ctrls->skip_mvp_tx = 1;
        ctrls->skip_tx_th = 70;
        ctrls->skip_inter_tx = 1;
        ctrls->use_neigh_coeff_info = 0;
        break;
#endif
#endif
    default:
        assert(0);
        break;
    }
}
#endif
#endif
#if SS_CLN_CFL_CTRLS
void set_cfl_ctrls(ModeDecisionContext* ctx, uint8_t cfl_level) {

    CflCtrls* ctrls = &ctx->cfl_ctrls;

    switch (cfl_level) {
    case 0:
        ctrls->enabled = 0;
        break;
    case 1:
        ctrls->enabled = 1;
        ctrls->itr_th = 2;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->itr_th = 1;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if CLN_RATE_EST_CTRLS
void set_rate_est_ctrls(ModeDecisionContext* ctx, uint8_t rate_est_level) {

    MdRateEstCtrls* ctrls = &ctx->rate_est_ctrls;

    switch (rate_est_level) {
    case 0:
        ctrls->update_skip_ctx_dc_sign_ctx = 0;
        ctrls->update_skip_coeff_ctx = 0;
        ctrls->coeff_rate_est_lvl = 0;
        ctrls->lpd0_qp_offset = 8;
        ctrls->pd0_fast_coeff_est_level = 2;
        break;
    case 1:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx = 1;
        ctrls->coeff_rate_est_lvl = 1;
        ctrls->lpd0_qp_offset = 0;
        ctrls->pd0_fast_coeff_est_level = 1;
        break;
    case 2:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx = 0;
        ctrls->coeff_rate_est_lvl = 1;
        ctrls->lpd0_qp_offset = 0;
        ctrls->pd0_fast_coeff_est_level = 2;
        break;
    case 3:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx = 0;
        ctrls->coeff_rate_est_lvl = 2;
        ctrls->lpd0_qp_offset = 0;
        ctrls->pd0_fast_coeff_est_level = 2;
        break;
    case 4:
        ctrls->update_skip_ctx_dc_sign_ctx = 0;
        ctrls->update_skip_coeff_ctx = 0;
        ctrls->coeff_rate_est_lvl = 2;
        ctrls->lpd0_qp_offset = 0;
        ctrls->pd0_fast_coeff_est_level = 2;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if CLN_INTRA_CTRLS
void set_intra_ctrls(PictureControlSet* pcs, ModeDecisionContext* ctx, uint8_t intra_level) {

    IntraCtrls* ctrls = &ctx->intra_ctrls;

    // If intra is disallowed at the pic level, must disallow at SB level
    if (pcs->skip_intra)
        intra_level = 0;

    assert(IMPLIES(pcs->slice_type == I_SLICE, intra_level > 0));

    switch (intra_level) {
    case 0:
        ctrls->enable_intra = 0;
        ctrls->intra_mode_end = DC_PRED;
        ctrls->angular_pred_level = 0;
        break;
    case 1:
        ctrls->enable_intra = 1;
        ctrls->intra_mode_end = PAETH_PRED;
        ctrls->angular_pred_level = 1;
        break;
    case 2:
        ctrls->enable_intra = 1;
        ctrls->intra_mode_end = PAETH_PRED;
        ctrls->angular_pred_level = 2;
        break;
    case 3:
        ctrls->enable_intra = 1;
        ctrls->intra_mode_end = SMOOTH_H_PRED;
        ctrls->angular_pred_level = 3;
        break;
    case 4:
        ctrls->enable_intra = 1;
        ctrls->intra_mode_end = SMOOTH_PRED;
        ctrls->angular_pred_level = 4;
        break;
    case 5:
        ctrls->enable_intra = 1;
        ctrls->intra_mode_end = DC_PRED;
        ctrls->angular_pred_level = 0;
        break;
    default:
        assert(0);
        break;
    }

    /* For PD1, the ability to skip intra must be set at the pic level to ensure all SBs
    perform inverse TX and generate the recon. */
    if (ctx->pd_pass == PD_PASS_1) {
        ctx->skip_intra = pcs->skip_intra;

        // Check user-defined settings
        if (pcs->parent_pcs_ptr->scs_ptr->static_config.enable_paeth == 0)
            ctrls->intra_mode_end = MIN(ctrls->intra_mode_end, SMOOTH_H_PRED);

        if (pcs->parent_pcs_ptr->scs_ptr->static_config.enable_smooth == 0)
            ctrls->intra_mode_end = MIN(ctrls->intra_mode_end, D67_PRED);

        if (pcs->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta == 0)
            ctrls->angular_pred_level = 0;
    }
    else {
        ctx->skip_intra = !(ctrls->enable_intra) || pcs->skip_intra;
    }
}
#endif
#if CLN_REG_PD1_TX_CTRLS
void set_tx_shortcut_ctrls(ModeDecisionContext* ctx, uint8_t tx_shortcut_level) {

    TxShortcutCtrls* ctrls = &ctx->tx_shortcut_ctrls;

    switch (tx_shortcut_level) {
    case 0:
        ctrls->bypass_tx_when_zcoeff = 0;
        ctrls->apply_pf_on_coeffs = 0;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->use_neighbour_info = 0;
        break;
    case 1:
        ctrls->bypass_tx_when_zcoeff = 1;
        ctrls->apply_pf_on_coeffs = 0;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->use_neighbour_info = 0;
        break;
    case 2:
        ctrls->bypass_tx_when_zcoeff = 2;
        ctrls->apply_pf_on_coeffs = 1;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->use_neighbour_info = 0;
        break;
    case 3:
        ctrls->bypass_tx_when_zcoeff = 2;
        ctrls->apply_pf_on_coeffs = 1;
        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info = 0;
        break;
    case 4:
        ctrls->bypass_tx_when_zcoeff = 2;
        ctrls->apply_pf_on_coeffs = 1;
        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info = 1;
        break;
    case 5:
        ctrls->bypass_tx_when_zcoeff = 2;
        ctrls->apply_pf_on_coeffs = 1;
        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info = 1;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if CLN_MDS0_CTRLS
void set_mds0_controls(PictureControlSet* pcs, ModeDecisionContext *ctx, uint8_t mds0_level) {

    Mds0Ctrls* ctrls = &ctx->mds0_ctrls;

    switch (mds0_level) {
    case 0:
        ctrls->mds0_dist_type = MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th = 0;
        break;
    case 1:
        ctrls->mds0_dist_type = MDS0_SSD;
        ctrls->enable_cost_based_early_exit = 0;
        ctrls->mds0_distortion_th = 0;
        break;
    case 2:
        ctrls->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 0;
        ctrls->mds0_distortion_th = 0;
        break;
    case 3:
        ctrls->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th = 50;
        break;
    case 4:
        ctrls->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th = 0;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#if FTR_VLPD0  && !CLN_MERGE_LPD0_VLPD0
// Set signals used for the very-light-pd0 path; only PD0 should call this function
// assumes NSQ OFF, no 4x4, no chroma, no mds3, no 64x64
EbErrorType signal_derivation_enc_dec_kernel_oq_very_light_pd0(
#if FTR_VLPD0_INTER_DEPTH
    PictureControlSet* pcs,
#endif
    ModeDecisionContext *ctx) {

    EbErrorType return_error = EB_ErrorNone;

    ctx->md_disallow_nsq = 1;
#if CLN_INTRA_CTRLS
    uint8_t intra_level = 0;
    if (pcs->slice_type == I_SLICE)
        intra_level = 1;
    else
        intra_level = 0;
    set_intra_ctrls(pcs, ctx, intra_level);
#else
    ctx->skip_intra = 1;
#endif
#if !CLN_INDEPTH
    set_in_depth_block_skip_ctrls(ctx, 1);
#endif
    ctx->shut_fast_rate = EB_TRUE;

#if FTR_VLPD0_INTER_DEPTH
    ctx->pd0_inter_depth_bias =  950 + pcs->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 1000)
        ctx->pd0_inter_depth_bias = ctx->pd0_inter_depth_bias - 150;
    else if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 500)
        ctx->pd0_inter_depth_bias = ctx->pd0_inter_depth_bias - 50;
    else if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 250)
        ctx->pd0_inter_depth_bias = ctx->pd0_inter_depth_bias + 50;
    else
        ctx->pd0_inter_depth_bias = ctx->pd0_inter_depth_bias + 150;
#else
    ctx->pd0_inter_depth_bias =  1000;
#endif
#if CLN_RATE_EST_CTRLS
    ctx->approx_inter_rate = 1;
#endif
    return return_error;
}
#endif

#if OPT_PD0_PATH
// Set signals used for light-pd0 path; only PD0 should call this function
// assumes NSQ OFF, no 4x4, no chroma, no TXT/TXS/RDOQ/SSSE, SB_64x64
#if CLN_LPD0_CTRL
#if CLN_MERGE_LPD0_VLPD0
EbErrorType signal_derivation_enc_dec_kernel_oq_light_pd0(
    SequenceControlSet* scs,
    PictureControlSet* pcs,
    ModeDecisionContext* ctx) {

    EbErrorType return_error = EB_ErrorNone;

    ctx->md_disallow_nsq = 1;
    ctx->inject_inter_candidates = 1;

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    ctx->shut_fast_rate = EB_TRUE;

    uint8_t intra_level = 0;
    if (pcs->slice_type == I_SLICE)
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
        ctx->pd0_inter_depth_bias = 950 + pcs->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
        if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 1000)
            ctx->pd0_inter_depth_bias = ctx->pd0_inter_depth_bias - 150;
        else if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 500)
            ctx->pd0_inter_depth_bias = ctx->pd0_inter_depth_bias - 50;
        else if (pcs->parent_pcs_ptr->me_8x8_cost_variance[ctx->sb_index] > 250)
            ctx->pd0_inter_depth_bias = ctx->pd0_inter_depth_bias + 50;
        else
            ctx->pd0_inter_depth_bias = ctx->pd0_inter_depth_bias + 150;
    }
    else {
        ctx->pd0_inter_depth_bias = 0;
    }
#if !CLN_CAND_REDUCTION_CTRLS
    if (ctx->pd0_level == VERY_LIGHT_PD0)
        ctx->reduce_unipred_candidates = 0;
    else
        ctx->reduce_unipred_candidates = ctx->pd0_level <= LIGHT_PD0_LVL1 ? 0 : 1;
#endif
    if (ctx->pd0_level == VERY_LIGHT_PD0)
        return return_error;

#if CLN_MDS0_CTRLS
    uint8_t mds0_level = 0;
    if (ctx->pd0_level <= LIGHT_PD0_LVL2)
        mds0_level = 4;
    else
        mds0_level = 0;

    set_mds0_controls(pcs, ctx, mds0_level);
#else
    if (ctx->pd0_level <= LIGHT_PD0_LVL2)
        ctx->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
    else
        ctx->mds0_dist_type = MDS0_SAD;
#endif
    set_chroma_controls(ctx, 0 /*chroma off*/);

    uint8_t pf_level = 1;
    if (pcs->slice_type != I_SLICE) {

        if (ctx->pd0_level <= LIGHT_PD0_LVL2) {
            pf_level = 1;
        }
        else {
            // Use ME distortion and variance detector to enable PF
            uint64_t use_pf_th = compute_pf_th(scs, pcs, ctx);
            uint32_t fast_lambda = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index]);

            if (ctx->pd0_level <= LIGHT_PD0_LVL3) {
                if (pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                    pf_level = 1;
                else
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
            }
            else {

                if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
                else {
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : (cost_64x64 < (2 * use_pf_th)) ? 2 : 1;

                    if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 16000 &&
                        ctx->depth_removal_ctrls.enabled &&
                        (ctx->depth_removal_ctrls.disallow_below_32x32 || ctx->depth_removal_ctrls.disallow_below_64x64)) {

                        if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 6 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index])
                            pf_level = 3;
                        else if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index])
                            pf_level = (cost_64x64 < (12 * use_pf_th)) ? 3 : 2;

                        pf_level = (cost_64x64 < (8 * use_pf_th)) ? 3 : MAX(2, pf_level);
                    }
                    else if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 16000 &&
                        ctx->depth_removal_ctrls.enabled &&
                        ctx->depth_removal_ctrls.disallow_below_16x16) {

                        if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 6 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index])
                            pf_level = (cost_64x64 < (4 * use_pf_th)) ? 3 : 2;

                        pf_level = (cost_64x64 < (2 * use_pf_th)) ? 3 : (cost_64x64 < (8 * use_pf_th)) ? MAX(2, pf_level) : pf_level;
                    }
                }
            }
        }
    }
    else {
        pf_level = 1;
    }
    set_pf_controls(ctx, pf_level);

    uint8_t subres_level;
    if (ctx->pd0_level <= LIGHT_PD0_LVL1) {
        subres_level = 0;
    }
    else {
        subres_level = 0;

        SbParams* sb_params_ptr = &pcs->parent_pcs_ptr->sb_params_array[ctx->sb_index];

        // The controls checks the deviation between: (1) the pred-to-src SAD of even rows and (2) the pred-to-src SAD of odd rows for each 64x64 to decide whether to use subres or not
        // then applies the result to the 64x64 block and to all children, therefore if incomplete 64x64 then shut subres
        if (sb_params_ptr->is_complete_sb) {
            // Use ME distortion and variance detector to enable subres
            uint64_t use_subres_th = compute_subres_th(scs, pcs, ctx);
            uint32_t fast_lambda = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index]);

            if (ctx->pd0_level <= LIGHT_PD0_LVL2) {
                if (pcs->slice_type == I_SLICE)
                    subres_level = 1;
                else
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
            }
            else if (ctx->pd0_level <= LIGHT_PD0_LVL3) {
                if (pcs->slice_type == I_SLICE)
                    subres_level = 1;
                else if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
                else
                    subres_level = 2;
            }
            else {
                if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    subres_level = (ctx->depth_removal_ctrls.enabled &&
                        (ctx->depth_removal_ctrls.disallow_below_16x16 ||
                            ctx->depth_removal_ctrls.disallow_below_32x32 ||
                            ctx->depth_removal_ctrls.disallow_below_64x64)) ? 2 : 1;
                else
                    subres_level = 2;
            }
#if !FIX_SUBRES_R2R
            ctx->subres_ctrls.odd_to_even_deviation_th = 5;
#endif
        }
        else {
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
#else
EbErrorType signal_derivation_enc_dec_kernel_oq_light_pd0(
    SequenceControlSet *scs,
    PictureControlSet *pcs,
    ModeDecisionContext *ctx) {

    EbErrorType return_error = EB_ErrorNone;

    set_chroma_controls(ctx, 0 /*chroma off*/);

    ctx->md_disallow_nsq = 1;
    ctx->inject_inter_candidates = 1;

#if CLN_INTRA_CTRLS
    uint8_t intra_level = 0;
    if (pcs->slice_type == I_SLICE)
        intra_level = 1;
    else if (ctx->pd0_level <= LIGHT_PD0_LVL1)
        intra_level = (pcs->temporal_layer_index == 0) ? 1 : 0;
    else
        intra_level = 0;
    set_intra_ctrls(pcs, ctx, intra_level);
#else
    if (pcs->slice_type == I_SLICE)
        ctx->skip_intra = 0;
     else if (ctx->pd0_level <= LIGHT_PD0_LVL1)
        ctx->skip_intra = (pcs->temporal_layer_index == 0) ? 0 : 1;
    else
        ctx->skip_intra = 1;
#endif

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    ctx->shut_fast_rate = EB_TRUE;

    if (ctx->pd0_level <= LIGHT_PD0_LVL2)
        ctx->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
    else
        ctx->mds0_dist_type = MDS0_SAD;
    uint8_t pf_level = 1;
    if (pcs->slice_type != I_SLICE) {

        if (ctx->pd0_level <= LIGHT_PD0_LVL2) {
            pf_level = 1;
        }
        else {
            // Use ME distortion and variance detector to enable PF
            uint64_t use_pf_th = compute_pf_th(scs, pcs, ctx);
            uint32_t fast_lambda = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index]);

            if (ctx->pd0_level <= LIGHT_PD0_LVL3) {
                if (pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                    pf_level = 1;
                else
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
            }
            else {

                if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
                else {
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : (cost_64x64 < (2 * use_pf_th)) ? 2 : 1;

                    if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 16000 &&
                        ctx->depth_removal_ctrls.enabled &&
                        (ctx->depth_removal_ctrls.disallow_below_32x32 || ctx->depth_removal_ctrls.disallow_below_64x64)) {

                        if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 6 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index])
                            pf_level = 3;
                        else if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index])
                            pf_level = (cost_64x64 < (12 * use_pf_th)) ? 3 : 2;

                        pf_level = (cost_64x64 < (8 * use_pf_th)) ? 3 : MAX(2, pf_level);
                    }
                    else if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 16000 &&
                        ctx->depth_removal_ctrls.enabled &&
                        ctx->depth_removal_ctrls.disallow_below_16x16) {

                        if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 6 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index])
                            pf_level = (cost_64x64 < (4 * use_pf_th)) ? 3 : 2;

                        pf_level = (cost_64x64 < (2 * use_pf_th)) ? 3 : (cost_64x64 < (8 * use_pf_th)) ? MAX(2, pf_level) : pf_level;
                    }
                }
            }
        }
    }
    else {
        pf_level = 1;
    }
    set_pf_controls(ctx, pf_level);

    uint8_t subres_level;
    if (ctx->pd0_level <= LIGHT_PD0_LVL1) {
        subres_level = 0;
    }
    else {
        subres_level = 0;

        SbParams *sb_params_ptr = &pcs->parent_pcs_ptr->sb_params_array[ctx->sb_index];

        // The controls checks the deviation between: (1) the pred-to-src SAD of even rows and (2) the pred-to-src SAD of odd rows for each 64x64 to decide whether to use subres or not
        // then applies the result to the 64x64 block and to all children, therefore if incomplete 64x64 then shut subres
        if (sb_params_ptr->is_complete_sb) {
            // Use ME distortion and variance detector to enable subres
            uint64_t use_subres_th = compute_subres_th(scs, pcs, ctx);
            uint32_t fast_lambda = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index]);

            if (ctx->pd0_level <= LIGHT_PD0_LVL2) {
                if (pcs->slice_type == I_SLICE)
                    subres_level = 1;
                else
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
            }
            else if (ctx->pd0_level <= LIGHT_PD0_LVL3) {
                if (pcs->slice_type == I_SLICE)
                    subres_level = 1;
                else if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
                else
                    subres_level = 2;
            }
            else {
                if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    subres_level = (ctx->depth_removal_ctrls.enabled &&
                    (ctx->depth_removal_ctrls.disallow_below_16x16 ||
                        ctx->depth_removal_ctrls.disallow_below_32x32 ||
                        ctx->depth_removal_ctrls.disallow_below_64x64)) ? 2 : 1;
                else
                    subres_level = 2;
            }
            ctx->subres_ctrls.odd_to_even_deviation_th = 5;
        }
        else {
            if (ctx->pd0_level <= LIGHT_PD0_LVL2)
                subres_level = 0;
            else
                subres_level = pcs->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 2;
        }
    }
    set_subres_controls(ctx, subres_level);

    ctx->pd0_inter_depth_bias = 0;

    ctx->reduce_unipred_candidates = ctx->pd0_level <= LIGHT_PD0_LVL1 ? 0 : 1;

    if (ctx->pd0_level <= LIGHT_PD0_LVL2)
        set_rate_est_ctrls(ctx, 2);
    else
        set_rate_est_ctrls(ctx, 0);

#if CLN_RATE_EST_CTRLS
    ctx->approx_inter_rate = 1;
#endif

    return return_error;
}
#endif
#else
EbErrorType signal_derivation_enc_dec_kernel_oq_light_pd0(
    SequenceControlSet *scs,
    PictureControlSet *pcs,
    ModeDecisionContext *ctx) {

    EbErrorType return_error = EB_ErrorNone;
    EbEncMode enc_mode = pcs->enc_mode;

    set_chroma_controls(ctx, 0 /*chroma off*/);

    ctx->md_disallow_nsq = 1;
    ctx->inject_inter_candidates = 1;

    if (pcs->slice_type == I_SLICE)
        ctx->skip_intra = 0;
    else if (enc_mode <= ENC_M1)
        ctx->skip_intra = 0;
#if TUNE_M10_M0
#if TUNE_M7_M8
    else if (enc_mode <= ENC_M8)
#else
    else if (enc_mode <= ENC_M7)
#endif
#else
    else if (enc_mode <= ENC_M6)
#endif
        ctx->skip_intra = (pcs->temporal_layer_index == 0) ? 0 : 1;
    else
        ctx->skip_intra = 1;
#if !CLN_INDEPTH
    uint8_t in_depth_block_skip_level = 0;
    if (pcs->parent_pcs_ptr->sc_class1)
        in_depth_block_skip_level = 0;
#if TUNE_M9_11_3
    else if (enc_mode <= ENC_M8)
        in_depth_block_skip_level = 0;
#endif
    else if (enc_mode <= ENC_M9)
#if TUNE_M9_11_3
        in_depth_block_skip_level = pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 0 : 1;
#else
        in_depth_block_skip_level = 0;
#endif
    else
        in_depth_block_skip_level = 1;

    set_in_depth_block_skip_ctrls(ctx, in_depth_block_skip_level);
#endif

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    ctx->shut_fast_rate = EB_TRUE;
#if !CLN_RATE_EST_CTRLS
    // Estimate the rate of the first(eob / fast_coeff_est_level) coeff(s), DC and last coeff only
    if (enc_mode <= ENC_M2)
        ctx->fast_coeff_est_level = 1;
#if TUNE_M8_M10 && !TUNE_M8_M11_MT
    else if (enc_mode <= ENC_M7)
#else
#if TUNE_M9_SLOWDOWN
    else if (enc_mode <= ENC_M9)
#else
    else if (enc_mode <= ENC_M8)
#endif
#endif
        ctx->fast_coeff_est_level = 2;
#if !TUNE_M10_M3_1
    else if (enc_mode <= ENC_M8)
        ctx->fast_coeff_est_level = (pcs->temporal_layer_index == 0) ? 2 : 4;
#endif
    else
        ctx->fast_coeff_est_level = (pcs->slice_type == I_SLICE) ? 2 : 4;
#endif
#if TUNE_MDS0_SAD
#if TUNE_M10_M0 && !TUNE_M7_11
    if (enc_mode <= ENC_M8)
#else
    if (enc_mode <= ENC_M9)
#endif
        ctx->mds0_dist_type = pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
    else
        ctx->mds0_dist_type = MDS0_SAD;
#else
#if SS_OPT_MDS0
    ctx->mds0_dist_type = (enc_mode <= ENC_MRS) ? MDS0_SAD : pcs->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
#else
    ctx->use_var_in_mds0 = (enc_mode <= ENC_MRS) ? 0 : pcs->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 0;
#endif
#endif
    uint8_t pf_level = 1;
    if (pcs->slice_type != I_SLICE) {
#if TUNE_M10_M3_1
#if TUNE_M9_SLOW && !TUNE_M9_M10
        if (enc_mode <= ENC_M9 || pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) {
#else
#if TUNE_M8_M11_MT
#if OPT_PD0_PF_LEVEL
        if (enc_mode <= ENC_M9) {
#else
        if (enc_mode <= ENC_M9 || pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) {
#endif
#else
        if (enc_mode <= ENC_M8 || pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) {
#endif
#endif
#else
        if (enc_mode <= ENC_M7 || pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) {
#endif
            pf_level = 1;
        }
        else {
            // Use ME distortion and variance detector to enable PF
            uint64_t use_pf_th = compute_pf_th(scs, pcs, ctx);
            uint32_t fast_lambda = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index]);
#if OPT_PD0_PF_LEVEL
#if TUNE_M11_SLOWDOWN
            if (enc_mode <= ENC_M11) {
#else
            if (enc_mode <= ENC_M10) {
#endif
                if (pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
                    pf_level = 1;
                else
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
            }
#if TUNE_NEW_M12
#if FTR_M13
            else if (enc_mode <= ENC_M13) {
#else
            else if (enc_mode <= ENC_M12) {
#endif
#else
            else if (enc_mode <= ENC_M11) {
#endif

                if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
                else {
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : (cost_64x64 < (2 * use_pf_th)) ? 2 : 1;

                    if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 16000 &&
                        ctx->depth_removal_ctrls.enabled &&
                        (ctx->depth_removal_ctrls.disallow_below_32x32 || ctx->depth_removal_ctrls.disallow_below_64x64)) {

                        if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 6 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index])
                            pf_level = 3;
                        else if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index])
                            pf_level = (cost_64x64 < (12 * use_pf_th)) ? 3 : 2;

                        pf_level = (cost_64x64 < (8 * use_pf_th)) ? 3 : MAX(2, pf_level);
                    }
                    else if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 16000 &&
                        ctx->depth_removal_ctrls.enabled &&
                        ctx->depth_removal_ctrls.disallow_below_16x16) {

                        if (pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index] < 6 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index] &&
                            pcs->parent_pcs_ptr->me_32x32_distortion[ctx->sb_index] < 4 * pcs->parent_pcs_ptr->me_16x16_distortion[ctx->sb_index])
                            pf_level = (cost_64x64 < (4 * use_pf_th)) ? 3 : 2;

                        pf_level = (cost_64x64 < (2 * use_pf_th)) ? 3 : (cost_64x64 < (8 * use_pf_th)) ? MAX(2, pf_level) : pf_level;
                    }
                }
            }
#else
            if (enc_mode <= ENC_M11) {
                pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
            }
#endif
            else {
                if (pcs->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                    pf_level = (cost_64x64 < ((use_pf_th * 3) >> 1)) ? 3 : 1;
                else
                    pf_level = (cost_64x64 < ((use_pf_th * 5) >> 1)) ? 3 : 1;
            }

        }
    }
    else {
        pf_level = 1;
    }
    set_pf_controls(ctx, pf_level);

    uint8_t subres_level;
    if (enc_mode <= ENC_M8) {
        subres_level = 0;
    }
    else {
        subres_level = 0;

        SbParams *sb_params_ptr = &pcs->parent_pcs_ptr->sb_params_array[ctx->sb_index];

        // The controls checks the deviation between: (1) the pred-to-src SAD of even rows and (2) the pred-to-src SAD of odd rows for each 64x64 to decide whether to use subres or not
        // then applies the result to the 64x64 block and to all children, therefore if incomplete 64x64 then shut subres
        if (sb_params_ptr->is_complete_sb) {
            // Use ME distortion and variance detector to enable subres
            uint64_t use_subres_th = compute_subres_th(scs, pcs, ctx);
            uint32_t fast_lambda = ctx->hbd_mode_decision ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs->parent_pcs_ptr->me_64x64_distortion[ctx->sb_index]);

#if TUNE_SUBRES_TX
#if CLN_RATE_EST_CTRLS
            if (enc_mode <= ENC_M9) {
#else
            if (enc_mode <= ENC_M10) {
#endif
                if (pcs->slice_type == I_SLICE)
                    subres_level = 1;
                else
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
            }
#if CLN_RATE_EST_CTRLS
            else if (enc_mode <= ENC_M10) {
                if (pcs->slice_type == I_SLICE)
                    subres_level = 1;
                else if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
                else
                    subres_level = 2;
            }
#endif
            else
#if TUNE_ADD_SUBRESS_FACTOR4
#if TUNE_SUBRES_LEVEL
#if CLN_RATE_EST_CTRLS
                if (pcs->parent_pcs_ptr->is_used_as_reference_flag)
                    subres_level = (ctx->depth_removal_ctrls.enabled &&
                    (ctx->depth_removal_ctrls.disallow_below_16x16 ||
                        ctx->depth_removal_ctrls.disallow_below_32x32 ||
                        ctx->depth_removal_ctrls.disallow_below_64x64)) ? 2 : 1;
                else
                    subres_level = 2;
#else
                subres_level = (ctx->depth_removal_ctrls.enabled &&
                (ctx->depth_removal_ctrls.disallow_below_16x16 ||
                    ctx->depth_removal_ctrls.disallow_below_32x32 ||
                    ctx->depth_removal_ctrls.disallow_below_64x64)) ? 2 : 1;
#endif
#else
                if (pcs->slice_type == I_SLICE)
                    subres_level = 1;
                else
                    subres_level = (ctx->depth_removal_ctrls.enabled &&
                    ((ctx->depth_removal_ctrls.disallow_below_16x16 && cost_64x64 < use_subres_th) ||
                        ctx->depth_removal_ctrls.disallow_below_32x32 ||
                        ctx->depth_removal_ctrls.disallow_below_64x64)) ? 2 : 1;
#endif
#else
                subres_level = 1;
#endif
#else
            if (pcs->slice_type == I_SLICE)
                subres_level = 1;
            else
                subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
#endif

            ctx->subres_ctrls.odd_to_even_deviation_th = 5;
        }
#if CLN_RATE_EST_CTRLS
        else {
            if (enc_mode <= ENC_M9)
                subres_level = 0;
            else
                subres_level = pcs->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 2;
        }
#endif
    }
    set_subres_controls(ctx, subres_level);
#if !CLN_REMOVE_UNUSED_FEATS
#if FTR_PD_EARLY_EXIT
    ctx->pd0_early_exit_th = enc_mode <= ENC_M8 ? 0 : 3;
#endif
#endif
#if FTR_PD_EARLY_EXIT
#if TUNE_M9_M11_OPTIMIZED_SUPER4KL
#if TUNE_M10_M11
    ctx->pd0_inter_depth_bias = enc_mode <= ENC_M10 ? 0 : 1003;
#else
    ctx->pd0_inter_depth_bias = enc_mode <= ENC_M9 ? 0 : 1003;
#endif
#else
    ctx->pd0_inter_depth_bias = enc_mode <= ENC_M8 ? 0 : 1003;
#endif
#endif
#if FTR_REDUCE_UNI_PRED
    ctx->reduce_unipred_candidates = enc_mode <= ENC_M8 ? 0 : 1;
#endif
#if CLN_RATE_EST_CTRLS
    if (enc_mode <= ENC_M9)
        set_rate_est_ctrls(ctx, 2);
    else
        set_rate_est_ctrls(ctx, 0);
#endif
    return return_error;
}
#endif
#endif

#if !CLN_REF_AREA
#if FTR_INTRA_DETECTOR
#if FTR_SELECTIVE_DLF
uint8_t ref_is_high_intra(PictureControlSet *pcs_ptr, uint8_t *dlf_th){
#else
uint8_t ref_is_high_intra(
    PictureControlSet *pcs_ptr) {
#endif
    if (pcs_ptr->slice_type == I_SLICE)
        return 100;
    uint8_t iperc = 0;

    EbReferenceObject *ref_obj_l0 =
            (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
    iperc = ref_obj_l0->intra_coded_area;
    if (pcs_ptr->slice_type == B_SLICE) {
        EbReferenceObject *ref_obj_l1 =
            (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
        iperc += ref_obj_l1->intra_coded_area;
    }
#if FTR_SELECTIVE_DLF
    if(dlf_th != NULL)
        *dlf_th = iperc;
#endif
    return (iperc > 70);
}
#endif
#if FTR_COEFF_DETECTOR
uint8_t ref_is_high_skip(PictureControlSet *pcs_ptr, uint8_t *skip_area) {
    uint8_t skip_perc = 0;

    if (pcs_ptr->temporal_layer_index == 0)
        return 0;

    EbReferenceObject *ref_obj_l0 =
        (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
    skip_perc = ref_obj_l0->skip_coded_area;
    if (pcs_ptr->slice_type == B_SLICE) {
        EbReferenceObject *ref_obj_l1 =
            (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
        skip_perc += ref_obj_l1->skip_coded_area;
    }
   * skip_area = skip_perc;

    return (skip_perc > 70);
}
#endif
#endif
#if CLN_LPD1_LVLS
void signal_derivation_enc_dec_kernel_oq_light_pd1(
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr) {

    EbEncMode lpd1_level = context_ptr->lpd1_ctrls.pd1_level;

    // Get ref info, used to set some feature levels
    const uint32_t picture_qp = pcs_ptr->picture_qp;
    uint32_t me_8x8_cost_variance = (uint32_t)~0;
    uint32_t me_64x64_distortion = (uint32_t)~0;
    uint8_t l0_was_skip = 0, l1_was_skip = 0;
    uint8_t l0_was_64x64_mvp = 0, l1_was_64x64_mvp = 0;

    if (pcs_ptr->slice_type != I_SLICE) {
        me_8x8_cost_variance = pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[context_ptr->sb_index];
        me_64x64_distortion = pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index];

        EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
        l0_was_skip = ref_obj_l0->sb_skip[context_ptr->sb_index], l1_was_skip = 1;
        l0_was_64x64_mvp = ref_obj_l0->sb_64x64_mvp[context_ptr->sb_index], l1_was_64x64_mvp = 1;
        if (pcs_ptr->slice_type == B_SLICE) {
            EbReferenceObject *ref_obj_l1 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
            l1_was_skip = ref_obj_l1->sb_skip[context_ptr->sb_index];
            l1_was_64x64_mvp = ref_obj_l1->sb_64x64_mvp[context_ptr->sb_index];
        }
    }
#if CLN_REF_AREA
    uint8_t ref_skip_perc = pcs_ptr->ref_skip_percentage;
#else
    uint8_t ref_intra_area = 0;
    ref_is_high_intra(pcs_ptr, &ref_intra_area);

    uint8_t ref_skip_area = 0;
    ref_is_high_skip(pcs_ptr, &ref_skip_area);
#endif

#if CLN_CAND_REDUCTION_CTRLS
    // Set candidate reduction levels
    uint8_t cand_reduction_level = 0;
    if (pcs_ptr->slice_type == I_SLICE)
        cand_reduction_level = 0;
    else if (lpd1_level <= LPD1_LVL_0)
        cand_reduction_level = 2;
    else if (lpd1_level <= LPD1_LVL_1)
        cand_reduction_level = 3;
    else if (lpd1_level <= LPD1_LVL_2)
        cand_reduction_level = 4;
    else
        cand_reduction_level = 5;

    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
        cand_reduction_level = 5;

    set_cand_reduction_ctrls(pcs_ptr, context_ptr, cand_reduction_level,
        picture_qp,
        me_8x8_cost_variance,
        me_64x64_distortion,
        l0_was_skip, l1_was_skip, ref_skip_perc);
#endif

#if !CLN_CAND_REDUCTION_CTRLS
    // Set feature levels

    uint8_t near_count_level = 0;
    if (lpd1_level <= LPD1_LVL_0)
        near_count_level = 2;
    else if (lpd1_level <= LPD1_LVL_2)
        near_count_level = 3;
    else
        near_count_level = 4;

#if FTR_OPT_MPASS_NEAR0
#if FTR_OP_TEST
    if (1)
#else
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
#endif
        near_count_level = 4;
#endif
    set_near_count_ctrls(context_ptr, near_count_level);

    // Can only use this feature when a single unipred ME candidate is selected, which is the case when the following conditions are true
    if (pcs_ptr->parent_pcs_ptr->ref_list0_count_try == 1 && pcs_ptr->parent_pcs_ptr->ref_list1_count_try == 1 && pcs_ptr->parent_pcs_ptr->use_best_me_unipred_cand_only) {
        if (lpd1_level <= LPD1_LVL_1)
            context_ptr->lpd1_mvp_best_me_list = 0;
        else
            context_ptr->lpd1_mvp_best_me_list = 1;
    }
    else {
        context_ptr->lpd1_mvp_best_me_list = 0;
    }
#endif
    if (lpd1_level <= LPD1_LVL_0)
        context_ptr->rdoq_level = 4;
    else
        context_ptr->rdoq_level = 5;
#if FTR_OPT_MPASS_RDOQ_OFF
#if FTR_OP_TEST
    if (1)
#else
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
#endif
        context_ptr->rdoq_level = 0;
#endif
    set_rdoq_controls(context_ptr, context_ptr->rdoq_level);
#if !CLN_INTRA_CTRLS
    // Set enable_paeth @ MD
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth == DEFAULT)
        if (lpd1_level <= LPD1_LVL_1)
            context_ptr->md_enable_paeth = 1;
        else
            context_ptr->md_enable_paeth = 0;
    else
        context_ptr->md_enable_paeth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth;

    // Set enable_smooth @ MD
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth == DEFAULT)
        if (lpd1_level <= LPD1_LVL_1)
            context_ptr->md_enable_smooth = 1;
        else
            context_ptr->md_enable_smooth = 0;
    else
        context_ptr->md_enable_smooth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth;
#endif
    if (lpd1_level <= LPD1_LVL_0)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 9 : 11) : 15;
    else if (lpd1_level <= LPD1_LVL_1)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13) : 15;
    else if (lpd1_level <= LPD1_LVL_2) {
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 14) : 15;
#if CLN_REF_AREA
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 50) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
#else
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 100) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
#endif
            me_8x8_cost_variance < (200 * picture_qp) &&
            me_64x64_distortion < (200 * picture_qp))
            context_ptr->md_subpel_me_level = 0;
    }
    else {
        if (pcs_ptr->temporal_layer_index != 0)
            context_ptr->md_subpel_me_level = 0;
        else {
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 14) : 15;
#if CLN_REF_AREA
            if (((l0_was_skip && l1_was_skip && ref_skip_perc > 50) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
#else
            if (((l0_was_skip && l1_was_skip && ref_skip_area > 100) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
#endif
                me_8x8_cost_variance < (200 * picture_qp) &&
                me_64x64_distortion < (200 * picture_qp))
                context_ptr->md_subpel_me_level = 0;
        }
    }

    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);
#if !CLN_INTRA_CTRLS
    // Set dc_cand_only_flag
    if (lpd1_level <= LPD1_LVL_1)
        context_ptr->dc_cand_only_flag = EB_FALSE;
    else
        context_ptr->dc_cand_only_flag = EB_TRUE;

    // Set intra_angle_delta @ MD
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta == DEFAULT)
        if (lpd1_level <= LPD1_LVL_1)
            context_ptr->md_intra_angle_delta = 1;
        else
            context_ptr->md_intra_angle_delta = 0;
    else
        context_ptr->md_intra_angle_delta = pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta;

    // Set disable_angle_z2_prediction_flag
    if (lpd1_level <= LPD1_LVL_1)
        context_ptr->disable_angle_z2_intra_flag = EB_FALSE;
    else
        context_ptr->disable_angle_z2_intra_flag = EB_TRUE;
#endif

#if CLN_MDS0_CTRLS
    uint8_t mds0_level = 0;
    if (lpd1_level <= LPD1_LVL_1)
        mds0_level = 4;
    else {
        mds0_level = 4;
#if CLN_REF_AREA
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 40) ||
#else
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 80) ||
#endif
            (me_8x8_cost_variance < (250 * picture_qp) &&
                me_64x64_distortion < (250 * picture_qp))))
            mds0_level = 0;
    }

    set_mds0_controls(pcs_ptr, context_ptr, mds0_level);
#else
    if (lpd1_level <= LPD1_LVL_1)
        context_ptr->mds0_dist_type = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
    else {
        context_ptr->mds0_dist_type = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;

        if (context_ptr->mds0_dist_type != MDS0_SAD &&
            ((l0_was_skip && l1_was_skip && ref_skip_area > 80) ||
            (me_8x8_cost_variance < (250 * picture_qp) &&
                me_64x64_distortion < (250 * picture_qp))))
            context_ptr->mds0_dist_type = MDS0_SAD;
    }
#endif
#if CLN_LPD1_TX_CTRLS
    uint8_t lpd1_tx_level = 0;
    if (lpd1_level <= LPD1_LVL_1)
        lpd1_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 3 : 2;
    else if (lpd1_level <= LPD1_LVL_3) {
        lpd1_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 4 : 2;
#if CLN_REF_AREA
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
#else
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 70) &&
#endif
            me_8x8_cost_variance < (800 * picture_qp) &&
            me_64x64_distortion < (800 * picture_qp)) ||
            (me_8x8_cost_variance < (100 * picture_qp) &&
                me_64x64_distortion < (100 * picture_qp)))
            lpd1_tx_level = 6;
    }
    else {
        lpd1_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 5 : 3;
#if CLN_REF_AREA
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
#else
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 70) &&
#endif
            me_8x8_cost_variance < (800 * picture_qp) &&
            me_64x64_distortion < (800 * picture_qp)) ||
            (me_8x8_cost_variance < (100 * picture_qp) &&
                me_64x64_distortion < (100 * picture_qp)))
            lpd1_tx_level = 6;
    }

    set_lpd1_tx_ctrls(context_ptr, lpd1_tx_level);
#else
    uint8_t mds1_skip_level = 0;
    if (lpd1_level <= LPD1_LVL_0)
        mds1_skip_level = 2;
    else if (lpd1_level <= LPD1_LVL_2)
        mds1_skip_level = 3;
    else if (lpd1_level <= LPD1_LVL_3)
        mds1_skip_level = 4;
    else
        mds1_skip_level = 5;

    set_mds1_skipping_controls(context_ptr, mds1_skip_level);

    uint8_t skip_tx_level = 0;
    if (lpd1_level <= LPD1_LVL_1)
        skip_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 3 : 2;
    else if (lpd1_level <= LPD1_LVL_3) {
        skip_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 4 : 2;
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 70) &&
            me_8x8_cost_variance < (800 * picture_qp) &&
            me_64x64_distortion < (800 * picture_qp)) ||
            (me_8x8_cost_variance < (100 * picture_qp) &&
                me_64x64_distortion < (100 * picture_qp)))
            skip_tx_level = 6;
    }
    else {
        skip_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 5 : 3;
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 70) &&
            me_8x8_cost_variance < (800 * picture_qp) &&
            me_64x64_distortion < (800 * picture_qp)) ||
            (me_8x8_cost_variance < (100 * picture_qp) &&
                me_64x64_distortion < (100 * picture_qp)))
            skip_tx_level = 6;
    }

    set_skip_tx_ctrls(context_ptr, skip_tx_level);
#endif
#if !CLN_CAND_REDUCTION_CTRLS
    if (lpd1_level <= LPD1_LVL_1)
    context_ptr->reduce_unipred_candidates = 1;
    else {
        context_ptr->reduce_unipred_candidates = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 3;
#if CLN_REF_AREA
        if ((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
#else
        if ((l0_was_skip && l1_was_skip && ref_skip_area > 70) &&
#endif
            me_8x8_cost_variance < (500 * picture_qp) &&
            me_64x64_distortion < (500 * picture_qp))
            context_ptr->reduce_unipred_candidates = 3;
    }
#endif
    uint8_t rate_est_level = 0;
    if (lpd1_level <= LPD1_LVL_0)
        rate_est_level = 4;
    else
        rate_est_level = pcs_ptr->slice_type == I_SLICE ? 4 : (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 0 : 4);

    set_rate_est_ctrls(context_ptr, rate_est_level);

    // If want to turn off approximating inter rate, must ensure that the approximation is also disabled
    // at the pic level (pcs_ptr->approx_inter_rate)
    context_ptr->approx_inter_rate = 1;

    uint8_t pf_level = 1;
    set_pf_controls(context_ptr, pf_level);
#if !CLN_LPD1_TX_CTRLS
    context_ptr->reduce_last_md_stage_candidate = 4;
#endif
#if !CLN_CAND_REDUCTION_CTRLS
    uint8_t neighbouring_mode_level = 1;
    set_use_neighbouring_mode_ctrls(context_ptr, neighbouring_mode_level);

    uint8_t redundant_cand_level = 1;
    set_redundant_cand_controls(context_ptr, redundant_cand_level);

    uint8_t eliminate_candidate_based_on_pme_me_results = 0;
    if (pcs_ptr->slice_type == I_SLICE) // not compatible with I_SLICE
        eliminate_candidate_based_on_pme_me_results = 0;
    else
        eliminate_candidate_based_on_pme_me_results = 2;

    set_cand_elimination_controls(context_ptr, eliminate_candidate_based_on_pme_me_results);
#endif
#if CLN_INTRA_CTRLS
    uint8_t intra_level = 0;
    if (lpd1_level <= LPD1_LVL_1)
        intra_level = 4;
    else
        intra_level = 5;

    set_intra_ctrls(pcs_ptr, context_ptr, intra_level);
#endif

    /* Set signals that have assumed values in the light-PD1 path (but need to be initialized as they may be checked) */

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    context_ptr->shut_fast_rate = EB_FALSE;
    context_ptr->uv_ctrls.enabled = 1;
    context_ptr->uv_ctrls.uv_mode = CHROMA_MODE_1;
    context_ptr->uv_ctrls.nd_uv_serach_mode = 0;
    set_cfl_ctrls(context_ptr, 0);
    context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;
    context_ptr->new_nearest_injection = 1;
    context_ptr->inject_inter_candidates = 1;
#if !CLN_INTRA_CTRLS
    context_ptr->skip_intra = pcs_ptr->skip_intra;
#endif
    context_ptr->blk_skip_decision = EB_TRUE;
    context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx = 0;
    context_ptr->rate_est_ctrls.update_skip_coeff_ctx = 0;
    context_ptr->subres_ctrls.odd_to_even_deviation_th = 0;
}
#else
#if FTR_VLPD1
void signal_derivation_enc_dec_kernel_oq_very_light_pd1(
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr) {

    EbEncMode enc_mode = pcs_ptr->enc_mode;

    // Get neighbour/colocated/pic-level info used in setting levels
#if !TUNE_VLPD1_LVLS
    const EbBool disallow_below_64x64 = (context_ptr->depth_removal_ctrls.enabled && context_ptr->depth_removal_ctrls.disallow_below_64x64);
#endif
    const uint32_t picture_qp = pcs_ptr->picture_qp;
#if FIX_LPD1_FOR_ISLICE
    uint32_t me_8x8_cost_variance = (uint32_t)~0;
    uint32_t me_64x64_distortion = (uint32_t)~0;
    uint8_t l0_was_skip = 0, l1_was_skip = 0;
    uint8_t l0_was_64x64_mvp = 0, l1_was_64x64_mvp = 0;

    if (pcs_ptr->slice_type != I_SLICE) {
        me_8x8_cost_variance = pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[context_ptr->sb_index];
        me_64x64_distortion = pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index];

        EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
        l0_was_skip = ref_obj_l0->sb_skip[context_ptr->sb_index], l1_was_skip = 1;
        l0_was_64x64_mvp = ref_obj_l0->sb_64x64_mvp[context_ptr->sb_index], l1_was_64x64_mvp = 1;
        if (pcs_ptr->slice_type == B_SLICE) {
            EbReferenceObject *ref_obj_l1 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
            l1_was_skip = ref_obj_l1->sb_skip[context_ptr->sb_index];
            l1_was_64x64_mvp = ref_obj_l1->sb_64x64_mvp[context_ptr->sb_index];
        }
    }
#else
    const uint32_t me_8x8_cost_variance = pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[context_ptr->sb_index];
    const uint32_t me_64x64_distortion = pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index];

    EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
    uint8_t l0_was_skip = ref_obj_l0->sb_skip[context_ptr->sb_index], l1_was_skip = 1;
    uint8_t l0_was_64x64_mvp = ref_obj_l0->sb_64x64_mvp[context_ptr->sb_index], l1_was_64x64_mvp = 1;
    if (pcs_ptr->slice_type == B_SLICE) {
        EbReferenceObject *ref_obj_l1 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
        l1_was_skip = ref_obj_l1->sb_skip[context_ptr->sb_index];
        l1_was_64x64_mvp = ref_obj_l1->sb_64x64_mvp[context_ptr->sb_index];
    }
#endif
    uint8_t ref_intra_area = 0;
    ref_is_high_intra(pcs_ptr, &ref_intra_area);

    uint8_t ref_skip_area = 0;
    ref_is_high_skip(pcs_ptr, &ref_skip_area);

    uint8_t near_count_level = 0;
    if (enc_mode <= ENC_M9)
        near_count_level = 1;
#if TUNE_M10_M11
    else if (enc_mode <= ENC_M11)
#else
    else if (enc_mode <= ENC_M10)
#endif
        near_count_level = 2;
    else
        near_count_level = 4;
#if FTR_OPT_MPASS_NEAR0
#if FTR_OP_TEST
    if (1)
#else
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
#endif
        near_count_level = 0;
#endif
    set_near_count_ctrls(context_ptr, near_count_level);

    // Can only use this feature when a single unipred ME candidate is selected, which is the case when the following conditions are true
    if (pcs_ptr->parent_pcs_ptr->ref_list0_count_try == 1 && pcs_ptr->parent_pcs_ptr->ref_list1_count_try == 1 && pcs_ptr->parent_pcs_ptr->use_best_me_unipred_cand_only) {
        if (enc_mode <= ENC_M10)
            context_ptr->lpd1_mvp_best_me_list = 0;
        else
            context_ptr->lpd1_mvp_best_me_list = 1;
    }
    else {
        context_ptr->lpd1_mvp_best_me_list = 0;
    }

    if (enc_mode <= ENC_M10)
        context_ptr->rdoq_level = 4;
    else
        context_ptr->rdoq_level = 5;
#if FTR_OPT_MPASS_RDOQ_OFF
#if FTR_OP_TEST
    if (1)
#else
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
#endif
        context_ptr->rdoq_level = 0;
#endif
    set_rdoq_controls(context_ptr, context_ptr->rdoq_level);

    uint8_t pf_level = 1;
    set_pf_controls(context_ptr, pf_level);

    if (enc_mode <= ENC_M8)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 7;
    else if (enc_mode <= ENC_M10)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 9 : 11) : 15;

#if OPT_M12_SPEED
    else if (enc_mode <= ENC_M11) {
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 14) : 15;
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 100) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
            me_8x8_cost_variance < (200 * picture_qp) &&
            me_64x64_distortion < (200 * picture_qp))
            context_ptr->md_subpel_me_level = 0;
    }
#if FTR_M13
    else if (enc_mode <= ENC_M12) {
#else
    else {
#endif
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 14) :
            pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->temporal_layer_index == 0 ? 12 : 15) : 15;
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 100) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
            me_8x8_cost_variance < (200 * picture_qp) &&
            me_64x64_distortion < (200 * picture_qp))
            context_ptr->md_subpel_me_level = 0;
    }
#else
    else {
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 14) : 15;
#if TUNE_VLPD1_LVLS
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 100) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
            me_8x8_cost_variance < (200 * picture_qp) &&
            me_64x64_distortion < (200 * picture_qp))
            context_ptr->md_subpel_me_level = 0;
#else
        if ((disallow_below_64x64 || (l0_was_skip && l1_was_skip && ref_skip_area > 80) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
            me_8x8_cost_variance < (300 * picture_qp) &&
            me_64x64_distortion < (300 * picture_qp))
            context_ptr->md_subpel_me_level = 0;
#endif
    }
#endif
#if FTR_M13
    else {
        if (pcs_ptr->temporal_layer_index != 0)
            context_ptr->md_subpel_me_level = 0;
        else {
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 14) : 15;

            if (((l0_was_skip && l1_was_skip && ref_skip_area > 100) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
                me_8x8_cost_variance < (200 * picture_qp) &&
                me_64x64_distortion < (200 * picture_qp))
                context_ptr->md_subpel_me_level = 0;
        }
    }
#endif

    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);
#if SS_FIX_MOVE_SKIP_INTRA_PIC
    context_ptr->skip_intra = pcs_ptr->skip_intra;
#else
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->skip_intra = 0;
    else {
        if (enc_mode <= ENC_M10) {
            context_ptr->skip_intra = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag || pcs_ptr->intra_percentage > 100 ? 0 : 1;
        }
        else
            context_ptr->skip_intra = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1) : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag || pcs_ptr->intra_percentage > 100 ? 0 : 1);
    }
#endif
#if !CLN_REMOVE_UNUSED_FEATS
    uint8_t early_cand_elimination = 0;
    if (pcs_ptr->slice_type == I_SLICE)
        early_cand_elimination = 0;
    else
        early_cand_elimination = 3;

    set_early_cand_elimination_controls(context_ptr, early_cand_elimination);
#endif
    /*reduce_last_md_stage_candidate
    0: OFF
    1: When MDS1 output has 0 coeffs, apply PFN4 and use DCT_DCT only
    2: 1 + disallow RDOQ and IFS when when MDS0 cand == MDS1 cand and
    the candidate does not belong to the best class
    3: 1 + 2 + remouve candidates when when MDS0 cand == MDS1 cand and they don't belong to the best class*/
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->reduce_last_md_stage_candidate = 0;
    else
        context_ptr->reduce_last_md_stage_candidate = 4;

    if (enc_mode <= ENC_M10)
        context_ptr->mds0_dist_type = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
    else {
        context_ptr->mds0_dist_type = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;

        if (context_ptr->mds0_dist_type != MDS0_SAD &&
            ((l0_was_skip && l1_was_skip && ref_skip_area > 80) ||
            (me_8x8_cost_variance < (250 * picture_qp) &&
             me_64x64_distortion  < (250 * picture_qp))))
            context_ptr->mds0_dist_type = MDS0_SAD;
    }

    uint8_t mds1_skip_level = 0;
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        mds1_skip_level = 0;
#if TUNE_M10_M11
    else if (enc_mode <= ENC_M11)
#else
    else if (enc_mode <= ENC_M10)
#endif
        mds1_skip_level = 2;
#if FTR_M13
    else if (enc_mode <= ENC_M12)
        mds1_skip_level = 5;
    else
        mds1_skip_level = 7;
#else
    else
        mds1_skip_level = 5;
#endif

    set_mds1_skipping_controls(context_ptr, mds1_skip_level);

    uint8_t skip_tx_level = 0;
    if (enc_mode <= ENC_M9)
        skip_tx_level = 1;
    else if (enc_mode <= ENC_M10)
        skip_tx_level = 3;
#if FTR_M13
    else if (enc_mode <= ENC_M12) {
#else
    else {
#endif
#if TUNE_M9_M11_OPTIMIZED_SUPER4KL
        skip_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 4 : 3;
#else
        skip_tx_level = 4;
#endif
#if TUNE_VLPD1_LVLS
        if (((l0_was_skip && l1_was_skip && ref_skip_area > 70) &&
            me_8x8_cost_variance < (800 * picture_qp) &&
            me_64x64_distortion < (800 * picture_qp)) ||
            (me_8x8_cost_variance < (100 * picture_qp) &&
             me_64x64_distortion < (100 * picture_qp)))
#if FTR_M13
            skip_tx_level = 7;
#else
            skip_tx_level = 5;
#endif
#else
        if ((disallow_below_64x64 || (l0_was_skip && l1_was_skip && ref_skip_area > 80)) &&
            me_8x8_cost_variance < (300 * picture_qp) &&
            me_64x64_distortion  < (300 * picture_qp))
            skip_tx_level = 5;
#endif
    }
#if FTR_M13
    else {
    skip_tx_level = 6;
    if (((l0_was_skip && l1_was_skip && ref_skip_area > 70) &&
        me_8x8_cost_variance < (800 * picture_qp) &&
        me_64x64_distortion < (800 * picture_qp)) ||
        (me_8x8_cost_variance < (100 * picture_qp) &&
            me_64x64_distortion < (100 * picture_qp)))
        skip_tx_level = 7;
    }
#endif

    set_skip_tx_ctrls(context_ptr, skip_tx_level);

    uint8_t redundant_cand_level = 0;
    if (enc_mode <= ENC_M9)
        redundant_cand_level = 0;
#if FTR_LESS_BI
    else if(enc_mode <= ENC_M11)
        redundant_cand_level = 2;
    else
        redundant_cand_level = 3;
#else
    else
        redundant_cand_level = 2;
#endif

    set_redundant_cand_controls(context_ptr, redundant_cand_level);

    if (enc_mode <= ENC_M10)
        context_ptr->reduce_unipred_candidates = 1;
    else {
        context_ptr->reduce_unipred_candidates = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 3;
#if TUNE_VLPD1_LVLS
        if ((l0_was_skip && l1_was_skip && ref_skip_area > 70) &&
            me_8x8_cost_variance < (500 * picture_qp) &&
            me_64x64_distortion < (500 * picture_qp))
            context_ptr->reduce_unipred_candidates = 3;
#else
        if ((disallow_below_64x64 || (l0_was_skip && l1_was_skip && ref_skip_area > 100)) &&
            me_8x8_cost_variance < (200 * picture_qp) &&
            me_64x64_distortion  < (200 * picture_qp))
            context_ptr->reduce_unipred_candidates = 3;
#endif
    }

    uint8_t use_neighbouring_mode = 1;
    set_use_neighbouring_mode_ctrls(context_ptr, use_neighbouring_mode);

#if CLN_RATE_EST_CTRLS
    uint8_t rate_est_level = 0;
    if (enc_mode <= ENC_M11)
        rate_est_level = 4;
    else
        rate_est_level = pcs_ptr->slice_type == I_SLICE ? 4 : (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 0 : 4);

    set_rate_est_ctrls(context_ptr, rate_est_level);

    context_ptr->approx_inter_rate = 1;
#endif
    /* Set signals that have assumed values in the light-PD1 path (but need to be initialized as they may be checked) */

    context_ptr->uv_ctrls.enabled = 1;
    context_ptr->uv_ctrls.uv_mode = CHROMA_MODE_1;
    context_ptr->uv_ctrls.nd_uv_serach_mode = 0;

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    context_ptr->shut_fast_rate = EB_FALSE;
#if SS_CLN_CFL_CTRLS
    set_cfl_ctrls(context_ptr, 0);
#else
    context_ptr->md_disable_cfl = 1;
#endif
    context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;

    context_ptr->new_nearest_injection = 1;
    context_ptr->inject_inter_candidates = 1;
    context_ptr->blk_skip_decision = EB_TRUE;
#if CLN_RATE_EST_CTRLS
    context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx = 0;
    context_ptr->rate_est_ctrls.update_skip_coeff_ctx = 0;
    context_ptr->subres_ctrls.odd_to_even_deviation_th = 0;
#else
    context_ptr->shut_skip_ctx_dc_sign_update = EB_TRUE;
    context_ptr->subres_ctrls.odd_to_even_deviation_th = 0;
    context_ptr->use_skip_coeff_context = 0;
#endif
    context_ptr->md_enable_paeth = 0;
    context_ptr->md_enable_smooth = 0;
#if OPT_LPD1_PME
    context_ptr->md_pme_level = 0;
    md_pme_search_controls(context_ptr, context_ptr->md_pme_level);
    context_ptr->md_subpel_pme_level = 0;
    md_subpel_pme_controls(context_ptr, context_ptr->md_subpel_pme_level);
#endif
    context_ptr->dc_cand_only_flag = EB_TRUE;
    context_ptr->md_intra_angle_delta = 0;
    context_ptr->disable_angle_z2_intra_flag = EB_TRUE;
    // no warp, no pme, already use DC only
    uint8_t eliminate_candidate_based_on_pme_me_results = 0;
    set_cand_elimination_controls(context_ptr, eliminate_candidate_based_on_pme_me_results);
}
#endif
#if LIGHT_PD1_MACRO
void signal_derivation_enc_dec_kernel_oq_light_pd1(
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr) {

EbEncMode enc_mode = pcs_ptr->enc_mode;

#if !OPT_REMOVE_TXT_LPD1
    uint8_t txt_level = 0;
    if (enc_mode <= ENC_M2)
        txt_level = 1;
    else if (enc_mode <= ENC_M4)
        txt_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 3;
    else if (enc_mode <= ENC_M9)
        txt_level = 5;
    else {
        txt_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? ((pcs_ptr->temporal_layer_index == 0) ? 6 : 7) : ((pcs_ptr->temporal_layer_index == 0) ? 6 : 8);
        if (pcs_ptr->intra_percentage < 170 && pcs_ptr->temporal_layer_index && pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_480p_RANGE && !pcs_ptr->parent_pcs_ptr->sc_class1) {
            txt_level = 0;
        }
    }

    set_txt_controls(context_ptr, txt_level, pcs_ptr->parent_pcs_ptr->input_resolution);
#endif

    uint8_t near_count_level = 0;
    if (enc_mode <= ENC_M9)
        near_count_level = 1;
#if TUNE_REDUCE_NR_LPD1
    else if (enc_mode <= ENC_M10)
        near_count_level = 2;
    else
        near_count_level = 3;
#else
    else
        near_count_level = 2;
#endif
#if FTR_OPT_MPASS_NEAR0
#if FTR_OP_TEST
    if (1)
#else
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
#endif
        near_count_level = 0;
#endif
    set_near_count_ctrls(context_ptr, near_count_level);

#if FTR_MVP_BEST_ME_LIST
    // Can only use this feature when a single unipred ME candidate is selected, which is the case when the following conditions are true
    if (pcs_ptr->parent_pcs_ptr->ref_list0_count_try == 1 && pcs_ptr->parent_pcs_ptr->ref_list1_count_try == 1 && pcs_ptr->parent_pcs_ptr->use_best_me_unipred_cand_only) {
#if TUNE_M10_M11
        if (enc_mode <= ENC_M11)
#else
        if (enc_mode <= ENC_M10)
#endif
            context_ptr->lpd1_mvp_best_me_list = 0;
        else
            context_ptr->lpd1_mvp_best_me_list = 1;
    }
    else {
        context_ptr->lpd1_mvp_best_me_list = 0;
    }
#endif

#if TUNE_M8_M10
    if (enc_mode <= ENC_M7)
#else
    if (enc_mode <= ENC_M9)
#endif
        context_ptr->rdoq_level = 1;
    else if (enc_mode <= ENC_M10)
        context_ptr->rdoq_level = 4;
    else
        context_ptr->rdoq_level = 5;
#if FTR_OPT_MPASS_RDOQ_OFF
#if FTR_OP_TEST
    if (1)
#else
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
#endif
        context_ptr->rdoq_level = 0;
#endif
    set_rdoq_controls(context_ptr, context_ptr->rdoq_level);

    // Set enable_paeth @ MD
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth == DEFAULT)
        context_ptr->md_enable_paeth = 1;
    else
        context_ptr->md_enable_paeth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth;

    // Set enable_smooth @ MD
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth == DEFAULT)
        context_ptr->md_enable_smooth = 1;
    else
        context_ptr->md_enable_smooth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth;

    uint8_t pf_level = 1;
    set_pf_controls(context_ptr, pf_level);

#if OPT_LPD1_PME
    // Set PME level
#if TUNE_PME_M0
    if (enc_mode <= ENC_M0)
        context_ptr->md_pme_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 2 : 1;
    else if (enc_mode <= ENC_M4)
        context_ptr->md_pme_level = 3;
#if TUNE_M8_M10_4K_SUPER
    else if (enc_mode <= ENC_M8)
#else
    else if (enc_mode <= ENC_M10)
#endif
        context_ptr->md_pme_level = 5;
    else
        context_ptr->md_pme_level = 0;
#else
    if (enc_mode <= ENC_M0)
        context_ptr->md_pme_level = 1;
    else if (enc_mode <= ENC_M4)
        context_ptr->md_pme_level = 2;
#if TUNE_M8_M10_4K_SUPER
#if TUNE_M7_M8_3
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M8)
#endif
#else
    else if (enc_mode <= ENC_M10)
#endif
        context_ptr->md_pme_level = 4;
    else
        context_ptr->md_pme_level = 0;
#endif
    md_pme_search_controls(context_ptr, context_ptr->md_pme_level);

    if (enc_mode <= ENC_M4)
        context_ptr->md_subpel_pme_level = 1;
#if TUNE_M8_M10_4K_SUPER
    else if (enc_mode <= ENC_M9)
        context_ptr->md_subpel_pme_level = 2;
    else if (enc_mode <= ENC_M10)
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 2 : 4;
#else
    else if (enc_mode <= ENC_M10)
        context_ptr->md_subpel_pme_level = 2;
#endif
#if FTR_M13
    else if (enc_mode <= ENC_M12)
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 3 : 4);
    else {
        if (pcs_ptr->temporal_layer_index != 0)
            context_ptr->md_subpel_pme_level = 0;
        else
            context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 3 : 4);
    }
#else
    else
#if TUNE_4K_M11
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 3 : 4);
#else
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : 3;
#endif
#endif

    md_subpel_pme_controls(context_ptr, context_ptr->md_subpel_pme_level);
#endif

    if (enc_mode <= ENC_M4)
        context_ptr->md_subpel_me_level = 1;
    else if (enc_mode <= ENC_M7)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 1 : 2;
    else if (enc_mode <= ENC_M8)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 7;
    else if (enc_mode <= ENC_M9)
#if TUNE_M9_11_3
#if FTR_VLPD1 // Fix levels to obey "onion ring"
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 7) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 10)) : 15;
#else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 7) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 10)) : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 15 : 16);
#endif
#else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 7) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 10);
#endif
    else if (enc_mode <= ENC_M10)
#if TUNE_M8_M10_4K_SUPER
#if FTR_VLPD1 // Fix levels to obey "onion ring"
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 9 : 11) : 15;
#else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 9 : 11) : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 15 : 16);
#endif
#else
        context_ptr->md_subpel_me_level = (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11;
#endif
#if OPT_M12_SPEED
    else if (enc_mode <= ENC_M11)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13) : 15;
#endif
#if OPT_SKIP_SUBPEL_ZZ
#if TUNE_NEW_M12
    else if (enc_mode <= ENC_M12)
#else
    else if (enc_mode <= ENC_M11)
#endif
#if TUNE_4K_M11
#if FTR_VLPD1 // Fix levels to obey "onion ring"
#if OPT_M12_SPEED
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13) :
        pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->temporal_layer_index == 0 ? 12 : 15) : 15;
#else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13) : 15;
#endif
#else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13) : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 15 : 16);
#endif
#else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13;
#endif
#if FTR_M13
    else {
        if (pcs_ptr->temporal_layer_index != 0)
            context_ptr->md_subpel_me_level = 0;
        else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13) : 15;
    }
#else
    else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 14;
#endif
#else
    else
#if TUNE_M11_SUBPEL
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 12 : 13;
#else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 12 : 13);
#endif
#endif
    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);

    // Set dc_cand_only_flag
#if FTR_M13
    context_ptr->dc_cand_only_flag = (enc_mode <= ENC_M12) ? EB_FALSE : (pcs_ptr->temporal_layer_index != 0) ? EB_TRUE : EB_FALSE;
#else
    context_ptr->dc_cand_only_flag = EB_FALSE;
#endif
    // Set intra_angle_delta @ MD
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta == DEFAULT)
        context_ptr->md_intra_angle_delta = 1;
    else
        context_ptr->md_intra_angle_delta = pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta;

    // Set disable_angle_z2_prediction_flag
    context_ptr->disable_angle_z2_intra_flag = EB_FALSE;

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    context_ptr->shut_fast_rate = EB_FALSE;

#if !OPT_REMOVE_TXT_LPD1
    // Estimate the rate of the first(eob / fast_coeff_est_level) coeff(s), DC and last coeff only
    if (enc_mode <= ENC_M10)
        context_ptr->fast_coeff_est_level = 1;
    else
        context_ptr->fast_coeff_est_level = 3;
#endif
#if SS_FIX_MOVE_SKIP_INTRA_PIC
    context_ptr->skip_intra = pcs_ptr->skip_intra;
#else
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->skip_intra = 0;
    else {
        if (enc_mode <= ENC_M8)
            context_ptr->skip_intra = 0;
        else if (enc_mode <= ENC_M10) {
            context_ptr->skip_intra = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag || pcs_ptr->intra_percentage > 100 ? 0 : 1;
        }
        else
#if TUNE_M9_11_3
            context_ptr->skip_intra = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1) : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag || pcs_ptr->intra_percentage > 100 ? 0 : 1);
#else
            context_ptr->skip_intra = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1;
#endif
    }
#endif
#if !CLN_REMOVE_UNUSED_FEATS
    uint8_t early_cand_elimination = 0;
    if (pcs_ptr->slice_type == I_SLICE)
        early_cand_elimination = 0;
    else if (enc_mode <= ENC_M6)
        early_cand_elimination = 0;
    else if (enc_mode <= ENC_M10)
        early_cand_elimination = 1;
    else
        early_cand_elimination = 3;

    set_early_cand_elimination_controls(context_ptr, early_cand_elimination);
#endif
    /*reduce_last_md_stage_candidate
    0: OFF
    1: When MDS1 output has 0 coeffs, apply PFN4 and use DCT_DCT only
    2: 1 + disallow RDOQ and IFS when when MDS0 cand == MDS1 cand and
    the candidate does not belong to the best class
    3: 1 + 2 + remouve candidates when when MDS0 cand == MDS1 cand and they don't belong to the best class*/
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->reduce_last_md_stage_candidate = 0;
    else
        if(enc_mode <= ENC_M6)
            context_ptr->reduce_last_md_stage_candidate = 0;
#if TUNE_M8_M10
        else if (enc_mode <= ENC_M7)
#else
        else if (enc_mode <= ENC_M8)
#endif
            context_ptr->reduce_last_md_stage_candidate = 3;
        else
            context_ptr->reduce_last_md_stage_candidate = 4;

    context_ptr->mds0_dist_type = (enc_mode <= ENC_MRS) ? MDS0_SAD : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;

    uint8_t eliminate_candidate_based_on_pme_me_results = 0;
    if (pcs_ptr->slice_type == I_SLICE)
        eliminate_candidate_based_on_pme_me_results = 0;
    else if (enc_mode <= ENC_M6)
        eliminate_candidate_based_on_pme_me_results = 0;
    else if (enc_mode <= ENC_M9)
        eliminate_candidate_based_on_pme_me_results = 1;
#if OPT_EARLY_ELIM_TH
#if TUNE_NEW_M12
#if FTR_M13
    else if (enc_mode <= ENC_M13)
#else
    else if (enc_mode <= ENC_M12)
#endif
#else
    else if (enc_mode <= ENC_M11)
#endif
        eliminate_candidate_based_on_pme_me_results = 2;
    else
        eliminate_candidate_based_on_pme_me_results = 3;
#else
    else
        eliminate_candidate_based_on_pme_me_results = 2;
#endif

    set_cand_elimination_controls(context_ptr, eliminate_candidate_based_on_pme_me_results);

    uint8_t mds1_skip_level = 0;
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        mds1_skip_level = 0;
#if TUNE_M9_M10
#if TUNE_M7_M8_3
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M8)
#endif
#else
    else if (enc_mode <= ENC_M9)
#endif
        mds1_skip_level = 0;
#if FTR_TX_NEIGH_INFO
    else if (enc_mode <= ENC_M10)
        mds1_skip_level = 2;
#if FTR_M13
    else if (enc_mode <= ENC_M12)
        mds1_skip_level = 3;
    else
        mds1_skip_level = 6;
#else
    else
        mds1_skip_level = 3;
#endif
#else
    else
#if CLN_DECPL_TX_FEATS
        mds1_skip_level = 2;
#else
        mds1_skip_level = 1;
#endif
#endif

    set_mds1_skipping_controls(context_ptr, mds1_skip_level);

#if FTR_SKIP_TX_LPD1
    uint8_t skip_tx_level = 0;
#if TUNE_M9_11_3
    if (enc_mode <= ENC_M8)
#else
    if (enc_mode <= ENC_M9)
#endif
        skip_tx_level = 0;
#if TUNE_M8_M10_4K_SUPER
    else if (enc_mode <= ENC_M9)
#else
    else if (enc_mode <= ENC_M10)
#endif
        skip_tx_level = 1;
#if FTR_M13
    else if (enc_mode <= ENC_M12)
#else
    else
#endif
#if OPT_TX_SKIP
        skip_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 3 : 2;
#else
        skip_tx_level = 2;
#endif
#if FTR_M13
    else
        skip_tx_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 5 : 2;
#endif

    set_skip_tx_ctrls(context_ptr, skip_tx_level);
#else
#if FTR_SKIP_COST_ZERO_COEFF
#if TUNE_M8_M11_MT
    if (enc_mode <= ENC_M9)
#else
    if (enc_mode <= ENC_M10)
#endif
        context_ptr->lpd1_zero_y_coeff_exit = 0;
    else
        context_ptr->lpd1_zero_y_coeff_exit = 1;
#endif
#endif
    uint8_t redundant_cand_level = 0;
    if (enc_mode <= ENC_M9)
        redundant_cand_level = 0;
    else
        redundant_cand_level = 1;

    set_redundant_cand_controls(context_ptr, redundant_cand_level);

#if TUNE_M9_M10 && !TUNE_M9_SLOWDOWN
    context_ptr->reduce_unipred_candidates = enc_mode <= ENC_M8 ? 0 : 1;
#else
    context_ptr->reduce_unipred_candidates = enc_mode <= ENC_M9 ? 0 : 1;
#endif
    uint8_t use_neighbouring_mode = 0;
#if TUNE_M8_M10_4K_SUPER
#if TUNE_M7_M8_3
    if (enc_mode <= ENC_M7)
#else
    if (enc_mode <= ENC_M8)
#endif
        use_neighbouring_mode = 0;
    else if (enc_mode <= ENC_M9)
        use_neighbouring_mode = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 0 : 1;
    else
        use_neighbouring_mode = 1;
#else
    if (enc_mode <= ENC_M9)
        use_neighbouring_mode = 0;
    else
        use_neighbouring_mode = 1;
#endif
    set_use_neighbouring_mode_ctrls(context_ptr, use_neighbouring_mode);

#if CLN_RATE_EST_CTRLS
    uint8_t rate_est_level = 0;
    if (enc_mode <= ENC_M11)
        rate_est_level = 4;
    else
        rate_est_level = pcs_ptr->slice_type == I_SLICE ? 4 : (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 0 : 4);

    set_rate_est_ctrls(context_ptr, rate_est_level);

    context_ptr->approx_inter_rate = 1;
#endif
    /* Set signals that have assumed values in the light-PD1 path (but need to be initialized as they may be checked) */

#if !OPT_UPDATE_MI_MAP
#if OPT_REMOVE_TXT_LPD1
    set_txt_controls(context_ptr, 0, pcs_ptr->parent_pcs_ptr->input_resolution);
#endif
    // mds1 not performed, so can't use coeff info
    context_ptr->bypass_tx_search_when_zcoef = 0;
#endif
    context_ptr->uv_ctrls.enabled = 1;
    context_ptr->uv_ctrls.uv_mode = CHROMA_MODE_1;
    context_ptr->uv_ctrls.nd_uv_serach_mode = 0;
#if SS_CLN_CFL_CTRLS
    set_cfl_ctrls(context_ptr, 0);
#else
    context_ptr->md_disable_cfl = 1;
#endif
    context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;
    context_ptr->new_nearest_injection = 1;
    context_ptr->inject_inter_candidates = 1;

#if !OPT_UPDATE_MI_MAP
    // compound off; sigs set in signal_derivation_block
    context_ptr->inter_compound_mode = 0;
#endif
    context_ptr->blk_skip_decision = EB_TRUE;
#if OPT_REDUCE_COPIES_LPD1
#if CLN_RATE_EST_CTRLS
    context_ptr->rate_est_ctrls.update_skip_ctx_dc_sign_ctx = 0;
    context_ptr->rate_est_ctrls.update_skip_coeff_ctx = 0;
#else
    context_ptr->shut_skip_ctx_dc_sign_update = EB_TRUE;
#endif
#else
    context_ptr->shut_skip_ctx_dc_sign_update = 0;
#endif
#if !OPT_UPDATE_MI_MAP
    set_spatial_sse_full_loop_level(context_ptr, 0);
    set_txs_controls(context_ptr, 0);

    // subres not used in MDS3 of PD1 (light-PD1 has no MDS1)
    set_subres_controls(context_ptr, 0);
#endif

    context_ptr->subres_ctrls.odd_to_even_deviation_th = 0;

#if !CLN_READ_REFINE_MVS
    context_ptr->md_sq_mv_search_level = 0;
    md_sq_motion_search_controls(context_ptr, context_ptr->md_sq_mv_search_level);

    context_ptr->md_nsq_mv_search_level = 0;
    md_nsq_motion_search_controls(context_ptr, context_ptr->md_nsq_mv_search_level);
#endif
#if !CLN_RATE_EST_CTRLS
#if FIX_SKIP_COEFF_CONTEXT
    context_ptr->use_skip_coeff_context = 0;
#endif
#endif
}
#endif
#endif
EbErrorType signal_derivation_enc_dec_kernel_oq(
    SequenceControlSet *sequence_control_set_ptr,
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    EbEncMode enc_mode = pcs_ptr->enc_mode;
    uint8_t pd_pass = context_ptr->pd_pass;
#if CLN_REG_PD_SIG_SET_0
    const uint32_t picture_qp = pcs_ptr->picture_qp;
    uint32_t me_8x8_cost_variance = (uint32_t)~0;
    uint32_t me_64x64_distortion = (uint32_t)~0;
    uint8_t l0_was_skip = 0, l1_was_skip = 0;
    uint8_t ref_skip_perc = pcs_ptr->ref_skip_percentage;

    set_cand_reduction_ctrls(
        pcs_ptr, context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->cand_reduction_level,
        picture_qp,
        me_8x8_cost_variance,
        me_64x64_distortion,
        l0_was_skip, l1_was_skip, ref_skip_perc);

    set_txt_controls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->txt_level,
        pcs_ptr->parent_pcs_ptr->input_resolution);

    set_tx_shortcut_ctrls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->tx_shortcut_level);

    set_interpolation_search_level_ctrls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->interpolation_search_level);

    set_chroma_controls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->chroma_level);

    set_cfl_ctrls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->cfl_level);
#else
#if CLN_CAND_REDUCTION_CTRLS
    const uint32_t picture_qp = pcs_ptr->picture_qp;
    uint32_t me_8x8_cost_variance = (uint32_t)~0;
    uint32_t me_64x64_distortion = (uint32_t)~0;
    uint8_t l0_was_skip = 0, l1_was_skip = 0;
    uint8_t ref_skip_perc = pcs_ptr->ref_skip_percentage;

    uint8_t cand_reduction_level = 0;
    if (pd_pass == PD_PASS_0)
        cand_reduction_level = 0;
    else if (pcs_ptr->slice_type == I_SLICE)
        cand_reduction_level = 0;
    else if (enc_mode <= ENC_M7)
        cand_reduction_level = 0;
    else if (enc_mode <= ENC_M9)
        cand_reduction_level = 1;
    else
        cand_reduction_level = 2;

    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
        cand_reduction_level = 5;

    set_cand_reduction_ctrls(pcs_ptr, context_ptr, cand_reduction_level,
        picture_qp,
        me_8x8_cost_variance,
        me_64x64_distortion,
        l0_was_skip, l1_was_skip, ref_skip_perc);
#endif
    uint8_t txt_level = 0;
    if (pd_pass == PD_PASS_0)
        txt_level = 0;
    else
        if (enc_mode <= ENC_M2)
            txt_level = 1;
#if TUNE_M0_M7_MEGA_FEB && !TUNE_MEGA_M9_M4
        else if (enc_mode <= ENC_M5)
#else
        else if (enc_mode <= ENC_M4)
#endif
            txt_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 3;
#if !TUNE_M0_M7_MEGA_FEB
        else if (enc_mode <= ENC_M5)
            txt_level = 3;
#endif
#if TUNE_NEW_TXT_LVLS
#if TUNE_M10_M0 && !TUNE_4K_M8_M11
        else if (enc_mode <= ENC_M9)
#else
        else if (enc_mode <= ENC_M8)
#endif
            txt_level = 5;
#if !TUNE_M1_M8
#if TUNE_4K_M8_M11
        else if (enc_mode <= ENC_M9) {
#if TUNE_REG_PD1
            txt_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 5 : ((pcs_ptr->temporal_layer_index == 0) ? 6 : 9);
#else
            txt_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 5 : ((pcs_ptr->temporal_layer_index == 0) ? 6 : 8);
#endif
#if CLN_REF_AREA
            if (pcs_ptr->ref_intra_percentage < 85 && pcs_ptr->temporal_layer_index && pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_1080p_RANGE && !pcs_ptr->parent_pcs_ptr->sc_class1) {
#else
            if (pcs_ptr->intra_percentage < 170 && pcs_ptr->temporal_layer_index && pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_1080p_RANGE && !pcs_ptr->parent_pcs_ptr->sc_class1) {
#endif
                txt_level = 0;
            }
        }
#endif
#endif
#if TUNE_NEW_M11
#if !TUNE_M10_M0
#if TUNE_M9_M10_MAR && !TUNE_M10_DEPTH_ME
        else if (enc_mode <= ENC_M9)
#else
#if TUNE_M11 && !TUNE_M10_M9_1
        else if (enc_mode <= ENC_M11)
#else
#if TUNE_M10_M3_1
        else if (enc_mode <= ENC_M9)
#else
        else if (enc_mode <= ENC_M10)
#endif
#endif
#endif
#if TUNE_TXT_LEVEL
            txt_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? ((pcs_ptr->temporal_layer_index == 0) ? 6 : 7) : ((pcs_ptr->temporal_layer_index == 0) ? 6 : 8);
#else
            txt_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 7;
#endif
#endif
#if TUNE_M10_M9_1
#if TUNE_REG_PD1
        else if (enc_mode <= ENC_M10)
#else
        else if (enc_mode <= ENC_M11)
#endif
#if TUNE_NEW_M11_2
        {
#if TUNE_REG_PD1
            txt_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? ((pcs_ptr->temporal_layer_index == 0) ? 6 : 8) : ((pcs_ptr->temporal_layer_index == 0) ? 6 : 9);
#else
            txt_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? ((pcs_ptr->temporal_layer_index == 0) ? 6 : 7) : ((pcs_ptr->temporal_layer_index == 0) ? 6 : 8);
#endif
#if TUNE_M11_2
#if CLN_REF_AREA
            if (pcs_ptr->ref_intra_percentage < 85 && pcs_ptr->temporal_layer_index && pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_480p_RANGE && !pcs_ptr->parent_pcs_ptr->sc_class1) {
#else
            if (pcs_ptr->intra_percentage < 170 && pcs_ptr->temporal_layer_index && pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_480p_RANGE && !pcs_ptr->parent_pcs_ptr->sc_class1) {
#endif
#else
            if (pcs_ptr->intra_percentage < 150 && pcs_ptr->temporal_layer_index && pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_480p_RANGE && !pcs_ptr->parent_pcs_ptr->sc_class1) {
#endif
                txt_level = 0;
            }
        }
#else
            txt_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? ((pcs_ptr->temporal_layer_index == 0) ? 6 : 7) : 0;
#endif
#endif
#if TUNE_REG_PD1
#if TUNE_NEW_M12
#if FTR_M13
        else if (enc_mode <= ENC_M13) {
#else
        else if (enc_mode <= ENC_M12) {
#endif
#else
        else if (enc_mode <= ENC_M11) {
#endif
            txt_level = (pcs_ptr->temporal_layer_index == 0) ? 7 : 10;
#if CLN_REF_AREA
            if (pcs_ptr->ref_intra_percentage < 85 && pcs_ptr->temporal_layer_index && pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_480p_RANGE && !pcs_ptr->parent_pcs_ptr->sc_class1) {
#else
            if (pcs_ptr->intra_percentage < 170 && pcs_ptr->temporal_layer_index && pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_480p_RANGE && !pcs_ptr->parent_pcs_ptr->sc_class1) {
#endif
                txt_level = 0;
            }
        }
#endif
        else
            txt_level = 0;
#else
        else
            txt_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 7;
#endif
        set_txt_controls(context_ptr, txt_level, pcs_ptr->parent_pcs_ptr->input_resolution);
#else
        else if (enc_mode <= ENC_M8)
            txt_level = 5;
        else
            txt_level = (pcs_ptr->temporal_layer_index == 0) ? 6 : 0;

    set_txt_controls(context_ptr, txt_level);
#endif

#if CLN_REG_PD1_TX_CTRLS
    uint8_t tx_shortcut_level = 0;
    if (pd_pass == PD_PASS_0)
        tx_shortcut_level = 0;
    else if (enc_mode <= ENC_M5)
        tx_shortcut_level = 0;
#if CLN_M6_M12_FEATURES
    else if (enc_mode <= ENC_M11)
        tx_shortcut_level = pcs_ptr->slice_type == I_SLICE ? 0 : 1;
#else
    else if (enc_mode <= ENC_M7)
        tx_shortcut_level = pcs_ptr->slice_type == I_SLICE ? 0 : 1;
    else if (enc_mode <= ENC_M9)
        tx_shortcut_level = pcs_ptr->slice_type == I_SLICE ? 0 : 2;
    else if (enc_mode <= ENC_M10)
        tx_shortcut_level = pcs_ptr->slice_type == I_SLICE ? 0 : 3;
#endif
    else if (enc_mode <= ENC_M12)
        tx_shortcut_level = pcs_ptr->slice_type == I_SLICE ? 0 : 4;
    else
        tx_shortcut_level = pcs_ptr->slice_type == I_SLICE ? 0 : 5;

    set_tx_shortcut_ctrls(context_ptr, tx_shortcut_level);
#else
    if (pd_pass == PD_PASS_0)
        context_ptr->bypass_tx_search_when_zcoef = 0;
    else
        if (pcs_ptr->slice_type == I_SLICE)
            context_ptr->bypass_tx_search_when_zcoef = 0;
        else
#if TUNE_M0_M7_MEGA_FEB
#if TUNE_MEGA_M9_M4
            context_ptr->bypass_tx_search_when_zcoef = (enc_mode <= ENC_M5) ? 0 : 1;
#else
            context_ptr->bypass_tx_search_when_zcoef = (enc_mode <= ENC_M4) ? 0 : 1;
#endif
#else
            context_ptr->bypass_tx_search_when_zcoef = (enc_mode <= ENC_M3) ? 0 : 1;
#endif
#endif
#if TUNE_BLOCK_SIZE
    uint8_t interpolation_search_level = 0;
    if (pd_pass == PD_PASS_0)
        interpolation_search_level = 0;
    else
        if (enc_mode <= ENC_MR)
            interpolation_search_level = 2;
#if FTR_COEFF_DETECTOR
#if TUNE_M10_M0
#if TUNE_4K_M8_M11
        else if (enc_mode <= ENC_M6)
#else
        else if (enc_mode <= ENC_M8)
#endif
#else
        else if (enc_mode <= ENC_M9)
#endif
            interpolation_search_level = 4;
#endif
#if TUNE_M8_M10_4K_SUPER
        else if (enc_mode <= ENC_M7 /* only for 4K*/) {
            interpolation_search_level = 4;
            if (pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_1080p_RANGE) {
#if CLN_REF_AREA
                const uint8_t th[INPUT_SIZE_COUNT] = { 100,85,85,55,50,45,40 };
                uint8_t skip_area = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : pcs_ptr->ref_skip_percentage;
#else
                uint8_t skip_area = 0;
                const uint8_t th[INPUT_SIZE_COUNT] = { 200,170,170,110,100,90,80 };
                ref_is_high_skip(pcs_ptr, &skip_area);
#endif
                if (skip_area > th[pcs_ptr->parent_pcs_ptr->input_resolution])
                    interpolation_search_level = 0;
            }
        }
#endif
#if TUNE_4K_M8_M11 && !TUNE_M7_M8_3
        else if (enc_mode <= ENC_M8) {
#if TUNE_M8_M10_4K_SUPER
            interpolation_search_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 4 : 6;
#else
            interpolation_search_level = 4;
#endif
            if (pcs_ptr->parent_pcs_ptr->input_resolution > INPUT_SIZE_1080p_RANGE) {
                uint8_t skip_area = 0;
#if TUNE_M8_M10_4K_SUPER
                const uint8_t th[INPUT_SIZE_COUNT] = { 200,170,170,110,100,30,30 };
#else
                const uint8_t th[INPUT_SIZE_COUNT] = { 200,170,170,110,100,90,80 };
#endif
                ref_is_high_skip(pcs_ptr, &skip_area);
                if (skip_area > th[pcs_ptr->parent_pcs_ptr->input_resolution])
                    interpolation_search_level = 0;
            }
        }
#endif

#if CLN_LIST0_ONLY_BASE_IFS_MFMV
        else {
            interpolation_search_level = 4;
            const uint8_t th[INPUT_SIZE_COUNT] = { 100,85,85,55,50,45,40 };
            uint8_t skip_area = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : pcs_ptr->ref_skip_percentage;
            if (skip_area > th[pcs_ptr->parent_pcs_ptr->input_resolution])
                interpolation_search_level = 0;
        }
    set_interpolation_search_level_ctrls(context_ptr, interpolation_search_level);
#else
#if TUNE_M7_M8_3
        else if (enc_mode <= ENC_M8) {
            interpolation_search_level = 4;
#if CLN_REF_AREA
            const uint8_t th[INPUT_SIZE_COUNT] = { 100,85,85,55,50,45,40 };
            uint8_t skip_area = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : pcs_ptr->ref_skip_percentage;
#else
            uint8_t skip_area = 0;
            const uint8_t th[INPUT_SIZE_COUNT] = { 200,170,170,110,100,90,80 };
            ref_is_high_skip(pcs_ptr, &skip_area);
#endif
            if (skip_area > th[pcs_ptr->parent_pcs_ptr->input_resolution])
                interpolation_search_level = 0;
        }
#endif
#if CLN_M10_M12_DIFFS
        else if (enc_mode <= ENC_M11) {
            interpolation_search_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 4 : 6;
#if CLN_REF_AREA
            const uint8_t th[INPUT_SIZE_COUNT] = { 100,85,85,55,50,15,15 };
            uint8_t skip_area = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : pcs_ptr->ref_skip_percentage;
#else
            uint8_t skip_area = 0;
            const uint8_t th[INPUT_SIZE_COUNT] = { 200,170,170,110,100,30,30 };
            ref_is_high_skip(pcs_ptr, &skip_area);
#endif
            if (skip_area > th[pcs_ptr->parent_pcs_ptr->input_resolution])
                interpolation_search_level = 0;
        }
#else
#if TUNE_M9_M10_MAR
#if TUNE_M7_M10_MT && !TUNE_M10_FASTER
        else if (enc_mode <= ENC_M10)
#else
#if TUNE_M10_M3_1 && !TUNE_4K_M8_M11
        else if (enc_mode <= ENC_M10)
#else
        else if (enc_mode <= ENC_M9)
#endif
#endif
#else
        else if (enc_mode <= ENC_M8)
#endif
#if FTR_COEFF_DETECTOR
        {
#if TUNE_M8_M10_4K_SUPER
            interpolation_search_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 4 : 6;
#if CLN_REF_AREA
            const uint8_t th[INPUT_SIZE_COUNT] = { 100,85,85,55,50,15,15 };
            uint8_t skip_area = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : pcs_ptr->ref_skip_percentage;
#else
            uint8_t skip_area = 0;
            const uint8_t th[INPUT_SIZE_COUNT] = {200,170,170,110,100,30,30};
#endif
#else
            interpolation_search_level = 4;
            uint8_t skip_area = 0;
            const uint8_t th[INPUT_SIZE_COUNT] = { 200,170,170,110,100,90,80 };
#endif
#if !CLN_REF_AREA
            ref_is_high_skip(pcs_ptr, &skip_area);
#endif
            if (skip_area > th[pcs_ptr->parent_pcs_ptr->input_resolution])
                interpolation_search_level = 0;
        }
#if TUNE_4K_M8_M11
        else if (enc_mode <= ENC_M10) {
            interpolation_search_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 4 : 6;
#if CLN_REF_AREA
            const uint8_t th[INPUT_SIZE_COUNT] = { 100,85,85,55,50,15,15 };
            uint8_t skip_area = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : pcs_ptr->ref_skip_percentage;
#else
            uint8_t skip_area = 0;
            const uint8_t th[INPUT_SIZE_COUNT] = { 200,170,170,110,100,30,30 };
            ref_is_high_skip(pcs_ptr, &skip_area);
#endif
            if (skip_area > th[pcs_ptr->parent_pcs_ptr->input_resolution])
                interpolation_search_level = 0;
        }
#endif
#else
            interpolation_search_level = 4;
#endif
#endif
#if TUNE_M10_DEPTH_ME
#if !TUNE_M7_M10_MT
#if TUNE_NEW_M10_M11 && !TUNE_M11
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M10)
#endif
            interpolation_search_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 5 : 0;
#endif
#if TUNE_M10_M9_1
#if TUNE_NEW_M11_2
#if TUNE_NEW_M12
#if FTR_M13
#if CLN_M10_M12_DIFFS
        else
#else
        else if (enc_mode <= ENC_M13)
#endif
#else
        else if (enc_mode <= ENC_M12)
#endif
#else
        else if (enc_mode <= ENC_M11)
#endif
#else
        else if (enc_mode <= ENC_M10)
#endif
#if OPT_IFS
#if FTR_COEFF_DETECTOR
        {
            interpolation_search_level = 6;
#if !CLN_REF_AREA
            uint8_t skip_area = 0;
#endif
#if TUNE_REG_PD1
            const uint8_t th[INPUT_SIZE_COUNT] = { 100,80,80,15,15,15,15 };
#if CLN_REF_AREA
            uint8_t skip_area = pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0 ? 0 : pcs_ptr->ref_skip_percentage;
#endif
#else
            const uint8_t th[INPUT_SIZE_COUNT] = { 200,170,170,30,30,30,30 };
#endif
#if !CLN_REF_AREA
            ref_is_high_skip(pcs_ptr, &skip_area);
#endif
            if (skip_area > th[pcs_ptr->parent_pcs_ptr->input_resolution])
                interpolation_search_level = 0;
        }
#else
            interpolation_search_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 6 : 7;
#endif
#else
            interpolation_search_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 4 : 5;
#endif
#endif
#if !CLN_M10_M12_DIFFS
        else
            interpolation_search_level = 0;
#endif
#else
        else
            interpolation_search_level = 5;
#endif
    set_interpolation_search_level_ctrls(context_ptr, interpolation_search_level);
#if !CLN_LIST0_ONLY_BASE_IFS_MFMV
#else
#endif
    if (pd_pass == PD_PASS_0)
        context_ptr->interpolation_search_level = IFS_OFF;
    else
        if (enc_mode <= ENC_MR)
            context_ptr->interpolation_search_level = IFS_MDS1;
#if TUNE_TXS_IFS_MFMV_DEPTH_M9
#if TUNE_M8_M9_FEB24
#if TUNE_MATCH_04_M8
        else if (enc_mode <= ENC_M9)
#else
        else if (enc_mode <= ENC_M8)
#endif
            context_ptr->interpolation_search_level = IFS_MDS3;
        else
            context_ptr->interpolation_search_level = (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? IFS_MDS3 : IFS_OFF;
#else
        else
            context_ptr->interpolation_search_level = IFS_MDS3;
#endif
#else
        else if (enc_mode <= ENC_M8)
            context_ptr->interpolation_search_level = IFS_MDS3;
        else
            context_ptr->interpolation_search_level = (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? IFS_MDS3 : IFS_OFF;
#endif
#endif
#endif

#if CHROMA_CLEANUP
    uint8_t chroma_level = 0;
    if (pd_pass == PD_PASS_0)
        chroma_level = 0;
    else if (sequence_control_set_ptr->static_config.set_chroma_mode == DEFAULT) {
        if (enc_mode <= ENC_MRS)
             chroma_level = 1;
        else if (enc_mode <= ENC_M1)
            chroma_level = 2;
        else if (enc_mode <= ENC_M5)
            chroma_level = 3;
        else
            chroma_level = 4;
    }
    else // use specified level
        chroma_level = sequence_control_set_ptr->static_config.set_chroma_mode;
    set_chroma_controls(context_ptr, chroma_level);
#else
    // Set Chroma Mode
    // Level                Settings
    // CHROMA_MODE_0  0     Full chroma search @ MD
    // CHROMA_MODE_1  1     Fast chroma search @ MD
    // CHROMA_MODE_2  2     Chroma blind @ MD + CFL @ EP
    if (pd_pass == PD_PASS_0)
        context_ptr->chroma_level = CHROMA_MODE_2;
    else if (sequence_control_set_ptr->static_config.set_chroma_mode == DEFAULT) {
        if (enc_mode <= ENC_M5)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else
            context_ptr->chroma_level = CHROMA_MODE_1;
    }
    else // use specified level
        context_ptr->chroma_level = sequence_control_set_ptr->static_config.set_chroma_mode;
#endif
    // Chroma independent modes search
    // Level                Settings
    // 0                    post first md_stage
    // 1                    post last md_stage
#if !CHROMA_CLEANUP
    if (enc_mode <= ENC_MRS) {
        context_ptr->chroma_at_last_md_stage = 0;
        context_ptr->chroma_at_last_md_stage_intra_th = (uint64_t)~0;
        context_ptr->chroma_at_last_md_stage_cfl_th = (uint64_t)~0;
    }
    else if (enc_mode <= ENC_M1) {
        context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
        context_ptr->chroma_at_last_md_stage_intra_th = 130;
        context_ptr->chroma_at_last_md_stage_cfl_th = 130;
    }
    else {
        context_ptr->chroma_at_last_md_stage = (context_ptr->chroma_level == CHROMA_MODE_0) ? 1 : 0;
        context_ptr->chroma_at_last_md_stage_intra_th = 100;
        context_ptr->chroma_at_last_md_stage_cfl_th = 100;
    }
#endif
#if SS_CLN_CFL_CTRLS
    uint8_t cfl_level = 0;
    if (pd_pass == PD_PASS_0)
        cfl_level = 0;
    else if (pcs_ptr->parent_pcs_ptr->sc_class1) {
        if (enc_mode <= ENC_M6)
            cfl_level = 1;
        else
            cfl_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 0;
    }
    else if (enc_mode <= ENC_M5)
        cfl_level = 1;
#if CLN_M6_M12_FEATURES
    else if (enc_mode <= ENC_M10)
        cfl_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 0;
#else
    else if (enc_mode <= ENC_M6)
        cfl_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 0;
    else if (enc_mode <= ENC_M10)
        cfl_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 0;
#endif
    else if (enc_mode <= ENC_M12)
        cfl_level = (pcs_ptr->slice_type == I_SLICE) ? 2 : 0;
    else
        cfl_level = 0;

    set_cfl_ctrls(context_ptr, cfl_level);
#else
    // Cfl level
    // Level                Settings
    // 0                    Allow cfl
    // 1                    Disable cfl
#if TUNE_M0_M7_MEGA_FEB
#if OPT_RECON_COPY
    if (pd_pass == PD_PASS_0)
        context_ptr->md_disable_cfl = EB_TRUE;
#if TUNE_SC_M7_M10
    else if (pcs_ptr->parent_pcs_ptr->sc_class1) {
        if (enc_mode <= ENC_M6)
            context_ptr->md_disable_cfl = EB_FALSE;
        else
            context_ptr->md_disable_cfl = (pcs_ptr->temporal_layer_index == 0) ? EB_FALSE : EB_TRUE;
    }
#endif
#if TUNE_M5_M6
    else if (enc_mode <= ENC_M5)
#else
    else if (enc_mode <= ENC_M6)
#endif
#else
    if (enc_mode <= ENC_M6)
#endif
#else
    if (enc_mode <= ENC_M5)
#endif
        context_ptr->md_disable_cfl = EB_FALSE;
#if TUNE_M9_MARCH && !TUNE_MEGA_M9_M4
    else if (enc_mode <= ENC_M8)
#else
#if (TUNE_M7_11 && !TUNE_M8_M10) || TUNE_M10_SLOWDOWN
    else if (enc_mode <= ENC_M10)
#else
    else if (enc_mode <= ENC_M9)
#endif
#endif
        context_ptr->md_disable_cfl = (pcs_ptr->temporal_layer_index == 0) ? EB_FALSE : EB_TRUE;
#if TUNE_M10_DEPTH_ME
#if TUNE_M11 && !TUNE_M10_M9_1
    else if (enc_mode <= ENC_M11)
#else
#if TUNE_NEW_M11_2
#if TUNE_NEW_M12
#if FTR_M13
    else if (enc_mode <= ENC_M13)
#else
    else if (enc_mode <= ENC_M12)
#endif
#else
    else if (enc_mode <= ENC_M11)
#endif
#else
    else if (enc_mode <= ENC_M10)
#endif
#endif
        context_ptr->md_disable_cfl = (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
#endif
#if TUNE_M10_M9_1 && !TUNE_NEW_M11_2
    else if (enc_mode <= ENC_M11)
        context_ptr->md_disable_cfl = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? ((pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE) : EB_TRUE;
#endif
    else
        context_ptr->md_disable_cfl = EB_TRUE;
#endif
#endif
#if !OPT_PD0_PATH
    // Set disallow_4x4
    context_ptr->disallow_4x4 = get_disallow_4x4 (enc_mode, pcs_ptr->slice_type);

#if !CLN_GEOM
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
#endif
    if (pd_pass == PD_PASS_0)
        context_ptr->md_disallow_nsq = enc_mode <= ENC_M0 ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
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
#if CLN_REG_PD_SIG_SET_1
    context_ptr->new_nearest_near_comb_injection = pd_pass == PD_PASS_0 ? 0 : pcs_ptr->new_nearest_near_comb_injection;
#else
    if (pd_pass == PD_PASS_0)
        context_ptr->new_nearest_near_comb_injection = 0;
    else
        if (sequence_control_set_ptr->static_config.new_nearest_comb_inject == DEFAULT)
            if (enc_mode <= ENC_M0)
                context_ptr->new_nearest_near_comb_injection = 1;
            else
                context_ptr->new_nearest_near_comb_injection = 0;
        else
            context_ptr->new_nearest_near_comb_injection = sequence_control_set_ptr->static_config.new_nearest_comb_inject;
#endif
#if !CLN_CAND_REDUCTION_CTRLS
#if CLN_REG_PD1_DIFFS
    uint8_t near_count_level = 0;
    if (pd_pass == PD_PASS_0)
        near_count_level = 0;
    else
        near_count_level = 1;
#else
    uint8_t near_count_level = 0;
    if (pd_pass == PD_PASS_0)
        near_count_level = 0;
    else
#if TUNE_MATCH_04_M8 && !TUNE_MEGA_M9_M4
        if (enc_mode <= ENC_M7)
#else
#if TUNE_M9_SLOW
#if TUNE_M9_M10
#if TUNE_M9_11_3
#if TUNE_REG_PD1
        if (enc_mode <= ENC_M10)
#else
        if (enc_mode <= ENC_M11)
#endif
#else
        if (enc_mode <= ENC_M10)
#endif
#else
        if (enc_mode <= ENC_M9)
#endif
#else
        if (enc_mode <= ENC_M8)
#endif
#endif
            near_count_level = 1;
#if !TUNE_M9_11_3
#if TUNE_M9_M10_MAR
#if TUNE_NEW_M11_2
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M10)
#endif
#else
        else if (enc_mode <= ENC_M9)
#endif
            near_count_level = 2;
#endif
        else
#if TUNE_NEW_M10_M11
#if TUNE_NEW_M11_2
#if TUNE_REG_PD1
            near_count_level = (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? 2 : 3;
#else
            near_count_level = (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
#endif
#else
            near_count_level = 4;
#endif
#else
            near_count_level = (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? 3 : 4;
#endif
#endif
#if FTR_OPT_MPASS_NEAR0
#if FTR_OP_TEST
    if (1)
#else
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
#endif
        near_count_level = 0;
#endif
    set_near_count_ctrls(context_ptr, near_count_level);
#endif


#if CLN_WM_SIG
    //set Warped-Motion controls from Picture level.
    set_wm_controls(context_ptr, pd_pass == PD_PASS_0 ? 0 : pcs_ptr->wm_level);
#else

#if FTR_NEW_WM_LVL
    // Set warped motion injection
    // Level                Settings
    // 0                    OFF
    // 1                    On
    uint8_t wm_level = 0;
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.allow_warped_motion) {
        if (pd_pass == PD_PASS_0)
            wm_level = 0;
#if !OPT_TXS_WM
#if TUNE_MEGA_M9_M4 && !TUNE_M7_M10_MT
        else if (enc_mode <= ENC_M9)
#else
#if TUNE_M7_11
        else if (enc_mode <= ENC_M7)
#else
        else if (enc_mode <= ENC_M8)
#endif
#endif
            wm_level = 1;
#endif
#if TUNE_M9_SLOW
#if TUNE_M10_M0 && !TUNE_M8_M11_MT
        else if (enc_mode <= ENC_M10)
#else
        else if (enc_mode <= ENC_M9)
#endif
#if OPT_TXS_WM
            wm_level = 1;
#else
            wm_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 1 : 0;
#endif
#endif
#if TUNE_NEW_M10_M11
#if TUNE_M10_FASTER && !TUNE_M10_M9_1
        else if (enc_mode <= ENC_M9)
#else
#if !TUNE_M10_M0
        else if (enc_mode <= ENC_M10)
#endif
#endif
#if TUNE_M10_M9_1
#if !TUNE_M10_M0
            wm_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 1 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 2 : 0);
#endif
#if OPT_TXS_WM
#if TUNE_M10_M11
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M10)
#endif
            wm_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE) ? 1 : 2;
        else
            wm_level = 2;
#else
        else
#if TUNE_M7_11
            wm_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE) ? 1 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 2 : 0);
#else
            wm_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 2 : 0;
#endif
#endif
#else
            wm_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 1 : 2;
        else
            wm_level = 2;
#endif
#else
        else
#if TUNE_M8_M9_FEB24
#if TUNE_MATCH_04_M8
            wm_level = 2;
#else
            wm_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 2 : 0;
#endif
#else
            wm_level = 2;
#endif
#endif
    }
    else {
        wm_level = 0;
    }

    set_wm_controls(context_ptr, wm_level);
#else
    // Set warped motion injection
    // Level                Settings
    // 0                    OFF
    // 1                    On
    if (pd_pass == PD_PASS_0) {
        context_ptr->warped_motion_injection = 0;
    }
    else
        context_ptr->warped_motion_injection = 1;
#endif

#endif

#if CLN_REG_PD_SIG_SET_1
    context_ptr->unipred3x3_injection = pd_pass == PD_PASS_0 ? 0 : pcs_ptr->unipred3x3_injection;
    context_ptr->bipred3x3_injection = pd_pass == PD_PASS_0 ? 0 : pcs_ptr->bipred3x3_injection;
    context_ptr->inject_inter_candidates = 1;
    context_ptr->inter_compound_mode = pd_pass == PD_PASS_0 ? 0 : pcs_ptr->inter_compound_mode;

    set_dist_based_ref_pruning_controls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->dist_based_ref_pruning);

    set_spatial_sse_full_loop_level(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->spatial_sse_full_loop_level);
#else
    // Set unipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
    if (pd_pass == PD_PASS_0) {
        context_ptr->unipred3x3_injection = 0;
    }
    else
        if (enc_mode <= ENC_M0)
            context_ptr->unipred3x3_injection = 1;
        else if (enc_mode <= ENC_M1)
            context_ptr->unipred3x3_injection = 2;
        else
            context_ptr->unipred3x3_injection = 0;

    // Set bipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
    if (pd_pass == PD_PASS_0) {
        context_ptr->bipred3x3_injection = 0;
    }
    else if (sequence_control_set_ptr->static_config.bipred_3x3_inject == DEFAULT) {
        if (enc_mode <= ENC_M1)
            context_ptr->bipred3x3_injection = 1;
        else if (enc_mode <= ENC_M5)
            context_ptr->bipred3x3_injection = 2;
        else
            context_ptr->bipred3x3_injection = 0;
    }
    else {
        context_ptr->bipred3x3_injection =
            sequence_control_set_ptr->static_config.bipred_3x3_inject;
    }

    context_ptr->inject_inter_candidates = 1;

    if (sequence_control_set_ptr->compound_mode) {
        if (sequence_control_set_ptr->static_config.compound_level == DEFAULT) {
            if (pd_pass == PD_PASS_0)
                context_ptr->inter_compound_mode = 0;
            else if (enc_mode <= ENC_MR)
                context_ptr->inter_compound_mode = 1;
#if TUNE_M0_M7_MEGA_FEB
#if TUNE_M3_M6_MEM_OPT
            else if (enc_mode <= ENC_M2)
#else
            else if (enc_mode <= ENC_M3)
#endif
#else
            else if (enc_mode <= ENC_M0)
#endif

                context_ptr->inter_compound_mode = 2;
#if !TUNE_M0_M7_MEGA_FEB
            else if (enc_mode <= ENC_M1)

                context_ptr->inter_compound_mode = (pcs_ptr->temporal_layer_index == 0) ? 2 : 3;
#endif
#if TUNE_M10_M3_1
            else if (enc_mode <= ENC_M3)
#else
            else if (enc_mode <= ENC_M4)
#endif

                context_ptr->inter_compound_mode = 3;
#if !TUNE_M10_M3_1
#if TUNE_MEGA_M9_M4
            else if (enc_mode <= ENC_M6)
#else
            else if (enc_mode <= ENC_M7)
#endif
                context_ptr->inter_compound_mode = (pcs_ptr->parent_pcs_ptr->input_resolution == INPUT_SIZE_240p_RANGE) ? 5 : 0;
#endif
            else
                context_ptr->inter_compound_mode = 0;

        }
        else {
            context_ptr->inter_compound_mode = sequence_control_set_ptr->static_config.compound_level;
        }
    }
    else {
        context_ptr->inter_compound_mode = 0;
    }

    // Set dist_based_ref_pruning
    if (pcs_ptr->parent_pcs_ptr->ref_list0_count_try > 1 || pcs_ptr->parent_pcs_ptr->ref_list1_count_try > 1) {
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->dist_based_ref_pruning = 0;
        else if (enc_mode <= ENC_MR)
            context_ptr->dist_based_ref_pruning = 1;
        else if (enc_mode <= ENC_M0)
            context_ptr->dist_based_ref_pruning = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
#if TUNE_M10_M0 && !TUNE_M8_M10_4K_SUPER
        else if (enc_mode <= ENC_M9)
#else
#if CLN_REF_PRUNING
        else
#else
        else if (enc_mode <= ENC_M8)
#endif
#endif
            context_ptr->dist_based_ref_pruning = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
#if !CLN_REF_PRUNING
#if TUNE_M10_M7
#if TUNE_M10_M0
        else if (enc_mode <= ENC_M10)
#else
        else if (enc_mode <= ENC_M9)
#endif
            context_ptr->dist_based_ref_pruning = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 4 : ((pcs_ptr->temporal_layer_index == 0) ? 2 : 4);
#endif
#endif
#if !CLN_REF_PRUNING
        else
            context_ptr->dist_based_ref_pruning = 4;
#endif
    }
    else {
        context_ptr->dist_based_ref_pruning = 0;
    }

    set_dist_based_ref_pruning_controls(context_ptr, context_ptr->dist_based_ref_pruning);
#if !CLN_MD_STAGING_CTRLS
    if (pd_pass == PD_PASS_0) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;
    }
    else
        if (enc_mode <= ENC_M4)
            context_ptr->md_staging_mode = MD_STAGING_MODE_2;
        else
            context_ptr->md_staging_mode = MD_STAGING_MODE_1;
#endif
    // spatial_sse_full_loop_level
#if TUNE_BLOCK_SIZE
    uint8_t spatial_sse_full_loop_level = 0;
    if (pd_pass == PD_PASS_0)
        spatial_sse_full_loop_level = 0;
    else if (sequence_control_set_ptr->static_config.spatial_sse_full_loop_level == DEFAULT)
#if TUNE_SC_M7_M10
        if (pcs_ptr->parent_pcs_ptr->sc_class1)
            spatial_sse_full_loop_level = 1;
#if CLN_FEAT_LEVEL
        else if (enc_mode <= ENC_M9)
#else
        else if (enc_mode <= ENC_M8)
#endif
#else
#if TUNE_M9_M10_MAR && !TUNE_M9_MARCH
        if (enc_mode <= ENC_M9)
#else
        if (enc_mode <= ENC_M8)
#endif
#endif
            spatial_sse_full_loop_level = 1;
#if TUNE_M10_DEPTH_ME
#if !CLN_FEAT_LEVEL
#if TUNE_M10_M7
        else if (enc_mode <= ENC_M9)
#else
        else if (enc_mode <= ENC_M10)
#endif
#if TUNE_M7_M10_MT
            spatial_sse_full_loop_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE) ? 2 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 0 : 1);
#else
            spatial_sse_full_loop_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE) ? 2 : 0;
#endif
#endif
        else
            spatial_sse_full_loop_level = 0;
#else
        else
            spatial_sse_full_loop_level = 2;
#endif
    else
        spatial_sse_full_loop_level =
        sequence_control_set_ptr->static_config.spatial_sse_full_loop_level;

    set_spatial_sse_full_loop_level(context_ptr, spatial_sse_full_loop_level);
#else
    if (pd_pass == PD_PASS_0)
        context_ptr->spatial_sse_full_loop_level = EB_FALSE;
    else if (sequence_control_set_ptr->static_config.spatial_sse_full_loop_level == DEFAULT)
#if TUNE_M8_M9_FEB24
#if TUNE_MATCH_04_M8
        if (enc_mode <= ENC_M9)
#else
        if (enc_mode <= ENC_M8)
#endif
            context_ptr->spatial_sse_full_loop_level = EB_TRUE;
        else if (enc_mode <= ENC_M10)
            context_ptr->spatial_sse_full_loop_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? EB_TRUE : EB_FALSE;
        else
            context_ptr->spatial_sse_full_loop_level = EB_FALSE;
#else
#if TUNE_M7_M10_PRESETS
        if (enc_mode <= ENC_M10)
#else
        if (enc_mode <= ENC_M9)
#endif
            context_ptr->spatial_sse_full_loop_level = EB_TRUE;
        else
            context_ptr->spatial_sse_full_loop_level = EB_FALSE;
#endif
    else
        context_ptr->spatial_sse_full_loop_level =
        sequence_control_set_ptr->static_config.spatial_sse_full_loop_level;
#endif
#endif
#if CHROMA_CLEANUP
    if (context_ptr->uv_ctrls.uv_mode <= CHROMA_MODE_1)
#else
    if (context_ptr->chroma_level <= CHROMA_MODE_1)
#endif
        context_ptr->blk_skip_decision = EB_TRUE;
    else
        context_ptr->blk_skip_decision = EB_FALSE;
    if (pd_pass == PD_PASS_0)
#if FIX_PD0_RDOQ
#if TUNE_M1_M8
        if (enc_mode <= ENC_M1)
#else
        if (enc_mode <= ENC_M0)
#endif
            context_ptr->rdoq_level = 1;
        else
            context_ptr->rdoq_level = 0;
#else
        context_ptr->rdoq_level = 0;
#endif
    else
#if TUNE_M10_M0 && !TUNE_M9_M10
        if (enc_mode <= ENC_M9)
#else
#if CLN_FEAT_LEVEL
        if (enc_mode <= ENC_M11)
#else
        if (enc_mode <= ENC_M8)
#endif
#endif
            context_ptr->rdoq_level = 1;
#if TUNE_NEW_M11
#if OPT_RDOQ_EOB
#if FTR_RDO_OPT
#if !TUNE_M10_M0
        else if (enc_mode <= ENC_M9)
            context_ptr->rdoq_level = 3;
#endif
#if !CLN_FEAT_LEVEL
#if CLN_M10_M12_DIFFS
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M10)
#endif
            context_ptr->rdoq_level = 4;
#endif
#if TUNE_NEW_M12
#if FTR_M13
#if CLN_M10_M12_DIFFS
        else
#else
        else if (enc_mode <= ENC_M13)
#endif
#else
        else if (enc_mode <= ENC_M12)
#endif
#else
        else if (enc_mode <= ENC_M11)
#endif
            context_ptr->rdoq_level = 5;
#else
        else if (enc_mode <= ENC_M10)
            context_ptr->rdoq_level = 3;
        else if (enc_mode <= ENC_M11)
            context_ptr->rdoq_level = 4;
#endif
#else
        else if (enc_mode <= ENC_M10)
            context_ptr->rdoq_level = 3;
#endif
#if !CLN_M10_M12_DIFFS
        else
            context_ptr->rdoq_level = 0;
#endif
#else
        else
#if FTR_RDOQ_ONLY_DCT_DCT
            context_ptr->rdoq_level = 3;
#else
            context_ptr->rdoq_level = 2;
#endif
#endif
#if FTR_OPT_MPASS_RDOQ_OFF
#if FTR_OP_TEST
    if (1)
#else
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->rc_stat_gen_pass_mode)
#endif
        context_ptr->rdoq_level = 0;
#endif
    set_rdoq_controls(context_ptr, context_ptr->rdoq_level);

    // Derive redundant block
#if SS_OPT_MD
    if (pd_pass == PD_PASS_0 || context_ptr->md_disallow_nsq)
#else
    if (pd_pass == PD_PASS_0)
#endif
        context_ptr->redundant_blk = EB_FALSE;
    else
        if (sequence_control_set_ptr->static_config.enable_redundant_blk ==
            DEFAULT)
#if TUNE_NEW_M11
            context_ptr->redundant_blk = EB_TRUE;
#else
            if (enc_mode <= ENC_M10)
                context_ptr->redundant_blk = EB_TRUE;
            else
                context_ptr->redundant_blk = EB_FALSE;
#endif
        else
            context_ptr->redundant_blk =
            sequence_control_set_ptr->static_config.enable_redundant_blk;
#if CLN_REG_PD_SIG_SET_2
    set_parent_sq_coeff_area_based_cycles_reduction_ctrls(
        context_ptr,
        pcs_ptr->parent_pcs_ptr->input_resolution,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->parent_sq_coeff_area_based_cycles_reduction_level);
#else
    uint8_t parent_sq_coeff_area_based_cycles_reduction_level = 0;
    if (pd_pass == PD_PASS_0)
        parent_sq_coeff_area_based_cycles_reduction_level = 0;
    else if (enc_mode <= ENC_MRS)
        parent_sq_coeff_area_based_cycles_reduction_level = 0;
    else if (enc_mode <= ENC_MR)
        parent_sq_coeff_area_based_cycles_reduction_level = pcs_ptr->slice_type == I_SLICE ? 0 : 1;
    else if (enc_mode <= ENC_M1)
        parent_sq_coeff_area_based_cycles_reduction_level = pcs_ptr->slice_type == I_SLICE ? 0 : 2;
    else if (enc_mode <= ENC_M2)

        parent_sq_coeff_area_based_cycles_reduction_level = pcs_ptr->slice_type == I_SLICE ? 0
        : (pcs_ptr->temporal_layer_index == 0) ? 2 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 4 : 7;


    else
        parent_sq_coeff_area_based_cycles_reduction_level = pcs_ptr->slice_type == I_SLICE ? 5 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 6 : 7;
    set_parent_sq_coeff_area_based_cycles_reduction_ctrls(context_ptr, pcs_ptr->parent_pcs_ptr->input_resolution, parent_sq_coeff_area_based_cycles_reduction_level);

#endif
#if CLN_REG_PD_SIG_SET_2
    context_ptr->sq_weight = pd_pass == PD_PASS_0 ? (uint32_t) ~0 : pcs_ptr->sq_weight;

    context_ptr->max_part0_to_part1_dev = pd_pass == PD_PASS_0 ? 0 : pcs_ptr->max_part0_to_part1_dev;
#else
    // Weighting (expressed as a percentage) applied to
    // square shape costs for determining if a and b
    // shapes should be skipped. Namely:
    // skip HA, HB, and H4 if h_cost > (weighted sq_cost)
    // skip VA, VB, and V4 if v_cost > (weighted sq_cost)
    if (pd_pass == PD_PASS_0)
        context_ptr->sq_weight = (uint32_t)~0;
    else
        if (enc_mode <= ENC_MRS)
            context_ptr->sq_weight = (uint32_t)~0;
        else
#if TUNE_M0_M7_MEGA_FEB
            if (enc_mode <= ENC_M0)
#else
            if (enc_mode <= ENC_MR)
#endif
                context_ptr->sq_weight = 105;
#if !TUNE_M0_M7_MEGA_FEB
            else
                if (enc_mode <= ENC_M0)
                    context_ptr->sq_weight = 102;
#endif
#if TUNE_M10_M3_1
                else
                    context_ptr->sq_weight = 95;
#else
#if TUNE_M0_M7_MEGA_FEB
#if TUNE_M3_M6_MEM_OPT
                else if (enc_mode <= ENC_M3)
#else
                else if (enc_mode <= ENC_M2)
#endif
#else
                else if (enc_mode <= ENC_M1)
#endif
                    context_ptr->sq_weight = 95;
                else
                    context_ptr->sq_weight = 90;
#endif


    // max_part0_to_part1_dev is used to:
    // (1) skip the H_Path if the deviation between the Parent-SQ src-to-recon distortion of (1st quadrant + 2nd quadrant) and the Parent-SQ src-to-recon distortion of (3rd quadrant + 4th quadrant) is less than TH,
    // (2) skip the V_Path if the deviation between the Parent-SQ src-to-recon distortion of (1st quadrant + 3rd quadrant) and the Parent-SQ src-to-recon distortion of (2nd quadrant + 4th quadrant) is less than TH.
    if (pd_pass == PD_PASS_0)
        context_ptr->max_part0_to_part1_dev = 0;
    else
#if TUNE_M0_M7_MEGA_FEB
        if (enc_mode <= ENC_M3)
#else
        if (enc_mode <= ENC_M2)
#endif
            context_ptr->max_part0_to_part1_dev = 0;
#if !TUNE_M0_M7_MEGA_FEB
        else if (enc_mode <= ENC_M3)
            context_ptr->max_part0_to_part1_dev = (pcs_ptr->temporal_layer_index == 0) ? 0 : 50;
#endif
        else
            context_ptr->max_part0_to_part1_dev = 100;
#endif



    // Set pic_obmc_level @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_pic_obmc_level = 0;
    else
        context_ptr->md_pic_obmc_level =
        pcs_ptr->parent_pcs_ptr->pic_obmc_level;

    set_obmc_controls(context_ptr, context_ptr->md_pic_obmc_level);

#if CLN_REG_PD_SIG_SET_2
    context_ptr->md_inter_intra_level = pcs_ptr->md_inter_intra_level;
    set_inter_intra_ctrls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->md_inter_intra_level);
#else
    // Set enable_inter_intra @ MD
    //Block level switch, has to follow the picture level
    // inter intra pred                      Settings
    // 0                                     OFF
    // 1                                     FULL
    // 2                                     FAST 1 : Do not inject for unipred3x3 or PME inter candidates
    // 3                                     FAST 2 : Level 1 + do not inject for non-closest ref frames or ref frames with high distortion
    if (pcs_ptr->parent_pcs_ptr->slice_type != I_SLICE && sequence_control_set_ptr->seq_header.enable_interintra_compound) {
        if (pd_pass == PD_PASS_0)
            context_ptr->md_inter_intra_level = 0;
        else if (enc_mode <= ENC_M2)
            context_ptr->md_inter_intra_level = 1;
        else
            context_ptr->md_inter_intra_level = 0;
    }
    else
        context_ptr->md_inter_intra_level = 0;

    set_inter_intra_ctrls(context_ptr, context_ptr->md_inter_intra_level);
#endif
#if CLN_REG_PD_SIG_SET_2
    set_txs_controls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->txs_level);

    set_nic_controls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 17 : pcs_ptr->nic_level);
#else

#if !CLN_INTRA_CTRLS
    // Set enable_paeth @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_enable_paeth = 1;
    else if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth == DEFAULT)
        context_ptr->md_enable_paeth = 1;
    else
        context_ptr->md_enable_paeth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_paeth;

    // Set enable_smooth @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_enable_smooth = 1;
    else if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth == DEFAULT)
        context_ptr->md_enable_smooth = 1;
    else
        context_ptr->md_enable_smooth = (uint8_t)pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.enable_smooth;
#endif
#if OPT_TXS_SEARCH
    // Set md_tx_size_search_mode @ MD
    uint8_t txs_level = 0;
    if (pcs_ptr->parent_pcs_ptr->tx_size_search_mode == 0)
        txs_level = 0;
    else if (pd_pass == PD_PASS_0)
        txs_level = 0;
    else if (enc_mode <= ENC_MRS)
        txs_level = 1;
    else if (enc_mode <= ENC_MR)
        txs_level = 2;
#if TUNE_M1_M8
    else if (enc_mode <= ENC_M2)
#else
    else if (enc_mode <= ENC_M1)
#endif
        txs_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 3;
#if TUNE_M9_M10_MAR
#if TUNE_M7_11
#if TUNE_M11_SLOWDOWN
    else if (enc_mode <= ENC_M11)
#else
    else if (enc_mode <= ENC_M10)
#endif
#else
    else if (enc_mode <= ENC_M9)
#endif
#else
    else if (enc_mode <= ENC_M8)
#endif
        txs_level = 3;
#if TUNE_M10_M3_1 && !TUNE_M7_11
    else if (enc_mode <= ENC_M10)
        txs_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 4 : 3;
#endif
#if TUNE_NEW_M11
#if TUNE_NEW_M10_M11
#if !TUNE_M10_M7
    else if (enc_mode <= ENC_M10)
        txs_level = 4;
#endif
#if TUNE_NEW_M12
#if FTR_M13
    else if (enc_mode <= ENC_M13)
#else
    else if (enc_mode <= ENC_M12)
#endif
#else
    else if (enc_mode <= ENC_M11)
#endif
#else
    else if (enc_mode <= ENC_M10)
#endif
#if TUNE_M8_M9_FEB24
#if TUNE_TXS_M11
#if OPT_TXS_WM
        txs_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_240p_RANGE) ? 5: (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 4 : 1;
#else
        txs_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 4 : 1;
#endif
#else
        txs_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 4 : 0;
#endif
#else
        txs_level = 4;
#endif
    else
        txs_level = 0;
#else
    else
        txs_level = 4;
#endif
    set_txs_controls(context_ptr, txs_level);
#else
    // Set md_tx_size_search_mode @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_tx_size_search_mode = 0;
    else
        context_ptr->md_tx_size_search_mode = pcs_ptr->parent_pcs_ptr->tx_size_search_mode;

    // Assign whether to use TXS in inter classes (if TXS is ON)
    // 0 OFF - Use TXS for intra candidates only
    // 1 ON  - Use TXS for all candidates
    // 2 ON  - INTER TXS restricted to max 1 depth
    if (enc_mode <= ENC_MRS)
        context_ptr->md_staging_tx_size_level = 1;

    else if (enc_mode <= ENC_MR)
        context_ptr->md_staging_tx_size_level = 2;

    else if (enc_mode <= ENC_M1)
        context_ptr->md_staging_tx_size_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 0;

    else
        context_ptr->md_staging_tx_size_level = 0;
#endif
#if CLN_NIC_SIGS
    uint8_t nic_level = get_nic_level(pd_pass, enc_mode, pcs_ptr->temporal_layer_index);
    set_nic_controls(context_ptr, nic_level);
#else
    uint8_t nic_scaling_level = get_nic_scaling_level (pd_pass , enc_mode, pcs_ptr->temporal_layer_index);

    set_nic_controls(context_ptr, nic_scaling_level);

    uint8_t nic_pruning_level;
    if (pd_pass == PD_PASS_0)
        nic_pruning_level = 0;
    else
        if (enc_mode <= ENC_MRS)
            nic_pruning_level = 0;
        else if (enc_mode <= ENC_MR)
            nic_pruning_level = 1;

        else if (enc_mode <= ENC_M0)
            nic_pruning_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 3;
#if TUNE_M0_M7_MEGA_FEB
#if TUNE_M5_M6
        else if (enc_mode <= ENC_M5)
#else
        else if (enc_mode <= ENC_M4)
#endif
#else
        else if (enc_mode <= ENC_M1)
#endif
            nic_pruning_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 4;
#if !TUNE_M0_M7_MEGA_FEB
        else if (enc_mode <= ENC_M2)
            nic_pruning_level = 4;
#endif
#if TUNE_M0_M7_MEGA_FEB
#if (TUNE_M10_M3_1 && !TUNE_M7_M8) || TUNE_M7_SLOWDOWN
        else if (enc_mode <= ENC_M7)
#else
        else if (enc_mode <= ENC_M6)
#endif
#else
        else if (enc_mode <= ENC_M4)
#endif
            nic_pruning_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 5;
#if !TUNE_M0_M7_MEGA_FEB
        else if (enc_mode <= ENC_M5)
            nic_pruning_level = 5;
        else if (enc_mode <= ENC_M6)
            nic_pruning_level = 6;
#endif
#if CLN_NIC_PRUNING
#if TUNE_M8_M9_FEB24
#if TUNE_M7_M10_MT
#if TUNE_M10_M0
#if CLN_M10_M12_DIFFS
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M10)
#endif
#else
        else if (enc_mode <= ENC_M9)
#endif
#else
        else if (enc_mode <= ENC_M8)
#endif
#else
        else if (enc_mode <= ENC_M7)
#endif
            nic_pruning_level = 7;
#if !TUNE_M10_FASTER
#if TUNE_M10_INITIAL_PRESET && (!TUNE_NEW_M10_M11 || TUNE_M7_M10_MT)
        else if (enc_mode <= ENC_M10)
#else
        else if (enc_mode <= ENC_M9)
#endif
            nic_pruning_level = 9;
#endif
#else
        else if (enc_mode <= ENC_M7)
            nic_pruning_level = (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 7 : 8;
        else if (enc_mode <= ENC_M9)
            nic_pruning_level = (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 9 : 10;
#endif
#if OPT_NIC_PRUNING
        else
            nic_pruning_level = 12;
#else
        else
            nic_pruning_level = 11;
#endif

    set_nic_pruning_controls(context_ptr, nic_pruning_level);
#endif
#endif
    // Set md_filter_intra_mode @ MD
    // md_filter_intra_level specifies whether filter intra would be active
    // for a given prediction candidate in mode decision.

    // md_filter_intra_level | Settings
    // 0                      | OFF
    // 1                      | ON
    if (pd_pass == PD_PASS_0)
        context_ptr->md_filter_intra_level = 0;
    else
        context_ptr->md_filter_intra_level =
        pcs_ptr->pic_filter_intra_level;
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

#if CLN_REMOVE_UNUSED_FEATS
    uint8_t pf_level = 1;
    set_pf_controls(context_ptr, pf_level);
#else
#if OPT_REDUCE_TX
    uint8_t pf_level = 1;

    // Only allow PF if using SB_64x64
    if (pcs_ptr->slice_type != I_SLICE && sequence_control_set_ptr->static_config.super_block_size != 128) {

        if (pd_pass == PD_PASS_0) {
#if TUNE_M3_M6_MEM_OPT
#if TUNE_MEGA_M9_M4
#if TUNE_M10_M7
            if (enc_mode <= ENC_M7 || pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) {
#else
            if (enc_mode <= ENC_M6 || pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) {
#endif
#else
            if (enc_mode <= ENC_M5 || pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) {
#endif
#else
            if (enc_mode <= ENC_M4 || pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) {
#endif
                pf_level = 1;
            }
            else {
                // Use ME distortion and variance detector to enable PF
                uint64_t use_pf_th = compute_pf_th(sequence_control_set_ptr, pcs_ptr, context_ptr);
                uint32_t fast_lambda = context_ptr->hbd_mode_decision ? context_ptr->fast_lambda_md[EB_10_BIT_MD] : context_ptr->fast_lambda_md[EB_8_BIT_MD];
                uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index]);
#if TUNE_M8_M9_FEB24
#if TUNE_M9_M10_MAR
#if TUNE_NEW_M10_M11 && !TUNE_M11
                if (enc_mode <= ENC_M11) {
#else
#if TUNE_NEW_M11_2
#if TUNE_NEW_M12
#if FTR_M13
                if (enc_mode <= ENC_M13) {
#else
                if (enc_mode <= ENC_M12) {
#endif
#else
                if (enc_mode <= ENC_M11) {
#endif
#else
                if (enc_mode <= ENC_M10) {
#endif
#endif
#else
                if (enc_mode <= ENC_M8) {
#endif
#else
                if (enc_mode <= ENC_M7) {
#endif
                    pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
                }
                else {
                    if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                        pf_level = (cost_64x64 < ((use_pf_th * 3) >> 1)) ? 3 : 1;
                    else
                        pf_level = (cost_64x64 < ((use_pf_th * 5) >> 1)) ? 3 : 1;
                }

            }
        }
        else { // PD1
#if TUNE_M9_M10_MAR
#if TUNE_M10_M9_1
#if TUNE_NEW_M11_2
#if TUNE_NEW_M12
#if FTR_M13
            if (enc_mode <= ENC_M13) {
#else
            if (enc_mode <= ENC_M12) {
#endif
#else
            if (enc_mode <= ENC_M11) {
#endif
#else
            if (enc_mode <= ENC_M10) {
#endif
#else
            if (enc_mode <= ENC_M9) {
#endif
#else
            if (enc_mode <= ENC_M8) {
#endif
                pf_level = 1;
            }
            else {
                // Use ME distortion and variance detector to enable PF
                uint64_t use_pf_th = compute_pf_th(sequence_control_set_ptr, pcs_ptr, context_ptr);
                uint32_t fast_lambda = context_ptr->hbd_mode_decision ? context_ptr->fast_lambda_md[EB_10_BIT_MD] : context_ptr->fast_lambda_md[EB_8_BIT_MD];
                uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index]);

                pf_level = (cost_64x64 < use_pf_th) ? 2 : 1;
            }
        }
    }
    else {
        pf_level = 1;
    }
    set_pf_controls(context_ptr, pf_level);
#else
    if (pcs_ptr->slice_type != I_SLICE) {

        if (pd_pass == PD_PASS_0) {

            // Only allow PF if using SB_64x64
#if TUNE_M0_M7_MEGA_FEB
            if (enc_mode <= ENC_M7 ||
#else
            if (enc_mode <= ENC_M4 ||
#endif
                pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ||
                sequence_control_set_ptr->static_config.super_block_size == 128) {
                context_ptr->pf_level = 1;
            }
            else {
                // Use ME distortion and variance detector to enable PF
                uint32_t fast_lambda = context_ptr->hbd_mode_decision ?
                    context_ptr->fast_lambda_md[EB_10_BIT_MD] :
                    context_ptr->fast_lambda_md[EB_8_BIT_MD];
                uint32_t sb_size = sequence_control_set_ptr->static_config.super_block_size * sequence_control_set_ptr->static_config.super_block_size;
                uint64_t cost_th_rate = 1 << 13;
                uint64_t use_pf_th = 0;

                if (pcs_ptr->parent_pcs_ptr->variance[context_ptr->sb_index][ME_TIER_ZERO_PU_64x64] <= 400)
                    use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
                else if (pcs_ptr->parent_pcs_ptr->variance[context_ptr->sb_index][ME_TIER_ZERO_PU_64x64] <= 800)
                    use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
                else
                    use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);

                uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index]);
#if !TUNE_M0_M7_MEGA_FEB
                if (enc_mode <= ENC_M7) {
                    context_ptr->pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
                }
                else {
#endif
                    if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                        context_ptr->pf_level = (cost_64x64 < ((use_pf_th * 3) >> 1)) ? 3 : 1;
                    else
                        context_ptr->pf_level = (cost_64x64 < ((use_pf_th * 5) >> 1)) ? 3 : 1;
#if !TUNE_M0_M7_MEGA_FEB
                }
#endif

            }
        }
        else
            context_ptr->pf_level = 1;
    }
    else {
        context_ptr->pf_level = 1;
    }
    set_pf_controls(context_ptr, context_ptr->pf_level);
#endif
#endif
#if !CLN_INDEPTH
    uint8_t in_depth_block_skip_level = 0;
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        in_depth_block_skip_level = 0;
    else if (context_ptr->pd_pass == PD_PASS_0)
#if TUNE_IN_DEPTH_LIST0_M9_LEVEL
#if TUNE_M9_M10_MAR
#if TUNE_NEW_M10_M11
        if (enc_mode <= ENC_M9)
#else
        if (enc_mode <= ENC_M10)
#endif
#else
        if (enc_mode <= ENC_M8)
#endif
            in_depth_block_skip_level = 0;
#if TUNE_NEW_M10_M11 && !TUNE_M10_M7
        else if (enc_mode <= ENC_M10)
            in_depth_block_skip_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? 1 : 0;
#endif
        else
            in_depth_block_skip_level = 1;
#else
        if (enc_mode <= ENC_M9)
            in_depth_block_skip_level = 0;
        else
            in_depth_block_skip_level = (pcs_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? 0 : 1;
#endif
    else
        in_depth_block_skip_level = 0;

    set_in_depth_block_skip_ctrls(context_ptr, in_depth_block_skip_level);
#endif
#if !CLN_REMOVE_UNUSED_FEATS
    uint8_t lower_depth_block_skip_level = 0;
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        lower_depth_block_skip_level = 0;
    else if (pd_pass == PD_PASS_0)
        lower_depth_block_skip_level = 0;
    else {
#if TUNE_M8_M9_FEB24
#if TUNE_M9_MARCH
#if TUNE_M10_M0
        if (enc_mode <= ENC_M10)
#else
        if (enc_mode <= ENC_M9)
#endif
#else
        if (enc_mode <= ENC_M8)
#endif
#else
        if (enc_mode <= ENC_M7)
#endif
            lower_depth_block_skip_level = 0;
#if TUNE_M7_M10_PRESETS
        else
            lower_depth_block_skip_level = 1;
#else
        else if (enc_mode <= ENC_M9)
            lower_depth_block_skip_level = 1;
        else
            lower_depth_block_skip_level = 2;
#endif
    }

    set_lower_depth_block_skip_ctrls(context_ptr, lower_depth_block_skip_level);
#endif
#if CLN_REG_PD_SIG_SET_2
    md_sq_motion_search_controls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->md_sq_mv_search_level);

    md_nsq_motion_search_controls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->md_nsq_mv_search_level);

    md_pme_search_controls(
        context_ptr,
        pd_pass == PD_PASS_0 ? 0 : pcs_ptr->md_pme_level);
#else
    if (pd_pass == PD_PASS_0)
        context_ptr->md_sq_mv_search_level = 0;
    else
#if TUNE_M0_M7_MEGA_FEB
        if (enc_mode <= ENC_M3)
#else
        if (enc_mode <= ENC_M1)
#endif
            context_ptr->md_sq_mv_search_level = 1;
#if !TUNE_M0_M7_MEGA_FEB
        else if (enc_mode <= ENC_M4)
            context_ptr->md_sq_mv_search_level = 4;
#endif
        else
            context_ptr->md_sq_mv_search_level = 0;
    md_sq_motion_search_controls(context_ptr, context_ptr->md_sq_mv_search_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->md_nsq_mv_search_level = 0;
    else
        if (enc_mode <= ENC_MRS)
            context_ptr->md_nsq_mv_search_level = 2;
        else
            context_ptr->md_nsq_mv_search_level = 4;

    md_nsq_motion_search_controls(context_ptr, context_ptr->md_nsq_mv_search_level);
    // Set PME level
    if (pd_pass == PD_PASS_0)
        context_ptr->md_pme_level = 0;
    else

        if (enc_mode <= ENC_M0)
#if TUNE_PME_M0
            context_ptr->md_pme_level = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 2 : 1;
#else
            context_ptr->md_pme_level = 1;
#endif
#if TUNE_MEGA_M9_M4
#if TUNE_M5_M6 && !TUNE_M1_M8
        else if (enc_mode <= ENC_M4)
#else
        else if (enc_mode <= ENC_M5)
#endif
#else
        else if (enc_mode <= ENC_M6)
#endif
#if TUNE_PME_M0
            context_ptr->md_pme_level = 3;
#else
            context_ptr->md_pme_level = 2;
#endif
#if TUNE_M9_11_3
#if TUNE_M1_M8
        else if (enc_mode <= ENC_M6)
#else
        else if (enc_mode <= ENC_M7)
#endif
            context_ptr->md_pme_level = 4;
#if TUNE_M9_M11_OPTIMIZED_SUPER4KL
#if TUNE_M10_M11
#if TUNE_IMPROVE_M11_M10
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M10)
#endif
#else
        else if (enc_mode <= ENC_M9)
#endif
#else
        else if (enc_mode <= ENC_M8)
#endif
            context_ptr->md_pme_level = 6;
#endif
#if !TUNE_M10_M11
#if TUNE_M7_M10_MT
#if TUNE_M10_M9_1
#if TUNE_M10_M0 && !TUNE_M9_11_3
        else if (enc_mode <= ENC_M10)
#else
        else if (enc_mode <= ENC_M9)
#endif
#else
        else if (enc_mode <= ENC_M8)
#endif
#else
        else if (enc_mode <= ENC_M7)
#endif
#if FTR_USE_PSAD
#if TUNE_M9_11_3
            context_ptr->md_pme_level = 8;
#else
            context_ptr->md_pme_level = 4;
#endif
#else
            context_ptr->md_pme_level = 3;
#endif
#endif
#if FTR_USE_PSAD
#if !TUNE_M7_M10_MT
        else if (enc_mode <= ENC_M8)
            context_ptr->md_pme_level = 5;
#endif
#if !TUNE_M10_M9_1
        else if (enc_mode <= ENC_M9)
            context_ptr->md_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE? 7 : 6;
#endif
#if !TUNE_M10_M0
        else if (enc_mode <= ENC_M10)
            context_ptr->md_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE? 7 : 8;
#endif
#else
#if TUNE_M7_M10_PRESETS && !TUNE_M8_M9_FEB24 || TUNE_M9_M10_MAR
        else if (enc_mode <= ENC_M9)
#else
        else if (enc_mode <= ENC_M8)
#endif
            context_ptr->md_pme_level = 4;
#if TUNE_NEW_M10_M11
        else if (enc_mode <= ENC_M10)
            context_ptr->md_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE? 4 : 5;
#endif
#endif
#if FTR_IMPROVE_PME
#if TUNE_M10_INITIAL_PRESET
#if TUNE_NEW_M10_M11
#if TUNE_REG_PD1
#if !TUNE_M10_M11
        else if (enc_mode <= ENC_M10)
            context_ptr->md_pme_level = 9;
#endif
#if TUNE_NEW_M12
        else if (enc_mode <= ENC_M12)
#else
        else if (enc_mode <= ENC_M11)
#endif
#else
        else if (enc_mode <= ENC_M11)
#endif
#else
        else if (enc_mode <= ENC_M10)
#endif
#else
        else if (enc_mode <= ENC_M9)
#endif
#if FTR_USE_PSAD
#if TUNE_M10_M3_1
#if OPT_M11_PME
#if TUNE_PME_M0
            context_ptr->md_pme_level = 10;
#else
            context_ptr->md_pme_level = 9;
#endif
#else
            context_ptr->md_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE? 0 : 8;
#endif
#else
            context_ptr->md_pme_level = 8;
#endif
#else
            context_ptr->md_pme_level = 5;
#endif
#endif
        else
            context_ptr->md_pme_level = 0;

    md_pme_search_controls(context_ptr, context_ptr->md_pme_level);
#endif
    if (pd_pass == PD_PASS_0)
#if TUNE_M0_M7_MEGA_FEB
        context_ptr->md_subpel_me_level = enc_mode <= ENC_M5 ? 3 : 0;
#else
        context_ptr->md_subpel_me_level = enc_mode <= ENC_M6 ? 3 : 0;
#endif
#if CLN_SUBPEL_LVL
    else if (enc_mode <= ENC_M0)
        context_ptr->md_subpel_me_level = 1;
#if TUNE_M1_M8
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M6)
#endif
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 1 : 2;
#else
    else
        if (enc_mode <= ENC_M4)
            context_ptr->md_subpel_me_level = 1;
#endif
#if TUNE_M8_M10 && !TUNE_4K_M8_M11
        else if (enc_mode <= ENC_M8)
#else
#if !CLN_SUBPEL_LVL
#if TUNE_M7_M8_3
        else if (enc_mode <= ENC_M6)
#else
        else if (enc_mode <= ENC_M7)
#endif
#endif
#endif
#if !CLN_SUBPEL_LVL
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 1 : 2;
#endif
#if TUNE_M7_M8_3
#if TUNE_M8_SLOWDOWN
        else if (enc_mode <= ENC_M8)
#else
        else if (enc_mode <= ENC_M7)
#endif
            context_ptr->md_subpel_me_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 4 : 7);
#endif
#if TUNE_4K_M8_M11
#if TUNE_M9_SLOWDOWN
#if CLN_M10_M12_DIFFS
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M9)
#endif
#else
        else if (enc_mode <= ENC_M8)
#endif
#if TUNE_M7_M8_3
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (4) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 10);
#else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 1 : (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 2 : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 10));
#endif
#endif
#if FTR_LOW_AC_SUBPEL && !TUNE_M8_M10
        else if (enc_mode <= ENC_M8)
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 7;
#endif
#if TUNE_CTR_QUARTER_PEL && !TUNE_M7_M10_MT
#if TUNE_M9_M10_MAR
        else if (enc_mode <= ENC_M9)
#else
        else if (enc_mode <= ENC_M8)
#endif
#if FTR_LOW_AC_SUBPEL
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 5 : 8;
#else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 6;
#endif
#endif
#if !TUNE_M9_M10_MAR
        else if (enc_mode <= ENC_M9)
#if TUNE_M8_M9_FEB24
#if TUNE_CTR_QUARTER_PEL
            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 6;
            else
                context_ptr->md_subpel_me_level = 7;
#else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 6;
#endif
#else
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 5;
#endif
#endif
#if TUNE_NEW_M10_M11 && !CLN_M10_M12_DIFFS
#if TUNE_M10_FASTER && !TUNE_M10_SLOWDOWN
        else if (enc_mode <= ENC_M9)
#else
#if TUNE_M10_M11
        else if (enc_mode <= ENC_M11)
#else
        else if (enc_mode <= ENC_M10)
#endif
#endif
#if FTR_LOW_AC_SUBPEL
#if OPT_SUBPEL
#if TUNE_M10_M0
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 7) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 10);
#else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 5 : 8) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 10);
#endif
#else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 5 : 8) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 9);
#endif
#else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 6) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 6 : 7);
#endif
#endif
#if !TUNE_M8_M10
#if OPT_M11_SUBPEL
        else if (enc_mode <= ENC_M10)
#else
        else
#endif
#endif
#if OPT_M12_SPEED
#if !TUNE_M10_M11
       else if (enc_mode <= ENC_M11)
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 12 : 13);
#endif
#endif
#if TUNE_CTR_QUARTER_PEL
#if FTR_LOW_AC_SUBPEL
#if OPT_SUBPEL
#if !TUNE_M8_M10
            context_ptr->md_subpel_me_level = (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11;
#endif
#if OPT_M11_SUBPEL
#if TUNE_4K_M11
#if TUNE_M9_11_3
#if TUNE_NEW_M12
       else if (enc_mode <= ENC_M12)
#else
       else if (enc_mode <= ENC_M11)
#endif
#else
       else if (enc_mode <= ENC_M10)
#endif
#if OPT_M12_SPEED
       context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11) :
       pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE ? ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 12 : 13) : ((pcs_ptr->temporal_layer_index == 0) ? 12 : 15);
#else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 12 : 13);
#endif
       else {
#if FTR_M13
           if (pcs_ptr->temporal_layer_index != 0)
               context_ptr->md_subpel_me_level = 0;
           else
               context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 12 : 13);
#else
           if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE)
               context_ptr->md_subpel_me_level = (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11;
           else if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
               context_ptr->md_subpel_me_level = (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 12 : 13;
           else
               context_ptr->md_subpel_me_level = (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 15 : 16;
#endif
       }
#else
        else
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 9 : 11) : ((pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 12 : 13);
#endif
#endif
#else
            context_ptr->md_subpel_me_level = (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 8 : 9;
#endif

#else
        context_ptr->md_subpel_me_level = (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) ? 6 : 7;
#endif
#else
        context_ptr->md_subpel_me_level = 6;
#endif
    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);
#if CLN_SUBPEL_LVL
    if (pd_pass == PD_PASS_0)
        context_ptr->md_subpel_pme_level = enc_mode <= ENC_M0 ? 3 : 0;
    else if (enc_mode <= ENC_M0)
        context_ptr->md_subpel_pme_level = 1;
    else
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 2 : 4;
#else
    if (pd_pass == PD_PASS_0)
        context_ptr->md_subpel_pme_level = enc_mode <= ENC_M4 ? 3 : 0;
#if TUNE_M10_M7
#if TUNE_M5_M6
#if TUNE_M10_M3_1
    else if (enc_mode <= ENC_M4)
#else
    else if (enc_mode <= ENC_M5)
#endif
#else
    else if (enc_mode <= ENC_M6)
#endif
#else
    else if (enc_mode <= ENC_M7)
#endif
        context_ptr->md_subpel_pme_level = 1;
#if TUNE_NEW_M10_M11
#if TUNE_M8_M10_4K_SUPER
#if TUNE_M9_11_3
    else if (enc_mode <= ENC_M8)
#else
    else if (enc_mode <= ENC_M9)
#endif
        context_ptr->md_subpel_pme_level = 2;
#if CLN_REG_PD1_DIFFS
    else
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 2 : 4;
#else
    else if (enc_mode <= ENC_M10)
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 2 : 4;
#endif
#else
#if TUNE_M11
#if OPT_M11_SUBPEL
    else if (enc_mode <= ENC_M10)
#else
    else if (enc_mode <= ENC_M11)
#endif
#else
    else if (enc_mode <= ENC_M10)
#endif
#if TUNE_M8_M10_4K_SUPER
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 2 : 4;
#else
        context_ptr->md_subpel_pme_level = 2;
#endif
#endif
#if !CLN_REG_PD1_DIFFS
#if FTR_M13
    else if (enc_mode <= ENC_M12)
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 3 : 4);
    else {
        if (pcs_ptr->temporal_layer_index != 0)
            context_ptr->md_subpel_pme_level = 0;
        else
            context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 3 : 4);
    }
#else
    else
#if OPT_M11_SUBPEL
#if TUNE_4K_M11
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : ((pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 3 : 4);
#else
        context_ptr->md_subpel_pme_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 2 : 3;
#endif
#else
        context_ptr->md_subpel_pme_level = 3;
#endif
#endif
#endif
#else
    else
        context_ptr->md_subpel_pme_level = 2;
#endif
#endif
    md_subpel_pme_controls(context_ptr, context_ptr->md_subpel_pme_level);
#if !CLN_INTRA_CTRLS
    // Set dc_cand_only_flag
    if (pd_pass == PD_PASS_0)
        context_ptr->dc_cand_only_flag = EB_TRUE;
    else
#if TUNE_M7_M10_PRESETS
#if TUNE_NEW_M10_M11
#if TUNE_NEW_M12
#if FTR_M13
        if (enc_mode <= ENC_M13)
#else
        if (enc_mode <= ENC_M12)
#endif
#else
        if (enc_mode <= ENC_M11)
#endif
#else
        if (enc_mode <= ENC_M10)
#endif
#else
        if (enc_mode <= ENC_M9)
#endif
            context_ptr->dc_cand_only_flag = EB_FALSE;
        else
            context_ptr->dc_cand_only_flag = EB_TRUE;

    // Set intra_angle_delta @ MD
    if (pd_pass == PD_PASS_0)
        context_ptr->md_intra_angle_delta = 0;
    else if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta == DEFAULT)
        context_ptr->md_intra_angle_delta = 1;
    else
        context_ptr->md_intra_angle_delta = pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.intra_angle_delta;

    // Set disable_angle_z2_prediction_flag
    if (pd_pass == PD_PASS_0)
        context_ptr->disable_angle_z2_intra_flag = EB_TRUE;
    else
        context_ptr->disable_angle_z2_intra_flag = EB_FALSE;
#endif
#if CLN_RATE_EST_CTRLS
    uint8_t rate_est_level = 0;
    if (pd_pass == PD_PASS_0) {
        if (enc_mode <= ENC_M2)
            rate_est_level = 1;
        else
            rate_est_level = 2;
    }
    else {
#if TUNE_M1_M8
        if (enc_mode <= ENC_M3)
#else
        if (enc_mode <= ENC_M2)
#endif
            rate_est_level = 1;
        else if (enc_mode <= ENC_M9)
            rate_est_level = 2;
        else if (enc_mode <= ENC_M11)
            rate_est_level = (pcs_ptr->slice_type == I_SLICE) ? 3 : 4;
        else
#if TUNE_IMPROVE_M12
            rate_est_level = pcs_ptr->slice_type == I_SLICE ? 3 : (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? 4 : 0;
#else
            rate_est_level = pcs_ptr->slice_type == I_SLICE ? 3 : (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? 0 : 4;
#endif
    }

    set_rate_est_ctrls(context_ptr, rate_est_level);

    // set at pic-level b/c feature depends on some pic-level initializations
    context_ptr->approx_inter_rate = pcs_ptr->approx_inter_rate;
#else
    // Shut skip_context and dc_sign update for rate estimation
    if (pd_pass == PD_PASS_0)
#if TUNE_M10_M0
        context_ptr->shut_skip_ctx_dc_sign_update = enc_mode <= ENC_M3 ? EB_FALSE : EB_TRUE;
#else
        context_ptr->shut_skip_ctx_dc_sign_update = enc_mode <= ENC_M4 ? EB_FALSE : EB_TRUE;
#endif
    else
#if TUNE_M8_M10_4K_SUPER
    {
#if TUNE_M9_SLOWDOWN
        if (enc_mode <= ENC_M9) {
#else
        if (enc_mode <= ENC_M8) {
#endif
            context_ptr->shut_skip_ctx_dc_sign_update = EB_FALSE;
        }
        else if (enc_mode <= ENC_M9) {
            context_ptr->shut_skip_ctx_dc_sign_update = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? ((pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE) : EB_FALSE;
        }
#if TUNE_NEW_M12
#if FTR_M13
        else if (enc_mode <= ENC_M13) {
#else
        else if (enc_mode <= ENC_M12) {
#endif
#else
        else if (enc_mode <= ENC_M11) {
#endif
            context_ptr->shut_skip_ctx_dc_sign_update = (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
        }
        else {
            context_ptr->shut_skip_ctx_dc_sign_update = EB_TRUE;
        }
    }
#else
#if TUNE_M10_M7
        context_ptr->shut_skip_ctx_dc_sign_update = enc_mode <= ENC_M8 ?
#else
        context_ptr->shut_skip_ctx_dc_sign_update = enc_mode <= ENC_M7 ?
#endif
        EB_FALSE :
#if TUNE_NEW_M10_M11
        enc_mode <= ENC_M11 ? ((pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE) :
        EB_TRUE;
#else
        (pcs_ptr->slice_type == I_SLICE) ?
        EB_FALSE :
        EB_TRUE;
#endif
#endif
#endif
    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    if (pd_pass == PD_PASS_0)
        context_ptr->shut_fast_rate = EB_TRUE;
    else
        context_ptr->shut_fast_rate = EB_FALSE;
#if !CLN_RATE_EST_CTRLS
#if FTR_ENHANCE_FAST_COEFF_RATE
    // Estimate the rate of the first(eob / fast_coeff_est_level) coeff(s), DC and last coeff only
    if (pd_pass == PD_PASS_0)
#if TUNE_M0_M7_MEGA_FEB
#if TUNE_M3_M6_MEM_OPT
        if (enc_mode <= ENC_M2)
#else
        if (enc_mode <= ENC_M3)
#endif
#else
        if (enc_mode <= ENC_M4)
#endif
            context_ptr->fast_coeff_est_level = 1;
#if TUNE_M10_M0
        else if (enc_mode <= ENC_M9)
#else
        else if (enc_mode <= ENC_M8)
#endif
            context_ptr->fast_coeff_est_level = 2;
#if !TUNE_M10_M0
#if TUNE_M8_M9_FEB24
        else if (enc_mode <= ENC_M8)
#else
        else if (enc_mode <= ENC_M9)
#endif
            context_ptr->fast_coeff_est_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
#endif
        else
#if TUNE_M9_TPL_ME_HME
            context_ptr->fast_coeff_est_level = (pcs_ptr->slice_type == I_SLICE) ? 2 : 4;
#else
            context_ptr->fast_coeff_est_level = 4;
#endif
    else
#if TUNE_NEW_M11_2
    {
        if (enc_mode <= ENC_M10)
            context_ptr->fast_coeff_est_level = 1;
#if TUNE_NEW_M12
#if FTR_M13
        else if (enc_mode <= ENC_M13)
#else
        else if (enc_mode <= ENC_M12)
#endif
#else
        else if (enc_mode <= ENC_M11)
#endif
            context_ptr->fast_coeff_est_level = 3;
#else
#if TUNE_M9_MARCH
#if TUNE_NEW_M10_M11
#if TUNE_M11
        if (enc_mode <= ENC_M11)
#else
        if (enc_mode <= ENC_M10)
#endif
#else
        if (enc_mode <= ENC_M9)
#endif
#else
        if (enc_mode <= ENC_M8)
#endif
            context_ptr->fast_coeff_est_level = 1;
#endif
        else
            context_ptr->fast_coeff_est_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
#if TUNE_NEW_M11_2
    }
#endif
#else
    // Estimate the rate of the first (eob/N) coeff(s) and last coeff only
    if (pd_pass == PD_PASS_0)
        if (enc_mode <= ENC_M4)
            context_ptr->fast_coeff_est_level = 0;
        else if (enc_mode <= ENC_M9)
            context_ptr->fast_coeff_est_level = 1;
        else
            context_ptr->fast_coeff_est_level = 2;
    else
        context_ptr->fast_coeff_est_level = 0;
#endif
#endif

#if SS_FIX_MOVE_SKIP_INTRA_PIC
#if CLN_INTRA_CTRLS
    // intra_level must be greater than 0 for I_SLICE
    uint8_t intra_level = 0;
    if (pd_pass == PD_PASS_0) {
        if (pcs_ptr->slice_type == I_SLICE)
            intra_level = 5;
        else if (enc_mode <= ENC_M1)
            intra_level = 5;
        else
            intra_level = (pcs_ptr->temporal_layer_index == 0) ? 5 : 0;
    }
    else if (enc_mode <= ENC_M0)
        intra_level = 1;
    else if (enc_mode <= ENC_M2)
        intra_level = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 4;
    else if (enc_mode <= ENC_M5)
        intra_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 4;
    else if (enc_mode <= ENC_M9)
        intra_level = pcs_ptr->slice_type == I_SLICE ? 1 : (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
    else
        intra_level = pcs_ptr->slice_type == I_SLICE ? 3 : 4;

    set_intra_ctrls(pcs_ptr, context_ptr, intra_level);
#else
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->skip_intra = 0;
    else if (pd_pass == PD_PASS_0){
        if (enc_mode <= ENC_M1)
            context_ptr->skip_intra = 0;
        else if (enc_mode <= ENC_M7)
            context_ptr->skip_intra = (pcs_ptr->temporal_layer_index == 0) ? 0 : 1;
        else
            context_ptr->skip_intra = 1;
    }
    else
        context_ptr->skip_intra = pcs_ptr->skip_intra;
#endif
#else
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->skip_intra = 0;
    else if (pd_pass == PD_PASS_0)
        if (enc_mode <= ENC_M1)
            context_ptr->skip_intra = 0;
#if TUNE_M0_M7_MEGA_FEB
#if TUNE_MEGA_M9_M4 && !LIGHT_PD0
        else if (enc_mode <= ENC_M8)
#else
#if TUNE_M10_M7 && !TUNE_M10_M0
        else if (enc_mode <= ENC_M6)
#else
        else if (enc_mode <= ENC_M7)
#endif
#endif
#else
        else if (enc_mode <= ENC_M4)
#endif
            context_ptr->skip_intra = (pcs_ptr->temporal_layer_index == 0) ? 0 : 1;
        else
            context_ptr->skip_intra = 1;
    else
#if TUNE_SKIP_INTRA_PD1
#if TUNE_M10_M9_1
#if TUNE_M9_SLOW
    {
        if (enc_mode <= ENC_M8)
            context_ptr->skip_intra = 0;
#if !TUNE_M10_M0
        else if (enc_mode <= ENC_M9)
            context_ptr->skip_intra = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1;
#endif
        else if (enc_mode <= ENC_M10) {
            context_ptr->skip_intra = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1;
            if (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag && pcs_ptr->intra_percentage > 100) {
                context_ptr->skip_intra = 0;
            }
        }
        else
#if TUNE_M9_11_3
            context_ptr->skip_intra = (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) ? (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1) : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : (pcs_ptr->intra_percentage > 100 ? 0 : 1));
#else
            context_ptr->skip_intra = pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1;
#endif
    }
#else
        context_ptr->skip_intra = (enc_mode <= ENC_M8) ? 0 : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1);
#endif
#else
        context_ptr->skip_intra = (enc_mode <= ENC_M9) ? 0 : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1);
#endif
#else
        context_ptr->skip_intra = 0;
#endif
#endif
#if !SS_OPT_MD
    if (pd_pass == PD_PASS_0)
        context_ptr->use_prev_mds_res = EB_FALSE;
    else
        context_ptr->use_prev_mds_res = EB_FALSE;
#endif
#if !CLN_MDS0_CTRLS
#if OPT_EARLY_CAND_ELIM
    uint8_t early_cand_elimination = 0;
    if (pd_pass == PD_PASS_0)
        early_cand_elimination = 0;
    else if (pcs_ptr->slice_type == I_SLICE)
        early_cand_elimination = 0;
#if TUNE_M7_M8
    else if (enc_mode <= ENC_M7)
#else
    else if (enc_mode <= ENC_M6)
#endif
        early_cand_elimination = 0;
#if TUNE_M9_SLOW
#if (TUNE_M7_11 && !TUNE_M9_M10) || TUNE_M10_SLOWDOWN
    else if (enc_mode <= ENC_M10)
#else
    else if (enc_mode <= ENC_M9)
#endif
#else
    else if (enc_mode <= ENC_M8)
#endif
        early_cand_elimination = 1;
#if !TUNE_M10_M0
    else if (enc_mode <= ENC_M10)
        early_cand_elimination = 2;
#endif
    else
        early_cand_elimination = 3;

    set_early_cand_elimination_controls(context_ptr, early_cand_elimination);
#else
    if (pd_pass == PD_PASS_0)
        context_ptr->early_cand_elimination = 0;
    else
        if (pcs_ptr->slice_type == I_SLICE)
            context_ptr->early_cand_elimination = 0;
        else
            if (enc_mode <= ENC_M6)
                context_ptr->early_cand_elimination = 0;
            else
                context_ptr->early_cand_elimination = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 120 : 102;
#endif
#endif
#if !CLN_REG_PD1_TX_CTRLS
#if OPT_REDUCE_TX
    /*reduce_last_md_stage_candidate
    0: OFF
    1: When MDS1 output has 0 coeffs, apply PFN4 and use DCT_DCT only
    2: 1 + disallow RDOQ and IFS when when MDS0 cand == MDS1 cand and
    the candidate does not belong to the best class
    3: 1 + 2 + remouve candidates when when MDS0 cand == MDS1 cand and they don't belong to the best class*/
#else
    /*reduce_last_md_stage_candidate
    0: OFF
    1: Aply PFN2 when the block is 0 coeff and PFN4 when MDS0 cand == MDS1 cand and
    the candidate does not belong to the best class
    2: 1 + disallow RDOQ and IFS when when MDS0 cand == MDS1 cand and
    the candidate does not belong to the best class
    3: 1 + 2 + remouve candidates when when MDS0 cand == MDS1 cand and they don't belong to the best class*/
#endif
    if (pd_pass == PD_PASS_0)
        context_ptr->reduce_last_md_stage_candidate = 0;
    else
        if (pcs_ptr->slice_type == I_SLICE)
            context_ptr->reduce_last_md_stage_candidate = 0;
        else
#if FTR_OPT_PF_PD1
#if TUNE_M7_MT
            if(enc_mode <= ENC_M7)
                context_ptr->reduce_last_md_stage_candidate = 0;
#else
            if(enc_mode <= ENC_M6)
                context_ptr->reduce_last_md_stage_candidate = 0;
#if TUNE_M8_M10
            else if (enc_mode <= ENC_M7)
#else
            else if (enc_mode <= ENC_M8)
#endif
                context_ptr->reduce_last_md_stage_candidate = 3;
#endif
            else
                context_ptr->reduce_last_md_stage_candidate = 4;
#else
#if TUNE_MEGA_M9_M4
            context_ptr->reduce_last_md_stage_candidate = (enc_mode <= ENC_M6) ? 0 : 3;
#else
            context_ptr->reduce_last_md_stage_candidate = (enc_mode <= ENC_M7) ? 0 : 3;
#endif
#endif
#endif
#if !CLN_CAND_REDUCTION_CTRLS
    if (pd_pass == PD_PASS_0)
        context_ptr->merge_inter_classes = 1;
    else
#if TUNE_MATCH_04_M8
#if TUNE_M7_11
        context_ptr->merge_inter_classes = (enc_mode <= ENC_M6) ? 0 : 1;
#else
        context_ptr->merge_inter_classes = (enc_mode <= ENC_M7) ? 0 : 1;
#endif
#else
        context_ptr->merge_inter_classes = (enc_mode <= ENC_M8) ? 0 : 1;
#endif
#endif
#if CLN_REG_PD_SIG_SET_2
    set_mds0_controls(
        pcs_ptr,
        context_ptr,
        pd_pass == PD_PASS_0 ? 2 : pcs_ptr->mds0_level);
#else
#if CLN_MDS0_CTRLS
    uint8_t mds0_level = 0;
    if (pd_pass == PD_PASS_0)
        mds0_level = 2;
#if CLN_M6_M12_FEATURES
    else if (enc_mode <= ENC_M11)
        mds0_level = 2;
#else
    else if (enc_mode <= ENC_M7)
        mds0_level = 2;
#if CLN_M10_M12_DIFFS
    else if (enc_mode <= ENC_M11)
#else
    else if (enc_mode <= ENC_M10)
#endif
        mds0_level = pcs_ptr->slice_type == I_SLICE ? 2 : 3;
#endif
    else
        mds0_level = pcs_ptr->slice_type == I_SLICE ? 2 : 4;

    set_mds0_controls(pcs_ptr, context_ptr, mds0_level);
#else
#if SS_OPT_MDS0
    context_ptr->mds0_dist_type = (enc_mode <= ENC_MRS) ? MDS0_SAD : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? MDS0_VAR : MDS0_SAD;
#else
    context_ptr->use_var_in_mds0 = (enc_mode <= ENC_MRS) ? 0 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 0;
#endif
#endif
#endif
#if !CLN_CAND_REDUCTION_CTRLS
    uint8_t eliminate_candidate_based_on_pme_me_results = 0;

    if (pd_pass == PD_PASS_0)
        eliminate_candidate_based_on_pme_me_results = 0;
    else
        if (pcs_ptr->slice_type == I_SLICE)
            eliminate_candidate_based_on_pme_me_results = 0;

        else if (enc_mode <= ENC_M6)
            eliminate_candidate_based_on_pme_me_results = 0;
#if CLN_REG_PD1_DIFFS
        else
            eliminate_candidate_based_on_pme_me_results = 1;
#else
#if TUNE_M10_M7 && !TUNE_M7_M10_MT
        else if (enc_mode <= ENC_M8)
#else
#if TUNE_M10_M0
#if TUNE_M9_M10
        else if (enc_mode <= ENC_M10)
#else
        else if (enc_mode <= ENC_M9)
#endif
#else
        else if (enc_mode <= ENC_M7)
#endif
#endif
            eliminate_candidate_based_on_pme_me_results = 1;
        else
            eliminate_candidate_based_on_pme_me_results = 2;
#endif
    set_cand_elimination_controls(context_ptr, eliminate_candidate_based_on_pme_me_results);
#endif
#if !TUNE_NEW_TXT_LVLS
    if (enc_mode <= ENC_M7)
        context_ptr->early_txt_search_exit_level = 0;
    else
        if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
            context_ptr->early_txt_search_exit_level = 1;
        else
            context_ptr->early_txt_search_exit_level = 2;
#endif
#if !OPT_PD0_PATH
#if TUNE_MEGA_M9_M4
    if (enc_mode <= ENC_M8)
#else
    if (enc_mode <= ENC_M7)
#endif
        context_ptr->ep_use_md_skip_decision = 0;
    else
        context_ptr->ep_use_md_skip_decision = 1;
#endif
#if CLN_SUBRES
    set_subres_controls(context_ptr, 0);
#else
#if FTR_SUBSAMPLE_RESIDUAL
    uint8_t subres_level;
#if TUNE_M10_M0 && !TUNE_4K_M8_M11
    if (enc_mode <= ENC_M9) {
#else
    if (enc_mode <= ENC_M8) {
#endif
        subres_level = 0;
    }
#if TUNE_4K_M8_M11
    else if (enc_mode <= ENC_M9 && pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE) {
        subres_level = 0;
    }
#endif
    else {
        subres_level = 0;
#if FTR_TX_SUB_PD0_INTRA
        if (sequence_control_set_ptr->static_config.super_block_size != 128 /* 128x128 not tested*/ && context_ptr->disallow_4x4 /* as no T4x2*/&& context_ptr->md_disallow_nsq /* as no T4x2*/) {
#else
        if (pcs_ptr->slice_type != I_SLICE && sequence_control_set_ptr->static_config.super_block_size != 128 /* 128x128 not tested*/ && context_ptr->disallow_4x4 /* as no T4x2*/&& context_ptr->md_disallow_nsq /* as no T4x2*/) {
#endif

            SbParams *sb_params_ptr = &pcs_ptr->parent_pcs_ptr->sb_params_array[context_ptr->sb_index];

            // The controls checks the deviation between: (1) the pred-to-src SAD of even rows and (2) the pred-to-src SAD of odd rows for each 64x64 to decide whether to use subres or not
            // then applies the result to the 64x64 block and to all children, therefore if incomplete 64x64 then shut subres
            if (sb_params_ptr->is_complete_sb) {

                // Use ME distortion and variance detector to enable subres
                uint64_t use_subres_th = compute_subres_th(sequence_control_set_ptr, pcs_ptr, context_ptr);
                uint32_t fast_lambda = context_ptr->hbd_mode_decision ? context_ptr->fast_lambda_md[EB_10_BIT_MD] : context_ptr->fast_lambda_md[EB_8_BIT_MD];
                uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index]);

                if (pd_pass == PD_PASS_0) {
                    if (pcs_ptr->slice_type == I_SLICE)
#if FTR_TX_SUB_PD0_INTRA
                        subres_level = 1;
#else
                        subres_level = 0;
#endif
                    else
                        subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;

                    context_ptr->subres_ctrls.odd_to_even_deviation_th = 5;
                }
                else {
#if PD1_SUBSAMPLE_RESIDUAL
                    if (pcs_ptr->temporal_layer_index == 0)
                        subres_level = 0;
#if TUNE_SUBRES_TX
#if CLN_REG_PD1_DIFFS
                    else
                        subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
#else
                    else if (enc_mode <= ENC_M10)
                        subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
                    else
#if TUNE_ADD_SUBRESS_FACTOR4
                        subres_level = (context_ptr->depth_removal_ctrls.enabled &&
                            (context_ptr->depth_removal_ctrls.disallow_below_32x32 ||
                                context_ptr->depth_removal_ctrls.disallow_below_64x64)) ? 2 : 1;
#else
                        subres_level = 1;
#endif
#endif
#else
                    else
                        subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
#endif
#else
                    subres_level = 0;
#endif
                    context_ptr->subres_ctrls.odd_to_even_deviation_th = 0;
                }
#if TUNE_ADD_SUBRESS_FACTOR4
                // If level 2 is used, 8x8 blocks must be disallowed
                if (subres_level == 2 &&
                    !(context_ptr->depth_removal_ctrls.enabled &&
                    (context_ptr->depth_removal_ctrls.disallow_below_16x16 ||
                        context_ptr->depth_removal_ctrls.disallow_below_32x32 ||
                        context_ptr->depth_removal_ctrls.disallow_below_64x64))) {
                    subres_level = 1;
                }
#endif
            }
        }
    }
    set_subres_controls(context_ptr, subres_level);
#endif
#endif
#if !LIGHT_PD1_MACRO
    context_ptr->use_best_mds0 = 0;
    if (pd_pass == PD_PASS_0) {
#if TUNE_M10_M0
        if (enc_mode <= ENC_M9)
#else
        if (enc_mode <= ENC_M8)
#endif
            context_ptr->use_best_mds0 = 0;
        else
            context_ptr->use_best_mds0 = 1;
    }
#endif
#if FTR_SKIP_MDS1 && !CLN_NIC_PRUNE_CTRLS
    uint8_t mds1_skip_level = 0;
    if (pd_pass == PD_PASS_0)
        mds1_skip_level = 0;
    else if (pcs_ptr->parent_pcs_ptr->sc_class1)
        mds1_skip_level = 0;
#if CLN_DECPL_TX_FEATS
#if TUNE_M9_11_3 && !TUNE_M9_M11_OPTIMIZED_SUPER4KL
    else if (enc_mode <= ENC_M8)
#else
    else if (enc_mode <= ENC_M9)
#endif
        mds1_skip_level = 0;
#if CLN_REG_PD1_TX_CTRLS
    else
        mds1_skip_level = 1;
#else
#if FTR_TX_NEIGH_INFO
    else if(enc_mode <= ENC_M10)
        mds1_skip_level = pcs_ptr->slice_type == I_SLICE ? 1 : 2;
#if FTR_M13
    else if (enc_mode <= ENC_M12)
        mds1_skip_level = pcs_ptr->slice_type == I_SLICE ? 1 : 3;
    else
#if CLN_LPD1_LVLS
        mds1_skip_level = pcs_ptr->slice_type == I_SLICE ? 1 : 4;
#else
        mds1_skip_level = pcs_ptr->slice_type == I_SLICE ? 1 : 6;
#endif
#else
    else
        mds1_skip_level = pcs_ptr->slice_type == I_SLICE ? 1 : 3;
#endif
#else
    else
        mds1_skip_level = pcs_ptr->slice_type == I_SLICE ? 1 : 2;
#endif
#endif
#else
#if TUNE_M10_M0
    else if (enc_mode <= ENC_M9)
#else
    else if (enc_mode <= ENC_M8)
#endif
        mds1_skip_level = 0;
#if TUNE_M10_M9_1
#if TUNE_NEW_M11_2
    else if (enc_mode <= ENC_M11)
#else
    else if (enc_mode <= ENC_M10)
#endif
#else
    else if (enc_mode <= ENC_M9)
#endif
        mds1_skip_level = 1;
#if !TUNE_M10_M9_1
    else if (enc_mode <= ENC_M10)
        mds1_skip_level = 2;
#endif
    else
        mds1_skip_level = 3;
#endif
    set_mds1_skipping_controls(context_ptr, mds1_skip_level);
#endif
#if !CLN_REMOVE_UNUSED_FEATS
#if FTR_PD_EARLY_EXIT
    if (pd_pass == PD_PASS_0)
        context_ptr->pd0_early_exit_th = enc_mode <= ENC_M8 ? 0 : 3;
    else
        context_ptr->pd0_early_exit_th = 0;
#endif
#endif
#if FTR_PD_EARLY_EXIT
#if CLN_REMOVE_UNUSED_FEATS
    context_ptr->pd0_inter_depth_bias = 0;
#else
    if (pd_pass == PD_PASS_0)
        context_ptr->pd0_inter_depth_bias = enc_mode <= ENC_M8 ? 0 : 1003;
    else
        context_ptr->pd0_inter_depth_bias = 0;
#endif
#endif
#if REMOVE_CLOSE_MVS
#if CLN_M10_M12_DIFFS
#if !CLN_CAND_REDUCTION_CTRLS
    uint8_t redundant_cand_level = 0;
    set_redundant_cand_controls(context_ptr, redundant_cand_level);
#endif
#else
    uint8_t redundant_cand_level = 0;
    if (pd_pass == PD_PASS_0)
        redundant_cand_level = 0;
    else
#if TUNE_M10_M0 && !TUNE_M7_11
        if (enc_mode <= ENC_M10)
#else
        if (enc_mode <= ENC_M9)
#endif
            redundant_cand_level = 0;
        else
            redundant_cand_level = 1;

    set_redundant_cand_controls(context_ptr, redundant_cand_level);
#endif
#endif
#if !SS_CLN_CFL_CTRLS
#if FTR_FASTER_CFL
#if TUNE_M10_M0
#if TUNE_M7_M8_3
    context_ptr->cfl_itr_th = enc_mode <= ENC_M6 ? 2 : 1;
#else
    context_ptr->cfl_itr_th = enc_mode <= ENC_M8 ? 2 : 1;
#endif
#else
    context_ptr->cfl_itr_th = enc_mode <= ENC_M9 ? 2 : 1;
#endif
#endif
#endif
#if CLN_REG_PD1_DIFFS
#if !CLN_CAND_REDUCTION_CTRLS
    context_ptr->reduce_unipred_candidates = 0;
#endif
#else
#if FTR_REDUCE_UNI_PRED
#if TUNE_M7_11 && !TUNE_M8_M10_4K_SUPER
    context_ptr->reduce_unipred_candidates = enc_mode <= ENC_M9 ? 0 : 1;
#else
    context_ptr->reduce_unipred_candidates = enc_mode <= ENC_M10 ? 0 : 1;
#endif
#endif
#endif
#if !CLN_ADD_LIST0_ONLY_CTRL
#if OPT_USE_INTRA_NEIGHBORING
    uint8_t use_neighbouring_mode = 0;
    if (pd_pass == PD_PASS_0)
        use_neighbouring_mode = 0;
#if TUNE_M8_M10_4K_SUPER
    else if (enc_mode <= ENC_M8)
        use_neighbouring_mode = 0;
    else if (enc_mode <= ENC_M9)
        use_neighbouring_mode = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE ? 0 : 1;
    else
        use_neighbouring_mode = 1;
#else
#if TUNE_M7_11
    else if (enc_mode <= ENC_M9)
#else
    else if (enc_mode <= ENC_M10)
#endif
        use_neighbouring_mode = 0;
    else
        use_neighbouring_mode = 1;
#endif
    set_use_neighbouring_mode_ctrls(context_ptr, use_neighbouring_mode);
#endif
#endif
#if !CLN_RATE_EST_CTRLS
#if FIX_SKIP_COEFF_CONTEXT
    if (pd_pass == PD_PASS_0)
        context_ptr->use_skip_coeff_context = (enc_mode <= ENC_M0) ? 1 : 0;
    else
        context_ptr->use_skip_coeff_context = (enc_mode <= ENC_M0) ? 1 : 0;
#endif
#endif
    return return_error;
}
#else

EbErrorType signal_derivation_enc_dec_kernel_oq(
    SequenceControlSet *sequence_control_set_ptr,
    PictureControlSet *pcs_ptr,
    ModeDecisionContext *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    EbEncMode enc_mode = pcs_ptr->enc_mode;
    uint8_t pd_pass = context_ptr->pd_pass;
    uint8_t txt_level = 0;
    if (pd_pass == PD_PASS_0)
        txt_level = 0;
    else if (pd_pass == PD_PASS_1)
        txt_level = 0;
    else
        if (enc_mode <= ENC_M2)
            txt_level = 1;
        else if (enc_mode <= ENC_M4)
            txt_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 3;
        else if (enc_mode <= ENC_M5)
            txt_level = 3;
        else if (enc_mode <= ENC_M8)
            txt_level = 5;
        else
            txt_level = (pcs_ptr->slice_type == I_SLICE) ? 5 : 0;
    set_txt_controls(context_ptr, txt_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->interpolation_search_level = IFS_OFF;
    else if (pd_pass == PD_PASS_1)
        context_ptr->interpolation_search_level = IFS_OFF;
    else
        if (enc_mode <= ENC_MR)
            context_ptr->interpolation_search_level = IFS_MDS1;
        else if (enc_mode <= ENC_M8)
            context_ptr->interpolation_search_level = IFS_MDS3;
        else
            context_ptr->interpolation_search_level = (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? IFS_MDS3 : IFS_OFF;
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
    if (enc_mode <= ENC_MRS) {
        context_ptr->chroma_at_last_md_stage = 0;
        context_ptr->chroma_at_last_md_stage_intra_th = (uint64_t)~0;
        context_ptr->chroma_at_last_md_stage_cfl_th = (uint64_t)~0;
    }
    else if (enc_mode <= ENC_M1) {
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
    if (enc_mode <= ENC_M5)
        context_ptr->md_disable_cfl = EB_FALSE;
    else
        context_ptr->md_disable_cfl = (pcs_ptr->temporal_layer_index == 0) ? EB_FALSE : EB_TRUE;
     // Set disallow_4x4
    context_ptr->disallow_4x4 = get_disallow_4x4 (enc_mode, pcs_ptr->slice_type);
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
         context_ptr->md_disallow_nsq = enc_mode <= ENC_M0 ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
     else if (pd_pass == PD_PASS_1)
         context_ptr->md_disallow_nsq = enc_mode <= ENC_MR ? pcs_ptr->parent_pcs_ptr->disallow_nsq : 1;
     else
         // Update nsq settings based on the sb_class
         context_ptr->md_disallow_nsq = pcs_ptr->parent_pcs_ptr->disallow_nsq;


     if (pd_pass == PD_PASS_0)
         context_ptr->global_mv_injection = 0;
     else if (pd_pass == PD_PASS_1)
         context_ptr->global_mv_injection = 0;
     else
         context_ptr->global_mv_injection = pcs_ptr->parent_pcs_ptr->gm_ctrls.enabled;

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
            if (enc_mode <= ENC_M0)
                    context_ptr->new_nearest_near_comb_injection = 1;
                else
                    context_ptr->new_nearest_near_comb_injection = 0;
        else
            context_ptr->new_nearest_near_comb_injection =
            sequence_control_set_ptr->static_config.new_nearest_comb_inject;

    uint8_t near_count_level = 0;

    if (pd_pass == PD_PASS_0)
        near_count_level = 0;
    if (pd_pass == PD_PASS_1)
        near_count_level = 0;
    else
        if (enc_mode <= ENC_M8)
            near_count_level = 1;
        else
            near_count_level = 2;
    set_near_count_ctrls(context_ptr, near_count_level);


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
        if (enc_mode <= ENC_M0)
            context_ptr->unipred3x3_injection = 1;
        else if (enc_mode <= ENC_M1)
            context_ptr->unipred3x3_injection = 2;
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
        if (enc_mode <= ENC_M1)
            context_ptr->bipred3x3_injection = 1;
        else if (enc_mode <= ENC_M5)
            context_ptr->bipred3x3_injection = 2;
        else
            context_ptr->bipred3x3_injection = 0;
        }
    else{
        context_ptr->bipred3x3_injection =
        sequence_control_set_ptr->static_config.bipred_3x3_inject;
        }

    context_ptr->inject_inter_candidates = 1;

    if (sequence_control_set_ptr->compound_mode) {
        if (sequence_control_set_ptr->static_config.compound_level == DEFAULT) {
            if (pd_pass == PD_PASS_0)
                context_ptr->inter_compound_mode = 0;
            else if (pd_pass == PD_PASS_1)
                context_ptr->inter_compound_mode = 0;
            else if (enc_mode <= ENC_MR)
                context_ptr->inter_compound_mode = 1;
            else if (enc_mode <= ENC_M3)

                context_ptr->inter_compound_mode = 2;
            else if (enc_mode <= ENC_M4)
                context_ptr->inter_compound_mode = 3;
            else if (enc_mode <= ENC_M5)
                context_ptr->inter_compound_mode = (pcs_ptr->parent_pcs_ptr->input_resolution == INPUT_SIZE_240p_RANGE) ? 5 : 0;
            else
                context_ptr->inter_compound_mode = 0;
        }
        else {
            context_ptr->inter_compound_mode = sequence_control_set_ptr->static_config.compound_level;
        }
    }
    else {
        context_ptr->inter_compound_mode = 0;
    }
    // Set dist_based_ref_pruning
    if (pcs_ptr->parent_pcs_ptr->ref_list0_count_try > 1 || pcs_ptr->parent_pcs_ptr->ref_list1_count_try > 1) {
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->dist_based_ref_pruning = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->dist_based_ref_pruning = 0;
        else if (enc_mode <= ENC_MR)
            context_ptr->dist_based_ref_pruning = 1;
        else if (enc_mode <= ENC_M0)
            context_ptr->dist_based_ref_pruning = (pcs_ptr->temporal_layer_index == 0) ? 1 : 2;
        else if (enc_mode <= ENC_M7)
            context_ptr->dist_based_ref_pruning = (pcs_ptr->temporal_layer_index == 0) ? 2 : 4;
        else
            context_ptr->dist_based_ref_pruning = 4;
    }
    else {
        context_ptr->dist_based_ref_pruning = 0;
    }

    set_dist_based_ref_pruning_controls(context_ptr, context_ptr->dist_based_ref_pruning);
    if (pd_pass == PD_PASS_0) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;
    }
    else if (pd_pass == PD_PASS_1) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    }
    else
        if (enc_mode <= ENC_M4)
            context_ptr->md_staging_mode = MD_STAGING_MODE_2;
        else
            context_ptr->md_staging_mode = MD_STAGING_MODE_1;


    // spatial_sse_full_loop_level | Default Encoder Settings            | Command Line Settings
    //             0               | OFF subject to possible constraints | OFF in PD_PASS_2
    //             1               | ON subject to possible constraints  | ON in PD_PASS_2
    if (pd_pass == PD_PASS_0)
        context_ptr->spatial_sse_full_loop_level = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->spatial_sse_full_loop_level = EB_FALSE;
    else if (sequence_control_set_ptr->static_config.spatial_sse_full_loop_level == DEFAULT)
        if (enc_mode <= ENC_M8)
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
        context_ptr->rdoq_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->rdoq_level = 0;
    else
        if (enc_mode <= ENC_M8)
            context_ptr->rdoq_level = 1;
        else
            context_ptr->rdoq_level = (pcs_ptr->parent_pcs_ptr->slice_type == I_SLICE) ? 2 : 3;
    set_rdoq_controls(context_ptr, context_ptr->rdoq_level);

    // Derive redundant block
    if (pd_pass == PD_PASS_0)
        context_ptr->redundant_blk = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->redundant_blk = EB_TRUE;
    else
        if (sequence_control_set_ptr->static_config.enable_redundant_blk ==
            DEFAULT)
            context_ptr->redundant_blk = EB_TRUE;
        else
            context_ptr->redundant_blk =
            sequence_control_set_ptr->static_config.enable_redundant_blk;
    uint8_t parent_sq_coeff_area_based_cycles_reduction_level = 0;
    if (pd_pass == PD_PASS_0)
        parent_sq_coeff_area_based_cycles_reduction_level = 0;
    else if (pd_pass == PD_PASS_1)
        parent_sq_coeff_area_based_cycles_reduction_level = 0;
    else if (enc_mode <= ENC_MRS)
        parent_sq_coeff_area_based_cycles_reduction_level = 0;
    else if (enc_mode <= ENC_MR)
        parent_sq_coeff_area_based_cycles_reduction_level = pcs_ptr->slice_type == I_SLICE ? 0 : 1;
    else if (enc_mode <= ENC_M1)
        parent_sq_coeff_area_based_cycles_reduction_level = pcs_ptr->slice_type == I_SLICE ? 0 : 2;
    else if (enc_mode <= ENC_M2)
        parent_sq_coeff_area_based_cycles_reduction_level = pcs_ptr->slice_type == I_SLICE ? 0
        : (pcs_ptr->temporal_layer_index == 0)? 2 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 4 : 7;
    else
        parent_sq_coeff_area_based_cycles_reduction_level = pcs_ptr->slice_type == I_SLICE ? 5 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 6 : 7;
    set_parent_sq_coeff_area_based_cycles_reduction_ctrls(context_ptr, pcs_ptr->parent_pcs_ptr->input_resolution, parent_sq_coeff_area_based_cycles_reduction_level);
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
                if (enc_mode <= ENC_M0)
                    context_ptr->sq_weight = 105;
                    else if (enc_mode <= ENC_M2)
                        context_ptr->sq_weight = 95;
                    else
                        context_ptr->sq_weight = 90;

        // max_part0_to_part1_dev is used to:
        // (1) skip the H_Path if the deviation between the Parent-SQ src-to-recon distortion of (1st quadrant + 2nd quadrant) and the Parent-SQ src-to-recon distortion of (3rd quadrant + 4th quadrant) is less than TH,
        // (2) skip the V_Path if the deviation between the Parent-SQ src-to-recon distortion of (1st quadrant + 3rd quadrant) and the Parent-SQ src-to-recon distortion of (2nd quadrant + 4th quadrant) is less than TH.
        if (pd_pass == PD_PASS_0)
            context_ptr->max_part0_to_part1_dev = 0;
        else if (pd_pass == PD_PASS_1)
            context_ptr->max_part0_to_part1_dev = 0;
        else
            if (enc_mode <= ENC_M3)
                context_ptr->max_part0_to_part1_dev = 0;
            else
                context_ptr->max_part0_to_part1_dev = 100;

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
    //Block level switch, has to follow the picture level
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
            context_ptr->md_inter_intra_level = 1;
        else
            context_ptr->md_inter_intra_level = 0;
    }
    else
        context_ptr->md_inter_intra_level = 0;

    set_inter_intra_ctrls(context_ptr, context_ptr->md_inter_intra_level);
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
        context_ptr->md_staging_tx_size_level = 1;
    else if (enc_mode <= ENC_MR)
        context_ptr->md_staging_tx_size_level = 2;
    else if (enc_mode <= ENC_M1)
        context_ptr->md_staging_tx_size_level = (pcs_ptr->temporal_layer_index == 0) ? 2 : 0;
    else
        context_ptr->md_staging_tx_size_level = 0;
    uint8_t nic_scaling_level = get_nic_scaling_level (pd_pass , enc_mode, pcs_ptr->temporal_layer_index);

    set_nic_controls(context_ptr, nic_scaling_level);

    uint8_t nic_pruning_level;
    if (pd_pass == PD_PASS_0)
        nic_pruning_level = 0;
    else if (pd_pass == PD_PASS_1)
        nic_pruning_level = 0;
    else
        if (enc_mode <= ENC_MRS)
            nic_pruning_level = 0;
        else if (enc_mode <= ENC_MR)
            nic_pruning_level = 1;
        else if (enc_mode <= ENC_M0)
            nic_pruning_level = (pcs_ptr->temporal_layer_index == 0) ? 1 : 3;
        else if (enc_mode <= ENC_M4)
            nic_pruning_level = (pcs_ptr->temporal_layer_index == 0) ? 3 : 4;
        else if (enc_mode <= ENC_M6)
            nic_pruning_level = (pcs_ptr->temporal_layer_index == 0) ? 4 : 5;
        else if (enc_mode <= ENC_M7)
            nic_pruning_level = (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 7 : 8;
        else if (enc_mode <= ENC_M8)
            nic_pruning_level = (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 9 : 10;
        else
            nic_pruning_level = 11;
    set_nic_pruning_controls(context_ptr, nic_pruning_level);
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
    if (pcs_ptr->slice_type != I_SLICE) {
        if (pd_pass == PD_PASS_0) {

            // Only allow PF if using SB_64x64
            if (enc_mode <= ENC_M6 ||
                pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ||
                sequence_control_set_ptr->static_config.super_block_size == 128) {
                context_ptr->pf_level = 1;
            }
            else {
                // Use ME distortion and variance detector to enable PF
                uint32_t fast_lambda = context_ptr->hbd_mode_decision ?
                    context_ptr->fast_lambda_md[EB_10_BIT_MD] :
                    context_ptr->fast_lambda_md[EB_8_BIT_MD];
                uint32_t sb_size = sequence_control_set_ptr->static_config.super_block_size * sequence_control_set_ptr->static_config.super_block_size;
                uint64_t cost_th_rate = 1 << 13;
                uint64_t use_pf_th = 0;

                if (pcs_ptr->parent_pcs_ptr->variance[context_ptr->sb_index][ME_TIER_ZERO_PU_64x64] <= 400)
                    use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 2);
                else if (pcs_ptr->parent_pcs_ptr->variance[context_ptr->sb_index][ME_TIER_ZERO_PU_64x64] <= 800)
                    use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size);
                else
                    use_pf_th = RDCOST(fast_lambda, cost_th_rate, sb_size >> 1);

                uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs_ptr->parent_pcs_ptr->me_64x64_distortion[context_ptr->sb_index]);
                if (enc_mode <= ENC_M8) {
                    context_ptr->pf_level = (cost_64x64 < use_pf_th) ? 3 : 1;
                }
                else {
                    if (pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
                        context_ptr->pf_level = (cost_64x64 < ((use_pf_th * 3) >> 1)) ? 3 : 1;
                    else
                        context_ptr->pf_level = (cost_64x64 < ((use_pf_th * 5) >> 1)) ? 3 : 1;
                }
            }
        }
        else if (pd_pass == PD_PASS_1)
            context_ptr->pf_level = 1;
        else
            context_ptr->pf_level = 1;
    } else {
        context_ptr->pf_level = 1;
    }
    set_pf_controls(context_ptr, context_ptr->pf_level);
#if !CLN_INDEPTH
    uint8_t in_depth_block_skip_level = 0;
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        in_depth_block_skip_level = 0;
    else if (context_ptr->pd_pass == PD_PASS_0)
        in_depth_block_skip_level = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        in_depth_block_skip_level = 0;
    else
        in_depth_block_skip_level = 0;

    set_in_depth_block_skip_ctrls(context_ptr, in_depth_block_skip_level);
#endif
    uint8_t lower_depth_block_skip_level = 0;
    if (pcs_ptr->parent_pcs_ptr->sc_class1)
        lower_depth_block_skip_level = 0;
    else if (pd_pass == PD_PASS_0)
        lower_depth_block_skip_level = 0;
    else if (pd_pass == PD_PASS_1)
        lower_depth_block_skip_level = 0;
    else {
        if (enc_mode <= ENC_M7)
            lower_depth_block_skip_level = 0;
        else
            lower_depth_block_skip_level = 1;
    }

    set_lower_depth_block_skip_ctrls(context_ptr, lower_depth_block_skip_level);

    if (pd_pass == PD_PASS_0)
        context_ptr->md_sq_mv_search_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_sq_mv_search_level = 0;
    else
        if (enc_mode <= ENC_M3)
            context_ptr->md_sq_mv_search_level = 1;
        else if (enc_mode <= ENC_M5)
            context_ptr->md_sq_mv_search_level = 4;
        else
            context_ptr->md_sq_mv_search_level = 0;
    md_sq_motion_search_controls(context_ptr, context_ptr->md_sq_mv_search_level);
    if (pd_pass == PD_PASS_0)
        context_ptr->md_nsq_mv_search_level = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_nsq_mv_search_level = 0;
    else
        if (enc_mode <= ENC_MRS)
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
        if (enc_mode <= ENC_M0)
            context_ptr->md_pme_level = 1;
        else if (enc_mode <= ENC_M5)
            context_ptr->md_pme_level = 2;
        else if (enc_mode <= ENC_M7)
            context_ptr->md_pme_level = 3;
        else if (enc_mode <= ENC_M8)
            context_ptr->md_pme_level = 4;
        else
            context_ptr->md_pme_level = 0;
    md_pme_search_controls(context_ptr, context_ptr->md_pme_level);

    if (pd_pass == PD_PASS_0)
        context_ptr->md_subpel_me_level = enc_mode <= ENC_M5 ? 3 : 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_subpel_me_level = 3;
    else
        if (enc_mode <= ENC_M4)
            context_ptr->md_subpel_me_level = 1;
        else if (enc_mode <= ENC_M6)
            context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 1 : 2;
        else if (enc_mode <= ENC_M8)
        context_ptr->md_subpel_me_level = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 4 : 5;
        else
        context_ptr->md_subpel_me_level = 6;

    md_subpel_me_controls(context_ptr, context_ptr->md_subpel_me_level);

    if (pd_pass == PD_PASS_0)
        context_ptr->md_subpel_pme_level = enc_mode <= ENC_M4 ? 3 : 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->md_subpel_pme_level = 3;
    else if (enc_mode <= ENC_M7)
        context_ptr->md_subpel_pme_level = 1;
    else
        context_ptr->md_subpel_pme_level = 2;
    md_subpel_pme_controls(context_ptr, context_ptr->md_subpel_pme_level);
    // Set dc_cand_only_flag
    if (pd_pass == PD_PASS_0)
        context_ptr->dc_cand_only_flag = EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->dc_cand_only_flag =
        (pcs_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
    else
        if (enc_mode <= ENC_M8)
            context_ptr->dc_cand_only_flag = EB_FALSE;
        else
            context_ptr->dc_cand_only_flag = EB_TRUE;

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
    // Shut skip_context and dc_sign update for rate estimation
    if (pd_pass == PD_PASS_0)
        context_ptr->shut_skip_ctx_dc_sign_update = enc_mode <= ENC_M4 ? EB_FALSE : EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->shut_skip_ctx_dc_sign_update = EB_TRUE;
    else
        context_ptr->shut_skip_ctx_dc_sign_update = enc_mode <= ENC_M7 ?
        EB_FALSE :
        (pcs_ptr->slice_type == I_SLICE) ?
        EB_FALSE :
        EB_TRUE;
    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    if (pd_pass == PD_PASS_0)
        context_ptr->shut_fast_rate = EB_TRUE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->shut_fast_rate = EB_FALSE;
    else
        context_ptr->shut_fast_rate = EB_FALSE;
    // Estimate the rate of the first (eob/N) coeff(s) and last coeff only
    if (pd_pass == PD_PASS_0)
        if (enc_mode <= ENC_M3)
            context_ptr->fast_coeff_est_level = 0;
        else if (enc_mode <= ENC_M8)
            context_ptr->fast_coeff_est_level = 1;
        else
            context_ptr->fast_coeff_est_level = 2;
    else if (pd_pass == PD_PASS_1)
        context_ptr->fast_coeff_est_level = 0;
    else
        context_ptr->fast_coeff_est_level = 0;
    if (pcs_ptr->slice_type == I_SLICE)
        context_ptr->skip_intra = 0;
    else if (pd_pass == PD_PASS_0)
        if (enc_mode <= ENC_M1)
            context_ptr->skip_intra = 0;
        else if (enc_mode <= ENC_M7)
            context_ptr->skip_intra = (pcs_ptr->temporal_layer_index == 0) ? 0 : 1;
        else
            context_ptr->skip_intra = 1;
    else
        context_ptr->skip_intra = 0;
    if (pd_pass == PD_PASS_0)
        context_ptr->use_prev_mds_res = EB_FALSE;
    else if (pd_pass == PD_PASS_1)
        context_ptr->use_prev_mds_res = EB_FALSE;
    else
        context_ptr->use_prev_mds_res = EB_FALSE;
    if (pd_pass == PD_PASS_0)
        context_ptr->early_cand_elimination = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->early_cand_elimination = 0;
    else
        if (pcs_ptr->slice_type == I_SLICE)
            context_ptr->early_cand_elimination = 0;
        else
            if (enc_mode <= ENC_M6)
                context_ptr->early_cand_elimination = 0;
            else
                context_ptr->early_cand_elimination = pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE ? 120 : 102;

    /*reduce_last_md_stage_candidate
    0: OFF
    1: Aply PFN2 when the block is 0 coeff and PFN4 when MDS0 cand == MDS1 cand and
    the candidate does not belong to the best class
    2: 1 + disallow RDOQ and IFS when when MDS0 cand == MDS1 cand and
    the candidate does not belong to the best class
    3: 1 + 2 + remouve candidates when when MDS0 cand == MDS1 cand and they don't belong to the best class*/
     if (pd_pass == PD_PASS_0)
        context_ptr->reduce_last_md_stage_candidate = 0;
    else if (pd_pass == PD_PASS_1)
        context_ptr->reduce_last_md_stage_candidate = 0;
    else
        if (pcs_ptr->slice_type == I_SLICE)
            context_ptr->reduce_last_md_stage_candidate = 0;
        else
            context_ptr->reduce_last_md_stage_candidate = (enc_mode <= ENC_M7) ? 0 : 3;
     if (pd_pass == PD_PASS_0)
         context_ptr->merge_inter_classes = 1;
     else if (pd_pass == PD_PASS_1)
         context_ptr->merge_inter_classes = 1;
     else
         context_ptr->merge_inter_classes = (enc_mode <= ENC_M8) ? 0 : 1;
     context_ptr->use_var_in_mds0 = (enc_mode <= ENC_MRS) ? 0 : pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 1 : 0;
     uint8_t eliminate_candidate_based_on_pme_me_results = 0;
     if (pd_pass == PD_PASS_0)
        eliminate_candidate_based_on_pme_me_results = 0;
    else if (pd_pass == PD_PASS_1)
        eliminate_candidate_based_on_pme_me_results = 0;
    else
        if (pcs_ptr->slice_type == I_SLICE)
            eliminate_candidate_based_on_pme_me_results = 0;
        else if (enc_mode <= ENC_M6)
            eliminate_candidate_based_on_pme_me_results = 0;
        else if (enc_mode <= ENC_M7)
            eliminate_candidate_based_on_pme_me_results = 1;
        else
            eliminate_candidate_based_on_pme_me_results = 2;
     set_cand_elimination_controls(context_ptr,eliminate_candidate_based_on_pme_me_results);
     if (pd_pass == PD_PASS_0)
         context_ptr->bypass_tx_search_when_zcoef = 0;
     else if (pd_pass == PD_PASS_1)
         context_ptr->bypass_tx_search_when_zcoef = 0;
     else
         if (pcs_ptr->slice_type == I_SLICE)
             context_ptr->bypass_tx_search_when_zcoef = 0;
         else
             context_ptr->bypass_tx_search_when_zcoef = (enc_mode <= ENC_M4) ? 0 : 1;
     if (enc_mode <= ENC_M8)
         context_ptr->early_txt_search_exit_level = 0;
     else
         if(pcs_ptr->parent_pcs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE)
             context_ptr->early_txt_search_exit_level = 1;
         else
             context_ptr->early_txt_search_exit_level = 2;
     if (enc_mode <= ENC_M7)
         context_ptr->ep_use_md_skip_decision = 0;
     else
         context_ptr->ep_use_md_skip_decision = 1;
     context_ptr->use_best_mds0 = 0;
     if (pd_pass == PD_PASS_0) {
         if (enc_mode <= ENC_M8)
             context_ptr->use_best_mds0 = 0;
         else
             context_ptr->use_best_mds0 = 1;
     }
    return return_error;
}
#endif
#if OPT_PD0_PATH
void copy_neighbour_arrays_light_pd0(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                                     uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds, uint32_t sb_org_x,
                                     uint32_t sb_org_y);
#endif
void copy_neighbour_arrays(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                           uint32_t src_idx, uint32_t dst_idx, uint32_t blk_mds, uint32_t sb_org_x,
                           uint32_t sb_org_y);

static void set_parent_to_be_considered(MdcSbData *results_ptr, uint32_t blk_index, int32_t sb_size,
                                        int8_t pred_depth,
                                        uint8_t pred_sq_idx,
                                        const uint8_t disallow_nsq,
                                        int8_t depth_step) {
    const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
    if (blk_geom->sq_size < ((sb_size == BLOCK_128X128) ? 128 : 64)) {
        //Set parent to be considered
#if LIGHT_PD0
        uint32_t parent_depth_idx_mds = blk_geom->parent_depth_idx_mds;
#else
        uint32_t parent_depth_idx_mds =
            (blk_geom->sqi_mds -
             (blk_geom->quadi - 3) * ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth]) -
            parent_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
#endif
        const BlockGeom *parent_blk_geom = get_blk_geom_mds(parent_depth_idx_mds);
        const uint32_t parent_tot_d1_blocks = disallow_nsq ? 1 :
            parent_blk_geom->sq_size == 128 ? 17 :
            parent_blk_geom->sq_size > 8 ? 25 :
            parent_blk_geom->sq_size == 8 ? 5 : 1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < parent_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[parent_depth_idx_mds + block_1d_idx] = 1;
        }

        if (depth_step < -1)
            set_parent_to_be_considered(results_ptr, parent_depth_idx_mds, sb_size, pred_depth, pred_sq_idx, disallow_nsq ,depth_step + 1);
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
#if CLN_GEOM
    if(blk_geom->geom_idx==GEOM_0)
        tot_d1_blocks = 1;
#endif

#if FIX_PRED_DEPTH_REFINEMNT_4X4
    if (blk_geom->sq_size == 8 && context_ptr->disallow_4x4) return;
#endif
    if (blk_geom->sq_size > 4) {
        for (uint32_t block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[blk_index + block_1d_idx]= 1;
            results_ptr->refined_split_flag[blk_index + block_1d_idx] = EB_TRUE;
        }
        //Set first child to be considered
#if LIGHT_PD0
        uint32_t child_block_idx_1 = blk_index + blk_geom->d1_depth_offset;
#else
        uint32_t child_block_idx_1 = blk_index +
            d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
#endif
        const BlockGeom *child1_blk_geom = get_blk_geom_mds(child_block_idx_1);
        const uint32_t child1_tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1 :
            child1_blk_geom->sq_size == 128 ? 17:
            child1_blk_geom->sq_size > 8 ? 25 :
            child1_blk_geom->sq_size == 8 ? 5 :
            1;

        for (uint32_t block_1d_idx = 0; block_1d_idx < child1_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[child_block_idx_1 + block_1d_idx] = 1;
            results_ptr->refined_split_flag[child_block_idx_1 + block_1d_idx] = EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_1])
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_1, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
        //Set second child to be considered
#if CLN_GEOM
        uint32_t child_block_idx_2 = child_block_idx_1 +
            ns_depth_offset[blk_geom->geom_idx][blk_geom->depth + 1];
#else
        uint32_t child_block_idx_2 = child_block_idx_1 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
#endif
        const BlockGeom *child2_blk_geom = get_blk_geom_mds(child_block_idx_2);
        const uint32_t child2_tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1 :
            child2_blk_geom->sq_size == 128 ? 17:
            child2_blk_geom->sq_size > 8 ? 25 :
            child2_blk_geom->sq_size == 8 ? 5 :
            1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < child2_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[child_block_idx_2 + block_1d_idx] = 1;
            results_ptr->refined_split_flag[child_block_idx_2 + block_1d_idx] = EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_2])
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_2, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
        //Set third child to be considered
#if CLN_GEOM
        uint32_t child_block_idx_3 = child_block_idx_2 +
            ns_depth_offset[blk_geom->geom_idx][blk_geom->depth + 1];
#else
        uint32_t child_block_idx_3 = child_block_idx_2 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
#endif
        const BlockGeom *child3_blk_geom = get_blk_geom_mds(child_block_idx_3);
        const uint32_t child3_tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1 :
            child3_blk_geom->sq_size == 128 ? 17:
            child3_blk_geom->sq_size > 8 ? 25 :
            child3_blk_geom->sq_size == 8 ? 5 :
            1;

        for (uint32_t block_1d_idx = 0; block_1d_idx < child3_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[child_block_idx_3 + block_1d_idx] = 1;
            results_ptr->refined_split_flag[child_block_idx_3 + block_1d_idx] = EB_FALSE;
        }

        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_3])
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_3, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
        //Set forth child to be considered
#if CLN_GEOM
        uint32_t child_block_idx_4 = child_block_idx_3 +
            ns_depth_offset[blk_geom->geom_idx][blk_geom->depth + 1];
#else
        uint32_t child_block_idx_4 = child_block_idx_3 +
            ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
#endif
        const BlockGeom *child4_blk_geom = get_blk_geom_mds(child_block_idx_4);
        const uint32_t child4_tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1 :
            child4_blk_geom->sq_size == 128 ? 17:
            child4_blk_geom->sq_size > 8 ? 25 :
            child4_blk_geom->sq_size == 8 ? 5 :
            1;
        for (uint32_t block_1d_idx = 0; block_1d_idx < child4_tot_d1_blocks; block_1d_idx++) {
            results_ptr->consider_block[child_block_idx_4 + block_1d_idx] = 1;
            results_ptr->refined_split_flag[child_block_idx_4 + block_1d_idx] = EB_FALSE;
        }
        // Add children blocks if more depth to consider (depth_step is > 1), or block not allowed (add next depth)
        if (depth_step > 1 || !pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[child_block_idx_4])
            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, child_block_idx_4, sb_index, sb_size, pred_depth, pred_sq_idx , depth_step > 1 ? depth_step - 1 : 1);
    }
}
#if OPT_MEM_PALETTE
uint32_t get_tot_1d_blks(struct PictureParentControlSet * ppcs, const int32_t sq_size, const uint8_t disallow_nsq) {
#else
static INLINE uint32_t get_tot_1d_blks(struct PictureParentControlSet * ppcs, const int32_t sq_size, const uint8_t disallow_nsq) {
#endif
    uint32_t tot_d1_blocks;

    tot_d1_blocks = (disallow_nsq) ||
        (sq_size >= 64 && ppcs->disallow_all_nsq_blocks_above_64x64) ||
        (sq_size >= 32 && ppcs->disallow_all_nsq_blocks_above_32x32) ||
        (sq_size >= 16 && ppcs->disallow_all_nsq_blocks_above_16x16) ||
        (sq_size <= 64 && ppcs->disallow_all_nsq_blocks_below_64x64) ||
        (sq_size <= 32 && ppcs->disallow_all_nsq_blocks_below_32x32) ||
        (sq_size <= 8 && ppcs->disallow_all_nsq_blocks_below_8x8) ||
        (sq_size <= 16 && ppcs->disallow_all_nsq_blocks_below_16x16) ? 1 :
        (sq_size == 16 && ppcs->disallow_all_non_hv_nsq_blocks_below_16x16) ? 5 :
        (sq_size == 16 && ppcs->disallow_all_h4_v4_blocks_below_16x16) ? 17 :
        sq_size == 128
        ? 17
        : sq_size > 8 ? 25 : sq_size == 8 ? 5 : 1;

    if (ppcs->disallow_HVA_HVB_HV4)
        tot_d1_blocks = MIN(5, tot_d1_blocks);

    if (ppcs->disallow_HV4)
        tot_d1_blocks = MIN(17, tot_d1_blocks);

    return tot_d1_blocks;
}

#if OPT_MEM_PALETTE
EbErrorType rtime_alloc_palette_info(BlkStruct * md_blk_arr_nsq) {

        EB_MALLOC_ARRAY(md_blk_arr_nsq->palette_info, 1);
        EB_MALLOC_ARRAY( md_blk_arr_nsq->palette_info->color_idx_map, MAX_PALETTE_SQUARE);

    return EB_ErrorNone;
}
#endif
// Initialize structures used to indicate which blocks will be tested at MD.
// MD data structures should be updated in init_block_data(), not here.
static void build_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
    uint32_t sb_index, EbBool is_complete_sb) {

    memset(context_ptr->tested_blk_flag, 0, sizeof(uint8_t) * scs_ptr->max_block_cnt);
#if !CLN_REMOVE_UNUSED_FEATS
    memset(context_ptr->do_not_process_blk, 0, sizeof(uint8_t) * scs_ptr->max_block_cnt);
#endif
#if OPT_PD0_PATH
    memset(context_ptr->avail_blk_flag, EB_FALSE, sizeof(uint8_t) * scs_ptr->max_block_cnt);
#endif

    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
    results_ptr->leaf_count = 0;
    uint32_t blk_index = 0;
#if !LIGHT_PD0
    const BlockSize sb_size = scs_ptr->seq_header.sb_size;
#endif
    const uint16_t max_block_cnt = scs_ptr->max_block_cnt;
    int32_t min_sq_size;
#if FTR_BYPASS_ENCDEC
    if (context_ptr->pred_depth_only)
        min_sq_size = (context_ptr->depth_removal_ctrls.enabled && context_ptr->depth_removal_ctrls.disallow_below_64x64)
        ? 64
        : (context_ptr->depth_removal_ctrls.enabled && context_ptr->depth_removal_ctrls.disallow_below_32x32)
        ? 32
        : (context_ptr->depth_removal_ctrls.enabled && context_ptr->depth_removal_ctrls.disallow_below_16x16)
        ? 16
        : context_ptr->disallow_4x4 ? 8 : 4;
    else
#endif
        min_sq_size = context_ptr->disallow_4x4 ? 8 : 4;

    while (blk_index < max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

        // Initialize here because may not be updated at inter-depth decision for incomplete SBs
        if (!is_complete_sb)
            context_ptr->md_blk_arr_nsq[blk_index].part = PARTITION_SPLIT;

        // SQ/NSQ block(s) filter based on the SQ size
        uint8_t is_block_tagged =
            (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE)
#if FTR_BYPASS_ENCDEC
            ||
            (context_ptr->pred_depth_only && (blk_geom->sq_size < min_sq_size))
#endif
            ? 0
            : 1;
#if FTR_M13
        if (context_ptr->skip_pd0)
            is_block_tagged = !context_ptr->depth_removal_ctrls.disallow_below_64x64 && context_ptr->depth_removal_ctrls.disallow_below_32x32 && (blk_geom->sq_size != 32) ? 0 : is_block_tagged;
#endif

        // SQ/NSQ block(s) filter based on the block validity
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_block_tagged) {

            const uint32_t tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1 :
                get_tot_1d_blks(pcs_ptr->parent_pcs_ptr, blk_geom->sq_size, context_ptr->md_disallow_nsq);

            for (uint32_t idx = blk_index; idx < (tot_d1_blocks + blk_index); ++idx) {

#if OPT_MEM_PALETTE
                //  MD palette info buffer
                if (pcs_ptr->parent_pcs_ptr->palette_level) {
                    if (context_ptr->md_blk_arr_nsq[idx].palette_mem == 0) {
                        rtime_alloc_palette_info(&context_ptr->md_blk_arr_nsq[idx]);
                            context_ptr->md_blk_arr_nsq[idx].palette_mem = 1;
                    }
                }

                    context_ptr->md_blk_arr_nsq[idx].palette_size[0] = 0;
                    context_ptr->md_blk_arr_nsq[idx].palette_size[1] = 0;
#endif

                if (results_ptr->consider_block[idx]) {
                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = idx;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = tot_d1_blocks;
                    results_ptr->split_flag[results_ptr->leaf_count++] = results_ptr->refined_split_flag[idx];
                }
            }
#if LIGHT_PD0
            blk_index += blk_geom->d1_depth_offset;
#else
            blk_index += d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
#endif
        }
        else {
#if LIGHT_PD0
            blk_index +=
                (blk_geom->sq_size > min_sq_size)
                ? blk_geom->d1_depth_offset
                : blk_geom->ns_depth_offset;
#else
            blk_index +=
                (blk_geom->sq_size > min_sq_size)
                ? d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
#endif
        }
    }
}
#if OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS
void update_pred_th_offset(ModeDecisionContext *ctx, const BlockGeom *blk_geom, int8_t *s_depth, int8_t *e_depth, int64_t *th_offset) {

    uint32_t full_lambda = ctx->hbd_mode_decision ?
        ctx->full_lambda_md[EB_10_BIT_MD] :
        ctx->full_lambda_md[EB_8_BIT_MD];

    // cost-band-based modulation
    uint64_t max_cost = RDCOST(full_lambda, 16, ctx->depth_refinement_ctrls.max_cost_multiplier * blk_geom->bwidth * blk_geom->bheight);

    if (ctx->md_local_blk_unit[blk_geom->sqi_mds].default_cost <= max_cost) {
        uint64_t band_size = max_cost / ctx->depth_refinement_ctrls.max_band_cnt;
        uint64_t band_idx = ctx->md_local_blk_unit[blk_geom->sqi_mds].default_cost / band_size;
        if (ctx->depth_refinement_ctrls.decrement_per_band[band_idx] == MAX_SIGNED_VALUE) {
            *s_depth = 0;
            *e_depth = 0;
        }
        else {
            *th_offset = -ctx->depth_refinement_ctrls.decrement_per_band[band_idx];
        }
    }
    else {
        *th_offset = 0;
    }
}
#else
void update_pred_th_offset(ModeDecisionContext *mdctxt, const BlockGeom *blk_geom, int8_t *s_depth, int8_t *e_depth, int64_t *th_offset) {

    uint32_t full_lambda = mdctxt->hbd_mode_decision ?
        mdctxt->full_lambda_md[EB_10_BIT_MD] :
        mdctxt->full_lambda_md[EB_8_BIT_MD];

    uint64_t cost_th_0 = (RDCOST(full_lambda, 16, 200 * blk_geom->bwidth * blk_geom->bheight) << (mdctxt->depth_refinement_ctrls.use_pred_block_cost - 1));
    uint64_t cost_th_1 = (RDCOST(full_lambda, 16, 300 * blk_geom->bwidth * blk_geom->bheight) << (mdctxt->depth_refinement_ctrls.use_pred_block_cost - 1));
    uint64_t cost_th_2 = (RDCOST(full_lambda, 16, 400 * blk_geom->bwidth * blk_geom->bheight) << (mdctxt->depth_refinement_ctrls.use_pred_block_cost - 1));

    if (mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost < cost_th_0) {
        *s_depth = 0;
        *e_depth = 0;
    }
    else if (mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost < cost_th_1) {
        *th_offset = -10;
    }
    else if (mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost < cost_th_2) {
        *th_offset = -5;
    }
}
#endif
#if LIGHT_PD0
uint8_t is_parent_to_current_deviation_small(ModeDecisionContext *mdctxt,
    const BlockGeom *blk_geom, int64_t th_offset) {
#else
uint8_t is_parent_to_current_deviation_small(SequenceControlSet *scs_ptr,
    ModeDecisionContext *mdctxt, const BlockGeom *blk_geom, int64_t th_offset) {
#endif
    if (mdctxt->depth_refinement_ctrls.parent_to_current_th == MIN_SIGNED_VALUE)
        return EB_FALSE;
    mdctxt->parent_to_current_deviation = MIN_SIGNED_VALUE;
    // block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
        // Get the parent of the current block
#if LIGHT_PD0
    uint32_t parent_depth_idx_mds = blk_geom->parent_depth_idx_mds;
#else
    uint32_t parent_depth_idx_mds =
        (blk_geom->sqi_mds -
        (blk_geom->quadi - 3) * ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]) -
        parent_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
#endif
    if (mdctxt->avail_blk_flag[parent_depth_idx_mds]) {
        mdctxt->parent_to_current_deviation =
            (int64_t)(((int64_t)MAX(mdctxt->md_local_blk_unit[parent_depth_idx_mds].default_cost, 1) - (int64_t)MAX((mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost * 4), 1)) * 100) /
            (int64_t)MAX((mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost * 4), 1);
    }
    if (mdctxt->parent_to_current_deviation <= (mdctxt->depth_refinement_ctrls.parent_to_current_th + th_offset))
        return EB_TRUE;

    return EB_FALSE;
}

uint8_t is_child_to_current_deviation_small(SequenceControlSet *scs_ptr,
    ModeDecisionContext *mdctxt, const BlockGeom *blk_geom, uint32_t blk_index, int64_t th_offset) {

    if (mdctxt->depth_refinement_ctrls.sub_to_current_th == MIN_SIGNED_VALUE)
        return EB_FALSE;
    mdctxt->child_to_current_deviation = MIN_SIGNED_VALUE;
#if LIGHT_PD0
    const uint32_t ns_d1_offset = blk_geom->d1_depth_offset;
#else
    const uint32_t ns_d1_offset = d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
#endif

#if CLN_GEOM
    (void)scs_ptr;
    const uint32_t ns_depth_plus1_offset = ns_depth_offset[blk_geom->geom_idx][blk_geom->depth + 1];
#else
    const uint32_t ns_depth_plus1_offset = ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth + 1];
#endif
    const uint32_t child_block_idx_1 = blk_index + ns_d1_offset;
    const uint32_t child_block_idx_2 = child_block_idx_1 + ns_depth_plus1_offset;
    const uint32_t child_block_idx_3 = child_block_idx_2 + ns_depth_plus1_offset;
    const uint32_t child_block_idx_4 = child_block_idx_3 + ns_depth_plus1_offset;

    uint64_t child_cost = 0;
    uint8_t child_cnt = 0;
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
            (int64_t)(((int64_t)MAX(child_cost, 1) - (int64_t)MAX(mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost, 1)) * 100) /
            (int64_t)(MAX(mdctxt->md_local_blk_unit[blk_geom->sqi_mds].default_cost, 1));
    }


    if (mdctxt->child_to_current_deviation <= (mdctxt->depth_refinement_ctrls.sub_to_current_th + th_offset))
        return EB_TRUE;

    return EB_FALSE;
}
static void perform_pred_depth_refinement(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                          ModeDecisionContext *context_ptr, uint32_t sb_index) {
    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
    uint32_t   blk_index   = 0;
    if (pcs_ptr->parent_pcs_ptr->disallow_nsq) {
        if (context_ptr->disallow_4x4) {
            memset(results_ptr->consider_block, 0, sizeof(uint8_t)*scs_ptr->max_block_cnt);
            memset(results_ptr->split_flag, 1, sizeof(uint8_t)*scs_ptr->max_block_cnt);
            memset(results_ptr->refined_split_flag, 1, sizeof(uint8_t)*scs_ptr->max_block_cnt);
        } else {
            while (blk_index < scs_ptr->max_block_cnt) {

                const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);

                    EbBool split_flag =
                    blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
                results_ptr->consider_block[blk_index] = 0;
                results_ptr->split_flag[blk_index] = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
                    results_ptr->refined_split_flag[blk_index] = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
#if LIGHT_PD0
                blk_index +=
                    split_flag
                    ? blk_geom->d1_depth_offset
                    : blk_geom->ns_depth_offset;
#else
                blk_index +=
                    split_flag
                    ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]
                    : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
#endif
            }
        }
    }
    else {
        // Reset mdc_sb_array data to defaults; it will be updated based on the predicted blocks (stored in md_blk_arr_nsq)
        while (blk_index < scs_ptr->max_block_cnt) {
            const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
            results_ptr->consider_block[blk_index] = 0;
            results_ptr->split_flag[blk_index] = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
            results_ptr->refined_split_flag[blk_index] = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
            blk_index++;
        }
    }
    results_ptr->leaf_count = 0;
    blk_index               = 0;
    while (blk_index < scs_ptr->max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
        const unsigned   tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1 :
            blk_geom->sq_size == 128 ? 17 :
            blk_geom->sq_size > 8 ? 25 :
            blk_geom->sq_size == 8 ? 5 :
            1;

        // if the parent square is inside inject this block
        uint8_t is_blk_allowed =
            pcs_ptr->slice_type != I_SLICE ? 1 : (blk_geom->sq_size < 128) ? 1 : 0;

        // derive split_flag
        EbBool split_flag = context_ptr->md_blk_arr_nsq[blk_index].split_flag;

        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_blk_arr_nsq[blk_index].split_flag == EB_FALSE) {
#if CLN_DEPTH_REFINEMENT
                    // Add current pred depth block(s)
                    for (unsigned block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
                        results_ptr->consider_block[blk_index + block_1d_idx] = 1;
                        results_ptr->refined_split_flag[blk_index + block_1d_idx] = EB_FALSE;
                    }

                    int8_t s_depth = context_ptr->depth_ctrls.s_depth;
                    int8_t e_depth = context_ptr->depth_ctrls.e_depth;

#if FTR_M13
                    if (context_ptr->skip_pd0) {
                        SbParams *sb_params = &pcs_ptr->parent_pcs_ptr->sb_params_array[context_ptr->sb_index];
                        if ((sb_params->width % 32 == 0) && (sb_params->height % 32 == 0)) {
                            s_depth = 0;
                            e_depth = 0;
                        }
                    }
#endif
                    // If multiple depths are selected, perform refinement
                    if (s_depth != 0 || e_depth != 0) {
                        // Check that the start and end depth are in allowed range, given other features
                        // which restrict allowable depths
                        if (context_ptr->disallow_4x4) {
                            e_depth = (blk_geom->sq_size == 8) ? 0
                                : (blk_geom->sq_size == 16) ? MIN(1, e_depth)
                                : (blk_geom->sq_size == 32) ? MIN(2, e_depth)
                                : e_depth;
                        }
                        if (context_ptr->depth_removal_ctrls.enabled) {
                            if (context_ptr->depth_removal_ctrls.disallow_below_64x64) {
                                e_depth = (blk_geom->sq_size <= 64) ? 0
                                    : (blk_geom->sq_size == 128) ? MIN(1, e_depth) : e_depth;
                            }
                            else if (context_ptr->depth_removal_ctrls.disallow_below_32x32) {
                                e_depth = (blk_geom->sq_size <= 32) ? 0
                                    : (blk_geom->sq_size == 64) ? MIN(1, e_depth)
                                    : (blk_geom->sq_size == 128) ? MIN(2, e_depth) : e_depth;
                            }
                            else if (context_ptr->depth_removal_ctrls.disallow_below_16x16) {

                                e_depth = (blk_geom->sq_size <= 16) ? 0
                                    : (blk_geom->sq_size == 32) ? MIN(1, e_depth)
                                    : (blk_geom->sq_size == 64) ? MIN(2, e_depth)
                                    : (blk_geom->sq_size == 128) ? MIN(3, e_depth) : e_depth;
                            }
                        }

                        uint8_t sq_size_idx = 7 - (uint8_t)svt_log2f((uint8_t)blk_geom->sq_size);
                        int64_t th_offset = 0;
                        if (context_ptr->depth_refinement_ctrls.enabled && context_ptr->depth_refinement_ctrls.cost_band_based_modulation && (s_depth != 0 || e_depth != 0)) {
                            update_pred_th_offset(context_ptr, blk_geom, &s_depth, &e_depth, &th_offset);
                        }

                        // Add block indices of upper depth(s)
                        // Block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
                        uint8_t add_parent_depth = 1;
                        if (context_ptr->depth_refinement_ctrls.enabled && s_depth == -1 && pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[blk_index] && blk_geom->sq_size < ((scs_ptr->seq_header.sb_size == BLOCK_128X128) ? 128 : 64)) {
#if LIGHT_PD0
                            add_parent_depth = is_parent_to_current_deviation_small(context_ptr, blk_geom, th_offset);
#else
                            add_parent_depth = is_parent_to_current_deviation_small(
                                scs_ptr, context_ptr, blk_geom, th_offset);
#endif
                        }

                        // Add block indices of lower depth(s)
                        // Block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
                        uint8_t add_sub_depth = 1;
                        if (context_ptr->depth_refinement_ctrls.enabled && e_depth == 1 && pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[blk_index] && blk_geom->sq_size > 4) {
                            add_sub_depth = is_child_to_current_deviation_small(
                                scs_ptr, context_ptr, blk_geom, blk_index, th_offset);
                        }

                        // Use a maximum of 2 depth per block (PRED+Parent or PRED+Sub)
                        if (context_ptr->depth_refinement_ctrls.enabled && context_ptr->depth_refinement_ctrls.up_to_2_depth) {
                            if ((s_depth == -1) && add_parent_depth && (e_depth == 1) && add_sub_depth) {
                                if (context_ptr->parent_to_current_deviation != MIN_SIGNED_VALUE && context_ptr->child_to_current_deviation != MIN_SIGNED_VALUE) {
                                    if (context_ptr->parent_to_current_deviation <= context_ptr->child_to_current_deviation) {
                                        add_sub_depth = 0;
                                    }
                                    else {
                                        add_parent_depth = 0;
                                    }
                                }
                            }
                        }

                        if (s_depth != 0 && add_parent_depth)
                            set_parent_to_be_considered(results_ptr, blk_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth, sq_size_idx,
                                pcs_ptr->parent_pcs_ptr->disallow_nsq, s_depth);

                        if (e_depth != 0 && add_sub_depth)
                            set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, blk_index, sb_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth, sq_size_idx, e_depth);
                    }
#else
                    int8_t s_depth = context_ptr->depth_ctrls.s_depth;
                    int8_t e_depth = context_ptr->depth_ctrls.e_depth;
                    // Check that the start and end depth are in allowed range, given other features
                    // which restrict allowable depths
                    if (context_ptr->disallow_4x4) {
                        e_depth = (blk_geom->sq_size == 8) ? 0
                                : (blk_geom->sq_size == 16) ? MIN(1, e_depth)
                                : (blk_geom->sq_size == 32) ? MIN(2, e_depth)
                                : e_depth;
                    }
                    if (context_ptr->depth_removal_ctrls.enabled) {
                        if (context_ptr->depth_removal_ctrls.disallow_below_64x64) {
                            e_depth = (blk_geom->sq_size <= 64) ? 0
                                : (blk_geom->sq_size == 128) ? MIN(1, e_depth) : e_depth;
                        }
                        else if (context_ptr->depth_removal_ctrls.disallow_below_32x32) {
                            e_depth = (blk_geom->sq_size <= 32) ? 0
                                : (blk_geom->sq_size == 64) ? MIN(1, e_depth)
                                : (blk_geom->sq_size == 128) ? MIN(2, e_depth) : e_depth;
                        }
                        else if (context_ptr->depth_removal_ctrls.disallow_below_16x16) {
                        e_depth = (blk_geom->sq_size <= 16) ? 0
                                : (blk_geom->sq_size ==  32) ? MIN(1, e_depth)
                                : (blk_geom->sq_size ==  64) ? MIN(2, e_depth)
                                : (blk_geom->sq_size == 128) ? MIN(3, e_depth) : e_depth;
                    }
                }
                    // Add current pred depth block(s)
                    for (unsigned block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
                        results_ptr->consider_block[blk_index + block_1d_idx] = 1;
                        results_ptr->refined_split_flag[blk_index + block_1d_idx] = EB_FALSE;
                    }

                    uint8_t sq_size_idx = 7 - (uint8_t)svt_log2f((uint8_t)blk_geom->sq_size);
#if OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS
                    int64_t th_offset = 0;
                    if (context_ptr->depth_refinement_ctrls.enabled && context_ptr->depth_refinement_ctrls.cost_band_based_modulation && (s_depth != 0 || e_depth != 0)) {
                        update_pred_th_offset(context_ptr, blk_geom, &s_depth, &e_depth, &th_offset);
                    }
#else
                    // Update pred and generate an offset to be used @ sub_to_current_th and parent_to_current_th derivation based on the cost range of the predicted block; use default ths for high cost(s) and more aggressive TH(s) or Pred only for low cost(s)
                    int64_t th_offset = 0;
                    if (context_ptr->depth_refinement_ctrls.enabled && context_ptr->depth_refinement_ctrls.use_pred_block_cost && (s_depth != 0 || e_depth != 0)) {
                        update_pred_th_offset(context_ptr, blk_geom, &s_depth, &e_depth, &th_offset);
                    }
#endif
                    // Add block indices of upper depth(s)
                    // Block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
                    uint8_t add_parent_depth = 1;
                    if (context_ptr->depth_refinement_ctrls.enabled && s_depth == -1 && pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[blk_index] && blk_geom->sq_size < ((scs_ptr->seq_header.sb_size == BLOCK_128X128) ? 128 : 64)) {
                        add_parent_depth = is_parent_to_current_deviation_small(
                            scs_ptr, context_ptr, blk_geom, th_offset);
                    }
                    // Add block indices of lower depth(s)
                    // Block-based depth refinement using cost is applicable for only [s_depth=-1, e_depth=1]
                    uint8_t add_sub_depth = 1;
                    if (context_ptr->depth_refinement_ctrls.enabled && e_depth == 1 && pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_allowed[blk_index] && blk_geom->sq_size > 4) {
                        add_sub_depth = is_child_to_current_deviation_small(
                            scs_ptr, context_ptr, blk_geom, blk_index, th_offset);
                    }
                    // Use a maximum of 2 depth per block (PRED+Parent or PRED+Sub)
                    if (context_ptr->depth_refinement_ctrls.enabled && context_ptr->depth_refinement_ctrls.up_to_2_depth) {
                        if ((s_depth == -1) && add_parent_depth && (e_depth == 1) && add_sub_depth) {
                            if (context_ptr->parent_to_current_deviation != MIN_SIGNED_VALUE && context_ptr->child_to_current_deviation != MIN_SIGNED_VALUE) {
                                if (context_ptr->parent_to_current_deviation <= context_ptr->child_to_current_deviation) {
                                    add_sub_depth = 0;
                                }
                                else {
                                    add_parent_depth = 0;
                                }
                            }
                        }
                    }
                    if (add_parent_depth)
                    if (s_depth != 0)
                        set_parent_to_be_considered(
                            results_ptr, blk_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth,sq_size_idx,
                            pcs_ptr->parent_pcs_ptr->disallow_nsq, s_depth);
                    if (add_sub_depth)
                    if (e_depth != 0)
                        set_child_to_be_considered(pcs_ptr, context_ptr, results_ptr, blk_index, sb_index, scs_ptr->seq_header.sb_size, (int8_t)blk_geom->depth,sq_size_idx, e_depth);
#endif
                }
            }
        }
#if LIGHT_PD0
        blk_index +=
            split_flag
            ? blk_geom->d1_depth_offset
            : blk_geom->ns_depth_offset;
#else
        blk_index +=
            split_flag
                ? d1_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[scs_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
#endif
    }
}
// Initialize structures used to indicate which blocks will be tested at MD.
// MD data structures should be updated in init_block_data(), not here.
#if OPT_MEM_PALETTE
 EbErrorType  build_starting_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr, uint32_t sb_index) {
#else
static void build_starting_cand_block_array(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr, uint32_t sb_index) {
#endif
    memset(context_ptr->tested_blk_flag, 0, sizeof(uint8_t) * scs_ptr->max_block_cnt);
#if !CLN_REMOVE_UNUSED_FEATS
    memset(context_ptr->do_not_process_blk, 0, sizeof(uint8_t) * scs_ptr->max_block_cnt);
#endif
#if OPT_PD0_PATH
    memset(context_ptr->avail_blk_flag, EB_FALSE, sizeof(uint8_t) * scs_ptr->max_block_cnt);
#endif
    MdcSbData *results_ptr = context_ptr->mdc_sb_array;
    results_ptr->leaf_count = 0;
    uint32_t blk_index = 0;
#if !LIGHT_PD0
    const BlockSize sb_size = scs_ptr->seq_header.sb_size;
#endif
    const uint16_t max_block_cnt = scs_ptr->max_block_cnt;
    const int32_t min_sq_size =
        (context_ptr->depth_removal_ctrls.enabled && context_ptr->depth_removal_ctrls.disallow_below_64x64)
        ? 64
        : (context_ptr->depth_removal_ctrls.enabled && context_ptr->depth_removal_ctrls.disallow_below_32x32)
        ? 32
        : (context_ptr->depth_removal_ctrls.enabled && context_ptr->depth_removal_ctrls.disallow_below_16x16)
        ? 16
        : context_ptr->disallow_4x4 ? 8 : 4;

    // Loop over all blocks to initialize data for partitions to be tested
    while (blk_index < max_block_cnt) {
        const BlockGeom *blk_geom = get_blk_geom_mds(blk_index);
        // SQ/NSQ block(s) filter based on the SQ size
#if FTR_M13
        uint8_t is_block_tagged = (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE) || (blk_geom->sq_size < min_sq_size) ? 0 : 1;
        if (context_ptr->skip_pd0)
            is_block_tagged = !context_ptr->depth_removal_ctrls.disallow_below_64x64 && context_ptr->depth_removal_ctrls.disallow_below_32x32 && (blk_geom->sq_size != 32) ? 0 : is_block_tagged;
#else
        const uint8_t is_block_tagged =
            (blk_geom->sq_size == 128 && pcs_ptr->slice_type == I_SLICE) ||
            (blk_geom->sq_size < min_sq_size)
            ? 0
            : 1;
#endif

        // SQ/NSQ block(s) filter based on the block validity
        if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_block_tagged) {

            const uint32_t tot_d1_blocks = pcs_ptr->parent_pcs_ptr->disallow_nsq ? 1 :
                get_tot_1d_blks(pcs_ptr->parent_pcs_ptr, blk_geom->sq_size, context_ptr->md_disallow_nsq);

            for (uint32_t idx = blk_index; idx < (tot_d1_blocks + blk_index); ++idx) {

                if (pcs_ptr->parent_pcs_ptr->sb_geom[sb_index].block_is_inside_md_scan[idx]) {
#if OPT_MEM_PALETTE
                    //  MD palette info buffer
                    if (pcs_ptr->parent_pcs_ptr->palette_level) {
                        if (context_ptr->md_blk_arr_nsq[idx].palette_mem == 0) {
                            rtime_alloc_palette_info(&context_ptr->md_blk_arr_nsq[idx]);
                            context_ptr->md_blk_arr_nsq[idx].palette_mem = 1;
                        }
                    }

                    context_ptr->md_blk_arr_nsq[idx].palette_size[0] = 0;
                    context_ptr->md_blk_arr_nsq[idx].palette_size[1] = 0;
#endif
                    results_ptr->leaf_data_array[results_ptr->leaf_count].mds_idx = idx;
                    results_ptr->leaf_data_array[results_ptr->leaf_count].tot_d1_blocks = tot_d1_blocks;
                    results_ptr->split_flag[results_ptr->leaf_count++] = (blk_geom->sq_size > min_sq_size) ? EB_TRUE : EB_FALSE;
                }
            }
#if LIGHT_PD0
            blk_index += blk_geom->d1_depth_offset;
#else
            blk_index += d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
#endif
        }
        else {
#if LIGHT_PD0
#if FTR_M13
            if (context_ptr->skip_pd0)
                context_ptr->md_blk_arr_nsq[blk_index].part = (blk_geom->sq_size > min_sq_size) ? PARTITION_SPLIT : PARTITION_NONE;
#endif
            blk_index +=
                (blk_geom->sq_size > min_sq_size)
                ? blk_geom->d1_depth_offset
                : blk_geom->ns_depth_offset;
#else
            blk_index +=
                (blk_geom->sq_size > min_sq_size)
                ? d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth]
                : ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
#endif
        }
    }

#if OPT_MEM_PALETTE
    return EB_ErrorNone;
#endif
}
void recode_loop_update_q(
    PictureParentControlSet *ppcs_ptr,
    int *const loop, int *const q, int *const q_low,
    int *const q_high, const int top_index, const int bottom_index,
    int *const undershoot_seen, int *const overshoot_seen,
    int *const low_cr_seen, const int loop_count);
void sb_qp_derivation_tpl_la(PictureControlSet *pcs_ptr);
void mode_decision_configuration_init_qp_update(PictureControlSet *pcs_ptr);
void init_enc_dec_segement(PictureParentControlSet *parentpicture_control_set_ptr);

static void recode_loop_decision_maker(PictureControlSet *pcs_ptr,
            SequenceControlSet *scs_ptr, EbBool *do_recode) {
    PictureParentControlSet *ppcs_ptr = pcs_ptr->parent_pcs_ptr;
    EncodeContext *const encode_context_ptr = ppcs_ptr->scs_ptr->encode_context_ptr;
    RATE_CONTROL *const rc = &(encode_context_ptr->rc);
    int32_t loop = 0;
    FrameHeader *frm_hdr = &ppcs_ptr->frm_hdr;
    int32_t q = frm_hdr->quantization_params.base_q_idx;
    if (ppcs_ptr->loop_count == 0) {
        ppcs_ptr->q_low = ppcs_ptr->bottom_index;
        ppcs_ptr->q_high = ppcs_ptr->top_index;
    }

    // Update q and decide whether to do a recode loop
    recode_loop_update_q(ppcs_ptr, &loop, &q,
            &ppcs_ptr->q_low, &ppcs_ptr->q_high,
            ppcs_ptr->top_index, ppcs_ptr->bottom_index,
            &ppcs_ptr->undershoot_seen, &ppcs_ptr->overshoot_seen,
            &ppcs_ptr->low_cr_seen, ppcs_ptr->loop_count);

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

        ppcs_ptr->picture_qp =
            (uint8_t)CLIP3((int32_t)scs_ptr->static_config.min_qp_allowed,
                    (int32_t)scs_ptr->static_config.max_qp_allowed,
                    (frm_hdr->quantization_params.base_q_idx + 2) >> 2);
        pcs_ptr->picture_qp = ppcs_ptr->picture_qp;

        // 2pass QPM with tpl_la
        if (scs_ptr->static_config.enable_adaptive_quantization == 2 &&
#if !FTR_RC_CAP
            !use_output_stat(scs_ptr) &&
            (use_input_stat(scs_ptr) || scs_ptr->lap_enabled) &&
#endif
            scs_ptr->static_config.enable_tpl_la &&
            ppcs_ptr->r0 != 0)
            sb_qp_derivation_tpl_la(pcs_ptr);
        else
        {
            ppcs_ptr->frm_hdr.delta_q_params.delta_q_present = 0;
#if !RFCTR_RC_P2
            ppcs_ptr->average_qp = 0;
#endif
            for (int sb_addr = 0; sb_addr < pcs_ptr->sb_total_count_pix; ++sb_addr) {
                SuperBlock * sb_ptr = pcs_ptr->sb_ptr_array[sb_addr];
                sb_ptr->qindex   = quantizer_to_qindex[pcs_ptr->picture_qp];
#if !RFCTR_RC_P2
                ppcs_ptr->average_qp += pcs_ptr->picture_qp;
#endif
            }
        }
    } else {
        ppcs_ptr->loop_count = 0;
    }
}
#if !OPT_PD0_PATH
static void init_avail_blk_flag(SequenceControlSet *scs_ptr, ModeDecisionContext *context_ptr) {
    // Initialize avail_blk_flag to false
    memset(context_ptr->avail_blk_flag, EB_FALSE, sizeof(uint8_t) * scs_ptr->max_block_cnt);
}
#endif

#if LIGHT_PD1_MACRO
/* for debug/documentation purposes: list all features assumed off for light pd1*/
void exaustive_light_pd1_features(
    ModeDecisionContext       *md_ctx,
    PictureParentControlSet   *ppcs,
    uint8_t                    use_light_pd1,
    uint8_t                    debug_lpd1_features
)
{

    if (debug_lpd1_features) {
        uint8_t light_pd1;

        // Use light-PD1 path if the assumed features are off
        if (
            md_ctx->obmc_ctrls.enabled == 0 &&
            md_ctx->md_allow_intrabc == 0 &&
            md_ctx->hbd_mode_decision == 0 &&
            md_ctx->ifs_ctrls.interpolation_search_level == IFS_OFF &&
            ppcs->frm_hdr.allow_warped_motion == 0 &&
            md_ctx->inter_intra_comp_ctrls.enabled == 0 &&
#if !CLN_REG_PD1_TX_CTRLS
            md_ctx->reduce_last_md_stage_candidate == 4 &&
#endif
#if CLN_RATE_EST_CTRLS
            md_ctx->rate_est_ctrls.update_skip_ctx_dc_sign_ctx == 0 &&
#else
            md_ctx->shut_skip_ctx_dc_sign_update == 1 &&
#endif
            md_ctx->spatial_sse_ctrls.spatial_sse_full_loop_level == 0 &&
            md_ctx->md_sq_me_ctrls.enabled == 0 &&
#if !OPT_LPD1_PME
            md_ctx->md_pme_ctrls.enabled == 0 &&
#endif
#if OPT_REMOVE_TXT_LPD1
            md_ctx->txt_ctrls.enabled == 0 &&
#endif
#if CLN_MDS0_CTRLS
            md_ctx->mds0_ctrls.mds0_dist_type != MDS0_SSD &&
#else
            md_ctx->mds0_dist_type != MDS0_SSD &&
#endif
            md_ctx->unipred3x3_injection == 0 &&
            md_ctx->bipred3x3_injection == 0 &&
            md_ctx->inter_compound_mode == 0 &&
            md_ctx->md_pic_obmc_level == 0 &&
            md_ctx->md_filter_intra_level == 0 &&
            md_ctx->new_nearest_near_comb_injection == 0 &&
            md_ctx->md_palette_level == 0 &&
#if CLN_CAND_REDUCTION_CTRLS
            md_ctx->cand_reduction_ctrls.merge_inter_classes&&
#else
            md_ctx->merge_inter_classes &&
#endif
            ppcs->gm_ctrls.enabled == 0 &&
#if OPT_REDUCE_COPIES_LPD1
            // If TXS enabled at picture level, there are necessary context updates that must be added to LPD1
            ppcs->frm_hdr.tx_mode != TX_MODE_SELECT &&
#endif
            md_ctx->txs_ctrls.enabled == 0 &&
            md_ctx->pred_depth_only &&
            md_ctx->md_disallow_nsq == EB_TRUE &&
            md_ctx->disallow_4x4 == EB_TRUE &&
            ppcs->scs_ptr->static_config.super_block_size == 64 &&
            ppcs->ref_list0_count_try == 1 &&
            ppcs->ref_list1_count_try == 1 &&
#if SS_CLN_CFL_CTRLS
            md_ctx->cfl_ctrls.enabled == 0 &&
#else
            md_ctx->md_disable_cfl == EB_TRUE &&
#endif
            md_ctx->uv_ctrls.nd_uv_serach_mode == 0 &&
            md_ctx->uv_ctrls.uv_mode == CHROMA_MODE_1
            ) {
            light_pd1 = 1;
        }
        else {
            light_pd1 = 0;
        }

        assert_err( light_pd1 == use_light_pd1, "Warning: light PD1 feature assumption is broken \n");

    }
}
#endif
#if CLN_LPD1_LVLS
/* Light-PD1 classifier used when cost/coeff info is available.  If PD0 is skipped, or the trasnsform is
not performed, a separate detector (lpd1_detector_skip_pd0) is used. */
void lpd1_detector_post_pd0(PictureControlSet* pcs, ModeDecisionContext* md_ctx) {

    for (int pd1_lvl = LPD1_LEVELS - 1; pd1_lvl > REGULAR_PD1; pd1_lvl--) {
        if (md_ctx->lpd1_ctrls.pd1_level == pd1_lvl) {
            if (md_ctx->lpd1_ctrls.use_lpd1_detector[pd1_lvl]) {

                // Use info from ref frames (if available)
                if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl] && pcs->slice_type != I_SLICE) {

                    EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                    uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
                    if (pcs->slice_type == B_SLICE) {
                        EbReferenceObject *ref_obj_l1 = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                        l1_was_intra = ref_obj_l1->sb_intra[md_ctx->sb_index];
                    }
                    if (l0_was_intra && l1_was_intra) {
                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        continue;
                    }
                    else if (l0_was_intra || l1_was_intra) {
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
                const uint32_t nz_coeffs = md_ctx->avail_blk_flag[0] == EB_TRUE ? md_ctx->md_local_blk_unit[0].count_non_zero_coeffs : (uint32_t)~0;

                const uint32_t lambda = md_ctx->full_sb_lambda_md[EB_8_BIT_MD]; // light-PD1 assumes 8-bit MD
                const uint32_t rate = md_ctx->lpd1_ctrls.coeff_th[pd1_lvl];
                const uint32_t dist = md_ctx->lpd1_ctrls.cost_th_dist[pd1_lvl];
                /* dist << 14 is equivalent to 64 * 64 * 4 * dist (64 * 64 so the distortion is the per-pixel SSD) and 4 because
                the distortion of the 64x64 block is shifted by 2 (same as multiplying by 4) in perform_tx_light_pd0. */
                const uint64_t low_th = RDCOST(lambda, 6000 + rate * 500, (uint64_t)dist << 14);
                const uint16_t coeff_th = md_ctx->lpd1_ctrls.coeff_th[pd1_lvl];

                // If the PD0 cost is very high and the number of non-zero coeffs is high, the block is difficult, so should use regular PD1
                if (pd0_cost > low_th && nz_coeffs > coeff_th) {
                    md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                }

                // If the best PD0 mode was INTER, check the MV length
                if (md_ctx->md_blk_arr_nsq[0].prediction_mode_flag == INTER_MODE && md_ctx->lpd1_ctrls.max_mv_length[pd1_lvl] != (uint16_t)~0 &&
                    md_ctx->avail_blk_flag[0] == EB_TRUE) {


                    PredictionUnit* pu_ptr = md_ctx->md_blk_arr_nsq[0].prediction_unit_array;
                    const uint16_t max_mv_length = md_ctx->lpd1_ctrls.max_mv_length[pd1_lvl];

                    if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0)
                    {
                        if (pu_ptr->mv[REF_LIST_0].x > max_mv_length || pu_ptr->mv[REF_LIST_0].y > max_mv_length)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                    }
                    else if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1)
                    {
                        if (pu_ptr->mv[REF_LIST_1].x > max_mv_length || pu_ptr->mv[REF_LIST_1].y > max_mv_length)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                    }
                    else //if (pu_ptr->inter_pred_direction_index == BI_PRED)
                    {
                        assert(pu_ptr->inter_pred_direction_index == BI_PRED);
                        if (pu_ptr->mv[REF_LIST_0].x > max_mv_length || pu_ptr->mv[REF_LIST_0].y > max_mv_length)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        if (pu_ptr->mv[REF_LIST_1].x > max_mv_length || pu_ptr->mv[REF_LIST_1].y > max_mv_length)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                    }
                }

                if (pcs->slice_type != I_SLICE) {
                    /* me_8x8_cost_variance_th is shifted by 5 then mulitplied by the pic QP (max 63).  Therefore, the TH must be less than
                       (((uint32_t)~0) >> 1) to avoid overflow issues from the multiplication. */
                    if (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] < (((uint32_t)~0) >> 1) &&
                        pcs->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] > (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] >> 5) * pcs->picture_qp)
                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                }

            }
        }
    }
}

/* Light-PD1 classifier used when cost/coeff info is unavailable.  If PD0 is skipped, or the trasnsform is
not performed, this detector is used (else lpd1_detector_post_pd0() is used). */
void lpd1_detector_skip_pd0(PictureControlSet* pcs, ModeDecisionContext* md_ctx, uint32_t pic_width_in_sb) {

    const int16_t left_sb_index = (int16_t)md_ctx->sb_index - 1;
    const int16_t top_sb_index = (int16_t)md_ctx->sb_index - (int16_t)pic_width_in_sb;

    for (int pd1_lvl = LPD1_LEVELS - 1; pd1_lvl > REGULAR_PD1; pd1_lvl--) {
        if (md_ctx->lpd1_ctrls.pd1_level == pd1_lvl) {
            if (md_ctx->lpd1_ctrls.use_lpd1_detector[pd1_lvl]) {

                // Use info from ref. frames (if available)
                if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl] && pcs->slice_type != I_SLICE) {

                    EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                    uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
                    uint8_t l0_was_skip = ref_obj_l0->sb_skip[md_ctx->sb_index], l1_was_skip = 1;
                    if (pcs->slice_type == B_SLICE) {
                        EbReferenceObject *ref_obj_l1 = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                        l1_was_intra = ref_obj_l1->sb_intra[md_ctx->sb_index];
                        l1_was_skip = ref_obj_l1->sb_skip[md_ctx->sb_index];
                    }

                    if (l0_was_intra && l1_was_intra) {
                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        continue;
                    }
                    else if (!l0_was_skip && !l1_was_skip && (l0_was_intra || l1_was_intra)) {
                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        continue;
                    }
                }

                // I_SLICE doesn't have ME info
                if (pcs->slice_type != I_SLICE) {

                    // If the SB origin of one dimension is zero, then this SB is the first block in a row/column, so won't have neighbours
                    if (md_ctx->sb_origin_x == 0 || md_ctx->sb_origin_y == 0) {
                        if (pcs->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] > md_ctx->lpd1_ctrls.skip_pd0_edge_dist_th[pd1_lvl])
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        /* me_8x8_cost_variance_th is shifted by 5 then mulitplied by the pic QP (max 63).  Therefore, the TH must be less than
                           (((uint32_t)~0) >> 1) to avoid overflow issues from the multiplication. */
                        else if (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] < (((uint32_t)~0) >> 1) &&
                            pcs->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] > (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] >> 5) * pcs->picture_qp)
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                    }
                    else {
                        if (pcs->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] >
                            ((pcs->parent_pcs_ptr->me_64x64_distortion[left_sb_index] + pcs->parent_pcs_ptr->me_64x64_distortion[top_sb_index]) << md_ctx->lpd1_ctrls.skip_pd0_me_shift[pd1_lvl]))
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        else if (pcs->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] >
                            ((pcs->parent_pcs_ptr->me_8x8_cost_variance[left_sb_index] + pcs->parent_pcs_ptr->me_8x8_cost_variance[top_sb_index]) << md_ctx->lpd1_ctrls.skip_pd0_me_shift[pd1_lvl])) {
                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                        }
                        else if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl]) {
                            // Use info from neighbouring SBs
                            if (pcs->sb_intra[left_sb_index] && pcs->sb_intra[top_sb_index]) {
                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                            }
                            else if (!pcs->sb_skip[left_sb_index] && !pcs->sb_intra[top_sb_index] &&
                                (pcs->sb_intra[left_sb_index] || pcs->sb_intra[top_sb_index])) {
                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif

#if CLN_MERGE_LPD0_VLPD0
/*
* Check whether vlpd0 is safe or not
*/
uint8_t is_vlpd0_safe(PictureControlSet* pcs_ptr, ModeDecisionContext* md_ctx) {

    uint8_t is_vlpd0_safe = EB_TRUE;

    EbReferenceObject* ref_obj_l0 = (EbReferenceObject*)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
    uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
    if (pcs_ptr->slice_type == B_SLICE) {
        EbReferenceObject* ref_obj_l1 = (EbReferenceObject*)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
        l1_was_intra = ref_obj_l1->sb_intra[md_ctx->sb_index];
    }
    if (l0_was_intra || l1_was_intra) {
        return EB_FALSE;
    }

    uint32_t me_8x8_cost_variance_th = 250000;
    if (pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] > (me_8x8_cost_variance_th >> 5) * pcs_ptr->picture_qp)
        return EB_FALSE;

    return is_vlpd0_safe;

}
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
#if LIGHT_PD1_MACRO || FIX_DO_NOT_TEST_CORRUPTED_MVS
        ModeDecisionContext *md_ctx = context_ptr->md_context;
#endif
#if LIGHT_PD1_MACRO
        struct PictureParentControlSet *ppcs = pcs_ptr->parent_pcs_ptr;
#endif
#if FTR_10BIT_MDS3_REG_PD1
        md_ctx->encoder_bit_depth = (uint8_t)scs_ptr->static_config.encoder_bit_depth;
#endif
#if FIX_DO_NOT_TEST_CORRUPTED_MVS
        md_ctx->corrupted_mv_check = (pcs_ptr->parent_pcs_ptr->aligned_width >= (1 << (MV_IN_USE_BITS - 3))) || (pcs_ptr->parent_pcs_ptr->aligned_height >= (1 << (MV_IN_USE_BITS - 3)));
#endif
        context_ptr->tile_group_index = enc_dec_tasks_ptr->tile_group_index;
        context_ptr->coded_sb_count   = 0;
        segments_ptr = pcs_ptr->enc_dec_segment_ctrl[context_ptr->tile_group_index];
        // SB Constants
        uint8_t sb_sz      = (uint8_t)scs_ptr->sb_size_pix;
        uint8_t sb_size_log2 = (uint8_t)svt_log2f(sb_sz);
        context_ptr->sb_sz = sb_sz;
        uint32_t pic_width_in_sb = (pcs_ptr->parent_pcs_ptr->aligned_width + sb_sz - 1) >>
            sb_size_log2;
        uint16_t tile_group_width_in_sb = pcs_ptr->parent_pcs_ptr
                                              ->tile_group_info[context_ptr->tile_group_index]
                                              .tile_group_width_in_sb;
#if  FTR_INTRA_DETECTOR
        context_ptr->tot_intra_coded_area       = 0;
#endif
#if FTR_COEFF_DETECTOR
        context_ptr->tot_skip_coded_area = 0;
#endif
        // Bypass encdec for the first pass
#if  FTR_OPT_MPASS_BYPASS_FRAMES
#if FTR_OP_TEST
        if (use_output_stat(scs_ptr) || (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag && 1 && !pcs_ptr->parent_pcs_ptr->first_frame_in_minigop)) {
#else
        if (use_output_stat(scs_ptr) || (!pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag && scs_ptr->rc_stat_gen_pass_mode  && !pcs_ptr->parent_pcs_ptr->first_frame_in_minigop)) {
#endif
#else
        if (use_output_stat(scs_ptr)) {
#endif
            svt_release_object(pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr);
            pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr = (EbObjectWrapper *)NULL;
            pcs_ptr->parent_pcs_ptr->pa_me_data = NULL;
            // Get Empty EncDec Results
            svt_get_empty_object(context_ptr->enc_dec_output_fifo_ptr, &enc_dec_results_wrapper_ptr);
            enc_dec_results_ptr = (EncDecResults *)enc_dec_results_wrapper_ptr->object_ptr;
            enc_dec_results_ptr->pcs_wrapper_ptr = enc_dec_tasks_ptr->pcs_wrapper_ptr;
            enc_dec_results_ptr->completed_sb_row_index_start = 0;
            enc_dec_results_ptr->completed_sb_row_count =
                ((pcs_ptr->parent_pcs_ptr->aligned_height + scs_ptr->sb_size_pix - 1) >> sb_size_log2);
            // Post EncDec Results
            svt_post_full_object(enc_dec_results_wrapper_ptr);
        }
        else{
            if (enc_dec_tasks_ptr->input_type == ENCDEC_TASKS_SUPERRES_INPUT) {
                // do as dorecode do
                pcs_ptr->enc_dec_coded_sb_count = 0;
                // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
                if (context_ptr->is_md_rate_estimation_ptr_owner) {
                    EB_FREE_ARRAY(context_ptr->md_rate_estimation_ptr);
                    context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
                }
                context_ptr->md_rate_estimation_ptr = pcs_ptr->md_rate_estimation_array;
                // re-init mode decision configuration for qp update for re-encode frame
                mode_decision_configuration_init_qp_update(pcs_ptr);
                // init segment for re-encode frame
                init_enc_dec_segement(pcs_ptr->parent_pcs_ptr);

                // post tile based encdec task
                EbObjectWrapper* enc_dec_re_encode_tasks_wrapper_ptr;
                uint16_t tg_count =
                    pcs_ptr->parent_pcs_ptr->tile_group_cols * pcs_ptr->parent_pcs_ptr->tile_group_rows;
                for (uint16_t tile_group_idx = 0; tile_group_idx < tg_count; tile_group_idx++) {
                    svt_get_empty_object(context_ptr->enc_dec_feedback_fifo_ptr,
                        &enc_dec_re_encode_tasks_wrapper_ptr);

                    EncDecTasks* enc_dec_re_encode_tasks_ptr = (EncDecTasks*)enc_dec_re_encode_tasks_wrapper_ptr->object_ptr;
                    enc_dec_re_encode_tasks_ptr->pcs_wrapper_ptr = enc_dec_tasks_ptr->pcs_wrapper_ptr;
                    enc_dec_re_encode_tasks_ptr->input_type = ENCDEC_TASKS_MDC_INPUT;
                    enc_dec_re_encode_tasks_ptr->tile_group_index = tile_group_idx;

                    // Post the Full Results Object
                    svt_post_full_object(enc_dec_re_encode_tasks_wrapper_ptr);
                }

                svt_release_object(enc_dec_tasks_wrapper_ptr);
                continue;
            }

#if CLN_UPDATE_CDF
        if (pcs_ptr->cdf_ctrl.enabled) {
#endif
            if (!pcs_ptr->cdf_ctrl.update_mv)
                copy_mv_rate(pcs_ptr, &context_ptr->md_context->rate_est_table);
            if (!pcs_ptr->cdf_ctrl.update_se)

#if OPT8_MDC
                av1_estimate_syntax_rate(&context_ptr->md_context->rate_est_table,
                    pcs_ptr->slice_type == I_SLICE ? EB_TRUE : EB_FALSE,
                    pcs_ptr->pic_filter_intra_level,
                    pcs_ptr->parent_pcs_ptr->frm_hdr.allow_screen_content_tools,
                    scs_ptr->seq_header.enable_restoration,
                    pcs_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc,
#if OPT9_RATE_ESTIMATION
                    pcs_ptr->parent_pcs_ptr->partition_contexts,
#endif
                    &pcs_ptr->md_frame_context);
#else
                av1_estimate_syntax_rate(&context_ptr->md_context->rate_est_table,
                    pcs_ptr->slice_type == I_SLICE ? EB_TRUE : EB_FALSE,
                    &pcs_ptr->md_frame_context);
#endif
            if (!pcs_ptr->cdf_ctrl.update_coef)
                av1_estimate_coefficients_rate(&context_ptr->md_context->rate_est_table,
                    &pcs_ptr->md_frame_context);
#if CLN_UPDATE_CDF
        }
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
#if !FIX_SANITIZER_RACE_CONDS
            if (pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL)
                ((EbReferenceObject *)
                     pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                    ->average_intensity = pcs_ptr->parent_pcs_ptr->average_intensity[0];
#endif
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
                    sb_index = context_ptr->md_context->sb_index =(uint16_t)((y_sb_index + tile_group_y_sb_start) * pic_width_in_sb +
                        x_sb_index + tile_group_x_sb_start);
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
                    mdc_ptr = context_ptr->md_context->mdc_sb_array;
                    context_ptr->sb_index = sb_index;
                    if (pcs_ptr->cdf_ctrl.enabled) {
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
                        if (pcs_ptr->cdf_ctrl.update_se)
#if OPT8_MDC
                            av1_estimate_syntax_rate(&context_ptr->md_context->rate_est_table,
                                pcs_ptr->slice_type == I_SLICE,
                                pcs_ptr->pic_filter_intra_level,
                                pcs_ptr->parent_pcs_ptr->frm_hdr.allow_screen_content_tools,
                                scs_ptr->seq_header.enable_restoration,
                                pcs_ptr->parent_pcs_ptr->frm_hdr.allow_intrabc,
#if OPT9_RATE_ESTIMATION
                                pcs_ptr->parent_pcs_ptr->partition_contexts,
#endif
                                &pcs_ptr->ec_ctx_array[sb_index]);
#else
                        av1_estimate_syntax_rate(&context_ptr->md_context->rate_est_table,
                            pcs_ptr->slice_type == I_SLICE,
                            &pcs_ptr->ec_ctx_array[sb_index]);
#endif
                        // Initial Rate Estimation of the Motion vectors
                        if (pcs_ptr->cdf_ctrl.update_mv)
                        av1_estimate_mv_rate(pcs_ptr,
                            &context_ptr->md_context->rate_est_table,
                            &pcs_ptr->ec_ctx_array[sb_index]);

                        if (pcs_ptr->cdf_ctrl.update_coef)
                        av1_estimate_coefficients_rate(&context_ptr->md_context->rate_est_table,
                            &pcs_ptr->ec_ctx_array[sb_index]);
                        context_ptr->md_context->md_rate_estimation_ptr =
                            &context_ptr->md_context->rate_est_table;
                    }
                    // Configure the SB
                    mode_decision_configure_sb(
                        context_ptr->md_context, pcs_ptr, (uint8_t)sb_ptr->qindex);
                    // signals set once per SB (i.e. not per PD)
                    signal_derivation_enc_dec_kernel_common(scs_ptr, pcs_ptr, context_ptr->md_context);


#if OPT_MEM_PALETTE
                if(pcs_ptr->parent_pcs_ptr->palette_level)
                    // Status of palette info alloc
                    for (int i = 0;i< scs_ptr->max_block_cnt; ++i)
                        context_ptr->md_context->md_blk_arr_nsq[i].palette_mem   = 0;
#endif

#if FIX_REMOVE_PD1

#if FTR_SUBSAMPLE_RESIDUAL
                    // Initialize is_subres_safe
                    context_ptr->md_context->is_subres_safe = (uint8_t) ~0;
#endif
#if FTR_10BIT_MDS3_REG_PD1
                    // Signal initialized here; if needed, will be set in md_encode_block before MDS3
                    md_ctx->need_hbd_comp_mds3 = 0;
#endif
                    uint8_t skip_pd_pass_0 = (scs_ptr->static_config.super_block_size == 64 && context_ptr->md_context->depth_removal_ctrls.disallow_below_64x64)
                        ? 1
                        : 0;
#if FTR_M13
                    if (context_ptr->md_context->skip_pd0)
                        if (context_ptr->md_context->depth_removal_ctrls.disallow_below_32x32)
                            skip_pd_pass_0 = 1;
#endif
#if FTR_VLPD0
                    if (context_ptr->md_context->pd0_level == VERY_LIGHT_PD0) {

#if CLN_MERGE_LPD0_VLPD0
                        // Use the next conservative level if not safe to use VLPD0
                        if(!is_vlpd0_safe(pcs_ptr,md_ctx))
                            context_ptr->md_context->pd0_level = VERY_LIGHT_PD0 - 1;
#else
                        EbReferenceObject* ref_obj_l0 = (EbReferenceObject*)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                        uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
                        if (pcs_ptr->slice_type == B_SLICE) {
                            EbReferenceObject* ref_obj_l1 = (EbReferenceObject*)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                            l1_was_intra = ref_obj_l1->sb_intra[md_ctx->sb_index];
                        }
                        if (l0_was_intra || l1_was_intra) {
                            context_ptr->md_context->pd0_level = VERY_LIGHT_PD0 - 1;
                        }

                        uint32_t me_8x8_cost_variance_th = 250000;
                        if (pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[sb_index] > (me_8x8_cost_variance_th >> 5) * pcs_ptr->picture_qp)
                            context_ptr->md_context->pd0_level = VERY_LIGHT_PD0 - 1;
#endif
                    }
#endif
#if FTR_LPD1_DETECTOR
                    // PD0 is only skipped if there is a single depth to test
                    if (skip_pd_pass_0)
                        md_ctx->pred_depth_only = 1;
#endif
                    // Multi-Pass PD
                    if (!skip_pd_pass_0 && pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_ON) {
#else
                    uint8_t pd_pass_2_only = (scs_ptr->static_config.super_block_size == 64 && context_ptr->md_context->depth_removal_ctrls.disallow_below_64x64)
                        ? 1
                        : 0;

                    // Multi-Pass PD
                    if (!pd_pass_2_only &&
                        (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_0 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_1 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4)
                        ) {
#endif
#if OPT_PD0_PATH
                        // [PD_PASS_0]
                        // Input : mdc_blk_ptr built @ mdc process (up to 4421)
                        // Output: md_blk_arr_nsq reduced set of block(s)
                        context_ptr->md_context->pd_pass = PD_PASS_0;
#else
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
                        signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);
                        // [PD_PASS_0]
                        // Input : mdc_blk_ptr built @ mdc process (up to 4421)
                        // Output: md_blk_arr_nsq reduced set of block(s)
                        // Build the t=0 cand_block_array
                        build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
                        // Initialize avail_blk_flag to false
                        init_avail_blk_flag(scs_ptr, context_ptr->md_context);
#endif
#if LIGHT_PD0
                        // skip_intra much be TRUE for non-I_SLICE pictures to use light_pd0 path
#if FTR_VLPD0 && !CLN_MERGE_LPD0_VLPD0
                        if (context_ptr->md_context->pd0_level == VERY_LIGHT_PD0) {

                            // [PD_PASS_0] Signal(s) derivation
#if FTR_VLPD0_INTER_DEPTH
                            signal_derivation_enc_dec_kernel_oq_very_light_pd0(pcs_ptr, context_ptr->md_context);
#else
                            signal_derivation_enc_dec_kernel_oq_very_light_pd0(context_ptr->md_context);
#endif
                            // Save a clean copy of the neighbor arrays
                            if (!context_ptr->md_context->skip_intra)
                                copy_neighbour_arrays_light_pd0(pcs_ptr,
                                    context_ptr->md_context,
                                    MD_NEIGHBOR_ARRAY_INDEX,
                                    MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                    0,
                                    sb_origin_x,
                                    sb_origin_y);

                            // Build the t=0 cand_block_array
                            build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);

                            mode_decision_sb_very_light_pd0(scs_ptr,
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
                                copy_neighbour_arrays_light_pd0(pcs_ptr,
                                    context_ptr->md_context,
                                    MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                    MD_NEIGHBOR_ARRAY_INDEX,
                                    0,
                                    sb_origin_x,
                                    sb_origin_y);
                        } else if(context_ptr->md_context->pd0_level != REGULAR_PD0) {
#else
#if CLN_MERGE_LPD0_VLPD0
                        if (context_ptr->md_context->pd0_level != REGULAR_PD0) {
#else
                        if (context_ptr->md_context->use_light_pd0) {
#endif
#endif
#if OPT_PD0_PATH

                            // [PD_PASS_0] Signal(s) derivation
                            signal_derivation_enc_dec_kernel_oq_light_pd0(scs_ptr, pcs_ptr, context_ptr->md_context);

                            // Save a clean copy of the neighbor arrays
                            if (!context_ptr->md_context->skip_intra)
                                copy_neighbour_arrays_light_pd0(pcs_ptr,
                                    context_ptr->md_context,
                                    MD_NEIGHBOR_ARRAY_INDEX,
                                    MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                    0,
                                    sb_origin_x,
                                    sb_origin_y);

                            // Build the t=0 cand_block_array
                            build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif
                            mode_decision_sb_light_pd0(scs_ptr,
                                pcs_ptr,
                                mdc_ptr,
                                sb_ptr,
                                sb_origin_x,
                                sb_origin_y,
                                sb_index,
                                context_ptr->md_context);
#if OPT_PD0_PATH
                            // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                            // Reset neighnor information to current SB @ position (0,0)
                            if (!context_ptr->md_context->skip_intra)
                                copy_neighbour_arrays_light_pd0(pcs_ptr,
                                    context_ptr->md_context,
                                    MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                    MD_NEIGHBOR_ARRAY_INDEX,
                                    0,
                                    sb_origin_x,
                                    sb_origin_y);
#endif
                        }
                        else {
#endif
#if OPT_PD0_PATH

                            // [PD_PASS_0] Signal(s) derivation
                            signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);

                            // Save a clean copy of the neighbor arrays
                            copy_neighbour_arrays(pcs_ptr,
                                context_ptr->md_context,
                                MD_NEIGHBOR_ARRAY_INDEX,
                                MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                0,
                                sb_origin_x,
                                sb_origin_y);

                            // Build the t=0 cand_block_array
                            build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#endif
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
#if OPT_PD0_PATH
                            // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                            // Reset neighnor information to current SB @ position (0,0)
                            copy_neighbour_arrays(pcs_ptr,
                                context_ptr->md_context,
                                MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                MD_NEIGHBOR_ARRAY_INDEX,
                                0,
                                sb_origin_x,
                                sb_origin_y);
#endif
#if LIGHT_PD0
                        }
#endif
#if FTR_LPD1_DETECTOR
#if FTR_VLPD1
#if FTR_VLPD0
                        // This classifier is used for only pd0_level 0 and pd0_level 1
                        // where the count_non_zero_coeffs is derived @ PD0
                        if (context_ptr->md_context->pd0_level != VERY_LIGHT_PD0)
#if CLN_LPD1_LVLS
                            lpd1_detector_post_pd0(pcs_ptr, md_ctx);
#endif
#endif
#if !CLN_LPD1_LVLS
                        for (int pd1_lvl = LPD1_LEVELS - 1; pd1_lvl > REGULAR_PD1; pd1_lvl--) {
                            if (md_ctx->lpd1_ctrls.pd1_level == pd1_lvl) {
                                if (md_ctx->lpd1_ctrls.use_lpd1_detector[pd1_lvl]) {
#if FIX_LPD1_FOR_ISLICE
                                    // Use info from ref frames (if available)
                                    if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl] && pcs_ptr->slice_type != I_SLICE) {
#else
                                    // Use info from ref frames
                                    if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl]) {
#endif
                                        EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                                        uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
                                        if (pcs_ptr->slice_type == B_SLICE) {
                                            EbReferenceObject *ref_obj_l1 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                                            l1_was_intra = ref_obj_l1->sb_intra[md_ctx->sb_index];
                                        }
                                        if (l0_was_intra && l1_was_intra) {
                                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                            continue;
                                        }
                                        else if (l0_was_intra || l1_was_intra) {
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
                                    const uint32_t nz_coeffs = md_ctx->avail_blk_flag[0] == EB_TRUE ? md_ctx->md_local_blk_unit[0].count_non_zero_coeffs : (uint32_t)~0;

                                    const uint32_t lambda = md_ctx->full_sb_lambda_md[EB_8_BIT_MD]; // light-PD1 assumes 8-bit MD
                                    const uint32_t rate = md_ctx->lpd1_ctrls.coeff_th[pd1_lvl];
                                    const uint32_t dist = md_ctx->lpd1_ctrls.cost_th_dist[pd1_lvl];
                                    /* dist << 14 is equivalent to 64 * 64 * 4 * dist (64 * 64 so the distortion is the per-pixel SSD) and 4 because
                                    the distortion of the 64x64 block is shifted by 2 (same as multiplying by 4) in perform_tx_light_pd0. */
                                    const uint64_t low_th = RDCOST(lambda, 6000 + rate * 500, dist << 14);
                                    const uint16_t coeff_th = md_ctx->lpd1_ctrls.coeff_th[pd1_lvl];

                                    // If the PD0 cost is very high and the number of non-zero coeffs is high, the block is difficult, so should use regular PD1
                                    if (pd0_cost > low_th && nz_coeffs > coeff_th) {
                                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                    }

                                    // If the best PD0 mode was INTER, check the MV length
                                    if (md_ctx->md_blk_arr_nsq[0].prediction_mode_flag == INTER_MODE && md_ctx->lpd1_ctrls.max_mv_length[pd1_lvl] != (uint16_t)~0 &&
                                        md_ctx->avail_blk_flag[0] == EB_TRUE) {


                                        PredictionUnit* pu_ptr = md_ctx->md_blk_arr_nsq[0].prediction_unit_array;
                                        const uint16_t max_mv_length = md_ctx->lpd1_ctrls.max_mv_length[pd1_lvl];

                                        if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0)
                                        {
                                            if (pu_ptr->mv[REF_LIST_0].x > max_mv_length || pu_ptr->mv[REF_LIST_0].y > max_mv_length)
                                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                        }
                                        else if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1)
                                        {
                                            if (pu_ptr->mv[REF_LIST_1].x > max_mv_length || pu_ptr->mv[REF_LIST_1].y > max_mv_length)
                                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                        }
                                        else //if (pu_ptr->inter_pred_direction_index == BI_PRED)
                                        {
                                            assert(pu_ptr->inter_pred_direction_index == BI_PRED);
                                            if (pu_ptr->mv[REF_LIST_0].x > max_mv_length || pu_ptr->mv[REF_LIST_0].y > max_mv_length)
                                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                            if (pu_ptr->mv[REF_LIST_1].x > max_mv_length || pu_ptr->mv[REF_LIST_1].y > max_mv_length)
                                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                        }
                                    }
#if FIX_LPD1_FOR_ISLICE
                                    if (pcs_ptr->slice_type != I_SLICE) {
                                        if (pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[sb_index] > (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] >> 5) * pcs_ptr->picture_qp)
                                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                    }
#else
                                    if (pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[sb_index] > (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] >> 5) * pcs_ptr->picture_qp)
                                        md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
#endif
                                }
                            }
                        }
#endif
                        // Force pred depth only for modes where that is not the default
                        if (md_ctx->lpd1_ctrls.pd1_level > REGULAR_PD1 && !ppcs->sc_class1) {
                            set_depth_ctrls(md_ctx, 0);
                            md_ctx->pred_depth_only = 1;
                        }
#else
                        if (md_ctx->lpd1_ctrls.enabled) {
                            if (md_ctx->lpd1_ctrls.use_light_pd1) {
                                if (md_ctx->lpd1_ctrls.use_lpd1_detector) {
                                    /* Use the cost and coeffs of the 64x64 block to avoid looping over all tested blocks to find
                                    the selected partitioning. */
                                    const uint64_t pd0_cost = md_ctx->md_local_blk_unit[0].cost;
                                    // If block was not tested in PD0, won't have coeff info, so set to max and base detection on cost only (which is set
                                    // even if 64x64 block is not tested)
                                    const uint32_t nz_coeffs = md_ctx->avail_blk_flag[0] == EB_TRUE ? md_ctx->md_local_blk_unit[0].count_non_zero_coeffs : (uint32_t)~0;

                                    const uint32_t lambda = md_ctx->full_sb_lambda_md[EB_8_BIT_MD]; // light-PD1 assumes 8-bit MD
                                    const uint32_t rate = md_ctx->lpd1_ctrls.coeff_th;
                                    const uint32_t dist = md_ctx->lpd1_ctrls.cost_th_dist;
                                    /* dist << 14 is equivalent to 64 * 64 * 4 * dist (64 * 64 so the distortion is the per-pixel SSD) and 4 because
                                    the distortion of the 64x64 block is shifted by 2 (same as multiplying by 4) in perform_tx_light_pd0. */
                                    const uint64_t low_th = RDCOST(lambda, 6000 + rate * 500, dist << 14);
                                    const uint16_t coeff_th = md_ctx->lpd1_ctrls.coeff_th;

                                    // If the PD0 cost is very high and the number of non-zero coeffs is high, the block is difficult, so should use regular PD1
                                    if (pd0_cost > low_th && nz_coeffs > coeff_th) {
                                        md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                    }
#if OPT_USE_MVS_LPD1_DETECTOR
                                    // If the best PD0 mode was INTER, check the MV length
                                    if (md_ctx->md_blk_arr_nsq[0].prediction_mode_flag == INTER_MODE && md_ctx->lpd1_ctrls.max_mv_length != (uint16_t)~0 &&
                                        md_ctx->avail_blk_flag[0] == EB_TRUE) {


                                        PredictionUnit* pu_ptr = md_ctx->md_blk_arr_nsq[0].prediction_unit_array;
                                        const uint16_t max_mv_length = md_ctx->lpd1_ctrls.max_mv_length;

                                        if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0)
                                        {
                                            if (pu_ptr->mv[REF_LIST_0].x > max_mv_length || pu_ptr->mv[REF_LIST_0].y > max_mv_length)
                                                md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                        }
                                        else if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1)
                                        {
                                            if (pu_ptr->mv[REF_LIST_1].x > max_mv_length || pu_ptr->mv[REF_LIST_1].y > max_mv_length)
                                                md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                        }
                                        else //if (pu_ptr->inter_pred_direction_index == BI_PRED)
                                        {
                                            assert(pu_ptr->inter_pred_direction_index == BI_PRED);
                                            if (pu_ptr->mv[REF_LIST_0].x > max_mv_length || pu_ptr->mv[REF_LIST_0].y > max_mv_length)
                                                md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                            if (pu_ptr->mv[REF_LIST_1].x > max_mv_length || pu_ptr->mv[REF_LIST_1].y > max_mv_length)
                                                md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                        }
                                    }
#endif
#if OPT_LIGHT_PD1_USE_ME_DIST_VAR
                                    if (pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[sb_index] > md_ctx->lpd1_ctrls.me_8x8_cost_variance_th)
                                        md_ctx->lpd1_ctrls.use_light_pd1 = 0;
#endif
                                    // Force pred depth only for modes where that is not the default
                                    if (md_ctx->lpd1_ctrls.use_light_pd1 && !ppcs->sc_class1) {
                                        set_depth_ctrls(md_ctx, 0);
                                        md_ctx->pred_depth_only = 1;
                                    }
                                }
                            }
                            //else {
                                // TODO: add opportunistic detector here
                            //}
                        }
#endif
#endif
                        // Perform Pred_0 depth refinement - add depth(s) to be considered in the next stage(s)
                        perform_pred_depth_refinement(
                            scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#if !OPT_PD0_PATH
                        // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]
                        // Reset neighnor information to current SB @ position (0,0)
                        copy_neighbour_arrays(pcs_ptr,
                                              context_ptr->md_context,
                                              MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                                              MD_NEIGHBOR_ARRAY_INDEX,
                                              0,
                                              sb_origin_x,
                                              sb_origin_y);
#endif
#if !FIX_REMOVE_PD1
                        if (pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_1 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_2 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_3 ||
                            pcs_ptr->parent_pcs_ptr->multi_pass_pd_level == MULTI_PASS_PD_LEVEL_4) {
                            // [PD_PASS_1] Signal(s) derivation
                            context_ptr->md_context->pd_pass = PD_PASS_1;
                            signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);
                            // Re-build mdc_blk_ptr for the 2nd PD Pass [PD_PASS_1]

                            build_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index, pcs_ptr->parent_pcs_ptr->sb_params_array[sb_index].is_complete_sb);
                            // Initialize avail_blk_flag to false
                            init_avail_blk_flag(scs_ptr, context_ptr->md_context);

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
#endif
                    }
#if FIX_REMOVE_PD1
                    // [PD_PASS_1] Signal(s) derivation
                    context_ptr->md_context->pd_pass = PD_PASS_1;
#if LIGHT_PD1_MACRO
#if FTR_LPD1_DETECTOR
#if FTR_VLPD1
#if FTR_VLPD0
                    // This classifier is used for the case PD0 is bypassed and for pd0_level 2
                    // where the count_non_zero_coeffs is not derived @ PD0
                    if (skip_pd_pass_0 || context_ptr->md_context->pd0_level == VERY_LIGHT_PD0) {
#else
                    if (skip_pd_pass_0) {
#endif
#if CLN_LPD1_LVLS
                        lpd1_detector_skip_pd0(pcs_ptr, md_ctx, pic_width_in_sb);
#else
                        const int16_t left_sb_index = (int16_t)md_ctx->sb_index - 1;
                        const int16_t top_sb_index = (int16_t)md_ctx->sb_index - (int16_t)pic_width_in_sb;

                        for (int pd1_lvl = LPD1_LEVELS - 1; pd1_lvl > REGULAR_PD1; pd1_lvl--) {
                            if (md_ctx->lpd1_ctrls.pd1_level == pd1_lvl) {
                                if (md_ctx->lpd1_ctrls.use_lpd1_detector[pd1_lvl]) {
#if FIX_LPD1_FOR_ISLICE
                                    // Use info from ref. frames (if available)
                                    if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl] && pcs_ptr->slice_type != I_SLICE) {
#else
                                    // Use info from ref. frames
                                    if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl]) {
#endif
                                        EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                                        uint8_t l0_was_intra = ref_obj_l0->sb_intra[md_ctx->sb_index], l1_was_intra = 0;
                                        uint8_t l0_was_skip = ref_obj_l0->sb_skip[md_ctx->sb_index], l1_was_skip = 1;
                                        if (pcs_ptr->slice_type == B_SLICE) {
                                            EbReferenceObject *ref_obj_l1 = (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
                                            l1_was_intra = ref_obj_l1->sb_intra[md_ctx->sb_index];
                                            l1_was_skip = ref_obj_l1->sb_skip[md_ctx->sb_index];
                                        }

                                        if (l0_was_intra && l1_was_intra) {
                                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                            continue;
                                        }
                                        else if (!l0_was_skip && !l1_was_skip && (l0_was_intra || l1_was_intra)) {
                                            md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                            continue;
                                        }
                                    }
#if FIX_LPD1_FOR_ISLICE
                                    // I_SLICE doesn't have ME info
                                    if (pcs_ptr->slice_type != I_SLICE) {
#endif
                                        // If the SB origin of one dimension is zero, then this SB is the first block in a row/column, so won't have neighbours
                                        if (md_ctx->sb_origin_x == 0 || md_ctx->sb_origin_y == 0) {
                                            if (pcs_ptr->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] > md_ctx->lpd1_ctrls.skip_pd0_edge_dist_th[pd1_lvl])
                                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                            else if (pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[sb_index] > (md_ctx->lpd1_ctrls.me_8x8_cost_variance_th[pd1_lvl] >> 5) * pcs_ptr->picture_qp)
                                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                        }
                                        else {
                                            if (pcs_ptr->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] >
                                                ((pcs_ptr->parent_pcs_ptr->me_64x64_distortion[left_sb_index] + pcs_ptr->parent_pcs_ptr->me_64x64_distortion[top_sb_index]) << md_ctx->lpd1_ctrls.skip_pd0_me_shift[pd1_lvl]))
                                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                            else if (pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[md_ctx->sb_index] >
                                                ((pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[left_sb_index] + pcs_ptr->parent_pcs_ptr->me_8x8_cost_variance[top_sb_index]) << md_ctx->lpd1_ctrls.skip_pd0_me_shift[pd1_lvl])) {
                                                md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                            }
                                            else if (md_ctx->lpd1_ctrls.use_ref_info[pd1_lvl]) {
                                                // Use info from neighbouring SBs
                                                if (pcs_ptr->sb_intra[left_sb_index] && pcs_ptr->sb_intra[top_sb_index]) {
                                                    md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                                }
                                                else if (!pcs_ptr->sb_skip[left_sb_index] && !pcs_ptr->sb_intra[top_sb_index] &&
                                                    (pcs_ptr->sb_intra[left_sb_index] || pcs_ptr->sb_intra[top_sb_index])) {
                                                    md_ctx->lpd1_ctrls.pd1_level = pd1_lvl - 1;
                                                }
                                            }
                                        }
#if FIX_LPD1_FOR_ISLICE
                                    }
#endif
                                }
                            }
                        }
#endif
                    }

                    // Can only use light-PD1 under the following conditions
                    if (!(md_ctx->hbd_mode_decision == 0 &&
#if !OPT_LPD1_MRP // NB reference pruning not done in light-PD1
                        ppcs->frm_hdr.allow_warped_motion == 0 &&
                        ppcs->ref_list0_count_try == 1 &&
                        ppcs->ref_list1_count_try == 1 &&
#endif
                        md_ctx->pred_depth_only &&
                        ppcs->disallow_nsq == EB_TRUE &&
                        md_ctx->disallow_4x4 == EB_TRUE &&
                        scs_ptr->static_config.super_block_size == 64)) {

                        md_ctx->lpd1_ctrls.pd1_level = REGULAR_PD1;
                    }
#else
                    if (skip_pd_pass_0 && md_ctx->lpd1_ctrls.enabled) {
                        if (md_ctx->lpd1_ctrls.use_light_pd1) {
                            if (md_ctx->lpd1_ctrls.use_lpd1_detector) {

                                const int16_t left_sb_index = (int16_t)md_ctx->sb_index - 1;
                                const int16_t top_sb_index = (int16_t)md_ctx->sb_index - (int16_t)pic_width_in_sb;
#if TUNE_LPD1_DETECTOR
#if FIX_LPD1_DETECTOR
                                // If the SB origin of one dimension is zero, then this SB is the first block in a row/column, so won't have neighbours
                                if (md_ctx->sb_origin_x == 0 || md_ctx->sb_origin_y == 0) {
#else
                                if (left_sb_index < 0 || top_sb_index < 0) {
#endif
                                    if (pcs_ptr->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] > md_ctx->lpd1_ctrls.skip_pd0_edge_dist_th)
                                        md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                }
                                else if (pcs_ptr->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] >
                                    ((pcs_ptr->parent_pcs_ptr->me_64x64_distortion[left_sb_index] + pcs_ptr->parent_pcs_ptr->me_64x64_distortion[top_sb_index]) << md_ctx->lpd1_ctrls.skip_pd0_me_dist_shift)) {
                                    md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                }
#else
                                if (left_sb_index < 0 || top_sb_index < 0) {
                                    if (pcs_ptr->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] > 1024)
                                        md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                }
                                else if (pcs_ptr->parent_pcs_ptr->me_64x64_distortion[md_ctx->sb_index] >
                                    ((pcs_ptr->parent_pcs_ptr->me_64x64_distortion[left_sb_index] + pcs_ptr->parent_pcs_ptr->me_64x64_distortion[top_sb_index]) << 1)) {
                                    md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                                }
#endif
                            }
                        }
                        //else {
                            // TODO: add opportunistic detector here
                        //}

                    }

                    // Can only use light-PD1 under the following conditions
                    if (!(md_ctx->hbd_mode_decision == 0 &&
#if !OPT_LPD1_MRP // NB reference pruning not done in light-PD1
                        ppcs->frm_hdr.allow_warped_motion == 0 &&
                        ppcs->ref_list0_count_try == 1 &&
                        ppcs->ref_list1_count_try == 1 &&
#endif
                        md_ctx->pred_depth_only &&
                        ppcs->disallow_nsq == EB_TRUE &&
                        md_ctx->disallow_4x4 == EB_TRUE &&
                        scs_ptr->static_config.super_block_size == 64)) {

                        md_ctx->lpd1_ctrls.use_light_pd1 = 0;
                    }
#endif
#else
                    // Can only use light-PD1 under the following conditions
                    if (md_ctx->hbd_mode_decision == 0 &&
#if !OPT_LPD1_MRP // NB reference pruning not done in light-PD1
                        ppcs->frm_hdr.allow_warped_motion == 0 &&
                        ppcs->ref_list0_count_try == 1 &&
                        ppcs->ref_list1_count_try == 1 &&
#endif
                        md_ctx->pred_depth_only &&
                        ppcs->disallow_nsq == EB_TRUE &&
                        md_ctx->disallow_4x4 == EB_TRUE &&
                        scs_ptr->static_config.super_block_size == 64) {

                        if (pcs_ptr->enc_mode <= ENC_M10)
                            md_ctx->use_light_pd1 = 0;
                        else
                            md_ctx->use_light_pd1 = (scs_ptr->input_resolution <= INPUT_SIZE_720p_RANGE) ? (pcs_ptr->temporal_layer_index == 0 ? 0 : 1) : (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ? 0 : 1);
                    }
                    else {
                        md_ctx->use_light_pd1 = 0;
                    }
#endif
#if FTR_VLPD1
                    exaustive_light_pd1_features(md_ctx, ppcs, md_ctx->lpd1_ctrls.pd1_level > REGULAR_PD1, 0);
#if CLN_LPD1_LVLS
                    if (md_ctx->lpd1_ctrls.pd1_level > REGULAR_PD1)
                        signal_derivation_enc_dec_kernel_oq_light_pd1(pcs_ptr, context_ptr->md_context);
#else
                    if (md_ctx->lpd1_ctrls.pd1_level == VERY_LIGHT_PD1)
                        signal_derivation_enc_dec_kernel_oq_very_light_pd1(pcs_ptr, context_ptr->md_context);
                    else if (md_ctx->lpd1_ctrls.pd1_level == LIGHT_PD1)
                        signal_derivation_enc_dec_kernel_oq_light_pd1(pcs_ptr, context_ptr->md_context);
#endif
                    else
                        signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);
#else
                    exaustive_light_pd1_features(md_ctx, ppcs, md_ctx->lpd1_ctrls.use_light_pd1,0);

                    if (md_ctx->lpd1_ctrls.use_light_pd1)
                        signal_derivation_enc_dec_kernel_oq_light_pd1(pcs_ptr, context_ptr->md_context);
                    else
                        signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);
#endif
#else
                    signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);
#endif
#else
                    // [PD_PASS_2] Signal(s) derivation
                    context_ptr->md_context->pd_pass = PD_PASS_2;
                        signal_derivation_enc_dec_kernel_oq(scs_ptr, pcs_ptr, context_ptr->md_context);
                    // Re-build mdc_blk_ptr for the 3rd PD Pass [PD_PASS_2]
#endif
#if FIX_REMOVE_PD1
                    if (!skip_pd_pass_0 && pcs_ptr->parent_pcs_ptr->multi_pass_pd_level != MULTI_PASS_PD_OFF)
#else
                    if (!pd_pass_2_only && pcs_ptr->parent_pcs_ptr->multi_pass_pd_level != MULTI_PASS_PD_OFF)
#endif
                        build_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index, pcs_ptr->parent_pcs_ptr->sb_params_array[sb_index].is_complete_sb);
                    else
                        // Build the t=0 cand_block_array
                        build_starting_cand_block_array(scs_ptr, pcs_ptr, context_ptr->md_context, sb_index);
#if !OPT_PD0_PATH
                    // Initialize avail_blk_flag to false
                    init_avail_blk_flag(scs_ptr, context_ptr->md_context);
#endif
#if FIX_REMOVE_PD1
                    // [PD_PASS_1] Mode Decision - Obtain the final partitioning decision using more accurate info
                    // than previous stages.  Reduce the total number of partitions to 1.
                    // Input : mdc_blk_ptr built @ PD0 refinement
                    // Output: md_blk_arr_nsq reduced set of block(s)

                    // PD1 MD Tool(s): default MD Tool(s)
#else
                    // [PD_PASS_2] Mode Decision - Obtain the final partitioning decision using more accurate info
                    // than previous stages.  Reduce the total number of partitions to 1.
                    // Input : mdc_blk_ptr built @ PD1 refinement
                    // Output: md_blk_arr_nsq reduced set of block(s)

                    // PD2 MD Tool(s): default MD Tool(s)
#endif
#if LIGHT_PD1_MACRO
#if FTR_LPD1_DETECTOR
#if FTR_VLPD1
                    if (md_ctx->lpd1_ctrls.pd1_level > REGULAR_PD1)
#else
                    if (md_ctx->lpd1_ctrls.use_light_pd1)
#endif
#else
                    if (md_ctx->use_light_pd1)
#endif
#if FTR_10BIT_MDS3_LPD1
                        mode_decision_sb_light_pd1(scs_ptr,
                                        pcs_ptr,
#else
                        mode_decision_sb_light_pd1(pcs_ptr,
#endif
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
#else
                    mode_decision_sb(scs_ptr,
                                     pcs_ptr,
                                     mdc_ptr,
                                     sb_ptr,
                                     sb_origin_x,
                                     sb_origin_y,
                                     sb_index,
                                     context_ptr->md_context);
#endif
                    //if (/*ppcs->is_used_as_reference_flag &&*/ md_ctx->hbd_mode_decision == 0 && scs_ptr->static_config.encoder_bit_depth > EB_8BIT)
                    //    md_ctx->bypass_encdec = 0;
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
#if FTR_BYPASS_ENCDEC
                    if (!context_ptr->md_context->bypass_encdec) {
#endif
                        av1_encode_decode(
                            scs_ptr, pcs_ptr, sb_ptr, sb_index, sb_origin_x, sb_origin_y, context_ptr);
#if FTR_BYPASS_ENCDEC
                    }
#endif
#if REFCTR_SEP_ENCDEC
                    av1_encdec_update(scs_ptr, pcs_ptr, sb_ptr, sb_index, sb_origin_x, sb_origin_y, context_ptr);
#endif
#endif

                    context_ptr->coded_sb_count++;
                }
                x_sb_start_index = (x_sb_start_index > 0) ? x_sb_start_index - 1 : 0;
            }
        }

        svt_block_on_mutex(pcs_ptr->intra_mutex);
#if  FTR_INTRA_DETECTOR
        pcs_ptr->intra_coded_area += (uint32_t)context_ptr->tot_intra_coded_area;
#endif
#if FTR_COEFF_DETECTOR
        pcs_ptr->skip_coded_area += (uint32_t)context_ptr->tot_skip_coded_area;
#endif
        // Accumulate block selection
        pcs_ptr->enc_dec_coded_sb_count += (uint32_t)context_ptr->coded_sb_count;
        EbBool last_sb_flag = (pcs_ptr->sb_total_count_pix == pcs_ptr->enc_dec_coded_sb_count);
        svt_release_mutex(pcs_ptr->intra_mutex);

        if (last_sb_flag) {
            EbBool do_recode = EB_FALSE;
#if !FIX_SANITIZER_RACE_CONDS
            scs_ptr->encode_context_ptr->recode_loop = scs_ptr->static_config.recode_loop;
#endif
#if FTR_RC_CAP
            if ((use_input_stat(scs_ptr) || scs_ptr->lap_enabled || scs_ptr->static_config.max_bit_rate != 0) &&
#else
            if ((use_input_stat(scs_ptr) || scs_ptr->lap_enabled) &&
#endif
                scs_ptr->encode_context_ptr->recode_loop != DISALLOW_RECODE) {
                recode_loop_decision_maker(pcs_ptr, scs_ptr, &do_recode);
            }

            if (do_recode) {

                pcs_ptr->enc_dec_coded_sb_count = 0;
#if !CLN_MDCONTEXT
                // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
                if (context_ptr->is_md_rate_estimation_ptr_owner) {
                    EB_FREE_ARRAY(context_ptr->md_rate_estimation_ptr);
                    context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
                }
                context_ptr->md_rate_estimation_ptr = pcs_ptr->md_rate_estimation_array;
#endif
                // re-init mode decision configuration for qp update for re-encode frame
                mode_decision_configuration_init_qp_update(pcs_ptr);
                // init segment for re-encode frame
                init_enc_dec_segement(pcs_ptr->parent_pcs_ptr);
                EbObjectWrapper *enc_dec_re_encode_tasks_wrapper_ptr;
                uint16_t tg_count =
                    pcs_ptr->parent_pcs_ptr->tile_group_cols * pcs_ptr->parent_pcs_ptr->tile_group_rows;
                for (uint16_t tile_group_idx = 0; tile_group_idx < tg_count; tile_group_idx++) {
                    svt_get_empty_object(context_ptr->enc_dec_feedback_fifo_ptr,
                            &enc_dec_re_encode_tasks_wrapper_ptr);

                    EncDecTasks *enc_dec_re_encode_tasks_ptr = (EncDecTasks *)enc_dec_re_encode_tasks_wrapper_ptr->object_ptr;
                    enc_dec_re_encode_tasks_ptr->pcs_wrapper_ptr  = enc_dec_tasks_ptr->pcs_wrapper_ptr;
                    enc_dec_re_encode_tasks_ptr->input_type       = ENCDEC_TASKS_MDC_INPUT;
                    enc_dec_re_encode_tasks_ptr->tile_group_index = tile_group_idx;

                    // Post the Full Results Object
                    svt_post_full_object(enc_dec_re_encode_tasks_wrapper_ptr);
                }

            }
            else {
#if OPT_MEMORY_CDF
            EB_FREE_ARRAY(pcs_ptr->ec_ctx_array);
#endif
            // Copy film grain data from parent picture set to the reference object for further reference
            if (scs_ptr->seq_header.film_grain_params_present) {
                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE &&
                    pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr) {
                    ((EbReferenceObject *)
                         pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                        ->film_grain_params = pcs_ptr->parent_pcs_ptr->frm_hdr.film_grain_params;
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
                        ->global_motion[frame] = pcs_ptr->parent_pcs_ptr->global_motion[frame];
#if CLN_MDCONTEXT
            svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->sgrproj_restore_cost,
                pcs_ptr->md_rate_estimation_array->sgrproj_restore_fac_bits,
                2 * sizeof(int32_t));
            svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->switchable_restore_cost,
                pcs_ptr->md_rate_estimation_array->switchable_restore_fac_bits,
                3 * sizeof(int32_t));
            svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->wiener_restore_cost,
                pcs_ptr->md_rate_estimation_array->wiener_restore_fac_bits,
                2 * sizeof(int32_t));
#else
            svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->sgrproj_restore_cost,
                       context_ptr->md_rate_estimation_ptr->sgrproj_restore_fac_bits,
                       2 * sizeof(int32_t));
            svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->switchable_restore_cost,
                       context_ptr->md_rate_estimation_ptr->switchable_restore_fac_bits,
                       3 * sizeof(int32_t));
            svt_memcpy(pcs_ptr->parent_pcs_ptr->av1x->wiener_restore_cost,
                       context_ptr->md_rate_estimation_ptr->wiener_restore_fac_bits,
                       2 * sizeof(int32_t));
#endif
            pcs_ptr->parent_pcs_ptr->av1x->rdmult =
                context_ptr->pic_full_lambda[(context_ptr->bit_depth == EB_10BIT) ? EB_10_BIT_MD
                                                                                  : EB_8_BIT_MD];
            if (pcs_ptr->parent_pcs_ptr->superres_total_recode_loop == 0) {
                svt_release_object(pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr);
                pcs_ptr->parent_pcs_ptr->me_data_wrapper_ptr = (EbObjectWrapper*)NULL;
                pcs_ptr->parent_pcs_ptr->pa_me_data = NULL;
            }
            // Get Empty EncDec Results
            svt_get_empty_object(context_ptr->enc_dec_output_fifo_ptr, &enc_dec_results_wrapper_ptr);
            enc_dec_results_ptr = (EncDecResults *)enc_dec_results_wrapper_ptr->object_ptr;
            enc_dec_results_ptr->pcs_wrapper_ptr = enc_dec_tasks_ptr->pcs_wrapper_ptr;
            //CHKN these are not needed for DLF
            enc_dec_results_ptr->completed_sb_row_index_start = 0;
            enc_dec_results_ptr->completed_sb_row_count =
                ((pcs_ptr->parent_pcs_ptr->aligned_height + scs_ptr->sb_size_pix - 1) >> sb_size_log2);
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
