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

#include "EbEncHandle.h"
#include "EbEncDecTasks.h"
#include "EbEncDecResults.h"
#include "EbCodingLoop.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbUtility.h"
#include "grainSynthesis.h"


#if MULTI_PASS_PD
#define FC_SKIP_TX_SR_TH025 125 // Fast cost skip tx search threshold.
#define FC_SKIP_TX_SR_TH010 110 // Fast cost skip tx search threshold.
#endif

void eb_av1_cdef_search(
    EncDecContext                *context_ptr,
    SequenceControlSet           *sequence_control_set_ptr,
    PictureControlSet            *picture_control_set_ptr
);

void av1_cdef_frame16bit(
    uint8_t is16bit,
    SequenceControlSet           *sequence_control_set_ptr,
    PictureControlSet            *pCs
);

void eb_av1_add_film_grain(EbPictureBufferDesc *src,
    EbPictureBufferDesc *dst,
    aom_film_grain_t *film_grain_ptr);

void eb_av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm, int32_t after_cdef);
void eb_av1_pick_filter_restoration(const Yv12BufferConfig *src, Yv12BufferConfig * trial_frame_rst /*Av1Comp *cpi*/, Macroblock *x, Av1Common *const cm);
void eb_av1_loop_restoration_filter_frame(Yv12BufferConfig *frame, Av1Common *cm, int32_t optimized_lr);

static void enc_dec_context_dctor(EbPtr p)
{
    EbThreadContext* thread_context_ptr = (EbThreadContext*)p;
    EncDecContext* obj = (EncDecContext*)thread_context_ptr->priv;
    EB_DELETE(obj->md_context);
    EB_DELETE(obj->residual_buffer);
    EB_DELETE(obj->transform_buffer);
    EB_DELETE(obj->inverse_quant_buffer);
    EB_DELETE(obj->input_sample16bit_buffer);
    if (obj->is_md_rate_estimation_ptr_owner)
        EB_FREE(obj->md_rate_estimation_ptr);
    EB_FREE_ARRAY(obj->transform_inner_array_ptr);
    EB_FREE_ARRAY(obj);
}

/******************************************************
 * Enc Dec Context Constructor
 ******************************************************/
EbErrorType enc_dec_context_ctor(
    EbThreadContext     *thread_context_ptr,
    const EbEncHandle   *enc_handle_ptr,
    int                 index,
    int                 tasks_index,
    int                 demux_index)
{
    const EbSvtAv1EncConfiguration* static_config = &enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr->static_config;
    EbBool is16bit = (EbBool)(static_config->encoder_bit_depth > EB_8BIT);
    EbColorFormat color_format = static_config->encoder_color_format;
    uint8_t enable_hbd_mode_decision = static_config->enable_hbd_mode_decision;

    EncDecContext         *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv = context_ptr;
    thread_context_ptr->dctor = enc_dec_context_dctor;

    context_ptr->is16bit = is16bit;
    context_ptr->color_format = color_format;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->enc_dec_tasks_resource_ptr, index);
    context_ptr->enc_dec_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->enc_dec_results_resource_ptr, index);
    context_ptr->enc_dec_feedback_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->enc_dec_tasks_resource_ptr, tasks_index);
    context_ptr->picture_demux_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);

    // Trasform Scratch Memory
    EB_MALLOC_ARRAY(context_ptr->transform_inner_array_ptr, 3152); //refer to EbInvTransform_SSE2.as. case 32x32
    // MD rate Estimation tables
    EB_MALLOC(context_ptr->md_rate_estimation_ptr, sizeof(MdRateEstimationContext));
    context_ptr->is_md_rate_estimation_ptr_owner = EB_TRUE;

    // Prediction Buffer
    {
        EbPictureBufferDescInitData initData;

        initData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        initData.max_width = SB_STRIDE_Y;
        initData.max_height = SB_STRIDE_Y;
        initData.bit_depth = EB_8BIT;
        initData.left_padding = 0;
        initData.right_padding = 0;
        initData.top_padding = 0;
        initData.bot_padding = 0;
        initData.split_mode = EB_FALSE;
        initData.color_format = color_format;

        context_ptr->input_sample16bit_buffer = (EbPictureBufferDesc *)EB_NULL;
        if (is16bit) {
            initData.bit_depth = EB_16BIT;

            EB_NEW(
                context_ptr->input_sample16bit_buffer,
                eb_picture_buffer_desc_ctor,
                (EbPtr)&initData);
        }
    }

    // Scratch Coeff Buffer
    {
        EbPictureBufferDescInitData initData;

        initData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        initData.max_width = SB_STRIDE_Y;
        initData.max_height = SB_STRIDE_Y;
        initData.bit_depth = EB_16BIT;
        initData.color_format = color_format;
        initData.left_padding = 0;
        initData.right_padding = 0;
        initData.top_padding = 0;
        initData.bot_padding = 0;
        initData.split_mode = EB_FALSE;

        EbPictureBufferDescInitData init32BitData;

        init32BitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        init32BitData.max_width = SB_STRIDE_Y;
        init32BitData.max_height = SB_STRIDE_Y;
        init32BitData.bit_depth = EB_32BIT;
        init32BitData.color_format = color_format;
        init32BitData.left_padding = 0;
        init32BitData.right_padding = 0;
        init32BitData.top_padding = 0;
        init32BitData.bot_padding = 0;
        init32BitData.split_mode = EB_FALSE;
        EB_NEW(
            context_ptr->inverse_quant_buffer,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&init32BitData);
        EB_NEW(
            context_ptr->transform_buffer,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&init32BitData);
        EB_NEW(
            context_ptr->residual_buffer,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&initData);
    }

    // Mode Decision Context
#if PAL_SUP
    EB_NEW(
        context_ptr->md_context,
        mode_decision_context_ctor,
        color_format, 0, 0, enable_hbd_mode_decision , static_config->screen_content_mode );
#else
    EB_NEW(
        context_ptr->md_context,
        mode_decision_context_ctor,
        color_format, 0, 0, enable_hbd_mode_decision);
#endif
    if (enable_hbd_mode_decision)
        context_ptr->md_context->input_sample16bit_buffer = context_ptr->input_sample16bit_buffer;

    context_ptr->md_context->enc_dec_context_ptr = context_ptr;

    return EB_ErrorNone;
}

/**************************************************
 * Reset Segmentation Map
 *************************************************/
static void reset_segmentation_map(SegmentationNeighborMap *segmentation_map){
    if(segmentation_map->data!=NULL)
        EB_MEMSET(segmentation_map->data, ~0, segmentation_map->map_size);
}

/**************************************************
 * Reset Mode Decision Neighbor Arrays
 *************************************************/
static void reset_encode_pass_neighbor_arrays(PictureControlSet *picture_control_set_ptr)
{
    neighbor_array_unit_reset(picture_control_set_ptr->ep_intra_luma_mode_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_intra_chroma_mode_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_mv_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_skip_flag_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_mode_type_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_leaf_depth_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_luma_recon_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_cb_recon_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_cr_recon_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_luma_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_cb_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_cr_dc_sign_level_coeff_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->ep_partition_context_neighbor_array);
    // TODO(Joel): 8-bit ep_luma_recon_neighbor_array (Cb,Cr) when is16bit==0?
    EbBool is16bit = (EbBool)(picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    if (is16bit) {
        neighbor_array_unit_reset(picture_control_set_ptr->ep_luma_recon_neighbor_array16bit);
        neighbor_array_unit_reset(picture_control_set_ptr->ep_cb_recon_neighbor_array16bit);
        neighbor_array_unit_reset(picture_control_set_ptr->ep_cr_recon_neighbor_array16bit);
    }
    return;
}

/**************************************************
 * Reset Coding Loop
 **************************************************/
static void reset_enc_dec(
    EncDecContext         *context_ptr,
    PictureControlSet     *picture_control_set_ptr,
    SequenceControlSet    *sequence_control_set_ptr,
    uint32_t                   segment_index)
{
    context_ptr->is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

#if ADD_DELTA_QP_SUPPORT
    uint16_t picture_qp = picture_control_set_ptr->picture_qp;
    context_ptr->qp = picture_qp;
    context_ptr->qp_index = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present ? (uint8_t)quantizer_to_qindex[context_ptr->qp] : (uint8_t)picture_control_set_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
#else
    context_ptr->qp = picture_control_set_ptr->picture_qp;
#endif
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation

    context_ptr->chroma_qp = (uint8_t)context_ptr->qp;

    // Lambda Assignement
    context_ptr->qp_index = (uint8_t)picture_control_set_ptr->
        parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    (*av1_lambda_assignment_function_table[picture_control_set_ptr->parent_pcs_ptr->pred_structure])(
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
        &context_ptr->fast_chroma_lambda,
        &context_ptr->full_chroma_lambda,
        (uint8_t)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index,
        context_ptr->md_context->hbd_mode_decision);
    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    if (context_ptr->is_md_rate_estimation_ptr_owner) {
        EB_FREE(context_ptr->md_rate_estimation_ptr);
        context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
    }
    context_ptr->md_rate_estimation_ptr = picture_control_set_ptr->md_rate_estimation_array;
    if (segment_index == 0){
        reset_encode_pass_neighbor_arrays(picture_control_set_ptr);
        reset_segmentation_map(picture_control_set_ptr->segmentation_neighbor_map);
    }

    return;
}

/******************************************************
 * EncDec Configure LCU
 ******************************************************/
static void enc_dec_configure_lcu(
    EncDecContext         *context_ptr,
    SuperBlock            *sb_ptr,
    PictureControlSet     *picture_control_set_ptr,
    uint8_t                    sb_qp)
{
    context_ptr->qp = sb_qp;

    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation
    context_ptr->chroma_qp = (uint8_t)context_ptr->qp;
    /* Note(CHKN) : when Qp modulation varies QP on a sub-LCU(CU) basis,  Lamda has to change based on Cu->QP , and then this code has to move inside the CU loop in MD */
    (void)sb_ptr;
    context_ptr->qp_index = (uint8_t)picture_control_set_ptr->
        parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    (*av1_lambda_assignment_function_table[picture_control_set_ptr->parent_pcs_ptr->pred_structure])(
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
        &context_ptr->fast_chroma_lambda,
        &context_ptr->full_chroma_lambda,
        (uint8_t)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index,
        context_ptr->md_context->hbd_mode_decision);

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
EbBool assign_enc_dec_segments(
    EncDecSegments   *segmentPtr,
    uint16_t             *segmentInOutIndex,
    EncDecTasks      *taskPtr,
    EbFifo           *srmFifoPtr)
{
    EbBool continueProcessingFlag = EB_FALSE;
    EbObjectWrapper *wrapper_ptr;
    EncDecTasks *feedbackTaskPtr;

    uint32_t rowSegmentIndex = 0;
    uint32_t segment_index;
    uint32_t rightSegmentIndex;
    uint32_t bottomLeftSegmentIndex;

    int16_t feedbackRowIndex = -1;

    uint32_t selfAssigned = EB_FALSE;

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
        *segmentInOutIndex = segmentPtr->row_array[0].current_seg_index;
        taskPtr->input_type = ENCDEC_TASKS_CONTINUE;
        ++segmentPtr->row_array[0].current_seg_index;
        continueProcessingFlag = EB_TRUE;

        //fprintf(trace, "Start  Pic: %u Seg: %u\n",
        //    (unsigned) ((PictureControlSet*) taskPtr->picture_control_set_wrapper_ptr->object_ptr)->picture_number,
        //    *segmentInOutIndex);

        break;

    case ENCDEC_TASKS_ENCDEC_INPUT:

        // Setup rowSegmentIndex to release the in_progress token
        //rowSegmentIndex = taskPtr->encDecSegmentRowArray[0];

        // Start on the assigned row immediately
        *segmentInOutIndex = segmentPtr->row_array[taskPtr->enc_dec_segment_row].current_seg_index;
        taskPtr->input_type = ENCDEC_TASKS_CONTINUE;
        ++segmentPtr->row_array[taskPtr->enc_dec_segment_row].current_seg_index;
        continueProcessingFlag = EB_TRUE;

        //fprintf(trace, "Start  Pic: %u Seg: %u\n",
        //    (unsigned) ((PictureControlSet*) taskPtr->picture_control_set_wrapper_ptr->object_ptr)->picture_number,
        //    *segmentInOutIndex);

        break;

    case ENCDEC_TASKS_CONTINUE:

        // Update the Dependency List for Right and Bottom Neighbors
        segment_index = *segmentInOutIndex;
        rowSegmentIndex = segment_index / segmentPtr->segment_band_count;

        rightSegmentIndex = segment_index + 1;
        bottomLeftSegmentIndex = segment_index + segmentPtr->segment_band_count;

        // Right Neighbor
        if (segment_index < segmentPtr->row_array[rowSegmentIndex].ending_seg_index)
        {
            eb_block_on_mutex(segmentPtr->row_array[rowSegmentIndex].assignment_mutex);

            --segmentPtr->dep_map.dependency_map[rightSegmentIndex];

            if (segmentPtr->dep_map.dependency_map[rightSegmentIndex] == 0) {
                *segmentInOutIndex = segmentPtr->row_array[rowSegmentIndex].current_seg_index;
                ++segmentPtr->row_array[rowSegmentIndex].current_seg_index;
                selfAssigned = EB_TRUE;
                continueProcessingFlag = EB_TRUE;

                //fprintf(trace, "Start  Pic: %u Seg: %u\n",
                //    (unsigned) ((PictureControlSet*) taskPtr->picture_control_set_wrapper_ptr->object_ptr)->picture_number,
                //    *segmentInOutIndex);
            }

            eb_release_mutex(segmentPtr->row_array[rowSegmentIndex].assignment_mutex);
        }

        // Bottom-left Neighbor
        if (rowSegmentIndex < segmentPtr->segment_row_count - 1 && bottomLeftSegmentIndex >= segmentPtr->row_array[rowSegmentIndex + 1].starting_seg_index)
        {
            eb_block_on_mutex(segmentPtr->row_array[rowSegmentIndex + 1].assignment_mutex);

            --segmentPtr->dep_map.dependency_map[bottomLeftSegmentIndex];

            if (segmentPtr->dep_map.dependency_map[bottomLeftSegmentIndex] == 0) {
                if (selfAssigned == EB_TRUE)
                    feedbackRowIndex = (int16_t)rowSegmentIndex + 1;
                else {
                    *segmentInOutIndex = segmentPtr->row_array[rowSegmentIndex + 1].current_seg_index;
                    ++segmentPtr->row_array[rowSegmentIndex + 1].current_seg_index;
                    selfAssigned = EB_TRUE;
                    continueProcessingFlag = EB_TRUE;

                    //fprintf(trace, "Start  Pic: %u Seg: %u\n",
                    //    (unsigned) ((PictureControlSet*) taskPtr->picture_control_set_wrapper_ptr->object_ptr)->picture_number,
                    //    *segmentInOutIndex);
                }
            }
            eb_release_mutex(segmentPtr->row_array[rowSegmentIndex + 1].assignment_mutex);
        }

        if (feedbackRowIndex > 0) {
            eb_get_empty_object(
                srmFifoPtr,
                &wrapper_ptr);
            feedbackTaskPtr = (EncDecTasks*)wrapper_ptr->object_ptr;
            feedbackTaskPtr->input_type = ENCDEC_TASKS_ENCDEC_INPUT;
            feedbackTaskPtr->enc_dec_segment_row = feedbackRowIndex;
            feedbackTaskPtr->picture_control_set_wrapper_ptr = taskPtr->picture_control_set_wrapper_ptr;
            eb_post_full_object(wrapper_ptr);
        }

        break;

    default:
        break;
    }

    return continueProcessingFlag;
}
void recon_output(
    PictureControlSet    *picture_control_set_ptr,
    SequenceControlSet   *sequence_control_set_ptr) {
    EbObjectWrapper             *outputReconWrapperPtr;
    EbBufferHeaderType           *outputReconPtr;
    EncodeContext               *encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;
    EbBool is16bit = (sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    // The totalNumberOfReconFrames counter has to be write/read protected as
    //   it is used to determine the end of the stream.  If it is not protected
    //   the encoder might not properly terminate.
    eb_block_on_mutex(encode_context_ptr->total_number_of_recon_frame_mutex);

    if (!picture_control_set_ptr->parent_pcs_ptr->is_alt_ref) {
        // Get Recon Buffer
        eb_get_empty_object(
            sequence_control_set_ptr->encode_context_ptr->recon_output_fifo_ptr,
            &outputReconWrapperPtr);
        outputReconPtr = (EbBufferHeaderType*)outputReconWrapperPtr->object_ptr;
        outputReconPtr->flags = 0;

        // START READ/WRITE PROTECTED SECTION
        if (encode_context_ptr->total_number_of_recon_frames == encode_context_ptr->terminating_picture_number)
            outputReconPtr->flags = EB_BUFFERFLAG_EOS;

        encode_context_ptr->total_number_of_recon_frames++;

        //eb_release_mutex(encode_context_ptr->terminating_conditions_mutex);

        // STOP READ/WRITE PROTECTED SECTION
        outputReconPtr->n_filled_len = 0;

        // Copy the Reconstructed Picture to the Output Recon Buffer
        {
            uint32_t sampleTotalCount;
            uint8_t *reconReadPtr;
            uint8_t *reconWritePtr;

            EbPictureBufferDesc *recon_ptr;
            {
                if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    recon_ptr = is16bit ?
                    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit :
                    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
                else {
                    if (is16bit)
                        recon_ptr = picture_control_set_ptr->recon_picture16bit_ptr;
                    else
                        recon_ptr = picture_control_set_ptr->recon_picture_ptr;
                }
            }

            // FGN: Create a buffer if needed, copy the reconstructed picture and run the film grain synthesis algorithm

            if (sequence_control_set_ptr->seq_header.film_grain_params_present) {
                EbPictureBufferDesc  *intermediateBufferPtr;
                {
                    if (is16bit)
                        intermediateBufferPtr = picture_control_set_ptr->film_grain_picture16bit_ptr;
                    else
                        intermediateBufferPtr = picture_control_set_ptr->film_grain_picture_ptr;
                }

                aom_film_grain_t *film_grain_ptr;

                if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    film_grain_ptr = &((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->film_grain_params;
                else
                    film_grain_ptr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr.film_grain_params;

                eb_av1_add_film_grain(recon_ptr, intermediateBufferPtr, film_grain_ptr);
                recon_ptr = intermediateBufferPtr;
            }

            // End running the film grain
            // Y Recon Samples
            sampleTotalCount = ((recon_ptr->max_width - sequence_control_set_ptr->max_input_pad_right) * (recon_ptr->max_height - sequence_control_set_ptr->max_input_pad_bottom)) << is16bit;
            reconReadPtr = recon_ptr->buffer_y + (recon_ptr->origin_y << is16bit) * recon_ptr->stride_y + (recon_ptr->origin_x << is16bit);
            reconWritePtr = &(outputReconPtr->p_buffer[outputReconPtr->n_filled_len]);

            CHECK_REPORT_ERROR(
                (outputReconPtr->n_filled_len + sampleTotalCount <= outputReconPtr->n_alloc_len),
                encode_context_ptr->app_callback_ptr,
                EB_ENC_ROB_OF_ERROR);

            // Initialize Y recon buffer
            picture_copy_kernel(
                reconReadPtr,
                recon_ptr->stride_y,
                reconWritePtr,
                recon_ptr->max_width - sequence_control_set_ptr->max_input_pad_right,
                recon_ptr->width - sequence_control_set_ptr->pad_right,
                recon_ptr->height - sequence_control_set_ptr->pad_bottom,
                1 << is16bit);

            outputReconPtr->n_filled_len += sampleTotalCount;

            // U Recon Samples
            sampleTotalCount = ((recon_ptr->max_width - sequence_control_set_ptr->max_input_pad_right) * (recon_ptr->max_height - sequence_control_set_ptr->max_input_pad_bottom) >> 2) << is16bit;
            reconReadPtr = recon_ptr->buffer_cb + ((recon_ptr->origin_y << is16bit) >> 1) * recon_ptr->stride_cb + ((recon_ptr->origin_x << is16bit) >> 1);
            reconWritePtr = &(outputReconPtr->p_buffer[outputReconPtr->n_filled_len]);

            CHECK_REPORT_ERROR(
                (outputReconPtr->n_filled_len + sampleTotalCount <= outputReconPtr->n_alloc_len),
                encode_context_ptr->app_callback_ptr,
                EB_ENC_ROB_OF_ERROR);

            // Initialize U recon buffer
            picture_copy_kernel(
                reconReadPtr,
                recon_ptr->stride_cb,
                reconWritePtr,
                (recon_ptr->max_width - sequence_control_set_ptr->max_input_pad_right) >> 1,
                (recon_ptr->width - sequence_control_set_ptr->pad_right) >> 1,
                (recon_ptr->height - sequence_control_set_ptr->pad_bottom) >> 1,
                1 << is16bit);
            outputReconPtr->n_filled_len += sampleTotalCount;

            // V Recon Samples
            sampleTotalCount = ((recon_ptr->max_width - sequence_control_set_ptr->max_input_pad_right) * (recon_ptr->max_height - sequence_control_set_ptr->max_input_pad_bottom) >> 2) << is16bit;
            reconReadPtr = recon_ptr->buffer_cr + ((recon_ptr->origin_y << is16bit) >> 1) * recon_ptr->stride_cr + ((recon_ptr->origin_x << is16bit) >> 1);
            reconWritePtr = &(outputReconPtr->p_buffer[outputReconPtr->n_filled_len]);

            CHECK_REPORT_ERROR(
                (outputReconPtr->n_filled_len + sampleTotalCount <= outputReconPtr->n_alloc_len),
                encode_context_ptr->app_callback_ptr,
                EB_ENC_ROB_OF_ERROR);

            // Initialize V recon buffer

            picture_copy_kernel(
                reconReadPtr,
                recon_ptr->stride_cr,
                reconWritePtr,
                (recon_ptr->max_width - sequence_control_set_ptr->max_input_pad_right) >> 1,
                (recon_ptr->width - sequence_control_set_ptr->pad_right) >> 1,
                (recon_ptr->height - sequence_control_set_ptr->pad_bottom) >> 1,
                1 << is16bit);
            outputReconPtr->n_filled_len += sampleTotalCount;
            outputReconPtr->pts = picture_control_set_ptr->picture_number;
        }

        // Post the Recon object
        eb_post_full_object(outputReconWrapperPtr);
    }
    else {
        // Overlay and altref have 1 recon only, which is from overlay pictures. So the recon of the alt_ref is not sent to the application.
        // However, to hanlde the end of sequence properly, total_number_of_recon_frames is increamented
        encode_context_ptr->total_number_of_recon_frames++;
    }
    eb_release_mutex(encode_context_ptr->total_number_of_recon_frame_mutex);
}

void psnr_calculations(
    PictureControlSet    *picture_control_set_ptr,
    SequenceControlSet   *sequence_control_set_ptr){
    EbBool is16bit = (sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    if (!is16bit) {
        EbPictureBufferDesc *recon_ptr;

        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_ptr = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
        else
            recon_ptr = picture_control_set_ptr->recon_picture_ptr;

        EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;

        uint64_t sseTotal[3] = { 0 };
        uint32_t   columnIndex;
        uint32_t   row_index = 0;
        uint64_t   residualDistortion = 0;
        EbByte  inputBuffer;
        EbByte  reconCoeffBuffer;

        EbByte buffer_y;
        EbByte buffer_cb;
        EbByte buffer_cr;

        // if current source picture was temporally filtered, use an alternative buffer which stores
        // the original source picture
        if(picture_control_set_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE){
            buffer_y = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[0];
            buffer_cb = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[1];
            buffer_cr = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[2];
        }
        else {
            buffer_y = input_picture_ptr->buffer_y;
            buffer_cb = input_picture_ptr->buffer_cb;
            buffer_cr = input_picture_ptr->buffer_cr;
        }

        reconCoeffBuffer = &((recon_ptr->buffer_y)[recon_ptr->origin_x + recon_ptr->origin_y * recon_ptr->stride_y]);
        inputBuffer = &(buffer_y[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y]);

        residualDistortion = 0;

        while (row_index < sequence_control_set_ptr->seq_header.max_frame_height) {
            columnIndex = 0;
            while (columnIndex < sequence_control_set_ptr->seq_header.max_frame_width) {
                residualDistortion += (int64_t)SQR((int64_t)(inputBuffer[columnIndex]) - (reconCoeffBuffer[columnIndex]));
                ++columnIndex;
            }

            inputBuffer += input_picture_ptr->stride_y;
            reconCoeffBuffer += recon_ptr->stride_y;
            ++row_index;
        }

        sseTotal[0] = residualDistortion;

        reconCoeffBuffer = &((recon_ptr->buffer_cb)[recon_ptr->origin_x / 2 + recon_ptr->origin_y / 2 * recon_ptr->stride_cb]);
        inputBuffer = &(buffer_cb[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);

        residualDistortion = 0;
        row_index = 0;
        while (row_index < sequence_control_set_ptr->chroma_height) {
            columnIndex = 0;
            while (columnIndex < sequence_control_set_ptr->chroma_width) {
                residualDistortion += (int64_t)SQR((int64_t)(inputBuffer[columnIndex]) - (reconCoeffBuffer[columnIndex]));
                ++columnIndex;
            }

            inputBuffer += input_picture_ptr->stride_cb;
            reconCoeffBuffer += recon_ptr->stride_cb;
            ++row_index;
        }

        sseTotal[1] = residualDistortion;

        reconCoeffBuffer = &((recon_ptr->buffer_cr)[recon_ptr->origin_x / 2 + recon_ptr->origin_y / 2 * recon_ptr->stride_cr]);
        inputBuffer = &(buffer_cr[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
        residualDistortion = 0;
        row_index = 0;

        while (row_index < sequence_control_set_ptr->chroma_height) {
            columnIndex = 0;
            while (columnIndex < sequence_control_set_ptr->chroma_width) {
                residualDistortion += (int64_t)SQR((int64_t)(inputBuffer[columnIndex]) - (reconCoeffBuffer[columnIndex]));
                ++columnIndex;
            }

            inputBuffer += input_picture_ptr->stride_cr;
            reconCoeffBuffer += recon_ptr->stride_cr;
            ++row_index;
        }

        sseTotal[2] = residualDistortion;
        picture_control_set_ptr->parent_pcs_ptr->luma_sse = (uint32_t)sseTotal[0];
        picture_control_set_ptr->parent_pcs_ptr->cb_sse = (uint32_t)sseTotal[1];
        picture_control_set_ptr->parent_pcs_ptr->cr_sse = (uint32_t)sseTotal[2];

        if(picture_control_set_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
            EB_FREE_ARRAY(buffer_y);
            EB_FREE_ARRAY(buffer_cb);
            EB_FREE_ARRAY(buffer_cr);
        }
    }
    else {
        EbPictureBufferDesc *recon_ptr;

        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_ptr = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            recon_ptr = picture_control_set_ptr->recon_picture16bit_ptr;
        EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;

        uint64_t sseTotal[3] = { 0 };
        uint32_t   columnIndex;
        uint32_t   row_index = 0;
        uint64_t   residualDistortion = 0;
        EbByte  inputBuffer;
        EbByte  inputBufferBitInc;
        uint16_t*  reconCoeffBuffer;

        if (sequence_control_set_ptr->static_config.ten_bit_format == 1) {
            const uint32_t luma_width = sequence_control_set_ptr->seq_header.max_frame_width;
            const uint32_t luma_height = sequence_control_set_ptr->seq_header.max_frame_height;
            const uint32_t chroma_width = sequence_control_set_ptr->chroma_width;
            const uint32_t picture_width_in_sb = (luma_width + 64 - 1) / 64;
            const uint32_t pictureHeighInLcu = (luma_height + 64 - 1) / 64;
            const uint32_t luma2BitWidth = luma_width / 4;
            const uint32_t chroma_height = luma_height / 2;
            const uint32_t chroma2BitWidth = luma_width / 8;
            uint32_t lcuNumberInHeight, lcuNumberInWidth;

            EbByte  inputBufferOrg = &((input_picture_ptr->buffer_y)[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
            uint16_t*  reconBufferOrg = (uint16_t*)(&((recon_ptr->buffer_y)[(recon_ptr->origin_x << is16bit) + (recon_ptr->origin_y << is16bit) * recon_ptr->stride_y]));;

            EbByte  inputBufferOrgU = &((input_picture_ptr->buffer_cb)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);;
            uint16_t*  reconBufferOrgU = reconCoeffBuffer = (uint16_t*)(&((recon_ptr->buffer_cb)[(recon_ptr->origin_x << is16bit) / 2 + (recon_ptr->origin_y << is16bit) / 2 * recon_ptr->stride_cb]));;

            EbByte  inputBufferOrgV = &((input_picture_ptr->buffer_cr)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);;
            uint16_t*  reconBufferOrgV = reconCoeffBuffer = (uint16_t*)(&((recon_ptr->buffer_cr)[(recon_ptr->origin_x << is16bit) / 2 + (recon_ptr->origin_y << is16bit) / 2 * recon_ptr->stride_cr]));;

            residualDistortion = 0;
            uint64_t residualDistortionU = 0;
            uint64_t residualDistortionV = 0;

            for (lcuNumberInHeight = 0; lcuNumberInHeight < pictureHeighInLcu; ++lcuNumberInHeight)
            {
                for (lcuNumberInWidth = 0; lcuNumberInWidth < picture_width_in_sb; ++lcuNumberInWidth)
                {
                    uint32_t tbOriginX = lcuNumberInWidth * 64;
                    uint32_t tbOriginY = lcuNumberInHeight * 64;
                    uint32_t sb_width = (luma_width - tbOriginX) < 64 ? (luma_width - tbOriginX) : 64;
                    uint32_t sb_height = (luma_height - tbOriginY) < 64 ? (luma_height - tbOriginY) : 64;

                    inputBuffer = inputBufferOrg + tbOriginY * input_picture_ptr->stride_y + tbOriginX;
                    inputBufferBitInc = input_picture_ptr->buffer_bit_inc_y + tbOriginY * luma2BitWidth + (tbOriginX / 4)*sb_height;
                    reconCoeffBuffer = reconBufferOrg + tbOriginY * recon_ptr->stride_y + tbOriginX;

                    uint64_t   j, k;
                    uint16_t   outPixel;
                    uint8_t    nBitPixel;
                    uint8_t   four2bitPels;
                    uint32_t     inn_stride = sb_width / 4;

                    for (j = 0; j < sb_height; j++)
                    {
                        for (k = 0; k < sb_width / 4; k++)
                        {
                            four2bitPels = inputBufferBitInc[k + j * inn_stride];

                            nBitPixel = (four2bitPels >> 6) & 3;
                            outPixel = inputBuffer[k * 4 + 0 + j * input_picture_ptr->stride_y] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortion += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 0 + j * recon_ptr->stride_y]);

                            nBitPixel = (four2bitPels >> 4) & 3;
                            outPixel = inputBuffer[k * 4 + 1 + j * input_picture_ptr->stride_y] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortion += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 1 + j * recon_ptr->stride_y]);

                            nBitPixel = (four2bitPels >> 2) & 3;
                            outPixel = inputBuffer[k * 4 + 2 + j * input_picture_ptr->stride_y] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortion += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 2 + j * recon_ptr->stride_y]);

                            nBitPixel = (four2bitPels >> 0) & 3;
                            outPixel = inputBuffer[k * 4 + 3 + j * input_picture_ptr->stride_y] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortion += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 3 + j * recon_ptr->stride_y]);
                        }
                    }

                    //U+V

                    tbOriginX = lcuNumberInWidth * 32;
                    tbOriginY = lcuNumberInHeight * 32;
                    sb_width = (chroma_width - tbOriginX) < 32 ? (chroma_width - tbOriginX) : 32;
                    sb_height = (chroma_height - tbOriginY) < 32 ? (chroma_height - tbOriginY) : 32;

                    inn_stride = sb_width / 4;

                    inputBuffer = inputBufferOrgU + tbOriginY * input_picture_ptr->stride_cb + tbOriginX;

                    inputBufferBitInc = input_picture_ptr->buffer_bit_inc_cb + tbOriginY * chroma2BitWidth + (tbOriginX / 4)*sb_height;

                    reconCoeffBuffer = reconBufferOrgU + tbOriginY * recon_ptr->stride_cb + tbOriginX;

                    for (j = 0; j < sb_height; j++)
                    {
                        for (k = 0; k < sb_width / 4; k++)
                        {
                            four2bitPels = inputBufferBitInc[k + j * inn_stride];

                            nBitPixel = (four2bitPels >> 6) & 3;
                            outPixel = inputBuffer[k * 4 + 0 + j * input_picture_ptr->stride_cb] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortionU += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 0 + j * recon_ptr->stride_cb]);

                            nBitPixel = (four2bitPels >> 4) & 3;
                            outPixel = inputBuffer[k * 4 + 1 + j * input_picture_ptr->stride_cb] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortionU += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 1 + j * recon_ptr->stride_cb]);

                            nBitPixel = (four2bitPels >> 2) & 3;
                            outPixel = inputBuffer[k * 4 + 2 + j * input_picture_ptr->stride_cb] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortionU += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 2 + j * recon_ptr->stride_cb]);

                            nBitPixel = (four2bitPels >> 0) & 3;
                            outPixel = inputBuffer[k * 4 + 3 + j * input_picture_ptr->stride_cb] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortionU += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 3 + j * recon_ptr->stride_cb]);
                        }
                    }

                    inputBuffer = inputBufferOrgV + tbOriginY * input_picture_ptr->stride_cr + tbOriginX;
                    inputBufferBitInc = input_picture_ptr->buffer_bit_inc_cr + tbOriginY * chroma2BitWidth + (tbOriginX / 4)*sb_height;
                    reconCoeffBuffer = reconBufferOrgV + tbOriginY * recon_ptr->stride_cr + tbOriginX;

                    for (j = 0; j < sb_height; j++)
                    {
                        for (k = 0; k < sb_width / 4; k++)
                        {
                            four2bitPels = inputBufferBitInc[k + j * inn_stride];

                            nBitPixel = (four2bitPels >> 6) & 3;
                            outPixel = inputBuffer[k * 4 + 0 + j * input_picture_ptr->stride_cr] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortionV += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 0 + j * recon_ptr->stride_cr]);

                            nBitPixel = (four2bitPels >> 4) & 3;
                            outPixel = inputBuffer[k * 4 + 1 + j * input_picture_ptr->stride_cr] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortionV += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 1 + j * recon_ptr->stride_cr]);

                            nBitPixel = (four2bitPels >> 2) & 3;
                            outPixel = inputBuffer[k * 4 + 2 + j * input_picture_ptr->stride_cr] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortionV += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 2 + j * recon_ptr->stride_cr]);

                            nBitPixel = (four2bitPels >> 0) & 3;
                            outPixel = inputBuffer[k * 4 + 3 + j * input_picture_ptr->stride_cr] << 2;
                            outPixel = outPixel | nBitPixel;
                            residualDistortionV += (int64_t)SQR((int64_t)outPixel - (int64_t)reconCoeffBuffer[k * 4 + 3 + j * recon_ptr->stride_cr]);
                        }
                    }
                }
            }

            sseTotal[0] = residualDistortion;
            sseTotal[1] = residualDistortionU;
            sseTotal[2] = residualDistortionV;
        }
        else {
            reconCoeffBuffer = (uint16_t*)(&((recon_ptr->buffer_y)[(recon_ptr->origin_x << is16bit) + (recon_ptr->origin_y << is16bit) * recon_ptr->stride_y]));

            // if current source picture was temporally filtered, use an alternative buffer which stores
            // the original source picture
            EbByte buffer_y, buffer_bit_inc_y;
            EbByte buffer_cb, buffer_bit_inc_cb;
            EbByte buffer_cr, buffer_bit_inc_cr;

            if(picture_control_set_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE){
                buffer_y = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[0];
                buffer_bit_inc_y = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[0];
                buffer_cb = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[1];
                buffer_bit_inc_cb = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[1];
                buffer_cr = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_ptr[2];
                buffer_bit_inc_cr = picture_control_set_ptr->parent_pcs_ptr->save_enhanced_picture_bit_inc_ptr[2];
            }else{
                buffer_y = input_picture_ptr->buffer_y;
                buffer_bit_inc_y = input_picture_ptr->buffer_bit_inc_y;
                buffer_cb = input_picture_ptr->buffer_cb;
                buffer_bit_inc_cb = input_picture_ptr->buffer_bit_inc_cb;
                buffer_cr = input_picture_ptr->buffer_cr;
                buffer_bit_inc_cr = input_picture_ptr->buffer_bit_inc_cr;
            }

            inputBuffer = &((buffer_y)[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
            inputBufferBitInc = &((buffer_bit_inc_y)[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_bit_inc_y]);

            residualDistortion = 0;

            while (row_index < sequence_control_set_ptr->seq_header.max_frame_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->seq_header.max_frame_width) {
                    residualDistortion += (int64_t)SQR((int64_t)((((inputBuffer[columnIndex]) << 2) | ((inputBufferBitInc[columnIndex] >> 6) & 3))) - (reconCoeffBuffer[columnIndex]));

                    ++columnIndex;
                }

                inputBuffer += input_picture_ptr->stride_y;
                inputBufferBitInc += input_picture_ptr->stride_bit_inc_y;
                reconCoeffBuffer += recon_ptr->stride_y;
                ++row_index;
            }

            sseTotal[0] = residualDistortion;

            reconCoeffBuffer = (uint16_t*)(&((recon_ptr->buffer_cb)[(recon_ptr->origin_x << is16bit) / 2 + (recon_ptr->origin_y << is16bit) / 2 * recon_ptr->stride_cb]));
            inputBuffer = &((buffer_cb)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
            inputBufferBitInc = &((buffer_bit_inc_cb)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_bit_inc_cb]);

            residualDistortion = 0;
            row_index = 0;
            while (row_index < sequence_control_set_ptr->chroma_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->chroma_width) {
                    residualDistortion += (int64_t)SQR((int64_t)((((inputBuffer[columnIndex]) << 2) | ((inputBufferBitInc[columnIndex] >> 6) & 3))) - (reconCoeffBuffer[columnIndex]));
                    ++columnIndex;
                }

                inputBuffer += input_picture_ptr->stride_cb;
                inputBufferBitInc += input_picture_ptr->stride_bit_inc_cb;
                reconCoeffBuffer += recon_ptr->stride_cb;
                ++row_index;
            }

            sseTotal[1] = residualDistortion;

            reconCoeffBuffer = (uint16_t*)(&((recon_ptr->buffer_cr)[(recon_ptr->origin_x << is16bit) / 2 + (recon_ptr->origin_y << is16bit) / 2 * recon_ptr->stride_cr]));
            inputBuffer = &((buffer_cr)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            inputBufferBitInc = &((buffer_bit_inc_cr)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_bit_inc_cr]);

            residualDistortion = 0;
            row_index = 0;

            while (row_index < sequence_control_set_ptr->chroma_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->chroma_width) {
                    residualDistortion += (int64_t)SQR((int64_t)((((inputBuffer[columnIndex]) << 2) | ((inputBufferBitInc[columnIndex] >> 6) & 3))) - (reconCoeffBuffer[columnIndex]));
                    ++columnIndex;
                }

                inputBuffer += input_picture_ptr->stride_cr;
                inputBufferBitInc += input_picture_ptr->stride_bit_inc_cr;
                reconCoeffBuffer += recon_ptr->stride_cr;
                ++row_index;
            }

            sseTotal[2] = residualDistortion;

            if(picture_control_set_ptr->parent_pcs_ptr->temporal_filtering_on == EB_TRUE) {
                EB_FREE_ARRAY(buffer_y);
                EB_FREE_ARRAY(buffer_cb);
                EB_FREE_ARRAY(buffer_cr);
                EB_FREE_ARRAY(buffer_bit_inc_y);
                EB_FREE_ARRAY(buffer_bit_inc_cb);
                EB_FREE_ARRAY(buffer_bit_inc_cr);
            }
        }

        picture_control_set_ptr->parent_pcs_ptr->luma_sse = (uint32_t)sseTotal[0];
        picture_control_set_ptr->parent_pcs_ptr->cb_sse = (uint32_t)sseTotal[1];
        picture_control_set_ptr->parent_pcs_ptr->cr_sse = (uint32_t)sseTotal[2];
    }
}

void pad_ref_and_set_flags(
    PictureControlSet    *picture_control_set_ptr,
    SequenceControlSet   *sequence_control_set_ptr){

    EbReferenceObject   *referenceObject = (EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
    EbPictureBufferDesc *refPicPtr = (EbPictureBufferDesc*)referenceObject->reference_picture;
    EbPictureBufferDesc *refPic16BitPtr = (EbPictureBufferDesc*)referenceObject->reference_picture16bit;
    EbBool                is16bit = (sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    if (!is16bit) {
        // Y samples
        generate_padding(
            refPicPtr->buffer_y,
            refPicPtr->stride_y,
            refPicPtr->width,
            refPicPtr->height,
            refPicPtr->origin_x,
            refPicPtr->origin_y);

        // Cb samples
        generate_padding(
            refPicPtr->buffer_cb,
            refPicPtr->stride_cb,
            refPicPtr->width >> 1,
            refPicPtr->height >> 1,
            refPicPtr->origin_x >> 1,
            refPicPtr->origin_y >> 1);

        // Cr samples
        generate_padding(
            refPicPtr->buffer_cr,
            refPicPtr->stride_cr,
            refPicPtr->width >> 1,
            refPicPtr->height >> 1,
            refPicPtr->origin_x >> 1,
            refPicPtr->origin_y >> 1);
    }

    //We need this for MCP
    if (is16bit) {
        // Y samples
        generate_padding16_bit(
            refPic16BitPtr->buffer_y,
            refPic16BitPtr->stride_y << 1,
            refPic16BitPtr->width << 1,
            refPic16BitPtr->height,
            refPic16BitPtr->origin_x << 1,
            refPic16BitPtr->origin_y);

        // Cb samples
        generate_padding16_bit(
            refPic16BitPtr->buffer_cb,
            refPic16BitPtr->stride_cb << 1,
            refPic16BitPtr->width,
            refPic16BitPtr->height >> 1,
            refPic16BitPtr->origin_x,
            refPic16BitPtr->origin_y >> 1);

        // Cr samples
        generate_padding16_bit(
            refPic16BitPtr->buffer_cr,
            refPic16BitPtr->stride_cr << 1,
            refPic16BitPtr->width,
            refPic16BitPtr->height >> 1,
            refPic16BitPtr->origin_x,
            refPic16BitPtr->origin_y >> 1);

        // Hsan: unpack ref samples (to be used @ MD)
        un_pack2d(
            (uint16_t*) refPic16BitPtr->buffer_y,
            refPic16BitPtr->stride_y,
            refPicPtr->buffer_y,
            refPicPtr->stride_y,
            refPicPtr->buffer_bit_inc_y,
            refPicPtr->stride_bit_inc_y,
            refPic16BitPtr->width  + (refPicPtr->origin_x << 1),
            refPic16BitPtr->height + (refPicPtr->origin_y << 1));

        un_pack2d(
            (uint16_t*)refPic16BitPtr->buffer_cb,
            refPic16BitPtr->stride_cb,
            refPicPtr->buffer_cb,
            refPicPtr->stride_cb,
            refPicPtr->buffer_bit_inc_cb,
            refPicPtr->stride_bit_inc_cb,
            (refPic16BitPtr->width + (refPicPtr->origin_x << 1)) >> 1,
            (refPic16BitPtr->height + (refPicPtr->origin_y << 1)) >> 1);

        un_pack2d(
            (uint16_t*)refPic16BitPtr->buffer_cr,
            refPic16BitPtr->stride_cr,
            refPicPtr->buffer_cr,
            refPicPtr->stride_cr,
            refPicPtr->buffer_bit_inc_cr,
            refPicPtr->stride_bit_inc_cr,
            (refPic16BitPtr->width + (refPicPtr->origin_x << 1)) >> 1,
            (refPic16BitPtr->height + (refPicPtr->origin_y << 1)) >> 1);
    }
    // set up the ref POC
    referenceObject->ref_poc = picture_control_set_ptr->parent_pcs_ptr->picture_number;

    // set up the QP
    referenceObject->qp = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->picture_qp;

    // set up the Slice Type
    referenceObject->slice_type = picture_control_set_ptr->parent_pcs_ptr->slice_type;
#if TWO_PASS
    referenceObject->referenced_area_avg = picture_control_set_ptr->parent_pcs_ptr->referenced_area_avg;
#endif
}

void copy_statistics_to_ref_object(
    PictureControlSet    *picture_control_set_ptr,
    SequenceControlSet   *sequence_control_set_ptr
)
{
    picture_control_set_ptr->intra_coded_area = (100 * picture_control_set_ptr->intra_coded_area) / (sequence_control_set_ptr->seq_header.max_frame_width * sequence_control_set_ptr->seq_header.max_frame_height);
    if (picture_control_set_ptr->slice_type == I_SLICE)
        picture_control_set_ptr->intra_coded_area = 0;

    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->intra_coded_area = (uint8_t)(picture_control_set_ptr->intra_coded_area);

    uint32_t sb_index;
    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index)
        ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->non_moving_index_array[sb_index] = picture_control_set_ptr->parent_pcs_ptr->non_moving_index_array[sb_index];
    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->tmp_layer_idx = (uint8_t)picture_control_set_ptr->temporal_layer_index;
    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->is_scene_change = picture_control_set_ptr->parent_pcs_ptr->scene_change_flag;

    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->cdef_frame_strength = picture_control_set_ptr->parent_pcs_ptr->cdef_frame_strength;

    Av1Common* cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;
    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->sg_frame_ep = cm->sg_frame_ep;
    if (sequence_control_set_ptr->mfmv_enabled) {
        ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->frame_type = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.frame_type;
        ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->order_hint = picture_control_set_ptr->parent_pcs_ptr->cur_order_hint;
        memcpy(((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->ref_order_hint, picture_control_set_ptr->parent_pcs_ptr->ref_order_hint, 7 * sizeof(uint32_t));
    }
}

/******************************************************
* Derive EncDec Settings for OQ
Input   : encoder mode and pd pass
Output  : EncDec Kernel signal(s)
******************************************************/
EbErrorType signal_derivation_enc_dec_kernel_oq(
    SequenceControlSet    *sequence_control_set_ptr,
    PictureControlSet     *picture_control_set_ptr,
    ModeDecisionContext   *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

#if MULTI_PASS_PD
    // Tx_search Level                                Settings
    // 0                                              OFF
    // 1                                              Tx search at encdec
    // 2                                              Tx search at inter-depth
    // 3                                              Tx search at full loop
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->tx_search_level = TX_SEARCH_OFF;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
    else if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
        if (picture_control_set_ptr->enc_mode <= ENC_M6)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
            else
                context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    else
        if (picture_control_set_ptr->enc_mode <= ENC_M4)
            context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else if (picture_control_set_ptr->enc_mode <= ENC_M7) {
            if (picture_control_set_ptr->temporal_layer_index == 0)
                context_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
            else
                context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
        }
        else
            context_ptr->tx_search_level = TX_SEARCH_ENC_DEC;

    // Set tx search skip weights (MAX_MODE_COST: no skipping; 0: always skipping)
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->tx_weight = MAX_MODE_COST;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
    else
#if PRESETS_TUNE
        if (MR_MODE) // tx weight
            context_ptr->tx_weight = MAX_MODE_COST;
        else {
            if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
                context_ptr->tx_weight = MAX_MODE_COST;
            else if (!MR_MODE && picture_control_set_ptr->enc_mode <= ENC_M5)
                context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
            else if (!MR_MODE) {
                if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                    context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
                else
                    context_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
            }
        }
#else
        if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
            context_ptr->tx_weight = MAX_MODE_COST;
        else if (!MR_MODE && picture_control_set_ptr->enc_mode <= ENC_M1)
            context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else if (!MR_MODE) {
            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                context_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
            else
                context_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
        }
#endif

    // Set tx search reduced set falg (0: full tx set; 1: reduced tx set; 1: two tx))
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->tx_search_reduced_set = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->tx_search_reduced_set = 0;
    else if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
#if PRESETS_TUNE
        if (picture_control_set_ptr->enc_mode <= ENC_M5)
            context_ptr->tx_search_reduced_set = 0;
        else if (picture_control_set_ptr->enc_mode <= ENC_M6)
            if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
                context_ptr->tx_search_reduced_set = 0;
            else
                context_ptr->tx_search_reduced_set = 1;
        else if (picture_control_set_ptr->enc_mode <= ENC_M7)
            context_ptr->tx_search_reduced_set = 1;
        else
            context_ptr->tx_search_reduced_set = 2;
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
            context_ptr->tx_search_reduced_set = 0;
        else if (picture_control_set_ptr->enc_mode <= ENC_M6)
            if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
                context_ptr->tx_search_reduced_set = 0;
            else
                context_ptr->tx_search_reduced_set = 1;
        else if (picture_control_set_ptr->enc_mode <= ENC_M7)
            context_ptr->tx_search_reduced_set = 1;
        else
            context_ptr->tx_search_reduced_set = 2;
#endif
    else
        if (context_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
            context_ptr->tx_search_reduced_set = 0;
#if PRESETS_TUNE
        else if (picture_control_set_ptr->enc_mode <= ENC_M2)
            context_ptr->tx_search_reduced_set = 0;
        else
            context_ptr->tx_search_reduced_set = 1;
#else
        else if (picture_control_set_ptr->enc_mode <= ENC_M1)
            context_ptr->tx_search_reduced_set = 0;
        else if (picture_control_set_ptr->enc_mode <= ENC_M3)
            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                context_ptr->tx_search_reduced_set = 0;
            else
                context_ptr->tx_search_reduced_set = 1;
        else
            context_ptr->tx_search_reduced_set = 1;
#endif
#endif
    // Interpolation search Level                     Settings
    // 0                                              OFF
    // 1                                              Interpolation search at inter-depth
    // 2                                              Interpolation search at full loop
    // 3                                              Chroma blind interpolation search at fast loop
    // 4                                              Interpolation search at fast loop
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
    else if (context_ptr->pd_pass == PD_PASS_1) {
        if (picture_control_set_ptr->temporal_layer_index == 0)
            context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
        else
            context_ptr->interpolation_search_level = IT_SEARCH_OFF;
    }
    else
    if (MR_MODE)
        context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP;
    else if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#if PRESETS_TUNE
    else if (picture_control_set_ptr->enc_mode <= ENC_M2)
#else
    else if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
        context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
#if !PRESETS_TUNE
    else if (picture_control_set_ptr->enc_mode <= ENC_M3)
        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
            context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
        else
            context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#endif
    else if (picture_control_set_ptr->enc_mode <= ENC_M7)
        if (picture_control_set_ptr->temporal_layer_index == 0)
            context_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
        else
            context_ptr->interpolation_search_level = IT_SEARCH_OFF;
    else
        context_ptr->interpolation_search_level = IT_SEARCH_OFF;
#endif
    // Set Chroma Mode
    // Level                Settings
    // CHROMA_MODE_0  0     Full chroma search @ MD
    // CHROMA_MODE_1  1     Fast chroma search @ MD
    // CHROMA_MODE_2  2     Chroma blind @ MD + CFL @ EP
    // CHROMA_MODE_3  3     Chroma blind @ MD + no CFL @ EP
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->chroma_level = CHROMA_MODE_2; // or CHROMA_MODE_3
    else if (context_ptr->pd_pass == PD_PASS_1) {
        if (picture_control_set_ptr->temporal_layer_index == 0)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else
            context_ptr->chroma_level = CHROMA_MODE_1;
    }
    else
#endif
    if (sequence_control_set_ptr->static_config.set_chroma_mode == DEFAULT) {
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            if (picture_control_set_ptr->enc_mode <= ENC_M6)
                context_ptr->chroma_level = CHROMA_MODE_1;
            else
                if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)
                    context_ptr->chroma_level = CHROMA_MODE_1;
                else
                    context_ptr->chroma_level = (sequence_control_set_ptr->encoder_bit_depth == EB_8BIT) ?
                    CHROMA_MODE_2 :
                    CHROMA_MODE_3;
        else
#if PRESETS_TUNE
#if ENHANCED_M0_SETTINGS
            if (picture_control_set_ptr->enc_mode == ENC_M0)
                context_ptr->chroma_level = CHROMA_MODE_0;
#else
            if (MR_MODE)
                context_ptr->chroma_level = CHROMA_MODE_0;
#endif
            else
                if (picture_control_set_ptr->enc_mode <= ENC_M5 && picture_control_set_ptr->temporal_layer_index == 0)
                    context_ptr->chroma_level = CHROMA_MODE_0;
                else
                    if (picture_control_set_ptr->enc_mode <= ENC_M5)
                        context_ptr->chroma_level = CHROMA_MODE_1;
                    else
                        context_ptr->chroma_level = (sequence_control_set_ptr->encoder_bit_depth == EB_8BIT) ?
                        CHROMA_MODE_2 :
                        CHROMA_MODE_3;

#else
#if ENHANCED_M0_SETTINGS
            if (picture_control_set_ptr->enc_mode == ENC_M0)
                context_ptr->chroma_level = CHROMA_MODE_0;
#else
        if (MR_MODE)
            context_ptr->chroma_level = CHROMA_MODE_0;
        else
        if (picture_control_set_ptr->enc_mode == ENC_M0 && picture_control_set_ptr->temporal_layer_index == 0)
            context_ptr->chroma_level = CHROMA_MODE_0;
#endif
        else
        if (picture_control_set_ptr->enc_mode <= ENC_M4)
            context_ptr->chroma_level = CHROMA_MODE_1;
        else
            context_ptr->chroma_level = (sequence_control_set_ptr->encoder_bit_depth == EB_8BIT) ?
                CHROMA_MODE_2 :
                CHROMA_MODE_3 ;
#endif
    }
    else // use specified level
      context_ptr->chroma_level = sequence_control_set_ptr->static_config.set_chroma_mode;

    // Set the full loop escape level
    // Level                Settings
    // 0                    Off
    // 1                    On but only INTRA
    // 2                    On both INTRA and INTER
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->full_loop_escape = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->full_loop_escape = 0;
    else
#endif
    if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
#if PRESETS_TUNE
        if (picture_control_set_ptr->enc_mode <= ENC_M5)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            context_ptr->full_loop_escape = 0;
        else
            context_ptr->full_loop_escape = 2;
    else
    if (picture_control_set_ptr->enc_mode <= ENC_M5)
        context_ptr->full_loop_escape = 0;
    else
        context_ptr->full_loop_escape = 2;

    // Set global MV injection
    // Level                Settings
    // 0                    Injection off (Hsan: but not derivation as used by MV ref derivation)
    // 1                    On
#if GLOBAL_WARPED_MOTION
    if (sequence_control_set_ptr->static_config.enable_global_motion == EB_TRUE)
    {
#if MULTI_PASS_PD
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->global_mv_injection = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->global_mv_injection = 0;
        else
#endif
#if PRESETS_OPT
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#else
        if (picture_control_set_ptr->enc_mode == ENC_M0)
#endif
            context_ptr->global_mv_injection = 1;
        else
            context_ptr->global_mv_injection = 0;
    }
    else
        context_ptr->global_mv_injection = 0;
#else
    if (sequence_control_set_ptr->static_config.enable_global_warped_motion == EB_TRUE)
    {
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
        {
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
                context_ptr->global_mv_injection = 1;
            else
                context_ptr->global_mv_injection = 0;
        }
        else
        {
            if (picture_control_set_ptr->enc_mode <= ENC_M7)
                context_ptr->global_mv_injection = 1;
            else
                context_ptr->global_mv_injection = 0;
        }
    }
    else
        context_ptr->global_mv_injection = 0;
#endif
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->new_nearest_injection = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->new_nearest_injection = 0;
    else
        context_ptr->new_nearest_injection = 1;
#endif

#if MULTI_PASS_PD // Shut nx4 and 4xn if 1st pass
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->new_nearest_near_comb_injection = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->new_nearest_near_comb_injection = 0;
    else
#endif
    if (sequence_control_set_ptr->static_config.new_nearest_comb_inject == DEFAULT)
#if M0_OPT
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            context_ptr->new_nearest_near_comb_injection = 0;
        else
#endif
            if (picture_control_set_ptr->enc_mode <= ENC_M0 && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
                context_ptr->new_nearest_near_comb_injection = 1;
            else
                context_ptr->new_nearest_near_comb_injection = 0;
    else
        context_ptr->new_nearest_near_comb_injection = sequence_control_set_ptr->static_config.new_nearest_comb_inject;

#if MULTI_PASS_PD // Shut nx4 and 4xn if 1st pass
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->nx4_4xn_parent_mv_injection = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->nx4_4xn_parent_mv_injection = 0;
    else
#endif
    if (sequence_control_set_ptr->static_config.nx4_4xn_parent_mv_inject == DEFAULT)
#if PRESETS_TUNE
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            if (picture_control_set_ptr->enc_mode <= ENC_M2)
                context_ptr->nx4_4xn_parent_mv_injection = 1;
            else
                context_ptr->nx4_4xn_parent_mv_injection = 0;
        else
#if PRESETS_OPT
            if (picture_control_set_ptr->enc_mode <= ENC_M7)
#else
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
                context_ptr->nx4_4xn_parent_mv_injection = 1;
            else
                context_ptr->nx4_4xn_parent_mv_injection = 0;
#else
        if (picture_control_set_ptr->enc_mode == ENC_M0)
            context_ptr->nx4_4xn_parent_mv_injection = 1;
        else
            context_ptr->nx4_4xn_parent_mv_injection = 0;
#endif
    else
        context_ptr->nx4_4xn_parent_mv_injection = sequence_control_set_ptr->static_config.nx4_4xn_parent_mv_inject;

    // Set warped motion injection
    // Level                Settings
    // 0                    OFF
    // 1                    On
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->warped_motion_injection = 0;
    }
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->warped_motion_injection = 1;
    }
    else
#endif
    if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
        context_ptr->warped_motion_injection = 0;
    else
        context_ptr->warped_motion_injection = 1;

    // Set unipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->unipred3x3_injection = 0;
    }
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->unipred3x3_injection = 2;
    }
    else
#endif
    if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
#if PRESETS_TUNE
        if (picture_control_set_ptr->enc_mode <= ENC_M2)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
#if M0_OPT
            context_ptr->unipred3x3_injection = 2;
#else
            context_ptr->unipred3x3_injection = 1;
#endif
        else
            context_ptr->unipred3x3_injection = 0;
    else
#if PRESETS_OPT
    if (picture_control_set_ptr->enc_mode <= ENC_M7)
            context_ptr->unipred3x3_injection = 1;
#else
    if (picture_control_set_ptr->enc_mode <= ENC_M1)
        context_ptr->unipred3x3_injection = 1;
    else if (picture_control_set_ptr->enc_mode <= ENC_M4)
        context_ptr->unipred3x3_injection = 2;
#endif
    else
        context_ptr->unipred3x3_injection = 0;

    // Set bipred3x3 injection
    // Level                Settings
    // 0                    OFF
    // 1                    ON FULL
    // 2                    Reduced set
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->bipred3x3_injection = 0;
    }
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->bipred3x3_injection = 2;
    }
    else
#endif
    if (sequence_control_set_ptr->static_config.bipred_3x3_inject == DEFAULT)
#if PRESETS_TUNE
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            if (picture_control_set_ptr->enc_mode <= ENC_M4)
                context_ptr->bipred3x3_injection = 1;
            else
                context_ptr->bipred3x3_injection = 0;
        else
            if (picture_control_set_ptr->enc_mode <= ENC_M2)
                context_ptr->bipred3x3_injection = 1;
            else if (picture_control_set_ptr->enc_mode <= ENC_M4)
                context_ptr->bipred3x3_injection = 2;
            else
                context_ptr->bipred3x3_injection = 0;
#else
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
                context_ptr->bipred3x3_injection = 1;
            else
                context_ptr->bipred3x3_injection = 0;
        else
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
                context_ptr->bipred3x3_injection = 1;
            else if (picture_control_set_ptr->enc_mode <= ENC_M4)
                context_ptr->bipred3x3_injection = 2;
            else
                context_ptr->bipred3x3_injection = 0;
#endif
    else
        context_ptr->bipred3x3_injection = sequence_control_set_ptr->static_config.bipred_3x3_inject;

    // Level                Settings
    // 0                    Level 0: OFF
    // 1                    Level 1: 7x5 full-pel search + sub-pel refinement off
    // 2                    Level 2: 7x5 full-pel search +  (H + V) sub-pel refinement only = 4 half-pel + 4 quarter-pel = 8 positions + pred_me_distortion to pa_me_distortion deviation on
    // 3                    Level 3: 7x5 full-pel search +  (H + V + D only ~ the best) sub-pel refinement = up to 6 half-pel + up to 6  quarter-pel = up to 12 positions + pred_me_distortion to pa_me_distortion deviation on
    // 4                    Level 4: 7x5 full-pel search +  (H + V + D) sub-pel refinement = 8 half-pel + 8 quarter-pel = 16 positions + pred_me_distortion to pa_me_distortion deviation on
    // 5                    Level 5: 7x5 full-pel search +  (H + V + D) sub-pel refinement = 8 half-pel + 8 quarter-pel = 16 positions + pred_me_distortion to pa_me_distortion deviation off
    if (picture_control_set_ptr->slice_type != I_SLICE) {
#if MULTI_PASS_PD // Shut predictive if 1st pass
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->predictive_me_level = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->predictive_me_level = 2;
    else
#endif
        if (sequence_control_set_ptr->static_config.pred_me == DEFAULT) {
#if PRESETS_TUNE
            if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
                if (picture_control_set_ptr->enc_mode <= ENC_M1)
#if M0_OPT
                    context_ptr->predictive_me_level = (picture_control_set_ptr->enc_mode <= ENC_M0) ? 2 : 4;
#else
                    context_ptr->predictive_me_level = 4;
#endif
                else if (picture_control_set_ptr->enc_mode <= ENC_M4)
                    context_ptr->predictive_me_level = 2;
                else
                    context_ptr->predictive_me_level = 0;
            else
#if M0_OPT
                if (picture_control_set_ptr->enc_mode <= ENC_M2)
#if PRESETS_OPT
                    context_ptr->predictive_me_level = 5;
#else
                    context_ptr->predictive_me_level = (picture_control_set_ptr->enc_mode <= ENC_M0) ? 5 : 4;
#endif
#else
                if (picture_control_set_ptr->enc_mode <= ENC_M1)
                    context_ptr->predictive_me_level = 4;
#endif
                else if (picture_control_set_ptr->enc_mode <= ENC_M4)
                    context_ptr->predictive_me_level = 2;
                else
                    context_ptr->predictive_me_level = 0;
#else
#if M0_OPT
            if (picture_control_set_ptr->enc_mode <= ENC_M0)
                context_ptr->predictive_me_level = picture_control_set_ptr->parent_pcs_ptr->sc_content_detected ? 2 : 5;
            else if (picture_control_set_ptr->enc_mode <= ENC_M1)
#else
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
                context_ptr->predictive_me_level = 4;
            else if (picture_control_set_ptr->enc_mode <= ENC_M4)
                context_ptr->predictive_me_level = 2;
            else
                context_ptr->predictive_me_level = 0;
#endif
        } else
            context_ptr->predictive_me_level = sequence_control_set_ptr->static_config.pred_me;
    } else
        context_ptr->predictive_me_level = 0;


    // Derive md_staging_mode
    //
    // MD_STAGING_MODE_0
    // Default Parameters
    //
    // MD_STAGING_MODE_1
    //  __________________________________________________________________________________________________________________
    // |        | md_stage_0                  | md_stage_1                     | md_stage_2                              |
    // |________|_____________________________|________________________________|_________________________________________|
    // |CLASS_0 |Prediction for Luma & Chroma |T, Q, Q-1, T-1 for Luma Only    |T, Q, Q-1, T-1 or Luma & Chroma          |
    // |CLASS_6 |                             |No RDOQ                         |RDOQ                                     |
    // |        |                             |No Tx Type Search               |Tx Type Search                           |
    // |        |                             |No Tx Size Search               |Tx Size Search                           |
    // |        |                             |                                |CFL vs. Independent                      |
    // |________|_____________________________|________________________________|_________________________________________|
    // |CLASS_1 |Prediction for Luma Only     |T, Q, Q-1, T-1 for Luma Only    |T, Q, Q-1, T-1 for Luma & Chroma         |
    // |CLASS_2 |No Interpolation Search      |No RDOQ                         |RDOQ                                     |
    // |CLASS_3 |Bilinear Interpolation       |No Tx Type Search               |Tx Type Search                           |
    // |CLASS_4 |                             |No Tx Size Search               |Tx Size Search                           |
    // |CLASS_5 |                             |Interpolation Search            |                                         |
    // |________|_____________________________|________________________________|_________________________________________|
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;
    }
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    }
    else
#endif
    if (picture_control_set_ptr->enc_mode <= ENC_M4)
        context_ptr->md_staging_mode = MD_STAGING_MODE_1;
    else
        context_ptr->md_staging_mode = MD_STAGING_MODE_0;

#if MULTI_PASS_PD
    // Set md staging count level
    // Level 0              minimum count = 1
    // Level 1              set towards the best possible partitioning (to further optimize)
    // Level 2              HG: breack down or look up-table(s) are required !
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->md_staging_count_level = 0;
    }
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->md_staging_count_level = 1;
    }
    else {
        context_ptr->md_staging_count_level = 2;
    }
#endif

    // Combine MD Class1&2
    // 0                    OFF
    // 1                    ON
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->combine_class12 = 0;
    } else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->combine_class12 = 0;
    }
    else
#endif
    if (sequence_control_set_ptr->static_config.combine_class_12 == DEFAULT)
#if PRESETS_TUNE
        context_ptr->combine_class12 = (picture_control_set_ptr->enc_mode <= ENC_M4) ? 0 : 1;
#else
        context_ptr->combine_class12 = (picture_control_set_ptr->enc_mode == ENC_M0) ? 0 : 1;
#endif
    else
        context_ptr->combine_class12 = sequence_control_set_ptr->static_config.combine_class_12;

    // Set interpolation filter search blk size
    // Level                Settings
    // 0                    ON for 8x8 and above
    // 1                    ON for 16x16 and above
    // 2                    ON for 32x32 and above
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->interpolation_filter_search_blk_size = 0;
    }
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->interpolation_filter_search_blk_size = 0;
    }
    else
#endif
    if (picture_control_set_ptr->enc_mode <= ENC_M4)
        context_ptr->interpolation_filter_search_blk_size = 0;
    else
        context_ptr->interpolation_filter_search_blk_size = 1;

    // Set PF MD
    context_ptr->pf_md_mode = PF_OFF;
    // Derive Spatial SSE Flag
#if MULTI_PASS_PD // Shut spatial SSE @ full loop
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->spatial_sse_full_loop = EB_TRUE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->spatial_sse_full_loop = EB_TRUE;
    else
#endif
    if (sequence_control_set_ptr->static_config.spatial_sse_fl == DEFAULT)
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            if (picture_control_set_ptr->enc_mode <= ENC_M6)
                context_ptr->spatial_sse_full_loop = EB_TRUE;
            else
                context_ptr->spatial_sse_full_loop = EB_FALSE;
        else
            if (picture_control_set_ptr->enc_mode <= ENC_M4)
                context_ptr->spatial_sse_full_loop = EB_TRUE;
            else
                context_ptr->spatial_sse_full_loop = EB_FALSE;
    else
        context_ptr->spatial_sse_full_loop = sequence_control_set_ptr->static_config.spatial_sse_fl;

    if (context_ptr->chroma_level <= CHROMA_MODE_1)
        context_ptr->blk_skip_decision = EB_TRUE;
    else
        context_ptr->blk_skip_decision = EB_FALSE;

    // Derive Trellis Quant Coeff Optimization Flag
#if MULTI_PASS_PD // Shut spatial SSE @ full loop
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->trellis_quant_coeff_optimization = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->trellis_quant_coeff_optimization = EB_FALSE;
    else
#endif
    if (sequence_control_set_ptr->static_config.enable_trellis == DEFAULT)
#if PRESETS_TUNE
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            if (picture_control_set_ptr->enc_mode <= ENC_M2)
                context_ptr->trellis_quant_coeff_optimization = EB_TRUE;
            else
                context_ptr->trellis_quant_coeff_optimization = EB_FALSE;
        else
            if (picture_control_set_ptr->enc_mode <= ENC_M2)
                context_ptr->trellis_quant_coeff_optimization = EB_TRUE;
            else
                context_ptr->trellis_quant_coeff_optimization = EB_FALSE;
#else
        if (picture_control_set_ptr->enc_mode == ENC_M0)
            context_ptr->trellis_quant_coeff_optimization = EB_TRUE;
        else
            context_ptr->trellis_quant_coeff_optimization = EB_FALSE;
#endif
    else
        context_ptr->trellis_quant_coeff_optimization = sequence_control_set_ptr->static_config.enable_trellis;

    // Derive redundant block
#if MULTI_PASS_PD // Shut redundant_blk
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->redundant_blk = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->redundant_blk = EB_TRUE;
    else
#endif
    if (sequence_control_set_ptr->static_config.enable_redundant_blk == DEFAULT)
#if PRESETS_TUNE
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            if (picture_control_set_ptr->enc_mode <= ENC_M8)
                context_ptr->redundant_blk = EB_TRUE;
            else
                context_ptr->redundant_blk = EB_FALSE;
        else
            if (picture_control_set_ptr->enc_mode <= ENC_M5)
                context_ptr->redundant_blk = EB_TRUE;
            else
                context_ptr->redundant_blk = EB_FALSE;
#else
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
                context_ptr->redundant_blk = EB_TRUE;
            else
                context_ptr->redundant_blk = EB_FALSE;
        else
            if (picture_control_set_ptr->enc_mode <= ENC_M5)
                context_ptr->redundant_blk = EB_TRUE;
            else
                context_ptr->redundant_blk = EB_FALSE;
#endif
    else
        context_ptr->redundant_blk = sequence_control_set_ptr->static_config.enable_redundant_blk;

    if (sequence_control_set_ptr->static_config.encoder_bit_depth == EB_8BIT)
#if MULTI_PASS_PD
        if (context_ptr->pd_pass == PD_PASS_0)
            context_ptr->edge_based_skip_angle_intra = 0;
        else if (context_ptr->pd_pass == PD_PASS_1)
            context_ptr->edge_based_skip_angle_intra = 1;
        else
#endif
        if (sequence_control_set_ptr->static_config.edge_skp_angle_intra == DEFAULT) {
        if (MR_MODE)
            context_ptr->edge_based_skip_angle_intra = 0;
#if M0_OPT
#if PRESETS_TUNE
#if PRESETS_OPT
        else if (picture_control_set_ptr->enc_mode <= ENC_M7 && !picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
#else
        else if (picture_control_set_ptr->enc_mode <= ENC_M1 && !picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
#endif
            context_ptr->edge_based_skip_angle_intra = 0;
#else
        else if (picture_control_set_ptr->enc_mode <= ENC_M0)
            context_ptr->edge_based_skip_angle_intra = picture_control_set_ptr->parent_pcs_ptr->sc_content_detected ? 1 : 0;
#endif
        else
            context_ptr->edge_based_skip_angle_intra = 1;
#else
        else
            context_ptr->edge_based_skip_angle_intra = (picture_control_set_ptr->enc_mode == ENC_M0 && picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0) ? 0 : 1;
#endif
        } else
            context_ptr->edge_based_skip_angle_intra = sequence_control_set_ptr->static_config.edge_skp_angle_intra;
    else
        context_ptr->edge_based_skip_angle_intra = 0;

#if MULTI_PASS_PD // Shut spatial SSE @ full loop
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->prune_ref_frame_for_rec_partitions = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->prune_ref_frame_for_rec_partitions = 1;
    else
#endif
    if (sequence_control_set_ptr->static_config.prune_ref_rec_part == DEFAULT)
#if PRESETS_TUNE
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected || picture_control_set_ptr->enc_mode <= ENC_M1)
#else
        if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected || picture_control_set_ptr->enc_mode == ENC_M0)
#endif
            context_ptr->prune_ref_frame_for_rec_partitions = 0;
        else
            context_ptr->prune_ref_frame_for_rec_partitions = 1;
    else
        context_ptr->prune_ref_frame_for_rec_partitions = sequence_control_set_ptr->static_config.prune_ref_rec_part;


    // Derive INTER/INTER WEDGE variance TH
    // Phoenix: Active only when inter/inter compound is on
    if (picture_control_set_ptr->enc_mode <= ENC_M7)
        context_ptr->inter_inter_wedge_variance_th = 0;
    else
        context_ptr->inter_inter_wedge_variance_th = 100;

    // Derive MD Exit TH
#if MULTI_PASS_PD // Shut md_exit_th
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_exit_th = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_exit_th = (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected) ? 10 : 18;
    else
#endif
#if M0_OPT
    if (MR_MODE ||( picture_control_set_ptr->enc_mode == ENC_M0 && picture_control_set_ptr->parent_pcs_ptr->sc_content_detected == 0))
#else
    if (MR_MODE)
#endif
        context_ptr->md_exit_th = 0;
    else
        context_ptr->md_exit_th = (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected) ? 10 : 18;

    // md_stage_1_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_1_cand_prune_th
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
#if ENHANCED_M0_SETTINGS // lossless change - PD_PASS_2 and PD_PASS_1 should be completely decoupled: previous merge conflict
        context_ptr->md_stage_1_cand_prune_th = 75;
#else
        context_ptr->md_stage_1_cand_prune_th = sequence_control_set_ptr->static_config.md_stage_1_cand_prune_th;
#endif
    else
#endif
#if M0_OPT
    if (MR_MODE || (picture_control_set_ptr->enc_mode == ENC_M0 && (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected == 0)) || sequence_control_set_ptr->input_resolution == INPUT_SIZE_576p_RANGE_OR_LOWER)
#else
    if (MR_MODE || picture_control_set_ptr->enc_mode == ENC_M0 || sequence_control_set_ptr->input_resolution == INPUT_SIZE_576p_RANGE_OR_LOWER)
#endif
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;
    else if (picture_control_set_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_1_cand_prune_th = sequence_control_set_ptr->static_config.md_stage_1_cand_prune_th;
    else
        context_ptr->md_stage_1_cand_prune_th = (uint64_t)~0;



    // md_stage_1_class_prune_th (for class removal)
    // Remove class if deviation to the best higher than TH_C
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
#if ENHANCED_M0_SETTINGS // lossless change - PD_PASS_2 and PD_PASS_1 should be completely decoupled: previous merge conflict
        context_ptr->md_stage_1_class_prune_th = 100;
#else
        context_ptr->md_stage_1_class_prune_th = sequence_control_set_ptr->static_config.md_stage_1_class_prune_th;
#endif
    else
#endif
#if M0_OPT
    if (MR_MODE || (picture_control_set_ptr->enc_mode == ENC_M0 && (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected == 0)) || sequence_control_set_ptr->input_resolution == INPUT_SIZE_576p_RANGE_OR_LOWER)
#else
    if (MR_MODE || picture_control_set_ptr->enc_mode == ENC_M0 || sequence_control_set_ptr->input_resolution == INPUT_SIZE_576p_RANGE_OR_LOWER)
#endif
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;
    else if (picture_control_set_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_1_class_prune_th = sequence_control_set_ptr->static_config.md_stage_1_class_prune_th;
    else
        context_ptr->md_stage_1_class_prune_th = (uint64_t)~0;

    // md_stage_2_cand_prune_th (for single candidate removal per class)
    // Remove candidate if deviation to the best is higher than md_stage_2_cand_prune_th
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_cand_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_stage_2_cand_prune_th = sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE ? 5 : 3;
    else
#endif
#if M0_OPT || PRESETS_TUNE
#if PRESETS_TUNE && !PRESETS_OPT
        if (MR_MODE || picture_control_set_ptr->parent_pcs_ptr->sc_content_detected || picture_control_set_ptr->enc_mode <= ENC_M1)
#else
        if (MR_MODE || picture_control_set_ptr->parent_pcs_ptr->sc_content_detected || picture_control_set_ptr->enc_mode <= ENC_M0)
#endif
#else
    if (MR_MODE || picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
#endif
        context_ptr->md_stage_2_cand_prune_th = (uint64_t)~0;
#if !M0_OPT
    else if (picture_control_set_ptr->enc_mode <= ENC_M0)
        context_ptr->md_stage_2_cand_prune_th = sequence_control_set_ptr->static_config.md_stage_2_cand_prune_th;
#endif
    else if (picture_control_set_ptr->enc_mode <= ENC_M2)
        context_ptr->md_stage_2_cand_prune_th = sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE ? 15 : 12;
    else if (picture_control_set_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_2_cand_prune_th = sequence_control_set_ptr->input_resolution <= INPUT_SIZE_1080i_RANGE ? 5 : 3;
    else
        context_ptr->md_stage_2_cand_prune_th = (uint64_t)~0;

    // md_stage_2_class_prune_th (for class removal)
    // Remove class if deviation to the best is higher than md_stage_2_class_prune_th
#if MULTI_PASS_PD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_stage_2_class_prune_th = (uint64_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
#if ENHANCED_M0_SETTINGS // lossless change - PD_PASS_2 and PD_PASS_1 should be completely decoupled: previous merge conflict
        context_ptr->md_stage_2_class_prune_th = 25;
#else
        context_ptr->md_stage_2_class_prune_th = sequence_control_set_ptr->static_config.md_stage_2_class_prune_th;
#endif
    else
#endif
#if M0_OPT
    if (MR_MODE)
#else
    if (MR_MODE || picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
#endif
        context_ptr->md_stage_2_class_prune_th = (uint64_t)~0;
    else if (picture_control_set_ptr->enc_mode <= ENC_M4)
        context_ptr->md_stage_2_class_prune_th = sequence_control_set_ptr->static_config.md_stage_2_class_prune_th;
    else // to be tested for m5-m8
        context_ptr->md_stage_2_class_prune_th = (uint64_t)~0;


    // Weighting (expressed as a percentage) applied to
    // square shape costs for determining if a and b
    // shapes should be skipped. Namely:
    // skip HA and HB if h_cost > (weighted sq_cost)
    // skip VA and VB if v_cost > (weighted sq_cost)

    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->sq_weight = (uint32_t)~0;
    else if (context_ptr->pd_pass == PD_PASS_1)
#if ENHANCED_M0_SETTINGS // lossless change - PD_PASS_2 and PD_PASS_1 should be completely decoupled: previous merge conflict
        context_ptr->sq_weight = 100;
#else
        context_ptr->sq_weight = sequence_control_set_ptr->static_config.sq_weight;
#endif
    else

    if (MR_MODE)
        context_ptr->sq_weight = (uint32_t)~0;
    else
        context_ptr->sq_weight = sequence_control_set_ptr->static_config.sq_weight;


#if MULTI_PASS_PD
    // Set pred ME full search area
    if (context_ptr->pd_pass == PD_PASS_0) {
        context_ptr->full_pel_ref_window_width_th = FULL_PEL_REF_WINDOW_WIDTH;
        context_ptr->full_pel_ref_window_height_th = FULL_PEL_REF_WINDOW_HEIGHT;
    }
    else if (context_ptr->pd_pass == PD_PASS_1) {
        context_ptr->full_pel_ref_window_width_th = FULL_PEL_REF_WINDOW_WIDTH;
        context_ptr->full_pel_ref_window_height_th = FULL_PEL_REF_WINDOW_HEIGHT;
    }
    else {
        context_ptr->full_pel_ref_window_width_th = (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected == 0 && picture_control_set_ptr->enc_mode == ENC_M0) ? FULL_PEL_REF_WINDOW_WIDTH_EXTENDED : FULL_PEL_REF_WINDOW_WIDTH;
        context_ptr->full_pel_ref_window_height_th = (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected == 0 && picture_control_set_ptr->enc_mode == ENC_M0) ? FULL_PEL_REF_WINDOW_HEIGHT_EXTENDED : FULL_PEL_REF_WINDOW_HEIGHT;
    }

    // set compound_types_to_try
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->compound_types_to_try = MD_COMP_AVG;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->compound_types_to_try = MD_COMP_AVG;
    else {
        if (picture_control_set_ptr->parent_pcs_ptr->compound_mode)
            context_ptr->compound_types_to_try = picture_control_set_ptr->parent_pcs_ptr->compound_mode == 1 ? MD_COMP_DIFF0 : MD_COMP_WEDGE;
        else
            context_ptr->compound_types_to_try = MD_COMP_AVG;
    }
    // Set coeff_based_nsq_cand_reduction
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
#if ENHANCED_M0_SETTINGS
    else
        context_ptr->coeff_based_nsq_cand_reduction = EB_TRUE;
#else
    else if (MR_MODE)
        context_ptr->coeff_based_nsq_cand_reduction = EB_FALSE;
    else
        context_ptr->coeff_based_nsq_cand_reduction = EB_TRUE;
#endif

    // Set rdoq_quantize_fp @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->rdoq_quantize_fp = EB_FALSE;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->rdoq_quantize_fp = EB_FALSE;
    else
#if PRESETS_OPT
        context_ptr->rdoq_quantize_fp = (picture_control_set_ptr->enc_mode <= ENC_M7) ? EB_TRUE : EB_FALSE;
#elif PRESETS_TUNE
        context_ptr->rdoq_quantize_fp = (picture_control_set_ptr->enc_mode <= ENC_M2) ? EB_TRUE : EB_FALSE;
#else
        context_ptr->rdoq_quantize_fp = (picture_control_set_ptr->enc_mode == ENC_M0) ? EB_TRUE : EB_FALSE;
#endif

    // Set pic_obmc_mode @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_pic_obmc_mode = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_pic_obmc_mode = 0;
    else
        context_ptr->md_pic_obmc_mode = picture_control_set_ptr->parent_pcs_ptr->pic_obmc_mode;

    // Set enable_inter_intra @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_enable_inter_intra = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_enable_inter_intra = 0;
    else
        context_ptr->md_enable_inter_intra = picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra;

    // Set md_atb_mode @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_atb_mode = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_atb_mode = 0;
    else
        context_ptr->md_atb_mode = picture_control_set_ptr->parent_pcs_ptr->atb_mode;

    // Set md_filter_intra_mode @ MD
    if (context_ptr->pd_pass == PD_PASS_0)
        context_ptr->md_filter_intra_mode = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->md_filter_intra_mode = 0;
    else
        context_ptr->md_filter_intra_mode = picture_control_set_ptr->pic_filter_intra_mode;


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
        context_ptr->dc_cand_only_flag = (picture_control_set_ptr->slice_type == I_SLICE) ? EB_FALSE : EB_TRUE;
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
#endif
    // signal for enabling shortcut to skip search depths
    if (context_ptr->pd_pass == PD_PASS_0)
       context_ptr->enable_auto_max_partition = 0;
    else if (context_ptr->pd_pass == PD_PASS_1)
        context_ptr->enable_auto_max_partition = 0;
    else if (picture_control_set_ptr->enc_mode <= ENC_M0 || picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT)
        context_ptr->enable_auto_max_partition = 0;
    else
        context_ptr->enable_auto_max_partition = sequence_control_set_ptr->static_config.enable_auto_max_partition;


    return return_error;
}
void copy_neighbour_arrays(
    PictureControlSet   *picture_control_set_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t              src_idx,
    uint32_t              dst_idx,
    uint32_t              blk_mds,
    uint32_t              sb_org_x,
    uint32_t              sb_org_y);

static void set_parent_to_be_considered(
    MdcLcuData *resultsPtr,
    uint32_t    blk_index,
    int32_t     sb_size,
    int8_t      depth_step) {
    uint32_t  parent_depth_idx_mds, block_1d_idx;
    const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
    if (blk_geom->sq_size < ((sb_size == BLOCK_128X128) ? 128 : 64)) {
        //Set parent to be considered
        parent_depth_idx_mds = (blk_geom->sqi_mds - (blk_geom->quadi - 3) * ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth]) - parent_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
        const BlockGeom * parent_blk_geom = get_blk_geom_mds(parent_depth_idx_mds);
        uint32_t parent_tot_d1_blocks =
            parent_blk_geom->sq_size == 128 ? 17 :
            parent_blk_geom->sq_size > 8 ? 25 :
            parent_blk_geom->sq_size == 8 ? 5 : 1;
        for (block_1d_idx = 0; block_1d_idx < parent_tot_d1_blocks; block_1d_idx++) {
            resultsPtr->leaf_data_array[parent_depth_idx_mds + block_1d_idx].consider_block = 1;
        }

        if (depth_step < -1)
            set_parent_to_be_considered(
                resultsPtr,
                parent_depth_idx_mds,
                sb_size,
                depth_step + 1);
    }
}

static void set_child_to_be_considered(
    MdcLcuData *resultsPtr,
    uint32_t    blk_index,
    int32_t     sb_size,
    int8_t      depth_step) {
    uint32_t child_block_idx_1, child_block_idx_2, child_block_idx_3, child_block_idx_4;
    uint32_t tot_d1_blocks, block_1d_idx;
    const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
    tot_d1_blocks =
        blk_geom->sq_size == 128 ? 17 :
        blk_geom->sq_size > 8 ? 25 :
        blk_geom->sq_size == 8 ? 5 : 1;
    if (blk_geom->sq_size > 4) {
        for (block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
            resultsPtr->leaf_data_array[blk_index + block_1d_idx].consider_block = 1;
            resultsPtr->leaf_data_array[blk_index + block_1d_idx].refined_split_flag = EB_TRUE;
        }
        //Set first child to be considered
        child_block_idx_1 = blk_index + d1_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth];
        const BlockGeom * child1_blk_geom = get_blk_geom_mds(child_block_idx_1);
        uint32_t child1_tot_d1_blocks =
            child1_blk_geom->sq_size == 128 ? 17 :
            child1_blk_geom->sq_size > 8 ? 25 :
            child1_blk_geom->sq_size == 8 ? 5 : 1;

        for (block_1d_idx = 0; block_1d_idx < child1_tot_d1_blocks; block_1d_idx++) {
            resultsPtr->leaf_data_array[child_block_idx_1 + block_1d_idx].consider_block = 1;
            resultsPtr->leaf_data_array[child_block_idx_1 + block_1d_idx].refined_split_flag = EB_FALSE;
        }
        if (depth_step > 1)
            set_child_to_be_considered(
                resultsPtr,
                child_block_idx_1,
                sb_size,
                depth_step - 1);
        //Set second child to be considered
        child_block_idx_2 = child_block_idx_1 + ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom * child2_blk_geom = get_blk_geom_mds(child_block_idx_2);
        uint32_t child2_tot_d1_blocks =
            child2_blk_geom->sq_size == 128 ? 17 :
            child2_blk_geom->sq_size > 8 ? 25 :
            child2_blk_geom->sq_size == 8 ? 5 : 1;
        for (block_1d_idx = 0; block_1d_idx < child2_tot_d1_blocks; block_1d_idx++) {
            resultsPtr->leaf_data_array[child_block_idx_2 + block_1d_idx].consider_block = 1;
            resultsPtr->leaf_data_array[child_block_idx_2 + block_1d_idx].refined_split_flag = EB_FALSE;
        }
        if (depth_step > 1)
            set_child_to_be_considered(
                resultsPtr,
                child_block_idx_2,
                sb_size,
                depth_step - 1);
        //Set third child to be considered
        child_block_idx_3 = child_block_idx_2 + ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom * child3_blk_geom = get_blk_geom_mds(child_block_idx_3);
        uint32_t child3_tot_d1_blocks =
            child3_blk_geom->sq_size == 128 ? 17 :
            child3_blk_geom->sq_size > 8 ? 25 :
            child3_blk_geom->sq_size == 8 ? 5 : 1;

        for (block_1d_idx = 0; block_1d_idx < child3_tot_d1_blocks; block_1d_idx++) {
            resultsPtr->leaf_data_array[child_block_idx_3 + block_1d_idx].consider_block = 1;
            resultsPtr->leaf_data_array[child_block_idx_3 + block_1d_idx].refined_split_flag = EB_FALSE;
        }
        if (depth_step > 1)
            set_child_to_be_considered(
                resultsPtr,
                child_block_idx_3,
                sb_size,
                depth_step - 1);
        //Set forth child to be considered
        child_block_idx_4 = child_block_idx_3 + ns_depth_offset[sb_size == BLOCK_128X128][blk_geom->depth + 1];
        const BlockGeom * child4_blk_geom = get_blk_geom_mds(child_block_idx_4);
        uint32_t child4_tot_d1_blocks =
            child4_blk_geom->sq_size == 128 ? 17 :
            child4_blk_geom->sq_size > 8 ? 25 :
            child4_blk_geom->sq_size == 8 ? 5 : 1;
        for (block_1d_idx = 0; block_1d_idx < child4_tot_d1_blocks; block_1d_idx++) {
            resultsPtr->leaf_data_array[child_block_idx_4 + block_1d_idx].consider_block = 1;
            resultsPtr->leaf_data_array[child_block_idx_4 + block_1d_idx].refined_split_flag = EB_FALSE;
        }
        if (depth_step > 1)
            set_child_to_be_considered(
                resultsPtr,
                child_block_idx_4,
                sb_size,
                depth_step - 1);
    }
}

static void build_cand_block_array(
    SequenceControlSet *sequence_control_set_ptr,
    PictureControlSet  *picture_control_set_ptr,
    uint32_t            sb_index) {

    MdcLcuData *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
    resultsPtr->leaf_count = 0;
    uint32_t  blk_index = 0;
    uint32_t d1_blocks_accumlated, tot_d1_blocks = 0, d1_block_idx;

    EbBool split_flag;
    while (blk_index < sequence_control_set_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        split_flag = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        //if the parent sq is inside inject this block
        uint8_t is_blk_allowed = picture_control_set_ptr->slice_type != I_SLICE ? 1 : (blk_geom->sq_size < 128) ? 1 : 0;
        //init consider block flag
        if (sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed) {
            tot_d1_blocks =
                blk_geom->sq_size == 128 ? 17 :
                blk_geom->sq_size > 8 ? 25 :
                blk_geom->sq_size == 8 ? 5 : 1;
            d1_blocks_accumlated = 0;
            for (d1_block_idx = 0; d1_block_idx < tot_d1_blocks; d1_block_idx++)
                d1_blocks_accumlated += resultsPtr->leaf_data_array[blk_index + d1_block_idx].consider_block ? 1 : 0;

            for (uint32_t idx = 0; idx < tot_d1_blocks; ++idx) {
                if (resultsPtr->leaf_data_array[blk_index].consider_block) {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = d1_blocks_accumlated;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = 0;//valid only for square 85 world. will be removed.
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = blk_index;
                    split_flag = resultsPtr->leaf_data_array[blk_index].refined_split_flag;
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag;
                }
                blk_index++;
            }
        }
        blk_index += split_flag ? d1_depth_offset[sequence_control_set_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] - tot_d1_blocks : ns_depth_offset[sequence_control_set_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] - tot_d1_blocks;
    }
}

uint64_t  pd_level_tab[2][9][2][3] =
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

void derive_start_end_depth(
    PictureControlSet  *picture_control_set_ptr,
    SuperBlock         *sb_ptr,
    uint32_t sb_size,
    int8_t *s_depth,
    int8_t *e_depth,
    const BlockGeom * blk_geom) {

    uint8_t encode_mode = picture_control_set_ptr->parent_pcs_ptr->enc_mode;

    int8_t start_depth = sb_size == BLOCK_128X128 ? 0 : 1;
    int8_t end_depth = 5;
    int8_t depth = blk_geom->depth + start_depth;

    int8_t depthp1 = MIN(depth + 1, end_depth);
    int8_t depthp2 = MIN(depth + 2, end_depth);
    int8_t depthp3 = MIN(depth + 3, end_depth);

    uint8_t depthm1 = MAX(depth - 1, start_depth);
    uint8_t depthm2 = MAX(depth - 2, start_depth);
    uint8_t depthm3 = MAX(depth - 3, start_depth);

    uint64_t max_distance = 0xFFFFFFFFFFFFFFFF;

    uint64_t mth01 = pd_level_tab[!picture_control_set_ptr->parent_pcs_ptr->sc_content_detected && picture_control_set_ptr->slice_type != I_SLICE][encode_mode][0][0];
    uint64_t mth02 = pd_level_tab[!picture_control_set_ptr->parent_pcs_ptr->sc_content_detected && picture_control_set_ptr->slice_type != I_SLICE][encode_mode][0][1];
    uint64_t mth03 = pd_level_tab[!picture_control_set_ptr->parent_pcs_ptr->sc_content_detected && picture_control_set_ptr->slice_type != I_SLICE][encode_mode][0][2];
    uint64_t pth01 = pd_level_tab[!picture_control_set_ptr->parent_pcs_ptr->sc_content_detected && picture_control_set_ptr->slice_type != I_SLICE][encode_mode][1][0];
    uint64_t pth02 = pd_level_tab[!picture_control_set_ptr->parent_pcs_ptr->sc_content_detected && picture_control_set_ptr->slice_type != I_SLICE][encode_mode][1][1];
    uint64_t pth03 = pd_level_tab[!picture_control_set_ptr->parent_pcs_ptr->sc_content_detected && picture_control_set_ptr->slice_type != I_SLICE][encode_mode][1][2];

    uint64_t dist_001 =
        sb_ptr->depth_cost[depth] == 0 ?
        max_distance :
        sb_ptr->depth_cost[depthp1] <= sb_ptr->depth_cost[depth] ?
        0 :
        (((int64_t)sb_ptr->depth_cost[depthp1] - (int64_t)sb_ptr->depth_cost[depth]) * 100) / sb_ptr->depth_cost[depth];

    uint64_t dist_100 =
        sb_ptr->depth_cost[depth] == 0 ?
        max_distance :
        sb_ptr->depth_cost[depthm1] <= sb_ptr->depth_cost[depth] ?
        0 :
        (((int64_t)sb_ptr->depth_cost[depthm1] - (int64_t)sb_ptr->depth_cost[depth]) * 100) / sb_ptr->depth_cost[depth];

    uint64_t dist_002 =
        sb_ptr->depth_cost[depth] == 0 ?
        max_distance :
        sb_ptr->depth_cost[depthp2] <= sb_ptr->depth_cost[depth] ?
        0 :
        (((int64_t)sb_ptr->depth_cost[depthp2] - (int64_t)sb_ptr->depth_cost[depth]) * 100) / sb_ptr->depth_cost[depth];

    uint64_t dist_200 =
        sb_ptr->depth_cost[depth] == 0 ?
        max_distance :
        sb_ptr->depth_cost[depthm2] <= sb_ptr->depth_cost[depth] ?
        0 :
        (((int64_t)sb_ptr->depth_cost[depthm2] - (int64_t)sb_ptr->depth_cost[depth]) * 100) / sb_ptr->depth_cost[depth];

    uint64_t dist_003 =
        sb_ptr->depth_cost[depth] == 0 ?
        max_distance :
        sb_ptr->depth_cost[depthp3] <= sb_ptr->depth_cost[depth] ?
        0 :
        (((int64_t)sb_ptr->depth_cost[depthp3] - (int64_t)sb_ptr->depth_cost[depth]) * 100) / sb_ptr->depth_cost[depth];

    uint64_t dist_300 =
        sb_ptr->depth_cost[depth] == 0 ?
        max_distance :
        sb_ptr->depth_cost[depthm3] <= sb_ptr->depth_cost[depth] ?
        0 :
        (((int64_t)sb_ptr->depth_cost[depthm3] - (int64_t)sb_ptr->depth_cost[depth]) * 100) / sb_ptr->depth_cost[depth];

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

static void perform_pred_depth_refinement(
    SequenceControlSet  *sequence_control_set_ptr,
    PictureControlSet   *picture_control_set_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t             sb_index) {

    MdcLcuData *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
    uint32_t  blk_index = 0;

    // Reset mdc_sb_array data to defaults; it will be updated based on the predicted blocks (stored in md_cu_arr_nsq)
    while (blk_index < sequence_control_set_ptr->max_block_cnt) {
        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        resultsPtr->leaf_data_array[blk_index].consider_block = 0;
        resultsPtr->leaf_data_array[blk_index].split_flag = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        resultsPtr->leaf_data_array[blk_index].refined_split_flag = blk_geom->sq_size > 4 ? EB_TRUE : EB_FALSE;
        blk_index++;
    }

    resultsPtr->leaf_count = 0;
    blk_index = 0;

    SuperBlock  *sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];

    uint32_t tot_d1_blocks, block_1d_idx;
    EbBool split_flag;

    while (blk_index < sequence_control_set_ptr->max_block_cnt) {

        const BlockGeom * blk_geom = get_blk_geom_mds(blk_index);
        tot_d1_blocks =
            blk_geom->sq_size == 128 ? 17 :
            blk_geom->sq_size > 8 ? 25 :
            blk_geom->sq_size == 8 ? 5 : 1;

        // if the parent square is inside inject this block
        uint8_t is_blk_allowed = picture_control_set_ptr->slice_type != I_SLICE ? 1 : (blk_geom->sq_size < 128) ? 1 : 0;

        // derive split_flag
        split_flag = context_ptr->md_cu_arr_nsq[blk_index].split_flag;

        if (sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[blk_index] && is_blk_allowed) {
            if (blk_geom->shape == PART_N) {
                if (context_ptr->md_cu_arr_nsq[blk_index].split_flag == EB_FALSE) {

                    int8_t s_depth = 0;
                    int8_t e_depth = 0;

                    if (context_ptr->pd_pass == PD_PASS_0) {
                        derive_start_end_depth(
                            picture_control_set_ptr,
                            sb_ptr,
                            sequence_control_set_ptr->seq_header.sb_size,
                            &s_depth,
                            &e_depth,
                            blk_geom);
                    }
                    else if (context_ptr->pd_pass == PD_PASS_1) {

                        EbBool zero_coeff_present_flag = EB_FALSE;

                        if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2)
                            zero_coeff_present_flag = context_ptr->md_cu_arr_nsq[blk_index].block_has_coeff == 0;

                        else if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3) {
                            switch (blk_geom->bsize) {

                            case BLOCK_128X128:
                                zero_coeff_present_flag = (context_ptr->md_local_cu_unit[blk_index].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index].block_has_coeff == 0) || // SQ
                                    (context_ptr->md_local_cu_unit[blk_index + 1].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 1].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 2].block_has_coeff == 0) || // H
                                    (context_ptr->md_local_cu_unit[blk_index + 3].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 3].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 4].block_has_coeff == 0) || // V
                                    (context_ptr->md_local_cu_unit[blk_index + 5].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 5].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 6].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 7].block_has_coeff == 0) ||
                                    (context_ptr->md_local_cu_unit[blk_index + 8].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 8].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 9].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 10].block_has_coeff == 0) ||
                                    (context_ptr->md_local_cu_unit[blk_index + 11].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 11].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 12].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 13].block_has_coeff == 0) ||
                                    (context_ptr->md_local_cu_unit[blk_index + 14].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 14].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 15].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 16].block_has_coeff == 0);
                                break;

                            case BLOCK_64X64:
                            case BLOCK_32X32:
                            case BLOCK_16X16:
                                zero_coeff_present_flag = (context_ptr->md_local_cu_unit[blk_index].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index].block_has_coeff == 0) || // SQ
                                    (context_ptr->md_local_cu_unit[blk_index + 1].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 1].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 2].block_has_coeff == 0) || // H
                                    (context_ptr->md_local_cu_unit[blk_index + 3].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 3].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 4].block_has_coeff == 0) || // V
                                    (context_ptr->md_local_cu_unit[blk_index + 5].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 5].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 6].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 7].block_has_coeff == 0) ||
                                    (context_ptr->md_local_cu_unit[blk_index + 8].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 8].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 9].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 10].block_has_coeff == 0) ||
                                    (context_ptr->md_local_cu_unit[blk_index + 11].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 11].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 12].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 13].block_has_coeff == 0) ||
                                    (context_ptr->md_local_cu_unit[blk_index + 14].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 14].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 15].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 16].block_has_coeff == 0) ||
                                    (context_ptr->md_local_cu_unit[blk_index + 17].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 17].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 18].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 19].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 20].block_has_coeff == 0) ||
                                    (context_ptr->md_local_cu_unit[blk_index + 21].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 21].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 22].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 23].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 24].block_has_coeff == 0);
                                break;

                            case BLOCK_8X8:
                                zero_coeff_present_flag = (context_ptr->md_local_cu_unit[blk_index].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index].block_has_coeff == 0) || // SQ
                                    (context_ptr->md_local_cu_unit[blk_index + 1].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 1].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 2].block_has_coeff == 0) || // H
                                    (context_ptr->md_local_cu_unit[blk_index + 3].avail_blk_flag && context_ptr->md_cu_arr_nsq[blk_index + 3].block_has_coeff == 0 && context_ptr->md_cu_arr_nsq[blk_index + 4].block_has_coeff == 0);  // V
                                break;

                            case BLOCK_4X4:
                                zero_coeff_present_flag = (context_ptr->md_cu_arr_nsq[blk_index].block_has_coeff == 0); // SQ
                                break;

                            default:
                                assert(0);
                                break;
                            }
                        }

                        if (zero_coeff_present_flag) {
                            s_depth = 0;
                            e_depth = 0;
                        }
                        else

                            if (context_ptr->md_cu_arr_nsq[blk_index].best_d1_blk == blk_index) {
                                s_depth = -1;
                                e_depth = 0;
                            }
                            else {
                                s_depth = 0;
                                e_depth = 1;
                            }
                    }

                    // Add current pred depth block(s)
                    for (block_1d_idx = 0; block_1d_idx < tot_d1_blocks; block_1d_idx++) {
                        resultsPtr->leaf_data_array[blk_index + block_1d_idx].consider_block = 1;
                        resultsPtr->leaf_data_array[blk_index + block_1d_idx].refined_split_flag = EB_FALSE;
                    }

                    // Add block indices of upper depth(s)
                    if (s_depth != 0)
                        set_parent_to_be_considered(
                            resultsPtr,
                            blk_index,
                            sequence_control_set_ptr->seq_header.sb_size,
                            s_depth);

                    // Add block indices of lower depth(s)
                    if (e_depth != 0)
                        set_child_to_be_considered(
                            resultsPtr,
                            blk_index,
                            sequence_control_set_ptr->seq_header.sb_size,
                            e_depth);
                }
            }
        }
        blk_index += split_flag ? d1_depth_offset[sequence_control_set_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth] : ns_depth_offset[sequence_control_set_ptr->seq_header.sb_size == BLOCK_128X128][blk_geom->depth];
    }
}

/******************************************************
 * EncDec Kernel
 ******************************************************/
void* enc_dec_kernel(void *input_ptr)
{
    // Context & SCS & PCS
    EbThreadContext                       *thread_context_ptr = (EbThreadContext*)input_ptr;
    EncDecContext                         *context_ptr = (EncDecContext*)thread_context_ptr->priv;
    PictureControlSet                     *picture_control_set_ptr;
    SequenceControlSet                    *sequence_control_set_ptr;

    // Input
    EbObjectWrapper                       *encDecTasksWrapperPtr;
    EncDecTasks                           *encDecTasksPtr;

    // Output
    EbObjectWrapper                       *encDecResultsWrapperPtr;
    EncDecResults                         *encDecResultsPtr;
    // SB Loop variables
    SuperBlock                               *sb_ptr;
    uint16_t                                 sb_index;
    uint8_t                                  sb_sz;
    uint8_t                                  lcuSizeLog2;
    uint32_t                                 x_lcu_index;
    uint32_t                                 y_lcu_index;
    uint32_t                                 sb_origin_x;
    uint32_t                                 sb_origin_y;
    EbBool                                   lastLcuFlag;
    EbBool                                   endOfRowFlag;
    uint32_t                                 lcuRowIndexStart;
    uint32_t                                 lcuRowIndexCount;
    uint32_t                                 picture_width_in_sb;
    MdcLcuData                              *mdcPtr;

    // Variables
    EbBool                                   is16bit;

    // Segments
    //EbBool                                 initialProcessCall;
    uint16_t                                 segment_index;
    uint32_t                                 xLcuStartIndex;
    uint32_t                                 yLcuStartIndex;
    uint32_t                                 lcuStartIndex;
    uint32_t                                 lcuSegmentCount;
    uint32_t                                 lcuSegmentIndex;
    uint32_t                                 segmentRowIndex;
    uint32_t                                 segmentBandIndex;
    uint32_t                                 segmentBandSize;
    EncDecSegments                          *segments_ptr;

    segment_index = 0;

    for (;;) {
        // Get Mode Decision Results
        eb_get_full_object(
            context_ptr->mode_decision_input_fifo_ptr,
            &encDecTasksWrapperPtr);

        encDecTasksPtr = (EncDecTasks*)encDecTasksWrapperPtr->object_ptr;
        picture_control_set_ptr = (PictureControlSet*)encDecTasksPtr->picture_control_set_wrapper_ptr->object_ptr;
        sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
        segments_ptr = picture_control_set_ptr->enc_dec_segment_ctrl;
        lastLcuFlag = EB_FALSE;
        is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
        (void)is16bit;
        (void)endOfRowFlag;
#if !MULTI_PASS_PD
        // EncDec Kernel Signal(s) derivation

        signal_derivation_enc_dec_kernel_oq(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            context_ptr->md_context);
#endif
        // SB Constants
        sb_sz = (uint8_t)sequence_control_set_ptr->sb_size_pix;
        lcuSizeLog2 = (uint8_t)Log2f(sb_sz);
        context_ptr->sb_sz = sb_sz;
        picture_width_in_sb = (sequence_control_set_ptr->seq_header.max_frame_width + sb_sz - 1) >> lcuSizeLog2;
        endOfRowFlag = EB_FALSE;
        lcuRowIndexStart = lcuRowIndexCount = 0;
        context_ptr->tot_intra_coded_area = 0;

        // Segment-loop
        while (assign_enc_dec_segments(segments_ptr, &segment_index, encDecTasksPtr, context_ptr->enc_dec_feedback_fifo_ptr) == EB_TRUE)
        {
            xLcuStartIndex = segments_ptr->x_start_array[segment_index];
            yLcuStartIndex = segments_ptr->y_start_array[segment_index];
            lcuStartIndex = yLcuStartIndex * picture_width_in_sb + xLcuStartIndex;
            lcuSegmentCount = segments_ptr->valid_lcu_count_array[segment_index];

            segmentRowIndex = segment_index / segments_ptr->segment_band_count;
            segmentBandIndex = segment_index - segmentRowIndex * segments_ptr->segment_band_count;
            segmentBandSize = (segments_ptr->lcu_band_count * (segmentBandIndex + 1) + segments_ptr->segment_band_count - 1) / segments_ptr->segment_band_count;

            // Reset Coding Loop State
            reset_mode_decision(
                sequence_control_set_ptr,
                context_ptr->md_context,
                picture_control_set_ptr,
                segment_index);

            // Reset EncDec Coding State
            reset_enc_dec(    // HT done
                context_ptr,
                picture_control_set_ptr,
                sequence_control_set_ptr,
                segment_index);

            if (picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL)
                ((EbReferenceObject  *)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->average_intensity = picture_control_set_ptr->parent_pcs_ptr->average_intensity[0];
            for (y_lcu_index = yLcuStartIndex, lcuSegmentIndex = lcuStartIndex; lcuSegmentIndex < lcuStartIndex + lcuSegmentCount; ++y_lcu_index) {
                for (x_lcu_index = xLcuStartIndex; x_lcu_index < picture_width_in_sb && (x_lcu_index + y_lcu_index < segmentBandSize) && lcuSegmentIndex < lcuStartIndex + lcuSegmentCount; ++x_lcu_index, ++lcuSegmentIndex) {
                    sb_index = (uint16_t)(y_lcu_index * picture_width_in_sb + x_lcu_index);
                    sb_ptr = picture_control_set_ptr->sb_ptr_array[sb_index];
                    sb_origin_x = x_lcu_index << lcuSizeLog2;
                    sb_origin_y = y_lcu_index << lcuSizeLog2;
                    lastLcuFlag = (sb_index == sequence_control_set_ptr->sb_tot_cnt - 1) ? EB_TRUE : EB_FALSE;
                    endOfRowFlag = (x_lcu_index == picture_width_in_sb - 1) ? EB_TRUE : EB_FALSE;
                    lcuRowIndexStart = (x_lcu_index == picture_width_in_sb - 1 && lcuRowIndexCount == 0) ? y_lcu_index : lcuRowIndexStart;
                    lcuRowIndexCount = (x_lcu_index == picture_width_in_sb - 1) ? lcuRowIndexCount + 1 : lcuRowIndexCount;
                    mdcPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
                    context_ptr->sb_index = sb_index;

                    if (picture_control_set_ptr->update_cdf) {
                        picture_control_set_ptr->rate_est_array[sb_index] = *picture_control_set_ptr->md_rate_estimation_array;
#if CABAC_SERIAL
                        if (sb_index == 0)
                            picture_control_set_ptr->ec_ctx_array[sb_index] = *picture_control_set_ptr->coeff_est_entropy_coder_ptr->fc;
                        else
                            picture_control_set_ptr->ec_ctx_array[sb_index] = picture_control_set_ptr->ec_ctx_array[sb_index - 1];
#else
                        // Use the latest available CDF for the current SB
                        // Use the weighted average of left (3x) and top (1x) if available.
                        int8_t up_available     = ((int32_t)(sb_origin_y >> MI_SIZE_LOG2) > sb_ptr->tile_info.mi_row_start);
                        int8_t left_available   = ((int32_t)(sb_origin_x >> MI_SIZE_LOG2) > sb_ptr->tile_info.mi_col_start);
                        if (!left_available && !up_available)
                            picture_control_set_ptr->ec_ctx_array[sb_index] = *picture_control_set_ptr->coeff_est_entropy_coder_ptr->fc;
                        else if (!left_available)
                            picture_control_set_ptr->ec_ctx_array[sb_index] = picture_control_set_ptr->ec_ctx_array[sb_index - picture_width_in_sb];
                        else if (!up_available)
                            picture_control_set_ptr->ec_ctx_array[sb_index] = picture_control_set_ptr->ec_ctx_array[sb_index - 1];
                        else {
                            picture_control_set_ptr->ec_ctx_array[sb_index] = picture_control_set_ptr->ec_ctx_array[sb_index - 1];
                            avg_cdf_symbols(
                                &picture_control_set_ptr->ec_ctx_array[sb_index],
                                &picture_control_set_ptr->ec_ctx_array[sb_index - picture_width_in_sb],
                                AVG_CDF_WEIGHT_LEFT,
                                AVG_CDF_WEIGHT_TOP);
                        }

#endif
                        // Initial Rate Estimation of the syntax elements
                        av1_estimate_syntax_rate(
                            &picture_control_set_ptr->rate_est_array[sb_index],
                            picture_control_set_ptr->slice_type == I_SLICE,
                            &picture_control_set_ptr->ec_ctx_array[sb_index]);
                        // Initial Rate Estimation of the Motion vectors
                        av1_estimate_mv_rate(
                            picture_control_set_ptr,
                            &picture_control_set_ptr->rate_est_array[sb_index],
                            &picture_control_set_ptr->ec_ctx_array[sb_index]);

                        av1_estimate_coefficients_rate(
                            &picture_control_set_ptr->rate_est_array[sb_index],
                            &picture_control_set_ptr->ec_ctx_array[sb_index]);

                        //let the candidate point to the new rate table.
                        uint32_t  candidateIndex;
                        for (candidateIndex = 0; candidateIndex < MODE_DECISION_CANDIDATE_MAX_COUNT; ++candidateIndex)
                            context_ptr->md_context->fast_candidate_ptr_array[candidateIndex]->md_rate_estimation_ptr = &picture_control_set_ptr->rate_est_array[sb_index];
                        context_ptr->md_context->md_rate_estimation_ptr = &picture_control_set_ptr->rate_est_array[sb_index];
                    }
                    // Configure the LCU
                    mode_decision_configure_lcu(
                        context_ptr->md_context,
                        picture_control_set_ptr,
                        (uint8_t)sb_ptr->qp);

#if MULTI_PASS_PD
                    // Multi-Pass PD Path
                    // For each SB, all blocks are tested in PD0 (4421 blocks if 128x128 SB, and 1101 blocks if 64x64 SB).
                    // Then the PD0 predicted Partitioning Structure is refined by considering up to three refinements depths away from the predicted depth, both in the direction of smaller block sizes and in the direction of larger block sizes (up to Pred - 3 / Pred + 3 refinement). The selection of the refinement depth is performed using the cost
                    // deviation between the current depth cost and candidate depth cost. The generated blocks are used as input candidates to PD1.
                    // The PD1 predicted Partitioning Structure is also refined (up to Pred - 1 / Pred + 1 refinement) using the square (SQ) vs. non-square (NSQ) decision(s)
                    // inside the predicted depth and using coefficient information. The final set of blocks is evaluated in PD2 to output the final Partitioning Structure

                    if ((picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_0 ||
                         picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_1 ||
                         picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2 ||
                         picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3) && sequence_control_set_ptr->sb_geom[sb_index].is_complete_sb) {

                        // Save a clean copy of the neighbor arrays
                        copy_neighbour_arrays(
                            picture_control_set_ptr,
                            context_ptr->md_context,
                            MD_NEIGHBOR_ARRAY_INDEX,
                            MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                            0,
                            sb_origin_x,
                            sb_origin_y);

                        // [PD_PASS_0] Signal(s) derivation
                        context_ptr->md_context->pd_pass = PD_PASS_0;
                        signal_derivation_enc_dec_kernel_oq(
                            sequence_control_set_ptr,
                            picture_control_set_ptr,
                            context_ptr->md_context);

                        // [PD_PASS_0] Mode Decision - Reduce the total number of partitions to be tested in later stages.
                        // Input : mdc_cu_ptr built @ mdc process (up to 4421)
                        // Output: md_cu_arr_nsq reduced set of block(s)

                        // PD0 MD Tool(s) : Best ME candidate only as INTER candidate(s), DC only as INTRA candidate(s), Chroma blind, Spatial SSE,
                        // no MVP table generation, no fast rate @ full cost derivation, Md-Stage 0 and Md-Stage 2 using count=1 (i.e. only best md-stage-0 candidate)
                        mode_decision_sb(
                            sequence_control_set_ptr,
                            picture_control_set_ptr,
                            mdcPtr,
                            sb_ptr,
                            sb_origin_x,
                            sb_origin_y,
                            sb_index,
                            context_ptr->md_context);

                        // Perform Pred_0 depth refinement - Add blocks to be considered in the next stage(s) of PD based on depth cost.
                        perform_pred_depth_refinement(
                            sequence_control_set_ptr,
                            picture_control_set_ptr,
                            context_ptr->md_context,
                            sb_index);

                        // Re-build mdc_cu_ptr for the 2nd PD Pass [PD_PASS_1]
                        build_cand_block_array(
                            sequence_control_set_ptr,
                            picture_control_set_ptr,
                            sb_index);

                        // Reset neighnor information to current SB @ position (0,0)
                        copy_neighbour_arrays(
                            picture_control_set_ptr,
                            context_ptr->md_context,
                            MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX,
                            MD_NEIGHBOR_ARRAY_INDEX,
                            0,
                            sb_origin_x,
                            sb_origin_y);

                        if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_1 ||
                            picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2 ||
                            picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3 ){

                            // [PD_PASS_1] Signal(s) derivation
                            context_ptr->md_context->pd_pass = PD_PASS_1;
                            signal_derivation_enc_dec_kernel_oq(
                                sequence_control_set_ptr,
                                picture_control_set_ptr,
                                context_ptr->md_context);

                            // [PD_PASS_1] Mode Decision - Further reduce the number of
                            // partitions to be considered in later PD stages. This pass uses more accurate
                            // info than PD0 to give a better PD estimate.
                            // Input : mdc_cu_ptr built @ PD0 refinement
                            // Output: md_cu_arr_nsq reduced set of block(s)

                            // PD1 MD Tool(s) : ME and Predictive ME only as INTER candidate(s) but MRP blind (only reference index 0 for motion compensation),
                            // DC only as INTRA candidate(s)
                            mode_decision_sb(
                                sequence_control_set_ptr,
                                picture_control_set_ptr,
                                mdcPtr,
                                sb_ptr,
                                sb_origin_x,
                                sb_origin_y,
                                sb_index,
                                context_ptr->md_context);

                            // Perform Pred_1 depth refinement - Add blocks to be considered in the next stage(s) of PD based on depth cost.
                            perform_pred_depth_refinement(
                                sequence_control_set_ptr,
                                picture_control_set_ptr,
                                context_ptr->md_context,
                                sb_index);

                            // Re-build mdc_cu_ptr for the 3rd PD Pass [PD_PASS_2]
                            build_cand_block_array(
                                sequence_control_set_ptr,
                                picture_control_set_ptr,
                                sb_index);

                            // Reset neighnor information to current SB @ position (0,0)
                            copy_neighbour_arrays(
                                picture_control_set_ptr,
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
                    signal_derivation_enc_dec_kernel_oq(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        context_ptr->md_context);

                    // [PD_PASS_2] Mode Decision - Obtain the final partitioning decision using more accurate info
                    // than previous stages.  Reduce the total number of partitions to 1.
                    // Input : mdc_cu_ptr built @ PD1 refinement
                    // Output: md_cu_arr_nsq reduced set of block(s)

                    // PD2 MD Tool(s): default MD Tool(s)
#endif
                    mode_decision_sb(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        mdcPtr,
                        sb_ptr,
                        sb_origin_x,
                        sb_origin_y,
                        sb_index,
                        context_ptr->md_context);

                    // Configure the LCU
                    enc_dec_configure_lcu(
                        context_ptr,
                        sb_ptr,
                        picture_control_set_ptr,
                        (uint8_t)sb_ptr->qp);

#if NO_ENCDEC
                    no_enc_dec_pass(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        sb_ptr,
                        sb_index,
                        sb_origin_x,
                        sb_origin_y,
                        sb_ptr->qp,
                        context_ptr);
#else
                    // Encode Pass
                    av1_encode_pass(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        sb_ptr,
                        sb_index,
                        sb_origin_x,
                        sb_origin_y,
                        context_ptr);
#endif

                    if (picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL)
                        ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->intra_coded_area_sb[sb_index] = (uint8_t)((100 * context_ptr->intra_coded_area_sb[sb_index]) / (64 * 64));
                }
                xLcuStartIndex = (xLcuStartIndex > 0) ? xLcuStartIndex - 1 : 0;
            }
        }

        eb_block_on_mutex(picture_control_set_ptr->intra_mutex);
        picture_control_set_ptr->intra_coded_area += (uint32_t)context_ptr->tot_intra_coded_area;
        eb_release_mutex(picture_control_set_ptr->intra_mutex);

        if (lastLcuFlag) {
            // Copy film grain data from parent picture set to the reference object for further reference
            if (sequence_control_set_ptr->seq_header.film_grain_params_present)
            {
                if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE && picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr) {
                    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->film_grain_params
                        = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.film_grain_params;
                }
            }
            if (picture_control_set_ptr->parent_pcs_ptr->frame_end_cdf_update_mode && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE && picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr)
                for (int frame = LAST_FRAME; frame <= ALTREF_FRAME; ++frame)
                    ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->global_motion[frame]
                    = picture_control_set_ptr->parent_pcs_ptr->global_motion[frame];
            EB_MEMCPY(picture_control_set_ptr->parent_pcs_ptr->av1x->sgrproj_restore_cost, context_ptr->md_rate_estimation_ptr->sgrproj_restore_fac_bits, 2 * sizeof(int32_t));
            EB_MEMCPY(picture_control_set_ptr->parent_pcs_ptr->av1x->switchable_restore_cost, context_ptr->md_rate_estimation_ptr->switchable_restore_fac_bits, 3 * sizeof(int32_t));
            EB_MEMCPY(picture_control_set_ptr->parent_pcs_ptr->av1x->wiener_restore_cost, context_ptr->md_rate_estimation_ptr->wiener_restore_fac_bits, 2 * sizeof(int32_t));
            picture_control_set_ptr->parent_pcs_ptr->av1x->rdmult = context_ptr->full_lambda;
        }

        if (lastLcuFlag)
        {
            // Get Empty EncDec Results
            eb_get_empty_object(
                context_ptr->enc_dec_output_fifo_ptr,
                &encDecResultsWrapperPtr);
            encDecResultsPtr = (EncDecResults*)encDecResultsWrapperPtr->object_ptr;
            encDecResultsPtr->picture_control_set_wrapper_ptr = encDecTasksPtr->picture_control_set_wrapper_ptr;
            //CHKN these are not needed for DLF
            encDecResultsPtr->completed_lcu_row_index_start = 0;
            encDecResultsPtr->completed_lcu_row_count = ((sequence_control_set_ptr->seq_header.max_frame_height + sequence_control_set_ptr->sb_size_pix - 1) >> lcuSizeLog2);
            // Post EncDec Results
            eb_post_full_object(encDecResultsWrapperPtr);
        }
        // Release Mode Decision Results
        eb_release_object(encDecTasksWrapperPtr);
    }
    return EB_NULL;
}

void eb_av1_add_film_grain(EbPictureBufferDesc *src,
    EbPictureBufferDesc *dst,
    aom_film_grain_t *film_grain_ptr) {
    uint8_t *luma, *cb, *cr;
    int32_t height, width, luma_stride, chroma_stride;
    int32_t use_high_bit_depth = 0;
    int32_t chroma_subsamp_x = 0;
    int32_t chroma_subsamp_y = 0;

    aom_film_grain_t params = *film_grain_ptr;

    switch (src->bit_depth) {
    case EB_8BIT:
        params.bit_depth = 8;
        use_high_bit_depth = 0;
        chroma_subsamp_x = 1;
        chroma_subsamp_y = 1;
        break;
    case EB_10BIT:
        params.bit_depth = 10;
        use_high_bit_depth = 1;
        chroma_subsamp_x = 1;
        chroma_subsamp_y = 1;
        break;
    default:  //todo: Throw an error if unknown format?
        params.bit_depth = 10;
        use_high_bit_depth = 1;
        chroma_subsamp_x = 1;
        chroma_subsamp_y = 1;
    }

    dst->max_width = src->max_width;
    dst->max_height = src->max_height;

    fgn_copy_rect(src->buffer_y + ((src->origin_y * src->stride_y + src->origin_x) << use_high_bit_depth), src->stride_y,
        dst->buffer_y + ((dst->origin_y * dst->stride_y + dst->origin_x) << use_high_bit_depth), dst->stride_y,
        dst->width, dst->height, use_high_bit_depth);

    fgn_copy_rect(src->buffer_cb + ((src->stride_cb * (src->origin_y >> chroma_subsamp_y)
        + (src->origin_x >> chroma_subsamp_x)) << use_high_bit_depth), src->stride_cb,
        dst->buffer_cb + ((dst->stride_cb * (dst->origin_y >> chroma_subsamp_y)
            + (dst->origin_x >> chroma_subsamp_x)) << use_high_bit_depth), dst->stride_cb,
        dst->width >> chroma_subsamp_x, dst->height >> chroma_subsamp_y,
        use_high_bit_depth);

    fgn_copy_rect(src->buffer_cr + ((src->stride_cr * (src->origin_y >> chroma_subsamp_y)
        + (src->origin_x >> chroma_subsamp_x)) << use_high_bit_depth), src->stride_cr,
        dst->buffer_cr + ((dst->stride_cr * (dst->origin_y >> chroma_subsamp_y)
            + (dst->origin_x >> chroma_subsamp_x)) << use_high_bit_depth), dst->stride_cr,
        dst->width >> chroma_subsamp_x, dst->height >> chroma_subsamp_y,
        use_high_bit_depth);

    luma = dst->buffer_y + ((dst->origin_y * dst->stride_y + dst->origin_x) << use_high_bit_depth);
    cb = dst->buffer_cb + ((dst->stride_cb * (dst->origin_y >> chroma_subsamp_y)
        + (dst->origin_x >> chroma_subsamp_x)) << use_high_bit_depth);
    cr = dst->buffer_cr + ((dst->stride_cr * (dst->origin_y >> chroma_subsamp_y)
        + (dst->origin_x >> chroma_subsamp_x)) << use_high_bit_depth);

    luma_stride = dst->stride_y;
    chroma_stride = dst->stride_cb;

    width = dst->width;
    height = dst->height;

    eb_av1_add_film_grain_run(&params, luma, cb, cr, height, width, luma_stride,
        chroma_stride, use_high_bit_depth, chroma_subsamp_y,
        chroma_subsamp_x);
    return;
}
// clang-format on
