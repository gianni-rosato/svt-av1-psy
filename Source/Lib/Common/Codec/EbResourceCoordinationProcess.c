/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureBufferDesc.h"
#include "EbResourceCoordinationProcess.h"
#include "EbResourceCoordinationResults.h"
#include "EbTransforms.h"
#include "EbTime.h"

void resource_coordination_context_dctor(EbPtr p)
{
    ResourceCoordinationContext *obj = (ResourceCoordinationContext*)p;

    EB_FREE_ARRAY(obj->sequenceControlSetActiveArray);
    EB_FREE_ARRAY(obj->picture_number_array);
}

/************************************************
 * Resource Coordination Context Constructor
 ************************************************/
EbErrorType resource_coordination_context_ctor(
    ResourceCoordinationContext *context_ptr,
    EbFifo                        *inputBufferFifoPtr,
    EbFifo                        *resource_coordination_results_output_fifo_ptr,
    EbFifo                        **picture_control_set_fifo_ptr_array,
    EbSequenceControlSetInstance  **sequence_control_set_instance_array,
    EbFifo                         *sequence_control_set_empty_fifo_ptr,
    EbCallback                    **app_callback_ptr_array,
    uint32_t                       compute_segments_total_count_array,
    uint32_t                        encode_instances_total_count){

    context_ptr->dctor = resource_coordination_context_dctor;

    context_ptr->input_buffer_fifo_ptr = inputBufferFifoPtr;
    context_ptr->resource_coordination_results_output_fifo_ptr = resource_coordination_results_output_fifo_ptr;
    context_ptr->picture_control_set_fifo_ptr_array = picture_control_set_fifo_ptr_array;
    context_ptr->sequence_control_set_instance_array = sequence_control_set_instance_array;
    context_ptr->sequence_control_set_empty_fifo_ptr = sequence_control_set_empty_fifo_ptr;
    context_ptr->app_callback_ptr_array = app_callback_ptr_array;
    context_ptr->compute_segments_total_count_array = compute_segments_total_count_array;
    context_ptr->encode_instances_total_count = encode_instances_total_count;

    // Allocate SequenceControlSetActiveArray
    EB_CALLOC_ARRAY(context_ptr->sequenceControlSetActiveArray, context_ptr->encode_instances_total_count);

    EB_CALLOC_ARRAY(context_ptr->picture_number_array, context_ptr->encode_instances_total_count);

    context_ptr->average_enc_mod = 0;
    context_ptr->prev_enc_mod = 0;
    context_ptr->prev_enc_mode_delta = 0;
    context_ptr->cur_speed = 0; // speed x 1000
    context_ptr->previous_mode_change_buffer = 0;
    context_ptr->first_in_pic_arrived_time_seconds = 0;
    context_ptr->first_in_pic_arrived_timeu_seconds = 0;
    context_ptr->previous_frame_in_check1 = 0;
    context_ptr->previous_frame_in_check2 = 0;
    context_ptr->previous_frame_in_check3 = 0;
    context_ptr->previous_mode_change_frame_in = 0;
    context_ptr->prevs_time_seconds = 0;
    context_ptr->prevs_timeu_seconds = 0;
    context_ptr->prev_frame_out = 0;
    context_ptr->start_flag = EB_FALSE;

    context_ptr->previous_buffer_check1 = 0;
    context_ptr->prev_change_cond = 0;

    return EB_ErrorNone;
}

/******************************************************
* Derive Pre-Analysis settings for OQ
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
EbErrorType signal_derivation_pre_analysis_oq(
    SequenceControlSet       *sequence_control_set_ptr,
    PictureParentControlSet  *picture_control_set_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    uint8_t input_resolution = sequence_control_set_ptr->input_resolution;

    // HME Flags updated @ signal_derivation_multi_processes_oq
    uint8_t  hme_me_level = picture_control_set_ptr->enc_mode;

    // Derive HME Flag
    if (sequence_control_set_ptr->static_config.use_default_me_hme) {
        picture_control_set_ptr->enable_hme_flag = EB_TRUE;
        picture_control_set_ptr->enable_hme_level0_flag = enable_hme_level0_flag[0][input_resolution][hme_me_level] || enable_hme_level0_flag[1][input_resolution][hme_me_level];
        picture_control_set_ptr->enable_hme_level1_flag = enable_hme_level1_flag[0][input_resolution][hme_me_level] || enable_hme_level1_flag[1][input_resolution][hme_me_level];
        picture_control_set_ptr->enable_hme_level2_flag = enable_hme_level2_flag[0][input_resolution][hme_me_level] || enable_hme_level2_flag[1][input_resolution][hme_me_level];
    }
    else {
        picture_control_set_ptr->enable_hme_flag = sequence_control_set_ptr->static_config.enable_hme_flag;
        picture_control_set_ptr->enable_hme_level0_flag = sequence_control_set_ptr->static_config.enable_hme_level0_flag;
        picture_control_set_ptr->enable_hme_level1_flag = sequence_control_set_ptr->static_config.enable_hme_level1_flag;
        picture_control_set_ptr->enable_hme_level2_flag = sequence_control_set_ptr->static_config.enable_hme_level2_flag;
    }

    picture_control_set_ptr->tf_enable_hme_flag = EB_TRUE;
    picture_control_set_ptr->tf_enable_hme_level0_flag = tf_enable_hme_level0_flag[0][input_resolution][hme_me_level] || tf_enable_hme_level0_flag[1][input_resolution][hme_me_level];
    picture_control_set_ptr->tf_enable_hme_level1_flag = tf_enable_hme_level1_flag[0][input_resolution][hme_me_level] || tf_enable_hme_level1_flag[1][input_resolution][hme_me_level];
    picture_control_set_ptr->tf_enable_hme_level2_flag = tf_enable_hme_level2_flag[0][input_resolution][hme_me_level] || tf_enable_hme_level2_flag[1][input_resolution][hme_me_level];

    if (picture_control_set_ptr->enc_mode >= ENC_M8)
        sequence_control_set_ptr->seq_header.enable_restoration = 0;

    return return_error;
}

//******************************************************************************//
// Modify the Enc mode based on the buffer Status
// Inputs: TargetSpeed, Status of the SCbuffer
// Output: EncMod
//******************************************************************************//
void SpeedBufferControl(
    ResourceCoordinationContext   *context_ptr,
    PictureParentControlSet       *picture_control_set_ptr,
    SequenceControlSet            *sequence_control_set_ptr)
{
    uint64_t cursTimeSeconds = 0;
    uint64_t cursTimeuSeconds = 0;
    double overallDuration = 0.0;
    double instDuration = 0.0;
    int8_t  encoderModeDelta = 0;
    int64_t inputFramesCount = 0;
    int8_t changeCond = 0;
    int64_t targetFps = (sequence_control_set_ptr->static_config.injector_frame_rate >> 16);

    int64_t bufferTrshold1 = SC_FRAMES_INTERVAL_T1;
    int64_t bufferTrshold2 = SC_FRAMES_INTERVAL_T2;
    int64_t bufferTrshold3 = MIN(targetFps * 3, SC_FRAMES_INTERVAL_T3);
    eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->sc_buffer_mutex);

    if (sequence_control_set_ptr->encode_context_ptr->sc_frame_in == 0)
        EbStartTime(&context_ptr->first_in_pic_arrived_time_seconds, &context_ptr->first_in_pic_arrived_timeu_seconds);
    else if (sequence_control_set_ptr->encode_context_ptr->sc_frame_in == SC_FRAMES_TO_IGNORE)
        context_ptr->start_flag = EB_TRUE;
    // Compute duration since the start of the encode and since the previous checkpoint
    EbFinishTime(&cursTimeSeconds, &cursTimeuSeconds);

    EbComputeOverallElapsedTimeMs(
        context_ptr->first_in_pic_arrived_time_seconds,
        context_ptr->first_in_pic_arrived_timeu_seconds,
        cursTimeSeconds,
        cursTimeuSeconds,
        &overallDuration);

    EbComputeOverallElapsedTimeMs(
        context_ptr->prevs_time_seconds,
        context_ptr->prevs_timeu_seconds,
        cursTimeSeconds,
        cursTimeuSeconds,
        &instDuration);

    inputFramesCount = (int64_t)overallDuration *(sequence_control_set_ptr->static_config.injector_frame_rate >> 16) / 1000;
    sequence_control_set_ptr->encode_context_ptr->sc_buffer = inputFramesCount - sequence_control_set_ptr->encode_context_ptr->sc_frame_in;

    encoderModeDelta = 0;

    // Check every bufferTsshold1 for the changes (previous_frame_in_check1 variable)
    if ((sequence_control_set_ptr->encode_context_ptr->sc_frame_in > context_ptr->previous_frame_in_check1 + bufferTrshold1 && sequence_control_set_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)) {
        // Go to a slower mode based on the fullness and changes of the buffer
        if (sequence_control_set_ptr->encode_context_ptr->sc_buffer < targetFps && (context_ptr->prev_enc_mode_delta > -1 || (context_ptr->prev_enc_mode_delta < 0 && sequence_control_set_ptr->encode_context_ptr->sc_frame_in > context_ptr->previous_mode_change_frame_in + targetFps * 2))) {
            if (context_ptr->previous_buffer_check1 > sequence_control_set_ptr->encode_context_ptr->sc_buffer + bufferTrshold1) {
                encoderModeDelta += -1;
                changeCond = 2;
            }
            else if (context_ptr->previous_mode_change_buffer > bufferTrshold1 + sequence_control_set_ptr->encode_context_ptr->sc_buffer && sequence_control_set_ptr->encode_context_ptr->sc_buffer < bufferTrshold1) {
                encoderModeDelta += -1;
                changeCond = 4;
            }
        }

        // Go to a faster mode based on the fullness and changes of the buffer
        if (sequence_control_set_ptr->encode_context_ptr->sc_buffer > bufferTrshold1 + context_ptr->previous_buffer_check1) {
            encoderModeDelta += +1;
            changeCond = 1;
        }
        else if (sequence_control_set_ptr->encode_context_ptr->sc_buffer > bufferTrshold1 + context_ptr->previous_mode_change_buffer) {
            encoderModeDelta += +1;
            changeCond = 3;
        }

        // Update the encode mode based on the fullness of the buffer
        // If previous ChangeCond was the same, double the threshold2
        if (sequence_control_set_ptr->encode_context_ptr->sc_buffer > bufferTrshold3 &&
            (context_ptr->prev_change_cond != 7 || sequence_control_set_ptr->encode_context_ptr->sc_frame_in > context_ptr->previous_mode_change_frame_in + bufferTrshold2 * 2) &&
            sequence_control_set_ptr->encode_context_ptr->sc_buffer > context_ptr->previous_mode_change_buffer) {
            encoderModeDelta += 1;
            changeCond = 7;
        }
        encoderModeDelta = CLIP3(-1, 1, encoderModeDelta);
        sequence_control_set_ptr->encode_context_ptr->enc_mode = (EbEncMode)CLIP3(1, MAX_ENC_PRESET, (int8_t)sequence_control_set_ptr->encode_context_ptr->enc_mode + encoderModeDelta);

        // Update previous stats
        context_ptr->previous_frame_in_check1 = sequence_control_set_ptr->encode_context_ptr->sc_frame_in;
        context_ptr->previous_buffer_check1 = sequence_control_set_ptr->encode_context_ptr->sc_buffer;

        if (encoderModeDelta) {
            context_ptr->previous_mode_change_buffer = sequence_control_set_ptr->encode_context_ptr->sc_buffer;
            context_ptr->previous_mode_change_frame_in = sequence_control_set_ptr->encode_context_ptr->sc_frame_in;
            context_ptr->prev_enc_mode_delta = encoderModeDelta;
        }
    }

    // Check every bufferTrshold2 for the changes (previous_frame_in_check2 variable)
    if ((sequence_control_set_ptr->encode_context_ptr->sc_frame_in > context_ptr->previous_frame_in_check2 + bufferTrshold2 && sequence_control_set_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)) {
        encoderModeDelta = 0;

        // if no change in the encoder mode and buffer is low enough and level is not increasing, switch to a slower encoder mode
        // If previous ChangeCond was the same, double the threshold2
        if (encoderModeDelta == 0 && sequence_control_set_ptr->encode_context_ptr->sc_frame_in > context_ptr->previous_mode_change_frame_in + bufferTrshold2 &&
            (context_ptr->prev_change_cond != 8 || sequence_control_set_ptr->encode_context_ptr->sc_frame_in > context_ptr->previous_mode_change_frame_in + bufferTrshold2 * 2) &&
            ((sequence_control_set_ptr->encode_context_ptr->sc_buffer - context_ptr->previous_mode_change_buffer < (targetFps / 3)) || context_ptr->previous_mode_change_buffer == 0) &&
            sequence_control_set_ptr->encode_context_ptr->sc_buffer < bufferTrshold3) {
            encoderModeDelta = -1;
            changeCond = 8;
        }

        encoderModeDelta = CLIP3(-1, 1, encoderModeDelta);
        sequence_control_set_ptr->encode_context_ptr->enc_mode = (EbEncMode)CLIP3(1, MAX_ENC_PRESET, (int8_t)sequence_control_set_ptr->encode_context_ptr->enc_mode + encoderModeDelta);

        // Update previous stats
        context_ptr->previous_frame_in_check2 = sequence_control_set_ptr->encode_context_ptr->sc_frame_in;

        if (encoderModeDelta) {
            context_ptr->previous_mode_change_buffer = sequence_control_set_ptr->encode_context_ptr->sc_buffer;
            context_ptr->previous_mode_change_frame_in = sequence_control_set_ptr->encode_context_ptr->sc_frame_in;
            context_ptr->prev_enc_mode_delta = encoderModeDelta;
        }
    }
    // Check every SC_FRAMES_INTERVAL_SPEED frames for the speed calculation (previous_frame_in_check3 variable)
    if (context_ptr->start_flag || (sequence_control_set_ptr->encode_context_ptr->sc_frame_in > context_ptr->previous_frame_in_check3 + SC_FRAMES_INTERVAL_SPEED && sequence_control_set_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)) {
        if (context_ptr->start_flag)
            context_ptr->cur_speed = (uint64_t)(sequence_control_set_ptr->encode_context_ptr->sc_frame_out - 0) * 1000 / (uint64_t)(overallDuration);
        else {
            if (instDuration != 0)
                context_ptr->cur_speed = (uint64_t)(sequence_control_set_ptr->encode_context_ptr->sc_frame_out - context_ptr->prev_frame_out) * 1000 / (uint64_t)(instDuration);
        }
        context_ptr->start_flag = EB_FALSE;

        // Update previous stats
        context_ptr->previous_frame_in_check3 = sequence_control_set_ptr->encode_context_ptr->sc_frame_in;
        context_ptr->prevs_time_seconds = cursTimeSeconds;
        context_ptr->prevs_timeu_seconds = cursTimeuSeconds;
        context_ptr->prev_frame_out = sequence_control_set_ptr->encode_context_ptr->sc_frame_out;
    }
    else if (sequence_control_set_ptr->encode_context_ptr->sc_frame_in < SC_FRAMES_TO_IGNORE && (overallDuration != 0))
        context_ptr->cur_speed = (uint64_t)(sequence_control_set_ptr->encode_context_ptr->sc_frame_out - 0) * 1000 / (uint64_t)(overallDuration);
    if (changeCond)
        context_ptr->prev_change_cond = changeCond;
    sequence_control_set_ptr->encode_context_ptr->sc_frame_in++;
    if (sequence_control_set_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)
        context_ptr->average_enc_mod += sequence_control_set_ptr->encode_context_ptr->enc_mode;
    else
        context_ptr->average_enc_mod = 0;
    // Set the encoder level
    picture_control_set_ptr->enc_mode = sequence_control_set_ptr->encode_context_ptr->enc_mode;

    eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->sc_buffer_mutex);
    context_ptr->prev_enc_mod = sequence_control_set_ptr->encode_context_ptr->enc_mode;
}

void ResetPcsAv1(
    PictureParentControlSet       *picture_control_set_ptr) {
    FrameHeader *frm_hdr = &picture_control_set_ptr->frm_hdr;
    Av1Common *cm = picture_control_set_ptr->av1_cm;

    picture_control_set_ptr->is_skip_mode_allowed = 0;
    picture_control_set_ptr->skip_mode_flag = 0;
    frm_hdr->frame_type                 = KEY_FRAME;
    frm_hdr->show_frame = 1;
    frm_hdr->showable_frame = 1;  // frame can be used as show existing frame in future
    // Flag for a frame used as a reference - not written to the bitstream
    picture_control_set_ptr->is_reference_frame = 0;
    // Flag signaling that the frame is encoded using only INTRA modes.
    picture_control_set_ptr->intra_only = 0;
    // uint8_t last_intra_only;

    frm_hdr->disable_cdf_update = 0;
    frm_hdr->allow_high_precision_mv = 0;
    frm_hdr->force_integer_mv = 0;  // 0 the default in AOM, 1 only integer
    frm_hdr->allow_warped_motion = 0;

    /* profile settings */
#if CONFIG_ENTROPY_STATS
    int32_t coef_cdf_category;
#endif

    frm_hdr->quantization_params.base_q_idx = 31;
    frm_hdr->quantization_params.delta_q_y_dc = 0;
    frm_hdr->quantization_params.delta_q_u_dc = 0;
    frm_hdr->quantization_params.delta_q_v_dc = 0;
    frm_hdr->quantization_params.delta_q_u_ac = 0;
    frm_hdr->quantization_params.delta_q_v_ac = 0;

    picture_control_set_ptr->separate_uv_delta_q = 0;
    // Encoder
    frm_hdr->quantization_params.using_qmatrix = 0;
    frm_hdr->quantization_params.qm_y = 5;
    frm_hdr->quantization_params.qm_u = 5;
    frm_hdr->quantization_params.qm_v = 5;
    // Whether to use previous frame's motion vectors for prediction.
    frm_hdr->use_ref_frame_mvs = 0;
    frm_hdr->is_motion_mode_switchable = 0;
    // Flag signaling how frame contexts should be updated at the end of
    // a frame decode
    picture_control_set_ptr->refresh_frame_context = REFRESH_FRAME_CONTEXT_DISABLED;

    frm_hdr->loop_filter_params.filter_level[0] = 0;
    frm_hdr->loop_filter_params.filter_level[1] = 0;
    frm_hdr->loop_filter_params.filter_level_u = 0;
    frm_hdr->loop_filter_params.filter_level_v = 0;
    frm_hdr->loop_filter_params.sharpness_level = 0;

    frm_hdr->loop_filter_params.mode_ref_delta_enabled = 0;
    frm_hdr->loop_filter_params.mode_ref_delta_update = 0;
    frm_hdr->loop_filter_params.mode_deltas[0] = 0;
    frm_hdr->loop_filter_params.mode_deltas[1] = 0;

    frm_hdr->loop_filter_params.ref_deltas[0] = 1;
    frm_hdr->loop_filter_params.ref_deltas[1] = 0;
    frm_hdr->loop_filter_params.ref_deltas[2] = 0;
    frm_hdr->loop_filter_params.ref_deltas[3] = 0;
    frm_hdr->loop_filter_params.ref_deltas[4] = -1;
    frm_hdr->loop_filter_params.ref_deltas[5] = 0;
    frm_hdr->loop_filter_params.ref_deltas[6] = -1;
    frm_hdr->loop_filter_params.ref_deltas[7] = -1;

    frm_hdr->all_lossless = 0;
    frm_hdr->coded_lossless = 0;
    frm_hdr->reduced_tx_set= 0;
    frm_hdr->reference_mode = SINGLE_REFERENCE;
    picture_control_set_ptr->frame_context_idx = 0; /* Context to use/update */
    for (int32_t i = 0; i < REF_FRAMES; i++)
        picture_control_set_ptr->fb_of_context_type[i] = 0;
    frm_hdr->primary_ref_frame = PRIMARY_REF_NONE;
    picture_control_set_ptr->frame_offset = picture_control_set_ptr->picture_number;
    frm_hdr->error_resilient_mode = 0;
    cm->tiles_info.uniform_tile_spacing_flag = 1;
    picture_control_set_ptr->large_scale_tile = 0;
    picture_control_set_ptr->film_grain_params_present = 0;

    //cdef_pri_damping & cdef_sec_damping are consolidated to cdef_damping
    frm_hdr->CDEF_params.cdef_damping = 0;
    //picture_control_set_ptr->cdef_pri_damping = 0;
    //picture_control_set_ptr->cdef_sec_damping = 0;

    picture_control_set_ptr->nb_cdef_strengths = 1;
    for (int32_t i = 0; i < CDEF_MAX_STRENGTHS; i++) {
        frm_hdr->CDEF_params.cdef_y_strength[i] = 0;
        frm_hdr->CDEF_params.cdef_uv_strength[i] = 0;
    }
    frm_hdr->CDEF_params.cdef_bits = 0;

#if ADD_DELTA_QP_SUPPORT
    frm_hdr->delta_q_params.delta_q_present = 1;
    frm_hdr->delta_lf_params.delta_lf_present = 0;
    frm_hdr->delta_q_params.delta_q_res = DEFAULT_DELTA_Q_RES;
#else
    frm_hdr->delta_q_params.delta_q_present = 0;
#endif

    frm_hdr->delta_lf_params.delta_lf_present = 0;
    frm_hdr->delta_lf_params.delta_lf_res = 0;
    frm_hdr->delta_lf_params.delta_lf_multi = 0;

    frm_hdr->current_frame_id = 0;
    frm_hdr->frame_refs_short_signaling = 0;
    picture_control_set_ptr->allow_comp_inter_inter = 0;
    //  int32_t all_one_sided_refs;
}
/***********************************************
**** Copy the input buffer from the
**** sample application to the library buffers
************************************************/
static EbErrorType copy_frame_buffer(
    SequenceControlSet            *sequence_control_set_ptr,
    uint8_t                          *dst,
    uint8_t                          *src)
{
    EbSvtAv1EncConfiguration        *config = &sequence_control_set_ptr->static_config;
    EbErrorType                      return_error = EB_ErrorNone;

    EbPictureBufferDesc             *dst_picture_ptr = (EbPictureBufferDesc*)dst;
    EbPictureBufferDesc             *src_picture_ptr = (EbPictureBufferDesc*)src;
    uint16_t                         inputRowIndex;
    EbBool                           is16BitInput = (EbBool)(config->encoder_bit_depth > EB_8BIT);

    // Need to include for Interlacing on the fly with pictureScanType = 1

    if (!is16BitInput) {
        uint32_t     lumaBufferOffset = (dst_picture_ptr->stride_y*sequence_control_set_ptr->top_padding + sequence_control_set_ptr->left_padding) << is16BitInput;
        uint32_t     chromaBufferOffset = (dst_picture_ptr->stride_cr*(sequence_control_set_ptr->top_padding >> 1) + (sequence_control_set_ptr->left_padding >> 1)) << is16BitInput;
        uint16_t     lumaStride = dst_picture_ptr->stride_y << is16BitInput;
        uint16_t     chromaStride = dst_picture_ptr->stride_cb << is16BitInput;
        uint16_t     lumaWidth = (uint16_t)(dst_picture_ptr->width - sequence_control_set_ptr->max_input_pad_right) << is16BitInput;
        uint16_t     chromaWidth = (lumaWidth >> 1) << is16BitInput;
        uint16_t     lumaHeight = (uint16_t)(dst_picture_ptr->height - sequence_control_set_ptr->max_input_pad_bottom);

        //uint16_t     lumaHeight  = input_picture_ptr->max_height;
        // Y
        for (inputRowIndex = 0; inputRowIndex < lumaHeight; inputRowIndex++) {
            EB_MEMCPY((dst_picture_ptr->buffer_y + lumaBufferOffset + lumaStride * inputRowIndex),
                (src_picture_ptr->buffer_y + lumaBufferOffset + lumaStride * inputRowIndex),
                lumaWidth);
        }

        // U
        for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
            EB_MEMCPY((dst_picture_ptr->buffer_cb + chromaBufferOffset + chromaStride * inputRowIndex),
                (src_picture_ptr->buffer_cb + chromaBufferOffset + chromaStride * inputRowIndex),
                chromaWidth);
        }

        // V
        for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
            EB_MEMCPY((dst_picture_ptr->buffer_cr + chromaBufferOffset + chromaStride * inputRowIndex),
                (src_picture_ptr->buffer_cr + chromaBufferOffset + chromaStride * inputRowIndex),
                chromaWidth);
        }
    }
    else if (is16BitInput && config->compressed_ten_bit_format == 1)
    {
        {
            uint32_t  lumaBufferOffset = (dst_picture_ptr->stride_y*sequence_control_set_ptr->top_padding + sequence_control_set_ptr->left_padding);
            uint32_t  chromaBufferOffset = (dst_picture_ptr->stride_cr*(sequence_control_set_ptr->top_padding >> 1) + (sequence_control_set_ptr->left_padding >> 1));
            uint16_t  lumaStride = dst_picture_ptr->stride_y;
            uint16_t  chromaStride = dst_picture_ptr->stride_cb;
            uint16_t  lumaWidth = (uint16_t)(dst_picture_ptr->width - sequence_control_set_ptr->max_input_pad_right);
            uint16_t  chromaWidth = (lumaWidth >> 1);
            uint16_t  lumaHeight = (uint16_t)(dst_picture_ptr->height - sequence_control_set_ptr->max_input_pad_bottom);

            // Y 8bit
            for (inputRowIndex = 0; inputRowIndex < lumaHeight; inputRowIndex++) {
                EB_MEMCPY((dst_picture_ptr->buffer_y + lumaBufferOffset + lumaStride * inputRowIndex),
                    (src_picture_ptr->buffer_y + lumaBufferOffset + lumaStride * inputRowIndex),
                    lumaWidth);
            }

            // U 8bit
            for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
                EB_MEMCPY((dst_picture_ptr->buffer_cb + chromaBufferOffset + chromaStride * inputRowIndex),
                    (src_picture_ptr->buffer_cb + chromaBufferOffset + chromaStride * inputRowIndex),
                    chromaWidth);
            }

            // V 8bit
            for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
                EB_MEMCPY((dst_picture_ptr->buffer_cr + chromaBufferOffset + chromaStride * inputRowIndex),
                    (src_picture_ptr->buffer_cr + chromaBufferOffset + chromaStride * inputRowIndex),
                    chromaWidth);
            }
            // AMIR to update
            ////efficient copy - final
            ////compressed 2Bit in 1D format
            //{
            //    uint16_t luma2BitWidth = sequence_control_set_ptr->max_input_luma_width / 4;
            //    uint16_t lumaHeight = sequence_control_set_ptr->max_input_luma_height;

            //    uint16_t sourceLuma2BitStride = sourceLumaStride / 4;
            //    uint16_t sourceChroma2BitStride = sourceLuma2BitStride >> 1;

            //    for (inputRowIndex = 0; inputRowIndex < lumaHeight; inputRowIndex++) {
            //        EB_MEMCPY(input_picture_ptr->buffer_bit_inc_y + luma2BitWidth * inputRowIndex, inputPtr->luma_ext + sourceLuma2BitStride * inputRowIndex, luma2BitWidth);
            //    }
            //    for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
            //        EB_MEMCPY(input_picture_ptr->buffer_bit_inc_cb + (luma2BitWidth >> 1)*inputRowIndex, inputPtr->cb_ext + sourceChroma2BitStride * inputRowIndex, luma2BitWidth >> 1);
            //    }
            //    for (inputRowIndex = 0; inputRowIndex < lumaHeight >> 1; inputRowIndex++) {
            //        EB_MEMCPY(input_picture_ptr->buffer_bit_inc_cr + (luma2BitWidth >> 1)*inputRowIndex, inputPtr->cr_ext + sourceChroma2BitStride * inputRowIndex, luma2BitWidth >> 1);
            //    }
            //}
        }
    }
    else { // 10bit packed

        EB_MEMCPY(dst_picture_ptr->buffer_y,
            src_picture_ptr->buffer_y ,
            src_picture_ptr->luma_size);

        EB_MEMCPY(dst_picture_ptr->buffer_cb,
            src_picture_ptr->buffer_cb,
            src_picture_ptr->chroma_size);

        EB_MEMCPY(dst_picture_ptr->buffer_cr,
            src_picture_ptr->buffer_cr,
            src_picture_ptr->chroma_size);

        EB_MEMCPY(dst_picture_ptr->buffer_bit_inc_y,
            src_picture_ptr->buffer_bit_inc_y,
            src_picture_ptr->luma_size);

        EB_MEMCPY(dst_picture_ptr->buffer_bit_inc_cb,
            src_picture_ptr->buffer_bit_inc_cb,
            src_picture_ptr->chroma_size);

        EB_MEMCPY(dst_picture_ptr->buffer_bit_inc_cr,
            src_picture_ptr->buffer_bit_inc_cr,
            src_picture_ptr->chroma_size);
    }
    return return_error;
}
static void CopyInputBuffer(
    SequenceControlSet*    sequenceControlSet,
    EbBufferHeaderType*     dst,
    EbBufferHeaderType*     src
)
{
    // Copy the higher level structure
    dst->n_alloc_len = src->n_alloc_len;
    dst->n_filled_len = src->n_filled_len;
    dst->flags = src->flags;
    dst->pts = src->pts;
    dst->n_tick_count = src->n_tick_count;
    dst->size = src->size;
    dst->qp = src->qp;
    dst->pic_type = src->pic_type;

    // Copy the picture buffer
    if (src->p_buffer != NULL)
        copy_frame_buffer(sequenceControlSet, dst->p_buffer, src->p_buffer);
}
/***************************************
 * ResourceCoordination Kernel
 ***************************************/
void* resource_coordination_kernel(void *input_ptr)
{
    ResourceCoordinationContext   *context_ptr = (ResourceCoordinationContext*)input_ptr;

    EbObjectWrapper               *picture_control_set_wrapper_ptr;

    PictureParentControlSet       *picture_control_set_ptr;

    EbObjectWrapper               *previousSequenceControlSetWrapperPtr;
    SequenceControlSet            *sequence_control_set_ptr;

    EbObjectWrapper               *ebInputWrapperPtr;
    EbBufferHeaderType              *ebInputPtr;
    EbObjectWrapper               *outputWrapperPtr;
    ResourceCoordinationResults   *outputResultsPtr;

    EbObjectWrapper               *input_picture_wrapper_ptr;
    EbObjectWrapper               *reference_picture_wrapper_ptr;

    uint32_t                         instance_index;
    EbBool                           end_of_sequence_flag = EB_FALSE;
    uint32_t                         aspectRatio;

    uint32_t                         input_size = 0;
    EbObjectWrapper               *prevPictureControlSetWrapperPtr = 0;

    for (;;) {
        // Tie instance_index to zero for now...
        instance_index = 0;

        // Get the Next svt Input Buffer [BLOCKING]
        eb_get_full_object(
            context_ptr->input_buffer_fifo_ptr,
            &ebInputWrapperPtr);
        ebInputPtr = (EbBufferHeaderType*)ebInputWrapperPtr->object_ptr;
        sequence_control_set_ptr = context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr;

        // If config changes occured since the last picture began encoding, then
        //   prepare a new sequence_control_set_ptr containing the new changes and update the state
        //   of the previous Active SequenceControlSet
        eb_block_on_mutex(context_ptr->sequence_control_set_instance_array[instance_index]->config_mutex);
        if (context_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->initial_picture) {
            // Update picture width, picture height, cropping right offset, cropping bottom offset, and conformance windows
            if (context_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->initial_picture)

            {
                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->seq_header.max_frame_width = context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width;
                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->seq_header.max_frame_height = context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height;
                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->chroma_width = (context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_width >> 1);
                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->chroma_height = (context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_luma_height >> 1);

                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pad_right = context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_pad_right;
                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->cropping_right_offset = context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pad_right;
                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pad_bottom = context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->max_input_pad_bottom;
                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->cropping_bottom_offset = context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pad_bottom;

                if (context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pad_right != 0 || context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->pad_bottom != 0)
                    context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->conformance_window_flag = 1;
                else
                    context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->conformance_window_flag = 0;
                input_size = context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->seq_header.max_frame_width * context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr->seq_header.max_frame_height;
            }

            // Copy previous Active SequenceControlSetPtr to a place holder
            previousSequenceControlSetWrapperPtr = context_ptr->sequenceControlSetActiveArray[instance_index];

            // Get empty SequenceControlSet [BLOCKING]
            eb_get_empty_object(
                context_ptr->sequence_control_set_empty_fifo_ptr,
                &context_ptr->sequenceControlSetActiveArray[instance_index]);

            // Copy the contents of the active SequenceControlSet into the new empty SequenceControlSet
            copy_sequence_control_set(
                (SequenceControlSet*)context_ptr->sequenceControlSetActiveArray[instance_index]->object_ptr,
                context_ptr->sequence_control_set_instance_array[instance_index]->sequence_control_set_ptr);

            // Disable releaseFlag of new SequenceControlSet
            eb_object_release_disable(
                context_ptr->sequenceControlSetActiveArray[instance_index]);

            if (previousSequenceControlSetWrapperPtr != EB_NULL) {
                // Enable releaseFlag of old SequenceControlSet
                eb_object_release_enable(
                    previousSequenceControlSetWrapperPtr);

                // Check to see if previous SequenceControlSet is already inactive, if TRUE then release the SequenceControlSet
                if (previousSequenceControlSetWrapperPtr->live_count == 0) {
                    eb_release_object(
                        previousSequenceControlSetWrapperPtr);
                }
            }
        }
        eb_release_mutex(context_ptr->sequence_control_set_instance_array[instance_index]->config_mutex);
#if ENABLE_CDF_UPDATE
        // Seque Control Set is released by Rate Control after passing through MDC->MD->ENCDEC->Packetization->RateControl,
        // in the PictureManager after receiving the reference and in PictureManager after receiving the feedback
        eb_object_inc_live_count(
            context_ptr->sequenceControlSetActiveArray[instance_index],
            3);
#else
        // Sequence Control Set is released by Rate Control after passing through MDC->MD->ENCDEC->Packetization->RateControl
        //   and in the PictureManager
        eb_object_inc_live_count(
            context_ptr->sequenceControlSetActiveArray[instance_index],
            2);
#endif
        // Set the current SequenceControlSet
        sequence_control_set_ptr = (SequenceControlSet*)context_ptr->sequenceControlSetActiveArray[instance_index]->object_ptr;

        // Init SB Params
        if (context_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->initial_picture) {
            derive_input_resolution(
                sequence_control_set_ptr,
                input_size);

            sb_params_init(sequence_control_set_ptr);
            sb_geom_init(sequence_control_set_ptr);

            sequence_control_set_ptr->enable_altrefs =  sequence_control_set_ptr->static_config.enable_altrefs &&
                sequence_control_set_ptr->static_config.altref_nframes > 1 &&
                sequence_control_set_ptr->static_config.encoder_bit_depth == EB_8BIT ? EB_TRUE : EB_FALSE;

            // Sep PM mode
            sequence_control_set_ptr->pm_mode = sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE ?
                PM_MODE_2 :
                PM_MODE_1;

            // Construct PM Trans Coeff Shaping
            if (context_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->initial_picture) {
                if (sequence_control_set_ptr->pm_mode == PM_MODE_0)
                    construct_pm_trans_coeff_shaping(sequence_control_set_ptr);
            }
        }
        // Since at this stage we do not know the prediction structure and the location of ALT_REF pictures,
        // for every picture (except first picture), we allocate two: 1. original picture, 2. potential Overlay picture.
        // In Picture Decision Process, where the overlay frames are known, they extra pictures are released
        uint8_t has_overlay = (sequence_control_set_ptr->static_config.enable_overlays == EB_FALSE ||
            context_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->initial_picture) ? 0 : 1;
        for (uint8_t loop_index = 0; loop_index <= has_overlay && !end_of_sequence_flag; loop_index++) {
            //Get a New ParentPCS where we will hold the new inputPicture
            eb_get_empty_object(
                context_ptr->picture_control_set_fifo_ptr_array[instance_index],
                &picture_control_set_wrapper_ptr);

            // Parent PCS is released by the Rate Control after passing through MDC->MD->ENCDEC->Packetization
            eb_object_inc_live_count(
                picture_control_set_wrapper_ptr,
                1);

            picture_control_set_ptr = (PictureParentControlSet*)picture_control_set_wrapper_ptr->object_ptr;

            picture_control_set_ptr->p_pcs_wrapper_ptr = picture_control_set_wrapper_ptr;

            picture_control_set_ptr->overlay_ppcs_ptr = NULL;
            picture_control_set_ptr->is_alt_ref       = 0;
            if (loop_index) {
                picture_control_set_ptr->is_overlay = 1;
                // set the overlay_ppcs_ptr in the original (ALT_REF) ppcs to the current ppcs
                EbObjectWrapper               *alt_ref_picture_control_set_wrapper_ptr = (context_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->initial_picture) ?
                    picture_control_set_wrapper_ptr :
                    sequence_control_set_ptr->encode_context_ptr->previous_picture_control_set_wrapper_ptr;

                picture_control_set_ptr->alt_ref_ppcs_ptr = ((PictureParentControlSet*)alt_ref_picture_control_set_wrapper_ptr->object_ptr);
                picture_control_set_ptr->alt_ref_ppcs_ptr->overlay_ppcs_ptr = picture_control_set_ptr;
            }
            else {
                picture_control_set_ptr->is_overlay = 0;
                picture_control_set_ptr->alt_ref_ppcs_ptr = NULL;
            }
            // Set the Encoder mode
            picture_control_set_ptr->enc_mode = sequence_control_set_ptr->static_config.enc_mode;

            // Keep track of the previous input for the ZZ SADs computation
            picture_control_set_ptr->previous_picture_control_set_wrapper_ptr = (context_ptr->sequence_control_set_instance_array[instance_index]->encode_context_ptr->initial_picture) ?
                picture_control_set_wrapper_ptr :
                sequence_control_set_ptr->encode_context_ptr->previous_picture_control_set_wrapper_ptr;
            if (loop_index == 0)
                sequence_control_set_ptr->encode_context_ptr->previous_picture_control_set_wrapper_ptr = picture_control_set_wrapper_ptr;
            // Copy data from the svt buffer to the input frame
            // *Note - Assumes 4:2:0 planar
            input_picture_wrapper_ptr = ebInputWrapperPtr;
            picture_control_set_ptr->enhanced_picture_ptr = (EbPictureBufferDesc*)ebInputPtr->p_buffer;
            picture_control_set_ptr->input_ptr = ebInputPtr;
            end_of_sequence_flag = (picture_control_set_ptr->input_ptr->flags & EB_BUFFERFLAG_EOS) ? EB_TRUE : EB_FALSE;
            EbStartTime(&picture_control_set_ptr->start_time_seconds, &picture_control_set_ptr->start_time_u_seconds);

            picture_control_set_ptr->sequence_control_set_wrapper_ptr = context_ptr->sequenceControlSetActiveArray[instance_index];
            picture_control_set_ptr->sequence_control_set_ptr = sequence_control_set_ptr;
            picture_control_set_ptr->input_picture_wrapper_ptr = input_picture_wrapper_ptr;
            picture_control_set_ptr->end_of_sequence_flag = end_of_sequence_flag;

            if (loop_index == 1) {
                // Get a new input picture for overlay.
                EbObjectWrapper     *input_pic_wrapper_ptr;

                // Get a new input picture for overlay.
                eb_get_empty_object(
                    sequence_control_set_ptr->encode_context_ptr->overlay_input_picture_pool_fifo_ptr,
                    &input_pic_wrapper_ptr);

                // Copy from original picture (picture_control_set_ptr->input_picture_wrapper_ptr), which is shared between overlay and alt_ref up to this point, to the new input picture.
                if (picture_control_set_ptr->alt_ref_ppcs_ptr->input_picture_wrapper_ptr->object_ptr != NULL) {
                    CopyInputBuffer(
                        sequence_control_set_ptr,
                        (EbBufferHeaderType*)input_pic_wrapper_ptr->object_ptr,
                        (EbBufferHeaderType*)picture_control_set_ptr->alt_ref_ppcs_ptr->input_picture_wrapper_ptr->object_ptr);
                }
                // Assign the new picture to the new pointers
                picture_control_set_ptr->input_ptr = (EbBufferHeaderType*)input_pic_wrapper_ptr->object_ptr;
                picture_control_set_ptr->enhanced_picture_ptr = (EbPictureBufferDesc*)picture_control_set_ptr->input_ptr->p_buffer;
                picture_control_set_ptr->input_picture_wrapper_ptr = input_pic_wrapper_ptr;
            }
            // Set Picture Control Flags
            picture_control_set_ptr->idr_flag = sequence_control_set_ptr->encode_context_ptr->initial_picture || (picture_control_set_ptr->input_ptr->pic_type == EB_AV1_KEY_PICTURE);
            picture_control_set_ptr->cra_flag = (picture_control_set_ptr->input_ptr->pic_type == EB_AV1_INTRA_ONLY_PICTURE) ? EB_TRUE : EB_FALSE;
            picture_control_set_ptr->scene_change_flag = EB_FALSE;
            picture_control_set_ptr->qp_on_the_fly = EB_FALSE;
            picture_control_set_ptr->sb_total_count = sequence_control_set_ptr->sb_total_count;
            picture_control_set_ptr->eos_coming = (ebInputPtr->flags & (EB_BUFFERFLAG_EOS << 1)) ? EB_TRUE : EB_FALSE;

            if (sequence_control_set_ptr->static_config.speed_control_flag) {
                SpeedBufferControl(
                    context_ptr,
                    picture_control_set_ptr,
                    sequence_control_set_ptr);
            }
            else
                picture_control_set_ptr->enc_mode = (EbEncMode)sequence_control_set_ptr->static_config.enc_mode;
            aspectRatio = (sequence_control_set_ptr->seq_header.max_frame_width * 10) / sequence_control_set_ptr->seq_header.max_frame_height;
            aspectRatio = (aspectRatio <= ASPECT_RATIO_4_3) ? ASPECT_RATIO_CLASS_0 : (aspectRatio <= ASPECT_RATIO_16_9) ? ASPECT_RATIO_CLASS_1 : ASPECT_RATIO_CLASS_2;

            // Set the SCD Mode
            sequence_control_set_ptr->scd_mode = sequence_control_set_ptr->static_config.scene_change_detection == 0 ?
                SCD_MODE_0 :
                SCD_MODE_1;

            // Set the block mean calculation prec
            sequence_control_set_ptr->block_mean_calc_prec = BLOCK_MEAN_PREC_SUB;

            // Pre-Analysis Signal(s) derivation
            signal_derivation_pre_analysis_oq(
                sequence_control_set_ptr,
                picture_control_set_ptr);
#if QPS_TUNING
            picture_control_set_ptr->filtered_sse = 0;
            picture_control_set_ptr->filtered_sse_uv = 0;
#endif
            // Rate Control
            // Set the ME Distortion and OIS Historgrams to zero
            if (sequence_control_set_ptr->static_config.rate_control_mode) {
                EB_MEMSET(picture_control_set_ptr->me_distortion_histogram, 0, NUMBER_OF_SAD_INTERVALS * sizeof(uint16_t));
                EB_MEMSET(picture_control_set_ptr->ois_distortion_histogram, 0, NUMBER_OF_INTRA_SAD_INTERVALS * sizeof(uint16_t));
            }
            picture_control_set_ptr->full_sb_count = 0;

            if (sequence_control_set_ptr->static_config.use_qp_file == 1) {
                picture_control_set_ptr->qp_on_the_fly = EB_TRUE;
                if (picture_control_set_ptr->input_ptr->qp > MAX_QP_VALUE) {
                    SVT_LOG("SVT [WARNING]: INPUT QP OUTSIDE OF RANGE\n");
                    picture_control_set_ptr->qp_on_the_fly = EB_FALSE;
                    picture_control_set_ptr->picture_qp = (uint8_t)sequence_control_set_ptr->qp;
                }
                picture_control_set_ptr->picture_qp = (uint8_t)picture_control_set_ptr->input_ptr->qp;
            }
            else {
                picture_control_set_ptr->qp_on_the_fly = EB_FALSE;
                picture_control_set_ptr->picture_qp = (uint8_t)sequence_control_set_ptr->qp;
            }

            // Picture Stats
            if (loop_index == has_overlay || end_of_sequence_flag)
                picture_control_set_ptr->picture_number = context_ptr->picture_number_array[instance_index]++;
            else
                picture_control_set_ptr->picture_number = context_ptr->picture_number_array[instance_index];
            ResetPcsAv1(picture_control_set_ptr);

            sequence_control_set_ptr->encode_context_ptr->initial_picture = EB_FALSE;

            // Get Empty Reference Picture Object
            eb_get_empty_object(
                sequence_control_set_ptr->encode_context_ptr->pa_reference_picture_pool_fifo_ptr,
                &reference_picture_wrapper_ptr);

            picture_control_set_ptr->pa_reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
            // Since overlay pictures are not added to PA_Reference queue in PD and not released there, the life count is only set to 1
            if (picture_control_set_ptr->is_overlay)
                // Give the new Reference a nominal live_count of 1
                eb_object_inc_live_count(
                    picture_control_set_ptr->pa_reference_picture_wrapper_ptr,
                    1);
            else
                eb_object_inc_live_count(
                    picture_control_set_ptr->pa_reference_picture_wrapper_ptr,
                    2);
            ((EbPaReferenceObject*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->object_ptr)->input_padded_picture_ptr->buffer_y = picture_control_set_ptr->enhanced_picture_ptr->buffer_y;

            // Get Empty Output Results Object
            if (picture_control_set_ptr->picture_number > 0 && (prevPictureControlSetWrapperPtr != NULL))
            {
                ((PictureParentControlSet       *)prevPictureControlSetWrapperPtr->object_ptr)->end_of_sequence_flag = end_of_sequence_flag;
                eb_get_empty_object(
                    context_ptr->resource_coordination_results_output_fifo_ptr,
                    &outputWrapperPtr);
                outputResultsPtr = (ResourceCoordinationResults*)outputWrapperPtr->object_ptr;
                outputResultsPtr->picture_control_set_wrapper_ptr = prevPictureControlSetWrapperPtr;
                // since overlay frame has the end of sequence set properly, set the end of sequence to true in the alt ref picture
                if (((PictureParentControlSet       *)prevPictureControlSetWrapperPtr->object_ptr)->is_overlay && end_of_sequence_flag)
                    ((PictureParentControlSet       *)prevPictureControlSetWrapperPtr->object_ptr)->alt_ref_ppcs_ptr->end_of_sequence_flag = EB_TRUE;
                // Post the finished Results Object
                eb_post_full_object(outputWrapperPtr);
            }
            prevPictureControlSetWrapperPtr = picture_control_set_wrapper_ptr;
        }
    }

    return EB_NULL;
}
