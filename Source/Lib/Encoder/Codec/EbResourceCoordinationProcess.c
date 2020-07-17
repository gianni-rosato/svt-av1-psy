/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbEncHandle.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureBufferDesc.h"
#include "EbResourceCoordinationProcess.h"
#include "EbResourceCoordinationResults.h"
#include "EbTransforms.h"
#include "EbTime.h"
#include "EbObject.h"
#include "EbLog.h"
#include "common_dsp_rtcd.h"
typedef struct ResourceCoordinationContext {
    EbFifo *                       input_buffer_fifo_ptr;
    EbFifo *                       resource_coordination_results_output_fifo_ptr;
    EbFifo **                      picture_control_set_fifo_ptr_array;
    EbSequenceControlSetInstance **scs_instance_array;
    EbObjectWrapper **             sequence_control_set_active_array;
    EbFifo *                       sequence_control_set_empty_fifo_ptr;
    EbCallback **                  app_callback_ptr_array;

    // Compute Segments
    uint32_t compute_segments_total_count_array;
    uint32_t encode_instances_total_count;

    // Picture Number Array
    uint64_t *picture_number_array;

    uint64_t average_enc_mod;
    uint8_t  prev_enc_mod;
    int8_t   prev_enc_mode_delta;
    uint8_t  prev_change_cond;

    int64_t previous_mode_change_buffer;
    int64_t previous_mode_change_frame_in;
    int64_t previous_buffer_check1;
    int64_t previous_frame_in_check1;
    int64_t previous_frame_in_check2;
    int64_t previous_frame_in_check3;

    uint64_t cur_speed; // speed x 1000
    uint64_t prevs_time_seconds;
    uint64_t prevs_timeu_seconds;
    int64_t  prev_frame_out;

    uint64_t first_in_pic_arrived_time_seconds;
    uint64_t first_in_pic_arrived_timeu_seconds;
    EbBool   start_flag;
} ResourceCoordinationContext;

static void resource_coordination_context_dctor(EbPtr p) {
    EbThreadContext *thread_contxt_ptr = (EbThreadContext *)p;
    if (thread_contxt_ptr->priv) {
        ResourceCoordinationContext *obj = (ResourceCoordinationContext *)thread_contxt_ptr->priv;

        EB_FREE_ARRAY(obj->sequence_control_set_active_array);
        EB_FREE_ARRAY(obj->picture_number_array);
        EB_FREE_ARRAY(obj->picture_control_set_fifo_ptr_array);
        EB_FREE_ARRAY(obj);
    }
}

/************************************************
 * Resource Coordination Context Constructor
 ************************************************/
EbErrorType resource_coordination_context_ctor(EbThreadContext *thread_contxt_ptr,
                                               EbEncHandle *    enc_handle_ptr) {
    ResourceCoordinationContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_contxt_ptr->priv  = context_ptr;
    thread_contxt_ptr->dctor = resource_coordination_context_dctor;

    EB_MALLOC_ARRAY(context_ptr->picture_control_set_fifo_ptr_array,
                    enc_handle_ptr->encode_instance_total_count);
    for (uint32_t i = 0; i < enc_handle_ptr->encode_instance_total_count; i++) {
        //ResourceCoordination works with ParentPCS
        context_ptr->picture_control_set_fifo_ptr_array[i] = eb_system_resource_get_producer_fifo(
            enc_handle_ptr->picture_parent_control_set_pool_ptr_array[i], 0);
    }

    context_ptr->input_buffer_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->input_buffer_resource_ptr, 0);
    context_ptr->resource_coordination_results_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(
            enc_handle_ptr->resource_coordination_results_resource_ptr, 0);
    context_ptr->scs_instance_array = enc_handle_ptr->scs_instance_array;
    context_ptr->sequence_control_set_empty_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->scs_pool_ptr, 0);
    context_ptr->app_callback_ptr_array = enc_handle_ptr->app_callback_ptr_array;
    context_ptr->compute_segments_total_count_array =
        enc_handle_ptr->compute_segments_total_count_array;
    context_ptr->encode_instances_total_count = enc_handle_ptr->encode_instance_total_count;

    // Allocate SequenceControlSetActiveArray
    EB_CALLOC_ARRAY(context_ptr->sequence_control_set_active_array,
                    context_ptr->encode_instances_total_count);

    EB_CALLOC_ARRAY(context_ptr->picture_number_array, context_ptr->encode_instances_total_count);

    context_ptr->average_enc_mod                    = 0;
    context_ptr->prev_enc_mod                       = 0;
    context_ptr->prev_enc_mode_delta                = 0;
    context_ptr->cur_speed                          = 0; // speed x 1000
    context_ptr->previous_mode_change_buffer        = 0;
    context_ptr->first_in_pic_arrived_time_seconds  = 0;
    context_ptr->first_in_pic_arrived_timeu_seconds = 0;
    context_ptr->previous_frame_in_check1           = 0;
    context_ptr->previous_frame_in_check2           = 0;
    context_ptr->previous_frame_in_check3           = 0;
    context_ptr->previous_mode_change_frame_in      = 0;
    context_ptr->prevs_time_seconds                 = 0;
    context_ptr->prevs_timeu_seconds                = 0;
    context_ptr->prev_frame_out                     = 0;
    context_ptr->start_flag                         = EB_FALSE;

    context_ptr->previous_buffer_check1 = 0;
    context_ptr->prev_change_cond       = 0;

    return EB_ErrorNone;
}

/******************************************************
* Derive Pre-Analysis settings for OQ
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
EbErrorType signal_derivation_pre_analysis_oq(SequenceControlSet *     scs_ptr,
                                              PictureParentControlSet *pcs_ptr) {
    EbErrorType return_error     = EB_ErrorNone;
#if !REFACTOR_ME_HME
    uint8_t     input_resolution = scs_ptr->input_resolution;

    // HME Flags updated @ signal_derivation_multi_processes_oq
    uint8_t hme_me_level =
        scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#endif
    // Derive HME Flag
    if (scs_ptr->static_config.use_default_me_hme) {
#if REFACTOR_ME_HME
        // Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
        // Will be later updated for SC/NSC in PictureDecisionProcess
        pcs_ptr->enable_hme_flag        = 1;
        pcs_ptr->enable_hme_level0_flag = 1;
        pcs_ptr->enable_hme_level1_flag = 1;
        pcs_ptr->enable_hme_level2_flag = 1;
#else
        pcs_ptr->enable_hme_flag = enable_hme_flag[0][input_resolution][hme_me_level] ||
                                   enable_hme_flag[1][input_resolution][hme_me_level];
        pcs_ptr->enable_hme_level0_flag =
            enable_hme_level0_flag[0][input_resolution][hme_me_level] ||
            enable_hme_level0_flag[1][input_resolution][hme_me_level];
        pcs_ptr->enable_hme_level1_flag =
            enable_hme_level1_flag[0][input_resolution][hme_me_level] ||
            enable_hme_level1_flag[1][input_resolution][hme_me_level];
        pcs_ptr->enable_hme_level2_flag =
            enable_hme_level2_flag[0][input_resolution][hme_me_level] ||
            enable_hme_level2_flag[1][input_resolution][hme_me_level];
#endif
    } else {
        pcs_ptr->enable_hme_flag        = scs_ptr->static_config.enable_hme_flag;
        pcs_ptr->enable_hme_level0_flag = scs_ptr->static_config.enable_hme_level0_flag;
        pcs_ptr->enable_hme_level1_flag = scs_ptr->static_config.enable_hme_level1_flag;
        pcs_ptr->enable_hme_level2_flag = scs_ptr->static_config.enable_hme_level2_flag;
    }
#if REFACTOR_ME_HME
    // Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
    // Will be later updated for SC/NSC in PictureDecisionProcess
    pcs_ptr->tf_enable_hme_flag        = 1;
    pcs_ptr->tf_enable_hme_level0_flag = 1;
    pcs_ptr->tf_enable_hme_level1_flag = 1;
    pcs_ptr->tf_enable_hme_level2_flag = 1;
#else
    pcs_ptr->tf_enable_hme_flag = tf_enable_hme_flag[0][input_resolution][hme_me_level] ||
                                  tf_enable_hme_flag[1][input_resolution][hme_me_level];
    pcs_ptr->tf_enable_hme_level0_flag =
        tf_enable_hme_level0_flag[0][input_resolution][hme_me_level] ||
        tf_enable_hme_level0_flag[1][input_resolution][hme_me_level];
    pcs_ptr->tf_enable_hme_level1_flag =
        tf_enable_hme_level1_flag[0][input_resolution][hme_me_level] ||
        tf_enable_hme_level1_flag[1][input_resolution][hme_me_level];
    pcs_ptr->tf_enable_hme_level2_flag =
        tf_enable_hme_level2_flag[0][input_resolution][hme_me_level] ||
        tf_enable_hme_level2_flag[1][input_resolution][hme_me_level];
#endif
    if (scs_ptr->static_config.enable_intra_edge_filter == DEFAULT)
        scs_ptr->seq_header.enable_intra_edge_filter = 1;
    else
        scs_ptr->seq_header.enable_intra_edge_filter = (uint8_t)scs_ptr->static_config.enable_intra_edge_filter;

    if (scs_ptr->static_config.pic_based_rate_est == DEFAULT)
#if PIC_BASED_RE_OFF
        scs_ptr->seq_header.pic_based_rate_est = 0;
#else
        scs_ptr->seq_header.pic_based_rate_est = 1;
#endif
    else
        scs_ptr->seq_header.pic_based_rate_est = (uint8_t)scs_ptr->static_config.pic_based_rate_est;

    if (scs_ptr->static_config.enable_restoration_filtering == DEFAULT) {
#if !MAR10_ADOPTIONS
        if (pcs_ptr->enc_mode >= ENC_M8)
            scs_ptr->seq_header.enable_restoration = 0;
        else
#endif
#if M8_RESTORATION && !UPGRADE_M8
            scs_ptr->seq_header.enable_restoration = (pcs_ptr->enc_mode <= ENC_M5) ? 1 : 0;
#else
#if REVERT_BLUE
#if JUNE17_ADOPTIONS
            scs_ptr->seq_header.enable_restoration = (pcs_ptr->enc_mode <= ENC_M6) ? 1 : 0;
#else
#if PRESET_SHIFITNG
            scs_ptr->seq_header.enable_restoration = (pcs_ptr->enc_mode <= ENC_M5) ? 1 : 0;
#else
            scs_ptr->seq_header.enable_restoration = (pcs_ptr->enc_mode <= ENC_M7) ? 1 : 0;
#endif
#endif
#else
            scs_ptr->seq_header.enable_restoration = 1;
#endif
#endif
    } else
        scs_ptr->seq_header.enable_restoration =
            (uint8_t)scs_ptr->static_config.enable_restoration_filtering;

    if (scs_ptr->static_config.cdef_mode == DEFAULT)
        scs_ptr->seq_header.enable_cdef = 1;
    else
        scs_ptr->seq_header.enable_cdef = (uint8_t)(scs_ptr->static_config.cdef_mode>0);

#if !M8_CDF
#if MAR2_M7_ADOPTIONS
#if MAR10_ADOPTIONS
    scs_ptr->cdf_mode = (pcs_ptr->enc_mode <= ENC_M8) ? 0 : 1;
#else
    scs_ptr->cdf_mode = (pcs_ptr->enc_mode <= ENC_M7) ? 0 : 1;
#endif
#else
    scs_ptr->cdf_mode = (pcs_ptr->enc_mode <= ENC_M6) ? 0 : 1;
#endif
#endif
    if (scs_ptr->static_config.enable_warped_motion == DEFAULT) {
        scs_ptr->seq_header.enable_warped_motion = 1;
    } else
        scs_ptr->seq_header.enable_warped_motion = (uint8_t)scs_ptr->static_config.enable_warped_motion;

    return return_error;
}

//******************************************************************************//
// Modify the Enc mode based on the buffer Status
// Inputs: TargetSpeed, Status of the SCbuffer
// Output: EncMod
//******************************************************************************//
void speed_buffer_control(ResourceCoordinationContext *context_ptr,
                          PictureParentControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    uint64_t curs_time_seconds  = 0;
    uint64_t curs_time_useconds = 0;
    double   overall_duration   = 0.0;
    double   inst_duration      = 0.0;
    int8_t   encoder_mode_delta = 0;
    int64_t  input_frames_count = 0;
    int8_t   change_cond        = 0;
    int64_t  target_fps         = (scs_ptr->static_config.injector_frame_rate >> 16);

    int64_t buffer_threshold_1 = SC_FRAMES_INTERVAL_T1;
    int64_t buffer_threshold_2 = SC_FRAMES_INTERVAL_T2;
    int64_t buffer_threshold_3 = MIN(target_fps * 3, SC_FRAMES_INTERVAL_T3);
    eb_block_on_mutex(scs_ptr->encode_context_ptr->sc_buffer_mutex);

    if (scs_ptr->encode_context_ptr->sc_frame_in == 0)
        eb_start_time(&context_ptr->first_in_pic_arrived_time_seconds,
                      &context_ptr->first_in_pic_arrived_timeu_seconds);
    else if (scs_ptr->encode_context_ptr->sc_frame_in == SC_FRAMES_TO_IGNORE)
        context_ptr->start_flag = EB_TRUE;
    // Compute duration since the start of the encode and since the previous checkpoint
    eb_finish_time(&curs_time_seconds, &curs_time_useconds);

    eb_compute_overall_elapsed_time_ms(context_ptr->first_in_pic_arrived_time_seconds,
                                       context_ptr->first_in_pic_arrived_timeu_seconds,
                                       curs_time_seconds,
                                       curs_time_useconds,
                                       &overall_duration);

    eb_compute_overall_elapsed_time_ms(context_ptr->prevs_time_seconds,
                                       context_ptr->prevs_timeu_seconds,
                                       curs_time_seconds,
                                       curs_time_useconds,
                                       &inst_duration);

    input_frames_count =
        (int64_t)overall_duration * (scs_ptr->static_config.injector_frame_rate >> 16) / 1000;
    scs_ptr->encode_context_ptr->sc_buffer =
        input_frames_count - scs_ptr->encode_context_ptr->sc_frame_in;

    encoder_mode_delta = 0;

    // Check every bufferTsshold1 for the changes (previous_frame_in_check1 variable)
    if ((scs_ptr->encode_context_ptr->sc_frame_in >
             context_ptr->previous_frame_in_check1 + buffer_threshold_1 &&
         scs_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)) {
        // Go to a slower mode based on the fullness and changes of the buffer
        if (scs_ptr->encode_context_ptr->sc_buffer < target_fps &&
            (context_ptr->prev_enc_mode_delta > -1 ||
             (context_ptr->prev_enc_mode_delta < 0 &&
              scs_ptr->encode_context_ptr->sc_frame_in >
                  context_ptr->previous_mode_change_frame_in + target_fps * 2))) {
            if (context_ptr->previous_buffer_check1 >
                scs_ptr->encode_context_ptr->sc_buffer + buffer_threshold_1) {
                encoder_mode_delta += -1;
                change_cond = 2;
            } else if (context_ptr->previous_mode_change_buffer >
                           buffer_threshold_1 + scs_ptr->encode_context_ptr->sc_buffer &&
                       scs_ptr->encode_context_ptr->sc_buffer < buffer_threshold_1) {
                encoder_mode_delta += -1;
                change_cond = 4;
            }
        }

        // Go to a faster mode based on the fullness and changes of the buffer
        if (scs_ptr->encode_context_ptr->sc_buffer >
            buffer_threshold_1 + context_ptr->previous_buffer_check1) {
            encoder_mode_delta += +1;
            change_cond = 1;
        } else if (scs_ptr->encode_context_ptr->sc_buffer >
                   buffer_threshold_1 + context_ptr->previous_mode_change_buffer) {
            encoder_mode_delta += +1;
            change_cond = 3;
        }

        // Update the encode mode based on the fullness of the buffer
        // If previous ChangeCond was the same, double the threshold2
        if (scs_ptr->encode_context_ptr->sc_buffer > buffer_threshold_3 &&
            (context_ptr->prev_change_cond != 7 ||
             scs_ptr->encode_context_ptr->sc_frame_in >
                 context_ptr->previous_mode_change_frame_in + buffer_threshold_2 * 2) &&
            scs_ptr->encode_context_ptr->sc_buffer > context_ptr->previous_mode_change_buffer) {
            encoder_mode_delta += 1;
            change_cond = 7;
        }
        encoder_mode_delta                    = CLIP3(-1, 1, encoder_mode_delta);
        scs_ptr->encode_context_ptr->enc_mode = (EbEncMode)CLIP3(
            1, MAX_ENC_PRESET, (int8_t)scs_ptr->encode_context_ptr->enc_mode + encoder_mode_delta);

        // Update previous stats
        context_ptr->previous_frame_in_check1 = scs_ptr->encode_context_ptr->sc_frame_in;
        context_ptr->previous_buffer_check1   = scs_ptr->encode_context_ptr->sc_buffer;

        if (encoder_mode_delta) {
            context_ptr->previous_mode_change_buffer   = scs_ptr->encode_context_ptr->sc_buffer;
            context_ptr->previous_mode_change_frame_in = scs_ptr->encode_context_ptr->sc_frame_in;
            context_ptr->prev_enc_mode_delta           = encoder_mode_delta;
        }
    }

    // Check every buffer_threshold_2 for the changes (previous_frame_in_check2 variable)
    if ((scs_ptr->encode_context_ptr->sc_frame_in >
             context_ptr->previous_frame_in_check2 + buffer_threshold_2 &&
         scs_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)) {
        encoder_mode_delta = 0;

        // if no change in the encoder mode and buffer is low enough and level is not increasing, switch to a slower encoder mode
        // If previous ChangeCond was the same, double the threshold2
        if (encoder_mode_delta == 0 &&
            scs_ptr->encode_context_ptr->sc_frame_in >
                context_ptr->previous_mode_change_frame_in + buffer_threshold_2 &&
            (context_ptr->prev_change_cond != 8 ||
             scs_ptr->encode_context_ptr->sc_frame_in >
                 context_ptr->previous_mode_change_frame_in + buffer_threshold_2 * 2) &&
            ((scs_ptr->encode_context_ptr->sc_buffer - context_ptr->previous_mode_change_buffer <
              (target_fps / 3)) ||
             context_ptr->previous_mode_change_buffer == 0) &&
            scs_ptr->encode_context_ptr->sc_buffer < buffer_threshold_3) {
            encoder_mode_delta = -1;
            change_cond        = 8;
        }

        encoder_mode_delta                    = CLIP3(-1, 1, encoder_mode_delta);
        scs_ptr->encode_context_ptr->enc_mode = (EbEncMode)CLIP3(
            1, MAX_ENC_PRESET, (int8_t)scs_ptr->encode_context_ptr->enc_mode + encoder_mode_delta);

        // Update previous stats
        context_ptr->previous_frame_in_check2 = scs_ptr->encode_context_ptr->sc_frame_in;

        if (encoder_mode_delta) {
            context_ptr->previous_mode_change_buffer   = scs_ptr->encode_context_ptr->sc_buffer;
            context_ptr->previous_mode_change_frame_in = scs_ptr->encode_context_ptr->sc_frame_in;
            context_ptr->prev_enc_mode_delta           = encoder_mode_delta;
        }
    }
    // Check every SC_FRAMES_INTERVAL_SPEED frames for the speed calculation (previous_frame_in_check3 variable)
    if (context_ptr->start_flag ||
        (scs_ptr->encode_context_ptr->sc_frame_in >
             context_ptr->previous_frame_in_check3 + SC_FRAMES_INTERVAL_SPEED &&
         scs_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)) {
        if (context_ptr->start_flag)
            context_ptr->cur_speed = (uint64_t)(scs_ptr->encode_context_ptr->sc_frame_out - 0) *
                                     1000 / (uint64_t)(overall_duration);
        else {
            if (inst_duration != 0)
                context_ptr->cur_speed = (uint64_t)(scs_ptr->encode_context_ptr->sc_frame_out -
                                                    context_ptr->prev_frame_out) *
                                         1000 / (uint64_t)(inst_duration);
        }
        context_ptr->start_flag = EB_FALSE;

        // Update previous stats
        context_ptr->previous_frame_in_check3 = scs_ptr->encode_context_ptr->sc_frame_in;
        context_ptr->prevs_time_seconds       = curs_time_seconds;
        context_ptr->prevs_timeu_seconds      = curs_time_useconds;
        context_ptr->prev_frame_out           = scs_ptr->encode_context_ptr->sc_frame_out;
    } else if (scs_ptr->encode_context_ptr->sc_frame_in < SC_FRAMES_TO_IGNORE &&
               (overall_duration != 0))
        context_ptr->cur_speed = (uint64_t)(scs_ptr->encode_context_ptr->sc_frame_out - 0) * 1000 /
                                 (uint64_t)(overall_duration);
    if (change_cond) context_ptr->prev_change_cond = change_cond;
    scs_ptr->encode_context_ptr->sc_frame_in++;
    if (scs_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)
        context_ptr->average_enc_mod += scs_ptr->encode_context_ptr->enc_mode;
    else
        context_ptr->average_enc_mod = 0;
    // Set the encoder level
    pcs_ptr->enc_mode = scs_ptr->encode_context_ptr->enc_mode;

    eb_release_mutex(scs_ptr->encode_context_ptr->sc_buffer_mutex);
    context_ptr->prev_enc_mod = scs_ptr->encode_context_ptr->enc_mode;
}

void reset_pcs_av1(PictureParentControlSet *pcs_ptr) {
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;
    Av1Common *  cm      = pcs_ptr->av1_cm;

    pcs_ptr->is_skip_mode_allowed = 0;
    pcs_ptr->skip_mode_flag       = 0;
    frm_hdr->frame_type           = KEY_FRAME;
    frm_hdr->show_frame           = 1;
    frm_hdr->showable_frame       = 1; // frame can be used as show existing frame in future
    // Flag for a frame used as a reference - not written to the Bitstream
    pcs_ptr->is_reference_frame = 0;
    // Flag signaling that the frame is encoded using only INTRA modes.
    pcs_ptr->intra_only = 0;
    // uint8_t last_intra_only;

    frm_hdr->disable_cdf_update      = 0;
    frm_hdr->allow_high_precision_mv = 0;
    frm_hdr->force_integer_mv        = 0; // 0 the default in AOM, 1 only integer
    frm_hdr->allow_warped_motion     = 0;

    /* profile settings */
#if CONFIG_ENTROPY_STATS
    int32_t coef_cdf_category;
#endif

    frm_hdr->quantization_params.base_q_idx              = 31;
    frm_hdr->quantization_params.delta_q_dc[AOM_PLANE_Y] = 0;
    frm_hdr->quantization_params.delta_q_ac[AOM_PLANE_Y] = 0;
    frm_hdr->quantization_params.delta_q_ac[AOM_PLANE_U] = 0;
    frm_hdr->quantization_params.delta_q_dc[AOM_PLANE_U] = 0;
    frm_hdr->quantization_params.delta_q_ac[AOM_PLANE_V] = 0;
    frm_hdr->quantization_params.delta_q_dc[AOM_PLANE_V] = 0;

    pcs_ptr->separate_uv_delta_q = 0;
    // Encoder
    frm_hdr->quantization_params.using_qmatrix   = 0;
    frm_hdr->quantization_params.qm[AOM_PLANE_Y] = 5;
    frm_hdr->quantization_params.qm[AOM_PLANE_U] = 5;
    frm_hdr->quantization_params.qm[AOM_PLANE_V] = 5;
    frm_hdr->is_motion_mode_switchable           = 0;
    // Flag signaling how frame contexts should be updated at the end of
    // a frame decode
    pcs_ptr->refresh_frame_context = REFRESH_FRAME_CONTEXT_DISABLED;

    frm_hdr->loop_filter_params.filter_level[0] = 0;
    frm_hdr->loop_filter_params.filter_level[1] = 0;
    frm_hdr->loop_filter_params.filter_level_u  = 0;
    frm_hdr->loop_filter_params.filter_level_v  = 0;
    frm_hdr->loop_filter_params.sharpness_level = 0;

    frm_hdr->loop_filter_params.mode_ref_delta_enabled = 0;
    frm_hdr->loop_filter_params.mode_ref_delta_update  = 0;
    frm_hdr->loop_filter_params.mode_deltas[0]         = 0;
    frm_hdr->loop_filter_params.mode_deltas[1]         = 0;

    frm_hdr->loop_filter_params.ref_deltas[0] = 1;
    frm_hdr->loop_filter_params.ref_deltas[1] = 0;
    frm_hdr->loop_filter_params.ref_deltas[2] = 0;
    frm_hdr->loop_filter_params.ref_deltas[3] = 0;
    frm_hdr->loop_filter_params.ref_deltas[4] = -1;
    frm_hdr->loop_filter_params.ref_deltas[5] = 0;
    frm_hdr->loop_filter_params.ref_deltas[6] = -1;
    frm_hdr->loop_filter_params.ref_deltas[7] = -1;

    frm_hdr->all_lossless      = 0;
    frm_hdr->coded_lossless    = 0;
    frm_hdr->reduced_tx_set    = 0;
    frm_hdr->reference_mode    = SINGLE_REFERENCE;
    pcs_ptr->frame_context_idx = 0; /* Context to use/update */
    for (int32_t i = 0; i < REF_FRAMES; i++) pcs_ptr->fb_of_context_type[i] = 0;
    frm_hdr->primary_ref_frame               = PRIMARY_REF_NONE;
    pcs_ptr->frame_offset                    = pcs_ptr->picture_number;
    frm_hdr->error_resilient_mode            = 0;
    cm->tiles_info.uniform_tile_spacing_flag = 1;
    pcs_ptr->large_scale_tile                = 0;
    pcs_ptr->film_grain_params_present       = 0;

    //cdef_pri_damping & cdef_sec_damping are consolidated to cdef_damping
    frm_hdr->cdef_params.cdef_damping = 0;
    //pcs_ptr->cdef_pri_damping = 0;
    //pcs_ptr->cdef_sec_damping = 0;

    pcs_ptr->nb_cdef_strengths = 1;
    for (int32_t i = 0; i < CDEF_MAX_STRENGTHS; i++) {
        frm_hdr->cdef_params.cdef_y_strength[i]  = 0;
        frm_hdr->cdef_params.cdef_uv_strength[i] = 0;
    }
    frm_hdr->cdef_params.cdef_bits            = 0;
    frm_hdr->delta_q_params.delta_q_present   = 1;
    frm_hdr->delta_lf_params.delta_lf_present = 0;
    frm_hdr->delta_q_params.delta_q_res       = DEFAULT_DELTA_Q_RES;
    frm_hdr->delta_lf_params.delta_lf_present = 0;
    frm_hdr->delta_lf_params.delta_lf_res     = 0;
    frm_hdr->delta_lf_params.delta_lf_multi   = 0;

    frm_hdr->current_frame_id           = 0;
    frm_hdr->frame_refs_short_signaling = 0;
    pcs_ptr->allow_comp_inter_inter     = 0;
    //  int32_t all_one_sided_refs;
}
/***********************************************
**** Copy the input buffer from the
**** sample application to the library buffers
************************************************/
static EbErrorType copy_frame_buffer(SequenceControlSet *scs_ptr, uint8_t *dst, uint8_t *src) {
    EbSvtAv1EncConfiguration *config       = &scs_ptr->static_config;
    EbErrorType               return_error = EB_ErrorNone;

    EbPictureBufferDesc *dst_picture_ptr = (EbPictureBufferDesc *)dst;
    EbPictureBufferDesc *src_picture_ptr = (EbPictureBufferDesc *)src;
    uint16_t             input_row_index;
    EbBool               is_16bit_input = (EbBool)(config->encoder_bit_depth > EB_8BIT);

    // Need to include for Interlacing on the fly with pictureScanType = 1

    if (!is_16bit_input) {
        uint32_t luma_buffer_offset =
            (dst_picture_ptr->stride_y * scs_ptr->top_padding + scs_ptr->left_padding)
            << is_16bit_input;
        uint32_t chroma_buffer_offset = (dst_picture_ptr->stride_cr * (scs_ptr->top_padding >> 1) +
                                         (scs_ptr->left_padding >> 1))
                                        << is_16bit_input;
        uint16_t luma_stride   = dst_picture_ptr->stride_y << is_16bit_input;
        uint16_t chroma_stride = dst_picture_ptr->stride_cb << is_16bit_input;
        uint16_t luma_width    = (uint16_t)(dst_picture_ptr->width - scs_ptr->max_input_pad_right)
                              << is_16bit_input;
        uint16_t chroma_width = (luma_width >> 1) << is_16bit_input;
        uint16_t luma_height  = (uint16_t)(dst_picture_ptr->height - scs_ptr->max_input_pad_bottom);

        //uint16_t     luma_height  = input_picture_ptr->max_height;
        // Y
        for (input_row_index = 0; input_row_index < luma_height; input_row_index++) {
            eb_memcpy(
                (dst_picture_ptr->buffer_y + luma_buffer_offset + luma_stride * input_row_index),
                (src_picture_ptr->buffer_y + luma_buffer_offset + luma_stride * input_row_index),
                luma_width);
        }

        // U
        for (input_row_index = 0; input_row_index<(luma_height>> 1); input_row_index++) {
            eb_memcpy((dst_picture_ptr->buffer_cb + chroma_buffer_offset +
                       chroma_stride * input_row_index),
                      (src_picture_ptr->buffer_cb + chroma_buffer_offset +
                       chroma_stride * input_row_index),
                      chroma_width);
        }

        // V
        for (input_row_index = 0; input_row_index<(luma_height>> 1); input_row_index++) {
            eb_memcpy((dst_picture_ptr->buffer_cr + chroma_buffer_offset +
                       chroma_stride * input_row_index),
                      (src_picture_ptr->buffer_cr + chroma_buffer_offset +
                       chroma_stride * input_row_index),
                      chroma_width);
        }
    } else if (is_16bit_input && config->compressed_ten_bit_format == 1) {
        {
            uint32_t luma_buffer_offset =
                (dst_picture_ptr->stride_y * scs_ptr->top_padding + scs_ptr->left_padding);
            uint32_t chroma_buffer_offset =
                (dst_picture_ptr->stride_cr * (scs_ptr->top_padding >> 1) +
                 (scs_ptr->left_padding >> 1));
            uint16_t luma_stride   = dst_picture_ptr->stride_y;
            uint16_t chroma_stride = dst_picture_ptr->stride_cb;
            uint16_t luma_width = (uint16_t)(dst_picture_ptr->width - scs_ptr->max_input_pad_right);
            uint16_t chroma_width = (luma_width >> 1);
            uint16_t luma_height =
                (uint16_t)(dst_picture_ptr->height - scs_ptr->max_input_pad_bottom);

            // Y 8bit
            for (input_row_index = 0; input_row_index < luma_height; input_row_index++) {
                eb_memcpy((dst_picture_ptr->buffer_y + luma_buffer_offset +
                           luma_stride * input_row_index),
                          (src_picture_ptr->buffer_y + luma_buffer_offset +
                           luma_stride * input_row_index),
                          luma_width);
            }

            // U 8bit
            for (input_row_index = 0; input_row_index<(luma_height>> 1); input_row_index++) {
                eb_memcpy((dst_picture_ptr->buffer_cb + chroma_buffer_offset +
                           chroma_stride * input_row_index),
                          (src_picture_ptr->buffer_cb + chroma_buffer_offset +
                           chroma_stride * input_row_index),
                          chroma_width);
            }

            // V 8bit
            for (input_row_index = 0; input_row_index<(luma_height>> 1); input_row_index++) {
                eb_memcpy((dst_picture_ptr->buffer_cr + chroma_buffer_offset +
                           chroma_stride * input_row_index),
                          (src_picture_ptr->buffer_cr + chroma_buffer_offset +
                           chroma_stride * input_row_index),
                          chroma_width);
            }
            // AMIR to update
            ////efficient copy - final
            ////compressed 2Bit in 1D format
            //{
            //    uint16_t luma_2bit_width = scs_ptr->max_input_luma_width / 4;
            //    uint16_t luma_height = scs_ptr->max_input_luma_height;

            //    uint16_t source_luma_2bit_stride = source_luma_stride / 4;
            //    uint16_t source_chroma_2bit_stride = source_luma_2bit_stride >> 1;

            //    for (input_row_index = 0; input_row_index < luma_height; input_row_index++) {
            //        eb_memcpy(input_picture_ptr->buffer_bit_inc_y + luma_2bit_width * input_row_index, input_ptr->luma_ext + source_luma_2bit_stride * input_row_index, luma_2bit_width);
            //    }
            //    for (input_row_index = 0; input_row_index < luma_height >> 1; input_row_index++) {
            //        eb_memcpy(input_picture_ptr->buffer_bit_inc_cb + (luma_2bit_width >> 1)*input_row_index, input_ptr->cb_ext + source_chroma_2bit_stride * input_row_index, luma_2bit_width >> 1);
            //    }
            //    for (input_row_index = 0; input_row_index < luma_height >> 1; input_row_index++) {
            //        eb_memcpy(input_picture_ptr->buffer_bit_inc_cr + (luma_2bit_width >> 1)*input_row_index, input_ptr->cr_ext + source_chroma_2bit_stride * input_row_index, luma_2bit_width >> 1);
            //    }
            //}
        }
    } else { // 10bit packed

        eb_memcpy(dst_picture_ptr->buffer_y, src_picture_ptr->buffer_y, src_picture_ptr->luma_size);

        eb_memcpy(
            dst_picture_ptr->buffer_cb, src_picture_ptr->buffer_cb, src_picture_ptr->chroma_size);

        eb_memcpy(
            dst_picture_ptr->buffer_cr, src_picture_ptr->buffer_cr, src_picture_ptr->chroma_size);

        eb_memcpy(dst_picture_ptr->buffer_bit_inc_y,
                  src_picture_ptr->buffer_bit_inc_y,
                  src_picture_ptr->luma_size);

        eb_memcpy(dst_picture_ptr->buffer_bit_inc_cb,
                  src_picture_ptr->buffer_bit_inc_cb,
                  src_picture_ptr->chroma_size);

        eb_memcpy(dst_picture_ptr->buffer_bit_inc_cr,
                  src_picture_ptr->buffer_bit_inc_cr,
                  src_picture_ptr->chroma_size);
    }
    return return_error;
}
static void copy_input_buffer(SequenceControlSet *sequenceControlSet, EbBufferHeaderType *dst,
                              EbBufferHeaderType *src) {
    // Copy the higher level structure
    dst->n_alloc_len  = src->n_alloc_len;
    dst->n_filled_len = src->n_filled_len;
    dst->flags        = src->flags;
    dst->pts          = src->pts;
    dst->n_tick_count = src->n_tick_count;
    dst->size         = src->size;
    dst->qp           = src->qp;
    dst->pic_type     = src->pic_type;

    // Copy the picture buffer
    if (src->p_buffer != NULL) copy_frame_buffer(sequenceControlSet, dst->p_buffer, src->p_buffer);
}
/******************************************************
 * Read Stat from File
 * reads StatStruct per frame from the file and stores under pcs_ptr
 ******************************************************/
static void read_stat_from_file(PictureParentControlSet *pcs_ptr, SequenceControlSet *scs_ptr) {
    eb_block_on_mutex(scs_ptr->encode_context_ptr->stat_file_mutex);

    int32_t fseek_return_value = fseek(scs_ptr->static_config.input_stat_file,
                                       (long)pcs_ptr->picture_number * sizeof(StatStruct),
                                       SEEK_SET);

    if (fseek_return_value != 0) {
        SVT_LOG("Error in fseek  returnVal %i\n", (int)fseek_return_value);
    }
    size_t fread_return_value = fread(&pcs_ptr->stat_struct,
                                      (size_t)1,
                                      sizeof(StatStruct),
                                      scs_ptr->static_config.input_stat_file);
    if (fread_return_value != sizeof(StatStruct)) {
        SVT_LOG("Error in freed  returnVal %i\n", (int)fread_return_value);
    }

    uint64_t referenced_area_avg          = 0;
    uint64_t referenced_area_has_non_zero = 0;
    for (int sb_addr = 0; sb_addr < scs_ptr->sb_total_count; ++sb_addr) {
        referenced_area_avg +=
            (pcs_ptr->stat_struct.referenced_area[sb_addr] /
             pcs_ptr->sb_params_array[sb_addr].width / pcs_ptr->sb_params_array[sb_addr].height);
        referenced_area_has_non_zero += pcs_ptr->stat_struct.referenced_area[sb_addr];
    }
    referenced_area_avg /= scs_ptr->sb_total_count;
    // adjust the reference area based on the intra refresh
    if (scs_ptr->intra_period_length && scs_ptr->intra_period_length < TWO_PASS_IR_THRSHLD)
        referenced_area_avg =
            referenced_area_avg * (scs_ptr->intra_period_length + 1) / TWO_PASS_IR_THRSHLD;
    pcs_ptr->referenced_area_avg          = referenced_area_avg;
    pcs_ptr->referenced_area_has_non_zero = referenced_area_has_non_zero ? 1 : 0;
    eb_release_mutex(scs_ptr->encode_context_ptr->stat_file_mutex);
}


/* Resource Coordination Kernel */
/*********************************************************************************
*
* @brief
*  The Resource Coordination Process is the first stage that input pictures
*  this process is a single threaded, picture-based process that handles one picture at a time
*  in display order
*
* @par Description:
*  Input input picture samples are available once the input_buffer_fifo_ptr queue gets any items
*  The Resource Coordination Process assembles the input information and creates
*  the appropriate buffers that would travel with the input picture all along
*  the encoding pipeline and passes this data along with the current encoder settings
*  to the picture analysis process
*  Encoder settings include, but are not limited to QPs, picture type, encoding
*  parameters that change per picture sequence
*
* @param[in] EbBufferHeaderType
*  EbBufferHeaderType containing the input picture samples along with settings specific to that picture
*
* @param[out] Input picture in Picture buffers
*  Initialized picture level (PictureParentControlSet) / sequence level
*  (SequenceControlSet if it's the initial picture) structures
*
* @param[out] Settings
*  Encoder settings include picture timing and order settings (POC) resolution settings, sequence level
*  parameters (if it is the initial picture) and other encoding parameters such as QP, Bitrate, picture type ...
*
********************************************************************************/
void *resource_coordination_kernel(void *input_ptr) {
    EbThreadContext *            enc_contxt_ptr = (EbThreadContext *)input_ptr;
    ResourceCoordinationContext *context_ptr = (ResourceCoordinationContext *)enc_contxt_ptr->priv;

    EbObjectWrapper *pcs_wrapper_ptr;

    PictureParentControlSet *pcs_ptr;

    EbObjectWrapper *   prev_scs_wrapper_ptr;
    SequenceControlSet *scs_ptr;

    EbObjectWrapper *            eb_input_wrapper_ptr;
    EbBufferHeaderType *         eb_input_ptr;
    EbObjectWrapper *            output_wrapper_ptr;
    ResourceCoordinationResults *out_results_ptr;

    EbObjectWrapper *input_picture_wrapper_ptr;
    EbObjectWrapper *reference_picture_wrapper_ptr;

    uint32_t instance_index;
    EbBool   end_of_sequence_flag = EB_FALSE;

    uint32_t         input_size           = 0;
    EbObjectWrapper *prev_pcs_wrapper_ptr = 0;

    for (;;) {
        // Tie instance_index to zero for now...
        instance_index = 0;

        // Get the Next svt Input Buffer [BLOCKING]
        EB_GET_FULL_OBJECT(context_ptr->input_buffer_fifo_ptr, &eb_input_wrapper_ptr);

        eb_input_ptr = (EbBufferHeaderType *)eb_input_wrapper_ptr->object_ptr;
        scs_ptr      = context_ptr->scs_instance_array[instance_index]->scs_ptr;

        // If config changes occured since the last picture began encoding, then
        //   prepare a new scs_ptr containing the new changes and update the state
        //   of the previous Active SequenceControlSet
        eb_block_on_mutex(context_ptr->scs_instance_array[instance_index]->config_mutex);
        if (context_ptr->scs_instance_array[instance_index]->encode_context_ptr->initial_picture) {
            // Update picture width, picture height, cropping right offset, cropping bottom offset, and conformance windows
            if (context_ptr->scs_instance_array[instance_index]
                    ->encode_context_ptr->initial_picture)

            {
                context_ptr->scs_instance_array[instance_index]
                    ->scs_ptr->seq_header.max_frame_width =
                    context_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_width;
                context_ptr->scs_instance_array[instance_index]
                    ->scs_ptr->seq_header.max_frame_height =
                    context_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_luma_height;
                context_ptr->scs_instance_array[instance_index]->scs_ptr->chroma_width =
                    (context_ptr->scs_instance_array[instance_index]
                         ->scs_ptr->max_input_luma_width >>
                     1);
                context_ptr->scs_instance_array[instance_index]->scs_ptr->chroma_height =
                    (context_ptr->scs_instance_array[instance_index]
                         ->scs_ptr->max_input_luma_height >>
                     1);

                context_ptr->scs_instance_array[instance_index]->scs_ptr->pad_right =
                    context_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_pad_right;
                context_ptr->scs_instance_array[instance_index]->scs_ptr->pad_bottom =
                    context_ptr->scs_instance_array[instance_index]->scs_ptr->max_input_pad_bottom;
                input_size = context_ptr->scs_instance_array[instance_index]
                                 ->scs_ptr->seq_header.max_frame_width *
                             context_ptr->scs_instance_array[instance_index]
                                 ->scs_ptr->seq_header.max_frame_height;
            }

            // Copy previous Active SequenceControlSetPtr to a place holder
            prev_scs_wrapper_ptr = context_ptr->sequence_control_set_active_array[instance_index];

            // Get empty SequenceControlSet [BLOCKING]
            eb_get_empty_object(context_ptr->sequence_control_set_empty_fifo_ptr,
                                &context_ptr->sequence_control_set_active_array[instance_index]);

            // Copy the contents of the active SequenceControlSet into the new empty SequenceControlSet
            copy_sequence_control_set(
                (SequenceControlSet *)context_ptr->sequence_control_set_active_array[instance_index]
                    ->object_ptr,
                context_ptr->scs_instance_array[instance_index]->scs_ptr);

            // Disable releaseFlag of new SequenceControlSet
            eb_object_release_disable(
                context_ptr->sequence_control_set_active_array[instance_index]);

            if (prev_scs_wrapper_ptr != NULL) {
                // Enable releaseFlag of old SequenceControlSet
                eb_object_release_enable(prev_scs_wrapper_ptr);

                // Check to see if previous SequenceControlSet is already inactive, if TRUE then release the SequenceControlSet
                if (prev_scs_wrapper_ptr->live_count == 0) {
                    eb_release_object(prev_scs_wrapper_ptr);
                }
            }
        }
        eb_release_mutex(context_ptr->scs_instance_array[instance_index]->config_mutex);
        // Seque Control Set is released by Rate Control after passing through MDC->MD->ENCDEC->Packetization->RateControl,
        // in the PictureManager after receiving the reference and in PictureManager after receiving the feedback
        eb_object_inc_live_count(context_ptr->sequence_control_set_active_array[instance_index], 3);

        // Set the current SequenceControlSet
        scs_ptr =
            (SequenceControlSet *)context_ptr->sequence_control_set_active_array[instance_index]
                ->object_ptr;

        // Init SB Params
        if (context_ptr->scs_instance_array[instance_index]->encode_context_ptr->initial_picture) {
            derive_input_resolution(&scs_ptr->input_resolution, input_size);

            sb_params_init(scs_ptr);
            sb_geom_init(scs_ptr);
            scs_ptr->enable_altrefs = scs_ptr->static_config.enable_altrefs ? EB_TRUE : EB_FALSE;

            // initialize sequence level enable_superres
            scs_ptr->seq_header.enable_superres = 0;

            if (scs_ptr->static_config.inter_intra_compound == DEFAULT) {
                // Set inter-intra mode      Settings
                // 0                 OFF
                // 1                 ON
                scs_ptr->seq_header.enable_interintra_compound =
#if MAY12_ADOPTIONS
#if UNIFY_SC_NSC
                (scs_ptr->static_config.enc_mode <= ENC_M2) ? 1 : 0;
#else
#if PRESET_SHIFITNG
                (scs_ptr->static_config.enc_mode <= ENC_M2 &&
#else
                (scs_ptr->static_config.enc_mode <= ENC_M4 &&
#endif
                    scs_ptr->static_config.screen_content_mode != 1)
                    ? 1
                    : 0;
#endif
#else
#if APR24_M3_ADOPTIONS
                (scs_ptr->static_config.enc_mode <= ENC_M2 &&
                    scs_ptr->static_config.screen_content_mode != 1)
                    ? 1
                    : 0;
#else
#if APR23_ADOPTIONS
                (scs_ptr->static_config.enc_mode <= ENC_M5 &&
                    scs_ptr->static_config.screen_content_mode != 1)
                    ? 1
                    : 0;
#else
#if PRESETS_SHIFT
                    (scs_ptr->static_config.enc_mode <= ENC_M2 &&
                    scs_ptr->static_config.screen_content_mode != 1)
                    ? 1
                    : 0;
#else
#if MAR18_MR_TESTS_ADOPTIONS
                    (scs_ptr->static_config.enc_mode <= ENC_M3 &&
                    scs_ptr->static_config.screen_content_mode != 1)
                    ? 1
                    : 0;
#else
#if MAR3_M2_ADOPTIONS
#if MAR4_M3_ADOPTIONS
                    MR_MODE || (scs_ptr->static_config.enc_mode <= ENC_M3 &&
#else
                    MR_MODE || (scs_ptr->static_config.enc_mode <= ENC_M2 &&
#endif
#else
                    MR_MODE || (scs_ptr->static_config.enc_mode <= ENC_M1 &&
#endif
                                scs_ptr->static_config.screen_content_mode != 1)
                        ? 1
                        : 0;
#endif
#endif
#endif
#endif
#endif

            } else
                scs_ptr->seq_header.enable_interintra_compound =
                    scs_ptr->static_config.inter_intra_compound;
#if FILTER_INTRA_CLI
            // Enable/Disable Filter Intra
            // seq_header.filter_intra_level | Settings
            // 0                             | Disable
            // 1                             | Enable
            if (scs_ptr->static_config.filter_intra_level == DEFAULT)
#if MAY19_ADOPTIONS
                scs_ptr->seq_header.filter_intra_level =
#if JUNE17_ADOPTIONS
                (scs_ptr->static_config.enc_mode <= ENC_M6) ? 1 : 0;
#else
#if PRESET_SHIFITNG
                (scs_ptr->static_config.enc_mode <= ENC_M4) ? 1 : 0;
#else
                (scs_ptr->static_config.enc_mode <= ENC_M6) ? 1 : 0;
#endif
#endif
#else
#if APR23_ADOPTIONS_2
                scs_ptr->seq_header.enable_filter_intra =
                (scs_ptr->static_config.enc_mode <= ENC_M5) ? 1 : 0;
#else
#if PRESETS_SHIFT
                scs_ptr->seq_header.enable_filter_intra =
                (scs_ptr->static_config.enc_mode <= ENC_M4) ? 1 : 0;
#else
#if MAR10_ADOPTIONS
                if (scs_ptr->static_config.screen_content_mode == 1)
#if MAR17_ADOPTIONS
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M7) ? 1 : 0;
#else
#if MAR12_ADOPTIONS
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M3) ? 1 : 0;
#else
#if MAR11_ADOPTIONS
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M1) ? 1 : 0;
#else
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M2) ? 1 : 0;
#endif
#endif
#endif
                else
#endif
#if MAR17_ADOPTIONS
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M7) ? 1 : 0;
#else
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M4) ? 1 : 0;
#endif
#endif
#endif
#endif
            else
                scs_ptr->seq_header.filter_intra_level = (scs_ptr->static_config.filter_intra_level == 0) ? 0 : 1;
#else
            // Set filter intra mode      Settings
            // 0                 OFF
            // 1                 ON
            if (scs_ptr->static_config.enable_filter_intra)
#if MAY19_ADOPTIONS
                scs_ptr->seq_header.enable_filter_intra =
#if JUNE17_ADOPTIONS
                (scs_ptr->static_config.enc_mode <= ENC_M6) ? 1 : 0;
#else
#if PRESET_SHIFITNG
                (scs_ptr->static_config.enc_mode <= ENC_M4) ? 1 : 0;
#else
                (scs_ptr->static_config.enc_mode <= ENC_M6) ? 1 : 0;
#endif
#endif
#else
#if APR23_ADOPTIONS_2
                scs_ptr->seq_header.enable_filter_intra =
                (scs_ptr->static_config.enc_mode <= ENC_M5) ? 1 : 0;
#else
#if PRESETS_SHIFT
                scs_ptr->seq_header.enable_filter_intra =
                (scs_ptr->static_config.enc_mode <= ENC_M4) ? 1 : 0;
#else
#if MAR10_ADOPTIONS
                if (scs_ptr->static_config.screen_content_mode == 1)
#if MAR17_ADOPTIONS
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M7) ? 1 : 0;
#else
#if MAR12_ADOPTIONS
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M3) ? 1 : 0;
#else
#if MAR11_ADOPTIONS
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M1) ? 1 : 0;
#else
                    scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M2) ? 1 : 0;
#endif
#endif
#endif
                else
#endif
#if MAR17_ADOPTIONS
                scs_ptr->seq_header.enable_filter_intra =
                (scs_ptr->static_config.enc_mode <= ENC_M7) ? 1 : 0;
#else
                scs_ptr->seq_header.enable_filter_intra =
                    (scs_ptr->static_config.enc_mode <= ENC_M4) ? 1 : 0;
#endif
#endif
#endif
#endif
            else
                scs_ptr->seq_header.enable_filter_intra = 0;
#endif
            // Set compound mode      Settings
            // 0                 OFF: No compond mode search : AVG only
            // 1                 ON: full
            if (scs_ptr->static_config.compound_level == DEFAULT) {
#if MAR11_ADOPTIONS
                scs_ptr->compound_mode = (scs_ptr->static_config.enc_mode <= ENC_M8) ? 1 : 0;
#else
                scs_ptr->compound_mode = (scs_ptr->static_config.enc_mode <= ENC_M4) ? 1 : 0;
#endif
            } else
                scs_ptr->compound_mode = scs_ptr->static_config.compound_level;

#if M8_NEW_REF
                scs_ptr->compound_mode = 1 ;
#endif
            if (scs_ptr->compound_mode) {
                scs_ptr->seq_header.order_hint_info.enable_jnt_comp = 1; //DISTANCE
                scs_ptr->seq_header.enable_masked_compound          = 1; //DIFF+WEDGE
            } else {
                scs_ptr->seq_header.order_hint_info.enable_jnt_comp = 0;
                scs_ptr->seq_header.enable_masked_compound          = 0;
            }
        }
        // Since at this stage we do not know the prediction structure and the location of ALT_REF pictures,
        // for every picture (except first picture), we allocate two: 1. original picture, 2. potential Overlay picture.
        // In Picture Decision Process, where the overlay frames are known, they extra pictures are released
        uint8_t has_overlay =
            (scs_ptr->static_config.enable_overlays == EB_FALSE ||
             context_ptr->scs_instance_array[instance_index]->encode_context_ptr->initial_picture)
                ? 0
                : 1;
        for (uint8_t loop_index = 0; loop_index <= has_overlay && !end_of_sequence_flag;
             loop_index++) {
            //Get a New ParentPCS where we will hold the new input_picture
            eb_get_empty_object(context_ptr->picture_control_set_fifo_ptr_array[instance_index],
                                &pcs_wrapper_ptr);

            // Parent PCS is released by the Rate Control after passing through MDC->MD->ENCDEC->Packetization
            eb_object_inc_live_count(pcs_wrapper_ptr, 1);

            pcs_ptr = (PictureParentControlSet *)pcs_wrapper_ptr->object_ptr;

            pcs_ptr->p_pcs_wrapper_ptr = pcs_wrapper_ptr;

            pcs_ptr->sb_params_array = scs_ptr->sb_params_array;
            pcs_ptr->sb_geom = scs_ptr->sb_geom;
            pcs_ptr->input_resolution = scs_ptr->input_resolution;
            pcs_ptr->picture_sb_width = scs_ptr->pic_width_in_sb;
            pcs_ptr->picture_sb_height = scs_ptr->picture_height_in_sb;

            pcs_ptr->overlay_ppcs_ptr = NULL;
            pcs_ptr->is_alt_ref       = 0;
            if (loop_index) {
                pcs_ptr->is_overlay = 1;
                // set the overlay_ppcs_ptr in the original (ALT_REF) ppcs to the current ppcs
                EbObjectWrapper *alt_ref_picture_control_set_wrapper_ptr =
                    (context_ptr->scs_instance_array[instance_index]
                         ->encode_context_ptr->initial_picture)
                        ? pcs_wrapper_ptr
                        : scs_ptr->encode_context_ptr->previous_picture_control_set_wrapper_ptr;

                pcs_ptr->alt_ref_ppcs_ptr =
                    ((PictureParentControlSet *)
                         alt_ref_picture_control_set_wrapper_ptr->object_ptr);
                pcs_ptr->alt_ref_ppcs_ptr->overlay_ppcs_ptr = pcs_ptr;
            } else {
                pcs_ptr->is_overlay       = 0;
                pcs_ptr->alt_ref_ppcs_ptr = NULL;
            }
            // Set the Encoder mode
            pcs_ptr->enc_mode = scs_ptr->static_config.enc_mode;

            // Keep track of the previous input for the ZZ SADs computation
            pcs_ptr->previous_picture_control_set_wrapper_ptr =
                (context_ptr->scs_instance_array[instance_index]
                     ->encode_context_ptr->initial_picture)
                    ? pcs_wrapper_ptr
                    : scs_ptr->encode_context_ptr->previous_picture_control_set_wrapper_ptr;
            if (loop_index == 0)
                scs_ptr->encode_context_ptr->previous_picture_control_set_wrapper_ptr =
                    pcs_wrapper_ptr;
            // Copy data from the svt buffer to the input frame
            // *Note - Assumes 4:2:0 planar
            input_picture_wrapper_ptr     = eb_input_wrapper_ptr;
            pcs_ptr->enhanced_picture_ptr = (EbPictureBufferDesc *)eb_input_ptr->p_buffer;
            pcs_ptr->input_ptr            = eb_input_ptr;
            end_of_sequence_flag =
                (pcs_ptr->input_ptr->flags & EB_BUFFERFLAG_EOS) ? EB_TRUE : EB_FALSE;
            eb_start_time(&pcs_ptr->start_time_seconds, &pcs_ptr->start_time_u_seconds);

            pcs_ptr->scs_wrapper_ptr =
                context_ptr->sequence_control_set_active_array[instance_index];
            pcs_ptr->scs_ptr                   = scs_ptr;
            pcs_ptr->input_picture_wrapper_ptr = input_picture_wrapper_ptr;
            pcs_ptr->end_of_sequence_flag      = end_of_sequence_flag;

            if (loop_index == 1) {
                // Get a new input picture for overlay.
                EbObjectWrapper *input_pic_wrapper_ptr;

                // Get a new input picture for overlay.
                eb_get_empty_object(
                    scs_ptr->encode_context_ptr->overlay_input_picture_pool_fifo_ptr,
                    &input_pic_wrapper_ptr);

                // Copy from original picture (pcs_ptr->input_picture_wrapper_ptr), which is shared between overlay and alt_ref up to this point, to the new input picture.
                if (pcs_ptr->alt_ref_ppcs_ptr->input_picture_wrapper_ptr->object_ptr != NULL) {
                    copy_input_buffer(scs_ptr,
                                      (EbBufferHeaderType *)input_pic_wrapper_ptr->object_ptr,
                                      (EbBufferHeaderType *)pcs_ptr->alt_ref_ppcs_ptr
                                          ->input_picture_wrapper_ptr->object_ptr);
                }
                // Assign the new picture to the new pointers
                pcs_ptr->input_ptr = (EbBufferHeaderType *)input_pic_wrapper_ptr->object_ptr;
                pcs_ptr->enhanced_picture_ptr = (EbPictureBufferDesc *)pcs_ptr->input_ptr->p_buffer;
                pcs_ptr->input_picture_wrapper_ptr = input_pic_wrapper_ptr;
            }
            // Set Picture Control Flags
            pcs_ptr->idr_flag = scs_ptr->encode_context_ptr->initial_picture ||
                                (pcs_ptr->input_ptr->pic_type == EB_AV1_KEY_PICTURE);
            pcs_ptr->cra_flag =
                (pcs_ptr->input_ptr->pic_type == EB_AV1_INTRA_ONLY_PICTURE) ? EB_TRUE : EB_FALSE;
            pcs_ptr->scene_change_flag = EB_FALSE;
            pcs_ptr->qp_on_the_fly     = EB_FALSE;
            pcs_ptr->sb_total_count    = scs_ptr->sb_total_count;
            if (scs_ptr->static_config.speed_control_flag) {
                speed_buffer_control(context_ptr, pcs_ptr, scs_ptr);
            } else
                pcs_ptr->enc_mode = (EbEncMode)scs_ptr->static_config.enc_mode;
            //  If the mode of the second pass is not set from CLI, it is set to enc_mode
            pcs_ptr->snd_pass_enc_mode =
                (scs_ptr->use_output_stat_file &&
                 scs_ptr->static_config.snd_pass_enc_mode != MAX_ENC_PRESET + 1)
                    ? (EbEncMode)scs_ptr->static_config.snd_pass_enc_mode
                    : pcs_ptr->enc_mode;

            // Set the SCD Mode
            scs_ptr->scd_mode =
                scs_ptr->static_config.scene_change_detection == 0 ? SCD_MODE_0 : SCD_MODE_1;

            // Set the block mean calculation prec
            scs_ptr->block_mean_calc_prec = BLOCK_MEAN_PREC_SUB;

            // Pre-Analysis Signal(s) derivation
            signal_derivation_pre_analysis_oq(scs_ptr, pcs_ptr);
            pcs_ptr->filtered_sse    = 0;
            pcs_ptr->filtered_sse_uv = 0;
            // Rate Control
            // Set the ME Distortion and OIS Historgrams to zero
            if (scs_ptr->static_config.rate_control_mode) {
                EB_MEMSET(pcs_ptr->me_distortion_histogram,
                          0,
                          NUMBER_OF_SAD_INTERVALS * sizeof(uint16_t));
                EB_MEMSET(pcs_ptr->ois_distortion_histogram,
                          0,
                          NUMBER_OF_INTRA_SAD_INTERVALS * sizeof(uint16_t));
            }
            pcs_ptr->full_sb_count = 0;

            if (scs_ptr->static_config.use_qp_file == 1) {
                pcs_ptr->qp_on_the_fly = EB_TRUE;
                if (pcs_ptr->input_ptr->qp > MAX_QP_VALUE) {
                    SVT_LOG("SVT [WARNING]: INPUT QP OUTSIDE OF RANGE\n");
                    pcs_ptr->qp_on_the_fly = EB_FALSE;
                    pcs_ptr->picture_qp    = (uint8_t)scs_ptr->static_config.qp;
                }
                pcs_ptr->picture_qp = (uint8_t)pcs_ptr->input_ptr->qp;
            } else {
                pcs_ptr->qp_on_the_fly = EB_FALSE;
                pcs_ptr->picture_qp    = (uint8_t)scs_ptr->static_config.qp;
            }

            // Picture Stats
            if (loop_index == has_overlay || end_of_sequence_flag)
                pcs_ptr->picture_number = context_ptr->picture_number_array[instance_index]++;
            else
                pcs_ptr->picture_number = context_ptr->picture_number_array[instance_index];
            reset_pcs_av1(pcs_ptr);
            if (scs_ptr->use_input_stat_file && !end_of_sequence_flag)
                read_stat_from_file(pcs_ptr, scs_ptr);
            else {
                memset(&pcs_ptr->stat_struct, 0, sizeof(StatStruct));
            }
            scs_ptr->encode_context_ptr->initial_picture = EB_FALSE;

            // Get Empty Reference Picture Object
            eb_get_empty_object(scs_ptr->encode_context_ptr->pa_reference_picture_pool_fifo_ptr,
                                &reference_picture_wrapper_ptr);

            pcs_ptr->pa_reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
            // Since overlay pictures are not added to PA_Reference queue in PD and not released there, the life count is only set to 1
            if (pcs_ptr->is_overlay)
                // Give the new Reference a nominal live_count of 1
                eb_object_inc_live_count(pcs_ptr->pa_reference_picture_wrapper_ptr, 1);
            else
                eb_object_inc_live_count(pcs_ptr->pa_reference_picture_wrapper_ptr, 2);
            if (scs_ptr->static_config.unrestricted_motion_vector == 0) {
                struct PictureParentControlSet *ppcs_ptr = pcs_ptr;
                Av1Common *const                cm       = ppcs_ptr->av1_cm;
                uint8_t                         pic_width_in_sb =
                    (uint8_t)((pcs_ptr->aligned_width + scs_ptr->sb_size_pix - 1) /
                              scs_ptr->sb_size_pix);
                int       tile_row, tile_col;
                uint32_t  x_sb_index, y_sb_index;
                const int tile_cols = cm->tiles_info.tile_cols;
                const int tile_rows = cm->tiles_info.tile_rows;
                TileInfo  tile_info;
                int       sb_size_log2 = scs_ptr->seq_header.sb_size_log2;
                //Tile Loop
                for (tile_row = 0; tile_row < tile_rows; tile_row++) {
                    eb_av1_tile_set_row(&tile_info, &cm->tiles_info, cm->mi_rows, tile_row);

                    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
                        eb_av1_tile_set_col(&tile_info, &cm->tiles_info, cm->mi_cols, tile_col);

                        for ((y_sb_index =
                                  cm->tiles_info.tile_row_start_mi[tile_row] >> sb_size_log2);
                             (y_sb_index<((uint32_t)cm->tiles_info.tile_row_start_mi[tile_row + 1]>>
                              sb_size_log2));
                             y_sb_index++) {
                            for ((x_sb_index =
                                      cm->tiles_info.tile_col_start_mi[tile_col] >> sb_size_log2);
                                 (x_sb_index<((uint32_t)cm->tiles_info
                                                 .tile_col_start_mi[tile_col + 1]>> sb_size_log2));
                                 x_sb_index++) {
                                int sb_index =
                                    (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                                scs_ptr->sb_params_array[sb_index].tile_start_x =
                                    4 * tile_info.mi_col_start;
                                scs_ptr->sb_params_array[sb_index].tile_end_x =
                                    4 * tile_info.mi_col_end;
                                scs_ptr->sb_params_array[sb_index].tile_start_y =
                                    4 * tile_info.mi_row_start;
                                scs_ptr->sb_params_array[sb_index].tile_end_y =
                                    4 * tile_info.mi_row_end;
                            }
                        }
                    }
                }
            }

            // Get Empty Output Results Object
            if (pcs_ptr->picture_number > 0 && (prev_pcs_wrapper_ptr != NULL)) {
                ((PictureParentControlSet *)prev_pcs_wrapper_ptr->object_ptr)
                    ->end_of_sequence_flag = end_of_sequence_flag;
                eb_get_empty_object(context_ptr->resource_coordination_results_output_fifo_ptr,
                                    &output_wrapper_ptr);
                out_results_ptr = (ResourceCoordinationResults *)output_wrapper_ptr->object_ptr;
                out_results_ptr->pcs_wrapper_ptr = prev_pcs_wrapper_ptr;
                // since overlay frame has the end of sequence set properly, set the end of sequence to true in the alt ref picture
                if (((PictureParentControlSet *)prev_pcs_wrapper_ptr->object_ptr)->is_overlay &&
                    end_of_sequence_flag)
                    ((PictureParentControlSet *)prev_pcs_wrapper_ptr->object_ptr)
                        ->alt_ref_ppcs_ptr->end_of_sequence_flag = EB_TRUE;
                // Post the finished Results Object
                eb_post_full_object(output_wrapper_ptr);
            }
            prev_pcs_wrapper_ptr = pcs_wrapper_ptr;
        }
    }

    return NULL;
}
