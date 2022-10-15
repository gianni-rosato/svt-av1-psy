/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
#include "EbResize.h"
#include "EbTime.h"
#include "EbObject.h"
#include "EbLog.h"
#include "pass2_strategy.h"
#include "common_dsp_rtcd.h"
#include "EbResize.h"
#include "EbMetadataHandle.h"
typedef struct ResourceCoordinationContext {
    EbFifo                        *input_cmd_fifo_ptr;
    EbFifo                        *resource_coordination_results_output_fifo_ptr;
    EbFifo                       **picture_control_set_fifo_ptr_array;
    EbSequenceControlSetInstance **scs_instance_array;
    EbCallback                   **app_callback_ptr_array;

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
    Bool     start_flag;
} ResourceCoordinationContext;

static void resource_coordination_context_dctor(EbPtr p) {
    EbThreadContext *thread_contxt_ptr = (EbThreadContext *)p;
    if (thread_contxt_ptr->priv) {
        ResourceCoordinationContext *obj = (ResourceCoordinationContext *)thread_contxt_ptr->priv;
        EB_FREE_ARRAY(obj->picture_number_array);
        EB_FREE_ARRAY(obj->picture_control_set_fifo_ptr_array);
        EB_FREE_ARRAY(obj);
    }
}

/************************************************
 * Resource Coordination Context Constructor
 ************************************************/
EbErrorType resource_coordination_context_ctor(EbThreadContext *thread_contxt_ptr,
                                               EbEncHandle     *enc_handle_ptr) {
    ResourceCoordinationContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_contxt_ptr->priv  = context_ptr;
    thread_contxt_ptr->dctor = resource_coordination_context_dctor;

    EB_MALLOC_ARRAY(context_ptr->picture_control_set_fifo_ptr_array,
                    enc_handle_ptr->encode_instance_total_count);
    for (uint32_t i = 0; i < enc_handle_ptr->encode_instance_total_count; i++) {
        //ResourceCoordination works with ParentPCS
        context_ptr->picture_control_set_fifo_ptr_array[i] = svt_system_resource_get_producer_fifo(
            enc_handle_ptr->picture_parent_control_set_pool_ptr_array[i], 0);
    }
    context_ptr->input_cmd_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->input_cmd_resource_ptr, 0);
    context_ptr->resource_coordination_results_output_fifo_ptr =
        svt_system_resource_get_producer_fifo(
            enc_handle_ptr->resource_coordination_results_resource_ptr, 0);
    context_ptr->scs_instance_array     = enc_handle_ptr->scs_instance_array;
    context_ptr->app_callback_ptr_array = enc_handle_ptr->app_callback_ptr_array;
    context_ptr->compute_segments_total_count_array =
        enc_handle_ptr->compute_segments_total_count_array;
    context_ptr->encode_instances_total_count = enc_handle_ptr->encode_instance_total_count;

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
    context_ptr->start_flag                         = FALSE;

    context_ptr->previous_buffer_check1 = 0;
    context_ptr->prev_change_cond       = 0;

    return EB_ErrorNone;
}

uint8_t get_wn_filter_level(EncMode enc_mode, uint8_t input_resolution, Bool is_ref);
uint8_t get_sg_filter_level(EncMode enc_mode, Bool fast_decode, uint8_t input_resolution,
                            Bool is_base);
/*
* return true if restoration filtering is enabled; false otherwise
  Used by signal_derivation_pre_analysis_oq and memory allocation
*/
uint8_t get_enable_restoration(EncMode enc_mode, int8_t config_enable_restoration,
                               uint8_t input_resolution, Bool fast_decode) {
    if (config_enable_restoration != DEFAULT)
        return config_enable_restoration;

    uint8_t wn = 0;
    for (int is_ref = 0; is_ref < 2; is_ref++) {
        wn = get_wn_filter_level(enc_mode, input_resolution, is_ref);
        if (wn)
            break;
    }

    uint8_t sg = 0;
    for (int is_base = 0; is_base < 2; is_base++) {
        sg = get_sg_filter_level(enc_mode, fast_decode, input_resolution, is_base);
        if (sg)
            break;
    }

    return (sg > 0 || wn > 0);
}

/******************************************************
* Derive Pre-Analysis settings for OQ for pcs
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
EbErrorType signal_derivation_pre_analysis_oq_pcs(PictureParentControlSet *pcs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    // Derive HME Flag
    // Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
    // Will be later updated for SC/NSC in PictureDecisionProcess
    pcs_ptr->enable_hme_flag        = 1;
    pcs_ptr->enable_hme_level0_flag = 1;
    pcs_ptr->enable_hme_level1_flag = 1;
    pcs_ptr->enable_hme_level2_flag = 1;
    // Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
    // Will be later updated for SC/NSC in PictureDecisionProcess
    pcs_ptr->tf_enable_hme_flag        = 1;
    pcs_ptr->tf_enable_hme_level0_flag = 1;
    pcs_ptr->tf_enable_hme_level1_flag = 1;
    pcs_ptr->tf_enable_hme_level2_flag = 1;

    //if (scs_ptr->static_config.enable_tpl_la)
    //assert_err(pcs_ptr->is_720p_or_larger == (pcs_ptr->tpl_ctrls.synth_blk_size == 16), "TPL Synth Size Error");

    return return_error;
}

/******************************************************
* Derive Pre-Analysis settings for OQ for scs
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
EbErrorType signal_derivation_pre_analysis_oq_scs(SequenceControlSet *scs_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    // Set the SCD Mode
    scs_ptr->scd_mode = scs_ptr->static_config.scene_change_detection == 0 ? SCD_MODE_0
                                                                           : SCD_MODE_1;

    // initialize sequence level enable_superres
    scs_ptr->seq_header.enable_superres = scs_ptr->static_config.superres_mode > SUPERRES_NONE ? 1
                                                                                               : 0;
    if (scs_ptr->inter_intra_compound == DEFAULT)
        scs_ptr->seq_header.enable_interintra_compound = 1;
    else
        scs_ptr->seq_header.enable_interintra_compound = scs_ptr->inter_intra_compound;
    // Enable/Disable Filter Intra
    // seq_header.filter_intra_level | Settings
    // 0                             | Disable
    // 1                             | Enable
    if (scs_ptr->filter_intra_level == DEFAULT)
        scs_ptr->seq_header.filter_intra_level = 1;
    else
        scs_ptr->seq_header.filter_intra_level = scs_ptr->filter_intra_level;
    // Set compound mode      Settings
    // 0                 OFF: No compond mode search : AVG only
    // 1                 ON: full
    if (scs_ptr->compound_level == DEFAULT)
        scs_ptr->compound_mode = 1;
    else
        scs_ptr->compound_mode = scs_ptr->compound_level;
    if (scs_ptr->compound_mode) {
        scs_ptr->seq_header.order_hint_info.enable_jnt_comp = 1; //DISTANCE
        scs_ptr->seq_header.enable_masked_compound          = 1; //DIFF+WEDGE
    } else {
        scs_ptr->seq_header.order_hint_info.enable_jnt_comp = 0;
        scs_ptr->seq_header.enable_masked_compound          = 0;
    }

    if (scs_ptr->enable_intra_edge_filter == DEFAULT)
        scs_ptr->seq_header.enable_intra_edge_filter = 1;
    else
        scs_ptr->seq_header.enable_intra_edge_filter = (uint8_t)scs_ptr->enable_intra_edge_filter;

    if (scs_ptr->pic_based_rate_est == DEFAULT)
        scs_ptr->seq_header.pic_based_rate_est = 0;
    else
        scs_ptr->seq_header.pic_based_rate_est = (uint8_t)scs_ptr->pic_based_rate_est;

    if (scs_ptr->static_config.enable_restoration_filtering == DEFAULT) {
        scs_ptr->seq_header.enable_restoration = get_enable_restoration(
            scs_ptr->static_config.enc_mode,
            scs_ptr->static_config.enable_restoration_filtering,
            scs_ptr->input_resolution,
            scs_ptr->static_config.fast_decode);
    } else
        scs_ptr->seq_header.enable_restoration =
            (uint8_t)scs_ptr->static_config.enable_restoration_filtering;

    if (scs_ptr->static_config.cdef_level == DEFAULT)
        scs_ptr->seq_header.cdef_level = 1;
    else
        scs_ptr->seq_header.cdef_level = (uint8_t)(scs_ptr->static_config.cdef_level > 0);

    if (scs_ptr->enable_warped_motion == DEFAULT) {
        scs_ptr->seq_header.enable_warped_motion = 1;
    } else
        scs_ptr->seq_header.enable_warped_motion = (uint8_t)scs_ptr->enable_warped_motion;

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
    int64_t  target_fps         = (60 >> 16);

    int64_t buffer_threshold_1 = SC_FRAMES_INTERVAL_T1;
    int64_t buffer_threshold_2 = SC_FRAMES_INTERVAL_T2;
    int64_t buffer_threshold_3 = MIN(target_fps * 3, SC_FRAMES_INTERVAL_T3);
    svt_block_on_mutex(scs_ptr->encode_context_ptr->sc_buffer_mutex);

    if (scs_ptr->encode_context_ptr->sc_frame_in == 0)
        svt_av1_get_time(&context_ptr->first_in_pic_arrived_time_seconds,
                         &context_ptr->first_in_pic_arrived_timeu_seconds);
    else if (scs_ptr->encode_context_ptr->sc_frame_in == SC_FRAMES_TO_IGNORE)
        context_ptr->start_flag = TRUE;
    // Compute duration since the start of the encode and since the previous checkpoint
    svt_av1_get_time(&curs_time_seconds, &curs_time_useconds);

    overall_duration = svt_av1_compute_overall_elapsed_time_ms(
        context_ptr->first_in_pic_arrived_time_seconds,
        context_ptr->first_in_pic_arrived_timeu_seconds,
        curs_time_seconds,
        curs_time_useconds);

    inst_duration = svt_av1_compute_overall_elapsed_time_ms(context_ptr->prevs_time_seconds,
                                                            context_ptr->prevs_timeu_seconds,
                                                            curs_time_seconds,
                                                            curs_time_useconds);

    input_frames_count                     = (int64_t)overall_duration * (60 >> 16) / 1000;
    scs_ptr->encode_context_ptr->sc_buffer = input_frames_count -
        scs_ptr->encode_context_ptr->sc_frame_in;

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
        scs_ptr->encode_context_ptr->enc_mode = (EncMode)CLIP3(
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
        if (scs_ptr->encode_context_ptr->sc_frame_in >
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
        scs_ptr->encode_context_ptr->enc_mode = (EncMode)CLIP3(
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
        context_ptr->start_flag = FALSE;

        // Update previous stats
        context_ptr->previous_frame_in_check3 = scs_ptr->encode_context_ptr->sc_frame_in;
        context_ptr->prevs_time_seconds       = curs_time_seconds;
        context_ptr->prevs_timeu_seconds      = curs_time_useconds;
        context_ptr->prev_frame_out           = scs_ptr->encode_context_ptr->sc_frame_out;
    } else if (scs_ptr->encode_context_ptr->sc_frame_in < SC_FRAMES_TO_IGNORE &&
               (overall_duration != 0))
        context_ptr->cur_speed = (uint64_t)(scs_ptr->encode_context_ptr->sc_frame_out - 0) * 1000 /
            (uint64_t)(overall_duration);
    if (change_cond)
        context_ptr->prev_change_cond = change_cond;
    scs_ptr->encode_context_ptr->sc_frame_in++;
    if (scs_ptr->encode_context_ptr->sc_frame_in >= SC_FRAMES_TO_IGNORE)
        context_ptr->average_enc_mod += scs_ptr->encode_context_ptr->enc_mode;
    else
        context_ptr->average_enc_mod = 0;
    // Set the encoder level
    pcs_ptr->enc_mode = scs_ptr->encode_context_ptr->enc_mode;

    svt_release_mutex(scs_ptr->encode_context_ptr->sc_buffer_mutex);
    context_ptr->prev_enc_mod = scs_ptr->encode_context_ptr->enc_mode;
}
static EbErrorType reset_pcs_av1(PictureParentControlSet *pcs_ptr) {
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;
    Av1Common   *cm      = pcs_ptr->av1_cm;

    pcs_ptr->gf_interval = 0;

    pcs_ptr->reference_released   = 0;
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
    frm_hdr->quantization_params.delta_q_ac[AOM_PLANE_Y] = 0;
    frm_hdr->quantization_params.delta_q_dc[AOM_PLANE_Y] =
        pcs_ptr->scs_ptr->static_config.luma_y_dc_qindex_offset;
    frm_hdr->quantization_params.delta_q_ac[AOM_PLANE_U] =
        pcs_ptr->scs_ptr->static_config.chroma_u_ac_qindex_offset;
    frm_hdr->quantization_params.delta_q_dc[AOM_PLANE_U] =
        pcs_ptr->scs_ptr->static_config.chroma_u_dc_qindex_offset;
    frm_hdr->quantization_params.delta_q_ac[AOM_PLANE_V] =
        pcs_ptr->scs_ptr->static_config.chroma_v_ac_qindex_offset;
    frm_hdr->quantization_params.delta_q_dc[AOM_PLANE_V] =
        pcs_ptr->scs_ptr->static_config.chroma_v_dc_qindex_offset;

    // Encoder
    frm_hdr->quantization_params.using_qmatrix   = pcs_ptr->scs_ptr->static_config.enable_qm;
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
    frm_hdr->primary_ref_frame = PRIMARY_REF_NONE;
    if (pcs_ptr->scs_ptr->static_config.rate_control_mode == SVT_AV1_RC_MODE_CBR &&
        pcs_ptr->scs_ptr->static_config.intra_period_length != -1) {
        pcs_ptr->frame_offset = pcs_ptr->picture_number %
            (pcs_ptr->scs_ptr->static_config.intra_period_length + 1);
    } else
        pcs_ptr->frame_offset = pcs_ptr->picture_number;
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
    pcs_ptr->me_data_wrapper_ptr             = NULL;
    pcs_ptr->down_scaled_picture_wrapper_ptr = NULL;
    pcs_ptr->ds_pics.picture_ptr             = NULL;
    pcs_ptr->ds_pics.quarter_picture_ptr     = NULL;
    pcs_ptr->ds_pics.sixteenth_picture_ptr   = NULL;
    pcs_ptr->max_number_of_pus_per_sb        = SQUARE_PU_COUNT;

    atomic_set_u32(&pcs_ptr->pame_done, 0);

    svt_create_cond_var(&pcs_ptr->me_ready);

    SequenceControlSet *scs_ptr           = pcs_ptr->scs_ptr;
    pcs_ptr->me_segments_completion_count = 0;
    pcs_ptr->me_segments_column_count     = (uint8_t)(scs_ptr->me_segment_column_count_array[0]);
    pcs_ptr->me_segments_row_count        = (uint8_t)(scs_ptr->me_segment_row_count_array[0]);

    pcs_ptr->me_segments_total_count = (uint16_t)(pcs_ptr->me_segments_column_count *
                                                  pcs_ptr->me_segments_row_count);
    pcs_ptr->tpl_disp_coded_sb_count = 0;

    pcs_ptr->tpl_src_data_ready  = 0;
    pcs_ptr->tf_motion_direction = -1;

    return EB_ErrorNone;
}
/***********************************************
**** Copy the input buffer from the
**** sample application to the library buffers
************************************************/
static EbErrorType copy_frame_buffer_overlay(SequenceControlSet *scs_ptr, uint8_t *dst,
                                             uint8_t *src) {
    EbSvtAv1EncConfiguration *config       = &scs_ptr->static_config;
    EbErrorType               return_error = EB_ErrorNone;

    EbPictureBufferDesc *dst_picture_ptr = (EbPictureBufferDesc *)dst;
    EbPictureBufferDesc *src_picture_ptr = (EbPictureBufferDesc *)src;
    Bool                 is_16bit_input  = (Bool)(config->encoder_bit_depth > EB_EIGHT_BIT);

    // Need to include for Interlacing on the fly with pictureScanType = 1

    if (!is_16bit_input) {
        uint16_t input_row_index;
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
            svt_memcpy(
                (dst_picture_ptr->buffer_y + luma_buffer_offset + luma_stride * input_row_index),
                (src_picture_ptr->buffer_y + luma_buffer_offset + luma_stride * input_row_index),
                luma_width);
        }

        // U
        for (input_row_index = 0; input_row_index < (luma_height >> 1); input_row_index++) {
            svt_memcpy((dst_picture_ptr->buffer_cb + chroma_buffer_offset +
                        chroma_stride * input_row_index),
                       (src_picture_ptr->buffer_cb + chroma_buffer_offset +
                        chroma_stride * input_row_index),
                       chroma_width);
        }

        // V
        for (input_row_index = 0; input_row_index < (luma_height >> 1); input_row_index++) {
            svt_memcpy((dst_picture_ptr->buffer_cr + chroma_buffer_offset +
                        chroma_stride * input_row_index),
                       (src_picture_ptr->buffer_cr + chroma_buffer_offset +
                        chroma_stride * input_row_index),
                       chroma_width);
        }
    } else { // 10bit packed

        svt_memcpy(
            dst_picture_ptr->buffer_y, src_picture_ptr->buffer_y, src_picture_ptr->luma_size);

        svt_memcpy(
            dst_picture_ptr->buffer_cb, src_picture_ptr->buffer_cb, src_picture_ptr->chroma_size);

        svt_memcpy(
            dst_picture_ptr->buffer_cr, src_picture_ptr->buffer_cr, src_picture_ptr->chroma_size);

        svt_memcpy(dst_picture_ptr->buffer_bit_inc_y,
                   src_picture_ptr->buffer_bit_inc_y,
                   src_picture_ptr->luma_size >> 2);

        svt_memcpy(dst_picture_ptr->buffer_bit_inc_cb,
                   src_picture_ptr->buffer_bit_inc_cb,
                   src_picture_ptr->chroma_size >> 2);

        svt_memcpy(dst_picture_ptr->buffer_bit_inc_cr,
                   src_picture_ptr->buffer_bit_inc_cr,
                   src_picture_ptr->chroma_size >> 2);
    }
    return return_error;
}

/* overlay specific version of copy_input_buffer without passes specializations */
static void copy_input_buffer_overlay(SequenceControlSet *sequenceControlSet,
                                      EbBufferHeaderType *dst, EbBufferHeaderType *src) {
    // Copy the higher level structure
    dst->n_alloc_len  = src->n_alloc_len;
    dst->n_filled_len = src->n_filled_len;
    dst->flags        = src->flags;
    dst->pts          = src->pts;
    dst->n_tick_count = src->n_tick_count;
    dst->size         = src->size;
    dst->qp           = src->qp;
    dst->pic_type     = src->pic_type;

    // Copy the metadata array
    if (svt_aom_copy_metadata_buffer(dst, src->metadata) != EB_ErrorNone)
        dst->metadata = NULL;

    // Copy the picture buffer
    if (src->p_buffer != NULL)
        copy_frame_buffer_overlay(sequenceControlSet, dst->p_buffer, src->p_buffer);
}

/******************************************************
 * Read Stat from File
 ******************************************************/
void read_stat(SequenceControlSet *scs_ptr) {
    EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;

    encode_context_ptr->rc_stats_buffer = scs_ptr->static_config.rc_stats_buffer;
}
void setup_two_pass(SequenceControlSet *scs_ptr) {
    EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
    scs_ptr->twopass.passes           = scs_ptr->passes;
    scs_ptr->twopass.stats_buf_ctx    = &encode_context_ptr->stats_buf_context;
    scs_ptr->twopass.stats_in         = scs_ptr->twopass.stats_buf_ctx->stats_in_start;
    if (scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
        scs_ptr->static_config.pass == ENC_LAST_PASS) {
        const size_t packet_sz = sizeof(FIRSTPASS_STATS);
        const int    packets   = (int)(encode_context_ptr->rc_stats_buffer.sz / packet_sz);

        if (!scs_ptr->lap_rc) {
            /*Re-initialize to stats buffer, populated by application in the case of
                * two pass*/
            scs_ptr->twopass.stats_buf_ctx->stats_in_start =
                encode_context_ptr->rc_stats_buffer.buf;
            scs_ptr->twopass.stats_in = scs_ptr->twopass.stats_buf_ctx->stats_in_start;
            scs_ptr->twopass.stats_buf_ctx->stats_in_end_write =
                &scs_ptr->twopass.stats_buf_ctx->stats_in_start[packets - 1];
            scs_ptr->twopass.stats_buf_ctx->stats_in_end =
                &scs_ptr->twopass.stats_buf_ctx->stats_in_start[packets - 1];
            svt_av1_init_second_pass(scs_ptr);
            //less than 200 frames or gop_constraint_rc, used in VBR and set in multipass encode
            scs_ptr->is_short_clip = scs_ptr->twopass.stats_buf_ctx->total_stats->count < 200
                ? 1
                : scs_ptr->is_short_clip;
        }
    } else if (scs_ptr->lap_rc)
        svt_av1_init_single_pass_lap(scs_ptr);
}

extern EbErrorType first_pass_signal_derivation_pre_analysis_pcs(PictureParentControlSet *pcs_ptr);
extern EbErrorType first_pass_signal_derivation_pre_analysis_scs(SequenceControlSet *scs_ptr);

static EbErrorType realloc_sb_param(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {
    EB_FREE_ARRAY(pcs_ptr->b64_geom);
    EB_MALLOC_ARRAY(pcs_ptr->b64_geom, scs_ptr->b64_total_count);
    memcpy(pcs_ptr->b64_geom, scs_ptr->b64_geom, sizeof(B64Geom) * scs_ptr->b64_total_count);
    EB_FREE_ARRAY(pcs_ptr->sb_geom);
    EB_MALLOC_ARRAY(pcs_ptr->sb_geom, scs_ptr->sb_total_count);
    memcpy(pcs_ptr->sb_geom, scs_ptr->sb_geom, sizeof(SbGeom) * scs_ptr->sb_total_count);
    pcs_ptr->is_pcs_sb_params = TRUE;
    return EB_ErrorNone;
}

static void update_frame_event(PictureParentControlSet *pcs_ptr) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    EbPrivDataNode     *node    = (EbPrivDataNode *)pcs_ptr->input_ptr->p_app_private;
    while (node) {
        if (node->node_type == REF_FRAME_SCALING_EVENT) {
            // update resize denominator by input event
            assert_err(node->size == sizeof(EbRefFrameScale),
                       "private data size mismatch of REF_FRAME_SCALING_EVENT");
            // update scaling event for future pictures
            scs_ptr->encode_context_ptr->resize_evt = *(EbRefFrameScale *)node->data;
            // set reset flag of rate control
            pcs_ptr->rc_reset_flag = TRUE;
        }
        node = node->next;
    }
    // update current picture scaling event
    pcs_ptr->resize_evt = scs_ptr->encode_context_ptr->resize_evt;
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
    EbThreadContext             *enc_contxt_ptr = (EbThreadContext *)input_ptr;
    ResourceCoordinationContext *context_ptr = (ResourceCoordinationContext *)enc_contxt_ptr->priv;

    EbObjectWrapper *pcs_wrapper_ptr;

    PictureParentControlSet *pcs_ptr;
    SequenceControlSet      *scs_ptr;

    EbObjectWrapper             *eb_input_wrapper_ptr;
    EbBufferHeaderType          *eb_input_ptr;
    EbObjectWrapper             *output_wrapper_ptr;
    ResourceCoordinationResults *out_results_ptr;
    EbObjectWrapper             *eb_input_cmd_wrapper;
    InputCommand                *input_cmd_obj;
    EbObjectWrapper             *input_picture_wrapper_ptr;
    EbObjectWrapper             *reference_picture_wrapper_ptr;

    Bool             end_of_sequence_flag = FALSE;
    EbObjectWrapper *prev_pcs_wrapper_ptr = 0;

    for (;;) {
        // Tie instance_index to zero for now...
        uint32_t instance_index = 0;
        // Get the input command containing 2 input buffers: y8b & rest(uv8b+yuvbitInc)
        EB_GET_FULL_OBJECT(context_ptr->input_cmd_fifo_ptr, &eb_input_cmd_wrapper);

        input_cmd_obj = (InputCommand *)eb_input_cmd_wrapper->object_ptr;

        EbObjectWrapper    *eb_y8b_wrapper_ptr = input_cmd_obj->eb_y8b_wrapper_ptr;
        EbBufferHeaderType *y8b_header = (EbBufferHeaderType *)eb_y8b_wrapper_ptr->object_ptr;
        uint8_t            *buff_y8b   = ((EbPictureBufferDesc *)y8b_header->p_buffer)->buffer_y;
        eb_input_wrapper_ptr           = input_cmd_obj->eb_input_wrapper_ptr;
        eb_input_ptr                   = (EbBufferHeaderType *)eb_input_wrapper_ptr->object_ptr;

        // Set the SequenceControlSet
        scs_ptr = context_ptr->scs_instance_array[instance_index]->scs_ptr;

        if (context_ptr->scs_instance_array[instance_index]->encode_context_ptr->initial_picture) {
            // Update picture width, picture height, cropping right offset, cropping bottom offset, and conformance windows
            scs_ptr->chroma_width  = (scs_ptr->max_input_luma_width >> 1);
            scs_ptr->chroma_height = (scs_ptr->max_input_luma_height >> 1);

            scs_ptr->pad_right  = scs_ptr->max_input_pad_right;
            scs_ptr->pad_bottom = scs_ptr->max_input_pad_bottom;

            // Pre-Analysis Signal(s) derivation
            if (scs_ptr->static_config.pass == ENC_FIRST_PASS)
                first_pass_signal_derivation_pre_analysis_scs(scs_ptr);
            else
                signal_derivation_pre_analysis_oq_scs(scs_ptr);

            // Init SB Params
            const uint32_t input_size = scs_ptr->max_input_luma_width *
                scs_ptr->max_input_luma_height;
            derive_input_resolution(&scs_ptr->input_resolution, input_size);

            b64_geom_init(scs_ptr);
            sb_geom_init(scs_ptr);

            // sf_identity
            svt_av1_setup_scale_factors_for_frame(&scs_ptr->sf_identity,
                                                  scs_ptr->max_input_luma_width,
                                                  scs_ptr->max_input_luma_height,
                                                  scs_ptr->max_input_luma_width,
                                                  scs_ptr->max_input_luma_height);
        }
        // Since at this stage we do not know the prediction structure and the location of ALT_REF pictures,
        // for every picture (except first picture), we allocate two: 1. original picture, 2. potential Overlay picture.
        // In Picture Decision Process, where the overlay frames are known, they extra pictures are released
        uint8_t has_overlay =
            (scs_ptr->static_config.enable_overlays == FALSE ||
             context_ptr->scs_instance_array[instance_index]->encode_context_ptr->initial_picture)
            ? 0
            : 1;
        for (uint8_t loop_index = 0; loop_index <= has_overlay && !end_of_sequence_flag;
             loop_index++) {
            //Get a New ParentPCS where we will hold the new input_picture
            svt_get_empty_object(context_ptr->picture_control_set_fifo_ptr_array[instance_index],
                                 &pcs_wrapper_ptr);

            // Parent PCS is released by the Rate Control after passing through MDC->MD->ENCDEC->Packetization
            svt_object_inc_live_count(pcs_wrapper_ptr, 1);

            pcs_ptr = (PictureParentControlSet *)pcs_wrapper_ptr->object_ptr;

            // - p_pcs_wrapper_ptr is a direct copy of pcs_wrapper_ptr (live_count == 1).
            // - Most of p_pcs_wrapper_ptr in pre-allocated overlay candidates will be released & recycled to empty fifo
            //     by altref candidate's svt_release_object(pcs_ptr->overlay_ppcs_ptr->p_pcs_wrapper_ptr) in PictureDecision.
            // - The recycled ppcs may be assigned a new picture_number in ResourceCoordination.
            // - If the to-be-removed overlay candidate runs in picture_decision_kernel() after above release/recycle/assign,
            //     picture_decision_reorder_queue will update by the same picture_number (of the same ppcs ptr) twice and CHECK_REPORT_ERROR_NC occur.
            // - So need ppcs live_count + 1 before post ResourceCoordinationResults, and release ppcs before end of PictureDecision,
            //     to avoid recycling overlay candidate's ppcs to empty fifo too early.
            pcs_ptr->p_pcs_wrapper_ptr = pcs_wrapper_ptr;

            // reallocate sb_param_array and sb_geom for super-res or reference scaling mode on
            if (scs_ptr->static_config.superres_mode > SUPERRES_NONE ||
                scs_ptr->static_config.resize_mode > RESIZE_NONE)
                realloc_sb_param(scs_ptr, pcs_ptr);
            else {
                pcs_ptr->b64_geom         = scs_ptr->b64_geom;
                pcs_ptr->sb_geom          = scs_ptr->sb_geom;
                pcs_ptr->is_pcs_sb_params = FALSE;
            }
            pcs_ptr->input_resolution  = scs_ptr->input_resolution;
            pcs_ptr->picture_sb_width  = scs_ptr->pic_width_in_b64;
            pcs_ptr->picture_sb_height = scs_ptr->pic_height_in_b64;

            pcs_ptr->overlay_ppcs_ptr   = NULL;
            pcs_ptr->is_alt_ref         = 0;
            pcs_ptr->transition_present = -1;
            pcs_ptr->is_noise_level     = 0;
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
            //make pcs input buffer access the luma8bit part from the Luma8bit Pool
            pcs_ptr->enhanced_picture_ptr->buffer_y = buff_y8b;
            pcs_ptr->input_ptr                      = eb_input_ptr;
            end_of_sequence_flag = (pcs_ptr->input_ptr->flags & EB_BUFFERFLAG_EOS) ? TRUE : FALSE;
            // Check whether super-res is previously enabled in this recycled parent pcs and restore to non-scale-down default if so.
            if (pcs_ptr->frame_superres_enabled || pcs_ptr->frame_resize_enabled)
                reset_resized_picture(scs_ptr, pcs_ptr, pcs_ptr->enhanced_picture_ptr);
            pcs_ptr->superres_total_recode_loop = 0;
            pcs_ptr->superres_recode_loop       = 0;
            svt_av1_get_time(&pcs_ptr->start_time_seconds, &pcs_ptr->start_time_u_seconds);
            pcs_ptr->scs_ptr                   = scs_ptr;
            pcs_ptr->input_picture_wrapper_ptr = input_picture_wrapper_ptr;
            //store the y8b warapper to be used for release later
            pcs_ptr->eb_y8b_wrapper_ptr   = eb_y8b_wrapper_ptr;
            pcs_ptr->end_of_sequence_flag = end_of_sequence_flag;
            pcs_ptr->rc_reset_flag        = FALSE;
            update_frame_event(pcs_ptr);
            pcs_ptr->is_not_scaled = (scs_ptr->static_config.superres_mode == SUPERRES_NONE) &&
                scs_ptr->static_config.resize_mode == RESIZE_NONE;
            if (loop_index == 1) {
                // Get a new input picture for overlay.
                EbObjectWrapper *input_pic_wrapper_ptr;

                // Get a new input picture for overlay.
                svt_get_empty_object(
                    scs_ptr->encode_context_ptr->overlay_input_picture_pool_fifo_ptr,
                    &input_pic_wrapper_ptr);

                // Copy from original picture (pcs_ptr->input_picture_wrapper_ptr), which is shared between overlay and alt_ref up to this point, to the new input picture.
                if (pcs_ptr->alt_ref_ppcs_ptr->input_picture_wrapper_ptr->object_ptr != NULL) {
                    copy_input_buffer_overlay(
                        scs_ptr,
                        (EbBufferHeaderType *)input_pic_wrapper_ptr->object_ptr,
                        (EbBufferHeaderType *)
                            pcs_ptr->alt_ref_ppcs_ptr->input_picture_wrapper_ptr->object_ptr);
                }
                // Assign the new picture to the new pointers
                pcs_ptr->input_ptr = (EbBufferHeaderType *)input_pic_wrapper_ptr->object_ptr;
                pcs_ptr->enhanced_picture_ptr = (EbPictureBufferDesc *)pcs_ptr->input_ptr->p_buffer;
                pcs_ptr->input_picture_wrapper_ptr = input_pic_wrapper_ptr;

                // overlay does NOT use y8b buffer, set to NULL to avoid eb_y8b_wrapper_ptr->live_count disorder
                pcs_ptr->eb_y8b_wrapper_ptr = NULL;
            }
            // Set Picture Control Flags
            pcs_ptr->idr_flag          = scs_ptr->encode_context_ptr->initial_picture;
            pcs_ptr->cra_flag          = 0;
            pcs_ptr->scene_change_flag = FALSE;
            pcs_ptr->qp_on_the_fly     = FALSE;
            pcs_ptr->b64_total_count   = scs_ptr->b64_total_count;
            if (scs_ptr->speed_control_flag) {
                speed_buffer_control(context_ptr, pcs_ptr, scs_ptr);
            } else
                pcs_ptr->enc_mode = (EncMode)scs_ptr->static_config.enc_mode;
            //  If the mode of the second pass is not set from CLI, it is set to enc_mode

            // Pre-Analysis Signal(s) derivation
            if (scs_ptr->static_config.pass == ENC_FIRST_PASS)
                first_pass_signal_derivation_pre_analysis_pcs(pcs_ptr);
            else
                signal_derivation_pre_analysis_oq_pcs(pcs_ptr);
            // Rate Control

            // Picture Stats
            if (loop_index == has_overlay || end_of_sequence_flag)
                pcs_ptr->picture_number = context_ptr->picture_number_array[instance_index]++;
            else
                pcs_ptr->picture_number = context_ptr->picture_number_array[instance_index];
            if (pcs_ptr->picture_number == 0) {
                if (scs_ptr->static_config.pass == ENC_MIDDLE_PASS ||
                    scs_ptr->static_config.pass == ENC_LAST_PASS)
                    read_stat(scs_ptr);
                if (scs_ptr->static_config.pass != ENC_SINGLE_PASS || scs_ptr->lap_rc)
                    setup_two_pass(scs_ptr);
                else
                    set_rc_param(scs_ptr);
                if (scs_ptr->static_config.pass == ENC_MIDDLE_PASS)
                    find_init_qp_middle_pass(scs_ptr, pcs_ptr);
            }
            if (scs_ptr->passes == 3 && !end_of_sequence_flag &&
                scs_ptr->static_config.pass == ENC_LAST_PASS &&
                scs_ptr->static_config.rate_control_mode) {
                pcs_ptr->stat_struct = (scs_ptr->twopass.stats_buf_ctx->stats_in_start +
                                        pcs_ptr->picture_number)
                                           ->stat_struct;
                if (pcs_ptr->stat_struct.poc != pcs_ptr->picture_number)
                    SVT_LOG("Error reading data in multi pass encoding\n");
            }
            if (scs_ptr->static_config.use_qp_file == 1) {
                pcs_ptr->qp_on_the_fly = TRUE;
                if (pcs_ptr->input_ptr->qp > MAX_QP_VALUE) {
                    SVT_WARN("INPUT QP/CRF OUTSIDE OF RANGE\n");
                    pcs_ptr->qp_on_the_fly = FALSE;
                }
                pcs_ptr->picture_qp = (uint8_t)pcs_ptr->input_ptr->qp;
            } else {
                pcs_ptr->qp_on_the_fly = FALSE;
                pcs_ptr->picture_qp    = (uint8_t)scs_ptr->static_config.qp;
            }

            pcs_ptr->ts_duration = (double)10000000 * (1 << 16) / scs_ptr->frame_rate;
            scs_ptr->encode_context_ptr->initial_picture = FALSE;

            // Get Empty Reference Picture Object
            svt_get_empty_object(scs_ptr->encode_context_ptr->pa_reference_picture_pool_fifo_ptr,
                                 &reference_picture_wrapper_ptr);

            pcs_ptr->pa_reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
            //make pa_ref full sample buffer access the luma8bit part from the y8b Pool
            EbPaReferenceObject *pa_ref_obj =
                (EbPaReferenceObject *)pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
            EbPictureBufferDesc *input_padded_picture_ptr =
                (EbPictureBufferDesc *)pa_ref_obj->input_padded_picture_ptr;
            input_padded_picture_ptr->buffer_y = buff_y8b;
            // Since overlay pictures are not added to PA_Reference queue in PD and not released there, the life count is only set to 1
            if (pcs_ptr->is_overlay)
                // Give the new Reference a nominal live_count of 1
                svt_object_inc_live_count(pcs_ptr->pa_reference_picture_wrapper_ptr, 1);
            else
                svt_object_inc_live_count(pcs_ptr->pa_reference_picture_wrapper_ptr, 2);
            if (pcs_ptr->eb_y8b_wrapper_ptr) {
                // y8b follows longest life cycle of pa ref and input. so it needs to build on top of live count of pa ref
                svt_object_inc_live_count(pcs_ptr->eb_y8b_wrapper_ptr, 2);
            }
            if (scs_ptr->static_config.restricted_motion_vector) {
                struct PictureParentControlSet *ppcs_ptr = pcs_ptr;
                Av1Common *const                cm       = ppcs_ptr->av1_cm;
                uint8_t                         pic_width_in_sb =
                    (uint8_t)((pcs_ptr->aligned_width + scs_ptr->sb_size - 1) / scs_ptr->sb_size);
                int       tile_row, tile_col;
                uint32_t  x_sb_index, y_sb_index;
                const int tile_cols = cm->tiles_info.tile_cols;
                const int tile_rows = cm->tiles_info.tile_rows;
                TileInfo  tile_info;
                int       sb_size_log2 = scs_ptr->seq_header.sb_size_log2;
                //Tile Loop
                for (tile_row = 0; tile_row < tile_rows; tile_row++) {
                    svt_av1_tile_set_row(&tile_info, &cm->tiles_info, cm->mi_rows, tile_row);

                    for (tile_col = 0; tile_col < tile_cols; tile_col++) {
                        svt_av1_tile_set_col(&tile_info, &cm->tiles_info, cm->mi_cols, tile_col);

                        for ((y_sb_index =
                                  cm->tiles_info.tile_row_start_mi[tile_row] >> sb_size_log2);
                             (y_sb_index <
                              ((uint32_t)cm->tiles_info.tile_row_start_mi[tile_row + 1] >>
                               sb_size_log2));
                             y_sb_index++) {
                            for ((x_sb_index =
                                      cm->tiles_info.tile_col_start_mi[tile_col] >> sb_size_log2);
                                 (x_sb_index <
                                  ((uint32_t)cm->tiles_info.tile_col_start_mi[tile_col + 1] >>
                                   sb_size_log2));
                                 x_sb_index++) {
                                int sb_index                             = (uint16_t)(x_sb_index +
                                                          y_sb_index * pic_width_in_sb);
                                scs_ptr->b64_geom[sb_index].tile_start_x = 4 *
                                    tile_info.mi_col_start;
                                scs_ptr->b64_geom[sb_index].tile_end_x   = 4 * tile_info.mi_col_end;
                                scs_ptr->b64_geom[sb_index].tile_start_y = 4 *
                                    tile_info.mi_row_start;
                                scs_ptr->b64_geom[sb_index].tile_end_y = 4 * tile_info.mi_row_end;
                            }
                        }
                    }
                }
            }

            // Get Empty Output Results Object
            if (pcs_ptr->picture_number > 0 && (prev_pcs_wrapper_ptr != NULL)) {
                PictureParentControlSet *ppcs_out = (PictureParentControlSet *)
                                                        prev_pcs_wrapper_ptr->object_ptr;

                ppcs_out->end_of_sequence_flag = end_of_sequence_flag;
                // since overlay frame has the end of sequence set properly, set the end of sequence to true in the alt ref picture
                if (ppcs_out->is_overlay && end_of_sequence_flag)
                    ppcs_out->alt_ref_ppcs_ptr->end_of_sequence_flag = TRUE;

                reset_pcs_av1(ppcs_out);

                svt_get_empty_object(context_ptr->resource_coordination_results_output_fifo_ptr,
                                     &output_wrapper_ptr);
                out_results_ptr = (ResourceCoordinationResults *)output_wrapper_ptr->object_ptr;

                if (scs_ptr->static_config.enable_overlays == TRUE) {
                    // ppcs live_count + 1 for PictureAnalysis & PictureDecision, will svt_release_object(ppcs) at the end of picture_decision_kernel.
                    svt_object_inc_live_count(prev_pcs_wrapper_ptr, 1);
                }

                out_results_ptr->pcs_wrapper_ptr = prev_pcs_wrapper_ptr;
                // Post the finished Results Object
                svt_post_full_object(output_wrapper_ptr);
            }
            prev_pcs_wrapper_ptr = pcs_wrapper_ptr;
        }
        // Release the Input Command
        svt_release_object(eb_input_cmd_wrapper);
    }

    return NULL;
}
