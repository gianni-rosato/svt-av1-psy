/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#include "EbEncHandle.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbMotionEstimationResults.h"
#include "EbInitialRateControlProcess.h"
#include "EbInitialRateControlResults.h"
#include "EbMotionEstimationContext.h"
#include "EbUtility.h"
#include "EbReferenceObject.h"
#include "EbResize.h"
#include "common_dsp_rtcd.h"
/**************************************
 * Context
 **************************************/
typedef struct InitialRateControlContext {
    EbFifo *motion_estimation_results_input_fifo_ptr;
    EbFifo *initialrate_control_results_output_fifo_ptr;
} InitialRateControlContext;

/**************************************
* Macros
**************************************/
#define PAN_SB_PERCENTAGE 75
#define LOW_AMPLITUDE_TH 16

static void eb_get_mv(PictureParentControlSet *pcs_ptr, uint32_t sb_index, int32_t *x_current_mv,
            int32_t *y_current_mv) {
    uint32_t me_candidate_index;

    const MeSbResults *me_results       = pcs_ptr->me_results[sb_index];
    uint8_t            total_me_cnt     = me_results->total_me_candidate_index[0];
    const MeCandidate *me_block_results = me_results->me_candidate[0];
    for (me_candidate_index = 0; me_candidate_index < total_me_cnt; me_candidate_index++) {
        if (me_block_results->direction == UNI_PRED_LIST_0) {
            *x_current_mv = me_results->me_mv_array[0][0].x_mv;
            *y_current_mv = me_results->me_mv_array[0][0].y_mv;
            break;
        }
    }
}
EbBool check_mv_for_pan_high_amp(uint32_t hierarchical_levels, uint32_t temporal_layer_index,
                                 int32_t *x_current_mv, int32_t *x_candidate_mv) {
    if (*x_current_mv * *x_candidate_mv >
            0 // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*x_current_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*x_candidate_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*x_current_mv - *x_candidate_mv) < LOW_AMPLITUDE_TH) { // close amplitude

        return (EB_TRUE);
    }

    else
        return (EB_FALSE);
}

EbBool check_mv_for_tilt_high_amp(uint32_t hierarchical_levels, uint32_t temporal_layer_index,
                                  int32_t *y_current_mv, int32_t *y_candidate_mv) {
    if (*y_current_mv * *y_candidate_mv >
            0 // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*y_current_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*y_candidate_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*y_current_mv - *y_candidate_mv) < LOW_AMPLITUDE_TH) { // close amplitude

        return (EB_TRUE);
    }

    else
        return (EB_FALSE);
}

EbBool check_mv_for_pan(uint32_t hierarchical_levels, uint32_t temporal_layer_index,
                        int32_t *x_current_mv, int32_t *y_current_mv, int32_t *x_candidate_mv,
                        int32_t *y_candidate_mv) {
    if (*y_current_mv < LOW_AMPLITUDE_TH &&
        *y_candidate_mv<
            LOW_AMPLITUDE_TH && * x_current_mv * *
            x_candidate_mv > 0 // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*x_current_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*x_candidate_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*x_current_mv - *x_candidate_mv) < LOW_AMPLITUDE_TH) { // close amplitude

        return (EB_TRUE);
    }

    else
        return (EB_FALSE);
}

EbBool check_mv_for_tilt(uint32_t hierarchical_levels, uint32_t temporal_layer_index,
                         int32_t *x_current_mv, int32_t *y_current_mv, int32_t *x_candidate_mv,
                         int32_t *y_candidate_mv) {
    if (*x_current_mv < LOW_AMPLITUDE_TH &&
        *x_candidate_mv<
            LOW_AMPLITUDE_TH && * y_current_mv * *
            y_candidate_mv > 0 // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*y_current_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*y_candidate_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*y_current_mv - *y_candidate_mv) < LOW_AMPLITUDE_TH) { // close amplitude

        return (EB_TRUE);
    }

    else
        return (EB_FALSE);
}

EbBool check_mv_for_non_uniform_motion(int32_t *x_current_mv, int32_t *y_current_mv,
                                       int32_t *x_candidate_mv, int32_t *y_candidate_mv) {
    int32_t mv_threshold = 40; //LOW_AMPLITUDE_TH + 18;
    // Either the x or the y direction is greater than threshold
    if ((ABS(*x_current_mv - *x_candidate_mv) > mv_threshold) ||
        (ABS(*y_current_mv - *y_candidate_mv) > mv_threshold))
        return (EB_TRUE);
    else
        return (EB_FALSE);
}

void check_for_non_uniform_motion_vector_field(PictureParentControlSet *pcs_ptr) {
    uint32_t sb_count;
    uint32_t pic_width_in_sb =
        (pcs_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t sb_origin_x;
    uint32_t sb_origin_y;

    int32_t  x_current_mv                   = 0;
    int32_t  y_current_mv                   = 0;
    int32_t  x_left_mv                      = 0;
    int32_t  y_left_mv                      = 0;
    int32_t  x_top_mv                       = 0;
    int32_t  y_top_mv                       = 0;
    int32_t  x_right_mv                     = 0;
    int32_t  y_right_mv                       = 0;
    int32_t  x_bottom_mv                    = 0;
    int32_t  y_bottom_mv                    = 0;
    uint32_t count_of_non_uniform_neighbors = 0;

    for (sb_count = 0; sb_count < pcs_ptr->sb_total_count; ++sb_count) {
        count_of_non_uniform_neighbors = 0;

        sb_origin_x = (sb_count % pic_width_in_sb) * BLOCK_SIZE_64;
        sb_origin_y = (sb_count / pic_width_in_sb) * BLOCK_SIZE_64;

        if (((sb_origin_x + BLOCK_SIZE_64) <= pcs_ptr->enhanced_picture_ptr->width) &&
            ((sb_origin_y + BLOCK_SIZE_64) <= pcs_ptr->enhanced_picture_ptr->height)) {
            // Current MV
            eb_get_mv(pcs_ptr, sb_count, &x_current_mv, &y_current_mv);

            // Left MV
            if (sb_origin_x == 0) {
                x_left_mv = 0;
                y_left_mv = 0;
            } else
                eb_get_mv(pcs_ptr, sb_count - 1, &x_left_mv, &y_left_mv);
            count_of_non_uniform_neighbors += check_mv_for_non_uniform_motion(
                &x_current_mv, &y_current_mv, &x_left_mv, &y_left_mv);

            // Top MV
            if (sb_origin_y == 0) {
                x_top_mv = 0;
                y_top_mv = 0;
            } else
                eb_get_mv(pcs_ptr, sb_count - pic_width_in_sb, &x_top_mv, &y_top_mv);
            count_of_non_uniform_neighbors +=
                check_mv_for_non_uniform_motion(&x_current_mv, &y_current_mv, &x_top_mv, &y_top_mv);

            // Right MV
            if ((sb_origin_x + (BLOCK_SIZE_64 << 1)) > pcs_ptr->enhanced_picture_ptr->width) {
                x_right_mv = 0;
                y_right_mv   = 0;
            } else
                eb_get_mv(pcs_ptr, sb_count + 1, &x_right_mv, &y_right_mv);
            count_of_non_uniform_neighbors += check_mv_for_non_uniform_motion(
                &x_current_mv, &y_current_mv, &x_right_mv, &y_right_mv);

            // Bottom MV
            if ((sb_origin_y + (BLOCK_SIZE_64 << 1)) > pcs_ptr->enhanced_picture_ptr->height) {
                x_bottom_mv = 0;
                y_bottom_mv = 0;
            } else
                eb_get_mv(pcs_ptr, sb_count + pic_width_in_sb, &x_bottom_mv, &y_bottom_mv);
            count_of_non_uniform_neighbors += check_mv_for_non_uniform_motion(
                &x_current_mv, &y_current_mv, &x_bottom_mv, &y_bottom_mv);
        }
    }
}

void detect_global_motion(PictureParentControlSet *pcs_ptr) {
    //initilize global motion to be OFF for all references frames.
    memset(pcs_ptr->is_global_motion, EB_FALSE, MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);

    if (pcs_ptr->gm_level <= GM_DOWN) {
        uint32_t num_of_list_to_search =
            (pcs_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;

        for (uint32_t list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
            uint32_t num_of_ref_pic_to_search;
            if (pcs_ptr->is_alt_ref == EB_TRUE)
                num_of_ref_pic_to_search = 1;
            else
                num_of_ref_pic_to_search = pcs_ptr->slice_type == P_SLICE
                                               ? pcs_ptr->ref_list0_count
                                               : list_index == REF_LIST_0
                                                     ? pcs_ptr->ref_list0_count
                                                     : pcs_ptr->ref_list1_count;

            // Ref Picture Loop
            for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
                 ++ref_pic_index) {
                pcs_ptr->is_global_motion[list_index][ref_pic_index] = EB_FALSE;
                if (pcs_ptr->global_motion_estimation[list_index][ref_pic_index].wmtype >
                    TRANSLATION)
                    pcs_ptr->is_global_motion[list_index][ref_pic_index] = EB_TRUE;
            }
        }
    } else {
        uint32_t sb_count;
        uint32_t pic_width_in_sb =
            (pcs_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
        uint32_t sb_origin_x;
        uint32_t sb_origin_y;

        uint32_t total_checked_sbs = 0;
        uint32_t total_pan_sbs     = 0;

        int32_t  x_current_mv   = 0;
        int32_t  y_current_mv   = 0;
        int32_t  x_left_mv      = 0;
        int32_t  y_left_mv      = 0;
        int32_t  x_top_mv       = 0;
        int32_t  y_top_mv       = 0;
        int32_t  x_right_mv     = 0;
        int32_t  y_right_mv       = 0;
        int32_t  x_bottom_mv    = 0;
        int32_t  y_bottom_mv    = 0;
        int64_t  x_tile_mv_sum  = 0;
        int64_t  y_tilt_mv_sum  = 0;
        int64_t  x_pan_mv_sum   = 0;
        int64_t  y_pan_mv_sum   = 0;
        uint32_t total_tilt_sbs = 0;

        uint32_t total_tilt_high_amp_sbs = 0;
        uint32_t total_pan_high_amp_sbs  = 0;

        for (sb_count = 0; sb_count < pcs_ptr->sb_total_count; ++sb_count) {
            sb_origin_x = (sb_count % pic_width_in_sb) * BLOCK_SIZE_64;
            sb_origin_y = (sb_count / pic_width_in_sb) * BLOCK_SIZE_64;
            if (((sb_origin_x + BLOCK_SIZE_64) <= pcs_ptr->enhanced_picture_ptr->width) &&
                ((sb_origin_y + BLOCK_SIZE_64) <= pcs_ptr->enhanced_picture_ptr->height)) {
                // Current MV
                eb_get_mv(pcs_ptr, sb_count, &x_current_mv, &y_current_mv);

                // Left MV
                if (sb_origin_x == 0) {
                    x_left_mv = 0;
                    y_left_mv = 0;
                } else
                    eb_get_mv(pcs_ptr, sb_count - 1, &x_left_mv, &y_left_mv);
                // Top MV
                if (sb_origin_y == 0) {
                    x_top_mv = 0;
                    y_top_mv = 0;
                } else
                    eb_get_mv(pcs_ptr, sb_count - pic_width_in_sb, &x_top_mv, &y_top_mv);
                // Right MV
                if ((sb_origin_x + (BLOCK_SIZE_64 << 1)) > pcs_ptr->enhanced_picture_ptr->width) {
                    x_right_mv = 0;
                    y_right_mv   = 0;
                } else
                    eb_get_mv(pcs_ptr, sb_count + 1, &x_right_mv, &y_right_mv);
                // Bottom MV
                if ((sb_origin_y + (BLOCK_SIZE_64 << 1)) > pcs_ptr->enhanced_picture_ptr->height) {
                    x_bottom_mv = 0;
                    y_bottom_mv = 0;
                } else
                    eb_get_mv(pcs_ptr, sb_count + pic_width_in_sb, &x_bottom_mv, &y_bottom_mv);
                total_checked_sbs++;

                if ((EbBool)(check_mv_for_pan(pcs_ptr->hierarchical_levels,
                                              pcs_ptr->temporal_layer_index,
                                              &x_current_mv,
                                              &y_current_mv,
                                              &x_left_mv,
                                              &y_left_mv) ||
                             check_mv_for_pan(pcs_ptr->hierarchical_levels,
                                              pcs_ptr->temporal_layer_index,
                                              &x_current_mv,
                                              &y_current_mv,
                                              &x_top_mv,
                                              &y_top_mv) ||
                             check_mv_for_pan(pcs_ptr->hierarchical_levels,
                                              pcs_ptr->temporal_layer_index,
                                              &x_current_mv,
                                              &y_current_mv,
                                              &x_right_mv,
                                              &y_right_mv) ||
                             check_mv_for_pan(pcs_ptr->hierarchical_levels,
                                              pcs_ptr->temporal_layer_index,
                                              &x_current_mv,
                                              &y_current_mv,
                                              &x_bottom_mv,
                                              &y_bottom_mv))) {
                    total_pan_sbs++;

                    x_pan_mv_sum += x_current_mv;
                    y_pan_mv_sum += y_current_mv;
                }

                if ((EbBool)(check_mv_for_tilt(pcs_ptr->hierarchical_levels,
                                               pcs_ptr->temporal_layer_index,
                                               &x_current_mv,
                                               &y_current_mv,
                                               &x_left_mv,
                                               &y_left_mv) ||
                             check_mv_for_tilt(pcs_ptr->hierarchical_levels,
                                               pcs_ptr->temporal_layer_index,
                                               &x_current_mv,
                                               &y_current_mv,
                                               &x_top_mv,
                                               &y_top_mv) ||
                             check_mv_for_tilt(pcs_ptr->hierarchical_levels,
                                               pcs_ptr->temporal_layer_index,
                                               &x_current_mv,
                                               &y_current_mv,
                                               &x_right_mv,
                                               &y_right_mv) ||
                             check_mv_for_tilt(pcs_ptr->hierarchical_levels,
                                               pcs_ptr->temporal_layer_index,
                                               &x_current_mv,
                                               &y_current_mv,
                                               &x_bottom_mv,
                                               &y_bottom_mv))) {
                    total_tilt_sbs++;

                    x_tile_mv_sum += x_current_mv;
                    y_tilt_mv_sum += y_current_mv;
                }

                if ((EbBool)(check_mv_for_pan_high_amp(pcs_ptr->hierarchical_levels,
                                                       pcs_ptr->temporal_layer_index,
                                                       &x_current_mv,
                                                       &x_left_mv) ||
                             check_mv_for_pan_high_amp(pcs_ptr->hierarchical_levels,
                                                       pcs_ptr->temporal_layer_index,
                                                       &x_current_mv,
                                                       &x_top_mv) ||
                             check_mv_for_pan_high_amp(pcs_ptr->hierarchical_levels,
                                                       pcs_ptr->temporal_layer_index,
                                                       &x_current_mv,
                                                       &x_right_mv) ||
                             check_mv_for_pan_high_amp(pcs_ptr->hierarchical_levels,
                                                       pcs_ptr->temporal_layer_index,
                                                       &x_current_mv,
                                                       &x_bottom_mv))) {
                    total_pan_high_amp_sbs++;
                }

                if ((EbBool)(check_mv_for_tilt_high_amp(pcs_ptr->hierarchical_levels,
                                                        pcs_ptr->temporal_layer_index,
                                                        &y_current_mv,
                                                        &y_left_mv) ||
                             check_mv_for_tilt_high_amp(pcs_ptr->hierarchical_levels,
                                                        pcs_ptr->temporal_layer_index,
                                                        &y_current_mv,
                                                        &y_top_mv) ||
                             check_mv_for_tilt_high_amp(pcs_ptr->hierarchical_levels,
                                                        pcs_ptr->temporal_layer_index,
                                                        &y_current_mv,
                                                        &y_right_mv) ||
                             check_mv_for_tilt_high_amp(pcs_ptr->hierarchical_levels,
                                                        pcs_ptr->temporal_layer_index,
                                                        &y_current_mv,
                                                        &y_bottom_mv))) {
                    total_tilt_high_amp_sbs++;
                }
            }
        }
        pcs_ptr->is_pan  = EB_FALSE;
        pcs_ptr->is_tilt = EB_FALSE;

        pcs_ptr->pan_mvx  = 0;
        pcs_ptr->pan_mvy  = 0;
        pcs_ptr->tilt_mvx = 0;
        pcs_ptr->tilt_mvy = 0;

        // If more than PAN_SB_PERCENTAGE % of SBs are PAN
        if ((total_checked_sbs > 0) &&
            ((total_pan_sbs * 100 / total_checked_sbs) > PAN_SB_PERCENTAGE)) {
            pcs_ptr->is_pan = EB_TRUE;
            if (total_pan_sbs > 0) {
                pcs_ptr->pan_mvx = (int16_t)(x_pan_mv_sum / total_pan_sbs);
                pcs_ptr->pan_mvy = (int16_t)(y_pan_mv_sum / total_pan_sbs);
            }
        }

        if ((total_checked_sbs > 0) &&
            ((total_tilt_sbs * 100 / total_checked_sbs) > PAN_SB_PERCENTAGE)) {
            pcs_ptr->is_tilt = EB_TRUE;
            if (total_tilt_sbs > 0) {
                pcs_ptr->tilt_mvx = (int16_t)(x_tile_mv_sum / total_tilt_sbs);
                pcs_ptr->tilt_mvy = (int16_t)(y_tilt_mv_sum / total_tilt_sbs);
            }
        }
    }
}

static void initial_rate_control_context_dctor(EbPtr p) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)p;
    InitialRateControlContext *obj = (InitialRateControlContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}

/************************************************
* Initial Rate Control Context Constructor
************************************************/
EbErrorType initial_rate_control_context_ctor(EbThreadContext *  thread_context_ptr,
                                              const EbEncHandle *enc_handle_ptr) {
    InitialRateControlContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = initial_rate_control_context_dctor;

    context_ptr->motion_estimation_results_input_fifo_ptr = eb_system_resource_get_consumer_fifo(
        enc_handle_ptr->motion_estimation_results_resource_ptr, 0);
    context_ptr->initialrate_control_results_output_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->initial_rate_control_results_resource_ptr, 0);

    return EB_ErrorNone;
}

/************************************************
* Release Pa Reference Objects
** Check if reference pictures are needed
** release them when appropriate
************************************************/
void release_pa_reference_objects(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {
    // PA Reference Pictures
    uint32_t num_of_list_to_search;
    uint32_t list_index;
    uint32_t ref_pic_index;
    if (pcs_ptr->slice_type != I_SLICE) {
        num_of_list_to_search = (pcs_ptr->slice_type == P_SLICE) ? REF_LIST_0 : REF_LIST_1;

        // List Loop
        for (list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
            // Release PA Reference Pictures
            uint8_t num_of_ref_pic_to_search =
                (pcs_ptr->slice_type == P_SLICE)
                    ? MIN(pcs_ptr->ref_list0_count, scs_ptr->reference_count)
                    : (list_index == REF_LIST_0)
                          ? MIN(pcs_ptr->ref_list0_count, scs_ptr->reference_count)
                          : MIN(pcs_ptr->ref_list1_count, scs_ptr->reference_count);

            for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
                if (pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index] != NULL) {
                    eb_release_object(pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]);
                }
            }
        }
    }

    if (pcs_ptr->pa_reference_picture_wrapper_ptr != NULL) {
        eb_release_object(pcs_ptr->pa_reference_picture_wrapper_ptr);
    }

    return;
}

/************************************************
* Global Motion Detection Based on ME information
** Mark pictures for pan
** Mark pictures for tilt
** No lookahead information used in this function
************************************************/
void me_based_global_motion_detection(PictureParentControlSet *pcs_ptr) {
    // PAN Generation
    pcs_ptr->is_pan  = EB_FALSE;
    pcs_ptr->is_tilt = EB_FALSE;

    if (pcs_ptr->slice_type != I_SLICE) detect_global_motion(pcs_ptr);
    // Check if the motion vector field for temporal layer 0 pictures
    if (pcs_ptr->slice_type != I_SLICE && pcs_ptr->temporal_layer_index == 0)
        check_for_non_uniform_motion_vector_field(pcs_ptr);
    return;
}
/************************************************
* Global Motion Detection Based on Lookahead
** Mark pictures for pan
** Mark pictures for tilt
** LAD Window: min (8 or sliding window size)
************************************************/
void update_global_motion_detection_over_time(EncodeContext *          encode_context_ptr,
                                              SequenceControlSet *     scs_ptr,
                                              PictureParentControlSet *pcs_ptr) {
    InitialRateControlReorderEntry *temp_queue_entry_ptr;
    PictureParentControlSet *       temp_pcs_ptr;

    uint32_t total_pan_pictures     = 0;
    uint32_t total_checked_pictures = 0;
    uint32_t total_tilt_pictures    = 0;
    uint32_t update_is_pan_frames_to_check;
    uint32_t input_queue_index;
    uint32_t frames_to_check_index;

    (void)scs_ptr;

    // Determine number of frames to check (8 frames)
    update_is_pan_frames_to_check = MIN(8, pcs_ptr->frames_in_sw);

    // Walk the first N entries in the sliding window
    input_queue_index = encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    uint32_t update_frames_to_check = update_is_pan_frames_to_check;
    for (frames_to_check_index = 0; frames_to_check_index < update_frames_to_check;
         frames_to_check_index++) {
        temp_queue_entry_ptr =
            encode_context_ptr->initial_rate_control_reorder_queue[input_queue_index];
        temp_pcs_ptr =
            ((PictureParentControlSet *)(temp_queue_entry_ptr->parent_pcs_wrapper_ptr)->object_ptr);

        if (temp_pcs_ptr->slice_type != I_SLICE) {
            total_pan_pictures += (temp_pcs_ptr->is_pan == EB_TRUE);

            total_tilt_pictures += (temp_pcs_ptr->is_tilt == EB_TRUE);

            // Keep track of checked pictures
            total_checked_pictures++;
        }

        // Increment the input_queue_index Iterator
        input_queue_index = (input_queue_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                                ? 0
                                : input_queue_index + 1;
    }

    pcs_ptr->is_pan  = EB_FALSE;
    pcs_ptr->is_tilt = EB_FALSE;

    if (total_checked_pictures) {
        if (pcs_ptr->slice_type != I_SLICE) {
            if ((total_pan_pictures * 100 / total_checked_pictures) > 75) pcs_ptr->is_pan = EB_TRUE;
        }
    }
    return;
}

/************************************************
* Update BEA Information Based on Lookahead
** Average zzCost of Collocated SB throughout lookahead frames
** Set isMostOfPictureNonMoving based on number of non moving SBs
** LAD Window: min (2xmgpos+1 or sliding window size)
************************************************/

void update_bea_info_over_time(EncodeContext *          encode_context_ptr,
                               PictureParentControlSet *pcs_ptr) {
    InitialRateControlReorderEntry *temp_queue_entry_ptr;
    PictureParentControlSet *       temp_pcs_ptr;
    uint32_t                        update_non_moving_index_array_frames_to_check;
    uint16_t                        sb_idx;
    uint16_t                        frames_to_check_index;
    uint64_t                        non_moving_index_sum = 0;
    uint32_t                        input_queue_index;

    SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
    // Update motionIndexArray of the current picture by averaging the motionIndexArray of the N future pictures
    // Determine number of frames to check N
    update_non_moving_index_array_frames_to_check =
        MIN(MIN(((pcs_ptr->pred_struct_ptr->pred_struct_period << 1) + 1), pcs_ptr->frames_in_sw),
            scs_ptr->static_config.look_ahead_distance);
    uint64_t me_dist           = 0;
    uint8_t  me_dist_pic_count = 0;
    // SB Loop
    for (sb_idx = 0; sb_idx < pcs_ptr->sb_total_count; ++sb_idx) {
        uint16_t non_moving_index_over_sliding_window = pcs_ptr->non_moving_index_array[sb_idx];

        // Walk the first N entries in the sliding window starting picture + 1
        input_queue_index =
            (encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
             INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                ? 0
                : encode_context_ptr->initial_rate_control_reorder_queue_head_index + 1;
        for (frames_to_check_index = 0;
             frames_to_check_index < update_non_moving_index_array_frames_to_check - 1;
             frames_to_check_index++) {
            temp_queue_entry_ptr =
                encode_context_ptr->initial_rate_control_reorder_queue[input_queue_index];
            temp_pcs_ptr =
                ((PictureParentControlSet *)(temp_queue_entry_ptr->parent_pcs_wrapper_ptr)
                     ->object_ptr);

            if (temp_pcs_ptr->slice_type == I_SLICE || temp_pcs_ptr->end_of_sequence_flag) break;
            // Limit the distortion to lower layers 0, 1 and 2 only. Higher layers have close temporal distance and lower distortion that might contaminate the data
            if (temp_pcs_ptr->temporal_layer_index <
                MAX((int8_t)pcs_ptr->hierarchical_levels - 1, 2)) {
                if (sb_idx == 0) me_dist_pic_count++;
                me_dist += (temp_pcs_ptr->slice_type == I_SLICE)
                               ? 0
                               : (uint64_t)temp_pcs_ptr->rc_me_distortion[sb_idx];
            }
            // Store the filtered_sse of next ALT_REF picture in the I slice to be used in QP Scaling
            if (pcs_ptr->slice_type == I_SLICE && pcs_ptr->filtered_sse == 0 && sb_idx == 0 &&
                temp_pcs_ptr->temporal_layer_index == 0) {
                pcs_ptr->filtered_sse    = temp_pcs_ptr->filtered_sse;
                pcs_ptr->filtered_sse_uv = temp_pcs_ptr->filtered_sse_uv;
            }
            non_moving_index_over_sliding_window += temp_pcs_ptr->non_moving_index_array[sb_idx];

            // Increment the input_queue_index Iterator
            input_queue_index =
                (input_queue_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : input_queue_index + 1;
        }
        pcs_ptr->non_moving_index_array[sb_idx] =
            (uint8_t)(non_moving_index_over_sliding_window / (frames_to_check_index + 1));

        non_moving_index_sum += pcs_ptr->non_moving_index_array[sb_idx];
    }

    pcs_ptr->non_moving_index_average = (uint16_t)non_moving_index_sum / pcs_ptr->sb_total_count;
    me_dist_pic_count                 = MAX(me_dist_pic_count, 1);
    pcs_ptr->qp_scaling_average_complexity =
        (uint16_t)((uint64_t)me_dist / pcs_ptr->sb_total_count / 256 / me_dist_pic_count);
    return;
}

/****************************************
* Init ZZ Cost array to default values
** Used when no Lookahead is available
****************************************/
void init_zz_cost_info(PictureParentControlSet *pcs_ptr) {
    uint16_t sb_idx;
    pcs_ptr->non_moving_index_average = INVALID_ZZ_COST;

    // SB Loop
    for (sb_idx = 0; sb_idx < pcs_ptr->sb_total_count; ++sb_idx)
        pcs_ptr->non_moving_index_array[sb_idx] = INVALID_ZZ_COST;
    return;
}

/************************************************
* Update uniform motion field
** Update Uniformly moving SBs using
** collocated SBs infor in lookahead pictures
** LAD Window: min (2xmgpos+1 or sliding window size)
************************************************/
void update_motion_field_uniformity_over_time(EncodeContext *          encode_context_ptr,
                                              SequenceControlSet *     scs_ptr,
                                              PictureParentControlSet *pcs_ptr) {
    InitialRateControlReorderEntry *temp_queue_entry_ptr;
    PictureParentControlSet *       temp_pcs_ptr;
    uint32_t                        input_queue_index;
    uint32_t                        no_frames_to_check;
    uint32_t                        frames_to_check_index;
    //SVT_LOG("To update POC %d\tframesInSw = %d\n", pcs_ptr->picture_number, pcs_ptr->frames_in_sw);
    // Determine number of frames to check N
    no_frames_to_check =
        MIN(MIN(((pcs_ptr->pred_struct_ptr->pred_struct_period << 1) + 1), pcs_ptr->frames_in_sw),
            scs_ptr->static_config.look_ahead_distance);

    // Walk the first N entries in the sliding window starting picture + 1
    input_queue_index = (encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
                         INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? 0
                            : encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    for (frames_to_check_index = 0; frames_to_check_index < no_frames_to_check - 1;
         frames_to_check_index++) {
        temp_queue_entry_ptr =
            encode_context_ptr->initial_rate_control_reorder_queue[input_queue_index];
        temp_pcs_ptr =
            ((PictureParentControlSet *)(temp_queue_entry_ptr->parent_pcs_wrapper_ptr)->object_ptr);

        if (temp_pcs_ptr->end_of_sequence_flag) break;
        // Increment the input_queue_index Iterator
        input_queue_index = (input_queue_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                                ? 0
                                : input_queue_index + 1;
    }
    return;
}
InitialRateControlReorderEntry *determine_picture_offset_in_queue(
    EncodeContext *encode_context_ptr, PictureParentControlSet *pcs_ptr,
    MotionEstimationResults *in_results_ptr) {
    InitialRateControlReorderEntry *queue_entry_ptr;
    int32_t                         queue_entry_index;

    queue_entry_index =
        (int32_t)(pcs_ptr->picture_number -
                  encode_context_ptr
                      ->initial_rate_control_reorder_queue
                          [encode_context_ptr->initial_rate_control_reorder_queue_head_index]
                      ->picture_number);
    queue_entry_index += encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    queue_entry_index = (queue_entry_index > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? queue_entry_index - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                            : queue_entry_index;
    queue_entry_ptr = encode_context_ptr->initial_rate_control_reorder_queue[queue_entry_index];
    queue_entry_ptr->parent_pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;
    queue_entry_ptr->picture_number         = pcs_ptr->picture_number;

    return queue_entry_ptr;
}

void get_histogram_queue_data(SequenceControlSet *scs_ptr, EncodeContext *encode_context_ptr,
                              PictureParentControlSet *pcs_ptr) {
    HlRateControlHistogramEntry *histogram_queue_entry_ptr;
    int32_t                      histogram_queue_entry_index;

    // Determine offset from the Head Ptr for HLRC histogram queue
    eb_block_on_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
    histogram_queue_entry_index = (int32_t)(
        pcs_ptr->picture_number -
        encode_context_ptr
            ->hl_rate_control_historgram_queue[encode_context_ptr
                                                   ->hl_rate_control_historgram_queue_head_index]
            ->picture_number);
    histogram_queue_entry_index += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
    histogram_queue_entry_index =
        (histogram_queue_entry_index > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1)
            ? histogram_queue_entry_index - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
            : (histogram_queue_entry_index < 0)
                  ? histogram_queue_entry_index + HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
                  : histogram_queue_entry_index;
    histogram_queue_entry_ptr =
        encode_context_ptr->hl_rate_control_historgram_queue[histogram_queue_entry_index];

    //histogram_queue_entry_ptr->parent_pcs_wrapper_ptr  = in_results_ptr->pcs_wrapper_ptr;
    histogram_queue_entry_ptr->picture_number       = pcs_ptr->picture_number;
    histogram_queue_entry_ptr->end_of_sequence_flag = pcs_ptr->end_of_sequence_flag;
    histogram_queue_entry_ptr->slice_type           = pcs_ptr->slice_type;
    histogram_queue_entry_ptr->temporal_layer_index = pcs_ptr->temporal_layer_index;
    histogram_queue_entry_ptr->full_sb_count        = pcs_ptr->full_sb_count;
    histogram_queue_entry_ptr->life_count           = 0;
    histogram_queue_entry_ptr->passed_to_hlrc       = EB_FALSE;
    histogram_queue_entry_ptr->is_coded             = EB_FALSE;
    histogram_queue_entry_ptr->total_num_bits_coded = 0;
    histogram_queue_entry_ptr->frames_in_sw         = 0;
    eb_memcpy(histogram_queue_entry_ptr->me_distortion_histogram,
              pcs_ptr->me_distortion_histogram,
              sizeof(uint16_t) * NUMBER_OF_SAD_INTERVALS);

    eb_memcpy(histogram_queue_entry_ptr->ois_distortion_histogram,
              pcs_ptr->ois_distortion_histogram,
              sizeof(uint16_t) * NUMBER_OF_INTRA_SAD_INTERVALS);

    eb_release_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
    //SVT_LOG("Test1 POC: %d\t POC: %d\t LifeCount: %d\n", histogram_queue_entry_ptr->picture_number, pcs_ptr->picture_number,  histogram_queue_entry_ptr->life_count);

    return;
}

void update_histogram_queue_entry(SequenceControlSet *scs_ptr, EncodeContext *encode_context_ptr,
                                  PictureParentControlSet *pcs_ptr, uint32_t frames_in_sw) {
    HlRateControlHistogramEntry *histogram_queue_entry_ptr;
    int32_t                      histogram_queue_entry_index;

    eb_block_on_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);

    histogram_queue_entry_index = (int32_t)(
        pcs_ptr->picture_number -
        encode_context_ptr
            ->hl_rate_control_historgram_queue[encode_context_ptr
                                                   ->hl_rate_control_historgram_queue_head_index]
            ->picture_number);
    histogram_queue_entry_index += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
    histogram_queue_entry_index =
        (histogram_queue_entry_index > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1)
            ? histogram_queue_entry_index - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
            : (histogram_queue_entry_index < 0)
                  ? histogram_queue_entry_index + HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
                  : histogram_queue_entry_index;
    histogram_queue_entry_ptr =
        encode_context_ptr->hl_rate_control_historgram_queue[histogram_queue_entry_index];
    histogram_queue_entry_ptr->passed_to_hlrc = EB_TRUE;

    if (scs_ptr->static_config.rate_control_mode == 2)
        histogram_queue_entry_ptr->life_count +=
            (int16_t)(scs_ptr->static_config.intra_period_length + 1) -
            3; // FramelevelRC does not decrease the life count for first picture in each temporal layer
    else
        histogram_queue_entry_ptr->life_count += pcs_ptr->historgram_life_count;

    histogram_queue_entry_ptr->frames_in_sw = frames_in_sw;
    eb_release_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);

    return;
}

/* Initial Rate Control Kernel */

/*********************************************************************************
*
* @brief
*  The Initial Rate Control process determines the initial bit budget for each picture
*  depending on the data gathered in the Picture Analysis and Motion Estimation processes
*  as well as the settings determined in the Picture Decision process.
*
* @par Description:
*  The Initial Rate Control process employs a sliding window buffer to analyze
*  multiple pictures if a delay is allowed. Note that no reference picture data is
*  used in this process.
*
* @param[in] Picture
*  The Initial Rate Control Kernel takes a picture and determines the initial bit budget
*  for each picture depending on the data that was gathered in Picture Analysis and
*  Motion Estimation processes
*
* @param[out] Bit Budget
*  Bit Budget is the amount of budgetted bits for a picture
*
* @remarks
*  Temporal noise reduction is currently performed in Initial Rate Control Process.
*  In the future we might decide to move it to Motion Analysis Process.
*
********************************************************************************/
void *initial_rate_control_kernel(void *input_ptr) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)input_ptr;
    InitialRateControlContext *context_ptr = (InitialRateControlContext *)thread_context_ptr->priv;
    PictureParentControlSet *  pcs_ptr;
    PictureParentControlSet *  pcs_ptr_temp;
    EncodeContext *            encode_context_ptr;
    SequenceControlSet *       scs_ptr;

    EbObjectWrapper *        in_results_wrapper_ptr;
    MotionEstimationResults *in_results_ptr;

    EbObjectWrapper *          out_results_wrapper_ptr;
    InitialRateControlResults *out_results_ptr;

    // Queue variables
    uint32_t                        queue_entry_index_temp;
    uint32_t                        queue_entry_index_temp2;
    InitialRateControlReorderEntry *queue_entry_ptr;

    EbBool           move_slide_window_flag = EB_TRUE;
    EbBool           end_of_sequence_flag   = EB_TRUE;
    uint8_t          frames_in_sw;
    uint8_t          temporal_layer_index;
    EbObjectWrapper *reference_picture_wrapper_ptr;

    // Segments
    uint32_t segment_index;
    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->motion_estimation_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);

        in_results_ptr = (MotionEstimationResults *)in_results_wrapper_ptr->object_ptr;
        pcs_ptr        = (PictureParentControlSet *)in_results_ptr->pcs_wrapper_ptr->object_ptr;

        segment_index = in_results_ptr->segment_index;

        // Set the segment mask
        SEGMENT_COMPLETION_MASK_SET(pcs_ptr->me_segments_completion_mask, segment_index);

        // If the picture is complete, proceed
        if (SEGMENT_COMPLETION_MASK_TEST(pcs_ptr->me_segments_completion_mask,
                                         pcs_ptr->me_segments_total_count)) {
            scs_ptr            = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
            encode_context_ptr = (EncodeContext *)scs_ptr->encode_context_ptr;
            // Mark picture when global motion is detected using ME results
            //reset intra_coded_estimation_sb
            me_based_global_motion_detection(pcs_ptr);
            // Release Pa Ref pictures when not needed
            release_pa_reference_objects(scs_ptr, pcs_ptr);

            //****************************************************
            // Input Motion Analysis Results into Reordering Queue
            //****************************************************

            if (!pcs_ptr->is_overlay)
                // Determine offset from the Head Ptr
                queue_entry_ptr =
                    determine_picture_offset_in_queue(encode_context_ptr, pcs_ptr, in_results_ptr);

            if (scs_ptr->static_config.rate_control_mode) {
                if (scs_ptr->static_config.look_ahead_distance != 0) {
                    // Getting the Histogram Queue Data
                    get_histogram_queue_data(scs_ptr, encode_context_ptr, pcs_ptr);
                }
            }

            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS;
                 temporal_layer_index++)
                pcs_ptr->frames_in_interval[temporal_layer_index] = 0;
            pcs_ptr->frames_in_sw          = 0;
            pcs_ptr->historgram_life_count = 0;
            pcs_ptr->scene_change_in_gop   = EB_FALSE;
            move_slide_window_flag = EB_TRUE;
            while (move_slide_window_flag) {
                // Check if the sliding window condition is valid
                queue_entry_index_temp =
                    encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                if (encode_context_ptr->initial_rate_control_reorder_queue[queue_entry_index_temp]
                        ->parent_pcs_wrapper_ptr != NULL)
                    end_of_sequence_flag =
                        (((PictureParentControlSet
                               *)(encode_context_ptr
                                      ->initial_rate_control_reorder_queue[queue_entry_index_temp]
                                      ->parent_pcs_wrapper_ptr)
                              ->object_ptr))
                            ->end_of_sequence_flag;
                else
                    end_of_sequence_flag = EB_FALSE;
                frames_in_sw = 0;
                while (move_slide_window_flag && !end_of_sequence_flag &&
                       queue_entry_index_temp <=
                           encode_context_ptr->initial_rate_control_reorder_queue_head_index +
                               scs_ptr->static_config.look_ahead_distance) {
                    // frames_in_sw <= scs_ptr->static_config.look_ahead_distance){
                    frames_in_sw++;

                    queue_entry_index_temp2 =
                        (queue_entry_index_temp > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? queue_entry_index_temp - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                            : queue_entry_index_temp;

                    move_slide_window_flag =
                        (EbBool)(move_slide_window_flag &&
                                 (encode_context_ptr
                                      ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                                      ->parent_pcs_wrapper_ptr != NULL));
                    if (encode_context_ptr
                            ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                            ->parent_pcs_wrapper_ptr != NULL) {
                        // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                        end_of_sequence_flag =
                            ((PictureParentControlSet *)(encode_context_ptr
                                                             ->initial_rate_control_reorder_queue
                                                                 [queue_entry_index_temp2]
                                                             ->parent_pcs_wrapper_ptr)
                                 ->object_ptr)
                                ->end_of_sequence_flag;
                    } else
                        end_of_sequence_flag = EB_FALSE;
                    queue_entry_index_temp++;
                }

                if (move_slide_window_flag) {
                    //get a new entry spot
                    queue_entry_ptr =
                        encode_context_ptr->initial_rate_control_reorder_queue
                            [encode_context_ptr->initial_rate_control_reorder_queue_head_index];
                    pcs_ptr = ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr)
                                   ->object_ptr);
                    scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
                    // overlay picture was not added to the queue. For the alt_ref picture with an overlay picture, it loops on both alt ref and overlay pictures
                    uint8_t has_overlay = pcs_ptr->is_alt_ref ? 1 : 0;
                    for (uint8_t loop_index = 0; loop_index <= has_overlay; loop_index++) {
                        if (loop_index) pcs_ptr = pcs_ptr->overlay_ppcs_ptr;
                        pcs_ptr->frames_in_sw = frames_in_sw;
                        queue_entry_index_temp =
                            encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                        end_of_sequence_flag = EB_FALSE;
                        // find the frames_in_interval for the peroid I frames
                        while (
                            !end_of_sequence_flag &&
                            queue_entry_index_temp <=
                                encode_context_ptr->initial_rate_control_reorder_queue_head_index +
                                    scs_ptr->static_config.look_ahead_distance) {
                            queue_entry_index_temp2 =
                                (queue_entry_index_temp >
                                 INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                                    ? queue_entry_index_temp -
                                          INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                                    : queue_entry_index_temp;
                            pcs_ptr_temp = ((PictureParentControlSet
                                                 *)(encode_context_ptr
                                                        ->initial_rate_control_reorder_queue
                                                            [queue_entry_index_temp2]
                                                        ->parent_pcs_wrapper_ptr)
                                                ->object_ptr);
                            if (scs_ptr->intra_period_length != -1) {
                                if (pcs_ptr->picture_number %
                                        ((scs_ptr->intra_period_length + 1)) ==
                                    0) {
                                    pcs_ptr
                                        ->frames_in_interval[pcs_ptr_temp->temporal_layer_index]++;
                                    if (pcs_ptr_temp->scene_change_flag)
                                        pcs_ptr->scene_change_in_gop = EB_TRUE;
                                }
                            }

                            pcs_ptr_temp->historgram_life_count++;
                            end_of_sequence_flag = pcs_ptr_temp->end_of_sequence_flag;
                            queue_entry_index_temp++;
                        }

                        if ((scs_ptr->static_config.look_ahead_distance != 0) &&
                            (frames_in_sw < (scs_ptr->static_config.look_ahead_distance + 1)))
                            pcs_ptr->end_of_sequence_region = EB_TRUE;
                        else
                            pcs_ptr->end_of_sequence_region = EB_FALSE;

                        if (scs_ptr->static_config.rate_control_mode) {
                            // Determine offset from the Head Ptr for HLRC histogram queue and set the life count
                            if (scs_ptr->static_config.look_ahead_distance != 0) {
                                // Update Histogram Queue Entry Life count
                                update_histogram_queue_entry(
                                    scs_ptr, encode_context_ptr, pcs_ptr, frames_in_sw);
                            }
                        }

                        // Mark each input picture as PAN or not
                        // If a lookahead is present then check PAN for a period of time
                        if (!pcs_ptr->end_of_sequence_flag &&
                            scs_ptr->static_config.look_ahead_distance != 0) {
                            // Check for Pan,Tilt, Zoom and other global motion detectors over the future pictures in the lookahead
                            update_global_motion_detection_over_time(
                                encode_context_ptr, scs_ptr, pcs_ptr);
                        } else {
                            if (pcs_ptr->slice_type != I_SLICE) detect_global_motion(pcs_ptr);
                        }

                        // BACKGROUND ENHANCEMENT Part II
                        if (!pcs_ptr->end_of_sequence_flag &&
                            scs_ptr->static_config.look_ahead_distance != 0) {
                            // Update BEA information based on Lookahead information
                            update_bea_info_over_time(encode_context_ptr, pcs_ptr);
                        } else {
                            // Reset zzCost information to default When there's no lookahead available
                            init_zz_cost_info(pcs_ptr);
                        }

                        // Use the temporal layer 0 is_sb_motion_field_non_uniform array for all the other layer pictures in the mini GOP
                        if (!pcs_ptr->end_of_sequence_flag &&
                            scs_ptr->static_config.look_ahead_distance != 0) {
                            // Updat uniformly moving SBs based on Collocated SBs in LookAhead window
                            update_motion_field_uniformity_over_time(
                                encode_context_ptr, scs_ptr, pcs_ptr);
                        }
                        // Get Empty Reference Picture Object
                        eb_get_empty_object(
                            scs_ptr->encode_context_ptr->reference_picture_pool_fifo_ptr,
                            &reference_picture_wrapper_ptr);
                        if (loop_index) {
                            pcs_ptr->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                            // Give the new Reference a nominal live_count of 1
                            eb_object_inc_live_count(pcs_ptr->reference_picture_wrapper_ptr, 1);
                        } else {
                            ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr
                                                             ->object_ptr))
                                ->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                            // Give the new Reference a nominal live_count of 1
                            eb_object_inc_live_count(
                                ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr
                                                                 ->object_ptr))
                                    ->reference_picture_wrapper_ptr,
                                1);
                        }
                        pcs_ptr->stat_struct_first_pass_ptr =
                            pcs_ptr->is_used_as_reference_flag
                                ? &((EbReferenceObject *)
                                        pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                                       ->stat_struct
                                : &pcs_ptr->stat_struct;
                        if (scs_ptr->use_output_stat_file)
                            memset(pcs_ptr->stat_struct_first_pass_ptr, 0, sizeof(StatStruct));
                        // Get Empty Results Object
                        eb_get_empty_object(
                            context_ptr->initialrate_control_results_output_fifo_ptr,
                            &out_results_wrapper_ptr);

                        out_results_ptr =
                            (InitialRateControlResults *)out_results_wrapper_ptr->object_ptr;

                        if (loop_index)
                            out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
                        else
                            out_results_ptr->pcs_wrapper_ptr =
                                queue_entry_ptr->parent_pcs_wrapper_ptr;
                        // Post the Full Results Object
                        eb_post_full_object(out_results_wrapper_ptr);
                    }
                    // Reset the Reorder Queue Entry
                    queue_entry_ptr->picture_number += INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
                    queue_entry_ptr->parent_pcs_wrapper_ptr = (EbObjectWrapper *)NULL;

                    // Increment the Reorder Queue head Ptr
                    encode_context_ptr->initial_rate_control_reorder_queue_head_index =
                        (encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
                         INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? 0
                            : encode_context_ptr->initial_rate_control_reorder_queue_head_index + 1;

                    queue_entry_ptr =
                        encode_context_ptr->initial_rate_control_reorder_queue
                            [encode_context_ptr->initial_rate_control_reorder_queue_head_index];
                }
            }
        }

        // Release the Input Results
        eb_release_object(in_results_wrapper_ptr);
    }
    return NULL;
}
