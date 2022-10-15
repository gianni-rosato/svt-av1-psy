// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* Copyright(c) 2019 Netflix, Inc.
*
* This source code is subject to the terms of the BSD 3-Clause Clear License and
* the Alliance for Open Media Patent License 1.0. If the BSD 3-Clause Clear License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "EbPictureDecisionProcess.h"
#include "EbDefinitions.h"
#include "EbEncHandle.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureAnalysisProcess.h"
#include "EbPictureAnalysisResults.h"
#include "EbPictureDecisionResults.h"
#include "EbReferenceObject.h"
#include "EbSvtAv1ErrorCodes.h"
#include "EbTemporalFiltering.h"
#include "EbObject.h"
#include "EbUtility.h"
#include "EbLog.h"
#include "common_dsp_rtcd.h"
#include "EbResize.h"
#include "EbMalloc.h"
#include "EbInterPrediction.h"
#include "aom_dsp_rtcd.h"

#include "EbPictureOperators.h"
/************************************************
 * Defines
 ************************************************/
#define  LAY0_OFF  0
#define  LAY1_OFF  3
#define  LAY2_OFF  5
#define  LAY3_OFF  6
#define  LAY4_OFF  7
#define  LAY5_OFF  2
extern PredictionStructureConfigEntry flat_pred_struct[];
extern PredictionStructureConfigEntry two_level_hierarchical_pred_struct[];
extern PredictionStructureConfigEntry three_level_hierarchical_pred_struct[];
extern PredictionStructureConfigEntry four_level_hierarchical_pred_struct[];
extern PredictionStructureConfigEntry five_level_hierarchical_pred_struct[];
extern PredictionStructureConfigEntry six_level_hierarchical_pred_struct[];
void  get_max_allocated_me_refs(uint8_t ref_count_used_list0, uint8_t ref_count_used_list1, uint8_t* max_ref_to_alloc, uint8_t* max_cand_to_alloc);
void init_resize_picture(SequenceControlSet* scs_ptr, PictureParentControlSet* pcs_ptr);
uint64_t  get_ref_poc(PictureDecisionContext *context, uint64_t curr_picture_number, int32_t delta_poc)
{
    uint64_t ref_poc;
    uint64_t boundary_poc = MAX(context->key_poc, context->sframe_poc);
    if ((int64_t)curr_picture_number - (int64_t)delta_poc < (int64_t)boundary_poc) {
        ref_poc = boundary_poc;
    }
    else {
        ref_poc = (int64_t)curr_picture_number - (int64_t)delta_poc;
    }

    return ref_poc;
}

MvReferenceFrame svt_get_ref_frame_type(uint8_t list, uint8_t ref_idx);

static INLINE int get_relative_dist(const OrderHintInfo *oh, int a, int b) {
    if (!oh->enable_order_hint) return 0;

    const int bits = oh->order_hint_bits;

    assert(bits >= 1);
    assert(a >= 0 && a < (1 << bits));
    assert(b >= 0 && b < (1 << bits));

    int diff = a - b;
    const int m = 1 << (bits - 1);
    diff = (diff & (m - 1)) - (diff & m);
    return diff;
}

void svt_av1_setup_skip_mode_allowed(PictureParentControlSet  *parent_pcs_ptr) {

    FrameHeader *frm_hdr = &parent_pcs_ptr->frm_hdr;

    struct RefFrameInfo {
        int used;
        uint64_t poc;
    } ref_frame_arr_single[7];

    for (uint8_t i = 0; i < 7; ++i)
        ref_frame_arr_single[i].used = 1;

    for (uint8_t i = 0; i < 7; ++i) {
        ref_frame_arr_single[i].poc = parent_pcs_ptr->av1_ref_signal.ref_poc_array[i] % (uint64_t)(1 << (parent_pcs_ptr->scs_ptr->seq_header.order_hint_info.order_hint_bits));
    }

    OrderHintInfo order_hint_info_st;
    order_hint_info_st.enable_order_hint = 1;
    order_hint_info_st.order_hint_bits = 6+1;

    const OrderHintInfo *const order_hint_info = &order_hint_info_st;// cm->seq_params.order_hint_info;
    SkipModeInfo *const skip_mode_info = &frm_hdr->skip_mode_params;// cm->current_frame.skip_mode_info;

    skip_mode_info->skip_mode_allowed = 0;
    skip_mode_info->ref_frame_idx_0 = INVALID_IDX;
    skip_mode_info->ref_frame_idx_1 = INVALID_IDX;
    parent_pcs_ptr->cur_order_hint = parent_pcs_ptr->picture_number % (uint64_t)(1 << (parent_pcs_ptr->scs_ptr->seq_header.order_hint_info.order_hint_bits));
    for (uint8_t i = 0; i < 7; ++i)
        parent_pcs_ptr->ref_order_hint[i] = (uint32_t)ref_frame_arr_single[i].poc;
    if (/*!order_hint_info->enable_order_hint ||*/ parent_pcs_ptr->slice_type == I_SLICE /*frame_is_intra_only(cm)*/ ||
        frm_hdr->reference_mode == SINGLE_REFERENCE)
        return;

    const int cur_order_hint = parent_pcs_ptr->picture_number % (uint64_t)(1 << (parent_pcs_ptr->scs_ptr->seq_header.order_hint_info.order_hint_bits));
    int ref_order_hints[2] = { -1, INT_MAX };
    int ref_idx[2] = { INVALID_IDX, INVALID_IDX };

    // Identify the nearest forward and backward references.

    for (int i = 0; i < INTER_REFS_PER_FRAME; ++i) {
        //const RefCntBuffer *const buf = get_ref_frame_buf(cm, LAST_FRAME + i);
        //if (buf == NULL) continue;

        if (ref_frame_arr_single[i].used == 0) continue;

        const int ref_order_hint = (const int)ref_frame_arr_single[i].poc;// buf->order_hint;
        if (get_relative_dist(order_hint_info, ref_order_hint, cur_order_hint) <
            0) {
            // Forward reference
            if (ref_order_hints[0] == -1 ||
                get_relative_dist(order_hint_info, ref_order_hint,
                    ref_order_hints[0]) > 0) {
                ref_order_hints[0] = ref_order_hint;
                ref_idx[0] = i;
            }
        }
        else if (get_relative_dist(order_hint_info, ref_order_hint,
            cur_order_hint) > 0) {
            // Backward reference
            if (ref_order_hints[1] == INT_MAX ||
                get_relative_dist(order_hint_info, ref_order_hint,
                    ref_order_hints[1]) < 0) {
                ref_order_hints[1] = ref_order_hint;
                ref_idx[1] = i;
            }
        }
    }

    if (ref_idx[0] != INVALID_IDX && ref_idx[1] != INVALID_IDX) {
        // == Bi-directional prediction ==
        skip_mode_info->skip_mode_allowed = 1;
        skip_mode_info->ref_frame_idx_0 = AOMMIN(ref_idx[0], ref_idx[1]);
        skip_mode_info->ref_frame_idx_1 = AOMMAX(ref_idx[0], ref_idx[1]);
    }
    else if (ref_idx[0] != INVALID_IDX && ref_idx[1] == INVALID_IDX) {
        // == Forward prediction only ==
        // Identify the second nearest forward reference.
        ref_order_hints[1] = -1;
        for (int i = 0; i < INTER_REFS_PER_FRAME; ++i) {
            //const RefCntBuffer *const buf = get_ref_frame_buf(cm, LAST_FRAME + i);
            //if (buf == NULL) continue;
            if (ref_frame_arr_single[i].used == 0) continue;

            const int ref_order_hint = (const int)ref_frame_arr_single[i].poc;// buf->order_hint;
            if ((ref_order_hints[0] != -1 &&
                get_relative_dist(order_hint_info, ref_order_hint,
                    ref_order_hints[0]) < 0) &&
                    (ref_order_hints[1] == -1 ||
                        get_relative_dist(order_hint_info, ref_order_hint,
                            ref_order_hints[1]) > 0)) {
                // Second closest forward reference
                ref_order_hints[1] = ref_order_hint;
                ref_idx[1] = i;
            }
        }
        if (ref_order_hints[1] != -1) {
            skip_mode_info->skip_mode_allowed = 1;
            skip_mode_info->ref_frame_idx_0 = AOMMIN(ref_idx[0], ref_idx[1]);
            skip_mode_info->ref_frame_idx_1 = AOMMAX(ref_idx[0], ref_idx[1]);
        }
    }

    //output: idx
    //0 :LAST
    //1 :LAST2
    //2 :LAST3
    //3 :GOLD
    //4 :BWD
    //5 :ALT2
    //6 :ALT
}
uint8_t  circ_inc(uint8_t max, uint8_t off, uint8_t input)
{
    input++;
    if (input >= max)
        input = 0;

    if (off == 2)
    {
        input++;
        if (input >= max)
            input = 0;
    }

    return input;
}
#define POC_CIRCULAR_ADD(base, offset/*, bits*/)             (/*(((int32_t) (base)) + ((int32_t) (offset)) > ((int32_t) (1 << (bits))))   ? ((base) + (offset) - (1 << (bits))) : \
                                                             (((int32_t) (base)) + ((int32_t) (offset)) < 0)                           ? ((base) + (offset) + (1 << (bits))) : \
                                                                                                                                       */((base) + (offset)))
EbErrorType derive_tf_window_params(
    SequenceControlSet *scs_ptr,
    EncodeContext *encode_context_ptr,
    PictureParentControlSet *pcs_ptr,
    PictureDecisionContext *context_ptr,
    uint32_t out_stride_diff64);
#define FLASH_TH                            5
#define FADE_TH                             3
#define SCENE_TH                            3000
#define NOISY_SCENE_TH                      4500    // SCD TH in presence of noise
#define HIGH_PICTURE_VARIANCE_TH            1500
#define NUM64x64INPIC(w,h)          ((w*h)>> (svt_log2f(BLOCK_SIZE_64)<<1))
#define QUEUE_GET_PREVIOUS_SPOT(h)  ((h == 0) ? PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1 : h - 1)
#define QUEUE_GET_NEXT_SPOT(h,off)  (( (h+off) >= PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH) ? h+off - PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH  : h + off)

#define WTH 64
#define OTH 64
static void picture_decision_context_dctor(EbPtr p)
{
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    PictureDecisionContext* obj = (PictureDecisionContext*)thread_context_ptr->priv;

    if (obj->prev_picture_histogram) {
        for (int region_in_picture_width_index = 0; region_in_picture_width_index < MAX_NUMBER_OF_REGIONS_IN_WIDTH; region_in_picture_width_index++) {
            if (obj->prev_picture_histogram[region_in_picture_width_index]) {
                for (int region_in_picture_height_index = 0; region_in_picture_height_index < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; region_in_picture_height_index++) {
                    EB_FREE_ARRAY(obj->prev_picture_histogram[region_in_picture_width_index][region_in_picture_height_index]);
                }
            }
            EB_FREE_PTR_ARRAY(obj->prev_picture_histogram[region_in_picture_width_index], MAX_NUMBER_OF_REGIONS_IN_HEIGHT);
        }
        EB_FREE_PTR_ARRAY(obj->prev_picture_histogram, MAX_NUMBER_OF_REGIONS_IN_WIDTH);
    }
    EB_FREE_2D(obj->ahd_running_avg);
    EB_FREE_2D(obj->ahd_running_avg_cr);
    EB_FREE_2D(obj->ahd_running_avg_cb);
    EB_FREE_ARRAY(obj);
}

 /************************************************
  * Picture Analysis Context Constructor
  ************************************************/
EbErrorType picture_decision_context_ctor(
    EbThreadContext     *thread_context_ptr,
    const EbEncHandle   *enc_handle_ptr,
    uint8_t scene_change_detection)
{
    PictureDecisionContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv = context_ptr;
    thread_context_ptr->dctor = picture_decision_context_dctor;

     memset(context_ptr->tf_pic_array, 0, (1 << MAX_TEMPORAL_LAYERS) * sizeof(PictureParentControlSet *));
     context_ptr->tf_pic_arr_cnt = 0;
    context_ptr->picture_analysis_results_input_fifo_ptr =
        svt_system_resource_get_consumer_fifo(enc_handle_ptr->picture_analysis_results_resource_ptr, 0);
    context_ptr->picture_decision_results_output_fifo_ptr =
        svt_system_resource_get_producer_fifo(enc_handle_ptr->picture_decision_results_resource_ptr, 0);

    if (scene_change_detection) {
        EB_ALLOC_PTR_ARRAY(context_ptr->prev_picture_histogram, MAX_NUMBER_OF_REGIONS_IN_WIDTH);
        for (uint32_t region_in_picture_width_index = 0; region_in_picture_width_index < MAX_NUMBER_OF_REGIONS_IN_WIDTH; region_in_picture_width_index++) { // loop over horizontal regions
            EB_ALLOC_PTR_ARRAY(context_ptr->prev_picture_histogram[region_in_picture_width_index], MAX_NUMBER_OF_REGIONS_IN_HEIGHT);
            for (uint32_t region_in_picture_height_index = 0; region_in_picture_height_index < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; region_in_picture_height_index++) {
                EB_CALLOC_ARRAY(context_ptr->prev_picture_histogram[region_in_picture_width_index][region_in_picture_height_index], HISTOGRAM_NUMBER_OF_BINS * sizeof(uint32_t));
            }
        }

        EB_CALLOC_2D(context_ptr->ahd_running_avg, MAX_NUMBER_OF_REGIONS_IN_WIDTH * sizeof(uint32_t), MAX_NUMBER_OF_REGIONS_IN_HEIGHT * sizeof(uint32_t));
    }
    context_ptr->reset_running_avg = TRUE;
    context_ptr->me_fifo_ptr = svt_system_resource_get_producer_fifo(
            enc_handle_ptr->me_pool_ptr_array[0], 0);


    context_ptr->mg_progress_id = 0;
    context_ptr->last_i_noise_levels_log1p_fp16[0] = 0;
    context_ptr->transition_detected = -1;
    context_ptr->sframe_poc = 0;
    context_ptr->sframe_due = 0;
    return EB_ErrorNone;
}
static Bool scene_transition_detector(
    PictureDecisionContext* context_ptr,
    SequenceControlSet* scs_ptr,
    PictureParentControlSet** parent_pcs_window)
{
    PictureParentControlSet* current_pcs_ptr = parent_pcs_window[1];
    PictureParentControlSet* future_pcs_ptr = parent_pcs_window[2];

    // calculating the frame threshold based on the number of 64x64 blocks in the frame
    uint32_t  region_threshhold;

    Bool is_abrupt_change; // this variable signals an abrubt change (scene change or flash)
    Bool is_scene_change; // this variable signals a frame representing a scene change

    uint32_t** ahd_running_avg = context_ptr->ahd_running_avg;

    uint32_t  region_in_picture_width_index;
    uint32_t  region_in_picture_height_index;

    uint32_t  region_width;
    uint32_t  region_height;
    uint32_t  region_width_offset;
    uint32_t  region_height_offset;

    uint32_t  is_abrupt_change_count = 0;
    uint32_t  is_scene_change_count = 0;

    uint32_t  region_count_threshold = (uint32_t)(((float)((scs_ptr->picture_analysis_number_of_regions_per_width * scs_ptr->picture_analysis_number_of_regions_per_height) * 50) / 100) + 0.5);

    region_width = parent_pcs_window[1]->enhanced_picture_ptr->width / scs_ptr->picture_analysis_number_of_regions_per_width;
    region_height = parent_pcs_window[1]->enhanced_picture_ptr->height / scs_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (region_in_picture_width_index = 0; region_in_picture_width_index < scs_ptr->picture_analysis_number_of_regions_per_width; region_in_picture_width_index++) {  // loop over horizontal regions
        for (region_in_picture_height_index = 0; region_in_picture_height_index < scs_ptr->picture_analysis_number_of_regions_per_height; region_in_picture_height_index++) { // loop over vertical regions

            is_abrupt_change = FALSE;
            is_scene_change = FALSE;

            // accumulative histogram (absolute) differences between the past and current frame
            uint32_t ahd = 0;

            region_width_offset = (region_in_picture_width_index == scs_ptr->picture_analysis_number_of_regions_per_width - 1) ?
                parent_pcs_window[1]->enhanced_picture_ptr->width - (scs_ptr->picture_analysis_number_of_regions_per_width * region_width) :
                0;

            region_height_offset = (region_in_picture_height_index == scs_ptr->picture_analysis_number_of_regions_per_height - 1) ?
                parent_pcs_window[1]->enhanced_picture_ptr->height - (scs_ptr->picture_analysis_number_of_regions_per_height * region_height) :
                0;

            region_width += region_width_offset;
            region_height += region_height_offset;

            region_threshhold = SCENE_TH * NUM64x64INPIC(region_width, region_height);

            for (int bin = 0; bin < HISTOGRAM_NUMBER_OF_BINS; ++bin) {
                ahd += ABS((int32_t)current_pcs_ptr->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][bin] - (int32_t)context_ptr->prev_picture_histogram[region_in_picture_width_index][region_in_picture_height_index][bin]);
            }

            if (context_ptr->reset_running_avg) {
                ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] = ahd;
            }

            uint32_t ahd_error = ABS(
                (int32_t)
                ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] -
                (int32_t)ahd);

            if (ahd_error > region_threshhold && ahd >= ahd_error) {
                is_abrupt_change = TRUE;
            }
            if (is_abrupt_change)
            {
                // this variable denotes the average intensity difference between the next and the past frames
                uint8_t aid_future_past = (uint8_t)ABS(
                    (int16_t)future_pcs_ptr
                    ->average_intensity_per_region[region_in_picture_width_index]
                    [region_in_picture_height_index] -
                    (int16_t)context_ptr
                    ->prev_average_intensity_per_region[region_in_picture_width_index]
                    [region_in_picture_height_index]);
                uint8_t   aid_future_present = (uint8_t)ABS((int16_t)future_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index] - (int16_t)current_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index]);
                uint8_t   aid_present_past = (uint8_t)ABS((int16_t)current_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index] - (int16_t)context_ptr->prev_average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index]);

                if (aid_future_past < FLASH_TH && aid_future_present >= FLASH_TH && aid_present_past >= FLASH_TH) {
                    //SVT_LOG ("\nFlash in frame# %i , %i\n", current_pcs_ptr->picture_number,aid_future_past);
                }
                else if (aid_future_present < FADE_TH && aid_present_past < FADE_TH) {
                    //SVT_LOG ("\nFlash in frame# %i , %i\n", current_pcs_ptr->picture_number,aid_future_past);
                }
                else {
                    is_scene_change = TRUE;
                    //SVT_LOG ("\nScene Change in frame# %i , %i\n", current_pcs_ptr->picture_number,aid_future_past);
                }
            }
            else
                ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] = (3 * ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] + ahd) / 4;
            is_abrupt_change_count += is_abrupt_change;
            is_scene_change_count += is_scene_change;
        }
    }

    context_ptr->reset_running_avg = is_abrupt_change_count >= region_count_threshold;
    return is_scene_change_count >= region_count_threshold;
}
/***************************************************************************************************
* release_prev_picture_from_reorder_queue
***************************************************************************************************/
EbErrorType release_prev_picture_from_reorder_queue(
    EncodeContext                 *encode_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    PictureDecisionReorderEntry   *queue_previous_entry_ptr;
    int32_t                           previous_entry_index;

    // Get the previous entry from the Picture Decision Reordering Queue (Entry N-1)
    // P.S. The previous entry in display order is needed for Scene Change Detection
    previous_entry_index = (encode_context_ptr->picture_decision_reorder_queue_head_index == 0) ? PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1 : encode_context_ptr->picture_decision_reorder_queue_head_index - 1;
    queue_previous_entry_ptr = encode_context_ptr->picture_decision_reorder_queue[previous_entry_index];

    // SB activity classification based on (0,0) SAD & picture activity derivation
    if (queue_previous_entry_ptr->parent_pcs_wrapper_ptr) {
        // Reset the Picture Decision Reordering Queue Entry
        // P.S. The reset of the Picture Decision Reordering Queue Entry could not be done before running the Scene Change Detector
        queue_previous_entry_ptr->picture_number += PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH;
        queue_previous_entry_ptr->parent_pcs_wrapper_ptr = (EbObjectWrapper *)NULL;
    }

    return return_error;
}

/***************************************************************************************************
* Initializes mini GOP activity array
*
***************************************************************************************************/
EbErrorType initialize_mini_gop_activity_array(
    PictureDecisionContext        *context_ptr) {
    // Loop over all mini GOPs
    for (uint32_t gopindex = 0; gopindex < MINI_GOP_MAX_COUNT; ++gopindex)
        context_ptr->mini_gop_activity_array[gopindex] = get_mini_gop_stats(gopindex)->hierarchical_levels != 3; /*MIN_HIERARCHICAL_LEVEL*/
    return EB_ErrorNone;
}

/***************************************************************************************************
* Generates block picture map
*
*
***************************************************************************************************/
static EbErrorType generate_picture_window_split(
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr) {
    context_ptr->total_number_of_mini_gops = 0;
    // Loop over all mini GOPs
    for (uint32_t gopindex = 0; gopindex < MINI_GOP_MAX_COUNT; gopindex += context_ptr->mini_gop_activity_array[gopindex]
        ? 1
        : mini_gop_offset[get_mini_gop_stats(gopindex)->hierarchical_levels - MIN_HIERARCHICAL_LEVEL]) {
        // Only for a valid mini GOP
        if (get_mini_gop_stats(gopindex)->end_index < encode_context_ptr->pre_assignment_buffer_count && !context_ptr->mini_gop_activity_array[gopindex]) {
            context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(gopindex)->start_index;
            context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(gopindex)->end_index;
            context_ptr->mini_gop_length[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(gopindex)->lenght;
            context_ptr->mini_gop_hierarchical_levels[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(gopindex)->hierarchical_levels;
            context_ptr->mini_gop_intra_count[context_ptr->total_number_of_mini_gops] = 0;
            context_ptr->mini_gop_idr_count[context_ptr->total_number_of_mini_gops] = 0;
            context_ptr->total_number_of_mini_gops++;
        }
    }
    // Only in presence of at least 1 valid mini GOP
    if (context_ptr->total_number_of_mini_gops != 0) {
        context_ptr->mini_gop_intra_count[context_ptr->total_number_of_mini_gops - 1] = encode_context_ptr->pre_assignment_buffer_intra_count;
        context_ptr->mini_gop_idr_count[context_ptr->total_number_of_mini_gops - 1] = encode_context_ptr->pre_assignment_buffer_idr_count;
    }
    return EB_ErrorNone;
}

/***************************************************************************************************
* Handles an incomplete picture window map
*
*
***************************************************************************************************/
EbErrorType handle_incomplete_picture_window_map(
    uint32_t                       hierarchical_level,
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    if (context_ptr->total_number_of_mini_gops == 0) {
        //The trailing frames of minigop16 will try to use minigop8 first
        hierarchical_level = hierarchical_level >=3 ? 3 : hierarchical_level;
        context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops] = 0;
        context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - 1;
        context_ptr->mini_gop_length[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops];
        context_ptr->mini_gop_hierarchical_levels[context_ptr->total_number_of_mini_gops] = hierarchical_level;

        context_ptr->total_number_of_mini_gops++;
    }
    else if (context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops - 1] < encode_context_ptr->pre_assignment_buffer_count - 1) {
        context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops] = context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops - 1] + 1;
        context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - 1;
        context_ptr->mini_gop_length[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops];
        context_ptr->mini_gop_hierarchical_levels[context_ptr->total_number_of_mini_gops] = 3;// MIN_HIERARCHICAL_LEVEL;// AMIR
        context_ptr->mini_gop_intra_count[context_ptr->total_number_of_mini_gops - 1] = 0;
        context_ptr->mini_gop_idr_count[context_ptr->total_number_of_mini_gops - 1] = 0;

        context_ptr->total_number_of_mini_gops++;
    }

    context_ptr->mini_gop_intra_count[context_ptr->total_number_of_mini_gops - 1] = encode_context_ptr->pre_assignment_buffer_intra_count;
    context_ptr->mini_gop_idr_count[context_ptr->total_number_of_mini_gops - 1] = encode_context_ptr->pre_assignment_buffer_idr_count;

    return return_error;
}
/*
   This function searches if the target input picture
   is still in the pre-assign buffer, if yes it returns
   its ppcs, else it returns Null
*/
PictureParentControlSet * is_pic_still_in_pre_assign_buffer(
    EncodeContext                 *encode_context_ptr,
    PictureDecisionContext        *context_ptr,
    uint32_t                       mini_gop_index,
    uint64_t                       target_pic)
{
    for (uint32_t pic_i = context_ptr->mini_gop_start_index[mini_gop_index]; pic_i <= context_ptr->mini_gop_end_index[mini_gop_index]; ++pic_i) {
        PictureParentControlSet*pcs_i = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pic_i]->object_ptr;
        if (pcs_i->picture_number == target_pic)
            return pcs_i;
    }
    return NULL;
}
/*
   This function tells if a picture is part of a short
   mg in RA configuration
*/
uint8_t is_pic_cutting_short_ra_mg(PictureDecisionContext   *context_ptr, PictureParentControlSet *pcs_ptr, uint32_t mg_idx)
{
    //if the size < complete MG or if there is usage of closed GOP
    if ((context_ptr->mini_gop_length[mg_idx] < pcs_ptr->pred_struct_ptr->pred_struct_period || context_ptr->mini_gop_idr_count[mg_idx] > 0) &&
        pcs_ptr->pred_struct_ptr->pred_type == SVT_AV1_PRED_RANDOM_ACCESS &&
        pcs_ptr->idr_flag == FALSE &&
        pcs_ptr->cra_flag == FALSE) {

        return 1;
    }
    else {
        return 0;
    }
}
/***************************************************************************************************
* If a switch happens, then update the RPS of the base layer frame separating the 2 different prediction structures
* Clean up the reference queue dependant counts of the base layer frame separating the 2 different prediction structures
*
***************************************************************************************************/
EbErrorType update_base_layer_reference_queue_dependent_count(
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr,
    SequenceControlSet            *scs_ptr,
    uint32_t                         gopindex) {
    if (!context_ptr || !encode_context_ptr || !scs_ptr)
        return EB_ErrorBadParameter;

    EbErrorType return_error = EB_ErrorNone;

    PictureParentControlSet       *pcs_ptr;

    // Get the 1st PCS mini GOP
    pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_start_index[gopindex]]->object_ptr;

    PictureParentControlSet  *trig_pcs; //triggering dep-cnt clean up picture, should be the first in dec order that goes to PicMgr
    if ( is_pic_cutting_short_ra_mg(context_ptr, pcs_ptr, gopindex)) {
        trig_pcs = pcs_ptr;//this an LDP minigop so the first picture in minigop is the first in dec order
    }
    else {
        //this an RA minigop so the last picture in minigop is the first in dec order
        uint32_t  last_pic_in_mg = context_ptr->mini_gop_end_index[gopindex];
        trig_pcs = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[last_pic_in_mg]->object_ptr;
    }
    // Derive the temporal layer difference between the current mini GOP and the previous mini GOP
    pcs_ptr->hierarchical_layers_diff = (int32_t)encode_context_ptr->previous_mini_gop_hierarchical_levels - (int32_t)pcs_ptr->hierarchical_levels;

    // Set init_pred_struct_position_flag to TRUE if mini GOP switch
    pcs_ptr->init_pred_struct_position_flag = encode_context_ptr->is_mini_gop_changed = (pcs_ptr->hierarchical_layers_diff != 0) ?
        TRUE :
        FALSE;

    // If the current mini GOP is different than the previous mini GOP update then update the positive dependant counts of the reference entry separating the 2 mini GOPs
    if (pcs_ptr->hierarchical_layers_diff != 0) {
        uint32_t input_queue_index = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

        while (input_queue_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
            PaReferenceQueueEntry *input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[input_queue_index];

            int32_t diff_n = 0;
            // Find the reference entry separating the 2 mini GOPs  (pcs_ptr->picture_number is the POC of the first isput in the mini GOP)
            if (input_entry_ptr->picture_number == (pcs_ptr->picture_number - 1)) {
                // Update the positive dependant counts

                // 1st step: remove all positive entries from the dependant list0 and dependant list1
                uint32_t dependant_list_positive_entries = 0;
                for (uint32_t dep_idx = 0; dep_idx < input_entry_ptr->list0.list_count; ++dep_idx) {
                    if (input_entry_ptr->list0.list[dep_idx] >= 0)
                        dependant_list_positive_entries++;
                }
                input_entry_ptr->list0.list_count = input_entry_ptr->list0.list_count - dependant_list_positive_entries;
                dependant_list_positive_entries = 0;
                for (uint32_t dep_idx = 0; dep_idx < input_entry_ptr->list1.list_count; ++dep_idx) {
                    if (input_entry_ptr->list1.list[dep_idx] >= 0)
                        dependant_list_positive_entries++;
                }
                input_entry_ptr->list1.list_count = input_entry_ptr->list1.list_count - dependant_list_positive_entries;

                // 2nd step: inherit the positive dependant counts of the current mini GOP
                // Get the RPS set of the current mini GOP
                PredictionStructure *next_pred_struct_ptr = get_prediction_structure(
                    encode_context_ptr->prediction_structure_group_ptr,
                    pcs_ptr->pred_structure,
                    scs_ptr->reference_count,
                    pcs_ptr->hierarchical_levels);            // Number of temporal layer in the current mini GOP

                // Get the RPS of a base layer input
                PredictionStructureEntry *next_base_layer_pred_position_ptr = next_pred_struct_ptr->pred_struct_entry_ptr_array[next_pred_struct_ptr->pred_struct_entry_count - 1];

                for (uint32_t dep_idx = 0; dep_idx < next_base_layer_pred_position_ptr->dep_list0.list_count; ++dep_idx) {
                    if (next_base_layer_pred_position_ptr->dep_list0.list[dep_idx] >= 0)
                        input_entry_ptr->list0.list[input_entry_ptr->list0.list_count++] = next_base_layer_pred_position_ptr->dep_list0.list[dep_idx];
                }

                for (uint32_t dep_idx = 0; dep_idx < next_base_layer_pred_position_ptr->dep_list1.list_count; ++dep_idx) {
                    if (next_base_layer_pred_position_ptr->dep_list1.list[dep_idx] >= 0)
                        input_entry_ptr->list1.list[input_entry_ptr->list1.list_count++] = next_base_layer_pred_position_ptr->dep_list1.list[dep_idx];
                }

                diff_n =
                    (int32_t)(input_entry_ptr->list0.list_count + input_entry_ptr->list1.list_count) - //depCnt after clean up
                    (int32_t)(input_entry_ptr->dep_list0_count + input_entry_ptr->dep_list1_count) + //depCnt from org prediction struct
                    (input_entry_ptr->is_alt_ref ? 1 : 0);

                // explanation for adding above line "((input_entry_ptr->is_alt_ref) ? 1 : 0)":
                // int32_t dep_list_count_org = input_entry_ptr->dep_list0_count + input_entry_ptr->dep_list1_count;
                // int32_t dep_list_count_new = input_entry_ptr->list0.list_count + input_entry_ptr->list1.list_count + ((input_entry_ptr->is_alt_ref) ? 1 : 0); // refer to below line 665
                // diff_n = dep_list_count_new - dep_list_count_org;
                // so, additional 1 from is_alt_ref is added to diff_n.

                //these refs are defintely not in the pre-assignment buffer
                if (diff_n) {
                    //TODO: change triggereing picture from the firt pic to the last pic in the MG for RA ( first in dec-order = base layer)
                    trig_pcs->updated_links_arr[trig_pcs->other_updated_links_cnt].pic_num = input_entry_ptr->picture_number;
                    trig_pcs->updated_links_arr[trig_pcs->other_updated_links_cnt++].dep_cnt_diff = diff_n;
                }
                // 3rd step: update the dependant count
                uint32_t dependant_list_removed_entries = input_entry_ptr->dep_list0_count + input_entry_ptr->dep_list1_count - input_entry_ptr->dependent_count;
                input_entry_ptr->dep_list0_count = (input_entry_ptr->is_alt_ref) ? input_entry_ptr->list0.list_count + 1 : input_entry_ptr->list0.list_count;
                input_entry_ptr->dep_list1_count = input_entry_ptr->list1.list_count;
                input_entry_ptr->dependent_count = input_entry_ptr->dep_list0_count + input_entry_ptr->dep_list1_count - dependant_list_removed_entries;
            } else {
                // Modify Dependent List0
                uint32_t dep_list_count = input_entry_ptr->list0.list_count;
                for (uint32_t dep_idx = 0; dep_idx < dep_list_count; ++dep_idx) {
                    // Adjust the latest currentInputPoc in case we're in a POC rollover scenario
                    // currentInputPoc += (currentInputPoc < input_entry_ptr->pocNumber) ? (1 << scs_ptr->bitsForPictureOrderCount) : 0;
                    int64_t dep_poc = POC_CIRCULAR_ADD(
                        (int64_t)input_entry_ptr->picture_number, // can't use a value that gets reset
                        input_entry_ptr->list0.list[dep_idx]/*,
                                                         scs_ptr->bitsForPictureOrderCount*/);

                                                         // If Dependent POC is greater or equal to the IDR POC
                    if (dep_poc >= (int64_t)pcs_ptr->picture_number && input_entry_ptr->list0.list[dep_idx]) {
                        input_entry_ptr->list0.list[dep_idx] = 0;

                        // Decrement the Reference's reference_count
                        --input_entry_ptr->dependent_count;
                        diff_n--;
                        CHECK_REPORT_ERROR(
                            (input_entry_ptr->dependent_count != ~0u),
                            encode_context_ptr->app_callback_ptr,
                            EB_ENC_PD_ERROR3);
                    }
                }
                // Modify Dependent List1
                dep_list_count = input_entry_ptr->list1.list_count;
                for (uint32_t dep_idx = 0; dep_idx < dep_list_count; ++dep_idx) {
                    // Adjust the latest currentInputPoc in case we're in a POC rollover scenario
                    // currentInputPoc += (currentInputPoc < input_entry_ptr->pocNumber) ? (1 << scs_ptr->bitsForPictureOrderCount) : 0;
                    int64_t dep_poc = POC_CIRCULAR_ADD(
                        (int64_t)input_entry_ptr->picture_number,
                        input_entry_ptr->list1.list[dep_idx]/*,
                                                         scs_ptr->bitsForPictureOrderCount*/);

                    // If Dependent POC is greater or equal to the IDR POC
                    if ((dep_poc >= (int64_t)pcs_ptr->picture_number) && input_entry_ptr->list1.list[dep_idx]) {
                        input_entry_ptr->list1.list[dep_idx] = 0;
                        // Decrement the Reference's reference_count
                        --input_entry_ptr->dependent_count;
                        diff_n--;
                        CHECK_REPORT_ERROR(
                            (input_entry_ptr->dependent_count != ~0u),
                            encode_context_ptr->app_callback_ptr,
                            EB_ENC_PD_ERROR3);
                    }
                }

                //these refs are defintely not in the pre-ass buffer
                if (diff_n) {
                    trig_pcs->updated_links_arr[trig_pcs->other_updated_links_cnt].pic_num = input_entry_ptr->picture_number;
                    trig_pcs->updated_links_arr[trig_pcs->other_updated_links_cnt++].dep_cnt_diff = diff_n;
                }
            }
            // Increment the input_queue_index Iterator
            input_queue_index = (input_queue_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : input_queue_index + 1;
        }
    }

    return return_error;
}

/***************************************************************************************************
* Generates mini GOP RPSs
*
*
***************************************************************************************************/
EbErrorType generate_mini_gop_rps(
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    SequenceControlSet            *scs_ptr;
    uint32_t                         mini_gop_index;
    PictureParentControlSet    *pcs_ptr;
    uint32_t                         out_stride_diff64;

    // Loop over all mini GOPs
    for (mini_gop_index = 0; mini_gop_index < context_ptr->total_number_of_mini_gops; ++mini_gop_index) {
        // Loop over picture within the mini GOP
        for (out_stride_diff64 = context_ptr->mini_gop_start_index[mini_gop_index]; out_stride_diff64 <= context_ptr->mini_gop_end_index[mini_gop_index]; out_stride_diff64++) {
            pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[out_stride_diff64]->object_ptr;
            scs_ptr = pcs_ptr->scs_ptr;
            pcs_ptr->pred_structure = scs_ptr->static_config.pred_structure;
            pcs_ptr->hierarchical_levels = (uint8_t)context_ptr->mini_gop_hierarchical_levels[mini_gop_index];

            pcs_ptr->pred_struct_ptr = get_prediction_structure(
                encode_context_ptr->prediction_structure_group_ptr,
                pcs_ptr->pred_structure,
                scs_ptr->reference_count,
                pcs_ptr->hierarchical_levels);
        }
    }
    return return_error;
}

static INLINE void update_list0_only_base(SequenceControlSet* scs, PictureParentControlSet* pcs) {

    // If noise_variance_th is MAX, then always skip list 1, else compare to the avg variance (if available)
    if (pcs->list0_only_base_ctrls.noise_variance_th == (uint16_t)~0 ||
        (scs->calculate_variance && (pcs->pic_avg_variance < pcs->list0_only_base_ctrls.noise_variance_th))) {
        pcs->ref_list1_count_try = 0;
    }
}

/* Wrapper function to compute TPL Synthesizer block size: Used in init memory allocation and TPL Controls*/
uint8_t svt_aom_get_tpl_synthesizer_block_size(int8_t tpl_level, uint32_t picture_width,
    uint32_t picture_height) {
    uint8_t blk_size;
    if (tpl_level <= 5)
        blk_size = AOMMIN(picture_width, picture_height) >= 720 ? 16 : 8;
    if (tpl_level <= 6)
        blk_size = 16;
    else
        blk_size = AOMMIN(picture_width, picture_height) >= 720 ? 32 : 16;

    return blk_size;
}

/*************************************************************************************
Set the TPL controls that control TPL search complexity.
***************************************************************************************/
static void set_tpl_extended_controls(PictureParentControlSet *pcs, uint8_t tpl_level) {
    TplControls *       tpl_ctrls = &pcs->tpl_ctrls;
    SequenceControlSet *scs = pcs->scs_ptr;
    const uint8_t            is_islice = pcs->slice_type == I_SLICE;
    const EbInputResolution resolution = scs->input_resolution;

    switch (tpl_level) {
    case 0:
        tpl_ctrls->enable = 0;
        tpl_ctrls->compute_rate = 0;
        tpl_ctrls->enable_tpl_qps = 0;
        tpl_ctrls->disable_intra_pred_nref = 0;
        tpl_ctrls->intra_mode_end = DC_PRED;
        tpl_ctrls->reduced_tpl_group = -1;
        tpl_ctrls->pf_shape = DEFAULT_SHAPE;
        tpl_ctrls->use_pred_sad_in_intra_search = 0;
        tpl_ctrls->use_pred_sad_in_inter_search = 0;
        tpl_ctrls->dispenser_search_level = 0;
        tpl_ctrls->subsample_tx = 0;
        tpl_ctrls->subpel_depth = FULL_PEL;
        break;
    case 1:
        tpl_ctrls->enable = 1;
        tpl_ctrls->compute_rate = 1;
        tpl_ctrls->enable_tpl_qps = 1;
        tpl_ctrls->disable_intra_pred_nref = 0;
        tpl_ctrls->intra_mode_end = PAETH_PRED;
        tpl_ctrls->reduced_tpl_group = -1;
        tpl_ctrls->pf_shape = DEFAULT_SHAPE;
        tpl_ctrls->use_pred_sad_in_intra_search = 0;
        tpl_ctrls->use_pred_sad_in_inter_search = 0;
        tpl_ctrls->dispenser_search_level = 0;
        tpl_ctrls->subsample_tx = 0;
        tpl_ctrls->subpel_depth = QUARTER_PEL;
        break;
    case 2:
        tpl_ctrls->enable = 1;
        tpl_ctrls->compute_rate = 0;
        tpl_ctrls->enable_tpl_qps = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end = DC_PRED;
        tpl_ctrls->reduced_tpl_group = -1;
        tpl_ctrls->pf_shape = DEFAULT_SHAPE;
        tpl_ctrls->use_pred_sad_in_intra_search = 0;
        tpl_ctrls->use_pred_sad_in_inter_search = 0;
        tpl_ctrls->dispenser_search_level = 0;
        tpl_ctrls->subsample_tx = 0;
        tpl_ctrls->subpel_depth = QUARTER_PEL;
        break;
    case 3:
        tpl_ctrls->enable = 1;
        tpl_ctrls->compute_rate = 0;
        tpl_ctrls->enable_tpl_qps = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end = DC_PRED;
        tpl_ctrls->reduced_tpl_group = -1;
        tpl_ctrls->pf_shape = DEFAULT_SHAPE;
        tpl_ctrls->use_pred_sad_in_intra_search = 0;
        tpl_ctrls->use_pred_sad_in_inter_search = 0;
        tpl_ctrls->dispenser_search_level = 0;
        tpl_ctrls->subsample_tx = 0;
        tpl_ctrls->subpel_depth = FULL_PEL;
        break;
    case 4:
        tpl_ctrls->enable = 1;
        tpl_ctrls->compute_rate = 0;
        tpl_ctrls->enable_tpl_qps = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end = DC_PRED;
        tpl_ctrls->reduced_tpl_group = is_islice ? -1 : 3;
        tpl_ctrls->pf_shape = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_pred_sad_in_intra_search = 1;
        tpl_ctrls->use_pred_sad_in_inter_search = 1;
        tpl_ctrls->dispenser_search_level = 0;
        tpl_ctrls->subsample_tx = 0;
        tpl_ctrls->subpel_depth = FULL_PEL;
        break;
    case 5:
        tpl_ctrls->enable = 1;
        tpl_ctrls->compute_rate = 0;
        tpl_ctrls->enable_tpl_qps = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end = DC_PRED;
        tpl_ctrls->reduced_tpl_group = pcs->hierarchical_levels == 5
            ? (is_islice ? 4 : 3)
            : (is_islice ? 3 : 2);
        tpl_ctrls->pf_shape = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_pred_sad_in_intra_search = 1;
        tpl_ctrls->use_pred_sad_in_inter_search = 1;
        tpl_ctrls->dispenser_search_level = 0;
        tpl_ctrls->subsample_tx = 0;
        tpl_ctrls->subpel_depth = FULL_PEL;
        break;
    case 6:
        tpl_ctrls->enable = 1;
        tpl_ctrls->compute_rate = 0;
        tpl_ctrls->enable_tpl_qps = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end = DC_PRED;
        tpl_ctrls->reduced_tpl_group = pcs->hierarchical_levels == 5
            ? is_islice ? 4 : (resolution <= INPUT_SIZE_480p_RANGE ? 3 : 2)
            : is_islice ? 3 : (resolution <= INPUT_SIZE_480p_RANGE ? 2 : 1);
        tpl_ctrls->pf_shape = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_pred_sad_in_intra_search = 1;
        tpl_ctrls->use_pred_sad_in_inter_search = 1;
        tpl_ctrls->dispenser_search_level = 0;
        tpl_ctrls->subsample_tx = 1;
        tpl_ctrls->subpel_depth = FULL_PEL;
        break;
    case 7:
        tpl_ctrls->enable = 1;
        tpl_ctrls->compute_rate = 0;
        tpl_ctrls->enable_tpl_qps = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end = DC_PRED;
        tpl_ctrls->reduced_tpl_group = pcs->hierarchical_levels == 5
            ? is_islice ? 4 : (resolution <= INPUT_SIZE_480p_RANGE ? 3 : 1)
            : is_islice ? 3 : (resolution <= INPUT_SIZE_480p_RANGE ? 2 : 0);
        tpl_ctrls->pf_shape = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_pred_sad_in_intra_search = 1;
        tpl_ctrls->use_pred_sad_in_inter_search = 1;
        tpl_ctrls->dispenser_search_level = 1;
        tpl_ctrls->subsample_tx = 2;
        tpl_ctrls->subpel_depth = FULL_PEL;
        break;
    }

    // Check user-defined settings for MAX intra mode
    if (scs->enable_paeth == 0)
        tpl_ctrls->intra_mode_end = MIN(tpl_ctrls->intra_mode_end, SMOOTH_H_PRED);

    if (scs->enable_smooth == 0)
        tpl_ctrls->intra_mode_end = MIN(tpl_ctrls->intra_mode_end, D67_PRED);

    // Derive synthesizer block size from frame size and tpl level
    tpl_ctrls->synth_blk_size = svt_aom_get_tpl_synthesizer_block_size(tpl_level, pcs->aligned_width, pcs->aligned_height);

    if ((int)scs->static_config.hierarchical_levels <= tpl_ctrls->reduced_tpl_group)
        tpl_ctrls->reduced_tpl_group = -1;

    // TPL may only look at a subset of available pictures in tpl group, which may affect the r0 calcuation.
    // As a result, we defined a factor to adjust r0 (to compensate for TPL not using all available frames).
    if (tpl_ctrls->reduced_tpl_group >= 0) {
        switch ((pcs->hierarchical_levels - tpl_ctrls->reduced_tpl_group)) {
        case 0:
        default:
            tpl_ctrls->r0_adjust_factor = 0;
            break;
        case 1:
            tpl_ctrls->r0_adjust_factor = 0.1;
            break;
        case 2:
            tpl_ctrls->r0_adjust_factor = 0.3;
            break;
        case 3:
            tpl_ctrls->r0_adjust_factor = 0.7;
            break;
        case 4:
            tpl_ctrls->r0_adjust_factor = 2.0;
            break;
        case 5:
            tpl_ctrls->r0_adjust_factor = 3.0;
            break;
        }

        // Adjust r0 scaling factor based on GOP structure and lookahead
        if (!scs->tpl_lad_mg)
            tpl_ctrls->r0_adjust_factor *= 3;
    }
    else {
        // No r0 adjustment when all frames are used
        tpl_ctrls->r0_adjust_factor = 0;

        // If no lookahead, apply r0 scaling
        if (!scs->tpl_lad_mg) {
            tpl_ctrls->r0_adjust_factor = is_islice ? 0 : 0.1;
        }
    }

    // Calculated qindex based on r0 using qstep calculation
    tpl_ctrls->qstep_based_q_calc = 1;
}

uint8_t pf_gi[16] = { 0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60 };

void set_cdef_controls(PictureParentControlSet *pcs_ptr, uint8_t cdef_level, Bool fast_decode) {

    CdefControls *cdef_ctrls = &pcs_ptr->cdef_ctrls;
    int i, j, sf_idx, second_pass_fs_num;
    cdef_ctrls->use_reference_cdef_fs = 0;
    cdef_ctrls->use_skip_detector = 0;
    switch (cdef_level)
    {
        // OFF
    case 0:
        cdef_ctrls->enabled = 0;
        break;
    case 1:
        // pf_set {0,1,..,15}
        // sf_set {0,1,2,3}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 16;
        second_pass_fs_num = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[1];
        cdef_ctrls->default_first_pass_fs[2] = pf_gi[2];
        cdef_ctrls->default_first_pass_fs[3] = pf_gi[3];
        cdef_ctrls->default_first_pass_fs[4] = pf_gi[4];
        cdef_ctrls->default_first_pass_fs[5] = pf_gi[5];
        cdef_ctrls->default_first_pass_fs[6] = pf_gi[6];
        cdef_ctrls->default_first_pass_fs[7] = pf_gi[7];
        cdef_ctrls->default_first_pass_fs[8] = pf_gi[8];
        cdef_ctrls->default_first_pass_fs[9] = pf_gi[9];
        cdef_ctrls->default_first_pass_fs[10] = pf_gi[10];
        cdef_ctrls->default_first_pass_fs[11] = pf_gi[11];
        cdef_ctrls->default_first_pass_fs[12] = pf_gi[12];
        cdef_ctrls->default_first_pass_fs[13] = pf_gi[13];
        cdef_ctrls->default_first_pass_fs[14] = pf_gi[14];
        cdef_ctrls->default_first_pass_fs[15] = pf_gi[15];
        sf_idx = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 2: // N
        // pf_set {0,1,2,4,5,6,8,9,10,12,13,14}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 12;
        second_pass_fs_num = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[1];
        cdef_ctrls->default_first_pass_fs[2] = pf_gi[2];
        cdef_ctrls->default_first_pass_fs[3] = pf_gi[4];
        cdef_ctrls->default_first_pass_fs[4] = pf_gi[5];
        cdef_ctrls->default_first_pass_fs[5] = pf_gi[6];
        cdef_ctrls->default_first_pass_fs[6] = pf_gi[8];
        cdef_ctrls->default_first_pass_fs[7] = pf_gi[9];
        cdef_ctrls->default_first_pass_fs[8] = pf_gi[10];
        cdef_ctrls->default_first_pass_fs[9] = pf_gi[12];
        cdef_ctrls->default_first_pass_fs[10] = pf_gi[13];
        cdef_ctrls->default_first_pass_fs[11] = pf_gi[14];
        sf_idx = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = -1; // cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 3:
        // pf_set {0,2,4,6,8,10,12,14}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 8;
        second_pass_fs_num = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[2];
        cdef_ctrls->default_first_pass_fs[2] = pf_gi[4];
        cdef_ctrls->default_first_pass_fs[3] = pf_gi[6];
        cdef_ctrls->default_first_pass_fs[4] = pf_gi[8];
        cdef_ctrls->default_first_pass_fs[5] = pf_gi[10];
        cdef_ctrls->default_first_pass_fs[6] = pf_gi[12];
        cdef_ctrls->default_first_pass_fs[7] = pf_gi[14];
        sf_idx = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 4:
        // pf_set {0,4,8,12,15}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 5;
        second_pass_fs_num = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[4];
        cdef_ctrls->default_first_pass_fs[2] = pf_gi[8];
        cdef_ctrls->default_first_pass_fs[3] = pf_gi[12];
        cdef_ctrls->default_first_pass_fs[4] = pf_gi[15];
        sf_idx = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = -1; // cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 5:
        // pf_set {0,5,10,15}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 4;
        second_pass_fs_num = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[5];
        cdef_ctrls->default_first_pass_fs[2] = pf_gi[10];
        cdef_ctrls->default_first_pass_fs[3] = pf_gi[15];
        sf_idx = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 6:
        // pf_set {0,7,15}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 3;
        second_pass_fs_num = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[7];
        cdef_ctrls->default_first_pass_fs[2] = pf_gi[15];
        sf_idx = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = -1; // cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 7:
        // pf_set {0,7,15}
        // sf_set {0,1,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 3;
        second_pass_fs_num = 2;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[7];
        cdef_ctrls->default_first_pass_fs[2] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 1;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[2] = pf_gi[7] + 1;
        cdef_ctrls->default_second_pass_fs[3] = pf_gi[7] + 2;
        cdef_ctrls->default_second_pass_fs[4] = pf_gi[15] + 1;
        cdef_ctrls->default_second_pass_fs[5] = pf_gi[15] + 2;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 8:
        // pf_set {0,7,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 3;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[7];
        cdef_ctrls->default_first_pass_fs[2] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[7] + 2;
        cdef_ctrls->default_second_pass_fs[2] = pf_gi[15] + 2;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = -1; // cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
        case 9:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;
        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1; // cdef_ctrls->default_second_pass_fs[0];
        cdef_ctrls->default_second_pass_fs_uv[1] = -1; // cdef_ctrls->default_second_pass_fs[1];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 4;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 10:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// when using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;

        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 1;
        cdef_ctrls->subsampling_factor = 4;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 11:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;
        cdef_ctrls->use_reference_cdef_fs = 1;
        cdef_ctrls->search_best_ref_fs = 1;
        cdef_ctrls->subsampling_factor = 4;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
        case 12:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;
        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1; // cdef_ctrls->default_second_pass_fs[0];
        cdef_ctrls->default_second_pass_fs_uv[1] = -1; // cdef_ctrls->default_second_pass_fs[1];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 4;
        cdef_ctrls->use_skip_detector = 0;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 13:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// when using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;

        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 1;
        cdef_ctrls->subsampling_factor = 4;
        cdef_ctrls->use_skip_detector = 0;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 14:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;
        cdef_ctrls->use_reference_cdef_fs = 1;
        cdef_ctrls->search_best_ref_fs = 1;
        cdef_ctrls->subsampling_factor = 4;
        cdef_ctrls->use_skip_detector = 0;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs_ptr->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 15:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;
        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1; // cdef_ctrls->default_second_pass_fs[0];
        cdef_ctrls->default_second_pass_fs_uv[1] = -1; // cdef_ctrls->default_second_pass_fs[1];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 0;
        cdef_ctrls->subsampling_factor = 4;
        cdef_ctrls->zero_fs_cost_bias = 62;
        cdef_ctrls->use_skip_detector = 0;
        break;
    case 16:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// when using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;

        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs = 1;
        cdef_ctrls->subsampling_factor = 4;
        cdef_ctrls->zero_fs_cost_bias = 62;
        cdef_ctrls->use_skip_detector = 1;
        break;
    case 17:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 2;
        second_pass_fs_num = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1] = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1] = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2] = -1;// if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;
        cdef_ctrls->use_reference_cdef_fs = 1;
        cdef_ctrls->search_best_ref_fs = 1;
        cdef_ctrls->subsampling_factor = 4;
        cdef_ctrls->zero_fs_cost_bias = 62;
        cdef_ctrls->use_skip_detector = 1;
        break;
    case 18:
        // pf_set {0}
        // sf_set {0}
        cdef_ctrls->enabled = 1;
        cdef_ctrls->first_pass_fs_num = 1;
        second_pass_fs_num = 0;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num*second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0] = pf_gi[0];
        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->use_reference_cdef_fs = 1;
        cdef_ctrls->search_best_ref_fs = 1;
        cdef_ctrls->subsampling_factor = 4;
        cdef_ctrls->zero_fs_cost_bias = 62;
        cdef_ctrls->use_skip_detector = 1;
        break;
    default:
        assert(0);
        break;
    }
}

void set_wn_filter_ctrls(Av1Common* cm, uint8_t wn_filter_lvl) {
    WnFilterCtrls* ctrls = &cm->wn_filter_ctrls;

    switch (wn_filter_lvl) {
    case 0:
        ctrls->enabled = 0;
        break;
    case 1:
        ctrls->enabled = 1;
        ctrls->use_chroma = 1;
        ctrls->filter_tap_lvl = 1;
        ctrls->use_refinement = 1;
        ctrls->max_one_refinement_step = 0;
        ctrls->use_prev_frame_coeffs = 0;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->use_chroma = 1;
        ctrls->filter_tap_lvl = 1;
        ctrls->use_refinement = 1;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs = 0;
        break;
    case 3:
        ctrls->enabled = 1;
        ctrls->use_chroma = 1;
        ctrls->filter_tap_lvl = 2;
        ctrls->use_refinement = 1;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs = 0;
        break;
    case 4:
        ctrls->enabled = 1;
        ctrls->use_chroma = 1;
        ctrls->filter_tap_lvl = 2;
        ctrls->use_refinement = 0;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs = 0;
        break;
    case 5:
        ctrls->enabled = 1;
        ctrls->use_chroma = 0;
        ctrls->filter_tap_lvl = 2;
        ctrls->use_refinement = 0;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs = 0;
        break;
    case 6:
        ctrls->enabled = 1;
        ctrls->use_chroma = 0;
        ctrls->filter_tap_lvl = 2;
        ctrls->use_refinement = 0;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs = 1;
        break;
    default:
        assert(0);
        break;
    }
}

void set_sg_filter_ctrls(Av1Common* cm, uint8_t sg_filter_lvl) {
    SgFilterCtrls* ctrls = &cm->sg_filter_ctrls;

    switch (sg_filter_lvl) {
    case 0:
        ctrls->enabled = 0;
        break;
    case 1:
        ctrls->enabled = 1;
        ctrls->use_chroma = 1;
        ctrls->step_range = 16;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->use_chroma = 1;
        ctrls->step_range = 4;
        break;
    case 3:
        ctrls->enabled = 1;
        ctrls->use_chroma = 1;
        ctrls->step_range = 1;
        break;
    case 4:
        ctrls->enabled = 1;
        ctrls->use_chroma = 1;
        ctrls->step_range = 0;
        break;
    default:
        assert(0);
        break;
    }
}

// Returns the level for Wiener restoration filter
uint8_t get_wn_filter_level(EncMode enc_mode, uint8_t input_resolution, Bool is_ref) {
    uint8_t wn_filter_lvl = 0;
    if (enc_mode <= ENC_M4)
        wn_filter_lvl = 1;
    else if (enc_mode <= ENC_M5)
        wn_filter_lvl = 2;
    else if (enc_mode <= ENC_M9)
        wn_filter_lvl = is_ref ? 5 : 0;
    else
        wn_filter_lvl = 0;

    // higher resolutions will shut restoration to save memory
    if (input_resolution >= INPUT_SIZE_8K_RANGE)
        wn_filter_lvl = 0;

    return wn_filter_lvl;
}

// Returns the level for self-guided restoration filter
uint8_t get_sg_filter_level(EncMode enc_mode, Bool fast_decode, uint8_t input_resolution, Bool is_base) {
    uint8_t sg_filter_lvl = 0;
    if (fast_decode == 0) {
        if (enc_mode <= ENC_M3)
            sg_filter_lvl = 1;
        else if (enc_mode <= ENC_M4)
            sg_filter_lvl = is_base ? 1 : 4;
        else if (enc_mode <= ENC_M5)
            sg_filter_lvl = is_base ? 1 : 0;
        else
            sg_filter_lvl = 0;
    }
    else {
        if (enc_mode <= ENC_M2)
            sg_filter_lvl = input_resolution <= INPUT_SIZE_360p_RANGE ? 1 : 0;
        else if (enc_mode <= ENC_M4)
            sg_filter_lvl = input_resolution <= INPUT_SIZE_360p_RANGE ? (is_base ? 1 : 4)
                                                                       : 0;
        else if (enc_mode <= ENC_M5)
            sg_filter_lvl = input_resolution <= INPUT_SIZE_360p_RANGE ? (is_base ? 1 : 0) : 0;
        else
            sg_filter_lvl = 0;
    }

    // higher resolutions will shut restoration to save memory
    if (input_resolution >= INPUT_SIZE_8K_RANGE)
        sg_filter_lvl = 0;

    return sg_filter_lvl;
}

void set_list0_only_base(PictureParentControlSet* pcs_ptr, uint8_t list0_only_base) {
    List0OnlyBase* ctrls = &pcs_ptr->list0_only_base_ctrls;

    switch (list0_only_base) {
    case 0:
        ctrls->enabled = 0;
        ctrls->noise_variance_th = 0;
        break;
    case 1:
        ctrls->enabled = 1;
        ctrls->noise_variance_th = 1000;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->noise_variance_th = (uint16_t)~0;
        break;
    default:
        assert(0);
        break;
    }
}

/*
* return the DLF level
* Used by signal_derivation_multi_processes_oq and memory allocation
*
* 0 -- DLF OFF
* 1 -- Full Frame DLF
* 2 -- SB-Based DLF
* 3 -- SB-Based DLF + skip DLF if SB has 0 coeffs
*/
static uint8_t get_dlf_level(EncMode enc_mode, uint8_t is_used_as_reference_flag, uint8_t is_16bit, Bool fast_decode, EbInputResolution resolution) {

    uint8_t dlf_level;
    // Don't disable DLF for low resolutions when fast-decode is used
    if (fast_decode == 0 || resolution <= INPUT_SIZE_360p_RANGE) {
        if (enc_mode <= ENC_M5)
        dlf_level = 1;
        else if (enc_mode <= ENC_M7)
            dlf_level = 2;
        else if (enc_mode <= ENC_M8)
            dlf_level = resolution <= INPUT_SIZE_360p_RANGE ? 2 : (is_used_as_reference_flag ? 2 : 4);
        else if (enc_mode <= ENC_M11)
            dlf_level = resolution <= INPUT_SIZE_360p_RANGE ? 2 : (is_used_as_reference_flag ? 3 : 0);
        else if (enc_mode <= ENC_M12)
            dlf_level = is_used_as_reference_flag ? 4 : 0;
        else
            dlf_level = (is_16bit && is_used_as_reference_flag) ? 4 : 0;
    }
    else {
        dlf_level = is_used_as_reference_flag ? 4 : 0;
    }

    return dlf_level;
}
void svt_aom_set_dlf_controls(PictureParentControlSet* pcs_ptr, uint8_t dlf_level, uint8_t bit_depth) {
    DlfCtrls* ctrls = &pcs_ptr->dlf_ctrls;

    switch (dlf_level) {
    case 0:
        ctrls->enabled = 0;
        ctrls->sb_based_dlf = 0;
        break;
    case 1:
        ctrls->enabled = 1;
        ctrls->sb_based_dlf = 0;
        ctrls->min_filter_level = 0;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->sb_based_dlf = 1;
        ctrls->min_filter_level = 0;
        break;
    case 3:
        ctrls->enabled = 1;
        ctrls->sb_based_dlf = 1;
        ctrls->min_filter_level = (bit_depth == EB_EIGHT_BIT) ? 4 : 16;
        break;
    case 4:
        ctrls->enabled = 1;
        ctrls->sb_based_dlf = 1;
        ctrls->min_filter_level = (bit_depth == EB_EIGHT_BIT) ? 16 : 32;
        break;
    default:
        assert(0);
        break;
    }
}

uint16_t  get_max_can_count(EncMode enc_mode );

/*
    set controls for intra block copy
*/
static void set_intrabc_level(PictureParentControlSet* pcs, SequenceControlSet* scs, uint8_t ibc_level) {
    IntraBCCtrls* intraBC_ctrls = &pcs->intraBC_ctrls;

    switch (ibc_level) {
    case 0:
        intraBC_ctrls->enabled = 0;
        break;
    case 1:
        intraBC_ctrls->enabled = pcs->sc_class1;
        intraBC_ctrls->ibc_shift = 0;
        intraBC_ctrls->ibc_direction = 0;
        intraBC_ctrls->hash_4x4_blocks = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = scs->super_block_size;
        break;
    case 2:
        intraBC_ctrls->enabled = pcs->sc_class1;
        intraBC_ctrls->ibc_shift = 1;
        intraBC_ctrls->ibc_direction = 0;
        intraBC_ctrls->hash_4x4_blocks = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = scs->super_block_size;
        break;
    case 3:
        intraBC_ctrls->enabled = pcs->sc_class1;
        intraBC_ctrls->ibc_shift = 1;
        intraBC_ctrls->ibc_direction = 1;
        intraBC_ctrls->hash_4x4_blocks = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = scs->super_block_size;
        break;
    case 4:
        intraBC_ctrls->enabled = pcs->sc_class1;
        intraBC_ctrls->ibc_shift = 1;
        intraBC_ctrls->ibc_direction = 0;
        intraBC_ctrls->hash_4x4_blocks = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = block_size_wide[BLOCK_16X16];
        break;
    case 5:
        intraBC_ctrls->enabled = pcs->sc_class1;
        intraBC_ctrls->ibc_shift = 1;
        intraBC_ctrls->ibc_direction = 0;
        intraBC_ctrls->hash_4x4_blocks = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = block_size_wide[BLOCK_8X8];
        break;
    case 6:
        intraBC_ctrls->enabled = pcs->sc_class1;
        intraBC_ctrls->ibc_shift = 1;
        intraBC_ctrls->ibc_direction = 1;
        intraBC_ctrls->hash_4x4_blocks = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = block_size_wide[BLOCK_8X8];
        break;
    default:
        assert(0);
        break;
    }
}
/*
    set controls for Palette prediction
*/
void set_palette_level(PictureParentControlSet* pcs_ptr, uint8_t palette_level) {

    PaletteCtrls* palette_ctrls = &pcs_ptr->palette_ctrls;

    switch (palette_level) {
    case 0:
        palette_ctrls->enabled = 0;
        break;
    case 1:
        palette_ctrls->enabled = 1;
        palette_ctrls->dominant_color_step = 1;
        break;
    case 2:
        palette_ctrls->enabled = 1;
        palette_ctrls->dominant_color_step = 2;
        break;
    default:
        assert(0);
        break;
    }
}
/******************************************************
* GM controls
******************************************************/
void set_gm_controls(PictureParentControlSet *pcs_ptr, uint8_t gm_level)
{
    GmControls *gm_ctrls = &pcs_ptr->gm_ctrls;
    switch (gm_level)
    {
    case 0:
        gm_ctrls->enabled = 0;
        break;
    case 1:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 0;
        gm_ctrls->rotzoom_model_only = 0;
        gm_ctrls->bipred_only = 0;
        gm_ctrls->bypass_based_on_me = 0;
        gm_ctrls->use_stationary_block = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps = 5;
        gm_ctrls->downsample_level = GM_FULL;
        break;
    case 2:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 0;
        gm_ctrls->bipred_only = 0;
        gm_ctrls->bypass_based_on_me = 0;
        gm_ctrls->use_stationary_block = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps = 5;
        gm_ctrls->downsample_level = GM_FULL;
        break;
    case 3:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 1;
        gm_ctrls->bipred_only = 0;
        gm_ctrls->bypass_based_on_me = 0;
        gm_ctrls->use_stationary_block = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps = 5;
        gm_ctrls->downsample_level = GM_FULL;
        break;
    case 4:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 1;
        gm_ctrls->bipred_only = 0;
        gm_ctrls->bypass_based_on_me = 1;
        gm_ctrls->use_stationary_block = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps = 5;
        gm_ctrls->downsample_level = GM_ADAPT_0;
        break;
    case 5:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 1;
        gm_ctrls->bipred_only = 0;
        gm_ctrls->bypass_based_on_me = 1;
        gm_ctrls->use_stationary_block = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps = 5;
        gm_ctrls->downsample_level = GM_ADAPT_1;
        break;
    case 6:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 1;
        gm_ctrls->bipred_only = 1;
        gm_ctrls->bypass_based_on_me = 1;
        gm_ctrls->use_stationary_block = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps = 5;
        gm_ctrls->downsample_level = GM_DOWN16;
        break;
    case 7:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 1;
        gm_ctrls->bipred_only = 1;
        gm_ctrls->bypass_based_on_me = 1;
        gm_ctrls->use_stationary_block = 1;
        gm_ctrls->use_distance_based_active_th = 1;
        gm_ctrls->params_refinement_steps = 5;
        gm_ctrls->downsample_level = GM_DOWN16;
        break;
    case 8:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 1;
        gm_ctrls->bipred_only = 1;
        gm_ctrls->bypass_based_on_me = 1;
        gm_ctrls->use_stationary_block = 1;
        gm_ctrls->use_distance_based_active_th = 1;
        gm_ctrls->params_refinement_steps = 1;
        gm_ctrls->downsample_level = GM_DOWN16;
        break;
    default:
        assert(0);
        break;
    }
}
uint8_t derive_gm_level(PictureParentControlSet* pcs_ptr) {
    SequenceControlSet* scs_ptr = pcs_ptr->scs_ptr;
    uint8_t gm_level = 0;
    const EncMode enc_mode = pcs_ptr->enc_mode;
    const uint8_t is_ref = pcs_ptr->is_used_as_reference_flag;

    // disable global motion when reference scaling enabled,
    // even if current pic is not scaled, because its reference
    // pics might be scaled in different size
    // super-res is ok for its reference pics are always upscaled
    // to original size
    if (scs_ptr->enable_global_motion == TRUE &&
        pcs_ptr->frame_superres_enabled == FALSE &&
        scs_ptr->static_config.resize_mode == RESIZE_NONE) {
        if (pcs_ptr->enc_mode <= ENC_MRS)
            gm_level = 2;
        else if (pcs_ptr->enc_mode <= ENC_M2)
            gm_level = 3;
        else if (enc_mode <= ENC_M3)
            gm_level = is_ref ? 4 : 0;
        else if (pcs_ptr->enc_mode <= ENC_M5)
            gm_level = is_ref ? 5 : 0;
        else if (pcs_ptr->enc_mode <= ENC_M6)
            gm_level = is_ref ? 6 : 0;
        else
            gm_level = 0;
    }
    return gm_level;
}

Bool is_pic_skipped(PictureParentControlSet *pcs) {
    if (!pcs->is_used_as_reference_flag &&
        pcs->scs_ptr->rc_stat_gen_pass_mode &&
        !pcs->first_frame_in_minigop)
        return TRUE;
    return FALSE;
}

/******************************************************
* Derive Multi-Processes Settings for OQ
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/
EbErrorType signal_derivation_multi_processes_oq(
    SequenceControlSet *scs_ptr,
    PictureParentControlSet *pcs_ptr,
    PictureDecisionContext *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;
    EncMode                enc_mode = pcs_ptr->enc_mode;
    const uint8_t            is_islice = pcs_ptr->slice_type == I_SLICE;
    const uint8_t            is_ref = pcs_ptr->is_used_as_reference_flag;
    const uint8_t            is_base = pcs_ptr->temporal_layer_index == 0;
    const EbInputResolution  input_resolution = pcs_ptr->input_resolution;
    const uint32_t           bit_depth = scs_ptr->static_config.encoder_bit_depth;
    const Bool               fast_decode = scs_ptr->static_config.fast_decode;
    const uint32_t           hierarchical_levels = scs_ptr->static_config.hierarchical_levels;

    // If enabled here, the hme enable flags should also be enabled in ResourceCoordinationProcess
    // to ensure that resources are allocated for the downsampled pictures used in HME
    pcs_ptr->enable_hme_flag        = 1;
    pcs_ptr->enable_hme_level0_flag = 1;
    if (pcs_ptr->sc_class1) {
        pcs_ptr->enable_hme_level1_flag = 1;
        pcs_ptr->enable_hme_level2_flag = 1;
    }
    else if (enc_mode <= ENC_M7) {
        pcs_ptr->enable_hme_level1_flag = 1;
        pcs_ptr->enable_hme_level2_flag = 1;
    }
    else {
        pcs_ptr->enable_hme_level1_flag = 1;
        pcs_ptr->enable_hme_level2_flag = 0;
    }
    switch (pcs_ptr->tf_ctrls.hme_me_level) {
    case 0:
        pcs_ptr->tf_enable_hme_flag        = 1;
        pcs_ptr->tf_enable_hme_level0_flag = 1;
        pcs_ptr->tf_enable_hme_level1_flag = 1;
        pcs_ptr->tf_enable_hme_level2_flag = 1;
        break;

    case 1:
        pcs_ptr->tf_enable_hme_flag        = 1;
        pcs_ptr->tf_enable_hme_level0_flag = 1;
        pcs_ptr->tf_enable_hme_level1_flag = 1;
        pcs_ptr->tf_enable_hme_level2_flag = 0;
        break;

    case 2:
        pcs_ptr->tf_enable_hme_flag       =  1;
        pcs_ptr->tf_enable_hme_level0_flag = 1;
        pcs_ptr->tf_enable_hme_level1_flag = 0;
        pcs_ptr->tf_enable_hme_level2_flag = 0;
        break;

    default:
        assert(0);
        break;
    }
    // Set the Multi-Pass PD level
    pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_ON;
    pcs_ptr->disallow_nsq = svt_aom_get_disallow_nsq(enc_mode, is_islice);
    // Set disallow_all_nsq_blocks_below_8x8: 8x4, 4x8


    pcs_ptr->disallow_all_nsq_blocks_below_8x8 = FALSE;


    // Set disallow_all_nsq_blocks_below_16x16: 16x8, 8x16, 16x4, 4x16
    pcs_ptr->disallow_all_nsq_blocks_below_16x16 = FALSE;

    pcs_ptr->disallow_all_nsq_blocks_below_64x64 = FALSE;
    pcs_ptr->disallow_all_nsq_blocks_below_32x32 = FALSE;
    pcs_ptr->disallow_all_nsq_blocks_above_64x64= FALSE;
    // disallow_all_nsq_blocks_above_32x32
    pcs_ptr->disallow_all_nsq_blocks_above_32x32 = FALSE;
    // disallow_all_nsq_blocks_above_16x16
    pcs_ptr->disallow_all_nsq_blocks_above_16x16 = FALSE;

    // Set if HA/HB/VA/VB blocks should considered
    if (enc_mode <= ENC_M0)
        pcs_ptr->disallow_HVA_HVB = FALSE;
    else
        pcs_ptr->disallow_HVA_HVB = TRUE;

    // Set if H4/V4 blocks should considered
    if (enc_mode <= ENC_M3)
        pcs_ptr->disallow_HV4 = FALSE;
    else
        pcs_ptr->disallow_HV4 = is_base ? FALSE : TRUE;

    // Set disallow_all_non_hv_nsq_blocks_below_16x16
    pcs_ptr->disallow_all_non_hv_nsq_blocks_below_16x16 = FALSE;

    // Set disallow_all_h4_v4_blocks_below_16x16
    pcs_ptr->disallow_all_h4_v4_blocks_below_16x16 = FALSE;
    // Loop filter Level                            Settings
    // 0                                            OFF
    // 1                                            CU-BASED
    // 2                                            LIGHT FRAME-BASED
    // 3                                            FULL FRAME-BASED

    //for now only I frames are allowed to use sc tools.
    //TODO: we can force all frames in GOP with the same detection status of leading I frame.
    uint8_t intrabc_level = 0;
    if (is_islice) {
        frm_hdr->allow_screen_content_tools = pcs_ptr->sc_class1;
        if (scs_ptr->intrabc_mode == DEFAULT) {
            if (enc_mode <= ENC_M5)
                intrabc_level = 1;
            else if (enc_mode <= ENC_M7)
                intrabc_level = 4;
            else if (enc_mode <= ENC_M9)
                intrabc_level = 5;
            else if (enc_mode <= ENC_M10)
                intrabc_level = 6;
            else
                intrabc_level = 0;
        }
        else {
            intrabc_level = (uint8_t)scs_ptr->intrabc_mode;
        }
    }
    else {
        //this will enable sc tools for P frames. hence change Bitstream even if palette mode is OFF
        frm_hdr->allow_screen_content_tools = pcs_ptr->sc_class1;
        intrabc_level = 0;
    }
    set_intrabc_level(pcs_ptr, scs_ptr, intrabc_level);
    frm_hdr->allow_intrabc = pcs_ptr->intraBC_ctrls.enabled;
    // Set palette_level
    if (frm_hdr->allow_screen_content_tools) {
        if (scs_ptr->palette_level == DEFAULT) { //auto mode; if not set by cfg
            if (enc_mode <= ENC_M11)
                pcs_ptr->palette_level = is_base ? 2 : 0;
            else if (enc_mode <= ENC_M12)
                pcs_ptr->palette_level = is_islice ? 2 : 0;
            else
                pcs_ptr->palette_level = 0;
        }
        else
            pcs_ptr->palette_level = scs_ptr->palette_level;
    }
    else
        pcs_ptr->palette_level = 0;

    set_palette_level(pcs_ptr, pcs_ptr->palette_level);
    uint8_t dlf_level = 0;
    if (pcs_ptr->scs_ptr->static_config.enable_dlf_flag && frm_hdr->allow_intrabc == 0) {
        dlf_level = get_dlf_level(enc_mode, is_ref, bit_depth > EB_EIGHT_BIT, fast_decode, input_resolution);
    }
    svt_aom_set_dlf_controls(pcs_ptr, dlf_level, bit_depth);

    // Set CDEF controls
    if (scs_ptr->seq_header.cdef_level && frm_hdr->allow_intrabc == 0) {
        if (scs_ptr->static_config.cdef_level == DEFAULT) {
            if (enc_mode <= ENC_M1)
                pcs_ptr->cdef_level = 1;
            else if (enc_mode <= ENC_M4)
                pcs_ptr->cdef_level = 2;
            else if (enc_mode <= ENC_M5)
                pcs_ptr->cdef_level = 4;
            else if (enc_mode <= ENC_M10)
                pcs_ptr->cdef_level = is_base ? 8 : is_ref ? 9 : 10;
            else if (enc_mode <= ENC_M11)
                pcs_ptr->cdef_level = is_base ? 15 : is_ref ? 16 : 17;
            else if (enc_mode <= ENC_M12) {
                if (input_resolution <= INPUT_SIZE_1080p_RANGE)
                    pcs_ptr->cdef_level = is_base ? 15 : is_ref ? 16 : 17;
                else
                    pcs_ptr->cdef_level = is_islice ? 15 : is_ref ? 16 : 17;
            }
            else {
                if (hierarchical_levels <= 2) {
                    pcs_ptr->cdef_level = is_islice ? 15 : is_base ? 16 : 0;
                }
                else {
                if (input_resolution <= INPUT_SIZE_1080p_RANGE)
                    pcs_ptr->cdef_level = is_base ? 15 : 0;
                else
                    pcs_ptr->cdef_level = is_islice ? 15 : 0;
                }
            }
        }
        else
            pcs_ptr->cdef_level = (int8_t)(scs_ptr->static_config.cdef_level);
    }
    else
        pcs_ptr->cdef_level = 0;
    set_cdef_controls(pcs_ptr, pcs_ptr->cdef_level, fast_decode);

    uint8_t wn = 0, sg = 0;
    // If restoration filtering is enabled at the sequence level, derive the settings used for this frame
    if (scs_ptr->seq_header.enable_restoration) {
        wn = get_wn_filter_level(enc_mode, input_resolution, is_ref);
        sg = get_sg_filter_level(enc_mode, fast_decode, input_resolution, is_base);
    }

    Av1Common* cm = pcs_ptr->av1_cm;
    set_wn_filter_ctrls(cm, wn);
    set_sg_filter_ctrls(cm, sg);

    // Set whether restoration filtering is enabled for this frame
    pcs_ptr->enable_restoration = (wn > 0 || sg > 0);


    // Set frame end cdf update mode      Settings
    // 0                                     OFF
    // 1                                     ON
    if (scs_ptr->frame_end_cdf_update == DEFAULT)
        pcs_ptr->frame_end_cdf_update_mode = 1;
    else
        pcs_ptr->frame_end_cdf_update_mode = scs_ptr->frame_end_cdf_update;

    (void)context_ptr;

    // Tune TPL for better chroma.Only for 240P. 0 is OFF
#if TUNE_CHROMA_SSIM
    pcs_ptr->tune_tpl_for_chroma = 1;
#else
    pcs_ptr->tune_tpl_for_chroma = 0;
#endif
    uint8_t list0_only_base = 0;
    if (enc_mode <= ENC_M6)
        list0_only_base = 0;
    else if (enc_mode <= ENC_M7) {
        if (hierarchical_levels <= 3)
            list0_only_base = 1;
        else
            list0_only_base = 0;
    }
    else if (enc_mode <= ENC_M8) {
        if (hierarchical_levels <= 3)
            list0_only_base = 2;
        else
            list0_only_base = 0;
    }
    else  if (enc_mode <= ENC_M9) {
        if (hierarchical_levels <= 3)
            list0_only_base = 2;
        else
            list0_only_base = 1;
    }
    else
        list0_only_base = 2;

    set_list0_only_base(pcs_ptr, list0_only_base);

    if (scs_ptr->enable_hbd_mode_decision == DEFAULT)
        if (enc_mode <= ENC_MR)
            pcs_ptr->hbd_mode_decision = 1;
        else if (enc_mode <= ENC_M1)
            pcs_ptr->hbd_mode_decision = is_ref ? 1 : 2;
        else if (enc_mode <= ENC_M3)
            pcs_ptr->hbd_mode_decision = 2;
        else if (enc_mode <= ENC_M4)
            pcs_ptr->hbd_mode_decision = is_ref ? 2 : 0;
        else if (enc_mode <= ENC_M7)
            pcs_ptr->hbd_mode_decision = is_base ? 2 : 0;
        else
            pcs_ptr->hbd_mode_decision = is_islice ? 2 : 0;
    else
        pcs_ptr->hbd_mode_decision = scs_ptr->enable_hbd_mode_decision;

    pcs_ptr->max_can_count = get_max_can_count(enc_mode);

    if (enc_mode <= ENC_M7)
        pcs_ptr->use_best_me_unipred_cand_only = 0;
    else
        pcs_ptr->use_best_me_unipred_cand_only = 1;

    //TPL level should not be modified outside of this function
    set_tpl_extended_controls(pcs_ptr, scs_ptr->tpl_level);

    pcs_ptr->adjust_under_shoot_gf = 0;
    if (scs_ptr->passes == 1 && scs_ptr->static_config.rate_control_mode == SVT_AV1_RC_MODE_VBR)
        pcs_ptr->adjust_under_shoot_gf = enc_mode <= ENC_M11 ? 1 : 2;
    return return_error;
}

/******************************************************
* Derive Multi-Processes Settings for first pass
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_multi_processes(
    SequenceControlSet *scs_ptr,
    PictureParentControlSet *pcs_ptr);
//set the ref frame types used for this picture,
static void set_all_ref_frame_type(PictureParentControlSet  *parent_pcs_ptr, MvReferenceFrame ref_frame_arr[], uint8_t* tot_ref_frames)
{
    MvReferenceFrame rf[2];
    *tot_ref_frames = 0;

    //SVT_LOG("POC %i  totRef L0:%i   totRef L1: %i\n", parent_pcs_ptr->picture_number, parent_pcs_ptr->ref_list0_count, parent_pcs_ptr->ref_list1_count);

     //single ref - List0
    for (uint8_t ref_idx0 = 0; ref_idx0 < parent_pcs_ptr->ref_list0_count_try; ++ref_idx0) {
        rf[0] = svt_get_ref_frame_type(REF_LIST_0, ref_idx0);
        ref_frame_arr[(*tot_ref_frames)++] = rf[0];
    }

    //single ref - List1
    for (uint8_t ref_idx1 = 0; ref_idx1 < parent_pcs_ptr->ref_list1_count_try; ++ref_idx1) {
        rf[1] = svt_get_ref_frame_type(REF_LIST_1, ref_idx1);
        ref_frame_arr[(*tot_ref_frames)++] = rf[1];
    }

    //compound Bi-Dir
    for (uint8_t ref_idx0 = 0; ref_idx0 < parent_pcs_ptr->ref_list0_count_try; ++ref_idx0) {
        for (uint8_t ref_idx1 = 0; ref_idx1 < parent_pcs_ptr->ref_list1_count_try; ++ref_idx1) {
            rf[0] = svt_get_ref_frame_type(REF_LIST_0, ref_idx0);
            rf[1] = svt_get_ref_frame_type(REF_LIST_1, ref_idx1);
            ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
        }
    }
    if (parent_pcs_ptr->slice_type == B_SLICE)
    {

        //compound Uni-Dir
        if (parent_pcs_ptr->ref_list0_count_try > 1) {
            rf[0] = LAST_FRAME;
            rf[1] = LAST2_FRAME;
            ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
            if (parent_pcs_ptr->ref_list0_count_try > 2) {
                rf[1] = LAST3_FRAME;
                ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
                if (parent_pcs_ptr->ref_list0_count_try > 3) {
                    rf[1] = GOLDEN_FRAME;
                    ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
                }
            }
        }
        if (parent_pcs_ptr->ref_list1_count_try > 2) {
            rf[0] = BWDREF_FRAME;
            rf[1] = ALTREF_FRAME;
            ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
        }
    }

}

static void prune_refs(PredictionStructureEntry *pred_position_ptr, Av1RpsNode *av1_rps)
{
    if (pred_position_ptr->ref_list0.reference_list_count < 4) {
        av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];
    }
    if (pred_position_ptr->ref_list0.reference_list_count < 3) {
        av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
    }
    if (pred_position_ptr->ref_list0.reference_list_count < 2) {
        av1_rps->ref_dpb_index[LAST2] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_poc_array[LAST2] = av1_rps->ref_poc_array[LAST];
    }

    if (pred_position_ptr->ref_list1.reference_list_count < 3) {
        av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
        av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
    }
    if (pred_position_ptr->ref_list1.reference_list_count < 2) {
        av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
        av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
    }
}

// Set the show_frame and show_existing_frame for current picture if it's:
// 1)Low delay P, 2)Low delay b and 3)I frames of RA
// For b frames of RA, need to set it manually based on picture_index
static Bool set_frame_display_params(
        PictureParentControlSet       *pcs_ptr,
        PictureDecisionContext        *context_ptr,
        uint32_t                       mini_gop_index)
{
    Av1RpsNode *av1_rps = &pcs_ptr->av1_ref_signal;
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;

    if (pcs_ptr->pred_struct_ptr->pred_type == SVT_AV1_PRED_LOW_DELAY_P || pcs_ptr->is_overlay) {
        //P frames
        av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_poc_array[BWD] = av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[LAST];

        frm_hdr->show_frame = TRUE;
        pcs_ptr->has_show_existing = FALSE;
    } else if (pcs_ptr->pred_struct_ptr->pred_type == SVT_AV1_PRED_LOW_DELAY_B) {
        av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[LAST2];
        av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST3];
        av1_rps->ref_poc_array[BWD] = av1_rps->ref_poc_array[LAST];
        av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[LAST2];
        av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[LAST3];

        frm_hdr->show_frame = TRUE;
        pcs_ptr->has_show_existing = FALSE;
    } else if (pcs_ptr->pred_struct_ptr->pred_type == SVT_AV1_PRED_RANDOM_ACCESS) {
        //Decide on Show Mecanism
        if (pcs_ptr->slice_type == I_SLICE) {
            //3 cases for I slice:  1:Key Frame treated above.  2: broken MiniGop due to sc or intra refresh  3: complete miniGop due to sc or intra refresh
            if (context_ptr->mini_gop_length[mini_gop_index] < pcs_ptr->pred_struct_ptr->pred_struct_period) {
                //Scene Change that breaks the mini gop and switch to LDP (if I scene change happens to be aligned with a complete miniGop, then we do not break the pred structure)
                frm_hdr->show_frame = TRUE;
                pcs_ptr->has_show_existing = FALSE;
            } else {
                frm_hdr->show_frame = FALSE;
                pcs_ptr->has_show_existing = FALSE;
            }
        } else {
            if (context_ptr->mini_gop_length[mini_gop_index] != pcs_ptr->pred_struct_ptr->pred_struct_period) {
                SVT_LOG("Error in GOP indexing3\n");
            }
            // Handle b frame of Random Access out
            return FALSE;
        }
    }
    return TRUE;
}

static void set_key_frame_rps(PictureParentControlSet *pcs, PictureDecisionContext *context_ptr)
{
    FrameHeader *frm_hdr = &pcs->frm_hdr;
    context_ptr->lay0_toggle = 0;
    context_ptr->lay1_toggle = 0;
    context_ptr->lay2_toggle = 0;

    frm_hdr->show_frame = TRUE;
    pcs->has_show_existing = FALSE;
    return;
}

// Decide whether to make an inter frame into an S-Frame
static void set_sframe_type(PictureParentControlSet *ppcs, EncodeContext *encode_context_ptr, PictureDecisionContext *context_ptr)
{
    SequenceControlSet* scs = ppcs->scs_ptr;
    FrameHeader *frm_hdr = &ppcs->frm_hdr;
    const int sframe_dist = encode_context_ptr->sf_cfg.sframe_dist;
    const EbSFrameMode sframe_mode = encode_context_ptr->sf_cfg.sframe_mode;

    // s-frame supports low-delay
    assert_err(scs->static_config.pred_structure == 0 || scs->static_config.pred_structure == 1,
        "S-frame supports only low delay");
    // handle multiple hierarchical levels only, no flat IPPP support
    assert_err(ppcs->hierarchical_levels > 0, "S-frame doesn't support flat IPPP...");

    const int is_arf = ppcs->temporal_layer_index == 0 ? TRUE : FALSE;
    const uint64_t frames_since_key = ppcs->picture_number - context_ptr->key_poc;
    if (sframe_mode == SFRAME_STRICT_BASE) {
        // SFRAME_STRICT_ARF: insert sframe if it matches altref frame.
        if (is_arf && (frames_since_key % sframe_dist) == 0) {
            frm_hdr->frame_type = S_FRAME;
        }
    }
    else {
        // SFRAME_NEAREST_ARF: if sframe will be inserted at the next available altref frame
        if ((frames_since_key % sframe_dist) == 0) {
            context_ptr->sframe_due = 1;
        }
        if (context_ptr->sframe_due && is_arf) {
            frm_hdr->frame_type = S_FRAME;
            context_ptr->sframe_due = 0;
        }
    }

#if DEBUG_SFRAME
    if (frm_hdr->frame_type == S_FRAME) {
        fprintf(stderr, "\nFrame %d - set sframe\n", (int)ppcs->picture_number);
    }
#endif
    return;
}

// Update RPS info for S-Frame
static void set_sframe_rps(PictureParentControlSet *ppcs, EncodeContext *encode_context_ptr, PictureDecisionContext *context_ptr)
{
    ppcs->frm_hdr.error_resilient_mode = 1;
    ppcs->av1_ref_signal.refresh_frame_mask = 0xFF;

    context_ptr->lay0_toggle = 0;
    context_ptr->lay1_toggle = 0;
    context_ptr->lay2_toggle = 0;

    // Bookmark latest switch frame poc to prevent following frames referencing frames before the switch frame
    context_ptr->sframe_poc = ppcs->picture_number;
    // Reset pred_struct_position
    encode_context_ptr->elapsed_non_cra_count = 0;
    return;
}

/*************************************************
* AV1 Reference Picture Signalling:
* Stateless derivation of RPS info to be stored in
* Picture Header
*
* This function uses the picture index from the just
* collected miniGop to derive the RPS(refIndexes+refresh)
* the miniGop is always 4L but could be complete (8 pictures)
or non-complete (less than 8 pictures).
* We get to this function if the picture is:
* 1) first Key frame
* 2) part of a complete RA MiniGop where the last frame could be a regular I for open GOP
* 3) part of complete LDP MiniGop where the last frame could be Key frame for closed GOP
* 4) part of non-complete LDP MiniGop where the last frame is a regularI+SceneChange.
This miniGOP has P frames with predStruct=LDP, and the last frame=I with pred struct=RA.
* 5) part of non-complete LDP MiniGop at the end of the stream.This miniGOP has P frames with
predStruct=LDP, and the last frame=I with pred struct=RA.
*
*Note: the  SceneChange I has pred_type = SVT_AV1_PRED_RANDOM_ACCESS. if SChange is aligned on the miniGop,
we do not break the GOP.
*************************************************/
static void  av1_generate_rps_info(
    PictureParentControlSet       *pcs_ptr,
    EncodeContext                 *encode_context_ptr,
    PictureDecisionContext        *context_ptr,
    uint32_t                       picture_index,
    uint32_t                       mini_gop_index
)
{
    (void)encode_context_ptr;
    Av1RpsNode *av1_rps = &pcs_ptr->av1_ref_signal;
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;

    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    PredictionStructureEntry *pred_position_ptr = pcs_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[pcs_ptr->pred_struct_index];
    //Set frame type
    if (pcs_ptr->slice_type == I_SLICE) {
        frm_hdr->frame_type = pcs_ptr->idr_flag ? KEY_FRAME : INTRA_ONLY_FRAME;
        pcs_ptr->av1_ref_signal.refresh_frame_mask = 0xFF;
#if DEBUG_SFRAME
        fprintf(stderr, "\nFrame %d - key frame\n", (int)pcs_ptr->picture_number);
#endif
    }
    else {
        frm_hdr->frame_type = INTER_FRAME;

        // test s-frame on base layer inter frames
        if (encode_context_ptr->sf_cfg.sframe_dist > 0) {
            set_sframe_type(pcs_ptr, encode_context_ptr, context_ptr);
        }
    }


    pcs_ptr->intra_only = pcs_ptr->slice_type == I_SLICE ? 1 : 0;
    // When a hierarchical change happens, the references should be chosen properly. In the following case, we switch from 6L to 5L or 4L.
    // In this case we need to update lay0_tpggle to choose the right reference
    if (pcs_ptr->hierarchical_layers_diff != 0) {
        if ((pcs_ptr->hierarchical_layers_diff == 2 && pcs_ptr->hierarchical_levels == 3) ||
            (pcs_ptr->hierarchical_layers_diff == 1 && pcs_ptr->hierarchical_levels == 4)) {
            context_ptr->lay0_toggle = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 1 : 0;
        }
    }
    if (pcs_ptr->hierarchical_levels == 0) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(pcs_ptr, context_ptr);
            return;
        }

        const uint8_t  base0_idx = (context_ptr->lay0_toggle + 8 - 1) % 8;
        const uint8_t  base1_idx = (context_ptr->lay0_toggle + 8 - 2) % 8;
        const uint8_t  base2_idx = (context_ptr->lay0_toggle + 8 - 3) % 8;
        const uint8_t  base3_idx = (context_ptr->lay0_toggle + 8 - 4) % 8;
        const uint8_t  base4_idx = (context_ptr->lay0_toggle + 8 - 5) % 8;
        const uint8_t  base5_idx = (context_ptr->lay0_toggle + 8 - 6) % 8;
        //const uint8_t  base6_idx = (context_ptr->lay0_toggle + 8 - 7) % 8;
        const uint8_t  base7_idx = (context_ptr->lay0_toggle + 8 - 8) % 8;

       // {1, 3, 5, 7},    // GOP Index 0 - Ref List 0
       // { 2, 4, 6, 0 }     // GOP Index 0 - Ref List 1
        av1_rps->ref_dpb_index[LAST] = base0_idx;
        av1_rps->ref_dpb_index[LAST2] = base2_idx;
        av1_rps->ref_dpb_index[LAST3] = base4_idx;
        av1_rps->ref_dpb_index[GOLD] = base7_idx;

        av1_rps->ref_dpb_index[BWD] = base1_idx;
        av1_rps->ref_dpb_index[ALT2] = base3_idx;
        av1_rps->ref_dpb_index[ALT] = base5_idx;
        gop_i = 0;
        av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list0[0]);
        av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list0[1]);
        av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list0[2]);
        av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list0[3]);

        av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list1[0]);
        av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list1[1]);
        av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list1[2]);

        prune_refs(pred_position_ptr, av1_rps);
        av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;

        // Flat mode, output all frames
        set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index);
        frm_hdr->show_frame = TRUE;
        pcs_ptr->has_show_existing = FALSE;
        context_ptr->lay0_toggle = (1 + context_ptr->lay0_toggle) % 8;
    } else if (pcs_ptr->hierarchical_levels == 1) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(pcs_ptr, context_ptr);
            return;
        }

        const uint8_t  base_off0_idx = (context_ptr->lay0_toggle + 5 - 5) % 5;
        const uint8_t  base_off2_idx = (context_ptr->lay0_toggle + 5 - 1) % 5;
        const uint8_t  base_off4_idx = (context_ptr->lay0_toggle + 5 - 2) % 5;
        const uint8_t  base_off6_idx = (context_ptr->lay0_toggle + 5 - 3) % 5;
        const uint8_t  base_off8_idx = (context_ptr->lay0_toggle + 5 - 4) % 5;
        //const uint8_t  lay1_off0_idx = (context_ptr->lay1_toggle + 3 - 3) % 3 + 5;
        const uint8_t  lay1_off2_idx = (context_ptr->lay1_toggle + 3 - 1) % 3 + 5;
        //const uint8_t  lay1_off4_idx = (context_ptr->lay1_toggle + 3 - 2) % 3 + 5;

        switch (pcs_ptr->temporal_layer_index) {
        case 0:
            //{2, 6, 10, 0},     // GOP Index 0 - Ref List 0
            //{4, 8, 0, 0}       // GOP Index 0 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base_off2_idx;
            av1_rps->ref_dpb_index[LAST2] = base_off6_idx;
            av1_rps->ref_dpb_index[LAST3] = base_off0_idx;
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base_off4_idx;
            av1_rps->ref_dpb_index[ALT2] = base_off8_idx;
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
            gop_i = 0;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;
            break;

        case 1:
            //{ 1, 2, 3, 5},    // GOP Index 1 - Ref List 0
            //{ -1, 0, 0, 0 }     // GOP Index 1 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base_off2_idx;
            av1_rps->ref_dpb_index[LAST2] = lay1_off2_idx;
            av1_rps->ref_dpb_index[LAST3] = base_off4_idx;
            av1_rps->ref_dpb_index[GOLD] = base_off6_idx;

            av1_rps->ref_dpb_index[BWD] = base_off0_idx;
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];

            gop_i = 1;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[3]);;

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            av1_rps->refresh_frame_mask = 1 << (context_ptr->lay1_toggle+5);

            break;

        default:
            SVT_ERROR("unexpected picture mini Gop number\n");
            break;
        }

        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {
            if (pcs_ptr->is_used_as_reference_flag && picture_index != 0) {
                frm_hdr->show_frame = FALSE;
                pcs_ptr->has_show_existing = FALSE;
            } else {
                frm_hdr->show_frame = TRUE;
                pcs_ptr->has_show_existing = TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = base_off0_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", pcs_ptr->hierarchical_levels);
            }
        }

        // Toggle layer0 and layer1
        if ((picture_index == context_ptr->mini_gop_end_index[mini_gop_index] &&
                    !pcs_ptr->is_overlay &&
                    scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS) ||
                ((scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {
            context_ptr->lay0_toggle = (1 + context_ptr->lay0_toggle) % 5;
            context_ptr->lay1_toggle = (1 + context_ptr->lay1_toggle) % 3;
        }
    } else if (pcs_ptr->hierarchical_levels == 2) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(pcs_ptr, context_ptr);
            return;
        }
        // Slot 0,1,2 for base layer, 3,4 for layer 1, 5 for layer2
        const uint8_t  base0_idx = context_ptr->lay0_toggle == 0 ? 1 : context_ptr->lay0_toggle == 1 ? 2 : 0;
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 0 : 1;
        const uint8_t  base2_idx = context_ptr->lay0_toggle == 0 ? 0 : context_ptr->lay0_toggle == 1 ? 1 : 2;

        const uint8_t  lay1_0_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 1 : LAY1_OFF + 0;
        const uint8_t  lay1_1_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 0 : LAY1_OFF + 1;
        const uint8_t  lay2_idx = LAY2_OFF;

        switch (pcs_ptr->temporal_layer_index) {
        case 0:
            //{4, 12, 0, 0},     // GOP Index 0 - Ref List 0
            //{ 4, 8, 0, 0 }      // GOP Index 0 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = base2_idx;
            av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base1_idx;
            av1_rps->ref_dpb_index[ALT2] = base0_idx;
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
            gop_i = 0;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;
            break;

        case 1:
            //{ 2, 4, 6, 0}   // GOP Index 2 - Ref List 0
            //{-2, 0, 0, 0}   // GOP Index 2 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
            av1_rps->ref_dpb_index[LAST3] = base0_idx;
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base2_idx;
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];

            gop_i = 2;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 1 << (LAY1_OFF + context_ptr->lay1_toggle);

            break;

        case 2:
            // update RPS for the overlay frame.
            if (pcs_ptr->is_overlay) {
                //{ 0, 0, 0, 0}         // GOP Index 1 - Ref List 0
                //{ 0, 0, 0, 0 }       // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST]  = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD]  = base1_idx;
                av1_rps->ref_dpb_index[BWD]  = base1_idx;
                av1_rps->ref_dpb_index[ALT2]  = base1_idx;
                av1_rps->ref_dpb_index[ALT] = base1_idx;

                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
            } else {
                if (picture_index == 0) {
                    //{ 1, 3, 5, 0}      // GOP Index 1 - Ref List 0
                    //{ -1, -3, 0, 0}    // GOP Index 1 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = base1_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
                    av1_rps->ref_dpb_index[LAST3] = base0_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

                    av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                    av1_rps->ref_dpb_index[ALT2] = base2_idx;
                    av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                    gop_i = 1;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
                }
                else if (picture_index == 2) {
                    // { 1, 3, 2, 0},     // GOP Index 3 - Ref List 0
                    // { -1,0, 0, 0}      // GOP Index 3 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay2_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

                    av1_rps->ref_dpb_index[BWD] = base2_idx;
                    av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                    av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                    gop_i = 3;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
                    av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
                }
                else
                    SVT_LOG("Error in GOp indexing\n");
            }

            if (picture_index == 0 )
                av1_rps->refresh_frame_mask = 1 << (lay2_idx);
            else
                av1_rps->refresh_frame_mask = 0;
            break;

        default:
            SVT_ERROR("unexpected picture mini Gop number\n");
            break;
        }

        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {
            if (pcs_ptr->is_used_as_reference_flag && picture_index != 0) {
                frm_hdr->show_frame = FALSE;
                pcs_ptr->has_show_existing = FALSE;
            } else {
                frm_hdr->show_frame = TRUE;
                pcs_ptr->has_show_existing = TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = lay1_1_idx;
                else if (picture_index == 2)
                    frm_hdr->show_existing_frame = base2_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", pcs_ptr->hierarchical_levels);
            }
        }

        if ((picture_index == (context_ptr->mini_gop_end_index[mini_gop_index]) &&
                    !pcs_ptr->is_overlay &&
                    scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS) ||
                ((scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_P ||
                 scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {
            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
        }
    } else if (pcs_ptr->hierarchical_levels == 3) {

        uint8_t gop_i;
        Bool is_trailing_frames = FALSE;

        if (frm_hdr->frame_type == KEY_FRAME)
        {
            set_key_frame_rps(pcs_ptr, context_ptr);
            return;
        }

        //picture_index has this order:
        //         0     2    4      6
        //            1          5
        //                 3
        //                              7(could be an I)
        //

        //DPB: Loc7|Loc6|Loc5|Loc4|Loc3|Loc2|Loc1|Loc0
        //Layer 0 : circular move 0-1-2
        //Layer 1 : circular move 3-4
        //Layer 2 : circular move 5-6
        //Layer 3 : 7
        //pic_num
        //         1     3    5      7    9     11     13      15
        //            2          6           10            14
        //                 4                        12
        //
        //base0:0                   base1:8                          base2:16
        if (pcs_ptr->pred_struct_ptr->pred_type == SVT_AV1_PRED_LOW_DELAY_P &&
                context_ptr->mini_gop_length[mini_gop_index] < 8 &&
                pcs_ptr->scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS) {
            is_trailing_frames = TRUE;
        }


        const uint8_t  base0_idx = context_ptr->lay0_toggle == 0 ? 1 : context_ptr->lay0_toggle == 1 ? 2 : 0;   //the oldest L0 picture in the DPB
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 0 : 1;   //the middle L0 picture in the DPB
        const uint8_t  base2_idx = context_ptr->lay0_toggle == 0 ? 0 : context_ptr->lay0_toggle == 1 ? 1 : 2;   //the newest L0 picture in the DPB

        const uint8_t  lay1_0_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 1 : LAY1_OFF + 0;   //the oldest L1 picture in the DPB
        const uint8_t  lay1_1_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 0 : LAY1_OFF + 1;   //the newest L1 picture in the DPB

        const uint8_t  lay2_0_idx = picture_index < 4 ? LAY2_OFF + 1 : LAY2_OFF + 0;   //the oldest L2 picture in the DPB
        const uint8_t  lay2_1_idx = picture_index < 4 ? LAY2_OFF + 0 : LAY2_OFF + 1;   //the newest L2 picture in the DPB
        const uint8_t  lay3_idx = 7;    //the newest L3 picture in the DPB
        /*in 5L struct, we switich to  4L at the end of the seq.
        the current pred struct is reset, and the previous 5L minGop is out of reach!
        four_level_hierarchical_pred_struct will be set to follow 4L order, but this will generate some RPS mismatch for some frames.
        PKTization DPB can detect those  */

        switch (pcs_ptr->temporal_layer_index) {
        case 0:

            //{8, 0, 0, 0},     // GOP Index 0 - Ref List 0
            //{8, 0,  0, 0}      // GOP Index 0 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = base0_idx;
            av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base1_idx;
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
            gop_i = 0;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;

            break;

        case 1:
            //{ 4, 8, 12,  0},   // GOP Index 4 - Ref List 0
            //{-4,  0, 0,  0}     // GOP Index 4 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
            av1_rps->ref_dpb_index[LAST3] = base0_idx;
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base2_idx;
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];

            gop_i = 4;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            av1_rps->refresh_frame_mask = 1 << (LAY1_OFF + context_ptr->lay1_toggle);

            break;

        case 2:

            if (picture_index == 1) {
                //{  2,  4,  6,  10}    // GOP Index 2 - Ref List 0
                //{ -2, -6,  0,   0}    // GOP Index 2 - Ref List 1

                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_0_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = base0_idx;

                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 2;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 5) {
                //{ 2, 4, 6, 10}   // GOP Index 6 - Ref List 0
                //{ -2,  0, 0,  0 }    // GOP Index 6 - Ref List 1

                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_0_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;// av1_rps->ref_dpb_index[LAST];

                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 6;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }

            av1_rps->refresh_frame_mask = 1 << (LAY2_OFF + context_ptr->lay2_toggle);
            //toggle 3->4
            if (!is_trailing_frames ||
                    (context_ptr->mini_gop_length[mini_gop_index] >= 7 && is_trailing_frames &&
                     pcs_ptr->scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS)) {
                //For trailing frames, Only toggle it if we are sure we have 2 layer 2 frames in trailing frames
                context_ptr->lay2_toggle = 1 - context_ptr->lay2_toggle;
            }

            break;

        case 3:
            // update RPS for the overlay frame.
            if (pcs_ptr->is_overlay) {
                //{ 0, 0, 0, 0}         // GOP Index 1 - Ref List 0
                //{ 0, 0, 0, 0 }       // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST]  = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD]  = base1_idx;

                av1_rps->ref_dpb_index[BWD]  = base1_idx;
                av1_rps->ref_dpb_index[ALT2]  = base1_idx;
                av1_rps->ref_dpb_index[ALT] = base1_idx;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
            }
            else

            if (picture_index == 0) {
                //{ 1, 3, 5, 8}         // GOP Index 1 - Ref List 0
                //{ -1, -3, -7,  0 }    // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_0_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = lay3_idx;
                av1_rps->ref_dpb_index[BWD] = lay2_1_idx;
                av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 1;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 2) {
                // { 1, 3, 2, 5},        // GOP Index 3 - Ref List 0
                //{ -1,  -5, 0,  0 }     //     GOP Index 3 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay3_idx;
                av1_rps->ref_dpb_index[GOLD] = lay2_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 3;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);
                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 4) {
               // { 1, 3, 5, 4},    // GOP Index 5 - Ref List 0

                //{ -1,  -3, 0,  0 }     // GOP Index 5 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_0_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay3_idx;

                av1_rps->ref_dpb_index[BWD] = lay2_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 5;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 6) {

                //{ 1,  3, 5,  6},     //  GOP Index 7 - Ref List 0

                //{ -1,  0, 0,  0 }      // GOP Index 7 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay2_0_idx;
                av1_rps->ref_dpb_index[GOLD] = lay3_idx;
                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 7;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else
                SVT_LOG("Error in GOp indexing\n");
            if (pcs_ptr->scs_ptr->mrp_ctrls.referencing_scheme == 1) {
                av1_rps->refresh_frame_mask = pcs_ptr->is_overlay ? 0 : 1 << (lay3_idx);
            }
            else {
              if (picture_index == 0)
                av1_rps->refresh_frame_mask = 1 << (lay3_idx);
              else
                av1_rps->refresh_frame_mask = 0;
            }
            break;

        default:
            SVT_ERROR("unexpected picture mini Gop number\n");
            break;
        }

        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {

            uint8_t no_show;
            if (pcs_ptr->scs_ptr->mrp_ctrls.referencing_scheme == 1)
                no_show = pcs_ptr->is_used_as_reference_flag && picture_index != 0 && picture_index != 2 && picture_index != 4 && picture_index != 6;
            else
                no_show = pcs_ptr->is_used_as_reference_flag && picture_index != 0;

            if (no_show)
            {
                frm_hdr->show_frame = FALSE;
                pcs_ptr->has_show_existing = FALSE;
            } else {
                frm_hdr->show_frame = TRUE;
                pcs_ptr->has_show_existing = TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = lay2_1_idx;
                else if (picture_index == 2)
                    frm_hdr->show_existing_frame = lay1_1_idx;
                else if (picture_index == 4)
                    frm_hdr->show_existing_frame = lay2_1_idx;
                else if (picture_index == 6)
                    frm_hdr->show_existing_frame = base2_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", pcs_ptr->hierarchical_levels);
            }
        }

        //last pic in MiniGop: Base layer toggling
        //mini GOP toggling since last Key Frame.
        //a regular I keeps the toggling process and does not reset the toggle.  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
        //whoever needs a miniGOP Level toggling, this is the time
        if (((picture_index == (context_ptr->mini_gop_end_index[mini_gop_index] % 8) &&
                        !pcs_ptr->is_overlay &&
                        scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS)) ||
                ((scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {

            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;


            //layer2 toggle for low delay b/P case
            if ((scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) &&
                    scs_ptr->static_config.intra_period_length != -1 &&
                    pcs_ptr->cra_flag) {
                int32_t trailing_count = (scs_ptr->static_config.intra_period_length % 8);
                if (trailing_count >=2 && trailing_count <= 5) {
                    context_ptr->lay2_toggle = 1 - context_ptr->lay2_toggle;
                }
            }
        }
    } else if (pcs_ptr->hierarchical_levels == 4) {
        uint8_t gop_i;
        //Av1RpsNode_t *av1_rps = &pcs_ptr->av1RefSignal2;

        //Reset miniGop Toggling. The first miniGop after a KEY frame has toggle=0
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(pcs_ptr, context_ptr);
            return;
        }

        //picture_index has this order:
        //         0     2    4      6    8     10     12      14
        //            1          5           9            13
        //                 3                        11
        //                              7
        //                                                          15(could be an I)

        //DPB: Loc7|Loc6|Loc5|Loc4|Loc3|Loc2|Loc1|Loc0
        //Layer 0 : circular move 0-1-2
        //Layer 1 : circular move 3-4
        //Layer 2 : DPB Location 5
        //Layer 3 : DPB Location 6
        //Layer 4 : DPB Location 7
        //pic_num                  for poc 17
        //         1     3    5      7    9     11     13      15         17    19     21    23   25     27    29    31
        //            2          6           10            14                18           22          26          30
        //                 4                        12:L2_0                         20:L2_1                 28
        //                              8:L1_0                                                       24:L1_1
        //base0:0                                               base1:16                                           base2:32

        const uint8_t  base0_idx = context_ptr->lay0_toggle == 0 ? 1 : context_ptr->lay0_toggle == 1 ? 2 : 0;   //the oldest L0 picture in the DPB
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 0 : 1;   //the middle L0 picture in the DPB
        const uint8_t  base2_idx = context_ptr->lay0_toggle == 0 ? 0 : context_ptr->lay0_toggle == 1 ? 1 : 2;   //the newest L0 picture in the DPB

        const uint8_t  lay1_0_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 1 : LAY1_OFF + 0;   //the oldest L1 picture in the DPB
        const uint8_t  lay1_1_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 0 : LAY1_OFF + 1;   //the newest L1 picture in the DPB
        const uint8_t  lay2_1_idx = LAY2_OFF;   //the oldest L2 picture in the DPB
        const uint8_t  lay3_idx = LAY3_OFF;    //the newest L3 picture in the DPB
        const uint8_t  lay4_idx = LAY4_OFF;    //the newest L3 picture in the DPB
        switch (pcs_ptr->temporal_layer_index) {
        case 0:

            //{16, 48, 0, 0},      // GOP Index 0 - Ref List 0
           //{16, 32, 0, 0}       // GOP Index 0 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = base2_idx;
            av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base1_idx;
            av1_rps->ref_dpb_index[ALT2] = base0_idx;
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];

            gop_i = 0;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;

            break;

        case 1:
            //{  8, 16, 24, 0},   // GOP Index 8 - Ref List 0
            //{ -8, 0, 0, 0}      // GOP Index 8 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
            av1_rps->ref_dpb_index[LAST3] = base0_idx;
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base2_idx;
            av1_rps->ref_dpb_index[ALT2] = lay2_1_idx;
            av1_rps->ref_dpb_index[ALT] = lay3_idx;

            gop_i = 8;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
            av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            av1_rps->refresh_frame_mask = 1 << (LAY1_OFF + context_ptr->lay1_toggle);

            break;

        case 2:

            if (picture_index == 3) {
                //{  4,   8,  12,  20 },  // GOP Index 4 - Ref List 0
                //{ -4, -12,  0,  0 }     // GOP Index 4 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = base0_idx;

                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 4;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            }
            else if (picture_index == 11) {
                //{ 4, 8, 12, 0},       // GOP Index 12 - Ref List 0
                //{ -4,  0, 0,  0 }     // GOP Index 12 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay3_idx;

                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = lay4_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 12;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            av1_rps->refresh_frame_mask = 1 << (LAY2_OFF);
            //toggle 3->4
            context_ptr->lay2_toggle = 1 - context_ptr->lay2_toggle;

            break;

        case 3:

            if (picture_index == 1) {
                //{ 2, 4, 10, 18},        // GOP Index 2 - Ref List 0
                //{ -2, -6, -14,  0 }   // GOP Index 2 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = base0_idx;
                av1_rps->ref_dpb_index[BWD] = lay2_1_idx;
                av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 2;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 5) {
                //{ 2, 4, 6, 14},        // GOP Index 6 - Ref List 0
                //{ -2, -10,  0,  0 }   // GOP Index 6 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 6;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 9) {
                //{ 2, 4, 10, 18},       // GOP Index 10 - Ref List 0
                //{ -2, -6,  0,  0 }    // GOP Index 10 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay2_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 10;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 13) {
                //{ 2, 4, 6, 14},    // GOP Index 14 - Ref List 0
                //{ -2, 0,  0, 0 }   // GOP Index 14 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_1_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;

                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = lay4_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 14;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else
                SVT_LOG("Error in GOp indexing\n");
            av1_rps->refresh_frame_mask = 1 << (lay3_idx);
            break;

        case 4:
            // update RPS for the overlay frame.
            if (pcs_ptr->is_overlay) {
                //{ 0, 0, 0, 0}         // GOP Index 1 - Ref List 0
                //{ 0, 0, 0, 0 }       // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = base1_idx;
                av1_rps->ref_dpb_index[ALT2] = base1_idx;
                av1_rps->ref_dpb_index[ALT] = base1_idx;

                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
            }
            else
            if (picture_index == 0) {
                //{ 1, 9, 8, 17},  // GOP Index 1 - Ref List 0
                //{ -1, -3, -7, 0 }   // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
                av1_rps->ref_dpb_index[LAST3] = lay4_idx;
                av1_rps->ref_dpb_index[GOLD] = base0_idx;
                av1_rps->ref_dpb_index[BWD] = lay3_idx;
                av1_rps->ref_dpb_index[ALT2] = lay2_1_idx;
                av1_rps->ref_dpb_index[ALT] = lay1_1_idx;
                gop_i = 1;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 2) {
                //{ 1, 3, 2, 11},  // GOP Index 3 - Ref List 0
               //{ -1, -5, -13, 0 }   // GOP Index 3 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay4_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay2_1_idx;
                av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 3;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 4) {
                //{ 1, 5, 4, 13},  // GOP Index 5 - Ref List 0
               //{ -1, -3, -11, 0 }   // GOP Index 5 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay4_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay3_idx;
                av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 5;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 6) {
                //{ 1, 3, 6, 7},  // GOP Index 7 - Ref List 0
               //{ -1, -9, 0, 0 }   // GOP Index 7 - Ref List 1
               av1_rps->ref_dpb_index[LAST] = lay3_idx;
               av1_rps->ref_dpb_index[LAST2] = lay2_1_idx;
               av1_rps->ref_dpb_index[LAST3] = lay4_idx;
               av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 7;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 8) {
                //{ 1, 9, 8, 17},  // GOP Index 9 - Ref List 0
                //{ -1, -3, -7, 0 }   // GOP Index 9 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay4_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay3_idx;
                av1_rps->ref_dpb_index[ALT2] = lay2_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 9;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 10) {
                //{ 1, 3, 2, 11},  // GOP Index 11 - Ref List 0
                //{ -1, -5, 0, 0 }   // GOP Index 11 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay4_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = lay2_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 11;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 12) {
                //{ 1, 5, 4, 13},  // GOP Index 13 - Ref List 0
                //{ -1, -3, 0, 0 }   // GOP Index 13 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay4_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = lay3_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 13;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 14) {
                //{ 1, 3, 6, 7},  // GOP Index 15 - Ref List 0
                //{ -1, 0, 0, 0 }   // GOP Index 15 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST3] = lay4_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_1_idx;
                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = base1_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 15;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else
                SVT_LOG("Error in GOp indexing\n");

            uint8_t keep_lay4_as_ref;
            if (pcs_ptr->scs_ptr->mrp_ctrls.referencing_scheme == 1)
                keep_lay4_as_ref = (picture_index == 0 || picture_index == 8 || picture_index == 14 || picture_index == 2 || picture_index == 4 || picture_index == 6 || picture_index == 10 || picture_index == 12);
            else
                keep_lay4_as_ref = (picture_index == 0 || picture_index == 8);
            if (keep_lay4_as_ref)
                av1_rps->refresh_frame_mask = 1 << (lay4_idx);
            else
                av1_rps->refresh_frame_mask = 0;
            break;

        default:
            SVT_ERROR("unexpected picture mini Gop number\n");
            break;
        }


        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {

            uint8_t no_show;
            if (pcs_ptr->scs_ptr->mrp_ctrls.referencing_scheme == 1)
                no_show = (pcs_ptr->is_used_as_reference_flag && picture_index != 0 && picture_index != 8 && picture_index != 14 && picture_index != 2 && picture_index != 4 && picture_index != 6 && picture_index != 10 && picture_index != 12);
            else
                no_show = (pcs_ptr->is_used_as_reference_flag && picture_index != 0 && picture_index != 8);
            if (no_show)
            {
                frm_hdr->show_frame = FALSE;
                pcs_ptr->has_show_existing = FALSE;
            } else {
                frm_hdr->show_frame = TRUE;
                pcs_ptr->has_show_existing = TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = lay3_idx;
                else if (picture_index == 2)
                    frm_hdr->show_existing_frame = lay2_1_idx;
                else if (picture_index == 4)
                    frm_hdr->show_existing_frame = lay3_idx;
                else if (picture_index == 6)
                    frm_hdr->show_existing_frame = lay1_1_idx;
                else if (picture_index == 8)
                    frm_hdr->show_existing_frame = lay3_idx;
                else if (picture_index == 10)
                    frm_hdr->show_existing_frame = lay2_1_idx;
                else if (picture_index == 12)
                    frm_hdr->show_existing_frame = lay3_idx;
                else if (picture_index == 14)
                    frm_hdr->show_existing_frame = base2_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", pcs_ptr->hierarchical_levels);
            }
        }

        //last pic in MiniGop: Base layer toggling
        //mini GOP toggling since last Key Frame.
        //a regular I keeps the toggling process and does not reset the toggle.  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
        //whoever needs a miniGOP Level toggling, this is the time
        if (((picture_index == (context_ptr->mini_gop_end_index[mini_gop_index]) &&
                        !pcs_ptr->is_overlay &&
                        scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS)) ||
                ((scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {

            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
        }
    } else if (pcs_ptr->hierarchical_levels == 5) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(pcs_ptr, context_ptr);
            return;
        }
        //DPB: Loc7|Loc6|Loc5|Loc4|Loc3|Loc2|Loc1|Loc0
        //Layer 0 : circular move 0-1
        //Layer 1 : circular move 2-3
        //Layer 2 : DPB Location 4
        //Layer 3 : DPB Location 5
        //Layer 4 : DPB Location 6
        //Layer 5 : DPB Location 7
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 1 : 0;   //the oldest L0 picture in the DPB
        const uint8_t  base2_idx = context_ptr->lay0_toggle == 0 ? 0 : 1;   //the newest L0 picture in the DPB

        const uint8_t  lay1_0_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 1 : LAY1_OFF + 0;   //the oldest L1 picture in the DPB
        const uint8_t  lay1_1_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 0 : LAY1_OFF + 1;   //the newest L1 picture in the DPB
        const uint8_t  lay2_idx = LAY2_OFF;   //the oldest L2 picture in the DPB
        const uint8_t  lay3_idx = LAY3_OFF;    //the newest L3 picture in the DPB
        const uint8_t  lay4_idx = LAY4_OFF;    //the newest L4 picture in the DPB
        const uint8_t  lay5_idx = LAY5_OFF;    //the newest L5 picture in the DPB
        switch (pcs_ptr->temporal_layer_index) {
        case 0:
            //{32, 64, 48, 0}, // GOP Index 0 - Ref List 0
            //{32, 0, 0, 0}       // GOP Index 0 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = base2_idx;
            av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[BWD] = base1_idx;
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];

            gop_i = 0;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;
            break;

        case 1:
            //{16, 32, 0, 0}, // GOP Index 16 - Ref List 0
            //{-16, 24, 20, 0} // GOP Index 16 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
            av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base2_idx;
            av1_rps->ref_dpb_index[ALT2] = lay2_idx;
            av1_rps->ref_dpb_index[ALT] = lay3_idx;

            gop_i = 16;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
            av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            av1_rps->refresh_frame_mask = 1 << (LAY1_OFF + context_ptr->lay1_toggle);

            break;

        case 2:

            if (picture_index == 7) {
                //{8, 16, 24, 0}, // GOP Index 8 - Ref List 0
                //{-8, -24, 12, 0 } // GOP Index 8 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = lay3_idx;
                gop_i = 8;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];
                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 23) {
                //{8, 16, 24, 0}    // GOP Index 24 - Ref List 0
                //{-8, 10, 9, 0} // GOP Index 24 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = lay4_idx;
                av1_rps->ref_dpb_index[ALT] = lay5_idx;
                gop_i = 24;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            av1_rps->refresh_frame_mask = 1 << (LAY2_OFF);

            break;

        case 3:

            if (picture_index == 3) {
                //{4, 8, 20, 0}, // GOP Index 4 - Ref List 0
                //{-4, -12, -28, 0} // GOP Index 4 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                av1_rps->ref_dpb_index[BWD] = lay2_idx;
                av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 4;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 11) {
                //{4, 8, 12, 28}, // GOP Index 12 - Ref List 0
                //{-4, -20, 5, 0} // GOP Index 12 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = lay5_idx;
                gop_i = 12;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 19) {
                //{4, 8, 20, 0}, // GOP Index 20 - Ref List 0
                //{-4, -12, 0, 0} // GOP Index 20 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                av1_rps->ref_dpb_index[BWD] = lay2_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 20;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 27) {
                //{4, 8, 12, 28}, // GOP Index 28 - Ref List 0
                //{-4, 5, 0, 0} // GOP Index 28 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_1_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;

                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = lay5_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 28;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else
                SVT_LOG("Error in GOp indexing\n");
            av1_rps->refresh_frame_mask = 1 << (LAY3_OFF);
            break;

        case 4:
            // update RPS for the overlay frame.
            if (picture_index == 1) {
                //{2, 4, 18, 0}, // GOP Index 2 - Ref List 0
                //{-2, -6, -14, 0} // GOP Index 2 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                av1_rps->ref_dpb_index[BWD] = lay3_idx;
                av1_rps->ref_dpb_index[ALT2] = lay2_idx;
                av1_rps->ref_dpb_index[ALT] = lay1_1_idx;
                gop_i = 2;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 5) {
                //{2, 4, 6, 22}, // GOP Index 6 - Ref List 0
                //{-2, -10, -26, 0} // GOP Index 6 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay2_idx;
                av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 6;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 9) {
                //{2, 4, 10, 26}, // GOP Index 10 - Ref List 0
                //{-2, -6, -22, 0} // GOP Index 10 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay3_idx;
                av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 10;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 13) {
                //{2, 4, 6, 14}, // GOP Index 14 - Ref List 0
                //{-2, -18, 3, 0} // GOP Index 14 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = lay2_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = lay5_idx;
                gop_i = 14;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 17) {
                //{2, 4, 18,  34}, // GOP Index 18 - Ref List 0
                //{-2, -6, -14, 0} // GOP Index 18 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay3_idx;
                av1_rps->ref_dpb_index[ALT2] = lay2_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 18;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 21) {
                //{2, 4, 6, 22}, // GOP Index 22 - Ref List 0
                //{-2, -10, 0, 0} // GOP Index 22 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_1_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = lay2_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 22;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 25) {
                //{2, 4, 10, 26}, // GOP Index 26 - Ref List 0
                //{-2, -6, 0, 0} // GOP Index 26 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_1_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = lay3_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 26;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 29) {
                //{2, 4, 6, 14}, // GOP Index 30 - Ref List 0
                //{-2, 3, 0, 0} // GOP Index 30 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = lay2_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_1_idx;
                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = lay5_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 30;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else
                SVT_LOG("Error in GOp indexing\n");
            av1_rps->refresh_frame_mask = 1 << (LAY4_OFF);
            break;

        case 5:
            // update RPS for the overlay frame.
            if (pcs_ptr->is_overlay) {
                //{ 0, 0, 0, 0}         // GOP Index 1 - Ref List 0
                //{ 0, 0, 0, 0 }       // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = base1_idx;
                av1_rps->ref_dpb_index[ALT2] = base1_idx;
                av1_rps->ref_dpb_index[ALT] = base1_idx;

                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, 0);
            } else {
                if (picture_index == 0) {
                    //{1, 17, 2, 0}, // GOP Index 1 - Ref List 0
                    //{-1, -3, -7, 0} // GOP Index 1 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = base1_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT] = lay2_idx;
                    gop_i = 1;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 2) {
                    //{1, 3, 2, 19}, // GOP Index 3 - Ref List 0
                    //{-1, -5, -13, 0} // GOP Index 3 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                    av1_rps->ref_dpb_index[BWD] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay2_idx;
                    av1_rps->ref_dpb_index[ALT] = lay1_1_idx;
                    gop_i = 3;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 4) {
                    //{1, 5, 2, 21}, // GOP Index 5 - Ref List 0
                    //{-1, -3, -11, 0} // GOP Index 5 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay2_idx;
                    av1_rps->ref_dpb_index[ALT] = lay1_1_idx;
                    gop_i = 5;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 6) {
                    //{1, 3, 2, 7}, // GOP Index 7 - Ref List 0
                    //{-1, -9, -25, 0} // GOP Index 7 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay2_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[ALT] = base2_idx;
                    gop_i = 7;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 8) {
                    //{1, 9, 2, 25}, // GOP Index 9 - Ref List 0
                    //{-1, -3, -7, 0} // GOP Index 9 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT] = lay1_1_idx;
                    gop_i = 9;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 10) {
                    //{1, 3, 2, 11}, // GOP Index 11 - Ref List 0
                    //{-1, -5, -21, 0} // GOP Index 11 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[ALT] = base2_idx;
                    gop_i = 11;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 12) {
                    //{1, 5, 2, 13}, // GOP Index 13 - Ref List 0
                    //{-1, -3, -19, 0} // GOP Index 13 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[ALT] = base2_idx;
                    gop_i = 13;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 14) {
                    //{1, 3, 2, 7}, // GOP Index 15 - Ref List 0
                    //{ -1, -17, 15, 0 } // GOP Index 15 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay2_idx;
                    av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                    av1_rps->ref_dpb_index[ALT2] = base2_idx;
                    av1_rps->ref_dpb_index[ALT] = base1_idx;
                    gop_i = 15;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 16) {
                    //{1, 17, 2, 33}, // GOP Index 17 - Ref List 0
                    //{-1, -3, -7, 0} // GOP Index 17 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT] = lay2_idx;
                    gop_i = 17;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 18) {
                    //{1, 3, 2, 19}, // GOP Index 19 - Ref List 0
                    //{-1, -5, -13, 0} // GOP Index 19 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay2_idx;
                    av1_rps->ref_dpb_index[ALT] = base2_idx;
                    gop_i = 19;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 20) {
                    //{1, 5, 2, 21}, // GOP Index 21 - Ref List 0
                    //{-1, -3, -11, 0} // GOP Index 21 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay2_idx;
                    av1_rps->ref_dpb_index[ALT] = base2_idx;
                    gop_i = 21;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 22) {
                    //{1, 3, 2, 7}, // GOP Index 23 - Ref List 0
                    //{-1, -9, 0, 0} // GOP Index 23 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay2_idx;
                    av1_rps->ref_dpb_index[ALT2] = base2_idx;
                    av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                    gop_i = 23;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
                }
                else if (picture_index == 24) {
                    //{1, 9, 2, 25}, // GOP Index 25 - Ref List 0
                    //{-1, -3, -7, 0} // GOP Index 25 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT] = base2_idx;
                    gop_i = 25;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 26) {
                    //{1, 3, 2, 11}, // GOP Index 27 - Ref List 0
                    //{-1, -5, 0, 0} // GOP Index 27 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT2] = base2_idx;
                    av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                    gop_i = 27;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
                }
                else if (picture_index == 28) {
                    //{1, 5, 2, 13}, // GOP Index 29 - Ref List 0
                    //{-1, -3, 0, 0} // GOP Index 29 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = base2_idx;
                    av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                    gop_i = 29;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
                }
                else if (picture_index == 30) {
                    //{1, 3, 2, 7}, // GOP Index 31 - Ref List 0
                    //{ -1, 15, 0, 0 } // GOP Index 31 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay5_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay2_idx;
                    av1_rps->ref_dpb_index[BWD] = base2_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                    gop_i = 31;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
                }
                else
                    SVT_LOG("Error in GOp indexing\n");
            }
            if (picture_index % 2  == 0 && picture_index < 31 )
                av1_rps->refresh_frame_mask = 1 << (lay5_idx);
            break;

        default:
            SVT_ERROR("unexpected picture mini Gop number\n");
            break;
        }


        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {
            if (pcs_ptr->is_used_as_reference_flag && (picture_index % 2 != 0)) {
                frm_hdr->show_frame = FALSE;
                pcs_ptr->has_show_existing = FALSE;
            } else {
                frm_hdr->show_frame = TRUE;
                pcs_ptr->has_show_existing = TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = lay4_idx;
                else if (picture_index == 2)
                    frm_hdr->show_existing_frame = lay3_idx;
                else if (picture_index == 4)
                    frm_hdr->show_existing_frame = lay4_idx;
                else if (picture_index == 6)
                    frm_hdr->show_existing_frame = lay2_idx;
                else if (picture_index == 8)
                    frm_hdr->show_existing_frame = lay4_idx;
                else if (picture_index == 10)
                    frm_hdr->show_existing_frame = lay3_idx;
                else if (picture_index == 12)
                    frm_hdr->show_existing_frame = lay4_idx;
                else if (picture_index == 14)
                    frm_hdr->show_existing_frame = lay1_1_idx;
                else if (picture_index == 16)
                    frm_hdr->show_existing_frame = lay4_idx;
                else if (picture_index == 18)
                    frm_hdr->show_existing_frame = lay3_idx;
                else if (picture_index == 20)
                    frm_hdr->show_existing_frame = lay4_idx;
                else if (picture_index == 22)
                    frm_hdr->show_existing_frame = lay2_idx;
                else if (picture_index == 24)
                    frm_hdr->show_existing_frame = lay4_idx;
                else if (picture_index == 26)
                    frm_hdr->show_existing_frame = lay3_idx;
                else if (picture_index == 28)
                    frm_hdr->show_existing_frame = lay4_idx;
                else if (picture_index == 30)
                    frm_hdr->show_existing_frame = base2_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", pcs_ptr->hierarchical_levels);
            }
        }

        //last pic in MiniGop: Base layer toggling
        //mini GOP toggling since last Key Frame.
        //a regular I keeps the toggling process and does not reset the toggle.  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
        //whoever needs a miniGOP Level toggling, this is the time
        if (((picture_index == (context_ptr->mini_gop_end_index[mini_gop_index]) &&
                        !pcs_ptr->is_overlay &&
                        scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS)) ||
                ((scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {
            //Layer0 toggle 0->1
            context_ptr->lay0_toggle = 1 - context_ptr->lay0_toggle;
            //Layer1 toggle 2->3
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
        }
    }else {
        SVT_ERROR("Not supported GOP structure!");
        exit(0);
    }

    if (frm_hdr->frame_type == S_FRAME) {
        set_sframe_rps(pcs_ptr, encode_context_ptr, context_ptr);
    }
}

static EbErrorType av1_generate_rps_info_from_user_config(
    PictureParentControlSet *picture_control_set_ptr, EncodeContext *encode_context_ptr) {
    Av1RpsNode *              av1_rps = &picture_control_set_ptr->av1_ref_signal;
    FrameHeader *             frm_hdr = &picture_control_set_ptr->frm_hdr;
    PredictionStructureEntry *pred_position_ptr =
        picture_control_set_ptr->pred_struct_ptr
            ->pred_struct_entry_ptr_array[picture_control_set_ptr->pred_struct_index];
    DPBInfo *dpb_list_ptr = &encode_context_ptr->dpb_list[0];
    int32_t  dpb_list_idx = 0;

    if (frm_hdr->frame_type == KEY_FRAME) {
        frm_hdr->show_frame                        = TRUE;
        picture_control_set_ptr->has_show_existing = FALSE;
        EB_MEMSET(dpb_list_ptr, 0, sizeof(DPBInfo) * REF_FRAMES);
        encode_context_ptr->display_picture_number = picture_control_set_ptr->picture_number;
        dpb_list_ptr[dpb_list_idx].is_used         = TRUE;
        dpb_list_ptr[dpb_list_idx].picture_number  = picture_control_set_ptr->picture_number;
        dpb_list_ptr[dpb_list_idx].is_displayed    = TRUE;
        dpb_list_ptr[dpb_list_idx].temporal_layer_index =
            picture_control_set_ptr->temporal_layer_index;

        // Construct dependent lists
        dpb_list_ptr[dpb_list_idx].dep_list0.list_count = 0;
        for (uint32_t dep_idx = 0; dep_idx < pred_position_ptr->dep_list0.list_count; ++dep_idx)
            if (pred_position_ptr->dep_list0.list[dep_idx] >= 0)
                dpb_list_ptr[dpb_list_idx]
                    .dep_list0.list[dpb_list_ptr[dpb_list_idx].dep_list0.list_count++] =
                    pred_position_ptr->dep_list0.list[dep_idx];
        dpb_list_ptr[dpb_list_idx].dep_list1.list_count = pred_position_ptr->dep_list1.list_count;
        for (uint32_t dep_idx = 0; dep_idx < pred_position_ptr->dep_list1.list_count; ++dep_idx)
            dpb_list_ptr[dpb_list_idx].dep_list1.list[dep_idx] =
                pred_position_ptr->dep_list1.list[dep_idx];
        dpb_list_ptr[dpb_list_idx].dep_list0_count =
            dpb_list_ptr[dpb_list_idx].dep_list0.list_count;
        dpb_list_ptr[dpb_list_idx].dep_list1_count =
            dpb_list_ptr[dpb_list_idx].dep_list1.list_count;
        dpb_list_ptr[dpb_list_idx].dep_count = dpb_list_ptr[dpb_list_idx].dep_list0_count +
            dpb_list_ptr[dpb_list_idx].dep_list1_count;
        return EB_ErrorNone;
    }

    // If there was an I-frame or Scene Change, then cleanup the Reference Queue's Dependent Counts
    if (picture_control_set_ptr->slice_type == I_SLICE) {
        encode_context_ptr->is_i_slice_in_last_mini_gop = TRUE;
        encode_context_ptr->i_slice_picture_number_in_last_mini_gop =
            picture_control_set_ptr->picture_number;
    }

    // Construct dpb index mapping for ref list0
    if ((picture_control_set_ptr->slice_type == P_SLICE) ||
        (picture_control_set_ptr->slice_type == B_SLICE)) {
        uint8_t ref_idx = 0;
        for (ref_idx = LAST; ref_idx < LAST + picture_control_set_ptr->ref_list0_count; ++ref_idx) {
            uint64_t ref_poc = picture_control_set_ptr->picture_number - pred_position_ptr->ref_list0.reference_list[ref_idx-LAST];
            if (picture_control_set_ptr->is_overlay) {
                ref_poc = picture_control_set_ptr->picture_number;
            }

            dpb_list_idx = 0;
            do {
                if (dpb_list_ptr[dpb_list_idx].is_used == TRUE &&
                    dpb_list_ptr[dpb_list_idx].picture_number == ref_poc) {
                    if (dpb_list_ptr[dpb_list_idx].temporal_layer_index >
                        picture_control_set_ptr->temporal_layer_index)
                        return EB_Corrupt_Frame;
                    break;
                }
            } while (++dpb_list_idx < REF_FRAMES);
            if (dpb_list_idx < REF_FRAMES) {
                av1_rps->ref_dpb_index[ref_idx] = dpb_list_idx;
                --dpb_list_ptr[dpb_list_idx].dep_count;
                if (dpb_list_ptr[dpb_list_idx].dep_count < 0) {
                    SVT_ERROR("dep_count error in dpb list0\n");
                    return EB_Corrupt_Frame;
                }
            } else {
                SVT_ERROR("can't find ref frame in dpb list0\n");
                return EB_Corrupt_Frame;
            }
        }
        for (; ref_idx <= GOLD; ++ref_idx)
            av1_rps->ref_dpb_index[ref_idx] = av1_rps->ref_dpb_index[LAST];

        // Construct dpb index mapping for ref list1
        if (picture_control_set_ptr->slice_type == B_SLICE) {
            for (ref_idx = BWD; ref_idx < BWD + picture_control_set_ptr->ref_list1_count;
                 ++ref_idx) {
                uint64_t ref_poc = picture_control_set_ptr->picture_number -
                    pred_position_ptr->ref_list1.reference_list[ref_idx - BWD];
                dpb_list_idx = 0;
                do {
                    if (dpb_list_ptr[dpb_list_idx].is_used == TRUE &&
                        dpb_list_ptr[dpb_list_idx].picture_number == ref_poc) {
                        if (dpb_list_ptr[dpb_list_idx].temporal_layer_index >
                            picture_control_set_ptr->temporal_layer_index)
                            return EB_Corrupt_Frame;
                        break;
                    }
                } while (++dpb_list_idx < REF_FRAMES);

                if (dpb_list_idx < REF_FRAMES) {
                    av1_rps->ref_dpb_index[ref_idx] = dpb_list_idx;
                    --dpb_list_ptr[dpb_list_idx].dep_count;
                    if (dpb_list_ptr[dpb_list_idx].dep_count < 0) {
                        SVT_ERROR("dep_count error in dpb list1\n");
                        return EB_Corrupt_Frame;
                    }
                } else {
                    SVT_ERROR("can't find ref frame in dpb list1\n");
                    return EB_Corrupt_Frame;
                }
            }
            for (; ref_idx <= ALT; ++ref_idx)
                av1_rps->ref_dpb_index[ref_idx] = av1_rps->ref_dpb_index[BWD];
        } else
            av1_rps->ref_dpb_index[BWD]     = av1_rps->ref_dpb_index[ALT2] =
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
    }

    // Release unused positions in dpb list
    dpb_list_idx = 0;
    do {
        if (dpb_list_ptr[dpb_list_idx].is_used == TRUE &&
            dpb_list_ptr[dpb_list_idx].is_displayed == TRUE &&
            dpb_list_ptr[dpb_list_idx].dep_count == 0) {
            dpb_list_ptr[dpb_list_idx].is_used      = FALSE;
            dpb_list_ptr[dpb_list_idx].is_displayed = FALSE;
        }
    } while (++dpb_list_idx < REF_FRAMES);

    // Insert current frame into dpb list
    if (picture_control_set_ptr->picture_number == encode_context_ptr->display_picture_number + 1 &&
        picture_control_set_ptr->is_used_as_reference_flag == 0) {
        frm_hdr->show_frame = TRUE;
        ++encode_context_ptr->display_picture_number;
        av1_rps->refresh_frame_mask = 0;
    } else {
        frm_hdr->show_frame = FALSE;
        dpb_list_idx        = 0;
        do {
            if (dpb_list_ptr[dpb_list_idx].is_used == FALSE)
                break;
        } while (++dpb_list_idx < REF_FRAMES);
        if (dpb_list_idx < REF_FRAMES) {
            av1_rps->refresh_frame_mask               = 1 << dpb_list_idx;
            dpb_list_ptr[dpb_list_idx].is_used        = TRUE;
            dpb_list_ptr[dpb_list_idx].picture_number = picture_control_set_ptr->picture_number;
            dpb_list_ptr[dpb_list_idx].is_alt_ref     = picture_control_set_ptr->is_alt_ref;
            dpb_list_ptr[dpb_list_idx].temporal_layer_index =
                picture_control_set_ptr->temporal_layer_index;
            if (picture_control_set_ptr->is_alt_ref)
                dpb_list_ptr[dpb_list_idx].is_displayed = TRUE;
            // Construct dependent lists
            dpb_list_ptr[dpb_list_idx].dep_list0.list_count = 0;
            for (uint32_t dep_idx = 0; dep_idx < pred_position_ptr->dep_list0.list_count; ++dep_idx)
                if (pred_position_ptr->dep_list0.list[dep_idx] >= 0)
                    dpb_list_ptr[dpb_list_idx]
                        .dep_list0.list[dpb_list_ptr[dpb_list_idx].dep_list0.list_count++] =
                        pred_position_ptr->dep_list0.list[dep_idx];
            dpb_list_ptr[dpb_list_idx].dep_list1.list_count =
                pred_position_ptr->dep_list1.list_count;
            for (uint32_t dep_idx = 0; dep_idx < pred_position_ptr->dep_list1.list_count; ++dep_idx)
                dpb_list_ptr[dpb_list_idx].dep_list1.list[dep_idx] =
                    pred_position_ptr->dep_list1.list[dep_idx];
            dpb_list_ptr[dpb_list_idx].dep_list0_count = (picture_control_set_ptr->is_alt_ref)
                ? dpb_list_ptr[dpb_list_idx].dep_list0.list_count + 1
                : dpb_list_ptr[dpb_list_idx].dep_list0.list_count;
            dpb_list_ptr[dpb_list_idx].dep_list1_count =
                dpb_list_ptr[dpb_list_idx].dep_list1.list_count;
            dpb_list_ptr[dpb_list_idx].dep_count = dpb_list_ptr[dpb_list_idx].dep_list0_count +
                dpb_list_ptr[dpb_list_idx].dep_list1_count;
        } else {
            SVT_ERROR("can't find unused dpb to hold current frame\n");
            return EB_Corrupt_Frame;
        }
    }

    // Calculate output flag
    picture_control_set_ptr->has_show_existing = FALSE;
    do {
        dpb_list_idx = 0;
        do {
            if (dpb_list_ptr[dpb_list_idx].is_used == TRUE &&
                dpb_list_ptr[dpb_list_idx].is_displayed == FALSE &&
                dpb_list_ptr[dpb_list_idx].picture_number ==
                    encode_context_ptr->display_picture_number + 1)
                break;
        } while (++dpb_list_idx < REF_FRAMES);
        if (dpb_list_idx < REF_FRAMES) {
            dpb_list_ptr[dpb_list_idx].is_displayed = TRUE;
            ++encode_context_ptr->display_picture_number;
            if (encode_context_ptr->display_picture_number ==
                picture_control_set_ptr->picture_number)
                frm_hdr->show_frame = TRUE;
            else {
                picture_control_set_ptr->has_show_existing = TRUE;
                frm_hdr->show_existing_frame               = dpb_list_idx;
                break;
            }
        } else
            break;
    } while (1);

    return EB_ErrorNone;
}

static EbErrorType av1_generate_minigop_rps_info_from_user_config(
    EncodeContext                 *encode_context_ptr,
    PictureDecisionContext        *context_ptr,
    uint32_t                       mini_gop_index
)
{
    PictureParentControlSet       *picture_control_set_ptr = NULL;
    if (encode_context_ptr->is_mini_gop_changed) {
        PredictionStructure          *next_pred_struct_ptr;
        PredictionStructureEntry     *next_base_layer_pred_position_ptr;
        uint32_t                      dependant_list_removed_entries, dep_idx = 0;
        int32_t                       dpb_list_idx = 0;
        DPBInfo                      *dpb_list_ptr = &encode_context_ptr->dpb_list[0];

        picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_start_index[mini_gop_index]]->object_ptr;
        while (dpb_list_idx < REF_FRAMES) {
            DPBInfo *referenceEntryPtr = &dpb_list_ptr[dpb_list_idx];

            if (referenceEntryPtr->picture_number == (picture_control_set_ptr->picture_number - 1)) {

                next_pred_struct_ptr = picture_control_set_ptr->pred_struct_ptr;

                next_base_layer_pred_position_ptr = next_pred_struct_ptr->pred_struct_entry_ptr_array[next_pred_struct_ptr->pred_struct_entry_count - 1];

                // Remove all positive entries from the dependant lists
                int32_t dependant_list_positive_entries = 0;
                for (dep_idx = 0; dep_idx < referenceEntryPtr->dep_list0.list_count; ++dep_idx) {
                    if (referenceEntryPtr->dep_list0.list[dep_idx] >= 0) {
                        dependant_list_positive_entries++;
                    }
                }
                referenceEntryPtr->dep_list0.list_count = referenceEntryPtr->dep_list0.list_count - dependant_list_positive_entries;

                dependant_list_positive_entries = 0;
                for (dep_idx = 0; dep_idx < referenceEntryPtr->dep_list1.list_count; ++dep_idx) {
                    if (referenceEntryPtr->dep_list1.list[dep_idx] >= 0) {
                        dependant_list_positive_entries++;
                    }
                }
                referenceEntryPtr->dep_list1.list_count = referenceEntryPtr->dep_list1.list_count - dependant_list_positive_entries;

                for (dep_idx = 0; dep_idx < next_base_layer_pred_position_ptr->dep_list0.list_count; ++dep_idx) {
                    if (next_base_layer_pred_position_ptr->dep_list0.list[dep_idx] >= 0) {
                        referenceEntryPtr->dep_list0.list[referenceEntryPtr->dep_list0.list_count++] = next_base_layer_pred_position_ptr->dep_list0.list[dep_idx];
                    }
                }

                for (dep_idx = 0; dep_idx < next_base_layer_pred_position_ptr->dep_list1.list_count; ++dep_idx) {
                    if (next_base_layer_pred_position_ptr->dep_list1.list[dep_idx] >= 0) {
                        referenceEntryPtr->dep_list1.list[referenceEntryPtr->dep_list1.list_count++] = next_base_layer_pred_position_ptr->dep_list1.list[dep_idx];
                    }
                }

                // Update the dependant count update
                dependant_list_removed_entries = referenceEntryPtr->dep_list0_count + referenceEntryPtr->dep_list1_count - referenceEntryPtr->dep_count;
                referenceEntryPtr->dep_list0_count = (referenceEntryPtr->is_alt_ref) ? referenceEntryPtr->dep_list0.list_count + 1 : referenceEntryPtr->dep_list0.list_count;
                referenceEntryPtr->dep_list1_count = referenceEntryPtr->dep_list1.list_count;
                referenceEntryPtr->dep_count = referenceEntryPtr->dep_list0_count + referenceEntryPtr->dep_list1_count - dependant_list_removed_entries;
            }
            else {
                // Modify Dependent List0
                for (dep_idx = 0; dep_idx < referenceEntryPtr->dep_list0.list_count; ++dep_idx) {
                    uint64_t depPoc = POC_CIRCULAR_ADD(
                    referenceEntryPtr->picture_number,
                    referenceEntryPtr->dep_list0.list[dep_idx]);

                    if (depPoc >= picture_control_set_ptr->picture_number && referenceEntryPtr->dep_list0.list[dep_idx]) {
                        referenceEntryPtr->dep_list0.list[dep_idx] = 0;

                        --referenceEntryPtr->dep_count;

                        if (referenceEntryPtr->dep_count < 0) {
                            return EB_Corrupt_Frame;
                        }
                    }
                }

                // Modify Dependent List1
                for (dep_idx = 0; dep_idx < referenceEntryPtr->dep_list1.list_count; ++dep_idx) {
                    uint64_t depPoc = POC_CIRCULAR_ADD(
                    referenceEntryPtr->picture_number,
                    referenceEntryPtr->dep_list1.list[dep_idx]);

                    if ((depPoc >= picture_control_set_ptr->picture_number) && referenceEntryPtr->dep_list1.list[dep_idx]) {
                        referenceEntryPtr->dep_list1.list[dep_idx] = 0;

                        --referenceEntryPtr->dep_count;

                        if (referenceEntryPtr->dep_count < 0) {
                            return EB_Corrupt_Frame;
                        }
                    }
                }
            }
            ++dpb_list_idx;
        }
    }

    if (encode_context_ptr->is_i_slice_in_last_mini_gop == TRUE)
    {
        encode_context_ptr->is_i_slice_in_last_mini_gop = FALSE;
        int32_t dpb_list_idx = 0;
        uint32_t dep_idx = 0;
        DPBInfo *dpb_list_ptr = &encode_context_ptr->dpb_list[0];
        while (dpb_list_idx < REF_FRAMES) {
            DPBInfo *referenceEntryPtr = &dpb_list_ptr[dpb_list_idx];

            if (referenceEntryPtr->picture_number == encode_context_ptr->i_slice_picture_number_in_last_mini_gop) {
                ++dpb_list_idx;
                continue;
            }

            // Modify Dependent List0
            for (dep_idx = 0; dep_idx < referenceEntryPtr->dep_list0.list_count; ++dep_idx) {
                uint64_t depPoc = POC_CIRCULAR_ADD(
                referenceEntryPtr->picture_number,
                referenceEntryPtr->dep_list0.list[dep_idx]);

                if (depPoc >= encode_context_ptr->i_slice_picture_number_in_last_mini_gop && referenceEntryPtr->dep_list0.list[dep_idx]) {
                    referenceEntryPtr->dep_list0.list[dep_idx] = 0;
                    --referenceEntryPtr->dep_count;
                    if (referenceEntryPtr->dep_count < 0) {
                        return EB_Corrupt_Frame;
                    }
                }
            }

            // Modify Dependent List1
            for (dep_idx = 0; dep_idx < referenceEntryPtr->dep_list1.list_count; ++dep_idx) {
                uint64_t depPoc = POC_CIRCULAR_ADD(
                referenceEntryPtr->picture_number,
                referenceEntryPtr->dep_list1.list[dep_idx]);

                if ((depPoc >= encode_context_ptr->i_slice_picture_number_in_last_mini_gop)
                && referenceEntryPtr->dep_list1.list[dep_idx]) {
                    referenceEntryPtr->dep_list1.list[dep_idx] = 0;
                    --referenceEntryPtr->dep_count;
                    if (referenceEntryPtr->dep_count < 0) {
                        return EB_Corrupt_Frame;
                    }
                }
            }
            ++dpb_list_idx;
        }
    }

    // Add 1 to the loop for the overlay picture. If the last picture is alt ref, increase the loop by 1 to add the overlay picture
    uint32_t has_overlay = ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_end_index[mini_gop_index]]->object_ptr)->is_alt_ref ? 1 : 0;
    picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_start_index[mini_gop_index]]->object_ptr;
    SequenceControlSet *scs_ptr = picture_control_set_ptr->scs_ptr;
    for (uint32_t decode_order = 0,pictureIndex = context_ptr->mini_gop_start_index[mini_gop_index]; pictureIndex <= context_ptr->mini_gop_end_index[mini_gop_index]+has_overlay; ++pictureIndex, ++decode_order) {
        if (has_overlay && pictureIndex == context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay) {
            picture_control_set_ptr = ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_end_index[mini_gop_index]]->object_ptr)->overlay_ppcs_ptr;
        }
        else if ((context_ptr->mini_gop_length[mini_gop_index] == picture_control_set_ptr->pred_struct_ptr->pred_struct_period)
        &&(context_ptr->mini_gop_idr_count[mini_gop_index] == 0)
        &&(scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS)){
            uint32_t pic_idx = context_ptr->mini_gop_start_index[mini_gop_index];
            do {
                picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pic_idx]->object_ptr;
                if (decode_order == picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[picture_control_set_ptr->pred_struct_index]->decode_order) {
                    break;
                }
            } while (++pic_idx <= context_ptr->mini_gop_end_index[mini_gop_index]);
            if (pic_idx > context_ptr->mini_gop_end_index[mini_gop_index]) {
                return EB_Corrupt_Frame;
            }
        }
        else {
            picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pictureIndex]->object_ptr;
        }
        EbErrorType ret = av1_generate_rps_info_from_user_config(picture_control_set_ptr, encode_context_ptr);
        if (ret != EB_ErrorNone) {
            return EB_Corrupt_Frame;
        }
    }
    return EB_ErrorNone;
}

static void av1_generate_rps_ref_poc_from_user_config(
    PictureParentControlSet *picture_control_set_ptr) {
    Av1RpsNode *              av1_rps = &picture_control_set_ptr->av1_ref_signal;
    FrameHeader *             frm_hdr = &picture_control_set_ptr->frm_hdr;
    PredictionStructureEntry *pred_position_ptr =
        picture_control_set_ptr->pred_struct_ptr
            ->pred_struct_entry_ptr_array[picture_control_set_ptr->pred_struct_index];

    frm_hdr->frame_type = picture_control_set_ptr->slice_type == I_SLICE
        ? picture_control_set_ptr->idr_flag ? KEY_FRAME : INTRA_ONLY_FRAME
        : INTER_FRAME;

    picture_control_set_ptr->intra_only = picture_control_set_ptr->slice_type == I_SLICE;

    if (picture_control_set_ptr->is_overlay)
        for (int ref_idx = LAST; ref_idx <= ALT; ++ref_idx)
            av1_rps->ref_poc_array[ref_idx] = picture_control_set_ptr->picture_number;
    else if (picture_control_set_ptr->slice_type == P_SLICE ||
             picture_control_set_ptr->slice_type == B_SLICE) {
        uint32_t ref_idx;
        for (ref_idx = LAST; ref_idx < LAST + pred_position_ptr->ref_list0.reference_list_count;
             ++ref_idx)
            av1_rps->ref_poc_array[ref_idx] = picture_control_set_ptr->picture_number -
                pred_position_ptr->ref_list0.reference_list[ref_idx - LAST];
        for (; ref_idx <= GOLD; ++ref_idx)
            av1_rps->ref_poc_array[ref_idx] = av1_rps->ref_poc_array[LAST];
        if (picture_control_set_ptr->slice_type == B_SLICE) {
            for (ref_idx = BWD; ref_idx < BWD + pred_position_ptr->ref_list1.reference_list_count;
                 ++ref_idx)
                av1_rps->ref_poc_array[ref_idx] = picture_control_set_ptr->picture_number -
                    pred_position_ptr->ref_list1.reference_list[ref_idx - BWD];
            for (; ref_idx <= ALT; ++ref_idx)
                av1_rps->ref_poc_array[ref_idx] = av1_rps->ref_poc_array[BWD];
        } else
            av1_rps->ref_poc_array[BWD]     = av1_rps->ref_poc_array[ALT2] =
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[LAST];
    }
}

/***************************************************************************************************
// Perform Required Picture Analysis Processing for the Overlay frame
***************************************************************************************************/
    void perform_simple_picture_analysis_for_overlay(PictureParentControlSet     *pcs_ptr) {
        EbPictureBufferDesc           *input_padded_picture_ptr;
        EbPictureBufferDesc           *input_picture_ptr;
        EbPaReferenceObject           *pa_ref_obj_;

        SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
        input_picture_ptr = pcs_ptr->enhanced_picture_ptr;
        pa_ref_obj_ = (EbPaReferenceObject*)pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
        input_padded_picture_ptr = (EbPictureBufferDesc*)pa_ref_obj_->input_padded_picture_ptr;

        // Pad pictures to multiple min cu size
        pad_picture_to_multiple_of_min_blk_size_dimensions(
            scs_ptr,
            input_picture_ptr);

        // Pre processing operations performed on the input picture
        picture_pre_processing_operations(
            pcs_ptr,
            scs_ptr);

        if (input_picture_ptr->color_format >= EB_YUV422) {
            // Jing: Do the conversion of 422/444=>420 here since it's multi-threaded kernel
            //       Reuse the Y, only add cb/cr in the newly created buffer desc
            //       NOTE: since denoise may change the src, so this part is after picture_pre_processing_operations()
            pcs_ptr->chroma_downsampled_picture_ptr->buffer_y = input_picture_ptr->buffer_y;
            down_sample_chroma(input_picture_ptr, pcs_ptr->chroma_downsampled_picture_ptr);
        }
        else
            pcs_ptr->chroma_downsampled_picture_ptr = input_picture_ptr;

        // R2R FIX: copying input_picture_ptr to input_padded_picture_ptr for motion_estimate_sb needs it
        {
            uint8_t *pa = input_padded_picture_ptr->buffer_y + input_padded_picture_ptr->origin_x +
                input_padded_picture_ptr->origin_y * input_padded_picture_ptr->stride_y;
            uint8_t *in = input_picture_ptr->buffer_y + input_picture_ptr->origin_x +
                input_picture_ptr->origin_y * input_picture_ptr->stride_y;
            for (uint32_t row = 0; row < input_picture_ptr->height; row++)
                svt_memcpy(pa + row * input_padded_picture_ptr->stride_y, in + row * input_picture_ptr->stride_y, sizeof(uint8_t) * input_picture_ptr->width);
        }

        // Pad input picture to complete border SBs
        pad_picture_to_multiple_of_sb_dimensions(
            input_padded_picture_ptr);
        // 1/4 & 1/16 input picture downsampling through filtering
        if (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
            downsample_filtering_input_picture(
                pcs_ptr,
                input_padded_picture_ptr,
                (EbPictureBufferDesc*)pa_ref_obj_->quarter_downsampled_picture_ptr,
                (EbPictureBufferDesc*)pa_ref_obj_->sixteenth_downsampled_picture_ptr);
        }
        else {
            downsample_decimation_input_picture(
                pcs_ptr,
                input_padded_picture_ptr,
                (EbPictureBufferDesc*)pa_ref_obj_->quarter_downsampled_picture_ptr,
                (EbPictureBufferDesc*)pa_ref_obj_->sixteenth_downsampled_picture_ptr);
        }
        // Gathering statistics of input picture, including Variance Calculation, Histogram Bins
        gathering_picture_statistics(
            scs_ptr,
            pcs_ptr,
            input_padded_picture_ptr,
            pa_ref_obj_->sixteenth_downsampled_picture_ptr);

        pcs_ptr->sc_class0 = pcs_ptr->alt_ref_ppcs_ptr->sc_class0;
        pcs_ptr->sc_class1 = pcs_ptr->alt_ref_ppcs_ptr->sc_class1;
        pcs_ptr->sc_class2 = pcs_ptr->alt_ref_ppcs_ptr->sc_class2;
    }
/***************************************************************************************************
 * Initialize the overlay frame
***************************************************************************************************/
void initialize_overlay_frame(PictureParentControlSet     *pcs_ptr) {
    pcs_ptr->scene_change_flag = FALSE;
    pcs_ptr->cra_flag = FALSE;
    pcs_ptr->idr_flag = FALSE;
    pcs_ptr->last_idr_picture = pcs_ptr->alt_ref_ppcs_ptr->last_idr_picture;
    pcs_ptr->pred_structure = pcs_ptr->alt_ref_ppcs_ptr->pred_structure;
    pcs_ptr->pred_struct_ptr = pcs_ptr->alt_ref_ppcs_ptr->pred_struct_ptr;
    pcs_ptr->hierarchical_levels = pcs_ptr->alt_ref_ppcs_ptr->hierarchical_levels;
    pcs_ptr->hierarchical_layers_diff = 0;
    pcs_ptr->init_pred_struct_position_flag = FALSE;
    pcs_ptr->pre_assignment_buffer_count = pcs_ptr->alt_ref_ppcs_ptr->pre_assignment_buffer_count;

    perform_simple_picture_analysis_for_overlay(pcs_ptr);
 }

/*
  ret number of past picture(not including current) in mg buffer.

*/
int32_t avail_past_pictures(PictureParentControlSet**buf, uint32_t buf_size, uint64_t input_pic)
{
    //buffer has at least curr picture
    int32_t tot_past = 0;
    for (uint32_t pic = 0; pic < buf_size; pic++) {
        if (buf[pic]->picture_number < input_pic) {
            tot_past++;
        }
    }
    return tot_past;
}


/*
  searches a picture in a given pcs buffer
*/
int32_t search_this_pic(PictureParentControlSet**buf, uint32_t buf_size, uint64_t input_pic)
{
    int32_t index = -1;
    for (uint32_t pic = 0; pic < buf_size; pic++) {
        if (buf[pic]->picture_number == input_pic) {
            index = (int32_t)pic;
            break;
        }
    }
    return index;
}
/*
  Tells if an Intra picture should be delayed to get next mini-gop
*/
Bool is_delayed_intra(PictureParentControlSet *pcs) {


    if ((pcs->idr_flag || pcs->cra_flag) && pcs->pred_structure == SVT_AV1_PRED_RANDOM_ACCESS) {
        if (pcs->scs_ptr->static_config.intra_period_length == 0 || pcs->end_of_sequence_flag)
            return 0;
        else if (pcs->idr_flag || (pcs->cra_flag && pcs->pre_assignment_buffer_count < pcs->pred_struct_ptr->pred_struct_period))
            return 1;
        else
            return 0;
    }
    else
        return 0;
}
/*
  Performs first pass in ME process
*/
void process_first_pass_frame(
    SequenceControlSet      *scs_ptr, PictureParentControlSet *pcs_ptr, PictureDecisionContext  *context_ptr) {
    int16_t seg_idx;

    // Initialize Segments
    pcs_ptr->first_pass_seg_column_count = (uint8_t)(scs_ptr->me_segment_column_count_array[0]);
    pcs_ptr->first_pass_seg_row_count = (uint8_t)(scs_ptr->me_segment_row_count_array[0]);

    pcs_ptr->first_pass_seg_total_count = (uint16_t)(pcs_ptr->first_pass_seg_column_count  * pcs_ptr->first_pass_seg_row_count);
    pcs_ptr->first_pass_seg_acc = 0;
    first_pass_signal_derivation_multi_processes(scs_ptr, pcs_ptr);
    if (pcs_ptr->me_data_wrapper_ptr == NULL) {
        EbObjectWrapper               *me_wrapper;
        svt_get_empty_object(context_ptr->me_fifo_ptr, &me_wrapper);
        pcs_ptr->me_data_wrapper_ptr = me_wrapper;
        pcs_ptr->pa_me_data = (MotionEstimationData *)me_wrapper->object_ptr;
    }

    for (seg_idx = 0; seg_idx < pcs_ptr->first_pass_seg_total_count; ++seg_idx) {

        EbObjectWrapper               *out_results_wrapper_ptr;
        PictureDecisionResults        *out_results_ptr;
        svt_get_empty_object(
            context_ptr->picture_decision_results_output_fifo_ptr,
            &out_results_wrapper_ptr);
        out_results_ptr = (PictureDecisionResults*)out_results_wrapper_ptr->object_ptr;
        out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
        out_results_ptr->segment_index = seg_idx;
        out_results_ptr->task_type = TASK_FIRST_PASS_ME;
        svt_post_full_object(out_results_wrapper_ptr);
    }

    svt_block_on_semaphore(pcs_ptr->first_pass_done_semaphore);
    svt_release_object(pcs_ptr->me_data_wrapper_ptr);
    pcs_ptr->me_data_wrapper_ptr = (EbObjectWrapper *)NULL;
    pcs_ptr->pa_me_data = NULL;
}
/*
store this input  picture to be used for TF-ing of upcoming base
increment live count of the required ressources to be used by TF of upcoming base.
will be released once TF is done
*/
void low_delay_store_tf_pictures(
    SequenceControlSet      *scs,
    PictureParentControlSet *pcs,
    PictureDecisionContext  *ctx)
{
    const uint32_t mg_size = 1 << (scs->static_config.hierarchical_levels);
    const uint32_t tot_past = scs->tf_params_per_type[1].max_num_past_pics;
    if (pcs->temporal_layer_index != 0 && pcs->pic_index + 1 + tot_past >= mg_size)
    {
        //printf("Store:%lld \n", pcs->picture_number);
        //store this picture to be used for TF-ing upcoming base
        ctx->tf_pic_array[ctx->tf_pic_arr_cnt++] = pcs;

        //increment live count of these ressources to be used by TF of upcoming base. will be released once TF is done.
        svt_object_inc_live_count(pcs->p_pcs_wrapper_ptr, 1);
        svt_object_inc_live_count(pcs->input_picture_wrapper_ptr, 1);
        svt_object_inc_live_count(pcs->pa_reference_picture_wrapper_ptr, 1);
        if (pcs->eb_y8b_wrapper_ptr)
            svt_object_inc_live_count(pcs->eb_y8b_wrapper_ptr, 1);
    }
}
/*
 TF is done, release ressources and reset the tf picture buffer.
*/
void low_delay_release_tf_pictures(
    PictureDecisionContext  *ctx)
{
    for (uint32_t pic_it = 0; pic_it < ctx->tf_pic_arr_cnt; pic_it++) {

        PictureParentControlSet *past_pcs = ctx->tf_pic_array[pic_it];
        //printf("                   Release:%lld \n", past_pcs->picture_number);

        svt_release_object(past_pcs->input_picture_wrapper_ptr);

        if (past_pcs->eb_y8b_wrapper_ptr)
            svt_release_object(past_pcs->eb_y8b_wrapper_ptr);

        svt_release_object(past_pcs->pa_reference_picture_wrapper_ptr);
        //ppcs should be the last one to release
        svt_release_object(past_pcs->p_pcs_wrapper_ptr);
    }

    memset(ctx->tf_pic_array, 0, ctx->tf_pic_arr_cnt * sizeof(PictureParentControlSet *));
    ctx->tf_pic_arr_cnt = 0;
}

/*
  Performs Motion Compensated Temporal Filtering in ME process
*/
void mctf_frame(
    SequenceControlSet      *scs_ptr,
    PictureParentControlSet *pcs_ptr,
    PictureDecisionContext  *context_ptr,
    uint32_t               out_stride_diff64
)
{
    if (scs_ptr->static_config.pred_structure != SVT_AV1_PRED_RANDOM_ACCESS &&
        scs_ptr->tf_params_per_type[1].enabled)
        low_delay_store_tf_pictures(
            scs_ptr,
            pcs_ptr,
            context_ptr);
    if (pcs_ptr->tf_ctrls.enabled) {
        derive_tf_window_params(
            scs_ptr,
            scs_ptr->encode_context_ptr,
            pcs_ptr,
            context_ptr,
            out_stride_diff64);
        pcs_ptr->temp_filt_prep_done = 0;


        pcs_ptr->tf_tot_horz_blks = pcs_ptr->tf_tot_vert_blks = 0;

        // Start Filtering in ME processes
        {
            int16_t seg_idx;

            // Initialize Segments
            pcs_ptr->tf_segments_column_count = scs_ptr->tf_segment_column_count;
            pcs_ptr->tf_segments_row_count = scs_ptr->tf_segment_row_count;
            pcs_ptr->tf_segments_total_count = (uint16_t)(pcs_ptr->tf_segments_column_count  * pcs_ptr->tf_segments_row_count);
            pcs_ptr->temp_filt_seg_acc = 0;
            for (seg_idx = 0; seg_idx < pcs_ptr->tf_segments_total_count; ++seg_idx) {

                EbObjectWrapper               *out_results_wrapper_ptr;
                PictureDecisionResults        *out_results_ptr;

                svt_get_empty_object(
                    context_ptr->picture_decision_results_output_fifo_ptr,
                    &out_results_wrapper_ptr);
                out_results_ptr = (PictureDecisionResults*)out_results_wrapper_ptr->object_ptr;
                out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
                out_results_ptr->segment_index = seg_idx;
                out_results_ptr->task_type = 1;
                svt_post_full_object(out_results_wrapper_ptr);
            }

            svt_block_on_semaphore(pcs_ptr->temp_filt_done_semaphore);
        }


        if (pcs_ptr->tf_tot_horz_blks > pcs_ptr->tf_tot_vert_blks * 6 / 4){
            context_ptr->tf_motion_direction = 0;
        }
        else  if (pcs_ptr->tf_tot_vert_blks > pcs_ptr->tf_tot_horz_blks * 6 / 4) {
            context_ptr->tf_motion_direction = 1;
        }
        else {
            context_ptr->tf_motion_direction = -1;
        }

    }
    else
        pcs_ptr->temporal_filtering_on = FALSE; // set temporal filtering flag OFF for current picture

    pcs_ptr->is_noise_level = (context_ptr->last_i_noise_levels_log1p_fp16[0] >= VQ_NOISE_LVL_TH);

    if (scs_ptr->static_config.pred_structure != SVT_AV1_PRED_RANDOM_ACCESS &&
        scs_ptr->tf_params_per_type[1].enabled&&
        pcs_ptr->temporal_layer_index == 0)
        low_delay_release_tf_pictures(context_ptr);
}

/* this function sets up ME refs for a regular pic*/
void tpl_regular_setup_me_refs(
    PictureParentControlSet         *base_pcs,
    PictureParentControlSet         *cur_pcs)
{

    for (uint8_t list_index = REF_LIST_0; list_index < TOTAL_NUM_OF_REF_LISTS; list_index++) {

        uint8_t  ref_list_count = (list_index == REF_LIST_0) ?
                cur_pcs->ref_list0_count_try :
                cur_pcs->ref_list1_count_try;

        if (list_index == REF_LIST_0)
            cur_pcs->tpl_data.tpl_ref0_count = ref_list_count;
        else
            cur_pcs->tpl_data.tpl_ref1_count = ref_list_count;


        for (uint8_t ref_idx = 0; ref_idx < ref_list_count; ref_idx++) {

            uint64_t ref_poc = cur_pcs->ref_pic_poc_array[list_index][ref_idx];

            cur_pcs->tpl_data.ref_tpl_group_idx[list_index][ref_idx] = -1;
            for (uint32_t j = 0; j < base_pcs->tpl_group_size; j++) {
                if (ref_poc == base_pcs->tpl_group[j]->picture_number) {
                    cur_pcs->tpl_data.ref_in_slide_window[list_index][ref_idx] = TRUE;
                    cur_pcs->tpl_data.ref_tpl_group_idx[list_index][ref_idx] = j;
                    break;
                }
            }

            EbPaReferenceObject *ref_obj = (EbPaReferenceObject *)cur_pcs->ref_pa_pic_ptr_array[list_index][ref_idx]->object_ptr;

            cur_pcs->tpl_data.tpl_ref_ds_ptr_array[list_index][ref_idx].picture_number = ref_obj->picture_number;
            cur_pcs->tpl_data.tpl_ref_ds_ptr_array[list_index][ref_idx].picture_ptr = ref_obj->input_padded_picture_ptr;
            //not needed for TPL but could be linked.
            cur_pcs->tpl_data.tpl_ref_ds_ptr_array[list_index][ref_idx].sixteenth_picture_ptr = NULL;
            cur_pcs->tpl_data.tpl_ref_ds_ptr_array[list_index][ref_idx].quarter_picture_ptr = NULL;
        }
    }


}
/*
  prepare TPL data fields
*/
void tpl_prep_info(PictureParentControlSet    *pcs) {

     for (uint32_t pic_i = 0; pic_i < pcs->tpl_group_size; ++pic_i) {

         PictureParentControlSet* pcs_tpl = pcs->tpl_group[pic_i];

         {

             pcs_tpl->tpl_data.tpl_ref0_count = 0;
             pcs_tpl->tpl_data.tpl_ref1_count = 0;
             EB_MEMSET(pcs_tpl->tpl_data.ref_in_slide_window, 0, MAX_NUM_OF_REF_PIC_LIST*REF_LIST_MAX_DEPTH * sizeof(Bool));
             EB_MEMSET(pcs_tpl->tpl_data.tpl_ref_ds_ptr_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbDownScaledBufDescPtrArray));
             EB_MEMSET(pcs_tpl->tpl_data.tpl_ref_ds_ptr_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbDownScaledBufDescPtrArray));

             pcs_tpl->tpl_data.tpl_slice_type = pcs_tpl->slice_type;
             pcs_tpl->tpl_data.tpl_temporal_layer_index = pcs_tpl->temporal_layer_index;
             pcs_tpl->tpl_data.is_used_as_reference_flag = pcs_tpl->is_used_as_reference_flag;
             pcs_tpl->tpl_data.tpl_decode_order = pcs_tpl->decode_order;

             pcs_tpl->tpl_data.base_pcs = pcs;

             if (pcs_tpl->tpl_data.tpl_slice_type != I_SLICE) {
                tpl_regular_setup_me_refs(
                    pcs,
                    pcs_tpl);

             }
         }
     }
}

void send_picture_out(
    SequenceControlSet      *scs,
    PictureParentControlSet *pcs,
    PictureDecisionContext  *ctx)
{
    EbObjectWrapper               *me_wrapper;
    EbObjectWrapper               *out_results_wrapper;


    //every picture enherits latest motion direction from TF
    pcs->tf_motion_direction = ctx->tf_motion_direction;




        if (pcs->is_used_as_reference_flag) {
            EbObjectWrapper* reference_picture_wrapper;
            // Get Empty Reference Picture Object
            svt_get_empty_object(
                scs->encode_context_ptr->reference_picture_pool_fifo_ptr,
                &reference_picture_wrapper);
            pcs->reference_picture_wrapper_ptr = reference_picture_wrapper;
            // reset reference object in case of its members are altered by superres tool
            EbReferenceObject* ref =
                (EbReferenceObject*)reference_picture_wrapper->object_ptr;
            svt_reference_object_reset(ref, scs);
            // Give the new Reference a nominal live_count of 1
            svt_object_inc_live_count(pcs->reference_picture_wrapper_ptr, 1);
#if SRM_REPORT
    pcs->reference_picture_wrapper_ptr->pic_number= pcs->picture_number;
#endif
    }else {
            pcs->reference_picture_wrapper_ptr = NULL;
        }

        //get a new ME data buffer
        if (pcs->me_data_wrapper_ptr == NULL) {
            svt_get_empty_object(ctx->me_fifo_ptr, &me_wrapper);
            pcs->me_data_wrapper_ptr = me_wrapper;
            pcs->pa_me_data = (MotionEstimationData *)me_wrapper->object_ptr;
            //printf("[%ld]: Got me data [NORMAL] %p\n", pcs->picture_number, pcs->pa_me_data);
        }


        MrpCtrls* mrp_ctrl = &(scs->mrp_ctrls);
        uint8_t ref_count_used_list0 =
            MAX(mrp_ctrl->sc_base_ref_list0_count,
                MAX(mrp_ctrl->base_ref_list0_count,
                    MAX(mrp_ctrl->sc_non_base_ref_list0_count, mrp_ctrl->non_base_ref_list0_count)));

        uint8_t ref_count_used_list1 =
            MAX(mrp_ctrl->sc_base_ref_list1_count,
                MAX(mrp_ctrl->base_ref_list1_count,
                    MAX(mrp_ctrl->sc_non_base_ref_list1_count, mrp_ctrl->non_base_ref_list1_count)));

        uint8_t max_ref_to_alloc, max_cand_to_alloc;

        get_max_allocated_me_refs(ref_count_used_list0, ref_count_used_list1, &max_ref_to_alloc, &max_cand_to_alloc);

        pcs->pa_me_data->max_cand = max_cand_to_alloc;
        pcs->pa_me_data->max_refs = max_ref_to_alloc;
        pcs->pa_me_data->max_l0 = ref_count_used_list0;

    //****************************************************
    // Picture resizing for super-res tool
    //****************************************************

    // Scale picture if super-res is used
    // Handle SUPERRES_FIXED and SUPERRES_RANDOM modes here.
    // SUPERRES_QTHRESH and SUPERRES_AUTO modes are handled in rate control process because these modes depend on qindex
    if (scs->static_config.pass == ENC_SINGLE_PASS) {
        if (scs->static_config.resize_mode > RESIZE_NONE ||
            scs->static_config.superres_mode == SUPERRES_FIXED ||
            scs->static_config.superres_mode == SUPERRES_RANDOM) {
            init_resize_picture(scs, pcs);
        }
    }

    uint8_t gm_level = derive_gm_level(pcs);
    set_gm_controls(pcs, gm_level);

    for (uint32_t segment_index = 0; segment_index < pcs->me_segments_total_count; ++segment_index) {
        // Get Empty Results Object
        svt_get_empty_object(
            ctx->picture_decision_results_output_fifo_ptr,
            &out_results_wrapper);

        PictureDecisionResults* out_results = (PictureDecisionResults*)out_results_wrapper->object_ptr;
        out_results->pcs_wrapper_ptr = pcs->p_pcs_wrapper_ptr;
        out_results->segment_index = segment_index;
        out_results->task_type = TASK_PAME;
        //Post the Full Results Object
        svt_post_full_object(out_results_wrapper);
    }

}
/***************************************************************************************************
* Store the pcs pointers in the gf group, set the gf_interval and gf_update_due
***************************************************************************************************/
void store_gf_group(
    PictureParentControlSet *pcs,
    PictureDecisionContext  *ctx,
    uint32_t                 mg_size) {
    if (pcs->slice_type == I_SLICE || (!is_delayed_intra(pcs) && pcs->temporal_layer_index == 0) || pcs->slice_type == P_SLICE) {
        if (is_delayed_intra(pcs)) {
            pcs->gf_group[0] = (void*)pcs;
            EB_MEMCPY(&pcs->gf_group[1], ctx->mg_pictures_array, mg_size * sizeof(PictureParentControlSet*));
            pcs->gf_interval = 1 + mg_size;
        }
        else {
            if (pcs->slice_type == P_SLICE && mg_size > 0 && ctx->mg_pictures_array[mg_size - 1]->idr_flag)
                mg_size = MAX(0, (int)mg_size - 1);
            EB_MEMCPY(&pcs->gf_group[0], ctx->mg_pictures_array, mg_size * sizeof(PictureParentControlSet*));
            pcs->gf_interval = mg_size;
        }

        if (pcs->slice_type == I_SLICE && pcs->end_of_sequence_flag) {
            pcs->gf_interval = 1;
            pcs->gf_group[0] = (void*)pcs;
        }

        for (int pic_i = 0; pic_i < pcs->gf_interval; ++pic_i) {
            if (pcs->gf_group[pic_i]->slice_type == I_SLICE || (!is_delayed_intra(pcs) && pcs->gf_group[pic_i]->temporal_layer_index == 0) || pcs->gf_group[pic_i]->slice_type == P_SLICE)
                pcs->gf_group[pic_i]->gf_update_due = 1;
            else
                pcs->gf_group[pic_i]->gf_update_due = 0;

            // For P picture that come after I, we need to set the gf_group pictures. It is used later in RC
            if (pcs->slice_type == I_SLICE && pcs->gf_group[pic_i]->slice_type == P_SLICE && pcs->picture_number < pcs->gf_group[pic_i]->picture_number) {
                pcs->gf_group[pic_i]->gf_interval = pcs->gf_interval - 1;
                EB_MEMCPY(&pcs->gf_group[pic_i]->gf_group[0], &ctx->mg_pictures_array[1], pcs->gf_group[pic_i]->gf_interval * sizeof(PictureParentControlSet*));
                pcs->gf_group[pic_i]->gf_update_due = 0;
            }
        }
    }
}

#if LAD_MG_PRINT
/* prints content of pre-assignment buffer */
void print_pre_ass_buffer(EncodeContext *ctx, PictureParentControlSet *pcs_ptr, uint8_t log)
{
    if (log) {

        if (ctx->pre_assignment_buffer_intra_count > 0)
            SVT_LOG("PRE-ASSIGN INTRA   (%i pictures)  POC:%lld \n", ctx->pre_assignment_buffer_count, pcs_ptr->picture_number);
        if (ctx->pre_assignment_buffer_count == (uint32_t)(1 << pcs_ptr->scs_ptr->static_config.hierarchical_levels))
            SVT_LOG("PRE-ASSIGN COMPLETE   (%i pictures)  POC:%lld \n", ctx->pre_assignment_buffer_count, pcs_ptr->picture_number);
        if (ctx->pre_assignment_buffer_eos_flag == 1)
            SVT_LOG("PRE-ASSIGN EOS   (%i pictures)  POC:%lld \n", ctx->pre_assignment_buffer_count, pcs_ptr->picture_number);
        if (pcs_ptr->pred_structure == SVT_AV1_PRED_LOW_DELAY_P)
            SVT_LOG("PRE-ASSIGN LDP   (%i pictures)  POC:%lld \n", ctx->pre_assignment_buffer_count, pcs_ptr->picture_number);

        SVT_LOG("\n Pre-Assign(%i):  ", ctx->pre_assignment_buffer_count);
        for (uint32_t pic = 0; pic < ctx->pre_assignment_buffer_count; pic++) {
            PictureParentControlSet *pcs = (PictureParentControlSet*)ctx->pre_assignment_buffer[pic]->object_ptr;
            SVT_LOG("%ld ", pcs->picture_number);
        }
        SVT_LOG("\n");
    }
}
#endif

void pack_highbd_pic(const EbPictureBufferDesc *pic_ptr, uint16_t *buffer_16bit[3], uint32_t ss_x,
    uint32_t ss_y, Bool include_padding);

EbErrorType derive_tf_window_params(
    SequenceControlSet *scs_ptr,
    EncodeContext *encode_context_ptr,
    PictureParentControlSet *pcs_ptr,
    PictureDecisionContext *context_ptr,
    uint32_t out_stride_diff64) {
    PictureParentControlSet * picture_control_set_ptr_central = pcs_ptr;
    EbPictureBufferDesc * central_picture_ptr = picture_control_set_ptr_central->enhanced_picture_ptr;
    uint32_t encoder_bit_depth = picture_control_set_ptr_central->scs_ptr->static_config.encoder_bit_depth;
    Bool is_highbd = (encoder_bit_depth == 8) ? (uint8_t)FALSE : (uint8_t)TRUE;

    // chroma subsampling
    uint32_t ss_x = picture_control_set_ptr_central->scs_ptr->subsampling_x;
    uint32_t ss_y = picture_control_set_ptr_central->scs_ptr->subsampling_y;
    int32_t *noise_levels_log1p_fp16 = &(picture_control_set_ptr_central->noise_levels_log1p_fp16[0]);
    int32_t noise_level_fp16;
    double *noise_levels = &(picture_control_set_ptr_central->noise_levels[0]);


    uint8_t do_noise_est = pcs_ptr->tf_ctrls.use_intra_for_noise_est ? 0 : 1;
    if (picture_control_set_ptr_central->slice_type == I_SLICE)
        do_noise_est = 1;
    // allocate 16 bit buffer
    if (is_highbd) {
        EB_MALLOC_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_Y],
            central_picture_ptr->luma_size);
        if (pcs_ptr->tf_ctrls.do_chroma) {
            EB_MALLOC_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_U],
                central_picture_ptr->chroma_size);
            EB_MALLOC_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_V],
                central_picture_ptr->chroma_size);
        }

        // pack byte buffers to 16 bit buffer
        pack_highbd_pic(central_picture_ptr,
            picture_control_set_ptr_central->altref_buffer_highbd,
            ss_x,
            ss_y,
            TRUE);
        // Estimate source noise level
        uint16_t *altref_buffer_highbd_start[COLOR_CHANNELS];
        altref_buffer_highbd_start[C_Y] =
            picture_control_set_ptr_central->altref_buffer_highbd[C_Y] +
            central_picture_ptr->origin_y * central_picture_ptr->stride_y +
            central_picture_ptr->origin_x;

        if (pcs_ptr->tf_ctrls.do_chroma) {
            altref_buffer_highbd_start[C_U] =
                picture_control_set_ptr_central->altref_buffer_highbd[C_U] +
                (central_picture_ptr->origin_y >> ss_y) * central_picture_ptr->stride_bit_inc_cb +
                (central_picture_ptr->origin_x >> ss_x);

            altref_buffer_highbd_start[C_V] =
                picture_control_set_ptr_central->altref_buffer_highbd[C_V] +
                (central_picture_ptr->origin_y >> ss_y) * central_picture_ptr->stride_bit_inc_cr +
                (central_picture_ptr->origin_x >> ss_x);
        }
        else {
            altref_buffer_highbd_start[C_U] = NOT_USED_VALUE;
            altref_buffer_highbd_start[C_V] = NOT_USED_VALUE;
        }

    if (pcs_ptr->tf_ctrls.use_fixed_point || pcs_ptr->tf_ctrls.use_medium_filter) {
    if (do_noise_est)
    {
            noise_level_fp16 = svt_estimate_noise_highbd_fp16(altref_buffer_highbd_start[C_Y], // Y only
                                                              central_picture_ptr->width,
                                                              central_picture_ptr->height,
                                                              central_picture_ptr->stride_y,
                                                              encoder_bit_depth);
            noise_levels_log1p_fp16[C_Y] = noise_log1p_fp16(noise_level_fp16);
    }
    } else
            if (do_noise_est)
            noise_levels[0] = svt_estimate_noise_highbd(altref_buffer_highbd_start[C_Y], // Y only
            central_picture_ptr->width,
            central_picture_ptr->height,
            central_picture_ptr->stride_y,
            encoder_bit_depth);

        if (pcs_ptr->tf_ctrls.do_chroma) {
        if (pcs_ptr->tf_ctrls.use_fixed_point || pcs_ptr->tf_ctrls.use_medium_filter) {
            noise_level_fp16 = svt_estimate_noise_highbd_fp16(altref_buffer_highbd_start[C_U], // U only
                                                             (central_picture_ptr->width >> 1),
                                                             (central_picture_ptr->height >> 1),
                                                             central_picture_ptr->stride_cb,
                                                             encoder_bit_depth);
            noise_levels_log1p_fp16[C_U] = noise_log1p_fp16(noise_level_fp16);

            noise_level_fp16 = svt_estimate_noise_highbd_fp16(altref_buffer_highbd_start[C_V], // V only
                                                             (central_picture_ptr->width >> 1),
                                                             (central_picture_ptr->height >> 1),
                                                             central_picture_ptr->stride_cb,
                                                             encoder_bit_depth);
            noise_levels_log1p_fp16[C_V] = noise_log1p_fp16(noise_level_fp16);
        } else {
        noise_levels[1] = svt_estimate_noise_highbd(altref_buffer_highbd_start[C_U], // U only
            (central_picture_ptr->width >> 1),
            (central_picture_ptr->height >> 1),
            central_picture_ptr->stride_cb,
            encoder_bit_depth);

        noise_levels[2] = svt_estimate_noise_highbd(altref_buffer_highbd_start[C_V], // V only
            (central_picture_ptr->width >> 1),
            (central_picture_ptr->height >> 1),
            central_picture_ptr->stride_cb,
            encoder_bit_depth);
        }
        }
    }
    else {
        EbByte buffer_y = central_picture_ptr->buffer_y +
            central_picture_ptr->origin_y * central_picture_ptr->stride_y +
            central_picture_ptr->origin_x;
        EbByte buffer_u =
            central_picture_ptr->buffer_cb +
            (central_picture_ptr->origin_y >> ss_y) * central_picture_ptr->stride_cb +
            (central_picture_ptr->origin_x >> ss_x);
        EbByte buffer_v =
            central_picture_ptr->buffer_cr +
            (central_picture_ptr->origin_y >> ss_x) * central_picture_ptr->stride_cr +
            (central_picture_ptr->origin_x >> ss_x);

        if (pcs_ptr->tf_ctrls.use_fixed_point || pcs_ptr->tf_ctrls.use_medium_filter) {
        if (do_noise_est)
        {
            noise_level_fp16 = svt_estimate_noise_fp16(buffer_y, // Y
                                                       central_picture_ptr->width,
                                                       central_picture_ptr->height,
                                                       central_picture_ptr->stride_y);
            noise_levels_log1p_fp16[C_Y] = noise_log1p_fp16(noise_level_fp16);
        }
        } else
            if (do_noise_est)
            noise_levels[0] = svt_estimate_noise(buffer_y, // Y
            central_picture_ptr->width,
            central_picture_ptr->height,
            central_picture_ptr->stride_y);
        if (pcs_ptr->tf_ctrls.do_chroma) {
        if (pcs_ptr->tf_ctrls.use_fixed_point || pcs_ptr->tf_ctrls.use_medium_filter) {
            noise_level_fp16 = svt_estimate_noise_fp16(buffer_u, // U
                                                      (central_picture_ptr->width >> ss_x),
                                                      (central_picture_ptr->height >> ss_y),
                                                      central_picture_ptr->stride_cb);
            noise_levels_log1p_fp16[C_U] = noise_log1p_fp16(noise_level_fp16);

            noise_level_fp16 = svt_estimate_noise_fp16(buffer_v, // V
                                                      (central_picture_ptr->width >> ss_x),
                                                      (central_picture_ptr->height >> ss_y),
                                                      central_picture_ptr->stride_cr);
            noise_levels_log1p_fp16[C_V] = noise_log1p_fp16(noise_level_fp16);
        } else {
        noise_levels[1] = svt_estimate_noise(buffer_u, // U
            (central_picture_ptr->width >> ss_x),
            (central_picture_ptr->height >> ss_y),
            central_picture_ptr->stride_cb);

        noise_levels[2] = svt_estimate_noise(buffer_v, // V
            (central_picture_ptr->width >> ss_x),
            (central_picture_ptr->height >> ss_y),
            central_picture_ptr->stride_cr);
        }
        }
    }
    if(pcs_ptr->tf_ctrls.use_fixed_point || pcs_ptr->tf_ctrls.use_medium_filter) {
        if (do_noise_est) {
            context_ptr->last_i_noise_levels_log1p_fp16[0] = noise_levels_log1p_fp16[0];
        } else {
            noise_levels_log1p_fp16[0] = context_ptr->last_i_noise_levels_log1p_fp16[0];
        }
    } else
        if (do_noise_est)
            context_ptr->last_i_noise_levels[0] = noise_levels[0];
       else
           noise_levels[0] = context_ptr->last_i_noise_levels[0];

    // Set is_noise_level for the tf off case
    pcs_ptr->is_noise_level = (context_ptr->last_i_noise_levels_log1p_fp16[0] >= VQ_NOISE_LVL_TH);
    // Adjust number of filtering frames based on noise and quantization factor.
    // Basically, we would like to use more frames to filter low-noise frame such
    // that the filtered frame can provide better predictions for more frames.
    // Also, when the quantization factor is small enough (lossless compression),
    // we will not change the number of frames for key frame filtering, which is
    // to avoid visual quality drop.
    int adjust_num = 0;
    if (pcs_ptr->tf_ctrls.use_fixed_point || pcs_ptr->tf_ctrls.use_medium_filter) {
        if (noise_levels_log1p_fp16[0] < 26572 /*FLOAT2FP(log1p(0.5), 16, int32_t)*/) {
            adjust_num = 6;
        }
        else if (noise_levels_log1p_fp16[0] < 45426 /*FLOAT2FP(log1p(1.0), 16, int32_t)*/) {
            adjust_num = 4;
        }
        else if (noise_levels_log1p_fp16[0] < 71998 /*FLOAT2FP(log1p(2.0), 16, int32_t)*/) {
            adjust_num = 2;
        }
    } else
    if (noise_levels[0] < 0.5) {
        adjust_num = 6;
    }
    else if (noise_levels[0] < 1.0) {
        adjust_num = 4;
    }
    else if (noise_levels[0] < 2.0) {
        adjust_num = 2;
    }
    (void)out_stride_diff64;
    if (scs_ptr->static_config.pred_structure != SVT_AV1_PRED_RANDOM_ACCESS) {
        int num_past_pics = pcs_ptr->tf_ctrls.num_past_pics + (pcs_ptr->tf_ctrls.noise_adjust_past_pics ? (adjust_num >> 1) : 0);
        num_past_pics = MIN(pcs_ptr->tf_ctrls.max_num_past_pics, num_past_pics);

        int num_future_pics = pcs_ptr->tf_ctrls.num_future_pics + (pcs_ptr->tf_ctrls.noise_adjust_future_pics ? (adjust_num >> 1) : 0);
        num_future_pics = MIN(pcs_ptr->tf_ctrls.max_num_future_pics, num_future_pics);

        //initilize list
        for (int pic_itr = 0; pic_itr < ALTREF_MAX_NFRAMES; pic_itr++)
            pcs_ptr->temp_filt_pcs_list[pic_itr] = NULL;

        //get previous
        for (int pic_itr = 0; pic_itr < num_past_pics; pic_itr++) {
            int32_t idx = search_this_pic(context_ptr->tf_pic_array, context_ptr->tf_pic_arr_cnt, pcs_ptr->picture_number - num_past_pics + pic_itr);
            if (idx >= 0)
                pcs_ptr->temp_filt_pcs_list[pic_itr] = context_ptr->tf_pic_array[idx];
        }

        //get central
        pcs_ptr->temp_filt_pcs_list[num_past_pics] = pcs_ptr;

        int actual_past_pics = num_past_pics;
        int actual_future_pics = 0;
        int pic_i;
        //search reord-queue to get the future pictures
        for (pic_i = 0; pic_i < num_future_pics; pic_i++) {
            int32_t q_index = QUEUE_GET_NEXT_SPOT(pcs_ptr->pic_decision_reorder_queue_idx, pic_i + 1);
            if (encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr != NULL) {
                PictureParentControlSet* pcs_itr = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr->object_ptr;
                pcs_ptr->temp_filt_pcs_list[pic_i + num_past_pics + 1] = pcs_itr;
                actual_future_pics++;
            }
            else
                break;
        }

        //search in pre-ass if still short
        if (pic_i < num_future_pics) {
            actual_future_pics = 0;
            for (int pic_i_future = 0; pic_i_future < num_future_pics; pic_i_future++) {
                for (uint32_t pic_i_pa = 0; pic_i_pa < encode_context_ptr->pre_assignment_buffer_count; pic_i_pa++) {
                    PictureParentControlSet* pcs_itr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pic_i_pa]->object_ptr;
                    if (pcs_itr->picture_number == pcs_ptr->picture_number + pic_i_future + 1) {
                        pcs_ptr->temp_filt_pcs_list[pic_i_future + num_past_pics + 1] = pcs_itr;
                        actual_future_pics++;
                        break; //exist the pre-ass loop, go search the next
                    }
                }
            }
        }
        pcs_ptr->past_altref_nframes = actual_past_pics;
        pcs_ptr->future_altref_nframes = actual_future_pics;

        // adjust the temporal filtering pcs buffer to remove unused past pictures
        if (actual_past_pics != num_past_pics) {

            pic_i = 0;
            while (pcs_ptr->temp_filt_pcs_list[pic_i] != NULL) {
                pcs_ptr->temp_filt_pcs_list[pic_i] = pcs_ptr->temp_filt_pcs_list[pic_i + num_past_pics - actual_past_pics];
                pic_i++;
            }
        }
    }
    else {
        if (is_delayed_intra(pcs_ptr)) {
            //initilize list
            for (int pic_itr = 0; pic_itr < ALTREF_MAX_NFRAMES; pic_itr++)
                pcs_ptr->temp_filt_pcs_list[pic_itr] = NULL;

            pcs_ptr->temp_filt_pcs_list[0] = pcs_ptr;
            uint32_t num_future_pics = pcs_ptr->tf_ctrls.num_future_pics + (pcs_ptr->tf_ctrls.noise_adjust_future_pics ? adjust_num : 0);
            num_future_pics = MIN(pcs_ptr->tf_ctrls.max_num_future_pics, num_future_pics);

            uint32_t pic_i;
            for (pic_i = 0; pic_i < num_future_pics; pic_i++) {
                int32_t idx = search_this_pic(context_ptr->mg_pictures_array, context_ptr->mg_size, pcs_ptr->picture_number + pic_i + 1);
                if (idx >= 0)
                    pcs_ptr->temp_filt_pcs_list[pic_i + 1] = context_ptr->mg_pictures_array[idx];
                else
                    break;
            }

            pcs_ptr->past_altref_nframes = 0;
            pcs_ptr->future_altref_nframes = pic_i;
        }
        else
            if (pcs_ptr->idr_flag) {

                //initilize list
                for (int pic_itr = 0; pic_itr < ALTREF_MAX_NFRAMES; pic_itr++)
                    pcs_ptr->temp_filt_pcs_list[pic_itr] = NULL;

                pcs_ptr->temp_filt_pcs_list[0] = pcs_ptr;
                uint32_t num_future_pics = pcs_ptr->tf_ctrls.num_future_pics + (pcs_ptr->tf_ctrls.noise_adjust_future_pics ? adjust_num : 0);
                num_future_pics = MIN(pcs_ptr->tf_ctrls.max_num_future_pics, num_future_pics);
                uint32_t num_past_pics = 0;
                uint32_t pic_i;
                //search reord-queue to get the future pictures
                for (pic_i = 0; pic_i < num_future_pics; pic_i++) {
                    int32_t q_index = QUEUE_GET_NEXT_SPOT(pcs_ptr->pic_decision_reorder_queue_idx, pic_i + 1);
                    if (encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr != NULL) {
                        PictureParentControlSet* pcs_itr = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr->object_ptr;
                        pcs_ptr->temp_filt_pcs_list[pic_i + num_past_pics + 1] = pcs_itr;

                    }
                    else
                        break;
                }

                pcs_ptr->past_altref_nframes = 0;
                pcs_ptr->future_altref_nframes = pic_i;
            }
            else
            {
                int num_past_pics = pcs_ptr->tf_ctrls.num_past_pics + (pcs_ptr->tf_ctrls.noise_adjust_past_pics ? (adjust_num >> 1) : 0);
                num_past_pics = MIN(pcs_ptr->tf_ctrls.max_num_past_pics, num_past_pics);

                int num_future_pics = pcs_ptr->tf_ctrls.num_future_pics + (pcs_ptr->tf_ctrls.noise_adjust_future_pics ? (adjust_num >> 1) : 0);
                num_future_pics = MIN(pcs_ptr->tf_ctrls.max_num_future_pics, num_future_pics);

                //initilize list
                for (int pic_itr = 0; pic_itr < ALTREF_MAX_NFRAMES; pic_itr++)
                    pcs_ptr->temp_filt_pcs_list[pic_itr] = NULL;
                // limit the number of pictures to make sure there are enough pictures in the buffer. i.e. Intra CRA case
                // limit the number of pictures to make sure there are enough pictures in the buffer. i.e. Intra CRA case
                num_past_pics = MIN(num_past_pics, avail_past_pictures(context_ptr->mg_pictures_array, context_ptr->mg_size, pcs_ptr->picture_number));
                //get previous+current pictures from the the pre-assign buffer
                for (int pic_itr = 0; pic_itr <= num_past_pics; pic_itr++) {
                    int32_t idx = search_this_pic(context_ptr->mg_pictures_array, context_ptr->mg_size, pcs_ptr->picture_number - num_past_pics + pic_itr);
                    if (idx >= 0)
                        pcs_ptr->temp_filt_pcs_list[pic_itr] = context_ptr->mg_pictures_array[idx];
                }
                int actual_past_pics = num_past_pics;
                int actual_future_pics = 0;
                int pic_i;
                //search reord-queue to get the future pictures
                for (pic_i = 0; pic_i < num_future_pics; pic_i++) {
                    int32_t q_index = QUEUE_GET_NEXT_SPOT(pcs_ptr->pic_decision_reorder_queue_idx, pic_i + 1);
                    if (encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr != NULL) {
                        PictureParentControlSet* pcs_itr = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr->object_ptr;
                        pcs_ptr->temp_filt_pcs_list[pic_i + num_past_pics + 1] = pcs_itr;
                        actual_future_pics++;
                    }
                    else
                        break;
                }

                //search in pre-ass if still short
                if (pic_i < num_future_pics) {
                    actual_future_pics = 0;
                    for (int pic_i_future = 0; pic_i_future < num_future_pics; pic_i_future++) {
                        for (uint32_t pic_i_pa = 0; pic_i_pa < encode_context_ptr->pre_assignment_buffer_count; pic_i_pa++) {
                            PictureParentControlSet* pcs_itr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pic_i_pa]->object_ptr;
                            if (pcs_itr->picture_number == pcs_ptr->picture_number + pic_i_future + 1) {
                                pcs_ptr->temp_filt_pcs_list[pic_i_future + num_past_pics + 1] = pcs_itr;
                                actual_future_pics++;
                                break; //exist the pre-ass loop, go search the next
                            }
                        }
                    }
                }
                pcs_ptr->past_altref_nframes = actual_past_pics;
                pcs_ptr->future_altref_nframes = actual_future_pics;

                // adjust the temporal filtering pcs buffer to remove unused past pictures
                if (actual_past_pics != num_past_pics) {

                    pic_i = 0;
                    while (pcs_ptr->temp_filt_pcs_list[pic_i] != NULL) {
                        pcs_ptr->temp_filt_pcs_list[pic_i] = pcs_ptr->temp_filt_pcs_list[pic_i + num_past_pics - actual_past_pics];
                        pic_i++;
                    }
                }
            }

    }
    return EB_ErrorNone;
}

PaReferenceQueueEntry * search_ref_in_ref_queue_pa(
    EncodeContext *encode_context_ptr,
    uint64_t ref_poc)
{
    PaReferenceQueueEntry * ref_entry_ptr = NULL;
    uint32_t ref_queue_i = encode_context_ptr->picture_decision_pa_reference_queue_head_index;
    // Find the Reference in the Reference Queue
    do {
        ref_entry_ptr =
            encode_context_ptr->picture_decision_pa_reference_queue[ref_queue_i];
        if (ref_entry_ptr->picture_number == ref_poc)
            return ref_entry_ptr;

        // Increment the reference_queue_index Iterator
        ref_queue_i = (ref_queue_i == REFERENCE_QUEUE_MAX_DEPTH - 1)
            ? 0
            : ref_queue_i + 1;

    } while (ref_queue_i != encode_context_ptr->picture_decision_pa_reference_queue_tail_index);


    return NULL;
}
/*
 * Copy TF params: sps -> pcs
 */
void copy_tf_params(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {

    // Map TF settings sps -> pcs
   if (scs_ptr->static_config.pred_structure != SVT_AV1_PRED_RANDOM_ACCESS)
   {
        if (pcs_ptr->slice_type != I_SLICE && pcs_ptr->temporal_layer_index == 0)
            pcs_ptr->tf_ctrls = scs_ptr->tf_params_per_type[1];
        else
            pcs_ptr->tf_ctrls.enabled = 0;

        return;
   }
    if (is_delayed_intra(pcs_ptr))
        pcs_ptr->tf_ctrls = scs_ptr->tf_params_per_type[0];
    else if (pcs_ptr->temporal_layer_index == 0)  // BASE
        pcs_ptr->tf_ctrls = scs_ptr->tf_params_per_type[1];
    else if (pcs_ptr->temporal_layer_index == 1)  // L1
        pcs_ptr->tf_ctrls = scs_ptr->tf_params_per_type[2];
    else
        pcs_ptr->tf_ctrls.enabled = 0;
}
void is_screen_content(PictureParentControlSet *pcs_ptr);
/*
* Update the list0 count try and the list1 count try based on the Enc-Mode, whether BASE or not, whether SC or not
*/
void update_count_try(SequenceControlSet* scs_ptr, PictureParentControlSet* pcs_ptr) {
    MrpCtrls* mrp_ctrl = &scs_ptr->mrp_ctrls;
    if (pcs_ptr->sc_class1) {
        if (pcs_ptr->temporal_layer_index == 0) {
            pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, mrp_ctrl->sc_base_ref_list0_count);
            pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, mrp_ctrl->sc_base_ref_list1_count);
        }
        else {
            pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, mrp_ctrl->sc_non_base_ref_list0_count);
            pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, mrp_ctrl->sc_non_base_ref_list1_count);
        }
    }
    else {
        if (pcs_ptr->temporal_layer_index == 0) {
            pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, mrp_ctrl->base_ref_list0_count);
            pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, mrp_ctrl->base_ref_list1_count);
        }
        else {
            pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, mrp_ctrl->non_base_ref_list0_count);
            pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, mrp_ctrl->non_base_ref_list1_count);
        }
    }
}
/*
* Switch frame's pcs_ptr->dpb_order_hint[8] will be packed to uncompressed_header as ref_order_hint[8], ref to spec 5.9.2.
* Pictures are inputted in this process in display order and no need to consider reordering since the switch frame feature only supports low delay pred structure by design (not by spec).
*/
static void update_sframe_ref_order_hint(PictureParentControlSet *ppcs, PictureDecisionContext *context_ptr)
{
    assert(sizeof(ppcs->dpb_order_hint) == sizeof(context_ptr->ref_order_hint));
    memcpy(ppcs->dpb_order_hint, context_ptr->ref_order_hint, sizeof(ppcs->dpb_order_hint));
    if (ppcs->av1_ref_signal.refresh_frame_mask != 0) {
        const uint32_t cur_order_hint = ppcs->picture_number % ((uint64_t)1 << (ppcs->scs_ptr->seq_header.order_hint_info.order_hint_bits));
        for (int32_t i = 0; i < REF_FRAMES; i++) {
            if ((ppcs->av1_ref_signal.refresh_frame_mask >> i) & 1) {
                context_ptr->ref_order_hint[i] = cur_order_hint;
            }
        }
    }
}

/*****************************************************************
* Update the RC param queue
* Set the size of the previous Gop/param, Check if all the frames in gop are processed, if yes reset
* Increament the head index to assign a new spot in the queue for the new gop
*****************************************************************/
static void update_rc_param_queue(
    PictureParentControlSet *ppcs,
    EncodeContext           *enc_cxt) {
    if (ppcs->idr_flag == TRUE && ppcs->picture_number > 0) {
        svt_block_on_mutex(enc_cxt->rc_param_queue_mutex);
        // Set the size of the previous Gop/param
        enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index]->size =
            (int32_t)(ppcs->picture_number - enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index]->first_poc);
        // Check if all the frames in gop are processed, if yes reset
        if (enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index]->size ==
            enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index]->processed_frame_number) {
            enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index]->size = -1;
            enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index]->processed_frame_number = 0;
        }
        // Increament the head index to assign a new spot in the queue for the new gop
        enc_cxt->rc_param_queue_head_index = (enc_cxt->rc_param_queue_head_index == PARALLEL_GOP_MAX_NUMBER - 1) ?
            0 : enc_cxt->rc_param_queue_head_index + 1;
        assert_err(enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index]->size == -1, "The head in rc paramqueue is not empty");
        enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index]->first_poc = ppcs->picture_number;
        svt_release_mutex(enc_cxt->rc_param_queue_mutex);

    }
    // Store the pointer to the right spot in the RC param queue under PCS
    ppcs->rate_control_param_ptr = enc_cxt->rc_param_queue[enc_cxt->rc_param_queue_head_index];
}

/****************************************************************************************
* set_layer_depth()
* Set the layer depth per frame based on frame type, temporal layer
****************************************************************************************/
static void set_layer_depth(PictureParentControlSet *ppcs) {
   // SequenceControlSet *scs = ppcs->scs_ptr;
    if (ppcs->frm_hdr.frame_type == KEY_FRAME)
        ppcs->layer_depth = 0;
    else
        ppcs->layer_depth = ppcs->temporal_layer_index + 1;
}
/****************************************************************************************
* set_frame_update_type()
* Set the update type per frame based on frame type, temporal layer and prediction structure
* For Low delay, there is a special case where all non key frames are treated as LF_UPDATE.
* Every MAX_GF_INTERVAL frames, update type is set to GF_UPDATE
****************************************************************************************/
static void set_frame_update_type(PictureParentControlSet *ppcs) {
    SequenceControlSet *scs = ppcs->scs_ptr;
    if (ppcs->frm_hdr.frame_type == KEY_FRAME) {
        ppcs->update_type = KF_UPDATE;
    }
    else if (scs->max_temporal_layers > 0 && ppcs->pred_structure != SVT_AV1_PRED_LOW_DELAY_B) {
        if (ppcs->temporal_layer_index == 0) {
            ppcs->update_type = ARF_UPDATE;
        }
        else if (ppcs->temporal_layer_index == ppcs->hierarchical_levels) {
            ppcs->update_type = LF_UPDATE;
        }
        else {
            ppcs->update_type = INTNL_ARF_UPDATE;
        }
    }
    else if (ppcs->pred_structure == SVT_AV1_PRED_LOW_DELAY_B && (ppcs->frame_offset % MAX_GF_INTERVAL) == 0) {
        ppcs->update_type = GF_UPDATE;
    }
    else {
        ppcs->update_type = LF_UPDATE;
    }
}
static void set_gf_group_param(PictureParentControlSet *ppcs) {
    set_frame_update_type(ppcs);
    set_layer_depth(ppcs);
}

/* Picture Decision Kernel */

/***************************************************************************************************
*
* @brief
*  The Picture Decision process performs multi-picture level decisions, including setting of the prediction structure,
*  setting the picture type and performing scene change detection.
*
* @par Description:
*  Since the prior Picture Analysis process stage is multithreaded, inputs to the Picture Decision Process can arrive
*  out-of-display-order, so a reordering queue is used to enforce processing of pictures in display order. The algorithms
*  employed in the Picture Decision process are dependent on prior pictures’ statistics, so the order in which pictures are
*  processed must be strictly enforced. Additionally, the Picture Decision process uses the reorder queue to hold input pictures
*  until they can be started into the Motion Estimation process while following the proper prediction structure.
*
* @param[in] Pictures
*  The Picture Decision Process takes images spontaneously as they arive and perform multi-picture level decisions,
*  including setting of the picture structure, setting the picture type and scene change detection.
*
* @param[out] Picture Control Set
*  Picture Control Set with fully available Picture Analysis Reference List
*
* @remarks
*  For Low Delay Sequences, pictures are started into the encoder pipeline immediately.
*
*  For Random Access Sequences, pictures are held for up to a PredictionStructurePeriod
*    in order to determine if a Scene Change or Intra Frame is forthcoming. Either of
*    those events (and additionally a End of Sequence Flag) will change the expected
*    prediction structure.
*
*  Below is an example worksheet for how Intra Flags and Scene Change Flags interact
*    together to affect the prediction structure.
*
*  The base prediction structure for this example is a 3-Level Hierarchical Random Access,
*    Single Reference Prediction Structure:
*
*        b   b
*       / \ / \
*      /   b   \
*     /   / \   \
*    I-----------b
*
*  From this base structure, the following RPS positions are derived:
*
*    p   p       b   b       p   p
*     \   \     / \ / \     /   /
*      P   \   /   b   \   /   P
*       \   \ /   / \   \ /   /
*        ----I-----------b----
*
*    L L L   I  [ Normal ]   T T T
*    2 1 0   n               0 1 2
*            t
*            r
*            a
*
*  The RPS is composed of Leading Picture [L2-L0], Intra (CRA), Base/Normal Pictures,
*    and Trailing Pictures [T0-t2]. Generally speaking, Leading Pictures are useful
*    for handling scene changes without adding extraneous I-pictures and the Trailing
*    pictures are useful for terminating GOPs.
*
*  Here is a table of possible combinations of pictures needed to handle intra and
*    scene changes happening in quick succession.
*
*        Distance to scene change ------------>
*
*                  0              1                 2                3+
*   I
*   n
*   t   0        I   I           n/a               n/a              n/a
*   r
*   a              p              p
*                   \            /
*   P   1        I   I          I   I              n/a              n/a
*   e
*   r               p                               p
*   i                \                             /
*   o            p    \         p   p             /   p
*   d             \    \       /     \           /   /
*       2     I    -----I     I       I         I----    I          n/a
*   |
*   |            p   p           p   p            p   p            p   p
*   |             \   \         /     \          /     \          /   /
*   |              P   \       /   p   \        /   p   \        /   P
*   |               \   \     /     \   \      /   /     \      /   /
*   V   3+   I       ----I   I       ----I    I----       I    I----       I
*
*   The table is interpreted as follows:
*
*   If there are no SCs or Intras encountered for a PredPeriod, then the normal
*     prediction structure is applied.
*
*   If there is an intra in the PredPeriod, then one of the above combinations of
*     Leading and Trailing pictures is used.  If there is no scene change, the last
*     valid column consisting of Trailing Pictures only is used.  However, if there
*     is an upcoming scene change before the next intra, then one of the above patterns
*     is used. In the case of End of Sequence flags, only the last valid column of Trailing
*     Pictures is used. The intention here is that any combination of Intra Flag and Scene
*     Change flag can be coded.
***************************************************************************************************/
void* picture_decision_kernel(void *input_ptr)
{
    EbThreadContext               *thread_context_ptr = (EbThreadContext*)input_ptr;
    PictureDecisionContext        *context_ptr = (PictureDecisionContext*)thread_context_ptr->priv;

    PictureParentControlSet       *pcs_ptr;

    EncodeContext                 *encode_context_ptr;
    SequenceControlSet            *scs_ptr;
    EbObjectWrapper               *in_results_wrapper_ptr;
    PictureAnalysisResults        *in_results_ptr;

    PredictionStructureEntry      *pred_position_ptr;

    Bool                          pre_assignment_buffer_first_pass_flag;
    SliceType                         picture_type;

    PictureDecisionReorderEntry   *queue_entry_ptr;
    int32_t                           queue_entry_index;

    int32_t                           previous_entry_index;

    PaReferenceQueueEntry         *input_entry_ptr = (PaReferenceQueueEntry*)NULL;;
    uint32_t                           input_queue_index;

    PaReferenceQueueEntry         *pa_reference_entry_ptr;
    uint64_t                           ref_poc;

    uint32_t                           dep_idx;
    int64_t                            dep_poc;

    uint32_t                           dep_list_count;

    // Dynamic GOP
    uint32_t                           mini_gop_index;
    uint32_t                           out_stride_diff64;
    int64_t                            current_input_poc = -1;

    Bool                          window_avail, frame_passthrough;
    uint32_t                           window_index;
    uint32_t                           entry_index;
    // Debug
    uint64_t                           loop_count = 0;

    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(
            context_ptr->picture_analysis_results_input_fifo_ptr,
            &in_results_wrapper_ptr);

        in_results_ptr = (PictureAnalysisResults*)in_results_wrapper_ptr->object_ptr;
        pcs_ptr = (PictureParentControlSet*)in_results_ptr->pcs_wrapper_ptr->object_ptr;
        scs_ptr = pcs_ptr->scs_ptr;
        encode_context_ptr = (EncodeContext*)scs_ptr->encode_context_ptr;
        loop_count++;

        // Input Picture Analysis Results into the Picture Decision Reordering Queue
        // P.S. Since the prior Picture Analysis processes stage is multithreaded, inputs to the Picture Decision Process
        // can arrive out-of-display-order, so a the Picture Decision Reordering Queue is used to enforce processing of
        // pictures in display order
        if (!pcs_ptr->is_overlay ) {
            queue_entry_index = (int32_t)(pcs_ptr->picture_number - encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index]->picture_number);
            queue_entry_index += encode_context_ptr->picture_decision_reorder_queue_head_index;
            queue_entry_index = (queue_entry_index > PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1) ? queue_entry_index - PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH : queue_entry_index;
            queue_entry_ptr = encode_context_ptr->picture_decision_reorder_queue[queue_entry_index];
            if (queue_entry_ptr->parent_pcs_wrapper_ptr != NULL) {
                CHECK_REPORT_ERROR_NC(
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PD_ERROR8);
            }
            else {
                queue_entry_ptr->parent_pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;
                queue_entry_ptr->picture_number = pcs_ptr->picture_number;
            }

            pcs_ptr->pic_decision_reorder_queue_idx = queue_entry_index;
            pcs_ptr->first_pass_done = 0;
        }
        // Process the head of the Picture Decision Reordering Queue (Entry N)
        // P.S. The Picture Decision Reordering Queue should be parsed in the display order to be able to construct a pred structure
        queue_entry_ptr = encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index];

        while (queue_entry_ptr->parent_pcs_wrapper_ptr != NULL) {
            if (((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr->object_ptr))->end_of_sequence_flag == TRUE) {
                frame_passthrough = TRUE;
            }
            else
                frame_passthrough = FALSE;
            window_avail = TRUE;
            previous_entry_index = QUEUE_GET_PREVIOUS_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index);
            if (scs_ptr->static_config.pass == ENC_FIRST_PASS || scs_ptr->lap_rc) {
                for (window_index = 0; window_index < scs_ptr->scd_delay + 1; window_index++)
                {
                    entry_index = QUEUE_GET_NEXT_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index, window_index);
                    PictureDecisionReorderEntry   *first_pass_queue_entry = encode_context_ptr->picture_decision_reorder_queue[entry_index];
                    if (first_pass_queue_entry->parent_pcs_wrapper_ptr == NULL) {
                        break;
                    }
                    else {

                        PictureParentControlSet *first_pass_pcs_ptr = (PictureParentControlSet*)first_pass_queue_entry->parent_pcs_wrapper_ptr->object_ptr;
                        if (!first_pass_pcs_ptr->first_pass_done) {
                            int32_t temp_entry_index = QUEUE_GET_PREVIOUS_SPOT(entry_index);
                            first_pass_pcs_ptr->first_pass_ref_ppcs_ptr[0] = first_pass_queue_entry->picture_number > 0 ?
                                (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[temp_entry_index]->parent_pcs_wrapper_ptr->object_ptr : NULL;
                            temp_entry_index = QUEUE_GET_PREVIOUS_SPOT(temp_entry_index);
                            first_pass_pcs_ptr->first_pass_ref_ppcs_ptr[1] = first_pass_queue_entry->picture_number > 1 ?
                                (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[temp_entry_index]->parent_pcs_wrapper_ptr->object_ptr : NULL;
                            first_pass_pcs_ptr->first_pass_ref_count = first_pass_queue_entry->picture_number > 1 ? 2 : first_pass_queue_entry->picture_number > 0 ? 1 : 0;

                            process_first_pass_frame(scs_ptr, first_pass_pcs_ptr, context_ptr);
                            first_pass_pcs_ptr->first_pass_done = 1;
                        }
                    }
                }
            }

            pcs_ptr = (PictureParentControlSet*)queue_entry_ptr->parent_pcs_wrapper_ptr->object_ptr;
            memset(pcs_ptr->pd_window, 0, PD_WINDOW_SIZE * sizeof(PictureParentControlSet*));

            //for poc 0, ignore previous frame check
            if (queue_entry_ptr->picture_number > 0 && encode_context_ptr->picture_decision_reorder_queue[previous_entry_index]->parent_pcs_wrapper_ptr == NULL)
                window_avail = FALSE;
            else {

                //TODO: risk of a race condition accessing prev(pcs0 is released, and pcs1 still doing sc).
                //Actually we dont need to keep prev, just keep previous copy of histograms.
                pcs_ptr->pd_window[0] =
                    queue_entry_ptr->picture_number > 0 ? (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[previous_entry_index]->parent_pcs_wrapper_ptr->object_ptr : NULL;
                pcs_ptr->pd_window[1] =
                    (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index]->parent_pcs_wrapper_ptr->object_ptr;
                for (window_index = 0; window_index < scs_ptr->scd_delay; window_index++) {
                    entry_index = QUEUE_GET_NEXT_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index, window_index + 1);
                    if (encode_context_ptr->picture_decision_reorder_queue[entry_index]->parent_pcs_wrapper_ptr == NULL) {
                        window_avail = FALSE;
                        break;
                    }
                    else if (((PictureParentControlSet *)(encode_context_ptr->picture_decision_reorder_queue[entry_index]->parent_pcs_wrapper_ptr->object_ptr))->end_of_sequence_flag == TRUE) {
                        window_avail = FALSE;
                        frame_passthrough = TRUE;
                        break;
                    }else {
                        pcs_ptr->pd_window[2 + window_index] =
                            (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[entry_index]->parent_pcs_wrapper_ptr->object_ptr;
                    }
                }
            }

            pcs_ptr = (PictureParentControlSet*)queue_entry_ptr->parent_pcs_wrapper_ptr->object_ptr;
            if (pcs_ptr->idr_flag == TRUE)
                context_ptr->last_solid_color_frame_poc = 0xFFFFFFFF;
            if (window_avail == TRUE && queue_entry_ptr->picture_number > 0) {
                if (scs_ptr->static_config.scene_change_detection) {

                    pcs_ptr->scene_change_flag = scene_transition_detector(
                        context_ptr,
                        scs_ptr,
                        (PictureParentControlSet**)pcs_ptr->pd_window);

                }
                else if (scs_ptr->vq_ctrls.sharpness_ctrls.scene_transition && (context_ptr->transition_detected == -1 || context_ptr->transition_detected == 0)) {
                    context_ptr->transition_detected = scene_transition_detector(
                        context_ptr,
                        scs_ptr,
                        (PictureParentControlSet**)pcs_ptr->pd_window);
                }
                else
                    pcs_ptr->scene_change_flag = FALSE;

                pcs_ptr->cra_flag = (pcs_ptr->scene_change_flag == TRUE) ?
                    TRUE :
                    pcs_ptr->cra_flag;

                // Store scene change in context
                context_ptr->is_scene_change_detected = pcs_ptr->scene_change_flag;
            }

            if (window_avail == TRUE || frame_passthrough == TRUE)
            {
                // Place the PCS into the Pre-Assignment Buffer
                // P.S. The Pre-Assignment Buffer is used to store a whole pre-structure
                encode_context_ptr->pre_assignment_buffer[encode_context_ptr->pre_assignment_buffer_count] = queue_entry_ptr->parent_pcs_wrapper_ptr;

                // Setup the PCS & SCS
                pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[encode_context_ptr->pre_assignment_buffer_count]->object_ptr;
                // Set the POC Number
                pcs_ptr->picture_number = ++current_input_poc;

                pcs_ptr->pred_structure = scs_ptr->static_config.pred_structure;

                pcs_ptr->hierarchical_layers_diff = 0;

                pcs_ptr->init_pred_struct_position_flag = FALSE;

                pcs_ptr->self_updated_links = 0;
                pcs_ptr->other_updated_links_cnt = 0;

                pcs_ptr->tpl_group_size = 0;
                if (pcs_ptr->picture_number == 0)
                    context_ptr->prev_delayed_intra = NULL;

                release_prev_picture_from_reorder_queue(
                    encode_context_ptr);

                // If the Intra period length is 0, then introduce an intra for every picture
                if (scs_ptr->static_config.intra_period_length == 0)
                    pcs_ptr->cra_flag = TRUE;
                // If an #IntraPeriodLength has passed since the last Intra, then introduce a CRA or IDR based on Intra Refresh type
                else if (scs_ptr->static_config.intra_period_length != -1) {

                    pcs_ptr->cra_flag =
                        (scs_ptr->static_config.intra_refresh_type != SVT_AV1_FWDKF_REFRESH) ?
                        pcs_ptr->cra_flag :
                        ((encode_context_ptr->intra_period_position == (uint32_t)scs_ptr->static_config.intra_period_length) || (pcs_ptr->scene_change_flag == TRUE)) ?
                        TRUE :
                        pcs_ptr->cra_flag;

                    pcs_ptr->idr_flag =
                        (scs_ptr->static_config.intra_refresh_type != SVT_AV1_KF_REFRESH) ?
                        pcs_ptr->idr_flag :
                        ((encode_context_ptr->intra_period_position == (uint32_t)scs_ptr->static_config.intra_period_length) ||
                         (pcs_ptr->scene_change_flag == TRUE) ||
                         (scs_ptr->static_config.force_key_frames && pcs_ptr->input_ptr->pic_type == EB_AV1_KEY_PICTURE)) ?
                        TRUE :
                        pcs_ptr->idr_flag;
                }

                encode_context_ptr->pre_assignment_buffer_eos_flag = (pcs_ptr->end_of_sequence_flag) ? (uint32_t)TRUE : encode_context_ptr->pre_assignment_buffer_eos_flag;

                // Histogram data to be used at the next input (N + 1)
                if (scs_ptr->static_config.scene_change_detection || scs_ptr->vq_ctrls.sharpness_ctrls.scene_transition) {
                    for (int region_in_picture_width_index = 0; region_in_picture_width_index < MAX_NUMBER_OF_REGIONS_IN_WIDTH; region_in_picture_width_index++) {
                        for (int region_in_picture_height_index = 0; region_in_picture_height_index < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; region_in_picture_height_index++) {

                            svt_memcpy(
                                &(context_ptr->prev_picture_histogram[region_in_picture_width_index][region_in_picture_height_index][0]),
                                &(pcs_ptr->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][0]),
                                HISTOGRAM_NUMBER_OF_BINS * sizeof(uint32_t));

                            context_ptr->prev_average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index] = pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index];

                        }
                    }
                }

                // Increment the Pre-Assignment Buffer Intra Count
                encode_context_ptr->pre_assignment_buffer_intra_count += (pcs_ptr->idr_flag || pcs_ptr->cra_flag);
                encode_context_ptr->pre_assignment_buffer_idr_count += pcs_ptr->idr_flag;
                encode_context_ptr->pre_assignment_buffer_count += 1;

                if (scs_ptr->static_config.rate_control_mode)
                {
                    // Increment the Intra Period Position
                    encode_context_ptr->intra_period_position = ((encode_context_ptr->intra_period_position == (uint32_t)scs_ptr->static_config.intra_period_length) || (pcs_ptr->scene_change_flag == TRUE)) ? 0 : encode_context_ptr->intra_period_position + 1;
                }
                else
                {
                    // Increment the Intra Period Position
                    encode_context_ptr->intra_period_position =
                        ((encode_context_ptr->intra_period_position == (uint32_t)scs_ptr->static_config.intra_period_length) ||
                         (pcs_ptr->scene_change_flag == TRUE) ||
                         (scs_ptr->static_config.force_key_frames && pcs_ptr->input_ptr->pic_type == EB_AV1_KEY_PICTURE)) ?
                        0 : encode_context_ptr->intra_period_position + 1;
                }

#if LAD_MG_PRINT
                 print_pre_ass_buffer(encode_context_ptr, pcs_ptr, 1);
#endif


                // Determine if Pictures can be released from the Pre-Assignment Buffer
                if ((encode_context_ptr->pre_assignment_buffer_intra_count > 0) ||
                    (encode_context_ptr->pre_assignment_buffer_count == (uint32_t)(1 << scs_ptr->static_config.hierarchical_levels)) ||
                    (encode_context_ptr->pre_assignment_buffer_eos_flag == TRUE) ||
                    (pcs_ptr->pred_structure == SVT_AV1_PRED_LOW_DELAY_P) ||
                    (pcs_ptr->pred_structure == SVT_AV1_PRED_LOW_DELAY_B))
                {
#if LAD_MG_PRINT
                    print_pre_ass_buffer(encode_context_ptr, pcs_ptr,0);
#endif
                    // Initialize Picture Block Params
                    context_ptr->mini_gop_start_index[0] = 0;
                    context_ptr->mini_gop_end_index[0] = encode_context_ptr->pre_assignment_buffer_count - 1;
                    context_ptr->mini_gop_length[0] = encode_context_ptr->pre_assignment_buffer_count;

                    context_ptr->mini_gop_hierarchical_levels[0] = scs_ptr->static_config.hierarchical_levels;
                    context_ptr->mini_gop_intra_count[0] = encode_context_ptr->pre_assignment_buffer_intra_count;
                    context_ptr->mini_gop_idr_count[0] = encode_context_ptr->pre_assignment_buffer_idr_count;
                    context_ptr->total_number_of_mini_gops = 1;

                    encode_context_ptr->previous_mini_gop_hierarchical_levels = (pcs_ptr->picture_number == 0) ?
                        scs_ptr->static_config.hierarchical_levels :
                        encode_context_ptr->previous_mini_gop_hierarchical_levels;

                    {
                        if (scs_ptr->static_config.hierarchical_levels == 1) {
                            //minigop 2 case
                            context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops] = 0;
                            context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - 1;
                            context_ptr->mini_gop_length[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops];
                            context_ptr->mini_gop_hierarchical_levels[context_ptr->total_number_of_mini_gops] = 2;
                        } else {
                            //minigop 4,8,16,32
                            if (encode_context_ptr->pre_assignment_buffer_count > 1) {
                                initialize_mini_gop_activity_array(
                                        context_ptr);

                                if (encode_context_ptr->pre_assignment_buffer_count >= 32 &&
                                    !(encode_context_ptr->pre_assignment_buffer_count == 32 && pcs_ptr->idr_flag))
                                    context_ptr->mini_gop_activity_array[L6_INDEX] = FALSE;
                                if (encode_context_ptr->pre_assignment_buffer_count >= 16 &&
                                    !(encode_context_ptr->pre_assignment_buffer_count == 16 && pcs_ptr->idr_flag))
                                    context_ptr->mini_gop_activity_array[L5_0_INDEX] = FALSE;
                                if (encode_context_ptr->pre_assignment_buffer_count >= 8 &&
                                    !(encode_context_ptr->pre_assignment_buffer_count == 8 && pcs_ptr->idr_flag)) {
                                    context_ptr->mini_gop_activity_array[L4_0_INDEX] = FALSE;
                                    context_ptr->mini_gop_activity_array[L4_1_INDEX] = FALSE;
                                }

                                generate_picture_window_split(
                                        context_ptr,
                                        encode_context_ptr);

                                handle_incomplete_picture_window_map(
                                        scs_ptr->static_config.hierarchical_levels,
                                        context_ptr,
                                        encode_context_ptr);
                            }
                        }
                    }

                    generate_mini_gop_rps(
                        context_ptr,
                        encode_context_ptr);


                   // SVT_LOG(" pre-assign split into  %i MGs \n", context_ptr->total_number_of_mini_gops);

                    // Loop over Mini GOPs

                    for (mini_gop_index = 0; mini_gop_index < context_ptr->total_number_of_mini_gops; ++mini_gop_index) {
                        pre_assignment_buffer_first_pass_flag = TRUE;
                        encode_context_ptr->is_mini_gop_changed = FALSE;
                        {
                            update_base_layer_reference_queue_dependent_count(
                                context_ptr,
                                encode_context_ptr,
                                scs_ptr,
                                mini_gop_index);

                            // Keep track of the number of hierarchical levels of the latest implemented mini GOP
                            encode_context_ptr->previous_mini_gop_hierarchical_levels = context_ptr->mini_gop_hierarchical_levels[mini_gop_index];
                        }

                        // 1st Loop over Pictures in the Pre-Assignment Buffer
                        for (out_stride_diff64 = context_ptr->mini_gop_start_index[mini_gop_index]; out_stride_diff64 <= context_ptr->mini_gop_end_index[mini_gop_index]; ++out_stride_diff64) {
                            Bool is_trailing_frame = FALSE;
                            pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[out_stride_diff64]->object_ptr;
                            scs_ptr = pcs_ptr->scs_ptr;
                            // Keep track of the mini GOP size to which the input picture belongs - needed @ PictureManagerProcess()
                            pcs_ptr->pre_assignment_buffer_count = context_ptr->mini_gop_length[mini_gop_index];

                            // Update the Pred Structure if cutting short a Random Access period
                            if (is_pic_cutting_short_ra_mg(context_ptr, pcs_ptr, mini_gop_index))
                            {
                                // Correct the Pred Index before switching structures
                                if (pre_assignment_buffer_first_pass_flag == TRUE)
                                    encode_context_ptr->pred_struct_position -= pcs_ptr->pred_struct_ptr->init_pic_index;
                                pcs_ptr->pred_struct_ptr = get_prediction_structure(
                                    encode_context_ptr->prediction_structure_group_ptr,
                                    SVT_AV1_PRED_LOW_DELAY_P,
                                    scs_ptr->reference_count,
                                    pcs_ptr->hierarchical_levels);
                                picture_type = P_SLICE;
                                if (scs_ptr->static_config.hierarchical_levels == 1 &&
                                    encode_context_ptr->prediction_structure_group_ptr->ref_count_used < MAX_REF_IDX) {
                                    // Only works for 1B case
                                    is_trailing_frame = TRUE;
                                    encode_context_ptr->pred_struct_position =
                                        pcs_ptr->pred_struct_ptr->steady_state_index + out_stride_diff64 - context_ptr->mini_gop_start_index[mini_gop_index];
                                }
                            }
                            // Open GOP CRA - adjust the RPS
                            else if ((context_ptr->mini_gop_length[mini_gop_index] == pcs_ptr->pred_struct_ptr->pred_struct_period) &&

                                (pcs_ptr->pred_struct_ptr->pred_type == SVT_AV1_PRED_RANDOM_ACCESS || pcs_ptr->pred_struct_ptr->temporal_layer_count == 1) &&
                                pcs_ptr->idr_flag == FALSE &&
                                pcs_ptr->cra_flag == TRUE)
                            {
                                picture_type = I_SLICE;
                            }
                            else {
                                // Set the Picture Type
                                picture_type =
                                    (pcs_ptr->idr_flag) ? I_SLICE :
                                    (pcs_ptr->cra_flag) ? I_SLICE :
                                    (pcs_ptr->pred_structure == SVT_AV1_PRED_LOW_DELAY_P) ? P_SLICE :
                                    (pcs_ptr->pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? B_SLICE :
                                    (pcs_ptr->pre_assignment_buffer_count == pcs_ptr->pred_struct_ptr->pred_struct_period) ? ((out_stride_diff64 == context_ptr->mini_gop_end_index[mini_gop_index] && 0) ? P_SLICE : B_SLICE) :

                                    (encode_context_ptr->pre_assignment_buffer_eos_flag) ? P_SLICE :
                                    B_SLICE;
                            }
                            if (!is_trailing_frame) {
                            // If mini GOP switch, reset position
                            encode_context_ptr->pred_struct_position = (pcs_ptr->init_pred_struct_position_flag) ?
                                pcs_ptr->pred_struct_ptr->init_pic_index :
                                encode_context_ptr->pred_struct_position;

                            // If Intra, reset position
                            if (pcs_ptr->idr_flag == TRUE)
                                encode_context_ptr->pred_struct_position = pcs_ptr->pred_struct_ptr->init_pic_index;
                            else if (pcs_ptr->cra_flag == TRUE && context_ptr->mini_gop_length[mini_gop_index] < pcs_ptr->pred_struct_ptr->pred_struct_period)
                                encode_context_ptr->pred_struct_position = pcs_ptr->pred_struct_ptr->init_pic_index;
                            else if (encode_context_ptr->elapsed_non_cra_count == 0) {
                                // If we are the picture directly after a CRA, we have to not use references that violate the CRA
                                encode_context_ptr->pred_struct_position = pcs_ptr->pred_struct_ptr->init_pic_index + 1;
                            }
                            // Elif Scene Change, determine leading and trailing pictures
                            //else if (encode_context_ptr->pre_assignment_buffer_scene_change_count > 0) {
                            //    if(buffer_index < encode_context_ptr->pre_assignment_buffer_scene_change_index) {
                            //        ++encode_context_ptr->pred_struct_position;
                            //        picture_type = P_SLICE;
                            //    }
                            //    else {
                            //        encode_context_ptr->pred_struct_position = pcs_ptr->pred_struct_ptr->init_pic_index + encode_context_ptr->pre_assignment_buffer_count - buffer_index - 1;
                            //    }
                            //}
                            // Else, Increment the position normally
                            else
                                ++encode_context_ptr->pred_struct_position;
                            }
                            // The poc number of the latest IDR picture is stored so that last_idr_picture (present in PCS) for the incoming pictures can be updated.
                            // The last_idr_picture is used in reseting the poc (in entropy coding) whenever IDR is encountered.
                            // Note IMP: This logic only works when display and decode order are the same. Currently for Random Access, IDR is inserted (similar to CRA) by using trailing P pictures (low delay fashion) and breaking prediction structure.
                            // Note: When leading P pictures are implemented, this logic has to change..
                            if (pcs_ptr->idr_flag == TRUE)
                                encode_context_ptr->last_idr_picture = pcs_ptr->picture_number;
                            else
                                pcs_ptr->last_idr_picture = encode_context_ptr->last_idr_picture;
                            if (!is_trailing_frame) {
                            // Cycle the PredStructPosition if its overflowed
                            encode_context_ptr->pred_struct_position = (encode_context_ptr->pred_struct_position == pcs_ptr->pred_struct_ptr->pred_struct_entry_count) ?
                                encode_context_ptr->pred_struct_position - pcs_ptr->pred_struct_ptr->pred_struct_period :
                                encode_context_ptr->pred_struct_position;
                            }

                            pred_position_ptr = pcs_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[encode_context_ptr->pred_struct_position];
                            if (scs_ptr->static_config.enable_overlays == TRUE) {
                                // At this stage we know the prediction structure and the location of ALT_REF pictures.
                                // For every ALTREF picture, there is an overlay picture. They extra pictures are released
                                // is_alt_ref flag is set for non-slice base layer pictures
                                if (pred_position_ptr->temporal_layer_index == 0 && picture_type != I_SLICE) {
                                    pcs_ptr->is_alt_ref = 1;
                                    pcs_ptr->frm_hdr.show_frame = 0;
                                    ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[out_stride_diff64 - 1]->object_ptr)->has_show_existing = FALSE;
                                }
                                // release the overlay PCS for non alt ref pictures. First picture does not have overlay PCS
                                else if (pcs_ptr->picture_number) {
                                    svt_release_object(pcs_ptr->overlay_ppcs_ptr->input_picture_wrapper_ptr);
                                    // release the pa_reference_picture
                                    svt_release_object(pcs_ptr->overlay_ppcs_ptr->pa_reference_picture_wrapper_ptr);
                                    // release the parent pcs
                                    // Note: this release will recycle ppcs to empty fifo if not live_count+1 in ResourceCoordination.
                                    svt_release_object(pcs_ptr->overlay_ppcs_ptr->p_pcs_wrapper_ptr);
                                    pcs_ptr->overlay_ppcs_ptr = NULL;
                                }
                            }
                            PictureParentControlSet       *cur_picture_control_set_ptr = pcs_ptr;
                            for (uint8_t loop_index = 0; loop_index <= pcs_ptr->is_alt_ref; loop_index++) {
                                // Set missing parts in the overlay  pictures
                                if (loop_index == 1) {
                                    pcs_ptr = pcs_ptr->overlay_ppcs_ptr;
                                    initialize_overlay_frame(pcs_ptr);
                                    picture_type = P_SLICE;
                                }
                                // Set the Slice type
                                pcs_ptr->slice_type = picture_type;

                                switch (picture_type) {
                                case I_SLICE:

                                    // Reset Prediction Structure Position & Reference Struct Position
                                    if (pcs_ptr->picture_number == 0)
                                        encode_context_ptr->intra_period_position = 0;
                                    encode_context_ptr->elapsed_non_cra_count = 0;

                                    //-------------------------------
                                    // IDR
                                    //-------------------------------
                                    if (pcs_ptr->idr_flag == TRUE) {
                                        // Set CRA flag
                                        pcs_ptr->cra_flag = FALSE;

                                        // Reset the pictures since last IDR counter
                                        encode_context_ptr->elapsed_non_idr_count = 0;
                                        //log latest key frame poc
                                        context_ptr->key_poc = pcs_ptr->picture_number;
                                    }
                                    //-------------------------------
                                    // CRA
                                    //-------------------------------
                                    else {
                                        // Set a Random Access Point if not an IDR
                                        pcs_ptr->cra_flag = TRUE;
                                    }
                                    break;

                                case P_SLICE:
                                case B_SLICE:

                                    // Reset CRA and IDR Flag
                                    pcs_ptr->cra_flag = FALSE;
                                    pcs_ptr->idr_flag = FALSE;

                                    if (loop_index == 0)
                                    {
                                        // Increment & Clip the elapsed Non-IDR Counter. This is clipped rather than allowed to free-run
                                        // inorder to avoid rollover issues.  This assumes that any the GOP period is less than MAX_ELAPSED_IDR_COUNT
                                        encode_context_ptr->elapsed_non_idr_count = MIN(encode_context_ptr->elapsed_non_idr_count + 1, MAX_ELAPSED_IDR_COUNT);
                                        encode_context_ptr->elapsed_non_cra_count = MIN(encode_context_ptr->elapsed_non_cra_count + 1, MAX_ELAPSED_IDR_COUNT);
                                    }

                                    CHECK_REPORT_ERROR(
                                        (pcs_ptr->pred_struct_ptr->pred_struct_entry_count < MAX_ELAPSED_IDR_COUNT),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR1);

                                    break;

                                default:

                                    CHECK_REPORT_ERROR_NC(
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR2);

                                    break;
                                }
                                pcs_ptr->pred_struct_index = (uint8_t)encode_context_ptr->pred_struct_position;
                                if (pcs_ptr->is_overlay) {
                                    // set the overlay frame as non reference frame with max temporal layer index
                                    pcs_ptr->temporal_layer_index = (uint8_t)pcs_ptr->hierarchical_levels;
                                    pcs_ptr->is_used_as_reference_flag = FALSE;
                                }
                                else {
                                    pcs_ptr->temporal_layer_index = (uint8_t)pred_position_ptr->temporal_layer_index;
                                    pcs_ptr->is_used_as_reference_flag = pred_position_ptr->is_referenced;
                                }

                                pre_assignment_buffer_first_pass_flag = FALSE;

                                // Film grain (assigning the random-seed)
                                {
                                    uint16_t *fgn_random_seed_ptr = &pcs_ptr->scs_ptr->film_grain_random_seed;
                                    pcs_ptr->frm_hdr.film_grain_params.random_seed = *fgn_random_seed_ptr;
                                    *fgn_random_seed_ptr += 3381;  // Changing random seed for film grain
                                    if (!(*fgn_random_seed_ptr))     // Random seed should not be zero
                                        *fgn_random_seed_ptr += 7391;
                                }
                                uint32_t pic_index = 0;
                                if (scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS) {
                                    pic_index = out_stride_diff64 - context_ptr->mini_gop_start_index[mini_gop_index];
                                } else {
                                    // For low delay P or low delay b case, get the the picture_index by mini_gop size
                                    if (scs_ptr->static_config.intra_period_length >= 0) {
                                        pic_index = (pcs_ptr->picture_number == 0) ? 0 :
                                            (uint32_t)(((pcs_ptr->picture_number - 1) % (scs_ptr->static_config.intra_period_length+1)) % pcs_ptr->pred_struct_ptr->pred_struct_period);
                                    } else {
                                        // intra-period=-1 case, no gop
                                        pic_index = (pcs_ptr->picture_number == 0) ? 0 :
                                            (uint32_t)((pcs_ptr->picture_number - 1) % pcs_ptr->pred_struct_ptr->pred_struct_period);
                                    }
                                }
                                if(scs_ptr->static_config.enable_manual_pred_struct){
                                    av1_generate_rps_ref_poc_from_user_config(pcs_ptr);
                                }
                                else{
                                    av1_generate_rps_info(
                                        pcs_ptr,
                                        encode_context_ptr,
                                        context_ptr,
                                        pic_index,
                                        mini_gop_index);
                                }

                                pcs_ptr->pic_index = pic_index;

                                pcs_ptr->allow_comp_inter_inter = 0;
                                pcs_ptr->is_skip_mode_allowed = 0;

                                pcs_ptr->frm_hdr.reference_mode = (ReferenceMode)0xFF;

                                if (pcs_ptr->slice_type != I_SLICE) {
                                    pcs_ptr->allow_comp_inter_inter = 1;
                                    if (pcs_ptr->slice_type == P_SLICE) {
                                        pcs_ptr->is_skip_mode_allowed = 0;
                                        pcs_ptr->frm_hdr.reference_mode = SINGLE_REFERENCE;
                                        pcs_ptr->skip_mode_flag = 0;
                                    }
                                    else if (pcs_ptr->temporal_layer_index == 0) {
                                        pcs_ptr->frm_hdr.reference_mode = REFERENCE_MODE_SELECT;
                                        pcs_ptr->frm_hdr.skip_mode_params.skip_mode_flag = 0;
                                    }
                                    else {
                                        pcs_ptr->frm_hdr.reference_mode = REFERENCE_MODE_SELECT;
                                        pcs_ptr->is_skip_mode_allowed = 1;
                                        pcs_ptr->skip_mode_flag = 1;
                                    }
                                }

                                pcs_ptr->av1_cm->mi_cols = pcs_ptr->aligned_width >> MI_SIZE_LOG2;
                                pcs_ptr->av1_cm->mi_rows = pcs_ptr->aligned_height >> MI_SIZE_LOG2;

                                //Jing: For low delay b/P case, don't alter the bias
                                memset(pcs_ptr->av1_cm->ref_frame_sign_bias, 0, 8 * sizeof(int32_t));
                                if (pcs_ptr->frm_hdr.reference_mode == REFERENCE_MODE_SELECT &&
                                        pcs_ptr->temporal_layer_index &&
                                        scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS)
                                {
                                    pcs_ptr->av1_cm->ref_frame_sign_bias[ALTREF_FRAME] =
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[ALTREF2_FRAME] =
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[BWDREF_FRAME] = 1;
                                }
                                if (pcs_ptr->frm_hdr.reference_mode == REFERENCE_MODE_SELECT &&
                                    pcs_ptr->temporal_layer_index &&
                                    scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS)
                                {
                                    if (pcs_ptr->av1_ref_signal.ref_poc_array[ALT] >= pcs_ptr->picture_number)
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[ALTREF_FRAME] = 1;
                                    else
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[ALTREF_FRAME] = 0;

                                    if (pcs_ptr->av1_ref_signal.ref_poc_array[ALT2] >= pcs_ptr->picture_number)
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[ALTREF2_FRAME] = 1;
                                    else
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[ALTREF2_FRAME] = 0;

                                    if (pcs_ptr->av1_ref_signal.ref_poc_array[BWD] >= pcs_ptr->picture_number)
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[BWDREF_FRAME] = 1;
                                    else
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[BWDREF_FRAME] = 0;
                                }

                                if (pcs_ptr->slice_type == I_SLICE) {
                                    // If running multi-threaded mode, perform SC detection in picture_analysis_kernel, else in picture_decision_kernel
                                    if (scs_ptr->static_config.logical_processors == 1) {
                                            int copy_frame = 1;
                                        if (pcs_ptr->scs_ptr->ipp_pass_ctrls.skip_frame_first_pass)
                                            copy_frame = (((pcs_ptr->picture_number % 8) == 0) || ((pcs_ptr->picture_number % 8) == 6) || ((pcs_ptr->picture_number % 8) == 7));
                                        // Bypass copy for the unecessary picture in IPPP pass
                                        if ((scs_ptr->static_config.pass != ENC_FIRST_PASS || copy_frame) == 0) {
                                            pcs_ptr->sc_class0 = pcs_ptr->sc_class1 = pcs_ptr->sc_class2 = 0;
                                        }
                                        else if (scs_ptr->static_config.screen_content_mode == 2) // auto detect
                                        {
                                                // SC Detection is OFF for 4K and higher
                                             if (scs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE)
                                                is_screen_content(pcs_ptr);
                                            else
                                                pcs_ptr->sc_class0 = pcs_ptr->sc_class1 = pcs_ptr->sc_class2 = 0;
                                        }
                                        else
                                            pcs_ptr->sc_class0 = pcs_ptr->sc_class1 = pcs_ptr->sc_class2 = scs_ptr->static_config.screen_content_mode;
                                    }
                                    context_ptr->last_i_picture_sc_class0 = pcs_ptr->sc_class0;
                                    context_ptr->last_i_picture_sc_class1 = pcs_ptr->sc_class1;
                                    context_ptr->last_i_picture_sc_class2 = pcs_ptr->sc_class2;
                                }
                                else {
                                    pcs_ptr->sc_class0 = context_ptr->last_i_picture_sc_class0;
                                    pcs_ptr->sc_class1 = context_ptr->last_i_picture_sc_class1;
                                    pcs_ptr->sc_class2 = context_ptr->last_i_picture_sc_class2;
                                }
                                // TODO: put this in EbMotionEstimationProcess?
                                copy_tf_params(scs_ptr, pcs_ptr);
                                // TODO: put this in EbMotionEstimationProcess?
                                // ME Kernel Multi-Processes Signal(s) derivation
                                if (scs_ptr->static_config.pass == ENC_FIRST_PASS)
                                    first_pass_signal_derivation_multi_processes(scs_ptr, pcs_ptr);
                                else
                                    signal_derivation_multi_processes_oq(scs_ptr, pcs_ptr, context_ptr);
                                // Update the Dependant List Count - If there was an I-frame or Scene Change or S-frame, then cleanup the Picture Decision PA Reference Queue Dependent Counts
                                if (pcs_ptr->slice_type == I_SLICE || pcs_ptr->frm_hdr.frame_type == S_FRAME)
                                {
                                    input_queue_index = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                                    while (input_queue_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
                                        input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[input_queue_index];

                                        int32_t diff_n = 0;
                                        // Modify Dependent List0
                                        dep_list_count = input_entry_ptr->list0.list_count;
                                        for (dep_idx = 0; dep_idx < dep_list_count; ++dep_idx) {
                                            // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                                            // current_input_poc += (current_input_poc < input_entry_ptr->pocNumber) ? (1 << scs_ptr->bits_for_picture_order_count) : 0;

                                            dep_poc = POC_CIRCULAR_ADD(
                                                (int64_t)input_entry_ptr->picture_number, // can't use a value that gets reset
                                                input_entry_ptr->list0.list[dep_idx]/*,
                                                scs_ptr->bits_for_picture_order_count*/);

                                            // If Dependent POC is greater or equal to the IDR POC
                                            Bool cleanup_dep = FALSE;
                                            if (pcs_ptr->slice_type == I_SLICE) {
                                                if (dep_poc >= (int64_t)pcs_ptr->picture_number && input_entry_ptr->list0.list[dep_idx]) {
                                                    cleanup_dep = TRUE;
                                                }
                                            }
                                            else if (dep_poc > (int64_t)pcs_ptr->picture_number && input_entry_ptr->list0.list[dep_idx]) {  // s-frame
                                                // Reason for why using condition 'dep_poc > (int64_t)pcs_ptr->picture_number', but not '>=' like I-SLICE case:
                                                // s-frame is allowed to use forward reference frames while I_SLICE has no reference frame.
                                                cleanup_dep = TRUE;
                                            }
                                            if (cleanup_dep) {
                                                input_entry_ptr->list0.list[dep_idx] = 0;

                                                // Decrement the Reference's referenceCount
                                                --input_entry_ptr->dependent_count;

                                                diff_n--;

                                                CHECK_REPORT_ERROR(
                                                    (input_entry_ptr->dependent_count != ~0u),
                                                    encode_context_ptr->app_callback_ptr,
                                                    EB_ENC_PD_ERROR3);
                                            }
                                        }

                                        // Modify Dependent List1
                                        dep_list_count = input_entry_ptr->list1.list_count;
                                        for (dep_idx = 0; dep_idx < dep_list_count; ++dep_idx) {
                                            // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                                            // current_input_poc += (current_input_poc < input_entry_ptr->pocNumber) ? (1 << scs_ptr->bits_for_picture_order_count) : 0;

                                            dep_poc = POC_CIRCULAR_ADD(
                                                (int64_t)input_entry_ptr->picture_number,
                                                input_entry_ptr->list1.list[dep_idx]/*,
                                                scs_ptr->bits_for_picture_order_count*/);

                                            // If Dependent POC is greater or equal to the IDR POC
                                            Bool cleanup_dep = FALSE;
                                            if (pcs_ptr->slice_type == I_SLICE) {
                                                if (((dep_poc >= (int64_t)pcs_ptr->picture_number) || (((pcs_ptr->pre_assignment_buffer_count != pcs_ptr->pred_struct_ptr->pred_struct_period) || (pcs_ptr->idr_flag == TRUE)) && (dep_poc > ((int64_t)pcs_ptr->picture_number - (int64_t)pcs_ptr->pre_assignment_buffer_count)))) && input_entry_ptr->list1.list[dep_idx]) {
                                                    cleanup_dep = TRUE;
                                                }
                                            }
                                            else if (dep_poc > (int64_t)pcs_ptr->picture_number && input_entry_ptr->list1.list[dep_idx]) {  // s-frame
                                                cleanup_dep = TRUE;
                                            }

                                            if (cleanup_dep)
                                            {
                                                input_entry_ptr->list1.list[dep_idx] = 0;

                                                // Decrement the Reference's referenceCount
                                                --input_entry_ptr->dependent_count;

                                                diff_n--;

                                                CHECK_REPORT_ERROR(
                                                    (input_entry_ptr->dependent_count != ~0u),
                                                    encode_context_ptr->app_callback_ptr,
                                                    EB_ENC_PD_ERROR3);
                                            }
                                        }

                                        if (diff_n) {

                                            PictureParentControlSet * pcs_of_this_entry =
                                                is_pic_still_in_pre_assign_buffer(
                                                    encode_context_ptr,
                                                    context_ptr,
                                                    mini_gop_index,
                                                    input_entry_ptr->picture_number);

                                            if (pcs_of_this_entry != NULL) {
                                               //The reference picture that needs a cut link is still in PD, it can do self clean-up later in picMgr.
                                               //the dec order of this entry is usually > current pcs; so  if curr goes first to picMgr,
                                               //then it will not find this entry in picMgr ref Q to clear it; so the ref has to do a self cleaning
                                                pcs_of_this_entry->self_updated_links = diff_n;

                                            }
                                            else {
                                               //picture has left PD, so pcs_ptr has to do the cleanup later in PicMgr (and should find it in picMgr ref Q)
                                                //TODO: better to use first pic in MG as triggering picture
                                               pcs_ptr->updated_links_arr[pcs_ptr->other_updated_links_cnt].pic_num = input_entry_ptr->picture_number;
                                               pcs_ptr->updated_links_arr[pcs_ptr->other_updated_links_cnt++].dep_cnt_diff = diff_n;

                                            }
                                        }
                                        // Increment the input_queue_index Iterator
                                        input_queue_index = (input_queue_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : input_queue_index + 1;
                                    }
                                }
                                else if (pcs_ptr->idr_flag == TRUE) {
                                    // Set the Picture Decision PA Reference Entry pointer
                                    input_entry_ptr = (PaReferenceQueueEntry*)NULL;
                                }

                                // Place Picture in Picture Decision PA Reference Queue
                                // there is no need to add the overlay frames in the PA Reference Queue
                                if (!pcs_ptr->is_overlay)
                                {
                                    input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_tail_index];
                                    input_entry_ptr->input_object_ptr = pcs_ptr->pa_reference_picture_wrapper_ptr;
                                    input_entry_ptr->picture_number = pcs_ptr->picture_number;
                                    input_entry_ptr->is_alt_ref = pcs_ptr->is_alt_ref;
                                    input_entry_ptr->eb_y8b_wrapper_ptr = pcs_ptr->eb_y8b_wrapper_ptr;
                                    encode_context_ptr->picture_decision_pa_reference_queue_tail_index =
                                        (encode_context_ptr->picture_decision_pa_reference_queue_tail_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->picture_decision_pa_reference_queue_tail_index + 1;

                                    // Check if the Picture Decision PA Reference is full
                                    CHECK_REPORT_ERROR(
                                        (((encode_context_ptr->picture_decision_pa_reference_queue_head_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) ||
                                        (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->input_object_ptr == NULL))),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR4);
                                }
                                // Copy the reference lists into the inputEntry and
                                // set the Reference Counts Based on Temporal Layer and how many frames are active
                                pcs_ptr->ref_list0_count = (picture_type == I_SLICE) ? 0 :
                                                                            (pcs_ptr->is_overlay) ? 1 : (uint8_t)pred_position_ptr->ref_list0.reference_list_count;
                                pcs_ptr->ref_list1_count = (picture_type == I_SLICE || pcs_ptr->is_overlay) ? 0 : (uint8_t)pred_position_ptr->ref_list1.reference_list_count;

                                update_count_try(scs_ptr, pcs_ptr);

                                if (context_ptr->transition_detected == 1) {
                                    if (pcs_ptr->slice_type == P_SLICE) {
                                        pcs_ptr->transition_present = 1;
                                    }
                                    else if (pcs_ptr->temporal_layer_index == 0) {
                                        pcs_ptr->transition_present = 1;
                                        context_ptr->transition_detected = 0;
                                    }
                                }
                                if (picture_type == B_SLICE && pcs_ptr->temporal_layer_index == 0 && pcs_ptr->list0_only_base_ctrls.enabled) {
                                    update_list0_only_base(scs_ptr, pcs_ptr);
                                }
                                assert(pcs_ptr->ref_list0_count_try <= pcs_ptr->ref_list0_count);
                                assert(pcs_ptr->ref_list1_count_try <= pcs_ptr->ref_list1_count);
                                if (!pcs_ptr->is_overlay) {
                                    input_entry_ptr->list0_ptr = &pred_position_ptr->ref_list0;
                                    input_entry_ptr->list1_ptr = &pred_position_ptr->ref_list1;
                                    // Copy the Dependent Lists
                                    // *Note - we are removing any leading picture dependencies for now
                                    input_entry_ptr->list0.list_count = 0;
                                    for (dep_idx = 0; dep_idx < pred_position_ptr->dep_list0.list_count; ++dep_idx) {
                                        if (pred_position_ptr->dep_list0.list[dep_idx] >= 0)
                                            input_entry_ptr->list0.list[input_entry_ptr->list0.list_count++] = pred_position_ptr->dep_list0.list[dep_idx];
                                    }
                                    input_entry_ptr->list1.list_count = pred_position_ptr->dep_list1.list_count;
                                    for (dep_idx = 0; dep_idx < pred_position_ptr->dep_list1.list_count; ++dep_idx)
                                        input_entry_ptr->list1.list[dep_idx] = pred_position_ptr->dep_list1.list[dep_idx];
                                    // add the overlay picture to the dependant list
                                    input_entry_ptr->dep_list0_count = (pcs_ptr->is_alt_ref) ? input_entry_ptr->list0.list_count + 1 : input_entry_ptr->list0.list_count;

                                    input_entry_ptr->dep_list1_count = input_entry_ptr->list1.list_count;
                                    input_entry_ptr->dependent_count = input_entry_ptr->dep_list0_count + input_entry_ptr->dep_list1_count;


                                }

                                CHECK_REPORT_ERROR(
                                    (pcs_ptr->pred_struct_ptr->pred_struct_period * REF_LIST_MAX_DEPTH < MAX_ELAPSED_IDR_COUNT),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PD_ERROR5);

                                // Reset the PA Reference Lists
                                EB_MEMSET(pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                                EB_MEMSET(pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                                EB_MEMSET(pcs_ptr->ref_y8b_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                                EB_MEMSET(pcs_ptr->ref_y8b_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                                EB_MEMSET(pcs_ptr->ref_pic_poc_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                                EB_MEMSET(pcs_ptr->ref_pic_poc_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));

                                // Update the RC param queue
                                update_rc_param_queue(pcs_ptr, encode_context_ptr);
                            }
                            pcs_ptr = cur_picture_control_set_ptr;
                        }
                        if(scs_ptr->static_config.enable_manual_pred_struct){
                            EbErrorType ret = av1_generate_minigop_rps_info_from_user_config(encode_context_ptr,context_ptr,mini_gop_index);
                            if (ret != EB_ErrorNone) {
                                CHECK_REPORT_ERROR_NC(
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PD_ERROR9);
                            }
                        }
                        // Add 1 to the loop for the overlay picture. If the last picture is alt ref, increase the loop by 1 to add the overlay picture
                        uint32_t has_overlay = ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_end_index[mini_gop_index]]->object_ptr)->is_alt_ref ? 1 : 0;
                        for (out_stride_diff64 = context_ptr->mini_gop_start_index[mini_gop_index]; out_stride_diff64 <= context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay; ++out_stride_diff64) {
                            // 2nd Loop over Pictures in the Pre-Assignment Buffer
                            // Assign the overlay pcs. Since Overlay picture is not added to the picture_decision_pa_reference_queue, in the next stage, the loop finds the alt_ref picture. The reference for overlay frame is hardcoded later
                            if (has_overlay && out_stride_diff64 == context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay) {
                                pcs_ptr = ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_end_index[mini_gop_index]]->object_ptr)->overlay_ppcs_ptr;
                            }
                            else {
                                pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[out_stride_diff64]->object_ptr;

                                if (scs_ptr->static_config.sframe_dist != 0 || !pcs_ptr->is_not_scaled) {
                                    update_sframe_ref_order_hint(pcs_ptr, context_ptr);
                                }
                            }
                            pcs_ptr->picture_number_alt = encode_context_ptr->picture_number_alt++;

                            // Set the Decode Order
                            if ((context_ptr->mini_gop_idr_count[mini_gop_index] == 0) &&
                                (context_ptr->mini_gop_length[mini_gop_index] == pcs_ptr->pred_struct_ptr->pred_struct_period) &&
                                (scs_ptr->static_config.pred_structure == SVT_AV1_PRED_RANDOM_ACCESS) &&
                                !pcs_ptr->is_overlay) {
                                pcs_ptr->decode_order = encode_context_ptr->decode_base_number + pcs_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[pcs_ptr->pred_struct_index]->decode_order;
                            }
                            else
                                pcs_ptr->decode_order = pcs_ptr->picture_number_alt;

                            if (pcs_ptr->end_of_sequence_flag == TRUE) {
                                encode_context_ptr->terminating_sequence_flag_received = TRUE;
                                encode_context_ptr->terminating_picture_number = pcs_ptr->picture_number_alt;
                            }

                            // Find the Reference in the Picture Decision PA Reference Queue
                            input_queue_index = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                            do {
                                // Setup the Picture Decision PA Reference Queue Entry
                                input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[input_queue_index];

                                // Increment the reference_queue_index Iterator
                                input_queue_index = (input_queue_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : input_queue_index + 1;
                            } while ((input_queue_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) && (input_entry_ptr->picture_number != pcs_ptr->picture_number));

                            EB_MEMSET(pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                            EB_MEMSET(pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                            EB_MEMSET(pcs_ptr->ref_y8b_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                            EB_MEMSET(pcs_ptr->ref_y8b_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                            EB_MEMSET(pcs_ptr->ref_pic_poc_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                            EB_MEMSET(pcs_ptr->ref_pic_poc_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                            CHECK_REPORT_ERROR(
                                (input_entry_ptr->picture_number == pcs_ptr->picture_number),
                                encode_context_ptr->app_callback_ptr,
                                EB_ENC_PD_ERROR7);

                            // Reset the PA Reference Lists
                            EB_MEMSET(pcs_ptr->ref_pa_pic_ptr_array, 0, 2 * sizeof(EbObjectWrapper*));

                            EB_MEMSET(pcs_ptr->ref_pic_poc_array, 0, 2 * sizeof(uint64_t));

                            // Configure List0
                            if ((pcs_ptr->slice_type == P_SLICE) || (pcs_ptr->slice_type == B_SLICE)) {
                                uint8_t ref_pic_index;
                                for (ref_pic_index = 0; ref_pic_index < pcs_ptr->ref_list0_count; ++ref_pic_index) {
                                    if (pcs_ptr->ref_list0_count) {
                                        if (pcs_ptr->is_overlay) {
                                        // hardcode the reference for the overlay frame
                                            ref_poc = pcs_ptr->picture_number;
                                        } else {
                                            ref_poc = POC_CIRCULAR_ADD(
                                                (int64_t)pcs_ptr->picture_number,
                                                -input_entry_ptr->list0_ptr->reference_list[ref_pic_index]);
                                        }

                                        pa_reference_entry_ptr = search_ref_in_ref_queue_pa(encode_context_ptr, ref_poc);
                                        assert(pa_reference_entry_ptr != 0);
                                        CHECK_REPORT_ERROR((pa_reference_entry_ptr),
                                            encode_context_ptr->app_callback_ptr,
                                            EB_ENC_PM_ERROR10);
                                            // Set the Reference Object
                                        pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_0][ref_pic_index] = pa_reference_entry_ptr->input_object_ptr;
                                        pcs_ptr->ref_pic_poc_array[REF_LIST_0][ref_pic_index] = ref_poc;
                                        // Increment the PA Reference's liveCount by the number of tiles in the input picture
                                        //assert((int32_t)pa_reference_entry_ptr->input_object_ptr->live_count > 0);
                                        svt_object_inc_live_count(
                                            pa_reference_entry_ptr->input_object_ptr,
                                            1);

                                        pcs_ptr->ref_y8b_array[REF_LIST_0][ref_pic_index] = pa_reference_entry_ptr->eb_y8b_wrapper_ptr;

                                        if (pa_reference_entry_ptr->eb_y8b_wrapper_ptr) {
                                            //y8b follows longest life cycle of pa ref and input. so it needs to build on top of live count of pa ref
                                            svt_object_inc_live_count(
                                                pa_reference_entry_ptr->eb_y8b_wrapper_ptr,
                                                1);
                                        }
                                        --pa_reference_entry_ptr->dependent_count;
                                    }
                                }
                            }

                            // Configure List1
                            if (pcs_ptr->slice_type == B_SLICE) {
                                uint8_t ref_pic_index;
                                for (ref_pic_index = 0; ref_pic_index < pcs_ptr->ref_list1_count; ++ref_pic_index) {
                                    if (pcs_ptr->ref_list1_count) {
                                        ref_poc = POC_CIRCULAR_ADD(
                                            (int64_t)pcs_ptr->picture_number,
                                            -input_entry_ptr->list1_ptr->reference_list[ref_pic_index]);

                                        pa_reference_entry_ptr = search_ref_in_ref_queue_pa(encode_context_ptr, ref_poc);

                                        assert(pa_reference_entry_ptr != 0);
                                        CHECK_REPORT_ERROR((pa_reference_entry_ptr),
                                            encode_context_ptr->app_callback_ptr,
                                            EB_ENC_PM_ERROR10);
                                        // Set the Reference Object
                                        pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_1][ref_pic_index] = pa_reference_entry_ptr->input_object_ptr;
                                        pcs_ptr->ref_pic_poc_array[REF_LIST_1][ref_pic_index] = ref_poc;

                                        // Increment the PA Reference's liveCount by the number of tiles in the input picture
                                        //assert((int32_t)pa_reference_entry_ptr->input_object_ptr->live_count > 0);
                                        svt_object_inc_live_count(
                                            pa_reference_entry_ptr->input_object_ptr,
                                            1);
                                        pcs_ptr->ref_y8b_array[REF_LIST_1][ref_pic_index] = pa_reference_entry_ptr->eb_y8b_wrapper_ptr;

                                        if (pa_reference_entry_ptr->eb_y8b_wrapper_ptr) {
                                            //y8b follows longest life cycle of pa ref and input. so it needs to build on top of live count of pa ref
                                            svt_object_inc_live_count(
                                                pa_reference_entry_ptr->eb_y8b_wrapper_ptr,
                                                1);
                                        }
                                        --pa_reference_entry_ptr->dependent_count;
                                    }
                                }
                            }

                            svt_av1_setup_skip_mode_allowed(pcs_ptr);

                            pcs_ptr->is_skip_mode_allowed = pcs_ptr->frm_hdr.skip_mode_params.skip_mode_allowed;
                            pcs_ptr->skip_mode_flag = pcs_ptr->is_skip_mode_allowed;
                            //SVT_LOG("POC:%i  skip_mode_allowed:%i  REF_SKIP_0: %i   REF_SKIP_1: %i \n",pcs_ptr->picture_number, pcs_ptr->skip_mode_info.skip_mode_allowed, pcs_ptr->skip_mode_info.ref_frame_idx_0, pcs_ptr->skip_mode_info.ref_frame_idx_1);

                            {
                                if (pcs_ptr->ref_list0_count)
                                    pcs_ptr->scene_transition_flag[REF_LIST_0] = FALSE;
                                if (pcs_ptr->ref_list1_count)
                                    pcs_ptr->scene_transition_flag[REF_LIST_1] = FALSE;
                            }

                            //set the ref frame types used for this picture,
                            set_all_ref_frame_type(pcs_ptr, pcs_ptr->ref_frame_type_arr, &pcs_ptr->tot_ref_frame_types);

                            pcs_ptr->me_processed_b64_count = 0;

                            uint32_t pic_it = out_stride_diff64 - context_ptr->mini_gop_start_index[mini_gop_index];
                            context_ptr->mg_pictures_array[pic_it] = pcs_ptr;
                            if (out_stride_diff64 == context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay) {
                                // Increment the Decode Base Number
                                encode_context_ptr->decode_base_number += context_ptr->mini_gop_length[mini_gop_index] + has_overlay;
                            }
                        }
                        uint32_t mg_size = context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay - context_ptr->mini_gop_start_index[mini_gop_index]+1;
                        context_ptr->mg_size = mg_size;
                        EB_MEMCPY(context_ptr->mg_pictures_array_disp_order, context_ptr->mg_pictures_array, mg_size * sizeof(PictureParentControlSet*));
                        //sort based on decoder order
                        {
                            uint32_t  num_of_cand_to_sort = mg_size;
                            uint32_t i, j;
                            for (i = 0; i < num_of_cand_to_sort - 1; ++i) {
                                for (j = i + 1; j < num_of_cand_to_sort; ++j) {
                                    if (context_ptr->mg_pictures_array[j]->decode_order < context_ptr->mg_pictures_array[i]->decode_order) {
                                        PictureParentControlSet* temp = context_ptr->mg_pictures_array[i];
                                        context_ptr->mg_pictures_array[i] = context_ptr->mg_pictures_array[j];
                                        context_ptr->mg_pictures_array[j] = temp;
                                    }
                                }
                            }
                        }
                        for (uint32_t pic_i = 0; pic_i < mg_size; ++pic_i) {
                            if (pic_i == 0)
                                context_ptr->mg_pictures_array_disp_order[pic_i]->first_frame_in_minigop = 1;
                            else
                                context_ptr->mg_pictures_array_disp_order[pic_i]->first_frame_in_minigop = 0;
                        }

                        for (uint32_t pic_i = 0; pic_i < mg_size; ++pic_i) {
                            pcs_ptr = context_ptr->mg_pictures_array_disp_order[pic_i];
                            set_gf_group_param(pcs_ptr);
                        }

                        //Process previous delayed Intra if we have one
                        pcs_ptr->is_new_gf_group = 0;
                        if (context_ptr->prev_delayed_intra) {
                            pcs_ptr = context_ptr->prev_delayed_intra;
                            store_gf_group(pcs_ptr, context_ptr, mg_size);
                        }
                        else {
                            for (uint32_t pic_i = 0; pic_i < mg_size; ++pic_i) {
                                pcs_ptr = context_ptr->mg_pictures_array_disp_order[pic_i];
                                if (is_delayed_intra(pcs_ptr) == FALSE) {
                                    store_gf_group(pcs_ptr, context_ptr, mg_size);
                                }
                            }
                        }
                        //Process previous delayed Intra if we have one
                        if (context_ptr->prev_delayed_intra) {
                            pcs_ptr = context_ptr->prev_delayed_intra;
                            mctf_frame(scs_ptr, pcs_ptr, context_ptr, out_stride_diff64);
                        }

                        //Do TF loop in display order
                        for (uint32_t pic_i = 0; pic_i < mg_size; ++pic_i) {
                            pcs_ptr = context_ptr->mg_pictures_array_disp_order[pic_i];

                            if (is_delayed_intra(pcs_ptr) == FALSE) {
                                mctf_frame(scs_ptr, pcs_ptr, context_ptr, out_stride_diff64);
                            }
                        }

                        if (context_ptr->prev_delayed_intra) {
                            pcs_ptr = context_ptr->prev_delayed_intra;
                            context_ptr->prev_delayed_intra = NULL;
                            send_picture_out(scs_ptr, pcs_ptr, context_ptr);
                        }

                       //split MG into two for these two special cases
                       uint8_t ldp_delayi_mg = 0;
                       uint8_t ldp_i_eos_mg = 0;
                       if (pcs_ptr->pred_struct_ptr->pred_type == SVT_AV1_PRED_RANDOM_ACCESS &&
                           context_ptr->mg_pictures_array[0]->slice_type == P_SLICE) {
                           if (is_delayed_intra(context_ptr->mg_pictures_array[mg_size - 1]))
                               ldp_delayi_mg = 1;
                           else if (context_ptr->mg_pictures_array[mg_size - 1]->slice_type == I_SLICE &&
                               context_ptr->mg_pictures_array[mg_size - 1]->end_of_sequence_flag)
                               ldp_i_eos_mg = 1;
                       }


                        for (uint32_t pic_i = 0; pic_i < mg_size; ++pic_i){

                            pcs_ptr = context_ptr->mg_pictures_array[pic_i];
                            if (is_delayed_intra(pcs_ptr)) {
                                context_ptr->prev_delayed_intra = pcs_ptr;

                                if(ldp_delayi_mg )
                                   context_ptr->mg_progress_id++;

                                pcs_ptr->ext_mg_id = context_ptr->mg_progress_id;
                                pcs_ptr->ext_mg_size = 1;

                            }else{
                                pcs_ptr->ext_mg_id = context_ptr->mg_progress_id;
                                pcs_ptr->ext_mg_size = ldp_delayi_mg ? mg_size-1 : mg_size;

                                if (ldp_i_eos_mg) {
                                    if (pcs_ptr->slice_type == P_SLICE) {
                                        pcs_ptr->ext_mg_size = mg_size - 1;
                                    }
                                    else if (pcs_ptr->slice_type == I_SLICE){
                                        context_ptr->mg_progress_id++;
                                        pcs_ptr->ext_mg_id = context_ptr->mg_progress_id;
                                        pcs_ptr->ext_mg_size =  1 ;
                                    }
                                }


                                send_picture_out(scs_ptr, pcs_ptr, context_ptr);
                            }


                        }

                        context_ptr->mg_progress_id++;

                    } // End MINI GOPs loop
                    // Reset the Pre-Assignment Buffer
                    encode_context_ptr->pre_assignment_buffer_count = 0;
                    encode_context_ptr->pre_assignment_buffer_idr_count = 0;
                    encode_context_ptr->pre_assignment_buffer_intra_count = 0;
                    encode_context_ptr->pre_assignment_buffer_scene_change_count = 0;
                    encode_context_ptr->pre_assignment_buffer_eos_flag = FALSE;
                }

                // Walk the picture_decision_pa_reference_queue and remove entries that have been completely referenced.
                input_queue_index = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                while (input_queue_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
                    input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[input_queue_index];

                    // Remove the entry
                    if ((input_entry_ptr->dependent_count == 0) &&
                        (input_entry_ptr->input_object_ptr)) {
                        // Release the nominal live_count value
                        //assert((int32_t)input_entry_ptr->input_object_ptr->live_count > 0);
                        svt_release_object(input_entry_ptr->input_object_ptr);

                        if (input_entry_ptr->eb_y8b_wrapper_ptr) {
                            //y8b needs to get decremented at the same time of pa ref
                            svt_release_object(input_entry_ptr->eb_y8b_wrapper_ptr);
                        }

                        input_entry_ptr->input_object_ptr = (EbObjectWrapper*)NULL;
                    }

                    // Increment the head_index if the head is null
                    encode_context_ptr->picture_decision_pa_reference_queue_head_index =
                        (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->input_object_ptr) ? encode_context_ptr->picture_decision_pa_reference_queue_head_index :
                        (encode_context_ptr->picture_decision_pa_reference_queue_head_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0
                        : encode_context_ptr->picture_decision_pa_reference_queue_head_index + 1;

                    CHECK_REPORT_ERROR(
                        (((encode_context_ptr->picture_decision_pa_reference_queue_head_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) ||
                        (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->input_object_ptr == NULL))),
                        encode_context_ptr->app_callback_ptr,
                        EB_ENC_PD_ERROR4);

                    // Increment the input_queue_index Iterator
                    input_queue_index = (input_queue_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : input_queue_index + 1;
                }

                // Increment the Picture Decision Reordering Queue Head Ptr
                encode_context_ptr->picture_decision_reorder_queue_head_index = (encode_context_ptr->picture_decision_reorder_queue_head_index == PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->picture_decision_reorder_queue_head_index + 1;

                // Get the next entry from the Picture Decision Reordering Queue (Entry N+1)
                queue_entry_ptr = encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index];
            }
            else
                break;
        }

        if (scs_ptr->static_config.enable_overlays == TRUE) {
            // release ppcs, since live_count + 1 before post in ResourceCoordination
            svt_release_object(in_results_ptr->pcs_wrapper_ptr);
        }

        // Release the Input Results
        svt_release_object(in_results_wrapper_ptr);
    }

    return NULL;
}
// clang-format on
