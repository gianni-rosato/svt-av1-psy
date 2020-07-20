// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
#if NOISE_BASED_TF_FRAMES
#include "EbMalloc.h"
#endif

/************************************************
 * Defines
 ************************************************/
#define  LAY0_OFF  0
#define  LAY1_OFF  3
#define  LAY2_OFF  5
#define  LAY3_OFF  6
#define  LAY4_OFF  7
extern PredictionStructureConfigEntry flat_pred_struct[];
extern PredictionStructureConfigEntry two_level_hierarchical_pred_struct[];
extern PredictionStructureConfigEntry three_level_hierarchical_pred_struct[];
extern PredictionStructureConfigEntry four_level_hierarchical_pred_struct[];
extern PredictionStructureConfigEntry five_level_hierarchical_pred_struct[];
extern PredictionStructureConfigEntry six_level_hierarchical_pred_struct[];
#if TF_LEVELS
typedef struct  TfControls {
    uint8_t enabled;
    uint8_t window_size;
#if IMPROVED_TF_LEVELS
    uint8_t noise_based_window_adjust;
#endif
}TfControls;
#endif

/**************************************
 * Context
 **************************************/
typedef struct PictureDecisionContext
{
    EbDctor      dctor;
    EbFifo       *picture_analysis_results_input_fifo_ptr;
    EbFifo       *picture_decision_results_output_fifo_ptr;
#if DECOUPLE_ME_RES
    EbFifo       *me_fifo_ptr;
#endif
    uint64_t      last_solid_color_frame_poc;

    EbBool        reset_running_avg;

    uint32_t    **ahd_running_avg_cb;
    uint32_t    **ahd_running_avg_cr;
    uint32_t    **ahd_running_avg;
    EbBool        is_scene_change_detected;

    // Dynamic GOP
    uint32_t      ttl_region_activity_cost[MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT];

    uint32_t      total_number_of_mini_gops;

    uint32_t      mini_gop_start_index[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_end_index[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_length[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_intra_count[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_idr_count[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t      mini_gop_hierarchical_levels[MINI_GOP_WINDOW_MAX_COUNT];
    EbBool        mini_gop_activity_array[MINI_GOP_MAX_COUNT];
    uint32_t      mini_gop_region_activity_cost_array[MINI_GOP_MAX_COUNT][MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT];

    uint32_t      mini_gop_group_faded_in_pictures_count[MINI_GOP_MAX_COUNT];
    uint32_t      mini_gop_group_faded_out_pictures_count[MINI_GOP_MAX_COUNT];
    uint8_t       lay0_toggle; //3 way toggle 0->1->2
    uint8_t       lay1_toggle; //2 way toggle 0->1
    uint8_t       lay2_toggle; //2 way toggle 0->1
    EbBool        mini_gop_toggle;    //mini GOP toggling since last Key Frame  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
    uint8_t       last_i_picture_sc_detection;
    uint64_t      key_poc;
#if TF_LEVELS
    uint8_t tf_level;
    TfControls tf_ctrls;
#endif
#if DECOUPLE_ME_RES
    PictureParentControlSet* mg_pictures_array[1<<MAX_TEMPORAL_LAYERS];
    DepCntPicInfo updated_links_arr[UPDATED_LINKS];//if not empty, this picture is a depn-cnt-cleanUp triggering picture (I frame; or MG size change )
                                                      //this array will store all others pictures needing a dep-cnt clean up.
    uint32_t other_updated_links_cnt; //how many other pictures in the above array needing a dep-cnt clean-up
#endif
} PictureDecisionContext;

void init_resize_picture(SequenceControlSet* scs_ptr, PictureParentControlSet* pcs_ptr);

uint64_t  get_ref_poc(PictureDecisionContext *context, uint64_t curr_picture_number, int32_t delta_poc)
{
    uint64_t ref_poc;

    if ((int64_t)curr_picture_number - (int64_t)delta_poc < (int64_t)context->key_poc)
        ref_poc = context->key_poc;
    else
        ref_poc = curr_picture_number - delta_poc;

    return ref_poc;
}

typedef struct {
    MvReferenceFrame ref_type;
    int used;
    uint64_t poc;
} RefFrameInfo;

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

void eb_av1_setup_skip_mode_allowed(PictureParentControlSet  *parent_pcs_ptr) {

    FrameHeader *frm_hdr = &parent_pcs_ptr->frm_hdr;

    RefFrameInfo ref_frame_arr_single[7];

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
#if NOISE_BASED_TF_FRAMES
#define FUTURE_WINDOW_WIDTH                 12
#else
#define FUTURE_WINDOW_WIDTH                 6
#endif
#define FLASH_TH                            5
#define FADE_TH                             3
#define SCENE_TH                            3000
#define NOISY_SCENE_TH                      4500    // SCD TH in presence of noise
#define HIGH_PICTURE_VARIANCE_TH            1500
#define NUM64x64INPIC(w,h)          ((w*h)>> (eb_log2f(BLOCK_SIZE_64)<<1))
#define QUEUE_GET_PREVIOUS_SPOT(h)  ((h == 0) ? PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1 : h - 1)
#define QUEUE_GET_NEXT_SPOT(h,off)  (( (h+off) >= PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH) ? h+off - PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH  : h + off)

#define WTH 64
#define OTH 64
static void picture_decision_context_dctor(EbPtr p)
{
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    PictureDecisionContext* obj = (PictureDecisionContext*)thread_context_ptr->priv;

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
    const EbEncHandle   *enc_handle_ptr)
{
    uint32_t arr_row, arr_col;

    PictureDecisionContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv = context_ptr;
    thread_context_ptr->dctor = picture_decision_context_dctor;

    context_ptr->picture_analysis_results_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->picture_analysis_results_resource_ptr, 0);
    context_ptr->picture_decision_results_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->picture_decision_results_resource_ptr, 0);

    EB_MALLOC_2D(context_ptr->ahd_running_avg_cb,  MAX_NUMBER_OF_REGIONS_IN_WIDTH, MAX_NUMBER_OF_REGIONS_IN_HEIGHT);
    EB_MALLOC_2D(context_ptr->ahd_running_avg_cr, MAX_NUMBER_OF_REGIONS_IN_WIDTH, MAX_NUMBER_OF_REGIONS_IN_HEIGHT);
    EB_MALLOC_2D(context_ptr->ahd_running_avg, MAX_NUMBER_OF_REGIONS_IN_WIDTH, MAX_NUMBER_OF_REGIONS_IN_HEIGHT);

    for (arr_row = 0; arr_row < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; arr_row++)
    {
        for (arr_col = 0; arr_col < MAX_NUMBER_OF_REGIONS_IN_WIDTH; arr_col++) {
            context_ptr->ahd_running_avg_cb[arr_col][arr_row] = 0;
            context_ptr->ahd_running_avg_cr[arr_col][arr_row] = 0;
            context_ptr->ahd_running_avg[arr_col][arr_row] = 0;
        }
    }

    context_ptr->reset_running_avg = EB_TRUE;
#if DECOUPLE_ME_RES
    context_ptr->me_fifo_ptr = eb_system_resource_get_producer_fifo(
            enc_handle_ptr->me_pool_ptr_array[0], 0);
#endif
    return EB_ErrorNone;
}

EbBool scene_transition_detector(
    PictureDecisionContext *context_ptr,
    SequenceControlSet                 *scs_ptr,
    PictureParentControlSet           **parent_pcs_window,
    uint32_t                                windowWidthFuture)
{
    PictureParentControlSet       *prev_pcs_ptr = parent_pcs_window[0];
    PictureParentControlSet       *current_pcs_ptr = parent_pcs_window[1];
    PictureParentControlSet       *future_pcs_ptr = parent_pcs_window[2];

    // calculating the frame threshold based on the number of 64x64 blocks in the frame
    uint32_t  region_threshhold;
    uint32_t  region_threshhold_chroma;
    // this variable determines whether the running average should be reset to equal the ahd or not after detecting a scene change.
    //EbBool reset_running_avg = context_ptr->reset_running_avg;

    EbBool is_abrupt_change; // this variable signals an abrubt change (scene change or flash)
    EbBool is_scene_change; // this variable signals a frame representing a scene change
    EbBool is_flash; // this variable signals a frame that contains a flash
    EbBool is_fade; // this variable signals a frame that contains a fade
    EbBool gradual_change; // this signals the detection of a light scene change a small/localized flash or the start of a fade

    uint32_t  ahd; // accumulative histogram (absolute) differences between the past and current frame

    uint32_t  ahd_cb;
    uint32_t  ahd_cr;

    uint32_t  ahd_error_cb = 0;
    uint32_t  ahd_error_cr = 0;

    uint32_t **ahd_running_avg_cb = context_ptr->ahd_running_avg_cb;
    uint32_t **ahd_running_avg_cr = context_ptr->ahd_running_avg_cr;
    uint32_t **ahd_running_avg = context_ptr->ahd_running_avg;

    uint32_t  ahd_error = 0; // the difference between the ahd and the running average at the current frame.

    uint8_t   aid_future_past = 0; // this variable denotes the average intensity difference between the next and the past frames
    uint8_t   aid_future_present = 0;
    uint8_t   aid_present_past = 0;

    uint32_t  bin = 0; // variable used to iterate through the bins of the histograms

    uint32_t  region_in_picture_width_index;
    uint32_t  region_in_picture_height_index;

    uint32_t  region_width;
    uint32_t  region_height;
    uint32_t  region_width_offset;
    uint32_t  region_height_offset;

    uint32_t  is_abrupt_change_count = 0;
    uint32_t  is_scene_change_count = 0;

    uint32_t  region_count_threshold = (scs_ptr->scd_mode == SCD_MODE_2) ?
        (uint32_t)(((float)((scs_ptr->picture_analysis_number_of_regions_per_width * scs_ptr->picture_analysis_number_of_regions_per_height) * 75) / 100) + 0.5) :
        (uint32_t)(((float)((scs_ptr->picture_analysis_number_of_regions_per_width * scs_ptr->picture_analysis_number_of_regions_per_height) * 50) / 100) + 0.5);

    region_width = parent_pcs_window[1]->enhanced_picture_ptr->width / scs_ptr->picture_analysis_number_of_regions_per_width;
    region_height = parent_pcs_window[1]->enhanced_picture_ptr->height / scs_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (region_in_picture_width_index = 0; region_in_picture_width_index < scs_ptr->picture_analysis_number_of_regions_per_width; region_in_picture_width_index++) {  // loop over horizontal regions
        for (region_in_picture_height_index = 0; region_in_picture_height_index < scs_ptr->picture_analysis_number_of_regions_per_height; region_in_picture_height_index++) { // loop over vertical regions

            is_abrupt_change = EB_FALSE;
            is_scene_change = EB_FALSE;
            is_flash = EB_FALSE;
            gradual_change = EB_FALSE;

            // Reset accumulative histogram (absolute) differences between the past and current frame
            ahd = 0;
            ahd_cb = 0;
            ahd_cr = 0;

            region_width_offset = (region_in_picture_width_index == scs_ptr->picture_analysis_number_of_regions_per_width - 1) ?
                parent_pcs_window[1]->enhanced_picture_ptr->width - (scs_ptr->picture_analysis_number_of_regions_per_width * region_width) :
                0;

            region_height_offset = (region_in_picture_height_index == scs_ptr->picture_analysis_number_of_regions_per_height - 1) ?
                parent_pcs_window[1]->enhanced_picture_ptr->height - (scs_ptr->picture_analysis_number_of_regions_per_height * region_height) :
                0;

            region_width += region_width_offset;
            region_height += region_height_offset;

            region_threshhold = (
                // Noise insertion/removal detection
                ((ABS((int64_t)current_pcs_ptr->pic_avg_variance - (int64_t)prev_pcs_ptr->pic_avg_variance)) > NOISE_VARIANCE_TH) &&
                (current_pcs_ptr->pic_avg_variance > HIGH_PICTURE_VARIANCE_TH || prev_pcs_ptr->pic_avg_variance > HIGH_PICTURE_VARIANCE_TH)) ?
                NOISY_SCENE_TH * NUM64x64INPIC(region_width, region_height) : // SCD TH function of noise insertion/removal.
                SCENE_TH * NUM64x64INPIC(region_width, region_height);

            region_threshhold_chroma = region_threshhold / 4;

            for (bin = 0; bin < HISTOGRAM_NUMBER_OF_BINS; ++bin) {
                ahd += ABS((int32_t)current_pcs_ptr->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][0][bin] - (int32_t)prev_pcs_ptr->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][0][bin]);
                ahd_cb += ABS((int32_t)current_pcs_ptr->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][1][bin] - (int32_t)prev_pcs_ptr->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][1][bin]);
                ahd_cr += ABS((int32_t)current_pcs_ptr->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][2][bin] - (int32_t)prev_pcs_ptr->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][2][bin]);
            }

            if (context_ptr->reset_running_avg) {
                ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] = ahd;
                ahd_running_avg_cb[region_in_picture_width_index][region_in_picture_height_index] = ahd_cb;
                ahd_running_avg_cr[region_in_picture_width_index][region_in_picture_height_index] = ahd_cr;
            }

            ahd_error = ABS((int32_t)ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] - (int32_t)ahd);
            ahd_error_cb = ABS((int32_t)ahd_running_avg_cb[region_in_picture_width_index][region_in_picture_height_index] - (int32_t)ahd_cb);
            ahd_error_cr = ABS((int32_t)ahd_running_avg_cr[region_in_picture_width_index][region_in_picture_height_index] - (int32_t)ahd_cr);

            if ((ahd_error > region_threshhold       && ahd >= ahd_error) ||
                (ahd_error_cb > region_threshhold_chroma && ahd_cb >= ahd_error_cb) ||
                (ahd_error_cr > region_threshhold_chroma && ahd_cr >= ahd_error_cr)) {
                is_abrupt_change = EB_TRUE;
            }
            else if ((ahd_error > (region_threshhold >> 1)) && ahd >= ahd_error)
                gradual_change = EB_TRUE;
            if (is_abrupt_change)
            {
                aid_future_past = (uint8_t)ABS((int16_t)future_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index][0] - (int16_t)prev_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index][0]);
                aid_future_present = (uint8_t)ABS((int16_t)future_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index][0] - (int16_t)current_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index][0]);
                aid_present_past = (uint8_t)ABS((int16_t)current_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index][0] - (int16_t)prev_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index][0]);

                if (aid_future_past < FLASH_TH && aid_future_present >= FLASH_TH && aid_present_past >= FLASH_TH) {
                    is_flash = EB_TRUE;
                    //SVT_LOG ("\nFlash in frame# %i , %i\n", current_pcs_ptr->picture_number,aid_future_past);
                }
                else if (aid_future_present < FADE_TH && aid_present_past < FADE_TH) {
                    is_fade = EB_TRUE;
                    //SVT_LOG ("\nFlash in frame# %i , %i\n", current_pcs_ptr->picture_number,aid_future_past);
                }
                else {
                    is_scene_change = EB_TRUE;
                    //SVT_LOG ("\nScene Change in frame# %i , %i\n", current_pcs_ptr->picture_number,aid_future_past);
                }
            }
            else if (gradual_change) {
                aid_future_past = (uint8_t)ABS((int16_t)future_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index][0] - (int16_t)prev_pcs_ptr->average_intensity_per_region[region_in_picture_width_index][region_in_picture_height_index][0]);
                if (aid_future_past < FLASH_TH) {
                    // proper action to be signalled
                    //SVT_LOG ("\nLight Flash in frame# %i , %i\n", current_pcs_ptr->picture_number,aid_future_past);
                    ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] = (3 * ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] + ahd) / 4;
                }
                else {
                    // proper action to be signalled
                    //SVT_LOG ("\nLight Scene Change / fade detected in frame# %i , %i\n", current_pcs_ptr->picture_number,aid_future_past);
                    ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] = (3 * ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] + ahd) / 4;
                }
            }
            else
                ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] = (3 * ahd_running_avg[region_in_picture_width_index][region_in_picture_height_index] + ahd) / 4;
            is_abrupt_change_count += is_abrupt_change;
            is_scene_change_count += is_scene_change;
        }
    }

    (void)windowWidthFuture;
    (void)is_flash;
    (void)is_fade;

    if (is_abrupt_change_count >= region_count_threshold)
        context_ptr->reset_running_avg = EB_TRUE;
    else
        context_ptr->reset_running_avg = EB_FALSE;
    if ((is_scene_change_count >= region_count_threshold) && ((!parent_pcs_window[1]->fade_in_to_black) && (!parent_pcs_window[1]->fade_out_from_black)))
        return(EB_TRUE);
    else
        return(EB_FALSE);
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
    EbErrorType return_error = EB_ErrorNone;

    uint32_t MinigopIndex;

    // Loop over all mini GOPs
    for (MinigopIndex = 0; MinigopIndex < MINI_GOP_MAX_COUNT; ++MinigopIndex) {
        context_ptr->mini_gop_activity_array[MinigopIndex] = (get_mini_gop_stats(MinigopIndex)->hierarchical_levels == 3 /*MIN_HIERARCHICAL_LEVEL*/) ?
            EB_FALSE :
            EB_TRUE;
    }

    return return_error;
}

/***************************************************************************************************
* Generates block picture map
*
*
***************************************************************************************************/
EbErrorType generate_picture_window_split(
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    uint32_t    MinigopIndex;

    context_ptr->total_number_of_mini_gops = 0;

    // Loop over all mini GOPs
    MinigopIndex = 0;
    while (MinigopIndex < MINI_GOP_MAX_COUNT) {
        // Only for a valid mini GOP
        if (get_mini_gop_stats(MinigopIndex)->end_index < encode_context_ptr->pre_assignment_buffer_count && context_ptr->mini_gop_activity_array[MinigopIndex] == EB_FALSE) {
            context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(MinigopIndex)->start_index;
            context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(MinigopIndex)->end_index;
            context_ptr->mini_gop_length[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(MinigopIndex)->lenght;
            context_ptr->mini_gop_hierarchical_levels[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(MinigopIndex)->hierarchical_levels;
            context_ptr->mini_gop_intra_count[context_ptr->total_number_of_mini_gops] = 0;
            context_ptr->mini_gop_idr_count[context_ptr->total_number_of_mini_gops] = 0;

            context_ptr->total_number_of_mini_gops++;
        }

        MinigopIndex += context_ptr->mini_gop_activity_array[MinigopIndex] ?
            1 :
            mini_gop_offset[get_mini_gop_stats(MinigopIndex)->hierarchical_levels - MIN_HIERARCHICAL_LEVEL];
    }

    // Only in presence of at least 1 valid mini GOP
    if (context_ptr->total_number_of_mini_gops != 0) {
        context_ptr->mini_gop_intra_count[context_ptr->total_number_of_mini_gops - 1] = encode_context_ptr->pre_assignment_buffer_intra_count;
        context_ptr->mini_gop_idr_count[context_ptr->total_number_of_mini_gops - 1] = encode_context_ptr->pre_assignment_buffer_idr_count;
    }

    return return_error;
}

/***************************************************************************************************
* Handles an incomplete picture window map
*
*
***************************************************************************************************/
EbErrorType handle_incomplete_picture_window_map(
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    if (context_ptr->total_number_of_mini_gops == 0) {
        context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops] = 0;
        context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - 1;
        context_ptr->mini_gop_length[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops];
        context_ptr->mini_gop_hierarchical_levels[context_ptr->total_number_of_mini_gops] = 3;// MIN_HIERARCHICAL_LEVEL; // AMIR to be updated after other predictions are supported

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
#if DECOUPLE_ME_RES
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
        pcs_ptr->pred_struct_ptr->pred_type == EB_PRED_RANDOM_ACCESS &&
        pcs_ptr->idr_flag == EB_FALSE &&
        pcs_ptr->cra_flag == EB_FALSE) {

        return 1;
    }
    else {
        return 0;
    }
}
#endif
/***************************************************************************************************
* If a switch happens, then update the RPS of the base layer frame separating the 2 different prediction structures
* Clean up the reference queue dependant counts of the base layer frame separating the 2 different prediction structures
*
***************************************************************************************************/
EbErrorType update_base_layer_reference_queue_dependent_count(
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr,
    SequenceControlSet            *scs_ptr,
    uint32_t                         MinigopIndex) {
    if (!context_ptr || !encode_context_ptr || !scs_ptr)
        return EB_ErrorBadParameter;

    EbErrorType return_error = EB_ErrorNone;

    PaReferenceQueueEntry         *input_entry_ptr;
    uint32_t                         input_queue_index;

    PredictionStructure           *next_pred_struct_ptr;
    PredictionStructureEntry      *next_base_layer_pred_position_ptr;

    uint32_t                         dependant_list_positive_entries;
    uint32_t                         dependant_list_removed_entries;
    uint32_t                         dep_list_count;

    uint32_t                         dep_idx;
    uint64_t                         dep_poc;

    PictureParentControlSet       *pcs_ptr;

    // Get the 1st PCS mini GOP
    pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_start_index[MinigopIndex]]->object_ptr;

#if DECOUPLE_ME_RES
    PictureParentControlSet  *trig_pcs; //triggering dep-cnt clean up picture, should be the first in dec order that goes to PicMgr
    if ( is_pic_cutting_short_ra_mg(context_ptr, pcs_ptr, MinigopIndex)) {
        trig_pcs = pcs_ptr;//this an LDP minigop so the first picture in minigop is the first in dec order
    }
    else {
        //this an RA minigop so the last picture in minigop is the first in dec order
        uint32_t  last_pic_in_mg = context_ptr->mini_gop_end_index[MinigopIndex];
        trig_pcs = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[last_pic_in_mg]->object_ptr;
    }
#endif
    // Derive the temporal layer difference between the current mini GOP and the previous mini GOP
    pcs_ptr->hierarchical_layers_diff = (uint8_t)(encode_context_ptr->previous_mini_gop_hierarchical_levels - pcs_ptr->hierarchical_levels);

    // Set init_pred_struct_position_flag to TRUE if mini GOP switch
    pcs_ptr->init_pred_struct_position_flag = encode_context_ptr->is_mini_gop_changed = (pcs_ptr->hierarchical_layers_diff != 0) ?
        EB_TRUE :
        EB_FALSE;

    // If the current mini GOP is different than the previous mini GOP update then update the positive dependant counts of the reference entry separating the 2 mini GOPs
    if (pcs_ptr->hierarchical_layers_diff != 0) {
        input_queue_index = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

        while (input_queue_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
            input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[input_queue_index];

#if DECOUPLE_ME_RES
            int32_t diff_n = 0;
#endif
            // Find the reference entry separating the 2 mini GOPs  (pcs_ptr->picture_number is the POC of the first isput in the mini GOP)
            if (input_entry_ptr->picture_number == (pcs_ptr->picture_number - 1)) {
                // Update the positive dependant counts

                // 1st step: remove all positive entries from the dependant list0 and dependant list1
                dependant_list_positive_entries = 0;
                for (dep_idx = 0; dep_idx < input_entry_ptr->list0.list_count; ++dep_idx) {
                    if (input_entry_ptr->list0.list[dep_idx] >= 0)
                        dependant_list_positive_entries++;
                }
                input_entry_ptr->list0.list_count = input_entry_ptr->list0.list_count - dependant_list_positive_entries;
                dependant_list_positive_entries = 0;
                for (dep_idx = 0; dep_idx < input_entry_ptr->list1.list_count; ++dep_idx) {
                    if (input_entry_ptr->list1.list[dep_idx] >= 0)
                        dependant_list_positive_entries++;
                }
                input_entry_ptr->list1.list_count = input_entry_ptr->list1.list_count - dependant_list_positive_entries;

                // 2nd step: inherit the positive dependant counts of the current mini GOP
                // Get the RPS set of the current mini GOP
                next_pred_struct_ptr = get_prediction_structure(
                    encode_context_ptr->prediction_structure_group_ptr,
                    pcs_ptr->pred_structure,
                    scs_ptr->reference_count,
                    pcs_ptr->hierarchical_levels);            // Number of temporal layer in the current mini GOP

                // Get the RPS of a base layer input
                next_base_layer_pred_position_ptr = next_pred_struct_ptr->pred_struct_entry_ptr_array[next_pred_struct_ptr->pred_struct_entry_count - 1];

                for (dep_idx = 0; dep_idx < next_base_layer_pred_position_ptr->dep_list0.list_count; ++dep_idx) {
                    if (next_base_layer_pred_position_ptr->dep_list0.list[dep_idx] >= 0)
                        input_entry_ptr->list0.list[input_entry_ptr->list0.list_count++] = next_base_layer_pred_position_ptr->dep_list0.list[dep_idx];
                }

                for (dep_idx = 0; dep_idx < next_base_layer_pred_position_ptr->dep_list1.list_count; ++dep_idx) {
                    if (next_base_layer_pred_position_ptr->dep_list1.list[dep_idx] >= 0)
                        input_entry_ptr->list1.list[input_entry_ptr->list1.list_count++] = next_base_layer_pred_position_ptr->dep_list1.list[dep_idx];
                }

#if DECOUPLE_ME_RES
                diff_n =
                    (input_entry_ptr->list0.list_count + input_entry_ptr->list1.list_count) - //depCnt after clean up
                    (input_entry_ptr->dep_list0_count + input_entry_ptr->dep_list1_count); //depCnt from org prediction struct


                //these refs are defintely not in the pre-assignment buffer
                if (diff_n) {
                    //TODO: change triggereing picture from the firt pic to the last pic in the MG for RA ( first in dec-order = base layer)
                    trig_pcs->updated_links_arr[trig_pcs->other_updated_links_cnt].pic_num = input_entry_ptr->picture_number;
                    trig_pcs->updated_links_arr[trig_pcs->other_updated_links_cnt++].dep_cnt_diff = diff_n;
                }
#endif
                // 3rd step: update the dependant count
                dependant_list_removed_entries = input_entry_ptr->dep_list0_count + input_entry_ptr->dep_list1_count - input_entry_ptr->dependent_count;
                input_entry_ptr->dep_list0_count = (input_entry_ptr->is_alt_ref) ? input_entry_ptr->list0.list_count + 1 : input_entry_ptr->list0.list_count;
                input_entry_ptr->dep_list1_count = input_entry_ptr->list1.list_count;
                input_entry_ptr->dependent_count = input_entry_ptr->dep_list0_count + input_entry_ptr->dep_list1_count - dependant_list_removed_entries;
            }
            else {
                // Modify Dependent List0
                dep_list_count = input_entry_ptr->list0.list_count;
                for (dep_idx = 0; dep_idx < dep_list_count; ++dep_idx) {
                    // Adjust the latest currentInputPoc in case we're in a POC rollover scenario
                    // currentInputPoc += (currentInputPoc < input_entry_ptr->pocNumber) ? (1 << scs_ptr->bitsForPictureOrderCount) : 0;
                    dep_poc = POC_CIRCULAR_ADD(
                        input_entry_ptr->picture_number, // can't use a value that gets reset
                        input_entry_ptr->list0.list[dep_idx]/*,
                                                         scs_ptr->bitsForPictureOrderCount*/);

                                                         // If Dependent POC is greater or equal to the IDR POC
                    if (dep_poc >= pcs_ptr->picture_number && input_entry_ptr->list0.list[dep_idx]) {
                        input_entry_ptr->list0.list[dep_idx] = 0;

                        // Decrement the Reference's reference_count
                        --input_entry_ptr->dependent_count;
#if DECOUPLE_ME_RES
                        diff_n--;
#endif
                        CHECK_REPORT_ERROR(
                            (input_entry_ptr->dependent_count != ~0u),
                            encode_context_ptr->app_callback_ptr,
                            EB_ENC_PD_ERROR3);
                    }
                }
                // Modify Dependent List1
                dep_list_count = input_entry_ptr->list1.list_count;
                for (dep_idx = 0; dep_idx < dep_list_count; ++dep_idx) {
                    // Adjust the latest currentInputPoc in case we're in a POC rollover scenario
                    // currentInputPoc += (currentInputPoc < input_entry_ptr->pocNumber) ? (1 << scs_ptr->bitsForPictureOrderCount) : 0;
                    dep_poc = POC_CIRCULAR_ADD(
                        input_entry_ptr->picture_number,
                        input_entry_ptr->list1.list[dep_idx]/*,
                                                         scs_ptr->bitsForPictureOrderCount*/);

                    // If Dependent POC is greater or equal to the IDR POC
                    if ((dep_poc >= pcs_ptr->picture_number) && input_entry_ptr->list1.list[dep_idx]) {
                        input_entry_ptr->list1.list[dep_idx] = 0;
                        // Decrement the Reference's reference_count
                        --input_entry_ptr->dependent_count;
#if DECOUPLE_ME_RES
                        diff_n--;
#endif
                        CHECK_REPORT_ERROR(
                            (input_entry_ptr->dependent_count != ~0u),
                            encode_context_ptr->app_callback_ptr,
                            EB_ENC_PD_ERROR3);
                    }
                }

#if DECOUPLE_ME_RES
                //these refs are defintely not in the pre-ass buffer
                if (diff_n) {
                    trig_pcs->updated_links_arr[trig_pcs->other_updated_links_cnt].pic_num = input_entry_ptr->picture_number;
                    trig_pcs->updated_links_arr[trig_pcs->other_updated_links_cnt++].dep_cnt_diff = diff_n;
                }
#endif
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
            scs_ptr = (SequenceControlSet*)pcs_ptr->scs_wrapper_ptr->object_ptr;
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

#if TF_LEVELS
#if IMPROVED_TF_LEVELS
void set_tf_controls(PictureDecisionContext *context_ptr, uint8_t tf_level) {

    TfControls *tf_ctrls = &context_ptr->tf_ctrls;

    switch (tf_level)
    {
    case 0:
        tf_ctrls->enabled = 0;
        break;
    case 1:
        tf_ctrls->enabled = 1;
        tf_ctrls->window_size = 7;
        tf_ctrls->noise_based_window_adjust = 1;
        break;
    case 2:
        tf_ctrls->enabled = 1;
        tf_ctrls->window_size = 3;
        tf_ctrls->noise_based_window_adjust = 1;
        break;
    case 3:
        tf_ctrls->enabled = 1;
        tf_ctrls->window_size = 3;
        tf_ctrls->noise_based_window_adjust = 0;
        break;
    default:
        assert(0);
        break;
    }
}
#else
void set_tf_controls(PictureDecisionContext *context_ptr, uint8_t tf_level) {

    TfControls *tf_ctrls = &context_ptr->tf_ctrls;

    switch (tf_level)
    {
    case 0:
        tf_ctrls->enabled = 1;
        tf_ctrls->window_size = 7;
        break;
    case 1:
        tf_ctrls->enabled = 1;
        tf_ctrls->window_size = 5;
        break;
    case 2:
        tf_ctrls->enabled = 1;
        tf_ctrls->window_size = 3;
        break;
    case 3:
        tf_ctrls->enabled = 0;
        tf_ctrls->window_size = 3;
        break;
    default:
        assert(0);
        break;
    }
}
#endif
#endif
/******************************************************
* Derive Multi-Processes Settings for OQ
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/
#if TF_LEVELS
EbErrorType signal_derivation_multi_processes_oq(
    SequenceControlSet *scs_ptr,
    PictureParentControlSet *pcs_ptr,
    PictureDecisionContext *context_ptr) {
#else
EbErrorType signal_derivation_multi_processes_oq(
    SequenceControlSet        *scs_ptr,
    PictureParentControlSet   *pcs_ptr) {
#endif
    EbErrorType return_error = EB_ErrorNone;
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;

#if !UNIFY_SC_NSC
    uint8_t sc_content_detected = pcs_ptr->sc_content_detected;
#endif
#if !REFACTOR_ME_HME
    uint8_t enc_mode_hme = scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#endif
#if REFACTOR_ME_HME
    // If enabled here, the hme enable flags should also be enabled in ResourceCoordinationProcess
    // to ensure that resources are allocated for the downsampled pictures used in HME
    pcs_ptr->enable_hme_flag        = 1;
    pcs_ptr->enable_hme_level0_flag = 1;

#if ALLOW_HME_L1L2_REFINEMENT
    pcs_ptr->enable_hme_level1_flag = 1;
    pcs_ptr->enable_hme_level2_flag = 1;
#else
    if (sc_content_detected)
#if APR23_ADOPTIONS_2
        // HME refinement was turned on everywhere before the branch was used for M8 testing;
        // this restores the previous settings
        if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if PRESETS_SHIFT
        if (pcs_ptr->enc_mode <= ENC_M2) {
#else
        if (pcs_ptr->enc_mode <= ENC_M3) {
#endif
#endif
            pcs_ptr->enable_hme_level1_flag = 1;
            pcs_ptr->enable_hme_level2_flag = 1;
        }
        else {
#if 0//APR23_ADOPTIONS
            pcs_ptr->enable_hme_level1_flag = 1;
            pcs_ptr->enable_hme_level2_flag = 1;
#else
            pcs_ptr->enable_hme_level1_flag = 0;
            pcs_ptr->enable_hme_level2_flag = 0;
#endif
        }
    else {
#if UPGRADE_M6_M7_M8
        if (pcs_ptr->enc_mode <= ENC_M7) {
            pcs_ptr->enable_hme_level1_flag = 1;
            pcs_ptr->enable_hme_level2_flag = 1;
        }
        else {
            pcs_ptr->enable_hme_level1_flag = 1;
            pcs_ptr->enable_hme_level2_flag = 0;
        }
#else
        pcs_ptr->enable_hme_level1_flag = 1;
        pcs_ptr->enable_hme_level2_flag = 1;
#endif
    }
#endif

    pcs_ptr->tf_enable_hme_flag = 1;
    pcs_ptr->tf_enable_hme_level0_flag = 1;
#if PRESETS_SHIFT
    // Can enable everywhere b/c TF is off for SC anyway; remove fake diff
#if UPGRADE_M6_M7_M8
#if JUNE26_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if JUNE25_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M4) {
#else
#if PRESET_SHIFITNG
    if (pcs_ptr->enc_mode <= ENC_M5) {
#else
    if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
#endif
#endif
        pcs_ptr->tf_enable_hme_level1_flag = 1;
        pcs_ptr->tf_enable_hme_level2_flag = 1;
    }
    else {
        pcs_ptr->tf_enable_hme_level1_flag = 1;
        pcs_ptr->tf_enable_hme_level2_flag = 0;
    }
#else
    pcs_ptr->tf_enable_hme_level1_flag = 1;
    pcs_ptr->tf_enable_hme_level2_flag = 1;
#endif
#else
    if (sc_content_detected)
        if (pcs_ptr->enc_mode <= ENC_M3) {
            pcs_ptr->tf_enable_hme_level1_flag = 1;
            pcs_ptr->tf_enable_hme_level2_flag = 1;
        }
        else {
            pcs_ptr->tf_enable_hme_level1_flag = 0;
            pcs_ptr->tf_enable_hme_level2_flag = 0;
        }
    else {
        pcs_ptr->tf_enable_hme_level1_flag = 1;
        pcs_ptr->tf_enable_hme_level2_flag = 1;
    }
#endif
#else
    pcs_ptr->enable_hme_flag = enable_hme_flag[pcs_ptr->sc_content_detected][scs_ptr->input_resolution][enc_mode_hme];

    pcs_ptr->enable_hme_level0_flag = enable_hme_level0_flag[pcs_ptr->sc_content_detected][scs_ptr->input_resolution][enc_mode_hme];
    pcs_ptr->enable_hme_level1_flag = enable_hme_level1_flag[pcs_ptr->sc_content_detected][scs_ptr->input_resolution][enc_mode_hme];
    pcs_ptr->enable_hme_level2_flag = enable_hme_level2_flag[pcs_ptr->sc_content_detected][scs_ptr->input_resolution][enc_mode_hme];

    pcs_ptr->tf_enable_hme_flag = tf_enable_hme_flag[pcs_ptr->sc_content_detected][scs_ptr->input_resolution][enc_mode_hme];
    pcs_ptr->tf_enable_hme_level0_flag = tf_enable_hme_level0_flag[pcs_ptr->sc_content_detected][scs_ptr->input_resolution][enc_mode_hme];
    pcs_ptr->tf_enable_hme_level1_flag = tf_enable_hme_level1_flag[pcs_ptr->sc_content_detected][scs_ptr->input_resolution][enc_mode_hme];
    pcs_ptr->tf_enable_hme_level2_flag = tf_enable_hme_level2_flag[pcs_ptr->sc_content_detected][scs_ptr->input_resolution][enc_mode_hme];

#endif
#if DEPTH_PART_CLEAN_UP
    // Set the Multi-Pass PD level
#if ADD_NEW_MPPD_LEVEL
#if ADOPT_SKIPPING_PD1
    pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_LEVEL_0;
#else
#if MAR23_ADOPTIONS
    if (sc_content_detected)
#if MAR30_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M3)
            pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_LEVEL_1;
        else
            pcs_ptr->multi_pass_pd_level =
            (pcs_ptr->slice_type == I_SLICE)
            ? MULTI_PASS_PD_LEVEL_1
            : MULTI_PASS_PD_LEVEL_2;
    else if (pcs_ptr->enc_mode <= ENC_M0)
#else
        pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_LEVEL_1;
    else if (pcs_ptr->enc_mode <= ENC_M1)
#endif
        pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_LEVEL_1;
    else
        pcs_ptr->multi_pass_pd_level =
        (pcs_ptr->slice_type == I_SLICE)
        ? MULTI_PASS_PD_LEVEL_1
        : MULTI_PASS_PD_LEVEL_2;
#else
    if (sc_content_detected)
#if MAR19_ADOPTIONS
        // Use a single-stage PD if I_SLICE
        pcs_ptr->multi_pass_pd_level =
        (pcs_ptr->slice_type == I_SLICE)
        ? MULTI_PASS_PD_OFF
        : MULTI_PASS_PD_LEVEL_1;
#else
        pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_OFF;
#endif
    else if (pcs_ptr->enc_mode <= ENC_M1)
        // Use a single-stage PD if I_SLICE
        pcs_ptr->multi_pass_pd_level =
        (pcs_ptr->slice_type == I_SLICE)
        ? MULTI_PASS_PD_OFF
        : pcs_ptr->hierarchical_levels < 3
        ? MULTI_PASS_PD_LEVEL_0
        : MULTI_PASS_PD_LEVEL_1;
    else if (pcs_ptr->enc_mode <= ENC_M7)
        // Use a single-stage PD if I_SLICE
        pcs_ptr->multi_pass_pd_level =
        (pcs_ptr->slice_type == I_SLICE)
        ? MULTI_PASS_PD_OFF
        : pcs_ptr->hierarchical_levels < 3
        ? MULTI_PASS_PD_LEVEL_0
        : MULTI_PASS_PD_LEVEL_2;
    else
        // Use a single-stage PD if I_SLICE
        pcs_ptr->multi_pass_pd_level =
        (pcs_ptr->slice_type == I_SLICE)
        ? MULTI_PASS_PD_OFF
        : MULTI_PASS_PD_LEVEL_3;
#endif
#endif
#else
    if (sc_content_detected)
#if MAR17_ADOPTIONS
        pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_OFF;
#else
        if (pcs_ptr->enc_mode <= ENC_M4)
            pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_OFF;
        else
            // Use a single-stage PD if I_SLICE
            pcs_ptr->multi_pass_pd_level =
            (pcs_ptr->slice_type == I_SLICE)
            ? MULTI_PASS_PD_OFF
            : MULTI_PASS_PD_LEVEL_2;
#endif
    else if (MR_MODE)
        pcs_ptr->multi_pass_pd_level =
        (pcs_ptr->slice_type == I_SLICE)
        ? MULTI_PASS_PD_OFF
        : MULTI_PASS_PD_LEVEL_0;
#if MAR17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7)
#else
    else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
        // Use a single-stage PD if I_SLICE
        pcs_ptr->multi_pass_pd_level =
        (pcs_ptr->slice_type == I_SLICE)
        ? MULTI_PASS_PD_OFF
        : pcs_ptr->hierarchical_levels < 3
        ? MULTI_PASS_PD_LEVEL_0
        : MULTI_PASS_PD_LEVEL_1;
    else
        // Use a single-stage PD if I_SLICE
        pcs_ptr->multi_pass_pd_level =
        (pcs_ptr->slice_type == I_SLICE)
        ? MULTI_PASS_PD_OFF
        : MULTI_PASS_PD_LEVEL_2;
#if !DEPTH_PART_CLEAN_UP
    // If ADP then set multi_pass_pd_level to INVALID
    if(pcs_ptr->adp_level != ADP_OFF)
        pcs_ptr->multi_pass_pd_level = MULTI_PASS_PD_INVALID;
#endif
#endif


    // Set disallow_nsq
#if APR25_12AM_ADOPTIONS
#if JUNE11_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if APR25_1PM_ADOPTIONS
#if PRESET_SHIFITNG
#if PRESET_SHIFITNG
    if (pcs_ptr->enc_mode <= ENC_M3 || (pcs_ptr->enc_mode <= ENC_M4 && sc_content_detected)) {
#else
    if (pcs_ptr->enc_mode <= ENC_M3 || (pcs_ptr->enc_mode <= ENC_M6 && sc_content_detected)) {
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M5 || (pcs_ptr->enc_mode <= ENC_M6 && sc_content_detected)) {
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
#endif
        pcs_ptr->disallow_nsq = EB_FALSE;
    }
#if !JUNE11_ADOPTIONS
#if APR25_3AM_ADOPTIONS
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M4 && !sc_content_detected) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M6 && !sc_content_detected) {
#endif
#else
    else if (pcs_ptr->enc_mode <= ENC_M6) {
#endif
        pcs_ptr->disallow_nsq = pcs_ptr->is_used_as_reference_flag ? EB_FALSE : EB_TRUE;
    }
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
#if APR25_7PM_ADOPTIONS
        if (sc_content_detected)
            pcs_ptr->disallow_nsq = EB_FALSE;
        else
            pcs_ptr->disallow_nsq = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
#else
        pcs_ptr->disallow_nsq = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
#endif
    }
#endif
#if JUNE17_ADOPTIONS && !SOFT_CYCLES_M6M7
    else if (pcs_ptr->enc_mode <= ENC_M6) {
#if !UNIFY_SC_NSC
        if (sc_content_detected)
            pcs_ptr->disallow_nsq = EB_FALSE;
        else
#endif
            pcs_ptr->disallow_nsq = pcs_ptr->temporal_layer_index == 0 ? EB_FALSE : EB_TRUE;
    }
#endif
    else {
#if APR25_7PM_ADOPTIONS
#if !UNIFY_SC_NSC
        if (sc_content_detected)
            pcs_ptr->disallow_nsq = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
        else
#endif
#if JUNE26_ADOPTIONS
            pcs_ptr->disallow_nsq = EB_TRUE;
#else
#if SOFT_CYCLES_M6M7
            pcs_ptr->disallow_nsq = EB_FALSE;
#else
            pcs_ptr->disallow_nsq = EB_TRUE;
#endif
#endif
#else
        pcs_ptr->disallow_nsq = EB_TRUE;
#endif
    }
#else
#if M8_NSQ
#if M5_I_NSQ
#if ALLOW_NSQ_M6
    if (pcs_ptr->enc_mode <= ENC_M6) {
#else
    if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
        pcs_ptr->disallow_nsq = EB_FALSE;
    }
#if UPGRADE_M6_M7_M8
#if ALLOW_NSQ_M6
    else if (pcs_ptr->enc_mode <= ENC_M7) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M6) {
#endif
        if (pcs_ptr->sc_content_detected)
            pcs_ptr->disallow_nsq = pcs_ptr->is_used_as_reference_flag ? EB_FALSE : EB_TRUE;
        else
            pcs_ptr->disallow_nsq = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
    }
#endif
    else {
        pcs_ptr->disallow_nsq = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
    }
#else
    pcs_ptr->disallow_nsq = pcs_ptr->enc_mode <= ENC_M5 ? EB_FALSE : EB_TRUE;
#endif
#else
    pcs_ptr->disallow_nsq = EB_FALSE;
#endif
#endif
#if NSQ_REMOVAL_CODE_CLEAN_UP
    pcs_ptr->max_number_of_pus_per_sb = SQUARE_PU_COUNT;
#else
    if (!pcs_ptr->disallow_nsq)
        assert(scs_ptr->nsq_present == 1 && "use nsq_present 1");

    pcs_ptr->max_number_of_pus_per_sb =
        pcs_ptr->disallow_nsq
        ? SQUARE_PU_COUNT
        : MAX_ME_PU_COUNT;
#endif
#if !REMOVE_UNUSED_CODE_PH2
    // Set sb_64x64_simulated - only allow when SB size is not already 64x64
#if MAR23_ADOPTIONS
    if (scs_ptr->static_config.super_block_size != 64) {
        if (pcs_ptr->sc_content_detected)
#if APR08_ADOPTIONS
            pcs_ptr->sb_64x64_simulated = EB_FALSE;
#else
            pcs_ptr->sb_64x64_simulated = EB_TRUE;
#endif
        else
            pcs_ptr->sb_64x64_simulated = EB_FALSE;
    }
    else
        pcs_ptr->sb_64x64_simulated = EB_FALSE;
#else
    if (sc_content_detected) {
        pcs_ptr->sb_64x64_simulated = EB_TRUE;
    }
    else {
        pcs_ptr->sb_64x64_simulated = EB_FALSE;
    }
#endif
#endif
#if !M8_4x4
    // Set disallow_4x4
    pcs_ptr->disallow_4x4 = EB_FALSE;
#endif
    // Set disallow_all_nsq_blocks_below_8x8: 8x4, 4x8
#if UNIFY_SC_NSC
#if NEW_M8
    if (pcs_ptr->enc_mode <= ENC_M8)
#else
#if SOFT_CYCLES_M6M7
    if (pcs_ptr->enc_mode <= ENC_M7)
#else
    if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#endif
        pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_FALSE;
    else
        pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_TRUE;
#else
    if (sc_content_detected) {
        pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_FALSE;
    }
    else {
#if JUNE11_ADOPTIONS
#if SOFT_CYCLES_REDUCTION
#if SOFT_CYCLES_M6M7
    if (pcs_ptr->enc_mode <= ENC_M7)
#else
        if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
        if (pcs_ptr->enc_mode <= ENC_M3)
#endif
            pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_FALSE;
#else
#if MAY16_7PM_ADOPTIONS
#if PRESET_SHIFITNG
        if (pcs_ptr->enc_mode <= ENC_M0 ||
            (pcs_ptr->enc_mode <= ENC_M2 && pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE))
            pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_FALSE;
        else if (pcs_ptr->enc_mode <= ENC_M2)
            pcs_ptr->disallow_all_nsq_blocks_below_8x8 = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
#else
        if (pcs_ptr->enc_mode <= ENC_M0 ||
            (pcs_ptr->enc_mode <= ENC_M4 && pcs_ptr->input_resolution <= INPUT_SIZE_1080p_RANGE))
            pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_FALSE;
        else if (pcs_ptr->enc_mode <= ENC_M4)
            pcs_ptr->disallow_all_nsq_blocks_below_8x8 = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
#endif
#else
#if PRESETS_SHIFT
        if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if MAR19_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M7)
#else
        if (pcs_ptr->enc_mode <= ENC_M4)
#endif
#endif
            pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_FALSE;
#endif
#endif
        else
            pcs_ptr->disallow_all_nsq_blocks_below_8x8 = EB_TRUE;
    }
#endif

    // Set disallow_all_nsq_blocks_below_16x16: 16x8, 8x16, 16x4, 4x16
    pcs_ptr->disallow_all_nsq_blocks_below_16x16 = EB_FALSE;

#if NO_NSQ_B32
    pcs_ptr->disallow_all_nsq_blocks_below_64x64 = EB_FALSE;
    pcs_ptr->disallow_all_nsq_blocks_below_32x32 = EB_FALSE;
#endif
#if NO_NSQ_ABOVE
    pcs_ptr->disallow_all_nsq_blocks_above_64x64= EB_FALSE;
#if APR25_7PM_ADOPTIONS
    // disallow_all_nsq_blocks_above_32x32
#if UNIFY_SC_NSC
    pcs_ptr->disallow_all_nsq_blocks_above_32x32 = EB_FALSE;
#else
#if PRESET_SHIFITNG
#if SOFT_CYCLES_REDUCTION
#if SOFT_CYCLES_M6M7
    if (!sc_content_detected || pcs_ptr->enc_mode <= ENC_M7)
#else
    if (!sc_content_detected || pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
    if (!sc_content_detected || pcs_ptr->enc_mode <= ENC_M3)
#endif
#else
    if (!sc_content_detected || pcs_ptr->enc_mode <= ENC_M5)
#endif
        pcs_ptr->disallow_all_nsq_blocks_above_32x32 = EB_FALSE;
    else
        pcs_ptr->disallow_all_nsq_blocks_above_32x32 = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
#endif
    // disallow_all_nsq_blocks_above_16x16
#if JUNE17_ADOPTIONS
#if SOFT_CYCLES_REDUCTION
#if NEW_M8
    if (pcs_ptr->enc_mode <= ENC_M8)
#else
#if SOFT_CYCLES_M6M7
    if (pcs_ptr->enc_mode <= ENC_M7)
#else
    if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M4)
#endif
#else
#if PRESET_SHIFITNG
    if (!sc_content_detected || pcs_ptr->enc_mode <= ENC_M4)
#else
    if (!sc_content_detected || pcs_ptr->enc_mode <= ENC_M6)
#endif
#endif
        pcs_ptr->disallow_all_nsq_blocks_above_16x16 = EB_FALSE;
#if JUNE17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M6)
#else
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M5)
#else
    else if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#endif
        pcs_ptr->disallow_all_nsq_blocks_above_16x16 = pcs_ptr->slice_type == I_SLICE ? EB_FALSE : EB_TRUE;
    else
        pcs_ptr->disallow_all_nsq_blocks_above_16x16 = EB_TRUE;
#else
    pcs_ptr->disallow_all_nsq_blocks_above_32x32= EB_FALSE;
    pcs_ptr->disallow_all_nsq_blocks_above_16x16= EB_FALSE;
#endif
#endif
#if NO_AB_HV4
    pcs_ptr->disallow_HVA_HVB_HV4 = EB_FALSE;
    pcs_ptr->disallow_HV4 = EB_FALSE;
#endif

    // Set disallow_all_non_hv_nsq_blocks_below_16x16
    pcs_ptr->disallow_all_non_hv_nsq_blocks_below_16x16 = EB_FALSE;

    // Set disallow_all_h4_v4_blocks_below_16x16
    pcs_ptr->disallow_all_h4_v4_blocks_below_16x16 = EB_FALSE;
#else

        if (sc_content_detected)
#if MAR10_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M4)
#else
            if (pcs_ptr->enc_mode <= ENC_M2)
#endif
                pcs_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
#if DEPTH_PART_CLEAN_UP
        else
            // Use a single-stage PD if I_SLICE
            pcs_ptr->pic_depth_mode =
            (pcs_ptr->slice_type == I_SLICE)
            ? PIC_ALL_DEPTH_MODE
            : PIC_MULTI_PASS_PD_MODE_2;
#else
#if MAR10_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M8)
#else
            else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
                if (pcs_ptr->slice_type == I_SLICE)
                    pcs_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
                else
                    pcs_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
            else
                pcs_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
#endif
        else if (MR_MODE)
            pcs_ptr->pic_depth_mode =
            (pcs_ptr->slice_type == I_SLICE)
            ? PIC_ALL_DEPTH_MODE
            : PIC_MULTI_PASS_PD_MODE_0;
#if MAR4_M6_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M5)
#else
        else if (pcs_ptr->enc_mode <= ENC_M1)
#endif
            // Use a single-stage PD if I_SLICE
            pcs_ptr->pic_depth_mode =
            (pcs_ptr->slice_type == I_SLICE)
            ? PIC_ALL_DEPTH_MODE
            : pcs_ptr->hierarchical_levels < 3
            ? PIC_MULTI_PASS_PD_MODE_0
            : PIC_MULTI_PASS_PD_MODE_1;
#if DEPTH_PART_CLEAN_UP
        else
            // Use a single-stage PD if I_SLICE
            pcs_ptr->pic_depth_mode =
            (pcs_ptr->slice_type == I_SLICE)
            ? PIC_ALL_DEPTH_MODE
            : PIC_MULTI_PASS_PD_MODE_2;
#else
#if MAR2_M7_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M5)
#else
        else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
            // Use a single-stage PD if I_SLICE
            pcs_ptr->pic_depth_mode = (pcs_ptr->slice_type == I_SLICE) ? PIC_ALL_DEPTH_MODE : PIC_MULTI_PASS_PD_MODE_2;
        else
            pcs_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;

        if (pcs_ptr->pic_depth_mode < PIC_SQ_DEPTH_MODE)
            assert(scs_ptr->nsq_present == 1 && "use nsq_present 1");

#endif
        pcs_ptr->max_number_of_pus_per_sb = (pcs_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) ? MAX_ME_PU_COUNT : SQUARE_PU_COUNT;
#endif
#if !DEPTH_PART_CLEAN_UP
    // NSQ search Level                               Settings
    // NSQ_SEARCH_OFF                                 OFF
    // NSQ_SEARCH_LEVEL1                              Allow only NSQ Inter-NEAREST/NEAR/GLOBAL if parent SQ has no coeff + reordering nsq_table number and testing only 1 NSQ SHAPE
    // NSQ_SEARCH_LEVEL2                              Allow only NSQ Inter-NEAREST/NEAR/GLOBAL if parent SQ has no coeff + reordering nsq_table number and testing only 2 NSQ SHAPE
    // NSQ_SEARCH_LEVEL3                              Allow only NSQ Inter-NEAREST/NEAR/GLOBAL if parent SQ has no coeff + reordering nsq_table number and testing only 3 NSQ SHAPE
    // NSQ_SEARCH_LEVEL4                              Allow only NSQ Inter-NEAREST/NEAR/GLOBAL if parent SQ has no coeff + reordering nsq_table number and testing only 4 NSQ SHAPE
    // NSQ_SEARCH_LEVEL5                              Allow only NSQ Inter-NEAREST/NEAR/GLOBAL if parent SQ has no coeff + reordering nsq_table number and testing only 5 NSQ SHAPE
    // NSQ_SEARCH_LEVEL6                              Allow only NSQ Inter-NEAREST/NEAR/GLOBAL if parent SQ has no coeff + reordering nsq_table number and testing only 6 NSQ SHAPE
    // NSQ_SEARCH_FULL                                Allow NSQ Intra-FULL and Inter-FULL

        if (MR_MODE)
            pcs_ptr->nsq_search_level = NSQ_SEARCH_FULL;
        else if (pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_0 ||
                 pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_1 ||
                 pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2 ||
                 pcs_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3 ){

            pcs_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;
        }
        else if (pcs_ptr->enc_mode <= ENC_M1)
            pcs_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;
        else if (pcs_ptr->enc_mode <= ENC_M1)
            pcs_ptr->nsq_search_level = (pcs_ptr->is_used_as_reference_flag) ? NSQ_SEARCH_LEVEL6 : NSQ_SEARCH_LEVEL3;
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M8)
#else
        else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
            if (pcs_ptr->is_used_as_reference_flag)
                pcs_ptr->nsq_search_level = NSQ_SEARCH_LEVEL3;
            else
                pcs_ptr->nsq_search_level = NSQ_SEARCH_LEVEL1;
        else
            pcs_ptr->nsq_search_level = NSQ_SEARCH_OFF;
    if (pcs_ptr->nsq_search_level > NSQ_SEARCH_OFF)
        assert(scs_ptr->nsq_present == 1 && "use nsq_present 1");

    switch (pcs_ptr->nsq_search_level) {
    case NSQ_SEARCH_OFF:
        pcs_ptr->nsq_max_shapes_md = 0;
        break;
    case NSQ_SEARCH_LEVEL1:
        pcs_ptr->nsq_max_shapes_md = 1;
        break;
    case NSQ_SEARCH_LEVEL2:
        pcs_ptr->nsq_max_shapes_md = 2;
        break;
    case NSQ_SEARCH_LEVEL3:
        pcs_ptr->nsq_max_shapes_md = 3;
        break;
    case NSQ_SEARCH_LEVEL4:
        pcs_ptr->nsq_max_shapes_md = 4;
        break;
    case NSQ_SEARCH_LEVEL5:
        pcs_ptr->nsq_max_shapes_md = 5;
        break;
    case NSQ_SEARCH_LEVEL6:
        pcs_ptr->nsq_max_shapes_md = 6;
        break;
    case NSQ_SEARCH_FULL:
        pcs_ptr->nsq_max_shapes_md = 6;
        break;
    default:
        SVT_LOG("nsq_search_level is not supported\n");
        break;
    }

    if (pcs_ptr->nsq_search_level == NSQ_SEARCH_OFF)
        if (pcs_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) pcs_ptr->pic_depth_mode = PIC_SQ_DEPTH_MODE;
    if (pcs_ptr->pic_depth_mode > PIC_SQ_DEPTH_MODE)
        assert(pcs_ptr->nsq_search_level == NSQ_SEARCH_OFF);
#endif
    // Loop filter Level                            Settings
    // 0                                            OFF
    // 1                                            CU-BASED
    // 2                                            LIGHT FRAME-BASED
    // 3                                            FULL FRAME-BASED

    //for now only I frames are allowed to use sc tools.
    //TODO: we can force all frames in GOP with the same detection status of leading I frame.
    if (pcs_ptr->slice_type == I_SLICE) {
        frm_hdr->allow_screen_content_tools = pcs_ptr->sc_content_detected;
        if (scs_ptr->static_config.intrabc_mode == DEFAULT) {
#if IBC_SYNCH
            // ENABLE/DISABLE IBC
#if MAR3_M6_ADOPTIONS
#if MAR4_M8_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M8)
#else
            if (pcs_ptr->enc_mode <= ENC_M6)
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M5)
#endif
            {
                frm_hdr->allow_intrabc =  pcs_ptr->sc_content_detected;
            } else {
                frm_hdr->allow_intrabc =  0;
            }

            //IBC Modes:   0: OFF 1:Slow   2:Faster   3:Fastest
#if MAR4_M8_ADOPTIONS
#if M8_IBC
#if UPGRADE_M6_M7_M8
#if M5_I_IBC
            if (pcs_ptr->enc_mode <= ENC_M8)
#else
            if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M8)
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
            {
                pcs_ptr->ibc_mode = 1; // Slow
            } else {
                pcs_ptr->ibc_mode = 2; // Faster
            }

#else
            //IBC Modes:   0: OFF 1:Slow   2:Faster   3:Fastest
            if (pcs_ptr->enc_mode <= ENC_M5) {
                frm_hdr->allow_intrabc =  pcs_ptr->sc_content_detected;
                pcs_ptr->ibc_mode = 1; // Slow
            } else {
                frm_hdr->allow_intrabc =  0;
                pcs_ptr->ibc_mode = 2; // Faster
            }
#endif
        } else {
            frm_hdr->allow_intrabc =  (uint8_t)(scs_ptr->static_config.intrabc_mode > 0);
            pcs_ptr->ibc_mode = (uint8_t)scs_ptr->static_config.intrabc_mode;
        }
    }
    else {
        //this will enable sc tools for P frames. hence change Bitstream even if palette mode is OFF
        frm_hdr->allow_screen_content_tools = pcs_ptr->sc_content_detected;
        frm_hdr->allow_intrabc = 0;
        pcs_ptr->ibc_mode = 0; // OFF
    }

   /*Palette Modes:
        0:OFF
        1:Slow    NIC=7/4/4
        2:        NIC=7/2/2
        3:        NIC=7/2/2 + No K means for non ref
        4:        NIC=4/2/1
        5:        NIC=4/2/1 + No K means for Inter frame
        6:Fastest NIC=4/2/1 + No K means for non base + step for non base for most dominent

    */
    if (frm_hdr->allow_screen_content_tools)
        if (scs_ptr->static_config.enable_palette == -1)//auto mode; if not set by cfg
#if PALETTE_SYNCH
            pcs_ptr->palette_mode =
            (frm_hdr->allow_screen_content_tools) &&
#if MAR4_M3_ADOPTIONS
#if MAR10_ADOPTIONS
#if M8_PALETTE
#if UPGRADE_M6_M7_M8
#if M1_SC_ADOPTION
#if REVERT_WHITE // palette_mode
#if MAY19_ADOPTIONS
#if JUNE23_ADOPTIONS
            ((pcs_ptr->enc_mode <= ENC_M3) || (pcs_ptr->temporal_layer_index == 0 && pcs_ptr->enc_mode <= ENC_M8))
#else
#if NEW_M8
            ((pcs_ptr->enc_mode <= ENC_M3) || (pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M4) ||
            (pcs_ptr->temporal_layer_index == 0 && pcs_ptr->enc_mode <= ENC_M8))
#else
#if UNIFY_SC_NSC
            ((pcs_ptr->enc_mode <= ENC_M3) || (pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M4) ||
            (pcs_ptr->temporal_layer_index == 0 && pcs_ptr->enc_mode <= ENC_M7))
#else
#if JUNE17_ADOPTIONS
           ((pcs_ptr->enc_mode <= ENC_M3) || (pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M4) ||
           (pcs_ptr->temporal_layer_index == 0 && pcs_ptr->enc_mode <= ENC_M6))
#else
#if JUNE15_ADOPTIONS
           ((pcs_ptr->enc_mode <= ENC_M0) || (pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M4) ||
#else
#if PRESET_SHIFITNG
           (MR_MODE || (pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M4) ||
#else
           (MR_MODE || (pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M6) ||
#endif
#endif
#endif
#endif
#endif
#endif
#if !JUNE17_ADOPTIONS
#if PRESET_SHIFITNG
            (pcs_ptr->slice_type == I_SLICE && pcs_ptr->enc_mode <= ENC_M5))
#else
            (pcs_ptr->slice_type == I_SLICE && pcs_ptr->enc_mode <= ENC_M7))
#endif
#endif
#else
#if MAY16_7PM_ADOPTIONS
           ((pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M6) ||
            (pcs_ptr->slice_type == I_SLICE && pcs_ptr->enc_mode <= ENC_M7))
#else
#if APR25_3AM_ADOPTIONS
           (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M6) ||
            (pcs_ptr->slice_type == I_SLICE && pcs_ptr->enc_mode <= ENC_M7))
#else
           (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->is_used_as_reference_flag && pcs_ptr->enc_mode <= ENC_M7))
#endif
#endif
#endif
#else
           (pcs_ptr->enc_mode <= ENC_M0 || pcs_ptr->is_used_as_reference_flag)
#endif
#else
           (pcs_ptr->enc_mode <= ENC_M5 || pcs_ptr->is_used_as_reference_flag)
#endif
#else
            pcs_ptr->enc_mode <= ENC_M5
#endif
#else
            pcs_ptr->enc_mode <= ENC_M8
#endif
#else
            pcs_ptr->enc_mode <= ENC_M3
#endif
#else
            pcs_ptr->enc_mode <= ENC_M1
#endif

            ? 6 : 0;
#else
            pcs_ptr->palette_mode = pcs_ptr->enc_mode <= ENC_M1 ? 6 : 0;
#endif
        else
            pcs_ptr->palette_mode = scs_ptr->static_config.enable_palette;
    else
        pcs_ptr->palette_mode = 0;

    assert(pcs_ptr->palette_mode<7);

    if (!pcs_ptr->scs_ptr->static_config.disable_dlf_flag && frm_hdr->allow_intrabc == 0) {
#if MAR2_M8_ADOPTIONS
#if M8_LOOP_FILTER && !UPGRADE_M8 || REVERT_WHITE
#if REVERT_WHITE //  loop_filter_mode
#if PRESET_SHIFITNG
#if M6_LOOP_FILTER_MODE
#if UNIFY_SC_NSC
        if (pcs_ptr->enc_mode <= ENC_M6)
#else
        if (pcs_ptr->enc_mode <= ENC_M6 || pcs_ptr->sc_content_detected)
#endif
#else
        if (pcs_ptr->enc_mode <= ENC_M5 || pcs_ptr->sc_content_detected)
#endif
#else
        if (pcs_ptr->enc_mode <= ENC_M7 || pcs_ptr->sc_content_detected)
#endif
#else
        if (pcs_ptr->enc_mode <= ENC_M5)
#endif
            pcs_ptr->loop_filter_mode = 3;
        else
            pcs_ptr->loop_filter_mode =
            pcs_ptr->is_used_as_reference_flag ? 1 : 0;
#else
        pcs_ptr->loop_filter_mode = 3;
#endif
#else
        if (pcs_ptr->enc_mode <= ENC_M7)
            pcs_ptr->loop_filter_mode = 3;
        else
            pcs_ptr->loop_filter_mode = pcs_ptr->is_used_as_reference_flag ? 1 : 0;
#endif
    }
    else
        pcs_ptr->loop_filter_mode = 0;

    // CDEF Level                                   Settings
    // 0                                            OFF
    // 1                                            1 step refinement
    // 2                                            4 step refinement
    // 3                                            8 step refinement
    // 4                                            16 step refinement
    // 5                                            64 step refinement
    if (scs_ptr->seq_header.enable_cdef && frm_hdr->allow_intrabc == 0) {
        if (scs_ptr->static_config.cdef_mode == DEFAULT) {
#if MAR17_ADOPTIONS
#if M8_CDEF
#if !UNIFY_SC_NSC
        if (pcs_ptr->sc_content_detected)
#if UPGRADE_M6_M7_M8
#if PRESET_SHIFITNG
            if (pcs_ptr->enc_mode <= ENC_M5)
#else
            if (pcs_ptr->enc_mode <= ENC_M7)
#endif
                pcs_ptr->cdef_filter_mode = 5;
            else
#if M5_I_CDEF
                pcs_ptr->cdef_filter_mode = pcs_ptr->slice_type == I_SLICE ? 5 : 2;
#else
                pcs_ptr->cdef_filter_mode = 2;
#endif
#else
            pcs_ptr->cdef_filter_mode = 5;
#endif
        else
#endif
#if UPGRADE_M6_M7_M8
#if PRESET_SHIFITNG
            if (pcs_ptr->enc_mode <= ENC_M5)
#else
            if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M5)
#endif
                pcs_ptr->cdef_filter_mode = 5;
            else
#if M5_I_CDEF
                pcs_ptr->cdef_filter_mode = pcs_ptr->slice_type == I_SLICE ? 5 : 2;
#else
                pcs_ptr->cdef_filter_mode = 2;
#endif
#else
        pcs_ptr->cdef_filter_mode = 5;
#endif
#else
#if MAR10_ADOPTIONS
        if (pcs_ptr->sc_content_detected)
            pcs_ptr->cdef_filter_mode = 5;
        else
#endif
                if (pcs_ptr->enc_mode <= ENC_M7)
                    pcs_ptr->cdef_filter_mode = 5;
                else
                    pcs_ptr->cdef_filter_mode = 2;
#endif
        } else
            pcs_ptr->cdef_filter_mode = (int8_t)(scs_ptr->static_config.cdef_mode);
    }
    else
        pcs_ptr->cdef_filter_mode = 0;

    // SG Level                                    Settings
    // 0                                            OFF
    // 1                                            0 step refinement
    // 2                                            1 step refinement
    // 3                                            4 step refinement
    // 4                                            16 step refinement

    Av1Common* cm = pcs_ptr->av1_cm;
    if (scs_ptr->static_config.sg_filter_mode == DEFAULT) {
#if UNIFY_SC_NSC
        if (pcs_ptr->enc_mode <= ENC_M3)
            cm->sg_filter_mode = 4;
        else
            cm->sg_filter_mode = pcs_ptr->slice_type == I_SLICE ? 4 : 1;
#else
        if (sc_content_detected)
#if MAR12_M8_ADOPTIONS
#if M8_SG
#if UPGRADE_M6_M7_M8
#if JUNE17_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if PRESET_SHIFITNG
            if (pcs_ptr->enc_mode <= ENC_M5)
#else
            if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M5)
#endif
                cm->sg_filter_mode = 4;
            else
#if M5_I_SG
                cm->sg_filter_mode = pcs_ptr->slice_type == I_SLICE ? 4 : 1;
#else
                cm->sg_filter_mode = 1;
#endif
#else
            cm->sg_filter_mode = 4;
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M5)
                cm->sg_filter_mode = 4;
            else
                cm->sg_filter_mode = 0;
#endif
#if MAY19_ADOPTIONS
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M3)
#else
    else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#else
#if PRESETS_SHIFT
    else if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if MAR17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7)
#else
#if MAR10_ADOPTIONS
#if MAR11_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M4)
#else
    else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#endif
#endif
#endif
            cm->sg_filter_mode = 4;
#if MAY19_ADOPTIONS
    else
        cm->sg_filter_mode = pcs_ptr->slice_type == I_SLICE ? 4 : 1;
#else
#if MAR12_M8_ADOPTIONS
#if M8_SG
#if UPGRADE_M6_M7_M8
#if APR25_3AM_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M6)
#else
    else if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#else
    else if (pcs_ptr->enc_mode <= ENC_M5)
#endif
        cm->sg_filter_mode = 3;
    else
#if M5_I_SG
        cm->sg_filter_mode = pcs_ptr->slice_type == I_SLICE ? 3 : 1;
#else
        cm->sg_filter_mode = 1;
#endif
#else
    else
        cm->sg_filter_mode = 3;
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M6)
            cm->sg_filter_mode = 3;
        else
            cm->sg_filter_mode = 1;
#endif
#endif
#endif
    }
    else
        cm->sg_filter_mode = scs_ptr->static_config.sg_filter_mode;

    // WN Level                                     Settings
    // 0                                            OFF
    // 1                                            3-Tap luma/ 3-Tap chroma
    // 2                                            5-Tap luma/ 5-Tap chroma
    // 3                                            7-Tap luma/ 5-Tap chroma

    if (scs_ptr->static_config.wn_filter_mode == DEFAULT) {
#if !UNIFY_SC_NSC
        if (sc_content_detected)
#if MAR12_M8_ADOPTIONS
            cm->wn_filter_mode = 3;
#else
#if MAR3_M6_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M6)
#else
            if (pcs_ptr->enc_mode <= ENC_M5)
#endif
                cm->wn_filter_mode = 3;
            else
                cm->wn_filter_mode = 0;
#endif
        else
#endif
#if MAY19_ADOPTIONS
#if JUNE17_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M5)
#else
#if PRESET_SHIFITNG
        if (pcs_ptr->enc_mode <= ENC_M4)
#else
        if (pcs_ptr->enc_mode <= ENC_M6)
#endif
#endif
#else
#if PRESETS_SHIFT
        if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if MAR17_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M7)
#else
            if (pcs_ptr->enc_mode <= ENC_M5)
#endif
#endif
#endif
            cm->wn_filter_mode = 3;
#if MAR12_M8_ADOPTIONS
        else
            cm->wn_filter_mode = 2;
#else
        else if (pcs_ptr->enc_mode <= ENC_M7)
            cm->wn_filter_mode = 2;
        else
            cm->wn_filter_mode = 0;
#endif
    }
    else
        cm->wn_filter_mode = scs_ptr->static_config.wn_filter_mode;

    // Intra prediction modes                       Settings
    // 0                                            FULL
    // 1                                            LIGHT per block : disable_z2_prediction && disable_angle_refinement  for 64/32/4
    // 2                                            OFF per block : disable_angle_prediction for 64/32/4
    // 3                                            OFF : disable_angle_prediction
    // 4                                            OIS based Intra
    // 5                                            Light OIS based Intra

    if (pcs_ptr->slice_type == I_SLICE)
#if MAR2_M8_ADOPTIONS
        pcs_ptr->intra_pred_mode = 0;
#else
    if (sc_content_detected)
        if (pcs_ptr->enc_mode <= ENC_M6)
            pcs_ptr->intra_pred_mode = 0;
        else
            pcs_ptr->intra_pred_mode = 4;
    else
        if (pcs_ptr->enc_mode <= ENC_M7)
            pcs_ptr->intra_pred_mode = 0;
        else
            pcs_ptr->intra_pred_mode = 4;
#endif
    else {
#if !UNIFY_SC_NSC
    if (sc_content_detected)
#if JUNE17_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M3)
#else
#if JUNE8_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M2)
#else
#if APR24_M3_ADOPTIONS
#if PRESET_SHIFITNG
            if (pcs_ptr->enc_mode <= ENC_M1)
#else
            if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#else
#if APR23_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M5)
#else
#if PRESETS_SHIFT
            if (pcs_ptr->enc_mode <= ENC_M2)
#else
#if MAR12_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M3)
#else
#if MAR10_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M2)
#else
        if (pcs_ptr->enc_mode <= ENC_M1)
#endif
#endif
#endif
#endif
#endif
#endif
#endif
            pcs_ptr->intra_pred_mode = 0;
#if !MAR12_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M2)
            if (pcs_ptr->temporal_layer_index == 0)
                pcs_ptr->intra_pred_mode = 1;
            else
                pcs_ptr->intra_pred_mode = 2;
#endif
#if MAR17_ADOPTIONS
            else
                if (pcs_ptr->temporal_layer_index == 0)
                    pcs_ptr->intra_pred_mode = 2;
                else
                    pcs_ptr->intra_pred_mode = 3;
#else
        else if (pcs_ptr->enc_mode <= ENC_M6)
            if (pcs_ptr->temporal_layer_index == 0)
                pcs_ptr->intra_pred_mode = 2;
            else
                pcs_ptr->intra_pred_mode = 3;
        else
            pcs_ptr->intra_pred_mode = 4;
#endif
#if JUNE8_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M2)
#else
#if APR23_ADOPTIONS_2
#if PRESET_SHIFITNG
        else if (pcs_ptr->enc_mode <= ENC_M1)
#else
        else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#else
#if PRESETS_SHIFT
        else if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if MAR20_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M7)
#else
#if MAR10_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M4)
#else
    else
        if ((pcs_ptr->enc_mode <= ENC_M1) || (pcs_ptr->enc_mode <= ENC_M2 && pcs_ptr->temporal_layer_index == 0))
#endif
#endif
#endif
#endif
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M2)
#endif
            pcs_ptr->intra_pred_mode = 0;
#if MAR2_M8_ADOPTIONS
#if M8_INTRA_MODE && !UPGRADE_M8
        else if (pcs_ptr->enc_mode <= ENC_M5)
            if (pcs_ptr->temporal_layer_index == 0)
                pcs_ptr->intra_pred_mode = 1;
            else
                pcs_ptr->intra_pred_mode = 3;
        else
            pcs_ptr->intra_pred_mode = 3;
#else
        else
            if (pcs_ptr->temporal_layer_index == 0)
                pcs_ptr->intra_pred_mode = 1;
            else
                pcs_ptr->intra_pred_mode = 3;
#endif
#else
        else if (pcs_ptr->enc_mode <= ENC_M7)
            if (pcs_ptr->temporal_layer_index == 0)
                pcs_ptr->intra_pred_mode = 1;
            else
                pcs_ptr->intra_pred_mode = 3;
        else
            pcs_ptr->intra_pred_mode = 4;
#endif
    }

        if (MR_MODE)
            pcs_ptr->intra_pred_mode = 0;

        // Set tx size search mode      Settings
        // 0                 OFF: no transform partitioning
        // 1                 ON for INTRA blocks
#if !UNIFY_SC_NSC
#if MAR4_M6_ADOPTIONS
    if (pcs_ptr->sc_content_detected)
#if NEW_TXS_SETTINGS
#if JUNE17_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if JUNE11_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M3)
#else
        if (pcs_ptr->enc_mode <= ENC_M0)
#endif
#endif
            pcs_ptr->tx_size_search_mode = 1;
        else
#endif
#if MAR10_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M8)
#else
        if (pcs_ptr->enc_mode <= ENC_M5)
#endif
            pcs_ptr->tx_size_search_mode = (pcs_ptr->slice_type == I_SLICE) ? 1 : 0;
        else
            pcs_ptr->tx_size_search_mode = 0;
#if JUNE17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if JUNE11_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M3)
#else
#if PRESETS_SHIFT
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M2)
#else
    else if (pcs_ptr->enc_mode <= ENC_M4)
#endif
#else
#if MAR17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7)
#else
#if MAR10_ADOPTIONS
#if MAR11_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M4)
#else
    else if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#else
    else if (pcs_ptr->enc_mode <= ENC_M3)
#endif
#endif
#endif
#endif
#endif
        pcs_ptr->tx_size_search_mode = 1;
#else
        if (pcs_ptr->enc_mode <= ENC_M2)
            pcs_ptr->tx_size_search_mode = (pcs_ptr->sc_content_detected && pcs_ptr->slice_type != I_SLICE) ? 0 : 1;
#endif
#else
#if JUNE26_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if JUNE25_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M5)
#else
    if (pcs_ptr->enc_mode <= ENC_M4)
#endif
#endif
        pcs_ptr->tx_size_search_mode = 1;
#endif
#if NEW_TXS_SETTINGS
#if NEW_M8
    else if (pcs_ptr->enc_mode <= ENC_M8)
#else
#if JUNE17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7)
#else
    else if (pcs_ptr->enc_mode <= ENC_M6)
#endif
#endif
        pcs_ptr->tx_size_search_mode = (pcs_ptr->slice_type == I_SLICE) ? 1 : 0;
#endif
        else
            pcs_ptr->tx_size_search_mode = 0;

#if APR22_ADOPTIONS
    // Assign whether to use TXS in inter classes (if TXS is ON)
    // 0 OFF - TXS in intra classes only
    // 1 ON - TXS in all classes
    // 2 ON - INTER TXS restricted to max 1 depth
#if RESTRICT_INTER_TXS_DEPTH
#if UNIFY_SC_NSC
#if REMOVE_MR_MACRO
    if (pcs_ptr->enc_mode <= ENC_MRS)
#else
    if (MRS_MODE)
#endif
#else
#if JUNE15_ADOPTIONS
    if (MRS_MODE || (pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M0))
#else
    if (MR_MODE)
#endif
#endif
        pcs_ptr->txs_in_inter_classes = 1;
#if UNIFY_SC_NSC
    else if (pcs_ptr->enc_mode <= ENC_M0)
#else
#if JUNE17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->sc_content_detected && pcs_ptr->enc_mode <= ENC_M3))
#else
    else if (pcs_ptr->enc_mode <= ENC_M0)
#endif
#endif
        pcs_ptr->txs_in_inter_classes = 2;
    else
        pcs_ptr->txs_in_inter_classes = 0;
#else
#if MAY16_7PM_ADOPTIONS
    pcs_ptr->txs_in_inter_classes = 0;
#else
    if (pcs_ptr->sc_content_detected)
        pcs_ptr->txs_in_inter_classes = 0;
    else if (pcs_ptr->enc_mode <= ENC_M0)
        pcs_ptr->txs_in_inter_classes = (pcs_ptr->is_used_as_reference_flag && pcs_ptr->input_resolution <= INPUT_SIZE_480p_RANGE) ? 1 : 0;
    else
        pcs_ptr->txs_in_inter_classes = 0;
#endif
#endif
#endif
#if !INTER_COMP_REDESIGN
        // Set Wedge mode      Settings
        // 0                 FULL: Full search
        // 1                 Fast: If two predictors are very similar, skip wedge compound mode search
        // 2                 Fast: estimate Wedge sign
        // 3                 Fast: Mode 1 & Mode 2
#if WEDGE_SYNCH
#if MAR11_ADOPTIONS
    if (pcs_ptr->sc_content_detected)
#if MAR12_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M3)
#else
        if (pcs_ptr->enc_mode <= ENC_M1)
#endif
            pcs_ptr->wedge_mode = 0;
        else
            pcs_ptr->wedge_mode = 1;
#if MAR17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M7)
#else
    else if (pcs_ptr->enc_mode <= ENC_M4)
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M2)
#endif
        pcs_ptr->wedge_mode = 0;
    else
        pcs_ptr->wedge_mode = 1;
#else
        pcs_ptr->wedge_mode = 0;
#endif
#if M8_NEW_REF
    pcs_ptr->wedge_mode = 0;
#endif
#endif
        // inter intra pred                      Settings
        // 0                                     OFF
        // 1                                     ON
        //picture level switch,  has to follow the sequence level.
        if (pcs_ptr->slice_type != I_SLICE && scs_ptr->seq_header.enable_interintra_compound) {
#if INTRA_COMPOUND_OPT
    // inter intra pred                      Settings
    // 0                                     OFF
    // 1                                     FULL
    // 2                                     FAST 1 : Do not inject for non basic inter
    // 3                                     FAST 2 : 1 + MRP pruning/ similar based disable + NIC tuning
#if PRESETS_SHIFT
#if !UNIFY_SC_NSC
#if ENABLE_SC_DETECTOR
        if (pcs_ptr->sc_content_detected)
            pcs_ptr->enable_inter_intra = 0;
        else
#endif
#endif
#if JUNE8_ADOPTIONS
        pcs_ptr->enable_inter_intra = pcs_ptr->enc_mode <= ENC_M2 ? 2 : 3;
#else
#if PRESET_SHIFITNG
        pcs_ptr->enable_inter_intra = pcs_ptr->enc_mode <= ENC_M1 ? 2 : 3;
#else
        pcs_ptr->enable_inter_intra = pcs_ptr->enc_mode <= ENC_M2 ? 2 : 3;
#endif
#endif
#else
#if MAR30_ADOPTIONS
        pcs_ptr->enable_inter_intra = pcs_ptr->enc_mode <= ENC_M3 ? 2 : 3;
#else
        pcs_ptr->enable_inter_intra= pcs_ptr->enc_mode <= ENC_M2 ? 2 : 3;
#endif
#endif
#else
            pcs_ptr->enable_inter_intra = 1;//shut for sc , if needed.
#endif
        }
        else {
            pcs_ptr->enable_inter_intra = 0;
        }
        // Set compound mode      Settings
        // 0                 OFF: No compond mode search : AVG only
        // 1                 ON: compond mode search: AVG/DIST/DIFF
        // 2                 ON: AVG/DIST/DIFF/WEDGE
        if (scs_ptr->static_config.compound_level == DEFAULT) {
            if (scs_ptr->compound_mode)
#if MAR17_ADOPTIONS
#if INTER_COMP_REDESIGN
#if PRESETS_SHIFT
#if !UNIFY_SC_NSC
#if APR23_ADOPTIONS
            if (pcs_ptr->sc_content_detected)
#if JUNE17_ADOPTIONS
                pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M0 ? 1 :
                                            pcs_ptr->enc_mode <= ENC_M3 ? 3 : 0;
#else
#if JUNE15_ADOPTIONS
                pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M0 ? 1 :
                                            pcs_ptr->enc_mode <= ENC_M2 ? 3 : 0;
#else
#if NEW_MRP_SETTINGS
                pcs_ptr->compound_mode = MR_MODE ? 1 :
                                            pcs_ptr->enc_mode <= ENC_M0 ? 2 :
                                            pcs_ptr->enc_mode <= ENC_M2 ? 3 : 0;
#else
#if JUNE8_ADOPTIONS
                pcs_ptr->compound_mode = MR_MODE ? 1 : pcs_ptr->enc_mode <= ENC_M2 ? 2 : 0;
#else
#if SHIFT_M5_SC_TO_M3
#if PRESET_SHIFITNG
                pcs_ptr->compound_mode = MR_MODE ? 1 : pcs_ptr->enc_mode <= ENC_M1 ? 2 : 0;
#else
                pcs_ptr->compound_mode = MR_MODE ? 1 : pcs_ptr->enc_mode <= ENC_M2 ? 2 : 0;
#endif
#else
                pcs_ptr->compound_mode = MR_MODE ? 1 : pcs_ptr->enc_mode <= ENC_M4 ? 2 : 0;
#endif
#endif
#endif
#endif
#endif
            else
#endif
#endif
#if M2_COMBO_1 || M1_COMBO_2 || NEW_M1_CAND
#if UNIFY_SC_NSC
                pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M1 ? 1 :
                                            pcs_ptr->enc_mode <= ENC_M3 ? 3 : 0;
#else
#if JUNE17_ADOPTIONS
                pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M1 ? 1 :
                                            pcs_ptr->enc_mode <= ENC_M2 ? 3 : 0;
#else
#if NEW_MRP_SETTINGS
                pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M0 ? 1 :
                                            pcs_ptr->enc_mode <= ENC_M2 ? 3 : 0;
#else
                pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M0 ? 1 :
#if PRESET_SHIFITNG
                pcs_ptr->enc_mode <= ENC_M2 ? 2 : 0;
#else
                pcs_ptr->enc_mode <= ENC_M4 ? 2 : 0;
#endif
#endif
#endif
#endif
#else
            pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M3 ? 1 :
                                     pcs_ptr->enc_mode <= ENC_M4 ? 2 : 0;
#endif
#else
            pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M4 ? 1 :
                                     pcs_ptr->enc_mode <= ENC_M7 ? 2 : 0;
#endif
#else
            pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M7 ? 2 : 1;
#endif
#else
#if MAR4_M3_ADOPTIONS
#if MAR10_ADOPTIONS
            pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M4 ? 2 : 1;
#else
            pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M3 ? 2 : 1;
#endif
#else
                    pcs_ptr->compound_mode = pcs_ptr->enc_mode <= ENC_M1 ? 2 : 1;
#endif
#endif
            else
                pcs_ptr->compound_mode = 0;
        }
        else
            pcs_ptr->compound_mode = scs_ptr->static_config.compound_level;
#if !WEDGE_SYNCH
        if (pcs_ptr->wedge_mode > 0 && pcs_ptr->compound_mode != 2)
            SVT_LOG("wedge_mode set but will not be active\n");
#endif
        // Set frame end cdf update mode      Settings
        // 0                                     OFF
        // 1                                     ON
        if (scs_ptr->static_config.frame_end_cdf_update == DEFAULT)
            pcs_ptr->frame_end_cdf_update_mode = 1;
        else
            pcs_ptr->frame_end_cdf_update_mode = scs_ptr->static_config.frame_end_cdf_update;

#if !SHUT_ME_CAND_SORTING
        if (scs_ptr->static_config.prune_unipred_me == DEFAULT)
#if MAR4_M6_ADOPTIONS
#if MAR10_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M8)
#else
        if (pcs_ptr->enc_mode <= ENC_M5)
#endif
            pcs_ptr->prune_unipred_at_me = 1;
        else
            pcs_ptr->prune_unipred_at_me = 0;
#else
            if (pcs_ptr->enc_mode >= ENC_M4)
                pcs_ptr->prune_unipred_at_me = 0;
            else
                pcs_ptr->prune_unipred_at_me = 1;
#endif
        else
            pcs_ptr->prune_unipred_at_me = scs_ptr->static_config.prune_unipred_me;
#endif
        //CHKN: Temporal MVP should be disabled for pictures beloning to 4L MiniGop preceeded by 5L miniGOP. in this case the RPS is wrong(known issue). check RPS construction for more info.
#if !UNIFY_SC_NSC
#if ENABLE_SC_DETECTOR
        if (pcs_ptr->sc_content_detected)
            pcs_ptr->frm_hdr.use_ref_frame_mvs = 0;
        else
#endif
#endif
        if (pcs_ptr->slice_type == I_SLICE)
            pcs_ptr->frm_hdr.use_ref_frame_mvs = 0;
        else
            pcs_ptr->frm_hdr.use_ref_frame_mvs = scs_ptr->mfmv_enabled;

        // Global motion level                        Settings
        // GM_FULL                                    Exhaustive search mode.
        // GM_DOWN                                    Downsampled resolution with a downsampling factor of 2 in each dimension
        // GM_TRAN_ONLY                               Translation only using ME MV.
#if MAR17_ADOPTIONS
#if MAR20_M4_ADOPTIONS
#if MAR30_ADOPTIONS
#if APR08_ADOPTIONS
#if PRESETS_SHIFT
#if APR22_ADOPTIONS
#if SHIFT_M5_SC_TO_M3
#if JUNE17_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M2)
#else
#if PRESET_SHIFITNG
    if (pcs_ptr->enc_mode <= ENC_M1)
#else
    if (pcs_ptr->enc_mode <= ENC_M2)
#endif
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M2 || (pcs_ptr->enc_mode <= ENC_M4 && pcs_ptr->sc_content_detected))
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->enc_mode <= ENC_M4 && pcs_ptr->sc_content_detected))
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M0 || (pcs_ptr->enc_mode <= ENC_M7 && pcs_ptr->sc_content_detected))
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M3 || (pcs_ptr->enc_mode <= ENC_M7 && pcs_ptr->sc_content_detected))
#endif
#else
    if (pcs_ptr->enc_mode <= ENC_M3)
#endif
#else
#if MAR18_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M4)
#else
    if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#endif
        pcs_ptr->gm_level = GM_FULL;
#if JUNE26_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M4)
        pcs_ptr->gm_level = GM_DOWN;
    else
        pcs_ptr->gm_level = GM_DOWN16;
#else
    else
        pcs_ptr->gm_level = GM_DOWN;
#endif
#else
#if MAR12_ADOPTIONS
    if (pcs_ptr->sc_content_detected)
        if (pcs_ptr->enc_mode <= ENC_M3)
            pcs_ptr->gm_level = GM_FULL;
        else
            pcs_ptr->gm_level = GM_DOWN;
    else
#endif
#if MAR10_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M1)
#else
        if (pcs_ptr->enc_mode <= ENC_M0 ||
        (pcs_ptr->enc_mode <= ENC_M1 && pcs_ptr->sc_content_detected))
#endif
            pcs_ptr->gm_level = GM_FULL;
        else
            pcs_ptr->gm_level = GM_DOWN;
#endif
        //Exit TX size search when all coefficients are zero
        // 0: OFF
        // 1: ON
        pcs_ptr->tx_size_early_exit = 1;
        //Prune reference and reduce ME SR based on HME/ME distortion
        // 0: OFF
        // 1: ON
#if !ME_HME_PRUNING_CLEANUP
#if MAR30_ADOPTIONS
    if (MR_MODE || pcs_ptr->sc_content_detected)
        pcs_ptr->prune_ref_based_me = 0;
    else
        pcs_ptr->prune_ref_based_me = 1;
#else
#if MAR20_ADOPTIONS
    if (MR_MODE)
        pcs_ptr->prune_ref_based_me = 0;
    else
        pcs_ptr->prune_ref_based_me = 1;
#else
        if (pcs_ptr->sc_content_detected)
            pcs_ptr->prune_ref_based_me = 0;
        else
            pcs_ptr->prune_ref_based_me = 1;
#endif
#endif
#endif
#if TF_LEVELS
#if UNIFY_SC_NSC
    uint8_t perform_filtering =
        (scs_ptr->enable_altrefs == EB_TRUE && scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS && scs_ptr->static_config.hierarchical_levels >= 1)
        ? 1 : 0;
#else
#if UPGRADE_M6_M7_M8
    uint8_t perform_filtering =
        (scs_ptr->enable_altrefs == EB_TRUE && scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS && pcs_ptr->sc_content_detected == 0 && scs_ptr->static_config.hierarchical_levels >= 1)
        ? 1 : 0;
#else
    uint8_t perform_filtering =
        (scs_ptr->enable_altrefs == EB_TRUE && scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS && pcs_ptr->sc_content_detected == 0 && scs_ptr->static_config.hierarchical_levels >= 1) &&
        (pcs_ptr->temporal_layer_index == 0 ||(pcs_ptr->temporal_layer_index == 1 && scs_ptr->static_config.hierarchical_levels >= 3))
        ? 1 : 0;
#endif
#endif


#if IMPROVED_TF_LEVELS
    if (perform_filtering) {
#if JUNE26_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if JUNE25_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M4) {
#else
        if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
#endif
            if (pcs_ptr->temporal_layer_index == 0 || (pcs_ptr->temporal_layer_index == 1 && scs_ptr->static_config.hierarchical_levels >= 3))
                context_ptr->tf_level = 1;
            else
                context_ptr->tf_level = 0;
        }
        else if (pcs_ptr->enc_mode <= ENC_M6) {
            if (pcs_ptr->temporal_layer_index == 0 || (pcs_ptr->temporal_layer_index == 1 && scs_ptr->static_config.hierarchical_levels >= 3))
                context_ptr->tf_level = 2;
            else
                context_ptr->tf_level = 0;
        }
        else {
            if (pcs_ptr->temporal_layer_index == 0 || (pcs_ptr->temporal_layer_index == 1 && scs_ptr->static_config.hierarchical_levels >= 3))
                context_ptr->tf_level = 3;
            else
                context_ptr->tf_level = 0;
        }
    }
    else
        context_ptr->tf_level = 0;
#else
    if (perform_filtering) {
#if UPGRADE_M6_M7_M8
#if PRESET_SHIFITNG
        if (pcs_ptr->enc_mode <= ENC_M5) {
#else
        if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
            if (pcs_ptr->temporal_layer_index == 0 || (pcs_ptr->temporal_layer_index == 1 && scs_ptr->static_config.hierarchical_levels >= 3))
                context_ptr->tf_level = 0;
            else
                context_ptr->tf_level = 3;
        }
        else {
#if ADD_SKIP_INTRA_SIGNAL
            if (pcs_ptr->enc_mode <= ENC_M6) {
                if (pcs_ptr->temporal_layer_index == 0)
                    context_ptr->tf_level = 0;
                else
                    context_ptr->tf_level = 3;
            }
            else
                context_ptr->tf_level = 3;
#else
            if (pcs_ptr->temporal_layer_index == 0)
                context_ptr->tf_level = 0;
            else
                context_ptr->tf_level = 3;
#endif
        }
#else
#if UPGRADE_M8
        if (pcs_ptr->enc_mode <= ENC_M8) {
#else
        if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
            context_ptr->tf_level = 0;
        }
        else {
            context_ptr->tf_level = 2;
        }
#endif
    }
    else
        context_ptr->tf_level = 3;
#endif
    set_tf_controls(context_ptr, context_ptr->tf_level);
#endif
#if ON_OFF_FEATURE_MRP
    // MRP control
    // 0: OFF (1,1)  ; override features
    // 1: FULL (4,3) ; override features
    // 2: (4,3) ; No-override features
    // 3: (3,3) ; No-override features
    // 4: (3,2) ; No-override features
    // 5: (2,3) ; No-override features
    // 6: (2,2) ; No-override features
    // 7: (2,1) ; No-override features
    // 8: (1,2) ; No-override features
    // 9: (1,1) ; No-override features

    // Level 0 , 1  : set ref_list0_count_try and ref_list1_count_try and Override MRP-related features
    // Level 2 .. 9 : Only set ref_list0_count_try and ref_list1_count_try

    if (scs_ptr->static_config.mrp_level == DEFAULT) {
#if !UNIFY_SC_NSC
        if (pcs_ptr->sc_content_detected)
            if (MRS_MODE)
                pcs_ptr->mrp_level = 3;
            else if (MR_MODE)
                pcs_ptr->mrp_level = 4;
            else if (pcs_ptr->enc_mode <= ENC_M4)
                pcs_ptr->mrp_level = 6;
            else if (pcs_ptr->enc_mode <= ENC_M5)
                pcs_ptr->mrp_level = 7;
            else
                pcs_ptr->mrp_level = 9;
        else
#endif
#if JUNE26_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M6)
                pcs_ptr->mrp_level = 2;
#else
            if (pcs_ptr->enc_mode <= ENC_M5)
                pcs_ptr->mrp_level = 2;
#if NEW_M7_MRP
            else if (pcs_ptr->enc_mode <= ENC_M6)
#else
            else if (pcs_ptr->enc_mode <= ENC_M7)
#endif
                pcs_ptr->mrp_level = 6;
#endif
            else
                pcs_ptr->mrp_level = pcs_ptr->is_used_as_reference_flag  ? 6 : 9;
    }else
        pcs_ptr->mrp_level = scs_ptr->static_config.mrp_level;
#endif
#if TPL_OPT
    if (pcs_ptr->enc_mode <= ENC_M5)
        pcs_ptr->tpl_opt_flag = 0;
    else
        pcs_ptr->tpl_opt_flag = 1;
#endif
    return return_error;
}

int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
//set the ref frame types used for this picture,
#if  !REMOVE_MRP_MODE || !MRP_CTRL
static void set_all_ref_frame_type(SequenceControlSet *scs_ptr, PictureParentControlSet  *parent_pcs_ptr, MvReferenceFrame ref_frame_arr[], uint8_t* tot_ref_frames)
#else
static void set_all_ref_frame_type(PictureParentControlSet  *parent_pcs_ptr, MvReferenceFrame ref_frame_arr[], uint8_t* tot_ref_frames)
#endif
{
    MvReferenceFrame rf[2];
    *tot_ref_frames = 0;

    //SVT_LOG("POC %i  totRef L0:%i   totRef L1: %i\n", parent_pcs_ptr->picture_number, parent_pcs_ptr->ref_list0_count, parent_pcs_ptr->ref_list1_count);

     //single ref - List0
#if ON_OFF_FEATURE_MRP
    for (uint8_t ref_idx0 = 0; ref_idx0 < parent_pcs_ptr->mrp_ctrls.ref_list0_count_try; ++ref_idx0) {
#else
    for (uint8_t ref_idx0 = 0; ref_idx0 < parent_pcs_ptr->ref_list0_count_try; ++ref_idx0) {
#endif
        rf[0] = svt_get_ref_frame_type(REF_LIST_0, ref_idx0);
        ref_frame_arr[(*tot_ref_frames)++] = rf[0];
    }

    //single ref - List1
#if ON_OFF_FEATURE_MRP
    for (uint8_t ref_idx1 = 0; ref_idx1 < parent_pcs_ptr->mrp_ctrls.ref_list1_count_try; ++ref_idx1) {
#else
    for (uint8_t ref_idx1 = 0; ref_idx1 < parent_pcs_ptr->ref_list1_count_try; ++ref_idx1) {
#endif
        rf[1] = svt_get_ref_frame_type(REF_LIST_1, ref_idx1);
        ref_frame_arr[(*tot_ref_frames)++] = rf[1];
    }

    //compound Bi-Dir
#if ON_OFF_FEATURE_MRP
    for (uint8_t ref_idx0 = 0; ref_idx0 < parent_pcs_ptr->mrp_ctrls.ref_list0_count_try; ++ref_idx0) {
        for (uint8_t ref_idx1 = 0; ref_idx1 < parent_pcs_ptr->mrp_ctrls.ref_list1_count_try; ++ref_idx1) {
#else
    for (uint8_t ref_idx0 = 0; ref_idx0 < parent_pcs_ptr->ref_list0_count_try; ++ref_idx0) {
        for (uint8_t ref_idx1 = 0; ref_idx1 < parent_pcs_ptr->ref_list1_count_try; ++ref_idx1) {
#endif
            rf[0] = svt_get_ref_frame_type(REF_LIST_0, ref_idx0);
            rf[1] = svt_get_ref_frame_type(REF_LIST_1, ref_idx1);
            ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
        }
    }
#if  REMOVE_MRP_MODE
    if (parent_pcs_ptr->slice_type == B_SLICE)
#else
    if (scs_ptr->mrp_mode == 0 && parent_pcs_ptr->slice_type == B_SLICE)
#endif
    {

        //compound Uni-Dir
#if ON_OFF_FEATURE_MRP
        if (parent_pcs_ptr->mrp_ctrls.ref_list0_count_try > 1) {
#else
        if (parent_pcs_ptr->ref_list0_count_try > 1) {
#endif
            rf[0] = LAST_FRAME;
            rf[1] = LAST2_FRAME;
            ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
#if ON_OFF_FEATURE_MRP
            if (parent_pcs_ptr->mrp_ctrls.ref_list0_count_try > 2) {
#else
            if (parent_pcs_ptr->ref_list0_count_try > 2) {
#endif
                rf[1] = LAST3_FRAME;
                ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
#if ON_OFF_FEATURE_MRP
                if (parent_pcs_ptr->mrp_ctrls.ref_list0_count_try > 3) {
#else
                if (parent_pcs_ptr->ref_list0_count_try > 3) {
#endif
                    rf[1] = GOLDEN_FRAME;
                    ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
                }
            }
        }
#if ON_OFF_FEATURE_MRP
        if (parent_pcs_ptr->mrp_ctrls.ref_list0_count_try > 2) {
#else
        if (parent_pcs_ptr->ref_list1_count_try > 2) {
#endif
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
static EbBool set_frame_display_params(
        PictureParentControlSet       *pcs_ptr,
        PictureDecisionContext        *context_ptr,
        uint32_t                       mini_gop_index)
{
    Av1RpsNode *av1_rps = &pcs_ptr->av1_ref_signal;
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;

    if (pcs_ptr->pred_struct_ptr->pred_type == EB_PRED_LOW_DELAY_P || pcs_ptr->is_overlay) {
        //P frames
        av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_poc_array[BWD] = av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[LAST];

        frm_hdr->show_frame = EB_TRUE;
        pcs_ptr->has_show_existing = EB_FALSE;
    } else if (pcs_ptr->pred_struct_ptr->pred_type == EB_PRED_LOW_DELAY_B) {
        av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[LAST2];
        av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST3];
        av1_rps->ref_poc_array[BWD] = av1_rps->ref_poc_array[LAST];
        av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[LAST2];
        av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[LAST3];

        frm_hdr->show_frame = EB_TRUE;
        pcs_ptr->has_show_existing = EB_FALSE;
    } else if (pcs_ptr->pred_struct_ptr->pred_type == EB_PRED_RANDOM_ACCESS) {
        //Decide on Show Mecanism
        if (pcs_ptr->slice_type == I_SLICE) {
            //3 cases for I slice:  1:Key Frame treated above.  2: broken MiniGop due to sc or intra refresh  3: complete miniGop due to sc or intra refresh
            if (context_ptr->mini_gop_length[mini_gop_index] < pcs_ptr->pred_struct_ptr->pred_struct_period) {
                //Scene Change that breaks the mini gop and switch to LDP (if I scene change happens to be aligned with a complete miniGop, then we do not break the pred structure)
                frm_hdr->show_frame = EB_TRUE;
                pcs_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_FALSE;
                pcs_ptr->has_show_existing = EB_FALSE;
            }
        } else {
            if (context_ptr->mini_gop_length[mini_gop_index] != pcs_ptr->pred_struct_ptr->pred_struct_period) {
                SVT_LOG("Error in GOP indexing3\n");
            }
            // Handle b frame of Random Access out
            return EB_FALSE;
        }
    }
    return EB_TRUE;
}

static void set_key_frame_rps(PictureParentControlSet *pcs, PictureDecisionContext *context_ptr)
{
    FrameHeader *frm_hdr = &pcs->frm_hdr;
    context_ptr->lay0_toggle = 0;
    context_ptr->lay1_toggle = 0;
    context_ptr->lay2_toggle = 0;

    frm_hdr->show_frame = EB_TRUE;
    pcs->has_show_existing = EB_FALSE;
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
*Note: the  SceneChange I has pred_type = EB_PRED_RANDOM_ACCESS. if SChange is aligned on the miniGop,
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

    SequenceControlSet *scs_ptr = (SequenceControlSet*)pcs_ptr->scs_wrapper_ptr->object_ptr;
    PredictionStructureEntry *pred_position_ptr = pcs_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[pcs_ptr->pred_struct_index];
    //Set frame type
    if (pcs_ptr->slice_type == I_SLICE)
        frm_hdr->frame_type = pcs_ptr->idr_flag ? KEY_FRAME : INTRA_ONLY_FRAME;
    else
        frm_hdr->frame_type = INTER_FRAME;

    pcs_ptr->intra_only = pcs_ptr->slice_type == I_SLICE ? 1 : 0;

    if (pcs_ptr->hierarchical_levels == 0) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(pcs_ptr, context_ptr);
            return;
        }

        const uint8_t  base0_idx = context_ptr->lay0_toggle == 0 ? 1 : context_ptr->lay0_toggle == 1 ? 2 : 0;
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 0 : 1;

        //{1, 2, 0, 0},     // GOP Index 0 - Ref List 0
        //{1, 2, 0, 0},     // GOP Index 0 - Ref List 1
        av1_rps->ref_dpb_index[LAST] = base1_idx;
        av1_rps->ref_dpb_index[LAST2] = base0_idx;
        av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

        av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[LAST2];
        av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
        gop_i = 0;
        av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list0[0]);
        av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list0[1]);
        av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
        av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

        av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list1[0]);
        av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, flat_pred_struct[gop_i].ref_list1[1]);
        av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

        prune_refs(pred_position_ptr, av1_rps);
        av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;

        // Flat mode, output all frames
        set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index);
        frm_hdr->show_frame = EB_TRUE;
        pcs_ptr->has_show_existing = EB_FALSE;
        context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
    } else if (pcs_ptr->hierarchical_levels == 1) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(pcs_ptr, context_ptr);
            return;
        }

        const uint8_t  base0_idx = context_ptr->lay0_toggle == 0 ? 1 : context_ptr->lay0_toggle == 1 ? 2 : 0;
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 0 : 1;
        const uint8_t  base2_idx = context_ptr->lay0_toggle == 0 ? 0 : context_ptr->lay0_toggle == 1 ? 1 : 2;

        switch (pcs_ptr->temporal_layer_index) {
        case 0:

            //{2, 4, 0, 0},     // GOP Index 0 - Ref List 0
            //{2, 4, 0, 0},     // GOP Index 0 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = base0_idx;
            av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[LAST2];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
            gop_i = 0;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;
            break;

        case 1:
            //{ 1, 3, 0, 0}   // GOP Index 2 - Ref List 0
            //{-1, 0, 0, 0}   // GOP Index 2 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = base0_idx;
            av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base2_idx;
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];

            gop_i = 1;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 0;

            break;

        default:
            SVT_LOG("Error: unexpected picture mini Gop number\n");
            break;
        }

        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {
            if (pcs_ptr->is_used_as_reference_flag && picture_index != 0) {
                frm_hdr->show_frame = EB_FALSE;
                pcs_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                pcs_ptr->has_show_existing = EB_TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = base2_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", pcs_ptr->hierarchical_levels);
            }
        }

        // Toggle layer0 and layer1
        if ((picture_index == context_ptr->mini_gop_end_index[mini_gop_index] &&
                    !pcs_ptr->is_overlay &&
                    scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) ||
                ((scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {
            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
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

            //{4, 8, 0, 0},     // GOP Index 0 - Ref List 0
            //{4, 8, 0, 0},     // GOP Index 0 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = base0_idx;
            av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[LAST2];
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
            SVT_LOG("Error: unexpected picture mini Gop number\n");
            break;
        }

        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {
            if (pcs_ptr->is_used_as_reference_flag && picture_index != 0) {
                frm_hdr->show_frame = EB_FALSE;
                pcs_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                pcs_ptr->has_show_existing = EB_TRUE;

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
                    scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) ||
                ((scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                 scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {
            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
        }
    } else if (pcs_ptr->hierarchical_levels == 3) {

        uint8_t gop_i;
        EbBool is_trailing_frames = EB_FALSE;

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
        if (pcs_ptr->pred_struct_ptr->pred_type == EB_PRED_LOW_DELAY_P &&
                context_ptr->mini_gop_length[mini_gop_index] < 8 &&
                pcs_ptr->scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {
            is_trailing_frames = EB_TRUE;
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
            av1_rps->ref_dpb_index[LAST2] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[LAST3] = av1_rps->ref_dpb_index[LAST];
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base1_idx;
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
            gop_i = 0;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = av1_rps->ref_poc_array[LAST];
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
                     pcs_ptr->scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)) {
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
            if (picture_index == 0)
                av1_rps->refresh_frame_mask = 1 << (lay3_idx);
            else
                av1_rps->refresh_frame_mask = 0;
            break;

        default:
            SVT_LOG("Error: unexpected picture mini Gop number\n");
            break;
        }

        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {
            if (pcs_ptr->is_used_as_reference_flag && picture_index != 0)
            {
                frm_hdr->show_frame = EB_FALSE;
                pcs_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                pcs_ptr->has_show_existing = EB_TRUE;

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
                        scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)) ||
                ((scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {

            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;


            //layer2 toggle for low delay b/P case
            if ((scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                    scs_ptr->intra_period_length != -1 &&
                    pcs_ptr->cra_flag) {
                int32_t trailing_count = (scs_ptr->intra_period_length % 8);
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
            if (picture_index == 0 || picture_index == 8)
                av1_rps->refresh_frame_mask = 1 << (lay4_idx);
            else
                av1_rps->refresh_frame_mask = 0;
            break;

        default:
            SVT_LOG("Error: unexpected picture mini Gop number\n");
            break;
        }


        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {
            if (pcs_ptr->is_used_as_reference_flag && picture_index != 0 && picture_index != 8)
            {
                frm_hdr->show_frame = EB_FALSE;
                pcs_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                pcs_ptr->has_show_existing = EB_TRUE;

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
                        scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)) ||
                ((scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
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

        const uint8_t  base0_idx = context_ptr->lay0_toggle == 0 ? 1 : context_ptr->lay0_toggle == 1 ? 2 : 0;   //the oldest L0 picture in the DPB
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 0 : 1;   //the middle L0 picture in the DPB
        const uint8_t  base2_idx = context_ptr->lay0_toggle == 0 ? 0 : context_ptr->lay0_toggle == 1 ? 1 : 2;   //the newest L0 picture in the DPB

        const uint8_t  lay1_0_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 1 : LAY1_OFF + 0;   //the oldest L1 picture in the DPB
        const uint8_t  lay1_1_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 0 : LAY1_OFF + 1;   //the newest L1 picture in the DPB
        const uint8_t  lay2_idx = LAY2_OFF;   //the oldest L2 picture in the DPB
        const uint8_t  lay3_idx = LAY3_OFF;    //the newest L3 picture in the DPB
        const uint8_t  lay4_idx = LAY4_OFF;    //the newest L3 picture in the DPB

        switch (pcs_ptr->temporal_layer_index) {
        case 0:
            //{32, 64, 96, 0},      // GOP Index 0 - Ref List 0
            //{32, 0, 0, 0}       // GOP Index 0 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = base0_idx;
            av1_rps->ref_dpb_index[LAST3] = base2_idx;
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
            //{16, 32, 48, 0},   // GOP Index 16 - Ref List 0
            //{-16, 0, 0, 0}    // GOP Index 16 - Ref List 1
            av1_rps->ref_dpb_index[LAST] = base1_idx;
            av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
            av1_rps->ref_dpb_index[LAST3] = base0_idx;
            av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

            av1_rps->ref_dpb_index[BWD] = base2_idx;
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];

            gop_i = 16;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            av1_rps->refresh_frame_mask = 1 << (LAY1_OFF + context_ptr->lay1_toggle);

            break;

        case 2:

            if (picture_index == 7) {
                //{8, 16, 24, 40}, // GOP Index 8 - Ref List 0
                //{-8, -24, 0, 0} // GOP Index 8 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = base0_idx;

                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 8;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 23) {
                //{8, 16, 24, 0}    // GOP Index 24 - Ref List 0
                //{-8, 0, 0, 0}     // GOP Index 24 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 24;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            av1_rps->refresh_frame_mask = 1 << (LAY2_OFF);

            break;

        case 3:

            if (picture_index == 3) {
                //{4, 8, 20, 36}, // GOP Index 4 - Ref List 0
                //{-4, -12, -28, 0} // GOP Index 4 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                av1_rps->ref_dpb_index[GOLD] = base0_idx;
                av1_rps->ref_dpb_index[BWD] = lay2_idx;
                av1_rps->ref_dpb_index[ALT2] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT] = base2_idx;
                gop_i = 4;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
            }
            else if (picture_index == 11) {
                //{4, 8, 12, 28}, // GOP Index 12 - Ref List 0
                //{-4, -20, 0, 0} // GOP Index 12 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 12;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
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
                //{-4, 0, 0, 0} // GOP Index 28 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay2_idx;
                av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                av1_rps->ref_dpb_index[LAST3] = lay1_1_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;

                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 28;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
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
                //{-2, -18, 0, 0} // GOP Index 14 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = lay2_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                av1_rps->ref_dpb_index[ALT2] = base2_idx;
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 14;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
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
                //{-2, 0, 0, 0} // GOP Index 30 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay3_idx;
                av1_rps->ref_dpb_index[LAST2] = lay4_idx;
                av1_rps->ref_dpb_index[LAST3] = lay2_idx;
                av1_rps->ref_dpb_index[GOLD] = lay1_1_idx;
                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 30;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
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
                    //{1, 17, 33, 0}, // GOP Index 1 - Ref List 0
                    //{-1, -3, -7, 0} // GOP Index 1 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = base1_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_0_idx;
                    av1_rps->ref_dpb_index[LAST3] = base0_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT] = lay2_idx;
                    gop_i = 1;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 2) {
                    //{1, 3, 19, 0}, // GOP Index 3 - Ref List 0
                    //{-1, -5, -13, 0} // GOP Index 3 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                    av1_rps->ref_dpb_index[BWD] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay2_idx;
                    av1_rps->ref_dpb_index[ALT] = lay1_1_idx;
                    gop_i = 3;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 4) {
                    //{1, 5, 21, 0}, // GOP Index 5 - Ref List 0
                    //{-1, -3, -11, 0} // GOP Index 5 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay2_idx;
                    av1_rps->ref_dpb_index[ALT] = lay1_1_idx;
                    gop_i = 5;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 6) {
                    //{1, 3, 7, 23}, // GOP Index 7 - Ref List 0
                    //{-1, -9, -25, 0} // GOP Index 7 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST3] = base1_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
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
                    //{1, 9, 25, 0}, // GOP Index 9 - Ref List 0
                    //{-1, -3, -7, 0} // GOP Index 9 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT] = lay1_1_idx;
                    gop_i = 9;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 10) {
                    //{1, 3, 11, 27}, // GOP Index 11 - Ref List 0
                    //{-1, -5, -21, 0} // GOP Index 11 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST3] = base1_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
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
                    //{1, 5, 13, 29}, // GOP Index 13 - Ref List 0
                    //{-1, -3, -19, 0} // GOP Index 13 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST3] = base1_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
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
                    //{1, 3, 7, 15}, // GOP Index 15 - Ref List 0
                    //{-1, -17, 0, 0} // GOP Index 15 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay2_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
                    av1_rps->ref_dpb_index[BWD] = lay1_1_idx;
                    av1_rps->ref_dpb_index[ALT2] = base2_idx;
                    av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                    gop_i = 15;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
                }
                else if (picture_index == 16) {
                    //{1, 17, 33, 0}, // GOP Index 17 - Ref List 0
                    //{-1, -3, -7, 0} // GOP Index 17 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST2] = base1_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay1_0_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT] = lay2_idx;
                    gop_i = 17;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 18) {
                    //{1, 3, 19, 35}, // GOP Index 19 - Ref List 0
                    //{-1, -5, -13, 0} // GOP Index 19 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST3] = base1_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
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
                    //{1, 5, 21, 37}, // GOP Index 21 - Ref List 0
                    //{-1, -3, -11, 0} // GOP Index 21 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST3] = base1_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_0_idx;
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
                    //{1, 3, 7, 23}, // GOP Index 23 - Ref List 0
                    //{-1, -9, 0, 0} // GOP Index 23 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay1_1_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
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
                    //{1, 9, 25, 0}, // GOP Index 25 - Ref List 0
                    //{-1, -3, -7, 0} // GOP Index 25 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay1_1_idx;
                    av1_rps->ref_dpb_index[LAST3] = base1_idx;
                    av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];
                    av1_rps->ref_dpb_index[BWD] = lay4_idx;
                    av1_rps->ref_dpb_index[ALT2] = lay3_idx;
                    av1_rps->ref_dpb_index[ALT] = base2_idx;
                    gop_i = 25;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                    av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
                }
                else if (picture_index == 26) {
                    //{1, 3, 11, 27}, // GOP Index 27 - Ref List 0
                    //{-1, -5, 0, 0} // GOP Index 27 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay1_1_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
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
                    //{1, 5, 13, 29}, // GOP Index 29 - Ref List 0
                    //{-1, -3, 0, 0} // GOP Index 29 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay2_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay1_1_idx;
                    av1_rps->ref_dpb_index[GOLD] = base1_idx;
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
                    //{1, 3, 7, 15}, // GOP Index 31 - Ref List 0
                    //{-1, 0, 0, 0} // GOP Index 31 - Ref List 1
                    av1_rps->ref_dpb_index[LAST] = lay4_idx;
                    av1_rps->ref_dpb_index[LAST2] = lay3_idx;
                    av1_rps->ref_dpb_index[LAST3] = lay2_idx;
                    av1_rps->ref_dpb_index[GOLD] = lay1_1_idx;
                    av1_rps->ref_dpb_index[BWD] = base2_idx;
                    av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                    av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                    gop_i = 31;
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, pcs_ptr->picture_number, six_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
                    av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
                }
                else
                    SVT_LOG("Error in GOp indexing\n");
            }

            av1_rps->refresh_frame_mask = 0;
            break;

        default:
            SVT_LOG("Error: unexpected picture mini Gop number\n");
            break;
        }


        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(pcs_ptr, context_ptr, mini_gop_index)) {
            if (pcs_ptr->is_used_as_reference_flag && (picture_index % 8 != 0)) {
                frm_hdr->show_frame = EB_FALSE;
                pcs_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                pcs_ptr->has_show_existing = EB_TRUE;

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
                        scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)) ||
                ((scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  scs_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                 (pcs_ptr->temporal_layer_index == 0))) {

            //Layer0 toggle 0->1
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 2->3
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
        }
    }else {
        SVT_LOG("Error: Not supported GOP structure!");
        exit(0);
    }
}

static EbErrorType av1_generate_rps_info_from_user_config(
    PictureParentControlSet       *picture_control_set_ptr,
    EncodeContext                 *encode_context_ptr
)
{
    Av1RpsNode *av1_rps = &picture_control_set_ptr->av1_ref_signal;
    FrameHeader *frm_hdr = &picture_control_set_ptr->frm_hdr;
    PredictionStructureEntry *pred_position_ptr = picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[picture_control_set_ptr->pred_struct_index];
    DPBInfo *dpb_list_ptr = &encode_context_ptr->dpb_list[0];
    int32_t dpb_list_idx = 0;
    uint64_t ref_poc = 0;
    uint32_t dep_idx = 0;
    uint8_t  ref_idx = 0;

    if (frm_hdr->frame_type == KEY_FRAME) {
        frm_hdr->show_frame = EB_TRUE;
        picture_control_set_ptr->has_show_existing = EB_FALSE;
        EB_MEMSET(dpb_list_ptr, 0, sizeof(DPBInfo)*REF_FRAMES);
        encode_context_ptr->display_picture_number = picture_control_set_ptr->picture_number;
        dpb_list_ptr[dpb_list_idx].is_used = EB_TRUE;
        dpb_list_ptr[dpb_list_idx].picture_number = picture_control_set_ptr->picture_number;
        dpb_list_ptr[dpb_list_idx].is_displayed = EB_TRUE;
        dpb_list_ptr[dpb_list_idx].temporal_layer_index = picture_control_set_ptr->temporal_layer_index;

        // Construct dependent lists
        dpb_list_ptr[dpb_list_idx].dep_list0.list_count = 0;
        for (dep_idx = 0; dep_idx < pred_position_ptr->dep_list0.list_count; ++dep_idx) {
            if (pred_position_ptr->dep_list0.list[dep_idx] >= 0){
                dpb_list_ptr[dpb_list_idx].dep_list0.list[dpb_list_ptr[dpb_list_idx].dep_list0.list_count++] = pred_position_ptr->dep_list0.list[dep_idx];
            }
        }
        dpb_list_ptr[dpb_list_idx].dep_list1.list_count = pred_position_ptr->dep_list1.list_count;
        for (dep_idx = 0; dep_idx < pred_position_ptr->dep_list1.list_count; ++dep_idx) {
            dpb_list_ptr[dpb_list_idx].dep_list1.list[dep_idx] = pred_position_ptr->dep_list1.list[dep_idx];
        }
        dpb_list_ptr[dpb_list_idx].dep_list0_count = dpb_list_ptr[dpb_list_idx].dep_list0.list_count;
        dpb_list_ptr[dpb_list_idx].dep_list1_count = dpb_list_ptr[dpb_list_idx].dep_list1.list_count;
        dpb_list_ptr[dpb_list_idx].dep_count = dpb_list_ptr[dpb_list_idx].dep_list0_count + dpb_list_ptr[dpb_list_idx].dep_list1_count;
        return EB_ErrorNone;
    }

    // If there was an I-frame or Scene Change, then cleanup the Reference Queue's Dependent Counts
    if (picture_control_set_ptr->slice_type == I_SLICE)
    {
        encode_context_ptr->is_i_slice_in_last_mini_gop = EB_TRUE;
        encode_context_ptr->i_slice_picture_number_in_last_mini_gop = picture_control_set_ptr->picture_number;
    }

    // Construct dpb index mapping for ref list0
    if ((picture_control_set_ptr->slice_type == P_SLICE) || (picture_control_set_ptr->slice_type == B_SLICE)) {
        for (ref_idx = LAST; ref_idx < LAST + picture_control_set_ptr->ref_list0_count; ++ref_idx) {
            if (picture_control_set_ptr->is_overlay) {
                ref_poc = picture_control_set_ptr->picture_number;
            }
            else {
                ref_poc = picture_control_set_ptr->picture_number - pred_position_ptr->ref_list0.reference_list[ref_idx-LAST];
            }
            dpb_list_idx = 0;
            do {
                if (dpb_list_ptr[dpb_list_idx].is_used == EB_TRUE && dpb_list_ptr[dpb_list_idx].picture_number == ref_poc) {
                    if (dpb_list_ptr[dpb_list_idx].temporal_layer_index > picture_control_set_ptr->temporal_layer_index) {
                        return EB_Corrupt_Frame;
                    }
                    else {
                        break;
                    }
                }
            } while (++dpb_list_idx < REF_FRAMES);
            if (dpb_list_idx < REF_FRAMES) {
                av1_rps->ref_dpb_index[ref_idx] = dpb_list_idx;
                --dpb_list_ptr[dpb_list_idx].dep_count;
                if(dpb_list_ptr[dpb_list_idx].dep_count < 0){
                    SVT_LOG("Error: dep_count error in dpb list0\n");
                    return EB_Corrupt_Frame;
                }
            }
            else {
                SVT_LOG("Error: can't find ref frame in dpb list0\n");
                return EB_Corrupt_Frame;
            }
        }
        for (; ref_idx <= GOLD; ++ref_idx) {
            av1_rps->ref_dpb_index[ref_idx] = av1_rps->ref_dpb_index[LAST];
        }

        // Construct dpb index mapping for ref list1
        if (picture_control_set_ptr->slice_type == B_SLICE) {
            for (ref_idx = BWD; ref_idx < BWD + picture_control_set_ptr->ref_list1_count; ++ref_idx) {
                ref_poc = picture_control_set_ptr->picture_number - pred_position_ptr->ref_list1.reference_list[ref_idx-BWD];
                dpb_list_idx = 0;
                do {
                    if (dpb_list_ptr[dpb_list_idx].is_used == EB_TRUE && dpb_list_ptr[dpb_list_idx].picture_number == ref_poc) {
                        if (dpb_list_ptr[dpb_list_idx].temporal_layer_index > picture_control_set_ptr->temporal_layer_index) {
                            return EB_Corrupt_Frame;
                        }
                        else {
                            break;
                        }
                    }
                } while (++dpb_list_idx < REF_FRAMES);

                if (dpb_list_idx < REF_FRAMES) {
                    av1_rps->ref_dpb_index[ref_idx] = dpb_list_idx;
                    --dpb_list_ptr[dpb_list_idx].dep_count;
                    if(dpb_list_ptr[dpb_list_idx].dep_count < 0){
                        SVT_LOG("Error: dep_count error in dpb list1\n");
                        return EB_Corrupt_Frame;
                    }
                }
                else {
                    SVT_LOG("Error: can't find ref frame in dpb list1\n");
                    return EB_Corrupt_Frame;
                }
            }
            for (; ref_idx <= ALT; ++ref_idx) {
                av1_rps->ref_dpb_index[ref_idx] = av1_rps->ref_dpb_index[BWD];
            }
        }
        else{
            av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
        }
    }

    // Release unused positions in dpb list
    dpb_list_idx = 0;
    do {
        if (dpb_list_ptr[dpb_list_idx].is_used == EB_TRUE
        && dpb_list_ptr[dpb_list_idx].is_displayed == EB_TRUE
        && dpb_list_ptr[dpb_list_idx].dep_count == 0) {
            dpb_list_ptr[dpb_list_idx].is_used = EB_FALSE;
            dpb_list_ptr[dpb_list_idx].is_displayed = EB_FALSE;
        }
    } while (++dpb_list_idx < REF_FRAMES);

    // Insert current frame into dpb list
    if (picture_control_set_ptr->picture_number == encode_context_ptr->display_picture_number + 1
    && picture_control_set_ptr->is_used_as_reference_flag == 0 ) {
        frm_hdr->show_frame = EB_TRUE;
        ++encode_context_ptr->display_picture_number;
        av1_rps->refresh_frame_mask = 0;
    }
    else {
        frm_hdr->show_frame = EB_FALSE;
        dpb_list_idx = 0;
        do {
            if (dpb_list_ptr[dpb_list_idx].is_used == EB_FALSE) {
                break;
            }
        } while (++dpb_list_idx < REF_FRAMES);
        if (dpb_list_idx < REF_FRAMES) {
            av1_rps->refresh_frame_mask = 1 << dpb_list_idx;
            dpb_list_ptr[dpb_list_idx].is_used = EB_TRUE;
            dpb_list_ptr[dpb_list_idx].picture_number = picture_control_set_ptr->picture_number;
            dpb_list_ptr[dpb_list_idx].is_alt_ref = picture_control_set_ptr->is_alt_ref;
            dpb_list_ptr[dpb_list_idx].temporal_layer_index = picture_control_set_ptr->temporal_layer_index;
            if (picture_control_set_ptr->is_alt_ref) {
                dpb_list_ptr[dpb_list_idx].is_displayed = EB_TRUE;
            }
            // Construct dependent lists
            dpb_list_ptr[dpb_list_idx].dep_list0.list_count = 0;
            for (dep_idx = 0; dep_idx < pred_position_ptr->dep_list0.list_count; ++dep_idx) {
                if (pred_position_ptr->dep_list0.list[dep_idx] >= 0){
                    dpb_list_ptr[dpb_list_idx].dep_list0.list[dpb_list_ptr[dpb_list_idx].dep_list0.list_count++] = pred_position_ptr->dep_list0.list[dep_idx];
                }
            }
            dpb_list_ptr[dpb_list_idx].dep_list1.list_count = pred_position_ptr->dep_list1.list_count;
            for (dep_idx = 0; dep_idx < pred_position_ptr->dep_list1.list_count; ++dep_idx) {
                dpb_list_ptr[dpb_list_idx].dep_list1.list[dep_idx] = pred_position_ptr->dep_list1.list[dep_idx];
            }
            dpb_list_ptr[dpb_list_idx].dep_list0_count = (picture_control_set_ptr->is_alt_ref) ? dpb_list_ptr[dpb_list_idx].dep_list0.list_count + 1 : dpb_list_ptr[dpb_list_idx].dep_list0.list_count;
            dpb_list_ptr[dpb_list_idx].dep_list1_count = dpb_list_ptr[dpb_list_idx].dep_list1.list_count;
            dpb_list_ptr[dpb_list_idx].dep_count = dpb_list_ptr[dpb_list_idx].dep_list0_count + dpb_list_ptr[dpb_list_idx].dep_list1_count;
        }
        else {
            SVT_LOG("Error: can't find unused dpb to hold current frame\n");
            return EB_Corrupt_Frame;
        }
    }

    // Calculate output flag
    picture_control_set_ptr->has_show_existing = EB_FALSE;
    do {
        dpb_list_idx = 0;
        do {
            if (dpb_list_ptr[dpb_list_idx].is_used == EB_TRUE
            && dpb_list_ptr[dpb_list_idx].is_displayed == EB_FALSE
            && dpb_list_ptr[dpb_list_idx].picture_number == encode_context_ptr->display_picture_number + 1) {
                break;
            }
        } while (++dpb_list_idx < REF_FRAMES);
        if (dpb_list_idx < REF_FRAMES) {
            dpb_list_ptr[dpb_list_idx].is_displayed = EB_TRUE;
            ++encode_context_ptr->display_picture_number;
            if (encode_context_ptr->display_picture_number == picture_control_set_ptr->picture_number) {
                frm_hdr->show_frame = EB_TRUE;
            }
            else {
                picture_control_set_ptr->has_show_existing = EB_TRUE;
                frm_hdr->show_existing_frame = dpb_list_idx;
                break;
            }
        }
        else {
            break;
        }
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

    if (encode_context_ptr->is_i_slice_in_last_mini_gop == EB_TRUE)
    {
        encode_context_ptr->is_i_slice_in_last_mini_gop = EB_FALSE;
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
    SequenceControlSet* scs_ptr = (SequenceControlSet*)picture_control_set_ptr->scs_wrapper_ptr->object_ptr;
    for (uint32_t decode_order = 0,pictureIndex = context_ptr->mini_gop_start_index[mini_gop_index]; pictureIndex <= context_ptr->mini_gop_end_index[mini_gop_index]+has_overlay; ++pictureIndex, ++decode_order) {
        if (has_overlay && pictureIndex == context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay) {
            picture_control_set_ptr = ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_end_index[mini_gop_index]]->object_ptr)->overlay_ppcs_ptr;
        }
        else if ((context_ptr->mini_gop_length[mini_gop_index] == picture_control_set_ptr->pred_struct_ptr->pred_struct_period)
        &&(context_ptr->mini_gop_idr_count[mini_gop_index] == 0)
        &&(scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)){
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

static void av1_generate_rps_ref_poc_from_user_config(PictureParentControlSet *picture_control_set_ptr)
{
    Av1RpsNode *av1_rps = &picture_control_set_ptr->av1_ref_signal;
    FrameHeader *frm_hdr = &picture_control_set_ptr->frm_hdr;
    PredictionStructureEntry *pred_position_ptr = picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[picture_control_set_ptr->pred_struct_index];
    uint32_t ref_idx = 0;

    if (picture_control_set_ptr->slice_type == I_SLICE)
        frm_hdr->frame_type = picture_control_set_ptr->idr_flag ? KEY_FRAME : INTRA_ONLY_FRAME;
    else
        frm_hdr->frame_type = INTER_FRAME;

    picture_control_set_ptr->intra_only = picture_control_set_ptr->slice_type == I_SLICE ? 1 : 0;

    if (picture_control_set_ptr->is_overlay) {
        for (ref_idx = LAST; ref_idx <= ALT; ++ref_idx) {
            av1_rps->ref_poc_array[ref_idx] = picture_control_set_ptr->picture_number;
        }
    }
    else {
        if ((picture_control_set_ptr->slice_type == P_SLICE) || (picture_control_set_ptr->slice_type == B_SLICE)) {
            for (ref_idx = LAST; ref_idx < LAST + pred_position_ptr->ref_list0.reference_list_count; ++ref_idx) {
                av1_rps->ref_poc_array[ref_idx] = picture_control_set_ptr->picture_number - pred_position_ptr->ref_list0.reference_list[ref_idx-LAST];
            }
            for (; ref_idx <= GOLD; ++ref_idx) {
                av1_rps->ref_poc_array[ref_idx] = av1_rps->ref_poc_array[LAST];
            }
            if (picture_control_set_ptr->slice_type == B_SLICE) {
                for (ref_idx = BWD; ref_idx < BWD + pred_position_ptr->ref_list1.reference_list_count; ++ref_idx) {
                  av1_rps->ref_poc_array[ref_idx] = picture_control_set_ptr->picture_number - pred_position_ptr->ref_list1.reference_list[ref_idx-BWD];
                }
                for (; ref_idx <= ALT; ++ref_idx) {
                  av1_rps->ref_poc_array[ref_idx] = av1_rps->ref_poc_array[BWD];
                }
            }
            else {
                av1_rps->ref_poc_array[BWD] = av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[LAST];
            }
        }
    }
    return;
}

/***************************************************************************************************
// Perform Required Picture Analysis Processing for the Overlay frame
***************************************************************************************************/
void perform_simple_picture_analysis_for_overlay(PictureParentControlSet     *pcs_ptr) {
    EbPictureBufferDesc           *input_padded_picture_ptr;
    EbPictureBufferDesc           *input_picture_ptr;
    EbPaReferenceObject           *pa_ref_obj_;
    uint32_t                        pic_width_in_sb;
    uint32_t                        pic_height_in_sb;
    uint32_t                        sb_total_count;

    SequenceControlSet *scs_ptr = (SequenceControlSet*)pcs_ptr->scs_wrapper_ptr->object_ptr;
    input_picture_ptr               = pcs_ptr->enhanced_picture_ptr;
    pa_ref_obj_               = (EbPaReferenceObject*)pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
    input_padded_picture_ptr        = (EbPictureBufferDesc*)pa_ref_obj_->input_padded_picture_ptr;
    pic_width_in_sb = (pcs_ptr->aligned_width + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
    pic_height_in_sb   = (pcs_ptr->aligned_height + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
    sb_total_count      = pic_width_in_sb * pic_height_in_sb;

    // Pad pictures to multiple min cu size
    pad_picture_to_multiple_of_min_blk_size_dimensions(
        scs_ptr,
        input_picture_ptr);

    // Pre processing operations performed on the input picture
    picture_pre_processing_operations(
        pcs_ptr,
        scs_ptr,
        sb_total_count);

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
            eb_memcpy(pa + row * input_padded_picture_ptr->stride_y, in + row * input_picture_ptr->stride_y, sizeof(uint8_t) * input_picture_ptr->width);
    }

    // Pad input picture to complete border SBs
    pad_picture_to_multiple_of_sb_dimensions(
        input_padded_picture_ptr);

    // 1/4 & 1/16 input picture decimation
    downsample_decimation_input_picture(
        pcs_ptr,
        input_padded_picture_ptr,
        (EbPictureBufferDesc*)pa_ref_obj_->quarter_decimated_picture_ptr,
        (EbPictureBufferDesc*)pa_ref_obj_->sixteenth_decimated_picture_ptr);

    // 1/4 & 1/16 input picture downsampling through filtering
    if (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
        downsample_filtering_input_picture(
            pcs_ptr,
            input_padded_picture_ptr,
            (EbPictureBufferDesc*)pa_ref_obj_->quarter_filtered_picture_ptr,
            (EbPictureBufferDesc*)pa_ref_obj_->sixteenth_filtered_picture_ptr);
    }
    // Gathering statistics of input picture, including Variance Calculation, Histogram Bins
    gathering_picture_statistics(
        scs_ptr,
        pcs_ptr,
        pcs_ptr->chroma_downsampled_picture_ptr, //420 input_picture_ptr
        input_padded_picture_ptr,
        pa_ref_obj_->sixteenth_decimated_picture_ptr, // Hsan: always use decimated until studying the trade offs
        sb_total_count);

    pcs_ptr->sc_content_detected = pcs_ptr->alt_ref_ppcs_ptr->sc_content_detected;
}
/***************************************************************************************************
 * Initialize the overlay frame
***************************************************************************************************/
void initialize_overlay_frame(PictureParentControlSet     *pcs_ptr) {
    pcs_ptr->fade_out_from_black = pcs_ptr->alt_ref_ppcs_ptr->fade_out_from_black;
    pcs_ptr->fade_in_to_black = pcs_ptr->alt_ref_ppcs_ptr->fade_in_to_black;
    pcs_ptr->scene_change_flag = EB_FALSE;
    pcs_ptr->cra_flag = EB_FALSE;
    pcs_ptr->idr_flag = EB_FALSE;
    pcs_ptr->target_bit_rate = pcs_ptr->alt_ref_ppcs_ptr->target_bit_rate;
    pcs_ptr->last_idr_picture = pcs_ptr->alt_ref_ppcs_ptr->last_idr_picture;
    pcs_ptr->pred_structure = pcs_ptr->alt_ref_ppcs_ptr->pred_structure;
    pcs_ptr->pred_struct_ptr = pcs_ptr->alt_ref_ppcs_ptr->pred_struct_ptr;
    pcs_ptr->hierarchical_levels = pcs_ptr->alt_ref_ppcs_ptr->hierarchical_levels;
    pcs_ptr->hierarchical_layers_diff = 0;
    pcs_ptr->init_pred_struct_position_flag = EB_FALSE;
    pcs_ptr->pre_assignment_buffer_count = pcs_ptr->alt_ref_ppcs_ptr->pre_assignment_buffer_count;

    perform_simple_picture_analysis_for_overlay(pcs_ptr);
 }

/***************************************************************************************************
 * Helper function. Compare two frames: center frame and target frame. Return the summation of
 * absolute difference between the two frames from a histogram of luma values
***************************************************************************************************/

static __inline uint32_t compute_luma_sad_between_center_and_target_frame(
    int center_index,
    int target_frame_index,
    PictureParentControlSet *pcs_ptr,
    SequenceControlSet *scs_ptr) {

    int32_t center_sum = 0, altref_sum = 0;
    uint32_t ahd = 0;

    for (int bin = 0; bin < HISTOGRAM_NUMBER_OF_BINS; ++bin) {
        center_sum = 0, altref_sum = 0;
        for (uint32_t region_in_picture_width_index = 0; region_in_picture_width_index < scs_ptr->picture_analysis_number_of_regions_per_width; region_in_picture_width_index++) {
            for (uint32_t region_in_picture_height_index = 0; region_in_picture_height_index < scs_ptr->picture_analysis_number_of_regions_per_height; region_in_picture_height_index++) {
                center_sum += pcs_ptr->temp_filt_pcs_list[center_index]->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][0][bin];
                altref_sum += pcs_ptr->temp_filt_pcs_list[target_frame_index]->picture_histogram[region_in_picture_width_index][region_in_picture_height_index][0][bin];
            }
        }
        ahd += ABS(center_sum - altref_sum);
    }
    return ahd;
}

#if NOISE_BASED_TF_FRAMES
double estimate_noise(const uint8_t *src, uint16_t width, uint16_t height,
    uint16_t stride_y);

double estimate_noise_highbd(const uint16_t *src, int width, int height, int stride,
    int bd);

void pack_highbd_pic(const EbPictureBufferDesc *pic_ptr, uint16_t *buffer_16bit[3], uint32_t ss_x,
    uint32_t ss_y, EbBool include_padding);

EbErrorType derive_tf_window_params(
#else
// Derive past_altref_nframes and future_altref_nframes using the central_frame activity compared to default past_altref_nframes and future_altref_nframes
void derive_tf_window_params(
#endif
    SequenceControlSet *scs_ptr,
    EncodeContext *encode_context_ptr,
    PictureParentControlSet *pcs_ptr,
#if TF_LEVELS
    PictureDecisionContext *context_ptr,
#endif
    uint32_t out_stride_diff64) {
#if NOISE_BASED_TF_FRAMES
    PictureParentControlSet * picture_control_set_ptr_central = pcs_ptr;
    EbPictureBufferDesc * central_picture_ptr = picture_control_set_ptr_central->enhanced_picture_ptr;
    uint32_t encoder_bit_depth = picture_control_set_ptr_central->scs_ptr->static_config.encoder_bit_depth;
    EbBool is_highbd = (encoder_bit_depth == 8) ? (uint8_t)EB_FALSE : (uint8_t)EB_TRUE;

    // chroma subsampling
    uint32_t ss_x = picture_control_set_ptr_central->scs_ptr->subsampling_x;
    uint32_t ss_y = picture_control_set_ptr_central->scs_ptr->subsampling_y;
    double *noise_levels = &(picture_control_set_ptr_central->noise_levels[0]);

    // allocate 16 bit buffer
    if (is_highbd) {
        EB_MALLOC_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_Y],
            central_picture_ptr->luma_size);
        EB_MALLOC_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_U],
            central_picture_ptr->chroma_size);
        EB_MALLOC_ARRAY(picture_control_set_ptr_central->altref_buffer_highbd[C_V],
            central_picture_ptr->chroma_size);

        // pack byte buffers to 16 bit buffer
        pack_highbd_pic(central_picture_ptr,
            picture_control_set_ptr_central->altref_buffer_highbd,
            ss_x,
            ss_y,
            EB_TRUE);
    }

    // Estimate source noise level
    if (is_highbd) {
        uint16_t *altref_buffer_highbd_start[COLOR_CHANNELS];
        altref_buffer_highbd_start[C_Y] =
            picture_control_set_ptr_central->altref_buffer_highbd[C_Y] +
            central_picture_ptr->origin_y * central_picture_ptr->stride_y +
            central_picture_ptr->origin_x;

        altref_buffer_highbd_start[C_U] =
            picture_control_set_ptr_central->altref_buffer_highbd[C_U] +
            (central_picture_ptr->origin_y >> ss_y) * central_picture_ptr->stride_bit_inc_cb +
            (central_picture_ptr->origin_x >> ss_x);

        altref_buffer_highbd_start[C_V] =
            picture_control_set_ptr_central->altref_buffer_highbd[C_V] +
            (central_picture_ptr->origin_y >> ss_y) * central_picture_ptr->stride_bit_inc_cr +
            (central_picture_ptr->origin_x >> ss_x);

        noise_levels[0] = estimate_noise_highbd(altref_buffer_highbd_start[C_Y], // Y only
            central_picture_ptr->width,
            central_picture_ptr->height,
            central_picture_ptr->stride_y,
            encoder_bit_depth);

        noise_levels[1] = estimate_noise_highbd(altref_buffer_highbd_start[C_U], // U only
            (central_picture_ptr->width >> 1),
            (central_picture_ptr->height >> 1),
            central_picture_ptr->stride_cb,
            encoder_bit_depth);

        noise_levels[2] = estimate_noise_highbd(altref_buffer_highbd_start[C_V], // V only
            (central_picture_ptr->width >> 1),
            (central_picture_ptr->height >> 1),
            central_picture_ptr->stride_cb,
            encoder_bit_depth);

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

        noise_levels[0] = estimate_noise(buffer_y, // Y
            central_picture_ptr->width,
            central_picture_ptr->height,
            central_picture_ptr->stride_y);

        noise_levels[1] = estimate_noise(buffer_u, // U
            (central_picture_ptr->width >> ss_x),
            (central_picture_ptr->height >> ss_y),
            central_picture_ptr->stride_cb);

        noise_levels[2] = estimate_noise(buffer_v, // V
            (central_picture_ptr->width >> ss_x),
            (central_picture_ptr->height >> ss_y),
            central_picture_ptr->stride_cr);
    }

    // Adjust number of filtering frames based on noise and quantization factor.
    // Basically, we would like to use more frames to filter low-noise frame such
    // that the filtered frame can provide better predictions for more frames.
    // Also, when the quantization factor is small enough (lossless compression),
    // we will not change the number of frames for key frame filtering, which is
    // to avoid visual quality drop.
    int adjust_num = 0;
#if IMPROVED_TF_LEVELS
    if (context_ptr->tf_ctrls.noise_based_window_adjust) {
#endif
#if 0
    if (num_frames == 1) {  // `arnr_max_frames = 1` is used to disable filtering.
        adjust_num = 0;
    }
    else if (filter_frame_lookahead_idx < 0 && q <= 10) {
        adjust_num = 0;
    }
    else
#endif
    if (noise_levels[0] < 0.5) {
        adjust_num = 6;
    }
    else if (noise_levels[0] < 1.0) {
        adjust_num = 4;
    }
    else if (noise_levels[0] < 2.0) {
        adjust_num = 2;
    }
#endif
#if IMPROVED_TF_LEVELS
    }
#endif
#if TF_LEVELS
#if NOISE_BASED_TF_FRAMES
    int altref_nframes = MIN(scs_ptr->static_config.altref_nframes, context_ptr->tf_ctrls.window_size + adjust_num);
#else
    int altref_nframes = MIN(scs_ptr->static_config.altref_nframes, context_ptr->tf_ctrls.window_size);
#endif
#else
    int altref_nframes = scs_ptr->static_config.altref_nframes;
#endif
    if (pcs_ptr->idr_flag) {

        //initilize list
        for (int pic_itr = 0; pic_itr < ALTREF_MAX_NFRAMES; pic_itr++)
            pcs_ptr->temp_filt_pcs_list[pic_itr] = NULL;

        pcs_ptr->temp_filt_pcs_list[0] = pcs_ptr;
#if TF_LEVELS
        uint32_t num_future_pics = altref_nframes - 1;
#else
        uint32_t num_future_pics = 6;
#endif
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
        int index_center = 0;
        uint32_t actual_future_pics = pcs_ptr->future_altref_nframes;
        int pic_itr, ahd;

        int ahd_th = (((pcs_ptr->aligned_width * pcs_ptr->aligned_height) * AHD_TH_WEIGHT) / 100);

        // Accumulative histogram absolute differences between the central and future frame
        for (pic_itr = (index_center + actual_future_pics); pic_itr > index_center; pic_itr--) {
            ahd = compute_luma_sad_between_center_and_target_frame(index_center, pic_itr, pcs_ptr, scs_ptr);
            if (ahd < ahd_th)
                break;
        }
        pcs_ptr->future_altref_nframes = pic_itr - index_center;
        //SVT_LOG("\nPOC %d\t PAST %d\t FUTURE %d\n", pcs_ptr->picture_number, pcs_ptr->past_altref_nframes, pcs_ptr->future_altref_nframes);
    }
    else
    {
        int num_past_pics = altref_nframes / 2;
        int num_future_pics = altref_nframes - num_past_pics - 1;
        assert(altref_nframes <= ALTREF_MAX_NFRAMES);

        //initilize list
        for (int pic_itr = 0; pic_itr < ALTREF_MAX_NFRAMES; pic_itr++)
            pcs_ptr->temp_filt_pcs_list[pic_itr] = NULL;
        // limit the number of pictures to make sure there are enough pictures in the buffer. i.e. Intra CRA case
        num_past_pics = MIN(MIN(num_past_pics, (int)encode_context_ptr->pre_assignment_buffer_count - 1), (int)out_stride_diff64);

        //get previous+current pictures from the the pre-assign buffer
        for (int pic_itr = 0; pic_itr <= num_past_pics; pic_itr++) {
            PictureParentControlSet* pcs_itr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[out_stride_diff64 - num_past_pics + pic_itr]->object_ptr;
            pcs_ptr->temp_filt_pcs_list[pic_itr] = pcs_itr;
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
        int index_center = actual_past_pics;
        int pic_itr;
        int ahd;

        int ahd_th = (((pcs_ptr->aligned_width * pcs_ptr->aligned_height) * AHD_TH_WEIGHT) / 100);

        // Accumulative histogram absolute differences between the central and past frame
        for (pic_itr = index_center - actual_past_pics; pic_itr < index_center; pic_itr++) {
            ahd = compute_luma_sad_between_center_and_target_frame(index_center, pic_itr, pcs_ptr, scs_ptr);

            if (ahd < ahd_th)
                break;
        }
        pcs_ptr->past_altref_nframes = actual_past_pics = index_center - pic_itr;

        if (pcs_ptr->temporal_layer_index == 1)
            pcs_ptr->past_altref_nframes = actual_past_pics = MIN(1, pcs_ptr->past_altref_nframes);

        // Accumulative histogram absolute differences between the central and past frame
        for (pic_itr = (index_center + actual_future_pics); pic_itr > index_center; pic_itr--) {
            ahd = compute_luma_sad_between_center_and_target_frame(index_center, pic_itr, pcs_ptr, scs_ptr);
            if (ahd < ahd_th)
                break;
        }
        pcs_ptr->future_altref_nframes = pic_itr - index_center;

        if (pcs_ptr->temporal_layer_index == 1)
            pcs_ptr->future_altref_nframes = MIN(1, pcs_ptr->future_altref_nframes);

        //SVT_LOG("\nPOC %d\t PAST %d\t FUTURE %d\n", pcs_ptr->picture_number, pcs_ptr->past_altref_nframes, pcs_ptr->future_altref_nframes);

        // adjust the temporal filtering pcs buffer to remove unused past pictures
        if (actual_past_pics != num_past_pics) {

            pic_i = 0;
            while (pcs_ptr->temp_filt_pcs_list[pic_i] != NULL) {
                pcs_ptr->temp_filt_pcs_list[pic_i] = pcs_ptr->temp_filt_pcs_list[pic_i + num_past_pics - actual_past_pics];
                pic_i++;
            }
        }
    }
#if NOISE_BASED_TF_FRAMES
    return EB_ErrorNone;
#endif
}

#if DECOUPLE_ME_RES
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
#endif
#if ON_OFF_FEATURE_MRP

/***************************************************************************************************
// set number of references to try based on mrp level
***************************************************************************************************/
void set_mrp_controls(PictureParentControlSet *pcs_ptr) {

    MrpControls *mrp_ctrls = &pcs_ptr->mrp_ctrls;

    switch (pcs_ptr->mrp_level)
    {
    case 0:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 1);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 1);
        break;
    case 1:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 4);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
        break;
    case 2:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 4);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
        break;
    case 3:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 3);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
        break;
    case 4:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 3);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 2);
        break;
    case 5:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 2);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
        break;
    case 6:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 2);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 2);
        break;
    case 7:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 2);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 1);
        break;
    case 8:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 1);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 2);
        break;
    case 9:
        mrp_ctrls->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 1);
        mrp_ctrls->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 1);
        break;
    default:
        assert(0);
        break;
    }
    assert(mrp_ctrls->ref_list0_count_try <= pcs_ptr->ref_list0_count);
    assert(mrp_ctrls->ref_list1_count_try <= pcs_ptr->ref_list1_count);
}
#endif
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
    FrameHeader                   *frm_hdr;

    EncodeContext                 *encode_context_ptr;
    SequenceControlSet            *scs_ptr;
#if DECOUPLE_ME_RES
    EbObjectWrapper               *me_wrapper_ptr;
#endif
    EbObjectWrapper               *in_results_wrapper_ptr;
    PictureAnalysisResults        *in_results_ptr;

    EbObjectWrapper               *out_results_wrapper_ptr;
    PictureDecisionResults        *out_results_ptr;

    PredictionStructureEntry      *pred_position_ptr;

    EbBool                          pre_assignment_buffer_first_pass_flag;
    EB_SLICE                         picture_type;

    PictureDecisionReorderEntry   *queue_entry_ptr;
    int32_t                           queue_entry_index;

    int32_t                           previous_entry_index;

    PaReferenceQueueEntry         *input_entry_ptr = (PaReferenceQueueEntry*)NULL;;
    uint32_t                           input_queue_index;

    PaReferenceQueueEntry         *pa_reference_entry_ptr;
#if !DECOUPLE_ME_RES
    uint32_t                           pa_reference_queue_index;
#endif
    uint64_t                           ref_poc;

    uint32_t                           dep_idx;
    uint64_t                           dep_poc;

    uint32_t                           dep_list_count;

    // Dynamic GOP
    uint32_t                           mini_gop_index;
    uint32_t                           out_stride_diff64;

    EbBool                          window_avail, frame_passthrough;
    uint32_t                           window_index;
    uint32_t                           entry_index;
    PictureParentControlSet        *parent_pcs_window[FUTURE_WINDOW_WIDTH + 2];

    // Debug
    uint64_t                           loop_count = 0;

    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(
            context_ptr->picture_analysis_results_input_fifo_ptr,
            &in_results_wrapper_ptr);

        in_results_ptr = (PictureAnalysisResults*)in_results_wrapper_ptr->object_ptr;
        pcs_ptr = (PictureParentControlSet*)in_results_ptr->pcs_wrapper_ptr->object_ptr;
        scs_ptr = (SequenceControlSet*)pcs_ptr->scs_wrapper_ptr->object_ptr;
        frm_hdr = &pcs_ptr->frm_hdr;
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
        }
        // Process the head of the Picture Decision Reordering Queue (Entry N)
        // P.S. The Picture Decision Reordering Queue should be parsed in the display order to be able to construct a pred structure
        queue_entry_ptr = encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index];

        while (queue_entry_ptr->parent_pcs_wrapper_ptr != NULL) {
            if (((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr->object_ptr))->end_of_sequence_flag == EB_TRUE) {
                frame_passthrough = EB_TRUE;
            }
            else
                frame_passthrough = EB_FALSE;
            window_avail = EB_TRUE;
            previous_entry_index = QUEUE_GET_PREVIOUS_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index);
#if NOISE_BASED_TF_FRAMES
            parent_pcs_window[ 0] = parent_pcs_window[ 1] = parent_pcs_window[ 2] = parent_pcs_window[ 3] = parent_pcs_window[ 4] = parent_pcs_window[ 5] =
            parent_pcs_window[ 6] = parent_pcs_window[ 7] = parent_pcs_window[ 8] = parent_pcs_window[ 9] = parent_pcs_window[10] = parent_pcs_window[11] =
            parent_pcs_window[12] = parent_pcs_window[13] = NULL;
#else
            parent_pcs_window[0] = parent_pcs_window[1] = parent_pcs_window[2] = parent_pcs_window[3] = parent_pcs_window[4] = parent_pcs_window[5] = NULL;
#endif
            //for poc 0, ignore previous frame check
            if (queue_entry_ptr->picture_number > 0 && encode_context_ptr->picture_decision_reorder_queue[previous_entry_index]->parent_pcs_wrapper_ptr == NULL)
                window_avail = EB_FALSE;
            else {

                parent_pcs_window[0] = queue_entry_ptr->picture_number > 0 ? (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[previous_entry_index]->parent_pcs_wrapper_ptr->object_ptr : NULL;
                parent_pcs_window[1] = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index]->parent_pcs_wrapper_ptr->object_ptr;
                for (window_index = 0; window_index < scs_ptr->scd_delay; window_index++) {
                    entry_index = QUEUE_GET_NEXT_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index, window_index + 1);
                    if (encode_context_ptr->picture_decision_reorder_queue[entry_index]->parent_pcs_wrapper_ptr == NULL) {
                        window_avail = EB_FALSE;
                        break;
                    }
                    else if (((PictureParentControlSet *)(encode_context_ptr->picture_decision_reorder_queue[entry_index]->parent_pcs_wrapper_ptr->object_ptr))->end_of_sequence_flag == EB_TRUE) {
                        window_avail = EB_FALSE;
                        frame_passthrough = EB_TRUE;
                        break;
                    }else {
                        parent_pcs_window[2 + window_index] = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[entry_index]->parent_pcs_wrapper_ptr->object_ptr;
                    }
                }
            }

            pcs_ptr = (PictureParentControlSet*)queue_entry_ptr->parent_pcs_wrapper_ptr->object_ptr;
            frm_hdr = &pcs_ptr->frm_hdr;
            pcs_ptr->fade_out_from_black = 0;

            pcs_ptr->fade_in_to_black = 0;
            if (pcs_ptr->idr_flag == EB_TRUE)
                context_ptr->last_solid_color_frame_poc = 0xFFFFFFFF;
            if (window_avail == EB_TRUE && queue_entry_ptr->picture_number > 0) {
                if (scs_ptr->static_config.scene_change_detection) {
                    pcs_ptr->scene_change_flag = scene_transition_detector(
                        context_ptr,
                        scs_ptr,
                        parent_pcs_window,
                        FUTURE_WINDOW_WIDTH);
                }
                else
                    pcs_ptr->scene_change_flag = EB_FALSE;
                pcs_ptr->cra_flag = (pcs_ptr->scene_change_flag == EB_TRUE) ?
                    EB_TRUE :
                    pcs_ptr->cra_flag;

                // Store scene change in context
                context_ptr->is_scene_change_detected = pcs_ptr->scene_change_flag;
            }

            if (window_avail == EB_TRUE || frame_passthrough == EB_TRUE)
            {
                // Place the PCS into the Pre-Assignment Buffer
                // P.S. The Pre-Assignment Buffer is used to store a whole pre-structure
                encode_context_ptr->pre_assignment_buffer[encode_context_ptr->pre_assignment_buffer_count] = queue_entry_ptr->parent_pcs_wrapper_ptr;

                // Setup the PCS & SCS
                pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[encode_context_ptr->pre_assignment_buffer_count]->object_ptr;
                frm_hdr = &pcs_ptr->frm_hdr;
                // Set the POC Number
                pcs_ptr->picture_number = (encode_context_ptr->current_input_poc + 1) /*& ((1 << scs_ptr->bits_for_picture_order_count)-1)*/;
                encode_context_ptr->current_input_poc = pcs_ptr->picture_number;

                pcs_ptr->pred_structure = scs_ptr->static_config.pred_structure;

                pcs_ptr->hierarchical_layers_diff = 0;

                pcs_ptr->init_pred_struct_position_flag = EB_FALSE;

                pcs_ptr->target_bit_rate = scs_ptr->static_config.target_bit_rate;

#if DECOUPLE_ME_RES
                pcs_ptr->self_updated_links = 0;
                pcs_ptr->other_updated_links_cnt = 0;
#endif
                release_prev_picture_from_reorder_queue(
                    encode_context_ptr);

                // If the Intra period length is 0, then introduce an intra for every picture
                if (scs_ptr->intra_period_length == 0)
                    pcs_ptr->cra_flag = EB_TRUE;
                // If an #IntraPeriodLength has passed since the last Intra, then introduce a CRA or IDR based on Intra Refresh type
                else if (scs_ptr->intra_period_length != -1) {
                    pcs_ptr->cra_flag =
                        (scs_ptr->intra_refresh_type != CRA_REFRESH) ?
                        pcs_ptr->cra_flag :
                        (encode_context_ptr->intra_period_position == (uint32_t)scs_ptr->intra_period_length) ?
                        EB_TRUE :
                        pcs_ptr->cra_flag;

                    pcs_ptr->idr_flag =
                        (scs_ptr->intra_refresh_type != IDR_REFRESH) ?
                        pcs_ptr->idr_flag :
                        (encode_context_ptr->intra_period_position == (uint32_t)scs_ptr->intra_period_length) ?
                        EB_TRUE :
                        pcs_ptr->idr_flag;
                }

                encode_context_ptr->pre_assignment_buffer_eos_flag = (pcs_ptr->end_of_sequence_flag) ? (uint32_t)EB_TRUE : encode_context_ptr->pre_assignment_buffer_eos_flag;

                // Increment the Pre-Assignment Buffer Intra Count
                encode_context_ptr->pre_assignment_buffer_intra_count += (pcs_ptr->idr_flag || pcs_ptr->cra_flag);
                encode_context_ptr->pre_assignment_buffer_idr_count += pcs_ptr->idr_flag;
                encode_context_ptr->pre_assignment_buffer_count += 1;

                if (scs_ptr->static_config.rate_control_mode)
                {
                    // Increment the Intra Period Position
                    encode_context_ptr->intra_period_position = (encode_context_ptr->intra_period_position == (uint32_t)scs_ptr->intra_period_length) ? 0 : encode_context_ptr->intra_period_position + 1;
                }
                else
                {
                    // Increment the Intra Period Position
                    encode_context_ptr->intra_period_position = ((encode_context_ptr->intra_period_position == (uint32_t)scs_ptr->intra_period_length) || (pcs_ptr->scene_change_flag == EB_TRUE)) ? 0 : encode_context_ptr->intra_period_position + 1;
                }

                // Determine if Pictures can be released from the Pre-Assignment Buffer
                if ((encode_context_ptr->pre_assignment_buffer_intra_count > 0) ||
                    (encode_context_ptr->pre_assignment_buffer_count == (uint32_t)(1 << scs_ptr->static_config.hierarchical_levels)) ||
                    (encode_context_ptr->pre_assignment_buffer_eos_flag == EB_TRUE) ||
                    (pcs_ptr->pred_structure == EB_PRED_LOW_DELAY_P) ||
                    (pcs_ptr->pred_structure == EB_PRED_LOW_DELAY_B))
                {
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

                                if (encode_context_ptr->pre_assignment_buffer_count >= 32)
                                    context_ptr->mini_gop_activity_array[L6_INDEX] = EB_FALSE;
                                if (encode_context_ptr->pre_assignment_buffer_count >= 16)
                                    context_ptr->mini_gop_activity_array[L5_0_INDEX] = EB_FALSE;
                                if (encode_context_ptr->pre_assignment_buffer_count >= 8) {
                                    context_ptr->mini_gop_activity_array[L4_0_INDEX] = EB_FALSE;
                                    context_ptr->mini_gop_activity_array[L4_1_INDEX] = EB_FALSE;
                                }

                                generate_picture_window_split(
                                        context_ptr,
                                        encode_context_ptr);

                                handle_incomplete_picture_window_map(
                                        context_ptr,
                                        encode_context_ptr);
                            }
                        }
                    }

                    generate_mini_gop_rps(
                        context_ptr,
                        encode_context_ptr);

                    // Loop over Mini GOPs

                    for (mini_gop_index = 0; mini_gop_index < context_ptr->total_number_of_mini_gops; ++mini_gop_index) {
                        pre_assignment_buffer_first_pass_flag = EB_TRUE;
                        encode_context_ptr->is_mini_gop_changed = EB_FALSE;
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
                            pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[out_stride_diff64]->object_ptr;
                            scs_ptr = (SequenceControlSet*)pcs_ptr->scs_wrapper_ptr->object_ptr;
                            frm_hdr = &pcs_ptr->frm_hdr;
                            // Keep track of the mini GOP size to which the input picture belongs - needed @ PictureManagerProcess()
                            pcs_ptr->pre_assignment_buffer_count = context_ptr->mini_gop_length[mini_gop_index];

                            // Update the Pred Structure if cutting short a Random Access period
#if DECOUPLE_ME_RES
                            if (is_pic_cutting_short_ra_mg(context_ptr, pcs_ptr, mini_gop_index))
#else
                            if ((context_ptr->mini_gop_length[mini_gop_index] < pcs_ptr->pred_struct_ptr->pred_struct_period || context_ptr->mini_gop_idr_count[mini_gop_index] > 0) &&
                                pcs_ptr->pred_struct_ptr->pred_type == EB_PRED_RANDOM_ACCESS &&
                                pcs_ptr->idr_flag == EB_FALSE &&
                                pcs_ptr->cra_flag == EB_FALSE)
#endif
                            {
                                // Correct the Pred Index before switching structures
                                if (pre_assignment_buffer_first_pass_flag == EB_TRUE)
                                    encode_context_ptr->pred_struct_position -= pcs_ptr->pred_struct_ptr->init_pic_index;
                                pcs_ptr->pred_struct_ptr = get_prediction_structure(
                                    encode_context_ptr->prediction_structure_group_ptr,
                                    EB_PRED_LOW_DELAY_P,
                                    scs_ptr->reference_count,
                                    pcs_ptr->hierarchical_levels);
                                picture_type = P_SLICE;
                            }
                            // Open GOP CRA - adjust the RPS
                            else if ((context_ptr->mini_gop_length[mini_gop_index] == pcs_ptr->pred_struct_ptr->pred_struct_period) &&

                                (pcs_ptr->pred_struct_ptr->pred_type == EB_PRED_RANDOM_ACCESS || pcs_ptr->pred_struct_ptr->temporal_layer_count == 1) &&
                                pcs_ptr->idr_flag == EB_FALSE &&
                                pcs_ptr->cra_flag == EB_TRUE)
                            {
                                picture_type = I_SLICE;
                            }
                            else {
                                // Set the Picture Type
                                picture_type =
                                    (pcs_ptr->idr_flag) ? I_SLICE :
                                    (pcs_ptr->cra_flag) ? I_SLICE :
                                    (pcs_ptr->pred_structure == EB_PRED_LOW_DELAY_P) ? P_SLICE :
                                    (pcs_ptr->pred_structure == EB_PRED_LOW_DELAY_B) ? B_SLICE :
                                    (pcs_ptr->pre_assignment_buffer_count == pcs_ptr->pred_struct_ptr->pred_struct_period) ? ((out_stride_diff64 == context_ptr->mini_gop_end_index[mini_gop_index] && 0) ? P_SLICE : B_SLICE) :

                                    (encode_context_ptr->pre_assignment_buffer_eos_flag) ? P_SLICE :
                                    B_SLICE;
                            }
                            // If mini GOP switch, reset position
                            encode_context_ptr->pred_struct_position = (pcs_ptr->init_pred_struct_position_flag) ?
                                pcs_ptr->pred_struct_ptr->init_pic_index :
                                encode_context_ptr->pred_struct_position;

                            // If Intra, reset position
                            if (pcs_ptr->idr_flag == EB_TRUE)
                                encode_context_ptr->pred_struct_position = pcs_ptr->pred_struct_ptr->init_pic_index;
                            else if (pcs_ptr->cra_flag == EB_TRUE && context_ptr->mini_gop_length[mini_gop_index] < pcs_ptr->pred_struct_ptr->pred_struct_period)
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
                            // The poc number of the latest IDR picture is stored so that last_idr_picture (present in PCS) for the incoming pictures can be updated.
                            // The last_idr_picture is used in reseting the poc (in entropy coding) whenever IDR is encountered.
                            // Note IMP: This logic only works when display and decode order are the same. Currently for Random Access, IDR is inserted (similar to CRA) by using trailing P pictures (low delay fashion) and breaking prediction structure.
                            // Note: When leading P pictures are implemented, this logic has to change..
                            if (pcs_ptr->idr_flag == EB_TRUE)
                                encode_context_ptr->last_idr_picture = pcs_ptr->picture_number;
                            else
                                pcs_ptr->last_idr_picture = encode_context_ptr->last_idr_picture;
                            // Cycle the PredStructPosition if its overflowed
                            encode_context_ptr->pred_struct_position = (encode_context_ptr->pred_struct_position == pcs_ptr->pred_struct_ptr->pred_struct_entry_count) ?
                                encode_context_ptr->pred_struct_position - pcs_ptr->pred_struct_ptr->pred_struct_period :
                                encode_context_ptr->pred_struct_position;

                            pred_position_ptr = pcs_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[encode_context_ptr->pred_struct_position];
                            if (scs_ptr->static_config.enable_overlays == EB_TRUE) {
                                // At this stage we know the prediction structure and the location of ALT_REF pictures.
                                // For every ALTREF picture, there is an overlay picture. They extra pictures are released
                                // is_alt_ref flag is set for non-slice base layer pictures
                                if (pred_position_ptr->temporal_layer_index == 0 && picture_type != I_SLICE) {
                                    pcs_ptr->is_alt_ref = 1;
                                    frm_hdr->show_frame = 0;
                                    ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[out_stride_diff64 - 1]->object_ptr)->has_show_existing = EB_FALSE;
                                }
                                // release the overlay PCS for non alt ref pictures. First picture does not have overlay PCS
                                else if (pcs_ptr->picture_number) {
                                    eb_release_object(pcs_ptr->overlay_ppcs_ptr->input_picture_wrapper_ptr);
                                    // release the pa_reference_picture
                                    eb_release_object(pcs_ptr->overlay_ppcs_ptr->pa_reference_picture_wrapper_ptr);
                                    // release the parent pcs
                                    eb_release_object(pcs_ptr->overlay_ppcs_ptr->p_pcs_wrapper_ptr);
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
                                ((EbPaReferenceObject*)pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr)->slice_type = pcs_ptr->slice_type;

                                switch (picture_type) {
                                case I_SLICE:

                                    // Reset Prediction Structure Position & Reference Struct Position
                                    if (pcs_ptr->picture_number == 0)
                                        encode_context_ptr->intra_period_position = 0;
                                    encode_context_ptr->elapsed_non_cra_count = 0;

                                    //-------------------------------
                                    // IDR
                                    //-------------------------------
                                    if (pcs_ptr->idr_flag == EB_TRUE) {
                                        // Set CRA flag
                                        pcs_ptr->cra_flag = EB_FALSE;

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
                                        pcs_ptr->cra_flag = EB_TRUE;
                                    }
                                    break;

                                case P_SLICE:
                                case B_SLICE:

                                    // Reset CRA and IDR Flag
                                    pcs_ptr->cra_flag = EB_FALSE;
                                    pcs_ptr->idr_flag = EB_FALSE;

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
                                    pcs_ptr->is_used_as_reference_flag = EB_FALSE;
                                }
                                else {
                                    pcs_ptr->temporal_layer_index = (uint8_t)pred_position_ptr->temporal_layer_index;
                                    pcs_ptr->is_used_as_reference_flag = pred_position_ptr->is_referenced;
                                }

                                pre_assignment_buffer_first_pass_flag = EB_FALSE;

                                // Film grain (assigning the random-seed)
                                {
                                    uint16_t *fgn_random_seed_ptr = &pcs_ptr->scs_ptr->film_grain_random_seed;
                                    frm_hdr->film_grain_params.random_seed = *fgn_random_seed_ptr;
                                    *fgn_random_seed_ptr += 3381;  // Changing random seed for film grain
                                    if (!(*fgn_random_seed_ptr))     // Random seed should not be zero
                                        *fgn_random_seed_ptr += 7391;
                                }
                                uint32_t pic_index = 0;
                                if (scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {
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
                                pcs_ptr->allow_comp_inter_inter = 0;
                                pcs_ptr->is_skip_mode_allowed = 0;

                                frm_hdr->reference_mode = (ReferenceMode)0xFF;

                                if (pcs_ptr->slice_type != I_SLICE) {
                                    pcs_ptr->allow_comp_inter_inter = 1;
                                    if (pcs_ptr->slice_type == P_SLICE) {
                                        pcs_ptr->is_skip_mode_allowed = 0;
                                        frm_hdr->reference_mode = SINGLE_REFERENCE;
                                        pcs_ptr->skip_mode_flag = 0;
                                    }
                                    else if (pcs_ptr->temporal_layer_index == 0) {
                                        frm_hdr->reference_mode = REFERENCE_MODE_SELECT;
                                        frm_hdr->skip_mode_params.skip_mode_flag = 0;
                                    }
                                    else {
                                        frm_hdr->reference_mode = REFERENCE_MODE_SELECT;
                                        pcs_ptr->is_skip_mode_allowed = 1;
                                        pcs_ptr->skip_mode_flag = 1;
                                    }
                                }

                                pcs_ptr->av1_cm->mi_cols = pcs_ptr->aligned_width >> MI_SIZE_LOG2;
                                pcs_ptr->av1_cm->mi_rows = pcs_ptr->aligned_height >> MI_SIZE_LOG2;

                                //Jing: For low delay b/P case, don't alter the bias
                                memset(pcs_ptr->av1_cm->ref_frame_sign_bias, 0, 8 * sizeof(int32_t));
                                if (frm_hdr->reference_mode == REFERENCE_MODE_SELECT &&
                                        pcs_ptr->temporal_layer_index &&
                                        scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)
                                {
                                    pcs_ptr->av1_cm->ref_frame_sign_bias[ALTREF_FRAME] =
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[ALTREF2_FRAME] =
                                        pcs_ptr->av1_cm->ref_frame_sign_bias[BWDREF_FRAME] = 1;
                                }
                                if (frm_hdr->reference_mode == REFERENCE_MODE_SELECT &&
                                    pcs_ptr->temporal_layer_index &&
                                    scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)
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

                                if (pcs_ptr->slice_type == I_SLICE)
                                    context_ptr->last_i_picture_sc_detection = pcs_ptr->sc_content_detected;
                                else
                                    pcs_ptr->sc_content_detected = context_ptr->last_i_picture_sc_detection;

                                // TODO: put this in EbMotionEstimationProcess?
                                // ME Kernel Multi-Processes Signal(s) derivation
#if TF_LEVELS
                                signal_derivation_multi_processes_oq(
                                    scs_ptr,
                                    pcs_ptr,
                                    context_ptr);
#else
                                signal_derivation_multi_processes_oq(
                                scs_ptr,
                                    pcs_ptr);
#endif

                            // Set tx_mode
                            frm_hdr->tx_mode = (pcs_ptr->tx_size_search_mode) ?
                                TX_MODE_SELECT :
                                TX_MODE_LARGEST;
                                // Update the Dependant List Count - If there was an I-frame or Scene Change, then cleanup the Picture Decision PA Reference Queue Dependent Counts
                                if (pcs_ptr->slice_type == I_SLICE)
                                {
                                    input_queue_index = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                                    while (input_queue_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
                                        input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[input_queue_index];

#if DECOUPLE_ME_RES
                                        int32_t diff_n = 0;
#endif
                                        // Modify Dependent List0
                                        dep_list_count = input_entry_ptr->list0.list_count;
                                        for (dep_idx = 0; dep_idx < dep_list_count; ++dep_idx) {
                                            // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                                            // current_input_poc += (current_input_poc < input_entry_ptr->pocNumber) ? (1 << scs_ptr->bits_for_picture_order_count) : 0;

                                            dep_poc = POC_CIRCULAR_ADD(
                                                input_entry_ptr->picture_number, // can't use a value that gets reset
                                                input_entry_ptr->list0.list[dep_idx]/*,
                                                scs_ptr->bits_for_picture_order_count*/);

                                                // If Dependent POC is greater or equal to the IDR POC
                                            if (dep_poc >= pcs_ptr->picture_number && input_entry_ptr->list0.list[dep_idx]) {
                                                input_entry_ptr->list0.list[dep_idx] = 0;

                                                // Decrement the Reference's referenceCount
                                                --input_entry_ptr->dependent_count;

#if DECOUPLE_ME_RES
                                                diff_n--;
#endif

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
                                                input_entry_ptr->picture_number,
                                                input_entry_ptr->list1.list[dep_idx]/*,
                                                scs_ptr->bits_for_picture_order_count*/);

                                                // If Dependent POC is greater or equal to the IDR POC
                                            if (((dep_poc >= pcs_ptr->picture_number) || (((pcs_ptr->pre_assignment_buffer_count != pcs_ptr->pred_struct_ptr->pred_struct_period) || (pcs_ptr->idr_flag == EB_TRUE)) && (dep_poc > (pcs_ptr->picture_number - pcs_ptr->pre_assignment_buffer_count)))) && input_entry_ptr->list1.list[dep_idx]) {
                                                input_entry_ptr->list1.list[dep_idx] = 0;

                                                // Decrement the Reference's referenceCount
                                                --input_entry_ptr->dependent_count;

#if DECOUPLE_ME_RES
                                                diff_n--;
#endif

                                                CHECK_REPORT_ERROR(
                                                    (input_entry_ptr->dependent_count != ~0u),
                                                    encode_context_ptr->app_callback_ptr,
                                                    EB_ENC_PD_ERROR3);
                                            }
                                        }

#if DECOUPLE_ME_RES
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
#endif
                                        // Increment the input_queue_index Iterator
                                        input_queue_index = (input_queue_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : input_queue_index + 1;
                                    }
                                }
                                else if (pcs_ptr->idr_flag == EB_TRUE) {
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
#if !DECOUPLE_ME_RES
                                    input_entry_ptr->reference_entry_index = encode_context_ptr->picture_decision_pa_reference_queue_tail_index;
#endif
                                    input_entry_ptr->is_alt_ref = pcs_ptr->is_alt_ref;
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

#if ON_OFF_FEATURE_MRP // one function or under signal derivation ?
                                // set number of references to try based on mrp level
                                set_mrp_controls(pcs_ptr);
#else

#if MRP_CTRL
                                //set the number of references to try in ME/MD.Note: PicMgr will still use the original values to sync the references.
#if UPGRADE_M6_M7_M8
#if !SYNCH_MRP
                                if (pcs_ptr->sc_content_detected) {
#if NEW_MRP_SETTINGS
                                    if (MRS_MODE) {
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 3);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
                                    }
                                    else if (MR_MODE) {
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 3);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 2);
                                    }
#if JUNE17_ADOPTIONS
                                    else if (pcs_ptr->enc_mode <= ENC_M4) {
#else
                                    else if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 2);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 2);
                                    }
#if JUNE17_ADOPTIONS
                                    else if (pcs_ptr->enc_mode <= ENC_M5) {
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 2);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 1);
                                    }
#endif
                                    else {
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 1);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 1);
                                    }
#else
#if MAY12_ADOPTIONS
#if PRESET_SHIFITNG
                                    if (pcs_ptr->enc_mode <= ENC_M1) {
#else
                                    if (pcs_ptr->enc_mode <= ENC_M2) {
#endif
#else
#if SHIFT_M6_SC_TO_M5
                                    if (pcs_ptr->enc_mode <= ENC_M4) {
#else
#if APR25_3AM_ADOPTIONS
                                    if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if MRP_ADOPTIONS
                                    if (pcs_ptr->enc_mode <= ENC_M6) {
#else
                                    if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
#endif
#endif
#endif
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 4);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
                                    }
#if APR25_12AM_ADOPTIONS
#if MAY19_ADOPTIONS
#if PRESET_SHIFITNG
                                    else if (pcs_ptr->enc_mode <= ENC_M2) {
#else
                                    else if (pcs_ptr->enc_mode <= ENC_M4) {
#endif
#else
#if APR25_3AM_ADOPTIONS
#if APR23_4AM_M6_ADOPTIONS
                                    else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
                                    else if (pcs_ptr->enc_mode <= ENC_M6) {
#endif
#else
                                    else if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
#endif
#if APR25_10AM_ADOPTIONS
                                        pcs_ptr->ref_list0_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list0_count, 2) : MIN(pcs_ptr->ref_list0_count, 1);
                                        pcs_ptr->ref_list1_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list1_count, 2) : MIN(pcs_ptr->ref_list1_count, 1);
#else
                                        pcs_ptr->ref_list0_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list0_count, 4) : MIN(pcs_ptr->ref_list0_count, 1);
                                        pcs_ptr->ref_list1_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list1_count, 3) : MIN(pcs_ptr->ref_list1_count, 1);
#endif
                                    }
#endif
                                    else {
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 1);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 1);
                                    }
#endif
                                }
                                else {
#endif
#if MRP_ADOPTIONS
#if JUNE17_ADOPTIONS
#if SYNCH_MRP
                                    if (pcs_ptr->enc_mode <= ENC_M6) {
#else
                                    if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
#else
#if PRESET_SHIFITNG
                                    if (pcs_ptr->enc_mode <= ENC_M4) {
#else
                                    if (pcs_ptr->enc_mode <= ENC_M6) {
#endif
#endif
#else
                                    if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 4);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
                                    }
#if APR25_12AM_ADOPTIONS
#if NEW_MRP_SETTINGS
#if JUNE17_ADOPTIONS
#if !SYNCH_MRP
                                    else if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
#else
                                    else if (pcs_ptr->enc_mode <= ENC_M6) {
#endif
#if !SYNCH_MRP
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 2);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 2);
#endif
#else
#if PRESET_SHIFITNG
                                    else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
                                    else if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
#if MAY19_ADOPTIONS
                                        pcs_ptr->ref_list0_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list0_count, 4) : MIN(pcs_ptr->ref_list0_count, 2);
                                        pcs_ptr->ref_list1_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list1_count, 3) : MIN(pcs_ptr->ref_list1_count, 1);
#else
                                        pcs_ptr->ref_list0_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list0_count, 4) : MIN(pcs_ptr->ref_list0_count, 1);
                                        pcs_ptr->ref_list1_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list1_count, 3) : MIN(pcs_ptr->ref_list1_count, 1);
#endif
#endif
#if !SYNCH_MRP
                                    }
#endif
#endif
                                    else {
#if APR25_7PM_ADOPTIONS
                                        pcs_ptr->ref_list0_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list0_count, 2) : MIN(pcs_ptr->ref_list0_count, 1);
                                        pcs_ptr->ref_list1_count_try = pcs_ptr->is_used_as_reference_flag ? MIN(pcs_ptr->ref_list1_count, 2) : MIN(pcs_ptr->ref_list1_count, 1);
#else
                                        pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 1);
                                        pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 1);
#endif
                                    }
#if !SYNCH_MRP
                                }
#endif
#else
                                // Set the number of references to try in ME/MD. Note: PicMgr/RPS will still use the original values to sync the references.
                                if (pcs_ptr->sc_content_detected) {
                                    pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 4);
                                    pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
                                }
                                else {
                                    pcs_ptr->ref_list0_count_try = MIN(pcs_ptr->ref_list0_count, 4);
                                    pcs_ptr->ref_list1_count_try = MIN(pcs_ptr->ref_list1_count, 3);
                                }
#endif
                                assert(pcs_ptr->ref_list0_count_try <= pcs_ptr->ref_list0_count);
                                assert(pcs_ptr->ref_list1_count_try <= pcs_ptr->ref_list1_count);
#endif
#endif
                                if (!pcs_ptr->is_overlay) {
                                    input_entry_ptr->list0_ptr = &pred_position_ptr->ref_list0;
                                    input_entry_ptr->list1_ptr = &pred_position_ptr->ref_list1;
                                }
                                if (!pcs_ptr->is_overlay)
                                {
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

                                    ((EbPaReferenceObject*)pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr)->dependent_pictures_count = input_entry_ptr->dependent_count;
                                }

                                CHECK_REPORT_ERROR(
                                    (pcs_ptr->pred_struct_ptr->pred_struct_period * REF_LIST_MAX_DEPTH < MAX_ELAPSED_IDR_COUNT),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PD_ERROR5);

                                // Reset the PA Reference Lists
                                EB_MEMSET(pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                                EB_MEMSET(pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));

                                EB_MEMSET(pcs_ptr->ref_pic_poc_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                                EB_MEMSET(pcs_ptr->ref_pic_poc_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                            }
                            pcs_ptr = cur_picture_control_set_ptr;
#if TF_LEVELS
                            if(context_ptr->tf_ctrls.enabled) {
#else
                            uint8_t perform_filtering =
                                (scs_ptr->enable_altrefs == EB_TRUE && scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS &&
                                 scs_ptr->static_config.hierarchical_levels >= 2) &&  pcs_ptr->sc_content_detected == 0 &&
                                (pcs_ptr->temporal_layer_index == 0 ||
                                  (pcs_ptr->temporal_layer_index == 1 && scs_ptr->static_config.hierarchical_levels >= 3))
                                ? 1 : 0;

                            if (perform_filtering){
#endif
                                derive_tf_window_params(
                                    scs_ptr,
                                    encode_context_ptr,
                                    pcs_ptr,
#if TF_LEVELS
                                    context_ptr,
#endif
                                    out_stride_diff64);
                                pcs_ptr->temp_filt_prep_done = 0;

                                // Start Filtering in ME processes
                                {
                                    int16_t seg_idx;

                                    // Initialize Segments
                                    pcs_ptr->tf_segments_column_count = scs_ptr->tf_segment_column_count;
                                    pcs_ptr->tf_segments_row_count    = scs_ptr->tf_segment_row_count;
                                    pcs_ptr->tf_segments_total_count = (uint16_t)(pcs_ptr->tf_segments_column_count  * pcs_ptr->tf_segments_row_count);
                                    pcs_ptr->temp_filt_seg_acc = 0;
                                    if (pcs_ptr->temporal_layer_index == 0)
                                        pcs_ptr->altref_strength = scs_ptr->static_config.altref_strength;
                                    else
                                        pcs_ptr->altref_strength = 2;

                                    for (seg_idx = 0; seg_idx < pcs_ptr->tf_segments_total_count; ++seg_idx) {
                                        eb_get_empty_object(
                                            context_ptr->picture_decision_results_output_fifo_ptr,
                                            &out_results_wrapper_ptr);
                                        out_results_ptr = (PictureDecisionResults*)out_results_wrapper_ptr->object_ptr;
                                        out_results_ptr->pcs_wrapper_ptr = encode_context_ptr->pre_assignment_buffer[out_stride_diff64];
                                        out_results_ptr->segment_index = seg_idx;
                                        out_results_ptr->task_type = 1;
                                        eb_post_full_object(out_results_wrapper_ptr);
                                    }

                                    eb_block_on_semaphore(pcs_ptr->temp_filt_done_semaphore);
                                }

                            }else
                                pcs_ptr->temporal_filtering_on = EB_FALSE; // set temporal filtering flag OFF for current picture

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
                                frm_hdr = &pcs_ptr->frm_hdr;
                            }
                            else {
                                pcs_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[out_stride_diff64]->object_ptr;
                                frm_hdr = &pcs_ptr->frm_hdr;
                            }
                            pcs_ptr->picture_number_alt = encode_context_ptr->picture_number_alt++;

                            // Set the Decode Order
                            if ((context_ptr->mini_gop_idr_count[mini_gop_index] == 0) &&
                                (context_ptr->mini_gop_length[mini_gop_index] == pcs_ptr->pred_struct_ptr->pred_struct_period) &&
                                (scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) &&
                                !pcs_ptr->is_overlay) {
                                pcs_ptr->decode_order = encode_context_ptr->decode_base_number + pcs_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[pcs_ptr->pred_struct_index]->decode_order;
                            }
                            else
                                pcs_ptr->decode_order = pcs_ptr->picture_number_alt;
                            encode_context_ptr->terminating_sequence_flag_received = (pcs_ptr->end_of_sequence_flag == EB_TRUE) ?
                                EB_TRUE :
                                encode_context_ptr->terminating_sequence_flag_received;

                            encode_context_ptr->terminating_picture_number = (pcs_ptr->end_of_sequence_flag == EB_TRUE) ?
                                pcs_ptr->picture_number_alt :
                                encode_context_ptr->terminating_picture_number;

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
#if DECOUPLE_ME_RES
                                        if (pcs_ptr->is_overlay)
                                        // hardcode the reference for the overlay frame
                                            ref_poc = pcs_ptr->picture_number;
                                        else
                                            ref_poc = POC_CIRCULAR_ADD(
                                                pcs_ptr->picture_number,
                                                -input_entry_ptr->list0_ptr->reference_list[ref_pic_index]);

                                        pa_reference_entry_ptr = search_ref_in_ref_queue_pa(encode_context_ptr, ref_poc);
                                        assert(pa_reference_entry_ptr != 0);
                                        CHECK_REPORT_ERROR((pa_reference_entry_ptr),
                                            encode_context_ptr->app_callback_ptr,
                                            EB_ENC_PM_ERROR10);
#else
                                        // hardcode the reference for the overlay frame
                                        if(pcs_ptr->is_overlay){
                                            pa_reference_queue_index = (uint32_t)CIRCULAR_ADD(
                                                (int32_t)input_entry_ptr->reference_entry_index,
                                                PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                            pa_reference_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[pa_reference_queue_index];

                                            // Calculate the Ref POC
                                            ref_poc = POC_CIRCULAR_ADD(
                                                pcs_ptr->picture_number,
                                                0);
                                        }
                                        else {
                                            pa_reference_queue_index = (uint32_t)CIRCULAR_ADD(
                                                ((int32_t)input_entry_ptr->reference_entry_index) - input_entry_ptr->list0_ptr->reference_list[ref_pic_index],
                                                PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                            pa_reference_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[pa_reference_queue_index];

                                            // Calculate the Ref POC
                                            ref_poc = POC_CIRCULAR_ADD(
                                                pcs_ptr->picture_number,
                                                -input_entry_ptr->list0_ptr->reference_list[ref_pic_index] /*
                                                scs_ptr->bits_for_picture_order_count*/);
                                        }
#endif
                                            // Set the Reference Object
                                        pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_0][ref_pic_index] = pa_reference_entry_ptr->input_object_ptr;
                                        pcs_ptr->ref_pic_poc_array[REF_LIST_0][ref_pic_index] = ref_poc;
                                        // Increment the PA Reference's liveCount by the number of tiles in the input picture
                                        eb_object_inc_live_count(
                                            pa_reference_entry_ptr->input_object_ptr,
                                            1);
                                        --pa_reference_entry_ptr->dependent_count;
                                    }
                                }
                            }

                            // Configure List1
                            if (pcs_ptr->slice_type == B_SLICE) {
                                uint8_t ref_pic_index;
                                for (ref_pic_index = 0; ref_pic_index < pcs_ptr->ref_list1_count; ++ref_pic_index) {
                                    if (pcs_ptr->ref_list1_count) {
#if DECOUPLE_ME_RES
                                        ref_poc = POC_CIRCULAR_ADD(
                                            pcs_ptr->picture_number,
                                            -input_entry_ptr->list1_ptr->reference_list[ref_pic_index]);

                                        pa_reference_entry_ptr = search_ref_in_ref_queue_pa(encode_context_ptr, ref_poc);

                                        assert(pa_reference_entry_ptr != 0);
                                        CHECK_REPORT_ERROR((pa_reference_entry_ptr),
                                            encode_context_ptr->app_callback_ptr,
                                            EB_ENC_PM_ERROR10);
#else
                                        pa_reference_queue_index = (uint32_t)CIRCULAR_ADD(
                                            ((int32_t)input_entry_ptr->reference_entry_index) - input_entry_ptr->list1_ptr->reference_list[ref_pic_index],
                                            PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                        pa_reference_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[pa_reference_queue_index];

                                        // Calculate the Ref POC
                                        ref_poc = POC_CIRCULAR_ADD(
                                            pcs_ptr->picture_number,
                                            -input_entry_ptr->list1_ptr->reference_list[ref_pic_index]/*,
                                            scs_ptr->bits_for_picture_order_count*/);
#endif
                                        // Set the Reference Object
                                        pcs_ptr->ref_pa_pic_ptr_array[REF_LIST_1][ref_pic_index] = pa_reference_entry_ptr->input_object_ptr;
                                        pcs_ptr->ref_pic_poc_array[REF_LIST_1][ref_pic_index] = ref_poc;

                                        // Increment the PA Reference's liveCount by the number of tiles in the input picture
                                        eb_object_inc_live_count(
                                            pa_reference_entry_ptr->input_object_ptr,
                                            1);
                                        --pa_reference_entry_ptr->dependent_count;
                                    }
                                }
                            }

                            eb_av1_setup_skip_mode_allowed(pcs_ptr);

                            pcs_ptr->is_skip_mode_allowed = frm_hdr->skip_mode_params.skip_mode_allowed;
                            pcs_ptr->skip_mode_flag = pcs_ptr->is_skip_mode_allowed;
                            //SVT_LOG("POC:%i  skip_mode_allowed:%i  REF_SKIP_0: %i   REF_SKIP_1: %i \n",pcs_ptr->picture_number, pcs_ptr->skip_mode_info.skip_mode_allowed, pcs_ptr->skip_mode_info.ref_frame_idx_0, pcs_ptr->skip_mode_info.ref_frame_idx_1);

                            {
                                if (pcs_ptr->ref_list0_count)
                                    pcs_ptr->scene_transition_flag[REF_LIST_0] = EB_FALSE;
                                if (pcs_ptr->ref_list1_count)
                                    pcs_ptr->scene_transition_flag[REF_LIST_1] = EB_FALSE;
                            }

                            //set the ref frame types used for this picture,
#if  !REMOVE_MRP_MODE || !MRP_CTRL
                            set_all_ref_frame_type(scs_ptr, pcs_ptr, pcs_ptr->ref_frame_type_arr, &pcs_ptr->tot_ref_frame_types);
#else
                            set_all_ref_frame_type(pcs_ptr, pcs_ptr->ref_frame_type_arr, &pcs_ptr->tot_ref_frame_types);
#endif
                            // Initialize Segments
                            pcs_ptr->me_segments_column_count = (uint8_t)(scs_ptr->me_segment_column_count_array[pcs_ptr->temporal_layer_index]);
                            pcs_ptr->me_segments_row_count = (uint8_t)(scs_ptr->me_segment_row_count_array[pcs_ptr->temporal_layer_index]);
                            pcs_ptr->me_segments_total_count = (uint16_t)(pcs_ptr->me_segments_column_count  * pcs_ptr->me_segments_row_count);
                            pcs_ptr->me_segments_completion_mask = 0;

                            //****************************************************
                            // Picture resizing for super-res tool
                            //****************************************************

                            // Scale picture if super-res is used
                            if(scs_ptr->static_config.superres_mode > SUPERRES_NONE){
                                init_resize_picture(pcs_ptr->scs_ptr,
                                                    pcs_ptr);
                            }
#if DECOUPLE_ME_RES
                            uint32_t pic_it = out_stride_diff64 - context_ptr->mini_gop_start_index[mini_gop_index];
                            context_ptr->mg_pictures_array[pic_it] = pcs_ptr;
#else
                            //****************************************************

                            // Post the results to the ME processes
                            {
                                uint32_t segment_index;

                                for (segment_index = 0; segment_index < pcs_ptr->me_segments_total_count; ++segment_index) {
                                    // Get Empty Results Object
                                    eb_get_empty_object(
                                            context_ptr->picture_decision_results_output_fifo_ptr,
                                            &out_results_wrapper_ptr);

                                    out_results_ptr = (PictureDecisionResults*)out_results_wrapper_ptr->object_ptr;
                                    if (pcs_ptr->is_overlay)
                                        out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
                                    else
                                        out_results_ptr->pcs_wrapper_ptr = encode_context_ptr->pre_assignment_buffer[out_stride_diff64];

                                    out_results_ptr->segment_index = segment_index;
                                    out_results_ptr->task_type = 0;
                                    // Post the Full Results Object
                                    eb_post_full_object(out_results_wrapper_ptr);
                                }
                            }
#endif
                            if (out_stride_diff64 == context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay) {
                                // Increment the Decode Base Number
                                encode_context_ptr->decode_base_number += context_ptr->mini_gop_length[mini_gop_index] + has_overlay;
                            }
                            if (out_stride_diff64 == encode_context_ptr->pre_assignment_buffer_count - 1 + has_overlay) {

                                // Reset the Pre-Assignment Buffer
                                encode_context_ptr->pre_assignment_buffer_count = 0;
                                encode_context_ptr->pre_assignment_buffer_idr_count = 0;
                                encode_context_ptr->pre_assignment_buffer_intra_count = 0;
                                encode_context_ptr->pre_assignment_buffer_scene_change_count = 0;
                                encode_context_ptr->pre_assignment_buffer_eos_flag = EB_FALSE;
                            }
                        }
#if DECOUPLE_ME_RES
                        uint32_t mg_size = context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay - context_ptr->mini_gop_start_index[mini_gop_index]+1;

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


                        for (uint32_t pic_i = 0; pic_i < mg_size; ++pic_i){

                            pcs_ptr = context_ptr->mg_pictures_array[pic_i];

                            if (scs_ptr->static_config.look_ahead_distance == 0) {

                                EbObjectWrapper* reference_picture_wrapper_ptr;
                                // Get Empty Reference Picture Object
                                eb_get_empty_object(
                                    scs_ptr->encode_context_ptr->reference_picture_pool_fifo_ptr,
                                    &reference_picture_wrapper_ptr);
                                pcs_ptr->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                                // Give the new Reference a nominal live_count of 1
                                eb_object_inc_live_count(pcs_ptr->reference_picture_wrapper_ptr, 1);
                            }
                            //get a new ME data buffer
                            eb_get_empty_object(context_ptr->me_fifo_ptr,
                                &me_wrapper_ptr);
                            pcs_ptr->me_data_wrapper_ptr = me_wrapper_ptr;

                            pcs_ptr->pa_me_data = (MotionEstimationData *)me_wrapper_ptr->object_ptr;

                            for (uint32_t segment_index = 0; segment_index < pcs_ptr->me_segments_total_count; ++segment_index) {
                                // Get Empty Results Object
                                eb_get_empty_object(
                                    context_ptr->picture_decision_results_output_fifo_ptr,
                                    &out_results_wrapper_ptr);

                                out_results_ptr = (PictureDecisionResults*)out_results_wrapper_ptr->object_ptr;
                                //if (pcs_ptr->is_overlay)
                                    out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
                                //else
                               //     out_results_ptr->pcs_wrapper_ptr = encode_context_ptr->pre_assignment_buffer[out_stride_diff64];

                                out_results_ptr->segment_index = segment_index;
                                out_results_ptr->task_type = 0;
                                // Post the Full Results Object
                                eb_post_full_object(out_results_wrapper_ptr);
                            }


                        }
#endif
                    } // End MINI GOPs loop
                }

                // Walk the picture_decision_pa_reference_queue and remove entries that have been completely referenced.
                input_queue_index = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                while (input_queue_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
                    input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[input_queue_index];

                    // Remove the entry
                    if ((input_entry_ptr->dependent_count == 0) &&
                        (input_entry_ptr->input_object_ptr)) {
                        // Release the nominal live_count value
                        eb_release_object(input_entry_ptr->input_object_ptr);
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
            if (window_avail == EB_FALSE && frame_passthrough == EB_FALSE)
                break;
        }

        // Release the Input Results
        eb_release_object(in_results_wrapper_ptr);
    }

    return NULL;
}
// clang-format on
