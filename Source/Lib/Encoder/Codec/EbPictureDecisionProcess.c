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

/**************************************
 * Context
 **************************************/
typedef struct PictureDecisionContext
{
    EbDctor      dctor;
    EbFifo       *picture_analysis_results_input_fifo_ptr;
    EbFifo       *picture_decision_results_output_fifo_ptr;

    uint64_t      last_solid_color_frame_poc;

    EbBool        reset_running_avg;

    uint32_t    **ahd_running_avg_cb;
    uint32_t    **ahd_running_avg_cr;
    uint32_t    **ahd_running_avg;
    EbBool        is_scene_change_detected;

    // Dynamic GOP
    uint32_t      totalRegionActivityCost[MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT];

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
} PictureDecisionContext;

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
        ref_frame_arr_single[i].poc = parent_pcs_ptr->av1_ref_signal.ref_poc_array[i] % (uint64_t)(1 << (parent_pcs_ptr->sequence_control_set_ptr->seq_header.order_hint_info.order_hint_bits));
    }

    OrderHintInfo order_hint_info_st;
    order_hint_info_st.enable_order_hint = 1;
    order_hint_info_st.order_hint_bits = 6+1;

    const OrderHintInfo *const order_hint_info = &order_hint_info_st;// cm->seq_params.order_hint_info;
    SkipModeInfo *const skip_mode_info = &frm_hdr->skip_mode_params;// cm->current_frame.skip_mode_info;

    skip_mode_info->skip_mode_allowed = 0;
    skip_mode_info->ref_frame_idx_0 = INVALID_IDX;
    skip_mode_info->ref_frame_idx_1 = INVALID_IDX;
    parent_pcs_ptr->cur_order_hint = parent_pcs_ptr->picture_number % (uint64_t)(1 << (parent_pcs_ptr->sequence_control_set_ptr->seq_header.order_hint_info.order_hint_bits));
    for (uint8_t i = 0; i < 7; ++i)
        parent_pcs_ptr->ref_order_hint[i] = (uint32_t)ref_frame_arr_single[i].poc;
    if (/*!order_hint_info->enable_order_hint ||*/ parent_pcs_ptr->slice_type == I_SLICE /*frame_is_intra_only(cm)*/ ||
        frm_hdr->reference_mode == SINGLE_REFERENCE)
        return;

    const int cur_order_hint = parent_pcs_ptr->picture_number % (uint64_t)(1 << (parent_pcs_ptr->sequence_control_set_ptr->seq_header.order_hint_info.order_hint_bits));
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

uint8_t  circ_dec(uint8_t max, uint8_t off, uint8_t input)
{
    int8_t x = input;

    x--;
    if (x < 0)
        x = max;

    if (off == 2)
    {
        x--;
        if (x < 0)
            x = max;
    }

    return x;
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

#define FUTURE_WINDOW_WIDTH                 6
#define FLASH_TH                            5
#define FADE_TH                             3
#define SCENE_TH                            3000
#define NOISY_SCENE_TH                      4500    // SCD TH in presence of noise
#define HIGH_PICTURE_VARIANCE_TH            1500
#define NUM64x64INPIC(w,h)          ((w*h)>> (LOG2F(BLOCK_SIZE_64)<<1))
#define QUEUE_GET_PREVIOUS_SPOT(h)  ((h == 0) ? PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1 : h - 1)
#define QUEUE_GET_NEXT_SPOT(h,off)  (( (h+off) >= PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH) ? h+off - PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH  : h + off)

#define WTH 64
#define OTH 64
#if !MULTI_PASS_PD
#define FC_SKIP_TX_SR_TH025                     125 // Fast cost skip tx search threshold.
#define FC_SKIP_TX_SR_TH010                     110 // Fast cost skip tx search threshold.
#endif
void picture_decision_context_dctor(EbPtr p)
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
    uint32_t arrayRow, arrowColumn;

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

    for (arrayRow = 0; arrayRow < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; arrayRow++)
    {
        for (arrowColumn = 0; arrowColumn < MAX_NUMBER_OF_REGIONS_IN_WIDTH; arrowColumn++) {
            context_ptr->ahd_running_avg_cb[arrowColumn][arrayRow] = 0;
            context_ptr->ahd_running_avg_cr[arrowColumn][arrayRow] = 0;
            context_ptr->ahd_running_avg[arrowColumn][arrayRow] = 0;
        }
    }

    context_ptr->reset_running_avg = EB_TRUE;

    return EB_ErrorNone;
}

EbBool SceneTransitionDetector(
    PictureDecisionContext *context_ptr,
    SequenceControlSet                 *sequence_control_set_ptr,
    PictureParentControlSet           **ParentPcsWindow,
    uint32_t                                windowWidthFuture)
{
    PictureParentControlSet       *previousPictureControlSetPtr = ParentPcsWindow[0];
    PictureParentControlSet       *currentPictureControlSetPtr = ParentPcsWindow[1];
    PictureParentControlSet       *futurePictureControlSetPtr = ParentPcsWindow[2];

    // calculating the frame threshold based on the number of 64x64 blocks in the frame
    uint32_t  regionThreshHold;
    uint32_t  regionThreshHoldChroma;
    // this variable determines whether the running average should be reset to equal the ahd or not after detecting a scene change.
    //EbBool reset_running_avg = context_ptr->reset_running_avg;

    EbBool isAbruptChange; // this variable signals an abrubt change (scene change or flash)
    EbBool is_scene_change; // this variable signals a frame representing a scene change
    EbBool isFlash; // this variable signals a frame that contains a flash
    EbBool isFade; // this variable signals a frame that contains a fade
    EbBool gradualChange; // this signals the detection of a light scene change a small/localized flash or the start of a fade

    uint32_t  ahd; // accumulative histogram (absolute) differences between the past and current frame

    uint32_t  ahdCb;
    uint32_t  ahdCr;

    uint32_t  ahdErrorCb = 0;
    uint32_t  ahdErrorCr = 0;

    uint32_t **ahd_running_avg_cb = context_ptr->ahd_running_avg_cb;
    uint32_t **ahd_running_avg_cr = context_ptr->ahd_running_avg_cr;
    uint32_t **ahd_running_avg = context_ptr->ahd_running_avg;

    uint32_t  ahdError = 0; // the difference between the ahd and the running average at the current frame.

    uint8_t   aidFuturePast = 0; // this variable denotes the average intensity difference between the next and the past frames
    uint8_t   aidFuturePresent = 0;
    uint8_t   aidPresentPast = 0;

    uint32_t  bin = 0; // variable used to iterate through the bins of the histograms

    uint32_t  regionInPictureWidthIndex;
    uint32_t  regionInPictureHeightIndex;

    uint32_t  regionWidth;
    uint32_t  regionHeight;
    uint32_t  regionWidthOffset;
    uint32_t  regionHeightOffset;

    uint32_t  isAbruptChangeCount = 0;
    uint32_t  isSceneChangeCount = 0;

    uint32_t  regionCountThreshold = (sequence_control_set_ptr->scd_mode == SCD_MODE_2) ?
        (uint32_t)(((float)((sequence_control_set_ptr->picture_analysis_number_of_regions_per_width * sequence_control_set_ptr->picture_analysis_number_of_regions_per_height) * 75) / 100) + 0.5) :
        (uint32_t)(((float)((sequence_control_set_ptr->picture_analysis_number_of_regions_per_width * sequence_control_set_ptr->picture_analysis_number_of_regions_per_height) * 50) / 100) + 0.5);

    regionWidth = ParentPcsWindow[1]->enhanced_picture_ptr->width / sequence_control_set_ptr->picture_analysis_number_of_regions_per_width;
    regionHeight = ParentPcsWindow[1]->enhanced_picture_ptr->height / sequence_control_set_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_width; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_height; regionInPictureHeightIndex++) { // loop over vertical regions

            isAbruptChange = EB_FALSE;
            is_scene_change = EB_FALSE;
            isFlash = EB_FALSE;
            gradualChange = EB_FALSE;

            // Reset accumulative histogram (absolute) differences between the past and current frame
            ahd = 0;
            ahdCb = 0;
            ahdCr = 0;

            regionWidthOffset = (regionInPictureWidthIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_width - 1) ?
                ParentPcsWindow[1]->enhanced_picture_ptr->width - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_width * regionWidth) :
                0;

            regionHeightOffset = (regionInPictureHeightIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_height - 1) ?
                ParentPcsWindow[1]->enhanced_picture_ptr->height - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_height * regionHeight) :
                0;

            regionWidth += regionWidthOffset;
            regionHeight += regionHeightOffset;

            regionThreshHold = (
                // Noise insertion/removal detection
                ((ABS((int64_t)currentPictureControlSetPtr->pic_avg_variance - (int64_t)previousPictureControlSetPtr->pic_avg_variance)) > NOISE_VARIANCE_TH) &&
                (currentPictureControlSetPtr->pic_avg_variance > HIGH_PICTURE_VARIANCE_TH || previousPictureControlSetPtr->pic_avg_variance > HIGH_PICTURE_VARIANCE_TH)) ?
                NOISY_SCENE_TH * NUM64x64INPIC(regionWidth, regionHeight) : // SCD TH function of noise insertion/removal.
                SCENE_TH * NUM64x64INPIC(regionWidth, regionHeight);

            regionThreshHoldChroma = regionThreshHold / 4;

            for (bin = 0; bin < HISTOGRAM_NUMBER_OF_BINS; ++bin) {
                ahd += ABS((int32_t)currentPictureControlSetPtr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][bin] - (int32_t)previousPictureControlSetPtr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][bin]);
                ahdCb += ABS((int32_t)currentPictureControlSetPtr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1][bin] - (int32_t)previousPictureControlSetPtr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1][bin]);
                ahdCr += ABS((int32_t)currentPictureControlSetPtr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2][bin] - (int32_t)previousPictureControlSetPtr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2][bin]);
            }

            if (context_ptr->reset_running_avg) {
                ahd_running_avg[regionInPictureWidthIndex][regionInPictureHeightIndex] = ahd;
                ahd_running_avg_cb[regionInPictureWidthIndex][regionInPictureHeightIndex] = ahdCb;
                ahd_running_avg_cr[regionInPictureWidthIndex][regionInPictureHeightIndex] = ahdCr;
            }

            ahdError = ABS((int32_t)ahd_running_avg[regionInPictureWidthIndex][regionInPictureHeightIndex] - (int32_t)ahd);
            ahdErrorCb = ABS((int32_t)ahd_running_avg_cb[regionInPictureWidthIndex][regionInPictureHeightIndex] - (int32_t)ahdCb);
            ahdErrorCr = ABS((int32_t)ahd_running_avg_cr[regionInPictureWidthIndex][regionInPictureHeightIndex] - (int32_t)ahdCr);

            if ((ahdError > regionThreshHold       && ahd >= ahdError) ||
                (ahdErrorCb > regionThreshHoldChroma && ahdCb >= ahdErrorCb) ||
                (ahdErrorCr > regionThreshHoldChroma && ahdCr >= ahdErrorCr)) {
                isAbruptChange = EB_TRUE;
            }
            else if ((ahdError > (regionThreshHold >> 1)) && ahd >= ahdError)
                gradualChange = EB_TRUE;
            if (isAbruptChange)
            {
                aidFuturePast = (uint8_t)ABS((int16_t)futurePictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] - (int16_t)previousPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0]);
                aidFuturePresent = (uint8_t)ABS((int16_t)futurePictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] - (int16_t)currentPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0]);
                aidPresentPast = (uint8_t)ABS((int16_t)currentPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] - (int16_t)previousPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0]);

                if (aidFuturePast < FLASH_TH && aidFuturePresent >= FLASH_TH && aidPresentPast >= FLASH_TH) {
                    isFlash = EB_TRUE;
                    //SVT_LOG ("\nFlash in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                }
                else if (aidFuturePresent < FADE_TH && aidPresentPast < FADE_TH) {
                    isFade = EB_TRUE;
                    //SVT_LOG ("\nFlash in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                }
                else {
                    is_scene_change = EB_TRUE;
                    //SVT_LOG ("\nScene Change in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                }
            }
            else if (gradualChange) {
                aidFuturePast = (uint8_t)ABS((int16_t)futurePictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] - (int16_t)previousPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0]);
                if (aidFuturePast < FLASH_TH) {
                    // proper action to be signalled
                    //SVT_LOG ("\nLight Flash in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                    ahd_running_avg[regionInPictureWidthIndex][regionInPictureHeightIndex] = (3 * ahd_running_avg[regionInPictureWidthIndex][regionInPictureHeightIndex] + ahd) / 4;
                }
                else {
                    // proper action to be signalled
                    //SVT_LOG ("\nLight Scene Change / fade detected in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                    ahd_running_avg[regionInPictureWidthIndex][regionInPictureHeightIndex] = (3 * ahd_running_avg[regionInPictureWidthIndex][regionInPictureHeightIndex] + ahd) / 4;
                }
            }
            else
                ahd_running_avg[regionInPictureWidthIndex][regionInPictureHeightIndex] = (3 * ahd_running_avg[regionInPictureWidthIndex][regionInPictureHeightIndex] + ahd) / 4;
            isAbruptChangeCount += isAbruptChange;
            isSceneChangeCount += is_scene_change;
        }
    }

    (void)windowWidthFuture;
    (void)isFlash;
    (void)isFade;

    if (isAbruptChangeCount >= regionCountThreshold)
        context_ptr->reset_running_avg = EB_TRUE;
    else
        context_ptr->reset_running_avg = EB_FALSE;
    if ((isSceneChangeCount >= regionCountThreshold) && ((!ParentPcsWindow[1]->fade_in_to_black) && (!ParentPcsWindow[1]->fade_out_from_black)))
        return(EB_TRUE);
    else
        return(EB_FALSE);
}

/***************************************************************************************************
* ReleasePrevPictureFromReorderQueue
***************************************************************************************************/
EbErrorType ReleasePrevPictureFromReorderQueue(
    EncodeContext                 *encode_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    PictureDecisionReorderEntry   *queuePreviousEntryPtr;
    int32_t                           previousEntryIndex;

    // Get the previous entry from the Picture Decision Reordering Queue (Entry N-1)
    // P.S. The previous entry in display order is needed for Scene Change Detection
    previousEntryIndex = (encode_context_ptr->picture_decision_reorder_queue_head_index == 0) ? PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1 : encode_context_ptr->picture_decision_reorder_queue_head_index - 1;
    queuePreviousEntryPtr = encode_context_ptr->picture_decision_reorder_queue[previousEntryIndex];

    // SB activity classification based on (0,0) SAD & picture activity derivation
    if (queuePreviousEntryPtr->parent_pcs_wrapper_ptr) {
        // Reset the Picture Decision Reordering Queue Entry
        // P.S. The reset of the Picture Decision Reordering Queue Entry could not be done before running the Scene Change Detector
        queuePreviousEntryPtr->picture_number += PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH;
        queuePreviousEntryPtr->parent_pcs_wrapper_ptr = (EbObjectWrapper *)EB_NULL;
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

    uint32_t MiniGopIndex;

    // Loop over all mini GOPs
    for (MiniGopIndex = 0; MiniGopIndex < MINI_GOP_MAX_COUNT; ++MiniGopIndex) {
        context_ptr->mini_gop_activity_array[MiniGopIndex] = (get_mini_gop_stats(MiniGopIndex)->hierarchical_levels == MIN_HIERARCHICAL_LEVEL) ?
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

    uint32_t    MiniGopIndex;

    context_ptr->total_number_of_mini_gops = 0;

    // Loop over all mini GOPs
    MiniGopIndex = 0;
    while (MiniGopIndex < MINI_GOP_MAX_COUNT) {
        // Only for a valid mini GOP
        if (get_mini_gop_stats(MiniGopIndex)->end_index < encode_context_ptr->pre_assignment_buffer_count && context_ptr->mini_gop_activity_array[MiniGopIndex] == EB_FALSE) {
            context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(MiniGopIndex)->start_index;
            context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(MiniGopIndex)->end_index;
            context_ptr->mini_gop_length[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(MiniGopIndex)->lenght;
            context_ptr->mini_gop_hierarchical_levels[context_ptr->total_number_of_mini_gops] = get_mini_gop_stats(MiniGopIndex)->hierarchical_levels;
            context_ptr->mini_gop_intra_count[context_ptr->total_number_of_mini_gops] = 0;
            context_ptr->mini_gop_idr_count[context_ptr->total_number_of_mini_gops] = 0;

            context_ptr->total_number_of_mini_gops++;
        }

        MiniGopIndex += context_ptr->mini_gop_activity_array[MiniGopIndex] ?
            1 :
            mini_gop_offset[get_mini_gop_stats(MiniGopIndex)->hierarchical_levels - MIN_HIERARCHICAL_LEVEL];
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
/***************************************************************************************************
* If a switch happens, then update the RPS of the base layer frame separating the 2 different prediction structures
* Clean up the reference queue dependant counts of the base layer frame separating the 2 different prediction structures
*
***************************************************************************************************/
EbErrorType update_base_layer_reference_queue_dependent_count(
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr,
    SequenceControlSet            *sequence_control_set_ptr,
    uint32_t                         MiniGopIndex) {
    if (!context_ptr || !encode_context_ptr || !sequence_control_set_ptr)
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

    PictureParentControlSet       *picture_control_set_ptr;

    // Get the 1st PCS mini GOP
    picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_start_index[MiniGopIndex]]->object_ptr;

    // Derive the temporal layer difference between the current mini GOP and the previous mini GOP
    picture_control_set_ptr->hierarchical_layers_diff = (uint8_t)(encode_context_ptr->previous_mini_gop_hierarchical_levels - picture_control_set_ptr->hierarchical_levels);

    // Set init_pred_struct_position_flag to TRUE if mini GOP switch
    picture_control_set_ptr->init_pred_struct_position_flag = (picture_control_set_ptr->hierarchical_layers_diff != 0) ?
        EB_TRUE :
        EB_FALSE;

    // If the current mini GOP is different than the previous mini GOP update then update the positive dependant counts of the reference entry separating the 2 mini GOPs
    if (picture_control_set_ptr->hierarchical_layers_diff != 0) {
        input_queue_index = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

        while (input_queue_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
            input_entry_ptr = encode_context_ptr->picture_decision_pa_reference_queue[input_queue_index];

            // Find the reference entry separating the 2 mini GOPs  (picture_control_set_ptr->picture_number is the POC of the first isput in the mini GOP)
            if (input_entry_ptr->picture_number == (picture_control_set_ptr->picture_number - 1)) {
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
                    picture_control_set_ptr->pred_structure,
                    sequence_control_set_ptr->reference_count,
                    picture_control_set_ptr->hierarchical_levels);            // Number of temporal layer in the current mini GOP

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
                    // currentInputPoc += (currentInputPoc < input_entry_ptr->pocNumber) ? (1 << sequence_control_set_ptr->bitsForPictureOrderCount) : 0;
                    dep_poc = POC_CIRCULAR_ADD(
                        input_entry_ptr->picture_number, // can't use a value that gets reset
                        input_entry_ptr->list0.list[dep_idx]/*,
                                                         sequence_control_set_ptr->bitsForPictureOrderCount*/);

                                                         // If Dependent POC is greater or equal to the IDR POC
                    if (dep_poc >= picture_control_set_ptr->picture_number && input_entry_ptr->list0.list[dep_idx]) {
                        input_entry_ptr->list0.list[dep_idx] = 0;

                        // Decrement the Reference's reference_count
                        --input_entry_ptr->dependent_count;
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
                    // currentInputPoc += (currentInputPoc < input_entry_ptr->pocNumber) ? (1 << sequence_control_set_ptr->bitsForPictureOrderCount) : 0;
                    dep_poc = POC_CIRCULAR_ADD(
                        input_entry_ptr->picture_number,
                        input_entry_ptr->list1.list[dep_idx]/*,
                                                         sequence_control_set_ptr->bitsForPictureOrderCount*/);

                    // If Dependent POC is greater or equal to the IDR POC
                    if ((dep_poc >= picture_control_set_ptr->picture_number) && input_entry_ptr->list1.list[dep_idx]) {
                        input_entry_ptr->list1.list[dep_idx] = 0;
                        // Decrement the Reference's reference_count
                        --input_entry_ptr->dependent_count;

                        CHECK_REPORT_ERROR(
                            (input_entry_ptr->dependent_count != ~0u),
                            encode_context_ptr->app_callback_ptr,
                            EB_ENC_PD_ERROR3);
                    }
                }
            }
            // Increment the input_queue_index Iterator
            input_queue_index = (input_queue_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : input_queue_index + 1;
        }
    }

    return return_error;
}

EbBool is_supposedly_4L_reference_frame(
    PictureDecisionContext        *context_ptr,
    uint32_t                         MiniGopIndex,
    uint32_t                        picture_index) {
    if ((context_ptr->mini_gop_hierarchical_levels[MiniGopIndex] == 4 && context_ptr->mini_gop_length[MiniGopIndex] == 16 && (picture_index == 7 || picture_index == 23)) ||    // supposedly a 4L reference frame for 5L prediction structure
        (context_ptr->mini_gop_hierarchical_levels[MiniGopIndex] == 5 && context_ptr->mini_gop_length[MiniGopIndex] == 32 && (picture_index == 7 || picture_index == 23))) { // supposedly a 4L reference frame for 6L prediction structure
        return(EB_TRUE);
    }
    else
        return(EB_FALSE);
}

/***************************************************************************************************
* Generates mini GOP RPSs
*
*
***************************************************************************************************/
EbErrorType GenerateMiniGopRps(
    PictureDecisionContext        *context_ptr,
    EncodeContext                 *encode_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    SequenceControlSet            *sequence_control_set_ptr;
    uint32_t                         mini_gop_index;
    PictureParentControlSet    *picture_control_set_ptr;
    uint32_t                         pictureIndex;

    // Loop over all mini GOPs
    for (mini_gop_index = 0; mini_gop_index < context_ptr->total_number_of_mini_gops; ++mini_gop_index) {
        // Loop over picture within the mini GOP
        for (pictureIndex = context_ptr->mini_gop_start_index[mini_gop_index]; pictureIndex <= context_ptr->mini_gop_end_index[mini_gop_index]; pictureIndex++) {
            picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pictureIndex]->object_ptr;
            sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
            picture_control_set_ptr->pred_structure = sequence_control_set_ptr->static_config.pred_structure;
            picture_control_set_ptr->hierarchical_levels = (uint8_t)context_ptr->mini_gop_hierarchical_levels[mini_gop_index];

            picture_control_set_ptr->pred_struct_ptr = get_prediction_structure(
                encode_context_ptr->prediction_structure_group_ptr,
                picture_control_set_ptr->pred_structure,
                sequence_control_set_ptr->reference_count,
                picture_control_set_ptr->hierarchical_levels);
        }
    }
    return return_error;
}

/******************************************************
* Derive Multi-Processes Settings for OQ
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/
EbErrorType signal_derivation_multi_processes_oq(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    FrameHeader *frm_hdr = &picture_control_set_ptr->frm_hdr;

    uint8_t sc_content_detected = picture_control_set_ptr->sc_content_detected;
    uint8_t enc_mode_hme = sequence_control_set_ptr->use_output_stat_file ? picture_control_set_ptr->snd_pass_enc_mode : picture_control_set_ptr->enc_mode;
    picture_control_set_ptr->enable_hme_flag = enable_hme_flag[picture_control_set_ptr->sc_content_detected][sequence_control_set_ptr->input_resolution][enc_mode_hme];

    picture_control_set_ptr->enable_hme_level0_flag = enable_hme_level0_flag[picture_control_set_ptr->sc_content_detected][sequence_control_set_ptr->input_resolution][enc_mode_hme];
    picture_control_set_ptr->enable_hme_level1_flag = enable_hme_level1_flag[picture_control_set_ptr->sc_content_detected][sequence_control_set_ptr->input_resolution][enc_mode_hme];
    picture_control_set_ptr->enable_hme_level2_flag = enable_hme_level2_flag[picture_control_set_ptr->sc_content_detected][sequence_control_set_ptr->input_resolution][enc_mode_hme];

    picture_control_set_ptr->tf_enable_hme_flag = tf_enable_hme_flag[picture_control_set_ptr->sc_content_detected][sequence_control_set_ptr->input_resolution][enc_mode_hme];
    picture_control_set_ptr->tf_enable_hme_level0_flag = tf_enable_hme_level0_flag[picture_control_set_ptr->sc_content_detected][sequence_control_set_ptr->input_resolution][enc_mode_hme];
    picture_control_set_ptr->tf_enable_hme_level1_flag = tf_enable_hme_level1_flag[picture_control_set_ptr->sc_content_detected][sequence_control_set_ptr->input_resolution][enc_mode_hme];
    picture_control_set_ptr->tf_enable_hme_level2_flag = tf_enable_hme_level2_flag[picture_control_set_ptr->sc_content_detected][sequence_control_set_ptr->input_resolution][enc_mode_hme];


        if (sc_content_detected)
#if PRESETS_TUNE
            if (picture_control_set_ptr->enc_mode <= ENC_M2)
                picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
            else if (picture_control_set_ptr->enc_mode <= ENC_M5)
                if (picture_control_set_ptr->slice_type == I_SLICE)
                    picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
                else
                    picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
            else
                picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;

#else
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
                picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
            else if (picture_control_set_ptr->enc_mode <= ENC_M3)
                if (picture_control_set_ptr->temporal_layer_index == 0)
                    picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
                else if (picture_control_set_ptr->is_used_as_reference_flag)
                    picture_control_set_ptr->pic_depth_mode = PIC_ALL_C_DEPTH_MODE;
                else
                    picture_control_set_ptr->pic_depth_mode = PIC_SQ_DEPTH_MODE;
            else if (picture_control_set_ptr->enc_mode <= ENC_M4)
                if (picture_control_set_ptr->slice_type == I_SLICE)
                    picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
                else
                    picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
            else
                picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
#endif
#if MULTI_PASS_PD
        else if (MR_MODE)
            picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
        else if (picture_control_set_ptr->enc_mode == ENC_M0)
            // Use a single-stage PD if I_SLICE
            picture_control_set_ptr->pic_depth_mode = (picture_control_set_ptr->slice_type == I_SLICE) ? PIC_ALL_DEPTH_MODE : PIC_MULTI_PASS_PD_MODE_1;
        else if (picture_control_set_ptr->enc_mode <= ENC_M2)
            // Use a single-stage PD if I_SLICE
            picture_control_set_ptr->pic_depth_mode = (picture_control_set_ptr->slice_type == I_SLICE) ? PIC_ALL_DEPTH_MODE : PIC_MULTI_PASS_PD_MODE_2;
#endif

#if PRESETS_TUNE
#if !MULTI_PASS_PD
        else if (picture_control_set_ptr->enc_mode <= ENC_M1)
            picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
        else if (picture_control_set_ptr->enc_mode <= ENC_M2)
            if (picture_control_set_ptr->slice_type == I_SLICE)
                picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
            else
                picture_control_set_ptr->pic_depth_mode = PIC_ALL_C_DEPTH_MODE;
        //Jing: TODO:
        //In M3, sb_sz may be 128x128, and init_sq_non4_block only works for 64x64 sb size
#endif
        else if (picture_control_set_ptr->enc_mode <= ENC_M6)
            picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
        else
            if (picture_control_set_ptr->slice_type == I_SLICE)
                picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
            else
                picture_control_set_ptr->pic_depth_mode = PIC_SB_SWITCH_DEPTH_MODE;
#else
// These are the default settings to use
#if !MULTI_PASS_PD
        else if (picture_control_set_ptr->enc_mode <= ENC_M2)
            picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
#endif

        //Jing: TODO:
        //In M3, sb_sz may be 128x128, and init_sq_non4_block only works for 64x64 sb size
        else if (picture_control_set_ptr->enc_mode <= ENC_M3)
            if (picture_control_set_ptr->slice_type == I_SLICE)
                picture_control_set_ptr->pic_depth_mode = PIC_ALL_C_DEPTH_MODE;
            else
                if (sequence_control_set_ptr->seq_header.sb_size == BLOCK_128X128)
                    picture_control_set_ptr->pic_depth_mode = PIC_SQ_DEPTH_MODE;
                else
                    picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;

        else if (picture_control_set_ptr->enc_mode <= ENC_M5)
            picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;

        else
            if (picture_control_set_ptr->slice_type == I_SLICE)
                picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
            else
                picture_control_set_ptr->pic_depth_mode = PIC_SB_SWITCH_DEPTH_MODE;
#endif

        if (picture_control_set_ptr->pic_depth_mode < PIC_SQ_DEPTH_MODE)
            assert(sequence_control_set_ptr->nsq_present == 1 && "use nsq_present 1");

        picture_control_set_ptr->max_number_of_pus_per_sb = (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) ? MAX_ME_PU_COUNT : SQUARE_PU_COUNT;
#if MDC_ADAPTIVE_LEVEL
        // Adaptive Ol  Level                    Settings
        // 0                                     OFF
        // 1                                     ON
        //NM : Please note that the open_loop_partitioning is operational only when
        // pic_depth_mode is set to PIC_ALL_DEPTH_MODE or PIC_ALL_C_DEPTH_MODE when
        // the motion information information for NSQ is generated.
        if (sequence_control_set_ptr->static_config.olpd_refinement == -1) { //auto mode; if not set by cfg
            if (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
                if (MR_MODE || sc_content_detected)
                    picture_control_set_ptr->enable_adaptive_ol_partitioning = 0;
                else if (picture_control_set_ptr->enc_mode <= ENC_M0)
                    picture_control_set_ptr->enable_adaptive_ol_partitioning = 1;
                else
                    picture_control_set_ptr->enable_adaptive_ol_partitioning = 0;
            }
            else
                picture_control_set_ptr->enable_adaptive_ol_partitioning = 0;
        }
        else {
            picture_control_set_ptr->enable_adaptive_ol_partitioning = sequence_control_set_ptr->static_config.olpd_refinement;
        }
#else

#if PREDICT_NSQ_SHAPE
    // Depth Level                           Settings
    // 0                                     pred only
    // 1                                     pred + 1
    // 2                                     pred + 2
    // 3                                     pred + 3
    // 4                                     pred - 1 + 1
    // 5                                     pred - 1 + 2
    // 6                                     pred - 1 + 3
    // 7                                     All
    if (MR_MODE || sc_content_detected)
        picture_control_set_ptr->mdc_depth_level = MAX_MDC_LEVEL;
    else if (picture_control_set_ptr->enc_mode == ENC_M0)
        picture_control_set_ptr->mdc_depth_level = (sequence_control_set_ptr->input_resolution == INPUT_SIZE_576p_RANGE_OR_LOWER) ? MAX_MDC_LEVEL : 6;
    else
        picture_control_set_ptr->mdc_depth_level = MAX_MDC_LEVEL; // Not tuned yet.

#endif
#endif

#if MULTI_PASS_PD && MDC_ADAPTIVE_LEVEL
    picture_control_set_ptr->enable_adaptive_ol_partitioning = 0;
#endif

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
            picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_FULL;
#if MULTI_PASS_PD
        else if (picture_control_set_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_0 ||
                 picture_control_set_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_1 ||
                 picture_control_set_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_2 ||
                 picture_control_set_ptr->pic_depth_mode == PIC_MULTI_PASS_PD_MODE_3 ){

            picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;
        }
#endif
        else if (sc_content_detected)
#if PRESETS_TUNE
#if M0_OPT
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
#else
            if (picture_control_set_ptr->enc_mode == ENC_M0)
                picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;
            else if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
                picture_control_set_ptr->nsq_search_level = (picture_control_set_ptr->is_used_as_reference_flag) ?
                NSQ_SEARCH_LEVEL6 : NSQ_SEARCH_LEVEL3;
            else if (picture_control_set_ptr->enc_mode <= ENC_M2)
                if (picture_control_set_ptr->is_used_as_reference_flag)
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL5;
                else
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL2;
            else
                picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_OFF;

#else
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
#if M0_OPT
                picture_control_set_ptr->nsq_search_level = (picture_control_set_ptr->is_used_as_reference_flag) ?
                NSQ_SEARCH_LEVEL6 : NSQ_SEARCH_LEVEL3;
#else
                picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;
#endif
            else if (picture_control_set_ptr->enc_mode <= ENC_M2)
                if (picture_control_set_ptr->temporal_layer_index == 0)
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;
                else if (picture_control_set_ptr->is_used_as_reference_flag)
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL4;
                else
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_OFF;
            else if (picture_control_set_ptr->enc_mode <= ENC_M3)
                if (picture_control_set_ptr->temporal_layer_index == 0)
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;
                else if (picture_control_set_ptr->is_used_as_reference_flag)
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL1;
                else
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_OFF;
            else
                    picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_OFF;
#endif
#if PREDICT_NSQ_SHAPE && !MDC_ADAPTIVE_LEVEL
        else if (picture_control_set_ptr->mdc_depth_level == (MAX_MDC_LEVEL - 1))
            picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL7;
#endif
#if PRESETS_TUNE
        else if (picture_control_set_ptr->enc_mode <= ENC_M0)
            picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;
        else if (picture_control_set_ptr->enc_mode <= ENC_M1)
            picture_control_set_ptr->nsq_search_level = (picture_control_set_ptr->is_used_as_reference_flag) ? NSQ_SEARCH_LEVEL6 : NSQ_SEARCH_LEVEL3;
        else if (picture_control_set_ptr->enc_mode <= ENC_M2)
            if (picture_control_set_ptr->is_used_as_reference_flag)
                picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL3;
            else
                picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL1;
        else
            picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_OFF;
#else
        else if (picture_control_set_ptr->enc_mode <= ENC_M1)
            picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL6;

        else if (picture_control_set_ptr->enc_mode <= ENC_M2)
            if (picture_control_set_ptr->is_used_as_reference_flag)
                picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL5;
            else
                picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_LEVEL3;
        else
            picture_control_set_ptr->nsq_search_level = NSQ_SEARCH_OFF;
#endif

    if (picture_control_set_ptr->nsq_search_level > NSQ_SEARCH_OFF)
        assert(sequence_control_set_ptr->nsq_present == 1 && "use nsq_present 1");

    switch (picture_control_set_ptr->nsq_search_level) {
    case NSQ_SEARCH_OFF:
        picture_control_set_ptr->nsq_max_shapes_md = 0;
        break;
    case NSQ_SEARCH_LEVEL1:
        picture_control_set_ptr->nsq_max_shapes_md = 1;
        break;
    case NSQ_SEARCH_LEVEL2:
        picture_control_set_ptr->nsq_max_shapes_md = 2;
        break;
    case NSQ_SEARCH_LEVEL3:
        picture_control_set_ptr->nsq_max_shapes_md = 3;
        break;
    case NSQ_SEARCH_LEVEL4:
        picture_control_set_ptr->nsq_max_shapes_md = 4;
        break;
    case NSQ_SEARCH_LEVEL5:
        picture_control_set_ptr->nsq_max_shapes_md = 5;
        break;
    case NSQ_SEARCH_LEVEL6:
        picture_control_set_ptr->nsq_max_shapes_md = 6;
        break;
#if PREDICT_NSQ_SHAPE
    case NSQ_SEARCH_LEVEL7:
        picture_control_set_ptr->nsq_max_shapes_md = 7;
        break;
#endif
    case NSQ_SEARCH_FULL:
        picture_control_set_ptr->nsq_max_shapes_md = 6;
        break;
    default:
        SVT_LOG("nsq_search_level is not supported\n");
        break;
    }

    if (picture_control_set_ptr->nsq_search_level == NSQ_SEARCH_OFF)
        if (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) picture_control_set_ptr->pic_depth_mode = PIC_SQ_DEPTH_MODE;
    if (picture_control_set_ptr->pic_depth_mode > PIC_SQ_DEPTH_MODE)
        assert(picture_control_set_ptr->nsq_search_level == NSQ_SEARCH_OFF);
#if !MULTI_PASS_PD
    // Interpolation search Level                     Settings
    // 0                                              OFF
    // 1                                              Interpolation search at inter-depth
    // 2                                              Interpolation search at full loop
    // 3                                              Chroma blind interpolation search at fast loop
    // 4                                              Interpolation search at fast loop

        if (MR_MODE)
            picture_control_set_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP;
        else if (sc_content_detected)
            picture_control_set_ptr->interpolation_search_level = IT_SEARCH_OFF;
        else if (picture_control_set_ptr->enc_mode <= ENC_M1)
            picture_control_set_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
        else if (picture_control_set_ptr->enc_mode <= ENC_M3)
            if (picture_control_set_ptr->is_used_as_reference_flag)
                picture_control_set_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
            else
                picture_control_set_ptr->interpolation_search_level = IT_SEARCH_OFF;
        else if (picture_control_set_ptr->enc_mode <= ENC_M7)
            if (picture_control_set_ptr->temporal_layer_index == 0)
                picture_control_set_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP_UV_BLIND;
            else
                picture_control_set_ptr->interpolation_search_level = IT_SEARCH_OFF;
        else
            picture_control_set_ptr->interpolation_search_level = IT_SEARCH_OFF;
#endif
    // Loop filter Level                            Settings
    // 0                                            OFF
    // 1                                            CU-BASED
    // 2                                            LIGHT FRAME-BASED
    // 3                                            FULL FRAME-BASED

    //for now only I frames are allowed to use sc tools.
    //TODO: we can force all frames in GOP with the same detection status of leading I frame.
    if (picture_control_set_ptr->slice_type == I_SLICE) {
        frm_hdr->allow_screen_content_tools = picture_control_set_ptr->sc_content_detected;
        if (picture_control_set_ptr->enc_mode <= ENC_M5)
            frm_hdr->allow_intrabc =  picture_control_set_ptr->sc_content_detected;
        else
            frm_hdr->allow_intrabc =  0;

        //IBC Modes:   0:Slow   1:Fast   2:Faster
#if PRESETS_TUNE
        if (picture_control_set_ptr->enc_mode <= ENC_M5)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M2)
#endif
            picture_control_set_ptr->ibc_mode = 0;
        else
            picture_control_set_ptr->ibc_mode = 1;
    }
    else {
#if PAL_SUP
        //this will enable sc tools for P frames. hence change bitstream even if palette mode is OFF
        frm_hdr->allow_screen_content_tools = picture_control_set_ptr->sc_content_detected;
#else
        frm_hdr->allow_screen_content_tools = 0;
#endif
        frm_hdr->allow_intrabc = 0;
    }

#if PAL_SUP

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
        if (sequence_control_set_ptr->static_config.enable_palette == -1)//auto mode; if not set by cfg
            picture_control_set_ptr->palette_mode =
            (sequence_control_set_ptr->static_config.encoder_bit_depth == EB_8BIT ||
            (sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT && sequence_control_set_ptr->static_config.enable_hbd_mode_decision == 0)) &&
            picture_control_set_ptr->enc_mode == ENC_M0 ? 6 : 0;
        else
            picture_control_set_ptr->palette_mode = sequence_control_set_ptr->static_config.enable_palette;
    else
        picture_control_set_ptr->palette_mode = 0;

    assert(picture_control_set_ptr->palette_mode<7);
#endif
    if (!picture_control_set_ptr->sequence_control_set_ptr->static_config.disable_dlf_flag && frm_hdr->allow_intrabc == 0) {
    if (sc_content_detected)
#if M0_OPT
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#else
        if (picture_control_set_ptr->enc_mode == ENC_M0)
            picture_control_set_ptr->loop_filter_mode = 3;
        else if (picture_control_set_ptr->enc_mode == ENC_M1)
#endif
            picture_control_set_ptr->loop_filter_mode = picture_control_set_ptr->is_used_as_reference_flag ? 3 : 0;
        else
            picture_control_set_ptr->loop_filter_mode = 0;
    else
#if PRESETS_OPT
        if (picture_control_set_ptr->enc_mode <= ENC_M7)
            picture_control_set_ptr->loop_filter_mode = 3;
#else
#if PRESETS_TUNE
        if (picture_control_set_ptr->enc_mode <= ENC_M2)
#else
        if (picture_control_set_ptr->enc_mode == ENC_M0)
#endif
            picture_control_set_ptr->loop_filter_mode = 3;
        else if (picture_control_set_ptr->enc_mode <= ENC_M5)
            picture_control_set_ptr->loop_filter_mode = picture_control_set_ptr->is_used_as_reference_flag ? 3 : 0;
#endif
        else
            picture_control_set_ptr->loop_filter_mode = picture_control_set_ptr->is_used_as_reference_flag ? 1 : 0;
    }
    else
        picture_control_set_ptr->loop_filter_mode = 0;
    // CDEF Level                                   Settings
    // 0                                            OFF
    // 1                                            1 step refinement
    // 2                                            4 step refinement
    // 3                                            8 step refinement
    // 4                                            16 step refinement
    // 5                                            64 step refinement
    if (sequence_control_set_ptr->seq_header.enable_cdef && frm_hdr->allow_intrabc == 0) {
        if (sc_content_detected)
#if PRESETS_TUNE
            if (picture_control_set_ptr->enc_mode <= ENC_M5)
#else
            if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
                picture_control_set_ptr->cdef_filter_mode = 4;
            else
                picture_control_set_ptr->cdef_filter_mode = 0;
        else
        if (picture_control_set_ptr->enc_mode <= ENC_M7)
#if M0_OPT
#if PRESETS_TUNE
#if PRESETS_OPT
            picture_control_set_ptr->cdef_filter_mode = 5;
#else
            picture_control_set_ptr->cdef_filter_mode = (picture_control_set_ptr->enc_mode <= ENC_M1) ? 5 : 4;
#endif
#else
            picture_control_set_ptr->cdef_filter_mode = (picture_control_set_ptr->enc_mode <= ENC_M0) ? 5 : 4;
#endif
#else
            picture_control_set_ptr->cdef_filter_mode = 4;
#endif
        else
            picture_control_set_ptr->cdef_filter_mode = 2;
    }
    else
        picture_control_set_ptr->cdef_filter_mode = 0;

    // SG Level                                    Settings
    // 0                                            OFF
    // 1                                            0 step refinement
    // 2                                            1 step refinement
    // 3                                            4 step refinement
    // 4                                            16 step refinement

    Av1Common* cm = picture_control_set_ptr->av1_cm;
    if (sc_content_detected)
#if PRESETS_TUNE
        if (picture_control_set_ptr->enc_mode <= ENC_M5)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            cm->sg_filter_mode = 4;
        else
            cm->sg_filter_mode = 0;
    else
#if PRESETS_TUNE
        if (picture_control_set_ptr->enc_mode <= ENC_M2)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M4)
#endif
        cm->sg_filter_mode = 4;
    else if (picture_control_set_ptr->enc_mode <= ENC_M6)
        cm->sg_filter_mode = 3;
    else
        cm->sg_filter_mode = 1;

    // WN Level                                     Settings
    // 0                                            OFF
    // 1                                            3-Tap luma/ 3-Tap chroma
    // 2                                            5-Tap luma/ 5-Tap chroma
    // 3                                            7-Tap luma/ 5-Tap chroma

    if (sc_content_detected)
#if PRESETS_TUNE
        if (picture_control_set_ptr->enc_mode <= ENC_M5)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            cm->wn_filter_mode = 3;
        else
            cm->wn_filter_mode = 0;
    else

    if (picture_control_set_ptr->enc_mode <= ENC_M5)
        cm->wn_filter_mode = 3;
    else if (picture_control_set_ptr->enc_mode <= ENC_M7)
        cm->wn_filter_mode = 2;
    else
        cm->wn_filter_mode = 0;
#if !MULTI_PASS_PD
    // Tx_search Level                                Settings
    // 0                                              OFF
    // 1                                              Tx search at encdec
    // 2                                              Tx search at inter-depth
    // 3                                              Tx search at full loop
    if (sc_content_detected)
        if (picture_control_set_ptr->enc_mode <= ENC_M6)
            picture_control_set_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            if (picture_control_set_ptr->is_used_as_reference_flag)
                picture_control_set_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
            else
                picture_control_set_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    else
    if (picture_control_set_ptr->enc_mode <= ENC_M4)
        picture_control_set_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
    else if (picture_control_set_ptr->enc_mode <= ENC_M7) {
        if (picture_control_set_ptr->temporal_layer_index == 0)
            picture_control_set_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
        else
            picture_control_set_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    }
    else
        picture_control_set_ptr->tx_search_level = TX_SEARCH_ENC_DEC;

    // Set tx search skip weights (MAX_MODE_COST: no skipping; 0: always skipping)
    if (picture_control_set_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
        picture_control_set_ptr->tx_weight = MAX_MODE_COST;
    else if (!MR_MODE && picture_control_set_ptr->enc_mode <= ENC_M1)
        picture_control_set_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
    else if (!MR_MODE){
        if (picture_control_set_ptr->is_used_as_reference_flag)
            picture_control_set_ptr->tx_weight = FC_SKIP_TX_SR_TH025;
        else
            picture_control_set_ptr->tx_weight = FC_SKIP_TX_SR_TH010;
    }

    // Set tx search reduced set falg (0: full tx set; 1: reduced tx set; 1: two tx))
    if (sc_content_detected)
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
            picture_control_set_ptr->tx_search_reduced_set = 0;
        else if (picture_control_set_ptr->enc_mode <= ENC_M6)
            if (picture_control_set_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
                picture_control_set_ptr->tx_search_reduced_set = 0;
            else
                picture_control_set_ptr->tx_search_reduced_set = 1;
        else if (picture_control_set_ptr->enc_mode <= ENC_M7)
            picture_control_set_ptr->tx_search_reduced_set = 1;
        else
            picture_control_set_ptr->tx_search_reduced_set = 2;
    else

    if (picture_control_set_ptr->tx_search_level == TX_SEARCH_ENC_DEC)
        picture_control_set_ptr->tx_search_reduced_set = 0;
    else if (picture_control_set_ptr->enc_mode <= ENC_M1)
        picture_control_set_ptr->tx_search_reduced_set = 0;
    else if (picture_control_set_ptr->enc_mode <= ENC_M3)
        if (picture_control_set_ptr->is_used_as_reference_flag)
            picture_control_set_ptr->tx_search_reduced_set = 0;
        else
            picture_control_set_ptr->tx_search_reduced_set = 1;
    else
        picture_control_set_ptr->tx_search_reduced_set = 1;
#endif
    // Intra prediction modes                       Settings
    // 0                                            FULL
    // 1                                            LIGHT per block : disable_z2_prediction && disable_angle_refinement  for 64/32/4
    // 2                                            OFF per block : disable_angle_prediction for 64/32/4
    // 3                                            OFF : disable_angle_prediction
    // 4                                            OIS based Intra
    // 5                                            Light OIS based Intra

    if (picture_control_set_ptr->slice_type == I_SLICE)
    if (sc_content_detected)
        if (picture_control_set_ptr->enc_mode <= ENC_M6)
            picture_control_set_ptr->intra_pred_mode = 0;
        else
            picture_control_set_ptr->intra_pred_mode = 4;
    else
#if PRESETS_OPT
        if (picture_control_set_ptr->enc_mode <= ENC_M7)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M6)
#endif
            picture_control_set_ptr->intra_pred_mode = 0;
        else
            picture_control_set_ptr->intra_pred_mode = 4;
    else {
    if (sc_content_detected)
#if M0_OPT

        if (picture_control_set_ptr->enc_mode <= ENC_M2)
#else
        if (picture_control_set_ptr->enc_mode == ENC_M0)
            picture_control_set_ptr->intra_pred_mode = 0;
        else if (picture_control_set_ptr->enc_mode <= ENC_M2)
#endif
            if (picture_control_set_ptr->temporal_layer_index == 0)
                picture_control_set_ptr->intra_pred_mode = 1;
            else
                picture_control_set_ptr->intra_pred_mode = 2;
        else if (picture_control_set_ptr->enc_mode <= ENC_M6)
            if (picture_control_set_ptr->temporal_layer_index == 0)
                picture_control_set_ptr->intra_pred_mode = 2;
            else
                picture_control_set_ptr->intra_pred_mode = 3;
        else
            picture_control_set_ptr->intra_pred_mode = 4;
    else
#if PRESETS_TUNE
        if ((picture_control_set_ptr->enc_mode <= ENC_M1) || (picture_control_set_ptr->enc_mode <= ENC_M2 && picture_control_set_ptr->temporal_layer_index == 0))
            picture_control_set_ptr->intra_pred_mode = 0;
#else
        if (picture_control_set_ptr->enc_mode == ENC_M0)
            picture_control_set_ptr->intra_pred_mode = 0;
        else if (picture_control_set_ptr->enc_mode  <= ENC_M1)
            if (picture_control_set_ptr->temporal_layer_index == 0)
                picture_control_set_ptr->intra_pred_mode = 1;
            else
                picture_control_set_ptr->intra_pred_mode = 2;
#endif
#if PRESETS_OPT
        else if (picture_control_set_ptr->enc_mode <= ENC_M7)
#else
        else if(picture_control_set_ptr->enc_mode <= ENC_M6)
#endif
            if (picture_control_set_ptr->temporal_layer_index == 0)
                picture_control_set_ptr->intra_pred_mode = 1;
            else
                picture_control_set_ptr->intra_pred_mode = 3;
        else
            picture_control_set_ptr->intra_pred_mode = 4;
    }

    if (MR_MODE)
        picture_control_set_ptr->intra_pred_mode = 0;

    // Skip sub blk based on neighbors depth        Settings
    // 0                                            OFF
    // 1                                            ON
    picture_control_set_ptr->skip_sub_blks =   0;

        if (picture_control_set_ptr->sc_content_detected)
            picture_control_set_ptr->cu8x8_mode = (picture_control_set_ptr->temporal_layer_index > 0) ?
            CU_8x8_MODE_1 :
            CU_8x8_MODE_0;
        else
        if (picture_control_set_ptr->enc_mode <= ENC_M8)
            picture_control_set_ptr->cu8x8_mode = CU_8x8_MODE_0;
        else
            picture_control_set_ptr->cu8x8_mode = (picture_control_set_ptr->temporal_layer_index > 0) ?
            CU_8x8_MODE_1 :
            CU_8x8_MODE_0;

        // Set atb mode      Settings
        // 0                 OFF: no transform partitioning
        // 1                 ON for INTRA blocks
        if (picture_control_set_ptr->enc_mode <= ENC_M2)
            picture_control_set_ptr->atb_mode = (MR_MODE || picture_control_set_ptr->temporal_layer_index == 0) ? 1 : 0;
        else
            picture_control_set_ptr->atb_mode = 0;

        // Set skip atb                          Settings
        // 0                                     OFF
        // 1                                     ON
        // Phoenix: Active only when atb compound is on
        if (picture_control_set_ptr->enc_mode <= ENC_M7)
            picture_control_set_ptr->coeff_based_skip_atb = 0;
        else
            picture_control_set_ptr->coeff_based_skip_atb = 1;
        // Set Wedge mode      Settings
        // 0                 FULL: Full search
        // 1                 Fast: If two predictors are very similar, skip wedge compound mode search
        // 2                 Fast: estimate Wedge sign
        // 3                 Fast: Mode 1 & Mode 2

        picture_control_set_ptr->wedge_mode = 0;

        // inter intra pred                      Settings
        // 0                                     OFF
        // 1                                     ON
#if II_COMP_FLAG
        if (sequence_control_set_ptr->static_config.inter_intra_compound == DEFAULT)
            picture_control_set_ptr->enable_inter_intra = picture_control_set_ptr->slice_type != I_SLICE ? sequence_control_set_ptr->seq_header.enable_interintra_compound : 0;
        else
            picture_control_set_ptr->enable_inter_intra = sequence_control_set_ptr->static_config.inter_intra_compound;
#endif

        // Set compound mode      Settings
        // 0                 OFF: No compond mode search : AVG only
        // 1                 ON: compond mode search: AVG/DIST/DIFF
        // 2                 ON: AVG/DIST/DIFF/WEDGE
        if (sequence_control_set_ptr->static_config.compound_level == DEFAULT) {
            if (sequence_control_set_ptr->compound_mode)
#if PRESETS_TUNE
                if (picture_control_set_ptr->sc_content_detected)
                    picture_control_set_ptr->compound_mode = (picture_control_set_ptr->enc_mode <= ENC_M0) ? 2 : 0;
                else
                    picture_control_set_ptr->compound_mode = picture_control_set_ptr->enc_mode <= ENC_M1 ? 2 : 1;
#else
                picture_control_set_ptr->compound_mode = picture_control_set_ptr->sc_content_detected ? 0 :
                picture_control_set_ptr->enc_mode <= ENC_M1 ? 2 : 1;
#endif
            else
                picture_control_set_ptr->compound_mode = 0;
#if !MULTI_PASS_PD
            // set compound_types_to_try
            if (picture_control_set_ptr->compound_mode)
                picture_control_set_ptr->compound_types_to_try = picture_control_set_ptr->compound_mode == 1 ? MD_COMP_DIFF0 : MD_COMP_WEDGE;
            else
                picture_control_set_ptr->compound_types_to_try = MD_COMP_AVG;
#endif
        }
        else
            picture_control_set_ptr->compound_mode = sequence_control_set_ptr->static_config.compound_level;

        // Set frame end cdf update mode      Settings
        // 0                                     OFF
        // 1                                     ON
        if (sequence_control_set_ptr->static_config.frame_end_cdf_update == DEFAULT)
#if PRESETS_TUNE
            picture_control_set_ptr->frame_end_cdf_update_mode = 1;
#else
            if (picture_control_set_ptr->enc_mode == ENC_M0)
                picture_control_set_ptr->frame_end_cdf_update_mode = 1;
            else
                picture_control_set_ptr->frame_end_cdf_update_mode = 0;
#endif
        else
            picture_control_set_ptr->frame_end_cdf_update_mode = sequence_control_set_ptr->static_config.frame_end_cdf_update;

        if (sequence_control_set_ptr->static_config.prune_unipred_me == DEFAULT)
#if M0_OPT
            if (picture_control_set_ptr->sc_content_detected || picture_control_set_ptr->enc_mode >= ENC_M4)
#else
            if (picture_control_set_ptr->sc_content_detected || picture_control_set_ptr->enc_mode == ENC_M0 || picture_control_set_ptr->enc_mode >= ENC_M4)
#endif
                picture_control_set_ptr->prune_unipred_at_me = 0;
            else
                picture_control_set_ptr->prune_unipred_at_me = 1;
        else
            picture_control_set_ptr->prune_unipred_at_me = sequence_control_set_ptr->static_config.prune_unipred_me;

        //CHKN: Temporal MVP should be disabled for pictures beloning to 4L MiniGop preceeded by 5L miniGOP. in this case the RPS is wrong(known issue). check RPS construction for more info.
        if (picture_control_set_ptr->slice_type == I_SLICE)
            picture_control_set_ptr->frm_hdr.use_ref_frame_mvs = 0;
        else
            picture_control_set_ptr->frm_hdr.use_ref_frame_mvs = sequence_control_set_ptr->mfmv_enabled;

#if GM_OPT
        // Global motion level                        Settings
        // GM_FULL                                    Exhaustive search mode.
        // GM_DOWN                                    Downsampled resolution with a downsampling factor of 2 in each dimension
        // GM_TRAN_ONLY                               Translation only using ME MV.
        picture_control_set_ptr->gm_level = GM_FULL;
#endif
        //Exit TX size search when all coefficients are zero
        // 0: OFF
        // 1: ON
        picture_control_set_ptr->tx_size_early_exit = 1;
    return return_error;
}

int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
//set the ref frame types used for this picture,
static void set_all_ref_frame_type(SequenceControlSet *sequence_control_set_ptr, PictureParentControlSet  *parent_pcs_ptr, MvReferenceFrame ref_frame_arr[], uint8_t* tot_ref_frames)
{
    MvReferenceFrame rf[2];
    *tot_ref_frames = 0;

    //SVT_LOG("POC %i  totRef L0:%i   totRef L1: %i\n", parent_pcs_ptr->picture_number, parent_pcs_ptr->ref_list0_count, parent_pcs_ptr->ref_list1_count);

    //single ref - List0
    for (uint8_t ref_idx0 = 0; ref_idx0 < parent_pcs_ptr->ref_list0_count; ++ref_idx0) {
        rf[0] = svt_get_ref_frame_type(REF_LIST_0, ref_idx0);
        ref_frame_arr[(*tot_ref_frames)++] = rf[0];
    }

    //single ref - List1
    for (uint8_t ref_idx1 = 0; ref_idx1 < parent_pcs_ptr->ref_list1_count; ++ref_idx1) {
        rf[1] = svt_get_ref_frame_type(REF_LIST_1, ref_idx1);
        ref_frame_arr[(*tot_ref_frames)++] = rf[1];
    }

    //compound Bi-Dir
    for (uint8_t ref_idx0 = 0; ref_idx0 < parent_pcs_ptr->ref_list0_count; ++ref_idx0) {
        for (uint8_t ref_idx1 = 0; ref_idx1 < parent_pcs_ptr->ref_list1_count; ++ref_idx1) {
            rf[0] = svt_get_ref_frame_type(REF_LIST_0, ref_idx0);
            rf[1] = svt_get_ref_frame_type(REF_LIST_1, ref_idx1);
            ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
        }
    }

    if (sequence_control_set_ptr->mrp_mode == 0 && parent_pcs_ptr->slice_type == B_SLICE)
    {

        //compound Uni-Dir
        if (parent_pcs_ptr->ref_list0_count > 1) {
            rf[0] = LAST_FRAME;
            rf[1] = LAST2_FRAME;
            ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
            if (parent_pcs_ptr->ref_list0_count > 2) {
                rf[1] = LAST3_FRAME;
                ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
                if (parent_pcs_ptr->ref_list0_count > 3) {
                    rf[1] = GOLDEN_FRAME;
                    ref_frame_arr[(*tot_ref_frames)++] = av1_ref_frame_type(rf);
                }
            }
        }
        if (parent_pcs_ptr->ref_list1_count > 2) {
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
// 1)Low delay P, 2)Low delay B and 3)I frames of RA
// For B frames of RA, need to set it manually based on picture_index
static EbBool set_frame_display_params(
        PictureParentControlSet       *picture_control_set_ptr,
        PictureDecisionContext        *context_ptr,
        uint32_t                       mini_gop_index)
{
    Av1RpsNode *av1_rps = &picture_control_set_ptr->av1_ref_signal;
    FrameHeader *frm_hdr = &picture_control_set_ptr->frm_hdr;

    if (picture_control_set_ptr->pred_struct_ptr->pred_type == EB_PRED_LOW_DELAY_P || picture_control_set_ptr->is_overlay) {
        //P frames
        av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_poc_array[BWD] = av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[LAST];

        frm_hdr->show_frame = EB_TRUE;
        picture_control_set_ptr->has_show_existing = EB_FALSE;
    } else if (picture_control_set_ptr->pred_struct_ptr->pred_type == EB_PRED_LOW_DELAY_B) {
        av1_rps->ref_dpb_index[BWD] = av1_rps->ref_dpb_index[LAST];
        av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[LAST2];
        av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[LAST3];
        av1_rps->ref_poc_array[BWD] = av1_rps->ref_poc_array[LAST];
        av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[LAST2];
        av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[LAST3];

        frm_hdr->show_frame = EB_TRUE;
        picture_control_set_ptr->has_show_existing = EB_FALSE;
    } else if (picture_control_set_ptr->pred_struct_ptr->pred_type == EB_PRED_RANDOM_ACCESS) {
        //Decide on Show Mecanism
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            //3 cases for I slice:  1:Key Frame treated above.  2: broken MiniGop due to sc or intra refresh  3: complete miniGop due to sc or intra refresh
            if (context_ptr->mini_gop_length[mini_gop_index] < picture_control_set_ptr->pred_struct_ptr->pred_struct_period) {
                //Scene Change that breaks the mini gop and switch to LDP (if I scene change happens to be aligned with a complete miniGop, then we do not break the pred structure)
                frm_hdr->show_frame = EB_TRUE;
                picture_control_set_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_FALSE;
                picture_control_set_ptr->has_show_existing = EB_FALSE;
            }
        } else {
            if (context_ptr->mini_gop_length[0] != picture_control_set_ptr->pred_struct_ptr->pred_struct_period) {
                SVT_LOG("Error in GOP indexing3\n");
            }
            // Handle B frame of Random Access out
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
static void  Av1GenerateRpsInfo(
    PictureParentControlSet       *picture_control_set_ptr,
    EncodeContext                 *encode_context_ptr,
    PictureDecisionContext        *context_ptr,
    uint32_t                       picture_index,
    uint32_t                       mini_gop_index
)
{
    (void)encode_context_ptr;
    Av1RpsNode *av1_rps = &picture_control_set_ptr->av1_ref_signal;
    FrameHeader *frm_hdr = &picture_control_set_ptr->frm_hdr;

    SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    PredictionStructureEntry *pred_position_ptr = picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[picture_control_set_ptr->pred_struct_index];
    //Set frame type
    if (picture_control_set_ptr->slice_type == I_SLICE)
        frm_hdr->frame_type = picture_control_set_ptr->idr_flag ? KEY_FRAME : INTRA_ONLY_FRAME;
    else
        frm_hdr->frame_type = INTER_FRAME;

    picture_control_set_ptr->intra_only = picture_control_set_ptr->slice_type == I_SLICE ? 1 : 0;

    if (picture_control_set_ptr->hierarchical_levels == 0) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(picture_control_set_ptr, context_ptr);
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
        av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, flat_pred_struct[gop_i].ref_list0[0]);
        av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, flat_pred_struct[gop_i].ref_list0[1]);
        av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
        av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

        av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, flat_pred_struct[gop_i].ref_list1[0]);
        av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, flat_pred_struct[gop_i].ref_list1[1]);
        av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

        prune_refs(pred_position_ptr, av1_rps);
        av1_rps->refresh_frame_mask = 1 << context_ptr->lay0_toggle;

        // Flat mode, output all frames
        set_frame_display_params(picture_control_set_ptr, context_ptr, mini_gop_index);
        frm_hdr->show_frame = EB_TRUE;
        picture_control_set_ptr->has_show_existing = EB_FALSE;
        context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
    } else if (picture_control_set_ptr->hierarchical_levels == 1) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(picture_control_set_ptr, context_ptr);
            return;
        }

        const uint8_t  base0_idx = context_ptr->lay0_toggle == 0 ? 1 : context_ptr->lay0_toggle == 1 ? 2 : 0;
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 0 : 1;
        const uint8_t  base2_idx = context_ptr->lay0_toggle == 0 ? 0 : context_ptr->lay0_toggle == 1 ? 1 : 2;

        switch (picture_control_set_ptr->temporal_layer_index) {
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
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, two_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 0;

            break;

        default:
            SVT_LOG("Error: unexpected picture mini Gop number\n");
            break;
        }

        prune_refs(pred_position_ptr, av1_rps);

        if (!set_frame_display_params(picture_control_set_ptr, context_ptr, mini_gop_index)) {
            if (picture_control_set_ptr->is_used_as_reference_flag && picture_index != 0) {
                frm_hdr->show_frame = EB_FALSE;
                picture_control_set_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                picture_control_set_ptr->has_show_existing = EB_TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = base2_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", picture_control_set_ptr->hierarchical_levels);
            }
        }

        // Toggle layer0 and layer1
        if ((picture_index == context_ptr->mini_gop_end_index[mini_gop_index] &&
                    !picture_control_set_ptr->is_overlay &&
                    sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) ||
                ((sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                 (picture_control_set_ptr->temporal_layer_index == 0))) {
            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
        }
    } else if (picture_control_set_ptr->hierarchical_levels == 2) {
        uint8_t gop_i;
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(picture_control_set_ptr, context_ptr);
            return;
        }
        // Slot 0,1,2 for base layer, 3,4 for layer 1, 5 for layer2
        const uint8_t  base0_idx = context_ptr->lay0_toggle == 0 ? 1 : context_ptr->lay0_toggle == 1 ? 2 : 0;
        const uint8_t  base1_idx = context_ptr->lay0_toggle == 0 ? 2 : context_ptr->lay0_toggle == 1 ? 0 : 1;
        const uint8_t  base2_idx = context_ptr->lay0_toggle == 0 ? 0 : context_ptr->lay0_toggle == 1 ? 1 : 2;

        const uint8_t  lay1_0_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 1 : LAY1_OFF + 0;
        const uint8_t  lay1_1_idx = context_ptr->lay1_toggle == 0 ? LAY1_OFF + 0 : LAY1_OFF + 1;
        const uint8_t  lay2_idx = LAY2_OFF;

        switch (picture_control_set_ptr->temporal_layer_index) {
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
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];

            av1_rps->refresh_frame_mask = 1 << (LAY1_OFF + context_ptr->lay1_toggle);

            break;

        case 2:
            // update RPS for the overlay frame.
            if (picture_control_set_ptr->is_overlay) {
                //{ 0, 0, 0, 0}         // GOP Index 1 - Ref List 0
                //{ 0, 0, 0, 0 }       // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST]  = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD]  = base1_idx;
                av1_rps->ref_dpb_index[BWD]  = base1_idx;
                av1_rps->ref_dpb_index[ALT2]  = base1_idx;
                av1_rps->ref_dpb_index[ALT] = base1_idx;

                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
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
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                    av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                    av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                    av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                    av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                    av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                    av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, three_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
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

        if (!set_frame_display_params(picture_control_set_ptr, context_ptr, mini_gop_index)) {
            if (picture_control_set_ptr->is_used_as_reference_flag && picture_index != 0) {
                frm_hdr->show_frame = EB_FALSE;
                picture_control_set_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                picture_control_set_ptr->has_show_existing = EB_TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = lay1_1_idx;
                else if (picture_index == 2)
                    frm_hdr->show_existing_frame = base2_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", picture_control_set_ptr->hierarchical_levels);
            }
        }

        if ((picture_index == (context_ptr->mini_gop_end_index[mini_gop_index]) &&
                    !picture_control_set_ptr->is_overlay &&
                    sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) ||
                ((sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                 sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                 (picture_control_set_ptr->temporal_layer_index == 0))) {
            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
        }
    } else if (picture_control_set_ptr->hierarchical_levels == 3) {

        uint8_t gop_i;
        EbBool is_trailing_frames = EB_FALSE;

        if (frm_hdr->frame_type == KEY_FRAME)
        {
            set_key_frame_rps(picture_control_set_ptr, context_ptr);
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
        if (picture_control_set_ptr->pred_struct_ptr->pred_type == EB_PRED_LOW_DELAY_P &&
                context_ptr->mini_gop_length[mini_gop_index] < 8 &&
                picture_control_set_ptr->sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {
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

        switch (picture_control_set_ptr->temporal_layer_index) {
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
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
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
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }

            av1_rps->refresh_frame_mask = 1 << (LAY2_OFF + context_ptr->lay2_toggle);
            //toggle 3->4
            if (!is_trailing_frames ||
                    (context_ptr->mini_gop_length[mini_gop_index] >= 7 && is_trailing_frames &&
                     picture_control_set_ptr->sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)) {
                //For trailing frames, Only toggle it if we are sure we have 2 layer 2 frames in trailing frames
                context_ptr->lay2_toggle = 1 - context_ptr->lay2_toggle;
            }

            break;

        case 3:
            // update RPS for the overlay frame.
            if (picture_control_set_ptr->is_overlay) {
                //{ 0, 0, 0, 0}         // GOP Index 1 - Ref List 0
                //{ 0, 0, 0, 0 }       // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST]  = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD]  = base1_idx;

                av1_rps->ref_dpb_index[BWD]  = base1_idx;
                av1_rps->ref_dpb_index[ALT2]  = base1_idx;
                av1_rps->ref_dpb_index[ALT] = base1_idx;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);
                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, four_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
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

        if (!set_frame_display_params(picture_control_set_ptr, context_ptr, mini_gop_index)) {
            if (picture_control_set_ptr->is_used_as_reference_flag && picture_index != 0)
            {
                frm_hdr->show_frame = EB_FALSE;
                picture_control_set_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                picture_control_set_ptr->has_show_existing = EB_TRUE;

                if (picture_index == 0)
                    frm_hdr->show_existing_frame = lay2_1_idx;
                else if (picture_index == 2)
                    frm_hdr->show_existing_frame = lay1_1_idx;
                else if (picture_index == 4)
                    frm_hdr->show_existing_frame = lay2_1_idx;
                else if (picture_index == 6)
                    frm_hdr->show_existing_frame = base2_idx;
                else
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", picture_control_set_ptr->hierarchical_levels);
            }
        }

        //last pic in MiniGop: Base layer toggling
        //mini GOP toggling since last Key Frame.
        //a regular I keeps the toggling process and does not reset the toggle.  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
        //whoever needs a miniGOP Level toggling, this is the time
        if (((picture_index == (context_ptr->mini_gop_end_index[mini_gop_index] % 8) &&
                        !picture_control_set_ptr->is_overlay &&
                        sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)) ||
                ((sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                 (picture_control_set_ptr->temporal_layer_index == 0))) {

            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;


            //layer2 toggle for low delay B/P case
            if ((sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                    sequence_control_set_ptr->intra_period_length != -1 &&
                    picture_control_set_ptr->cra_flag) {
                int32_t trailing_count = (sequence_control_set_ptr->intra_period_length % 8);
                if (trailing_count >=2 && trailing_count <= 5) {
                    context_ptr->lay2_toggle = 1 - context_ptr->lay2_toggle;
                }
            }
        }
    } else if (picture_control_set_ptr->hierarchical_levels == 4) {
        uint8_t gop_i;
        //Av1RpsNode_t *av1_rps = &picture_control_set_ptr->av1RefSignal2;

        //Reset miniGop Toggling. The first miniGop after a KEY frame has toggle=0
        if (frm_hdr->frame_type == KEY_FRAME) {
            set_key_frame_rps(picture_control_set_ptr, context_ptr);
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
        switch (picture_control_set_ptr->temporal_layer_index) {
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
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = av1_rps->ref_poc_array[LAST];
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
            av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
            av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];

            gop_i = 8;
            av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
            av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
            av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
            av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

            av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
            av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
            av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else if (picture_index == 11) {
                //{ 4, 8, 12, 0},       // GOP Index 12 - Ref List 0
                //{ -4,  0, 0,  0 }     // GOP Index 12 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = lay1_1_idx;
                av1_rps->ref_dpb_index[LAST2] = lay2_1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = av1_rps->ref_dpb_index[LAST];

                av1_rps->ref_dpb_index[BWD] = base2_idx;
                av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 12;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = av1_rps->ref_poc_array[LAST];

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 14;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
                av1_rps->ref_poc_array[ALT] = av1_rps->ref_poc_array[BWD];
            }
            else
                SVT_LOG("Error in GOp indexing\n");
            av1_rps->refresh_frame_mask = 1 << (lay3_idx);
            break;

        case 4:
            // update RPS for the overlay frame.
            if (picture_control_set_ptr->is_overlay) {
                //{ 0, 0, 0, 0}         // GOP Index 1 - Ref List 0
                //{ 0, 0, 0, 0 }       // GOP Index 1 - Ref List 1
                av1_rps->ref_dpb_index[LAST] = base1_idx;
                av1_rps->ref_dpb_index[LAST2] = base1_idx;
                av1_rps->ref_dpb_index[LAST3] = base1_idx;
                av1_rps->ref_dpb_index[GOLD] = base1_idx;
                av1_rps->ref_dpb_index[BWD] = base1_idx;
                av1_rps->ref_dpb_index[ALT2] = base1_idx;
                av1_rps->ref_dpb_index[ALT] = base1_idx;

                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, 0);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
                av1_rps->ref_poc_array[ALT] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[2]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[1]);
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
                av1_rps->ref_dpb_index[ALT2] = av1_rps->ref_dpb_index[BWD];
                av1_rps->ref_dpb_index[ALT] = av1_rps->ref_dpb_index[BWD];
                gop_i = 15;
                av1_rps->ref_poc_array[LAST] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[0]);
                av1_rps->ref_poc_array[LAST2] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[1]);
                av1_rps->ref_poc_array[LAST3] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[2]);
                av1_rps->ref_poc_array[GOLD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list0[3]);

                av1_rps->ref_poc_array[BWD] = get_ref_poc(context_ptr, picture_control_set_ptr->picture_number, five_level_hierarchical_pred_struct[gop_i].ref_list1[0]);
                av1_rps->ref_poc_array[ALT2] = av1_rps->ref_poc_array[BWD];
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

        if (!set_frame_display_params(picture_control_set_ptr, context_ptr, mini_gop_index)) {
            if (picture_control_set_ptr->is_used_as_reference_flag && picture_index != 0 && picture_index != 8)
            {
                frm_hdr->show_frame = EB_FALSE;
                picture_control_set_ptr->has_show_existing = EB_FALSE;
            } else {
                frm_hdr->show_frame = EB_TRUE;
                picture_control_set_ptr->has_show_existing = EB_TRUE;

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
                    SVT_LOG("Error in GOP indexing for hierarchical level %d\n", picture_control_set_ptr->hierarchical_levels);
            }
        }

        //last pic in MiniGop: Base layer toggling
        //mini GOP toggling since last Key Frame.
        //a regular I keeps the toggling process and does not reset the toggle.  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
        //whoever needs a miniGOP Level toggling, this is the time
        if (((picture_index == (context_ptr->mini_gop_end_index[mini_gop_index]) &&
                        !picture_control_set_ptr->is_overlay &&
                        sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)) ||
                ((sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_P ||
                  sequence_control_set_ptr->static_config.pred_structure == EB_PRED_LOW_DELAY_B) &&
                 (picture_control_set_ptr->temporal_layer_index == 0))) {

            //Layer0 toggle 0->1->2
            context_ptr->lay0_toggle = circ_inc(3, 1, context_ptr->lay0_toggle);
            //Layer1 toggle 3->4
            context_ptr->lay1_toggle = 1 - context_ptr->lay1_toggle;
        }
    } else {
        SVT_LOG("Error: Not supported GOP structure!");
        exit(0);
    }
}

/***************************************************************************************************
// Perform Required Picture Analysis Processing for the Overlay frame
***************************************************************************************************/
void perform_simple_picture_analysis_for_overlay(PictureParentControlSet     *picture_control_set_ptr) {
    EbPictureBufferDesc           *input_padded_picture_ptr;
    EbPictureBufferDesc           *input_picture_ptr;
    EbPaReferenceObject           *paReferenceObject;
    uint32_t                        picture_width_in_sb;
    uint32_t                        pictureHeighInLcu;
    uint32_t                        sb_total_count;

    SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    input_picture_ptr               = picture_control_set_ptr->enhanced_picture_ptr;
    paReferenceObject               = (EbPaReferenceObject*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
    input_padded_picture_ptr        = (EbPictureBufferDesc*)paReferenceObject->input_padded_picture_ptr;
    picture_width_in_sb = (sequence_control_set_ptr->seq_header.max_frame_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
    pictureHeighInLcu   = (sequence_control_set_ptr->seq_header.max_frame_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
    sb_total_count      = picture_width_in_sb * pictureHeighInLcu;

    // Pad pictures to multiple min cu size
    PadPictureToMultipleOfMinCuSizeDimensions(
        sequence_control_set_ptr,
        input_picture_ptr);

    // Pre processing operations performed on the input picture
    PicturePreProcessingOperations(
        picture_control_set_ptr,
        sequence_control_set_ptr,
        sb_total_count);
    if (input_picture_ptr->color_format >= EB_YUV422) {
        // Jing: Do the conversion of 422/444=>420 here since it's multi-threaded kernel
        //       Reuse the Y, only add cb/cr in the newly created buffer desc
        //       NOTE: since denoise may change the src, so this part is after PicturePreProcessingOperations()
        picture_control_set_ptr->chroma_downsampled_picture_ptr->buffer_y = input_picture_ptr->buffer_y;
        DownSampleChroma(input_picture_ptr, picture_control_set_ptr->chroma_downsampled_picture_ptr);
    }
    else
        picture_control_set_ptr->chroma_downsampled_picture_ptr = input_picture_ptr;
    // Pad input picture to complete border LCUs
    PadPictureToMultipleOfLcuDimensions(
        input_padded_picture_ptr);
    // 1/4 & 1/16 input picture decimation
    DownsampleDecimationInputPicture(
        picture_control_set_ptr,
        input_padded_picture_ptr,
        (EbPictureBufferDesc*)paReferenceObject->quarter_decimated_picture_ptr,
        (EbPictureBufferDesc*)paReferenceObject->sixteenth_decimated_picture_ptr);

    // 1/4 & 1/16 input picture downsampling through filtering
    if (sequence_control_set_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
        DownsampleFilteringInputPicture(
            picture_control_set_ptr,
            input_padded_picture_ptr,
            (EbPictureBufferDesc*)paReferenceObject->quarter_filtered_picture_ptr,
            (EbPictureBufferDesc*)paReferenceObject->sixteenth_filtered_picture_ptr);
    }
    // Gathering statistics of input picture, including Variance Calculation, Histogram Bins
    GatheringPictureStatistics(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        picture_control_set_ptr->chroma_downsampled_picture_ptr, //420 input_picture_ptr
        input_padded_picture_ptr,
        paReferenceObject->sixteenth_decimated_picture_ptr, // Hsan: always use decimated until studying the trade offs
        sb_total_count);

    picture_control_set_ptr->sc_content_detected = picture_control_set_ptr->alt_ref_ppcs_ptr->sc_content_detected;
}
/***************************************************************************************************
 * Initialize the overlay frame
***************************************************************************************************/
void initialize_overlay_frame(PictureParentControlSet     *picture_control_set_ptr) {
    picture_control_set_ptr->fade_out_from_black = picture_control_set_ptr->alt_ref_ppcs_ptr->fade_out_from_black;
    picture_control_set_ptr->fade_in_to_black = picture_control_set_ptr->alt_ref_ppcs_ptr->fade_in_to_black;
    picture_control_set_ptr->scene_change_flag = EB_FALSE;
    picture_control_set_ptr->cra_flag = EB_FALSE;
    picture_control_set_ptr->idr_flag = EB_FALSE;
    picture_control_set_ptr->target_bit_rate = picture_control_set_ptr->alt_ref_ppcs_ptr->target_bit_rate;
    picture_control_set_ptr->last_idr_picture = picture_control_set_ptr->alt_ref_ppcs_ptr->last_idr_picture;
    picture_control_set_ptr->use_rps_in_sps = EB_FALSE;
    picture_control_set_ptr->open_gop_cra_flag = EB_FALSE;
    picture_control_set_ptr->pred_structure = picture_control_set_ptr->alt_ref_ppcs_ptr->pred_structure;
    picture_control_set_ptr->pred_struct_ptr = picture_control_set_ptr->alt_ref_ppcs_ptr->pred_struct_ptr;
    picture_control_set_ptr->hierarchical_levels = picture_control_set_ptr->alt_ref_ppcs_ptr->hierarchical_levels;
    picture_control_set_ptr->hierarchical_layers_diff = 0;
    picture_control_set_ptr->init_pred_struct_position_flag = EB_FALSE;
    picture_control_set_ptr->pre_assignment_buffer_count = picture_control_set_ptr->alt_ref_ppcs_ptr->pre_assignment_buffer_count;

    perform_simple_picture_analysis_for_overlay(picture_control_set_ptr);
 }

/***************************************************************************************************
 * Helper function. Compare two frames: center frame and target frame. Return the summation of
 * absolute difference between the two frames from a histogram of luma values
***************************************************************************************************/

static __inline uint32_t compute_luma_sad_between_center_and_target_frame(
    int center_index,
    int target_frame_index,
    PictureParentControlSet *picture_control_set_ptr,
    SequenceControlSet *sequence_control_set_ptr) {

    int32_t center_sum = 0, altref_sum = 0;
    uint32_t ahd = 0;

    for (int bin = 0; bin < HISTOGRAM_NUMBER_OF_BINS; ++bin) {
        center_sum = 0, altref_sum = 0;
        for (uint32_t regionInPictureWidthIndex = 0; regionInPictureWidthIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_width; regionInPictureWidthIndex++) {
            for (uint32_t regionInPictureHeightIndex = 0; regionInPictureHeightIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_height; regionInPictureHeightIndex++) {
                center_sum += picture_control_set_ptr->temp_filt_pcs_list[center_index]->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][bin];
                altref_sum += picture_control_set_ptr->temp_filt_pcs_list[target_frame_index]->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][bin];
            }
        }
        ahd += ABS(center_sum - altref_sum);
    }
    return ahd;
}

/***************************************************************************************************
 * Picture Decision Kernel
 *
 * Notes on the Picture Decision:
 *
 * The Picture Decision process performs multi-picture level decisions, including setting of the prediction structure,
 * setting the picture type and scene change detection.
 *
 * Inputs:
 * Input Picture
 *   -Input Picture Data
 *
 *  Outputs:
 *   -Picture Control Set with fully available PA Reference List
 *
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
 *      /   B   \
 *     /   / \   \
 *    I-----------B
 *
 *  From this base structure, the following RPS positions are derived:
 *
 *    p   p       b   b       p   p
 *     \   \     / \ / \     /   /
 *      P   \   /   B   \   /   P
 *       \   \ /   / \   \ /   /
 *        ----I-----------B----
 *
 *    L L L   I  [ Normal ]   T T T
 *    2 1 0   n               0 1 2
 *            t
 *            r
 *            a
 *
 *  The RPS is composed of Leading Picture [L2-L0], Intra (CRA), Base/Normal Pictures,
 *    and Trailing Pictures [T0-T2]. Generally speaking, Leading Pictures are useful
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
 *
 ***************************************************************************************************/
void* picture_decision_kernel(void *input_ptr)
{
    EbThreadContext               *thread_context_ptr = (EbThreadContext*)input_ptr;
    PictureDecisionContext        *context_ptr = (PictureDecisionContext*)thread_context_ptr->priv;

    PictureParentControlSet       *picture_control_set_ptr;
    FrameHeader                   *frm_hdr;

    EncodeContext                 *encode_context_ptr;
    SequenceControlSet            *sequence_control_set_ptr;

    EbObjectWrapper               *inputResultsWrapperPtr;
    PictureAnalysisResults        *inputResultsPtr;

    EbObjectWrapper               *outputResultsWrapperPtr;
    PictureDecisionResults        *outputResultsPtr;

    PredictionStructureEntry      *predPositionPtr;

    EbBool                          preAssignmentBufferFirstPassFlag;
    EB_SLICE                         picture_type;

    PictureDecisionReorderEntry   *queueEntryPtr;
    int32_t                           queueEntryIndex;

    int32_t                           previousEntryIndex;

    PaReferenceQueueEntry         *inputEntryPtr = (PaReferenceQueueEntry*)EB_NULL;;
    uint32_t                           inputQueueIndex;

    PaReferenceQueueEntry         *paReferenceEntryPtr;
    uint32_t                           paReferenceQueueIndex;

    uint64_t                           ref_poc;

    uint32_t                           depIdx;
    uint64_t                           depPoc;

    uint32_t                           depListCount;

    // Dynamic GOP
    uint32_t                           mini_gop_index;
    uint32_t                           pictureIndex;

    EbBool                          windowAvail, framePasseThru;
    uint32_t                           windowIndex;
    uint32_t                           entryIndex;
    PictureParentControlSet        *ParentPcsWindow[FUTURE_WINDOW_WIDTH + 2];

    // Debug
    uint64_t                           loopCount = 0;

    for (;;) {
        // Get Input Full Object
        eb_get_full_object(
            context_ptr->picture_analysis_results_input_fifo_ptr,
            &inputResultsWrapperPtr);

        inputResultsPtr = (PictureAnalysisResults*)inputResultsWrapperPtr->object_ptr;
        picture_control_set_ptr = (PictureParentControlSet*)inputResultsPtr->picture_control_set_wrapper_ptr->object_ptr;
        sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
        frm_hdr = &picture_control_set_ptr->frm_hdr;
        encode_context_ptr = (EncodeContext*)sequence_control_set_ptr->encode_context_ptr;
        loopCount++;

        // Input Picture Analysis Results into the Picture Decision Reordering Queue
        // P.S. Since the prior Picture Analysis processes stage is multithreaded, inputs to the Picture Decision Process
        // can arrive out-of-display-order, so a the Picture Decision Reordering Queue is used to enforce processing of
        // pictures in display order
        if (!picture_control_set_ptr->is_overlay ) {
            queueEntryIndex = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index]->picture_number);
            queueEntryIndex += encode_context_ptr->picture_decision_reorder_queue_head_index;
            queueEntryIndex = (queueEntryIndex > PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1) ? queueEntryIndex - PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH : queueEntryIndex;
            queueEntryPtr = encode_context_ptr->picture_decision_reorder_queue[queueEntryIndex];
            if (queueEntryPtr->parent_pcs_wrapper_ptr != NULL) {
                CHECK_REPORT_ERROR_NC(
                    encode_context_ptr->app_callback_ptr,
                    EB_ENC_PD_ERROR8);
            }
            else {
                queueEntryPtr->parent_pcs_wrapper_ptr = inputResultsPtr->picture_control_set_wrapper_ptr;
                queueEntryPtr->picture_number = picture_control_set_ptr->picture_number;
            }

            picture_control_set_ptr->pic_decision_reorder_queue_idx = queueEntryIndex;
        }
        // Process the head of the Picture Decision Reordering Queue (Entry N)
        // P.S. The Picture Decision Reordering Queue should be parsed in the display order to be able to construct a pred structure
        queueEntryPtr = encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index];

        while (queueEntryPtr->parent_pcs_wrapper_ptr != EB_NULL) {
            if (((PictureParentControlSet *)(queueEntryPtr->parent_pcs_wrapper_ptr->object_ptr))->end_of_sequence_flag == EB_TRUE) {
                framePasseThru = EB_TRUE;
            }
            else
                framePasseThru = EB_FALSE;
            windowAvail = EB_TRUE;
            previousEntryIndex = QUEUE_GET_PREVIOUS_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index);

            ParentPcsWindow[0] = ParentPcsWindow[1] = ParentPcsWindow[2] = ParentPcsWindow[3] = ParentPcsWindow[4] = ParentPcsWindow[5] = NULL;
            //for poc 0, ignore previous frame check
            if (queueEntryPtr->picture_number > 0 && encode_context_ptr->picture_decision_reorder_queue[previousEntryIndex]->parent_pcs_wrapper_ptr == NULL)
                windowAvail = EB_FALSE;
            else {

                ParentPcsWindow[0] = queueEntryPtr->picture_number > 0 ? (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[previousEntryIndex]->parent_pcs_wrapper_ptr->object_ptr : NULL;
                ParentPcsWindow[1] = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index]->parent_pcs_wrapper_ptr->object_ptr;
                for (windowIndex = 0; windowIndex < sequence_control_set_ptr->scd_delay; windowIndex++) {
                    entryIndex = QUEUE_GET_NEXT_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index, windowIndex + 1);
                    if (encode_context_ptr->picture_decision_reorder_queue[entryIndex]->parent_pcs_wrapper_ptr == NULL) {
                        windowAvail = EB_FALSE;
                        break;
                    }
                    else if (((PictureParentControlSet *)(encode_context_ptr->picture_decision_reorder_queue[entryIndex]->parent_pcs_wrapper_ptr->object_ptr))->end_of_sequence_flag == EB_TRUE) {
                        windowAvail = EB_FALSE;
                        framePasseThru = EB_TRUE;
                        break;
                    }else {
                        ParentPcsWindow[2 + windowIndex] = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[entryIndex]->parent_pcs_wrapper_ptr->object_ptr;
                    }
                }
            }

            picture_control_set_ptr = (PictureParentControlSet*)queueEntryPtr->parent_pcs_wrapper_ptr->object_ptr;
            frm_hdr = &picture_control_set_ptr->frm_hdr;
            picture_control_set_ptr->fade_out_from_black = 0;

            picture_control_set_ptr->fade_in_to_black = 0;
            if (picture_control_set_ptr->idr_flag == EB_TRUE)
                context_ptr->last_solid_color_frame_poc = 0xFFFFFFFF;
            if (windowAvail == EB_TRUE && queueEntryPtr->picture_number > 0) {
                if (sequence_control_set_ptr->static_config.scene_change_detection) {
                    picture_control_set_ptr->scene_change_flag = SceneTransitionDetector(
                        context_ptr,
                        sequence_control_set_ptr,
                        ParentPcsWindow,
                        FUTURE_WINDOW_WIDTH);
                }
                else
                    picture_control_set_ptr->scene_change_flag = EB_FALSE;
                picture_control_set_ptr->cra_flag = (picture_control_set_ptr->scene_change_flag == EB_TRUE) ?
                    EB_TRUE :
                    picture_control_set_ptr->cra_flag;

                // Store scene change in context
                context_ptr->is_scene_change_detected = picture_control_set_ptr->scene_change_flag;
            }

            if (windowAvail == EB_TRUE || framePasseThru == EB_TRUE)
            {
                // Place the PCS into the Pre-Assignment Buffer
                // P.S. The Pre-Assignment Buffer is used to store a whole pre-structure
                encode_context_ptr->pre_assignment_buffer[encode_context_ptr->pre_assignment_buffer_count] = queueEntryPtr->parent_pcs_wrapper_ptr;

                // Setup the PCS & SCS
                picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[encode_context_ptr->pre_assignment_buffer_count]->object_ptr;
                frm_hdr = &picture_control_set_ptr->frm_hdr;
                // Set the POC Number
                picture_control_set_ptr->picture_number = (encode_context_ptr->current_input_poc + 1) /*& ((1 << sequence_control_set_ptr->bits_for_picture_order_count)-1)*/;
                encode_context_ptr->current_input_poc = picture_control_set_ptr->picture_number;

                picture_control_set_ptr->pred_structure = sequence_control_set_ptr->static_config.pred_structure;

                picture_control_set_ptr->hierarchical_layers_diff = 0;

                picture_control_set_ptr->init_pred_struct_position_flag = EB_FALSE;

                picture_control_set_ptr->target_bit_rate = sequence_control_set_ptr->static_config.target_bit_rate;

                ReleasePrevPictureFromReorderQueue(
                    encode_context_ptr);

                // If the Intra period length is 0, then introduce an intra for every picture
                if (sequence_control_set_ptr->intra_period_length == 0)
                    picture_control_set_ptr->cra_flag = EB_TRUE;
                // If an #IntraPeriodLength has passed since the last Intra, then introduce a CRA or IDR based on Intra Refresh type
                else if (sequence_control_set_ptr->intra_period_length != -1) {
                    picture_control_set_ptr->cra_flag =
                        (sequence_control_set_ptr->intra_refresh_type != CRA_REFRESH) ?
                        picture_control_set_ptr->cra_flag :
                        (encode_context_ptr->intra_period_position == (uint32_t)sequence_control_set_ptr->intra_period_length) ?
                        EB_TRUE :
                        picture_control_set_ptr->cra_flag;

                    picture_control_set_ptr->idr_flag =
                        (sequence_control_set_ptr->intra_refresh_type != IDR_REFRESH) ?
                        picture_control_set_ptr->idr_flag :
                        (encode_context_ptr->intra_period_position == (uint32_t)sequence_control_set_ptr->intra_period_length) ?
                        EB_TRUE :
                        picture_control_set_ptr->idr_flag;
                }

                encode_context_ptr->pre_assignment_buffer_eos_flag = (picture_control_set_ptr->end_of_sequence_flag) ? (uint32_t)EB_TRUE : encode_context_ptr->pre_assignment_buffer_eos_flag;

                // Increment the Pre-Assignment Buffer Intra Count
                encode_context_ptr->pre_assignment_buffer_intra_count += (picture_control_set_ptr->idr_flag || picture_control_set_ptr->cra_flag);
                encode_context_ptr->pre_assignment_buffer_idr_count += picture_control_set_ptr->idr_flag;
                encode_context_ptr->pre_assignment_buffer_count += 1;

                if (sequence_control_set_ptr->static_config.rate_control_mode)
                {
                    // Increment the Intra Period Position
                    encode_context_ptr->intra_period_position = (encode_context_ptr->intra_period_position == (uint32_t)sequence_control_set_ptr->intra_period_length) ? 0 : encode_context_ptr->intra_period_position + 1;
                }
                else
                {
                    // Increment the Intra Period Position
                    encode_context_ptr->intra_period_position = ((encode_context_ptr->intra_period_position == (uint32_t)sequence_control_set_ptr->intra_period_length) || (picture_control_set_ptr->scene_change_flag == EB_TRUE)) ? 0 : encode_context_ptr->intra_period_position + 1;
                }

                // Determine if Pictures can be released from the Pre-Assignment Buffer
                if ((encode_context_ptr->pre_assignment_buffer_intra_count > 0) ||
                    (encode_context_ptr->pre_assignment_buffer_count == (uint32_t)(1 << sequence_control_set_ptr->static_config.hierarchical_levels)) ||
                    (encode_context_ptr->pre_assignment_buffer_eos_flag == EB_TRUE) ||
                    (picture_control_set_ptr->pred_structure == EB_PRED_LOW_DELAY_P) ||
                    (picture_control_set_ptr->pred_structure == EB_PRED_LOW_DELAY_B))
                {
                    // Initialize Picture Block Params
                    context_ptr->mini_gop_start_index[0] = 0;
                    context_ptr->mini_gop_end_index[0] = encode_context_ptr->pre_assignment_buffer_count - 1;
                    context_ptr->mini_gop_length[0] = encode_context_ptr->pre_assignment_buffer_count;

                    context_ptr->mini_gop_hierarchical_levels[0] = sequence_control_set_ptr->static_config.hierarchical_levels;
                    context_ptr->mini_gop_intra_count[0] = encode_context_ptr->pre_assignment_buffer_intra_count;
                    context_ptr->mini_gop_idr_count[0] = encode_context_ptr->pre_assignment_buffer_idr_count;
                    context_ptr->total_number_of_mini_gops = 1;

                    encode_context_ptr->previous_mini_gop_hierarchical_levels = (picture_control_set_ptr->picture_number == 0) ?
                        sequence_control_set_ptr->static_config.hierarchical_levels :
                        encode_context_ptr->previous_mini_gop_hierarchical_levels;

                    {
                        if (sequence_control_set_ptr->static_config.hierarchical_levels == 1) {
                            //minigop 2 case
                            context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops] = 0;
                            context_ptr->mini_gop_end_index[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - 1;
                            context_ptr->mini_gop_length[context_ptr->total_number_of_mini_gops] = encode_context_ptr->pre_assignment_buffer_count - context_ptr->mini_gop_start_index[context_ptr->total_number_of_mini_gops];
                            context_ptr->mini_gop_hierarchical_levels[context_ptr->total_number_of_mini_gops] = 2;
                        } else {
                            //minigop 4,8,16
                            if (encode_context_ptr->pre_assignment_buffer_count > 1) {
                                initialize_mini_gop_activity_array(
                                        context_ptr);

                                if (encode_context_ptr->pre_assignment_buffer_count == 16)
                                    context_ptr->mini_gop_activity_array[L5_0_INDEX] = EB_FALSE;
                                else if (sequence_control_set_ptr->static_config.hierarchical_levels >= 3) {
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

                    GenerateMiniGopRps(
                        context_ptr,
                        encode_context_ptr);

                    // Loop over Mini GOPs

                    for (mini_gop_index = 0; mini_gop_index < context_ptr->total_number_of_mini_gops; ++mini_gop_index) {
                        preAssignmentBufferFirstPassFlag = EB_TRUE;
                        {
                            update_base_layer_reference_queue_dependent_count(
                                context_ptr,
                                encode_context_ptr,
                                sequence_control_set_ptr,
                                mini_gop_index);

                            // Keep track of the number of hierarchical levels of the latest implemented mini GOP
                            encode_context_ptr->previous_mini_gop_hierarchical_levels = context_ptr->mini_gop_hierarchical_levels[mini_gop_index];
                        }

                        // 1st Loop over Pictures in the Pre-Assignment Buffer
                        for (pictureIndex = context_ptr->mini_gop_start_index[mini_gop_index]; pictureIndex <= context_ptr->mini_gop_end_index[mini_gop_index]; ++pictureIndex) {
                            picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pictureIndex]->object_ptr;
                            sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
                            frm_hdr = &picture_control_set_ptr->frm_hdr;
                            // Keep track of the mini GOP size to which the input picture belongs - needed @ PictureManagerProcess()
                            picture_control_set_ptr->pre_assignment_buffer_count = context_ptr->mini_gop_length[mini_gop_index];

                            // Update the Pred Structure if cutting short a Random Access period
                            if ((context_ptr->mini_gop_length[mini_gop_index] < picture_control_set_ptr->pred_struct_ptr->pred_struct_period || context_ptr->mini_gop_idr_count[mini_gop_index] > 0) &&
                                picture_control_set_ptr->pred_struct_ptr->pred_type == EB_PRED_RANDOM_ACCESS &&
                                picture_control_set_ptr->idr_flag == EB_FALSE &&
                                picture_control_set_ptr->cra_flag == EB_FALSE)
                            {
                                // Correct the Pred Index before switching structures
                                if (preAssignmentBufferFirstPassFlag == EB_TRUE)
                                    encode_context_ptr->pred_struct_position -= picture_control_set_ptr->pred_struct_ptr->init_pic_index;
                                picture_control_set_ptr->pred_struct_ptr = get_prediction_structure(
                                    encode_context_ptr->prediction_structure_group_ptr,
                                    EB_PRED_LOW_DELAY_P,
                                    sequence_control_set_ptr->reference_count,
                                    picture_control_set_ptr->hierarchical_levels);

                                // Set the RPS Override Flag - this current only will convert a Random Access structure to a Low Delay structure
                                picture_control_set_ptr->use_rps_in_sps = EB_FALSE;
                                picture_control_set_ptr->open_gop_cra_flag = EB_FALSE;

                                picture_type = P_SLICE;
                            }
                            // Open GOP CRA - adjust the RPS
                            else if ((context_ptr->mini_gop_length[mini_gop_index] == picture_control_set_ptr->pred_struct_ptr->pred_struct_period) &&

                                (picture_control_set_ptr->pred_struct_ptr->pred_type == EB_PRED_RANDOM_ACCESS || picture_control_set_ptr->pred_struct_ptr->temporal_layer_count == 1) &&
                                picture_control_set_ptr->idr_flag == EB_FALSE &&
                                picture_control_set_ptr->cra_flag == EB_TRUE)
                            {
                                picture_control_set_ptr->use_rps_in_sps = EB_FALSE;
                                picture_control_set_ptr->open_gop_cra_flag = EB_TRUE;

                                picture_type = I_SLICE;
                            }
                            else {
                                picture_control_set_ptr->use_rps_in_sps = EB_FALSE;
                                picture_control_set_ptr->open_gop_cra_flag = EB_FALSE;

                                // Set the Picture Type
                                picture_type =
                                    (picture_control_set_ptr->idr_flag) ? I_SLICE :
                                    (picture_control_set_ptr->cra_flag) ? I_SLICE :
                                    (picture_control_set_ptr->pred_structure == EB_PRED_LOW_DELAY_P) ? P_SLICE :
                                    (picture_control_set_ptr->pred_structure == EB_PRED_LOW_DELAY_B) ? B_SLICE :
                                    (picture_control_set_ptr->pre_assignment_buffer_count == picture_control_set_ptr->pred_struct_ptr->pred_struct_period) ? ((pictureIndex == context_ptr->mini_gop_end_index[mini_gop_index] && sequence_control_set_ptr->static_config.base_layer_switch_mode) ? P_SLICE : B_SLICE) :

                                    (encode_context_ptr->pre_assignment_buffer_eos_flag) ? P_SLICE :
                                    B_SLICE;
                            }
                            // If mini GOP switch, reset position
                            encode_context_ptr->pred_struct_position = (picture_control_set_ptr->init_pred_struct_position_flag) ?
                                picture_control_set_ptr->pred_struct_ptr->init_pic_index :
                                encode_context_ptr->pred_struct_position;

                            // If Intra, reset position
                            if (picture_control_set_ptr->idr_flag == EB_TRUE)
                                encode_context_ptr->pred_struct_position = picture_control_set_ptr->pred_struct_ptr->init_pic_index;
                            else if (picture_control_set_ptr->cra_flag == EB_TRUE && context_ptr->mini_gop_length[mini_gop_index] < picture_control_set_ptr->pred_struct_ptr->pred_struct_period)
                                encode_context_ptr->pred_struct_position = picture_control_set_ptr->pred_struct_ptr->init_pic_index;
                            else if (encode_context_ptr->elapsed_non_cra_count == 0) {
                                // If we are the picture directly after a CRA, we have to not use references that violate the CRA
                                encode_context_ptr->pred_struct_position = picture_control_set_ptr->pred_struct_ptr->init_pic_index + 1;
                            }
                            // Elif Scene Change, determine leading and trailing pictures
                            //else if (encode_context_ptr->pre_assignment_buffer_scene_change_count > 0) {
                            //    if(bufferIndex < encode_context_ptr->pre_assignment_buffer_scene_change_index) {
                            //        ++encode_context_ptr->pred_struct_position;
                            //        picture_type = P_SLICE;
                            //    }
                            //    else {
                            //        encode_context_ptr->pred_struct_position = picture_control_set_ptr->pred_struct_ptr->init_pic_index + encode_context_ptr->pre_assignment_buffer_count - bufferIndex - 1;
                            //    }
                            //}
                            // Else, Increment the position normally
                            else
                                ++encode_context_ptr->pred_struct_position;
                            // The poc number of the latest IDR picture is stored so that last_idr_picture (present in PCS) for the incoming pictures can be updated.
                            // The last_idr_picture is used in reseting the poc (in entropy coding) whenever IDR is encountered.
                            // Note IMP: This logic only works when display and decode order are the same. Currently for Random Access, IDR is inserted (similar to CRA) by using trailing P pictures (low delay fashion) and breaking prediction structure.
                            // Note: When leading P pictures are implemented, this logic has to change..
                            if (picture_control_set_ptr->idr_flag == EB_TRUE)
                                encode_context_ptr->last_idr_picture = picture_control_set_ptr->picture_number;
                            else
                                picture_control_set_ptr->last_idr_picture = encode_context_ptr->last_idr_picture;
                            // Cycle the PredStructPosition if its overflowed
                            encode_context_ptr->pred_struct_position = (encode_context_ptr->pred_struct_position == picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_count) ?
                                encode_context_ptr->pred_struct_position - picture_control_set_ptr->pred_struct_ptr->pred_struct_period :
                                encode_context_ptr->pred_struct_position;

                            predPositionPtr = picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[encode_context_ptr->pred_struct_position];
                            if (sequence_control_set_ptr->static_config.enable_overlays == EB_TRUE) {
                                // At this stage we know the prediction structure and the location of ALT_REF pictures.
                                // For every ALTREF picture, there is an overlay picture. They extra pictures are released
                                // is_alt_ref flag is set for non-slice base layer pictures
                                if (predPositionPtr->temporal_layer_index == 0 && picture_type != I_SLICE) {
                                    picture_control_set_ptr->is_alt_ref = 1;
                                    frm_hdr->show_frame = 0;
                                    ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pictureIndex - 1]->object_ptr)->has_show_existing = EB_FALSE;
                                }
                                // release the overlay PCS for non alt ref pictures. First picture does not have overlay PCS
                                else if (picture_control_set_ptr->picture_number) {
                                    eb_release_object(picture_control_set_ptr->overlay_ppcs_ptr->input_picture_wrapper_ptr);
                                    // release the pa_reference_picture
                                    eb_release_object(picture_control_set_ptr->overlay_ppcs_ptr->pa_reference_picture_wrapper_ptr);
                                    // release the parent pcs
                                    eb_release_object(picture_control_set_ptr->overlay_ppcs_ptr->p_pcs_wrapper_ptr);
                                    picture_control_set_ptr->overlay_ppcs_ptr = EB_NULL;
                                }
                            }
                            PictureParentControlSet       *cur_picture_control_set_ptr = picture_control_set_ptr;
                            for (uint8_t loop_index = 0; loop_index <= picture_control_set_ptr->is_alt_ref; loop_index++) {
                                // Set missing parts in the overlay  pictures
                                if (loop_index == 1) {
                                    picture_control_set_ptr = picture_control_set_ptr->overlay_ppcs_ptr;
                                    initialize_overlay_frame(picture_control_set_ptr);
                                    picture_type = P_SLICE;
                                }
                                // Set the Slice type
                                picture_control_set_ptr->slice_type = picture_type;
                                ((EbPaReferenceObject*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->object_ptr)->slice_type = picture_control_set_ptr->slice_type;

                                switch (picture_type) {
                                case I_SLICE:

                                    // Reset Prediction Structure Position & Reference Struct Position
                                    if (picture_control_set_ptr->picture_number == 0)
                                        encode_context_ptr->intra_period_position = 0;
                                    encode_context_ptr->elapsed_non_cra_count = 0;

                                    //-------------------------------
                                    // IDR
                                    //-------------------------------
                                    if (picture_control_set_ptr->idr_flag == EB_TRUE) {
                                        // Set CRA flag
                                        picture_control_set_ptr->cra_flag = EB_FALSE;

                                        // Reset the pictures since last IDR counter
                                        encode_context_ptr->elapsed_non_idr_count = 0;
                                        //log latest key frame poc
                                        context_ptr->key_poc = picture_control_set_ptr->picture_number;
                                    }
                                    //-------------------------------
                                    // CRA
                                    //-------------------------------
                                    else {
                                        // Set a Random Access Point if not an IDR
                                        picture_control_set_ptr->cra_flag = EB_TRUE;
                                    }
                                    break;

                                case P_SLICE:
                                case B_SLICE:

                                    // Reset CRA and IDR Flag
                                    picture_control_set_ptr->cra_flag = EB_FALSE;
                                    picture_control_set_ptr->idr_flag = EB_FALSE;

                                    if (loop_index == 0)
                                    {
                                        // Increment & Clip the elapsed Non-IDR Counter. This is clipped rather than allowed to free-run
                                        // inorder to avoid rollover issues.  This assumes that any the GOP period is less than MAX_ELAPSED_IDR_COUNT
                                        encode_context_ptr->elapsed_non_idr_count = MIN(encode_context_ptr->elapsed_non_idr_count + 1, MAX_ELAPSED_IDR_COUNT);
                                        encode_context_ptr->elapsed_non_cra_count = MIN(encode_context_ptr->elapsed_non_cra_count + 1, MAX_ELAPSED_IDR_COUNT);
                                    }

                                    CHECK_REPORT_ERROR(
                                        (picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_count < MAX_ELAPSED_IDR_COUNT),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR1);

                                    break;

                                default:

                                    CHECK_REPORT_ERROR_NC(
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR2);

                                    break;
                                }
                                picture_control_set_ptr->pred_struct_index = (uint8_t)encode_context_ptr->pred_struct_position;
                                if (picture_control_set_ptr->is_overlay) {
                                    // set the overlay frame as non reference frame with max temporal layer index
                                    picture_control_set_ptr->temporal_layer_index = (uint8_t)picture_control_set_ptr->hierarchical_levels;
                                    picture_control_set_ptr->is_used_as_reference_flag = EB_FALSE;
                                }
                                else {
                                    picture_control_set_ptr->temporal_layer_index = (uint8_t)predPositionPtr->temporal_layer_index;
                                    picture_control_set_ptr->is_used_as_reference_flag = predPositionPtr->is_referenced;
                                }

                                preAssignmentBufferFirstPassFlag = EB_FALSE;

                                // Film grain (assigning the random-seed)
                                {
                                    uint16_t *fgn_random_seed_ptr = &picture_control_set_ptr->sequence_control_set_ptr->film_grain_random_seed;
                                    frm_hdr->film_grain_params.random_seed = *fgn_random_seed_ptr;
                                    *fgn_random_seed_ptr += 3381;  // Changing random seed for film grain
                                    if (!(*fgn_random_seed_ptr))     // Random seed should not be zero
                                        *fgn_random_seed_ptr += 7391;
                                }
                                uint32_t pic_index = 0;
                                if (sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {
                                    pic_index = pictureIndex - context_ptr->mini_gop_start_index[mini_gop_index];
                                } else {
                                    // For low delay P or low delay B case, get the the picture_index by mini_gop size
                                    if (sequence_control_set_ptr->static_config.intra_period_length >= 0) {
                                        pic_index = (picture_control_set_ptr->picture_number == 0) ? 0 :
                                            (uint32_t)(((picture_control_set_ptr->picture_number - 1) % (sequence_control_set_ptr->static_config.intra_period_length+1)) % picture_control_set_ptr->pred_struct_ptr->pred_struct_period);
                                    } else {
                                        // intra-period=-1 case, no gop
                                        pic_index = (picture_control_set_ptr->picture_number == 0) ? 0 :
                                            (uint32_t)((picture_control_set_ptr->picture_number - 1) % picture_control_set_ptr->pred_struct_ptr->pred_struct_period);
                                    }
                                }

                                Av1GenerateRpsInfo(
                                    picture_control_set_ptr,
                                    encode_context_ptr,
                                    context_ptr,
                                    pic_index,
                                    mini_gop_index);
                                picture_control_set_ptr->allow_comp_inter_inter = 0;
                                picture_control_set_ptr->is_skip_mode_allowed = 0;

                                frm_hdr->reference_mode = (ReferenceMode)0xFF;

                                if (picture_control_set_ptr->slice_type != I_SLICE) {
                                    picture_control_set_ptr->allow_comp_inter_inter = 1;
                                    if (picture_control_set_ptr->slice_type == P_SLICE) {
                                        picture_control_set_ptr->is_skip_mode_allowed = 0;
                                        frm_hdr->reference_mode = SINGLE_REFERENCE;
                                        picture_control_set_ptr->skip_mode_flag = 0;
                                    }
                                    else if (picture_control_set_ptr->temporal_layer_index == 0) {
                                        frm_hdr->reference_mode = REFERENCE_MODE_SELECT;
                                        frm_hdr->skip_mode_params.skip_mode_flag = 0;
                                    }
                                    else {
                                        frm_hdr->reference_mode = REFERENCE_MODE_SELECT;
                                        picture_control_set_ptr->is_skip_mode_allowed = 1;
                                        picture_control_set_ptr->skip_mode_flag = 1;
                                    }
                                }

                                picture_control_set_ptr->av1_cm->mi_cols = picture_control_set_ptr->sequence_control_set_ptr->seq_header.max_frame_width >> MI_SIZE_LOG2;
                                picture_control_set_ptr->av1_cm->mi_rows = picture_control_set_ptr->sequence_control_set_ptr->seq_header.max_frame_height >> MI_SIZE_LOG2;

                                //Jing: For low delay B/P case, don't alter the bias
                                memset(picture_control_set_ptr->av1_cm->ref_frame_sign_bias, 0, 8 * sizeof(int32_t));
                                if (frm_hdr->reference_mode == REFERENCE_MODE_SELECT &&
                                        picture_control_set_ptr->temporal_layer_index &&
                                        sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS)
                                {
                                    picture_control_set_ptr->av1_cm->ref_frame_sign_bias[ALTREF_FRAME] =
                                        picture_control_set_ptr->av1_cm->ref_frame_sign_bias[ALTREF2_FRAME] =
                                        picture_control_set_ptr->av1_cm->ref_frame_sign_bias[BWDREF_FRAME] = 1;
                                }

                                if (picture_control_set_ptr->slice_type == I_SLICE)
                                    context_ptr->last_i_picture_sc_detection = picture_control_set_ptr->sc_content_detected;
                                else
                                    picture_control_set_ptr->sc_content_detected = context_ptr->last_i_picture_sc_detection;

                                // TODO: put this in EbMotionEstimationProcess?
                                // ME Kernel Multi-Processes Signal(s) derivation
                                signal_derivation_multi_processes_oq(
                                sequence_control_set_ptr,
                                    picture_control_set_ptr);

                            // Set tx_mode
                            frm_hdr->tx_mode = (picture_control_set_ptr->atb_mode) ?
                                TX_MODE_SELECT :
                                TX_MODE_LARGEST;

                                picture_control_set_ptr->use_src_ref = EB_FALSE;
                                picture_control_set_ptr->limit_ois_to_dc_mode_flag = EB_FALSE;

                                // Update the Dependant List Count - If there was an I-frame or Scene Change, then cleanup the Picture Decision PA Reference Queue Dependent Counts
                                if (picture_control_set_ptr->slice_type == I_SLICE)
                                {
                                    inputQueueIndex = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                                    while (inputQueueIndex != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
                                        inputEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[inputQueueIndex];

                                        // Modify Dependent List0
                                        depListCount = inputEntryPtr->list0.list_count;
                                        for (depIdx = 0; depIdx < depListCount; ++depIdx) {
                                            // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                                            // current_input_poc += (current_input_poc < inputEntryPtr->pocNumber) ? (1 << sequence_control_set_ptr->bits_for_picture_order_count) : 0;

                                            depPoc = POC_CIRCULAR_ADD(
                                                inputEntryPtr->picture_number, // can't use a value that gets reset
                                                inputEntryPtr->list0.list[depIdx]/*,
                                                sequence_control_set_ptr->bits_for_picture_order_count*/);

                                                // If Dependent POC is greater or equal to the IDR POC
                                            if (depPoc >= picture_control_set_ptr->picture_number && inputEntryPtr->list0.list[depIdx]) {
                                                inputEntryPtr->list0.list[depIdx] = 0;

                                                // Decrement the Reference's referenceCount
                                                --inputEntryPtr->dependent_count;

                                                CHECK_REPORT_ERROR(
                                                    (inputEntryPtr->dependent_count != ~0u),
                                                    encode_context_ptr->app_callback_ptr,
                                                    EB_ENC_PD_ERROR3);
                                            }
                                        }

                                        // Modify Dependent List1
                                        depListCount = inputEntryPtr->list1.list_count;
                                        for (depIdx = 0; depIdx < depListCount; ++depIdx) {
                                            // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                                            // current_input_poc += (current_input_poc < inputEntryPtr->pocNumber) ? (1 << sequence_control_set_ptr->bits_for_picture_order_count) : 0;

                                            depPoc = POC_CIRCULAR_ADD(
                                                inputEntryPtr->picture_number,
                                                inputEntryPtr->list1.list[depIdx]/*,
                                                sequence_control_set_ptr->bits_for_picture_order_count*/);

                                                // If Dependent POC is greater or equal to the IDR POC
                                            if (((depPoc >= picture_control_set_ptr->picture_number) || (((picture_control_set_ptr->pre_assignment_buffer_count != picture_control_set_ptr->pred_struct_ptr->pred_struct_period) || (picture_control_set_ptr->idr_flag == EB_TRUE)) && (depPoc > (picture_control_set_ptr->picture_number - picture_control_set_ptr->pre_assignment_buffer_count)))) && inputEntryPtr->list1.list[depIdx]) {
                                                inputEntryPtr->list1.list[depIdx] = 0;

                                                // Decrement the Reference's referenceCount
                                                --inputEntryPtr->dependent_count;

                                                CHECK_REPORT_ERROR(
                                                    (inputEntryPtr->dependent_count != ~0u),
                                                    encode_context_ptr->app_callback_ptr,
                                                    EB_ENC_PD_ERROR3);
                                            }
                                        }

                                        // Increment the inputQueueIndex Iterator
                                        inputQueueIndex = (inputQueueIndex == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
                                    }
                                }
                                else if (picture_control_set_ptr->idr_flag == EB_TRUE) {
                                    // Set the Picture Decision PA Reference Entry pointer
                                    inputEntryPtr = (PaReferenceQueueEntry*)EB_NULL;
                                }

                                // Place Picture in Picture Decision PA Reference Queue
                                // there is no need to add the overlay frames in the PA Reference Queue
                                if (!picture_control_set_ptr->is_overlay)
                                {
                                    inputEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_tail_index];
                                    inputEntryPtr->input_object_ptr = picture_control_set_ptr->pa_reference_picture_wrapper_ptr;
                                    inputEntryPtr->picture_number = picture_control_set_ptr->picture_number;
                                    inputEntryPtr->reference_entry_index = encode_context_ptr->picture_decision_pa_reference_queue_tail_index;
                                    inputEntryPtr->is_alt_ref = picture_control_set_ptr->is_alt_ref;
                                    encode_context_ptr->picture_decision_pa_reference_queue_tail_index =
                                        (encode_context_ptr->picture_decision_pa_reference_queue_tail_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->picture_decision_pa_reference_queue_tail_index + 1;

                                    // Check if the Picture Decision PA Reference is full
                                    CHECK_REPORT_ERROR(
                                        (((encode_context_ptr->picture_decision_pa_reference_queue_head_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) ||
                                        (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->input_object_ptr == EB_NULL))),
                                        encode_context_ptr->app_callback_ptr,
                                        EB_ENC_PD_ERROR4);
                                }
                                // Copy the reference lists into the inputEntry and
                                // set the Reference Counts Based on Temporal Layer and how many frames are active
                                picture_control_set_ptr->ref_list0_count = (picture_type == I_SLICE) ? 0 :
                                                                            (picture_control_set_ptr->is_overlay) ? 1 : (uint8_t)predPositionPtr->ref_list0.reference_list_count;
                                picture_control_set_ptr->ref_list1_count = (picture_type == I_SLICE || picture_control_set_ptr->is_overlay) ? 0 : (uint8_t)predPositionPtr->ref_list1.reference_list_count;
                                if (!picture_control_set_ptr->is_overlay) {
                                    inputEntryPtr->list0_ptr = &predPositionPtr->ref_list0;
                                    inputEntryPtr->list1_ptr = &predPositionPtr->ref_list1;
                                }
                                if (!picture_control_set_ptr->is_overlay)
                                {
                                    // Copy the Dependent Lists
                                    // *Note - we are removing any leading picture dependencies for now
                                    inputEntryPtr->list0.list_count = 0;
                                    for (depIdx = 0; depIdx < predPositionPtr->dep_list0.list_count; ++depIdx) {
                                        if (predPositionPtr->dep_list0.list[depIdx] >= 0)
                                            inputEntryPtr->list0.list[inputEntryPtr->list0.list_count++] = predPositionPtr->dep_list0.list[depIdx];
                                    }
                                    inputEntryPtr->list1.list_count = predPositionPtr->dep_list1.list_count;
                                    for (depIdx = 0; depIdx < predPositionPtr->dep_list1.list_count; ++depIdx)
                                        inputEntryPtr->list1.list[depIdx] = predPositionPtr->dep_list1.list[depIdx];
                                    // add the overlay picture to the dependant list
                                    inputEntryPtr->dep_list0_count = (picture_control_set_ptr->is_alt_ref) ? inputEntryPtr->list0.list_count + 1 : inputEntryPtr->list0.list_count;

                                    inputEntryPtr->dep_list1_count = inputEntryPtr->list1.list_count;
                                    inputEntryPtr->dependent_count = inputEntryPtr->dep_list0_count + inputEntryPtr->dep_list1_count;

                                    ((EbPaReferenceObject*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->object_ptr)->dependent_pictures_count = inputEntryPtr->dependent_count;
                                }

                                CHECK_REPORT_ERROR(
                                    (picture_control_set_ptr->pred_struct_ptr->pred_struct_period * REF_LIST_MAX_DEPTH < MAX_ELAPSED_IDR_COUNT),
                                    encode_context_ptr->app_callback_ptr,
                                    EB_ENC_PD_ERROR5);

                                // Reset the PA Reference Lists
                                EB_MEMSET(picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                                EB_MEMSET(picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));

                                EB_MEMSET(picture_control_set_ptr->ref_pic_poc_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                                EB_MEMSET(picture_control_set_ptr->ref_pic_poc_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                            }
                            picture_control_set_ptr = cur_picture_control_set_ptr;

                            if ((sequence_control_set_ptr->enable_altrefs == EB_TRUE &&
                                    sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS &&
                                    sequence_control_set_ptr->static_config.hierarchical_levels >=2) &&

                                ((picture_control_set_ptr->slice_type == I_SLICE && picture_control_set_ptr->sc_content_detected == 0) ||
                                  (picture_control_set_ptr->slice_type != I_SLICE && picture_control_set_ptr->temporal_layer_index == 0)
#if TWO_PASS
                                    || (sequence_control_set_ptr->use_input_stat_file && picture_control_set_ptr->temporal_layer_index == 1 && picture_control_set_ptr->sc_content_detected == 0)
#endif
                                    ) ) {
                                int altref_nframes = picture_control_set_ptr->sequence_control_set_ptr->static_config.altref_nframes;
                                if (picture_control_set_ptr->idr_flag) {

                                    //initilize list
                                    for (int pic_itr = 0; pic_itr < ALTREF_MAX_NFRAMES; pic_itr++)
                                        picture_control_set_ptr->temp_filt_pcs_list[pic_itr] = NULL;

                                    picture_control_set_ptr->temp_filt_pcs_list[0] = picture_control_set_ptr;

                                    uint32_t num_future_pics = 6;
                                    uint32_t num_past_pics = 0;
                                    uint32_t pic_i;
                                    //search reord-queue to get the future pictures
                                    for (pic_i = 0; pic_i < num_future_pics; pic_i++) {
                                        int32_t q_index = QUEUE_GET_NEXT_SPOT(picture_control_set_ptr->pic_decision_reorder_queue_idx, pic_i + 1);
                                        if (encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr != NULL) {
                                            PictureParentControlSet* pcs_itr = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr->object_ptr;
                                            picture_control_set_ptr->temp_filt_pcs_list[pic_i + num_past_pics + 1] = pcs_itr;

                                        }
                                        else
                                           break;
                                    }

                                    picture_control_set_ptr->past_altref_nframes = 0;
                                    picture_control_set_ptr->future_altref_nframes = pic_i;
                                    int index_center = 0;
                                    uint32_t actual_future_pics = picture_control_set_ptr->future_altref_nframes;
                                    int pic_itr, ahd;

                                    int ahd_th = (((sequence_control_set_ptr->seq_header.max_frame_width * sequence_control_set_ptr->seq_header.max_frame_height) * AHD_TH_WEIGHT) / 100);

                                    // Accumulative histogram absolute differences between the central and future frame
                                    for (pic_itr = (index_center + actual_future_pics); pic_itr > index_center; pic_itr--) {
                                        ahd = compute_luma_sad_between_center_and_target_frame(index_center, pic_itr, picture_control_set_ptr, sequence_control_set_ptr);
                                        if (ahd < ahd_th)
                                            break;
                                    }
                                    picture_control_set_ptr->future_altref_nframes = pic_itr - index_center;
                                    //SVT_LOG("\nPOC %d\t PAST %d\t FUTURE %d\n", picture_control_set_ptr->picture_number, picture_control_set_ptr->past_altref_nframes, picture_control_set_ptr->future_altref_nframes);
                                }
                                else
                                {
                                int num_past_pics = altref_nframes / 2;
                                int num_future_pics = altref_nframes - num_past_pics - 1;
                                assert(altref_nframes <= ALTREF_MAX_NFRAMES);

                                //initilize list
                                for (int pic_itr = 0; pic_itr < ALTREF_MAX_NFRAMES; pic_itr++)
                                    picture_control_set_ptr->temp_filt_pcs_list[pic_itr] = NULL;
                                // Jing: Intra CRA case
                                num_past_pics = MIN(num_past_pics, (int)encode_context_ptr->pre_assignment_buffer_count -1);

                                //get previous+current pictures from the the pre-assign buffer
                                for (int pic_itr = 0; pic_itr <= num_past_pics; pic_itr++) {
                                    PictureParentControlSet* pcs_itr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pictureIndex - num_past_pics + pic_itr]->object_ptr;
                                    picture_control_set_ptr->temp_filt_pcs_list[pic_itr] = pcs_itr;
                                }
                                int actual_past_pics = num_past_pics;
                                int actual_future_pics = 0;
                                int pic_i;
                                //search reord-queue to get the future pictures
                                for (pic_i = 0; pic_i < num_future_pics; pic_i++) {
                                    int32_t q_index = QUEUE_GET_NEXT_SPOT(picture_control_set_ptr->pic_decision_reorder_queue_idx, pic_i + 1);
                                    if (encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr != NULL) {
                                        PictureParentControlSet* pcs_itr = (PictureParentControlSet *)encode_context_ptr->picture_decision_reorder_queue[q_index]->parent_pcs_wrapper_ptr->object_ptr;
                                        picture_control_set_ptr->temp_filt_pcs_list[pic_i + num_past_pics + 1] = pcs_itr;
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
                                            if (pcs_itr->picture_number == picture_control_set_ptr->picture_number + pic_i_future + 1) {
                                                picture_control_set_ptr->temp_filt_pcs_list[pic_i_future + num_past_pics + 1] = pcs_itr;
                                                actual_future_pics++;
                                                break; //exist the pre-ass loop, go search the next
                                            }
                                        }
                                    }
                                }
                                int index_center = actual_past_pics;
                                int pic_itr;
                                int ahd;

                                int ahd_th = (((sequence_control_set_ptr->seq_header.max_frame_width * sequence_control_set_ptr->seq_header.max_frame_height) * AHD_TH_WEIGHT) / 100);

                                // Accumulative histogram absolute differences between the central and past frame
                                for (pic_itr = index_center - actual_past_pics; pic_itr < index_center; pic_itr++) {
                                    ahd = compute_luma_sad_between_center_and_target_frame(index_center, pic_itr, picture_control_set_ptr, sequence_control_set_ptr);

                                    if (ahd < ahd_th)
                                        break;
                                }
                                picture_control_set_ptr->past_altref_nframes = actual_past_pics = index_center - pic_itr;

                                // Accumulative histogram absolute differences between the central and past frame
                                for (pic_itr = (index_center + actual_future_pics); pic_itr > index_center; pic_itr--) {
                                    ahd = compute_luma_sad_between_center_and_target_frame(index_center, pic_itr, picture_control_set_ptr, sequence_control_set_ptr);
                                    if (ahd < ahd_th)
                                        break;
                                }
                                picture_control_set_ptr->future_altref_nframes = pic_itr - index_center;
                                //SVT_LOG("\nPOC %d\t PAST %d\t FUTURE %d\n", picture_control_set_ptr->picture_number, picture_control_set_ptr->past_altref_nframes, picture_control_set_ptr->future_altref_nframes);

                                // adjust the temporal filtering pcs buffer to remove unused past pictures
                                if(actual_past_pics != num_past_pics) {

                                    pic_i = 0;
                                    while (picture_control_set_ptr->temp_filt_pcs_list[pic_i] != NULL){
                                        picture_control_set_ptr->temp_filt_pcs_list[pic_i] = picture_control_set_ptr->temp_filt_pcs_list[pic_i + num_past_pics - actual_past_pics];
                                        pic_i++;
                                    }

                                }
                                }

                                picture_control_set_ptr->temp_filt_prep_done = 0;

                                // Start Filtering in ME processes
                                {
                                    int16_t seg_idx;

                                    // Initialize Segments
                                    picture_control_set_ptr->tf_segments_column_count = sequence_control_set_ptr->tf_segment_column_count;
                                    picture_control_set_ptr->tf_segments_row_count    = sequence_control_set_ptr->tf_segment_row_count;
                                    picture_control_set_ptr->tf_segments_total_count = (uint16_t)(picture_control_set_ptr->tf_segments_column_count  * picture_control_set_ptr->tf_segments_row_count);

                                    picture_control_set_ptr->temp_filt_seg_acc = 0;
#if  TWO_PASS
                                    if (picture_control_set_ptr->temporal_layer_index == 0)
                                        picture_control_set_ptr->altref_strength = sequence_control_set_ptr->static_config.altref_strength;
                                    else
                                        picture_control_set_ptr->altref_strength = 2;
#else
                                    picture_control_set_ptr->altref_strength = sequence_control_set_ptr->static_config.altref_strength;
#endif

                                    for (seg_idx = 0; seg_idx < picture_control_set_ptr->tf_segments_total_count; ++seg_idx) {
                                        eb_get_empty_object(
                                            context_ptr->picture_decision_results_output_fifo_ptr,
                                            &outputResultsWrapperPtr);
                                        outputResultsPtr = (PictureDecisionResults*)outputResultsWrapperPtr->object_ptr;
                                        outputResultsPtr->picture_control_set_wrapper_ptr = encode_context_ptr->pre_assignment_buffer[pictureIndex];
                                        outputResultsPtr->segment_index = seg_idx;
                                        outputResultsPtr->task_type = 1;
                                        eb_post_full_object(outputResultsWrapperPtr);
                                    }

                                    eb_block_on_semaphore(picture_control_set_ptr->temp_filt_done_semaphore);
                                }

                            }else
                                picture_control_set_ptr->temporal_filtering_on = EB_FALSE; // set temporal filtering flag OFF for current picture

                        }

                        // Add 1 to the loop for the overlay picture. If the last picture is alt ref, increase the loop by 1 to add the overlay picture
                        uint32_t has_overlay = ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_end_index[mini_gop_index]]->object_ptr)->is_alt_ref ? 1 : 0;
                        for (pictureIndex = context_ptr->mini_gop_start_index[mini_gop_index]; pictureIndex <= context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay; ++pictureIndex) {
                            // 2nd Loop over Pictures in the Pre-Assignment Buffer
                            // Assign the overlay pcs. Since Overlay picture is not added to the picture_decision_pa_reference_queue, in the next stage, the loop finds the alt_ref picture. The reference for overlay frame is hardcoded later
                            if (has_overlay && pictureIndex == context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay) {
                                picture_control_set_ptr = ((PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[context_ptr->mini_gop_end_index[mini_gop_index]]->object_ptr)->overlay_ppcs_ptr;
                                frm_hdr = &picture_control_set_ptr->frm_hdr;
                            }
                            else {
                                picture_control_set_ptr = (PictureParentControlSet*)encode_context_ptr->pre_assignment_buffer[pictureIndex]->object_ptr;
                                frm_hdr = &picture_control_set_ptr->frm_hdr;
                            }
                            picture_control_set_ptr->picture_number_alt = encode_context_ptr->picture_number_alt++;

                            // Set the Decode Order
                            if ((context_ptr->mini_gop_idr_count[mini_gop_index] == 0) &&
                                (context_ptr->mini_gop_length[mini_gop_index] == picture_control_set_ptr->pred_struct_ptr->pred_struct_period) &&
                                (sequence_control_set_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) &&
                                !picture_control_set_ptr->is_overlay) {
                                picture_control_set_ptr->decode_order = encode_context_ptr->decode_base_number + picture_control_set_ptr->pred_struct_ptr->pred_struct_entry_ptr_array[picture_control_set_ptr->pred_struct_index]->decode_order;
                            }
                            else
                                picture_control_set_ptr->decode_order = picture_control_set_ptr->picture_number_alt;
                            encode_context_ptr->terminating_sequence_flag_received = (picture_control_set_ptr->end_of_sequence_flag == EB_TRUE) ?
                                EB_TRUE :
                                encode_context_ptr->terminating_sequence_flag_received;

                            encode_context_ptr->terminating_picture_number = (picture_control_set_ptr->end_of_sequence_flag == EB_TRUE) ?
                                picture_control_set_ptr->picture_number_alt :
                                encode_context_ptr->terminating_picture_number;

                            // Find the Reference in the Picture Decision PA Reference Queue
                            inputQueueIndex = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                            do {
                                // Setup the Picture Decision PA Reference Queue Entry
                                inputEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[inputQueueIndex];

                                // Increment the referenceQueueIndex Iterator
                                inputQueueIndex = (inputQueueIndex == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
                            } while ((inputQueueIndex != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) && (inputEntryPtr->picture_number != picture_control_set_ptr->picture_number));

                            EB_MEMSET(picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                            EB_MEMSET(picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(EbObjectWrapper*));
                            EB_MEMSET(picture_control_set_ptr->ref_pic_poc_array[REF_LIST_0], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                            EB_MEMSET(picture_control_set_ptr->ref_pic_poc_array[REF_LIST_1], 0, REF_LIST_MAX_DEPTH * sizeof(uint64_t));
                            CHECK_REPORT_ERROR(
                                (inputEntryPtr->picture_number == picture_control_set_ptr->picture_number),
                                encode_context_ptr->app_callback_ptr,
                                EB_ENC_PD_ERROR7);

                            // Reset the PA Reference Lists
                            EB_MEMSET(picture_control_set_ptr->ref_pa_pic_ptr_array, 0, 2 * sizeof(EbObjectWrapper*));

                            EB_MEMSET(picture_control_set_ptr->ref_pic_poc_array, 0, 2 * sizeof(uint64_t));

                            // Configure List0
                            if ((picture_control_set_ptr->slice_type == P_SLICE) || (picture_control_set_ptr->slice_type == B_SLICE)) {
                                uint8_t ref_pic_index;
                                for (ref_pic_index = 0; ref_pic_index < picture_control_set_ptr->ref_list0_count; ++ref_pic_index) {
                                    if (picture_control_set_ptr->ref_list0_count) {
                                        // hardcode the reference for the overlay frame
                                        if(picture_control_set_ptr->is_overlay){
                                            paReferenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                                (int32_t)inputEntryPtr->reference_entry_index,
                                                PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                            paReferenceEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[paReferenceQueueIndex];

                                            // Calculate the Ref POC
                                            ref_poc = POC_CIRCULAR_ADD(
                                                picture_control_set_ptr->picture_number,
                                                0);
                                        }
                                        else {
                                            paReferenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                                ((int32_t)inputEntryPtr->reference_entry_index) - inputEntryPtr->list0_ptr->reference_list[ref_pic_index],
                                                PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                            paReferenceEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[paReferenceQueueIndex];

                                            // Calculate the Ref POC
                                            ref_poc = POC_CIRCULAR_ADD(
                                                picture_control_set_ptr->picture_number,
                                                -inputEntryPtr->list0_ptr->reference_list[ref_pic_index] /*
                                                sequence_control_set_ptr->bits_for_picture_order_count*/);
                                        }
                                            // Set the Reference Object
                                        picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_0][ref_pic_index] = paReferenceEntryPtr->input_object_ptr;
                                        picture_control_set_ptr->ref_pic_poc_array[REF_LIST_0][ref_pic_index] = ref_poc;
                                        // Increment the PA Reference's liveCount by the number of tiles in the input picture
                                        eb_object_inc_live_count(
                                            paReferenceEntryPtr->input_object_ptr,
                                            1);
                                        --paReferenceEntryPtr->dependent_count;
                                    }
                                }
                            }

                            // Configure List1
                            if (picture_control_set_ptr->slice_type == B_SLICE) {
                                uint8_t ref_pic_index;
                                for (ref_pic_index = 0; ref_pic_index < picture_control_set_ptr->ref_list1_count; ++ref_pic_index) {
                                    if (picture_control_set_ptr->ref_list1_count) {
                                        paReferenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                            ((int32_t)inputEntryPtr->reference_entry_index) - inputEntryPtr->list1_ptr->reference_list[ref_pic_index],
                                            PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                        paReferenceEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[paReferenceQueueIndex];

                                        // Calculate the Ref POC
                                        ref_poc = POC_CIRCULAR_ADD(
                                            picture_control_set_ptr->picture_number,
                                            -inputEntryPtr->list1_ptr->reference_list[ref_pic_index]/*,
                                            sequence_control_set_ptr->bits_for_picture_order_count*/);
                                        // Set the Reference Object
                                        picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_1][ref_pic_index] = paReferenceEntryPtr->input_object_ptr;
                                        picture_control_set_ptr->ref_pic_poc_array[REF_LIST_1][ref_pic_index] = ref_poc;

                                        // Increment the PA Reference's liveCount by the number of tiles in the input picture
                                        eb_object_inc_live_count(
                                            paReferenceEntryPtr->input_object_ptr,
                                            1);
                                        --paReferenceEntryPtr->dependent_count;
                                    }
                                }
                            }

                            eb_av1_setup_skip_mode_allowed(picture_control_set_ptr);

                            picture_control_set_ptr->is_skip_mode_allowed = frm_hdr->skip_mode_params.skip_mode_allowed;
                            picture_control_set_ptr->skip_mode_flag = picture_control_set_ptr->is_skip_mode_allowed;
                            //SVT_LOG("POC:%i  skip_mode_allowed:%i  REF_SKIP_0: %i   REF_SKIP_1: %i \n",picture_control_set_ptr->picture_number, picture_control_set_ptr->skip_mode_info.skip_mode_allowed, picture_control_set_ptr->skip_mode_info.ref_frame_idx_0, picture_control_set_ptr->skip_mode_info.ref_frame_idx_1);

                            {
                                picture_control_set_ptr->intensity_transition_flag = EB_FALSE;
                                if (picture_control_set_ptr->ref_list0_count)
                                    picture_control_set_ptr->scene_transition_flag[REF_LIST_0] = EB_FALSE;
                                if (picture_control_set_ptr->ref_list1_count)
                                    picture_control_set_ptr->scene_transition_flag[REF_LIST_1] = EB_FALSE;
                            }

                            //set the ref frame types used for this picture,
                            set_all_ref_frame_type(sequence_control_set_ptr, picture_control_set_ptr, picture_control_set_ptr->ref_frame_type_arr, &picture_control_set_ptr->tot_ref_frame_types);
                            // Initialize Segments
                            picture_control_set_ptr->me_segments_column_count = (uint8_t)(sequence_control_set_ptr->me_segment_column_count_array[picture_control_set_ptr->temporal_layer_index]);
                            picture_control_set_ptr->me_segments_row_count = (uint8_t)(sequence_control_set_ptr->me_segment_row_count_array[picture_control_set_ptr->temporal_layer_index]);
                            picture_control_set_ptr->me_segments_total_count = (uint16_t)(picture_control_set_ptr->me_segments_column_count  * picture_control_set_ptr->me_segments_row_count);
                            picture_control_set_ptr->me_segments_completion_mask = 0;

                            // Post the results to the ME processes
                            {
                                uint32_t segment_index;

                                for (segment_index = 0; segment_index < picture_control_set_ptr->me_segments_total_count; ++segment_index) {
                                    // Get Empty Results Object
                                    eb_get_empty_object(
                                            context_ptr->picture_decision_results_output_fifo_ptr,
                                            &outputResultsWrapperPtr);

                                    outputResultsPtr = (PictureDecisionResults*)outputResultsWrapperPtr->object_ptr;
                                    if (picture_control_set_ptr->is_overlay)
                                        outputResultsPtr->picture_control_set_wrapper_ptr = picture_control_set_ptr->p_pcs_wrapper_ptr;
                                    else
                                        outputResultsPtr->picture_control_set_wrapper_ptr = encode_context_ptr->pre_assignment_buffer[pictureIndex];

                                    outputResultsPtr->segment_index = segment_index;
                                    outputResultsPtr->task_type = 0;
                                    // Post the Full Results Object
                                    eb_post_full_object(outputResultsWrapperPtr);
                                }
                            }

                            if (pictureIndex == context_ptr->mini_gop_end_index[mini_gop_index] + has_overlay) {
                                // Increment the Decode Base Number
                                encode_context_ptr->decode_base_number += context_ptr->mini_gop_length[mini_gop_index] + has_overlay;
                            }
                            if (pictureIndex == encode_context_ptr->pre_assignment_buffer_count - 1 + has_overlay) {

                                // Reset the Pre-Assignment Buffer
                                encode_context_ptr->pre_assignment_buffer_count = 0;
                                encode_context_ptr->pre_assignment_buffer_idr_count = 0;
                                encode_context_ptr->pre_assignment_buffer_intra_count = 0;
                                encode_context_ptr->pre_assignment_buffer_scene_change_count = 0;
                                encode_context_ptr->pre_assignment_buffer_eos_flag = EB_FALSE;
                            }
                        }
                    } // End MINI GOPs loop
                }

                // Walk the picture_decision_pa_reference_queue and remove entries that have been completely referenced.
                inputQueueIndex = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                while (inputQueueIndex != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {
                    inputEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[inputQueueIndex];

                    // Remove the entry
                    if ((inputEntryPtr->dependent_count == 0) &&
                        (inputEntryPtr->input_object_ptr)) {
                        // Release the nominal live_count value
                        eb_release_object(inputEntryPtr->input_object_ptr);
                        inputEntryPtr->input_object_ptr = (EbObjectWrapper*)EB_NULL;
                    }

                    // Increment the head_index if the head is null
                    encode_context_ptr->picture_decision_pa_reference_queue_head_index =
                        (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->input_object_ptr) ? encode_context_ptr->picture_decision_pa_reference_queue_head_index :
                        (encode_context_ptr->picture_decision_pa_reference_queue_head_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0
                        : encode_context_ptr->picture_decision_pa_reference_queue_head_index + 1;

                    CHECK_REPORT_ERROR(
                        (((encode_context_ptr->picture_decision_pa_reference_queue_head_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) ||
                        (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->input_object_ptr == EB_NULL))),
                        encode_context_ptr->app_callback_ptr,
                        EB_ENC_PD_ERROR4);

                    // Increment the inputQueueIndex Iterator
                    inputQueueIndex = (inputQueueIndex == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
                }

                // Increment the Picture Decision Reordering Queue Head Ptr
                encode_context_ptr->picture_decision_reorder_queue_head_index = (encode_context_ptr->picture_decision_reorder_queue_head_index == PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->picture_decision_reorder_queue_head_index + 1;

                // Get the next entry from the Picture Decision Reordering Queue (Entry N+1)
                queueEntryPtr = encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index];
            }
            if (windowAvail == EB_FALSE && framePasseThru == EB_FALSE)
                break;
        }

        // Release the Input Results
        eb_release_object(inputResultsWrapperPtr);
    }

    return EB_NULL;
}
// clang-format on
