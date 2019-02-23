/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureAnalysisResults.h"
#include "EbPictureDecisionProcess.h"
#include "EbPictureDecisionResults.h"
#include "EbReferenceObject.h"
#include "EbErrorCodes.h"

/************************************************
 * Defines
 ************************************************/
#define POC_CIRCULAR_ADD(base, offset/*, bits*/)             (/*(((int32_t) (base)) + ((int32_t) (offset)) > ((int32_t) (1 << (bits))))   ? ((base) + (offset) - (1 << (bits))) : \
                                                             (((int32_t) (base)) + ((int32_t) (offset)) < 0)                           ? ((base) + (offset) + (1 << (bits))) : \
                                                                                                                                       */((base) + (offset)))
#define FUTURE_WINDOW_WIDTH                 4
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
#if TUNED_SETTINGS_FOR_M0
#define FC_SKIP_TX_SR_TH                       125 // Fast cost skip tx search threshold.
#endif
#if TUNED_SETTINGS_FOR_M1
#define FC_SKIP_TX_SR_TH_M1                    110 // Fast cost skip tx search threshold.
#endif
 /************************************************
  * Picture Analysis Context Constructor
  ************************************************/
EbErrorType PictureDecisionContextCtor(
    PictureDecisionContext_t **context_dbl_ptr,
    EbFifo_t *pictureAnalysisResultsInputFifoPtr,
    EbFifo_t *pictureDecisionResultsOutputFifoPtr)
{
    PictureDecisionContext_t *context_ptr;
    uint32_t arrayIndex;
    uint32_t arrayRow, arrowColumn;
    EB_MALLOC(PictureDecisionContext_t*, context_ptr, sizeof(PictureDecisionContext_t), EB_N_PTR);
    *context_dbl_ptr = context_ptr;

    context_ptr->pictureAnalysisResultsInputFifoPtr = pictureAnalysisResultsInputFifoPtr;
    context_ptr->pictureDecisionResultsOutputFifoPtr = pictureDecisionResultsOutputFifoPtr;

    EB_MALLOC(uint32_t**, context_ptr->ahdRunningAvgCb, sizeof(uint32_t*) * MAX_NUMBER_OF_REGIONS_IN_WIDTH, EB_N_PTR);

    EB_MALLOC(uint32_t**, context_ptr->ahdRunningAvgCr, sizeof(uint32_t*) * MAX_NUMBER_OF_REGIONS_IN_WIDTH, EB_N_PTR);

    EB_MALLOC(uint32_t**, context_ptr->ahdRunningAvg, sizeof(uint32_t*) * MAX_NUMBER_OF_REGIONS_IN_WIDTH, EB_N_PTR);

    for (arrayIndex = 0; arrayIndex < MAX_NUMBER_OF_REGIONS_IN_WIDTH; arrayIndex++)
    {
        EB_MALLOC(uint32_t*, context_ptr->ahdRunningAvgCb[arrayIndex], sizeof(uint32_t) * MAX_NUMBER_OF_REGIONS_IN_HEIGHT, EB_N_PTR);

        EB_MALLOC(uint32_t*, context_ptr->ahdRunningAvgCr[arrayIndex], sizeof(uint32_t) * MAX_NUMBER_OF_REGIONS_IN_HEIGHT, EB_N_PTR);

        EB_MALLOC(uint32_t*, context_ptr->ahdRunningAvg[arrayIndex], sizeof(uint32_t) * MAX_NUMBER_OF_REGIONS_IN_HEIGHT, EB_N_PTR);
    }

    for (arrayRow = 0; arrayRow < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; arrayRow++)
    {
        for (arrowColumn = 0; arrowColumn < MAX_NUMBER_OF_REGIONS_IN_WIDTH; arrowColumn++) {
            context_ptr->ahdRunningAvgCb[arrowColumn][arrayRow] = 0;
            context_ptr->ahdRunningAvgCr[arrowColumn][arrayRow] = 0;
            context_ptr->ahdRunningAvg[arrowColumn][arrayRow] = 0;
        }
    }


    context_ptr->resetRunningAvg = EB_TRUE;

    context_ptr->isSceneChangeDetected = EB_FALSE;


    return EB_ErrorNone;
}

EbBool SceneTransitionDetector(
    PictureDecisionContext_t *context_ptr,
    SequenceControlSet_t                 *sequence_control_set_ptr,
    PictureParentControlSet_t           **ParentPcsWindow,
    uint32_t                                windowWidthFuture)
{
    PictureParentControlSet_t       *previousPictureControlSetPtr = ParentPcsWindow[0];
    PictureParentControlSet_t       *currentPictureControlSetPtr = ParentPcsWindow[1];
    PictureParentControlSet_t       *futurePictureControlSetPtr = ParentPcsWindow[2];

    // calculating the frame threshold based on the number of 64x64 blocks in the frame
    uint32_t  regionThreshHold;
    uint32_t  regionThreshHoldChroma;
    // this variable determines whether the running average should be reset to equal the ahd or not after detecting a scene change.
    //EbBool resetRunningAvg = context_ptr->resetRunningAvg;

    EbBool isAbruptChange; // this variable signals an abrubt change (scene change or flash)
    EbBool isSceneChange; // this variable signals a frame representing a scene change
    EbBool isFlash; // this variable signals a frame that contains a flash
    EbBool isFade; // this variable signals a frame that contains a fade
    EbBool gradualChange; // this signals the detection of a light scene change a small/localized flash or the start of a fade

    uint32_t  ahd; // accumulative histogram (absolute) differences between the past and current frame

    uint32_t  ahdCb;
    uint32_t  ahdCr;

    uint32_t  ahdErrorCb = 0;
    uint32_t  ahdErrorCr = 0;

    uint32_t **ahdRunningAvgCb = context_ptr->ahdRunningAvgCb;
    uint32_t **ahdRunningAvgCr = context_ptr->ahdRunningAvgCr;
    uint32_t **ahdRunningAvg = context_ptr->ahdRunningAvg;

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
            isSceneChange = EB_FALSE;
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

            if (context_ptr->resetRunningAvg) {
                ahdRunningAvg[regionInPictureWidthIndex][regionInPictureHeightIndex] = ahd;
                ahdRunningAvgCb[regionInPictureWidthIndex][regionInPictureHeightIndex] = ahdCb;
                ahdRunningAvgCr[regionInPictureWidthIndex][regionInPictureHeightIndex] = ahdCr;
            }

            ahdError = ABS((int32_t)ahdRunningAvg[regionInPictureWidthIndex][regionInPictureHeightIndex] - (int32_t)ahd);
            ahdErrorCb = ABS((int32_t)ahdRunningAvgCb[regionInPictureWidthIndex][regionInPictureHeightIndex] - (int32_t)ahdCb);
            ahdErrorCr = ABS((int32_t)ahdRunningAvgCr[regionInPictureWidthIndex][regionInPictureHeightIndex] - (int32_t)ahdCr);


            if ((ahdError > regionThreshHold       && ahd >= ahdError) ||
                (ahdErrorCb > regionThreshHoldChroma && ahdCb >= ahdErrorCb) ||
                (ahdErrorCr > regionThreshHoldChroma && ahdCr >= ahdErrorCr)) {

                isAbruptChange = EB_TRUE;

            }
            else if ((ahdError > (regionThreshHold >> 1)) && ahd >= ahdError) {
                gradualChange = EB_TRUE;
            }

            if (isAbruptChange)
            {
                aidFuturePast = (uint8_t)ABS((int16_t)futurePictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] - (int16_t)previousPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0]);
                aidFuturePresent = (uint8_t)ABS((int16_t)futurePictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] - (int16_t)currentPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0]);
                aidPresentPast = (uint8_t)ABS((int16_t)currentPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] - (int16_t)previousPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0]);

                if (aidFuturePast < FLASH_TH && aidFuturePresent >= FLASH_TH && aidPresentPast >= FLASH_TH) {
                    isFlash = EB_TRUE;
                    //printf ("\nFlash in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                }
                else if (aidFuturePresent < FADE_TH && aidPresentPast < FADE_TH) {
                    isFade = EB_TRUE;
                    //printf ("\nFlash in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                }
                else {
                    isSceneChange = EB_TRUE;
                    //printf ("\nScene Change in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                }

            }
            else if (gradualChange) {

                aidFuturePast = (uint8_t)ABS((int16_t)futurePictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] - (int16_t)previousPictureControlSetPtr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0]);
                if (aidFuturePast < FLASH_TH) {
                    // proper action to be signalled
                    //printf ("\nLight Flash in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                    ahdRunningAvg[regionInPictureWidthIndex][regionInPictureHeightIndex] = (3 * ahdRunningAvg[regionInPictureWidthIndex][regionInPictureHeightIndex] + ahd) / 4;
                }
                else {
                    // proper action to be signalled
                    //printf ("\nLight Scene Change / fade detected in frame# %i , %i\n", currentPictureControlSetPtr->picture_number,aidFuturePast);
                    ahdRunningAvg[regionInPictureWidthIndex][regionInPictureHeightIndex] = (3 * ahdRunningAvg[regionInPictureWidthIndex][regionInPictureHeightIndex] + ahd) / 4;
                }

            }
            else {
                ahdRunningAvg[regionInPictureWidthIndex][regionInPictureHeightIndex] = (3 * ahdRunningAvg[regionInPictureWidthIndex][regionInPictureHeightIndex] + ahd) / 4;
            }

            isAbruptChangeCount += isAbruptChange;
            isSceneChangeCount += isSceneChange;
        }
    }

    (void)windowWidthFuture;
    (void)isFlash;
    (void)isFade;

    if (isAbruptChangeCount >= regionCountThreshold) {
        context_ptr->resetRunningAvg = EB_TRUE;
    }
    else {
        context_ptr->resetRunningAvg = EB_FALSE;
    }

    if ((isSceneChangeCount >= regionCountThreshold) && ((!ParentPcsWindow[1]->fade_in_to_black) && (!ParentPcsWindow[1]->fade_out_from_black))) {
        return(EB_TRUE);
    }
    else {
        return(EB_FALSE);
    }

}

/***************************************************************************************************
* ReleasePrevPictureFromReorderQueue
***************************************************************************************************/
EbErrorType ReleasePrevPictureFromReorderQueue(
    EncodeContext_t                 *encode_context_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    PictureDecisionReorderEntry_t   *queuePreviousEntryPtr;
    int32_t                           previousEntryIndex;


    // Get the previous entry from the Picture Decision Reordering Queue (Entry N-1)
    // P.S. The previous entry in display order is needed for Scene Change Detection
    previousEntryIndex = (encode_context_ptr->picture_decision_reorder_queue_head_index == 0) ? PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1 : encode_context_ptr->picture_decision_reorder_queue_head_index - 1;
    queuePreviousEntryPtr = encode_context_ptr->picture_decision_reorder_queue[previousEntryIndex];

    // SB activity classification based on (0,0) SAD & picture activity derivation
    if (queuePreviousEntryPtr->parentPcsWrapperPtr) {

        // Reset the Picture Decision Reordering Queue Entry
        // P.S. The reset of the Picture Decision Reordering Queue Entry could not be done before running the Scene Change Detector
        queuePreviousEntryPtr->picture_number += PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH;
        queuePreviousEntryPtr->parentPcsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;
    }

    return return_error;
}
#if NEW_PRED_STRUCT

/***************************************************************************************************
* Initializes mini GOP activity array
*
***************************************************************************************************/
EbErrorType initialize_mini_gop_activity_array(
    PictureDecisionContext_t        *context_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    uint32_t mini_gop_index;

    // Loop over all mini GOPs
    for (mini_gop_index = 0; mini_gop_index < MINI_GOP_MAX_COUNT; ++mini_gop_index) {

        context_ptr->miniGopActivityArray[mini_gop_index] = (GetMiniGopStats(mini_gop_index)->hierarchical_levels == MIN_HIERARCHICAL_LEVEL) ?
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
    PictureDecisionContext_t        *context_ptr,
    EncodeContext_t                 *encode_context_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    uint32_t	mini_gop_index;

    context_ptr->totalNumberOfMiniGops = 0;

    // Loop over all mini GOPs
    mini_gop_index = 0;
    while (mini_gop_index < MINI_GOP_MAX_COUNT) {

        // Only for a valid mini GOP
        if (GetMiniGopStats(mini_gop_index)->endIndex < encode_context_ptr->pre_assignment_buffer_count && context_ptr->miniGopActivityArray[mini_gop_index] == EB_FALSE) {

            context_ptr->miniGopStartIndex[context_ptr->totalNumberOfMiniGops] = GetMiniGopStats(mini_gop_index)->startIndex;
            context_ptr->miniGopEndIndex[context_ptr->totalNumberOfMiniGops] = GetMiniGopStats(mini_gop_index)->endIndex;
            context_ptr->miniGopLength[context_ptr->totalNumberOfMiniGops] = GetMiniGopStats(mini_gop_index)->lenght;
            context_ptr->miniGopHierarchicalLevels[context_ptr->totalNumberOfMiniGops] = GetMiniGopStats(mini_gop_index)->hierarchical_levels;
            context_ptr->miniGopIntraCount[context_ptr->totalNumberOfMiniGops] = 0;
            context_ptr->miniGopIdrCount[context_ptr->totalNumberOfMiniGops] = 0;

            context_ptr->totalNumberOfMiniGops++;
        }

        mini_gop_index += context_ptr->miniGopActivityArray[mini_gop_index] ?
            1 :
            MiniGopOffset[GetMiniGopStats(mini_gop_index)->hierarchical_levels - MIN_HIERARCHICAL_LEVEL];

    }

    // Only in presence of at least 1 valid mini GOP
    if (context_ptr->totalNumberOfMiniGops != 0) {
        context_ptr->miniGopIntraCount[context_ptr->totalNumberOfMiniGops - 1] = encode_context_ptr->pre_assignment_buffer_intra_count;
        context_ptr->miniGopIdrCount[context_ptr->totalNumberOfMiniGops - 1] = encode_context_ptr->pre_assignment_buffer_idr_count;
    }

    return return_error;
}

/***************************************************************************************************
* Handles an incomplete picture window map
*
*
***************************************************************************************************/
EbErrorType handle_incomplete_picture_window_map(
    PictureDecisionContext_t        *context_ptr,
    EncodeContext_t                 *encode_context_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    if (context_ptr->totalNumberOfMiniGops == 0) {

        context_ptr->miniGopStartIndex[context_ptr->totalNumberOfMiniGops] = 0;
        context_ptr->miniGopEndIndex[context_ptr->totalNumberOfMiniGops] = encode_context_ptr->pre_assignment_buffer_count - 1;
        context_ptr->miniGopLength[context_ptr->totalNumberOfMiniGops] = encode_context_ptr->pre_assignment_buffer_count - context_ptr->miniGopStartIndex[context_ptr->totalNumberOfMiniGops];
        context_ptr->miniGopHierarchicalLevels[context_ptr->totalNumberOfMiniGops] = 3;// MIN_HIERARCHICAL_LEVEL; // AMIR to be updated after other predictions are supported

        context_ptr->totalNumberOfMiniGops++;

    }
    else if (context_ptr->miniGopEndIndex[context_ptr->totalNumberOfMiniGops - 1] < encode_context_ptr->pre_assignment_buffer_count - 1) {

        context_ptr->miniGopStartIndex[context_ptr->totalNumberOfMiniGops] = context_ptr->miniGopEndIndex[context_ptr->totalNumberOfMiniGops - 1] + 1;
        context_ptr->miniGopEndIndex[context_ptr->totalNumberOfMiniGops] = encode_context_ptr->pre_assignment_buffer_count - 1;
        context_ptr->miniGopLength[context_ptr->totalNumberOfMiniGops] = encode_context_ptr->pre_assignment_buffer_count - context_ptr->miniGopStartIndex[context_ptr->totalNumberOfMiniGops];
        context_ptr->miniGopHierarchicalLevels[context_ptr->totalNumberOfMiniGops] = 3;// MIN_HIERARCHICAL_LEVEL;// AMIR
        context_ptr->miniGopIntraCount[context_ptr->totalNumberOfMiniGops - 1] = 0;
        context_ptr->miniGopIdrCount[context_ptr->totalNumberOfMiniGops - 1] = 0;

        context_ptr->totalNumberOfMiniGops++;
    }

    context_ptr->miniGopIntraCount[context_ptr->totalNumberOfMiniGops - 1] = encode_context_ptr->pre_assignment_buffer_intra_count;
    context_ptr->miniGopIdrCount[context_ptr->totalNumberOfMiniGops - 1] = encode_context_ptr->pre_assignment_buffer_idr_count;

    return return_error;
}
/***************************************************************************************************
* If a switch happens, then update the RPS of the base layer frame separating the 2 different prediction structures
* Clean up the reference queue dependant counts of the base layer frame separating the 2 different prediction structures
*
***************************************************************************************************/
EbErrorType update_base_layer_reference_queue_dependent_count(
    PictureDecisionContext_t        *context_ptr,
    EncodeContext_t                 *encode_context_ptr,
    SequenceControlSet_t            *sequence_control_set_ptr,
    uint32_t                         mini_gop_index) {

    if (!context_ptr || !encode_context_ptr || !sequence_control_set_ptr)
        return EB_ErrorBadParameter;

    EbErrorType return_error = EB_ErrorNone;

    PaReferenceQueueEntry_t         *input_entry_ptr;
    uint32_t                         input_queue_index;

    PredictionStructure_t           *next_pred_struct_ptr;
    PredictionStructureEntry_t      *next_base_layer_pred_position_ptr;

    uint32_t                         dependant_list_positive_entries;
    uint32_t                         dependant_list_removed_entries;
    uint32_t                         dep_list_count;

    uint32_t                         dep_idx;
    uint64_t                         dep_poc;

    PictureParentControlSet_t       *picture_control_set_ptr;

    // Get the 1st PCS mini GOP
    picture_control_set_ptr = (PictureParentControlSet_t*)encode_context_ptr->pre_assignment_buffer[context_ptr->miniGopStartIndex[mini_gop_index]]->objectPtr;

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
                for (dep_idx = 0; dep_idx < input_entry_ptr->list0.listCount; ++dep_idx) {
                    if (input_entry_ptr->list0.list[dep_idx] >= 0) {
                        dependant_list_positive_entries++;
                    }
                }
                input_entry_ptr->list0.listCount = input_entry_ptr->list0.listCount - dependant_list_positive_entries;
                dependant_list_positive_entries = 0;
                for (dep_idx = 0; dep_idx < input_entry_ptr->list1.listCount; ++dep_idx) {
                    if (input_entry_ptr->list1.list[dep_idx] >= 0) {
                        dependant_list_positive_entries++;
                    }
                }
                input_entry_ptr->list1.listCount = input_entry_ptr->list1.listCount - dependant_list_positive_entries;

                // 2nd step: inherit the positive dependant counts of the current mini GOP
                // Get the RPS set of the current mini GOP
                next_pred_struct_ptr = GetPredictionStructure(
                    encode_context_ptr->prediction_structure_group_ptr,
                    picture_control_set_ptr->pred_structure,
                    1,
                    picture_control_set_ptr->hierarchical_levels);			// Number of temporal layer in the current mini GOP  

                // Get the RPS of a base layer input
                next_base_layer_pred_position_ptr = next_pred_struct_ptr->predStructEntryPtrArray[next_pred_struct_ptr->predStructEntryCount - 1];

                for (dep_idx = 0; dep_idx < next_base_layer_pred_position_ptr->depList0.listCount; ++dep_idx) {
                    if (next_base_layer_pred_position_ptr->depList0.list[dep_idx] >= 0) {
                        input_entry_ptr->list0.list[input_entry_ptr->list0.listCount++] = next_base_layer_pred_position_ptr->depList0.list[dep_idx];
                    }
                }


                for (dep_idx = 0; dep_idx < next_base_layer_pred_position_ptr->depList1.listCount; ++dep_idx) {
                    if (next_base_layer_pred_position_ptr->depList1.list[dep_idx] >= 0) {
                        input_entry_ptr->list1.list[input_entry_ptr->list1.listCount++] = next_base_layer_pred_position_ptr->depList1.list[dep_idx];
                    }
                }

                // 3rd step: update the dependant count
                dependant_list_removed_entries = input_entry_ptr->depList0Count + input_entry_ptr->depList1Count - input_entry_ptr->dependentCount;
                input_entry_ptr->depList0Count = input_entry_ptr->list0.listCount;
                input_entry_ptr->depList1Count = input_entry_ptr->list1.listCount;
                input_entry_ptr->dependentCount = input_entry_ptr->depList0Count + input_entry_ptr->depList1Count - dependant_list_removed_entries;

            }
            else {
                // Modify Dependent List0
                dep_list_count = input_entry_ptr->list0.listCount;
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
                        --input_entry_ptr->dependentCount;
                        CHECK_REPORT_ERROR(
                            (input_entry_ptr->dependentCount != ~0u),
                            encode_context_ptr->app_callback_ptr,
                            EB_ENC_PD_ERROR3);
                    }
                }
                // Modify Dependent List1
                dep_list_count = input_entry_ptr->list1.listCount;
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
                        --input_entry_ptr->dependentCount;

                        CHECK_REPORT_ERROR(
                            (input_entry_ptr->dependentCount != ~0u),
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
    PictureDecisionContext_t        *context_ptr,
    uint32_t                         mini_gop_index,
    uint32_t				        picture_index) {

    if ((context_ptr->miniGopHierarchicalLevels[mini_gop_index] == 4 && context_ptr->miniGopLength[mini_gop_index] == 16 && (picture_index == 7 || picture_index == 23)) ||	// supposedly a 4L reference frame for 5L prediction structure 
        (context_ptr->miniGopHierarchicalLevels[mini_gop_index] == 5 && context_ptr->miniGopLength[mini_gop_index] == 32 && (picture_index == 7 || picture_index == 23))) { // supposedly a 4L reference frame for 6L prediction structure
        return(EB_TRUE);
    }
    else {
        return(EB_FALSE);
    }
}

#endif

/***************************************************************************************************
* Generates mini GOP RPSs
*
*
***************************************************************************************************/
EbErrorType GenerateMiniGopRps(
    PictureDecisionContext_t        *context_ptr,
    EncodeContext_t                 *encode_context_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    uint32_t                         miniGopIndex;
    PictureParentControlSet_t    *picture_control_set_ptr;
    uint32_t                         pictureIndex;

    // Loop over all mini GOPs
    for (miniGopIndex = 0; miniGopIndex < context_ptr->totalNumberOfMiniGops; ++miniGopIndex) {

        // Loop over picture within the mini GOP
        for (pictureIndex = context_ptr->miniGopStartIndex[miniGopIndex]; pictureIndex <= context_ptr->miniGopEndIndex[miniGopIndex]; pictureIndex++) {

            picture_control_set_ptr = (PictureParentControlSet_t*)encode_context_ptr->pre_assignment_buffer[pictureIndex]->objectPtr;

            picture_control_set_ptr->pred_structure = EB_PRED_RANDOM_ACCESS;
            picture_control_set_ptr->hierarchical_levels = (uint8_t)context_ptr->miniGopHierarchicalLevels[miniGopIndex];

            picture_control_set_ptr->pred_struct_ptr = GetPredictionStructure(
                encode_context_ptr->prediction_structure_group_ptr,
                picture_control_set_ptr->pred_structure,
                1,
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
    PictureParentControlSet_t   *picture_control_set_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    //  MDC Partitioning Method              Settings
    //  PIC_ALL_DEPTH_MODE                   ALL sq and nsq: SB size -> 4x4
    //  PIC_ALL_C_DEPTH_MODE                 ALL sq and nsq: SB size -> 4x4  (No 4xN ; Nx4)
    //  PIC_SQ_DEPTH_MODE                    ONLY sq: SB size -> 4x4
    //  PIC_SQ_NON4_DEPTH_MODE               ONLY sq: SB size -> 8x8  (No 4x4)

    if (picture_control_set_ptr->enc_mode == ENC_M0)
        picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;

    else if (picture_control_set_ptr->enc_mode == ENC_M1) {
#if TUNED_SETTINGS_FOR_M1
        picture_control_set_ptr->pic_depth_mode = PIC_ALL_C_DEPTH_MODE;
#else
        if (picture_control_set_ptr->is_used_as_reference_flag == EB_TRUE)
            picture_control_set_ptr->pic_depth_mode = PIC_ALL_DEPTH_MODE;
        else
            picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
#endif
    }
    else if (picture_control_set_ptr->enc_mode == ENC_M2) {
        if (picture_control_set_ptr->is_used_as_reference_flag == EB_TRUE)
            picture_control_set_ptr->pic_depth_mode = PIC_ALL_C_DEPTH_MODE;
        else
            picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
    }
    else {
        if (picture_control_set_ptr->is_used_as_reference_flag == EB_TRUE)
            picture_control_set_ptr->pic_depth_mode = PIC_SQ_DEPTH_MODE;
        else
            picture_control_set_ptr->pic_depth_mode = PIC_SQ_NON4_DEPTH_MODE;
    }

    picture_control_set_ptr->max_number_of_pus_per_sb = (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) ? MAX_ME_PU_COUNT : SQUARE_PU_COUNT;
#if NSQ_SEARCH_LEVELS
    // NSQ search Level                               Settings
    // 0                                              OFF
    // 1                                              Allow only NSQ Intra-FULL if parent SQ is intra-coded and vice versa.
    // 2                                              Allow only NSQ Inter-NEAREST/NEAR/GLOBAL if parent SQ has no coeff
    // 3                                              Allow only NSQ Intra-FULL and Inter-NEWMV if parent SQ is NEWMV
    // 4                                              Allow only NSQ Inter-FULL and Intra-Z3 if parent SQ is intra-coded
    // 5                                              Allow NSQ Intra-FULL and Inter-FULL
#if TUNED_SETTINGS_FOR_M0 || TUNED_SETTINGS_FOR_M1
    if (!MR_MODE)
    picture_control_set_ptr->nsq_search_level        = NSQ_SEARCH_BASE_ON_SQ_COEFF;
    else
#endif
    picture_control_set_ptr->nsq_search_level        = NSQ_SEARCH_FULL;


    if (picture_control_set_ptr->nsq_search_level == NSQ_SEARCH_OFF) {
        if (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) picture_control_set_ptr->pic_depth_mode = PIC_SQ_DEPTH_MODE;
    }
    if (picture_control_set_ptr->pic_depth_mode > PIC_SQ_DEPTH_MODE) {
        assert(picture_control_set_ptr->nsq_search_level != NSQ_SEARCH_OFF);
    }
#endif
#if INTERPOLATION_SEARCH_LEVELS
    // Interpolation search Level                     Settings
    // 0                                              OFF
    // 1                                              Interpolation search at inter-depth
    // 2                                              Interpolation search at full loop
    // 3                                              Interpolation search at fast loop
    if (picture_control_set_ptr->enc_mode == ENC_M0) {
        picture_control_set_ptr->interpolation_search_level = IT_SEARCH_FAST_LOOP;
    }
    else {
        picture_control_set_ptr->interpolation_search_level = IT_SEARCH_OFF;
    }
#else
    // Interpolation filter search Level MD         Settings
    // 0                                            OFF
    // 1                                            FAST
    // 2                                            ON
    if (picture_control_set_ptr->enc_mode == ENC_M0)
        picture_control_set_ptr->interpolation_filter_search_mode = 1;
    else
        picture_control_set_ptr->interpolation_filter_search_mode = 0;

#endif

    // Loop filter Level                            Settings
    // 0                                            OFF
    // 1                                            CU-BASED
    // 2                                            LIGHT FRAME-BASED
    // 3                                            FULL FRAME-BASED

    if (!picture_control_set_ptr->sequence_control_set_ptr->static_config.disable_dlf_flag) {
        if (picture_control_set_ptr->enc_mode >= ENC_M2)
            picture_control_set_ptr->loop_filter_mode = 1;
        else  if (picture_control_set_ptr->enc_mode == ENC_M1)
#if TUNED_SETTINGS_FOR_M1
            picture_control_set_ptr->loop_filter_mode = 3;
#else
            picture_control_set_ptr->loop_filter_mode = 2;
#endif
        else  if (picture_control_set_ptr->enc_mode == ENC_M0)
            picture_control_set_ptr->loop_filter_mode = 3;
    }
    else {
        picture_control_set_ptr->loop_filter_mode = 0;
    }
#if FAST_CDEF
    // CDEF Level                                   Settings
    // 0                                            OFF
    // 1                                            4 step refinement
    // 2                                            8 step refinement
    // 3                                            16 step refinement
    SequenceControlSet_t                    *sequence_control_set_ptr;
    sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    if (sequence_control_set_ptr->enable_cdef) {
        if (picture_control_set_ptr->enc_mode >= ENC_M3)
            picture_control_set_ptr->cdef_filter_mode = 1;
        else  if (picture_control_set_ptr->enc_mode == ENC_M2)
            picture_control_set_ptr->cdef_filter_mode = 2;
        else  if (picture_control_set_ptr->enc_mode <= ENC_M1)
            picture_control_set_ptr->cdef_filter_mode = 3;
    }
    else {
        picture_control_set_ptr->cdef_filter_mode = 0;
    }
#endif
#if FAST_SG
    // SG Level                                    Settings
    // 0                                            OFF
    // 1                                            0 step refinement
    // 2                                            1 step refinement
    // 3                                            4 step refinement
    // 4                                            16 step refinement

    Av1Common* cm = picture_control_set_ptr->av1_cm;

    if (picture_control_set_ptr->enc_mode >= ENC_M3)
        cm->sg_filter_mode = 1;
    else  if (picture_control_set_ptr->enc_mode == ENC_M2)
        cm->sg_filter_mode = 2;
#if TUNED_SETTINGS_FOR_M0
    else  if (picture_control_set_ptr->enc_mode == ENC_M1)
        cm->sg_filter_mode = 3;
    else  if (picture_control_set_ptr->enc_mode == ENC_M0)
        cm->sg_filter_mode = 4;
#else
    else  if (picture_control_set_ptr->enc_mode <= ENC_M1)
        cm->sg_filter_mode = 3;
#endif

#endif

#if FAST_WN
    // WN Level                                     Settings
    // 1                                            5-Tap luma/ 5-Tap chroma
    // 2                                            7-Tap luma/ 5-Tap chroma

    if (picture_control_set_ptr->enc_mode >= ENC_M1)
        cm->wn_filter_mode = 1;
    else
        cm->wn_filter_mode = 2;
#endif

#if TX_SEARCH_LEVELS
    // Tx_search Level                                Settings
    // 0                                              OFF
    // 1                                              Tx search at encdec
    // 2                                              Tx search at inter-depth
    // 3                                              Tx search at full loop

    if (picture_control_set_ptr->enc_mode > ENC_M1) {
        picture_control_set_ptr->tx_search_level = TX_SEARCH_ENC_DEC;
    }
    else {
        picture_control_set_ptr->tx_search_level = TX_SEARCH_FULL_LOOP;
    }

    // Set tx search skip weights (MAX_MODE_COST: no skipping; 0: always skipping)
#if TUNED_SETTINGS_FOR_M0
    if (!MR_MODE && picture_control_set_ptr->enc_mode == ENC_M0)
        picture_control_set_ptr->tx_weight = FC_SKIP_TX_SR_TH;
    else
#endif
#if TUNED_SETTINGS_FOR_M1
    if (!MR_MODE && picture_control_set_ptr->enc_mode == ENC_M1)
        picture_control_set_ptr->tx_weight = FC_SKIP_TX_SR_TH_M1;
    else
#endif
        picture_control_set_ptr->tx_weight = MAX_MODE_COST;

    // Set tx search reduced set falg (0: full tx set; 1: reduced tx set)
#if TUNED_SETTINGS_FOR_M1
    if (picture_control_set_ptr->enc_mode == ENC_M2) {
#else
    if (picture_control_set_ptr->enc_mode == ENC_M1) {
#endif
        picture_control_set_ptr->tx_search_reduced_set = 1;
    }
    else {
        picture_control_set_ptr->tx_search_reduced_set = 0;
    }
#endif

    // Intra prediction mode                       Settings
    // 0                                            OFF : disable_angle_prediction
    // 1                                            OFF per block : disable_angle_prediction for 64/32/4
    // 2                                            LIGHT: disable_z2_prediction && disable_angle_refinement
    // 3                                            LIGHT per block : disable_z2_prediction && disable_angle_refinement  for 64/32/4
    // 4                                            FULL  

    // The intra prediction level is a combination of the intra prediction modes.

    // Intra prediction levels                      Settings
    // 0                                            OFF : disable_angle_prediction
    // 1                                            Disable_angle_prediction for 64/32/4 (mode 1) @ BASE AND OFF (mode 0) Otherwise
    // 2                                            Disable_z2_prediction && disable_angle_refinement for 64/32/4 (mode 3) @ BASE AND OFF (mode 0) Otherwise
    // 3                                            Full (mode 4) @ BASE AND Disable_z2_prediction && disable_angle_refinement (mode 2) Otherwise
    // 4                                            FULL 

    uint8_t intra_pred_level = 4;

    switch (picture_control_set_ptr->enc_mode) {
    case 0:
        intra_pred_level = 3; //ENC_M0
        break;
    case 1:
        intra_pred_level = 2; //ENC_M1
        break;
    case 2:
        intra_pred_level = 2; //ENC_M2
        break;
    case 3:
        intra_pred_level = 2; //ENC_M3
        break;
    default:
        intra_pred_level = 4; //MR_MODE
        break;
    }

    if (intra_pred_level == 4) {
        picture_control_set_ptr->intra_pred_mode = 4;

    }else if (intra_pred_level == 3) { 
        if (picture_control_set_ptr->temporal_layer_index == 0)
            picture_control_set_ptr->intra_pred_mode = 4;
        else
            picture_control_set_ptr->intra_pred_mode = 2;
    }else if (intra_pred_level == 2) {  
        if (picture_control_set_ptr->temporal_layer_index == 0)
            picture_control_set_ptr->intra_pred_mode = 3;
        else
            picture_control_set_ptr->intra_pred_mode = 0;
    }
    else if (intra_pred_level == 1) { 
        if (picture_control_set_ptr->temporal_layer_index == 0)
            picture_control_set_ptr->intra_pred_mode = 1;
        else
            picture_control_set_ptr->intra_pred_mode = 0;
    }

    return return_error;
}

/***************************************************************************
* Set the default subPel enble/disable flag for each frame
****************************************************************************/
uint8_t PictureLevelSubPelSettings(
    uint8_t   input_resolution,
    uint8_t   enc_mode,
    uint8_t   temporal_layer_index,
    EbBool    is_used_as_reference_flag){

    // Set Subpel Flag
    uint8_t subPelMode = 0;
#if ENCODER_MODE_CLEANUP
    UNUSED(input_resolution);
    UNUSED(enc_mode);
    UNUSED(temporal_layer_index);
    UNUSED(is_used_as_reference_flag);
    subPelMode =  1;
#else
    if (input_resolution >= INPUT_SIZE_4K_RANGE) {
        subPelMode = (enc_mode <= ENC_M1) ? 1 : 0;
    }
    else {

        if (enc_mode <= ENC_M2) {

            subPelMode = 1;

        }
        else if (enc_mode <= ENC_M3) {

            subPelMode = is_used_as_reference_flag ? 1 : 0;
        }

        else if (enc_mode == ENC_M5) {

            subPelMode = temporal_layer_index == 0 ? 1 : 0;
        }
        else {
            subPelMode = 0;
        }
    }
#endif
    return subPelMode;
}
#if !CHROMA_BLIND
/***************************************************************************
* Set the default chroma mode for each frame
****************************************************************************/
EbChromaMode PictureLevelChromaSettings(
    uint8_t   input_resolution,
    uint8_t   enc_mode,
    uint8_t   slice_type,
    uint8_t   temporal_layer_index,
    EbBool is_used_as_reference_flag){

    EbChromaMode chroma_mode = CHROMA_MODE_FULL;
    UNUSED(input_resolution);
    UNUSED(enc_mode);
    UNUSED(slice_type);
    UNUSED(temporal_layer_index);
    UNUSED(is_used_as_reference_flag);

#if !ENCODER_MODE_CLEANUP
    if ((enc_mode >= ENC_M3 && input_resolution >= INPUT_SIZE_4K_RANGE) || enc_mode > ENC_M3)
    {
        if (enc_mode == ENC_M6 && input_resolution >= INPUT_SIZE_4K_RANGE)
            chroma_mode = CHROMA_MODE_BEST;
        else
            chroma_mode = (temporal_layer_index > 0 || slice_type == I_SLICE) ? CHROMA_MODE_BEST : CHROMA_MODE_FULL;
    }

    if (enc_mode == ENC_M3 || enc_mode == ENC_M4)
        chroma_mode = (is_used_as_reference_flag) ? CHROMA_MODE_FULL : chroma_mode;
#endif
    return chroma_mode;
}
#endif

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
*Note: the  SceneChange I has predType = EB_PRED_RANDOM_ACCESS. if SChange is aligned on the miniGop,
we do not break the GOP.
*************************************************/
void  Av1GenerateRpsInfo(
    PictureParentControlSet_t       *picture_control_set_ptr,
    EncodeContext_t                 *encode_context_ptr,
    PictureDecisionContext_t        *context_ptr,
    uint32_t                           pictureIndex
)
{
    (void)encode_context_ptr;
    Av1RpsNode_t *av1Rps = &picture_control_set_ptr->av1RefSignal;

    //set Frame Type
    if (picture_control_set_ptr->slice_type == I_SLICE)
        picture_control_set_ptr->av1FrameType = picture_control_set_ptr->idr_flag ? KEY_FRAME : INTRA_ONLY_FRAME;
    else
        picture_control_set_ptr->av1FrameType = INTER_FRAME;


    picture_control_set_ptr->intra_only = picture_control_set_ptr->slice_type == I_SLICE ? 1 : 0;

    //RPS for Flat GOP
    if (picture_control_set_ptr->hierarchical_levels == 0)
    {

        memset(av1Rps->refDpbIndex, 0, 7);
        av1Rps->refreshFrameMask = 1;
        picture_control_set_ptr->showFrame = EB_TRUE;
        picture_control_set_ptr->hasShowExisting = EB_FALSE;

    }
    else if (picture_control_set_ptr->hierarchical_levels == 3)//RPS for 4L GOP
    {

        //Reset miniGop Toggling. The first miniGop after a KEY frame has toggle=0
        if (picture_control_set_ptr->av1FrameType == KEY_FRAME)
        {
            context_ptr->miniGopToggle = 0;
            picture_control_set_ptr->showFrame = EB_TRUE;
            picture_control_set_ptr->hasShowExisting = EB_FALSE;
            return;
        }

        //pictureIndex has this order:
        //        0      2    4      6
        //            1          5
        //                 3
        //                             8(could be an I)


        //DPB: Loc7|Loc6|Loc5|Loc4|Loc3|Loc2|Loc1|Loc0
        //Layer 0 : toggling bwteween DPB Location 0, and  locations 3-4-5-6-7
        //Layer 1 : DPB Location 1
        //Layer 2 : DPB Location 2


        //         1     3    5      7
        //            2          6
        //                 4
        //base0:0                      base1:8
        const uint8_t  base0_idx = context_ptr->miniGopToggle ? 0 : 3; //Base layer for prediction from past
        const uint8_t  base1_idx = context_ptr->miniGopToggle ? 3 : 0; //Base layer for prediction from future
        const uint8_t  layer1_idx = 1;
        const uint8_t  layer2_idx = 2;


        switch (picture_control_set_ptr->temporal_layer_index) {

        case 0:

            av1Rps->refDpbIndex[0] = base0_idx;
            av1Rps->refDpbIndex[6] = base0_idx;
            av1Rps->refreshFrameMask = context_ptr->miniGopToggle ? 248 : 1;
            break;
        case 1:
            av1Rps->refDpbIndex[0] = base0_idx;
            av1Rps->refDpbIndex[6] = base1_idx;
            av1Rps->refreshFrameMask = 2;
            break;
        case 2:

            if (pictureIndex == 1) {
                av1Rps->refDpbIndex[0] = base0_idx;
                av1Rps->refDpbIndex[6] = layer1_idx;
            }
            else if (pictureIndex == 5) {
                av1Rps->refDpbIndex[0] = layer1_idx;
                av1Rps->refDpbIndex[6] = base1_idx;
            }
            else {
                printf("Error in GOp indexing\n");
            }
            av1Rps->refreshFrameMask = 4;
            break;
        case 3:
            if (pictureIndex == 0) {
                av1Rps->refDpbIndex[0] = base0_idx;
                av1Rps->refDpbIndex[6] = layer2_idx;
            }
            else if (pictureIndex == 2) {
                av1Rps->refDpbIndex[0] = layer2_idx;
                av1Rps->refDpbIndex[6] = layer1_idx;
            }
            else if (pictureIndex == 4) {
                av1Rps->refDpbIndex[0] = layer1_idx;
                av1Rps->refDpbIndex[6] = layer2_idx;
            }
            else if (pictureIndex == 6) {
                av1Rps->refDpbIndex[0] = layer2_idx;
                av1Rps->refDpbIndex[6] = base1_idx;
            }
            else {
                printf("Error in GOp indexing\n");
            }
            av1Rps->refreshFrameMask = 0;
            break;
        default:
            printf("Error: unexpected picture mini Gop number\n");
            break;
        }

        if (picture_control_set_ptr->pred_struct_ptr->predType == EB_PRED_LOW_DELAY_P)
        {
            //P frames.
            av1Rps->refDpbIndex[1] = av1Rps->refDpbIndex[2] = av1Rps->refDpbIndex[3] = av1Rps->refDpbIndex[0];
            av1Rps->refDpbIndex[4] = av1Rps->refDpbIndex[5] = av1Rps->refDpbIndex[6] = av1Rps->refDpbIndex[0];
            picture_control_set_ptr->showFrame = EB_TRUE;
            picture_control_set_ptr->hasShowExisting = EB_FALSE;
        }
        else if (picture_control_set_ptr->pred_struct_ptr->predType == EB_PRED_RANDOM_ACCESS)
        {
            av1Rps->refDpbIndex[1] = av1Rps->refDpbIndex[2] = av1Rps->refDpbIndex[3] = av1Rps->refDpbIndex[0];
            av1Rps->refDpbIndex[4] = av1Rps->refDpbIndex[5] = av1Rps->refDpbIndex[6];

            //Decide on Show Mecanism
            if (picture_control_set_ptr->slice_type == I_SLICE)
            {
                //3 cases for I slice:  1:Key Frame treated above.  2: broken MiniGop due to sc or intra refresh  3: complete miniGop due to sc or intra refresh
                if (context_ptr->miniGopLength[0] < picture_control_set_ptr->pred_struct_ptr->predStructPeriod)
                {
                    //Scene Change that breaks the mini gop and switch to LDP (if I scene change happens to be aligned with a complete miniGop, then we do not break the pred structure)
                    picture_control_set_ptr->showFrame = EB_TRUE;
                    picture_control_set_ptr->hasShowExisting = EB_FALSE;
                }
                else
                {
                    picture_control_set_ptr->showFrame = EB_FALSE;
                    picture_control_set_ptr->hasShowExisting = EB_FALSE;
                }
            }
            else//B pic
            {
                if (context_ptr->miniGopLength[0] != picture_control_set_ptr->pred_struct_ptr->predStructPeriod)
                    printf("Error in GOp indexing3\n");

                if (picture_control_set_ptr->is_used_as_reference_flag)
                {
                    picture_control_set_ptr->showFrame = EB_FALSE;
                    picture_control_set_ptr->hasShowExisting = EB_FALSE;
                }
                else
                {
                    picture_control_set_ptr->showFrame = EB_TRUE;
                    picture_control_set_ptr->hasShowExisting = EB_TRUE;

                    if (pictureIndex == 0) {
                        picture_control_set_ptr->showExistingLoc = layer2_idx;
                    }
                    else if (pictureIndex == 2) {
                        picture_control_set_ptr->showExistingLoc = layer1_idx;
                    }
                    else if (pictureIndex == 4) {
                        picture_control_set_ptr->showExistingLoc = layer2_idx;
                    }
                    else if (pictureIndex == 6) {
                        picture_control_set_ptr->showExistingLoc = base1_idx;
                    }
                    else {
                        printf("Error in GOp indexing2\n");
                    }

                }

            }

        }
        else {
            printf("Error: Not supported GOP structure!");
            exit(0);
        }

        //last pic in MiniGop: mGop Toggling
        //mini GOP toggling since last Key Frame.
        //a regular I keeps the toggling process and does not reset the toggle.  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
        if (pictureIndex == context_ptr->miniGopEndIndex[0])
            context_ptr->miniGopToggle = 1 - context_ptr->miniGopToggle;
    }
#if NEW_PRED_STRUCT
    else if (picture_control_set_ptr->hierarchical_levels == 4)//RPS for 4L GOP
    {

    //Reset miniGop Toggling. The first miniGop after a KEY frame has toggle=0
    if (picture_control_set_ptr->av1FrameType == KEY_FRAME)
    {
        context_ptr->miniGopToggle = 0;
        picture_control_set_ptr->showFrame = EB_TRUE;
        picture_control_set_ptr->hasShowExisting = EB_FALSE;
        return;
    }


    //         0     2    4      6    8     10     12      14
    //            1          5           9            13
    //                 3                        11
    //                              7

    //DPB: Loc7|Loc6|Loc5|Loc4|Loc3|Loc2|Loc1|Loc0
    //Layer 0 : toggling bwteween DPB Location 0, and  locations 3-4-5-6-7
    //Layer 1 : DPB Location 1
    //Layer 2 : DPB Location 2
    //Layer 3 : DPB Location 3

    //         1     3    5      7    9     11     13      15
    //            2          6           10            14
    //                 4                        12
    //                              8
    //base0:0                                               base1:16
    const uint8_t  base0_idx = context_ptr->miniGopToggle ? 0 : 3; //Base layer for prediction from past
    const uint8_t  base1_idx = context_ptr->miniGopToggle ? 3 : 0; //Base layer for prediction from future
    const uint8_t  layer1_idx = 1;
    const uint8_t  layer2_idx = 2;
    const uint8_t  layer3_idx1 = 4;
    const uint8_t  layer3_idx2 = 5;

    switch (picture_control_set_ptr->temporal_layer_index) {

    case 0:

        av1Rps->refDpbIndex[0] = base0_idx;
        av1Rps->refDpbIndex[6] = base0_idx;
        av1Rps->refreshFrameMask = context_ptr->miniGopToggle ? 200 : 1;
        break;
    case 1:
        av1Rps->refDpbIndex[0] = base0_idx;
        av1Rps->refDpbIndex[6] = base1_idx;
        av1Rps->refreshFrameMask = 2;
        break;
    case 2:

        if (pictureIndex == 3) {
            av1Rps->refDpbIndex[0] = base0_idx;
            av1Rps->refDpbIndex[6] = layer1_idx;
    }
        else if (pictureIndex == 11) {
            av1Rps->refDpbIndex[0] = layer1_idx;
            av1Rps->refDpbIndex[6] = base1_idx;
        }
        av1Rps->refreshFrameMask = 4;
        break;
    case 3:

        if (pictureIndex == 1) {
            av1Rps->refDpbIndex[0] = base0_idx;
            av1Rps->refDpbIndex[6] = layer2_idx;
            av1Rps->refreshFrameMask = 16;
        }
        else if (pictureIndex == 5) {
            av1Rps->refDpbIndex[0] = layer2_idx;
            av1Rps->refDpbIndex[6] = layer1_idx;
            av1Rps->refreshFrameMask = 32;
        }
        else if (pictureIndex == 9) {
            av1Rps->refDpbIndex[0] = layer1_idx;
            av1Rps->refDpbIndex[6] = layer2_idx;
            av1Rps->refreshFrameMask = 16;
        }
        else if (pictureIndex == 13) {
            av1Rps->refDpbIndex[0] = layer2_idx;
            av1Rps->refDpbIndex[6] = base1_idx;
            av1Rps->refreshFrameMask = 32;
        }
        else {
            printf("Error in GOp indexing\n");
        }
        break;
    case 4:
        if (pictureIndex == 0) {
            av1Rps->refDpbIndex[0] = base0_idx;
            av1Rps->refDpbIndex[6] = layer3_idx1;
        }
        else if (pictureIndex == 2) {
            av1Rps->refDpbIndex[0] = layer3_idx1;
            av1Rps->refDpbIndex[6] = layer2_idx;
        }
        else if (pictureIndex == 4) {
            av1Rps->refDpbIndex[0] = layer2_idx;
            av1Rps->refDpbIndex[6] = layer3_idx2;
        }
        else if (pictureIndex == 6) {
            av1Rps->refDpbIndex[0] = layer3_idx2;
            av1Rps->refDpbIndex[6] = layer1_idx;
        }
        else if (pictureIndex == 8) {
            av1Rps->refDpbIndex[0] = layer1_idx;
            av1Rps->refDpbIndex[6] = layer3_idx1;
        }
        else if (pictureIndex == 10) {
            av1Rps->refDpbIndex[0] = layer3_idx1;
            av1Rps->refDpbIndex[6] = layer2_idx;
        }
        else if (pictureIndex == 12) {
            av1Rps->refDpbIndex[0] = layer2_idx;
            av1Rps->refDpbIndex[6] = layer3_idx2;
        }
        else if (pictureIndex == 14) {
            av1Rps->refDpbIndex[0] = layer3_idx2;
            av1Rps->refDpbIndex[6] = base1_idx;
        }
        else {
            printf("Error in GOp indexing\n");
        }
        av1Rps->refreshFrameMask = 0;
        break;
    default:
        printf("Error: unexpected picture mini Gop number\n");
        break;
    }

    if (picture_control_set_ptr->pred_struct_ptr->predType == EB_PRED_LOW_DELAY_P)
    {
        //P frames.
        av1Rps->refDpbIndex[1] = av1Rps->refDpbIndex[2] = av1Rps->refDpbIndex[3] = av1Rps->refDpbIndex[0];
        av1Rps->refDpbIndex[4] = av1Rps->refDpbIndex[5] = av1Rps->refDpbIndex[6] = av1Rps->refDpbIndex[0];
        picture_control_set_ptr->showFrame = EB_TRUE;
        picture_control_set_ptr->hasShowExisting = EB_FALSE;
    }
    else if (picture_control_set_ptr->pred_struct_ptr->predType == EB_PRED_RANDOM_ACCESS)
    {
        av1Rps->refDpbIndex[1] = av1Rps->refDpbIndex[2] = av1Rps->refDpbIndex[3] = av1Rps->refDpbIndex[0];
        av1Rps->refDpbIndex[4] = av1Rps->refDpbIndex[5] = av1Rps->refDpbIndex[6];

        //Decide on Show Mecanism
        if (picture_control_set_ptr->slice_type == I_SLICE)
        {
            //3 cases for I slice:  1:Key Frame treated above.  2: broken MiniGop due to sc or intra refresh  3: complete miniGop due to sc or intra refresh
            if (context_ptr->miniGopLength[0] < picture_control_set_ptr->pred_struct_ptr->predStructPeriod)
            {
                //Scene Change that breaks the mini gop and switch to LDP (if I scene change happens to be aligned with a complete miniGop, then we do not break the pred structure)
                picture_control_set_ptr->showFrame = EB_TRUE;
                picture_control_set_ptr->hasShowExisting = EB_FALSE;
            }
    else
    {
                picture_control_set_ptr->showFrame = EB_FALSE;
                picture_control_set_ptr->hasShowExisting = EB_FALSE;
            }
        }
        else//B pic
        {
            if (context_ptr->miniGopLength[0] != picture_control_set_ptr->pred_struct_ptr->predStructPeriod)
                printf("Error in GOp indexing3\n");

            if (picture_control_set_ptr->is_used_as_reference_flag)
            {
                picture_control_set_ptr->showFrame = EB_FALSE;
                picture_control_set_ptr->hasShowExisting = EB_FALSE;
            }
            else
            {
                picture_control_set_ptr->showFrame = EB_TRUE;
                picture_control_set_ptr->hasShowExisting = EB_TRUE;

                if (pictureIndex == 0) {
                    picture_control_set_ptr->showExistingLoc = layer3_idx1;
                }
                else if (pictureIndex == 2) {
                    picture_control_set_ptr->showExistingLoc = layer2_idx;
                }
                else if (pictureIndex == 4) {
                    picture_control_set_ptr->showExistingLoc = layer3_idx2;
                }
                else if (pictureIndex == 6) {
                    picture_control_set_ptr->showExistingLoc = layer1_idx;
                }
                else if (pictureIndex == 8) {
                    picture_control_set_ptr->showExistingLoc = layer3_idx1;
                }
                else if (pictureIndex == 10) {
                    picture_control_set_ptr->showExistingLoc = layer2_idx;
                }
                else if (pictureIndex == 12) {
                    picture_control_set_ptr->showExistingLoc = layer3_idx2;
                }
                else if (pictureIndex == 14) {
                    picture_control_set_ptr->showExistingLoc = base1_idx;
                }
                else {
                    printf("Error in GOp indexing2\n");
                }

            }

        }

    }
    else {
        printf("Error: Not supported GOP structure!");
        exit(0);
    }

    //last pic in MiniGop: mGop Toggling
    //mini GOP toggling since last Key Frame.
    //a regular I keeps the toggling process and does not reset the toggle.  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
    if (pictureIndex == context_ptr->miniGopEndIndex[0])
        context_ptr->miniGopToggle = 1 - context_ptr->miniGopToggle;

    }
#endif
    else
    {
        printf("Error: Not supported GOP structure!");
        exit(0);
    }
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
void* PictureDecisionKernel(void *input_ptr)
{
    PictureDecisionContext_t        *context_ptr = (PictureDecisionContext_t*)input_ptr;

    PictureParentControlSet_t       *picture_control_set_ptr;

    EncodeContext_t                 *encode_context_ptr;
    SequenceControlSet_t            *sequence_control_set_ptr;

    EbObjectWrapper_t               *inputResultsWrapperPtr;
    PictureAnalysisResults_t        *inputResultsPtr;

    EbObjectWrapper_t               *outputResultsWrapperPtr;
    PictureDecisionResults_t        *outputResultsPtr;

    PredictionStructureEntry_t      *predPositionPtr;

    EbBool                          preAssignmentBufferFirstPassFlag;
    EB_SLICE                         pictureType;

    PictureDecisionReorderEntry_t   *queueEntryPtr;
    int32_t                           queueEntryIndex;

    int32_t                           previousEntryIndex;

    PaReferenceQueueEntry_t         *inputEntryPtr;
    uint32_t                           inputQueueIndex;

    PaReferenceQueueEntry_t         *paReferenceEntryPtr;
    uint32_t                           paReferenceQueueIndex;

    uint64_t                           refPoc;

    uint32_t                           depIdx;
    uint64_t                           depPoc;

    uint32_t                           depListCount;

    // Dynamic GOP
    uint32_t                           miniGopIndex;
    uint32_t                           pictureIndex;

    // Initialization
    uint32_t                           picture_width_in_sb;

    EbBool                          windowAvail, framePasseThru;
    uint32_t                           windowIndex;
    uint32_t                           entryIndex;
    PictureParentControlSet_t        *ParentPcsWindow[FUTURE_WINDOW_WIDTH + 2];

    // Debug
    uint64_t                           loopCount = 0;

    for (;;) {

        // Get Input Full Object
        EbGetFullObject(
            context_ptr->pictureAnalysisResultsInputFifoPtr,
            &inputResultsWrapperPtr);

        inputResultsPtr = (PictureAnalysisResults_t*)inputResultsWrapperPtr->objectPtr;
        picture_control_set_ptr = (PictureParentControlSet_t*)inputResultsPtr->pictureControlSetWrapperPtr->objectPtr;
        sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
        encode_context_ptr = (EncodeContext_t*)sequence_control_set_ptr->encode_context_ptr;

        loopCount++;

        // Input Picture Analysis Results into the Picture Decision Reordering Queue
        // P.S. Since the prior Picture Analysis processes stage is multithreaded, inputs to the Picture Decision Process
        // can arrive out-of-display-order, so a the Picture Decision Reordering Queue is used to enforce processing of
        // pictures in display order

        queueEntryIndex = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index]->picture_number);
        queueEntryIndex += encode_context_ptr->picture_decision_reorder_queue_head_index;
        queueEntryIndex = (queueEntryIndex > PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH - 1) ? queueEntryIndex - PICTURE_DECISION_REORDER_QUEUE_MAX_DEPTH : queueEntryIndex;
        queueEntryPtr = encode_context_ptr->picture_decision_reorder_queue[queueEntryIndex];
        if (queueEntryPtr->parentPcsWrapperPtr != NULL) {
            CHECK_REPORT_ERROR_NC(
                encode_context_ptr->app_callback_ptr,
                EB_ENC_PD_ERROR8);
        }
        else {
            queueEntryPtr->parentPcsWrapperPtr = inputResultsPtr->pictureControlSetWrapperPtr;
            queueEntryPtr->picture_number = picture_control_set_ptr->picture_number;
        }
        // Process the head of the Picture Decision Reordering Queue (Entry N)
        // P.S. The Picture Decision Reordering Queue should be parsed in the display order to be able to construct a pred structure
        queueEntryPtr = encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index];

        while (queueEntryPtr->parentPcsWrapperPtr != EB_NULL) {

            if (queueEntryPtr->picture_number == 0 ||
                ((PictureParentControlSet_t *)(queueEntryPtr->parentPcsWrapperPtr->objectPtr))->end_of_sequence_flag == EB_TRUE){
                framePasseThru = EB_TRUE;
            }
            else {
                framePasseThru = EB_FALSE;
            }
            windowAvail = EB_TRUE;
            previousEntryIndex = QUEUE_GET_PREVIOUS_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index);

            if (encode_context_ptr->picture_decision_reorder_queue[previousEntryIndex]->parentPcsWrapperPtr == NULL) {
                windowAvail = EB_FALSE;
            }
            else {
                ParentPcsWindow[0] = (PictureParentControlSet_t *)encode_context_ptr->picture_decision_reorder_queue[previousEntryIndex]->parentPcsWrapperPtr->objectPtr;
                ParentPcsWindow[1] = (PictureParentControlSet_t *)encode_context_ptr->picture_decision_reorder_queue[encode_context_ptr->picture_decision_reorder_queue_head_index]->parentPcsWrapperPtr->objectPtr;
                for (windowIndex = 0; windowIndex < FUTURE_WINDOW_WIDTH; windowIndex++) {
                    entryIndex = QUEUE_GET_NEXT_SPOT(encode_context_ptr->picture_decision_reorder_queue_head_index, windowIndex + 1);
                    if (encode_context_ptr->picture_decision_reorder_queue[entryIndex]->parentPcsWrapperPtr == NULL) {
                        windowAvail = EB_FALSE;
                        break;
                    }
                    else if (((PictureParentControlSet_t *)(encode_context_ptr->picture_decision_reorder_queue[entryIndex]->parentPcsWrapperPtr->objectPtr))->end_of_sequence_flag == EB_TRUE) {
                        windowAvail = EB_FALSE;
                        framePasseThru = EB_TRUE;
                        break;
                    }else {
                        ParentPcsWindow[2 + windowIndex] = (PictureParentControlSet_t *)encode_context_ptr->picture_decision_reorder_queue[entryIndex]->parentPcsWrapperPtr->objectPtr;
                    }
                }
            }
            picture_control_set_ptr = (PictureParentControlSet_t*)queueEntryPtr->parentPcsWrapperPtr->objectPtr;

            picture_control_set_ptr->fade_out_from_black = 0;

            picture_control_set_ptr->fade_in_to_black = 0;
            if (picture_control_set_ptr->idr_flag == EB_TRUE)
                context_ptr->lastSolidColorFramePoc = 0xFFFFFFFF;

            if (windowAvail == EB_TRUE) {
                if (sequence_control_set_ptr->static_config.scene_change_detection) {

                    picture_control_set_ptr->scene_change_flag = SceneTransitionDetector(
                        context_ptr,
                        sequence_control_set_ptr,
                        ParentPcsWindow,
                        FUTURE_WINDOW_WIDTH);


                }
                else {
                    picture_control_set_ptr->scene_change_flag = EB_FALSE;
                }
                picture_control_set_ptr->cra_flag = (picture_control_set_ptr->scene_change_flag == EB_TRUE) ?
                    EB_TRUE :
                    picture_control_set_ptr->cra_flag;

                // Store scene change in context
                context_ptr->isSceneChangeDetected = picture_control_set_ptr->scene_change_flag;

            }

            if (windowAvail == EB_TRUE || framePasseThru == EB_TRUE)
            {
                // Place the PCS into the Pre-Assignment Buffer
                // P.S. The Pre-Assignment Buffer is used to store a whole pre-structure
                encode_context_ptr->pre_assignment_buffer[encode_context_ptr->pre_assignment_buffer_count] = queueEntryPtr->parentPcsWrapperPtr;

                // Setup the PCS & SCS
                picture_control_set_ptr = (PictureParentControlSet_t*)encode_context_ptr->pre_assignment_buffer[encode_context_ptr->pre_assignment_buffer_count]->objectPtr;

                // Set the POC Number
                picture_control_set_ptr->picture_number = (encode_context_ptr->current_input_poc + 1) /*& ((1 << sequence_control_set_ptr->bits_for_picture_order_count)-1)*/;
                encode_context_ptr->current_input_poc = picture_control_set_ptr->picture_number;


                picture_control_set_ptr->pred_structure = EB_PRED_RANDOM_ACCESS;

#if NEW_PRED_STRUCT
                picture_control_set_ptr->hierarchical_layers_diff = 0;

                picture_control_set_ptr->init_pred_struct_position_flag = EB_FALSE;
#endif
                picture_control_set_ptr->target_bit_rate = sequence_control_set_ptr->static_config.target_bit_rate;

                ReleasePrevPictureFromReorderQueue(
                    encode_context_ptr);

                // If the Intra period length is 0, then introduce an intra for every picture
                if (sequence_control_set_ptr->intra_period_length == 0) {
                    picture_control_set_ptr->cra_flag = EB_TRUE;
                }
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
                    context_ptr->miniGopStartIndex[0] = 0;
                    context_ptr->miniGopEndIndex[0] = encode_context_ptr->pre_assignment_buffer_count - 1;
                    context_ptr->miniGopLength[0] = encode_context_ptr->pre_assignment_buffer_count;

                    context_ptr->miniGopHierarchicalLevels[0] = sequence_control_set_ptr->static_config.hierarchical_levels; 
                    context_ptr->miniGopIntraCount[0] = encode_context_ptr->pre_assignment_buffer_intra_count;
                    context_ptr->miniGopIdrCount[0] = encode_context_ptr->pre_assignment_buffer_idr_count;
                    context_ptr->totalNumberOfMiniGops = 1;

                    encode_context_ptr->previous_mini_gop_hierarchical_levels = (picture_control_set_ptr->picture_number == 0) ?
                        sequence_control_set_ptr->static_config.hierarchical_levels : 
                        encode_context_ptr->previous_mini_gop_hierarchical_levels;

#if NEW_PRED_STRUCT
                    {
                        if (encode_context_ptr->pre_assignment_buffer_count > 1)
                        {
                            initialize_mini_gop_activity_array(
                                context_ptr);

                            if (encode_context_ptr->pre_assignment_buffer_count == 16) 
                                context_ptr->miniGopActivityArray[L5_0_INDEX] = EB_FALSE;
                            else {
                                context_ptr->miniGopActivityArray[L4_0_INDEX] = EB_FALSE;
                                context_ptr->miniGopActivityArray[L4_1_INDEX] = EB_FALSE;
                            }

                            generate_picture_window_split(
                                context_ptr,
                                encode_context_ptr);

                            handle_incomplete_picture_window_map(
                                context_ptr,
                                encode_context_ptr);
                        }
                    }
#endif
                    GenerateMiniGopRps(
                        context_ptr,
                        encode_context_ptr);

                    // Loop over Mini GOPs

                    for (miniGopIndex = 0; miniGopIndex < context_ptr->totalNumberOfMiniGops; ++miniGopIndex) {

                        preAssignmentBufferFirstPassFlag = EB_TRUE;
#if NEW_PRED_STRUCT
                        {
                            update_base_layer_reference_queue_dependent_count(
                                context_ptr,
                                encode_context_ptr,
                                sequence_control_set_ptr,
                                miniGopIndex);

                            // Keep track of the number of hierarchical levels of the latest implemented mini GOP
                            encode_context_ptr->previous_mini_gop_hierarchical_levels = context_ptr->miniGopHierarchicalLevels[miniGopIndex];
                        }
#endif
                        // 1st Loop over Pictures in the Pre-Assignment Buffer
                        for (pictureIndex = context_ptr->miniGopStartIndex[miniGopIndex]; pictureIndex <= context_ptr->miniGopEndIndex[miniGopIndex]; ++pictureIndex) {

                            picture_control_set_ptr = (PictureParentControlSet_t*)encode_context_ptr->pre_assignment_buffer[pictureIndex]->objectPtr;
                            sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;

                            // Keep track of the mini GOP size to which the input picture belongs - needed @ PictureManagerProcess()
                            picture_control_set_ptr->pre_assignment_buffer_count = context_ptr->miniGopLength[miniGopIndex];

                            // Update the Pred Structure if cutting short a Random Access period
                            if ((context_ptr->miniGopLength[miniGopIndex] < picture_control_set_ptr->pred_struct_ptr->predStructPeriod || context_ptr->miniGopIdrCount[miniGopIndex] > 0) &&

                                picture_control_set_ptr->pred_struct_ptr->predType == EB_PRED_RANDOM_ACCESS &&
                                picture_control_set_ptr->idr_flag == EB_FALSE &&
                                picture_control_set_ptr->cra_flag == EB_FALSE)
                            {
                                // Correct the Pred Index before switching structures
                                if (preAssignmentBufferFirstPassFlag == EB_TRUE) {
                                    encode_context_ptr->pred_struct_position -= picture_control_set_ptr->pred_struct_ptr->initPicIndex;
                                }

                                picture_control_set_ptr->pred_struct_ptr = GetPredictionStructure(
                                    encode_context_ptr->prediction_structure_group_ptr,
                                    EB_PRED_LOW_DELAY_P,
                                    1,
                                    picture_control_set_ptr->hierarchical_levels);

                                // Set the RPS Override Flag - this current only will convert a Random Access structure to a Low Delay structure
                                picture_control_set_ptr->use_rps_in_sps = EB_FALSE;
                                picture_control_set_ptr->open_gop_cra_flag = EB_FALSE;

                                pictureType = P_SLICE;

                            }
                            // Open GOP CRA - adjust the RPS
                            else if ((context_ptr->miniGopLength[miniGopIndex] == picture_control_set_ptr->pred_struct_ptr->predStructPeriod) &&

                                (picture_control_set_ptr->pred_struct_ptr->predType == EB_PRED_RANDOM_ACCESS || picture_control_set_ptr->pred_struct_ptr->temporalLayerCount == 1) &&
                                picture_control_set_ptr->idr_flag == EB_FALSE &&
                                picture_control_set_ptr->cra_flag == EB_TRUE)
                            {
                                picture_control_set_ptr->use_rps_in_sps = EB_FALSE;
                                picture_control_set_ptr->open_gop_cra_flag = EB_TRUE;

                                pictureType = I_SLICE;
                            }
                            else {

                                picture_control_set_ptr->use_rps_in_sps = EB_FALSE;
                                picture_control_set_ptr->open_gop_cra_flag = EB_FALSE;

                                // Set the Picture Type
                                pictureType =
                                    (picture_control_set_ptr->idr_flag) ? I_SLICE :
                                    (picture_control_set_ptr->cra_flag) ? I_SLICE :
                                    (picture_control_set_ptr->pred_structure == EB_PRED_LOW_DELAY_P) ? P_SLICE :
                                    (picture_control_set_ptr->pred_structure == EB_PRED_LOW_DELAY_B) ? B_SLICE :
                                    (picture_control_set_ptr->pre_assignment_buffer_count == picture_control_set_ptr->pred_struct_ptr->predStructPeriod) ? ((pictureIndex == context_ptr->miniGopEndIndex[miniGopIndex] && sequence_control_set_ptr->static_config.base_layer_switch_mode) ? P_SLICE : B_SLICE) :

                                    (encode_context_ptr->pre_assignment_buffer_eos_flag) ? P_SLICE :
                                    B_SLICE;
                            }
#if NEW_PRED_STRUCT
                            // If mini GOP switch, reset position
                            encode_context_ptr->pred_struct_position = (picture_control_set_ptr->init_pred_struct_position_flag) ?
                                picture_control_set_ptr->pred_struct_ptr->initPicIndex :
                                encode_context_ptr->pred_struct_position;
#endif

                            // If Intra, reset position
                            if (picture_control_set_ptr->idr_flag == EB_TRUE) {
                                encode_context_ptr->pred_struct_position = picture_control_set_ptr->pred_struct_ptr->initPicIndex;
                            }

                            else if (picture_control_set_ptr->cra_flag == EB_TRUE && context_ptr->miniGopLength[miniGopIndex] < picture_control_set_ptr->pred_struct_ptr->predStructPeriod) {

                                encode_context_ptr->pred_struct_position = picture_control_set_ptr->pred_struct_ptr->initPicIndex;
                            }
                            else if (encode_context_ptr->elapsed_non_cra_count == 0) {
                                // If we are the picture directly after a CRA, we have to not use references that violate the CRA
                                encode_context_ptr->pred_struct_position = picture_control_set_ptr->pred_struct_ptr->initPicIndex + 1;
                            }
                            // Elif Scene Change, determine leading and trailing pictures
                            //else if (encode_context_ptr->pre_assignment_buffer_scene_change_count > 0) {
                            //    if(bufferIndex < encode_context_ptr->pre_assignment_buffer_scene_change_index) {
                            //        ++encode_context_ptr->pred_struct_position;
                            //        pictureType = P_SLICE;
                            //    }
                            //    else {
                            //        encode_context_ptr->pred_struct_position = picture_control_set_ptr->pred_struct_ptr->initPicIndex + encode_context_ptr->pre_assignment_buffer_count - bufferIndex - 1;
                            //    }
                            //}
                            // Else, Increment the position normally
                            else {
                                ++encode_context_ptr->pred_struct_position;
                            }

                            // The poc number of the latest IDR picture is stored so that last_idr_picture (present in PCS) for the incoming pictures can be updated.
                            // The last_idr_picture is used in reseting the poc (in entropy coding) whenever IDR is encountered.
                            // Note IMP: This logic only works when display and decode order are the same. Currently for Random Access, IDR is inserted (similar to CRA) by using trailing P pictures (low delay fashion) and breaking prediction structure.
                            // Note: When leading P pictures are implemented, this logic has to change..
                            if (picture_control_set_ptr->idr_flag == EB_TRUE) {
                                encode_context_ptr->last_idr_picture = picture_control_set_ptr->picture_number;
                            }
                            else {
                                picture_control_set_ptr->last_idr_picture = encode_context_ptr->last_idr_picture;
                            }


                            // Cycle the PredStructPosition if its overflowed
                            encode_context_ptr->pred_struct_position = (encode_context_ptr->pred_struct_position == picture_control_set_ptr->pred_struct_ptr->predStructEntryCount) ?
                                encode_context_ptr->pred_struct_position - picture_control_set_ptr->pred_struct_ptr->predStructPeriod :
                                encode_context_ptr->pred_struct_position;

                            predPositionPtr = picture_control_set_ptr->pred_struct_ptr->predStructEntryPtrArray[encode_context_ptr->pred_struct_position];

                            // Set the Slice type
                            picture_control_set_ptr->slice_type = pictureType;
                            ((EbPaReferenceObject_t*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->objectPtr)->slice_type = picture_control_set_ptr->slice_type;

                            switch (pictureType) {

                            case I_SLICE:

                                // Reset Prediction Structure Position & Reference Struct Position
                                if (picture_control_set_ptr->picture_number == 0) {
                                    encode_context_ptr->intra_period_position = 0;
                                }
                                encode_context_ptr->elapsed_non_cra_count = 0;

                                //-------------------------------
                                // IDR
                                //-------------------------------
                                if (picture_control_set_ptr->idr_flag == EB_TRUE) {

                                    // Set CRA flag
                                    picture_control_set_ptr->cra_flag = EB_FALSE;

                                    // Reset the pictures since last IDR counter
                                    encode_context_ptr->elapsed_non_idr_count = 0;

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

                                // Increment & Clip the elapsed Non-IDR Counter. This is clipped rather than allowed to free-run
                                // inorder to avoid rollover issues.  This assumes that any the GOP period is less than MAX_ELAPSED_IDR_COUNT
                                encode_context_ptr->elapsed_non_idr_count = MIN(encode_context_ptr->elapsed_non_idr_count + 1, MAX_ELAPSED_IDR_COUNT);
                                encode_context_ptr->elapsed_non_cra_count = MIN(encode_context_ptr->elapsed_non_cra_count + 1, MAX_ELAPSED_IDR_COUNT);

                                CHECK_REPORT_ERROR(
                                    (picture_control_set_ptr->pred_struct_ptr->predStructEntryCount < MAX_ELAPSED_IDR_COUNT),
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
                            picture_control_set_ptr->temporal_layer_index = (uint8_t)predPositionPtr->temporal_layer_index;
                            picture_control_set_ptr->is_used_as_reference_flag = predPositionPtr->isReferenced;

                            // Set the Decode Order
                            if ((context_ptr->miniGopIdrCount[miniGopIndex] == 0) &&
                                (context_ptr->miniGopLength[miniGopIndex] == picture_control_set_ptr->pred_struct_ptr->predStructPeriod))

                            {
                                picture_control_set_ptr->decode_order = encode_context_ptr->decode_base_number + predPositionPtr->decode_order;
                            }
                            else {
                                picture_control_set_ptr->decode_order = picture_control_set_ptr->picture_number;
                            }

                            encode_context_ptr->terminating_sequence_flag_received = (picture_control_set_ptr->end_of_sequence_flag == EB_TRUE) ?
                                EB_TRUE :
                                encode_context_ptr->terminating_sequence_flag_received;

                            encode_context_ptr->terminating_picture_number = (picture_control_set_ptr->end_of_sequence_flag == EB_TRUE) ?
                                picture_control_set_ptr->picture_number :
                                encode_context_ptr->terminating_picture_number;

                            preAssignmentBufferFirstPassFlag = EB_FALSE;

                            // Film grain (assigning the random-seed)
                            {
                                uint16_t *fgn_random_seed_ptr = &picture_control_set_ptr->sequence_control_set_ptr->film_grain_random_seed;
                                picture_control_set_ptr->film_grain_params.random_seed = *fgn_random_seed_ptr;
                                *fgn_random_seed_ptr += 3381;  // Changing random seed for film grain
                                if (!(*fgn_random_seed_ptr))     // Random seed should not be zero
                                    *fgn_random_seed_ptr += 7391;
                            }

                            Av1GenerateRpsInfo(
                                picture_control_set_ptr,
                                encode_context_ptr,
                                context_ptr,
#if NEW_PRED_STRUCT
                                pictureIndex - context_ptr->miniGopStartIndex[miniGopIndex]);
#else
                                pictureIndex);
#endif
                            picture_control_set_ptr->allow_comp_inter_inter = 0;
                            picture_control_set_ptr->is_skip_mode_allowed = 0;

                            picture_control_set_ptr->reference_mode = (ReferenceMode)0xFF;

                            if (picture_control_set_ptr->slice_type != I_SLICE) {

                                if (picture_control_set_ptr->temporal_layer_index == 0 || picture_control_set_ptr->slice_type == P_SLICE) {

                                    picture_control_set_ptr->allow_comp_inter_inter = 1;

                                    picture_control_set_ptr->reference_mode = SINGLE_REFERENCE;
                                    picture_control_set_ptr->is_skip_mode_allowed = 0;
                                    picture_control_set_ptr->skip_mode_flag = 0;
                                }
                                else {
                                    picture_control_set_ptr->allow_comp_inter_inter = 1;
                                    picture_control_set_ptr->reference_mode = REFERENCE_MODE_SELECT;
                                    picture_control_set_ptr->is_skip_mode_allowed = 1;
                                    picture_control_set_ptr->skip_mode_flag = 1;
                                }
                            }

                            picture_control_set_ptr->av1_cm->mi_cols = picture_control_set_ptr->sequence_control_set_ptr->luma_width >> MI_SIZE_LOG2;
                            picture_control_set_ptr->av1_cm->mi_rows = picture_control_set_ptr->sequence_control_set_ptr->luma_height >> MI_SIZE_LOG2;

                            memset(picture_control_set_ptr->av1_cm->ref_frame_sign_bias, 0, 8 * sizeof(int32_t));
                            if (picture_control_set_ptr->reference_mode == REFERENCE_MODE_SELECT)
                            {
                                picture_control_set_ptr->av1_cm->ref_frame_sign_bias[ALTREF_FRAME] =
                                    picture_control_set_ptr->av1_cm->ref_frame_sign_bias[ALTREF2_FRAME] =
                                    picture_control_set_ptr->av1_cm->ref_frame_sign_bias[BWDREF_FRAME] = 1;
                            }

                            // ME Kernel Multi-Processes Signal(s) derivation
                            signal_derivation_multi_processes_oq(
                                picture_control_set_ptr);

                            // Set the default settings of  subpel
                            picture_control_set_ptr->use_subpel_flag = PictureLevelSubPelSettings(
                                sequence_control_set_ptr->input_resolution,
                                picture_control_set_ptr->enc_mode,
                                picture_control_set_ptr->temporal_layer_index,
                                picture_control_set_ptr->is_used_as_reference_flag);
#if !CHROMA_BLIND
                            // Set the default settings of  chroma
                            picture_control_set_ptr->chroma_mode = PictureLevelChromaSettings(
                                sequence_control_set_ptr->input_resolution,
                                picture_control_set_ptr->enc_mode,
                                picture_control_set_ptr->slice_type,
                                picture_control_set_ptr->temporal_layer_index,
                                picture_control_set_ptr->is_used_as_reference_flag);
#endif

                            picture_control_set_ptr->use_src_ref = EB_FALSE;
#if DISABLE_IN_LOOP_ME
                            picture_control_set_ptr->enable_in_loop_motion_estimation_flag = EB_FALSE;
#else
                            picture_control_set_ptr->enable_in_loop_motion_estimation_flag = sequence_control_set_ptr->static_config.in_loop_me_flag && picture_control_set_ptr->slice_type != I_SLICE ? EB_TRUE : EB_FALSE;
#endif
#if ENCODER_MODE_CLEANUP
                            picture_control_set_ptr->limit_ois_to_dc_mode_flag = EB_FALSE;
#else
                            picture_control_set_ptr->limit_ois_to_dc_mode_flag = picture_control_set_ptr->enc_mode >= ENC_M6 &&
                                picture_control_set_ptr->slice_type != I_SLICE ? EB_TRUE : EB_FALSE;
#endif
#if ENCODER_MODE_CLEANUP
                            picture_control_set_ptr->cu8x8_mode = CU_8x8_MODE_0;
#else
                            if ((picture_control_set_ptr->enc_mode > ENC_M1 && sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE) || (picture_control_set_ptr->enc_mode > ENC_M3 && sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE)) {

                                if (picture_control_set_ptr->enc_mode == ENC_M2 && sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE)
                                    picture_control_set_ptr->cu8x8_mode = (picture_control_set_ptr->is_used_as_reference_flag) ? CU_8x8_MODE_0 : CU_8x8_MODE_1;
                                else
                                    picture_control_set_ptr->cu8x8_mode = (picture_control_set_ptr->temporal_layer_index == 0) ? CU_8x8_MODE_0 : CU_8x8_MODE_1;

                            }
                            else {
                                picture_control_set_ptr->cu8x8_mode = CU_8x8_MODE_0;
                            }
#endif


                            // Update the Dependant List Count - If there was an I-frame or Scene Change, then cleanup the Picture Decision PA Reference Queue Dependent Counts
                            if (picture_control_set_ptr->slice_type == I_SLICE)
                            {

                                inputQueueIndex = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                                while (inputQueueIndex != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) {

                                    inputEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[inputQueueIndex];

                                    // Modify Dependent List0
                                    depListCount = inputEntryPtr->list0.listCount;
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
                                            --inputEntryPtr->dependentCount;

                                            CHECK_REPORT_ERROR(
                                                (inputEntryPtr->dependentCount != ~0u),
                                                encode_context_ptr->app_callback_ptr,
                                                EB_ENC_PD_ERROR3);
                                        }
                                    }

                                    // Modify Dependent List1
                                    depListCount = inputEntryPtr->list1.listCount;
                                    for (depIdx = 0; depIdx < depListCount; ++depIdx) {

                                        // Adjust the latest current_input_poc in case we're in a POC rollover scenario
                                        // current_input_poc += (current_input_poc < inputEntryPtr->pocNumber) ? (1 << sequence_control_set_ptr->bits_for_picture_order_count) : 0;

                                        depPoc = POC_CIRCULAR_ADD(
                                            inputEntryPtr->picture_number,
                                            inputEntryPtr->list1.list[depIdx]/*,
                                            sequence_control_set_ptr->bits_for_picture_order_count*/);

                                            // If Dependent POC is greater or equal to the IDR POC
                                        if (((depPoc >= picture_control_set_ptr->picture_number) || (((picture_control_set_ptr->pre_assignment_buffer_count != picture_control_set_ptr->pred_struct_ptr->predStructPeriod) || (picture_control_set_ptr->idr_flag == EB_TRUE)) && (depPoc > (picture_control_set_ptr->picture_number - picture_control_set_ptr->pre_assignment_buffer_count)))) && inputEntryPtr->list1.list[depIdx]) {
                                            inputEntryPtr->list1.list[depIdx] = 0;

                                            // Decrement the Reference's referenceCount
                                            --inputEntryPtr->dependentCount;

                                            CHECK_REPORT_ERROR(
                                                (inputEntryPtr->dependentCount != ~0u),
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
                                inputEntryPtr = (PaReferenceQueueEntry_t*)EB_NULL;
                            }

                            // Place Picture in Picture Decision PA Reference Queue
                            inputEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_tail_index];
                            inputEntryPtr->inputObjectPtr = picture_control_set_ptr->pa_reference_picture_wrapper_ptr;
                            inputEntryPtr->picture_number = picture_control_set_ptr->picture_number;
                            inputEntryPtr->referenceEntryIndex = encode_context_ptr->picture_decision_pa_reference_queue_tail_index;
                            inputEntryPtr->pPcsPtr = picture_control_set_ptr;
                            encode_context_ptr->picture_decision_pa_reference_queue_tail_index =
                                (encode_context_ptr->picture_decision_pa_reference_queue_tail_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->picture_decision_pa_reference_queue_tail_index + 1;

                            // Check if the Picture Decision PA Reference is full
                            CHECK_REPORT_ERROR(
                                (((encode_context_ptr->picture_decision_pa_reference_queue_head_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) ||
                                (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->inputObjectPtr == EB_NULL))),
                                encode_context_ptr->app_callback_ptr,
                                EB_ENC_PD_ERROR4);

                            // Copy the reference lists into the inputEntry and
                            // set the Reference Counts Based on Temporal Layer and how many frames are active
                            picture_control_set_ptr->ref_list0_count = (pictureType == I_SLICE) ? 0 : (uint8_t)predPositionPtr->refList0.referenceListCount;
                            picture_control_set_ptr->ref_list1_count = (pictureType == I_SLICE) ? 0 : (uint8_t)predPositionPtr->refList1.referenceListCount;

                            inputEntryPtr->list0Ptr = &predPositionPtr->refList0;
                            inputEntryPtr->list1Ptr = &predPositionPtr->refList1;

                            {

                                // Copy the Dependent Lists
                                // *Note - we are removing any leading picture dependencies for now
                                inputEntryPtr->list0.listCount = 0;
                                for (depIdx = 0; depIdx < predPositionPtr->depList0.listCount; ++depIdx) {
                                    if (predPositionPtr->depList0.list[depIdx] >= 0) {
                                        inputEntryPtr->list0.list[inputEntryPtr->list0.listCount++] = predPositionPtr->depList0.list[depIdx];
                                    }
                                }

                                inputEntryPtr->list1.listCount = predPositionPtr->depList1.listCount;
                                for (depIdx = 0; depIdx < predPositionPtr->depList1.listCount; ++depIdx) {
                                    inputEntryPtr->list1.list[depIdx] = predPositionPtr->depList1.list[depIdx];
                                }

                                inputEntryPtr->depList0Count = inputEntryPtr->list0.listCount;
                                inputEntryPtr->depList1Count = inputEntryPtr->list1.listCount;
                                inputEntryPtr->dependentCount = inputEntryPtr->depList0Count + inputEntryPtr->depList1Count;

                            }

                            ((EbPaReferenceObject_t*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->objectPtr)->dependentPicturesCount = inputEntryPtr->dependentCount;

                            /* uint32_t depCnt = ((EbPaReferenceObject_t*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->objectPtr)->dependentPicturesCount;
                            if (picture_control_set_ptr->picture_number>0 && picture_control_set_ptr->slice_type==I_SLICE && depCnt!=8 )
                            printf("depCnt Error1  POC:%i  TL:%i   is needed:%i\n",picture_control_set_ptr->picture_number,picture_control_set_ptr->temporal_layer_index,inputEntryPtr->dependentCount);
                            else if (picture_control_set_ptr->slice_type==B_SLICE && picture_control_set_ptr->temporal_layer_index == 0 && depCnt!=8)
                            printf("depCnt Error2  POC:%i  TL:%i   is needed:%i\n",picture_control_set_ptr->picture_number,picture_control_set_ptr->temporal_layer_index,inputEntryPtr->dependentCount);
                            else if (picture_control_set_ptr->slice_type==B_SLICE && picture_control_set_ptr->temporal_layer_index == 1 && depCnt!=4)
                            printf("depCnt Error3  POC:%i  TL:%i   is needed:%i\n",picture_control_set_ptr->picture_number,picture_control_set_ptr->temporal_layer_index,inputEntryPtr->dependentCount);
                            else if (picture_control_set_ptr->slice_type==B_SLICE && picture_control_set_ptr->temporal_layer_index == 2 && depCnt!=2)
                            printf("depCnt Error4  POC:%i  TL:%i   is needed:%i\n",picture_control_set_ptr->picture_number,picture_control_set_ptr->temporal_layer_index,inputEntryPtr->dependentCount);
                            else if (picture_control_set_ptr->slice_type==B_SLICE && picture_control_set_ptr->temporal_layer_index == 3 && depCnt!=0)
                            printf("depCnt Error5  POC:%i  TL:%i   is needed:%i\n",picture_control_set_ptr->picture_number,picture_control_set_ptr->temporal_layer_index,inputEntryPtr->dependentCount);*/
                            //if (picture_control_set_ptr->slice_type==P_SLICE )
                            //     printf("POC:%i  TL:%i   is needed:%i\n",picture_control_set_ptr->picture_number,picture_control_set_ptr->temporal_layer_index,inputEntryPtr->dependentCount);

                            CHECK_REPORT_ERROR(
                                (picture_control_set_ptr->pred_struct_ptr->predStructPeriod < MAX_ELAPSED_IDR_COUNT),
                                encode_context_ptr->app_callback_ptr,
                                EB_ENC_PD_ERROR5);

                            // Reset the PA Reference Lists
                            EB_MEMSET(picture_control_set_ptr->ref_pa_pic_ptr_array, 0, 2 * sizeof(EbObjectWrapper_t*));

                            EB_MEMSET(picture_control_set_ptr->ref_pa_pic_ptr_array, 0, 2 * sizeof(uint32_t));

                        }

                        // 2nd Loop over Pictures in the Pre-Assignment Buffer
                        for (pictureIndex = context_ptr->miniGopStartIndex[miniGopIndex]; pictureIndex <= context_ptr->miniGopEndIndex[miniGopIndex]; ++pictureIndex) {

                            picture_control_set_ptr = (PictureParentControlSet_t*)encode_context_ptr->pre_assignment_buffer[pictureIndex]->objectPtr;

                            // Find the Reference in the Picture Decision PA Reference Queue
                            inputQueueIndex = encode_context_ptr->picture_decision_pa_reference_queue_head_index;

                            do {

                                // Setup the Picture Decision PA Reference Queue Entry
                                inputEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[inputQueueIndex];

                                // Increment the referenceQueueIndex Iterator
                                inputQueueIndex = (inputQueueIndex == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;

                            } while ((inputQueueIndex != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) && (inputEntryPtr->picture_number != picture_control_set_ptr->picture_number));

                            CHECK_REPORT_ERROR(
                                (inputEntryPtr->picture_number == picture_control_set_ptr->picture_number),
                                encode_context_ptr->app_callback_ptr,
                                EB_ENC_PD_ERROR7);

                            // Reset the PA Reference Lists
                            EB_MEMSET(picture_control_set_ptr->ref_pa_pic_ptr_array, 0, 2 * sizeof(EbObjectWrapper_t*));

                            EB_MEMSET(picture_control_set_ptr->ref_pic_poc_array, 0, 2 * sizeof(uint64_t));


                            // Configure List0
                            if ((picture_control_set_ptr->slice_type == P_SLICE) || (picture_control_set_ptr->slice_type == B_SLICE)) {

                                if (picture_control_set_ptr->ref_list0_count) {
                                    paReferenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                        ((int32_t)inputEntryPtr->referenceEntryIndex) - inputEntryPtr->list0Ptr->referenceList,
                                        PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                    paReferenceEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[paReferenceQueueIndex];

                                    // Calculate the Ref POC
                                    refPoc = POC_CIRCULAR_ADD(
                                        picture_control_set_ptr->picture_number,
                                        -inputEntryPtr->list0Ptr->referenceList/*,
                                        sequence_control_set_ptr->bits_for_picture_order_count*/);

                                        // Set the Reference Object
                                    picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_0] = paReferenceEntryPtr->inputObjectPtr;
                                    picture_control_set_ptr->ref_pic_poc_array[REF_LIST_0] = refPoc;
                                    picture_control_set_ptr->ref_pa_pcs_array[REF_LIST_0] = paReferenceEntryPtr->pPcsPtr;

                                    // Increment the PA Reference's liveCount by the number of tiles in the input picture
                                    EbObjectIncLiveCount(
                                        paReferenceEntryPtr->inputObjectPtr,
                                        1);

                                    ((EbPaReferenceObject_t*)picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_0]->objectPtr)->pPcsPtr = paReferenceEntryPtr->pPcsPtr;

                                    EbObjectIncLiveCount(
                                        paReferenceEntryPtr->pPcsPtr->p_pcs_wrapper_ptr,
                                        1);

                                    --paReferenceEntryPtr->dependentCount;
                                }
                            }

                            // Configure List1
                            if (picture_control_set_ptr->slice_type == B_SLICE) {

                                if (picture_control_set_ptr->ref_list1_count) {
                                    paReferenceQueueIndex = (uint32_t)CIRCULAR_ADD(
                                        ((int32_t)inputEntryPtr->referenceEntryIndex) - inputEntryPtr->list1Ptr->referenceList,
                                        PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH);                                                                                             // Max

                                    paReferenceEntryPtr = encode_context_ptr->picture_decision_pa_reference_queue[paReferenceQueueIndex];

                                    // Calculate the Ref POC
                                    refPoc = POC_CIRCULAR_ADD(
                                        picture_control_set_ptr->picture_number,
                                        -inputEntryPtr->list1Ptr->referenceList/*,
                                        sequence_control_set_ptr->bits_for_picture_order_count*/);
                                    picture_control_set_ptr->ref_pa_pcs_array[REF_LIST_1] = paReferenceEntryPtr->pPcsPtr;
                                    // Set the Reference Object
                                    picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_1] = paReferenceEntryPtr->inputObjectPtr;
                                    picture_control_set_ptr->ref_pic_poc_array[REF_LIST_1] = refPoc;

                                    // Increment the PA Reference's liveCount by the number of tiles in the input picture
                                    EbObjectIncLiveCount(
                                        paReferenceEntryPtr->inputObjectPtr,
                                        1);

                                    ((EbPaReferenceObject_t*)picture_control_set_ptr->ref_pa_pic_ptr_array[REF_LIST_1]->objectPtr)->pPcsPtr = paReferenceEntryPtr->pPcsPtr;

                                    EbObjectIncLiveCount(
                                        paReferenceEntryPtr->pPcsPtr->p_pcs_wrapper_ptr,
                                        1);

                                    --paReferenceEntryPtr->dependentCount;
                                }
                            }

                            // SB Loop to reset similarColocatedLcu Array
                            uint16_t *variancePtr;
                            uint32_t  NullVarCnt = 0;
                            uint32_t  varLcuCnt = 0;
                            uint32_t  lcuCodingOrder;
                            uint32_t  sb_origin_x;
                            uint32_t  sb_origin_y;

                            if ((picture_control_set_ptr->slice_type == P_SLICE) || (picture_control_set_ptr->slice_type == B_SLICE)) {
                                picture_width_in_sb = (sequence_control_set_ptr->luma_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
                                for (lcuCodingOrder = 0; lcuCodingOrder < picture_control_set_ptr->sb_total_count; ++lcuCodingOrder) {
                                    sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
                                    sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
                                    variancePtr = picture_control_set_ptr->variance[lcuCodingOrder];
                                    //do it only for complete 64x64 blocks
                                    if (sb_origin_x + 64 <= picture_control_set_ptr->enhanced_picture_ptr->width && sb_origin_y + 64 <= picture_control_set_ptr->enhanced_picture_ptr->height)
                                    {
                                        NullVarCnt += (variancePtr[0] == 0) ? 1 : 0;
                                        varLcuCnt++;
                                    }
                                }
                            }

                            {

                                picture_control_set_ptr->intensity_transition_flag = EB_FALSE;
                                if (picture_control_set_ptr->ref_list0_count) {
                                    picture_control_set_ptr->scene_transition_flag[REF_LIST_0] = EB_FALSE;

                                }

                                if (picture_control_set_ptr->ref_list1_count) {

                                    picture_control_set_ptr->scene_transition_flag[REF_LIST_1] = EB_FALSE;
                                }
                            }


                            // Initialize Segments
                            picture_control_set_ptr->me_segments_column_count = (uint8_t)(sequence_control_set_ptr->me_segment_column_count_array[picture_control_set_ptr->temporal_layer_index]);
                            picture_control_set_ptr->me_segments_row_count = (uint8_t)(sequence_control_set_ptr->me_segment_row_count_array[picture_control_set_ptr->temporal_layer_index]);
                            picture_control_set_ptr->me_segments_total_count = (uint16_t)(picture_control_set_ptr->me_segments_column_count  * picture_control_set_ptr->me_segments_row_count);
                            picture_control_set_ptr->me_segments_completion_mask = 0;

                            // Post the results to the ME processes
                            {
                                uint32_t segment_index;

                                for (segment_index = 0; segment_index < picture_control_set_ptr->me_segments_total_count; ++segment_index)
                                {
                                    // Get Empty Results Object
                                    EbGetEmptyObject(
                                        context_ptr->pictureDecisionResultsOutputFifoPtr,
                                        &outputResultsWrapperPtr);

                                    outputResultsPtr = (PictureDecisionResults_t*)outputResultsWrapperPtr->objectPtr;

                                    outputResultsPtr->pictureControlSetWrapperPtr = encode_context_ptr->pre_assignment_buffer[pictureIndex];

                                    outputResultsPtr->segment_index = segment_index;

                                    // Post the Full Results Object
                                    EbPostFullObject(outputResultsWrapperPtr);
                                }
                            }

                            if (pictureIndex == context_ptr->miniGopEndIndex[miniGopIndex]) {

                                // Increment the Decode Base Number
                                encode_context_ptr->decode_base_number += context_ptr->miniGopLength[miniGopIndex];
                            }

                            if (pictureIndex == encode_context_ptr->pre_assignment_buffer_count - 1) {

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
                    if ((inputEntryPtr->dependentCount == 0) &&
                        (inputEntryPtr->inputObjectPtr)) {
                        EbReleaseObject(inputEntryPtr->pPcsPtr->p_pcs_wrapper_ptr);
                        // Release the nominal liveCount value
                        EbReleaseObject(inputEntryPtr->inputObjectPtr);
                        inputEntryPtr->inputObjectPtr = (EbObjectWrapper_t*)EB_NULL;
                    }

                    // Increment the HeadIndex if the head is null
                    encode_context_ptr->picture_decision_pa_reference_queue_head_index =
                        (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->inputObjectPtr) ? encode_context_ptr->picture_decision_pa_reference_queue_head_index :
                        (encode_context_ptr->picture_decision_pa_reference_queue_head_index == PICTURE_DECISION_PA_REFERENCE_QUEUE_MAX_DEPTH - 1) ? 0
                        : encode_context_ptr->picture_decision_pa_reference_queue_head_index + 1;

                    CHECK_REPORT_ERROR(
                        (((encode_context_ptr->picture_decision_pa_reference_queue_head_index != encode_context_ptr->picture_decision_pa_reference_queue_tail_index) ||
                        (encode_context_ptr->picture_decision_pa_reference_queue[encode_context_ptr->picture_decision_pa_reference_queue_head_index]->inputObjectPtr == EB_NULL))),
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
        EbReleaseObject(inputResultsWrapperPtr);
    }

    return EB_NULL;
}
