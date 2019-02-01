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

#include "EbDefinitions.h"
#include "EbRateControlProcess.h"
#include "EbSystemResourceManager.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbUtility.h"
#include "EbErrorCodes.h"

#include "EbRateControlResults.h"
#include "EbRateControlTasks.h"
#include "RateControlModel.h"

#if V2_QP_SCALING
static const uint8_t DEFAULT_QP_OFFSET_LAYER_ARRAY[MAX_TEMPORAL_LAYERS] = {
     1, 2, 3, 4, 5, 6
};
#endif
#if NEW_QP_SCALING || P2_NEW_QP_SCALING
static const uint8_t QP_OFFSET_LAYER_ARRAY_BDRATE[MAX_TEMPORAL_LAYERS] = {
    1, 3, 5, 5, 6, 7
};
#endif

static uint8_t QP_OFFSET_LAYER_ARRAY[MAX_TEMPORAL_LAYERS] =
{
    1, 2, 4, 5, 6, 7
};

/*****************************
* Internal Typedefs
*****************************/
void RateControlLayerReset(
    RateControlLayerContext_t   *rateControlLayerPtr,
    PictureControlSet_t         *picture_control_set_ptr,
    RateControlContext_t        *rateControlContextPtr,
    uint32_t                       pictureAreaInPixel,
    EbBool                      wasUsed)
{

    SequenceControlSet_t *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    uint32_t                sliceNum;
    uint32_t                temporal_layer_index;
    uint64_t                totalFrameInInterval;
    uint64_t                sumBitsPerSw = 0;

    rateControlLayerPtr->target_bit_rate = picture_control_set_ptr->parent_pcs_ptr->target_bit_rate*RATE_PERCENTAGE_LAYER_ARRAY[sequence_control_set_ptr->static_config.hierarchical_levels][rateControlLayerPtr->temporalIndex] / 100;
    // update this based on temporal layers
    rateControlLayerPtr->frame_rate = sequence_control_set_ptr->frame_rate;

    totalFrameInInterval = sequence_control_set_ptr->static_config.intra_period_length + 1;

    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
        if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
            totalFrameInInterval = 0;
            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
                rateControlContextPtr->frames_in_interval[temporal_layer_index] = picture_control_set_ptr->parent_pcs_ptr->frames_in_interval[temporal_layer_index];
                totalFrameInInterval += picture_control_set_ptr->parent_pcs_ptr->frames_in_interval[temporal_layer_index];
                sumBitsPerSw += picture_control_set_ptr->parent_pcs_ptr->bits_per_sw_per_layer[temporal_layer_index];
            }
#if ADAPTIVE_PERCENTAGE
            rateControlLayerPtr->target_bit_rate = picture_control_set_ptr->parent_pcs_ptr->target_bit_rate* picture_control_set_ptr->parent_pcs_ptr->bits_per_sw_per_layer[rateControlLayerPtr->temporalIndex] / sumBitsPerSw;
#endif
        }
    }


    if (rateControlLayerPtr->temporalIndex == 0) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 5;
    }
    else if (rateControlLayerPtr->temporalIndex == 1) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 5;
    }
    else if (rateControlLayerPtr->temporalIndex == 2) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 4;
    }
    else if (rateControlLayerPtr->temporalIndex == 3) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 3;
    }
    else if (rateControlLayerPtr->temporalIndex == 4) {
        rateControlLayerPtr->coeffAveragingWeight1 = 5;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 2;
    }
    else if (rateControlLayerPtr->temporalIndex == 5) {
        rateControlLayerPtr->coeffAveragingWeight1 = 3;
        rateControlLayerPtr->frame_rate = rateControlLayerPtr->frame_rate >> 1;
    }
    if (sequence_control_set_ptr->static_config.intra_period_length != -1) {
        rateControlLayerPtr->frame_rate = sequence_control_set_ptr->frame_rate * rateControlContextPtr->frames_in_interval[rateControlLayerPtr->temporalIndex] / totalFrameInInterval;
    }

    rateControlLayerPtr->coeffAveragingWeight2 = 16 - rateControlLayerPtr->coeffAveragingWeight1;
    if (rateControlLayerPtr->frame_rate == 0) { // no frame in that layer
        rateControlLayerPtr->frame_rate = 1 << RC_PRECISION;
    }
    rateControlLayerPtr->channelBitRate = (((rateControlLayerPtr->target_bit_rate << (2 * RC_PRECISION)) / rateControlLayerPtr->frame_rate) + RC_PRECISION_OFFSET) >> RC_PRECISION;
    rateControlLayerPtr->channelBitRate = (uint64_t)MAX((int64_t)1, (int64_t)rateControlLayerPtr->channelBitRate);
    rateControlLayerPtr->ecBitConstraint = rateControlLayerPtr->channelBitRate;


    // This is only for the initial frame, because the feedback is from packetization now and all of these are considered
    // considering the bits for slice header
    // *Note - only one-slice-per picture is supported for UHD
    sliceNum = 1;

    rateControlLayerPtr->ecBitConstraint -= SLICE_HEADER_BITS_NUM * sliceNum;

    rateControlLayerPtr->ecBitConstraint = MAX(1, rateControlLayerPtr->ecBitConstraint);

    rateControlLayerPtr->previousBitConstraint = rateControlLayerPtr->channelBitRate;
    rateControlLayerPtr->bitConstraint = rateControlLayerPtr->channelBitRate;
    rateControlLayerPtr->difTotalAndEcBits = 0;

    rateControlLayerPtr->frameSameSADMinQpCount = 0;
    rateControlLayerPtr->maxQp = picture_control_set_ptr->picture_qp;

    rateControlLayerPtr->alpha = 1 << (RC_PRECISION - 1);
    {
        if (!wasUsed) {


            rateControlLayerPtr->sameSADCount = 0;

            rateControlLayerPtr->kCoeff = 3 << RC_PRECISION;
            rateControlLayerPtr->previousKCoeff = 3 << RC_PRECISION;

            rateControlLayerPtr->cCoeff = (rateControlLayerPtr->channelBitRate << (2 * RC_PRECISION)) / pictureAreaInPixel / CCOEFF_INIT_FACT;
            rateControlLayerPtr->previousCCoeff = (rateControlLayerPtr->channelBitRate << (2 * RC_PRECISION)) / pictureAreaInPixel / CCOEFF_INIT_FACT;

            // These are for handling Pred structure 2, when for higher temporal layer, frames can arrive in different orders
            // They should be modifed in a way that gets these from previous layers
            rateControlLayerPtr->previousFrameQp = 32;
            rateControlLayerPtr->previousFrameBitActual = 1200;
            rateControlLayerPtr->previousFrameQuantizedCoeffBitActual = 1000;
            rateControlLayerPtr->previousFrameSadMe = 10000000;
            rateControlLayerPtr->previousFrameQp = picture_control_set_ptr->picture_qp;
            rateControlLayerPtr->deltaQpFraction = 0;
            rateControlLayerPtr->previousFrameAverageQp = picture_control_set_ptr->picture_qp;
            rateControlLayerPtr->previousCalculatedFrameQp = picture_control_set_ptr->picture_qp;
            rateControlLayerPtr->calculatedFrameQp = picture_control_set_ptr->picture_qp;
            rateControlLayerPtr->criticalStates = 0;
        }
        else {
            rateControlLayerPtr->sameSADCount = 0;
            rateControlLayerPtr->criticalStates = 0;
        }
    }
}


void RateControlLayerResetPart2(
    RateControlLayerContext_t   *rateControlLayerPtr,
    PictureControlSet_t         *picture_control_set_ptr)
{

    // update this based on temporal layers

    rateControlLayerPtr->maxQp = (uint32_t)CLIP3(0, 63, (int32_t)(picture_control_set_ptr->picture_qp + QP_OFFSET_LAYER_ARRAY[rateControlLayerPtr->temporalIndex]));;
    {

        // These are for handling Pred structure 2, when for higher temporal layer, frames can arrive in different orders
        // They should be modifed in a way that gets these from previous layers
        rateControlLayerPtr->previousFrameQp = rateControlLayerPtr->maxQp;
        rateControlLayerPtr->previousFrameAverageQp = rateControlLayerPtr->maxQp;
        rateControlLayerPtr->previousCalculatedFrameQp = rateControlLayerPtr->maxQp;
        rateControlLayerPtr->calculatedFrameQp = rateControlLayerPtr->maxQp;
    }
}

EbErrorType HighLevelRateControlContextCtor(
    HighLevelRateControlContext_t   **entryDblPtr) {

    HighLevelRateControlContext_t *entryPtr;
    EB_MALLOC(HighLevelRateControlContext_t*, entryPtr, sizeof(HighLevelRateControlContext_t), EB_N_PTR);
    *entryDblPtr = entryPtr;

    return EB_ErrorNone;
}


EbErrorType RateControlLayerContextCtor(
    RateControlLayerContext_t   **entryDblPtr) {

    RateControlLayerContext_t *entryPtr;
    EB_MALLOC(RateControlLayerContext_t*, entryPtr, sizeof(RateControlLayerContext_t), EB_N_PTR);

    *entryDblPtr = entryPtr;

    entryPtr->firstFrame = 1;
    entryPtr->firstNonIntraFrame = 1;
    entryPtr->feedbackArrived = EB_FALSE;

    return EB_ErrorNone;
}



EbErrorType RateControlIntervalParamContextCtor(
    RateControlIntervalParamContext_t   **entryDblPtr) {

    uint32_t temporalIndex;
    EbErrorType return_error = EB_ErrorNone;
    RateControlIntervalParamContext_t *entryPtr;
    EB_MALLOC(RateControlIntervalParamContext_t*, entryPtr, sizeof(RateControlIntervalParamContext_t), EB_N_PTR);

    *entryDblPtr = entryPtr;

    entryPtr->inUse = EB_FALSE;
    entryPtr->wasUsed = EB_FALSE;
    entryPtr->lastGop = EB_FALSE;
    entryPtr->processedFramesNumber = 0;
    EB_MALLOC(RateControlLayerContext_t**, entryPtr->rateControlLayerArray, sizeof(RateControlLayerContext_t*)*EB_MAX_TEMPORAL_LAYERS, EB_N_PTR);

    for (temporalIndex = 0; temporalIndex < EB_MAX_TEMPORAL_LAYERS; temporalIndex++) {
        return_error = RateControlLayerContextCtor(&entryPtr->rateControlLayerArray[temporalIndex]);
        entryPtr->rateControlLayerArray[temporalIndex]->temporalIndex = temporalIndex;
        entryPtr->rateControlLayerArray[temporalIndex]->frame_rate = 1 << RC_PRECISION;
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    entryPtr->min_target_rate_assigned = EB_FALSE;

    entryPtr->intraFramesQp = 0;
    entryPtr->nextGopIntraFrameQp = 0;
    entryPtr->firstPicPredBits = 0;
    entryPtr->firstPicActualBits = 0;
    entryPtr->firstPicPredQp = 0;
    entryPtr->firstPicActualQp = 0;
    entryPtr->firstPicActualQpAssigned = EB_FALSE;
    entryPtr->scene_change_in_gop = EB_FALSE;
    entryPtr->extraApBitRatioI = 0;

    return EB_ErrorNone;
}

EbErrorType RateControlCodedFramesStatsContextCtor(
    CodedFramesStatsEntry_t   **entryDblPtr,
    uint64_t                      picture_number) {

    CodedFramesStatsEntry_t *entryPtr;
    EB_MALLOC(CodedFramesStatsEntry_t*, entryPtr, sizeof(CodedFramesStatsEntry_t), EB_N_PTR);

    *entryDblPtr = entryPtr;

    entryPtr->picture_number = picture_number;
    entryPtr->frameTotalBitActual = -1;

    return EB_ErrorNone;
}


EbErrorType RateControlContextCtor(
    RateControlContext_t   **context_dbl_ptr,
    EbFifo_t                *rateControlInputTasksFifoPtr,
    EbFifo_t                *rateControlOutputResultsFifoPtr,
    int32_t                   intra_period_length)
{
    uint32_t temporalIndex;
    uint32_t intervalIndex;

#if OVERSHOOT_STAT_PRINT
    uint32_t pictureIndex;
#endif

    EbErrorType return_error = EB_ErrorNone;
    RateControlContext_t *context_ptr;
    EB_MALLOC(RateControlContext_t*, context_ptr, sizeof(RateControlContext_t), EB_N_PTR);

    *context_dbl_ptr = context_ptr;

    context_ptr->rateControlInputTasksFifoPtr = rateControlInputTasksFifoPtr;
    context_ptr->rateControlOutputResultsFifoPtr = rateControlOutputResultsFifoPtr;

    // High level RC
    return_error = HighLevelRateControlContextCtor(
        &context_ptr->highLevelRateControlPtr);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    for (temporalIndex = 0; temporalIndex < EB_MAX_TEMPORAL_LAYERS; temporalIndex++) {
        context_ptr->frames_in_interval[temporalIndex] = 0;
    }

    EB_MALLOC(RateControlIntervalParamContext_t**, context_ptr->rateControlParamQueue, sizeof(RateControlIntervalParamContext_t*)*PARALLEL_GOP_MAX_NUMBER, EB_N_PTR);

    context_ptr->rateControlParamQueueHeadIndex = 0;
    for (intervalIndex = 0; intervalIndex < PARALLEL_GOP_MAX_NUMBER; intervalIndex++) {
        return_error = RateControlIntervalParamContextCtor(
            &context_ptr->rateControlParamQueue[intervalIndex]);
        context_ptr->rateControlParamQueue[intervalIndex]->firstPoc = (intervalIndex*(uint32_t)(intra_period_length + 1));
        context_ptr->rateControlParamQueue[intervalIndex]->lastPoc = ((intervalIndex + 1)*(uint32_t)(intra_period_length + 1)) - 1;
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

#if OVERSHOOT_STAT_PRINT
    context_ptr->codedFramesStatQueueHeadIndex = 0;
    context_ptr->codedFramesStatQueueTailIndex = 0;
    EB_MALLOC(CodedFramesStatsEntry_t**, context_ptr->codedFramesStatQueue, sizeof(CodedFramesStatsEntry_t*)*CODED_FRAMES_STAT_QUEUE_MAX_DEPTH, EB_N_PTR);

    for (pictureIndex = 0; pictureIndex < CODED_FRAMES_STAT_QUEUE_MAX_DEPTH; ++pictureIndex) {
        return_error = RateControlCodedFramesStatsContextCtor(
            &context_ptr->codedFramesStatQueue[pictureIndex],
            pictureIndex);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    context_ptr->maxBitActualPerSw = 0;
    context_ptr->maxBitActualPerGop = 0;
#endif

    context_ptr->baseLayerFramesAvgQp = 0;
    context_ptr->baseLayerIntraFramesAvgQp = 0;


    context_ptr->intraCoefRate = 4;
    context_ptr->extraBits = 0;
    context_ptr->extraBitsGen = 0;
    context_ptr->maxRateAdjustDeltaQP = 0;


    return EB_ErrorNone;
}
void HighLevelRcInputPictureMode2(
    PictureParentControlSet_t         *picture_control_set_ptr,
    SequenceControlSet_t              *sequence_control_set_ptr,
    EncodeContext_t                   *encode_context_ptr,
    RateControlContext_t              *context_ptr,
    HighLevelRateControlContext_t     *highLevelRateControlPtr)
{

    EbBool                             end_of_sequence_flag = EB_TRUE;

    HlRateControlHistogramEntry_t      *hlRateControlHistogramPtrTemp;
    // Queue variables
    uint32_t                             queueEntryIndexTemp;
    uint32_t                             queueEntryIndexTemp2;
    uint32_t                             queueEntryIndexHeadTemp;


    uint64_t                              minLaBitDistance;
    uint32_t                              selectedRefQpTableIndex;
    uint32_t                              selectedRefQp;
#if RC_UPDATE_TARGET_RATE
    uint32_t                              selectedOrgRefQp;
#endif
    uint32_t                                previous_selected_ref_qp = encode_context_ptr->previous_selected_ref_qp;
    uint64_t                                max_coded_poc = encode_context_ptr->max_coded_poc;
    uint32_t                                max_coded_poc_selected_ref_qp = encode_context_ptr->max_coded_poc_selected_ref_qp;


    uint32_t                              refQpIndex;
    uint32_t                              refQpIndexTemp;
    uint32_t                              refQpTableIndex;

    uint32_t                              areaInPixel;
    uint32_t                              numOfFullLcus;
    uint32_t                              qpSearchMin;
    uint32_t                              qpSearchMax;
    int32_t                              qpStep = 1;
    EbBool                             bestQpFound;
    uint32_t                              temporal_layer_index;
    EbBool                             tables_updated;

    uint64_t                              bitConstraintPerSw = 0;

    RateControlTables_t                    *rateControlTablesPtr;
    EB_Bit_Number                        *sadBitsArrayPtr;
    EB_Bit_Number                        *intraSadBitsArrayPtr;
    uint32_t                               pred_bits_ref_qp;

    for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
        picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = 0;
    }
    picture_control_set_ptr->total_bits_per_gop = 0;

    areaInPixel = sequence_control_set_ptr->luma_width * sequence_control_set_ptr->luma_height;;

    EbBlockOnMutex(sequence_control_set_ptr->encode_context_ptr->rate_table_update_mutex);

    tables_updated = sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array_updated;
    picture_control_set_ptr->percentage_updated = EB_FALSE;

    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {

        // Increamenting the head of the hl_rate_control_historgram_queue and clean up the entores
        hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]);
        while ((hlRateControlHistogramPtrTemp->lifeCount == 0) && hlRateControlHistogramPtrTemp->passedToHlrc) {

            EbBlockOnMutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            // Reset the Reorder Queue Entry
            hlRateControlHistogramPtrTemp->picture_number += INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
            hlRateControlHistogramPtrTemp->lifeCount = -1;
            hlRateControlHistogramPtrTemp->passedToHlrc = EB_FALSE;
            hlRateControlHistogramPtrTemp->isCoded = EB_FALSE;
            hlRateControlHistogramPtrTemp->totalNumBitsCoded = 0;

            // Increment the Reorder Queue head Ptr
            encode_context_ptr->hl_rate_control_historgram_queue_head_index =
                (encode_context_ptr->hl_rate_control_historgram_queue_head_index == HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->hl_rate_control_historgram_queue_head_index + 1;
            EbReleaseMutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            hlRateControlHistogramPtrTemp = encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index];

        }
        // For the case that number of frames in the sliding window is less than size of the look ahead or intra Refresh. i.e. end of sequence
        if ((picture_control_set_ptr->frames_in_sw < MIN(sequence_control_set_ptr->static_config.look_ahead_distance + 1, (uint32_t)sequence_control_set_ptr->intra_period_length + 1))) {

            selectedRefQp = max_coded_poc_selected_ref_qp;

            // Update the QP for the sliding window based on the status of RC
            if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize << 3))) {
                selectedRefQp = (uint32_t)MAX((int32_t)selectedRefQp - 2, 0);
            }
            else if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize << 2))) {
                selectedRefQp = (uint32_t)MAX((int32_t)selectedRefQp - 1, 0);
            }
            if ((context_ptr->extraBitsGen < -(int64_t)(context_ptr->virtualBufferSize << 2))) {
                selectedRefQp += 2;
            }
            else if ((context_ptr->extraBitsGen < -(int64_t)(context_ptr->virtualBufferSize << 1))) {
                selectedRefQp += 1;
            }

            if ((picture_control_set_ptr->frames_in_sw < (uint32_t)(sequence_control_set_ptr->intra_period_length + 1)) &&
                (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0)) {
                selectedRefQp = (uint32_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    selectedRefQp + 1);
            }

            queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
            queueEntryIndexHeadTemp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
            queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                queueEntryIndexHeadTemp;

            queueEntryIndexTemp = queueEntryIndexHeadTemp;
            {

                hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexTemp]);
                refQpIndexTemp = selectedRefQp + QP_OFFSET_LAYER_ARRAY[hlRateControlHistogramPtrTemp->temporal_layer_index];
                refQpIndexTemp = (uint32_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    refQpIndexTemp);

                if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                    refQpIndexTemp = (uint32_t)MAX((int32_t)refQpIndexTemp + RC_INTRA_QP_OFFSET, 0);
                }

                hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = 0;
                rateControlTablesPtr = &encode_context_ptr->rate_control_tables_array[refQpIndexTemp];
                sadBitsArrayPtr = rateControlTablesPtr->sadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                intraSadBitsArrayPtr = rateControlTablesPtr->intraSadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                pred_bits_ref_qp = 0;
                numOfFullLcus = 0;

                if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                    // Loop over block in the frame and calculated the predicted bits at reg QP
                    {
                        unsigned i;
                        uint32_t accum = 0;
                        for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                        {
                            accum += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);
                        }

                        pred_bits_ref_qp = accum;
                        numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                    }
                    hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                }

                else {
                    {
                        unsigned i;
                        uint32_t accum = 0;
                        for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
                        {
                            accum += (uint32_t)(hlRateControlHistogramPtrTemp->me_distortion_histogram[i] * sadBitsArrayPtr[i]);
                        }

                        pred_bits_ref_qp = accum;
                        numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;

                    }
                    hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                }

                // Scale for in complete
                //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the picture boundries
                hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] * (uint64_t)areaInPixel / (numOfFullLcus << 12);

                // Store the pred_bits_ref_qp for the first frame in the window to PCS
                picture_control_set_ptr->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];

            }
        }
        else {
            // Loop over the QPs and find the best QP
            minLaBitDistance = MAX_UNSIGNED_VALUE;
            qpSearchMin = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                (uint32_t)MAX((int32_t)sequence_control_set_ptr->qp - 20, 0));

            qpSearchMax = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                sequence_control_set_ptr->qp + 20);

            for (refQpTableIndex = qpSearchMin; refQpTableIndex < qpSearchMax; refQpTableIndex++) {
                highLevelRateControlPtr->predBitsRefQpPerSw[refQpTableIndex] = 0;
            }

            bitConstraintPerSw = highLevelRateControlPtr->bitConstraintPerSw * picture_control_set_ptr->frames_in_sw / (sequence_control_set_ptr->static_config.look_ahead_distance + 1);

            // Update the target rate for the sliding window based on the status of RC
            if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize * 10))) {
                bitConstraintPerSw = bitConstraintPerSw * 130 / 100;
            }
            else if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize << 3))) {
                bitConstraintPerSw = bitConstraintPerSw * 120 / 100;
            }
            else if ((context_ptr->extraBitsGen > (int64_t)(context_ptr->virtualBufferSize << 2))) {
                bitConstraintPerSw = bitConstraintPerSw * 110 / 100;
            }
            if ((context_ptr->extraBitsGen < -(int64_t)(context_ptr->virtualBufferSize << 3))) {
                bitConstraintPerSw = bitConstraintPerSw * 80 / 100;
            }
            else if ((context_ptr->extraBitsGen < -(int64_t)(context_ptr->virtualBufferSize << 2))) {
                bitConstraintPerSw = bitConstraintPerSw * 90 / 100;
            }

            // Loop over proper QPs and find the Predicted bits for that QP. Find the QP with the closest total predicted rate to target bits for the sliding window.
            previous_selected_ref_qp = CLIP3(
                qpSearchMin,
                qpSearchMax,
                previous_selected_ref_qp);
            refQpTableIndex = previous_selected_ref_qp;
            selectedRefQpTableIndex = refQpTableIndex;
            selectedRefQp = refQpListTable[selectedRefQpTableIndex];
            bestQpFound = EB_FALSE;
            while (refQpTableIndex >= qpSearchMin && refQpTableIndex <= qpSearchMax && !bestQpFound) {

                refQpIndex = CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    refQpListTable[refQpTableIndex]);
                highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] = 0;

                // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
                queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queueEntryIndexHeadTemp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queueEntryIndexHeadTemp;

                queueEntryIndexTemp = queueEntryIndexHeadTemp;
                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    queueEntryIndexTemp <= queueEntryIndexHeadTemp + sequence_control_set_ptr->static_config.look_ahead_distance) {

                    queueEntryIndexTemp2 = (queueEntryIndexTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queueEntryIndexTemp;
                    hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexTemp2]);

                    refQpIndexTemp = refQpIndex + QP_OFFSET_LAYER_ARRAY[hlRateControlHistogramPtrTemp->temporal_layer_index];
                    refQpIndexTemp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        refQpIndexTemp);

                    if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                        refQpIndexTemp = (uint32_t)MAX((int32_t)refQpIndexTemp + RC_INTRA_QP_OFFSET, 0);
                    }

                    hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = 0;

                    if (refQpTableIndex == previous_selected_ref_qp) {
                        hlRateControlHistogramPtrTemp->lifeCount--;
                    }
                    if (hlRateControlHistogramPtrTemp->isCoded) {
                        // If the frame is already coded, use the actual number of bits
                        hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->totalNumBitsCoded;
                    }
                    else {
                        rateControlTablesPtr = &encode_context_ptr->rate_control_tables_array[refQpIndexTemp];
                        sadBitsArrayPtr = rateControlTablesPtr->sadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                        intraSadBitsArrayPtr = rateControlTablesPtr->intraSadBitsArray[0];
                        pred_bits_ref_qp = 0;
                        numOfFullLcus = 0;

                        if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                            // Loop over block in the frame and calculated the predicted bits at reg QP
                            unsigned i;
                            uint32_t accum = 0;
                            for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                            {
                                accum += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);
                            }

                            pred_bits_ref_qp = accum;
                            numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                            hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                        }
                        else {
                            unsigned i;
                            uint32_t accum = 0;
                            uint32_t accumIntra = 0;
                            for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
                            {
                                accum += (uint32_t)(hlRateControlHistogramPtrTemp->me_distortion_histogram[i] * sadBitsArrayPtr[i]);
                                accumIntra += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);

                            }
                            if (accum > accumIntra * 3)
                                pred_bits_ref_qp = accumIntra;
                            else
                                pred_bits_ref_qp = accum;
                            numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                            hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                        }

                        // Scale for in complete LCSs
                        //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the picture boundries
                        hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] * (uint64_t)areaInPixel / (numOfFullLcus << 12);

                    }
                    highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] += hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    // Store the pred_bits_ref_qp for the first frame in the window to PCS
                    if (queueEntryIndexHeadTemp == queueEntryIndexTemp2)
                        picture_control_set_ptr->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];

                    end_of_sequence_flag = hlRateControlHistogramPtrTemp->end_of_sequence_flag;
                    queueEntryIndexTemp++;
                }

                if (minLaBitDistance >= (uint64_t)ABS((int64_t)highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] - (int64_t)bitConstraintPerSw)) {
                    minLaBitDistance = (uint64_t)ABS((int64_t)highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] - (int64_t)bitConstraintPerSw);
                    selectedRefQpTableIndex = refQpTableIndex;
                    selectedRefQp = refQpIndex;
                }
                else {
                    bestQpFound = EB_TRUE;
                }

                if (refQpTableIndex == previous_selected_ref_qp) {
                    if (highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] > bitConstraintPerSw) {
                        qpStep = +1;
                    }
                    else {
                        qpStep = -1;
                    }
                }
                refQpTableIndex = (uint32_t)(refQpTableIndex + qpStep);

            }
        }

#if RC_UPDATE_TARGET_RATE
        selectedOrgRefQp = selectedRefQp;
        if (sequence_control_set_ptr->intra_period_length != -1 && picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0 &&
            (int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) {
            if (picture_control_set_ptr->picture_number > 0) {
                picture_control_set_ptr->intra_selected_org_qp = (uint8_t)selectedRefQp;
            }
            else {
                selectedOrgRefQp = selectedRefQp + 1;
                selectedRefQp = selectedRefQp + 1;
            }
            refQpIndex = selectedRefQp;
            highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] = 0;

            if (highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] == 0) {

                // Finding the predicted bits for each frame in the sliding window at the reference Qp(s)
                //queueEntryIndexTemp = encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queueEntryIndexHeadTemp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queueEntryIndexHeadTemp;

                queueEntryIndexTemp = queueEntryIndexHeadTemp;

                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    //queueEntryIndexTemp <= encode_context_ptr->hl_rate_control_historgram_queue_head_index+sequence_control_set_ptr->static_config.look_ahead_distance){
                    queueEntryIndexTemp <= queueEntryIndexHeadTemp + sequence_control_set_ptr->static_config.look_ahead_distance) {

                    queueEntryIndexTemp2 = (queueEntryIndexTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queueEntryIndexTemp;
                    hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexTemp2]);


                    refQpIndexTemp = refQpIndex + QP_OFFSET_LAYER_ARRAY[hlRateControlHistogramPtrTemp->temporal_layer_index];
                    refQpIndexTemp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        refQpIndexTemp);

                    if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                        refQpIndexTemp = (uint32_t)MAX((int32_t)refQpIndexTemp + RC_INTRA_QP_OFFSET, 0);
                    }

                    hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = 0;

                    if (hlRateControlHistogramPtrTemp->isCoded) {
                        hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->totalNumBitsCoded;
                    }
                    else {
                        rateControlTablesPtr = &encode_context_ptr->rate_control_tables_array[refQpIndexTemp];
                        sadBitsArrayPtr = rateControlTablesPtr->sadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                        intraSadBitsArrayPtr = rateControlTablesPtr->intraSadBitsArray[hlRateControlHistogramPtrTemp->temporal_layer_index];
                        pred_bits_ref_qp = 0;

                        numOfFullLcus = 0;

                        if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                            // Loop over block in the frame and calculated the predicted bits at reg QP

                            {
                                unsigned i;
                                uint32_t accum = 0;
                                for (i = 0; i < NUMBER_OF_INTRA_SAD_INTERVALS; ++i)
                                {
                                    accum += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);
                                }

                                pred_bits_ref_qp = accum;
                                numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                            }
                            hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                        }

                        else {
                            unsigned i;
                            uint32_t accum = 0;
                            uint32_t accumIntra = 0;
                            for (i = 0; i < NUMBER_OF_SAD_INTERVALS; ++i)
                            {
                                accum += (uint32_t)(hlRateControlHistogramPtrTemp->me_distortion_histogram[i] * sadBitsArrayPtr[i]);
                                accumIntra += (uint32_t)(hlRateControlHistogramPtrTemp->ois_distortion_histogram[i] * intraSadBitsArrayPtr[i]);

                            }
                            if (accum > accumIntra * 3)
                                pred_bits_ref_qp = accumIntra;
                            else
                                pred_bits_ref_qp = accum;
                            numOfFullLcus = hlRateControlHistogramPtrTemp->full_sb_count;
                            hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] += pred_bits_ref_qp;
                        }

                        // Scale for in complete
                        //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the picture boundries
                        hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp] * (uint64_t)areaInPixel / (numOfFullLcus << 12);

                    }
                    highLevelRateControlPtr->predBitsRefQpPerSw[refQpIndex] += hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    // Store the pred_bits_ref_qp for the first frame in the window to PCS
                    //  if(encode_context_ptr->hl_rate_control_historgram_queue_head_index == queueEntryIndexTemp2)
                    if (queueEntryIndexHeadTemp == queueEntryIndexTemp2)
                        picture_control_set_ptr->pred_bits_ref_qp[refQpIndexTemp] = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];

                    end_of_sequence_flag = hlRateControlHistogramPtrTemp->end_of_sequence_flag;
                    queueEntryIndexTemp++;
                }
            }
        }
#endif
        picture_control_set_ptr->tables_updated = tables_updated;
        EbBool expensiveISlice = EB_FALSE;
        // Looping over the window to find the percentage of bit allocation in each layer
        if ((sequence_control_set_ptr->intra_period_length != -1) &&
            ((int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length) &&
            ((int32_t)picture_control_set_ptr->frames_in_sw > sequence_control_set_ptr->intra_period_length)) {
            uint64_t iSliceBits = 0;

            if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {

                queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
                queueEntryIndexHeadTemp += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
                queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                    queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                    queueEntryIndexHeadTemp;

                queueEntryIndexTemp = queueEntryIndexHeadTemp;

                // This is set to false, so the last frame would go inside the loop
                end_of_sequence_flag = EB_FALSE;

                while (!end_of_sequence_flag &&
                    queueEntryIndexTemp <= queueEntryIndexHeadTemp + sequence_control_set_ptr->intra_period_length) {

                    queueEntryIndexTemp2 = (queueEntryIndexTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH : queueEntryIndexTemp;
                    hlRateControlHistogramPtrTemp = (encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexTemp2]);

                    refQpIndexTemp = selectedRefQp + QP_OFFSET_LAYER_ARRAY[hlRateControlHistogramPtrTemp->temporal_layer_index];
                    refQpIndexTemp = (uint32_t)CLIP3(
                        sequence_control_set_ptr->static_config.min_qp_allowed,
                        sequence_control_set_ptr->static_config.max_qp_allowed,
                        refQpIndexTemp);

                    if (hlRateControlHistogramPtrTemp->slice_type == I_SLICE) {
                        refQpIndexTemp = (uint32_t)MAX((int32_t)refQpIndexTemp + RC_INTRA_QP_OFFSET, 0);
                    }
                    if (queueEntryIndexTemp == queueEntryIndexHeadTemp) {
                        iSliceBits = hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    }
                    picture_control_set_ptr->total_bits_per_gop += hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    picture_control_set_ptr->bits_per_sw_per_layer[hlRateControlHistogramPtrTemp->temporal_layer_index] += hlRateControlHistogramPtrTemp->pred_bits_ref_qp[refQpIndexTemp];
                    picture_control_set_ptr->percentage_updated = EB_TRUE;

                    end_of_sequence_flag = hlRateControlHistogramPtrTemp->end_of_sequence_flag;
                    queueEntryIndexTemp++;
                }
                if (iSliceBits * 100 > 85 * picture_control_set_ptr->total_bits_per_gop) {
                    expensiveISlice = EB_TRUE;
                }
                if (picture_control_set_ptr->total_bits_per_gop == 0) {
                    for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
                        picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = RATE_PERCENTAGE_LAYER_ARRAY[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_index];
                    }
                }
            }
        }
        else {
            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
                picture_control_set_ptr->bits_per_sw_per_layer[temporal_layer_index] = RATE_PERCENTAGE_LAYER_ARRAY[sequence_control_set_ptr->static_config.hierarchical_levels][temporal_layer_index];
            }
        }
        if (expensiveISlice) {
            if (tables_updated) {
                selectedRefQp = (uint32_t)MAX((int32_t)selectedRefQp - 1, 0);
            }
            else {
                selectedRefQp = (uint32_t)MAX((int32_t)selectedRefQp - 3, 0);
            }
            selectedRefQp = (uint32_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                selectedRefQp);
        }
        // Set the QP
        previous_selected_ref_qp = selectedRefQp;
        if (picture_control_set_ptr->picture_number > max_coded_poc && picture_control_set_ptr->temporal_layer_index < 2 && !picture_control_set_ptr->end_of_sequence_region) {

            max_coded_poc = picture_control_set_ptr->picture_number;
            max_coded_poc_selected_ref_qp = previous_selected_ref_qp;
            encode_context_ptr->previous_selected_ref_qp = previous_selected_ref_qp;
            encode_context_ptr->max_coded_poc = max_coded_poc;
            encode_context_ptr->max_coded_poc_selected_ref_qp = max_coded_poc_selected_ref_qp;

        }
        picture_control_set_ptr->best_pred_qp = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            sequence_control_set_ptr->static_config.max_qp_allowed,
            selectedRefQp + QP_OFFSET_LAYER_ARRAY[picture_control_set_ptr->temporal_layer_index]);
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            picture_control_set_ptr->best_pred_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->best_pred_qp + RC_INTRA_QP_OFFSET, 0);
        }
#if RC_UPDATE_TARGET_RATE
        if (picture_control_set_ptr->picture_number == 0) {
            highLevelRateControlPtr->prevIntraSelectedRefQp = selectedRefQp;
            highLevelRateControlPtr->prevIntraOrgSelectedRefQp = selectedRefQp;
        }
        if (sequence_control_set_ptr->intra_period_length != -1) {
            if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                highLevelRateControlPtr->prevIntraSelectedRefQp = selectedRefQp;
                highLevelRateControlPtr->prevIntraOrgSelectedRefQp = selectedOrgRefQp;
            }
        }
#endif
        picture_control_set_ptr->target_bits_best_pred_qp = picture_control_set_ptr->pred_bits_ref_qp[picture_control_set_ptr->best_pred_qp];
        //if (picture_control_set_ptr->slice_type == 2)
        // {
        //printf("\nTID: %d\t", picture_control_set_ptr->temporal_layer_index);
        //printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
        //    picture_control_set_ptr->picture_number,
        //    picture_control_set_ptr->best_pred_qp,
        //    (int32_t)picture_control_set_ptr->target_bits_best_pred_qp,
        //    (int32_t)highLevelRateControlPtr->predBitsRefQpPerSw[selectedRefQp - 1],
        //    (int32_t)highLevelRateControlPtr->predBitsRefQpPerSw[selectedRefQp],
        //    (int32_t)highLevelRateControlPtr->predBitsRefQpPerSw[selectedRefQp + 1],
        //    (int32_t)highLevelRateControlPtr->bitConstraintPerSw,
        //    (int32_t)bitConstraintPerSw,
        //    (int32_t)highLevelRateControlPtr->virtualBufferLevel);
        //}
    }
    EbReleaseMutex(sequence_control_set_ptr->encode_context_ptr->rate_table_update_mutex);
}
void FrameLevelRcInputPictureMode2(
    PictureControlSet_t               *picture_control_set_ptr,
    SequenceControlSet_t              *sequence_control_set_ptr,
    RateControlContext_t              *context_ptr,
    RateControlLayerContext_t         *rateControlLayerPtr,
    RateControlIntervalParamContext_t *rateControlParamPtr,
    uint32_t                             bestOisCuIndex)
{

    RateControlLayerContext_t   *rateControlLayerTempPtr;

    // Tiles
    uint32_t                       pictureAreaInPixel;
    uint32_t                       areaInPixel;

    // SB Loop variables
    uint32_t                       sb_total_count;
    LargestCodingUnit_t         *sb_ptr;
    uint32_t                       lcuCodingOrder;
    uint32_t                       sb_height;
    uint32_t                       sb_width;
    uint64_t                       tempQp;
    uint32_t                       areaInLcus;
    (void)bestOisCuIndex;
    sb_total_count = picture_control_set_ptr->sb_total_count;
    pictureAreaInPixel = sequence_control_set_ptr->luma_height*sequence_control_set_ptr->luma_width;

    if (rateControlLayerPtr->firstFrame == 1) {
        rateControlLayerPtr->firstFrame = 0;
        picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer = 1;
    }
    else {
        picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer = 0;
    }
    if (picture_control_set_ptr->slice_type != I_SLICE) {
        if (rateControlLayerPtr->firstNonIntraFrame == 1) {
            rateControlLayerPtr->firstNonIntraFrame = 0;
            picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer = 1;
        }
        else {
            picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer = 0;
        }
    }

    picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = 0;

    // ***Rate Control***
    areaInLcus = 0;
    areaInPixel = 0;

    for (lcuCodingOrder = 0; lcuCodingOrder < picture_control_set_ptr->sb_total_count; ++lcuCodingOrder) {

        sb_ptr = picture_control_set_ptr->sb_ptr_array[lcuCodingOrder];
        sb_width = (sequence_control_set_ptr->luma_width - sb_ptr->origin_x >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_width - sb_ptr->origin_x;
        sb_height = (sequence_control_set_ptr->luma_height - sb_ptr->origin_y >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_height - sb_ptr->origin_y;

        // This is because of the tile boundry LCUs which do not have correct SAD from ME.
        if ((sb_width == BLOCK_SIZE_64) && (sb_height == BLOCK_SIZE_64)) {
            // add the area of one SB (64x64=4096) to the area of the tile
            areaInPixel += 4096;
            areaInLcus++;
        }
        else {
            // add the area of the SB to the area of the tile
            areaInPixel += sb_width * sb_height;
        }
    }
    rateControlLayerPtr->areaInPixel = areaInPixel;

    if (picture_control_set_ptr->parent_pcs_ptr->first_frame_in_temporal_layer || (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc)) {
        if (sequence_control_set_ptr->static_config.enable_qp_scaling_flag && (picture_control_set_ptr->picture_number != rateControlParamPtr->firstPoc)) {
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (int32_t)sequence_control_set_ptr->static_config.min_qp_allowed,
                (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed,
                (int32_t)(rateControlParamPtr->intraFramesQp + QP_OFFSET_LAYER_ARRAY[picture_control_set_ptr->temporal_layer_index]) - 1 - (RC_INTRA_QP_OFFSET));

        }

        if (picture_control_set_ptr->picture_number == 0) {
            rateControlParamPtr->intraFramesQp = sequence_control_set_ptr->qp;
        }

        if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc) {
            uint32_t temporalLayerIdex;
            rateControlParamPtr->previousVirtualBufferLevel = context_ptr->virtualBufferLevelInitialValue;
            rateControlParamPtr->virtualBufferLevel = context_ptr->virtualBufferLevelInitialValue;
            rateControlParamPtr->extraApBitRatioI = 0;
            if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                rateControlParamPtr->lastPoc = MAX(rateControlParamPtr->firstPoc + picture_control_set_ptr->parent_pcs_ptr->frames_in_sw - 1, rateControlParamPtr->firstPoc);
                rateControlParamPtr->lastGop = EB_TRUE;
            }

            if ((context_ptr->extraBits > (int64_t)(context_ptr->virtualBufferSize >> 8)) ||
                (context_ptr->extraBits < -(int64_t)(context_ptr->virtualBufferSize >> 8))) {

                int64_t extraBitsPerGop = 0;

                if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                    if ((context_ptr->extraBits > (int64_t)(context_ptr->virtualBufferSize << 4)) ||
                        (context_ptr->extraBits < -(int64_t)(context_ptr->virtualBufferSize << 4))) {
                        extraBitsPerGop = context_ptr->extraBits;
                        extraBitsPerGop = CLIP3(
                            -(int64_t)(context_ptr->vbFillThreshold2 << 3),
                            (int64_t)(context_ptr->vbFillThreshold2 << 3),
                            extraBitsPerGop);
                    }
                    else
                        if ((context_ptr->extraBits > (int64_t)(context_ptr->virtualBufferSize << 3)) ||
                            (context_ptr->extraBits < -(int64_t)(context_ptr->virtualBufferSize << 3))) {
                            extraBitsPerGop = context_ptr->extraBits;
                            extraBitsPerGop = CLIP3(
                                -(int64_t)(context_ptr->vbFillThreshold2 << 2),
                                (int64_t)(context_ptr->vbFillThreshold2 << 2),
                                extraBitsPerGop);
                        }
                        else if ((context_ptr->extraBits > (int64_t)(context_ptr->virtualBufferSize << 2)) ||
                            (context_ptr->extraBits < -(int64_t)(context_ptr->virtualBufferSize << 2))) {

                            extraBitsPerGop = CLIP3(
                                -(int64_t)context_ptr->vbFillThreshold2 << 1,
                                (int64_t)context_ptr->vbFillThreshold2 << 1,
                                extraBitsPerGop);
                        }
                        else {
                            extraBitsPerGop = CLIP3(
                                -(int64_t)context_ptr->vbFillThreshold1,
                                (int64_t)context_ptr->vbFillThreshold1,
                                extraBitsPerGop);
                        }
                }
                else {
                    if ((context_ptr->extraBits > (int64_t)(context_ptr->virtualBufferSize << 3)) ||
                        (context_ptr->extraBits < -(int64_t)(context_ptr->virtualBufferSize << 3))) {
                        extraBitsPerGop = context_ptr->extraBits;
                        extraBitsPerGop = CLIP3(
                            -(int64_t)(context_ptr->vbFillThreshold2 << 2),
                            (int64_t)(context_ptr->vbFillThreshold2 << 2),
                            extraBitsPerGop);
                    }
                    else if ((context_ptr->extraBits > (int64_t)(context_ptr->virtualBufferSize << 2)) ||
                        (context_ptr->extraBits < -(int64_t)(context_ptr->virtualBufferSize << 2))) {

                        extraBitsPerGop = CLIP3(
                            -(int64_t)context_ptr->vbFillThreshold2 << 1,
                            (int64_t)context_ptr->vbFillThreshold2 << 1,
                            extraBitsPerGop);
                    }
                }

                rateControlParamPtr->virtualBufferLevel -= extraBitsPerGop;
                rateControlParamPtr->previousVirtualBufferLevel -= extraBitsPerGop;
                context_ptr->extraBits -= extraBitsPerGop;
            }

            for (temporalLayerIdex = 0; temporalLayerIdex < EB_MAX_TEMPORAL_LAYERS; temporalLayerIdex++) {

                rateControlLayerTempPtr = rateControlParamPtr->rateControlLayerArray[temporalLayerIdex];
                RateControlLayerReset(
                    rateControlLayerTempPtr,
                    picture_control_set_ptr,
                    context_ptr,
                    pictureAreaInPixel,
                    rateControlParamPtr->wasUsed);
            }
        }

        picture_control_set_ptr->parent_pcs_ptr->sad_me = 0;
        // Finding the QP of the Intra frame by using variance tables
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            uint32_t         selectedRefQp;

            if (sequence_control_set_ptr->static_config.look_ahead_distance == 0) {
                uint32_t         selectedRefQpTableIndex;
                uint32_t         intra_sad_interval_index;
                uint32_t         refQpIndex;
                uint32_t         refQpTableIndex;
                uint32_t         qpSearchMin;
                uint32_t         qpSearchMax;
                uint32_t         numOfFullLcus;
                uint64_t         minLaBitDistance;

                minLaBitDistance = MAX_UNSIGNED_VALUE;
                selectedRefQpTableIndex = 0;
                selectedRefQp = refQpListTable[selectedRefQpTableIndex];
                qpSearchMin = (uint32_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    (uint32_t)MAX((int32_t)sequence_control_set_ptr->qp - 30, 0));

                qpSearchMax = CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    sequence_control_set_ptr->qp + 30);

                if (!sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array_updated) {
                    context_ptr->intraCoefRate = CLIP3(
                        1,
                        2,
                        (uint32_t)(rateControlLayerPtr->frame_rate >> 16) / 4);
                }
                else {
                    if (context_ptr->baseLayerFramesAvgQp > context_ptr->baseLayerIntraFramesAvgQp + 3)
                        context_ptr->intraCoefRate--;
                    else if (context_ptr->baseLayerFramesAvgQp <= context_ptr->baseLayerIntraFramesAvgQp + 2)
                        context_ptr->intraCoefRate += 2;
                    else if (context_ptr->baseLayerFramesAvgQp <= context_ptr->baseLayerIntraFramesAvgQp)
                        context_ptr->intraCoefRate++;

                    context_ptr->intraCoefRate = CLIP3(
                        1,
                        15,//(uint32_t) (context_ptr->frames_in_interval[0]+1) / 2,
                        context_ptr->intraCoefRate);
                }

                // Loop over proper QPs and find the Predicted bits for that QP. Find the QP with the closest total predicted rate to target bits for the sliding window.
                for (refQpTableIndex = qpSearchMin; refQpTableIndex < qpSearchMax; refQpTableIndex++) {
                    refQpIndex = refQpListTable[refQpTableIndex];

                    picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex] = 0;

                    numOfFullLcus = 0;
                    // Loop over block in the frame and calculated the predicted bits at reg QP
                    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

                        sb_ptr = picture_control_set_ptr->sb_ptr_array[lcuCodingOrder];
                        sb_width = (sequence_control_set_ptr->luma_width - sb_ptr->origin_x >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_width - sb_ptr->origin_x;
                        sb_height = (sequence_control_set_ptr->luma_height - sb_ptr->origin_y >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_height - sb_ptr->origin_y;
                        // This is because of the tile boundry LCUs which do not have correct SAD from ME.
                        // ME doesn't know about Tile Boundries
                        if ((sb_width == BLOCK_SIZE_64) && (sb_height == BLOCK_SIZE_64)) {
                            numOfFullLcus++;
                            intra_sad_interval_index = picture_control_set_ptr->parent_pcs_ptr->intra_sad_interval_index[lcuCodingOrder];
                            picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex] += sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array[refQpIndex].intraSadBitsArray[picture_control_set_ptr->temporal_layer_index][intra_sad_interval_index];

                        }
                    }

                    // Scale for in complete LCUs
                    //  pred_bits_ref_qp is normalized based on the area because of the LCUs at the tile boundries
                    picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex] = picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex] * rateControlLayerPtr->areaInPixel / (numOfFullLcus << 12);

                    if (minLaBitDistance > (uint64_t)ABS((int64_t)rateControlLayerPtr->ecBitConstraint*context_ptr->intraCoefRate - (int64_t)picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex])) {
                        minLaBitDistance = (uint64_t)ABS((int64_t)rateControlLayerPtr->ecBitConstraint*context_ptr->intraCoefRate - (int64_t)picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex]);

                        selectedRefQpTableIndex = refQpTableIndex;
                        selectedRefQp = refQpIndex;
                    }

                }

                if (!sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array_updated) {
                    picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)selectedRefQp - (int32_t)1, 0);
                    rateControlLayerPtr->calculatedFrameQp = (uint8_t)MAX((int32_t)selectedRefQp - (int32_t)1, 0);
                    picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
                }
                else {
                    picture_control_set_ptr->picture_qp = (uint8_t)selectedRefQp;
                    rateControlLayerPtr->calculatedFrameQp = (uint8_t)selectedRefQp;
                    picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                        (uint32_t)MAX((int32_t)context_ptr->baseLayerFramesAvgQp - (int32_t)3, 0),
                        context_ptr->baseLayerFramesAvgQp,
                        picture_control_set_ptr->picture_qp);
                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                        (uint32_t)MAX((int32_t)context_ptr->baseLayerIntraFramesAvgQp - (int32_t)5, 0),
                        context_ptr->baseLayerIntraFramesAvgQp + 2,
                        picture_control_set_ptr->picture_qp);
                }
            }
            else {
                selectedRefQp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                picture_control_set_ptr->picture_qp = (uint8_t)selectedRefQp;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
                if (rateControlParamPtr->firstPoc == 0) {
                    picture_control_set_ptr->picture_qp++;
                }
            }

            // Update the QP based on the VB
            if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                if (rateControlParamPtr->virtualBufferLevel >= context_ptr->vbFillThreshold2 << 1) {
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 2;
                }
                else if (rateControlParamPtr->virtualBufferLevel >= context_ptr->vbFillThreshold2) {
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE;
                }
                else if (rateControlParamPtr->virtualBufferLevel >= context_ptr->vbFillThreshold1 &&
                    rateControlParamPtr->virtualBufferLevel < context_ptr->vbFillThreshold2) {
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD1QPINCREASE;
                }
                if (rateControlParamPtr->virtualBufferLevel <= -(context_ptr->vbFillThreshold2 << 2))
                    picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - (int32_t)2, 0);
                else
                    if (rateControlParamPtr->virtualBufferLevel <= -(context_ptr->vbFillThreshold2 << 1))
                        picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - (int32_t)1, 0);
                    else if (rateControlParamPtr->virtualBufferLevel <= 0)
                        picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
            }
            else {

                if (rateControlParamPtr->virtualBufferLevel >= context_ptr->vbFillThreshold2) {
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE;
                }
                if (rateControlParamPtr->virtualBufferLevel <= -(context_ptr->vbFillThreshold2 << 2))
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp - (uint8_t)THRESHOLD2QPINCREASE - (int32_t)2;
                else if (rateControlParamPtr->virtualBufferLevel <= -(context_ptr->vbFillThreshold2 << 1))
                    picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp - (uint8_t)THRESHOLD2QPINCREASE - (int32_t)1;
                else
                    if (rateControlParamPtr->virtualBufferLevel <= 0)
                        picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
            }
            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                picture_control_set_ptr->picture_qp);
        }
        else {

            // SB Loop
            for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

                sb_ptr = picture_control_set_ptr->sb_ptr_array[lcuCodingOrder];
                sb_width = (sequence_control_set_ptr->luma_width - sb_ptr->origin_x >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_width - sb_ptr->origin_x;
                sb_height = (sequence_control_set_ptr->luma_height - sb_ptr->origin_y >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_height - sb_ptr->origin_y;
                // This is because of the tile boundry LCUs which do not have correct SAD from ME.
                // ME doesn't know about Tile Boundries
                if ((sb_width == BLOCK_SIZE_64) && (sb_height == BLOCK_SIZE_64)) {
                    picture_control_set_ptr->parent_pcs_ptr->sad_me += picture_control_set_ptr->parent_pcs_ptr->rc_me_distortion[lcuCodingOrder];
                }
            }

            //  tileSadMe is normalized based on the area because of the LCUs at the tile boundries
            picture_control_set_ptr->parent_pcs_ptr->sad_me = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me*rateControlLayerPtr->areaInPixel / (areaInLcus << 12)), 1);

            // totalSquareMad has RC_PRECISION precision
            picture_control_set_ptr->parent_pcs_ptr->sad_me <<= RC_PRECISION;

        }

        tempQp = picture_control_set_ptr->picture_qp;

        if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc) {
            uint32_t temporalLayerIdex;
            for (temporalLayerIdex = 0; temporalLayerIdex < EB_MAX_TEMPORAL_LAYERS; temporalLayerIdex++) {
                rateControlLayerTempPtr = rateControlParamPtr->rateControlLayerArray[temporalLayerIdex];
                RateControlLayerResetPart2(
                    rateControlLayerTempPtr,
                    picture_control_set_ptr);
            }
        }

        if (picture_control_set_ptr->picture_number == 0) {
            context_ptr->baseLayerFramesAvgQp = picture_control_set_ptr->picture_qp + 1;
            context_ptr->baseLayerIntraFramesAvgQp = picture_control_set_ptr->picture_qp;
        }
    }
    else {

        picture_control_set_ptr->parent_pcs_ptr->sad_me = 0;

        // if the pixture is an I slice, for now we set the QP as the QP of the previous frame
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            uint32_t         selectedRefQp;

            if (sequence_control_set_ptr->static_config.look_ahead_distance == 0)
            {
                uint32_t         selectedRefQpTableIndex;
                uint32_t         intra_sad_interval_index;
                uint32_t         refQpIndex;
                uint32_t         refQpTableIndex;
                uint32_t         qpSearchMin;
                uint32_t         qpSearchMax;
                uint32_t         numOfFullLcus;
                uint64_t         minLaBitDistance;

                minLaBitDistance = MAX_UNSIGNED_VALUE;
                selectedRefQpTableIndex = 0;
                selectedRefQp = refQpListTable[selectedRefQpTableIndex];
                qpSearchMin = (uint8_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    (uint32_t)MAX((int32_t)sequence_control_set_ptr->qp - 20, 0));

                qpSearchMax = (uint8_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    sequence_control_set_ptr->qp + 20);

                context_ptr->intraCoefRate = CLIP3(
                    1,
                    (uint32_t)(rateControlLayerPtr->frame_rate >> 16) / 4,
                    context_ptr->intraCoefRate);
                // Loop over proper QPs and find the Predicted bits for that QP. Find the QP with the closest total predicted rate to target bits for the sliding window.
                for (refQpTableIndex = qpSearchMin; refQpTableIndex < qpSearchMax; refQpTableIndex++) {
                    refQpIndex = refQpListTable[refQpTableIndex];
                    picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex] = 0;
                    numOfFullLcus = 0;
                    // Loop over block in the frame and calculated the predicted bits at reg QP
                    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

                        sb_ptr = picture_control_set_ptr->sb_ptr_array[lcuCodingOrder];
                        sb_width = (sequence_control_set_ptr->luma_width - sb_ptr->origin_x >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_width - sb_ptr->origin_x;
                        sb_height = (sequence_control_set_ptr->luma_height - sb_ptr->origin_y >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_height - sb_ptr->origin_y;
                        // This is because of the tile boundry LCUs which do not have correct SAD from ME.
                        // ME doesn't know about Tile Boundries
                        if ((sb_width == BLOCK_SIZE_64) && (sb_height == BLOCK_SIZE_64)) {
                            numOfFullLcus++;
                            intra_sad_interval_index = picture_control_set_ptr->parent_pcs_ptr->intra_sad_interval_index[lcuCodingOrder];
                            picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex] += sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array[refQpIndex].intraSadBitsArray[picture_control_set_ptr->temporal_layer_index][intra_sad_interval_index];
                        }
                    }

                    // Scale for in complete LCUs
                    // pred_bits_ref_qp is normalized based on the area because of the LCUs at the tile boundries
                    picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex] = picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex] * rateControlLayerPtr->areaInPixel / (numOfFullLcus << 12);
                    if (minLaBitDistance > (uint64_t)ABS((int64_t)rateControlLayerPtr->ecBitConstraint*context_ptr->intraCoefRate - (int64_t)picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex])) {
                        minLaBitDistance = (uint64_t)ABS((int64_t)rateControlLayerPtr->ecBitConstraint*context_ptr->intraCoefRate - (int64_t)picture_control_set_ptr->parent_pcs_ptr->pred_bits_ref_qp[refQpIndex]);

                        selectedRefQpTableIndex = refQpTableIndex;
                        selectedRefQp = refQpIndex;
                    }
                }
                if (!sequence_control_set_ptr->encode_context_ptr->rate_control_tables_array_updated) {
                    picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)selectedRefQp - (int32_t)1, 0);
                    rateControlLayerPtr->calculatedFrameQp = (uint8_t)MAX((int32_t)selectedRefQp - (int32_t)1, 0);
                    picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
                }
                else {
                    picture_control_set_ptr->picture_qp = (uint8_t)selectedRefQp;
                    rateControlLayerPtr->calculatedFrameQp = (uint8_t)selectedRefQp;
                    picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                        (uint32_t)MAX((int32_t)context_ptr->baseLayerFramesAvgQp - (int32_t)3, 0),
                        context_ptr->baseLayerFramesAvgQp + 1,
                        picture_control_set_ptr->picture_qp);
                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                        (uint32_t)MAX((int32_t)context_ptr->baseLayerIntraFramesAvgQp - (int32_t)5, 0),
                        context_ptr->baseLayerIntraFramesAvgQp + 2,
                        picture_control_set_ptr->picture_qp);
                }
            }
            else {
                selectedRefQp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                picture_control_set_ptr->picture_qp = (uint8_t)selectedRefQp;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->picture_qp;
            }

            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                sequence_control_set_ptr->static_config.min_qp_allowed,
                sequence_control_set_ptr->static_config.max_qp_allowed,
                picture_control_set_ptr->picture_qp);

            tempQp = picture_control_set_ptr->picture_qp;

        }

        else { // Not an I slice
            // combining the target rate from initial RC and frame level RC
            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = rateControlLayerPtr->bitConstraint;
                rateControlLayerPtr->ecBitConstraint = (rateControlLayerPtr->alpha * picture_control_set_ptr->parent_pcs_ptr->target_bits_best_pred_qp +
                    ((1 << RC_PRECISION) - rateControlLayerPtr->alpha) * picture_control_set_ptr->parent_pcs_ptr->target_bits_rc + RC_PRECISION_OFFSET) >> RC_PRECISION;

                rateControlLayerPtr->ecBitConstraint = (uint64_t)MAX((int64_t)rateControlLayerPtr->ecBitConstraint - (int64_t)rateControlLayerPtr->difTotalAndEcBits, 1);

                picture_control_set_ptr->parent_pcs_ptr->target_bits_rc = rateControlLayerPtr->ecBitConstraint;
            }

            // SB Loop
            for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {

                sb_ptr = picture_control_set_ptr->sb_ptr_array[lcuCodingOrder];
                sb_width = (sequence_control_set_ptr->luma_width - sb_ptr->origin_x >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_width - sb_ptr->origin_x;
                sb_height = (sequence_control_set_ptr->luma_height - sb_ptr->origin_y >= (uint16_t)BLOCK_SIZE_64) ? sb_ptr->size : sequence_control_set_ptr->luma_height - sb_ptr->origin_y;
                // This is because of the tile boundry LCUs which do not have correct SAD from ME.
                // ME doesn't know about Tile Boundries
                if ((sb_width == BLOCK_SIZE_64) && (sb_height == BLOCK_SIZE_64)) {
                    picture_control_set_ptr->parent_pcs_ptr->sad_me += picture_control_set_ptr->parent_pcs_ptr->rc_me_distortion[lcuCodingOrder];
                }
            }

            //  tileSadMe is normalized based on the area because of the LCUs at the tile boundries
            picture_control_set_ptr->parent_pcs_ptr->sad_me = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me*rateControlLayerPtr->areaInPixel / (areaInLcus << 12)), 1);
            picture_control_set_ptr->parent_pcs_ptr->sad_me <<= RC_PRECISION;

            rateControlLayerPtr->totalMad = MAX((picture_control_set_ptr->parent_pcs_ptr->sad_me / rateControlLayerPtr->areaInPixel), 1);

            if (!rateControlLayerPtr->feedbackArrived) {
                rateControlLayerPtr->previousFrameSadMe = picture_control_set_ptr->parent_pcs_ptr->sad_me;
            }

            {
                uint64_t qpCalcTemp1, qpCalcTemp2, qpCalcTemp3;

                qpCalcTemp1 = picture_control_set_ptr->parent_pcs_ptr->sad_me *rateControlLayerPtr->totalMad;
                qpCalcTemp2 =
                    MAX((int64_t)(rateControlLayerPtr->ecBitConstraint << (2 * RC_PRECISION)) - (int64_t)rateControlLayerPtr->cCoeff*(int64_t)rateControlLayerPtr->areaInPixel,
                    (int64_t)(rateControlLayerPtr->ecBitConstraint << (2 * RC_PRECISION - 2)));

                // This is a more complex but with higher precision implementation
                if (qpCalcTemp1 > qpCalcTemp2)
                    qpCalcTemp3 = (uint64_t)((qpCalcTemp1 / qpCalcTemp2)*rateControlLayerPtr->kCoeff);
                else
                    qpCalcTemp3 = (uint64_t)(qpCalcTemp1*rateControlLayerPtr->kCoeff / qpCalcTemp2);
                tempQp = (uint64_t)(Log2fHighPrecision(MAX(((qpCalcTemp3 + RC_PRECISION_OFFSET) >> RC_PRECISION)*((qpCalcTemp3 + RC_PRECISION_OFFSET) >> RC_PRECISION)*((qpCalcTemp3 + RC_PRECISION_OFFSET) >> RC_PRECISION), 1), RC_PRECISION));

                rateControlLayerPtr->calculatedFrameQp = (uint8_t)(CLIP3(1, 63, (uint32_t)(tempQp + RC_PRECISION_OFFSET) >> RC_PRECISION));
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp = (uint8_t)(CLIP3(1, 63, (uint32_t)(tempQp + RC_PRECISION_OFFSET) >> RC_PRECISION));
            }

            tempQp += rateControlLayerPtr->deltaQpFraction;
            picture_control_set_ptr->picture_qp = (uint8_t)((tempQp + RC_PRECISION_OFFSET) >> RC_PRECISION);
            // Use the QP of HLRC instead of calculated one in FLRC
            picture_control_set_ptr->picture_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
            picture_control_set_ptr->parent_pcs_ptr->calculated_qp = picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
            if (rateControlParamPtr->firstPoc == 0) {
                picture_control_set_ptr->parent_pcs_ptr->best_pred_qp++;
                picture_control_set_ptr->picture_qp++;
                picture_control_set_ptr->parent_pcs_ptr->calculated_qp++;
            }
        }
        if (picture_control_set_ptr->parent_pcs_ptr->first_non_intra_frame_in_temporal_layer && picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {
            picture_control_set_ptr->picture_qp = (uint8_t)rateControlParamPtr->intraFramesQp + 1;
        }

        if (!rateControlLayerPtr->feedbackArrived && picture_control_set_ptr->slice_type != I_SLICE) {

            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                (int32_t)sequence_control_set_ptr->static_config.min_qp_allowed,
                (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed,
                (int32_t)(rateControlParamPtr->intraFramesQp + QP_OFFSET_LAYER_ARRAY[picture_control_set_ptr->temporal_layer_index] - RC_INTRA_QP_OFFSET));

        }

        if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
            if (rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold2 << 2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 4;
            }
            else if (rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold2 << 1) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 3;
            }
            else if (rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 2;
            }
            else if (rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold1 &&
                rateControlParamPtr->virtualBufferLevel < context_ptr->vbFillThreshold2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD1QPINCREASE + 2;
            }
        }
        else {


            if (rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold2 << 2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 2;

            }
            else if (rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold2 << 1) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 1;


            }
            else if (rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD2QPINCREASE + 1;
            }
            else if (rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold1 &&
                rateControlParamPtr->virtualBufferLevel < context_ptr->vbFillThreshold2) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + (uint8_t)THRESHOLD1QPINCREASE;
            }

        }
        if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
            if (rateControlParamPtr->virtualBufferLevel < -(context_ptr->vbFillThreshold2 << 2))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 2, 0);
            else if (rateControlParamPtr->virtualBufferLevel < -(context_ptr->vbFillThreshold2 << 1))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 1, 0);
            else if (rateControlParamPtr->virtualBufferLevel < 0)
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
        }
        else {

            if (rateControlParamPtr->virtualBufferLevel < -(context_ptr->vbFillThreshold2 << 2))
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE - 1, 0);
            else if (rateControlParamPtr->virtualBufferLevel < -context_ptr->vbFillThreshold2)
                picture_control_set_ptr->picture_qp = (uint8_t)MAX((int32_t)picture_control_set_ptr->picture_qp - (int32_t)THRESHOLD2QPINCREASE, 0);
        }
        // limiting the QP based on the predicted QP
        if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
            if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)MAX((int32_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp - 8, 0),
                    (uint32_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp + 8,
                    (uint32_t)picture_control_set_ptr->picture_qp);
            }
            else {
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    (uint32_t)MAX((int32_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp - 8, 0),
                    (uint32_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp + 8,
                    (uint32_t)picture_control_set_ptr->picture_qp);

            }
        }
        if (picture_control_set_ptr->picture_number != rateControlParamPtr->firstPoc &&
            picture_control_set_ptr->picture_qp == picture_control_set_ptr->parent_pcs_ptr->best_pred_qp && rateControlParamPtr->virtualBufferLevel > context_ptr->vbFillThreshold1) {
            if (rateControlParamPtr->extraApBitRatioI > 200) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 3;
            }
            else if (rateControlParamPtr->extraApBitRatioI > 100) {
                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 2;
            }
            else if (rateControlParamPtr->extraApBitRatioI > 50) {
                picture_control_set_ptr->picture_qp++;
            }
        }
        //Limiting the QP based on the QP of the Reference frame
        {
            uint32_t ref_qp;

            if (!picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                if (context_ptr->frames_in_interval[picture_control_set_ptr->temporal_layer_index] < 5) {
                    if ((int32_t)picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {
                        if (picture_control_set_ptr->ref_slice_type_array[0] == I_SLICE) {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0],
                                (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0] + 2,
                                picture_control_set_ptr->picture_qp);
                        }
                        else {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)picture_control_set_ptr->ref_pic_qp_array[0] - 1, 0),
                                (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0] + 1,
                                picture_control_set_ptr->picture_qp);

                        }
                    }
                    if ((int32_t)picture_control_set_ptr->temporal_layer_index == 1) {
                        ref_qp = 0;
                        if (picture_control_set_ptr->ref_slice_type_array[0] != I_SLICE) {
                            ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[0]);
                        }

                        if ((picture_control_set_ptr->slice_type == B_SLICE) && (picture_control_set_ptr->ref_slice_type_array[1] != I_SLICE)) {
                            ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[1]);
                        }
                        if (ref_qp > 0) {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)ref_qp - 1,
                                (uint32_t)ref_qp + 1,
                                picture_control_set_ptr->picture_qp);
                        }
                    }
                }
                else {
                    if ((int32_t)picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {
                        if (picture_control_set_ptr->ref_slice_type_array[0] == I_SLICE) {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0],
                                (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0] + 2,
                                picture_control_set_ptr->picture_qp);
                        }
                        else {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)picture_control_set_ptr->ref_pic_qp_array[0] - 1, 0),
                                (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0] + 1,
                                picture_control_set_ptr->picture_qp);
                        }
                    }
                    if ((int32_t)picture_control_set_ptr->temporal_layer_index == 1) {
                        ref_qp = 0;
                        if (picture_control_set_ptr->ref_slice_type_array[0] != I_SLICE) {
                            ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[0]);
                        }
                        if ((picture_control_set_ptr->slice_type == B_SLICE) && (picture_control_set_ptr->ref_slice_type_array[1] != I_SLICE)) {
                            ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[1]);
                        }
                        if (ref_qp > 0) {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)ref_qp - 1,
                                (uint32_t)ref_qp + 1,
                                picture_control_set_ptr->picture_qp);
                        }
                    }
                }

                if ((int32_t)picture_control_set_ptr->temporal_layer_index == 2) {
                    ref_qp = 0;
                    if (picture_control_set_ptr->ref_slice_type_array[0] != I_SLICE) {
                        ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[0]);
                    }

                    if ((picture_control_set_ptr->slice_type == B_SLICE) && (picture_control_set_ptr->ref_slice_type_array[1] != I_SLICE)) {
                        ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[1]);
                    }
                    if (ref_qp > 0) {
                        if (picture_control_set_ptr->slice_type == P_SLICE) {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)ref_qp - 1, 0),
                                (uint32_t)ref_qp + 1,
                                picture_control_set_ptr->picture_qp);
                        }
                        else {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)ref_qp - 1, 0),
                                (uint32_t)ref_qp + 2,
                                picture_control_set_ptr->picture_qp);
                        }
                    }
                }

                if ((int32_t)picture_control_set_ptr->temporal_layer_index >= 3) {

                    ref_qp = 0;
                    if (picture_control_set_ptr->ref_slice_type_array[0] != I_SLICE) {
                        ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[0]);
                    }
                    if ((picture_control_set_ptr->slice_type == B_SLICE) && (picture_control_set_ptr->ref_slice_type_array[1] != I_SLICE)) {
                        ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[1]);
                    }
                    if (ref_qp > 0) {
                        if (picture_control_set_ptr->slice_type == P_SLICE) {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)ref_qp - 1, 0),
                                (uint32_t)ref_qp + 2,
                                picture_control_set_ptr->picture_qp);
                        }
                        else {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)ref_qp - 1, 0),
                                (uint32_t)ref_qp + 3,
                                picture_control_set_ptr->picture_qp);
                        }
                    }
                }
            }
            else {
                if ((int32_t)picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {
                    if (picture_control_set_ptr->ref_slice_type_array[0] == I_SLICE) {
                        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                            (uint32_t)MAX((int32_t)picture_control_set_ptr->ref_pic_qp_array[0] - 2, 0),
                            (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0] + 3,
                            picture_control_set_ptr->picture_qp);
                    }
                    else {
                        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                            (uint32_t)MAX((int32_t)picture_control_set_ptr->ref_pic_qp_array[0] - 2, 0),
                            (uint32_t)picture_control_set_ptr->ref_pic_qp_array[0] + 3,
                            picture_control_set_ptr->picture_qp);
                    }
                }

                if ((int32_t)picture_control_set_ptr->temporal_layer_index >= 1) {
                    ref_qp = 0;
                    if (picture_control_set_ptr->ref_slice_type_array[0] != I_SLICE) {
                        ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[0]);
                    }
                    if ((picture_control_set_ptr->slice_type == B_SLICE) && (picture_control_set_ptr->ref_slice_type_array[1] != I_SLICE)) {
                        ref_qp = MAX(ref_qp, picture_control_set_ptr->ref_pic_qp_array[1]);
                    }
                    if (ref_qp > 0) {
                        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                            (uint32_t)MAX((int32_t)ref_qp - 2, 0),
                            (uint32_t)ref_qp + 3,
                            picture_control_set_ptr->picture_qp);
                    }
                }
            }
        }
        // limiting the QP between min Qp allowed and max Qp allowed
        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
            sequence_control_set_ptr->static_config.min_qp_allowed,
            sequence_control_set_ptr->static_config.max_qp_allowed,
            picture_control_set_ptr->picture_qp);

        rateControlLayerPtr->deltaQpFraction = CLIP3(-RC_PRECISION_OFFSET, RC_PRECISION_OFFSET, -((int64_t)tempQp - (int64_t)(picture_control_set_ptr->picture_qp << RC_PRECISION)));

        if (picture_control_set_ptr->parent_pcs_ptr->sad_me == rateControlLayerPtr->previousFrameSadMe &&
            (rateControlLayerPtr->previousFrameSadMe != 0))
            rateControlLayerPtr->sameSADCount++;
        else
            rateControlLayerPtr->sameSADCount = 0;
    }

    rateControlLayerPtr->previousCCoeff = rateControlLayerPtr->cCoeff;
    rateControlLayerPtr->previousKCoeff = rateControlLayerPtr->kCoeff;
    rateControlLayerPtr->previousCalculatedFrameQp = rateControlLayerPtr->calculatedFrameQp;
}
void FrameLevelRcFeedbackPictureMode2(
    PictureParentControlSet_t         *parentPictureControlSetPtr,
    SequenceControlSet_t              *sequence_control_set_ptr,
    RateControlContext_t              *context_ptr)
{



    RateControlLayerContext_t           *rateControlLayerTempPtr;
    RateControlIntervalParamContext_t   *rateControlParamPtr;
    RateControlLayerContext_t           *rateControlLayerPtr;
    // SB Loop variables
    uint32_t                       sliceNum;
    uint64_t                       previousFrameBitActual;

    if (sequence_control_set_ptr->intra_period_length == -1)
        rateControlParamPtr = context_ptr->rateControlParamQueue[0];
    else {
        uint32_t intervalIndexTemp = 0;
        while ((!(parentPictureControlSetPtr->picture_number >= context_ptr->rateControlParamQueue[intervalIndexTemp]->firstPoc &&
            parentPictureControlSetPtr->picture_number <= context_ptr->rateControlParamQueue[intervalIndexTemp]->lastPoc)) &&
            (intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER)) {
            intervalIndexTemp++;
        }
        CHECK_REPORT_ERROR(
            intervalIndexTemp != PARALLEL_GOP_MAX_NUMBER,
            sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
            EB_ENC_RC_ERROR2);
        rateControlParamPtr = context_ptr->rateControlParamQueue[intervalIndexTemp];
    }

    rateControlLayerPtr = rateControlParamPtr->rateControlLayerArray[parentPictureControlSetPtr->temporal_layer_index];

    rateControlLayerPtr->maxQp = 0;

    rateControlLayerPtr->feedbackArrived = EB_TRUE;
    rateControlLayerPtr->maxQp = MAX(rateControlLayerPtr->maxQp, parentPictureControlSetPtr->picture_qp);

    rateControlLayerPtr->previousFrameQp = parentPictureControlSetPtr->picture_qp;
    rateControlLayerPtr->previousFrameBitActual = parentPictureControlSetPtr->total_num_bits;
    if (parentPictureControlSetPtr->quantized_coeff_num_bits == 0)
        parentPictureControlSetPtr->quantized_coeff_num_bits = 1;
    rateControlLayerPtr->previousFrameQuantizedCoeffBitActual = parentPictureControlSetPtr->quantized_coeff_num_bits;

    // Setting Critical states for adjusting the averaging weights on C and K
    if ((parentPictureControlSetPtr->sad_me > (3 * rateControlLayerPtr->previousFrameSadMe) >> 1) &&
        (rateControlLayerPtr->previousFrameSadMe != 0)) {
        rateControlLayerPtr->criticalStates = 3;
    }
    else if (rateControlLayerPtr->criticalStates) {
        rateControlLayerPtr->criticalStates--;
    }
    else {
        rateControlLayerPtr->criticalStates = 0;
    }

    if (parentPictureControlSetPtr->slice_type != I_SLICE) {
        // Updating CCoeff
        rateControlLayerPtr->cCoeff = (((int64_t)rateControlLayerPtr->previousFrameBitActual - (int64_t)rateControlLayerPtr->previousFrameQuantizedCoeffBitActual) << (2 * RC_PRECISION))
            / rateControlLayerPtr->areaInPixel;
        rateControlLayerPtr->cCoeff = MAX(rateControlLayerPtr->cCoeff, 1);

        // Updating KCoeff
        if ((parentPictureControlSetPtr->sad_me + RC_PRECISION_OFFSET) >> RC_PRECISION != 0) {
            {
                uint64_t test1, test2, test3;
                test1 = rateControlLayerPtr->previousFrameQuantizedCoeffBitActual*(TWO_TO_POWER_QP_OVER_THREE[parentPictureControlSetPtr->picture_qp]);
                test2 = MAX(parentPictureControlSetPtr->sad_me / rateControlLayerPtr->areaInPixel, 1);
                test3 = test1 * 65536 / test2 * 65536 / parentPictureControlSetPtr->sad_me;

                rateControlLayerPtr->kCoeff = test3;
            }
        }

        if (rateControlLayerPtr->criticalStates) {
            rateControlLayerPtr->kCoeff = (8 * rateControlLayerPtr->kCoeff + 8 * rateControlLayerPtr->previousKCoeff + 8) >> 4;
            rateControlLayerPtr->cCoeff = (8 * rateControlLayerPtr->cCoeff + 8 * rateControlLayerPtr->previousCCoeff + 8) >> 4;
        }
        else {
            rateControlLayerPtr->kCoeff = (rateControlLayerPtr->coeffAveragingWeight1*rateControlLayerPtr->kCoeff + rateControlLayerPtr->coeffAveragingWeight2*rateControlLayerPtr->previousKCoeff + 8) >> 4;
            rateControlLayerPtr->cCoeff = (rateControlLayerPtr->coeffAveragingWeight1*rateControlLayerPtr->cCoeff + rateControlLayerPtr->coeffAveragingWeight2*rateControlLayerPtr->previousCCoeff + 8) >> 4;
        }

        if (parentPictureControlSetPtr->slice_type != I_SLICE) {
            rateControlLayerPtr->previousFrameSadMe = parentPictureControlSetPtr->sad_me;
        }
        else {
            rateControlLayerPtr->previousFrameSadMe = 0;
        }
    }

    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
        if (parentPictureControlSetPtr->slice_type == I_SLICE) {
            if (parentPictureControlSetPtr->total_num_bits < parentPictureControlSetPtr->target_bits_best_pred_qp << 1)
                context_ptr->baseLayerIntraFramesAvgQp = (3 * context_ptr->baseLayerIntraFramesAvgQp + parentPictureControlSetPtr->picture_qp + 2) >> 2;
            else if (parentPictureControlSetPtr->total_num_bits > parentPictureControlSetPtr->target_bits_best_pred_qp << 2)
                context_ptr->baseLayerIntraFramesAvgQp = (3 * context_ptr->baseLayerIntraFramesAvgQp + parentPictureControlSetPtr->picture_qp + 4 + 2) >> 2;
            else if (parentPictureControlSetPtr->total_num_bits > parentPictureControlSetPtr->target_bits_best_pred_qp << 1)
                context_ptr->baseLayerIntraFramesAvgQp = (3 * context_ptr->baseLayerIntraFramesAvgQp + parentPictureControlSetPtr->picture_qp + 2 + 2) >> 2;
        }
    }


    {
        uint64_t previousFrameEcBits = 0;
        EbBool pictureMinQpAllowed = EB_TRUE;
        rateControlLayerPtr->previousFrameAverageQp = 0;
        rateControlLayerPtr->previousFrameAverageQp += rateControlLayerPtr->previousFrameQp;
        previousFrameEcBits += rateControlLayerPtr->previousFrameBitActual;
        if (rateControlLayerPtr->sameSADCount == 0 ||
            parentPictureControlSetPtr->picture_qp != sequence_control_set_ptr->static_config.min_qp_allowed) {
            pictureMinQpAllowed = EB_FALSE;
        }
        if (pictureMinQpAllowed)
            rateControlLayerPtr->frameSameSADMinQpCount++;
        else
            rateControlLayerPtr->frameSameSADMinQpCount = 0;

        rateControlLayerPtr->previousEcBits = previousFrameEcBits;
        previousFrameBitActual = parentPictureControlSetPtr->total_num_bits;
        if (parentPictureControlSetPtr->first_frame_in_temporal_layer) {
            rateControlLayerPtr->difTotalAndEcBits = (previousFrameBitActual - previousFrameEcBits);
        }
        else {
            rateControlLayerPtr->difTotalAndEcBits = ((previousFrameBitActual - previousFrameEcBits) + rateControlLayerPtr->difTotalAndEcBits) >> 1;
        }

        // update bitrate of different layers in the interval based on the rate of the I frame
        if (parentPictureControlSetPtr->picture_number == rateControlParamPtr->firstPoc &&
            (parentPictureControlSetPtr->slice_type == I_SLICE) &&
            sequence_control_set_ptr->static_config.intra_period_length != -1) {
            uint32_t temporalLayerIdex;
            uint64_t target_bit_rate;
            uint64_t channelBitRate;
            uint64_t sumBitsPerSw = 0;
#if ADAPTIVE_PERCENTAGE
            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                if (parentPictureControlSetPtr->tables_updated && parentPictureControlSetPtr->percentage_updated) {
                    parentPictureControlSetPtr->bits_per_sw_per_layer[0] =
                        (uint64_t)MAX((int64_t)parentPictureControlSetPtr->bits_per_sw_per_layer[0] + (int64_t)parentPictureControlSetPtr->total_num_bits - (int64_t)parentPictureControlSetPtr->target_bits_best_pred_qp, 1);
                }
            }
#endif

            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
                for (temporalLayerIdex = 0; temporalLayerIdex < EB_MAX_TEMPORAL_LAYERS; temporalLayerIdex++) {
                    sumBitsPerSw += parentPictureControlSetPtr->bits_per_sw_per_layer[temporalLayerIdex];
                }
            }

            for (temporalLayerIdex = 0; temporalLayerIdex < EB_MAX_TEMPORAL_LAYERS; temporalLayerIdex++) {
                rateControlLayerTempPtr = rateControlParamPtr->rateControlLayerArray[temporalLayerIdex];

                target_bit_rate = (uint64_t)((int64_t)parentPictureControlSetPtr->target_bit_rate -
                    MIN((int64_t)parentPictureControlSetPtr->target_bit_rate * 3 / 4, (int64_t)(parentPictureControlSetPtr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION))
                    *RATE_PERCENTAGE_LAYER_ARRAY[sequence_control_set_ptr->static_config.hierarchical_levels][temporalLayerIdex] / 100;

#if ADAPTIVE_PERCENTAGE
                if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->intra_period_length != -1) {
                    target_bit_rate = (uint64_t)((int64_t)parentPictureControlSetPtr->target_bit_rate -
                        MIN((int64_t)parentPictureControlSetPtr->target_bit_rate * 3 / 4, (int64_t)(parentPictureControlSetPtr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION))
                        *parentPictureControlSetPtr->bits_per_sw_per_layer[temporalLayerIdex] / sumBitsPerSw;
                }
#endif
                // update this based on temporal layers
                if (temporalLayerIdex == 0)
                    channelBitRate = (((target_bit_rate << (2 * RC_PRECISION)) / MAX(1, rateControlLayerTempPtr->frame_rate - (1 * context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)))) + RC_PRECISION_OFFSET) >> RC_PRECISION;
                else
                    channelBitRate = (((target_bit_rate << (2 * RC_PRECISION)) / rateControlLayerTempPtr->frame_rate) + RC_PRECISION_OFFSET) >> RC_PRECISION;
                channelBitRate = (uint64_t)MAX((int64_t)1, (int64_t)channelBitRate);
                rateControlLayerTempPtr->ecBitConstraint = channelBitRate;

                sliceNum = 1;
                rateControlLayerTempPtr->ecBitConstraint -= SLICE_HEADER_BITS_NUM * sliceNum;

                rateControlLayerTempPtr->previousBitConstraint = channelBitRate;
                rateControlLayerTempPtr->bitConstraint = channelBitRate;
                rateControlLayerTempPtr->channelBitRate = channelBitRate;
            }
            if ((int64_t)parentPictureControlSetPtr->target_bit_rate * 3 / 4 < (int64_t)(parentPictureControlSetPtr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION) {
                rateControlParamPtr->previousVirtualBufferLevel += (int64_t)((parentPictureControlSetPtr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION) - (int64_t)parentPictureControlSetPtr->target_bit_rate * 3 / 4;
                context_ptr->extraBitsGen -= (int64_t)((parentPictureControlSetPtr->total_num_bits*context_ptr->frame_rate / (sequence_control_set_ptr->static_config.intra_period_length + 1)) >> RC_PRECISION) - (int64_t)parentPictureControlSetPtr->target_bit_rate * 3 / 4;
            }
        }

        if (previousFrameBitActual) {

            uint64_t bitChangesRate;
            // Updating virtual buffer level and it can be negative
            if ((parentPictureControlSetPtr->picture_number == rateControlParamPtr->firstPoc) &&
                (parentPictureControlSetPtr->slice_type == I_SLICE) &&
                (rateControlParamPtr->lastGop == EB_FALSE) &&
                sequence_control_set_ptr->static_config.intra_period_length != -1) {
                rateControlParamPtr->virtualBufferLevel =
                    (int64_t)rateControlParamPtr->previousVirtualBufferLevel;
            }
            else {
                rateControlParamPtr->virtualBufferLevel =
                    (int64_t)rateControlParamPtr->previousVirtualBufferLevel +
                    (int64_t)previousFrameBitActual - (int64_t)rateControlLayerPtr->channelBitRate;
                context_ptr->extraBitsGen -= (int64_t)previousFrameBitActual - (int64_t)rateControlLayerPtr->channelBitRate;
            }

            if (rateControlLayerPtr->frameSameSADMinQpCount > 10) {
                rateControlLayerPtr->previousBitConstraint = (int64_t)rateControlLayerPtr->channelBitRate;
                rateControlParamPtr->virtualBufferLevel = ((int64_t)context_ptr->virtualBufferSize >> 1);
            }
            // Updating bit difference
            rateControlLayerPtr->bitDiff = (int64_t)rateControlParamPtr->virtualBufferLevel
                //- ((int64_t)context_ptr->virtualBufferSize>>1);
                - ((int64_t)rateControlLayerPtr->channelBitRate >> 1);

            // Limit the bit difference
            rateControlLayerPtr->bitDiff = CLIP3(-(int64_t)(rateControlLayerPtr->channelBitRate), (int64_t)(rateControlLayerPtr->channelBitRate >> 1), rateControlLayerPtr->bitDiff);
            bitChangesRate = rateControlLayerPtr->frame_rate;

            // Updating bit Constraint
            rateControlLayerPtr->bitConstraint = MAX((int64_t)rateControlLayerPtr->previousBitConstraint - ((rateControlLayerPtr->bitDiff << RC_PRECISION) / ((int64_t)bitChangesRate)), 1);

            // Limiting the bitConstraint
            if (parentPictureControlSetPtr->temporal_layer_index == 0) {
                rateControlLayerPtr->bitConstraint = CLIP3(rateControlLayerPtr->channelBitRate >> 2,
                    rateControlLayerPtr->channelBitRate * 200 / 100,
                    rateControlLayerPtr->bitConstraint);
            }
            else {
                rateControlLayerPtr->bitConstraint = CLIP3(rateControlLayerPtr->channelBitRate >> 1,
                    rateControlLayerPtr->channelBitRate * 200 / 100,
                    rateControlLayerPtr->bitConstraint);
            }
            rateControlLayerPtr->ecBitConstraint = (uint64_t)MAX((int64_t)rateControlLayerPtr->bitConstraint - (int64_t)rateControlLayerPtr->difTotalAndEcBits, 1);
            rateControlParamPtr->previousVirtualBufferLevel = rateControlParamPtr->virtualBufferLevel;
            rateControlLayerPtr->previousBitConstraint = rateControlLayerPtr->bitConstraint;
        }

        rateControlParamPtr->processedFramesNumber++;
        rateControlParamPtr->inUse = EB_TRUE;
        // check if all the frames in the interval have arrived
        if (rateControlParamPtr->processedFramesNumber == (rateControlParamPtr->lastPoc - rateControlParamPtr->firstPoc + 1) &&
            sequence_control_set_ptr->intra_period_length != -1) {

            uint32_t temporalIndex;
            int64_t extraBits;
            rateControlParamPtr->firstPoc += PARALLEL_GOP_MAX_NUMBER * (uint32_t)(sequence_control_set_ptr->intra_period_length + 1);
            rateControlParamPtr->lastPoc += PARALLEL_GOP_MAX_NUMBER * (uint32_t)(sequence_control_set_ptr->intra_period_length + 1);
            rateControlParamPtr->processedFramesNumber = 0;
            rateControlParamPtr->virtualBufferLevel = rateControlParamPtr->virtualBufferLevel;
            rateControlParamPtr->previousVirtualBufferLevel = rateControlParamPtr->previousVirtualBufferLevel;
            rateControlParamPtr->extraApBitRatioI = 0;
            rateControlParamPtr->inUse = EB_FALSE;
            rateControlParamPtr->wasUsed = EB_TRUE;
            rateControlParamPtr->lastGop = EB_FALSE;
            rateControlParamPtr->firstPicActualQpAssigned = EB_FALSE;
            for (temporalIndex = 0; temporalIndex < EB_MAX_TEMPORAL_LAYERS; temporalIndex++) {
                rateControlParamPtr->rateControlLayerArray[temporalIndex]->firstFrame = 1;
                rateControlParamPtr->rateControlLayerArray[temporalIndex]->firstNonIntraFrame = 1;
                rateControlParamPtr->rateControlLayerArray[temporalIndex]->feedbackArrived = EB_FALSE;
            }
            extraBits = ((int64_t)context_ptr->virtualBufferSize >> 1) - (int64_t)rateControlParamPtr->virtualBufferLevel;

            rateControlParamPtr->virtualBufferLevel = context_ptr->virtualBufferSize >> 1;
            context_ptr->extraBits += extraBits;

        }
        // Allocate the extraBits among other GOPs
        if ((parentPictureControlSetPtr->temporal_layer_index <= 2) &&
            ((context_ptr->extraBits > (int64_t)(context_ptr->virtualBufferSize >> 8)) ||
            (context_ptr->extraBits < -(int64_t)(context_ptr->virtualBufferSize >> 8)))) {
            uint32_t intervalIndexTemp, intervalInUseCount;
            int64_t extraBitsPerGop;
            int64_t extraBits = context_ptr->extraBits;
            int32_t clipCoef1, clipCoef2;
            if (parentPictureControlSetPtr->end_of_sequence_region) {
                clipCoef1 = -1;
                clipCoef2 = -1;
            }
            else {
                if (context_ptr->extraBits > (int64_t)(context_ptr->virtualBufferSize << 3) ||
                    context_ptr->extraBits < -(int64_t)(context_ptr->virtualBufferSize << 3)) {
                    clipCoef1 = 0;
                    clipCoef2 = 0;
                }
                else {
                    clipCoef1 = 2;
                    clipCoef2 = 4;
                }
            }

            intervalInUseCount = 0;

            if (extraBits > 0) {
                // Extra bits to be distributed
                // Distribute it among those that are consuming more
                for (intervalIndexTemp = 0; intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER; intervalIndexTemp++) {
                    if (context_ptr->rateControlParamQueue[intervalIndexTemp]->inUse &&
                        context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel > ((int64_t)context_ptr->virtualBufferSize >> 1)) {
                        intervalInUseCount++;
                    }
                }
                // Distribute the rate among them
                if (intervalInUseCount) {
                    extraBitsPerGop = extraBits / intervalInUseCount;
                    if (clipCoef1 > 0)
                        extraBitsPerGop = CLIP3(
                            -(int64_t)context_ptr->virtualBufferSize >> clipCoef1,
                            (int64_t)context_ptr->virtualBufferSize >> clipCoef1,
                            extraBitsPerGop);
                    else
                        extraBitsPerGop = CLIP3(
                            -(int64_t)context_ptr->virtualBufferSize << (-clipCoef1),
                            (int64_t)context_ptr->virtualBufferSize << (-clipCoef1),
                            extraBitsPerGop);

                    for (intervalIndexTemp = 0; intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER; intervalIndexTemp++) {
                        if (context_ptr->rateControlParamQueue[intervalIndexTemp]->inUse &&
                            context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel > ((int64_t)context_ptr->virtualBufferSize >> 1)) {
                            context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel -= extraBitsPerGop;
                            context_ptr->rateControlParamQueue[intervalIndexTemp]->previousVirtualBufferLevel -= extraBitsPerGop;
                            context_ptr->extraBits -= extraBitsPerGop;
                        }
                    }
                }
                // if no interval with more consuming was found, allocate it to ones with consuming less
                else {
                    intervalInUseCount = 0;
                    // Distribute it among those that are consuming less
                    for (intervalIndexTemp = 0; intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER; intervalIndexTemp++) {

                        if (context_ptr->rateControlParamQueue[intervalIndexTemp]->inUse &&
                            context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel <= ((int64_t)context_ptr->virtualBufferSize >> 1)) {
                            intervalInUseCount++;
                        }
                    }
                    if (intervalInUseCount) {
                        extraBitsPerGop = extraBits / intervalInUseCount;
                        if (clipCoef2 > 0)
                            extraBitsPerGop = CLIP3(
                                -(int64_t)context_ptr->virtualBufferSize >> clipCoef2,
                                (int64_t)context_ptr->virtualBufferSize >> clipCoef2,
                                extraBitsPerGop);
                        else
                            extraBitsPerGop = CLIP3(
                                -(int64_t)context_ptr->virtualBufferSize << (-clipCoef2),
                                (int64_t)context_ptr->virtualBufferSize << (-clipCoef2),
                                extraBitsPerGop);
                        // Distribute the rate among them
                        for (intervalIndexTemp = 0; intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER; intervalIndexTemp++) {

                            if (context_ptr->rateControlParamQueue[intervalIndexTemp]->inUse &&
                                context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel <= ((int64_t)context_ptr->virtualBufferSize >> 1)) {
                                context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel -= extraBitsPerGop;
                                context_ptr->rateControlParamQueue[intervalIndexTemp]->previousVirtualBufferLevel -= extraBitsPerGop;
                                context_ptr->extraBits -= extraBitsPerGop;
                            }
                        }
                    }
                }
            }
            else {
                // Distribute it among those that are consuming less
                for (intervalIndexTemp = 0; intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER; intervalIndexTemp++) {

                    if (context_ptr->rateControlParamQueue[intervalIndexTemp]->inUse &&
                        context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel < ((int64_t)context_ptr->virtualBufferSize >> 1)) {
                        intervalInUseCount++;
                    }
                }
                if (intervalInUseCount) {
                    extraBitsPerGop = extraBits / intervalInUseCount;
                    if (clipCoef1 > 0)
                        extraBitsPerGop = CLIP3(
                            -(int64_t)context_ptr->virtualBufferSize >> clipCoef1,
                            (int64_t)context_ptr->virtualBufferSize >> clipCoef1,
                            extraBitsPerGop);
                    else
                        extraBitsPerGop = CLIP3(
                            -(int64_t)context_ptr->virtualBufferSize << (-clipCoef1),
                            (int64_t)context_ptr->virtualBufferSize << (-clipCoef1),
                            extraBitsPerGop);
                    // Distribute the rate among them
                    for (intervalIndexTemp = 0; intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER; intervalIndexTemp++) {
                        if (context_ptr->rateControlParamQueue[intervalIndexTemp]->inUse &&
                            context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel < ((int64_t)context_ptr->virtualBufferSize >> 1)) {
                            context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel -= extraBitsPerGop;
                            context_ptr->rateControlParamQueue[intervalIndexTemp]->previousVirtualBufferLevel -= extraBitsPerGop;
                            context_ptr->extraBits -= extraBitsPerGop;
                        }
                    }
                }
                // if no interval with less consuming was found, allocate it to ones with consuming more
                else {
                    intervalInUseCount = 0;
                    for (intervalIndexTemp = 0; intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER; intervalIndexTemp++) {

                        if (context_ptr->rateControlParamQueue[intervalIndexTemp]->inUse &&
                            context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel < (int64_t)(context_ptr->virtualBufferSize)) {
                            intervalInUseCount++;
                        }
                    }
                    if (intervalInUseCount) {
                        extraBitsPerGop = extraBits / intervalInUseCount;
                        if (clipCoef2 > 0)
                            extraBitsPerGop = CLIP3(
                                -(int64_t)context_ptr->virtualBufferSize >> clipCoef2,
                                (int64_t)context_ptr->virtualBufferSize >> clipCoef2,
                                extraBitsPerGop);
                        else
                            extraBitsPerGop = CLIP3(
                                -(int64_t)context_ptr->virtualBufferSize << (-clipCoef2),
                                (int64_t)context_ptr->virtualBufferSize << (-clipCoef2),
                                extraBitsPerGop);
                        // Distribute the rate among them
                        for (intervalIndexTemp = 0; intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER; intervalIndexTemp++) {

                            if (context_ptr->rateControlParamQueue[intervalIndexTemp]->inUse &&
                                context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel < (int64_t)(context_ptr->virtualBufferSize)) {
                                context_ptr->rateControlParamQueue[intervalIndexTemp]->virtualBufferLevel -= extraBitsPerGop;
                                context_ptr->rateControlParamQueue[intervalIndexTemp]->previousVirtualBufferLevel -= extraBitsPerGop;
                                context_ptr->extraBits -= extraBitsPerGop;
                            }
                        }
                    }
                }
            }
        }
    }
}

void HighLevelRcFeedBackPicture(
    PictureParentControlSet_t         *picture_control_set_ptr,
    SequenceControlSet_t              *sequence_control_set_ptr)
{

    // Queue variables
    HlRateControlHistogramEntry_t      *hlRateControlHistogramPtrTemp;
    uint32_t                             queueEntryIndexHeadTemp;


    //printf("\nOut:%d Slidings: ",picture_control_set_ptr->picture_number);
    if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {

        // Update the coded rate in the histogram queue
        if (picture_control_set_ptr->picture_number >= sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue[sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number) {
            queueEntryIndexHeadTemp = (int32_t)(picture_control_set_ptr->picture_number - sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue[sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
            queueEntryIndexHeadTemp += sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_head_index;
            queueEntryIndexHeadTemp = (queueEntryIndexHeadTemp > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
                queueEntryIndexHeadTemp - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
                queueEntryIndexHeadTemp;

            hlRateControlHistogramPtrTemp = (sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue[queueEntryIndexHeadTemp]);
            if (hlRateControlHistogramPtrTemp->picture_number == picture_control_set_ptr->picture_number &&
                hlRateControlHistogramPtrTemp->passedToHlrc) {
                EbBlockOnMutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
                hlRateControlHistogramPtrTemp->totalNumBitsCoded = picture_control_set_ptr->total_num_bits;
                hlRateControlHistogramPtrTemp->isCoded = EB_TRUE;
                EbReleaseMutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
            }
        }

    }
}
#if ADD_DELTA_QP_SUPPORT ||  NEW_QPS

static const uint8_t quantizer_to_qindex[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48,
    52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100,
    104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152,
    156, 160, 164, 168, 172, 176, 180, 184, 188, 192, 196, 200, 204,
    208, 212, 216, 220, 224, 228, 232, 236, 240, 244, 249, 255,
};
#endif

#if NEW_QPS
#define MAX_Q_INDEX 255
#define MIN_Q_INDEX 0

extern int16_t av1_ac_quant_Q3(int32_t qindex, int32_t delta, aom_bit_depth_t bit_depth);
// These functions use formulaic calculations to make playing with the
// quantizer tables easier. If necessary they can be replaced by lookup
// tables if and when things settle down in the experimental bitstream

double av1_convert_qindex_to_q(int32_t qindex, aom_bit_depth_t bit_depth) {
    // Convert the index to a real Q value (scaled down to match old Q values)
    switch (bit_depth) {
    case AOM_BITS_8: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 4.0;
    case AOM_BITS_10: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 16.0;
    case AOM_BITS_12: return av1_ac_quant_Q3(qindex, 0, bit_depth) / 64.0;
    default:
        assert(0 && "bit_depth should be AOM_BITS_8, AOM_BITS_10 or AOM_BITS_12");
        return -1.0;
    }
}
int32_t av1_compute_qdelta(double qstart, double qtarget,
    aom_bit_depth_t bit_depth) {
    int32_t start_index = MAX_Q_INDEX;
    int32_t target_index = MAX_Q_INDEX;
    int32_t i;

    // Convert the average q value to an index.
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        start_index = i;
        if (av1_convert_qindex_to_q(i, bit_depth) >= qstart) break;
    }

    // Convert the q target to an index
    for (i = MIN_Q_INDEX; i < MAX_Q_INDEX; ++i) {
        target_index = i;
        if (av1_convert_qindex_to_q(i, bit_depth) >= qtarget) break;
    }

    return target_index - start_index;
}
#endif
void* RateControlKernel(void *input_ptr)
{
    // Context
    RateControlContext_t        *context_ptr = (RateControlContext_t*)input_ptr;
    // EncodeContext_t             *encode_context_ptr;

    RateControlIntervalParamContext_t *rateControlParamPtr;

    RateControlIntervalParamContext_t *prevGopRateControlParamPtr;
    RateControlIntervalParamContext_t *nextGopRateControlParamPtr;

    PictureControlSet_t         *picture_control_set_ptr;
    PictureParentControlSet_t   *parentPictureControlSetPtr;

    // Config
    SequenceControlSet_t        *sequence_control_set_ptr;

    // Input
    EbObjectWrapper_t           *rateControlTasksWrapperPtr;
    RateControlTasks_t          *rateControlTasksPtr;

    // Output
    EbObjectWrapper_t           *rateControlResultsWrapperPtr;
    RateControlResults_t        *rateControlResultsPtr;

    RateControlLayerContext_t   *rateControlLayerPtr;

    // SB Loop variables
#if !MEM_RED
    uint32_t                       sb_total_count;
#endif
    LargestCodingUnit_t         *sb_ptr;
    uint32_t                       lcuCodingOrder;
    uint64_t                       totalNumberOfFbFrames = 0;
    uint32_t                       bestOisCuIndex = 0;

    RATE_CONTROL_TASKTYPES       taskType;
    RateControlModel_t          *rc_model_ptr;

    rate_control_model_ctor(&rc_model_ptr);

    for (;;) {

        // Get RateControl Task
        EbGetFullObject(
            context_ptr->rateControlInputTasksFifoPtr,
            &rateControlTasksWrapperPtr);

        rateControlTasksPtr = (RateControlTasks_t*)rateControlTasksWrapperPtr->objectPtr;
        taskType = rateControlTasksPtr->taskType;

        // Modify these for different temporal layers later
        switch (taskType) {

        case RC_PICTURE_MANAGER_RESULT:

            picture_control_set_ptr = (PictureControlSet_t*)rateControlTasksPtr->pictureControlSetWrapperPtr->objectPtr;
            sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;

            // High level RC
            if (picture_control_set_ptr->picture_number == 0) {

                rate_control_model_init(rc_model_ptr, sequence_control_set_ptr);
                context_ptr->highLevelRateControlPtr->target_bit_rate = sequence_control_set_ptr->static_config.target_bit_rate;
                context_ptr->highLevelRateControlPtr->frame_rate = sequence_control_set_ptr->frame_rate;
                context_ptr->highLevelRateControlPtr->channelBitRatePerFrame = (uint64_t)MAX((int64_t)1, (int64_t)((context_ptr->highLevelRateControlPtr->target_bit_rate << RC_PRECISION) / context_ptr->highLevelRateControlPtr->frame_rate));

                context_ptr->highLevelRateControlPtr->channelBitRatePerSw = context_ptr->highLevelRateControlPtr->channelBitRatePerFrame * (sequence_control_set_ptr->static_config.look_ahead_distance + 1);
                context_ptr->highLevelRateControlPtr->bitConstraintPerSw = context_ptr->highLevelRateControlPtr->channelBitRatePerSw;

#if RC_UPDATE_TARGET_RATE
                context_ptr->highLevelRateControlPtr->previousUpdatedBitConstraintPerSw = context_ptr->highLevelRateControlPtr->channelBitRatePerSw;
#endif

                int32_t totalFrameInInterval = sequence_control_set_ptr->intra_period_length;
                uint32_t gopPeriod = (1 << picture_control_set_ptr->parent_pcs_ptr->hierarchical_levels);
                context_ptr->frame_rate = sequence_control_set_ptr->frame_rate;
                while (totalFrameInInterval >= 0) {
                    if (totalFrameInInterval % (gopPeriod) == 0)
                        context_ptr->frames_in_interval[0] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 1) == 0)
                        context_ptr->frames_in_interval[1] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 2) == 0)
                        context_ptr->frames_in_interval[2] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 3) == 0)
                        context_ptr->frames_in_interval[3] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 4) == 0)
                        context_ptr->frames_in_interval[4] ++;
                    else if (totalFrameInInterval % (gopPeriod >> 5) == 0)
                        context_ptr->frames_in_interval[5] ++;
                    totalFrameInInterval--;
                }
                context_ptr->virtualBufferSize = (((uint64_t)sequence_control_set_ptr->static_config.target_bit_rate * 3) << RC_PRECISION) / (context_ptr->frame_rate);
                context_ptr->rateAveragePeriodinFrames = (uint64_t)sequence_control_set_ptr->static_config.intra_period_length + 1;
                context_ptr->virtualBufferLevelInitialValue = context_ptr->virtualBufferSize >> 1;
                context_ptr->virtualBufferLevel = context_ptr->virtualBufferSize >> 1;
                context_ptr->previousVirtualBufferLevel = context_ptr->virtualBufferSize >> 1;
                context_ptr->vbFillThreshold1 = (context_ptr->virtualBufferSize * 6) >> 3;
                context_ptr->vbFillThreshold2 = (context_ptr->virtualBufferSize << 3) >> 3;
                context_ptr->baseLayerFramesAvgQp = sequence_control_set_ptr->qp;
                context_ptr->baseLayerIntraFramesAvgQp = sequence_control_set_ptr->qp;
            }
            if (sequence_control_set_ptr->static_config.rate_control_mode)
            {
                picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp = 0;
                HighLevelRcInputPictureMode2(
                    picture_control_set_ptr->parent_pcs_ptr,
                    sequence_control_set_ptr,
                    sequence_control_set_ptr->encode_context_ptr,
                    context_ptr,
                    context_ptr->highLevelRateControlPtr);


            }

            // Frame level RC
            if (sequence_control_set_ptr->intra_period_length == -1 || sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                rateControlParamPtr = context_ptr->rateControlParamQueue[0];
                prevGopRateControlParamPtr = context_ptr->rateControlParamQueue[0];
                nextGopRateControlParamPtr = context_ptr->rateControlParamQueue[0];
            }
            else {
                uint32_t intervalIndexTemp = 0;
                EbBool intervalFound = EB_FALSE;
                while ((intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER) && !intervalFound) {

                    if (picture_control_set_ptr->picture_number >= context_ptr->rateControlParamQueue[intervalIndexTemp]->firstPoc &&
                        picture_control_set_ptr->picture_number <= context_ptr->rateControlParamQueue[intervalIndexTemp]->lastPoc) {
                        intervalFound = EB_TRUE;
                    }
                    else {
                        intervalIndexTemp++;
                    }
                }
                CHECK_REPORT_ERROR(
                    intervalIndexTemp != PARALLEL_GOP_MAX_NUMBER,
                    sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
                    EB_ENC_RC_ERROR2);

                rateControlParamPtr = context_ptr->rateControlParamQueue[intervalIndexTemp];

                prevGopRateControlParamPtr = (intervalIndexTemp == 0) ?
                    context_ptr->rateControlParamQueue[PARALLEL_GOP_MAX_NUMBER - 1] :
                    context_ptr->rateControlParamQueue[intervalIndexTemp - 1];
                nextGopRateControlParamPtr = (intervalIndexTemp == PARALLEL_GOP_MAX_NUMBER - 1) ?
                    context_ptr->rateControlParamQueue[0] :
                    context_ptr->rateControlParamQueue[intervalIndexTemp + 1];
            }

            rateControlLayerPtr = rateControlParamPtr->rateControlLayerArray[picture_control_set_ptr->temporal_layer_index];
#if !MEM_RED
            sb_total_count = picture_control_set_ptr->sb_total_count;
#endif

            if (sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                // if RC mode is 0,  fixed QP is used
                // QP scaling based on POC number for Flat IPPP structure
#if NEW_QPS
                picture_control_set_ptr->parent_pcs_ptr->base_qindex = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#endif
                if (sequence_control_set_ptr->static_config.enable_qp_scaling_flag && picture_control_set_ptr->parent_pcs_ptr->qp_on_the_fly == EB_FALSE) {
#if NEW_QPS
                    const int32_t qindex = quantizer_to_qindex[(uint8_t)sequence_control_set_ptr->qp];
                    const double q_val = av1_convert_qindex_to_q(qindex, (aom_bit_depth_t)sequence_control_set_ptr->static_config.encoder_bit_depth);
                    if (picture_control_set_ptr->slice_type == I_SLICE) {

                        const int32_t delta_qindex = av1_compute_qdelta(
                            q_val,
                            q_val * 0.25,
                            (aom_bit_depth_t)sequence_control_set_ptr->static_config.encoder_bit_depth);
                        picture_control_set_ptr->parent_pcs_ptr->base_qindex =
                            (uint8_t)CLIP3(
                            (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.min_qp_allowed],
                                (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.max_qp_allowed],
                                (int32_t)(qindex + delta_qindex));

                    }
                    else {
                        // const double delta_rate[FIXED_GF_INTERVAL] = { 0.50, 1.0, 0.85, 1.0,
                        //     0.70, 1.0, 0.85, 1.0 };
                        const double delta_rate_new[6] = { 0.40, 0.7, 0.85, 1.0, 1.0, 1.0 };
                        const int32_t delta_qindex = av1_compute_qdelta(
                            q_val,
                            q_val * delta_rate_new[picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index],
                            (aom_bit_depth_t)sequence_control_set_ptr->static_config.encoder_bit_depth);

                        picture_control_set_ptr->parent_pcs_ptr->base_qindex =
                            (uint8_t)CLIP3(
                            (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.min_qp_allowed],
                                (int32_t)quantizer_to_qindex[sequence_control_set_ptr->static_config.max_qp_allowed],
                                (int32_t)(qindex + delta_qindex));

                    }
#endif

                    if (picture_control_set_ptr->slice_type == I_SLICE) {
                        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, (int32_t)(sequence_control_set_ptr->qp) + context_ptr->maxRateAdjustDeltaQP);
                    }
                    else {

#if V2_QP_SCALING || NEW_QP_SCALING
#if V2_QP_SCALING
                        if (/*picture_control_set_ptr->enc_mode == ENC_M0 &&*/ 0/*sequence_control_set_ptr->static_config.tune == TUNE_VMAF*/) {
#endif
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, (int32_t)(sequence_control_set_ptr->qp + DEFAULT_QP_OFFSET_LAYER_ARRAY[picture_control_set_ptr->temporal_layer_index]));
                        }
                        else {
#if NEW_QP_SCALING
#if ENCODER_MODE_CLEANUP
                            if(1){
#else
                            if (picture_control_set_ptr->enc_mode <= ENC_M1 /*&&sequence_control_set_ptr->static_config.tune != TUNE_VQ*/) {
#endif
                                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, (int32_t)(sequence_control_set_ptr->qp + QP_OFFSET_LAYER_ARRAY_BDRATE[picture_control_set_ptr->temporal_layer_index]));
                            }
                            else {
                                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, (int32_t)(sequence_control_set_ptr->qp + QP_OFFSET_LAYER_ARRAY[picture_control_set_ptr->temporal_layer_index]));
                            }
#else
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, (int32_t)(sequence_control_set_ptr->qp + QP_OFFSET_LAYER_ARRAY[picture_control_set_ptr->temporal_layer_index]));
#endif
#if V2_QP_SCALING
                        }
#endif
#else
                        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, (int32_t)(sequence_control_set_ptr->qp + QP_OFFSET_LAYER_ARRAY[picture_control_set_ptr->temporal_layer_index]));
#endif
                    }
                }

                else if (picture_control_set_ptr->parent_pcs_ptr->qp_on_the_fly == EB_TRUE) {

                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3((int32_t)sequence_control_set_ptr->static_config.min_qp_allowed, (int32_t)sequence_control_set_ptr->static_config.max_qp_allowed, picture_control_set_ptr->parent_pcs_ptr->picture_qp);
#if NEW_QPS
                    picture_control_set_ptr->parent_pcs_ptr->base_qindex = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#endif
                }

                picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp;

            }
            else {
                FrameLevelRcInputPictureMode2(
                    picture_control_set_ptr,
                    sequence_control_set_ptr,
                    context_ptr,
                    rateControlLayerPtr,
                    rateControlParamPtr,
                    bestOisCuIndex);

                picture_control_set_ptr->picture_qp = rate_control_get_quantizer(rc_model_ptr, picture_control_set_ptr->parent_pcs_ptr);

                if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc && picture_control_set_ptr->picture_number != 0 && !prevGopRateControlParamPtr->scene_change_in_gop) {
                    int16_t deltaApQp = (int16_t)prevGopRateControlParamPtr->firstPicActualQp - (int16_t)prevGopRateControlParamPtr->firstPicPredQp;
                    int64_t extraApBitRatio = (prevGopRateControlParamPtr->firstPicPredBits != 0) ?
                        (((int64_t)prevGopRateControlParamPtr->firstPicActualBits - (int64_t)prevGopRateControlParamPtr->firstPicPredBits) * 100) / ((int64_t)prevGopRateControlParamPtr->firstPicPredBits) :
                        0;
                    extraApBitRatio += (int64_t)deltaApQp * 15;
                    if (extraApBitRatio > 200) {
                        picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 3;
                    }
                    else if (extraApBitRatio > 100) {
                        picture_control_set_ptr->picture_qp = picture_control_set_ptr->picture_qp + 2;
                    }
                    else if (extraApBitRatio > 50) {
                        picture_control_set_ptr->picture_qp++;
                    }
                }

                if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc && picture_control_set_ptr->picture_number != 0) {
                    uint8_t qpIncAllowed = 3;
                    uint8_t qpDecAllowed = 4;
                    if (picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp + 10 <= prevGopRateControlParamPtr->firstPicActualQp)
                    {
                        qpDecAllowed = (uint8_t)(prevGopRateControlParamPtr->firstPicActualQp - picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp) >> 1;
                    }

                    if (picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp >= prevGopRateControlParamPtr->firstPicActualQp + 10)
                    {
                        qpIncAllowed = (uint8_t)(picture_control_set_ptr->parent_pcs_ptr->intra_selected_org_qp - prevGopRateControlParamPtr->firstPicActualQp) * 2 / 3;
                        if (prevGopRateControlParamPtr->firstPicActualQp <= 15)
                            qpIncAllowed += 5;
                        else if (prevGopRateControlParamPtr->firstPicActualQp <= 20)
                            qpIncAllowed += 4;
                        else if (prevGopRateControlParamPtr->firstPicActualQp <= 25)
                            qpIncAllowed += 3;
                    }
                    else if (prevGopRateControlParamPtr->scene_change_in_gop) {
                        qpIncAllowed = 5;
                    }
                    if (picture_control_set_ptr->parent_pcs_ptr->end_of_sequence_region) {
                        qpIncAllowed += 2;
                        qpDecAllowed += 4;
                    }
                    picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                        (uint32_t)MAX((int32_t)prevGopRateControlParamPtr->firstPicActualQp - (int32_t)qpDecAllowed, 0),
                        (uint32_t)prevGopRateControlParamPtr->firstPicActualQp + qpIncAllowed,
                        picture_control_set_ptr->picture_qp);
                }

                // Scene change
                if (picture_control_set_ptr->slice_type == I_SLICE && picture_control_set_ptr->picture_number != rateControlParamPtr->firstPoc) {
                    if (nextGopRateControlParamPtr->firstPicActualQpAssigned) {

                        picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                            (uint32_t)MAX((int32_t)nextGopRateControlParamPtr->firstPicActualQp - (int32_t)1, 0),
                            (uint32_t)nextGopRateControlParamPtr->firstPicActualQp + 8,
                            picture_control_set_ptr->picture_qp);
                    }
                    else {
                        if (rateControlParamPtr->firstPicActualQp < 20) {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)rateControlParamPtr->firstPicActualQp - (int32_t)4, 0),
                                (uint32_t)rateControlParamPtr->firstPicActualQp + 10,
                                picture_control_set_ptr->picture_qp);
                        }
                        else {
                            picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                                (uint32_t)MAX((int32_t)rateControlParamPtr->firstPicActualQp - (int32_t)4, 0),
                                (uint32_t)rateControlParamPtr->firstPicActualQp + 8,
                                picture_control_set_ptr->picture_qp);

                        }

                    }
                }
                picture_control_set_ptr->picture_qp = (uint8_t)CLIP3(
                    sequence_control_set_ptr->static_config.min_qp_allowed,
                    sequence_control_set_ptr->static_config.max_qp_allowed,
                    picture_control_set_ptr->picture_qp);
#if NEW_QPS
                picture_control_set_ptr->parent_pcs_ptr->base_qindex = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#endif
            }

            picture_control_set_ptr->parent_pcs_ptr->picture_qp = picture_control_set_ptr->picture_qp;
            if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0 && sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                context_ptr->baseLayerFramesAvgQp = (3 * context_ptr->baseLayerFramesAvgQp + picture_control_set_ptr->picture_qp + 2) >> 2;
            }
            if (picture_control_set_ptr->slice_type == I_SLICE) {
                if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc) {
                    rateControlParamPtr->firstPicPredQp = (uint16_t)picture_control_set_ptr->parent_pcs_ptr->best_pred_qp;
                    rateControlParamPtr->firstPicActualQp = (uint16_t)picture_control_set_ptr->picture_qp;
                    rateControlParamPtr->scene_change_in_gop = picture_control_set_ptr->parent_pcs_ptr->scene_change_in_gop;
                    rateControlParamPtr->firstPicActualQpAssigned = EB_TRUE;
                }
                {
                    if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc) {
                        if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                            context_ptr->baseLayerIntraFramesAvgQp = (3 * context_ptr->baseLayerIntraFramesAvgQp + picture_control_set_ptr->picture_qp + 2) >> 2;
                        }
                    }

                    if (picture_control_set_ptr->picture_number == rateControlParamPtr->firstPoc) {
                        rateControlParamPtr->intraFramesQp = picture_control_set_ptr->picture_qp;
                        rateControlParamPtr->nextGopIntraFrameQp = picture_control_set_ptr->picture_qp;

                    }
                }
            }
            picture_control_set_ptr->parent_pcs_ptr->average_qp = 0;
#if MEM_RED
            for (lcuCodingOrder = 0; lcuCodingOrder < sequence_control_set_ptr->sb_tot_cnt; ++lcuCodingOrder) {
#else

            for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
#endif

                sb_ptr = picture_control_set_ptr->sb_ptr_array[lcuCodingOrder];
#if ADD_DELTA_QP_SUPPORT

                sb_ptr->qp = quantizer_to_qindex[picture_control_set_ptr->picture_qp];
#else
                sb_ptr->qp = (uint8_t)picture_control_set_ptr->picture_qp;
#endif
                picture_control_set_ptr->parent_pcs_ptr->average_qp += sb_ptr->qp;
            }

            // Get Empty Rate Control Results Buffer
            EbGetEmptyObject(
                context_ptr->rateControlOutputResultsFifoPtr,
                &rateControlResultsWrapperPtr);
            rateControlResultsPtr = (RateControlResults_t*)rateControlResultsWrapperPtr->objectPtr;
            rateControlResultsPtr->pictureControlSetWrapperPtr = rateControlTasksPtr->pictureControlSetWrapperPtr;

            // Post Full Rate Control Results
            EbPostFullObject(rateControlResultsWrapperPtr);

            // Release Rate Control Tasks
            EbReleaseObject(rateControlTasksWrapperPtr);

            break;

        case RC_PACKETIZATION_FEEDBACK_RESULT:

            parentPictureControlSetPtr = (PictureParentControlSet_t*)rateControlTasksPtr->pictureControlSetWrapperPtr->objectPtr;
            sequence_control_set_ptr = (SequenceControlSet_t*)parentPictureControlSetPtr->sequence_control_set_wrapper_ptr->objectPtr;

            if (sequence_control_set_ptr->static_config.rate_control_mode) {
                rate_control_update_model(rc_model_ptr, parentPictureControlSetPtr);
            }

            // Frame level RC
            if (sequence_control_set_ptr->intra_period_length == -1 || sequence_control_set_ptr->static_config.rate_control_mode == 0) {
                rateControlParamPtr = context_ptr->rateControlParamQueue[0];
                prevGopRateControlParamPtr = context_ptr->rateControlParamQueue[0];
                if (parentPictureControlSetPtr->slice_type == I_SLICE) {

                    if (parentPictureControlSetPtr->total_num_bits > MAX_BITS_PER_FRAME) {
                        context_ptr->maxRateAdjustDeltaQP++;
                    }
                    else if (context_ptr->maxRateAdjustDeltaQP > 0 && parentPictureControlSetPtr->total_num_bits < MAX_BITS_PER_FRAME * 85 / 100) {
                        context_ptr->maxRateAdjustDeltaQP--;
                    }
                    context_ptr->maxRateAdjustDeltaQP = CLIP3(0, 63, context_ptr->maxRateAdjustDeltaQP);
                    context_ptr->maxRateAdjustDeltaQP = 0;
                }
            }
            else {
                uint32_t intervalIndexTemp = 0;
                EbBool intervalFound = EB_FALSE;
                while ((intervalIndexTemp < PARALLEL_GOP_MAX_NUMBER) && !intervalFound) {

                    if (parentPictureControlSetPtr->picture_number >= context_ptr->rateControlParamQueue[intervalIndexTemp]->firstPoc &&
                        parentPictureControlSetPtr->picture_number <= context_ptr->rateControlParamQueue[intervalIndexTemp]->lastPoc) {
                        intervalFound = EB_TRUE;
                    }
                    else {
                        intervalIndexTemp++;
                    }
                }
                CHECK_REPORT_ERROR(
                    intervalIndexTemp != PARALLEL_GOP_MAX_NUMBER,
                    sequence_control_set_ptr->encode_context_ptr->app_callback_ptr,
                    EB_ENC_RC_ERROR2);

                rateControlParamPtr = context_ptr->rateControlParamQueue[intervalIndexTemp];

                prevGopRateControlParamPtr = (intervalIndexTemp == 0) ?
                    context_ptr->rateControlParamQueue[PARALLEL_GOP_MAX_NUMBER - 1] :
                    context_ptr->rateControlParamQueue[intervalIndexTemp - 1];

            }
            if (sequence_control_set_ptr->static_config.rate_control_mode != 0) {

                context_ptr->previousVirtualBufferLevel = context_ptr->virtualBufferLevel;

                context_ptr->virtualBufferLevel =
                    (int64_t)context_ptr->previousVirtualBufferLevel +
                    (int64_t)parentPictureControlSetPtr->total_num_bits - (int64_t)context_ptr->highLevelRateControlPtr->channelBitRatePerFrame;

                HighLevelRcFeedBackPicture(
                    parentPictureControlSetPtr,
                    sequence_control_set_ptr);
                FrameLevelRcFeedbackPictureMode2(
                    parentPictureControlSetPtr,
                    sequence_control_set_ptr,
                    context_ptr);
                if (parentPictureControlSetPtr->picture_number == rateControlParamPtr->firstPoc) {
                    rateControlParamPtr->firstPicPredBits = parentPictureControlSetPtr->target_bits_best_pred_qp;
                    rateControlParamPtr->firstPicActualBits = parentPictureControlSetPtr->total_num_bits;
                    {
                        int16_t deltaApQp = (int16_t)rateControlParamPtr->firstPicActualQp - (int16_t)rateControlParamPtr->firstPicPredQp;
                        rateControlParamPtr->extraApBitRatioI = (rateControlParamPtr->firstPicPredBits != 0) ?
                            (((int64_t)rateControlParamPtr->firstPicActualBits - (int64_t)rateControlParamPtr->firstPicPredBits) * 100) / ((int64_t)rateControlParamPtr->firstPicPredBits) :
                            0;
                        rateControlParamPtr->extraApBitRatioI += (int64_t)deltaApQp * 15;
                    }


                }

            }

            // Queue variables
#if OVERSHOOT_STAT_PRINT
            if (sequence_control_set_ptr->intra_period_length != -1) {

                int32_t                       queueEntryIndex;
                uint32_t                       queueEntryIndexTemp;
                uint32_t                       queueEntryIndexTemp2;
                CodedFramesStatsEntry_t     *queueEntryPtr;
                EbBool                      moveSlideWondowFlag = EB_TRUE;
                EbBool                      end_of_sequence_flag = EB_TRUE;
                uint32_t                       frames_in_sw;

                // Determine offset from the Head Ptr
                queueEntryIndex = (int32_t)(parentPictureControlSetPtr->picture_number - context_ptr->codedFramesStatQueue[context_ptr->codedFramesStatQueueHeadIndex]->picture_number);
                queueEntryIndex += context_ptr->codedFramesStatQueueHeadIndex;
                queueEntryIndex = (queueEntryIndex > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queueEntryIndex - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queueEntryIndex;
                queueEntryPtr = context_ptr->codedFramesStatQueue[queueEntryIndex];

                queueEntryPtr->frameTotalBitActual = (uint64_t)parentPictureControlSetPtr->total_num_bits;
                queueEntryPtr->picture_number = parentPictureControlSetPtr->picture_number;
                queueEntryPtr->end_of_sequence_flag = parentPictureControlSetPtr->end_of_sequence_flag;
                context_ptr->rateAveragePeriodinFrames = (uint64_t)sequence_control_set_ptr->static_config.intra_period_length + 1;

                //printf("\n0_POC: %d\n",
                //    queueEntryPtr->picture_number);
                moveSlideWondowFlag = EB_TRUE;
                while (moveSlideWondowFlag) {
                    //  printf("\n1_POC: %d\n",
                    //      queueEntryPtr->picture_number);
                      // Check if the sliding window condition is valid
                    queueEntryIndexTemp = context_ptr->codedFramesStatQueueHeadIndex;
                    if (context_ptr->codedFramesStatQueue[queueEntryIndexTemp]->frameTotalBitActual != -1) {
                        end_of_sequence_flag = context_ptr->codedFramesStatQueue[queueEntryIndexTemp]->end_of_sequence_flag;
                    }
                    else {
                        end_of_sequence_flag = EB_FALSE;
                    }
                    while (moveSlideWondowFlag && !end_of_sequence_flag &&
                        queueEntryIndexTemp < context_ptr->codedFramesStatQueueHeadIndex + context_ptr->rateAveragePeriodinFrames) {
                        // printf("\n2_POC: %d\n",
                        //     queueEntryPtr->picture_number);

                        queueEntryIndexTemp2 = (queueEntryIndexTemp > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queueEntryIndexTemp;

                        moveSlideWondowFlag = (EbBool)(moveSlideWondowFlag && (context_ptr->codedFramesStatQueue[queueEntryIndexTemp2]->frameTotalBitActual != -1));

                        if (context_ptr->codedFramesStatQueue[queueEntryIndexTemp2]->frameTotalBitActual != -1) {
                            // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                            end_of_sequence_flag = context_ptr->codedFramesStatQueue[queueEntryIndexTemp]->end_of_sequence_flag;
                        }
                        else {
                            end_of_sequence_flag = EB_FALSE;
                        }
                        queueEntryIndexTemp =
                            (queueEntryIndexTemp == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : queueEntryIndexTemp + 1;

                    }

                    if (moveSlideWondowFlag) {
                        //get a new entry spot
                        queueEntryPtr = (context_ptr->codedFramesStatQueue[context_ptr->codedFramesStatQueueHeadIndex]);
                        queueEntryIndexTemp = context_ptr->codedFramesStatQueueHeadIndex;
                        // This is set to false, so the last frame would go inside the loop
                        end_of_sequence_flag = EB_FALSE;
                        frames_in_sw = 0;
                        context_ptr->totalBitActualPerSw = 0;

                        while (!end_of_sequence_flag &&
                            queueEntryIndexTemp < context_ptr->codedFramesStatQueueHeadIndex + context_ptr->rateAveragePeriodinFrames) {
                            frames_in_sw++;

                            queueEntryIndexTemp2 = (queueEntryIndexTemp > CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - CODED_FRAMES_STAT_QUEUE_MAX_DEPTH : queueEntryIndexTemp;

                            context_ptr->totalBitActualPerSw += context_ptr->codedFramesStatQueue[queueEntryIndexTemp2]->frameTotalBitActual;
                            end_of_sequence_flag = context_ptr->codedFramesStatQueue[queueEntryIndexTemp2]->end_of_sequence_flag;

                            queueEntryIndexTemp =
                                (queueEntryIndexTemp == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : queueEntryIndexTemp + 1;

                        }
                        //

                        //if(frames_in_sw == context_ptr->rateAveragePeriodinFrames)
                        //    printf("POC:%d\t %.3f\n", queueEntryPtr->picture_number, (double)context_ptr->totalBitActualPerSw*(sequence_control_set_ptr->frame_rate>> RC_PRECISION)/(double)frames_in_sw/1000);
                        if (frames_in_sw == (uint32_t)sequence_control_set_ptr->intra_period_length + 1) {
                            context_ptr->maxBitActualPerSw = MAX(context_ptr->maxBitActualPerSw, context_ptr->totalBitActualPerSw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / 1000);
                            if (queueEntryPtr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                                context_ptr->maxBitActualPerGop = MAX(context_ptr->maxBitActualPerGop, context_ptr->totalBitActualPerSw*(sequence_control_set_ptr->frame_rate >> RC_PRECISION) / frames_in_sw / 1000);
                                if (context_ptr->totalBitActualPerSw > sequence_control_set_ptr->static_config.maxBufferSize) {
                                    printf("\nPOC:%d\tOvershoot:%.0f%% \n",
                                        (int32_t)queueEntryPtr->picture_number,
                                        (double)((int64_t)context_ptr->totalBitActualPerSw * 100 / (int64_t)sequence_control_set_ptr->static_config.maxBufferSize - 100));
                                }
                            }
                        }
                        if (frames_in_sw == context_ptr->rateAveragePeriodinFrames - 1) {
                            printf("\n%d MAX\n", (int32_t)context_ptr->maxBitActualPerSw);
                            printf("\n%d GopMa\n", (int32_t)context_ptr->maxBitActualPerGop);
                        }

                        // Reset the Queue Entry
                        queueEntryPtr->picture_number += CODED_FRAMES_STAT_QUEUE_MAX_DEPTH;
                        queueEntryPtr->frameTotalBitActual = -1;

                        // Increment the Reorder Queue head Ptr
                        context_ptr->codedFramesStatQueueHeadIndex =
                            (context_ptr->codedFramesStatQueueHeadIndex == CODED_FRAMES_STAT_QUEUE_MAX_DEPTH - 1) ? 0 : context_ptr->codedFramesStatQueueHeadIndex + 1;

                        queueEntryPtr = (context_ptr->codedFramesStatQueue[context_ptr->codedFramesStatQueueHeadIndex]);

                    }
                }
            }
#endif
            totalNumberOfFbFrames++;

            // Release the SequenceControlSet
            EbReleaseObject(parentPictureControlSetPtr->sequence_control_set_wrapper_ptr);

            // Release the input buffer 
            EbReleaseObject(parentPictureControlSetPtr->input_picture_wrapper_ptr);

            // Release the ParentPictureControlSet
            EbReleaseObject(rateControlTasksPtr->pictureControlSetWrapperPtr);

            // Release Rate Control Tasks
            EbReleaseObject(rateControlTasksWrapperPtr);

            break;

        case RC_ENTROPY_CODING_ROW_FEEDBACK_RESULT:

            // Extract bits-per-lcu-row

            // Release Rate Control Tasks
            EbReleaseObject(rateControlTasksWrapperPtr);

            break;

        default:
            picture_control_set_ptr = (PictureControlSet_t*)rateControlTasksPtr->pictureControlSetWrapperPtr->objectPtr;
            sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
            //encode_context_ptr            = sequence_control_set_ptr->encode_context_ptr;
            //CHECK_REPORT_ERROR_NC(
            //             encode_context_ptr->app_callback_ptr,
            //             EB_ENC_RC_ERROR1);

            break;
        }
    }
    return EB_NULL;
}
