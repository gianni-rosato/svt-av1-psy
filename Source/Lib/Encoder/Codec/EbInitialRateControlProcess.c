// clang-format off
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

/**************************************
 * Context
 **************************************/
typedef struct InitialRateControlContext
{
    EbFifo                    *motion_estimation_results_input_fifo_ptr;
    EbFifo                    *initialrate_control_results_output_fifo_ptr;
}InitialRateControlContext;


/**************************************
* Macros
**************************************/
#define PAN_LCU_PERCENTAGE                    75
#define LOW_AMPLITUDE_TH                      16

void GetMv(
    PictureParentControlSet    *picture_control_set_ptr,
    uint32_t                         sb_index,
    int32_t                        *xCurrentMv,
    int32_t                        *yCurrentMv)
{
    uint32_t             meCandidateIndex;

    const MeLcuResults *me_results = picture_control_set_ptr->me_results[sb_index];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[0];
    const MeCandidate *me_block_results = me_results->me_candidate[0];
    for (meCandidateIndex = 0; meCandidateIndex < total_me_cnt; meCandidateIndex++) {
        if (me_block_results->direction == UNI_PRED_LIST_0) {
            *xCurrentMv = me_results->me_mv_array[0][0].x_mv;
            *yCurrentMv = me_results->me_mv_array[0][0].y_mv;
            break;
        }
    }
}

void GetMeDist(
    PictureParentControlSet    *picture_control_set_ptr,
    uint32_t                         sb_index,
    uint32_t                      *distortion)
{
    *distortion = (uint32_t)picture_control_set_ptr->me_results[sb_index]->me_candidate[0][0].distortion;
}

EbBool CheckMvForPanHighAmp(
    uint32_t   hierarchical_levels,
    uint32_t     temporal_layer_index,
    int32_t    *xCurrentMv,
    int32_t    *xCandidateMv)
{
    if (*xCurrentMv * *xCandidateMv > 0                        // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*xCurrentMv) >= global_motion_threshold[hierarchical_levels][temporal_layer_index]    // high amplitude
        && ABS(*xCandidateMv) >= global_motion_threshold[hierarchical_levels][temporal_layer_index]    // high amplitude
        && ABS(*xCurrentMv - *xCandidateMv) < LOW_AMPLITUDE_TH) {    // close amplitude

        return(EB_TRUE);
    }

    else
        return(EB_FALSE);
}

EbBool CheckMvForTiltHighAmp(
    uint32_t   hierarchical_levels,
    uint32_t     temporal_layer_index,
    int32_t    *yCurrentMv,
    int32_t    *yCandidateMv)
{
    if (*yCurrentMv * *yCandidateMv > 0                        // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*yCurrentMv) >= global_motion_threshold[hierarchical_levels][temporal_layer_index]    // high amplitude
        && ABS(*yCandidateMv) >= global_motion_threshold[hierarchical_levels][temporal_layer_index]    // high amplitude
        && ABS(*yCurrentMv - *yCandidateMv) < LOW_AMPLITUDE_TH) {    // close amplitude

        return(EB_TRUE);
    }

    else
        return(EB_FALSE);
}

EbBool CheckMvForPan(
    uint32_t   hierarchical_levels,
    uint32_t     temporal_layer_index,
    int32_t    *xCurrentMv,
    int32_t    *yCurrentMv,
    int32_t    *xCandidateMv,
    int32_t    *yCandidateMv)
{
    if (*yCurrentMv < LOW_AMPLITUDE_TH
        && *yCandidateMv < LOW_AMPLITUDE_TH
        && *xCurrentMv * *xCandidateMv        > 0                        // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*xCurrentMv) >= global_motion_threshold[hierarchical_levels][temporal_layer_index]    // high amplitude
        && ABS(*xCandidateMv) >= global_motion_threshold[hierarchical_levels][temporal_layer_index]    // high amplitude
        && ABS(*xCurrentMv - *xCandidateMv) < LOW_AMPLITUDE_TH) {    // close amplitude

        return(EB_TRUE);
    }

    else
        return(EB_FALSE);
}

EbBool CheckMvForTilt(
    uint32_t   hierarchical_levels,
    uint32_t     temporal_layer_index,
    int32_t    *xCurrentMv,
    int32_t    *yCurrentMv,
    int32_t    *xCandidateMv,
    int32_t    *yCandidateMv)
{
    if (*xCurrentMv < LOW_AMPLITUDE_TH
        && *xCandidateMv < LOW_AMPLITUDE_TH
        && *yCurrentMv * *yCandidateMv        > 0                        // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*yCurrentMv) >= global_motion_threshold[hierarchical_levels][temporal_layer_index]    // high amplitude
        && ABS(*yCandidateMv) >= global_motion_threshold[hierarchical_levels][temporal_layer_index]    // high amplitude
        && ABS(*yCurrentMv - *yCandidateMv) < LOW_AMPLITUDE_TH) {    // close amplitude

        return(EB_TRUE);
    }

    else
        return(EB_FALSE);
}

EbBool CheckMvForNonUniformMotion(
    int32_t    *xCurrentMv,
    int32_t    *yCurrentMv,
    int32_t    *xCandidateMv,
    int32_t    *yCandidateMv)
{
    int32_t mvThreshold = 40;//LOW_AMPLITUDE_TH + 18;
    // Either the x or the y direction is greater than threshold
    if ((ABS(*xCurrentMv - *xCandidateMv) > mvThreshold) || (ABS(*yCurrentMv - *yCandidateMv) > mvThreshold))
        return(EB_TRUE);
    else
        return(EB_FALSE);
}

void CheckForNonUniformMotionVectorField(
    PictureParentControlSet    *picture_control_set_ptr)
{
    uint32_t    sb_count;
    uint32_t    picture_width_in_sb = (picture_control_set_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t    sb_origin_x;
    uint32_t    sb_origin_y;

    int32_t    xCurrentMv = 0;
    int32_t    yCurrentMv = 0;
    int32_t    xLeftMv = 0;
    int32_t    yLeftMv = 0;
    int32_t    xTopMv = 0;
    int32_t    yTopMv = 0;
    int32_t    xRightMv = 0;
    int32_t    yRightMv = 0;
    int32_t    xBottomMv = 0;
    int32_t    yBottomMv = 0;
    uint32_t countOfNonUniformNeighbors = 0;

    for (sb_count = 0; sb_count < picture_control_set_ptr->sb_total_count; ++sb_count) {
        countOfNonUniformNeighbors = 0;

        sb_origin_x = (sb_count % picture_width_in_sb) * BLOCK_SIZE_64;
        sb_origin_y = (sb_count / picture_width_in_sb) * BLOCK_SIZE_64;

        if (((sb_origin_x + BLOCK_SIZE_64) <= picture_control_set_ptr->enhanced_picture_ptr->width) &&
            ((sb_origin_y + BLOCK_SIZE_64) <= picture_control_set_ptr->enhanced_picture_ptr->height)) {
            // Current MV
            GetMv(picture_control_set_ptr, sb_count, &xCurrentMv, &yCurrentMv);

            // Left MV
            if (sb_origin_x == 0) {
                xLeftMv = 0;
                yLeftMv = 0;
            }
            else
                GetMv(picture_control_set_ptr, sb_count - 1, &xLeftMv, &yLeftMv);
            countOfNonUniformNeighbors += CheckMvForNonUniformMotion(&xCurrentMv, &yCurrentMv, &xLeftMv, &yLeftMv);

            // Top MV
            if (sb_origin_y == 0) {
                xTopMv = 0;
                yTopMv = 0;
            }
            else
                GetMv(picture_control_set_ptr, sb_count - picture_width_in_sb, &xTopMv, &yTopMv);
            countOfNonUniformNeighbors += CheckMvForNonUniformMotion(&xCurrentMv, &yCurrentMv, &xTopMv, &yTopMv);

            // Right MV
            if ((sb_origin_x + (BLOCK_SIZE_64 << 1)) > picture_control_set_ptr->enhanced_picture_ptr->width) {
                xRightMv = 0;
                yRightMv = 0;
            }
            else
                GetMv(picture_control_set_ptr, sb_count + 1, &xRightMv, &yRightMv);
            countOfNonUniformNeighbors += CheckMvForNonUniformMotion(&xCurrentMv, &yCurrentMv, &xRightMv, &yRightMv);

            // Bottom MV
            if ((sb_origin_y + (BLOCK_SIZE_64 << 1)) > picture_control_set_ptr->enhanced_picture_ptr->height) {
                xBottomMv = 0;
                yBottomMv = 0;
            }
            else
                GetMv(picture_control_set_ptr, sb_count + picture_width_in_sb, &xBottomMv, &yBottomMv);
            countOfNonUniformNeighbors += CheckMvForNonUniformMotion(&xCurrentMv, &yCurrentMv, &xBottomMv, &yBottomMv);
        }
    }
}


void DetectGlobalMotion(
    PictureParentControlSet    *picture_control_set_ptr)
{

#if  INIT_GM_FIX
    //initilize global motion to be OFF for all references frames.
    memset(picture_control_set_ptr->is_global_motion, EB_FALSE, MAX_NUM_OF_REF_PIC_LIST*REF_LIST_MAX_DEPTH);
#endif

#if GLOBAL_WARPED_MOTION
#if GM_OPT
    if (picture_control_set_ptr->gm_level <= GM_DOWN) {
#endif
    uint32_t numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE)
        ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;

    for (uint32_t listIndex = REF_LIST_0; listIndex <= numOfListToSearch; ++listIndex) {

        uint32_t num_of_ref_pic_to_search;
        if (picture_control_set_ptr->is_alt_ref == EB_TRUE)
            num_of_ref_pic_to_search = 1;
        else
            num_of_ref_pic_to_search = picture_control_set_ptr->slice_type == P_SLICE
                ? picture_control_set_ptr->ref_list0_count
                : listIndex == REF_LIST_0
                    ? picture_control_set_ptr->ref_list0_count
                    : picture_control_set_ptr->ref_list1_count;

        // Ref Picture Loop
        for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
             ++ref_pic_index)
        {
            picture_control_set_ptr->is_global_motion[listIndex][ref_pic_index] = EB_FALSE;
            if (picture_control_set_ptr->global_motion_estimation[listIndex][ref_pic_index].wmtype > TRANSLATION)
                picture_control_set_ptr->is_global_motion[listIndex][ref_pic_index] = EB_TRUE;
        }
    }
#endif
#if GM_OPT && GLOBAL_WARPED_MOTION || !GLOBAL_WARPED_MOTION
#if GM_OPT && GLOBAL_WARPED_MOTION
    }
    else {
#endif
    uint32_t    sb_count;
    uint32_t    picture_width_in_sb = (picture_control_set_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t    sb_origin_x;
    uint32_t    sb_origin_y;

    uint32_t  totalCheckedLcus = 0;
    uint32_t  totalPanLcus = 0;

    int32_t    xCurrentMv = 0;
    int32_t    yCurrentMv = 0;
    int32_t    xLeftMv = 0;
    int32_t    yLeftMv = 0;
    int32_t    xTopMv = 0;
    int32_t    yTopMv = 0;
    int32_t    xRightMv = 0;
    int32_t    yRightMv = 0;
    int32_t    xBottomMv = 0;
    int32_t    yBottomMv = 0;
    int64_t  xTiltMvSum = 0;
    int64_t  yTiltMvSum = 0;
    int64_t xPanMvSum = 0;
    int64_t yPanMvSum = 0;
    uint32_t  totalTiltLcus = 0;

    uint32_t  totalTiltHighAmpLcus = 0;
    uint32_t  totalPanHighAmpLcus = 0;

    for (sb_count = 0; sb_count < picture_control_set_ptr->sb_total_count; ++sb_count) {
        sb_origin_x = (sb_count % picture_width_in_sb) * BLOCK_SIZE_64;
        sb_origin_y = (sb_count / picture_width_in_sb) * BLOCK_SIZE_64;
        if (((sb_origin_x + BLOCK_SIZE_64) <= picture_control_set_ptr->enhanced_picture_ptr->width) &&
            ((sb_origin_y + BLOCK_SIZE_64) <= picture_control_set_ptr->enhanced_picture_ptr->height)) {
            // Current MV
            GetMv(picture_control_set_ptr, sb_count, &xCurrentMv, &yCurrentMv);

            // Left MV
            if (sb_origin_x == 0) {
                xLeftMv = 0;
                yLeftMv = 0;
            }
            else
                GetMv(picture_control_set_ptr, sb_count - 1, &xLeftMv, &yLeftMv);
            // Top MV
            if (sb_origin_y == 0) {
                xTopMv = 0;
                yTopMv = 0;
            }
            else
                GetMv(picture_control_set_ptr, sb_count - picture_width_in_sb, &xTopMv, &yTopMv);
            // Right MV
            if ((sb_origin_x + (BLOCK_SIZE_64 << 1)) > picture_control_set_ptr->enhanced_picture_ptr->width) {
                xRightMv = 0;
                yRightMv = 0;
            }
            else
                GetMv(picture_control_set_ptr, sb_count + 1, &xRightMv, &yRightMv);
            // Bottom MV
            if ((sb_origin_y + (BLOCK_SIZE_64 << 1)) > picture_control_set_ptr->enhanced_picture_ptr->height) {
                xBottomMv = 0;
                yBottomMv = 0;
            }
            else
                GetMv(picture_control_set_ptr, sb_count + picture_width_in_sb, &xBottomMv, &yBottomMv);
            totalCheckedLcus++;

            if ((EbBool)(CheckMvForPan(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &yCurrentMv, &xLeftMv, &yLeftMv) ||
                CheckMvForPan(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &yCurrentMv, &xTopMv, &yTopMv) ||
                CheckMvForPan(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &yCurrentMv, &xRightMv, &yRightMv) ||
                CheckMvForPan(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &yCurrentMv, &xBottomMv, &yBottomMv))) {
                totalPanLcus++;

                xPanMvSum += xCurrentMv;
                yPanMvSum += yCurrentMv;
            }

            if ((EbBool)(CheckMvForTilt(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &yCurrentMv, &xLeftMv, &yLeftMv) ||
                CheckMvForTilt(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &yCurrentMv, &xTopMv, &yTopMv) ||
                CheckMvForTilt(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &yCurrentMv, &xRightMv, &yRightMv) ||
                CheckMvForTilt(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &yCurrentMv, &xBottomMv, &yBottomMv))) {
                totalTiltLcus++;

                xTiltMvSum += xCurrentMv;
                yTiltMvSum += yCurrentMv;
            }

            if ((EbBool)(CheckMvForPanHighAmp(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &xLeftMv) ||
                CheckMvForPanHighAmp(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &xTopMv) ||
                CheckMvForPanHighAmp(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &xRightMv) ||
                CheckMvForPanHighAmp(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &xCurrentMv, &xBottomMv))) {
                totalPanHighAmpLcus++;
            }

            if ((EbBool)(CheckMvForTiltHighAmp(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &yCurrentMv, &yLeftMv) ||
                CheckMvForTiltHighAmp(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &yCurrentMv, &yTopMv) ||
                CheckMvForTiltHighAmp(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &yCurrentMv, &yRightMv) ||
                CheckMvForTiltHighAmp(picture_control_set_ptr->hierarchical_levels, picture_control_set_ptr->temporal_layer_index, &yCurrentMv, &yBottomMv))) {
                totalTiltHighAmpLcus++;
            }
        }
    }
    picture_control_set_ptr->is_pan = EB_FALSE;
    picture_control_set_ptr->is_tilt = EB_FALSE;

    picture_control_set_ptr->panMvx = 0;
    picture_control_set_ptr->panMvy = 0;
    picture_control_set_ptr->tiltMvx = 0;
    picture_control_set_ptr->tiltMvy = 0;

    // If more than PAN_LCU_PERCENTAGE % of LCUs are PAN
    if ((totalCheckedLcus > 0) && ((totalPanLcus * 100 / totalCheckedLcus) > PAN_LCU_PERCENTAGE)) {
        picture_control_set_ptr->is_pan = EB_TRUE;
        if (totalPanLcus > 0) {
            picture_control_set_ptr->panMvx = (int16_t)(xPanMvSum / totalPanLcus);
            picture_control_set_ptr->panMvy = (int16_t)(yPanMvSum / totalPanLcus);
        }
    }

    if ((totalCheckedLcus > 0) && ((totalTiltLcus * 100 / totalCheckedLcus) > PAN_LCU_PERCENTAGE)) {
        picture_control_set_ptr->is_tilt = EB_TRUE;
        if (totalTiltLcus > 0) {
            picture_control_set_ptr->tiltMvx = (int16_t)(xTiltMvSum / totalTiltLcus);
            picture_control_set_ptr->tiltMvy = (int16_t)(yTiltMvSum / totalTiltLcus);
        }
    }
#if GM_OPT && GLOBAL_WARPED_MOTION
    }
#endif
#endif
}

static void initial_rate_control_context_dctor(EbPtr p)
{
    EbThreadContext   *thread_context_ptr = (EbThreadContext*)p;
    InitialRateControlContext* obj = (InitialRateControlContext*)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}


/************************************************
* Initial Rate Control Context Constructor
************************************************/
EbErrorType initial_rate_control_context_ctor(
    EbThreadContext     *thread_context_ptr,
    const EbEncHandle   *enc_handle_ptr)
{
    InitialRateControlContext  *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv = context_ptr;
    thread_context_ptr->dctor = initial_rate_control_context_dctor;

    context_ptr->motion_estimation_results_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->motion_estimation_results_resource_ptr, 0);
    context_ptr->initialrate_control_results_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->initial_rate_control_results_resource_ptr, 0);

    return EB_ErrorNone;
}

/************************************************
* Release Pa Reference Objects
** Check if reference pictures are needed
** release them when appropriate
************************************************/
void ReleasePaReferenceObjects(
    SequenceControlSet              *sequence_control_set_ptr,
    PictureParentControlSet         *picture_control_set_ptr)
{
    // PA Reference Pictures
    uint32_t                             numOfListToSearch;
    uint32_t                             listIndex;
    uint32_t                             ref_pic_index;
    if (picture_control_set_ptr->slice_type != I_SLICE) {
        numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE) ? REF_LIST_0 : REF_LIST_1;

        // List Loop
        for (listIndex = REF_LIST_0; listIndex <= numOfListToSearch; ++listIndex) {
            // Release PA Reference Pictures
            uint8_t num_of_ref_pic_to_search = (picture_control_set_ptr->slice_type == P_SLICE) ?
                MIN(picture_control_set_ptr->ref_list0_count, sequence_control_set_ptr->reference_count) :
                (listIndex == REF_LIST_0) ?
                MIN(picture_control_set_ptr->ref_list0_count, sequence_control_set_ptr->reference_count) :
                MIN(picture_control_set_ptr->ref_list1_count, sequence_control_set_ptr->reference_count);

            for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
                if (picture_control_set_ptr->ref_pa_pic_ptr_array[listIndex][ref_pic_index] != EB_NULL) {
                    eb_release_object(picture_control_set_ptr->ref_pa_pic_ptr_array[listIndex][ref_pic_index]);
                }
            }
        }
    }

    if (picture_control_set_ptr->pa_reference_picture_wrapper_ptr != EB_NULL) {
        eb_release_object(picture_control_set_ptr->pa_reference_picture_wrapper_ptr);
    }

    return;
}

/************************************************
* Global Motion Detection Based on ME information
** Mark pictures for pan
** Mark pictures for tilt
** No lookahead information used in this function
************************************************/
void MeBasedGlobalMotionDetection(
    PictureParentControlSet         *picture_control_set_ptr)
{
    // PAN Generation
    picture_control_set_ptr->is_pan = EB_FALSE;
    picture_control_set_ptr->is_tilt = EB_FALSE;

    if (picture_control_set_ptr->slice_type != I_SLICE)
        DetectGlobalMotion(picture_control_set_ptr);
    // Check if the motion vector field for temporal layer 0 pictures
    if (picture_control_set_ptr->slice_type != I_SLICE && picture_control_set_ptr->temporal_layer_index == 0)
        CheckForNonUniformMotionVectorField(picture_control_set_ptr);
    return;
}

void StationaryEdgeCountLcu(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    PictureParentControlSet   *temporalPictureControlSetPtr,
    uint32_t                       totalLcuCount)
{
    uint32_t               sb_index;
    for (sb_index = 0; sb_index < totalLcuCount; sb_index++) {
        SbParams sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
        SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];
        if (sb_params.potential_logo_sb &&sb_params.is_complete_sb && sb_stat_ptr->check1_for_logo_stationary_edge_over_time_flag && sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag) {
            SbStat *tempLcuStatPtr = &temporalPictureControlSetPtr->sb_stat_array[sb_index];
            uint32_t rasterScanCuIndex;

            if (tempLcuStatPtr->check1_for_logo_stationary_edge_over_time_flag)
            {
                for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++)
                    sb_stat_ptr->cu_stat_array[rasterScanCuIndex].similar_edge_count += tempLcuStatPtr->cu_stat_array[rasterScanCuIndex].edge_cu;
            }
        }

        if (sb_params.potential_logo_sb &&sb_params.is_complete_sb && sb_stat_ptr->pm_check1_for_logo_stationary_edge_over_time_flag && sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag) {
            SbStat *tempLcuStatPtr = &temporalPictureControlSetPtr->sb_stat_array[sb_index];
            uint32_t rasterScanCuIndex;

            if (tempLcuStatPtr->pm_check1_for_logo_stationary_edge_over_time_flag)
            {
                for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++)
                    sb_stat_ptr->cu_stat_array[rasterScanCuIndex].pm_similar_edge_count += tempLcuStatPtr->cu_stat_array[rasterScanCuIndex].edge_cu;
            }
        }
    }
}

void StationaryEdgeOverUpdateOverTimeLcuPart1(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr)
{
    uint32_t               sb_index;
    int32_t                 xCurrentMv = 0;
    int32_t                 yCurrentMv = 0;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; sb_index++) {
        SbParams sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
        SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

        if (sb_params.potential_logo_sb &&sb_params.is_complete_sb) {
            // Current MV
            if (picture_control_set_ptr->temporal_layer_index > 0)
                GetMv(picture_control_set_ptr, sb_index, &xCurrentMv, &yCurrentMv);

            EbBool lowMotion = picture_control_set_ptr->temporal_layer_index == 0 ? EB_TRUE : (ABS(xCurrentMv) < 16) && (ABS(yCurrentMv) < 16) ? EB_TRUE : EB_FALSE;
            uint16_t *yVariancePtr = picture_control_set_ptr->variance[sb_index];
            uint64_t var0 = yVariancePtr[ME_TIER_ZERO_PU_32x32_0];
            uint64_t var1 = yVariancePtr[ME_TIER_ZERO_PU_32x32_1];
            uint64_t var2 = yVariancePtr[ME_TIER_ZERO_PU_32x32_2];
            uint64_t var3 = yVariancePtr[ME_TIER_ZERO_PU_32x32_3];

            uint64_t averageVar = (var0 + var1 + var2 + var3) >> 2;
            uint64_t varOfVar = (((int32_t)(var0 - averageVar) * (int32_t)(var0 - averageVar)) +
                ((int32_t)(var1 - averageVar) * (int32_t)(var1 - averageVar)) +
                ((int32_t)(var2 - averageVar) * (int32_t)(var2 - averageVar)) +
                ((int32_t)(var3 - averageVar) * (int32_t)(var3 - averageVar))) >> 2;

            if ((varOfVar <= 50000) || !lowMotion)
                sb_stat_ptr->check1_for_logo_stationary_edge_over_time_flag = 0;
            else
                sb_stat_ptr->check1_for_logo_stationary_edge_over_time_flag = 1;
            if ((varOfVar <= 1000))
                sb_stat_ptr->pm_check1_for_logo_stationary_edge_over_time_flag = 0;
            else
                sb_stat_ptr->pm_check1_for_logo_stationary_edge_over_time_flag = 1;
        }
        else {
            sb_stat_ptr->check1_for_logo_stationary_edge_over_time_flag = 0;

            sb_stat_ptr->pm_check1_for_logo_stationary_edge_over_time_flag = 0;
        }
    }
}
void StationaryEdgeOverUpdateOverTimeLcuPart2(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr)
{
    uint32_t               sb_index;

    uint32_t               lowSadTh = (sequence_control_set_ptr->input_resolution < INPUT_SIZE_1080p_RANGE) ? 5 : 2;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; sb_index++) {
        SbParams sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
        SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

        if (sb_params.potential_logo_sb &&sb_params.is_complete_sb) {
            uint32_t meDist = 0;

            EbBool lowSad = EB_FALSE;

            if (picture_control_set_ptr->slice_type == B_SLICE)
                GetMeDist(picture_control_set_ptr, sb_index, &meDist);
            lowSad = (picture_control_set_ptr->slice_type != B_SLICE) ?

                EB_FALSE : (meDist < 64 * 64 * lowSadTh) ? EB_TRUE : EB_FALSE;

            if (lowSad) {
                sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag = 0;
                sb_stat_ptr->low_dist_logo = 1;
            }
            else {
                sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag = 1;

                sb_stat_ptr->low_dist_logo = 0;
            }
        }
        else {
            sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag = 0;

            sb_stat_ptr->low_dist_logo = 0;
        }
        sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag = 1;
    }
}

void StationaryEdgeOverUpdateOverTimeLcu(
    SequenceControlSet        *sequence_control_set_ptr,
    uint32_t                        totalCheckedPictures,
    PictureParentControlSet   *picture_control_set_ptr,
    uint32_t                       totalLcuCount)
{
    uint32_t               sb_index;
    const uint32_t         slideWindowTh = ((totalCheckedPictures / 4) - 1);

    for (sb_index = 0; sb_index < totalLcuCount; sb_index++) {
        SbParams sb_params = sequence_control_set_ptr->sb_params_array[sb_index];

        SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];
        sb_stat_ptr->stationary_edge_over_time_flag = EB_FALSE;
        if (sb_params.potential_logo_sb &&sb_params.is_complete_sb && sb_stat_ptr->check1_for_logo_stationary_edge_over_time_flag && sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag) {
            uint32_t rasterScanCuIndex;
            uint32_t similarEdgeCountLcu = 0;
            // CU Loop
            for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++)
                similarEdgeCountLcu += (sb_stat_ptr->cu_stat_array[rasterScanCuIndex].similar_edge_count > slideWindowTh) ? 1 : 0;
            sb_stat_ptr->stationary_edge_over_time_flag = (similarEdgeCountLcu >= 4) ? EB_TRUE : EB_FALSE;
        }

        sb_stat_ptr->pm_stationary_edge_over_time_flag = EB_FALSE;
        if (sb_params.potential_logo_sb &&sb_params.is_complete_sb && sb_stat_ptr->pm_check1_for_logo_stationary_edge_over_time_flag && sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag) {
            uint32_t rasterScanCuIndex;
            uint32_t similarEdgeCountLcu = 0;
            // CU Loop
            for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++)
                similarEdgeCountLcu += (sb_stat_ptr->cu_stat_array[rasterScanCuIndex].pm_similar_edge_count > slideWindowTh) ? 1 : 0;
            sb_stat_ptr->pm_stationary_edge_over_time_flag = (similarEdgeCountLcu >= 4) ? EB_TRUE : EB_FALSE;
        }
    }
    {
        uint32_t sb_index;
        uint32_t sb_x, sb_y;
        uint32_t countOfNeighbors = 0;

        uint32_t countOfNeighborsPm = 0;

        int32_t lcuHor, lcuVer, lcuVerOffset;
        int32_t lcuHorS, lcuVerS, lcuHorE, lcuVerE;
        uint32_t picture_width_in_sb = sequence_control_set_ptr->picture_width_in_sb;
        uint32_t picture_height_in_sb = sequence_control_set_ptr->picture_height_in_sb;

        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
            SbParams                sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
            SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            sb_x = sb_params.horizontal_index;
            sb_y = sb_params.vertical_index;
            if (sb_params.potential_logo_sb &&sb_params.is_complete_sb && sb_stat_ptr->check1_for_logo_stationary_edge_over_time_flag && sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag) {
                {
                    lcuHorS = (sb_x > 0) ? -1 : 0;
                    lcuHorE = (sb_x < picture_width_in_sb - 1) ? 1 : 0;
                    lcuVerS = (sb_y > 0) ? -1 : 0;
                    lcuVerE = (sb_y < picture_height_in_sb - 1) ? 1 : 0;
                    countOfNeighbors = 0;
                    for (lcuVer = lcuVerS; lcuVer <= lcuVerE; lcuVer++) {
                        lcuVerOffset = lcuVer * (int32_t)picture_width_in_sb;
                        for (lcuHor = lcuHorS; lcuHor <= lcuHorE; lcuHor++)
                            countOfNeighbors += (picture_control_set_ptr->sb_stat_array[sb_index + lcuVerOffset + lcuHor].stationary_edge_over_time_flag == 1);
                    }
                    if (countOfNeighbors == 1)
                        sb_stat_ptr->stationary_edge_over_time_flag = 0;
                }
            }

            if (sb_params.potential_logo_sb &&sb_params.is_complete_sb && sb_stat_ptr->pm_check1_for_logo_stationary_edge_over_time_flag && sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag) {
                {
                    lcuHorS = (sb_x > 0) ? -1 : 0;
                    lcuHorE = (sb_x < picture_width_in_sb - 1) ? 1 : 0;
                    lcuVerS = (sb_y > 0) ? -1 : 0;
                    lcuVerE = (sb_y < picture_height_in_sb - 1) ? 1 : 0;
                    countOfNeighborsPm = 0;
                    for (lcuVer = lcuVerS; lcuVer <= lcuVerE; lcuVer++) {
                        lcuVerOffset = lcuVer * (int32_t)picture_width_in_sb;
                        for (lcuHor = lcuHorS; lcuHor <= lcuHorE; lcuHor++)
                            countOfNeighborsPm += (picture_control_set_ptr->sb_stat_array[sb_index + lcuVerOffset + lcuHor].pm_stationary_edge_over_time_flag == 1);
                    }

                    if (countOfNeighborsPm == 1)
                        sb_stat_ptr->pm_stationary_edge_over_time_flag = 0;
                }
            }
        }
    }

    {
        uint32_t sb_index;
        uint32_t sb_x, sb_y;
        uint32_t countOfNeighbors = 0;
        int32_t lcuHor, lcuVer, lcuVerOffset;
        int32_t lcuHorS, lcuVerS, lcuHorE, lcuVerE;
        uint32_t picture_width_in_sb = sequence_control_set_ptr->picture_width_in_sb;
        uint32_t picture_height_in_sb = sequence_control_set_ptr->picture_height_in_sb;

        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
            SbParams                sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
            SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            sb_x = sb_params.horizontal_index;
            sb_y = sb_params.vertical_index;

            {
                if (sb_stat_ptr->stationary_edge_over_time_flag == 0 && sb_params.potential_logo_sb && (sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag || !sb_params.is_complete_sb)) {
                    lcuHorS = (sb_x > 0) ? -1 : 0;
                    lcuHorE = (sb_x < picture_width_in_sb - 1) ? 1 : 0;
                    lcuVerS = (sb_y > 0) ? -1 : 0;
                    lcuVerE = (sb_y < picture_height_in_sb - 1) ? 1 : 0;
                    countOfNeighbors = 0;
                    for (lcuVer = lcuVerS; lcuVer <= lcuVerE; lcuVer++) {
                        lcuVerOffset = lcuVer * (int32_t)picture_width_in_sb;
                        for (lcuHor = lcuHorS; lcuHor <= lcuHorE; lcuHor++)
                            countOfNeighbors += (picture_control_set_ptr->sb_stat_array[sb_index + lcuVerOffset + lcuHor].stationary_edge_over_time_flag == 1);
                    }
                    if (countOfNeighbors > 0)
                        sb_stat_ptr->stationary_edge_over_time_flag = 2;
                }
            }
        }

        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
            SbParams                sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
            SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            sb_x = sb_params.horizontal_index;
            sb_y = sb_params.vertical_index;

            {
                if (sb_stat_ptr->stationary_edge_over_time_flag == 0 && sb_params.potential_logo_sb && (sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag || !sb_params.is_complete_sb)) {
                    lcuHorS = (sb_x > 0) ? -1 : 0;
                    lcuHorE = (sb_x < picture_width_in_sb - 1) ? 1 : 0;
                    lcuVerS = (sb_y > 0) ? -1 : 0;
                    lcuVerE = (sb_y < picture_height_in_sb - 1) ? 1 : 0;
                    countOfNeighbors = 0;
                    for (lcuVer = lcuVerS; lcuVer <= lcuVerE; lcuVer++) {
                        lcuVerOffset = lcuVer * (int32_t)picture_width_in_sb;
                        for (lcuHor = lcuHorS; lcuHor <= lcuHorE; lcuHor++)
                            countOfNeighbors += (picture_control_set_ptr->sb_stat_array[sb_index + lcuVerOffset + lcuHor].stationary_edge_over_time_flag == 2);
                    }
                    if (countOfNeighbors > 3)
                        sb_stat_ptr->stationary_edge_over_time_flag = 3;
                }
            }
        }
    }

    {
        uint32_t sb_index;
        uint32_t sb_x, sb_y;
        uint32_t countOfNeighbors = 0;
        int32_t lcuHor, lcuVer, lcuVerOffset;
        int32_t lcuHorS, lcuVerS, lcuHorE, lcuVerE;
        uint32_t picture_width_in_sb = sequence_control_set_ptr->picture_width_in_sb;
        uint32_t picture_height_in_sb = sequence_control_set_ptr->picture_height_in_sb;

        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
            SbParams                sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
            SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            sb_x = sb_params.horizontal_index;
            sb_y = sb_params.vertical_index;

            {
                if (sb_stat_ptr->pm_stationary_edge_over_time_flag == 0 && sb_params.potential_logo_sb && (sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag || !sb_params.is_complete_sb)) {
                    lcuHorS = (sb_x > 0) ? -1 : 0;
                    lcuHorE = (sb_x < picture_width_in_sb - 1) ? 1 : 0;
                    lcuVerS = (sb_y > 0) ? -1 : 0;
                    lcuVerE = (sb_y < picture_height_in_sb - 1) ? 1 : 0;
                    countOfNeighbors = 0;
                    for (lcuVer = lcuVerS; lcuVer <= lcuVerE; lcuVer++) {
                        lcuVerOffset = lcuVer * (int32_t)picture_width_in_sb;
                        for (lcuHor = lcuHorS; lcuHor <= lcuHorE; lcuHor++)
                            countOfNeighbors += (picture_control_set_ptr->sb_stat_array[sb_index + lcuVerOffset + lcuHor].pm_stationary_edge_over_time_flag == 1);
                    }
                    if (countOfNeighbors > 0)
                        sb_stat_ptr->pm_stationary_edge_over_time_flag = 2;
                }
            }
        }

        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
            SbParams                sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
            SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            sb_x = sb_params.horizontal_index;
            sb_y = sb_params.vertical_index;

            {
                if (sb_stat_ptr->pm_stationary_edge_over_time_flag == 0 && sb_params.potential_logo_sb && (sb_stat_ptr->check2_for_logo_stationary_edge_over_time_flag || !sb_params.is_complete_sb)) {
                    lcuHorS = (sb_x > 0) ? -1 : 0;
                    lcuHorE = (sb_x < picture_width_in_sb - 1) ? 1 : 0;
                    lcuVerS = (sb_y > 0) ? -1 : 0;
                    lcuVerE = (sb_y < picture_height_in_sb - 1) ? 1 : 0;
                    countOfNeighbors = 0;
                    for (lcuVer = lcuVerS; lcuVer <= lcuVerE; lcuVer++) {
                        lcuVerOffset = lcuVer * (int32_t)picture_width_in_sb;
                        for (lcuHor = lcuHorS; lcuHor <= lcuHorE; lcuHor++)
                            countOfNeighbors += (picture_control_set_ptr->sb_stat_array[sb_index + lcuVerOffset + lcuHor].pm_stationary_edge_over_time_flag == 2);
                    }
                    if (countOfNeighbors > 3)
                        sb_stat_ptr->pm_stationary_edge_over_time_flag = 3;
                }
            }
        }
    }
}

/************************************************
* Global Motion Detection Based on Lookahead
** Mark pictures for pan
** Mark pictures for tilt
** LAD Window: min (8 or sliding window size)
************************************************/
void UpdateGlobalMotionDetectionOverTime(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *sequence_control_set_ptr,
    PictureParentControlSet         *picture_control_set_ptr)
{
    InitialRateControlReorderEntry   *temporaryQueueEntryPtr;
    PictureParentControlSet          *temporaryPictureControlSetPtr;

    uint32_t                                totalPanPictures = 0;
    uint32_t                                totalCheckedPictures = 0;
    uint32_t                                totalTiltPictures = 0;
    uint32_t                                updateIsPanFramesToCheck;
    uint32_t                                inputQueueIndex;
    uint32_t                                framesToCheckIndex;

    (void)sequence_control_set_ptr;

    // Determine number of frames to check (8 frames)
    updateIsPanFramesToCheck = MIN(8, picture_control_set_ptr->frames_in_sw);

    // Walk the first N entries in the sliding window
    inputQueueIndex = encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    uint32_t updateFramesToCheck = updateIsPanFramesToCheck;
    for (framesToCheckIndex = 0; framesToCheckIndex < updateFramesToCheck; framesToCheckIndex++) {
        temporaryQueueEntryPtr = encode_context_ptr->initial_rate_control_reorder_queue[inputQueueIndex];
        temporaryPictureControlSetPtr = ((PictureParentControlSet*)(temporaryQueueEntryPtr->parent_pcs_wrapper_ptr)->object_ptr);

        if (temporaryPictureControlSetPtr->slice_type != I_SLICE) {
            totalPanPictures += (temporaryPictureControlSetPtr->is_pan == EB_TRUE);

            totalTiltPictures += (temporaryPictureControlSetPtr->is_tilt == EB_TRUE);

            // Keep track of checked pictures
            totalCheckedPictures++;
        }

        // Increment the inputQueueIndex Iterator
        inputQueueIndex = (inputQueueIndex == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
    }

    picture_control_set_ptr->is_pan = EB_FALSE;
    picture_control_set_ptr->is_tilt = EB_FALSE;

    if (totalCheckedPictures) {
        if (picture_control_set_ptr->slice_type != I_SLICE) {
            if ((totalPanPictures * 100 / totalCheckedPictures) > 75)
                picture_control_set_ptr->is_pan = EB_TRUE;
        }
    }
    return;
}

/************************************************
* Update BEA Information Based on Lookahead
** Average zzCost of Collocated SB throughout lookahead frames
** Set isMostOfPictureNonMoving based on number of non moving LCUs
** LAD Window: min (2xmgpos+1 or sliding window size)
************************************************/

void UpdateBeaInfoOverTime(
    EncodeContext                   *encode_context_ptr,
    PictureParentControlSet         *picture_control_set_ptr)
{
    InitialRateControlReorderEntry   *temporaryQueueEntryPtr;
    PictureParentControlSet          *temporaryPictureControlSetPtr;
    uint32_t                                updateNonMovingIndexArrayFramesToCheck;
    uint16_t                              lcuIdx;
    uint16_t                                framesToCheckIndex;
    uint64_t                                nonMovingIndexSum = 0;
    uint32_t                                inputQueueIndex;

    SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    // Update motionIndexArray of the current picture by averaging the motionIndexArray of the N future pictures
    // Determine number of frames to check N
    updateNonMovingIndexArrayFramesToCheck = MIN(MIN(((picture_control_set_ptr->pred_struct_ptr->pred_struct_period << 1) + 1), picture_control_set_ptr->frames_in_sw), sequence_control_set_ptr->static_config.look_ahead_distance);
    uint64_t me_dist = 0;
    uint8_t me_dist_pic_count = 0;
    // SB Loop
    for (lcuIdx = 0; lcuIdx < picture_control_set_ptr->sb_total_count; ++lcuIdx) {
        uint16_t nonMovingIndexOverSlidingWindow = picture_control_set_ptr->non_moving_index_array[lcuIdx];

        // Walk the first N entries in the sliding window starting picture + 1
        inputQueueIndex = (encode_context_ptr->initial_rate_control_reorder_queue_head_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->initial_rate_control_reorder_queue_head_index + 1;
        for (framesToCheckIndex = 0; framesToCheckIndex < updateNonMovingIndexArrayFramesToCheck - 1; framesToCheckIndex++) {
            temporaryQueueEntryPtr = encode_context_ptr->initial_rate_control_reorder_queue[inputQueueIndex];
            temporaryPictureControlSetPtr = ((PictureParentControlSet*)(temporaryQueueEntryPtr->parent_pcs_wrapper_ptr)->object_ptr);

            if (temporaryPictureControlSetPtr->slice_type == I_SLICE || temporaryPictureControlSetPtr->end_of_sequence_flag)
                break;
            // Limit the distortion to lower layers 0, 1 and 2 only. Higher layers have close temporal distance and lower distortion that might contaminate the data
            if (temporaryPictureControlSetPtr->temporal_layer_index < MAX((int8_t)picture_control_set_ptr->hierarchical_levels - 1, 2) ) {
                if (lcuIdx == 0)
                    me_dist_pic_count++;
                me_dist += (temporaryPictureControlSetPtr->slice_type == I_SLICE) ? 0 : (uint64_t)temporaryPictureControlSetPtr->rc_me_distortion[lcuIdx];
            }
            // Store the filtered_sse of next ALT_REF picture in the I slice to be used in QP Scaling
            if (picture_control_set_ptr->slice_type == I_SLICE && picture_control_set_ptr->filtered_sse == 0 && lcuIdx == 0 && temporaryPictureControlSetPtr->temporal_layer_index == 0) {
                picture_control_set_ptr->filtered_sse = temporaryPictureControlSetPtr->filtered_sse;
                picture_control_set_ptr->filtered_sse_uv = temporaryPictureControlSetPtr->filtered_sse_uv;
            }
            nonMovingIndexOverSlidingWindow += temporaryPictureControlSetPtr->non_moving_index_array[lcuIdx];

            // Increment the inputQueueIndex Iterator
            inputQueueIndex = (inputQueueIndex == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
        }
        picture_control_set_ptr->non_moving_index_array[lcuIdx] = (uint8_t)(nonMovingIndexOverSlidingWindow / (framesToCheckIndex + 1));

        nonMovingIndexSum += picture_control_set_ptr->non_moving_index_array[lcuIdx];
    }

    picture_control_set_ptr->non_moving_index_average = (uint16_t)nonMovingIndexSum / picture_control_set_ptr->sb_total_count;
    me_dist_pic_count = MAX(me_dist_pic_count, 1);
    picture_control_set_ptr->qp_scaling_average_complexity = (uint16_t)((uint64_t)me_dist / picture_control_set_ptr->sb_total_count / 256 / me_dist_pic_count);
    return;
}

/****************************************
* Init ZZ Cost array to default values
** Used when no Lookahead is available
****************************************/
void InitZzCostInfo(
    PictureParentControlSet         *picture_control_set_ptr)
{
    uint16_t lcuIdx;
    picture_control_set_ptr->non_moving_index_average = INVALID_ZZ_COST;

    // SB Loop
    for (lcuIdx = 0; lcuIdx < picture_control_set_ptr->sb_total_count; ++lcuIdx)
        picture_control_set_ptr->non_moving_index_array[lcuIdx] = INVALID_ZZ_COST;
    return;
}

/************************************************
* Update uniform motion field
** Update Uniformly moving LCUs using
** collocated LCUs infor in lookahead pictures
** LAD Window: min (2xmgpos+1 or sliding window size)
************************************************/
void UpdateMotionFieldUniformityOverTime(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *sequence_control_set_ptr,
    PictureParentControlSet         *picture_control_set_ptr)
{
    InitialRateControlReorderEntry   *temporaryQueueEntryPtr;
    PictureParentControlSet          *temporaryPictureControlSetPtr;
    uint32_t                                inputQueueIndex;
    uint32_t                              NoFramesToCheck;
    uint32_t                                framesToCheckIndex;
    //SVT_LOG("To update POC %d\tframesInSw = %d\n", picture_control_set_ptr->picture_number, picture_control_set_ptr->frames_in_sw);

    //Check conditions for statinary edge over time
    StationaryEdgeOverUpdateOverTimeLcuPart2(
        sequence_control_set_ptr,
        picture_control_set_ptr);

    // Determine number of frames to check N
    NoFramesToCheck = MIN(MIN(((picture_control_set_ptr->pred_struct_ptr->pred_struct_period << 1) + 1), picture_control_set_ptr->frames_in_sw), sequence_control_set_ptr->static_config.look_ahead_distance);

    // Walk the first N entries in the sliding window starting picture + 1
    inputQueueIndex = (encode_context_ptr->initial_rate_control_reorder_queue_head_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    for (framesToCheckIndex = 0; framesToCheckIndex < NoFramesToCheck - 1; framesToCheckIndex++) {
        temporaryQueueEntryPtr = encode_context_ptr->initial_rate_control_reorder_queue[inputQueueIndex];
        temporaryPictureControlSetPtr = ((PictureParentControlSet*)(temporaryQueueEntryPtr->parent_pcs_wrapper_ptr)->object_ptr);

        if (temporaryPictureControlSetPtr->end_of_sequence_flag)
            break;
        // The values are calculated for every 4th frame
        if ((temporaryPictureControlSetPtr->picture_number & 3) == 0) {
            StationaryEdgeCountLcu(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                temporaryPictureControlSetPtr,
                picture_control_set_ptr->sb_total_count);
        }
        // Increment the inputQueueIndex Iterator
        inputQueueIndex = (inputQueueIndex == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
    }
    StationaryEdgeOverUpdateOverTimeLcu(
        sequence_control_set_ptr,
        NoFramesToCheck,
        picture_control_set_ptr,
        picture_control_set_ptr->sb_total_count);
    return;
}
InitialRateControlReorderEntry  * DeterminePictureOffsetInQueue(
    EncodeContext                   *encode_context_ptr,
    PictureParentControlSet         *picture_control_set_ptr,
    MotionEstimationResults         *inputResultsPtr)
{
    InitialRateControlReorderEntry  *queueEntryPtr;
    int32_t                             queueEntryIndex;

    queueEntryIndex = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->initial_rate_control_reorder_queue[encode_context_ptr->initial_rate_control_reorder_queue_head_index]->picture_number);
    queueEntryIndex += encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    queueEntryIndex = (queueEntryIndex > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? queueEntryIndex - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH : queueEntryIndex;
    queueEntryPtr = encode_context_ptr->initial_rate_control_reorder_queue[queueEntryIndex];
    queueEntryPtr->parent_pcs_wrapper_ptr = inputResultsPtr->picture_control_set_wrapper_ptr;
    queueEntryPtr->picture_number = picture_control_set_ptr->picture_number;

    return queueEntryPtr;
}

void GetHistogramQueueData(
    SequenceControlSet              *sequence_control_set_ptr,
    EncodeContext                   *encode_context_ptr,
    PictureParentControlSet         *picture_control_set_ptr)
{
    HlRateControlHistogramEntry     *histogramQueueEntryPtr;
    int32_t                             histogramQueueEntryIndex;

    // Determine offset from the Head Ptr for HLRC histogram queue
    eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
    histogramQueueEntryIndex = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
    histogramQueueEntryIndex += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
    histogramQueueEntryIndex = (histogramQueueEntryIndex > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
        histogramQueueEntryIndex - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
        (histogramQueueEntryIndex < 0) ?
        histogramQueueEntryIndex + HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
        histogramQueueEntryIndex;
    histogramQueueEntryPtr = encode_context_ptr->hl_rate_control_historgram_queue[histogramQueueEntryIndex];

    //histogramQueueEntryPtr->parent_pcs_wrapper_ptr  = inputResultsPtr->picture_control_set_wrapper_ptr;
    histogramQueueEntryPtr->picture_number = picture_control_set_ptr->picture_number;
    histogramQueueEntryPtr->end_of_sequence_flag = picture_control_set_ptr->end_of_sequence_flag;
    histogramQueueEntryPtr->slice_type = picture_control_set_ptr->slice_type;
    histogramQueueEntryPtr->temporal_layer_index = picture_control_set_ptr->temporal_layer_index;
    histogramQueueEntryPtr->full_sb_count = picture_control_set_ptr->full_sb_count;
    histogramQueueEntryPtr->life_count = 0;
    histogramQueueEntryPtr->passed_to_hlrc = EB_FALSE;
    histogramQueueEntryPtr->is_coded = EB_FALSE;
    histogramQueueEntryPtr->total_num_bits_coded = 0;
    histogramQueueEntryPtr->frames_in_sw = 0;
    EB_MEMCPY(
        histogramQueueEntryPtr->me_distortion_histogram,
        picture_control_set_ptr->me_distortion_histogram,
        sizeof(uint16_t) * NUMBER_OF_SAD_INTERVALS);

    EB_MEMCPY(
        histogramQueueEntryPtr->ois_distortion_histogram,
        picture_control_set_ptr->ois_distortion_histogram,
        sizeof(uint16_t) * NUMBER_OF_INTRA_SAD_INTERVALS);

    eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
    //SVT_LOG("Test1 POC: %d\t POC: %d\t LifeCount: %d\n", histogramQueueEntryPtr->picture_number, picture_control_set_ptr->picture_number,  histogramQueueEntryPtr->life_count);

    return;
}

void UpdateHistogramQueueEntry(
    SequenceControlSet              *sequence_control_set_ptr,
    EncodeContext                   *encode_context_ptr,
    PictureParentControlSet         *picture_control_set_ptr,
    uint32_t                           frames_in_sw)
{
    HlRateControlHistogramEntry     *histogramQueueEntryPtr;
    int32_t                             histogramQueueEntryIndex;

    eb_block_on_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);

    histogramQueueEntryIndex = (int32_t)(picture_control_set_ptr->picture_number - encode_context_ptr->hl_rate_control_historgram_queue[encode_context_ptr->hl_rate_control_historgram_queue_head_index]->picture_number);
    histogramQueueEntryIndex += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
    histogramQueueEntryIndex = (histogramQueueEntryIndex > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1) ?
        histogramQueueEntryIndex - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
        (histogramQueueEntryIndex < 0) ?
        histogramQueueEntryIndex + HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH :
        histogramQueueEntryIndex;
    histogramQueueEntryPtr = encode_context_ptr->hl_rate_control_historgram_queue[histogramQueueEntryIndex];
    histogramQueueEntryPtr->passed_to_hlrc = EB_TRUE;
    if (sequence_control_set_ptr->static_config.rate_control_mode == 2)
        histogramQueueEntryPtr->life_count += (int16_t)(sequence_control_set_ptr->static_config.intra_period_length + 1) - 3; // FramelevelRC does not decrease the life count for first picture in each temporal layer
    else
        histogramQueueEntryPtr->life_count += picture_control_set_ptr->historgram_life_count;
    histogramQueueEntryPtr->frames_in_sw = frames_in_sw;
    eb_release_mutex(sequence_control_set_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);

    return;
}


/************************************************
* Initial Rate Control Kernel
* The Initial Rate Control Process determines the initial bit budget for each
* picture depending on the data gathered in the Picture Analysis and Motion
* Analysis processes as well as the settings determined in the Picture Decision process.
* The Initial Rate Control process also employs a sliding window buffer to analyze multiple
* pictures if the delay is allowed.  Note that through this process, until the subsequent
* Picture Manager process, no reference picture data has been used.
* P.S. Temporal noise reduction is now performed in Initial Rate Control Process.
* In future we might decide to move it to Motion Analysis Process.
************************************************/
void* initial_rate_control_kernel(void *input_ptr)
{
    EbThreadContext                 *thread_context_ptr = (EbThreadContext*)input_ptr;
    InitialRateControlContext       *context_ptr = (InitialRateControlContext*)thread_context_ptr->priv;
    PictureParentControlSet         *picture_control_set_ptr;
    PictureParentControlSet         *pictureControlSetPtrTemp;
    EncodeContext                   *encode_context_ptr;
    SequenceControlSet              *sequence_control_set_ptr;

    EbObjectWrapper                 *inputResultsWrapperPtr;
    MotionEstimationResults         *inputResultsPtr;

    EbObjectWrapper                 *outputResultsWrapperPtr;
    InitialRateControlResults       *outputResultsPtr;

    // Queue variables
    uint32_t                             queueEntryIndexTemp;
    uint32_t                             queueEntryIndexTemp2;
    InitialRateControlReorderEntry  *queueEntryPtr;

    EbBool                            moveSlideWondowFlag = EB_TRUE;
    EbBool                            end_of_sequence_flag = EB_TRUE;
    uint8_t                               frames_in_sw;
    uint8_t                               temporal_layer_index;
    EbObjectWrapper                  *reference_picture_wrapper_ptr;

    // Segments
    uint32_t                              segment_index;
    for (;;) {
        // Get Input Full Object
        eb_get_full_object(
            context_ptr->motion_estimation_results_input_fifo_ptr,
            &inputResultsWrapperPtr);

        inputResultsPtr = (MotionEstimationResults*)inputResultsWrapperPtr->object_ptr;
        picture_control_set_ptr = (PictureParentControlSet*)inputResultsPtr->picture_control_set_wrapper_ptr->object_ptr;

        segment_index = inputResultsPtr->segment_index;

        // Set the segment mask
        SEGMENT_COMPLETION_MASK_SET(picture_control_set_ptr->me_segments_completion_mask, segment_index);

        // If the picture is complete, proceed
        if (SEGMENT_COMPLETION_MASK_TEST(picture_control_set_ptr->me_segments_completion_mask, picture_control_set_ptr->me_segments_total_count)) {
            sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
            encode_context_ptr = (EncodeContext*)sequence_control_set_ptr->encode_context_ptr;
            // Mark picture when global motion is detected using ME results
            //reset intraCodedEstimationLcu
            MeBasedGlobalMotionDetection(
                picture_control_set_ptr);
            // Release Pa Ref pictures when not needed
            ReleasePaReferenceObjects(
                sequence_control_set_ptr,
                picture_control_set_ptr);

            //****************************************************
            // Input Motion Analysis Results into Reordering Queue
            //****************************************************

            if(!picture_control_set_ptr->is_overlay)
            // Determine offset from the Head Ptr
            queueEntryPtr = DeterminePictureOffsetInQueue(
                encode_context_ptr,
                picture_control_set_ptr,
                inputResultsPtr);

            if (sequence_control_set_ptr->static_config.rate_control_mode)
            {
                if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                    // Getting the Histogram Queue Data
                    GetHistogramQueueData(
                        sequence_control_set_ptr,
                        encode_context_ptr,
                        picture_control_set_ptr);
                }
            }

            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS; temporal_layer_index++)
                picture_control_set_ptr->frames_in_interval[temporal_layer_index] = 0;
            picture_control_set_ptr->frames_in_sw = 0;
            picture_control_set_ptr->historgram_life_count = 0;
            picture_control_set_ptr->scene_change_in_gop = EB_FALSE;

            //Check conditions for statinary edge over time

            StationaryEdgeOverUpdateOverTimeLcuPart1(
                sequence_control_set_ptr,
                picture_control_set_ptr);

            moveSlideWondowFlag = EB_TRUE;
            while (moveSlideWondowFlag) {
                // Check if the sliding window condition is valid
                queueEntryIndexTemp = encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                if (encode_context_ptr->initial_rate_control_reorder_queue[queueEntryIndexTemp]->parent_pcs_wrapper_ptr != EB_NULL)
                    end_of_sequence_flag = (((PictureParentControlSet*)(encode_context_ptr->initial_rate_control_reorder_queue[queueEntryIndexTemp]->parent_pcs_wrapper_ptr)->object_ptr))->end_of_sequence_flag;
                else
                    end_of_sequence_flag = EB_FALSE;
                frames_in_sw = 0;
                while (moveSlideWondowFlag && !end_of_sequence_flag &&
                    queueEntryIndexTemp <= encode_context_ptr->initial_rate_control_reorder_queue_head_index + sequence_control_set_ptr->static_config.look_ahead_distance) {
                    // frames_in_sw <= sequence_control_set_ptr->static_config.look_ahead_distance){
                    frames_in_sw++;

                    queueEntryIndexTemp2 = (queueEntryIndexTemp > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH : queueEntryIndexTemp;

                    moveSlideWondowFlag = (EbBool)(moveSlideWondowFlag && (encode_context_ptr->initial_rate_control_reorder_queue[queueEntryIndexTemp2]->parent_pcs_wrapper_ptr != EB_NULL));
                    if (encode_context_ptr->initial_rate_control_reorder_queue[queueEntryIndexTemp2]->parent_pcs_wrapper_ptr != EB_NULL) {
                        // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                        end_of_sequence_flag = ((PictureParentControlSet*)(encode_context_ptr->initial_rate_control_reorder_queue[queueEntryIndexTemp2]->parent_pcs_wrapper_ptr)->object_ptr)->end_of_sequence_flag;
                    }
                    else
                        end_of_sequence_flag = EB_FALSE;
                    queueEntryIndexTemp++;
                }

                if (moveSlideWondowFlag) {
                    //get a new entry spot
                    queueEntryPtr = encode_context_ptr->initial_rate_control_reorder_queue[encode_context_ptr->initial_rate_control_reorder_queue_head_index];
                    picture_control_set_ptr = ((PictureParentControlSet*)(queueEntryPtr->parent_pcs_wrapper_ptr)->object_ptr);
                    sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
                    // overlay picture was not added to the queue. For the alt_ref picture with an overlay picture, it loops on both alt ref and overlay pictures
                    uint8_t has_overlay = picture_control_set_ptr->is_alt_ref ? 1 : 0;
                    for (uint8_t loop_index = 0; loop_index <= has_overlay; loop_index++) {
                        if (loop_index)
                            picture_control_set_ptr = picture_control_set_ptr->overlay_ppcs_ptr;
                        picture_control_set_ptr->frames_in_sw = frames_in_sw;
                        queueEntryIndexTemp = encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                        end_of_sequence_flag = EB_FALSE;
                        // find the frames_in_interval for the peroid I frames
                        while (!end_of_sequence_flag &&
                            queueEntryIndexTemp <= encode_context_ptr->initial_rate_control_reorder_queue_head_index + sequence_control_set_ptr->static_config.look_ahead_distance) {
                            queueEntryIndexTemp2 = (queueEntryIndexTemp > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? queueEntryIndexTemp - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH : queueEntryIndexTemp;
                            pictureControlSetPtrTemp = ((PictureParentControlSet*)(encode_context_ptr->initial_rate_control_reorder_queue[queueEntryIndexTemp2]->parent_pcs_wrapper_ptr)->object_ptr);
                            if (sequence_control_set_ptr->intra_period_length != -1) {
                                if (picture_control_set_ptr->picture_number % ((sequence_control_set_ptr->intra_period_length + 1)) == 0) {
                                    picture_control_set_ptr->frames_in_interval[pictureControlSetPtrTemp->temporal_layer_index] ++;
                                    if (pictureControlSetPtrTemp->scene_change_flag)
                                        picture_control_set_ptr->scene_change_in_gop = EB_TRUE;
                                }
                            }

                            pictureControlSetPtrTemp->historgram_life_count++;
                            end_of_sequence_flag = pictureControlSetPtrTemp->end_of_sequence_flag;
                            queueEntryIndexTemp++;
                        }

                        if ((sequence_control_set_ptr->static_config.look_ahead_distance != 0) && (frames_in_sw < (sequence_control_set_ptr->static_config.look_ahead_distance + 1)))
                            picture_control_set_ptr->end_of_sequence_region = EB_TRUE;
                        else
                            picture_control_set_ptr->end_of_sequence_region = EB_FALSE;

                        if (sequence_control_set_ptr->static_config.rate_control_mode)
                        {
                            // Determine offset from the Head Ptr for HLRC histogram queue and set the life count
                            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                                // Update Histogram Queue Entry Life count
                                UpdateHistogramQueueEntry(
                                    sequence_control_set_ptr,
                                    encode_context_ptr,
                                    picture_control_set_ptr,
                                    frames_in_sw);
                            }
                        }

                        // Mark each input picture as PAN or not
                        // If a lookahead is present then check PAN for a period of time
                        if (!picture_control_set_ptr->end_of_sequence_flag && sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                            // Check for Pan,Tilt, Zoom and other global motion detectors over the future pictures in the lookahead
                            UpdateGlobalMotionDetectionOverTime(
                                encode_context_ptr,
                                sequence_control_set_ptr,
                                picture_control_set_ptr);
                        }
                        else {
                            if (picture_control_set_ptr->slice_type != I_SLICE)
                                DetectGlobalMotion(picture_control_set_ptr);
                        }

                        // BACKGROUND ENHANCEMENT PART II
                        if (!picture_control_set_ptr->end_of_sequence_flag && sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                            // Update BEA information based on Lookahead information
                            UpdateBeaInfoOverTime(
                                encode_context_ptr,
                                picture_control_set_ptr);
                        }
                        else {
                            // Reset zzCost information to default When there's no lookahead available
                            InitZzCostInfo(
                                picture_control_set_ptr);
                        }

                        // Use the temporal layer 0 isLcuMotionFieldNonUniform array for all the other layer pictures in the mini GOP
                        if (!picture_control_set_ptr->end_of_sequence_flag && sequence_control_set_ptr->static_config.look_ahead_distance != 0) {
                            // Updat uniformly moving LCUs based on Collocated LCUs in LookAhead window
                            UpdateMotionFieldUniformityOverTime(
                                encode_context_ptr,
                                sequence_control_set_ptr,
                                picture_control_set_ptr);
                        }
                        // Get Empty Reference Picture Object
                        eb_get_empty_object(
                            sequence_control_set_ptr->encode_context_ptr->reference_picture_pool_fifo_ptr,
                            &reference_picture_wrapper_ptr);
                        if (loop_index) {
                            picture_control_set_ptr->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                            // Give the new Reference a nominal live_count of 1
                            eb_object_inc_live_count(
                                picture_control_set_ptr->reference_picture_wrapper_ptr,
                                1);
                        }
                        else {
                            ((PictureParentControlSet*)(queueEntryPtr->parent_pcs_wrapper_ptr->object_ptr))->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                            // Give the new Reference a nominal live_count of 1
                            eb_object_inc_live_count(
                                ((PictureParentControlSet*)(queueEntryPtr->parent_pcs_wrapper_ptr->object_ptr))->reference_picture_wrapper_ptr,
                                1);
                        }
#if TWO_PASS
                        picture_control_set_ptr->stat_struct_first_pass_ptr = picture_control_set_ptr->is_used_as_reference_flag ? &((EbReferenceObject*)picture_control_set_ptr->reference_picture_wrapper_ptr->object_ptr)->stat_struct : &picture_control_set_ptr->stat_struct;
                        if (sequence_control_set_ptr->use_output_stat_file)
                            memset(picture_control_set_ptr->stat_struct_first_pass_ptr, 0, sizeof(stat_struct_t));
#endif
                        // Get Empty Results Object
                        eb_get_empty_object(
                            context_ptr->initialrate_control_results_output_fifo_ptr,
                            &outputResultsWrapperPtr);

                        outputResultsPtr = (InitialRateControlResults*)outputResultsWrapperPtr->object_ptr;

                        if (loop_index)
                            outputResultsPtr->picture_control_set_wrapper_ptr = picture_control_set_ptr->p_pcs_wrapper_ptr;
                        else
                        outputResultsPtr->picture_control_set_wrapper_ptr = queueEntryPtr->parent_pcs_wrapper_ptr;
                        // Post the Full Results Object
                        eb_post_full_object(outputResultsWrapperPtr);
                    }
                    // Reset the Reorder Queue Entry
                    queueEntryPtr->picture_number += INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
                    queueEntryPtr->parent_pcs_wrapper_ptr = (EbObjectWrapper *)EB_NULL;

                    // Increment the Reorder Queue head Ptr
                    encode_context_ptr->initial_rate_control_reorder_queue_head_index =
                        (encode_context_ptr->initial_rate_control_reorder_queue_head_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : encode_context_ptr->initial_rate_control_reorder_queue_head_index + 1;

                    queueEntryPtr = encode_context_ptr->initial_rate_control_reorder_queue[encode_context_ptr->initial_rate_control_reorder_queue_head_index];
                }
            }
        }

        // Release the Input Results
        eb_release_object(inputResultsWrapperPtr);
    }
    return EB_NULL;
}
// clang-format on
