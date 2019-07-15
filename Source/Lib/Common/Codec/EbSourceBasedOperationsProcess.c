/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"

#include "EbSourceBasedOperationsProcess.h"
#include "EbInitialRateControlResults.h"
#include "EbPictureDemuxResults.h"
#include "EbMotionEstimationContext.h"
#include "emmintrin.h"

/**************************************
* Macros
**************************************/
#define VARIANCE_PRECISION          16
#define PAN_LCU_PERCENTAGE 75
#define LOW_AMPLITUDE_TH   16

#define Y_MEAN_RANGE_02                  70
#define Y_MEAN_RANGE_01                  130
#define CB_MEAN_RANGE_02                 115
#define CR_MEAN_RANGE_00                 110
#define CR_MEAN_RANGE_02                 135

#define DARK_FRM_TH                      45
#define CB_MEAN_RANGE_00                 80

#define SAD_DEVIATION_LCU_TH             15
#define SAD_DEVIATION_LCU_NON_M4_TH      20

#define MAX_DELTA_QP_SHAPE_TH            4
#define MIN_DELTA_QP_SHAPE_TH            1

#define MIN_BLACK_AREA_PERCENTAGE        20
#define LOW_MEAN_THLD                    25

#define MIN_WHITE_AREA_PERCENTAGE         1
#define LOW_MEAN_THLD1                   40
#define HIGH_MEAN_THLD1                  210
#define NORM_FACTOR                      10 // Used ComplexityClassifier32x32
const int8_t  MinDeltaQPdefault[3] = {
    -4, -3, -2
};
const uint8_t MaxDeltaQPdefault[3] = {
    4, 5, 6
};

/************************************************
* Initial Rate Control Context Constructor
************************************************/
EbErrorType source_based_operations_context_ctor(
    SourceBasedOperationsContext  *context_ptr,
    EbFifo                        *initialRateControlResultsInputFifoPtr,
    EbFifo                        *picture_demux_results_output_fifo_ptr,
    SequenceControlSet            *sequence_control_set_ptr)
{
    UNUSED(sequence_control_set_ptr);

    context_ptr->initial_rate_control_results_input_fifo_ptr = initialRateControlResultsInputFifoPtr;
    context_ptr->picture_demux_results_output_fifo_ptr       = picture_demux_results_output_fifo_ptr;
    return EB_ErrorNone;
}

/***************************************************
* Derives BEA statistics and set activity flags
***************************************************/
void DerivePictureActivityStatistics(
    SequenceControlSet            *sequence_control_set_ptr,
    PictureParentControlSet       *picture_control_set_ptr)

{
    uint64_t               nonMovingIndexMin = ~0u;
    uint64_t               nonMovingIndexMax = 0;
    uint64_t               nonMovingIndexSum = 0;
    uint32_t               complete_sb_count = 0;
    uint32_t               non_moving_sb_count = 0;
    uint32_t               sb_total_count = picture_control_set_ptr->sb_total_count;
    uint32_t                 totNmvIdx = 0;

    uint32_t               sb_index;
    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
        SbParams *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
        if (sb_params->is_complete_sb)
        {
            nonMovingIndexMin = picture_control_set_ptr->non_moving_index_array[sb_index] < nonMovingIndexMin ?
                picture_control_set_ptr->non_moving_index_array[sb_index] :
                nonMovingIndexMin;

            nonMovingIndexMax = picture_control_set_ptr->non_moving_index_array[sb_index] > nonMovingIndexMax ?
                picture_control_set_ptr->non_moving_index_array[sb_index] :
                nonMovingIndexMax;
            if (picture_control_set_ptr->non_moving_index_array[sb_index] < NON_MOVING_SCORE_1)
                non_moving_sb_count++;
            complete_sb_count++;

            nonMovingIndexSum += picture_control_set_ptr->non_moving_index_array[sb_index];

            if (picture_control_set_ptr->non_moving_index_array[sb_index] < NON_MOVING_SCORE_1)
                totNmvIdx++;
        }
    }

    if (complete_sb_count > 0) {
        picture_control_set_ptr->non_moving_index_average = (uint16_t)(nonMovingIndexSum / complete_sb_count);
        picture_control_set_ptr->kf_zeromotion_pct = (non_moving_sb_count * 100) / complete_sb_count;
    }

    return;
}

void GrassLcu(
    SourceBasedOperationsContext        *context_ptr,
    SequenceControlSet                *sequence_control_set_ptr,
    PictureParentControlSet            *picture_control_set_ptr,
    uint32_t                                 sb_index) {
    uint32_t                  childIndex;

    EbBool                 lcuGrassFlag = EB_FALSE;

    uint32_t grassLcuInrange;
    uint32_t processedCus;
    uint32_t  rasterScanCuIndex;

    SbParams *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    SbStat *sb_stat_ptr = &(picture_control_set_ptr->sb_stat_array[sb_index]);

    _mm_prefetch((const char*)sb_stat_ptr, _MM_HINT_T0);

    lcuGrassFlag = EB_FALSE;
    grassLcuInrange = 0;
    processedCus = 0;

    for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++) {
        if (sb_params->raster_scan_cu_validity[rasterScanCuIndex]) {
            const uint32_t mdScanCuIndex = raster_scan_to_md_scan[rasterScanCuIndex];
            const uint32_t rasterScanParentCuIndex = raster_scan_cu_parent_index[rasterScanCuIndex];
            const uint32_t mdScanParentCuIndex = raster_scan_to_md_scan[rasterScanParentCuIndex];
            CuStat *cuStatPtr = &(sb_stat_ptr->cu_stat_array[mdScanCuIndex]);

            const uint32_t perfectCondition = 7;
            const uint8_t y_mean = context_ptr->y_mean_ptr[rasterScanCuIndex];
            const uint8_t cbMean = context_ptr->cb_mean_ptr[rasterScanCuIndex];
            const uint8_t crMean = context_ptr->cr_mean_ptr[rasterScanCuIndex];
            uint32_t grassCondition = 0;
            uint32_t skinCondition = 0;

            // GRASS
            grassCondition += (y_mean > Y_MEAN_RANGE_02 && y_mean < Y_MEAN_RANGE_01) ? 1 : 0;
            grassCondition += (cbMean > CB_MEAN_RANGE_00 && cbMean < CB_MEAN_RANGE_02) ? 2 : 0;
            grassCondition += (crMean > CR_MEAN_RANGE_00 && crMean < CR_MEAN_RANGE_02) ? 4 : 0;

            grassLcuInrange += (grassCondition == perfectCondition) ? 1 : 0;
            processedCus++;

            lcuGrassFlag = grassCondition == perfectCondition ? EB_TRUE : lcuGrassFlag;

            cuStatPtr->grass_area = (EbBool)(grassCondition == perfectCondition);
            // SKIN
            skinCondition += (y_mean > Y_MEAN_RANGE_02 && y_mean < Y_MEAN_RANGE_01) ? 1 : 0;
            skinCondition += (cbMean > 100 && cbMean < 120) ? 2 : 0;
            skinCondition += (crMean > 135 && crMean < 160) ? 4 : 0;

            cuStatPtr->skin_area = (EbBool)(skinCondition == perfectCondition);
            for (childIndex = mdScanCuIndex + 1; childIndex < mdScanCuIndex + 5; ++childIndex) {
                sb_stat_ptr->cu_stat_array[childIndex].grass_area = cuStatPtr->grass_area;
                sb_stat_ptr->cu_stat_array[childIndex].skin_area = cuStatPtr->skin_area;
            }

            sb_stat_ptr->cu_stat_array[mdScanParentCuIndex].grass_area = cuStatPtr->grass_area ? EB_TRUE :
                sb_stat_ptr->cu_stat_array[mdScanParentCuIndex].grass_area;
            sb_stat_ptr->cu_stat_array[0].grass_area = cuStatPtr->grass_area ? EB_TRUE :
                sb_stat_ptr->cu_stat_array[0].grass_area;
            sb_stat_ptr->cu_stat_array[mdScanParentCuIndex].skin_area = cuStatPtr->skin_area ? EB_TRUE :
                sb_stat_ptr->cu_stat_array[mdScanParentCuIndex].skin_area;
            sb_stat_ptr->cu_stat_array[0].skin_area = cuStatPtr->skin_area ? EB_TRUE :
                sb_stat_ptr->cu_stat_array[0].skin_area;
        }
    }

    context_ptr->picture_num_grass_sb += lcuGrassFlag ? 1 : 0;
}

void GrassSkinPicture(
    SourceBasedOperationsContext        *context_ptr,
    PictureParentControlSet            *picture_control_set_ptr) {
    picture_control_set_ptr->grass_percentage_in_picture = (uint8_t)(context_ptr->picture_num_grass_sb * 100 / picture_control_set_ptr->sb_total_count);
}

void SetDefaultDeltaQpRange(
    SourceBasedOperationsContext    *context_ptr,
    PictureParentControlSet        *picture_control_set_ptr,
    EbBool                             scene_transition_flag) {
    int8_t    min_delta_qp;
    uint8_t    max_delta_qp;
    if (picture_control_set_ptr->temporal_layer_index == 0) {
        min_delta_qp = MinDeltaQPdefault[0];
        max_delta_qp = MaxDeltaQPdefault[0];
    }
    else if (picture_control_set_ptr->is_used_as_reference_flag) {
        min_delta_qp = MinDeltaQPdefault[1];
        max_delta_qp = MaxDeltaQPdefault[1];
    }
    else {
        min_delta_qp = MinDeltaQPdefault[2];
        max_delta_qp = MaxDeltaQPdefault[2];
    }

    // Shape the min degrade
    min_delta_qp = (((int8_t)(min_delta_qp + MIN_DELTA_QP_SHAPE_TH) > 0) ? 0 : (min_delta_qp + MIN_DELTA_QP_SHAPE_TH));

    // Shape the max degrade
    max_delta_qp = (((int8_t)(max_delta_qp - MAX_DELTA_QP_SHAPE_TH) < 0) ? 0 : (max_delta_qp - MAX_DELTA_QP_SHAPE_TH));

    // Check on Scene Transition Flag
    max_delta_qp = scene_transition_flag ? 0 : max_delta_qp;

    context_ptr->min_delta_qp = min_delta_qp;
    context_ptr->max_delta_qp = max_delta_qp;
}

void DetermineMorePotentialAuraAreas(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet    *picture_control_set_ptr)
{
    uint16_t sb_index;
    int32_t lcuHor, lcuVer, lcuVerOffset;
    uint8_t  sb_x, sb_y;
    uint32_t countOfEdgeBlocks = 0, countOfNonEdgeBlocks = 0;

    uint32_t lightLumaValue = 150;

    uint16_t sb_total_count = picture_control_set_ptr->sb_total_count;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        SbParams *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

        sb_x = sb_params->horizontal_index;
        sb_y = sb_params->vertical_index;
        // For all the internal LCUs
        if ((sb_x > 0) && (sb_x < sequence_control_set_ptr->picture_width_in_sb - 1) && (sb_y > 0) && (sb_y < sequence_control_set_ptr->picture_height_in_sb - 1)) {
            countOfNonEdgeBlocks = 0;
            if (picture_control_set_ptr->edge_results_ptr[sb_index].edge_block_num
                && picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_64x64] >= lightLumaValue) {
                for (lcuVer = -1; lcuVer <= 1; lcuVer++) {
                    lcuVerOffset = lcuVer * (int32_t)sequence_control_set_ptr->picture_width_in_sb;
                    for (lcuHor = -1; lcuHor <= 1; lcuHor++) {
                        countOfNonEdgeBlocks += (!picture_control_set_ptr->edge_results_ptr[sb_index + lcuVerOffset + lcuHor].edge_block_num) &&
                            (picture_control_set_ptr->non_moving_index_array[sb_index + lcuVerOffset + lcuHor] < 30);
                    }
                }
            }

            if (countOfNonEdgeBlocks > 1)
                countOfEdgeBlocks++;
        }
    }

    // To check the percentage of potential aura in the picture.. If a large area is detected then this is not isolated
    picture_control_set_ptr->percentage_of_edgein_light_background = (uint8_t)(countOfEdgeBlocks * 100 / sb_total_count);
}

/***************************************************
* Detects the presence of dark area
***************************************************/
void DeriveHighDarkAreaDensityFlag(
    SequenceControlSet                *sequence_control_set_ptr,
    PictureParentControlSet           *picture_control_set_ptr) {
    uint32_t    regionInPictureWidthIndex;
    uint32_t    regionInPictureHeightIndex;
    uint32_t    lumaHistogramBin;
    uint32_t    blackSamplesCount = 0;
    uint32_t    blackAreaPercentage;
    // Loop over regions inside the picture
    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_width; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_height; regionInPictureHeightIndex++) { // loop over vertical regions
            for (lumaHistogramBin = 0; lumaHistogramBin < LOW_MEAN_THLD; lumaHistogramBin++) { // loop over the 1st LOW_MEAN_THLD bins
                blackSamplesCount += picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][lumaHistogramBin];
            }
        }
    }

    blackAreaPercentage = (blackSamplesCount * 100) / (sequence_control_set_ptr->seq_header.max_frame_width * sequence_control_set_ptr->seq_header.max_frame_height);
    picture_control_set_ptr->high_dark_area_density_flag = (EbBool)(blackAreaPercentage >= MIN_BLACK_AREA_PERCENTAGE);

    blackSamplesCount = 0;
    uint32_t    whiteSamplesCount = 0;
    uint32_t    whiteAreaPercentage;
    // Loop over regions inside the picture
    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_width; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_height; regionInPictureHeightIndex++) { // loop over vertical regions
            for (lumaHistogramBin = 0; lumaHistogramBin < LOW_MEAN_THLD1; lumaHistogramBin++) { // loop over the 1st LOW_MEAN_THLD bins
                blackSamplesCount += picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][lumaHistogramBin];
            }
            for (lumaHistogramBin = HIGH_MEAN_THLD1; lumaHistogramBin < HISTOGRAM_NUMBER_OF_BINS; lumaHistogramBin++)
                whiteSamplesCount += picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][lumaHistogramBin];
        }
    }

    blackAreaPercentage = (blackSamplesCount * 100) / (sequence_control_set_ptr->seq_header.max_frame_width * sequence_control_set_ptr->seq_header.max_frame_height);
    whiteAreaPercentage = (whiteSamplesCount * 100) / (sequence_control_set_ptr->seq_header.max_frame_width * sequence_control_set_ptr->seq_header.max_frame_height);
    picture_control_set_ptr->high_dark_low_light_area_density_flag = (EbBool)(blackAreaPercentage >= MIN_BLACK_AREA_PERCENTAGE) && (whiteAreaPercentage >= MIN_WHITE_AREA_PERCENTAGE);
}
#define NORM_FACTOR    10
#define VAR_MIN        10
#define VAR_MAX        300
#define MIN_Y        70
#define MAX_Y        145
#define MID_CB        140
#define MID_CR        115
#define TH_CB        10
#define TH_CR        15
void SpatialHighContrastClassifier(
    SourceBasedOperationsContext    *context_ptr,
    PictureParentControlSet       *picture_control_set_ptr,
    uint32_t                           sb_index)
{
    uint32_t                 blkIt;

    context_ptr->high_contrast_num = 0;
    context_ptr->high_contrast_num_ii = 0;
    //16x16 blocks
    for (blkIt = 0; blkIt < 16; blkIt++) {
        uint8_t y_mean = context_ptr->y_mean_ptr[5 + blkIt];
        uint8_t umean = context_ptr->cb_mean_ptr[5 + blkIt];
        uint8_t vmean = context_ptr->cr_mean_ptr[5 + blkIt];

        uint16_t var = picture_control_set_ptr->variance[sb_index][5 + blkIt];

        if (var > VAR_MIN && var<VAR_MAX            &&  //medium texture
            y_mean>MIN_Y && y_mean < MAX_Y            &&  //medium brightness(not too dark and not too bright)
            ABS((int64_t)umean - MID_CB) < TH_CB &&  //middle of the color plane
            ABS((int64_t)vmean - MID_CR) < TH_CR     //middle of the color plane
            )
        {
            context_ptr->high_contrast_num++;
        }

        if (
            y_mean < 30 &&  //medium brightness(not too dark and not too bright)
            ABS((int64_t)umean - 128) < 5 &&  //middle of the color plane
            ABS((int64_t)vmean - 128) < 5     //middle of the color plane
            )
        {
            context_ptr->high_contrast_num_ii++;
        }
    }
}
/************************************************
 * Source Based Operations Kernel
 ************************************************/
void* source_based_operations_kernel(void *input_ptr)
{
    SourceBasedOperationsContext    *context_ptr = (SourceBasedOperationsContext*)input_ptr;
    PictureParentControlSet       *picture_control_set_ptr;
    SequenceControlSet            *sequence_control_set_ptr;
    EbObjectWrapper               *inputResultsWrapperPtr;
    InitialRateControlResults        *inputResultsPtr;
    EbObjectWrapper               *outputResultsWrapperPtr;
    PictureDemuxResults           *outputResultsPtr;

    for (;;) {
        // Get Input Full Object
        eb_get_full_object(
            context_ptr->initial_rate_control_results_input_fifo_ptr,
            &inputResultsWrapperPtr);

        inputResultsPtr = (InitialRateControlResults*)inputResultsWrapperPtr->object_ptr;
        picture_control_set_ptr = (PictureParentControlSet*)inputResultsPtr->picture_control_set_wrapper_ptr->object_ptr;
        sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;

        picture_control_set_ptr->dark_back_groundlight_fore_ground = EB_FALSE;
        context_ptr->picture_num_grass_sb = 0;

        context_ptr->sb_cmplx_contrast_count = 0;
        context_ptr->sb_high_contrast_count = 0;
        context_ptr->complete_sb_count = 0;
        uint32_t sb_total_count = picture_control_set_ptr->sb_total_count;
        uint32_t sb_index;

        /***********************************************LCU-based operations************************************************************/
        for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
            SbParams *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            EbBool is_complete_sb = sb_params->is_complete_sb;
            uint8_t  *y_mean_ptr = picture_control_set_ptr->y_mean[sb_index];

            _mm_prefetch((const char*)y_mean_ptr, _MM_HINT_T0);
            uint8_t  *cr_mean_ptr = picture_control_set_ptr->crMean[sb_index];
            uint8_t  *cb_mean_ptr = picture_control_set_ptr->cbMean[sb_index];

            _mm_prefetch((const char*)cr_mean_ptr, _MM_HINT_T0);
            _mm_prefetch((const char*)cb_mean_ptr, _MM_HINT_T0);

            context_ptr->y_mean_ptr = y_mean_ptr;
            context_ptr->cr_mean_ptr = cr_mean_ptr;
            context_ptr->cb_mean_ptr = cb_mean_ptr;

            // Grass detection
            GrassLcu(
                context_ptr,
                sequence_control_set_ptr,
                picture_control_set_ptr,
                sb_index);

            // Spatial high contrast classifier
            if (is_complete_sb) {
                SpatialHighContrastClassifier(
                    context_ptr,
                    picture_control_set_ptr,
                    sb_index);
            }

            if (is_complete_sb) {
                context_ptr->complete_sb_count++;
            }
        }

        /*********************************************Picture-based operations**********************************************************/
        // Delta QP range adjustments
        SetDefaultDeltaQpRange(
            context_ptr,
            picture_control_set_ptr,
            picture_control_set_ptr->slice_type == I_SLICE ? EB_FALSE : picture_control_set_ptr->scene_transition_flag[REF_LIST_0]);
        // Dark density derivation (histograms not available when no SCD)

        DeriveHighDarkAreaDensityFlag(
            sequence_control_set_ptr,
            picture_control_set_ptr);
        // Detect aura areas in lighter background when subject is moving similar to background
        DetermineMorePotentialAuraAreas(
            sequence_control_set_ptr,
            picture_control_set_ptr);

        // Activity statistics derivation
        DerivePictureActivityStatistics(
            sequence_control_set_ptr,
            picture_control_set_ptr);

        // Skin & Grass detection
        GrassSkinPicture(
            context_ptr,
            picture_control_set_ptr);
        // Get Empty Results Object
        eb_get_empty_object(
            context_ptr->picture_demux_results_output_fifo_ptr,
            &outputResultsWrapperPtr);

        outputResultsPtr = (PictureDemuxResults*)outputResultsWrapperPtr->object_ptr;
        outputResultsPtr->picture_control_set_wrapper_ptr = inputResultsPtr->picture_control_set_wrapper_ptr;
        outputResultsPtr->picture_type = EB_PIC_INPUT;

        // Release the Input Results
        eb_release_object(inputResultsWrapperPtr);

        // Post the Full Results Object
        eb_post_full_object(outputResultsWrapperPtr);
    }
    return EB_NULL;
}
