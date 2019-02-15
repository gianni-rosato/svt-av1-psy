/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureBufferDesc.h"

#include "EbSourceBasedOperationsProcess.h"
#include "EbInitialRateControlResults.h"
#include "EbPictureDemuxResults.h"
#include "EbPictureOperators.h"
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
#define CB_MEAN_RANGE_00                80


#define SAD_DEVIATION_LCU_TH        15
#define SAD_DEVIATION_LCU_NON_M4_TH 20

#define MAX_DELTA_QP_SHAPE_TH            4
#define MIN_DELTA_QP_SHAPE_TH           1

#define MIN_BLACK_AREA_PERCENTAGE        20
#define LOW_MEAN_THLD                    25

#define MIN_WHITE_AREA_PERCENTAGE         1
#define LOW_MEAN_THLD1                   40
#define HIGH_MEAN_THLD1                  210
#define NORM_FACTOR                      10 // Used ComplexityClassifier32x32

const uint32_t    THRESHOLD_NOISE[MAX_TEMPORAL_LAYERS] = { 33, 28, 27, 26, 26, 26 }; // [Temporal Layer Index]  // Used ComplexityClassifier32x32
// Outlier removal threshold per depth {2%, 2%, 4%, 4%}
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
    SourceBasedOperationsContext_t **context_dbl_ptr,
    EbFifo_t                        *initialRateControlResultsInputFifoPtr,
    EbFifo_t                        *picture_demux_results_output_fifo_ptr,
    SequenceControlSet_t            *sequence_control_set_ptr)
{
    SourceBasedOperationsContext_t *context_ptr;

    uint32_t  pictureLcuWidth = (sequence_control_set_ptr->max_input_luma_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
    uint32_t    pictureLcuHeight = (sequence_control_set_ptr->max_input_luma_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
    uint32_t    sb_total_count = pictureLcuWidth * pictureLcuHeight;

    EB_MALLOC(SourceBasedOperationsContext_t*, context_ptr, sizeof(SourceBasedOperationsContext_t), EB_N_PTR);
    *context_dbl_ptr = context_ptr;
    context_ptr->initial_rate_control_results_input_fifo_ptr = initialRateControlResultsInputFifoPtr;
    context_ptr->picture_demux_results_output_fifo_ptr = picture_demux_results_output_fifo_ptr;

    EB_MALLOC(uint8_t*, context_ptr->sb_high_contrast_array, sizeof(uint8_t) * sb_total_count, EB_N_PTR);

    return EB_ErrorNone;
}


/****************************************
* Init BEA QPM array to 0
** Used when no Lookahead is available
** and or zz_cost_array is invalid
****************************************/
void InitBeaQpmInfo(
    PictureParentControlSet_t        *picture_control_set_ptr,
    SequenceControlSet_t            *sequence_control_set_ptr)
{
    uint32_t sb_index;
    uint32_t zz_cost_average = 0, zzSum = 0;
    uint32_t complete_sb_count = 0;
    picture_control_set_ptr->low_motion_content_flag = EB_FALSE;
    picture_control_set_ptr->zz_cost_average = INVALID_ZZ_COST;


    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
        if (sb_params->is_complete_sb) {
            zzSum += picture_control_set_ptr->zz_cost_array[sb_index];
            complete_sb_count++;
        }
    }

    zz_cost_average = zzSum / complete_sb_count;
    picture_control_set_ptr->low_motion_content_flag = zz_cost_average == 0 ? EB_TRUE : EB_FALSE;
    picture_control_set_ptr->zz_cost_average = zz_cost_average;

    return;
}
/***************************************************
* Derives BEA statistics and set activity flags
***************************************************/
void DerivePictureActivityStatistics(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr)

{

    uint64_t               nonMovingIndexMin = ~0u;
    uint64_t               nonMovingIndexMax = 0;
    uint64_t               nonMovingIndexSum = 0;
#if NEW_PRED_STRUCT
    uint32_t               complete_sb_count = 0;
    uint32_t               non_moving_sb_count = 0;
#endif
    uint32_t               sb_total_count = picture_control_set_ptr->sb_total_count;
    uint32_t                 totNmvIdx = 0;

    uint32_t               sb_index;
    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {


        nonMovingIndexMin = picture_control_set_ptr->non_moving_index_array[sb_index] < nonMovingIndexMin ?
            picture_control_set_ptr->non_moving_index_array[sb_index] :
            nonMovingIndexMin;

        nonMovingIndexMax = picture_control_set_ptr->non_moving_index_array[sb_index] > nonMovingIndexMax ?
            picture_control_set_ptr->non_moving_index_array[sb_index] :
            nonMovingIndexMax;
#if NEW_PRED_STRUCT
        if (picture_control_set_ptr->non_moving_index_array[sb_index] < NON_MOVING_SCORE_1) {
            non_moving_sb_count++;
        }
        complete_sb_count++;
#endif
        nonMovingIndexSum += picture_control_set_ptr->non_moving_index_array[sb_index];

        if (picture_control_set_ptr->non_moving_index_array[sb_index] < 10)
            totNmvIdx++;
    }
    picture_control_set_ptr->non_moving_index_average = (uint16_t)(nonMovingIndexSum / sb_total_count);
#if NEW_PRED_STRUCT
    picture_control_set_ptr->kf_zeromotion_pct = (non_moving_sb_count * 100) / complete_sb_count;
#endif
    InitBeaQpmInfo(
        picture_control_set_ptr,
        sequence_control_set_ptr);

    return;
}
/***************************************************
* complexity Classification
***************************************************/
void ComplexityClassifier32x32(
    SequenceControlSet_t      *sequence_control_set_ptr,
    PictureParentControlSet_t *picture_control_set_ptr) {


    //No need to have Threshold depending on TempLayer. Higher Temoral Layer would have smaller classes
    //TODO: add more classes using above threshold when needed.


    uint32_t          sb_index;
    SbParams_t     *sb_params;
    uint8_t           noiseClass;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {


        picture_control_set_ptr->cmplx_status_sb[sb_index] = CMPLX_LOW;

        if (picture_control_set_ptr->temporal_layer_index >= 1 && picture_control_set_ptr->slice_type == B_SLICE) {

            sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];


            if (sb_params->is_complete_sb) {

                noiseClass = CMPLX_LOW;

                uint32_t         blkIt;

                for (blkIt = 0; blkIt < 4; blkIt++) {


                    uint32_t distortion = picture_control_set_ptr->me_results[sb_index][1 + blkIt].distortionDirection[0].distortion;

                    if ((((uint32_t)(distortion)) >> NORM_FACTOR) > THRESHOLD_NOISE[picture_control_set_ptr->temporal_layer_index])


                        noiseClass++;
                }


                picture_control_set_ptr->cmplx_status_sb[sb_index] = noiseClass > 0 ? CMPLX_NOISE : CMPLX_LOW;
            }
        }
    }
}





/******************************************************
* Pre-MD Uncovered Area Detection
******************************************************/
void FailingMotionLcu(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t        *picture_control_set_ptr,
    uint32_t                             sb_index) {

    uint32_t rasterScanCuIndex;

    // SB Loop : Failing motion detector for L2 only
    SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    // Detection variables
    uint64_t                  sortedcuOisSAD = 0;
    uint64_t                  cuMeSAD = 0;
    int64_t                  meToOisSadDeviation = 0;
    // SB loop variables

    int64_t failing_motion_sb_flag = 0;

    if (picture_control_set_ptr->slice_type != I_SLICE && sb_params->is_complete_sb && (!picture_control_set_ptr->similar_colocated_sb_array[sb_index])) {
        for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_64x64; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_32x32_3; rasterScanCuIndex++) {

            meToOisSadDeviation = 0;

            // Get ME SAD

            cuMeSAD = picture_control_set_ptr->me_results[sb_index][rasterScanCuIndex].distortionDirection[0].distortion;




            OisCu32Cu16Results_t *oisResultsPtr = picture_control_set_ptr->ois_cu32_cu16_results[sb_index];
            if (RASTER_SCAN_CU_SIZE[rasterScanCuIndex] > 32) {
                sortedcuOisSAD = oisResultsPtr->sorted_ois_candidate[1][0].distortion +
                    oisResultsPtr->sorted_ois_candidate[2][0].distortion +
                    oisResultsPtr->sorted_ois_candidate[3][0].distortion +
                    oisResultsPtr->sorted_ois_candidate[4][0].distortion;
            }
            else { //32x32
                sortedcuOisSAD = oisResultsPtr->sorted_ois_candidate[rasterScanCuIndex][0].distortion;
            }



            int64_t  meToOisSadDiff = (int32_t)cuMeSAD - (int32_t)sortedcuOisSAD;
            meToOisSadDeviation = (sortedcuOisSAD == 0) || (meToOisSadDiff < 0) ? 0 : (meToOisSadDiff * 100) / sortedcuOisSAD;

            if (meToOisSadDeviation > SAD_DEVIATION_LCU_TH) {
                failing_motion_sb_flag += 1;
            }
        }

        // Update failing motion flag
        picture_control_set_ptr->failing_motion_sb_flag[sb_index] = (failing_motion_sb_flag && (picture_control_set_ptr->intensity_transition_flag == EB_FALSE)) ? EB_TRUE : EB_FALSE;
    }

}

/******************************************************
* Pre-MD Uncovered Area Detection
******************************************************/
void DetectUncoveredLcu(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t        *picture_control_set_ptr,
    uint32_t                             sb_index) {

    uint32_t rasterScanCuIndex;

    // SB Loop : Uncovered area detector -- ON only for 4k
    SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    // Detection variables
    uint64_t                  sortedcuOisSAD = 0;
    uint64_t                  cuMeSAD = 0;
    int64_t                  meToOisSadDeviation = 0;
    // SB loop variables

    int64_t uncovered_area_sb_flag = 0;


    if (picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {
        if (sb_params->is_complete_sb && (!picture_control_set_ptr->similar_colocated_sb_array[sb_index])) {


            for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_64x64; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_32x32_3; rasterScanCuIndex++) {


                meToOisSadDeviation = 0;

                // Get ME SAD

                cuMeSAD = picture_control_set_ptr->me_results[sb_index][rasterScanCuIndex].distortionDirection[0].distortion;




                OisCu32Cu16Results_t *oisResultsPtr = picture_control_set_ptr->ois_cu32_cu16_results[sb_index];
                if (RASTER_SCAN_CU_SIZE[rasterScanCuIndex] > 32) {
                    sortedcuOisSAD = oisResultsPtr->sorted_ois_candidate[1][0].distortion +
                        oisResultsPtr->sorted_ois_candidate[2][0].distortion +
                        oisResultsPtr->sorted_ois_candidate[3][0].distortion +
                        oisResultsPtr->sorted_ois_candidate[4][0].distortion;
                }
                else if (RASTER_SCAN_CU_SIZE[rasterScanCuIndex] == 32) {
                    sortedcuOisSAD = oisResultsPtr->sorted_ois_candidate[rasterScanCuIndex][0].distortion;
                }
                else {
                    sortedcuOisSAD = oisResultsPtr->sorted_ois_candidate[rasterScanCuIndex][0].distortion;
                }



                int64_t  meToOisSadDiff = (int32_t)cuMeSAD - (int32_t)sortedcuOisSAD;
                meToOisSadDeviation = (sortedcuOisSAD == 0) || (meToOisSadDiff < 0) ? 0 : (meToOisSadDiff * 100) / sortedcuOisSAD;


                if (RASTER_SCAN_CU_SIZE[rasterScanCuIndex] > 16) {

                    if (meToOisSadDeviation > SAD_DEVIATION_LCU_NON_M4_TH) {
                        uncovered_area_sb_flag += 1;
                    }

                }
            }

            // Update Uncovered area flag
            picture_control_set_ptr->uncovered_area_sb_flag[sb_index] = (uncovered_area_sb_flag && (picture_control_set_ptr->intensity_transition_flag == EB_FALSE)) ? EB_TRUE : EB_FALSE;
        }

    }
}


/******************************************************
* Calculates AC Energy
******************************************************/
void CalculateAcEnergy(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t        *picture_control_set_ptr,
    uint32_t                             sb_index) {

    EbPictureBufferDesc_t    *inputPicturePtr = picture_control_set_ptr->enhanced_picture_ptr;
    uint32_t                     inputLumaStride = inputPicturePtr->strideY;
    uint32_t                   inputOriginIndex;
    SbParams_t  *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

    uint8_t       *meanPtr = picture_control_set_ptr->yMean[sb_index];
    inputOriginIndex = (sb_params->origin_y + inputPicturePtr->origin_y) * inputLumaStride + (sb_params->origin_x + inputPicturePtr->origin_x);

    if (sb_params->is_complete_sb && picture_control_set_ptr->slice_type == I_SLICE) {

        uint32_t inputCuOriginIndex;
        uint32_t cuNum, cu_size;
        uint16_t cuH, cuW;


        picture_control_set_ptr->sb_y_src_energy_cu_array[sb_index][0] = ComputeNxMSatdSadLCU(
            &(inputPicturePtr->bufferY[inputOriginIndex]),
            inputPicturePtr->strideY,
            sb_params->width,
            sb_params->height,
            sequence_control_set_ptr->encode_context_ptr->asm_type);

        //64 x 64
        picture_control_set_ptr->sb_y_src_mean_cu_array[sb_index][0] = meanPtr[0];

        // 32x32
        cu_size = 32;
        cuNum = 64 / cu_size;
        for (cuH = 0; cuH < cuNum; cuH++) {
            for (cuW = 0; cuW < cuNum; cuW++) {
                inputCuOriginIndex = inputOriginIndex + cuH * (64 / cuNum)*inputLumaStride + cuW * (64 / cuNum);

                picture_control_set_ptr->sb_y_src_energy_cu_array[sb_index][1 + cuH * cuNum + cuW] = ComputeNxMSatdSadLCU(
                    &(inputPicturePtr->bufferY[inputCuOriginIndex]),
                    inputPicturePtr->strideY,
                    cu_size,
                    cu_size,
                    sequence_control_set_ptr->encode_context_ptr->asm_type);
                picture_control_set_ptr->sb_y_src_mean_cu_array[sb_index][1 + cuH * cuNum + cuW] = meanPtr[1 + cuH * cuNum + cuW];
            }
        }


    }
    else {
        picture_control_set_ptr->sb_y_src_energy_cu_array[sb_index][0] = 100000000;
        picture_control_set_ptr->sb_y_src_energy_cu_array[sb_index][1] = 100000000;
        picture_control_set_ptr->sb_y_src_energy_cu_array[sb_index][2] = 100000000;
        picture_control_set_ptr->sb_y_src_energy_cu_array[sb_index][3] = 100000000;
        picture_control_set_ptr->sb_y_src_energy_cu_array[sb_index][4] = 100000000;
        picture_control_set_ptr->sb_y_src_mean_cu_array[sb_index][0] = 100000000;
        picture_control_set_ptr->sb_y_src_mean_cu_array[sb_index][1] = 100000000;
        picture_control_set_ptr->sb_y_src_mean_cu_array[sb_index][2] = 100000000;
        picture_control_set_ptr->sb_y_src_mean_cu_array[sb_index][3] = 100000000;
        picture_control_set_ptr->sb_y_src_mean_cu_array[sb_index][4] = 100000000;
    }
}

void LumaContrastDetectorLcu(
    SourceBasedOperationsContext_t *context_ptr,
    SequenceControlSet_t           *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    uint32_t                            sb_index) {

    uint64_t                  cuOisSAD = 0;
    uint64_t                  cuMeSAD = 0;

    // Calculate Luma mean of the frame by averaging the mean of LCUs to Detect Dark Frames (On only for 4k and BQMode)
    uint8_t  *y_mean_ptr = context_ptr->y_mean_ptr;
    SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    if (sb_params->is_complete_sb) {
        if (picture_control_set_ptr->slice_type != I_SLICE && picture_control_set_ptr->temporal_layer_index == 0) {




            OisCu32Cu16Results_t *oisResultsPtr = picture_control_set_ptr->ois_cu32_cu16_results[sb_index];
            cuOisSAD = oisResultsPtr->sorted_ois_candidate[1][0].distortion +
                oisResultsPtr->sorted_ois_candidate[2][0].distortion +
                oisResultsPtr->sorted_ois_candidate[3][0].distortion +
                oisResultsPtr->sorted_ois_candidate[4][0].distortion;



            cuMeSAD = picture_control_set_ptr->me_results[sb_index][0].distortionDirection[0].distortion;


            context_ptr->to_be_intra_coded_probability += cuOisSAD < cuMeSAD ? 1 : 0;
            context_ptr->depth1_block_num++;
        }
    }

    if (picture_control_set_ptr->non_moving_index_array[sb_index] < 10)
    {
        context_ptr->y_non_moving_mean += y_mean_ptr[0];
        context_ptr->countOfNonMovingLcus++;
    }
    else {
        context_ptr->y_moving_mean += y_mean_ptr[0];
        context_ptr->count_of_moving_sbs++;
    }
}

void LumaContrastDetectorPicture(
    SourceBasedOperationsContext_t        *context_ptr,
    PictureParentControlSet_t            *picture_control_set_ptr) {

    context_ptr->y_non_moving_mean = (context_ptr->countOfNonMovingLcus != 0) ? (context_ptr->y_non_moving_mean / context_ptr->countOfNonMovingLcus) : 0;
    context_ptr->y_moving_mean = (context_ptr->count_of_moving_sbs != 0) ? (context_ptr->y_moving_mean / context_ptr->count_of_moving_sbs) : 0;

    picture_control_set_ptr->dark_back_groundlight_fore_ground = ((context_ptr->y_moving_mean > (2 * context_ptr->y_non_moving_mean)) && (context_ptr->y_non_moving_mean < DARK_FRM_TH)) ?
        EB_TRUE :
        EB_FALSE;

    picture_control_set_ptr->intra_coded_block_probability = 0;

    if (picture_control_set_ptr->slice_type != I_SLICE && picture_control_set_ptr->temporal_layer_index == 0) {
        picture_control_set_ptr->intra_coded_block_probability = (uint8_t)(context_ptr->depth1_block_num != 0 ? context_ptr->to_be_intra_coded_probability * 100 / context_ptr->depth1_block_num : 0);
    }
}

void GrassLcu(
    SourceBasedOperationsContext_t        *context_ptr,
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureParentControlSet_t            *picture_control_set_ptr,
    uint32_t                                 sb_index) {

    uint32_t                  childIndex;

    EbBool                 lcuGrassFlag = EB_FALSE;

    uint32_t grassLcuInrange;
    uint32_t processedCus;
    uint32_t  rasterScanCuIndex;

    SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    SbStat_t *sb_stat_ptr = &(picture_control_set_ptr->sb_stat_array[sb_index]);

    _mm_prefetch((const char*)sb_stat_ptr, _MM_HINT_T0);

    lcuGrassFlag = EB_FALSE;
    grassLcuInrange = 0;
    processedCus = 0;


    for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++) {
        if (sb_params->raster_scan_cu_validity[rasterScanCuIndex]) {
            const uint32_t mdScanCuIndex = RASTER_SCAN_TO_MD_SCAN[rasterScanCuIndex];
            const uint32_t rasterScanParentCuIndex = RASTER_SCAN_CU_PARENT_INDEX[rasterScanCuIndex];
            const uint32_t mdScanParentCuIndex = RASTER_SCAN_TO_MD_SCAN[rasterScanParentCuIndex];
            CuStat_t *cuStatPtr = &(sb_stat_ptr->cu_stat_array[mdScanCuIndex]);


            const uint32_t perfectCondition = 7;
            const uint8_t yMean = context_ptr->y_mean_ptr[rasterScanCuIndex];
            const uint8_t cbMean = context_ptr->cb_mean_ptr[rasterScanCuIndex];
            const uint8_t crMean = context_ptr->cr_mean_ptr[rasterScanCuIndex];
            uint32_t grassCondition = 0;
            uint32_t skinCondition = 0;


            // GRASS
            grassCondition += (yMean > Y_MEAN_RANGE_02 && yMean < Y_MEAN_RANGE_01) ? 1 : 0;
            grassCondition += (cbMean > CB_MEAN_RANGE_00 && cbMean < CB_MEAN_RANGE_02) ? 2 : 0;
            grassCondition += (crMean > CR_MEAN_RANGE_00 && crMean < CR_MEAN_RANGE_02) ? 4 : 0;


            grassLcuInrange += (grassCondition == perfectCondition) ? 1 : 0;
            processedCus++;

            lcuGrassFlag = grassCondition == perfectCondition ? EB_TRUE : lcuGrassFlag;

            cuStatPtr->grass_area = (EbBool)(grassCondition == perfectCondition);
            // SKIN
            skinCondition += (yMean > Y_MEAN_RANGE_02 && yMean < Y_MEAN_RANGE_01) ? 1 : 0;
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
    SourceBasedOperationsContext_t        *context_ptr,
    PictureParentControlSet_t            *picture_control_set_ptr) {
    picture_control_set_ptr->grass_percentage_in_picture = (uint8_t)(context_ptr->picture_num_grass_sb * 100 / picture_control_set_ptr->sb_total_count);
}

/******************************************************
* Detect and mark SB and 32x32 CUs which belong to an isolated non-homogeneous region surrounding a homogenous and flat region
******************************************************/
void DetermineIsolatedNonHomogeneousRegionInPicture(
    SequenceControlSet_t            *sequence_control_set_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr)
{
    uint32_t sb_index;
    uint32_t cuuIndex;
    int32_t lcuHor, lcuVer, lcuVerOffset;
    uint32_t sb_total_count = picture_control_set_ptr->sb_total_count;
    uint32_t picture_width_in_sb = sequence_control_set_ptr->picture_width_in_sb;
    uint32_t picture_height_in_sb = sequence_control_set_ptr->picture_height_in_sb;

    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {

        SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
        // Initialize
        picture_control_set_ptr->sb_isolated_non_homogeneous_area_array[sb_index] = EB_FALSE;
        if ((sb_params->horizontal_index > 0) && (sb_params->horizontal_index < picture_width_in_sb - 1) && (sb_params->vertical_index > 0) && (sb_params->vertical_index < picture_height_in_sb - 1)) {
            uint32_t countOfMedVarianceLcu;
            countOfMedVarianceLcu = 0;

            // top neighbors
            countOfMedVarianceLcu += ((picture_control_set_ptr->variance[sb_index - picture_width_in_sb - 1][ME_TIER_ZERO_PU_64x64]) <= MEDIUM_LCU_VARIANCE) ? 1 : 0;
            countOfMedVarianceLcu += (picture_control_set_ptr->variance[sb_index - picture_width_in_sb][ME_TIER_ZERO_PU_64x64] <= MEDIUM_LCU_VARIANCE) ? 1 : 0;
            countOfMedVarianceLcu += (picture_control_set_ptr->variance[sb_index - picture_width_in_sb + 1][ME_TIER_ZERO_PU_64x64] <= MEDIUM_LCU_VARIANCE && sequence_control_set_ptr->sb_params_array[sb_index - picture_width_in_sb + 1].is_complete_sb) ? 1 : 0;
            // bottom
            countOfMedVarianceLcu += (picture_control_set_ptr->variance[sb_index + picture_width_in_sb - 1][ME_TIER_ZERO_PU_64x64] <= MEDIUM_LCU_VARIANCE && sequence_control_set_ptr->sb_params_array[sb_index + picture_width_in_sb - 1].is_complete_sb) ? 1 : 0;
            countOfMedVarianceLcu += (picture_control_set_ptr->variance[sb_index + picture_width_in_sb][ME_TIER_ZERO_PU_64x64] <= MEDIUM_LCU_VARIANCE && sequence_control_set_ptr->sb_params_array[sb_index + picture_width_in_sb].is_complete_sb) ? 1 : 0;
            countOfMedVarianceLcu += (picture_control_set_ptr->variance[sb_index + picture_width_in_sb + 1][ME_TIER_ZERO_PU_64x64] <= MEDIUM_LCU_VARIANCE && sequence_control_set_ptr->sb_params_array[sb_index + picture_width_in_sb + 1].is_complete_sb) ? 1 : 0;
            // left right
            countOfMedVarianceLcu += (picture_control_set_ptr->variance[sb_index + 1][ME_TIER_ZERO_PU_64x64] <= MEDIUM_LCU_VARIANCE && sequence_control_set_ptr->sb_params_array[sb_index + 1].is_complete_sb) ? 1 : 0;
            countOfMedVarianceLcu += (picture_control_set_ptr->variance[sb_index - 1][ME_TIER_ZERO_PU_64x64] <= MEDIUM_LCU_VARIANCE) ? 1 : 0;

            // At least two neighbors are flat
            if ((countOfMedVarianceLcu > 2) || countOfMedVarianceLcu > 1)
            {
                // Search within an SB if any of the 32x32 CUs is non-homogeneous
                uint32_t count32x32NonhomCusInLcu = 0;
                for (cuuIndex = 0; cuuIndex < 4; cuuIndex++)
                {
                    if (picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][cuuIndex] > VAR_BASED_DETAIL_PRESERVATION_SELECTOR_THRSLHD)
                        count32x32NonhomCusInLcu++;
                }
                // If atleast one is non-homogeneous, then check all its neighbors (top left, top, top right, left, right, btm left, btm, btm right)
                uint32_t countOfHomogeneousNeighborLcus = 0;
                if (count32x32NonhomCusInLcu > 0) {
                    for (lcuVer = -1; lcuVer <= 1; lcuVer++) {
                        lcuVerOffset = lcuVer * (int32_t)picture_width_in_sb;
                        for (lcuHor = -1; lcuHor <= 1; lcuHor++) {
                            if (lcuVer != 0 && lcuHor != 0)
                                countOfHomogeneousNeighborLcus += (picture_control_set_ptr->sb_homogeneous_area_array[sb_index + lcuVerOffset + lcuHor] == EB_TRUE);

                        }
                    }
                }

                // To determine current lcu is isolated non-homogeneous, at least 2 neighbors must be homogeneous
                if (countOfHomogeneousNeighborLcus >= 2) {
                    for (cuuIndex = 0; cuuIndex < 4; cuuIndex++)
                    {
                        if (picture_control_set_ptr->var_of_var32x32_based_sb_array[sb_index][cuuIndex] > VAR_BASED_DETAIL_PRESERVATION_SELECTOR_THRSLHD)
                        {
                            picture_control_set_ptr->sb_isolated_non_homogeneous_area_array[sb_index] = EB_TRUE;
                        }
                    }
                }

            }
        }
    }
    return;
}


void SetDefaultDeltaQpRange(
    SourceBasedOperationsContext_t    *context_ptr,
    PictureParentControlSet_t        *picture_control_set_ptr,
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
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t    *picture_control_set_ptr)
{
    uint16_t sb_index;
    int32_t lcuHor, lcuVer, lcuVerOffset;
    uint8_t  sb_x, sb_y;
    uint32_t countOfEdgeBlocks = 0, countOfNonEdgeBlocks = 0;

    uint32_t lightLumaValue = 150;

    uint16_t sb_total_count = picture_control_set_ptr->sb_total_count;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

        sb_x = sb_params->horizontal_index;
        sb_y = sb_params->vertical_index;
        // For all the internal LCUs
        if ((sb_x > 0) && (sb_x < sequence_control_set_ptr->picture_width_in_sb - 1) && (sb_y > 0) && (sb_y < sequence_control_set_ptr->picture_height_in_sb - 1)) {
            countOfNonEdgeBlocks = 0;
            if (picture_control_set_ptr->edge_results_ptr[sb_index].edge_block_num
                && picture_control_set_ptr->yMean[sb_index][ME_TIER_ZERO_PU_64x64] >= lightLumaValue) {

                for (lcuVer = -1; lcuVer <= 1; lcuVer++) {
                    lcuVerOffset = lcuVer * (int32_t)sequence_control_set_ptr->picture_width_in_sb;
                    for (lcuHor = -1; lcuHor <= 1; lcuHor++) {
                        countOfNonEdgeBlocks += (!picture_control_set_ptr->edge_results_ptr[sb_index + lcuVerOffset + lcuHor].edge_block_num) &&
                            (picture_control_set_ptr->non_moving_index_array[sb_index + lcuVerOffset + lcuHor] < 30);
                    }
                }
            }

            if (countOfNonEdgeBlocks > 1) {
                countOfEdgeBlocks++;
            }
        }
    }


    // To check the percentage of potential aura in the picture.. If a large area is detected then this is not isolated
    picture_control_set_ptr->percentage_of_edgein_light_background = (uint8_t)(countOfEdgeBlocks * 100 / sb_total_count);

}

/***************************************************
* Detects the presence of dark area
***************************************************/
void DeriveHighDarkAreaDensityFlag(
    SequenceControlSet_t                *sequence_control_set_ptr,
    PictureParentControlSet_t           *picture_control_set_ptr) {


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

    blackAreaPercentage = (blackSamplesCount * 100) / (sequence_control_set_ptr->luma_width * sequence_control_set_ptr->luma_height);
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
            for (lumaHistogramBin = HIGH_MEAN_THLD1; lumaHistogramBin < HISTOGRAM_NUMBER_OF_BINS; lumaHistogramBin++) {
                whiteSamplesCount += picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][lumaHistogramBin];
            }
        }
    }

    blackAreaPercentage = (blackSamplesCount * 100) / (sequence_control_set_ptr->luma_width * sequence_control_set_ptr->luma_height);
    whiteAreaPercentage = (whiteSamplesCount * 100) / (sequence_control_set_ptr->luma_width * sequence_control_set_ptr->luma_height);
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

/******************************************************
* High  contrast classifier
******************************************************/
void TemporalHighContrastClassifier(
    SourceBasedOperationsContext_t    *context_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    uint32_t                           sb_index)
{

    uint32_t blkIt;
    uint32_t nsadTable[] = { 10, 5, 5, 5, 5, 5 };
    uint32_t thRes = 0;
    uint32_t nsad;
    uint32_t meDist = 0;

    if (picture_control_set_ptr->slice_type == B_SLICE) {


        for (blkIt = 0; blkIt < 4; blkIt++) {

            nsad = ((uint32_t)picture_control_set_ptr->me_results[sb_index][1 + blkIt].distortionDirection[0].distortion) >> NORM_FACTOR;

            if (nsad >= nsadTable[picture_control_set_ptr->temporal_layer_index] + thRes)
                meDist++;
        }

    }
    context_ptr->high_dist = meDist > 0 ? EB_TRUE : EB_FALSE;
}

void SpatialHighContrastClassifier(
    SourceBasedOperationsContext_t    *context_ptr,
    PictureParentControlSet_t       *picture_control_set_ptr,
    uint32_t                           sb_index)
{

    uint32_t                 blkIt;

    context_ptr->high_contrast_num = 0;
    context_ptr->high_contrast_num_ii = 0;
    //16x16 blocks
    for (blkIt = 0; blkIt < 16; blkIt++) {

        uint8_t ymean = context_ptr->y_mean_ptr[5 + blkIt];
        uint8_t umean = context_ptr->cb_mean_ptr[5 + blkIt];
        uint8_t vmean = context_ptr->cr_mean_ptr[5 + blkIt];

        uint16_t var = picture_control_set_ptr->variance[sb_index][5 + blkIt];


        if (var > VAR_MIN && var<VAR_MAX            &&  //medium texture
            ymean>MIN_Y && ymean < MAX_Y            &&  //medium brightness(not too dark and not too bright)
            ABS((int64_t)umean - MID_CB) < TH_CB &&  //middle of the color plane
            ABS((int64_t)vmean - MID_CR) < TH_CR     //middle of the color plane
            )
        {
            context_ptr->high_contrast_num++;

        }


        if (
            ymean < 30 &&  //medium brightness(not too dark and not too bright)
            ABS((int64_t)umean - 128) < 5 &&  //middle of the color plane
            ABS((int64_t)vmean - 128) < 5     //middle of the color plane
            )
        {
            context_ptr->high_contrast_num_ii++;

        }
    }
}

void DeriveComplexityContrastPicture(
    SourceBasedOperationsContext_t    *context_ptr,
    SequenceControlSet_t         *sequence_control_set_ptr,
    PictureParentControlSet_t    *picture_control_set_ptr)

{
    uint32_t    sb_index;

    //look only for isolated shapes.
    if ((context_ptr->sb_cmplx_contrast_count * 100) / context_ptr->complete_sb_count > 10) {
        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
            picture_control_set_ptr->sb_cmplx_contrast_array[sb_index] = 0;
        }
    }

    //look only for isolated shapes.
    if ((context_ptr->sb_high_contrast_count * 100) / context_ptr->complete_sb_count <= 10) {

        for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index)
        {
            SbParams_t     *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            EbBool isCentralArea = EB_FALSE;
            isCentralArea = (EbBool)(sb_params->origin_y > 16 * BLOCK_SIZE_64 && sb_params->origin_y < 21 * BLOCK_SIZE_64 && sb_params->origin_x> 12 * BLOCK_SIZE_64 && sb_params->origin_x < 32 * 64);
            if (context_ptr->sb_high_contrast_array[sb_index] > 0 && isCentralArea) {
                picture_control_set_ptr->sb_high_contrast_array_dialated[sb_index] = 4;
                int32_t i, j;
                uint8_t * ptr = &picture_control_set_ptr->sb_high_contrast_array_dialated[(int32_t)sb_index - (int32_t)sequence_control_set_ptr->picture_width_in_sb * 1 - 1];

                if (sb_index == (uint32_t)(picture_control_set_ptr->sb_total_count - 1))
                {
                    ptr[0] = 4;
                    ptr[1] = 4;
                    ptr[2] = 4;
                    ptr[sequence_control_set_ptr->picture_width_in_sb] = 4;
                    ptr[1 + sequence_control_set_ptr->picture_width_in_sb] = 4;
                }
                else
                {
                    for (j = 0; j < 2; j++) {
                        for (i = 0; i < 3; i++) {

                            ptr[i + j * sequence_control_set_ptr->picture_width_in_sb] = 4;
                        }
                    }
                }
            }
        }
    }

    return;
}

/******************************************************
* Detect Cu32x32 Clean Sparse Array
******************************************************/
void DetectCu32x32CleanSparseLcu(
    SequenceControlSet_t        *sequence_control_set_ptr,
    PictureParentControlSet_t    *picture_control_set_ptr,
    uint32_t                         sb_index)
{
    int32_t  blockIndex, blockIndexX, blockIndexY, cu32x32Index;
    uint16_t * variancePtr;
    uint8_t * meanPtr;
    SbParams_t     *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

    if (sb_params->is_complete_sb) {
        variancePtr = picture_control_set_ptr->variance[sb_index];
        meanPtr = picture_control_set_ptr->yMean[sb_index];
        for (blockIndexY = 0; blockIndexY < 2; ++blockIndexY) {
            for (blockIndexX = 0; blockIndexX < 2; ++blockIndexX) {
                cu32x32Index = blockIndexX + (blockIndexY << 1) + 1;
                blockIndex = ((sb_params->horizontal_index << 1) + blockIndexX) + ((sb_params->vertical_index << 1) + blockIndexY) * picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array_stride;

                if ((meanPtr[cu32x32Index]) > 130 && (meanPtr[cu32x32Index]) < 220 && (variancePtr[cu32x32Index]) < 30) {

                    picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array[blockIndex] = 2;
                }
                else if ((meanPtr[cu32x32Index]) > 130 && (meanPtr[cu32x32Index]) < 220) {
                    picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array[blockIndex] = 1;
                }
            }
        }
    }

    return;
}

void DetectCu32x32CleanSparsePicture(
    PictureParentControlSet_t    *picture_control_set_ptr)
{

    int32_t  blockIndex, blockIndexX, blockIndexY;
    int32_t  blockIndexTemp, blockIndexXTemp, blockIndexYTemp;
    uint32_t  rowNumber = picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array_size / picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array_stride;

    for (blockIndexY = 1; blockIndexY < (int32_t)rowNumber - 1; ++blockIndexY) {
        for (blockIndexX = 1; blockIndexX < picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array_stride - 1; ++blockIndexX) {
            blockIndex = (blockIndexX)+(blockIndexY)* picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array_stride;
            uint8_t neighCount = 0;
            if (picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array[blockIndex] == 1)
            {
                for (blockIndexYTemp = -1; blockIndexYTemp <= 1; ++blockIndexYTemp) {
                    for (blockIndexXTemp = -1; blockIndexXTemp <= 1; ++blockIndexXTemp) {
                        blockIndexTemp = blockIndex + blockIndexXTemp + blockIndexYTemp * picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array_stride;
                        if (picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array[blockIndexTemp] == 2) {
                            neighCount++;

                        }
                    }
                }
                if (neighCount >= 6)
                    picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array[blockIndex] = 2;
            }
        }
    }

    return;
}

/************************************************
 * Source Based Operations Kernel
 ************************************************/
void* source_based_operations_kernel(void *input_ptr)
{
    SourceBasedOperationsContext_t    *context_ptr = (SourceBasedOperationsContext_t*)input_ptr;
    PictureParentControlSet_t       *picture_control_set_ptr;
    SequenceControlSet_t            *sequence_control_set_ptr;
    EbObjectWrapper_t               *inputResultsWrapperPtr;
    InitialRateControlResults_t        *inputResultsPtr;
    EbObjectWrapper_t               *outputResultsWrapperPtr;
    PictureDemuxResults_t           *outputResultsPtr;

    for (;;) {

        // Get Input Full Object
        EbGetFullObject(
            context_ptr->initial_rate_control_results_input_fifo_ptr,
            &inputResultsWrapperPtr);

        inputResultsPtr = (InitialRateControlResults_t*)inputResultsWrapperPtr->objectPtr;
        picture_control_set_ptr = (PictureParentControlSet_t*)inputResultsPtr->pictureControlSetWrapperPtr->objectPtr;
        sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;

        picture_control_set_ptr->dark_back_groundlight_fore_ground = EB_FALSE;
        context_ptr->picture_num_grass_sb = 0;

        context_ptr->sb_cmplx_contrast_count = 0;
        context_ptr->sb_high_contrast_count = 0;
        context_ptr->complete_sb_count = 0;

        context_ptr->count_of_moving_sbs = 0;
        context_ptr->countOfNonMovingLcus = 0;
        context_ptr->y_non_moving_mean = 0;
        context_ptr->y_moving_mean = 0;
        context_ptr->to_be_intra_coded_probability = 0;
        context_ptr->depth1_block_num = 0;

        // Reset the cu 32x32 array for Clean Sparse flag
        EB_MEMSET(picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array, 0, picture_control_set_ptr->cu32x32_clean_sparse_coeff_map_array_size);

        uint32_t sb_total_count = picture_control_set_ptr->sb_total_count;
        uint32_t sb_index;

        /***********************************************LCU-based operations************************************************************/
        for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
            SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            picture_control_set_ptr->sb_cmplx_contrast_array[sb_index] = 0;
            context_ptr->sb_high_contrast_array[sb_index] = 0;
            picture_control_set_ptr->sb_high_contrast_array_dialated[sb_index] = 0;
            EbBool is_complete_sb = sb_params->is_complete_sb;
            uint8_t  *y_mean_ptr = picture_control_set_ptr->yMean[sb_index];

            _mm_prefetch((const char*)y_mean_ptr, _MM_HINT_T0);

            // 32x32 spare coefficient detection
            if (picture_control_set_ptr->slice_type == I_SLICE) {
                DetectCu32x32CleanSparseLcu(
                    sequence_control_set_ptr,
                    picture_control_set_ptr,
                    sb_index);
            }
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


            // Luma Contrast detection
            LumaContrastDetectorLcu(
                context_ptr,
                sequence_control_set_ptr,
                picture_control_set_ptr,
                sb_index);

            // AC energy computation
            CalculateAcEnergy(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                sb_index);

            // Failing Motion Detection
            picture_control_set_ptr->failing_motion_sb_flag[sb_index] = EB_FALSE;

            if (picture_control_set_ptr->slice_type != I_SLICE && is_complete_sb) {

                FailingMotionLcu(
                    sequence_control_set_ptr,
                    picture_control_set_ptr,
                    sb_index);
            }

            picture_control_set_ptr->uncovered_area_sb_flag[sb_index] = EB_FALSE;
            if (picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->slice_type != I_SLICE) {

                if (is_complete_sb && (!picture_control_set_ptr->similar_colocated_sb_array[sb_index])) {
                    DetectUncoveredLcu(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        sb_index);
                }
            }

            // Uncovered area detection II

            // Temporal high contrast classifier
            if (is_complete_sb) {
                TemporalHighContrastClassifier(
                    context_ptr,
                    picture_control_set_ptr,
                    sb_index);

                if (context_ptr->high_contrast_num > 0 && context_ptr->high_dist == EB_TRUE) {
                    picture_control_set_ptr->sb_cmplx_contrast_array[sb_index] = 4;
                    context_ptr->sb_cmplx_contrast_count++;
                }

                if ((context_ptr->high_dist == EB_TRUE && context_ptr->high_contrast_num_ii > 0) || picture_control_set_ptr->sb_cmplx_contrast_array[sb_index] == 4) {
                    context_ptr->sb_high_contrast_array[sb_index] = 4;
                    context_ptr->sb_high_contrast_count++;
                }

                context_ptr->complete_sb_count++;
            }

        }

        /*********************************************Picture-based operations**********************************************************/
        LumaContrastDetectorPicture(
            context_ptr,
            picture_control_set_ptr);

        if (picture_control_set_ptr->slice_type == I_SLICE) {
            DetectCu32x32CleanSparsePicture(
                picture_control_set_ptr);
        }

        DeriveComplexityContrastPicture(
            context_ptr,
            sequence_control_set_ptr,
            picture_control_set_ptr);

        // Delta QP range adjustments
        SetDefaultDeltaQpRange(
            context_ptr,
            picture_control_set_ptr,
            picture_control_set_ptr->slice_type == I_SLICE ? EB_FALSE : picture_control_set_ptr->scene_transition_flag[REF_LIST_0]);
        // Dark density derivation (histograms not available when no SCD)

        DeriveHighDarkAreaDensityFlag(
            sequence_control_set_ptr,
            picture_control_set_ptr);

        // Detect and mark SB and 32x32 CUs which belong to an isolated non-homogeneous region surrounding a homogenous and flat region.
        DetermineIsolatedNonHomogeneousRegionInPicture(
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

        // Complexity Classification
        ComplexityClassifier32x32(
            sequence_control_set_ptr,
            picture_control_set_ptr);

        // Get Empty Results Object
        EbGetEmptyObject(
            context_ptr->picture_demux_results_output_fifo_ptr,
            &outputResultsWrapperPtr);

        outputResultsPtr = (PictureDemuxResults_t*)outputResultsWrapperPtr->objectPtr;
        outputResultsPtr->pictureControlSetWrapperPtr = inputResultsPtr->pictureControlSetWrapperPtr;
        outputResultsPtr->pictureType = EB_PIC_INPUT;

        // Release the Input Results
        EbReleaseObject(inputResultsWrapperPtr);

        // Post the Full Results Object
        EbPostFullObject(outputResultsWrapperPtr);

    }
    return EB_NULL;
}
