/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/


#include "EbModeDecisionConfiguration.h"
#include "EbRateDistortionCost.h"
#include "EbUtility.h"
#include "EbModeDecisionProcess.h"
#include "EbDefinitions.h"
#if !OPEN_LOOP_EARLY_PARTITION
static const uint32_t me2Nx2NOffset[4] = { 0, 1, 5, 21 };
#endif
/********************************************
 * Constants
 ********************************************/
#if OPEN_LOOP_EARLY_PARTITION
int pa_to_ep_block_index[85] = {

    0    ,

    25   ,

    50   ,
    75   ,    84   ,    93   ,    102  ,
    111  ,
    136  ,    145  ,    154  ,    163  ,
    172  ,
    197  ,    206  ,    215  ,    224  ,
    233  ,
    258  ,    267  ,    276  ,    285  ,

    294  ,

    319  ,
    344  ,    353  ,    362  ,    371  ,
    380  ,
    405  ,    414  ,    423  ,    432  ,
    441  ,
    466  ,    475  ,    484  ,   493  ,
    502  ,
    527  ,    536  ,    545  ,    554  ,

    563  ,

    588  ,
    613  ,    622  ,    631  ,    640  ,
    649  ,
    674  ,    683  ,    692  ,    701  ,
    710  ,
    735  ,    744  ,    753  ,    762  ,
    771  ,
    796  ,    805  ,    814  ,    823  ,

    832  ,

    857  ,
    882  ,    891  ,    900  ,    909  ,
    918  ,
    943  ,    952  ,    961  ,    970  ,
    979  ,
    1004 ,    1013 ,    1022 ,    1031 ,
    1040 ,
    1065 ,    1074 ,    1083 ,    1092
};
#endif
#define ADD_CU_STOP_SPLIT             0   // Take into account & Stop Splitting
#define ADD_CU_CONTINUE_SPLIT         1   // Take into account & Continue Splitting
#define DO_NOT_ADD_CU_CONTINUE_SPLIT  2   // Do not take into account & Continue Splitting

#define DEPTH_64                      0   // Depth corresponding to the CU size
#define DEPTH_32                      1   // Depth corresponding to the CU size
#define DEPTH_16                      2   // Depth corresponding to the CU size
#define DEPTH_8                       3   // Depth corresponding to the CU size

static const uint8_t parentCuIndex[85] =
{
    0,
    0, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    21, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    42, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
    36, 0, 0, 1, 2, 3, 5, 0, 1, 2, 3, 10, 0, 1, 2, 3, 15, 0, 1, 2, 3,
};

const uint8_t incrementalCount[85] = {

    //64x64
    0,
    //32x32
    4, 4,
    4, 4,
    //16x16
    0, 0, 0, 0,
    0, 4, 0, 4,
    0, 0, 0, 0,
    0, 4, 0, 4,
    //8x8
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 4, 0, 0, 0, 4,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 4, 0, 0, 0, 4

};

/*******************************************
mdcSetDepth : set depth to be tested
*******************************************/
#define REFINEMENT_P        0x01
#define REFINEMENT_Pp1      0x02
#define REFINEMENT_Pp2      0x04
#define REFINEMENT_Pp3      0x08
#define REFINEMENT_Pm1      0x10
#define REFINEMENT_Pm2      0x20
#define REFINEMENT_Pm3      0x40



EbErrorType MdcRefinement(
    MdcpLocalCodingUnit_t                   *localCuArray,
    uint32_t                                  cu_index,
    uint32_t                                  depth,
    uint8_t                                   refinementLevel,
    uint8_t                                   lowestLevel)
{
    EbErrorType return_error = EB_ErrorNone;


    if (refinementLevel & REFINEMENT_P) {
        if (lowestLevel == REFINEMENT_P) {
            localCuArray[cu_index].stopSplit = EB_TRUE;
        }

    }
    else {
        localCuArray[cu_index].slectedCu = EB_FALSE;
    }

    if (refinementLevel & REFINEMENT_Pp1) {

        if (depth < 3 && cu_index < 81) {
            localCuArray[cu_index + 1].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + DepthOffset[depth + 1]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1]].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pp1) {
            if (depth < 3 && cu_index < 81) {
                localCuArray[cu_index + 1].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + DepthOffset[depth + 1]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1]].stopSplit = EB_TRUE;
            }
        }
    }

    if (refinementLevel & REFINEMENT_Pp2) {
        if (depth < 2 && cu_index < 65) {
            localCuArray[cu_index + 1 + 1].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 1 + DepthOffset[depth + 2]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 1 + 2 * DepthOffset[depth + 2]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 1 + 3 * DepthOffset[depth + 2]].slectedCu = EB_TRUE;

            localCuArray[cu_index + 1 + DepthOffset[depth + 1] + 1].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + DepthOffset[depth + 1] + 1 + DepthOffset[depth + 2]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + DepthOffset[depth + 1] + 1 + 2 * DepthOffset[depth + 2]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + DepthOffset[depth + 1] + 1 + 3 * DepthOffset[depth + 2]].slectedCu = EB_TRUE;

            localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1] + 1].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1] + 1 + DepthOffset[depth + 2]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1] + 1 + 2 * DepthOffset[depth + 2]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1] + 1 + 3 * DepthOffset[depth + 2]].slectedCu = EB_TRUE;

            localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1] + 1].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1] + 1 + DepthOffset[depth + 2]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1] + 1 + 2 * DepthOffset[depth + 2]].slectedCu = EB_TRUE;
            localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1] + 1 + 3 * DepthOffset[depth + 2]].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pp2) {
            if (depth < 2 && cu_index < 65) {
                localCuArray[cu_index + 1 + 1].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 1 + DepthOffset[depth + 2]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 1 + 2 * DepthOffset[depth + 2]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 1 + 3 * DepthOffset[depth + 2]].stopSplit = EB_TRUE;

                localCuArray[cu_index + 1 + DepthOffset[depth + 1] + 1].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + DepthOffset[depth + 1] + 1 + DepthOffset[depth + 2]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + DepthOffset[depth + 1] + 1 + 2 * DepthOffset[depth + 2]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + DepthOffset[depth + 1] + 1 + 3 * DepthOffset[depth + 2]].stopSplit = EB_TRUE;

                localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1] + 1].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1] + 1 + DepthOffset[depth + 2]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1] + 1 + 2 * DepthOffset[depth + 2]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 2 * DepthOffset[depth + 1] + 1 + 3 * DepthOffset[depth + 2]].stopSplit = EB_TRUE;

                localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1] + 1].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1] + 1 + DepthOffset[depth + 2]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1] + 1 + 2 * DepthOffset[depth + 2]].stopSplit = EB_TRUE;
                localCuArray[cu_index + 1 + 3 * DepthOffset[depth + 1] + 1 + 3 * DepthOffset[depth + 2]].stopSplit = EB_TRUE;
            }
        }
    }

    if (refinementLevel & REFINEMENT_Pp3) {
        uint8_t inLoop;
        uint8_t outLoop;
        uint8_t cu_index = 2;
        if (depth == 0) {

            for (outLoop = 0; outLoop < 16; ++outLoop) {
                for (inLoop = 0; inLoop < 4; ++inLoop) {
                    localCuArray[++cu_index].slectedCu = EB_TRUE;

                }
                cu_index += cu_index == 21 ? 2 : cu_index == 42 ? 2 : cu_index == 63 ? 2 : 1;

            }
            if (lowestLevel == REFINEMENT_Pp3) {
                cu_index = 2;
                for (outLoop = 0; outLoop < 16; ++outLoop) {
                    for (inLoop = 0; inLoop < 4; ++inLoop) {
                        localCuArray[++cu_index].stopSplit = EB_TRUE;
                    }
                    cu_index += cu_index == 21 ? 2 : cu_index == 42 ? 2 : cu_index == 63 ? 2 : 1;
                }
            }
        }

    }

    if (refinementLevel & REFINEMENT_Pm1) {
        if (depth > 0) {
            localCuArray[cu_index - 1 - parentCuIndex[cu_index]].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pm1) {
            if (depth > 0) {
                localCuArray[cu_index - 1 - parentCuIndex[cu_index]].stopSplit = EB_TRUE;
            }
        }
    }

    if (refinementLevel & REFINEMENT_Pm2) {
        if (depth == 2) {
            localCuArray[0].slectedCu = EB_TRUE;
        }
        if (depth == 3) {
            localCuArray[1].slectedCu = EB_TRUE;
            localCuArray[22].slectedCu = EB_TRUE;
            localCuArray[43].slectedCu = EB_TRUE;
            localCuArray[64].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pm2) {
            if (depth == 2) {
                localCuArray[0].stopSplit = EB_TRUE;
            }
            if (depth == 3) {
                localCuArray[1].stopSplit = EB_TRUE;
                localCuArray[22].stopSplit = EB_TRUE;
                localCuArray[43].stopSplit = EB_TRUE;
                localCuArray[64].stopSplit = EB_TRUE;
            }
        }
    }

    if (refinementLevel & REFINEMENT_Pm3) {
        if (depth == 3) {
            localCuArray[0].slectedCu = EB_TRUE;
        }
        if (lowestLevel == REFINEMENT_Pm2) {
            if (depth == 3) {
                localCuArray[0].stopSplit = EB_TRUE;
            }
        }
    }

    return return_error;
}
/*******************************************
Cost Computation for intra CU
*******************************************/
// TO be updated with AV1 functions
EbErrorType MdcIntraCuRate(
    uint32_t                                  cu_depth,
    MdRateEstimationContext_t               *md_rate_estimation_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    uint64_t                                  *mdcIntraRate)
{
    EbErrorType return_error = EB_ErrorNone;

    (void)mdcIntraRate;
    (void)picture_control_set_ptr;
    (void)md_rate_estimation_ptr;
    (void)cu_depth;
    //uint32_t chroma_mode = EB_INTRA_CHROMA_DM;
    //EB_PART_MODE partitionMode = SIZE_2Nx2N;
    //int32_t prediction_index = -1;

    //// Number of bits for each synatax element
    //uint64_t partSizeIntraBitsNum = 0;
    //uint64_t intraLumaModeBitsNum = 0;

    //// Luma and chroma rate
    //uint64_t lumaRate;
    //uint64_t chromaRate;

    //EncodeContext_t *encode_context_ptr = ((SequenceControlSet_t*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr))->encode_context_ptr;

    //CHECK_REPORT_ERROR(
    //    (partitionMode == SIZE_2Nx2N),
    //    encode_context_ptr->app_callback_ptr,
    //    EB_ENC_RD_COST_ERROR3);

    //// Estimate Chroma Mode Bits
    //chromaRate = md_rate_estimation_ptr->intraChromaBits[chroma_mode];

    //// Estimate Partition Size Bits :
    //// *Note - Intra is implicitly 2Nx2N
    //partSizeIntraBitsNum = ((uint8_t)cu_depth == (((SequenceControlSet_t *)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr)->max_sb_depth - 1)) ?
    //    md_rate_estimation_ptr->intraPartSizeBits[partitionMode] :
    //    ZERO_COST;

    //// Estimate Luma Mode Bits for Intra
    ///*intra_luma_mode_context(
    //cu_ptr,
    //luma_mode,
    //&prediction_index);*/

    //intraLumaModeBitsNum = (prediction_index != -1) ?
    //    md_rate_estimation_ptr->intraLumaBits[prediction_index] :
    //    md_rate_estimation_ptr->intraLumaBits[3];

    //// Rate of the candiadate mode is equal to the sum of the rate of independent syntax element
    //lumaRate =
    //    partSizeIntraBitsNum +
    //    intraLumaModeBitsNum;

    //// Assign fast cost
    //*(mdcIntraRate) = chromaRate + lumaRate;

    return return_error;
}

/*******************************************
Cost Computation for inter CU
*******************************************/
// TO be updated with AV1 functions
uint64_t MdcInterCuRate(
    uint32_t                                   direction,
    int16_t                                   xMvL0,
    int16_t                                   yMvL0,
    int16_t                                   xMvL1,
    int16_t                                   yMvL1)
{

    int32_t MVs_0, MVs_1, MVs_2, MVs_3;

    uint64_t rate = 86440;


    switch (direction)
    {
    default:
    case UNI_PRED_LIST_0:

        // Estimate the Motion Vector Prediction Index Bits
        rate += 23196;

        // Estimate the Motion Vector Difference Bits
        MVs_0 = ABS(xMvL0);
        MVs_1 = ABS(yMvL0);
        MVs_0 = MVs_0 > 499 ? 499 : MVs_0;
        MVs_1 = MVs_1 > 499 ? 499 : MVs_1;
        rate += mvBitTable[MVs_0][MVs_1];
        break;

    case UNI_PRED_LIST_1:

        // Estimate the Motion Vector Prediction Index Bits
        rate += 23196;

        // Estimate the Motion Vector Difference Bits

        MVs_2 = ABS(xMvL1);
        MVs_3 = ABS(yMvL1);
        MVs_2 = MVs_2 > 499 ? 499 : MVs_2;
        MVs_3 = MVs_3 > 499 ? 499 : MVs_3;


        rate += mvBitTable[MVs_2][MVs_3];




        break;

    case BI_PRED:
        //LIST 0 RATE ESTIMATION

        // Estimate the Motion Vector Prediction Index Bits

        rate += 46392;

        // Estimate the Motion Vector Difference Bits

        MVs_0 = ABS(xMvL0);
        MVs_1 = ABS(yMvL0);
        MVs_0 = MVs_0 > 499 ? 499 : MVs_0;
        MVs_1 = MVs_1 > 499 ? 499 : MVs_1;


        rate += mvBitTable[MVs_0][MVs_1];

        // Estimate the Motion Vector Difference Bits
        MVs_2 = ABS(xMvL1);
        MVs_3 = ABS(yMvL1);
        MVs_2 = MVs_2 > 499 ? 499 : MVs_2;
        MVs_3 = MVs_3 > 499 ? 499 : MVs_3;


        rate += mvBitTable[MVs_2][MVs_3];

        break;

    }

    return rate;
}


/*******************************************
Derive the contouring class
If (AC energy < 32 * 32) then apply aggressive action (Class 1),
else if (AC energy < 32 * 32 * 1.6) OR (32 * 32 * 3.5 < AC energy < 32 * 32 * 4.5 AND non-8x8) then moderate action (Class 2),
else no action
*******************************************/
uint8_t DeriveContouringClass(
    PictureParentControlSet_t   *parentPcsPtr,
    uint16_t                       sb_index,
    uint8_t                        leaf_index)
{
    uint8_t contouringClass = 0;

    SequenceControlSet_t *sequence_control_set_ptr = (SequenceControlSet_t*)parentPcsPtr->sequence_control_set_wrapper_ptr->object_ptr;

    if (parentPcsPtr->is_sb_homogeneous_over_time[sb_index]) {
        if (leaf_index > 0) {
            SbParams_t            *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            if (sb_params->is_edge_sb) {

                if (parentPcsPtr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_1) {
                    contouringClass = 2;
                }
                else if (parentPcsPtr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_2) {
                    contouringClass = 3;
                }
                else if (parentPcsPtr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < (ANTI_CONTOURING_TH_1 + ANTI_CONTOURING_TH_2)) {
                    contouringClass = 3;
                }
            }
            else {
                if (parentPcsPtr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_0) {
                    contouringClass = 1;
                }
                else if (parentPcsPtr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_1) {
                    contouringClass = 2;
                }
                else if (parentPcsPtr->sb_y_src_energy_cu_array[sb_index][(leaf_index - 1) / 21 + 1] < ANTI_CONTOURING_TH_2) {
                    contouringClass = 3;
                }
            }
        }
    }
    return(contouringClass);
}


void RefinementPredictionLoop(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    LargestCodingUnit_t                    *sb_ptr,
    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext_t     *context_ptr)
{

    MdcpLocalCodingUnit_t    *localCuArray         = context_ptr->localCuArray;
    SbParams_t               *sb_params            = &sequence_control_set_ptr->sb_params_array[sb_index];
    uint32_t                  temporal_layer_index = picture_control_set_ptr->temporal_layer_index;
    uint32_t                  cu_index             = 0;
#if !MDC_FIX_1
    uint8_t                   stationary_edge_over_time_flag = (&(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]))->stationary_edge_over_time_flag;
#endif
    sb_ptr->pred64 = EB_FALSE;
    while (cu_index < CU_MAX_COUNT)
    {
        if (sb_params->raster_scan_cu_validity[MD_SCAN_TO_RASTER_SCAN[cu_index]] && (localCuArray[cu_index].earlySplitFlag == EB_FALSE))
        {
            localCuArray[cu_index].slectedCu = EB_TRUE;
            sb_ptr->pred64 = (cu_index == 0) ? EB_TRUE : sb_ptr->pred64;
            uint32_t depth = GetCodedUnitStats(cu_index)->depth;
            uint8_t refinementLevel;
#if !MDC_FIX_1
            if (sb_ptr->picture_control_set_ptr->slice_type == I_SLICE) {

                {
                    uint8_t lowestLevel = 0x00;

                    if (sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE)
                        refinementLevel = NdpRefinementControl_ISLICE_M4[0][depth]; // HG: why always 0
                    else
                        refinementLevel = NdpRefinementControl_ISLICE_1080P_M4[0][depth]; // HG: why always 0

                    if (depth <= 1 && stationary_edge_over_time_flag > 0) {
                        if (depth == 0)
                            refinementLevel = Predp1 + Predp2 + Predp3;
                        else
                            refinementLevel = Pred + Predp1 + Predp2;

                    }
                    lowestLevel = (refinementLevel & REFINEMENT_Pp3) ? REFINEMENT_Pp3 : (refinementLevel & REFINEMENT_Pp2) ? REFINEMENT_Pp2 : (refinementLevel & REFINEMENT_Pp1) ? REFINEMENT_Pp1 :
                        (refinementLevel & REFINEMENT_P) ? REFINEMENT_P :
                        (refinementLevel & REFINEMENT_Pm1) ? REFINEMENT_Pm1 : (refinementLevel & REFINEMENT_Pm2) ? REFINEMENT_Pm2 : (refinementLevel & REFINEMENT_Pm3) ? REFINEMENT_Pm3 : 0x00;

                    MdcRefinement(
                        &(*context_ptr->localCuArray),
                        cu_index,
                        depth,
                        refinementLevel,
                        lowestLevel);
                }
            }
            else 
#endif    
            {
#if ADAPTIVE_DEPTH_PARTITIONING
                if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && (picture_control_set_ptr->parent_pcs_ptr->sb_depth_mode_array[sb_index] == SB_PRED_OPEN_LOOP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_depth_mode_array[sb_index] == SB_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE)) {
#else
                if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && (picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_PRED_OPEN_LOOP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE)) {
#endif
                    refinementLevel = Pred;
                }
                else

#if MDC_FIX_1
                    refinementLevel = NdpRefinementControl[temporal_layer_index][depth];
#else
                    if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_OPEN_LOOP_DEPTH_MODE ||
#if ADAPTIVE_DEPTH_PARTITIONING
                        (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && (picture_control_set_ptr->parent_pcs_ptr->sb_depth_mode_array[sb_index] == SB_OPEN_LOOP_DEPTH_MODE)))
#else
                        (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && (picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_OPEN_LOOP_DEPTH_MODE || picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_AVC_DEPTH_MODE)))
#endif
                        refinementLevel = NdpRefinementControlNREF[temporal_layer_index][depth];
                    else
                        refinementLevel = NdpRefinementControl_FAST[temporal_layer_index][depth];
#endif
                if (picture_control_set_ptr->parent_pcs_ptr->cu8x8_mode == CU_8x8_MODE_1) {
                    refinementLevel = ((refinementLevel & REFINEMENT_Pp1) && depth == 2) ? refinementLevel - REFINEMENT_Pp1 :
                        ((refinementLevel & REFINEMENT_Pp2) && depth == 1) ? refinementLevel - REFINEMENT_Pp2 :
                        ((refinementLevel & REFINEMENT_Pp3) && depth == 0) ? refinementLevel - REFINEMENT_Pp3 : refinementLevel;
                }

                uint8_t lowestLevel = 0x00;

                lowestLevel = (refinementLevel & REFINEMENT_Pp3) ? REFINEMENT_Pp3 : (refinementLevel & REFINEMENT_Pp2) ? REFINEMENT_Pp2 : (refinementLevel & REFINEMENT_Pp1) ? REFINEMENT_Pp1 :
                    (refinementLevel & REFINEMENT_P) ? REFINEMENT_P :
                    (refinementLevel & REFINEMENT_Pm1) ? REFINEMENT_Pm1 : (refinementLevel & REFINEMENT_Pm2) ? REFINEMENT_Pm2 : (refinementLevel & REFINEMENT_Pm3) ? REFINEMENT_Pm3 : 0x00;

                MdcRefinement(
                    &(*context_ptr->localCuArray),
                    cu_index,
                    depth,
                    refinementLevel,
                    lowestLevel);
            }

            cu_index += DepthOffset[depth];

        }
        else {

            cu_index++;
        }
    } // End while 1 CU Loop
}


void PrePredictionRefinement(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    LargestCodingUnit_t                    *sb_ptr,
    uint32_t                                  sb_index,
    uint32_t                                 *startDepth,
    uint32_t                                 *endDepth
)
{
    SbParams_t    *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

    EB_SLICE        slice_type = picture_control_set_ptr->slice_type;

    uint8_t           edge_block_num = picture_control_set_ptr->parent_pcs_ptr->edge_results_ptr[sb_index].edge_block_num;

    SbStat_t      *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);
    uint8_t           stationary_edge_over_time_flag = sb_stat_ptr->stationary_edge_over_time_flag;

    uint8_t           aura_status_iii = sb_ptr->aura_status_iii;






    if (picture_control_set_ptr->parent_pcs_ptr->high_dark_low_light_area_density_flag && picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0 && picture_control_set_ptr->parent_pcs_ptr->sharp_edge_sb_flag[sb_index] && !picture_control_set_ptr->parent_pcs_ptr->similar_colocated_sb_array_ii[sb_index]) {
        *startDepth = DEPTH_16;
    }

    if ((slice_type != I_SLICE && picture_control_set_ptr->high_intra_slection == 0) && (sb_params->is_complete_sb)) {

        if (picture_control_set_ptr->scene_caracteristic_id == EB_FRAME_CARAC_0) {

            if (picture_control_set_ptr->parent_pcs_ptr->grass_percentage_in_picture > 60 && aura_status_iii) {
                *startDepth = DEPTH_16;
            }
        }
    }

    if (picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag && edge_block_num)
    {
        *startDepth = DEPTH_16;
    }


    // S-LOGO

    if (stationary_edge_over_time_flag > 0) {

        *startDepth = DEPTH_16;
        *endDepth = DEPTH_16;

    }

    if (picture_control_set_ptr->parent_pcs_ptr->complex_sb_array[sb_ptr->index] == SB_COMPLEXITY_STATUS_2) {
        *startDepth = DEPTH_16;
    }
}


void ForwardCuToModeDecision(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,

    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext_t     *context_ptr
)
{

    uint8_t                   cu_index = 0;
    uint32_t                  cuClass = DO_NOT_ADD_CU_CONTINUE_SPLIT;
    EbBool                 split_flag = EB_TRUE;
    MdcLcuData_t           *resultsPtr = &picture_control_set_ptr->mdc_sb_array[sb_index];
    SbParams_t            *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    MdcpLocalCodingUnit_t  *localCuArray = context_ptr->localCuArray;
    EB_SLICE                slice_type = picture_control_set_ptr->slice_type;


    // CU Loop
    const CodedUnitStats_t *cuStatsPtr = GetCodedUnitStats(0);

    SbStat_t *sb_stat_ptr = &(picture_control_set_ptr->parent_pcs_ptr->sb_stat_array[sb_index]);

    EbBool    testAllDepthIntraSliceFlag = EB_FALSE;


    testAllDepthIntraSliceFlag = slice_type == I_SLICE &&
        (sb_stat_ptr->stationary_edge_over_time_flag || picture_control_set_ptr->parent_pcs_ptr->logo_pic_flag ||
        (picture_control_set_ptr->parent_pcs_ptr->very_low_var_pic_flag && picture_control_set_ptr->parent_pcs_ptr->low_motion_content_flag)) ?
        EB_TRUE : testAllDepthIntraSliceFlag;


    resultsPtr->leaf_count = 0;
#if OPEN_LOOP_EARLY_PARTITION
    uint8_t   enable_blk_4x4 = 0;
#endif
    cu_index = 0;

    while (cu_index < CU_MAX_COUNT)
    {
        split_flag = EB_TRUE;
        if (sb_params->raster_scan_cu_validity[MD_SCAN_TO_RASTER_SCAN[cu_index]])
        {
            cuStatsPtr = GetCodedUnitStats(cu_index);

            switch (cuStatsPtr->depth) {

            case 0:
            case 1:
            case 2:

                cuClass = DO_NOT_ADD_CU_CONTINUE_SPLIT;


                if (slice_type == I_SLICE) {
                    if (testAllDepthIntraSliceFlag) {
                        cuClass = ADD_CU_CONTINUE_SPLIT;
                    }
                    else {
                        cuClass = localCuArray[cu_index].slectedCu == EB_TRUE ? ADD_CU_CONTINUE_SPLIT : cuClass;
                        cuClass = localCuArray[cu_index].stopSplit == EB_TRUE ? ADD_CU_STOP_SPLIT : cuClass;
                    }
                }
                else {
                    cuClass = localCuArray[cu_index].slectedCu == EB_TRUE ? ADD_CU_CONTINUE_SPLIT : cuClass;
                    cuClass = localCuArray[cu_index].stopSplit == EB_TRUE ? ADD_CU_STOP_SPLIT : cuClass;

                }

                // Take into account MAX CU size & MAX intra size (from the API)
                cuClass = (cuStatsPtr->size > sequence_control_set_ptr->max_cu_size || (slice_type == I_SLICE && cuStatsPtr->size > sequence_control_set_ptr->max_intra_size)) ?
                    DO_NOT_ADD_CU_CONTINUE_SPLIT :
                    cuClass;

                // Take into account MIN CU size & Min intra size(from the API)
                cuClass = (cuStatsPtr->size == sequence_control_set_ptr->min_cu_size || (slice_type == I_SLICE && cuStatsPtr->size == sequence_control_set_ptr->min_intra_size)) ?
                    ADD_CU_STOP_SPLIT :
                    cuClass;

                switch (cuClass) {

                case ADD_CU_STOP_SPLIT:
                    // Stop
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
#if OPEN_LOOP_EARLY_PARTITION
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = pa_to_ep_block_index[cu_index];
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = 1;
#endif
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;

                    break;

                case ADD_CU_CONTINUE_SPLIT:
                    // Go Down + consider the current CU as candidate
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
#if OPEN_LOOP_EARLY_PARTITION
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = pa_to_ep_block_index[cu_index];
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = 1;
#endif
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;

                case DO_NOT_ADD_CU_CONTINUE_SPLIT:
                    // Go Down + do not consider the current CU as candidate
                    split_flag = EB_TRUE;

                    break;

                default:
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
#if OPEN_LOOP_EARLY_PARTITION
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = pa_to_ep_block_index[cu_index];
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = 1;
#endif
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    break;
                }

                break;
            case 3:

                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
#if OPEN_LOOP_EARLY_PARTITION
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = pa_to_ep_block_index[cu_index];
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = 1;

                if (enable_blk_4x4) {
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;

                    int first_4_index = pa_to_ep_block_index[cu_index] + d1_depth_offset[sequence_control_set_ptr->sb_size == BLOCK_128X128][cuStatsPtr->depth];
                    for (int i = 0; i < 4; ++i) {

                        resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;

                        resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = first_4_index + i;
                        resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = 1;

                        resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;
                    }
                }else
                    resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;
#else
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_FALSE;
#endif

                break;

            default:
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].leaf_index = cu_index;
#if OPEN_LOOP_EARLY_PARTITION
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].mds_idx = pa_to_ep_block_index[cu_index];
                resultsPtr->leaf_data_array[resultsPtr->leaf_count].tot_d1_blocks = 1;
#endif
                resultsPtr->leaf_data_array[resultsPtr->leaf_count++].split_flag = split_flag = EB_TRUE;
                break;
            }
        }

        cu_index += (split_flag == EB_TRUE) ? 1 : DepthOffset[cuStatsPtr->depth];

    } // End CU Loop

}



void MdcInterDepthDecision(
    ModeDecisionConfigurationContext_t     *context_ptr,
    uint32_t                                 origin_x,
    uint32_t                                 origin_y,
    uint32_t                                 endDepth,
#if !MDC_FIX_1
    uint32_t                                 splitFlagBits0,
    uint32_t                                 splitFlagBits1,
    uint64_t                                 lambda,
#endif
    uint32_t                                 cu_index)
{

    uint32_t               leftCuIndex;
    uint32_t               topCuIndex;
    uint32_t               topLeftCuIndex;
    uint32_t               depthZeroCandidateCuIndex;
    uint32_t               depthOneCandidateCuIndex = cu_index;
    uint32_t               depthTwoCandidateCuIndex = cu_index;
    uint64_t               depthNRate = 0;
    uint64_t               depthNPlusOneRate = 0;
    uint64_t               depthNCost = 0;
    uint64_t               depthNPlusOneCost = 0;
    MdcpLocalCodingUnit_t *localCuArray = context_ptr->localCuArray;
    /*** Stage 0: Inter depth decision: depth 2 vs depth 3 ***/
    // Walks to the last coded 8x8 block for merging
    uint8_t  group_of8x8_blocks_count = context_ptr->group_of8x8_blocks_count;
    uint8_t  group_of16x16_blocks_count = context_ptr->group_of16x16_blocks_count;
    if ((GROUP_OF_4_8x8_BLOCKS(origin_x, origin_y))) {

        group_of8x8_blocks_count++;

        // From the last coded cu index, get the indices of the left, top, and top left cus
        leftCuIndex = cu_index - DEPTH_THREE_STEP;
        topCuIndex = leftCuIndex - DEPTH_THREE_STEP;
        topLeftCuIndex = topCuIndex - DEPTH_THREE_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depthTwoCandidateCuIndex = topLeftCuIndex - 1;

        // Compute depth N cost
        localCuArray[depthTwoCandidateCuIndex].splitContext = 0;
#if !MDC_FIX_1
        depthNRate = (((lambda  * splitFlagBits0) + MD_OFFSET) >> MD_SHIFT);
#endif
        depthNCost = (localCuArray[depthTwoCandidateCuIndex]).earlyCost + depthNRate;

        if (endDepth < 3) {

            (localCuArray[depthTwoCandidateCuIndex]).earlySplitFlag = EB_FALSE;
            (localCuArray[depthTwoCandidateCuIndex]).earlyCost = depthNCost;

        }
        else {
#if !MDC_FIX_1
            // Assign rate
            depthNPlusOneRate = (((lambda  * splitFlagBits1) + MD_OFFSET) >> MD_SHIFT);
#endif
            depthNPlusOneCost = (localCuArray[cu_index]).earlyCost + (localCuArray[leftCuIndex]).earlyCost + (localCuArray[topCuIndex]).earlyCost + (localCuArray[topLeftCuIndex]).earlyCost + depthNPlusOneRate;

            if (depthNCost <= depthNPlusOneCost) {

                // If the cost is low enough to warrant not spliting further:
                // 1. set the split flag of the candidate pu for merging to false
                // 2. update the last pu index
                (localCuArray[depthTwoCandidateCuIndex]).earlySplitFlag = EB_FALSE;
                (localCuArray[depthTwoCandidateCuIndex]).earlyCost = depthNCost;

            }
            else {
                // If the cost is not low enough:
                // update the cost of the candidate pu for merging
                // this update is required for the next inter depth decision
                (&localCuArray[depthTwoCandidateCuIndex])->earlyCost = depthNPlusOneCost;
            }

        }
    }

    // Walks to the last coded 16x16 block for merging
    if (GROUP_OF_4_16x16_BLOCKS(GetCodedUnitStats(depthTwoCandidateCuIndex)->origin_x, GetCodedUnitStats(depthTwoCandidateCuIndex)->origin_y) &&
        (group_of8x8_blocks_count == 4)) {

        group_of8x8_blocks_count = 0;
        group_of16x16_blocks_count++;

        // From the last coded pu index, get the indices of the left, top, and top left pus
        leftCuIndex = depthTwoCandidateCuIndex - DEPTH_TWO_STEP;
        topCuIndex = leftCuIndex - DEPTH_TWO_STEP;
        topLeftCuIndex = topCuIndex - DEPTH_TWO_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depthOneCandidateCuIndex = topLeftCuIndex - 1;

        if (GetCodedUnitStats(depthOneCandidateCuIndex)->depth == 1) {
#if !MDC_FIX_1
            depthNRate = (((lambda  *splitFlagBits0) + MD_OFFSET) >> MD_SHIFT);
#endif
            depthNCost = localCuArray[depthOneCandidateCuIndex].earlyCost + depthNRate;
            if (endDepth < 2) {

                localCuArray[depthOneCandidateCuIndex].earlySplitFlag = EB_FALSE;
                localCuArray[depthOneCandidateCuIndex].earlyCost = depthNCost;

            }
            else {
                // Compute depth N+1 cost
#if !MDC_FIX_1
                // Assign rate
                depthNPlusOneRate = (((lambda  *splitFlagBits1) + MD_OFFSET) >> MD_SHIFT);
#endif
                depthNPlusOneCost = localCuArray[depthTwoCandidateCuIndex].earlyCost +
                    localCuArray[leftCuIndex].earlyCost +
                    localCuArray[topCuIndex].earlyCost +
                    localCuArray[topLeftCuIndex].earlyCost +
                    depthNPlusOneRate;

                // Inter depth comparison: depth 1 vs depth 2
                if (depthNCost <= depthNPlusOneCost) {
                    // If the cost is low enough to warrant not spliting further:
                    // 1. set the split flag of the candidate pu for merging to false
                    // 2. update the last pu index
                    localCuArray[depthOneCandidateCuIndex].earlySplitFlag = EB_FALSE;
                    localCuArray[depthOneCandidateCuIndex].earlyCost = depthNCost;
                }
                else {
                    // If the cost is not low enough:
                    // update the cost of the candidate pu for merging
                    // this update is required for the next inter depth decision
                    localCuArray[depthOneCandidateCuIndex].earlyCost = depthNPlusOneCost;
                }
            }
        }
    }

    // Stage 2: Inter depth decision: depth 0 vs depth 1

    // Walks to the last coded 32x32 block for merging
    // Stage 2 isn't performed in I slices since the abcense of 64x64 candidates
    if (GROUP_OF_4_32x32_BLOCKS(GetCodedUnitStats(depthOneCandidateCuIndex)->origin_x, GetCodedUnitStats(depthOneCandidateCuIndex)->origin_y) &&
        (group_of16x16_blocks_count == 4)) {

        group_of16x16_blocks_count = 0;

        // From the last coded pu index, get the indices of the left, top, and top left pus
        leftCuIndex = depthOneCandidateCuIndex - DEPTH_ONE_STEP;
        topCuIndex = leftCuIndex - DEPTH_ONE_STEP;
        topLeftCuIndex = topCuIndex - DEPTH_ONE_STEP;

        // From the top left index, get the index of the candidate pu for merging
        depthZeroCandidateCuIndex = topLeftCuIndex - 1;

        if (GetCodedUnitStats(depthZeroCandidateCuIndex)->depth == 0) {

            // Compute depth N cost
#if !MDC_FIX_1
            depthNRate = (((lambda  *splitFlagBits0) + MD_OFFSET) >> MD_SHIFT);
#endif
            depthNCost = (&localCuArray[depthZeroCandidateCuIndex])->earlyCost + depthNRate;
            if (endDepth < 1) {

                (&localCuArray[depthZeroCandidateCuIndex])->earlySplitFlag = EB_FALSE;

            }
            else {
                // Compute depth N+1 cost
#if !MDC_FIX_1
                // Assign rate
                depthNPlusOneRate = (((lambda  *splitFlagBits1) + MD_OFFSET) >> MD_SHIFT);
#endif
                depthNPlusOneCost = localCuArray[depthOneCandidateCuIndex].earlyCost +
                    localCuArray[leftCuIndex].earlyCost +
                    localCuArray[topCuIndex].earlyCost +
                    localCuArray[topLeftCuIndex].earlyCost +
                    depthNPlusOneRate;

                // Inter depth comparison: depth 0 vs depth 1
                if (depthNCost <= depthNPlusOneCost) {
                    // If the cost is low enough to warrant not spliting further:
                    // 1. set the split flag of the candidate pu for merging to false
                    // 2. update the last pu index
                    (&localCuArray[depthZeroCandidateCuIndex])->earlySplitFlag = EB_FALSE;
                }
            }
        }
    }

    context_ptr->group_of8x8_blocks_count = group_of8x8_blocks_count;
    context_ptr->group_of16x16_blocks_count = group_of16x16_blocks_count;
}

void PredictionPartitionLoop(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    uint32_t                                sb_index,
    uint32_t                                tbOriginX,
    uint32_t                                tbOriginY,
    uint32_t                                startDepth,
    uint32_t                                endDepth,
    ModeDecisionConfigurationContext_t     *context_ptr){

#if !OPEN_LOOP_EARLY_PARTITION
    MdRateEstimationContext_t *md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
#endif
    MdcpLocalCodingUnit_t *localCuArray = context_ptr->localCuArray;
    MdcpLocalCodingUnit_t   *cu_ptr;

    SbParams_t *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
#if !MDC_FIX_0
    uint32_t      cuInterSad = 0;
    uint64_t      cuInterRate = 0;
    uint64_t      cuInterCost = 0;
#endif
#if !OPEN_LOOP_EARLY_PARTITION
    uint32_t      cuIntraSad = 0;
    uint64_t      cuIntraRate = 0;
    uint64_t      cuIntraCost = 0;
#endif
    uint32_t      cuIndexInRaterScan;
    uint32_t      cu_index = 0;
    uint32_t      startIndex = 0;

    (void)tbOriginX;
    (void)tbOriginY;

    const CodedUnitStats_t *cuStatsPtr;

    for (cu_index = startIndex; cu_index < CU_MAX_COUNT; ++cu_index)

    {

        localCuArray[cu_index].slectedCu = EB_FALSE;
        localCuArray[cu_index].stopSplit = EB_FALSE;

        cu_ptr = &localCuArray[cu_index];
        cuIndexInRaterScan = MD_SCAN_TO_RASTER_SCAN[cu_index];
        if (sb_params->raster_scan_cu_validity[cuIndexInRaterScan])
        {
            uint32_t depth;
#if !OPEN_LOOP_EARLY_PARTITION
            uint32_t size;
#endif
            cuStatsPtr = GetCodedUnitStats(cu_index);

            depth = cuStatsPtr->depth;
#if !OPEN_LOOP_EARLY_PARTITION
            size = cuStatsPtr->size;
#endif
            cu_ptr->earlySplitFlag = (depth < endDepth) ? EB_TRUE : EB_FALSE;

            if (depth >= startDepth && depth <= endDepth) {
                //reset the flags here:   all CU splitFalg=TRUE. default: we always split. interDepthDecision will select where  to stop splitting(ie setting the flag to False)
#if !OPEN_LOOP_EARLY_PARTITION
                {
                    MdcIntraCuRate(
                        depth,
                        md_rate_estimation_ptr,
                        picture_control_set_ptr,
                        &cuIntraRate);


                    OisCu32Cu16Results_t            *oisCu32Cu16ResultsPtr = picture_control_set_ptr->parent_pcs_ptr->ois_cu32_cu16_results[sb_index];
                    OisCu8Results_t                   *oisCu8ResultsPtr = picture_control_set_ptr->parent_pcs_ptr->ois_cu8_results[sb_index];
                    OisCandidate_t * OisCuPtr = cuIndexInRaterScan < RASTER_SCAN_CU_INDEX_8x8_0 ?
                        oisCu32Cu16ResultsPtr->sorted_ois_candidate[cuIndexInRaterScan] : oisCu8ResultsPtr->sorted_ois_candidate[cuIndexInRaterScan - RASTER_SCAN_CU_INDEX_8x8_0];


                    if (size > 32) {
                        cuIntraSad =
                            oisCu32Cu16ResultsPtr->sorted_ois_candidate[1][0].distortion +
                            oisCu32Cu16ResultsPtr->sorted_ois_candidate[2][0].distortion +
                            oisCu32Cu16ResultsPtr->sorted_ois_candidate[3][0].distortion +
                            oisCu32Cu16ResultsPtr->sorted_ois_candidate[4][0].distortion;
                    }
                    else if (size == 32) {
                        cuIntraSad = oisCu32Cu16ResultsPtr->sorted_ois_candidate[cuIndexInRaterScan][0].distortion;
                    }
                    else {
                        if (size > 8) {
                            cuIntraSad = oisCu32Cu16ResultsPtr->sorted_ois_candidate[cuIndexInRaterScan][0].distortion;
                        }
                        else {
                            if (OisCuPtr[0].valid_distortion) {
                                cuIntraSad = OisCuPtr[0].distortion;
                            }
                            else {

                                const CodedUnitStats_t  *cu_stats = GetCodedUnitStats(ParentBlockIndex[cu_index]);
                                const uint32_t me2Nx2NTableOffset = cu_stats->cuNumInDepth + me2Nx2NOffset[cu_stats->depth];


                                if (oisCu8ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].valid_distortion) {
                                    cuIntraSad = oisCu8ResultsPtr->sorted_ois_candidate[me2Nx2NTableOffset][0].distortion;
                                }
                                else {
                                    cuIntraSad = 0;
                                }
                            }
                        }
                    }


                    cuIntraCost = (cuIntraSad << COST_PRECISION) + ((context_ptr->lambda * cuIntraRate + MD_OFFSET) >> MD_SHIFT);
                    cu_ptr->earlyCost = cuIntraCost;

                }
#endif
                if (picture_control_set_ptr->slice_type != I_SLICE) {

                    MeCuResults_t * mePuResult = &picture_control_set_ptr->parent_pcs_ptr->me_results[sb_index][cuIndexInRaterScan];
#if MDC_FIX_0            
                    // Initialize the mdc candidate (only av1 rate estimation inputs)
                    // Hsan: mode, direction, .. could be modified toward better early inter depth decision (e.g. NEARESTMV instead of NEWMV)
                    context_ptr->mdc_candidate_ptr->md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
                    context_ptr->mdc_candidate_ptr->type = INTER_MODE;
                    context_ptr->mdc_candidate_ptr->merge_flag = EB_FALSE;
                    context_ptr->mdc_candidate_ptr->merge_index = 0;
#if MDC_FIX_1
                    context_ptr->mdc_candidate_ptr->prediction_direction[0] = (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0) ?
                        UNI_PRED_LIST_0 :
                        mePuResult->distortionDirection[0].direction;
#else
                    context_ptr->mdc_candidate_ptr->prediction_direction[0] = UNI_PRED_LIST_0;
#endif
                    context_ptr->mdc_candidate_ptr->is_skip_mode_flag = 0;
                    // Hsan: what's the best mode for rate simulation
#if MDC_FIX_1
                    context_ptr->mdc_candidate_ptr->inter_mode = NEARESTMV;
                    context_ptr->mdc_candidate_ptr->pred_mode = NEARESTMV;
#else
                    context_ptr->mdc_candidate_ptr->inter_mode = NEWMV;
                    context_ptr->mdc_candidate_ptr->pred_mode = NEWMV;
#endif
                    context_ptr->mdc_candidate_ptr->motion_mode = SIMPLE_TRANSLATION;
#if !MDC_FIX_1
                    context_ptr->mdc_candidate_ptr->is_compound = 0;
#endif
                    context_ptr->mdc_candidate_ptr->is_new_mv = 1;
                    context_ptr->mdc_candidate_ptr->is_zero_mv = 0;
                    context_ptr->mdc_candidate_ptr->drl_index = 0;
                    context_ptr->mdc_candidate_ptr->motionVector_x_L0 = mePuResult->xMvL0 << 1;
                    context_ptr->mdc_candidate_ptr->motionVector_y_L0 = mePuResult->yMvL0 << 1;
#if MDC_FIX_1
                    context_ptr->mdc_candidate_ptr->motionVector_x_L1 = mePuResult->xMvL1 << 1;
                    context_ptr->mdc_candidate_ptr->motionVector_y_L1 = mePuResult->yMvL1 << 1;
#endif
                    context_ptr->mdc_candidate_ptr->ref_mv_index = 0;
                    context_ptr->mdc_candidate_ptr->pred_mv_weight = 0;
#if MDC_FIX_1
                    if (context_ptr->mdc_candidate_ptr->prediction_direction[0] == BI_PRED) {
                        context_ptr->mdc_candidate_ptr->ref_frame_type = LAST_BWD_FRAME;
                        context_ptr->mdc_candidate_ptr->is_compound = 1;
                    }
                    else if (context_ptr->mdc_candidate_ptr->prediction_direction[0] == UNI_PRED_LIST_0) {
                        context_ptr->mdc_candidate_ptr->ref_frame_type = LAST_FRAME;
                        context_ptr->mdc_candidate_ptr->is_compound = 0;
                    }
                    else { // context_ptr->mdc_candidate_ptr->prediction_direction[0]
                        context_ptr->mdc_candidate_ptr->ref_frame_type = BWDREF_FRAME;
                        context_ptr->mdc_candidate_ptr->is_compound = 0;
                    }
#else
                    context_ptr->mdc_candidate_ptr->ref_frame_type = LAST_FRAME;
#endif
                    context_ptr->mdc_candidate_ptr->motion_vector_pred_x[REF_LIST_0] = 0;
                    context_ptr->mdc_candidate_ptr->motion_vector_pred_y[REF_LIST_0] = 0;
                    // Initialize the ref mv
                    memset(context_ptr->mdc_ref_mv_stack,0,sizeof(CandidateMv));
                    context_ptr->blk_geom = get_blk_geom_mds(pa_to_ep_block_index[cu_index]);
                    // Initialize mdc cu (only av1 rate estimation inputs)
                    context_ptr->mdc_cu_ptr->is_inter_ctx = 0;
                    context_ptr->mdc_cu_ptr->skip_flag_context = 0;
#if MDC_FIX_1
                    context_ptr->mdc_cu_ptr->inter_mode_ctx[context_ptr->mdc_candidate_ptr->ref_frame_type] = 0;
                    context_ptr->mdc_cu_ptr->reference_mode_context = 0;
                    context_ptr->mdc_cu_ptr->compoud_reference_type_context = 0;
#endif
                    av1_zero(context_ptr->mdc_cu_ptr->av1xd->neighbors_ref_counts); // Hsan: neighbor not generated @ open loop partitioning => assumes always (0,0)

                    // Fast Cost Calc
                    cu_ptr->earlyCost = av1_inter_fast_cost(
                        context_ptr->mdc_cu_ptr,
                        context_ptr->mdc_candidate_ptr,
                        context_ptr->qp,
                        mePuResult->distortionDirection[0].distortion,
                        (uint64_t) 0,
                        context_ptr->lambda,
                        picture_control_set_ptr,
                        context_ptr->mdc_ref_mv_stack,
                        context_ptr->blk_geom,
                        (tbOriginY + context_ptr->blk_geom->origin_y) >> MI_SIZE_LOG2,
                        (tbOriginX + context_ptr->blk_geom->origin_x) >> MI_SIZE_LOG2,
                        DC_PRED,        // Hsan: neighbor not generated @ open loop partitioning
                        DC_PRED);       // Hsan: neighbor not generated @ open loop partitioning
#else
                    cuInterRate = MdcInterCuRate(
                        mePuResult->distortionDirection[0].direction,
                        mePuResult->xMvL0,
                        mePuResult->yMvL0,
                        mePuResult->xMvL1,
                        mePuResult->yMvL1);
                    cuInterSad = mePuResult->distortionDirection[0].distortion;


                    cuInterCost = (cuInterSad << COST_PRECISION) + ((context_ptr->lambda * cuInterRate + MD_OFFSET) >> MD_SHIFT);
                    cu_ptr->earlyCost = cuInterCost;
#endif
                }
#if !OPEN_LOOP_EARLY_PARTITION
                cu_ptr->earlyCost = picture_control_set_ptr->slice_type == I_SLICE ? cuIntraCost : MIN(cuInterCost, cuIntraCost);
#endif
                if (endDepth == 2) {
                    context_ptr->group_of8x8_blocks_count = depth == 2 ? incrementalCount[cuIndexInRaterScan] : 0;
                }

                if (endDepth == 1) {
                    context_ptr->group_of16x16_blocks_count = depth == 1 ? incrementalCount[cuIndexInRaterScan] : 0;
                }
                MdcInterDepthDecision(
                    context_ptr,
                    cuStatsPtr->origin_x,
                    cuStatsPtr->origin_y,
                    endDepth,
#if !MDC_FIX_1
                    md_rate_estimation_ptr->splitFlagBits[0],
                    md_rate_estimation_ptr->splitFlagBits[3],
                    context_ptr->lambda,
#endif
                    cu_index);
            }
            else {
                cu_ptr->earlyCost = ~0u;
            }

        }

    }// End CU Loop

}

EbErrorType EarlyModeDecisionLcu(
    SequenceControlSet_t                   *sequence_control_set_ptr,
    PictureControlSet_t                    *picture_control_set_ptr,
    LargestCodingUnit_t                    *sb_ptr,
    uint32_t                                  sb_index,
    ModeDecisionConfigurationContext_t     *context_ptr){

    EbErrorType    return_error = EB_ErrorNone;
    uint32_t       tbOriginX    = sb_ptr->origin_x;
    uint32_t       tbOriginY    = sb_ptr->origin_y;
#if !OPEN_LOOP_EARLY_PARTITION  
    EB_SLICE       slice_type   = picture_control_set_ptr->slice_type;
#endif

#if ADAPTIVE_DEPTH_PARTITIONING
    uint32_t      startDepth = DEPTH_64;
#else
    uint32_t      startDepth = (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode == PIC_SB_SWITCH_DEPTH_MODE && picture_control_set_ptr->parent_pcs_ptr->sb_md_mode_array[sb_index] == LCU_AVC_DEPTH_MODE) ?
        DEPTH_16 :
        DEPTH_64;
#endif

#if OPEN_LOOP_EARLY_PARTITION        
    uint32_t      endDepth =  DEPTH_8 ;
#else
    uint32_t      endDepth = (slice_type == I_SLICE) ? DEPTH_8 : DEPTH_16;
#endif
    context_ptr->group_of8x8_blocks_count = 0;
    context_ptr->group_of16x16_blocks_count = 0;

#if !OPEN_LOOP_EARLY_PARTITION
    // The MDC refinements had been taken into account at the budgeting algorithm & therefore could be skipped for the ADP case
    if (picture_control_set_ptr->parent_pcs_ptr->pic_depth_mode != PIC_SB_SWITCH_DEPTH_MODE) {
        PrePredictionRefinement(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_ptr,
            sb_index,
            &startDepth,
            &endDepth);
    }
#endif
    PredictionPartitionLoop(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sb_index,
        tbOriginX,
        tbOriginY,
        startDepth,
        endDepth,
        context_ptr
    );

    RefinementPredictionLoop(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sb_ptr,
        sb_index,
        context_ptr);

    ForwardCuToModeDecision(
        sequence_control_set_ptr,
        picture_control_set_ptr,

        sb_index,
        context_ptr);

    return return_error;

}



