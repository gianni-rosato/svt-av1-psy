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

/***************************************
* Includes
***************************************/
#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbSequenceControlSet.h"

#include "EbModeDecision.h"
#include "EbAdaptiveMotionVectorPrediction.h"
#include "EbTransformUnit.h"
#include "EbModeDecisionProcess.h"
#include "EbMotionEstimation.h"

#if ICOPY
#include "av1me.h"
#include "hash.h"
#endif

/********************************************
* Constants
********************************************/
//static uint32_t  AntiContouringIntraMode[11] = { EB_INTRA_PLANAR, EB_INTRA_DC, EB_INTRA_HORIZONTAL, EB_INTRA_VERTICAL,
//EB_INTRA_MODE_2, EB_INTRA_MODE_6, EB_INTRA_MODE_14, EB_INTRA_MODE_18, EB_INTRA_MODE_22, EB_INTRA_MODE_30, EB_INTRA_MODE_34 };

const uint32_t parentIndex[85] = { 0, 0, 0, 2, 2, 2, 2, 0, 7, 7, 7, 7, 0, 12, 12, 12, 12, 0, 17, 17, 17, 17, 0, 0,
23, 23, 23, 23, 0, 28, 28, 28, 28, 0, 33, 33, 33, 33, 0, 38, 38, 38, 38, 0, 0,
44, 44, 44, 44, 0, 49, 49, 49, 49, 0, 54, 54, 54, 54, 0, 59, 59, 59, 59, 0, 0,
65, 65, 65, 65, 0, 70, 70, 70, 70, 0, 75, 75, 75, 75, 0, 80, 80, 80, 80 };
#if !OIS_BASED_INTRA
uint8_t GetNumOfIntraModesFromOisPoint(
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                       meSad,
    uint32_t                       oisDcSad
);
#endif
extern uint32_t stage1ModesArray[];

uint8_t GetMaxDrlIndex(uint8_t  refmvCnt, PredictionMode   mode);
int32_t av1_mv_bit_cost(const MV *mv, const MV *ref, const int32_t *mvjcost,
    int32_t *mvcost[2], int32_t weight);
#define MV_COST_WEIGHT 108

void ChooseBestAv1MvPred(
    ModeDecisionContext_t            *context_ptr,
    struct MdRateEstimationContext_s      *md_rate_estimation_ptr,
    CodingUnit_t      *cu_ptr,
    MvReferenceFrame ref_frame,
    uint8_t              is_compound,
    PredictionMode    mode,              //NEW or NEW_NEW
    int16_t             mv0x,
    int16_t             mv0y,
    int16_t             mv1x,
    int16_t             mv1y,
    uint8_t             *bestDrlIndex,      // output
    IntMv             bestPredmv[2]      // output
)
{
    uint8_t              drli, maxDrlIndex;
    IntMv             nearestmv[2];
    IntMv             nearmv[2];
    IntMv             ref_mv[2];
    uint32_t             bestmvCost = 0xFFFFFFFF;
    MV                 mv;

    maxDrlIndex = GetMaxDrlIndex(cu_ptr->av1xd->ref_mv_count[ref_frame], mode);
    // maxDrlIndex = 1;

    for (drli = 0; drli < maxDrlIndex; drli++) {

        get_av1_mv_pred_drl(
            context_ptr,
            cu_ptr,
            ref_frame,
            is_compound,
            mode,
            drli,
            nearestmv,
            nearmv,
            ref_mv);

        //compute the rate for this drli Cand
        mv.row = mv0y;
        mv.col = mv0x;

        uint32_t mvRate = (uint32_t)av1_mv_bit_cost(
            &mv,
            &(ref_mv[0].as_mv),
            md_rate_estimation_ptr->nmv_vec_cost,
            md_rate_estimation_ptr->nmvcoststack,
            MV_COST_WEIGHT);

        if (is_compound) {

            mv.row = mv1y;
            mv.col = mv1x;

            mvRate += (uint32_t)av1_mv_bit_cost(
                &mv,
                &(ref_mv[1].as_mv),
                md_rate_estimation_ptr->nmv_vec_cost,
                md_rate_estimation_ptr->nmvcoststack,
                MV_COST_WEIGHT);
        }

        if (mvRate < bestmvCost) {

            bestmvCost = mvRate;
            *bestDrlIndex = drli;
            bestPredmv[0] = ref_mv[0];
            bestPredmv[1] = ref_mv[1];
        }

    }

}


/***************************************
* Mode Decision Candidate Ctor
***************************************/
#if INTRA_INTER_FAST_LOOP
EbErrorType mode_decision_candidate_buffer_ctor(
    ModeDecisionCandidateBuffer_t **buffer_dbl_ptr,
    uint64_t                       *fast_cost_ptr,
    uint64_t                       *full_cost_ptr,
    uint64_t                       *full_cost_skip_ptr,
    uint64_t                       *full_cost_merge_ptr)
#else
EbErrorType mode_decision_candidate_buffer_ctor(
    ModeDecisionCandidateBuffer_t **buffer_dbl_ptr,
    uint16_t                          sb_max_size,
    EB_BITDEPTH                     max_bitdepth,
    uint64_t                         *fast_cost_ptr,
    uint64_t                         *full_cost_ptr,
    uint64_t                         *full_cost_skip_ptr,
    uint64_t                         *full_cost_merge_ptr
)
#endif
{
    EbPictureBufferDescInitData_t pictureBufferDescInitData;
    EbPictureBufferDescInitData_t doubleWidthPictureBufferDescInitData;

    EbPictureBufferDescInitData_t ThirtyTwoWidthPictureBufferDescInitData;


    EbErrorType return_error = EB_ErrorNone;
    // Allocate Buffer
    ModeDecisionCandidateBuffer_t *bufferPtr;
    EB_MALLOC(ModeDecisionCandidateBuffer_t*, bufferPtr, sizeof(ModeDecisionCandidateBuffer_t), EB_N_PTR);
    *buffer_dbl_ptr = bufferPtr;

    // Init Picture Data
#if INTRA_INTER_FAST_LOOP
    pictureBufferDescInitData.maxWidth = MAX_SB_SIZE;
    pictureBufferDescInitData.maxHeight = MAX_SB_SIZE;
    pictureBufferDescInitData.bit_depth = EB_8BIT;
#else
    pictureBufferDescInitData.maxWidth = sb_max_size;
    pictureBufferDescInitData.maxHeight = sb_max_size;
    pictureBufferDescInitData.bit_depth = max_bitdepth;
#endif
    pictureBufferDescInitData.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;
    pictureBufferDescInitData.left_padding = 0;
    pictureBufferDescInitData.right_padding = 0;
    pictureBufferDescInitData.top_padding = 0;
    pictureBufferDescInitData.bot_padding = 0;
    pictureBufferDescInitData.splitMode = EB_FALSE;
#if INTRA_INTER_FAST_LOOP
    doubleWidthPictureBufferDescInitData.maxWidth = MAX_SB_SIZE;
    doubleWidthPictureBufferDescInitData.maxHeight = MAX_SB_SIZE;
#else
    doubleWidthPictureBufferDescInitData.maxWidth = sb_max_size;
    doubleWidthPictureBufferDescInitData.maxHeight = sb_max_size;
#endif
    doubleWidthPictureBufferDescInitData.bit_depth = EB_16BIT;
    doubleWidthPictureBufferDescInitData.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;
    doubleWidthPictureBufferDescInitData.left_padding = 0;
    doubleWidthPictureBufferDescInitData.right_padding = 0;
    doubleWidthPictureBufferDescInitData.top_padding = 0;
    doubleWidthPictureBufferDescInitData.bot_padding = 0;
    doubleWidthPictureBufferDescInitData.splitMode = EB_FALSE;

#if INTRA_INTER_FAST_LOOP
    ThirtyTwoWidthPictureBufferDescInitData.maxWidth = MAX_SB_SIZE;
    ThirtyTwoWidthPictureBufferDescInitData.maxHeight = MAX_SB_SIZE;
#else
    ThirtyTwoWidthPictureBufferDescInitData.maxWidth = sb_max_size;
    ThirtyTwoWidthPictureBufferDescInitData.maxHeight = sb_max_size;
#endif
    ThirtyTwoWidthPictureBufferDescInitData.bit_depth = EB_32BIT;
    ThirtyTwoWidthPictureBufferDescInitData.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;
    ThirtyTwoWidthPictureBufferDescInitData.left_padding = 0;
    ThirtyTwoWidthPictureBufferDescInitData.right_padding = 0;
    ThirtyTwoWidthPictureBufferDescInitData.top_padding = 0;
    ThirtyTwoWidthPictureBufferDescInitData.bot_padding = 0;
    ThirtyTwoWidthPictureBufferDescInitData.splitMode = EB_FALSE;

    // Candidate Ptr
    bufferPtr->candidate_ptr = (ModeDecisionCandidate_t*)EB_NULL;

    // Video Buffers
    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*)&(bufferPtr->prediction_ptr),
        (EbPtr)&pictureBufferDescInitData);

    // Video Buffers
    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*)&(bufferPtr->predictionPtrTemp),
        (EbPtr)&pictureBufferDescInitData);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*)&(bufferPtr->cflTempPredictionPtr),
        (EbPtr)&pictureBufferDescInitData);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }


    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*)&(bufferPtr->residual_ptr),
        (EbPtr)&doubleWidthPictureBufferDescInitData);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*)&(bufferPtr->residualQuantCoeffPtr),
        (EbPtr)&ThirtyTwoWidthPictureBufferDescInitData);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*)&(bufferPtr->reconCoeffPtr),
        (EbPtr)&ThirtyTwoWidthPictureBufferDescInitData);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*)&(bufferPtr->recon_ptr),
        (EbPtr)&pictureBufferDescInitData);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    //Distortion
    bufferPtr->residual_luma_sad = 0;

    bufferPtr->full_lambda_rate = 0;


    // Costs
    bufferPtr->fast_cost_ptr = fast_cost_ptr;
    bufferPtr->full_cost_ptr = full_cost_ptr;
    bufferPtr->full_cost_skip_ptr = full_cost_skip_ptr;
    bufferPtr->full_cost_merge_ptr = full_cost_merge_ptr;
    return EB_ErrorNone;
}


// Function Declarations
void RoundMv(
    ModeDecisionCandidate_t    *candidateArray,
    uint32_t                   canTotalCnt)
{

    candidateArray[canTotalCnt].motionVector_x_L0 = (candidateArray[canTotalCnt].motionVector_x_L0 + 2)&~0x03;
    candidateArray[canTotalCnt].motionVector_y_L0 = (candidateArray[canTotalCnt].motionVector_y_L0 + 2)&~0x03;

    candidateArray[canTotalCnt].motionVector_x_L1 = (candidateArray[canTotalCnt].motionVector_x_L1 + 2)&~0x03;
    candidateArray[canTotalCnt].motionVector_y_L1 = (candidateArray[canTotalCnt].motionVector_y_L1 + 2)&~0x03;

    return;

}

#if REMOVED_DUPLICATE_INTER
/***************************************
* return true if the MV candidate is already injected
***************************************/
EbBool is_already_injected_mv_l0(
    ModeDecisionContext_t *context_ptr,
    int16_t                mv_x,
    int16_t                mv_y) {

    for (int inter_candidate_index = 0; inter_candidate_index < context_ptr->injected_mv_count_l0; inter_candidate_index++) {
        if (context_ptr->injected_mv_x_l0_array[inter_candidate_index] == mv_x &&
            context_ptr->injected_mv_y_l0_array[inter_candidate_index] == mv_y) {
            return(EB_TRUE);
        }
    }

    return(EB_FALSE);
}

EbBool is_already_injected_mv_l1(
    ModeDecisionContext_t *context_ptr,
    int16_t                mv_x,
    int16_t                mv_y) {

    for (int inter_candidate_index = 0; inter_candidate_index < context_ptr->injected_mv_count_l1; inter_candidate_index++) {
        if (context_ptr->injected_mv_x_l1_array[inter_candidate_index] == mv_x &&
            context_ptr->injected_mv_y_l1_array[inter_candidate_index] == mv_y) {
            return(EB_TRUE);
        }
    }

    return(EB_FALSE);
}

EbBool is_already_injected_mv_bipred(
    ModeDecisionContext_t *context_ptr,
    int16_t                mv_x_l0,
    int16_t                mv_y_l0,
    int16_t                mv_x_l1,
    int16_t                mv_y_l1) {

    for (int inter_candidate_index = 0; inter_candidate_index < context_ptr->injected_mv_count_bipred; inter_candidate_index++) {
        if (context_ptr->injected_mv_x_bipred_l0_array[inter_candidate_index] == mv_x_l0 &&
            context_ptr->injected_mv_y_bipred_l0_array[inter_candidate_index] == mv_y_l0 &&
            context_ptr->injected_mv_x_bipred_l1_array[inter_candidate_index] == mv_x_l1 &&
            context_ptr->injected_mv_y_bipred_l1_array[inter_candidate_index] == mv_y_l1) {
            return(EB_TRUE);
        }
    }
    return(EB_FALSE);
}

#endif

EbErrorType SetMvpClipMVs(
    ModeDecisionCandidate_t  *candidate_ptr,
    uint32_t                    cu_origin_x,
    uint32_t                    cu_origin_y,
    uint32_t                    pu_index,
    uint32_t                    tbSize,
    PictureControlSet_t      *picture_control_set_ptr)
{
    EbErrorType  return_error = EB_ErrorNone;

    uint32_t        picture_width = ((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr)->luma_width;
    uint32_t        picture_height = ((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr)->luma_height;

    candidate_ptr->motion_vector_pred_idx[REF_LIST_0] = 0;
    candidate_ptr->motion_vector_pred_x[REF_LIST_0] = 0;
    candidate_ptr->motion_vector_pred_y[REF_LIST_0] = 0;
    candidate_ptr->motion_vector_pred_idx[REF_LIST_1] = 0;
    candidate_ptr->motion_vector_pred_x[REF_LIST_1] = 0;
    candidate_ptr->motion_vector_pred_y[REF_LIST_1] = 0;

    switch (candidate_ptr->prediction_direction[pu_index]) {
    case UNI_PRED_LIST_0:
        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motionVector_x_L0,
            &candidate_ptr->motionVector_y_L0,
            picture_width,
            picture_height,
            tbSize);

        break;

    case UNI_PRED_LIST_1:

        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motionVector_x_L1,
            &candidate_ptr->motionVector_y_L1,
            picture_width,
            picture_height,
            tbSize);

        break;

    case BI_PRED:

        // Choose the MVP in list0
        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motionVector_x_L0,
            &candidate_ptr->motionVector_y_L0,
            picture_width,
            picture_height,
            tbSize);

        // Choose the MVP in list1
        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motionVector_x_L1,
            &candidate_ptr->motionVector_y_L1,
            picture_width,
            picture_height,
            tbSize);
        break;

    default:
        break;
    }

    return return_error;
}


void LimitMvOverBound(
    int16_t *mvx,
    int16_t *mvy,
    ModeDecisionContext_t     *ctxtPtr,
    const SequenceControlSet_t      *sCSet)
{
    int32_t mvxF, mvyF;

    //L0
    mvxF = (*mvx) >> 2;
    mvyF = (*mvy) >> 2;

    if ((int32_t)ctxtPtr->cu_origin_x + mvxF + (int32_t)ctxtPtr->blk_geom->bwidth > (int32_t)sCSet->luma_width) {
        *mvx = (int16_t)(sCSet->luma_width - ctxtPtr->blk_geom->bwidth - ctxtPtr->cu_origin_x);
    }

    if ((int32_t)ctxtPtr->cu_origin_y + mvyF + (int32_t)ctxtPtr->blk_geom->bheight > (int32_t)sCSet->luma_height) {
        *mvy = (int16_t)(sCSet->luma_height - ctxtPtr->blk_geom->bheight - ctxtPtr->cu_origin_y);
    }

    if ((int32_t)ctxtPtr->cu_origin_x + mvxF < 0) {
        *mvx = -(int16_t)ctxtPtr->cu_origin_x;
    }

    if ((int32_t)ctxtPtr->cu_origin_y + mvyF < 0) {
        *mvy = -(int16_t)ctxtPtr->cu_origin_y;
    }


}


#if INTRA_INTER_FAST_LOOP
void sort_fast_loop_candidates(
    struct ModeDecisionContext_s   *context_ptr,
    uint32_t                        buffer_total_count,
    ModeDecisionCandidateBuffer_t **buffer_ptr_array,
    uint8_t                        *best_candidate_index_array,
    uint8_t                        *sorted_candidate_index_array,
    uint64_t                       *ref_fast_cost) {

    uint32_t fullReconCandidateCount = context_ptr->full_recon_search_count;

    //  move the scratch candidates (MAX_CU_COST) to the last spots (if any)
    uint32_t best_candidate_start_index = 0;
    uint32_t best_candidate_end_index = buffer_total_count - 1;
    for (uint8_t full_buffer_index = 0; full_buffer_index < buffer_total_count; full_buffer_index++) {
        if (*(buffer_ptr_array[full_buffer_index]->fast_cost_ptr) == MAX_CU_COST) {
            best_candidate_index_array[best_candidate_end_index--] = full_buffer_index;
        }
        else {
            best_candidate_index_array[best_candidate_start_index++] = full_buffer_index;
        }
    }

    // fl escape: inter then intra
    uint32_t i, j, index;
    for (i = 0; i < fullReconCandidateCount - 1; ++i) {
        for (j = i + 1; j < fullReconCandidateCount; ++j) {
            if ((buffer_ptr_array[best_candidate_index_array[i]]->candidate_ptr->type == INTRA_MODE) &&
                (buffer_ptr_array[best_candidate_index_array[j]]->candidate_ptr->type == INTER_MODE)) {
                index = best_candidate_index_array[i];
                best_candidate_index_array[i] = (uint8_t)best_candidate_index_array[j];
                best_candidate_index_array[j] = (uint8_t)index;
            }
        }
    }

    // tx search
    for (i = 0; i < fullReconCandidateCount; i++) {
        if (*(buffer_ptr_array[i]->fast_cost_ptr) < *ref_fast_cost) {
            *ref_fast_cost = *(buffer_ptr_array[i]->fast_cost_ptr);
        }
    }
    for (i = 0; i < MAX_NFL; ++i) {
        sorted_candidate_index_array[i] = best_candidate_index_array[i];
    }
    for (i = 0; i < fullReconCandidateCount - 1; ++i) {
        for (j = i + 1; j < fullReconCandidateCount; ++j) {
            if (*(buffer_ptr_array[j]->fast_cost_ptr) < *(buffer_ptr_array[i]->fast_cost_ptr)) {
                index = sorted_candidate_index_array[i];
                sorted_candidate_index_array[i] = (uint8_t)sorted_candidate_index_array[j];
                sorted_candidate_index_array[j] = (uint8_t)index;
            }
        }
    }
}
#else
/***************************************
* Pre-Mode Decision
*   Selects which fast cost modes to
*   do full reconstruction on.
***************************************/
EbErrorType PreModeDecision(
    CodingUnit_t                   *cu_ptr,
    uint32_t                          buffer_total_count,
    ModeDecisionCandidateBuffer_t **buffer_ptr_array,
    uint32_t                         *full_candidate_total_count_ptr,
    uint8_t                          *best_candidate_index_array,
#if USED_NFL_FEATURE_BASED
    uint8_t                          *sorted_candidate_index_array,
#endif
    uint8_t                          *disable_merge_index,
    uint64_t                         *ref_fast_cost,
    EbBool                           same_fast_full_candidate){

    UNUSED(cu_ptr);
    EbErrorType return_error = EB_ErrorNone;
    uint32_t fullCandidateIndex;
    uint32_t fullReconCandidateCount;
    uint32_t                          highestCostIndex;
    uint64_t                          highestCost;
    uint32_t                          candIndx = 0, i, j, index;

    *full_candidate_total_count_ptr = buffer_total_count;

    //Second,  We substract one, because with N buffers we can determine the best N-1 candidates.
    //Note/TODO: in the case number of fast candidate is less or equal to the number of buffers, N buffers would be enough
    if (same_fast_full_candidate)
        fullReconCandidateCount = MAX(1, (*full_candidate_total_count_ptr));
    else
        fullReconCandidateCount = MAX(1, (*full_candidate_total_count_ptr) - 1);

    //With N buffers, we get here with the best N-1, plus the last candidate. We need to exclude the worst, and keep the best N-1.
    highestCost = *(buffer_ptr_array[0]->fast_cost_ptr);
    highestCostIndex = 0;

    if (buffer_total_count > 1) {
        if (same_fast_full_candidate) {
            for (i = 0; i < buffer_total_count; i++) {
                best_candidate_index_array[candIndx++] = (uint8_t)i;
            }
        }
        else {
            for (i = 1; i < buffer_total_count; i++) {

                if (*(buffer_ptr_array[i]->fast_cost_ptr) >= highestCost) {
                    highestCost = *(buffer_ptr_array[i]->fast_cost_ptr);
                    highestCostIndex = i;
                }
            }
            for (i = 0; i < buffer_total_count; i++) {

                if (i != highestCostIndex) {
                    best_candidate_index_array[candIndx++] = (uint8_t)i;
                }
            }
        }

    }
    else
        best_candidate_index_array[0] = 0;
    for (i = 0; i < fullReconCandidateCount - 1; ++i) {
        for (j = i + 1; j < fullReconCandidateCount; ++j) {
            if ((buffer_ptr_array[best_candidate_index_array[i]]->candidate_ptr->type == INTRA_MODE) &&
                (buffer_ptr_array[best_candidate_index_array[j]]->candidate_ptr->type == INTER_MODE)) {
                index = best_candidate_index_array[i];
                best_candidate_index_array[i] = (uint8_t)best_candidate_index_array[j];
                best_candidate_index_array[j] = (uint8_t)index;
            }
        }
    }
    for (i = 0; i < fullReconCandidateCount; i++) {
        if (*(buffer_ptr_array[i]->fast_cost_ptr) < *ref_fast_cost) {
            *ref_fast_cost = *(buffer_ptr_array[i]->fast_cost_ptr);
        }
    }
#if USED_NFL_FEATURE_BASED
    for (i = 0; i < MAX_NFL; ++i) {
        sorted_candidate_index_array[i] = best_candidate_index_array[i];
    }

    for (i = 0; i < fullReconCandidateCount - 1; ++i) {
        for (j = i + 1; j < fullReconCandidateCount; ++j) {
            if (*(buffer_ptr_array[j]->fast_cost_ptr) < *(buffer_ptr_array[i]->fast_cost_ptr)) {
                index = sorted_candidate_index_array[i];
                sorted_candidate_index_array[i] = (uint8_t)sorted_candidate_index_array[j];
                sorted_candidate_index_array[j] = (uint8_t)index;
            }
        }
    }

#endif
    // Set (*full_candidate_total_count_ptr) to fullReconCandidateCount
    (*full_candidate_total_count_ptr) = fullReconCandidateCount;

    for (i = 0; i < fullReconCandidateCount; ++i) {
        fullCandidateIndex = best_candidate_index_array[i];

        // Set disable_merge_index
        *disable_merge_index = buffer_ptr_array[fullCandidateIndex]->candidate_ptr->type == INTER_MODE ? 1 : *disable_merge_index;
    }

    return return_error;
}
#endif
#if IMPROVED_BIPRED_INJECTION || IMPROVED_UNIPRED_INJECTION

#define BIPRED_3x3_REFINMENT_POSITIONS 8

int8_t ALLOW_REFINEMENT_FLAG[BIPRED_3x3_REFINMENT_POSITIONS] = {  1, 0, 1, 0, 1,  0,  1, 0 };
int8_t BIPRED_3x3_X_POS[BIPRED_3x3_REFINMENT_POSITIONS] = { -1, -1, 0, 1, 1, 1, 0, -1 };
int8_t BIPRED_3x3_Y_POS[BIPRED_3x3_REFINMENT_POSITIONS] = { 0, 1, 1, 1, 0, -1, -1, -1 };
#endif
#if IMPROVED_UNIPRED_INJECTION
void Unipred3x3CandidatesInjection(
    PictureControlSet_t            *picture_control_set_ptr,
    ModeDecisionContext_t          *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                        me_sb_addr,
    SsMeContext_t                  *inloop_me_context,
    EbBool                          use_close_loop_me,
    uint32_t                        close_loop_me_index,
    const uint32_t                  me2Nx2NTableOffset,
    uint32_t                       *candidateTotalCnt){

    UNUSED(sb_ptr);
    uint32_t                   bipredIndex;
    uint32_t                   canTotalCnt = (*candidateTotalCnt);
    MeCuResults_t * mePuResult = &picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr][me2Nx2NTableOffset];
    ModeDecisionCandidate_t    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    IntMv  bestPredmv[2] = { {0}, {0} };

    // (8 Best_L0 neighbors)
    for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
    {
        /**************
        NEWMV L0
        ************* */
        if (context_ptr->unipred3x3_injection >= 2)
            if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                continue;
#if REMOVED_DUPLICATE_INTER
        int16_t to_inject_mv_x = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (mePuResult->xMvL0 + BIPRED_3x3_X_POS[bipredIndex]) << 1;
        int16_t to_inject_mv_y = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (mePuResult->yMvL0 + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
        if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
            candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
            candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
            candidateArray[canTotalCnt].use_intrabc = 0;
#endif
            candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].merge_index = 0;
            candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif
            candidateArray[canTotalCnt].inter_mode = NEWMV;
            candidateArray[canTotalCnt].pred_mode = NEWMV;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

            candidateArray[canTotalCnt].is_compound = 0;
            candidateArray[canTotalCnt].is_new_mv = 1;
            candidateArray[canTotalCnt].is_zero_mv = 0;

            candidateArray[canTotalCnt].drl_index = 0;

            // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER
            candidateArray[canTotalCnt].motionVector_x_L0 = to_inject_mv_x;
            candidateArray[canTotalCnt].motionVector_y_L0 = to_inject_mv_y;
#else
            candidateArray[canTotalCnt].motionVector_x_L0 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (mePuResult->xMvL0 + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            candidateArray[canTotalCnt].motionVector_y_L0 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (mePuResult->yMvL0 + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
#endif

            // will be needed later by the rate estimation
            candidateArray[canTotalCnt].ref_mv_index = 0;
            candidateArray[canTotalCnt].pred_mv_weight = 0;
            candidateArray[canTotalCnt].ref_frame_type = LAST_FRAME;


            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

            ChooseBestAv1MvPred(
                context_ptr,
                candidateArray[canTotalCnt].md_rate_estimation_ptr,
                context_ptr->cu_ptr,
                candidateArray[canTotalCnt].ref_frame_type,
                candidateArray[canTotalCnt].is_compound,
                candidateArray[canTotalCnt].pred_mode,
                candidateArray[canTotalCnt].motionVector_x_L0,
                candidateArray[canTotalCnt].motionVector_y_L0,
                0, 0,
                &candidateArray[canTotalCnt].drl_index,
                bestPredmv);

            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;

            ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER
            context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
            context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
            ++context_ptr->injected_mv_count_l0;
        }
#endif
    }
    // (8 Best_L1 neighbors)
    for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
    {
        if (isCompoundEnabled) {
            /**************
            NEWMV L1
            ************* */
        if (context_ptr->unipred3x3_injection >= 2)
            if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                continue;
#if REMOVED_DUPLICATE_INTER_L1
            int16_t to_inject_mv_x = use_close_loop_me ? (inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (mePuResult->xMvL1 + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            int16_t to_inject_mv_y = use_close_loop_me ? (inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (mePuResult->yMvL1 + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
            if (context_ptr->injected_mv_count_l1 == 0 || is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
                candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
                candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
                candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
                candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
                candidateArray[canTotalCnt].use_intrabc = 0;
#endif
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canTotalCnt].merge_index = 0;
                candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)1;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif
                candidateArray[canTotalCnt].inter_mode = NEWMV;
                candidateArray[canTotalCnt].pred_mode = NEWMV;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

                candidateArray[canTotalCnt].is_compound = 0;
                candidateArray[canTotalCnt].is_new_mv = 1;
                candidateArray[canTotalCnt].is_zero_mv = 0;

                candidateArray[canTotalCnt].drl_index = 0;

                // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER_L1
                candidateArray[canTotalCnt].motionVector_x_L1 = to_inject_mv_x;
                candidateArray[canTotalCnt].motionVector_y_L1 = to_inject_mv_y;
#else
                candidateArray[canTotalCnt].motionVector_x_L1 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (mePuResult->xMvL1 + BIPRED_3x3_X_POS[bipredIndex]) << 1;
                candidateArray[canTotalCnt].motionVector_y_L1 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (mePuResult->yMvL1 + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
#endif
                // will be needed later by the rate estimation
                candidateArray[canTotalCnt].ref_mv_index = 0;
                candidateArray[canTotalCnt].pred_mv_weight = 0;
                candidateArray[canTotalCnt].ref_frame_type = BWDREF_FRAME;


                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

                ChooseBestAv1MvPred(
                    context_ptr,
                    candidateArray[canTotalCnt].md_rate_estimation_ptr,
                    context_ptr->cu_ptr,
                    candidateArray[canTotalCnt].ref_frame_type,
                    candidateArray[canTotalCnt].is_compound,
                    candidateArray[canTotalCnt].pred_mode,
                    candidateArray[canTotalCnt].motionVector_x_L1,
                    candidateArray[canTotalCnt].motionVector_y_L1,
                    0, 0,
                    &candidateArray[canTotalCnt].drl_index,
                    bestPredmv);

                candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[0].as_mv.col;
                candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[0].as_mv.row;

                ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER_L1
                context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
                context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
                ++context_ptr->injected_mv_count_l1;
            }
#endif
        }
    }

    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;

    return;
}
#endif
#if IMPROVED_BIPRED_INJECTION
void Bipred3x3CandidatesInjection(
    PictureControlSet_t            *picture_control_set_ptr,
    ModeDecisionContext_t          *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                        me_sb_addr,
    SsMeContext_t                  *inloop_me_context,
    EbBool                          use_close_loop_me,
    uint32_t                        close_loop_me_index,
    const uint32_t                  me2Nx2NTableOffset,
    uint32_t                       *candidateTotalCnt){

    UNUSED(sb_ptr);
    uint32_t                   bipredIndex;
    uint32_t                   canTotalCnt = (*candidateTotalCnt);
    MeCuResults_t * mePuResult = &picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr][me2Nx2NTableOffset];
    ModeDecisionCandidate_t    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    IntMv  bestPredmv[2] = { {0}, {0} };

    if (isCompoundEnabled) {
        /**************
       NEW_NEWMV
       ************* */
       // (Best_L0, 8 Best_L1 neighbors)
        for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
        {
        if (context_ptr->bipred3x3_injection >= 2)
            if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                continue;
#if REMOVED_DUPLICATE_INTER_BIPRED
            int16_t to_inject_mv_x_l0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : mePuResult->xMvL0 << 1;
            int16_t to_inject_mv_y_l0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : mePuResult->yMvL0 << 1;
            int16_t to_inject_mv_x_l1 = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1) : (mePuResult->xMvL1 + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            int16_t to_inject_mv_y_l1 = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1) : (mePuResult->yMvL1 + BIPRED_3x3_Y_POS[bipredIndex]) << 1;

            if (context_ptr->injected_mv_count_bipred == 0 || is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1) == EB_FALSE) {
#endif
            candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
            candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
            candidateArray[canTotalCnt].use_intrabc = 0;
#endif
            candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].merge_index = 0;
            candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
            candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif


            candidateArray[canTotalCnt].is_new_mv = 1;
            candidateArray[canTotalCnt].is_zero_mv = 0;

            candidateArray[canTotalCnt].drl_index = 0;

            // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER_BIPRED
            candidateArray[canTotalCnt].motionVector_x_L0 = to_inject_mv_x_l0;
            candidateArray[canTotalCnt].motionVector_y_L0 = to_inject_mv_y_l0;
            candidateArray[canTotalCnt].motionVector_x_L1 = to_inject_mv_x_l1;
            candidateArray[canTotalCnt].motionVector_y_L1 = to_inject_mv_y_l1;
#else
            candidateArray[canTotalCnt].motionVector_x_L0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : mePuResult->xMvL0 << 1;
            candidateArray[canTotalCnt].motionVector_y_L0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : mePuResult->yMvL0 << 1;

            candidateArray[canTotalCnt].motionVector_x_L1 = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1) : (mePuResult->xMvL1 + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            candidateArray[canTotalCnt].motionVector_y_L1 = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1) : (mePuResult->yMvL1 + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
#endif
            // will be needed later by the rate estimation
            candidateArray[canTotalCnt].ref_mv_index = 0;
            candidateArray[canTotalCnt].pred_mv_weight = 0;

            candidateArray[canTotalCnt].inter_mode = NEW_NEWMV;
            candidateArray[canTotalCnt].pred_mode = NEW_NEWMV;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            candidateArray[canTotalCnt].is_compound = 1;
            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
            candidateArray[canTotalCnt].ref_frame_type = LAST_BWD_FRAME;

            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

            ChooseBestAv1MvPred(
                context_ptr,
                candidateArray[canTotalCnt].md_rate_estimation_ptr,
                context_ptr->cu_ptr,
                candidateArray[canTotalCnt].ref_frame_type,
                candidateArray[canTotalCnt].is_compound,
                candidateArray[canTotalCnt].pred_mode,
                candidateArray[canTotalCnt].motionVector_x_L0,
                candidateArray[canTotalCnt].motionVector_y_L0,
                candidateArray[canTotalCnt].motionVector_x_L1,
                candidateArray[canTotalCnt].motionVector_y_L1,
                &candidateArray[canTotalCnt].drl_index,
                bestPredmv);

            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;

            ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER_BIPRED
            context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
            context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
            context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
            context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
            ++context_ptr->injected_mv_count_bipred;
            }
#endif
        }

        // (8 Best_L0 neighbors, Best_L1) :
        for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
        {
        if (context_ptr->bipred3x3_injection >= 2)
            if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                continue;
#if REMOVED_DUPLICATE_INTER_BIPRED
            int16_t to_inject_mv_x_l0 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (mePuResult->xMvL0 + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            int16_t to_inject_mv_y_l0 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (mePuResult->yMvL0 + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
            int16_t to_inject_mv_x_l1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : mePuResult->xMvL1 << 1;
            int16_t to_inject_mv_y_l1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : mePuResult->yMvL1 << 1;

            if (context_ptr->injected_mv_count_bipred == 0 || is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1) == EB_FALSE) {
#endif
            candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
            candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
            candidateArray[canTotalCnt].use_intrabc = 0;
#endif
            candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].merge_index = 0;
            candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
            candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif

            candidateArray[canTotalCnt].is_new_mv = 1;
            candidateArray[canTotalCnt].is_zero_mv = 0;

            candidateArray[canTotalCnt].drl_index = 0;

            // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER_BIPRED
            candidateArray[canTotalCnt].motionVector_x_L0 = to_inject_mv_x_l0;
            candidateArray[canTotalCnt].motionVector_y_L0 = to_inject_mv_y_l0;
            candidateArray[canTotalCnt].motionVector_x_L1 = to_inject_mv_x_l1;
            candidateArray[canTotalCnt].motionVector_y_L1 = to_inject_mv_y_l1;
#else
            candidateArray[canTotalCnt].motionVector_x_L0 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (mePuResult->xMvL0 + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            candidateArray[canTotalCnt].motionVector_y_L0 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (mePuResult->yMvL0 + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
            candidateArray[canTotalCnt].motionVector_x_L1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : mePuResult->xMvL1 << 1;
            candidateArray[canTotalCnt].motionVector_y_L1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : mePuResult->yMvL1 << 1;
#endif
            // will be needed later by the rate estimation
            candidateArray[canTotalCnt].ref_mv_index = 0;
            candidateArray[canTotalCnt].pred_mv_weight = 0;

            candidateArray[canTotalCnt].inter_mode = NEW_NEWMV;
            candidateArray[canTotalCnt].pred_mode = NEW_NEWMV;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            candidateArray[canTotalCnt].is_compound = 1;
            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
            candidateArray[canTotalCnt].ref_frame_type = LAST_BWD_FRAME;

            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

            ChooseBestAv1MvPred(
                context_ptr,
                candidateArray[canTotalCnt].md_rate_estimation_ptr,
                context_ptr->cu_ptr,
                candidateArray[canTotalCnt].ref_frame_type,
                candidateArray[canTotalCnt].is_compound,
                candidateArray[canTotalCnt].pred_mode,
                candidateArray[canTotalCnt].motionVector_x_L0,
                candidateArray[canTotalCnt].motionVector_y_L0,
                candidateArray[canTotalCnt].motionVector_x_L1,
                candidateArray[canTotalCnt].motionVector_y_L1,
                &candidateArray[canTotalCnt].drl_index,
                bestPredmv);

            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;

            ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER_BIPRED
            context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
            context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
            context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
            context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
            ++context_ptr->injected_mv_count_bipred;
        }
#endif
        }

    }

    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;

    return;
}
#endif

uint8_t GetMaxDrlIndex(uint8_t  refmvCnt, PredictionMode   mode)
{
    uint8_t maxDrl = 0;

    if (mode == NEWMV || mode == NEW_NEWMV) {

        if (refmvCnt < 2)
            maxDrl = 1;
        else if (refmvCnt == 2)
            maxDrl = 2;
        else
            maxDrl = 3;
    }

    if (mode == NEARMV || mode == NEAR_NEARMV || mode == NEAR_NEWMV || mode == NEW_NEARMV) {

        if (refmvCnt < 3)
            maxDrl = 1;
        else if (refmvCnt == 3)
            maxDrl = 2;
        else
            maxDrl = 3;
    }

    return maxDrl;
}
/*********************************************************************
**********************************************************************
        Upto 12 inter Candidated injected
        Min 6 inter Candidated injected
UniPred L0 : NEARST         + upto 3x NEAR
UniPred L1 : NEARST         + upto 3x NEAR
BIPred     : NEARST_NEARST  + upto 3x NEAR_NEAR
**********************************************************************
**********************************************************************/
void InjectAv1MvpCandidates(
    struct ModeDecisionContext_s     *context_ptr,
    CodingUnit_t                     *cu_ptr,
    MvReferenceFrame               *refFrames,
    PictureControlSet_t              *picture_control_set_ptr,
    uint32_t                            lcuAddr,
#if !M8_SKIP_BLK    
    uint32_t                            leaf_index,
#endif    
    EbBool                           allow_bipred,
    uint32_t                           *candTotCnt)
{
#if !M8_SKIP_BLK  
    (void)leaf_index;
#endif
    (void)lcuAddr;
    (void)refFrames;
    uint32_t                   canIdx = *candTotCnt;
    ModeDecisionCandidate_t    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
#if  !BASE_LAYER_REF
    isCompoundEnabled = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : isCompoundEnabled;
#endif
    MacroBlockD  *xd = cu_ptr->av1xd;
    uint8_t drli, maxDrlIndex;
    IntMv    nearestmv[2], nearmv[2], ref_mv[2];

    //NEAREST_L0
#if REMOVED_DUPLICATE_INTER
    int16_t to_inject_mv_x = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.col;
    int16_t to_inject_mv_y = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.row;
    if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
    candidateArray[canIdx].type = INTER_MODE;
    candidateArray[canIdx].inter_mode = NEARESTMV;
    candidateArray[canIdx].pred_mode = NEARESTMV;
    candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
    candidateArray[canIdx].is_compound = 0;
#if TWO_FAST_LOOP 
    candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
    candidateArray[canIdx].distortion_ready = 0;
#endif
#else
    candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
    candidateArray[canIdx].use_intrabc = 0;
#endif
    candidateArray[canIdx].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
    candidateArray[canIdx].merge_index = 0;
    candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
    candidateArray[canIdx].prediction_direction[0] = UNI_PRED_LIST_0;
#if !INTRA_INTER_FAST_LOOP
    candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
    candidateArray[canIdx].is_new_mv = 0;
    candidateArray[canIdx].is_zero_mv = 0;
#if REMOVED_DUPLICATE_INTER
    candidateArray[canIdx].motionVector_x_L0 = to_inject_mv_x;
    candidateArray[canIdx].motionVector_y_L0 = to_inject_mv_y;
#else
    candidateArray[canIdx].motionVector_x_L0 = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.col;
    candidateArray[canIdx].motionVector_y_L0 = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.row;
#endif
    candidateArray[canIdx].drl_index = 0;
    candidateArray[canIdx].ref_mv_index = 0;
    candidateArray[canIdx].pred_mv_weight = 0;
    candidateArray[canIdx].ref_frame_type = LAST_FRAME;
    candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
    candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;
    ++canIdx;
#if REMOVED_DUPLICATE_INTER
    context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
    context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
    ++context_ptr->injected_mv_count_l0;
    }
#endif
    //NEAR_L0
    maxDrlIndex = GetMaxDrlIndex(xd->ref_mv_count[LAST_FRAME], NEARMV);
    //maxDrlIndex = 1;
    for (drli = 0; drli < maxDrlIndex; drli++) {

        get_av1_mv_pred_drl(
            context_ptr,
            cu_ptr,
            LAST_FRAME,
            0,
            NEARMV,
            drli,
            nearestmv,
            nearmv,
            ref_mv);

#if REMOVED_DUPLICATE_INTER
        int16_t to_inject_mv_x = nearmv[0].as_mv.col;
        int16_t to_inject_mv_y = nearmv[0].as_mv.row;
        if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
        candidateArray[canIdx].type = INTER_MODE;
        candidateArray[canIdx].inter_mode = NEARMV;
        candidateArray[canIdx].pred_mode = NEARMV;
        candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
        candidateArray[canIdx].is_compound = 0;
#if TWO_FAST_LOOP 
        candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
        candidateArray[canIdx].distortion_ready = 0;
#endif
#else
        candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
        candidateArray[canIdx].use_intrabc = 0;
#endif
        candidateArray[canIdx].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canIdx].merge_index = 0;
        candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
        candidateArray[canIdx].prediction_direction[0] = UNI_PRED_LIST_0;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
        candidateArray[canIdx].is_new_mv = 0;
        candidateArray[canIdx].is_zero_mv = 0;
#if REMOVED_DUPLICATE_INTER
        candidateArray[canIdx].motionVector_x_L0 = to_inject_mv_x;
        candidateArray[canIdx].motionVector_y_L0 = to_inject_mv_y;
#else
        candidateArray[canIdx].motionVector_x_L0 = nearmv[0].as_mv.col;
        candidateArray[canIdx].motionVector_y_L0 = nearmv[0].as_mv.row;
#endif
        candidateArray[canIdx].drl_index = drli;
        candidateArray[canIdx].ref_mv_index = 0;
        candidateArray[canIdx].pred_mv_weight = 0;
        candidateArray[canIdx].ref_frame_type = LAST_FRAME;
        candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
        candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;
        ++canIdx;
#if REMOVED_DUPLICATE_INTER
        context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
        context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
        ++context_ptr->injected_mv_count_l0;
        }
#endif
    }


    if (isCompoundEnabled) {
#if REMOVED_DUPLICATE_INTER_L1
        int16_t to_inject_mv_x = context_ptr->cu_ptr->ref_mvs[BWDREF_FRAME][0].as_mv.col;
        int16_t to_inject_mv_y = context_ptr->cu_ptr->ref_mvs[BWDREF_FRAME][0].as_mv.row;
        if (context_ptr->injected_mv_count_l1 == 0 || is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
        //NEAREST_L1
        candidateArray[canIdx].type = INTER_MODE;
        candidateArray[canIdx].inter_mode = NEARESTMV;
        candidateArray[canIdx].pred_mode = NEARESTMV;
        candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
        candidateArray[canIdx].is_compound = 0;
#if TWO_FAST_LOOP 
        candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
        candidateArray[canIdx].distortion_ready = 0;
#endif
#else
        candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
        candidateArray[canIdx].use_intrabc = 0;
#endif
        candidateArray[canIdx].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canIdx].merge_index = 0;
        candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
        candidateArray[canIdx].prediction_direction[0] = UNI_PRED_LIST_1;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
        candidateArray[canIdx].is_new_mv = 0;
        candidateArray[canIdx].is_zero_mv = 0;
#if REMOVED_DUPLICATE_INTER_L1
        candidateArray[canIdx].motionVector_x_L1 = to_inject_mv_x;
        candidateArray[canIdx].motionVector_y_L1 = to_inject_mv_y;
#else
        candidateArray[canIdx].motionVector_x_L1 = context_ptr->cu_ptr->ref_mvs[BWDREF_FRAME][0].as_mv.col;
        candidateArray[canIdx].motionVector_y_L1 = context_ptr->cu_ptr->ref_mvs[BWDREF_FRAME][0].as_mv.row;
#endif
        candidateArray[canIdx].drl_index = 0;
        candidateArray[canIdx].ref_mv_index = 0;
        candidateArray[canIdx].pred_mv_weight = 0;
        candidateArray[canIdx].ref_frame_type = BWDREF_FRAME;
        candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
        candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;
        ++canIdx;
#if REMOVED_DUPLICATE_INTER_L1
        context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
        context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
        ++context_ptr->injected_mv_count_l1;
    }
#endif
        //NEAR_L1
        maxDrlIndex = GetMaxDrlIndex(xd->ref_mv_count[BWDREF_FRAME], NEARMV);

        for (drli = 0; drli < maxDrlIndex; drli++) {

            get_av1_mv_pred_drl(
                context_ptr,
                cu_ptr,
                BWDREF_FRAME,
                0,
                NEARMV,
                drli,
                nearestmv,
                nearmv,
                ref_mv);
#if REMOVED_DUPLICATE_INTER_L1
            int16_t to_inject_mv_x = nearmv[0].as_mv.col;
            int16_t to_inject_mv_y = nearmv[0].as_mv.row;
            if (context_ptr->injected_mv_count_l1 == 0 || is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
            candidateArray[canIdx].type = INTER_MODE;
            candidateArray[canIdx].inter_mode = NEARMV;
            candidateArray[canIdx].pred_mode = NEARMV;
            candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
            candidateArray[canIdx].is_compound = 0;
#if TWO_FAST_LOOP 
            candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
            candidateArray[canIdx].distortion_ready = 0;
#endif
#else
            candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
            candidateArray[canIdx].use_intrabc = 0;
#endif
            candidateArray[canIdx].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canIdx].merge_index = 0;
            candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
            candidateArray[canIdx].prediction_direction[0] = UNI_PRED_LIST_1;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
            candidateArray[canIdx].is_new_mv = 0;
            candidateArray[canIdx].is_zero_mv = 0;
#if REMOVED_DUPLICATE_INTER_L1
            candidateArray[canIdx].motionVector_x_L1 = to_inject_mv_x;
            candidateArray[canIdx].motionVector_y_L1 = to_inject_mv_y;
#else
            candidateArray[canIdx].motionVector_x_L1 = nearmv[0].as_mv.col;
            candidateArray[canIdx].motionVector_y_L1 = nearmv[0].as_mv.row;
#endif
            candidateArray[canIdx].drl_index = drli;
            candidateArray[canIdx].ref_mv_index = 0;
            candidateArray[canIdx].pred_mv_weight = 0;
            candidateArray[canIdx].ref_frame_type = BWDREF_FRAME;
            candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
            candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;
            ++canIdx;
#if REMOVED_DUPLICATE_INTER_L1
            context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
            context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
            ++context_ptr->injected_mv_count_l1;
            }
#endif
        }

#if BASE_LAYER_REF
        if (allow_bipred)
#endif
        {
            //SKIP (NEAREST_NEAREST with LAST_BWD_FRAME)
#if REMOVED_DUPLICATE_INTER_BIPRED
            int16_t to_inject_mv_x_l0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[LAST_BWD_FRAME][0].this_mv.as_mv.col;
            int16_t to_inject_mv_y_l0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[LAST_BWD_FRAME][0].this_mv.as_mv.row;
            int16_t to_inject_mv_x_l1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[LAST_BWD_FRAME][0].comp_mv.as_mv.col;
            int16_t to_inject_mv_y_l1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[LAST_BWD_FRAME][0].comp_mv.as_mv.row;
            if (context_ptr->injected_mv_count_bipred == 0 || is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1) == EB_FALSE) {
#endif
                candidateArray[canIdx].type = INTER_MODE;
                candidateArray[canIdx].inter_mode = NEAREST_NEARESTMV;
                candidateArray[canIdx].pred_mode = NEAREST_NEARESTMV;
                candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canIdx].is_compound = 1;
#if TWO_FAST_LOOP 
                candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
                candidateArray[canIdx].distortion_ready = 0;
#endif
#else
                candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
                candidateArray[canIdx].use_intrabc = 0;
#endif
                candidateArray[canIdx].merge_flag = EB_TRUE;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canIdx].merge_index = 0;
                candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
                candidateArray[canIdx].prediction_direction[0] = BI_PRED;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
                candidateArray[canIdx].is_new_mv = 0;
                candidateArray[canIdx].is_zero_mv = 0;
#if REMOVED_DUPLICATE_INTER_BIPRED
                candidateArray[canIdx].motionVector_x_L0 = to_inject_mv_x_l0;
                candidateArray[canIdx].motionVector_y_L0 = to_inject_mv_y_l0;
                candidateArray[canIdx].motionVector_x_L1 = to_inject_mv_x_l1;
                candidateArray[canIdx].motionVector_y_L1 = to_inject_mv_y_l1;
#else
                candidateArray[canIdx].motionVector_x_L0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[LAST_BWD_FRAME][0].this_mv.as_mv.col;
                candidateArray[canIdx].motionVector_y_L0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[LAST_BWD_FRAME][0].this_mv.as_mv.row;
                candidateArray[canIdx].motionVector_x_L1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[LAST_BWD_FRAME][0].comp_mv.as_mv.col;
                candidateArray[canIdx].motionVector_y_L1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[LAST_BWD_FRAME][0].comp_mv.as_mv.row;
#endif
                candidateArray[canIdx].drl_index = 0;
                candidateArray[canIdx].ref_mv_index = 0;
                candidateArray[canIdx].pred_mv_weight = 0;
                candidateArray[canIdx].ref_frame_type = LAST_BWD_FRAME;
                candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;
                ++canIdx;
#if REMOVED_DUPLICATE_INTER_BIPRED
                context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                ++context_ptr->injected_mv_count_bipred;
            }
#endif
        }
        //NEAR_NEAR
        if (allow_bipred) {

            maxDrlIndex = GetMaxDrlIndex(xd->ref_mv_count[LAST_BWD_FRAME], NEAR_NEARMV);
            //maxDrlIndex = 1;
            for (drli = 0; drli < maxDrlIndex; drli++) {

                get_av1_mv_pred_drl(
                    context_ptr,
                    cu_ptr,
                    LAST_BWD_FRAME,
                    1,
                    NEAR_NEARMV,
                    drli,
                    nearestmv,
                    nearmv,
                    ref_mv);
#if REMOVED_DUPLICATE_INTER_BIPRED
                int16_t to_inject_mv_x_l0 = nearmv[0].as_mv.col;
                int16_t to_inject_mv_y_l0 = nearmv[0].as_mv.row;
                int16_t to_inject_mv_x_l1 = nearmv[1].as_mv.col;
                int16_t to_inject_mv_y_l1 = nearmv[1].as_mv.row;
                if (context_ptr->injected_mv_count_bipred == 0 || is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1) == EB_FALSE) {
#endif
                candidateArray[canIdx].type = INTER_MODE;
                candidateArray[canIdx].inter_mode = NEAR_NEARMV;
                candidateArray[canIdx].pred_mode = NEAR_NEARMV;
                candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canIdx].is_compound = 1;
#if TWO_FAST_LOOP 
                candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
                candidateArray[canIdx].distortion_ready = 0;
#endif
#else
                candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
                candidateArray[canIdx].use_intrabc = 0;
#endif
                candidateArray[canIdx].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canIdx].merge_index = 0;
                candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
                candidateArray[canIdx].prediction_direction[0] = BI_PRED;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
                candidateArray[canIdx].is_new_mv = 0;
                candidateArray[canIdx].is_zero_mv = 0;
#if REMOVED_DUPLICATE_INTER_BIPRED
                candidateArray[canIdx].motionVector_x_L0 = to_inject_mv_x_l0;
                candidateArray[canIdx].motionVector_y_L0 = to_inject_mv_y_l0;
                candidateArray[canIdx].motionVector_x_L1 = to_inject_mv_x_l1;
                candidateArray[canIdx].motionVector_y_L1 = to_inject_mv_y_l1;
#else
                candidateArray[canIdx].motionVector_x_L0 = nearmv[0].as_mv.col;
                candidateArray[canIdx].motionVector_y_L0 = nearmv[0].as_mv.row;
                candidateArray[canIdx].motionVector_x_L1 = nearmv[1].as_mv.col;
                candidateArray[canIdx].motionVector_y_L1 = nearmv[1].as_mv.row;
#endif
                candidateArray[canIdx].drl_index = drli;
                candidateArray[canIdx].ref_mv_index = 0;
                candidateArray[canIdx].pred_mv_weight = 0;
                candidateArray[canIdx].ref_frame_type = LAST_BWD_FRAME;
                candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;
                ++canIdx;
#if REMOVED_DUPLICATE_INTER_BIPRED
                context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                ++context_ptr->injected_mv_count_bipred;
                }
#endif
            }
        }
    }
    //update tot Candidate count
    *candTotCnt = canIdx;
}

void inject_warped_motion_candidates(
    PictureControlSet_t              *picture_control_set_ptr,
    struct ModeDecisionContext_s     *context_ptr,
    CodingUnit_t                     *cu_ptr,
    uint32_t                         *candTotCnt,
    SsMeContext_t                    *ss_mecontext,
    MeCuResults_t                    *mePuResult,
    EbBool                            use_close_loop_me,
    uint32_t                          close_loop_me_index)
{
    uint32_t canIdx = *candTotCnt;
    ModeDecisionCandidate_t *candidateArray = context_ptr->fast_candidate_array;
    MacroBlockD  *xd = cu_ptr->av1xd;
    uint8_t drli, maxDrlIndex;
    IntMv nearestmv[2], nearmv[2], ref_mv[2];

    //NEAREST_L0
#if REMOVED_DUPLICATE_INTER
    int16_t to_inject_mv_x = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.col;
    int16_t to_inject_mv_y = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.row;
    if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
        candidateArray[canIdx].type = INTER_MODE;
        candidateArray[canIdx].inter_mode = NEARESTMV;
        candidateArray[canIdx].pred_mode = NEARESTMV;
        candidateArray[canIdx].motion_mode = WARPED_CAUSAL;
        candidateArray[canIdx].wm_params.wmtype = AFFINE;
        candidateArray[canIdx].is_compound = 0;
#if TWO_FAST_LOOP 
        candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
        candidateArray[canIdx].distortion_ready = 0;
#endif
#else
        candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
        candidateArray[canIdx].use_intrabc = 0;
#endif
        candidateArray[canIdx].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canIdx].merge_index = 0;
        candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
        candidateArray[canIdx].prediction_direction[0] = UNI_PRED_LIST_0;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
        candidateArray[canIdx].is_new_mv = 0;
        candidateArray[canIdx].is_zero_mv = 0;
#if REMOVED_DUPLICATE_INTER
        candidateArray[canIdx].motionVector_x_L0 = to_inject_mv_x;
        candidateArray[canIdx].motionVector_y_L0 = to_inject_mv_y;
#else
        candidateArray[canIdx].motionVector_x_L0 = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.col;
        candidateArray[canIdx].motionVector_y_L0 = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.row;
#endif
        candidateArray[canIdx].drl_index = 0;
        candidateArray[canIdx].ref_mv_index = 0;
        candidateArray[canIdx].pred_mv_weight = 0;
        candidateArray[canIdx].ref_frame_type = LAST_FRAME;
        candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
        candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;

        MvUnit_t mv_unit;
        mv_unit.mv[0].x = candidateArray[canIdx].motionVector_x_L0;
        mv_unit.mv[0].y = candidateArray[canIdx].motionVector_y_L0;
        candidateArray[canIdx].local_warp_valid = warped_motion_parameters(
            picture_control_set_ptr,
            context_ptr->cu_ptr,
            &mv_unit,
            context_ptr->blk_geom,
            context_ptr->cu_origin_x,
            context_ptr->cu_origin_y,
            candidateArray[canIdx].ref_frame_type,
            &candidateArray[canIdx].wm_params,
            &candidateArray[canIdx].num_proj_ref);

        if (candidateArray[canIdx].local_warp_valid)
            ++canIdx;
#if REMOVED_DUPLICATE_INTER
        context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
        context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
        ++context_ptr->injected_mv_count_l0;
    }
#endif

    //NEAR_L0
    maxDrlIndex = GetMaxDrlIndex(xd->ref_mv_count[LAST_FRAME], NEARMV);
    for (drli = 0; drli < maxDrlIndex; drli++) {
        get_av1_mv_pred_drl(
            context_ptr,
            cu_ptr,
            LAST_FRAME,
            0,
            NEARMV,
            drli,
            nearestmv,
            nearmv,
            ref_mv);
#if REMOVED_DUPLICATE_INTER
        int16_t to_inject_mv_x = nearmv[0].as_mv.col;
        int16_t to_inject_mv_y = nearmv[0].as_mv.row;
        if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
            candidateArray[canIdx].type = INTER_MODE;
            candidateArray[canIdx].inter_mode = NEARMV;
            candidateArray[canIdx].pred_mode = NEARMV;
            candidateArray[canIdx].motion_mode = WARPED_CAUSAL;
            candidateArray[canIdx].wm_params.wmtype = AFFINE;
            candidateArray[canIdx].is_compound = 0;
#if TWO_FAST_LOOP 
            candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
            candidateArray[canIdx].distortion_ready = 0;
#endif
#else
            candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
            candidateArray[canIdx].use_intrabc = 0;
#endif
            candidateArray[canIdx].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canIdx].merge_index = 0;
            candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
            candidateArray[canIdx].prediction_direction[0] = UNI_PRED_LIST_0;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
            candidateArray[canIdx].is_new_mv = 0;
            candidateArray[canIdx].is_zero_mv = 0;
#if REMOVED_DUPLICATE_INTER
            candidateArray[canIdx].motionVector_x_L0 = to_inject_mv_x;
            candidateArray[canIdx].motionVector_y_L0 = to_inject_mv_y;
#else
            candidateArray[canIdx].motionVector_x_L0 = nearmv[0].as_mv.col;
            candidateArray[canIdx].motionVector_y_L0 = nearmv[0].as_mv.row;
#endif
            candidateArray[canIdx].drl_index = drli;
            candidateArray[canIdx].ref_mv_index = 0;
            candidateArray[canIdx].pred_mv_weight = 0;
            candidateArray[canIdx].ref_frame_type = LAST_FRAME;
            candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
            candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;

            MvUnit_t mv_unit;
            mv_unit.mv[0].x = candidateArray[canIdx].motionVector_x_L0;
            mv_unit.mv[0].y = candidateArray[canIdx].motionVector_y_L0;
            candidateArray[canIdx].local_warp_valid = warped_motion_parameters(
                picture_control_set_ptr,
                context_ptr->cu_ptr,
                &mv_unit,
                context_ptr->blk_geom,
                context_ptr->cu_origin_x,
                context_ptr->cu_origin_y,
                candidateArray[canIdx].ref_frame_type,
                &candidateArray[canIdx].wm_params,
                &candidateArray[canIdx].num_proj_ref);

            if (candidateArray[canIdx].local_warp_valid)
                ++canIdx;
#if REMOVED_DUPLICATE_INTER
            context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
            context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
            ++context_ptr->injected_mv_count_l0;
        }
#endif
    }

    // NEWMV L0
    const MV neighbors[9] = { { 0, 0 },
        { 0, -2 }, { 2, 0 }, { 0, 2 }, { -2, 0 } ,
        { 2, -2 }, { 2, 2 }, { 2, 2 }, { -2, 2 } };

    IntMv  bestPredmv[2] = { {0}, {0} };
    for (int i=0; i<9; i++){
#if REMOVED_DUPLICATE_INTER
        int16_t to_inject_mv_x = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : mePuResult->xMvL0 << 1; // context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.col;
        int16_t to_inject_mv_y = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : mePuResult->yMvL0 << 1; // context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.row;
        to_inject_mv_x += neighbors[i].col;
        to_inject_mv_y += neighbors[i].row;
        if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
        candidateArray[canIdx].type = INTER_MODE;
#if TWO_FAST_LOOP 
        candidateArray[canIdx].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
        candidateArray[canIdx].distortion_ready = 0;
#endif
#else
        candidateArray[canIdx].distortion_ready = 0;
#endif
#if ICOPY
        candidateArray[canIdx].use_intrabc = 0;
#endif
        candidateArray[canIdx].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canIdx].merge_index = 0;
        candidateArray[canIdx].mpm_flag = EB_FALSE;
#endif
        candidateArray[canIdx].prediction_direction[0] = (EbPredDirection)0;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canIdx].is_skip_mode_flag = 0;
#endif
        candidateArray[canIdx].inter_mode = NEWMV;
        candidateArray[canIdx].pred_mode = NEWMV;
        candidateArray[canIdx].motion_mode = WARPED_CAUSAL;
        candidateArray[canIdx].wm_params.wmtype = AFFINE;

        candidateArray[canIdx].is_compound = 0;
        candidateArray[canIdx].is_new_mv = 1;
        candidateArray[canIdx].is_zero_mv = 0;

        candidateArray[canIdx].drl_index = 0;

        // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER
        candidateArray[canIdx].motionVector_x_L0 = to_inject_mv_x;
        candidateArray[canIdx].motionVector_y_L0 = to_inject_mv_y;
#else
        candidateArray[canIdx].motionVector_x_L0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : mePuResult->xMvL0 << 1; // context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.col;
        candidateArray[canIdx].motionVector_y_L0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : mePuResult->yMvL0 << 1; // context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.row;
        candidateArray[canIdx].motionVector_x_L0 += neighbors[i].col;
        candidateArray[canIdx].motionVector_y_L0 += neighbors[i].row;
#endif
        candidateArray[canIdx].ref_mv_index = 0;
        candidateArray[canIdx].pred_mv_weight = 0;
        candidateArray[canIdx].ref_frame_type = LAST_FRAME;

        candidateArray[canIdx].transform_type[PLANE_TYPE_Y] = DCT_DCT;
        candidateArray[canIdx].transform_type[PLANE_TYPE_UV] = DCT_DCT;

        ChooseBestAv1MvPred(
            context_ptr,
            candidateArray[canIdx].md_rate_estimation_ptr,
            context_ptr->cu_ptr,
            candidateArray[canIdx].ref_frame_type,
            candidateArray[canIdx].is_compound,
            candidateArray[canIdx].pred_mode,
            candidateArray[canIdx].motionVector_x_L0,
            candidateArray[canIdx].motionVector_y_L0,
            0, 0,
            &candidateArray[canIdx].drl_index,
            bestPredmv);

        candidateArray[canIdx].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
        candidateArray[canIdx].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;

        MvUnit_t mv_unit;
        mv_unit.mv[0].x = candidateArray[canIdx].motionVector_x_L0;
        mv_unit.mv[0].y = candidateArray[canIdx].motionVector_y_L0;
        candidateArray[canIdx].local_warp_valid = warped_motion_parameters(
            picture_control_set_ptr,
            context_ptr->cu_ptr,
            &mv_unit,
            context_ptr->blk_geom,
            context_ptr->cu_origin_x,
            context_ptr->cu_origin_y,
            candidateArray[canIdx].ref_frame_type,
            &candidateArray[canIdx].wm_params,
            &candidateArray[canIdx].num_proj_ref);

        if (candidateArray[canIdx].local_warp_valid)
            ++canIdx;
#if REMOVED_DUPLICATE_INTER
        context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
        context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
        ++context_ptr->injected_mv_count_l0;
    }
#endif
    }

    *candTotCnt = canIdx;
}


// END of Function Declarations
void  inject_inter_candidates(
    PictureControlSet_t            *picture_control_set_ptr,
    ModeDecisionContext_t          *context_ptr,
    SsMeContext_t                  *ss_mecontext,
    const SequenceControlSet_t     *sequence_control_set_ptr,
    LargestCodingUnit_t            *sb_ptr,
#if M8_SKIP_BLK
    uint32_t                       *candidateTotalCnt){
#else
    uint32_t                       *candidateTotalCnt,
    const uint32_t                  leaf_index){
#endif

    (void)sequence_control_set_ptr;
    uint32_t                   canTotalCnt = *candidateTotalCnt;
    const uint32_t             lcuAddr = sb_ptr->index;
    ModeDecisionCandidate_t    *candidateArray = context_ptr->fast_candidate_array;
    static MvReferenceFrame refFrames[] = { LAST_FRAME, BWDREF_FRAME, LAST_BWD_FRAME };
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    uint32_t me_sb_addr;
    uint32_t geom_offset_x = 0;
    uint32_t geom_offset_y = 0;

    if (sequence_control_set_ptr->sb_size == BLOCK_128X128) {

        uint32_t me_sb_size = sequence_control_set_ptr->sb_sz;
        uint32_t me_pic_width_in_sb = (sequence_control_set_ptr->luma_width + sequence_control_set_ptr->sb_sz - 1) / me_sb_size;
        uint32_t me_sb_x = (context_ptr->cu_origin_x / me_sb_size);
        uint32_t me_sb_y = (context_ptr->cu_origin_y / me_sb_size);

        me_sb_addr = me_sb_x + me_sb_y * me_pic_width_in_sb;

        geom_offset_x = (me_sb_x & 0x1) * me_sb_size;
        geom_offset_y = (me_sb_y & 0x1) * me_sb_size;

    }
    else {
        me_sb_addr = lcuAddr;
    }

    uint32_t max_number_of_pus_per_sb;
#if DISABLE_NSQ_FOR_NON_REF || DISABLE_NSQ

    max_number_of_pus_per_sb = picture_control_set_ptr->parent_pcs_ptr->max_number_of_pus_per_sb;
    
#if DISABLE_IN_LOOP_ME
#if TEST5_DISABLE_NSQ_ME
    const uint32_t me2Nx2NTableOffset = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4 || context_ptr->blk_geom->bwidth == 128 || context_ptr->blk_geom->bheight == 128 || context_ptr->blk_geom->shape != PART_N) ? 0 :
        get_me_info_index(max_number_of_pus_per_sb, context_ptr->blk_geom, geom_offset_x, geom_offset_y);

#else
    uint32_t me2Nx2NTableOffset;

        me2Nx2NTableOffset = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4 || context_ptr->blk_geom->bwidth == 128 || context_ptr->blk_geom->bheight == 128) ? 0 :
            get_me_info_index(max_number_of_pus_per_sb, context_ptr->blk_geom, geom_offset_x, geom_offset_y);
#endif
#else
    const uint32_t me2Nx2NTableOffset = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? 0 :
        get_me_info_index(max_number_of_pus_per_sb, context_ptr->blk_geom, geom_offset_x, geom_offset_y);
#endif
#else
#if DISABLE_IN_LOOP_ME
#if TEST5_DISABLE_NSQ_ME
    const uint32_t me2Nx2NTableOffset = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4 || context_ptr->blk_geom->bwidth == 128 || context_ptr->blk_geom->bheight == 128 || context_ptr->blk_geom->shape != PART_N) ? 0 :
        get_me_info_index(picture_control_set_ptr->parent_pcs_ptr->max_number_of_pus_per_sb, context_ptr->blk_geom, geom_offset_x, geom_offset_y);
#else
    const uint32_t me2Nx2NTableOffset = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4 || context_ptr->blk_geom->bwidth == 128 || context_ptr->blk_geom->bheight == 128) ? 0 :
        get_me_info_index(picture_control_set_ptr->parent_pcs_ptr->max_number_of_pus_per_sb, context_ptr->blk_geom, geom_offset_x, geom_offset_y);
#endif
#else
    const uint32_t me2Nx2NTableOffset = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? 0 :
        get_me_info_index(picture_control_set_ptr->parent_pcs_ptr->max_number_of_pus_per_sb, context_ptr->blk_geom, geom_offset_x, geom_offset_y);
#endif
#endif

    MeCuResults_t * mePuResult = &picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr][me2Nx2NTableOffset];
    EbBool use_close_loop_me = picture_control_set_ptr->parent_pcs_ptr->enable_in_loop_motion_estimation_flag &&
        ((context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) || (context_ptr->blk_geom->bwidth > 64 || context_ptr->blk_geom->bheight > 64)) ? EB_TRUE : EB_FALSE;

    uint32_t close_loop_me_index = use_close_loop_me ? get_in_loop_me_info_index(MAX_SS_ME_PU_COUNT, sequence_control_set_ptr->sb_size == BLOCK_128X128 ? 1 : 0, context_ptr->blk_geom) : 0;
#if BASE_LAYER_REF
    EbBool allow_bipred = (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0 || context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : EB_TRUE;
#else    
    EbBool allow_bipred = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : EB_TRUE;
#endif
    IntMv  bestPredmv[2] = { {0}, {0} };
    uint8_t sq_index = LOG2F(context_ptr->blk_geom->sq_size) - 2;
    uint8_t inject_newmv_candidate = 1;
#if !NSQ_OPTIMASATION
    if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level == NSQ_INTER_SEARCH_BASE_ON_SQ_MVMODE) {
        inject_newmv_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
            context_ptr->parent_sq_pred_mode[sq_index] == NEWMV || context_ptr->parent_sq_pred_mode[sq_index] == NEW_NEWMV ? inject_newmv_candidate : 0;
    }
#endif
#if NSQ_OPTIMASATION
    if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level >= NSQ_SEARCH_LEVEL1 &&
        picture_control_set_ptr->parent_pcs_ptr->nsq_search_level < NSQ_SEARCH_FULL) {
        inject_newmv_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
            context_ptr->parent_sq_has_coeff[sq_index] != 0 ? inject_newmv_candidate : 0;
    }
#else
    if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level == NSQ_SEARCH_BASE_ON_SQ_COEFF) {
        inject_newmv_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
            context_ptr->parent_sq_has_coeff[sq_index] != 0 ? inject_newmv_candidate : 0;
    }
#endif

    generate_av1_mvp_table(
#if TILES
        &sb_ptr->tile_info,
#endif
        context_ptr,
        context_ptr->cu_ptr,
        context_ptr->blk_geom,
        context_ptr->cu_origin_x,
        context_ptr->cu_origin_y,
        refFrames,
        (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 1 : 3,
        picture_control_set_ptr);

    uint32_t mi_row = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t mi_col = context_ptr->cu_origin_x >> MI_SIZE_LOG2;
    av1_count_overlappable_neighbors(
        picture_control_set_ptr,
        context_ptr->cu_ptr,
        context_ptr->blk_geom->bsize,
        mi_row,
        mi_col);

    /**************
         MVP
    ************* */
    InjectAv1MvpCandidates(
        context_ptr,
        context_ptr->cu_ptr,
        refFrames,
        picture_control_set_ptr,
        lcuAddr,
#if !M8_SKIP_BLK
        leaf_index,
#endif
        allow_bipred,
        &canTotalCnt);

    if (inject_newmv_candidate) {
        /**************
            NEWMV L0
        ************* */
#if REMOVED_DUPLICATE_INTER
        int16_t to_inject_mv_x = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : mePuResult->xMvL0 << 1;
        int16_t to_inject_mv_y = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : mePuResult->yMvL0 << 1;
        if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
        candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
        candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
        candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
        candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
        candidateArray[canTotalCnt].use_intrabc = 0;
#endif
        candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canTotalCnt].merge_index = 0;
        candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
        candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;
#if !INTRA_INTER_FAST_LOOP
        candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif
        candidateArray[canTotalCnt].inter_mode = NEWMV;
        candidateArray[canTotalCnt].pred_mode = NEWMV;
        candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

        candidateArray[canTotalCnt].is_compound = 0;
        candidateArray[canTotalCnt].is_new_mv = 1;
        candidateArray[canTotalCnt].is_zero_mv = 0;

        candidateArray[canTotalCnt].drl_index = 0;

        // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER
        candidateArray[canTotalCnt].motionVector_x_L0 = to_inject_mv_x;
        candidateArray[canTotalCnt].motionVector_y_L0 = to_inject_mv_y;
#else
        candidateArray[canTotalCnt].motionVector_x_L0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : mePuResult->xMvL0 << 1;
        candidateArray[canTotalCnt].motionVector_y_L0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : mePuResult->yMvL0 << 1;
#endif
        // will be needed later by the rate estimation
        candidateArray[canTotalCnt].ref_mv_index = 0;
        candidateArray[canTotalCnt].pred_mv_weight = 0;
        candidateArray[canTotalCnt].ref_frame_type = LAST_FRAME;


        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

        ChooseBestAv1MvPred(
            context_ptr,
            candidateArray[canTotalCnt].md_rate_estimation_ptr,
            context_ptr->cu_ptr,
            candidateArray[canTotalCnt].ref_frame_type,
            candidateArray[canTotalCnt].is_compound,
            candidateArray[canTotalCnt].pred_mode,
            candidateArray[canTotalCnt].motionVector_x_L0,
            candidateArray[canTotalCnt].motionVector_y_L0,
            0, 0,
            &candidateArray[canTotalCnt].drl_index,
            bestPredmv);

        candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
        candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;

        ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER
        context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
        context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
        ++context_ptr->injected_mv_count_l0;
        }
#endif
        if (isCompoundEnabled) {
            /**************
               NEWMV L1
           ************* */
#if REMOVED_DUPLICATE_INTER_L1
            int16_t to_inject_mv_x = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : mePuResult->xMvL1 << 1;
            int16_t to_inject_mv_y = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : mePuResult->yMvL1 << 1;
            if (context_ptr->injected_mv_count_l1 == 0 || is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
            candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
            candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
            candidateArray[canTotalCnt].use_intrabc = 0;
#endif
            candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].merge_index = 0;
            candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)1;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif
            candidateArray[canTotalCnt].inter_mode = NEWMV;
            candidateArray[canTotalCnt].pred_mode = NEWMV;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

            candidateArray[canTotalCnt].is_compound = 0;
            candidateArray[canTotalCnt].is_new_mv = 1;
            candidateArray[canTotalCnt].is_zero_mv = 0;

            candidateArray[canTotalCnt].drl_index = 0;

            // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER_L1
            candidateArray[canTotalCnt].motionVector_x_L1 = to_inject_mv_x;
            candidateArray[canTotalCnt].motionVector_y_L1 = to_inject_mv_y;
#else
            candidateArray[canTotalCnt].motionVector_x_L1 = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : mePuResult->xMvL1 << 1;
            candidateArray[canTotalCnt].motionVector_y_L1 = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : mePuResult->yMvL1 << 1;
#endif
            // will be needed later by the rate estimation
            candidateArray[canTotalCnt].ref_mv_index = 0;
            candidateArray[canTotalCnt].pred_mv_weight = 0;
            candidateArray[canTotalCnt].ref_frame_type = BWDREF_FRAME;


            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

            ChooseBestAv1MvPred(
                context_ptr,
                candidateArray[canTotalCnt].md_rate_estimation_ptr,
                context_ptr->cu_ptr,
                candidateArray[canTotalCnt].ref_frame_type,
                candidateArray[canTotalCnt].is_compound,
                candidateArray[canTotalCnt].pred_mode,
                candidateArray[canTotalCnt].motionVector_x_L1,
                candidateArray[canTotalCnt].motionVector_y_L1,
                0, 0,
                &candidateArray[canTotalCnt].drl_index,
                bestPredmv);

            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[0].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[0].as_mv.row;

            ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER_L1
            context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
            context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
            ++context_ptr->injected_mv_count_l1;
            }
#endif
            /**************
               NEW_NEWMV
            ************* */
            if (allow_bipred) {

#if REMOVED_DUPLICATE_INTER_BIPRED
                int16_t to_inject_mv_x_l0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : mePuResult->xMvL0 << 1;
                int16_t to_inject_mv_y_l0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : mePuResult->yMvL0 << 1;
                int16_t to_inject_mv_x_l1 = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : mePuResult->xMvL1 << 1;
                int16_t to_inject_mv_y_l1 = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : mePuResult->yMvL1 << 1;
                if (context_ptr->injected_mv_count_bipred == 0 || is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1) == EB_FALSE) {
#endif
                candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
                candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
                candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
                candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
                candidateArray[canTotalCnt].use_intrabc = 0;
#endif
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canTotalCnt].merge_index = 0;
                candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
                candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif

                candidateArray[canTotalCnt].is_new_mv = 1;
                candidateArray[canTotalCnt].is_zero_mv = 0;

                candidateArray[canTotalCnt].drl_index = 0;

                // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER_BIPRED
                candidateArray[canTotalCnt].motionVector_x_L0 = to_inject_mv_x_l0;
                candidateArray[canTotalCnt].motionVector_y_L0 = to_inject_mv_y_l0;
                candidateArray[canTotalCnt].motionVector_x_L1 = to_inject_mv_x_l1;
                candidateArray[canTotalCnt].motionVector_y_L1 = to_inject_mv_y_l1;
#else
                candidateArray[canTotalCnt].motionVector_x_L0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : mePuResult->xMvL0 << 1;
                candidateArray[canTotalCnt].motionVector_y_L0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : mePuResult->yMvL0 << 1;
                candidateArray[canTotalCnt].motionVector_x_L1 = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : mePuResult->xMvL1 << 1;
                candidateArray[canTotalCnt].motionVector_y_L1 = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : mePuResult->yMvL1 << 1;
#endif
                // will be needed later by the rate estimation
                candidateArray[canTotalCnt].ref_mv_index = 0;
                candidateArray[canTotalCnt].pred_mv_weight = 0;

                candidateArray[canTotalCnt].inter_mode = NEW_NEWMV;
                candidateArray[canTotalCnt].pred_mode = NEW_NEWMV;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canTotalCnt].is_compound = 1;
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
                candidateArray[canTotalCnt].ref_frame_type = LAST_BWD_FRAME;

                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

                ChooseBestAv1MvPred(
                    context_ptr,
                    candidateArray[canTotalCnt].md_rate_estimation_ptr,
                    context_ptr->cu_ptr,
                    candidateArray[canTotalCnt].ref_frame_type,
                    candidateArray[canTotalCnt].is_compound,
                    candidateArray[canTotalCnt].pred_mode,
                    candidateArray[canTotalCnt].motionVector_x_L0,
                    candidateArray[canTotalCnt].motionVector_y_L0,
                    candidateArray[canTotalCnt].motionVector_x_L1,
                    candidateArray[canTotalCnt].motionVector_y_L1,
                    &candidateArray[canTotalCnt].drl_index,
                    bestPredmv);

                candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
                candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
                candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;

                ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER_BIPRED
                context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                ++context_ptr->injected_mv_count_bipred;
                }
#endif
            }

        }
#if M0_ME_SEARCH_BASE && !BASE_LAYER_REF
        /**************
        inject NewMv from L1 as a candidate of NEWMV L0 in base layer frame (Only single reference support in base layer)
        ************* */
        if (picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0) {
#if REMOVED_DUPLICATE_INTER
            int16_t to_inject_mv_x = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : mePuResult->xMvL1 << 1;
            int16_t to_inject_mv_y = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : mePuResult->yMvL1 << 1;
            if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
            candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
            candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
            candidateArray[canTotalCnt].use_intrabc = 0;
#endif
            candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].merge_index = 0;
            candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif
            candidateArray[canTotalCnt].inter_mode = NEWMV;
            candidateArray[canTotalCnt].pred_mode = NEWMV;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

            candidateArray[canTotalCnt].is_compound = 0;
            candidateArray[canTotalCnt].is_new_mv = 1;
            candidateArray[canTotalCnt].is_zero_mv = 0;

            candidateArray[canTotalCnt].drl_index = 0;

            // Set the MV to ME result
#if REMOVED_DUPLICATE_INTER
            candidateArray[canTotalCnt].motionVector_x_L0 = to_inject_mv_x;
            candidateArray[canTotalCnt].motionVector_y_L0 = to_inject_mv_y;
#else
            candidateArray[canTotalCnt].motionVector_x_L0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : mePuResult->xMvL1 << 1;
            candidateArray[canTotalCnt].motionVector_y_L0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : mePuResult->yMvL1 << 1;
#endif

            // will be needed later by the rate estimation
            candidateArray[canTotalCnt].ref_mv_index = 0;
            candidateArray[canTotalCnt].pred_mv_weight = 0;
            candidateArray[canTotalCnt].ref_frame_type = LAST_FRAME;


            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

            ChooseBestAv1MvPred(
                context_ptr,
                candidateArray[canTotalCnt].md_rate_estimation_ptr,
                context_ptr->cu_ptr,
                candidateArray[canTotalCnt].ref_frame_type,
                candidateArray[canTotalCnt].is_compound,
                candidateArray[canTotalCnt].pred_mode,
                candidateArray[canTotalCnt].motionVector_x_L0,
                candidateArray[canTotalCnt].motionVector_y_L0,
                0, 0,
                &candidateArray[canTotalCnt].drl_index,
                bestPredmv);

            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;

            ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER
            context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
            context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
            ++context_ptr->injected_mv_count_l0;
            }
#endif
        }
#endif
    }

#if SHUT_GLOBAL_MV
    if (context_ptr->global_mv_injection) {
#endif
        /**************
         GLOBALMV L0
        ************* */
        {

#if REMOVED_DUPLICATE_INTER
            int16_t to_inject_mv_x = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_y = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
            if (context_ptr->injected_mv_count_l0 == 0 || is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y) == EB_FALSE) {
#endif
                candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
                candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
                candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
                candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
                candidateArray[canTotalCnt].use_intrabc = 0;
#endif
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canTotalCnt].merge_index = 0;
                candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif
                candidateArray[canTotalCnt].inter_mode = GLOBALMV;
                candidateArray[canTotalCnt].pred_mode = GLOBALMV;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canTotalCnt].is_compound = 0;
                candidateArray[canTotalCnt].is_new_mv = 0;
                candidateArray[canTotalCnt].is_zero_mv = 0;
                candidateArray[canTotalCnt].drl_index = 0;

                // will be needed later by the rate estimation
                candidateArray[canTotalCnt].ref_mv_index = 0;
                candidateArray[canTotalCnt].pred_mv_weight = 0;
                candidateArray[canTotalCnt].ref_frame_type = LAST_FRAME;

                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

                // Set the MV to frame MV
#if REMOVED_DUPLICATE_INTER
                candidateArray[canTotalCnt].motionVector_x_L0 = to_inject_mv_x;
                candidateArray[canTotalCnt].motionVector_y_L0 = to_inject_mv_y;
#else
                candidateArray[canTotalCnt].motionVector_y_L0 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
                candidateArray[canTotalCnt].motionVector_x_L0 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
#endif

                ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER
                context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
                context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
                ++context_ptr->injected_mv_count_l0;
            }
#endif
        }

        if (isCompoundEnabled && allow_bipred) {

            /**************
            GLOBAL_GLOBALMV
            ************* */
#if REMOVED_DUPLICATE_INTER_BIPRED
            int16_t to_inject_mv_x_l0 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_y_l0 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_x_l1 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_y_l1 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
            if (context_ptr->injected_mv_count_bipred == 0 || is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1) == EB_FALSE) {
#endif
                candidateArray[canTotalCnt].type = INTER_MODE;
#if TWO_FAST_LOOP 
                candidateArray[canTotalCnt].enable_two_fast_loops = 0;
#if OIS_BASED_INTRA
                candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
                candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
                candidateArray[canTotalCnt].use_intrabc = 0;
#endif
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canTotalCnt].merge_index = 0;
                candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
#if !INTRA_INTER_FAST_LOOP
                candidateArray[canTotalCnt].is_skip_mode_flag = 0;
#endif
                candidateArray[canTotalCnt].inter_mode = GLOBAL_GLOBALMV;
                candidateArray[canTotalCnt].pred_mode = GLOBAL_GLOBALMV;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canTotalCnt].is_compound = 1;
                candidateArray[canTotalCnt].is_new_mv = 0;
                candidateArray[canTotalCnt].is_zero_mv = 0;
                candidateArray[canTotalCnt].drl_index = 0;

                // will be needed later by the rate estimation
                candidateArray[canTotalCnt].ref_mv_index = 0;
                candidateArray[canTotalCnt].pred_mv_weight = 0;
                candidateArray[canTotalCnt].ref_frame_type = LAST_BWD_FRAME;

                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

                // Set the MV to frame MV
#if REMOVED_DUPLICATE_INTER_BIPRED
                candidateArray[canTotalCnt].motionVector_x_L0 = to_inject_mv_x_l0;
                candidateArray[canTotalCnt].motionVector_y_L0 = to_inject_mv_y_l0;
                candidateArray[canTotalCnt].motionVector_x_L1 = to_inject_mv_x_l1;
                candidateArray[canTotalCnt].motionVector_y_L1 = to_inject_mv_y_l1;
#else
                candidateArray[canTotalCnt].motionVector_y_L0 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
                candidateArray[canTotalCnt].motionVector_x_L0 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);


                // Set the MV to frame MV
                candidateArray[canTotalCnt].motionVector_y_L1 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
                candidateArray[canTotalCnt].motionVector_x_L1 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
#endif

                ++canTotalCnt;
#if REMOVED_DUPLICATE_INTER_BIPRED
                context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                ++context_ptr->injected_mv_count_bipred;
            }
#endif
        }
#if SHUT_GLOBAL_MV
    }
#endif
    // Warped Motion
    if (picture_control_set_ptr->parent_pcs_ptr->allow_warped_motion &&
        has_overlappable_candidates(context_ptr->cu_ptr) &&
        context_ptr->blk_geom->bwidth >= 8 &&
        context_ptr->blk_geom->bheight >= 8 &&
        context_ptr->warped_motion_injection){

        inject_warped_motion_candidates(
            picture_control_set_ptr,
            context_ptr,
            context_ptr->cu_ptr,
            &canTotalCnt,
            ss_mecontext,
            mePuResult,
            use_close_loop_me,
            close_loop_me_index);
    }

    if (inject_newmv_candidate) {
#if BASE_LAYER_REF
        if (isCompoundEnabled) {
            if (allow_bipred) {
#else
        if (allow_bipred) {
#endif
#if IMPROVED_BIPRED_INJECTION
            //----------------------
            // Bipred2Nx2N
            //----------------------
            if (context_ptr->bipred3x3_injection > 0)
                if (picture_control_set_ptr->slice_type == B_SLICE)
                    Bipred3x3CandidatesInjection(
                        picture_control_set_ptr,
                        context_ptr,
                        sb_ptr,
                        me_sb_addr,
                        ss_mecontext,
                        use_close_loop_me,
                        close_loop_me_index,
                        me2Nx2NTableOffset,
                        &canTotalCnt);
#endif
#if BASE_LAYER_REF
            }
#endif
#if IMPROVED_UNIPRED_INJECTION
            //----------------------
            // Unipred2Nx2N
            //----------------------
            if (context_ptr->unipred3x3_injection > 0)
                if (picture_control_set_ptr->slice_type != I_SLICE)
                    Unipred3x3CandidatesInjection(
                        picture_control_set_ptr,
                        context_ptr,
                        sb_ptr,
                        me_sb_addr,
                        ss_mecontext,
                        use_close_loop_me,
                        close_loop_me_index,
                        me2Nx2NTableOffset,
                        &canTotalCnt);
#endif
        }
    }
    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;


    return;
}




static INLINE PredictionMode get_uv_mode(UV_PredictionMode mode) {
    assert(mode < UV_INTRA_MODES);
    static const PredictionMode uv2y[] = {
        DC_PRED,        // UV_DC_PRED
        V_PRED,         // UV_V_PRED
        H_PRED,         // UV_H_PRED
        D45_PRED,       // UV_D45_PRED
        D135_PRED,      // UV_D135_PRED
        D113_PRED,      // UV_D113_PRED
        D157_PRED,      // UV_D157_PRED
        D203_PRED,      // UV_D203_PRED
        D67_PRED,       // UV_D67_PRED
        SMOOTH_PRED,    // UV_SMOOTH_PRED
        SMOOTH_V_PRED,  // UV_SMOOTH_V_PRED
        SMOOTH_H_PRED,  // UV_SMOOTH_H_PRED
        PAETH_PRED,     // UV_PAETH_PRED
        DC_PRED,        // UV_CFL_PRED
        INTRA_INVALID,  // UV_INTRA_MODES
        INTRA_INVALID,  // UV_MODE_INVALID
    };
    return uv2y[mode];
}
static TxType intra_mode_to_tx_type(const MbModeInfo *mbmi,
    PLANE_TYPE plane_type) {
    static const TxType _intra_mode_to_tx_type[INTRA_MODES] = {
        DCT_DCT,    // DC
        ADST_DCT,   // V
        DCT_ADST,   // H
        DCT_DCT,    // D45
        ADST_ADST,  // D135
        ADST_DCT,   // D117
        DCT_ADST,   // D153
        DCT_ADST,   // D207
        ADST_DCT,   // D63
        ADST_ADST,  // SMOOTH
        ADST_DCT,   // SMOOTH_V
        DCT_ADST,   // SMOOTH_H
        ADST_ADST,  // PAETH
    };
    const PredictionMode mode =
        (plane_type == PLANE_TYPE_Y) ? mbmi->mode : get_uv_mode(mbmi->uv_mode);
    assert(mode < INTRA_MODES);
    return _intra_mode_to_tx_type[mode];
}

static INLINE TxType av1_get_tx_type(
    block_size  sb_type,
    int32_t   is_inter,
    PredictionMode pred_mode,
    UV_PredictionMode pred_mode_uv,
    PLANE_TYPE plane_type,
    const MacroBlockD *xd, int32_t blk_row,
    int32_t blk_col, TxSize tx_size,
    int32_t reduced_tx_set)
{
    UNUSED(sb_type);
    UNUSED(*xd);
    UNUSED(blk_row);
    UNUSED(blk_col);

    // block_size  sb_type = BLOCK_8X8;

    MbModeInfo  mbmi;
    mbmi.mode = pred_mode;
    mbmi.uv_mode = pred_mode_uv;


    // const MbModeInfo *const mbmi = xd->mi[0];
    // const struct MacroblockdPlane *const pd = &xd->plane[plane_type];
    const TxSetType tx_set_type =
        /*av1_*/get_ext_tx_set_type(tx_size, is_inter, reduced_tx_set);

    TxType tx_type = DCT_DCT;
    if ( /*xd->lossless[mbmi->segment_id] ||*/ txsize_sqr_up_map[tx_size] > TX_32X32) {
        tx_type = DCT_DCT;
    }
    else {
        if (plane_type == PLANE_TYPE_Y) {
            //const int32_t txk_type_idx =
            //    av1_get_txk_type_index(/*mbmi->*/sb_type, blk_row, blk_col);
            //tx_type = mbmi->txk_type[txk_type_idx];
        }
        else if (is_inter /*is_inter_block(mbmi)*/) {
            // scale back to y plane's coordinate
            //blk_row <<= pd->subsampling_y;
            //blk_col <<= pd->subsampling_x;
            //const int32_t txk_type_idx =
            //    av1_get_txk_type_index(mbmi->sb_type, blk_row, blk_col);
            //tx_type = mbmi->txk_type[txk_type_idx];
        }
        else {
            // In intra mode, uv planes don't share the same prediction mode as y
            // plane, so the tx_type should not be shared
            tx_type = intra_mode_to_tx_type(&mbmi, PLANE_TYPE_UV);
        }
    }
    ASSERT(tx_type < TX_TYPES);
    if (!av1_ext_tx_used[tx_set_type][tx_type]) return DCT_DCT;
    return tx_type;
}


#if OIS_BASED_INTRA

void  inject_intra_candidates_ois(
    PictureControlSet_t            *picture_control_set_ptr,
    ModeDecisionContext_t          *context_ptr,
    LargestCodingUnit_t            *sb_ptr,
    uint32_t                       *candidate_total_cnt){

    uint8_t                     intra_candidate_counter;
    uint8_t                     intra_mode;
    uint32_t                    can_total_cnt = 0;
    EbBool                      use_angle_delta = (context_ptr->blk_geom->bsize >= BLOCK_8X8);
    ModeDecisionCandidate_t    *candidate_array = context_ptr->fast_candidate_array;

    EbBool                      disable_cfl_flag = (context_ptr->blk_geom->sq_size > 32 || 
                                                    context_ptr->blk_geom->bwidth == 4  ||   
                                                    context_ptr->blk_geom->bheight == 4)    ? EB_TRUE : EB_FALSE;


    ois_sb_results_t    *ois_sb_results_ptr = picture_control_set_ptr->parent_pcs_ptr->ois_sb_results[sb_ptr->index];
    ois_candidate_t     *ois_blk_ptr = ois_sb_results_ptr->ois_candidate_array[ep_to_pa_block_index[context_ptr->blk_geom->blkidx_mds]];
    uint8_t              total_intra_luma_mode = ois_sb_results_ptr-> total_ois_intra_candidate[ep_to_pa_block_index[context_ptr->blk_geom->blkidx_mds]];
    
#if TWO_FAST_LOOP 
    uint8_t enable_two_fast_loops = picture_control_set_ptr->parent_pcs_ptr->enable_two_fast_loops && (context_ptr->blk_geom->sq_size > 4 && context_ptr->blk_geom->shape == PART_N);          
#endif

    for (intra_candidate_counter = 0; intra_candidate_counter < total_intra_luma_mode; ++intra_candidate_counter) {
                

        intra_mode = ois_blk_ptr[can_total_cnt].intra_mode;
        if (av1_is_directional_mode((PredictionMode)intra_mode)) {

            int32_t angle_delta = ois_blk_ptr[can_total_cnt].angle_delta ;
            candidate_array[can_total_cnt].type = INTRA_MODE;
            candidate_array[can_total_cnt].intra_luma_mode = intra_mode;
#if TWO_FAST_LOOP 
            candidate_array[can_total_cnt].enable_two_fast_loops =  enable_two_fast_loops;
#endif
            candidate_array[can_total_cnt].distortion_ready =  1;
            candidate_array[can_total_cnt].me_distortion = ois_blk_ptr[can_total_cnt].distortion;
#if ICOPY
            candidate_array[can_total_cnt].use_intrabc = 0;
#endif
            candidate_array[can_total_cnt].is_directional_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)intra_mode);
            candidate_array[can_total_cnt].use_angle_delta = use_angle_delta ? candidate_array[can_total_cnt].is_directional_mode_flag : 0;
            candidate_array[can_total_cnt].angle_delta[PLANE_TYPE_Y] = angle_delta;
#if CHROMA_BLIND 
            candidate_array[can_total_cnt].intra_chroma_mode = disable_cfl_flag ? intra_luma_to_chroma[intra_mode] : 
                                                               context_ptr->chroma_level == CHROMA_MODE_0 ? UV_CFL_PRED : UV_DC_PRED;
#else
            candidate_array[can_total_cnt].intra_chroma_mode = disable_cfl_flag ? intra_luma_to_chroma[openLoopIntraCandidate] : UV_CFL_PRED;
#endif
            candidate_array[can_total_cnt].cfl_alpha_signs = 0;
            candidate_array[can_total_cnt].cfl_alpha_idx = 0;
            candidate_array[can_total_cnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidate_array[can_total_cnt].intra_chroma_mode);
            candidate_array[can_total_cnt].angle_delta[PLANE_TYPE_UV] = 0;
            candidate_array[can_total_cnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;

            if (candidate_array[can_total_cnt].intra_chroma_mode == UV_CFL_PRED)
                candidate_array[can_total_cnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;
            else
                candidate_array[can_total_cnt].transform_type[PLANE_TYPE_UV] =
                av1_get_tx_type(
                    context_ptr->blk_geom->bsize,
                    0,
                    (PredictionMode)candidate_array[can_total_cnt].intra_luma_mode,
                    (UV_PredictionMode)candidate_array[can_total_cnt].intra_chroma_mode,
                    PLANE_TYPE_UV,
                    0,
                    0,
                    0,
                    context_ptr->blk_geom->txsize_uv[0],
                    picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
#if !INTRA_INTER_FAST_LOOP
            // candidate_array[can_total_cnt].transform_type[PLANE_TYPE_UV]            = context_ptr->blk_geom->bwidth_uv >= 32 || context_ptr->blk_geom->bheight_uv >= 32  ? DCT_DCT : chroma_transform_type[candidate_array[can_total_cnt].intra_chroma_mode];
            candidate_array[can_total_cnt].mpm_flag = EB_FALSE;
#endif
            candidate_array[can_total_cnt].ref_frame_type = INTRA_FRAME;
            candidate_array[can_total_cnt].pred_mode = (PredictionMode)intra_mode;
            candidate_array[can_total_cnt].motion_mode = SIMPLE_TRANSLATION;
            ++can_total_cnt;

        }
        else {
            candidate_array[can_total_cnt].type = INTRA_MODE;
            candidate_array[can_total_cnt].intra_luma_mode = intra_mode;
#if TWO_FAST_LOOP 
            candidate_array[can_total_cnt].enable_two_fast_loops =  enable_two_fast_loops;
#endif
            candidate_array[can_total_cnt].distortion_ready =  1;
            candidate_array[can_total_cnt].me_distortion = ois_blk_ptr[can_total_cnt].distortion;
#if ICOPY
            candidate_array[can_total_cnt].use_intrabc = 0;
#endif
            candidate_array[can_total_cnt].is_directional_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)intra_mode);
            candidate_array[can_total_cnt].use_angle_delta = candidate_array[can_total_cnt].is_directional_mode_flag;
            candidate_array[can_total_cnt].angle_delta[PLANE_TYPE_Y] = 0;
#if TURN_OFF_CFL
            candidate_array[can_total_cnt].intra_chroma_mode = intra_luma_to_chroma[openLoopIntraCandidate];
#else
#if CHROMA_BLIND
            candidate_array[can_total_cnt].intra_chroma_mode =  disable_cfl_flag ? intra_luma_to_chroma[intra_mode] : 
                                                                context_ptr->chroma_level == CHROMA_MODE_0 ? UV_CFL_PRED : UV_DC_PRED;

#else
            candidate_array[can_total_cnt].intra_chroma_mode = disable_cfl_flag ? intra_luma_to_chroma[openLoopIntraCandidate] : UV_CFL_PRED;
#endif
#endif
            candidate_array[can_total_cnt].cfl_alpha_signs = 0;
            candidate_array[can_total_cnt].cfl_alpha_idx = 0;
            candidate_array[can_total_cnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidate_array[can_total_cnt].intra_chroma_mode);
            candidate_array[can_total_cnt].angle_delta[PLANE_TYPE_UV] = 0;
            candidate_array[can_total_cnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;

            if (candidate_array[can_total_cnt].intra_chroma_mode == UV_CFL_PRED)
                candidate_array[can_total_cnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;
            else
                candidate_array[can_total_cnt].transform_type[PLANE_TYPE_UV] =
                av1_get_tx_type(
                    context_ptr->blk_geom->bsize,
                    0,
                    (PredictionMode)candidate_array[can_total_cnt].intra_luma_mode,
                    (UV_PredictionMode)candidate_array[can_total_cnt].intra_chroma_mode,
                    PLANE_TYPE_UV,
                    0,
                    0,
                    0,
                    context_ptr->blk_geom->txsize_uv[0],
                    picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
#if !INTRA_INTER_FAST_LOOP
            candidate_array[can_total_cnt].mpm_flag = EB_FALSE;
#endif
            candidate_array[can_total_cnt].ref_frame_type = INTRA_FRAME;
            candidate_array[can_total_cnt].pred_mode = (PredictionMode)intra_mode;
            candidate_array[can_total_cnt].motion_mode = SIMPLE_TRANSLATION;
            ++can_total_cnt;
        }
    }


    // update the total number of candidates injected
    (*candidate_total_cnt) = can_total_cnt;


    return;
}

#endif

#if ICOPY
double av1_convert_qindex_to_q(int32_t qindex, aom_bit_depth_t bit_depth);

static INLINE void setup_pred_plane(struct buf_2d *dst, block_size bsize,
    uint8_t *src, int width, int height,
    int stride, int mi_row, int mi_col,
    int subsampling_x, int subsampling_y) {
    // Offset the buffer pointer
    if (subsampling_y && (mi_row & 0x01) && (mi_size_high[bsize] == 1))
        mi_row -= 1;
    if (subsampling_x && (mi_col & 0x01) && (mi_size_wide[bsize] == 1))
        mi_col -= 1;

    const int x = (MI_SIZE * mi_col) >> subsampling_x;
    const int y = (MI_SIZE * mi_row) >> subsampling_y;
    dst->buf = src + (y * stride + x);// scaled_buffer_offset(x, y, stride, scale);
    dst->buf0 = src;
    dst->width = width;
    dst->height = height;
    dst->stride = stride;
}
void av1_setup_pred_block(block_size sb_type,
    struct buf_2d dst[MAX_MB_PLANE],
    const Yv12BufferConfig *src, int mi_row, int mi_col) {
    int i;

    dst[0].buf = src->y_buffer;
    dst[0].stride = src->y_stride;
    dst[1].buf = src->u_buffer;
    dst[2].buf = src->v_buffer;
    dst[1].stride = dst[2].stride = src->uv_stride;

    i = 0;
    setup_pred_plane(dst + i, sb_type, dst[i].buf,
        i ? src->uv_crop_width : src->y_crop_width,
        i ? src->uv_crop_height : src->y_crop_height,
        dst[i].stride, mi_row, mi_col,
        0, 0);

}
// Values are now correlated to quantizer.
static int sad_per_bit16lut_8[QINDEX_RANGE];
static int sad_per_bit4lut_8[QINDEX_RANGE];

static void init_me_luts_bd(int *bit16lut, int *bit4lut, int range,
    aom_bit_depth_t bit_depth) {
    int i;
    // Initialize the sad lut tables using a formulaic calculation for now.
    // This is to make it easier to resolve the impact of experimental changes
    // to the quantizer tables.
    for (i = 0; i < range; i++) {
        const double q = av1_convert_qindex_to_q(i, bit_depth);
        bit16lut[i] = (int)(0.0418 * q + 2.4107);
        bit4lut[i] = (int)(0.063 * q + 2.742);
    }
}

void av1_init_me_luts(void) {
    init_me_luts_bd(sad_per_bit16lut_8, sad_per_bit4lut_8, QINDEX_RANGE,
        AOM_BITS_8);
}

static INLINE int mv_check_bounds(const MvLimits *mv_limits, const MV *mv) {
    return (mv->row >> 3) < mv_limits->row_min ||
        (mv->row >> 3) > mv_limits->row_max ||
        (mv->col >> 3) < mv_limits->col_min ||
        (mv->col >> 3) > mv_limits->col_max;
}
void assert_release(int statement)
{
    if (statement == 0)
        printf("ASSERT_ERRRR\n");
}

void  intra_bc_search(
    PictureControlSet_t            *pcs,
    ModeDecisionContext_t          *context_ptr,
    const SequenceControlSet_t     *scs,
    LargestCodingUnit_t            *sb_ptr,
    CodingUnit_t                   *cu_ptr,
    MV                             *dv_cand,
    uint8_t                        *num_dv_cand)
{

    IntraBcContext  x_st;
    IntraBcContext  *x = &x_st;
    //fill x with what needed.
#if IBC_EARLY_0
    x->is_exhaustive_allowed =  context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4 ? 1 : 0;
#endif
#if HASH_X
    //CHKN crc calculator could be moved to mdContext and these init at init time.
    av1_crc_calculator_init(&x->crc_calculator1, 24, 0x5D6DCB);
    av1_crc_calculator_init(&x->crc_calculator2, 24, 0x864CFB);
#endif

    x->xd = cu_ptr->av1xd;
    x->nmv_vec_cost = context_ptr->md_rate_estimation_ptr->nmv_vec_cost;
    x->mv_cost_stack = context_ptr->md_rate_estimation_ptr->nmvcoststack;
    block_size bsize = context_ptr->blk_geom->bsize;
    const Av1Common *const cm = pcs->parent_pcs_ptr->av1_cm;
    MvReferenceFrame ref_frame = INTRA_FRAME;
    generate_av1_mvp_table(
#if TILES
        &sb_ptr->tile_info,
#endif
        context_ptr,
        context_ptr->cu_ptr,
        context_ptr->blk_geom,
        context_ptr->cu_origin_x,
        context_ptr->cu_origin_y,
        &ref_frame,
        1,
        pcs);

    const int num_planes = 3;
    MacroBlockD * xd = cu_ptr->av1xd;
    const TileInfo *tile = &xd->tile;
    const int mi_row = -xd->mb_to_top_edge / (8 * MI_SIZE);
    const int mi_col = -xd->mb_to_left_edge / (8 * MI_SIZE);
    const int w = block_size_wide[bsize];
    const int h = block_size_high[bsize];
    const int sb_row = mi_row >> scs->mib_size_log2;
    const int sb_col = mi_col >> scs->mib_size_log2;

    // Set up limit values for MV components.
    // Mv beyond the range do not produce new/different prediction block.
    const int mi_width = mi_size_wide[bsize];
    const int mi_height = mi_size_high[bsize];
    x->mv_limits.row_min =
        -(((mi_row + mi_height) * MI_SIZE) + AOM_INTERP_EXTEND);
    x->mv_limits.col_min = -(((mi_col + mi_width) * MI_SIZE) + AOM_INTERP_EXTEND);
    x->mv_limits.row_max = (cm->mi_rows - mi_row) * MI_SIZE + AOM_INTERP_EXTEND;
    x->mv_limits.col_max = (cm->mi_cols - mi_col) * MI_SIZE + AOM_INTERP_EXTEND;
    //set search paramters
    x->sadperbit16 = sad_per_bit16lut_8[pcs->parent_pcs_ptr->base_qindex];
    x->errorperbit = context_ptr->full_lambda >> RD_EPB_SHIFT;
    x->errorperbit += (x->errorperbit == 0);
    //temp buffer for hash me
    for (int xi = 0; xi < 2; xi++)
        for (int yj = 0; yj < 2; yj++)
            x->hash_value_buffer[xi][yj] = (uint32_t*)malloc(AOM_BUFFER_SIZE_FOR_BLOCK_HASH * sizeof(uint32_t));

    IntMv nearestmv, nearmv;
    av1_find_best_ref_mvs_from_stack(0, context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack /*mbmi_ext*/, xd, ref_frame, &nearestmv, &nearmv,
        0);
    if (nearestmv.as_int == INVALID_MV) {
        nearestmv.as_int = 0;
    }
    if (nearmv.as_int == INVALID_MV) {
        nearmv.as_int = 0;
    }
    IntMv dv_ref = nearestmv.as_int == 0 ? nearmv : nearestmv;
    if (dv_ref.as_int == 0)
        av1_find_ref_dv(&dv_ref, tile, scs->mib_size, mi_row, mi_col);
    // Ref DV should not have sub-pel.
    assert((dv_ref.as_mv.col & 7) == 0);
    assert((dv_ref.as_mv.row & 7) == 0);
    context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[INTRA_FRAME][0].this_mv = dv_ref;

    /* pointer to current frame */
    Yv12BufferConfig cur_buf;
#if ICOPY_10B
    link_Eb_to_aom_buffer_desc_8bit(
        pcs->parent_pcs_ptr->enhanced_picture_ptr,
        &cur_buf);
#else
    LinkEbToAomBufferDesc(
        pcs->parent_pcs_ptr->enhanced_picture_ptr,
        &cur_buf);
#endif
    struct buf_2d yv12_mb[MAX_MB_PLANE];
    av1_setup_pred_block(bsize, yv12_mb, &cur_buf, mi_row, mi_col);
    for (int i = 0; i < num_planes; ++i) {
        x->xdplane[i].pre[0] = yv12_mb[i];  //ref in ME
    }
    //setup src for DV search same as ref
    x->plane[0].src = x->xdplane[0].pre[0];


    enum IntrabcMotionDirection {
        IBC_MOTION_ABOVE,
        IBC_MOTION_LEFT,
        IBC_MOTION_DIRECTIONS
    };

    //up to two dv candidates will be generated
#if IBC_MODES
    enum IntrabcMotionDirection max_dir = pcs->parent_pcs_ptr->ibc_mode > 1 ? IBC_MOTION_LEFT : IBC_MOTION_DIRECTIONS;

    for (enum IntrabcMotionDirection dir = IBC_MOTION_ABOVE;
        dir < max_dir; ++dir) {
#else
    for (enum IntrabcMotionDirection dir = IBC_MOTION_ABOVE;
        dir < IBC_MOTION_DIRECTIONS; ++dir) {
#endif

        const MvLimits tmp_mv_limits = x->mv_limits;

        switch (dir) {
        case IBC_MOTION_ABOVE:
            x->mv_limits.col_min = (tile->mi_col_start - mi_col) * MI_SIZE;
            x->mv_limits.col_max = (tile->mi_col_end - mi_col) * MI_SIZE - w;
            x->mv_limits.row_min = (tile->mi_row_start - mi_row) * MI_SIZE;
            x->mv_limits.row_max =
                (sb_row * scs->mib_size - mi_row) * MI_SIZE - h;
            break;
        case IBC_MOTION_LEFT:
            x->mv_limits.col_min = (tile->mi_col_start - mi_col) * MI_SIZE;
            x->mv_limits.col_max =
                (sb_col * scs->mib_size - mi_col) * MI_SIZE - w;
            // TODO: Minimize the overlap between above and
            // left areas.
            x->mv_limits.row_min = (tile->mi_row_start - mi_row) * MI_SIZE;
            int bottom_coded_mi_edge =
                AOMMIN((sb_row + 1) * scs->mib_size, tile->mi_row_end);
            x->mv_limits.row_max = (bottom_coded_mi_edge - mi_row) * MI_SIZE - h;
            break;
        default: assert(0);
        }
        assert_release(x->mv_limits.col_min >= tmp_mv_limits.col_min);
        assert_release(x->mv_limits.col_max <= tmp_mv_limits.col_max);
        assert_release(x->mv_limits.row_min >= tmp_mv_limits.row_min);
        assert_release(x->mv_limits.row_max <= tmp_mv_limits.row_max);

        av1_set_mv_search_range(&x->mv_limits, &dv_ref.as_mv);

        if (x->mv_limits.col_max < x->mv_limits.col_min ||
            x->mv_limits.row_max < x->mv_limits.row_min) {
            x->mv_limits = tmp_mv_limits;
            continue;
        }

        int step_param = 0;
        MV mvp_full = dv_ref.as_mv;
        mvp_full.col >>= 3;
        mvp_full.row >>= 3;
        const int sadpb = x->sadperbit16;
        x->best_mv.as_int = 0;

#define INT_VAR_MAX  2147483647    // maximum (signed) int value

        const int bestsme = av1_full_pixel_search(
            pcs, x, bsize, &mvp_full, step_param, 1, 0,
            sadpb, NULL, &dv_ref.as_mv, INT_VAR_MAX, 1,
            (MI_SIZE * mi_col), (MI_SIZE * mi_row), 1);

        x->mv_limits = tmp_mv_limits;
        if (bestsme == INT_VAR_MAX) continue;
        mvp_full = x->best_mv.as_mv;

        const MV dv = { .row = mvp_full.row * 8,.col = mvp_full.col * 8 };
        if (mv_check_bounds(&x->mv_limits, &dv)) continue;
        if (!av1_is_dv_valid(dv, xd, mi_row, mi_col, bsize,
            scs->mib_size_log2))
            continue;

        // DV should not have sub-pel.
        assert_release((dv.col & 7) == 0);
        assert_release((dv.row & 7) == 0);

        //store output
        dv_cand[*num_dv_cand] = dv;
        (*num_dv_cand)++;

    }

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            free(x->hash_value_buffer[i][j]);

}

void  inject_intra_bc_candidates(
    PictureControlSet_t            *picture_control_set_ptr,
    ModeDecisionContext_t          *context_ptr,
    const SequenceControlSet_t     *sequence_control_set_ptr,
    LargestCodingUnit_t            *sb_ptr,
    CodingUnit_t                   *cu_ptr,
    uint32_t                       *cand_cnt)
{
    MV dv_cand[2];
    uint8_t num_dv_cand = 0;

    //perform dv-pred + search up to 2 dv(s)
    intra_bc_search(
        picture_control_set_ptr,
        context_ptr,
        sequence_control_set_ptr,
        sb_ptr,
        cu_ptr,
        dv_cand,
        &num_dv_cand);

    ModeDecisionCandidate_t    *candidateArray = context_ptr->fast_candidate_array;
    uint32_t dv_i;

    for (dv_i = 0; dv_i < num_dv_cand; dv_i++)
    {
        candidateArray[*cand_cnt].type = INTRA_MODE;
        candidateArray[*cand_cnt].intra_luma_mode = DC_PRED;
#if TWO_FAST_LOOP 
        candidateArray[*cand_cnt].enable_two_fast_loops = 0;
#else
        candidateArray[*cand_cnt].distortion_ready = 0;
#endif
        candidateArray[*cand_cnt].use_intrabc = 1;
        candidateArray[*cand_cnt].is_directional_mode_flag = 0;
        candidateArray[*cand_cnt].use_angle_delta = 0;
        candidateArray[*cand_cnt].angle_delta[PLANE_TYPE_Y] = 0;
        candidateArray[*cand_cnt].intra_chroma_mode = UV_DC_PRED;
        candidateArray[*cand_cnt].cfl_alpha_signs = 0;
        candidateArray[*cand_cnt].cfl_alpha_idx = 0;
        candidateArray[*cand_cnt].is_directional_chroma_mode_flag = 0;
        candidateArray[*cand_cnt].angle_delta[PLANE_TYPE_UV] = 0;
        candidateArray[*cand_cnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
        candidateArray[*cand_cnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;
        candidateArray[*cand_cnt].ref_frame_type = INTRA_FRAME;
        candidateArray[*cand_cnt].pred_mode = DC_PRED;
        candidateArray[*cand_cnt].motion_mode = SIMPLE_TRANSLATION;
        //inter ralated
        candidateArray[*cand_cnt].is_compound = 0;
        candidateArray[*cand_cnt].merge_flag = EB_FALSE;
        candidateArray[*cand_cnt].prediction_direction[0] = UNI_PRED_LIST_0;
        candidateArray[*cand_cnt].is_new_mv = 0;
        candidateArray[*cand_cnt].is_zero_mv = 0;
        candidateArray[*cand_cnt].motionVector_x_L0 = dv_cand[dv_i].col;
        candidateArray[*cand_cnt].motionVector_y_L0 = dv_cand[dv_i].row;
        candidateArray[*cand_cnt].motion_vector_pred_x[REF_LIST_0] = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[INTRA_FRAME][0].this_mv.as_mv.col;
        candidateArray[*cand_cnt].motion_vector_pred_y[REF_LIST_0] = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[INTRA_FRAME][0].this_mv.as_mv.row;
        candidateArray[*cand_cnt].drl_index = 0;
        candidateArray[*cand_cnt].ref_mv_index = 0;
        candidateArray[*cand_cnt].pred_mv_weight = 0;
        candidateArray[*cand_cnt].interp_filters = av1_broadcast_interp_filter(BILINEAR);
        ++(*cand_cnt);

    }

}
#endif
// END of Function Declarations
void  inject_intra_candidates(
    PictureControlSet_t            *picture_control_set_ptr,
    ModeDecisionContext_t          *context_ptr,
    const SequenceControlSet_t     *sequence_control_set_ptr,
    LargestCodingUnit_t            *sb_ptr,
#if M8_SKIP_BLK
    uint32_t                       *candidateTotalCnt){
#else
    uint32_t                       *candidateTotalCnt,
    const uint32_t                  leaf_index){

    (void)leaf_index;
#endif
    (void)sequence_control_set_ptr;
    (void)sb_ptr;
#if ENABLE_PAETH
    uint8_t                     is16bit = (sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
#endif
    uint8_t                     intra_mode_start = DC_PRED;
#if ENABLE_PAETH
    uint8_t                     intra_mode_end   = is16bit ? SMOOTH_H_PRED : PAETH_PRED;
#else    
    uint8_t                     intra_mode_end   = SMOOTH_H_PRED;
#endif
    uint8_t                     openLoopIntraCandidate;
    uint32_t                    canTotalCnt = 0;
    uint8_t                     angleDeltaCounter = 0;
    EbBool                      use_angle_delta = (context_ptr->blk_geom->bsize >= BLOCK_8X8);
    uint8_t                     angleDeltaCandidateCount = use_angle_delta ? 7 : 1;
    ModeDecisionCandidate_t    *candidateArray = context_ptr->fast_candidate_array;

    EbBool                      disable_cfl_flag = (context_ptr->blk_geom->sq_size > 32 || 
                                                    context_ptr->blk_geom->bwidth == 4  ||   
                                                    context_ptr->blk_geom->bheight == 4)    ? EB_TRUE : EB_FALSE;

    uint8_t                     disable_z2_prediction;
    uint8_t                     disable_angle_refinement;
    uint8_t                     disable_angle_prediction;

    if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode == 3){
        disable_z2_prediction       = 0;
        disable_angle_refinement    = 0;
        disable_angle_prediction    = 1;
    } else if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode == 2) {
        disable_z2_prediction       = 0;
        disable_angle_refinement    = 0 ;
        disable_angle_prediction    = (context_ptr->blk_geom->sq_size > 16 ||
                                       context_ptr->blk_geom->bwidth == 4 ||
                                       context_ptr->blk_geom->bheight == 4) ? 1 : 0;
    } else if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode == 1) {
        disable_z2_prediction       = (context_ptr->blk_geom->sq_size > 16 ||
                                       context_ptr->blk_geom->bwidth == 4 ||
                                       context_ptr->blk_geom->bheight == 4) ? 1 : 0;
        disable_angle_refinement    = (context_ptr->blk_geom->sq_size > 16 ||
                                       context_ptr->blk_geom->bwidth == 4 ||
                                       context_ptr->blk_geom->bheight == 4) ? 1 : 0;
        disable_angle_prediction    = 0;
    } else {
        disable_z2_prediction       = 0;
        disable_angle_refinement    = 0;
        disable_angle_prediction    = 0;

    }
#if MR_MODE
    disable_z2_prediction       = 0;
    disable_angle_refinement    = 0;
    disable_angle_prediction    = 0;
#endif
    angleDeltaCandidateCount = disable_angle_refinement ? 1: angleDeltaCandidateCount;
#if !NSQ_OPTIMASATION
    uint8_t sq_index = LOG2F(context_ptr->blk_geom->sq_size) - 2;
    if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level == NSQ_INTER_SEARCH_BASE_ON_SQ_INTRAMODE) {
        disable_z2_prediction = context_ptr->blk_geom->shape == PART_N ? disable_z2_prediction :
            context_ptr->parent_sq_type[sq_index] == INTRA_MODE ? disable_z2_prediction : 0;
    }
#endif
#if TWO_FAST_LOOP 
    uint8_t enable_two_fast_loops = picture_control_set_ptr->parent_pcs_ptr->enable_two_fast_loops && (context_ptr->blk_geom->sq_size > 4 && context_ptr->blk_geom->shape == PART_N);          
#endif
#if !DIS_EDGE_FIL
    const int32_t disable_ang_uv = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) && context_ptr->blk_geom->has_uv ? 1 : 0;
#endif
    for (openLoopIntraCandidate = intra_mode_start; openLoopIntraCandidate <= intra_mode_end ; ++openLoopIntraCandidate) {

        if (av1_is_directional_mode((PredictionMode)openLoopIntraCandidate)) {

            if (!disable_angle_prediction) {
                for (angleDeltaCounter = 0; angleDeltaCounter < angleDeltaCandidateCount; ++angleDeltaCounter) {
                    int32_t angle_delta = angleDeltaCandidateCount == 1 ? 0 : angleDeltaCounter - (angleDeltaCandidateCount >> 1);
                    int32_t  p_angle = mode_to_angle_map[(PredictionMode)openLoopIntraCandidate] + angle_delta * ANGLE_STEP;
                    if (!disable_z2_prediction || (p_angle <= 90 || p_angle >= 180)) {
                        candidateArray[canTotalCnt].type = INTRA_MODE;
                        candidateArray[canTotalCnt].intra_luma_mode = openLoopIntraCandidate;
#if TWO_FAST_LOOP 
                        candidateArray[canTotalCnt].enable_two_fast_loops = enable_two_fast_loops;
#if OIS_BASED_INTRA
                        candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
                        candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
                        candidateArray[canTotalCnt].use_intrabc = 0;
#endif
                        candidateArray[canTotalCnt].is_directional_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)openLoopIntraCandidate);
                        candidateArray[canTotalCnt].use_angle_delta = use_angle_delta ? candidateArray[canTotalCnt].is_directional_mode_flag : 0;
                        candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y] = angle_delta;
#if CHROMA_BLIND 
                        candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ? 
                            intra_luma_to_chroma[openLoopIntraCandidate] : 
                            (context_ptr->chroma_level == CHROMA_MODE_0) ?
                                UV_CFL_PRED :
                                UV_DC_PRED;
#else
                        candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ? intra_luma_to_chroma[openLoopIntraCandidate] : UV_CFL_PRED;
#endif
#if !DIS_EDGE_FIL
                        candidateArray[canTotalCnt].intra_chroma_mode = disable_ang_uv && av1_is_directional_mode(candidateArray[canTotalCnt].intra_chroma_mode) ?
                            UV_DC_PRED : candidateArray[canTotalCnt].intra_chroma_mode;
#endif
                        candidateArray[canTotalCnt].cfl_alpha_signs = 0;
                        candidateArray[canTotalCnt].cfl_alpha_idx = 0;
                        candidateArray[canTotalCnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidateArray[canTotalCnt].intra_chroma_mode);
                        candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_UV] = 0;
                        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;

                        if (candidateArray[canTotalCnt].intra_chroma_mode == UV_CFL_PRED)
                            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;
                        else
                            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] =
                            av1_get_tx_type(
                                context_ptr->blk_geom->bsize,
                                0,
                                (PredictionMode)candidateArray[canTotalCnt].intra_luma_mode,
                                (UV_PredictionMode)candidateArray[canTotalCnt].intra_chroma_mode,
                                PLANE_TYPE_UV,
                                0,
                                0,
                                0,
                                context_ptr->blk_geom->txsize_uv[0],
                                picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
#if !INTRA_INTER_FAST_LOOP
                        // candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV]            = context_ptr->blk_geom->bwidth_uv >= 32 || context_ptr->blk_geom->bheight_uv >= 32  ? DCT_DCT : chroma_transform_type[candidateArray[canTotalCnt].intra_chroma_mode];
                        candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
                        candidateArray[canTotalCnt].ref_frame_type = INTRA_FRAME;
                        candidateArray[canTotalCnt].pred_mode = (PredictionMode)openLoopIntraCandidate;
                        candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                        ++canTotalCnt;
                    }
            }
        }
        }
        else {
            candidateArray[canTotalCnt].type = INTRA_MODE;
            candidateArray[canTotalCnt].intra_luma_mode = openLoopIntraCandidate;
#if TWO_FAST_LOOP 
            candidateArray[canTotalCnt].enable_two_fast_loops =  enable_two_fast_loops;
#if OIS_BASED_INTRA
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#else
            candidateArray[canTotalCnt].distortion_ready = 0;
#endif
#if ICOPY
            candidateArray[canTotalCnt].use_intrabc = 0;
#endif
            candidateArray[canTotalCnt].is_directional_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)openLoopIntraCandidate);
            candidateArray[canTotalCnt].use_angle_delta = candidateArray[canTotalCnt].is_directional_mode_flag;
            candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y] = 0;
#if TURN_OFF_CFL
            candidateArray[canTotalCnt].intra_chroma_mode = intra_luma_to_chroma[openLoopIntraCandidate];
#else
#if CHROMA_BLIND
            candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ? 
                intra_luma_to_chroma[openLoopIntraCandidate] : 
                (context_ptr->chroma_level == CHROMA_MODE_0) ?
                    UV_CFL_PRED :
                    UV_DC_PRED;

#else
            candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ? intra_luma_to_chroma[openLoopIntraCandidate] : UV_CFL_PRED;
#endif
#endif
#if !DIS_EDGE_FIL
            candidateArray[canTotalCnt].intra_chroma_mode = disable_ang_uv && av1_is_directional_mode(candidateArray[canTotalCnt].intra_chroma_mode) ?
                UV_DC_PRED : candidateArray[canTotalCnt].intra_chroma_mode;
#endif
            candidateArray[canTotalCnt].cfl_alpha_signs = 0;
            candidateArray[canTotalCnt].cfl_alpha_idx = 0;
            candidateArray[canTotalCnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidateArray[canTotalCnt].intra_chroma_mode);
            candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_UV] = 0;
            candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;

            if (candidateArray[canTotalCnt].intra_chroma_mode == UV_CFL_PRED)
                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;
            else
                candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] =
                av1_get_tx_type(
                    context_ptr->blk_geom->bsize,
                    0,
                    (PredictionMode)candidateArray[canTotalCnt].intra_luma_mode,
                    (UV_PredictionMode)candidateArray[canTotalCnt].intra_chroma_mode,
                    PLANE_TYPE_UV,
                    0,
                    0,
                    0,
                    context_ptr->blk_geom->txsize_uv[0],
                    picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
#if !INTRA_INTER_FAST_LOOP
            candidateArray[canTotalCnt].mpm_flag = EB_FALSE;
#endif
            candidateArray[canTotalCnt].ref_frame_type = INTRA_FRAME;
            candidateArray[canTotalCnt].pred_mode = (PredictionMode)openLoopIntraCandidate;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            ++canTotalCnt;
        }
    }


    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;


    return;
}

void ProductInitMdCandInjection(
    ModeDecisionContext_t          *context_ptr,
    uint32_t                         *candidateTotalCnt)

{

    *candidateTotalCnt = 0;
    context_ptr->generate_mvp = EB_FALSE;

    return;
}
/***************************************
* ProductGenerateMdCandidatesCu
*   Creates list of initial modes to
*   perform fast cost search on.
***************************************/
EbErrorType ProductGenerateMdCandidatesCu(
    LargestCodingUnit_t                 *sb_ptr,
    ModeDecisionContext_t             *context_ptr,
    SsMeContext_t                    *ss_mecontext,
#if !M8_SKIP_BLK
    const uint32_t                      leaf_index,
#endif
    const uint32_t                      lcuAddr,
#if !INTRA_INTER_FAST_LOOP
    uint32_t                           *bufferTotalCountPtr,
#endif
    uint32_t                           *candidateTotalCountPtr,
    EbPtr                              interPredContextPtr,
    PictureControlSet_t              *picture_control_set_ptr)
{

    (void)lcuAddr;
    (void)interPredContextPtr;
    const SequenceControlSet_t *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    const EB_SLICE slice_type = picture_control_set_ptr->slice_type;
    uint32_t       canTotalCnt;

#if REMOVED_DUPLICATE_INTER
    // Reset duplicates variables
    context_ptr->injected_mv_count_l0 = 0;
    context_ptr->injected_mv_count_l1 = 0;
    context_ptr->injected_mv_count_bipred = 0;
#endif

    ProductInitMdCandInjection(
        context_ptr,
        &canTotalCnt);

    uint8_t sq_index = LOG2F(context_ptr->blk_geom->sq_size) - 2;
    uint8_t inject_intra_candidate = 1;
    uint8_t inject_inter_candidate = 1;

    if (slice_type != I_SLICE) {
#if !NSQ_OPTIMASATION
        if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level == NSQ_SEARCH_BASE_ON_SQ_TYPE) {
            inject_intra_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
                context_ptr->parent_sq_type[sq_index] == INTRA_MODE ? inject_intra_candidate : 0;
            inject_inter_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
                context_ptr->parent_sq_type[sq_index] == INTER_MODE ? inject_inter_candidate : 0;
        }
#endif
#if NSQ_OPTIMASATION
        if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level >= NSQ_SEARCH_LEVEL1 &&
            picture_control_set_ptr->parent_pcs_ptr->nsq_search_level < NSQ_SEARCH_FULL) {
            inject_intra_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
                context_ptr->parent_sq_has_coeff[sq_index] != 0 ? inject_intra_candidate : 0;
        }
#else
        if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level == NSQ_SEARCH_BASE_ON_SQ_COEFF) {
            inject_intra_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
                context_ptr->parent_sq_has_coeff[sq_index] != 0 ? inject_intra_candidate : 0;
        }
#endif
}
    //----------------------
    // Intra
    if (context_ptr->blk_geom->sq_size < 128) {
#if OIS_BASED_INTRA
        if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode >= 4 && context_ptr->blk_geom->sq_size > 4 && context_ptr->blk_geom->shape == PART_N)
            inject_intra_candidates_ois(
                picture_control_set_ptr,
                context_ptr,
                sb_ptr,
                &canTotalCnt);
        else 
#endif
            if (inject_intra_candidate)
            inject_intra_candidates(
                picture_control_set_ptr,
                context_ptr,
                sequence_control_set_ptr,
                sb_ptr,
#if M8_SKIP_BLK
                &canTotalCnt);
#else
                &canTotalCnt,
                leaf_index);
#endif
    }

#if 0//INTRA_INTER_FAST_LOOP
    // Track the total number of fast intra candidates
    context_ptr->fast_candidate_intra_count = canTotalCnt;
#endif

#if ICOPY
    if (picture_control_set_ptr->parent_pcs_ptr->allow_intrabc)
        inject_intra_bc_candidates(
            picture_control_set_ptr,
            context_ptr,
            sequence_control_set_ptr,
            sb_ptr,
            context_ptr->cu_ptr,
            &canTotalCnt
        );
#endif

#if INTRA_INTER_FAST_LOOP
    // Track the total number of fast intra candidates
    context_ptr->fast_candidate_intra_count = canTotalCnt;
#endif

    if (slice_type != I_SLICE) {
        if (inject_inter_candidate)
            inject_inter_candidates(
                picture_control_set_ptr,
                context_ptr,
                ss_mecontext,
                sequence_control_set_ptr,
                sb_ptr,
#if M8_SKIP_BLK
                &canTotalCnt);
#else
                &canTotalCnt,
                leaf_index);
#endif
    }
#if !INTRA_INTER_FAST_LOOP
    // Set BufferTotalCount: determines the number of candidates to fully reconstruct
    *bufferTotalCountPtr = context_ptr->full_recon_search_count;
#endif
    *candidateTotalCountPtr = canTotalCnt;
#if !INTRA_INTER_FAST_LOOP
    // Make sure buffer_total_count is not larger than the number of fast modes
    *bufferTotalCountPtr = MIN(*candidateTotalCountPtr, *bufferTotalCountPtr);
#endif
    return EB_ErrorNone;
}

/***************************************
* Full Mode Decision
***************************************/
uint8_t product_full_mode_decision(
    struct ModeDecisionContext_s   *context_ptr,
    CodingUnit_t                   *cu_ptr,
    uint8_t                           bwidth,
    uint8_t                           bheight,
    ModeDecisionCandidateBuffer_t **buffer_ptr_array,
    uint32_t                          candidate_total_count,
    uint8_t                          *best_candidate_index_array,
    uint32_t                           *best_intra_mode)
{
    UNUSED(bwidth);
    UNUSED(bheight);

    uint8_t                   candidateIndex;
    uint64_t                  lowestCost = 0xFFFFFFFFFFFFFFFFull;
    uint64_t                  lowestIntraCost = 0xFFFFFFFFFFFFFFFFull;
    uint8_t                   lowestCostIndex = 0;
    PredictionUnit_t       *pu_ptr;
    uint32_t                   i;
    ModeDecisionCandidate_t       *candidate_ptr;

    lowestCostIndex = best_candidate_index_array[0];

    // Find the candidate with the lowest cost
    for (i = 0; i < candidate_total_count; ++i) {
        candidateIndex = best_candidate_index_array[i];

        // Compute fullCostBis
        if ((*(buffer_ptr_array[candidateIndex]->full_cost_ptr) < lowestIntraCost) && buffer_ptr_array[candidateIndex]->candidate_ptr->type == INTRA_MODE) {
            *best_intra_mode = buffer_ptr_array[candidateIndex]->candidate_ptr->pred_mode;
            lowestIntraCost = *(buffer_ptr_array[candidateIndex]->full_cost_ptr);
        }

        if (*(buffer_ptr_array[candidateIndex]->full_cost_ptr) < lowestCost) {
            lowestCostIndex = candidateIndex;
            lowestCost = *(buffer_ptr_array[candidateIndex]->full_cost_ptr);
        }
    }

    candidate_ptr = buffer_ptr_array[lowestCostIndex]->candidate_ptr;

    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost = *(buffer_ptr_array[lowestCostIndex]->full_cost_ptr);
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost = (context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost - buffer_ptr_array[lowestCostIndex]->candidate_ptr->chroma_distortion) + buffer_ptr_array[lowestCostIndex]->candidate_ptr->chroma_distortion_inter_depth;

    if (candidate_ptr->type == INTRA_MODE)
        context_ptr->md_local_cu_unit[cu_ptr->mds_idx].cost_luma = buffer_ptr_array[lowestCostIndex]->full_cost_luma;

    context_ptr->md_ep_pipe_sb[cu_ptr->mds_idx].merge_cost = *buffer_ptr_array[lowestCostIndex]->full_cost_merge_ptr;
    context_ptr->md_ep_pipe_sb[cu_ptr->mds_idx].skip_cost = *buffer_ptr_array[lowestCostIndex]->full_cost_skip_ptr;


    if (candidate_ptr->type == INTER_MODE && candidate_ptr->merge_flag == EB_TRUE) {
        context_ptr->md_ep_pipe_sb[cu_ptr->mds_idx].chroma_distortion = buffer_ptr_array[lowestCostIndex]->candidate_ptr->chroma_distortion;
    }

    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].full_distortion = buffer_ptr_array[lowestCostIndex]->candidate_ptr->full_distortion;
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].chroma_distortion = (uint32_t)buffer_ptr_array[lowestCostIndex]->candidate_ptr->chroma_distortion;
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].chroma_distortion_inter_depth = (uint32_t)buffer_ptr_array[lowestCostIndex]->candidate_ptr->chroma_distortion_inter_depth;

    cu_ptr->prediction_mode_flag = candidate_ptr->type;
    cu_ptr->skip_flag = candidate_ptr->skip_flag; // note, the skip flag is re-checked in the ENCDEC process
    cu_ptr->block_has_coeff = ((candidate_ptr->block_has_coeff) > 0) ? EB_TRUE : EB_FALSE;
    cu_ptr->quantized_dc[0] = buffer_ptr_array[lowestCostIndex]->candidate_ptr->quantized_dc[0];
    cu_ptr->quantized_dc[1] = buffer_ptr_array[lowestCostIndex]->candidate_ptr->quantized_dc[1];
    cu_ptr->quantized_dc[2] = buffer_ptr_array[lowestCostIndex]->candidate_ptr->quantized_dc[2];

    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].count_non_zero_coeffs = candidate_ptr->count_non_zero_coeffs;

#if ICOPY    
    cu_ptr->av1xd->use_intrabc = candidate_ptr->use_intrabc;
#endif

    // Set the PU level variables
    cu_ptr->interp_filters = candidate_ptr->interp_filters;
    {
        pu_ptr = cu_ptr->prediction_unit_array;
        // Intra Prediction
        pu_ptr->intra_luma_mode = 0x1F;
        if (cu_ptr->prediction_mode_flag == INTRA_MODE)
        {
            pu_ptr->intra_luma_mode = candidate_ptr->intra_luma_mode;

            pu_ptr->is_directional_mode_flag = candidate_ptr->is_directional_mode_flag;
            pu_ptr->use_angle_delta = candidate_ptr->use_angle_delta;
            pu_ptr->angle_delta[PLANE_TYPE_Y] = candidate_ptr->angle_delta[PLANE_TYPE_Y];


            pu_ptr->cfl_alpha_idx = candidate_ptr->cfl_alpha_idx;
            pu_ptr->cfl_alpha_signs = candidate_ptr->cfl_alpha_signs;

            pu_ptr->intra_chroma_mode = candidate_ptr->intra_chroma_mode;
            pu_ptr->is_directional_chroma_mode_flag = candidate_ptr->is_directional_chroma_mode_flag;
            pu_ptr->angle_delta[PLANE_TYPE_UV] = candidate_ptr->angle_delta[PLANE_TYPE_UV];
        }

        // Inter Prediction
        pu_ptr->inter_pred_direction_index = candidate_ptr->prediction_direction[0];
        pu_ptr->merge_flag = candidate_ptr->merge_flag;
#if ICOPY
        if (cu_ptr->prediction_mode_flag != INTER_MODE && cu_ptr->av1xd->use_intrabc == 0)
#else
        if (cu_ptr->prediction_mode_flag != INTER_MODE)
#endif
        {
            pu_ptr->inter_pred_direction_index = 0x03;
            pu_ptr->merge_flag = EB_FALSE;
        }
#if !INTRA_INTER_FAST_LOOP
        pu_ptr->merge_index = candidate_ptr->merge_index;
#endif
        pu_ptr->mv[REF_LIST_0].x = 0;
        pu_ptr->mv[REF_LIST_0].y = 0;

        pu_ptr->mv[REF_LIST_1].x = 0;
        pu_ptr->mv[REF_LIST_1].y = 0;

        cu_ptr->pred_mode = candidate_ptr->pred_mode;
        cu_ptr->drl_index = candidate_ptr->drl_index;

        pu_ptr->inter_mode = candidate_ptr->inter_mode;
        pu_ptr->is_compound = candidate_ptr->is_compound;
        pu_ptr->pred_mv_weight = candidate_ptr->pred_mv_weight;
        pu_ptr->ref_frame_type = candidate_ptr->ref_frame_type;
        pu_ptr->ref_mv_index = candidate_ptr->ref_mv_index;
#if !INTRA_INTER_FAST_LOOP
        pu_ptr->is_skip_mode_flag = candidate_ptr->is_skip_mode_flag;
#endif
        pu_ptr->is_new_mv = candidate_ptr->is_new_mv;
        pu_ptr->is_zero_mv = candidate_ptr->is_zero_mv;

        if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0)
        {
            //EB_MEMCPY(&pu_ptr->mv[REF_LIST_0].x,&candidate_ptr->MVsL0,4);
            pu_ptr->mv[REF_LIST_0].x = candidate_ptr->motionVector_x_L0;
            pu_ptr->mv[REF_LIST_0].y = candidate_ptr->motionVector_y_L0;
        }

        if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1)
        {
            //EB_MEMCPY(&pu_ptr->mv[REF_LIST_1].x,&candidate_ptr->MVsL1,4);
            pu_ptr->mv[REF_LIST_1].x = candidate_ptr->motionVector_x_L1;
            pu_ptr->mv[REF_LIST_1].y = candidate_ptr->motionVector_y_L1;
        }

        if (pu_ptr->inter_pred_direction_index == BI_PRED)
        {
            //EB_MEMCPY(&pu_ptr->mv[REF_LIST_0].x,&candidate_ptr->MVs,8);
            pu_ptr->mv[REF_LIST_0].x = candidate_ptr->motionVector_x_L0;
            pu_ptr->mv[REF_LIST_0].y = candidate_ptr->motionVector_y_L0;
            pu_ptr->mv[REF_LIST_1].x = candidate_ptr->motionVector_x_L1;
            pu_ptr->mv[REF_LIST_1].y = candidate_ptr->motionVector_y_L1;
        }

        // The MV prediction indicies are recalcated by the EncDec.
        pu_ptr->mvd[REF_LIST_0].predIdx = 0;
        pu_ptr->mvd[REF_LIST_1].predIdx = 0;

        pu_ptr->overlappable_neighbors[0] = context_ptr->cu_ptr->prediction_unit_array[0].overlappable_neighbors[0];
        pu_ptr->overlappable_neighbors[1] = context_ptr->cu_ptr->prediction_unit_array[0].overlappable_neighbors[1];
        pu_ptr->motion_mode = candidate_ptr->motion_mode;
        pu_ptr->num_proj_ref = candidate_ptr->num_proj_ref;
        if (pu_ptr->motion_mode == WARPED_CAUSAL)
            EB_MEMCPY(&pu_ptr->wm_params, &candidate_ptr->wm_params, sizeof(EbWarpedMotionParams));
    }

    TransformUnit_t        *txb_ptr;
    uint32_t                  txb_itr;
    uint32_t                  tu_index;
    uint32_t                  tuTotalCount;
    uint32_t  cu_size_log2 = context_ptr->cu_size_log2;

    {
        tuTotalCount = context_ptr->blk_geom->txb_count;
        tu_index = 0;
        txb_itr = 0;
    }

#if NO_ENCDEC
    int32_t txb_1d_offset = 0, txb_1d_offset_uv = 0;

    cu_ptr->block_has_coeff = 0;
#endif

    //cu_ptr->forceSmallTu = candidate_ptr->forceSmallTu;

    // Set TU
    do {

        txb_ptr = &cu_ptr->transform_unit_array[tu_index];

        txb_ptr->split_flag = EB_FALSE;
        txb_ptr->y_has_coeff = (EbBool)(((candidate_ptr->y_has_coeff)  & (1 << tu_index)) > 0);
        txb_ptr->u_has_coeff = (EbBool)(((candidate_ptr->u_has_coeff) & (1 << (tu_index))) > 0);
        txb_ptr->v_has_coeff = (EbBool)(((candidate_ptr->v_has_coeff) & (1 << (tu_index))) > 0);
        txb_ptr->transform_type[PLANE_TYPE_Y] = candidate_ptr->transform_type[PLANE_TYPE_Y];
        txb_ptr->transform_type[PLANE_TYPE_UV] = candidate_ptr->transform_type[PLANE_TYPE_UV];


#if NO_ENCDEC

        if (context_ptr->blk_geom->has_uv) {
            cu_ptr->block_has_coeff |= txb_ptr->y_has_coeff;
            cu_ptr->block_has_coeff |= txb_ptr->u_has_coeff;
            cu_ptr->block_has_coeff |= txb_ptr->v_has_coeff;
        }
        else {
            cu_ptr->block_has_coeff |= txb_ptr->y_has_coeff;
        }


        cu_ptr->cand_buff_index = lowestCostIndex;

        cu_ptr->skip_flag = 0;   //SKIP is turned OFF for this case!!
        txb_ptr->nz_coef_count[0] = candidate_ptr->eob[0][tu_index];
        txb_ptr->nz_coef_count[1] = candidate_ptr->eob[1][tu_index];
        txb_ptr->nz_coef_count[2] = candidate_ptr->eob[2][tu_index];

        if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0) {
            cu_ptr->predmv[0].as_mv.col = candidate_ptr->motion_vector_pred_x[REF_LIST_0];
            cu_ptr->predmv[0].as_mv.row = candidate_ptr->motion_vector_pred_y[REF_LIST_0];
        }
        else if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1) {
            cu_ptr->predmv[0].as_mv.col = candidate_ptr->motion_vector_pred_x[REF_LIST_1];
            cu_ptr->predmv[0].as_mv.row = candidate_ptr->motion_vector_pred_y[REF_LIST_1];
        }
        else if (pu_ptr->inter_pred_direction_index == BI_PRED) {
            cu_ptr->predmv[0].as_mv.col = candidate_ptr->motion_vector_pred_x[REF_LIST_0];
            cu_ptr->predmv[0].as_mv.row = candidate_ptr->motion_vector_pred_y[REF_LIST_0];
            cu_ptr->predmv[1].as_mv.col = candidate_ptr->motion_vector_pred_x[REF_LIST_1];
            cu_ptr->predmv[1].as_mv.row = candidate_ptr->motion_vector_pred_y[REF_LIST_1];
        }
#endif
#if NO_ENCDEC
        //copy coeff
        {

            uint32_t  bwidth = context_ptr->blk_geom->tx_width[txb_itr] < 64 ? context_ptr->blk_geom->tx_width[txb_itr] : 32;
            uint32_t  bheight = context_ptr->blk_geom->tx_height[txb_itr] < 64 ? context_ptr->blk_geom->tx_height[txb_itr] : 32;

            int32_t* src_ptr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residualQuantCoeffPtr->buffer_y)[txb_1d_offset]);
            int32_t* dst_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->buffer_y)[txb_1d_offset]);

            uint32_t j;

            for (j = 0; j < bheight; j++)
            {
                memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
            }

            if (context_ptr->blk_geom->has_uv)
            {
                // Cb
                bwidth = context_ptr->blk_geom->tx_width_uv[txb_itr];
                bheight = context_ptr->blk_geom->tx_height_uv[txb_itr];

                src_ptr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residualQuantCoeffPtr->bufferCb)[txb_1d_offset_uv]);
                dst_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->bufferCb)[txb_1d_offset_uv]);

                for (j = 0; j < bheight; j++)
                {
                    memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                }

                src_ptr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residualQuantCoeffPtr->bufferCr)[txb_1d_offset_uv]);
                dst_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->bufferCr)[txb_1d_offset_uv]);

                for (j = 0; j < bheight; j++)
                {
                    memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                }
            }

            txb_1d_offset += context_ptr->blk_geom->tx_width[txb_itr] * context_ptr->blk_geom->tx_height[txb_itr];
            if (context_ptr->blk_geom->has_uv)
                txb_1d_offset_uv += context_ptr->blk_geom->tx_width_uv[txb_itr] * context_ptr->blk_geom->tx_height_uv[txb_itr];

        }

#endif

        ++tu_index;
        ++txb_itr;

    } while (txb_itr < tuTotalCount);
    UNUSED(cu_size_log2);
    return lowestCostIndex;
}

