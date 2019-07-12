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

#include "av1me.h"
#include "hash.h"

#define  INCRMENT_CAND_TOTAL_COUNT(cnt) cnt++; if(cnt>=MODE_DECISION_CANDIDATE_MAX_COUNT) printf(" ERROR: reaching limit for MODE_DECISION_CANDIDATE_MAX_COUNT %i\n",cnt);
int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
/********************************************
* Constants
********************************************/
//static uint32_t  AntiContouringIntraMode[11] = { EB_INTRA_PLANAR, EB_INTRA_DC, EB_INTRA_HORIZONTAL, EB_INTRA_VERTICAL,
//EB_INTRA_MODE_2, EB_INTRA_MODE_6, EB_INTRA_MODE_14, EB_INTRA_MODE_18, EB_INTRA_MODE_22, EB_INTRA_MODE_30, EB_INTRA_MODE_34 };

const uint32_t parentIndex[85] = { 0, 0, 0, 2, 2, 2, 2, 0, 7, 7, 7, 7, 0, 12, 12, 12, 12, 0, 17, 17, 17, 17, 0, 0,
23, 23, 23, 23, 0, 28, 28, 28, 28, 0, 33, 33, 33, 33, 0, 38, 38, 38, 38, 0, 0,
44, 44, 44, 44, 0, 49, 49, 49, 49, 0, 54, 54, 54, 54, 0, 59, 59, 59, 59, 0, 0,
65, 65, 65, 65, 0, 70, 70, 70, 70, 0, 75, 75, 75, 75, 0, 80, 80, 80, 80 };
/*
  NORMAL ORDER
  |-------------------------------------------------------------|
  | ref_idx          0            1           2            3       |
  | List0            LAST        LAST2        LAST3        GOLD    |
  | List1            BWD            ALT2        ALT                 |
  |-------------------------------------------------------------|
*/
#define INVALID_REF 0xF

uint8_t get_list_idx(uint8_t ref_type) {
    if (ref_type == LAST_FRAME || ref_type == LAST2_FRAME || ref_type == LAST3_FRAME || ref_type == GOLDEN_FRAME)
        return 0;
    else if (ref_type == BWDREF_FRAME || ref_type == ALTREF_FRAME || ref_type == ALTREF2_FRAME)
        return 1;
    else
        return (INVALID_REF);
};

uint8_t get_ref_frame_idx(uint8_t ref_type) {
    if (ref_type == LAST_FRAME || ref_type == BWDREF_FRAME)
        return 0;
    else if (ref_type == LAST2_FRAME || ref_type == ALTREF2_FRAME)
        return 1;
    else if (ref_type == LAST3_FRAME || ref_type == ALTREF_FRAME)
        return 2;
    else if (ref_type == GOLDEN_FRAME)
        return 3;
    else
        return (INVALID_REF);
};
MvReferenceFrame svt_get_ref_frame_type(uint8_t list, uint8_t ref_idx) {
    switch (list) {
    case 0:
        return (ref_idx == 0 ? LAST_FRAME : ref_idx == 1 ? LAST2_FRAME : ref_idx == 2 ? LAST3_FRAME : ref_idx == 3 ? GOLDEN_FRAME : INVALID_REF);
    case 1:
        return (ref_idx == 0 ? BWDREF_FRAME : ref_idx == 1 ? ALTREF2_FRAME : ref_idx == 2 ? ALTREF_FRAME : INVALID_REF);
    default:
        return (INVALID_REF);
    }
};
extern uint32_t stage1ModesArray[];

uint8_t GetMaxDrlIndex(uint8_t  refmvCnt, PredictionMode   mode);
int32_t av1_mv_bit_cost(const MV *mv, const MV *ref, const int32_t *mvjcost,
    int32_t *mvcost[2], int32_t weight);
#define MV_COST_WEIGHT 108

void ChooseBestAv1MvPred(
    ModeDecisionContext            *context_ptr,
    struct MdRateEstimationContext      *md_rate_estimation_ptr,
    CodingUnit      *cu_ptr,
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

static void mode_decision_candidate_buffer_dctor(EbPtr p)
{
    ModeDecisionCandidateBuffer *obj = (ModeDecisionCandidateBuffer*)p;
    EB_DELETE(obj->prediction_ptr);
    EB_DELETE(obj->prediction_ptr_temp);
    EB_DELETE(obj->cfl_temp_prediction_ptr);
    EB_DELETE(obj->residual_ptr);
    EB_DELETE(obj->residual_quant_coeff_ptr);
    EB_DELETE(obj->recon_coeff_ptr);
    EB_DELETE(obj->recon_ptr);
}
/***************************************
* Mode Decision Candidate Ctor
***************************************/
EbErrorType mode_decision_candidate_buffer_ctor(
    ModeDecisionCandidateBuffer    *buffer_ptr,
    uint64_t                       *fast_cost_ptr,
    uint64_t                       *full_cost_ptr,
    uint64_t                       *full_cost_skip_ptr,
    uint64_t                       *full_cost_merge_ptr)
{
    EbPictureBufferDescInitData pictureBufferDescInitData;
    EbPictureBufferDescInitData doubleWidthPictureBufferDescInitData;

    EbPictureBufferDescInitData ThirtyTwoWidthPictureBufferDescInitData;


    buffer_ptr->dctor = mode_decision_candidate_buffer_dctor;

    // Init Picture Data
    pictureBufferDescInitData.max_width = MAX_SB_SIZE;
    pictureBufferDescInitData.max_height = MAX_SB_SIZE;
    pictureBufferDescInitData.bit_depth = EB_8BIT;
    pictureBufferDescInitData.color_format = EB_YUV420;
    pictureBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    pictureBufferDescInitData.left_padding = 0;
    pictureBufferDescInitData.right_padding = 0;
    pictureBufferDescInitData.top_padding = 0;
    pictureBufferDescInitData.bot_padding = 0;
    pictureBufferDescInitData.split_mode = EB_FALSE;
    doubleWidthPictureBufferDescInitData.max_width = MAX_SB_SIZE;
    doubleWidthPictureBufferDescInitData.max_height = MAX_SB_SIZE;
    doubleWidthPictureBufferDescInitData.bit_depth = EB_16BIT;
    doubleWidthPictureBufferDescInitData.color_format = EB_YUV420;
    doubleWidthPictureBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    doubleWidthPictureBufferDescInitData.left_padding = 0;
    doubleWidthPictureBufferDescInitData.right_padding = 0;
    doubleWidthPictureBufferDescInitData.top_padding = 0;
    doubleWidthPictureBufferDescInitData.bot_padding = 0;
    doubleWidthPictureBufferDescInitData.split_mode = EB_FALSE;

    ThirtyTwoWidthPictureBufferDescInitData.max_width = MAX_SB_SIZE;
    ThirtyTwoWidthPictureBufferDescInitData.max_height = MAX_SB_SIZE;
    ThirtyTwoWidthPictureBufferDescInitData.bit_depth = EB_32BIT;
    ThirtyTwoWidthPictureBufferDescInitData.color_format = EB_YUV420;
    ThirtyTwoWidthPictureBufferDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    ThirtyTwoWidthPictureBufferDescInitData.left_padding = 0;
    ThirtyTwoWidthPictureBufferDescInitData.right_padding = 0;
    ThirtyTwoWidthPictureBufferDescInitData.top_padding = 0;
    ThirtyTwoWidthPictureBufferDescInitData.bot_padding = 0;
    ThirtyTwoWidthPictureBufferDescInitData.split_mode = EB_FALSE;

    // Candidate Ptr
    buffer_ptr->candidate_ptr = (ModeDecisionCandidate*)EB_NULL;

    // Video Buffers
    EB_NEW(
        buffer_ptr->prediction_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&pictureBufferDescInitData);

    // Video Buffers
    EB_NEW(
        buffer_ptr->prediction_ptr_temp,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&pictureBufferDescInitData);

    EB_NEW(
        buffer_ptr->cfl_temp_prediction_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&pictureBufferDescInitData);

    EB_NEW(
        buffer_ptr->residual_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&doubleWidthPictureBufferDescInitData);

    EB_NEW(
        buffer_ptr->residual_quant_coeff_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&ThirtyTwoWidthPictureBufferDescInitData);

    EB_NEW(
        buffer_ptr->recon_coeff_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&ThirtyTwoWidthPictureBufferDescInitData);

    EB_NEW(
        buffer_ptr->recon_ptr,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&pictureBufferDescInitData);

    //Distortion
    buffer_ptr->residual_luma_sad = 0;

    buffer_ptr->full_lambda_rate = 0;

    // Costs
    buffer_ptr->fast_cost_ptr = fast_cost_ptr;
    buffer_ptr->full_cost_ptr = full_cost_ptr;
    buffer_ptr->full_cost_skip_ptr = full_cost_skip_ptr;
    buffer_ptr->full_cost_merge_ptr = full_cost_merge_ptr;
    return EB_ErrorNone;
}

// Function Declarations
void RoundMv(
    ModeDecisionCandidate    *candidateArray,
    uint32_t                   canTotalCnt)
{
    candidateArray[canTotalCnt].motion_vector_xl0 = (candidateArray[canTotalCnt].motion_vector_xl0 + 2)&~0x03;
    candidateArray[canTotalCnt].motion_vector_yl0 = (candidateArray[canTotalCnt].motion_vector_yl0 + 2)&~0x03;

    candidateArray[canTotalCnt].motion_vector_xl1 = (candidateArray[canTotalCnt].motion_vector_xl1 + 2)&~0x03;
    candidateArray[canTotalCnt].motion_vector_yl1 = (candidateArray[canTotalCnt].motion_vector_yl1 + 2)&~0x03;

    return;
}

/***************************************
* return true if the MV candidate is already injected
***************************************/
EbBool mrp_is_already_injected_mv_l0(
    ModeDecisionContext *context_ptr,
    int16_t                mv_x,
    int16_t                mv_y,
    uint8_t                ref_type) {
    for (int inter_candidate_index = 0; inter_candidate_index < context_ptr->injected_mv_count_l0; inter_candidate_index++) {
        if (context_ptr->injected_mv_x_l0_array[inter_candidate_index] == mv_x &&
            context_ptr->injected_mv_y_l0_array[inter_candidate_index] == mv_y &&
            context_ptr->injected_ref_type_l0_array[inter_candidate_index] == ref_type) {
            return(EB_TRUE);
        }
    }

    return(EB_FALSE);
}

EbBool mrp_is_already_injected_mv_l1(
    ModeDecisionContext *context_ptr,
    int16_t                mv_x,
    int16_t                mv_y,
    uint8_t                ref_type) {
    for (int inter_candidate_index = 0; inter_candidate_index < context_ptr->injected_mv_count_l1; inter_candidate_index++) {
        if (context_ptr->injected_mv_x_l1_array[inter_candidate_index] == mv_x &&
            context_ptr->injected_mv_y_l1_array[inter_candidate_index] == mv_y &&
            context_ptr->injected_ref_type_l1_array[inter_candidate_index] == ref_type) {
            return(EB_TRUE);
        }
    }

    return(EB_FALSE);
}

EbBool mrp_is_already_injected_mv_bipred(
    ModeDecisionContext *context_ptr,
    int16_t                mv_x_l0,
    int16_t                mv_y_l0,
    int16_t                mv_x_l1,
    int16_t                mv_y_l1,
    uint8_t                ref_type) {
    for (int inter_candidate_index = 0; inter_candidate_index < context_ptr->injected_mv_count_bipred; inter_candidate_index++) {
        if (context_ptr->injected_mv_x_bipred_l0_array[inter_candidate_index] == mv_x_l0 &&
            context_ptr->injected_mv_y_bipred_l0_array[inter_candidate_index] == mv_y_l0 &&
            context_ptr->injected_mv_x_bipred_l1_array[inter_candidate_index] == mv_x_l1 &&
            context_ptr->injected_mv_y_bipred_l1_array[inter_candidate_index] == mv_y_l1 &&
            context_ptr->injected_ref_type_bipred_array[inter_candidate_index] == ref_type) {
            return(EB_TRUE);
        }
    }
    return(EB_FALSE);
}

EbErrorType SetMvpClipMVs(
    ModeDecisionCandidate  *candidate_ptr,
    uint32_t                    cu_origin_x,
    uint32_t                    cu_origin_y,
    uint32_t                    pu_index,
    uint32_t                    tb_size,
    PictureControlSet      *picture_control_set_ptr)
{
    EbErrorType  return_error = EB_ErrorNone;

    uint32_t        picture_width = ((SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr)->seq_header.max_frame_width;
    uint32_t        picture_height = ((SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr)->seq_header.max_frame_height;

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
            &candidate_ptr->motion_vector_xl0,
            &candidate_ptr->motion_vector_yl0,
            picture_width,
            picture_height,
            tb_size);

        break;

    case UNI_PRED_LIST_1:

        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motion_vector_xl1,
            &candidate_ptr->motion_vector_yl1,
            picture_width,
            picture_height,
            tb_size);

        break;

    case BI_PRED:

        // Choose the MVP in list0
        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motion_vector_xl0,
            &candidate_ptr->motion_vector_yl0,
            picture_width,
            picture_height,
            tb_size);

        // Choose the MVP in list1
        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motion_vector_xl1,
            &candidate_ptr->motion_vector_yl1,
            picture_width,
            picture_height,
            tb_size);
        break;

    default:
        break;
    }

    return return_error;
}

void LimitMvOverBound(
    int16_t *mvx,
    int16_t *mvy,
    ModeDecisionContext     *ctxtPtr,
    const SequenceControlSet      *sCSet)
{
    int32_t mvxF, mvyF;

    //L0
    mvxF = (*mvx) >> 2;
    mvyF = (*mvy) >> 2;

    if ((int32_t)ctxtPtr->cu_origin_x + mvxF + (int32_t)ctxtPtr->blk_geom->bwidth > (int32_t)sCSet->seq_header.max_frame_width)
        *mvx = (int16_t)(sCSet->seq_header.max_frame_width - ctxtPtr->blk_geom->bwidth - ctxtPtr->cu_origin_x);
    if ((int32_t)ctxtPtr->cu_origin_y + mvyF + (int32_t)ctxtPtr->blk_geom->bheight > (int32_t)sCSet->seq_header.max_frame_height)
        *mvy = (int16_t)(sCSet->seq_header.max_frame_height - ctxtPtr->blk_geom->bheight - ctxtPtr->cu_origin_y);
    if ((int32_t)ctxtPtr->cu_origin_x + mvxF < 0)
        *mvx = -(int16_t)ctxtPtr->cu_origin_x;
    if ((int32_t)ctxtPtr->cu_origin_y + mvyF < 0)
        *mvy = -(int16_t)ctxtPtr->cu_origin_y;
}

void sort_fast_loop_candidates(
    struct ModeDecisionContext   *context_ptr,
    uint32_t                        buffer_total_count,
    ModeDecisionCandidateBuffer **buffer_ptr_array,
    uint8_t                        *best_candidate_index_array,
    uint8_t                        *sorted_candidate_index_array,
    uint64_t                       *ref_fast_cost) {
    uint32_t fullReconCandidateCount = context_ptr->full_recon_search_count;

    //  move the scratch candidates (MAX_CU_COST) to the last spots (if any)
    uint32_t best_candidate_start_index = 0;
    uint32_t best_candidate_end_index = buffer_total_count - 1;
    for (uint8_t full_buffer_index = 0; full_buffer_index < buffer_total_count; full_buffer_index++) {
        if (*(buffer_ptr_array[full_buffer_index]->fast_cost_ptr) == MAX_CU_COST)
            best_candidate_index_array[best_candidate_end_index--] = full_buffer_index;
        else
            best_candidate_index_array[best_candidate_start_index++] = full_buffer_index;
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

    // fl escape level 2: inter then intra
    for (i = 0; i < fullReconCandidateCount; ++i)
        sorted_candidate_index_array[i] = best_candidate_index_array[i];
    for (i = 0; i < fullReconCandidateCount - 1; ++i) {
        for (j = i + 1; j < fullReconCandidateCount; ++j) {
            if (*(buffer_ptr_array[sorted_candidate_index_array[j]]->fast_cost_ptr) < *(buffer_ptr_array[sorted_candidate_index_array[i]]->fast_cost_ptr)) {
                index = sorted_candidate_index_array[i];
                sorted_candidate_index_array[i] = (uint8_t)sorted_candidate_index_array[j];
                sorted_candidate_index_array[j] = (uint8_t)index;
            }
        }
    }
    // tx search
    *ref_fast_cost = *(buffer_ptr_array[sorted_candidate_index_array[0]]->fast_cost_ptr);
}

#define BIPRED_3x3_REFINMENT_POSITIONS 8

int8_t ALLOW_REFINEMENT_FLAG[BIPRED_3x3_REFINMENT_POSITIONS] = {  1, 0, 1, 0, 1,  0,  1, 0 };
int8_t BIPRED_3x3_X_POS[BIPRED_3x3_REFINMENT_POSITIONS] = { -1, -1, 0, 1, 1, 1, 0, -1 };
int8_t BIPRED_3x3_Y_POS[BIPRED_3x3_REFINMENT_POSITIONS] = { 0, 1, 1, 1, 0, -1, -1, -1 };

void Unipred3x3CandidatesInjection(
    const SequenceControlSet  *sequence_control_set_ptr,
    PictureControlSet         *picture_control_set_ptr,
    ModeDecisionContext       *context_ptr,
    LargestCodingUnit         *sb_ptr,
    uint32_t                   me_sb_addr,
    SsMeContext               *inloop_me_context,
    EbBool                     use_close_loop_me,
    uint32_t                   close_loop_me_index,
    uint32_t                  *candidateTotalCnt){
    UNUSED(sb_ptr);
    uint32_t                   bipredIndex;
    uint32_t                   canTotalCnt = (*candidateTotalCnt);
    const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[context_ptr->me_block_offset];
    const MeCandidate *me_block_results = me_results->me_candidate[context_ptr->me_block_offset];
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    IntMv  bestPredmv[2] = { {0}, {0} };

    // (8 Best_L0 neighbors)
    //const MeLcuResults_t *meResults = pictureControlSetPtr->ParentPcsPtr->meResultsPtr[lcuAddr];
    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
    {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;

        if (inter_direction == 0) {
    for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
    {
        /**************
        NEWMV L0
        ************* */
        if (context_ptr->unipred3x3_injection >= 2){
            if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                continue;
        }
        int16_t to_inject_mv_x = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv + BIPRED_3x3_X_POS[bipredIndex]) << 1;
        int16_t to_inject_mv_y = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
        uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
        if (context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {
            candidateArray[canTotalCnt].type = INTER_MODE;
            candidateArray[canTotalCnt].distortion_ready = 0;
            candidateArray[canTotalCnt].use_intrabc = 0;
            candidateArray[canTotalCnt].merge_flag = EB_FALSE;
            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;
            candidateArray[canTotalCnt].inter_mode = NEWMV;
            candidateArray[canTotalCnt].pred_mode = NEWMV;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

            candidateArray[canTotalCnt].is_compound = 0;
            candidateArray[canTotalCnt].is_new_mv = 1;
            candidateArray[canTotalCnt].is_zero_mv = 0;

            candidateArray[canTotalCnt].drl_index = 0;

            // Set the MV to ME result
            candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x;
            candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y;

            // will be needed later by the rate estimation
            candidateArray[canTotalCnt].ref_mv_index = 0;
            candidateArray[canTotalCnt].pred_mv_weight = 0;
            candidateArray[canTotalCnt].ref_frame_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
            candidateArray[canTotalCnt].ref_frame_index_l0 = list0_ref_index;
            candidateArray[canTotalCnt].ref_frame_index_l1 = -1;

            candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
            candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;

            ChooseBestAv1MvPred(
                context_ptr,
                candidateArray[canTotalCnt].md_rate_estimation_ptr,
                context_ptr->cu_ptr,
                candidateArray[canTotalCnt].ref_frame_type,
                candidateArray[canTotalCnt].is_compound,
                candidateArray[canTotalCnt].pred_mode,
                candidateArray[canTotalCnt].motion_vector_xl0,
                candidateArray[canTotalCnt].motion_vector_yl0,
                0, 0,
                &candidateArray[canTotalCnt].drl_index,
                bestPredmv);

            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
            INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
            context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
            context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
            context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = to_inject_ref_type;
            ++context_ptr->injected_mv_count_l0;
        }
           }
        }
    }

    // (8 Best_L1 neighbors)
//const MeLcuResults_t *meResults = pictureControlSetPtr->ParentPcsPtr->meResultsPtr[lcuAddr];
    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
    {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;
        if (inter_direction == 1) {
    for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
    {
        if (isCompoundEnabled) {
            /**************
            NEWMV L1
            ************* */
            if (context_ptr->unipred3x3_injection >= 2) {
                if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                    continue;
            }
            int16_t to_inject_mv_x = use_close_loop_me ? (inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].x_mv + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            int16_t to_inject_mv_y = use_close_loop_me ? (inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].y_mv + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
            if (context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {
                candidateArray[canTotalCnt].type = INTER_MODE;
                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)1;
                candidateArray[canTotalCnt].inter_mode = NEWMV;
                candidateArray[canTotalCnt].pred_mode = NEWMV;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

                candidateArray[canTotalCnt].is_compound = 0;
                candidateArray[canTotalCnt].is_new_mv = 1;
                candidateArray[canTotalCnt].is_zero_mv = 0;

                candidateArray[canTotalCnt].drl_index = 0;

                // Set the MV to ME result
                candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x;
                candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y;
                // will be needed later by the rate estimation
                candidateArray[canTotalCnt].ref_mv_index = 0;
                candidateArray[canTotalCnt].pred_mv_weight = 0;
                candidateArray[canTotalCnt].ref_frame_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                candidateArray[canTotalCnt].ref_frame_index_l0 = -1;
                candidateArray[canTotalCnt].ref_frame_index_l1 = list1_ref_index;
                candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
                ChooseBestAv1MvPred(
                    context_ptr,
                    candidateArray[canTotalCnt].md_rate_estimation_ptr,
                    context_ptr->cu_ptr,
                    candidateArray[canTotalCnt].ref_frame_type,
                    candidateArray[canTotalCnt].is_compound,
                    candidateArray[canTotalCnt].pred_mode,
                    candidateArray[canTotalCnt].motion_vector_xl1,
                    candidateArray[canTotalCnt].motion_vector_yl1,
                    0, 0,
                    &candidateArray[canTotalCnt].drl_index,
                    bestPredmv);

                candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[0].as_mv.col;
                candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[0].as_mv.row;
                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
                context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
                context_ptr->injected_ref_type_l1_array[context_ptr->injected_mv_count_l1] = to_inject_ref_type;
                ++context_ptr->injected_mv_count_l1;
            }
        }
    }
        }
    }

    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;

    return;
}

void Bipred3x3CandidatesInjection(
    const SequenceControlSet *sequence_control_set_ptr,
    PictureControlSet        *picture_control_set_ptr,
    ModeDecisionContext      *context_ptr,
    LargestCodingUnit        *sb_ptr,
    uint32_t                  me_sb_addr,
    SsMeContext              *inloop_me_context,
    EbBool                    use_close_loop_me,
    uint32_t                  close_loop_me_index,
    uint32_t                 *candidateTotalCnt){
    UNUSED(sb_ptr);
    uint32_t                   bipredIndex;
    uint32_t                   canTotalCnt = (*candidateTotalCnt);
    const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[context_ptr->me_block_offset];
    const MeCandidate *me_block_results = me_results->me_candidate[context_ptr->me_block_offset];
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    IntMv  bestPredmv[2] = { {0}, {0} };

    if (isCompoundEnabled) {
        /**************
       NEW_NEWMV
       ************* */
       //const MeLcuResults_t *meResults = pictureControlSetPtr->ParentPcsPtr->meResultsPtr[lcuAddr];
        for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
        {
            const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
            const uint8_t inter_direction = me_block_results_ptr->direction;
            const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
            const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;

            if (inter_direction == 2) {
       // (Best_L0, 8 Best_L1 neighbors)
        for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
        {
        if (context_ptr->bipred3x3_injection >= 2){
            if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                continue;
        }
        int16_t to_inject_mv_x_l0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1;
        int16_t to_inject_mv_y_l0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1;
        int16_t to_inject_mv_x_l1 = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1) : (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv + BIPRED_3x3_X_POS[bipredIndex]) << 1;
        int16_t to_inject_mv_y_l1 = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1) : (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv + BIPRED_3x3_Y_POS[bipredIndex]) << 1;

        MvReferenceFrame rf[2];
        rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
        rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
        uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
        if (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE) {
                candidateArray[canTotalCnt].type = INTER_MODE;
                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                candidateArray[canTotalCnt].is_new_mv = 1;
                candidateArray[canTotalCnt].is_zero_mv = 0;

            candidateArray[canTotalCnt].drl_index = 0;

            // Set the MV to ME result
            candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x_l0;
            candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y_l0;
            candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x_l1;
            candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y_l1;
            // will be needed later by the rate estimation
            candidateArray[canTotalCnt].ref_mv_index = 0;
            candidateArray[canTotalCnt].pred_mv_weight = 0;

            candidateArray[canTotalCnt].inter_mode = NEW_NEWMV;
            candidateArray[canTotalCnt].pred_mode = NEW_NEWMV;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            candidateArray[canTotalCnt].is_compound = 1;
            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
            MvReferenceFrame rf[2];
            rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
            rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
            candidateArray[canTotalCnt].ref_frame_type = av1_ref_frame_type(rf);
            candidateArray[canTotalCnt].ref_frame_index_l0 = list0_ref_index;
            candidateArray[canTotalCnt].ref_frame_index_l1 = list1_ref_index;
            candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
            candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
            ChooseBestAv1MvPred(
                context_ptr,
                candidateArray[canTotalCnt].md_rate_estimation_ptr,
                context_ptr->cu_ptr,
                candidateArray[canTotalCnt].ref_frame_type,
                candidateArray[canTotalCnt].is_compound,
                candidateArray[canTotalCnt].pred_mode,
                candidateArray[canTotalCnt].motion_vector_xl0,
                candidateArray[canTotalCnt].motion_vector_yl0,
                candidateArray[canTotalCnt].motion_vector_xl1,
                candidateArray[canTotalCnt].motion_vector_yl1,
                &candidateArray[canTotalCnt].drl_index,
                bestPredmv);

            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
            candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
            candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;
            INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
            context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
            context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
            context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
            context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
            context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = to_inject_ref_type;
            ++context_ptr->injected_mv_count_bipred;
        }
        }

        // (8 Best_L0 neighbors, Best_L1) :
        for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
        {
            if (context_ptr->bipred3x3_injection >= 2){
                if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                    continue;
            }
            int16_t to_inject_mv_x_l0 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] + BIPRED_3x3_X_POS[bipredIndex]) << 1 : (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            int16_t to_inject_mv_y_l0 = use_close_loop_me ? (inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] + BIPRED_3x3_Y_POS[bipredIndex]) << 1 : (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
            int16_t to_inject_mv_x_l1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv << 1;
            int16_t to_inject_mv_y_l1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv << 1;

            MvReferenceFrame rf[2];
            rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
            rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
            uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
            if (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE) {
                candidateArray[canTotalCnt].type = INTER_MODE;
                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;

                candidateArray[canTotalCnt].is_new_mv = 1;
                candidateArray[canTotalCnt].is_zero_mv = 0;

                candidateArray[canTotalCnt].drl_index = 0;

                // Set the MV to ME result
                candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x_l0;
                candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y_l0;
                candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x_l1;
                candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y_l1;
                // will be needed later by the rate estimation
                candidateArray[canTotalCnt].ref_mv_index = 0;
                candidateArray[canTotalCnt].pred_mv_weight = 0;

                candidateArray[canTotalCnt].inter_mode = NEW_NEWMV;
                candidateArray[canTotalCnt].pred_mode = NEW_NEWMV;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canTotalCnt].is_compound = 1;
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
                MvReferenceFrame rf[2];
                rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
                rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
                candidateArray[canTotalCnt].ref_frame_type = av1_ref_frame_type(rf);
                candidateArray[canTotalCnt].ref_frame_index_l0 = list0_ref_index;
                candidateArray[canTotalCnt].ref_frame_index_l1 = list1_ref_index;
                candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
                ChooseBestAv1MvPred(
                    context_ptr,
                    candidateArray[canTotalCnt].md_rate_estimation_ptr,
                    context_ptr->cu_ptr,
                    candidateArray[canTotalCnt].ref_frame_type,
                    candidateArray[canTotalCnt].is_compound,
                    candidateArray[canTotalCnt].pred_mode,
                    candidateArray[canTotalCnt].motion_vector_xl0,
                    candidateArray[canTotalCnt].motion_vector_yl0,
                    candidateArray[canTotalCnt].motion_vector_xl1,
                    candidateArray[canTotalCnt].motion_vector_yl1,
                    &candidateArray[canTotalCnt].drl_index,
                    bestPredmv);

                candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
                candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
                candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;
                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = to_inject_ref_type;
                ++context_ptr->injected_mv_count_bipred;
            }
        }
            }
        }
     }

    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;

    return;
}
#if EIGTH_PEL_MV
void eighth_pel_unipred_refinement(
    const SequenceControlSet  *sequence_control_set_ptr,
    PictureControlSet         *picture_control_set_ptr,
    ModeDecisionContext       *context_ptr,
    uint32_t                   me_sb_addr,
    SsMeContext               *inloop_me_context,
    EbBool                     use_close_loop_me,
    uint32_t                   close_loop_me_index,
    uint32_t                  *candidateTotalCnt) {
    uint32_t                   bipredIndex;
    uint32_t                   canTotalCnt = (*candidateTotalCnt);
    const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[context_ptr->me_block_offset];
    const MeCandidate *me_block_results = me_results->me_candidate[context_ptr->me_block_offset];
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    IntMv  bestPredmv[2] = { {0}, {0} };

    // (8 Best_L0 neighbors)
    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index) {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
        if (inter_direction == 0) {
            for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex) {
                /**************
                NEWMV L0
                ************* */
                if (context_ptr->unipred3x3_injection >= 2) {
                    if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                        continue;
                }

                int16_t to_inject_mv_x = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] << 1) + BIPRED_3x3_X_POS[bipredIndex])  : ((me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1) + BIPRED_3x3_X_POS[bipredIndex]);
                int16_t to_inject_mv_y = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] << 1) + BIPRED_3x3_Y_POS[bipredIndex])  : ((me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1) + BIPRED_3x3_Y_POS[bipredIndex]);
                uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);

                if (context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {
                    candidateArray[canTotalCnt].type = INTER_MODE;
                    candidateArray[canTotalCnt].distortion_ready = 0;
                    candidateArray[canTotalCnt].use_intrabc = 0;
                    candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                    candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;
                    candidateArray[canTotalCnt].inter_mode = NEWMV;
                    candidateArray[canTotalCnt].pred_mode = NEWMV;
                    candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                    candidateArray[canTotalCnt].is_compound = 0;
                    candidateArray[canTotalCnt].is_new_mv = 1;
                    candidateArray[canTotalCnt].is_zero_mv = 0;
                    candidateArray[canTotalCnt].drl_index = 0;
                    // Set the MV to ME result
                    candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x;
                    candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y;
                    // will be needed later by the rate estimation
                    candidateArray[canTotalCnt].ref_mv_index = 0;
                    candidateArray[canTotalCnt].pred_mv_weight = 0;
                    candidateArray[canTotalCnt].ref_frame_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
                    candidateArray[canTotalCnt].ref_frame_index_l0 = list0_ref_index;
                    candidateArray[canTotalCnt].ref_frame_index_l1 = -1;
                    candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                    candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

                    ChooseBestAv1MvPred(
                        context_ptr,
                        candidateArray[canTotalCnt].md_rate_estimation_ptr,
                        context_ptr->cu_ptr,
                        candidateArray[canTotalCnt].ref_frame_type,
                        candidateArray[canTotalCnt].is_compound,
                        candidateArray[canTotalCnt].pred_mode,
                        candidateArray[canTotalCnt].motion_vector_xl0,
                        candidateArray[canTotalCnt].motion_vector_yl0,
                        0, 0,
                        &candidateArray[canTotalCnt].drl_index,
                        bestPredmv);

                    candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                    candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
                    INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                    context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
                    context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
                    context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = to_inject_ref_type;

                    ++context_ptr->injected_mv_count_l0;
                }
            }
        }
    }

    // (8 Best_L1 neighbors)
    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index) {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;
        if (inter_direction == 1) {
            for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex) {
                if (isCompoundEnabled) {
                    /**************
                    NEWMV L1
                    ************* */
                    if (context_ptr->unipred3x3_injection >= 2) {
                        if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                            continue;
                    }
                    int16_t to_inject_mv_x = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] << 1) + BIPRED_3x3_X_POS[bipredIndex])  : ((me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].x_mv << 1) + BIPRED_3x3_X_POS[bipredIndex]);
                    int16_t to_inject_mv_y = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] << 1) + BIPRED_3x3_Y_POS[bipredIndex])  : ((me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].y_mv << 1) + BIPRED_3x3_Y_POS[bipredIndex]);
                    uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                    if (context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {
                        candidateArray[canTotalCnt].type = INTER_MODE;
                        candidateArray[canTotalCnt].distortion_ready = 0;
                        candidateArray[canTotalCnt].use_intrabc = 0;
                        candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                        candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)1;
                        candidateArray[canTotalCnt].inter_mode = NEWMV;
                        candidateArray[canTotalCnt].pred_mode = NEWMV;
                        candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                        candidateArray[canTotalCnt].is_compound = 0;
                        candidateArray[canTotalCnt].is_new_mv = 1;
                        candidateArray[canTotalCnt].is_zero_mv = 0;
                        candidateArray[canTotalCnt].drl_index = 0;
                        // Set the MV to ME result
                        candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x;
                        candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y;
                        // will be needed later by the rate estimation
                        candidateArray[canTotalCnt].ref_mv_index = 0;
                        candidateArray[canTotalCnt].pred_mv_weight = 0;
                        candidateArray[canTotalCnt].ref_frame_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                        candidateArray[canTotalCnt].ref_frame_index_l0 = -1;
                        candidateArray[canTotalCnt].ref_frame_index_l1 = list1_ref_index;
                        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

                        ChooseBestAv1MvPred(
                            context_ptr,
                            candidateArray[canTotalCnt].md_rate_estimation_ptr,
                            context_ptr->cu_ptr,
                            candidateArray[canTotalCnt].ref_frame_type,
                            candidateArray[canTotalCnt].is_compound,
                            candidateArray[canTotalCnt].pred_mode,
                            candidateArray[canTotalCnt].motion_vector_xl1,
                            candidateArray[canTotalCnt].motion_vector_yl1,
                            0, 0,
                            &candidateArray[canTotalCnt].drl_index,
                            bestPredmv);

                        candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[0].as_mv.col;
                        candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[0].as_mv.row;
                        INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                        context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
                        context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
                        context_ptr->injected_ref_type_l1_array[context_ptr->injected_mv_count_l1] = to_inject_ref_type;
                        ++context_ptr->injected_mv_count_l1;
                    }
                }
            }
        }
    }
    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;
    return;
}

void eighth_pel_bipred_refinement(
    const SequenceControlSet *sequence_control_set_ptr,
    PictureControlSet        *picture_control_set_ptr,
    ModeDecisionContext      *context_ptr,
    uint32_t                  me_sb_addr,
    SsMeContext              *inloop_me_context,
    EbBool                    use_close_loop_me,
    uint32_t                  close_loop_me_index,
    uint32_t                 *candidateTotalCnt) {
    uint32_t                   bipredIndex;
    uint32_t                   canTotalCnt = (*candidateTotalCnt);
    const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[context_ptr->me_block_offset];
    const MeCandidate *me_block_results = me_results->me_candidate[context_ptr->me_block_offset];
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    IntMv  bestPredmv[2] = { {0}, {0} };
    MvReferenceFrame rf[2];
    if (isCompoundEnabled) {
        /**************
       NEW_NEWMV
       ************* */
        for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index) {
            const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
            const uint8_t inter_direction = me_block_results_ptr->direction;
            const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
            const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;

            if (inter_direction == 2) {
                // (Best_L0, 8 Best_L1 neighbors)
                for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex) {
                    if (context_ptr->bipred3x3_injection >= 2) {
                        if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                            continue;
                    }
                    int16_t to_inject_mv_x_l0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1;
                    int16_t to_inject_mv_y_l0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1;
                    int16_t to_inject_mv_x_l1 = use_close_loop_me ? (((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] << 1) + BIPRED_3x3_X_POS[bipredIndex])) : ((me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv << 1) + BIPRED_3x3_X_POS[bipredIndex]);
                    int16_t to_inject_mv_y_l1 = use_close_loop_me ? (((inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] << 1) + BIPRED_3x3_Y_POS[bipredIndex])) : ((me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv << 1) + BIPRED_3x3_Y_POS[bipredIndex]);
                    rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
                    rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
                    uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
                    if (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE) {
                        candidateArray[canTotalCnt].type = INTER_MODE;
                        candidateArray[canTotalCnt].distortion_ready = 0;
                        candidateArray[canTotalCnt].use_intrabc = 0;
                        candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                        candidateArray[canTotalCnt].is_new_mv = 1;
                        candidateArray[canTotalCnt].is_zero_mv = 0;
                        candidateArray[canTotalCnt].drl_index = 0;
                        // Set the MV to ME result
                        candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x_l0;
                        candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y_l0;
                        candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x_l1;
                        candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y_l1;
                        // will be needed later by the rate estimation
                        candidateArray[canTotalCnt].ref_mv_index = 0;
                        candidateArray[canTotalCnt].pred_mv_weight = 0;
                        candidateArray[canTotalCnt].inter_mode = NEW_NEWMV;
                        candidateArray[canTotalCnt].pred_mode = NEW_NEWMV;
                        candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                        candidateArray[canTotalCnt].is_compound = 1;
                        candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
                        rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
                        rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
                        candidateArray[canTotalCnt].ref_frame_type = av1_ref_frame_type(rf);
                        candidateArray[canTotalCnt].ref_frame_index_l0 = list0_ref_index;
                        candidateArray[canTotalCnt].ref_frame_index_l1 = list1_ref_index;
                        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

                        ChooseBestAv1MvPred(
                            context_ptr,
                            candidateArray[canTotalCnt].md_rate_estimation_ptr,
                            context_ptr->cu_ptr,
                            candidateArray[canTotalCnt].ref_frame_type,
                            candidateArray[canTotalCnt].is_compound,
                            candidateArray[canTotalCnt].pred_mode,
                            candidateArray[canTotalCnt].motion_vector_xl0,
                            candidateArray[canTotalCnt].motion_vector_yl0,
                            candidateArray[canTotalCnt].motion_vector_xl1,
                            candidateArray[canTotalCnt].motion_vector_yl1,
                            &candidateArray[canTotalCnt].drl_index,
                            bestPredmv);

                        candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                        candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
                        candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
                        candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;
                        INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                        context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                        context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                        context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                        context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                        context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = to_inject_ref_type;
                        ++context_ptr->injected_mv_count_bipred;
                    }
                }

                // (8 Best_L0 neighbors, Best_L1) :
                for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
                {
                    if (context_ptr->bipred3x3_injection >= 2) {
                        if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                            continue;
                    }
                    int16_t to_inject_mv_x_l0 = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] << 1) + BIPRED_3x3_X_POS[bipredIndex]) : ((me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1) + BIPRED_3x3_X_POS[bipredIndex]);
                    int16_t to_inject_mv_y_l0 = use_close_loop_me ? ((inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] << 1) + BIPRED_3x3_Y_POS[bipredIndex]) : ((me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1) + BIPRED_3x3_Y_POS[bipredIndex]);
                    int16_t to_inject_mv_x_l1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv << 1;
                    int16_t to_inject_mv_y_l1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv << 1;
                    rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
                    rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
                    uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
                    if (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE) {
                        candidateArray[canTotalCnt].type = INTER_MODE;
                        candidateArray[canTotalCnt].distortion_ready = 0;
                        candidateArray[canTotalCnt].use_intrabc = 0;
                        candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                        candidateArray[canTotalCnt].is_new_mv = 1;
                        candidateArray[canTotalCnt].is_zero_mv = 0;
                        candidateArray[canTotalCnt].drl_index = 0;
                        // Set the MV to ME result
                        candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x_l0;
                        candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y_l0;
                        candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x_l1;
                        candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y_l1;
                        // will be needed later by the rate estimation
                        candidateArray[canTotalCnt].ref_mv_index = 0;
                        candidateArray[canTotalCnt].pred_mv_weight = 0;
                        candidateArray[canTotalCnt].inter_mode = NEW_NEWMV;
                        candidateArray[canTotalCnt].pred_mode = NEW_NEWMV;
                        candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                        candidateArray[canTotalCnt].is_compound = 1;
                        candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
                        rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
                        rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
                        candidateArray[canTotalCnt].ref_frame_type = av1_ref_frame_type(rf);
                        candidateArray[canTotalCnt].ref_frame_index_l0 = list0_ref_index;
                        candidateArray[canTotalCnt].ref_frame_index_l1 = list1_ref_index;
                        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_Y] = DCT_DCT;
                        candidateArray[canTotalCnt].transform_type[PLANE_TYPE_UV] = DCT_DCT;

                        ChooseBestAv1MvPred(
                            context_ptr,
                            candidateArray[canTotalCnt].md_rate_estimation_ptr,
                            context_ptr->cu_ptr,
                            candidateArray[canTotalCnt].ref_frame_type,
                            candidateArray[canTotalCnt].is_compound,
                            candidateArray[canTotalCnt].pred_mode,
                            candidateArray[canTotalCnt].motion_vector_xl0,
                            candidateArray[canTotalCnt].motion_vector_yl0,
                            candidateArray[canTotalCnt].motion_vector_xl1,
                            candidateArray[canTotalCnt].motion_vector_yl1,
                            &candidateArray[canTotalCnt].drl_index,
                            bestPredmv);

                        candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                        candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
                        candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
                        candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;
                        INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                        context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                        context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                        context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                        context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                        context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = to_inject_ref_type;
                        ++context_ptr->injected_mv_count_bipred;
                    }
                }
            }
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
void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);
void inject_mvp_candidates_II(
    struct ModeDecisionContext     *context_ptr,
    PictureControlSet              *picture_control_set_ptr,
    CodingUnit                     *cu_ptr,
    MvReferenceFrame                 ref_pair,
    uint32_t                         *candTotCnt)
{
    EbBool allow_compound = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE || context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : EB_TRUE;
    uint8_t inj_mv;
    uint32_t                   canIdx = *candTotCnt;
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    MacroBlockD  *xd = cu_ptr->av1xd;
    uint8_t drli, maxDrlIndex;
    IntMv    nearestmv[2], nearmv[2], ref_mv[2];

    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_pair);

    //single ref/list
    if (rf[1] == NONE_FRAME)
    {
        MvReferenceFrame frame_type = rf[0];
        uint8_t list_idx = get_list_idx(rf[0]);
        uint8_t ref_idx = get_ref_frame_idx(rf[0]);

        //NEAREST
        int16_t to_inject_mv_x = context_ptr->cu_ptr->ref_mvs[frame_type][0].as_mv.col;
        int16_t to_inject_mv_y = context_ptr->cu_ptr->ref_mvs[frame_type][0].as_mv.row;

        inj_mv = list_idx == 0 ?
            context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, frame_type) == EB_FALSE :
            context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, frame_type) == EB_FALSE;

        if (inj_mv) {
            candidateArray[canIdx].type = INTER_MODE;
            candidateArray[canIdx].inter_mode = NEARESTMV;
            candidateArray[canIdx].pred_mode = NEARESTMV;
            candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
            candidateArray[canIdx].is_compound = 0;
            candidateArray[canIdx].distortion_ready = 0;
            candidateArray[canIdx].use_intrabc = 0;
            candidateArray[canIdx].merge_flag = EB_FALSE;
            candidateArray[canIdx].prediction_direction[0] = list_idx;
            candidateArray[canIdx].is_new_mv = 0;
            candidateArray[canIdx].is_zero_mv = 0;

            candidateArray[canIdx].drl_index = 0;
            candidateArray[canIdx].ref_mv_index = 0;
            candidateArray[canIdx].pred_mv_weight = 0;
            candidateArray[canIdx].ref_frame_type = frame_type;

            candidateArray[canIdx].ref_frame_index_l0 = (list_idx == 0) ? ref_idx : -1;
            candidateArray[canIdx].ref_frame_index_l1 = (list_idx == 1) ? ref_idx : -1;
            candidateArray[canIdx].transform_type[0] = DCT_DCT;
            candidateArray[canIdx].transform_type_uv = DCT_DCT;
            if (list_idx == 0) {
                candidateArray[canIdx].motion_vector_xl0 = to_inject_mv_x;
                candidateArray[canIdx].motion_vector_yl0 = to_inject_mv_y;
                context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
                context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
                context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = frame_type;
                ++context_ptr->injected_mv_count_l0;
            }
            else {
                candidateArray[canIdx].motion_vector_xl1 = to_inject_mv_x;
                candidateArray[canIdx].motion_vector_yl1 = to_inject_mv_y;
                context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
                context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
                context_ptr->injected_ref_type_l1_array[context_ptr->injected_mv_count_l1] = frame_type;
                ++context_ptr->injected_mv_count_l1;
            }

            INCRMENT_CAND_TOTAL_COUNT(canIdx);
        }

        //NEAR
        maxDrlIndex = GetMaxDrlIndex(xd->ref_mv_count[frame_type], NEARMV);

        for (drli = 0; drli < maxDrlIndex; drli++)
        {
            get_av1_mv_pred_drl(
                context_ptr,
                cu_ptr,
                frame_type,
                0,
                NEARMV,
                drli,
                nearestmv,
                nearmv,
                ref_mv);

            int16_t to_inject_mv_x = nearmv[0].as_mv.col;
            int16_t to_inject_mv_y = nearmv[0].as_mv.row;

            inj_mv = list_idx == 0 ?
                context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, frame_type) == EB_FALSE :
                context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, frame_type) == EB_FALSE;

            if (inj_mv) {
                candidateArray[canIdx].type = INTER_MODE;
                candidateArray[canIdx].inter_mode = NEARMV;
                candidateArray[canIdx].pred_mode = NEARMV;
                candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canIdx].is_compound = 0;
                candidateArray[canIdx].distortion_ready = 0;
                candidateArray[canIdx].use_intrabc = 0;
                candidateArray[canIdx].merge_flag = EB_FALSE;
                candidateArray[canIdx].prediction_direction[0] = list_idx;
                candidateArray[canIdx].is_new_mv = 0;
                candidateArray[canIdx].is_zero_mv = 0;
                candidateArray[canIdx].drl_index = drli;
                candidateArray[canIdx].ref_mv_index = 0;
                candidateArray[canIdx].pred_mv_weight = 0;
                candidateArray[canIdx].ref_frame_type = frame_type;

                candidateArray[canIdx].ref_frame_index_l0 = (list_idx == 0) ? ref_idx : -1;
                candidateArray[canIdx].ref_frame_index_l1 = (list_idx == 1) ? ref_idx : -1;

                candidateArray[canIdx].transform_type[0] = DCT_DCT;
                candidateArray[canIdx].transform_type_uv = DCT_DCT;
                if (list_idx == 0) {
                    candidateArray[canIdx].motion_vector_xl0 = to_inject_mv_x;
                    candidateArray[canIdx].motion_vector_yl0 = to_inject_mv_y;
                    context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
                    context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
                    context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = frame_type;
                    ++context_ptr->injected_mv_count_l0;
                }
                else {
                    candidateArray[canIdx].motion_vector_xl1 = to_inject_mv_x;
                    candidateArray[canIdx].motion_vector_yl1 = to_inject_mv_y;
                    context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
                    context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
                    context_ptr->injected_ref_type_l1_array[context_ptr->injected_mv_count_l1] = frame_type;
                    ++context_ptr->injected_mv_count_l1;
                }
                INCRMENT_CAND_TOTAL_COUNT(canIdx);
            }
        }
    }
    else if (allow_compound)
    {
        uint8_t ref_idx_0 = get_ref_frame_idx(rf[0]);
        uint8_t ref_idx_1 = get_ref_frame_idx(rf[1]);

        {
            //NEAREST_NEAREST
            int16_t to_inject_mv_x_l0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].this_mv.as_mv.col;
            int16_t to_inject_mv_y_l0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].this_mv.as_mv.row;
            int16_t to_inject_mv_x_l1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].comp_mv.as_mv.col;
            int16_t to_inject_mv_y_l1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].comp_mv.as_mv.row;

            inj_mv = context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, ref_pair) == EB_FALSE;

            if (inj_mv) {
                candidateArray[canIdx].type = INTER_MODE;
                candidateArray[canIdx].inter_mode = NEAREST_NEARESTMV;
                candidateArray[canIdx].pred_mode = NEAREST_NEARESTMV;
                candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canIdx].is_compound = 1;
                candidateArray[canIdx].distortion_ready = 0;
                candidateArray[canIdx].use_intrabc = 0;

                candidateArray[canIdx].merge_flag =
                    picture_control_set_ptr->parent_pcs_ptr->is_skip_mode_allowed &&
                    (rf[0] == picture_control_set_ptr->parent_pcs_ptr->skip_mode_info.ref_frame_idx_0 + 1) &&
                    (rf[1] == picture_control_set_ptr->parent_pcs_ptr->skip_mode_info.ref_frame_idx_1 + 1) ? EB_TRUE : EB_FALSE;

                candidateArray[canIdx].prediction_direction[0] = BI_PRED;
                candidateArray[canIdx].is_new_mv = 0;
                candidateArray[canIdx].is_zero_mv = 0;
                candidateArray[canIdx].motion_vector_xl0 = to_inject_mv_x_l0;
                candidateArray[canIdx].motion_vector_yl0 = to_inject_mv_y_l0;
                candidateArray[canIdx].motion_vector_xl1 = to_inject_mv_x_l1;
                candidateArray[canIdx].motion_vector_yl1 = to_inject_mv_y_l1;
                candidateArray[canIdx].drl_index = 0;
                candidateArray[canIdx].ref_mv_index = 0;
                candidateArray[canIdx].pred_mv_weight = 0;
                candidateArray[canIdx].ref_frame_type = ref_pair;
                candidateArray[canIdx].ref_frame_index_l0 = ref_idx_0;
                candidateArray[canIdx].ref_frame_index_l1 = ref_idx_1;

                candidateArray[canIdx].transform_type[0] = DCT_DCT;
                candidateArray[canIdx].transform_type_uv = DCT_DCT;

                context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = ref_pair;
                ++context_ptr->injected_mv_count_bipred;

                INCRMENT_CAND_TOTAL_COUNT(canIdx);
            }

            //NEAR_NEAR
            maxDrlIndex = GetMaxDrlIndex(xd->ref_mv_count[ref_pair], NEAR_NEARMV);
            for (drli = 0; drli < maxDrlIndex; drli++) {
                get_av1_mv_pred_drl(
                    context_ptr,
                    cu_ptr,
                    ref_pair,
                    1,
                    NEAR_NEARMV,
                    drli,
                    nearestmv,
                    nearmv,
                    ref_mv);

                int16_t to_inject_mv_x_l0 = nearmv[0].as_mv.col;
                int16_t to_inject_mv_y_l0 = nearmv[0].as_mv.row;
                int16_t to_inject_mv_x_l1 = nearmv[1].as_mv.col;
                int16_t to_inject_mv_y_l1 = nearmv[1].as_mv.row;

                inj_mv = context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, ref_pair) == EB_FALSE;

                if (inj_mv) {
                    candidateArray[canIdx].type = INTER_MODE;
                    candidateArray[canIdx].inter_mode = NEAR_NEARMV;
                    candidateArray[canIdx].pred_mode = NEAR_NEARMV;
                    candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                    candidateArray[canIdx].is_compound = 1;
                    candidateArray[canIdx].distortion_ready = 0;
                    candidateArray[canIdx].use_intrabc = 0;
                    candidateArray[canIdx].merge_flag = EB_FALSE;
                    candidateArray[canIdx].prediction_direction[0] = BI_PRED;
                    candidateArray[canIdx].is_new_mv = 0;
                    candidateArray[canIdx].is_zero_mv = 0;

                    candidateArray[canIdx].motion_vector_xl0 = to_inject_mv_x_l0;
                    candidateArray[canIdx].motion_vector_yl0 = to_inject_mv_y_l0;
                    candidateArray[canIdx].motion_vector_xl1 = to_inject_mv_x_l1;
                    candidateArray[canIdx].motion_vector_yl1 = to_inject_mv_y_l1;

                    candidateArray[canIdx].drl_index = drli;
                    candidateArray[canIdx].ref_mv_index = 0;
                    candidateArray[canIdx].pred_mv_weight = 0;

                    candidateArray[canIdx].ref_frame_type = ref_pair;

                    candidateArray[canIdx].ref_frame_index_l0 = ref_idx_0;
                    candidateArray[canIdx].ref_frame_index_l1 = ref_idx_1;

                    candidateArray[canIdx].transform_type[0] = DCT_DCT;
                    candidateArray[canIdx].transform_type_uv = DCT_DCT;

                    context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                    context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                    context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                    context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                    context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = ref_pair;
                    ++context_ptr->injected_mv_count_bipred;

                    INCRMENT_CAND_TOTAL_COUNT(canIdx);
                }
            }
        }
    }

    //update tot Candidate count
    *candTotCnt = canIdx;
}

void inject_new_nearest_new_comb_candidates(
    const SequenceControlSet       *sequence_control_set_ptr,
    struct ModeDecisionContext     *context_ptr,
    PictureControlSet              *picture_control_set_ptr,
    MvReferenceFrame                ref_pair,
    uint32_t                       *candTotCnt)
{
    uint8_t inj_mv;
    uint32_t                  canIdx = *candTotCnt;
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;

    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_pair);

    {
        uint8_t ref_idx_0 = get_ref_frame_idx(rf[0]);
        uint8_t ref_idx_1 = get_ref_frame_idx(rf[1]);

        if (rf[1] != NONE_FRAME)
        {
            {
                //NEAREST_NEWMV
                const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[context_ptr->me_sb_addr];

                int16_t to_inject_mv_x_l0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].this_mv.as_mv.col;
                int16_t to_inject_mv_y_l0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].this_mv.as_mv.row;
                int16_t to_inject_mv_x_l1 = me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (get_list_idx(rf[1]) << 2) : (get_list_idx(rf[1]) << 1)) + ref_idx_1].x_mv << 1;
                int16_t to_inject_mv_y_l1 = me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (get_list_idx(rf[1]) << 2) : (get_list_idx(rf[1]) << 1)) + ref_idx_1].y_mv << 1;

                inj_mv = context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, ref_pair) == EB_FALSE;

                if (inj_mv) {

                    candidateArray[canIdx].type = INTER_MODE;
                    candidateArray[canIdx].inter_mode = NEAREST_NEWMV;
                    candidateArray[canIdx].pred_mode = NEAREST_NEWMV;
                    candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                    candidateArray[canIdx].is_compound = 1;
                    candidateArray[canIdx].distortion_ready = 0;
                    candidateArray[canIdx].use_intrabc = 0;

                    candidateArray[canIdx].merge_flag = EB_FALSE;

                    candidateArray[canIdx].prediction_direction[0] = BI_PRED;
                    candidateArray[canIdx].is_new_mv = 0;
                    candidateArray[canIdx].is_zero_mv = 0;
                    candidateArray[canIdx].motion_vector_xl0 = to_inject_mv_x_l0;
                    candidateArray[canIdx].motion_vector_yl0 = to_inject_mv_y_l0;
                    candidateArray[canIdx].motion_vector_xl1 = to_inject_mv_x_l1;
                    candidateArray[canIdx].motion_vector_yl1 = to_inject_mv_y_l1;
                    candidateArray[canIdx].drl_index = 0;
                    candidateArray[canIdx].ref_mv_index = 0;
                    candidateArray[canIdx].pred_mv_weight = 0;
                    candidateArray[canIdx].ref_frame_type = ref_pair;
                    candidateArray[canIdx].ref_frame_index_l0 = ref_idx_0;
                    candidateArray[canIdx].ref_frame_index_l1 = ref_idx_1;
                    candidateArray[canIdx].transform_type[0] = DCT_DCT;
                    candidateArray[canIdx].transform_type_uv = DCT_DCT;

                    IntMv  bestPredmv[2] = { {0}, {0} };

                    ChooseBestAv1MvPred(
                        context_ptr,
                        candidateArray[canIdx].md_rate_estimation_ptr,
                        context_ptr->cu_ptr,
                        candidateArray[canIdx].ref_frame_type,
                        candidateArray[canIdx].is_compound,
                        candidateArray[canIdx].pred_mode,
                        candidateArray[canIdx].motion_vector_xl0,
                        candidateArray[canIdx].motion_vector_yl0,
                        candidateArray[canIdx].motion_vector_xl1,
                        candidateArray[canIdx].motion_vector_yl1,
                        &candidateArray[canIdx].drl_index,
                        bestPredmv);

                    candidateArray[canIdx].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                    candidateArray[canIdx].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
                    candidateArray[canIdx].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
                    candidateArray[canIdx].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;

                    context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                    context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                    context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                    context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                    context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = ref_pair;
                    ++context_ptr->injected_mv_count_bipred;

                    INCRMENT_CAND_TOTAL_COUNT(canIdx);
                }
            }

            {
                //NEW_NEARESTMV
                const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[context_ptr->me_sb_addr];

                int16_t to_inject_mv_x_l0 = me_results->me_mv_array[context_ptr->me_block_offset][ref_idx_0].x_mv << 1;
                int16_t to_inject_mv_y_l0 = me_results->me_mv_array[context_ptr->me_block_offset][ref_idx_0].y_mv << 1;
                int16_t to_inject_mv_x_l1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].comp_mv.as_mv.col;
                int16_t to_inject_mv_y_l1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].comp_mv.as_mv.row;

                inj_mv = context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, ref_pair) == EB_FALSE;

                if (inj_mv)
                {

                    candidateArray[canIdx].type = INTER_MODE;
                    candidateArray[canIdx].inter_mode = NEW_NEARESTMV;
                    candidateArray[canIdx].pred_mode = NEW_NEARESTMV;
                    candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                    candidateArray[canIdx].is_compound = 1;
                    candidateArray[canIdx].distortion_ready = 0;
                    candidateArray[canIdx].use_intrabc = 0;

                    candidateArray[canIdx].merge_flag = EB_FALSE;

                    candidateArray[canIdx].prediction_direction[0] = BI_PRED;
                    candidateArray[canIdx].is_new_mv = 0;
                    candidateArray[canIdx].is_zero_mv = 0;
                    candidateArray[canIdx].motion_vector_xl0 = to_inject_mv_x_l0;
                    candidateArray[canIdx].motion_vector_yl0 = to_inject_mv_y_l0;
                    candidateArray[canIdx].motion_vector_xl1 = to_inject_mv_x_l1;
                    candidateArray[canIdx].motion_vector_yl1 = to_inject_mv_y_l1;
                    candidateArray[canIdx].drl_index = 0;
                    candidateArray[canIdx].ref_mv_index = 0;
                    candidateArray[canIdx].pred_mv_weight = 0;
                    candidateArray[canIdx].ref_frame_type = ref_pair;
                    candidateArray[canIdx].ref_frame_index_l0 = ref_idx_0;
                    candidateArray[canIdx].ref_frame_index_l1 = ref_idx_1;
                    candidateArray[canIdx].transform_type[0] = DCT_DCT;
                    candidateArray[canIdx].transform_type_uv = DCT_DCT;

                    IntMv  bestPredmv[2] = { {0}, {0} };

                    ChooseBestAv1MvPred(
                        context_ptr,
                        candidateArray[canIdx].md_rate_estimation_ptr,
                        context_ptr->cu_ptr,
                        candidateArray[canIdx].ref_frame_type,
                        candidateArray[canIdx].is_compound,
                        candidateArray[canIdx].pred_mode,
                        candidateArray[canIdx].motion_vector_xl0,
                        candidateArray[canIdx].motion_vector_yl0,
                        candidateArray[canIdx].motion_vector_xl1,
                        candidateArray[canIdx].motion_vector_yl1,
                        &candidateArray[canIdx].drl_index,
                        bestPredmv);

                    candidateArray[canIdx].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                    candidateArray[canIdx].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
                    candidateArray[canIdx].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
                    candidateArray[canIdx].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;

                    context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                    context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                    context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                    context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                    context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = ref_pair;
                    ++context_ptr->injected_mv_count_bipred;

                    INCRMENT_CAND_TOTAL_COUNT(canIdx);
                }
            }

        }
    }

    //update tot Candidate count
    *candTotCnt = canIdx;
}

void inject_warped_motion_candidates(
    PictureControlSet              *picture_control_set_ptr,
    struct ModeDecisionContext     *context_ptr,
    CodingUnit                     *cu_ptr,
    uint32_t                       *candTotCnt,
    SsMeContext                    *ss_mecontext,
    MeLcuResults                   *meResult,
    EbBool                          use_close_loop_me,
    uint32_t                        close_loop_me_index)
{
    uint32_t canIdx = *candTotCnt;
    ModeDecisionCandidate *candidateArray = context_ptr->fast_candidate_array;
    MacroBlockD  *xd = cu_ptr->av1xd;
    uint8_t drli, maxDrlIndex;
    IntMv nearestmv[2], nearmv[2], ref_mv[2];

    //NEAREST_L0
        candidateArray[canIdx].type = INTER_MODE;
        candidateArray[canIdx].inter_mode = NEARESTMV;
        candidateArray[canIdx].pred_mode = NEARESTMV;
        candidateArray[canIdx].motion_mode = WARPED_CAUSAL;
        candidateArray[canIdx].wm_params.wmtype = AFFINE;
        candidateArray[canIdx].is_compound = 0;
        candidateArray[canIdx].distortion_ready = 0;
        candidateArray[canIdx].use_intrabc = 0;
        candidateArray[canIdx].merge_flag = EB_FALSE;
        candidateArray[canIdx].prediction_direction[0] = UNI_PRED_LIST_0;
        candidateArray[canIdx].is_new_mv = 0;
        candidateArray[canIdx].is_zero_mv = 0;
        candidateArray[canIdx].motion_vector_xl0 = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.col;
        candidateArray[canIdx].motion_vector_yl0 = context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.row;
        candidateArray[canIdx].drl_index = 0;
        candidateArray[canIdx].ref_mv_index = 0;
        candidateArray[canIdx].pred_mv_weight = 0;
        candidateArray[canIdx].ref_frame_type = LAST_FRAME;
        candidateArray[canIdx].ref_frame_index_l0 = 0;
        candidateArray[canIdx].ref_frame_index_l1 = -1;
        candidateArray[canIdx].transform_type[0] = DCT_DCT;
        candidateArray[canIdx].transform_type_uv = DCT_DCT;

        Mv mv_0;
        mv_0.x = candidateArray[canIdx].motion_vector_xl0;
        mv_0.y = candidateArray[canIdx].motion_vector_yl0;
        MvUnit mv_unit;
        mv_unit.mv[0] = mv_0;
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
            INCRMENT_CAND_TOTAL_COUNT(canIdx);

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
            candidateArray[canIdx].type = INTER_MODE;
            candidateArray[canIdx].inter_mode = NEARMV;
            candidateArray[canIdx].pred_mode = NEARMV;
            candidateArray[canIdx].motion_mode = WARPED_CAUSAL;
            candidateArray[canIdx].wm_params.wmtype = AFFINE;
            candidateArray[canIdx].is_compound = 0;
            candidateArray[canIdx].distortion_ready = 0;
            candidateArray[canIdx].use_intrabc = 0;
            candidateArray[canIdx].merge_flag = EB_FALSE;
            candidateArray[canIdx].prediction_direction[0] = UNI_PRED_LIST_0;
            candidateArray[canIdx].is_new_mv = 0;
            candidateArray[canIdx].is_zero_mv = 0;
            candidateArray[canIdx].motion_vector_xl0 = nearmv[0].as_mv.col;
            candidateArray[canIdx].motion_vector_yl0 = nearmv[0].as_mv.row;
            candidateArray[canIdx].drl_index = drli;
            candidateArray[canIdx].ref_mv_index = 0;
            candidateArray[canIdx].pred_mv_weight = 0;
            candidateArray[canIdx].ref_frame_type = LAST_FRAME;
            candidateArray[canIdx].ref_frame_index_l0 = 0;
            candidateArray[canIdx].ref_frame_index_l1 = -1;
            candidateArray[canIdx].transform_type[0] = DCT_DCT;
            candidateArray[canIdx].transform_type_uv = DCT_DCT;

            Mv mv_0;
            mv_0.x = candidateArray[canIdx].motion_vector_xl0;
            mv_0.y = candidateArray[canIdx].motion_vector_yl0;
            MvUnit mv_unit;
            mv_unit.mv[0] = mv_0;
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
                INCRMENT_CAND_TOTAL_COUNT(canIdx);
    }

    // NEWMV L0
    const MV neighbors[9] = { { 0, 0 },
        { 0, -2 }, { 2, 0 }, { 0, 2 }, { -2, 0 } ,
        { 2, -2 }, { 2, 2 }, { 2, 2 }, { -2, 2 } };

    IntMv  bestPredmv[2] = { {0}, {0} };

    uint8_t total_me_cnt = meResult->total_me_candidate_index[context_ptr->me_block_offset];
    const MeCandidate *me_block_results = meResult->me_candidate[context_ptr->me_block_offset];
    //const MeLcuResults_t *meResults = pictureControlSetPtr->ParentPcsPtr->meResultsPtr[lcuAddr];
    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
    {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
        if (inter_direction == 0) {
    for (int i=0; i<9; i++){

        candidateArray[canIdx].type = INTER_MODE;
        candidateArray[canIdx].distortion_ready = 0;
        candidateArray[canIdx].use_intrabc = 0;
        candidateArray[canIdx].merge_flag = EB_FALSE;
        candidateArray[canIdx].prediction_direction[0] = (EbPredDirection)0;
        candidateArray[canIdx].inter_mode = NEWMV;
        candidateArray[canIdx].pred_mode = NEWMV;
        candidateArray[canIdx].motion_mode = WARPED_CAUSAL;
        candidateArray[canIdx].wm_params.wmtype = AFFINE;

        candidateArray[canIdx].is_compound = 0;
        candidateArray[canIdx].is_new_mv = 1;
        candidateArray[canIdx].is_zero_mv = 0;

        candidateArray[canIdx].drl_index = 0;

        // Set the MV to ME result
        candidateArray[canIdx].motion_vector_xl0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : meResult->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1;
        candidateArray[canIdx].motion_vector_yl0 = use_close_loop_me ? ss_mecontext->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : meResult->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1;
        candidateArray[canIdx].motion_vector_xl0 += neighbors[i].col;
        candidateArray[canIdx].motion_vector_yl0 += neighbors[i].row;
        candidateArray[canIdx].ref_mv_index = 0;
        candidateArray[canIdx].pred_mv_weight = 0;
        candidateArray[canIdx].ref_frame_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
        candidateArray[canIdx].ref_frame_index_l0 = list0_ref_index;
        candidateArray[canIdx].ref_frame_index_l1 = -1;

        candidateArray[canIdx].transform_type[0] = DCT_DCT;
        candidateArray[canIdx].transform_type_uv = DCT_DCT;

        ChooseBestAv1MvPred(
            context_ptr,
            candidateArray[canIdx].md_rate_estimation_ptr,
            context_ptr->cu_ptr,
            candidateArray[canIdx].ref_frame_type,
            candidateArray[canIdx].is_compound,
            candidateArray[canIdx].pred_mode,
            candidateArray[canIdx].motion_vector_xl0,
            candidateArray[canIdx].motion_vector_yl0,
            0, 0,
            &candidateArray[canIdx].drl_index,
            bestPredmv);

        candidateArray[canIdx].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
        candidateArray[canIdx].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;

        Mv mv_0;
        mv_0.x = candidateArray[canIdx].motion_vector_xl0;
        mv_0.y = candidateArray[canIdx].motion_vector_yl0;
        MvUnit mv_unit;
        mv_unit.mv[0] = mv_0;
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
            INCRMENT_CAND_TOTAL_COUNT(canIdx);
    }
        }
    }

    *candTotCnt = canIdx;
}




void inject_new_candidates(
    const SequenceControlSet   *sequence_control_set_ptr,
    struct ModeDecisionContext *context_ptr,
    PictureControlSet          *picture_control_set_ptr,
    EbBool                      isCompoundEnabled,
    EbBool                      allow_bipred,
    uint32_t                    me_sb_addr,
    SsMeContext                *inloop_me_context,
    EbBool                      use_close_loop_me,
    uint32_t                    close_loop_me_index,
    uint32_t                    me_block_offset,
    uint32_t                   *candidateTotalCnt) {

    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    IntMv  bestPredmv[2] = { {0}, {0} };
    uint32_t canTotalCnt = (*candidateTotalCnt);

    const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[me_block_offset];
    const MeCandidate *me_block_results = me_results->me_candidate[me_block_offset];

    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
    {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
        const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;

        /**************
            NEWMV L0
        ************* */
        if (inter_direction == 0) {

            int16_t to_inject_mv_x = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : me_results->me_mv_array[me_block_offset][list0_ref_index].x_mv << 1;
            int16_t to_inject_mv_y = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : me_results->me_mv_array[me_block_offset][list0_ref_index].y_mv << 1;
            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
            if (context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {


                candidateArray[canTotalCnt].type = INTER_MODE;
                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;
                candidateArray[canTotalCnt].inter_mode = NEWMV;
                candidateArray[canTotalCnt].pred_mode = NEWMV;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

                candidateArray[canTotalCnt].is_compound = 0;
                candidateArray[canTotalCnt].is_new_mv = 1;
                candidateArray[canTotalCnt].is_zero_mv = 0;

                candidateArray[canTotalCnt].drl_index = 0;

                // Set the MV to ME result
                candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x;
                candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y;

                // will be needed later by the rate estimation
                candidateArray[canTotalCnt].ref_mv_index = 0;
                candidateArray[canTotalCnt].pred_mv_weight = 0;
                candidateArray[canTotalCnt].ref_frame_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
                candidateArray[canTotalCnt].ref_frame_index_l0 = list0_ref_index;
                candidateArray[canTotalCnt].ref_frame_index_l1 = -1;

                candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;

                ChooseBestAv1MvPred(
                    context_ptr,
                    candidateArray[canTotalCnt].md_rate_estimation_ptr,
                    context_ptr->cu_ptr,
                    candidateArray[canTotalCnt].ref_frame_type,
                    candidateArray[canTotalCnt].is_compound,
                    candidateArray[canTotalCnt].pred_mode,
                    candidateArray[canTotalCnt].motion_vector_xl0,
                    candidateArray[canTotalCnt].motion_vector_yl0,
                    0, 0,
                    &candidateArray[canTotalCnt].drl_index,
                    bestPredmv);

                candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;

                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

                context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
                context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
                context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = to_inject_ref_type;
                ++context_ptr->injected_mv_count_l0;
            }

            }

        if (isCompoundEnabled) {
            /**************
               NEWMV L1
           ************* */
            if (inter_direction == 1) {
                int16_t to_inject_mv_x = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : me_results->me_mv_array[me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].x_mv << 1;
                int16_t to_inject_mv_y = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : me_results->me_mv_array[me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].y_mv << 1;
                uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                if (context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {

                    candidateArray[canTotalCnt].type = INTER_MODE;
                    candidateArray[canTotalCnt].distortion_ready = 0;
                    candidateArray[canTotalCnt].use_intrabc = 0;
                    candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                    candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)1;

                    candidateArray[canTotalCnt].inter_mode = NEWMV;
                    candidateArray[canTotalCnt].pred_mode = NEWMV;
                    candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

                    candidateArray[canTotalCnt].is_compound = 0;
                    candidateArray[canTotalCnt].is_new_mv = 1;
                    candidateArray[canTotalCnt].is_zero_mv = 0;

                    candidateArray[canTotalCnt].drl_index = 0;

                    // Set the MV to ME result
                    candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x;
                    candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y;

                    // will be needed later by the rate estimation
                    candidateArray[canTotalCnt].ref_mv_index = 0;
                    candidateArray[canTotalCnt].pred_mv_weight = 0;
                    candidateArray[canTotalCnt].ref_frame_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                    candidateArray[canTotalCnt].ref_frame_index_l0 = -1;
                    candidateArray[canTotalCnt].ref_frame_index_l1 = list1_ref_index;

                    candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                    candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;

                    ChooseBestAv1MvPred(
                        context_ptr,
                        candidateArray[canTotalCnt].md_rate_estimation_ptr,
                        context_ptr->cu_ptr,
                        candidateArray[canTotalCnt].ref_frame_type,
                        candidateArray[canTotalCnt].is_compound,
                        candidateArray[canTotalCnt].pred_mode,
                        candidateArray[canTotalCnt].motion_vector_xl1,
                        candidateArray[canTotalCnt].motion_vector_yl1,
                        0, 0,
                        &candidateArray[canTotalCnt].drl_index,
                        bestPredmv);

                    candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[0].as_mv.col;
                    candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[0].as_mv.row;
                    INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

                    context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
                    context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
                    context_ptr->injected_ref_type_l1_array[context_ptr->injected_mv_count_l1] = to_inject_ref_type;
                    ++context_ptr->injected_mv_count_l1;
                }

                }
            /**************
               NEW_NEWMV
            ************* */
            if (allow_bipred) {

                if (inter_direction == 2) {
                    int16_t to_inject_mv_x_l0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][0] << 1 : me_results->me_mv_array[me_block_offset][list0_ref_index].x_mv << 1;
                    int16_t to_inject_mv_y_l0 = use_close_loop_me ? inloop_me_context->inloop_me_mv[0][0][close_loop_me_index][1] << 1 : me_results->me_mv_array[me_block_offset][list0_ref_index].y_mv << 1;
                    int16_t to_inject_mv_x_l1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][0] << 1 : me_results->me_mv_array[me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv << 1;
                    int16_t to_inject_mv_y_l1 = use_close_loop_me ? inloop_me_context->inloop_me_mv[1][0][close_loop_me_index][1] << 1 : me_results->me_mv_array[me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv << 1;
                    MvReferenceFrame rf[2];
                    rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
                    rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
                    uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
                    if (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE) {

                        candidateArray[canTotalCnt].type = INTER_MODE;

                        candidateArray[canTotalCnt].distortion_ready = 0;
                        candidateArray[canTotalCnt].use_intrabc = 0;

                        candidateArray[canTotalCnt].merge_flag = EB_FALSE;

                        candidateArray[canTotalCnt].is_new_mv = 1;
                        candidateArray[canTotalCnt].is_zero_mv = 0;

                        candidateArray[canTotalCnt].drl_index = 0;

                        // Set the MV to ME result

                        candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x_l0;
                        candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y_l0;
                        candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x_l1;
                        candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y_l1;

                        // will be needed later by the rate estimation
                        candidateArray[canTotalCnt].ref_mv_index = 0;
                        candidateArray[canTotalCnt].pred_mv_weight = 0;

                        candidateArray[canTotalCnt].inter_mode = NEW_NEWMV;
                        candidateArray[canTotalCnt].pred_mode = NEW_NEWMV;
                        candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                        candidateArray[canTotalCnt].is_compound = 1;
                        candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;
                        MvReferenceFrame rf[2];
                        rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
                        rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
                        candidateArray[canTotalCnt].ref_frame_type = av1_ref_frame_type(rf);

                        candidateArray[canTotalCnt].ref_frame_index_l0 = list0_ref_index;
                        candidateArray[canTotalCnt].ref_frame_index_l1 = list1_ref_index;
                        candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                        candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
                        ChooseBestAv1MvPred(
                            context_ptr,
                            candidateArray[canTotalCnt].md_rate_estimation_ptr,
                            context_ptr->cu_ptr,
                            candidateArray[canTotalCnt].ref_frame_type,
                            candidateArray[canTotalCnt].is_compound,
                            candidateArray[canTotalCnt].pred_mode,
                            candidateArray[canTotalCnt].motion_vector_xl0,
                            candidateArray[canTotalCnt].motion_vector_yl0,
                            candidateArray[canTotalCnt].motion_vector_xl1,
                            candidateArray[canTotalCnt].motion_vector_yl1,
                            &candidateArray[canTotalCnt].drl_index,
                            bestPredmv);

                        candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
                        candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
                        candidateArray[canTotalCnt].motion_vector_pred_x[REF_LIST_1] = bestPredmv[1].as_mv.col;
                        candidateArray[canTotalCnt].motion_vector_pred_y[REF_LIST_1] = bestPredmv[1].as_mv.row;
                        INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

                        context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                        context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                        context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                        context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                        context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = to_inject_ref_type;
                        ++context_ptr->injected_mv_count_bipred;
                    }

                    }
                }
            }
            }
    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;
        }

void  inject_inter_candidates(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    SsMeContext                  *ss_mecontext,
    const SequenceControlSet     *sequence_control_set_ptr,
    LargestCodingUnit            *sb_ptr,
    uint32_t                       *candidateTotalCnt) {

    (void)sequence_control_set_ptr;
    uint32_t                   canTotalCnt = *candidateTotalCnt;
    const uint32_t             lcuAddr = sb_ptr->index;
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    uint32_t geom_offset_x = 0;
    uint32_t geom_offset_y = 0;

    if (sequence_control_set_ptr->seq_header.sb_size == BLOCK_128X128) {
        uint32_t me_sb_size = sequence_control_set_ptr->sb_sz;
        uint32_t me_pic_width_in_sb = (sequence_control_set_ptr->seq_header.max_frame_width + sequence_control_set_ptr->sb_sz - 1) / me_sb_size;
        uint32_t me_sb_x = (context_ptr->cu_origin_x / me_sb_size);
        uint32_t me_sb_y = (context_ptr->cu_origin_y / me_sb_size);
        context_ptr->me_sb_addr = me_sb_x + me_sb_y * me_pic_width_in_sb;
        geom_offset_x = (me_sb_x & 0x1) * me_sb_size;
        geom_offset_y = (me_sb_y & 0x1) * me_sb_size;
    }
    else
        context_ptr->me_sb_addr = lcuAddr;
    uint32_t max_number_of_pus_per_sb;

    max_number_of_pus_per_sb = picture_control_set_ptr->parent_pcs_ptr->max_number_of_pus_per_sb;
    context_ptr->me_block_offset =
        (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4 || context_ptr->blk_geom->bwidth == 128 || context_ptr->blk_geom->bheight == 128) ?
            0 :
            get_me_info_index(max_number_of_pus_per_sb, context_ptr->blk_geom, geom_offset_x, geom_offset_y);

    MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[context_ptr->me_sb_addr];

    EbBool use_close_loop_me = picture_control_set_ptr->parent_pcs_ptr->enable_in_loop_motion_estimation_flag &&
        ((context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) || (context_ptr->blk_geom->bwidth > 64 || context_ptr->blk_geom->bheight > 64)) ? EB_TRUE : EB_FALSE;

    uint32_t close_loop_me_index = use_close_loop_me ? get_in_loop_me_info_index(MAX_SS_ME_PU_COUNT, sequence_control_set_ptr->seq_header.sb_size == BLOCK_128X128 ? 1 : 0, context_ptr->blk_geom) : 0;
    EbBool allow_bipred = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : EB_TRUE;
    uint8_t sq_index = LOG2F(context_ptr->blk_geom->sq_size) - 2;
    uint8_t inject_newmv_candidate = 1;

    if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level >= NSQ_SEARCH_LEVEL1 &&
        picture_control_set_ptr->parent_pcs_ptr->nsq_search_level < NSQ_SEARCH_FULL) {
        inject_newmv_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
            context_ptr->parent_sq_has_coeff[sq_index] != 0 ? inject_newmv_candidate : 0;
    }

    generate_av1_mvp_table(
        &sb_ptr->tile_info,
        context_ptr,
        context_ptr->cu_ptr,
        context_ptr->blk_geom,
        context_ptr->cu_origin_x,
        context_ptr->cu_origin_y,
        picture_control_set_ptr->parent_pcs_ptr->ref_frame_type_arr,
        picture_control_set_ptr->parent_pcs_ptr->tot_ref_frame_types,
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

    uint32_t refIt;
    //all of ref pairs: (1)single-ref List0  (2)single-ref List1  (3)compound Bi-Dir List0-List1  (4)compound Uni-Dir List0-List0  (5)compound Uni-Dir List1-List1
    for (refIt = 0; refIt < picture_control_set_ptr->parent_pcs_ptr->tot_ref_frame_types; ++refIt) {
        MvReferenceFrame ref_frame_pair = picture_control_set_ptr->parent_pcs_ptr->ref_frame_type_arr[refIt];
        inject_mvp_candidates_II(
            context_ptr,
            picture_control_set_ptr,
            context_ptr->cu_ptr,
            ref_frame_pair,
            &canTotalCnt);
    }

    //----------------------
    //    NEAREST_NEWMV, NEW_NEARESTMV, NEAR_NEWMV, NEW_NEARMV.
    //----------------------
    if (context_ptr->new_nearest_near_comb_injection) {
        EbBool allow_compound = (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE || context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : EB_TRUE;
        if (allow_compound) {
            //all of ref pairs: (1)single-ref List0  (2)single-ref List1  (3)compound Bi-Dir List0-List1  (4)compound Uni-Dir List0-List0  (5)compound Uni-Dir List1-List1
            for (refIt = 0; refIt < picture_control_set_ptr->parent_pcs_ptr->tot_ref_frame_types; ++refIt) {
                MvReferenceFrame ref_frame_pair = picture_control_set_ptr->parent_pcs_ptr->ref_frame_type_arr[refIt];
                inject_new_nearest_new_comb_candidates(
                    sequence_control_set_ptr,
                    context_ptr,
                    picture_control_set_ptr,
                    ref_frame_pair,
                    &canTotalCnt);
            }
        }
    }

    if (inject_newmv_candidate) {

        inject_new_candidates(
            sequence_control_set_ptr,
            context_ptr,
            picture_control_set_ptr,
            isCompoundEnabled,
            allow_bipred,
            context_ptr->me_sb_addr,
            ss_mecontext,
            use_close_loop_me,
            close_loop_me_index,
            context_ptr->me_block_offset,
            &canTotalCnt);

        if (context_ptr->nx4_4xn_parent_mv_injection) {
            // If Nx4 or 4xN the inject the MV of the aprent block


            // Derive whether if current block would need to have offsets made
            uint32_t bwidth_offset_to_8 = (context_ptr->blk_geom->bwidth == 4) << 2;
            uint32_t bheight_offset_to_8 = (context_ptr->blk_geom->bheight == 4) << 2;

            // if there is an offset needed to set either dimension to 8
            if (bwidth_offset_to_8 || bheight_offset_to_8) {

                // Align parent block has dimensions inherited by current block, if current block has a dimension of 4
                // add 4 so the resulting block follows an 8x8 basis
                uint32_t bwidth_to_search = context_ptr->blk_geom->bwidth + bwidth_offset_to_8;
                uint32_t bheight_to_search = context_ptr->blk_geom->bheight + bheight_offset_to_8;

                // Align parent block has origin inherited by current block
                uint32_t x_to_search = context_ptr->blk_geom->origin_x - (geom_offset_x + ((context_ptr->blk_geom->origin_x & 0x7) ? 4 : 0));
                uint32_t y_to_search = context_ptr->blk_geom->origin_y - (geom_offset_y + ((context_ptr->blk_geom->origin_y & 0x7) ? 4 : 0));

                // Search the me_info_index of the parent block
                uint32_t me_info_index = 0;
                for (uint32_t block_index = 0; block_index < max_number_of_pus_per_sb; block_index++) {

                    if (
                        (bwidth_to_search == partition_width[block_index]) &&
                        (bheight_to_search == partition_height[block_index]) &&
                        (x_to_search == pu_search_index_map[block_index][0]) &&
                        (y_to_search == pu_search_index_map[block_index][1]))
                    {
                        me_info_index = block_index;
                        break;
                    }
                }

                inject_new_candidates(
                    sequence_control_set_ptr,
                    context_ptr,
                    picture_control_set_ptr,
                    isCompoundEnabled,
                    allow_bipred,
                    context_ptr->me_sb_addr,
                    ss_mecontext,
                    use_close_loop_me,
                    close_loop_me_index,
                    me_info_index,
                    &canTotalCnt);
            }
        }
    }

    if (context_ptr->global_mv_injection) {
        /**************
         GLOBALMV L0
        ************* */
        {
            int16_t to_inject_mv_x = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_y = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, 0/*list0_ref_index*/);
            if (context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {

                candidateArray[canTotalCnt].type = INTER_MODE;

                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;

                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;

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
                candidateArray[canTotalCnt].ref_frame_index_l0 = 0;
                candidateArray[canTotalCnt].ref_frame_index_l1 = -1;

                candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
                // Set the MV to frame MV
                candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x;
                candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y;

                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

                context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
                context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
                context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = to_inject_ref_type;
                ++context_ptr->injected_mv_count_l0;
            }
            }

        if (isCompoundEnabled && allow_bipred) {
            /**************
            GLOBAL_GLOBALMV
            ************* */

            int16_t to_inject_mv_x_l0 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_y_l0 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_x_l1 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_y_l1 = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[BWDREF_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
            MvReferenceFrame rf[2];
            rf[0] = svt_get_ref_frame_type(REF_LIST_0, 0/*list0_ref_index*/);
            rf[1] = svt_get_ref_frame_type(REF_LIST_1, 0/*list1_ref_index*/);
            uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
            if (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE) {

                candidateArray[canTotalCnt].type = INTER_MODE;
                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;

                candidateArray[canTotalCnt].merge_flag = EB_FALSE;

                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;

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
                candidateArray[canTotalCnt].ref_frame_index_l0 = 0;
                candidateArray[canTotalCnt].ref_frame_index_l1 = 0;
                candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
                // Set the MV to frame MV

                candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x_l0;
                candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y_l0;
                candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x_l1;
                candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y_l1;

                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

                context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = to_inject_ref_type;
                ++context_ptr->injected_mv_count_bipred;
            }
            }
        }

    // Warped Motion
    if (picture_control_set_ptr->parent_pcs_ptr->allow_warped_motion &&
        has_overlappable_candidates(context_ptr->cu_ptr) &&
        context_ptr->blk_geom->bwidth >= 8 &&
        context_ptr->blk_geom->bheight >= 8 &&
        context_ptr->warped_motion_injection) {
        inject_warped_motion_candidates(
            picture_control_set_ptr,
            context_ptr,
            context_ptr->cu_ptr,
            &canTotalCnt,
            ss_mecontext,
            me_results,
            use_close_loop_me,
            close_loop_me_index);
    }

    if (inject_newmv_candidate) {
        if (isCompoundEnabled) {
            if (allow_bipred) {

            //----------------------
            // Bipred2Nx2N
            //----------------------
            if (context_ptr->bipred3x3_injection > 0)
                if (picture_control_set_ptr->slice_type == B_SLICE)
                    Bipred3x3CandidatesInjection(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        context_ptr,
                        sb_ptr,
                        context_ptr->me_sb_addr,
                        ss_mecontext,
                        use_close_loop_me,
                        close_loop_me_index,
                        &canTotalCnt);

        }

        //----------------------
        // Unipred2Nx2N
        //----------------------
        if (context_ptr->unipred3x3_injection > 0)
            if (picture_control_set_ptr->slice_type != I_SLICE)
                Unipred3x3CandidatesInjection(
                    sequence_control_set_ptr,
                    picture_control_set_ptr,
                    context_ptr,
                    sb_ptr,
                    context_ptr->me_sb_addr,
                    ss_mecontext,
                    use_close_loop_me,
                    close_loop_me_index,
                    &canTotalCnt);
            }
        }

#if EIGTH_PEL_MV
        //----------------------
        // Eighth-pel refinement
        //----------------------
        if (inject_newmv_candidate && picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv) {
            if (allow_bipred) {
                //----------------------
                // Inject eight-pel bi-pred
                //----------------------
                if (context_ptr->bipred3x3_injection > 0)
                    if (picture_control_set_ptr->slice_type == B_SLICE)
                        eighth_pel_bipred_refinement(
                            sequence_control_set_ptr,
                            picture_control_set_ptr,
                            context_ptr,
                            me_sb_addr,
                            ss_mecontext,
                            use_close_loop_me,
                            close_loop_me_index,
                            &canTotalCnt);

            //----------------------
            // Inject eight-pel uni-pred
            //----------------------
            if (context_ptr->unipred3x3_injection > 0)
                if (picture_control_set_ptr->slice_type != I_SLICE)
                    eighth_pel_unipred_refinement(
                        sequence_control_set_ptr,
                        picture_control_set_ptr,
                        context_ptr,
                        me_sb_addr,
                        ss_mecontext,
                        use_close_loop_me,
                        close_loop_me_index,
                        &canTotalCnt);
                }
            }
#endif

// update the total number of candidates injected
(*candidateTotalCnt) = canTotalCnt;

return;
    }

 extern PredictionMode get_uv_mode(UvPredictionMode mode) {
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
    PlaneType plane_type) {
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
    BlockSize  sb_type,
    int32_t   is_inter,
    PredictionMode pred_mode,
    UvPredictionMode pred_mode_uv,
    PlaneType plane_type,
    const MacroBlockD *xd, int32_t blk_row,
    int32_t blk_col, TxSize tx_size,
    int32_t reduced_tx_set)
{
    UNUSED(sb_type);
    UNUSED(*xd);
    UNUSED(blk_row);
    UNUSED(blk_col);

    // BlockSize  sb_type = BLOCK_8X8;

    MbModeInfo  mbmi;
    mbmi.mode = pred_mode;
    mbmi.uv_mode = pred_mode_uv;

    // const MbModeInfo *const mbmi = xd->mi[0];
    // const struct MacroblockdPlane *const pd = &xd->plane[plane_type];
    const TxSetType tx_set_type =
        /*av1_*/get_ext_tx_set_type(tx_size, is_inter, reduced_tx_set);

    TxType tx_type = DCT_DCT;
    if ( /*xd->lossless[mbmi->segment_id] ||*/ txsize_sqr_up_map[tx_size] > TX_32X32)
        tx_type = DCT_DCT;
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

void  inject_intra_candidates_ois(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    LargestCodingUnit            *sb_ptr,
    uint32_t                       *candidate_total_cnt){
    uint8_t                     intra_candidate_counter;
    uint8_t                     intra_mode;
    uint32_t                    can_total_cnt = 0;
    ModeDecisionCandidate    *candidate_array = context_ptr->fast_candidate_array;
    EbBool                      disable_cfl_flag = (MAX(context_ptr->blk_geom->bheight, context_ptr->blk_geom->bwidth) > 32) ? EB_TRUE : EB_FALSE;

    OisSbResults    *ois_sb_results_ptr = picture_control_set_ptr->parent_pcs_ptr->ois_sb_results[sb_ptr->index];
    OisCandidate     *ois_blk_ptr = ois_sb_results_ptr->ois_candidate_array[ep_to_pa_block_index[context_ptr->blk_geom->blkidx_mds]];
    uint8_t              total_intra_luma_mode = ois_sb_results_ptr-> total_ois_intra_candidate[ep_to_pa_block_index[context_ptr->blk_geom->blkidx_mds]];

    for (intra_candidate_counter = 0; intra_candidate_counter < total_intra_luma_mode; ++intra_candidate_counter) {
        intra_mode = ois_blk_ptr[can_total_cnt].intra_mode;
        assert(intra_mode < INTRA_MODES);
        if (av1_is_directional_mode((PredictionMode)intra_mode)) {
            int32_t angle_delta = ois_blk_ptr[can_total_cnt].angle_delta ;
            candidate_array[can_total_cnt].type = INTRA_MODE;
            candidate_array[can_total_cnt].intra_luma_mode = intra_mode;
            candidate_array[can_total_cnt].distortion_ready =  1;
            candidate_array[can_total_cnt].me_distortion = ois_blk_ptr[can_total_cnt].distortion;
            candidate_array[can_total_cnt].use_intrabc = 0;
            candidate_array[can_total_cnt].is_directional_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)intra_mode);
            candidate_array[can_total_cnt].angle_delta[PLANE_TYPE_Y] = angle_delta;
            candidate_array[can_total_cnt].intra_chroma_mode = disable_cfl_flag ? intra_luma_to_chroma[intra_mode] :
                                                               context_ptr->chroma_level <= CHROMA_MODE_1 ? UV_CFL_PRED : UV_DC_PRED;

            candidate_array[can_total_cnt].cfl_alpha_signs = 0;
            candidate_array[can_total_cnt].cfl_alpha_idx = 0;
            candidate_array[can_total_cnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidate_array[can_total_cnt].intra_chroma_mode);
            candidate_array[can_total_cnt].angle_delta[PLANE_TYPE_UV] = 0;

            candidate_array[can_total_cnt].transform_type[0] = DCT_DCT;

            if (candidate_array[can_total_cnt].intra_chroma_mode == UV_CFL_PRED)
                candidate_array[can_total_cnt].transform_type_uv = DCT_DCT;
            else
                candidate_array[can_total_cnt].transform_type_uv =
                av1_get_tx_type(
                    context_ptr->blk_geom->bsize,
                    0,
                    (PredictionMode)candidate_array[can_total_cnt].intra_luma_mode,
                    (UvPredictionMode)candidate_array[can_total_cnt].intra_chroma_mode,
                    PLANE_TYPE_UV,
                    0,
                    0,
                    0,
                    context_ptr->blk_geom->txsize_uv[0][0],
                    picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
            candidate_array[can_total_cnt].ref_frame_type = INTRA_FRAME;
            candidate_array[can_total_cnt].pred_mode = (PredictionMode)intra_mode;
            candidate_array[can_total_cnt].motion_mode = SIMPLE_TRANSLATION;
            INCRMENT_CAND_TOTAL_COUNT(can_total_cnt);
        }
        else {
            candidate_array[can_total_cnt].type = INTRA_MODE;
            candidate_array[can_total_cnt].intra_luma_mode = intra_mode;
            candidate_array[can_total_cnt].distortion_ready =  1;
            candidate_array[can_total_cnt].me_distortion = ois_blk_ptr[can_total_cnt].distortion;
            candidate_array[can_total_cnt].use_intrabc = 0;
            candidate_array[can_total_cnt].is_directional_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)intra_mode);
            candidate_array[can_total_cnt].angle_delta[PLANE_TYPE_Y] = 0;
            candidate_array[can_total_cnt].intra_chroma_mode =  disable_cfl_flag ? intra_luma_to_chroma[intra_mode] :
                                                                context_ptr->chroma_level <= CHROMA_MODE_1 ? UV_CFL_PRED : UV_DC_PRED;

            candidate_array[can_total_cnt].cfl_alpha_signs = 0;
            candidate_array[can_total_cnt].cfl_alpha_idx = 0;
            candidate_array[can_total_cnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidate_array[can_total_cnt].intra_chroma_mode);
            candidate_array[can_total_cnt].angle_delta[PLANE_TYPE_UV] = 0;
            candidate_array[can_total_cnt].transform_type[0] = DCT_DCT;

            if (candidate_array[can_total_cnt].intra_chroma_mode == UV_CFL_PRED)
                candidate_array[can_total_cnt].transform_type_uv = DCT_DCT;
            else
                candidate_array[can_total_cnt].transform_type_uv =
                av1_get_tx_type(
                    context_ptr->blk_geom->bsize,
                    0,
                    (PredictionMode)candidate_array[can_total_cnt].intra_luma_mode,
                    (UvPredictionMode)candidate_array[can_total_cnt].intra_chroma_mode,
                    PLANE_TYPE_UV,
                    0,
                    0,
                    0,
                    context_ptr->blk_geom->txsize_uv[0][0],
                    picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
            candidate_array[can_total_cnt].ref_frame_type = INTRA_FRAME;
            candidate_array[can_total_cnt].pred_mode = (PredictionMode)intra_mode;
            candidate_array[can_total_cnt].motion_mode = SIMPLE_TRANSLATION;
            INCRMENT_CAND_TOTAL_COUNT(can_total_cnt);
        }
    }

    // update the total number of candidates injected
    (*candidate_total_cnt) = can_total_cnt;

    return;
}

double av1_convert_qindex_to_q(int32_t qindex, AomBitDepth bit_depth);

static INLINE void setup_pred_plane(struct Buf2D *dst, BlockSize bsize,
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
void av1_setup_pred_block(BlockSize sb_type,
    struct Buf2D dst[MAX_MB_PLANE],
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
    AomBitDepth bit_depth) {
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
    PictureControlSet            *pcs,
    ModeDecisionContext          *context_ptr,
    const SequenceControlSet     *scs,
    LargestCodingUnit            *sb_ptr,
    CodingUnit                   *cu_ptr,
    MV                             *dv_cand,
    uint8_t                        *num_dv_cand)
{
    IntraBcContext  x_st;
    IntraBcContext  *x = &x_st;
    //fill x with what needed.
    x->is_exhaustive_allowed =  context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4 ? 1 : 0;
    //CHKN crc calculator could be moved to mdContext and these init at init time.
    av1_crc_calculator_init(&x->crc_calculator1, 24, 0x5D6DCB);
    av1_crc_calculator_init(&x->crc_calculator2, 24, 0x864CFB);

    x->xd = cu_ptr->av1xd;
    x->nmv_vec_cost = context_ptr->md_rate_estimation_ptr->nmv_vec_cost;
    x->mv_cost_stack = context_ptr->md_rate_estimation_ptr->nmvcoststack;
    BlockSize bsize = context_ptr->blk_geom->bsize;
    assert(bsize < BlockSizeS_ALL);
    const Av1Common *const cm = pcs->parent_pcs_ptr->av1_cm;
    MvReferenceFrame ref_frame = INTRA_FRAME;
    generate_av1_mvp_table(
        &sb_ptr->tile_info,
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
    const int sb_row = mi_row >> scs->seq_header.sb_size_log2;
    const int sb_col = mi_col >> scs->seq_header.sb_size_log2;

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
    if (nearestmv.as_int == INVALID_MV)
        nearestmv.as_int = 0;
    if (nearmv.as_int == INVALID_MV)
        nearmv.as_int = 0;
    IntMv dv_ref = nearestmv.as_int == 0 ? nearmv : nearestmv;
    if (dv_ref.as_int == 0)
        av1_find_ref_dv(&dv_ref, tile, scs->seq_header.sb_mi_size, mi_row, mi_col);
    // Ref DV should not have sub-pel.
    assert((dv_ref.as_mv.col & 7) == 0);
    assert((dv_ref.as_mv.row & 7) == 0);
    context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[INTRA_FRAME][0].this_mv = dv_ref;

    /* pointer to current frame */
    Yv12BufferConfig cur_buf;
    link_Eb_to_aom_buffer_desc_8bit(
        pcs->parent_pcs_ptr->enhanced_picture_ptr,
        &cur_buf);
    struct Buf2D yv12_mb[MAX_MB_PLANE];
    av1_setup_pred_block(bsize, yv12_mb, &cur_buf, mi_row, mi_col);
    for (int i = 0; i < num_planes; ++i)
        x->xdplane[i].pre[0] = yv12_mb[i];  //ref in ME
    //setup src for DV search same as ref
    x->plane[0].src = x->xdplane[0].pre[0];

    enum IntrabcMotionDirection {
        IBC_MOTION_ABOVE,
        IBC_MOTION_LEFT,
        IBC_MOTION_DIRECTIONS
    };

    //up to two dv candidates will be generated
    enum IntrabcMotionDirection max_dir = pcs->parent_pcs_ptr->ibc_mode > 1 ? IBC_MOTION_LEFT : IBC_MOTION_DIRECTIONS;

    for (enum IntrabcMotionDirection dir = IBC_MOTION_ABOVE;
        dir < max_dir; ++dir) {
        const MvLimits tmp_mv_limits = x->mv_limits;

        switch (dir) {
        case IBC_MOTION_ABOVE:
            x->mv_limits.col_min = (tile->mi_col_start - mi_col) * MI_SIZE;
            x->mv_limits.col_max = (tile->mi_col_end - mi_col) * MI_SIZE - w;
            x->mv_limits.row_min = (tile->mi_row_start - mi_row) * MI_SIZE;
            x->mv_limits.row_max =
                (sb_row * scs->seq_header.sb_mi_size - mi_row) * MI_SIZE - h;
            break;
        case IBC_MOTION_LEFT:
            x->mv_limits.col_min = (tile->mi_col_start - mi_col) * MI_SIZE;
            x->mv_limits.col_max =
                (sb_col * scs->seq_header.sb_mi_size - mi_col) * MI_SIZE - w;
            // TODO: Minimize the overlap between above and
            // left areas.
            x->mv_limits.row_min = (tile->mi_row_start - mi_row) * MI_SIZE;
            int bottom_coded_mi_edge =
                AOMMIN((sb_row + 1) * scs->seq_header.sb_mi_size, tile->mi_row_end);
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
            scs->seq_header.sb_size_log2))
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
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    const SequenceControlSet     *sequence_control_set_ptr,
    LargestCodingUnit            *sb_ptr,
    CodingUnit                   *cu_ptr,
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

    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    uint32_t dv_i;

    for (dv_i = 0; dv_i < num_dv_cand; dv_i++)
    {
        candidateArray[*cand_cnt].type = INTRA_MODE;
        candidateArray[*cand_cnt].intra_luma_mode = DC_PRED;
        candidateArray[*cand_cnt].distortion_ready = 0;
        candidateArray[*cand_cnt].use_intrabc = 1;
        candidateArray[*cand_cnt].is_directional_mode_flag = 0;
        candidateArray[*cand_cnt].angle_delta[PLANE_TYPE_Y] = 0;
        candidateArray[*cand_cnt].intra_chroma_mode = UV_DC_PRED;
        candidateArray[*cand_cnt].cfl_alpha_signs = 0;
        candidateArray[*cand_cnt].cfl_alpha_idx = 0;
        candidateArray[*cand_cnt].is_directional_chroma_mode_flag = 0;
        candidateArray[*cand_cnt].angle_delta[PLANE_TYPE_UV] = 0;
        candidateArray[*cand_cnt].transform_type[0] = DCT_DCT;
        candidateArray[*cand_cnt].transform_type_uv = DCT_DCT;
        candidateArray[*cand_cnt].ref_frame_type = INTRA_FRAME;
        candidateArray[*cand_cnt].pred_mode = DC_PRED;
        candidateArray[*cand_cnt].motion_mode = SIMPLE_TRANSLATION;
        //inter ralated
        candidateArray[*cand_cnt].is_compound = 0;
        candidateArray[*cand_cnt].merge_flag = EB_FALSE;
        candidateArray[*cand_cnt].prediction_direction[0] = UNI_PRED_LIST_0;
        candidateArray[*cand_cnt].is_new_mv = 0;
        candidateArray[*cand_cnt].is_zero_mv = 0;
        candidateArray[*cand_cnt].motion_vector_xl0 = dv_cand[dv_i].col;
        candidateArray[*cand_cnt].motion_vector_yl0 = dv_cand[dv_i].row;
        candidateArray[*cand_cnt].motion_vector_pred_x[REF_LIST_0] = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[INTRA_FRAME][0].this_mv.as_mv.col;
        candidateArray[*cand_cnt].motion_vector_pred_y[REF_LIST_0] = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[INTRA_FRAME][0].this_mv.as_mv.row;
        candidateArray[*cand_cnt].drl_index = 0;
        candidateArray[*cand_cnt].ref_mv_index = 0;
        candidateArray[*cand_cnt].pred_mv_weight = 0;
        candidateArray[*cand_cnt].interp_filters = av1_broadcast_interp_filter(BILINEAR);
        INCRMENT_CAND_TOTAL_COUNT( (*cand_cnt) );
    }
}
// END of Function Declarations
void  inject_intra_candidates(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    const SequenceControlSet     *sequence_control_set_ptr,
    LargestCodingUnit            *sb_ptr,
    uint32_t                       *candidateTotalCnt){
    (void)sequence_control_set_ptr;
    (void)sb_ptr;
    uint8_t                     is16bit = (sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    uint8_t                     intra_mode_start = DC_PRED;
    uint8_t                     intra_mode_end   = is16bit ? SMOOTH_H_PRED : PAETH_PRED;
    uint8_t                     openLoopIntraCandidate;
    uint32_t                    canTotalCnt = 0;
    uint8_t                     angleDeltaCounter = 0;
    EbBool                      use_angle_delta = (context_ptr->blk_geom->bsize >= BLOCK_8X8);
    uint8_t                     angleDeltaCandidateCount = use_angle_delta ? 7 : 1;
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool                      disable_cfl_flag = (MAX(context_ptr->blk_geom->bheight, context_ptr->blk_geom->bwidth) > 32) ? EB_TRUE : EB_FALSE;

    uint8_t                     disable_z2_prediction;
    uint8_t                     disable_angle_refinement;
    uint8_t                     disable_angle_prediction;
        uint8_t     angle_delta_shift = 1;
    if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode == 4) {
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            intra_mode_end = is16bit ? SMOOTH_H_PRED : PAETH_PRED;
            angleDeltaCandidateCount = use_angle_delta ? 5 : 1;
            disable_angle_prediction = 0;
            angle_delta_shift = 2;
            disable_z2_prediction = 0;
        }
        else {
            intra_mode_end = DC_PRED;
            disable_angle_prediction = 1;
            angleDeltaCandidateCount = 1;
            angle_delta_shift = 1;
            disable_z2_prediction = 0;
        }
    }else
    if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode == 3){
        disable_z2_prediction       = 0;
        disable_angle_refinement    = 0;
        disable_angle_prediction    = 1;
        angleDeltaCandidateCount = disable_angle_refinement ? 1: angleDeltaCandidateCount;
    } else if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode == 2) {
        disable_z2_prediction       = 0;
        disable_angle_refinement    = 0 ;
        disable_angle_prediction    = (context_ptr->blk_geom->sq_size > 16 ||
                                       context_ptr->blk_geom->bwidth == 4 ||
                                       context_ptr->blk_geom->bheight == 4) ? 1 : 0;
        angleDeltaCandidateCount = disable_angle_refinement ? 1: angleDeltaCandidateCount;
    } else if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode == 1) {
        disable_z2_prediction       = (context_ptr->blk_geom->sq_size > 16 ||
                                       context_ptr->blk_geom->bwidth == 4 ||
                                       context_ptr->blk_geom->bheight == 4) ? 1 : 0;
        disable_angle_refinement    = (context_ptr->blk_geom->sq_size > 16 ||
                                       context_ptr->blk_geom->bwidth == 4 ||
                                       context_ptr->blk_geom->bheight == 4) ? 1 : 0;
        disable_angle_prediction    = 0;
        angleDeltaCandidateCount = disable_angle_refinement ? 1: angleDeltaCandidateCount;
    } else {
        disable_z2_prediction       = 0;
        disable_angle_refinement    = 0;
        disable_angle_prediction    = 0;
        angleDeltaCandidateCount = disable_angle_refinement ? 1: angleDeltaCandidateCount;
    }
#if MR_MODE
    disable_z2_prediction       = 0;
    disable_angle_refinement    = 0;
    disable_angle_prediction    = 0;
#endif
    for (openLoopIntraCandidate = intra_mode_start; openLoopIntraCandidate <= intra_mode_end ; ++openLoopIntraCandidate) {
        if (av1_is_directional_mode((PredictionMode)openLoopIntraCandidate)) {
            if (!disable_angle_prediction) {
                for (angleDeltaCounter = 0; angleDeltaCounter < angleDeltaCandidateCount; ++angleDeltaCounter) {
                    int32_t angle_delta = CLIP( angle_delta_shift * (angleDeltaCandidateCount == 1 ? 0 : angleDeltaCounter - (angleDeltaCandidateCount >> 1)), -3 , 3);
                    int32_t  p_angle = mode_to_angle_map[(PredictionMode)openLoopIntraCandidate] + angle_delta * ANGLE_STEP;
                    if (!disable_z2_prediction || (p_angle <= 90 || p_angle >= 180)) {
                        candidateArray[canTotalCnt].type = INTRA_MODE;
                        candidateArray[canTotalCnt].intra_luma_mode = openLoopIntraCandidate;
                        candidateArray[canTotalCnt].distortion_ready = 0;
                        candidateArray[canTotalCnt].use_intrabc = 0;
                        candidateArray[canTotalCnt].is_directional_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)openLoopIntraCandidate);
                        candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y] = angle_delta;
                        // Search the best independent intra chroma mode
                        if (context_ptr->chroma_level == CHROMA_MODE_0) {
                            candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ?
                                context_ptr->best_uv_mode[openLoopIntraCandidate][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]] :
                                UV_CFL_PRED ;
                            candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_UV] = disable_cfl_flag ?
                                context_ptr->best_uv_angle[candidateArray[canTotalCnt].intra_luma_mode][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]] : 0;
                            candidateArray[canTotalCnt].is_directional_chroma_mode_flag = disable_cfl_flag ?
                                (uint8_t)av1_is_directional_mode((PredictionMode)(context_ptr->best_uv_mode[candidateArray[canTotalCnt].intra_luma_mode][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]])) : 0;
                        }
                        else {
                            // Hsan/Omar: why the restriction below ? (i.e. disable_ang_uv)
                            const int32_t disable_ang_uv = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) && context_ptr->blk_geom->has_uv ? 1 : 0;
                            candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ?
                                intra_luma_to_chroma[openLoopIntraCandidate] :
                                (context_ptr->chroma_level == CHROMA_MODE_1) ?
                                UV_CFL_PRED :
                                UV_DC_PRED;
                            candidateArray[canTotalCnt].intra_chroma_mode = disable_ang_uv && av1_is_directional_mode(candidateArray[canTotalCnt].intra_chroma_mode) ?
                                UV_DC_PRED : candidateArray[canTotalCnt].intra_chroma_mode;
                            candidateArray[canTotalCnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidateArray[canTotalCnt].intra_chroma_mode);
                            candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_UV] = 0;
                        }
                        candidateArray[canTotalCnt].cfl_alpha_signs = 0;
                        candidateArray[canTotalCnt].cfl_alpha_idx = 0;
                        candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;

                        if (candidateArray[canTotalCnt].intra_chroma_mode == UV_CFL_PRED)
                            candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
                        else
                            candidateArray[canTotalCnt].transform_type_uv =
                            av1_get_tx_type(
                                context_ptr->blk_geom->bsize,
                                0,
                                (PredictionMode)candidateArray[canTotalCnt].intra_luma_mode,
                                (UvPredictionMode)candidateArray[canTotalCnt].intra_chroma_mode,
                                PLANE_TYPE_UV,
                                0,
                                0,
                                0,
                                context_ptr->blk_geom->txsize_uv[0][0],
                                picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
                        candidateArray[canTotalCnt].ref_frame_type = INTRA_FRAME;
                        candidateArray[canTotalCnt].pred_mode = (PredictionMode)openLoopIntraCandidate;
                        candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                        INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                    }
            }
        }
        }
        else {
            candidateArray[canTotalCnt].type = INTRA_MODE;
            candidateArray[canTotalCnt].intra_luma_mode = openLoopIntraCandidate;
            candidateArray[canTotalCnt].distortion_ready = 0;
            candidateArray[canTotalCnt].use_intrabc = 0;
            candidateArray[canTotalCnt].is_directional_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)openLoopIntraCandidate);
            candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y] = 0;
            // Search the best independent intra chroma mode
            if (context_ptr->chroma_level == CHROMA_MODE_0) {
                candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ?
                    context_ptr->best_uv_mode[openLoopIntraCandidate][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]] :
                    UV_CFL_PRED;
                candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_UV] = disable_cfl_flag ?
                    context_ptr->best_uv_angle[candidateArray[canTotalCnt].intra_luma_mode][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]] : 0;
                candidateArray[canTotalCnt].is_directional_chroma_mode_flag = disable_cfl_flag ?
                    (uint8_t)av1_is_directional_mode((PredictionMode)(context_ptr->best_uv_mode[candidateArray[canTotalCnt].intra_luma_mode][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]])) : 0;
            }
            else {
                // Hsan/Omar: why the restriction below ? (i.e. disable_ang_uv)
                const int32_t disable_ang_uv = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) && context_ptr->blk_geom->has_uv ? 1 : 0;
                candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ?
                    intra_luma_to_chroma[openLoopIntraCandidate] :
                    (context_ptr->chroma_level == CHROMA_MODE_1) ?
                        UV_CFL_PRED :
                        UV_DC_PRED;

                candidateArray[canTotalCnt].intra_chroma_mode = disable_ang_uv && av1_is_directional_mode(candidateArray[canTotalCnt].intra_chroma_mode) ?
                    UV_DC_PRED : candidateArray[canTotalCnt].intra_chroma_mode;

                candidateArray[canTotalCnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidateArray[canTotalCnt].intra_chroma_mode);
                candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_UV] = 0;

            }
            candidateArray[canTotalCnt].cfl_alpha_signs = 0;
            candidateArray[canTotalCnt].cfl_alpha_idx = 0;
            candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;

            if (candidateArray[canTotalCnt].intra_chroma_mode == UV_CFL_PRED)
                candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
            else
                candidateArray[canTotalCnt].transform_type_uv =
                av1_get_tx_type(
                    context_ptr->blk_geom->bsize,
                    0,
                    (PredictionMode)candidateArray[canTotalCnt].intra_luma_mode,
                    (UvPredictionMode)candidateArray[canTotalCnt].intra_chroma_mode,
                    PLANE_TYPE_UV,
                    0,
                    0,
                    0,
                    context_ptr->blk_geom->txsize_uv[0][0],
                    picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used);
            candidateArray[canTotalCnt].ref_frame_type = INTRA_FRAME;
            candidateArray[canTotalCnt].pred_mode = (PredictionMode)openLoopIntraCandidate;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
        }
    }

    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;

    return;
}
/***************************************
* ProductGenerateMdCandidatesCu
*   Creates list of initial modes to
*   perform fast cost search on.
***************************************/
EbErrorType ProductGenerateMdCandidatesCu(
    LargestCodingUnit                 *sb_ptr,
    ModeDecisionContext             *context_ptr,
    SsMeContext                    *ss_mecontext,
    const uint32_t                      lcuAddr,
    uint32_t                           *candidateTotalCountPtr,
    EbPtr                              interPredContextPtr,
    PictureControlSet              *picture_control_set_ptr)
{
    (void)lcuAddr;
    (void)interPredContextPtr;
    const SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    const EB_SLICE slice_type = picture_control_set_ptr->slice_type;
    uint32_t canTotalCnt = 0;
    // Reset duplicates variables
    context_ptr->injected_mv_count_l0 = 0;
    context_ptr->injected_mv_count_l1 = 0;
    context_ptr->injected_mv_count_bipred = 0;
    uint8_t sq_index = LOG2F(context_ptr->blk_geom->sq_size) - 2;
    uint8_t inject_intra_candidate = 1;
    uint8_t inject_inter_candidate = 1;

    if (slice_type != I_SLICE) {
        if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level >= NSQ_SEARCH_LEVEL1 &&
            picture_control_set_ptr->parent_pcs_ptr->nsq_search_level < NSQ_SEARCH_FULL) {
            inject_intra_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
                context_ptr->parent_sq_has_coeff[sq_index] != 0 ? inject_intra_candidate : 0;
        }
}
    //----------------------
    // Intra
    if (context_ptr->blk_geom->sq_size < 128) {
        if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode >= 5 && context_ptr->blk_geom->sq_size > 4 && context_ptr->blk_geom->shape == PART_N)
            inject_intra_candidates_ois(
                picture_control_set_ptr,
                context_ptr,
                sb_ptr,
                &canTotalCnt);
        else
            if (inject_intra_candidate)
            inject_intra_candidates(
                picture_control_set_ptr,
                context_ptr,
                sequence_control_set_ptr,
                sb_ptr,
                &canTotalCnt);
    }

    if (picture_control_set_ptr->parent_pcs_ptr->allow_intrabc)
        inject_intra_bc_candidates(
            picture_control_set_ptr,
            context_ptr,
            sequence_control_set_ptr,
            sb_ptr,
            context_ptr->cu_ptr,
            &canTotalCnt
        );

    // Track the total number of fast intra candidates
    context_ptr->fast_candidate_intra_count = canTotalCnt;

    if (slice_type != I_SLICE) {
        if (inject_inter_candidate)
            inject_inter_candidates(
                picture_control_set_ptr,
                context_ptr,
                ss_mecontext,
                sequence_control_set_ptr,
                sb_ptr,
                &canTotalCnt);
    }
    *candidateTotalCountPtr = canTotalCnt;
    return EB_ErrorNone;
}

/***************************************
* Full Mode Decision
***************************************/
uint8_t product_full_mode_decision(
    struct ModeDecisionContext   *context_ptr,
    CodingUnit                   *cu_ptr,
    uint8_t                           bwidth,
    uint8_t                           bheight,
    ModeDecisionCandidateBuffer **buffer_ptr_array,
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
    PredictionUnit       *pu_ptr;
    uint32_t                   i;
    ModeDecisionCandidate       *candidate_ptr;

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

    if (candidate_ptr->type == INTER_MODE && candidate_ptr->merge_flag == EB_TRUE)
        context_ptr->md_ep_pipe_sb[cu_ptr->mds_idx].chroma_distortion = buffer_ptr_array[lowestCostIndex]->candidate_ptr->chroma_distortion;
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].full_distortion = buffer_ptr_array[lowestCostIndex]->candidate_ptr->full_distortion;
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].chroma_distortion = (uint32_t)buffer_ptr_array[lowestCostIndex]->candidate_ptr->chroma_distortion;
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].chroma_distortion_inter_depth = (uint32_t)buffer_ptr_array[lowestCostIndex]->candidate_ptr->chroma_distortion_inter_depth;

    cu_ptr->prediction_mode_flag = candidate_ptr->type;
    cu_ptr->tx_depth = candidate_ptr->tx_depth;
    cu_ptr->skip_flag = candidate_ptr->skip_flag; // note, the skip flag is re-checked in the ENCDEC process
    cu_ptr->block_has_coeff = ((candidate_ptr->block_has_coeff) > 0) ? EB_TRUE : EB_FALSE;
    cu_ptr->quantized_dc[1][0] = buffer_ptr_array[lowestCostIndex]->candidate_ptr->quantized_dc[1][0];
    cu_ptr->quantized_dc[2][0] = buffer_ptr_array[lowestCostIndex]->candidate_ptr->quantized_dc[2][0];
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].count_non_zero_coeffs = candidate_ptr->count_non_zero_coeffs;

    cu_ptr->av1xd->use_intrabc = candidate_ptr->use_intrabc;

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
        if (cu_ptr->prediction_mode_flag != INTER_MODE && cu_ptr->av1xd->use_intrabc == 0)
        {
            pu_ptr->inter_pred_direction_index = 0x03;
            pu_ptr->merge_flag = EB_FALSE;
        }
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
        pu_ptr->ref_frame_index_l0 = candidate_ptr->ref_frame_index_l0;
        pu_ptr->ref_frame_index_l1 = candidate_ptr->ref_frame_index_l1;
        pu_ptr->ref_mv_index = candidate_ptr->ref_mv_index;
        pu_ptr->is_new_mv = candidate_ptr->is_new_mv;
        pu_ptr->is_zero_mv = candidate_ptr->is_zero_mv;

        if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_0)
        {
            //EB_MEMCPY(&pu_ptr->mv[REF_LIST_0].x,&candidate_ptr->mvs_l0,4);
            pu_ptr->mv[REF_LIST_0].x = candidate_ptr->motion_vector_xl0;
            pu_ptr->mv[REF_LIST_0].y = candidate_ptr->motion_vector_yl0;
        }

        if (pu_ptr->inter_pred_direction_index == UNI_PRED_LIST_1)
        {
            //EB_MEMCPY(&pu_ptr->mv[REF_LIST_1].x,&candidate_ptr->mvs_l1,4);
            pu_ptr->mv[REF_LIST_1].x = candidate_ptr->motion_vector_xl1;
            pu_ptr->mv[REF_LIST_1].y = candidate_ptr->motion_vector_yl1;
        }

        if (pu_ptr->inter_pred_direction_index == BI_PRED)
        {
            //EB_MEMCPY(&pu_ptr->mv[REF_LIST_0].x,&candidate_ptr->mvs,8);
            pu_ptr->mv[REF_LIST_0].x = candidate_ptr->motion_vector_xl0;
            pu_ptr->mv[REF_LIST_0].y = candidate_ptr->motion_vector_yl0;
            pu_ptr->mv[REF_LIST_1].x = candidate_ptr->motion_vector_xl1;
            pu_ptr->mv[REF_LIST_1].y = candidate_ptr->motion_vector_yl1;
        }

        // The MV prediction indicies are recalcated by the EncDec.
        pu_ptr->mvd[REF_LIST_0].pred_idx = 0;
        pu_ptr->mvd[REF_LIST_1].pred_idx = 0;

        pu_ptr->overlappable_neighbors[0] = context_ptr->cu_ptr->prediction_unit_array[0].overlappable_neighbors[0];
        pu_ptr->overlappable_neighbors[1] = context_ptr->cu_ptr->prediction_unit_array[0].overlappable_neighbors[1];
        pu_ptr->motion_mode = candidate_ptr->motion_mode;
        pu_ptr->num_proj_ref = candidate_ptr->num_proj_ref;
        if (pu_ptr->motion_mode == WARPED_CAUSAL)
            EB_MEMCPY(&pu_ptr->wm_params, &candidate_ptr->wm_params, sizeof(EbWarpedMotionParams));
    }

    TransformUnit *txb_ptr;
    uint32_t txb_itr;
    uint32_t tu_index;
    uint32_t tuTotalCount;
    uint32_t cu_size_log2 = context_ptr->cu_size_log2;
    tuTotalCount = context_ptr->blk_geom->txb_count[cu_ptr->tx_depth];
    tu_index = 0;
    txb_itr = 0;
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
        txb_ptr->transform_type[PLANE_TYPE_Y] = candidate_ptr->transform_type[tu_index];
        txb_ptr->transform_type[PLANE_TYPE_UV] = candidate_ptr->transform_type_uv;

        cu_ptr->quantized_dc[0][tu_index] = candidate_ptr->quantized_dc[0][tu_index];

#if NO_ENCDEC

        if (context_ptr->blk_geom->has_uv) {
            cu_ptr->block_has_coeff |= txb_ptr->y_has_coeff;
            cu_ptr->block_has_coeff |= txb_ptr->u_has_coeff;
            cu_ptr->block_has_coeff |= txb_ptr->v_has_coeff;
        }
        else
            cu_ptr->block_has_coeff |= txb_ptr->y_has_coeff;
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

            int32_t* src_ptr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residual_quant_coeff_ptr->buffer_y)[txb_1d_offset]);
            int32_t* dst_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->buffer_y)[txb_1d_offset]);

            uint32_t j;

            for (j = 0; j < bheight; j++)
                memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
            if (context_ptr->blk_geom->has_uv)
            {
                // Cb
                bwidth = context_ptr->blk_geom->tx_width_uv[txb_itr];
                bheight = context_ptr->blk_geom->tx_height_uv[txb_itr];

                src_ptr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residual_quant_coeff_ptr->buffer_cb)[txb_1d_offset_uv]);
                dst_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->buffer_cb)[txb_1d_offset_uv]);

                for (j = 0; j < bheight; j++)
                    memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
                src_ptr = &(((int32_t*)buffer_ptr_array[lowestCostIndex]->residual_quant_coeff_ptr->buffer_cr)[txb_1d_offset_uv]);
                dst_ptr = &(((int32_t*)context_ptr->cu_ptr->coeff_tmp->buffer_cr)[txb_1d_offset_uv]);

                for (j = 0; j < bheight; j++)
                    memcpy(dst_ptr + j * bwidth, src_ptr + j * bwidth, bwidth * sizeof(int32_t));
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
