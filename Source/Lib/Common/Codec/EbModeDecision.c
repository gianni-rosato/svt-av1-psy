// clang-format off
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

#include "EbCommonUtils.h"
#include "EbSequenceControlSet.h"
#include "EbModeDecision.h"
#include "EbTransformUnit.h"
#include "EbModeDecisionProcess.h"
#include "EbMotionEstimation.h"


#include "av1me.h"
#include "hash.h"
#include "EbInterPrediction.h"
#if II_COMP_FLAG
#include "EbRateDistortionCost.h"
#endif
#include "aom_dsp_rtcd.h"
#define  INCRMENT_CAND_TOTAL_COUNT(cnt) MULTI_LINE_MACRO_BEGIN cnt++; if(cnt>=MODE_DECISION_CANDIDATE_MAX_COUNT) printf(" ERROR: reaching limit for MODE_DECISION_CANDIDATE_MAX_COUNT %i\n",cnt); MULTI_LINE_MACRO_END

int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
int av1_filter_intra_allowed_bsize(uint8_t enable_filter_intra, BlockSize bs);
#if OBMC_FLAG
#define INT_MAX       2147483647    // maximum (signed) int value
#endif

static EB_AV1_INTER_PREDICTION_FUNC_PTR   av1_inter_prediction_function_table[2] =
{
    av1_inter_prediction,
    av1_inter_prediction_hbd
};

#if II_COMP_FLAG
void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);

static INLINE int is_interintra_allowed_bsize(const BlockSize bsize) {
    return (bsize >= BLOCK_8X8) && (bsize <= BLOCK_32X32);
}

static INLINE int is_interintra_allowed_mode(const PredictionMode mode) {
    return (mode >= SINGLE_INTER_MODE_START) && (mode < SINGLE_INTER_MODE_END);
}

static INLINE int is_interintra_allowed_ref(const MvReferenceFrame rf[2]) {
    return (rf[0] > INTRA_FRAME) && (rf[1] <= INTRA_FRAME);
}
int svt_is_interintra_allowed(
    uint8_t enable_inter_intra,
    BlockSize sb_type,
    PredictionMode mode,
    MvReferenceFrame ref_frame[2])
{
    return
        enable_inter_intra &&
        is_interintra_allowed_bsize((const BlockSize)sb_type) &&
        is_interintra_allowed_mode(mode)  &&
        is_interintra_allowed_ref(ref_frame);

}
#endif
/********************************************
* Constants
********************************************/
#if II_COMP_FLAG
// 1 - Regular uni-pred ,
// 2 - Regular uni-pred + Wedge compound Inter Intra
// 3 - Regular uni-pred + Wedge compound Inter Intra + Smooth compound Inter Intra

#define II_COUNT                3
#endif
#if OBMC_FLAG
#if MULTI_PASS_PD
static INLINE int is_inter_mode(PredictionMode mode)
{
    return mode >= SINGLE_INTER_MODE_START && mode < SINGLE_INTER_MODE_END;
}
MotionMode obmc_motion_mode_allowed(
    const PictureControlSet    *picture_control_set_ptr,
    struct ModeDecisionContext *context_ptr,
    const BlockSize             bsize,
    MvReferenceFrame            rf0,
    MvReferenceFrame            rf1,
    PredictionMode              mode)
{

    if (!context_ptr->md_pic_obmc_mode)
        return SIMPLE_TRANSLATION;

    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

    if (!frm_hdr->is_motion_mode_switchable)
        return SIMPLE_TRANSLATION;

    if (frm_hdr->force_integer_mv == 0) {
        const TransformationType gm_type =
            picture_control_set_ptr->parent_pcs_ptr->global_motion[rf0].wmtype;
        if (is_global_mv_block(mode, bsize, gm_type))
            return SIMPLE_TRANSLATION;
    }

    if (is_motion_variation_allowed_bsize(bsize) &&
        is_inter_mode(mode) &&
        rf1 != INTRA_FRAME &&
        !(rf1 > INTRA_FRAME)) // is_motion_variation_allowed_compound
    {
        if (!has_overlappable_candidates(context_ptr->cu_ptr)) // check_num_overlappable_neighbors
            return SIMPLE_TRANSLATION;

        return OBMC_CAUSAL;
    }
    else
        return SIMPLE_TRANSLATION;
}
#else
MotionMode obmc_motion_mode_allowed(
    const PictureControlSet       *picture_control_set_ptr,
    const CodingUnit              *cu_ptr,
    const BlockSize                 bsize,
    MvReferenceFrame                rf0,
    MvReferenceFrame                rf1,
    PredictionMode                  mode);
#endif
void precompute_obmc_data(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr);
#endif
//static uint32_t  AntiContouringIntraMode[11] = { EB_INTRA_PLANAR, EB_INTRA_DC, EB_INTRA_HORIZONTAL, EB_INTRA_VERTICAL,
//EB_INTRA_MODE_2, EB_INTRA_MODE_6, EB_INTRA_MODE_14, EB_INTRA_MODE_18, EB_INTRA_MODE_22, EB_INTRA_MODE_30, EB_INTRA_MODE_34 };
int32_t have_newmv_in_inter_mode(PredictionMode mode) {

    return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAREST_NEWMV ||
        mode == NEW_NEARESTMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}
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
int32_t eb_av1_mv_bit_cost(const MV *mv, const MV *ref, const int32_t *mvjcost,
    int32_t *mvcost[2], int32_t weight);
#define MV_COST_WEIGHT 108

#if II_COMP_FLAG
#define MAX_INTERINTRA_SB_SQUARE 32 * 32

EbErrorType  intra_luma_prediction_for_interintra(
    ModeDecisionContext         *md_context_ptr,
    PictureControlSet           *picture_control_set_ptr,
    INTERINTRA_MODE              interintra_mode,
    EbPictureBufferDesc         *prediction_ptr);
int64_t pick_wedge_fixed_sign(
    ModeDecisionCandidate        *candidate_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    ModeDecisionContext                  *context_ptr,
    const BlockSize bsize,
    const int16_t *const residual1,
    const int16_t *const diff10,
    const int8_t wedge_sign,
    int8_t *const best_wedge_index);
void model_rd_for_sb_with_curvfit(
    PictureControlSet      *picture_control_set_ptr,
    ModeDecisionContext                  *context_ptr,
    BlockSize bsize, int bw, int bh,
    uint8_t* src_buf, uint32_t src_stride, uint8_t* pred_buf, uint32_t pred_stride,
    int plane_from, int plane_to, int mi_row, int mi_col, int *out_rate_sum,
    int64_t *out_dist_sum, int *skip_txfm_sb, int64_t *skip_sse_sb,
    int *plane_rate, int64_t *plane_sse, int64_t *plane_dist);

static int64_t pick_interintra_wedge(
    ModeDecisionCandidate        *candidate_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    ModeDecisionContext                  *context_ptr,
    const BlockSize bsize,
    const uint8_t *const p0,
    const uint8_t *const p1,
    uint8_t * src_buf,
    uint32_t  src_stride,
    int32_t *wedge_index_out
    )
{

    assert(is_interintra_wedge_used(bsize));
   // assert(cpi->common.seq_params.enable_interintra_compound);

    const int bw = block_size_wide[bsize];
    const int bh = block_size_high[bsize];
    DECLARE_ALIGNED(32, int16_t, residual1[MAX_SB_SQUARE]);  // src - pred1
    DECLARE_ALIGNED(32, int16_t, diff10[MAX_SB_SQUARE]);     // pred1 - pred0
    if (context_ptr->hbd_mode_decision)
    {
        aom_highbd_subtract_block(bh, bw, residual1, bw,  src_buf, src_stride, p1, bw, EB_10BIT);
        aom_highbd_subtract_block(bh, bw, diff10, bw, p1, bw,  p0, bw, EB_10BIT);

    }else
    {
        aom_subtract_block(bh, bw, residual1, bw, src_buf, src_stride, p1, bw);
        aom_subtract_block(bh, bw, diff10, bw, p1, bw, p0, bw);
    }

    int8_t wedge_index = -1;
    int64_t rd =
        pick_wedge_fixed_sign(candidate_ptr,picture_control_set_ptr, context_ptr, bsize, residual1, diff10, 0, &wedge_index);

    *wedge_index_out = wedge_index;

    return rd;
}
#if II_COMP_FLAG
//for every CU, perform DC/V/H/S intra prediction to be used later in inter-intra search
void precompute_intra_pred_for_inter_intra(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr)
{
    uint32_t j;
    EbPictureBufferDesc  pred_desc;
    pred_desc.origin_x = pred_desc.origin_y = 0;
    pred_desc.stride_y = context_ptr->blk_geom->bwidth;

    for (j = 0; j < INTERINTRA_MODES; ++j)
    {

        INTERINTRA_MODE interintra_mode = (INTERINTRA_MODE)j;
        pred_desc.buffer_y = context_ptr->intrapred_buf[j];
        intra_luma_prediction_for_interintra(
            context_ptr,
            picture_control_set_ptr,
            interintra_mode,
            &pred_desc);
    }
}
#endif
 void combine_interintra(INTERINTRA_MODE mode,
    int8_t use_wedge_interintra, int wedge_index,
    int wedge_sign, BlockSize bsize,
    BlockSize plane_bsize, uint8_t *comppred,
    int compstride, const uint8_t *interpred,
    int interstride, const uint8_t *intrapred,
    int intrastride);
void inter_intra_search(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    ModeDecisionCandidate        *candidate_ptr)
{
    DECLARE_ALIGNED(16, uint8_t, tmp_buf[ 2* MAX_INTERINTRA_SB_SQUARE]);
#if !II_COMP_FLAG
    DECLARE_ALIGNED(16, uint8_t, intrapred[ MAX_INTERINTRA_SB_SQUARE]);
#endif
    DECLARE_ALIGNED(16, uint8_t, ii_pred_buf[2*MAX_INTERINTRA_SB_SQUARE]);
    //get inter pred for ref0
    EbPictureBufferDesc   *src_pic = context_ptr->hbd_mode_decision ? picture_control_set_ptr->input_frame16bit :  picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    uint16_t *src_buf_hbd = (uint16_t*)src_pic->buffer_y + (context_ptr->cu_origin_x + src_pic->origin_x) + (context_ptr->cu_origin_y + src_pic->origin_y) * src_pic->stride_y;
    uint8_t *src_buf = src_pic->buffer_y + (context_ptr->cu_origin_x + src_pic->origin_x) + (context_ptr->cu_origin_y + src_pic->origin_y) * src_pic->stride_y;

    uint8_t bit_depth = context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT;

    uint32_t  bwidth = context_ptr->blk_geom->bwidth;
    uint32_t  bheight = context_ptr->blk_geom->bheight;
    EbPictureBufferDesc  pred_desc;
    pred_desc.origin_x = pred_desc.origin_y = 0;
    pred_desc.stride_y = bwidth;

    EbPictureBufferDesc  *ref_pic_list0;
    EbPictureBufferDesc  *ref_pic_list1 = NULL;
    Mv mv_0;
    Mv mv_1;
    mv_0.x = candidate_ptr->motion_vector_xl0;
    mv_0.y = candidate_ptr->motion_vector_yl0;
    mv_1.x = candidate_ptr->motion_vector_xl1;
    mv_1.y = candidate_ptr->motion_vector_yl1;
    MvUnit mv_unit;
    mv_unit.mv[0] = mv_0;
    mv_unit.mv[1] = mv_1;
    int8_t ref_idx_l0 = candidate_ptr->ref_frame_index_l0;
    int8_t ref_idx_l1 = candidate_ptr->ref_frame_index_l1;
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, candidate_ptr->ref_frame_type);
    uint8_t list_idx0, list_idx1;
    list_idx0 = get_list_idx(rf[0]);
    if (rf[1] == NONE_FRAME)
        list_idx1 = get_list_idx(rf[0]);
    else
        list_idx1 = get_list_idx(rf[1]);
    assert(list_idx0 < MAX_NUM_OF_REF_PIC_LIST);
    assert(list_idx1 < MAX_NUM_OF_REF_PIC_LIST);
//
    if (ref_idx_l0 >= 0)
        ref_pic_list0 = context_ptr->hbd_mode_decision ?
                        ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[list_idx0][ref_idx_l0]->object_ptr)->reference_picture16bit :
                        ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[list_idx0][ref_idx_l0]->object_ptr)->reference_picture;
    else
        ref_pic_list0 = (EbPictureBufferDesc*)EB_NULL;

    if (ref_idx_l1 >= 0)
        ref_pic_list1 = context_ptr->hbd_mode_decision ?
                        ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[list_idx1][ref_idx_l1]->object_ptr)->reference_picture16bit :
                        ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[list_idx1][ref_idx_l1]->object_ptr)->reference_picture;
    else
        ref_pic_list1 = (EbPictureBufferDesc*)EB_NULL;

    mv_unit.pred_direction = candidate_ptr->prediction_direction[0];

    pred_desc.buffer_y = tmp_buf;

    //we call the regular inter prediction path here(no compound)
    av1_inter_prediction_function_table[context_ptr->hbd_mode_decision > EB_8_BIT_MD](
        picture_control_set_ptr,
        0,//ASSUMPTION: fixed interpolation filter.
        context_ptr->cu_ptr,
        candidate_ptr->ref_frame_type,
        &mv_unit,
        0,//use_intrabc,
#if OBMC_FLAG
        SIMPLE_TRANSLATION,
        0,
        0,
#endif
        1,//compound_idx not used
        NULL,// interinter_comp not used
#if II_COMP_FLAG
        NULL,
        NULL,
        NULL,
        NULL,
        0,
        0,
        0,
        0,
#endif
        context_ptr->cu_origin_x,
        context_ptr->cu_origin_y,
        bwidth,
        bheight,
        ref_pic_list0,
        ref_pic_list1,
        &pred_desc, //output
        0,          //output origin_x,
        0,          //output origin_y,
        0, //do chroma
        context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT);

    assert(is_interintra_wedge_used(context_ptr->blk_geom->bsize));//if not I need to add nowedge path!!

    int64_t rd = INT64_MAX;
    int64_t best_interintra_rd = INT64_MAX;
    int rmode = 0, rate_sum;
    int64_t dist_sum;
    int tmp_rate_mv = 0;


    INTERINTRA_MODE best_interintra_mode = INTERINTRA_MODES;
#if !II_COMP_FLAG
    EbPictureBufferDesc  inra_pred_desc;
    inra_pred_desc.origin_x = inra_pred_desc.origin_y = 0;
    inra_pred_desc.stride_y = bwidth;
    inra_pred_desc.buffer_y = intrapred;
#endif
    int8_t enable_smooth_interintra =1;
      //if (cpi->oxcf.enable_smooth_interintra &&
      //!cpi->sf.disable_smooth_interintra) {
    if (enable_smooth_interintra) {
        int j = 0;
        if (/*cpi->sf.reuse_inter_intra_mode == 0 ||*/
            best_interintra_mode == INTERINTRA_MODES) {
            for (j = 0; j < INTERINTRA_MODES; ++j) {
                //if ((!cpi->oxcf.enable_smooth_intra || cpi->sf.disable_smooth_intra) &&
                //    (INTERINTRA_MODE)j == II_SMOOTH_PRED)
                //  continue;
                INTERINTRA_MODE interintra_mode = (INTERINTRA_MODE)j;
                //rmode = interintra_mode_cost[mbmi->interintra_mode];
                const int bsize_group = size_group_lookup[context_ptr->blk_geom->bsize];
                rmode  = candidate_ptr->md_rate_estimation_ptr->inter_intra_mode_fac_bits[bsize_group][interintra_mode];
#if !II_COMP_FLAG
                inra_pred_desc.buffer_y = context_ptr->intrapred_buf[j];

                //av1_build_intra_predictors_for_interintra(cm, xd, bsize, 0, orig_dst,
                //                                          intrapred, bw);
                intra_luma_prediction_for_interintra(
                    context_ptr,
                    picture_control_set_ptr,
                    interintra_mode,
                    &inra_pred_desc);
#endif

                //av1_combine_interintra(xd, bsize, 0, tmp_buf, bw, intrapred, bw);
                if (context_ptr->hbd_mode_decision)
                    combine_interintra_highbd(
                        interintra_mode,//mode,
                        0,//use_wedge_interintra,
                        0,//candidate_ptr->interintra_wedge_index,
                        0,//int wedge_sign,
                        context_ptr->blk_geom->bsize,
                        context_ptr->blk_geom->bsize,// plane_bsize,
                         ii_pred_buf, bwidth, /*uint8_t *comppred, int compstride,*/
                        tmp_buf, bwidth,  /*const uint8_t *interpred, int interstride,*/
                        context_ptr->intrapred_buf[j], bwidth /*const uint8_t *intrapred,   int intrastride*/,bit_depth);
                else

                combine_interintra(
                    interintra_mode,//mode,
                    0,//use_wedge_interintra,
                    0,//candidate_ptr->interintra_wedge_index,
                    0,//int wedge_sign,
                    context_ptr->blk_geom->bsize,
                    context_ptr->blk_geom->bsize,// plane_bsize,
                    ii_pred_buf, bwidth, /*uint8_t *comppred, int compstride,*/
                    tmp_buf, bwidth,  /*const uint8_t *interpred, int interstride,*/
                    context_ptr->intrapred_buf[j], bwidth /*const uint8_t *intrapred,   int intrastride*/);

                //model_rd_sb_fn[MODELRD_TYPE_INTERINTRA](
                //    cpi, bsize, x, xd, 0, 0, mi_row, mi_col, &rate_sum, &dist_sum,
                //    &tmp_skip_txfm_sb, &tmp_skip_sse_sb, NULL, NULL, NULL);
                if (context_ptr->hbd_mode_decision) {

                    model_rd_for_sb_with_curvfit(picture_control_set_ptr, context_ptr, context_ptr->blk_geom->bsize, bwidth, bheight,
                        (uint8_t*)src_buf_hbd, src_pic->stride_y, ii_pred_buf, bwidth,
                        0, 0, 0, 0, &rate_sum, &dist_sum, NULL, NULL, NULL, NULL, NULL);
                }else
                {

                model_rd_for_sb_with_curvfit(picture_control_set_ptr, context_ptr, context_ptr->blk_geom->bsize, bwidth, bheight,
                    src_buf, src_pic->stride_y, ii_pred_buf, bwidth,
                    0, 0, 0, 0, &rate_sum, &dist_sum, NULL, NULL, NULL, NULL, NULL);
                }
                // rd = RDCOST(x->rdmult, tmp_rate_mv + rate_sum + rmode, dist_sum);
                rd = RDCOST(context_ptr->full_lambda, tmp_rate_mv + rate_sum + rmode, dist_sum);

                if (rd < best_interintra_rd) {
                    best_interintra_rd = rd;
                    candidate_ptr->interintra_mode = best_interintra_mode = interintra_mode;
                }
            }

            /* best_interintra_rd_wedge =
                 pick_interintra_wedge(cpi, x, bsize, intrapred_, tmp_buf_);*/

            //CHKN need to re-do intra pred using the winner, or have a separate intra serch for wedge

            pick_interintra_wedge(
                candidate_ptr,
                picture_control_set_ptr,
                context_ptr,
                context_ptr->blk_geom->bsize,
#if II_COMP_FLAG
                context_ptr->intrapred_buf[best_interintra_mode],
#else
                intrapred,
#endif
                tmp_buf,
                context_ptr->hbd_mode_decision ? (uint8_t*) src_buf_hbd : src_buf,
                src_pic->stride_y,
                &candidate_ptr->interintra_wedge_index
            );

            //if (best_interintra_rd_wedge < best_interintra_rd) {

                //candidate_ptr->use_wedge_interintra = 1;
                //candidate_ptr->ii_wedge_sign = 0;
            //}
            //args->inter_intra_mode[mbmi->ref_frame[0]] = best_interintra_mode;
        }
    }
    // Enable wedge search if source variance and edge strength are above the thresholds.



}
#endif

COMPOUND_TYPE to_av1_compound_lut[] = {
    COMPOUND_AVERAGE,
    COMPOUND_DISTWTD,
    COMPOUND_DIFFWTD,
    COMPOUND_WEDGE
};

void determine_compound_mode(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    ModeDecisionCandidate        *candidatePtr,
    MD_COMP_TYPE                 cur_type) {


    candidatePtr->interinter_comp.type = to_av1_compound_lut[cur_type];

    if (cur_type == MD_COMP_AVG) {

        candidatePtr->comp_group_idx = 0;
        candidatePtr->compound_idx = 1;
#if 0
        if (candidatePtr->merge_flag == 0)
            search_compound_avg_dist(
                picture_control_set_ptr,
                context_ptr,
                candidatePtr);
#endif
    }
    else if (cur_type == MD_COMP_DIST) {

        candidatePtr->comp_group_idx = 0;
        candidatePtr->compound_idx = 0;
    }
    else if (cur_type == MD_COMP_DIFF0) {
        candidatePtr->comp_group_idx = 1;
        candidatePtr->compound_idx = 1;
        candidatePtr->interinter_comp.mask_type = 55;
        search_compound_diff_wedge(
            picture_control_set_ptr,
            context_ptr,
            candidatePtr
        );

    }
    //else if (cur_type == MD_COMP_DIFF1) {
    //    candidatePtr->comp_group_idx = 1;
    //    candidatePtr->compound_idx = 1;
    //    candidatePtr->interinter_comp.mask_type = 1;
    //}
    else if (cur_type == MD_COMP_WEDGE) {

        candidatePtr->comp_group_idx = 1;
        candidatePtr->compound_idx = 1;
        search_compound_diff_wedge(
            picture_control_set_ptr,
            context_ptr,
            candidatePtr
        );

        candidatePtr->interinter_comp.wedge_index = candidatePtr->interinter_comp.wedge_index;
        candidatePtr->interinter_comp.wedge_sign = candidatePtr->interinter_comp.wedge_sign;
    }
    else {
        printf("ERROR: not used comp type\n");
    }
}

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

        uint32_t mvRate = (uint32_t)eb_av1_mv_bit_cost(
            &mv,
            &(ref_mv[0].as_mv),
            md_rate_estimation_ptr->nmv_vec_cost,
            md_rate_estimation_ptr->nmvcoststack,
            MV_COST_WEIGHT);

        if (is_compound) {
            mv.row = mv1y;
            mv.col = mv1x;

            mvRate += (uint32_t)eb_av1_mv_bit_cost(
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
static void mode_decision_scratch_candidate_buffer_dctor(EbPtr p)
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
    EbBitDepthEnum                  max_bitdepth,
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
    pictureBufferDescInitData.bit_depth = max_bitdepth;
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

    // Costs
    buffer_ptr->fast_cost_ptr = fast_cost_ptr;
    buffer_ptr->full_cost_ptr = full_cost_ptr;
    buffer_ptr->full_cost_skip_ptr = full_cost_skip_ptr;
    buffer_ptr->full_cost_merge_ptr = full_cost_merge_ptr;
    return EB_ErrorNone;
}
EbErrorType mode_decision_scratch_candidate_buffer_ctor(
    ModeDecisionCandidateBuffer    *buffer_ptr,
    EbBitDepthEnum                  max_bitdepth)
{

    EbPictureBufferDescInitData pictureBufferDescInitData;
    EbPictureBufferDescInitData doubleWidthPictureBufferDescInitData;
    EbPictureBufferDescInitData ThirtyTwoWidthPictureBufferDescInitData;


    buffer_ptr->dctor = mode_decision_scratch_candidate_buffer_dctor;

    // Init Picture Data
    pictureBufferDescInitData.max_width = MAX_SB_SIZE;
    pictureBufferDescInitData.max_height = MAX_SB_SIZE;
    pictureBufferDescInitData.bit_depth = max_bitdepth;
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
    return EB_ErrorNone;
}

uint8_t check_ref_beackout(
    struct ModeDecisionContext *context_ptr,
    uint8_t                     ref_frame_type,
    PART                        shape)
{
    uint8_t skip_candidate = 0;
    uint8_t ref_cnt = 0;
    uint8_t allowed_nsq_ref_th = (uint8_t)PRUNE_REC_TH;
    if (context_ptr->prune_ref_frame_for_rec_partitions) {
        if (shape != PART_N) {
            uint8_t ref_idx;
            assert(ref_frame_type < 30);
            ref_cnt = 0;
            for (ref_idx = 0; ref_idx < allowed_nsq_ref_th; ref_idx++) {
                if (ref_frame_type == context_ptr->ref_best_ref_sq_table[ref_idx]) {
                    ref_cnt++;
                }
            }
            skip_candidate = ref_cnt ? 0 : 1;
        }
    }
    return skip_candidate;
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

#define BIPRED_3x3_REFINMENT_POSITIONS 8

int8_t ALLOW_REFINEMENT_FLAG[BIPRED_3x3_REFINMENT_POSITIONS] = {  1, 0, 1, 0, 1,  0,  1, 0 };
int8_t BIPRED_3x3_X_POS[BIPRED_3x3_REFINMENT_POSITIONS] = { -1, -1, 0, 1, 1, 1, 0, -1 };
int8_t BIPRED_3x3_Y_POS[BIPRED_3x3_REFINMENT_POSITIONS] = { 0, 1, 1, 1, 0, -1, -1, -1 };

void Unipred3x3CandidatesInjection(
    const SequenceControlSet  *sequence_control_set_ptr,
    PictureControlSet         *picture_control_set_ptr,
    ModeDecisionContext       *context_ptr,
    SuperBlock                *sb_ptr,
    uint32_t                   me_sb_addr,
    uint32_t                  *candidateTotalCnt){
    UNUSED(sb_ptr);
    uint32_t                   bipredIndex;
    uint32_t                   canTotalCnt = (*candidateTotalCnt);
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[context_ptr->me_block_offset];
    const MeCandidate *me_block_results = me_results->me_candidate[context_ptr->me_block_offset];
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (frm_hdr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    IntMv  bestPredmv[2] = { {0}, {0} };
    int inside_tile = 1;
    MacroBlockD  *xd = context_ptr->cu_ptr->av1xd;
    int umv0tile = (sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0);
    uint32_t mi_row = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t mi_col = context_ptr->cu_origin_x >> MI_SIZE_LOG2;

    // (8 Best_L0 neighbors)
    //const MeLcuResults_t *meResults = pictureControlSetPtr->ParentPcsPtr->meResultsPtr[lcuAddr];
    total_me_cnt = MIN(total_me_cnt, BEST_CANDIDATE_COUNT);
    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
    {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
#if MULTI_PASS_PD
        if (list0_ref_index > context_ptr->md_max_ref_count - 1)
            continue;
#endif
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
        int16_t to_inject_mv_x;
        int16_t to_inject_mv_y;
        if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv) {
            to_inject_mv_x = (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1) + BIPRED_3x3_X_POS[bipredIndex];
            to_inject_mv_y = (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1) + BIPRED_3x3_Y_POS[bipredIndex];
        }
        else {
            to_inject_mv_x = (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            to_inject_mv_y = (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
        }
        uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
        uint8_t skip_cand = check_ref_beackout(
            context_ptr,
            to_inject_ref_type,
            context_ptr->blk_geom->shape);

        inside_tile = 1;
        if(umv0tile)
            inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x, to_inject_mv_y, mi_col, mi_row, context_ptr->blk_geom->bsize);
        skip_cand = skip_cand || (!inside_tile);
        if (!skip_cand && (context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE)) {

#if II_COMP_FLAG // 3x3  L0
             //MvReferenceFrame rf[2];
             //rf[0] = to_inject_ref_type;
             //rf[1] = -1;

             uint8_t inter_type;
             uint8_t is_ii_allowed = 0;//svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, NEWMV, rf);
             uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
             //uint8_t is_obmc_allowed =  obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
             //tot_inter_types = is_obmc_allowed ? tot_inter_types+1 : tot_inter_types;

            for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
            {
#endif
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
#if II_COMP_FLAG
            if (inter_type == 0) {
                candidateArray[canTotalCnt].is_interintra_used = 0;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            }
            else {
                if (is_ii_allowed) {
                    if (inter_type == 1) {
                        inter_intra_search(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canTotalCnt]);
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].use_wedge_interintra = 1;
                        candidateArray[canTotalCnt].ii_wedge_sign = 0;
                    }
                    else if (inter_type == 2) {
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].interintra_mode = candidateArray[canTotalCnt - 1].interintra_mode;
                        candidateArray[canTotalCnt].use_wedge_interintra = 0;
                    }
                }
                //if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                //    candidateArray[canTotalCnt].is_interintra_used = 0;
                //    candidateArray[canTotalCnt].motion_mode = OBMC_CAUSAL;
                //}
            }
#endif
            INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

#if II_COMP_FLAG
            }
#endif
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
    total_me_cnt = MIN(total_me_cnt, BEST_CANDIDATE_COUNT);
    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
    {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;
#if MULTI_PASS_PD
        if (list1_ref_index > context_ptr->md_max_ref_count - 1)
            continue;
#endif
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
            int16_t to_inject_mv_x;
            int16_t to_inject_mv_y;
            if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv) {
                to_inject_mv_x = (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].x_mv << 1) + BIPRED_3x3_X_POS[bipredIndex];
                to_inject_mv_y = (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].y_mv << 1) + BIPRED_3x3_Y_POS[bipredIndex];
            }
            else {
                to_inject_mv_x = (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].x_mv + BIPRED_3x3_X_POS[bipredIndex]) << 1;
                to_inject_mv_y = (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].y_mv + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
            }
            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
            uint8_t skip_cand = check_ref_beackout(
                context_ptr,
                to_inject_ref_type,
                context_ptr->blk_geom->shape);

            inside_tile = 1;
            if(umv0tile)
                inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x, to_inject_mv_y, mi_col, mi_row, context_ptr->blk_geom->bsize);
            skip_cand = skip_cand || (!inside_tile);
            if (!skip_cand && (context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE)) {
#if II_COMP_FLAG // 3x3  L1
             //MvReferenceFrame rf[2];
             //rf[0] = to_inject_ref_type;
             //rf[1] = -1;
             uint8_t inter_type;
             uint8_t is_ii_allowed = 0;//svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, NEWMV, rf);
             uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
             //uint8_t is_obmc_allowed =  obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
             //tot_inter_types = is_obmc_allowed ? tot_inter_types+1 : tot_inter_types;
            for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
            {
#endif
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
#if II_COMP_FLAG
            if (inter_type == 0) {
                candidateArray[canTotalCnt].is_interintra_used = 0;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            }
            else {
                if (is_ii_allowed) {
                    if (inter_type == 1) {
                        inter_intra_search(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canTotalCnt]);
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].use_wedge_interintra = 1;
                        candidateArray[canTotalCnt].ii_wedge_sign = 0;
                    }
                    else if (inter_type == 2) {
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].interintra_mode = candidateArray[canTotalCnt - 1].interintra_mode;
                        candidateArray[canTotalCnt].use_wedge_interintra = 0;
                    }
                }
                //if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                //    candidateArray[canTotalCnt].is_interintra_used = 0;
                //    candidateArray[canTotalCnt].motion_mode = OBMC_CAUSAL;
                //}
            }
#endif
                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
#if II_COMP_FLAG
            }
#endif
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
    SuperBlock               *sb_ptr,
    uint32_t                  me_sb_addr,
    uint32_t                 *candidateTotalCnt){
    UNUSED(sb_ptr);
    uint32_t                   bipredIndex;
    uint32_t                   canTotalCnt = (*candidateTotalCnt);
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[context_ptr->me_block_offset];
    const MeCandidate *me_block_results = me_results->me_candidate[context_ptr->me_block_offset];
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (frm_hdr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
    IntMv  bestPredmv[2] = { {0}, {0} };
    int inside_tile = 1;
    MacroBlockD  *xd = context_ptr->cu_ptr->av1xd;
    int umv0tile = (sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0);
    uint32_t mi_row = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t mi_col = context_ptr->cu_origin_x >> MI_SIZE_LOG2;
    MD_COMP_TYPE cur_type; //BIP 3x3
    BlockSize bsize = context_ptr->blk_geom->bsize;
    MD_COMP_TYPE tot_comp_types =
        (picture_control_set_ptr->parent_pcs_ptr->compound_mode == 1 || context_ptr->compound_types_to_try == MD_COMP_AVG) ?
            MD_COMP_AVG :
            (bsize >= BLOCK_8X8 && bsize <= BLOCK_32X32) ?
                context_ptr->compound_types_to_try :
                    context_ptr->compound_types_to_try == MD_COMP_WEDGE ? MD_COMP_DIFF0 :
                    context_ptr->compound_types_to_try;

    if (context_ptr->source_variance < context_ptr->inter_inter_wedge_variance_th)
        tot_comp_types = MIN(tot_comp_types, MD_COMP_DIFF0);

    if (isCompoundEnabled) {
        /**************
       NEW_NEWMV
       ************* */
       //const MeLcuResults_t *meResults = pictureControlSetPtr->ParentPcsPtr->meResultsPtr[lcuAddr];
        total_me_cnt = MIN(total_me_cnt, BEST_CANDIDATE_COUNT);
        for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
        {
            const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
            const uint8_t inter_direction = me_block_results_ptr->direction;
            const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
            const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;

#if MULTI_PASS_PD
            if (list0_ref_index > context_ptr->md_max_ref_count - 1 || list1_ref_index > context_ptr->md_max_ref_count - 1)
                continue;
#endif

            if (inter_direction == 2) {
       // (Best_L0, 8 Best_L1 neighbors)
        for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
        {
        if (context_ptr->bipred3x3_injection >= 2){
            if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                continue;
        }
        int16_t to_inject_mv_x_l0 =  me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1;
        int16_t to_inject_mv_y_l0 =  me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1;
        int16_t to_inject_mv_x_l1;
        int16_t to_inject_mv_y_l1;
        if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv) {
            to_inject_mv_x_l1 = (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) :
                (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv << 1) + BIPRED_3x3_X_POS[bipredIndex];
            to_inject_mv_y_l1 = (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) :
                (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv << 1) + BIPRED_3x3_Y_POS[bipredIndex];
        }
        else {
            to_inject_mv_x_l1 = (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) :
                (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv + BIPRED_3x3_X_POS[bipredIndex]) << 1;
            to_inject_mv_y_l1 = (me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) :
                (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
        }
        MvReferenceFrame rf[2];
        rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
        rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
        uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
        uint8_t skip_cand = check_ref_beackout(
            context_ptr,
            to_inject_ref_type,
            context_ptr->blk_geom->shape);

        inside_tile = 1;
        if(umv0tile) {
            inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l0, mi_col, mi_row, context_ptr->blk_geom->bsize) &&
                          is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l1, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize);
        }
        skip_cand = skip_cand || (!inside_tile);
        if (!skip_cand && (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE)) {
            context_ptr->variance_ready = 0;
            for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
            {
                if (cur_type == MD_COMP_WEDGE && wedge_params_lookup[context_ptr->blk_geom->bsize].bits == 0) continue;
                // If two predictors are very similar, skip wedge compound mode search
                if (context_ptr->variance_ready)
                    if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEW_NEWMV) && context_ptr->prediction_mse < 64))
                        continue;

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
#if II_COMP_FLAG
            candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
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
            //BIP 3x3
            determine_compound_mode(
                picture_control_set_ptr,
                context_ptr,
                &candidateArray[canTotalCnt],
                cur_type);
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

        // (8 Best_L0 neighbors, Best_L1) :
        for (bipredIndex = 0; bipredIndex < BIPRED_3x3_REFINMENT_POSITIONS; ++bipredIndex)
        {
            if (context_ptr->bipred3x3_injection >= 2){
                if (ALLOW_REFINEMENT_FLAG[bipredIndex] == 0)
                    continue;
            }
            int16_t to_inject_mv_x_l0;
            int16_t to_inject_mv_y_l0;
            if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv) {
                to_inject_mv_x_l0 = (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1) + BIPRED_3x3_X_POS[bipredIndex];
                to_inject_mv_y_l0 = (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1) + BIPRED_3x3_Y_POS[bipredIndex];
            }
            else {
                to_inject_mv_x_l0 = (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv + BIPRED_3x3_X_POS[bipredIndex]) << 1;
                to_inject_mv_y_l0 = (me_results->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv + BIPRED_3x3_Y_POS[bipredIndex]) << 1;
            }
            int16_t to_inject_mv_x_l1 =  me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv << 1;
            int16_t to_inject_mv_y_l1 =  me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv << 1;

            MvReferenceFrame rf[2];
            rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
            rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
            uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
            uint8_t skip_cand = check_ref_beackout(
                context_ptr,
                to_inject_ref_type,
                context_ptr->blk_geom->shape);

            inside_tile = 1;
            if(umv0tile) {
                inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l0, mi_col, mi_row, context_ptr->blk_geom->bsize) &&
                              is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l1, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize);
            }
            skip_cand = skip_cand || (!inside_tile);
            if (!skip_cand && (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE)) {
                context_ptr->variance_ready = 0;
                for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
                {
                    if (cur_type == MD_COMP_WEDGE && wedge_params_lookup[context_ptr->blk_geom->bsize].bits == 0) continue;
                    // If two predictors are very similar, skip wedge compound mode search
                    if (context_ptr->variance_ready)
                        if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEW_NEWMV) && context_ptr->prediction_mse < 64))
                            continue;
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
#if II_COMP_FLAG
                candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
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
                //BIP 3x3
                determine_compound_mode(
                    picture_control_set_ptr,
                    context_ptr,
                    &candidateArray[canTotalCnt],
                    cur_type);
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
     }

    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;

    return;
}

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
#if !II_COMP_FLAG
void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);
#endif
void inject_mvp_candidates_II(
    struct ModeDecisionContext     *context_ptr,
    PictureControlSet              *picture_control_set_ptr,
    CodingUnit                     *cu_ptr,
    MvReferenceFrame                 ref_pair,
    uint32_t                         *candTotCnt)
{
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    EbBool allow_compound = (frm_hdr->reference_mode == SINGLE_REFERENCE || context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : EB_TRUE;
    uint8_t inj_mv;
    uint32_t                   canIdx = *candTotCnt;
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    MacroBlockD  *xd = cu_ptr->av1xd;
    uint8_t drli, maxDrlIndex;
    IntMv    nearestmv[2], nearmv[2], ref_mv[2];

    MvReferenceFrame rf[2];
    int inside_tile = 1;
    SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet *)picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    int umv0tile = (sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0);
    uint32_t mi_row = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t mi_col = context_ptr->cu_origin_x >> MI_SIZE_LOG2;
    av1_set_ref_frame(rf, ref_pair);
    MD_COMP_TYPE cur_type; //MVP
    BlockSize bsize = context_ptr->blk_geom->bsize;                       // bloc size
    MD_COMP_TYPE tot_comp_types =
            (bsize >= BLOCK_8X8 && bsize <= BLOCK_32X32) ? context_ptr->compound_types_to_try :
                context_ptr->compound_types_to_try == MD_COMP_WEDGE ? MD_COMP_DIFF0 :
                context_ptr->compound_types_to_try;
    if (context_ptr->source_variance < context_ptr->inter_inter_wedge_variance_th)
        tot_comp_types = MIN(tot_comp_types, MD_COMP_DIFF0);
    //single ref/list
    if (rf[1] == NONE_FRAME)
    {
        MvReferenceFrame frame_type = rf[0];
        uint8_t list_idx = get_list_idx(rf[0]);
        uint8_t ref_idx = get_ref_frame_idx(rf[0]);
#if MULTI_PASS_PD
        if (ref_idx > context_ptr->md_max_ref_count - 1)
            return;
#endif
        //NEAREST
        int16_t to_inject_mv_x = context_ptr->cu_ptr->ref_mvs[frame_type][0].as_mv.col;
        int16_t to_inject_mv_y = context_ptr->cu_ptr->ref_mvs[frame_type][0].as_mv.row;

        inj_mv = list_idx == 0 ?
            context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, frame_type) == EB_FALSE :
            context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, frame_type) == EB_FALSE;

        if(umv0tile)
            inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x, to_inject_mv_y, mi_col, mi_row, context_ptr->blk_geom->bsize);
        inj_mv = inj_mv && inside_tile;
        if (inj_mv) {
#if II_COMP_FLAG // NEARESTMV
            uint8_t inter_type;
#if MULTI_PASS_PD
            uint8_t is_ii_allowed = svt_is_interintra_allowed(context_ptr->md_enable_inter_intra, bsize, NEARESTMV, rf);
#else
            uint8_t is_ii_allowed = svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, NEARESTMV, rf);
#endif
            uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
#if OBMC_FLAG
#if MULTI_PASS_PD
            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr, bsize, rf[0], rf[1], NEARESTMV) == OBMC_CAUSAL;
#else
            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEARESTMV) == OBMC_CAUSAL;
#endif
            tot_inter_types = is_obmc_allowed ? tot_inter_types + 1 : tot_inter_types;
#endif
            for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
            {
#endif
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
#if II_COMP_FLAG
            if (inter_type == 0) {
                candidateArray[canIdx].is_interintra_used = 0;
                candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
            }
            else {
                if (is_ii_allowed) {
                    if (inter_type == 1) {
                        inter_intra_search(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canIdx]);
                        candidateArray[canIdx].is_interintra_used = 1;
                        candidateArray[canIdx].use_wedge_interintra = 1;
                        candidateArray[canIdx].ii_wedge_sign = 0;
                    }
                    else if (inter_type == 2) {
                        candidateArray[canIdx].is_interintra_used = 1;
                        candidateArray[canIdx].interintra_mode = candidateArray[canIdx - 1].interintra_mode;
                        candidateArray[canIdx].use_wedge_interintra = 0;
                    }
                }
#if OBMC_FLAG
                if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                    candidateArray[canIdx].is_interintra_used = 0;
                    candidateArray[canIdx].motion_mode = OBMC_CAUSAL;
                }
#endif
            }
#endif
            INCRMENT_CAND_TOTAL_COUNT(canIdx);
#if II_COMP_FLAG
            }
#endif
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

            if(umv0tile)
                inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x, to_inject_mv_y, mi_col, mi_row, context_ptr->blk_geom->bsize);
            inj_mv = inj_mv && inside_tile;
            if (inj_mv) {
#if II_COMP_FLAG // NEARMV
            uint8_t inter_type;
#if MULTI_PASS_PD
            uint8_t is_ii_allowed = svt_is_interintra_allowed(context_ptr->md_enable_inter_intra, bsize, NEARMV, rf);
#else
            uint8_t is_ii_allowed = svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, NEARMV, rf);
#endif
            uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
#if OBMC_FLAG
#if MULTI_PASS_PD
            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr, bsize, rf[0], rf[1], NEARMV) == OBMC_CAUSAL;
#else
            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEARMV) == OBMC_CAUSAL;
#endif
            tot_inter_types = is_obmc_allowed ? tot_inter_types+1 : tot_inter_types;
#endif

            for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
            {
#endif
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
#if II_COMP_FLAG
            if (inter_type == 0) {
                candidateArray[canIdx].is_interintra_used = 0;
                candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
            }
            else {
                if (is_ii_allowed) {
                    if (inter_type == 1) {
                        inter_intra_search(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canIdx]);
                        candidateArray[canIdx].is_interintra_used = 1;
                        candidateArray[canIdx].use_wedge_interintra = 1;
                        candidateArray[canIdx].ii_wedge_sign = 0;
                    }
                    else if (inter_type == 2) {
                        candidateArray[canIdx].is_interintra_used = 1;
                        candidateArray[canIdx].interintra_mode = candidateArray[canIdx - 1].interintra_mode;
                        candidateArray[canIdx].use_wedge_interintra = 0;
                    }
                }
#if OBMC_FLAG
                if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                    candidateArray[canIdx].is_interintra_used = 0;
                    candidateArray[canIdx].motion_mode = OBMC_CAUSAL;
                }
#endif
            }
#endif
                INCRMENT_CAND_TOTAL_COUNT(canIdx);
#if II_COMP_FLAG
            }
#endif
            }
        }
    }
    else if (allow_compound)
    {
        uint8_t ref_idx_0 = get_ref_frame_idx(rf[0]);
        uint8_t ref_idx_1 = get_ref_frame_idx(rf[1]);
#if MULTI_PASS_PD
        if (ref_idx_0 > context_ptr->md_max_ref_count - 1 || ref_idx_1 > context_ptr->md_max_ref_count - 1)
            return;
#endif
        {
            //NEAREST_NEAREST
            int16_t to_inject_mv_x_l0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].this_mv.as_mv.col;
            int16_t to_inject_mv_y_l0 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].this_mv.as_mv.row;
            int16_t to_inject_mv_x_l1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].comp_mv.as_mv.col;
            int16_t to_inject_mv_y_l1 = context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].comp_mv.as_mv.row;

            inj_mv = context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, ref_pair) == EB_FALSE;

            if(umv0tile) {
                inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l0, mi_col, mi_row, context_ptr->blk_geom->bsize) &&
                              is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l1, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize);
            }
            inj_mv = inj_mv && inside_tile;
            if (inj_mv) {

                context_ptr->variance_ready = 0;
                for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
                {
                    if (cur_type == MD_COMP_WEDGE && wedge_params_lookup[context_ptr->blk_geom->bsize].bits == 0) continue;
                    // If two predictors are very similar, skip wedge compound mode search
                    if (context_ptr->variance_ready)
                        if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEAREST_NEARESTMV) && context_ptr->prediction_mse < 64))
                            continue;
                candidateArray[canIdx].type = INTER_MODE;
                candidateArray[canIdx].inter_mode = NEAREST_NEARESTMV;
                candidateArray[canIdx].pred_mode = NEAREST_NEARESTMV;
                candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canIdx].is_compound = 1;
#if II_COMP_FLAG
                candidateArray[canIdx].is_interintra_used = 0;
#endif
                candidateArray[canIdx].distortion_ready = 0;
                candidateArray[canIdx].use_intrabc = 0;

                candidateArray[canIdx].merge_flag =
                    cur_type == MD_COMP_AVG &&
                    picture_control_set_ptr->parent_pcs_ptr->is_skip_mode_allowed &&
                    (rf[0] == frm_hdr->skip_mode_params.ref_frame_idx_0 + 1) &&
                    (rf[1] == frm_hdr->skip_mode_params.ref_frame_idx_1 + 1) ? EB_TRUE : EB_FALSE;

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
                //NRST-NRST
                determine_compound_mode(
                    picture_control_set_ptr,
                    context_ptr,
                    &candidateArray[canIdx],
                    cur_type);
                INCRMENT_CAND_TOTAL_COUNT(canIdx);
                }
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

                if(umv0tile) {
                    inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l0, mi_col, mi_row, context_ptr->blk_geom->bsize) &&
                                  is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l1, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize);
                }
                inj_mv = inj_mv && inside_tile;
                if (inj_mv) {

                    context_ptr->variance_ready = 0;
                    for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
                    {
                        if (cur_type == MD_COMP_WEDGE && wedge_params_lookup[context_ptr->blk_geom->bsize].bits == 0) continue;
                        // If two predictors are very similar, skip wedge compound mode search
                        if (context_ptr->variance_ready)
                            if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEAR_NEARMV) && context_ptr->prediction_mse < 64))
                                continue;
                    candidateArray[canIdx].type = INTER_MODE;
                    candidateArray[canIdx].inter_mode = NEAR_NEARMV;
                    candidateArray[canIdx].pred_mode = NEAR_NEARMV;
                    candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                    candidateArray[canIdx].is_compound = 1;
#if II_COMP_FLAG
                    candidateArray[canIdx].is_interintra_used = 0;
#endif
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
                    //NR-NR
                    determine_compound_mode(
                        picture_control_set_ptr,
                        context_ptr,
                        &candidateArray[canIdx],
                        cur_type);
                    INCRMENT_CAND_TOTAL_COUNT(canIdx);
                    }
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
    MacroBlockD  *xd = context_ptr->cu_ptr->av1xd;
    IntMv    nearestmv[2], nearmv[2], ref_mv[2];
    uint8_t drli, maxDrlIndex;
    int inside_tile = 1;
    int umv0tile = (sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0);
    uint32_t mi_row = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t mi_col = context_ptr->cu_origin_x >> MI_SIZE_LOG2;

    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_pair);
    MD_COMP_TYPE cur_type; //N_NR N_NRST
    BlockSize bsize = context_ptr->blk_geom->bsize;                       // bloc size
    MD_COMP_TYPE tot_comp_types =
        (bsize >= BLOCK_8X8 && bsize <= BLOCK_32X32) ?
            context_ptr->compound_types_to_try :
            context_ptr->compound_types_to_try == MD_COMP_WEDGE ?
                MD_COMP_DIFF0 :
                context_ptr->compound_types_to_try;

    if (context_ptr->source_variance < context_ptr->inter_inter_wedge_variance_th)
        tot_comp_types = MIN(tot_comp_types, MD_COMP_DIFF0);
    {
        uint8_t ref_idx_0 = get_ref_frame_idx(rf[0]);
        uint8_t ref_idx_1 = get_ref_frame_idx(rf[1]);


#if MULTI_PASS_PD
        if (ref_idx_0 > context_ptr->md_max_ref_count - 1 || ref_idx_1 > context_ptr->md_max_ref_count - 1)
            return;
#endif

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

                if(umv0tile) {
                    inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l0, mi_col, mi_row, context_ptr->blk_geom->bsize) &&
                                  is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l1, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize);
                }
                inj_mv = inj_mv && inside_tile;
                if (inj_mv) {

                    context_ptr->variance_ready = 0;
                    for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
                    {
                        if (cur_type == MD_COMP_WEDGE && wedge_params_lookup[context_ptr->blk_geom->bsize].bits == 0) continue;
                        // If two predictors are very similar, skip wedge compound mode search
                        if (context_ptr->variance_ready)
                            if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEAREST_NEWMV) && context_ptr->prediction_mse < 64))
                                continue;
                    candidateArray[canIdx].type = INTER_MODE;
                    candidateArray[canIdx].inter_mode = NEAREST_NEWMV;
                    candidateArray[canIdx].pred_mode = NEAREST_NEWMV;
                    candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                    candidateArray[canIdx].is_compound = 1;
#if II_COMP_FLAG
                    candidateArray[canIdx].is_interintra_used = 0;
#endif
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
                    get_av1_mv_pred_drl(
                        context_ptr,
                        context_ptr->cu_ptr,
                        candidateArray[canIdx].ref_frame_type,
                        candidateArray[canIdx].is_compound,
                        NEAREST_NEWMV,
                        0,//not needed drli,
                        nearestmv,
                        nearmv,
                        ref_mv);
                    candidateArray[canIdx].motion_vector_pred_x[REF_LIST_1] = ref_mv[1].as_mv.col;
                    candidateArray[canIdx].motion_vector_pred_y[REF_LIST_1] = ref_mv[1].as_mv.row;
                    context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                    context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                    context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                    context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                    context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = ref_pair;
                    ++context_ptr->injected_mv_count_bipred;
                    //NRST_N
                    determine_compound_mode(
                        picture_control_set_ptr,
                        context_ptr,
                        &candidateArray[canIdx],
                        cur_type);
                    INCRMENT_CAND_TOTAL_COUNT(canIdx);
                    }
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

                if(umv0tile) {
                    inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l0, mi_col, mi_row, context_ptr->blk_geom->bsize) &&
                                  is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l1, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize);
                }
                inj_mv = inj_mv && inside_tile;
                if (inj_mv)
                {

                    context_ptr->variance_ready = 0;
                    for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
                    {
                        if (cur_type == MD_COMP_WEDGE && wedge_params_lookup[context_ptr->blk_geom->bsize].bits == 0) continue;
                        // If two predictors are very similar, skip wedge compound mode search
                        if (context_ptr->variance_ready)
                            if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEW_NEARESTMV) && context_ptr->prediction_mse < 64))
                                continue;
                    candidateArray[canIdx].type = INTER_MODE;
                    candidateArray[canIdx].inter_mode = NEW_NEARESTMV;
                    candidateArray[canIdx].pred_mode = NEW_NEARESTMV;
                    candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                    candidateArray[canIdx].is_compound = 1;
#if II_COMP_FLAG
                    candidateArray[canIdx].is_interintra_used = 0;
#endif
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
                    get_av1_mv_pred_drl(
                        context_ptr,
                        context_ptr->cu_ptr,
                        candidateArray[canIdx].ref_frame_type,
                        candidateArray[canIdx].is_compound,
                        NEW_NEARESTMV,
                        0,//not needed drli,
                        nearestmv,
                        nearmv,
                        ref_mv);
                    candidateArray[canIdx].motion_vector_pred_x[REF_LIST_0] = ref_mv[0].as_mv.col;
                    candidateArray[canIdx].motion_vector_pred_y[REF_LIST_0] = ref_mv[0].as_mv.row;
                    context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                    context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                    context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                    context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                    context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = ref_pair;
                    ++context_ptr->injected_mv_count_bipred;
                    //N_NRST
                    determine_compound_mode(
                        picture_control_set_ptr,
                        context_ptr,
                        &candidateArray[canIdx],
                        cur_type);
                    INCRMENT_CAND_TOTAL_COUNT(canIdx);
                    }
                }
            }
             //NEW_NEARMV
            {
                maxDrlIndex = GetMaxDrlIndex(xd->ref_mv_count[ref_pair], NEW_NEARMV);

                for (drli = 0; drli < maxDrlIndex; drli++) {
                    get_av1_mv_pred_drl(
                        context_ptr,
                        context_ptr->cu_ptr,
                        ref_pair,
                        1,
                        NEW_NEARMV,
                        drli,
                        nearestmv,
                        nearmv,
                        ref_mv);

                        //NEW_NEARMV
                        const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[context_ptr->me_sb_addr];

                        int16_t to_inject_mv_x_l0 = me_results->me_mv_array[context_ptr->me_block_offset][ref_idx_0].x_mv << 1;
                        int16_t to_inject_mv_y_l0 = me_results->me_mv_array[context_ptr->me_block_offset][ref_idx_0].y_mv << 1;
                        int16_t to_inject_mv_x_l1 = nearmv[1].as_mv.col;
                        int16_t to_inject_mv_y_l1 = nearmv[1].as_mv.row;

                        inj_mv = context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, ref_pair) == EB_FALSE;

                        if (inj_mv){
                            context_ptr->variance_ready = 0 ;
                            for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++){
                                // If two predictors are very similar, skip wedge compound mode search
                                if (context_ptr->variance_ready)
                                    if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEW_NEARMV) && context_ptr->prediction_mse < 64))
                                        continue;

                                candidateArray[canIdx].type = INTER_MODE;
                                candidateArray[canIdx].inter_mode = NEW_NEARMV;
                                candidateArray[canIdx].pred_mode = NEW_NEARMV;
                                candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                                candidateArray[canIdx].is_compound = 1;
                                candidateArray[canIdx].is_interintra_used = 0;

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

                                candidateArray[canIdx].motion_vector_pred_x[REF_LIST_0] = ref_mv[0].as_mv.col;
                                candidateArray[canIdx].motion_vector_pred_y[REF_LIST_0] = ref_mv[0].as_mv.row;

                                context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                                context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                                context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                                context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                                context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = ref_pair;
                                ++context_ptr->injected_mv_count_bipred;

                                //NEW_NEARMV
                                determine_compound_mode(
                                    picture_control_set_ptr,
                                    context_ptr,
                                    &candidateArray[canIdx],
                                    cur_type);

                                INCRMENT_CAND_TOTAL_COUNT(canIdx);
                        }
                    }
                }
            }
            //NEAR_NEWMV
            {
               maxDrlIndex = GetMaxDrlIndex(xd->ref_mv_count[ref_pair], NEAR_NEWMV);

               for (drli = 0; drli < maxDrlIndex; drli++) {
                   get_av1_mv_pred_drl(
                       context_ptr,
                       context_ptr->cu_ptr,
                       ref_pair,
                       1,
                       NEAR_NEWMV,
                       drli,
                       nearestmv,
                       nearmv,
                       ref_mv);

                   //NEAR_NEWMV
                   const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[context_ptr->me_sb_addr];

                   int16_t to_inject_mv_x_l0 = nearmv[0].as_mv.col;
                   int16_t to_inject_mv_y_l0 = nearmv[0].as_mv.row;
                   int16_t to_inject_mv_x_l1 = me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (get_list_idx(rf[1]) << 2) : (get_list_idx(rf[1]) << 1)) + ref_idx_1].x_mv << 1;//context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].comp_mv.as_mv.col;
                   int16_t to_inject_mv_y_l1 = me_results->me_mv_array[context_ptr->me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (get_list_idx(rf[1]) << 2) : (get_list_idx(rf[1]) << 1)) + ref_idx_1].y_mv << 1;//context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[ref_pair][0].comp_mv.as_mv.row;

                   inj_mv = context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, ref_pair) == EB_FALSE;

                   if (inj_mv) {

                       context_ptr->variance_ready = 0 ;
                       for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types ; cur_type++){
                           // If two predictors are very similar, skip wedge compound mode search
                           if (context_ptr->variance_ready)
                               if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEAR_NEWMV) && context_ptr->prediction_mse  < 64))
                                   continue;

                           candidateArray[canIdx].type = INTER_MODE;
                           candidateArray[canIdx].inter_mode = NEAR_NEWMV;
                           candidateArray[canIdx].pred_mode = NEAR_NEWMV;
                           candidateArray[canIdx].motion_mode = SIMPLE_TRANSLATION;
                           candidateArray[canIdx].is_compound = 1;
                           candidateArray[canIdx].is_interintra_used = 0;

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

                           candidateArray[canIdx].motion_vector_pred_x[REF_LIST_1] = ref_mv[1].as_mv.col;
                           candidateArray[canIdx].motion_vector_pred_y[REF_LIST_1] = ref_mv[1].as_mv.row;

                           context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                           context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                           context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                           context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                           context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = ref_pair;
                           ++context_ptr->injected_mv_count_bipred;

                           //NEAR_NEWMV
                           determine_compound_mode(
                               picture_control_set_ptr,
                               context_ptr,
                               &candidateArray[canIdx],
                               cur_type);

                           INCRMENT_CAND_TOTAL_COUNT(canIdx);
                        }
                   }
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
    MeLcuResults                   *meResult)
{
    uint32_t canIdx = *candTotCnt;
    ModeDecisionCandidate *candidateArray = context_ptr->fast_candidate_array;
    MacroBlockD  *xd = cu_ptr->av1xd;
    uint8_t drli, maxDrlIndex;
    IntMv nearestmv[2], nearmv[2], ref_mv[2];

    int inside_tile = 1;
    SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet *)picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    int umv0tile = (sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0);
    uint32_t mi_row = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t mi_col = context_ptr->cu_origin_x >> MI_SIZE_LOG2;

    if(umv0tile)
        inside_tile = is_inside_tile_boundary(&(xd->tile), context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.col, context_ptr->cu_ptr->ref_mvs[LAST_FRAME][0].as_mv.row, mi_col, mi_row, context_ptr->blk_geom->bsize);
    if(inside_tile)
    {
    //NEAREST_L0
        candidateArray[canIdx].type = INTER_MODE;
        candidateArray[canIdx].inter_mode = NEARESTMV;
        candidateArray[canIdx].pred_mode = NEARESTMV;
        candidateArray[canIdx].motion_mode = WARPED_CAUSAL;
        candidateArray[canIdx].wm_params_l0.wmtype = AFFINE;
        candidateArray[canIdx].is_compound = 0;
#if II_COMP_FLAG
        candidateArray[canIdx].is_interintra_used = 0;
#endif
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
            &candidateArray[canIdx].wm_params_l0,
            &candidateArray[canIdx].num_proj_ref);

        if (candidateArray[canIdx].local_warp_valid)
            INCRMENT_CAND_TOTAL_COUNT(canIdx);
    }

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
            candidateArray[canIdx].wm_params_l0.wmtype = AFFINE;
            candidateArray[canIdx].is_compound = 0;
#if II_COMP_FLAG
            candidateArray[canIdx].is_interintra_used = 0;
#endif
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
            if(umv0tile)
                inside_tile = is_inside_tile_boundary(&(xd->tile), mv_0.x, mv_0.y, mi_col, mi_row, context_ptr->blk_geom->bsize);
            if(inside_tile)
            {
            candidateArray[canIdx].local_warp_valid = warped_motion_parameters(
                picture_control_set_ptr,
                context_ptr->cu_ptr,
                &mv_unit,
                context_ptr->blk_geom,
                context_ptr->cu_origin_x,
                context_ptr->cu_origin_y,
                candidateArray[canIdx].ref_frame_type,
                &candidateArray[canIdx].wm_params_l0,
                &candidateArray[canIdx].num_proj_ref);

            if (candidateArray[canIdx].local_warp_valid)
                INCRMENT_CAND_TOTAL_COUNT(canIdx);
            }
    }

    // NEWMV L0
    const MV neighbors[9] = { { 0, 0 },
        { 0, -1 }, { 1, 0 }, { 0, 1 }, { -1, 0 } ,
        { 0, -2 }, { 2, 0 }, { 0, 2 }, { -2, 0 } };
    IntMv  bestPredmv[2] = { {0}, {0} };

    uint8_t total_me_cnt = meResult->total_me_candidate_index[context_ptr->me_block_offset];
    const MeCandidate *me_block_results = meResult->me_candidate[context_ptr->me_block_offset];
    //const MeLcuResults_t *meResults = pictureControlSetPtr->ParentPcsPtr->meResultsPtr[lcuAddr];
    for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index)
    {
        const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
        const uint8_t inter_direction = me_block_results_ptr->direction;
        const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;

#if MULTI_PASS_PD
        if (list0_ref_index > context_ptr->md_max_ref_count - 1)
            continue;
#endif

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
        candidateArray[canIdx].wm_params_l0.wmtype = AFFINE;

        candidateArray[canIdx].is_compound = 0;
#if II_COMP_FLAG
        candidateArray[canIdx].is_interintra_used = 0;
#endif
        candidateArray[canIdx].is_new_mv = 1;
        candidateArray[canIdx].is_zero_mv = 0;

        candidateArray[canIdx].drl_index = 0;

        // Set the MV to ME result
        candidateArray[canIdx].motion_vector_xl0 =  meResult->me_mv_array[context_ptr->me_block_offset][list0_ref_index].x_mv << 1;
        candidateArray[canIdx].motion_vector_yl0 =  meResult->me_mv_array[context_ptr->me_block_offset][list0_ref_index].y_mv << 1;
        if (picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv) {
            candidateArray[canIdx].motion_vector_xl0 += (neighbors[i].col);
            candidateArray[canIdx].motion_vector_yl0 += (neighbors[i].row);
        }
        else {
            candidateArray[canIdx].motion_vector_xl0 += (neighbors[i].col << 1);
            candidateArray[canIdx].motion_vector_yl0 += (neighbors[i].row << 1);
        }
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
        if(umv0tile)
            inside_tile = is_inside_tile_boundary(&(xd->tile), mv_0.x, mv_0.y, mi_col, mi_row, context_ptr->blk_geom->bsize);
        if(inside_tile)
        {
        candidateArray[canIdx].local_warp_valid = warped_motion_parameters(
            picture_control_set_ptr,
            context_ptr->cu_ptr,
            &mv_unit,
            context_ptr->blk_geom,
            context_ptr->cu_origin_x,
            context_ptr->cu_origin_y,
            candidateArray[canIdx].ref_frame_type,
            &candidateArray[canIdx].wm_params_l0,
            &candidateArray[canIdx].num_proj_ref);

        if (candidateArray[canIdx].local_warp_valid)
            INCRMENT_CAND_TOTAL_COUNT(canIdx);
        }
    }
        }
    }

    *candTotCnt = canIdx;
}



#if OBMC_FLAG

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
void eb_av1_setup_pred_block(BlockSize sb_type,
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

extern aom_variance_fn_ptr_t mefn_ptr[BlockSizeS_ALL];

int av1_find_best_obmc_sub_pixel_tree_up(
    ModeDecisionContext *context_ptr,IntraBcContext *x, const AV1_COMMON *const cm, int mi_row, int mi_col,
    MV *bestmv, const MV *ref_mv, int allow_hp, int error_per_bit,
    const aom_variance_fn_ptr_t *vfp, int forced_stop, int iters_per_step,
    int *mvjcost, int *mvcost[2], int *distortion, unsigned int *sse1,
    int is_second, int use_accurate_subpel_search)  ;

int av1_obmc_full_pixel_search(
    ModeDecisionContext *context_ptr,
    IntraBcContext *x,
    MV *mvp_full,
    int sadpb,
    const aom_variance_fn_ptr_t *fn_ptr,
    const MV *ref_mv,
    MV *dst_mv,
    int is_second);

static void single_motion_search(
    PictureControlSet *pcs,
    ModeDecisionContext *context_ptr,
    ModeDecisionCandidate        *candidate_ptr,
    MvReferenceFrame *rf,
    IntMv  bestPredmv,
    IntraBcContext  *x,
    BlockSize bsize,
    MV*      ref_mv,
    int ref_idx,
    int *rate_mv) {

    (void)ref_idx;
    const Av1Common *const cm = pcs->parent_pcs_ptr->av1_cm;
    FrameHeader *frm_hdr = &pcs->parent_pcs_ptr->frm_hdr;

    x->xd = context_ptr->cu_ptr->av1xd;
    const int mi_row = -x->xd->mb_to_top_edge / (8 * MI_SIZE);
    const int mi_col = -x->xd->mb_to_left_edge / (8 * MI_SIZE);

    x->nmv_vec_cost = context_ptr->md_rate_estimation_ptr->nmv_vec_cost;
    x->mv_cost_stack = context_ptr->md_rate_estimation_ptr->nmvcoststack;
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
    x->sadperbit16 = sad_per_bit16lut_8[frm_hdr->quantization_params.base_q_idx];
    x->errorperbit = context_ptr->full_lambda >> RD_EPB_SHIFT;
    x->errorperbit += (x->errorperbit == 0);


    int bestsme = INT_MAX;
    int sadpb = x->sadperbit16;
    MV mvp_full;

    MvLimits tmp_mv_limits = x->mv_limits;

    // Note: MV limits are modified here. Always restore the original values
    // after full-pixel motion search.
    eb_av1_set_mv_search_range(&x->mv_limits, ref_mv);

    mvp_full = bestPredmv.as_mv; // mbmi->mv[0].as_mv;


    mvp_full.col >>= 3;
    mvp_full.row >>= 3;

    x->best_mv.as_int = x->second_best_mv.as_int = INVALID_MV; //D

    switch (candidate_ptr->motion_mode) {

    case OBMC_CAUSAL:
        bestsme = av1_obmc_full_pixel_search(
            context_ptr,
            x,
            &mvp_full,
            sadpb,
            &mefn_ptr[bsize],
            ref_mv,
            &(x->best_mv.as_mv),
            0);
        break;
    default: assert(0 && "Invalid motion mode!\n");
    }

    x->mv_limits = tmp_mv_limits;

    const int use_fractional_mv =
        bestsme < INT_MAX && frm_hdr->force_integer_mv == 0;
    if (use_fractional_mv) {
        int dis; /* TODO: use dis in distortion calculation later. */
        switch (candidate_ptr->motion_mode) {
        case OBMC_CAUSAL:
            av1_find_best_obmc_sub_pixel_tree_up(
                context_ptr,
                x,
                cm,
                mi_row,
                mi_col,
                &x->best_mv.as_mv,
                ref_mv,
                frm_hdr->allow_high_precision_mv,
                x->errorperbit,
                &mefn_ptr[bsize],
                0, // mv.subpel_force_stop
                2, //  mv.subpel_iters_per_step
                x->nmv_vec_cost,
                x->mv_cost_stack,
                &dis,
                &context_ptr->pred_sse[rf[0]],
                0,
                USE_8_TAPS );
            break;
        default: assert(0 && "Invalid motion mode!\n");
        }
    }
    *rate_mv = eb_av1_mv_bit_cost(&x->best_mv.as_mv, ref_mv, x->nmv_vec_cost,
        x->mv_cost_stack, MV_COST_WEIGHT);

}

void obmc_motion_refinement(
    PictureControlSet          *picture_control_set_ptr,
    struct ModeDecisionContext *context_ptr,
    ModeDecisionCandidate    *candidate,
    uint8_t                   ref_list_idx)
{
    IntMv  bestPredmv[2] = { {0}, {0} };
    IntraBcContext  x_st;
    IntraBcContext  *x = &x_st;

    MacroBlockD * xd;
    xd = x->xd = context_ptr->cu_ptr->av1xd;
    const int mi_row = -xd->mb_to_top_edge / (8 * MI_SIZE);
    const int mi_col = -xd->mb_to_left_edge / (8 * MI_SIZE);

    {
        uint8_t ref_idx = get_ref_frame_idx(candidate->ref_frame_type);
        uint8_t list_idx = get_list_idx(candidate->ref_frame_type);
        EbPictureBufferDesc  *reference_picture = ((EbReferenceObject*)picture_control_set_ptr->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr)->reference_picture;
        Yv12BufferConfig ref_buf;
        link_Eb_to_aom_buffer_desc_8bit(
            reference_picture,
            &ref_buf);

        struct Buf2D yv12_mb[MAX_MB_PLANE];
        eb_av1_setup_pred_block(context_ptr->blk_geom->bsize, yv12_mb, &ref_buf, mi_row, mi_col);
        for (int i = 0; i < 1; ++i)
            x->xdplane[i].pre[0] = yv12_mb[i];  //ref in ME

        x->plane[0].src.buf = 0;// x->xdplane[0].pre[0];
        x->plane[0].src.buf0 = 0;
    }

    IntMv  best_mv ;
    best_mv.as_int = 0;
    if (ref_list_idx == 0) {
        best_mv.as_mv.col = candidate->motion_vector_xl0;// to_inject_mv_x;
        best_mv.as_mv.row = candidate->motion_vector_yl0; //to_inject_mv_y;
    }
    else
    {
        best_mv.as_mv.col = candidate->motion_vector_xl1;// to_inject_mv_x;
        best_mv.as_mv.row = candidate->motion_vector_yl1; //to_inject_mv_y;
    }
    int tmp_rate_mv;

    MvReferenceFrame rf[2];
    rf[0] = candidate->ref_frame_type;
    rf[1] = -1;

    MV      ref_mv;
    ref_mv.col = candidate->motion_vector_pred_x[ref_list_idx];
    ref_mv.row = candidate->motion_vector_pred_y[ref_list_idx];

    single_motion_search(
        picture_control_set_ptr,
        context_ptr,
        candidate,
        rf,
        best_mv,
        x,
        context_ptr->blk_geom->bsize,
        &ref_mv,
        0,
        &tmp_rate_mv);

    if (ref_list_idx == 0) {
        candidate->motion_vector_xl0 = x->best_mv.as_mv.col;
        candidate->motion_vector_yl0 = x->best_mv.as_mv.row;
    }
    else {
        candidate->motion_vector_xl1 = x->best_mv.as_mv.col;
        candidate->motion_vector_yl1 = x->best_mv.as_mv.row;
    }

    ChooseBestAv1MvPred(
        context_ptr,
        candidate->md_rate_estimation_ptr,
        context_ptr->cu_ptr,
        candidate->ref_frame_type,
        candidate->is_compound,
        candidate->pred_mode,
        ref_list_idx == 0 ? candidate->motion_vector_xl0 : candidate->motion_vector_xl1,
        ref_list_idx == 0 ? candidate->motion_vector_yl0 : candidate->motion_vector_yl1,
        0, 0,
        &candidate->drl_index,
        bestPredmv);

    if (ref_list_idx == 0) {
        candidate->motion_vector_pred_x[REF_LIST_0] = bestPredmv[0].as_mv.col;
        candidate->motion_vector_pred_y[REF_LIST_0] = bestPredmv[0].as_mv.row;
    }
    else {
        candidate->motion_vector_pred_x[REF_LIST_1] = bestPredmv[0].as_mv.col;
        candidate->motion_vector_pred_y[REF_LIST_1] = bestPredmv[0].as_mv.row;
    }
}
#endif

void inject_new_candidates(
    const SequenceControlSet   *sequence_control_set_ptr,
    struct ModeDecisionContext *context_ptr,
    PictureControlSet          *picture_control_set_ptr,
    EbBool                      isCompoundEnabled,
    EbBool                      allow_bipred,
    uint32_t                    me_sb_addr,
    uint32_t                    me_block_offset,
    uint32_t                   *candidateTotalCnt) {

    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    IntMv  bestPredmv[2] = { {0}, {0} };
    uint32_t canTotalCnt = (*candidateTotalCnt);

    const MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[me_sb_addr];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[me_block_offset];
    const MeCandidate *me_block_results = me_results->me_candidate[me_block_offset];
    MacroBlockD  *xd = context_ptr->cu_ptr->av1xd;
    int inside_tile = 1;
    int umv0tile = (sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0);
    uint32_t mi_row = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t mi_col = context_ptr->cu_origin_x >> MI_SIZE_LOG2;
    BlockSize bsize = context_ptr->blk_geom->bsize;                       // bloc size
    MD_COMP_TYPE cur_type; //NN
    MD_COMP_TYPE tot_comp_types =
        (picture_control_set_ptr->parent_pcs_ptr->compound_mode == 1 || context_ptr->compound_types_to_try == MD_COMP_AVG) ?
            MD_COMP_AVG :
            (bsize >= BLOCK_8X8 && bsize <= BLOCK_32X32) ?
                context_ptr->compound_types_to_try :
                context_ptr->compound_types_to_try == MD_COMP_WEDGE ?
                     MD_COMP_DIFF0 :
                     context_ptr->compound_types_to_try;
    if (context_ptr->source_variance < context_ptr->inter_inter_wedge_variance_th)
        tot_comp_types = MIN(tot_comp_types, MD_COMP_DIFF0);
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
#if MULTI_PASS_PD
            if (list0_ref_index > context_ptr->md_max_ref_count - 1)
                continue;
#endif
            int16_t to_inject_mv_x = me_results->me_mv_array[me_block_offset][list0_ref_index].x_mv << 1;
            int16_t to_inject_mv_y = me_results->me_mv_array[me_block_offset][list0_ref_index].y_mv << 1;
            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
            uint8_t skip_cand = check_ref_beackout(
                context_ptr,
                to_inject_ref_type,
                context_ptr->blk_geom->shape);

            inside_tile = 1;
            if(umv0tile)
                inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x, to_inject_mv_y, mi_col, mi_row, context_ptr->blk_geom->bsize);
            skip_cand = skip_cand || (!inside_tile);

            if (!skip_cand && (context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE)) {

#if II_COMP_FLAG    // NEWMV L0
             MvReferenceFrame rf[2];
             rf[0] = to_inject_ref_type;
             rf[1] = -1;

            uint8_t inter_type;
#if MULTI_PASS_PD
            uint8_t is_ii_allowed = svt_is_interintra_allowed(context_ptr->md_enable_inter_intra, bsize, NEWMV, rf);
#else
            uint8_t is_ii_allowed = svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, NEWMV, rf);
#endif
            uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
#if OBMC_FLAG
#if MULTI_PASS_PD
            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
            tot_inter_types = is_obmc_allowed && context_ptr->md_pic_obmc_mode <= 2 ? tot_inter_types + 1 : tot_inter_types;
#else
             uint8_t is_obmc_allowed =  obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
             tot_inter_types = is_obmc_allowed && picture_control_set_ptr->parent_pcs_ptr->pic_obmc_mode <= 2 ? tot_inter_types + 1 : tot_inter_types;
#endif
#endif
            for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
            {
#endif
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

#if II_COMP_FLAG
            if (inter_type == 0) {
                candidateArray[canTotalCnt].is_interintra_used = 0;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            }
            else {
                if (is_ii_allowed) {
                    if (inter_type == 1) {
                        inter_intra_search(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canTotalCnt]);
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].use_wedge_interintra = 1;
                        candidateArray[canTotalCnt].ii_wedge_sign = 0;
                    }
                    else if (inter_type == 2) {
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].interintra_mode = candidateArray[canTotalCnt - 1].interintra_mode;
                        candidateArray[canTotalCnt].use_wedge_interintra = 0;
                    }
                }
 #if OBMC_FLAG
                if (is_obmc_allowed && inter_type == tot_inter_types-1){
                        candidateArray[canTotalCnt].is_interintra_used = 0;
                        candidateArray[canTotalCnt].motion_mode = OBMC_CAUSAL;

                        obmc_motion_refinement(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canTotalCnt],
                            REF_LIST_0);
                    }
#endif
            }
#endif
                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

#if II_COMP_FLAG
            }
#endif
                context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
                context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
                context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = to_inject_ref_type;
                ++context_ptr->injected_mv_count_l0;
#if MULTI_PASS_PD // Test 1 inter only if 1st pass
                if (context_ptr->best_me_cand_only_flag)
                    break;
#endif
            }

            }

        if (isCompoundEnabled) {
            /**************
               NEWMV L1
           ************* */
            if (inter_direction == 1) {
#if MULTI_PASS_PD
                if (list1_ref_index > context_ptr->md_max_ref_count - 1)
                    continue;
#endif
                int16_t to_inject_mv_x =  me_results->me_mv_array[me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].x_mv << 1;
                int16_t to_inject_mv_y =  me_results->me_mv_array[me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? 4 : 2) + list1_ref_index].y_mv << 1;
                uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                uint8_t skip_cand = check_ref_beackout(
                    context_ptr,
                    to_inject_ref_type,
                    context_ptr->blk_geom->shape);

                inside_tile = 1;
                if(umv0tile)
                    inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x, to_inject_mv_y, mi_col, mi_row, context_ptr->blk_geom->bsize);
                skip_cand = skip_cand || !inside_tile;

                if (!skip_cand && (context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE)) {
  #if II_COMP_FLAG // NEWMV L1
             MvReferenceFrame rf[2];
             rf[0] = to_inject_ref_type;
             rf[1] = -1;

            uint8_t inter_type;
#if MULTI_PASS_PD
            uint8_t is_ii_allowed = svt_is_interintra_allowed(context_ptr->md_enable_inter_intra, bsize, NEWMV, rf);
#else
            uint8_t is_ii_allowed = svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, NEWMV, rf);
#endif
            uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
#if OBMC_FLAG
#if MULTI_PASS_PD
            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
#else
            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
#endif
            tot_inter_types = is_obmc_allowed  && picture_control_set_ptr->parent_pcs_ptr->pic_obmc_mode <= 2 ? tot_inter_types + 1 : tot_inter_types;
#endif
            for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
            {
#endif
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
#if II_COMP_FLAG
            if (inter_type == 0) {
                candidateArray[canTotalCnt].is_interintra_used = 0;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            }
            else {
                if (is_ii_allowed) {
                    if (inter_type == 1) {
                        inter_intra_search(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canTotalCnt]);
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].use_wedge_interintra = 1;
                        candidateArray[canTotalCnt].ii_wedge_sign = 0;
                    }
                    else if (inter_type == 2) {
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].interintra_mode = candidateArray[canTotalCnt - 1].interintra_mode;
                        candidateArray[canTotalCnt].use_wedge_interintra = 0;
                    }
                }
#if OBMC_FLAG
                if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                    candidateArray[canTotalCnt].is_interintra_used = 0;
                    candidateArray[canTotalCnt].motion_mode = OBMC_CAUSAL;

                    obmc_motion_refinement(
                        picture_control_set_ptr,
                        context_ptr,
                        &candidateArray[canTotalCnt],
                        REF_LIST_1);
                }
#endif

            }
#endif
                    INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

#if II_COMP_FLAG
            }
#endif
                    context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
                    context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
                    context_ptr->injected_ref_type_l1_array[context_ptr->injected_mv_count_l1] = to_inject_ref_type;
                    ++context_ptr->injected_mv_count_l1;
#if MULTI_PASS_PD // Test 1 inter only if 1st pass
                    if (context_ptr->best_me_cand_only_flag)
                        break;
#endif
                }

                }
            /**************
               NEW_NEWMV
            ************* */
            if (allow_bipred) {
#if MULTI_PASS_PD
                if (list0_ref_index > context_ptr->md_max_ref_count - 1 || list1_ref_index > context_ptr->md_max_ref_count - 1)
                    continue;
#endif
                if (inter_direction == 2) {
                    int16_t to_inject_mv_x_l0 = me_results->me_mv_array[me_block_offset][list0_ref_index].x_mv << 1;
                    int16_t to_inject_mv_y_l0 = me_results->me_mv_array[me_block_offset][list0_ref_index].y_mv << 1;
                    int16_t to_inject_mv_x_l1 = me_results->me_mv_array[me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].x_mv << 1;
                    int16_t to_inject_mv_y_l1 = me_results->me_mv_array[me_block_offset][((sequence_control_set_ptr->mrp_mode == 0) ? (me_block_results_ptr->ref1_list << 2) : (me_block_results_ptr->ref1_list << 1)) + list1_ref_index].y_mv << 1;
                    MvReferenceFrame rf[2];
                    rf[0] = svt_get_ref_frame_type(me_block_results_ptr->ref0_list, list0_ref_index);
                    rf[1] = svt_get_ref_frame_type(me_block_results_ptr->ref1_list, list1_ref_index);
                    uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
                    uint8_t skip_cand = check_ref_beackout(
                        context_ptr,
                        to_inject_ref_type,
                        context_ptr->blk_geom->shape);

                    inside_tile = 1;
                    if(umv0tile) {
                        inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l0, mi_col, mi_row, context_ptr->blk_geom->bsize) &&
                                      is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l1, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize);
                    }
                    skip_cand = skip_cand || (!inside_tile);
                    if (!skip_cand && (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE)) {
                        context_ptr->variance_ready = 0;
                        for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
                        {
                            if (cur_type == MD_COMP_WEDGE && wedge_params_lookup[context_ptr->blk_geom->bsize].bits == 0) continue;
                            // If two predictors are very similar, skip wedge compound mode search
                            if (context_ptr->variance_ready)
                                if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEW_NEWMV) && context_ptr->prediction_mse < 64))
                                    continue;
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
#if II_COMP_FLAG
                        candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
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
                        //NEW_NEW
                        determine_compound_mode(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canTotalCnt],
                            cur_type);
                        INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

                        context_ptr->injected_mv_x_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l0;
                        context_ptr->injected_mv_y_bipred_l0_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l0;
                        context_ptr->injected_mv_x_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_x_l1;
                        context_ptr->injected_mv_y_bipred_l1_array[context_ptr->injected_mv_count_bipred] = to_inject_mv_y_l1;
                        context_ptr->injected_ref_type_bipred_array[context_ptr->injected_mv_count_bipred] = to_inject_ref_type;
                        ++context_ptr->injected_mv_count_bipred;
#if MULTI_PASS_PD // Test 1 inter only if 1st pass
                        if (context_ptr->best_me_cand_only_flag)
                            break;
#endif
                        }
                    }

                    }
                }
            }
            }
    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;
        }
        void inject_predictive_me_candidates(
            //const SequenceControlSet   *sequence_control_set_ptr,
            struct ModeDecisionContext *context_ptr,
            PictureControlSet          *picture_control_set_ptr,
            EbBool                      isCompoundEnabled,
            EbBool                      allow_bipred,
            uint32_t                   *candidateTotalCnt) {

            ModeDecisionCandidate *candidateArray = context_ptr->fast_candidate_array;
            IntMv  bestPredmv[2] = { {0}, {0} };
            uint32_t canTotalCnt = (*candidateTotalCnt);
            BlockSize bsize = context_ptr->blk_geom->bsize;                       // bloc size

            MD_COMP_TYPE cur_type; //BIP 3x3 MiSize >= BLOCK_8X8 && MiSize <= BLOCK_32X32)
            MD_COMP_TYPE tot_comp_types =
                (bsize >= BLOCK_8X8 && bsize <= BLOCK_32X32) ?
                    context_ptr->compound_types_to_try :
                        context_ptr->compound_types_to_try == MD_COMP_WEDGE ? MD_COMP_DIFF0 :
                        context_ptr->compound_types_to_try;
            if (context_ptr->source_variance < context_ptr->inter_inter_wedge_variance_th)
                tot_comp_types = MIN(tot_comp_types, MD_COMP_DIFF0);
            uint8_t listIndex;
            uint8_t ref_pic_index;
            listIndex = REF_LIST_0;
            {
                // Ref Picture Loop
                for (ref_pic_index = 0; ref_pic_index < 4; ++ref_pic_index) {
                    if (context_ptr->valid_refined_mv[listIndex][ref_pic_index]) {
                        int16_t to_inject_mv_x = context_ptr->best_spatial_pred_mv[listIndex][ref_pic_index][0];
                        int16_t to_inject_mv_y = context_ptr->best_spatial_pred_mv[listIndex][ref_pic_index][1];
                        uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, ref_pic_index);
                        if (context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {
#if OBMC_FLAG
                            MvReferenceFrame rf[2];
                            rf[0] = to_inject_ref_type;
                            rf[1] = -1;
                            uint8_t inter_type;
                            uint8_t is_ii_allowed = 0;// svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, NEWMV, rf);
                            uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
#if MULTI_PASS_PD
                            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
#else
                            uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
#endif
                            tot_inter_types = is_obmc_allowed ? tot_inter_types + 1 : tot_inter_types;
                            for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
                            {
#endif
                            candidateArray[canTotalCnt].type = INTER_MODE;
                            candidateArray[canTotalCnt].distortion_ready = 0;
                            candidateArray[canTotalCnt].use_intrabc = 0;
                            candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;
                            candidateArray[canTotalCnt].inter_mode = NEWMV;
                            candidateArray[canTotalCnt].pred_mode = NEWMV;
                            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                            candidateArray[canTotalCnt].is_compound = 0;
#if II_COMP_FLAG // PME OFF   L0
                            candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
                            candidateArray[canTotalCnt].is_new_mv = 1;
                            candidateArray[canTotalCnt].is_zero_mv = 0;
                            candidateArray[canTotalCnt].drl_index = 0;
                            candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x;
                            candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y;
                            candidateArray[canTotalCnt].ref_mv_index = 0;
                            candidateArray[canTotalCnt].pred_mv_weight = 0;
                            candidateArray[canTotalCnt].ref_frame_type = svt_get_ref_frame_type(REF_LIST_0, ref_pic_index);
                            candidateArray[canTotalCnt].ref_frame_index_l0 = ref_pic_index;
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
#if OBMC_FLAG
                            if (inter_type == 0) {
                                candidateArray[canTotalCnt].is_interintra_used = 0;
                                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                            }
                            else {
                                if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                                    candidateArray[canTotalCnt].is_interintra_used = 0;
                                    candidateArray[canTotalCnt].motion_mode = OBMC_CAUSAL;

                                    obmc_motion_refinement(
                                        picture_control_set_ptr,
                                        context_ptr,
                                        &candidateArray[canTotalCnt],
                                        REF_LIST_0);
                                }
                            }
#endif
                            INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                            context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x;
                            context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y;
                            context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = to_inject_ref_type;
                            ++context_ptr->injected_mv_count_l0;
#if OBMC_FLAG
                        }
#endif
                        }
                    }
                }
            }
            if (isCompoundEnabled) {
                /**************
                   NEWMV L1
               ************* */
                listIndex = REF_LIST_1;
                {
                    // Ref Picture Loop
                    for (ref_pic_index = 0; ref_pic_index < 3; ++ref_pic_index) {
                        if (context_ptr->valid_refined_mv[listIndex][ref_pic_index]) {
                            int16_t to_inject_mv_x = context_ptr->best_spatial_pred_mv[listIndex][ref_pic_index][0];
                            int16_t to_inject_mv_y = context_ptr->best_spatial_pred_mv[listIndex][ref_pic_index][1];
                            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_1, ref_pic_index);
                            if (context_ptr->injected_mv_count_l1 == 0 || mrp_is_already_injected_mv_l1(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE) {
#if OBMC_FLAG
                                MvReferenceFrame rf[2];
                                rf[0] = to_inject_ref_type;
                                rf[1] = -1;
                                uint8_t inter_type;
                                uint8_t is_ii_allowed = 0;// svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, NEWMV, rf);
                                uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
#if MULTI_PASS_PD
                                uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
#else
                                uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
#endif
                                tot_inter_types = is_obmc_allowed ? tot_inter_types + 1 : tot_inter_types;
                                for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
                                {
#endif
                                candidateArray[canTotalCnt].type = INTER_MODE;
                                candidateArray[canTotalCnt].distortion_ready = 0;
                                candidateArray[canTotalCnt].use_intrabc = 0;
                                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)1;
                                candidateArray[canTotalCnt].inter_mode = NEWMV;
                                candidateArray[canTotalCnt].pred_mode = NEWMV;
                                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                                candidateArray[canTotalCnt].is_compound = 0;
#if II_COMP_FLAG // PME OFF   L1
                                candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
                                candidateArray[canTotalCnt].is_new_mv = 1;
                                candidateArray[canTotalCnt].is_zero_mv = 0;
                                candidateArray[canTotalCnt].drl_index = 0;
                                candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x;
                                candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y;
                                candidateArray[canTotalCnt].ref_mv_index = 0;
                                candidateArray[canTotalCnt].pred_mv_weight = 0;
                                candidateArray[canTotalCnt].ref_frame_type = svt_get_ref_frame_type(REF_LIST_1, ref_pic_index);
                                candidateArray[canTotalCnt].ref_frame_index_l0 = -1;
                                candidateArray[canTotalCnt].ref_frame_index_l1 = ref_pic_index;
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
#if OBMC_FLAG
                                if (inter_type == 0) {
                                    candidateArray[canTotalCnt].is_interintra_used = 0;
                                    candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                                }
                                else {
                                    if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                                        candidateArray[canTotalCnt].is_interintra_used = 0;
                                        candidateArray[canTotalCnt].motion_mode = OBMC_CAUSAL;

                                        obmc_motion_refinement(
                                            picture_control_set_ptr,
                                            context_ptr,
                                            &candidateArray[canTotalCnt],
                                            REF_LIST_1);
                                    }
                                }
#endif
                                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                                context_ptr->injected_mv_x_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_x;
                                context_ptr->injected_mv_y_l1_array[context_ptr->injected_mv_count_l1] = to_inject_mv_y;
                                context_ptr->injected_ref_type_l1_array[context_ptr->injected_mv_count_l1] = to_inject_ref_type;
                                ++context_ptr->injected_mv_count_l1;
#if OBMC_FLAG
                            }
#endif
                            }
                        }
                    }
                }
            }
            /**************
                NEW_NEWMV
            ************* */
            if (allow_bipred) {
                uint8_t ref_pic_index_l0;
                uint8_t ref_pic_index_l1;
                {
                    // Ref Picture Loop
                    for (ref_pic_index_l0 = 0; ref_pic_index_l0 < 4; ++ref_pic_index_l0) {
                        for (ref_pic_index_l1 = 0; ref_pic_index_l1 < 4; ++ref_pic_index_l1) {
                            if (context_ptr->valid_refined_mv[REF_LIST_0][ref_pic_index_l0] && context_ptr->valid_refined_mv[REF_LIST_1][ref_pic_index_l1]) {
                                int16_t to_inject_mv_x_l0 = context_ptr->best_spatial_pred_mv[REF_LIST_0][ref_pic_index_l0][0];
                                int16_t to_inject_mv_y_l0 = context_ptr->best_spatial_pred_mv[REF_LIST_0][ref_pic_index_l0][1];
                                int16_t to_inject_mv_x_l1 = context_ptr->best_spatial_pred_mv[REF_LIST_1][ref_pic_index_l1][0];
                                int16_t to_inject_mv_y_l1 = context_ptr->best_spatial_pred_mv[REF_LIST_1][ref_pic_index_l1][1];

                                MvReferenceFrame rf[2];
                                rf[0] = svt_get_ref_frame_type(REF_LIST_0, ref_pic_index_l0);
                                rf[1] = svt_get_ref_frame_type(REF_LIST_1, ref_pic_index_l1);
                                uint8_t to_inject_ref_type = av1_ref_frame_type(rf);
                                if (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE) {

                                    context_ptr->variance_ready = 0;
                                    for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
                                    {
                                        // If two predictors are very similar, skip wedge compound mode search
                                        if (context_ptr->variance_ready)
                                            if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(NEW_NEWMV) && context_ptr->prediction_mse < 64))
                                                continue;

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
#if II_COMP_FLAG
                                        candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
                                        candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;

                                        MvReferenceFrame rf[2];
                                        rf[0] = svt_get_ref_frame_type(REF_LIST_0, ref_pic_index_l0);
                                        rf[1] = svt_get_ref_frame_type(REF_LIST_1, ref_pic_index_l1);
                                        candidateArray[canTotalCnt].ref_frame_type = av1_ref_frame_type(rf);
                                        candidateArray[canTotalCnt].ref_frame_index_l0 = ref_pic_index_l0;
                                        candidateArray[canTotalCnt].ref_frame_index_l1 = ref_pic_index_l1;

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

                                        //MVP REFINE
                                        determine_compound_mode(
                                            picture_control_set_ptr,
                                            context_ptr,
                                            &candidateArray[canTotalCnt],
                                            cur_type);
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
                }
            }

            (*candidateTotalCnt) = canTotalCnt;
        }

void  inject_inter_candidates(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    const SequenceControlSet     *sequence_control_set_ptr,
    SuperBlock                   *sb_ptr,
#if ENHANCED_M0_SETTINGS
    EbBool                        coeff_based_nsq_cand_reduction,
#endif
    uint32_t                       *candidateTotalCnt) {

    (void)sequence_control_set_ptr;

    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    uint32_t                   canTotalCnt = *candidateTotalCnt;
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool isCompoundEnabled = (frm_hdr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
#if GM_OPT
    uint8_t inj_mv = 1;
#endif
    int inside_tile = 1;
    MacroBlockD  *xd = context_ptr->cu_ptr->av1xd;
    int umv0tile = (sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0);
    MeLcuResults *me_results = picture_control_set_ptr->parent_pcs_ptr->me_results[context_ptr->me_sb_addr];
    EbBool allow_bipred = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : EB_TRUE;
#if !ENHANCED_M0_SETTINGS
    uint8_t sq_index = LOG2F(context_ptr->blk_geom->sq_size) - 2;
    uint8_t inject_newmv_candidate = 1;
#if MULTI_PASS_PD // Shut coef-based inter skip if 1st pass
    if (context_ptr->coeff_based_nsq_cand_reduction) {
#else
    if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level >= NSQ_SEARCH_LEVEL1 &&
        picture_control_set_ptr->parent_pcs_ptr->nsq_search_level < NSQ_SEARCH_FULL) {
#endif
        if (context_ptr->md_local_cu_unit[context_ptr->blk_geom->sqi_mds].avail_blk_flag)
            inject_newmv_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
                context_ptr->parent_sq_has_coeff[sq_index] != 0 ? inject_newmv_candidate : 0;
    }
#endif
    BlockSize bsize = context_ptr->blk_geom->bsize;                       // bloc size
#if MULTI_PASS_PD

#else

#endif
    MD_COMP_TYPE cur_type; //GG

    MD_COMP_TYPE tot_comp_types =
        (picture_control_set_ptr->parent_pcs_ptr->compound_mode == 1 || context_ptr->compound_types_to_try == MD_COMP_AVG) ?
            MD_COMP_AVG :
            (bsize >= BLOCK_8X8 && bsize <= BLOCK_32X32) ?
            context_ptr->compound_types_to_try :
            context_ptr->compound_types_to_try == MD_COMP_WEDGE ?
                MD_COMP_DIFF0 :
                context_ptr->compound_types_to_try;

    if (context_ptr->source_variance < context_ptr->inter_inter_wedge_variance_th)
        tot_comp_types = MIN(tot_comp_types, MD_COMP_DIFF0);

    uint32_t mi_row = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t mi_col = context_ptr->cu_origin_x >> MI_SIZE_LOG2;
    eb_av1_count_overlappable_neighbors(
        picture_control_set_ptr,
        context_ptr->cu_ptr,
        context_ptr->blk_geom->bsize,
        mi_row,
        mi_col);
    uint8_t is_obmc_allowed = obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr, context_ptr->blk_geom->bsize, LAST_FRAME, -1, NEWMV) == OBMC_CAUSAL;
    if (is_obmc_allowed)
        precompute_obmc_data(
            picture_control_set_ptr,
            context_ptr);
    /**************
         MVP
    ************* */

    uint32_t refIt;
#if MULTI_PASS_PD
    if (context_ptr->new_nearest_injection)
#endif
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
        EbBool allow_compound = (frm_hdr->reference_mode == SINGLE_REFERENCE || context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) ? EB_FALSE : EB_TRUE;
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
#if !ENHANCED_M0_SETTINGS
    if (inject_newmv_candidate) {
#endif

        inject_new_candidates(
            sequence_control_set_ptr,
            context_ptr,
            picture_control_set_ptr,
            isCompoundEnabled,
            allow_bipred,
            context_ptr->me_sb_addr,
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
                uint32_t x_to_search = context_ptr->blk_geom->origin_x - (context_ptr->geom_offset_x + ((context_ptr->blk_geom->origin_x & 0x7) ? 4 : 0));
                uint32_t y_to_search = context_ptr->blk_geom->origin_y - (context_ptr->geom_offset_y + ((context_ptr->blk_geom->origin_y & 0x7) ? 4 : 0));

                // Search the me_info_index of the parent block
                uint32_t me_info_index = 0;
                for (uint32_t block_index = 0; block_index < picture_control_set_ptr->parent_pcs_ptr->max_number_of_pus_per_sb; block_index++) {

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
                    me_info_index,
                    &canTotalCnt);
            }
        }
#if !ENHANCED_M0_SETTINGS
    }
#endif
    if (context_ptr->global_mv_injection) {
#if GLOBAL_WARPED_MOTION
#if GM_OPT
        if (picture_control_set_ptr->parent_pcs_ptr->gm_level <= GM_DOWN) {
#endif
        for (unsigned list_ref_index_l0 = 0; list_ref_index_l0 < 1; ++list_ref_index_l0)
        for (unsigned list_ref_index_l1 = 0; list_ref_index_l1 < 1; ++list_ref_index_l1) {

            /**************
             GLOBALMV
            ************* */

            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, list_ref_index_l0);
            EbWarpedMotionParams *params_l0 = &picture_control_set_ptr->parent_pcs_ptr->global_motion[to_inject_ref_type];

            IntMv mv_l0 = gm_get_motion_vector_enc(
                params_l0,
                picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv,
                context_ptr->blk_geom->bsize,
                mi_col, mi_row,
                0 /* force_integer_mv */);

            int16_t to_inject_mv_x_l0 = mv_l0.as_mv.col;
            int16_t to_inject_mv_y_l0 = mv_l0.as_mv.row;

            if(umv0tile)
                inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l0,
                                                      mi_col, mi_row, context_ptr->blk_geom->bsize);

            if (inside_tile
                && (((params_l0->wmtype > TRANSLATION
                      && context_ptr->blk_geom->bwidth >= 8
                      && context_ptr->blk_geom->bheight >= 8)
                     || params_l0->wmtype <= TRANSLATION))) {

#if II_COMP_FLAG      // GLOBALMV L0
                 MvReferenceFrame rf[2];
                 rf[0] = to_inject_ref_type;
                 rf[1] = -1;

                uint8_t inter_type;
#if MULTI_PASS_PD
                uint8_t is_ii_allowed = svt_is_interintra_allowed(context_ptr->md_enable_inter_intra, bsize, GLOBALMV, rf);
#else
                uint8_t is_ii_allowed = svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, GLOBALMV, rf);
#endif
                uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
                //uint8_t is_obmc_allowed =  obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
                //tot_inter_types = is_obmc_allowed ? tot_inter_types+1 : tot_inter_types;

                for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
                {
#endif
                candidateArray[canTotalCnt].type = INTER_MODE;

                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;

                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)0;

                candidateArray[canTotalCnt].inter_mode = GLOBALMV;
                candidateArray[canTotalCnt].pred_mode = GLOBALMV;
                candidateArray[canTotalCnt].motion_mode = params_l0->wmtype > TRANSLATION ? WARPED_CAUSAL : SIMPLE_TRANSLATION;

                candidateArray[canTotalCnt].wm_params_l0 = *params_l0;

                candidateArray[canTotalCnt].is_compound = 0;
                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;
                candidateArray[canTotalCnt].merge_flag = EB_FALSE;
                candidateArray[canTotalCnt].prediction_direction[0] = UNI_PRED_LIST_0;
                candidateArray[canTotalCnt].is_new_mv = 0;
                candidateArray[canTotalCnt].is_zero_mv = 0;
                candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x_l0;
                candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y_l0;
                candidateArray[canTotalCnt].drl_index = 0;
                candidateArray[canTotalCnt].ref_mv_index = 0;
                candidateArray[canTotalCnt].pred_mv_weight = 0;
                candidateArray[canTotalCnt].ref_frame_type = av1_ref_frame_type(rf);
                candidateArray[canTotalCnt].ref_frame_index_l0 = 0;
                candidateArray[canTotalCnt].ref_frame_index_l1 = -1;
                candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;

#if II_COMP_FLAG
                if (inter_type == 0) {
                    candidateArray[canTotalCnt].is_interintra_used = 0;
                }
                else {
                    if (is_ii_allowed) {
                        if (inter_type == 1) {
                            inter_intra_search(
                                picture_control_set_ptr,
                                context_ptr,
                                &candidateArray[canTotalCnt]);
                            candidateArray[canTotalCnt].is_interintra_used = 1;
                            candidateArray[canTotalCnt].use_wedge_interintra = 1;
                            candidateArray[canTotalCnt].ii_wedge_sign = 0;
                        }
                        else if (inter_type == 2) {
                            candidateArray[canTotalCnt].is_interintra_used = 1;
                            candidateArray[canTotalCnt].interintra_mode = candidateArray[canTotalCnt - 1].interintra_mode;
                            candidateArray[canTotalCnt].use_wedge_interintra = 0;
                        }
                    }
                    //if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                    //    candidateArray[canTotalCnt].is_interintra_used = 0;
                    //    candidateArray[canTotalCnt].motion_mode = OBMC_CAUSAL;
                    //}
                }
#endif

                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

#if II_COMP_FLAG
                }
#endif

                context_ptr->injected_mv_x_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_x_l0;
                context_ptr->injected_mv_y_l0_array[context_ptr->injected_mv_count_l0] = to_inject_mv_y_l0;
                context_ptr->injected_ref_type_l0_array[context_ptr->injected_mv_count_l0] = to_inject_ref_type;
                ++context_ptr->injected_mv_count_l0;

                EbWarpedMotionParams *params_l1 = &picture_control_set_ptr->parent_pcs_ptr
                        ->global_motion[svt_get_ref_frame_type(REF_LIST_1, list_ref_index_l1)];

                if (isCompoundEnabled && allow_bipred
                    && (params_l0->wmtype > TRANSLATION && params_l1->wmtype > TRANSLATION)) {
                    /**************
                    GLOBAL_GLOBALMV
                    ************* */

                    IntMv mv_l1 = gm_get_motion_vector_enc(
                        params_l1,
                        picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv,
                        context_ptr->blk_geom->bsize,
                        mi_col, mi_row,
                        0 /* force_integer_mv */);

                    int16_t to_inject_mv_x_l1 = mv_l1.as_mv.col;
                    int16_t to_inject_mv_y_l1 = mv_l1.as_mv.row;

                    inside_tile = 1;
                    if (umv0tile) {
                        inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l1,
                                                              mi_col, mi_row, context_ptr->blk_geom->bsize)
                                      && is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l1,
                                                                 mi_col, mi_row, context_ptr->blk_geom->bsize);
                    }

                    if (inside_tile) {
                        MvReferenceFrame rf[2];
                        rf[0] = svt_get_ref_frame_type(REF_LIST_0, list_ref_index_l0);
                        rf[1] = svt_get_ref_frame_type(REF_LIST_1, list_ref_index_l1);
                        uint8_t to_inject_ref_type = av1_ref_frame_type(rf);

                        // Warped prediction is only compatible with MD_COMP_AVG and MD_COMP_DIST.
                        for (cur_type = MD_COMP_AVG; cur_type <= MIN(MD_COMP_DIST, tot_comp_types); cur_type++)
                        {
                            candidateArray[canTotalCnt].type = INTER_MODE;
                            candidateArray[canTotalCnt].distortion_ready = 0;
                            candidateArray[canTotalCnt].use_intrabc = 0;

                            candidateArray[canTotalCnt].merge_flag = EB_FALSE;

                            candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;

                            candidateArray[canTotalCnt].inter_mode = GLOBAL_GLOBALMV;
                            candidateArray[canTotalCnt].pred_mode = GLOBAL_GLOBALMV;
                            candidateArray[canTotalCnt].motion_mode = params_l0->wmtype > TRANSLATION ? WARPED_CAUSAL : SIMPLE_TRANSLATION;
                            candidateArray[canTotalCnt].wm_params_l0 = *params_l0;
                            candidateArray[canTotalCnt].wm_params_l1 = *params_l1;
                            candidateArray[canTotalCnt].is_compound = 1;

#if II_COMP_FLAG
                            candidateArray[canTotalCnt].is_interintra_used = 0;
#endif

                            candidateArray[canTotalCnt].is_new_mv = 0;
                            candidateArray[canTotalCnt].is_zero_mv = 0;
                            candidateArray[canTotalCnt].drl_index = 0;

                            // will be needed later by the rate estimation
                            candidateArray[canTotalCnt].ref_mv_index = 0;
                            candidateArray[canTotalCnt].pred_mv_weight = 0;
                            candidateArray[canTotalCnt].ref_frame_type = to_inject_ref_type;
                            candidateArray[canTotalCnt].ref_frame_index_l0 = 0;
                            candidateArray[canTotalCnt].ref_frame_index_l1 = 0;
                            candidateArray[canTotalCnt].transform_type[0] = DCT_DCT;
                            candidateArray[canTotalCnt].transform_type_uv = DCT_DCT;
                            // Set the MV to frame MV

                            candidateArray[canTotalCnt].motion_vector_xl0 = to_inject_mv_x_l0;
                            candidateArray[canTotalCnt].motion_vector_yl0 = to_inject_mv_y_l0;
                            candidateArray[canTotalCnt].motion_vector_xl1 = to_inject_mv_x_l1;
                            candidateArray[canTotalCnt].motion_vector_yl1 = to_inject_mv_y_l1;
                            //GLOB-GLOB
                            determine_compound_mode(
                                picture_control_set_ptr,
                                context_ptr,
                                &candidateArray[canTotalCnt],
                                cur_type);
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
#endif
#if GM_OPT && GLOBAL_WARPED_MOTION || !GLOBAL_WARPED_MOTION
#if GM_OPT && GLOBAL_WARPED_MOTION
    }else{
#endif
        /**************
         GLOBALMV L0
        ************* */
        {
            int16_t to_inject_mv_x = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
            int16_t to_inject_mv_y = (int16_t)(picture_control_set_ptr->parent_pcs_ptr->global_motion[LAST_FRAME].wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
            uint8_t to_inject_ref_type = svt_get_ref_frame_type(REF_LIST_0, 0/*list0_ref_index*/);
            inj_mv = context_ptr->injected_mv_count_l0 == 0 || mrp_is_already_injected_mv_l0(context_ptr, to_inject_mv_x, to_inject_mv_y, to_inject_ref_type) == EB_FALSE;
            if(umv0tile)
                inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x, to_inject_mv_y, mi_col, mi_row, context_ptr->blk_geom->bsize);
            inj_mv = inj_mv && inside_tile;
            if (inj_mv) {


#if II_COMP_FLAG      // GLOBALMV L0
             MvReferenceFrame rf[2];
             rf[0] = to_inject_ref_type;
             rf[1] = -1;

            uint8_t inter_type;
            uint8_t is_ii_allowed = svt_is_interintra_allowed(picture_control_set_ptr->parent_pcs_ptr->enable_inter_intra, bsize, GLOBALMV, rf);
            uint8_t tot_inter_types = is_ii_allowed ? II_COUNT : 1;
            //uint8_t is_obmc_allowed =  obmc_motion_mode_allowed(picture_control_set_ptr, context_ptr->cu_ptr, bsize, rf[0], rf[1], NEWMV) == OBMC_CAUSAL;
            //tot_inter_types = is_obmc_allowed ? tot_inter_types+1 : tot_inter_types;

            for (inter_type = 0; inter_type < tot_inter_types; inter_type++)
            {
#endif
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

#if II_COMP_FLAG
            if (inter_type == 0) {
                candidateArray[canTotalCnt].is_interintra_used = 0;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
            }
            else {
                if (is_ii_allowed) {
                    if (inter_type == 1) {
                        inter_intra_search(
                            picture_control_set_ptr,
                            context_ptr,
                            &candidateArray[canTotalCnt]);
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].use_wedge_interintra = 1;
                        candidateArray[canTotalCnt].ii_wedge_sign = 0;
                    }
                    else if (inter_type == 2) {
                        candidateArray[canTotalCnt].is_interintra_used = 1;
                        candidateArray[canTotalCnt].interintra_mode = candidateArray[canTotalCnt - 1].interintra_mode;
                        candidateArray[canTotalCnt].use_wedge_interintra = 0;
                    }
                }
                //if (is_obmc_allowed && inter_type == tot_inter_types - 1) {
                //    candidateArray[canTotalCnt].is_interintra_used = 0;
                //    candidateArray[canTotalCnt].motion_mode = OBMC_CAUSAL;
                //}
            }
#endif
                INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);

#if II_COMP_FLAG
            }
#endif
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
            inside_tile = 1;
            if(umv0tile) {
                inside_tile = is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize) &&
                              is_inside_tile_boundary(&(xd->tile), to_inject_mv_x_l0, to_inject_mv_y_l1, mi_col, mi_row, context_ptr->blk_geom->bsize);
            }
            if (inside_tile && (context_ptr->injected_mv_count_bipred == 0 || mrp_is_already_injected_mv_bipred(context_ptr, to_inject_mv_x_l0, to_inject_mv_y_l0, to_inject_mv_x_l1, to_inject_mv_y_l1, to_inject_ref_type) == EB_FALSE)) {
                context_ptr->variance_ready = 0;
                for (cur_type = MD_COMP_AVG; cur_type <= tot_comp_types; cur_type++)
                {
                    if (cur_type == MD_COMP_WEDGE && wedge_params_lookup[context_ptr->blk_geom->bsize].bits == 0) continue;
                    // If two predictors are very similar, skip wedge compound mode search
                    if (context_ptr->variance_ready)
                        if (context_ptr->prediction_mse < 8 || (!have_newmv_in_inter_mode(GLOBAL_GLOBALMV) && context_ptr->prediction_mse < 64))
                            continue;
                candidateArray[canTotalCnt].type = INTER_MODE;
                candidateArray[canTotalCnt].distortion_ready = 0;
                candidateArray[canTotalCnt].use_intrabc = 0;

                candidateArray[canTotalCnt].merge_flag = EB_FALSE;

                candidateArray[canTotalCnt].prediction_direction[0] = (EbPredDirection)2;

                candidateArray[canTotalCnt].inter_mode = GLOBAL_GLOBALMV;
                candidateArray[canTotalCnt].pred_mode = GLOBAL_GLOBALMV;
                candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
                candidateArray[canTotalCnt].is_compound = 1;
#if II_COMP_FLAG
                candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
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
                //GLOB-GLOB
                determine_compound_mode(
                    picture_control_set_ptr,
                    context_ptr,
                    &candidateArray[canTotalCnt],
                    cur_type);
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
#if GM_OPT && GLOBAL_WARPED_MOTION
    }
#endif
#endif
    }

    // Warped Motion
    if (frm_hdr->allow_warped_motion &&
        has_overlappable_candidates(context_ptr->cu_ptr) &&
        context_ptr->blk_geom->bwidth >= 8 &&
        context_ptr->blk_geom->bheight >= 8 &&
        context_ptr->warped_motion_injection) {
        inject_warped_motion_candidates(
            picture_control_set_ptr,
            context_ptr,
            context_ptr->cu_ptr,
            &canTotalCnt,
            me_results);
    }
#if ENHANCED_M0_SETTINGS
    if (!coeff_based_nsq_cand_reduction) {
#else
    if (inject_newmv_candidate) {
#endif
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
                    &canTotalCnt);
            }
        }


        if (context_ptr->predictive_me_level)
            inject_predictive_me_candidates(
                context_ptr,
                picture_control_set_ptr,
                isCompoundEnabled,
                allow_bipred,
                &canTotalCnt);
// update the total number of candidates injected
(*candidateTotalCnt) = canTotalCnt;

return;
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
    mbmi.block_mi.mode = pred_mode;
    mbmi.block_mi.uv_mode = pred_mode_uv;

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
            tx_type = intra_mode_to_tx_type(&mbmi.block_mi, PLANE_TYPE_UV);
        }
    }
    ASSERT(tx_type < TX_TYPES);
    if (!av1_ext_tx_used[tx_set_type][tx_type]) return DCT_DCT;
    return tx_type;
}

void  inject_intra_candidates_ois(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    SuperBlock                   *sb_ptr,
    uint32_t                       *candidate_total_cnt){
    uint8_t                     intra_candidate_counter;
    uint8_t                     intra_mode;
    uint32_t                    can_total_cnt = 0;
    ModeDecisionCandidate    *candidate_array = context_ptr->fast_candidate_array;
    EbBool                      disable_cfl_flag = (MAX(context_ptr->blk_geom->bheight, context_ptr->blk_geom->bwidth) > 32) ? EB_TRUE : EB_FALSE;

    OisSbResults    *ois_sb_results_ptr = picture_control_set_ptr->parent_pcs_ptr->ois_sb_results[sb_ptr->index];
    OisCandidate     *ois_blk_ptr = ois_sb_results_ptr->ois_candidate_array[ep_to_pa_block_index[context_ptr->blk_geom->blkidx_mds]];
    uint8_t              total_intra_luma_mode = ois_sb_results_ptr-> total_ois_intra_candidate[ep_to_pa_block_index[context_ptr->blk_geom->blkidx_mds]];
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

    for (intra_candidate_counter = 0; intra_candidate_counter < total_intra_luma_mode; ++intra_candidate_counter) {
        intra_mode = ois_blk_ptr[can_total_cnt].intra_mode;
        assert(intra_mode < INTRA_MODES);
        if (av1_is_directional_mode((PredictionMode)intra_mode)) {
            int32_t angle_delta = ois_blk_ptr[can_total_cnt].angle_delta ;
            candidate_array[can_total_cnt].type = INTRA_MODE;
#if PAL_SUP
            candidate_array[can_total_cnt].palette_info.pmi.palette_size[0] = 0;
            candidate_array[can_total_cnt].palette_info.pmi.palette_size[1] = 0;
#endif
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
                    frm_hdr->reduced_tx_set);
            candidate_array[can_total_cnt].ref_frame_type = INTRA_FRAME;
            candidate_array[can_total_cnt].pred_mode = (PredictionMode)intra_mode;
            candidate_array[can_total_cnt].motion_mode = SIMPLE_TRANSLATION;
            INCRMENT_CAND_TOTAL_COUNT(can_total_cnt);
        }
        else {
            candidate_array[can_total_cnt].type = INTRA_MODE;
#if PAL_SUP
            candidate_array[can_total_cnt].palette_info.pmi.palette_size[0] = 0;
            candidate_array[can_total_cnt].palette_info.pmi.palette_size[1] = 0;
#endif
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
                    frm_hdr->reduced_tx_set);
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

double eb_av1_convert_qindex_to_q(int32_t qindex, AomBitDepth bit_depth);

#if !OBMC_FLAG
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
void eb_av1_setup_pred_block(BlockSize sb_type,
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
#endif
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
        const double q = eb_av1_convert_qindex_to_q(i, bit_depth);
        bit16lut[i] = (int)(0.0418 * q + 2.4107);
        bit4lut[i] = (int)(0.063 * q + 2.742);
    }
}

void eb_av1_init_me_luts(void) {
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
    FrameHeader *frm_hdr = &pcs->parent_pcs_ptr->frm_hdr;
    const Av1Common *const cm = pcs->parent_pcs_ptr->av1_cm;
    MvReferenceFrame ref_frame = INTRA_FRAME;
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
    x->sadperbit16 = sad_per_bit16lut_8[frm_hdr->quantization_params.base_q_idx];
    x->errorperbit = context_ptr->full_lambda >> RD_EPB_SHIFT;
    x->errorperbit += (x->errorperbit == 0);
    //temp buffer for hash me
    for (int xi = 0; xi < 2; xi++)
        for (int yj = 0; yj < 2; yj++)
            x->hash_value_buffer[xi][yj] = (uint32_t*)malloc(AOM_BUFFER_SIZE_FOR_BLOCK_HASH * sizeof(uint32_t));

    IntMv nearestmv, nearmv;
    eb_av1_find_best_ref_mvs_from_stack(0, context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack /*mbmi_ext*/, xd, ref_frame, &nearestmv, &nearmv,
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
    eb_av1_setup_pred_block(bsize, yv12_mb, &cur_buf, mi_row, mi_col);
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

        eb_av1_set_mv_search_range(&x->mv_limits, &dv_ref.as_mv);

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

        const int bestsme = eb_av1_full_pixel_search(
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
        cu_ptr,
        dv_cand,
        &num_dv_cand);

    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    uint32_t dv_i;

    for (dv_i = 0; dv_i < num_dv_cand; dv_i++)
    {
#if PAL_SUP
        candidateArray[*cand_cnt].palette_info.pmi.palette_size[0] = 0;
        candidateArray[*cand_cnt].palette_info.pmi.palette_size[1] = 0;
#endif
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
#if II_COMP_FLAG
        candidateArray[*cand_cnt].is_interintra_used = 0;
#endif
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
        candidateArray[*cand_cnt].filter_intra_mode = FILTER_INTRA_MODES;
        INCRMENT_CAND_TOTAL_COUNT( (*cand_cnt) );
    }
}
// Indices are sign, integer, and fractional part of the gradient value
static const uint8_t gradient_to_angle_bin[2][7][16] = {
  {
      { 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 },
      { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
      { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
      { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
      { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 },
      { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 },
  },
  {
      { 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4 },
      { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3 },
      { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 },
      { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 },
      { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 },
      { 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2 },
      { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 },
  },
};

/* clang-format off */
static const uint8_t mode_to_angle_bin[INTRA_MODES] = {
  0, 2, 6, 0, 4, 3, 5, 7, 1, 0,
  0,
};
void av1_get_gradient_hist_c(const uint8_t *src, int src_stride, int rows,
    int cols, uint64_t *hist) {
    src += src_stride;
    for (int r = 1; r < rows; ++r) {
        for (int c = 1; c < cols; ++c) {
            int dx = src[c] - src[c - 1];
            int dy = src[c] - src[c - src_stride];
            int index;
            const int temp = dx * dx + dy * dy;
            if (dy == 0) {
                index = 2;
            }
            else {
                const int sn = (dx > 0) ^ (dy > 0);
                dx = abs(dx);
                dy = abs(dy);
                const int remd = (dx % dy) * 16 / dy;
                const int quot = dx / dy;
                index = gradient_to_angle_bin[sn][AOMMIN(quot, 6)][AOMMIN(remd, 15)];
            }
            hist[index] += temp;
        }
        src += src_stride;
    }
}
static void angle_estimation(
    const uint8_t *src,
    int src_stride,
    int rows,
    int cols,
    //BLOCK_SIZE bsize,
    uint8_t *directional_mode_skip_mask)
{
    // Check if angle_delta is used
    //if (!av1_use_angle_delta(bsize)) return;

    uint64_t hist[DIRECTIONAL_MODES] = { 0 };
    //if (is_hbd)
    //    get_highbd_gradient_hist(src, src_stride, rows, cols, hist);
    //else
    av1_get_gradient_hist(src, src_stride, rows, cols, hist);

    int i;
    uint64_t hist_sum = 0;
    for (i = 0; i < DIRECTIONAL_MODES; ++i) hist_sum += hist[i];
    for (i = 0; i < INTRA_MODES; ++i) {
        if (av1_is_directional_mode(i)) {
            const uint8_t angle_bin = mode_to_angle_bin[i];
            uint64_t score = 2 * hist[angle_bin];
            int weight = 2;
            if (angle_bin > 0) {
                score += hist[angle_bin - 1];
                ++weight;
            }
            if (angle_bin < DIRECTIONAL_MODES - 1) {
                score += hist[angle_bin + 1];
                ++weight;
            }
            const int thresh = 10;
            if (score * thresh < hist_sum * weight) directional_mode_skip_mask[i] = 1;
        }
    }
}
// END of Function Declarations
void  inject_intra_candidates(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    const SequenceControlSet     *sequence_control_set_ptr,
    SuperBlock                   *sb_ptr,
#if ENHANCED_M0_SETTINGS
    EbBool                        dc_cand_only_flag,
#endif
    uint32_t                     *candidateTotalCnt){
    (void)sequence_control_set_ptr;
    (void)sb_ptr;
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    uint8_t                     intra_mode_start = DC_PRED;
#if MULTI_PASS_PD
#if ENHANCED_M0_SETTINGS
    uint8_t                     intra_mode_end = dc_cand_only_flag ? DC_PRED : PAETH_PRED;
#else
    uint8_t                     intra_mode_end = context_ptr->dc_cand_only_flag ? DC_PRED : PAETH_PRED;
#endif
#else
    uint8_t                     intra_mode_end   =  PAETH_PRED;
#endif

    uint8_t                     openLoopIntraCandidate;
    uint32_t                    canTotalCnt = 0;
    uint8_t                     angleDeltaCounter = 0;
    EbBool                      use_angle_delta = av1_use_angle_delta(context_ptr->blk_geom->bsize);
    uint8_t                     angleDeltaCandidateCount = use_angle_delta ? 7 : 1;
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool                      disable_cfl_flag = (MAX(context_ptr->blk_geom->bheight, context_ptr->blk_geom->bwidth) > 32) ? EB_TRUE : EB_FALSE;

    uint8_t                     disable_z2_prediction;
    uint8_t                     disable_angle_refinement;
    uint8_t                     disable_angle_prediction;
    uint8_t directional_mode_skip_mask[INTRA_MODES] = { 0 };

    if (context_ptr->edge_based_skip_angle_intra && use_angle_delta)
    {
        EbPictureBufferDesc   *src_pic = picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
        uint8_t               *src_buf = src_pic->buffer_y + (context_ptr->cu_origin_x + src_pic->origin_x) + (context_ptr->cu_origin_y + src_pic->origin_y) * src_pic->stride_y;
        const int rows = block_size_high[context_ptr->blk_geom->bsize];
        const int cols = block_size_wide[context_ptr->blk_geom->bsize];
        angle_estimation(src_buf, src_pic->stride_y, rows, cols, /*context_ptr->blk_geom->bsize,*/directional_mode_skip_mask);
    }
    uint8_t     angle_delta_shift = 1;

#if MULTI_PASS_PD
    if (context_ptr->disable_angle_z2_intra_flag) {
        disable_angle_prediction = 1;
        angleDeltaCandidateCount = 1;
        angle_delta_shift = 1;
        disable_z2_prediction = 1;
    }
    else
#endif
    if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode == 4) {
        if (picture_control_set_ptr->slice_type == I_SLICE) {
            intra_mode_end =  PAETH_PRED;
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
            if (!disable_angle_prediction &&
                directional_mode_skip_mask[(PredictionMode)openLoopIntraCandidate] == 0) {
                for (angleDeltaCounter = 0; angleDeltaCounter < angleDeltaCandidateCount; ++angleDeltaCounter) {
                    int32_t angle_delta = CLIP( angle_delta_shift * (angleDeltaCandidateCount == 1 ? 0 : angleDeltaCounter - (angleDeltaCandidateCount >> 1)), -3 , 3);
                    int32_t  p_angle = mode_to_angle_map[(PredictionMode)openLoopIntraCandidate] + angle_delta * ANGLE_STEP;
                    if (!disable_z2_prediction || (p_angle <= 90 || p_angle >= 180)) {
                        candidateArray[canTotalCnt].type = INTRA_MODE;
                        candidateArray[canTotalCnt].palette_info.pmi.palette_size[0] = 0;
                        candidateArray[canTotalCnt].palette_info.pmi.palette_size[1] = 0;
                        candidateArray[canTotalCnt].intra_luma_mode = openLoopIntraCandidate;
                        candidateArray[canTotalCnt].distortion_ready = 0;
                        candidateArray[canTotalCnt].use_intrabc = 0;
                        candidateArray[canTotalCnt].filter_intra_mode = FILTER_INTRA_MODES;
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
                                frm_hdr->reduced_tx_set);
                        candidateArray[canTotalCnt].ref_frame_type = INTRA_FRAME;
                        candidateArray[canTotalCnt].pred_mode = (PredictionMode)openLoopIntraCandidate;
                        candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
#if II_COMP_FLAG
                        candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
                        INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
                    }
            }
        }
        }
        else {
            candidateArray[canTotalCnt].type = INTRA_MODE;
#if PAL_SUP
            candidateArray[canTotalCnt].palette_info.pmi.palette_size[0] = 0;
            candidateArray[canTotalCnt].palette_info.pmi.palette_size[1] = 0;
#endif
            candidateArray[canTotalCnt].intra_luma_mode = openLoopIntraCandidate;
            candidateArray[canTotalCnt].distortion_ready = 0;
            candidateArray[canTotalCnt].use_intrabc = 0;
            candidateArray[canTotalCnt].filter_intra_mode = FILTER_INTRA_MODES;
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
                    frm_hdr->reduced_tx_set);
            candidateArray[canTotalCnt].ref_frame_type = INTRA_FRAME;
            candidateArray[canTotalCnt].pred_mode = (PredictionMode)openLoopIntraCandidate;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;
#if II_COMP_FLAG
            candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
            INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
        }
    }

    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;

    return;
}
// END of Function Declarations
void  inject_filter_intra_candidates(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    uint32_t                     *candidateTotalCnt){

    FilterIntraMode             intra_mode_start = FILTER_DC_PRED;
    FilterIntraMode             intra_mode_end   = FILTER_INTRA_MODES;
    FilterIntraMode             filter_intra_mode;
    uint32_t                    canTotalCnt = *candidateTotalCnt;
    ModeDecisionCandidate      *candidateArray = context_ptr->fast_candidate_array;

    EbBool                      disable_cfl_flag = (MAX(context_ptr->blk_geom->bheight, context_ptr->blk_geom->bwidth) > 32) ? EB_TRUE : EB_FALSE;

    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

    for (filter_intra_mode = intra_mode_start; filter_intra_mode < intra_mode_end ; ++filter_intra_mode) {

            candidateArray[canTotalCnt].type = INTRA_MODE;
            candidateArray[canTotalCnt].intra_luma_mode = DC_PRED;
            candidateArray[canTotalCnt].distortion_ready = 0;
            candidateArray[canTotalCnt].use_intrabc = 0;
            candidateArray[canTotalCnt].filter_intra_mode = filter_intra_mode;
            candidateArray[canTotalCnt].is_directional_mode_flag = 0;

#if PAL_SUP
            candidateArray[canTotalCnt].palette_info.pmi.palette_size[0] = 0;
            candidateArray[canTotalCnt].palette_info.pmi.palette_size[1] = 0;
#endif

            candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y] = 0;

            // Search the best independent intra chroma mode
            if (context_ptr->chroma_level == CHROMA_MODE_0) {
                candidateArray[canTotalCnt].intra_chroma_mode  = disable_cfl_flag ? context_ptr->best_uv_mode[fimode_to_intramode[filter_intra_mode]][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]] : UV_CFL_PRED;

                candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_UV] = disable_cfl_flag ? context_ptr->best_uv_angle[fimode_to_intramode[filter_intra_mode]][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]] : 0;
                candidateArray[canTotalCnt].is_directional_chroma_mode_flag = disable_cfl_flag ? (uint8_t)av1_is_directional_mode((PredictionMode)(context_ptr->best_uv_mode[fimode_to_intramode[filter_intra_mode]][MAX_ANGLE_DELTA + candidateArray[canTotalCnt].angle_delta[PLANE_TYPE_Y]])) : 0;

            }
            else {
                // Hsan/Omar: why the restriction below ? (i.e. disable_ang_uv)
                const int32_t disable_ang_uv = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) && context_ptr->blk_geom->has_uv ? 1 : 0;
                candidateArray[canTotalCnt].intra_chroma_mode = disable_cfl_flag ? intra_luma_to_chroma[fimode_to_intramode[filter_intra_mode]] :
                    (context_ptr->chroma_level == CHROMA_MODE_1) ?
                        UV_CFL_PRED :
                        UV_DC_PRED;

                candidateArray[canTotalCnt].intra_chroma_mode =  disable_ang_uv && av1_is_directional_mode(candidateArray[canTotalCnt].intra_chroma_mode) ?
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
                    frm_hdr->reduced_tx_set);

            candidateArray[canTotalCnt].ref_frame_type = INTRA_FRAME;
            candidateArray[canTotalCnt].pred_mode = DC_PRED;
            candidateArray[canTotalCnt].motion_mode = SIMPLE_TRANSLATION;

#if II_COMP_FLAG
            candidateArray[canTotalCnt].is_interintra_used = 0;
#endif
            INCRMENT_CAND_TOTAL_COUNT(canTotalCnt);
    }

    // update the total number of candidates injected
    (*candidateTotalCnt) = canTotalCnt;

    return;
}

#if PAL_SUP
int svt_av1_allow_palette(int allow_palette,
    BlockSize sb_type) {
    assert(sb_type < BlockSizeS_ALL);
    return allow_palette && block_size_wide[sb_type] <= 64 &&
        block_size_high[sb_type] <= 64 && sb_type >= BLOCK_8X8;
}
void  search_palette_luma(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    PaletteInfo                 *palette_cand,
    uint32_t                     *tot_palette_cands);

void  inject_palette_candidates(
    PictureControlSet            *picture_control_set_ptr,
    ModeDecisionContext          *context_ptr,
    uint32_t                       *candidate_total_cnt) {



    uint32_t                  can_total_cnt = *candidate_total_cnt;
    ModeDecisionCandidate    *candidateArray = context_ptr->fast_candidate_array;
    EbBool                    disable_cfl_flag = (MAX(context_ptr->blk_geom->bheight, context_ptr->blk_geom->bwidth) > 32) ? EB_TRUE : EB_FALSE;
    uint32_t cand_i;
    uint32_t tot_palette_cands = 0;
    PaletteInfo    *palette_cand_array = context_ptr->palette_cand_array;

    search_palette_luma(
        picture_control_set_ptr,
        context_ptr,
        palette_cand_array,
        &tot_palette_cands);

    for (cand_i = 0; cand_i < tot_palette_cands; ++cand_i) {

        palette_cand_array[cand_i].pmi.palette_size[1] = 0;
        memcpy(candidateArray[can_total_cnt].palette_info.color_idx_map, palette_cand_array[cand_i].color_idx_map, 64 * 64);
        memcpy(&candidateArray[can_total_cnt].palette_info.pmi, &palette_cand_array[cand_i].pmi, sizeof(PaletteModeInfo));
        assert(palette_cand_array[cand_i].pmi.palette_size[0] < 9);
        //to re check these fields
        candidateArray[can_total_cnt].type = INTRA_MODE;
        candidateArray[can_total_cnt].intra_luma_mode = DC_PRED;
        candidateArray[can_total_cnt].distortion_ready = 0;
        candidateArray[can_total_cnt].use_intrabc = 0;

        candidateArray[can_total_cnt].filter_intra_mode = FILTER_INTRA_MODES;

        candidateArray[can_total_cnt].is_directional_mode_flag = 0;

        candidateArray[can_total_cnt].angle_delta[PLANE_TYPE_Y] = 0;

        // Search the best independent intra chroma mode
        if (context_ptr->chroma_level == CHROMA_MODE_0) {
            candidateArray[can_total_cnt].intra_chroma_mode = disable_cfl_flag ?
                context_ptr->best_uv_mode[DC_PRED][MAX_ANGLE_DELTA + candidateArray[can_total_cnt].angle_delta[PLANE_TYPE_Y]] :
                UV_CFL_PRED;

            candidateArray[can_total_cnt].angle_delta[PLANE_TYPE_UV] = disable_cfl_flag ?
                context_ptr->best_uv_angle[candidateArray[can_total_cnt].intra_luma_mode][MAX_ANGLE_DELTA + candidateArray[can_total_cnt].angle_delta[PLANE_TYPE_Y]] : 0;
            candidateArray[can_total_cnt].is_directional_chroma_mode_flag = disable_cfl_flag ?
                (uint8_t)av1_is_directional_mode((PredictionMode)(context_ptr->best_uv_mode[candidateArray[can_total_cnt].intra_luma_mode][MAX_ANGLE_DELTA + candidateArray[can_total_cnt].angle_delta[PLANE_TYPE_Y]])) : 0;

        }
        else {
            // Hsan/Omar: why the restriction below ? (i.e. disable_ang_uv)
            const int32_t disable_ang_uv = (context_ptr->blk_geom->bwidth == 4 || context_ptr->blk_geom->bheight == 4) && context_ptr->blk_geom->has_uv ? 1 : 0;
            candidateArray[can_total_cnt].intra_chroma_mode = disable_cfl_flag ?
                intra_luma_to_chroma[DC_PRED] :
                (context_ptr->chroma_level == CHROMA_MODE_1) ?
                UV_CFL_PRED :
                UV_DC_PRED;

            candidateArray[can_total_cnt].intra_chroma_mode = disable_ang_uv && av1_is_directional_mode(candidateArray[can_total_cnt].intra_chroma_mode) ?
                UV_DC_PRED : candidateArray[can_total_cnt].intra_chroma_mode;

            candidateArray[can_total_cnt].is_directional_chroma_mode_flag = (uint8_t)av1_is_directional_mode((PredictionMode)candidateArray[can_total_cnt].intra_chroma_mode);
            candidateArray[can_total_cnt].angle_delta[PLANE_TYPE_UV] = 0;

        }

        candidateArray[can_total_cnt].cfl_alpha_signs = 0;
        candidateArray[can_total_cnt].cfl_alpha_idx = 0;
        candidateArray[can_total_cnt].transform_type[0] = DCT_DCT;

        if (candidateArray[can_total_cnt].intra_chroma_mode == UV_CFL_PRED)
            candidateArray[can_total_cnt].transform_type_uv = DCT_DCT;
        else
            candidateArray[can_total_cnt].transform_type_uv =

            av1_get_tx_type(
                context_ptr->blk_geom->bsize,
                0,
                (PredictionMode)candidateArray[can_total_cnt].intra_luma_mode,
                (UvPredictionMode)candidateArray[can_total_cnt].intra_chroma_mode,
                PLANE_TYPE_UV,
                0,
                0,
                0,
                context_ptr->blk_geom->txsize_uv[0][0],
                picture_control_set_ptr->parent_pcs_ptr->frm_hdr.reduced_tx_set);

        candidateArray[can_total_cnt].ref_frame_type = INTRA_FRAME;
        candidateArray[can_total_cnt].pred_mode = (PredictionMode)DC_PRED;
        candidateArray[can_total_cnt].motion_mode = SIMPLE_TRANSLATION;
        INCRMENT_CAND_TOTAL_COUNT(can_total_cnt);
    }

    // update the total number of candidates injected
    (*candidate_total_cnt) = can_total_cnt;

    return;
}
#endif
EbErrorType generate_md_stage_0_cand(
    SuperBlock          *sb_ptr,
    ModeDecisionContext *context_ptr,
    uint32_t            *candidate_total_count_ptr,
    PictureControlSet   *picture_control_set_ptr)
{

    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    const SequenceControlSet *sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    const EB_SLICE slice_type = picture_control_set_ptr->slice_type;
    uint32_t canTotalCnt = 0;
    // Reset duplicates variables
    context_ptr->injected_mv_count_l0 = 0;
    context_ptr->injected_mv_count_l1 = 0;
    context_ptr->injected_mv_count_bipred = 0;
    uint8_t sq_index = LOG2F(context_ptr->blk_geom->sq_size) - 2;
#if ENHANCED_M0_SETTINGS
    EbBool coeff_based_nsq_cand_reduction = EB_FALSE;
#else
    uint8_t inject_intra_candidate = 1;
    uint8_t inject_inter_candidate = 1;
#endif
    if (slice_type != I_SLICE) {
#if MULTI_PASS_PD // Shut coef-based inter skip if 1st pass
        if (context_ptr->coeff_based_nsq_cand_reduction) {
#else
        if (picture_control_set_ptr->parent_pcs_ptr->nsq_search_level >= NSQ_SEARCH_LEVEL1 &&
            picture_control_set_ptr->parent_pcs_ptr->nsq_search_level < NSQ_SEARCH_FULL) {
#endif
            if (context_ptr->md_local_cu_unit[context_ptr->blk_geom->sqi_mds].avail_blk_flag)
#if ENHANCED_M0_SETTINGS
                coeff_based_nsq_cand_reduction = context_ptr->blk_geom->shape == PART_N || context_ptr->parent_sq_has_coeff[sq_index] != 0 ? EB_FALSE : EB_TRUE;
#else
                inject_intra_candidate = context_ptr->blk_geom->shape == PART_N ? 1 :
                    context_ptr->parent_sq_has_coeff[sq_index] != 0 ? inject_intra_candidate : 0;
#endif
        }
}
    //----------------------
    // Intra
    if (context_ptr->blk_geom->sq_size < 128) {
#if MULTI_PASS_PD
#if ENHANCED_M0_SETTINGS
        if (!context_ptr->dc_cand_only_flag && !coeff_based_nsq_cand_reduction && picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode >= 5 && context_ptr->blk_geom->sq_size > 4 && context_ptr->blk_geom->shape == PART_N)
#else
        if (!context_ptr->dc_cand_only_flag && picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode >= 5 && context_ptr->blk_geom->sq_size > 4 && context_ptr->blk_geom->shape == PART_N)
#endif
#else
        if (picture_control_set_ptr->parent_pcs_ptr->intra_pred_mode >= 5 && context_ptr->blk_geom->sq_size > 4 && context_ptr->blk_geom->shape == PART_N)
#endif
            inject_intra_candidates_ois(
                picture_control_set_ptr,
                context_ptr,
                sb_ptr,
                &canTotalCnt);
        else
#if !ENHANCED_M0_SETTINGS
            if (inject_intra_candidate)
#endif
                inject_intra_candidates(
                    picture_control_set_ptr,
                    context_ptr,
                    sequence_control_set_ptr,
                    sb_ptr,
#if ENHANCED_M0_SETTINGS
                    context_ptr->dc_cand_only_flag || coeff_based_nsq_cand_reduction,
#endif
                &canTotalCnt);
    }
#if MULTI_PASS_PD
#if ENHANCED_M0_SETTINGS
    if (!coeff_based_nsq_cand_reduction)
#endif
       if (context_ptr->md_filter_intra_mode > 0 && av1_filter_intra_allowed_bsize(sequence_control_set_ptr->seq_header.enable_filter_intra, context_ptr->blk_geom->bsize))
#else
       if (picture_control_set_ptr->pic_filter_intra_mode > 0 && av1_filter_intra_allowed_bsize(sequence_control_set_ptr->seq_header.enable_filter_intra, context_ptr->blk_geom->bsize))
#endif
            inject_filter_intra_candidates(
                picture_control_set_ptr,
                context_ptr,
                &canTotalCnt);


    if (frm_hdr->allow_intrabc)
        inject_intra_bc_candidates(
            picture_control_set_ptr,
            context_ptr,
            sequence_control_set_ptr,
            context_ptr->cu_ptr,
            &canTotalCnt
        );

#if PAL_SUP
    //can be removed later if need be
    for (uint16_t i = 0; i < canTotalCnt; i++) {
        assert(context_ptr->fast_candidate_array[i].palette_info.pmi.palette_size[0] == 0);
        assert(context_ptr->fast_candidate_array[i].palette_info.pmi.palette_size[1] == 0);
    }
    if (svt_av1_allow_palette(picture_control_set_ptr->parent_pcs_ptr->palette_mode, context_ptr->blk_geom->bsize)) {
        inject_palette_candidates(
            picture_control_set_ptr,
            context_ptr,
            &canTotalCnt);
    }
    for (uint16_t i = 0; i < canTotalCnt; i++) {
        assert(context_ptr->fast_candidate_array[i].palette_info.pmi.palette_size[0] < 9);
        assert(context_ptr->fast_candidate_array[i].palette_info.pmi.palette_size[1] == 0);
    }
#endif

    if (slice_type != I_SLICE) {
#if !ENHANCED_M0_SETTINGS
        if (inject_inter_candidate)
#endif
            inject_inter_candidates(
                picture_control_set_ptr,
                context_ptr,
                sequence_control_set_ptr,
                sb_ptr,
#if ENHANCED_M0_SETTINGS
                coeff_based_nsq_cand_reduction,
#endif
                &canTotalCnt);
    }

    *candidate_total_count_ptr = canTotalCnt;
    CAND_CLASS  cand_class_it;
    memset(context_ptr->md_stage_0_count, 0, CAND_CLASS_TOTAL * sizeof(uint32_t));

    uint32_t cand_i;
    for (cand_i = 0; cand_i < canTotalCnt; cand_i++)
    {
        ModeDecisionCandidate * cand_ptr = &context_ptr->fast_candidate_array[cand_i];

        if (cand_ptr->type == INTRA_MODE) {
            // Intra prediction
                if (cand_ptr->filter_intra_mode == FILTER_INTRA_MODES) {
                  if (cand_ptr->palette_info.pmi.palette_size[0] == 0) {
                    cand_ptr->cand_class = CAND_CLASS_0;
                    context_ptr->md_stage_0_count[CAND_CLASS_0]++;
                  }
                  else {
                     cand_ptr->cand_class = CAND_CLASS_7;
                     context_ptr->md_stage_0_count[CAND_CLASS_7]++;
                  }
                }
                else {
                    cand_ptr->cand_class = CAND_CLASS_6;
                    context_ptr->md_stage_0_count[CAND_CLASS_6]++;
                }

        } else if (cand_ptr->inter_mode == GLOBALMV || cand_ptr->inter_mode == GLOBAL_GLOBALMV) {
            cand_ptr->cand_class = CAND_CLASS_8;
            context_ptr->md_stage_0_count[CAND_CLASS_8]++;
        }
#if OBMC_FLAG
        else {
            // Inter pred
            if (cand_ptr->motion_mode == OBMC_CAUSAL) {
                // OBMC
                cand_ptr->cand_class = CAND_CLASS_5;
                context_ptr->md_stage_0_count[CAND_CLASS_5]++;
            }
            else if (cand_ptr->is_compound == 0 ||
                    (cand_ptr->is_compound == 1 && cand_ptr->interinter_comp.type == COMPOUND_AVERAGE)) {
#else
            else if ((cand_ptr->type == INTER_MODE && cand_ptr->is_compound == 0) ||
                (cand_ptr->type == INTER_MODE && cand_ptr->is_compound == 1 && cand_ptr->interinter_comp.type == COMPOUND_AVERAGE)) {
#endif
#if II_COMP_FLAG
                if (cand_ptr->is_interintra_used && cand_ptr->is_compound == 0) {
                    // InterIntra
                    cand_ptr->cand_class = CAND_CLASS_4;
                    context_ptr->md_stage_0_count[CAND_CLASS_4]++;

                }
                else
#endif
                if (context_ptr->combine_class12) {
                    cand_ptr->cand_class = CAND_CLASS_1;
                    context_ptr->md_stage_0_count[CAND_CLASS_1]++;
                }
                else {
                    if (cand_ptr->is_new_mv) {
                        // ME pred
                        cand_ptr->cand_class = CAND_CLASS_1;
                        context_ptr->md_stage_0_count[CAND_CLASS_1]++;
                    }
                    else {
                        // MV pred
                        cand_ptr->cand_class = CAND_CLASS_2;
                        context_ptr->md_stage_0_count[CAND_CLASS_2]++;
                    }
                }
            }
            else {
#if !OBMC_FLAG
                if (context_ptr->combine_class12) {
                    cand_ptr->cand_class = CAND_CLASS_2;
                    context_ptr->md_stage_0_count[CAND_CLASS_2]++;
                }
                else {
#endif
                // InterInter
                cand_ptr->cand_class = CAND_CLASS_3;
                context_ptr->md_stage_0_count[CAND_CLASS_3]++;
#if !OBMC_FLAG
                }
#endif
            }
#if OBMC_FLAG
        }
#endif
    }
    uint32_t fast_accum = 0;
    for (cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++) {
        fast_accum += context_ptr->md_stage_0_count[cand_class_it];
    }
    assert(fast_accum == canTotalCnt);

    return EB_ErrorNone;
}

/***************************************
* Full Mode Decision
***************************************/
uint32_t product_full_mode_decision(
    struct ModeDecisionContext   *context_ptr,
    CodingUnit                   *cu_ptr,
    ModeDecisionCandidateBuffer **buffer_ptr_array,
    uint32_t                      candidate_total_count,
    uint32_t                     *best_candidate_index_array,
    uint8_t                       prune_ref_frame_for_rec_partitions,
    uint32_t                     *best_intra_mode)
{
    uint32_t                  candidateIndex;
    uint64_t                  lowestCost = 0xFFFFFFFFFFFFFFFFull;
    uint64_t                  lowestIntraCost = 0xFFFFFFFFFFFFFFFFull;
    uint32_t                  lowestCostIndex = 0;
    if (prune_ref_frame_for_rec_partitions) {
        if (context_ptr->blk_geom->shape == PART_N) {
            for (uint32_t i = 0; i < candidate_total_count; ++i) {
                candidateIndex = best_candidate_index_array[i];
                ModeDecisionCandidate *candidate_ptr = buffer_ptr_array[candidateIndex]->candidate_ptr;
                EbBool is_inter = (candidate_ptr->pred_mode >= NEARESTMV) ? EB_TRUE : EB_FALSE;
                EbBool is_simple_translation = (candidate_ptr->motion_mode != WARPED_CAUSAL) ? EB_TRUE : EB_FALSE;
                if (is_inter && is_simple_translation) {
                    uint8_t ref_frame_type = candidate_ptr->ref_frame_type;
                    assert(ref_frame_type < MAX_REF_TYPE_CAND);
                    context_ptr->ref_best_cost_sq_table[ref_frame_type] = *(buffer_ptr_array[candidateIndex]->full_cost_ptr);
                }

            }
            //Sort ref_frame by sq cost.
            uint32_t i, j, index;
            for (i = 0; i < MAX_REF_TYPE_CAND; ++i) {
                context_ptr->ref_best_ref_sq_table[i] = i;
            }
            for (i = 0; i < MAX_REF_TYPE_CAND - 1; ++i) {
                for (j = i + 1; j < MAX_REF_TYPE_CAND; ++j) {
                    if (context_ptr->ref_best_cost_sq_table[context_ptr->ref_best_ref_sq_table[j]] < context_ptr->ref_best_cost_sq_table[context_ptr->ref_best_ref_sq_table[i]]) {
                        index = context_ptr->ref_best_ref_sq_table[i];
                        context_ptr->ref_best_ref_sq_table[i] = context_ptr->ref_best_ref_sq_table[j];
                        context_ptr->ref_best_ref_sq_table[j] = index;
                    }
                }
            }
        }
    }


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
    if (cu_ptr->prediction_mode_flag == INTER_MODE && candidate_ptr->is_compound)
    {
        cu_ptr->interinter_comp.type = candidate_ptr->interinter_comp.type;
        cu_ptr->interinter_comp.mask_type = candidate_ptr->interinter_comp.mask_type;
        cu_ptr->interinter_comp.wedge_index = candidate_ptr->interinter_comp.wedge_index;
        cu_ptr->interinter_comp.wedge_sign = candidate_ptr->interinter_comp.wedge_sign;
        cu_ptr->compound_idx = candidate_ptr->compound_idx;
        cu_ptr->comp_group_idx = candidate_ptr->comp_group_idx;
        if (cu_ptr->interinter_comp.type == COMPOUND_AVERAGE){
            if (cu_ptr->comp_group_idx != 0 || cu_ptr->compound_idx != 1)
                printf("Error: Compound combination not allowed\n");
        }
    }
#if II_COMP_FLAG

    cu_ptr->is_interintra_used          = candidate_ptr->is_interintra_used;
    cu_ptr->interintra_mode             = candidate_ptr->interintra_mode;
    cu_ptr->use_wedge_interintra        = candidate_ptr->use_wedge_interintra;
    cu_ptr->interintra_wedge_index      = candidate_ptr->interintra_wedge_index;
    cu_ptr->ii_wedge_sign               = candidate_ptr->ii_wedge_sign;

#endif

    // Set the PU level variables
    cu_ptr->interp_filters = candidate_ptr->interp_filters;
    {
        pu_ptr = cu_ptr->prediction_unit_array;
        // Intra Prediction
        pu_ptr->intra_luma_mode = 0x1F;
        if (cu_ptr->prediction_mode_flag == INTRA_MODE)
        {
            cu_ptr->filter_intra_mode= candidate_ptr->filter_intra_mode;
            pu_ptr->intra_luma_mode = candidate_ptr->intra_luma_mode;

            pu_ptr->is_directional_mode_flag = candidate_ptr->is_directional_mode_flag;
            pu_ptr->angle_delta[PLANE_TYPE_Y] = candidate_ptr->angle_delta[PLANE_TYPE_Y];

            pu_ptr->cfl_alpha_idx = candidate_ptr->cfl_alpha_idx;
            pu_ptr->cfl_alpha_signs = candidate_ptr->cfl_alpha_signs;

            pu_ptr->intra_chroma_mode = candidate_ptr->intra_chroma_mode;
            pu_ptr->is_directional_chroma_mode_flag = candidate_ptr->is_directional_chroma_mode_flag;
            pu_ptr->angle_delta[PLANE_TYPE_UV] = candidate_ptr->angle_delta[PLANE_TYPE_UV];
        }

#if PAL_SUP
        if (cu_ptr->prediction_mode_flag == INTRA_MODE)
        {
            memcpy(&cu_ptr->palette_info.pmi, &candidate_ptr->palette_info.pmi, sizeof(PaletteModeInfo));
            if(svt_av1_allow_palette(context_ptr->sb_ptr->picture_control_set_ptr->parent_pcs_ptr->palette_mode, context_ptr->blk_geom->bsize))
               memcpy(cu_ptr->palette_info.color_idx_map, candidate_ptr->palette_info.color_idx_map, MAX_PALETTE_SQUARE);
        }
        else {
            cu_ptr->palette_info.pmi.palette_size[0] = cu_ptr->palette_info.pmi.palette_size[1] = 0;
        }
#endif
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
        pu_ptr->compound_idx = candidate_ptr->compound_idx;
        pu_ptr->interinter_comp = candidate_ptr->interinter_comp;
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
        // The MV prediction indicies are recalcated by the EncDec.
        pu_ptr->mvd[REF_LIST_0].pred_idx = 0;
        pu_ptr->mvd[REF_LIST_1].pred_idx = 0;

        pu_ptr->overlappable_neighbors[0] = context_ptr->cu_ptr->prediction_unit_array[0].overlappable_neighbors[0];
        pu_ptr->overlappable_neighbors[1] = context_ptr->cu_ptr->prediction_unit_array[0].overlappable_neighbors[1];
        pu_ptr->motion_mode = candidate_ptr->motion_mode;
        pu_ptr->num_proj_ref = candidate_ptr->num_proj_ref;
        if (pu_ptr->motion_mode == WARPED_CAUSAL) {
            EB_MEMCPY(&pu_ptr->wm_params_l0, &candidate_ptr->wm_params_l0, sizeof(EbWarpedMotionParams));
            EB_MEMCPY(&pu_ptr->wm_params_l1, &candidate_ptr->wm_params_l1, sizeof(EbWarpedMotionParams));
        }
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
// clang-format on
