/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbRateDistortionCost_h
#define EbRateDistortionCost_h

/***************************************
 * Includes
 ***************************************/
#include "EbModeDecisionProcess.h"
#include "EbEncIntraPrediction.h"
#include "EbEncInterPrediction.h"
#include "EbLambdaRateTables.h"
#include "EbTransforms.h"
#include "EbEncDecProcess.h"
#include "EbEntropyCoding.h"

#ifdef __cplusplus
extern "C" {
#endif
extern uint64_t svt_av1_cost_coeffs_txb(struct ModeDecisionContext *ctx, uint8_t allow_update_cdf,
                                        FRAME_CONTEXT *ec_ctx, struct ModeDecisionCandidateBuffer *cand_bf,
                                        const TranLow *const qcoeff, uint16_t eob, PlaneType plane_type,
                                        TxSize transform_size, TxType transform_type, int16_t txb_skip_ctx,
                                        int16_t dc_sign_ctx, Bool reduced_transform_set_flag);
void            svt_aom_coding_loop_context_generation(PictureControlSet *pcs, ModeDecisionContext *ctx);
#define RDDIV_BITS 7

#define RDCOST(RM, R, D)                                                         \
    (ROUND_POWER_OF_TWO(((int64_t)(R)) * ((int64_t)(RM)), AV1_PROB_COST_SHIFT) + \
     ((int64_t)(D) * ((int64_t)1 << RDDIV_BITS)))

extern uint64_t svt_aom_partition_rate_cost(PictureParentControlSet *pcs, ModeDecisionContext *ctx,
                                            uint32_t blk_mds_idx, PartitionType p, uint64_t lambda,
                                            bool use_accurate_part_ctx, MdRateEstimationContext *md_rate_est_ctx);
uint64_t        svt_aom_get_intra_uv_fast_rate(PictureControlSet *pcs, struct ModeDecisionContext *ctx,
                                               ModeDecisionCandidateBuffer *cand_bf, bool use_accurate_cfl);
uint64_t        svt_aom_intra_fast_cost(PictureControlSet *pcs, struct ModeDecisionContext *ctx,
                                        ModeDecisionCandidateBuffer *cand_bf, uint64_t lambda, uint64_t luma_distortion,
                                        uint64_t chroma_distortion);
uint64_t        svt_aom_inter_fast_cost(PictureControlSet *pcs, struct ModeDecisionContext *ctx,
                                        ModeDecisionCandidateBuffer *cand_bf, uint64_t lambda, uint64_t luma_distortion,
                                        uint64_t chroma_distortion);
EbErrorType     svt_aom_full_cost_light_pd0(ModeDecisionContext *ctx, struct ModeDecisionCandidateBuffer *cand_bf,
                                            uint64_t *y_distortion, uint64_t lambda, uint64_t *y_coeff_bits);
void svt_aom_full_cost(PictureControlSet *pcs, ModeDecisionContext *ctx, struct ModeDecisionCandidateBuffer *cand_bf,
                       uint64_t lambda, uint64_t y_distortion[DIST_TOTAL][DIST_CALC_TOTAL],
                       uint64_t cb_distortion[DIST_TOTAL][DIST_CALC_TOTAL],
                       uint64_t cr_distortion[DIST_TOTAL][DIST_CALC_TOTAL], uint64_t *y_coeff_bits,
                       uint64_t *cb_coeff_bits, uint64_t *cr_coeff_bits);

uint64_t svt_aom_tx_size_bits(MdRateEstimationContext *md_rate_est_ctx, MacroBlockD *xd, const MbModeInfo *mbmi,
                              TxSize tx_size, TxMode tx_mode, BlockSize bsize, uint8_t skip, FRAME_CONTEXT *ec_ctx,
                              uint8_t allow_update_cdf);

uint64_t svt_aom_get_tx_size_bits(ModeDecisionCandidateBuffer *candidateBuffer, ModeDecisionContext *ctx,
                                  PictureControlSet *pcs, uint8_t tx_depth, Bool block_has_coeff);

MvJointType svt_av1_get_mv_joint(const MV *mv);
int32_t svt_av1_mv_bit_cost(const MV *mv, const MV *ref, const int32_t *mvjcost, int32_t *mvcost[2], int32_t weight);
int32_t svt_av1_mv_bit_cost_light(const MV *mv, const MV *ref);
static INLINE uint8_t av1_drl_ctx(const CandidateMv *ref_mv_stack, int32_t ref_idx) {
    return ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL   ? ref_mv_stack[ref_idx + 1].weight >= REF_CAT_LEVEL ? 0 : 1
        : ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL ? 2
                                                           : 0;
}
#ifdef __cplusplus
}
#endif
#endif //EbRateDistortionCost_h
