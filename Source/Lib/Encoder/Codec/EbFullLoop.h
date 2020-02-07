/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbFullLoop_h
#define EbFullLoop_h

#include "EbModeDecisionProcess.h"
#include "EbCommonUtils.h"
#include "EbInvTransforms.h"
#include "EbTransforms.h"

#ifdef __cplusplus
extern "C" {
#endif

void full_loop_r(SuperBlock *sb_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
                 ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                 PictureControlSet *pcs_ptr, uint32_t component_mask, uint32_t cb_qp,
                 uint32_t cr_qp, uint32_t *cb_count_non_zero_coeffs,
                 uint32_t *cr_count_non_zero_coeffs);

void cu_full_distortion_fast_txb_mode_r(
    SuperBlock *sb_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
    ModeDecisionContext *context_ptr, ModeDecisionCandidate *candidate_ptr,
    PictureControlSet *pcs_ptr, EbPictureBufferDesc *input_picture_ptr,
    uint64_t cb_full_distortion[DIST_CALC_TOTAL], uint64_t cr_full_distortion[DIST_CALC_TOTAL],
    uint32_t count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU], COMPONENT_TYPE component_type,
    uint64_t *cb_coeff_bits, uint64_t *cr_coeff_bits, EbBool is_full_loop);

void product_full_loop(ModeDecisionCandidateBuffer *candidate_buffer,
                       ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                       EbPictureBufferDesc *input_picture_ptr, uint32_t qp,
                       uint32_t *y_count_non_zero_coeffs, uint64_t *y_coeff_bits,
                       uint64_t *y_full_distortion);

void inv_transform_recon_wrapper(uint8_t *pred_buffer, uint32_t pred_offset, uint32_t pred_stride,
                                 uint8_t *rec_buffer, uint32_t rec_offset, uint32_t rec_stride,
                                 int32_t *rec_coeff_buffer, uint32_t coeff_offset, EbBool hbd,
                                 TxSize txsize, TxType transform_type, PlaneType component_type,
                                 uint32_t eob);

extern uint32_t d2_inter_depth_block_decision(ModeDecisionContext *context_ptr, uint32_t blk_mds,
                                              SuperBlock *tb_ptr, uint32_t sb_addr,
                                              uint32_t tb_origin_x, uint32_t tb_origin_y,
                                              uint64_t                 full_lambda,
                                              MdRateEstimationContext *md_rate_estimation_ptr,
                                              PictureControlSet *      pcs_ptr);

// compute the cost of curr depth, and the depth above
extern void compute_depth_costs_md_skip(ModeDecisionContext *context_ptr,
                                        SequenceControlSet *scs_ptr,
                                        PictureParentControlSet *pcs_ptr,
                                        uint32_t above_depth_mds,
                                        uint32_t step, uint64_t *above_depth_cost,
                                        uint64_t *curr_depth_cost);
uint64_t    d1_non_square_block_decision(ModeDecisionContext *context_ptr, uint32_t d1_block_itr);

#ifdef __cplusplus
}
#endif
#endif // EbFullLoop_h
