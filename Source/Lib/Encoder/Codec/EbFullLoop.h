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

#ifndef EbFullLoop_h
#define EbFullLoop_h

#include "EbModeDecisionProcess.h"
#include "EbCommonUtils.h"
#include "EbInvTransforms.h"
#include "EbTransforms.h"
#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

#if OPT_CHROMA_PATH
void full_loop_r(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr,
                 ModeDecisionCandidateBuffer *candidate_buffer,
                 EbPictureBufferDesc *input_picture_ptr,
                 COMPONENT_TYPE component_type, uint32_t chroma_qindex,
                 uint32_t count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU],
                 uint64_t cb_full_distortion[DIST_CALC_TOTAL],
                 uint64_t cr_full_distortion[DIST_CALC_TOTAL],
                 uint64_t *cb_coeff_bits,
                 uint64_t *cr_coeff_bits,
                 EbBool is_full_loop);
#else
void full_loop_r(SuperBlock *sb_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
                 ModeDecisionContext *context_ptr, EbPictureBufferDesc *input_picture_ptr,
                 PictureControlSet *pcs_ptr, uint32_t component_mask, uint32_t cb_qindex,
                 uint32_t cr_qindex, uint32_t *cb_count_non_zero_coeffs,
                 uint32_t *cr_count_non_zero_coeffs);
void cu_full_distortion_fast_txb_mode_r(
    SuperBlock *sb_ptr, ModeDecisionCandidateBuffer *candidate_buffer,
    ModeDecisionContext *context_ptr, ModeDecisionCandidate *candidate_ptr,
    PictureControlSet *pcs_ptr, EbPictureBufferDesc *input_picture_ptr,
    uint64_t cb_full_distortion[DIST_CALC_TOTAL], uint64_t cr_full_distortion[DIST_CALC_TOTAL],
    uint32_t count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU], COMPONENT_TYPE component_type,
    uint64_t *cb_coeff_bits, uint64_t *cr_coeff_bits, EbBool is_full_loop);
#endif
void inv_transform_recon_wrapper(uint8_t *pred_buffer, uint32_t pred_offset, uint32_t pred_stride,
                                 uint8_t *rec_buffer, uint32_t rec_offset, uint32_t rec_stride,
                                 int32_t *rec_coeff_buffer, uint32_t coeff_offset, EbBool hbd,
                                 TxSize txsize, TxType transform_type, PlaneType component_type,
                                 uint32_t eob);

extern uint32_t d2_inter_depth_block_decision(SequenceControlSet* scs_ptr,
                                              PictureControlSet* pcs_ptr,
                                              ModeDecisionContext* context_ptr,
                                              uint32_t blk_mds,
                                              uint32_t sb_addr);
// compute the cost of curr depth, and the depth above
extern void compute_depth_costs_md_skip(ModeDecisionContext *    context_ptr,
                                        SequenceControlSet *     scs_ptr,
                                        PictureParentControlSet *pcs_ptr, uint32_t above_depth_mds,
                                        uint32_t step, uint64_t *above_depth_cost,
                                        uint64_t *curr_depth_cost);
#if LIGHT_PD0_2
void compute_depth_costs_md_skip_light_pd0(ModeDecisionContext *context_ptr,
                                           uint32_t above_depth_mds, uint32_t step,
                                           uint64_t *above_depth_cost, uint64_t *curr_depth_cost);
#endif
uint64_t    d1_non_square_block_decision(ModeDecisionContext *context_ptr, uint32_t d1_block_itr);

static const int av1_get_tx_scale_tab[TX_SIZES_ALL] = {
    0, 0, 0, 1, 2, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 1, 1};

static const TxSize get_txsize_entropy_ctx_tab[TX_SIZES_ALL] = {
    0, 1, 2, 3, 4, 1, 1, 2, 2, 3, 3, 4, 4, 1, 1, 2, 2, 3, 3};

static const int get_txb_bwl_tab[TX_SIZES_ALL] = {
    2, 3, 4, 5, 5, 2, 3, 3, 4, 4, 5, 5, 5, 2, 4, 3, 5, 4, 5};

static const int get_txb_wide_tab[TX_SIZES_ALL] = {
    4, 8, 16, 32, 32, 4, 8, 8, 16, 16, 32, 32, 32, 4, 16, 8, 32, 16, 32};

static const int get_txb_high_tab[TX_SIZES_ALL] = {
    4, 8, 16, 32, 32, 8, 4, 16, 8, 32, 16, 32, 32, 16, 4, 32, 8, 32, 16};
#ifdef __cplusplus
}
#endif
#endif // EbFullLoop_h
