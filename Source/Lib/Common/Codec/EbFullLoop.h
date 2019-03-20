/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbFullLoop_h
#define EbFullLoop_h

#include "EbModeDecisionProcess.h"
#ifdef __cplusplus
extern "C" {
#endif

    void FullLoop_R(
        LargestCodingUnit_t            *sb_ptr,
        ModeDecisionCandidateBuffer_t  *candidateBuffer,
        ModeDecisionContext_t          *context_ptr,
        EbPictureBufferDesc_t          *input_picture_ptr,
        PictureControlSet_t            *picture_control_set_ptr,
        uint32_t                          component_mask,
        uint32_t                          cbQp,
        uint32_t                          crQp,
        uint32_t                          *cb_count_non_zero_coeffs,
        uint32_t                          *cr_count_non_zero_coeffs);

    void CuFullDistortionFastTuMode_R(
        LargestCodingUnit_t            *sb_ptr,
        ModeDecisionCandidateBuffer_t  *candidateBuffer,
        ModeDecisionContext_t            *context_ptr,
        ModeDecisionCandidate_t           *candidate_ptr,
        PictureControlSet_t            *picture_control_set_ptr,
        uint64_t                          cbFullDistortion[DIST_CALC_TOTAL],
        uint64_t                          crFullDistortion[DIST_CALC_TOTAL],
        uint32_t                          count_non_zero_coeffs[3][MAX_NUM_OF_TU_PER_CU],

        COMPONENT_TYPE                  component_type,

        uint64_t                         *cb_coeff_bits,
        uint64_t                         *cr_coeff_bits,
        EbAsm                            asm_type);

    void ProductFullLoop(
        ModeDecisionCandidateBuffer_t  *candidateBuffer,
        ModeDecisionContext_t          *context_ptr,
        PictureControlSet_t            *picture_control_set_ptr,
        uint32_t                          qp,
        uint32_t                           *y_count_non_zero_coeffs,
        uint64_t                         *y_coeff_bits,
        uint64_t                         *y_full_distortion);


    void ProductFullLoopTxSearch(
        ModeDecisionCandidateBuffer_t  *candidateBuffer,
        ModeDecisionContext_t          *context_ptr,
        PictureControlSet_t            *picture_control_set_ptr);
    extern uint32_t d2_inter_depth_block_decision(
        ModeDecisionContext_t          *context_ptr,
        uint32_t                        blk_mds,
        LargestCodingUnit_t            *tbPtr,
        uint32_t                          lcuAddr,
        uint32_t                          tbOriginX,
        uint32_t                          tbOriginY,
        uint64_t                          full_lambda,
        MdRateEstimationContext_t      *md_rate_estimation_ptr,
        PictureControlSet_t            *picture_control_set_ptr);


    void  d1_non_square_block_decision(
        ModeDecisionContext_t               *context_ptr
    );



#ifdef __cplusplus
}
#endif
#endif // EbFullLoop_h