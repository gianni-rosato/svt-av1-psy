/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbCodingLoop_h
#define EbCodingLoop_h

#include "EbDefinitions.h"
#include "EbCodingUnit.h"
#include "EbSequenceControlSet.h"
#include "EbModeDecisionProcess.h"
#include "EbEncDecProcess.h"
#ifdef __cplusplus
extern "C" {
#endif
    /*******************************************
     * ModeDecisionLCU
     *   performs CL (LCU)
     *******************************************/

    typedef EbErrorType(*EB_MODE_DECISION)(
        SequenceControlSet                *sequence_control_set_ptr,
        PictureControlSet                 *picture_control_set_ptr,
        const MdcLcuData * const           mdcResultTbPtr,
        LargestCodingUnit                 *sb_ptr,
        uint32_t                             sb_origin_x,
        uint32_t                             sb_origin_y,
        uint32_t                             lcuAddr,
        ModeDecisionContext               *context_ptr);

    extern EbErrorType AV1ProductModeDecisionLcu(
        SequenceControlSet                *sequence_control_set_ptr,
        PictureControlSet                 *picture_control_set_ptr,
        const MdcLcuData * const           mdcResultTbPtr,
        LargestCodingUnit                 *sb_ptr,
        uint32_t                             sb_origin_x,
        uint32_t                             sb_origin_y,
        uint32_t                             lcuAddr,
        ModeDecisionContext               *context_ptr);

    extern EbErrorType in_loop_motion_estimation_sblock(
        PictureControlSet                 *picture_control_set_ptr,  // input parameter, Picture Control Set Ptr
        uint32_t                             sb_origin_x,            // input parameter, SB Origin X
        uint32_t                             sb_origin_y,            // input parameter, SB Origin X
        int16_t                              x_mv_l0,
        int16_t                              y_mv_l0,
        int16_t                              x_mv_l1,
        int16_t                              y_mv_l1,
        SsMeContext                        *context_ptr);          // input parameter, ME Context Ptr, used to store decimated/interpolated LCU/SR

    extern EbErrorType mode_decision_sb(
        SequenceControlSet                *sequence_control_set_ptr,
        PictureControlSet                 *picture_control_set_ptr,
        const MdcLcuData * const           mdcResultTbPtr,
        LargestCodingUnit                 *sb_ptr,
        uint16_t                             sb_origin_x,
        uint16_t                             sb_origin_y,
        uint32_t                             lcuAddr,
        SsMeContext                       *ss_mecontext,
        ModeDecisionContext               *context_ptr);

    extern EbErrorType ModeDecisionRefinementLcu(
        SequenceControlSet                *sequence_control_set_ptr,
        PictureControlSet                 *picture_control_set_ptr,
        LargestCodingUnit                 *sb_ptr,
        uint32_t                               sb_origin_x,
        uint32_t                               sb_origin_y,
        ModeDecisionContext               *context_ptr);
#if !MEMORY_FOOTPRINT_OPT
    extern EbErrorType QpmDeriveWeightsMinAndMax(
        PictureControlSet                    *picture_control_set_ptr,
        EncDecContext                        *context_ptr);
#endif
    uint8_t get_skip_tx_search_flag(
#if BYPASS_USELESS_TX_SEARCH
        const BlockGeom         *blk_geom,
#else
        int32_t                  sq_size,
#endif
        uint64_t                 ref_fast_cost,
        uint64_t                 cu_cost,
        uint64_t                 weight);

    extern void av1_encode_pass(
        SequenceControlSet    *sequence_control_set_ptr,
        PictureControlSet     *picture_control_set_ptr,
        LargestCodingUnit     *sb_ptr,
        uint32_t                   tbAddr,
        uint32_t                   sb_origin_x,
        uint32_t                   sb_origin_y,
#if !MEMORY_FOOTPRINT_OPT
        uint32_t                   sb_qp,
#endif
        EncDecContext         *context_ptr);

#if NO_ENCDEC

    void no_enc_dec_pass(
        SequenceControlSet    *sequence_control_set_ptr,
        PictureControlSet     *picture_control_set_ptr,
        LargestCodingUnit     *sb_ptr,
        uint32_t                   tbAddr,
        uint32_t                   sb_origin_x,
        uint32_t                   sb_origin_y,
        uint32_t                   sb_qp,
        EncDecContext         *context_ptr);
#endif

#ifdef __cplusplus
}
#endif
#endif // EbCodingLoop_h
