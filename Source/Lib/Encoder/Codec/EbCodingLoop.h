// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbCodingLoop_h
#define EbCodingLoop_h

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
    extern EbErrorType mode_decision_sb(
        SequenceControlSet                *sequence_control_set_ptr,
        PictureControlSet                 *picture_control_set_ptr,
        const MdcLcuData * const           mdcResultTbPtr,
        SuperBlock                        *sb_ptr,
        uint16_t                             sb_origin_x,
        uint16_t                             sb_origin_y,
        uint32_t                             lcuAddr,
        ModeDecisionContext               *context_ptr);

    uint8_t get_skip_tx_search_flag(
        int32_t                  sq_size,
        uint64_t                 ref_fast_cost,
        uint64_t                 cu_cost,
        uint64_t                 weight);

    extern void av1_encode_pass(
        SequenceControlSet    *sequence_control_set_ptr,
        PictureControlSet     *picture_control_set_ptr,
        SuperBlock            *sb_ptr,
        uint32_t                   tbAddr,
        uint32_t                   sb_origin_x,
        uint32_t                   sb_origin_y,
        EncDecContext         *context_ptr);

#if NO_ENCDEC
    void no_enc_dec_pass(
        SequenceControlSet    *sequence_control_set_ptr,
        PictureControlSet     *picture_control_set_ptr,
        SuperBlock            *sb_ptr,
        uint32_t                   tbAddr,
        uint32_t                   sb_origin_x,
        uint32_t                   sb_origin_y,
        uint32_t                   sb_qp,
        EncDecContext         *context_ptr);
#endif

    void store16bit_input_src(
        EbPictureBufferDesc     *input_sample16bit_buffer,
        PictureControlSet       *picture_control_set_ptr,
        uint32_t                 lcuX,
        uint32_t                 lcuY,
        uint32_t                 lcuW,
        uint32_t                 lcuH);

    void residual_kernel(
        uint8_t   *input,
        uint32_t   input_offset,
        uint32_t   input_stride,
        uint8_t   *pred,
        uint32_t   pred_offset,
        uint32_t   pred_stride,
        int16_t   *residual,
        uint32_t   residual_offset,
        uint32_t   residual_stride,
        EbBool     hbd,
        uint32_t   area_width,
        uint32_t   area_height);

#ifdef __cplusplus
}
#endif
#endif // EbCodingLoop_h
// clang-format on
