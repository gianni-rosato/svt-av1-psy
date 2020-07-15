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

#ifndef EbEntropyCoding_h
#define EbEntropyCoding_h

#include "EbDefinitions.h"
#include "EbEntropyCodingObject.h"
#include "EbCodingUnit.h"
#include "EbPredictionUnit.h"
#include "EbPictureBufferDesc.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbCabacContextModel.h"
#include "EbModeDecision.h"
#include "EbEncIntraPrediction.h"
#include "EbBitstreamUnit.h"
#include "EbPacketizationProcess.h"
#include "EbModeDecisionProcess.h"
#include "EbInterPrediction.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ModeDecisionCandidateBuffer;
struct ModeDecisionCandidate;

/**************************************
     * Extern Function Declarations
     **************************************/
struct EntropyCodingContext;
extern EbErrorType write_sb(struct EntropyCodingContext *context_ptr, SuperBlock *tb_ptr,
                            PictureControlSet *pcs_ptr, uint16_t tile_idx,
                            EntropyCoder *entropy_coder_ptr,
                            EbPictureBufferDesc *coeff_ptr);

extern int get_wedge_params_bits(BlockSize sb_type);

extern EbErrorType encode_slice_finish(EntropyCoder *entropy_coder_ptr);

extern EbErrorType reset_entropy_coder(EncodeContext *encode_context_ptr,
                                       EntropyCoder *entropy_coder_ptr, uint32_t qp,
                                       EB_SLICE slice_type);
#if MD_FRAME_CONTEXT_MEM_OPT
EbErrorType av1_txb_estimate_coeff_bits(
    struct ModeDecisionContext *md_context, uint8_t allow_update_cdf, FRAME_CONTEXT *ec_ctx,
    PictureControlSet *pcs_ptr, struct ModeDecisionCandidateBuffer *candidate_buffer_ptr,
    uint32_t txb_origin_index, uint32_t txb_chroma_origin_index,
    EbPictureBufferDesc *coeff_buffer_sb, uint32_t y_eob, uint32_t cb_eob, uint32_t cr_eob,
    uint64_t *y_txb_coeff_bits, uint64_t *cb_txb_coeff_bits, uint64_t *cr_txb_coeff_bits,
    TxSize txsize, TxSize txsize_uv, TxType tx_type, TxType tx_type_uv,
    COMPONENT_TYPE component_type);
#else
EbErrorType av1_txb_estimate_coeff_bits(
    struct ModeDecisionContext *md_context, uint8_t allow_update_cdf, FRAME_CONTEXT *ec_ctx,
    PictureControlSet *pcs_ptr, struct ModeDecisionCandidateBuffer *candidate_buffer_ptr,
    uint32_t txb_origin_index, uint32_t txb_chroma_origin_index, EntropyCoder *entropy_coder_ptr,
    EbPictureBufferDesc *coeff_buffer_sb, uint32_t y_eob, uint32_t cb_eob, uint32_t cr_eob,
    uint64_t *y_txb_coeff_bits, uint64_t *cb_txb_coeff_bits, uint64_t *cr_txb_coeff_bits,
    TxSize txsize, TxSize txsize_uv, TxType tx_type, TxType tx_type_uv,
    COMPONENT_TYPE component_type);
#endif

//**********************************************************************************************************//
//onyxc_int.h
static INLINE int32_t frame_is_intra_only(const PictureParentControlSet *const pcs_ptr) {
    return pcs_ptr->frm_hdr.frame_type == KEY_FRAME ||
           pcs_ptr->frm_hdr.frame_type == INTRA_ONLY_FRAME;
}

static INLINE int32_t frame_is_sframe(const PictureParentControlSet *pcs_ptr) {
    return pcs_ptr->frm_hdr.frame_type == S_FRAME;
}

// Returns 1 if this frame might allow mvs from some reference frame.

static INLINE int32_t frame_might_allow_ref_frame_mvs(const PictureParentControlSet *pcs_ptr,
                                                      SequenceControlSet *           scs_ptr) {
    return !pcs_ptr->frm_hdr.error_resilient_mode &&
           scs_ptr->seq_header.order_hint_info.enable_ref_frame_mvs &&
           scs_ptr->seq_header.order_hint_info.enable_order_hint && !frame_is_intra_only(pcs_ptr);
}

// Returns 1 if this frame might use warped_motion
static INLINE int32_t frame_might_allow_warped_motion(const PictureParentControlSet *pcs_ptr,
                                                      SequenceControlSet *           scs_ptr) {
    return !pcs_ptr->frm_hdr.error_resilient_mode && !frame_is_intra_only(pcs_ptr) &&
           scs_ptr->static_config.enable_warped_motion;
}

static INLINE uint8_t major_minor_to_seq_level_idx(BitstreamLevel bl) {
    assert(bl.major >= LEVEL_MAJOR_MIN && bl.major <= LEVEL_MAJOR_MAX);
    return ((bl.major - LEVEL_MAJOR_MIN) << LEVEL_MINOR_BITS) + bl.minor;
}

//**********************************************************************************************************//
//encoder.h
static INLINE int32_t get_ref_frame_map_idx(const PictureParentControlSet *pcs_ptr,
                                            MvReferenceFrame               ref_frame) {
    return pcs_ptr->av1_ref_signal
            .ref_dpb_index[ref_frame - LAST_FRAME]; //LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
}

//*******************************************************************************************//
// bitwriter_buffer.h
struct AomWriteBitBuffer {
    uint8_t *bit_buffer;
    uint32_t bit_offset;
};

int32_t  eb_aom_wb_is_byte_aligned(const struct AomWriteBitBuffer *wb);
uint32_t eb_aom_wb_bytes_written(const struct AomWriteBitBuffer *wb);

void eb_aom_wb_write_bit(struct AomWriteBitBuffer *wb, int32_t bit);
void eb_aom_wb_write_literal(struct AomWriteBitBuffer *wb, int32_t data, int32_t bits);

void eb_aom_wb_write_inv_signed_literal(struct AomWriteBitBuffer *wb, int32_t data, int32_t bits);
//*******************************************************************************************//
// Bitstream.h
struct AomWriteBitBuffer;

void write_sequence_header(SequenceControlSet *      scs_ptr /*Av1Comp *cpi*/,
                           struct AomWriteBitBuffer *wb);

uint32_t write_obu_header(ObuType ObuType, int32_t obuExtension, uint8_t *const dst);

int32_t write_uleb_obu_size(uint32_t obu_header_size, uint32_t obu_payload_size, uint8_t *dest);

//*******************************************************************************************//
// blockd.h

void get_txb_ctx(PictureControlSet *pcs_ptr, const int32_t plane,
                 NeighborArrayUnit *dc_sign_level_coeff_neighbor_array, uint32_t blk_origin_x,
                 uint32_t blk_origin_y, const BlockSize plane_bsize, const TxSize tx_size,
                 int16_t *const txb_skip_ctx, int16_t *const dc_sign_ctx);

extern int32_t eb_av1_get_reference_mode_context(uint32_t blk_origin_x, uint32_t blk_origin_y,
                                                 NeighborArrayUnit *mode_type_neighbor_array,
                                                 NeighborArrayUnit *inter_pred_dir_neighbor_array);

extern int32_t eb_av1_get_comp_reference_type_context(
    uint32_t blk_origin_x, uint32_t blk_origin_y, NeighborArrayUnit *mode_type_neighbor_array,
    NeighborArrayUnit *inter_pred_dir_neighbor_array);

extern void av1_collect_neighbors_ref_counts_new(MacroBlockD *const xd);
// Obtain contexts to signal a reference frame be either BWDREF/ALTREF2, or
// ALTREF.
//extern int32_t get_pred_context_brfarf2_or_arf(const MacroBlockD *xd);
// Obtain contexts to signal a reference frame be either BWDREF or ALTREF2.
//extern int32_t get_pred_context_brf_or_arf2(const MacroBlockD *xd);
// == Context functions for comp ref ==
//
// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode be either
// GOLDEN/LAST3, or LAST/LAST2.
extern int32_t eb_av1_get_pred_context_comp_ref_p(const MacroBlockD *xd);

// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode be LAST,
// conditioning on that it is known either LAST/LAST2.
extern int32_t eb_av1_get_pred_context_comp_ref_p1(const MacroBlockD *xd);

// Returns a context number for the given MB prediction signal
// Signal the first reference frame for a compound mode be GOLDEN,
// conditioning on that it is known either GOLDEN or LAST3.
extern int32_t eb_av1_get_pred_context_comp_ref_p2(const MacroBlockD *xd);

// Signal the 2nd reference frame for a compound mode be either
// ALTREF, or ALTREF2/BWDREF.
extern int32_t eb_av1_get_pred_context_comp_bwdref_p(const MacroBlockD *xd);

// Signal the 2nd reference frame for a compound mode be either
// ALTREF2 or BWDREF.
extern int32_t eb_av1_get_pred_context_comp_bwdref_p1(const MacroBlockD *xd);
// == Context functions for single ref ==
//
// For the bit to signal whether the single reference is a forward reference
// frame or a backward reference frame.
extern int32_t eb_av1_get_pred_context_single_ref_p1(const MacroBlockD *xd);

// For the bit to signal whether the single reference is ALTREF_FRAME or
// non-ALTREF backward reference frame, knowing that it shall be either of
// these 2 choices.
extern int32_t eb_av1_get_pred_context_single_ref_p2(const MacroBlockD *xd);

// For the bit to signal whether the single reference is LAST3/GOLDEN or
// LAST2/LAST, knowing that it shall be either of these 2 choices.
extern int32_t eb_av1_get_pred_context_single_ref_p3(const MacroBlockD *xd);

// For the bit to signal whether the single reference is LAST2_FRAME or
// LAST_FRAME, knowing that it shall be either of these 2 choices.
extern int32_t eb_av1_get_pred_context_single_ref_p4(const MacroBlockD *xd);

// For the bit to signal whether the single reference is GOLDEN_FRAME or
// LAST3_FRAME, knowing that it shall be either of these 2 choices.
extern int32_t eb_av1_get_pred_context_single_ref_p5(const MacroBlockD *xd);

// For the bit to signal whether the single reference is ALTREF2_FRAME or
// BWDREF_FRAME, knowing that it shall be either of these 2 choices.
extern int32_t eb_av1_get_pred_context_single_ref_p6(const MacroBlockD *xd);

extern EbErrorType write_frame_header_av1(Bitstream *bitstream_ptr, SequenceControlSet *scs_ptr,
                                          PictureControlSet *pcs_ptr, uint8_t show_existing);
extern EbErrorType encode_td_av1(uint8_t *bitstream_ptr);
extern EbErrorType encode_sps_av1(Bitstream *bitstream_ptr, SequenceControlSet *scs_ptr);

//*******************************************************************************************//

MotionMode motion_mode_allowed(const PictureControlSet *pcs_ptr, const BlkStruct *blk_ptr,
                               const BlockSize bsize, MvReferenceFrame rf0, MvReferenceFrame rf1,
                               PredictionMode mode);

int is_masked_compound_type(COMPOUND_TYPE type);


int32_t eb_aom_count_primitive_subexpfin(uint16_t n, uint16_t k, uint16_t v);
int32_t eb_aom_count_primitive_refsubexpfin(uint16_t n, uint16_t k, uint16_t ref, uint16_t v);

#ifdef __cplusplus
}
#endif
#endif //EbEntropyCoding_h
