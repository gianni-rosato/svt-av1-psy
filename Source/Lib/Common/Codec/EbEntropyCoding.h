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
#include "EbEntropyCodingUtil.h"
#include "EbEntropyCodingProcess.h"

#include "EbCodingUnit.h"
#include "EbPredictionUnit.h"
#include "EbPictureBufferDesc.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "EbCabacContextModel.h"
#include "EbModeDecision.h"
#include "EbIntraPrediction.h"
#include "EbBitstreamUnit.h"
#include "EbPacketizationProcess.h"
#include "EbModeDecisionProcess.h"
#ifdef __cplusplus
extern "C" {
#endif

#define MAX_TILE_WIDTH (4096)        // Max Tile width in pixels
#define MAX_TILE_AREA (4096 * 2304)  // Maximum tile area in pixels

#define CHECK_BACKWARD_REFS(ref_frame) \
  (((ref_frame) >= BWDREF_FRAME) && ((ref_frame) <= ALTREF_FRAME))
#define IS_BACKWARD_REF_FRAME(ref_frame) CHECK_BACKWARD_REFS(ref_frame)

    /**************************************
     * Extern Function Declarations
     **************************************/
    struct EntropyCodingContext;
    extern EbErrorType write_sb(
        struct EntropyCodingContext   *context_ptr,
        LargestCodingUnit     *tb_ptr,
        PictureControlSet     *picture_control_set_ptr,
        EntropyCoder          *entropy_coder_ptr,
        EbPictureBufferDesc   *coeff_ptr);

    extern EbErrorType encode_slice_finish(
        EntropyCoder        *entropy_coder_ptr);

    extern EbErrorType reset_bitstream(
        EbPtr                 bitstream_ptr);

    extern EbErrorType reset_entropy_coder(
        EncodeContext       *encode_context_ptr,
        EntropyCoder        *entropy_coder_ptr,
        uint32_t                 qp,
        EB_SLICE               slice_type);

    extern EbErrorType av1_tu_estimate_coeff_bits(
        struct ModeDecisionContext         *md_context,
        uint8_t                             allow_update_cdf,
        FRAME_CONTEXT                      *ec_ctx,
        PictureControlSet                  *picture_control_set_ptr,
        struct ModeDecisionCandidateBuffer *candidate_buffer_ptr,
        uint32_t                            tu_origin_index,
        uint32_t                            tu_chroma_origin_index,
        EntropyCoder                       *entropy_coder_ptr,
        EbPictureBufferDesc                *coeff_buffer_sb,
        uint32_t                            y_eob,
        uint32_t                            cb_eob,
        uint32_t                            cr_eob,
        uint64_t                            *y_tu_coeff_bits,
        uint64_t                            *cb_tu_coeff_bits,
        uint64_t                            *cr_tu_coeff_bits,
        TxSize                               txsize,
        TxSize                               txsize_uv,
        TxType                               tx_type,
        TxType                               tx_type_uv,
        COMPONENT_TYPE                       component_type);

    extern EbErrorType copy_rbsp_bitstream_to_payload(
        Bitstream *bitstream_ptr,
        EbByte      output_buffer,
        uint32_t      *output_buffer_index,
        uint32_t      *output_buffer_size,
        EncodeContext         *encode_context_ptr);

    //**********************************************************************************************************//
    //onyxc_int.h
    static INLINE int32_t frame_is_intra_only(const PictureParentControlSet *const pcs_ptr) {
        return pcs_ptr->frm_hdr.frame_type == KEY_FRAME || pcs_ptr->frm_hdr.frame_type == INTRA_ONLY_FRAME;
    }

    static INLINE int32_t frame_is_sframe(const PictureParentControlSet *pcs_ptr) {
        return pcs_ptr->frm_hdr.frame_type == S_FRAME;
    }

    // Returns 1 if this frame might allow mvs from some reference frame.

    static INLINE int32_t frame_might_allow_ref_frame_mvs(const PictureParentControlSet *pcs_ptr,
        SequenceControlSet    *scs_ptr) {
        return !pcs_ptr->frm_hdr.error_resilient_mode &&
            scs_ptr->seq_header.order_hint_info.enable_ref_frame_mvs &&
            scs_ptr->seq_header.order_hint_info.enable_order_hint && !frame_is_intra_only(pcs_ptr);
    }

    // Returns 1 if this frame might use warped_motion
    static INLINE int32_t frame_might_allow_warped_motion(const PictureParentControlSet *pcs_ptr,
        SequenceControlSet    *scs_ptr) {
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
        MvReferenceFrame ref_frame) {
        return pcs_ptr->av1_ref_signal.ref_dpb_index[ref_frame - LAST_FRAME];//LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
    }

    static INLINE int is_intrabc_block(const BlockModeInfo *block_mi) {
        return block_mi->use_intrabc;
    }

    static INLINE int bsize_to_max_depth(BlockSize bsize) {
        TxSize tx_size = max_txsize_rect_lookup[bsize];
        int depth = 0;
        while (depth < MAX_TX_DEPTH && tx_size != TX_4X4) {
            depth++;
            tx_size = sub_tx_size_map[tx_size];
        }
        return depth;
    }

    static AomCdfProb cdf_element_prob(const AomCdfProb *const cdf, size_t element) {
        assert(cdf != NULL);
        return (element > 0 ? cdf[element - 1] : CDF_PROB_TOP) - cdf[element];
    }

    static INLINE int bsize_to_tx_size_cat(BlockSize bsize) {
        TxSize tx_size = max_txsize_rect_lookup[bsize];
        assert(tx_size != TX_4X4);
        int depth = 0;
        while (tx_size != TX_4X4) {
            depth++;
            tx_size = sub_tx_size_map[tx_size];
            assert(depth < 10);
        }
        assert(depth <= MAX_TX_CATS);
        return depth - 1;
    }

    static INLINE void partition_gather_horz_alike(AomCdfProb *out,
        const AomCdfProb *const in, BlockSize bsize)
    {
        out[0] = CDF_PROB_TOP;
        out[0] -= cdf_element_prob(in, PARTITION_HORZ);
        out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
        out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
        out[0] -= cdf_element_prob(in, PARTITION_HORZ_B);
        out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
        if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_HORZ_4);
        out[0] = AOM_ICDF(out[0]);
        out[1] = AOM_ICDF(CDF_PROB_TOP);
        out[2] = 0;
    }

    static INLINE void partition_gather_vert_alike(AomCdfProb *out,
        const AomCdfProb *const in, BlockSize bsize)
    {
        out[0] = CDF_PROB_TOP;
        out[0] -= cdf_element_prob(in, PARTITION_VERT);
        out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
        out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
        out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
        out[0] -= cdf_element_prob(in, PARTITION_VERT_B);
        if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_VERT_4);
        out[0] = AOM_ICDF(out[0]);
        out[1] = AOM_ICDF(CDF_PROB_TOP);
        out[2] = 0;
    }

    static INLINE int get_palette_bsize_ctx(BlockSize bsize) {
        return num_pels_log2_lookup[bsize] - num_pels_log2_lookup[BLOCK_8X8];
    }

    //*******************************************************************************************//
    // bitwriter_buffer.h
    struct AomWriteBitBuffer
    {
        uint8_t *bit_buffer;
        uint32_t bit_offset;
    };

    int32_t eb_aom_wb_is_byte_aligned(const struct AomWriteBitBuffer *wb);
    uint32_t eb_aom_wb_bytes_written(const struct AomWriteBitBuffer *wb);

    void eb_aom_wb_write_bit(struct AomWriteBitBuffer *wb, int32_t bit);

    void eb_aom_wb_overwrite_bit(struct AomWriteBitBuffer *wb, int32_t bit);

    void eb_aom_wb_write_literal(struct AomWriteBitBuffer *wb, int32_t data, int32_t bits);

    void eb_aom_wb_write_inv_signed_literal(struct AomWriteBitBuffer *wb, int32_t data,
        int32_t bits);
    //*******************************************************************************************//
    // bitstream.h
    struct AomWriteBitBuffer;

    void write_sequence_header(SequenceControlSet *scs_ptr/*Av1Comp *cpi*/, struct AomWriteBitBuffer *wb);

    uint32_t WriteObuHeader(obuType obuType, int32_t obuExtension,
        uint8_t *const dst);

    int32_t WriteUlebObuSize(uint32_t obuHeaderSize, uint32_t obuPayloadSize,
        uint8_t *dest);

    //*******************************************************************************************//
    // blockd.h
    static INLINE uint32_t have_nearmv_in_inter_mode(PredictionMode mode) {
        return (mode == NEARMV || mode == NEAR_NEARMV || mode == NEAR_NEWMV ||
            mode == NEW_NEARMV);
    }

    void get_txb_ctx(
        SequenceControlSet *sequence_control_set_ptr,
        const int32_t               plane,
        NeighborArrayUnit     *dc_sign_level_coeff_neighbor_array,
        uint32_t                  cu_origin_x,
        uint32_t                  cu_origin_y,
        const BlockSize        plane_bsize,
        const TxSize           tx_size,
        int16_t *const           txb_skip_ctx,
        int16_t *const           dc_sign_ctx);

    extern int32_t eb_av1_get_reference_mode_context(
        uint32_t                  cu_origin_x,
        uint32_t                  cu_origin_y,
        NeighborArrayUnit    *mode_type_neighbor_array,
        NeighborArrayUnit    *inter_pred_dir_neighbor_array);

    extern int32_t eb_av1_get_comp_reference_type_context(
        uint32_t                  cu_origin_x,
        uint32_t                  cu_origin_y,
        NeighborArrayUnit    *mode_type_neighbor_array,
        NeighborArrayUnit     *inter_pred_dir_neighbor_array);

    extern void av1_collect_neighbors_ref_counts(
        CodingUnit            *cu_ptr,
        uint32_t                   cu_origin_x,
        uint32_t                   cu_origin_y,
        NeighborArrayUnit     *mode_type_neighbor_array,
        NeighborArrayUnit     *inter_pred_dir_neighbor_array,
        NeighborArrayUnit     *ref_frame_type_neighbor_array);
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

    extern EbErrorType write_frame_header_av1(
        Bitstream *bitstream_ptr,
        SequenceControlSet *scs_ptr,
        PictureControlSet *pcs_ptr,
        uint8_t showExisting);
    extern EbErrorType encode_td_av1(
        uint8_t *bitstream_ptr);
    extern EbErrorType encode_sps_av1(
        Bitstream *bitstream_ptr,
        SequenceControlSet *scs_ptr);

    //*******************************************************************************************//

    MotionMode motion_mode_allowed(
        const PictureControlSet       *picture_control_set_ptr,
        const CodingUnit              *cu_ptr,
        const BlockSize                 bsize,
        MvReferenceFrame                rf0,
        MvReferenceFrame                rf1,
        PredictionMode                  mode);

    int is_masked_compound_type(COMPOUND_TYPE type);

    static INLINE int32_t is_comp_ref_allowed(BlockSize bsize) {
        return AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
    }

    static INLINE int is_interinter_compound_used(CompoundType type, BlockSize sb_type) {
        const int comp_allowed = is_comp_ref_allowed(sb_type);
        switch (type) {
        case COMPOUND_AVERAGE:
        case COMPOUND_DISTWTD:
        case COMPOUND_DIFFWTD: return comp_allowed;
        case COMPOUND_WEDGE:
            return comp_allowed && wedge_params_lookup[sb_type].bits > 0;
        default: assert(0); return 0;
        }
    }

    static INLINE int is_any_masked_compound_used(BlockSize sb_type) {
        CompoundType comp_type;
        int i;
        if (!is_comp_ref_allowed(sb_type)) return 0;
        for (i = 0; i < COMPOUND_TYPES; i++) {
            comp_type = (CompoundType)i;
            if (is_masked_compound_type(comp_type) &&
                is_interinter_compound_used(comp_type, sb_type))
                return 1;
        }
        return 0;
    }

    static INLINE int32_t tile_log2(int32_t blk_size, int32_t target) {
        int32_t k;
        for (k = 0; (blk_size << k) < target; k++) {
        }
        return k;
}
    void eb_av1_tile_set_col(TileInfo *tile, const TilesInfo *tiles_info,
        int32_t mi_cols, int col);
    void eb_av1_tile_set_row(TileInfo *tile, TilesInfo *tiles_info,
        int32_t mi_rows, int row);
#ifdef __cplusplus
}
#endif
#endif //EbEntropyCoding_h
