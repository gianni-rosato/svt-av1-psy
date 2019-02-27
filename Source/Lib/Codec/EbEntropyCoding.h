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

#ifdef __cplusplus
extern "C" {
#endif
    /*!\brief OBU types. */
    typedef enum ATTRIBUTE_PACKED {
        OBU_SEQUENCE_HEADER = 1,
        OBU_TEMPORAL_DELIMITER = 2,
        OBU_FRAME_HEADER = 3,
        OBU_TILE_GROUP = 4,
        OBU_METADATA = 5,
        OBU_FRAME = 6,
        OBU_REDUNDANT_FRAME_HEADER = 7,
        OBU_PADDING = 15,
    } obuType;

    /**************************************
     * Extern Function Declarations
     **************************************/
    struct EntropyCodingContext_s;
    extern EbErrorType write_sb(
        struct EntropyCodingContext_s   *context_ptr,
        LargestCodingUnit_t     *tbPtr,
        PictureControlSet_t     *picture_control_set_ptr,
        EntropyCoder_t          *entropy_coder_ptr,
        EbPictureBufferDesc_t   *coeffPtr);


    extern EbErrorType EncodeSliceFinish(
        EntropyCoder_t        *entropy_coder_ptr);

    extern EbErrorType ResetBitstream(
        EbPtr                 bitstreamPtr);

    extern EbErrorType ResetEntropyCoder(
        EncodeContext_t       *encode_context_ptr,
        EntropyCoder_t        *entropy_coder_ptr,
        uint32_t                 qp,
        EB_SLICE               slice_type);

    extern EbErrorType Av1TuEstimateCoeffBits(
        PictureControlSet_t                    *picture_control_set_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        CodingUnit_t                           *cu_ptr,
        uint32_t                                  tuOriginIndex,
        uint32_t                                  tuChromaOriginIndex,
        EntropyCoder_t                         *entropy_coder_ptr,
        EbPictureBufferDesc_t                  *coeff_buffer_sb,
        uint32_t                                 yEob,
        uint32_t                                 cbEob,
        uint32_t                                 crEob,
        uint64_t                                 *y_tu_coeff_bits,
        uint64_t                                 *cb_tu_coeff_bits,
        uint64_t                                 *cr_tu_coeff_bits,
        TxSize                                 txsize,
        TxSize                                 txsize_uv,
        COMPONENT_TYPE                          component_type,
        EbAsm                                  asm_type);

    extern EbErrorType CopyRbspBitstreamToPayload(
        Bitstream_t *bitstreamPtr,
        EbByte      outputBuffer,
        uint32_t      *outputBufferIndex,
        uint32_t      *outputBufferSize,
        EncodeContext_t         *encode_context_ptr);


    //**********************************************************************************************************//
    //onyxc_int.h
    static INLINE int32_t frame_is_intra_only(const PictureParentControlSet_t *const pcsPtr) {
        return pcsPtr->av1FrameType == KEY_FRAME || pcsPtr->av1FrameType == INTRA_ONLY_FRAME;
    }

    static INLINE int32_t frame_is_sframe(const PictureParentControlSet_t *pcsPtr) {
        return pcsPtr->av1FrameType == S_FRAME;
    }

    // Returns 1 if this frame might allow mvs from some reference frame.

    static INLINE int32_t frame_might_allow_ref_frame_mvs(const PictureParentControlSet_t *pcsPtr,
        SequenceControlSet_t    *scsPtr) {
#if AV1_UPGRADE
        return !pcsPtr->error_resilient_mode &&
#else
        return !pcsPtr->error_resilient_mode && !pcsPtr->large_scale_tile &&
#endif
            scsPtr->enable_ref_frame_mvs &&
            scsPtr->enable_order_hint && !frame_is_intra_only(pcsPtr);
    }

    // Returns 1 if this frame might use warped_motion
    static INLINE int32_t frame_might_allow_warped_motion(const PictureParentControlSet_t *pcsPtr,
        SequenceControlSet_t    *scsPtr) {
        return !pcsPtr->error_resilient_mode && !frame_is_intra_only(pcsPtr) &&
            scsPtr->static_config.enable_warped_motion;
    }

    static INLINE uint8_t major_minor_to_seq_level_idx(BitstreamLevel bl) {
        assert(bl.major >= LEVEL_MAJOR_MIN && bl.major <= LEVEL_MAJOR_MAX);
        //assert(bl.minor >= LEVEL_MINOR_MIN && bl.minor <= LEVEL_MINOR_MAX);
        return ((bl.major - LEVEL_MAJOR_MIN) << LEVEL_MINOR_BITS) + bl.minor;
    }

    //**********************************************************************************************************//
    //encoder.h
    static INLINE int32_t get_ref_frame_map_idx(const PictureParentControlSet_t *pcsPtr,
        MvReferenceFrame ref_frame) {
        // (void)(*pcsPtr);
        // (void)ref_frame;
        // return 0;

        return pcsPtr->av1RefSignal.refDpbIndex[ref_frame - LAST_FRAME];//LAST-LAST2-LAST3-GOLDEN-BWD-ALT2-ALT
        //if (ref_frame >= LAST_FRAME && ref_frame <= LAST3_FRAME)
        //    return pcsPtr->lst_fb_idxes[ref_frame - 1];
        //else if (ref_frame == GOLDEN_FRAME)
        //    return pcsPtr->gld_fb_idx;
        //else if (ref_frame == BWDREF_FRAME)
        //    return pcsPtr->bwd_fb_idx;
        //else if (ref_frame == ALTREF2_FRAME)
        //    return pcsPtr->alt2_fb_idx;
        //else
        //    return pcsPtr->alt_fb_idx;
    }

    //*******************************************************************************************//
    // bitwriter_buffer.h
    struct aom_write_bit_buffer {
        uint8_t *bit_buffer;
        uint32_t bit_offset;
    };

    int32_t aom_wb_is_byte_aligned(const struct aom_write_bit_buffer *wb);
    uint32_t aom_wb_bytes_written(const struct aom_write_bit_buffer *wb);

    void aom_wb_write_bit(struct aom_write_bit_buffer *wb, int32_t bit);

    void aom_wb_overwrite_bit(struct aom_write_bit_buffer *wb, int32_t bit);

    void aom_wb_write_literal(struct aom_write_bit_buffer *wb, int32_t data, int32_t bits);

    void aom_wb_write_inv_signed_literal(struct aom_write_bit_buffer *wb, int32_t data,
        int32_t bits);
    //*******************************************************************************************//
    // bitstream.h
    struct aom_write_bit_buffer;

    //void WriteSequenceHeader(/*AV1_COMP *cpi, */struct aom_write_bit_buffer *wb);
    void WriteSequenceHeader(SequenceControlSet_t *scsPtr/*AV1_COMP *cpi*/, struct aom_write_bit_buffer *wb);

    uint32_t WriteObuHeader(obuType obuType, int32_t obuExtension,
        uint8_t *const dst);

    int32_t WriteUlebObuSize(uint32_t obuHeaderSize, uint32_t obuPayloadSize,
        uint8_t *dest);

    /*int32_t av1_pack_bitstream(AV1_COMP *const cpi, uint8_t *dest, size_t *size);

    static INLINE int32_t av1_preserve_existing_gf(AV1_COMP *cpi) {
    // Do not swap gf and arf indices for internal overlay frames
    return !cpi->multi_arf_allowed && cpi->rc.is_src_frame_alt_ref &&
    !cpi->rc.is_src_frame_ext_arf;
    }

    void av1_write_tx_type(const Av1Common *const cm, const MacroBlockD *xd,
    int32_t blk_row, int32_t blk_col, int32_t plane, TxSize tx_size,
    aom_writer *w);
    */

    //*******************************************************************************************//
    // blockd.h
    static INLINE uint32_t have_nearmv_in_inter_mode(PredictionMode mode) {
        return (mode == NEARMV || mode == NEAR_NEARMV || mode == NEAR_NEWMV ||
            mode == NEW_NEARMV);
    }

    void GetTxbCtx(
        const int32_t               plane,
        NeighborArrayUnit_t     *dcSignLevelCoeffNeighborArray,
        uint32_t                  cu_origin_x,
        uint32_t                  cu_origin_y,
        const block_size        plane_bsize,
        const TxSize           tx_size,
        int16_t *const           txb_skip_ctx,
        int16_t *const           dc_sign_ctx);

    extern int32_t Av1GetReferenceModeContext(
        uint32_t                  cu_origin_x,
        uint32_t                  cu_origin_y,
        NeighborArrayUnit_t    *mode_type_neighbor_array,
        NeighborArrayUnit_t    *inter_pred_dir_neighbor_array);

    extern int32_t Av1GetCompReferenceTypeContext(
        uint32_t                  cu_origin_x,
        uint32_t                  cu_origin_y,
        NeighborArrayUnit_t    *mode_type_neighbor_array,
        NeighborArrayUnit_t     *inter_pred_dir_neighbor_array);

    extern void Av1CollectNeighborsRefCounts(
        CodingUnit_t            *cu_ptr,
        uint32_t                   cu_origin_x,
        uint32_t                   cu_origin_y,
        NeighborArrayUnit_t     *mode_type_neighbor_array,
        NeighborArrayUnit_t     *inter_pred_dir_neighbor_array,
        NeighborArrayUnit_t     *ref_frame_type_neighbor_array);

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
    extern int32_t av1_get_pred_context_comp_ref_p(const MacroBlockD *xd);

    // Returns a context number for the given MB prediction signal
    // Signal the first reference frame for a compound mode be LAST,
    // conditioning on that it is known either LAST/LAST2.
    extern int32_t av1_get_pred_context_comp_ref_p1(const MacroBlockD *xd);

    // Returns a context number for the given MB prediction signal
    // Signal the first reference frame for a compound mode be GOLDEN,
    // conditioning on that it is known either GOLDEN or LAST3.
    extern int32_t av1_get_pred_context_comp_ref_p2(const MacroBlockD *xd);

    // Signal the 2nd reference frame for a compound mode be either
    // ALTREF, or ALTREF2/BWDREF.
    extern int32_t av1_get_pred_context_comp_bwdref_p(const MacroBlockD *xd);

    // Signal the 2nd reference frame for a compound mode be either
    // ALTREF2 or BWDREF.
    extern int32_t av1_get_pred_context_comp_bwdref_p1(const MacroBlockD *xd);
    // == Context functions for single ref ==
    //
    // For the bit to signal whether the single reference is a forward reference
    // frame or a backward reference frame.
    extern int32_t av1_get_pred_context_single_ref_p1(const MacroBlockD *xd);

    // For the bit to signal whether the single reference is ALTREF_FRAME or
    // non-ALTREF backward reference frame, knowing that it shall be either of
    // these 2 choices.
    extern int32_t av1_get_pred_context_single_ref_p2(const MacroBlockD *xd);

    // For the bit to signal whether the single reference is LAST3/GOLDEN or
    // LAST2/LAST, knowing that it shall be either of these 2 choices.
    extern int32_t av1_get_pred_context_single_ref_p3(const MacroBlockD *xd);

    // For the bit to signal whether the single reference is LAST2_FRAME or
    // LAST_FRAME, knowing that it shall be either of these 2 choices.
    extern int32_t av1_get_pred_context_single_ref_p4(const MacroBlockD *xd);

    // For the bit to signal whether the single reference is GOLDEN_FRAME or
    // LAST3_FRAME, knowing that it shall be either of these 2 choices.
    extern int32_t av1_get_pred_context_single_ref_p5(const MacroBlockD *xd);

    // For the bit to signal whether the single reference is ALTREF2_FRAME or
    // BWDREF_FRAME, knowing that it shall be either of these 2 choices.
    extern int32_t av1_get_pred_context_single_ref_p6(const MacroBlockD *xd);


    extern EbErrorType WriteFrameHeaderAv1(
        Bitstream_t *bitstreamPtr,
        SequenceControlSet_t *scsPtr,
        PictureControlSet_t *pcsPtr,
        uint8_t showExisting);
    extern EbErrorType encode_td_av1(
        uint8_t *bitstreamPtr);
    extern EbErrorType EncodeSPSAv1(
        Bitstream_t *bitstreamPtr,
        SequenceControlSet_t *scsPtr);

    //*******************************************************************************************//

    MOTION_MODE motion_mode_allowed(
        const PictureControlSet_t       *picture_control_set_ptr,
        const CodingUnit_t              *cu_ptr,
        const block_size                 bsize,
        MvReferenceFrame                rf0,
        MvReferenceFrame                rf1,
        PredictionMode                  mode);

#ifdef __cplusplus
}
#endif
#endif //EbEntropyCoding_h