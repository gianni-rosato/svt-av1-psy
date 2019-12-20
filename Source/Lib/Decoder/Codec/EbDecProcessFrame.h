/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecProcessFrame_h
#define EbDecProcessFrame_h

#ifdef __cplusplus
extern "C" {
#endif

#include "../../Encoder/Codec/EbIntraPrediction.h"
#include "EbDecObmc.h"

typedef struct DecModCtxt {
    /** Decoder Handle */
    void *dec_handle_ptr;

    SeqHeader *seq_header;

    FrameHeader *frame_header;

    int32_t *sb_iquant_ptr;

    int32_t *iquant_cur_ptr;

#if !FRAME_MI_MAP
    /* TODO: cur SB row idx. Should be moved out */
    int32_t         sb_row_mi;
    /* TODO: cur SB col idx. Should be moved out */
    int32_t         sb_col_mi;

    /* Left and above SBInfo pointers */
    SBInfo  *left_sb_info;
    SBInfo  *above_sb_info;
#endif
    /* TODO: Points to the cur coeff_buf in SB */
    int32_t *cur_coeff[MAX_MB_PLANE];

    /* Current tile info */
    TileInfo    cur_tile_info;

    /* CFL context */
    CflCtx  cfl_ctx;

    /*OBMC context*/
    ObmcCtx obmc_ctx;

    /* TODO: IntraRef Scratch buf! Should be moved to thrd ctxt */
    uint16_t    topNeighArray[64 * 2 + 1];
    uint16_t    leftNeighArray[64 * 2 + 1];

    /* Dequantization context */
    Dequant                dequants;

    /* This need to be moved to thread context */
    Dequant                *dequants_delta_q;

    /* Inverse Quantization Matrix */
    const QmVal          *giqmatrix[NUM_QM_LEVELS][3][TX_SIZES_ALL];

    /*Mask for Comp mode blending*/
    DECLARE_ALIGNED(16, uint8_t, seg_mask[2 * MAX_SB_SQUARE]);
} DecModCtxt;

typedef struct LRCtxt {
    /** Decoder Handle */
    void *dec_handle_ptr;

    /* Wiener and SGR Filter holder */
    RestorationUnitInfo    *lr_unit[MAX_MB_PLANE];

    int32_t lr_stride[MAX_MB_PLANE];

    /* Buffer to store deblocked line buffer around stripe boundary */
    RestorationStripeBoundaries boundaries[MAX_MB_PLANE];

    /* Used to store CDEF line buffer around stripe boundary */
    RestorationLineBuffers *rlbs;

    /* Scratch buffer to hold LR output */
    uint8_t *dst;
    uint16_t dst_stride;

    /* Pointer to a scratch buffer used by self-guided restoration */
    int32_t *rst_tmpbuf;
}LRCtxt;


void decode_super_block(DecModCtxt *dec_mod_ctxt,
                        uint32_t mi_row, uint32_t mi_col,
                        SBInfo *sbInfo);

EbErrorType start_decode_tile(EbDecHandle *dec_handle_ptr,
    DecModCtxt *dec_mod_ctxt, TilesInfo *tiles_info, int32_t tile_num);

#if !ENABLE_ROW_MT_DECODE
EbErrorType decode_tile(
    DecModCtxt *dec_mod_ctxt, TilesInfo *tile_info,
    DecMTParseReconTileInfo *parse_recon_tile_info_array,
    int32_t tile_row, int32_t tile_col);
#else
EbErrorType decode_tile(
    DecModCtxt *dec_mod_ctxt, TilesInfo *tile_info,
    DecMTParseReconTileInfo *parse_recon_tile_info_array,
    int32_t tile_col);
#endif

/* TODO: Should be moved out once decode tile is moved out from parse_tile */
void cfl_init(CflCtx *cfl, EbColorConfig *cc);

#ifdef __cplusplus
}
#endif

#endif // EbDecProcessFrame_h
