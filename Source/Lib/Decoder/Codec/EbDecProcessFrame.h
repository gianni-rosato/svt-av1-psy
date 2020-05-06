/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecProcessFrame_h
#define EbDecProcessFrame_h

#ifdef __cplusplus
extern "C" {
#endif

#include "EbIntraCommon.h"
#include "EbDecObmc.h"

typedef struct DecModCtxt {
    /** Decoder Handle */
    void *dec_handle_ptr;

    SeqHeader *seq_header;

    FrameHeader *frame_header;

    int32_t *sb_iquant_ptr;

    int32_t *iquant_cur_ptr;

    /* TODO: Points to the cur coeff_buf in SB */
    int32_t *cur_coeff[MAX_MB_PLANE];

    /* Current tile info */
    TileInfo cur_tile_info;

    /* CFL context */
    CflCtx cfl_ctx;

    /*OBMC context*/
    ObmcCtx obmc_ctx;

    /* TODO: IntraRef Scratch buf! Should be moved to thrd ctxt */
    uint16_t top_neigh_array[64 * 2 + 1];
    uint16_t left_neigh_array[64 * 2 + 1];

    /* Dequantization context */
    Dequant dequants;

    /* This need to be moved to thread context */
    Dequant *dequants_delta_q;

    /* Inverse Quantization Matrix */
    const QmVal *giqmatrix[NUM_QM_LEVELS][3][TX_SIZES_ALL];

    /*Mask for Comp mode blending*/
    DECLARE_ALIGNED(16, uint8_t, seg_mask[2 * MAX_SB_SQUARE]);
    /*MC temp buff for dynamic padding*/
    uint8_t *mc_buf[2];
} DecModCtxt;

typedef struct LrCtxt {
    /** Decoder Handle */
    void *dec_handle_ptr;

    /* Wiener and SGR Filter holder */
    RestorationUnitInfo *lr_unit[MAX_MB_PLANE];

    int32_t lr_stride[MAX_MB_PLANE];

    /* Buffer to store deblocked line buffer around stripe boundary */
    RestorationStripeBoundaries boundaries[MAX_MB_PLANE];

    /* Used to store CDEF line buffer around stripe boundary */
    RestorationLineBuffers ***rlbs;

    /* Scratch buffer to hold LR output */
    //uint8_t *dst;
    //uint16_t dst_stride;

    /* Temporary block level scratch buffer to store
       LR output of [SB_Size x 64] block */
    uint8_t *dst;

    /* Pointer to a scratch buffer used by self-guided restoration */
    int32_t **rst_tmpbuf;

    /* Flag to indicate if the access of buffers must be
       based on thread index or SB row index */
    EbBool is_thread_min;

} LrCtxt;


void decode_super_block(DecModCtxt *dec_mod_ctxt, uint32_t mi_row, uint32_t mi_col,
                        SBInfo *sb_info);

EbErrorType start_decode_tile(EbDecHandle *dec_handle_ptr, DecModCtxt *dec_mod_ctxt,
                              TilesInfo *tiles_info, int32_t tile_num);
EbErrorType decode_tile(DecModCtxt *dec_mod_ctxt, TilesInfo *tile_info,
                        DecMtParseReconTileInfo *parse_recon_tile_info_array, int32_t tile_col);

/* TODO: Should be moved out once decode tile is moved out from parse_tile */
void cfl_init(CflCtx *cfl, EbColorConfig *cc);

#ifdef __cplusplus
}
#endif

#endif // EbDecProcessFrame_h
