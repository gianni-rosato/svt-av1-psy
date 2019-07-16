/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecProcessFrame_h
#define EbDecProcessFrame_h

#ifdef __cplusplus
extern "C" {
#endif

#include "EbIntraPrediction.h"

typedef struct DecModCtxt {
    /** Decoder Handle */
    void *dec_handle_ptr;

    int32_t *sb_iquant_ptr;
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
    TileInfo    *cur_tile_info;

    /* CFL context */
    CflCtx  cfl_ctx;

    /* TODO: IntraRef Scratch buf! Should be moved to thrd ctxt */
    uint16_t    topNeighArray[64 * 2 + 1];
    uint16_t    leftNeighArray[64 * 2 + 1];

    /* Dequantization context */
    Dequant                dequants;

    /* This need to be moved to thread context */
    Dequant                *dequants_delta_q;

    /* Inverse Quantization Matrix */
    const QmVal          *giqmatrix[NUM_QM_LEVELS][3][TX_SIZES_ALL];

} DecModCtxt;

void decode_super_block(DecModCtxt *dec_mod_ctxt,
                        uint32_t mi_row, uint32_t mi_col,
                        SBInfo *sbInfo);

/* TODO: Should be moved out once decode tile is moved out from parse_tile */
void cfl_init(CflCtx *cfl, EbColorConfig *cc);

#ifdef __cplusplus
}
#endif

#endif // EbDecProcessFrame_h
