/*
* Copyright(c) 2019 Netflix, Inc.
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbDecParseFrame_h
#define EbDecParseFrame_h

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>

typedef struct ParseAboveNbr4x4Ctxt {
    /* Buffer holding the transform sizes of the previous 4x4 block row. */
    uint8_t *above_tx_wd;

    /* Buffer holding the partition context of the previous 4x4 block row. */
    uint8_t *above_part_wd;

    /* Buffer holding the sign of the DC coefficients and the cumulative sum of
       the coefficient levels of the previous 4x4 block row. */
    uint8_t *above_ctx[MAX_MB_PLANE];

    /* Buffer holding the seg_id_predicted of the left 4x4 blocks corresponding
       to the current super block row. */
    uint8_t *above_seg_pred_ctx;

    /* Value of base colors for Y, U, and V */
    uint16_t *above_palette_colors[MAX_MB_PLANE];

    int8_t *above_comp_grp_idx;

} ParseAboveNbr4x4Ctxt;

typedef struct ParseLeftNbr4x4Ctxt {
    /* Buffer holding the transform sizes of the left 4x4 blocks corresponding
     to the current super block row. */
    uint8_t *left_tx_ht;

    /* Buffer holding the partition context of the left 4x4 blocks corresponding
     to the current super block row. */
    uint8_t *left_part_ht;

    /* Buffer holding the sign of the DC coefficients and the
       cumulative sum of the coefficient levels of the left 4x4
       blocks corresponding to the current super block row. */
    uint8_t *left_ctx[MAX_MB_PLANE];

    /* Buffer holding the seg_id_predicted of the previous 4x4 block row. */
    uint8_t *left_seg_pred_ctx;

    /* Value of base colors for Y, U, and V */
    uint16_t *left_palette_colors[MAX_MB_PLANE];

    int8_t *left_comp_grp_idx;
} ParseLeftNbr4x4Ctxt;

/* Bhavna : Add comment */
typedef struct ParseTileData {
    uint8_t *data;
    uint8_t *data_end;
    size_t   tile_size;
} ParseTileData;

typedef struct ParseCtxt {
    /** Decoder Handle */
    void *dec_handle_ptr;

    /** symbol decoder Handle */
    SvtReader r;

    SeqHeader *  seq_header;
    FrameHeader *frame_header;

    ParseAboveNbr4x4Ctxt *parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * parse_left_nbr4x4_ctxt;

    //FRAME_CONTEXT   frm_ctx[DEC_MAX_NUM_FRM_PRLL];
    FRAME_CONTEXT cur_tile_ctx;

    TileInfo cur_tile_info;

    /* Stored here for current block and should be updated to next block modeinfo */
    /*!< Offset of first Luma and Chroma transform info from strat of SB pointer */
    uint16_t first_txb_offset[MAX_MB_PLANE - 1];

    /* TODO: Points to the cur ModeInfo_t in SB. Should be moved out */
    BlockModeInfo *cur_mode_info;
    /* TODO: Points to the cur ModeInfo_t in SB. Should be moved out */
    int32_t cur_mode_info_cnt;
    /* TODO: cur SB row idx. Should be moved out */
    int32_t sb_row_mi;
    /* TODO: cur SB col idx. Should be moved out */
    int32_t sb_col_mi;

    /* TODO: Points to the cur coeff_buf in SB. Should be moved out */
    int32_t *cur_coeff_buf[MAX_MB_PLANE];

    /* Points to the cur luma_trans_info in a block */
    TransformInfo_t *cur_luma_trans_info;

    /* Count of cur luma_trans_info */
    uint8_t cur_blk_luma_count;
    /*!< Chroma mode indo state acroos sub 8x8 blocks
     * if Prev block does not have chroma info then this state is remembered in this variable to be used in next block
    */
    int32_t prev_blk_has_chroma;

    TransformInfo_t *inter_trans_chroma;

    /*!< Number of TUs in block or force split block */
    uint8_t num_tus[MAX_MB_PLANE][4 /*Max force TU split*/];

    /*!< Reference Loop Restoration Unit  */
    RestorationUnitInfo ref_lr_unit[MAX_MB_PLANE];

    EbBool read_deltas;

    /* Buffer holding the delta LF values*/
    int32_t delta_lf[FRAME_LF_COUNT];

    /* Place holder for palette color information */
    uint16_t palette_colors[MAX_MB_PLANE][PALETTE_MAX_SIZE];

    /* Place holder for the current q index*/
    int32_t cur_q_ind;
} ParseCtxt;

typedef struct MasterParseCtxt {
    /* Array holding the context for each of the tile.*/
    ParseCtxt *tile_parse_ctxt;

    /* Frame level above neighbour parse context.*/
    ParseAboveNbr4x4Ctxt *parse_above_nbr4x4_ctxt;

    /* Frame level left neighbour parse context.*/
    ParseLeftNbr4x4Ctxt *parse_left_nbr4x4_ctxt;

    /* Frame level CDF value holder.*/
    FRAME_CONTEXT init_frm_ctx;

    /* Curent number of instances of tile context.*/
    int32_t context_count;

    /* Curent number of Tiles.*/
    int32_t num_tiles;

    /* Array of ParseTileData for each Tile */
    ParseTileData *parse_tile_data;
} MasterParseCtxt;

void parse_super_block(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, uint32_t blk_row,
                       uint32_t blk_col, SBInfo *sb_info);

void svt_tile_init(TileInfo *cur_tile_info, FrameHeader *frame_header, int32_t tile_row,
                   int32_t tile_col);

static int read_is_valid(const uint8_t *start, size_t len, const uint8_t *end) {
    return len != 0 && len <= (size_t)(end - start);
}

static INLINE EbErrorType init_svt_reader(SvtReader *r, const uint8_t *data, const uint8_t *data_end,
                            const size_t read_size, uint8_t allow_update_cdf) {
    if (read_is_valid(data, read_size, data_end) && !svt_reader_init(r, data, read_size))
        r->allow_update_cdf = allow_update_cdf;
    else
        return EB_Corrupt_Frame;
    return EB_ErrorNone;
}

EbErrorType start_parse_tile(EbDecHandle *dec_handle_ptr, ParseCtxt *parse_ctxt,
                             TilesInfo *tiles_info, int tile_num, int is_mt);

EbErrorType parse_tile(EbDecHandle *dec_handle_ptr, ParseCtxt *parse_ctx, TilesInfo *tile_info,
                       int tile_num, int32_t tile_row, int32_t tile_col, int32_t is_mt);

#ifdef __cplusplus
}
#endif

#endif // EbDecParseFrame_h
