/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecObuParser_h
#define EbDecObuParser_h

#include "assert.h"
#include "EbCabacContextModel.h" //for ENTROPY_CONTEXT
#include "EbDecBitstream.h"
#include "EbDecBitReader.h"
#include "EbCodingUnit.h"
#include "EbEntropyCoding.h"

#define HEADER_DUMP 0

#if HEADER_DUMP
#define PRINT_NL printf("\n");
#define PRINT(name, val) printf("\n%s :\t%X", name, val);
#define PRINT_NAME(name) printf("\n%s :\t", name);
#define PRINT_FRAME(name, val) printf("\n%s :\t%X", name, val);
#else
#define PRINT_NL
#define PRINT(name, val)
#define PRINT_NAME(name)
#define PRINT_FRAME(name, val)
#endif

#define ZERO_ARRAY(dest, n) memset(dest, 0, n * sizeof(*(dest)))

#define MAX_NUM_TEMPORAL_LAYERS 8
#define MAX_NUM_SPATIAL_LAYERS 4
#define MAX_NUM_OPERATING_POINTS \
  MAX_NUM_TEMPORAL_LAYERS * MAX_NUM_SPATIAL_LAYERS
#define OP_POINTS_CNT_MINUS_1_BITS 5
#define OP_POINTS_IDC_BITS 12

#define SELECT_SCREEN_CONTENT_TOOLS 2
#define SELECT_INTEGER_MV 2
#define NUM_REF_FRAMES 8

enum {
  SEQ_LEVEL_2_0,
  SEQ_LEVEL_2_1,
  SEQ_LEVEL_2_2,
  SEQ_LEVEL_2_3,
  SEQ_LEVEL_3_0,
  SEQ_LEVEL_3_1,
  SEQ_LEVEL_3_2,
  SEQ_LEVEL_3_3,
  SEQ_LEVEL_4_0,
  SEQ_LEVEL_4_1,
  SEQ_LEVEL_4_2,
  SEQ_LEVEL_4_3,
  SEQ_LEVEL_5_0,
  SEQ_LEVEL_5_1,
  SEQ_LEVEL_5_2,
  SEQ_LEVEL_5_3,
  SEQ_LEVEL_6_0,
  SEQ_LEVEL_6_1,
  SEQ_LEVEL_6_2,
  SEQ_LEVEL_6_3,
  SEQ_LEVEL_7_0,
  SEQ_LEVEL_7_1,
  SEQ_LEVEL_7_2,
  SEQ_LEVEL_7_3,
  SEQ_LEVELS,
  SEQ_LEVEL_MAX = 31
} UENUM1BYTE(AV1_LEVEL);


typedef struct ParseNbr4x4Ctxt {

    /* Buffer holding the transform sizes of the previous 4x4 block row. */
    uint8_t *above_tx_wd;

    /* Buffer holding the transform sizes of the left 4x4 blocks corresponding
     to the current super block row. */
    uint8_t *left_tx_ht;

    /* Buffer holding the partition context of the previous 4x4 block row. */
    uint8_t *above_part_wd;

    /* Buffer holding the partition context of the left 4x4 blocks corresponding
     to the current super block row. */
    uint8_t *left_part_ht;

    /* Buffer holding the sign of the DC coefficients and the cumulative
       sum of the coefficient levels of the previous 4x4 block row. */
    int8_t *above_ctx[MAX_MB_PLANE];

    /* Buffer holding the sign of the DC coefficients and the cumulative
       sum of the coefficient levels of the left 4x4 blocks
       corresponding to the current super block row. */
    int8_t *left_ctx[MAX_MB_PLANE];

    /* Number of mi columns with respect to the aligned width. */
    uint32_t num_mi_col;

    /* Buffer holding the seg_id_predicted of the previous 4x4 block row. */
    uint8_t *left_seg_pred_ctx;

    /* Buffer holding the seg_id_predicted of the left 4x4 blocks corresponding
     to the current super block row. */
    uint8_t *above_seg_pred_ctx;

    /* Value of base colors for Y, U, and V */
    uint16_t *above_palette_colors[MAX_MB_PLANE];
    uint16_t *left_palette_colors[MAX_MB_PLANE];

    /* Buffer holding the delta LF values*/
    int32_t delta_lf[FRAME_LF_COUNT];

    /* Place holder for the current q index*/
    int32_t cur_q_ind;

    /* Place holder for palette color information */
    uint16_t palette_colors[MAX_MB_PLANE][PALETTE_MAX_SIZE];

    int8_t *above_comp_grp_idx;

    int8_t *left_comp_grp_idx;

} ParseNbr4x4Ctxt;

typedef struct ParseCtxt {
    /** Decoder Handle */
    void *dec_handle_ptr;

    /** bitstream engine Handle */
    bitstrm_t   *bs;

    /** symbol decoder Handle */
    SvtReader   r;

    ParseNbr4x4Ctxt parse_nbr4x4_ctxt;

    //FRAME_CONTEXT   frm_ctx[DEC_MAX_NUM_FRM_PRLL];
    FRAME_CONTEXT   cur_tile_ctx;
    FRAME_CONTEXT   init_frm_ctx;

    TileInfo        cur_tile_info;

    /* Stored here for current block and should be updated to next block modeinfo */
    /*!< Offset of first transform info from strat of SB pointer for each plane */
    uint16_t        first_tu_offset[MAX_MB_PLANE - 1];
    /* TODO: Points to the cur ModeInfo_t in SB. Should be moved out */
    BlockModeInfo   *cur_mode_info;
    /* TODO: Points to the cur ModeInfo_t in SB. Should be moved out */
    int32_t         cur_mode_info_cnt;
    /* TODO: cur SB row idx. Should be moved out */
    int32_t         sb_row_mi;
    /* TODO: cur SB col idx. Should be moved out */
    int32_t         sb_col_mi;

    /* TODO: Points to the cur coeff_buf in SB. Should be moved out */
    int32_t *cur_coeff_buf[MAX_MB_PLANE];

    /* Points to the cur luma_trans_info in a block */
    TransformInfo_t *cur_luma_trans_info;

    /* Count of cur luma_trans_info */
    uint8_t cur_blk_luma_count;
#if !FRAME_MI_MAP
    /* Left and above SBInfo pointers */
    SBInfo  *left_sb_info;
    SBInfo  *above_sb_info;
#endif

    TransformInfo_t *inter_trans_chroma;

    /*!< Number of TUs in block or force split block */
    uint8_t         num_tus[MAX_MB_PLANE][4 /*Max force TU split*/];

    /*!< Reference Loop Restoration Unit  */
    RestorationUnitInfo ref_lr_unit[MAX_MB_PLANE];

    EbBool  read_deltas;
} ParseCtxt;

int get_qindex(SegmentationParams *seg_params, int segment_id, int base_q_idx);
void parse_super_block(EbDecHandle *dec_handle, uint32_t blk_row,
                       uint32_t blk_col, SBInfo *sbInfo);

void svt_setup_motion_field(EbDecHandle *dec_handle);

EbErrorType decode_obu(EbDecHandle *dec_handle_ptr, uint8_t *data, uint32_t data_size);
EbErrorType decode_multiple_obu(EbDecHandle *dec_handle_ptr, uint8_t **data,
    size_t data_size, uint32_t is_annexb);

static INLINE int allow_intrabc(const EbDecHandle *dec_handle) {
    return  (dec_handle->frame_header.frame_type == KEY_FRAME
        || dec_handle->frame_header.frame_type == INTRA_ONLY_FRAME)
        && dec_handle->seq_header.seq_force_screen_content_tools
        && dec_handle->frame_header.allow_intrabc;
}

#endif  // EbDecObuParser_h
