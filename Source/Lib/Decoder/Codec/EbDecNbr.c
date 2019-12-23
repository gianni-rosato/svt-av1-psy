/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbObuParse.h"
#include "EbDecParseFrame.h"
#include "EbDecNbr.h"

void update_block_nbrs(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, int mi_row, int mi_col,
                       BlockSize subsize) {
    FrameMiMap *frame_mi_map = &dec_handle->master_frame_buf.frame_mi_map;

    int32_t offset = parse_ctx->cur_mode_info_cnt;
    // int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;
    int bw4 = mi_size_wide[subsize];
    int bh4 = mi_size_high[subsize];
    /*TODO : Can remove later*/
    assert(mi_row >= 0);
    assert(mi_row + bh4 <= frame_mi_map->mi_rows_algnsb);
    assert(mi_col >= 0);
    assert(mi_col + bw4 <= frame_mi_map->mi_cols_algnsb);
    /* Update 4x4 nbr offset map */
    for (int i = mi_row; i < mi_row + bh4; i++) {
        uint16_t *p_mi_offset = frame_mi_map->p_mi_offset + i * frame_mi_map->mi_cols_algnsb;
        for (int j = mi_col; j < mi_col + bw4; j++) p_mi_offset[j] = offset;
    }
}
/* TODO : Should remove dec_mod_ctxt dependency */
BlockModeInfo *get_cur_mode_info(void *pv_dec_handle, int mi_row, int mi_col, SBInfo *sb_info) {
    EbDecHandle *dec_handle   = (EbDecHandle *)pv_dec_handle;
    FrameMiMap * frame_mi_map = &dec_handle->master_frame_buf.frame_mi_map;

    BlockModeInfo *cur_mi = NULL;
    (void)sb_info;
    int32_t cur_sb_row = mi_row >> (frame_mi_map->sb_size_log2 - MI_SIZE_LOG2);
    int32_t cur_sb_col = mi_col >> (frame_mi_map->sb_size_log2 - MI_SIZE_LOG2);
    SBInfo *cur_sb_info =
        frame_mi_map->pps_sb_info[cur_sb_row * frame_mi_map->sb_cols + cur_sb_col];
    int32_t offset = *(frame_mi_map->p_mi_offset + mi_row * frame_mi_map->mi_cols_algnsb + mi_col);
    cur_mi         = &cur_sb_info->sb_mode_info[offset];
    return cur_mi;
}

/* TODO : Should remove parse_ctx dependency */
BlockModeInfo *get_left_mode_info(EbDecHandle *dec_handle, int mi_row, int mi_col,
                                  SBInfo *sb_info) {
    (void)sb_info;
    return get_cur_mode_info(dec_handle, mi_row, mi_col - 1, NULL);
}

/* TODO : Should remove parse_ctx dependency */
BlockModeInfo *get_top_mode_info(EbDecHandle *dec_handle, int mi_row, int mi_col, SBInfo *sb_info) {
    (void)sb_info;
    return get_cur_mode_info(dec_handle, mi_row - 1, mi_col, NULL);
}
