/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbObuParse.h"
#include "EbDecProcessFrame.h"

#include "EbDecNbr.h"
#if !FRAME_MI_MAP
void update_nbrs_before_sb(FrameMiMap *frame_mi_map, int32_t sb_col) {
    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    /* Update from top sbrow to current SB top */
    memcpy(&frame_mi_map->cur_sb_mi_map[0][1], &frame_mi_map->top_sbrow_mi_map[sb_col*num_mis_in_sb_wd],
        num_mis_in_sb_wd * sizeof(int16_t));
}

void update_nbrs_after_sb(FrameMiMap *frame_mi_map, int32_t sb_col) {
    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    /* Update cur right 4x4 col as next left 4x4 */
    for (int i = 1; i < num_mis_in_sb_wd + 1; i++)
        frame_mi_map->cur_sb_mi_map[i][0] = frame_mi_map->cur_sb_mi_map[i][num_mis_in_sb_wd];
    /* Update bottom 4x4 of cur as top 4x4 for next row */
    memcpy(&frame_mi_map->top_sbrow_mi_map[sb_col*num_mis_in_sb_wd], &frame_mi_map->cur_sb_mi_map[num_mis_in_sb_wd][1],
        num_mis_in_sb_wd * sizeof(int16_t));
}
#endif

void update_block_nbrs(EbDecHandle *dec_handle,
    ParseCtxt   *parse_ctx, int mi_row, int mi_col,
    BlockSize subsize)
{
    FrameMiMap  *frame_mi_map = &dec_handle->master_frame_buf.frame_mi_map;

    int32_t offset = parse_ctx->cur_mode_info_cnt;
    // int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;
    int bw4 = mi_size_wide[subsize];
    int bh4 = mi_size_high[subsize];
#if FRAME_MI_MAP
    /*TODO : Can remove later*/
    assert(mi_row >= 0); assert(mi_row+ bh4 <= frame_mi_map->mi_rows_algnsb);
    assert(mi_col >= 0); assert(mi_col + bw4 <= frame_mi_map->mi_cols_algnsb);
    /* Update 4x4 nbr offset map */
    for (int i = mi_row; i < mi_row + bh4; i++) {
        uint16_t *p_mi_offset = frame_mi_map->p_mi_offset + i *
                                    frame_mi_map->mi_cols_algnsb;
        for (int j = mi_col; j < mi_col + bw4; j++)
            p_mi_offset[j] = offset;
    }
#else
    int32_t blk_mi_row_start = (mi_row - parse_ctx->sb_row_mi) + 1;
    int32_t blk_mi_col_start = (mi_col - parse_ctx->sb_col_mi) + 1;

    /*TODO : Can remove later*/
    assert(blk_mi_row_start > 0);
    assert(blk_mi_col_start > 0);

    /* Update 4x4 nbr map */
    for (int i = blk_mi_row_start; i < blk_mi_row_start + bh4; i++)
        for (int j = blk_mi_col_start; j < blk_mi_col_start + bw4; j++)
            frame_mi_map->cur_sb_mi_map[i][j] = offset;
#endif
}

#if !FRAME_MI_MAP
/* Should be called within same SB */
#endif
/* TODO : Should remove dec_mod_ctxt dependency */
BlockModeInfo* get_cur_mode_info(void *pv_dec_handle,
                              int mi_row, int mi_col, SBInfo *sb_info)
{
    EbDecHandle *dec_handle     = (EbDecHandle *)pv_dec_handle;
    FrameMiMap  *frame_mi_map   = &dec_handle->master_frame_buf.frame_mi_map;

    BlockModeInfo *cur_mi = NULL;
    (void)sb_info;
#if FRAME_MI_MAP
    int32_t cur_sb_row = mi_row >> (frame_mi_map->sb_size_log2 - MI_SIZE_LOG2);
    int32_t cur_sb_col = mi_col >> (frame_mi_map->sb_size_log2 - MI_SIZE_LOG2);
    SBInfo *cur_sb_info = frame_mi_map->pps_sb_info[cur_sb_row *
                                            frame_mi_map->sb_cols + cur_sb_col];
    int32_t offset = *(frame_mi_map->p_mi_offset +
                        mi_row * frame_mi_map->mi_cols_algnsb + mi_col);
    cur_mi = &cur_sb_info->sb_mode_info[offset];
#else
    DecModCtxt  *dec_mod_ctxt = (DecModCtxt *)dec_handle->pv_dec_mod_ctxt;

    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    int32_t cur_blk_mi_row = (mi_row - dec_mod_ctxt->sb_row_mi) + 1;
    int32_t cur_blk_mi_col = (mi_col - dec_mod_ctxt->sb_col_mi) + 1;

    /* Can be removed later */
    assert(cur_blk_mi_row >= 1 && cur_blk_mi_row <= num_mis_in_sb_wd);
    assert(cur_blk_mi_col >= 1 && cur_blk_mi_col <= num_mis_in_sb_wd);

    (void)num_mis_in_sb_wd;
    int32_t offset = frame_mi_map->cur_sb_mi_map[cur_blk_mi_row][cur_blk_mi_col];
    cur_mi = &sb_info->sb_mode_info[offset];
#endif
    return cur_mi;
}

/* TODO : Should remove parse_ctx dependency */
BlockModeInfo * get_left_mode_info(EbDecHandle *dec_handle,
    int mi_row, int mi_col, SBInfo *sb_info)
{
    (void)sb_info;
#if FRAME_MI_MAP
    return get_cur_mode_info(dec_handle, mi_row, mi_col - 1, NULL);
#else
    ParseCtxt   *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameMiMap  *frame_mi_map = &dec_handle->master_frame_buf.frame_mi_map;
    BlockModeInfo  *left_mi = NULL;

    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    int32_t left_blk_mi_row = (mi_row - parse_ctx->sb_row_mi) + 1;
    int32_t left_blk_mi_col = (mi_col - parse_ctx->sb_col_mi) + 1 - 1; // +1 -1 is only for better readibility of logic

    /* Can be removed later */
    assert(left_blk_mi_row >= 1 && left_blk_mi_row <= num_mis_in_sb_wd);
    assert(left_blk_mi_col >= 0 && left_blk_mi_col <= num_mis_in_sb_wd);
    (void)num_mis_in_sb_wd;

    int32_t offset = frame_mi_map->cur_sb_mi_map[left_blk_mi_row][left_blk_mi_col];
    /* From Left SB */
    if (left_blk_mi_col < 1)
        left_mi = &parse_ctx->left_sb_info->sb_mode_info[offset];
    else
        left_mi = &sb_info->sb_mode_info[offset];

    return left_mi;
#endif
}

/* TODO : Should remove parse_ctx dependency */
BlockModeInfo* get_top_mode_info(EbDecHandle *dec_handle,
    int mi_row, int mi_col, SBInfo *sb_info)
{
    (void)sb_info;
#if FRAME_MI_MAP
    return get_cur_mode_info(dec_handle, mi_row - 1, mi_col, NULL);
#else
    ParseCtxt   *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameMiMap  *frame_mi_map = &dec_handle->master_frame_buf.frame_mi_map;
    BlockModeInfo  *top_mi = NULL;

    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    int32_t top_blk_mi_row = (mi_row - parse_ctx->sb_row_mi) + 1 - 1; // +1 -1 is only for better readibility of logic
    int32_t top_blk_mi_col = (mi_col - parse_ctx->sb_col_mi) + 1;

    /* Can be removed later */
    assert(top_blk_mi_row >= 0 && top_blk_mi_row <= num_mis_in_sb_wd);
    assert(top_blk_mi_col >= 1 && top_blk_mi_col <= num_mis_in_sb_wd);
    (void)num_mis_in_sb_wd;

    int32_t offset = frame_mi_map->cur_sb_mi_map[top_blk_mi_row][top_blk_mi_col];
    /* From Top SB */
    if (top_blk_mi_row < 1)
        top_mi = &parse_ctx->above_sb_info->sb_mode_info[offset];
    else
        top_mi = &sb_info->sb_mode_info[offset];

    return top_mi;
#endif
}
