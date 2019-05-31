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

/* TODO : Should be moved to nbr file */
void update_nbrs_before_sb(FrameMiMap *frame_mi_map, int32_t sb_col) {
    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    /* Update from top sbrow to current SB top */
    memcpy(&frame_mi_map->cur_sb_mi_map[0][1], &frame_mi_map->top_sbrow_mi_map[sb_col*num_mis_in_sb_wd],
        num_mis_in_sb_wd * sizeof(int16_t));
}

/* TODO : Should be moved to nbr file */
void update_nbrs_after_sb(FrameMiMap *frame_mi_map, int32_t sb_col) {
    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    /* Update cur right 4x4 col as next left 4x4 */
    for (int i = 1; i < num_mis_in_sb_wd + 1; i++)
        frame_mi_map->cur_sb_mi_map[i][0] = frame_mi_map->cur_sb_mi_map[i][num_mis_in_sb_wd];
    /* Update bottom 4x4 of cur as top 4x4 for next row */
    memcpy(&frame_mi_map->top_sbrow_mi_map[sb_col*num_mis_in_sb_wd], &frame_mi_map->cur_sb_mi_map[num_mis_in_sb_wd][1],
        num_mis_in_sb_wd * sizeof(int16_t));
}

/* TODO : Should be moved to nbr file */
void update_block_nbrs(EbDecHandle *dec_handle,
    int mi_row, int mi_col,
    BlockSize subsize)
{
    ParseCtxt   *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameMiMap  *frame_mi_map = &dec_handle->master_frame_buf.frame_mi_map;

    int32_t offset = parse_ctx->cur_mode_info_cnt;
    // int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;
    int bw4 = mi_size_wide[subsize];
    int bh4 = mi_size_high[subsize];

    int32_t blk_mi_row_start = (mi_row - parse_ctx->sb_row_mi) + 1;
    int32_t blk_mi_col_start = (mi_col - parse_ctx->sb_col_mi) + 1;

    /*TODO : Can remove later*/
    assert(blk_mi_row_start > 0);
    assert(blk_mi_col_start > 0);

    /* Update 4x4 nbr map */
    for (int i = blk_mi_row_start; i < blk_mi_row_start + bh4; i++)
        for (int j = blk_mi_col_start; j < blk_mi_col_start + bw4; j++)
            frame_mi_map->cur_sb_mi_map[i][j] = offset;
}

/* TODO : Should be moved to nbr file */
ModeInfo_t* get_cur_mode_info(void *pv_dec_handle,
                              int mi_row, int mi_col, SBInfo *sb_info)
{
    EbDecHandle *dec_handle     = (EbDecHandle *)pv_dec_handle;
    FrameMiMap  *frame_mi_map   = &dec_handle->master_frame_buf.frame_mi_map;
    DecModCtxt *dec_mod_ctxt    = (DecModCtxt *)dec_handle->pv_dec_mod_ctxt;

    ModeInfo_t *cur_mi = NULL;

    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    int32_t cur_blk_mi_row = (mi_row - dec_mod_ctxt->sb_row_mi) + 1;
    int32_t cur_blk_mi_col = (mi_col - dec_mod_ctxt->sb_col_mi) + 1;

    /* Can be removed later */
    assert(cur_blk_mi_row >= 1 && cur_blk_mi_row <= num_mis_in_sb_wd);
    assert(cur_blk_mi_col >= 1 && cur_blk_mi_col <= num_mis_in_sb_wd);

    int32_t offset = frame_mi_map->cur_sb_mi_map[cur_blk_mi_row][cur_blk_mi_col];
    cur_mi = &sb_info->sb_mode_info[offset];

    return cur_mi;
}

/* TODO : Should remove parse_ctx dependency */
ModeInfo_t * get_left_mode_info(EbDecHandle *dec_handle,
    int mi_row, int mi_col, SBInfo *sb_info)
{
    ParseCtxt   *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameMiMap  *frame_mi_map = &dec_handle->master_frame_buf.frame_mi_map;
    ModeInfo_t  *left_mi = NULL;

    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    int32_t left_blk_mi_row = (mi_row - parse_ctx->sb_row_mi) + 1;
    int32_t left_blk_mi_col = (mi_col - parse_ctx->sb_col_mi) + 1 - 1; // +1 -1 is only for better readibility of logic

    /* Can be removed later */
    assert(left_blk_mi_row >= 1 && left_blk_mi_row <= num_mis_in_sb_wd);
    assert(left_blk_mi_col >= 0 && left_blk_mi_col <= num_mis_in_sb_wd);

    int32_t offset = frame_mi_map->cur_sb_mi_map[left_blk_mi_row][left_blk_mi_col];
    /* From Left SB */
    if (left_blk_mi_col < 1)
        left_mi = &parse_ctx->left_sb_info->sb_mode_info[offset];
    else
        left_mi = &sb_info->sb_mode_info[offset];

    return left_mi;
}

/* TODO : Should remove parse_ctx dependency */
ModeInfo_t* get_top_mode_info(EbDecHandle *dec_handle,
    int mi_row, int mi_col, SBInfo *sb_info)
{
    ParseCtxt   *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameMiMap  *frame_mi_map = &dec_handle->master_frame_buf.frame_mi_map;
    ModeInfo_t  *top_mi = NULL;

    int32_t num_mis_in_sb_wd = frame_mi_map->num_mis_in_sb_wd;

    int32_t top_blk_mi_row = (mi_row - parse_ctx->sb_row_mi) + 1 - 1; // +1 -1 is only for better readibility of logic
    int32_t top_blk_mi_col = (mi_col - parse_ctx->sb_col_mi) + 1;

    /* Can be removed later */
    assert(top_blk_mi_row >= 0 && top_blk_mi_row <= num_mis_in_sb_wd);
    assert(top_blk_mi_col >= 1 && top_blk_mi_col <= num_mis_in_sb_wd);

    int32_t offset = frame_mi_map->cur_sb_mi_map[top_blk_mi_row][top_blk_mi_col];
    /* From Top SB */
    if (top_blk_mi_row < 1)
        top_mi = &parse_ctx->above_sb_info->sb_mode_info[offset];
    else
        top_mi = &sb_info->sb_mode_info[offset];

    return top_mi;
}
