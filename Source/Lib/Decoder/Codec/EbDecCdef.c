/*
 * Copyright(c) 2019 Netflix, Inc.
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "EbDecUtils.h"
#include "EbDecNbr.h"
#include "EbUtility.h"
#include "EbDecCdef.h"

/*Compute's whether 8x8 block is skip or not skip block*/
static INLINE int32_t dec_is_8x8_block_skip(BlockModeInfo *mbmi) {
    int32_t is_skip = mbmi->skip;
    /* To get mode info for special cases wx4, 4xh & 4x4 blocksize */
    /* Skip is set to(Skips[r][c] && Skips[r + 1][c] && Skips[r][c + 1] &&
       Skips[r + 1][c + 1]) as per the Spec sec. 7.15.1 */
    if (BLOCK_4X4 == mbmi->sb_type)
        is_skip = mbmi[0].skip && mbmi[1].skip && mbmi[2].skip && mbmi[3].skip;
    else if (1 == mi_size_wide[mbmi->sb_type] || 1 == mi_size_high[mbmi->sb_type]) {
        is_skip = mbmi[0].skip && mbmi[1].skip;
    }
    return is_skip;
}

/*Compute's no. of cdef blocks in units of 8x8 manner in a 64x64 block */
static INLINE int32_t dec_sb_compute_cdef_list(EbDecHandle *dec_handle, SBInfo *sb_info,
                                               FrameHeader *frame_info, int32_t mi_row,
                                               int32_t mi_col, CdefList *dlist, BlockSize bs) {
    int32_t maxc = frame_info->mi_cols - mi_col;
    int32_t maxr = frame_info->mi_rows - mi_row;

    if (bs == BLOCK_128X128 || bs == BLOCK_128X64)
        maxc = AOMMIN(maxc, MI_SIZE_128X128);
    else
        maxc = AOMMIN(maxc, MI_SIZE_64X64);
    if (bs == BLOCK_128X128 || bs == BLOCK_64X128)
        maxr = AOMMIN(maxr, MI_SIZE_128X128);
    else
        maxr = AOMMIN(maxr, MI_SIZE_64X64);

    const int32_t r_step  = mi_size_high[BLOCK_8X8];
    const int32_t c_step  = mi_size_wide[BLOCK_8X8];
    const int32_t r_shift = (r_step == 2);
    const int32_t c_shift = (c_step == 2);

    assert(r_step == 1 || r_step == 2);
    assert(c_step == 1 || c_step == 2);

    int32_t count = 0;
    for (int32_t r = 0; r < maxr; r += r_step) {
        for (int32_t c = 0; c < maxc; c += c_step) {
            BlockModeInfo *mbmi =
                get_cur_mode_info(dec_handle, (mi_row + r), (mi_col + c), sb_info);
            if (!dec_is_8x8_block_skip(mbmi)) {
                dlist[count].by   = (uint8_t)(r >> r_shift);
                dlist[count].bx   = (uint8_t)(c >> c_shift);
                dlist[count].skip = 0;
                count++;
            }
        }
    }
    return count;
}

void svt_cdef_block(EbDecHandle *dec_handle, int32_t *mi_wide_l2, int32_t *mi_high_l2,
                    uint16_t **colbuf, uint8_t *prev_row_cdef, uint8_t *curr_row_cdef, int32_t fbr,
                    int32_t fbc, uint32_t *cdef_left, int32_t num_planes, uint16_t *src,
                    int32_t *curr_recon_stride, uint8_t **curr_blk_recon_buf,
                    uint16_t **linebuf_above, uint16_t **linebuf_curr, int32_t stride) {
    MasterFrameBuf *     master_frame_buf  = &dec_handle->master_frame_buf;
    CurFrameBuf *        frame_buf         = &master_frame_buf->cur_frame_bufs[0];
    EbPictureBufferDesc *recon_picture_ptr = dec_handle->cur_pic_buf[0]->ps_pic_buf;

    int8_t use_highbd = (dec_handle->seq_header.color_config.bit_depth > EB_8BIT ||
        dec_handle->is_16bit_pipeline);
    const int32_t cdef_mask  = 1;
    uint32_t      cdef_count;
    int32_t       coeff_shift = AOMMAX(recon_picture_ptr->bit_depth - 8, 0);
    FrameHeader * frame_info  = &dec_handle->frame_header;
    /*const int32_t stride = (frame_info->mi_cols << MI_SIZE_LOG2) +
        2 * CDEF_HBORDER;*/
    const int32_t nvfb = (frame_info->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    const int32_t nhfb = (frame_info->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    CdefList      dlist[MI_SIZE_64X64 * MI_SIZE_64X64];
    int32_t       dir[CDEF_NBLOCKS][CDEF_NBLOCKS] = {{0}};
    int32_t       var[CDEF_NBLOCKS][CDEF_NBLOCKS] = {{0}};
    /* Logic for getting SBinfo,
    SbInfo points to every super block.*/
    SBInfo *sb_info = NULL;
    if (dec_handle->seq_header.sb_size == BLOCK_128X128) {
        sb_info = frame_buf->sb_info + ((fbr >> 1) * master_frame_buf->sb_cols) + (fbc >> 1);
    } else {
        sb_info = frame_buf->sb_info + ((fbr)*master_frame_buf->sb_cols) + (fbc);
    }

    /*Logic for consuming cdef values from super block,
    Index will vary from 0 to 3 based on position of 64x64 block
    in Superblock.*/
    const int32_t index = dec_handle->seq_header.sb_size == BLOCK_128X128
                              ? (!!(fbc & cdef_mask) + 2 * !!(fbr & cdef_mask))
                              : 0;

    int32_t level, sec_strength;
    int32_t uv_level, uv_sec_strength;
    int32_t nhb, nvb;
    int32_t cstart     = 0;
    curr_row_cdef[fbc] = 0;
    if (sb_info == NULL || sb_info->sb_cdef_strength[index] == -1) {
        *cdef_left = 0;
        return;
    }
    if (!*cdef_left) cstart = -CDEF_HBORDER;
    nhb = AOMMIN(MI_SIZE_64X64, frame_info->mi_cols - MI_SIZE_64X64 * fbc);
    nvb = AOMMIN(MI_SIZE_64X64, frame_info->mi_rows - MI_SIZE_64X64 * fbr);
    int32_t frame_top, frame_left, frame_bottom, frame_right;
    int32_t row_ofset = MI_SIZE_64X64 * fbr;
    int32_t col_ofset = MI_SIZE_64X64 * fbc;

    /*For the current filter block, it's top left corner mi structure (mi_tl)
    is first accessed to check whether the top and left boundaries are
    frame boundaries. Then bottom-left and top-right mi structures are
    accessed to check whether the bottom and right boundaries
    (respectively) are frame boundaries.

    Note that we can't just check the bottom-right mi structure - eg. if
    we're at the right-hand edge of the frame but not the bottom, then
    the bottom-right mi is NULL but the bottom-left is not.  */

    frame_top  = (row_ofset == 0) ? 1 : 0;
    frame_left = (col_ofset == 0) ? 1 : 0;

    if (fbr != nvfb - 1) {
        frame_bottom = ((uint32_t)row_ofset + MI_SIZE_64X64 == frame_info->mi_rows) ? 1 : 0;
    } else
        frame_bottom = 1;

    if (fbc != nhfb - 1) {
        frame_right = ((uint32_t)col_ofset + MI_SIZE_64X64 == frame_info->mi_cols) ? 1 : 0;
    } else
        frame_right = 1;

    const int32_t cdef_strength = sb_info->sb_cdef_strength[index];
    level        = frame_info->cdef_params.cdef_y_strength[cdef_strength] / CDEF_SEC_STRENGTHS;
    sec_strength = frame_info->cdef_params.cdef_y_strength[cdef_strength] % CDEF_SEC_STRENGTHS;
    sec_strength += sec_strength == 3;
    uv_level        = frame_info->cdef_params.cdef_uv_strength[cdef_strength] / CDEF_SEC_STRENGTHS;
    uv_sec_strength = frame_info->cdef_params.cdef_uv_strength[cdef_strength] % CDEF_SEC_STRENGTHS;
    uv_sec_strength += uv_sec_strength == 3;

    if ((level == 0 && sec_strength == 0 && uv_level == 0 && uv_sec_strength == 0) ||
        (cdef_count = dec_sb_compute_cdef_list(dec_handle,
                                               sb_info,
                                               frame_info,
                                               (fbr * MI_SIZE_64X64),
                                               (fbc * MI_SIZE_64X64),
                                               dlist,
                                               BLOCK_64X64)) == 0) {
        *cdef_left = 0;
        return;
    }
    curr_row_cdef[fbc] = 1;
    /*Cdef loop for each plane*/
    for (int32_t pli = 0; pli < num_planes; pli++) {
        uint32_t coffset;
        int32_t  rend, cend;
        int32_t  pri_damping = frame_info->cdef_params.cdef_damping;
        int32_t  sec_damping = pri_damping;
        int32_t  hsize       = nhb << mi_wide_l2[pli];
        int32_t  vsize       = nvb << mi_high_l2[pli];
        int32_t  sub_x       = (pli == 0) ? 0 : dec_handle->seq_header.color_config.subsampling_x;
        int32_t  sub_y       = (pli == 0) ? 0 : dec_handle->seq_header.color_config.subsampling_y;
        if (pli) {
            level        = uv_level;
            sec_strength = uv_sec_strength;
        }

        if (fbc == nhfb - 1)
            cend = hsize;
        else
            cend = hsize + CDEF_HBORDER;

        if (fbr == nvfb - 1)
            rend = vsize;
        else
            rend = vsize + CDEF_VBORDER;

        coffset = fbc * MI_SIZE_64X64 << mi_wide_l2[pli];
        if (fbc == nhfb - 1) {
            /* On the last superblock column, fill in the right border with
               CDEF_VERY_LARGE to avoid filtering with the outside. */
            fill_rect(&src[cend + CDEF_HBORDER],
                      CDEF_BSTRIDE,
                      rend + CDEF_VBORDER,
                      hsize + CDEF_HBORDER - cend,
                      CDEF_VERY_LARGE);
        }
        if (fbr == nvfb - 1) {
            /* On the last superblock row, fill in the bottom border with
               CDEF_VERY_LARGE to avoid filtering with the outside. */
            fill_rect(&src[(rend + CDEF_VBORDER) * CDEF_BSTRIDE],
                      CDEF_BSTRIDE,
                      CDEF_VBORDER,
                      hsize + 2 * CDEF_HBORDER,
                      CDEF_VERY_LARGE);
        }
        uint8_t *rec_buff   = 0;
        uint32_t rec_stride = 0;
        switch (pli) {
        case 0:
            rec_buff   = curr_blk_recon_buf[0];
            rec_stride = curr_recon_stride[0];
            break;
        case 1:
            rec_buff   = curr_blk_recon_buf[1];
            rec_stride = curr_recon_stride[1];
            break;
        case 2:
            rec_buff   = curr_blk_recon_buf[2];
            rec_stride = curr_recon_stride[2];
            break;
        }
        /* Copy in the pixels we need from the current superblock for
           deringing.*/
        if (use_highbd)
            copy_sb16_16(&src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER + cstart],
                         CDEF_BSTRIDE,
                         (uint16_t *)rec_buff /*xd->plane[pli].dst.buf*/,
                         (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr,
                         coffset + cstart,
                         rec_stride /*xd->plane[pli].dst.stride*/,
                         rend,
                         cend - cstart);
        else
            copy_sb8_16(&src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER + cstart],
                        CDEF_BSTRIDE,
                        rec_buff /*xd->plane[pli].dst.buf*/,
                        (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr,
                        coffset + cstart,
                        rec_stride /*xd->plane[pli].dst.stride*/,
                        rend,
                        cend - cstart);

        if (!prev_row_cdef[fbc]) {
            if (use_highbd)
                copy_sb16_16( //cm,
                    &src[CDEF_HBORDER],
                    CDEF_BSTRIDE,
                    (uint16_t *)rec_buff /*xd->plane[pli].dst.buf*/,
                    (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                    coffset,
                    rec_stride /*xd->plane[pli].dst.stride*/,
                    CDEF_VBORDER,
                    hsize);
            else
                copy_sb8_16( //cm,
                    &src[CDEF_HBORDER],
                    CDEF_BSTRIDE,
                    rec_buff /*xd->plane[pli].dst.buf*/,
                    (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                    coffset,
                    rec_stride /*xd->plane[pli].dst.stride*/,
                    CDEF_VBORDER,
                    hsize);
        } else if (fbr > 0) {
            copy_rect(&src[CDEF_HBORDER],
                      CDEF_BSTRIDE,
                      &linebuf_above[pli][coffset],
                      stride,
                      CDEF_VBORDER,
                      hsize);
        } else {
            fill_rect(&src[CDEF_HBORDER], CDEF_BSTRIDE, CDEF_VBORDER, hsize, CDEF_VERY_LARGE);
        }

        if (!prev_row_cdef[fbc - 1]) {
            if (use_highbd)
                copy_sb16_16( //cm,
                    src,
                    CDEF_BSTRIDE,
                    (uint16_t *)rec_buff /*xd->plane[pli].dst.buf*/,
                    (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                    coffset - CDEF_HBORDER,
                    rec_stride /*xd->plane[pli].
                    dst.stride*/
                    ,
                    CDEF_VBORDER,
                    CDEF_HBORDER);
            else
                copy_sb8_16( //cm,
                    src,
                    CDEF_BSTRIDE,
                    rec_buff /*xd->plane[pli].dst.buf*/,
                    (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                    coffset - CDEF_HBORDER,
                    rec_stride /*xd->plane[pli].
                    dst.stride*/
                    ,
                    CDEF_VBORDER,
                    CDEF_HBORDER);
        } else if (fbr > 0 && fbc > 0) {
            copy_rect(src,
                      CDEF_BSTRIDE,
                      &linebuf_above[pli][coffset - CDEF_HBORDER],
                      stride,
                      CDEF_VBORDER,
                      CDEF_HBORDER);
        } else {
            fill_rect(src, CDEF_BSTRIDE, CDEF_VBORDER, CDEF_HBORDER, CDEF_VERY_LARGE);
        }

        if (!prev_row_cdef[fbc + 1]) {
            if (use_highbd)
                copy_sb16_16( //cm,
                    &src[CDEF_HBORDER + (nhb << mi_wide_l2[pli])],
                    CDEF_BSTRIDE,
                    (uint16_t *)rec_buff /*xd->plane[pli].dst.buf*/,
                    (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                    coffset + hsize,
                    rec_stride /*xd->plane[pli].dst.stride*/,
                    CDEF_VBORDER,
                    CDEF_HBORDER);
            else
                copy_sb8_16( //cm,
                    &src[CDEF_HBORDER + (nhb << mi_wide_l2[pli])],
                    CDEF_BSTRIDE,
                    rec_buff /*xd->plane[pli].dst.buf*/,
                    (MI_SIZE_64X64 << mi_high_l2[pli]) * fbr - CDEF_VBORDER,
                    coffset + hsize,
                    rec_stride /*xd->plane[pli].dst.stride*/,
                    CDEF_VBORDER,
                    CDEF_HBORDER);
        } else if (fbr > 0 && fbc < nhfb - 1) {
            copy_rect(&src[hsize + CDEF_HBORDER],
                      CDEF_BSTRIDE,
                      &linebuf_above[pli][coffset + hsize],
                      stride,
                      CDEF_VBORDER,
                      CDEF_HBORDER);
        } else {
            fill_rect(&src[hsize + CDEF_HBORDER],
                      CDEF_BSTRIDE,
                      CDEF_VBORDER,
                      CDEF_HBORDER,
                      CDEF_VERY_LARGE);
        }

        if (*cdef_left) {
            /* If we deringed the superblock on the left
               then we need to copy in saved pixels. */
            copy_rect(
                src, CDEF_BSTRIDE, colbuf[pli], CDEF_HBORDER, rend + CDEF_VBORDER, CDEF_HBORDER);
        }

        /* Saving pixels in case we need to dering the superblock
            on the right. */
        if (fbc < nhfb - 1)
            copy_rect(colbuf[pli],
                      CDEF_HBORDER,
                      src + hsize,
                      CDEF_BSTRIDE,
                      rend + CDEF_VBORDER,
                      CDEF_HBORDER);

        if (fbr < nvfb - 1) {
            if (use_highbd)
                copy_sb16_16(&linebuf_curr[pli][coffset],
                             stride,
                             (uint16_t *)rec_buff,
                             (MI_SIZE_64X64 << mi_high_l2[pli]) * (fbr + 1) - CDEF_VBORDER,
                             coffset,
                             rec_stride,
                             CDEF_VBORDER,
                             hsize);
            else
                copy_sb8_16(&linebuf_curr[pli][coffset],
                            stride,
                            rec_buff,
                            (MI_SIZE_64X64 << mi_high_l2[pli]) * (fbr + 1) - CDEF_VBORDER,
                            coffset,
                            rec_stride,
                            CDEF_VBORDER,
                            hsize);
        }

        if (frame_top) {
            fill_rect(src, CDEF_BSTRIDE, CDEF_VBORDER, hsize + 2 * CDEF_HBORDER, CDEF_VERY_LARGE);
        }
        if (frame_left) {
            fill_rect(src, CDEF_BSTRIDE, vsize + 2 * CDEF_VBORDER, CDEF_HBORDER, CDEF_VERY_LARGE);
        }
        if (frame_bottom) {
            fill_rect(&src[(vsize + CDEF_VBORDER) * CDEF_BSTRIDE],
                      CDEF_BSTRIDE,
                      CDEF_VBORDER,
                      hsize + 2 * CDEF_HBORDER,
                      CDEF_VERY_LARGE);
        }
        if (frame_right) {
            fill_rect(&src[hsize + CDEF_HBORDER],
                      CDEF_BSTRIDE,
                      vsize + 2 * CDEF_VBORDER,
                      CDEF_HBORDER,
                      CDEF_VERY_LARGE);
        }
        /*Cdef filter calling function for high bit depth */
        if (use_highbd) {
            uint16_t *tmp_buff = (uint16_t *)rec_buff;
            svt_cdef_filter_fb(NULL,
                               &tmp_buff[rec_stride * (MI_SIZE_64X64 * fbr << mi_high_l2[pli]) +
                               (fbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
                               rec_stride,
                               &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER],
                               sub_x,
                               sub_y,
                               dir,
                               NULL,
                               var,
                               pli,
                               dlist,
                               cdef_count,
                               level,
                               sec_strength,
                               pri_damping,
                               sec_damping,
                               coeff_shift);
        }
        /*Cdef filter calling function for 8 bit depth */
        else
            svt_cdef_filter_fb(&rec_buff[rec_stride * (MI_SIZE_64X64 * fbr << mi_high_l2[pli]) +
                               (fbc * MI_SIZE_64X64 << mi_wide_l2[pli])],
                               NULL,
                               rec_stride,
                               &src[CDEF_VBORDER * CDEF_BSTRIDE + CDEF_HBORDER],
                               sub_x,
                               sub_y,
                               dir,
                               NULL,
                               var,
                               pli,
                               dlist,
                               cdef_count,
                               level,
                               sec_strength,
                               pri_damping,
                               sec_damping,
                               coeff_shift);
    } /*cdef plane loop ending*/
    //CHKN filtered data is written back directy to recFrame.
    *cdef_left = 1;
}

void svt_cdef_sb_row_mt(EbDecHandle *dec_handle, int32_t *mi_wide_l2, int32_t *mi_high_l2,
                        uint16_t **colbuf, int32_t sb_fbr, uint16_t *src,
                        int32_t *curr_recon_stride, uint8_t **curr_blk_recon_buf) {
    FrameHeader * frame_info = &dec_handle->frame_header;
    const int32_t nvfb       = (frame_info->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    const int32_t nhfb       = (frame_info->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    /*For SB SIZE 128x128, we need two colbuf because, because we do cdef for
     each 64x64 in SB block in raster scan order,
     i.e for transversing across 0 - 3 64x64s in SB block*/
    uint16_t *colbuf_64[2][3];
    int32_t   sb_size_w     = block_size_wide[dec_handle->seq_header.sb_size];
    int32_t pic_width_in_sb = (dec_handle->seq_header.max_frame_width + sb_size_w - 1) / sb_size_w;

    const int32_t num_planes = av1_num_planes(&dec_handle->seq_header.color_config);

    EbBool sb_128 = dec_handle->seq_header.sb_size == BLOCK_128X128;
    for (int32_t pli = 0; pli < num_planes; pli++) {
        const int32_t block_height = (MI_SIZE_64X64 << mi_high_l2[pli]) + 2 * CDEF_VBORDER;
        /*Filling the colbuff's with some values.*/
        fill_rect(colbuf[pli], CDEF_HBORDER, block_height, CDEF_HBORDER, CDEF_VERY_LARGE);
        colbuf_64[0][pli] = colbuf[pli];
        /*For SB SIZE 128x128, we are allocating extra colbuff for
        transversing across 0 - 3 64x64s in SB block */
        if (sb_128) {
            fill_rect(colbuf[pli + 3], CDEF_HBORDER, block_height, CDEF_HBORDER, CDEF_VERY_LARGE);
            colbuf_64[1][pli] = colbuf[3 + pli];
        }
    }
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    volatile uint32_t *cdef_completed_in_prev_row = NULL;
    uint32_t *         cdef_completed_in_row, nsync = 1;

    uint8_t *curr_row_cdef_map = dec_mt_frame_data->row_cdef_map;
    uint32_t cdef_map_stride   = dec_mt_frame_data->cdef_map_stride;

    if (sb_fbr) {
        cdef_completed_in_prev_row =
            (volatile uint32_t *)&dec_mt_frame_data->cdef_completed_in_row[sb_fbr - 1];
    }
    cdef_completed_in_row = &dec_mt_frame_data->cdef_completed_in_row[sb_fbr];

    int32_t fbc_64, fbr_64;
    int32_t cnt_64           = sb_128 ? 4 : 1;
    int32_t sb_64_shift      = sb_128 ? 1 : 0;
    curr_row_cdef_map        = curr_row_cdef_map + 1 + cdef_map_stride;
    uint32_t cdef_left_64[2] = {1, 1};
    uint8_t *curr_row_cdef[2], *prev_row_cdef[2];

    curr_row_cdef[0] = curr_row_cdef_map + ((sb_fbr << sb_64_shift) * cdef_map_stride);
    prev_row_cdef[0] = curr_row_cdef[0] - cdef_map_stride;
    if (sb_128) {
        curr_row_cdef[1] = curr_row_cdef[0] + cdef_map_stride;
        prev_row_cdef[1] = curr_row_cdef[0];
    }
    /*for loop of SB cols*/
    for (int32_t sb_fbc = 0; sb_fbc < pic_width_in_sb; sb_fbc++) {
        /* Top-Right Sync*/
        if (sb_fbr) {
            if (sb_fbc == pic_width_in_sb - 1) nsync = 0;
            while (*cdef_completed_in_prev_row < (sb_fbc + nsync))
                ;
            //Sleep(5); /* ToDo : Change */
        }
        /*Curr multi thread implementation of cdef goes through every SB SIZE row*/
        /*If SB SIZE is 128x128, as cdef excepts top right sync,
        to process Bottom Right 64x64 block in (n) th SB BLOCK 128x128,
        we must have to process Top Left 64x64 block in (n+1) th  SB BLOCK 128x128 in the current SB ROW,
        below logic helps in do so*/
        if (0 == sb_fbc && sb_128) {
            int32_t row_cnt = 0;
            fbr_64          = (sb_fbr << sb_64_shift) + row_cnt;
            fbc_64          = row_cnt;
            svt_cdef_block(
                dec_handle,
                mi_wide_l2,
                mi_high_l2,
                colbuf_64[row_cnt],
                prev_row_cdef[row_cnt],
                curr_row_cdef[row_cnt],
                fbr_64,
                fbc_64,
                &cdef_left_64[row_cnt],
                num_planes,
                src,
                curr_recon_stride,
                curr_blk_recon_buf,
                dec_mt_frame_data->cdef_linebuf[fbr_64], /*above*/
                dec_mt_frame_data->cdef_linebuf[AOMMIN(fbr_64 + 1, nvfb - 1)], /*current*/
                dec_mt_frame_data->cdef_linebuf_stride);
        }

        for (int32_t cnt = 0; cnt < cnt_64; cnt++) {
            int32_t row_cnt = (cnt & 2) >> 1;
            /*transversing across 0 - 3 64x64s or 1 128x128*/
            fbr_64 = (sb_fbr << sb_64_shift) + (row_cnt);
            fbc_64 = (sb_fbc << sb_64_shift) + (cnt & 1) + (sb_128 ? !row_cnt : 0);
            if (fbc_64 >= nhfb) continue;
            if (fbr_64 >= nvfb) continue;

            svt_cdef_block(
                dec_handle,
                mi_wide_l2,
                mi_high_l2,
                colbuf_64[row_cnt],
                prev_row_cdef[row_cnt],
                curr_row_cdef[row_cnt],
                fbr_64,
                fbc_64,
                &cdef_left_64[row_cnt],
                num_planes,
                src,
                curr_recon_stride,
                curr_blk_recon_buf,
                dec_mt_frame_data->cdef_linebuf[fbr_64], /*above*/
                dec_mt_frame_data->cdef_linebuf[AOMMIN(fbr_64 + 1, nvfb - 1)], /*current*/
                dec_mt_frame_data->cdef_linebuf_stride);
        }
        /* Update Top-Right Sync*/
        *cdef_completed_in_row = sb_fbc;
    }
}

/* Frame level call, for CDEF */
void svt_cdef_frame(EbDecHandle *dec_handle, int enable_flag) {
    if (!enable_flag) return;

    EbPictureBufferDesc *recon_picture_ptr = dec_handle->cur_pic_buf[0]->ps_pic_buf;

    uint8_t *     curr_blk_recon_buf[MAX_MB_PLANE];
    int32_t       curr_recon_stride[MAX_MB_PLANE];
    FrameHeader * frame_info = &dec_handle->frame_header;
    const int32_t num_planes = av1_num_planes(&dec_handle->seq_header.color_config);

    DECLARE_ALIGNED(16, uint16_t, src[CDEF_INBUF_SIZE]);
    uint16_t *    linebuf[3];
    uint16_t *    colbuf[3];
    uint8_t *     row_cdef, *prev_row_cdef, *curr_row_cdef;
    int32_t       mi_wide_l2[3];
    int32_t       mi_high_l2[3];
    const int32_t nvfb        = (frame_info->mi_rows + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    const int32_t nhfb        = (frame_info->mi_cols + MI_SIZE_64X64 - 1) / MI_SIZE_64X64;
    row_cdef                  = (uint8_t *)svt_aom_malloc(sizeof(*row_cdef) * (nhfb + 2) * 2);

    assert(row_cdef != NULL);
    memset(row_cdef, 1, sizeof(*row_cdef) * (nhfb + 2) * 2);
    prev_row_cdef = row_cdef + 1;
    curr_row_cdef = prev_row_cdef + nhfb + 2;

    const int32_t stride = (frame_info->mi_cols << MI_SIZE_LOG2) + 2 * CDEF_HBORDER;

    for (int32_t pli = 0; pli < num_planes; pli++) {
        int32_t sub_x   = (pli == 0) ? 0 : dec_handle->seq_header.color_config.subsampling_x;
        int32_t sub_y   = (pli == 0) ? 0 : dec_handle->seq_header.color_config.subsampling_y;

        mi_wide_l2[pli] = MI_SIZE_LOG2 - sub_x;
        mi_high_l2[pli] = MI_SIZE_LOG2 - sub_y;

        /*Deriveing  recon pict buffer ptr's*/
        derive_blk_pointers(recon_picture_ptr,
                            pli,
                            0,
                            0,
                            (void *)&curr_blk_recon_buf[pli],
                            &curr_recon_stride[pli],
                            sub_x,
                            sub_y);
        /*Allocating memory for line buffes->to fill from src if needed*/
        linebuf[pli] = (uint16_t *)svt_aom_malloc(sizeof(*linebuf) * CDEF_VBORDER * stride);
        /*Allocating memory for col buffes->to fill from src if needed*/
        colbuf[pli] = (uint16_t *)svt_aom_malloc(
            sizeof(*colbuf) * ((CDEF_BLOCKSIZE << mi_high_l2[pli]) + 2 * CDEF_VBORDER) *
            CDEF_HBORDER);
    }

    /*Loop for 64x64 block wise, along col wise for frame size*/
    for (int32_t fbr = 0; fbr < nvfb; fbr++) {
        for (int32_t pli = 0; pli < num_planes; pli++) {
            const int32_t block_height = (MI_SIZE_64X64 << mi_high_l2[pli]) + 2 * CDEF_VBORDER;
            /*Filling the colbuff's with some values.*/
            fill_rect(colbuf[pli], CDEF_HBORDER, block_height, CDEF_HBORDER, CDEF_VERY_LARGE);
        }

        uint32_t cdef_left = 1;
        /*Loop for 64x64 block wise, along row wise for frame size*/
        for (int32_t fbc = 0; fbc < nhfb; fbc++) {
            svt_cdef_block(dec_handle,
                           mi_wide_l2,
                           mi_high_l2,
                           colbuf,
                           prev_row_cdef,
                           curr_row_cdef,
                           fbr,
                           fbc,
                           &cdef_left,
                           num_planes,
                           src,
                           curr_recon_stride,
                           curr_blk_recon_buf,
                           linebuf,
                           linebuf,
                           stride);
        }
        uint8_t *tmp  = prev_row_cdef;
        prev_row_cdef = curr_row_cdef;
        curr_row_cdef = tmp;
    }
    svt_aom_free(row_cdef);
    for (int32_t pli = 0; pli < num_planes; pli++) {
        svt_aom_free(linebuf[pli]);
        svt_aom_free(colbuf[pli]);
    }
}
