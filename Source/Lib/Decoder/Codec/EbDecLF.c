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

/*SUMMARY
Contains the Decoder Loop Filtering related functions*/

#include "EbDefinitions.h"
#include "EbUtility.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbObuParse.h"
#include "EbDecUtils.h"
#include "EbDeblockingCommon.h"
#include "EbDecNbr.h"
#include "EbDecLF.h"
#include "common_dsp_rtcd.h"
#define FILTER_LEN 4

/*Filter_length is mapped to int indx for the filter tap arrays*/
/*4-> 0   6-> 1   8-> 2   14->3 */
static int8_t filter_map[15] = {-1, -1, -1, -1, 0, -1, 1, -1, 2, -1, -1, -1, -1, -1, 3};

SvtLbdFilterTapFn lbd_vert_filter_tap[FILTER_LEN];
SvtHbdFilterTapFn hbd_vert_filter_tap[FILTER_LEN];
SvtLbdFilterTapFn lbd_horz_filter_tap[FILTER_LEN];
SvtHbdFilterTapFn hbd_horz_filter_tap[FILTER_LEN];

void set_lbd_lf_filter_tap_functions(void) {
    lbd_horz_filter_tap[0] = svt_aom_lpf_horizontal_4;
    lbd_horz_filter_tap[1] = svt_aom_lpf_horizontal_6;
    lbd_horz_filter_tap[2] = svt_aom_lpf_horizontal_8;
    lbd_horz_filter_tap[3] = svt_aom_lpf_horizontal_14;

    lbd_vert_filter_tap[0] = svt_aom_lpf_vertical_4;
    lbd_vert_filter_tap[1] = svt_aom_lpf_vertical_6;
    lbd_vert_filter_tap[2] = svt_aom_lpf_vertical_8;
    lbd_vert_filter_tap[3] = svt_aom_lpf_vertical_14;
}

void set_hbd_lf_filter_tap_functions(void) {
    hbd_horz_filter_tap[0] = svt_aom_highbd_lpf_horizontal_4;
    hbd_horz_filter_tap[1] = svt_aom_highbd_lpf_horizontal_6;
    hbd_horz_filter_tap[2] = svt_aom_highbd_lpf_horizontal_8;
    hbd_horz_filter_tap[3] = svt_aom_highbd_lpf_horizontal_14;

    hbd_vert_filter_tap[0] = svt_aom_highbd_lpf_vertical_4;
    hbd_vert_filter_tap[1] = svt_aom_highbd_lpf_vertical_6;
    hbd_vert_filter_tap[2] = svt_aom_highbd_lpf_vertical_8;
    hbd_vert_filter_tap[3] = svt_aom_highbd_lpf_vertical_14;
}

/*Population of neighbour block lf params for each 4x4 block*/
void fill_4x4_lf_param(LfCtxt* lf_ctxt,
                       int32_t tu_x, int32_t tu_y,
                       int32_t stride,
                       TxSize tx_size,
                       int32_t sub_x, int32_t sub_y,
                       int plane) {
    /*Population of chroma info is done at 4x4
    to be sync with all other luma population*/
    const int txw = tx_size_wide_unit[tx_size] << sub_x;
    const int txh = tx_size_high_unit[tx_size] << sub_y;
    int tx_offset = tu_y * stride + tu_x;
    if (plane == 0) {
        TxSize *tx_size_l = lf_ctxt->tx_size_l + tx_offset;

        *tx_size_l = tx_size;
        for (int b4c = 1; b4c < txw; b4c++)
            *(tx_size_l + b4c) = *tx_size_l;

        for (int b4r = 1; b4r < txh; b4r++) {
            for (int b4c = 0; b4c < txw; b4c++)
                *(tx_size_l + b4r * stride + b4c) = *tx_size_l;
        }
    } else {
        assert(plane == 1);
        TxSize *tx_size_uv = lf_ctxt->tx_size_uv + tx_offset;

        *tx_size_uv = tx_size;
        for (int b4c = 1; b4c < txw; b4c++)
            *(tx_size_uv + b4c) = *tx_size_uv;

        for (int b4r = 1; b4r < txh; b4r++) {
            for (int b4c = 0; b4c < txw; b4c++)
                *(tx_size_uv + b4r * stride + b4c) = *tx_size_uv;
        }
    }
}

/*Function to get transform size*/
static INLINE TxSize dec_get_transform_size(const EdgeDir edge_dir, TxSize tx_size) {
    /*since in case of chrominance or non-square transorm need to convert
    transform size into transform size in particular direction.
    for vertical edge, filter direction is horizontal, for horizontal
    edge, filter direction is vertical.*/
    tx_size = (VERT_EDGE == edge_dir) ? txsize_horz_map[tx_size] : txsize_vert_map[tx_size];
    return tx_size;
}

/*Return TxSize from get_transform_size(), so it is plane and direction
 awared*/
static AOM_FORCE_INLINE TxSize dec_set_lpf_parameters(Av1DeblockingParameters *const params,
                                                      FrameHeader *frm_hdr,
                                                      EbColorConfig *color_config,
                                                      BlockModeInfo *mi_cur,
                                                      BlockModeInfo *mi_near,
                                                      LfCtxt* lf_ctxt,
                                                      const EdgeDir edge_dir,
                                                      const uint32_t x, const uint32_t y,
                                                      const int32_t plane,
                                                      int32_t *sb_delta_lf,
                                                      int32_t *sb_delta_lf_prev,
                                                      int32_t *min_tx_dim) {
    LoopFilterInfoN *lf_info = &lf_ctxt->lf_info;

    /*reset to initial values*/
    params->filter_length = 0;

    /*no deblocking is required*/
    const uint32_t width  = frm_hdr->frame_size.frame_width;
    const uint32_t height = frm_hdr->frame_size.frame_height;
    if ((width <= x) || (height <= y)) {
        /*just return the smallest transform unit size*/
        return TX_4X4;
    }
    const uint32_t sub_x     = (plane > 0) ? color_config->subsampling_x : 0;
    const uint32_t sub_y     = (plane > 0) ? color_config->subsampling_y : 0;
    const int32_t  mi_row    = y >> MI_SIZE_LOG2;
    const int32_t  mi_col    = x >> MI_SIZE_LOG2;
    const int32_t  lf_stride = frm_hdr->mi_stride;

    /*Alligned at 8x8 boundary :
    case 1 : if 8x8 block, last 4x4 block property is populated for all
            the TL, Left and Top 4x4 block
    case 2 : if 8xH & Wx8 block , then 2nd 4xH or Wx4 block property is
            populated for the left or rigth block */
    int32_t lf_offset = (mi_row | sub_y) * lf_stride + (mi_col | sub_x);

    TxSize *lf_tx_l_cur  = lf_ctxt->tx_size_l + lf_offset;
    TxSize *lf_tx_uv_cur = lf_ctxt->tx_size_uv + lf_offset;

    if (lf_tx_l_cur  == NULL) return TX_INVALID;
    if (lf_tx_uv_cur == NULL) return TX_INVALID;

    /*TODO : to be written get_trnasfrom_size function*/
    TxSize ts_cur = plane == 0 ? *lf_tx_l_cur : *lf_tx_uv_cur;
    const TxSize ts = dec_get_transform_size(edge_dir, ts_cur);
    /* For rectangulare blocks,as we need min_hiegt/width only,
        so there is no need of converting to square transfom*/
    *min_tx_dim = (VERT_EDGE == edge_dir) ? tx_size_high_unit[ts_cur] :
                                            tx_size_wide_unit[ts_cur];

    {
        const uint32_t coord = (VERT_EDGE == edge_dir) ? (x >> sub_x) : (y >> sub_y);
        const uint32_t transform_masks = edge_dir == VERT_EDGE ?
                                         tx_size_wide[ts] - 1 :
                                         tx_size_high[ts] - 1;
        const int32_t txb_edge = (coord & transform_masks) ? (0) : (1);

        if (!txb_edge) return ts;
        /*prepare outer edge parameters. deblock the edge if it's an edge of a TU*/
        {
            uint32_t curr_level; // Added to address 4x4 problem
            PredictionMode mode = (mi_cur->mode == INTRA_MODE_4x4)
                                  ? DC_PRED : mi_cur->mode;
            if (frm_hdr->delta_lf_params.delta_lf_present)
                curr_level = get_filter_level_delta_lf(frm_hdr, edge_dir,
                                                       plane,
                                                       sb_delta_lf,
                                                       mi_cur->segment_id,
                                                       mode,
                                                       mi_cur->ref_frame[0]);
            else
                curr_level = lf_info->lvl[plane]
                                         [mi_cur->segment_id]
                                         [edge_dir]
                                         [mi_cur->ref_frame[0]]
                                         [mode_lf_lut[mode]];

            const int32_t curr_skipped =
                mi_cur->skip && is_inter_block_no_intrabc(mi_cur->ref_frame[0]);

            uint32_t level = curr_level;
            if (coord) {
                /* Since block are alligned at 8x8 boundary , so prev block
                    will be cur - 1<<sub_x */
                TxSize *lf_tx_l_prev =  (VERT_EDGE == edge_dir) ?
                                        lf_tx_l_cur - (uint64_t)(1 << sub_x) :
                                        lf_tx_l_cur - (lf_stride << sub_y);
                TxSize *lf_tx_uv_prev = (VERT_EDGE == edge_dir) ?
                                        lf_tx_uv_cur - (uint64_t)(1 << sub_x) :
                                        lf_tx_uv_cur - (lf_stride << sub_y);

                if (lf_tx_l_prev  == NULL) return TX_INVALID;
                if (lf_tx_uv_prev == NULL) return TX_INVALID;

                TxSize pv_ts_act = plane == 0 ? *lf_tx_l_prev : *lf_tx_uv_prev;
                const TxSize pv_ts = dec_get_transform_size(edge_dir, pv_ts_act);

                /* For rectangulare blocks,as we need min_hiegt/width only,
                    so there is no need of converting to square transfom*/
                int32_t prev_ts_dim = (VERT_EDGE == edge_dir) ?
                                       tx_size_high_unit[pv_ts_act] :
                                       tx_size_wide_unit[pv_ts_act];
                *min_tx_dim = AOMMIN(*min_tx_dim, prev_ts_dim);

                uint32_t pv_lvl;
                mode = (mi_near->mode == INTRA_MODE_4x4)
                        ? DC_PRED : mi_near->mode;
                if (frm_hdr->delta_lf_params.delta_lf_present)
                    pv_lvl = get_filter_level_delta_lf(frm_hdr,
                                                       edge_dir,
                                                       plane,
                                                       sb_delta_lf_prev,
                                                       mi_near->segment_id,
                                                       mode,
                                                       mi_near->ref_frame[0]);
                else
                    pv_lvl = lf_info->lvl[plane][mi_near->segment_id]
                                         [edge_dir]
                                         [mi_near->ref_frame[0]]
                                         [mode_lf_lut[mode]];

                const int32_t pv_skip = mi_near->skip &&
                                        is_inter_block_no_intrabc(mi_near->ref_frame[0]);
                const BlockSize bsize = get_plane_block_size(mi_near->sb_type,
                                                             sub_x, sub_y);
                assert(bsize < BlockSizeS_ALL);
                const int32_t prediction_masks = edge_dir == VERT_EDGE
                                                 ? block_size_wide[bsize] - 1
                                                 : block_size_high[bsize] - 1;
                const int32_t pu_edge = !(coord & prediction_masks);
                /*if the current and the previous blocks are skipped,
                deblock the edge if the edge belongs to a PU's edge only.*/
                if ((curr_level || pv_lvl) && (!pv_skip || !curr_skipped || pu_edge)) {
                    const TxSize min_ts = AOMMIN(ts, pv_ts);
                    if (TX_4X4 >= min_ts)
                        params->filter_length = 4;
                    else if (TX_8X8 == min_ts) {
                        if (plane != 0)
                            params->filter_length = 6;
                        else
                            params->filter_length = 8;
                    } else {
                        params->filter_length = 14;
                        /* No wide filtering for chroma plane*/
                        if (plane != 0) params->filter_length = 6;
                    }

                    /* update the level if the current block is skipped,
                    but the previous one is not*/
                    level = (curr_level) ? (curr_level) : (pv_lvl);
                }
            }
            /* prepare common parameters*/
            if (params->filter_length) {
                const LoopFilterThresh *const limits = lf_info->lfthr + level;
                params->lim                          = limits->lim;
                params->mblim                        = limits->mblim;
                params->hev_thr                      = limits->hev_thr;
            }
        }
    }
    return ts;
}

/*It applies Vertical Loop Filtering in a superblock*/
void dec_av1_filter_block_plane_vert(EbDecHandle *dec_handle,
                                     SBInfo *sb_info,
                                     EbPictureBufferDesc *recon_picture_buf,
                                     LfCtxt *lf_ctxt,
                                     const int32_t num_planes,
                                     const int32_t sb_mi_row,
                                     const int32_t sb_mi_col,
                                     int32_t *sb_delta_lf)
{
    FrameHeader *frm_hdr = &dec_handle->frame_header;
    EbColorConfig *color_config = &dec_handle->seq_header.color_config;
    EbBitDepthEnum is16bit = (recon_picture_buf->bit_depth > EB_8BIT ||
        dec_handle->is_16bit_pipeline);
    int32_t sub_x = color_config->subsampling_x;
    int32_t sub_y = color_config->subsampling_y;
    uint8_t no_lf_luma = !(frm_hdr->loop_filter_params.filter_level[0]) &&
                         !(frm_hdr->loop_filter_params.filter_level[1]);
    uint8_t no_lf_u = !(frm_hdr->loop_filter_params.filter_level_u);
    uint8_t no_lf_v = !(frm_hdr->loop_filter_params.filter_level_v);

    void *blk_recon_buf;
    int32_t recon_stride;
    int32_t *sb_delta_lf_left;

    TransformInfo_t *trans_info = NULL;
    uint32_t num_tu;

    uint32_t mi_cols = (&dec_handle->frame_header)->mi_cols;
    uint32_t mi_rows = (&dec_handle->frame_header)->mi_rows;

    BlockModeInfo *mode_info = get_cur_mode_info(dec_handle,
        sb_mi_row, sb_mi_col, sb_info);

    int n_blocks = sb_info->num_block;
    for (int sub_blck = 0; sub_blck < n_blocks; sub_blck++) {

        int32_t blk_mi_row = sb_mi_row + mode_info->mi_row_in_sb;
        int32_t blk_mi_col = sb_mi_col + mode_info->mi_col_in_sb;
        BlockSize bsize = mode_info->sb_type;

        int32_t bw4 = mi_size_wide[bsize];
        int32_t bh4 = mi_size_high[bsize];
        int32_t mb_to_bottom_edge = ((mi_rows - bh4 - blk_mi_row) * MI_SIZE) * 8;
        int32_t mb_to_right_edge  = ((mi_cols - bw4 - blk_mi_col) * MI_SIZE) * 8;

        uint8_t lossless   = frm_hdr->lossless_array[mode_info->segment_id];
        int lossless_block = (lossless && ((bsize >= BLOCK_64X64) &&
                             (bsize <= BLOCK_128X128)));

        int max_blocks_wide = block_size_wide[bsize];
        int max_blocks_high = block_size_high[bsize];
        if (mb_to_right_edge < 0)
            max_blocks_wide += gcc_right_shift(mb_to_right_edge, 3);
        if (mb_to_bottom_edge < 0)
            max_blocks_high += gcc_right_shift(mb_to_right_edge, 3);
        max_blocks_wide = max_blocks_wide >> tx_size_wide_log2[0];
        max_blocks_high = max_blocks_high >> tx_size_high_log2[0];

        int num_chroma_tus = lossless_block ? (max_blocks_wide * max_blocks_high)
            >> (sub_x + sub_y) : mode_info->num_tus[AOM_PLANE_U];

        for (int plane = 0; plane < num_planes; ++plane) {
            if ((plane == 0) && no_lf_luma)
                break;
            else if ((plane == 1) && no_lf_u)
                continue;
            else if ((plane == 2) && no_lf_v)
                continue;

            sub_x = (plane > 0) ? color_config->subsampling_x : 0;
            sub_y = (plane > 0) ? color_config->subsampling_y : 0;

            trans_info = (plane == 2) ? (sb_info->sb_trans_info[plane - 1] +
                          mode_info->first_txb_offset[plane - 1] + num_chroma_tus) :
                          (sb_info->sb_trans_info[plane] +
                          mode_info->first_txb_offset[plane]);

            if (lossless_block)
            {
                assert(trans_info->tx_size == TX_4X4);
                num_tu = (max_blocks_wide * max_blocks_high) >> (sub_x + sub_y);

            }
            else
                num_tu = mode_info->num_tus[!!plane];

            derive_blk_pointers(recon_picture_buf, plane,
                               ((blk_mi_col >> sub_x) * MI_SIZE),
                               ((blk_mi_row >> sub_y) * MI_SIZE),
                               &blk_recon_buf, &recon_stride,
                               sub_x, sub_y);

            BlockModeInfo *mi_left = mode_info;

            Av1DeblockingParameters params;
            memset(&params, 0, sizeof(params));

            for (uint32_t tu = 0; tu < num_tu; tu++)
            {
                void *tu_recon_buf;
                int32_t tu_offset;
                TxSize cur_tx_size;
                int32_t cur_txh, min_txh;

                cur_tx_size = trans_info->tx_size;
                cur_txh = tx_size_high_unit[cur_tx_size];

                tu_offset = (trans_info->txb_y_offset * recon_stride +
                            trans_info->txb_x_offset) << MI_SIZE_LOG2;
                tu_recon_buf = (void*)((uint8_t*)blk_recon_buf
                                + (tu_offset << is16bit));

                uint8_t *p = (uint8_t*)tu_recon_buf;

                /*inner loop always filter vertical edges in a MI block. If MI size
                is 8x8, it will filter the vertical edge aligned with a 8x8 block.
                If 4x4 trasnform is used, it will then filter the internal edge
                aligned with a 4x4 block*/
                /*When there are two 2xX, we combine the them do the LF in the first 2xX edge
                  Since the first 2xX was skipped because of num_tu 0*/

                uint32_t curr_luma_x = (blk_mi_col & (~sub_x)) * MI_SIZE +
                                       ((trans_info->txb_x_offset << sub_x) * MI_SIZE);
                uint32_t curr_luma_y = (blk_mi_row & (~sub_y)) * MI_SIZE +
                                       ((trans_info->txb_y_offset << sub_y)* MI_SIZE);

                int32_t left_mi_row = (blk_mi_row | sub_y) +
                                      (trans_info->txb_y_offset << sub_y);
                int32_t left_mi_col = (blk_mi_col & (~sub_x)) +
                                      (trans_info->txb_x_offset << sub_x);

                /* For VERT_EDGE edge and x_range is for SB scan */
                sb_delta_lf_left = (left_mi_col == sb_mi_col) ?
                                   sb_delta_lf - FRAME_LF_COUNT : sb_delta_lf;

                min_txh = cur_txh;
                for (int32_t temp_txh = 0; temp_txh < cur_txh;
                    temp_txh += min_txh) {

                    if(params.filter_length)
                        memset(&params, 0, sizeof(params));

                    if (left_mi_col > 0) {
                        mi_left = get_left_mode_info(dec_handle,
                                                     left_mi_row,
                                                     left_mi_col,
                                                     sb_info);
                    }

                    cur_tx_size = dec_set_lpf_parameters(&params,
                                                         frm_hdr, color_config,
                                                         mode_info,
                                                         mi_left,
                                                         lf_ctxt,
                                                         VERT_EDGE,
                                                         curr_luma_x,
                                                         curr_luma_y,
                                                         plane,
                                                         sb_delta_lf,
                                                         sb_delta_lf_left,
                                                         &min_txh);

                    if (cur_tx_size == TX_INVALID) {
                        params.filter_length = 0;
                        cur_tx_size = TX_4X4;
                    }

                    /*Do the filtering for only for the actual frame boundry*/
                    int32_t min_high = min_txh << MI_SIZE_LOG2;
                    int32_t frame_height = frm_hdr->frame_size.frame_height;
                    int32_t ext_height = curr_luma_y + (min_high << sub_y);
                    if (frame_height < ext_height)
                        min_high = ((frame_height - (int32_t)curr_luma_y) + sub_y) >> sub_y;

                    for (int32_t h = 0; h < min_high; h += 4) {
                        int8_t filter_idx = filter_map[params.filter_length];
                        if (filter_idx != -1) {
                            if (is16bit)
                                hbd_vert_filter_tap[filter_idx]((uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                                                                recon_stride,
                                                                params.mblim,
                                                                params.lim,
                                                                params.hev_thr,
                                                                recon_picture_buf->bit_depth);
                            else
                                lbd_vert_filter_tap[filter_idx](p,
                                                                recon_stride,
                                                                params.mblim,
                                                                params.lim,
                                                                params.hev_thr);
                        }
                        p += ((4 * recon_stride) << is16bit);
                    }
                    curr_luma_y += (min_txh << (sub_y + MI_SIZE_LOG2));
                    left_mi_row += (min_txh << sub_y);
                }
                trans_info++;
            }
        }
        mode_info++;
    }
}

void dec_av1_filter_block_plane_horz(EbDecHandle *dec_handle, SBInfo *sb_info,
                                     EbPictureBufferDesc *recon_picture_buf,
                                     LfCtxt *lf_ctxt,
                                     const int32_t num_planes,
                                     const int32_t sb_mi_row,
                                     const uint32_t sb_mi_col,
                                     int32_t *sb_delta_lf) {
    FrameHeader *frm_hdr        = &dec_handle->frame_header;
    EbColorConfig *color_config = &dec_handle->seq_header.color_config;

    EbBool is16bit = (recon_picture_buf->bit_depth > EB_8BIT ||
        dec_handle->is_16bit_pipeline);

    int32_t sub_x = color_config->subsampling_x;
    int32_t sub_y = color_config->subsampling_y;

    uint8_t no_lf_luma = !(frm_hdr->loop_filter_params.filter_level[0]) &&
                         !(frm_hdr->loop_filter_params.filter_level[1]);
    uint8_t no_lf_u = !(frm_hdr->loop_filter_params.filter_level_u);
    uint8_t no_lf_v = !(frm_hdr->loop_filter_params.filter_level_v);

    void *blk_recon_buf;
    int32_t recon_stride;
    int32_t *sb_delta_lf_above;

    TransformInfo_t *trans_info = NULL;
    uint32_t num_tu;

    uint32_t mi_cols = (&dec_handle->frame_header)->mi_cols;
    uint32_t mi_rows = (&dec_handle->frame_header)->mi_rows;

    BlockModeInfo *mode_info = get_cur_mode_info(dec_handle, sb_mi_row,
                                                 sb_mi_col, sb_info);
    int n_blocks = sb_info->num_block;
    for (int sub_blck = 0; sub_blck < n_blocks; sub_blck++) {
        int32_t blk_mi_row = sb_mi_row + mode_info->mi_row_in_sb;
        int32_t blk_mi_col = sb_mi_col + mode_info->mi_col_in_sb;
        BlockSize bsize = mode_info->sb_type;

        int32_t bw4 = mi_size_wide[bsize];
        int32_t bh4 = mi_size_high[bsize];
        int32_t mb_to_bottom_edge = ((mi_rows - bh4 - blk_mi_row)
                                    * MI_SIZE) * 8;
        int32_t mb_to_right_edge = ((mi_cols - bw4 - blk_mi_col)
                                    * MI_SIZE) * 8;

        uint8_t lossless = frm_hdr->lossless_array[mode_info->segment_id];
        int lossless_block = (lossless && ((bsize >= BLOCK_64X64) &&
                             (bsize <= BLOCK_128X128)));

        int max_blocks_wide = block_size_wide[bsize];
        int max_blocks_high = block_size_high[bsize];
        if (mb_to_right_edge < 0)
            max_blocks_wide += gcc_right_shift(mb_to_right_edge, 3);
        if (mb_to_bottom_edge < 0)
            max_blocks_high += gcc_right_shift(mb_to_bottom_edge, 3);
        max_blocks_wide = max_blocks_wide >> tx_size_wide_log2[0];
        max_blocks_high = max_blocks_high >> tx_size_high_log2[0];

        int num_chroma_tus = lossless_block ? (max_blocks_wide * max_blocks_high)
            >> (color_config->subsampling_x + color_config->subsampling_y) :
            mode_info->num_tus[AOM_PLANE_U];

        for (int plane = 0; plane < num_planes; ++plane) {
            if ((plane == 0) && no_lf_luma)
                break;
            else if ((plane == 1) && no_lf_u)
                continue;
            else if ((plane == 2) && no_lf_v)
                continue;

            sub_x = (plane > 0) ? color_config->subsampling_x : 0;
            sub_y = (plane > 0) ? color_config->subsampling_y : 0;

            trans_info = (plane == 2) ? (sb_info->sb_trans_info[plane - 1] +
                mode_info->first_txb_offset[plane - 1] + num_chroma_tus) :
                (sb_info->sb_trans_info[plane]
                    + mode_info->first_txb_offset[plane]);

            if (lossless_block)
            {
                assert(trans_info->tx_size == TX_4X4);
                num_tu = (max_blocks_wide * max_blocks_high) >> (sub_x + sub_y);
            }
            else
                num_tu = mode_info->num_tus[!!plane];

            derive_blk_pointers(recon_picture_buf, plane,
                               ((blk_mi_col >> sub_x) * MI_SIZE),
                               ((blk_mi_row >> sub_y) * MI_SIZE),
                               &blk_recon_buf, &recon_stride,
                               sub_x, sub_y);

            BlockModeInfo *mi_top = mode_info;

            Av1DeblockingParameters params;
            memset(&params, 0, sizeof(params));

            for (uint32_t tu = 0; tu < num_tu; tu++)
            {
                void *tu_recon_buf;
                int32_t tu_offset;
                TxSize cur_tx_size;
                int cur_txw, min_txw;

                cur_tx_size = trans_info->tx_size;
                cur_txw = tx_size_wide_unit[cur_tx_size];

                tu_offset = (trans_info->txb_y_offset * recon_stride +
                            trans_info->txb_x_offset) << MI_SIZE_LOG2;
                tu_recon_buf = (void*)((uint8_t*)blk_recon_buf
                               + (tu_offset << is16bit));

                uint8_t *p = (uint8_t*)tu_recon_buf;

                /*inner loop always filter vertical edges in a MI block. If MI size
                is 8x8, it will filter the vertical edge aligned with a 8x8 block.
                If 4x4 trasnform is used, it will then filter the internal edge
                aligned with a 4x4 block*/
                /*When there are two 2xX, we combine the them do the LF in the first 2xX edge
                 Since the first 2xX was skipped because of num_tu 0*/
                uint32_t curr_luma_x = (blk_mi_col & (~sub_x)) * MI_SIZE +
                                       ((trans_info->txb_x_offset << sub_x) * MI_SIZE);
                uint32_t curr_luma_y = (blk_mi_row & (~sub_y)) * MI_SIZE +
                                       ((trans_info->txb_y_offset << sub_y)* MI_SIZE);

                int32_t left_mi_row = (blk_mi_row & (~sub_y)) +
                                      (trans_info->txb_y_offset << sub_y);
                int32_t left_mi_col = (blk_mi_col | sub_x) +
                                      (trans_info->txb_x_offset << sub_x);

                /* For VERT_EDGE edge and x_range is for SB scan */
                sb_delta_lf_above = (left_mi_row == sb_mi_row) ?
                                    sb_delta_lf - lf_ctxt->delta_lf_stride :
                                    sb_delta_lf;

                min_txw = cur_txw;
                for (int temp_txw = 0; temp_txw < cur_txw; temp_txw += min_txw) {

                    if(params.filter_length)
                        memset(&params, 0, sizeof(params));

                    if (left_mi_row > 0)
                        mi_top = get_top_mode_info(dec_handle,
                                                   left_mi_row,
                                                   left_mi_col,
                                                   sb_info);

                    cur_tx_size = dec_set_lpf_parameters(&params,
                                                         frm_hdr, color_config,
                                                         mode_info,
                                                         mi_top,
                                                         lf_ctxt,
                                                         HORZ_EDGE,
                                                         curr_luma_x,
                                                         curr_luma_y,
                                                         plane,
                                                         sb_delta_lf,
                                                         sb_delta_lf_above,
                                                         &min_txw);

                    if (cur_tx_size == TX_INVALID) {
                        params.filter_length = 0;
                        cur_tx_size = TX_4X4;
                    }

                    /*Do the filtering for only for the actual frame boundry*/
                    int32_t min_width = min_txw << MI_SIZE_LOG2;
                    int32_t frame_width = frm_hdr->frame_size.frame_width;
                    int32_t ext_width = curr_luma_x + (min_width << sub_x);
                    if (frame_width < ext_width)
                        min_width = ((frame_width - (int32_t)curr_luma_x) + sub_x) >> sub_x;

                    for (uint8_t w = 0; w < min_width; w += 4) {
                        int filter_idx = filter_map[params.filter_length];

                        if (filter_idx != -1) {
                            if (is16bit)
                                hbd_horz_filter_tap[filter_idx]((uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                                                                recon_stride,
                                                                params.mblim,
                                                                params.lim,
                                                                params.hev_thr,
                                                                recon_picture_buf->bit_depth);
                            else
                                lbd_horz_filter_tap[filter_idx](p,
                                                                recon_stride,
                                                                params.mblim,
                                                                params.lim,
                                                                params.hev_thr);
                        }
                        p += (4 << is16bit);
                    }
                    curr_luma_x += (min_txw << (sub_x + MI_SIZE_LOG2));
                    left_mi_col += (min_txw << sub_x);
                }
                trans_info++;
            }
        }
        mode_info++;
    }
}

/*LF function to filter each SB*/
void dec_loop_filter_sb(EbDecHandle *dec_handle,
                        SBInfo *sb_info,
                        FrameHeader *frm_hdr,
                        SeqHeader *seq_header,
                        EbPictureBufferDesc *recon_picture_buf,
                        LfCtxt *lf_ctxt,
                        const int32_t mi_row, const int32_t mi_col,
                        int32_t plane_start,
                        int32_t plane_end, uint8_t last_col,
                        int32_t *sb_delta_lf) {

    int num_planes = plane_end - plane_start;
    if (frm_hdr->loop_filter_params.combine_vert_horz_lf) {
        /*filter all vertical and horizontal edges in every 64x64 super block
         filter vertical edges*/
        dec_av1_filter_block_plane_vert(dec_handle, sb_info,
                                        recon_picture_buf,
                                        lf_ctxt,
                                        num_planes,
                                        mi_row,
                                        mi_col,
                                        sb_delta_lf);

        /*filter horizontal edges*/
        int32_t max_mib_size =
            seq_header->sb_size == BLOCK_128X128 ? MAX_MIB_SIZE : SB64_MIB_SIZE;

        if ((int32_t)mi_col - max_mib_size >= 0) {
            dec_av1_filter_block_plane_horz(dec_handle,
                                            (sb_info - 1),
                                            recon_picture_buf,
                                            lf_ctxt,
                                            num_planes,
                                            mi_row,
                                            mi_col - max_mib_size,
                                            (sb_delta_lf - FRAME_LF_COUNT));
        }

        /*Filter the horizontal edges of the last sb in each row*/
        if (last_col) {
            dec_av1_filter_block_plane_horz(dec_handle, sb_info,
                                            recon_picture_buf,
                                            lf_ctxt,
                                            num_planes,
                                            mi_row,
                                            mi_col,
                                            sb_delta_lf);
        }
    } else {
        /*filter all vertical edges in every 64x64 super block*/
        dec_av1_filter_block_plane_vert(dec_handle, sb_info,
                                        recon_picture_buf,
                                        lf_ctxt,
                                        num_planes,
                                        mi_row,
                                        mi_col,
                                        sb_delta_lf);

        /*filter all horizontal edges in every 64x64 super block*/
        dec_av1_filter_block_plane_horz(dec_handle, sb_info,
                                        recon_picture_buf,
                                        lf_ctxt,
                                        num_planes,
                                        mi_row,
                                        mi_col,
                                        sb_delta_lf);
    }
}

/* Row level function to trigger loop filter for each superblock*/
void dec_loop_filter_row(EbDecHandle *dec_handle_ptr,
                         EbPictureBufferDesc *recon_picture_buf,
                         LfCtxt *lf_ctxt,
                         uint32_t y_sb_index,
                         int32_t plane_start, int32_t plane_end) {
    MasterFrameBuf *master_frame_buf = &dec_handle_ptr->master_frame_buf;
    CurFrameBuf *   frame_buf        = &master_frame_buf->cur_frame_bufs[0];
    FrameHeader *   frm_hdr          = &dec_handle_ptr->frame_header;
    SeqHeader *     seq_header       = &dec_handle_ptr->seq_header;
    uint8_t         sb_size_log2     = seq_header->sb_size_log2;
    int32_t         sb_size_w        = block_size_wide[seq_header->sb_size];
    int32_t         pic_width_in_sb  = (frm_hdr->frame_size.frame_width + sb_size_w - 1) / sb_size_w;
    uint32_t        sb_origin_y      = y_sb_index << sb_size_log2;

    volatile int32_t *sb_lf_completed_in_prev_row = NULL;
    DecMtlfFrameInfo *lf_frame_info =
        &dec_handle_ptr->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data.lf_frame_info;
    if (y_sb_index) {
        sb_lf_completed_in_prev_row =
            (volatile int32_t *)&lf_frame_info->sb_lf_completed_in_row[y_sb_index - 1];
    }
    int32_t *sb_lf_completed_in_row = &lf_frame_info->sb_lf_completed_in_row[y_sb_index];

    for (int32_t x_sb_index = 0; x_sb_index < pic_width_in_sb; ++x_sb_index) {
        int32_t sb_origin_x     = x_sb_index << sb_size_log2;
        int32_t end_of_row_flag = (x_sb_index == pic_width_in_sb - 1) ? EB_TRUE : EB_FALSE;

        SBInfo *sb_info =
            frame_buf->sb_info + (((y_sb_index * master_frame_buf->sb_cols) + x_sb_index));

        /* Top-Right Sync*/
        if (y_sb_index) {
            while (*sb_lf_completed_in_prev_row < MIN((x_sb_index + 2), pic_width_in_sb - 1))
                ;
        }
        /*LF function for a SB*/
        dec_loop_filter_sb(dec_handle_ptr,
                           sb_info,
                           frm_hdr,
                           seq_header,
                           recon_picture_buf,
                           lf_ctxt,
                           sb_origin_y >> 2,
                           sb_origin_x >> 2,
                           plane_start,
                           plane_end,
                           end_of_row_flag,
                           sb_info->sb_delta_lf);
        /* Update Top-Right Sync*/
        *sb_lf_completed_in_row = x_sb_index;
    }
}

/*Frame level function to trigger loop filter for each superblock*/
void dec_av1_loop_filter_frame(EbDecHandle *dec_handle_ptr,
                               EbPictureBufferDesc *recon_picture_buf,
                               LfCtxt *lf_ctxt,
                               int32_t plane_start,
                               int32_t plane_end,
                               int32_t is_mt, int enable_flag) {
    if (!enable_flag) return;

    FrameHeader *frm_hdr      = &dec_handle_ptr->frame_header;
    SeqHeader *  seq_header   = &dec_handle_ptr->seq_header;
    uint8_t      sb_size_log2 = seq_header->sb_size_log2;

    LoopFilterInfoN *lf_info = &lf_ctxt->lf_info;
    lf_ctxt->delta_lf_stride = dec_handle_ptr->master_frame_buf.sb_cols * FRAME_LF_COUNT;

    int32_t  sb_size_w            = block_size_wide[seq_header->sb_size];
    int32_t  sb_size_h            = block_size_high[seq_header->sb_size];
    uint32_t pic_width_in_sb      = (frm_hdr->frame_size.frame_width + sb_size_w - 1) / sb_size_w;
    uint32_t picture_height_in_sb = (frm_hdr->frame_size.frame_height + sb_size_h - 1) / sb_size_h;

    frm_hdr->loop_filter_params.combine_vert_horz_lf = 1;
    /*init hev threshold const vectors*/
    for (int lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++)
        memset(lf_info->lfthr[lvl].hev_thr, (lvl >> 4), SIMD_WIDTH);

    svt_av1_loop_filter_frame_init(frm_hdr, lf_info, plane_start, plane_end);

    set_lbd_lf_filter_tap_functions();
    set_hbd_lf_filter_tap_functions();

    if (is_mt) {
        for (uint32_t y_sb_index = 0; y_sb_index < picture_height_in_sb; ++y_sb_index) {
            dec_loop_filter_row(dec_handle_ptr,
                                recon_picture_buf,
                                lf_ctxt,
                                y_sb_index,
                                plane_start,
                                plane_end);
        }
    } else {
        /*Loop over a frame : tregger dec_loop_filter_sb for each SB*/
        for (uint32_t y_sb_index = 0; y_sb_index < picture_height_in_sb; ++y_sb_index) {
            for (uint32_t x_sb_index = 0; x_sb_index < pic_width_in_sb; ++x_sb_index) {
                uint32_t sb_origin_x     = x_sb_index << sb_size_log2;
                uint32_t sb_origin_y     = y_sb_index << sb_size_log2;
                EbBool end_of_row_flag = x_sb_index == pic_width_in_sb - 1;

                MasterFrameBuf *master_frame_buf = &dec_handle_ptr->master_frame_buf;
                CurFrameBuf *   frame_buf        = &master_frame_buf->cur_frame_bufs[0];

                SBInfo *sb_info =
                    frame_buf->sb_info + (((y_sb_index * master_frame_buf->sb_cols) + x_sb_index));

                /*sb_info->sb_delta_lf = frame_buf->delta_lf + (FRAME_LF_COUNT *
                    ((y_sb_index * master_frame_buf->sb_cols) + x_sb_index));*/

                /*LF function for a SB*/
                dec_loop_filter_sb(dec_handle_ptr,
                                   sb_info,
                                   frm_hdr,
                                   seq_header,
                                   recon_picture_buf,
                                   lf_ctxt,
                                   sb_origin_y >> 2,
                                   sb_origin_x >> 2,
                                   plane_start,
                                   plane_end,
                                   end_of_row_flag,
                                   sb_info->sb_delta_lf);
            }
        }
    }
}
