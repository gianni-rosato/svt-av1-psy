/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

/*SUMMARY
Contains the Decoder Loop Filtering related functions*/

#include "EbDefinitions.h"
#include "EbUtility.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbObuParse.h"
#include "EbDecProcessFrame.h"
#include "EbDecUtils.h"
#include "EbDeblockingFilter.h"
#include "EbDecLF.h"

/*TO check if block is INTER or not*/
static INLINE int32_t dec_is_inter_block(const LFBlockParamL *lf_block_l) {
    return /*is_intrabc_block(mbmi) ||*/ lf_block_l->ref_frame_0 > INTRA_FRAME;
}

/*Population of neighbour block LUMA params for each 4x4 block*/
void fill_4x4_param_luma(LFBlockParamL* lf_block_l,
    int32_t tu_x, int32_t tu_y, int32_t stride,
    TxSize tx_size, ModeInfo_t *mode_info)
{
    const int txw = tx_size_wide_unit[tx_size];
    const int txh = tx_size_high_unit[tx_size];

    lf_block_l += tu_y * stride + tu_x;

    for (int b4r = 0; b4r < txh; b4r++) {
        for (int b4c = 0; b4c < txw; b4c++) {
            if (b4r == 0 && b4c == 0) {
                lf_block_l->skip = mode_info->skip;
                lf_block_l->segment_id = mode_info->segment_id;
                lf_block_l->bsize = mode_info->sb_type;
                lf_block_l->tx_size_l = tx_size;
                lf_block_l->ref_frame_0 = mode_info->ref_frame[0];
                lf_block_l->mode = mode_info->mode;
            }
            else {
                memcpy(lf_block_l + b4r * stride + b4c,
                    lf_block_l, sizeof(LFBlockParamL));
            }
        }
    }
}

/*Population of neighbour block UV params for each 4x4 block*/
void fill_4x4_param_uv(LFBlockParamUV* lf_block_uv, int32_t tu_x, int32_t tu_y,
        int32_t stride, TxSize tx_size, int32_t sub_x, int32_t sub_y)
{
    /*Population of chroma info is done at 4x4
    to be sync with all other luma population*/
    const int txw = tx_size_wide_unit[tx_size] << sub_x;
    const int txh = tx_size_high_unit[tx_size] << sub_y;

    lf_block_uv += tu_y * stride + tu_x;

    for (int b4r = 0; b4r < txh; b4r++) {
        for (int b4c = 0; b4c < txw; b4c++) {

            if (b4r == 0 && b4c == 0) {
                lf_block_uv->tx_size_uv = tx_size;
            }
            else {
                memcpy(lf_block_uv + b4r * stride + b4c,
                    lf_block_uv, sizeof(LFBlockParamUV));
            }
        }
    }
}

/*Function for calculating the levels for each block*/
static uint8_t dec_get_filter_level(FrameHeader*frm_hdr,
    const LoopFilterInfoN *lfi_n,const int32_t dir_idx, int32_t plane,
    LFBlockParamL* lf_block_l)
{
    const int32_t segment_id = lf_block_l->segment_id;
    /*Added to address 4x4 problem*/
    PredictionMode mode;
    mode = (lf_block_l->mode == INTRA_MODE_4x4) ? DC_PRED : lf_block_l->mode;

    if (frm_hdr->delta_lf_params.delta_lf_present) {
        assert(0); int32_t delta_lf = -1;
        if (frm_hdr->delta_lf_params.delta_lf_multi) {
            /*const int32_t delta_lf_idx = delta_lf_id_lut[plane][dir_idx];*/
            assert(0);//delta_lf = mbmi->curr_delta_lf[delta_lf_idx];
        }
        else {
            assert(0); //delta_lf = mbmi->current_delta_lf_from_base;
        }
        int32_t base_level;
        if (plane == 0)
            base_level = frm_hdr->loop_filter_params.filter_level[dir_idx];
        else if (plane == 1)
            base_level = frm_hdr->loop_filter_params.filter_level_u;
        else
            base_level = frm_hdr->loop_filter_params.filter_level_v;
        int32_t lvl_seg = clamp(delta_lf + base_level, 0, MAX_LOOP_FILTER);
        assert(plane >= 0 && plane <= 2);
        /*const int32_t seg_lf_feature_id = seg_lvl_lf_lut[plane][dir_idx];*/
        assert(0);/*if (segfeature_active(&cm->seg, segment_id, seg_lf_feature_id)) {
            const int32_t data =
                get_segdata(&cm->seg, segment_id, seg_lf_feature_id);
            lvl_seg = clamp(lvl_seg + data, 0, MAX_LOOP_FILTER);
        }*/

        if (frm_hdr->loop_filter_params.mode_ref_delta_enabled) {
            const int32_t scale = 1 << (lvl_seg >> 5);
            lvl_seg += frm_hdr->loop_filter_params.ref_deltas
                [lf_block_l->ref_frame_0] * scale;
            if (lf_block_l->ref_frame_0 > INTRA_FRAME)
                lvl_seg += frm_hdr->loop_filter_params.mode_deltas
                    [mode_lf_lut[mode]] * scale;
            lvl_seg = clamp(lvl_seg, 0, MAX_LOOP_FILTER);
        }
        return lvl_seg;
    }
    else {
        ASSERT(mode < MB_MODE_COUNT);
        return lfi_n->lvl[plane][segment_id][dir_idx][lf_block_l->ref_frame_0]
            [mode_lf_lut[mode]];
    }
}

/*Function to get transform size*/
static TxSize dec_get_transform_size(const EDGE_DIR edge_dir, TxSize tx_size) {
    /*since in case of chrominance or non-square transorm need to convert
    transform size into transform size in particular direction.
    for vertical edge, filter direction is horizontal, for horizontal
    edge, filter direction is vertical.*/
    tx_size = (VERT_EDGE == edge_dir) ? txsize_horz_map[tx_size]
        : txsize_vert_map[tx_size];
    return tx_size;
}

/*Return TxSize from get_transform_size(), so it is plane and direction
 awared*/
static TxSize dec_set_lpf_parameters(AV1_DEBLOCKING_PARAMETERS *const params,
    EbPictureBufferDesc *recon_picture_buf, FrameHeader *frm_hdr,
    EbColorConfig *color_config, LFCtxt *lf_ctxt, LoopFilterInfoN *lf_info,
    const EDGE_DIR edge_dir, const uint32_t x,
    const uint32_t y, const int32_t plane)
{
    /*reset to initial values*/
    params->filter_length = 0;

    /*no deblocking is required*/
    const uint32_t width = recon_picture_buf->width;
    const uint32_t height = recon_picture_buf->height;
    if ((width <= x) || (height <= y)) {
        /*just return the smallest transform unit size*/
        return TX_4X4;
    }
    const uint32_t sub_x = (plane > 0) ? color_config->subsampling_x : 0;
    const uint32_t sub_y = (plane > 0) ? color_config->subsampling_y : 0;
    const int32_t mi_row = y >> MI_SIZE_LOG2;
    const int32_t mi_col = x >> MI_SIZE_LOG2;
    const int32_t lf_stride = frm_hdr->mi_stride;

    /*Alligned at 8x8 boundary :
    case 1 : if 8x8 block, last 4x4 block property is populated for all
            the TL, Left and Top 4x4 block
    case 2 : if 8xH & Wx8 block , then 2nd 4xH or Wx4 block property is
            populated for the left or rigth block */
    int32_t lf_offset = (mi_row | sub_y) * lf_stride + (mi_col | sub_x);

    LFBlockParamL* lf_block_l_cur = lf_ctxt->lf_block_luma + lf_offset;
    LFBlockParamUV* lf_block_uv_cur = lf_ctxt->lf_block_uv + lf_offset;

    if (lf_block_l_cur == NULL) return TX_INVALID;
    if (lf_block_uv_cur == NULL) return TX_INVALID;

    /*TODO : to be written get_trnasfrom_size function*/
    const TxSize ts = plane == 0 ?
        dec_get_transform_size(edge_dir, lf_block_l_cur->tx_size_l) :
        dec_get_transform_size(edge_dir, lf_block_uv_cur->tx_size_uv);
    {
        const uint32_t coord = (VERT_EDGE == edge_dir) ?
            (x >> sub_x) : (y >> sub_y);
        const uint32_t transform_masks =
            edge_dir == VERT_EDGE ? tx_size_wide[ts] - 1 : tx_size_high[ts] - 1;
        const int32_t tu_edge = (coord & transform_masks) ? (0) : (1);

        if (!tu_edge) return ts;
        /*prepare outer edge parameters. deblock the edge if it's an edge of a TU*/
        {
            const uint32_t curr_level = dec_get_filter_level(frm_hdr, lf_info,
                edge_dir, plane, lf_block_l_cur);
            const int32_t curr_skipped = lf_block_l_cur->skip
                && dec_is_inter_block(lf_block_l_cur);
            uint32_t level = curr_level;
            if (coord) {
                /* Since block are alligned at 8x8 boundary , so prev block
                    will be cur - 1<<sub_x */
                LFBlockParamL* lf_block_l_prev = (VERT_EDGE == edge_dir) ?
                    lf_block_l_cur - (uint64_t)(1 << sub_x) :
                    lf_block_l_cur - (lf_stride << sub_y);
                LFBlockParamUV* lf_block_uv_prev = (VERT_EDGE == edge_dir) ?
                    lf_block_uv_cur - (uint64_t)(1 << sub_x) :
                    lf_block_uv_cur - (lf_stride << sub_y);

                if (lf_block_l_prev == NULL) return TX_INVALID;
                if (lf_block_uv_prev == NULL) return TX_INVALID;

                const TxSize pv_ts = plane == 0 ?
                    dec_get_transform_size(edge_dir,
                        lf_block_l_prev->tx_size_l) :
                    dec_get_transform_size(edge_dir,
                        lf_block_uv_prev->tx_size_uv);

                const uint32_t pv_lvl =
                    dec_get_filter_level(frm_hdr, lf_info,
                        edge_dir, plane, lf_block_l_prev);
                const int32_t pv_skip = lf_block_l_prev->skip &&
                    dec_is_inter_block(lf_block_l_prev);
                const BlockSize bsize =
                    get_plane_block_size(lf_block_l_prev->bsize,
                        sub_x, sub_y);
                assert(bsize < BlockSizeS_ALL);
                const int32_t prediction_masks = edge_dir == VERT_EDGE
                    ? block_size_wide[bsize] - 1
                    : block_size_high[bsize] - 1;
                const int32_t pu_edge = !(coord & prediction_masks);
                /*if the current and the previous blocks are skipped,
                deblock the edge if the edge belongs to a PU's edge only.*/
                if ((curr_level || pv_lvl) &&
                    (!pv_skip || !curr_skipped || pu_edge)) {
                    const TxSize min_ts = AOMMIN(ts, pv_ts);
                    if (TX_4X4 >= min_ts)
                        params->filter_length = 4;
                    else if (TX_8X8 == min_ts) {
                        if (plane != 0)
                            params->filter_length = 6;
                        else
                            params->filter_length = 8;
                    }
                    else {
                        params->filter_length = 14;
                        /* No wide filtering for chroma plane*/
                        if (plane != 0)
                            params->filter_length = 6;
                    }

                    /* update the level if the current block is skipped,
                    but the previous one is not*/
                    level = (curr_level) ? (curr_level) : (pv_lvl);
                }
            }
            /* prepare common parameters*/
            if (params->filter_length) {
                const LoopFilterThresh *const limits = lf_info->lfthr + level;
                params->lim = limits->lim;
                params->mblim = limits->mblim;
                params->hev_thr = limits->hev_thr;
            }
        }
    }
    return ts;
}

/*It applies Vertical Loop Filtering in a superblock*/
void dec_av1_filter_block_plane_vert(
    FrameHeader *frm_hdr,EbColorConfig *color_config,
    EbPictureBufferDesc *recon_picture_buf, LFCtxt *lf_ctxt,
    LoopFilterInfoN *lf_info, const int32_t plane, BlockSize sb_size,
    const uint32_t mi_row, const uint32_t mi_col)
{
    EbBitDepthEnum is16bit = recon_picture_buf->bit_depth > 8;
    const int32_t row_step = MI_SIZE >> MI_SIZE_LOG2;
    const uint32_t sub_x = (plane > 0) ? color_config->subsampling_x : 0;
    const uint32_t sub_y = (plane > 0) ? color_config->subsampling_y : 0;
    const int32_t y_range = sb_size == BLOCK_128X128
        ? (MAX_MIB_SIZE >> sub_y) : (SB64_MIB_SIZE >> sub_y);
    const int32_t x_range = sb_size == BLOCK_128X128
        ? (MAX_MIB_SIZE >> sub_x) : (SB64_MIB_SIZE >> sub_x);
    void *blk_recon_buf;
    int32_t recon_stride;

    derive_blk_pointers(recon_picture_buf, plane,
        (mi_col*MI_SIZE >> sub_x), (mi_row*MI_SIZE >> sub_x),
        &blk_recon_buf, &recon_stride, sub_x, sub_y);

    for (int32_t y = 0; y < y_range; y += row_step) {
        uint8_t *p = (uint8_t*)blk_recon_buf +
                     ((y * MI_SIZE * recon_stride) << is16bit);

        for (int32_t x = 0; x < x_range;) {
            /*inner loop always filter vertical edges in a MI block. If MI size
            is 8x8, it will filter the vertical edge aligned with a 8x8 block.
            If 4x4 trasnform is used, it will then filter the internal edge
             aligned with a 4x4 block*/
            const uint32_t curr_luma_x = mi_col * MI_SIZE +
                                        ((x << sub_x) * MI_SIZE);
            const uint32_t curr_luma_y = mi_row * MI_SIZE +
                                        ((y << sub_y )* MI_SIZE);
            uint32_t advance_units;
            TxSize tx_size;
            AV1_DEBLOCKING_PARAMETERS params;
            memset(&params, 0, sizeof(params));
            tx_size =
                dec_set_lpf_parameters(&params, recon_picture_buf,
                    frm_hdr,color_config, lf_ctxt, lf_info, VERT_EDGE,
                    curr_luma_x, curr_luma_y, plane);

            if (tx_size == TX_INVALID) {
                params.filter_length = 0;
                tx_size = TX_4X4;
            }

            switch (params.filter_length) {
                /*apply 4-tap filtering*/
            case 4:
                if (is16bit)
                    aom_highbd_lpf_vertical_4(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        recon_picture_buf->bit_depth);
                else
                    aom_lpf_vertical_4(
                        p,
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
            case 6:  /* apply 6-tap filter for chroma plane only*/
                assert(plane != 0);
                if (is16bit)
                    aom_highbd_lpf_vertical_6(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        recon_picture_buf->bit_depth);
                else
                    aom_lpf_vertical_6(
                        p,
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                /*apply 8-tap filtering*/
            case 8:
                if (is16bit)
                    aom_highbd_lpf_vertical_8(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        recon_picture_buf->bit_depth);
                else
                    aom_lpf_vertical_8(
                        p,
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                /*apply 14-tap filtering*/
            case 14:
                if (is16bit)
                    aom_highbd_lpf_vertical_14(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        recon_picture_buf->bit_depth);
                else
                    aom_lpf_vertical_14(
                        p,
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                /*no filtering*/
            default: break;
            }
            /*advance the destination pointer*/
            assert(tx_size < TX_SIZES_ALL);
            advance_units = tx_size_wide_unit[tx_size];
            x += advance_units;
            p += ((advance_units * MI_SIZE) << is16bit);
        }
    }
}

/*It applies Horizonatal Loop Filtering in a superblock*/
void dec_av1_filter_block_plane_horz(
    FrameHeader *frm_hdr,EbColorConfig *color_config,
    EbPictureBufferDesc *recon_picture_buf,LFCtxt *lf_ctxt,
    LoopFilterInfoN *lf_info, const int32_t plane, BlockSize sb_size,
    const uint32_t mi_row, const uint32_t mi_col)
{
    EbBool is16bit = recon_picture_buf->bit_depth > 8;
    const int32_t col_step = MI_SIZE >> MI_SIZE_LOG2;
    const uint32_t sub_x = (plane > 0) ? color_config->subsampling_x : 0;
    const uint32_t sub_y = (plane > 0) ? color_config->subsampling_y : 0;
    const int32_t y_range = sb_size == BLOCK_128X128
        ? (MAX_MIB_SIZE >> sub_y) : (SB64_MIB_SIZE >> sub_y);
    const int32_t x_range = sb_size == BLOCK_128X128
        ? (MAX_MIB_SIZE >> sub_x) : (SB64_MIB_SIZE >> sub_x);
    void *blk_recon_buf;
    int32_t recon_stride;

    derive_blk_pointers(recon_picture_buf, plane,
        (mi_col*MI_SIZE >> sub_x), (mi_row*MI_SIZE >> sub_y),
        &blk_recon_buf, &recon_stride, sub_x, sub_y);

    for (int32_t x = 0; x < x_range; x += col_step) {
        uint8_t *p = (uint8_t*)blk_recon_buf + ((x * MI_SIZE) << is16bit);
        for (int32_t y = 0; y < y_range;) {
            /*inner loop always filter vertical edges in a MI block.If MI size
            is 8x8, it will first filter the vertical edge aligned with a 8x8
            block. If 4x4 trasnform is used, it will then filter the internal
            edge aligned with a 4x4 block*/
            const uint32_t curr_luma_x = mi_col * MI_SIZE +
                ((x << sub_x) * MI_SIZE);
            const uint32_t curr_luma_y = mi_row * MI_SIZE +
                ((y << sub_y)* MI_SIZE);
            uint32_t advance_units;
            TxSize tx_size;
            AV1_DEBLOCKING_PARAMETERS params;
            memset(&params, 0, sizeof(params));

            tx_size = dec_set_lpf_parameters(&params, recon_picture_buf, frm_hdr,
                    color_config, lf_ctxt, lf_info, HORZ_EDGE,
                    curr_luma_x,  curr_luma_y, plane);

            if (tx_size == TX_INVALID) {
                params.filter_length = 0;
                tx_size = TX_4X4;
            }

            switch (params.filter_length) {
                /*apply 4-tap filtering*/
            case 4:
                if (is16bit)
                    aom_highbd_lpf_horizontal_4(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        recon_picture_buf->bit_depth);
                else
                    aom_lpf_horizontal_4(
                        p,
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                /*apply 6-tap filtering*/
            case 6:
                assert(plane != 0);
                if (is16bit)
                    aom_highbd_lpf_horizontal_6(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        recon_picture_buf->bit_depth);
                else
                    aom_lpf_horizontal_6(
                        p,
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                /*apply 8-tap filtering*/
            case 8:
                if (is16bit)
                    aom_highbd_lpf_horizontal_8(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        recon_picture_buf->bit_depth);
                else
                    aom_lpf_horizontal_8(
                        p,
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                /*apply 14-tap filtering*/
            case 14:
                if (is16bit)
                    aom_highbd_lpf_horizontal_14(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        recon_picture_buf->bit_depth);
                else
                    aom_lpf_horizontal_14(
                        p,
                        recon_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                /* no filtering */
            default: break;
            }
            /*advance the destination pointer*/
            assert(tx_size < TX_SIZES_ALL);
            advance_units = tx_size_high_unit[tx_size];
            y += advance_units;
            p += ((advance_units * recon_stride * MI_SIZE) << is16bit);
        }
    }
}

/*LF function to filter each SB*/
void dec_loop_filter_sb(
    FrameHeader *frm_hdr, SeqHeader *seq_header,
    EbPictureBufferDesc *recon_picture_buf, LFCtxt *lf_ctxt,
    LoopFilterInfoN *lf_info, const uint32_t mi_row, const uint32_t mi_col,
    int32_t plane_start, int32_t plane_end, uint8_t LastCol)
{
    int32_t plane;
    for (plane = plane_start; plane < plane_end; plane++) {
        if (plane == 0 && !(frm_hdr->loop_filter_params.filter_level[0]) &&
            !(frm_hdr->loop_filter_params.filter_level[1]))
            break;
        else if (plane == 1 && !(frm_hdr->loop_filter_params.filter_level_u))
            continue;
        else if (plane == 2 && !(frm_hdr->loop_filter_params.filter_level_v))
            continue;

        if (frm_hdr->loop_filter_params.combine_vert_horz_lf) {
            /*filter all vertical and horizontal edges in every 64x64 super block
             filter vertical edges*/
            dec_av1_filter_block_plane_vert(frm_hdr, &seq_header->color_config,
                recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                mi_row, mi_col);

            /*filter horizontal edges*/
            int32_t max_mib_size = seq_header->sb_size == BLOCK_128X128 ?
                MAX_MIB_SIZE : SB64_MIB_SIZE;

            if ((int32_t)mi_col - max_mib_size >= 0) {
                dec_av1_filter_block_plane_horz(frm_hdr, &seq_header->color_config,
                    recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                    mi_row, mi_col - max_mib_size);
            }

            /*Filter the horizontal edges of the last lcu in each row*/
            if (LastCol) {
                dec_av1_filter_block_plane_horz(frm_hdr, &seq_header->color_config,
                    recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                    mi_row, mi_col);
            }
        }
        else {
            /*filter all vertical edges in every 64x64 super block*/
            dec_av1_filter_block_plane_vert(frm_hdr, &seq_header->color_config,
                recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                mi_row, mi_col);

            /*filter all horizontal edges in every 64x64 super block*/
            dec_av1_filter_block_plane_horz(frm_hdr, &seq_header->color_config,
                recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                mi_row, mi_col);
        }
    }
}

/*To update the LoopFilterInfo paramters*/
static void dec_update_sharpness(LoopFilterInfoN *lfi, int32_t sharpness_lvl) {
    int32_t lvl;

    // For each possible value for the loop filter fill out limits
    for (lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++) {
        // Set loop filter parameters that control sharpness.
        int32_t block_inside_limit = lvl >> ((sharpness_lvl > 0)
                                     + (sharpness_lvl > 4));

        if (sharpness_lvl > 0) {
            if (block_inside_limit > (9 - sharpness_lvl))
                block_inside_limit = (9 - sharpness_lvl);
        }

        if (block_inside_limit < 1) block_inside_limit = 1;

        memset(lfi->lfthr[lvl].lim, block_inside_limit, SIMD_WIDTH);
        memset(lfi->lfthr[lvl].mblim, (2 * (lvl + 2) + block_inside_limit),
            SIMD_WIDTH);
    }
}

/*Update the loop filter for the current frame */
static void dec_av1_loop_filter_frame_init(FrameHeader *frm_hdr,
    LoopFilterInfoN *lf_info, int32_t plane_start, int32_t plane_end)
{
    int32_t filt_lvl[MAX_MB_PLANE], filt_lvl_r[MAX_MB_PLANE];
    int32_t plane;
    int32_t seg_id;

    /*n_shift is the multiplier for lf_deltas
    the multiplier is 1 for when filter_lvl is between 0 and 31;
    2 when filter_lvl is between 32 and 63*/
    struct LoopFilter *const lf = &frm_hdr->loop_filter_params;

    lf->combine_vert_horz_lf = 1;

    /*update sharpness limits*/
    dec_update_sharpness(lf_info, lf->sharpness_level);

    /*init hev threshold const vectors*/
    for (int lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++)
        memset(lf_info->lfthr[lvl].hev_thr, (lvl >> 4), SIMD_WIDTH);


    filt_lvl[0] = frm_hdr->loop_filter_params.filter_level[0];
    filt_lvl[1] = frm_hdr->loop_filter_params.filter_level_u;
    filt_lvl[2] = frm_hdr->loop_filter_params.filter_level_v;

    filt_lvl_r[0] = frm_hdr->loop_filter_params.filter_level[1];
    filt_lvl_r[1] = frm_hdr->loop_filter_params.filter_level_u;
    filt_lvl_r[2] = frm_hdr->loop_filter_params.filter_level_v;

    for (plane = plane_start; plane < plane_end; plane++) {
        if (plane == 0 && !filt_lvl[0] && !filt_lvl_r[0])
            break;
        else if (plane == 1 && !filt_lvl[1])
            continue;
        else if (plane == 2 && !filt_lvl[2])
            continue;

        for (seg_id = 0; seg_id < MAX_SEGMENTS; seg_id++) {
            for (int32_t dir = 0; dir < 2; ++dir) {

                int32_t lvl_seg = (dir == 0) ? filt_lvl[plane]
                                  : filt_lvl_r[plane];
                assert(plane >= 0 && plane <= 2);

                //const int32_t seg_lf_feature_id = seg_lvl_lf_lut[plane][dir];
                /*if (segfeature_active(seg, seg_id, seg_lf_feature_id)) {
                    const int32_t data = get_segdata(&pcs_ptr->parent_pcs_ptr->seg,
                                            seg_id, seg_lf_feature_id);
                    lvl_seg = clamp(lvl_seg + data, 0, MAX_LOOP_FILTER);
                }*/

                if (!lf->mode_ref_delta_enabled) {
                    /*we could get rid of this if we assume that deltas are
                    set to zero when not in use; decoder always uses deltas*/
                    memset(lf_info->lvl[plane][seg_id][dir], lvl_seg,
                        sizeof(lf_info->lvl[plane][seg_id][dir]));
                }
                else {
                    int32_t ref, mode;
                    const int32_t scale = 1 << (lvl_seg >> 5);
                    const int32_t intra_lvl =
                        lvl_seg + lf->ref_deltas[INTRA_FRAME] * scale;
                    lf_info->lvl[plane][seg_id][dir][INTRA_FRAME][0] =
                        (uint8_t)clamp(intra_lvl, 0, MAX_LOOP_FILTER);

                    for (ref = LAST_FRAME; ref < REF_FRAMES; ++ref) {
                        for (mode = 0; mode < MAX_MODE_LF_DELTAS; ++mode) {
                            const int32_t inter_lvl =
                                lvl_seg + lf->ref_deltas[ref] * scale +
                                lf->mode_deltas[mode] * scale;
                            lf_info->lvl[plane][seg_id][dir][ref][mode] =
                                (uint8_t)clamp(inter_lvl, 0, MAX_LOOP_FILTER);
                        }
                    }
                }
            }
        }
    }
}

/*Frame level function to trigger loop filter for each superblock*/
void dec_av1_loop_filter_frame(
    FrameHeader *frm_hdr, SeqHeader *seq_header,
    EbPictureBufferDesc *recon_picture_buf, LFCtxt *lf_ctxt,
    int32_t plane_start, int32_t plane_end) {

    uint8_t     sb_size_Log2 = seq_header->sb_size_log2;
    uint32_t    x_lcu_index;
    uint32_t    y_lcu_index;
    uint32_t    sb_origin_x;
    uint32_t    sb_origin_y;
    EbBool      endOfRowFlag;

    LoopFilterInfoN *lf_info = &lf_ctxt->lf_info;

    int32_t sb_size_w = block_size_wide[seq_header->sb_size];
    int32_t sb_size_h = block_size_high[seq_header->sb_size];
    uint32_t picture_width_in_sb    =
        (seq_header->max_frame_width + sb_size_w - 1) / sb_size_w;
    uint32_t picture_height_in_sb   =
        (seq_header->max_frame_height + sb_size_h- 1) / sb_size_h;

    dec_av1_loop_filter_frame_init(frm_hdr, lf_info, plane_start, plane_end);

    /*Loop over a frame : tregger dec_loop_filter_sb for each SB*/
    for (y_lcu_index = 0; y_lcu_index < picture_height_in_sb; ++y_lcu_index) {
        for (x_lcu_index = 0; x_lcu_index < picture_width_in_sb; ++x_lcu_index) {
            sb_origin_x = x_lcu_index << sb_size_Log2;
            sb_origin_y = y_lcu_index << sb_size_Log2;
            endOfRowFlag = (x_lcu_index == picture_width_in_sb - 1) ?
                EB_TRUE : EB_FALSE;

            /*LF function for a SB*/
            dec_loop_filter_sb(frm_hdr, seq_header, recon_picture_buf,
                lf_ctxt, lf_info, sb_origin_y >> 2, sb_origin_x >> 2,
                plane_start, plane_end, endOfRowFlag);
        }
    }
}
