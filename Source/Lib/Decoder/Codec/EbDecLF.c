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
#include "EbDecParseHelper.h"

#define FILTER_LEN 4

/*Filter_length is mapped to int indx for the filter tap arrays*/
/*4-> 0   6-> 1   8-> 2   14->3 */
static int8_t filter_map[15] = { -1, -1, -1, -1, 0,-1, 1, -1, 2, -1,
                                -1, -1, -1, -1, 3 };

svt_lbd_filter_tap_fn_t lbd_vert_filter_tap[FILTER_LEN];
svt_hbd_filter_tap_fn_t hbd_vert_filter_tap[FILTER_LEN];
svt_lbd_filter_tap_fn_t lbd_horz_filter_tap[FILTER_LEN];
svt_hbd_filter_tap_fn_t hbd_horz_filter_tap[FILTER_LEN];

void set_lbd_lf_filter_tap_functions(void)
{
    lbd_horz_filter_tap[0] = aom_lpf_horizontal_4;
    lbd_horz_filter_tap[1] = aom_lpf_horizontal_6;
    lbd_horz_filter_tap[2] = aom_lpf_horizontal_8;
    lbd_horz_filter_tap[3] = aom_lpf_horizontal_14;

    lbd_vert_filter_tap[0] = aom_lpf_vertical_4;
    lbd_vert_filter_tap[1] = aom_lpf_vertical_6;
    lbd_vert_filter_tap[2] = aom_lpf_vertical_8;
    lbd_vert_filter_tap[3] = aom_lpf_vertical_14;
}

void set_hbd_lf_filter_tap_functions(void)
{
    hbd_horz_filter_tap[0] = aom_highbd_lpf_horizontal_4;
    hbd_horz_filter_tap[1] = aom_highbd_lpf_horizontal_6;
    hbd_horz_filter_tap[2] = aom_highbd_lpf_horizontal_8;
    hbd_horz_filter_tap[3] = aom_highbd_lpf_horizontal_14;

    hbd_vert_filter_tap[0] = aom_highbd_lpf_vertical_4;
    hbd_vert_filter_tap[1] = aom_highbd_lpf_vertical_6;
    hbd_vert_filter_tap[2] = aom_highbd_lpf_vertical_8;
    hbd_vert_filter_tap[3] = aom_highbd_lpf_vertical_14;
}

/*Population of neighbour block LUMA params for each 4x4 block*/
void fill_4x4_param_luma(LFBlockParamL* lf_block_l,
    int32_t tu_x, int32_t tu_y, int32_t stride,
    TxSize tx_size, BlockModeInfo *mode_info)
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

/*Function to get transform size*/
static INLINE TxSize dec_get_transform_size(const EDGE_DIR edge_dir, TxSize tx_size) {
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
static INLINE TxSize dec_set_lpf_parameters(AV1_DEBLOCKING_PARAMETERS *const params,
    EbPictureBufferDesc *recon_picture_buf, FrameHeader *frm_hdr,
    EbColorConfig *color_config, LFCtxt *lf_ctxt, LoopFilterInfoN *lf_info,
    const EDGE_DIR edge_dir, const uint32_t x,
    const uint32_t y, const int32_t plane,
    int32_t *sb_delta_lf, int32_t *sb_delta_lf_prev)
{
    UNUSED(recon_picture_buf);
    /*reset to initial values*/
    params->filter_length = 0;

    /*no deblocking is required*/
    const uint32_t width = frm_hdr->frame_size.frame_width;
    const uint32_t height = frm_hdr->frame_size.frame_height;
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
            const uint32_t curr_level = get_filter_level(frm_hdr, lf_info,
                edge_dir, plane, sb_delta_lf,
                lf_block_l_cur->segment_id, lf_block_l_cur->mode,
                lf_block_l_cur->ref_frame_0);

            const int32_t curr_skipped = lf_block_l_cur->skip
                && is_inter_block_no_intrabc(lf_block_l_cur->ref_frame_0);

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
                    get_filter_level(frm_hdr, lf_info,
                        edge_dir, plane, sb_delta_lf_prev,
                        lf_block_l_prev->segment_id, lf_block_l_prev->mode,
                        lf_block_l_prev->ref_frame_0);

                const int32_t pv_skip = lf_block_l_prev->skip &&
                    is_inter_block_no_intrabc(lf_block_l_prev->ref_frame_0);
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
    FrameHeader *frm_hdr, EbColorConfig *color_config,
    EbPictureBufferDesc *recon_picture_buf, LFCtxt *lf_ctxt,
    LoopFilterInfoN *lf_info, const int32_t plane, BlockSize sb_size,
    const uint32_t mi_row, const uint32_t mi_col, int32_t *sb_delta_lf)
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
    int32_t *sb_delta_lf_left;

    derive_blk_pointers(recon_picture_buf, plane,
        (mi_col*MI_SIZE >> sub_x), (mi_row*MI_SIZE >> sub_y),
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
                ((y << sub_y)* MI_SIZE);
            uint32_t advance_units;
            TxSize tx_size;
            AV1_DEBLOCKING_PARAMETERS params;
            memset(&params, 0, sizeof(params));

            /* For VERT_EDGE edge and x_range is for SB scan */
            sb_delta_lf_left = x == 0 ? sb_delta_lf - FRAME_LF_COUNT : sb_delta_lf;

            tx_size =
                dec_set_lpf_parameters(&params, recon_picture_buf,
                    frm_hdr, color_config, lf_ctxt, lf_info, VERT_EDGE,
                    curr_luma_x, curr_luma_y, plane,
                    sb_delta_lf, sb_delta_lf_left);

            if (tx_size == TX_INVALID) {
                params.filter_length = 0;
                tx_size = TX_4X4;
            }

            int8_t filter_idx = filter_map[params.filter_length];

            if (filter_idx != -1) {
                if (is16bit)
                    hbd_vert_filter_tap[filter_idx](
                        (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
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
    FrameHeader *frm_hdr, EbColorConfig *color_config,
    EbPictureBufferDesc *recon_picture_buf, LFCtxt *lf_ctxt,
    LoopFilterInfoN *lf_info, const int32_t plane, BlockSize sb_size,
    const uint32_t mi_row, const uint32_t mi_col, int32_t *sb_delta_lf)
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
    int32_t *sb_delta_lf_above;

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

            /* For HORZ_EDGE edge and y_range is for SB scan */
            sb_delta_lf_above = y == 0 ? sb_delta_lf - lf_ctxt->delta_lf_stride :
                sb_delta_lf;

            tx_size = dec_set_lpf_parameters(&params, recon_picture_buf, frm_hdr,
                color_config, lf_ctxt, lf_info, HORZ_EDGE,
                curr_luma_x, curr_luma_y, plane,
                sb_delta_lf, sb_delta_lf_above);

            if (tx_size == TX_INVALID) {
                params.filter_length = 0;
                tx_size = TX_4X4;
            }

            int filter_idx = filter_map[params.filter_length];

            if (filter_idx != -1) {
                if (is16bit)
                    hbd_horz_filter_tap[filter_idx](
                        (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
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
    int32_t plane_start, int32_t plane_end, uint8_t LastCol, int32_t *sb_delta_lf)
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
                mi_row, mi_col, sb_delta_lf);

            /*filter horizontal edges*/
            int32_t max_mib_size = seq_header->sb_size == BLOCK_128X128 ?
                MAX_MIB_SIZE : SB64_MIB_SIZE;

            if ((int32_t)mi_col - max_mib_size >= 0) {
                dec_av1_filter_block_plane_horz(frm_hdr, &seq_header->color_config,
                    recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                    mi_row, mi_col - max_mib_size, (sb_delta_lf-FRAME_LF_COUNT));
            }

            /*Filter the horizontal edges of the last lcu in each row*/
            if (LastCol) {
                dec_av1_filter_block_plane_horz(frm_hdr, &seq_header->color_config,
                    recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                    mi_row, mi_col, sb_delta_lf);
            }
        }
        else {
            /*filter all vertical edges in every 64x64 super block*/
            dec_av1_filter_block_plane_vert(frm_hdr, &seq_header->color_config,
                recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                mi_row, mi_col, sb_delta_lf);

            /*filter all horizontal edges in every 64x64 super block*/
            dec_av1_filter_block_plane_horz(frm_hdr, &seq_header->color_config,
                recon_picture_buf, lf_ctxt, lf_info, plane, seq_header->sb_size,
                mi_row, mi_col, sb_delta_lf);
        }
    }
}

/*Frame level function to trigger loop filter for each superblock*/
void dec_av1_loop_filter_frame(EbDecHandle *dec_handle_ptr,
    EbPictureBufferDesc *recon_picture_buf, LFCtxt *lf_ctxt,
    int32_t plane_start, int32_t plane_end)
{
    FrameHeader *frm_hdr = &dec_handle_ptr->frame_header;
    SeqHeader *seq_header = &dec_handle_ptr->seq_header;
    uint8_t     sb_size_Log2 = seq_header->sb_size_log2;
    uint32_t    x_lcu_index;
    uint32_t    y_lcu_index;
    uint32_t    sb_origin_x;
    uint32_t    sb_origin_y;
    EbBool      endOfRowFlag;

    LoopFilterInfoN *lf_info = &lf_ctxt->lf_info;
    lf_ctxt->delta_lf_stride = dec_handle_ptr->master_frame_buf.sb_cols *
                               FRAME_LF_COUNT;

    int32_t sb_size_w = block_size_wide[seq_header->sb_size];
    int32_t sb_size_h = block_size_high[seq_header->sb_size];
    uint32_t picture_width_in_sb    =
        (seq_header->max_frame_width + sb_size_w - 1) / sb_size_w;
    uint32_t picture_height_in_sb   =
        (seq_header->max_frame_height + sb_size_h- 1) / sb_size_h;

    frm_hdr->loop_filter_params.combine_vert_horz_lf = 1;
    /*init hev threshold const vectors*/
    for (int lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++)
        memset(lf_info->lfthr[lvl].hev_thr, (lvl >> 4), SIMD_WIDTH);

    eb_av1_loop_filter_frame_init(frm_hdr, lf_info, plane_start, plane_end);

    set_lbd_lf_filter_tap_functions();
    set_hbd_lf_filter_tap_functions();

    /*Loop over a frame : tregger dec_loop_filter_sb for each SB*/
    for (y_lcu_index = 0; y_lcu_index < picture_height_in_sb; ++y_lcu_index) {
        for (x_lcu_index = 0; x_lcu_index < picture_width_in_sb; ++x_lcu_index) {
            sb_origin_x = x_lcu_index << sb_size_Log2;
            sb_origin_y = y_lcu_index << sb_size_Log2;
            endOfRowFlag = (x_lcu_index == picture_width_in_sb - 1) ?
                EB_TRUE : EB_FALSE;

            MasterFrameBuf *master_frame_buf = &dec_handle_ptr->master_frame_buf;
            CurFrameBuf    *frame_buf = &master_frame_buf->cur_frame_bufs[0];

            SBInfo  *sb_info = frame_buf->sb_info + (
                ((y_lcu_index * master_frame_buf->sb_cols) + x_lcu_index));

            /*sb_info->sb_delta_lf = frame_buf->delta_lf + (FRAME_LF_COUNT *
                ((y_lcu_index * master_frame_buf->sb_cols) + x_lcu_index));*/

            /*LF function for a SB*/
            dec_loop_filter_sb(frm_hdr, seq_header, recon_picture_buf,
                lf_ctxt, lf_info, sb_origin_y >> 2, sb_origin_x >> 2,
                plane_start, plane_end, endOfRowFlag, sb_info->sb_delta_lf);
        }
    }
}
