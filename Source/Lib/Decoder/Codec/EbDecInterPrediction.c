/*
* Copyright(c) 2019 Intel Corporation
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

#include <stdlib.h>
#include <string.h>

#include "EbInterPrediction.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbObuParse.h"
#include "EbDecParseHelper.h"

#include "EbDecPicMgr.h"
#include "EbDecNbr.h"
#include "EbDecUtils.h"

static INLINE void dec_clamp_mv(MV *mv, int32_t min_col, int32_t max_col, int32_t min_row,
    int32_t max_row) {
    mv->col = (int16_t)clamp(mv->col, min_col, max_col);
    mv->row = (int16_t)clamp(mv->row, min_row, max_row);
}

static INLINE MV dec_clamp_mv_to_umv_border_sb (
    int32_t mb_to_left_edge, int32_t mb_to_right_edge, int32_t mb_to_top_edge,
    int32_t mb_to_bottom_edge, const MV *src_mv, int32_t bw, int32_t bh,
    int32_t ss_x, int32_t ss_y)
{
    // If the MV points so far into the UMV border that no visible pixels
    // are used for reconstruction, the subpel part of the MV can be
    // discarded and the MV limited to 16 pixels with equivalent results.
    const int32_t spel_left = (AOM_INTERP_EXTEND + bw) << SUBPEL_BITS;
    const int32_t spel_right = spel_left - SUBPEL_SHIFTS;
    const int32_t spel_top = (AOM_INTERP_EXTEND + bh) << SUBPEL_BITS;
    const int32_t spel_bottom = spel_top - SUBPEL_SHIFTS;
    MV clamped_mv = { (int16_t)(src_mv->row * (1 << (1 - ss_y))),
        (int16_t)(src_mv->col * (1 << (1 - ss_x))) };
    assert(ss_x <= 1);
    assert(ss_y <= 1);

    dec_clamp_mv(&clamped_mv,
        mb_to_left_edge   * (1 << (1 - ss_x)) - spel_left,
        mb_to_right_edge  * (1 << (1 - ss_x)) + spel_right,
        mb_to_top_edge    * (1 << (1 - ss_y)) - spel_top,
        mb_to_bottom_edge * (1 << (1 - ss_y)) + spel_bottom);

    return clamped_mv;
}


void svtav1_predict_inter_block_plane(
    EbDecHandle *dec_hdl, PartitionInfo_t *part_info, int32_t plane,
    int32_t build_for_obmc, int32_t mi_x, int32_t mi_y,
    void *dst, int32_t dst_stride,
    int32_t some_use_intra, int32_t bit_depth)
{
    const ModeInfo_t *mi = part_info->mi;
    const FrameHeader *cur_frm_hdr = &dec_hdl->frame_header;
    int32_t is_compound = has_second_ref(mi);
    int32_t ref;
    const int32_t is_intrabc = is_intrabc_block(mi);

    //temporary buffer for joint compound, move this to context if stack does not hold.
    DECLARE_ALIGNED(32, uint16_t, tmp_dst[128 * 128]);

    int32_t highbd = bit_depth > EB_8BIT;

    const BlockSize bsize = mi->sb_type;
    const int32_t ss_x = plane ? part_info->subsampling_x : 0;
    const int32_t ss_y = plane ? part_info->subsampling_y : 0;
    int32_t bw = part_info->wpx[0] >> ss_x;
    int32_t bh = part_info->hpx[0] >> ss_y;
    int32_t row_start = 0;
    int32_t col_start = 0;
    // For sub8x8 chroma blocks, we may be covering more than one luma block's
    // worth of pixels. Thus (mi_x, mi_y) may not be the correct coordinates for
    // the top-left corner of the prediction source - the correct top-left corner
    // is at (pre_x, pre_y).
    if (some_use_intra && plane && !build_for_obmc) {
        bool sub8_w = (block_size_wide[bsize] == 4) && ss_x;
        bool sub8_h = (block_size_high[bsize] == 4) && ss_y;

        if (sub8_w) {
            if (part_info->mi_col & 0x1) {
                bw = bw << 1;
                col_start = -1;
            }
            else
                return;
        }

        if (sub8_h) {
            if (part_info->mi_row & 0x1) {
                bh = bh << 1;
                row_start = -1;
            }
            else
                return;
        }
    }
    const int32_t pre_x = (mi_x + MI_SIZE * col_start) >> ss_x;
    const int32_t pre_y = (mi_y + MI_SIZE * row_start) >> ss_y;

    int32_t dst_offset = ((MI_SIZE * col_start) >> ss_x) +
            ((MI_SIZE * row_start * dst_stride) >> ss_y);
    void *dst_mod = (void *)((uint8_t *)dst + (dst_offset << highbd));

    assert(IMPLIES(is_intrabc, !is_compound));
    {

        ConvolveParams conv_params = get_conv_params_no_round(0, 0,
            plane, tmp_dst, MAX_SB_SIZE, is_compound, bit_depth);

        /* TODO: support Distantance WTD compound inter prediction later
        av1_dist_wtd_comp_weight_assign(
            cm, mi, 0, &conv_params.fwd_offset, &conv_params.bck_offset,
            &conv_params.use_dist_wtd_comp_avg, is_compound);
        */

        for (ref = 0; ref < 1 + is_compound; ++ref) {
            ScaleFactors *const sf = NULL; //TODO: Add reference scaling support later
                /* is_intrabc ? &cm->sf_identity : xd->block_ref_scale_factors[ref]; */
            const int32_t mode = mi->mode;
            const MV mv = mi->mv[ref].as_mv;
            MV mv_q4;
            const EbWarpedMotionParams *const wm_global =
                &part_info->ps_global_motion[mi->ref_frame[ref]];
            const EbWarpedMotionParams *const wm_local =
                &part_info->local_warp_params;
            EbDecPicBuf *ref_buf  = get_ref_frame_buf(dec_hdl, mi->ref_frame[ref]);
            EbPictureBufferDesc *ps_ref_pic_buf = ref_buf->ps_pic_buf;
            SubpelParams subpel_params;

            int32_t do_warp = (bw >= 8 && bh >= 8 && !build_for_obmc &&
                (cur_frm_hdr->force_integer_mv == 0) &&
                (((mode == GLOBALMV || mode == GLOBAL_GLOBALMV) &&
                (wm_global->wmtype > TRANSLATION)) ||
                    (mi->motion_mode == WARPED_CAUSAL)));

            void   *src;
            int32_t src_stride;

            derive_blk_pointers(ps_ref_pic_buf, plane, 0, 0, &src, &src_stride, ss_x, ss_y);

            mv_q4 = dec_clamp_mv_to_umv_border_sb(
                    part_info->mb_to_left_edge,
                    part_info->mb_to_right_edge,
                    part_info->mb_to_top_edge,
                    part_info->mb_to_bottom_edge,
                    &mv, bw, bh, ss_x, ss_y);

            subpel_params.xs = 0;
            subpel_params.ys = 0;
            subpel_params.subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_params.subpel_y = mv_q4.row & SUBPEL_MASK;

            conv_params.do_average = ref;
            /* TODO: support masked inter prediction based on WEDGE / DIFFWTD compound type later
            if (is_masked_compound_type(mi->inter_compound.type)) {
                // masked compound type has its own average mechanism
                conv_params.do_average = 0;
            }

            if (ref && is_masked_compound_type(mi->inter_compound.type))
                av1_make_masked_inter_predictor(
                    pre[ref], src_stride[ref], dst_mod, dst_buf->stride,
                    &subpel_params[ref], sf, bw, bh, &conv_params, mi->interp_filters,
                    plane, &warp_types, pre_x,
                    pre_y, ref, xd, cm->allow_warped_motion);
            else
             */

            assert(IMPLIES(is_intrabc, !do_warp));

            if (do_warp) {
                const EbWarpedMotionParams *wm_params = &default_warp_params;

                wm_params = (mi->motion_mode == WARPED_CAUSAL) ? wm_local : wm_global;

                av1_warp_plane((EbWarpedMotionParams *)wm_params,
                    highbd, bit_depth,
                    src, ref_buf->ps_pic_buf->width >> ss_x,
                    ref_buf->ps_pic_buf->height >> ss_y,
                    src_stride, dst_mod, pre_x, pre_y, bw, bh, dst_stride,
                    ss_x, ss_y, &conv_params);
            }
            else if (highbd) {
                uint16_t *src16 = (uint16_t *)src +
                    (((pre_y) + (mv_q4.row >> SUBPEL_BITS)) * src_stride) +
                    (pre_x) + (mv_q4.col >> SUBPEL_BITS);

                svt_highbd_inter_predictor(src16, src_stride, dst_mod, dst_stride,
                    &subpel_params, sf, bw, bh, &conv_params,
                    mi->interp_filters, is_intrabc, bit_depth);
            }
            else {
                uint8_t *src8 = (uint8_t *)src +
                    (((pre_y) + (mv_q4.row >> SUBPEL_BITS)) * src_stride) +
                    (pre_x) + (mv_q4.col >> SUBPEL_BITS);

                svt_inter_predictor(src8, src_stride, dst_mod, dst_stride,
                    &subpel_params, sf, bw, bh, &conv_params,
                    mi->interp_filters, is_intrabc);
            }
        }
    }
}

void svtav1_predict_inter_block(
    EbDecHandle *dec_hdl, PartitionInfo_t *part_info,
    int32_t mi_row, int32_t mi_col, int32_t num_planes)
{
    void *blk_recon_buf;
    int32_t recon_strd;
    int32_t sub_x, sub_y;
    int32_t some_use_intra;

    EbPictureBufferDesc *recon_picture_buf = dec_hdl->cur_pic_buf[0]->ps_pic_buf;

    /* scan through sub 8 blocks and see if anyof them is intra */
    some_use_intra = 0;
    const BlockSize bsize = part_info->mi->sb_type;
    sub_x = part_info->subsampling_x;
    sub_y = part_info->subsampling_y;
    bool sub8_w = (block_size_wide[bsize] == 4) && sub_x;
    bool sub8_h = (block_size_high[bsize] == 4) && sub_y;
    if (sub8_h || sub8_w) {
        int32_t i, j;
        /* Floor and Ceil to nearest 8x8 blks */
        const int row_start = sub8_h ? (mi_row & (~1)) : mi_row;
        const int row_end   = sub8_h ? (mi_row | (1)) : mi_row;
        const int col_start = sub8_w ? (mi_col & (~1)) : mi_col;
        const int col_end   = sub8_w ? (mi_col | (1)) : mi_col;

        for (i = row_start; i <= row_end; i++) {
            for (j = col_start; j <= col_end; j++) {
                ModeInfo_t *mode_info = get_cur_mode_info(dec_hdl,
                    i, j, part_info->sb_info);
                if (mode_info->ref_frame[0] == INTRA_FRAME)
                    some_use_intra = 1;
            }
        }
    }

    for (int plane = 0; plane < num_planes; ++plane) {
        sub_x = (plane > 0) ? part_info->subsampling_x : 0;
        sub_y = (plane > 0) ? part_info->subsampling_y : 0;

        derive_blk_pointers(recon_picture_buf, plane,
            mi_col*MI_SIZE >> sub_x, mi_row*MI_SIZE >> sub_y,
            &blk_recon_buf, &recon_strd, sub_x, sub_y);

        svtav1_predict_inter_block_plane(dec_hdl, part_info, plane,
            0, mi_col*MI_SIZE, mi_row*MI_SIZE, blk_recon_buf, recon_strd,
            some_use_intra, recon_picture_buf->bit_depth);
    }

    return;
}
