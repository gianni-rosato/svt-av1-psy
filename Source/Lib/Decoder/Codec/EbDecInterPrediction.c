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
#include <stdbool.h>

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbObuParse.h"

#include "EbDecPicMgr.h"
#include "EbDecNbr.h"
#include "EbDecUtils.h"
#include "EbDecObmc.h"

#include "common_dsp_rtcd.h"

#include "EbDecProcessFrame.h"
#include "EbDecIntraPrediction.h"
#include "EbDecInterPrediction.h"
#include "EbUtility.h"
#include "EbDefinitions.h"
#include "EbWarpedMotion.h"

static INLINE void dec_clamp_mv(MV *mv, int32_t min_col, int32_t max_col, int32_t min_row,
                                int32_t max_row) {
    mv->col = (int16_t)clamp(mv->col, min_col, max_col);
    mv->row = (int16_t)clamp(mv->row, min_row, max_row);
}

static INLINE MV dec_clamp_mv_to_umv_border_sb(int32_t mb_to_left_edge, int32_t mb_to_right_edge,
                                               int32_t mb_to_top_edge, int32_t mb_to_bottom_edge,
                                               const MV *src_mv, int32_t bw, int32_t bh,
                                               int32_t ss_x, int32_t ss_y) {
    // If the MV points so far into the UMV border that no visible pixels
    // are used for reconstruction, the subpel part of the MV can be
    // discarded and the MV limited to 16 pixels with equivalent results.
    const int32_t spel_left   = (AOM_INTERP_EXTEND + bw) << SUBPEL_BITS;
    const int32_t spel_right  = spel_left - SUBPEL_SHIFTS;
    const int32_t spel_top    = (AOM_INTERP_EXTEND + bh) << SUBPEL_BITS;
    const int32_t spel_bottom = spel_top - SUBPEL_SHIFTS;
    MV            clamped_mv  = {(int16_t)(src_mv->row * (1 << (1 - ss_y))),
                     (int16_t)(src_mv->col * (1 << (1 - ss_x)))};
    assert(ss_x <= 1);
    assert(ss_y <= 1);

    dec_clamp_mv(&clamped_mv,
                 mb_to_left_edge * (1 << (1 - ss_x)) - spel_left,
                 mb_to_right_edge * (1 << (1 - ss_x)) + spel_right,
                 mb_to_top_edge * (1 << (1 - ss_y)) - spel_top,
                 mb_to_bottom_edge * (1 << (1 - ss_y)) + spel_bottom);

    return clamped_mv;
}

#define OPT_DYN_PAD        0
void *aom_memset16(void *dest, int val, size_t length) {
    size_t i;
    uint16_t *dest16 = (uint16_t *)dest;
    for (i = 0; i < length; i++) *dest16++ = val;
    return dest;
}

static void highbd_build_mc_border(const uint8_t *src8, int32_t src_stride,
    uint8_t *dst8, int32_t dst_stride, int32_t x, int32_t y,
    int32_t b_w, int32_t b_h, int32_t w, int32_t h)
{
    const uint16_t *ref_row = (uint16_t *) src8;
    uint16_t *dst = (uint16_t *) dst8;

#if OPT_DYN_PAD
    int32_t x0 = x, y0 = y;
    uint16_t *dst_temp = dst;
    const uint16_t *ref_temp = ref_row;
    int32_t blk_fill_itt = 0;

    int32_t left = x0 < 0 ? -x0 : 0;
    if (left > b_w) left = b_w;

    int32_t right = b_w - left;

    int32_t top = y0 < 0 ? -y0 : 0;
    if (top > b_h) top = b_h;

    int32_t bottom = b_h - top;

    if (y0 + b_h <= 0) {
        blk_fill_itt = top;
        top = 0;
    }
    else if (y0 < 0) {
        dst_temp = dst + (top * dst_stride);
        blk_fill_itt = bottom;
        bottom = 0;
    }
    else if (y0 > h) {
        ref_temp = ref_row + (h - 1)* src_stride;
        blk_fill_itt = bottom;
        bottom = 0;
}
    else {
        ref_temp = ref_row + (y0 * src_stride);
        if (y0 + b_h > h)
            blk_fill_itt = h - y0;
        else {
            blk_fill_itt = bottom;
            bottom = 0;
        }
    }

    if (left) {
        y0 = y0 + top;
        for (int32_t i = 0; i < blk_fill_itt; i++) {
            aom_memset16((void *)dst_temp, ref_temp[0], left);
            if (right) memcpy(dst_temp + left, ref_temp, right * sizeof(uint16_t));
            dst_temp += dst_stride;
            ++y0;
            if (y0 > 0 && y0 < h) ref_temp += src_stride;
        }
        right = 0;
    }

    if (right) {
        int32_t copy = 0, set = 0, ptr_ext = 0;
        y0 = y0 + top;
        if (x0 + b_w < w) copy = right;
        else if (x0 > w) set = right;
        else {
            copy = w - x0;
            set = b_w - copy;
            ptr_ext = copy;
        }
        for (int32_t i = 0; i < blk_fill_itt; i++) {
            if (copy) memcpy(dst_temp, ref_temp + x0, copy * sizeof(uint16_t));
            if (set) aom_memset16((void *)(dst_temp + ptr_ext), ref_temp[w - 1], set);
            dst_temp += dst_stride;
            ++y0;
            if (y0 > 0 && y0 < h) ref_temp += src_stride;
        }
    }

    if (top) {
        dst_temp = dst + (top * dst_stride);
        for (int32_t i = 0; i < top; i++) {
            memcpy(dst, dst_temp, b_w * sizeof(uint16_t));
            dst += dst_stride;
        }
    }

    if (bottom) {
        int32_t itt = b_h - (h - y);
        dst = dst_temp;
        dst_temp -= dst_stride;
        for (int32_t i = 0; i < itt; i++) {
            memcpy(dst, dst_temp, b_w * sizeof(uint16_t));
            dst += dst_stride;
        }
    }

#else //OPT_DYN_PAD
    if (y >= h)
        ref_row += (h - 1) * src_stride;
    else if (y > 0)
        ref_row += y * src_stride;

    do {
        int right = 0, copy;
        int left = x < 0 ? -x : 0;

        if (left > b_w) left = b_w;

        if (x + b_w > w) right = x + b_w - w;

        if (right > b_w) right = b_w;

        copy = b_w - left - right;

        if (left) aom_memset16((void*)dst, ref_row[0], left);

        if (copy) memcpy((void*)(dst + left), ref_row + x + left, copy * sizeof(uint16_t));

        if (right) aom_memset16((void*)(dst + left + copy), ref_row[w - 1], right);

        dst += dst_stride;
        ++y;

        if (y > 0 && y < h) ref_row += src_stride;
    } while (--b_h);
#endif
}

static void build_mc_border(const uint8_t *ref_row, int32_t src_stride,
    uint8_t *dst, int32_t dst_stride, int32_t x, int32_t y, int32_t b_w,
    int32_t b_h, int32_t w, int32_t h)
{
    //x,y = block start postion in ref frame
    //b_w, b_h  = block size in ref frame
    //w, h = ref frame width  & height

#if OPT_DYN_PAD
    int32_t x0 = x, y0 = y;
    uint8_t *dst_temp = dst;
    const uint8_t *ref_temp = ref_row;
    int32_t blk_fill_itt = 0;

    int32_t left = x0 < 0 ? -x0 : 0;
    if (left > b_w) left = b_w;

    int32_t right = b_w - left;

    int32_t top = y0 < 0 ? -y0 : 0;
    if (top > b_h) top = b_h;

    int32_t bottom = b_h - top;

    if (y0 + b_h <= 0) {
        blk_fill_itt = top;
        top = 0;
    }
    else if (y0 < 0) {
        dst_temp = dst + (top * dst_stride);
        blk_fill_itt = bottom;
        bottom = 0;
    }
    else if (y0 > h) {
        ref_temp = ref_row + (h - 1)* src_stride;
        blk_fill_itt = bottom;
        bottom = 0;
    }
    else {
        ref_temp = ref_row + (y0 * src_stride);
        if (y0 + b_h > h)
            blk_fill_itt = h - y0;
        else {
            blk_fill_itt = bottom;
            bottom = 0;
        }
    }

    if (left) {
        y0 = y0 + top;
        for (int32_t i = 0; i < blk_fill_itt; i++) {
            memset(dst_temp, ref_temp[0], left);
            if (right) memcpy(dst_temp + left, ref_temp, right);
            dst_temp += dst_stride;
            ++y0;
            if (y0 > 0 && y0 < h) ref_temp += src_stride;
        }
        right = 0;
    }

    if (right) {
        int32_t copy = 0, set = 0, ptr_ext = 0;
        y0 = y0 + top;
        if (x0 + b_w < w) copy = right;
        else if (x0 > w) set = right;
        else {
            copy = w - x0;
            set = b_w - copy;
            ptr_ext = copy;
        }
        for (int32_t i = 0; i < blk_fill_itt; i++) {
            if (copy) memcpy(dst_temp, ref_temp + x0, copy);
            if (set) memset(dst_temp + ptr_ext, ref_temp[w - 1], set);
            dst_temp += dst_stride;
            ++y0;
            if (y0 > 0 && y0 < h) ref_temp +=  src_stride;
        }
    }

    if (top) {
        dst_temp = dst + (top * dst_stride);
        for (int32_t i = 0; i < top; i++) {
            memcpy(dst, dst_temp, b_w);
            dst += dst_stride;
        }
    }

    if (bottom) {
        int32_t itt = b_h - (h - y);
        dst = dst_temp;
        dst_temp -= dst_stride;
        for (int32_t i = 0; i < itt; i++) {
            memcpy(dst, dst_temp, b_w);
            dst += dst_stride;
        }
    }

#else //OPT_DYN_PAD
    if (y >= h)
        ref_row += (h - 1) * src_stride;
    else if (y > 0)
        ref_row += y * src_stride;

    do {
        int right = 0, copy;
        int left = x < 0 ? -x : 0;

        if (left > b_w) left = b_w;

        if (x + b_w > w) right = x + b_w - w;

        if (right > b_w) right = b_w;

        copy = b_w - left - right;

        if (left) memset(dst, ref_row[0], left);

        if (copy) memcpy(dst + left, ref_row + x + left, copy);

        if (right) memset(dst + left + copy, ref_row[w - 1], right);

        dst += dst_stride;
        ++y;

        if (y > 0 && y < h) ref_row += src_stride;
    } while (--b_h);
#endif //OPT_DYN_PAD
}

static INLINE int update_extend_mc_border_params(PadBlock *block,
     EbDecPicBuf *ref_buf, MV32 scaled_mv, const ScaleFactors *sf,
    int32_t *x_pad, int32_t *y_pad,
    int32_t ss_x, int32_t ss_y)
{
    int frame_width = (ref_buf->superres_upscaled_width + ss_x) >> ss_x;
    int frame_height = (ref_buf->frame_height + ss_y) >> ss_y;


    // Skip border extension if block is inside the frame.
    if (block->x0 < -(DYNIMIC_PAD_VALUE >> ss_x) ||
        block->x1 >((frame_width - 1) + (DYNIMIC_PAD_VALUE >> ss_x)) ||
        block->y0 < -(DYNIMIC_PAD_VALUE >> ss_y) ||
        block->y1 >((frame_height - 1) + (DYNIMIC_PAD_VALUE >> ss_y)))
    {
        if ((scaled_mv.col & SUBPEL_MASK) || (sf->x_step_q4 != SUBPEL_SHIFTS)) {
            block->x0 -= AOM_INTERP_EXTEND;
            block->x1 += AOM_INTERP_EXTEND;
            *x_pad = 1;
        }

        if ((scaled_mv.row & SUBPEL_MASK) || (sf->y_step_q4 != SUBPEL_SHIFTS)) {
            block->y0 -= AOM_INTERP_EXTEND;
            block->y1 += AOM_INTERP_EXTEND;
            *y_pad = 1;
        }

        return 1;
    }

    return 0;
}

static INLINE void extend_mc_border(void *src, int32_t *src_stride,
    PadBlock *block, MV32 scaled_mv, const ScaleFactors *sf, int32_t highbd,
    uint8_t *mc_buf, EbDecPicBuf *ref_buf,
    void **src_mod, int32_t ss_x, int32_t ss_y)
    {
    int32_t x_pad = 0, y_pad = 0;
    if (update_extend_mc_border_params(block, ref_buf, scaled_mv, sf, &x_pad,
        &y_pad, ss_x, ss_y))
    {
        int32_t buf_stride = *src_stride;
        *src_mod = NULL;
        const int32_t b_w = block->x1 - block->x0;
        const int32_t b_h = block->y1 - block->y0;

        // Extend the border.
        if (highbd) {
            highbd_build_mc_border((uint8_t*)src, buf_stride, mc_buf, b_w,
                block->x0, block->y0, b_w, b_h,
                ref_buf->superres_upscaled_width >> ss_x,
                ref_buf->frame_height >> ss_y);
        }
        else {
            build_mc_border((uint8_t*)src, buf_stride, mc_buf, b_w, block->x0,
                block->y0, b_w, b_h,
                ref_buf->superres_upscaled_width >> ss_x ,
                ref_buf->frame_height >> ss_y);
        }
        mc_buf = mc_buf + ((y_pad * AOM_INTERP_EXTEND * b_w +
            x_pad * AOM_INTERP_EXTEND) << highbd);
        *src_mod = (void*)mc_buf;
        *src_stride = b_w;
    }
}

void svt_make_inter_predictor(PartitionInfo *part_info, int32_t ref, void *src, int32_t src_stride,
                              void *dst_mod, int32_t dst_stride, EbDecPicBuf *ref_buf,
                              int32_t pre_x, int32_t pre_y, int32_t bw, int32_t bh,
                              ConvolveParams *conv_params, int32_t plane, int32_t do_warp, EbBool is_16bit) {
    const BlockModeInfo *mi         = part_info->mi;
    const int32_t        is_intrabc = is_intrabc_block(mi);

    const int32_t ss_x      = plane ? part_info->subsampling_x : 0;
    const int32_t ss_y      = plane ? part_info->subsampling_y : 0;
    int32_t       bit_depth = ref_buf->ps_pic_buf->bit_depth;
    int32_t       highbd = bit_depth > EB_8BIT || is_16bit;

    /*ScaleFactor*/
    const struct ScaleFactors *const sf =
        is_intrabc ? part_info->sf_identity : part_info->block_ref_sf[ref];

    const MV     mv = mi->mv[ref].as_mv;
    MV           mv_q4;
    void *       src_mod;
    SubpelParams subpel_params;
    do_warp = do_warp && !av1_is_scaled(sf);
    PadBlock block;
    MV32 scaled_mv;

    const int32_t is_scaled = av1_is_scaled(sf);
    if (is_scaled) {
        int orig_pos_y = (pre_y + 0) << SUBPEL_BITS;
        orig_pos_y += mv.row * (1 << (1 - ss_y));
        int orig_pos_x = (pre_x + 0) << SUBPEL_BITS;
        orig_pos_x += mv.col * (1 << (1 - ss_x));
        int pos_y = sf->scale_value_y(orig_pos_y, sf);
        int pos_x = sf->scale_value_x(orig_pos_x, sf);
        pos_x += SCALE_EXTRA_OFF;
        pos_y += SCALE_EXTRA_OFF;

        const int top    = -AOM_LEFT_TOP_MARGIN_SCALED(ss_y);
        const int left   = -AOM_LEFT_TOP_MARGIN_SCALED(ss_x);
        const int bottom = ((ref_buf->frame_height >> ss_y) + AOM_INTERP_EXTEND)
                           << SCALE_SUBPEL_BITS;
        const int right = ((ref_buf->superres_upscaled_width >> ss_x) + AOM_INTERP_EXTEND)
                          << SCALE_SUBPEL_BITS;

        pos_y = clamp(pos_y, top, bottom);
        pos_x = clamp(pos_x, left, right);

        subpel_params.subpel_x = pos_x & SCALE_SUBPEL_MASK;
        subpel_params.subpel_y = pos_y & SCALE_SUBPEL_MASK;
        subpel_params.xs       = sf->x_step_q4;
        subpel_params.ys       = sf->y_step_q4;

        // Get reference block top left coordinate.
        block.x0 = pos_x >> SCALE_SUBPEL_BITS;;
        block.y0 = pos_y >> SCALE_SUBPEL_BITS;;

        // Get reference block bottom right coordinate.
        block.x1 =
            ((pos_x + (bw - 1) * subpel_params.xs) >> SCALE_SUBPEL_BITS) + 1;
        block.y1 =
            ((pos_y + (bh - 1) * subpel_params.ys) >> SCALE_SUBPEL_BITS) + 1;

        MV temp_mv;
        temp_mv = dec_clamp_mv_to_umv_border_sb(part_info->mb_to_left_edge,
                                                part_info->mb_to_right_edge,
                                                part_info->mb_to_top_edge,
                                                part_info->mb_to_bottom_edge,
                                                &mv,
                                                bw,
                                                bh,
                                                ss_x,
                                                ss_y);
        scaled_mv = av1_scale_mv(&temp_mv, (pre_x + 0), (pre_y + 0), sf);
        scaled_mv.row += SCALE_EXTRA_OFF;
        scaled_mv.col += SCALE_EXTRA_OFF;
        int32_t src_offset = (block.y0 * src_stride ) + block.x0;
        src_mod            = (void *)((uint8_t *)src + (src_offset << highbd));
    } else {
        mv_q4 = dec_clamp_mv_to_umv_border_sb(part_info->mb_to_left_edge,
                                              part_info->mb_to_right_edge,
                                              part_info->mb_to_top_edge,
                                              part_info->mb_to_bottom_edge,
                                              &mv,
                                              bw,
                                              bh,
                                              ss_x,
                                              ss_y);


        // Get block position in current frame.
        int pos_x = (pre_x + 0) << SUBPEL_BITS;
        int pos_y = (pre_y + 0) << SUBPEL_BITS;
        // Get reference block top left coordinate.
        pos_x += mv_q4.col;
        pos_y += mv_q4.row;
        block.x0 = pos_x >> SUBPEL_BITS;
        block.y0 = pos_y >> SUBPEL_BITS;

        // Get reference block bottom right coordinate.
        block.x1 = (pos_x >> SUBPEL_BITS) + (bw - 1) + 1;
        block.y1 = (pos_y >> SUBPEL_BITS) + (bh - 1) + 1;
        scaled_mv.row =(int32_t) mv_q4.row;
        scaled_mv.col = (int32_t)mv_q4.col;
        int32_t src_offset = (((pre_y) + (mv_q4.row >> SUBPEL_BITS)) * src_stride) + (pre_x) +
                             (mv_q4.col >> SUBPEL_BITS);
        src_mod = (void *)((uint8_t *)src + (src_offset << highbd));

        subpel_params.xs       = SCALE_SUBPEL_SHIFTS;
        subpel_params.ys       = SCALE_SUBPEL_SHIFTS;
        subpel_params.subpel_x = (mv_q4.col & SUBPEL_MASK) << SCALE_EXTRA_BITS;
        subpel_params.subpel_y = (mv_q4.row & SUBPEL_MASK) << SCALE_EXTRA_BITS;
    }

    if ((!do_warp && !is_intrabc) || (is_scaled && !do_warp && !is_intrabc)) {
        extend_mc_border(src, &src_stride, &block, scaled_mv, sf, highbd,
            part_info->mc_buf[ref], ref_buf, &src_mod, ss_x, ss_y);
    }
    assert(IMPLIES(is_intrabc, !do_warp));

    if (do_warp) {
        const EbWarpedMotionParams *wm_params = &default_warp_params;

        const EbWarpedMotionParams *const wm_global =
            &part_info->ps_global_motion[mi->ref_frame[ref]];
        const EbWarpedMotionParams *const wm_local = &part_info->local_warp_params;

        wm_params = (mi->motion_mode == WARPED_CAUSAL) ? wm_local : wm_global;

        eb_av1_warp_plane((EbWarpedMotionParams *)wm_params,
                          highbd,
                          bit_depth,
                          src,
                          ref_buf->ps_pic_buf->width >> ss_x,
                          ref_buf->ps_pic_buf->height >> ss_y,
                          src_stride,
                          dst_mod,
                          pre_x,
                          pre_y,
                          bw,
                          bh,
                          dst_stride,
                          ss_x,
                          ss_y,
                          conv_params);
    } else if (highbd) {
        uint16_t *src16 = (uint16_t *)src_mod;

        svt_highbd_inter_predictor(src16,
                                   src_stride,
                                   dst_mod,
                                   dst_stride,
                                   &subpel_params,
                                   sf,
                                   bw,
                                   bh,
                                   conv_params,
                                   mi->interp_filters,
                                   is_intrabc,
                                   bit_depth);
    } else {
        svt_inter_predictor(src_mod,
                            src_stride,
                            dst_mod,
                            dst_stride,
                            &subpel_params,
                            sf,
                            bw,
                            bh,
                            conv_params,
                            mi->interp_filters,
                            is_intrabc);
    }
}

void svt_make_masked_inter_predictor(PartitionInfo *part_info, int32_t ref, void *src,
                                     int32_t src_stride, void *dst_ptr, int32_t dst_stride,
                                     EbDecPicBuf *ref_buf, int32_t pre_x, int32_t pre_y, int32_t bw,
                                     int32_t bh, ConvolveParams *conv_params, int32_t plane,
                                     uint8_t *seg_mask, int32_t do_warp, EbBool is_16bit) {
    InterInterCompoundData *comp_data = &part_info->mi->inter_inter_compound;
    const BlockSize         bsize     = part_info->mi->sb_type;
    int32_t                 bit_depth = ref_buf->ps_pic_buf->bit_depth;
    //We come here when we have a prediction done using regular path for the ref0 stored in conv_param.dst.
    //use regular path to generate a prediction for ref1 into  a temporary buffer,
    //then  blend that temporary buffer with that from  the first reference.

#define INTER_PRED_BYTES_PER_PIXEL 2
    DECLARE_ALIGNED(32, uint8_t, tmp_buf[INTER_PRED_BYTES_PER_PIXEL * MAX_SB_SQUARE]);
#undef INTER_PRED_BYTES_PER_PIXEL
    //uint8_t *tmp_dst =  tmp_buf;
    const int tmp_buf_stride = MAX_SB_SIZE;

    CONV_BUF_TYPE *org_dst        = conv_params->dst; //save the ref0 prediction pointer
    int            org_dst_stride = conv_params->dst_stride;
    CONV_BUF_TYPE *tmp_buf16      = (CONV_BUF_TYPE *)tmp_buf;
    conv_params->dst              = tmp_buf16;
    conv_params->dst_stride       = tmp_buf_stride;
    assert(conv_params->do_average == 0);
    assert(conv_params->is_compound == 1);

    svt_make_inter_predictor(part_info,
                             ref,
                             src,
                             src_stride,
                             dst_ptr,
                             dst_stride,
                             ref_buf,
                             pre_x,
                             pre_y,
                             bw,
                             bh,
                             conv_params,
                             plane,
                             do_warp,
                             is_16bit);

    if (!plane && comp_data->type == COMPOUND_DIFFWTD) {
        //CHKN  for DIFF: need to compute the mask  comp_data->seg_mask is
        //the output computed from the two preds org_dst and tmp_buf16
        //for WEDGE the mask is fixed from the table based on wedge_sign/index
        av1_build_compound_diffwtd_mask_d16(seg_mask,
                                            comp_data->mask_type,
                                            org_dst,
                                            org_dst_stride,
                                            tmp_buf16,
                                            tmp_buf_stride,
                                            bh,
                                            bw,
                                            conv_params,
                                            bit_depth);
    }

    build_masked_compound_no_round((uint8_t *)dst_ptr,
                                   dst_stride,
                                   org_dst,
                                   org_dst_stride,
                                   tmp_buf16,
                                   tmp_buf_stride,
                                   comp_data,
                                   seg_mask,
                                   bsize,
                                   bh,
                                   bw,
                                   conv_params,
                                   (uint8_t)bit_depth,
                                   is_16bit);
}

void av1_combine_interintra(PartitionInfo *part_info, BlockSize bsize, int plane,
                            uint8_t *inter_pred, int inter_stride, uint8_t *intra_pred,
                            int intra_stride, EbBitDepthEnum bit_depth, EbBool is_16bit) {
    BlockModeInfo * mi          = part_info->mi;
    int32_t         sub_x       = (plane > 0) ? part_info->subsampling_x : 0;
    int32_t         sub_y       = (plane > 0) ? part_info->subsampling_y : 0;
    const BlockSize plane_bsize = get_plane_block_size(bsize, sub_x, sub_y);

    if (bit_depth > EB_8BIT || is_16bit) {
        /*As per spec we r considering interitra_wedge_sign is always "zero"*/
        /*Check buffers, Aom  2nd time inter_pred buffer plane is plane independent */
        combine_interintra_highbd(mi->interintra_mode_params.interintra_mode,
                                  mi->interintra_mode_params.wedge_interintra,
                                  mi->interintra_mode_params.interintra_wedge_index,
                                  0 /*interintra_wedgesign*/,
                                  bsize,
                                  plane_bsize,
                                  inter_pred,
                                  inter_stride,
                                  inter_pred,
                                  inter_stride,
                                  intra_pred,
                                  intra_stride,
                                  bit_depth);
        return;
    }

    /*Check buffers, Aom  2nd time inter_pred buffer plane is plane independent */
    combine_interintra(mi->interintra_mode_params.interintra_mode,
                       mi->interintra_mode_params.wedge_interintra,
                       mi->interintra_mode_params.interintra_wedge_index,
                       0 /*interintra_wedgesign*/,
                       bsize,
                       plane_bsize,
                       inter_pred,
                       inter_stride,
                       inter_pred,
                       inter_stride,
                       intra_pred,
                       intra_stride);
}

void av1_build_intra_predictors_for_interintra(DecModCtxt *dec_mod_ctxt, PartitionInfo *part_info,
                                               void *pv_blk_recon_buf, int32_t recon_stride,
                                               BlockSize bsize, int32_t plane, uint8_t *dst,
                                               int dst_stride, EbBitDepthEnum bit_depth) {
    EbDecHandle *dec_handle = (EbDecHandle *)dec_mod_ctxt->dec_handle_ptr;
    EbBool is16b = dec_handle->is_16bit_pipeline;
    BlockModeInfo *mi          = part_info->mi;
    int32_t        sub_x       = (plane > 0) ? part_info->subsampling_x : 0;
    int32_t        sub_y       = (plane > 0) ? part_info->subsampling_y : 0;
    BlockSize      plane_bsize = get_plane_block_size(bsize, sub_x, sub_y);
    PredictionMode mode = interintra_to_intra_mode[mi->interintra_mode_params.interintra_mode];
    assert(mi->angle_delta[PLANE_TYPE_Y] == 0);
    assert(mi->angle_delta[PLANE_TYPE_UV] == 0);
    assert(mi->filter_intra_mode_info.use_filter_intra == 0);
    assert(mi->use_intrabc == 0);
    assert(mi->palette_size[plane != 0] == 0);

    void *pv_top_neighbor_array, *pv_left_neighbor_array;

    if (bit_depth == EB_8BIT && !is16b) {
        EbByte buf = (EbByte)pv_blk_recon_buf;

        pv_top_neighbor_array  = (void *)(buf - recon_stride);
        pv_left_neighbor_array = (void *)(buf - 1);
    } else { //16bit
        uint16_t *buf          = (uint16_t *)pv_blk_recon_buf;
        pv_top_neighbor_array  = (void *)(buf - recon_stride);
        pv_left_neighbor_array = (void *)(buf - 1);
    }

    /*Calling Intra prediction */
    svtav1_predict_intra_block(part_info,
                               plane,
                               max_txsize_rect_lookup[plane_bsize],
                               &dec_mod_ctxt->cur_tile_info,
                               (void *)dst,
                               dst_stride,
                               pv_top_neighbor_array,
                               pv_left_neighbor_array,
                               recon_stride,
                               dec_mod_ctxt->seq_header,
                               mode,
                               0,
                               0,
                               bit_depth,
                               is16b);
}

/* Build interintra_predictors */
void av1_build_interintra_predictors(DecModCtxt *dec_mod_ctxt, PartitionInfo *part_info, void *pred,
                                     int32_t stride, int plane, BlockSize bsize,
                                     EbBitDepthEnum bit_depth, EbBool is_16bit) {
    if (bit_depth > EB_8BIT || is_16bit) {
        DECLARE_ALIGNED(16, uint16_t, intrapredictor[MAX_SB_SQUARE]);
        av1_build_intra_predictors_for_interintra(dec_mod_ctxt,
                                                  part_info,
                                                  pred,
                                                  stride,
                                                  bsize,
                                                  plane,
                                                  (uint8_t *)intrapredictor,
                                                  MAX_SB_SIZE,
                                                  bit_depth);
        av1_combine_interintra(part_info,
                               bsize,
                               plane,
                               pred,
                               stride,
                               (uint8_t *)intrapredictor,
                               MAX_SB_SIZE,
                               bit_depth,
                               is_16bit);
    } else {
        DECLARE_ALIGNED(16, uint8_t, intrapredictor[MAX_SB_SQUARE]);
        av1_build_intra_predictors_for_interintra(dec_mod_ctxt,
                                                  part_info,
                                                  pred,
                                                  stride,
                                                  bsize,
                                                  plane,
                                                  intrapredictor,
                                                  MAX_SB_SIZE,
                                                  bit_depth);
        av1_combine_interintra(
            part_info, bsize, plane, pred, stride, intrapredictor, MAX_SB_SIZE, bit_depth, is_16bit);
    }
}

void svtav1_predict_inter_block_plane(DecModCtxt *dec_mod_ctx, EbDecHandle *dec_hdl,
                                      PartitionInfo *part_info, int32_t plane,
                                      int32_t build_for_obmc, int32_t mi_x, int32_t mi_y, void *dst,
                                      int32_t dst_stride, int32_t some_use_intra,
                                      int32_t bit_depth) {
    const BlockModeInfo *mi          = part_info->mi;
    const FrameHeader *  cur_frm_hdr = dec_mod_ctx->frame_header;
    SeqHeader *          seq_header  = dec_mod_ctx->seq_header;
    int32_t              is_compound = has_second_ref(mi);
    int32_t              ref;
    const int32_t        is_intrabc = is_intrabc_block(mi);

    //temporary buffer for joint compound, move this to context if stack does not hold.
    DECLARE_ALIGNED(32, uint16_t, tmp_dst[128 * 128]);

    EbBool is16b = (bit_depth > EB_8BIT) || dec_hdl->is_16bit_pipeline;

    const BlockSize bsize     = mi->sb_type;
    assert(bsize < BlockSizeS_ALL);
    const int32_t   ss_x      = plane ? part_info->subsampling_x : 0;
    const int32_t   ss_y      = plane ? part_info->subsampling_y : 0;
    int32_t         bw        = part_info->wpx[0] >> ss_x;
    int32_t         bh        = part_info->hpx[0] >> ss_y;
    int32_t         row_start = 0;
    int32_t         col_start = 0;
    // For sub8x8 chroma blocks, we may be covering more than one luma block's
    // worth of pixels. Thus (mi_x, mi_y) may not be the correct coordinates for
    // the top-left corner of the prediction source - the correct top-left corner
    // is at (pre_x, pre_y).
    if (some_use_intra && plane && !build_for_obmc) {
        bool sub8_w = (block_size_wide[bsize] == 4) && ss_x;
        bool sub8_h = (block_size_high[bsize] == 4) && ss_y;

        if (sub8_w) {
            if (part_info->mi_col & 0x1) {
                bw        = bw << 1;
                col_start = -1;
            } else
                return;
        }

        if (sub8_h) {
            if (part_info->mi_row & 0x1) {
                bh        = bh << 1;
                row_start = -1;
            } else
                return;
        }
    }
    const int32_t pre_x = (mi_x + MI_SIZE * col_start) >> ss_x;
    const int32_t pre_y = (mi_y + MI_SIZE * row_start) >> ss_y;

    int32_t dst_offset =
        ((MI_SIZE * col_start) >> ss_x) + ((MI_SIZE * row_start * dst_stride) >> ss_y);
    void *dst_mod = (void *)((uint8_t *)dst + (dst_offset << is16b));

    assert(IMPLIES(is_intrabc, !is_compound));
    {
        ConvolveParams conv_params =
            get_conv_params_no_round(0, 0, plane, tmp_dst, MAX_SB_SIZE, is_compound, bit_depth);

        int bck_frame_index = 0, fwd_frame_index = 0;
        int cur_frame_index = cur_frm_hdr->order_hint;

        EbDecPicBuf *bck_buf = get_ref_frame_buf(dec_hdl, mi->ref_frame[0]);
        EbDecPicBuf *fwd_buf = get_ref_frame_buf(dec_hdl, mi->ref_frame[1]);

        if (bck_buf != NULL) bck_frame_index = bck_buf->order_hint;
        if (fwd_buf != NULL) fwd_frame_index = fwd_buf->order_hint;

        /*Distantance WTD compound inter prediction */
        av1_dist_wtd_comp_weight_assign(seq_header,
                                        cur_frame_index,
                                        bck_frame_index,
                                        fwd_frame_index,
                                        (int)mi->compound_idx,
                                        0,
                                        &conv_params.fwd_offset,
                                        &conv_params.bck_offset,
                                        &conv_params.use_dist_wtd_comp_avg,
                                        is_compound);
        conv_params.use_jnt_comp_avg = conv_params.use_dist_wtd_comp_avg;

        for (ref = 0; ref < 1 + is_compound; ++ref) {
            const int32_t                     mode = mi->mode;
            const EbWarpedMotionParams *const wm_global =
                &part_info->ps_global_motion[mi->ref_frame[ref]];

            EbDecPicBuf *ref_buf = is_intrabc ? dec_hdl->cur_pic_buf[0]
                                              : get_ref_frame_buf(dec_hdl, mi->ref_frame[ref]);
            EbPictureBufferDesc *ps_ref_pic_buf = ref_buf->ps_pic_buf;

            int32_t do_warp =
                (bw >= 8 && bh >= 8 && !build_for_obmc && (cur_frm_hdr->force_integer_mv == 0) &&
                 (((mode == GLOBALMV || mode == GLOBAL_GLOBALMV) &&
                   (wm_global->wmtype > TRANSLATION)) ||
                  (mi->motion_mode == WARPED_CAUSAL)));

            void *  src;
            int32_t src_stride;

            derive_blk_pointers(ps_ref_pic_buf, plane, 0, 0, &src, &src_stride, ss_x, ss_y);

            conv_params.do_average = ref;
            /*support masked inter prediction based on WEDGE / DIFFWTD compound type */
            if (is_masked_compound_type(mi->inter_inter_compound.type)) {
                // masked compound type has its own average mechanism
                conv_params.do_average = 0;
            }

            if (ref && is_masked_compound_type(mi->inter_inter_compound.type)) {
                svt_make_masked_inter_predictor(part_info,
                                                ref,
                                                src,
                                                src_stride,
                                                dst_mod,
                                                dst_stride,
                                                ref_buf,
                                                pre_x,
                                                pre_y,
                                                bw,
                                                bh,
                                                &conv_params,
                                                plane,
                                                dec_mod_ctx->seg_mask,
                                                do_warp,
                                                is16b);
            } else {
                svt_make_inter_predictor(part_info,
                                         ref,
                                         src,
                                         src_stride,
                                         dst_mod,
                                         dst_stride,
                                         ref_buf,
                                         pre_x,
                                         pre_y,
                                         bw,
                                         bh,
                                         &conv_params,
                                         plane,
                                         do_warp,
                                         is16b);
            }
        }
    }
}

void svtav1_predict_inter_block(DecModCtxt *dec_mod_ctxt, EbDecHandle *dec_hdl,
                                PartitionInfo *part_info, int32_t mi_row, int32_t mi_col,
                                int32_t num_planes) {
    void *  blk_recon_buf;
    int32_t recon_stride;
    int32_t sub_x, sub_y;
    int32_t some_use_intra;

    EbPictureBufferDesc *recon_picture_buf = dec_hdl->cur_pic_buf[0]->ps_pic_buf;

    /* scan through sub 8 blocks and see if anyof them is intra */
    some_use_intra        = 0;
    const BlockSize bsize = part_info->mi->sb_type;
    sub_x                 = part_info->subsampling_x;
    sub_y                 = part_info->subsampling_y;
    bool sub8_w           = (block_size_wide[bsize] == 4) && sub_x;
    bool sub8_h           = (block_size_high[bsize] == 4) && sub_y;
    if (sub8_h || sub8_w) {
        int32_t i, j;
        /* Floor and Ceil to nearest 8x8 blks */
        const int row_start = sub8_h ? (mi_row & (~1)) : mi_row;
        const int row_end   = sub8_h ? (mi_row | (1)) : mi_row;
        const int col_start = sub8_w ? (mi_col & (~1)) : mi_col;
        const int col_end   = sub8_w ? (mi_col | (1)) : mi_col;

        for (i = row_start; i <= row_end; i++) {
            for (j = col_start; j <= col_end; j++) {
                BlockModeInfo *mode_info = get_cur_mode_info(dec_hdl, i, j, part_info->sb_info);
                if (mode_info->ref_frame[0] == INTRA_FRAME) some_use_intra = 1;
            }
        }
    }

    for (int plane = 0; plane < num_planes; ++plane) {
        sub_x = (plane > 0) ? part_info->subsampling_x : 0;
        sub_y = (plane > 0) ? part_info->subsampling_y : 0;

        derive_blk_pointers(recon_picture_buf,
                            plane,
                            mi_col * MI_SIZE >> sub_x,
                            mi_row * MI_SIZE >> sub_y,
                            &blk_recon_buf,
                            &recon_stride,
                            sub_x,
                            sub_y);

        svtav1_predict_inter_block_plane(dec_mod_ctxt,
                                         dec_hdl,
                                         part_info,
                                         plane,
                                         0,
                                         mi_col * MI_SIZE,
                                         mi_row * MI_SIZE,
                                         blk_recon_buf,
                                         recon_stride,
                                         some_use_intra,
                                         recon_picture_buf->bit_depth);

        if (is_interintra_pred(part_info->mi)) {
            /*Inter prd is done in above function, In the below function Intra prd happens follwed by interintra blending */
            av1_build_interintra_predictors(dec_mod_ctxt,
                                            part_info,
                                            blk_recon_buf,
                                            recon_stride,
                                            plane,
                                            bsize,
                                            recon_picture_buf->bit_depth,
                                            dec_hdl->is_16bit_pipeline);
        }
    }
    if (part_info->mi->motion_mode == OBMC_CAUSAL) {
        dec_build_obmc_inter_predictors_sb(
            (void *)dec_mod_ctxt, dec_hdl, part_info, mi_row, mi_col);
    }

    return;
}
