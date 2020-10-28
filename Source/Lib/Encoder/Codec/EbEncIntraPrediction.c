// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>
#include <string.h>
#include "EbEncIntraPrediction.h"
#include "EbModeDecisionProcess.h"
#include "common_dsp_rtcd.h"

static PartitionType from_shape_to_part[] = {
        PARTITION_NONE,
        PARTITION_HORZ,
        PARTITION_VERT,
        PARTITION_HORZ_A,
        PARTITION_HORZ_B,
        PARTITION_VERT_A,
        PARTITION_VERT_B,
        PARTITION_HORZ_4,
        PARTITION_VERT_4,
        PARTITION_SPLIT
};


static int get_filt_type(const MacroBlockD *xd, int plane) {
    int ab_sm, le_sm;

    if (plane == 0) {
        const MbModeInfo *ab = xd->above_mbmi;
        const MbModeInfo *le = xd->left_mbmi;
        ab_sm = ab ? is_smooth(&ab->block_mi, plane) : 0;
        le_sm = le ? is_smooth(&le->block_mi, plane) : 0;
    }
    else {
        const MbModeInfo *ab = xd->chroma_above_mbmi;
        const MbModeInfo *le = xd->chroma_left_mbmi;
        ab_sm = ab ? is_smooth(&ab->block_mi, plane) : 0;
        le_sm = le ? is_smooth(&le->block_mi, plane) : 0;
    }

    return (ab_sm || le_sm) ? 1 : 0;
}

////////////#####################...........Recurssive intra prediction ending...........#####################////////////

static void build_intra_predictors(
        const MacroBlockD *xd,
        uint8_t* top_neigh_array,
        uint8_t* left_neigh_array,
        // const uint8_t *ref,    int32_t ref_stride,
        uint8_t *dst, int32_t dst_stride,
        PredictionMode mode, int32_t angle_delta,
        FilterIntraMode filter_intra_mode,
        TxSize tx_size, int32_t disable_edge_filter,
        int32_t n_top_px, int32_t n_topright_px,
        int32_t n_left_px, int32_t n_bottomleft_px,
        int32_t plane)
{
    (void)xd;
    int32_t i;

    int32_t ref_stride = 1;
    const uint8_t *above_ref = top_neigh_array;//CHKN ref - ref_stride;
    const uint8_t *left_ref = left_neigh_array;//CHKN ref - 1;
    DECLARE_ALIGNED(32, uint8_t, left_data[MAX_TX_SIZE * 2 + 48]);
    DECLARE_ALIGNED(32, uint8_t, above_data[MAX_TX_SIZE * 2 + 48]);
    uint8_t *const above_row = above_data + 32;
    uint8_t *const left_col = left_data + 32;

    const int32_t txwpx = tx_size_wide[tx_size];
    const int32_t txhpx = tx_size_high[tx_size];
    int32_t need_left = extend_modes[mode] & NEED_LEFT;
    int32_t need_above = extend_modes[mode] & NEED_ABOVE;
    int32_t need_above_left = extend_modes[mode] & NEED_ABOVELEFT;
    int32_t p_angle = 0;
    const int32_t is_dr_mode = av1_is_directional_mode(mode);
    const int32_t use_filter_intra = filter_intra_mode != FILTER_INTRA_MODES;

    if (is_dr_mode) {
        p_angle = mode_to_angle_map[mode] + angle_delta * ANGLE_STEP;
        if (p_angle <= 90)
            need_above = 1, need_left = 0, need_above_left = 1;
        else if (p_angle < 180)
            need_above = 1, need_left = 1, need_above_left = 1;
        else
            need_above = 0, need_left = 1, need_above_left = 1;
    }
    if (use_filter_intra) need_left = need_above = need_above_left = 1;

    assert(n_top_px >= 0);
    assert(n_topright_px >= 0);
    assert(n_left_px >= 0);
    assert(n_bottomleft_px >= 0);

    if ((!need_above && n_left_px == 0) || (!need_left && n_top_px == 0)) {
        int32_t val;
        if (need_left)
            val = (n_top_px > 0) ? above_ref[0] : 129;
        else
            val = (n_left_px > 0) ? left_ref[0] : 127;
        for (i = 0; i < txhpx; ++i) {
            memset(dst, val, txwpx);
            dst += dst_stride;
        }
        return;
    }

    // NEED_LEFT
    if (need_left) {
        int32_t need_bottom = !!(extend_modes[mode] & NEED_BOTTOMLEFT);
        if (use_filter_intra) need_bottom = 0;
        if (is_dr_mode) need_bottom = p_angle > 180;
        const int32_t num_left_pixels_needed = txhpx + (need_bottom ? txwpx : 0);
        i = 0;
        if (n_left_px > 0) {
            for (; i < n_left_px; i++) left_col[i] = left_ref[i * ref_stride];
            if (need_bottom && n_bottomleft_px > 0) {
                assert(i == txhpx);
                for (; i < txhpx + n_bottomleft_px; i++)
                    left_col[i] = left_ref[i * ref_stride];
            }
            if (i < num_left_pixels_needed)
                memset(&left_col[i], left_col[i - 1], num_left_pixels_needed - i);
        }
        else {
            if (n_top_px > 0)
                memset(left_col, above_ref[0], num_left_pixels_needed);
            else
                memset(left_col, 129, num_left_pixels_needed);
        }
    }

    // NEED_ABOVE
    if (need_above) {
        int32_t need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
        if (use_filter_intra) need_right = 0;
        if (is_dr_mode) need_right = p_angle < 90;
        const int32_t num_top_pixels_needed = txwpx + (need_right ? txhpx : 0);
        if (n_top_px > 0) {
            svt_memcpy(above_row, above_ref, n_top_px);
            i = n_top_px;
            if (need_right && n_topright_px > 0) {
                assert(n_top_px == txwpx);
                svt_memcpy(above_row + txwpx, above_ref + txwpx, n_topright_px);
                i += n_topright_px;
            }
            if (i < num_top_pixels_needed)
                memset(&above_row[i], above_row[i - 1], num_top_pixels_needed - i);
        }
        else {
            if (n_left_px > 0)
                memset(above_row, left_ref[0], num_top_pixels_needed);
            else
                memset(above_row, 127, num_top_pixels_needed);
        }
    }

    if (need_above_left) {
        if (n_top_px > 0 && n_left_px > 0)
            above_row[-1] = above_ref[-1];
        else if (n_top_px > 0)
            above_row[-1] = above_ref[0];
        else if (n_left_px > 0)
            above_row[-1] = left_ref[0];
        else
            above_row[-1] = 128;
        left_col[-1] = above_row[-1];
    }
    if (use_filter_intra) {
        svt_av1_filter_intra_predictor(dst, dst_stride, tx_size, above_row, left_col,
                                       filter_intra_mode);
        return;
    }

    if (is_dr_mode) {
        int32_t upsample_above = 0;
        int32_t upsample_left = 0;
        if (!disable_edge_filter) {
            const int32_t need_right = p_angle < 90;
            const int32_t need_bottom = p_angle > 180;
            const int32_t filt_type = get_filt_type(xd, plane);

            if (p_angle != 90 && p_angle != 180) {
                const int32_t ab_le = need_above_left ? 1 : 0;
                if (need_above && need_left && (txwpx + txhpx >= 24))
                    filter_intra_edge_corner(above_row, left_col);
                if (need_above && n_top_px > 0) {
                    const int32_t strength =
                            intra_edge_filter_strength(txwpx, txhpx, p_angle - 90, filt_type);
                    const int32_t n_px = n_top_px + ab_le + (need_right ? txhpx : 0);
                    svt_av1_filter_intra_edge(above_row - ab_le, n_px, strength);
                }
                if (need_left && n_left_px > 0) {
                    const int32_t strength = intra_edge_filter_strength(
                            txhpx, txwpx, p_angle - 180, filt_type);
                    const int32_t n_px = n_left_px + ab_le + (need_bottom ? txwpx : 0);
                    svt_av1_filter_intra_edge(left_col - ab_le, n_px, strength);
                }
            }
            upsample_above =
                    use_intra_edge_upsample(txwpx, txhpx, p_angle - 90, filt_type);
            if (need_above && upsample_above) {
                const int32_t n_px = txwpx + (need_right ? txhpx : 0);
                svt_av1_upsample_intra_edge(above_row, n_px);
            }
            upsample_left =
                    use_intra_edge_upsample(txhpx, txwpx, p_angle - 180, filt_type);
            if (need_left && upsample_left) {
                const int32_t n_px = txhpx + (need_bottom ? txwpx : 0);
                svt_av1_upsample_intra_edge(left_col, n_px);
            }
        }
        dr_predictor(dst, dst_stride, tx_size, above_row, left_col, upsample_above,
                     upsample_left, p_angle);
        return;
    }

    // predict
    if (mode == DC_PRED) {
        dc_pred[n_left_px > 0][n_top_px > 0][tx_size](dst, dst_stride, above_row,
                                                      left_col);
    }
    else
        eb_pred[mode][tx_size](dst, dst_stride, above_row, left_col);
}
static void build_intra_predictors_high(
        const MacroBlockD *xd,
        uint16_t* top_neigh_array, // int8_t
        uint16_t* left_neigh_array, // int8_t
        //const uint8_t *ref8, int32_t ref_stride,
        uint16_t *dst,//uint8_t *dst8
        int32_t dst_stride, PredictionMode mode, int32_t angle_delta,
        FilterIntraMode filter_intra_mode, TxSize tx_size,
        int32_t disable_edge_filter, int32_t n_top_px, int32_t n_topright_px, int32_t n_left_px,
        int32_t n_bottomleft_px, int32_t plane, int32_t bd) {
    (void)xd;
    int32_t i;
    //uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
    //uint16_t *ref = CONVERT_TO_SHORTPTR(ref8);

    DECLARE_ALIGNED(16, uint16_t, left_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(16, uint16_t, above_data[MAX_TX_SIZE * 2 + 32]);
    uint16_t *const above_row = above_data + 16;
    uint16_t *const left_col = left_data + 16;
    const int32_t txwpx = tx_size_wide[tx_size];
    const int32_t txhpx = tx_size_high[tx_size];
    int32_t need_left = extend_modes[mode] & NEED_LEFT;
    int32_t need_above = extend_modes[mode] & NEED_ABOVE;
    int32_t need_above_left = extend_modes[mode] & NEED_ABOVELEFT;

    int32_t ref_stride = 1;
    const uint16_t *above_ref = top_neigh_array;
    const uint16_t *left_ref = left_neigh_array;
    //const uint16_t *above_ref = ref - ref_stride;
    //const uint16_t *left_ref = ref - 1;
    int32_t p_angle = 0;
    const int32_t is_dr_mode = av1_is_directional_mode(mode);
    const int32_t use_filter_intra = filter_intra_mode != FILTER_INTRA_MODES;
    int32_t base = 128 << (bd - 8);

    // The default values if ref pixels are not available:
    // base-1 base-1 base-1 .. base-1 base-1 base-1 base-1 base-1 base-1
    // base+1   A      b  ..     Y      Z
    // base+1   C      D  ..     W      X
    // base+1   E      F  ..     U      V
    // base+1   G      H  ..     S      T      T      T      T      T

    if (is_dr_mode) {
        p_angle = mode_to_angle_map[mode] + angle_delta * ANGLE_STEP;
        if (p_angle <= 90)
            need_above = 1, need_left = 0, need_above_left = 1;
        else if (p_angle < 180)
            need_above = 1, need_left = 1, need_above_left = 1;
        else
            need_above = 0, need_left = 1, need_above_left = 1;
    }
    if (use_filter_intra) need_left = need_above = need_above_left = 1;

    assert(n_top_px >= 0);
    assert(n_topright_px >= 0);
    assert(n_left_px >= 0);
    assert(n_bottomleft_px >= 0);

    if ((!need_above && n_left_px == 0) || (!need_left && n_top_px == 0)) {
        int32_t val;
        if (need_left)
            val = (n_top_px > 0) ? above_ref[0] : base + 1;
        else
            val = (n_left_px > 0) ? left_ref[0] : base - 1;
        for (i = 0; i < txhpx; ++i) {
            svt_aom_memset16(dst, val, txwpx);
            dst += dst_stride;
        }
        return;
    }

    // NEED_LEFT
    if (need_left) {
        int32_t need_bottom = !!(extend_modes[mode] & NEED_BOTTOMLEFT);
        if (use_filter_intra) need_bottom = 0;
        if (is_dr_mode) need_bottom = p_angle > 180;
        const int32_t num_left_pixels_needed = txhpx + (need_bottom ? txwpx : 0);
        i = 0;
        if (n_left_px > 0) {
            for (; i < n_left_px; i++) left_col[i] = left_ref[i * ref_stride];
            if (need_bottom && n_bottomleft_px > 0) {
                assert(i == txhpx);
                for (; i < txhpx + n_bottomleft_px; i++)
                    left_col[i] = left_ref[i * ref_stride];
            }
            if (i < num_left_pixels_needed)
                svt_aom_memset16(&left_col[i], left_col[i - 1], num_left_pixels_needed - i);
        }
        else {
            if (n_top_px > 0)
                svt_aom_memset16(left_col, above_ref[0], num_left_pixels_needed);
            else
                svt_aom_memset16(left_col, base + 1, num_left_pixels_needed);
        }
    }

    // NEED_ABOVE
    if (need_above) {
        int32_t need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
        if (use_filter_intra) need_right = 0;
        if (is_dr_mode) need_right = p_angle < 90;
        const int32_t num_top_pixels_needed = txwpx + (need_right ? txhpx : 0);
        if (n_top_px > 0) {
            svt_memcpy(above_row, above_ref, n_top_px * sizeof(above_ref[0]));
            i = n_top_px;
            if (need_right && n_topright_px > 0) {
                assert(n_top_px == txwpx);
                svt_memcpy(above_row + txwpx, above_ref + txwpx,
                       n_topright_px * sizeof(above_ref[0]));
                i += n_topright_px;
            }
            if (i < num_top_pixels_needed)
                svt_aom_memset16(&above_row[i], above_row[i - 1],
                                num_top_pixels_needed - i);
        }
        else {
            if (n_left_px > 0)
                svt_aom_memset16(above_row, left_ref[0], num_top_pixels_needed);
            else
                svt_aom_memset16(above_row, base - 1, num_top_pixels_needed);
        }
    }

    if (need_above_left) {
        if (n_top_px > 0 && n_left_px > 0)
            above_row[-1] = above_ref[-1];
        else if (n_top_px > 0)
            above_row[-1] = above_ref[0];
        else if (n_left_px > 0)
            above_row[-1] = left_ref[0];
        else
            above_row[-1] = (uint16_t)base;
        left_col[-1] = above_row[-1];
    }
    if (use_filter_intra) {
        highbd_filter_intra_predictor(dst, dst_stride, tx_size, above_row, left_col,
                                      filter_intra_mode, bd);
        return;
    }
    if (is_dr_mode) {
        int32_t upsample_above = 0;
        int32_t upsample_left = 0;
        if (!disable_edge_filter) {
            const int32_t need_right = p_angle < 90;
            const int32_t need_bottom = p_angle > 180;
            //const int32_t filt_type = get_filt_type(xd, plane);
            const int32_t filt_type = get_filt_type(xd, plane);
            if (p_angle != 90 && p_angle != 180) {
                const int32_t ab_le = need_above_left ? 1 : 0;
                if (need_above && need_left && (txwpx + txhpx >= 24))
                    filter_intra_edge_corner_high(above_row, left_col);
                if (need_above && n_top_px > 0) {
                    const int32_t strength =
                            intra_edge_filter_strength(txwpx, txhpx, p_angle - 90, filt_type);
                    const int32_t n_px = n_top_px + ab_le + (need_right ? txhpx : 0);
                    svt_av1_filter_intra_edge_high(above_row - ab_le, n_px, strength);
                }
                if (need_left && n_left_px > 0) {
                    const int32_t strength = intra_edge_filter_strength(
                            txhpx, txwpx, p_angle - 180, filt_type);
                    const int32_t n_px = n_left_px + ab_le + (need_bottom ? txwpx : 0);

                    svt_av1_filter_intra_edge_high(left_col - ab_le, n_px, strength);
                }
            }
            upsample_above =
                    use_intra_edge_upsample(txwpx, txhpx, p_angle - 90, filt_type);
            if (need_above && upsample_above) {
                const int32_t n_px = txwpx + (need_right ? txhpx : 0);
                //av1_upsample_intra_edge_high(above_row, n_px, bd);// AMIR : to be replaced by optimized code
                svt_av1_upsample_intra_edge_high_c(above_row, n_px, bd);
            }
            upsample_left =
                    use_intra_edge_upsample(txhpx, txwpx, p_angle - 180, filt_type);
            if (need_left && upsample_left) {
                const int32_t n_px = txhpx + (need_bottom ? txwpx : 0);
                //av1_upsample_intra_edge_high(left_col, n_px, bd);// AMIR: to be replaced by optimized code
                svt_av1_upsample_intra_edge_high_c(left_col, n_px, bd);
            }
        }
        highbd_dr_predictor(dst, dst_stride, tx_size, above_row, left_col,
                            upsample_above, upsample_left, p_angle, bd);
        return;
    }

    // predict
    if (mode == DC_PRED) {
        dc_pred_high[n_left_px > 0][n_top_px > 0][tx_size](
                dst, dst_stride, above_row, left_col, bd);
    }
    else
        pred_high[mode][tx_size](dst, dst_stride, above_row, left_col, bd);
}


void svt_av1_predict_intra_block(
        TileInfo * tile,
        STAGE       stage,
        const BlockGeom * blk_geom,
        const Av1Common *cm,
        int32_t wpx,
        int32_t hpx,
        TxSize tx_size,
        PredictionMode mode,
        int32_t angle_delta,
        int32_t use_palette,
        PaletteInfo  *palette_info,
        FilterIntraMode filter_intra_mode,
        uint8_t* top_neigh_array,
        uint8_t* left_neigh_array,
        EbPictureBufferDesc  *recon_buffer,
        int32_t col_off,
        int32_t row_off,
        int32_t plane,
        BlockSize bsize,
        uint32_t txb_org_x_pict,
        uint32_t txb_org_y_pict,
        uint32_t bl_org_x_pict,
        uint32_t bl_org_y_pict,
        uint32_t bl_org_x_mb,
        uint32_t bl_org_y_mb,
        ModeInfo **mi_grid_base,
        SeqHeader *seq_header_ptr)
{
    MacroBlockD xd_s;
    MacroBlockD *xd = &xd_s;

    xd_s.chroma_left_mbmi = 0;
    xd_s.chroma_above_mbmi = 0;

    uint32_t  pred_buf_x_offest;
    uint32_t  pred_buf_y_offest;

    if (stage == ED_STAGE) { // EncDec
        pred_buf_x_offest = plane ? ((bl_org_x_pict >> 3) << 3) >> 1 : txb_org_x_pict;
        pred_buf_y_offest = plane ? ((bl_org_y_pict >> 3) << 3) >> 1 : txb_org_y_pict;
    }
    else { // MD
        pred_buf_x_offest = bl_org_x_mb;
        pred_buf_y_offest = bl_org_y_mb;
    }

    // Adjust mirow , micol ;
    // All plane have the same values

    int32_t mirow = bl_org_y_pict >> 2;
    int32_t micol = bl_org_x_pict >> 2;
    xd->up_available   = (mirow > tile->mi_row_start);
    xd->left_available = (micol > tile->mi_col_start);
    const int32_t bw = mi_size_wide[bsize];
    const int32_t bh = mi_size_high[bsize];

    xd->mb_to_top_edge = -((mirow * MI_SIZE) * 8);
    xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mirow) * MI_SIZE) * 8;
    xd->mb_to_left_edge = -((micol * MI_SIZE) * 8);
    xd->mb_to_right_edge = ((cm->mi_cols - bw - micol) * MI_SIZE) * 8;
    xd->tile.mi_col_start = tile->mi_col_start;
    xd->tile.mi_col_end = tile->mi_col_end;
    xd->tile.mi_row_start = tile->mi_row_start;
    xd->tile.mi_row_end = tile->mi_row_end;
    xd->n8_h = bh;
    xd->n8_w = bw;
    xd->is_sec_rect = 0;
    if (xd->n8_w < xd->n8_h) {
        // Only mark is_sec_rect as 1 for the last block.
        // For PARTITION_VERT_4, it would be (0, 0, 0, 1);
        // For other partitions, it would be (0, 1).
        if (!((micol + xd->n8_w) & (xd->n8_h - 1))) xd->is_sec_rect = 1;
    }

    if (xd->n8_w > xd->n8_h)
        if (mirow & (xd->n8_w - 1)) xd->is_sec_rect = 1;
    uint8_t  *dst;
    int32_t dst_stride;
    if (plane == 0) {
        dst = recon_buffer->buffer_y + pred_buf_x_offest + recon_buffer->origin_x + (pred_buf_y_offest + recon_buffer->origin_y)*recon_buffer->stride_y;
        dst_stride = recon_buffer->stride_y;
    }
    else if (plane == 1) {
        dst = recon_buffer->buffer_cb + (pred_buf_x_offest + recon_buffer->origin_x / 2 + (pred_buf_y_offest + recon_buffer->origin_y / 2)*recon_buffer->stride_cb);
        dst_stride = recon_buffer->stride_cb;
    }
    else {
        dst = recon_buffer->buffer_cr + (pred_buf_x_offest + recon_buffer->origin_x / 2 + (pred_buf_y_offest + recon_buffer->origin_y / 2)*recon_buffer->stride_cr);
        dst_stride = recon_buffer->stride_cr;
    }

    int32_t chroma_up_available = xd->up_available;
    int32_t chroma_left_available = xd->left_available;
    const int32_t ss_x = plane == 0 ? 0 : 1; //CHKN
    const int32_t ss_y = plane == 0 ? 0 : 1;

    if (ss_x && bw < mi_size_wide[BLOCK_8X8])
        chroma_left_available = (micol - 1) > tile->mi_col_start;
    if (ss_y && bh < mi_size_high[BLOCK_8X8])
        chroma_up_available = (mirow - 1) > tile->mi_row_start;

    int mi_stride = cm->mi_stride;
    const int32_t offset = mirow * mi_stride + micol;
    xd->mi = mi_grid_base + offset;
    ModeInfo *mi_ptr = *xd->mi;

    if (xd->up_available) {
        // xd->above_mbmi = xd->mi[-xd->mi_stride].mbmi;
        xd->above_mbmi = &mi_ptr[-mi_stride].mbmi;
    }
    else
        xd->above_mbmi = NULL;
    if (xd->left_available) {
        //xd->left_mbmi = xd->mi[-1].mbmi;
        xd->left_mbmi = &mi_ptr[-1].mbmi;
    }
    else
        xd->left_mbmi = NULL;
    const int chroma_ref = ((mirow & 0x01) || !(bh & 0x01) || !ss_y) &&
                           ((micol & 0x01) || !(bw & 0x01) || !ss_x);
    if (chroma_ref) {
        // To help calculate the "above" and "left" chroma blocks, note that the
        // current block may cover multiple luma blocks (eg, if partitioned into
        // 4x4 luma blocks).
        // First, find the top-left-most luma block covered by this chroma block

        mi_ptr = xd->mi[-(mirow & ss_y) * mi_stride - (micol & ss_x)];

        // Then, we consider the luma region covered by the left or above 4x4 chroma
        // prediction. We want to point to the chroma reference block in that
        // region, which is the bottom-right-most mi unit.
        // This leads to the following offsets:
        xd->chroma_above_mbmi = chroma_up_available ? &mi_ptr[-mi_stride + ss_x].mbmi : NULL;

        xd->chroma_left_mbmi = chroma_left_available ? &mi_ptr[ss_y * mi_stride - 1].mbmi : NULL;
    }

    //CHKN  const MbModeInfo *const mbmi = xd->mi[0];
    const int32_t txwpx = tx_size_wide[tx_size];
    const int32_t txhpx = tx_size_high[tx_size];
    const int32_t x = col_off << tx_size_wide_log2[0];
    const int32_t y = row_off << tx_size_high_log2[0];
    if (use_palette) {
        const uint8_t *const map = palette_info->color_idx_map;
        const uint16_t *const palette =
                palette_info->pmi.palette_colors + plane * PALETTE_MAX_SIZE;
        for (int32_t r = 0; r < txhpx; ++r)
            for (int32_t c = 0; c < txwpx; ++c)
                dst[r * dst_stride + c] =
                        (uint8_t)palette[map[(r + y) * wpx + c + x]];
        return;
    }

    //CHKN BlockSize bsize = mbmi->sb_type;
    struct MacroblockdPlane  pd_s;
    struct MacroblockdPlane * pd = &pd_s;
    if (plane == 0)
        pd->subsampling_x = pd->subsampling_y = 0;
    else
        pd->subsampling_x = pd->subsampling_y = 1;
    const int32_t txw = tx_size_wide_unit[tx_size];
    const int32_t txh = tx_size_high_unit[tx_size];
    const int32_t have_top = row_off || (pd->subsampling_y ? /*xd->*/chroma_up_available
                                                           : xd->up_available);
    const int32_t have_left =
            col_off ||
            (pd->subsampling_x ? /*xd->*/chroma_left_available : xd->left_available);
    const int32_t mi_row = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
    const int32_t mi_col = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);
    const int32_t xr_chr_offset = 0;
    const int32_t yd_chr_offset = 0;

    // Distance between the right edge of this prediction block to
    // the frame right edge
    const int32_t xr = (xd->mb_to_right_edge >> (3 + pd->subsampling_x)) +
                       (wpx - x - txwpx) - xr_chr_offset;
    // Distance between the bottom edge of this prediction block to
    // the frame bottom edge
    const int32_t yd = (xd->mb_to_bottom_edge >> (3 + pd->subsampling_y)) +
                       (hpx - y - txhpx) - yd_chr_offset;
    const int32_t right_available =
            mi_col + ((col_off + txw) << pd->subsampling_x) < xd->tile.mi_col_end;
    const int32_t bottom_available =
            (yd > 0) &&
            (mi_row + ((row_off + txh) << pd->subsampling_y) < xd->tile.mi_row_end);

    const PartitionType partition = from_shape_to_part[blk_geom->shape]; //blk_ptr->part;// PARTITION_NONE;//CHKN this is good enough as the avail functions need to know if VERT part is used or not mbmi->partition;

    // force 4x4 chroma component block size.
    bsize = scale_chroma_bsize(bsize, pd->subsampling_x, pd->subsampling_y);

    const int32_t have_top_right = intra_has_top_right(
            seq_header_ptr->sb_size, bsize,
            mi_row, mi_col, have_top, right_available, partition, tx_size,
            row_off, col_off, pd->subsampling_x, pd->subsampling_y);
    const int32_t have_bottom_left = intra_has_bottom_left(
            seq_header_ptr->sb_size, bsize,
            mi_row, mi_col, bottom_available, have_left, partition,
            tx_size, row_off, col_off, pd->subsampling_x, pd->subsampling_y);

    const int32_t disable_edge_filter = !(seq_header_ptr->enable_intra_edge_filter);

    //if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    //  build_intra_predictors_high(
    //      xd, ref, ref_stride, dst, dst_stride, mode, angle_delta,
    //      filter_intra_mode, tx_size, disable_edge_filter,
    //      have_top ? AOMMIN(txwpx, xr + txwpx) : 0,
    //      have_top_right ? AOMMIN(txwpx, xr) : 0,
    //      have_left ? AOMMIN(txhpx, yd + txhpx) : 0,
    //      have_bottom_left ? AOMMIN(txhpx, yd) : 0, plane);
    //  return;
    //}

    build_intra_predictors(
            xd,
            top_neigh_array,
            left_neigh_array,
            // ref, ref_stride,
            dst, dst_stride, mode,
            angle_delta, filter_intra_mode, tx_size,
            disable_edge_filter,
            have_top ? AOMMIN(txwpx, xr + txwpx) : 0,
            have_top_right ? AOMMIN(txwpx, xr) : 0,
            have_left ? AOMMIN(txhpx, yd + txhpx) : 0,
            have_bottom_left ? AOMMIN(txhpx, yd) : 0, plane);
}

void svt_av1_predict_intra_block_16bit(
        EbBitDepthEnum bit_depth,
        TileInfo * tile,
        STAGE       stage,
        const BlockGeom * blk_geom,
        const Av1Common *cm,
        int32_t wpx,
        int32_t hpx,
        TxSize tx_size,
        PredictionMode mode,
        int32_t angle_delta,
        int32_t use_palette,
        PaletteInfo  *palette_info,
        FilterIntraMode filter_intra_mode,
        uint16_t* top_neigh_array,
        uint16_t* left_neigh_array,
        EbPictureBufferDesc  *recon_buffer,
        int32_t col_off,
        int32_t row_off,
        int32_t plane,
        BlockSize bsize,
        uint32_t txb_org_x_pict,
        uint32_t txb_org_y_pict,
        uint32_t bl_org_x_pict,
        uint32_t bl_org_y_pict,
        uint32_t bl_org_x_mb,
        uint32_t bl_org_y_mb,
        ModeInfo **mi_grid_base,
        SeqHeader *seq_header_ptr)
{
    MacroBlockD xd_s;
    MacroBlockD *xd = &xd_s;

    xd_s.chroma_left_mbmi = 0;
    xd_s.chroma_above_mbmi = 0;

    uint32_t  pred_buf_x_offest;
    uint32_t  pred_buf_y_offest;

    if (stage == ED_STAGE) { // EncDec
        pred_buf_x_offest = plane ? ((bl_org_x_pict >> 3) << 3) >> 1 : txb_org_x_pict;
        pred_buf_y_offest = plane ? ((bl_org_y_pict >> 3) << 3) >> 1 : txb_org_y_pict;
    } else { // MD
        pred_buf_x_offest = bl_org_x_mb;
        pred_buf_y_offest = bl_org_y_mb;
    }

    int32_t mirow = bl_org_y_pict >> 2;
    int32_t micol = bl_org_x_pict >> 2;
    xd->up_available = (mirow > tile->mi_row_start);
    xd->left_available = (micol > tile->mi_col_start);

    const int32_t bw = mi_size_wide[bsize];
    const int32_t bh = mi_size_high[bsize];

    xd->mb_to_top_edge = -((mirow * MI_SIZE) * 8);
    xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mirow) * MI_SIZE) * 8;
    xd->mb_to_left_edge = -((micol * MI_SIZE) * 8);
    xd->mb_to_right_edge = ((cm->mi_cols - bw - micol) * MI_SIZE) * 8;
    xd->tile.mi_col_start = tile->mi_col_start;
    xd->tile.mi_col_end = tile->mi_col_end;
    xd->tile.mi_row_start = tile->mi_row_start;
    xd->tile.mi_row_end = tile->mi_row_end;
    xd->n8_h = bh;
    xd->n8_w = bw;
    xd->is_sec_rect = 0;
    if (xd->n8_w < xd->n8_h) {
        // Only mark is_sec_rect as 1 for the last block.
        // For PARTITION_VERT_4, it would be (0, 0, 0, 1);
        // For other partitions, it would be (0, 1).
        if (!((micol + xd->n8_w) & (xd->n8_h - 1))) xd->is_sec_rect = 1;
    }

    if (xd->n8_w > xd->n8_h)
        if (mirow & (xd->n8_w - 1)) xd->is_sec_rect = 1;

    // Adjust prediction pointers
    uint16_t *dst;
    int32_t dst_stride;
    if (plane == 0) {
        dst = (uint16_t*)(recon_buffer->buffer_y) + pred_buf_x_offest + recon_buffer->origin_x + (pred_buf_y_offest + recon_buffer->origin_y)*recon_buffer->stride_y;
        dst_stride = recon_buffer->stride_y;
    }
    else if (plane == 1) {
        dst = (uint16_t*)(recon_buffer->buffer_cb) + (pred_buf_x_offest + recon_buffer->origin_x / 2 + (pred_buf_y_offest + recon_buffer->origin_y / 2)*recon_buffer->stride_cb);
        dst_stride = recon_buffer->stride_cb;
    }
    else {
        dst = (uint16_t*)(recon_buffer->buffer_cr) + (pred_buf_x_offest + recon_buffer->origin_x / 2 + (pred_buf_y_offest + recon_buffer->origin_y / 2)*recon_buffer->stride_cr);
        dst_stride = recon_buffer->stride_cr;
    }

    int32_t chroma_up_available = xd->up_available;
    int32_t chroma_left_available = xd->left_available;

    const int32_t ss_x = plane == 0 ? 0 : 1;
    const int32_t ss_y = plane == 0 ? 0 : 1;

    if (ss_x && bw < mi_size_wide[BLOCK_8X8])
        chroma_left_available = (micol - 1) > tile->mi_col_start;
    if (ss_y && bh < mi_size_high[BLOCK_8X8])
        chroma_up_available = (mirow - 1) > tile->mi_row_start;

    int mi_stride = cm->mi_stride;
    const int32_t offset = mirow * mi_stride + micol;
    xd->mi = mi_grid_base + offset;
    ModeInfo *mi_ptr = *xd->mi;

    if (xd->up_available) {
        // xd->above_mbmi = xd->mi[-xd->mi_stride].mbmi;
        xd->above_mbmi = &mi_ptr[-mi_stride].mbmi;
    }
    else
        xd->above_mbmi = NULL;
    if (xd->left_available) {
        //xd->left_mbmi = xd->mi[-1].mbmi;
        xd->left_mbmi = &mi_ptr[-1].mbmi;
    }
    else
        xd->left_mbmi = NULL;
    const int chroma_ref = ((mirow & 0x01) || !(bh & 0x01) || !ss_y) &&
                           ((micol & 0x01) || !(bw & 0x01) || !ss_x);
    if (chroma_ref) {
        // To help calculate the "above" and "left" chroma blocks, note that the
        // current block may cover multiple luma blocks (eg, if partitioned into
        // 4x4 luma blocks).
        // First, find the top-left-most luma block covered by this chroma block

        mi_ptr = xd->mi[-(mirow & ss_y) * mi_stride - (micol & ss_x)];

        // Then, we consider the luma region covered by the left or above 4x4 chroma
        // prediction. We want to point to the chroma reference block in that
        // region, which is the bottom-right-most mi unit.
        // This leads to the following offsets:
        xd->chroma_above_mbmi = chroma_up_available ? &mi_ptr[-mi_stride + ss_x].mbmi : NULL;
        xd->chroma_left_mbmi = chroma_left_available ? &mi_ptr[ss_y * mi_stride - 1].mbmi : NULL;
    }

    //CHKN  const MbModeInfo *const mbmi = xd->mi[0];
    const int32_t txwpx = tx_size_wide[tx_size];
    const int32_t txhpx = tx_size_high[tx_size];
    const int32_t x = col_off << tx_size_wide_log2[0];
    const int32_t y = row_off << tx_size_high_log2[0];
    if (use_palette) {
        const uint8_t *const map = palette_info->color_idx_map;
        const uint16_t *const palette =
            palette_info->pmi.palette_colors + plane * PALETTE_MAX_SIZE;
        uint16_t              max_val = (bit_depth == EB_8BIT) ? 0xFF : 0xFFFF;
        for (int32_t r = 0; r < txhpx; ++r)
            for (int32_t c = 0; c < txwpx; ++c)
                dst[r * dst_stride + c] = palette[map[(r + y) * wpx + c + x]] > max_val ? max_val : palette[map[(r + y) * wpx + c + x]];
        return;
    }

    //CHKN BlockSize bsize = mbmi->sb_type;

    struct MacroblockdPlane  pd_s;
    struct MacroblockdPlane * pd = &pd_s;
    if (plane == 0)
        pd->subsampling_x = pd->subsampling_y = 0;
    else
        pd->subsampling_x = pd->subsampling_y = 1;
    const int32_t txw = tx_size_wide_unit[tx_size];
    const int32_t txh = tx_size_high_unit[tx_size];
    const int32_t have_top = row_off || (pd->subsampling_y ? /*xd->*/chroma_up_available
                                                           : xd->up_available);
    const int32_t have_left =
            col_off ||
            (pd->subsampling_x ? /*xd->*/chroma_left_available : xd->left_available);
    const int32_t mi_row = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
    const int32_t mi_col = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);
    const int32_t xr_chr_offset = 0;
    const int32_t yd_chr_offset = 0;

    // Distance between the right edge of this prediction block to
    // the frame right edge
    const int32_t xr = (xd->mb_to_right_edge >> (3 + pd->subsampling_x)) +
                       (wpx - x - txwpx) - xr_chr_offset;
    // Distance between the bottom edge of this prediction block to
    // the frame bottom edge
    const int32_t yd = (xd->mb_to_bottom_edge >> (3 + pd->subsampling_y)) +
                       (hpx - y - txhpx) - yd_chr_offset;
    const int32_t right_available =
            mi_col + ((col_off + txw) << pd->subsampling_x) < xd->tile.mi_col_end;
    const int32_t bottom_available =
            (yd > 0) &&
            (mi_row + ((row_off + txh) << pd->subsampling_y) < xd->tile.mi_row_end);

    const PartitionType partition = from_shape_to_part[blk_geom->shape]; //blk_ptr->part;// PARTITION_NONE;//CHKN this is good enough as the avail functions need to know if VERT part is used or not mbmi->partition;

    // force 4x4 chroma component block size.
    bsize = scale_chroma_bsize(bsize, pd->subsampling_x, pd->subsampling_y);

    const int32_t have_top_right = intra_has_top_right(
            seq_header_ptr->sb_size, bsize,
            mi_row, mi_col, have_top, right_available, partition, tx_size,
            row_off, col_off, pd->subsampling_x, pd->subsampling_y);
    const int32_t have_bottom_left = intra_has_bottom_left(
            seq_header_ptr->sb_size, bsize,
            mi_row, mi_col, bottom_available, have_left, partition,
            tx_size, row_off, col_off, pd->subsampling_x, pd->subsampling_y);

    const int32_t disable_edge_filter = !(seq_header_ptr->enable_intra_edge_filter);

    build_intra_predictors_high(
            xd,
            top_neigh_array,
            left_neigh_array,
            // ref, ref_stride,
            dst, dst_stride, mode,
            angle_delta, filter_intra_mode, tx_size,
            disable_edge_filter,
            have_top ? AOMMIN(txwpx, xr + txwpx) : 0,
            have_top_right ? AOMMIN(txwpx, xr) : 0,
            have_left ? AOMMIN(txhpx, yd + txhpx) : 0,
        have_bottom_left ? AOMMIN(txhpx, yd) : 0, plane, bit_depth);
}

/** IntraPrediction()
is the main function to compute intra prediction for a PU
*/
EbErrorType svt_av1_intra_prediction_cl(
        uint8_t                              hbd_mode_decision,
        ModeDecisionContext                  *md_context_ptr,
        PictureControlSet                    *pcs_ptr,
        ModeDecisionCandidateBuffer           *candidate_buffer_ptr)
{
    (void) hbd_mode_decision;
    EbErrorType return_error = EB_ErrorNone;

    uint32_t mode_type_left_neighbor_index = get_neighbor_array_unit_left_index(
            md_context_ptr->mode_type_neighbor_array,
            md_context_ptr->blk_origin_y);
    uint32_t mode_type_top_neighbor_index = get_neighbor_array_unit_top_index(
            md_context_ptr->mode_type_neighbor_array,
            md_context_ptr->blk_origin_x);
    uint32_t intra_luma_mode_left_neighbor_index = get_neighbor_array_unit_left_index(
            md_context_ptr->intra_luma_mode_neighbor_array,
            md_context_ptr->blk_origin_y);
    uint32_t intra_luma_mode_top_neighbor_index = get_neighbor_array_unit_top_index(
            md_context_ptr->intra_luma_mode_neighbor_array,
            md_context_ptr->blk_origin_x);

    uint32_t intra_chroma_mode_left_neighbor_index = get_neighbor_array_unit_left_index(
            md_context_ptr->intra_chroma_mode_neighbor_array,
            md_context_ptr->round_origin_y >> 1);
    uint32_t intra_chroma_mode_top_neighbor_index = get_neighbor_array_unit_top_index(
            md_context_ptr->intra_chroma_mode_neighbor_array,
            md_context_ptr->round_origin_x >> 1);

    md_context_ptr->intra_luma_left_mode = (uint32_t)(
            (md_context_ptr->mode_type_neighbor_array->left_array[mode_type_left_neighbor_index] != INTRA_MODE) ? DC_PRED/*EB_INTRA_DC*/ :
            (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->left_array[intra_luma_mode_left_neighbor_index]);

    md_context_ptr->intra_luma_top_mode = (uint32_t)(
            (md_context_ptr->mode_type_neighbor_array->top_array[mode_type_top_neighbor_index] != INTRA_MODE) ? DC_PRED/*EB_INTRA_DC*/ :
            (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->top_array[intra_luma_mode_top_neighbor_index]);       //   use DC. This seems like we could use a SB-width

    md_context_ptr->intra_chroma_left_mode = (uint32_t)(
            (md_context_ptr->mode_type_neighbor_array->left_array[mode_type_left_neighbor_index] != INTRA_MODE) ? UV_DC_PRED :
            (uint32_t)md_context_ptr->intra_chroma_mode_neighbor_array->left_array[intra_chroma_mode_left_neighbor_index]);

    md_context_ptr->intra_chroma_top_mode = (uint32_t)(
            (md_context_ptr->mode_type_neighbor_array->top_array[mode_type_top_neighbor_index] != INTRA_MODE) ? UV_DC_PRED :
            (uint32_t)md_context_ptr->intra_chroma_mode_neighbor_array->top_array[intra_chroma_mode_top_neighbor_index]);       //   use DC. This seems like we could use a SB-width
    TxSize  tx_size = md_context_ptr->blk_geom->txsize[candidate_buffer_ptr->candidate_ptr->tx_depth][0]; // Nader - Intra 128x128 not supported
    TxSize  tx_size_chroma = md_context_ptr->blk_geom->txsize_uv[candidate_buffer_ptr->candidate_ptr->tx_depth][0]; //Nader - Intra 128x128 not supported

    if(!md_context_ptr->hbd_mode_decision) {
        uint8_t    top_neigh_array[64 * 2 + 1];
        uint8_t    left_neigh_array[64 * 2 + 1];
        PredictionMode mode;
        // Hsan: plane should be derived @ an earlier stage (e.g. @ the call of perform_fast_loop())
        int32_t start_plane = (md_context_ptr->uv_intra_comp_only) ? 1 : 0;
        int32_t end_plane = (md_context_ptr->blk_geom->has_uv && md_context_ptr->chroma_level <= CHROMA_MODE_1 && !md_context_ptr->md_staging_skip_chroma_pred) ? (int)MAX_MB_PLANE : 1;
        for (int32_t plane = start_plane; plane < end_plane; ++plane) {
            if (plane == 0) {
                if (md_context_ptr->blk_origin_y != 0)
                    svt_memcpy(top_neigh_array + 1, md_context_ptr->luma_recon_neighbor_array->top_array + md_context_ptr->blk_origin_x, md_context_ptr->blk_geom->bwidth * 2);
                if (md_context_ptr->blk_origin_x != 0)
                    svt_memcpy(left_neigh_array + 1, md_context_ptr->luma_recon_neighbor_array->left_array + md_context_ptr->blk_origin_y, md_context_ptr->blk_geom->bheight * 2);
                if (md_context_ptr->blk_origin_y != 0 && md_context_ptr->blk_origin_x != 0)
                    top_neigh_array[0] = left_neigh_array[0] = md_context_ptr->luma_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE + md_context_ptr->blk_origin_x - md_context_ptr->blk_origin_y];
            }

            else if (plane == 1) {
                if (md_context_ptr->round_origin_y != 0)
                    svt_memcpy(top_neigh_array + 1, md_context_ptr->cb_recon_neighbor_array->top_array + md_context_ptr->round_origin_x / 2, md_context_ptr->blk_geom->bwidth_uv * 2);

                if (md_context_ptr->round_origin_x != 0)
                    svt_memcpy(left_neigh_array + 1, md_context_ptr->cb_recon_neighbor_array->left_array + md_context_ptr->round_origin_y / 2, md_context_ptr->blk_geom->bheight_uv * 2);

                if (md_context_ptr->round_origin_y != 0 && md_context_ptr->round_origin_x != 0)
                    top_neigh_array[0] = left_neigh_array[0] = md_context_ptr->cb_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE / 2 + md_context_ptr->round_origin_x / 2 - md_context_ptr->round_origin_y / 2];
            }
            else {
                if (md_context_ptr->round_origin_y != 0)
                    svt_memcpy(top_neigh_array + 1, md_context_ptr->cr_recon_neighbor_array->top_array + md_context_ptr->round_origin_x / 2, md_context_ptr->blk_geom->bwidth_uv * 2);

                if (md_context_ptr->round_origin_x != 0)
                    svt_memcpy(left_neigh_array + 1, md_context_ptr->cr_recon_neighbor_array->left_array + md_context_ptr->round_origin_y / 2, md_context_ptr->blk_geom->bheight_uv * 2);

                if (md_context_ptr->round_origin_y != 0 && md_context_ptr->round_origin_x != 0)
                    top_neigh_array[0] = left_neigh_array[0] = md_context_ptr->cr_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE / 2 + md_context_ptr->round_origin_x / 2 - md_context_ptr->round_origin_y / 2];
            }

            if (plane)
                mode = (candidate_buffer_ptr->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode) UV_DC_PRED : (PredictionMode) candidate_buffer_ptr->candidate_ptr->intra_chroma_mode;
            else
                mode = candidate_buffer_ptr->candidate_ptr->pred_mode;

            svt_av1_predict_intra_block(
                    &md_context_ptr->sb_ptr->tile_info,
                    !ED_STAGE,
                    md_context_ptr->blk_geom,
                    pcs_ptr->parent_pcs_ptr->av1_cm,                                      //const Av1Common *cm,
                    plane ? md_context_ptr->blk_geom->bwidth_uv : md_context_ptr->blk_geom->bwidth,          //int32_t wpx,
                    plane ? md_context_ptr->blk_geom->bheight_uv : md_context_ptr->blk_geom->bheight,          //int32_t hpx,
                    plane ? tx_size_chroma : tx_size,                                               //TxSize tx_size,
                    mode,                                                                           //PredictionMode mode,
                    plane ? candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_UV] : candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
                    plane==0 ? (candidate_buffer_ptr->candidate_ptr->palette_info ?
                                    candidate_buffer_ptr->candidate_ptr->palette_info->pmi.palette_size[0]>0 : 0) : 0,
                    plane==0 ? candidate_buffer_ptr->candidate_ptr->palette_info : NULL,    //MD
                    plane ? FILTER_INTRA_MODES : candidate_buffer_ptr->candidate_ptr->filter_intra_mode,
                    top_neigh_array + 1,
                    left_neigh_array + 1,
                    candidate_buffer_ptr->prediction_ptr,                                              //uint8_t *dst,
                    //int32_t dst_stride,
                    0,                                                                              //int32_t col_off,
                    0,                                                                              //int32_t row_off,
                    plane,                                                                          //int32_t plane,
                    md_context_ptr->blk_geom->bsize,       //uint32_t puSize,
                    md_context_ptr->blk_origin_x,
                    md_context_ptr->blk_origin_y,
                    md_context_ptr->blk_origin_x,                  //uint32_t cuOrgX,
                    md_context_ptr->blk_origin_y,                  //uint32_t cuOrgY
                    plane ? ((md_context_ptr->blk_geom->origin_x >> 3) << 3) / 2 : md_context_ptr->blk_geom->origin_x,  //uint32_t cuOrgX used only for prediction Ptr
                    plane ? ((md_context_ptr->blk_geom->origin_y >> 3) << 3) / 2 : md_context_ptr->blk_geom->origin_y,   //uint32_t cuOrgY used only for prediction Ptr
                    pcs_ptr->mi_grid_base,
                    &((SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr)->seq_header
            );
        }
    } else {
        uint16_t    top_neigh_array[64 * 2 + 1];
        uint16_t    left_neigh_array[64 * 2 + 1];
        PredictionMode mode;
        // Hsan: plane should be derived @ an earlier stage (e.g. @ the call of perform_fast_loop())
        int32_t start_plane = (md_context_ptr->uv_intra_comp_only) ? 1 : 0;
        int32_t end_plane = (md_context_ptr->blk_geom->has_uv && md_context_ptr->chroma_level <= CHROMA_MODE_1 && !md_context_ptr->md_staging_skip_chroma_pred) ? (int)MAX_MB_PLANE : 1;
        for (int32_t plane = start_plane; plane < end_plane; ++plane) {
            if (plane == 0) {
                if (md_context_ptr->blk_origin_y != 0)
                    svt_memcpy(top_neigh_array + 1, (uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->top_array) + md_context_ptr->blk_origin_x, md_context_ptr->blk_geom->bwidth * 2 * sizeof(uint16_t));

                if (md_context_ptr->blk_origin_x != 0)
                    svt_memcpy(left_neigh_array + 1, (uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->left_array) + md_context_ptr->blk_origin_y, md_context_ptr->blk_geom->bheight * 2 * sizeof(uint16_t));

                if (md_context_ptr->blk_origin_y != 0 && md_context_ptr->blk_origin_x != 0)
                    top_neigh_array[0] = left_neigh_array[0] = ((uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->top_left_array) + MAX_PICTURE_HEIGHT_SIZE + md_context_ptr->blk_origin_x - md_context_ptr->blk_origin_y)[0];
            }
            else if (plane == 1) {
                if (md_context_ptr->round_origin_y != 0)
                    svt_memcpy(top_neigh_array + 1, (uint16_t*)(md_context_ptr->cb_recon_neighbor_array16bit->top_array) + md_context_ptr->round_origin_x / 2, md_context_ptr->blk_geom->bwidth_uv * 2 * sizeof(uint16_t));

                if (md_context_ptr->round_origin_x != 0)
                    svt_memcpy(left_neigh_array + 1, (uint16_t*)(md_context_ptr->cb_recon_neighbor_array16bit->left_array) + md_context_ptr->round_origin_y / 2, md_context_ptr->blk_geom->bheight_uv * 2 * sizeof(uint16_t));

                if (md_context_ptr->round_origin_y != 0 && md_context_ptr->round_origin_x != 0)
                    top_neigh_array[0] = left_neigh_array[0] = ((uint16_t*) (md_context_ptr->cb_recon_neighbor_array16bit->top_left_array) + MAX_PICTURE_HEIGHT_SIZE / 2 + md_context_ptr->round_origin_x / 2 - md_context_ptr->round_origin_y / 2)[0];
            }
            else {
                if (md_context_ptr->round_origin_y != 0)
                    svt_memcpy(top_neigh_array + 1, (uint16_t*)(md_context_ptr->cr_recon_neighbor_array16bit->top_array) + md_context_ptr->round_origin_x / 2, md_context_ptr->blk_geom->bwidth_uv * 2 * sizeof(uint16_t));

                if (md_context_ptr->round_origin_x != 0)
                    svt_memcpy(left_neigh_array + 1, (uint16_t*)(md_context_ptr->cr_recon_neighbor_array16bit->left_array) + md_context_ptr->round_origin_y / 2, md_context_ptr->blk_geom->bheight_uv * 2 * sizeof(uint16_t));

                if (md_context_ptr->round_origin_y != 0 && md_context_ptr->round_origin_x != 0)
                    top_neigh_array[0] = left_neigh_array[0] = ((uint16_t*) (md_context_ptr->cr_recon_neighbor_array16bit->top_left_array) + MAX_PICTURE_HEIGHT_SIZE / 2 + md_context_ptr->round_origin_x / 2 - md_context_ptr->round_origin_y / 2)[0];
            }

            if (plane)
                mode = (candidate_buffer_ptr->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode) UV_DC_PRED : (PredictionMode) candidate_buffer_ptr->candidate_ptr->intra_chroma_mode;
            else
                mode = candidate_buffer_ptr->candidate_ptr->pred_mode;

            svt_av1_predict_intra_block_16bit(
                    EB_10BIT,
                    &md_context_ptr->sb_ptr->tile_info,
                    !ED_STAGE,
                    md_context_ptr->blk_geom,
                    pcs_ptr->parent_pcs_ptr->av1_cm,                                      //const Av1Common *cm,
                    plane ? md_context_ptr->blk_geom->bwidth_uv : md_context_ptr->blk_geom->bwidth,          //int32_t wpx,
                    plane ? md_context_ptr->blk_geom->bheight_uv : md_context_ptr->blk_geom->bheight,          //int32_t hpx,
                    plane ? tx_size_chroma : tx_size,                                               //TxSize tx_size,
                    mode,                                                                           //PredictionMode mode,
                    plane ? candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_UV] : candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
                    plane==0 ? (candidate_buffer_ptr->candidate_ptr->palette_info ?
                                    candidate_buffer_ptr->candidate_ptr->palette_info->pmi.palette_size[0]>0 : 0) : 0,
                    plane==0 ? candidate_buffer_ptr->candidate_ptr->palette_info : NULL,    //MD
                    plane ? FILTER_INTRA_MODES : candidate_buffer_ptr->candidate_ptr->filter_intra_mode,
                    top_neigh_array + 1,
                    left_neigh_array + 1,
                    candidate_buffer_ptr->prediction_ptr,                                              //uint8_t *dst,
                    0,                                                                              //int32_t col_off,
                    0,                                                                              //int32_t row_off,
                    plane,                                                                          //int32_t plane,
                    md_context_ptr->blk_geom->bsize,       //uint32_t puSize,
                    md_context_ptr->blk_origin_x,
                    md_context_ptr->blk_origin_y,
                    md_context_ptr->blk_origin_x,                  //uint32_t cuOrgX,
                    md_context_ptr->blk_origin_y,                  //uint32_t cuOrgY
                    plane ? ((md_context_ptr->blk_geom->origin_x >> 3) << 3) / 2 : md_context_ptr->blk_geom->origin_x,  //uint32_t cuOrgX used only for prediction Ptr
                    plane ? ((md_context_ptr->blk_geom->origin_y >> 3) << 3) / 2 : md_context_ptr->blk_geom->origin_y,   //uint32_t cuOrgY used only for prediction Ptr
                    pcs_ptr->mi_grid_base,
                    &((SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr)->seq_header
            );
        }
    }

    return return_error;
}

EbErrorType  intra_luma_prediction_for_interintra(
        ModeDecisionContext         *md_context_ptr,
        PictureControlSet           *pcs_ptr,
        InterIntraMode              interintra_mode,
        EbPictureBufferDesc         *prediction_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    uint8_t is_inter = 0; // set to 0 b/c this is an intra path

    uint32_t mode_type_left_neighbor_index = get_neighbor_array_unit_left_index(
            md_context_ptr->mode_type_neighbor_array,
            md_context_ptr->blk_origin_y);
    uint32_t mode_type_top_neighbor_index = get_neighbor_array_unit_top_index(
            md_context_ptr->mode_type_neighbor_array,
            md_context_ptr->blk_origin_x);
    uint32_t intra_luma_mode_left_neighbor_index = get_neighbor_array_unit_left_index(
            md_context_ptr->intra_luma_mode_neighbor_array,
            md_context_ptr->blk_origin_y);
    uint32_t intra_luma_mode_top_neighbor_index = get_neighbor_array_unit_top_index(
            md_context_ptr->intra_luma_mode_neighbor_array,
            md_context_ptr->blk_origin_x);

    md_context_ptr->intra_luma_left_mode = (uint32_t)(
            (md_context_ptr->mode_type_neighbor_array->left_array[mode_type_left_neighbor_index] != INTRA_MODE) ? DC_PRED:
            (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->left_array[intra_luma_mode_left_neighbor_index]);

    md_context_ptr->intra_luma_top_mode = (uint32_t)(
            (md_context_ptr->mode_type_neighbor_array->top_array[mode_type_top_neighbor_index] != INTRA_MODE) ? DC_PRED:
            (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->top_array[intra_luma_mode_top_neighbor_index]);       //   use DC. This seems like we could use a SB-width

    TxSize  tx_size = md_context_ptr->blk_geom->txsize[0][0];  //CHKN  TOcheck
    PredictionMode mode = interintra_to_intra_mode[interintra_mode];

    if (!md_context_ptr->hbd_mode_decision) {
        uint8_t    top_neigh_array[64 * 2 + 1];
        uint8_t    left_neigh_array[64 * 2 + 1];

        if (md_context_ptr->blk_origin_y != 0)
            svt_memcpy(top_neigh_array + 1, md_context_ptr->luma_recon_neighbor_array->top_array + md_context_ptr->blk_origin_x, md_context_ptr->blk_geom->bwidth * 2);
        if (md_context_ptr->blk_origin_x != 0)
            svt_memcpy(left_neigh_array + 1, md_context_ptr->luma_recon_neighbor_array->left_array + md_context_ptr->blk_origin_y, md_context_ptr->blk_geom->bheight * 2);
        if (md_context_ptr->blk_origin_y != 0 && md_context_ptr->blk_origin_x != 0)
            top_neigh_array[0] = left_neigh_array[0] = md_context_ptr->luma_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE + md_context_ptr->blk_origin_x - md_context_ptr->blk_origin_y];

        svt_av1_predict_intra_block(
                &md_context_ptr->sb_ptr->tile_info,
                !ED_STAGE,
                md_context_ptr->blk_geom,
                pcs_ptr->parent_pcs_ptr->av1_cm,        //const Av1Common *cm,
                md_context_ptr->blk_geom->bwidth,                       //int32_t wpx,
                md_context_ptr->blk_geom->bheight,                      //int32_t hpx,
                tx_size,                                                //TxSize tx_size,
                mode,                                                   //PredictionMode mode,
                0,                                                      //candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
                0,                                                      //int32_t use_palette,
                NULL,  //Inter-Intra
                FILTER_INTRA_MODES,                                     //CHKN FilterIntraMode filter_intra_mode,
                top_neigh_array + 1,
                left_neigh_array + 1,
                prediction_ptr,                                         //uint8_t *dst,
                (md_context_ptr->blk_geom->tx_org_x[is_inter][0][0] - md_context_ptr->blk_geom->origin_x) >> 2,
                (md_context_ptr->blk_geom->tx_org_y[is_inter][0][0] - md_context_ptr->blk_geom->origin_y) >> 2,
                PLANE_TYPE_Y,                                           //int32_t plane,
                md_context_ptr->blk_geom->bsize,                        //uint32_t puSize,
                md_context_ptr->blk_origin_x,
                md_context_ptr->blk_origin_y,
                md_context_ptr->blk_origin_x,                            //uint32_t cuOrgX,
                md_context_ptr->blk_origin_y,                            //uint32_t cuOrgY
                0,                                                      //cuOrgX used only for prediction Ptr
                0,                                                       //cuOrgY used only for prediction Ptr
                pcs_ptr->mi_grid_base,
                &((SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr)->seq_header
        );
    } else {
        uint16_t top_neigh_array[64 * 2 + 1];
        uint16_t left_neigh_array[64 * 2 + 1];

        if (md_context_ptr->blk_origin_y != 0)
            svt_memcpy(top_neigh_array + 1, (uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->top_array) + md_context_ptr->blk_origin_x, md_context_ptr->blk_geom->bwidth * 2 * sizeof(uint16_t));
        if (md_context_ptr->blk_origin_x != 0)
            svt_memcpy(left_neigh_array + 1, (uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->left_array) + md_context_ptr->blk_origin_y, md_context_ptr->blk_geom->bheight * 2 * sizeof(uint16_t));
        if (md_context_ptr->blk_origin_y != 0 && md_context_ptr->blk_origin_x != 0)
            top_neigh_array[0] = left_neigh_array[0] = ((uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->top_left_array) + MAX_PICTURE_HEIGHT_SIZE + md_context_ptr->blk_origin_x - md_context_ptr->blk_origin_y)[0];

        svt_av1_predict_intra_block_16bit(
                EB_10BIT,
                &md_context_ptr->sb_ptr->tile_info,
                !ED_STAGE,
                md_context_ptr->blk_geom,
                pcs_ptr->parent_pcs_ptr->av1_cm,        //const Av1Common *cm,
                md_context_ptr->blk_geom->bwidth,                       //int32_t wpx,
                md_context_ptr->blk_geom->bheight,                      //int32_t hpx,
                tx_size,                                                //TxSize tx_size,
                mode,                                                   //PredictionMode mode,
                0,                                                      //candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
                0,                                                      //int32_t use_palette,
                NULL,  //Inter-Intra
                FILTER_INTRA_MODES,                                     //CHKN FilterIntraMode filter_intra_mode,
                top_neigh_array + 1,
                left_neigh_array + 1,
                prediction_ptr,                                         //uint8_t *dst,
                (md_context_ptr->blk_geom->tx_org_x[is_inter][0][0] - md_context_ptr->blk_geom->origin_x) >> 2,
                (md_context_ptr->blk_geom->tx_org_y[is_inter][0][0] - md_context_ptr->blk_geom->origin_y) >> 2,
                PLANE_TYPE_Y,                                           //int32_t plane,
                md_context_ptr->blk_geom->bsize,                        //uint32_t puSize,
                md_context_ptr->blk_origin_x,
                md_context_ptr->blk_origin_y,
                md_context_ptr->blk_origin_x,                            //uint32_t cuOrgX,
                md_context_ptr->blk_origin_y,                            //uint32_t cuOrgY
                0,                                                      //cuOrgX used only for prediction Ptr
                0,                                                      //cuOrgY used only for prediction Ptr
                pcs_ptr->mi_grid_base,
                &((SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr)->seq_header
        );
    }

    return return_error;
}


#define USE_PADDING_FIX 1
EbErrorType update_neighbor_samples_array_open_loop_mb(
        uint8_t                            *above_ref,
        uint8_t                            *left_ref,
        EbPictureBufferDesc                *input_ptr,
        uint32_t                            stride,
        uint32_t                            src_origin_x,
        uint32_t                            src_origin_y,
        uint8_t                             bwidth,
        uint8_t                             bheight)
{
    EbErrorType    return_error = EB_ErrorNone;

    uint32_t idx;
    uint8_t  *src_ptr;
    uint8_t  *read_ptr;
    uint32_t count;

    uint32_t width = input_ptr->width;
    uint32_t height = input_ptr->height;
    uint32_t block_size_half = bwidth << 1;

    // Adjust the Source ptr to start at the origin of the block being updated
    src_ptr = input_ptr->buffer_y + (((src_origin_y + input_ptr->origin_y) * stride) + (src_origin_x + input_ptr->origin_x));

    //Initialise the Luma Intra Reference Array to the mid range value 128 (for CUs at the picture boundaries)
    EB_MEMSET(above_ref, 127, (bwidth << 1) + 1);
    EB_MEMSET(left_ref, 129, (bheight << 1) + 1);

    // Get the upper left sample
    if (src_origin_x != 0 && src_origin_y != 0) {
        read_ptr = src_ptr - stride - 1;
        *above_ref = *read_ptr;
        *left_ref = *read_ptr;
        left_ref++;
        above_ref++;
    }else {
        *above_ref = *left_ref = 128;
        left_ref++;
        above_ref++;
    }
    // Get the left-column
    count = block_size_half;
    if (src_origin_x != 0) {
        read_ptr = src_ptr - 1;
#if USE_PADDING_FIX
        if(src_origin_y == 0)
            *(left_ref - 1) = *read_ptr;
#endif
        count = ((src_origin_y + count) > height) ? count - ((src_origin_y + count) - height) : count;
        for (idx = 0; idx < count; ++idx) {
            *left_ref = *read_ptr;
            read_ptr += stride;
            left_ref++;
        }
        left_ref += (block_size_half - count);
#if USE_PADDING_FIX
        // pading unknown left bottom pixels with value at(-1, -15)
        for(idx = 0; idx < bheight; idx++)
            *(left_ref - bheight + idx) = *(left_ref - bheight - 1);
#endif
#if USE_PADDING_FIX
    } else if (src_origin_y != 0 ) {
        count = ((src_origin_y + count) > height) ? count - ((src_origin_y + count) - height) : count;
        EB_MEMSET(left_ref - 1, *(src_ptr - stride), count + 1);
        *(above_ref - 1) = *(src_ptr - stride);
#endif
    }else
        left_ref += count;

    // Get the top-row
    count = block_size_half;
    if (src_origin_y != 0) {
        read_ptr = src_ptr - stride;
        count = ((src_origin_x + count) > width) ? count - ((src_origin_x + count) - width) : count;
        EB_MEMCPY(above_ref, read_ptr, count);
#if USE_PADDING_FIX
        // pading unknown top right pixels with value at(15, -1)
        if(src_origin_x != 0)
            for(idx = 0; idx < bwidth; idx++)
                *(above_ref + bwidth + idx) = *(above_ref + bwidth - 1);
#else
        above_ref += (block_size_half - count);
#endif
#if USE_PADDING_FIX
    } else if (src_origin_x != 0 ) {
        count = ((src_origin_x + count) > width) ? count - ((src_origin_x + count) - width) : count;
        EB_MEMSET(above_ref - 1, *(left_ref - count), count + 1);
#endif
    }

    return return_error;
}

EbErrorType update_neighbor_samples_array_open_loop_mb_recon(
    uint8_t *above_ref, uint8_t *left_ref, uint8_t *recon_ptr, uint32_t stride,
    uint32_t src_origin_x, uint32_t src_origin_y, uint8_t bwidth, uint8_t bheight, uint32_t width,
    uint32_t height)
{
    EbErrorType    return_error = EB_ErrorNone;

    uint32_t idx;
    uint8_t  *src_ptr;
    uint8_t  *read_ptr;
    uint32_t count;

    uint32_t block_size_half = bwidth << 1;

    // Adjust the Source ptr to start at the origin of the block being updated
    src_ptr = recon_ptr + (src_origin_y  * stride + src_origin_x );

    //Initialise the Luma Intra Reference Array to the mid range value 128 (for CUs at the picture boundaries)
    EB_MEMSET(above_ref, 127, (bwidth << 1) + 1);
    EB_MEMSET(left_ref, 129, (bheight << 1) + 1);

    // Get the upper left sample
    if (src_origin_x != 0 && src_origin_y != 0) {
        read_ptr = src_ptr - stride - 1;
        *above_ref = *read_ptr;
        *left_ref = *read_ptr;
        left_ref++;
        above_ref++;
    }
    else {
        *above_ref = *left_ref = 128;
        left_ref++;
        above_ref++;
    }
    // Get the left-column
    count = block_size_half;
    if (src_origin_x != 0) {
        read_ptr = src_ptr - 1;
#if USE_PADDING_FIX
        if (src_origin_y == 0)
            *(left_ref - 1) = *read_ptr;
#endif
        count = ((src_origin_y + count) > height) ? count - ((src_origin_y + count) - height) : count;
        for (idx = 0; idx < count; ++idx) {
            *left_ref = *read_ptr;
            read_ptr += stride;
            left_ref++;
        }
        left_ref += (block_size_half - count);
#if USE_PADDING_FIX
        // pading unknown left bottom pixels with value at(-1, -15)
        for (idx = 0; idx < bheight; idx++)
            *(left_ref - bheight + idx) = *(left_ref - bheight - 1);
#endif
#if USE_PADDING_FIX
    }
    else
        if (src_origin_y != 0) {
            count = ((src_origin_y + count) > height) ? count - ((src_origin_y + count) - height) : count;
            EB_MEMSET(left_ref - 1, *(src_ptr - stride), count + 1);
            *(above_ref - 1) = *(src_ptr - stride);
#endif
        }
        else
            left_ref += count;

    // Get the top-row
    count = block_size_half;
    if (src_origin_y != 0) {
        read_ptr = src_ptr - stride;
        count = ((src_origin_x + count) > width) ? count - ((src_origin_x + count) - width) : count;
        EB_MEMCPY(above_ref, read_ptr, count);
#if USE_PADDING_FIX
        // pading unknown top right pixels with value at(15, -1)
        if (src_origin_x != 0)
            for (idx = 0; idx < bwidth; idx++)
                *(above_ref + bwidth + idx) = *(above_ref + bwidth - 1);
#else
        above_ref += (block_size_half - count);
#endif
#if USE_PADDING_FIX
    }
    else
        if (src_origin_x != 0) {
            count = ((src_origin_x + count) > width) ? count - ((src_origin_x + count) - width) : count;
            EB_MEMSET(above_ref - 1, *(left_ref - count), count + 1);
#endif
        }

    return return_error;
}

