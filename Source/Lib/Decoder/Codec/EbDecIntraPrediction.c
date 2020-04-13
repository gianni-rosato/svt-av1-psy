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

#include <string.h>

#include "common_dsp_rtcd.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbDecParseHelper.h"
#include "EbDecProcessFrame.h"
#include "EbCommonUtils.h"
#include "EbIntraPrediction.h"

/* Avoid 12-bit output mismatch by intra pred intrinsic kernel */
void dec_init_intra_predictors_12b_internal(void) {
    eb_av1_highbd_dr_prediction_z2 = eb_av1_highbd_dr_prediction_z2_c;
}

/*TODO: Remove replication and harmonize with encoder after data str. harmonization */
int32_t dec_get_filt_type(const PartitionInfo *part_info, int32_t plane) {
    int ab_sm, le_sm;

    if (plane == 0) {
        const BlockModeInfo *ab = part_info->above_mbmi;
        const BlockModeInfo *le = part_info->left_mbmi;
        ab_sm                   = ab ? is_smooth(ab, plane) : 0;
        le_sm                   = le ? is_smooth(le, plane) : 0;
    } else {
        const BlockModeInfo *ab = part_info->chroma_above_mbmi;
        const BlockModeInfo *le = part_info->chroma_left_mbmi;
        ab_sm                   = ab ? is_smooth(ab, plane) : 0;
        le_sm                   = le ? is_smooth(le, plane) : 0;
    }

    return (ab_sm || le_sm) ? 1 : 0;
}

void cfl_init(CflCtx *cfl, EbColorConfig *cc) {
    assert(block_size_wide[CFL_MAX_BlockSize] == CFL_BUF_LINE);
    assert(block_size_high[CFL_MAX_BlockSize] == CFL_BUF_LINE);

    memset(&cfl->recon_buf_q3, 0, sizeof(cfl->recon_buf_q3));
    cfl->subsampling_x           = cc->subsampling_x;
    cfl->subsampling_y           = cc->subsampling_y;
    cfl->are_parameters_computed = 0;
}

// Due to frame boundary issues, it is possible that the total area covered by
// chroma exceeds that of luma. When this happens, we fill the missing pixels by
// repeating the last columns and/or rows.

static INLINE void cfl_pad(CflCtx *cfl, int32_t width, int32_t height) {
    const int32_t diff_width  = width - cfl->buf_width;
    const int32_t diff_height = height - cfl->buf_height;

    if (diff_width > 0) {
        const int min_height   = height - diff_height;
        int16_t * recon_buf_q3 = cfl->recon_buf_q3 + (width - diff_width);
        for (int j = 0; j < min_height; j++) {
            const int16_t last_pixel = recon_buf_q3[-1];
            assert(recon_buf_q3 + diff_width <= cfl->recon_buf_q3 + CFL_BUF_SQUARE);
            for (int i = 0; i < diff_width; i++) recon_buf_q3[i] = last_pixel;
            recon_buf_q3 += CFL_BUF_LINE;
        }
        cfl->buf_width = width;
    }
    if (diff_height > 0) {
        int16_t *recon_buf_q3 = cfl->recon_buf_q3 + ((height - diff_height) * CFL_BUF_LINE);
        for (int j = 0; j < diff_height; j++) {
            const int16_t *last_row_q3 = recon_buf_q3 - CFL_BUF_LINE;
            assert(recon_buf_q3 + width <= cfl->recon_buf_q3 + CFL_BUF_SQUARE);
            for (int i = 0; i < width; i++) recon_buf_q3[i] = last_row_q3[i];
            recon_buf_q3 += CFL_BUF_LINE;
        }
        cfl->buf_height = height;
    }
}

void cfl_luma_subsampling_422_lbd_c(const uint8_t *input, int32_t input_stride, int16_t *output_q3,
                                    int32_t width, int32_t height) {
    assert((height - 1) * CFL_BUF_LINE + width <= CFL_BUF_SQUARE);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i += 2) output_q3[i >> 1] = (input[i] + input[i + 1]) << 2;
        input += input_stride;
        output_q3 += CFL_BUF_LINE;
    }
}

void cfl_luma_subsampling_444_lbd_c(const uint8_t *input, int32_t input_stride, int16_t *output_q3,
                                    int32_t width, int32_t height) {
    assert((height - 1) * CFL_BUF_LINE + width <= CFL_BUF_SQUARE);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) output_q3[i] = input[i] << 3;
        input += input_stride;
        output_q3 += CFL_BUF_LINE;
    }
}

void cfl_luma_subsampling_422_hbd_c(const uint16_t *input, int32_t input_stride, int16_t *output_q3,
                                    int32_t width, int32_t height) {
    assert((height - 1) * CFL_BUF_LINE + width <= CFL_BUF_SQUARE);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i += 2) output_q3[i >> 1] = (input[i] + input[i + 1]) << 2;
        input += input_stride;
        output_q3 += CFL_BUF_LINE;
    }
}

void cfl_luma_subsampling_444_hbd_c(const uint16_t *input, int32_t input_stride, int16_t *output_q3,
                                    int32_t width, int32_t height) {
    assert((height - 1) * CFL_BUF_LINE + width <= CFL_BUF_SQUARE);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) output_q3[i] = input[i] << 3;
        input += input_stride;
        output_q3 += CFL_BUF_LINE;
    }
}

static void cfl_subsampling_highbd(TxSize tx_size, int32_t sub_x, int32_t sub_y,
                                  const uint16_t *input, int32_t input_stride,
                                  int16_t *recon_buf_q3)
{
    int32_t width  = tx_size_wide[tx_size];
    int32_t height = tx_size_high[tx_size];
    assert(width != 64 || height != 64);
    if (sub_x == 1 && sub_y == 1) {
        cfl_luma_subsampling_420_hbd(input, input_stride, recon_buf_q3,
            width, height);
    }
    else if (sub_x == 1 && sub_y == 0) {
        cfl_luma_subsampling_422_hbd_c(input, input_stride, recon_buf_q3,
            width, height);
}
    else {
        cfl_luma_subsampling_444_hbd_c(input, input_stride, recon_buf_q3,
            width, height);
    }
}

static void cfl_subsampling_lowbd(TxSize tx_size, int32_t sub_x, int32_t sub_y,
                                  const uint8_t *input, int32_t input_stride,
                                  int16_t *recon_buf_q3)
{
    int32_t width  = tx_size_wide[tx_size];
    int32_t height = tx_size_high[tx_size];
    assert(width != 64 || height != 64);
    if (sub_x == 1 && sub_y == 1) {
        cfl_luma_subsampling_420_lbd(input, input_stride, recon_buf_q3,
            width, height);
    }
    else if (sub_x == 1 && sub_y == 0) {
        cfl_luma_subsampling_422_lbd_c(input, input_stride, recon_buf_q3,
            width, height);
    }
    else {
        cfl_luma_subsampling_444_lbd_c(input, input_stride, recon_buf_q3,
            width, height);
    }
}
//######...........Ending for CFL.................#####//

//####...Wrapper funtion calling CFL leaf level functions...####//
static INLINE CflAllowedType is_cfl_allowed_with_frame_header(const PartitionInfo *xd,
                                                              EbColorConfig *cc, FrameHeader *fh)

{
    const BlockModeInfo *mbmi  = xd->mi;
    const BlockSize      bsize = mbmi->sb_type;
    assert(bsize < BlockSizeS_ALL);
    if (fh->lossless_array[mbmi->segment_id]) {
        // In lossless, CfL is available when the partition size is equal to the
        // transform size.
        const int ssx         = cc->subsampling_x;
        const int ssy         = cc->subsampling_y;
        const int plane_bsize = get_plane_block_size(bsize, ssx, ssy);
        return (CflAllowedType)(plane_bsize == BLOCK_4X4);
    }
    // Spec: CfL is available to luma partitions lesser than or equal to 32x32
    return (CflAllowedType)(block_size_wide[bsize] <= 32 && block_size_high[bsize] <= 32);
}

void cfl_compute_parameters(CflCtx *cfl_ctx, TxSize tx_size) {
    //CFL_CTX *const cfl = &xd->cfl;
    // Do not call cfl_compute_parameters multiple time on the same values.
    assert(cfl_ctx->are_parameters_computed == 0);
    cfl_pad(cfl_ctx, tx_size_wide[tx_size], tx_size_high[tx_size]);
    get_subtract_average_fn(tx_size)(cfl_ctx->recon_buf_q3);
    cfl_ctx->are_parameters_computed = 1;
}

void cfl_predict_block(PartitionInfo *xd, CflCtx *cfl_ctx, uint8_t *dst, int32_t dst_stride,
                       TxSize tx_size, int32_t plane, EbColorConfig *cc, FrameHeader *fh, EbBool is_16bit) {
    BlockModeInfo *mbmi                = xd->mi;
    CflAllowedType is_cfl_allowed_flag = is_cfl_allowed_with_frame_header(xd, cc, fh);
    assert(is_cfl_allowed_flag == CFL_ALLOWED);
    (void)is_cfl_allowed_flag;

    if (!cfl_ctx->are_parameters_computed) cfl_compute_parameters(cfl_ctx, tx_size);

    const int32_t alpha_q3 =
        cfl_idx_to_alpha(mbmi->cfl_alpha_idx, mbmi->cfl_alpha_signs, plane - 1);
    assert((tx_size_high[tx_size] - 1) * CFL_BUF_LINE + tx_size_wide[tx_size] <= CFL_BUF_SQUARE);

    if ((cc->bit_depth != EB_8BIT) || is_16bit) {
        eb_cfl_predict_hbd(cfl_ctx->recon_buf_q3,
                           (uint16_t *)dst,
                           dst_stride,
                           (uint16_t *)dst,
                           dst_stride,
                           alpha_q3,
                           cc->bit_depth,
                           tx_size_wide[tx_size],
                           tx_size_high[tx_size]);
        return;
    }

    eb_cfl_predict_lbd(cfl_ctx->recon_buf_q3,
                       dst,
                       dst_stride,
                       dst,
                       dst_stride,
                       alpha_q3,
                       cc->bit_depth,
                       tx_size_wide[tx_size],
                       tx_size_high[tx_size]);
}

static void cfl_store(CflCtx *cfl_ctx, const uint8_t *input, int input_stride, int row, int col,
                      TxSize tx_size, uint8_t use_hbd) {
    const int width        = tx_size_wide[tx_size];
    const int height       = tx_size_high[tx_size];
    const int tx_off_log2  = tx_size_wide_log2[0];
    const int sub_x        = cfl_ctx->subsampling_x;
    const int sub_y        = cfl_ctx->subsampling_y;
    const int store_row    = row << (tx_off_log2 - sub_y);
    const int store_col    = col << (tx_off_log2 - sub_x);
    const int store_height = height >> sub_y;
    const int store_width  = width >> sub_x;

    // Invalidate current parameters
    cfl_ctx->are_parameters_computed = 0;

    // Store the surface of the pixel buffer that was written to, this way we
    // can manage chroma overrun (e.g. when the chroma surfaces goes beyond the
    // frame boundary)
    if (col == 0 && row == 0) {
        cfl_ctx->buf_width  = store_width;
        cfl_ctx->buf_height = store_height;
    } else {
        cfl_ctx->buf_width  = AOMMAX(store_col + store_width, cfl_ctx->buf_width);
        cfl_ctx->buf_height = AOMMAX(store_row + store_height, cfl_ctx->buf_height);
    }

    // Check that we will remain inside the pixel buffer.
    assert(store_row + store_height <= CFL_BUF_LINE);
    assert(store_col + store_width <= CFL_BUF_LINE);

    // Store the input into the CfL pixel buffer
    int16_t *recon_buf_q3 = cfl_ctx->recon_buf_q3 + (store_row * CFL_BUF_LINE + store_col);

    if (use_hbd) {
        cfl_subsampling_highbd(tx_size, sub_x, sub_y, (uint16_t *)input,
            input_stride, recon_buf_q3);
    } else {
        cfl_subsampling_lowbd(tx_size, sub_x, sub_y ,input,
            input_stride, recon_buf_q3);
    }
}

// Adjust the row and column of blocks smaller than 8X8, as chroma-referenced
// and non-chroma-referenced blocks are stored together in the CfL buffer.
static INLINE void sub8x8_adjust_offset(PartitionInfo *xd, const CflCtx *cfl_ctx, int *row_out,
                                        int *col_out) {
    // Increment row index for bottom: 8x4, 16x4 or both bottom 4x4s.
    if ((xd->mi_row & 0x01) && cfl_ctx->subsampling_y) {
        assert(*row_out == 0);
        (*row_out)++;
    }

    // Increment col index for right: 4x8, 4x16 or both right 4x4s.
    if ((xd->mi_col & 0x01) && cfl_ctx->subsampling_x) {
        assert(*col_out == 0);
        (*col_out)++;
    }
}

void cfl_store_tx(PartitionInfo *xd, CflCtx *cfl_ctx, int row, int col, TxSize tx_size,
                  BlockSize bsize, EbColorConfig *cc, uint8_t *dst_buff, uint32_t dst_stride, EbBool is_16bit) {
    if (block_size_high[bsize] == 4 || block_size_wide[bsize] == 4) {
        // Only dimensions of size 4 can have an odd offset.
        assert(!((col & 1) && tx_size_wide[tx_size] != 4));
        assert(!((row & 1) && tx_size_high[tx_size] != 4));
        sub8x8_adjust_offset(xd, cfl_ctx, &row, &col);
    }

    cfl_store(cfl_ctx, dst_buff, dst_stride, row, col, tx_size, ((cc->bit_depth != EB_8BIT) || is_16bit));
}
//#####.....................Ending for wrapper of CFL...............................####//

/* TODO : Harmonize with Encoder! */
static void decode_build_intra_predictors(PartitionInfo *part_info, uint8_t *top_neigh_array,
                                          uint8_t *left_neigh_array, int32_t ref_stride,
                                          uint8_t *dst, int32_t dst_stride, PredictionMode mode,
                                          int32_t angle_delta, FilterIntraMode filter_intra_mode,
                                          TxSize tx_size, int32_t disable_edge_filter,
                                          int32_t n_top_px, int32_t n_topright_px,
                                          int32_t n_left_px, int32_t n_bottomleft_px,
                                          int32_t plane) {
    int32_t i;

    const uint8_t *above_ref = top_neigh_array; //CHKN ref - ref_stride;
    const uint8_t *left_ref  = left_neigh_array; //CHKN ref - 1;

    DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);
    uint8_t *const above_row = above_data + 16;
    uint8_t *const left_col  = left_data + 16;

    const int32_t txwpx            = tx_size_wide[tx_size];
    const int32_t txhpx            = tx_size_high[tx_size];
    int32_t       need_left        = extend_modes[mode] & NEED_LEFT;
    int32_t       need_above       = extend_modes[mode] & NEED_ABOVE;
    int32_t       need_above_left  = extend_modes[mode] & NEED_ABOVELEFT;
    int32_t       p_angle          = 0;
    const int32_t is_dr_mode       = av1_is_directional_mode(mode);
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
        i                                    = 0;
        if (n_left_px > 0) {
            for (; i < n_left_px; i++) left_col[i] = left_ref[i * ref_stride];
            if (need_bottom && n_bottomleft_px > 0) {
                assert(i == txhpx);
                for (; i < txhpx + n_bottomleft_px; i++) left_col[i] = left_ref[i * ref_stride];
            }
            if (i < num_left_pixels_needed)
                memset(&left_col[i], left_col[i - 1], num_left_pixels_needed - i);
        } else {
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
            memcpy(above_row, above_ref, n_top_px);
            i = n_top_px;
            if (need_right && n_topright_px > 0) {
                assert(n_top_px == txwpx);
                memcpy(above_row + txwpx, above_ref + txwpx, n_topright_px);
                i += n_topright_px;
            }
            if (i < num_top_pixels_needed)
                memset(&above_row[i], above_row[i - 1], num_top_pixels_needed - i);
        } else {
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
        eb_av1_filter_intra_predictor(
            dst, dst_stride, tx_size, above_row, left_col, filter_intra_mode);
        return;
    }

    if (is_dr_mode) {
        int32_t upsample_above = 0;
        int32_t upsample_left  = 0;
        if (!disable_edge_filter) {
            const int32_t need_right  = p_angle < 90;
            const int32_t need_bottom = p_angle > 180;

            const int32_t filt_type = dec_get_filt_type(part_info, plane);

            if (p_angle != 90 && p_angle != 180) {
                const int32_t ab_le = need_above_left ? 1 : 0;
                if (need_above && need_left && (txwpx + txhpx >= 24))
                    filter_intra_edge_corner(above_row, left_col);
                if (need_above && n_top_px > 0) {
                    const int32_t strength =
                        intra_edge_filter_strength(txwpx, txhpx, p_angle - 90, filt_type);
                    const int32_t n_px = n_top_px + ab_le + (need_right ? txhpx : 0);
                    eb_av1_filter_intra_edge(above_row - ab_le, n_px, strength);
                }
                if (need_left && n_left_px > 0) {
                    const int32_t strength =
                        intra_edge_filter_strength(txhpx, txwpx, p_angle - 180, filt_type);
                    const int32_t n_px = n_left_px + ab_le + (need_bottom ? txwpx : 0);
                    eb_av1_filter_intra_edge(left_col - ab_le, n_px, strength);
                }
            }
            upsample_above = use_intra_edge_upsample(txwpx, txhpx, p_angle - 90, filt_type);
            if (need_above && upsample_above) {
                const int32_t n_px = txwpx + (need_right ? txhpx : 0);

                eb_av1_upsample_intra_edge(above_row, n_px);
            }
            upsample_left = use_intra_edge_upsample(txhpx, txwpx, p_angle - 180, filt_type);
            if (need_left && upsample_left) {
                const int32_t n_px = txhpx + (need_bottom ? txwpx : 0);

                eb_av1_upsample_intra_edge(left_col, n_px);
            }
        }
        dr_predictor(
            dst, dst_stride, tx_size, above_row, left_col, upsample_above, upsample_left, p_angle);
        return;
    }

    // predict
    if (mode == DC_PRED) {
        dc_pred[n_left_px > 0][n_top_px > 0][tx_size](dst, dst_stride, above_row, left_col);
    } else
        pred[mode][tx_size](dst, dst_stride, above_row, left_col);
}

/* TODO : Harmonize with Encoder! */
static void decode_build_intra_predictors_high(PartitionInfo *part_info, uint16_t *top_neigh_array,
                                               uint16_t *left_neigh_array, int32_t ref_stride,
                                               uint16_t *dst, int32_t dst_stride,
                                               PredictionMode mode, int32_t angle_delta,
                                               FilterIntraMode filter_intra_mode, TxSize tx_size,
                                               int32_t disable_edge_filter, int32_t n_top_px,
                                               int32_t n_topright_px, int32_t n_left_px,
                                               int32_t n_bottomleft_px, int32_t plane, int32_t bd) {
    int32_t i;

    DECLARE_ALIGNED(16, uint16_t, left_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(16, uint16_t, above_data[MAX_TX_SIZE * 2 + 32]);
    uint16_t *const above_row       = above_data + 16;
    uint16_t *const left_col        = left_data + 16;
    const int32_t   txwpx           = tx_size_wide[tx_size];
    const int32_t   txhpx           = tx_size_high[tx_size];
    int32_t         need_left       = extend_modes[mode] & NEED_LEFT;
    int32_t         need_above      = extend_modes[mode] & NEED_ABOVE;
    int32_t         need_above_left = extend_modes[mode] & NEED_ABOVELEFT;

    const uint16_t *above_ref = top_neigh_array; //CHKN ref - ref_stride;
    const uint16_t *left_ref  = left_neigh_array; //CHKN ref - 1;

    int32_t       p_angle          = 0;
    const int32_t is_dr_mode       = av1_is_directional_mode(mode);
    const int32_t use_filter_intra = filter_intra_mode != FILTER_INTRA_MODES;
    int32_t       base             = 128 << (bd - 8);

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
            eb_aom_memset16(dst, val, txwpx);
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
        i                                    = 0;
        if (n_left_px > 0) {
            for (; i < n_left_px; i++) left_col[i] = left_ref[i * ref_stride];
            if (need_bottom && n_bottomleft_px > 0) {
                assert(i == txhpx);
                for (; i < txhpx + n_bottomleft_px; i++) left_col[i] = left_ref[i * ref_stride];
            }
            if (i < num_left_pixels_needed)
                eb_aom_memset16(&left_col[i], left_col[i - 1], num_left_pixels_needed - i);
        } else {
            if (n_top_px > 0)
                eb_aom_memset16(left_col, above_ref[0], num_left_pixels_needed);
            else
                eb_aom_memset16(left_col, base + 1, num_left_pixels_needed);
        }
    }

    // NEED_ABOVE
    if (need_above) {
        int32_t need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
        if (use_filter_intra) need_right = 0;
        if (is_dr_mode) need_right = p_angle < 90;
        const int32_t num_top_pixels_needed = txwpx + (need_right ? txhpx : 0);
        if (n_top_px > 0) {
            memcpy(above_row, above_ref, n_top_px * sizeof(above_ref[0]));
            i = n_top_px;
            if (need_right && n_topright_px > 0) {
                assert(n_top_px == txwpx);
                memcpy(above_row + txwpx, above_ref + txwpx, n_topright_px * sizeof(above_ref[0]));
                i += n_topright_px;
            }
            if (i < num_top_pixels_needed)
                eb_aom_memset16(&above_row[i], above_row[i - 1], num_top_pixels_needed - i);
        } else {
            if (n_left_px > 0)
                eb_aom_memset16(above_row, left_ref[0], num_top_pixels_needed);
            else
                eb_aom_memset16(above_row, base - 1, num_top_pixels_needed);
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
        highbd_filter_intra_predictor(
            dst, dst_stride, tx_size, above_row, left_col, filter_intra_mode, bd);
        return;
    }

    if (is_dr_mode) {
        int32_t upsample_above = 0;
        int32_t upsample_left  = 0;
        if (!disable_edge_filter) {
            const int32_t need_right  = p_angle < 90;
            const int32_t need_bottom = p_angle > 180;
            const int32_t filt_type   = dec_get_filt_type(part_info, plane);
            if (p_angle != 90 && p_angle != 180) {
                const int32_t ab_le = need_above_left ? 1 : 0;
                if (need_above && need_left && (txwpx + txhpx >= 24))
                    filter_intra_edge_corner_high(above_row, left_col);
                if (need_above && n_top_px > 0) {
                    const int32_t strength =
                        intra_edge_filter_strength(txwpx, txhpx, p_angle - 90, filt_type);
                    const int32_t n_px = n_top_px + ab_le + (need_right ? txhpx : 0);
                    eb_av1_filter_intra_edge_high(above_row - ab_le, n_px, strength);
                }
                if (need_left && n_left_px > 0) {
                    const int32_t strength =
                        intra_edge_filter_strength(txhpx, txwpx, p_angle - 180, filt_type);
                    const int32_t n_px = n_left_px + ab_le + (need_bottom ? txwpx : 0);

                    eb_av1_filter_intra_edge_high(left_col - ab_le, n_px, strength);
                }
            }
            upsample_above = use_intra_edge_upsample(txwpx, txhpx, p_angle - 90, filt_type);
            if (need_above && upsample_above) {
                const int32_t n_px = txwpx + (need_right ? txhpx : 0);
                //av1_upsample_intra_edge_high(above_row, n_px, bd);// AMIR : to be replaced by optimized code
                eb_av1_upsample_intra_edge_high_c(above_row, n_px, bd);
            }
            upsample_left = use_intra_edge_upsample(txhpx, txwpx, p_angle - 180, filt_type);
            if (need_left && upsample_left) {
                const int32_t n_px = txhpx + (need_bottom ? txwpx : 0);
                //av1_upsample_intra_edge_high(left_col, n_px, bd);// AMIR: to be replaced by optimized code
                eb_av1_upsample_intra_edge_high_c(left_col, n_px, bd);
            }
        }
        highbd_dr_predictor(dst,
                            dst_stride,
                            tx_size,
                            above_row,
                            left_col,
                            upsample_above,
                            upsample_left,
                            p_angle,
                            bd);
        return;
    }

    // predict
    if (mode == DC_PRED) {
        dc_pred_high[n_left_px > 0][n_top_px > 0][tx_size](
            dst, dst_stride, above_row, left_col, bd);
    } else
        pred_high[mode][tx_size](dst, dst_stride, above_row, left_col, bd);
}

void svtav1_predict_intra_block(PartitionInfo *xd, int32_t plane, TxSize tx_size, TileInfo *td,
                                void *pv_pred_buf, int32_t pred_stride, void *top_neigh_array,
                                void *left_neigh_array, int32_t ref_stride, SeqHeader *seq_header,
                                const PredictionMode mode, int32_t blk_mi_col_off,
                                int32_t blk_mi_row_off, EbBitDepthEnum bit_depth, EbBool is_16bit) {
    //ToDo:are_parameters_computed variable for CFL so that cal part for V plane we can skip,
    //once we compute for U plane, this parameter is block level parameter.
    const EbColorConfig *cc    = &seq_header->color_config;
    int32_t              sub_x = plane ? cc->subsampling_x : 0;
    int32_t              sub_y = plane ? cc->subsampling_y : 0;

    const BlockModeInfo *const mbmi = xd->mi;

    const int txwpx = tx_size_wide[tx_size];
    const int txhpx = tx_size_high[tx_size];

    const int use_palette = mbmi->palette_size[plane != 0] > 0;

    if (use_palette) return;
    const FilterIntraMode filter_intra_mode =
        (plane == AOM_PLANE_Y && mbmi->filter_intra_mode_info.use_filter_intra)
            ? mbmi->filter_intra_mode_info.filter_intra_mode
            : FILTER_INTRA_MODES;

    const int angle_delta = mbmi->angle_delta[plane != AOM_PLANE_Y];

    BlockSize bsize = mbmi->sb_type;
    //const struct macroblockd_plane *const pd = &xd->plane[plane];
    const int txw      = tx_size_wide_unit[tx_size];
    const int txh      = tx_size_high_unit[tx_size];
    const int have_top = blk_mi_row_off || (sub_y ? xd->chroma_up_available : xd->up_available);
    const int have_left =
        blk_mi_col_off || (sub_x ? xd->chroma_left_available : xd->left_available);

    const int mi_row        = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
    const int mi_col        = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);
    const int xr_chr_offset = 0;
    const int yd_chr_offset = 0;

    // Distance between right edge of this pred block to frame right edge
    const int xr = (xd->mb_to_right_edge >> (3 + sub_x)) +
                   (xd->wpx[plane] - (blk_mi_col_off << MI_SIZE_LOG2) - txwpx) - xr_chr_offset;
    // Distance between bottom edge of this pred block to frame bottom edge
    const int yd = (xd->mb_to_bottom_edge >> (3 + sub_y)) +
                   (xd->hpx[plane] - (blk_mi_row_off << MI_SIZE_LOG2) - txhpx) - yd_chr_offset;
    const int right_available = mi_col + ((blk_mi_col_off + txw) << sub_x) < td->mi_col_end;
    const int bottom_available =
        (yd > 0) && (mi_row + ((blk_mi_row_off + txh) << sub_y) < td->mi_row_end);

    const PartitionType partition = mbmi->partition;

    // force 4x4 chroma component block size.
    bsize = scale_chroma_bsize(bsize, sub_x, sub_y);

    const int have_top_right   = intra_has_top_right(seq_header->sb_size,
                                                   bsize,
                                                   mi_row,
                                                   mi_col,
                                                   have_top,
                                                   right_available,
                                                   partition,
                                                   tx_size,
                                                   blk_mi_row_off,
                                                   blk_mi_col_off,
                                                   sub_x,
                                                   sub_y);
    const int have_bottom_left = intra_has_bottom_left(seq_header->sb_size,
                                                       bsize,
                                                       mi_row,
                                                       mi_col,
                                                       bottom_available,
                                                       have_left,
                                                       partition,
                                                       tx_size,
                                                       blk_mi_row_off,
                                                       blk_mi_col_off,
                                                       sub_x,
                                                       sub_y);

    const int32_t disable_edge_filter = !seq_header->enable_intra_edge_filter;

    //###..Calling all other intra predictors except CFL & pallate...//
    if (bit_depth == EB_8BIT && !is_16bit) {
        decode_build_intra_predictors(xd,
                                      (uint8_t *)top_neigh_array, /*As per SVT Enc*/
                                      (uint8_t *)left_neigh_array,
                                      /*As per SVT Enc*/ ref_stride,
                                      (uint8_t *)pv_pred_buf,
                                      pred_stride,
                                      mode,
                                      angle_delta,
                                      filter_intra_mode,
                                      tx_size,
                                      disable_edge_filter,
                                      have_top ? AOMMIN(txwpx, xr + txwpx) : 0,
                                      have_top_right ? AOMMIN(txwpx, xr) : 0,
                                      have_left ? AOMMIN(txhpx, yd + txhpx) : 0,
                                      have_bottom_left ? AOMMIN(txhpx, yd) : 0,
                                      plane);
    } else { //16bit
        decode_build_intra_predictors_high(xd,
                                           (uint16_t *)top_neigh_array, /*As per SVT Enc*/
                                           (uint16_t *)left_neigh_array,
                                           /*As per SVT Enc*/ ref_stride,
                                           (uint16_t *)pv_pred_buf,
                                           pred_stride,
                                           mode,
                                           angle_delta,
                                           filter_intra_mode,
                                           tx_size,
                                           disable_edge_filter,
                                           have_top ? AOMMIN(txwpx, xr + txwpx) : 0,
                                           have_top_right ? AOMMIN(txwpx, xr) : 0,
                                           have_left ? AOMMIN(txhpx, yd + txhpx) : 0,
                                           have_bottom_left ? AOMMIN(txhpx, yd) : 0,
                                           plane,
                                           bit_depth);
    }
}

void svt_av1_predict_intra(DecModCtxt *dec_mod_ctxt, PartitionInfo *part_info, int32_t plane,
                           TxSize tx_size, TileInfo *td, void *pv_blk_recon_buf,
                           int32_t recon_stride, EbBitDepthEnum bit_depth, int32_t blk_mi_col_off,
                           int32_t blk_mi_row_off) {
    void *pv_top_neighbor_array, *pv_left_neighbor_array;

    EbDecHandle *dec_handle = (EbDecHandle *)dec_mod_ctxt->dec_handle_ptr;
    EbBool is16b = dec_handle->is_16bit_pipeline;
    const PredictionMode mode =
        (plane == AOM_PLANE_Y) ? part_info->mi->mode : get_uv_mode(part_info->mi->uv_mode);

    if (bit_depth == EB_8BIT && !is16b) {
        EbByte buf             = (EbByte)pv_blk_recon_buf;
        pv_top_neighbor_array  = (void *)(buf - recon_stride);
        pv_left_neighbor_array = (void *)(buf - 1);
    } else { //16bit
        uint16_t *buf          = (uint16_t *)pv_blk_recon_buf;
        pv_top_neighbor_array  = (void *)(buf - recon_stride);
        pv_left_neighbor_array = (void *)(buf - 1);
    }

    if (plane != AOM_PLANE_Y && part_info->mi->uv_mode == UV_CFL_PRED) {
        svtav1_predict_intra_block(part_info,
                                   plane,
                                   tx_size,
                                   td,
                                   pv_blk_recon_buf,
                                   recon_stride,
                                   pv_top_neighbor_array,
                                   pv_left_neighbor_array,
                                   recon_stride,
                                   dec_mod_ctxt->seq_header,
                                   mode,
                                   blk_mi_col_off,
                                   blk_mi_row_off,
                                   bit_depth,
                                   is16b);

        cfl_predict_block(part_info,
                          part_info->pv_cfl_ctxt,
                          pv_blk_recon_buf,
                          recon_stride,
                          tx_size,
                          plane,
                          &dec_mod_ctxt->seq_header->color_config,
                          dec_mod_ctxt->frame_header,
                          is16b);

        return;
    }

    svtav1_predict_intra_block(part_info,
                               plane,
                               tx_size,
                               td,
                               pv_blk_recon_buf,
                               recon_stride,
                               pv_top_neighbor_array,
                               pv_left_neighbor_array,
                               recon_stride,
                               dec_mod_ctxt->seq_header,
                               mode,
                               blk_mi_col_off,
                               blk_mi_row_off,
                               bit_depth,
                               is16b);
}
