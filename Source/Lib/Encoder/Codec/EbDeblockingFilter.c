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

#include "EbDeblockingFilter.h"
#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbCodingUnit.h"
#include "EbSequenceControlSet.h"
#include "EbReferenceObject.h"
#include "EbCommonUtils.h"
//#include "EbLog.h"

void eb_av1_loop_filter_init(PictureControlSet *pcs_ptr) {
    //assert(MB_MODE_COUNT == n_elements(mode_lf_lut));
    LoopFilterInfoN *  lfi = &pcs_ptr->parent_pcs_ptr->lf_info;
    struct LoopFilter *lf  = &pcs_ptr->parent_pcs_ptr->frm_hdr.loop_filter_params;
    int32_t            lvl;

    lf->combine_vert_horz_lf = 1;

    // init limits for given sharpness
    update_sharpness(lfi, lf->sharpness_level);

    // init hev threshold const vectors
    for (lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++)
        memset(lfi->lfthr[lvl].hev_thr, (lvl >> 4), SIMD_WIDTH);
}


//***************************************************************************************************//

static INLINE int32_t scaled_buffer_offset(int32_t x_offset, int32_t y_offset, int32_t stride/*,
    const struct scale_factors *sf*/) {
    const int32_t x =
        /*sf ? sf->scale_value_x(x_offset, sf) >> SCALE_EXTRA_BITS :*/ x_offset;
    const int32_t y =
        /*sf ? sf->scale_value_y(y_offset, sf) >> SCALE_EXTRA_BITS :*/ y_offset;
    return y * stride + x;
}
static INLINE void setup_pred_plane(struct Buf2D *dst, BlockSize bsize, uint8_t *src, int32_t width,
                                    int32_t height, int32_t stride, int32_t mi_row, int32_t mi_col,
                                    /*const struct scale_factors *scale,*/
                                    int32_t subsampling_x, int32_t subsampling_y,
                                    int32_t is_16bit) {
    // Offset the buffer pointer
    if (subsampling_y && (mi_row & 0x01) && (mi_size_high[bsize] == 1)) mi_row -= 1;
    if (subsampling_x && (mi_col & 0x01) && (mi_size_wide[bsize] == 1)) mi_col -= 1;

    const int32_t x = (MI_SIZE * mi_col) >> subsampling_x;
    const int32_t y = (MI_SIZE * mi_row) >> subsampling_y;
    dst->buf        = src + (scaled_buffer_offset(x, y, stride /*, scale*/) << is_16bit);
    dst->buf0       = src;
    dst->width      = width;
    dst->height     = height;
    dst->stride     = stride;
}
void eb_av1_setup_dst_planes(struct MacroblockdPlane *planes, BlockSize bsize,
                             //const Yv12BufferConfig *src,
                             const EbPictureBufferDesc *src, int32_t mi_row, int32_t mi_col,
                             const int32_t plane_start, const int32_t plane_end) {
    // We use AOMMIN(num_planes, MAX_MB_PLANE) instead of num_planes to quiet
    // the static analysis warnings.
    //for (int32_t i = plane_start; i < AOMMIN(plane_end, MAX_MB_PLANE); ++i) {
    //    struct MacroblockdPlane *const pd = &planes[i];
    //    const int32_t is_uv = i > 0;
    //    setup_pred_plane(&pd->dst, bsize, src->buffers[i], src->crop_widths[is_uv],
    //        src->crop_heights[is_uv], src->strides[is_uv], mi_row,
    //        mi_col, NULL, pd->subsampling_x, pd->subsampling_y);
    //}
    for (int32_t i = plane_start; i < AOMMIN(plane_end, 3); ++i) {
        if (i == 0) {
            struct MacroblockdPlane *const pd = &planes[0];
            setup_pred_plane(
                &pd->dst,
                bsize,
                &src->buffer_y[(src->origin_x + src->origin_y * src->stride_y) << pd->is_16bit],
                src->width,
                src->height,
                src->stride_y,
                mi_row,
                mi_col,
                /*NULL,*/ pd->subsampling_x,
                pd->subsampling_y,
                pd->is_16bit); //AMIR: Updated to point to the right location
        } else if (i == 1) {
            struct MacroblockdPlane *const pd = &planes[1];
            setup_pred_plane(
                &pd->dst,
                bsize,
                &src->buffer_cb[((src->origin_x + src->origin_y * src->stride_cb) << pd->is_16bit) /
                                2],
                src->width / 2,
                src->height / 2,
                src->stride_cb,
                mi_row,
                mi_col,
                /*NULL,*/ pd->subsampling_x,
                pd->subsampling_y,
                pd->is_16bit);
        } else if (i == 2) {
            struct MacroblockdPlane *const pd = &planes[2];
            setup_pred_plane(
                &pd->dst,
                bsize,
                &src->buffer_cr[((src->origin_x + src->origin_y * src->stride_cr) << pd->is_16bit) /
                                2],
                src->width / 2,
                src->height / 2,
                src->stride_cr,
                mi_row,
                mi_col,
                /* NULL,*/ pd->subsampling_x,
                pd->subsampling_y,
                pd->is_16bit);
        }
    }
}


//***************************************************************************************************//

static TxSize get_transform_size(const MacroBlockD *const xd, const MbModeInfo *const mbmi,
                                 const EdgeDir edge_dir, const int32_t mi_row, const int32_t mi_col,
                                 const int32_t plane, const struct MacroblockdPlane *plane_ptr) {
    assert(mbmi != NULL);
    (void)mi_row;
    (void)mi_col;
    (void)xd;
    //if (xd->lossless[mbmi->segment_id]) return TX_4X4;

    TxSize tx_size =
        (plane == COMPONENT_LUMA)
            ? (is_inter_block_no_intrabc(mbmi->block_mi.ref_frame[0])
                   ? tx_depth_to_tx_size[0][mbmi->block_mi.sb_type]
                   : tx_depth_to_tx_size[mbmi->tx_depth][mbmi->block_mi.sb_type]) // use max_tx_size
            : av1_get_max_uv_txsize(
                  mbmi->block_mi.sb_type, plane_ptr->subsampling_x, plane_ptr->subsampling_y);
    assert(tx_size < TX_SIZES_ALL);
    if (((plane == COMPONENT_LUMA) && is_inter_block_no_intrabc(mbmi->block_mi.ref_frame[0]) &&
         !mbmi->block_mi.skip)) { // if split tx is used

        const TxSize mb_tx_size =
            tx_depth_to_tx_size[mbmi->tx_depth][mbmi->block_mi.sb_type]; // tx_size
        assert(mb_tx_size < TX_SIZES_ALL);
        tx_size = mb_tx_size;
    }
    // since in case of chrominance or non-square transorm need to convert
    // transform size into transform size in particular direction.
    // for vertical edge, filter direction is horizontal, for horizontal
    // edge, filter direction is vertical.
    tx_size = (VERT_EDGE == edge_dir) ? txsize_horz_map[tx_size] : txsize_vert_map[tx_size];
    return tx_size;
}

// Return TxSize from get_transform_size(), so it is plane and direction
// awared
static TxSize set_lpf_parameters(Av1DeblockingParameters *const params, const uint64_t mode_step,
                                 const PictureControlSet *const pcs_ptr,
                                 const MacroBlockD *const xd, const EdgeDir edge_dir,
                                 const uint32_t x, const uint32_t y, const int32_t plane,
                                 const struct MacroblockdPlane *const plane_ptr) {
    FrameHeader* frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    const LoopFilterInfoN *lfi_n = &pcs_ptr->parent_pcs_ptr->lf_info;

    // reset to initial values
    params->filter_length = 0;

    // no deblocking is required
    const uint32_t width  = plane_ptr->dst.width;
    const uint32_t height = plane_ptr->dst.height;
    if ((width <= x) || (height <= y)) {
        // just return the smallest transform unit size
        return TX_4X4;
    }

    const uint32_t scale_horz = plane_ptr->subsampling_x;
    const uint32_t scale_vert = plane_ptr->subsampling_y;
    // for sub8x8 block, chroma prediction mode is obtained from the bottom/right
    // mi structure of the co-located 8x8 luma block. so for chroma plane, mi_row
    // and mi_col should map to the bottom/right mi structure, i.e, both mi_row
    // and mi_col should be odd number for chroma plane.

    const int32_t mi_row    = scale_vert | ((y << scale_vert) >> MI_SIZE_LOG2);
    const int32_t mi_col    = scale_horz | ((x << scale_horz) >> MI_SIZE_LOG2);
    uint32_t      mi_stride = pcs_ptr->mi_stride;
    const int32_t offset    = mi_row * mi_stride + mi_col;
    ModeInfo **   mi        = (pcs_ptr->mi_grid_base + offset);
    //MbModeInfo **mi = cm->mi_grid_visible + mi_row * cm->mi_stride + mi_col;
    const MbModeInfo *mbmi = &mi[0]->mbmi;

    // If current mbmi is not correctly setup, return an invalid value to stop
    // filtering. One example is that if this tile is not coded, then its mbmi
    // it not set up.
    if (mbmi == NULL) return TX_INVALID;

    const TxSize ts = get_transform_size(xd, mbmi /*mi[0]*/, edge_dir,
                                         mi_row, mi_col, plane, plane_ptr);
    assert(ts < TX_SIZES_ALL);

    {
        const uint32_t coord = (VERT_EDGE == edge_dir) ? (x) : (y);
        const uint32_t transform_masks =
            edge_dir == VERT_EDGE ? tx_size_wide[ts] - 1 : tx_size_high[ts] - 1;
        const int32_t txb_edge = (coord & transform_masks) ? (0) : (1);

        if (!txb_edge) return ts;

        // prepare outer edge parameters. deblock the edge if it's an edge of a TU
        {
            uint32_t curr_level; // Added to address 4x4 problem
            PredictionMode mode = (mbmi->block_mi.mode == INTRA_MODE_4x4)
                                  ? DC_PRED : mbmi->block_mi.mode;
            if (frm_hdr->delta_lf_params.delta_lf_present)
                curr_level = get_filter_level_delta_lf(frm_hdr, edge_dir, plane,
                                                       pcs_ptr->parent_pcs_ptr->curr_delta_lf,
                                                       0 /*segment_id*/,
                                                       mode, mbmi->block_mi.ref_frame[0]);
            else
                curr_level = lfi_n->lvl[plane][0/*segment_id*/][edge_dir]
                             [mbmi->block_mi.ref_frame[0]][mode_lf_lut[mode]];

            const int32_t curr_skipped =
                mbmi->block_mi.skip && is_inter_block_no_intrabc(mbmi->block_mi.ref_frame[0]);

            uint32_t level = curr_level;
            if (coord) {
                {
                    //const ModeInfo *const mi_prev = *(mi - mode_step);
                    const ModeInfo *const   mi_prev_temp = *(mi - mode_step);
                    const MbModeInfo *const mi_prev      = &mi_prev_temp[0].mbmi;
                    //
                    if (mi_prev == NULL) return TX_INVALID;
                    const int32_t pv_row =
                        (VERT_EDGE == edge_dir) ? (mi_row) : (mi_row - (1 << scale_vert));
                    const int32_t pv_col =
                        (VERT_EDGE == edge_dir) ? (mi_col - (1 << scale_horz)) : (mi_col);
                    const TxSize pv_ts =
                        get_transform_size(xd, mi_prev, edge_dir, pv_row, pv_col, plane, plane_ptr);

                    uint32_t pv_lvl;
                    mode = (mi_prev->block_mi.mode == INTRA_MODE_4x4)
                           ? DC_PRED : mi_prev->block_mi.mode;
                    if (frm_hdr->delta_lf_params.delta_lf_present)
                        pv_lvl = get_filter_level_delta_lf(frm_hdr,
                                                           edge_dir, plane,
                                                           pcs_ptr->parent_pcs_ptr->curr_delta_lf,
                                                           0 /*segment_id*/,
                                                           mi_prev->block_mi.mode,
                                                           mi_prev->block_mi.ref_frame[0]);
                    else
                        pv_lvl = lfi_n->lvl[plane][0/*segment_id*/][edge_dir]
                                [mi_prev->block_mi.ref_frame[0]][mode_lf_lut[mode]];

                    const int32_t pv_skip =
                        mi_prev->block_mi.skip &&
                        is_inter_block_no_intrabc(mi_prev->block_mi.ref_frame[0]);

                    const BlockSize bsize = get_plane_block_size(
                        mbmi->block_mi.sb_type, plane_ptr->subsampling_x, plane_ptr->subsampling_y);
                    assert(bsize < BlockSizeS_ALL);
                    const int32_t prediction_masks = edge_dir == VERT_EDGE
                                                         ? block_size_wide[bsize] - 1
                                                         : block_size_high[bsize] - 1;
                    const int32_t pu_edge = !(coord & prediction_masks);
                    // if the current and the previous blocks are skipped,
                    // deblock the edge if the edge belongs to a PU's edge only.
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
                            // No wide filtering for chroma plane
                            if (plane != 0) params->filter_length = 6;
                        }

                        // update the level if the current block is skipped,
                        // but the previous one is not
                        level = (curr_level) ? (curr_level) : (pv_lvl);
                    }
                }
            }
            // prepare common parameters
            if (params->filter_length) {
                const LoopFilterThresh *const limits =
                    pcs_ptr->parent_pcs_ptr->lf_info.lfthr + level;
                params->lim     = limits->lim;
                params->mblim   = limits->mblim;
                params->hev_thr = limits->hev_thr;
            }
        }
    }

    return ts;
}

void eb_av1_filter_block_plane_vert(const PictureControlSet *const pcs_ptr,
                                    const MacroBlockD *const xd, const int32_t plane,
                                    const MacroblockdPlane *const plane_ptr, const uint32_t mi_row,
                                    const uint32_t mi_col) {
    SequenceControlSet *scs_ptr =
        (SequenceControlSet *)pcs_ptr->parent_pcs_ptr->scs_wrapper_ptr->object_ptr;
    EbBool         is_16bit   = scs_ptr->static_config.encoder_bit_depth > 8;
    // TODO
    // when loop_filter_mode = 1, dblk is processed in encdec
    // 16 bit dblk for loop_filter_mode = 1 needs to enabled after 16bit encdec is done
    if (scs_ptr->static_config.is_16bit_pipeline)
        is_16bit = EB_TRUE;
    const int32_t  row_step   = MI_SIZE >> MI_SIZE_LOG2;
    const uint32_t scale_horz = plane_ptr->subsampling_x;
    const uint32_t scale_vert = plane_ptr->subsampling_y;
    uint8_t *const dst_ptr    = plane_ptr->dst.buf;
    const int32_t  dst_stride = plane_ptr->dst.stride;
    const int32_t  y_range    = scs_ptr->seq_header.sb_size == BLOCK_128X128
                                ? (MAX_MIB_SIZE >> scale_vert)
                                : (SB64_MIB_SIZE >> scale_vert);
    const int32_t x_range = scs_ptr->seq_header.sb_size == BLOCK_128X128
                                ? (MAX_MIB_SIZE >> scale_horz)
                                : (SB64_MIB_SIZE >> scale_horz);
    for (int32_t y = 0; y < y_range; y += row_step) {
        uint8_t *p = dst_ptr + ((y * MI_SIZE * dst_stride) << plane_ptr->is_16bit);
        for (int32_t x = 0; x < x_range;) {
            // inner loop always filter vertical edges in a MI block. If MI size
            // is 8x8, it will filter the vertical edge aligned with a 8x8 block.
            // If 4x4 trasnform is used, it will then filter the internal edge
            //  aligned with a 4x4 block
            const uint32_t          curr_x = ((mi_col * MI_SIZE) >> scale_horz) + x * MI_SIZE;
            const uint32_t          curr_y = ((mi_row * MI_SIZE) >> scale_vert) + y * MI_SIZE;
            uint32_t                advance_units;
            TxSize                  tx_size;
            Av1DeblockingParameters params;
            memset(&params, 0, sizeof(params));

            tx_size = set_lpf_parameters(&params,
                                         ((uint64_t)1 << scale_horz),
                                         pcs_ptr,
                                         xd,
                                         VERT_EDGE,
                                         curr_x,
                                         curr_y,
                                         plane,
                                         plane_ptr);
            if (tx_size == TX_INVALID) {
                params.filter_length = 0;
                tx_size              = TX_4X4;
            }

            switch (params.filter_length) {
                // apply 4-tap filtering
            case 4:
                if (is_16bit)
                    aom_highbd_lpf_vertical_4((uint16_t *)(p), //CONVERT_TO_SHORTPTR(p),
                                              dst_stride,
                                              params.mblim,
                                              params.lim,
                                              params.hev_thr,
                                              scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_vertical_4(p, dst_stride, params.mblim, params.lim, params.hev_thr);
                break;
            case 6: // apply 6-tap filter for chroma plane only
                assert(plane != 0);
                if (is_16bit)
                    aom_highbd_lpf_vertical_6((uint16_t *)(p), //CONVERT_TO_SHORTPTR(p),
                                              dst_stride,
                                              params.mblim,
                                              params.lim,
                                              params.hev_thr,
                                              scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_vertical_6(p, dst_stride, params.mblim, params.lim, params.hev_thr);
                break;
                // apply 8-tap filtering
            case 8:
                if (is_16bit)
                    aom_highbd_lpf_vertical_8((uint16_t *)(p), //CONVERT_TO_SHORTPTR(p),
                                              dst_stride,
                                              params.mblim,
                                              params.lim,
                                              params.hev_thr,
                                              scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_vertical_8(p, dst_stride, params.mblim, params.lim, params.hev_thr);
                break;
                // apply 14-tap filtering
            case 14:
                if (is_16bit)
                    aom_highbd_lpf_vertical_14((uint16_t *)(p), //CONVERT_TO_SHORTPTR(p),
                                               dst_stride,
                                               params.mblim,
                                               params.lim,
                                               params.hev_thr,
                                               scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_vertical_14(p, dst_stride, params.mblim, params.lim, params.hev_thr);
                break;
                // no filtering
            default: break;
            }
            // advance the destination pointer
            assert(tx_size < TX_SIZES_ALL);
            advance_units = tx_size_wide_unit[tx_size];
            x += advance_units;
            p += ((advance_units * MI_SIZE) << plane_ptr->is_16bit);
        }
    }
}

void eb_av1_filter_block_plane_horz(const PictureControlSet *const pcs_ptr,
                                    const MacroBlockD *const xd, const int32_t plane,
                                    const MacroblockdPlane *const plane_ptr, const uint32_t mi_row,
                                    const uint32_t mi_col) {
    SequenceControlSet *scs_ptr =
        (SequenceControlSet *)pcs_ptr->parent_pcs_ptr->scs_wrapper_ptr->object_ptr;
    EbBool         is_16bit   = scs_ptr->static_config.encoder_bit_depth > 8;
    // when loop_filter_mode = 1, dblk is processed in encdec
    // 16 bit dblk for loop_filter_mode = 1 needs to enabled after 16bit encdec is done
    if (scs_ptr->static_config.is_16bit_pipeline)
        is_16bit = EB_TRUE;
    const int32_t  col_step   = MI_SIZE >> MI_SIZE_LOG2;
    const uint32_t scale_horz = plane_ptr->subsampling_x;
    const uint32_t scale_vert = plane_ptr->subsampling_y;
    uint8_t *const dst_ptr    = plane_ptr->dst.buf;
    const int32_t  dst_stride = plane_ptr->dst.stride;
    const int32_t  y_range    = scs_ptr->seq_header.sb_size == BLOCK_128X128
                                ? (MAX_MIB_SIZE >> scale_vert)
                                : (SB64_MIB_SIZE >> scale_vert);
    const int32_t x_range = scs_ptr->seq_header.sb_size == BLOCK_128X128
                                ? (MAX_MIB_SIZE >> scale_horz)
                                : (SB64_MIB_SIZE >> scale_horz);
    uint32_t mi_stride = pcs_ptr->mi_stride;
    for (int32_t x = 0; x < x_range; x += col_step) {
        uint8_t *p = dst_ptr + ((x * MI_SIZE) << plane_ptr->is_16bit);
        for (int32_t y = 0; y < y_range;) {
            // inner loop always filter vertical edges in a MI block. If MI size
            // is 8x8, it will first filter the vertical edge aligned with a 8x8
            // block. If 4x4 trasnform is used, it will then filter the internal
            // edge aligned with a 4x4 block
            const uint32_t          curr_x = ((mi_col * MI_SIZE) >> scale_horz) + x * MI_SIZE;
            const uint32_t          curr_y = ((mi_row * MI_SIZE) >> scale_vert) + y * MI_SIZE;
            uint32_t                advance_units;
            TxSize                  tx_size;
            Av1DeblockingParameters params;
            memset(&params, 0, sizeof(params));

            tx_size = set_lpf_parameters(&params,
                                         //(pcs_ptr->parent_pcs_ptr->av1_cm->mi_stride << scale_vert),
                                         (mi_stride << scale_vert),
                                         pcs_ptr,
                                         xd,
                                         HORZ_EDGE,
                                         curr_x,
                                         curr_y,
                                         plane,
                                         plane_ptr);
            if (tx_size == TX_INVALID) {
                params.filter_length = 0;
                tx_size              = TX_4X4;
            }

            switch (params.filter_length) {
                // apply 4-tap filtering
            case 4:
                if (is_16bit)
                    aom_highbd_lpf_horizontal_4((uint16_t *)(p), //CONVERT_TO_SHORTPTR(p),
                                                dst_stride,
                                                params.mblim,
                                                params.lim,
                                                params.hev_thr,
                                                scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_horizontal_4(p, dst_stride, params.mblim, params.lim, params.hev_thr);
                break;
                // apply 6-tap filtering
            case 6:
                assert(plane != 0);
                if (is_16bit)
                    aom_highbd_lpf_horizontal_6((uint16_t *)(p), //CONVERT_TO_SHORTPTR(p),
                                                dst_stride,
                                                params.mblim,
                                                params.lim,
                                                params.hev_thr,
                                                scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_horizontal_6(p, dst_stride, params.mblim, params.lim, params.hev_thr);
                break;
                // apply 8-tap filtering
            case 8:
                if (is_16bit)
                    aom_highbd_lpf_horizontal_8((uint16_t *)(p), //CONVERT_TO_SHORTPTR(p),
                                                dst_stride,
                                                params.mblim,
                                                params.lim,
                                                params.hev_thr,
                                                scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_horizontal_8(p, dst_stride, params.mblim, params.lim, params.hev_thr);
                break;
                // apply 14-tap filtering
            case 14:
                if (is_16bit)
                    aom_highbd_lpf_horizontal_14((uint16_t *)(p), //CONVERT_TO_SHORTPTR(p),
                                                 dst_stride,
                                                 params.mblim,
                                                 params.lim,
                                                 params.hev_thr,
                                                 scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_horizontal_14(p, dst_stride, params.mblim, params.lim, params.hev_thr);
                break;
                // no filtering
            default: break;
            }

            // advance the destination pointer
            assert(tx_size < TX_SIZES_ALL);
            advance_units = tx_size_high_unit[tx_size];
            y += advance_units;
            p += ((advance_units * dst_stride * MI_SIZE) << plane_ptr->is_16bit);
        }
    }
}

// New function to filter each sb (64x64)
void loop_filter_sb(EbPictureBufferDesc *frame_buffer, //reconpicture,
                    //Yv12BufferConfig *frame_buffer,
                    PictureControlSet *pcs_ptr, MacroBlockD *xd, int32_t mi_row, int32_t mi_col,
                    int32_t plane_start, int32_t plane_end, uint8_t last_col) {
    FrameHeader *           frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    struct MacroblockdPlane pd[3];
    int32_t                 plane;

    pd[0].subsampling_x = 0;
    pd[0].subsampling_y = 0;
    pd[0].plane_type    = PLANE_TYPE_Y;
    pd[0].is_16bit      = frame_buffer->bit_depth > 8;
    pd[1].subsampling_x = 1;
    pd[1].subsampling_y = 1;
    pd[1].plane_type    = PLANE_TYPE_UV;
    pd[1].is_16bit      = frame_buffer->bit_depth > 8;
    pd[2].subsampling_x = 1;
    pd[2].subsampling_y = 1;
    pd[2].plane_type    = PLANE_TYPE_UV;
    pd[2].is_16bit      = frame_buffer->bit_depth > 8;

    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline)
        pd[0].is_16bit = pd[1].is_16bit = pd[2].is_16bit = EB_TRUE;

    for (plane = plane_start; plane < plane_end; plane++) {
        if (plane == 0 && !(frm_hdr->loop_filter_params.filter_level[0]) &&
            !(frm_hdr->loop_filter_params.filter_level[1]))
            break;
        else if (plane == 1 && !(frm_hdr->loop_filter_params.filter_level_u))
            continue;
        else if (plane == 2 && !(frm_hdr->loop_filter_params.filter_level_v))
            continue;

        if (frm_hdr->loop_filter_params.combine_vert_horz_lf) {
            // filter all vertical and horizontal edges in every 64x64 super block
            // filter vertical edges
            eb_av1_setup_dst_planes(pd,
                                    pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.sb_size,
                                    frame_buffer,
                                    mi_row,
                                    mi_col,
                                    plane,
                                    plane + 1);
            eb_av1_filter_block_plane_vert(pcs_ptr, xd, plane, &pd[plane], mi_row, mi_col);
            // filter horizontal edges
            int32_t max_mib_size =
                pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.sb_size == BLOCK_128X128
                    ? MAX_MIB_SIZE
                    : SB64_MIB_SIZE;

            if (mi_col - max_mib_size >= 0) {
                eb_av1_setup_dst_planes(pd,
                                        pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.sb_size,
                                        frame_buffer,
                                        mi_row,
                                        mi_col - max_mib_size,
                                        plane,
                                        plane + 1);
                eb_av1_filter_block_plane_horz(
                    pcs_ptr, xd, plane, &pd[plane], mi_row, mi_col - max_mib_size);
            }
            // Filter the horizontal edges of the last sb in each row
            if (last_col) {
                eb_av1_setup_dst_planes(pd,
                                        pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.sb_size,
                                        frame_buffer,
                                        mi_row,
                                        mi_col,
                                        plane,
                                        plane + 1);
                eb_av1_filter_block_plane_horz(pcs_ptr, xd, plane, &pd[plane], mi_row, mi_col);
            }
        } else {
            // filter all vertical edges in every 64x64 super block
            eb_av1_setup_dst_planes(pd,
                                    pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.sb_size,
                                    frame_buffer,
                                    mi_row,
                                    mi_col,
                                    plane,
                                    plane + 1);

            eb_av1_filter_block_plane_vert(pcs_ptr, xd, plane, &pd[plane], mi_row, mi_col);

            // filter all horizontal edges in every 64x64 super block
            eb_av1_setup_dst_planes(pd,
                                    pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.sb_size,
                                    frame_buffer,
                                    mi_row,
                                    mi_col,
                                    plane,
                                    plane + 1);
            eb_av1_filter_block_plane_horz(pcs_ptr, xd, plane, &pd[plane], mi_row, mi_col);
        }
    }
}

void eb_av1_loop_filter_frame(EbPictureBufferDesc *frame_buffer, PictureControlSet *pcs_ptr,
                              int32_t plane_start, int32_t plane_end) {
    SequenceControlSet *scs_ptr =
        (SequenceControlSet *)pcs_ptr->parent_pcs_ptr->scs_wrapper_ptr->object_ptr;
    //SuperBlock                     *sb_ptr;
    //uint16_t                                   sb_index;
    uint8_t  sb_size_log2 = (uint8_t)eb_log2f(scs_ptr->sb_size_pix);
    uint32_t x_sb_index;
    uint32_t y_sb_index;
    uint32_t sb_origin_x;
    uint32_t sb_origin_y;
    EbBool   end_of_row_flag;

    uint32_t pic_width_in_sb =
        (pcs_ptr->parent_pcs_ptr->aligned_width + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;
    uint32_t picture_height_in_sb =
        (pcs_ptr->parent_pcs_ptr->aligned_height + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;

    eb_av1_loop_filter_frame_init(&pcs_ptr->parent_pcs_ptr->frm_hdr,
                                  &pcs_ptr->parent_pcs_ptr->lf_info,
                                  plane_start,
                                  plane_end);

    for (y_sb_index = 0; y_sb_index < picture_height_in_sb; ++y_sb_index) {
        for (x_sb_index = 0; x_sb_index < pic_width_in_sb; ++x_sb_index) {
            //sb_index        = (uint16_t)(y_sb_index * pic_width_in_sb + x_sb_index);
            //sb_ptr          = pcs_ptr->sb_ptr_array[sb_index];
            sb_origin_x     = x_sb_index << sb_size_log2;
            sb_origin_y     = y_sb_index << sb_size_log2;
            end_of_row_flag = (x_sb_index == pic_width_in_sb - 1) ? EB_TRUE : EB_FALSE;

            loop_filter_sb(frame_buffer,
                           pcs_ptr,
                           NULL,
                           sb_origin_y >> 2,
                           sb_origin_x >> 2,
                           plane_start,
                           plane_end,
                           end_of_row_flag);
        }
    }
}
extern int16_t eb_av1_ac_quant_q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);

void eb_copy_buffer(EbPictureBufferDesc *srcBuffer, EbPictureBufferDesc *dstBuffer,
                    PictureControlSet *pcs_ptr, uint8_t plane) {
    EbBool is_16bit =
        (EbBool)(pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline)
        is_16bit = EB_TRUE;
    dstBuffer->origin_x    = srcBuffer->origin_x;
    dstBuffer->origin_y    = srcBuffer->origin_y;
    dstBuffer->width       = srcBuffer->width;
    dstBuffer->height      = srcBuffer->height;
    dstBuffer->max_width   = srcBuffer->max_width;
    dstBuffer->max_height  = srcBuffer->max_height;
    dstBuffer->bit_depth   = srcBuffer->bit_depth;
    dstBuffer->luma_size   = srcBuffer->luma_size;
    dstBuffer->chroma_size = srcBuffer->chroma_size;
    dstBuffer->packed_flag = srcBuffer->packed_flag;

    uint32_t luma_buffer_offset = (srcBuffer->origin_x + srcBuffer->origin_y * srcBuffer->stride_y)
                                  << is_16bit;
    uint16_t luma_width = (uint16_t)(srcBuffer->width - pcs_ptr->parent_pcs_ptr->scs_ptr->pad_right)
                          << is_16bit;
    uint16_t luma_height =
        (uint16_t)(srcBuffer->height - pcs_ptr->parent_pcs_ptr->scs_ptr->pad_bottom);
    uint16_t chroma_width = (luma_width >> 1);
    if (plane == 0) {
        uint16_t stride_y = srcBuffer->stride_y << is_16bit;

        dstBuffer->stride_y         = srcBuffer->stride_y;
        dstBuffer->stride_bit_inc_y = srcBuffer->stride_bit_inc_y;

        for (int32_t input_row_index = 0; input_row_index < luma_height; input_row_index++) {
            eb_memcpy((dstBuffer->buffer_y + luma_buffer_offset + stride_y * input_row_index),
                      (srcBuffer->buffer_y + luma_buffer_offset + stride_y * input_row_index),
                      luma_width);
        }
    } else if (plane == 1) {
        uint16_t stride_cb           = srcBuffer->stride_cb << is_16bit;
        dstBuffer->stride_cb         = srcBuffer->stride_cb;
        dstBuffer->stride_bit_inc_cb = srcBuffer->stride_bit_inc_cb;

        uint32_t chroma_buffer_offset =
            (srcBuffer->origin_x / 2 + srcBuffer->origin_y / 2 * srcBuffer->stride_cb) << is_16bit;

        for (int32_t input_row_index = 0; input_row_index < luma_height / 2; input_row_index++) {
            eb_memcpy((dstBuffer->buffer_cb + chroma_buffer_offset + stride_cb * input_row_index),
                      (srcBuffer->buffer_cb + chroma_buffer_offset + stride_cb * input_row_index),
                      chroma_width);
        }
    } else if (plane == 2) {
        uint16_t stride_cr = srcBuffer->stride_cr << is_16bit;

        dstBuffer->stride_cr         = srcBuffer->stride_cr;
        dstBuffer->stride_bit_inc_cr = srcBuffer->stride_bit_inc_cr;

        uint32_t chroma_buffer_offset =
            (srcBuffer->origin_x / 2 + srcBuffer->origin_y / 2 * srcBuffer->stride_cr) << is_16bit;

        for (int32_t input_row_index = 0; input_row_index< luma_height/2; input_row_index++) {
            eb_memcpy((dstBuffer->buffer_cr + chroma_buffer_offset + stride_cr * input_row_index),
                      (srcBuffer->buffer_cr + chroma_buffer_offset + stride_cr * input_row_index),
                      chroma_width);
        }
    }
}

//int32_t av1_get_max_filter_level(const Av1Comp *cpi) {
//    if (cpi->oxcf.pass == 2) {
//        return cpi->twopass.section_intra_rating > 8 ? MAX_LOOP_FILTER * 3 / 4
//            : MAX_LOOP_FILTER;
//    }
//    else {
//        return MAX_LOOP_FILTER;
//    }
//}

uint64_t picture_sse_calculations(PictureControlSet *pcs_ptr, EbPictureBufferDesc *recon_ptr,
                                  int32_t plane)

{
    SequenceControlSet *scs_ptr  = pcs_ptr->parent_pcs_ptr->scs_ptr;
    EbBool              is_16bit = scs_ptr->static_config.is_16bit_pipeline ||
                                   (scs_ptr->static_config.encoder_bit_depth > EB_8BIT);

    const uint32_t ss_x = scs_ptr->subsampling_x;
    const uint32_t ss_y = scs_ptr->subsampling_y;

    if (!is_16bit) {
        EbPictureBufferDesc *input_picture_ptr =
            (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;

        uint32_t column_index;
        uint32_t row_index           = 0;
        uint64_t residual_distortion = 0;
        EbByte   input_buffer;
        EbByte   recon_coeff_buffer;
        if (plane == 0) {
            recon_coeff_buffer = &(
                (recon_ptr
                     ->buffer_y)[recon_ptr->origin_x + recon_ptr->origin_y * recon_ptr->stride_y]);
            input_buffer =
                &((input_picture_ptr
                       ->buffer_y)[input_picture_ptr->origin_x +
                                   input_picture_ptr->origin_y * input_picture_ptr->stride_y]);

            residual_distortion = 0;

            while (row_index < input_picture_ptr->height) {
                column_index = 0;
                while (column_index < input_picture_ptr->width) {
                    residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                        (recon_coeff_buffer[column_index]));
                    ++column_index;
                }
                input_buffer += input_picture_ptr->stride_y;
                recon_coeff_buffer += recon_ptr->stride_y;
                ++row_index;
            }

            return residual_distortion;
        }

        else if (plane == 1) {
            recon_coeff_buffer =
                &((recon_ptr->buffer_cb)[recon_ptr->origin_x / 2 +
                                         recon_ptr->origin_y / 2 * recon_ptr->stride_cb]);
            input_buffer = &((input_picture_ptr->buffer_cb)[input_picture_ptr->origin_x / 2 +
                                                            input_picture_ptr->origin_y / 2 *
                                                                input_picture_ptr->stride_cb]);

            residual_distortion = 0;
            row_index           = 0;
            while (row_index < (uint32_t)(input_picture_ptr->height >> ss_y)) {
                column_index = 0;
                while (column_index < (uint32_t)(input_picture_ptr->width >> ss_x)) {
                    residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                        (recon_coeff_buffer[column_index]));
                    ++column_index;
                }

                input_buffer += input_picture_ptr->stride_cb;
                recon_coeff_buffer += recon_ptr->stride_cb;
                ++row_index;
            }

            return residual_distortion;
        } else if (plane == 2) {
            recon_coeff_buffer =
                &((recon_ptr->buffer_cr)[recon_ptr->origin_x / 2 +
                                         recon_ptr->origin_y / 2 * recon_ptr->stride_cr]);
            input_buffer        = &((input_picture_ptr->buffer_cr)[input_picture_ptr->origin_x / 2 +
                                                            input_picture_ptr->origin_y / 2 *
                                                                input_picture_ptr->stride_cr]);

            residual_distortion = 0;
            row_index           = 0;

            while (row_index < (uint32_t)(input_picture_ptr->height >> ss_y)) {
                column_index = 0;
                while (column_index < (uint32_t)(input_picture_ptr->width >> ss_x)) {
                    residual_distortion += (int64_t)SQR((int64_t)(input_buffer[column_index]) -
                                                        (recon_coeff_buffer[column_index]));
                    ++column_index;
                }

                input_buffer += input_picture_ptr->stride_cr;
                recon_coeff_buffer += recon_ptr->stride_cr;
                ++row_index;
            }

            return residual_distortion;
        }
        return 0;
    } else {
        EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc *)pcs_ptr->input_frame16bit;

        uint32_t  column_index;
        uint32_t  row_index           = 0;
        uint64_t  residual_distortion = 0;
        uint16_t *input_buffer;
        uint16_t *recon_coeff_buffer;
        if (plane == 0) {
            recon_coeff_buffer = (uint16_t *)&(
                (recon_ptr
                     ->buffer_y)[(recon_ptr->origin_x + recon_ptr->origin_y * recon_ptr->stride_y)
                                 << is_16bit]);
            input_buffer = (uint16_t *)&(
                (input_picture_ptr
                     ->buffer_y)[(input_picture_ptr->origin_x +
                                  input_picture_ptr->origin_y * input_picture_ptr->stride_y)
                                 << is_16bit]);

            residual_distortion = 0;

            while (row_index < input_picture_ptr->height) {
                column_index = 0;
                while (column_index < input_picture_ptr->width) {
                    residual_distortion +=
                        (int64_t)SQR(((int64_t)input_buffer[column_index]) -
                                     (int64_t)(recon_coeff_buffer[column_index]));
                    ++column_index;
                }

                input_buffer += input_picture_ptr->stride_y;
                recon_coeff_buffer += recon_ptr->stride_y;
                ++row_index;
            }

            return residual_distortion;
        }

        else if (plane == 1) {
            recon_coeff_buffer = (uint16_t *)&(
                (recon_ptr->buffer_cb)[(recon_ptr->origin_x / 2 +
                                        recon_ptr->origin_y / 2 * recon_ptr->stride_cb)
                                       << is_16bit]);
            input_buffer = (uint16_t *)&(
                (input_picture_ptr
                     ->buffer_cb)[(input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb)
                                  << is_16bit]);

            residual_distortion = 0;
            row_index           = 0;
            while (row_index < (uint32_t)(input_picture_ptr->height >> ss_y)) {
                column_index = 0;
                while (column_index < (uint32_t)(input_picture_ptr->width >> ss_x)) {
                    residual_distortion +=
                        (int64_t)SQR(((int64_t)input_buffer[column_index]) -
                                     (int64_t)(recon_coeff_buffer[column_index]));
                    ++column_index;
                }

                input_buffer += input_picture_ptr->stride_cb;
                recon_coeff_buffer += recon_ptr->stride_cb;
                ++row_index;
            }

            return residual_distortion;
        } else if (plane == 2) {
            recon_coeff_buffer = (uint16_t *)&(
                (recon_ptr->buffer_cr)[(recon_ptr->origin_x / 2 +
                                        recon_ptr->origin_y / 2 * recon_ptr->stride_cr)
                                       << is_16bit]);
            input_buffer = (uint16_t *)&(
                (input_picture_ptr
                     ->buffer_cr)[(input_picture_ptr->origin_x / 2 +
                                   input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr)
                                  << is_16bit]);
            residual_distortion = 0;
            row_index           = 0;

            while (row_index < (uint32_t)(input_picture_ptr->height >> ss_y)) {
                column_index = 0;
                while (column_index < (uint32_t)(input_picture_ptr->width >> ss_x)) {
                    residual_distortion +=
                        (int64_t)SQR(((int64_t)input_buffer[column_index]) -
                                     (int64_t)(recon_coeff_buffer[column_index]));
                    ++column_index;
                }

                input_buffer += input_picture_ptr->stride_cr;
                recon_coeff_buffer += recon_ptr->stride_cr;
                ++row_index;
            }

            return residual_distortion;
        }

        return 0;
    }
}

static int64_t try_filter_frame(
    //const Yv12BufferConfig *sd,
    //Av1Comp *const cpi,
    const EbPictureBufferDesc *sd, EbPictureBufferDesc *temp_lf_recon_buffer,
    PictureControlSet *pcs_ptr, int32_t filt_level, int32_t partial_frame, int32_t plane,
    int32_t dir) {
    (void)sd;
    (void)partial_frame;
    (void)sd;
    int64_t      filt_err;
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    assert(plane >= 0 && plane <= 2);
    int32_t filter_level[2] = {filt_level, filt_level};
    if (plane == 0 && dir == 0) filter_level[1] = frm_hdr->loop_filter_params.filter_level[1];
    if (plane == 0 && dir == 1) filter_level[0] = frm_hdr->loop_filter_params.filter_level[0];

    EbBool is_16bit =
        (EbBool)(pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    EbPictureBufferDesc *recon_buffer =
        is_16bit ? pcs_ptr->recon_picture16bit_ptr : pcs_ptr->recon_picture_ptr;
    if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE) {
        //get the 16bit form of the input SB
        if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline || is_16bit)
            recon_buffer = ((EbReferenceObject *)
                                pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                               ->reference_picture16bit;
        else
            recon_buffer = ((EbReferenceObject *)
                                pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                               ->reference_picture;
    } else { // non ref pictures
        recon_buffer = (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline ||
                        is_16bit) ? pcs_ptr->recon_picture16bit_ptr : pcs_ptr->recon_picture_ptr;
    }

    // set base filters for use of get_filter_level when in DELTA_Q_LF mode
    switch (plane) {
    case 0:
        frm_hdr->loop_filter_params.filter_level[0] = filter_level[0];
        frm_hdr->loop_filter_params.filter_level[1] = filter_level[1];
        break;
    case 1: frm_hdr->loop_filter_params.filter_level_u = filter_level[0]; break;
    case 2: frm_hdr->loop_filter_params.filter_level_v = filter_level[0]; break;
    }

    eb_av1_loop_filter_frame(recon_buffer, pcs_ptr, plane, plane + 1);

    filt_err = picture_sse_calculations(pcs_ptr, recon_buffer, plane);

    // Re-instate the unfiltered frame
    eb_copy_buffer(temp_lf_recon_buffer /*cpi->last_frame_uf*/,
                   recon_buffer /*cm->frame_to_show*/,
                   pcs_ptr,
                   (uint8_t)plane);

    return filt_err;
}
static int32_t search_filter_level(
    //const Yv12BufferConfig *sd, Av1Comp *cpi,
    EbPictureBufferDesc *sd, // source
    EbPictureBufferDesc *temp_lf_recon_buffer, PictureControlSet *pcs_ptr, int32_t partial_frame,
    const int32_t *last_frame_filter_level, double *best_cost_ret, int32_t plane, int32_t dir) {
    const int32_t min_filter_level = 0;
    const int32_t max_filter_level = MAX_LOOP_FILTER; // av1_get_max_filter_level(cpi);
    int32_t       filt_direction   = 0;
    int64_t       best_err;
    int32_t       filt_best;
    FrameHeader * frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    //Macroblock *x = &cpi->td.mb;

    // Start the search at the previous frame filter level unless it is now out of
    // range.
    int32_t lvl;
    switch (plane) {
    case 0: lvl = last_frame_filter_level[dir]; break;
    case 1: lvl = last_frame_filter_level[2]; break;
    case 2: lvl = last_frame_filter_level[3]; break;
    default: assert(plane >= 0 && plane <= 2); return 0;
    }
    int32_t filt_mid    = clamp(lvl, min_filter_level, max_filter_level);
    int32_t filter_step = filt_mid < 16 ? 4 : filt_mid / 4;

    EbBool is_16bit =
        (EbBool)(pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
    EbPictureBufferDesc *recon_buffer =
        is_16bit ? pcs_ptr->recon_picture16bit_ptr : pcs_ptr->recon_picture_ptr;

    if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE) {
        //get the 16bit form of the input SB
        if (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline || is_16bit)
            recon_buffer = ((EbReferenceObject *)
                                pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                               ->reference_picture16bit;
        else
            recon_buffer = ((EbReferenceObject *)
                                pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                               ->reference_picture;
    } else { // non ref pictures
        recon_buffer = (pcs_ptr->parent_pcs_ptr->scs_ptr->static_config.is_16bit_pipeline ||
                        is_16bit) ? pcs_ptr->recon_picture16bit_ptr : pcs_ptr->recon_picture_ptr;
    }
    // Sum squared error at each filter level
    int64_t ss_err[MAX_LOOP_FILTER + 1];

    // Set each entry to -1
    memset(ss_err, 0xFF, sizeof(ss_err));
    // make a copy of recon_buffer
    eb_copy_buffer(recon_buffer /*cm->frame_to_show*/,
                   temp_lf_recon_buffer /*&cpi->last_frame_uf*/,
                   pcs_ptr,
                   (uint8_t)plane);

    best_err =
        try_filter_frame(sd, temp_lf_recon_buffer, pcs_ptr, filt_mid, partial_frame, plane, dir);
    filt_best        = filt_mid;
    ss_err[filt_mid] = best_err;

    if (pcs_ptr->parent_pcs_ptr->loop_filter_mode <= 2) {
        filter_step             = 2;
        const int32_t filt_high = AOMMIN(filt_mid + filter_step, max_filter_level);
        const int32_t filt_low  = AOMMAX(filt_mid - filter_step, min_filter_level);

        // Bias against raising loop filter in favor of lowering it.
        int64_t bias = (best_err >> (15 - (filt_mid / 8))) * filter_step;

        //if ((cpi->oxcf.pass == 2) && (cpi->twopass.section_intra_rating < 20))
        //    bias = (bias * cpi->twopass.section_intra_rating) / 20;

        // yx, bias less for large block size
        if (frm_hdr->tx_mode != ONLY_4X4) bias >>= 1;

        if (filt_direction <= 0 && filt_low != filt_mid) {
            // Get Low filter error score
            if (ss_err[filt_low] < 0) {
                ss_err[filt_low] = try_filter_frame(
                    sd, temp_lf_recon_buffer, pcs_ptr, filt_low, partial_frame, plane, dir);
            }
            // If value is close to the best so far then bias towards a lower loop
            // filter value.
            if (ss_err[filt_low] < (best_err + bias)) {
                // Was it actually better than the previous best?
                if (ss_err[filt_low] < best_err) best_err = ss_err[filt_low];
                filt_best = filt_low;
            }
        }

        // Now look at filt_high
        if (filt_direction >= 0 && filt_high != filt_mid) {
            if (ss_err[filt_high] < 0) {
                ss_err[filt_high] = try_filter_frame(
                    sd, temp_lf_recon_buffer, pcs_ptr, filt_high, partial_frame, plane, dir);
            }
            // If value is significantly better than previous best, bias added against
            // raising filter value
            if (ss_err[filt_high] < (best_err - bias)) {
                best_err  = ss_err[filt_high];
                filt_best = filt_high;
            }
        }
    } else {
        while (filter_step > 0) {
            const int32_t filt_high = AOMMIN(filt_mid + filter_step, max_filter_level);
            const int32_t filt_low  = AOMMAX(filt_mid - filter_step, min_filter_level);

            // Bias against raising loop filter in favor of lowering it.
            int64_t bias = (best_err >> (15 - (filt_mid / 8))) * filter_step;

            //if ((cpi->oxcf.pass == 2) && (cpi->twopass.section_intra_rating < 20))
            //    bias = (bias * cpi->twopass.section_intra_rating) / 20;

            // yx, bias less for large block size
            if (frm_hdr->tx_mode != ONLY_4X4) bias >>= 1;

            if (filt_direction <= 0 && filt_low != filt_mid) {
                // Get Low filter error score
                if (ss_err[filt_low] < 0) {
                    ss_err[filt_low] = try_filter_frame(
                        sd, temp_lf_recon_buffer, pcs_ptr, filt_low, partial_frame, plane, dir);
                }
                // If value is close to the best so far then bias towards a lower loop
                // filter value.
                if (ss_err[filt_low] < (best_err + bias)) {
                    // Was it actually better than the previous best?
                    if (ss_err[filt_low] < best_err) best_err = ss_err[filt_low];
                    filt_best = filt_low;
                }
            }

            // Now look at filt_high
            if (filt_direction >= 0 && filt_high != filt_mid) {
                if (ss_err[filt_high] < 0) {
                    ss_err[filt_high] = try_filter_frame(
                        sd, temp_lf_recon_buffer, pcs_ptr, filt_high, partial_frame, plane, dir);
                }
                // If value is significantly better than previous best, bias added against
                // raising filter value
                if (ss_err[filt_high] < (best_err - bias)) {
                    best_err  = ss_err[filt_high];
                    filt_best = filt_high;
                }
            }

            // Half the step distance if the best filter value was the same as last time
            if (filt_best == filt_mid) {
                filter_step /= 2;
                filt_direction = 0;
            } else {
                filt_direction = (filt_best < filt_mid) ? -1 : 1;
                filt_mid       = filt_best;
            }
        }
    }
    // Update best error
    best_err = ss_err[filt_best];

    if (best_cost_ret) *best_cost_ret = (double)best_err; //RDCOST_DBL(x->rdmult, 0, best_err);
    return filt_best;
}

void eb_av1_pick_filter_level(DlfContext *         context_ptr,
                              EbPictureBufferDesc *srcBuffer, // source input
                              PictureControlSet *pcs_ptr, LpfPickMethod method) {
    SequenceControlSet *scs_ptr =
        (SequenceControlSet *)pcs_ptr->parent_pcs_ptr->scs_wrapper_ptr->object_ptr;
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;

    const int32_t num_planes = 3;
    (void)srcBuffer;
    struct LoopFilter *const lf = &frm_hdr->loop_filter_params;
    lf->sharpness_level         = frm_hdr->frame_type == KEY_FRAME ? 0 : 0;

    if (method == LPF_PICK_MINIMAL_LPF) {
        lf->filter_level[0] = 0;
        lf->filter_level[1] = 0;
    } else if (method >= LPF_PICK_FROM_Q) {
        const int32_t min_filter_level = 0;
        const int32_t max_filter_level = MAX_LOOP_FILTER; // av1_get_max_filter_level(cpi);
        const int32_t q                = eb_av1_ac_quant_q3(frm_hdr->quantization_params.base_q_idx,
                                             0,
                                             (AomBitDepth)scs_ptr->static_config.encoder_bit_depth);
        // These values were determined by linear fitting the result of the
        // searched level for 8 bit depth:
        // Keyframes: filt_guess = q * 0.06699 - 1.60817
        // Other frames: filt_guess = q * 0.02295 + 2.48225
        //
        // And high bit depth separately:
        // filt_guess = q * 0.316206 + 3.87252
        int32_t filt_guess;
        switch (scs_ptr->static_config.encoder_bit_depth) {
        case EB_8BIT:
            filt_guess = (frm_hdr->frame_type == KEY_FRAME)
                             ? ROUND_POWER_OF_TWO(q * 17563 - 421574, 18)
                             : ROUND_POWER_OF_TWO(q * 6017 + 650707, 18);
            break;
        case EB_10BIT: filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 4060632, 20); break;
        case EB_12BIT: filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 16242526, 22); break;
        default:
            assert(0 &&
                   "bit_depth should be AOM_BITS_8, AOM_BITS_10 "
                   "or AOM_BITS_12");
            return;
        }
        if (scs_ptr->static_config.encoder_bit_depth != EB_8BIT && frm_hdr->frame_type == KEY_FRAME)
            filt_guess -= 4;

        filt_guess = filt_guess > 2 ? filt_guess - 2 : filt_guess > 1 ? filt_guess - 1 : filt_guess;
        int32_t filt_guess_chroma = filt_guess > 1 ? filt_guess / 2 : filt_guess;

        // TODO(chengchen): retrain the model for Y, U, V filter levels
        lf->filter_level[0] = clamp(filt_guess, min_filter_level, max_filter_level);
        lf->filter_level[1] = clamp(filt_guess, min_filter_level, max_filter_level);
        lf->filter_level_u  = clamp(filt_guess_chroma, min_filter_level, max_filter_level);
        lf->filter_level_v  = clamp(filt_guess_chroma, min_filter_level, max_filter_level);
    } else {
        const int32_t last_frame_filter_level[4] = {
            lf->filter_level[0], lf->filter_level[1], lf->filter_level_u, lf->filter_level_v};
        EbPictureBufferDesc *temp_lf_recon_buffer =
            (scs_ptr->static_config.is_16bit_pipeline ||
             scs_ptr->static_config.encoder_bit_depth != EB_8BIT)
                ? context_ptr->temp_lf_recon_picture16bit_ptr
                : context_ptr->temp_lf_recon_picture_ptr;

        lf->filter_level[0] = lf->filter_level[1] =
            search_filter_level(srcBuffer,
                                temp_lf_recon_buffer,
                                pcs_ptr,
                                method == LPF_PICK_FROM_SUBIMAGE,
                                last_frame_filter_level,
                                NULL,
                                0,
                                2);

        if (num_planes > 1) {
            lf->filter_level_u = search_filter_level(srcBuffer,
                                                     temp_lf_recon_buffer,
                                                     pcs_ptr,
                                                     method == LPF_PICK_FROM_SUBIMAGE,
                                                     last_frame_filter_level,
                                                     NULL,
                                                     1,
                                                     0);
            lf->filter_level_v = search_filter_level(srcBuffer,
                                                     temp_lf_recon_buffer,
                                                     pcs_ptr,
                                                     method == LPF_PICK_FROM_SUBIMAGE,
                                                     last_frame_filter_level,
                                                     NULL,
                                                     2,
                                                     0);
        }
    }
}
