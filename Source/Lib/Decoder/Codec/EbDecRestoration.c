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

#include "EbDefinitions.h"
#include "EbDecHandle.h"
#include "EbDecUtils.h"

#include "EbDecInverseQuantize.h"
#include "EbDecProcessFrame.h"
#include "EbRestoration.h"
#include "EbDecRestoration.h"

void save_tile_row_boundary_lines(uint8_t *src, int32_t src_stride,
    int32_t src_width, int32_t src_height, int32_t use_highbd, int32_t plane,
    Av1Common *cm, int32_t after_cdef, RestorationStripeBoundaries *boundaries);

void dec_av1_loop_restoration_filter_frame(EbDecHandle *dec_handle, int optimized_lr)
{
    assert(!dec_handle->frame_header.all_lossless);

    FrameHeader *frame_header = &dec_handle->frame_header;
    LRCtxt *lr_ctxt = (LRCtxt *)dec_handle->pv_lr_ctxt;
    MasterFrameBuf *master_frame_buf = &dec_handle->master_frame_buf;
    CurFrameBuf    *frame_buf = &master_frame_buf->cur_frame_bufs[0];
    RestorationTileLimits tile_limit;
    AV1PixelRect tile_rect;
    EbPictureBufferDesc *cur_pic_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;
    RestorationUnitInfo *lr_unit;

    lr_ctxt->lr_unit[AOM_PLANE_Y] = frame_buf->lr_unit[AOM_PLANE_Y];
    lr_ctxt->lr_unit[AOM_PLANE_U] = frame_buf->lr_unit[AOM_PLANE_U];
    lr_ctxt->lr_unit[AOM_PLANE_V] = frame_buf->lr_unit[AOM_PLANE_V];

    int num_plane = av1_num_planes(&dec_handle->seq_header.color_config);
    int use_highbd = (dec_handle->seq_header.color_config.bit_depth > 8);
    int bit_depth = dec_handle->seq_header.color_config.bit_depth;
    int h = 0, w = 0, x, y, unit_row, unit_col;
    int src_stride, dst_stride, tile_stripe0 = 0;
    uint8_t *src, *dst;

    for (int plane = 0; plane < num_plane; plane++)
    {
        LRParams *lr_params = &frame_header->lr_params[plane];
        int ext_size = lr_params->loop_restoration_size * 3 / 2;
        int is_uv = plane > 0;
        int sx = 0, sy = 0;

        if (plane) {
            sx = dec_handle->seq_header.color_config.subsampling_x;
            sy = dec_handle->seq_header.color_config.subsampling_y;
        }

        // src points to frame start
        derive_blk_pointers(cur_pic_buf, plane, 0, 0, (void *)&src,
                            &src_stride, sx, sy);

        dst = lr_ctxt->dst;
        dst_stride = lr_ctxt->dst_stride;

        tile_rect = whole_frame_rect(&dec_handle->frame_header.frame_size,
            dec_handle->seq_header.color_config.subsampling_x,
            dec_handle->seq_header.color_config.subsampling_y, is_uv);

        int tile_h = tile_rect.bottom - tile_rect.top;
        int tile_w = tile_rect.right - tile_rect.left;

        if (lr_params->frame_restoration_type == RESTORE_NONE)
            continue;

        for (y = 0, unit_row = 0; y < tile_h; y += h, unit_row++)
        {
            int remaining_h = tile_h - y;
            h = (remaining_h < ext_size) ? remaining_h :
                                           lr_params->loop_restoration_size;
            tile_limit.v_start = tile_rect.top + y;
            tile_limit.v_end   = tile_rect.top + y + h;
            assert(tile_limit.v_end <= tile_rect.bottom);

            // Offset the tile upwards to align with the restoration processing stripe
            const int voffset = RESTORATION_UNIT_OFFSET >> sy;
            tile_limit.v_start = AOMMAX(tile_rect.top, tile_limit.v_start - voffset);
            if (tile_limit.v_end < tile_rect.bottom) tile_limit.v_end -= voffset;

            for (x = 0, unit_col = 0; x < tile_w; x += w, unit_col++)
            {
                int remaining_w = tile_w - x;
                w = (remaining_w < ext_size) ? remaining_w :
                                               lr_params->loop_restoration_size;

                tile_limit.h_start = tile_rect.left + x;
                tile_limit.h_end = tile_rect.left + x + w;

                lr_unit = lr_ctxt->lr_unit[plane] +
                    unit_row * lr_ctxt->lr_stride[plane] + unit_col;

                if (!use_highbd)
                    eb_av1_loop_restoration_filter_unit(1, &tile_limit, lr_unit,
                        &lr_ctxt->boundaries[plane], lr_ctxt->rlbs, &tile_rect,
                        tile_stripe0, sx, sy, use_highbd, bit_depth, src,
                        src_stride, dst, dst_stride, lr_ctxt->rst_tmpbuf, optimized_lr);
                else
                    eb_av1_loop_restoration_filter_unit(1, &tile_limit, lr_unit,
                        &lr_ctxt->boundaries[plane], lr_ctxt->rlbs, &tile_rect,
                        tile_stripe0, sx, sy, use_highbd, bit_depth,
                        CONVERT_TO_BYTEPTR(src), src_stride, CONVERT_TO_BYTEPTR(dst),
                        dst_stride, lr_ctxt->rst_tmpbuf, optimized_lr);
            }
        }
        for (y = 0; y < tile_h; y++) {
            memcpy(src, dst, dst_stride * sizeof(*dst) << use_highbd);
            src += src_stride << use_highbd;
            dst += dst_stride << use_highbd;
        }
    }
}

void dec_av1_loop_restoration_save_boundary_lines(EbDecHandle *dec_handle,
    int after_cdef)
{
    const int num_planes = av1_num_planes(&dec_handle->seq_header.color_config);
    const int use_highbd = (dec_handle->seq_header.color_config.bit_depth > 8);

    for (int p = 0; p < num_planes; ++p) {
        LRCtxt *lr_ctxt = (LRCtxt *)dec_handle->pv_lr_ctxt;
        FrameSize *frame_size = &dec_handle->frame_header.frame_size;
        int32_t sx = 0, sy = 0;
        uint8_t *src;
        int32_t stride;
        if (p) {
            sx = dec_handle->seq_header.color_config.subsampling_x;
            sy = dec_handle->seq_header.color_config.subsampling_y;
        }

        int32_t crop_width = frame_size->frame_width >> sx;
        int32_t crop_height = frame_size->frame_height >> sy;
        EbPictureBufferDesc *cur_pic_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;
        derive_blk_pointers(cur_pic_buf, p, 0, 0, (void *)&src, &stride, sx, sy);
        uint8_t *src_buf = REAL_PTR(use_highbd, use_highbd ?
                                    CONVERT_TO_BYTEPTR(src) : src);
        int32_t src_stride = stride;
        RestorationStripeBoundaries *boundaries = &lr_ctxt->boundaries[p];

        save_tile_row_boundary_lines(src_buf, src_stride, crop_width, crop_height,
            use_highbd, p, &dec_handle->cm, after_cdef, boundaries);
    }
}
