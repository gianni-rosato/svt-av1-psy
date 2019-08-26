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

AV1PixelRect av1_whole_frame_rect(SeqHeader *seq_header, int is_uv) {
    AV1PixelRect rect;

    int ss_x = is_uv && seq_header->color_config.subsampling_x;
    int ss_y = is_uv && seq_header->color_config.subsampling_y;

    rect.top = 0;
    rect.bottom = ROUND_POWER_OF_TWO(seq_header->max_frame_height, ss_y);
    rect.left = 0;
    rect.right = ROUND_POWER_OF_TWO(seq_header->max_frame_width, ss_x);
    return rect;
}

int av1_superres_scaled(FrameSize *frame_size) {
    return !(frame_size->frame_width == frame_size->superres_upscaled_width);
}

void dec_av1_loop_restoration_filter_frame(EbDecHandle *dec_handle, int optimized_lr)
{
    assert(!dec_handle->frame_header.all_lossless);

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
    int sb_log2 = dec_handle->seq_header.sb_size_log2;
    int h = 0, w = 0, x, y, src_stride, dst_stride, tile_stripe0 = 0;
    uint8_t *src, *dst;

    for (int plane = 0; plane < num_plane; plane++)
    {
        LRParams *lr_params = &dec_handle->frame_header.lr_params[plane];
        int ext_size = lr_params->loop_restoration_size * 3 / 2;
        int is_uv = plane > 0;
        int sx = 0, sy = 0;
        int master_col = dec_handle->master_frame_buf.sb_cols;

        if (plane) {
            sx = dec_handle->seq_header.color_config.subsampling_x;
            sy = dec_handle->seq_header.color_config.subsampling_y;
        }

        // src points to frame start
        derive_blk_pointers(cur_pic_buf, plane, 0, 0, (void *)&src,
                            &src_stride, sx, sy);

        dst = lr_ctxt->dst;
        dst_stride = lr_ctxt->dst_stride;

        tile_rect = av1_whole_frame_rect(&dec_handle->seq_header, is_uv);
        int tile_h = tile_rect.bottom - tile_rect.top;
        int tile_w = tile_rect.right - tile_rect.left;

        if (lr_params->frame_restoration_type == RESTORE_NONE)
            continue;

        for (y = 0; y < tile_h; y += h)
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

            for (x = 0; x < tile_w; x += w)
            {
                int remaining_w = tile_w - x;
                w = (remaining_w < ext_size) ? remaining_w :
                                               lr_params->loop_restoration_size;
                tile_limit.h_start = tile_rect.left + x;
                tile_limit.h_end   = tile_rect.left + x + w;

                lr_unit = lr_ctxt->lr_unit[plane] +
                    ((y >> sb_log2) * master_col) + ((x >> sb_log2));

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

static void dec_save_deblock_boundary_lines(EbDecHandle *dec_handle,
    int plane, int row, int stripe, int use_highbd, int is_above,
    RestorationStripeBoundaries *boundaries)
{
    uint8_t *src;
    int32_t stride;
    int sx = 0, sy = 0;
    int frame_width = dec_handle->frame_header.frame_size.frame_width;
    int frame_height = dec_handle->frame_header.frame_size.frame_height;

    if (plane) {
        sx = dec_handle->seq_header.color_config.subsampling_x;
        sy = dec_handle->seq_header.color_config.subsampling_y;
    }

    EbPictureBufferDesc *cur_pic_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;
    // src points to frame start
    derive_blk_pointers(cur_pic_buf, plane, 0, 0, (void *)&src, &stride, sx, sy);

    const uint8_t *src_buf = REAL_PTR(use_highbd, use_highbd ?
                                      CONVERT_TO_BYTEPTR(src) : src);

    const int src_stride = stride << use_highbd;
    const uint8_t *src_rows = src_buf + row * src_stride;

    uint8_t *bdry_buf = is_above ? boundaries->stripe_boundary_above
        : boundaries->stripe_boundary_below;
    uint8_t *bdry_start = bdry_buf + (RESTORATION_EXTRA_HORZ << use_highbd);
    const int bdry_stride = boundaries->stripe_boundary_stride << use_highbd;
    uint8_t *bdry_rows = bdry_start + RESTORATION_CTX_VERT * stripe * bdry_stride;

    // There is a rare case in which a processing stripe can end 1px above the
    // crop border. In this case, we do want to use deblocked pixels from below
    // the stripe (hence why we ended up in this function), but instead of
    // fetching 2 "below" rows we need to fetch one and duplicate it.
    // This is equivalent to clamping the sample locations against the crop border
    const int lines_to_save =
        AOMMIN(RESTORATION_CTX_VERT, (frame_height >> sx) - row);
    assert(lines_to_save == 1 || lines_to_save == 2);

    int upscaled_width = frame_width >> sx;;
    int line_bytes = 0;

    if (av1_superres_scaled(&dec_handle->frame_header.frame_size))
        assert(0);
    else {
        line_bytes = upscaled_width << use_highbd;
        for (int i = 0; i < lines_to_save; i++) {
            memcpy(bdry_rows + i * bdry_stride, src_rows + i * src_stride,
                line_bytes);
        }
    }

    // If we only saved one line, then copy it into the second line buffer
    if (lines_to_save == 1)
        memcpy(bdry_rows + bdry_stride, bdry_rows, line_bytes);

    extend_lines(bdry_rows, upscaled_width, RESTORATION_CTX_VERT, bdry_stride,
        RESTORATION_EXTRA_HORZ, use_highbd);
}

void dec_save_cdef_boundary_lines(EbDecHandle *dec_handle,
    int plane, int row, int stripe, int use_highbd, int is_above,
    RestorationStripeBoundaries *boundaries)
{
    FrameHeader *frame_info = &dec_handle->frame_header;
    int frame_width = frame_info->frame_size.frame_width;
    int sx = 0, sy = 0;
    uint8_t *src;
    int32_t stride;
    if (plane) {
        sx = dec_handle->seq_header.color_config.subsampling_x;
        sy = dec_handle->seq_header.color_config.subsampling_y;
    }

    EbPictureBufferDesc *cur_pic_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;
    // src points to frame start
    derive_blk_pointers(cur_pic_buf, plane, 0, 0, (void *)&src, &stride, sx, sy);

    const int is_uv = plane > 0;
    const uint8_t *src_buf = REAL_PTR(use_highbd, use_highbd ?
                                      CONVERT_TO_BYTEPTR(src) : src);
    const int src_stride = stride << use_highbd;
    const uint8_t *src_rows = src_buf + row * src_stride;

    uint8_t *bdry_buf = is_above ? boundaries->stripe_boundary_above
        : boundaries->stripe_boundary_below;
    uint8_t *bdry_start = bdry_buf + (RESTORATION_EXTRA_HORZ << use_highbd);
    const int bdry_stride = boundaries->stripe_boundary_stride << use_highbd;
    uint8_t *bdry_rows = bdry_start + RESTORATION_CTX_VERT * stripe * bdry_stride;
    const int src_width = frame_width >> sx;

    // At the point where this function is called, we've already applied
    // superres. So we don't need to extend the lines here, we can just
    // pull directly from the topmost row of the upscaled frame.
    const int ss_x = is_uv && sx;
    const int upscaled_width =
        av1_superres_scaled(&dec_handle->frame_header.frame_size) ?
        (frame_info->frame_size.superres_upscaled_width + ss_x) >> ss_x : src_width;

    const int line_bytes = upscaled_width << use_highbd;
    for (int i = 0; i < RESTORATION_CTX_VERT; i++) {
        // Copy the line at 'row' into both context lines. This is because
        // we want to (effectively) extend the outermost row of CDEF data
        // from this tile to produce a border, rather than using deblocked
        // pixels from the tile above/below.
        memcpy(bdry_rows + i * bdry_stride, src_rows, line_bytes);
    }
    extend_lines(bdry_rows, upscaled_width, RESTORATION_CTX_VERT, bdry_stride,
        RESTORATION_EXTRA_HORZ, use_highbd);
}

void dec_save_tile_row_boundary_lines(EbDecHandle *dec_handle, int use_highbd,
                                      int plane, int after_cdef)
{
    const int is_uv = plane > 0;
    const int ss_y = is_uv && dec_handle->seq_header.color_config.subsampling_y;
    const int stripe_height = RESTORATION_PROC_UNIT_SIZE >> ss_y;
    const int stripe_off = RESTORATION_UNIT_OFFSET >> ss_y;
    int frame_height = dec_handle->frame_header.frame_size.frame_height;
    LRCtxt *lr_ctxt = (LRCtxt *)dec_handle->pv_lr_ctxt;

    // Get the tile rectangle, with height rounded up to the next multiple of 8
    // luma pixels (only relevant for the bottom tile of the frame)
    const AV1PixelRect tile_rect =
        av1_whole_frame_rect(&dec_handle->seq_header, is_uv);
    const int stripe0 = 0;

    RestorationStripeBoundaries *boundaries = &lr_ctxt->boundaries[plane];

    const int plane_height = ROUND_POWER_OF_TWO(frame_height, ss_y);

    int tile_stripe;
    for (tile_stripe = 0;; ++tile_stripe) {
        const int rel_y0 = AOMMAX(0, tile_stripe * stripe_height - stripe_off);
        const int y0 = tile_rect.top + rel_y0;
        if (y0 >= tile_rect.bottom)
            break;

        const int rel_y1 = (tile_stripe + 1) * stripe_height - stripe_off;
        const int y1 = AOMMIN(tile_rect.top + rel_y1, tile_rect.bottom);

        const int frame_stripe = stripe0 + tile_stripe;

        // In this case, we should only use CDEF pixels at the top
        // and bottom of the frame as a whole; internal tile boundaries
        // can use deblocked pixels from adjacent tiles for context.
        const int use_deblock_above = (frame_stripe > 0);
        const int use_deblock_below = (y1 < plane_height);

        if (!after_cdef) {
            // Save deblocked context where needed.
            if (use_deblock_above) {
                dec_save_deblock_boundary_lines(dec_handle, plane,
                    y0 - RESTORATION_CTX_VERT, frame_stripe, use_highbd, 1, boundaries);
            }
            if (use_deblock_below) {
                dec_save_deblock_boundary_lines(dec_handle, plane, y1, frame_stripe,
                    use_highbd, 0, boundaries);
            }
        }
        else {
            // Save CDEF context where needed. Note that we need to save the CDEF
            // context for a particular boundary iff we *didn't* save deblocked
            // context for that boundary.
            //
            // In addition, we need to save copies of the outermost line within
            // the tile, rather than using data from outside the tile.
            if (!use_deblock_above) {
                dec_save_cdef_boundary_lines(dec_handle, plane, y0,
                    frame_stripe, use_highbd, 1, boundaries);
            }
            if (!use_deblock_below) {
                dec_save_cdef_boundary_lines(dec_handle, plane, y1 - 1,
                    frame_stripe, use_highbd, 0, boundaries);
            }
        }
    }
}

void dec_av1_loop_restoration_save_boundary_lines(EbDecHandle *dec_handle,
                                                  int after_cdef)
{
    const int num_planes = av1_num_planes(&dec_handle->seq_header.color_config);
    const int use_highbd = (dec_handle->seq_header.color_config.bit_depth > 8);
    for (int p = 0; p < num_planes; ++p)
        dec_save_tile_row_boundary_lines(dec_handle, use_highbd, p, after_cdef);
}
