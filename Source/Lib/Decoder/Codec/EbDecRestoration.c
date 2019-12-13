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
#include "EbDecUtils.h"
#include "EbDecProcessFrame.h"
#include "EbRestoration.h"

#define LR_PAD_SIDE 3
#define LR_PAD_MAX  (LR_PAD_SIDE << 1)

void save_tile_row_boundary_lines(uint8_t *src, int32_t src_stride,
    int32_t src_width, int32_t src_height, int32_t use_highbd, int32_t plane,
    Av1Common *cm, int32_t after_cdef, RestorationStripeBoundaries *boundaries);

void lr_generate_padding(
    EbByte  src_pic,                    //output paramter, pointer to the source picture(0,0).
    uint32_t   src_stride,              //input paramter, the stride of the source picture to be padded.
    uint32_t   original_src_width,      //input paramter, the width of the source picture which excludes the padding.
    uint32_t   original_src_height)     //input paramter, the heigth of the source picture which excludes the padding.
{
    uint32_t   verticalIdx;
    EbByte  tempSrcPic0;
    EbByte  tempSrcPic1;
    EbByte  tempSrcPic2;
    EbByte  tempSrcPic3;

    tempSrcPic0 = src_pic;
    for (verticalIdx = original_src_height; verticalIdx > 0; --verticalIdx)
    {
        // horizontal padding
        EB_MEMSET(tempSrcPic0 - LR_PAD_SIDE, *tempSrcPic0, LR_PAD_SIDE);
        EB_MEMSET(tempSrcPic0 + original_src_width, *(tempSrcPic0 + original_src_width - 1), LR_PAD_SIDE);
        tempSrcPic0 += src_stride;
    }

    // vertical padding
    tempSrcPic0 = src_pic - LR_PAD_SIDE;
    tempSrcPic1 = src_pic + (original_src_height - 1) * src_stride - LR_PAD_SIDE;
    tempSrcPic2 = tempSrcPic0;
    tempSrcPic3 = tempSrcPic1;

    for (verticalIdx = LR_PAD_SIDE; verticalIdx > 0; --verticalIdx)
    {
        // top part data copy
        tempSrcPic2 -= src_stride;
        EB_MEMCPY(tempSrcPic2, tempSrcPic0, sizeof(uint8_t) * (original_src_width + LR_PAD_MAX));
        // bottom part data copy
        tempSrcPic3 += src_stride;
        EB_MEMCPY(tempSrcPic3, tempSrcPic1, sizeof(uint8_t) * (original_src_width + LR_PAD_MAX));
    }
    return;
}

void lr_generate_padding16_bit(
    EbByte  src_pic,                       //output paramter, pointer to the source picture to be padded.
    uint32_t   src_stride,                 //input paramter, the stride of the source picture to be padded.
    uint32_t   original_src_width,         //input paramter, the width of the source picture which excludes the padding.
    uint32_t   original_src_height)        //input paramter, the height of the source picture which excludes the padding.
{
    uint32_t   verticalIdx;
    EbByte  tempSrcPic0;
    EbByte  tempSrcPic1;
    EbByte  tempSrcPic2;
    EbByte  tempSrcPic3;
    uint8_t use_highbd = 1;

    tempSrcPic0 = src_pic;
    for (verticalIdx = original_src_height; verticalIdx > 0; --verticalIdx)
    {
        // horizontal padding
        memset16bit((uint16_t*)(tempSrcPic0 - (LR_PAD_SIDE << use_highbd)), ((uint16_t*)(tempSrcPic0))[0], LR_PAD_SIDE);
        memset16bit((uint16_t*)(tempSrcPic0 + original_src_width),
            ((uint16_t*)(tempSrcPic0 + original_src_width - 2/*1*/))[0], LR_PAD_SIDE);
        tempSrcPic0 += src_stride;
    }

    // vertical padding
    tempSrcPic0 = src_pic - (LR_PAD_SIDE << use_highbd);
    tempSrcPic1 = src_pic + (original_src_height - 1) * src_stride - (LR_PAD_SIDE << use_highbd);
    tempSrcPic2 = tempSrcPic0;
    tempSrcPic3 = tempSrcPic1;
    for (verticalIdx = LR_PAD_SIDE; verticalIdx > 0; --verticalIdx)
    {
        // top part data copy
        tempSrcPic2 -= src_stride;
        EB_MEMCPY(tempSrcPic2, tempSrcPic0, sizeof(uint8_t) * (original_src_width + (LR_PAD_MAX << use_highbd)));
        // bottom part data copy
        tempSrcPic3 += src_stride;
        EB_MEMCPY(tempSrcPic3, tempSrcPic1, sizeof(uint8_t) * (original_src_width + (LR_PAD_MAX << use_highbd)));
    }

    return;
}

void lr_pad_pic(EbPictureBufferDesc *recon_picture_buf, FrameHeader *frame_hdr, EbColorConfig *color_cfg) {

    FrameSize *frame_size = &frame_hdr->frame_size;
    uint8_t sx = color_cfg->subsampling_x;
    uint8_t sy = color_cfg->subsampling_y;

    if (recon_picture_buf->bit_depth == EB_8BIT) {
        // Y samples
        lr_generate_padding(
            recon_picture_buf->buffer_y + recon_picture_buf->origin_x +
            recon_picture_buf->stride_y * recon_picture_buf->origin_y,
            recon_picture_buf->stride_y,
            frame_size->superres_upscaled_width,
            frame_size->frame_height);

        if (recon_picture_buf->color_format != EB_YUV400) {
            // Cb samples
            lr_generate_padding(
                recon_picture_buf->buffer_cb + (recon_picture_buf->origin_x >> sx) +
                recon_picture_buf->stride_cb * (recon_picture_buf->origin_y >> sy),
                recon_picture_buf->stride_cb,
                (frame_size->superres_upscaled_width + sx) >> sx,
                (frame_size->frame_height + sy) >> sy);

            // Cr samples
            lr_generate_padding(
                recon_picture_buf->buffer_cr + (recon_picture_buf->origin_x >> sx) +
                recon_picture_buf->stride_cr * (recon_picture_buf->origin_y >> sy),
                recon_picture_buf->stride_cr,
                (frame_size->superres_upscaled_width + sx) >> sx,
                (frame_size->frame_height + sy) >> sy);
        }
    }
    else {
        // Y samples
        lr_generate_padding16_bit(
            recon_picture_buf->buffer_y + (recon_picture_buf->origin_x << 1) +
            (recon_picture_buf->stride_y << 1) * recon_picture_buf->origin_y,
            recon_picture_buf->stride_y << 1,
            frame_size->superres_upscaled_width << 1,
            frame_size->frame_height);

        if (recon_picture_buf->color_format != EB_YUV400) {
            // Cb samples
            lr_generate_padding16_bit(
                recon_picture_buf->buffer_cb + ((recon_picture_buf->origin_x >> sx) << 1) +
                (recon_picture_buf->stride_cb << 1) * (recon_picture_buf->origin_y >> sy),
                recon_picture_buf->stride_cb << 1,
                ((frame_size->superres_upscaled_width + sx) >> sx) << 1,
                (frame_size->frame_height + sy) >> sy);

            // Cr samples
            lr_generate_padding16_bit(
                recon_picture_buf->buffer_cr + (recon_picture_buf->origin_x >> sx << 1) +
                (recon_picture_buf->stride_cr << 1) * (recon_picture_buf->origin_y >> sy),
                recon_picture_buf->stride_cr << 1,
                ((frame_size->superres_upscaled_width + sx) >> sx) << 1,
                (frame_size->frame_height + sy) >> sy);
        }
    }
}

void dec_av1_loop_restoration_filter_row(EbDecHandle *dec_handle, int32_t plane,
    int32_t row, int *h, int32_t sx, int32_t sy, int src_stride, int dst_stride,
    uint8_t *src, uint8_t *dst, int optimized_lr, int unit_row)
{
    AV1PixelRect tile_rect;
    RestorationTileLimits tile_limit;
    RestorationUnitInfo *lr_unit;
    int use_highbd = (dec_handle->seq_header.color_config.bit_depth > 8);
    int bit_depth = dec_handle->seq_header.color_config.bit_depth;
    LRCtxt *lr_ctxt = (LRCtxt *)dec_handle->pv_lr_ctxt;
    LRParams *lr_params = &dec_handle->frame_header.lr_params[plane];
    int ext_size = lr_params->loop_restoration_size * 3 / 2;
    int w = 0, tile_stripe0 = 0;

    tile_rect = whole_frame_rect(&dec_handle->frame_header.frame_size,
        dec_handle->seq_header.color_config.subsampling_x,
        dec_handle->seq_header.color_config.subsampling_y, plane > 0);
    int tile_h = tile_rect.bottom - tile_rect.top;
    int tile_w = tile_rect.right - tile_rect.left;

    int remaining_h = tile_h - row;
    *h = (remaining_h < ext_size) ? remaining_h :
        lr_params->loop_restoration_size;
    tile_limit.v_start = tile_rect.top + row;
    tile_limit.v_end = tile_rect.top + row + (*h);
    assert(tile_limit.v_end <= tile_rect.bottom);

    // Offset the tile upwards to align with the restoration processing stripe
    const int voffset = RESTORATION_UNIT_OFFSET >> sy;
    tile_limit.v_start = AOMMAX(tile_rect.top, tile_limit.v_start - voffset);
    if (tile_limit.v_end < tile_rect.bottom) tile_limit.v_end -= voffset;

    for (int x = 0, unit_col = 0; x < tile_w; x += w, unit_col++)
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

void dec_av1_loop_restoration_filter_frame(EbDecHandle *dec_handle,
    int optimized_lr, int enable_flag)
{
    if (!enable_flag) return;

    assert(!dec_handle->frame_header.all_lossless);

    FrameHeader *frame_header = &dec_handle->frame_header;
    LRCtxt *lr_ctxt = (LRCtxt *)dec_handle->pv_lr_ctxt;
    MasterFrameBuf *master_frame_buf = &dec_handle->master_frame_buf;
    CurFrameBuf    *frame_buf = &master_frame_buf->cur_frame_bufs[0];
    AV1PixelRect tile_rect;
    EbPictureBufferDesc *cur_pic_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;

    lr_pad_pic(cur_pic_buf, frame_header, &dec_handle->seq_header.color_config);

    lr_ctxt->lr_unit[AOM_PLANE_Y] = frame_buf->lr_unit[AOM_PLANE_Y];
    lr_ctxt->lr_unit[AOM_PLANE_U] = frame_buf->lr_unit[AOM_PLANE_U];
    lr_ctxt->lr_unit[AOM_PLANE_V] = frame_buf->lr_unit[AOM_PLANE_V];

    int num_plane = av1_num_planes(&dec_handle->seq_header.color_config);
    int use_highbd = (dec_handle->seq_header.color_config.bit_depth > 8);
    int h = 0, y, unit_row;
    int src_stride, dst_stride;
    uint8_t *src, *dst;

    for (int plane = 0; plane < num_plane; plane++)
    {
        LRParams *lr_params = &frame_header->lr_params[plane];
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
        if (lr_params->frame_restoration_type == RESTORE_NONE)
            continue;

        for (y = 0, unit_row = 0; y < tile_h; y += h, unit_row++)
        {
            dec_av1_loop_restoration_filter_row(dec_handle, plane, y, &h,
                sx, sy, src_stride, dst_stride, src, dst, optimized_lr, unit_row);
        }
        for (y = 0; y < tile_h; y++) {
            memcpy(src, dst, dst_stride * sizeof(*dst) << use_highbd);
            src += src_stride << use_highbd;
            dst += dst_stride << use_highbd;
        }
    }
}

void dec_av1_loop_restoration_save_boundary_lines(EbDecHandle *dec_handle,
        int after_cdef, int enable_flag)
{
    if (!enable_flag) return;

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
