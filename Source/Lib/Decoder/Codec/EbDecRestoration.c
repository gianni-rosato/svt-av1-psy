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
#include "EbDecRestoration.h"
#include "EbPictureOperators.h"
#include "EbRestoration.h"
#include "common_dsp_rtcd.h"
void save_tile_row_boundary_lines(uint8_t *src, int32_t src_stride, int32_t src_width,
                                  int32_t src_height, int32_t use_highbd, int32_t plane,
                                  Av1Common *cm, int32_t after_cdef,
                                  RestorationStripeBoundaries *boundaries);

void lr_generate_padding(
    EbByte   src_pic, //output paramter, pointer to the source picture(0,0).
    uint32_t src_stride, //input paramter, the stride of the source picture to be padded.
    uint32_t
        original_src_width, //input paramter, the width of the source picture which excludes the padding.
    uint32_t
        original_src_height) //input paramter, the heigth of the source picture which excludes the padding.
{
    uint32_t vertical_idx;
    EbByte   temp_src_pic0;
    EbByte   temp_src_pic1;
    EbByte   temp_src_pic2;
    EbByte   temp_src_pic3;

    temp_src_pic0 = src_pic;
    for (vertical_idx = original_src_height; vertical_idx > 0; --vertical_idx) {
        // horizontal padding
        EB_MEMSET(temp_src_pic0 - LR_PAD_SIDE, *temp_src_pic0, LR_PAD_SIDE);
        EB_MEMSET(temp_src_pic0 + original_src_width,
                  *(temp_src_pic0 + original_src_width - 1),
                  LR_PAD_SIDE);
        temp_src_pic0 += src_stride;
    }

    // vertical padding
    temp_src_pic0 = src_pic - LR_PAD_SIDE;
    temp_src_pic1 = src_pic + (original_src_height - 1) * src_stride - LR_PAD_SIDE;
    temp_src_pic2 = temp_src_pic0;
    temp_src_pic3 = temp_src_pic1;

    for (vertical_idx = LR_PAD_SIDE; vertical_idx > 0; --vertical_idx) {
        // top part data copy
        temp_src_pic2 -= src_stride;
        eb_memcpy(
            temp_src_pic2, temp_src_pic0, sizeof(uint8_t) * (original_src_width + LR_PAD_MAX));
        // bottom part data copy
        temp_src_pic3 += src_stride;
        eb_memcpy(
            temp_src_pic3, temp_src_pic1, sizeof(uint8_t) * (original_src_width + LR_PAD_MAX));
    }
    return;
}

void lr_generate_padding16_bit(
    EbByte   src_pic, //output paramter, pointer to the source picture to be padded.
    uint32_t src_stride, //input paramter, the stride of the source picture to be padded.
    uint32_t
        original_src_width, //input paramter, the width of the source picture which excludes the padding.
    uint32_t
        original_src_height) //input paramter, the height of the source picture which excludes the padding.
{
    uint32_t vertical_idx;
    EbByte   temp_src_pic0;
    EbByte   temp_src_pic1;
    EbByte   temp_src_pic2;
    EbByte   temp_src_pic3;
    uint8_t  use_highbd = 1;

    temp_src_pic0 = src_pic;
    for (vertical_idx = original_src_height; vertical_idx > 0; --vertical_idx) {
        // horizontal padding
        memset16bit((uint16_t *)(temp_src_pic0 - (LR_PAD_SIDE << use_highbd)),
                    ((uint16_t *)(temp_src_pic0))[0],
                    LR_PAD_SIDE);
        memset16bit((uint16_t *)(temp_src_pic0 + original_src_width),
                    ((uint16_t *)(temp_src_pic0 + original_src_width - 2 /*1*/))[0],
                    LR_PAD_SIDE);
        temp_src_pic0 += src_stride;
    }

    // vertical padding
    temp_src_pic0 = src_pic - (LR_PAD_SIDE << use_highbd);
    temp_src_pic1 = src_pic + (original_src_height - 1) * src_stride - (LR_PAD_SIDE << use_highbd);
    temp_src_pic2 = temp_src_pic0;
    temp_src_pic3 = temp_src_pic1;
    for (vertical_idx = LR_PAD_SIDE; vertical_idx > 0; --vertical_idx) {
        // top part data copy
        temp_src_pic2 -= src_stride;
        eb_memcpy(temp_src_pic2,
                  temp_src_pic0,
                  sizeof(uint8_t) * (original_src_width + (LR_PAD_MAX << use_highbd)));
        // bottom part data copy
        temp_src_pic3 += src_stride;
        eb_memcpy(temp_src_pic3,
                  temp_src_pic1,
                  sizeof(uint8_t) * (original_src_width + (LR_PAD_MAX << use_highbd)));
    }

    return;
}

void lr_pad_pic(EbPictureBufferDesc *recon_picture_buf, FrameHeader *frame_hdr,
                EbColorConfig *color_cfg)
{
    FrameSize *frame_size = &frame_hdr->frame_size;
    uint8_t    sx         = color_cfg->subsampling_x;
    uint8_t    sy         = color_cfg->subsampling_y;

    if (recon_picture_buf->bit_depth == EB_8BIT &&
        (!recon_picture_buf->is_16bit_pipeline)) {
        // Y samples
        lr_generate_padding(recon_picture_buf->buffer_y + recon_picture_buf->origin_x +
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
    } else {
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

static const StripeFilterFun stripe_filters[NUM_STRIPE_FILTERS] = { wiener_filter_stripe,
                                                                   sgrproj_filter_stripe,
                                                                   wiener_filter_stripe_highbd,
                                                                   sgrproj_filter_stripe_highbd };

// Filter one restoration unit
// Duplicated to avoid frame level buffer copy ( frame to block level copy)
// and unnecessary block copy based on LR_Type
void eb_dec_av1_loop_restoration_filter_unit(uint8_t need_bounadaries,
                                            const RestorationTileLimits *limits,
                                            const RestorationUnitInfo *rui,
                                            const RestorationStripeBoundaries *rsb,
                                            RestorationLineBuffers *rlbs,
                                            const Av1PixelRect *tile_rect, int32_t tile_stripe0,
                                            int32_t ss_x, int32_t ss_y, int32_t highbd,
                                            int32_t bit_depth, uint8_t *data8,
                                            int32_t stride, uint8_t *dst8, int32_t dst_stride,
                                            int32_t *tmpbuf, int32_t optimized_lr)
{
    RestorationType unit_rtype = rui->restoration_type;

    int32_t unit_h = limits->v_end - limits->v_start;
    int32_t unit_w = limits->h_end - limits->h_start;
    uint8_t *data8_tl = data8 + limits->v_start * stride + limits->h_start;
    uint8_t *dst8_tl = dst8;

    if (unit_rtype == RESTORE_NONE)
        return;

    const int32_t filter_idx = 2 * highbd + (unit_rtype == RESTORE_SGRPROJ);
    assert(filter_idx < NUM_STRIPE_FILTERS);
    const StripeFilterFun stripe_filter = stripe_filters[filter_idx];

    const int32_t procunit_width = RESTORATION_PROC_UNIT_SIZE >> ss_x;

    // Convolve the whole tile one stripe at a time
    RestorationTileLimits remaining_stripes = *limits;
    int32_t i = 0;
    while (i < unit_h) {
        int32_t copy_above, copy_below;
        remaining_stripes.v_start = limits->v_start + i;

        get_stripe_boundary_info(&remaining_stripes, tile_rect,
            ss_y, &copy_above, &copy_below);

        const int32_t full_stripe_height = RESTORATION_PROC_UNIT_SIZE >> ss_y;
        const int32_t runit_offset = RESTORATION_UNIT_OFFSET >> ss_y;

        // Work out where this stripe's boundaries are within
        // rsb->stripe_boundary_{above,below}
        const int32_t tile_stripe =
            (remaining_stripes.v_start - tile_rect->top + runit_offset) /
            full_stripe_height;
        const int32_t frame_stripe = tile_stripe0 + tile_stripe;
        const int32_t rsb_row = RESTORATION_CTX_VERT * frame_stripe;

        /* Calculate this stripe's height, based on two rules:
           The topmost stripe in each tile is 8 luma pixels shorter
           than usual. We can't extend past the end of the current
           restoration unit */
        const int32_t nominal_stripe_height =
            full_stripe_height - ((tile_stripe == 0) ? runit_offset : 0);
        /*In wiener filter leaf level function assumes always h to be multiple of 2.
          we can see assert related to h in this function->eb_av1_wiener_convolve_add_src_avx2*/
        const int32_t h = AOMMIN(nominal_stripe_height,
            ((remaining_stripes.v_end - remaining_stripes.v_start) + 1) & ~1);

        if (need_bounadaries)
            setup_processing_stripe_boundary(&remaining_stripes,
                rsb, rsb_row, highbd,
                h, data8, stride, rlbs, copy_above,
                copy_below, optimized_lr);

        stripe_filter(rui, unit_w, h, procunit_width, data8_tl + i * stride,
            stride, dst8_tl + i * dst_stride, dst_stride, tmpbuf, bit_depth);
        if (need_bounadaries)
            restore_processing_stripe_boundary(&remaining_stripes,
                rlbs, highbd, h, data8, stride, copy_above, copy_below,
                optimized_lr);

        i += h;
    }
    // copy block level dst buffer to src buffer
    copy_tile(unit_w, unit_h, dst8_tl, dst_stride, data8_tl, stride, highbd);
}

void dec_av1_loop_restoration_filter_row(EbDecHandle *dec_handle, int32_t sb_row,
                                         uint8_t **rec_buff, int *rec_stride, Av1PixelRect *tile_rect,
                                         int optimized_lr, uint8_t *dst, int thread_cnt)
{
    RestorationTileLimits tile_limit;
    RestorationUnitInfo *lr_unit;
    LrCtxt *lr_ctxt = (LrCtxt *)dec_handle->pv_lr_ctxt;
    const int32_t num_planes = av1_num_planes(&dec_handle->seq_header.
        color_config);
    int bit_depth = dec_handle->seq_header.color_config.bit_depth;
    int use_highbd = (bit_depth > EB_8BIT ||
        dec_handle->is_16bit_pipeline);
    int w_y = 0, tile_stripe0 = 0;
    int tile_w_y = tile_rect[AOM_PLANE_Y].right - tile_rect[AOM_PLANE_Y].left;
    int tile_h_y = tile_rect[AOM_PLANE_Y].bottom - tile_rect[AOM_PLANE_Y].top;
    uint8_t dst_stride = RESTORATION_PROC_UNIT_SIZE;

    volatile int32_t *sb_lr_completed_in_prev_row = NULL;
    int32_t *sb_lr_completed_in_row, nsync = 1;
    EbBool is_mt = dec_handle->dec_config.threads > 1;

    int32_t sb_row_idx = (is_mt == 0) ? 0 : sb_row;
    int32_t index = lr_ctxt->is_thread_min ? thread_cnt : sb_row_idx;

    if (is_mt) {
        DecMtFrameData *dec_mt_frame_data = &dec_handle->master_frame_buf.
            cur_frame_bufs[0].dec_mt_frame_data;

        if (sb_row) {
            sb_lr_completed_in_prev_row = (volatile int32_t *)
                &dec_mt_frame_data->sb_lr_completed_in_row[sb_row - 1];
        }
        sb_lr_completed_in_row =
            &dec_mt_frame_data->sb_lr_completed_in_row[sb_row];
    }

    for (int col_y = 0, sb_col_y = 0; col_y < tile_w_y; col_y += w_y, sb_col_y++) {
        int remaining_w_y = tile_w_y - col_y;
        int proc_width_y = RESTORATION_PROC_UNIT_SIZE;
        w_y = (remaining_w_y < proc_width_y) ? remaining_w_y : proc_width_y;

        /* Top-Right Sync*/
        if (is_mt) {
            if (sb_row) {
                if (col_y >= tile_w_y - w_y)
                    nsync = 0;
                while (*sb_lr_completed_in_prev_row < (sb_col_y + nsync));
            }
        }
        int sx = 0, sy = 0;
        uint8_t* src = NULL;
        uint32_t src_stride = 0;

        for(int32_t plane = 0; plane < num_planes; plane++) {

            LrParams *lr_params = &dec_handle->frame_header.lr_params[plane];
            if (lr_params->frame_restoration_type == RESTORE_NONE)
                continue;

            uint16_t lr_size = lr_params->loop_restoration_size;
            src = rec_buff[plane];
            src_stride = rec_stride[plane];

            if (plane) {
                sx = dec_handle->seq_header.color_config.subsampling_x;
                sy = dec_handle->seq_header.color_config.subsampling_y;
            }

            int plane_tile_w = (tile_w_y +sx)>> sx;
            int plane_tile_h = (tile_h_y+sy) >> sy;
            int col = (col_y+sx) >> sx;

            int sb_log2 = dec_handle->seq_header.sb_size_log2 - sy;
            int16_t sb_size = (1 << sb_log2);
            int row = sb_row << sb_log2;

            int remaining_h = plane_tile_h - row;
            int h = remaining_h < sb_size ? remaining_h : sb_size;
            tile_limit.v_start = tile_rect[plane].top + row;
            tile_limit.v_end = tile_rect[plane].top + row + h;
            assert(tile_limit.v_end <= tile_rect[plane].bottom);

            /* Offset the tile upwards to align with
               the restoration processing stripe */
            const int voffset = RESTORATION_UNIT_OFFSET >> sy;
            tile_limit.v_start = AOMMAX(tile_rect[plane].top,
                tile_limit.v_start - voffset);
            if (tile_limit.v_end < tile_rect[plane].bottom)
                tile_limit.v_end -= voffset;

            int remaining_w = plane_tile_w - col;
            //remaining_w = remaining_w >> sx;
            int proc_width = RESTORATION_PROC_UNIT_SIZE >> sx;
            int block_w = (remaining_w < proc_width) ?
                remaining_w : proc_width;
            tile_limit.h_start = tile_rect[plane].left + col;
            tile_limit.h_end = tile_rect[plane].left + col + block_w;

            uint8_t lr_size_log2 = lr_params->lr_size_log2;
            uint8_t unit_row = row >> lr_size_log2;
            uint8_t unit_col = col >> lr_size_log2;
            uint8_t lr_size_by_2 = lr_size >> 1;

            /* To use previous LR_info, remaining width/height
            should be less than 0.5 times of lr_size */
            if ((plane_tile_w - (unit_col << lr_size_log2)) < lr_size_by_2)
                unit_col = (col < lr_size_by_2) ?
                0 : (col - lr_size_by_2) >> lr_size_log2;
            if (plane_tile_h - (unit_row << lr_size_log2) < lr_size_by_2)
                unit_row = (row < lr_size_by_2) ?
                0 : (row - lr_size_by_2) >> lr_size_log2;

            // lr_info retrieval
            lr_unit = lr_ctxt->lr_unit[plane] +
                (unit_row * lr_ctxt->lr_stride[plane]) + unit_col;

            uint8_t *bdry_cdef =
                (uint8_t *)lr_ctxt->rlbs[index][plane]->tmp_save_cdef;
            uint8_t *bdry_cdef_ptr = bdry_cdef;
            uint8_t *bdry_lr =
                (uint8_t *)lr_ctxt->rlbs[index][plane]->tmp_save_lr;
            uint8_t *bdry_lr_ptr = bdry_lr;
            int width = RESTORATION_EXTRA_HORZ << use_highbd;
            int height = tile_limit.v_end - tile_limit.v_start;
            assert(height <= MAX_SB_SIZE);
            int stride = use_highbd ? src_stride << use_highbd : src_stride;
            proc_width = proc_width << use_highbd;

            uint8_t *src_proc = src + tile_limit.v_start * stride +
                (tile_limit.h_start << use_highbd) + proc_width - width;
            uint8_t *src_ptr = src_proc;

            for (int proc = 0; proc < height; proc++)
            {
                /* save cdef_data of current block to temp
                _buf before LR processing */
                memcpy(bdry_cdef_ptr, src_ptr, width);
                src_ptr += stride;
                bdry_cdef_ptr += width;
            }

            if (!use_highbd)
                eb_dec_av1_loop_restoration_filter_unit(1, &tile_limit, lr_unit,
                    &lr_ctxt->boundaries[plane], lr_ctxt->rlbs[index][plane],
                    &tile_rect[plane], tile_stripe0, sx, sy, use_highbd, bit_depth,
                    src, src_stride, dst, dst_stride,
                    lr_ctxt->rst_tmpbuf[index], optimized_lr);
            else
                eb_dec_av1_loop_restoration_filter_unit(1, &tile_limit, lr_unit,
                    &lr_ctxt->boundaries[plane], lr_ctxt->rlbs[index][plane],
                    &tile_rect[plane],
                    tile_stripe0, sx, sy, use_highbd, bit_depth,
                    CONVERT_TO_BYTEPTR(src), src_stride, CONVERT_TO_BYTEPTR(dst),
                    dst_stride, lr_ctxt->rst_tmpbuf[index], optimized_lr);
            src_ptr = src_proc - proc_width;
            // restore LR_data of previous block
            if (col)
                for (int proc = 0; proc < height; proc++) {
                    memcpy(src_ptr, bdry_lr_ptr, width);
                    src_ptr += stride;
                    bdry_lr_ptr += width;
                }

            if ((col + block_w) < plane_tile_w)
                for (int proc = 0; proc < height; proc++) {
                    // save lr_data of current block to temp_buf
                    memcpy(bdry_lr, src_proc, width);
                    // save cdef_data of current block to src
                    memcpy(src_proc, bdry_cdef, width);

                    src_proc += stride;
                    bdry_cdef += width;
                    bdry_lr += width;
                }
        }

        if (is_mt) {
            *sb_lr_completed_in_row = sb_col_y;
        }
    }
}

void dec_av1_loop_restoration_filter_frame(EbDecHandle *dec_handle, int optimized_lr,
                                           int enable_flag) {
    if (!enable_flag) return;

    assert(!dec_handle->frame_header.all_lossless);
    int32_t y,sb_row;
    uint8_t *curr_blk_recon_buf[MAX_MB_PLANE];
    int32_t curr_recon_stride[MAX_MB_PLANE];
    Av1PixelRect tile_rect[MAX_MB_PLANE];

    EbPictureBufferDesc *recon_picture_ptr =
        dec_handle->cur_pic_buf[0]->ps_pic_buf;
    const int32_t num_planes = av1_num_planes(&dec_handle->seq_header.
        color_config);

    lr_pad_pic(recon_picture_ptr, &dec_handle->frame_header,
        &dec_handle->seq_header.color_config);

    LrCtxt *lr_ctxt = (LrCtxt *)dec_handle->pv_lr_ctxt;
    uint8_t *dst    = lr_ctxt->dst;

    for (int32_t pli = 0; pli < num_planes; pli++) {
        int32_t sub_x = (pli == 0) ? 0 :
            dec_handle->seq_header.color_config.subsampling_x;
        int32_t sub_y = (pli == 0) ? 0 :
            dec_handle->seq_header.color_config.subsampling_y;

        /*Deriveing  recon pict buffer ptr's*/
        derive_blk_pointers(recon_picture_ptr,
                            pli,
                            0,
                            0,
                            (void *)&curr_blk_recon_buf[pli],
                            &curr_recon_stride[pli],
                            sub_x,
                            sub_y);

        tile_rect[pli] = whole_frame_rect(&dec_handle->frame_header.frame_size,
                                          sub_x,
                                          sub_y,
                                          pli > 0);
    }

    int tile_h = tile_rect[AOM_PLANE_Y].bottom - tile_rect[AOM_PLANE_Y].top;
    int sb_size = 1 << dec_handle->seq_header.sb_size_log2;

    for (y = 0, sb_row = 0; y < tile_h; y += sb_size, sb_row++)
    {
        dec_av1_loop_restoration_filter_row(dec_handle,
                                            sb_row,
                                            &curr_blk_recon_buf[AOM_PLANE_Y],
                                            &curr_recon_stride[AOM_PLANE_Y],
                                            tile_rect,
                                            optimized_lr,
                                            dst,
                                            0);
    }
}

void dec_av1_loop_restoration_save_boundary_lines(EbDecHandle *dec_handle,
                                                  int after_cdef) {
    const int num_planes = av1_num_planes(&dec_handle->seq_header.color_config);
    const int use_highbd = (dec_handle->seq_header.color_config.bit_depth > EB_8BIT ||
        dec_handle->is_16bit_pipeline);

    for (int p = 0; p < num_planes; ++p) {
        LrCtxt *   lr_ctxt    = (LrCtxt *)dec_handle->pv_lr_ctxt;
        FrameSize *frame_size = &dec_handle->frame_header.frame_size;
        int32_t    sx = 0, sy = 0;
        uint8_t *  src;
        int32_t    stride;
        if (p) {
            sx = dec_handle->seq_header.color_config.subsampling_x;
            sy = dec_handle->seq_header.color_config.subsampling_y;
        }

        int32_t              crop_width  = frame_size->frame_width >> sx;
        int32_t              crop_height = frame_size->frame_height >> sy;
        EbPictureBufferDesc *cur_pic_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;
        derive_blk_pointers(cur_pic_buf, p, 0, 0, (void *)&src, &stride, sx, sy);
        uint8_t *src_buf    = REAL_PTR(use_highbd, use_highbd ? CONVERT_TO_BYTEPTR(src) : src);
        int32_t  src_stride = stride;
        RestorationStripeBoundaries *boundaries = &lr_ctxt->boundaries[p];

        save_tile_row_boundary_lines(src_buf,
                                     src_stride,
                                     crop_width,
                                     crop_height,
                                     use_highbd,
                                     p,
                                     &dec_handle->cm,
                                     after_cdef,
                                     boundaries);
    }
}
