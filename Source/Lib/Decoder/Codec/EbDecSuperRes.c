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

#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "EbDecRestoration.h"
#include "EbDecUtils.h"
#include "EbDecMemInit.h"


void av1_upscale_normative_rows(const Av1Common *cm, const uint8_t *src, int src_stride,
                                uint8_t *dst, int dst_stride, int rows, int sub_x, int bd, EbBool is_16bit_pipeline);

void av1_upscale_normative_and_extend_frame(struct Av1Common *cm, FrameHeader *frm_hdr,
                                            SeqHeader *seq_hdr, EbPictureBufferDesc *src,
                                            EbPictureBufferDesc *dst) {
    const int num_planes = seq_hdr->color_config.mono_chrome ? 1 : MAX_MB_PLANE;

    for (int plane = 0; plane < num_planes; ++plane) {
        uint8_t *src_buf, *dst_buf;
        int32_t  src_stride, dst_stride;

        int sub_x = plane ? seq_hdr->color_config.subsampling_x : 0;
        int sub_y = plane ? seq_hdr->color_config.subsampling_y : 0;
        derive_blk_pointers(src, plane, 0, 0, (void *)&src_buf, &src_stride, sub_x, sub_y);
        derive_blk_pointers(dst, plane, 0, 0, (void *)&dst_buf, &dst_stride, sub_x, sub_y);

        av1_upscale_normative_rows(cm,
                                   (const uint8_t *)src_buf,
                                   src_stride,
                                   dst_buf,
                                   dst_stride,
                                   frm_hdr->frame_size.frame_height >> sub_x,
                                   sub_x,
                                   src->bit_depth,
                                   dst->is_16bit_pipeline);
    }
}

EbErrorType copy_recon(SeqHeader *seq_hdr, EbPictureBufferDesc *recon_picture_src,
                       EbPictureBufferDesc *recon_picture_dst, int num_planes) {
    recon_picture_dst->origin_x     = recon_picture_src->origin_x;
    recon_picture_dst->origin_y     = recon_picture_src->origin_y;
    recon_picture_dst->width        = recon_picture_src->width;
    recon_picture_dst->height       = recon_picture_src->height;
    recon_picture_dst->max_width    = recon_picture_src->max_width;
    recon_picture_dst->max_height   = recon_picture_src->max_height;
    recon_picture_dst->bit_depth    = recon_picture_src->bit_depth;
    recon_picture_dst->color_format = recon_picture_src->color_format;

    recon_picture_dst->stride_y  = recon_picture_src->stride_y;
    recon_picture_dst->stride_cb = recon_picture_src->stride_cb;
    recon_picture_dst->stride_cr = recon_picture_src->stride_cr;

    recon_picture_dst->luma_size   = recon_picture_src->luma_size;
    recon_picture_dst->chroma_size = recon_picture_src->chroma_size;
    recon_picture_dst->packed_flag = recon_picture_src->packed_flag;

    recon_picture_dst->stride_bit_inc_y  = recon_picture_src->stride_bit_inc_y;
    recon_picture_dst->stride_bit_inc_cb = recon_picture_src->stride_bit_inc_cb;
    recon_picture_dst->stride_bit_inc_cr = recon_picture_src->stride_bit_inc_cr;

    recon_picture_dst->buffer_enable_mask = seq_hdr->color_config.mono_chrome
                                                ? PICTURE_BUFFER_DESC_LUMA_MASK
                                                : PICTURE_BUFFER_DESC_FULL_MASK;

    recon_picture_dst->is_16bit_pipeline = recon_picture_src->is_16bit_pipeline;
    uint32_t bytes_per_pixel = (recon_picture_dst->bit_depth > EB_8BIT ||
        recon_picture_dst->is_16bit_pipeline) ? 2 : 1;

    // Allocate the Picture Buffers (luma & chroma)
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte,
                             recon_picture_dst->buffer_y,
                             recon_picture_dst->luma_size * bytes_per_pixel,
                             EB_A_PTR);
        memset(recon_picture_dst->buffer_y, 0, recon_picture_dst->luma_size * bytes_per_pixel);
    } else
        recon_picture_dst->buffer_y = 0;
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte,
                             recon_picture_dst->buffer_cb,
                             recon_picture_dst->chroma_size * bytes_per_pixel,
                             EB_A_PTR);
        memset(recon_picture_dst->buffer_cb, 0, recon_picture_dst->chroma_size * bytes_per_pixel);
    } else
        recon_picture_dst->buffer_cb = 0;
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {
        EB_ALLIGN_MALLOC_DEC(EbByte,
                             recon_picture_dst->buffer_cr,
                             recon_picture_dst->chroma_size * bytes_per_pixel,
                             EB_A_PTR);
        memset(recon_picture_dst->buffer_cr, 0, recon_picture_dst->chroma_size * bytes_per_pixel);
    } else
        recon_picture_dst->buffer_cr = 0;

    int use_highbd = (seq_hdr->color_config.bit_depth > EB_8BIT ||
        recon_picture_src->is_16bit_pipeline);

    for (int plane = 0; plane < num_planes; ++plane) {
        uint8_t *src_buf, *dst_buf;
        int32_t  src_stride, dst_stride;

        int sub_x = plane ? seq_hdr->color_config.subsampling_x : 0;
        int sub_y = plane ? seq_hdr->color_config.subsampling_y : 0;

        derive_blk_pointers(
            recon_picture_src, plane, 0, 0, (void *)&src_buf, &src_stride, sub_x, sub_y);
        derive_blk_pointers(
            recon_picture_dst, plane, 0, 0, (void *)&dst_buf, &dst_stride, sub_x, sub_y);

        int height = (recon_picture_src->height >> sub_y);
        for (int row = 0; row < height; ++row) {
            memcpy(dst_buf,
                   src_buf,
                   (recon_picture_src->width >> sub_x) * sizeof(*src_buf) << use_highbd);
            src_buf += src_stride << use_highbd;
            dst_buf += dst_stride << use_highbd;
        }
    }
    return EB_ErrorNone;
}

void av1_superres_upscale(Av1Common *cm, FrameHeader *frm_hdr, SeqHeader *seq_hdr,
                          EbPictureBufferDesc *recon_picture_src, int enable_flag) {
    if (!enable_flag) return;

    const int num_planes = seq_hdr->color_config.mono_chrome ? 1 : MAX_MB_PLANE;
    if (av1_superres_unscaled(&frm_hdr->frame_size)) return;

    EbPictureBufferDesc  recon_pic_temp;
    EbPictureBufferDesc *ps_recon_pic_temp;
    ps_recon_pic_temp = &recon_pic_temp;

    EbErrorType return_error =
        copy_recon(seq_hdr, recon_picture_src, ps_recon_pic_temp, num_planes);

    if (return_error != EB_ErrorNone) {
        ps_recon_pic_temp = NULL;
        assert(0);
    }

    uint32_t bytes_per_pixel = (recon_picture_src->bit_depth > EB_8BIT ||
        recon_picture_src->is_16bit_pipeline) ? 2 : 1;

    memset(recon_picture_src->buffer_y, 0, recon_picture_src->luma_size * bytes_per_pixel);
    memset(recon_picture_src->buffer_cb, 0, recon_picture_src->chroma_size * bytes_per_pixel);
    memset(recon_picture_src->buffer_cr, 0, recon_picture_src->chroma_size * bytes_per_pixel);

    recon_picture_src->width = frm_hdr->frame_size.superres_upscaled_width;

    av1_upscale_normative_and_extend_frame(
        cm, frm_hdr, seq_hdr, ps_recon_pic_temp, recon_picture_src);
}
