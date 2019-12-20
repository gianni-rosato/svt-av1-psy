/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

// SUMMARY
//   Contains the Decoder Utility functions

/**************************************
 * Includes
 **************************************/
#include <stdlib.h>
#include "EbDecUtils.h"
#include "EbDefinitions.h"
#include "../../Encoder/Codec/EbEntropyCoding.h"
#include "EbDecBlock.h"
#include "EbDecMemInit.h"

EbErrorType check_add_tplmv_buf(EbDecHandle *dec_handle_ptr) {

    FrameHeader *ps_frm_hdr = &dec_handle_ptr->frame_header;
    const int32_t tpl_size = ((ps_frm_hdr->mi_rows + MAX_MIB_SIZE) >> 1) *
                              (ps_frm_hdr->mi_stride >> 1);

    int32_t realloc = (dec_handle_ptr->master_frame_buf.tpl_mvs == NULL) ||
                (dec_handle_ptr->master_frame_buf.tpl_mvs_size < tpl_size);

    if (realloc) {
        /*TODO: Add free now itself */
        EB_MALLOC_DEC(TemporalMvRef *, dec_handle_ptr->master_frame_buf.tpl_mvs,
         tpl_size * sizeof(*dec_handle_ptr->master_frame_buf.tpl_mvs), EB_N_PTR);
        dec_handle_ptr->master_frame_buf.tpl_mvs_size = tpl_size;
    }
    return EB_ErrorNone;
}

void derive_blk_pointers(EbPictureBufferDesc *recon_picture_buf, int32_t plane,
                         int32_t blk_col_px, int32_t blk_row_px,
                         void **pp_blk_recon_buf, int32_t *recon_stride,
                         int32_t sub_x, int32_t sub_y)
{
    int32_t block_offset;

    if (plane == 0) {
        block_offset = (recon_picture_buf->origin_y + blk_row_px) *
            recon_picture_buf->stride_y + (recon_picture_buf->origin_x +
                blk_col_px);
        *recon_stride = recon_picture_buf->stride_y;
    }
    else if (plane == 1) {
        block_offset = ((recon_picture_buf->origin_y >> sub_y) +
            blk_row_px) * recon_picture_buf->stride_cb +
            ((recon_picture_buf->origin_x >> sub_x) + blk_col_px);
        *recon_stride = recon_picture_buf->stride_cb;
    }
    else {
        block_offset = ((recon_picture_buf->origin_y >> sub_y) +
            blk_row_px) * recon_picture_buf->stride_cr +
            ((recon_picture_buf->origin_x >> sub_x) + blk_col_px);
        *recon_stride = recon_picture_buf->stride_cr;
    }

    if (recon_picture_buf->bit_depth != EB_8BIT) {//16bit
        if (plane == 0)
            *pp_blk_recon_buf = (void *)((uint16_t*)recon_picture_buf->buffer_y
                                        + block_offset);
        else if (plane == 1)
            *pp_blk_recon_buf = (void *)((uint16_t*)recon_picture_buf->buffer_cb
                                        + block_offset);
        else
            *pp_blk_recon_buf = (void *)((uint16_t*)recon_picture_buf->buffer_cr
                                        + block_offset);
    }
    else {
        if (plane == 0)
            *pp_blk_recon_buf = (void *)((uint8_t*)recon_picture_buf->buffer_y
                                        + block_offset);
        else if (plane == 1)
            *pp_blk_recon_buf = (void *)((uint8_t*)recon_picture_buf->buffer_cb
                                        + block_offset);
        else
            *pp_blk_recon_buf = (void *)((uint8_t*)recon_picture_buf->buffer_cr
                                        + block_offset);
    }
}

void pad_pic(EbPictureBufferDesc *recon_picture_buf,
    FrameHeader *frame_hdr, int enable_flag)
{
    if (!enable_flag) return;

    FrameSize *frame_size = &frame_hdr->frame_size;
    int32_t sx = 0, sy = 0;

    switch (recon_picture_buf->color_format) {
        case EB_YUV400:
            sx = -1;
            sy = -1;
            break;
        case EB_YUV420:
            sx = 1;
            sy = 1;
            break;
        case EB_YUV422:
            sx = 1;
            sy = 0;
            break;
        case EB_YUV444:
            sx = 0;
            sy = 0;
            break;
        default:
            assert(0);
    }

    if (recon_picture_buf->bit_depth == EB_8BIT) {
        // Y samples
        generate_padding(
            recon_picture_buf->buffer_y,
            recon_picture_buf->stride_y,
            frame_size->superres_upscaled_width,
            frame_size->frame_height,
            recon_picture_buf->origin_x,
            recon_picture_buf->origin_y);

        if (recon_picture_buf->color_format != EB_YUV400) {
            // Cb samples
            generate_padding(
                recon_picture_buf->buffer_cb,
                recon_picture_buf->stride_cb,
                (frame_size->superres_upscaled_width + sx) >> sx,
                (frame_size->frame_height + sy) >> sy,
                recon_picture_buf->origin_x >> sx,
                recon_picture_buf->origin_y >> sy);

            // Cr samples
            generate_padding(
                recon_picture_buf->buffer_cr,
                recon_picture_buf->stride_cr,
                (frame_size->superres_upscaled_width + sx) >> sx,
                (frame_size->frame_height + sy) >> sy,
                recon_picture_buf->origin_x >> sx,
                recon_picture_buf->origin_y >> sy);
        }
    }
    else {
        // Y samples
        generate_padding16_bit(
            recon_picture_buf->buffer_y,
            recon_picture_buf->stride_y << 1,
            frame_size->superres_upscaled_width << 1,
            frame_size->frame_height,
            recon_picture_buf->origin_x << 1,
            recon_picture_buf->origin_y);

        if (recon_picture_buf->color_format != EB_YUV400) {
            // Cb samples
            generate_padding16_bit(
                recon_picture_buf->buffer_cb,
                recon_picture_buf->stride_cb << 1,
                ((frame_size->superres_upscaled_width + sx) >> sx) << 1,
                (frame_size->frame_height + sy) >> sy,
                recon_picture_buf->origin_x >> sx << 1,
                recon_picture_buf->origin_y >> sy);

            // Cr samples
            generate_padding16_bit(
                recon_picture_buf->buffer_cr,
                recon_picture_buf->stride_cr << 1,
                ((frame_size->superres_upscaled_width + sx) >> sx) << 1,
                (frame_size->frame_height + sy) >> sy,
                recon_picture_buf->origin_x >> sx << 1,
                recon_picture_buf->origin_y >> sy);
        }
    }
}

int inverse_recenter(int r, int v)
{
    if (v > 2 * r)
        return v;
    else if (v & 1)
        return r - ((v + 1) >> 1);
    else
        return r + (v >> 1);
}
