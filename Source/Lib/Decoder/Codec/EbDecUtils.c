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
#include "EbMcp.h"
#include "EbDecBlock.h"
#include "EbDecMemInit.h"

EbErrorType check_add_tplmv_buf(EbDecHandle *dec_handle_ptr) {
    FrameHeader * ps_frm_hdr = &dec_handle_ptr->frame_header;
    const int32_t tpl_size =
        ((ps_frm_hdr->mi_rows + MAX_MIB_SIZE) >> 1) * (ps_frm_hdr->mi_stride >> 1);

    int32_t realloc = (dec_handle_ptr->master_frame_buf.tpl_mvs == NULL) ||
                      (dec_handle_ptr->master_frame_buf.tpl_mvs_size < tpl_size);

    if (realloc) {
        /*TODO: Add free now itself */
        EB_MALLOC_DEC(TemporalMvRef *,
                      dec_handle_ptr->master_frame_buf.tpl_mvs,
                      tpl_size * sizeof(*dec_handle_ptr->master_frame_buf.tpl_mvs),
                      EB_N_PTR);
        dec_handle_ptr->master_frame_buf.tpl_mvs_size = tpl_size;
    }
    return EB_ErrorNone;
}

void derive_blk_pointers(EbPictureBufferDesc *recon_picture_buf, int32_t plane, int32_t blk_col_px,
                         int32_t blk_row_px, void **pp_blk_recon_buf, int32_t *recon_stride,
                         int32_t sub_x, int32_t sub_y) {
    int32_t block_offset;

    if (plane == 0) {
        block_offset = (recon_picture_buf->origin_y + blk_row_px) * recon_picture_buf->stride_y +
                       (recon_picture_buf->origin_x + blk_col_px);
        *recon_stride = recon_picture_buf->stride_y;
    } else if (plane == 1) {
        block_offset =
            ((recon_picture_buf->origin_y >> sub_y) + blk_row_px) * recon_picture_buf->stride_cb +
            ((recon_picture_buf->origin_x >> sub_x) + blk_col_px);
        *recon_stride = recon_picture_buf->stride_cb;
    } else {
        block_offset =
            ((recon_picture_buf->origin_y >> sub_y) + blk_row_px) * recon_picture_buf->stride_cr +
            ((recon_picture_buf->origin_x >> sub_x) + blk_col_px);
        *recon_stride = recon_picture_buf->stride_cr;
    }

    if (recon_picture_buf->bit_depth != EB_8BIT || recon_picture_buf->is_16bit_pipeline) { //16bit
        if (plane == 0)
            *pp_blk_recon_buf = (void *)((uint16_t *)recon_picture_buf->buffer_y + block_offset);
        else if (plane == 1)
            *pp_blk_recon_buf = (void *)((uint16_t *)recon_picture_buf->buffer_cb + block_offset);
        else
            *pp_blk_recon_buf = (void *)((uint16_t *)recon_picture_buf->buffer_cr + block_offset);
    } else {
        if (plane == 0)
            *pp_blk_recon_buf = (void *)((uint8_t *)recon_picture_buf->buffer_y + block_offset);
        else if (plane == 1)
            *pp_blk_recon_buf = (void *)((uint8_t *)recon_picture_buf->buffer_cb + block_offset);
        else
            *pp_blk_recon_buf = (void *)((uint8_t *)recon_picture_buf->buffer_cr + block_offset);
    }
}

void pad_row(EbPictureBufferDesc *recon_picture_buf,
    EbByte buf_y, EbByte buf_cb, EbByte buf_cr,
    uint32_t row_width, uint32_t row_height,
    uint32_t pad_width, uint32_t pad_height,
    uint32_t sx, uint32_t sy, PadDir flags)
{
    uint16_t stride_y = recon_picture_buf->stride_y;
    uint16_t stride_cb = recon_picture_buf->stride_cb;
    uint16_t stride_cr = recon_picture_buf->stride_cr;
    int32_t shift = 0;
    assert(!(pad_width & 1));
    assert(!(pad_height & 1));
    if ((recon_picture_buf->bit_depth == EB_8BIT) &&
        (!recon_picture_buf->is_16bit_pipeline))
    {
        if (flags & LEFT) {
            generate_padding_l(buf_y, stride_y,
                row_height, pad_width);
            if (recon_picture_buf->color_format != EB_YUV400) {
                generate_padding_l(buf_cb, stride_cb,
                    (row_height + sy) >> sy, pad_width >> sx);
                generate_padding_l(buf_cr, stride_cr,
                    (row_height + sy) >> sy, pad_width >> sx);
            }
        }
        if (flags & RIGHT) {
            generate_padding_r(buf_y, stride_y,
                row_width, row_height, pad_width);
            if (recon_picture_buf->color_format != EB_YUV400) {
                generate_padding_r(buf_cb, stride_cb,
                    (row_width + sx) >> sx, (row_height + sy) >> sy,
                    pad_width >> sx);
                generate_padding_r(buf_cr, stride_cr,
                    (row_width + sx) >> sx, (row_height + sy) >> sy,
                    pad_width >> sx);
            }
        }
    }
    else {
        shift = 1;
        stride_y = stride_y << shift;
        stride_cb = stride_cb << shift;
        stride_cr = stride_cr << shift;
        if (flags & LEFT) {
            generate_padding_l_hbd(buf_y, stride_y,
                row_height, pad_width << shift);
            if (recon_picture_buf->color_format != EB_YUV400) {
                generate_padding_l_hbd(buf_cb, stride_cb,
                    (row_height + sy) >> sy, pad_width >> sx << shift);
                generate_padding_l_hbd(buf_cr, stride_cr,
                    (row_height + sy) >> sy, pad_width >> sx << shift);
            }
        }
        if (flags & RIGHT) {
            generate_padding_r_hbd(buf_y, stride_y,
                row_width << shift, row_height, pad_width << shift);
            if (recon_picture_buf->color_format != EB_YUV400) {
                generate_padding_r_hbd(buf_cb, stride_cb,
                    (row_width + sx) >> sx << shift,
                    (row_height + sy) >> sy,
                    pad_width >> sx << shift);
                generate_padding_r_hbd(buf_cr, stride_cr,
                    (row_width + sx) >> sx << shift,
                    (row_height + sy) >> sy,
                    pad_width >> sx << shift);
            }
        }
    }
    if (flags & TOP) {
        generate_padding_t(buf_y, stride_y,
            row_width << shift, pad_height);
        if (recon_picture_buf->color_format != EB_YUV400) {
            generate_padding_t(buf_cb, stride_cb,
                (row_width + sx) >> sx << shift, pad_height >> sy);
            generate_padding_t(buf_cr, stride_cr,
                (row_width + sx) >> sx << shift, pad_height >> sy);
        }

        if (flags & TOP_LEFT) {
            generate_padding_t(buf_y - (pad_width << shift),
                stride_y, pad_width << shift, pad_height);
            if (recon_picture_buf->color_format != EB_YUV400) {
                generate_padding_t(buf_cb - (pad_width >> sx << shift),
                    stride_cb, (pad_width + sx) >> sx << shift,
                    pad_height >> sy);
                generate_padding_t(buf_cr - (pad_width >> sx << shift),
                    stride_cr, (pad_width + sx) >> sx << shift,
                    pad_height >> sy);
            }
        }

        if (flags & TOP_RIGHT) {
            generate_padding_t(buf_y + (row_width << shift),
                stride_y, pad_width << shift, pad_height);
            if (recon_picture_buf->color_format != EB_YUV400) {
                generate_padding_t(buf_cb + (row_width >> sx << shift),
                    stride_cb, (pad_width + sx) >> sx << shift,
                    pad_height >> sy);
                generate_padding_t(buf_cr + (row_width >> sx << shift),
                    stride_cr, (pad_width + sx) >> sx << shift,
                    pad_height >> sy);
            }
        }
    }

    if (flags & BOTTOM) {
        generate_padding_b(buf_y, stride_y,
            row_width << shift, row_height, pad_height);
        if (recon_picture_buf->color_format != EB_YUV400) {
            generate_padding_b(buf_cb, stride_cb,
                (row_width + sx) >> sx << shift,
                (row_height + sy) >> sy,
                pad_height >> sy);
            generate_padding_b(buf_cr, stride_cr,
                (row_width + sx) >> sx << shift,
                (row_height + sy) >> sy,
                pad_height >> sy);
        }

        if (flags & BOTTOM_LEFT) {
            generate_padding_b(buf_y - (pad_width << shift), stride_y,
                pad_width << shift, row_height, pad_height);
            if (recon_picture_buf->color_format != EB_YUV400) {
                generate_padding_b(buf_cb - (pad_width >> sx << shift),
                    stride_cb, pad_width >> sx << shift,
                    (row_height + sy) >> sy, pad_height >> sy);
                generate_padding_b(buf_cr - (pad_width >> sx << shift),
                    stride_cr, pad_width >> sx << shift,
                    (row_height + sy) >> sy, pad_height >> sy);
            }
        }

        if (flags & BOTTOM_RIGHT) {
            generate_padding_b(buf_y + (row_width << shift), stride_y,
                pad_width << shift, row_height, pad_height);
            if (recon_picture_buf->color_format != EB_YUV400) {
                generate_padding_b(buf_cb + (row_width >> sx << shift),
                    stride_cb, pad_width >> sx << shift,
                    (row_height + sy) >> sy, pad_height >> sy);
                generate_padding_b(buf_cr + (row_width >> sx << shift),
                    stride_cr, pad_width >> sx << shift,
                    (row_height + sy) >> sy, pad_height >> sy);
            }
        }
    }
}

PadDir get_neighbour_flags(int32_t row, int32_t col,
    int32_t num_rows, int32_t num_cols)
{
    PadDir flags = 0;
    if (col == 0) flags |= LEFT;
    if (col == num_cols - 1) flags |= RIGHT;
    if (row == 0) {
        flags |= TOP;
        if (col == 0) flags |= TOP_LEFT;
        if (col == num_cols - 1) flags |= TOP_RIGHT;
    }
    if (row == num_rows - 1) {
        flags |= BOTTOM;
        if (col == 0) flags |= BOTTOM_LEFT;
        if (col == num_cols - 1) flags |= BOTTOM_RIGHT;
    }
    return flags;
}


void pad_pic(EbDecHandle *dec_handle_ptr)
{
    EbPictureBufferDesc *recon_picture_buf =
        dec_handle_ptr->cur_pic_buf[0]->ps_pic_buf;
    uint32_t frame_width = dec_handle_ptr->frame_header.
        frame_size.superres_upscaled_width;
    uint32_t frame_height = dec_handle_ptr->frame_header.
        frame_size.frame_height;

    int sx = dec_handle_ptr->seq_header.color_config.subsampling_x;
    int sy = dec_handle_ptr->seq_header.color_config.subsampling_y;

    int32_t sb_size = dec_handle_ptr->seq_header.
        use_128x128_superblock ? 128 : 64;
    int32_t sb_size_log2 = dec_handle_ptr->seq_header.sb_size_log2;
    int32_t sb_aligned_height = ALIGN_POWER_OF_TWO(frame_height, sb_size_log2);

    assert(recon_picture_buf->color_format <= EB_YUV444);
    int32_t shift = 0;
    if (recon_picture_buf->bit_depth != EB_8BIT || recon_picture_buf->is_16bit_pipeline)
        shift = 1;

    uint16_t stride_y = recon_picture_buf->stride_y << shift;
    uint16_t stride_cb = recon_picture_buf->stride_cb << shift;
    uint16_t stride_cr = recon_picture_buf->stride_cr << shift;

    uint32_t pad_width = recon_picture_buf->origin_x;
    uint32_t pad_height = recon_picture_buf->origin_y;

    /* To be adjusted according to the granularity of operation.
       Assuming a superblock row here. */
    int32_t row_id = 0;
    for (uint32_t row = 0; row < frame_height; row += sb_size, row_id++) {
        EbByte src_y = recon_picture_buf->buffer_y + (pad_width << shift) +
            (pad_height + row ) * stride_y;
        EbByte src_cb = recon_picture_buf->buffer_cb +
            (pad_width >> sx << shift) +
            ((pad_height + row) >> sy) * stride_cb;
        EbByte src_cr = recon_picture_buf->buffer_cr +
            (pad_width >> sx << shift) +
            ((pad_height + row) >> sy) * stride_cr;

        int32_t row_height = AOMMIN(sb_size, (int32_t)(frame_height - row));
        PadDir flags = get_neighbour_flags(row_id, 0,
            sb_aligned_height >> sb_size_log2, 1);

        pad_row(recon_picture_buf, src_y, src_cb, src_cr, frame_width,
            row_height, pad_width, pad_height, sx, sy, flags);
    }
}

int inverse_recenter(int r, int v) {
    if (v > 2 * r)
        return v;
    else if (v & 1)
        return r - ((v + 1) >> 1);
    else
        return r + (v >> 1);
}
