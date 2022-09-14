/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <stdlib.h>

#include "EbEncHandle.h"
#include "EbRestProcess.h"
#include "EbEncDecResults.h"
#include "EbThreads.h"
#include "EbPictureDemuxResults.h"
#include "EbReferenceObject.h"
#include "EbPictureControlSet.h"
#include "EbResourceCoordinationProcess.h"
#include "EbResize.h"

/**************************************
 * Rest Context
 **************************************/
typedef struct RestContext {
    EbDctor dctor;
    EbFifo *rest_input_fifo_ptr;
    EbFifo *rest_output_fifo_ptr;
    EbFifo *picture_demux_fifo_ptr;

    EbPictureBufferDesc *trial_frame_rst;

    EbPictureBufferDesc *
        org_rec_frame; // while doing the filtering recon gets updated uisng setup/restore processing_stripe_bounadaries
    // many threads doing the above will result in race condition.
    // each thread will hence have his own copy of recon to work on.
    // later we can have a search version that does not need the exact right recon
    int32_t *rst_tmpbuf;
} RestContext;

void pack_highbd_pic(const EbPictureBufferDesc *pic_ptr, uint16_t *buffer_16bit[3], uint32_t ss_x,
                     uint32_t ss_y, Bool include_padding);
void copy_buffer_info(EbPictureBufferDesc *src_ptr, EbPictureBufferDesc *dst_ptr);
void recon_output(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr);
void svt_av1_loop_restoration_filter_frame(Yv12BufferConfig *frame, Av1Common *cm,
                                           int32_t optimized_lr);
EbErrorType psnr_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr,
                              Bool free_memory);
EbErrorType ssim_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr,
                              Bool free_memory);
void        pad_ref_and_set_flags(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr);
void        generate_padding(EbByte src_pic, uint32_t src_stride, uint32_t original_src_width,
                             uint32_t original_src_height, uint32_t padding_width,
                             uint32_t padding_height);
void        restoration_seg_search(int32_t *rst_tmpbuf, Yv12BufferConfig *org_fts,
                                   const Yv12BufferConfig *src, Yv12BufferConfig *trial_frame_rst,
                                   PictureControlSet *pcs_ptr, uint32_t segment_index);
void        rest_finish_search(PictureControlSet *pcs_ptr);
void        svt_av1_upscale_normative_rows(const Av1Common *cm, const uint8_t *src, int src_stride,
                                           uint8_t *dst, int dst_stride, int rows, int sub_x, int bd,
                                           Bool is_16bit_pipeline);
#if DEBUG_UPSCALING
void save_YUV_to_file(char *filename, EbByte buffer_y, EbByte buffer_u, EbByte buffer_v,
                      uint16_t width, uint16_t height, uint16_t stride_y, uint16_t stride_u,
                      uint16_t stride_v, uint16_t origin_y, uint16_t origin_x, uint32_t ss_x,
                      uint32_t ss_y);
#endif

static void rest_context_dctor(EbPtr p) {
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    RestContext     *obj                = (RestContext *)thread_context_ptr->priv;
    EB_DELETE(obj->trial_frame_rst);
    // buffer only malloc'd if boundaries are used in rest. search.
    // see scs_ptr->seq_header.use_boundaries_in_rest_search
    if (obj->org_rec_frame)
        EB_DELETE(obj->org_rec_frame);
    EB_FREE_ALIGNED(obj->rst_tmpbuf);
    EB_FREE_ARRAY(obj);
}

/******************************************************
 * Rest Context Constructor
 ******************************************************/
EbErrorType rest_context_ctor(EbThreadContext   *thread_context_ptr,
                              const EbEncHandle *enc_handle_ptr, EbPtr object_init_data_ptr,
                              int index, int demux_index) {
    const SequenceControlSet       *scs_ptr       = enc_handle_ptr->scs_instance_array[0]->scs_ptr;
    const EbSvtAv1EncConfiguration *config        = &scs_ptr->static_config;
    EbColorFormat                   color_format  = config->encoder_color_format;
    EbPictureBufferDescInitData    *init_data_ptr = (EbPictureBufferDescInitData *)
        object_init_data_ptr;
    RestContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = rest_context_dctor;

    // Input/Output System Resource Manager FIFOs
    context_ptr->rest_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->cdef_results_resource_ptr, index);
    context_ptr->rest_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->rest_results_resource_ptr, index);
    context_ptr->picture_demux_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);

    Bool is_16bit = scs_ptr->is_16bit_pipeline;

    if (get_enable_restoration(init_data_ptr->enc_mode,
                               config->enable_restoration_filtering,
                               scs_ptr->input_resolution,
                               scs_ptr->static_config.fast_decode)) {
        EbPictureBufferDescInitData init_data;

        init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        init_data.max_width          = (uint16_t)scs_ptr->max_input_luma_width;
        init_data.max_height         = (uint16_t)scs_ptr->max_input_luma_height;
        init_data.bit_depth          = is_16bit ? EB_SIXTEEN_BIT : EB_EIGHT_BIT;
        init_data.color_format       = color_format;
        init_data.left_padding       = AOM_RESTORATION_FRAME_BORDER;
        init_data.right_padding      = AOM_RESTORATION_FRAME_BORDER;
        init_data.top_padding        = AOM_RESTORATION_FRAME_BORDER;
        init_data.bot_padding        = AOM_RESTORATION_FRAME_BORDER;
        init_data.split_mode         = FALSE;
        init_data.is_16bit_pipeline  = is_16bit;

        EB_NEW(context_ptr->trial_frame_rst, svt_picture_buffer_desc_ctor, (EbPtr)&init_data);
        if (scs_ptr->use_boundaries_in_rest_search)
            EB_NEW(context_ptr->org_rec_frame, svt_picture_buffer_desc_ctor, (EbPtr)&init_data);
        else
            context_ptr->org_rec_frame = NULL;
        if (!is_16bit) {
            context_ptr->trial_frame_rst->bit_depth = EB_EIGHT_BIT;
            if (scs_ptr->use_boundaries_in_rest_search)
                context_ptr->org_rec_frame->bit_depth = EB_EIGHT_BIT;
        }

        EB_MALLOC_ALIGNED(context_ptr->rst_tmpbuf, RESTORATION_TMPBUF_SIZE);
    }

    return EB_ErrorNone;
}
extern void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr,
                          Bool is_highbd) {
    if (!is_highbd) {
        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == TRUE)
            *recon_ptr = ((EbReferenceObject *)
                              pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                             ->reference_picture;
        else
            *recon_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr; //OMK
    } else {
        *recon_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr;
    }

    // recon buffer is created in full resolution, it is resized to difference size
    // when reference scaling enabled. recon width and height should be adjusted to
    // upscaled render size
    if (*recon_ptr &&
        (pcs_ptr->parent_pcs_ptr->render_width != (*recon_ptr)->width ||
         pcs_ptr->parent_pcs_ptr->render_height != (*recon_ptr)->height)) {
        (*recon_ptr)->width  = pcs_ptr->parent_pcs_ptr->render_width;
        (*recon_ptr)->height = pcs_ptr->parent_pcs_ptr->render_height;
    }
}

// If using boundaries during the filter search, copy the recon pic to a new buffer (to
// avoid race conidition from many threads modifying the same recon pic).
//
// If not using boundaries during the filter search, return the input recon picture location
// to be used in restoration search (save cycles/memory of copying pic to a new buffer).
// The recon pic should not be modified during the search, otherwise there will be a race
// condition between threads.
//
// Return a pointer to the recon pic to be used during the restoration search.
EbPictureBufferDesc *get_own_recon(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                                   RestContext *context_ptr, Bool is_16bit) {
    const uint32_t ss_x = scs_ptr->subsampling_x;
    const uint32_t ss_y = scs_ptr->subsampling_y;

    EbPictureBufferDesc *recon_picture_ptr;
    get_recon_pic(pcs_ptr, &recon_picture_ptr, is_16bit);
    // if boundaries are not used, don't need to copy pic to new buffer, as the
    // search will not modify the pic
    if (!scs_ptr->use_boundaries_in_rest_search) {
        return recon_picture_ptr;
    }
    if (is_16bit) {
        uint16_t *rec_ptr = (uint16_t *)recon_picture_ptr->buffer_y + recon_picture_ptr->origin_x +
            recon_picture_ptr->origin_y * recon_picture_ptr->stride_y;
        uint16_t *rec_ptr_cb = (uint16_t *)recon_picture_ptr->buffer_cb +
            recon_picture_ptr->origin_x / 2 +
            recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb;
        uint16_t *rec_ptr_cr = (uint16_t *)recon_picture_ptr->buffer_cr +
            recon_picture_ptr->origin_x / 2 +
            recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr;

        EbPictureBufferDesc *org_rec = context_ptr->org_rec_frame;
        uint16_t            *org_ptr = (uint16_t *)org_rec->buffer_y + org_rec->origin_x +
            org_rec->origin_y * org_rec->stride_y;
        uint16_t *org_ptr_cb = (uint16_t *)org_rec->buffer_cb + org_rec->origin_x / 2 +
            org_rec->origin_y / 2 * org_rec->stride_cb;
        uint16_t *org_ptr_cr = (uint16_t *)org_rec->buffer_cr + org_rec->origin_x / 2 +
            org_rec->origin_y / 2 * org_rec->stride_cr;

        for (int r = 0; r < recon_picture_ptr->height; ++r)
            svt_memcpy(org_ptr + r * org_rec->stride_y,
                       rec_ptr + r * recon_picture_ptr->stride_y,
                       recon_picture_ptr->width << 1);

        for (int r = 0; r < (recon_picture_ptr->height >> ss_y); ++r) {
            svt_memcpy(org_ptr_cb + r * org_rec->stride_cb,
                       rec_ptr_cb + r * recon_picture_ptr->stride_cb,
                       (recon_picture_ptr->width >> ss_x) << 1);
            svt_memcpy(org_ptr_cr + r * org_rec->stride_cr,
                       rec_ptr_cr + r * recon_picture_ptr->stride_cr,
                       (recon_picture_ptr->width >> ss_x) << 1);
        }
    } else {
        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == TRUE)
            recon_picture_ptr = ((EbReferenceObject *)pcs_ptr->parent_pcs_ptr
                                     ->reference_picture_wrapper_ptr->object_ptr)
                                    ->reference_picture;
        else
            recon_picture_ptr = pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr; //OMK
        // if boundaries are not used, don't need to copy pic to new buffer, as the
        // search will not modify the pic
        if (!scs_ptr->use_boundaries_in_rest_search) {
            return recon_picture_ptr;
        }
        uint8_t *rec_ptr    = &((recon_picture_ptr->buffer_y)[recon_picture_ptr->origin_x +
                                                           recon_picture_ptr->origin_y *
                                                               recon_picture_ptr->stride_y]);
        uint8_t *rec_ptr_cb = &((recon_picture_ptr->buffer_cb)[recon_picture_ptr->origin_x / 2 +
                                                               recon_picture_ptr->origin_y / 2 *
                                                                   recon_picture_ptr->stride_cb]);
        uint8_t *rec_ptr_cr = &((recon_picture_ptr->buffer_cr)[recon_picture_ptr->origin_x / 2 +
                                                               recon_picture_ptr->origin_y / 2 *
                                                                   recon_picture_ptr->stride_cr]);

        EbPictureBufferDesc *org_rec = context_ptr->org_rec_frame;
        uint8_t             *org_ptr = &(
            (org_rec->buffer_y)[org_rec->origin_x + org_rec->origin_y * org_rec->stride_y]);
        uint8_t *org_ptr_cb = &((org_rec->buffer_cb)[org_rec->origin_x / 2 +
                                                     org_rec->origin_y / 2 * org_rec->stride_cb]);
        uint8_t *org_ptr_cr = &((org_rec->buffer_cr)[org_rec->origin_x / 2 +
                                                     org_rec->origin_y / 2 * org_rec->stride_cr]);

        for (int r = 0; r < recon_picture_ptr->height; ++r)
            svt_memcpy(org_ptr + r * org_rec->stride_y,
                       rec_ptr + r * recon_picture_ptr->stride_y,
                       recon_picture_ptr->width);

        for (int r = 0; r < (recon_picture_ptr->height >> ss_y); ++r) {
            svt_memcpy(org_ptr_cb + r * org_rec->stride_cb,
                       rec_ptr_cb + r * recon_picture_ptr->stride_cb,
                       (recon_picture_ptr->width >> ss_x));
            svt_memcpy(org_ptr_cr + r * org_rec->stride_cr,
                       rec_ptr_cr + r * recon_picture_ptr->stride_cr,
                       (recon_picture_ptr->width >> ss_x));
        }
    }
    return context_ptr->org_rec_frame;
}

void svt_convert_pic_8bit_to_16bit(EbPictureBufferDesc *src_8bit, EbPictureBufferDesc *dst_16bit,
                                   uint16_t ss_x, uint16_t ss_y) {
    //copy input from 8bit to 16bit
    uint8_t  *buffer_8bit;
    int32_t   stride_8bit;
    uint16_t *buffer_16bit;
    int32_t   stride_16bit;
    // Y
    buffer_16bit = (uint16_t *)(dst_16bit->buffer_y) + dst_16bit->origin_x +
        dst_16bit->origin_y * dst_16bit->stride_y;
    stride_16bit = dst_16bit->stride_y;
    buffer_8bit = src_8bit->buffer_y + src_8bit->origin_x + src_8bit->origin_y * src_8bit->stride_y;
    stride_8bit = src_8bit->stride_y;

    svt_convert_8bit_to_16bit(
        buffer_8bit, stride_8bit, buffer_16bit, stride_16bit, src_8bit->width, src_8bit->height);

    // Cb
    buffer_16bit = (uint16_t *)(dst_16bit->buffer_cb) + (dst_16bit->origin_x >> ss_x) +
        (dst_16bit->origin_y >> ss_y) * dst_16bit->stride_cb;
    stride_16bit = dst_16bit->stride_cb;
    buffer_8bit  = src_8bit->buffer_cb + (src_8bit->origin_x >> ss_x) +
        (src_8bit->origin_y >> ss_y) * src_8bit->stride_cb;
    stride_8bit = src_8bit->stride_cb;

    svt_convert_8bit_to_16bit(buffer_8bit,
                              stride_8bit,
                              buffer_16bit,
                              stride_16bit,
                              src_8bit->width >> ss_x,
                              src_8bit->height >> ss_y);

    // Cr
    buffer_16bit = (uint16_t *)(dst_16bit->buffer_cr) + (dst_16bit->origin_x >> ss_x) +
        (dst_16bit->origin_y >> ss_y) * dst_16bit->stride_cr;
    stride_16bit = dst_16bit->stride_cr;
    buffer_8bit  = src_8bit->buffer_cr + (src_8bit->origin_x >> ss_x) +
        (src_8bit->origin_y >> ss_y) * src_8bit->stride_cr;
    stride_8bit = src_8bit->stride_cr;

    svt_convert_8bit_to_16bit(buffer_8bit,
                              stride_8bit,
                              buffer_16bit,
                              stride_16bit,
                              src_8bit->width >> ss_x,
                              src_8bit->height >> ss_y);

    dst_16bit->width  = src_8bit->width;
    dst_16bit->height = src_8bit->height;
}

extern void pack_2d_pic(EbPictureBufferDesc *input_picture, uint16_t *packed[3]);

void set_unscaled_input_16bit(PictureControlSet *pcs_ptr) {
    EbPictureBufferDesc *input_pic  = pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;
    EbPictureBufferDesc *output_pic = pcs_ptr->input_frame16bit;
    uint16_t             ss_x       = pcs_ptr->parent_pcs_ptr->scs_ptr->subsampling_x;
    uint16_t             ss_y       = pcs_ptr->parent_pcs_ptr->scs_ptr->subsampling_y;
    copy_buffer_info(input_pic, pcs_ptr->input_frame16bit);
    if (input_pic->bit_depth == EB_EIGHT_BIT)
        svt_convert_pic_8bit_to_16bit(input_pic, output_pic, ss_x, ss_y);
    else {
        uint16_t *planes[3] = {(uint16_t *)output_pic->buffer_y +
                                   (output_pic->origin_y * output_pic->stride_y) +
                                   (output_pic->origin_x),
                               (uint16_t *)output_pic->buffer_cb +
                                   (((output_pic->origin_y) >> ss_y) * output_pic->stride_cb) +
                                   ((output_pic->origin_x) >> ss_x),
                               (uint16_t *)output_pic->buffer_cr +
                                   (((output_pic->origin_y) >> ss_y) * output_pic->stride_cr) +
                                   ((output_pic->origin_x) >> ss_x)};
        pack_2d_pic(input_pic, planes);
    }
}

void derive_blk_pointers_enc(EbPictureBufferDesc *recon_picture_buf, int32_t plane,
                             int32_t blk_col_px, int32_t blk_row_px, void **pp_blk_recon_buf,
                             int32_t *recon_stride, int32_t sub_x, int32_t sub_y, Bool use_highbd) {
    int32_t block_offset;

    if (plane == 0) {
        block_offset = (recon_picture_buf->origin_y + blk_row_px) * recon_picture_buf->stride_y +
            (recon_picture_buf->origin_x + blk_col_px);
        *recon_stride = recon_picture_buf->stride_y;
    } else if (plane == 1) {
        block_offset = ((recon_picture_buf->origin_y >> sub_y) + blk_row_px) *
                recon_picture_buf->stride_cb +
            ((recon_picture_buf->origin_x >> sub_x) + blk_col_px);
        *recon_stride = recon_picture_buf->stride_cb;
    } else {
        block_offset = ((recon_picture_buf->origin_y >> sub_y) + blk_row_px) *
                recon_picture_buf->stride_cr +
            ((recon_picture_buf->origin_x >> sub_x) + blk_col_px);
        *recon_stride = recon_picture_buf->stride_cr;
    }

    if (use_highbd) { //16bit
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

EbErrorType copy_recon_enc(SequenceControlSet *scs_ptr, EbPictureBufferDesc *recon_picture_src,
                           EbPictureBufferDesc *recon_picture_dst, int num_planes, int skip_copy) {
    recon_picture_dst->origin_x     = recon_picture_src->origin_x;
    recon_picture_dst->origin_y     = recon_picture_src->origin_y;
    recon_picture_dst->origin_bot_y = recon_picture_src->origin_bot_y;
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

    recon_picture_dst->buffer_enable_mask = scs_ptr->seq_header.color_config.mono_chrome
        ? PICTURE_BUFFER_DESC_LUMA_MASK
        : PICTURE_BUFFER_DESC_FULL_MASK;

    uint32_t bytesPerPixel = scs_ptr->is_16bit_pipeline ? 2 : 1;

    // Allocate the Picture Buffers (luma & chroma)
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_y,
                          recon_picture_dst->luma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_y, 0, recon_picture_dst->luma_size * bytesPerPixel);
    } else
        recon_picture_dst->buffer_y = 0;
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_cb,
                          recon_picture_dst->chroma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_cb, 0, recon_picture_dst->chroma_size * bytesPerPixel);
    } else
        recon_picture_dst->buffer_cb = 0;
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_cr,
                          recon_picture_dst->chroma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_cr, 0, recon_picture_dst->chroma_size * bytesPerPixel);
    } else
        recon_picture_dst->buffer_cr = 0;

    int use_highbd = scs_ptr->is_16bit_pipeline;

    if (!skip_copy) {
        for (int plane = 0; plane < num_planes; ++plane) {
            uint8_t *src_buf, *dst_buf;
            int32_t  src_stride, dst_stride;

            int sub_x = plane ? scs_ptr->subsampling_x : 0;
            int sub_y = plane ? scs_ptr->subsampling_y : 0;

            derive_blk_pointers_enc(recon_picture_src,
                                    plane,
                                    0,
                                    0,
                                    (void *)&src_buf,
                                    &src_stride,
                                    sub_x,
                                    sub_y,
                                    use_highbd);
            derive_blk_pointers_enc(recon_picture_dst,
                                    plane,
                                    0,
                                    0,
                                    (void *)&dst_buf,
                                    &dst_stride,
                                    sub_x,
                                    sub_y,
                                    use_highbd);

            int height = ((recon_picture_src->height + sub_y) >> sub_y);
            for (int row = 0; row < height; ++row) {
                svt_memcpy(dst_buf,
                           src_buf,
                           ((recon_picture_src->width + sub_x) >> sub_x) * sizeof(*src_buf)
                               << use_highbd);
                src_buf += src_stride << use_highbd;
                dst_buf += dst_stride << use_highbd;
            }
        }
    }

    return EB_ErrorNone;
}

void svt_av1_superres_upscale_frame(struct Av1Common *cm, PictureControlSet *pcs_ptr,
                                    SequenceControlSet *scs_ptr) {
    // Set these parameters for testing since they are not correctly populated yet
    EbPictureBufferDesc *recon_ptr;

    Bool is_16bit = scs_ptr->is_16bit_pipeline;

    get_recon_pic(pcs_ptr, &recon_ptr, is_16bit);

    uint16_t  ss_x       = scs_ptr->subsampling_x;
    uint16_t  ss_y       = scs_ptr->subsampling_y;
    const int num_planes = scs_ptr->seq_header.color_config.mono_chrome ? 1 : MAX_MB_PLANE;

    EbPictureBufferDesc  recon_pic_temp;
    EbPictureBufferDesc *ps_recon_pic_temp;
    ps_recon_pic_temp = &recon_pic_temp;

    EbErrorType return_error = copy_recon_enc(scs_ptr, recon_ptr, ps_recon_pic_temp, num_planes, 0);

    if (return_error != EB_ErrorNone) {
        ps_recon_pic_temp = NULL;
        assert(0);
    }

    EbPictureBufferDesc *src = ps_recon_pic_temp;
    EbPictureBufferDesc *dst = recon_ptr;

    // get the bit-depth from the encoder config instead of from the recon ptr
    int bit_depth = scs_ptr->static_config.encoder_bit_depth;

    for (int plane = 0; plane < num_planes; ++plane) {
        uint8_t *src_buf, *dst_buf;
        int32_t  src_stride, dst_stride;

        int sub_x = plane ? ss_x : 0;
        int sub_y = plane ? ss_y : 0;
        derive_blk_pointers_enc(
            src, plane, 0, 0, (void *)&src_buf, &src_stride, sub_x, sub_y, is_16bit);
        derive_blk_pointers_enc(
            dst, plane, 0, 0, (void *)&dst_buf, &dst_stride, sub_x, sub_y, is_16bit);

        svt_av1_upscale_normative_rows(cm,
                                       (const uint8_t *)src_buf,
                                       src_stride,
                                       dst_buf,
                                       dst_stride,
                                       (src->height + sub_y) >> sub_y,
                                       sub_x,
                                       bit_depth,
                                       is_16bit);
    }

    // free the memory
    EB_FREE_ALIGNED_ARRAY(ps_recon_pic_temp->buffer_y);
    EB_FREE_ALIGNED_ARRAY(ps_recon_pic_temp->buffer_cb);
    EB_FREE_ALIGNED_ARRAY(ps_recon_pic_temp->buffer_cr);
}

static void copy_statistics_to_ref_obj_ect(PictureControlSet *pcs, SequenceControlSet *scs) {
    EbReferenceObject *obj = (EbReferenceObject *)
                                 pcs->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;

    pcs->intra_coded_area = (100 * pcs->intra_coded_area) /
        (pcs->parent_pcs_ptr->aligned_width * pcs->parent_pcs_ptr->aligned_height);
    pcs->skip_coded_area = (100 * pcs->skip_coded_area) /
        (pcs->parent_pcs_ptr->aligned_width * pcs->parent_pcs_ptr->aligned_height);

    if (pcs->slice_type == I_SLICE)
        pcs->intra_coded_area = 0;
    obj->intra_coded_area                   = (uint8_t)(pcs->intra_coded_area);
    obj->skip_coded_area                    = (uint8_t)(pcs->skip_coded_area);
    struct PictureParentControlSet *ppcs    = pcs->parent_pcs_ptr;
    FrameHeader                    *frm_hdr = &ppcs->frm_hdr;

    struct LoopFilter *const lf = &frm_hdr->loop_filter_params;

    obj->filter_level[0] = lf->filter_level[0];
    obj->filter_level[1] = lf->filter_level[1];
    obj->filter_level_u  = lf->filter_level_u;
    obj->filter_level_v  = lf->filter_level_v;

    obj->ref_cdef_strengths_num = ppcs->nb_cdef_strengths;
    for (int i = 0; i < ppcs->nb_cdef_strengths; i++) {
        obj->ref_cdef_strengths[0][i] = frm_hdr->cdef_params.cdef_y_strength[i];
        obj->ref_cdef_strengths[1][i] = frm_hdr->cdef_params.cdef_uv_strength[i];
    }
    uint32_t sb_index;
    for (sb_index = 0; sb_index < pcs->b64_total_count; ++sb_index) {
        obj->sb_intra[sb_index]           = pcs->sb_intra[sb_index];
        obj->sb_skip[sb_index]            = pcs->sb_skip[sb_index];
        obj->sb_64x64_mvp[sb_index]       = pcs->sb_64x64_mvp[sb_index];
        obj->sb_me_64x64_dist[sb_index]   = pcs->parent_pcs_ptr->me_64x64_distortion[sb_index];
        obj->sb_me_8x8_cost_var[sb_index] = pcs->parent_pcs_ptr->me_8x8_cost_variance[sb_index];
    }
    obj->tmp_layer_idx   = (uint8_t)pcs->temporal_layer_index;
    obj->is_scene_change = pcs->parent_pcs_ptr->scene_change_flag;

    Av1Common *cm    = pcs->parent_pcs_ptr->av1_cm;
    obj->sg_frame_ep = cm->sg_frame_ep;
    if (scs->mfmv_enabled || !pcs->parent_pcs_ptr->is_not_scaled) {
        obj->frame_type = pcs->parent_pcs_ptr->frm_hdr.frame_type;
        obj->order_hint = pcs->parent_pcs_ptr->cur_order_hint;
        svt_memcpy(obj->ref_order_hint, pcs->parent_pcs_ptr->ref_order_hint, 7 * sizeof(uint32_t));
    }
    // Copy the prev frame wn filter coeffs
    if (cm->wn_filter_ctrls.enabled && cm->wn_filter_ctrls.use_prev_frame_coeffs) {
        for (int32_t plane = 0; plane < MAX_MB_PLANE; ++plane) {
            int32_t ntiles = pcs->rst_info[plane].units_per_tile;
            for (int32_t u = 0; u < ntiles; ++u) {
                obj->unit_info[plane][u].restoration_type =
                    pcs->rst_info[plane].unit_info[u].restoration_type;
                if (pcs->rst_info[plane].unit_info[u].restoration_type == RESTORE_WIENER)
                    obj->unit_info[plane][u].wiener_info =
                        pcs->rst_info[plane].unit_info[u].wiener_info;
            }
        }
    }
}

/******************************************************
 * Rest Kernel
 ******************************************************/
void *rest_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext    *thread_context_ptr = (EbThreadContext *)input_ptr;
    RestContext        *context_ptr        = (RestContext *)thread_context_ptr->priv;
    PictureControlSet  *pcs_ptr;
    SequenceControlSet *scs_ptr;

    //// Input
    EbObjectWrapper *cdef_results_wrapper_ptr;
    CdefResults     *cdef_results_ptr;

    //// Output
    EbObjectWrapper     *rest_results_wrapper_ptr;
    RestResults         *rest_results_ptr;
    EbObjectWrapper     *picture_demux_results_wrapper_ptr;
    PictureDemuxResults *picture_demux_results_rtr;
    // SB Loop variables
    uint8_t tile_cols;
    uint8_t tile_rows;

    Bool superres_recode = FALSE;

    for (;;) {
        // Get Cdef Results
        EB_GET_FULL_OBJECT(context_ptr->rest_input_fifo_ptr, &cdef_results_wrapper_ptr);

        cdef_results_ptr = (CdefResults *)cdef_results_wrapper_ptr->object_ptr;
        pcs_ptr          = (PictureControlSet *)cdef_results_ptr->pcs_wrapper_ptr->object_ptr;
        PictureParentControlSet *ppcs = pcs_ptr->parent_pcs_ptr;
        scs_ptr                       = pcs_ptr->scs_ptr;
        FrameHeader *frm_hdr          = &pcs_ptr->parent_pcs_ptr->frm_hdr;
        Bool         is_16bit         = scs_ptr->is_16bit_pipeline;
        Av1Common   *cm               = pcs_ptr->parent_pcs_ptr->av1_cm;
        if (ppcs->enable_restoration && frm_hdr->allow_intrabc == 0) {
            // If using boundaries during the filter search, copy the recon pic to a new buffer (to
            // avoid race condition from many threads modifying the same recon pic).
            //
            // If not using boundaries during the filter search, copy the input recon picture location
            // to be used in restoration search (save cycles/memory of copying pic to a new buffer).
            // The recon pic should not be modified during the search, otherwise there will be a race
            // condition between threads.
            EbPictureBufferDesc *recon_picture_ptr = get_own_recon(
                scs_ptr, pcs_ptr, context_ptr, is_16bit);
            EbPictureBufferDesc *input_picture_ptr = is_16bit
                ? pcs_ptr->input_frame16bit
                : pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr;

            EbPictureBufferDesc *scaled_input_picture_ptr = NULL;
            // downscale input picture if recon is resized
            Bool is_resized = recon_picture_ptr->width != input_picture_ptr->width ||
                recon_picture_ptr->height != input_picture_ptr->height;
            if (is_resized) {
                scaled_input_picture_ptr = pcs_ptr->scaled_input_picture_ptr;
                input_picture_ptr        = scaled_input_picture_ptr;
            }

            // there are padding pixels if input pics are not 8 pixel aligned
            // but there is no extra padding after input pics are resized for
            // reference scaling
            Yv12BufferConfig cpi_source;
            link_eb_to_aom_buffer_desc(input_picture_ptr,
                                       &cpi_source,
                                       is_resized ? 0 : scs_ptr->max_input_pad_right,
                                       is_resized ? 0 : scs_ptr->max_input_pad_bottom,
                                       is_16bit);

            Yv12BufferConfig trial_frame_rst;
            link_eb_to_aom_buffer_desc(context_ptr->trial_frame_rst,
                                       &trial_frame_rst,
                                       is_resized ? 0 : scs_ptr->max_input_pad_right,
                                       is_resized ? 0 : scs_ptr->max_input_pad_bottom,
                                       is_16bit);

            Yv12BufferConfig org_fts;
            link_eb_to_aom_buffer_desc(recon_picture_ptr,
                                       &org_fts,
                                       is_resized ? 0 : scs_ptr->max_input_pad_right,
                                       is_resized ? 0 : scs_ptr->max_input_pad_bottom,
                                       is_16bit);

            if (pcs_ptr->parent_pcs_ptr->slice_type != I_SLICE && cm->wn_filter_ctrls.enabled &&
                cm->wn_filter_ctrls.use_prev_frame_coeffs) {
                EbReferenceObject *ref_obj_l0 =
                    (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                for (int32_t plane = 0; plane < MAX_MB_PLANE; ++plane) {
                    int32_t ntiles = pcs_ptr->rst_info[plane].units_per_tile;
                    for (int32_t u = 0; u < ntiles; ++u) {
                        pcs_ptr->rst_info[plane].unit_info[u].restoration_type =
                            ref_obj_l0->unit_info[plane][u].restoration_type;
                        if (ref_obj_l0->unit_info[plane][u].restoration_type == RESTORE_WIENER)
                            pcs_ptr->rst_info[plane].unit_info[u].wiener_info =
                                ref_obj_l0->unit_info[plane][u].wiener_info;
                    }
                }
            }
            restoration_seg_search(context_ptr->rst_tmpbuf,
                                   &org_fts,
                                   &cpi_source,
                                   &trial_frame_rst,
                                   pcs_ptr,
                                   cdef_results_ptr->segment_index);
        }

        //all seg based search is done. update total processed segments. if all done, finish the search and perfrom application.
        svt_block_on_mutex(pcs_ptr->rest_search_mutex);

        pcs_ptr->tot_seg_searched_rest++;
        if (pcs_ptr->tot_seg_searched_rest == pcs_ptr->rest_segments_total_count) {
            if (ppcs->enable_restoration && frm_hdr->allow_intrabc == 0) {
                rest_finish_search(pcs_ptr);

                // Only need recon if REF pic or recon is output
                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag ||
                    scs_ptr->static_config.recon_enabled) {
                    if (pcs_ptr->rst_info[0].frame_restoration_type != RESTORE_NONE ||
                        pcs_ptr->rst_info[1].frame_restoration_type != RESTORE_NONE ||
                        pcs_ptr->rst_info[2].frame_restoration_type != RESTORE_NONE) {
                        svt_av1_loop_restoration_filter_frame(cm->frame_to_show, cm, 0);
                    }
                }

                if (cm->sg_filter_ctrls.enabled) {
                    uint8_t best_ep_cnt = 0;
                    uint8_t best_ep     = 0;
                    for (uint8_t i = 0; i < SGRPROJ_PARAMS; i++) {
                        if (cm->sg_frame_ep_cnt[i] > best_ep_cnt) {
                            best_ep     = i;
                            best_ep_cnt = cm->sg_frame_ep_cnt[i];
                        }
                    }
                    cm->sg_frame_ep = best_ep;
                }
            } else {
                pcs_ptr->rst_info[0].frame_restoration_type = RESTORE_NONE;
                pcs_ptr->rst_info[1].frame_restoration_type = RESTORE_NONE;
                pcs_ptr->rst_info[2].frame_restoration_type = RESTORE_NONE;
            }

            // delete scaled_input_picture_ptr after lr finished
            EB_DELETE(pcs_ptr->scaled_input_picture_ptr);
            if (pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL) {
                // copy stat to ref object (intra_coded_area, Luminance, Scene change detection flags)
                copy_statistics_to_ref_obj_ect(pcs_ptr, scs_ptr);
            }

            superres_recode = pcs_ptr->parent_pcs_ptr->superres_total_recode_loop > 0 ? TRUE
                                                                                      : FALSE;

            // Pad the reference picture and set ref POC
            if (scs_ptr->static_config.pass != ENC_FIRST_PASS) {
                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == TRUE)
                    pad_ref_and_set_flags(pcs_ptr, scs_ptr);
                else {
                    // convert non-reference frame buffer from 16-bit to 8-bit, to export recon and psnr/ssim calculation
                    if (is_16bit && scs_ptr->static_config.encoder_bit_depth == EB_EIGHT_BIT) {
                        EbPictureBufferDesc *ref_pic_ptr =
                            pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture_ptr;
                        EbPictureBufferDesc *ref_pic_16bit_ptr =
                            pcs_ptr->parent_pcs_ptr->enc_dec_ptr->recon_picture16bit_ptr;
                        //Y
                        uint16_t *buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_y);
                        uint8_t  *buf_8bit  = ref_pic_ptr->buffer_y;
                        svt_convert_16bit_to_8bit(
                            buf_16bit,
                            ref_pic_16bit_ptr->stride_y,
                            buf_8bit,
                            ref_pic_ptr->stride_y,
                            ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1),
                            ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1));

                        //CB
                        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cb);
                        buf_8bit  = ref_pic_ptr->buffer_cb;
                        svt_convert_16bit_to_8bit(
                            buf_16bit,
                            ref_pic_16bit_ptr->stride_cb,
                            buf_8bit,
                            ref_pic_ptr->stride_cb,
                            (ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1)) >>
                                scs_ptr->subsampling_x,
                            (ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1)) >>
                                scs_ptr->subsampling_y);

                        //CR
                        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cr);
                        buf_8bit  = ref_pic_ptr->buffer_cr;
                        svt_convert_16bit_to_8bit(
                            buf_16bit,
                            ref_pic_16bit_ptr->stride_cr,
                            buf_8bit,
                            ref_pic_ptr->stride_cr,
                            (ref_pic_16bit_ptr->width + (ref_pic_ptr->origin_x << 1)) >>
                                scs_ptr->subsampling_x,
                            (ref_pic_16bit_ptr->height + (ref_pic_ptr->origin_y << 1)) >>
                                scs_ptr->subsampling_y);
                    }
                }
            }

            // PSNR and SSIM Calculation.
            if (superres_recode) { // superres needs psnr to compute rdcost
                // Note: if superres recode is actived, memory needs to be freed in packetization process by calling free_temporal_filtering_buffer()
                EbErrorType return_error = psnr_calculations(pcs_ptr, scs_ptr, FALSE);
                if (return_error != EB_ErrorNone) {
                    assert_err(0,
                               "Couldn't allocate memory for uncompressed 10bit buffers for PSNR "
                               "calculations");
                }
            } else if (scs_ptr->static_config.stat_report) {
                // Note: if temporal_filtering is used, memory needs to be freed in the last of these calls
                EbErrorType return_error = psnr_calculations(pcs_ptr, scs_ptr, FALSE);
                if (return_error != EB_ErrorNone) {
                    assert_err(0,
                               "Couldn't allocate memory for uncompressed 10bit buffers for PSNR "
                               "calculations");
                }
                return_error = ssim_calculations(pcs_ptr, scs_ptr, TRUE /* free memory here */);
                if (return_error != EB_ErrorNone) {
                    assert_err(0,
                               "Couldn't allocate memory for uncompressed 10bit buffers for SSIM "
                               "calculations");
                }
            }

            if (!superres_recode) {
                if (scs_ptr->static_config.recon_enabled) {
                    recon_output(pcs_ptr, scs_ptr);
                }
                // post reference picture task in packetization process if it's superres_recode
                if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) {
                    // Get Empty PicMgr Results
                    svt_get_empty_object(context_ptr->picture_demux_fifo_ptr,
                                         &picture_demux_results_wrapper_ptr);

                    picture_demux_results_rtr = (PictureDemuxResults *)
                                                    picture_demux_results_wrapper_ptr->object_ptr;
                    picture_demux_results_rtr->reference_picture_wrapper_ptr =
                        pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr;
                    picture_demux_results_rtr->scs_ptr        = pcs_ptr->scs_ptr;
                    picture_demux_results_rtr->picture_number = pcs_ptr->picture_number;
                    picture_demux_results_rtr->picture_type   = EB_PIC_REFERENCE;

                    // Post Reference Picture
                    svt_post_full_object(picture_demux_results_wrapper_ptr);
                }
            }

            tile_cols = pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_cols;
            tile_rows = pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_rows;

            for (int tile_row_idx = 0; tile_row_idx < tile_rows; tile_row_idx++) {
                for (int tile_col_idx = 0; tile_col_idx < tile_cols; tile_col_idx++) {
                    const int tile_idx = tile_row_idx * tile_cols + tile_col_idx;
                    svt_get_empty_object(context_ptr->rest_output_fifo_ptr,
                                         &rest_results_wrapper_ptr);
                    rest_results_ptr = (struct RestResults *)rest_results_wrapper_ptr->object_ptr;
                    rest_results_ptr->pcs_wrapper_ptr = cdef_results_ptr->pcs_wrapper_ptr;
                    rest_results_ptr->tile_index      = tile_idx;
                    // Post Rest Results
                    svt_post_full_object(rest_results_wrapper_ptr);
                }
            }
        }
        svt_release_mutex(pcs_ptr->rest_search_mutex);

        // Release input Results
        svt_release_object(cdef_results_wrapper_ptr);
    }

    return NULL;
}
