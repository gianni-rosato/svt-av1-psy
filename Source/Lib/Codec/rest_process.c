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

#include "enc_handle.h"
#include "rest_process.h"
#include "enc_dec_results.h"
#include "svt_threads.h"
#include "pic_demux_results.h"
#include "reference_object.h"
#include "pcs.h"
#include "resource_coordination_process.h"
#include "resize.h"
#include "enc_mode_config.h"

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

void        svt_aom_pack_highbd_pic(const EbPictureBufferDesc *pic_ptr, uint16_t *buffer_16bit[3], uint32_t ss_x,
                                    uint32_t ss_y, Bool include_padding);
void        svt_aom_copy_buffer_info(EbPictureBufferDesc *src_ptr, EbPictureBufferDesc *dst_ptr);
void        svt_aom_recon_output(PictureControlSet *pcs, SequenceControlSet *scs);
void        svt_av1_loop_restoration_filter_frame(int32_t *rst_tmpbuf, Yv12BufferConfig *frame, Av1Common *cm,
                                                  int32_t optimized_lr);
EbErrorType psnr_calculations(PictureControlSet *pcs, SequenceControlSet *scs, Bool free_memory);
EbErrorType svt_aom_ssim_calculations(PictureControlSet *pcs, SequenceControlSet *scs, Bool free_memory);
void        pad_ref_and_set_flags(PictureControlSet *pcs, SequenceControlSet *scs);
void        restoration_seg_search(int32_t *rst_tmpbuf, Yv12BufferConfig *org_fts, const Yv12BufferConfig *src,
                                   Yv12BufferConfig *trial_frame_rst, PictureControlSet *pcs, uint32_t segment_index);
void        rest_finish_search(PictureControlSet *pcs);
void        svt_av1_upscale_normative_rows(const Av1Common *cm, const uint8_t *src, int src_stride, uint8_t *dst,
                                           int dst_stride, int rows, int sub_x, int bd, Bool is_16bit_pipeline);
#if DEBUG_UPSCALING
void save_YUV_to_file(char *filename, EbByte buffer_y, EbByte buffer_u, EbByte buffer_v, uint16_t width,
                      uint16_t height, uint16_t stride_y, uint16_t stride_u, uint16_t stride_v, uint16_t org_y,
                      uint16_t org_x, uint32_t ss_x, uint32_t ss_y);
#endif

static void rest_context_dctor(EbPtr p) {
    EbThreadContext *thread_ctx = (EbThreadContext *)p;
    RestContext     *obj        = (RestContext *)thread_ctx->priv;
    EB_DELETE(obj->trial_frame_rst);
    // buffer only malloc'd if boundaries are used in rest. search.
    // see scs->seq_header.use_boundaries_in_rest_search
    if (obj->org_rec_frame)
        EB_DELETE(obj->org_rec_frame);
    EB_FREE_ALIGNED(obj->rst_tmpbuf);
    EB_FREE_ARRAY(obj);
}

/******************************************************
 * Rest Context Constructor
 ******************************************************/
EbErrorType svt_aom_rest_context_ctor(EbThreadContext *thread_ctx, const EbEncHandle *enc_handle_ptr,
                                      EbPtr object_init_data_ptr, int index, int demux_index) {
    const SequenceControlSet       *scs           = enc_handle_ptr->scs_instance_array[0]->scs;
    const EbSvtAv1EncConfiguration *config        = &scs->static_config;
    EbColorFormat                   color_format  = config->encoder_color_format;
    EbPictureBufferDescInitData    *init_data_ptr = (EbPictureBufferDescInitData *)object_init_data_ptr;
    RestContext                    *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_ctx->priv  = context_ptr;
    thread_ctx->dctor = rest_context_dctor;

    // Input/Output System Resource Manager FIFOs
    context_ptr->rest_input_fifo_ptr  = svt_system_resource_get_consumer_fifo(enc_handle_ptr->cdef_results_resource_ptr,
                                                                             index);
    context_ptr->rest_output_fifo_ptr = svt_system_resource_get_producer_fifo(enc_handle_ptr->rest_results_resource_ptr,
                                                                              index);
    context_ptr->picture_demux_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);

    Bool is_16bit = scs->is_16bit_pipeline;
    if (svt_aom_get_enable_restoration(init_data_ptr->enc_mode,
                                       config->enable_restoration_filtering,
                                       scs->input_resolution,
                                       config->fast_decode)) {
        EbPictureBufferDescInitData init_data;

        init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        init_data.max_width          = (uint16_t)scs->max_input_luma_width;
        init_data.max_height         = (uint16_t)scs->max_input_luma_height;
        init_data.bit_depth          = is_16bit ? EB_SIXTEEN_BIT : EB_EIGHT_BIT;
        init_data.color_format       = color_format;
        init_data.left_padding       = AOM_RESTORATION_FRAME_BORDER;
        init_data.right_padding      = AOM_RESTORATION_FRAME_BORDER;
        init_data.top_padding        = AOM_RESTORATION_FRAME_BORDER;
        init_data.bot_padding        = AOM_RESTORATION_FRAME_BORDER;
        init_data.split_mode         = FALSE;
        init_data.is_16bit_pipeline  = is_16bit;

        EB_NEW(context_ptr->trial_frame_rst, svt_picture_buffer_desc_ctor, (EbPtr)&init_data);
        if (scs->use_boundaries_in_rest_search)
            EB_NEW(context_ptr->org_rec_frame, svt_picture_buffer_desc_ctor, (EbPtr)&init_data);
        else
            context_ptr->org_rec_frame = NULL;
        if (!is_16bit) {
            context_ptr->trial_frame_rst->bit_depth = EB_EIGHT_BIT;
            if (scs->use_boundaries_in_rest_search)
                context_ptr->org_rec_frame->bit_depth = EB_EIGHT_BIT;
        }
        context_ptr->rst_tmpbuf = NULL;
        if (svt_aom_get_enable_sg(init_data_ptr->enc_mode, scs->input_resolution, config->fast_decode))
            EB_MALLOC_ALIGNED(context_ptr->rst_tmpbuf, RESTORATION_TMPBUF_SIZE);
    }

    return EB_ErrorNone;
}
extern void svt_aom_get_recon_pic(PictureControlSet *pcs, EbPictureBufferDesc **recon_ptr, Bool is_highbd) {
    if (!is_highbd) {
        if (pcs->ppcs->is_ref == TRUE)
            *recon_ptr = ((EbReferenceObject *)pcs->ppcs->ref_pic_wrapper->object_ptr)->reference_picture;
        else
            *recon_ptr = pcs->ppcs->enc_dec_ptr->recon_pic; // OMK
    } else {
        *recon_ptr = pcs->ppcs->enc_dec_ptr->recon_pic_16bit;
    }

    // recon buffer is created in full resolution, it is resized to difference size
    // when reference scaling enabled. recon width and height should be adjusted to
    // upscaled render size
    if (*recon_ptr &&
        (pcs->ppcs->render_width != (*recon_ptr)->width || pcs->ppcs->render_height != (*recon_ptr)->height)) {
        (*recon_ptr)->width  = pcs->ppcs->render_width;
        (*recon_ptr)->height = pcs->ppcs->render_height;
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
static EbPictureBufferDesc *get_own_recon(SequenceControlSet *scs, PictureControlSet *pcs, RestContext *context_ptr,
                                          Bool is_16bit) {
    const uint32_t ss_x = scs->subsampling_x;
    const uint32_t ss_y = scs->subsampling_y;

    EbPictureBufferDesc *recon_pic;
    svt_aom_get_recon_pic(pcs, &recon_pic, is_16bit);
    // if boundaries are not used, don't need to copy pic to new buffer, as the
    // search will not modify the pic
    if (!scs->use_boundaries_in_rest_search) {
        return recon_pic;
    }
    if (is_16bit) {
        uint16_t *rec_ptr = (uint16_t *)recon_pic->buffer_y + recon_pic->org_x + recon_pic->org_y * recon_pic->stride_y;
        uint16_t *rec_ptr_cb = (uint16_t *)recon_pic->buffer_cb + recon_pic->org_x / 2 +
            recon_pic->org_y / 2 * recon_pic->stride_cb;
        uint16_t *rec_ptr_cr = (uint16_t *)recon_pic->buffer_cr + recon_pic->org_x / 2 +
            recon_pic->org_y / 2 * recon_pic->stride_cr;

        EbPictureBufferDesc *org_rec = context_ptr->org_rec_frame;
        uint16_t *org_ptr    = (uint16_t *)org_rec->buffer_y + org_rec->org_x + org_rec->org_y * org_rec->stride_y;
        uint16_t *org_ptr_cb = (uint16_t *)org_rec->buffer_cb + org_rec->org_x / 2 +
            org_rec->org_y / 2 * org_rec->stride_cb;
        uint16_t *org_ptr_cr = (uint16_t *)org_rec->buffer_cr + org_rec->org_x / 2 +
            org_rec->org_y / 2 * org_rec->stride_cr;

        for (int r = 0; r < recon_pic->height; ++r)
            svt_memcpy(org_ptr + r * org_rec->stride_y, rec_ptr + r * recon_pic->stride_y, recon_pic->width << 1);

        for (int r = 0; r < (recon_pic->height >> ss_y); ++r) {
            svt_memcpy(org_ptr_cb + r * org_rec->stride_cb,
                       rec_ptr_cb + r * recon_pic->stride_cb,
                       (recon_pic->width >> ss_x) << 1);
            svt_memcpy(org_ptr_cr + r * org_rec->stride_cr,
                       rec_ptr_cr + r * recon_pic->stride_cr,
                       (recon_pic->width >> ss_x) << 1);
        }
    } else {
        if (pcs->ppcs->is_ref == TRUE)
            recon_pic = ((EbReferenceObject *)pcs->ppcs->ref_pic_wrapper->object_ptr)->reference_picture;
        else
            recon_pic = pcs->ppcs->enc_dec_ptr->recon_pic; // OMK
        // if boundaries are not used, don't need to copy pic to new buffer, as the
        // search will not modify the pic
        if (!scs->use_boundaries_in_rest_search) {
            return recon_pic;
        }
        uint8_t *rec_ptr    = &((recon_pic->buffer_y)[recon_pic->org_x + recon_pic->org_y * recon_pic->stride_y]);
        uint8_t *rec_ptr_cb = &(
            (recon_pic->buffer_cb)[recon_pic->org_x / 2 + recon_pic->org_y / 2 * recon_pic->stride_cb]);
        uint8_t *rec_ptr_cr = &(
            (recon_pic->buffer_cr)[recon_pic->org_x / 2 + recon_pic->org_y / 2 * recon_pic->stride_cr]);

        EbPictureBufferDesc *org_rec = context_ptr->org_rec_frame;
        uint8_t             *org_ptr = &((org_rec->buffer_y)[org_rec->org_x + org_rec->org_y * org_rec->stride_y]);
        uint8_t *org_ptr_cb = &((org_rec->buffer_cb)[org_rec->org_x / 2 + org_rec->org_y / 2 * org_rec->stride_cb]);
        uint8_t *org_ptr_cr = &((org_rec->buffer_cr)[org_rec->org_x / 2 + org_rec->org_y / 2 * org_rec->stride_cr]);

        for (int r = 0; r < recon_pic->height; ++r)
            svt_memcpy(org_ptr + r * org_rec->stride_y, rec_ptr + r * recon_pic->stride_y, recon_pic->width);

        for (int r = 0; r < (recon_pic->height >> ss_y); ++r) {
            svt_memcpy(
                org_ptr_cb + r * org_rec->stride_cb, rec_ptr_cb + r * recon_pic->stride_cb, (recon_pic->width >> ss_x));
            svt_memcpy(
                org_ptr_cr + r * org_rec->stride_cr, rec_ptr_cr + r * recon_pic->stride_cr, (recon_pic->width >> ss_x));
        }
    }
    return context_ptr->org_rec_frame;
}

void svt_convert_pic_8bit_to_16bit(EbPictureBufferDesc *src_8bit, EbPictureBufferDesc *dst_16bit, uint16_t ss_x,
                                   uint16_t ss_y) {
    //copy input from 8bit to 16bit
    uint8_t  *buffer_8bit;
    int32_t   stride_8bit;
    uint16_t *buffer_16bit;
    int32_t   stride_16bit;
    // Y
    buffer_16bit = (uint16_t *)(dst_16bit->buffer_y) + dst_16bit->org_x + dst_16bit->org_y * dst_16bit->stride_y;
    stride_16bit = dst_16bit->stride_y;
    buffer_8bit  = src_8bit->buffer_y + src_8bit->org_x + src_8bit->org_y * src_8bit->stride_y;
    stride_8bit  = src_8bit->stride_y;

    svt_convert_8bit_to_16bit(buffer_8bit, stride_8bit, buffer_16bit, stride_16bit, src_8bit->width, src_8bit->height);

    // Cb
    buffer_16bit = (uint16_t *)(dst_16bit->buffer_cb) + (dst_16bit->org_x >> ss_x) +
        (dst_16bit->org_y >> ss_y) * dst_16bit->stride_cb;
    stride_16bit = dst_16bit->stride_cb;
    buffer_8bit  = src_8bit->buffer_cb + (src_8bit->org_x >> ss_x) + (src_8bit->org_y >> ss_y) * src_8bit->stride_cb;
    stride_8bit  = src_8bit->stride_cb;

    svt_convert_8bit_to_16bit(
        buffer_8bit, stride_8bit, buffer_16bit, stride_16bit, src_8bit->width >> ss_x, src_8bit->height >> ss_y);

    // Cr
    buffer_16bit = (uint16_t *)(dst_16bit->buffer_cr) + (dst_16bit->org_x >> ss_x) +
        (dst_16bit->org_y >> ss_y) * dst_16bit->stride_cr;
    stride_16bit = dst_16bit->stride_cr;
    buffer_8bit  = src_8bit->buffer_cr + (src_8bit->org_x >> ss_x) + (src_8bit->org_y >> ss_y) * src_8bit->stride_cr;
    stride_8bit  = src_8bit->stride_cr;

    svt_convert_8bit_to_16bit(
        buffer_8bit, stride_8bit, buffer_16bit, stride_16bit, src_8bit->width >> ss_x, src_8bit->height >> ss_y);

    dst_16bit->width  = src_8bit->width;
    dst_16bit->height = src_8bit->height;
}

extern void svt_aom_pack_2d_pic(EbPictureBufferDesc *input_picture, uint16_t *packed[3]);

void set_unscaled_input_16bit(PictureControlSet *pcs) {
    EbPictureBufferDesc *input_pic  = pcs->ppcs->enhanced_unscaled_pic;
    EbPictureBufferDesc *output_pic = pcs->input_frame16bit;
    uint16_t             ss_x       = pcs->ppcs->scs->subsampling_x;
    uint16_t             ss_y       = pcs->ppcs->scs->subsampling_y;
    svt_aom_copy_buffer_info(input_pic, pcs->input_frame16bit);
    if (input_pic->bit_depth == EB_EIGHT_BIT)
        svt_convert_pic_8bit_to_16bit(input_pic, output_pic, ss_x, ss_y);
    else {
        uint16_t *planes[3] = {
            (uint16_t *)output_pic->buffer_y + (output_pic->org_y * output_pic->stride_y) + (output_pic->org_x),
            (uint16_t *)output_pic->buffer_cb + (((output_pic->org_y) >> ss_y) * output_pic->stride_cb) +
                ((output_pic->org_x) >> ss_x),
            (uint16_t *)output_pic->buffer_cr + (((output_pic->org_y) >> ss_y) * output_pic->stride_cr) +
                ((output_pic->org_x) >> ss_x)};
        svt_aom_pack_2d_pic(input_pic, planes);
    }
}

static void derive_blk_pointers_enc(EbPictureBufferDesc *recon_picture_buf, int32_t plane, int32_t blk_col_px,
                                    int32_t blk_row_px, void **pp_blk_recon_buf, int32_t *recon_stride, int32_t sub_x,
                                    int32_t sub_y, Bool use_highbd) {
    int32_t block_offset;

    if (plane == 0) {
        block_offset = (recon_picture_buf->org_y + blk_row_px) * recon_picture_buf->stride_y +
            (recon_picture_buf->org_x + blk_col_px);
        *recon_stride = recon_picture_buf->stride_y;
    } else if (plane == 1) {
        block_offset = ((recon_picture_buf->org_y >> sub_y) + blk_row_px) * recon_picture_buf->stride_cb +
            ((recon_picture_buf->org_x >> sub_x) + blk_col_px);
        *recon_stride = recon_picture_buf->stride_cb;
    } else {
        block_offset = ((recon_picture_buf->org_y >> sub_y) + blk_row_px) * recon_picture_buf->stride_cr +
            ((recon_picture_buf->org_x >> sub_x) + blk_col_px);
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

static EbErrorType copy_recon_enc(SequenceControlSet *scs, EbPictureBufferDesc *recon_picture_src,
                                  EbPictureBufferDesc *recon_picture_dst, int num_planes, int skip_copy) {
    recon_picture_dst->org_x        = recon_picture_src->org_x;
    recon_picture_dst->org_y        = recon_picture_src->org_y;
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

    recon_picture_dst->buffer_enable_mask = scs->seq_header.color_config.mono_chrome ? PICTURE_BUFFER_DESC_LUMA_MASK
                                                                                     : PICTURE_BUFFER_DESC_FULL_MASK;

    uint32_t bytesPerPixel = scs->is_16bit_pipeline ? 2 : 1;

    // Allocate the Picture Buffers (luma & chroma)
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_y, recon_picture_dst->luma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_y, 0, recon_picture_dst->luma_size * bytesPerPixel);
    } else
        recon_picture_dst->buffer_y = 0;
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_cb, recon_picture_dst->chroma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_cb, 0, recon_picture_dst->chroma_size * bytesPerPixel);
    } else
        recon_picture_dst->buffer_cb = 0;
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_cr, recon_picture_dst->chroma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_cr, 0, recon_picture_dst->chroma_size * bytesPerPixel);
    } else
        recon_picture_dst->buffer_cr = 0;

    int use_highbd = scs->is_16bit_pipeline;

    if (!skip_copy) {
        for (int plane = 0; plane < num_planes; ++plane) {
            uint8_t *src_buf, *dst_buf;
            int32_t  src_stride, dst_stride;

            int sub_x = plane ? scs->subsampling_x : 0;
            int sub_y = plane ? scs->subsampling_y : 0;

            derive_blk_pointers_enc(
                recon_picture_src, plane, 0, 0, (void *)&src_buf, &src_stride, sub_x, sub_y, use_highbd);
            derive_blk_pointers_enc(
                recon_picture_dst, plane, 0, 0, (void *)&dst_buf, &dst_stride, sub_x, sub_y, use_highbd);

            int height = ((recon_picture_src->height + sub_y) >> sub_y);
            for (int row = 0; row < height; ++row) {
                svt_memcpy(
                    dst_buf, src_buf, ((recon_picture_src->width + sub_x) >> sub_x) * sizeof(*src_buf) << use_highbd);
                src_buf += src_stride << use_highbd;
                dst_buf += dst_stride << use_highbd;
            }
        }
    }

    return EB_ErrorNone;
}

void svt_av1_superres_upscale_frame(struct Av1Common *cm, PictureControlSet *pcs, SequenceControlSet *scs) {
    // Set these parameters for testing since they are not correctly populated yet
    EbPictureBufferDesc *recon_ptr;

    Bool is_16bit = scs->is_16bit_pipeline;

    svt_aom_get_recon_pic(pcs, &recon_ptr, is_16bit);

    uint16_t  ss_x       = scs->subsampling_x;
    uint16_t  ss_y       = scs->subsampling_y;
    const int num_planes = scs->seq_header.color_config.mono_chrome ? 1 : MAX_MB_PLANE;

    EbPictureBufferDesc  recon_pic_temp;
    EbPictureBufferDesc *ps_recon_pic_temp;
    ps_recon_pic_temp = &recon_pic_temp;

    EbErrorType return_error = copy_recon_enc(scs, recon_ptr, ps_recon_pic_temp, num_planes, 0);

    if (return_error != EB_ErrorNone) {
        ps_recon_pic_temp = NULL;
        assert(0);
    }

    EbPictureBufferDesc *src = ps_recon_pic_temp;
    EbPictureBufferDesc *dst = recon_ptr;

    // get the bit-depth from the encoder config instead of from the recon ptr
    int bit_depth = scs->static_config.encoder_bit_depth;

    for (int plane = 0; plane < num_planes; ++plane) {
        uint8_t *src_buf, *dst_buf;
        int32_t  src_stride, dst_stride;

        int sub_x = plane ? ss_x : 0;
        int sub_y = plane ? ss_y : 0;
        derive_blk_pointers_enc(src, plane, 0, 0, (void *)&src_buf, &src_stride, sub_x, sub_y, is_16bit);
        derive_blk_pointers_enc(dst, plane, 0, 0, (void *)&dst_buf, &dst_stride, sub_x, sub_y, is_16bit);

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
    EbReferenceObject *obj = (EbReferenceObject *)pcs->ppcs->ref_pic_wrapper->object_ptr;

    pcs->intra_coded_area = (100 * pcs->intra_coded_area) / (pcs->ppcs->aligned_width * pcs->ppcs->aligned_height);
    pcs->skip_coded_area  = (100 * pcs->skip_coded_area) / (pcs->ppcs->aligned_width * pcs->ppcs->aligned_height);
    pcs->hp_coded_area    = (100 * pcs->hp_coded_area) / (pcs->ppcs->aligned_width * pcs->ppcs->aligned_height);
    if (pcs->slice_type == I_SLICE)
        pcs->intra_coded_area = 0;
    obj->intra_coded_area                   = (uint8_t)(pcs->intra_coded_area);
    obj->skip_coded_area                    = (uint8_t)(pcs->skip_coded_area);
    obj->hp_coded_area                      = (uint8_t)(pcs->hp_coded_area);
    struct PictureParentControlSet *ppcs    = pcs->ppcs;
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
        obj->sb_me_64x64_dist[sb_index]   = pcs->ppcs->me_64x64_distortion[sb_index];
        obj->sb_me_8x8_cost_var[sb_index] = pcs->ppcs->me_8x8_cost_variance[sb_index];
    }
    obj->tmp_layer_idx   = (uint8_t)pcs->temporal_layer_index;
    obj->is_scene_change = pcs->ppcs->scene_change_flag;

    Av1Common *cm    = pcs->ppcs->av1_cm;
    obj->sg_frame_ep = cm->sg_frame_ep;
    if (scs->mfmv_enabled || !pcs->ppcs->is_not_scaled) {
        obj->frame_type = pcs->ppcs->frm_hdr.frame_type;
        obj->order_hint = pcs->ppcs->cur_order_hint;
        svt_memcpy(obj->ref_order_hint, pcs->ppcs->ref_order_hint, 7 * sizeof(uint32_t));
    }
    // Copy the prev frame wn filter coeffs
    if (cm->wn_filter_ctrls.enabled && cm->wn_filter_ctrls.use_prev_frame_coeffs) {
        for (int32_t plane = 0; plane < MAX_MB_PLANE; ++plane) {
            int32_t ntiles = pcs->rst_info[plane].units_per_tile;
            for (int32_t u = 0; u < ntiles; ++u) {
                obj->unit_info[plane][u].restoration_type = pcs->rst_info[plane].unit_info[u].restoration_type;
                if (pcs->rst_info[plane].unit_info[u].restoration_type == RESTORE_WIENER)
                    obj->unit_info[plane][u].wiener_info = pcs->rst_info[plane].unit_info[u].wiener_info;
            }
        }
    }
}

/******************************************************
 * Rest Kernel
 ******************************************************/
void *svt_aom_rest_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext    *thread_ctx  = (EbThreadContext *)input_ptr;
    RestContext        *context_ptr = (RestContext *)thread_ctx->priv;
    PictureControlSet  *pcs;
    SequenceControlSet *scs;

    //// Input
    EbObjectWrapper *cdef_results_wrapper;
    CdefResults     *cdef_results;

    //// Output
    EbObjectWrapper     *rest_results_wrapper;
    RestResults         *rest_results;
    EbObjectWrapper     *picture_demux_results_wrapper_ptr;
    PictureDemuxResults *picture_demux_results_rtr;
    // SB Loop variables
    uint8_t tile_cols;
    uint8_t tile_rows;

    Bool superres_recode = FALSE;

    for (;;) {
        // Get Cdef Results
        EB_GET_FULL_OBJECT(context_ptr->rest_input_fifo_ptr, &cdef_results_wrapper);

        cdef_results                  = (CdefResults *)cdef_results_wrapper->object_ptr;
        pcs                           = (PictureControlSet *)cdef_results->pcs_wrapper->object_ptr;
        PictureParentControlSet *ppcs = pcs->ppcs;
        scs                           = pcs->scs;
        FrameHeader *frm_hdr          = &pcs->ppcs->frm_hdr;
        Bool         is_16bit         = scs->is_16bit_pipeline;
        Av1Common   *cm               = pcs->ppcs->av1_cm;
        if (ppcs->enable_restoration && frm_hdr->allow_intrabc == 0) {
            // If using boundaries during the filter search, copy the recon pic to a new buffer (to
            // avoid race condition from many threads modifying the same recon pic).
            //
            // If not using boundaries during the filter search, copy the input recon picture
            // location to be used in restoration search (save cycles/memory of copying pic to a new
            // buffer). The recon pic should not be modified during the search, otherwise there will
            // be a race condition between threads.
            EbPictureBufferDesc *recon_pic = get_own_recon(scs, pcs, context_ptr, is_16bit);
            EbPictureBufferDesc *input_pic = is_16bit ? pcs->input_frame16bit : pcs->ppcs->enhanced_unscaled_pic;

            EbPictureBufferDesc *scaled_input_pic = NULL;
            // downscale input picture if recon is resized
            Bool is_resized = recon_pic->width != input_pic->width || recon_pic->height != input_pic->height;
            if (is_resized) {
                scaled_input_pic = pcs->scaled_input_pic;
                input_pic        = scaled_input_pic;
            }

            // there are padding pixels if input pics are not 8 pixel aligned
            // but there is no extra padding after input pics are resized for
            // reference scaling
            Yv12BufferConfig cpi_source;
            svt_aom_link_eb_to_aom_buffer_desc(input_pic,
                                               &cpi_source,
                                               is_resized ? 0 : scs->max_input_pad_right,
                                               is_resized ? 0 : scs->max_input_pad_bottom,
                                               is_16bit);

            Yv12BufferConfig trial_frame_rst;
            svt_aom_link_eb_to_aom_buffer_desc(context_ptr->trial_frame_rst,
                                               &trial_frame_rst,
                                               is_resized ? 0 : scs->max_input_pad_right,
                                               is_resized ? 0 : scs->max_input_pad_bottom,
                                               is_16bit);

            Yv12BufferConfig org_fts;
            svt_aom_link_eb_to_aom_buffer_desc(recon_pic,
                                               &org_fts,
                                               is_resized ? 0 : scs->max_input_pad_right,
                                               is_resized ? 0 : scs->max_input_pad_bottom,
                                               is_16bit);

            if (pcs->ppcs->slice_type != I_SLICE && cm->wn_filter_ctrls.enabled &&
                cm->wn_filter_ctrls.use_prev_frame_coeffs) {
                EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
                for (int32_t plane = 0; plane < MAX_MB_PLANE; ++plane) {
                    int32_t ntiles = pcs->rst_info[plane].units_per_tile;
                    for (int32_t u = 0; u < ntiles; ++u) {
                        pcs->rst_info[plane].unit_info[u].restoration_type =
                            ref_obj_l0->unit_info[plane][u].restoration_type;
                        if (ref_obj_l0->unit_info[plane][u].restoration_type == RESTORE_WIENER)
                            pcs->rst_info[plane].unit_info[u].wiener_info = ref_obj_l0->unit_info[plane][u].wiener_info;
                    }
                }
            }
            restoration_seg_search(
                context_ptr->rst_tmpbuf, &org_fts, &cpi_source, &trial_frame_rst, pcs, cdef_results->segment_index);
        }

        //all seg based search is done. update total processed segments. if all done, finish the search and perfrom application.
        svt_block_on_mutex(pcs->rest_search_mutex);

        pcs->tot_seg_searched_rest++;
        if (pcs->tot_seg_searched_rest == pcs->rest_segments_total_count) {
            if (ppcs->enable_restoration && frm_hdr->allow_intrabc == 0) {
                rest_finish_search(pcs);

                // Only need recon if REF pic or recon is output
                if (pcs->ppcs->is_ref || scs->static_config.recon_enabled) {
                    if (pcs->rst_info[0].frame_restoration_type != RESTORE_NONE ||
                        pcs->rst_info[1].frame_restoration_type != RESTORE_NONE ||
                        pcs->rst_info[2].frame_restoration_type != RESTORE_NONE) {
                        svt_av1_loop_restoration_filter_frame(context_ptr->rst_tmpbuf, cm->frame_to_show, cm, 0);
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
                pcs->rst_info[0].frame_restoration_type = RESTORE_NONE;
                pcs->rst_info[1].frame_restoration_type = RESTORE_NONE;
                pcs->rst_info[2].frame_restoration_type = RESTORE_NONE;
            }

            // delete scaled_input_pic after lr finished
            EB_DELETE(pcs->scaled_input_pic);
            if (pcs->ppcs->ref_pic_wrapper != NULL) {
                // copy stat to ref object (intra_coded_area, Luminance, Scene change detection
                // flags)
                copy_statistics_to_ref_obj_ect(pcs, scs);
            }

            superres_recode = pcs->ppcs->superres_total_recode_loop > 0 ? TRUE : FALSE;

            // Pad the reference picture and set ref POC
            {
                if (pcs->ppcs->is_ref == TRUE)
                    pad_ref_and_set_flags(pcs, scs);
                else {
                    // convert non-reference frame buffer from 16-bit to 8-bit, to export recon and
                    // psnr/ssim calculation
                    if (is_16bit && scs->static_config.encoder_bit_depth == EB_EIGHT_BIT) {
                        EbPictureBufferDesc *ref_pic_ptr       = pcs->ppcs->enc_dec_ptr->recon_pic;
                        EbPictureBufferDesc *ref_pic_16bit_ptr = pcs->ppcs->enc_dec_ptr->recon_pic_16bit;
                        // Y
                        uint16_t *buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_y);
                        uint8_t  *buf_8bit  = ref_pic_ptr->buffer_y;
                        svt_convert_16bit_to_8bit(buf_16bit,
                                                  ref_pic_16bit_ptr->stride_y,
                                                  buf_8bit,
                                                  ref_pic_ptr->stride_y,
                                                  ref_pic_16bit_ptr->width + (ref_pic_ptr->org_x << 1),
                                                  ref_pic_16bit_ptr->height + (ref_pic_ptr->org_y << 1));

                        //CB
                        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cb);
                        buf_8bit  = ref_pic_ptr->buffer_cb;
                        svt_convert_16bit_to_8bit(
                            buf_16bit,
                            ref_pic_16bit_ptr->stride_cb,
                            buf_8bit,
                            ref_pic_ptr->stride_cb,
                            (ref_pic_16bit_ptr->width + (ref_pic_ptr->org_x << 1)) >> scs->subsampling_x,
                            (ref_pic_16bit_ptr->height + (ref_pic_ptr->org_y << 1)) >> scs->subsampling_y);

                        //CR
                        buf_16bit = (uint16_t *)(ref_pic_16bit_ptr->buffer_cr);
                        buf_8bit  = ref_pic_ptr->buffer_cr;
                        svt_convert_16bit_to_8bit(
                            buf_16bit,
                            ref_pic_16bit_ptr->stride_cr,
                            buf_8bit,
                            ref_pic_ptr->stride_cr,
                            (ref_pic_16bit_ptr->width + (ref_pic_ptr->org_x << 1)) >> scs->subsampling_x,
                            (ref_pic_16bit_ptr->height + (ref_pic_ptr->org_y << 1)) >> scs->subsampling_y);
                    }
                }
            }

            // PSNR and SSIM Calculation.
            if (superres_recode) { // superres needs psnr to compute rdcost
                // Note: if superres recode is actived, memory needs to be freed in packetization process by calling free_temporal_filtering_buffer()
                EbErrorType return_error = psnr_calculations(pcs, scs, FALSE);
                if (return_error != EB_ErrorNone) {
                    svt_aom_assert_err(0,
                                       "Couldn't allocate memory for uncompressed 10bit buffers for PSNR "
                                       "calculations");
                }
            } else if (scs->static_config.stat_report) {
                // Note: if temporal_filtering is used, memory needs to be freed in the last of these calls
                EbErrorType return_error = psnr_calculations(pcs, scs, FALSE);
                if (return_error != EB_ErrorNone) {
                    svt_aom_assert_err(0,
                                       "Couldn't allocate memory for uncompressed 10bit buffers for PSNR "
                                       "calculations");
                }
                return_error = svt_aom_ssim_calculations(pcs, scs, TRUE /* free memory here */);
                if (return_error != EB_ErrorNone) {
                    svt_aom_assert_err(0,
                                       "Couldn't allocate memory for uncompressed 10bit buffers for SSIM "
                                       "calculations");
                }
            }

            if (!superres_recode) {
                if (scs->static_config.recon_enabled) {
                    svt_aom_recon_output(pcs, scs);
                }
                // post reference picture task in packetization process if it's superres_recode
                if (pcs->ppcs->is_ref) {
                    // Get Empty PicMgr Results
                    svt_get_empty_object(context_ptr->picture_demux_fifo_ptr, &picture_demux_results_wrapper_ptr);

                    picture_demux_results_rtr = (PictureDemuxResults *)picture_demux_results_wrapper_ptr->object_ptr;
                    picture_demux_results_rtr->ref_pic_wrapper = pcs->ppcs->ref_pic_wrapper;
                    picture_demux_results_rtr->scs             = pcs->scs;
                    picture_demux_results_rtr->picture_number  = pcs->picture_number;
                    picture_demux_results_rtr->picture_type    = EB_PIC_REFERENCE;

                    // Post Reference Picture
                    svt_post_full_object(picture_demux_results_wrapper_ptr);
                }
            }

            tile_cols = pcs->ppcs->av1_cm->tiles_info.tile_cols;
            tile_rows = pcs->ppcs->av1_cm->tiles_info.tile_rows;

            for (int tile_row_idx = 0; tile_row_idx < tile_rows; tile_row_idx++) {
                for (int tile_col_idx = 0; tile_col_idx < tile_cols; tile_col_idx++) {
                    const int tile_idx = tile_row_idx * tile_cols + tile_col_idx;
                    svt_get_empty_object(context_ptr->rest_output_fifo_ptr, &rest_results_wrapper);
                    rest_results              = (struct RestResults *)rest_results_wrapper->object_ptr;
                    rest_results->pcs_wrapper = cdef_results->pcs_wrapper;
                    rest_results->tile_index  = tile_idx;
                    // Post Rest Results
                    svt_post_full_object(rest_results_wrapper);
                }
            }
        }
        svt_release_mutex(pcs->rest_search_mutex);

        // Release input Results
        svt_release_object(cdef_results_wrapper);
    }

    return NULL;
}
