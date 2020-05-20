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

#include <stdlib.h>

#include "EbEncHandle.h"
#include "EbRestProcess.h"
#include "EbEncDecResults.h"
#include "EbThreads.h"
#include "EbPictureDemuxResults.h"
#include "EbReferenceObject.h"
#include "EbPictureControlSet.h"

#define DEBUG_UPSCALING 0

/**************************************
 * Rest Context
 **************************************/
typedef struct RestContext {
    EbDctor dctor;
    EbFifo *rest_input_fifo_ptr;
    EbFifo *rest_output_fifo_ptr;
    EbFifo *picture_demux_fifo_ptr;

    EbPictureBufferDesc *trial_frame_rst;

    EbPictureBufferDesc *temp_lf_recon_picture_ptr;
    EbPictureBufferDesc *temp_lf_recon_picture16bit_ptr;

    EbPictureBufferDesc *
        org_rec_frame; // while doing the filtering recon gets updated uisng setup/restore processing_stripe_bounadaries
    // many threads doing the above will result in race condition.
    // each thread will hence have his own copy of recon to work on.
    // later we can have a search version that does not need the exact right recon
    int32_t *rst_tmpbuf;
} RestContext;

void pack_highbd_pic(const EbPictureBufferDesc *pic_ptr, uint16_t *buffer_16bit[3], uint32_t ss_x,
                     uint32_t ss_y, EbBool include_padding);
void copy_buffer_info(EbPictureBufferDesc *src_ptr, EbPictureBufferDesc *dst_ptr);
void recon_output(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr);
void eb_av1_loop_restoration_filter_frame(Yv12BufferConfig *frame, Av1Common *cm,
                                          int32_t optimized_lr);
void copy_statistics_to_ref_obj_ect(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr);
void psnr_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory);
void ssim_calculations(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr, EbBool free_memory);
void pad_ref_and_set_flags(PictureControlSet *pcs_ptr, SequenceControlSet *scs_ptr);
void generate_padding(EbByte src_pic, uint32_t src_stride, uint32_t original_src_width,
                      uint32_t original_src_height, uint32_t padding_width,
                      uint32_t padding_height);
void restoration_seg_search(int32_t *rst_tmpbuf, Yv12BufferConfig *org_fts,
                            const Yv12BufferConfig *src, Yv12BufferConfig *trial_frame_rst,
                            PictureControlSet *pcs_ptr, uint32_t segment_index);
void rest_finish_search(PictureParentControlSet *p_pcs_ptr, Macroblock *x, Av1Common *const cm);

void av1_upscale_normative_rows(const Av1Common *cm, const uint8_t *src,
                                int src_stride, uint8_t *dst, int dst_stride, int rows, int sub_x, int bd, EbBool is_16bit_pipeline);

#if DEBUG_UPSCALING
void save_YUV_to_file(char *filename, EbByte buffer_y, EbByte buffer_u, EbByte buffer_v,
                      uint16_t width, uint16_t height,
                      uint16_t stride_y, uint16_t stride_u, uint16_t stride_v,
                      uint16_t origin_y, uint16_t origin_x,
                      uint32_t ss_x, uint32_t ss_y);
#endif

static void rest_context_dctor(EbPtr p)
{
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    RestContext *    obj                = (RestContext *)thread_context_ptr->priv;
    EB_DELETE(obj->temp_lf_recon_picture_ptr);
    EB_DELETE(obj->temp_lf_recon_picture16bit_ptr);
    EB_DELETE(obj->trial_frame_rst);
    EB_DELETE(obj->org_rec_frame);
    EB_FREE_ALIGNED(obj->rst_tmpbuf);
    EB_FREE_ARRAY(obj);
}

/******************************************************
 * Rest Context Constructor
 ******************************************************/
EbErrorType rest_context_ctor(EbThreadContext *  thread_context_ptr,
                              const EbEncHandle *enc_handle_ptr, int index, int demux_index) {
    const SequenceControlSet *      scs_ptr      = enc_handle_ptr->scs_instance_array[0]->scs_ptr;
    const EbSvtAv1EncConfiguration *config       = &scs_ptr->static_config;
    EbBool                          is_16bit     = (EbBool)(config->encoder_bit_depth > EB_8BIT);
    EbColorFormat                   color_format = config->encoder_color_format;

    RestContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = rest_context_dctor;

    // Input/Output System Resource Manager FIFOs
    context_ptr->rest_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->cdef_results_resource_ptr, index);
    context_ptr->rest_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->rest_results_resource_ptr, index);
    context_ptr->picture_demux_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);

    {
        EbPictureBufferDescInitData init_data;

        init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        init_data.max_width          = (uint16_t)scs_ptr->max_input_luma_width;
        init_data.max_height         = (uint16_t)scs_ptr->max_input_luma_height;
        init_data.bit_depth          = config->is_16bit_pipeline || is_16bit
                                     ? EB_16BIT : EB_8BIT;
        init_data.color_format       = color_format;
        init_data.left_padding       = AOM_BORDER_IN_PIXELS;
        init_data.right_padding      = AOM_BORDER_IN_PIXELS;
        init_data.top_padding        = AOM_BORDER_IN_PIXELS;
        init_data.bot_padding        = AOM_BORDER_IN_PIXELS;
        init_data.split_mode         = EB_FALSE;
        init_data.is_16bit_pipeline = config->is_16bit_pipeline;

        EB_NEW(context_ptr->trial_frame_rst, eb_picture_buffer_desc_ctor, (EbPtr)&init_data);

        EB_NEW(context_ptr->org_rec_frame, eb_picture_buffer_desc_ctor, (EbPtr)&init_data);
        if (!is_16bit)
        {
            context_ptr->trial_frame_rst->bit_depth = EB_8BIT;
            context_ptr->org_rec_frame->bit_depth = EB_8BIT;
        }

        EB_MALLOC_ALIGNED(context_ptr->rst_tmpbuf, RESTORATION_TMPBUF_SIZE);
    }

    EbPictureBufferDescInitData temp_lf_recon_desc_init_data;
    temp_lf_recon_desc_init_data.max_width          = (uint16_t)scs_ptr->max_input_luma_width;
    temp_lf_recon_desc_init_data.max_height         = (uint16_t)scs_ptr->max_input_luma_height;
    temp_lf_recon_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

    temp_lf_recon_desc_init_data.left_padding  = PAD_VALUE;
    temp_lf_recon_desc_init_data.right_padding = PAD_VALUE;
    temp_lf_recon_desc_init_data.top_padding   = PAD_VALUE;
    temp_lf_recon_desc_init_data.bot_padding   = PAD_VALUE;
    temp_lf_recon_desc_init_data.split_mode    = EB_FALSE;
    temp_lf_recon_desc_init_data.color_format  = color_format;

    if (config->is_16bit_pipeline || is_16bit) {
        temp_lf_recon_desc_init_data.bit_depth = EB_16BIT;
        EB_NEW(context_ptr->temp_lf_recon_picture16bit_ptr,
               eb_recon_picture_buffer_desc_ctor,
               (EbPtr)&temp_lf_recon_desc_init_data);
    } else {
        temp_lf_recon_desc_init_data.bit_depth = EB_8BIT;
        EB_NEW(context_ptr->temp_lf_recon_picture_ptr,
               eb_recon_picture_buffer_desc_ctor,
               (EbPtr)&temp_lf_recon_desc_init_data);
    }

    return EB_ErrorNone;
}
void get_own_recon(SequenceControlSet *scs_ptr, PictureControlSet *pcs_ptr,
                   RestContext *context_ptr, EbBool is_16bit) {

    const uint32_t ss_x = scs_ptr->subsampling_x;
    const uint32_t ss_y = scs_ptr->subsampling_y;

    EbPictureBufferDesc *recon_picture_ptr;
    if (is_16bit) {
        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_picture_ptr =
                ((EbReferenceObject *)
                     pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                    ->reference_picture16bit;
        else
            recon_picture_ptr = pcs_ptr->recon_picture16bit_ptr;

        uint16_t *rec_ptr = (uint16_t *)recon_picture_ptr->buffer_y + recon_picture_ptr->origin_x +
                            recon_picture_ptr->origin_y * recon_picture_ptr->stride_y;
        uint16_t *rec_ptr_cb = (uint16_t *)recon_picture_ptr->buffer_cb +
                               recon_picture_ptr->origin_x / 2 +
                               recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb;
        uint16_t *rec_ptr_cr = (uint16_t *)recon_picture_ptr->buffer_cr +
                               recon_picture_ptr->origin_x / 2 +
                               recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr;

        EbPictureBufferDesc *org_rec = context_ptr->org_rec_frame;
        uint16_t *           org_ptr = (uint16_t *)org_rec->buffer_y + org_rec->origin_x +
                            org_rec->origin_y * org_rec->stride_y;
        uint16_t *org_ptr_cb = (uint16_t *)org_rec->buffer_cb + org_rec->origin_x / 2 +
                               org_rec->origin_y / 2 * org_rec->stride_cb;
        uint16_t *org_ptr_cr = (uint16_t *)org_rec->buffer_cr + org_rec->origin_x / 2 +
                               org_rec->origin_y / 2 * org_rec->stride_cr;

        for (int r = 0; r < recon_picture_ptr->height; ++r)
            memcpy(org_ptr + r * org_rec->stride_y,
                   rec_ptr + r * recon_picture_ptr->stride_y,
                   recon_picture_ptr->width << 1);

        for (int r = 0; r < (recon_picture_ptr->height >> ss_y); ++r) {
            memcpy(org_ptr_cb + r * org_rec->stride_cb,
                   rec_ptr_cb + r * recon_picture_ptr->stride_cb,
                   (recon_picture_ptr->width >> ss_x) << 1);
            memcpy(org_ptr_cr + r * org_rec->stride_cr,
                   rec_ptr_cr + r * recon_picture_ptr->stride_cr,
                   (recon_picture_ptr->width >> ss_x) << 1);
        }
    } else {
        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_picture_ptr =
                ((EbReferenceObject *)
                     pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                    ->reference_picture;
        else
            recon_picture_ptr = pcs_ptr->recon_picture_ptr;

        uint8_t *rec_ptr =
            &((recon_picture_ptr
                   ->buffer_y)[recon_picture_ptr->origin_x +
                               recon_picture_ptr->origin_y * recon_picture_ptr->stride_y]);
        uint8_t *rec_ptr_cb =
            &((recon_picture_ptr
                   ->buffer_cb)[recon_picture_ptr->origin_x / 2 +
                                recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb]);
        uint8_t *rec_ptr_cr =
            &((recon_picture_ptr
                   ->buffer_cr)[recon_picture_ptr->origin_x / 2 +
                                recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr]);

        EbPictureBufferDesc *org_rec = context_ptr->org_rec_frame;
        uint8_t *            org_ptr =
            &((org_rec->buffer_y)[org_rec->origin_x + org_rec->origin_y * org_rec->stride_y]);
        uint8_t *org_ptr_cb = &((org_rec->buffer_cb)[org_rec->origin_x / 2 +
                                                     org_rec->origin_y / 2 * org_rec->stride_cb]);
        uint8_t *org_ptr_cr = &((org_rec->buffer_cr)[org_rec->origin_x / 2 +
                                                     org_rec->origin_y / 2 * org_rec->stride_cr]);

        for (int r = 0; r < recon_picture_ptr->height; ++r)
            memcpy(org_ptr + r * org_rec->stride_y,
                   rec_ptr + r * recon_picture_ptr->stride_y,
                   recon_picture_ptr->width);

        for (int r = 0; r < (recon_picture_ptr->height >> ss_y); ++r) {
            memcpy(org_ptr_cb + r * org_rec->stride_cb,
                   rec_ptr_cb + r * recon_picture_ptr->stride_cb,
                   (recon_picture_ptr->width >> ss_x));
            memcpy(org_ptr_cr + r * org_rec->stride_cr,
                   rec_ptr_cr + r * recon_picture_ptr->stride_cr,
                   (recon_picture_ptr->width >> ss_x));
        }
    }
}

void set_unscaled_input_16bit(PictureControlSet *pcs_ptr){
    uint16_t *unscaled_input_frame16bit[MAX_MB_PLANE];
    unscaled_input_frame16bit[0] = (uint16_t *)pcs_ptr->input_frame16bit->buffer_y;
    unscaled_input_frame16bit[1] = (uint16_t *)pcs_ptr->input_frame16bit->buffer_cb;
    unscaled_input_frame16bit[2] = (uint16_t *)pcs_ptr->input_frame16bit->buffer_cr;

    pack_highbd_pic(pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr, unscaled_input_frame16bit, 1, 1, EB_TRUE);
    copy_buffer_info(pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr, pcs_ptr->input_frame16bit);
}

void derive_blk_pointers_enc(EbPictureBufferDesc *recon_picture_buf, int32_t plane,
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

EbErrorType copy_recon_enc(SequenceControlSet *scs_ptr,
                           EbPictureBufferDesc*recon_picture_src,
                           EbPictureBufferDesc *recon_picture_dst,
                           int num_planes,
                           int skip_copy){

    recon_picture_dst->origin_x      = recon_picture_src->origin_x;
    recon_picture_dst->origin_y      = recon_picture_src->origin_y;
    recon_picture_dst->width         = recon_picture_src->width;
    recon_picture_dst->height        = recon_picture_src->height;
    recon_picture_dst->max_width     = recon_picture_src->max_width;
    recon_picture_dst->max_height    = recon_picture_src->max_height;
    recon_picture_dst->bit_depth     = recon_picture_src->bit_depth;
    recon_picture_dst->color_format  = recon_picture_src->color_format;

    recon_picture_dst->stride_y  = recon_picture_src->stride_y;
    recon_picture_dst->stride_cb = recon_picture_src->stride_cb;
    recon_picture_dst->stride_cr = recon_picture_src->stride_cr;

    recon_picture_dst->luma_size    = recon_picture_src->luma_size;
    recon_picture_dst->chroma_size  = recon_picture_src->chroma_size;
    recon_picture_dst->packed_flag   = recon_picture_src->packed_flag;

    recon_picture_dst->stride_bit_inc_y = recon_picture_src->stride_bit_inc_y;
    recon_picture_dst->stride_bit_inc_cb = recon_picture_src->stride_bit_inc_cb;
    recon_picture_dst->stride_bit_inc_cr = recon_picture_src->stride_bit_inc_cr;

    recon_picture_dst->buffer_enable_mask = scs_ptr->seq_header.color_config.mono_chrome ?
                                            PICTURE_BUFFER_DESC_LUMA_MASK : PICTURE_BUFFER_DESC_FULL_MASK;

    uint32_t bytesPerPixel = (recon_picture_dst->bit_depth == EB_8BIT) ? 1 : 2;

    // Allocate the Picture Buffers (luma & chroma)
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_y, recon_picture_dst->luma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_y, 0,
               recon_picture_dst->luma_size * bytesPerPixel);
    }
    else
        recon_picture_dst->buffer_y = 0;
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_cb, recon_picture_dst->chroma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_cb, 0,
               recon_picture_dst->chroma_size * bytesPerPixel);
    }
    else
        recon_picture_dst->buffer_cb = 0;
    if (recon_picture_dst->buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {
        EB_MALLOC_ALIGNED(recon_picture_dst->buffer_cr, recon_picture_dst->chroma_size * bytesPerPixel);
        memset(recon_picture_dst->buffer_cr, 0,
               recon_picture_dst->chroma_size * bytesPerPixel);
    }
    else
        recon_picture_dst->buffer_cr = 0;

    int use_highbd = (scs_ptr->static_config.encoder_bit_depth > 8);

    if(!skip_copy){
        for (int plane = 0; plane < num_planes; ++plane) {
            uint8_t *src_buf, *dst_buf;
            int32_t src_stride, dst_stride;

            int sub_x = plane ? scs_ptr->subsampling_x : 0;
            int sub_y = plane ? scs_ptr->subsampling_y : 0;

            derive_blk_pointers_enc(recon_picture_src, plane, 0, 0, (void *)&src_buf,
                                    &src_stride, sub_x, sub_y);
            derive_blk_pointers_enc(recon_picture_dst, plane, 0, 0, (void *)&dst_buf,
                                    &dst_stride, sub_x, sub_y);

            int height = (recon_picture_src->height >> sub_y);
            for (int row = 0; row < height; ++row) {
                memcpy(dst_buf, src_buf, (recon_picture_src->width >> sub_x) *
                                         sizeof(*src_buf) << use_highbd);
                src_buf += src_stride << use_highbd;
                dst_buf += dst_stride << use_highbd;
            }
        }
    }

    return EB_ErrorNone;
}

void get_recon_pic(PictureControlSet *pcs_ptr,
                   EbPictureBufferDesc **recon_ptr,
                   EbBool is_highbd){
    if(!is_highbd){
        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            *recon_ptr = ((EbReferenceObject*)pcs_ptr->parent_pcs_ptr->
                    reference_picture_wrapper_ptr->object_ptr)->reference_picture;
        else
            *recon_ptr = pcs_ptr->recon_picture_ptr;
    }else {
        if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            *recon_ptr = ((EbReferenceObject *) pcs_ptr->parent_pcs_ptr->
                    reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            *recon_ptr = pcs_ptr->recon_picture16bit_ptr;
    }
}

void eb_av1_superres_upscale_frame(struct Av1Common *cm,
                                   PictureControlSet *pcs_ptr,
                                   SequenceControlSet *scs_ptr)
{
    // Set these parameters for testing since they are not correctly populated yet
    EbPictureBufferDesc *recon_ptr;

    EbBool is_16bit = (EbBool)(scs_ptr->static_config.encoder_bit_depth > EB_8BIT) || (scs_ptr->static_config.is_16bit_pipeline);

    get_recon_pic(pcs_ptr,
                  &recon_ptr,
                  is_16bit);

    uint16_t ss_x = scs_ptr->subsampling_x;
    uint16_t ss_y = scs_ptr->subsampling_y;
    const int num_planes = scs_ptr->seq_header.color_config.mono_chrome ? 1 : MAX_MB_PLANE;

    EbPictureBufferDesc recon_pic_temp;
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
        int32_t src_stride, dst_stride;

        int sub_x = plane ? ss_x : 0;
        int sub_y = plane ? ss_y : 0;
        derive_blk_pointers_enc(src, plane, 0, 0, (void *) &src_buf, &src_stride,
                                sub_x, sub_y);
        derive_blk_pointers_enc(dst, plane, 0, 0, (void *) &dst_buf, &dst_stride,
                                sub_x, sub_y);

        av1_upscale_normative_rows(cm, (const uint8_t *) src_buf, src_stride, dst_buf,
                                   dst_stride, src->height >> sub_x,
                                   sub_x, bit_depth, is_16bit);
    }

    // free the memory
    EB_FREE_ALIGNED_ARRAY(ps_recon_pic_temp->buffer_y);
    EB_FREE_ALIGNED_ARRAY(ps_recon_pic_temp->buffer_cb);
    EB_FREE_ALIGNED_ARRAY(ps_recon_pic_temp->buffer_cr);

}

/******************************************************
 * Rest Kernel
 ******************************************************/
void *rest_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext *   thread_context_ptr = (EbThreadContext *)input_ptr;
    RestContext *       context_ptr        = (RestContext *)thread_context_ptr->priv;
    PictureControlSet * pcs_ptr;
    SequenceControlSet *scs_ptr;
    FrameHeader *       frm_hdr;

    //// Input
    EbObjectWrapper *cdef_results_wrapper_ptr;
    CdefResults *    cdef_results_ptr;

    //// Output
    EbObjectWrapper *    rest_results_wrapper_ptr;
    RestResults *        rest_results_ptr;
    EbObjectWrapper *    picture_demux_results_wrapper_ptr;
    PictureDemuxResults *picture_demux_results_rtr;
    // SB Loop variables

    for (;;) {
        // Get Cdef Results
        EB_GET_FULL_OBJECT(context_ptr->rest_input_fifo_ptr, &cdef_results_wrapper_ptr);

        cdef_results_ptr = (CdefResults *)cdef_results_wrapper_ptr->object_ptr;
        pcs_ptr          = (PictureControlSet *)cdef_results_ptr->pcs_wrapper_ptr->object_ptr;
        scs_ptr          = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
        frm_hdr          = &pcs_ptr->parent_pcs_ptr->frm_hdr;
        EbBool     is_16bit = (EbBool)(scs_ptr->static_config.encoder_bit_depth > EB_8BIT);
        Av1Common *cm       = pcs_ptr->parent_pcs_ptr->av1_cm;

        if (scs_ptr->seq_header.enable_restoration && frm_hdr->allow_intrabc == 0) {

            // ------- start: Normative upscaling - super-resolution tool
            if(!av1_superres_unscaled(&cm->frm_size)) {
                eb_av1_superres_upscale_frame(cm,
                                              pcs_ptr,
                                              scs_ptr);

                if(scs_ptr->static_config.is_16bit_pipeline || is_16bit){
                    set_unscaled_input_16bit(pcs_ptr);
                }
            }
            // ------- end: Normative upscaling - super-resolution tool
            get_own_recon(scs_ptr, pcs_ptr, context_ptr,
                scs_ptr->static_config.is_16bit_pipeline || is_16bit);
            Yv12BufferConfig cpi_source;
            pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr->is_16bit_pipeline = scs_ptr->static_config.is_16bit_pipeline;
            link_eb_to_aom_buffer_desc(scs_ptr->static_config.is_16bit_pipeline || is_16bit
                                       ? pcs_ptr->input_frame16bit
                                       : pcs_ptr->parent_pcs_ptr->enhanced_unscaled_picture_ptr,
                                       &cpi_source,
                                       scs_ptr->max_input_pad_right,
                                       scs_ptr->max_input_pad_bottom,
                                       scs_ptr->static_config.is_16bit_pipeline || is_16bit);

            Yv12BufferConfig trial_frame_rst;
            link_eb_to_aom_buffer_desc(context_ptr->trial_frame_rst, &trial_frame_rst,
                                       scs_ptr->max_input_pad_right,
                                       scs_ptr->max_input_pad_bottom,
                                       scs_ptr->static_config.is_16bit_pipeline || is_16bit);

            Yv12BufferConfig org_fts;
            link_eb_to_aom_buffer_desc(context_ptr->org_rec_frame, &org_fts,
                                       scs_ptr->max_input_pad_right,
                                       scs_ptr->max_input_pad_bottom,
                                       scs_ptr->static_config.is_16bit_pipeline || is_16bit);

            restoration_seg_search(context_ptr->rst_tmpbuf,
                                   &org_fts,
                                   &cpi_source,
                                   &trial_frame_rst,
                                   pcs_ptr,
                                   cdef_results_ptr->segment_index);
        }

        //all seg based search is done. update total processed segments. if all done, finish the search and perfrom application.
        eb_block_on_mutex(pcs_ptr->rest_search_mutex);

        pcs_ptr->tot_seg_searched_rest++;
        if (pcs_ptr->tot_seg_searched_rest == pcs_ptr->rest_segments_total_count) {
            if (scs_ptr->seq_header.enable_restoration && frm_hdr->allow_intrabc == 0) {
                rest_finish_search(pcs_ptr->parent_pcs_ptr, pcs_ptr->parent_pcs_ptr->av1x, pcs_ptr->parent_pcs_ptr->av1_cm);

                if (cm->rst_info[0].frame_restoration_type != RESTORE_NONE ||
                    cm->rst_info[1].frame_restoration_type != RESTORE_NONE ||
                    cm->rst_info[2].frame_restoration_type != RESTORE_NONE) {
                    eb_av1_loop_restoration_filter_frame(cm->frame_to_show, cm, 0);
                }
            } else {
                cm->rst_info[0].frame_restoration_type = RESTORE_NONE;
                cm->rst_info[1].frame_restoration_type = RESTORE_NONE;
                cm->rst_info[2].frame_restoration_type = RESTORE_NONE;
            }

            uint8_t best_ep_cnt = 0;
            uint8_t best_ep     = 0;
            for (uint8_t i = 0; i < SGRPROJ_PARAMS; i++) {
                if (cm->sg_frame_ep_cnt[i] > best_ep_cnt) {
                    best_ep     = i;
                    best_ep_cnt = cm->sg_frame_ep_cnt[i];
                }
            }
            cm->sg_frame_ep = best_ep;

            if (pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL) {
                // copy stat to ref object (intra_coded_area, Luminance, Scene change detection flags)
                copy_statistics_to_ref_obj_ect(pcs_ptr, scs_ptr);
            }

            // PSNR and SSIM Calculation.
            // Note: if temporal_filtering is used, memory needs to be freed in the last of these calls
            if (scs_ptr->static_config.stat_report) {
                psnr_calculations(pcs_ptr, scs_ptr, EB_FALSE);
                ssim_calculations(pcs_ptr, scs_ptr, EB_TRUE /* free memory here */);
            }

            // Pad the reference picture and set ref POC
            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                pad_ref_and_set_flags(pcs_ptr, scs_ptr);
            if (scs_ptr->static_config.recon_enabled) { recon_output(pcs_ptr, scs_ptr); }

            if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag) {
                // Get Empty PicMgr Results
                eb_get_empty_object(context_ptr->picture_demux_fifo_ptr,
                                    &picture_demux_results_wrapper_ptr);

                picture_demux_results_rtr =
                    (PictureDemuxResults *)picture_demux_results_wrapper_ptr->object_ptr;
                picture_demux_results_rtr->reference_picture_wrapper_ptr =
                    pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr;
                picture_demux_results_rtr->scs_wrapper_ptr = pcs_ptr->scs_wrapper_ptr;
                picture_demux_results_rtr->picture_number  = pcs_ptr->picture_number;
                picture_demux_results_rtr->picture_type    = EB_PIC_REFERENCE;

                // Post Reference Picture
                eb_post_full_object(picture_demux_results_wrapper_ptr);
            }
            //Jing: TODO
            //Consider to add parallelism here, sending line by line, not waiting for a full frame
            int sb_size_log2 = scs_ptr->seq_header.sb_size_log2;
            for (int tile_row_idx = 0;
                 tile_row_idx < pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_rows;
                 tile_row_idx++) {
                uint16_t tile_height_in_sb =
                    (cm->tiles_info.tile_row_start_mi[tile_row_idx + 1] -
                     cm->tiles_info.tile_row_start_mi[tile_row_idx] + (1 << sb_size_log2) - 1)
                     >> sb_size_log2;
                for (int tile_col_idx = 0;
                     tile_col_idx < pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_cols;
                     tile_col_idx++) {
                    const int tile_idx =
                        tile_row_idx * pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_cols +
                        tile_col_idx;
                    eb_get_empty_object(context_ptr->rest_output_fifo_ptr,
                                        &rest_results_wrapper_ptr);
                    rest_results_ptr = (struct RestResults *)rest_results_wrapper_ptr->object_ptr;
                    rest_results_ptr->pcs_wrapper_ptr = cdef_results_ptr->pcs_wrapper_ptr;
                    rest_results_ptr->completed_sb_row_index_start = 0;
                    // Set to tile rows
                    rest_results_ptr->completed_sb_row_count = tile_height_in_sb;
                    rest_results_ptr->tile_index             = tile_idx;
                    // Post Rest Results
                    eb_post_full_object(rest_results_wrapper_ptr);
                }
            }
        }
        eb_release_mutex(pcs_ptr->rest_search_mutex);

        // Release input Results
        eb_release_object(cdef_results_wrapper_ptr);
    }

    return NULL;
}
