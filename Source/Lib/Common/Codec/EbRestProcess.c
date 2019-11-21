// clang-format off
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
#include "EbPsnr.h"
#include "EbReferenceObject.h"
#include "EbPictureControlSet.h"

/**************************************
 * Rest Context
 **************************************/
typedef struct RestContext
{
    EbDctor                       dctor;
    EbFifo                       *rest_input_fifo_ptr;
    EbFifo                       *rest_output_fifo_ptr;
    EbFifo                       *picture_demux_fifo_ptr;

    EbPictureBufferDesc          *trial_frame_rst;

    EbPictureBufferDesc          *temp_lf_recon_picture_ptr;
    EbPictureBufferDesc          *temp_lf_recon_picture16bit_ptr;

    EbPictureBufferDesc           *org_rec_frame; // while doing the filtering recon gets updated uisng setup/restore processing_stripe_bounadaries
                                                    // many threads doing the above will result in race condition.
                                                    // each thread will hence have his own copy of recon to work on.
                                                    // later we can have a search version that does not need the exact right recon
    int32_t *rst_tmpbuf;
} RestContext;

void recon_output(
    PictureControlSet    *picture_control_set_ptr,
    SequenceControlSet   *sequence_control_set_ptr);
void eb_av1_loop_restoration_filter_frame(Yv12BufferConfig *frame,
    Av1Common *cm, int32_t optimized_lr);
void copy_statistics_to_ref_object(
    PictureControlSet    *picture_control_set_ptr,
    SequenceControlSet   *sequence_control_set_ptr);
void psnr_calculations(
    PictureControlSet    *picture_control_set_ptr,
    SequenceControlSet   *sequence_control_set_ptr);
void pad_ref_and_set_flags(
    PictureControlSet    *picture_control_set_ptr,
    SequenceControlSet   *sequence_control_set_ptr);
void generate_padding(
    EbByte              src_pic,
    uint32_t            src_stride,
    uint32_t            original_src_width,
    uint32_t            original_src_height,
    uint32_t            padding_width,
    uint32_t            padding_height);
void restoration_seg_search(
    int32_t                 *rst_tmpbuf,
    Yv12BufferConfig       *org_fts,
    const Yv12BufferConfig *src,
    Yv12BufferConfig       *trial_frame_rst,
    PictureControlSet    *pcs_ptr,
    uint32_t                segment_index);
void rest_finish_search(Macroblock *x, Av1Common *const cm);

static void rest_context_dctor(EbPtr p)
{
    EbThreadContext* thread_context_ptr = (EbThreadContext*)p;
    RestContext *obj = (RestContext*)thread_context_ptr->priv;
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
EbErrorType rest_context_ctor(
    EbThreadContext     *thread_context_ptr,
    const EbEncHandle   *enc_handle_ptr,
    int                 index,
    int                 demux_index)
{
    const SequenceControlSet* sequence_control_set_ptr = enc_handle_ptr->sequence_control_set_instance_array[0]->sequence_control_set_ptr;
    const EbSvtAv1EncConfiguration* config = &sequence_control_set_ptr->static_config;
    EbBool is16bit = (EbBool)(config->encoder_bit_depth > EB_8BIT);
    EbColorFormat color_format = config->encoder_color_format;

    RestContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv = context_ptr;
    thread_context_ptr->dctor = rest_context_dctor;

    // Input/Output System Resource Manager FIFOs
    context_ptr->rest_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(enc_handle_ptr->cdef_results_resource_ptr, index);
    context_ptr->rest_output_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->rest_results_resource_ptr, index);
    context_ptr->picture_demux_fifo_ptr =
        eb_system_resource_get_producer_fifo(enc_handle_ptr->picture_demux_results_resource_ptr, demux_index);

    {
        EbPictureBufferDescInitData initData;

        initData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
        initData.max_width = (uint16_t)sequence_control_set_ptr->max_input_luma_width;
        initData.max_height = (uint16_t)sequence_control_set_ptr->max_input_luma_height;
        initData.bit_depth = is16bit ? EB_16BIT : EB_8BIT;
        initData.color_format = color_format;
        initData.left_padding = AOM_BORDER_IN_PIXELS;
        initData.right_padding = AOM_BORDER_IN_PIXELS;
        initData.top_padding = AOM_BORDER_IN_PIXELS;
        initData.bot_padding = AOM_BORDER_IN_PIXELS;
        initData.split_mode = EB_FALSE;

        EB_NEW(
            context_ptr->trial_frame_rst,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&initData);

        EB_NEW(
            context_ptr->org_rec_frame,
            eb_picture_buffer_desc_ctor,
            (EbPtr)&initData);

         EB_MALLOC_ALIGNED(context_ptr->rst_tmpbuf, RESTORATION_TMPBUF_SIZE);
    }

    EbPictureBufferDescInitData tempLfReconDescInitData;
    tempLfReconDescInitData.max_width = (uint16_t)sequence_control_set_ptr->max_input_luma_width;
    tempLfReconDescInitData.max_height = (uint16_t)sequence_control_set_ptr->max_input_luma_height;
    tempLfReconDescInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

    tempLfReconDescInitData.left_padding = PAD_VALUE;
    tempLfReconDescInitData.right_padding = PAD_VALUE;
    tempLfReconDescInitData.top_padding = PAD_VALUE;
    tempLfReconDescInitData.bot_padding = PAD_VALUE;
    tempLfReconDescInitData.split_mode = EB_FALSE;
    tempLfReconDescInitData.color_format = color_format;

    if (is16bit) {
        tempLfReconDescInitData.bit_depth = EB_16BIT;
        EB_NEW(
            context_ptr->temp_lf_recon_picture16bit_ptr,
            eb_recon_picture_buffer_desc_ctor,
            (EbPtr)&tempLfReconDescInitData);
    }
    else {
        tempLfReconDescInitData.bit_depth = EB_8BIT;
        EB_NEW(
            context_ptr->temp_lf_recon_picture_ptr,
            eb_recon_picture_buffer_desc_ctor,
            (EbPtr)&tempLfReconDescInitData);
    }

    return EB_ErrorNone;
}
void   get_own_recon(
    SequenceControlSet                    *sequence_control_set_ptr,
    PictureControlSet                     *picture_control_set_ptr,
    RestContext                            *context_ptr,
    EbBool  is16bit)
{
    EbPictureBufferDesc  * recon_picture_ptr;
    if (is16bit) {
        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_picture_ptr = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            recon_picture_ptr = picture_control_set_ptr->recon_picture16bit_ptr;

        uint16_t*  rec_ptr = (uint16_t*)recon_picture_ptr->buffer_y + recon_picture_ptr->origin_x + recon_picture_ptr->origin_y     * recon_picture_ptr->stride_y;
        uint16_t*  rec_ptr_cb = (uint16_t*)recon_picture_ptr->buffer_cb + recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb;
        uint16_t*  rec_ptr_cr = (uint16_t*)recon_picture_ptr->buffer_cr + recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr;

        EbPictureBufferDesc *org_rec = context_ptr->org_rec_frame;
        uint16_t*  org_ptr = (uint16_t*)org_rec->buffer_y + org_rec->origin_x + org_rec->origin_y     * org_rec->stride_y;
        uint16_t*  org_ptr_cb = (uint16_t*)org_rec->buffer_cb + org_rec->origin_x / 2 + org_rec->origin_y / 2 * org_rec->stride_cb;
        uint16_t*  org_ptr_cr = (uint16_t*)org_rec->buffer_cr + org_rec->origin_x / 2 + org_rec->origin_y / 2 * org_rec->stride_cr;

        for (int r = 0; r < sequence_control_set_ptr->seq_header.max_frame_height; ++r)
            memcpy(org_ptr + r * org_rec->stride_y, rec_ptr + r * recon_picture_ptr->stride_y, sequence_control_set_ptr->seq_header.max_frame_width << 1);

        for (int r = 0; r < sequence_control_set_ptr->seq_header.max_frame_height / 2; ++r) {
            memcpy(org_ptr_cb + r * org_rec->stride_cb, rec_ptr_cb + r * recon_picture_ptr->stride_cb, (sequence_control_set_ptr->seq_header.max_frame_width / 2) << 1);
            memcpy(org_ptr_cr + r * org_rec->stride_cr, rec_ptr_cr + r * recon_picture_ptr->stride_cr, (sequence_control_set_ptr->seq_header.max_frame_width / 2) << 1);
        }
    }
    else {
        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_picture_ptr = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
        else
            recon_picture_ptr = picture_control_set_ptr->recon_picture_ptr;

        uint8_t * rec_ptr = &((recon_picture_ptr->buffer_y)[recon_picture_ptr->origin_x + recon_picture_ptr->origin_y * recon_picture_ptr->stride_y]);
        uint8_t *  rec_ptr_cb = &((recon_picture_ptr->buffer_cb)[recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb]);
        uint8_t *  rec_ptr_cr = &((recon_picture_ptr->buffer_cr)[recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr]);

        EbPictureBufferDesc *org_rec = context_ptr->org_rec_frame;
        uint8_t *  org_ptr = &((org_rec->buffer_y)[org_rec->origin_x + org_rec->origin_y * org_rec->stride_y]);
        uint8_t *  org_ptr_cb = &((org_rec->buffer_cb)[org_rec->origin_x / 2 + org_rec->origin_y / 2 * org_rec->stride_cb]);
        uint8_t *  org_ptr_cr = &((org_rec->buffer_cr)[org_rec->origin_x / 2 + org_rec->origin_y / 2 * org_rec->stride_cr]);

        for (int r = 0; r < sequence_control_set_ptr->seq_header.max_frame_height; ++r)
            memcpy(org_ptr + r * org_rec->stride_y, rec_ptr + r * recon_picture_ptr->stride_y, sequence_control_set_ptr->seq_header.max_frame_width);

        for (int r = 0; r < sequence_control_set_ptr->seq_header.max_frame_height / 2; ++r) {
            memcpy(org_ptr_cb + r * org_rec->stride_cb, rec_ptr_cb + r * recon_picture_ptr->stride_cb, (sequence_control_set_ptr->seq_header.max_frame_width / 2));
            memcpy(org_ptr_cr + r * org_rec->stride_cr, rec_ptr_cr + r * recon_picture_ptr->stride_cr, (sequence_control_set_ptr->seq_header.max_frame_width / 2));
        }
    }
}

/******************************************************
 * Rest Kernel
 ******************************************************/
void* rest_kernel(void *input_ptr)
{
    // Context & SCS & PCS
    EbThreadContext                       *thread_context_ptr = (EbThreadContext*)input_ptr;
    RestContext                           *context_ptr = (RestContext*)thread_context_ptr->priv;
    PictureControlSet                     *picture_control_set_ptr;
    SequenceControlSet                    *sequence_control_set_ptr;
    FrameHeader                           *frm_hdr;

    //// Input
    EbObjectWrapper                       *cdef_results_wrapper_ptr;
    CdefResults                         *cdef_results_ptr;

    //// Output
    EbObjectWrapper                       *rest_results_wrapper_ptr;
    RestResults*                          rest_results_ptr;
    EbObjectWrapper                       *picture_demux_results_wrapper_ptr;
    PictureDemuxResults                   *picture_demux_results_rtr;
    // SB Loop variables

    for (;;) {
        // Get Cdef Results
        eb_get_full_object(
            context_ptr->rest_input_fifo_ptr,
            &cdef_results_wrapper_ptr);

        cdef_results_ptr = (CdefResults*)cdef_results_wrapper_ptr->object_ptr;
        picture_control_set_ptr = (PictureControlSet*)cdef_results_ptr->picture_control_set_wrapper_ptr->object_ptr;
        sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
        frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
        uint8_t lcuSizeLog2 = (uint8_t)Log2f(sequence_control_set_ptr->sb_size_pix);
        EbBool  is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
        Av1Common* cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;

        if (sequence_control_set_ptr->seq_header.enable_restoration && frm_hdr->allow_intrabc == 0)
        {
            get_own_recon(sequence_control_set_ptr, picture_control_set_ptr, context_ptr, is16bit);

            Yv12BufferConfig cpi_source;
            link_eb_to_aom_buffer_desc(
                is16bit ? picture_control_set_ptr->input_frame16bit : picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                &cpi_source);

            Yv12BufferConfig trial_frame_rst;
            link_eb_to_aom_buffer_desc(
                context_ptr->trial_frame_rst,
                &trial_frame_rst);

            Yv12BufferConfig org_fts;
            link_eb_to_aom_buffer_desc(
                context_ptr->org_rec_frame,
                &org_fts);

            restoration_seg_search(
                context_ptr->rst_tmpbuf,
                &org_fts,
                &cpi_source,
                &trial_frame_rst,
                picture_control_set_ptr,
                cdef_results_ptr->segment_index);
        }

        //all seg based search is done. update total processed segments. if all done, finish the search and perfrom application.
        eb_block_on_mutex(picture_control_set_ptr->rest_search_mutex);

        picture_control_set_ptr->tot_seg_searched_rest++;
        if (picture_control_set_ptr->tot_seg_searched_rest == picture_control_set_ptr->rest_segments_total_count)
        {
            if (sequence_control_set_ptr->seq_header.enable_restoration && frm_hdr->allow_intrabc == 0) {
                rest_finish_search(
                    picture_control_set_ptr->parent_pcs_ptr->av1x,
                    picture_control_set_ptr->parent_pcs_ptr->av1_cm);

                if (cm->rst_info[0].frame_restoration_type != RESTORE_NONE ||
                    cm->rst_info[1].frame_restoration_type != RESTORE_NONE ||
                    cm->rst_info[2].frame_restoration_type != RESTORE_NONE)
                {
                    eb_av1_loop_restoration_filter_frame(
                        cm->frame_to_show,
                        cm,
                        0);
                }
            }
            else {
                cm->rst_info[0].frame_restoration_type = RESTORE_NONE;
                cm->rst_info[1].frame_restoration_type = RESTORE_NONE;
                cm->rst_info[2].frame_restoration_type = RESTORE_NONE;
            }

            uint8_t best_ep_cnt = 0;
            uint8_t best_ep = 0;
            for (uint8_t i = 0; i < SGRPROJ_PARAMS; i++) {
                if (cm->sg_frame_ep_cnt[i] > best_ep_cnt) {
                    best_ep = i;
                    best_ep_cnt = cm->sg_frame_ep_cnt[i];
                }
            }
            cm->sg_frame_ep = best_ep;

            if (picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL) {
                // copy stat to ref object (intra_coded_area, Luminance, Scene change detection flags)
                copy_statistics_to_ref_object(
                    picture_control_set_ptr,
                    sequence_control_set_ptr);
            }

            // PSNR Calculation
            if (sequence_control_set_ptr->static_config.stat_report)
                psnr_calculations(
                    picture_control_set_ptr,
                    sequence_control_set_ptr);

            // Pad the reference picture and set ref POC
            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                pad_ref_and_set_flags(
                    picture_control_set_ptr,
                    sequence_control_set_ptr);
            if (sequence_control_set_ptr->static_config.recon_enabled) {
                recon_output(
                    picture_control_set_ptr,
                    sequence_control_set_ptr);
            }

            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
            {
                // Get Empty PicMgr Results
                eb_get_empty_object(
                    context_ptr->picture_demux_fifo_ptr,
                    &picture_demux_results_wrapper_ptr);

                picture_demux_results_rtr = (PictureDemuxResults*)picture_demux_results_wrapper_ptr->object_ptr;
                picture_demux_results_rtr->reference_picture_wrapper_ptr = picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr;
                picture_demux_results_rtr->sequence_control_set_wrapper_ptr = picture_control_set_ptr->sequence_control_set_wrapper_ptr;
                picture_demux_results_rtr->picture_number = picture_control_set_ptr->picture_number;
                picture_demux_results_rtr->picture_type = EB_PIC_REFERENCE;

                // Post Reference Picture
                eb_post_full_object(picture_demux_results_wrapper_ptr);
            }

            // Get Empty rest Results to EC
            eb_get_empty_object(
                context_ptr->rest_output_fifo_ptr,
                &rest_results_wrapper_ptr);
            rest_results_ptr = (struct RestResults*)rest_results_wrapper_ptr->object_ptr;
            rest_results_ptr->picture_control_set_wrapper_ptr = cdef_results_ptr->picture_control_set_wrapper_ptr;
            rest_results_ptr->completed_lcu_row_index_start = 0;
            rest_results_ptr->completed_lcu_row_count = ((sequence_control_set_ptr->seq_header.max_frame_height + sequence_control_set_ptr->sb_size_pix - 1) >> lcuSizeLog2);
            // Post Rest Results
            eb_post_full_object(rest_results_wrapper_ptr);
        }
        eb_release_mutex(picture_control_set_ptr->rest_search_mutex);

        // Release input Results
        eb_release_object(cdef_results_wrapper_ptr);
    }

    return EB_NULL;
}
// clang-format on
