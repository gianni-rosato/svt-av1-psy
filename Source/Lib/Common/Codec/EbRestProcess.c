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
#include "EbDefinitions.h"
#if FILT_PROC
#include "EbRestProcess.h"
#include "EbEncDecResults.h"

#include "EbEncDecTasks.h"
#include "EbPictureDemuxResults.h"
#include "EbReferenceObject.h"


void ReconOutput(
    PictureControlSet_t    *picture_control_set_ptr,
    SequenceControlSet_t   *sequence_control_set_ptr);
void av1_loop_restoration_filter_frame(Yv12BufferConfig *frame,
    Av1Common *cm, int32_t optimized_lr);
void CopyStatisticsToRefObject(
    PictureControlSet_t    *picture_control_set_ptr,
    SequenceControlSet_t   *sequence_control_set_ptr);
void PsnrCalculations(
    PictureControlSet_t    *picture_control_set_ptr,
    SequenceControlSet_t   *sequence_control_set_ptr);
void PadRefAndSetFlags(
    PictureControlSet_t    *picture_control_set_ptr,
    SequenceControlSet_t   *sequence_control_set_ptr);
void generate_padding(
    EbByte              src_pic,
    uint32_t            src_stride,
    uint32_t            original_src_width,
    uint32_t            original_src_height,
    uint32_t            padding_width,
    uint32_t            padding_height);
#if REST_M
void restoration_seg_search(
    RestContext_t          *context_ptr,
    Yv12BufferConfig       *org_fts,
    const Yv12BufferConfig *src,
    Yv12BufferConfig       *trial_frame_rst,
    PictureControlSet_t    *pcs_ptr,
    uint32_t                segment_index);
void rest_finish_search(Macroblock *x, Av1Common *const cm);
#endif
/******************************************************
 * Rest Context Constructor
 ******************************************************/
EbErrorType rest_context_ctor(
    RestContext_t **context_dbl_ptr,
    EbFifo_t                *rest_input_fifo_ptr,
    EbFifo_t                *rest_output_fifo_ptr ,
    EbFifo_t                *picture_demux_fifo_ptr,
    EbBool                  is16bit,
    uint32_t                max_input_luma_width,
    uint32_t                max_input_luma_height
   )
{
    EbErrorType return_error = EB_ErrorNone;
    RestContext_t *context_ptr;
    EB_MALLOC(RestContext_t*, context_ptr, sizeof(RestContext_t), EB_N_PTR);
    *context_dbl_ptr = context_ptr;

    // Input/Output System Resource Manager FIFOs
    context_ptr->rest_input_fifo_ptr = rest_input_fifo_ptr;
    context_ptr->rest_output_fifo_ptr = rest_output_fifo_ptr;
    context_ptr->picture_demux_fifo_ptr = picture_demux_fifo_ptr;


    {
        EbPictureBufferDescInitData_t initData;

        initData.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;
        initData.maxWidth = (uint16_t)max_input_luma_width;
        initData.maxHeight = (uint16_t)max_input_luma_height;
        initData.bit_depth = is16bit ? EB_16BIT : EB_8BIT;
        initData.left_padding = AOM_BORDER_IN_PIXELS;
        initData.right_padding = AOM_BORDER_IN_PIXELS;
        initData.top_padding = AOM_BORDER_IN_PIXELS;
        initData.bot_padding = AOM_BORDER_IN_PIXELS;
        initData.splitMode = EB_FALSE;

        return_error = eb_picture_buffer_desc_ctor(
            (EbPtr*)&context_ptr->trial_frame_rst,
            (EbPtr)&initData);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

#if REST_M
         return_error = eb_picture_buffer_desc_ctor(
            (EbPtr*)&context_ptr->org_rec_frame,
                (EbPtr)&initData);

         EB_MALLOC(int32_t *, context_ptr->rst_tmpbuf, RESTORATION_TMPBUF_SIZE, EB_N_PTR);
#endif

    }

    context_ptr->temp_lf_recon_picture16bit_ptr = (EbPictureBufferDesc_t *)EB_NULL;
    context_ptr->temp_lf_recon_picture_ptr = (EbPictureBufferDesc_t *)EB_NULL;
    EbPictureBufferDescInitData_t tempLfReconDescInitData;
    tempLfReconDescInitData.maxWidth = (uint16_t)max_input_luma_width;
    tempLfReconDescInitData.maxHeight = (uint16_t)max_input_luma_height;
    tempLfReconDescInitData.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;

    tempLfReconDescInitData.left_padding = PAD_VALUE;
    tempLfReconDescInitData.right_padding = PAD_VALUE;
    tempLfReconDescInitData.top_padding = PAD_VALUE;
    tempLfReconDescInitData.bot_padding = PAD_VALUE;

    tempLfReconDescInitData.splitMode = EB_FALSE;

    if (is16bit) {
        tempLfReconDescInitData.bit_depth = EB_16BIT;
        return_error = eb_recon_picture_buffer_desc_ctor(
            (EbPtr*)&(context_ptr->temp_lf_recon_picture16bit_ptr),
            (EbPtr)&tempLfReconDescInitData);
    }
    else {
        tempLfReconDescInitData.bit_depth = EB_8BIT;
        return_error = eb_recon_picture_buffer_desc_ctor(
            (EbPtr*)&(context_ptr->temp_lf_recon_picture_ptr),
            (EbPtr)&tempLfReconDescInitData);
    }


    return EB_ErrorNone;
}
#if REST_M
void   get_own_recon(
    SequenceControlSet_t                    *sequence_control_set_ptr,
    PictureControlSet_t                     *picture_control_set_ptr,
    RestContext_t                            *context_ptr,
    EbBool  is16bit)
{
    EbPictureBufferDesc_t  * recon_picture_ptr;
    if (is16bit) {
        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_picture_ptr = ((EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->referencePicture16bit;
        else
            recon_picture_ptr = picture_control_set_ptr->recon_picture16bit_ptr;

        uint16_t*  rec_ptr = (uint16_t*)recon_picture_ptr->buffer_y + recon_picture_ptr->origin_x + recon_picture_ptr->origin_y     * recon_picture_ptr->stride_y;
        uint16_t*  rec_ptr_cb = (uint16_t*)recon_picture_ptr->bufferCb + recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->strideCb;
        uint16_t*  rec_ptr_cr = (uint16_t*)recon_picture_ptr->bufferCr + recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->strideCr;

        EbPictureBufferDesc_t *org_rec = context_ptr->org_rec_frame;
        uint16_t*  org_ptr = (uint16_t*)org_rec->buffer_y + org_rec->origin_x + org_rec->origin_y     * org_rec->stride_y;
        uint16_t*  org_ptr_cb = (uint16_t*)org_rec->bufferCb + org_rec->origin_x / 2 + org_rec->origin_y / 2 * org_rec->strideCb;
        uint16_t*  org_ptr_cr = (uint16_t*)org_rec->bufferCr + org_rec->origin_x / 2 + org_rec->origin_y / 2 * org_rec->strideCr;

        for (int r = 0; r < sequence_control_set_ptr->luma_height; ++r) {
            memcpy(org_ptr + r * org_rec->stride_y, rec_ptr + r * recon_picture_ptr->stride_y, sequence_control_set_ptr->luma_width << 1);
        }

        for (int r = 0; r < sequence_control_set_ptr->luma_height / 2; ++r) {
            memcpy(org_ptr_cb + r * org_rec->strideCb, rec_ptr_cb + r * recon_picture_ptr->strideCb, (sequence_control_set_ptr->luma_width / 2) << 1);
            memcpy(org_ptr_cr + r * org_rec->strideCr, rec_ptr_cr + r * recon_picture_ptr->strideCr, (sequence_control_set_ptr->luma_width / 2) << 1);
        }
    }
    else {
        if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
            recon_picture_ptr = ((EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->referencePicture;
        else
            recon_picture_ptr = picture_control_set_ptr->recon_picture_ptr;

        uint8_t * rec_ptr = &((recon_picture_ptr->buffer_y)[recon_picture_ptr->origin_x + recon_picture_ptr->origin_y * recon_picture_ptr->stride_y]);
        uint8_t *  rec_ptr_cb = &((recon_picture_ptr->bufferCb)[recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->strideCb]);
        uint8_t *  rec_ptr_cr = &((recon_picture_ptr->bufferCr)[recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->strideCr]);

        EbPictureBufferDesc_t *org_rec = context_ptr->org_rec_frame;
        uint8_t *  org_ptr = &((org_rec->buffer_y)[org_rec->origin_x + org_rec->origin_y * org_rec->stride_y]);
        uint8_t *  org_ptr_cb = &((org_rec->bufferCb)[org_rec->origin_x / 2 + org_rec->origin_y / 2 * org_rec->strideCb]);
        uint8_t *  org_ptr_cr = &((org_rec->bufferCr)[org_rec->origin_x / 2 + org_rec->origin_y / 2 * org_rec->strideCr]);

        for (int r = 0; r < sequence_control_set_ptr->luma_height; ++r) {
            memcpy(org_ptr + r * org_rec->stride_y, rec_ptr + r * recon_picture_ptr->stride_y, sequence_control_set_ptr->luma_width);
        }

        for (int r = 0; r < sequence_control_set_ptr->luma_height / 2; ++r) {
            memcpy(org_ptr_cb + r * org_rec->strideCb, rec_ptr_cb + r * recon_picture_ptr->strideCb, (sequence_control_set_ptr->luma_width / 2));
            memcpy(org_ptr_cr + r * org_rec->strideCr, rec_ptr_cr + r * recon_picture_ptr->strideCr, (sequence_control_set_ptr->luma_width / 2));
        }
    }
}
#endif

/******************************************************
 * Rest Kernel
 ******************************************************/
void* rest_kernel(void *input_ptr)
{
    // Context & SCS & PCS
    RestContext_t                            *context_ptr = (RestContext_t*)input_ptr;
    PictureControlSet_t                     *picture_control_set_ptr;
    SequenceControlSet_t                    *sequence_control_set_ptr;

    //// Input
    EbObjectWrapper_t                       *cdef_results_wrapper_ptr;
    CdefResults_t                         *cdef_results_ptr;

    //// Output
    EbObjectWrapper_t                       *rest_results_wrapper_ptr;
    RestResults_t*                          rest_results_ptr;
    EbObjectWrapper_t                       *picture_demux_results_wrapper_ptr;
    PictureDemuxResults_t                   *picture_demux_results_rtr;
    // SB Loop variables


    for (;;) {

        // Get Cdef Results
        eb_get_full_object(
            context_ptr->rest_input_fifo_ptr,
            &cdef_results_wrapper_ptr);

        cdef_results_ptr = (CdefResults_t*)cdef_results_wrapper_ptr->object_ptr;
        picture_control_set_ptr = (PictureControlSet_t*)cdef_results_ptr->picture_control_set_wrapper_ptr->object_ptr;
        sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
        uint8_t lcuSizeLog2 = (uint8_t)Log2f(sequence_control_set_ptr->sb_size_pix);
        EbBool  is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
        Av1Common* cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;

#if  REST_M

        if (sequence_control_set_ptr->enable_restoration)
        {
            get_own_recon(sequence_control_set_ptr, picture_control_set_ptr, context_ptr, is16bit);

            Yv12BufferConfig cpi_source;
            LinkEbToAomBufferDesc(
                is16bit ? picture_control_set_ptr->input_frame16bit : picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                &cpi_source);

            Yv12BufferConfig trial_frame_rst;
            LinkEbToAomBufferDesc(
                context_ptr->trial_frame_rst,
                &trial_frame_rst);

            Yv12BufferConfig org_fts;
            LinkEbToAomBufferDesc(
                context_ptr->org_rec_frame,
                &org_fts);

            restoration_seg_search(
                context_ptr,
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

#endif



#if REST_REF_ONLY
            if (sequence_control_set_ptr->enable_restoration && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag) {
#else
            if (sequence_control_set_ptr->enable_restoration) {
#endif

#if  !REST_M
                av1_loop_restoration_save_boundary_lines(
                    cm->frame_to_show,
                    cm,
                    1);

                Yv12BufferConfig cpi_source;
                LinkEbToAomBufferDesc(
                    is16bit ? picture_control_set_ptr->input_frame16bit : picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                    &cpi_source);

                Yv12BufferConfig trial_frame_rst;
                LinkEbToAomBufferDesc(
                    context_ptr->trial_frame_rst,
                    &trial_frame_rst);


#endif


#if REST_M
                rest_finish_search(
                    picture_control_set_ptr->parent_pcs_ptr->av1x,
                    picture_control_set_ptr->parent_pcs_ptr->av1_cm);
#else
                av1_pick_filter_restoration(
                    &cpi_source,
                    &trial_frame_rst,
                    picture_control_set_ptr->parent_pcs_ptr->av1x,
                    picture_control_set_ptr->parent_pcs_ptr->av1_cm);
#endif

                if (cm->rst_info[0].frame_restoration_type != RESTORE_NONE ||
                    cm->rst_info[1].frame_restoration_type != RESTORE_NONE ||
                    cm->rst_info[2].frame_restoration_type != RESTORE_NONE)
                {
                    av1_loop_restoration_filter_frame(
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

#if FAST_SG
            uint8_t best_ep_cnt = 0;
            uint8_t best_ep = 0;
            for (uint8_t i = 0; i < SGRPROJ_PARAMS; i++) {
                if (cm->sg_frame_ep_cnt[i] > best_ep_cnt) {
                    best_ep = i;
                    best_ep_cnt = cm->sg_frame_ep_cnt[i];
                }
            }
            cm->sg_frame_ep = best_ep;
#endif

            if (picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr != NULL) {
                // copy stat to ref object (intra_coded_area, Luminance, Scene change detection flags)
                CopyStatisticsToRefObject(
                    picture_control_set_ptr,
                    sequence_control_set_ptr);
            }

            // PSNR Calculation
            if (sequence_control_set_ptr->static_config.stat_report) {
                PsnrCalculations(
                    picture_control_set_ptr,
                    sequence_control_set_ptr);
            }

            // Pad the reference picture and set up TMVP flag and ref POC
            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                PadRefAndSetFlags(
                    picture_control_set_ptr,
                    sequence_control_set_ptr);

            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE && picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr)
            {
                EbPictureBufferDesc_t *input_picture_ptr = (EbPictureBufferDesc_t*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
                const uint32_t  SrclumaOffSet = input_picture_ptr->origin_x + input_picture_ptr->origin_y    *input_picture_ptr->stride_y;
                const uint32_t  SrccbOffset = (input_picture_ptr->origin_x >> 1) + (input_picture_ptr->origin_y >> 1)*input_picture_ptr->strideCb;
                const uint32_t  SrccrOffset = (input_picture_ptr->origin_x >> 1) + (input_picture_ptr->origin_y >> 1)*input_picture_ptr->strideCr;

                EbReferenceObject_t   *referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
                EbPictureBufferDesc_t *refDenPic = referenceObject->refDenSrcPicture;
                const uint32_t           ReflumaOffSet = refDenPic->origin_x + refDenPic->origin_y    *refDenPic->stride_y;
                const uint32_t           RefcbOffset = (refDenPic->origin_x >> 1) + (refDenPic->origin_y >> 1)*refDenPic->strideCb;
                const uint32_t           RefcrOffset = (refDenPic->origin_x >> 1) + (refDenPic->origin_y >> 1)*refDenPic->strideCr;

                uint16_t  verticalIdx;

                for (verticalIdx = 0; verticalIdx < refDenPic->height; ++verticalIdx)
                {
                    EB_MEMCPY(refDenPic->buffer_y + ReflumaOffSet + verticalIdx * refDenPic->stride_y,
                        input_picture_ptr->buffer_y + SrclumaOffSet + verticalIdx * input_picture_ptr->stride_y,
                        input_picture_ptr->width);
                }

                for (verticalIdx = 0; verticalIdx < input_picture_ptr->height / 2; ++verticalIdx)
                {
                    EB_MEMCPY(refDenPic->bufferCb + RefcbOffset + verticalIdx * refDenPic->strideCb,
                        input_picture_ptr->bufferCb + SrccbOffset + verticalIdx * input_picture_ptr->strideCb,
                        input_picture_ptr->width / 2);

                    EB_MEMCPY(refDenPic->bufferCr + RefcrOffset + verticalIdx * refDenPic->strideCr,
                        input_picture_ptr->bufferCr + SrccrOffset + verticalIdx * input_picture_ptr->strideCr,
                        input_picture_ptr->width / 2);
                }

                generate_padding(
                    refDenPic->buffer_y,
                    refDenPic->stride_y,
                    refDenPic->width,
                    refDenPic->height,
                    refDenPic->origin_x,
                    refDenPic->origin_y);

                generate_padding(
                    refDenPic->bufferCb,
                    refDenPic->strideCb,
                    refDenPic->width >> 1,
                    refDenPic->height >> 1,
                    refDenPic->origin_x >> 1,
                    refDenPic->origin_y >> 1);

                generate_padding(
                    refDenPic->bufferCr,
                    refDenPic->strideCr,
                    refDenPic->width >> 1,
                    refDenPic->height >> 1,
                    refDenPic->origin_x >> 1,
                    refDenPic->origin_y >> 1);
            }
            if (sequence_control_set_ptr->static_config.recon_enabled) {
                ReconOutput(
                    picture_control_set_ptr,
                    sequence_control_set_ptr);
            }


            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag)
            {

                // Get Empty PicMgr Results
                eb_get_empty_object(
                    context_ptr->picture_demux_fifo_ptr,
                    &picture_demux_results_wrapper_ptr);

                picture_demux_results_rtr = (PictureDemuxResults_t*)picture_demux_results_wrapper_ptr->object_ptr;
                picture_demux_results_rtr->reference_picture_wrapper_ptr = picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr;
                picture_demux_results_rtr->sequence_control_set_wrapper_ptr = picture_control_set_ptr->sequence_control_set_wrapper_ptr;
                picture_demux_results_rtr->picture_number = picture_control_set_ptr->picture_number;
                picture_demux_results_rtr->pictureType = EB_PIC_REFERENCE;

                // Post Reference Picture
                eb_post_full_object(picture_demux_results_wrapper_ptr);
            }



            // Get Empty rest Results to EC
            eb_get_empty_object(
                context_ptr->rest_output_fifo_ptr,
                &rest_results_wrapper_ptr);
            rest_results_ptr = (struct RestResults_s*)rest_results_wrapper_ptr->object_ptr;
            rest_results_ptr->picture_control_set_wrapper_ptr = cdef_results_ptr->picture_control_set_wrapper_ptr;
            rest_results_ptr->completed_lcu_row_index_start = 0;
            rest_results_ptr->completed_lcu_row_count = ((sequence_control_set_ptr->luma_height + sequence_control_set_ptr->sb_size_pix - 1) >> lcuSizeLog2);
            // Post Rest Results
            eb_post_full_object(rest_results_wrapper_ptr);

#if REST_M
        }
        eb_release_mutex(picture_control_set_ptr->rest_search_mutex);
#endif


        // Release input Results
        eb_release_object(cdef_results_wrapper_ptr);

    }

    return EB_NULL;
}
#endif