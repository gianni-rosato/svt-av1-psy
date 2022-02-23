/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>
#include "EbEncHandle.h"
#include "EbDlfProcess.h"
#include "EbEncDecResults.h"
#include "EbReferenceObject.h"
#include "EbDeblockingFilter.h"
#include "EbDefinitions.h"
#include "EbSequenceControlSet.h"
#include "EbPictureControlSet.h"
#include "aom_dsp_rtcd.h"
void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr, EbBool is_highbd);
void svt_av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm,
                                                  int32_t after_cdef);
void svt_convert_pic_8bit_to_16bit(EbPictureBufferDesc *src_8bit, EbPictureBufferDesc *dst_16bit,
                                   uint16_t ss_x, uint16_t ss_y);

extern void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr,
                          EbBool is_highbd);

static void dlf_context_dctor(EbPtr p) {
    EbThreadContext *thread_context_ptr = (EbThreadContext *)p;
    DlfContext      *obj                = (DlfContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}
/******************************************************
 * Dlf Context Constructor
 ******************************************************/
EbErrorType dlf_context_ctor(EbThreadContext *thread_context_ptr, const EbEncHandle *enc_handle_ptr,
                             int index) {
    DlfContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = dlf_context_dctor;

    // Input/Output System Resource Manager FIFOs
    context_ptr->dlf_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->enc_dec_results_resource_ptr, index);
    context_ptr->dlf_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->dlf_results_resource_ptr, index);
    return EB_ErrorNone;
}

/******************************************************
 * Dlf Kernel
 ******************************************************/
void *dlf_kernel(void *input_ptr) {
    // Context & SCS & PCS
    EbThreadContext    *thread_context_ptr = (EbThreadContext *)input_ptr;
    DlfContext         *context_ptr        = (DlfContext *)thread_context_ptr->priv;
    PictureControlSet  *pcs_ptr;
    SequenceControlSet *scs_ptr;

    //// Input
    EbObjectWrapper *enc_dec_results_wrapper_ptr;
    EncDecResults   *enc_dec_results_ptr;

    //// Output
    EbObjectWrapper   *dlf_results_wrapper_ptr;
    struct DlfResults *dlf_results_ptr;

    // SB Loop variables
    for (;;) {
        // Get EncDec Results
        EB_GET_FULL_OBJECT(context_ptr->dlf_input_fifo_ptr, &enc_dec_results_wrapper_ptr);

        enc_dec_results_ptr = (EncDecResults *)enc_dec_results_wrapper_ptr->object_ptr;
        pcs_ptr             = (PictureControlSet *)enc_dec_results_ptr->pcs_wrapper_ptr->object_ptr;
        scs_ptr             = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;

        EbBool is_16bit = scs_ptr->is_16bit_pipeline;
        if (is_16bit && scs_ptr->static_config.encoder_bit_depth == EB_8BIT) {
            svt_convert_pic_8bit_to_16bit(pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                                          pcs_ptr->input_frame16bit,
                                          pcs_ptr->parent_pcs_ptr->scs_ptr->subsampling_x,
                                          pcs_ptr->parent_pcs_ptr->scs_ptr->subsampling_y);
            // convert 8-bit recon to 16-bit for it bypass encdec process
            if (pcs_ptr->pic_bypass_encdec) {
                EbPictureBufferDesc *recon_picture_ptr;
                EbPictureBufferDesc *recon_picture_16bit_ptr;
                get_recon_pic(pcs_ptr, &recon_picture_ptr, 0);
                get_recon_pic(pcs_ptr, &recon_picture_16bit_ptr, 1);
                svt_convert_pic_8bit_to_16bit(recon_picture_ptr,
                                              recon_picture_16bit_ptr,
                                              pcs_ptr->parent_pcs_ptr->scs_ptr->subsampling_x,
                                              pcs_ptr->parent_pcs_ptr->scs_ptr->subsampling_y);
            }
        }
        EbBool         dlf_enable_flag = (EbBool)pcs_ptr->parent_pcs_ptr->dlf_ctrls.enabled;
        const uint16_t tg_count        = pcs_ptr->parent_pcs_ptr->tile_group_cols *
            pcs_ptr->parent_pcs_ptr->tile_group_rows;
        // Move sb level lf to here if tile_parallel
        if ((dlf_enable_flag && !pcs_ptr->parent_pcs_ptr->dlf_ctrls.sb_based_dlf) ||
            (dlf_enable_flag && pcs_ptr->parent_pcs_ptr->dlf_ctrls.sb_based_dlf && tg_count > 1)) {
            EbPictureBufferDesc *recon_buffer;
            get_recon_pic(pcs_ptr, &recon_buffer, is_16bit);
            svt_av1_loop_filter_init(pcs_ptr);
            svt_av1_pick_filter_level(
                (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                pcs_ptr,
                LPF_PICK_FROM_FULL_IMAGE);

            svt_av1_loop_filter_frame(recon_buffer, pcs_ptr, 0, 3);
        }

        //pre-cdef prep
        {
            Av1Common           *cm = pcs_ptr->parent_pcs_ptr->av1_cm;
            EbPictureBufferDesc *recon_picture_ptr;
            get_recon_pic(pcs_ptr, &recon_picture_ptr, is_16bit);
            link_eb_to_aom_buffer_desc(recon_picture_ptr,
                                       cm->frame_to_show,
                                       scs_ptr->max_input_pad_right,
                                       scs_ptr->max_input_pad_bottom,
                                       is_16bit);
            if (scs_ptr->seq_header.enable_restoration)
                svt_av1_loop_restoration_save_boundary_lines(cm->frame_to_show, cm, 0);
            if (scs_ptr->seq_header.cdef_level && pcs_ptr->parent_pcs_ptr->cdef_level) {
                if (is_16bit) {
                    pcs_ptr->src[0] = (uint16_t *)recon_picture_ptr->buffer_y +
                        (recon_picture_ptr->origin_x +
                         recon_picture_ptr->origin_y * recon_picture_ptr->stride_y);
                    pcs_ptr->src[1] = (uint16_t *)recon_picture_ptr->buffer_cb +
                        (recon_picture_ptr->origin_x / 2 +
                         recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb);
                    pcs_ptr->src[2] = (uint16_t *)recon_picture_ptr->buffer_cr +
                        (recon_picture_ptr->origin_x / 2 +
                         recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr);

                    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->input_frame16bit;
                    pcs_ptr->ref_coeff[0] = (uint16_t *)input_picture_ptr->buffer_y +
                        (input_picture_ptr->origin_x +
                         input_picture_ptr->origin_y * input_picture_ptr->stride_y);
                    pcs_ptr->ref_coeff[1] = (uint16_t *)input_picture_ptr->buffer_cb +
                        (input_picture_ptr->origin_x / 2 +
                         input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb);
                    pcs_ptr->ref_coeff[2] = (uint16_t *)input_picture_ptr->buffer_cr +
                        (input_picture_ptr->origin_x / 2 +
                         input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr);
                } else {
                    EbByte rec_ptr    = &((
                        recon_picture_ptr
                            ->buffer_y)[recon_picture_ptr->origin_x +
                                        recon_picture_ptr->origin_y * recon_picture_ptr->stride_y]);
                    EbByte rec_ptr_cb = &(
                        (recon_picture_ptr->buffer_cb)[recon_picture_ptr->origin_x / 2 +
                                                       recon_picture_ptr->origin_y / 2 *
                                                           recon_picture_ptr->stride_cb]);
                    EbByte rec_ptr_cr = &(
                        (recon_picture_ptr->buffer_cr)[recon_picture_ptr->origin_x / 2 +
                                                       recon_picture_ptr->origin_y / 2 *
                                                           recon_picture_ptr->stride_cr]);

                    EbPictureBufferDesc *input_picture_ptr =
                        (EbPictureBufferDesc *)pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
                    EbByte enh_ptr    = &((
                        input_picture_ptr
                            ->buffer_y)[input_picture_ptr->origin_x +
                                        input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
                    EbByte enh_ptr_cb = &(
                        (input_picture_ptr->buffer_cb)[input_picture_ptr->origin_x / 2 +
                                                       input_picture_ptr->origin_y / 2 *
                                                           input_picture_ptr->stride_cb]);
                    EbByte enh_ptr_cr = &(
                        (input_picture_ptr->buffer_cr)[input_picture_ptr->origin_x / 2 +
                                                       input_picture_ptr->origin_y / 2 *
                                                           input_picture_ptr->stride_cr]);

                    pcs_ptr->src[0] = (uint16_t *)rec_ptr;
                    pcs_ptr->src[1] = (uint16_t *)rec_ptr_cb;
                    pcs_ptr->src[2] = (uint16_t *)rec_ptr_cr;

                    pcs_ptr->ref_coeff[0] = (uint16_t *)enh_ptr;
                    pcs_ptr->ref_coeff[1] = (uint16_t *)enh_ptr_cb;
                    pcs_ptr->ref_coeff[2] = (uint16_t *)enh_ptr_cr;
                }
            }
        }

        pcs_ptr->cdef_segments_column_count = scs_ptr->cdef_segment_column_count;
        pcs_ptr->cdef_segments_row_count    = scs_ptr->cdef_segment_row_count;
        pcs_ptr->cdef_segments_total_count  = (uint16_t)(pcs_ptr->cdef_segments_column_count *
                                                        pcs_ptr->cdef_segments_row_count);
        pcs_ptr->tot_seg_searched_cdef      = 0;
        uint32_t segment_index;

        for (segment_index = 0; segment_index < pcs_ptr->cdef_segments_total_count;
             ++segment_index) {
            // Get Empty DLF Results to Cdef
            svt_get_empty_object(context_ptr->dlf_output_fifo_ptr, &dlf_results_wrapper_ptr);
            dlf_results_ptr = (struct DlfResults *)dlf_results_wrapper_ptr->object_ptr;
            dlf_results_ptr->pcs_wrapper_ptr = enc_dec_results_ptr->pcs_wrapper_ptr;
            dlf_results_ptr->segment_index   = segment_index;
            // Post DLF Results
            svt_post_full_object(dlf_results_wrapper_ptr);
        }

        // Release EncDec Results
        svt_release_object(enc_dec_results_wrapper_ptr);
    }

    return NULL;
}
