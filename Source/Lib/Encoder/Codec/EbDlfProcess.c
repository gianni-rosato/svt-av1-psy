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
void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr, Bool is_highbd);
void svt_av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm,
                                                  int32_t after_cdef);
void svt_convert_pic_8bit_to_16bit(EbPictureBufferDesc *src_8bit, EbPictureBufferDesc *dst_16bit,
                                   uint16_t ss_x, uint16_t ss_y);

extern void get_recon_pic(PictureControlSet *pcs_ptr, EbPictureBufferDesc **recon_ptr,
                          Bool is_highbd);

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
        PictureParentControlSet *ppcs = pcs_ptr->parent_pcs_ptr;
        scs_ptr                       = pcs_ptr->scs_ptr;

        Bool is_16bit = scs_ptr->is_16bit_pipeline;
        if (is_16bit && scs_ptr->static_config.encoder_bit_depth == EB_EIGHT_BIT) {
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
        Bool           dlf_enable_flag = (Bool)pcs_ptr->parent_pcs_ptr->dlf_ctrls.enabled;
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
            EbPictureBufferDesc *recon_pic;
            get_recon_pic(pcs_ptr, &recon_pic, is_16bit);

            Av1Common *cm = pcs_ptr->parent_pcs_ptr->av1_cm;
            if (ppcs->enable_restoration) {
                link_eb_to_aom_buffer_desc(recon_pic,
                                           cm->frame_to_show,
                                           scs_ptr->max_input_pad_right,
                                           scs_ptr->max_input_pad_bottom,
                                           is_16bit);
                svt_av1_loop_restoration_save_boundary_lines(cm->frame_to_show, cm, 0);
            }

            if (scs_ptr->seq_header.cdef_level && pcs_ptr->parent_pcs_ptr->cdef_level) {
                const uint32_t offset_y = recon_pic->origin_x +
                    recon_pic->origin_y * recon_pic->stride_y;
                pcs_ptr->cdef_input_recon[0] = recon_pic->buffer_y + (offset_y << is_16bit);
                const uint32_t offset_cb     = (recon_pic->origin_x +
                                            recon_pic->origin_y * recon_pic->stride_cb) >>
                    1;
                pcs_ptr->cdef_input_recon[1] = recon_pic->buffer_cb + (offset_cb << is_16bit);
                const uint32_t offset_cr     = (recon_pic->origin_x +
                                            recon_pic->origin_y * recon_pic->stride_cr) >>
                    1;
                pcs_ptr->cdef_input_recon[2] = recon_pic->buffer_cr + (offset_cr << is_16bit);

                EbPictureBufferDesc *input_pic      = is_16bit
                         ? pcs_ptr->input_frame16bit
                         : pcs_ptr->parent_pcs_ptr->enhanced_picture_ptr;
                const uint32_t       input_offset_y = input_pic->origin_x +
                    input_pic->origin_y * input_pic->stride_y;
                pcs_ptr->cdef_input_source[0]  = input_pic->buffer_y + (input_offset_y << is_16bit);
                const uint32_t input_offset_cb = (input_pic->origin_x +
                                                  input_pic->origin_y * input_pic->stride_cb) >>
                    1;
                pcs_ptr->cdef_input_source[1] = input_pic->buffer_cb +
                    (input_offset_cb << is_16bit);
                const uint32_t input_offset_cr = (input_pic->origin_x +
                                                  input_pic->origin_y * input_pic->stride_cr) >>
                    1;
                pcs_ptr->cdef_input_source[2] = input_pic->buffer_cr +
                    (input_offset_cr << is_16bit);
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
