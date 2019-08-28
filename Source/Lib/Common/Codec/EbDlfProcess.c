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
#include "EbDlfProcess.h"
#include "EbEncDecResults.h"
#include "EbThreads.h"
#include "EbReferenceObject.h"

#include "EbDeblockingFilter.h"

void eb_av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm, int32_t after_cdef);

static void dlf_context_dctor(EbPtr p)
{
    DlfContext *obj = (DlfContext*)p;
    EB_DELETE(obj->temp_lf_recon_picture_ptr);
    EB_DELETE(obj->temp_lf_recon_picture16bit_ptr);
}
/******************************************************
 * Dlf Context Constructor
 ******************************************************/
EbErrorType dlf_context_ctor(
    DlfContext            *context_ptr,
    EbFifo                *dlf_input_fifo_ptr,
    EbFifo                *dlf_output_fifo_ptr ,
    EbBool                  is16bit,
    EbColorFormat           color_format,
    uint32_t                max_input_luma_width,
    uint32_t                max_input_luma_height
   )
{
    EbErrorType return_error = EB_ErrorNone;
    context_ptr->dctor = dlf_context_dctor;

    // Input/Output System Resource Manager FIFOs
    context_ptr->dlf_input_fifo_ptr = dlf_input_fifo_ptr;
    context_ptr->dlf_output_fifo_ptr = dlf_output_fifo_ptr;

    context_ptr->temp_lf_recon_picture16bit_ptr = (EbPictureBufferDesc *)EB_NULL;
    context_ptr->temp_lf_recon_picture_ptr = (EbPictureBufferDesc *)EB_NULL;
    EbPictureBufferDescInitData temp_lf_recon_desc_init_data;
    temp_lf_recon_desc_init_data.max_width = (uint16_t)max_input_luma_width;
    temp_lf_recon_desc_init_data.max_height = (uint16_t)max_input_luma_height;
    temp_lf_recon_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

    temp_lf_recon_desc_init_data.left_padding = PAD_VALUE;
    temp_lf_recon_desc_init_data.right_padding = PAD_VALUE;
    temp_lf_recon_desc_init_data.top_padding = PAD_VALUE;
    temp_lf_recon_desc_init_data.bot_padding = PAD_VALUE;

    temp_lf_recon_desc_init_data.split_mode = EB_FALSE;
    temp_lf_recon_desc_init_data.color_format = color_format;

    if (is16bit) {
        temp_lf_recon_desc_init_data.bit_depth = EB_16BIT;
        EB_NEW(
            context_ptr->temp_lf_recon_picture16bit_ptr,
            eb_recon_picture_buffer_desc_ctor,
            (EbPtr)&temp_lf_recon_desc_init_data);
    }
    else {
        temp_lf_recon_desc_init_data.bit_depth = EB_8BIT;
        EB_NEW(
            context_ptr->temp_lf_recon_picture_ptr,
            eb_recon_picture_buffer_desc_ctor,
            (EbPtr)&temp_lf_recon_desc_init_data);
    }

    return return_error;
}

/******************************************************
 * Dlf Kernel
 ******************************************************/
void* dlf_kernel(void *input_ptr)
{
    // Context & SCS & PCS
    DlfContext                            *context_ptr = (DlfContext*)input_ptr;
    PictureControlSet                     *picture_control_set_ptr;
    SequenceControlSet                    *sequence_control_set_ptr;

    //// Input
    EbObjectWrapper                       *enc_dec_results_wrapper_ptr;
    EncDecResults                         *enc_dec_results_ptr;

    //// Output
    EbObjectWrapper                       *dlf_results_wrapper_ptr;
    struct DlfResults*                     dlf_results_ptr;

    // SB Loop variables
    for (;;) {
        // Get EncDec Results
        eb_get_full_object(
            context_ptr->dlf_input_fifo_ptr,
            &enc_dec_results_wrapper_ptr);

        enc_dec_results_ptr         = (EncDecResults*)enc_dec_results_wrapper_ptr->object_ptr;
        picture_control_set_ptr     = (PictureControlSet*)enc_dec_results_ptr->picture_control_set_wrapper_ptr->object_ptr;
        sequence_control_set_ptr    = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;

        EbBool is16bit       = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

        EbBool dlfEnableFlag = (EbBool) picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode;
        if (dlfEnableFlag && picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode >= 2) {
            EbPictureBufferDesc  *recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;

            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                //get the 16bit form of the input LCU
                if (is16bit)
                    recon_buffer = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
                else
                    recon_buffer = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
            else  // non ref pictures
                recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;

            eb_av1_loop_filter_init(picture_control_set_ptr);

            if (picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode == 2) {
                eb_av1_pick_filter_level(
                    context_ptr,
                    (EbPictureBufferDesc*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                    picture_control_set_ptr,
                    LPF_PICK_FROM_Q);
            }

            eb_av1_pick_filter_level(
                context_ptr,
                (EbPictureBufferDesc*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                picture_control_set_ptr,
                LPF_PICK_FROM_FULL_IMAGE);

#if NO_ENCDEC
            //NO DLF
            picture_control_set_ptr->parent_pcs_ptr->lf.filter_level[0] = 0;
            picture_control_set_ptr->parent_pcs_ptr->lf.filter_level[1] = 0;
            picture_control_set_ptr->parent_pcs_ptr->lf.filter_level_u = 0;
            picture_control_set_ptr->parent_pcs_ptr->lf.filter_level_v = 0;
#endif
                eb_av1_loop_filter_frame(
                    recon_buffer,
                    picture_control_set_ptr,
                    0,
                    3);
            }

        //pre-cdef prep
        {
            Av1Common* cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;
            EbPictureBufferDesc  * recon_picture_ptr;
            if (is16bit) {
                if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    recon_picture_ptr = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
                else
                    recon_picture_ptr = picture_control_set_ptr->recon_picture16bit_ptr;
            }
            else {
                if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    recon_picture_ptr = ((EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
                else
                    recon_picture_ptr = picture_control_set_ptr->recon_picture_ptr;
            }

            link_eb_to_aom_buffer_desc(
                recon_picture_ptr,
                cm->frame_to_show);

            if (sequence_control_set_ptr->seq_header.enable_restoration)
                eb_av1_loop_restoration_save_boundary_lines(cm->frame_to_show, cm, 0);
            if (sequence_control_set_ptr->seq_header.enable_cdef && picture_control_set_ptr->parent_pcs_ptr->cdef_filter_mode)
            {
                if (is16bit)
                {
                    picture_control_set_ptr->src[0] = (uint16_t*)recon_picture_ptr->buffer_y + (recon_picture_ptr->origin_x + recon_picture_ptr->origin_y     * recon_picture_ptr->stride_y);
                    picture_control_set_ptr->src[1] = (uint16_t*)recon_picture_ptr->buffer_cb + (recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb);
                    picture_control_set_ptr->src[2] = (uint16_t*)recon_picture_ptr->buffer_cr + (recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr);

                    EbPictureBufferDesc *input_picture_ptr = picture_control_set_ptr->input_frame16bit;
                    picture_control_set_ptr->ref_coeff[0] = (uint16_t*)input_picture_ptr->buffer_y + (input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y);
                    picture_control_set_ptr->ref_coeff[1] = (uint16_t*)input_picture_ptr->buffer_cb + (input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb);
                    picture_control_set_ptr->ref_coeff[2] = (uint16_t*)input_picture_ptr->buffer_cr + (input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr);
                }
                else
                {
                    EbByte  rec_ptr = &((recon_picture_ptr->buffer_y)[recon_picture_ptr->origin_x + recon_picture_ptr->origin_y * recon_picture_ptr->stride_y]);
                    EbByte  rec_ptr_cb = &((recon_picture_ptr->buffer_cb)[recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cb]);
                    EbByte  rec_ptr_cr = &((recon_picture_ptr->buffer_cr)[recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->stride_cr]);

                    EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
                    EbByte  enh_ptr = &((input_picture_ptr->buffer_y)[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y]);
                    EbByte  enh_ptr_cb = &((input_picture_ptr->buffer_cb)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);
                    EbByte  enh_ptr_cr = &((input_picture_ptr->buffer_cr)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);

                    picture_control_set_ptr->src[0] = (uint16_t*)rec_ptr;
                    picture_control_set_ptr->src[1] = (uint16_t*)rec_ptr_cb;
                    picture_control_set_ptr->src[2] = (uint16_t*)rec_ptr_cr;

                    picture_control_set_ptr->ref_coeff[0] = (uint16_t*)enh_ptr;
                    picture_control_set_ptr->ref_coeff[1] = (uint16_t*)enh_ptr_cb;
                    picture_control_set_ptr->ref_coeff[2] = (uint16_t*)enh_ptr_cr;

                }
            }
        }

        picture_control_set_ptr->cdef_segments_column_count =  sequence_control_set_ptr->cdef_segment_column_count;
        picture_control_set_ptr->cdef_segments_row_count    = sequence_control_set_ptr->cdef_segment_row_count;
        picture_control_set_ptr->cdef_segments_total_count  = (uint16_t)(picture_control_set_ptr->cdef_segments_column_count  * picture_control_set_ptr->cdef_segments_row_count);
        picture_control_set_ptr->tot_seg_searched_cdef      = 0;
        uint32_t segment_index;

        for (segment_index = 0; segment_index < picture_control_set_ptr->cdef_segments_total_count; ++segment_index)
        {
            // Get Empty DLF Results to Cdef
            eb_get_empty_object(
                context_ptr->dlf_output_fifo_ptr,
                &dlf_results_wrapper_ptr);
            dlf_results_ptr = (struct DlfResults*)dlf_results_wrapper_ptr->object_ptr;
            dlf_results_ptr->picture_control_set_wrapper_ptr = enc_dec_results_ptr->picture_control_set_wrapper_ptr;
            dlf_results_ptr->segment_index = segment_index;
            // Post DLF Results
            eb_post_full_object(dlf_results_wrapper_ptr);
        }

            // Release EncDec Results
            eb_release_object(enc_dec_results_wrapper_ptr);
        }

    return EB_NULL;
}
