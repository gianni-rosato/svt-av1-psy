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
#include "EbDlfProcess.h"
#include "EbEncDecResults.h"
#include "EbEncDecTasks.h"
#include "EbPictureDemuxResults.h"
#include "EbReferenceObject.h"

#include "EbDeblockingFilter.h"

void av1_loop_restoration_save_boundary_lines(const Yv12BufferConfig *frame, Av1Common *cm, int32_t after_cdef);

/******************************************************
 * Dlf Context Constructor
 ******************************************************/
EbErrorType dlf_context_ctor(
    DlfContext_t **context_dbl_ptr,
    EbFifo_t                *dlf_input_fifo_ptr,
    EbFifo_t                *dlf_output_fifo_ptr ,
    EbBool                  is16bit,
    uint32_t                max_input_luma_width,
    uint32_t                max_input_luma_height
   )
{
    EbErrorType return_error = EB_ErrorNone;
    DlfContext_t *context_ptr;
    EB_MALLOC(DlfContext_t*, context_ptr, sizeof(DlfContext_t), EB_N_PTR);
    *context_dbl_ptr = context_ptr;

    // Input/Output System Resource Manager FIFOs
    context_ptr->dlf_input_fifo_ptr = dlf_input_fifo_ptr;
    context_ptr->dlf_output_fifo_ptr = dlf_output_fifo_ptr;



    context_ptr->temp_lf_recon_picture16bit_ptr = (EbPictureBufferDesc_t *)EB_NULL;
    context_ptr->temp_lf_recon_picture_ptr = (EbPictureBufferDesc_t *)EB_NULL;
    EbPictureBufferDescInitData_t temp_lf_recon_desc_init_data;
    temp_lf_recon_desc_init_data.maxWidth = (uint16_t)max_input_luma_width;
    temp_lf_recon_desc_init_data.maxHeight = (uint16_t)max_input_luma_height;
    temp_lf_recon_desc_init_data.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;

    temp_lf_recon_desc_init_data.left_padding = PAD_VALUE;
    temp_lf_recon_desc_init_data.right_padding = PAD_VALUE;
    temp_lf_recon_desc_init_data.top_padding = PAD_VALUE;
    temp_lf_recon_desc_init_data.bot_padding = PAD_VALUE;

    temp_lf_recon_desc_init_data.splitMode = EB_FALSE;

    if (is16bit) {
        temp_lf_recon_desc_init_data.bit_depth = EB_16BIT;
        return_error = EbReconPictureBufferDescCtor(
            (EbPtr*)&(context_ptr->temp_lf_recon_picture16bit_ptr),
            (EbPtr)&temp_lf_recon_desc_init_data);
    }
    else {
        temp_lf_recon_desc_init_data.bit_depth = EB_8BIT;
        return_error = EbReconPictureBufferDescCtor(
            (EbPtr*)&(context_ptr->temp_lf_recon_picture_ptr),
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
    DlfContext_t                            *context_ptr = (DlfContext_t*)input_ptr;
    PictureControlSet_t                     *picture_control_set_ptr;
    SequenceControlSet_t                    *sequence_control_set_ptr;

    //// Input
    EbObjectWrapper_t                       *enc_dec_results_wrapper_ptr;
    EncDecResults_t                         *enc_dec_results_ptr;

    //// Output
    EbObjectWrapper_t                       *dlf_results_wrapper_ptr;
    struct DlfResults_s*                     dlf_results_ptr;

    // SB Loop variables


    for (;;) {

        // Get EncDec Results
        EbGetFullObject(
            context_ptr->dlf_input_fifo_ptr,
            &enc_dec_results_wrapper_ptr);

        enc_dec_results_ptr = (EncDecResults_t*)enc_dec_results_wrapper_ptr->objectPtr;
        picture_control_set_ptr = (PictureControlSet_t*)enc_dec_results_ptr->pictureControlSetWrapperPtr->objectPtr;
        sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;

        EbBool  is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
        EbBool dlfEnableFlag = (EbBool)(picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode &&
            (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag ||
                sequence_control_set_ptr->static_config.recon_enabled ||
                sequence_control_set_ptr->static_config.stat_report));

        if (dlfEnableFlag && picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode >= 2) {

            EbPictureBufferDesc_t  *recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;
            if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE) {

                //get the 16bit form of the input LCU
                if (is16bit) {
                    recon_buffer = ((EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->objectPtr)->referencePicture16bit;
                }
                else {
                    recon_buffer = ((EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->objectPtr)->referencePicture;
                }
            }
            else { // non ref pictures
                recon_buffer = is16bit ? picture_control_set_ptr->recon_picture16bit_ptr : picture_control_set_ptr->recon_picture_ptr;
            }

            av1_loop_filter_init(picture_control_set_ptr);

            if (picture_control_set_ptr->parent_pcs_ptr->loop_filter_mode == 2) {

                av1_pick_filter_level(
                    context_ptr,
                    (EbPictureBufferDesc_t*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                    picture_control_set_ptr,
                    LPF_PICK_FROM_Q);

            }

            av1_pick_filter_level(
                context_ptr,
                (EbPictureBufferDesc_t*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr,
                picture_control_set_ptr,
                LPF_PICK_FROM_FULL_IMAGE);

#if NO_ENCDEC
            //NO DLF
            picture_control_set_ptr->parent_pcs_ptr->lf.filter_level[0] = 0;
            picture_control_set_ptr->parent_pcs_ptr->lf.filter_level[1] = 0;
            picture_control_set_ptr->parent_pcs_ptr->lf.filter_level_u = 0;
            picture_control_set_ptr->parent_pcs_ptr->lf.filter_level_v = 0;
#endif
                av1_loop_filter_frame(
                    recon_buffer,
                    picture_control_set_ptr,
                    0,
                    3);
            }


#if CDEF_M

        //pre-cdef prep
        {
            Av1Common* cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;
            EbPictureBufferDesc_t  * recon_picture_ptr;
            if (is16bit) {
                if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    recon_picture_ptr = ((EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->objectPtr)->referencePicture16bit;
                else
                    recon_picture_ptr = picture_control_set_ptr->recon_picture16bit_ptr;
            }
            else {
                if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
                    recon_picture_ptr = ((EbReferenceObject_t*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->objectPtr)->referencePicture;
                else
                    recon_picture_ptr = picture_control_set_ptr->recon_picture_ptr;
            }

            LinkEbToAomBufferDesc(
                recon_picture_ptr,
                cm->frame_to_show);

            if (sequence_control_set_ptr->enable_restoration) {
                av1_loop_restoration_save_boundary_lines(cm->frame_to_show, cm, 0);
            }

#if CDEF_M
            if (sequence_control_set_ptr->enable_cdef && picture_control_set_ptr->parent_pcs_ptr->cdef_filter_mode)
            {
#endif
                if (is16bit)
                {
                    picture_control_set_ptr->src[0] = (uint16_t*)recon_picture_ptr->bufferY + (recon_picture_ptr->origin_x + recon_picture_ptr->origin_y     * recon_picture_ptr->strideY);
                    picture_control_set_ptr->src[1] = (uint16_t*)recon_picture_ptr->bufferCb + (recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->strideCb);
                    picture_control_set_ptr->src[2] = (uint16_t*)recon_picture_ptr->bufferCr + (recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->strideCr);

                    EbPictureBufferDesc_t *inputPicturePtr = picture_control_set_ptr->input_frame16bit;
                    picture_control_set_ptr->ref_coeff[0] = (uint16_t*)inputPicturePtr->bufferY + (inputPicturePtr->origin_x + inputPicturePtr->origin_y * inputPicturePtr->strideY);
                    picture_control_set_ptr->ref_coeff[1] = (uint16_t*)inputPicturePtr->bufferCb + (inputPicturePtr->origin_x / 2 + inputPicturePtr->origin_y / 2 * inputPicturePtr->strideCb);
                    picture_control_set_ptr->ref_coeff[2] = (uint16_t*)inputPicturePtr->bufferCr + (inputPicturePtr->origin_x / 2 + inputPicturePtr->origin_y / 2 * inputPicturePtr->strideCr);

                }
                else
                {
                    //these copies should go!
                EbByte  rec_ptr = &((recon_picture_ptr->bufferY)[recon_picture_ptr->origin_x + recon_picture_ptr->origin_y * recon_picture_ptr->strideY]);
                    EbByte  rec_ptr_cb = &((recon_picture_ptr->bufferCb)[recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->strideCb]);
                    EbByte  rec_ptr_cr = &((recon_picture_ptr->bufferCr)[recon_picture_ptr->origin_x / 2 + recon_picture_ptr->origin_y / 2 * recon_picture_ptr->strideCr]);

                    EbPictureBufferDesc_t *inputPicturePtr = (EbPictureBufferDesc_t*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
                    EbByte  enh_ptr = &((inputPicturePtr->bufferY)[inputPicturePtr->origin_x + inputPicturePtr->origin_y * inputPicturePtr->strideY]);
                    EbByte  enh_ptr_cb = &((inputPicturePtr->bufferCb)[inputPicturePtr->origin_x / 2 + inputPicturePtr->origin_y / 2 * inputPicturePtr->strideCb]);
                    EbByte  enh_ptr_cr = &((inputPicturePtr->bufferCr)[inputPicturePtr->origin_x / 2 + inputPicturePtr->origin_y / 2 * inputPicturePtr->strideCr]);

                    for (int r = 0; r < sequence_control_set_ptr->luma_height; ++r) {
                        for (int c = 0; c < sequence_control_set_ptr->luma_width; ++c) {
                        picture_control_set_ptr->src[0]      [r * sequence_control_set_ptr->luma_width + c] = rec_ptr[r * recon_picture_ptr->strideY + c];
                        picture_control_set_ptr->ref_coeff[0][r * sequence_control_set_ptr->luma_width + c] = enh_ptr[r * inputPicturePtr->strideY + c];
                        }
                    }

                for (int r = 0; r < sequence_control_set_ptr->luma_height/2; ++r) {
                    for (int c = 0; c < sequence_control_set_ptr->luma_width/2; ++c) {
                        picture_control_set_ptr->src[1][r * sequence_control_set_ptr->luma_width/2 + c] = rec_ptr_cb[r * recon_picture_ptr->strideCb + c];
                        picture_control_set_ptr->ref_coeff[1][r * sequence_control_set_ptr->luma_width/2 + c] = enh_ptr_cb[r * inputPicturePtr->strideCb + c];
                            picture_control_set_ptr->src[2][r * sequence_control_set_ptr->luma_width / 2 + c] = rec_ptr_cr[r * recon_picture_ptr->strideCr + c];
                            picture_control_set_ptr->ref_coeff[2][r * sequence_control_set_ptr->luma_width / 2 + c] = enh_ptr_cr[r * inputPicturePtr->strideCr + c];
                        }
                    }

                }
#if CDEF_M
            }
#endif

        }

        picture_control_set_ptr->cdef_segments_column_count =  sequence_control_set_ptr->cdef_segment_column_count;
        picture_control_set_ptr->cdef_segments_row_count = sequence_control_set_ptr->cdef_segment_row_count;
        picture_control_set_ptr->cdef_segments_total_count  = (uint16_t)(picture_control_set_ptr->cdef_segments_column_count  * picture_control_set_ptr->cdef_segments_row_count);
        picture_control_set_ptr->tot_seg_searched_cdef = 0;
        uint32_t segment_index;

        for (segment_index = 0; segment_index < picture_control_set_ptr->cdef_segments_total_count; ++segment_index)
        {
            // Get Empty DLF Results to Cdef
            EbGetEmptyObject(
                context_ptr->dlf_output_fifo_ptr,
                &dlf_results_wrapper_ptr);
            dlf_results_ptr = (struct DlfResults_s*)dlf_results_wrapper_ptr->objectPtr;
            dlf_results_ptr->picture_control_set_wrapper_ptr = enc_dec_results_ptr->pictureControlSetWrapperPtr;

            dlf_results_ptr->segment_index = segment_index;
            // Post DLF Results
            EbPostFullObject(dlf_results_wrapper_ptr);
        }
#else

            // Get Empty DLF Results to Cdef
            EbGetEmptyObject(
                context_ptr->dlf_output_fifo_ptr,
                &dlf_results_wrapper_ptr);
            dlf_results_ptr = (struct DlfResults_s*)dlf_results_wrapper_ptr->objectPtr;
            dlf_results_ptr->pictureControlSetWrapperPtr = enc_dec_results_ptr->pictureControlSetWrapperPtr;
            dlf_results_ptr->completedLcuRowIndexStart = 0;
            dlf_results_ptr->completedLcuRowCount = ((sequence_control_set_ptr->luma_height + sequence_control_set_ptr->sb_size_pix - 1) >> lcuSizeLog2);
            // Post DLF Results
            EbPostFullObject(dlf_results_wrapper_ptr);
#endif

            // Release EncDec Results
            EbReleaseObject(enc_dec_results_wrapper_ptr);

        }

    return EB_NULL;
}
#endif