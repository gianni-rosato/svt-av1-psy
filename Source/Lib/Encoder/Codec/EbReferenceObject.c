/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>
#include <string.h>

#include "EbThreads.h"
#include "EbReferenceObject.h"
#include "EbPictureBufferDesc.h"
#include "EbUtility.h"
#include "EncModeConfig.h"

void initialize_samples_neighboring_reference_picture_8bit(EbByte recon_samples_buffer_ptr, uint16_t stride,
                                                           uint16_t recon_width, uint16_t recon_height,
                                                           uint16_t left_padding, uint16_t top_padding) {
    uint8_t *recon_samples_ptr;
    uint16_t sample_count;

    // 1. zero out the top row
    recon_samples_ptr = recon_samples_buffer_ptr + (top_padding - 1) * stride + left_padding - 1;
    EB_MEMSET(recon_samples_ptr, 0, sizeof(uint8_t) * (1 + recon_width + 1));

    // 2. zero out the bottom row
    recon_samples_ptr = recon_samples_buffer_ptr + (top_padding + recon_height) * stride + left_padding - 1;
    EB_MEMSET(recon_samples_ptr, 0, sizeof(uint8_t) * (1 + recon_width + 1));

    // 3. zero out the left column
    recon_samples_ptr = recon_samples_buffer_ptr + top_padding * stride + left_padding - 1;
    for (sample_count = 0; sample_count < recon_height; sample_count++) recon_samples_ptr[sample_count * stride] = 0;
    // 4. zero out the right column
    recon_samples_ptr = recon_samples_buffer_ptr + top_padding * stride + left_padding + recon_width;
    for (sample_count = 0; sample_count < recon_height; sample_count++) recon_samples_ptr[sample_count * stride] = 0;
}

static void initialize_samples_neighboring_reference_picture(
    EbReferenceObject *ref_object, EbPictureBufferDescInitData *picture_buffer_desc_init_data_ptr,
    EbBitDepth bit_depth) {
    UNUSED(bit_depth);
    {
        initialize_samples_neighboring_reference_picture_8bit(ref_object->reference_picture->buffer_y,
                                                              ref_object->reference_picture->stride_y,
                                                              ref_object->reference_picture->width,
                                                              ref_object->reference_picture->height,
                                                              picture_buffer_desc_init_data_ptr->left_padding,
                                                              picture_buffer_desc_init_data_ptr->top_padding);

        initialize_samples_neighboring_reference_picture_8bit(ref_object->reference_picture->buffer_cb,
                                                              ref_object->reference_picture->stride_cb,
                                                              ref_object->reference_picture->width >> 1,
                                                              ref_object->reference_picture->height >> 1,
                                                              picture_buffer_desc_init_data_ptr->left_padding >> 1,
                                                              picture_buffer_desc_init_data_ptr->top_padding >> 1);

        initialize_samples_neighboring_reference_picture_8bit(ref_object->reference_picture->buffer_cr,
                                                              ref_object->reference_picture->stride_cr,
                                                              ref_object->reference_picture->width >> 1,
                                                              ref_object->reference_picture->height >> 1,
                                                              picture_buffer_desc_init_data_ptr->left_padding >> 1,
                                                              picture_buffer_desc_init_data_ptr->top_padding >> 1);
    }
}

static void svt_reference_object_dctor(EbPtr p) {
    EbReferenceObject *obj = (EbReferenceObject *)p;

    EB_DELETE(obj->reference_picture);
    EB_FREE_2D(obj->unit_info);
    EB_FREE_ALIGNED_ARRAY(obj->mvs);
    EB_FREE_ARRAY(obj->sb_intra);
    EB_FREE_ARRAY(obj->sb_skip);
    EB_FREE_ARRAY(obj->sb_64x64_mvp);
    EB_FREE_ARRAY(obj->sb_me_64x64_dist);
    EB_FREE_ARRAY(obj->sb_me_8x8_cost_var);
    for (uint8_t sr_denom_idx = 0; sr_denom_idx < NUM_SR_SCALES + 1; sr_denom_idx++) {
        for (uint8_t resize_denom_idx = 0; resize_denom_idx < NUM_RESIZE_SCALES + 1; resize_denom_idx++) {
            if (obj->downscaled_reference_picture[sr_denom_idx][resize_denom_idx] != NULL) {
                EB_DELETE(obj->downscaled_reference_picture[sr_denom_idx][resize_denom_idx]);
            }
            EB_DESTROY_MUTEX(obj->resize_mutex[sr_denom_idx][resize_denom_idx]);
        }
    }
    EB_DELETE(obj->quarter_reference_picture);
    EB_DELETE(obj->sixteenth_reference_picture);
    EB_DELETE(obj->input_picture);
    EB_DELETE(obj->quarter_input_picture);
    EB_DELETE(obj->sixteenth_input_picture);
}

/*****************************************
 * svt_picture_buffer_desc_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType svt_reference_object_ctor(EbReferenceObject *ref_object, EbPtr object_init_data_ptr) {
    EbReferenceObjectDescInitData *ref_init_ptr = (EbReferenceObjectDescInitData *)object_init_data_ptr;
    EbPictureBufferDescInitData   *picture_buffer_desc_init_data_ptr = &ref_init_ptr->reference_picture_desc_init_data;
    EbPictureBufferDescInitData    picture_buffer_desc_init_data_16bit_ptr = *picture_buffer_desc_init_data_ptr;

    ref_object->dctor = svt_reference_object_dctor;
    //TODO:12bit
    if (picture_buffer_desc_init_data_16bit_ptr.bit_depth == EB_TEN_BIT) {
        // Hsan: set split_mode to 0 to construct the packed reference buffer (used @ EP)
        // Use 10bit here to use in MD
        picture_buffer_desc_init_data_16bit_ptr.split_mode = TRUE;
        picture_buffer_desc_init_data_16bit_ptr.bit_depth  = EB_TEN_BIT;
        EB_NEW(ref_object->reference_picture,
               svt_picture_buffer_desc_ctor,
               (EbPtr)&picture_buffer_desc_init_data_16bit_ptr);
    } else {
        // Hsan: set split_mode to 0 to as 8BIT input
        picture_buffer_desc_init_data_ptr->split_mode = FALSE;
        EB_NEW(ref_object->reference_picture, svt_picture_buffer_desc_ctor, (EbPtr)picture_buffer_desc_init_data_ptr);

        initialize_samples_neighboring_reference_picture(
            ref_object, picture_buffer_desc_init_data_ptr, picture_buffer_desc_init_data_16bit_ptr.bit_depth);
    }
    ref_object->input_picture = NULL;
    const bool gm_ref_info    = svt_aom_need_gm_ref_info(
#if TUNE_M6
        ref_init_ptr->static_config->enc_mode, ref_init_ptr->static_config->resize_mode == RESIZE_NONE);
#else
        ref_init_ptr->static_config->enc_mode, true, ref_init_ptr->static_config->resize_mode == RESIZE_NONE);
#endif
    if (gm_ref_info)
        EB_NEW(ref_object->input_picture, svt_picture_buffer_desc_ctor, (EbPtr)picture_buffer_desc_init_data_ptr);
    uint32_t mi_rows = ref_object->reference_picture->height >> MI_SIZE_LOG2;
    uint32_t mi_cols = ref_object->reference_picture->width >> MI_SIZE_LOG2;
    // there should be one unit info per plane and per rest unit
    EB_MALLOC_2D(ref_object->unit_info, MAX_MB_PLANE, picture_buffer_desc_init_data_ptr->rest_units_per_tile);

    if (picture_buffer_desc_init_data_ptr->mfmv) {
        //MFMV map is 8x8 based.
        const int mem_size = ((mi_rows + 1) >> 1) * ((mi_cols + 1) >> 1);
        EB_CALLOC_ALIGNED_ARRAY(ref_object->mvs, mem_size);
    }
    ref_object->quarter_reference_picture   = NULL;
    ref_object->sixteenth_reference_picture = NULL;
    if (gm_ref_info) {
        EbPictureBufferDescInitData buf_desc;
        buf_desc.max_width          = picture_buffer_desc_init_data_ptr->max_width >> 1;
        buf_desc.max_height         = picture_buffer_desc_init_data_ptr->max_height >> 1;
        buf_desc.bit_depth          = 8; //Should be 8bit
        buf_desc.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        buf_desc.left_padding       = 32;
        buf_desc.right_padding      = 32;
        buf_desc.top_padding        = 32;
        buf_desc.bot_padding        = 32;
        buf_desc.color_format       = EB_YUV420;
        buf_desc.split_mode         = FALSE;

        EB_NEW(ref_object->quarter_reference_picture, svt_picture_buffer_desc_ctor, (EbPtr)&buf_desc);

        buf_desc.max_width          = picture_buffer_desc_init_data_ptr->max_width >> 2;
        buf_desc.max_height         = picture_buffer_desc_init_data_ptr->max_height >> 2;
        buf_desc.bit_depth          = 8; //Should be 8bit
        buf_desc.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
        buf_desc.left_padding       = 16;
        buf_desc.right_padding      = 16;
        buf_desc.top_padding        = 16;
        buf_desc.bot_padding        = 16;
        buf_desc.color_format       = EB_YUV420;
        buf_desc.split_mode         = FALSE;

        EB_NEW(ref_object->sixteenth_reference_picture, svt_picture_buffer_desc_ctor, (EbPtr)&buf_desc);
    }
    ref_object->ds_pics.picture_ptr           = ref_object->reference_picture;
    ref_object->ds_pics.quarter_picture_ptr   = ref_object->quarter_reference_picture;
    ref_object->ds_pics.sixteenth_picture_ptr = ref_object->sixteenth_reference_picture;
    memset(&ref_object->film_grain_params, 0, sizeof(ref_object->film_grain_params));
    // set all supplemental downscaled reference picture pointers to NULL
    for (uint8_t sr_denom_idx = 0; sr_denom_idx < NUM_SR_SCALES + 1; sr_denom_idx++) {
        for (uint8_t resize_denom_idx = 0; resize_denom_idx < NUM_RESIZE_SCALES + 1; resize_denom_idx++) {
            ref_object->downscaled_reference_picture[sr_denom_idx][resize_denom_idx] = NULL;
            ref_object->downscaled_picture_number[sr_denom_idx][resize_denom_idx]    = (uint64_t)~0;
            EB_CREATE_MUTEX(ref_object->resize_mutex[sr_denom_idx][resize_denom_idx]);
        }
    }

    ref_object->mi_rows = mi_rows;
    ref_object->mi_cols = mi_cols;
    EB_MALLOC_ARRAY(ref_object->sb_intra, picture_buffer_desc_init_data_ptr->sb_total_count);
    EB_MALLOC_ARRAY(ref_object->sb_skip, picture_buffer_desc_init_data_ptr->sb_total_count);
    EB_MALLOC_ARRAY(ref_object->sb_64x64_mvp, picture_buffer_desc_init_data_ptr->sb_total_count);
    EB_MALLOC_ARRAY(ref_object->sb_me_64x64_dist, picture_buffer_desc_init_data_ptr->sb_total_count);
    EB_MALLOC_ARRAY(ref_object->sb_me_8x8_cost_var, picture_buffer_desc_init_data_ptr->sb_total_count);
    return EB_ErrorNone;
}

EbErrorType svt_reference_object_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    EbReferenceObject *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, svt_reference_object_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

EbErrorType svt_reference_object_reset(EbReferenceObject *ref_object, SequenceControlSet *scs) {
    ref_object->mi_rows = scs->max_input_luma_height >> MI_SIZE_LOG2;
    ref_object->mi_cols = scs->max_input_luma_width >> MI_SIZE_LOG2;

    return EB_ErrorNone;
}

static void svt_pa_reference_object_dctor(EbPtr p) {
    EbPaReferenceObject *obj = (EbPaReferenceObject *)p;
    if (obj->dummy_obj)
        return;
    EB_DELETE(obj->input_padded_pic);
    EB_DELETE(obj->quarter_downsampled_picture_ptr);
    EB_DELETE(obj->sixteenth_downsampled_picture_ptr);
    for (uint8_t sr_denom_idx = 0; sr_denom_idx < NUM_SR_SCALES + 1; sr_denom_idx++) {
        for (uint8_t resize_denom_idx = 0; resize_denom_idx < NUM_RESIZE_SCALES + 1; resize_denom_idx++) {
            if (obj->downscaled_input_padded_picture_ptr[sr_denom_idx][resize_denom_idx] != NULL) {
                EB_DELETE(obj->downscaled_input_padded_picture_ptr[sr_denom_idx][resize_denom_idx]);
                EB_DELETE(obj->downscaled_quarter_downsampled_picture_ptr[sr_denom_idx][resize_denom_idx]);
                EB_DELETE(obj->downscaled_sixteenth_downsampled_picture_ptr[sr_denom_idx][resize_denom_idx]);
            }
            EB_DESTROY_MUTEX(obj->resize_mutex[sr_denom_idx][resize_denom_idx]);
        }
    }
}

static void svt_tpl_reference_object_dctor(EbPtr p) {
    EbTplReferenceObject *obj = (EbTplReferenceObject *)p;
    EB_DELETE(obj->ref_picture_ptr);
}

/*****************************************
 * svt_pa_reference_object_ctor
 *  Initializes the Buffer Descriptor's
 *  values that are fixed for the life of
 *  the descriptor.
 *****************************************/
EbErrorType svt_pa_reference_object_ctor(EbPaReferenceObject *pa_ref_obj_, EbPtr object_init_data_ptr) {
    EbPictureBufferDescInitData *picture_buffer_desc_init_data_ptr = (EbPictureBufferDescInitData *)
        object_init_data_ptr;

    pa_ref_obj_->dctor = svt_pa_reference_object_dctor;

    // Reference picture constructor
    EB_NEW(pa_ref_obj_->input_padded_pic, svt_picture_buffer_desc_ctor, (EbPtr)picture_buffer_desc_init_data_ptr);
    // Downsampled reference picture constructor
    EB_NEW(pa_ref_obj_->quarter_downsampled_picture_ptr,
           svt_picture_buffer_desc_ctor,
           (EbPtr)(picture_buffer_desc_init_data_ptr + 1));
    EB_NEW(pa_ref_obj_->sixteenth_downsampled_picture_ptr,
           svt_picture_buffer_desc_ctor,
           (EbPtr)(picture_buffer_desc_init_data_ptr + 2));
    // set all supplemental downscaled reference picture pointers to NULL
    for (uint8_t sr_down_idx = 0; sr_down_idx < NUM_SR_SCALES + 1; sr_down_idx++) {
        for (uint8_t resize_down_idx = 0; resize_down_idx < NUM_RESIZE_SCALES + 1; resize_down_idx++) {
            pa_ref_obj_->downscaled_input_padded_picture_ptr[sr_down_idx][resize_down_idx]          = NULL;
            pa_ref_obj_->downscaled_quarter_downsampled_picture_ptr[sr_down_idx][resize_down_idx]   = NULL;
            pa_ref_obj_->downscaled_sixteenth_downsampled_picture_ptr[sr_down_idx][resize_down_idx] = NULL;
            pa_ref_obj_->downscaled_picture_number[sr_down_idx][resize_down_idx]                    = (uint64_t)~0;
            EB_CREATE_MUTEX(pa_ref_obj_->resize_mutex[sr_down_idx][resize_down_idx]);
        }
    }

    return EB_ErrorNone;
}

EbErrorType svt_pa_reference_object_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    EbPaReferenceObject *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, svt_pa_reference_object_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

EbErrorType svt_tpl_reference_object_ctor(EbTplReferenceObject *tpl_ref_obj_, EbPtr object_init_data_ptr) {
    EbPictureBufferDescInitData *picture_buffer_desc_init_data_ptr = (EbPictureBufferDescInitData *)
        object_init_data_ptr;

    tpl_ref_obj_->dctor = svt_tpl_reference_object_dctor;

    // Reference picture constructor
    EB_NEW(tpl_ref_obj_->ref_picture_ptr, svt_picture_buffer_desc_ctor, (EbPtr)picture_buffer_desc_init_data_ptr);

    return EB_ErrorNone;
}
EbErrorType svt_tpl_reference_object_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    EbTplReferenceObject *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, svt_tpl_reference_object_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}

/************************************************
* Release Pa Reference Objects
** Check if reference pictures are needed
** release them when appropriate
************************************************/
void svt_aom_release_pa_reference_objects(SequenceControlSet *scs, PictureParentControlSet *pcs) {
    (void)scs;
    // PA Reference Pictures
    if (pcs->slice_type != I_SLICE) {
        const uint32_t num_of_list_to_search = (pcs->slice_type == P_SLICE) ? 1 /*List 0 only*/ : 2 /*List 0 + 1*/;

        // List Loop
        for (uint32_t list_index = REF_LIST_0; list_index < num_of_list_to_search; ++list_index) {
            // Release PA Reference Pictures
            uint8_t num_of_ref_pic_to_search = (list_index == REF_LIST_0) ? pcs->ref_list0_count : pcs->ref_list1_count;

            for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
                if (pcs->ref_pa_pic_ptr_array[list_index][ref_pic_index] != NULL) {
                    //assert((int32_t)pcs->ref_pa_pic_ptr_array[list_index][ref_pic_index]->live_count > 0);
                    svt_release_object(pcs->ref_pa_pic_ptr_array[list_index][ref_pic_index]);

                    if (pcs->ref_y8b_array[list_index][ref_pic_index]) {
                        //y8b  needs to get decremented at the same time of pa ref
                        svt_release_object(pcs->ref_y8b_array[list_index][ref_pic_index]);
                    }
                }
            }
        }
    }

    if (pcs->pa_ref_pic_wrapper != NULL) {
        //assert((int32_t)pcs->pa_ref_pic_wrapper->live_count > 0);
        svt_release_object(pcs->pa_ref_pic_wrapper);

        if (pcs->y8b_wrapper) {
            //y8b needs to get decremented at the same time of pa ref
            svt_release_object(pcs->y8b_wrapper);
        }
    }
    // Mark that the PCS released PA references
    pcs->reference_released = 1;
    return;
}
