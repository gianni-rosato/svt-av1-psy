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
#include "EbPictureControlSet.h"
#include "EbPictureBufferDesc.h"

#if CDEF_M
void *aom_memalign(size_t align, size_t size);
void aom_free(void *memblk);
void *aom_malloc(size_t size);

#endif
EbErrorType av1_alloc_restoration_buffers(Av1Common *cm);

#if ICOPY
EbErrorType av1_hash_table_create(hash_table *p_hash_table);
#endif

static void set_restoration_unit_size(int32_t width, int32_t height, int32_t sx, int32_t sy,
    RestorationInfo *rst) {
    (void)width;
    (void)height;
    (void)sx;
    (void)sy;

    int32_t s = 0;


    if (width * height > 352 * 288)
        rst[0].restoration_unit_size = RESTORATION_UNITSIZE_MAX;
    else
        rst[0].restoration_unit_size = (RESTORATION_UNITSIZE_MAX >> 1);
    rst[1].restoration_unit_size = rst[0].restoration_unit_size >> s;
    rst[2].restoration_unit_size = rst[1].restoration_unit_size;
}


EbErrorType picture_control_set_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    PictureControlSet_t *object_ptr;
    PictureControlSetInitData_t *initDataPtr = (PictureControlSetInitData_t*)object_init_data_ptr;

    EbPictureBufferDescInitData_t input_picture_buffer_desc_init_data;
    EbPictureBufferDescInitData_t coeffBufferDescInitData;

    // Max/Min CU Sizes
    const uint32_t maxCuSize = initDataPtr->sb_sz;

    // LCUs
    const uint16_t pictureLcuWidth = (uint16_t)((initDataPtr->picture_width + initDataPtr->sb_sz - 1) / initDataPtr->sb_sz);
    const uint16_t pictureLcuHeight = (uint16_t)((initDataPtr->picture_height + initDataPtr->sb_sz - 1) / initDataPtr->sb_sz);
    uint16_t sb_index;
    uint16_t sb_origin_x;
    uint16_t sb_origin_y;
    EbErrorType return_error = EB_ErrorNone;

    EbBool is16bit = initDataPtr->is16bit;

    EB_MALLOC(PictureControlSet_t*, object_ptr, sizeof(PictureControlSet_t), EB_N_PTR);

    // Init Picture Init data
    input_picture_buffer_desc_init_data.maxWidth = initDataPtr->picture_width;
    input_picture_buffer_desc_init_data.maxHeight = initDataPtr->picture_height;
    input_picture_buffer_desc_init_data.bit_depth = initDataPtr->bit_depth;
    input_picture_buffer_desc_init_data.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;

    input_picture_buffer_desc_init_data.left_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.right_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.top_padding = PAD_VALUE;
    input_picture_buffer_desc_init_data.bot_padding = PAD_VALUE;

    input_picture_buffer_desc_init_data.splitMode = EB_FALSE;

    coeffBufferDescInitData.maxWidth = initDataPtr->picture_width;
    coeffBufferDescInitData.maxHeight = initDataPtr->picture_height;
    coeffBufferDescInitData.bit_depth = EB_16BIT;
    coeffBufferDescInitData.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;

    coeffBufferDescInitData.left_padding = PAD_VALUE;
    coeffBufferDescInitData.right_padding = PAD_VALUE;
    coeffBufferDescInitData.top_padding = PAD_VALUE;
    coeffBufferDescInitData.bot_padding = PAD_VALUE;

    coeffBufferDescInitData.splitMode = EB_FALSE;


    *object_dbl_ptr = (EbPtr)object_ptr;

    object_ptr->sequence_control_set_wrapper_ptr = (EbObjectWrapper_t *)EB_NULL;

    object_ptr->recon_picture16bit_ptr = (EbPictureBufferDesc_t *)EB_NULL;
    object_ptr->recon_picture_ptr = (EbPictureBufferDesc_t *)EB_NULL;

    EbPictureBufferDescInitData_t coeffBufferDes32bitInitData;
    coeffBufferDes32bitInitData.maxWidth = initDataPtr->picture_width;
    coeffBufferDes32bitInitData.maxHeight = initDataPtr->picture_height;
    coeffBufferDes32bitInitData.bit_depth = EB_32BIT;
    coeffBufferDes32bitInitData.bufferEnableMask = PICTURE_BUFFER_DESC_FULL_MASK;
    coeffBufferDes32bitInitData.left_padding = 0;
    coeffBufferDes32bitInitData.right_padding = 0;
    coeffBufferDes32bitInitData.top_padding = 0;
    coeffBufferDes32bitInitData.bot_padding = 0;
    coeffBufferDes32bitInitData.splitMode = EB_FALSE;

    object_ptr->recon_picture32bit_ptr = (EbPictureBufferDesc_t *)EB_NULL;
    return_error = eb_recon_picture_buffer_desc_ctor(
        (EbPtr*)&(object_ptr->recon_picture32bit_ptr),
        (EbPtr)&coeffBufferDes32bitInitData);

    // Reconstructed Picture Buffer
    if (initDataPtr->is16bit == EB_TRUE) {
        return_error = eb_recon_picture_buffer_desc_ctor(
            (EbPtr*) &(object_ptr->recon_picture16bit_ptr),
            (EbPtr)&coeffBufferDescInitData);
    }
    else
    {

        return_error = eb_recon_picture_buffer_desc_ctor(
            (EbPtr*) &(object_ptr->recon_picture_ptr),
            (EbPtr)&input_picture_buffer_desc_init_data);
    }

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }


    // Film Grain Picture Buffer
    if (initDataPtr->film_grain_noise_level) {
        if (initDataPtr->is16bit == EB_TRUE) {
            return_error = eb_recon_picture_buffer_desc_ctor(
                (EbPtr*) &(object_ptr->film_grain_picture16bit_ptr),
                (EbPtr)&coeffBufferDescInitData);
        }
        else
        {
            return_error = eb_recon_picture_buffer_desc_ctor(
                (EbPtr*) &(object_ptr->film_grain_picture_ptr),
                (EbPtr)&input_picture_buffer_desc_init_data);
        }

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    if (initDataPtr->is16bit == EB_TRUE) {
        return_error = eb_picture_buffer_desc_ctor(
            (EbPtr*)&(object_ptr->input_frame16bit),
            (EbPtr)&coeffBufferDescInitData);
    }
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }


    // Entropy Coder
    return_error = EntropyCoderCtor(
        &object_ptr->entropy_coder_ptr,
        SEGMENT_ENTROPY_BUFFER_SIZE);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    // Packetization process Bitstream
    return_error = BitstreamCtor(
        &object_ptr->bitstreamPtr,
        PACKETIZATION_PROCESS_BUFFER_SIZE);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Rate estimation entropy coder
    return_error = EntropyCoderCtor(
        &object_ptr->coeff_est_entropy_coder_ptr,
        SEGMENT_ENTROPY_BUFFER_SIZE);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // GOP
    object_ptr->picture_number = 0;
    object_ptr->temporal_layer_index = 0;

    // SB Array
    object_ptr->sb_max_depth = (uint8_t)initDataPtr->max_depth;
    object_ptr->sb_total_count = pictureLcuWidth * pictureLcuHeight;
    EB_MALLOC(LargestCodingUnit_t**, object_ptr->sb_ptr_array, sizeof(LargestCodingUnit_t*) * object_ptr->sb_total_count, EB_N_PTR);

    sb_origin_x = 0;
    sb_origin_y = 0;

    const uint16_t picture_sb_w   = (uint16_t)((initDataPtr->picture_width  + initDataPtr->sb_size_pix - 1) / initDataPtr->sb_size_pix); 
    const uint16_t picture_sb_h   = (uint16_t)((initDataPtr->picture_height + initDataPtr->sb_size_pix - 1) / initDataPtr->sb_size_pix);
    const uint16_t all_sb = picture_sb_w * picture_sb_h;

    for (sb_index = 0; sb_index < all_sb; ++sb_index) {

        return_error = largest_coding_unit_ctor(
            &(object_ptr->sb_ptr_array[sb_index]),
            (uint8_t)initDataPtr->sb_size_pix,
            (uint16_t)(sb_origin_x * maxCuSize),
            (uint16_t)(sb_origin_y * maxCuSize),
            (uint16_t)sb_index,
            object_ptr);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        // Increment the Order in coding order (Raster Scan Order)
        sb_origin_y = (sb_origin_x == picture_sb_w - 1) ? sb_origin_y + 1 : sb_origin_y;
        sb_origin_x = (sb_origin_x == picture_sb_w - 1) ? 0 : sb_origin_x + 1;

    }

    // Copy SB array map
    EB_MALLOC(LargestCodingUnit_t**, object_ptr->sb_ptr_array_copy, sizeof(LargestCodingUnit_t*) * object_ptr->sb_total_count, EB_N_PTR);

    EB_MEMCPY(object_ptr->sb_ptr_array_copy, object_ptr->sb_ptr_array, object_ptr->sb_total_count * sizeof(sizeof(LargestCodingUnit_t*)));

    // Mode Decision Control config
    EB_MALLOC(MdcLcuData_t*, object_ptr->mdc_sb_array, object_ptr->sb_total_count * sizeof(MdcLcuData_t), EB_N_PTR);
    object_ptr->qp_array_stride = (uint16_t)((initDataPtr->picture_width + MIN_BLOCK_SIZE - 1) / MIN_BLOCK_SIZE);
    object_ptr->qp_array_size = ((initDataPtr->picture_width + MIN_BLOCK_SIZE - 1) / MIN_BLOCK_SIZE) *
        ((initDataPtr->picture_height + MIN_BLOCK_SIZE - 1) / MIN_BLOCK_SIZE);


    // Allocate memory for qp array (used by DLF)
    EB_MALLOC(uint8_t*, object_ptr->qp_array, sizeof(uint8_t) * object_ptr->qp_array_size, EB_N_PTR);

    EB_MALLOC(uint8_t*, object_ptr->entropy_qp_array, sizeof(uint8_t) * object_ptr->qp_array_size, EB_N_PTR);

    // Allocate memory for cbf array (used by DLF)
    EB_MALLOC(uint8_t*, object_ptr->cbf_map_array, sizeof(uint8_t) * ((initDataPtr->picture_width >> 2) * (initDataPtr->picture_height >> 2)), EB_N_PTR);

    // Mode Decision Neighbor Arrays
    uint8_t depth;
    for (depth = 0; depth < NEIGHBOR_ARRAY_TOTAL_COUNT; depth++) {
        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_intra_luma_mode_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_intra_chroma_mode_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE >> 1,
            MAX_PICTURE_HEIGHT_SIZE >> 1,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_mv_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(MvUnit_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_skip_flag_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_mode_type_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_leaf_depth_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

        return_error = neighbor_array_unit_ctor(
            &object_ptr->mdleaf_partition_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(struct PartitionContext),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_luma_recon_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_cb_recon_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE >> 1,
            MAX_PICTURE_HEIGHT_SIZE >> 1,
            sizeof(uint8_t),
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_cr_recon_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE >> 1,
            MAX_PICTURE_HEIGHT_SIZE >> 1,
            sizeof(uint8_t),
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }


        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_skip_coeff_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        // for each 4x4
        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        // for each 4x4
        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        // for each 4x4
        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_cb_dc_sign_level_coeff_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_inter_pred_dir_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }

        return_error = neighbor_array_unit_ctor(
            &object_ptr->md_ref_frame_type_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint8_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);


        return_error = neighbor_array_unit_ctor32(
            &object_ptr->md_interpolation_type_neighbor_array[depth],
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint32_t),
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            PU_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);



        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }



    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->md_refinement_intra_luma_mode_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->md_refinement_mode_type_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->md_refinement_luma_recon_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    // Encode Pass Neighbor Arrays
    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_intra_luma_mode_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    // Encode Pass Neighbor Arrays
    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_intra_chroma_mode_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE >> 1,
        MAX_PICTURE_HEIGHT_SIZE >> 1,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_mv_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(MvUnit_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_skip_flag_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        CU_NEIGHBOR_ARRAY_GRANULARITY,
        CU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_mode_type_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_leaf_depth_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,

        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_luma_recon_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_cb_recon_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE >> 1,
        MAX_PICTURE_HEIGHT_SIZE >> 1,
        sizeof(uint8_t),
        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = neighbor_array_unit_ctor(
        &object_ptr->ep_cr_recon_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE >> 1,
        MAX_PICTURE_HEIGHT_SIZE >> 1,
        sizeof(uint8_t),
        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
        SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    if (is16bit) {
        return_error = neighbor_array_unit_ctor(
            &object_ptr->ep_luma_recon_neighbor_array16bit,
            MAX_PICTURE_WIDTH_SIZE,
            MAX_PICTURE_HEIGHT_SIZE,
            sizeof(uint16_t),
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        return_error = neighbor_array_unit_ctor(
            &object_ptr->ep_cb_recon_neighbor_array16bit,
            MAX_PICTURE_WIDTH_SIZE >> 1,
            MAX_PICTURE_HEIGHT_SIZE >> 1,
            sizeof(uint16_t),
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
        return_error = neighbor_array_unit_ctor(
            &object_ptr->ep_cr_recon_neighbor_array16bit,
            MAX_PICTURE_WIDTH_SIZE >> 1,
            MAX_PICTURE_HEIGHT_SIZE >> 1,
            sizeof(uint16_t),
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            SAMPLE_NEIGHBOR_ARRAY_GRANULARITY,
            NEIGHBOR_ARRAY_UNIT_FULL_MASK);
        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }
    else {
        object_ptr->ep_luma_recon_neighbor_array16bit = 0;
        object_ptr->ep_cb_recon_neighbor_array16bit = 0;
        object_ptr->ep_cr_recon_neighbor_array16bit = 0;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->amvp_mv_merge_mv_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(MvUnit_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return_error = neighbor_array_unit_ctor(
        &object_ptr->amvp_mv_merge_mode_type_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_FULL_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Entropy Coding Neighbor Arrays
    return_error = neighbor_array_unit_ctor(
        &object_ptr->mode_type_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->partition_context_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(struct PartitionContext),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->skip_flag_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->skip_coeff_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // for each 4x4
    return_error = neighbor_array_unit_ctor(
        &object_ptr->luma_dc_sign_level_coeff_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // for each 4x4
    return_error = neighbor_array_unit_ctor(
        &object_ptr->cr_dc_sign_level_coeff_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // for each 4x4
    return_error = neighbor_array_unit_ctor(
        &object_ptr->cb_dc_sign_level_coeff_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->inter_pred_dir_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    return_error = neighbor_array_unit_ctor(
        &object_ptr->ref_frame_type_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    return_error = neighbor_array_unit_ctor32(
        &object_ptr->interpolation_type_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint32_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }



    return_error = neighbor_array_unit_ctor(
        &object_ptr->intra_luma_mode_neighbor_array,
        MAX_PICTURE_WIDTH_SIZE,
        MAX_PICTURE_HEIGHT_SIZE,
        sizeof(uint8_t),
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        PU_NEIGHBOR_ARRAY_GRANULARITY,
        NEIGHBOR_ARRAY_UNIT_TOP_AND_LEFT_ONLY_MASK);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    // Note - non-zero offsets are not supported (to be fixed later in DLF chroma filtering)
    object_ptr->cb_qp_offset = 0;
    object_ptr->cr_qp_offset = 0;

    object_ptr->slice_level_chroma_qp_flag = EB_TRUE;
    // slice level chroma QP offsets
    object_ptr->slice_cb_qp_offset = 0;
    object_ptr->slice_cr_qp_offset = 0;


    //object_ptr->total_num_bits = 0;

    // Error Resilience
    object_ptr->constrained_intra_flag = EB_FALSE;

    // Segments
    return_error = EncDecSegmentsCtor(
        &object_ptr->enc_dec_segment_ctrl,
        initDataPtr->enc_dec_segment_col,
        initDataPtr->enc_dec_segment_row);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    // Entropy Rows
    EB_CREATEMUTEX(EbHandle, object_ptr->entropy_coding_mutex, sizeof(EbHandle), EB_MUTEX);

    EB_CREATEMUTEX(EbHandle, object_ptr->intra_mutex, sizeof(EbHandle), EB_MUTEX);

#if CDEF_M
    EB_CREATEMUTEX(EbHandle, object_ptr->cdef_search_mutex, sizeof(EbHandle), EB_MUTEX);

    //object_ptr->mse_seg[0] = (uint64_t(*)[64])aom_malloc(sizeof(**object_ptr->mse_seg) *  pictureLcuWidth * pictureLcuHeight);
   // object_ptr->mse_seg[1] = (uint64_t(*)[64])aom_malloc(sizeof(**object_ptr->mse_seg) *  pictureLcuWidth * pictureLcuHeight);
   
    EB_MALLOC(uint64_t(*)[64], object_ptr->mse_seg[0], sizeof(**object_ptr->mse_seg) *  pictureLcuWidth * pictureLcuHeight, EB_N_PTR);
    EB_MALLOC(uint64_t(*)[64], object_ptr->mse_seg[1], sizeof(**object_ptr->mse_seg) *  pictureLcuWidth * pictureLcuHeight, EB_N_PTR);

    if (is16bit == 0)
    {
        EB_MALLOC(uint16_t*, object_ptr->src[0],sizeof(*object_ptr->src)       * initDataPtr->picture_width * initDataPtr->picture_height,EB_N_PTR);
        EB_MALLOC(uint16_t*, object_ptr->ref_coeff[0],sizeof(*object_ptr->ref_coeff) * initDataPtr->picture_width * initDataPtr->picture_height, EB_N_PTR);
        EB_MALLOC(uint16_t*, object_ptr->src[1],sizeof(*object_ptr->src)       * initDataPtr->picture_width * initDataPtr->picture_height * 3 / 2, EB_N_PTR);
        EB_MALLOC(uint16_t*, object_ptr->ref_coeff[1],sizeof(*object_ptr->ref_coeff) * initDataPtr->picture_width * initDataPtr->picture_height * 3 / 2, EB_N_PTR);
        EB_MALLOC(uint16_t*, object_ptr->src[2],sizeof(*object_ptr->src)       * initDataPtr->picture_width * initDataPtr->picture_height * 3 / 2, EB_N_PTR);
        EB_MALLOC(uint16_t*,object_ptr->ref_coeff[2],sizeof(*object_ptr->ref_coeff) * initDataPtr->picture_width * initDataPtr->picture_height * 3 / 2, EB_N_PTR);
    }
#endif

#if REST_M
    EB_CREATEMUTEX(EbHandle, object_ptr->rest_search_mutex, sizeof(EbHandle), EB_MUTEX);
     
#endif

    object_ptr->cu32x32_quant_coeff_num_map_array_stride = (uint16_t)((initDataPtr->picture_width + 32 - 1) / 32);
    uint16_t cu32x32QuantCoeffNumMapArraySize = (uint16_t)(((initDataPtr->picture_width + 32 - 1) / 32) * ((initDataPtr->picture_height + 32 - 1) / 32));
    EB_MALLOC(int8_t*, object_ptr->cu32x32_quant_coeff_num_map_array, sizeof(int8_t) * cu32x32QuantCoeffNumMapArraySize, EB_N_PTR);

    //the granularity is 4x4
    EB_MALLOC(ModeInfo**, object_ptr->mi_grid_base, sizeof(ModeInfo*) * object_ptr->sb_total_count*(BLOCK_SIZE_64 / 4)*(BLOCK_SIZE_64 / 4), EB_N_PTR);
    EB_MALLOC(ModeInfo*, object_ptr->mip, sizeof(ModeInfo) * object_ptr->sb_total_count*(BLOCK_SIZE_64 / 4)*(BLOCK_SIZE_64 / 4), EB_N_PTR);

    memset(object_ptr->mip, 0, sizeof(ModeInfo) * object_ptr->sb_total_count*(BLOCK_SIZE_64 / 4)*(BLOCK_SIZE_64 / 4));
    // pictureLcuWidth * pictureLcuHeight

    uint32_t miIdx;
    for (miIdx = 0; miIdx < object_ptr->sb_total_count*(BLOCK_SIZE_64 >> MI_SIZE_LOG2)*(BLOCK_SIZE_64 >> MI_SIZE_LOG2); ++miIdx) {

        object_ptr->mi_grid_base[miIdx] = object_ptr->mip + miIdx;
    }

    object_ptr->mi_stride = pictureLcuWidth * (BLOCK_SIZE_64 / 4);
#if ICOPY
    object_ptr->hash_table.p_lookup_table = NULL;
    av1_hash_table_create(&object_ptr->hash_table);
#endif
    return EB_ErrorNone;
}


EbErrorType picture_parent_control_set_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    PictureParentControlSet_t   *object_ptr;
    PictureControlSetInitData_t *initDataPtr = (PictureControlSetInitData_t*)object_init_data_ptr;
    EbErrorType return_error = EB_ErrorNone;
    const uint16_t pictureLcuWidth = (uint16_t)((initDataPtr->picture_width + initDataPtr->sb_sz - 1) / initDataPtr->sb_sz);
    const uint16_t pictureLcuHeight = (uint16_t)((initDataPtr->picture_height + initDataPtr->sb_sz - 1) / initDataPtr->sb_sz);
    uint16_t sb_index;
    uint32_t regionInPictureWidthIndex;
    uint32_t regionInPictureHeightIndex;

    EB_MALLOC(PictureParentControlSet_t*, object_ptr, sizeof(PictureParentControlSet_t), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)object_ptr;
    object_ptr->sequence_control_set_wrapper_ptr = (EbObjectWrapper_t *)EB_NULL;
    object_ptr->input_picture_wrapper_ptr = (EbObjectWrapper_t *)EB_NULL;
    object_ptr->reference_picture_wrapper_ptr = (EbObjectWrapper_t *)EB_NULL;

    object_ptr->enhanced_picture_ptr = (EbPictureBufferDesc_t *)EB_NULL;

    // GOP
    object_ptr->pred_struct_index = 0;
    object_ptr->picture_number = 0;
    object_ptr->idr_flag = EB_FALSE;
    object_ptr->temporal_layer_index = 0;
    object_ptr->total_num_bits = 0;
    object_ptr->last_idr_picture = 0;
    object_ptr->sb_total_count = pictureLcuWidth * pictureLcuHeight;

    object_ptr->data_ll_head_ptr = (EbLinkedListNode *)EB_NULL;
    object_ptr->app_out_data_ll_head_ptr = (EbLinkedListNode *)EB_NULL;
    EB_MALLOC(uint16_t**, object_ptr->variance, sizeof(uint16_t*) * object_ptr->sb_total_count, EB_N_PTR);
    EB_MALLOC(uint8_t**, object_ptr->yMean, sizeof(uint8_t*) * object_ptr->sb_total_count, EB_N_PTR);
    EB_MALLOC(uint8_t**, object_ptr->cbMean, sizeof(uint8_t*) * object_ptr->sb_total_count, EB_N_PTR);
    EB_MALLOC(uint8_t**, object_ptr->crMean, sizeof(uint8_t*) * object_ptr->sb_total_count, EB_N_PTR);
    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
        EB_MALLOC(uint16_t*, object_ptr->variance[sb_index], sizeof(uint16_t) * MAX_ME_PU_COUNT, EB_N_PTR);
        EB_MALLOC(uint8_t*, object_ptr->yMean[sb_index], sizeof(uint8_t) * MAX_ME_PU_COUNT, EB_N_PTR);
        EB_MALLOC(uint8_t*, object_ptr->cbMean[sb_index], sizeof(uint8_t) * 21, EB_N_PTR);
        EB_MALLOC(uint8_t*, object_ptr->crMean[sb_index], sizeof(uint8_t) * 21, EB_N_PTR);
    }
    // Histograms
    uint32_t videoComponent;

    EB_MALLOC(uint32_t****, object_ptr->picture_histogram, sizeof(uint32_t***) * MAX_NUMBER_OF_REGIONS_IN_WIDTH, EB_N_PTR);

    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < MAX_NUMBER_OF_REGIONS_IN_WIDTH; regionInPictureWidthIndex++) {  // loop over horizontal regions
        EB_MALLOC(uint32_t***, object_ptr->picture_histogram[regionInPictureWidthIndex], sizeof(uint32_t**) * MAX_NUMBER_OF_REGIONS_IN_HEIGHT, EB_N_PTR);
    }

    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < MAX_NUMBER_OF_REGIONS_IN_WIDTH; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; regionInPictureHeightIndex++) { // loop over vertical regions
            EB_MALLOC(uint32_t**, object_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex], sizeof(uint32_t*) * 3, EB_N_PTR);
        }
    }


    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < MAX_NUMBER_OF_REGIONS_IN_WIDTH; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < MAX_NUMBER_OF_REGIONS_IN_HEIGHT; regionInPictureHeightIndex++) { // loop over vertical regions
            for (videoComponent = 0; videoComponent < 3; ++videoComponent) {
                EB_MALLOC(uint32_t*, object_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][videoComponent], sizeof(uint32_t) * HISTOGRAM_NUMBER_OF_BINS, EB_N_PTR);
            }
        }
    }
#if OIS_BASED_INTRA
        EB_MALLOC(ois_sb_results_t**, object_ptr->ois_sb_results, sizeof(ois_sb_results_t*) * object_ptr->sb_total_count, EB_N_PTR);
        
    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {

        EB_MALLOC(ois_sb_results_t*, object_ptr->ois_sb_results[sb_index], sizeof(ois_sb_results_t), EB_N_PTR);

        ois_candidate_t* contigousCand;
        EB_MALLOC(ois_candidate_t*, contigousCand, sizeof(ois_candidate_t) * MAX_OIS_CANDIDATES * CU_MAX_COUNT, EB_N_PTR);

        uint32_t cuIdx;
        for (cuIdx = 0; cuIdx < CU_MAX_COUNT; ++cuIdx) {
            object_ptr->ois_sb_results[sb_index]->ois_candidate_array[cuIdx] = &contigousCand[cuIdx*MAX_OIS_CANDIDATES];
        }
    }
#else
    uint32_t maxOisCand = MAX_OPEN_LOOP_INTRA_CANDIDATES ;

    EB_MALLOC(OisCu32Cu16Results_t**, object_ptr->ois_cu32_cu16_results, sizeof(OisCu32Cu16Results_t*) * object_ptr->sb_total_count, EB_N_PTR);

    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {

        EB_MALLOC(OisCu32Cu16Results_t*, object_ptr->ois_cu32_cu16_results[sb_index], sizeof(OisCu32Cu16Results_t), EB_N_PTR);

        OisCandidate_t* contigousCand;
        EB_MALLOC(OisCandidate_t*, contigousCand, sizeof(OisCandidate_t) * maxOisCand * 21, EB_N_PTR);

        uint32_t cuIdx;
        for (cuIdx = 0; cuIdx < 21; ++cuIdx) {
            object_ptr->ois_cu32_cu16_results[sb_index]->sorted_ois_candidate[cuIdx] = &contigousCand[cuIdx*maxOisCand];
        }
    }

    EB_MALLOC(OisCu8Results_t**, object_ptr->ois_cu8_results, sizeof(OisCu8Results_t*) * object_ptr->sb_total_count, EB_N_PTR);
    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {

        EB_MALLOC(OisCu8Results_t*, object_ptr->ois_cu8_results[sb_index], sizeof(OisCu8Results_t), EB_N_PTR);

        OisCandidate_t* contigousCand;
        EB_MALLOC(OisCandidate_t*, contigousCand, sizeof(OisCandidate_t) * maxOisCand * 64, EB_N_PTR);

        uint32_t cuIdx;
        for (cuIdx = 0; cuIdx < 64; ++cuIdx) {
            object_ptr->ois_cu8_results[sb_index]->sorted_ois_candidate[cuIdx] = &contigousCand[cuIdx*maxOisCand];
        }
    }
#endif
    // Motion Estimation Results
    object_ptr->max_number_of_pus_per_sb = (initDataPtr->ext_block_flag) ? MAX_ME_PU_COUNT : SQUARE_PU_COUNT;
    EB_MALLOC(MeCuResults_t**, object_ptr->me_results, sizeof(MeCuResults_t*) * object_ptr->sb_total_count, EB_N_PTR);

    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
        EB_MALLOC(MeCuResults_t*, object_ptr->me_results[sb_index], sizeof(MeCuResults_t) * MAX_ME_PU_COUNT, EB_N_PTR);
    }
    EB_MALLOC(uint32_t*, object_ptr->rc_me_distortion, sizeof(uint32_t) * object_ptr->sb_total_count, EB_N_PTR);
    // ME and OIS Distortion Histograms
    EB_MALLOC(uint16_t*, object_ptr->me_distortion_histogram, sizeof(uint16_t) * NUMBER_OF_SAD_INTERVALS, EB_N_PTR);
    EB_MALLOC(uint16_t*, object_ptr->ois_distortion_histogram, sizeof(uint16_t) * NUMBER_OF_INTRA_SAD_INTERVALS, EB_N_PTR);
    EB_MALLOC(uint32_t*, object_ptr->intra_sad_interval_index, sizeof(uint32_t) * object_ptr->sb_total_count, EB_N_PTR);
    EB_MALLOC(uint32_t*, object_ptr->inter_sad_interval_index, sizeof(uint32_t) * object_ptr->sb_total_count, EB_N_PTR);

    // Enhance background for base layer frames:  zz SAD array
    EB_MALLOC(uint8_t*, object_ptr->zz_cost_array, sizeof(uint8_t) * object_ptr->sb_total_count * 64, EB_N_PTR);

    // Non moving index array
    EB_MALLOC(uint8_t*, object_ptr->non_moving_index_array, sizeof(uint8_t) * object_ptr->sb_total_count, EB_N_PTR);

    // similar Colocated Lcu array
    EB_MALLOC(EbBool*, object_ptr->similar_colocated_sb_array, sizeof(EbBool) * object_ptr->sb_total_count, EB_N_PTR);

    // similar Colocated Lcu array
    EB_MALLOC(EbBool*, object_ptr->similar_colocated_sb_array_ii, sizeof(EbBool) * object_ptr->sb_total_count, EB_N_PTR);

    // SB noise variance array
    EB_MALLOC(uint8_t*, object_ptr->sb_flat_noise_array, sizeof(uint8_t) * object_ptr->sb_total_count, EB_N_PTR);
    EB_MALLOC(uint64_t*, object_ptr->sb_variance_of_variance_over_time, sizeof(uint64_t) * object_ptr->sb_total_count, EB_N_PTR);
    EB_MALLOC(EbBool*, object_ptr->is_sb_homogeneous_over_time, sizeof(EbBool) * object_ptr->sb_total_count, EB_N_PTR);
    EB_MALLOC(EdgeLcuResults_t*, object_ptr->edge_results_ptr, sizeof(EdgeLcuResults_t) * object_ptr->sb_total_count, EB_N_PTR);

    EB_MALLOC(uint8_t*, object_ptr->sharp_edge_sb_flag, sizeof(EbBool) * object_ptr->sb_total_count, EB_N_PTR);

    EB_MALLOC(uint8_t*, object_ptr->failing_motion_sb_flag, sizeof(EbBool) * object_ptr->sb_total_count, EB_N_PTR);

    EB_MALLOC(EbBool*, object_ptr->uncovered_area_sb_flag, sizeof(EbBool) * object_ptr->sb_total_count, EB_N_PTR);

    EB_MALLOC(EbBool*, object_ptr->sb_homogeneous_area_array, sizeof(EbBool) * object_ptr->sb_total_count, EB_N_PTR);

    EB_MALLOC(uint64_t**, object_ptr->var_of_var32x32_based_sb_array, sizeof(uint64_t*) * object_ptr->sb_total_count, EB_N_PTR);
    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
        EB_MALLOC(uint64_t*, object_ptr->var_of_var32x32_based_sb_array[sb_index], sizeof(uint64_t) * 4, EB_N_PTR);
    }
#if !INTRA_INTER_FAST_LOOP
    EB_MALLOC(uint8_t*, object_ptr->cmplx_status_sb, sizeof(uint8_t) * object_ptr->sb_total_count, EB_N_PTR);
#endif
    EB_MALLOC(EbBool*, object_ptr->sb_isolated_non_homogeneous_area_array, sizeof(EbBool) * object_ptr->sb_total_count, EB_N_PTR);

    EB_MALLOC(uint64_t**, object_ptr->sb_y_src_energy_cu_array, sizeof(uint64_t*) * object_ptr->sb_total_count, EB_N_PTR);
    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
        EB_MALLOC(uint64_t*, object_ptr->sb_y_src_energy_cu_array[sb_index], sizeof(uint64_t) * 5, EB_N_PTR);
    }

    EB_MALLOC(uint64_t**, object_ptr->sb_y_src_mean_cu_array, sizeof(uint64_t*) * object_ptr->sb_total_count, EB_N_PTR);
    for (sb_index = 0; sb_index < object_ptr->sb_total_count; ++sb_index) {
        EB_MALLOC(uint64_t*, object_ptr->sb_y_src_mean_cu_array[sb_index], sizeof(uint64_t) * 5, EB_N_PTR);
    }

    EB_MALLOC(uint8_t*, object_ptr->sb_cmplx_contrast_array, sizeof(uint8_t) * object_ptr->sb_total_count, EB_N_PTR);
    EB_MALLOC(uint8_t*, object_ptr->sb_high_contrast_array_dialated, sizeof(uint8_t) * object_ptr->sb_total_count, EB_N_PTR);

    object_ptr->cu32x32_clean_sparse_coeff_map_array_stride = (uint16_t)((initDataPtr->picture_width + 32 - 1) / 32);
    object_ptr->cu32x32_clean_sparse_coeff_map_array_size = (uint16_t)(((initDataPtr->picture_width + 32 - 1) / 32) * ((initDataPtr->picture_height + 32 - 1) / 32));
    // Allocate memory for 32x32 cu array for clean sparse coeff
    EB_MALLOC(uint8_t*, object_ptr->cu32x32_clean_sparse_coeff_map_array, sizeof(uint8_t) * object_ptr->cu32x32_clean_sparse_coeff_map_array_size, EB_N_PTR);

    EB_MALLOC(SbStat_t*, object_ptr->sb_stat_array, sizeof(SbStat_t) * object_ptr->sb_total_count, EB_N_PTR);

    EB_MALLOC(EbSbComplexityStatus*, object_ptr->complex_sb_array, sizeof(EbSbComplexityStatus) * object_ptr->sb_total_count, EB_N_PTR);

    EB_CREATEMUTEX(EbHandle, object_ptr->rc_distortion_histogram_mutex, sizeof(EbHandle), EB_MUTEX);

#if ADAPTIVE_DEPTH_PARTITIONING
    EB_MALLOC(EB_SB_DEPTH_MODE*, object_ptr->sb_depth_mode_array, sizeof(EB_SB_DEPTH_MODE) * object_ptr->sb_total_count, EB_N_PTR);
#else
    EB_MALLOC(EbLcuDepthMode*, object_ptr->sb_md_mode_array, sizeof(EbLcuDepthMode) * object_ptr->sb_total_count, EB_N_PTR);
#endif

    EB_MALLOC(Av1Common*, object_ptr->av1_cm, sizeof(Av1Common), EB_N_PTR);


    object_ptr->av1_cm->interp_filter = SWITCHABLE;

    object_ptr->av1_cm->mi_stride = pictureLcuWidth * (BLOCK_SIZE_64 / 4);

    object_ptr->av1_cm->p_pcs_ptr = object_ptr;

    EB_MALLOC(Yv12BufferConfig*, object_ptr->av1_cm->frame_to_show, sizeof(Yv12BufferConfig), EB_N_PTR);

    object_ptr->av1_cm->use_highbitdepth = initDataPtr->is16bit;
    object_ptr->av1_cm->bit_depth = initDataPtr->is16bit ? EB_10BIT : EB_8BIT;
    object_ptr->av1_cm->subsampling_x = 1;
    object_ptr->av1_cm->subsampling_y = 1;
    object_ptr->av1_cm->width = initDataPtr->picture_width;
    object_ptr->av1_cm->height = initDataPtr->picture_height;
    object_ptr->av1_cm->superres_upscaled_width = initDataPtr->picture_width;
    object_ptr->av1_cm->superres_upscaled_height = initDataPtr->picture_height;

    object_ptr->av1_cm->mi_cols = initDataPtr->picture_width >> MI_SIZE_LOG2;
    object_ptr->av1_cm->mi_rows = initDataPtr->picture_height >> MI_SIZE_LOG2;

    object_ptr->av1_cm->byte_alignment = 0;

    set_restoration_unit_size(initDataPtr->picture_width, initDataPtr->picture_height, 1, 1, object_ptr->av1_cm->rst_info);

    return_error = av1_alloc_restoration_buffers(object_ptr->av1_cm);

    memset(&object_ptr->av1_cm->rst_frame, 0, sizeof(Yv12BufferConfig));

#if REST_M
    int32_t ntiles[2];
    for (int32_t is_uv = 0; is_uv < 2; ++is_uv)
        ntiles[is_uv] = object_ptr->av1_cm->rst_info[is_uv].units_per_tile; //CHKN res_tiles_in_plane

    assert(ntiles[1] <= ntiles[0]);

    EB_MALLOC(RestUnitSearchInfo*, object_ptr->rusi_picture[0], sizeof(RestUnitSearchInfo) * ntiles[0], EB_N_PTR);
    EB_MALLOC(RestUnitSearchInfo*, object_ptr->rusi_picture[1], sizeof(RestUnitSearchInfo) * ntiles[1], EB_N_PTR);
    EB_MALLOC(RestUnitSearchInfo*, object_ptr->rusi_picture[2], sizeof(RestUnitSearchInfo) * ntiles[1], EB_N_PTR);

    //object_ptr->rusi_picture[0] = (RestUnitSearchInfo *)malloc(sizeof(RestUnitSearchInfo) * ntiles[0]);
    //object_ptr->rusi_picture[1] = (RestUnitSearchInfo *)malloc(sizeof(RestUnitSearchInfo) * ntiles[1]);
    //object_ptr->rusi_picture[2] = (RestUnitSearchInfo *)malloc(sizeof(RestUnitSearchInfo) * ntiles[1]);

    memset(object_ptr->rusi_picture[0], 0, sizeof(RestUnitSearchInfo) * ntiles[0]);
    memset(object_ptr->rusi_picture[1], 0, sizeof(RestUnitSearchInfo) * ntiles[1]);
    memset(object_ptr->rusi_picture[2], 0, sizeof(RestUnitSearchInfo) * ntiles[1]);

#endif
    EB_MALLOC(Macroblock*, object_ptr->av1x, sizeof(Macroblock), EB_N_PTR);

    // Film grain noise model if film grain is applied
    if (initDataPtr->film_grain_noise_level) {
        denoise_and_model_init_data_t fg_init_data;
        fg_init_data.encoder_bit_depth = initDataPtr->encoder_bit_depth;
        fg_init_data.noise_level = initDataPtr->film_grain_noise_level;
        fg_init_data.width = initDataPtr->picture_width;
        fg_init_data.height = initDataPtr->picture_height;
        fg_init_data.stride_y = initDataPtr->picture_width + initDataPtr->left_padding + initDataPtr->right_padding;
        fg_init_data.stride_cb = fg_init_data.stride_cr = fg_init_data.stride_y >> 1;

        return_error = denoise_and_model_ctor((EbPtr*)&(object_ptr->denoise_and_model),
            (EbPtr)&fg_init_data);

        if (return_error == EB_ErrorInsufficientResources) {
            return EB_ErrorInsufficientResources;
        }
    }

    return return_error;
}

