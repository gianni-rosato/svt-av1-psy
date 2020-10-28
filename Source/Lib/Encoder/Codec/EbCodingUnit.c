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
#include <stdint.h>

#include "EbTransformUnit.h"
#include "EbPictureControlSet.h"

void largest_coding_unit_dctor(EbPtr p) {
    SuperBlock *obj = (SuperBlock *)p;
    EB_DELETE(obj->quantized_coeff);
    EB_FREE_ARRAY(obj->av1xd);
    EB_FREE_ARRAY(obj->final_blk_arr);
    EB_FREE_ARRAY(obj->cu_partition_array);
}
/*
Tasks & Questions
    -Need a GetEmptyChain function for testing sub partitions.  Tie it to an Itr?
    -How many empty CUs do we need?  We need to have enough for the max CU count AND
       enough for temp storage when calculating a subdivision.
    -Where do we store temp reconstructed picture data while deciding to subdivide or not?
    -Need a ReconPicture for each candidate.
    -I don't see a way around doing the copies in temp memory and then copying it in...
*/
EbErrorType largest_coding_unit_ctor(SuperBlock *larget_coding_unit_ptr, uint8_t sb_size_pix,
                                     uint16_t sb_origin_x, uint16_t sb_origin_y, uint16_t sb_index,
                                     PictureControlSet *picture_control_set)

{
    EbPictureBufferDescInitData coeff_init_data;

    larget_coding_unit_ptr->dctor = largest_coding_unit_dctor;

    // ************ SB ***************
    // Which borderLargestCuSize is not a power of two

    // Which borderLargestCuSize is not a power of two
    larget_coding_unit_ptr->pcs_ptr = picture_control_set;

    larget_coding_unit_ptr->origin_x = sb_origin_x;
    larget_coding_unit_ptr->origin_y = sb_origin_y;

    larget_coding_unit_ptr->index = sb_index;
    uint32_t tot_blk_num                    = sb_size_pix == 128 ? 1024 : 256;
    EB_MALLOC_ARRAY(larget_coding_unit_ptr->final_blk_arr, tot_blk_num);
    EB_MALLOC_ARRAY(larget_coding_unit_ptr->av1xd, 1);
    // Do NOT initialize the final_blk_arr here
    // Malloc maximum but only initialize it only when actually used.
    // This will help to same actually memory usage
    uint32_t max_block_count = sb_size_pix == 128 ? BLOCK_MAX_COUNT_SB_128 : BLOCK_MAX_COUNT_SB_64;

    EB_MALLOC_ARRAY(larget_coding_unit_ptr->cu_partition_array, max_block_count);

    coeff_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    coeff_init_data.max_width          = sb_size_pix;
    coeff_init_data.max_height         = sb_size_pix;
    coeff_init_data.bit_depth          = EB_32BIT;
    coeff_init_data.color_format       = picture_control_set->color_format;
    coeff_init_data.left_padding       = 0;
    coeff_init_data.right_padding      = 0;
    coeff_init_data.top_padding        = 0;
    coeff_init_data.bot_padding        = 0;
    coeff_init_data.split_mode         = EB_FALSE;

    EB_NEW(larget_coding_unit_ptr->quantized_coeff,
           svt_picture_buffer_desc_ctor,
           (EbPtr)&coeff_init_data);

    return EB_ErrorNone;
}
