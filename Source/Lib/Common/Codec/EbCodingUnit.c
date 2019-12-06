/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <stdint.h>

#include "EbCodingUnit.h"
#include "EbDefinitions.h"
#include "EbTransformUnit.h"
#include "EbPictureControlSet.h"

void largest_coding_unit_dctor(EbPtr p)
{
    SuperBlock* obj = (SuperBlock*)p;
    EB_DELETE(obj->quantized_coeff);
    EB_FREE_ARRAY(obj->av1xd);
    EB_FREE_ARRAY(obj->final_cu_arr);
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
EbErrorType largest_coding_unit_ctor(
    SuperBlock                     *larget_coding_unit_ptr,
    uint8_t                        sb_size_pix,
    uint16_t                       sb_origin_x,
    uint16_t                       sb_origin_y,
    uint16_t                       sb_index,
    PictureControlSet  *picture_control_set)

{
    uint32_t tu_index;
    EbPictureBufferDescInitData coeffInitData;

    larget_coding_unit_ptr->dctor = largest_coding_unit_dctor;

    // ************ SB ***************
        // Which borderLargestCuSize is not a power of two

        // Which borderLargestCuSize is not a power of two
    larget_coding_unit_ptr->picture_control_set_ptr = picture_control_set;

    larget_coding_unit_ptr->origin_x = sb_origin_x;
    larget_coding_unit_ptr->origin_y = sb_origin_y;

    larget_coding_unit_ptr->index = sb_index;

    uint32_t cu_i;
    uint32_t  tot_cu_num = sb_size_pix == 128 ? 1024 : 256;
    larget_coding_unit_ptr->final_cu_count = tot_cu_num;

    EB_MALLOC_ARRAY(larget_coding_unit_ptr->final_cu_arr, tot_cu_num);
    EB_MALLOC_ARRAY(larget_coding_unit_ptr->av1xd, tot_cu_num);

    for (cu_i = 0; cu_i < tot_cu_num; ++cu_i) {
        for (tu_index = 0; tu_index < TRANSFORM_UNIT_MAX_COUNT; ++tu_index)
            larget_coding_unit_ptr->final_cu_arr[cu_i].transform_unit_array[tu_index].tu_index = tu_index;
        larget_coding_unit_ptr->final_cu_arr[cu_i].leaf_index = cu_i;
        larget_coding_unit_ptr->final_cu_arr[cu_i].av1xd = larget_coding_unit_ptr->av1xd + cu_i;
    }

    uint32_t  max_block_count = sb_size_pix == 128 ? BLOCK_MAX_COUNT_SB_128 : BLOCK_MAX_COUNT_SB_64;

    EB_MALLOC_ARRAY(larget_coding_unit_ptr->cu_partition_array, max_block_count);

    coeffInitData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    coeffInitData.max_width = SB_STRIDE_Y;
    coeffInitData.max_height = SB_STRIDE_Y;
    coeffInitData.bit_depth = EB_32BIT;
    coeffInitData.color_format = picture_control_set->color_format;
    coeffInitData.left_padding = 0;
    coeffInitData.right_padding = 0;
    coeffInitData.top_padding = 0;
    coeffInitData.bot_padding = 0;
    coeffInitData.split_mode = EB_FALSE;

    EB_NEW(
        larget_coding_unit_ptr->quantized_coeff,
        eb_picture_buffer_desc_ctor,
        (EbPtr)&coeffInitData);

    return EB_ErrorNone;
}
