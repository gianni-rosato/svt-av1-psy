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
    LargestCodingUnit        **larget_coding_unit_dbl_ptr,
    uint8_t                        sb_size_pix,
    uint16_t                       sb_origin_x,
    uint16_t                       sb_origin_y,
    uint16_t                       sb_index,
    PictureControlSet  *picture_control_set)

{
    EbErrorType return_error = EB_ErrorNone;
    uint32_t tu_index;
    EbPictureBufferDescInitData coeffInitData;

    LargestCodingUnit *largestCodingUnitPtr;
    EB_MALLOC(LargestCodingUnit*, largestCodingUnitPtr, sizeof(LargestCodingUnit), EB_N_PTR);

    *larget_coding_unit_dbl_ptr = largestCodingUnitPtr;

    // ************ SB ***************
        // Which borderLargestCuSize is not a power of two

        // Which borderLargestCuSize is not a power of two
    largestCodingUnitPtr->picture_control_set_ptr = picture_control_set;

    largestCodingUnitPtr->origin_x = sb_origin_x;
    largestCodingUnitPtr->origin_y = sb_origin_y;

    largestCodingUnitPtr->index = sb_index;

    uint32_t cu_i;
    uint32_t  tot_cu_num = sb_size_pix == 128 ? 1024 : 256;

    EB_MALLOC(CodingUnit*, largestCodingUnitPtr->final_cu_arr, sizeof(CodingUnit) * tot_cu_num, EB_N_PTR);

    for (cu_i = 0; cu_i < tot_cu_num; ++cu_i) {
        for (tu_index = 0; tu_index < TRANSFORM_UNIT_MAX_COUNT; ++tu_index)
            largestCodingUnitPtr->final_cu_arr[cu_i].transform_unit_array[tu_index].tu_index = tu_index;
        largestCodingUnitPtr->final_cu_arr[cu_i].leaf_index = cu_i;

        EB_MALLOC(MacroBlockD*, largestCodingUnitPtr->final_cu_arr[cu_i].av1xd, sizeof(MacroBlockD), EB_N_PTR);
    }

    uint32_t  max_block_count = sb_size_pix == 128 ? BLOCK_MAX_COUNT_SB_128 : BLOCK_MAX_COUNT_SB_64;

    EB_MALLOC(PartitionType*, largestCodingUnitPtr->cu_partition_array, sizeof(PartitionType) * max_block_count, EB_N_PTR);

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

    return_error = eb_picture_buffer_desc_ctor(
        (EbPtr*) &(largestCodingUnitPtr->quantized_coeff),
        (EbPtr)&coeffInitData);

    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
    return EB_ErrorNone;
}
