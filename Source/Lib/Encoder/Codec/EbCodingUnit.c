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
#include "EbUtility.h"

void largest_coding_unit_dctor(EbPtr p) {
    SuperBlock *obj = (SuperBlock *)p;
    EB_FREE_ARRAY(obj->av1xd);
    EB_FREE_ARRAY(obj->final_blk_arr);
    EB_FREE_ARRAY(obj->cu_partition_array);
}
/*
* return the NSQ level
  Used by signal_derivation_multi_processes_oq and memory allocation
*/
bool svt_aom_get_disallow_nsq(EncMode enc_mode, bool is_islice) {
    // Set disallow_nsq
    if (enc_mode <= ENC_M4)
        return false;
    else if (enc_mode <= ENC_M5)
        return (is_islice ? false : true);
    else
        return true;
}

/*
* return the 4x4 level
Used by signal_derivation_enc_dec_kernel_oq and memory allocation
*/
bool svt_aom_get_disallow_4x4(EncMode enc_mode, SliceType slice_type) {
    (void)slice_type;
    if (enc_mode <= ENC_M5)
        return false;
    else
        return true;
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
                                     uint8_t enc_mode, uint16_t max_block_cnt,
                                     PictureControlSet *picture_control_set)

{
    larget_coding_unit_ptr->dctor = largest_coding_unit_dctor;

    // ************ SB ***************
    // Which borderLargestCuSize is not a power of two

    // Which borderLargestCuSize is not a power of two
    larget_coding_unit_ptr->pcs_ptr = picture_control_set;

    larget_coding_unit_ptr->origin_x = sb_origin_x;
    larget_coding_unit_ptr->origin_y = sb_origin_y;

    larget_coding_unit_ptr->index = sb_index;
    bool disallow_nsq             = true;
    for (uint8_t is_islice = 0; is_islice <= 1; is_islice++)
        disallow_nsq = MIN(disallow_nsq, svt_aom_get_disallow_nsq(enc_mode, is_islice));
    bool disallow_4x4 = true;
    for (SliceType slice_type = 0; slice_type < IDR_SLICE + 1; slice_type++)
        disallow_4x4 = MIN(disallow_4x4, svt_aom_get_disallow_4x4(enc_mode, slice_type));
    uint32_t tot_blk_num;
    if (sb_size_pix == 128)
        if (disallow_4x4 && disallow_nsq)
            tot_blk_num = 260;
        else if (disallow_4x4)
            tot_blk_num = 512;
        else
            tot_blk_num = 1024;
    else if (disallow_4x4 && disallow_nsq)
        tot_blk_num = 65;
    else if (disallow_4x4)
        tot_blk_num = 128;
    else
        tot_blk_num = 256;
    EB_MALLOC_ARRAY(larget_coding_unit_ptr->final_blk_arr, tot_blk_num);
    EB_MALLOC_ARRAY(larget_coding_unit_ptr->av1xd, 1);
    // Do NOT initialize the final_blk_arr here
    // Malloc maximum but only initialize it only when actually used.
    // This will help to same actually memory usage
    EB_MALLOC_ARRAY(larget_coding_unit_ptr->cu_partition_array, max_block_cnt);
    return EB_ErrorNone;
}
