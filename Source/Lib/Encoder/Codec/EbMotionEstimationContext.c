/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbMotionEstimationContext.h"
#include "EbUtility.h"

void motion_estimation_pred_unit_ctor(MePredUnit *pu) {
    pu->distortion = 0xFFFFFFFFull;

    pu->prediction_direction = UNI_PRED_LIST_0;
    return;
}

static void me_context_dctor(EbPtr p) {
    MeContext *obj = (MeContext *)p;
    uint32_t   list_index;
    uint32_t   ref_pic_index;
    EB_FREE_ALIGNED_ARRAY(obj->quarter_sb_buffer);

    EB_FREE_ARRAY(obj->mvd_bits_array);

    for (list_index = 0; list_index < MAX_NUM_OF_REF_PIC_LIST; list_index++) {
        for (ref_pic_index = 0; ref_pic_index < MAX_REF_IDX; ref_pic_index++) {
            EB_FREE_ARRAY(obj->pos_b_buffer[list_index][ref_pic_index]);
            EB_FREE_ARRAY(obj->pos_h_buffer[list_index][ref_pic_index]);
            EB_FREE_ARRAY(obj->pos_j_buffer[list_index][ref_pic_index]);
        }
    }

    EB_FREE_ARRAY(obj->one_d_intermediate_results_buf0);
    EB_FREE_ARRAY(obj->one_d_intermediate_results_buf1);
    EB_FREE_ARRAY(obj->me_candidate);
    EB_FREE_ARRAY(obj->avctemp_buffer);
    EB_FREE_ARRAY(obj->p_eight_pos_sad16x16);
    EB_FREE_ALIGNED_ARRAY(obj->sixteenth_sb_buffer);
    EB_FREE_ALIGNED_ARRAY(obj->sb_buffer);
}
EbErrorType me_context_ctor(MeContext *object_ptr, uint16_t max_input_luma_width,
                            uint16_t max_input_luma_height, uint8_t nsq_present, uint8_t mrp_mode) {
    uint32_t list_index;
    uint32_t ref_pic_index;
    uint32_t pu_index;
    uint32_t me_candidate_index;

    object_ptr->dctor = me_context_dctor;

    // Intermediate SB-sized buffer to retain the input samples
    object_ptr->sb_buffer_stride = BLOCK_SIZE_64;
    EB_MALLOC_ALIGNED_ARRAY(object_ptr->sb_buffer, BLOCK_SIZE_64 * object_ptr->sb_buffer_stride);

    object_ptr->quarter_sb_buffer_stride = (BLOCK_SIZE_64 >> 1);
    EB_MALLOC_ALIGNED_ARRAY(object_ptr->quarter_sb_buffer,
                            (BLOCK_SIZE_64 >> 1) * object_ptr->quarter_sb_buffer_stride);

    object_ptr->sixteenth_sb_buffer_stride = (BLOCK_SIZE_64 >> 2);
    EB_MALLOC_ALIGNED_ARRAY(object_ptr->sixteenth_sb_buffer,
                            (BLOCK_SIZE_64 >> 2) * object_ptr->sixteenth_sb_buffer_stride);
    object_ptr->interpolated_stride =
        MIN((uint16_t)MAX_SEARCH_AREA_WIDTH, (uint16_t)(max_input_luma_width + (PAD_VALUE << 1)));
#if !REMOVE_ME_SUBPEL_CODE
#if ME_MEM_OPT2
    uint16_t max_search_area_height = MIN((uint16_t)MAX_SEARCH_AREA_HEIGHT,
        (uint16_t)(max_input_luma_height + (PAD_VALUE << 1)));
#else
    uint16_t max_search_area_height = MIN((uint16_t)MAX_PICTURE_HEIGHT_SIZE,
                                          (uint16_t)(max_input_luma_height + (PAD_VALUE << 1)));
#endif
#endif
    EB_MEMSET(
        object_ptr->sb_buffer, 0, sizeof(uint8_t) * BLOCK_SIZE_64 * object_ptr->sb_buffer_stride);
    EB_MALLOC_ARRAY(object_ptr->mvd_bits_array, NUMBER_OF_MVD_CASES);
    // 15 intermediate buffers to retain the interpolated reference samples

    //      0    1    2    3
    // 0    A    a    b    c
    // 1    d    e    f    g
    // 2    h    i    j    k
    // 3    n    p    q    r

    //                  _____________
    //                 |             |
    // --I samples --> |Interpolation|-- O samples -->
    //                 | ____________|

    // Before Interpolation: 2 x 3
    //   I   I
    //   I   I
    //   I   I

    // After 1-D Horizontal Interpolation: (2 + 1) x 3 - a, b, and c
    // O I O I O
    // O I O I O
    // O I O I O

    // After 1-D Vertical Interpolation: 2 x (3 + 1) - d, h, and n
    //   O   O
    //   I   I
    //   O   O
    //   I   I
    //   O   O
    //   I   I
    //   O   O

    // After 2-D (Horizontal/Vertical) Interpolation: (2 + 1) x (3 + 1) - e, f, g, i, j, k, n, p, q, and r
    // O   O   O
    //   I   I
    // O   O   O
    //   I   I
    // O   O   O
    //   I   I
    // O   O   O

    for (list_index = 0; list_index < MAX_NUM_OF_REF_PIC_LIST; list_index++) {
        for (ref_pic_index = 0; ref_pic_index < MAX_REF_IDX; ref_pic_index++) {
            EB_MALLOC_ARRAY(object_ptr->pos_b_buffer[list_index][ref_pic_index],
                            object_ptr->interpolated_stride * max_search_area_height);
            EB_MALLOC_ARRAY(object_ptr->pos_h_buffer[list_index][ref_pic_index],
                            object_ptr->interpolated_stride * max_search_area_height);
            EB_MALLOC_ARRAY(object_ptr->pos_j_buffer[list_index][ref_pic_index],
                            object_ptr->interpolated_stride * max_search_area_height);
        }
    }

    EB_MALLOC_ARRAY(object_ptr->one_d_intermediate_results_buf0, BLOCK_SIZE_64 * BLOCK_SIZE_64);

    EB_MALLOC_ARRAY(object_ptr->one_d_intermediate_results_buf1, BLOCK_SIZE_64 * BLOCK_SIZE_64);

    EB_MALLOC_ARRAY(object_ptr->me_candidate,
                    ((mrp_mode == 0) ? ME_RES_CAND_MRP_MODE_0 : ME_RES_CAND_MRP_MODE_1));
    for (pu_index = 0; pu_index < (uint32_t)(nsq_present ? MAX_ME_PU_COUNT : SQUARE_PU_COUNT);
         pu_index++) {
        for (me_candidate_index = 0;
             me_candidate_index <
             (uint32_t)((mrp_mode == 0) ? ME_RES_CAND_MRP_MODE_0 : ME_RES_CAND_MRP_MODE_1);
             me_candidate_index++) {
            motion_estimation_pred_unit_ctor(
                &(object_ptr->me_candidate[me_candidate_index]).pu[pu_index]);
        }
    }

    EB_MALLOC_ARRAY(object_ptr->avctemp_buffer,
                    object_ptr->interpolated_stride * max_search_area_height);
    EB_MALLOC_ARRAY(object_ptr->p_eight_pos_sad16x16,
                    8 * 16); //16= 16 16x16 blocks in a SB.       8=8search points

    // Initialize Alt-Ref parameters
    object_ptr->me_alt_ref = EB_FALSE;

    return EB_ErrorNone;
}
