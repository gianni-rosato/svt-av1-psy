/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbMotionEstimationContext.h"


void MotionEstimetionPredUnitCtor(
    MePredUnit_t   *pu)
{

    pu->distortion = 0xFFFFFFFFull;

    pu->prediction_direction = UNI_PRED_LIST_0;

    pu->mv[0] = 0;

    pu->mv[1] = 0;

    return;
}


EbErrorType MeContextCtor(
    MeContext_t     **object_dbl_ptr)
{
    uint32_t                   listIndex;
    uint32_t                   refPicIndex;
    uint32_t                   pu_index;
    uint32_t                   meCandidateIndex;

    EB_MALLOC(MeContext_t*, *object_dbl_ptr, sizeof(MeContext_t), EB_N_PTR);

    // Intermediate LCU-sized buffer to retain the input samples
    (*object_dbl_ptr)->sb_buffer_stride = BLOCK_SIZE_64;
    EB_ALLIGN_MALLOC(uint8_t *, (*object_dbl_ptr)->sb_buffer, sizeof(uint8_t) * BLOCK_SIZE_64 * (*object_dbl_ptr)->sb_buffer_stride, EB_A_PTR);

    (*object_dbl_ptr)->quarter_sb_buffer_stride = (BLOCK_SIZE_64 >> 1);
    EB_MALLOC(uint8_t *, (*object_dbl_ptr)->quarter_sb_buffer, sizeof(uint8_t) * (BLOCK_SIZE_64 >> 1) * (*object_dbl_ptr)->quarter_sb_buffer_stride, EB_N_PTR);

    (*object_dbl_ptr)->sixteenth_sb_buffer_stride = (BLOCK_SIZE_64 >> 2);
    EB_ALLIGN_MALLOC(uint8_t *, (*object_dbl_ptr)->sixteenth_sb_buffer, sizeof(uint8_t) * (BLOCK_SIZE_64 >> 2) * (*object_dbl_ptr)->sixteenth_sb_buffer_stride, EB_A_PTR);

    (*object_dbl_ptr)->interpolated_stride = MAX_SEARCH_AREA_WIDTH;

    EB_MEMSET((*object_dbl_ptr)->sb_buffer, 0, sizeof(uint8_t) * BLOCK_SIZE_64 * (*object_dbl_ptr)->sb_buffer_stride);
    EB_MALLOC(EB_BitFraction *, (*object_dbl_ptr)->mvd_bits_array, sizeof(EB_BitFraction) * NUMBER_OF_MVD_CASES, EB_N_PTR);
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

    for (listIndex = 0; listIndex < MAX_NUM_OF_REF_PIC_LIST; listIndex++) {

        for (refPicIndex = 0; refPicIndex < MAX_REF_IDX; refPicIndex++) {

            EB_MALLOC(uint8_t *, (*object_dbl_ptr)->pos_b_buffer[listIndex][refPicIndex], sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

            EB_MALLOC(uint8_t *, (*object_dbl_ptr)->pos_h_buffer[listIndex][refPicIndex], sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

            EB_MALLOC(uint8_t *, (*object_dbl_ptr)->pos_j_buffer[listIndex][refPicIndex], sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

        }

    }

    EB_MALLOC(EbByte, (*object_dbl_ptr)->one_d_intermediate_results_buf0, sizeof(uint8_t)*BLOCK_SIZE_64*BLOCK_SIZE_64, EB_N_PTR);

    EB_MALLOC(EbByte, (*object_dbl_ptr)->one_d_intermediate_results_buf1, sizeof(uint8_t)*BLOCK_SIZE_64*BLOCK_SIZE_64, EB_N_PTR);

    for (pu_index = 0; pu_index < MAX_ME_PU_COUNT; pu_index++) {
        for (meCandidateIndex = 0; meCandidateIndex < MAX_ME_CANDIDATE_PER_PU; meCandidateIndex++) {
            MotionEstimetionPredUnitCtor(&((*object_dbl_ptr)->me_candidate[meCandidateIndex]).pu[pu_index]);
        }
    }

    EB_MALLOC(uint8_t *, (*object_dbl_ptr)->avctemp_buffer, sizeof(uint8_t) * (*object_dbl_ptr)->interpolated_stride * MAX_SEARCH_AREA_HEIGHT, EB_N_PTR);

    EB_MALLOC(uint16_t *, (*object_dbl_ptr)->p_eight_pos_sad16x16, sizeof(uint16_t) * 8 * 16, EB_N_PTR);//16= 16 16x16 blocks in a LCU.       8=8search points

    return EB_ErrorNone;
}
