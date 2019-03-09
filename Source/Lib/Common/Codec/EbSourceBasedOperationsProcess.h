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

#ifndef EbSourceBasedOperations_h
#define EbSourceBasedOperations_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbNoiseExtractAVX2.h"
/**************************************
 * Context
 **************************************/

typedef struct SourceBasedOperationsContext_s
{
    EbFifo_t  *initial_rate_control_results_input_fifo_ptr;
    EbFifo_t  *picture_demux_results_output_fifo_ptr;

    // Delta QP Map
    int8_t      min_delta_qp;
    uint8_t     max_delta_qp;

    int16_t     min_delta_qp_weight[3][4];
    int16_t     max_delta_qp_weight[3][4];
                
    // Skin     
    uint8_t     grass_percentage_in_picture;

    // Variance
    uint8_t    *sb_high_contrast_array;

    // local zz cost array
    uint32_t    picture_num_grass_sb;
    uint32_t    sb_high_contrast_count;
    uint32_t    complete_sb_count;
    uint32_t    sb_cmplx_contrast_count;
    uint32_t    high_contrast_num;
    uint32_t    high_contrast_num_ii;
    EbBool      high_dist;
    uint32_t    count_of_moving_sbs;
    uint32_t    countOfNonMovingLcus;
    uint64_t    y_non_moving_mean;
    uint64_t    y_moving_mean;
    uint32_t    to_be_intra_coded_probability;
    uint32_t    depth1_block_num;
    uint8_t    *y_mean_ptr;
    uint8_t    *cr_mean_ptr;
    uint8_t    *cb_mean_ptr;

} SourceBasedOperationsContext_t;

/***************************************
 * Extern Function Declaration
 ***************************************/
extern EbErrorType source_based_operations_context_ctor(
    SourceBasedOperationsContext_t **context_dbl_ptr,
    EbFifo_t                        *initial_rate_control_results_input_fifo_ptr,
    EbFifo_t                        *picture_demux_results_output_fifo_ptr,
    SequenceControlSet_t            *sequence_control_set_ptr);

extern void* source_based_operations_kernel(void *input_ptr);


#endif // EbSourceBasedOperations_h