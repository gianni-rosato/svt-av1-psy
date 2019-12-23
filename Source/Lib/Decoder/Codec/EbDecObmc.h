/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef EbDecObmc_h
#define EbDecObmc_h

#include "EbDecHandle.h"

#define MAX_OBMC_LEN 32

typedef struct ObmcCtx {
    /*we can reduce memory size by fixing one side to 32, because
    for left or above pred one side max is 32 only*/
    /*Note : If we change one side DIM to 32,make sure strides also changed appropriately*/
    uint8_t tmp_obmc_bufs[MAX_MB_PLANE][MAX_SB_SIZE * MAX_SB_SIZE * 2];

    int dst_stride[MAX_MB_PLANE];

} ObmcCtx;

static const int max_neighbor_obmc[6] = {0, 1, 2, 3, 4, 4};

void dec_build_obmc_inter_predictors_sb(void *pv_dec_mod_ctx, EbDecHandle *dec_handle,
                                        PartitionInfo *pi, int mi_row, int mi_col);

#endif // EbDecObmc_h
