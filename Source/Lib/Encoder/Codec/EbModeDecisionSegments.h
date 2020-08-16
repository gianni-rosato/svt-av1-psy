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

#ifndef EbModeDecisionSegments_h
#define EbModeDecisionSegments_h

#include "EbDefinitions.h"

/**************************************
 * Mode Decision Segments
 **************************************/
typedef struct {
    uint64_t completion_mask;
    EbHandle write_lock_mutex;

    uint32_t total_count;
    uint32_t column_count;
    uint32_t row_count;

    EbBool   in_progress;
    uint32_t current_row_idx;
} MdSegments_t;
#endif // EbModeDecisionSegments_h
