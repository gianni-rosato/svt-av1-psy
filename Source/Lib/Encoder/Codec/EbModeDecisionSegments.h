/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
