/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbEncDecSegments.h"
#include "EbThreads.h"

EbErrorType enc_dec_segments_ctor(
    EncDecSegments **segments_dbl_ptr,
    uint32_t             segment_col_count,
    uint32_t             segment_row_count)
{
    uint32_t row_index;
    EncDecSegments *segments_ptr;
    EB_MALLOC(EncDecSegments*, segments_ptr, sizeof(EncDecSegments), EB_N_PTR);

    *segments_dbl_ptr = segments_ptr;

    segments_ptr->segment_max_row_count = segment_row_count;
    segments_ptr->segment_max_band_count = segment_row_count + segment_col_count;
    segments_ptr->segment_max_total_count = segments_ptr->segment_max_row_count * segments_ptr->segment_max_band_count;

    // Start Arrays
    EB_MALLOC(uint16_t*, segments_ptr->x_start_array, sizeof(uint16_t) * segments_ptr->segment_max_total_count, EB_N_PTR);

    EB_MALLOC(uint16_t*, segments_ptr->y_start_array, sizeof(uint16_t) * segments_ptr->segment_max_total_count, EB_N_PTR);

    EB_MALLOC(uint16_t*, segments_ptr->valid_lcu_count_array, sizeof(uint16_t) * segments_ptr->segment_max_total_count, EB_N_PTR);

    // Dependency map
    EB_MALLOC(uint8_t*, segments_ptr->dep_map.dependency_map, sizeof(uint8_t) * segments_ptr->segment_max_total_count, EB_N_PTR);

    EB_CREATEMUTEX(EbHandle, segments_ptr->dep_map.update_mutex, sizeof(EbHandle), EB_MUTEX);

    // Segment rows
    EB_MALLOC(EncDecSegSegmentRow*, segments_ptr->row_array, sizeof(EncDecSegSegmentRow) * segments_ptr->segment_max_row_count, EB_N_PTR)

        for (row_index = 0; row_index < segments_ptr->segment_max_row_count; ++row_index) {
            EB_CREATEMUTEX(EbHandle, segments_ptr->row_array[row_index].assignment_mutex, sizeof(EbHandle), EB_MUTEX);
        }

    return EB_ErrorNone;
}

void enc_dec_segments_init(
    EncDecSegments *segments_ptr,
    uint32_t            segColCount,
    uint32_t            segRowCount,
    uint32_t            pic_width_lcu,
    uint32_t            pic_height_lcu)
{
    unsigned x, y, yLast;
    unsigned row_index, band_index, segment_index;

    segments_ptr->lcu_row_count = pic_height_lcu;
    segments_ptr->lcu_band_count = BAND_TOTAL_COUNT(pic_height_lcu, pic_width_lcu);
    segments_ptr->segment_row_count = segRowCount;
    segments_ptr->segment_band_count = BAND_TOTAL_COUNT(segRowCount, segColCount);
    segments_ptr->segmentTotalCount = segments_ptr->segment_row_count * segments_ptr->segment_band_count;

    //EB_MEMSET(segments_ptr->inputMap.inputDependencyMap, 0, sizeof(uint16_t) * segments_ptr->segmentTotalCount);
    EB_MEMSET(segments_ptr->valid_lcu_count_array, 0, sizeof(uint16_t) * segments_ptr->segmentTotalCount);
    EB_MEMSET(segments_ptr->x_start_array, -1, sizeof(uint16_t) * segments_ptr->segmentTotalCount);
    EB_MEMSET(segments_ptr->y_start_array, -1, sizeof(uint16_t) * segments_ptr->segmentTotalCount);

    // Initialize the per-LCU input availability map & Start Arrays
    for (y = 0; y < pic_height_lcu; ++y) {
        for (x = 0; x < pic_width_lcu; ++x) {
            band_index = BAND_INDEX(x, y, segments_ptr->segment_band_count, segments_ptr->lcu_band_count);
            row_index = ROW_INDEX(y, segments_ptr->segment_row_count, segments_ptr->lcu_row_count);
            segment_index = SEGMENT_INDEX(row_index, band_index, segments_ptr->segment_band_count);

            //++segments_ptr->inputMap.inputDependencyMap[segment_index];
            ++segments_ptr->valid_lcu_count_array[segment_index];
            segments_ptr->x_start_array[segment_index] = (segments_ptr->x_start_array[segment_index] == (uint16_t)-1) ?
                (uint16_t)x :
                segments_ptr->x_start_array[segment_index];
            segments_ptr->y_start_array[segment_index] = (segments_ptr->y_start_array[segment_index] == (uint16_t)-1) ?
                (uint16_t)y :
                segments_ptr->y_start_array[segment_index];
        }
    }

    // Initialize the row-based controls
    for (row_index = 0; row_index < segments_ptr->segment_row_count; ++row_index) {
        y = ((row_index * segments_ptr->lcu_row_count) + (segments_ptr->segment_row_count - 1)) / segments_ptr->segment_row_count;
        yLast = ((((row_index + 1) * segments_ptr->lcu_row_count) + (segments_ptr->segment_row_count - 1)) / segments_ptr->segment_row_count) - 1;
        band_index = BAND_INDEX(0, y, segments_ptr->segment_band_count, segments_ptr->lcu_band_count);

        segments_ptr->row_array[row_index].starting_seg_index = (uint16_t)SEGMENT_INDEX(row_index, band_index, segments_ptr->segment_band_count);
        band_index = BAND_INDEX(pic_width_lcu - 1, yLast, segments_ptr->segment_band_count, segments_ptr->lcu_band_count);
        segments_ptr->row_array[row_index].ending_seg_index = (uint16_t)SEGMENT_INDEX(row_index, band_index, segments_ptr->segment_band_count);
        segments_ptr->row_array[row_index].current_seg_index = segments_ptr->row_array[row_index].starting_seg_index;
    }

    // Initialize the per-segment dependency map
    EB_MEMSET(segments_ptr->dep_map.dependency_map, 0, sizeof(uint8_t) * segments_ptr->segmentTotalCount);
    for (row_index = 0; row_index < segments_ptr->segment_row_count; ++row_index) {
        for (segment_index = segments_ptr->row_array[row_index].starting_seg_index; segment_index <= segments_ptr->row_array[row_index].ending_seg_index; ++segment_index) {
            // Check that segment is valid
            if (segments_ptr->valid_lcu_count_array[segment_index]) {
                // Right Neighbor
                if (segment_index < segments_ptr->row_array[row_index].ending_seg_index)
                    ++segments_ptr->dep_map.dependency_map[segment_index + 1];
                // Bottom Neighbor
                if (row_index < segments_ptr->segment_row_count - 1 && segment_index + segments_ptr->segment_band_count >= segments_ptr->row_array[row_index + 1].starting_seg_index)
                    ++segments_ptr->dep_map.dependency_map[segment_index + segments_ptr->segment_band_count];
            }
        }
    }

    return;
}
