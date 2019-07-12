/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEncDecSegments_h
#define EbEncDecSegments_h

#include "EbDefinitions.h"
#include "EbThreads.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * Defines
     **************************************/
#define ENCDEC_SEGMENTS_MAX_COL_COUNT  60
#define ENCDEC_SEGMENTS_MAX_ROW_COUNT  37
#define ENCDEC_SEGMENTS_MAX_BAND_COUNT ENCDEC_SEGMENTS_MAX_COL_COUNT + ENCDEC_SEGMENTS_MAX_ROW_COUNT
#define ENCDEC_SEGMENTS_MAX_COUNT      ENCDEC_SEGMENTS_MAX_BAND_COUNT * ENCDEC_SEGMENTS_MAX_ROW_COUNT
#define ENCDEC_SEGMENT_INVALID         0xFFFF

     /**************************************
      * Macros
      **************************************/
#define BAND_TOTAL_COUNT(lcu_row_total_count, lcu_col_total_count) \
    ((lcu_row_total_count) + (lcu_col_total_count) - 1)
#define ROW_INDEX(y_lcu_index, segment_row_count, lcu_row_total_count) \
    (((y_lcu_index) * (segment_row_count)) / (lcu_row_total_count))
#define BAND_INDEX(x_lcu_index, y_lcu_index, segment_band_count, lcu_band_total_count) \
    ((((x_lcu_index) + (y_lcu_index)) * (segment_band_count)) / (lcu_band_total_count))
#define SEGMENT_INDEX(row_index, band_index, segment_band_count) \
    (((row_index) * (segment_band_count)) + (band_index))

      /**************************************
       * Member definitions
       **************************************/
    typedef struct EncDecSegDependencyMap
    {
        uint8_t      *dependency_map;
        EbHandle   update_mutex;
    } EncDecSegDependencyMap;

    typedef struct EncDecSegSegmentRow
    {
        uint16_t      starting_seg_index;
        uint16_t      ending_seg_index;
        uint16_t      current_seg_index;
        EbHandle   assignment_mutex;
    } EncDecSegSegmentRow;

    /**************************************
     * ENCDEC Segments
     **************************************/
    typedef struct EncDecSegments
    {
        EbDctor                 dctor;
        EncDecSegDependencyMap  dep_map;
        EncDecSegSegmentRow    *row_array;

        uint16_t                   *x_start_array;
        uint16_t                   *y_start_array;
        uint16_t                   *valid_lcu_count_array;

        uint32_t                    segment_band_count;
        uint32_t                    segment_row_count;
        uint32_t                    segmentTotalCount;
        uint32_t                    lcu_band_count;
        uint32_t                    lcu_row_count;

        uint32_t                    segment_max_band_count;
        uint32_t                    segment_max_row_count;
        uint32_t                    segment_max_total_count;
    } EncDecSegments;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType enc_dec_segments_ctor(
        EncDecSegments      *segments_dbl_ptr,
        uint32_t             segment_col_count,
        uint32_t             segment_row_count);

    extern void enc_dec_segments_init(
        EncDecSegments *segments_ptr,
        uint32_t            col_count,
        uint32_t            row_count,
        uint32_t            pic_width_lcu,
        uint32_t            pic_height_lcu);
#ifdef __cplusplus
}
#endif
#endif // EbEncDecResults_h
