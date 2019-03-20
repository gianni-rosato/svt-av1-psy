/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEncDecSegments_h
#define EbEncDecSegments_h

#include "EbDefinitions.h"
#include "EbThreads.h"
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
#define BAND_TOTAL_COUNT(lcuRowTotalCount, lcuColTotalCount) \
    ((lcuRowTotalCount) + (lcuColTotalCount) - 1)
#define ROW_INDEX(yLcuIndex, segmentRowCount, lcuRowTotalCount) \
    (((yLcuIndex) * (segmentRowCount)) / (lcuRowTotalCount))
#define BAND_INDEX(xLcuIndex, yLcuIndex, segmentBandCount, lcuBandTotalCount) \
    ((((xLcuIndex) + (yLcuIndex)) * (segmentBandCount)) / (lcuBandTotalCount))
#define SEGMENT_INDEX(row_index, bandIndex, segmentBandCount) \
    (((row_index) * (segmentBandCount)) + (bandIndex))

      /**************************************
       * Member definitions
       **************************************/
    typedef struct {
        uint8_t      *dependencyMap;
        EbHandle   updateMutex;
    } EncDecSegDependencyMap_t;

    typedef struct {
        uint16_t      startingSegIndex;
        uint16_t      endingSegIndex;
        uint16_t      currentSegIndex;
        EbHandle   assignmentMutex;
    } EncDecSegSegmentRow_t;

    /**************************************
     * ENCDEC Segments
     **************************************/
    typedef struct
    {
        EncDecSegDependencyMap_t  depMap;
        EncDecSegSegmentRow_t    *rowArray;

        uint16_t                   *xStartArray;
        uint16_t                   *yStartArray;
        uint16_t                   *validLcuCountArray;

        uint32_t                    segmentBandCount;
        uint32_t                    segmentRowCount;
        uint32_t                    segmentTotalCount;
        uint32_t                    lcuBandCount;
        uint32_t                    lcuRowCount;

        uint32_t                    segmentMaxBandCount;
        uint32_t                    segmentMaxRowCount;
        uint32_t                    segmentMaxTotalCount;

    } EncDecSegments_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType EncDecSegmentsCtor(
        EncDecSegments_t **segmentsDblPtr,
        uint32_t             segmentColCount,
        uint32_t             segmentRowCount);


    extern void EncDecSegmentsInit(
        EncDecSegments_t *segmentsPtr,
        uint32_t            colCount,
        uint32_t            row_count,
        uint32_t            picWidthLcu,
        uint32_t            picHeightLcu);
#ifdef __cplusplus
}
#endif
#endif // EbEncDecResults_h