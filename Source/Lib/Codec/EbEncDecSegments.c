/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "EbEncDecSegments.h"
#include "EbThreads.h"

EbErrorType EncDecSegmentsCtor(
    EncDecSegments_t **segmentsDblPtr,
    uint32_t             segmentColCount,
    uint32_t             segmentRowCount)
{
    uint32_t rowIndex;
    EncDecSegments_t *segmentsPtr;
    EB_MALLOC(EncDecSegments_t*, segmentsPtr, sizeof(EncDecSegments_t), EB_N_PTR);

    *segmentsDblPtr = segmentsPtr;

    segmentsPtr->segmentMaxRowCount = segmentRowCount;
    segmentsPtr->segmentMaxBandCount = segmentRowCount + segmentColCount;
    segmentsPtr->segmentMaxTotalCount = segmentsPtr->segmentMaxRowCount * segmentsPtr->segmentMaxBandCount;

    // Start Arrays
    EB_MALLOC(uint16_t*, segmentsPtr->xStartArray, sizeof(uint16_t) * segmentsPtr->segmentMaxTotalCount, EB_N_PTR);

    EB_MALLOC(uint16_t*, segmentsPtr->yStartArray, sizeof(uint16_t) * segmentsPtr->segmentMaxTotalCount, EB_N_PTR);

    EB_MALLOC(uint16_t*, segmentsPtr->validLcuCountArray, sizeof(uint16_t) * segmentsPtr->segmentMaxTotalCount, EB_N_PTR);

    // Dependency map
    EB_MALLOC(uint8_t*, segmentsPtr->depMap.dependencyMap, sizeof(uint8_t) * segmentsPtr->segmentMaxTotalCount, EB_N_PTR);

    EB_CREATEMUTEX(EbHandle, segmentsPtr->depMap.updateMutex, sizeof(EbHandle), EB_MUTEX);

    // Segment rows
    EB_MALLOC(EncDecSegSegmentRow_t*, segmentsPtr->rowArray, sizeof(EncDecSegSegmentRow_t) * segmentsPtr->segmentMaxRowCount, EB_N_PTR)

        for (rowIndex = 0; rowIndex < segmentsPtr->segmentMaxRowCount; ++rowIndex) {
            EB_CREATEMUTEX(EbHandle, segmentsPtr->rowArray[rowIndex].assignmentMutex, sizeof(EbHandle), EB_MUTEX);
        }

    return EB_ErrorNone;
}



void EncDecSegmentsInit(
    EncDecSegments_t *segmentsPtr,
    uint32_t            segColCount,
    uint32_t            segRowCount,
    uint32_t            picWidthLcu,
    uint32_t            picHeightLcu)
{
    unsigned x, y, yLast;
    unsigned rowIndex, bandIndex, segment_index;

    segmentsPtr->lcuRowCount = picHeightLcu;
    segmentsPtr->lcuBandCount = BAND_TOTAL_COUNT(picHeightLcu, picWidthLcu);
    segmentsPtr->segmentRowCount = segRowCount;
    segmentsPtr->segmentBandCount = BAND_TOTAL_COUNT(segRowCount, segColCount);
    segmentsPtr->segmentTotalCount = segmentsPtr->segmentRowCount * segmentsPtr->segmentBandCount;

    //EB_MEMSET(segmentsPtr->inputMap.inputDependencyMap, 0, sizeof(uint16_t) * segmentsPtr->segmentTotalCount);
    EB_MEMSET(segmentsPtr->validLcuCountArray, 0, sizeof(uint16_t) * segmentsPtr->segmentTotalCount);
    EB_MEMSET(segmentsPtr->xStartArray, -1, sizeof(uint16_t) * segmentsPtr->segmentTotalCount);
    EB_MEMSET(segmentsPtr->yStartArray, -1, sizeof(uint16_t) * segmentsPtr->segmentTotalCount);

    // Initialize the per-LCU input availability map & Start Arrays
    for (y = 0; y < picHeightLcu; ++y) {
        for (x = 0; x < picWidthLcu; ++x) {
            bandIndex = BAND_INDEX(x, y, segmentsPtr->segmentBandCount, segmentsPtr->lcuBandCount);
            rowIndex = ROW_INDEX(y, segmentsPtr->segmentRowCount, segmentsPtr->lcuRowCount);
            segment_index = SEGMENT_INDEX(rowIndex, bandIndex, segmentsPtr->segmentBandCount);

            //++segmentsPtr->inputMap.inputDependencyMap[segment_index];
            ++segmentsPtr->validLcuCountArray[segment_index];
            segmentsPtr->xStartArray[segment_index] = (segmentsPtr->xStartArray[segment_index] == (uint16_t)-1) ?
                (uint16_t)x :
                segmentsPtr->xStartArray[segment_index];
            segmentsPtr->yStartArray[segment_index] = (segmentsPtr->yStartArray[segment_index] == (uint16_t)-1) ?
                (uint16_t)y :
                segmentsPtr->yStartArray[segment_index];
        }
    }

    // Initialize the row-based controls
    for (rowIndex = 0; rowIndex < segmentsPtr->segmentRowCount; ++rowIndex) {
        y = ((rowIndex * segmentsPtr->lcuRowCount) + (segmentsPtr->segmentRowCount - 1)) / segmentsPtr->segmentRowCount;
        yLast = ((((rowIndex + 1) * segmentsPtr->lcuRowCount) + (segmentsPtr->segmentRowCount - 1)) / segmentsPtr->segmentRowCount) - 1;
        bandIndex = BAND_INDEX(0, y, segmentsPtr->segmentBandCount, segmentsPtr->lcuBandCount);

        segmentsPtr->rowArray[rowIndex].startingSegIndex = (uint16_t)SEGMENT_INDEX(rowIndex, bandIndex, segmentsPtr->segmentBandCount);
        bandIndex = BAND_INDEX(picWidthLcu - 1, yLast, segmentsPtr->segmentBandCount, segmentsPtr->lcuBandCount);
        segmentsPtr->rowArray[rowIndex].endingSegIndex = (uint16_t)SEGMENT_INDEX(rowIndex, bandIndex, segmentsPtr->segmentBandCount);
        segmentsPtr->rowArray[rowIndex].currentSegIndex = segmentsPtr->rowArray[rowIndex].startingSegIndex;
    }

    // Initialize the per-segment dependency map
    EB_MEMSET(segmentsPtr->depMap.dependencyMap, 0, sizeof(uint8_t) * segmentsPtr->segmentTotalCount);
    for (rowIndex = 0; rowIndex < segmentsPtr->segmentRowCount; ++rowIndex) {
        for (segment_index = segmentsPtr->rowArray[rowIndex].startingSegIndex; segment_index <= segmentsPtr->rowArray[rowIndex].endingSegIndex; ++segment_index) {

            // Check that segment is valid
            if (segmentsPtr->validLcuCountArray[segment_index]) {
                // Right Neighbor
                if (segment_index < segmentsPtr->rowArray[rowIndex].endingSegIndex) {
                    ++segmentsPtr->depMap.dependencyMap[segment_index + 1];
                }
                // Bottom Neighbor
                if (rowIndex < segmentsPtr->segmentRowCount - 1 && segment_index + segmentsPtr->segmentBandCount >= segmentsPtr->rowArray[rowIndex + 1].startingSegIndex) {
                    ++segmentsPtr->depMap.dependencyMap[segment_index + segmentsPtr->segmentBandCount];
                }
            }
        }
    }

    return;
}

