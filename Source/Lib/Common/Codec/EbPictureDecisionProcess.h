/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureDecision_h
#define EbPictureDecision_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"

/**************************************
 * Context
 **************************************/
typedef struct PictureDecisionContext_s
{
    EbFifo_t                     *picture_analysis_results_input_fifo_ptr;
    EbFifo_t                     *picture_decision_results_output_fifo_ptr;

    uint64_t                       lastSolidColorFramePoc;

    EbBool                         resetRunningAvg;

    uint32_t **ahdRunningAvgCb;
    uint32_t **ahdRunningAvgCr;
    uint32_t **ahdRunningAvg;
    EbBool        isSceneChangeDetected;

    // Dynamic GOP
    uint32_t        totalRegionActivityCost[MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT];

    uint32_t        totalNumberOfMiniGops;

    uint32_t        miniGopStartIndex[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t        miniGopEndIndex[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t        miniGopLength[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t        miniGopIntraCount[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t        miniGopIdrCount[MINI_GOP_WINDOW_MAX_COUNT];
    uint32_t        miniGopHierarchicalLevels[MINI_GOP_WINDOW_MAX_COUNT];
    EbBool        miniGopActivityArray[MINI_GOP_MAX_COUNT];
    uint32_t        miniGopRegionActivityCostArray[MINI_GOP_MAX_COUNT][MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT];

    uint32_t        miniGopGroupFadedInPicturesCount[MINI_GOP_MAX_COUNT];
    uint32_t        miniGopGroupFadedOutPicturesCount[MINI_GOP_MAX_COUNT];


    EbBool miniGopToggle;    //mini GOP toggling since last Key Frame  K-0-1-0-1-0-K-0-1-0-1-K-0-1.....
} PictureDecisionContext_t;

/***************************************
 * Extern Function Declaration
 ***************************************/
extern EbErrorType picture_decision_context_ctor(
    PictureDecisionContext_t **context_dbl_ptr,
    EbFifo_t                  *picture_analysis_results_input_fifo_ptr,
    EbFifo_t                  *picture_decision_results_output_fifo_ptr);


extern void* picture_decision_kernel(void *input_ptr);

#endif // EbPictureDecision_h