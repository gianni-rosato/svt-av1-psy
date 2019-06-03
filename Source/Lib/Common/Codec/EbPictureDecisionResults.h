/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureDecisionResults_h
#define EbPictureDecisionResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"

/**************************************
 * Process Results
 **************************************/
typedef struct PictureDecisionResults
{
    EbObjectWrapper   *picture_control_set_wrapper_ptr;
    uint32_t               segment_index;
#if ALTREF_FILTERING_SUPPORT
    uint8_t           task_type; //0:ME   1:Temporal Filtering
#endif
} PictureDecisionResults;

typedef struct PictureDecisionResultInitData {
    int32_t junk;
} PictureDecisionResultInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType picture_decision_result_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);

#endif //EbPictureDecisionResults_h
