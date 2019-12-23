/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureDecisionResults_h
#define EbPictureDecisionResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"

/**************************************
 * Process Results
 **************************************/
typedef struct PictureDecisionResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
    uint32_t         segment_index;
    uint8_t          task_type; //0:ME   1:Temporal Filtering
} PictureDecisionResults;

typedef struct PictureDecisionResultInitData {
    int32_t junk;
} PictureDecisionResultInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType picture_decision_result_creator(EbPtr *object_dbl_ptr,
                                                   EbPtr  object_init_data_ptr);

#endif //EbPictureDecisionResults_h
