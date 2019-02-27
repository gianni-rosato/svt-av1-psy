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
typedef struct PictureDecisionResults_s
{
    EbObjectWrapper_t   *pictureControlSetWrapperPtr;
    uint32_t               segment_index;
} PictureDecisionResults_t;

typedef struct PictureDecisionResultInitData_s {
    int32_t junk;
} PictureDecisionResultInitData_t;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType picture_decision_result_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);


#endif //EbPictureDecisionResults_h