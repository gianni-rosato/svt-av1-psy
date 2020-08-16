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
