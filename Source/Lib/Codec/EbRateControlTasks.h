/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbRateControlTasks_h
#define EbRateControlTasks_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"


/**************************************
 * Tasks Types
 **************************************/
typedef enum RATE_CONTROL_TASKTYPES {
    RC_PICTURE_MANAGER_RESULT,
    RC_PACKETIZATION_FEEDBACK_RESULT,
    RC_ENTROPY_CODING_ROW_FEEDBACK_RESULT,
    RC_INVALID_TASK
} RATE_CONTROL_TASKTYPES;

/**************************************
 * Process Results
 **************************************/
typedef struct RateControlTasks_s
{
    RATE_CONTROL_TASKTYPES              taskType;
    EbObjectWrapper_t                  *pictureControlSetWrapperPtr;
    uint32_t                              segment_index;

    // Following are valid for RC_ENTROPY_CODING_ROW_FEEDBACK_RESULT only
    uint64_t                              picture_number;
    uint32_t                              rowNumber;
    uint32_t                              bitCount;

} RateControlTasks_t;

typedef struct RateControlTasksInitData_s {
    int32_t junk;
} RateControlTasksInitData_t;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType rate_control_tasks_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);


#endif // EbRateControlTasks_h