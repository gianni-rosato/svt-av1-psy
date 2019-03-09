/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMotionEstimationResults_h
#define EbMotionEstimationResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * Process Results
     **************************************/
    typedef struct MotionEstimationResults_s
    {
        EbObjectWrapper_t   *pictureControlSetWrapperPtr;
        uint32_t               segment_index;
    } MotionEstimationResults_t;

    typedef struct MotionEstimationResultsInitData_s
    {
        int32_t junk;
    } MotionEstimationResultsInitData_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType MotionEstimationResultsCtor(
        EbPtr *object_dbl_ptr,
        EbPtr object_init_data_ptr);


#ifdef __cplusplus
}
#endif
#endif // EbMotionEstimationResults_h