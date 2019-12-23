/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMotionEstimationResults_h
#define EbMotionEstimationResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
/**************************************
     * Process Results
     **************************************/
typedef struct MotionEstimationResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
    uint32_t         segment_index;
} MotionEstimationResults;

typedef struct MotionEstimationResultsInitData {
    int32_t junk;
} MotionEstimationResultsInitData;

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType motion_estimation_results_creator(EbPtr *object_dbl_ptr,
                                                     EbPtr  object_init_data_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbMotionEstimationResults_h
