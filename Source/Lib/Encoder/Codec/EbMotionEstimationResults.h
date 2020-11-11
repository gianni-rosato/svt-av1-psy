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
#if FEATURE_INL_ME
extern EbErrorType ime_results_creator(EbPtr *object_dbl_ptr,
    EbPtr  object_init_data_ptr);
#endif

#ifdef __cplusplus
}
#endif
#endif // EbMotionEstimationResults_h
