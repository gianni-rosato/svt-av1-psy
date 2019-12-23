/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbInitialRateControlResults_h
#define EbInitialRateControlResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"

/**************************************
 * Process Results
 **************************************/
typedef struct InitialRateControlResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
} InitialRateControlResults;

typedef struct InitialRateControlResultInitData {
    int32_t junk;
} InitialRateControlResultInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType initial_rate_control_results_creator(EbPtr *object_dbl_ptr,
                                                        EbPtr  object_init_data_ptr);

#endif //EbInitialRateControlResults_h
