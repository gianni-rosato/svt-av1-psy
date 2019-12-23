/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbRateControlResults_h
#define EbRateControlResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
/**************************************
 * Process Results
 **************************************/
typedef struct RateControlResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
} RateControlResults;

typedef struct RateControlResultsInitData {
    int32_t junk;
} RateControlResultsInitData;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType rate_control_results_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbRateControlResults_h
