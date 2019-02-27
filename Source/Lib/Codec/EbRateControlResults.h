/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbRateControlResults_h
#define EbRateControlResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#ifdef __cplusplus
extern "C" {
#endif
/**************************************
 * Process Results
 **************************************/
typedef struct RateControlResults_s {
    EbObjectWrapper_t *pictureControlSetWrapperPtr;
} RateControlResults_t;

typedef struct RateControlResultsInitData_s {
    int32_t junk;
} RateControlResultsInitData_t;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType rate_control_results_ctor(
    EbPtr *object_dbl_ptr, 
    EbPtr object_init_data_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbRateControlResults_h