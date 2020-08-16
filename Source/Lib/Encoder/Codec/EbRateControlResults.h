/*
* Copyright (c) 2019, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
