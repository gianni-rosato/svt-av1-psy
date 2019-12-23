/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbResourceCoordinationResults_h
#define EbResourceCoordinationResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif
/**************************************
     * Process Results
     **************************************/
typedef struct ResourceCoordinationResults {
    EbDctor          dctor;
    EbObjectWrapper *pcs_wrapper_ptr;
} ResourceCoordinationResults;

typedef struct ResourceCoordinationResultInitData {
    int32_t junk;
} ResourceCoordinationResultInitData;

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType resource_coordination_result_creator(EbPtr *object_dbl_ptr,
                                                        EbPtr  object_init_data_ptr);

#ifdef __cplusplus
}
#endif
#endif //EbResourceCoordinationResults_h
