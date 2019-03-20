/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbResourceCoordinationResults_h
#define EbResourceCoordinationResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * Process Results
     **************************************/
    typedef struct ResourceCoordinationResults_s {
        EbObjectWrapper_t *pictureControlSetWrapperPtr;
    } ResourceCoordinationResults_t;

    typedef struct ResourceCoordinationResultInitData_s {
        int32_t junk;
    } ResourceCoordinationResultInitData_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType resource_coordination_result_ctor(
        EbPtr *object_dbl_ptr,
        EbPtr  object_init_data_ptr);


#ifdef __cplusplus
}
#endif
#endif //EbResourceCoordinationResults_h