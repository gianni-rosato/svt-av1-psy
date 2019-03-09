/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbInitialRateControlResults_h
#define EbInitialRateControlResults_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"

/**************************************
 * Process Results
 **************************************/
typedef struct InitialRateControlResults_s
{
    EbObjectWrapper_t   *pictureControlSetWrapperPtr;
} InitialRateControlResults_t;

typedef struct InitialRateControlResultInitData_s
{
    int32_t junk;
} InitialRateControlResultInitData_t;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType InitialRateControlResultsCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr);


#endif //EbInitialRateControlResults_h