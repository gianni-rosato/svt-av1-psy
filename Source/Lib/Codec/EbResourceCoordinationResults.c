/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbResourceCoordinationResults.h"

EbErrorType ResourceCoordinationResultCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    ResourceCoordinationResults_t *objectPtr;
    EB_MALLOC(ResourceCoordinationResults_t*, objectPtr, sizeof(ResourceCoordinationResults_t), EB_N_PTR);

    *object_dbl_ptr = objectPtr;

    object_init_data_ptr = 0;
    (void)object_init_data_ptr;


    return EB_ErrorNone;
}

