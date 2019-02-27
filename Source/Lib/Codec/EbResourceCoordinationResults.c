/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbResourceCoordinationResults.h"

EbErrorType resource_coordination_result_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    ResourceCoordinationResults_t *object_ptr;
    EB_MALLOC(ResourceCoordinationResults_t*, object_ptr, sizeof(ResourceCoordinationResults_t), EB_N_PTR);

    *object_dbl_ptr = object_ptr;

    object_init_data_ptr = 0;
    (void)object_init_data_ptr;


    return EB_ErrorNone;
}

