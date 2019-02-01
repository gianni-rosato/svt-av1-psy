/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#include <stdlib.h>

#include "EbInitialRateControlResults.h"

EbErrorType InitialRateControlResultsCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    InitialRateControlResults_t *objectPtr;
    EB_MALLOC(InitialRateControlResults_t *, objectPtr, sizeof(InitialRateControlResults_t), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)objectPtr;
    object_init_data_ptr = 0;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

