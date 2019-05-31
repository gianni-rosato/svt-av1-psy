/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#include <stdlib.h>

#include "EbInitialRateControlResults.h"

EbErrorType initial_rate_control_results_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    InitialRateControlResults *object_ptr;
    EB_MALLOC(InitialRateControlResults *, object_ptr, sizeof(InitialRateControlResults), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)object_ptr;
    object_init_data_ptr = 0;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}
