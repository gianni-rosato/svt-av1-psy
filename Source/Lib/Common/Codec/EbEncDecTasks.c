/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbEncDecTasks.h"

EbErrorType enc_dec_tasks_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    EncDecTasks *context_ptr;
    EB_MALLOC(EncDecTasks*, context_ptr, sizeof(EncDecTasks), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)context_ptr;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}
