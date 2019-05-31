/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbDefinitions.h"
#include "EbEncDecResults.h"

EbErrorType enc_dec_results_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    EncDecResults *context_ptr;
    EB_MALLOC(EncDecResults*, context_ptr, sizeof(EncDecResults), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)context_ptr;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}
