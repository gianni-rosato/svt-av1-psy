/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbEntropyCodingResults.h"


EbErrorType EntropyCodingResultsCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    EntropyCodingResults_t *context_ptr;
    EB_MALLOC(EntropyCodingResults_t*, context_ptr, sizeof(EntropyCodingResults_t), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)context_ptr;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

