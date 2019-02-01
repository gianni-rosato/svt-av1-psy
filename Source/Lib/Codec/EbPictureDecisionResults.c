/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbPictureDecisionResults.h"

EbErrorType PictureDecisionResultCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    PictureDecisionResults_t *objectPtr;
    EB_MALLOC(PictureDecisionResults_t *, objectPtr, sizeof(PictureDecisionResults_t), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)objectPtr;
    object_init_data_ptr = 0;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}


