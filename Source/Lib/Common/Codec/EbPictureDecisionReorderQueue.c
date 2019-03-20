/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureManagerReorderQueue.h"

EbErrorType picture_manager_reorder_entry_ctor(
    PictureManagerReorderEntry_t   **entry_dbl_ptr,
    uint32_t                           picture_number)
{
    EB_MALLOC(PictureManagerReorderEntry_t*, *entry_dbl_ptr, sizeof(PictureManagerReorderEntry_t), EB_N_PTR);

    (*entry_dbl_ptr)->picture_number = picture_number;
    (*entry_dbl_ptr)->parentPcsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;

    return EB_ErrorNone;
}
