/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureDecisionReorderQueue.h"

EbErrorType picture_decision_reorder_entry_ctor(
    PictureDecisionReorderEntry_t   **entry_dbl_ptr,
    uint32_t                            picture_number)
{
    EB_MALLOC(PictureDecisionReorderEntry_t*, *entry_dbl_ptr, sizeof(PictureDecisionReorderEntry_t), EB_N_PTR);

    (*entry_dbl_ptr)->picture_number = picture_number;
    (*entry_dbl_ptr)->parentPcsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;

    return EB_ErrorNone;
}
