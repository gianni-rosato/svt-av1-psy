/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureManagerReorderQueue.h"

EbErrorType PictureManagerReorderEntryCtor(
    PictureManagerReorderEntry_t   **entryDblPtr,
    uint32_t                           picture_number)
{
    EB_MALLOC(PictureManagerReorderEntry_t*, *entryDblPtr, sizeof(PictureManagerReorderEntry_t), EB_N_PTR);

    (*entryDblPtr)->picture_number = picture_number;
    (*entryDblPtr)->parentPcsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;

    return EB_ErrorNone;
}
