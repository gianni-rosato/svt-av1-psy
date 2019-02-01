/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureDemuxResults.h"

EbErrorType PictureResultsCtor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    PictureDemuxResults_t *objectPtr;
    EB_MALLOC(PictureDemuxResults_t*, objectPtr, sizeof(PictureDemuxResults_t), EB_N_PTR);

    *object_dbl_ptr = objectPtr;

    objectPtr->pictureType = EB_PIC_INVALID;
    objectPtr->pictureControlSetWrapperPtr = 0;
    objectPtr->reference_picture_wrapper_ptr = 0;
    objectPtr->picture_number = 0;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}


