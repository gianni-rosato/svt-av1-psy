/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureDemuxResults.h"

EbErrorType picture_results_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    PictureDemuxResults *object_ptr;
    EB_MALLOC(PictureDemuxResults*, object_ptr, sizeof(PictureDemuxResults), EB_N_PTR);

    *object_dbl_ptr = object_ptr;

    object_ptr->picture_type = EB_PIC_INVALID;
    object_ptr->picture_control_set_wrapper_ptr = 0;
    object_ptr->reference_picture_wrapper_ptr = 0;
    object_ptr->picture_number = 0;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}
