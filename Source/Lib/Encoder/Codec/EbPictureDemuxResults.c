/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureDemuxResults.h"

EbErrorType picture_results_ctor(PictureDemuxResults *object_ptr, EbPtr object_init_data_ptr) {
    object_ptr->picture_type = EB_PIC_INVALID;

    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType picture_results_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    PictureDemuxResults *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, picture_results_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
