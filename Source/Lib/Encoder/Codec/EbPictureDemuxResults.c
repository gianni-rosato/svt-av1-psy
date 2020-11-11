/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
#if FEATURE_INL_ME
EbErrorType picture_manager_result_ctor(PictureManagerResults *object_ptr,
    EbPtr                   object_init_data_ptr) {
    (void)object_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType picture_manager_result_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    PictureManagerResults *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, picture_manager_result_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
#endif
