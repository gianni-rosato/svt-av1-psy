/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbPictureAnalysisResults.h"

EbErrorType picture_analysis_result_ctor(PictureAnalysisResults *object_ptr,
                                         EbPtr                   object_init_data_ptr) {
    (void)object_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType picture_analysis_result_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    PictureAnalysisResults *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, picture_analysis_result_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
