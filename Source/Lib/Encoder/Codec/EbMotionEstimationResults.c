/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbMotionEstimationResults.h"

EbErrorType motion_estimation_results_ctor(MotionEstimationResults *context_ptr,
                                           EbPtr                    object_init_data_ptr) {
    (void)context_ptr;
    (void)(object_init_data_ptr);
    return EB_ErrorNone;
}

EbErrorType motion_estimation_results_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    MotionEstimationResults *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, motion_estimation_results_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
