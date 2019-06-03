/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbMotionEstimationResults.h"

EbErrorType motion_estimation_results_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    MotionEstimationResults *context_ptr;
    EB_MALLOC(MotionEstimationResults*, context_ptr, sizeof(MotionEstimationResults), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)context_ptr;
    object_init_data_ptr = 0;
    //ASSERT(object_init_data_ptr == 0);

    (void)(object_init_data_ptr);
    return EB_ErrorNone;
}
