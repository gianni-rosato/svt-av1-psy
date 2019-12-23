/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbRateControlResults.h"

EbErrorType rate_control_results_ctor(RateControlResults *context_ptr, EbPtr object_init_data_ptr) {
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType rate_control_results_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    RateControlResults *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, rate_control_results_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
