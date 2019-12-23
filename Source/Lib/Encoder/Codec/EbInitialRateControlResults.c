/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#include <stdlib.h>

#include "EbInitialRateControlResults.h"

EbErrorType initial_rate_control_results_ctor(InitialRateControlResults *object_ptr,
                                              EbPtr                      object_init_data_ptr) {
    (void)object_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType initial_rate_control_results_creator(EbPtr *object_dbl_ptr,
                                                 EbPtr  object_init_data_ptr) {
    InitialRateControlResults *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, initial_rate_control_results_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
