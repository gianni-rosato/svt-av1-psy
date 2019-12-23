/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbResourceCoordinationResults.h"

EbErrorType resource_coordination_result_ctor(ResourceCoordinationResults *object_ptr,
                                              EbPtr                        object_init_data_ptr) {
    (void)object_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}

EbErrorType resource_coordination_result_creator(EbPtr *object_dbl_ptr,
                                                 EbPtr  object_init_data_ptr) {
    ResourceCoordinationResults *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, resource_coordination_result_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
