/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbEncDecTasks.h"

EbErrorType enc_dec_tasks_ctor(EncDecTasks *context_ptr, EbPtr object_init_data_ptr) {
    (void)context_ptr;
    (void)object_init_data_ptr;

    return EB_ErrorNone;
}
EbErrorType enc_dec_tasks_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    EncDecTasks *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, enc_dec_tasks_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;

    return EB_ErrorNone;
}
