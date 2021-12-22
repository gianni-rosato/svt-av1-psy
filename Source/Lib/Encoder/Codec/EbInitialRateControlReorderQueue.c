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
#include "EbInitialRateControlReorderQueue.h"

EbErrorType initial_rate_control_reorder_entry_ctor(InitialRateControlReorderEntry *entry_ptr,
                                                    uint32_t picture_number) {
    entry_ptr->picture_number         = picture_number;
    entry_ptr->parent_pcs_wrapper_ptr = (EbObjectWrapper *)NULL;

    return EB_ErrorNone;
}
