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
#include "EbPictureDecisionQueue.h"

static void pa_reference_queue_entry_dctor(EbPtr p) {
    PaReferenceQueueEntry* obj = (PaReferenceQueueEntry*)p;
    EB_FREE_ARRAY(obj->list0.list);
    EB_FREE_ARRAY(obj->list1.list);
}

EbErrorType pa_reference_queue_entry_ctor(PaReferenceQueueEntry* entry_ptr) {
    entry_ptr->dctor = pa_reference_queue_entry_dctor;
    EB_MALLOC_ARRAY(entry_ptr->list0.list, (1 << MAX_TEMPORAL_LAYERS));
    EB_MALLOC_ARRAY(entry_ptr->list1.list, (1 << MAX_TEMPORAL_LAYERS));

    return EB_ErrorNone;
}
