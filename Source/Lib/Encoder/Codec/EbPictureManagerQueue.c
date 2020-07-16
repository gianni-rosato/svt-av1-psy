/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureManagerQueue.h"

EbErrorType input_queue_entry_ctor(InputQueueEntry *entry_ptr) {
    (void)entry_ptr;
    return EB_ErrorNone;
}

static void reference_queue_entry_dctor(EbPtr p) {
    ReferenceQueueEntry *obj = (ReferenceQueueEntry *)p;
    EB_FREE_ARRAY(obj->list0.list);
    EB_FREE_ARRAY(obj->list1.list);
}

EbErrorType reference_queue_entry_ctor(ReferenceQueueEntry *entry_ptr) {
    entry_ptr->dctor                = reference_queue_entry_dctor;
    entry_ptr->picture_number       = ~0u;

    EB_MALLOC_ARRAY(entry_ptr->list0.list, (1 << MAX_TEMPORAL_LAYERS));
    EB_MALLOC_ARRAY(entry_ptr->list1.list, (1 << MAX_TEMPORAL_LAYERS));

    return EB_ErrorNone;
}
#if DECOUPLE_ME_RES
EbErrorType dep_cnt_queue_entry_ctor(PicQueueEntry *entry_ptr) {

    entry_ptr->pic_num = ~0u;
    entry_ptr->is_done = 1;
    return EB_ErrorNone;
}
#endif
