// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureManagerQueue.h"

EbErrorType input_queue_entry_ctor(
    InputQueueEntry      *entryPtr)
{
    (void)entryPtr;
    return EB_ErrorNone;
}

void reference_queue_entry_dctor(EbPtr p)
{
    ReferenceQueueEntry* obj = (ReferenceQueueEntry*)p;
    EB_FREE(obj->list0.list);
    EB_FREE(obj->list1.list);
}

EbErrorType reference_queue_entry_ctor(
    ReferenceQueueEntry  *entryPtr)
{
    entryPtr->dctor = reference_queue_entry_dctor;
    entryPtr->reference_object_ptr = (EbObjectWrapper*)EB_NULL;
    entryPtr->picture_number = ~0u;
    entryPtr->dependent_count = 0;
    entryPtr->reference_available = EB_FALSE;

    EB_MALLOC(entryPtr->list0.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS));

    EB_MALLOC(entryPtr->list1.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS));

    return EB_ErrorNone;
}
// clang-format on
