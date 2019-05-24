/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureManagerQueue.h"

EbErrorType input_queue_entry_ctor(
    InputQueueEntry      **entry_dbl_ptr)
{
    InputQueueEntry *entryPtr;
    EB_MALLOC(InputQueueEntry*, entryPtr, sizeof(InputQueueEntry), EB_N_PTR);
    *entry_dbl_ptr = entryPtr;

    entryPtr->input_object_ptr = (EbObjectWrapper*)EB_NULL;
    entryPtr->reference_entry_index = 0;
    entryPtr->dependent_count = 0;
#if BASE_LAYER_REF
    EB_MALLOC(ReferenceList*, entryPtr->list0_ptr, sizeof(ReferenceList), EB_N_PTR);
    EB_MALLOC(ReferenceList*, entryPtr->list1_ptr, sizeof(ReferenceList), EB_N_PTR);
#else
    entryPtr->list0_ptr = (ReferenceList*)EB_NULL;
    entryPtr->list1_ptr = (ReferenceList*)EB_NULL;
#endif

    return EB_ErrorNone;
}



EbErrorType reference_queue_entry_ctor(
    ReferenceQueueEntry  **entry_dbl_ptr)
{
    ReferenceQueueEntry *entryPtr;
    EB_MALLOC(ReferenceQueueEntry*, entryPtr, sizeof(ReferenceQueueEntry), EB_N_PTR);
    *entry_dbl_ptr = entryPtr;

    entryPtr->reference_object_ptr = (EbObjectWrapper*)EB_NULL;
    entryPtr->picture_number = ~0u;
    entryPtr->dependent_count = 0;
    entryPtr->reference_available = EB_FALSE;

    EB_MALLOC(int32_t*, entryPtr->list0.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS), EB_N_PTR);

    EB_MALLOC(int32_t*, entryPtr->list1.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS), EB_N_PTR);

    return EB_ErrorNone;
}
