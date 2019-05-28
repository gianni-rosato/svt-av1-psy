/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPictureDecisionQueue.h"


EbErrorType pa_reference_queue_entry_ctor(
    PaReferenceQueueEntry   **entry_dbl_ptr)
{
    PaReferenceQueueEntry *entryPtr;
    EB_MALLOC(PaReferenceQueueEntry*, entryPtr, sizeof(PaReferenceQueueEntry), EB_N_PTR);
    *entry_dbl_ptr = entryPtr;

    entryPtr->input_object_ptr = (EbObjectWrapper*)EB_NULL;
    entryPtr->picture_number = 0;
    entryPtr->reference_entry_index = 0;
    entryPtr->dependent_count = 0;
#if BASE_LAYER_REF
    EB_MALLOC(ReferenceList*, entryPtr->list0_ptr, sizeof(ReferenceList), EB_N_PTR);
    EB_MALLOC(ReferenceList*, entryPtr->list1_ptr, sizeof(ReferenceList), EB_N_PTR);
    entryPtr->list0_ptr->reference_list = 0;
    entryPtr->list0_ptr->reference_list_count = 0;
    entryPtr->list1_ptr->reference_list = 0;
    entryPtr->list1_ptr->reference_list_count = 0;
#else
    entryPtr->list0_ptr = (ReferenceList*)EB_NULL;
    entryPtr->list1_ptr = (ReferenceList*)EB_NULL;
#endif
    EB_MALLOC(int32_t*, entryPtr->list0.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS), EB_N_PTR);
    EB_MALLOC(int32_t*, entryPtr->list1.list, sizeof(int32_t) * (1 << MAX_TEMPORAL_LAYERS), EB_N_PTR);

    return EB_ErrorNone;
}


