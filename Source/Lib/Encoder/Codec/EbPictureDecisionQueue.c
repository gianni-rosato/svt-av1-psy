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
#if OPT_REPLACE_DEP_CNT_CL
    UNUSED(p);
#else
    PaReferenceQueueEntry* obj = (PaReferenceQueueEntry*)p;
    EB_FREE_ARRAY(obj->list0.list);
    EB_FREE_ARRAY(obj->list1.list);
#endif
}

#if CLN_PD_REF_Q
EbErrorType pa_reference_queue_entry_ctor(PaReferenceEntry* entry_ptr) {
#else
EbErrorType pa_reference_queue_entry_ctor(PaReferenceQueueEntry* entry_ptr) {
#endif
    entry_ptr->dctor = pa_reference_queue_entry_dctor;
#if OPT_PD_REF_QUEUE
    entry_ptr->is_valid = 0;
#endif
#if !OPT_REPLACE_DEP_CNT_CL
    EB_MALLOC_ARRAY(entry_ptr->list0.list, (1 << MAX_TEMPORAL_LAYERS));
    EB_MALLOC_ARRAY(entry_ptr->list1.list, (1 << MAX_TEMPORAL_LAYERS));
#endif

    return EB_ErrorNone;
}
