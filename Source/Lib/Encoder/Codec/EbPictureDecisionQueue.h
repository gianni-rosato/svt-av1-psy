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

#ifndef EbPictureDecisionQueue_h
#define EbPictureDecisionQueue_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPredictionStructure.h"
#include "EbObject.h"
/************************************************
 * PA Reference Queue Entry
 ************************************************/
typedef struct PaReferenceQueueEntry {
    EbDctor          dctor;
    EbObjectWrapper *input_object_ptr;
    uint64_t         picture_number;
    uint32_t         dependent_count;
#if !DECOUPLE_ME_RES
    uint32_t         reference_entry_index;
#endif
    ReferenceList *  list0_ptr;
    ReferenceList *  list1_ptr;
    uint32_t         dep_list0_count;
    uint32_t         dep_list1_count;
    DependentList    list0;
    DependentList    list1;
    uint8_t          is_alt_ref;
} PaReferenceQueueEntry;

extern EbErrorType pa_reference_queue_entry_ctor(PaReferenceQueueEntry *entry_ptr);

#endif // EbPictureDecisionQueue_h
