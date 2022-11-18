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
#if CLN_PD_REF_Q
/************************************************
 * PA Reference Queue Entry
 ************************************************/
typedef struct PaReferenceEntry {
    EbDctor          dctor;
    EbObjectWrapper *input_object_ptr;
    EbObjectWrapper *eb_y8b_wrapper_ptr;
    uint64_t         picture_number;
#if OPT_REPLACE_DEP_CNT
#if OPT_PD_REF_QUEUE
    /* clang-format off */
    bool is_valid; // The entry will be valid when it represents a valid DPB entry.
                   // This is used in case the DPB is accessed before being populated,
                   // and for when the DPB is cleared at EOS.
    /* clang-format on */
#else
    uint8_t refresh_frame_mask;
#endif
    uint64_t decode_order;
#else
    uint32_t dependent_count;
#endif
#if !OPT_REPLACE_DEP_CNT_CL
    ReferenceList *list0_ptr;
    ReferenceList *list1_ptr;
    uint32_t       dep_list0_count;
    uint32_t       dep_list1_count;
    DependentList  list0;
    DependentList  list1;
#endif
    uint8_t is_alt_ref;
} PaReferenceEntry;

extern EbErrorType pa_reference_queue_entry_ctor(PaReferenceEntry *entry_ptr);
#else
/************************************************
 * PA Reference Queue Entry
 ************************************************/
typedef struct PaReferenceQueueEntry {
    EbDctor          dctor;
    EbObjectWrapper *input_object_ptr;
    EbObjectWrapper *eb_y8b_wrapper_ptr;
    uint64_t         picture_number;
#if OPT_REPLACE_DEP_CNT
#if OPT_PD_REF_QUEUE
    bool     is_valid; // The entry will be valid when it represents a usable picture in the DPB
#else
    uint8_t refresh_frame_mask;
#endif
    uint64_t decode_order;
#else
    uint32_t dependent_count;
#endif
#if !OPT_REPLACE_DEP_CNT_CL
    ReferenceList *list0_ptr;
    ReferenceList *list1_ptr;
    uint32_t       dep_list0_count;
    uint32_t       dep_list1_count;
    DependentList  list0;
    DependentList  list1;
#endif
    uint8_t        is_alt_ref;
} PaReferenceQueueEntry;

extern EbErrorType pa_reference_queue_entry_ctor(PaReferenceQueueEntry *entry_ptr);
#endif

#endif // EbPictureDecisionQueue_h
