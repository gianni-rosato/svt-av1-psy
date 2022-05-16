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

#ifndef EbPictureManagerQueue_h
#define EbPictureManagerQueue_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPredictionStructure.h"
#include "EbObject.h"
#include "EbCabacContextModel.h"

#ifdef __cplusplus
extern "C" {
#endif
/************************************************
     * Input Queue Entry
     ************************************************/
struct ReferenceQueueEntry; // empty struct definition

typedef struct InputQueueEntry {
    EbDctor          dctor;
    EbObjectWrapper *input_object_ptr;
    uint32_t         dependent_count;
    ReferenceList   *list0_ptr;
    ReferenceList   *list1_ptr;
    uint32_t         use_count;
    Bool             memory_mgmt_loop_done;
    Bool             rate_control_loop_done;
    Bool             encoding_has_begun;
} InputQueueEntry;

/************************************************
     * Reference Queue Entry
     ************************************************/
typedef struct ReferenceQueueEntry {
    EbDctor dctor;

    uint64_t         picture_number;
    uint64_t         decode_order;
    EbObjectWrapper *reference_object_ptr;
    EbObjectWrapper *ref_wraper;
    uint32_t         dependent_count;
    Bool             release_enable;
    Bool             reference_available;
    uint32_t         dep_list0_count;
    uint32_t         dep_list1_count;
    DependentList    list0;
    DependentList    list1;
    Bool             is_used_as_reference_flag;
    uint64_t         rc_group_index;
    Bool             is_alt_ref;
    Bool             feedback_arrived;
    SliceType        slice_type;
    uint8_t          temporal_layer_index;
    Bool             frame_context_updated;
} ReferenceQueueEntry;

typedef struct PicQueueEntry {
    EbDctor dctor;

    uint64_t pic_num;
    int32_t
        dep_cnt_diff; //increase(e.g 4L->5L) or decrease of dep cnt . not including the run-time decrease
    uint8_t is_done;
} PicQueueEntry;

extern EbErrorType input_queue_entry_ctor(InputQueueEntry *entry_ptr);

extern EbErrorType reference_queue_entry_ctor(ReferenceQueueEntry *entry_ptr);
extern EbErrorType dep_cnt_queue_entry_ctor(PicQueueEntry *entry_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbPictureManagerQueue_h
