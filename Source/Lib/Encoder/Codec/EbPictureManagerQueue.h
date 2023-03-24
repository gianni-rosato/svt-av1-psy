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
    Bool             release_enable;
    Bool             reference_available;
    Bool             is_ref;
    uint64_t         rc_group_index;
    Bool             is_alt_ref;
    Bool             feedback_arrived;
    SliceType        slice_type;
    uint8_t          temporal_layer_index;
    Bool             frame_context_updated;
    uint8_t          refresh_frame_mask;
    // decode order of the last frame to use the current entry as a reference
    uint64_t dec_order_of_last_ref;
    // True if frame_end_cdf_update_mode is enabled for this frame
    bool frame_end_cdf_update_required;
    bool is_valid;
} ReferenceQueueEntry;

extern EbErrorType svt_aom_input_queue_entry_ctor(InputQueueEntry *entry_ptr);

extern EbErrorType svt_aom_reference_queue_entry_ctor(ReferenceQueueEntry *entry_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbPictureManagerQueue_h
