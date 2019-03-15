/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureManagerQueue_h
#define EbPictureManagerQueue_h

#include "EbDefinitions.h"
#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPredictionStructure.h"
#include "EbApiSei.h"

#ifdef __cplusplus
extern "C" {
#endif
    /************************************************
     * Input Queue Entry
     ************************************************/
    struct ReferenceQueueEntry_s;   // empty struct definition

    typedef struct InputQueueEntry_s {
        EbObjectWrapper_t              *inputObjectPtr;
        uint32_t                          dependentCount;
        uint32_t                          referenceEntryIndex;
        ReferenceList_t                *list0Ptr;
        ReferenceList_t                *list1Ptr;
        uint32_t                          useCount;
        EbBool                         memoryMgmtLoopDone;
        EbBool                         rateControlLoopDone;
        EbBool                         encodingHasBegun;

    } InputQueueEntry_t;

    /************************************************
     * Reference Queue Entry
     ************************************************/
    typedef struct ReferenceQueueEntry_s {

        uint64_t                          picture_number;
        uint64_t                          decode_order;
        EbObjectWrapper_t              *referenceObjectPtr;
        uint32_t                          dependentCount;
        EbBool                         releaseEnable;
        EbBool                         referenceAvailable;
        uint32_t                          depList0Count;
        uint32_t                          depList1Count;
        DependentList_t                 list0;
        DependentList_t                 list1;
        EbBool                         is_used_as_reference_flag;

        uint64_t                          rcGroupIndex;
#if BASE_LAYER_REF
        EB_SLICE                        slice_type;
        uint8_t                         temporal_layer_index;
        uint64_t                        last_islice_picture_number;
#endif

    } ReferenceQueueEntry_t;

    /************************************************
     * Rate Control Input Queue Entry
     ************************************************/

    typedef struct RcInputQueueEntry_s {
        uint64_t                          picture_number;
        EbObjectWrapper_t              *inputObjectPtr;

        EbBool                         isPassed;
        EbBool                         releaseEnabled;
        uint64_t                          groupId;
        uint64_t                          gopFirstPoc;
        uint32_t                          gopIndex;


    } RcInputQueueEntry_t;

    /************************************************
     * Rate Control FeedBack  Queue Entry
     ************************************************/
    typedef struct RcFeedbackQueueEntry_s {

        uint64_t                          picture_number;
        EbObjectWrapper_t              *feedbackObjectPtr;

        EbBool                         isAvailable;
        EbBool                         isUpdated;
        EbBool                         releaseEnabled;
        uint64_t                          groupId;
        uint64_t                          gopFirstPoc;
        uint32_t                          gopIndex;

    } RcFeedbackQueueEntry_t;

    extern EbErrorType input_queue_entry_ctor(
        InputQueueEntry_t **entry_dbl_ptr);



    extern EbErrorType reference_queue_entry_ctor(
        ReferenceQueueEntry_t  **entry_dbl_ptr);


#ifdef __cplusplus
}
#endif
#endif // EbPictureManagerQueue_h