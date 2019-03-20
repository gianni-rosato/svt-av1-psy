/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureDecisionQueue_h
#define EbPictureDecisionQueue_h

#include "EbDefinitions.h"
#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPredictionStructure.h"
#include "EbApiSei.h"
#include"EbPictureControlSet.h"
/************************************************
 * PA Reference Queue Entry
 ************************************************/
typedef struct PaReferenceQueueEntry_s 
{
    EbObjectWrapper_t              *inputObjectPtr;
    uint64_t                          picture_number;
    uint32_t                          dependentCount;
    uint32_t                          referenceEntryIndex;
    ReferenceList_t                *list0Ptr;
    ReferenceList_t                *list1Ptr;
    uint32_t                          depList0Count;
    uint32_t                          depList1Count;
    DependentList_t                 list0;
    DependentList_t                 list1;

    PictureParentControlSet_t       *pPcsPtr;
} PaReferenceQueueEntry_t;

extern EbErrorType pa_reference_queue_entry_ctor(
    PaReferenceQueueEntry_t  **entry_dbl_ptr);


#endif // EbPictureDecisionQueue_h