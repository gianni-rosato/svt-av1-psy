/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureDecisionReorderQueue_h
#define EbPictureDecisionReorderQueue_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbObject.h"

/************************************************
 * Packetization Reorder Queue Entry
 ************************************************/
typedef struct PictureDecisionReorderEntry {
    EbDctor          dctor;
    uint64_t         picture_number;
    EbObjectWrapper *parent_pcs_wrapper_ptr;
} PictureDecisionReorderEntry;

extern EbErrorType picture_decision_reorder_entry_ctor(PictureDecisionReorderEntry *entry_ptr,
                                                       uint32_t                     picture_number);

#endif //EbPictureDecisionReorderQueue_h
