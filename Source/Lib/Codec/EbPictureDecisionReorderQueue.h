/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureDecisionReorderQueue_h
#define EbPictureDecisionReorderQueue_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"

/************************************************
 * Packetization Reorder Queue Entry
 ************************************************/
typedef struct PictureDecisionReorderEntry_s {
    uint64_t                                picture_number;
    EbObjectWrapper_t                    *parentPcsWrapperPtr;
} PictureDecisionReorderEntry_t;

extern EbErrorType PictureDecisionReorderEntryCtor(
    PictureDecisionReorderEntry_t       **entryDblPtr,
    uint32_t                                picture_number);

#endif //EbPictureDecisionReorderQueue_h
