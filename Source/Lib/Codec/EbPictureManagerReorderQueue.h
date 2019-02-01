/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureManagerReorderQueue_h
#define EbPictureManagerReorderQueue_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"

/************************************************
 * Packetization Reorder Queue Entry
 ************************************************/
typedef struct PictureManagerReorderEntry_s {
    uint64_t                                picture_number;
    EbObjectWrapper_t                    *parentPcsWrapperPtr;
} PictureManagerReorderEntry_t;

extern EbErrorType PictureManagerReorderEntryCtor(
    PictureManagerReorderEntry_t       **entryDblPtr,
    uint32_t                                picture_number);

#endif //EbPictureManagerReorderQueue_h
