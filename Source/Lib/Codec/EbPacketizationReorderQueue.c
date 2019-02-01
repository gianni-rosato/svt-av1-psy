/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPacketizationReorderQueue.h"

EbErrorType PacketizationReorderEntryCtor(
    PacketizationReorderEntry_t   **entryDblPtr,
    uint32_t                          picture_number)
{
    EB_MALLOC(PacketizationReorderEntry_t*, *entryDblPtr, sizeof(PacketizationReorderEntry_t), EB_N_PTR);

    (*entryDblPtr)->picture_number = picture_number;
    (*entryDblPtr)->output_stream_wrapper_ptr = (EbObjectWrapper_t *)EB_NULL;
    (*entryDblPtr)->outputStatisticsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;
    (*entryDblPtr)->outMetaData = (EbLinkedListNode*)EB_NULL;

    return EB_ErrorNone;
}

