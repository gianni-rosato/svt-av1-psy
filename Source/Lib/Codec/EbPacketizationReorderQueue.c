/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPacketizationReorderQueue.h"

EbErrorType packetization_reorder_entry_ctor(
    PacketizationReorderEntry_t   **entry_dbl_ptr,
    uint32_t                          picture_number)
{
    EB_MALLOC(PacketizationReorderEntry_t*, *entry_dbl_ptr, sizeof(PacketizationReorderEntry_t), EB_N_PTR);

    (*entry_dbl_ptr)->picture_number = picture_number;
    (*entry_dbl_ptr)->output_stream_wrapper_ptr = (EbObjectWrapper_t *)EB_NULL;
    (*entry_dbl_ptr)->outputStatisticsWrapperPtr = (EbObjectWrapper_t *)EB_NULL;
    (*entry_dbl_ptr)->outMetaData = (EbLinkedListNode*)EB_NULL;

    return EB_ErrorNone;
}

