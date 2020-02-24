/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPacketizationReorderQueue.h"

static void packetization_reorder_entry_dctor(EbPtr p) {
    PacketizationReorderEntry* obj = (PacketizationReorderEntry*)p;
    EB_DELETE(obj->bitstream_ptr);
}

EbErrorType packetization_reorder_entry_ctor(PacketizationReorderEntry *entry_ptr,
                                             uint32_t                   picture_number) {
    entry_ptr->dctor = packetization_reorder_entry_dctor;
    entry_ptr->picture_number = picture_number;
    //16 should enough for show existing frame
    EB_NEW(entry_ptr->bitstream_ptr, bitstream_ctor, 16);
    return EB_ErrorNone;
}
