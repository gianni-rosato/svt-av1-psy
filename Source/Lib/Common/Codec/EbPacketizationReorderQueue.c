// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include "EbPacketizationReorderQueue.h"

EbErrorType packetization_reorder_entry_ctor(
    PacketizationReorderEntry   *entry_dbl_ptr,
    uint32_t                          picture_number)
{
    entry_dbl_ptr->picture_number = picture_number;

    return EB_ErrorNone;
}
// clang-format on
