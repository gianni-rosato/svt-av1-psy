/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPictureManagerReorderQueue.h"

EbErrorType picture_manager_reorder_entry_ctor(PictureManagerReorderEntry *entry_ptr,
                                               uint32_t                    picture_number) {
    entry_ptr->picture_number         = picture_number;
    return EB_ErrorNone;
}
