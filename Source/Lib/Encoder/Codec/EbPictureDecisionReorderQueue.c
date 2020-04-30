/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbPictureDecisionReorderQueue.h"

EbErrorType picture_decision_reorder_entry_ctor(PictureDecisionReorderEntry *entry_ptr,
                                                uint32_t                     picture_number) {
    entry_ptr->picture_number         = picture_number;
    return EB_ErrorNone;
}

