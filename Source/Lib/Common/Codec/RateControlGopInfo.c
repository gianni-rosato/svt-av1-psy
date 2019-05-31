/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "RateControlModel.h"
#include "RateControlGopInfo.h"

EbRateControlGopInfo *get_gop_infos(EbRateControlGopInfo *gop_info,
                                    uint64_t position)
{
    EbRateControlGopInfo    *current;

    while (1) { // First frame is always guaranteed to exist
        current = &gop_info[position];

        if (current->exists)
            return current;
        if (position == 0)
            return EB_NULL;
        position--;
    }

    return EB_NULL;
}
