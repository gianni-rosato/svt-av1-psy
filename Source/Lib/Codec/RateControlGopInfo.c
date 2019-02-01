/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "RateControlModel.h"
#include "RateControlGopInfo.h"

RateControlGopInfo_t *get_gop_infos(RateControlGopInfo_t *gop_info,
                                    uint64_t position)
{
    RateControlGopInfo_t    *current;

    while (position >= 0) {
        current = &gop_info[position];

        if (current->exists) {
            return current;
        }

        if (position == 0) {
            return EB_NULL;
        }
        position--;
    }

    return EB_NULL;
}
