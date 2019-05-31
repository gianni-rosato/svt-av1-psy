/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"

static const uint32_t PartitionCountLut[SIZE_PART_MODE] = {
    1, //SIZE_2Nx2N
    2, //SIZE_2NxN
    2, //SIZE_Nx2N
    4, //SIZE_NxN
    2, //SIZE_2NxnU
    2, //SIZE_2NxnD
    2, //SIZE_nLx2N
    2  //SIZE_nRx2N
};

uint32_t PartitionCount(const EbPartMode partMode)
{
    return PartitionCountLut[partMode];
    //uint32_t count = 2;
    //count -= (partMode == SIZE_2Nx2N);
    //count += (partMode == SIZE_NxN);
    //return count;
}
