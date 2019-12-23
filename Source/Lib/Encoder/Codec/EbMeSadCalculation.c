/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#include "EbMeSadCalculation.h"

void initialize_buffer_32bits_c(uint32_t* pointer, uint32_t count128, uint32_t count32,
                                uint32_t value) {
    uint32_t counter;
    for (counter = 0; counter < count128 * 4 + count32; ++counter) { pointer[counter] = value; }
}
