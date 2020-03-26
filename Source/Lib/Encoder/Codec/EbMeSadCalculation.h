/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMeSadCalculation_h
#define EbMeSadCalculation_h


#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

void initialize_buffer_32bits_c(uint32_t* pointer, uint32_t count128, uint32_t count32,
                                uint32_t value);
#ifdef __cplusplus
}
#endif
#endif // EbMeSadCalculation_h
