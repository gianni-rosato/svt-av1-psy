/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbComputeMean_SS2_h
#define EbComputeMean_SS2_h
#ifdef __cplusplus
extern "C" {
#endif

#include "EbDefinitions.h"


uint64_t compute_subd_mean_of_squared_values8x8_sse2_intrin(
    uint8_t* input_samples, // input parameter, input samples Ptr
    uint16_t input_stride);

#ifdef __cplusplus
}
#endif

#endif
