/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbComputeVariance_C_h
#define EbComputeVariance_C_h
#ifdef __cplusplus
extern "C" {
#endif

#include "EbDefinitions.h"

uint32_t variance_highbd_c(const uint16_t *a, int a_stride, const uint16_t *b, int b_stride, int w,
                           int h, uint32_t *sse);

#ifdef __cplusplus
}
#endif
#endif // EbComputeVariance_C_h
