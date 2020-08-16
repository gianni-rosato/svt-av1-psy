/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
