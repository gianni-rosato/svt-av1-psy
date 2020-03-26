/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBMCP_SSE2_H
#define EBMCP_SSE2_H

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif
#define USE_PRE_COMPUTE 0
extern EB_ALIGN(16) const int16_t intra_prediction_const_sse2[344];
/**************************************************
    * Assembly Declarations
    **************************************************/
extern void picture_copy_kernel_sse2(EbByte src, uint32_t src_stride, EbByte dst,
                                     uint32_t dst_stride, uint32_t area_width,
                                     uint32_t area_height);

void picture_average_kernel_sse2(EbByte src0, uint32_t src0_stride, EbByte src1,
                                 uint32_t src1_stride, EbByte dst, uint32_t dst_stride,
                                 uint32_t area_width, uint32_t area_height);

#ifdef __cplusplus
}
#endif
#endif //EBMCP_SSE2_H
