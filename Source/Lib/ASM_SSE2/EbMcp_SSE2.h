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
#define USE_PRE_COMPUTE             0
    extern EB_ALIGN(16) const int16_t IntraPredictionConst_SSE2[344];
    /**************************************************
    * Assembly Declarations
    **************************************************/
    extern void PictureCopyKernel_SSE2(EbByte src, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t areaWidth, uint32_t areaHeight);
    void PictureAverageKernel_SSE2(EbByte src0, uint32_t src0Stride, EbByte src1, uint32_t src1Stride, EbByte dst, uint32_t dst_stride, uint32_t areaWidth, uint32_t areaHeight);
    void PictureAverageKernel_SSE2_INTRIN(EbByte src0, uint32_t src0Stride, EbByte src1, uint32_t src1Stride, EbByte dst, uint32_t dst_stride, uint32_t areaWidth, uint32_t areaHeight);
    void PictureAverageKernel1Line_SSE2_INTRIN(
        EbByte                  src0,
        EbByte                  src1,
        EbByte                  dst,
        uint32_t                   areaWidth);


#ifdef __cplusplus
}
#endif
#endif //EBMCP_SSE2_H
