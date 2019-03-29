/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_SSE4_1_h
#define EbPictureOperators_SSE4_1_h

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    uint64_t Compute8x8Satd_SSE4(
        int16_t *diff);       // input parameter, diff samples Ptr

    uint64_t Compute8x8Satd_U8_SSE4(
        uint8_t  *src,       // input parameter, diff samples Ptr
        uint64_t *dc_value,
        uint32_t  src_stride);

    uint64_t SpatialFullDistortionKernel4x4_SSSE3_INTRIN(
        uint8_t   *input,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t SpatialFullDistortionKernel8x8_SSSE3_INTRIN(
        uint8_t   *input,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t SpatialFullDistortionKernel16MxN_SSSE3_INTRIN(
        uint8_t   *input,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);


#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_SSE4_1_h

