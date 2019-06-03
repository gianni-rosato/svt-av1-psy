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

    uint64_t compute8x8_satd_sse4(
        int16_t *diff);       // input parameter, diff samples Ptr

    uint64_t compute8x8_satd_u8_sse4(
        uint8_t  *src,       // input parameter, diff samples Ptr
        uint64_t *dc_value,
        uint32_t  src_stride);

#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_SSE4_1_h
