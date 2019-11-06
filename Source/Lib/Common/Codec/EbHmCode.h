/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbHmCode_h
#define EbHmCode_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif
    uint64_t compute4x4_satd(
        int16_t *diff);

    uint64_t compute4x4_satd_u8(
        uint8_t  *src,
        uint64_t *dc_value,
        uint32_t  src_stride);

    uint64_t compute8x8_satd_u8_c(
        uint8_t* src,
        uint64_t* dcValue,
        uint32_t  srcStride);
#ifdef __cplusplus
}
#endif
#endif // EbHmCode_h
