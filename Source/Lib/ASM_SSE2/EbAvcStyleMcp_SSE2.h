/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EBAVCSTYLEMCP_SSE2_H
#define EBAVCSTYLEMCP_SSE2_H

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    void AvcStyleCopy_SSE2(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);

#ifdef __cplusplus
}
#endif
#endif //EBAVCSTYLEMCP_H