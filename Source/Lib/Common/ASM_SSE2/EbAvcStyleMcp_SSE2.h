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

#ifndef EBAVCSTYLEMCP_SSE2_H
#define EBAVCSTYLEMCP_SSE2_H

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

void avc_style_copy_sse2(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride,
                         uint32_t pu_width, uint32_t pu_height, EbByte temp_buf,
                         uint32_t frac_pos);

#ifdef __cplusplus
}
#endif
#endif //EBAVCSTYLEMCP_H
