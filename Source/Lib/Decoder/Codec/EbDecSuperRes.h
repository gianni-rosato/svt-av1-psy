/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef EbDecSuperRes_h
#define EbDecSuperRes_h

#include "EbAv1Structs.h"
#include "EbDecUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

void av1_upscale_normative_rows(FrameHeader *frm_hdr, const uint8_t *src,
    int src_stride, uint8_t *dst, int dst_stride, int rows, int sub_x, int bd);

void av1_upscale_normative_and_extend_frame(FrameHeader *frm_hdr,
    SeqHeader*seq_hdr, EbPictureBufferDesc *src, EbPictureBufferDesc *dst);

void av1_superres_upscale(FrameHeader *frm_hdr, SeqHeader*seq_hdr,
    EbPictureBufferDesc *recon_picture_src);

#ifdef __cplusplus
}
#endif
#endif // EbDecSuperRes_h
