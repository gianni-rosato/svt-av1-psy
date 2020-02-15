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

#ifndef EbDecCdef_h
#define EbDecCdef_h

#include "EbDecHandle.h"
#include "EbDecUtils.h"
#include "EbCdef.h"

#ifdef __cplusplus
extern "C" {
#endif

void svt_cdef_frame(EbDecHandle *dec_handle, int enable_flag);

void svt_cdef_sb_row_mt(EbDecHandle *dec_handle, int32_t *mi_wide_l2, int32_t *mi_high_l2,
                        uint16_t **colbuf, int32_t sb_fbr, uint16_t *src,
                        int32_t *curr_recon_stride, uint8_t **curr_blk_recon_buf);

#ifdef __cplusplus
}
#endif
#endif // EbDecCdef_h_
