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
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbResize_h
#define EbResize_h

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#include "EbInterPrediction.h"
#include "EbSequenceControlSet.h"
#include "EbSuperRes.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint16_t encoding_width;
    uint16_t encoding_height;
    uint8_t  superres_denom;
} superres_params_type;

EbErrorType av1_resize_and_extend_frame(const EbPictureBufferDesc* src, EbPictureBufferDesc* dst,
                                        int bd, const int num_planes, const uint32_t ss_x,
                                        const uint32_t ss_y);

void init_resize_picture(SequenceControlSet* scs_ptr, PictureParentControlSet* pcs_ptr);

#define filteredinterp_filters1000 av1_resize_filter_normative

#ifdef __cplusplus
}
#endif
#endif // EbResize_h
