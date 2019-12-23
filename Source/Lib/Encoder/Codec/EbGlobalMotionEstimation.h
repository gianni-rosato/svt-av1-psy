/*
* Copyright(c) 2019 Intel Corporation
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

#ifndef EbGlobalMotionEstimation_h
#define EbGlobalMotionEstimation_h

#include "EbPictureBufferDesc.h"
#include "EbMotionEstimationContext.h"

void global_motion_estimation(PictureParentControlSet *pcs_ptr, MeContext *context_ptr,
                              EbPictureBufferDesc *input_picture_ptr);
void compute_global_motion(EbPictureBufferDesc *input_pic, EbPictureBufferDesc *ref_pic,
                           EbWarpedMotionParams *bestWarpedMotion, int allow_high_precision_mv);

#endif // EbGlobalMotionEstimation_h
