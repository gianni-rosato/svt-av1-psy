/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbGlobalMotionEstimation_h
#define EbGlobalMotionEstimation_h

#include "EbPictureBufferDesc.h"
#include "EbMotionEstimationContext.h"

void svt_aom_global_motion_estimation(PictureParentControlSet *pcs, EbPictureBufferDesc *input_pic);

void compute_global_motion(PictureParentControlSet *pcs, int *frm_corners, int num_frm_corners,
                           EbPictureBufferDesc *det_input_pic, //src frame for detection
                           EbPictureBufferDesc *det_ref_pic, //ref frame for detection
                           EbPictureBufferDesc *input_pic, //src frame for refinement
                           EbPictureBufferDesc *ref_pic, //ref frame for refinement
                           uint8_t              sf, //downsacle factor between det and refinement
                           uint8_t chess_refn, EbWarpedMotionParams *best_wm, int allow_high_precision_mv);

void                    svt_aom_upscale_wm_params(EbWarpedMotionParams *wm_params, uint8_t scale_factor);
extern MvReferenceFrame svt_get_ref_frame_type(uint8_t list, uint8_t ref_idx);
#endif // EbGlobalMotionEstimation_h
