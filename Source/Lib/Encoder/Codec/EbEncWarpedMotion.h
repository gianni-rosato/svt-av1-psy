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

#ifndef EbEncWarpedMotion_h
#define EbEncWarpedMotion_h

#include "EbWarpedMotion.h"

#ifdef __cplusplus
extern "C" {
#endif

static INLINE int error_measure(int err) { return error_measure_lut[255 + err]; }

// Returns the error between the result of applying motion 'wm' to the frame
// described by 'ref' and the frame described by 'dst'.
int64_t eb_av1_warp_error(EbWarpedMotionParams *wm, int use_hbd, int bd, const uint8_t *ref,
                          int width, int height, int stride, uint8_t *dst, int p_col, int p_row,
                          int p_width, int p_height, int p_stride, int subsampling_x,
                          int subsampling_y, int64_t best_error);

// Returns the error between the frame described by 'ref' and the frame
// described by 'dst'.
int64_t eb_av1_frame_error(int use_hbd, int bd, const uint8_t *ref, int stride, uint8_t *dst,
                           int p_width, int p_height, int p_stride);


#ifdef __cplusplus
}
#endif
#endif // EbEncWarpedMotion_h
