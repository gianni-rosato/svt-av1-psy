/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include "aom_dsp_rtcd.h"
#include "synonyms.h"

////////////////////////////////////////////////////////////////////////////////
// 8 bit
////////////////////////////////////////////////////////////////////////////////
// 2 tap bilinear filters
#define BIL_SUBPEL_BITS 3
#define BIL_SUBPEL_SHIFTS (1 << BIL_SUBPEL_BITS)
static const uint8_t bilinear_filters_2t[BIL_SUBPEL_SHIFTS][2] = {
    {128, 0},
    {112, 16},
    {96, 32},
    {80, 48},
    {64, 64},
    {48, 80},
    {32, 96},
    {16, 112},
};

void svt_aom_var_filter_block2d_bil_first_pass_ssse3(const uint8_t *a, uint16_t *b,
                                                 unsigned int src_pixels_per_line,
                                                 unsigned int pixel_step,
                                                 unsigned int output_height,
                                                 unsigned int output_width, const uint8_t *filter);

void svt_aom_var_filter_block2d_bil_second_pass_ssse3(const uint16_t *a, uint8_t *b,
                                                  unsigned int src_pixels_per_line,
                                                  unsigned int pixel_step,
                                                  unsigned int output_height,
                                                  unsigned int output_width, const uint8_t *filter);

#define OBMC_SUBPIX_VAR(W, H)                                                         \
    uint32_t svt_aom_obmc_sub_pixel_variance##W##x##H##_sse4_1(const uint8_t *pre,    \
                                                           int            pre_stride, \
                                                           int            xoffset,    \
                                                           int            yoffset,    \
                                                           const int32_t *wsrc,       \
                                                           const int32_t *mask,       \
                                                           unsigned int * sse) {      \
        uint16_t fdata3[(H + 1) * W];                                                 \
        uint8_t  temp2[H * W];                                                        \
                                                                                      \
        svt_aom_var_filter_block2d_bil_first_pass_ssse3(                              \
            pre, fdata3, pre_stride, 1, H + 1, W, bilinear_filters_2t[xoffset]);      \
        svt_aom_var_filter_block2d_bil_second_pass_ssse3(                             \
            fdata3, temp2, W, W, H, W, bilinear_filters_2t[yoffset]);                 \
                                                                                      \
        return svt_aom_obmc_variance##W##x##H##_avx2(temp2, W, wsrc, mask, sse);      \
    }

OBMC_SUBPIX_VAR(128, 128)
OBMC_SUBPIX_VAR(128, 64)
OBMC_SUBPIX_VAR(64, 128)
OBMC_SUBPIX_VAR(64, 64)
OBMC_SUBPIX_VAR(64, 32)
OBMC_SUBPIX_VAR(32, 64)
OBMC_SUBPIX_VAR(32, 32)
OBMC_SUBPIX_VAR(32, 16)
OBMC_SUBPIX_VAR(16, 32)
OBMC_SUBPIX_VAR(16, 16)
OBMC_SUBPIX_VAR(16, 8)
OBMC_SUBPIX_VAR(8, 16)
OBMC_SUBPIX_VAR(8, 8)
OBMC_SUBPIX_VAR(8, 4)
OBMC_SUBPIX_VAR(4, 8)
OBMC_SUBPIX_VAR(4, 4)
OBMC_SUBPIX_VAR(4, 16)
OBMC_SUBPIX_VAR(16, 4)
OBMC_SUBPIX_VAR(8, 32)
OBMC_SUBPIX_VAR(32, 8)
OBMC_SUBPIX_VAR(16, 64)
OBMC_SUBPIX_VAR(64, 16)
