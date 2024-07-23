/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved.
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>

#include "aom_dsp_rtcd.h"
#include "coding_unit.h"
#include "definitions.h"
#include "filter.h"

DECLARE_ALIGNED(256, static const InterpKernel, av1_bilinear_filters[SUBPEL_SHIFTS]) = {{0, 0, 0, 128, 0, 0, 0, 0},
                                                                                        {0, 0, 0, 120, 8, 0, 0, 0},
                                                                                        {0, 0, 0, 112, 16, 0, 0, 0},
                                                                                        {0, 0, 0, 104, 24, 0, 0, 0},
                                                                                        {0, 0, 0, 96, 32, 0, 0, 0},
                                                                                        {0, 0, 0, 88, 40, 0, 0, 0},
                                                                                        {0, 0, 0, 80, 48, 0, 0, 0},
                                                                                        {0, 0, 0, 72, 56, 0, 0, 0},
                                                                                        {0, 0, 0, 64, 64, 0, 0, 0},
                                                                                        {0, 0, 0, 56, 72, 0, 0, 0},
                                                                                        {0, 0, 0, 48, 80, 0, 0, 0},
                                                                                        {0, 0, 0, 40, 88, 0, 0, 0},
                                                                                        {0, 0, 0, 32, 96, 0, 0, 0},
                                                                                        {0, 0, 0, 24, 104, 0, 0, 0},
                                                                                        {0, 0, 0, 16, 112, 0, 0, 0},
                                                                                        {0, 0, 0, 8, 120, 0, 0, 0}};

DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_4[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 0, -4, 126, 8, -2, 0, 0},
    {0, 0, -8, 122, 18, -4, 0, 0},
    {0, 0, -10, 116, 28, -6, 0, 0},
    {0, 0, -12, 110, 38, -8, 0, 0},
    {0, 0, -12, 102, 48, -10, 0, 0},
    {0, 0, -14, 94, 58, -10, 0, 0},
    {0, 0, -12, 84, 66, -10, 0, 0},
    {0, 0, -12, 76, 76, -12, 0, 0},
    {0, 0, -10, 66, 84, -12, 0, 0},
    {0, 0, -10, 58, 94, -14, 0, 0},
    {0, 0, -10, 48, 102, -12, 0, 0},
    {0, 0, -8, 38, 110, -12, 0, 0},
    {0, 0, -6, 28, 116, -10, 0, 0},
    {0, 0, -4, 18, 122, -8, 0, 0},
    {0, 0, -2, 8, 126, -4, 0, 0}};
DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_4smooth[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 0, 30, 62, 34, 2, 0, 0},
    {0, 0, 26, 62, 36, 4, 0, 0},
    {0, 0, 22, 62, 40, 4, 0, 0},
    {0, 0, 20, 60, 42, 6, 0, 0},
    {0, 0, 18, 58, 44, 8, 0, 0},
    {0, 0, 16, 56, 46, 10, 0, 0},
    {0, 0, 14, 54, 48, 12, 0, 0},
    {0, 0, 12, 52, 52, 12, 0, 0},
    {0, 0, 12, 48, 54, 14, 0, 0},
    {0, 0, 10, 46, 56, 16, 0, 0},
    {0, 0, 8, 44, 58, 18, 0, 0},
    {0, 0, 6, 42, 60, 20, 0, 0},
    {0, 0, 4, 40, 62, 22, 0, 0},
    {0, 0, 4, 36, 62, 26, 0, 0},
    {0, 0, 2, 34, 62, 30, 0, 0}};
DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_8[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 2, -6, 126, 8, -2, 0, 0},
    {0, 2, -10, 122, 18, -4, 0, 0},
    {0, 2, -12, 116, 28, -8, 2, 0},
    {0, 2, -14, 110, 38, -10, 2, 0},
    {0, 2, -14, 102, 48, -12, 2, 0},
    {0, 2, -16, 94, 58, -12, 2, 0},
    {0, 2, -14, 84, 66, -12, 2, 0},
    {0, 2, -14, 76, 76, -14, 2, 0},
    {0, 2, -12, 66, 84, -14, 2, 0},
    {0, 2, -12, 58, 94, -16, 2, 0},
    {0, 2, -12, 48, 102, -14, 2, 0},
    {0, 2, -10, 38, 110, -14, 2, 0},
    {0, 2, -8, 28, 116, -12, 2, 0},
    {0, 0, -4, 18, 122, -10, 2, 0},
    {0, 0, -2, 8, 126, -6, 2, 0}};

DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_8sharp[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {-2, 2, -6, 126, 8, -2, 2, 0},
    {-2, 6, -12, 124, 16, -6, 4, -2},
    {-2, 8, -18, 120, 26, -10, 6, -2},
    {-4, 10, -22, 116, 38, -14, 6, -2},
    {-4, 10, -22, 108, 48, -18, 8, -2},
    {-4, 10, -24, 100, 60, -20, 8, -2},
    {-4, 10, -24, 90, 70, -22, 10, -2},
    {-4, 12, -24, 80, 80, -24, 12, -4},
    {-2, 10, -22, 70, 90, -24, 10, -4},
    {-2, 8, -20, 60, 100, -24, 10, -4},
    {-2, 8, -18, 48, 108, -22, 10, -4},
    {-2, 6, -14, 38, 116, -22, 10, -4},
    {-2, 6, -10, 26, 120, -18, 8, -2},
    {-2, 4, -6, 16, 124, -12, 6, -2},
    {0, 2, -2, 8, 126, -6, 2, -2}};

DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_8smooth[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 2, 28, 62, 34, 2, 0, 0},
    {0, 0, 26, 62, 36, 4, 0, 0},
    {0, 0, 22, 62, 40, 4, 0, 0},
    {0, 0, 20, 60, 42, 6, 0, 0},
    {0, 0, 18, 58, 44, 8, 0, 0},
    {0, 0, 16, 56, 46, 10, 0, 0},
    {0, -2, 16, 54, 48, 12, 0, 0},
    {0, -2, 14, 52, 52, 14, -2, 0},
    {0, 0, 12, 48, 54, 16, -2, 0},
    {0, 0, 10, 46, 56, 16, 0, 0},
    {0, 0, 8, 44, 58, 18, 0, 0},
    {0, 0, 6, 42, 60, 20, 0, 0},
    {0, 0, 4, 40, 62, 22, 0, 0},
    {0, 0, 4, 36, 62, 26, 0, 0},
    {0, 0, 2, 34, 62, 28, 2, 0}};
// For w<=4, MULTITAP_SHARP is the same as EIGHTTAP_REGULAR
static const InterpFilterParams av1_interp_4tap[SWITCHABLE_FILTERS + 1] = {
    {(const int16_t *)av1_sub_pel_filters_4, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_REGULAR},
    {(const int16_t *)av1_sub_pel_filters_4smooth, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_SMOOTH},
    {(const int16_t *)av1_sub_pel_filters_4, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_REGULAR},
    {(const int16_t *)av1_bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS, BILINEAR},
};
static const InterpFilterParams av1_interp_filter_params_list[SWITCHABLE_FILTERS + 1] = {
    {(const int16_t *)av1_sub_pel_filters_8, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_REGULAR},
    {(const int16_t *)av1_sub_pel_filters_8smooth, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_SMOOTH},
    {(const int16_t *)av1_sub_pel_filters_8sharp, SUBPEL_TAPS, SUBPEL_SHIFTS, MULTITAP_SHARP},
    {(const int16_t *)av1_bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS, BILINEAR}};
static INLINE const InterpFilterParams *get_4tap_interp_filter_params(const InterpFilter interp_filter) {
    return &av1_interp_4tap[interp_filter];
}
static INLINE const InterpFilterParams *av1_get_filter(int subpel_search) {
    assert(subpel_search >= USE_2_TAPS);

    switch (subpel_search) {
    case USE_2_TAPS: return get_4tap_interp_filter_params(BILINEAR);
    case USE_4_TAPS: return get_4tap_interp_filter_params(EIGHTTAP_REGULAR);
    case USE_8_TAPS: return &av1_interp_filter_params_list[EIGHTTAP_REGULAR];
    default: assert(0); return NULL;
    }
}

// Get pred block from up-sampled reference.
void svt_aom_upsampled_pred_neon(MacroBlockD                  *xd,
                                 const struct AV1Common *const cm, //const AV1_COMMON *const cm,
                                 int mi_row, int mi_col, const MV *const mv, uint8_t *comp_pred, int width, int height,
                                 int subpel_x_q3, int subpel_y_q3, const uint8_t *ref, int ref_stride,
                                 int subpel_search) {
    (void)xd;
    (void)cm;
    (void)mi_row;
    (void)mi_col;
    (void)mv;
    const InterpFilterParams *filter = av1_get_filter(subpel_search);
    assert(filter != NULL);
    if (!subpel_x_q3 && !subpel_y_q3) {
        for (int i = 0; i < height; i++) {
            svt_memcpy(comp_pred, ref, width * sizeof(*comp_pred));
            comp_pred += width;
            ref += ref_stride;
        }
    } else if (!subpel_y_q3) {
        const int16_t *const kernel = av1_get_interp_filter_subpel_kernel(*filter, subpel_x_q3 << 1);
        svt_aom_convolve8_horiz_neon(ref, ref_stride, comp_pred, width, kernel, 16, NULL, -1, width, height);
    } else if (!subpel_x_q3) {
        const int16_t *const kernel = av1_get_interp_filter_subpel_kernel(*filter, subpel_y_q3 << 1);
        svt_aom_convolve8_vert_neon(ref, ref_stride, comp_pred, width, NULL, -1, kernel, 16, width, height);
    } else {
        DECLARE_ALIGNED(16, uint8_t, temp[((MAX_SB_SIZE * 2 + 16) + 16) * MAX_SB_SIZE]);
        const int16_t *const kernel_x            = av1_get_interp_filter_subpel_kernel(*filter, subpel_x_q3 << 1);
        const int16_t *const kernel_y            = av1_get_interp_filter_subpel_kernel(*filter, subpel_y_q3 << 1);
        const int            intermediate_height = (((height - 1) * 8 + subpel_y_q3) >> 3) + filter->taps;
        assert(intermediate_height <= (MAX_SB_SIZE * 2 + 16) + 16);
        svt_aom_convolve8_horiz_neon(ref - ref_stride * ((filter->taps >> 1) - 1),
                                     ref_stride,
                                     temp,
                                     MAX_SB_SIZE,
                                     kernel_x,
                                     16,
                                     NULL,
                                     -1,
                                     width,
                                     intermediate_height);
        svt_aom_convolve8_vert_neon(temp + MAX_SB_SIZE * ((filter->taps >> 1) - 1),
                                    MAX_SB_SIZE,
                                    comp_pred,
                                    width,
                                    NULL,
                                    -1,
                                    kernel_y,
                                    16,
                                    width,
                                    height);
    }
}
