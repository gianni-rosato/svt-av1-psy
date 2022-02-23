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

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "EbResize.h"

#include "aom_dsp_rtcd.h"
#include "EbModeDecisionProcess.h"

EbPictureBufferDesc *get_ref_pic_buffer(PictureControlSet *pcs_ptr, uint8_t is_highbd,
                                        uint8_t list_idx, uint8_t ref_idx);

#define DIVIDE_AND_ROUND(x, y) (((x) + ((y) >> 1)) / (y))

// Filters for factor of 2 downsampling.
static const int16_t av1_down2_symeven_half_filter[] = {56, 12, -3, -1};
static const int16_t av1_down2_symodd_half_filter[]  = {64, 35, 0, -3};

// Filters for interpolation (0.5-band) - note this also filters integer pels.
static const InterpKernel filteredinterp_filters500[(1 << RS_SUBPEL_BITS)] = {
    {-3, 0, 35, 64, 35, 0, -3, 0},    {-3, 0, 34, 64, 36, 0, -3, 0},
    {-3, -1, 34, 64, 36, 1, -3, 0},   {-3, -1, 33, 64, 37, 1, -3, 0},
    {-3, -1, 32, 64, 38, 1, -3, 0},   {-3, -1, 31, 64, 39, 1, -3, 0},
    {-3, -1, 31, 63, 39, 2, -3, 0},   {-2, -2, 30, 63, 40, 2, -3, 0},
    {-2, -2, 29, 63, 41, 2, -3, 0},   {-2, -2, 29, 63, 41, 3, -4, 0},
    {-2, -2, 28, 63, 42, 3, -4, 0},   {-2, -2, 27, 63, 43, 3, -4, 0},
    {-2, -3, 27, 63, 43, 4, -4, 0},   {-2, -3, 26, 62, 44, 5, -4, 0},
    {-2, -3, 25, 62, 45, 5, -4, 0},   {-2, -3, 25, 62, 45, 5, -4, 0},
    {-2, -3, 24, 62, 46, 5, -4, 0},   {-2, -3, 23, 61, 47, 6, -4, 0},
    {-2, -3, 23, 61, 47, 6, -4, 0},   {-2, -3, 22, 61, 48, 7, -4, -1},
    {-2, -3, 21, 60, 49, 7, -4, 0},   {-1, -4, 20, 60, 49, 8, -4, 0},
    {-1, -4, 20, 60, 50, 8, -4, -1},  {-1, -4, 19, 59, 51, 9, -4, -1},
    {-1, -4, 19, 59, 51, 9, -4, -1},  {-1, -4, 18, 58, 52, 10, -4, -1},
    {-1, -4, 17, 58, 52, 11, -4, -1}, {-1, -4, 16, 58, 53, 11, -4, -1},
    {-1, -4, 16, 57, 53, 12, -4, -1}, {-1, -4, 15, 57, 54, 12, -4, -1},
    {-1, -4, 15, 56, 54, 13, -4, -1}, {-1, -4, 14, 56, 55, 13, -4, -1},
    {-1, -4, 14, 55, 55, 14, -4, -1}, {-1, -4, 13, 55, 56, 14, -4, -1},
    {-1, -4, 13, 54, 56, 15, -4, -1}, {-1, -4, 12, 54, 57, 15, -4, -1},
    {-1, -4, 12, 53, 57, 16, -4, -1}, {-1, -4, 11, 53, 58, 16, -4, -1},
    {-1, -4, 11, 52, 58, 17, -4, -1}, {-1, -4, 10, 52, 58, 18, -4, -1},
    {-1, -4, 9, 51, 59, 19, -4, -1},  {-1, -4, 9, 51, 59, 19, -4, -1},
    {-1, -4, 8, 50, 60, 20, -4, -1},  {0, -4, 8, 49, 60, 20, -4, -1},
    {0, -4, 7, 49, 60, 21, -3, -2},   {-1, -4, 7, 48, 61, 22, -3, -2},
    {0, -4, 6, 47, 61, 23, -3, -2},   {0, -4, 6, 47, 61, 23, -3, -2},
    {0, -4, 5, 46, 62, 24, -3, -2},   {0, -4, 5, 45, 62, 25, -3, -2},
    {0, -4, 5, 45, 62, 25, -3, -2},   {0, -4, 5, 44, 62, 26, -3, -2},
    {0, -4, 4, 43, 63, 27, -3, -2},   {0, -4, 3, 43, 63, 27, -2, -2},
    {0, -4, 3, 42, 63, 28, -2, -2},   {0, -4, 3, 41, 63, 29, -2, -2},
    {0, -3, 2, 41, 63, 29, -2, -2},   {0, -3, 2, 40, 63, 30, -2, -2},
    {0, -3, 2, 39, 63, 31, -1, -3},   {0, -3, 1, 39, 64, 31, -1, -3},
    {0, -3, 1, 38, 64, 32, -1, -3},   {0, -3, 1, 37, 64, 33, -1, -3},
    {0, -3, 1, 36, 64, 34, -1, -3},   {0, -3, 0, 36, 64, 34, 0, -3},
};

// Filters for interpolation (0.625-band) - note this also filters integer pels.
static const InterpKernel filteredinterp_filters625[(1 << RS_SUBPEL_BITS)] = {
    {-1, -8, 33, 80, 33, -8, -1, 0}, {-1, -8, 31, 80, 34, -8, -1, 1},
    {-1, -8, 30, 80, 35, -8, -1, 1}, {-1, -8, 29, 80, 36, -7, -2, 1},
    {-1, -8, 28, 80, 37, -7, -2, 1}, {-1, -8, 27, 80, 38, -7, -2, 1},
    {0, -8, 26, 79, 39, -7, -2, 1},  {0, -8, 25, 79, 40, -7, -2, 1},
    {0, -8, 24, 79, 41, -7, -2, 1},  {0, -8, 23, 78, 42, -6, -2, 1},
    {0, -8, 22, 78, 43, -6, -2, 1},  {0, -8, 21, 78, 44, -6, -2, 1},
    {0, -8, 20, 78, 45, -5, -3, 1},  {0, -8, 19, 77, 47, -5, -3, 1},
    {0, -8, 18, 77, 48, -5, -3, 1},  {0, -8, 17, 77, 49, -5, -3, 1},
    {0, -8, 16, 76, 50, -4, -3, 1},  {0, -8, 15, 76, 51, -4, -3, 1},
    {0, -8, 15, 75, 52, -3, -4, 1},  {0, -7, 14, 74, 53, -3, -4, 1},
    {0, -7, 13, 74, 54, -3, -4, 1},  {0, -7, 12, 73, 55, -2, -4, 1},
    {0, -7, 11, 73, 56, -2, -4, 1},  {0, -7, 10, 72, 57, -1, -4, 1},
    {1, -7, 10, 71, 58, -1, -5, 1},  {0, -7, 9, 71, 59, 0, -5, 1},
    {1, -7, 8, 70, 60, 0, -5, 1},    {1, -7, 7, 69, 61, 1, -5, 1},
    {1, -6, 6, 68, 62, 1, -5, 1},    {0, -6, 6, 68, 62, 2, -5, 1},
    {1, -6, 5, 67, 63, 2, -5, 1},    {1, -6, 5, 66, 64, 3, -6, 1},
    {1, -6, 4, 65, 65, 4, -6, 1},    {1, -6, 3, 64, 66, 5, -6, 1},
    {1, -5, 2, 63, 67, 5, -6, 1},    {1, -5, 2, 62, 68, 6, -6, 0},
    {1, -5, 1, 62, 68, 6, -6, 1},    {1, -5, 1, 61, 69, 7, -7, 1},
    {1, -5, 0, 60, 70, 8, -7, 1},    {1, -5, 0, 59, 71, 9, -7, 0},
    {1, -5, -1, 58, 71, 10, -7, 1},  {1, -4, -1, 57, 72, 10, -7, 0},
    {1, -4, -2, 56, 73, 11, -7, 0},  {1, -4, -2, 55, 73, 12, -7, 0},
    {1, -4, -3, 54, 74, 13, -7, 0},  {1, -4, -3, 53, 74, 14, -7, 0},
    {1, -4, -3, 52, 75, 15, -8, 0},  {1, -3, -4, 51, 76, 15, -8, 0},
    {1, -3, -4, 50, 76, 16, -8, 0},  {1, -3, -5, 49, 77, 17, -8, 0},
    {1, -3, -5, 48, 77, 18, -8, 0},  {1, -3, -5, 47, 77, 19, -8, 0},
    {1, -3, -5, 45, 78, 20, -8, 0},  {1, -2, -6, 44, 78, 21, -8, 0},
    {1, -2, -6, 43, 78, 22, -8, 0},  {1, -2, -6, 42, 78, 23, -8, 0},
    {1, -2, -7, 41, 79, 24, -8, 0},  {1, -2, -7, 40, 79, 25, -8, 0},
    {1, -2, -7, 39, 79, 26, -8, 0},  {1, -2, -7, 38, 80, 27, -8, -1},
    {1, -2, -7, 37, 80, 28, -8, -1}, {1, -2, -7, 36, 80, 29, -8, -1},
    {1, -1, -8, 35, 80, 30, -8, -1}, {1, -1, -8, 34, 80, 31, -8, -1},
};

// Filters for interpolation (0.75-band) - note this also filters integer pels.
static const InterpKernel filteredinterp_filters750[(1 << RS_SUBPEL_BITS)] = {
    {2, -11, 25, 96, 25, -11, 2, 0}, {2, -11, 24, 96, 26, -11, 2, 0},
    {2, -11, 22, 96, 28, -11, 2, 0}, {2, -10, 21, 96, 29, -12, 2, 0},
    {2, -10, 19, 96, 31, -12, 2, 0}, {2, -10, 18, 95, 32, -11, 2, 0},
    {2, -10, 17, 95, 34, -12, 2, 0}, {2, -9, 15, 95, 35, -12, 2, 0},
    {2, -9, 14, 94, 37, -12, 2, 0},  {2, -9, 13, 94, 38, -12, 2, 0},
    {2, -8, 12, 93, 40, -12, 1, 0},  {2, -8, 11, 93, 41, -12, 1, 0},
    {2, -8, 9, 92, 43, -12, 1, 1},   {2, -8, 8, 92, 44, -12, 1, 1},
    {2, -7, 7, 91, 46, -12, 1, 0},   {2, -7, 6, 90, 47, -12, 1, 1},
    {2, -7, 5, 90, 49, -12, 1, 0},   {2, -6, 4, 89, 50, -12, 1, 0},
    {2, -6, 3, 88, 52, -12, 0, 1},   {2, -6, 2, 87, 54, -12, 0, 1},
    {2, -5, 1, 86, 55, -12, 0, 1},   {2, -5, 0, 85, 57, -12, 0, 1},
    {2, -5, -1, 84, 58, -11, 0, 1},  {2, -5, -2, 83, 60, -11, 0, 1},
    {2, -4, -2, 82, 61, -11, -1, 1}, {1, -4, -3, 81, 63, -10, -1, 1},
    {2, -4, -4, 80, 64, -10, -1, 1}, {1, -4, -4, 79, 66, -10, -1, 1},
    {1, -3, -5, 77, 67, -9, -1, 1},  {1, -3, -6, 76, 69, -9, -1, 1},
    {1, -3, -6, 75, 70, -8, -2, 1},  {1, -2, -7, 74, 71, -8, -2, 1},
    {1, -2, -7, 72, 72, -7, -2, 1},  {1, -2, -8, 71, 74, -7, -2, 1},
    {1, -2, -8, 70, 75, -6, -3, 1},  {1, -1, -9, 69, 76, -6, -3, 1},
    {1, -1, -9, 67, 77, -5, -3, 1},  {1, -1, -10, 66, 79, -4, -4, 1},
    {1, -1, -10, 64, 80, -4, -4, 2}, {1, -1, -10, 63, 81, -3, -4, 1},
    {1, -1, -11, 61, 82, -2, -4, 2}, {1, 0, -11, 60, 83, -2, -5, 2},
    {1, 0, -11, 58, 84, -1, -5, 2},  {1, 0, -12, 57, 85, 0, -5, 2},
    {1, 0, -12, 55, 86, 1, -5, 2},   {1, 0, -12, 54, 87, 2, -6, 2},
    {1, 0, -12, 52, 88, 3, -6, 2},   {0, 1, -12, 50, 89, 4, -6, 2},
    {0, 1, -12, 49, 90, 5, -7, 2},   {1, 1, -12, 47, 90, 6, -7, 2},
    {0, 1, -12, 46, 91, 7, -7, 2},   {1, 1, -12, 44, 92, 8, -8, 2},
    {1, 1, -12, 43, 92, 9, -8, 2},   {0, 1, -12, 41, 93, 11, -8, 2},
    {0, 1, -12, 40, 93, 12, -8, 2},  {0, 2, -12, 38, 94, 13, -9, 2},
    {0, 2, -12, 37, 94, 14, -9, 2},  {0, 2, -12, 35, 95, 15, -9, 2},
    {0, 2, -12, 34, 95, 17, -10, 2}, {0, 2, -11, 32, 95, 18, -10, 2},
    {0, 2, -12, 31, 96, 19, -10, 2}, {0, 2, -12, 29, 96, 21, -10, 2},
    {0, 2, -11, 28, 96, 22, -11, 2}, {0, 2, -11, 26, 96, 24, -11, 2},
};

// Filters for interpolation (0.875-band) - note this also filters integer pels.
static const InterpKernel filteredinterp_filters875[(1 << RS_SUBPEL_BITS)] = {
    {3, -8, 13, 112, 13, -8, 3, 0},   {2, -7, 12, 112, 15, -8, 3, -1},
    {3, -7, 10, 112, 17, -9, 3, -1},  {2, -6, 8, 112, 19, -9, 3, -1},
    {2, -6, 7, 112, 21, -10, 3, -1},  {2, -5, 6, 111, 22, -10, 3, -1},
    {2, -5, 4, 111, 24, -10, 3, -1},  {2, -4, 3, 110, 26, -11, 3, -1},
    {2, -4, 1, 110, 28, -11, 3, -1},  {2, -4, 0, 109, 30, -12, 4, -1},
    {1, -3, -1, 108, 32, -12, 4, -1}, {1, -3, -2, 108, 34, -13, 4, -1},
    {1, -2, -4, 107, 36, -13, 4, -1}, {1, -2, -5, 106, 38, -13, 4, -1},
    {1, -1, -6, 105, 40, -14, 4, -1}, {1, -1, -7, 104, 42, -14, 4, -1},
    {1, -1, -7, 103, 44, -15, 4, -1}, {1, 0, -8, 101, 46, -15, 4, -1},
    {1, 0, -9, 100, 48, -15, 4, -1},  {1, 0, -10, 99, 50, -15, 4, -1},
    {1, 1, -11, 97, 53, -16, 4, -1},  {0, 1, -11, 96, 55, -16, 4, -1},
    {0, 1, -12, 95, 57, -16, 4, -1},  {0, 2, -13, 93, 59, -16, 4, -1},
    {0, 2, -13, 91, 61, -16, 4, -1},  {0, 2, -14, 90, 63, -16, 4, -1},
    {0, 2, -14, 88, 65, -16, 4, -1},  {0, 2, -15, 86, 67, -16, 4, 0},
    {0, 3, -15, 84, 69, -17, 4, 0},   {0, 3, -16, 83, 71, -17, 4, 0},
    {0, 3, -16, 81, 73, -16, 3, 0},   {0, 3, -16, 79, 75, -16, 3, 0},
    {0, 3, -16, 77, 77, -16, 3, 0},   {0, 3, -16, 75, 79, -16, 3, 0},
    {0, 3, -16, 73, 81, -16, 3, 0},   {0, 4, -17, 71, 83, -16, 3, 0},
    {0, 4, -17, 69, 84, -15, 3, 0},   {0, 4, -16, 67, 86, -15, 2, 0},
    {-1, 4, -16, 65, 88, -14, 2, 0},  {-1, 4, -16, 63, 90, -14, 2, 0},
    {-1, 4, -16, 61, 91, -13, 2, 0},  {-1, 4, -16, 59, 93, -13, 2, 0},
    {-1, 4, -16, 57, 95, -12, 1, 0},  {-1, 4, -16, 55, 96, -11, 1, 0},
    {-1, 4, -16, 53, 97, -11, 1, 1},  {-1, 4, -15, 50, 99, -10, 0, 1},
    {-1, 4, -15, 48, 100, -9, 0, 1},  {-1, 4, -15, 46, 101, -8, 0, 1},
    {-1, 4, -15, 44, 103, -7, -1, 1}, {-1, 4, -14, 42, 104, -7, -1, 1},
    {-1, 4, -14, 40, 105, -6, -1, 1}, {-1, 4, -13, 38, 106, -5, -2, 1},
    {-1, 4, -13, 36, 107, -4, -2, 1}, {-1, 4, -13, 34, 108, -2, -3, 1},
    {-1, 4, -12, 32, 108, -1, -3, 1}, {-1, 4, -12, 30, 109, 0, -4, 2},
    {-1, 3, -11, 28, 110, 1, -4, 2},  {-1, 3, -11, 26, 110, 3, -4, 2},
    {-1, 3, -10, 24, 111, 4, -5, 2},  {-1, 3, -10, 22, 111, 6, -5, 2},
    {-1, 3, -10, 21, 112, 7, -6, 2},  {-1, 3, -9, 19, 112, 8, -6, 2},
    {-1, 3, -9, 17, 112, 10, -7, 3},  {-1, 3, -8, 15, 112, 12, -7, 2},
};

void downsample_decimation_input_picture(PictureParentControlSet *pcs_ptr,
                                         EbPictureBufferDesc     *input_padded_picture_ptr,
                                         EbPictureBufferDesc     *quarter_decimated_picture_ptr,
                                         EbPictureBufferDesc     *sixteenth_decimated_picture_ptr);

void downsample_filtering_input_picture(PictureParentControlSet *pcs_ptr,
                                        EbPictureBufferDesc     *input_padded_picture_ptr,
                                        EbPictureBufferDesc     *quarter_picture_ptr,
                                        EbPictureBufferDesc     *sixteenth_picture_ptr);

void calculate_scaled_size_helper(uint16_t *dim, uint8_t denom);

void pad_and_decimate_filtered_pic(PictureParentControlSet *picture_control_set_ptr_central);

static int get_down2_length(int length, int steps) {
    for (int s = 0; s < steps; ++s) length = (length + 1) >> 1;
    return length;
}

static int get_down2_steps(int in_length, int out_length) {
    int steps = 0;
    int proj_in_length;
    while ((proj_in_length = get_down2_length(in_length, 1)) >= out_length) {
        ++steps;
        in_length = proj_in_length;
        if (in_length == 1) {
            // Special case: we break because any further calls to get_down2_length()
            // with be with length == 1, which return 1, resulting in an infinite
            // loop.
            break;
        }
    }
    return steps;
}

static void down2_symeven(const uint8_t *const input, int length, uint8_t *output) {
    // Actual filter len = 2 * filter_len_half.
    const int16_t *filter          = av1_down2_symeven_half_filter;
    const int      filter_len_half = sizeof(av1_down2_symeven_half_filter) / 2;
    int            i, j;
    uint8_t       *optr = output;
    int            l1   = filter_len_half;
    int            l2   = (length - filter_len_half);
    l1 += (l1 & 1);
    l2 += (l2 & 1);
    if (l1 > l2) {
        // Short input length.
        for (i = 0; i < length; i += 2) {
            int sum = (1 << (FILTER_BITS - 1));
            for (j = 0; j < filter_len_half; ++j) {
                sum += (input[AOMMAX(i - j, 0)] + input[AOMMIN(i + 1 + j, length - 1)]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel(sum);
        }
    } else {
        // Initial part.
        for (i = 0; i < l1; i += 2) {
            int sum = (1 << (FILTER_BITS - 1));
            for (j = 0; j < filter_len_half; ++j) {
                sum += (input[AOMMAX(i - j, 0)] + input[i + 1 + j]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel(sum);
        }
        // Middle part.
        for (; i < l2; i += 2) {
            int sum = (1 << (FILTER_BITS - 1));
            for (j = 0; j < filter_len_half; ++j) {
                sum += (input[i - j] + input[i + 1 + j]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel(sum);
        }
        // End part.
        for (; i < length; i += 2) {
            int sum = (1 << (FILTER_BITS - 1));
            for (j = 0; j < filter_len_half; ++j) {
                sum += (input[i - j] + input[AOMMIN(i + 1 + j, length - 1)]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel(sum);
        }
    }
}

static void down2_symodd(const uint8_t *const input, int length, uint8_t *output) {
    // Actual filter len = 2 * filter_len_half - 1.
    const int16_t *filter          = av1_down2_symodd_half_filter;
    const int      filter_len_half = sizeof(av1_down2_symodd_half_filter) / 2;
    int            i, j;
    uint8_t       *optr = output;
    int            l1   = filter_len_half - 1;
    int            l2   = (length - filter_len_half + 1);
    l1 += (l1 & 1);
    l2 += (l2 & 1);
    if (l1 > l2) {
        // Short input length.
        for (i = 0; i < length; i += 2) {
            int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
            for (j = 1; j < filter_len_half; ++j) {
                sum += (input[(i - j < 0 ? 0 : i - j)] +
                        input[(i + j >= length ? length - 1 : i + j)]) *
                    filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel(sum);
        }
    } else {
        // Initial part.
        for (i = 0; i < l1; i += 2) {
            int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
            for (j = 1; j < filter_len_half; ++j) {
                sum += (input[(i - j < 0 ? 0 : i - j)] + input[i + j]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel(sum);
        }
        // Middle part.
        for (; i < l2; i += 2) {
            int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
            for (j = 1; j < filter_len_half; ++j) {
                sum += (input[i - j] + input[i + j]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel(sum);
        }
        // End part.
        for (; i < length; i += 2) {
            int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
            for (j = 1; j < filter_len_half; ++j) {
                sum += (input[i - j] + input[(i + j >= length ? length - 1 : i + j)]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel(sum);
        }
    }
}

static const InterpKernel *choose_interp_filter(int in_length, int out_length) {
    int out_length16 = out_length * 16;
    if (out_length16 >= in_length * 16)
        return filteredinterp_filters1000;
    else if (out_length16 >= in_length * 13)
        return filteredinterp_filters875;
    else if (out_length16 >= in_length * 11)
        return filteredinterp_filters750;
    else if (out_length16 >= in_length * 9)
        return filteredinterp_filters625;
    else
        return filteredinterp_filters500;
}

static void interpolate_core(const uint8_t *const input, int in_length, uint8_t *output,
                             int out_length, const int16_t *interp_filters, int interp_taps) {
    const int32_t delta = (((uint32_t)in_length << RS_SCALE_SUBPEL_BITS) + out_length / 2) /
        out_length;
    const int32_t offset = in_length > out_length
        ? (((int32_t)(in_length - out_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
            out_length
        : -(((int32_t)(out_length - in_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
            out_length;
    uint8_t      *optr   = output;
    int           x, x1, x2, sum, k, int_pel, sub_pel;
    int32_t       y;

    x = 0;
    y = offset + RS_SCALE_EXTRA_OFF;
    while ((y >> RS_SCALE_SUBPEL_BITS) < (interp_taps / 2 - 1)) {
        x++;
        y += delta;
    }
    x1 = x;
    x  = out_length - 1;
    y  = delta * x + offset + RS_SCALE_EXTRA_OFF;
    while ((y >> RS_SCALE_SUBPEL_BITS) + (int32_t)(interp_taps / 2) >= in_length) {
        x--;
        y -= delta;
    }
    x2 = x;
    if (x1 > x2) {
        for (x = 0, y = offset + RS_SCALE_EXTRA_OFF; x < out_length; ++x, y += delta) {
            int_pel               = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel               = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
            const int16_t *filter = &interp_filters[sub_pel * interp_taps];
            sum                   = 0;
            for (k = 0; k < interp_taps; ++k) {
                const int pk = int_pel - interp_taps / 2 + 1 + k;
                sum += filter[k] * input[AOMMAX(AOMMIN(pk, in_length - 1), 0)];
            }
            *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
        }
    } else {
        // Initial part.
        for (x = 0, y = offset + RS_SCALE_EXTRA_OFF; x < x1; ++x, y += delta) {
            int_pel               = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel               = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
            const int16_t *filter = &interp_filters[sub_pel * interp_taps];
            sum                   = 0;
            for (k = 0; k < interp_taps; ++k)
                sum += filter[k] * input[AOMMAX(int_pel - interp_taps / 2 + 1 + k, 0)];
            *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
        }
        // Middle part.
        for (; x <= x2; ++x, y += delta) {
            int_pel               = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel               = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
            const int16_t *filter = &interp_filters[sub_pel * interp_taps];
            sum                   = 0;
            for (k = 0; k < interp_taps; ++k)
                sum += filter[k] * input[int_pel - interp_taps / 2 + 1 + k];
            *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
        }
        // End part.
        for (; x < out_length; ++x, y += delta) {
            int_pel               = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel               = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
            const int16_t *filter = &interp_filters[sub_pel * interp_taps];
            sum                   = 0;
            for (k = 0; k < interp_taps; ++k)
                sum += filter[k] * input[AOMMIN(int_pel - interp_taps / 2 + 1 + k, in_length - 1)];
            *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
        }
    }
}

static void interpolate(const uint8_t *const input, int in_length, uint8_t *output,
                        int out_length) {
    const InterpKernel *interp_filters = choose_interp_filter(in_length, out_length);

    interpolate_core(input, in_length, output, out_length, &interp_filters[0][0], SUBPEL_TAPS);
}

static void resize_multistep(const uint8_t *const input, int length, uint8_t *output, int olength,
                             uint8_t *otmp) {
    if (length == olength) {
        svt_memcpy(output, input, sizeof(output[0]) * length);
        return;
    }
    const int steps = get_down2_steps(length, olength);

    if (steps > 0) {
        uint8_t *out            = NULL;
        int      filteredlength = length;

        assert(otmp != NULL);
        uint8_t *otmp2 = otmp + get_down2_length(length, 1);
        for (int s = 0; s < steps; ++s) {
            const int            proj_filteredlength = get_down2_length(filteredlength, 1);
            const uint8_t *const in                  = (s == 0 ? input : out);
            if (s == steps - 1 && proj_filteredlength == olength)
                out = output;
            else
                out = (s & 1) ? otmp2 : otmp;
            if (filteredlength & 1)
                down2_symodd(in, filteredlength, out);
            else
                down2_symeven(in, filteredlength, out);
            filteredlength = proj_filteredlength;
        }
        if (filteredlength != olength) {
            interpolate(out, filteredlength, output, olength);
        }
    } else {
        interpolate(input, length, output, olength);
    }
}

static void fill_arr_to_col(uint8_t *img, int stride, int len, uint8_t *arr) {
    int      i;
    uint8_t *iptr = img;
    uint8_t *aptr = arr;
    for (i = 0; i < len; ++i, iptr += stride) { *iptr = *aptr++; }
}

static void fill_col_to_arr(uint8_t *img, int stride, int len, uint8_t *arr) {
    int      i;
    uint8_t *iptr = img;
    uint8_t *aptr = arr;
    for (i = 0; i < len; ++i, iptr += stride) { *aptr++ = *iptr; }
}

static EbErrorType av1_resize_plane(const uint8_t *const input, int height, int width,
                                    int in_stride, uint8_t *output, int height2, int width2,
                                    int out_stride) {
    int      i;
    uint8_t *intbuf, *tmpbuf, *arrbuf, *arrbuf2;

    assert(width > 0);
    assert(height > 0);
    assert(width2 > 0);
    assert(height2 > 0);

    EB_MALLOC_ARRAY(intbuf, width2 * height);
    EB_MALLOC_ARRAY(tmpbuf, AOMMAX(width, height));
    EB_MALLOC_ARRAY(arrbuf, height);
    EB_MALLOC_ARRAY(arrbuf2, height2);
    if (intbuf == NULL || tmpbuf == NULL || arrbuf == NULL || arrbuf2 == NULL) {
        EB_FREE_ARRAY(intbuf);
        EB_FREE_ARRAY(tmpbuf);
        EB_FREE_ARRAY(arrbuf);
        EB_FREE_ARRAY(arrbuf2);
        return EB_ErrorInsufficientResources;
    }
    for (i = 0; i < height; ++i)
        resize_multistep(input + in_stride * i, width, intbuf + width2 * i, width2, tmpbuf);

    for (i = 0; i < width2; ++i) {
        fill_col_to_arr(intbuf + i, width2, height, arrbuf);
        resize_multistep(arrbuf, height, arrbuf2, height2, tmpbuf);
        fill_arr_to_col(output + i, out_stride, height2, arrbuf2);
    }

    EB_FREE_ARRAY(intbuf);
    EB_FREE_ARRAY(tmpbuf);
    EB_FREE_ARRAY(arrbuf);
    EB_FREE_ARRAY(arrbuf2);

    return EB_ErrorNone;
}

static EbErrorType av1_resize_plane_horizontal(const uint8_t *const input, int height, int width,
                                               int in_stride, uint8_t *output, int height2,
                                               int width2, int out_stride) {
    int      i;
    uint8_t *tmpbuf;

    assert(width > 0);
    assert(height > 0);
    assert(width2 > 0);
    assert(height2 == height);
    UNUSED(height2);

    EB_MALLOC_ARRAY(tmpbuf, AOMMAX(width, height));
    if (tmpbuf == NULL) {
        return EB_ErrorInsufficientResources;
    }
    for (i = 0; i < height; ++i)
        resize_multistep(input + in_stride * i, width, output + out_stride * i, width2, tmpbuf);

    EB_FREE_ARRAY(tmpbuf);

    return EB_ErrorNone;
}

static void highbd_interpolate_core(const uint16_t *const input, int in_length, uint16_t *output,
                                    int out_length, int bd, const int16_t *interp_filters,
                                    int interp_taps) {
    const int32_t delta = (((uint32_t)in_length << RS_SCALE_SUBPEL_BITS) + out_length / 2) /
        out_length;
    const int32_t offset = in_length > out_length
        ? (((int32_t)(in_length - out_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
            out_length
        : -(((int32_t)(out_length - in_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
            out_length;
    uint16_t     *optr   = output;
    int           x, x1, x2, sum, k, int_pel, sub_pel;
    int32_t       y;

    x = 0;
    y = offset + RS_SCALE_EXTRA_OFF;
    while ((y >> RS_SCALE_SUBPEL_BITS) < (interp_taps / 2 - 1)) {
        x++;
        y += delta;
    }
    x1 = x;
    x  = out_length - 1;
    y  = delta * x + offset + RS_SCALE_EXTRA_OFF;
    while ((y >> RS_SCALE_SUBPEL_BITS) + (int32_t)(interp_taps / 2) >= in_length) {
        x--;
        y -= delta;
    }
    x2 = x;
    if (x1 > x2) {
        for (x = 0, y = offset + RS_SCALE_EXTRA_OFF; x < out_length; ++x, y += delta) {
            int_pel               = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel               = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
            const int16_t *filter = &interp_filters[sub_pel * interp_taps];
            sum                   = 0;
            for (k = 0; k < interp_taps; ++k) {
                const int pk = int_pel - interp_taps / 2 + 1 + k;
                sum += filter[k] * input[AOMMAX(AOMMIN(pk, in_length - 1), 0)];
            }
            *optr++ = clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd);
        }
    } else {
        // Initial part.
        for (x = 0, y = offset + RS_SCALE_EXTRA_OFF; x < x1; ++x, y += delta) {
            int_pel               = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel               = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
            const int16_t *filter = &interp_filters[sub_pel * interp_taps];
            sum                   = 0;
            for (k = 0; k < interp_taps; ++k)
                sum += filter[k] * input[AOMMAX(int_pel - interp_taps / 2 + 1 + k, 0)];
            *optr++ = clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd);
        }
        // Middle part.
        for (; x <= x2; ++x, y += delta) {
            int_pel               = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel               = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
            const int16_t *filter = &interp_filters[sub_pel * interp_taps];
            sum                   = 0;
            for (k = 0; k < interp_taps; ++k)
                sum += filter[k] * input[int_pel - interp_taps / 2 + 1 + k];
            *optr++ = clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd);
        }
        // End part.
        for (; x < out_length; ++x, y += delta) {
            int_pel               = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel               = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
            const int16_t *filter = &interp_filters[sub_pel * interp_taps];
            sum                   = 0;
            for (k = 0; k < interp_taps; ++k)
                sum += filter[k] * input[AOMMIN(int_pel - interp_taps / 2 + 1 + k, in_length - 1)];
            *optr++ = clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd);
        }
    }
}

static void highbd_interpolate(const uint16_t *const input, int in_length, uint16_t *output,
                               int out_length, int bd) {
    const InterpKernel *interp_filters = choose_interp_filter(in_length, out_length);

    highbd_interpolate_core(
        input, in_length, output, out_length, bd, &interp_filters[0][0], SUBPEL_TAPS);
}

static void highbd_down2_symeven(const uint16_t *const input, int length, uint16_t *output,
                                 int bd) {
    // Actual filter len = 2 * filter_len_half.
    static const int16_t *filter          = av1_down2_symeven_half_filter;
    const int             filter_len_half = sizeof(av1_down2_symeven_half_filter) / 2;
    int                   i, j;
    uint16_t             *optr = output;
    int                   l1   = filter_len_half;
    int                   l2   = (length - filter_len_half);
    l1 += (l1 & 1);
    l2 += (l2 & 1);
    if (l1 > l2) {
        // Short input length.
        for (i = 0; i < length; i += 2) {
            int sum = (1 << (FILTER_BITS - 1));
            for (j = 0; j < filter_len_half; ++j) {
                sum += (input[AOMMAX(0, i - j)] + input[AOMMIN(i + 1 + j, length - 1)]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel_highbd(sum, bd);
        }
    } else {
        // Initial part.
        for (i = 0; i < l1; i += 2) {
            int sum = (1 << (FILTER_BITS - 1));
            for (j = 0; j < filter_len_half; ++j) {
                sum += (input[AOMMAX(0, i - j)] + input[i + 1 + j]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel_highbd(sum, bd);
        }
        // Middle part.
        for (; i < l2; i += 2) {
            int sum = (1 << (FILTER_BITS - 1));
            for (j = 0; j < filter_len_half; ++j) {
                sum += (input[i - j] + input[i + 1 + j]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel_highbd(sum, bd);
        }
        // End part.
        for (; i < length; i += 2) {
            int sum = (1 << (FILTER_BITS - 1));
            for (j = 0; j < filter_len_half; ++j) {
                sum += (input[i - j] + input[AOMMIN(i + 1 + j, length - 1)]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel_highbd(sum, bd);
        }
    }
}

static void highbd_down2_symodd(const uint16_t *const input, int length, uint16_t *output, int bd) {
    // Actual filter len = 2 * filter_len_half - 1.
    static const int16_t *filter          = av1_down2_symodd_half_filter;
    const int             filter_len_half = sizeof(av1_down2_symodd_half_filter) / 2;
    int                   i, j;
    uint16_t             *optr = output;
    int                   l1   = filter_len_half - 1;
    int                   l2   = (length - filter_len_half + 1);
    l1 += (l1 & 1);
    l2 += (l2 & 1);
    if (l1 > l2) {
        // Short input length.
        for (i = 0; i < length; i += 2) {
            int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
            for (j = 1; j < filter_len_half; ++j) {
                sum += (input[AOMMAX(i - j, 0)] + input[AOMMIN(i + j, length - 1)]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel_highbd(sum, bd);
        }
    } else {
        // Initial part.
        for (i = 0; i < l1; i += 2) {
            int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
            for (j = 1; j < filter_len_half; ++j) {
                sum += (input[AOMMAX(i - j, 0)] + input[i + j]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel_highbd(sum, bd);
        }
        // Middle part.
        for (; i < l2; i += 2) {
            int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
            for (j = 1; j < filter_len_half; ++j) {
                sum += (input[i - j] + input[i + j]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel_highbd(sum, bd);
        }
        // End part.
        for (; i < length; i += 2) {
            int sum = (1 << (FILTER_BITS - 1)) + input[i] * filter[0];
            for (j = 1; j < filter_len_half; ++j) {
                sum += (input[i - j] + input[AOMMIN(i + j, length - 1)]) * filter[j];
            }
            sum >>= FILTER_BITS;
            *optr++ = clip_pixel_highbd(sum, bd);
        }
    }
}

static void highbd_resize_multistep(const uint16_t *const input, int length, uint16_t *output,
                                    int olength, uint16_t *otmp, int bd) {
    if (length == olength) {
        svt_memcpy(output, input, sizeof(output[0]) * length);
        return;
    }
    const int steps = get_down2_steps(length, olength);

    if (steps > 0) {
        uint16_t *out            = NULL;
        int       filteredlength = length;

        assert(otmp != NULL);
        uint16_t *otmp2 = otmp + get_down2_length(length, 1);
        for (int s = 0; s < steps; ++s) {
            const int             proj_filteredlength = get_down2_length(filteredlength, 1);
            const uint16_t *const in                  = (s == 0 ? input : out);
            if (s == steps - 1 && proj_filteredlength == olength)
                out = output;
            else
                out = (s & 1) ? otmp2 : otmp;
            if (filteredlength & 1)
                highbd_down2_symodd(in, filteredlength, out, bd);
            else
                highbd_down2_symeven(in, filteredlength, out, bd);
            filteredlength = proj_filteredlength;
        }
        if (filteredlength != olength) {
            highbd_interpolate(out, filteredlength, output, olength, bd);
        }
    } else {
        highbd_interpolate(input, length, output, olength, bd);
    }
}

static void highbd_fill_col_to_arr(uint16_t *img, int stride, int len, uint16_t *arr) {
    int       i;
    uint16_t *iptr = img;
    uint16_t *aptr = arr;
    for (i = 0; i < len; ++i, iptr += stride) { *aptr++ = *iptr; }
}

static void highbd_fill_arr_to_col(uint16_t *img, int stride, int len, uint16_t *arr) {
    int       i;
    uint16_t *iptr = img;
    uint16_t *aptr = arr;
    for (i = 0; i < len; ++i, iptr += stride) { *iptr = *aptr++; }
}

static EbErrorType av1_highbd_resize_plane(const uint16_t *const input, int height, int width,
                                           int in_stride, uint16_t *output, int height2, int width2,
                                           int out_stride, int bd) {
    int       i;
    uint16_t *intbuf;
    uint16_t *tmpbuf;
    uint16_t *arrbuf;
    uint16_t *arrbuf2;

    EB_MALLOC_ARRAY(intbuf, sizeof(uint16_t) * width2 * height);
    EB_MALLOC_ARRAY(tmpbuf, sizeof(uint16_t) * AOMMAX(width, height));
    EB_MALLOC_ARRAY(arrbuf, sizeof(uint16_t) * height);
    EB_MALLOC_ARRAY(arrbuf2, sizeof(uint16_t) * height2);
    if (intbuf == NULL || tmpbuf == NULL || arrbuf == NULL || arrbuf2 == NULL) {
        EB_FREE(intbuf);
        EB_FREE(tmpbuf);
        EB_FREE(arrbuf);
        EB_FREE(arrbuf2);
        return EB_ErrorInsufficientResources;
    }
    for (i = 0; i < height; ++i) {
        highbd_resize_multistep(
            input + in_stride * i, width, intbuf + width2 * i, width2, tmpbuf, bd);
    }
    for (i = 0; i < width2; ++i) {
        highbd_fill_col_to_arr(intbuf + i, width2, height, arrbuf);
        highbd_resize_multistep(arrbuf, height, arrbuf2, height2, tmpbuf, bd);
        highbd_fill_arr_to_col(output + i, out_stride, height2, arrbuf2);
    }

    EB_FREE(intbuf);
    EB_FREE(tmpbuf);
    EB_FREE(arrbuf);
    EB_FREE(arrbuf2);

    return EB_ErrorNone;
}

static EbErrorType av1_highbd_resize_plane_horizontal(const uint16_t *const input, int height,
                                                      int width, int in_stride, uint16_t *output,
                                                      int height2, int width2, int out_stride,
                                                      int bd) {
    int       i;
    uint16_t *tmpbuf;

    UNUSED(height2);
    assert(height2 == height);

    EB_MALLOC_ARRAY(tmpbuf, sizeof(uint16_t) * AOMMAX(width, height));
    if (tmpbuf == NULL) {
        return EB_ErrorInsufficientResources;
    }
    for (i = 0; i < height; ++i) {
        highbd_resize_multistep(
            input + in_stride * i, width, output + out_stride * i, width2, tmpbuf, bd);
    }

    EB_FREE_ARRAY(tmpbuf);
    return EB_ErrorNone;
}

void pack_highbd_pic(const EbPictureBufferDesc *pic_ptr, uint16_t *buffer_16bit[3], uint32_t ss_x,
                     uint32_t ss_y, EbBool include_padding);

void unpack_highbd_pic(uint16_t *buffer_highbd[3], EbPictureBufferDesc *pic_ptr, uint32_t ss_x,
                       uint32_t ss_y, EbBool include_padding);
#if DEBUG_SCALING
void save_YUV_to_file(char *filename, EbByte buffer_y, EbByte buffer_u, EbByte buffer_v,
                      uint16_t width, uint16_t height, uint16_t stride_y, uint16_t stride_u,
                      uint16_t stride_v, uint16_t origin_y, uint16_t origin_x, uint32_t ss_x,
                      uint32_t ss_y);

void save_YUV_to_file_highbd(char *filename, uint16_t *buffer_y, uint16_t *buffer_u,
                             uint16_t *buffer_v, uint16_t width, uint16_t height, uint16_t stride_y,
                             uint16_t stride_u, uint16_t stride_v, uint16_t origin_y,
                             uint16_t origin_x, uint32_t ss_x, uint32_t ss_y);
#endif
typedef EbErrorType (*Av1HighbdResizePlane)(const uint16_t *const input, int height, int width,
                                            int in_stride, uint16_t *output, int height2,
                                            int width2, int out_stride, int bd);
typedef EbErrorType (*Av1ResizePlane)(const uint8_t *const input, int height, int width,
                                      int in_stride, uint8_t *output, int height2, int width2,
                                      int out_stride);

/*
 * Resize frame according to dst resolution.
 * Supports 8-bit / 10-bit and either packed or unpacked buffers
 */
static EbErrorType av1_resize_frame(const EbPictureBufferDesc *src, EbPictureBufferDesc *dst,
                                    int bd, const int num_planes, const uint32_t ss_x,
                                    const uint32_t ss_y, uint8_t is_packed,
                                    uint32_t buffer_enable_mask) {
    uint16_t *src_buffer_highbd[MAX_MB_PLANE];
    uint16_t *dst_buffer_highbd[MAX_MB_PLANE];

    if (bd > 8 && !is_packed) {
        EB_MALLOC_ARRAY(src_buffer_highbd[0], src->luma_size);
        EB_MALLOC_ARRAY(src_buffer_highbd[1], src->chroma_size);
        EB_MALLOC_ARRAY(src_buffer_highbd[2], src->chroma_size);
        EB_MALLOC_ARRAY(dst_buffer_highbd[0], dst->luma_size);
        EB_MALLOC_ARRAY(dst_buffer_highbd[1], dst->chroma_size);
        EB_MALLOC_ARRAY(dst_buffer_highbd[2], dst->chroma_size);
        pack_highbd_pic(src, src_buffer_highbd, ss_x, ss_y, EB_TRUE);
    } else {
        src_buffer_highbd[0] = (uint16_t *)src->buffer_y;
        src_buffer_highbd[1] = (uint16_t *)src->buffer_cb;
        src_buffer_highbd[2] = (uint16_t *)src->buffer_cr;
        dst_buffer_highbd[0] = (uint16_t *)dst->buffer_y;
        dst_buffer_highbd[1] = (uint16_t *)dst->buffer_cb;
        dst_buffer_highbd[2] = (uint16_t *)dst->buffer_cr;
    }
#if DEBUG_SCALING
    if (bd > 8)
        save_YUV_to_file_highbd("unscaled_pic_highbd.yuv",
                                src_buffer_highbd[0],
                                src_buffer_highbd[1],
                                src_buffer_highbd[2],
                                src->width + src->origin_x * 2,
                                src->height + src->origin_y * 2,
                                src->stride_y,
                                src->stride_cb,
                                src->stride_cr,
                                0,
                                0,
                                1,
                                1);
    else
        save_YUV_to_file("unscaled_pic.yuv",
                         src->buffer_y,
                         src->buffer_cb,
                         src->buffer_cr,
                         src->width + src->origin_x * 2,
                         src->height + src->origin_y * 2,
                         src->stride_y,
                         src->stride_cb,
                         src->stride_cr,
                         0,
                         0,
                         1,
                         1);
#endif

    for (int plane = 0; plane <= AOMMIN(num_planes, MAX_MB_PLANE - 1); ++plane) {
        if (bd > 8) {
            Av1HighbdResizePlane resize_plane_func = (src->height == dst->height)
                ? av1_highbd_resize_plane_horizontal
                : av1_highbd_resize_plane;
            switch (plane) {
            case 0:
                if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) && src_buffer_highbd[0] &&
                    dst_buffer_highbd[0])
                    resize_plane_func(
                        src_buffer_highbd[0] + src->origin_y * src->stride_y + src->origin_x,
                        src->height,
                        src->width,
                        src->stride_y,
                        dst_buffer_highbd[0] + dst->origin_y * dst->stride_y + dst->origin_x,
                        dst->height,
                        dst->width,
                        dst->stride_y,
                        bd);
                break;
            case 1:
                if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) && src_buffer_highbd[1] &&
                    dst_buffer_highbd[1])
                    resize_plane_func(
                        src_buffer_highbd[1] + (src->origin_y >> ss_y) * src->stride_cb +
                            (src->origin_x >> ss_x),
                        src->height >> ss_y,
                        src->width >> ss_x,
                        src->stride_cb,
                        dst_buffer_highbd[1] + (dst->origin_y >> ss_y) * dst->stride_cb +
                            (dst->origin_x >> ss_x),
                        (dst->height + ss_y) >> ss_y,
                        (dst->width + ss_x) >> ss_x,
                        dst->stride_cb,
                        bd);
                break;
            case 2:
                if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) && src_buffer_highbd[2] &&
                    dst_buffer_highbd[2])
                    resize_plane_func(
                        src_buffer_highbd[2] + (src->origin_y >> ss_y) * src->stride_cr +
                            (src->origin_x >> ss_x),
                        src->height >> ss_y,
                        src->width >> ss_x,
                        src->stride_cr,
                        dst_buffer_highbd[2] + (dst->origin_y >> ss_y) * dst->stride_cr +
                            (dst->origin_x >> ss_x),
                        (dst->height + ss_y) >> ss_y,
                        (dst->width + ss_x) >> ss_x,
                        dst->stride_cr,
                        bd);
                break;
            default: break;
            }
        } else {
            Av1ResizePlane resize_plane_func = (src->height == dst->height)
                ? av1_resize_plane_horizontal
                : av1_resize_plane;
            switch (plane) {
            case 0:
                if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) && src->buffer_y &&
                    dst->buffer_y)
                    resize_plane_func(src->buffer_y + src->origin_y * src->stride_y + src->origin_x,
                                      src->height,
                                      src->width,
                                      src->stride_y,
                                      dst->buffer_y + dst->origin_y * dst->stride_y + dst->origin_x,
                                      dst->height,
                                      dst->width,
                                      dst->stride_y);
                break;
            case 1:
                if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) && src->buffer_cb &&
                    dst->buffer_cb)
                    resize_plane_func(src->buffer_cb + (src->origin_y >> ss_y) * src->stride_cb +
                                          (src->origin_x >> ss_x),
                                      src->height >> ss_y,
                                      src->width >> ss_x,
                                      src->stride_cb,
                                      dst->buffer_cb + (dst->origin_y >> ss_y) * dst->stride_cb +
                                          (dst->origin_x >> ss_x),
                                      (dst->height + ss_y) >> ss_y,
                                      (dst->width + ss_x) >> ss_x,
                                      dst->stride_cb);
                break;
            case 2:
                if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) && src->buffer_cr &&
                    dst->buffer_cr)
                    resize_plane_func(src->buffer_cr + (src->origin_y >> ss_y) * src->stride_cr +
                                          (src->origin_x >> ss_x),
                                      src->height >> ss_y,
                                      src->width >> ss_x,
                                      src->stride_cr,
                                      dst->buffer_cr + (dst->origin_y >> ss_y) * dst->stride_cr +
                                          (dst->origin_x >> ss_x),
                                      (dst->height + ss_y) >> ss_y,
                                      (dst->width + ss_x) >> ss_x,
                                      dst->stride_cr);
                break;
            default: break;
            }
        }
    }

    // padding before unpack to support 10-bit with 2b compressed format
    if (bd > 8) {
        if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) && dst_buffer_highbd[0])
            generate_padding16_bit(dst_buffer_highbd[0],
                                   dst->stride_y,
                                   dst->width,
                                   dst->height,
                                   dst->origin_x,
                                   dst->origin_y);
        if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) && dst_buffer_highbd[1])
            generate_padding16_bit(dst_buffer_highbd[1],
                                   dst->stride_cb,
                                   (dst->width + ss_x) >> ss_x,
                                   (dst->height + ss_y) >> ss_y,
                                   (dst->origin_x + ss_x) >> ss_x,
                                   (dst->origin_y + ss_y) >> ss_y);
        if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) && dst_buffer_highbd[2])
            generate_padding16_bit(dst_buffer_highbd[2],
                                   dst->stride_cb,
                                   (dst->width + ss_x) >> ss_x,
                                   (dst->height + ss_y) >> ss_y,
                                   (dst->origin_x + ss_x) >> ss_x,
                                   (dst->origin_y + ss_y) >> ss_y);
    } else {
        if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Y_FLAG) && dst->buffer_y)
            generate_padding(dst->buffer_y,
                             dst->stride_y,
                             dst->width,
                             dst->height,
                             dst->origin_x,
                             dst->origin_y);
        if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Cb_FLAG) && dst->buffer_cb)
            generate_padding(dst->buffer_cb,
                             dst->stride_cb,
                             (dst->width + ss_x) >> ss_x,
                             (dst->height + ss_y) >> ss_y,
                             (dst->origin_x + ss_x) >> ss_x,
                             (dst->origin_y + ss_y) >> ss_y);
        if ((buffer_enable_mask & PICTURE_BUFFER_DESC_Cr_FLAG) && dst->buffer_cr)
            generate_padding(dst->buffer_cr,
                             dst->stride_cr,
                             (dst->width + ss_x) >> ss_x,
                             (dst->height + ss_y) >> ss_y,
                             (dst->origin_x + ss_x) >> ss_x,
                             (dst->origin_y + ss_y) >> ss_y);
    }

#if DEBUG_SCALING
    if (bd > 8)
        save_YUV_to_file_highbd("scaled_pic_highbd.yuv",
                                dst_buffer_highbd[0],
                                dst_buffer_highbd[1],
                                dst_buffer_highbd[2],
                                dst->width + dst->origin_x * 2,
                                dst->height + dst->origin_y * 2,
                                dst->stride_y,
                                dst->stride_cb,
                                dst->stride_cr,
                                0,
                                0,
                                1,
                                1);
    else
        save_YUV_to_file("scaled_pic.yuv",
                         dst->buffer_y,
                         dst->buffer_cb,
                         dst->buffer_cr,
                         dst->width + dst->origin_x * 2,
                         dst->height + dst->origin_y * 2,
                         dst->stride_y,
                         dst->stride_cb,
                         dst->stride_cr,
                         0,
                         0,
                         1,
                         1);
#endif
    if (bd > 8 && !is_packed) {
        unpack_highbd_pic(dst_buffer_highbd, dst, ss_x, ss_y, EB_TRUE);

        EB_FREE(src_buffer_highbd[0]);
        EB_FREE(src_buffer_highbd[1]);
        EB_FREE(src_buffer_highbd[2]);
        EB_FREE(dst_buffer_highbd[0]);
        EB_FREE(dst_buffer_highbd[1]);
        EB_FREE(dst_buffer_highbd[2]);
    }

    return EB_ErrorNone;
}

// Generate a random number in the range [0, 32768).
static INLINE unsigned int lcg_rand16(unsigned int *state) {
    *state = (unsigned int)(*state * 1103515245ULL + 12345);
    return *state / 65536 % 32768;
}

// Compute the horizontal frequency components' energy in a frame
// by calculuating the 16x4 Horizontal DCT. This is to be used to
// decide the superresolution parameters.
static void analyze_hor_freq(PictureParentControlSet *pcs_ptr, double *energy) {
    uint64_t freq_energy[16] = {0};

    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_picture_ptr;
    uint8_t             *in = input_picture_ptr->buffer_y + input_picture_ptr->origin_x +
        input_picture_ptr->origin_y * input_picture_ptr->stride_y;
    const int width  = input_picture_ptr->width;
    const int height = input_picture_ptr->height;

    DECLARE_ALIGNED(64, int32_t, coeff[16 * 4]);
    int n = 0;
    memset(freq_energy, 0, sizeof(freq_energy));

    // All treated as 8-bit input. For 10-bit input, uses high 8 bit (pointed by input_picture_ptr->buffer_y)
    DECLARE_ALIGNED(64, int16_t, src16[16 * 4]);
    for (int i = 0; i < height - 4; i += 4) {
        for (int j = 0; j < width - 16; j += 16) {
            for (int ii = 0; ii < 4; ++ii)
                for (int jj = 0; jj < 16; ++jj)
                    src16[ii * 16 + jj] = in[(i + ii) * input_picture_ptr->stride_y + (j + jj)];
            svt_av1_fwd_txfm2d_16x4(src16, coeff, 16, H_DCT, EB_8BIT);
            for (int k = 1; k < 16; ++k) {
                const uint64_t this_energy = ((int64_t)coeff[k] * coeff[k]) +
                    ((int64_t)coeff[k + 16] * coeff[k + 16]) +
                    ((int64_t)coeff[k + 32] * coeff[k + 32]) +
                    ((int64_t)coeff[k + 48] * coeff[k + 48]);
                freq_energy[k] += ROUND_POWER_OF_TWO(this_energy, 2);
            }
            n++;
        }
    }

    if (n) {
        for (int k = 1; k < 16; ++k) energy[k] = (double)freq_energy[k] / n;
        // Convert to cumulative energy
        for (int k = 14; k > 0; --k) energy[k] += energy[k + 1];
    } else {
        for (int k = 1; k < 16; ++k) energy[k] = 1e+20;
    }
}

#define SUPERRES_ENERGY_BY_Q2_THRESH_KEYFRAME_SOLO 0.012
#define SUPERRES_ENERGY_BY_Q2_THRESH_KEYFRAME 0.008
#define SUPERRES_ENERGY_BY_Q2_THRESH_ARFFRAME 0.008
#define SUPERRES_ENERGY_BY_AC_THRESH 0.2

static double get_energy_by_q2_thresh(const RATE_CONTROL *rc, int frame_update_type) {
    // TODO(now): Return keyframe thresh * factor based on frame type / pyramid
    // level.
    if (frame_update_type == ARF_UPDATE) {
        return SUPERRES_ENERGY_BY_Q2_THRESH_ARFFRAME;
    } else if (frame_update_type == KF_UPDATE) {
        if (rc->frames_to_key <= 1)
            return SUPERRES_ENERGY_BY_Q2_THRESH_KEYFRAME_SOLO;
        else
            return SUPERRES_ENERGY_BY_Q2_THRESH_KEYFRAME;
    } else {
        assert(0);
    }
    return 0;
}

static int av1_superres_in_recode_allowed(SequenceControlSet *scs_ptr) {
    const EbSvtAv1EncConfiguration *const static_config = &scs_ptr->static_config;
    // Empirically found to not be beneficial for image coding.
    return static_config->superres_mode == SUPERRES_AUTO &&
        static_config->superres_auto_search_type != SUPERRES_AUTO_SOLO /* &&
        scs_ptr->encode_context_ptr->rc.frames_to_key > 1*/
        ;
}

static uint8_t get_superres_denom_from_qindex_energy(int qindex, double *energy, double threshq,
                                                     double threshp) {
    const double q      = svt_av1_convert_qindex_to_q(qindex, AOM_BITS_8);
    const double tq     = threshq * q * q;
    const double tp     = threshp * energy[1];
    const double thresh = AOMMIN(tq, tp);
    int          k;
    for (k = SCALE_NUMERATOR * 2; k > SCALE_NUMERATOR; --k) {
        if (energy[k - 1] > thresh)
            break;
    }
    return 3 * SCALE_NUMERATOR - k;
}

int32_t get_frame_update_type(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {
    // Reasons why not use gf_group->update_type:
    //   1. It is valid only in 2nd pass of 2-pass encoding or lap_enabled is true. E.g. It's invalid in 1-pass CQP mode.
    //   2. It is set in RC process, so can't use it in processes before RC.
    if (pcs_ptr->frm_hdr.frame_type == KEY_FRAME) {
        return KF_UPDATE;
    }

    if (scs_ptr->max_temporal_layers > 0) {
        if (pcs_ptr->temporal_layer_index == 0) {
            return ARF_UPDATE;
        } else if (pcs_ptr->temporal_layer_index == pcs_ptr->hierarchical_levels) {
            return LF_UPDATE;
        } else {
            return INTNL_ARF_UPDATE;
        }
    } else {
        return LF_UPDATE;
    }
}

static uint8_t get_superres_denom_for_qindex(SequenceControlSet      *scs_ptr,
                                             PictureParentControlSet *pcs_ptr, int qindex,
                                             int sr_kf, int sr_arf) {
    // Use superres for Key-frames and Alt-ref frames only.
    int32_t update_type = get_frame_update_type(scs_ptr, pcs_ptr);
    if (update_type != KF_UPDATE && update_type != ARF_UPDATE) {
        return SCALE_NUMERATOR;
    }
    if (update_type == KF_UPDATE && !sr_kf) {
        return SCALE_NUMERATOR;
    }
    if (update_type == ARF_UPDATE && !sr_arf) {
        return SCALE_NUMERATOR;
    }

    double energy[16];
    analyze_hor_freq(pcs_ptr, energy);

    const double energy_by_q2_thresh = get_energy_by_q2_thresh(&scs_ptr->encode_context_ptr->rc,
                                                               update_type);
    int          denom               = get_superres_denom_from_qindex_energy(
        qindex, energy, energy_by_q2_thresh, SUPERRES_ENERGY_BY_AC_THRESH);

#if DEBUG_SUPERRES_ENERGY
    printf("\nFrame %d. energy_by_q2_thresh %.3f, energy = [",
           (int)pcs_ptr->picture_number,
           energy_by_q2_thresh);
    for (int k = 1; k < 16; ++k) printf("%f, ", energy[k]);
    printf("]\n");
    printf("boost = %d\n",
           (update_type == KF_UPDATE) ? scs_ptr->encode_context_ptr->rc.kf_boost
                                      : scs_ptr->encode_context_ptr->rc.gfu_boost);
    printf("denom = %d\n", denom);
#endif

    if (av1_superres_in_recode_allowed(scs_ptr)) {
        assert(scs_ptr->static_config.superres_mode != SUPERRES_NONE);
        // Force superres to be tried in the recode loop, as full-res is also going
        // to be tried anyway.
        denom = AOMMAX(denom, SCALE_NUMERATOR + 1);
    }
    return denom;
}

/*
 * Given the superres configurations and the frame type, determine the denominator and
 * encoding resolution
 */
static void calc_superres_params(superres_params_type *spr_params, SequenceControlSet *scs_ptr,
                                 PictureParentControlSet *pcs_ptr) {
    pcs_ptr->superres_total_recode_loop = 0;
    pcs_ptr->superres_recode_loop       = 0;
    spr_params->superres_denom          = SCALE_NUMERATOR;
    static unsigned int seed            = 34567;
    FrameHeader        *frm_hdr         = &pcs_ptr->frm_hdr;

    uint8_t superres_mode      = scs_ptr->static_config.superres_mode;
    uint8_t cfg_denom          = scs_ptr->static_config.superres_denom;
    uint8_t cfg_kf_denom       = scs_ptr->static_config.superres_kf_denom;
    uint8_t superres_qthres    = quantizer_to_qindex[scs_ptr->static_config.superres_qthres];
    uint8_t superres_kf_qthres = quantizer_to_qindex[scs_ptr->static_config.superres_kf_qthres];

    // super-res can only be enabled in case allow_intrabc is disabled
    if (frm_hdr->allow_intrabc)
        return;

    switch (superres_mode) {
    case SUPERRES_NONE: spr_params->superres_denom = SCALE_NUMERATOR; break;
    case SUPERRES_FIXED:
        if (frm_hdr->frame_type == KEY_FRAME)
            spr_params->superres_denom = cfg_kf_denom;
        else
            spr_params->superres_denom = cfg_denom;
        break;
    case SUPERRES_RANDOM: spr_params->superres_denom = (uint8_t)(lcg_rand16(&seed) % 9 + 8); break;
    case SUPERRES_QTHRESH: {
        // Do not use superres when screen content tools are used.
        if (frm_hdr->allow_screen_content_tools)
            break;

        const int q       = quantizer_to_qindex[pcs_ptr->picture_qp];
        const int qthresh = frame_is_intra_only(pcs_ptr) ? superres_kf_qthres : superres_qthres;
        if (q <= qthresh) {
            spr_params->superres_denom = SCALE_NUMERATOR;
        } else {
            spr_params->superres_denom = get_superres_denom_for_qindex(scs_ptr, pcs_ptr, q, 1, 1);
        }
        break;
    }
    case SUPERRES_AUTO: {
        const int q = quantizer_to_qindex[pcs_ptr->picture_qp];

        const int sr_search_type = scs_ptr->static_config.superres_auto_search_type;
        const int qthresh        = (sr_search_type == SUPERRES_AUTO_SOLO) ? 128 : 0;
        if (q <= qthresh) {
            spr_params->superres_denom = SCALE_NUMERATOR; // Don't use superres.
        } else {
            if (sr_search_type == SUPERRES_AUTO_SOLO) {
                spr_params->superres_denom = get_superres_denom_for_qindex(
                    scs_ptr, pcs_ptr, q, 1, 1);
            } else if (sr_search_type == SUPERRES_AUTO_DUAL) {
                pcs_ptr->superres_denom_array[0] = get_superres_denom_for_qindex(
                    scs_ptr, pcs_ptr, q, 1, 1);
                if (pcs_ptr->superres_denom_array[0] != SCALE_NUMERATOR) {
                    pcs_ptr->superres_denom_array[1]    = SCALE_NUMERATOR;
                    spr_params->superres_denom          = pcs_ptr->superres_denom_array[0];
                    pcs_ptr->superres_total_recode_loop = 2;
                }
            } else { // SUPERRES_AUTO_ALL
                assert(sr_search_type == SUPERRES_AUTO_ALL);
                int32_t update_type = get_frame_update_type(scs_ptr, pcs_ptr);
                if (update_type == KF_UPDATE || update_type == ARF_UPDATE) {
                    for (int i = 0; i < SCALE_NUMERATOR + 1; i++) {
                        if (i < SCALE_NUMERATOR) {
                            pcs_ptr->superres_denom_array[i] = SCALE_NUMERATOR + 1 + i;
                        } else {
                            pcs_ptr->superres_denom_array[i] = SCALE_NUMERATOR;
                        }
                    }
                    spr_params->superres_denom          = pcs_ptr->superres_denom_array[0];
                    pcs_ptr->superres_total_recode_loop = SCALE_NUMERATOR + 1;
                }
            }
        }
        break;
    }
    default: break;
    }

    // only encoding width is adjusted
    calculate_scaled_size_helper(&spr_params->encoding_width, spr_params->superres_denom);
}

static EbErrorType downscaled_source_buffer_desc_ctor(
    EbPictureBufferDesc **picture_ptr, EbPictureBufferDesc *picture_ptr_for_reference,
    superres_params_type spr_params) {
    EbPictureBufferDescInitData initData;

    initData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    initData.max_width          = spr_params.encoding_width;
    initData.max_height         = spr_params.encoding_height;
    initData.bit_depth          = picture_ptr_for_reference->bit_depth;
    initData.color_format       = picture_ptr_for_reference->color_format;
    initData.split_mode    = (picture_ptr_for_reference->bit_depth > EB_8BIT) ? EB_TRUE : EB_FALSE;
    initData.left_padding  = picture_ptr_for_reference->origin_x;
    initData.right_padding = picture_ptr_for_reference->origin_x;
    initData.top_padding   = picture_ptr_for_reference->origin_y;
    initData.bot_padding   = picture_ptr_for_reference->origin_bot_y;

    EB_NEW(*picture_ptr, svt_picture_buffer_desc_ctor, (EbPtr)&initData);

    return EB_ErrorNone;
}

EbErrorType sb_geom_init_pcs(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr);

EbErrorType sb_params_init_pcs(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr);

static uint8_t get_denom_idx(uint8_t superres_denom) {
    uint8_t denom_idx = (uint8_t)(superres_denom - SCALE_NUMERATOR - 1);
    return denom_idx;
}

/*
 * Modify encoder parameters and structures that depend on picture resolution
 * Performed after a source picture has been scaled
 */
void scale_pcs_params(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
                      superres_params_type spr_params, uint16_t source_width,
                      uint16_t source_height) {
    Av1Common *cm = pcs_ptr->av1_cm;

    // frame sizes
    cm->frm_size.frame_width          = spr_params.encoding_width;
    cm->frm_size.frame_height         = spr_params.encoding_height;
    cm->frm_size.render_width         = source_width;
    cm->frm_size.render_height        = source_height;
    cm->frm_size.superres_denominator = spr_params.superres_denom;

    // align width and height to be a multiple of 8
    uint16_t aligned_width  = (uint16_t)ALIGN_POWER_OF_TWO(spr_params.encoding_width, 3);
    uint16_t aligned_height = (uint16_t)ALIGN_POWER_OF_TWO(spr_params.encoding_height, 3);

    // remove this assertion because decoder specify allowing the width can be multiple of 8
    // encoder should fix the issue if downscaled pictures are not 8 pixel aligned
    //assert((aligned_width == spr_params.encoding_width) &&
    //       "Downscaled width needs to be a multiple of 8 "
    //       "(otherwise not yet implemented)");

    // change frame width and height params in pcs
    pcs_ptr->frame_width  = spr_params.encoding_width;
    pcs_ptr->frame_height = spr_params.encoding_height;

    pcs_ptr->aligned_width  = aligned_width;
    pcs_ptr->aligned_height = aligned_height;

    // number of SBs
    const uint16_t picture_sb_width  = (uint16_t)((aligned_width + scs_ptr->sb_size_pix - 1) /
                                                 scs_ptr->sb_size_pix);
    const uint16_t picture_sb_height = (uint16_t)((aligned_height + scs_ptr->sb_size_pix - 1) /
                                                  scs_ptr->sb_size_pix);

    pcs_ptr->picture_sb_width  = picture_sb_width; // TODO: use this instead of re-computing
    pcs_ptr->picture_sb_height = picture_sb_height;

    pcs_ptr->sb_total_count = picture_sb_width * picture_sb_height;

    // mi params
    cm->mi_stride = picture_sb_width * (scs_ptr->sb_size_pix >> MI_SIZE_LOG2);
    cm->mi_cols   = aligned_width >> MI_SIZE_LOG2;
    cm->mi_rows   = aligned_height >> MI_SIZE_LOG2;

    derive_input_resolution(&pcs_ptr->input_resolution,
                            spr_params.encoding_width * spr_params.encoding_height);

    // create new picture level sb_params and sb_geom
    sb_params_init_pcs(scs_ptr, pcs_ptr);
    sb_geom_init_pcs(scs_ptr, pcs_ptr);

    if (pcs_ptr->frame_superres_enabled == EB_TRUE) {
        pcs_ptr->frm_hdr.use_ref_frame_mvs = 0;
    } else {
        if (pcs_ptr->slice_type == I_SLICE)
            pcs_ptr->frm_hdr.use_ref_frame_mvs = 0;
        else
            pcs_ptr->frm_hdr.use_ref_frame_mvs = scs_ptr->mfmv_enabled;
    }
}

/*
 * Memory allocation for donwscaled reconstructed reference pictures
 */
static EbErrorType allocate_downscaled_reference_pics(
    EbPictureBufferDesc **downscaled_reference_picture_ptr,
    EbPictureBufferDesc *picture_ptr_for_reference, PictureParentControlSet *pcs_ptr) {
    EbPictureBufferDescInitData ref_pic_buf_desc_init_data;

    // Initialize the various Picture types
    ref_pic_buf_desc_init_data.max_width    = pcs_ptr->frame_width; // aligned_width;
    ref_pic_buf_desc_init_data.max_height   = pcs_ptr->frame_height; // aligned_height;
    ref_pic_buf_desc_init_data.bit_depth    = pcs_ptr->scs_ptr->encoder_bit_depth;
    ref_pic_buf_desc_init_data.color_format = pcs_ptr->scs_ptr->static_config.encoder_color_format;
    ref_pic_buf_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;

    ref_pic_buf_desc_init_data.left_padding  = picture_ptr_for_reference->origin_x;
    ref_pic_buf_desc_init_data.right_padding = picture_ptr_for_reference->origin_x;
    ref_pic_buf_desc_init_data.top_padding   = picture_ptr_for_reference->origin_y;
    ref_pic_buf_desc_init_data.bot_padding   = picture_ptr_for_reference->origin_bot_y;
    ref_pic_buf_desc_init_data.mfmv          = pcs_ptr->scs_ptr->mfmv_enabled;

    //TODO:12bit
    if (ref_pic_buf_desc_init_data.bit_depth == EB_10BIT) {
        // Hsan: set split_mode to 0 to construct the packed reference buffer (used @ EP)

        // align with svt_reference_object_ctor()
        // Use 10bit here to use in MD
        ref_pic_buf_desc_init_data.split_mode = EB_TRUE;
        ref_pic_buf_desc_init_data.bit_depth  = EB_10BIT;
        EB_NEW(*downscaled_reference_picture_ptr,
               svt_picture_buffer_desc_ctor,
               (EbPtr)&ref_pic_buf_desc_init_data);
    } else {
        // Hsan: set split_mode to 0 to as 8BIT input
        ref_pic_buf_desc_init_data.split_mode = EB_FALSE;
        EB_NEW(*downscaled_reference_picture_ptr,
               svt_picture_buffer_desc_ctor,
               (EbPtr)&ref_pic_buf_desc_init_data);
    }

    return EB_ErrorNone;
}

/*
 * Memory allocation for donwscaled source reference pictures
 */
static EbErrorType allocate_downscaled_source_reference_pics(
    EbPictureBufferDesc **input_padded_picture_ptr,
    EbPictureBufferDesc **quarter_downsampled_picture_ptr,
    EbPictureBufferDesc **sixteenth_downsampled_picture_ptr,
    EbPictureBufferDesc *picture_ptr_for_reference, superres_params_type spr_params) {
    EbPictureBufferDescInitData initData;

    initData.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
    initData.max_width          = spr_params.encoding_width;
    initData.max_height         = spr_params.encoding_height;
    initData.bit_depth          = picture_ptr_for_reference->bit_depth;
    initData.color_format       = picture_ptr_for_reference->color_format;
    initData.split_mode    = (picture_ptr_for_reference->bit_depth > EB_8BIT) ? EB_TRUE : EB_FALSE;
    initData.left_padding  = picture_ptr_for_reference->origin_x;
    initData.right_padding = picture_ptr_for_reference->origin_x;
    initData.top_padding   = picture_ptr_for_reference->origin_y;
    initData.bot_padding   = picture_ptr_for_reference->origin_bot_y;

    EB_NEW(*input_padded_picture_ptr, svt_picture_buffer_desc_ctor, (EbPtr)&initData);

    initData.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
    initData.max_width          = spr_params.encoding_width >> 1;
    initData.max_height         = spr_params.encoding_height >> 1;
    initData.bit_depth          = picture_ptr_for_reference->bit_depth;
    initData.color_format       = picture_ptr_for_reference->color_format;
    initData.split_mode    = (picture_ptr_for_reference->bit_depth > EB_8BIT) ? EB_TRUE : EB_FALSE;
    initData.left_padding  = picture_ptr_for_reference->origin_x >> 1;
    initData.right_padding = picture_ptr_for_reference->origin_x >> 1;
    initData.top_padding   = picture_ptr_for_reference->origin_y >> 1;
    initData.bot_padding   = picture_ptr_for_reference->origin_bot_y >> 1;

    EB_NEW(*quarter_downsampled_picture_ptr, svt_picture_buffer_desc_ctor, (EbPtr)&initData);

    initData.buffer_enable_mask = PICTURE_BUFFER_DESC_LUMA_MASK;
    initData.max_width          = spr_params.encoding_width >> 2;
    initData.max_height         = spr_params.encoding_height >> 2;
    initData.bit_depth          = picture_ptr_for_reference->bit_depth;
    initData.color_format       = picture_ptr_for_reference->color_format;
    initData.split_mode    = (picture_ptr_for_reference->bit_depth > EB_8BIT) ? EB_TRUE : EB_FALSE;
    initData.left_padding  = picture_ptr_for_reference->origin_x >> 2;
    initData.right_padding = picture_ptr_for_reference->origin_x >> 2;
    initData.top_padding   = picture_ptr_for_reference->origin_y >> 2;
    initData.bot_padding   = picture_ptr_for_reference->origin_bot_y >> 2;

    EB_NEW(*sixteenth_downsampled_picture_ptr, svt_picture_buffer_desc_ctor, (EbPtr)&initData);
    return EB_ErrorNone;
}

/*
 * Allocate memory and perform scaling of the source reference picture (references to the current picture)
 * and its decimated/filtered versions to match with the input picture resolution
 * This is used in the open-loop stage.
 */
void scale_source_references(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
                             EbPictureBufferDesc *input_picture_ptr) {
    EbPaReferenceObject *reference_object;

    uint8_t        denom_idx  = get_denom_idx(pcs_ptr->superres_denom);
    const int32_t  num_planes = 0; // Y only
    const uint32_t ss_x       = scs_ptr->subsampling_x;
    const uint32_t ss_y       = scs_ptr->subsampling_y;
    uint32_t       num_of_list_to_search =
        (pcs_ptr->slice_type == P_SLICE) ? 1 /*List 0 only*/ : 2 /*List 0 + 1*/;

    for (uint8_t list_index = REF_LIST_0; list_index < num_of_list_to_search; ++list_index) {
        uint8_t ref_pic_index;
        uint8_t num_of_ref_pic_to_search = (pcs_ptr->slice_type == P_SLICE)
            ? pcs_ptr->ref_list0_count
            : (list_index == REF_LIST_0) ? pcs_ptr->ref_list0_count
                                         : pcs_ptr->ref_list1_count;
        for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            reference_object = (EbPaReferenceObject *)pcs_ptr
                                   ->ref_pa_pic_ptr_array[list_index][ref_pic_index]
                                   ->object_ptr;

            uint64_t ref_picture_number = pcs_ptr->ref_pic_poc_array[list_index][ref_pic_index];

            EbPictureBufferDesc *ref_pic_ptr = reference_object->input_padded_picture_ptr;

            // if the size of the reference pic is different than the size of the input pic, then scale references
            if (ref_pic_ptr->width != input_picture_ptr->width) {
                EbBool do_resize = EB_FALSE;

                svt_block_on_mutex(reference_object->resize_mutex[denom_idx]);
                if (reference_object->downscaled_input_padded_picture_ptr[denom_idx] == NULL) {
                    superres_params_type spr_params = {pcs_ptr->frame_width, // aligned_width
                                                       pcs_ptr->frame_height, // aligned_height
                                                       scs_ptr->static_config.superres_mode};

                    // Allocate downsampled reference picture buffer descriptors
                    allocate_downscaled_source_reference_pics(
                        &reference_object->downscaled_input_padded_picture_ptr[denom_idx],
                        &reference_object->downscaled_quarter_downsampled_picture_ptr[denom_idx],
                        &reference_object->downscaled_sixteenth_downsampled_picture_ptr[denom_idx],
                        ref_pic_ptr,
                        spr_params);

                    do_resize = EB_TRUE;
                }

                if (reference_object->downscaled_picture_number[denom_idx] != ref_picture_number) {
                    do_resize = EB_TRUE;
                }

                if (do_resize) {
                    EbPictureBufferDesc *down_ref_pic_ptr =
                        reference_object->downscaled_input_padded_picture_ptr[denom_idx];

                    // downsample input padded picture buffer
                    av1_resize_frame(ref_pic_ptr,
                                     down_ref_pic_ptr,
                                     8, // only 8-bit buffer needed for open-loop processing
                                     num_planes,
                                     ss_x,
                                     ss_y,
                                     0, // is_packed
                                     PICTURE_BUFFER_DESC_LUMA_MASK); // buffer_enable_mask

                    // 1/4 & 1/16 input picture downsampling
                    if (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
                        downsample_filtering_input_picture(
                            pcs_ptr,
                            down_ref_pic_ptr,
                            reference_object->downscaled_quarter_downsampled_picture_ptr[denom_idx],
                            reference_object
                                ->downscaled_sixteenth_downsampled_picture_ptr[denom_idx]);
                    } else {
                        downsample_decimation_input_picture(
                            pcs_ptr,
                            down_ref_pic_ptr,
                            reference_object->downscaled_quarter_downsampled_picture_ptr[denom_idx],
                            reference_object
                                ->downscaled_sixteenth_downsampled_picture_ptr[denom_idx]);
                    }

                    reference_object->downscaled_picture_number[denom_idx] = ref_picture_number;
                }

                svt_release_mutex(reference_object->resize_mutex[denom_idx]);
            }
        }
    }
}

/*
 * Allocate memory and perform scaling of the input reference picture (current picture)
 * and its decimated/filtered versions to match with the input picture resolution
 * This is used in the open-loop stage.
 */
static void scale_input_references(PictureParentControlSet *pcs_ptr,
                                   superres_params_type     superres_params) {
    uint8_t denom_idx = get_denom_idx(superres_params.superres_denom);

    // reference structures (padded pictures + downsampled versions)
    EbPaReferenceObject *src_object = (EbPaReferenceObject *)
                                          pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
    EbPictureBufferDesc *padded_pic_ptr = src_object->input_padded_picture_ptr;
    EbBool               do_resize      = EB_FALSE;

    svt_block_on_mutex(src_object->resize_mutex[denom_idx]);

    if (src_object->downscaled_input_padded_picture_ptr[denom_idx] == NULL) {
        // Allocate downsampled reference picture buffer descriptors
        allocate_downscaled_source_reference_pics(
            &src_object->downscaled_input_padded_picture_ptr[denom_idx],
            &src_object->downscaled_quarter_downsampled_picture_ptr[denom_idx],
            &src_object->downscaled_sixteenth_downsampled_picture_ptr[denom_idx],
            padded_pic_ptr,
            superres_params);
        do_resize = EB_TRUE;
    }

    if (src_object->downscaled_picture_number[denom_idx] != pcs_ptr->picture_number) {
        do_resize = EB_TRUE;
    }

    if (do_resize) {
        padded_pic_ptr = src_object->downscaled_input_padded_picture_ptr[denom_idx];
        EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_picture_ptr;

        for (uint32_t row = 0;
             row < (uint32_t)(input_picture_ptr->height + 2 * input_picture_ptr->origin_y);
             row++)
            EB_MEMCPY(padded_pic_ptr->buffer_y + row * padded_pic_ptr->stride_y,
                      input_picture_ptr->buffer_y + row * input_picture_ptr->stride_y,
                      sizeof(uint8_t) * input_picture_ptr->stride_y);

        // 1/4 & 1/16 downsampled input picture
        if (pcs_ptr->scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
            downsample_filtering_input_picture(
                pcs_ptr,
                padded_pic_ptr,
                src_object->downscaled_quarter_downsampled_picture_ptr[denom_idx],
                src_object->downscaled_sixteenth_downsampled_picture_ptr[denom_idx]);
        } else {
            downsample_decimation_input_picture(
                pcs_ptr,
                padded_pic_ptr,
                src_object->downscaled_quarter_downsampled_picture_ptr[denom_idx],
                src_object->downscaled_sixteenth_downsampled_picture_ptr[denom_idx]);
        }

        src_object->downscaled_picture_number[denom_idx] = pcs_ptr->picture_number;
    }

    svt_release_mutex(src_object->resize_mutex[denom_idx]);
}

/*
 * Allocate memory and perform scaling of the reconstructed reference pictures
 * to match with the input picture resolution
 */
void scale_rec_references(PictureControlSet *pcs_ptr, EbPictureBufferDesc *input_picture_ptr,
                          uint8_t hbd_mode_decision) {
    EbReferenceObject *reference_object;

    PictureParentControlSet *ppcs_ptr = pcs_ptr->parent_pcs_ptr;
    SequenceControlSet      *scs_ptr  = ppcs_ptr->scs_ptr;

    uint8_t        denom_idx  = get_denom_idx(ppcs_ptr->superres_denom);
    const int32_t  num_planes = av1_num_planes(&scs_ptr->seq_header.color_config);
    const uint32_t ss_x       = scs_ptr->subsampling_x;
    const uint32_t ss_y       = scs_ptr->subsampling_y;
    uint32_t       num_of_list_to_search =
        (pcs_ptr->slice_type == P_SLICE) ? 1 /*List 0 only*/ : 2 /*List 0 + 1*/;

    for (uint8_t list_index = REF_LIST_0; list_index < num_of_list_to_search; ++list_index) {
        uint8_t ref_pic_index;
        uint8_t num_of_ref_pic_to_search = (ppcs_ptr->slice_type == P_SLICE)
            ? ppcs_ptr->ref_list0_count
            : (list_index == REF_LIST_0) ? ppcs_ptr->ref_list0_count
                                         : ppcs_ptr->ref_list1_count;
        for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            reference_object = (EbReferenceObject *)pcs_ptr
                                   ->ref_pic_ptr_array[list_index][ref_pic_index]
                                   ->object_ptr;

            uint64_t ref_picture_number = ppcs_ptr->ref_pic_poc_array[list_index][ref_pic_index];

            EbPictureBufferDesc *ref_pic_ptr = get_ref_pic_buffer(
                pcs_ptr, hbd_mode_decision, list_index, ref_pic_index);
            // if the size of the reference pic is different than the size of the input pic, then scale references
            if (ref_pic_ptr->width != input_picture_ptr->width) {
                EbBool do_resize = EB_FALSE;

                svt_block_on_mutex(reference_object->resize_mutex[denom_idx]);

                EbPictureBufferDesc *down_ref_pic8bit =
                    reference_object->downscaled_reference_picture[denom_idx];

                EbPictureBufferDesc *down_ref_pic_ptr = down_ref_pic8bit;

                if (down_ref_pic_ptr != NULL) {
                    if (ref_pic_ptr->bit_depth != down_ref_pic_ptr->bit_depth) {
                        EB_DELETE(reference_object->downscaled_reference_picture[denom_idx]);
                        down_ref_pic8bit = NULL;
                        down_ref_pic_ptr = NULL;
                    }
                }

                if (down_ref_pic_ptr == NULL) {
                    // Allocate downsampled reference picture buffer descriptors
                    allocate_downscaled_reference_pics(
                        &reference_object->downscaled_reference_picture[denom_idx],
                        ref_pic_ptr,
                        ppcs_ptr);

                    down_ref_pic8bit = reference_object->downscaled_reference_picture[denom_idx];
                    do_resize        = EB_TRUE;
                }

                if (reference_object->downscaled_picture_number[denom_idx] != ref_picture_number) {
                    do_resize = EB_TRUE;
                }

                if (do_resize) {
                    // downsample input padded picture buffer
                    av1_resize_frame(reference_object->reference_picture,
                                     down_ref_pic8bit,
                                     down_ref_pic8bit->bit_depth,
                                     num_planes,
                                     ss_x,
                                     ss_y,
                                     down_ref_pic8bit->packed_flag,
                                     PICTURE_BUFFER_DESC_FULL_MASK); // buffer_enable_mask

                    reference_object->downscaled_picture_number[denom_idx] = ref_picture_number;
                    //printf("rescaled reference picture %d\n", (int)ref_picture_number);
                }

                svt_release_mutex(reference_object->resize_mutex[denom_idx]);
            }
        }
    }
}

/*
 * Check if width of input and reference <reconstructed> pictures match.
 * if not, return pointers to downscaled reference picture of the correct resolution
 */
void use_scaled_rec_refs_if_needed(PictureControlSet   *pcs_ptr,
                                   EbPictureBufferDesc *input_picture_ptr,
                                   EbReferenceObject *ref_obj, EbPictureBufferDesc **ref_pic,
                                   uint8_t hbd_mode_decision) {
    if ((*ref_pic)->width != input_picture_ptr->width) {
        uint8_t denom_idx = get_denom_idx(pcs_ptr->parent_pcs_ptr->superres_denom);
        UNUSED(hbd_mode_decision);
        {
            assert(ref_obj->downscaled_reference_picture[denom_idx] != NULL);
            *ref_pic = ref_obj->downscaled_reference_picture[denom_idx];
        }
    }
    assert((*ref_pic)->width == input_picture_ptr->width);
}

/*
 * Check if width of input and reference <source> pictures match.
 * if not, return pointers to downscaled reference picture of the correct resolution.
 * These are used in the open-loop stage.
 */
void use_scaled_source_refs_if_needed(PictureParentControlSet *pcs_ptr,
                                      EbPictureBufferDesc     *input_picture_ptr,
                                      EbPaReferenceObject     *ref_obj,
                                      EbPictureBufferDesc    **ref_pic_ptr,
                                      EbPictureBufferDesc    **quarter_ref_pic_ptr,
                                      EbPictureBufferDesc    **sixteenth_ref_pic_ptr) {
    if ((*ref_pic_ptr)->width != input_picture_ptr->width) {
        uint8_t denom_idx = get_denom_idx(pcs_ptr->superres_denom);

        assert(ref_obj->downscaled_input_padded_picture_ptr[denom_idx] != NULL);

        *ref_pic_ptr           = ref_obj->downscaled_input_padded_picture_ptr[denom_idx];
        *quarter_ref_pic_ptr   = ref_obj->downscaled_quarter_downsampled_picture_ptr[denom_idx];
        *sixteenth_ref_pic_ptr = ref_obj->downscaled_sixteenth_downsampled_picture_ptr[denom_idx];
    }
    assert((*ref_pic_ptr)->width == input_picture_ptr->width);
}

void reset_resized_picture(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
                           EbPictureBufferDesc *input_picture_ptr) {
    superres_params_type spr_params = {input_picture_ptr->width, // encoding_width
                                       input_picture_ptr->height, // encoding_height
                                       SCALE_NUMERATOR};
    pcs_ptr->superres_denom         = spr_params.superres_denom;
    pcs_ptr->frame_superres_enabled = EB_FALSE;
    // restore frame size (width) which was changed by super-res tool
    scale_pcs_params(
        scs_ptr, pcs_ptr, spr_params, input_picture_ptr->width, input_picture_ptr->height);
    // delete picture buffer allocated by super-res tool
    // TODO: reuse the buffer if current picture's denominator is the same as previous one's.
    EB_DELETE(pcs_ptr->enhanced_downscaled_picture_ptr);
}

/*
 * If super-res is ON, determine super-res denominator for current picture,
 * perform resizing of source picture and
 * adjust resolution related parameters
 */
void init_resize_picture(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {
    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_unscaled_picture_ptr;

    superres_params_type spr_params = {input_picture_ptr->width, // encoding_width
                                       input_picture_ptr->height, // encoding_height
                                       scs_ptr->static_config.superres_denom};

    EbBool first_loop = !(scs_ptr->static_config.superres_mode == SUPERRES_AUTO &&
                          pcs_ptr->superres_recode_loop > 0);
    if (first_loop) { // first loop of multiple coding loop (auto-dual or auto-all mode) or the only loop (all the other modes)
        // determine super-res denom
        calc_superres_params(&spr_params, scs_ptr, pcs_ptr);
    } else {
        if (pcs_ptr->superres_recode_loop < pcs_ptr->superres_total_recode_loop) {
            spr_params.superres_denom =
                pcs_ptr->superres_denom_array[pcs_ptr->superres_recode_loop];
        } else { // extra loop to pick up a scaled recode
            // denom is set by downstream packetization process
            spr_params.superres_denom = pcs_ptr->superres_denom;
        }
        calculate_scaled_size_helper(&spr_params.encoding_width, spr_params.superres_denom);
    }

    // delete picture buffer allocated by superres tool
    // TODO: reuse the buffer if current picture's denom is the same as previous one's.
    EB_DELETE(pcs_ptr->enhanced_downscaled_picture_ptr);

    if (spr_params.superres_denom != SCALE_NUMERATOR) {
        pcs_ptr->superres_denom = spr_params.superres_denom;

        // Allocate downsampled picture buffer descriptor
        downscaled_source_buffer_desc_ctor(
            &pcs_ptr->enhanced_downscaled_picture_ptr, input_picture_ptr, spr_params);

        const int32_t  num_planes = av1_num_planes(&scs_ptr->seq_header.color_config);
        const uint32_t ss_x       = scs_ptr->subsampling_x;
        const uint32_t ss_y       = scs_ptr->subsampling_y;

        // downsample picture buffer
        assert(pcs_ptr->enhanced_downscaled_picture_ptr);
        av1_resize_frame(input_picture_ptr,
                         pcs_ptr->enhanced_downscaled_picture_ptr,
                         pcs_ptr->enhanced_downscaled_picture_ptr->bit_depth,
                         num_planes,
                         ss_x,
                         ss_y,
                         pcs_ptr->enhanced_downscaled_picture_ptr->packed_flag,
                         PICTURE_BUFFER_DESC_FULL_MASK); // buffer_enable_mask

        // use downscaled picture instead of original res for mode decision, encoding loop etc
        // after temporal filtering and motion estimation
        pcs_ptr->enhanced_picture_ptr = pcs_ptr->enhanced_downscaled_picture_ptr;

        pcs_ptr->frame_superres_enabled = EB_TRUE;

        scale_pcs_params(
            scs_ptr, pcs_ptr, spr_params, input_picture_ptr->width, input_picture_ptr->height);

        scale_input_references(pcs_ptr, spr_params);

        if (pcs_ptr->slice_type != I_SLICE) {
            scale_source_references(scs_ptr, pcs_ptr, pcs_ptr->enhanced_picture_ptr);
        }
    } else {
        // pcs_ptr previously might be used and dirty in params
        // clean up if current frame doesn't need scaling
        pcs_ptr->enhanced_picture_ptr = pcs_ptr->enhanced_unscaled_picture_ptr;
        reset_resized_picture(scs_ptr, pcs_ptr, input_picture_ptr);
    }
}
