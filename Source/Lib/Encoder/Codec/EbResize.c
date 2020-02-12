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

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "EbResize.h"

#define DEBUG_SCALING 0
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

void calculate_scaled_size_helper(uint16_t *dim, uint8_t denom);

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
    uint8_t *      optr = output;
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
    uint8_t *      optr = output;
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
    // TODO: use original filter in libaom
    if (out_length16 >= in_length * 16) return filteredinterp_filters1000;
    if (out_length16 >= in_length * 16)
        return filteredinterp_filters875; // wrong
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
    const int32_t delta =
        (((uint32_t)in_length << RS_SCALE_SUBPEL_BITS) + out_length / 2) / out_length;
    const int32_t offset =
        in_length > out_length
            ? (((int32_t)(in_length - out_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
                  out_length
            : -(((int32_t)(out_length - in_length) << (RS_SCALE_SUBPEL_BITS - 1)) +
                out_length / 2) /
                  out_length;
    uint8_t *optr = output;
    int      x, x1, x2, sum, k, int_pel, sub_pel;
    int32_t  y;

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
        memcpy(output, input, sizeof(output[0]) * length);
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
                out = (s & 1 ? otmp2 : otmp);
            if (filteredlength & 1)
                down2_symodd(in, filteredlength, out);
            else
                down2_symeven(in, filteredlength, out);
            filteredlength = proj_filteredlength;
        }
        if (filteredlength != olength) { interpolate(out, filteredlength, output, olength); }
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

EbErrorType av1_resize_plane(const uint8_t *const input, int height, int width, int in_stride,
                             uint8_t *output, int height2, int width2, int out_stride) {
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

static void highbd_interpolate_core(const uint16_t *const input, int in_length, uint16_t *output,
                                    int out_length, int bd, const int16_t *interp_filters,
                                    int interp_taps) {
    const int32_t delta =
        (((uint32_t)in_length << RS_SCALE_SUBPEL_BITS) + out_length / 2) / out_length;
    const int32_t offset =
        in_length > out_length
            ? (((int32_t)(in_length - out_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
                  out_length
            : -(((int32_t)(out_length - in_length) << (RS_SCALE_SUBPEL_BITS - 1)) +
                out_length / 2) /
                  out_length;
    uint16_t *optr = output;
    int       x, x1, x2, sum, k, int_pel, sub_pel;
    int32_t   y;

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
    uint16_t *            optr = output;
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
    uint16_t *            optr = output;
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
        memcpy(output, input, sizeof(output[0]) * length);
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
                out = (s & 1 ? otmp2 : otmp);
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

EbErrorType av1_highbd_resize_plane(const uint16_t *const input, int height, int width,
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

void pack_highbd_pic(const EbPictureBufferDesc *pic_ptr, uint16_t *buffer_16bit[3], uint32_t ss_x,
                     uint32_t ss_y, EbBool include_padding);

void unpack_highbd_pic(uint16_t *buffer_highbd[3], EbPictureBufferDesc *pic_ptr, uint32_t ss_x,
                       uint32_t ss_y, EbBool include_padding);

void save_YUV_to_file(char *filename, EbByte buffer_y, EbByte buffer_u, EbByte buffer_v,
                      uint16_t width, uint16_t height, uint16_t stride_y, uint16_t stride_u,
                      uint16_t stride_v, uint16_t origin_y, uint16_t origin_x, uint32_t ss_x,
                      uint32_t ss_y);

void save_YUV_to_file_highbd(char *filename, uint16_t *buffer_y, uint16_t *buffer_u,
                             uint16_t *buffer_v, uint16_t width, uint16_t height, uint16_t stride_y,
                             uint16_t stride_u, uint16_t stride_v, uint16_t origin_y,
                             uint16_t origin_x, uint32_t ss_x, uint32_t ss_y);

EbErrorType av1_resize_and_extend_frame(const EbPictureBufferDesc *src, EbPictureBufferDesc *dst,
                                        int bd, const int num_planes, const uint32_t ss_x,
                                        const uint32_t ss_y) {
    uint16_t *src_buffer_highbd[MAX_MB_PLANE];
    uint16_t *dst_buffer_highbd[MAX_MB_PLANE];

    if (bd > 8) {
        EB_MALLOC_ARRAY(src_buffer_highbd[0], src->luma_size);
        EB_MALLOC_ARRAY(src_buffer_highbd[1], src->chroma_size);
        EB_MALLOC_ARRAY(src_buffer_highbd[2], src->chroma_size);
        EB_MALLOC_ARRAY(dst_buffer_highbd[0], dst->luma_size);
        EB_MALLOC_ARRAY(dst_buffer_highbd[1], dst->chroma_size);
        EB_MALLOC_ARRAY(dst_buffer_highbd[2], dst->chroma_size);
        pack_highbd_pic(src, src_buffer_highbd, ss_x, ss_y, EB_TRUE);
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

    for (int plane = 0; plane < AOMMIN(num_planes, MAX_MB_PLANE); ++plane) {
        if (bd > 8) {
            switch (plane) {
            case 0:
                av1_highbd_resize_plane(
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
                av1_highbd_resize_plane(
                    src_buffer_highbd[1] + (src->origin_y >> ss_y) * src->stride_cb +
                        (src->origin_x >> ss_x),
                    src->height >> ss_y,
                    src->width >> ss_x,
                    src->stride_cb,
                    dst_buffer_highbd[1] + (dst->origin_y >> ss_y) * dst->stride_cb +
                        (dst->origin_x >> ss_x),
                    dst->height >> ss_y,
                    dst->width >> ss_x,
                    dst->stride_cb,
                    bd);
                break;
            case 2:
                av1_highbd_resize_plane(
                    src_buffer_highbd[2] + (src->origin_y >> ss_y) * src->stride_cr +
                        (src->origin_x >> ss_x),
                    src->height >> ss_y,
                    src->width >> ss_x,
                    src->stride_cr,
                    dst_buffer_highbd[2] + (dst->origin_y >> ss_y) * dst->stride_cr +
                        (dst->origin_x >> ss_x),
                    dst->height >> ss_y,
                    dst->width >> ss_x,
                    dst->stride_cr,
                    bd);
                break;
            default: break;
            }
        } else {
            switch (plane) {
            case 0:
                av1_resize_plane(src->buffer_y + src->origin_y * src->stride_y + src->origin_x,
                                 src->height,
                                 src->width,
                                 src->stride_y,
                                 dst->buffer_y + dst->origin_y * dst->stride_y + dst->origin_x,
                                 dst->height,
                                 dst->width,
                                 dst->stride_y);
                break;
            case 1:
                av1_resize_plane(src->buffer_cb + (src->origin_y >> ss_y) * src->stride_cb +
                                     (src->origin_x >> ss_x),
                                 src->height >> ss_y,
                                 src->width >> ss_x,
                                 src->stride_cb,
                                 dst->buffer_cb + (dst->origin_y >> ss_y) * dst->stride_cb +
                                     (dst->origin_x >> ss_x),
                                 dst->height >> ss_y,
                                 dst->width >> ss_x,
                                 dst->stride_cb);
                break;
            case 2:
                av1_resize_plane(src->buffer_cr + (src->origin_y >> ss_y) * src->stride_cr +
                                     (src->origin_x >> ss_x),
                                 src->height >> ss_y,
                                 src->width >> ss_x,
                                 src->stride_cr,
                                 dst->buffer_cr + (dst->origin_y >> ss_y) * dst->stride_cr +
                                     (dst->origin_x >> ss_x),
                                 dst->height >> ss_y,
                                 dst->width >> ss_x,
                                 dst->stride_cr);
                break;
            default: break;
            }
        }
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

    if (bd > 8) {
        unpack_highbd_pic(dst_buffer_highbd, dst, ss_x, ss_y, EB_TRUE);

        EB_FREE(src_buffer_highbd[0]);
        EB_FREE(src_buffer_highbd[1]);
        EB_FREE(src_buffer_highbd[2]);
        EB_FREE(dst_buffer_highbd[0]);
        EB_FREE(dst_buffer_highbd[1]);
        EB_FREE(dst_buffer_highbd[2]);
    }

    // TODO: extend frame borders
    // use eb_extend_frame() instead
    // aom_extend_frame_borders(dst, num_planes);

    return EB_ErrorNone;
}

// Generate a random number in the range [0, 32768).
static INLINE unsigned int lcg_rand16(unsigned int *state) {
    *state = (unsigned int)(*state * 1103515245ULL + 12345);
    return *state / 65536 % 32768;
}

// Given the superres configurations and the frame type, determine the denominator and
// encoding resolution
void calc_superres_params(superres_params_type *spr_params, SequenceControlSet *scs_ptr,
                          PictureParentControlSet *pcs_ptr) {
    spr_params->superres_denom = SCALE_NUMERATOR;
    static unsigned int seed = 34567;
    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;

    uint8_t superres_mode = scs_ptr->static_config.superres_mode;
    uint8_t cfg_denom     = scs_ptr->static_config.superres_denom;
    uint8_t cfg_kf_denom  = scs_ptr->static_config.superres_kf_denom;
    //uint8_t superres_qthres = scs_ptr->static_config.superres_qthres;

    // For now, super-resolution can only be enabled for key frames or intra only frames
    // In addition, it can only be enabled in case allow_intrabc is disabled and
    // loop restoration is enabled
    if ((frm_hdr->frame_type != KEY_FRAME &&
        frm_hdr->frame_type != INTRA_ONLY_FRAME) ||
        frm_hdr->allow_intrabc ||
        !scs_ptr->seq_header.enable_restoration) { return; }

    // remove assertion when rest of the modes are implemented
    assert(superres_mode <= SUPERRES_RANDOM);

    switch (superres_mode) {
    case SUPERRES_NONE: spr_params->superres_denom = SCALE_NUMERATOR; break;
    case SUPERRES_FIXED:
        if (frm_hdr->frame_type == KEY_FRAME)
            spr_params->superres_denom = cfg_kf_denom;
        else
            spr_params->superres_denom = cfg_denom;
        break;
    case SUPERRES_RANDOM: spr_params->superres_denom = (uint8_t)(lcg_rand16(&seed) % 9 + 8); break;
    //SUPERRES_QTHRESH and SUPERRES_AUTO are not yet implemented
    case SUPERRES_QTHRESH: break;
    case SUPERRES_AUTO: break;
    default: break;
    }

    // only encoding width is adjusted
    calculate_scaled_size_helper(&spr_params->encoding_width, spr_params->superres_denom);
}

EbErrorType downscaled_source_buffer_desc_ctor(EbPictureBufferDesc **picture_ptr,
                                               EbPictureBufferDesc * picture_ptr_for_reference,
                                               superres_params_type  spr_params) {
    EbPictureBufferDescInitData initData;

    initData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    initData.max_width          = spr_params.encoding_width;
    initData.max_height         = spr_params.encoding_height;
    initData.bit_depth          = picture_ptr_for_reference->bit_depth;
    initData.color_format       = picture_ptr_for_reference->color_format;
    initData.split_mode         = EB_TRUE;
    initData.left_padding       = picture_ptr_for_reference->origin_x;
    initData.right_padding      = picture_ptr_for_reference->origin_x;
    initData.top_padding        = picture_ptr_for_reference->origin_y;
    initData.bot_padding        = picture_ptr_for_reference->origin_y;

    EB_NEW(*picture_ptr, eb_picture_buffer_desc_ctor, (EbPtr)&initData);

    return EB_ErrorNone;
}

EbErrorType sb_geom_init_pcs(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr);

EbErrorType sb_params_init_pcs(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr);

EbErrorType scale_pcs_params(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
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

    assert((aligned_width == spr_params.encoding_width) &&
           "Downscaled width needs to be a multiple of 8 "
           "(otherwise not yet implemented)");

    // change frame width and height params in pcs
    pcs_ptr->frame_width  = spr_params.encoding_width;
    pcs_ptr->frame_height = spr_params.encoding_height;

    pcs_ptr->aligned_width  = aligned_width;
    pcs_ptr->aligned_height = aligned_height;

    // number of SBs
    const uint16_t picture_sb_width =
        (uint16_t)((aligned_width + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz);
    const uint16_t picture_sb_height =
        (uint16_t)((aligned_height + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz);

    pcs_ptr->picture_sb_width  = picture_sb_width; // TODO: use this instead of re-computing
    pcs_ptr->picture_sb_height = picture_sb_height;

    pcs_ptr->sb_total_count = picture_sb_width * picture_sb_height;

    // mi params
    cm->mi_stride = picture_sb_width * (BLOCK_SIZE_64 / 4);
    cm->mi_cols   = aligned_width >> MI_SIZE_LOG2;
    cm->mi_rows   = aligned_height >> MI_SIZE_LOG2;

    if (cm->frm_size.superres_denominator != SCALE_NUMERATOR) {
        derive_input_resolution(&pcs_ptr->input_resolution,
                                spr_params.encoding_width * spr_params.encoding_height);

        // create new picture level sb_params and sb_geom
        sb_params_init_pcs(scs_ptr, pcs_ptr);

        sb_geom_init_pcs(scs_ptr, pcs_ptr);
    }

    return EB_ErrorNone;
}

void init_resize_picture(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {
    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_picture_ptr;

    superres_params_type spr_params = {input_picture_ptr->width, // encoding_width
                                       input_picture_ptr->height, // encoding_height
                                       scs_ptr->static_config.superres_mode};

    // determine super-resolution parameters - encoding resolution
    // given configs and frame type
    calc_superres_params(&spr_params, scs_ptr, pcs_ptr);

    if (spr_params.superres_denom != SCALE_NUMERATOR) {

        scs_ptr->seq_header.enable_superres = 1; // enable sequence level super-res flag
                                                 // if super-res is ON for any frame

        // Allocate downsampled picture buffer descriptor
        downscaled_source_buffer_desc_ctor(
            &pcs_ptr->enhanced_downscaled_picture_ptr, input_picture_ptr, spr_params);

        const int32_t  num_planes = av1_num_planes(&scs_ptr->seq_header.color_config);
        const uint32_t ss_x       = scs_ptr->subsampling_x;
        const uint32_t ss_y       = scs_ptr->subsampling_y;

        // downsample picture buffer
        av1_resize_and_extend_frame(input_picture_ptr,
                                    pcs_ptr->enhanced_downscaled_picture_ptr,
                                    pcs_ptr->enhanced_downscaled_picture_ptr->bit_depth,
                                    num_planes,
                                    ss_x,
                                    ss_y);

        // use downscaled picture instead of original res for mode decision, encoding loop etc
        // after temporal filtering and motion estimation
        pcs_ptr->enhanced_picture_ptr = pcs_ptr->enhanced_downscaled_picture_ptr;

        pcs_ptr->frame_superres_enabled = EB_TRUE;

        scale_pcs_params(
            scs_ptr, pcs_ptr, spr_params, input_picture_ptr->width, input_picture_ptr->height);
    }
}
