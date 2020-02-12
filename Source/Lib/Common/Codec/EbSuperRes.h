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

#ifndef EbSuperRes_h
#define EbSuperRes_h


#ifdef __cplusplus
extern "C" {
#endif


#define UPSCALE_NORMATIVE_TAPS 8

static const int16_t av1_resize_filter_normative[(1 << RS_SUBPEL_BITS)][UPSCALE_NORMATIVE_TAPS] = {
#if UPSCALE_NORMATIVE_TAPS == 8
{0, 0, 0, 128, 0, 0, 0, 0},        {0, 0, -1, 128, 2, -1, 0, 0},
    {0, 1, -3, 127, 4, -2, 1, 0},      {0, 1, -4, 127, 6, -3, 1, 0},
    {0, 2, -6, 126, 8, -3, 1, 0},      {0, 2, -7, 125, 11, -4, 1, 0},
    {-1, 2, -8, 125, 13, -5, 2, 0},    {-1, 3, -9, 124, 15, -6, 2, 0},
    {-1, 3, -10, 123, 18, -6, 2, -1},  {-1, 3, -11, 122, 20, -7, 3, -1},
    {-1, 4, -12, 121, 22, -8, 3, -1},  {-1, 4, -13, 120, 25, -9, 3, -1},
    {-1, 4, -14, 118, 28, -9, 3, -1},  {-1, 4, -15, 117, 30, -10, 4, -1},
    {-1, 5, -16, 116, 32, -11, 4, -1}, {-1, 5, -16, 114, 35, -12, 4, -1},
    {-1, 5, -17, 112, 38, -12, 4, -1}, {-1, 5, -18, 111, 40, -13, 5, -1},
    {-1, 5, -18, 109, 43, -14, 5, -1}, {-1, 6, -19, 107, 45, -14, 5, -1},
    {-1, 6, -19, 105, 48, -15, 5, -1}, {-1, 6, -19, 103, 51, -16, 5, -1},
    {-1, 6, -20, 101, 53, -16, 6, -1}, {-1, 6, -20, 99, 56, -17, 6, -1},
    {-1, 6, -20, 97, 58, -17, 6, -1},  {-1, 6, -20, 95, 61, -18, 6, -1},
    {-2, 7, -20, 93, 64, -18, 6, -2},  {-2, 7, -20, 91, 66, -19, 6, -1},
    {-2, 7, -20, 88, 69, -19, 6, -1},  {-2, 7, -20, 86, 71, -19, 6, -1},
    {-2, 7, -20, 84, 74, -20, 7, -2},  {-2, 7, -20, 81, 76, -20, 7, -1},
    {-2, 7, -20, 79, 79, -20, 7, -2},  {-1, 7, -20, 76, 81, -20, 7, -2},
    {-2, 7, -20, 74, 84, -20, 7, -2},  {-1, 6, -19, 71, 86, -20, 7, -2},
    {-1, 6, -19, 69, 88, -20, 7, -2},  {-1, 6, -19, 66, 91, -20, 7, -2},
    {-2, 6, -18, 64, 93, -20, 7, -2},  {-1, 6, -18, 61, 95, -20, 6, -1},
    {-1, 6, -17, 58, 97, -20, 6, -1},  {-1, 6, -17, 56, 99, -20, 6, -1},
    {-1, 6, -16, 53, 101, -20, 6, -1}, {-1, 5, -16, 51, 103, -19, 6, -1},
    {-1, 5, -15, 48, 105, -19, 6, -1}, {-1, 5, -14, 45, 107, -19, 6, -1},
    {-1, 5, -14, 43, 109, -18, 5, -1}, {-1, 5, -13, 40, 111, -18, 5, -1},
    {-1, 4, -12, 38, 112, -17, 5, -1}, {-1, 4, -12, 35, 114, -16, 5, -1},
    {-1, 4, -11, 32, 116, -16, 5, -1}, {-1, 4, -10, 30, 117, -15, 4, -1},
    {-1, 3, -9, 28, 118, -14, 4, -1},  {-1, 3, -9, 25, 120, -13, 4, -1},
    {-1, 3, -8, 22, 121, -12, 4, -1},  {-1, 3, -7, 20, 122, -11, 3, -1},
    {-1, 2, -6, 18, 123, -10, 3, -1},  {0, 2, -6, 15, 124, -9, 3, -1},
    {0, 2, -5, 13, 125, -8, 2, -1},    {0, 1, -4, 11, 125, -7, 2, 0},
    {0, 1, -3, 8, 126, -6, 2, 0},      {0, 1, -3, 6, 127, -4, 1, 0},
    {0, 1, -2, 4, 127, -3, 1, 0},      {0, 0, -1, 2, 128, -1, 0, 0},
#else
#error "Invalid value of UPSCALE_NORMATIVE_TAPS"
#endif // UPSCALE_NORMATIVE_TAPS == 8
};
// Filters for interpolation (full-band) - no filtering for integer pixels

#ifdef __cplusplus
}
#endif
#endif // EbSuperRes_h
