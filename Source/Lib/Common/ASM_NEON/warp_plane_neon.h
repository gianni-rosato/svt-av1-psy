/*
* Copyright (c) 2023, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#ifndef AOM_AV1_COMMON_ARM_WARP_PLANE_NEON_H_
#define AOM_AV1_COMMON_ARM_WARP_PLANE_NEON_H_

#include <assert.h>
#include <arm_neon.h>
#include <memory.h>
#include <math.h>

#include "EbDefinitions.h"
#include "transpose_neon.h"
#include "common_dsp_rtcd.h"
#include "EbWarpedMotion.h"
#include "convolve.h"
#include "Source/Lib/Encoder/ASM_NEON/sum_neon.h"

// For warping, we really use a 6-tap filter, but we do blocks of 8 pixels
// at a time. The zoom/rotation/shear in the model are applied to the
// "fractional" position of each pixel, which therefore varies within
// [-1, 2) * WARPEDPIXEL_PREC_SHIFTS.
// We need an extra 2 taps to fit this in, for a total of 8 taps.
/* clang-format off */
const int16_t av1_warped_filter[WARPEDPIXEL_PREC_SHIFTS * 3 + 1][8] = {
  // [-1, 0)
  { 0,   0, 127,   1,   0, 0, 0, 0 }, { 0, - 1, 127,   2,   0, 0, 0, 0 },
  { 1, - 3, 127,   4, - 1, 0, 0, 0 }, { 1, - 4, 126,   6, - 2, 1, 0, 0 },
  { 1, - 5, 126,   8, - 3, 1, 0, 0 }, { 1, - 6, 125,  11, - 4, 1, 0, 0 },
  { 1, - 7, 124,  13, - 4, 1, 0, 0 }, { 2, - 8, 123,  15, - 5, 1, 0, 0 },
  { 2, - 9, 122,  18, - 6, 1, 0, 0 }, { 2, -10, 121,  20, - 6, 1, 0, 0 },
  { 2, -11, 120,  22, - 7, 2, 0, 0 }, { 2, -12, 119,  25, - 8, 2, 0, 0 },
  { 3, -13, 117,  27, - 8, 2, 0, 0 }, { 3, -13, 116,  29, - 9, 2, 0, 0 },
  { 3, -14, 114,  32, -10, 3, 0, 0 }, { 3, -15, 113,  35, -10, 2, 0, 0 },
  { 3, -15, 111,  37, -11, 3, 0, 0 }, { 3, -16, 109,  40, -11, 3, 0, 0 },
  { 3, -16, 108,  42, -12, 3, 0, 0 }, { 4, -17, 106,  45, -13, 3, 0, 0 },
  { 4, -17, 104,  47, -13, 3, 0, 0 }, { 4, -17, 102,  50, -14, 3, 0, 0 },
  { 4, -17, 100,  52, -14, 3, 0, 0 }, { 4, -18,  98,  55, -15, 4, 0, 0 },
  { 4, -18,  96,  58, -15, 3, 0, 0 }, { 4, -18,  94,  60, -16, 4, 0, 0 },
  { 4, -18,  91,  63, -16, 4, 0, 0 }, { 4, -18,  89,  65, -16, 4, 0, 0 },
  { 4, -18,  87,  68, -17, 4, 0, 0 }, { 4, -18,  85,  70, -17, 4, 0, 0 },
  { 4, -18,  82,  73, -17, 4, 0, 0 }, { 4, -18,  80,  75, -17, 4, 0, 0 },
  { 4, -18,  78,  78, -18, 4, 0, 0 }, { 4, -17,  75,  80, -18, 4, 0, 0 },
  { 4, -17,  73,  82, -18, 4, 0, 0 }, { 4, -17,  70,  85, -18, 4, 0, 0 },
  { 4, -17,  68,  87, -18, 4, 0, 0 }, { 4, -16,  65,  89, -18, 4, 0, 0 },
  { 4, -16,  63,  91, -18, 4, 0, 0 }, { 4, -16,  60,  94, -18, 4, 0, 0 },
  { 3, -15,  58,  96, -18, 4, 0, 0 }, { 4, -15,  55,  98, -18, 4, 0, 0 },
  { 3, -14,  52, 100, -17, 4, 0, 0 }, { 3, -14,  50, 102, -17, 4, 0, 0 },
  { 3, -13,  47, 104, -17, 4, 0, 0 }, { 3, -13,  45, 106, -17, 4, 0, 0 },
  { 3, -12,  42, 108, -16, 3, 0, 0 }, { 3, -11,  40, 109, -16, 3, 0, 0 },
  { 3, -11,  37, 111, -15, 3, 0, 0 }, { 2, -10,  35, 113, -15, 3, 0, 0 },
  { 3, -10,  32, 114, -14, 3, 0, 0 }, { 2, - 9,  29, 116, -13, 3, 0, 0 },
  { 2, - 8,  27, 117, -13, 3, 0, 0 }, { 2, - 8,  25, 119, -12, 2, 0, 0 },
  { 2, - 7,  22, 120, -11, 2, 0, 0 }, { 1, - 6,  20, 121, -10, 2, 0, 0 },
  { 1, - 6,  18, 122, - 9, 2, 0, 0 }, { 1, - 5,  15, 123, - 8, 2, 0, 0 },
  { 1, - 4,  13, 124, - 7, 1, 0, 0 }, { 1, - 4,  11, 125, - 6, 1, 0, 0 },
  { 1, - 3,   8, 126, - 5, 1, 0, 0 }, { 1, - 2,   6, 126, - 4, 1, 0, 0 },
  { 0, - 1,   4, 127, - 3, 1, 0, 0 }, { 0,   0,   2, 127, - 1, 0, 0, 0 },

  // [0, 1)
  { 0,  0,   0, 127,   1,   0,  0,  0}, { 0,  0,  -1, 127,   2,   0,  0,  0},
  { 0,  1,  -3, 127,   4,  -2,  1,  0}, { 0,  1,  -5, 127,   6,  -2,  1,  0},
  { 0,  2,  -6, 126,   8,  -3,  1,  0}, {-1,  2,  -7, 126,  11,  -4,  2, -1},
  {-1,  3,  -8, 125,  13,  -5,  2, -1}, {-1,  3, -10, 124,  16,  -6,  3, -1},
  {-1,  4, -11, 123,  18,  -7,  3, -1}, {-1,  4, -12, 122,  20,  -7,  3, -1},
  {-1,  4, -13, 121,  23,  -8,  3, -1}, {-2,  5, -14, 120,  25,  -9,  4, -1},
  {-1,  5, -15, 119,  27, -10,  4, -1}, {-1,  5, -16, 118,  30, -11,  4, -1},
  {-2,  6, -17, 116,  33, -12,  5, -1}, {-2,  6, -17, 114,  35, -12,  5, -1},
  {-2,  6, -18, 113,  38, -13,  5, -1}, {-2,  7, -19, 111,  41, -14,  6, -2},
  {-2,  7, -19, 110,  43, -15,  6, -2}, {-2,  7, -20, 108,  46, -15,  6, -2},
  {-2,  7, -20, 106,  49, -16,  6, -2}, {-2,  7, -21, 104,  51, -16,  7, -2},
  {-2,  7, -21, 102,  54, -17,  7, -2}, {-2,  8, -21, 100,  56, -18,  7, -2},
  {-2,  8, -22,  98,  59, -18,  7, -2}, {-2,  8, -22,  96,  62, -19,  7, -2},
  {-2,  8, -22,  94,  64, -19,  7, -2}, {-2,  8, -22,  91,  67, -20,  8, -2},
  {-2,  8, -22,  89,  69, -20,  8, -2}, {-2,  8, -22,  87,  72, -21,  8, -2},
  {-2,  8, -21,  84,  74, -21,  8, -2}, {-2,  8, -22,  82,  77, -21,  8, -2},
  {-2,  8, -21,  79,  79, -21,  8, -2}, {-2,  8, -21,  77,  82, -22,  8, -2},
  {-2,  8, -21,  74,  84, -21,  8, -2}, {-2,  8, -21,  72,  87, -22,  8, -2},
  {-2,  8, -20,  69,  89, -22,  8, -2}, {-2,  8, -20,  67,  91, -22,  8, -2},
  {-2,  7, -19,  64,  94, -22,  8, -2}, {-2,  7, -19,  62,  96, -22,  8, -2},
  {-2,  7, -18,  59,  98, -22,  8, -2}, {-2,  7, -18,  56, 100, -21,  8, -2},
  {-2,  7, -17,  54, 102, -21,  7, -2}, {-2,  7, -16,  51, 104, -21,  7, -2},
  {-2,  6, -16,  49, 106, -20,  7, -2}, {-2,  6, -15,  46, 108, -20,  7, -2},
  {-2,  6, -15,  43, 110, -19,  7, -2}, {-2,  6, -14,  41, 111, -19,  7, -2},
  {-1,  5, -13,  38, 113, -18,  6, -2}, {-1,  5, -12,  35, 114, -17,  6, -2},
  {-1,  5, -12,  33, 116, -17,  6, -2}, {-1,  4, -11,  30, 118, -16,  5, -1},
  {-1,  4, -10,  27, 119, -15,  5, -1}, {-1,  4,  -9,  25, 120, -14,  5, -2},
  {-1,  3,  -8,  23, 121, -13,  4, -1}, {-1,  3,  -7,  20, 122, -12,  4, -1},
  {-1,  3,  -7,  18, 123, -11,  4, -1}, {-1,  3,  -6,  16, 124, -10,  3, -1},
  {-1,  2,  -5,  13, 125,  -8,  3, -1}, {-1,  2,  -4,  11, 126,  -7,  2, -1},
  { 0,  1,  -3,   8, 126,  -6,  2,  0}, { 0,  1,  -2,   6, 127,  -5,  1,  0},
  { 0,  1,  -2,   4, 127,  -3,  1,  0}, { 0,  0,   0,   2, 127,  -1,  0,  0},

  // [1, 2)
  { 0, 0, 0,   1, 127,   0,   0, 0 }, { 0, 0, 0, - 1, 127,   2,   0, 0 },
  { 0, 0, 1, - 3, 127,   4, - 1, 0 }, { 0, 0, 1, - 4, 126,   6, - 2, 1 },
  { 0, 0, 1, - 5, 126,   8, - 3, 1 }, { 0, 0, 1, - 6, 125,  11, - 4, 1 },
  { 0, 0, 1, - 7, 124,  13, - 4, 1 }, { 0, 0, 2, - 8, 123,  15, - 5, 1 },
  { 0, 0, 2, - 9, 122,  18, - 6, 1 }, { 0, 0, 2, -10, 121,  20, - 6, 1 },
  { 0, 0, 2, -11, 120,  22, - 7, 2 }, { 0, 0, 2, -12, 119,  25, - 8, 2 },
  { 0, 0, 3, -13, 117,  27, - 8, 2 }, { 0, 0, 3, -13, 116,  29, - 9, 2 },
  { 0, 0, 3, -14, 114,  32, -10, 3 }, { 0, 0, 3, -15, 113,  35, -10, 2 },
  { 0, 0, 3, -15, 111,  37, -11, 3 }, { 0, 0, 3, -16, 109,  40, -11, 3 },
  { 0, 0, 3, -16, 108,  42, -12, 3 }, { 0, 0, 4, -17, 106,  45, -13, 3 },
  { 0, 0, 4, -17, 104,  47, -13, 3 }, { 0, 0, 4, -17, 102,  50, -14, 3 },
  { 0, 0, 4, -17, 100,  52, -14, 3 }, { 0, 0, 4, -18,  98,  55, -15, 4 },
  { 0, 0, 4, -18,  96,  58, -15, 3 }, { 0, 0, 4, -18,  94,  60, -16, 4 },
  { 0, 0, 4, -18,  91,  63, -16, 4 }, { 0, 0, 4, -18,  89,  65, -16, 4 },
  { 0, 0, 4, -18,  87,  68, -17, 4 }, { 0, 0, 4, -18,  85,  70, -17, 4 },
  { 0, 0, 4, -18,  82,  73, -17, 4 }, { 0, 0, 4, -18,  80,  75, -17, 4 },
  { 0, 0, 4, -18,  78,  78, -18, 4 }, { 0, 0, 4, -17,  75,  80, -18, 4 },
  { 0, 0, 4, -17,  73,  82, -18, 4 }, { 0, 0, 4, -17,  70,  85, -18, 4 },
  { 0, 0, 4, -17,  68,  87, -18, 4 }, { 0, 0, 4, -16,  65,  89, -18, 4 },
  { 0, 0, 4, -16,  63,  91, -18, 4 }, { 0, 0, 4, -16,  60,  94, -18, 4 },
  { 0, 0, 3, -15,  58,  96, -18, 4 }, { 0, 0, 4, -15,  55,  98, -18, 4 },
  { 0, 0, 3, -14,  52, 100, -17, 4 }, { 0, 0, 3, -14,  50, 102, -17, 4 },
  { 0, 0, 3, -13,  47, 104, -17, 4 }, { 0, 0, 3, -13,  45, 106, -17, 4 },
  { 0, 0, 3, -12,  42, 108, -16, 3 }, { 0, 0, 3, -11,  40, 109, -16, 3 },
  { 0, 0, 3, -11,  37, 111, -15, 3 }, { 0, 0, 2, -10,  35, 113, -15, 3 },
  { 0, 0, 3, -10,  32, 114, -14, 3 }, { 0, 0, 2, - 9,  29, 116, -13, 3 },
  { 0, 0, 2, - 8,  27, 117, -13, 3 }, { 0, 0, 2, - 8,  25, 119, -12, 2 },
  { 0, 0, 2, - 7,  22, 120, -11, 2 }, { 0, 0, 1, - 6,  20, 121, -10, 2 },
  { 0, 0, 1, - 6,  18, 122, - 9, 2 }, { 0, 0, 1, - 5,  15, 123, - 8, 2 },
  { 0, 0, 1, - 4,  13, 124, - 7, 1 }, { 0, 0, 1, - 4,  11, 125, - 6, 1 },
  { 0, 0, 1, - 3,   8, 126, - 5, 1 }, { 0, 0, 1, - 2,   6, 126, - 4, 1 },
  { 0, 0, 0, - 1,   4, 127, - 3, 1 }, { 0, 0, 0,   0,   2, 127, - 1, 0 },
  // dummy (replicate row index 191)
  { 0, 0, 0,   0,   2, 127, - 1, 0 },
};

static INLINE int16x8_t horizontal_filter_4x1_f4(const uint8x16_t in, int sx,
                                                 int alpha);

static INLINE int16x8_t horizontal_filter_8x1_f8(const uint8x16_t in, int sx,
                                                 int alpha);

static INLINE int16x8_t horizontal_filter_4x1_f1(const uint8x16_t in, int sx);

static INLINE int16x8_t horizontal_filter_8x1_f1(const uint8x16_t in, int sx);

static INLINE void vertical_filter_4x1_f1(const int16x8_t *src, int32x4_t *res,
                                          int sy);

static INLINE void vertical_filter_4x1_f4(const int16x8_t *src, int32x4_t *res,
                                          int sy, int gamma);

static INLINE void vertical_filter_8x1_f1(const int16x8_t *src,
                                          int32x4_t *res_low,
                                          int32x4_t *res_high, int sy);

static INLINE void vertical_filter_8x1_f8(const int16x8_t *src,
                                          int32x4_t *res_low,
                                          int32x4_t *res_high, int sy,
                                          int gamma);

static INLINE void load_filters_4(int16x8_t out[], int offset, int stride) {
  out[0] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 0 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[1] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 1 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[2] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 2 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[3] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 3 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
}

static INLINE void load_filters_8(int16x8_t out[], int offset, int stride) {
  out[0] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 0 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[1] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 1 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[2] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 2 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[3] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 3 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[4] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 4 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[5] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 5 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[6] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 6 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
  out[7] = vld1q_s16((int16_t *)(av1_warped_filter + ((offset + 7 * stride) >>
                                                      WARPEDDIFF_PREC_BITS)));
}

static INLINE int clamp_iy(int iy, int height) {
  return clamp(iy, 0, height - 1);
}

static INLINE void warp_affine_horizontal(
    const uint8_t *ref, int width, int height, int stride, int p_width,
    int p_height, int16_t alpha, int16_t beta, const int64_t x4,
    const int64_t y4, const int i, int16x8_t tmp[], const uint8x16_t indx_vec) {
  const int bd = 8;
  const int reduce_bits_horiz = ROUND0_BITS;
  const int height_limit = AOMMIN(8, p_height - i) + 7;

  int32_t ix4 = (int32_t)(x4 >> WARPEDMODEL_PREC_BITS);
  int32_t iy4 = (int32_t)(y4 >> WARPEDMODEL_PREC_BITS);

  int32_t sx4 = x4 & ((1 << WARPEDMODEL_PREC_BITS) - 1);
  sx4 += alpha * (-4) + beta * (-4) + (1 << (WARPEDDIFF_PREC_BITS - 1)) +
         (WARPEDPIXEL_PREC_SHIFTS << WARPEDDIFF_PREC_BITS);
  sx4 &= ~((1 << WARP_PARAM_REDUCE_BITS) - 1);

  if (ix4 <= -7) {
    for (int k = 0; k < height_limit; ++k) {
      int iy = clamp_iy(iy4 + k - 7, height);
      int16_t dup_val =
          (1 << (bd + FILTER_BITS - reduce_bits_horiz - 1)) +
          ref[iy * stride] * (1 << (FILTER_BITS - reduce_bits_horiz));
      tmp[k] = vdupq_n_s16(dup_val);
    }
    return;
  } else if (ix4 >= width + 6) {
    for (int k = 0; k < height_limit; ++k) {
      int iy = clamp_iy(iy4 + k - 7, height);
      int16_t dup_val = (1 << (bd + FILTER_BITS - reduce_bits_horiz - 1)) +
                        ref[iy * stride + (width - 1)] *
                            (1 << (FILTER_BITS - reduce_bits_horiz));
      tmp[k] = vdupq_n_s16(dup_val);
    }
    return;
  }

  uint8x16_t in[15];
  if (((ix4 - 7) < 0) || ((ix4 + 9) > width)) {
    const int out_of_boundary_left = -(ix4 - 6);
    const int out_of_boundary_right = (ix4 + 8) - width;

    for (int k = 0; k < height_limit; ++k) {
      const int iy = clamp_iy(iy4 + k - 7, height);
      const uint8_t *src = ref + iy * stride + ix4 - 7;
      uint8x16_t src_1 = vld1q_u8(src);

      if (out_of_boundary_left >= 0) {
        int limit = out_of_boundary_left + 1;
        uint8x16_t cmp_vec = vdupq_n_u8(out_of_boundary_left);
        uint8x16_t vec_dup = vdupq_n_u8(*(src + limit));
        uint8x16_t mask_val = vcleq_u8(indx_vec, cmp_vec);
        src_1 = vbslq_u8(mask_val, vec_dup, src_1);
      }
      if (out_of_boundary_right >= 0) {
        int limit = 15 - (out_of_boundary_right + 1);
        uint8x16_t cmp_vec = vdupq_n_u8(15 - out_of_boundary_right);
        uint8x16_t vec_dup = vdupq_n_u8(*(src + limit));
        uint8x16_t mask_val = vcgeq_u8(indx_vec, cmp_vec);
        src_1 = vbslq_u8(mask_val, vec_dup, src_1);
      }
      in[k] = src_1;
    }
  } else {
    for (int k = 0; k < height_limit; ++k) {
      const int iy = clamp_iy(iy4 + k - 7, height);
      const uint8_t *src = ref + iy * stride + ix4 - 7;
      in[k] = vld1q_u8(src);
    }
  }

  if (p_width == 4) {
    if (beta == 0) {
      if (alpha == 0) {
        for (int k = 0; k < height_limit; ++k) {
          tmp[k] = horizontal_filter_4x1_f1(in[k], sx4);
        }
      } else {
        for (int k = 0; k < height_limit; ++k) {
          tmp[k] = horizontal_filter_4x1_f4(in[k], sx4, alpha);
        }
      }
    } else {
      if (alpha == 0) {
        for (int k = 0; k < height_limit; ++k) {
          const int sx = sx4 + beta * (k - 3);
          tmp[k] = horizontal_filter_4x1_f1(in[k], sx);
        }
      } else {
        for (int k = 0; k < height_limit; ++k) {
          const int sx = sx4 + beta * (k - 3);
          tmp[k] = horizontal_filter_4x1_f4(in[k], sx, alpha);
        }
      }
    }
  } else {
    if (beta == 0) {
      if (alpha == 0) {
        for (int k = 0; k < height_limit; ++k) {
          tmp[k] = horizontal_filter_8x1_f1(in[k], sx4);
        }
      } else {
        for (int k = 0; k < height_limit; ++k) {
          tmp[k] = horizontal_filter_8x1_f8(in[k], sx4, alpha);
        }
      }
    } else {
      if (alpha == 0) {
        for (int k = 0; k < height_limit; ++k) {
          const int sx = sx4 + beta * (k - 3);
          tmp[k] = horizontal_filter_8x1_f1(in[k], sx);
        }
      } else {
        for (int k = 0; k < height_limit; ++k) {
          const int sx = sx4 + beta * (k - 3);
          tmp[k] = horizontal_filter_8x1_f8(in[k], sx, alpha);
        }
      }
    }
  }
}

static INLINE void warp_affine_vertical(
    uint8_t *pred, int p_width, int p_height, int p_stride, int is_compound,
    uint16_t *dst, int dst_stride, int do_average, int use_dist_wtd_comp_avg,
    int16_t gamma, int16_t delta, const int64_t y4, const int i, const int j,
    int16x8_t tmp[], const int fwd, const int bwd) {
  const int bd = 8;
  const int reduce_bits_horiz = ROUND0_BITS;
  const int offset_bits_vert = bd + 2 * FILTER_BITS - reduce_bits_horiz;
  int add_const_vert;
  if (is_compound) {
    add_const_vert =
        (1 << offset_bits_vert) + (1 << (COMPOUND_ROUND1_BITS - 1));
  } else {
    add_const_vert =
        (1 << offset_bits_vert) + (1 << (2 * FILTER_BITS - ROUND0_BITS - 1));
  }
  const int sub_constant = (1 << (bd - 1)) + (1 << bd);

  const int offset_bits = bd + 2 * FILTER_BITS - ROUND0_BITS;
  const int res_sub_const =
      (1 << (2 * FILTER_BITS - ROUND0_BITS - COMPOUND_ROUND1_BITS - 1)) -
      (1 << (offset_bits - COMPOUND_ROUND1_BITS)) -
      (1 << (offset_bits - COMPOUND_ROUND1_BITS - 1));

  int32_t sy4 = y4 & ((1 << WARPEDMODEL_PREC_BITS) - 1);
  sy4 += gamma * (-4) + delta * (-4) + (1 << (WARPEDDIFF_PREC_BITS - 1)) +
         (WARPEDPIXEL_PREC_SHIFTS << WARPEDDIFF_PREC_BITS);
  sy4 &= ~((1 << WARP_PARAM_REDUCE_BITS) - 1);

  if (p_width > 4) {
    for (int k = -4; k < AOMMIN(4, p_height - i - 4); ++k) {
      int sy = sy4 + delta * (k + 4);
      const int16x8_t *v_src = tmp + (k + 4);

      int32x4_t res_lo, res_hi;
      if (gamma == 0) {
        vertical_filter_8x1_f1(v_src, &res_lo, &res_hi, sy);
      } else {
        vertical_filter_8x1_f8(v_src, &res_lo, &res_hi, sy, gamma);
      }

      res_lo = vaddq_s32(res_lo, vdupq_n_s32(add_const_vert));
      res_hi = vaddq_s32(res_hi, vdupq_n_s32(add_const_vert));

      if (is_compound) {
        uint16_t *const p = (uint16_t *)&dst[(i + k + 4) * dst_stride + j];
        int16x8_t res_s16 =
            vcombine_s16(vshrn_n_s32(res_lo, COMPOUND_ROUND1_BITS),
                         vshrn_n_s32(res_hi, COMPOUND_ROUND1_BITS));
        if (do_average) {
          int16x8_t tmp16 = vreinterpretq_s16_u16(vld1q_u16(p));
          if (use_dist_wtd_comp_avg) {
            int32x4_t tmp32_lo = vmull_n_s16(vget_low_s16(tmp16), fwd);
            int32x4_t tmp32_hi = vmull_n_s16(vget_high_s16(tmp16), fwd);
            tmp32_lo = vmlal_n_s16(tmp32_lo, vget_low_s16(res_s16), bwd);
            tmp32_hi = vmlal_n_s16(tmp32_hi, vget_high_s16(res_s16), bwd);
            tmp16 = vcombine_s16(vshrn_n_s32(tmp32_lo, DIST_PRECISION_BITS),
                                 vshrn_n_s32(tmp32_hi, DIST_PRECISION_BITS));
          } else {
            tmp16 = vhaddq_s16(tmp16, res_s16);
          }
          int16x8_t res = vaddq_s16(tmp16, vdupq_n_s16(res_sub_const));
          uint8x8_t res8 = vqshrun_n_s16(
              res, 2 * FILTER_BITS - ROUND0_BITS - COMPOUND_ROUND1_BITS);
          vst1_u8(&pred[(i + k + 4) * p_stride + j], res8);
        } else {
          vst1q_u16(p, vreinterpretq_u16_s16(res_s16));
        }
      } else {
        int16x8_t res16 =
            vcombine_s16(vshrn_n_s32(res_lo, 2 * FILTER_BITS - ROUND0_BITS),
                         vshrn_n_s32(res_hi, 2 * FILTER_BITS - ROUND0_BITS));
        res16 = vsubq_s16(res16, vdupq_n_s16(sub_constant));

        uint8_t *const p = (uint8_t *)&pred[(i + k + 4) * p_stride + j];
        vst1_u8(p, vqmovun_s16(res16));
      }
    }
  } else {
    // p_width == 4
    for (int k = -4; k < AOMMIN(4, p_height - i - 4); ++k) {
      int sy = sy4 + delta * (k + 4);
      const int16x8_t *v_src = tmp + (k + 4);

      int32x4_t res_lo;
      if (gamma == 0) {
        vertical_filter_4x1_f1(v_src, &res_lo, sy);
      } else {
        vertical_filter_4x1_f4(v_src, &res_lo, sy, gamma);
      }

      res_lo = vaddq_s32(res_lo, vdupq_n_s32(add_const_vert));

      if (is_compound) {
        uint16_t *const p = (uint16_t *)&dst[(i + k + 4) * dst_stride + j];

        int16x4_t res_lo_s16 = vshrn_n_s32(res_lo, COMPOUND_ROUND1_BITS);
        if (do_average) {
          uint8_t *const dst8 = &pred[(i + k + 4) * p_stride + j];
          int16x4_t tmp16_lo = vreinterpret_s16_u16(vld1_u16(p));
          if (use_dist_wtd_comp_avg) {
            int32x4_t tmp32_lo = vmull_n_s16(tmp16_lo, fwd);
            tmp32_lo = vmlal_n_s16(tmp32_lo, res_lo_s16, bwd);
            tmp16_lo = vshrn_n_s32(tmp32_lo, DIST_PRECISION_BITS);
          } else {
            tmp16_lo = vhadd_s16(tmp16_lo, res_lo_s16);
          }
          int16x4_t res = vadd_s16(tmp16_lo, vdup_n_s16(res_sub_const));
          uint8x8_t res8 = vqshrun_n_s16(
              vcombine_s16(res, vdup_n_s16(0)),
              2 * FILTER_BITS - ROUND0_BITS - COMPOUND_ROUND1_BITS);
          vst1_lane_u32((uint32_t *)dst8, vreinterpret_u32_u8(res8), 0);
        } else {
          uint16x4_t res_u16_low = vreinterpret_u16_s16(res_lo_s16);
          vst1_u16(p, res_u16_low);
        }
      } else {
        int16x4_t res16 = vshrn_n_s32(res_lo, 2 * FILTER_BITS - ROUND0_BITS);
        res16 = vsub_s16(res16, vdup_n_s16(sub_constant));

        uint8_t *const p = (uint8_t *)&pred[(i + k + 4) * p_stride + j];
        uint8x8_t val = vqmovun_s16(vcombine_s16(res16, vdup_n_s16(0)));
        vst1_lane_u32((uint32_t *)p, vreinterpret_u32_u8(val), 0);
      }
    }
  }
}

#endif  // AOM_AV1_COMMON_ARM_WARP_PLANE_NEON_H_
