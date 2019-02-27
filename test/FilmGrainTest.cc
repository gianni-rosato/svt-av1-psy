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

#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "grainSynthesis.h"

static aom_film_grain_t film_grain_test_vectors[3] = {
  /* Test 1 */
  {
      1 /* apply_grain */,
      1 /* update_parameters */,
      { { 16, 0 },
        { 25, 136 },
        { 33, 144 },
        { 41, 160 },
        { 48, 168 },
        { 56, 136 },
        { 67, 128 },
        { 82, 144 },
        { 97, 152 },
        { 113, 144 },
        { 128, 176 },
        { 143, 168 },
        { 158, 176 },
        { 178, 184 } },
      14 /* num_points_y */,
      { { 16, 0 },
        { 20, 64 },
        { 28, 88 },
        { 60, 104 },
        { 90, 136 },
        { 105, 160 },
        { 134, 168 },
        { 168, 208 } },
      8 /* num_cb_points */,
      { { 16, 0 },
        { 28, 96 },
        { 56, 80 },
        { 66, 96 },
        { 80, 104 },
        { 108, 96 },
        { 122, 112 },
        { 137, 112 },
        { 169, 176 } },
      9 /* num_cr_points */,
      11 /* scaling_shift */,
      2 /* ar_coeff_lag */,
      { 0, 0, -58, 0, 0, 0, -76, 100, -43, 0, -51, 82 },
      { 0, 0, -49, 0, 0, 0, -36, 22, -30, 0, -38, 7, 39 },
      { 0, 0, -47, 0, 0, 0, -31, 31, -25, 0, -32, 13, -100 },
      8 /* ar_coeff_shift */,
      247 /* cb_mult */,
      192 /* cb_luma_mult */,
      18 /* cb_offset */,
      229 /* cr_mult */,
      192 /* cr_luma_mult */,
      54 /* cr_offset */,
      0 /* overlap_flag */,
      1 /* clip_to_restricted_range */,
      8 /* bit_depth */,
      0 /* chroma_scaling_from_luma*/,
      0 /* grain_scale_shift*/,
      45231 /* random_seed */
  },
  /* Test 2 */
  {
      1 /* apply_grain */,
      0 /* update_parameters */,
      { { 16, 0 },
        { 25, 136 },
        { 33, 144 },
        { 41, 160 },
        { 48, 168 },
        { 56, 136 },
        { 67, 128 },
        { 82, 144 },
        { 97, 152 },
        { 113, 144 },
        { 128, 176 },
        { 143, 168 },
        { 158, 176 },
        { 178, 184 } },
      14 /* num_points_y */,
      { { 16, 0 },
        { 20, 64 },
        { 28, 88 },
        { 60, 104 },
        { 90, 136 },
        { 105, 160 },
        { 134, 168 },
        { 168, 208 } },
      8 /* num_cb_points */,
      { { 16, 0 },
        { 28, 96 },
        { 56, 80 },
        { 66, 96 },
        { 80, 104 },
        { 108, 96 },
        { 122, 112 },
        { 137, 112 },
        { 169, 176 } },
      9 /* num_cr_points */,
      11 /* scaling_shift */,
      2 /* ar_coeff_lag */,
      { 0, 0, -58, 0, 0, 0, -76, 100, -43, 0, -51, 82 },
      { 0, 0, -49, 0, 0, 0, -36, 22, -30, 0, -38, 7, 39 },
      { 0, 0, -47, 0, 0, 0, -31, 31, -25, 0, -32, 13, -100 },
      8 /* ar_coeff_shift */,
      247 /* cb_mult */,
      192 /* cb_luma_mult */,
      18 /* cb_offset */,
      229 /* cr_mult */,
      192 /* cr_luma_mult */,
      54 /* cr_offset */,
      0 /* overlap_flag */,
      1 /* clip_to_restricted_range */,
      8 /* bit_depth */,
      0 /* chroma_scaling_from_luma*/,
      0 /* grain_scale_shift*/,
      36007 /* random_seed */
  },  
  
  /* Test 3 */
  {
      1 /* apply_grain */,
      1 /* update_parameters */,
      { { 0, 96 }, { 255, 96 } },
      2 /* num_points_y */,
      { { 0, 64 }, { 255, 64 } },
      2 /* num_cb_points */,
      { { 0, 64 }, { 255, 64 } },
      2 /* num_cr_points */,
      11 /* scaling_shift */,
      3 /* ar_coeff_lag */,
      {
          4, 1,   3, 0,   1,  -3, 8,  -3, 7,  -23, 1, -25,
          0, -10, 6, -17, -4, 53, 36, 5,  -5, -17, 8, 66,
      },
      {
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 127,
      },
      {
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 127,
      },
      7 /* ar_coeff_shift */,
      128 /* cb_mult */,
      192 /* cb_luma_mult */,
      256 /* cb_offset */,
      128 /* cr_mult */,
      192 /* cr_luma_mult */,
      256 /* cr_offset */,
      1 /* overlap_flag */,
      0 /* clip_to_restricted_range */,
      8 /* bit_depth */,
      0 /*chroma_scaling_from_luma*/,
      0 /* grain_scale_shift*/,
      45231 /* random_seed */
  },
};

TEST(FilmGrain, parameters_equality)
{
    /* Film grain parameters equality should not depend on random_seed and update_parameters values */
    EXPECT_EQ( film_grain_params_equal( film_grain_test_vectors, film_grain_test_vectors+1), 1);
    
    /* These two instances of film grain parameters are different */
    EXPECT_EQ( film_grain_params_equal( film_grain_test_vectors, film_grain_test_vectors+2), 0);

}


