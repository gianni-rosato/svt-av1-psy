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

#include <stdlib.h>
#include "aom_dsp_rtcd.h"
#include "EbEncWarpedMotion.h"
//#include "EbPictureControlSet.h"
#include "EbBitstreamUnit.h"
#include "EbDefinitions.h"
#include "convolve.h"

#define WARP_ERROR_BLOCK 32

/* clang-format off */
const int svt_aom_error_measure_lut[512] = {
  // pow 0.7
  16384, 16339, 16294, 16249, 16204, 16158, 16113, 16068,
  16022, 15977, 15932, 15886, 15840, 15795, 15749, 15703,
  15657, 15612, 15566, 15520, 15474, 15427, 15381, 15335,
  15289, 15242, 15196, 15149, 15103, 15056, 15010, 14963,
  14916, 14869, 14822, 14775, 14728, 14681, 14634, 14587,
  14539, 14492, 14445, 14397, 14350, 14302, 14254, 14206,
  14159, 14111, 14063, 14015, 13967, 13918, 13870, 13822,
  13773, 13725, 13676, 13628, 13579, 13530, 13481, 13432,
  13383, 13334, 13285, 13236, 13187, 13137, 13088, 13038,
  12988, 12939, 12889, 12839, 12789, 12739, 12689, 12639,
  12588, 12538, 12487, 12437, 12386, 12335, 12285, 12234,
  12183, 12132, 12080, 12029, 11978, 11926, 11875, 11823,
  11771, 11719, 11667, 11615, 11563, 11511, 11458, 11406,
  11353, 11301, 11248, 11195, 11142, 11089, 11036, 10982,
  10929, 10875, 10822, 10768, 10714, 10660, 10606, 10552,
  10497, 10443, 10388, 10333, 10279, 10224, 10168, 10113,
  10058, 10002,  9947,  9891,  9835,  9779,  9723,  9666,
  9610, 9553, 9497, 9440, 9383, 9326, 9268, 9211,
  9153, 9095, 9037, 8979, 8921, 8862, 8804, 8745,
  8686, 8627, 8568, 8508, 8449, 8389, 8329, 8269,
  8208, 8148, 8087, 8026, 7965, 7903, 7842, 7780,
  7718, 7656, 7593, 7531, 7468, 7405, 7341, 7278,
  7214, 7150, 7086, 7021, 6956, 6891, 6826, 6760,
  6695, 6628, 6562, 6495, 6428, 6361, 6293, 6225,
  6157, 6089, 6020, 5950, 5881, 5811, 5741, 5670,
  5599, 5527, 5456, 5383, 5311, 5237, 5164, 5090,
  5015, 4941, 4865, 4789, 4713, 4636, 4558, 4480,
  4401, 4322, 4242, 4162, 4080, 3998, 3916, 3832,
  3748, 3663, 3577, 3490, 3402, 3314, 3224, 3133,
  3041, 2948, 2854, 2758, 2661, 2562, 2461, 2359,
  2255, 2148, 2040, 1929, 1815, 1698, 1577, 1452,
  1323, 1187, 1045,  894,  731,  550,  339,    0,
  339,  550,  731,  894, 1045, 1187, 1323, 1452,
  1577, 1698, 1815, 1929, 2040, 2148, 2255, 2359,
  2461, 2562, 2661, 2758, 2854, 2948, 3041, 3133,
  3224, 3314, 3402, 3490, 3577, 3663, 3748, 3832,
  3916, 3998, 4080, 4162, 4242, 4322, 4401, 4480,
  4558, 4636, 4713, 4789, 4865, 4941, 5015, 5090,
  5164, 5237, 5311, 5383, 5456, 5527, 5599, 5670,
  5741, 5811, 5881, 5950, 6020, 6089, 6157, 6225,
  6293, 6361, 6428, 6495, 6562, 6628, 6695, 6760,
  6826, 6891, 6956, 7021, 7086, 7150, 7214, 7278,
  7341, 7405, 7468, 7531, 7593, 7656, 7718, 7780,
  7842, 7903, 7965, 8026, 8087, 8148, 8208, 8269,
  8329, 8389, 8449, 8508, 8568, 8627, 8686, 8745,
  8804, 8862, 8921, 8979, 9037, 9095, 9153, 9211,
  9268, 9326, 9383, 9440, 9497, 9553, 9610, 9666,
  9723,  9779,  9835,  9891,  9947, 10002, 10058, 10113,
  10168, 10224, 10279, 10333, 10388, 10443, 10497, 10552,
  10606, 10660, 10714, 10768, 10822, 10875, 10929, 10982,
  11036, 11089, 11142, 11195, 11248, 11301, 11353, 11406,
  11458, 11511, 11563, 11615, 11667, 11719, 11771, 11823,
  11875, 11926, 11978, 12029, 12080, 12132, 12183, 12234,
  12285, 12335, 12386, 12437, 12487, 12538, 12588, 12639,
  12689, 12739, 12789, 12839, 12889, 12939, 12988, 13038,
  13088, 13137, 13187, 13236, 13285, 13334, 13383, 13432,
  13481, 13530, 13579, 13628, 13676, 13725, 13773, 13822,
  13870, 13918, 13967, 14015, 14063, 14111, 14159, 14206,
  14254, 14302, 14350, 14397, 14445, 14492, 14539, 14587,
  14634, 14681, 14728, 14775, 14822, 14869, 14916, 14963,
  15010, 15056, 15103, 15149, 15196, 15242, 15289, 15335,
  15381, 15427, 15474, 15520, 15566, 15612, 15657, 15703,
  15749, 15795, 15840, 15886, 15932, 15977, 16022, 16068,
  16113, 16158, 16204, 16249, 16294, 16339, 16384, 16384,
};
/* clang-format on */

int64_t svt_av1_calc_frame_error_c(const uint8_t *const ref, int stride, const uint8_t *const dst, int p_width,
                                   int p_height, int p_stride) {
    int64_t sum_error = 0;
    for (int i = 0; i < p_height; ++i) {
        for (int j = 0; j < p_width; ++j) {
            sum_error += (int64_t)error_measure(dst[j + i * p_stride] - ref[j + i * stride]);
        }
    }
    return sum_error;
}

static int64_t warp_error(EbWarpedMotionParams *wm, const uint8_t *const ref, int width, int height, int stride,
                          const uint8_t *const dst, int p_col, int p_row, int p_width, int p_height, int p_stride,
                          int subsampling_x, int subsampling_y, uint8_t chess_refn, int64_t best_error) {
    int64_t        gm_sumerr = 0;
    int            warp_w, warp_h;
    int            error_bsize_w = AOMMIN(p_width, WARP_ERROR_BLOCK);
    int            error_bsize_h = AOMMIN(p_height, WARP_ERROR_BLOCK);
    uint8_t        tmp[WARP_ERROR_BLOCK * WARP_ERROR_BLOCK];
    ConvolveParams conv_params   = get_conv_params(0, 0, 0, 8);
    conv_params.use_jnt_comp_avg = 0;

    int i_itr = 0;
    for (int i = p_row; i < p_row + p_height; i += WARP_ERROR_BLOCK) {
        int jstart = (i_itr & 1) ? p_col : p_col + WARP_ERROR_BLOCK;
        int jstep  = 2;

        if (chess_refn == 0) {
            jstart = p_col;
            jstep  = 1;
        }

        for (int j = jstart; j < p_col + p_width; j += jstep * WARP_ERROR_BLOCK) {
            // avoid warping extra 8x8 blocks in the padded region of the frame
            // when p_width and p_height are not multiples of WARP_ERROR_BLOCK
            warp_w = AOMMIN(error_bsize_w, p_col + p_width - j);
            warp_h = AOMMIN(error_bsize_h, p_row + p_height - i);
            svt_warp_plane(wm,
                           ref,
                           width,
                           height,
                           stride,
                           tmp,
                           j,
                           i,
                           warp_w,
                           warp_h,
                           WARP_ERROR_BLOCK,
                           subsampling_x,
                           subsampling_y,
                           &conv_params);

            gm_sumerr += svt_av1_calc_frame_error(
                tmp, WARP_ERROR_BLOCK, dst + j + i * p_stride, warp_w, warp_h, p_stride);
            if (gm_sumerr > best_error)
                return gm_sumerr;
        }

        i_itr++;
    }

    gm_sumerr = chess_refn ? gm_sumerr * 2 : gm_sumerr;

    return gm_sumerr;
}

int64_t svt_av1_warp_error(EbWarpedMotionParams *wm, const uint8_t *ref, int width, int height, int stride,
                           uint8_t *dst, int p_col, int p_row, int p_width, int p_height, int p_stride,
                           int subsampling_x, int subsampling_y, uint8_t chess_refn, int64_t best_error) {
    if (wm->wmtype <= AFFINE)
        if (!svt_get_shear_params(wm))
            return 1;
    return warp_error(wm,
                      ref,
                      width,
                      height,
                      stride,
                      dst,
                      p_col,
                      p_row,
                      p_width,
                      p_height,
                      p_stride,
                      subsampling_x,
                      subsampling_y,
                      chess_refn,
                      best_error);
}
