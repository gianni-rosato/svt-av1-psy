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

/*!\file
  * \brief Describes film grain parameters and film grain synthesis
  *
  */
#ifndef AOM_AOM_GRAIN_SYNTHESIS_H_
#define AOM_AOM_GRAIN_SYNTHESIS_H_

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

/*!\brief Structure containing film grain synthesis parameters for a frame
     *
     * This structure contains input parameters for film grain synthesis
     */
typedef struct {
    int32_t apply_grain;

    int32_t update_parameters;

    // 8 bit values
    int32_t scaling_points_y[14][2];
    int32_t num_y_points; // value: 0..14

    // 8 bit values
    int32_t scaling_points_cb[10][2];
    int32_t num_cb_points; // value: 0..10

    // 8 bit values
    int32_t scaling_points_cr[10][2];
    int32_t num_cr_points; // value: 0..10

    int32_t scaling_shift; // values : 8..11

    int32_t ar_coeff_lag; // values:  0..3

    // 8 bit values
    int32_t ar_coeffs_y[24];
    int32_t ar_coeffs_cb[25];
    int32_t ar_coeffs_cr[25];

    // Shift value: AR coeffs range
    // 6: [-2, 2)
    // 7: [-1, 1)
    // 8: [-0.5, 0.5)
    // 9: [-0.25, 0.25)
    int32_t ar_coeff_shift; // values : 6..9

    int32_t cb_mult; // 8 bits
    int32_t cb_luma_mult; // 8 bits
    int32_t cb_offset; // 9 bits

    int32_t cr_mult; // 8 bits
    int32_t cr_luma_mult; // 8 bits
    int32_t cr_offset; // 9 bits

    int32_t overlap_flag;

    int32_t clip_to_restricted_range;

    int32_t bit_depth; // video bit depth

    int32_t chroma_scaling_from_luma;

    int32_t grain_scale_shift;

    uint16_t random_seed;
} AomFilmGrain;

int32_t film_grain_params_equal(AomFilmGrain *pars_a, AomFilmGrain *pars_b);

/*!\brief Add film grain
     *
     * Add film grain to an image
     *
     * \param[in]    grain_params     Grain parameters
     * \param[in]    luma             luma plane
     * \param[in]    cb               cb plane
     * \param[in]    cr               cr plane
     * \param[in]    height           luma plane height
     * \param[in]    width            luma plane width
     * \param[in]    luma_stride      luma plane stride
     * \param[in]    chroma_stride    chroma plane stride
     */
void eb_av1_add_film_grain_run(AomFilmGrain *grain_params, uint8_t *luma, uint8_t *cb, uint8_t *cr,
                               int32_t height, int32_t width, int32_t luma_stride,
                               int32_t chroma_stride, int32_t use_high_bit_depth,
                               int32_t chroma_subsamp_y, int32_t chroma_subsamp_x);

/*!\brief Add film grain
     *
     * Add film grain to an image
     *
     * \param[in]    grain_params     Grain parameters
     * \param[in]    src              Source image
     * \param[in]    dst              Resulting image with grain
     */

//void eb_av1_add_film_grain(AomFilmGrain *grain_params, EbPictureBufferDesc *src,
//        EbPictureBufferDesc *dst);

//void av1_film_grain_write_updated(const AomFilmGrain *pars,
//  int32_t monochrome, struct AomWriteBitBuffer *wb);
//
//void av1_film_grain_read_updated(AomFilmGrain *pars,
//                                 int32_t monochrome, struct aom_read_bit_buffer *wb,
//                                 struct aom_internal_error_info *error);

void fgn_copy_rect(uint8_t *src, int32_t src_stride, uint8_t *dst, int32_t dst_stride,
                   int32_t width, int32_t height, int32_t use_high_bit_depth);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // AOM_AOM_GRAIN_SYNTHESIS_H_
