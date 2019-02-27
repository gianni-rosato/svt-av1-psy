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

#include <stdlib.h>

#include "EbPictureControlSet.h"
#include "EbReferenceObject.h"

#include "EbInterPrediction.h"
#include "EbMcp.h"
#include "EbAvcStyleMcp.h"
#include "EbAdaptiveMotionVectorPrediction.h"

#include "EbModeDecisionProcess.h"
#include "EbErrorCodes.h"
#include <smmintrin.h> /* SSE4.1 */


#include "convolve.h"
#include "aom_dsp_rtcd.h"

#define MVBOUNDLOW    36    //  (80-71)<<2 // 80 = ReferencePadding ; minus 71 is derived from the expression -64 + 1 - 8, and plus 7 is derived from expression -1 + 8
#define MVBOUNDHIGH   348   //  (80+7)<<2
#define REFPADD_QPEL  320   //  (16+64)<<2


#define UNUSED_VARIABLE(x) (void)(x)

#define AOM_INTERP_EXTEND 4

#define SCALE_SUBPEL_BITS 10
#define SCALE_SUBPEL_SHIFTS (1 << SCALE_SUBPEL_BITS)
#define SCALE_SUBPEL_MASK (SCALE_SUBPEL_SHIFTS - 1)
#define SCALE_EXTRA_BITS (SCALE_SUBPEL_BITS - SUBPEL_BITS)
#define SCALE_EXTRA_OFF ((1 << SCALE_EXTRA_BITS) / 2)

#define RS_SUBPEL_BITS 6
#define RS_SUBPEL_MASK ((1 << RS_SUBPEL_BITS) - 1)
#define RS_SCALE_SUBPEL_BITS 14
#define RS_SCALE_SUBPEL_MASK ((1 << RS_SCALE_SUBPEL_BITS) - 1)
#define RS_SCALE_EXTRA_BITS (RS_SCALE_SUBPEL_BITS - RS_SUBPEL_BITS)

#define BIL_SUBPEL_BITS 3
#define BIL_SUBPEL_SHIFTS (1 << BIL_SUBPEL_BITS)

#define ROUND0_BITS 3
#define COMPOUND_ROUND1_BITS 7

typedef struct SubpelParams {
    int32_t xs;
    int32_t ys;
    int32_t subpel_x;
    int32_t subpel_y;
} SubpelParams;

//extern INLINE void clamp_mv(MV *mv, int32_t min_col, int32_t max_col, int32_t min_row,int32_t max_row);

static INLINE void clamp_mv(MV *mv, int32_t min_col, int32_t max_col, int32_t min_row,
    int32_t max_row) {
    mv->col = (int16_t)clamp(mv->col, min_col, max_col);
    mv->row = (int16_t)clamp(mv->row, min_row, max_row);
}

extern void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);

static INLINE MV clamp_mv_to_umv_border_sb(const MacroBlockD *xd,
    const MV *src_mv, int32_t bw, int32_t bh,
    int32_t ss_x, int32_t ss_y) {
    // If the MV points so far into the UMV border that no visible pixels
    // are used for reconstruction, the subpel part of the MV can be
    // discarded and the MV limited to 16 pixels with equivalent results.
    const int32_t spel_left = (AOM_INTERP_EXTEND + bw) << SUBPEL_BITS;
    const int32_t spel_right = spel_left - SUBPEL_SHIFTS;
    const int32_t spel_top = (AOM_INTERP_EXTEND + bh) << SUBPEL_BITS;
    const int32_t spel_bottom = spel_top - SUBPEL_SHIFTS;
    MV clamped_mv = { (int16_t)(src_mv->row * (1 << (1 - ss_y))),
        (int16_t)(src_mv->col * (1 << (1 - ss_x))) };
    assert(ss_x <= 1);
    assert(ss_y <= 1);

    clamp_mv(&clamped_mv,
        xd->mb_to_left_edge   * (1 << (1 - ss_x)) - spel_left,
        xd->mb_to_right_edge  * (1 << (1 - ss_x)) + spel_right,
        xd->mb_to_top_edge    * (1 << (1 - ss_y)) - spel_top,
        xd->mb_to_bottom_edge * (1 << (1 - ss_y)) + spel_bottom);

    return clamped_mv;
}



DECLARE_ALIGNED(256, static const InterpKernel,
sub_pel_filters_8[SUBPEL_SHIFTS]) = {
    { 0, 0, 0, 128, 0, 0, 0, 0 },{ 0, 2, -6, 126, 8, -2, 0, 0 },
    { 0, 2, -10, 122, 18, -4, 0, 0 },{ 0, 2, -12, 116, 28, -8, 2, 0 },
    { 0, 2, -14, 110, 38, -10, 2, 0 },{ 0, 2, -14, 102, 48, -12, 2, 0 },
    { 0, 2, -16, 94, 58, -12, 2, 0 },{ 0, 2, -14, 84, 66, -12, 2, 0 },
    { 0, 2, -14, 76, 76, -14, 2, 0 },{ 0, 2, -12, 66, 84, -14, 2, 0 },
    { 0, 2, -12, 58, 94, -16, 2, 0 },{ 0, 2, -12, 48, 102, -14, 2, 0 },
    { 0, 2, -10, 38, 110, -14, 2, 0 },{ 0, 2, -8, 28, 116, -12, 2, 0 },
    { 0, 0, -4, 18, 122, -10, 2, 0 },{ 0, 0, -2, 8, 126, -6, 2, 0 }
};
DECLARE_ALIGNED(256, static const InterpKernel,
sub_pel_filters_4[SUBPEL_SHIFTS]) = {
    { 0, 0, 0, 128, 0, 0, 0, 0 },{ 0, 0, -4, 126, 8, -2, 0, 0 },
    { 0, 0, -8, 122, 18, -4, 0, 0 },{ 0, 0, -10, 116, 28, -6, 0, 0 },
    { 0, 0, -12, 110, 38, -8, 0, 0 },{ 0, 0, -12, 102, 48, -10, 0, 0 },
    { 0, 0, -14, 94, 58, -10, 0, 0 },{ 0, 0, -12, 84, 66, -10, 0, 0 },
    { 0, 0, -12, 76, 76, -12, 0, 0 },{ 0, 0, -10, 66, 84, -12, 0, 0 },
    { 0, 0, -10, 58, 94, -14, 0, 0 },{ 0, 0, -10, 48, 102, -12, 0, 0 },
    { 0, 0, -8, 38, 110, -12, 0, 0 },{ 0, 0, -6, 28, 116, -10, 0, 0 },
    { 0, 0, -4, 18, 122, -8, 0, 0 },{ 0, 0, -2, 8, 126, -4, 0, 0 }
};

#define MAX_FILTER_TAP 8

void av1_convolve_2d_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst,
    int32_t dst_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params)
{
    int16_t im_block[(MAX_SB_SIZE + MAX_FILTER_TAP - 1) * MAX_SB_SIZE];
    int32_t im_h = h + filter_params_y->taps - 1;
    int32_t im_stride = w;
    const int32_t fo_vert = filter_params_y->taps / 2 - 1;
    const int32_t fo_horiz = filter_params_x->taps / 2 - 1;
    const int32_t bd = 8;
    const int32_t bits =
        FILTER_BITS * 2 - conv_params->round_0 - conv_params->round_1;

    // horizontal filter
    const uint8_t *src_horiz = src - fo_vert * src_stride;
    const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < im_h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t sum = (1 << (bd + FILTER_BITS - 1));
            for (int32_t k = 0; k < filter_params_x->taps; ++k) {
                sum += x_filter[k] * src_horiz[y * src_stride + x - fo_horiz + k];
            }
            assert(0 <= sum && sum < (1 << (bd + FILTER_BITS + 1)));
            im_block[y * im_stride + x] =
                (int16_t)ROUND_POWER_OF_TWO(sum, conv_params->round_0);
        }
    }

    // vertical filter
    int16_t *src_vert = im_block + fo_vert * im_stride;
    const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_y, subpel_y_q4 & SUBPEL_MASK);
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t sum = 1 << offset_bits;
            for (int32_t k = 0; k < filter_params_y->taps; ++k) {
                sum += y_filter[k] * src_vert[(y - fo_vert + k) * im_stride + x];
            }
            assert(0 <= sum && sum < (1 << (offset_bits + 2)));
            int16_t res = (CONV_BUF_TYPE)(ROUND_POWER_OF_TWO(sum, conv_params->round_1) -
                ((1 << (offset_bits - conv_params->round_1)) +
                (1 << (offset_bits - conv_params->round_1 - 1))));
            dst[y * dst_stride + x] = (uint8_t)clip_pixel_highbd(ROUND_POWER_OF_TWO(res, bits), 8);
        }
    }
}

void av1_convolve_y_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst,
    int32_t dst_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params)
{
    const int32_t fo_vert = filter_params_y->taps / 2 - 1;
    (void)filter_params_x;
    (void)subpel_x_q4;
    (void)conv_params;

    assert(conv_params->round_0 <= FILTER_BITS);
    assert(((conv_params->round_0 + conv_params->round_1) <= (FILTER_BITS + 1)) ||
        ((conv_params->round_0 + conv_params->round_1) == (2 * FILTER_BITS)));

    // vertical filter
    const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_y, subpel_y_q4 & SUBPEL_MASK);

    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t res = 0;
            for (int32_t k = 0; k < filter_params_y->taps; ++k) {
                res += y_filter[k] * src[(y - fo_vert + k) * src_stride + x];
            }
            dst[y * dst_stride + x] =
                (uint8_t)clip_pixel_highbd(ROUND_POWER_OF_TWO(res, FILTER_BITS), 8);
        }
    }
}

void av1_convolve_x_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst,
    int32_t dst_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params)
{
    const int32_t fo_horiz = filter_params_x->taps / 2 - 1;
    const int32_t bits = FILTER_BITS - conv_params->round_0;
    (void)filter_params_y;
    (void)subpel_y_q4;
    (void)conv_params;

    assert(bits >= 0);
    assert((FILTER_BITS - conv_params->round_1) >= 0 ||
        ((conv_params->round_0 + conv_params->round_1) == 2 * FILTER_BITS));

    // horizontal filter
    const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_x, subpel_x_q4 & SUBPEL_MASK);

    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t res = 0;
            for (int32_t k = 0; k < filter_params_x->taps; ++k) {
                res += x_filter[k] * src[y * src_stride + x - fo_horiz + k];
            }
            res = ROUND_POWER_OF_TWO(res, conv_params->round_0);
            dst[y * dst_stride + x] = (uint8_t)clip_pixel_highbd(ROUND_POWER_OF_TWO(res, bits), 8);
        }
    }
}

void av1_convolve_2d_copy_sr_c(const uint8_t *src, int32_t src_stride, uint8_t *dst,
    int32_t dst_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params) {
    (void)filter_params_x;
    (void)filter_params_y;
    (void)subpel_x_q4;
    (void)subpel_y_q4;
    (void)conv_params;

    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            dst[y * dst_stride + x] = src[y * src_stride + x];
        }
    }
}

void av1_jnt_convolve_2d_c(const uint8_t *src, int32_t src_stride, uint8_t *dst8,
    int32_t dst8_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params)
{
    CONV_BUF_TYPE *dst = conv_params->dst;
    int32_t dst_stride = conv_params->dst_stride;
    int16_t im_block[(MAX_SB_SIZE + MAX_FILTER_TAP - 1) * MAX_SB_SIZE];
    int32_t im_h = h + filter_params_y->taps - 1;
    int32_t im_stride = w;
    const int32_t fo_vert = filter_params_y->taps / 2 - 1;
    const int32_t fo_horiz = filter_params_x->taps / 2 - 1;
    const int32_t bd = 8;
    const int32_t round_bits =
        2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;

    // horizontal filter
    const uint8_t *src_horiz = src - fo_vert * src_stride;
    const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < im_h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t sum = (1 << (bd + FILTER_BITS - 1));
            for (int32_t k = 0; k < filter_params_x->taps; ++k) {
                sum += x_filter[k] * src_horiz[y * src_stride + x - fo_horiz + k];
            }
            assert(0 <= sum && sum < (1 << (bd + FILTER_BITS + 1)));
            im_block[y * im_stride + x] =
                (int16_t)ROUND_POWER_OF_TWO(sum, conv_params->round_0);
        }
    }

    // vertical filter
    int16_t *src_vert = im_block + fo_vert * im_stride;
    const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_y, subpel_y_q4 & SUBPEL_MASK);
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t sum = 1 << offset_bits;
            for (int32_t k = 0; k < filter_params_y->taps; ++k) {
                sum += y_filter[k] * src_vert[(y - fo_vert + k) * im_stride + x];
            }
            assert(0 <= sum && sum < (1 << (offset_bits + 2)));
            CONV_BUF_TYPE res = (CONV_BUF_TYPE)ROUND_POWER_OF_TWO(sum, conv_params->round_1);
            if (conv_params->do_average) {
                int32_t tmp = dst[y * dst_stride + x];
                if (conv_params->use_jnt_comp_avg) {
                    tmp = tmp * conv_params->fwd_offset + res * conv_params->bck_offset;
                    tmp = tmp >> DIST_PRECISION_BITS;
                    printf("here");
                }
                else {
                    tmp += res;
                    tmp = tmp >> 1;
                }
                tmp -= (1 << (offset_bits - conv_params->round_1)) +
                    (1 << (offset_bits - conv_params->round_1 - 1));
                dst8[y * dst8_stride + x] =
                    (uint8_t)clip_pixel_highbd(ROUND_POWER_OF_TWO(tmp, round_bits), 8);
            }
            else {
                dst[y * dst_stride + x] = res;
            }

        }
    }
}

void av1_jnt_convolve_y_c(const uint8_t *src, int32_t src_stride, uint8_t *dst8,
    int32_t dst8_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params)
{
    CONV_BUF_TYPE *dst = conv_params->dst;
    int32_t dst_stride = conv_params->dst_stride;
    const int32_t fo_vert = filter_params_y->taps / 2 - 1;
    const int32_t bits = FILTER_BITS - conv_params->round_0;
    const int32_t bd = 8;
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    const int32_t round_offset = (1 << (offset_bits - conv_params->round_1)) +
        (1 << (offset_bits - conv_params->round_1 - 1));
    const int32_t round_bits =
        2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    (void)filter_params_x;
    (void)subpel_x_q4;

    // vertical filter
    const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_y, subpel_y_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t res = 0;
            for (int32_t k = 0; k < filter_params_y->taps; ++k) {
                res += y_filter[k] * src[(y - fo_vert + k) * src_stride + x];
            }
            res *= (1 << bits);
            res = ROUND_POWER_OF_TWO(res, conv_params->round_1) + round_offset;

            if (conv_params->do_average) {
                int32_t tmp = dst[y * dst_stride + x];
                if (conv_params->use_jnt_comp_avg) {
                    tmp = tmp * conv_params->fwd_offset + res * conv_params->bck_offset;
                    tmp = tmp >> DIST_PRECISION_BITS;
                }
                else {
                    tmp += res;
                    tmp = tmp >> 1;
                }
                tmp -= round_offset;
                dst8[y * dst8_stride + x] =
                    (uint8_t)clip_pixel_highbd(ROUND_POWER_OF_TWO(tmp, round_bits), 8);
            }
            else {
                dst[y * dst_stride + x] = (CONV_BUF_TYPE)res;
            }

        }
    }
}

void av1_jnt_convolve_x_c(const uint8_t *src, int32_t src_stride, uint8_t *dst8,
    int32_t dst8_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params)
{
    CONV_BUF_TYPE *dst = conv_params->dst;
    int32_t dst_stride = conv_params->dst_stride;
    const int32_t fo_horiz = filter_params_x->taps / 2 - 1;
    const int32_t bits = FILTER_BITS - conv_params->round_1;
    const int32_t bd = 8;
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    const int32_t round_offset = (1 << (offset_bits - conv_params->round_1)) +
        (1 << (offset_bits - conv_params->round_1 - 1));
    const int32_t round_bits =
        2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    (void)filter_params_y;
    (void)subpel_y_q4;

    // horizontal filter
    const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t res = 0;
            for (int32_t k = 0; k < filter_params_x->taps; ++k) {
                res += x_filter[k] * src[y * src_stride + x - fo_horiz + k];
            }
            res = (1 << bits) * ROUND_POWER_OF_TWO(res, conv_params->round_0);
            res += round_offset;

            if (conv_params->do_average) {
                int32_t tmp = dst[y * dst_stride + x];
                if (conv_params->use_jnt_comp_avg) {

                    tmp = tmp * conv_params->fwd_offset + res * conv_params->bck_offset;
                    tmp = tmp >> DIST_PRECISION_BITS;
                }
                else {
                    tmp += res;
                    tmp = tmp >> 1;
                }
                tmp -= round_offset;
                dst8[y * dst8_stride + x] =
                    (uint8_t)clip_pixel_highbd(ROUND_POWER_OF_TWO(tmp, round_bits), 8);
            }
            else {
                dst[y * dst_stride + x] = (CONV_BUF_TYPE)res;
            }
        }
    }
}

void av1_jnt_convolve_2d_copy_c(const uint8_t *src, int32_t src_stride,
    uint8_t *dst8, int32_t dst8_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params)
{
    CONV_BUF_TYPE *dst = conv_params->dst;
    int32_t dst_stride = conv_params->dst_stride;
    const int32_t bits =
        FILTER_BITS * 2 - conv_params->round_1 - conv_params->round_0;
    const int32_t bd = 8;
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    const int32_t round_offset = (1 << (offset_bits - conv_params->round_1)) +
        (1 << (offset_bits - conv_params->round_1 - 1));
    (void)filter_params_x;
    (void)filter_params_y;
    (void)subpel_x_q4;
    (void)subpel_y_q4;

    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            CONV_BUF_TYPE res = src[y * src_stride + x] << bits;
            res += (CONV_BUF_TYPE)round_offset;

            if (conv_params->do_average) {
                int32_t tmp = dst[y * dst_stride + x];
                if (conv_params->use_jnt_comp_avg) {
                    tmp = tmp * conv_params->fwd_offset + res * conv_params->bck_offset;
                    tmp = tmp >> DIST_PRECISION_BITS;
                }
                else {
                    tmp += res;
                    tmp = tmp >> 1;
                }
                tmp -= round_offset;
                dst8[y * dst8_stride + x] = (uint8_t)clip_pixel_highbd(ROUND_POWER_OF_TWO(tmp, bits), 8);
            }
            else {
                dst[y * dst_stride + x] = res;
            }
        }
    }
}

void av1_highbd_convolve_2d_copy_sr_c(
    const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w,
    int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd) {
    (void)filter_params_x;
    (void)filter_params_y;
    (void)subpel_x_q4;
    (void)subpel_y_q4;
    (void)conv_params;
    (void)bd;

    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            dst[y * dst_stride + x] = src[y * src_stride + x];
        }
    }
}

void av1_highbd_convolve_x_sr_c(const uint16_t *src, int32_t src_stride,
    uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h,
    const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params, int32_t bd) {
    const int32_t fo_horiz = filter_params_x->taps / 2 - 1;
    const int32_t bits = FILTER_BITS - conv_params->round_0;
    (void)filter_params_y;
    (void)subpel_y_q4;

    assert(bits >= 0);
    assert((FILTER_BITS - conv_params->round_1) >= 0 ||
        ((conv_params->round_0 + conv_params->round_1) == 2 * FILTER_BITS));

    // horizontal filter
    const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t res = 0;
            for (int32_t k = 0; k < filter_params_x->taps; ++k) {
                res += x_filter[k] * src[y * src_stride + x - fo_horiz + k];
            }
            res = ROUND_POWER_OF_TWO(res, conv_params->round_0);
            dst[y * dst_stride + x] =
                clip_pixel_highbd(ROUND_POWER_OF_TWO(res, bits), bd);
        }
    }
}

void av1_highbd_convolve_y_sr_c(const uint16_t *src, int32_t src_stride,
    uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h,
    const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params, int32_t bd) {
    const int32_t fo_vert = filter_params_y->taps / 2 - 1;
    (void)filter_params_x;
    (void)subpel_x_q4;
    (void)conv_params;

    assert(conv_params->round_0 <= FILTER_BITS);
    assert(((conv_params->round_0 + conv_params->round_1) <= (FILTER_BITS + 1)) ||
        ((conv_params->round_0 + conv_params->round_1) == (2 * FILTER_BITS)));
    // vertical filter
    const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_y, subpel_y_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t res = 0;
            for (int32_t k = 0; k < filter_params_y->taps; ++k) {
                res += y_filter[k] * src[(y - fo_vert + k) * src_stride + x];
            }
            dst[y * dst_stride + x] =
                clip_pixel_highbd(ROUND_POWER_OF_TWO(res, FILTER_BITS), bd);
        }
    }
}

void av1_highbd_convolve_2d_sr_c(const uint16_t *src, int32_t src_stride,
    uint16_t *dst, int32_t dst_stride, int32_t w, int32_t h,
    const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params, int32_t bd) {
    int16_t im_block[(MAX_SB_SIZE + MAX_FILTER_TAP - 1) * MAX_SB_SIZE];
    int32_t im_h = h + filter_params_y->taps - 1;
    int32_t im_stride = w;
    const int32_t fo_vert = filter_params_y->taps / 2 - 1;
    const int32_t fo_horiz = filter_params_x->taps / 2 - 1;
    const int32_t bits =
        FILTER_BITS * 2 - conv_params->round_0 - conv_params->round_1;
    assert(bits >= 0);

    // horizontal filter
    const uint16_t *src_horiz = src - fo_vert * src_stride;
    const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < im_h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t sum = (1 << (bd + FILTER_BITS - 1));
            for (int32_t k = 0; k < filter_params_x->taps; ++k) {
                sum += x_filter[k] * src_horiz[y * src_stride + x - fo_horiz + k];
            }
            assert(0 <= sum && sum < (1 << (bd + FILTER_BITS + 1)));
            im_block[y * im_stride + x] = (CONV_BUF_TYPE)
                ROUND_POWER_OF_TWO(sum, conv_params->round_0);
        }
    }

    // vertical filter
    int16_t *src_vert = im_block + fo_vert * im_stride;
    const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_y, subpel_y_q4 & SUBPEL_MASK);
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t sum = 1 << offset_bits;
            for (int32_t k = 0; k < filter_params_y->taps; ++k) {
                sum += y_filter[k] * src_vert[(y - fo_vert + k) * im_stride + x];
            }
            assert(0 <= sum && sum < (1 << (offset_bits + 2)));
            int32_t res = ROUND_POWER_OF_TWO(sum, conv_params->round_1) -
                ((1 << (offset_bits - conv_params->round_1)) +
                (1 << (offset_bits - conv_params->round_1 - 1)));
            dst[y * dst_stride + x] =
                clip_pixel_highbd(ROUND_POWER_OF_TWO(res, bits), bd);
        }
    }
}

void av1_highbd_jnt_convolve_x_c(const uint16_t *src, int32_t src_stride,
    uint16_t *dst16, int32_t dst16_stride, int32_t w,
    int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params, int32_t bd) {
    CONV_BUF_TYPE *dst = conv_params->dst;
    int32_t dst_stride = conv_params->dst_stride;
    const int32_t fo_horiz = filter_params_x->taps / 2 - 1;
    const int32_t bits = FILTER_BITS - conv_params->round_1;
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    const int32_t round_offset = (1 << (offset_bits - conv_params->round_1)) +
        (1 << (offset_bits - conv_params->round_1 - 1));
    const int32_t round_bits =
        2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    assert(round_bits >= 0);
    (void)filter_params_y;
    (void)subpel_y_q4;
    assert(bits >= 0);
    // horizontal filter
    const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t res = 0;
            for (int32_t k = 0; k < filter_params_x->taps; ++k) {
                res += x_filter[k] * src[y * src_stride + x - fo_horiz + k];
            }
            res = (1 << bits) * ROUND_POWER_OF_TWO(res, conv_params->round_0);
            res += round_offset;

            if (conv_params->do_average) {
                int32_t tmp = dst[y * dst_stride + x];
                if (conv_params->use_jnt_comp_avg) {
                    tmp = tmp * conv_params->fwd_offset + res * conv_params->bck_offset;
                    tmp = tmp >> DIST_PRECISION_BITS;
                }
                else {
                    tmp += res;
                    tmp = tmp >> 1;
                }
                tmp -= round_offset;
                dst16[y * dst16_stride + x] =
                    clip_pixel_highbd(ROUND_POWER_OF_TWO(tmp, round_bits), bd);
            }
            else {
                dst[y * dst_stride + x] = (CONV_BUF_TYPE)res;
            }
        }
    }
}

void av1_highbd_jnt_convolve_y_c(const uint16_t *src, int32_t src_stride,
    uint16_t *dst16, int32_t dst16_stride, int32_t w,
    int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params, int32_t bd) {
    CONV_BUF_TYPE *dst = conv_params->dst;
    int32_t dst_stride = conv_params->dst_stride;
    const int32_t fo_vert = filter_params_y->taps / 2 - 1;
    const int32_t bits = FILTER_BITS - conv_params->round_0;
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    const int32_t round_offset = (1 << (offset_bits - conv_params->round_1)) +
        (1 << (offset_bits - conv_params->round_1 - 1));
    const int32_t round_bits =
        2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    assert(round_bits >= 0);
    (void)filter_params_x;
    (void)subpel_x_q4;
    assert(bits >= 0);
    // vertical filter
    const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_y, subpel_y_q4 & SUBPEL_MASK);
    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            int32_t res = 0;
            for (int32_t k = 0; k < filter_params_y->taps; ++k) {
                res += y_filter[k] * src[(y - fo_vert + k) * src_stride + x];
            }
            res *= (1 << bits);
            res = ROUND_POWER_OF_TWO(res, conv_params->round_1) + round_offset;

            if (conv_params->do_average) {
                int32_t tmp = dst[y * dst_stride + x];
                if (conv_params->use_jnt_comp_avg) {
                    tmp = tmp * conv_params->fwd_offset + res * conv_params->bck_offset;
                    tmp = tmp >> DIST_PRECISION_BITS;
                }
                else {
                    tmp += res;
                    tmp = tmp >> 1;
                }
                tmp -= round_offset;
                dst16[y * dst16_stride + x] =
                    clip_pixel_highbd(ROUND_POWER_OF_TWO(tmp, round_bits), bd);
            }
            else {
                dst[y * dst_stride + x] = (CONV_BUF_TYPE)res;
            }
        }
    }
}

void av1_highbd_jnt_convolve_2d_copy_c(
    const uint16_t *src, int32_t src_stride, uint16_t *dst16, int32_t dst16_stride,
    int32_t w, int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
    const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd) {
    CONV_BUF_TYPE *dst = conv_params->dst;
    int32_t dst_stride = conv_params->dst_stride;
    const int32_t bits =
        FILTER_BITS * 2 - conv_params->round_1 - conv_params->round_0;
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    const int32_t round_offset = (1 << (offset_bits - conv_params->round_1)) +
        (1 << (offset_bits - conv_params->round_1 - 1));
    assert(bits >= 0);
    (void)filter_params_x;
    (void)filter_params_y;
    (void)subpel_x_q4;
    (void)subpel_y_q4;

    for (int32_t y = 0; y < h; ++y) {
        for (int32_t x = 0; x < w; ++x) {
            CONV_BUF_TYPE res = src[y * src_stride + x] << bits;
            res += (CONV_BUF_TYPE)round_offset;
            if (conv_params->do_average) {
                int32_t tmp = dst[y * dst_stride + x];
                if (conv_params->use_jnt_comp_avg) {
                    tmp = tmp * conv_params->fwd_offset + res * conv_params->bck_offset;
                    tmp = tmp >> DIST_PRECISION_BITS;
                }
                else {
                    tmp += res;
                    tmp = tmp >> 1;
                }
                tmp -= round_offset;
                dst16[y * dst16_stride + x] =
                    clip_pixel_highbd(ROUND_POWER_OF_TWO(tmp, bits), bd);
            }
            else {
                dst[y * dst_stride + x] = res;
            }
        }
    }
}

void av1_highbd_jnt_convolve_2d_c(const uint16_t *src, int32_t src_stride,
    uint16_t *dst16, int32_t dst16_stride, int32_t w,
    int32_t h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params, int32_t bd)

{
    int32_t x, y, k;
    int16_t im_block[(MAX_SB_SIZE + MAX_FILTER_TAP - 1) * MAX_SB_SIZE];
    CONV_BUF_TYPE *dst = conv_params->dst;
    int32_t dst_stride = conv_params->dst_stride;
    int32_t im_h = h + filter_params_y->taps - 1;
    int32_t im_stride = w;
    const int32_t fo_vert = filter_params_y->taps / 2 - 1;
    const int32_t fo_horiz = filter_params_x->taps / 2 - 1;

    const int32_t round_bits =
        2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;
    assert(round_bits >= 0);


    // horizontal filter
    const uint16_t *src_horiz = src - fo_vert * src_stride;
    const int16_t *x_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_x, subpel_x_q4 & SUBPEL_MASK);
    for (y = 0; y < im_h; ++y) {
        for (x = 0; x < w; ++x) {
            int32_t sum = (1 << (bd + FILTER_BITS - 1));
            for (k = 0; k < filter_params_x->taps; ++k) {
                sum += x_filter[k] * src_horiz[y * src_stride + x - fo_horiz + k];
            }
            assert(0 <= sum && sum < (1 << (bd + FILTER_BITS + 1)));
            (void)bd;
            im_block[y * im_stride + x] =
                (int16_t)ROUND_POWER_OF_TWO(sum, conv_params->round_0);
        }
    }

    // vertical filter
    int16_t *src_vert = im_block + fo_vert * im_stride;
    const int32_t offset_bits = bd + 2 * FILTER_BITS - conv_params->round_0;
    const int16_t *y_filter = av1_get_interp_filter_subpel_kernel(
        *filter_params_y, subpel_y_q4 & SUBPEL_MASK);
    for (y = 0; y < h; ++y) {
        for (x = 0; x < w; ++x) {
            int32_t sum = 1 << offset_bits;
            for (k = 0; k < filter_params_y->taps; ++k) {
                sum += y_filter[k] * src_vert[(y - fo_vert + k) * im_stride + x];
            }
            assert(0 <= sum && sum < (1 << (offset_bits + 2)));
            CONV_BUF_TYPE res = (CONV_BUF_TYPE)ROUND_POWER_OF_TWO(sum, conv_params->round_1);
            if (conv_params->do_average) {
                int32_t tmp = dst[y * dst_stride + x];
                if (conv_params->use_jnt_comp_avg) {
                    tmp = tmp * conv_params->fwd_offset + res * conv_params->bck_offset;
                    tmp = tmp >> DIST_PRECISION_BITS;
                }
                else {
                    tmp += res;
                    tmp = tmp >> 1;
                }
                tmp -= (1 << (offset_bits - conv_params->round_1)) +
                    (1 << (offset_bits - conv_params->round_1 - 1));
                dst16[y * dst16_stride + x] =
                    clip_pixel_highbd(ROUND_POWER_OF_TWO(tmp, round_bits), bd);
            }
            else {
                dst[y * dst_stride + x] = res;
            }
        }
    }
}

aom_highbd_convolve_fn_t convolveHbd[/*subX*/2][/*subY*/2][/*bi*/2];
void asmSetConvolveHbdAsmTable(void)
{
    convolveHbd[0][0][0] = av1_highbd_convolve_2d_copy_sr;
    convolveHbd[0][0][1] = av1_highbd_jnt_convolve_2d_copy;

    convolveHbd[0][1][0] = av1_highbd_convolve_y_sr;
    convolveHbd[0][1][1] = av1_highbd_jnt_convolve_y;

    convolveHbd[1][0][0] = av1_highbd_convolve_x_sr;
    convolveHbd[1][0][1] = av1_highbd_jnt_convolve_x;

    convolveHbd[1][1][0] = av1_highbd_convolve_2d_sr;
    convolveHbd[1][1][1] = av1_highbd_jnt_convolve_2d;
}

aom_convolve_fn_t convolve[/*subX*/2][/*subY*/2][/*bi*/2];
void asmSetConvolveAsmTable(void)
{
    convolve[0][0][0] = av1_convolve_2d_copy_sr;
    convolve[0][0][1] = av1_jnt_convolve_2d_copy;

    convolve[0][1][0] = av1_convolve_y_sr;
    convolve[0][1][1] = av1_jnt_convolve_y;

    convolve[1][0][0] = av1_convolve_x_sr;
    convolve[1][0][1] = av1_jnt_convolve_x;

    convolve[1][1][0] = av1_convolve_2d_sr;
    convolve[1][1][1] = av1_jnt_convolve_2d;
}

InterpFilterParams av1RegularFilter = { (const int16_t *)sub_pel_filters_8, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_REGULAR };
InterpFilterParams av1RegularFilterW4 = { (const int16_t *)sub_pel_filters_4, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_REGULAR };

DECLARE_ALIGNED(256, static const InterpKernel,
sub_pel_filters_8sharp[SUBPEL_SHIFTS]) = {
{ 0, 0, 0, 128, 0, 0, 0, 0 },         { -2, 2, -6, 126, 8, -2, 2, 0 },
{ -2, 6, -12, 124, 16, -6, 4, -2 },   { -2, 8, -18, 120, 26, -10, 6, -2 },
{ -4, 10, -22, 116, 38, -14, 6, -2 }, { -4, 10, -22, 108, 48, -18, 8, -2 },
{ -4, 10, -24, 100, 60, -20, 8, -2 }, { -4, 10, -24, 90, 70, -22, 10, -2 },
{ -4, 12, -24, 80, 80, -24, 12, -4 }, { -2, 10, -22, 70, 90, -24, 10, -4 },
{ -2, 8, -20, 60, 100, -24, 10, -4 }, { -2, 8, -18, 48, 108, -22, 10, -4 },
{ -2, 6, -14, 38, 116, -22, 10, -4 }, { -2, 6, -10, 26, 120, -18, 8, -2 },
{ -2, 4, -6, 16, 124, -12, 6, -2 },   { 0, 2, -2, 8, 126, -6, 2, -2 }
};

DECLARE_ALIGNED(256, static const InterpKernel,
sub_pel_filters_8smooth[SUBPEL_SHIFTS]) = {
{ 0, 0, 0, 128, 0, 0, 0, 0 },     { 0, 2, 28, 62, 34, 2, 0, 0 },
{ 0, 0, 26, 62, 36, 4, 0, 0 },    { 0, 0, 22, 62, 40, 4, 0, 0 },
{ 0, 0, 20, 60, 42, 6, 0, 0 },    { 0, 0, 18, 58, 44, 8, 0, 0 },
{ 0, 0, 16, 56, 46, 10, 0, 0 },   { 0, -2, 16, 54, 48, 12, 0, 0 },
{ 0, -2, 14, 52, 52, 14, -2, 0 }, { 0, 0, 12, 48, 54, 16, -2, 0 },
{ 0, 0, 10, 46, 56, 16, 0, 0 },   { 0, 0, 8, 44, 58, 18, 0, 0 },
{ 0, 0, 6, 42, 60, 20, 0, 0 },    { 0, 0, 4, 40, 62, 22, 0, 0 },
{ 0, 0, 4, 36, 62, 26, 0, 0 },    { 0, 0, 2, 34, 62, 28, 2, 0 }
};
DECLARE_ALIGNED(256, static const InterpKernel,
bilinear_filters[SUBPEL_SHIFTS]) = {
{ 0, 0, 0, 128, 0, 0, 0, 0 },  { 0, 0, 0, 120, 8, 0, 0, 0 },
{ 0, 0, 0, 112, 16, 0, 0, 0 }, { 0, 0, 0, 104, 24, 0, 0, 0 },
{ 0, 0, 0, 96, 32, 0, 0, 0 },  { 0, 0, 0, 88, 40, 0, 0, 0 },
{ 0, 0, 0, 80, 48, 0, 0, 0 },  { 0, 0, 0, 72, 56, 0, 0, 0 },
{ 0, 0, 0, 64, 64, 0, 0, 0 },  { 0, 0, 0, 56, 72, 0, 0, 0 },
{ 0, 0, 0, 48, 80, 0, 0, 0 },  { 0, 0, 0, 40, 88, 0, 0, 0 },
{ 0, 0, 0, 32, 96, 0, 0, 0 },  { 0, 0, 0, 24, 104, 0, 0, 0 },
{ 0, 0, 0, 16, 112, 0, 0, 0 }, { 0, 0, 0, 8, 120, 0, 0, 0 }
};
DECLARE_ALIGNED(256, static const InterpKernel,
sub_pel_filters_4smooth[SUBPEL_SHIFTS]) = {
{ 0, 0, 0, 128, 0, 0, 0, 0 },   { 0, 0, 30, 62, 34, 2, 0, 0 },
{ 0, 0, 26, 62, 36, 4, 0, 0 },  { 0, 0, 22, 62, 40, 4, 0, 0 },
{ 0, 0, 20, 60, 42, 6, 0, 0 },  { 0, 0, 18, 58, 44, 8, 0, 0 },
{ 0, 0, 16, 56, 46, 10, 0, 0 }, { 0, 0, 14, 54, 48, 12, 0, 0 },
{ 0, 0, 12, 52, 52, 12, 0, 0 }, { 0, 0, 12, 48, 54, 14, 0, 0 },
{ 0, 0, 10, 46, 56, 16, 0, 0 }, { 0, 0, 8, 44, 58, 18, 0, 0 },
{ 0, 0, 6, 42, 60, 20, 0, 0 },  { 0, 0, 4, 40, 62, 22, 0, 0 },
{ 0, 0, 4, 36, 62, 26, 0, 0 },  { 0, 0, 2, 34, 62, 30, 0, 0 }
};
static const InterpFilterParams
av1_interp_filter_params_list[SWITCHABLE_FILTERS + 1] = {
  { (const int16_t *)sub_pel_filters_8, SUBPEL_TAPS, SUBPEL_SHIFTS,
    EIGHTTAP_REGULAR },
  { (const int16_t *)sub_pel_filters_8smooth, SUBPEL_TAPS, SUBPEL_SHIFTS,
    EIGHTTAP_SMOOTH },
  { (const int16_t *)sub_pel_filters_8sharp, SUBPEL_TAPS, SUBPEL_SHIFTS,
    MULTITAP_SHARP },
  { (const int16_t *)bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS,
    BILINEAR }
};
static const InterpFilterParams av1_interp_4tap[2] = {
  { (const int16_t *)sub_pel_filters_4, SUBPEL_TAPS, SUBPEL_SHIFTS,
    EIGHTTAP_REGULAR },
  { (const int16_t *)sub_pel_filters_4smooth, SUBPEL_TAPS, SUBPEL_SHIFTS,
    EIGHTTAP_SMOOTH },
};
InterpFilterParams av1_get_interp_filter_params_with_block_size(
    const InterpFilter interp_filter, const int32_t w) {
    if (w <= 4 &&
        (interp_filter == MULTITAP_SHARP || interp_filter == EIGHTTAP_REGULAR))
        return av1_interp_4tap[0];
    else if (w <= 4 && interp_filter == EIGHTTAP_SMOOTH)
        return av1_interp_4tap[1];

    return av1_interp_filter_params_list[interp_filter];
}

static void av1_get_convolve_filter_params( uint32_t interp_filters,
    InterpFilterParams *params_x, InterpFilterParams *params_y,
    int32_t w, int32_t h)
{
    InterpFilter filter_x = av1_extract_interp_filter(interp_filters, 1);
    InterpFilter filter_y = av1_extract_interp_filter(interp_filters, 0);
    *params_x = av1_get_interp_filter_params_with_block_size(filter_x, w);
    *params_y = av1_get_interp_filter_params_with_block_size(filter_y, h);
}


int32_t is_inter_block(const MbModeInfo *mbmi);
block_size scale_chroma_bsize(block_size bsize, int32_t subsampling_x,
    int32_t subsampling_y);

EbErrorType av1_inter_prediction(
    PictureControlSet_t                    *picture_control_set_ptr,
    uint32_t                                interp_filters,
    CodingUnit_t                           *cu_ptr,
    uint8_t                                 ref_frame_type,
    MvUnit_t                               *mv_unit,
    uint16_t                                pu_origin_x,
    uint16_t                                pu_origin_y,
    uint8_t                                 bwidth,
    uint8_t                                 bheight,
    EbPictureBufferDesc_t                  *ref_pic_list0,
    EbPictureBufferDesc_t                  *ref_pic_list1,
    EbPictureBufferDesc_t                  *prediction_ptr,
    uint16_t                                dst_origin_x,
    uint16_t                                dst_origin_y,
#if CHROMA_BLIND
    EbBool                                  perform_chroma,
#endif
    EbAsm                                   asm_type)
{
    (void)asm_type;
    EbErrorType  return_error = EB_ErrorNone;
    uint8_t         is_compound = (mv_unit->predDirection == BI_PRED) ? 1 : 0;
    DECLARE_ALIGNED(32, uint16_t, tmp_dstY[128 * 128]);//move this to context if stack does not hold.
    DECLARE_ALIGNED(32, uint16_t, tmp_dstCb[64 * 64]);
    DECLARE_ALIGNED(32, uint16_t, tmp_dstCr[64 * 64]);

    MV  mv, mv_q4;

    int32_t subpel_x, subpel_y;
    uint8_t * src_ptr;
    uint8_t * dst_ptr;
    int32_t src_stride;
    int32_t dst_stride;
    ConvolveParams conv_params;

    InterpFilterParams filter_params_x, filter_params_y;

    const BlockGeom * blk_geom = get_blk_geom_mds(cu_ptr->mds_idx);


    //special treatment for chroma in 4XN/NX4 blocks
    //if one of the neighbour blocks of the parent square is intra the chroma prediction will follow the normal path using the luma MV of the current nsq block which is the latest sub8x8.
    //for this case: only uniPred is allowed.

    int32_t sub8x8_inter = 0;
#if CHROMA_BLIND
    if(perform_chroma && (blk_geom->has_uv && (blk_geom->bwidth == 4 || blk_geom->bheight == 4)))
#else
    if (blk_geom->has_uv &&
        (blk_geom->bwidth == 4 || blk_geom->bheight == 4)
        )
#endif
    {

        //CHKN setup input param

        int32_t bw = blk_geom->bwidth_uv;
        int32_t bh = blk_geom->bheight_uv;
        UNUSED_VARIABLE(bw);
        UNUSED_VARIABLE(bh);

        uint32_t mi_x = pu_origin_x;       //these are luma picture wise
        uint32_t mi_y = pu_origin_y;

        MacroBlockD  *xd = cu_ptr->av1xd;
        xd->mi_stride = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->picture_width_in_sb*(BLOCK_SIZE_64 / 4);
        const int32_t offset = (mi_y >> MI_SIZE_LOG2) * xd->mi_stride + (mi_x >> MI_SIZE_LOG2);
        xd->mi = picture_control_set_ptr->mi_grid_base + offset;

        //CHKN fill current mi from current block
        {
            ModeInfo *miPtr = *xd->mi;
            uint8_t  miX, miY;
            MvReferenceFrame rf[2];
            av1_set_ref_frame(rf, ref_frame_type);
            for (miY = 0; miY < (blk_geom->bheight >> MI_SIZE_LOG2); miY++) {
                for (miX = 0; miX < (blk_geom->bwidth >> MI_SIZE_LOG2); miX++) {
                    miPtr[miX + miY * xd->mi_stride].mbmi.ref_frame[0] = rf[0];
                    if (mv_unit->predDirection == UNI_PRED_LIST_0) {
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_0].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_0].y;
                    }
                    else if (mv_unit->predDirection == UNI_PRED_LIST_1) {
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_1].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_1].y;
                    }
                    else {

                        // printf("ERRRRRRR");

                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_0].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_0].y;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[1].as_mv.col = mv_unit->mv[REF_LIST_1].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[1].as_mv.row = mv_unit->mv[REF_LIST_1].y;
                    }

                }
            }
        }

        int32_t build_for_obmc = 0;

        const block_size bsize = blk_geom->bsize;//mi->sb_type;
        ASSERT(bsize < BlockSizeS_ALL);
        const int32_t ss_x = 1;// pd->subsampling_x;
        const int32_t ss_y = 1;//pd->subsampling_y;
        sub8x8_inter = (block_size_wide[bsize] < 8 && ss_x) ||
            (block_size_high[bsize] < 8 && ss_y);

        //if (is_intrabc) sub8x8_inter = 0;

        // For sub8x8 chroma blocks, we may be covering more than one luma block's
        // worth of pixels. Thus (mi_x, mi_y) may not be the correct coordinates for
        // the top-left corner of the prediction source - the correct top-left corner
        // is at (pre_x, pre_y).
        const int32_t row_start =
            (block_size_high[bsize] == 4) && ss_y && !build_for_obmc ? -1 : 0;
        const int32_t col_start =
            (block_size_wide[bsize] == 4) && ss_x && !build_for_obmc ? -1 : 0;

        const int32_t pre_x = (mi_x + MI_SIZE * col_start) >> ss_x;
        const int32_t pre_y = (mi_y + MI_SIZE * row_start) >> ss_y;
        UNUSED_VARIABLE(pre_x);
        UNUSED_VARIABLE(pre_y);

        sub8x8_inter = sub8x8_inter && !build_for_obmc;
        if (sub8x8_inter) {
            for (int32_t row = row_start; row <= 0 && sub8x8_inter; ++row) {
                for (int32_t col = col_start; col <= 0; ++col) {

                    ModeInfo *miPtr = *xd->mi;
                    const MbModeInfo *this_mbmi = &miPtr[row * xd->mi_stride + col].mbmi;

                    if (!is_inter_block(this_mbmi)) sub8x8_inter = 0;
                    //if (is_intrabc_block(this_mbmi)) sub8x8_inter = 0;
                }
            }
        }

        if (sub8x8_inter) {
            // block size
            const int32_t b4_w = block_size_wide[bsize] >> ss_x;
            const int32_t b4_h = block_size_high[bsize] >> ss_y;
            const block_size plane_bsize = scale_chroma_bsize(bsize, ss_x, ss_y);
            ASSERT(plane_bsize < BlockSizeS_ALL);
            const int32_t b8_w = block_size_wide[plane_bsize] >> ss_x;
            const int32_t b8_h = block_size_high[plane_bsize] >> ss_y;

            assert(!is_compound);

            if (is_compound)
                printf("ETTTT");

            //const struct Buf2d orig_pred_buf[2] = { pd->pre[0], pd->pre[1] };

            int32_t row = row_start;
            int32_t src_stride;
            for (int32_t y = 0; y < b8_h; y += b4_h) {
                int32_t col = col_start;
                for (int32_t x = 0; x < b8_w; x += b4_w) {

                    ModeInfo *miPtr = *xd->mi;
                    const MbModeInfo *this_mbmi = &miPtr[row * xd->mi_stride + col].mbmi;

                    // MbModeInfo *this_mbmi = xd->mi[row * xd->mi_stride + col];
                     //is_compound = has_second_ref(this_mbmi);
                    int32_t tmp_dst_stride = 8;
                    UNUSED_VARIABLE(tmp_dst_stride);
                    assert(bw < 8 || bh < 8);

                    // ConvolveParams conv_params = get_conv_params_no_round(
                     //    0, plane, xd->tmp_conv_dst, tmp_dst_stride, is_compound, xd->bd);
                    conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, BLOCK_SIZE_64, is_compound, EB_8BIT);
                    conv_params.use_jnt_comp_avg = 0;

                    EbPictureBufferDesc_t                  *ref_pic = this_mbmi->ref_frame[0] == LAST_FRAME ? ref_pic_list0 : ref_pic_list1;
                    ASSERT(ref_pic != NULL);
                    src_ptr = ref_pic->bufferCb + (ref_pic->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic->strideCb;
                    dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
                    src_stride = ref_pic->strideCb;
                    dst_stride = prediction_ptr->strideCb;
                    src_ptr = src_ptr + x + y * ref_pic->strideCb;
                    dst_ptr = dst_ptr + x + y * prediction_ptr->strideCb;

                    const MV mv = this_mbmi->mv[0].as_mv;
                    mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
                    subpel_x = mv_q4.col & SUBPEL_MASK;
                    subpel_y = mv_q4.row & SUBPEL_MASK;
                    src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);

                    av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                        &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);


                    convolve[subpel_x != 0][subpel_y != 0][is_compound](
                        src_ptr,
                        src_stride,
                        dst_ptr,
                        dst_stride,
                        b4_w,
                        b4_h,
                        &filter_params_x,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                        &filter_params_y,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                        subpel_x,
                        subpel_y,
                        &conv_params);

                    //Cr
                    conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, BLOCK_SIZE_64, is_compound, EB_8BIT);
                    conv_params.use_jnt_comp_avg = 0;

                    src_ptr = ref_pic->bufferCr + (ref_pic->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic->strideCr;
                    dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;

                    src_stride = ref_pic->strideCr;
                    dst_stride = prediction_ptr->strideCr;
                    src_ptr = src_ptr + x + y * ref_pic->strideCr;
                    dst_ptr = dst_ptr + x + y * prediction_ptr->strideCr;

                    // const MV mv = this_mbmi->mv[0].as_mv;
                    mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
                    subpel_x = mv_q4.col & SUBPEL_MASK;
                    subpel_y = mv_q4.row & SUBPEL_MASK;
                    src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);

                    av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                        &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

                    convolve[subpel_x != 0][subpel_y != 0][is_compound](
                        src_ptr,
                        src_stride,
                        dst_ptr,
                        dst_stride,
                        b4_w,
                        b4_h,
                        &filter_params_x,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                        &filter_params_y,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                        subpel_x,
                        subpel_y,
                        &conv_params);


                    ++col;
                }
                ++row;
            }

            //for (ref = 0; ref < 2; ++ref) pd->pre[ref] = orig_pred_buf[ref];

            //return;
        }


    }


    if (mv_unit->predDirection == UNI_PRED_LIST_0 || mv_unit->predDirection == BI_PRED) {

        //List0-Y
        mv.col = mv_unit->mv[REF_LIST_0].x;
        mv.row = mv_unit->mv[REF_LIST_0].y;

        src_ptr = ref_pic_list0->buffer_y + ref_pic_list0->origin_x + pu_origin_x + (ref_pic_list0->origin_y + pu_origin_y) * ref_pic_list0->stride_y;
        dst_ptr = prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
        src_stride = ref_pic_list0->stride_y;
        dst_stride = prediction_ptr->stride_y;

        mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, bwidth, bheight, 0, 0);//mv_q4 has 1 extra bit for fractionnal to accomodate chroma when accessing filter coeffs.
        subpel_x = mv_q4.col & SUBPEL_MASK;
        subpel_y = mv_q4.row & SUBPEL_MASK;
        src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
        conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstY, 128, is_compound, EB_8BIT);
        av1_get_convolve_filter_params(interp_filters, &filter_params_x,
            &filter_params_y, bwidth, bheight);

        convolve[subpel_x != 0][subpel_y != 0][is_compound](
            src_ptr,
            src_stride,
            dst_ptr,
            dst_stride,
            bwidth,
            bheight,
            &filter_params_x,//av1RegularFilter,
            &filter_params_y,//av1RegularFilter,
            subpel_x,
            subpel_y,
            &conv_params);
#if CHROMA_BLIND
        if (perform_chroma && blk_geom->has_uv && sub8x8_inter == 0) {
#else
        if (blk_geom->has_uv && sub8x8_inter == 0) {
#endif
            //List0-Cb
            src_ptr = ref_pic_list0->bufferCb + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCb;
            dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            src_stride = ref_pic_list0->strideCb;
            dst_stride = prediction_ptr->strideCb;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, 64, is_compound, EB_8BIT);


            av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                &filter_params_x,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                &filter_params_y,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,

                subpel_x,
                subpel_y,
                &conv_params);

            //List0-Cr
            src_ptr = ref_pic_list0->bufferCr + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCr;
            dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            src_stride = ref_pic_list0->strideCr;
            dst_stride = prediction_ptr->strideCr;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, 64, is_compound, EB_8BIT);
            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                &filter_params_x,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                &filter_params_y,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                subpel_x,
                subpel_y,
                &conv_params);
        }
    }

    if ((mv_unit->predDirection == UNI_PRED_LIST_1 || mv_unit->predDirection == BI_PRED) ) {

        //List0-Y
        mv.col = mv_unit->mv[REF_LIST_1].x;
        mv.row = mv_unit->mv[REF_LIST_1].y;
        ASSERT(ref_pic_list1 != NULL);
        src_ptr = ref_pic_list1->buffer_y + ref_pic_list1->origin_x + pu_origin_x + (ref_pic_list1->origin_y + pu_origin_y) * ref_pic_list1->stride_y;
        dst_ptr = prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
        src_stride = ref_pic_list1->stride_y;
        dst_stride = prediction_ptr->stride_y;

        mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, bwidth, bheight, 0, 0);//mv_q4 has 1 extra bit for fractionnal to accomodate chroma when accessing filter coeffs.
        subpel_x = mv_q4.col & SUBPEL_MASK;
        subpel_y = mv_q4.row & SUBPEL_MASK;

        src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
        conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstY, 128, is_compound, EB_8BIT);
        av1_get_convolve_filter_params(interp_filters, &filter_params_x,
            &filter_params_y, bwidth, bheight);

        convolve[subpel_x != 0][subpel_y != 0][is_compound](
            src_ptr,
            src_stride,
            dst_ptr,
            dst_stride,
            bwidth,
            bheight,
            &filter_params_x,//&av1RegularFilter,
            &filter_params_y,//&av1RegularFilter,

            subpel_x,
            subpel_y,
            &conv_params);
#if CHROMA_BLIND
        if (perform_chroma && blk_geom->has_uv && sub8x8_inter == 0) {
#else
        if (blk_geom->has_uv && sub8x8_inter == 0) {
#endif
            //List0-Cb
            src_ptr = ref_pic_list1->bufferCb + (ref_pic_list1->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list1->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list1->strideCb;
            dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            src_stride = ref_pic_list1->strideCb;
            dst_stride = prediction_ptr->strideCb;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstCb, 64, is_compound, EB_8BIT);

            av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);


            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                &filter_params_x,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                &filter_params_y,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                subpel_x,
                subpel_y,
                &conv_params);

            //List0-Cr
            src_ptr = ref_pic_list1->bufferCr + (ref_pic_list1->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list1->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list1->strideCr;
            dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            src_stride = ref_pic_list1->strideCr;
            dst_stride = prediction_ptr->strideCr;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstCr, 64, is_compound, EB_8BIT);
            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                &filter_params_x,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                &filter_params_y,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                subpel_x,
                subpel_y,
                &conv_params);
        }


    }


    return return_error;
}


/***************************************************
*  PreLoad Reference Block  for 16bit mode
***************************************************/
void Av1UnPackReferenceBlock(
    uint16_t                *SrcBuffer,
    uint32_t                 SrcStride,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    uint8_t                  *eightBitBuffer,
    uint32_t                 eightBitStride,
    EbBool                sub_pred,
    EbAsm                 asm_type,
    uint8_t                  Tap)
{

    pu_width += Tap;
    pu_height += Tap;
    uint16_t *ptr16 = (uint16_t *)SrcBuffer - (Tap >> 1) - ((Tap >> 1)*SrcStride);

    extract8_bitdata_safe_sub(
        ptr16,
        SrcStride << sub_pred,
        eightBitBuffer,
        eightBitStride << sub_pred,
        pu_width,
        pu_height >> sub_pred,
        sub_pred,
        asm_type
    );
}


#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
EbErrorType AV1InterPrediction10BitMD(
    uint32_t                                interp_filters,
#else
EbErrorType AV1MDInterPrediction(
#endif
    PictureControlSet_t                     *picture_control_set_ptr,
    uint8_t                                  ref_frame_type,
    ModeDecisionContext_t                   *md_context_ptr,
    CodingUnit_t                            *cu_ptr,
    MvUnit_t                                *mv_unit,
    uint16_t                                 pu_origin_x,
    uint16_t                                 pu_origin_y,
    uint8_t                                  bwidth,
    uint8_t                                  bheight,
    EbPictureBufferDesc_t                   *ref_pic_list0,
    EbPictureBufferDesc_t                   *ref_pic_list1,
    EbPictureBufferDesc_t                   *prediction_ptr,
    uint16_t                                 dst_origin_x,
    uint16_t                                 dst_origin_y,
#if CHROMA_BLIND
    EbBool                                  perform_chroma,
#endif
    EbAsm                                    asm_type)
{
    EbErrorType  return_error = EB_ErrorNone;
    InterPredictionContext_t *context_ptr = (InterPredictionContext_t*)(md_context_ptr->inter_prediction_context);
    uint8_t         is_compound = (mv_unit->predDirection == BI_PRED) ? 1 : 0;
    DECLARE_ALIGNED(32, uint16_t, tmp_dstY[128 * 128]);//move this to context if stack does not hold.
    DECLARE_ALIGNED(32, uint16_t, tmp_dstCb[64 * 64]);
    DECLARE_ALIGNED(32, uint16_t, tmp_dstCr[64 * 64]);



    MV  mv, mv_q4;

    int32_t subpel_x, subpel_y;
    uint16_t * src_ptr;
    uint8_t * dst_ptr;
    int32_t src_stride;
    int32_t dst_stride;
    ConvolveParams conv_params;
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
    InterpFilterParams filter_params_x, filter_params_y;
#endif

    const BlockGeom * blk_geom = get_blk_geom_mds(cu_ptr->mds_idx);

    //special treatment for chroma in 4XN/NX4 blocks
   //if one of the neighbour blocks of the parent square is intra the chroma prediction will follow the normal path using the luma MV of the current nsq block which is the latest sub8x8.
   //for this case: only uniPred is allowed.

    int32_t sub8x8_inter = 0;
#if CHROMA_BLIND
    if (perform_chroma && blk_geom->has_uv && (blk_geom->bwidth == 4 || blk_geom->bheight == 4))
#else
    if (blk_geom->has_uv &&
        (blk_geom->bwidth == 4 || blk_geom->bheight == 4)
        )
#endif
    {

        //CHKN setup input param

        int32_t bw = blk_geom->bwidth_uv;
        int32_t bh = blk_geom->bheight_uv;
        UNUSED_VARIABLE(bw);
        UNUSED_VARIABLE(bh);

        uint32_t mi_x = pu_origin_x;       //these are luma picture wise
        uint32_t mi_y = pu_origin_y;

        MacroBlockD  *xd = cu_ptr->av1xd;
        xd->mi_stride = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->picture_width_in_sb*(BLOCK_SIZE_64 / 4);
        const int32_t offset = (mi_y >> MI_SIZE_LOG2) * xd->mi_stride + (mi_x >> MI_SIZE_LOG2);
        xd->mi = picture_control_set_ptr->mi_grid_base + offset;

        //CHKN fill current mi from current block
        {
            ModeInfo *miPtr = *xd->mi;
            uint8_t  miX, miY;
            MvReferenceFrame rf[2];
            av1_set_ref_frame(rf, ref_frame_type);
            for (miY = 0; miY < (blk_geom->bheight >> MI_SIZE_LOG2); miY++) {
                for (miX = 0; miX < (blk_geom->bwidth >> MI_SIZE_LOG2); miX++) {
                    miPtr[miX + miY * xd->mi_stride].mbmi.ref_frame[0] = rf[0];
                    if (mv_unit->predDirection == UNI_PRED_LIST_0) {
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_0].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_0].y;
                    }
                    else if (mv_unit->predDirection == UNI_PRED_LIST_1) {
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_1].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_1].y;
                    }
                    else {

                        // printf("ERRRRRRR");

                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_0].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_0].y;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[1].as_mv.col = mv_unit->mv[REF_LIST_1].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[1].as_mv.row = mv_unit->mv[REF_LIST_1].y;
                    }

                }
            }
        }

        int32_t build_for_obmc = 0;

        const block_size bsize = blk_geom->bsize;//mi->sb_type;
        ASSERT(bsize < BlockSizeS_ALL);
        const int32_t ss_x = 1;// pd->subsampling_x;
        const int32_t ss_y = 1;//pd->subsampling_y;
        sub8x8_inter = (block_size_wide[bsize] < 8 && ss_x) ||
            (block_size_high[bsize] < 8 && ss_y);

        //if (is_intrabc) sub8x8_inter = 0;

        // For sub8x8 chroma blocks, we may be covering more than one luma block's
        // worth of pixels. Thus (mi_x, mi_y) may not be the correct coordinates for
        // the top-left corner of the prediction source - the correct top-left corner
        // is at (pre_x, pre_y).
        const int32_t row_start =
            (block_size_high[bsize] == 4) && ss_y && !build_for_obmc ? -1 : 0;
        const int32_t col_start =
            (block_size_wide[bsize] == 4) && ss_x && !build_for_obmc ? -1 : 0;

        const int32_t pre_x = (mi_x + MI_SIZE * col_start) >> ss_x;
        const int32_t pre_y = (mi_y + MI_SIZE * row_start) >> ss_y;
        UNUSED_VARIABLE(pre_x);
        UNUSED_VARIABLE(pre_y);

        sub8x8_inter = sub8x8_inter && !build_for_obmc;
        if (sub8x8_inter) {
            for (int32_t row = row_start; row <= 0 && sub8x8_inter; ++row) {
                for (int32_t col = col_start; col <= 0; ++col) {

                    ModeInfo *miPtr = *xd->mi;
                    const MbModeInfo *this_mbmi = &miPtr[row * xd->mi_stride + col].mbmi;

                    if (!is_inter_block(this_mbmi)) sub8x8_inter = 0;
                    //if (is_intrabc_block(this_mbmi)) sub8x8_inter = 0;
                }
            }
        }

        if (sub8x8_inter) {
            // block size
            const int32_t b4_w = block_size_wide[bsize] >> ss_x;
            const int32_t b4_h = block_size_high[bsize] >> ss_y;
            const block_size plane_bsize = scale_chroma_bsize(bsize, ss_x, ss_y);
            ASSERT(plane_bsize < BlockSizeS_ALL);
            const int32_t b8_w = block_size_wide[plane_bsize] >> ss_x;
            const int32_t b8_h = block_size_high[plane_bsize] >> ss_y;

            assert(!is_compound);
            if (is_compound)
                printf("ETTTT");

            //const struct Buf2d orig_pred_buf[2] = { pd->pre[0], pd->pre[1] };

            int32_t row = row_start;
            int32_t src_stride;
            for (int32_t y = 0; y < b8_h; y += b4_h) {
                int32_t col = col_start;
                for (int32_t x = 0; x < b8_w; x += b4_w) {

                    ModeInfo *miPtr = *xd->mi;
                    const MbModeInfo *this_mbmi = &miPtr[row * xd->mi_stride + col].mbmi;

                    // MbModeInfo *this_mbmi = xd->mi[row * xd->mi_stride + col];
                     //is_compound = has_second_ref(this_mbmi);
                    int32_t tmp_dst_stride = 8;
                    UNUSED_VARIABLE(tmp_dst_stride);
                    assert(bw < 8 || bh < 8);

                    // ConvolveParams conv_params = get_conv_params_no_round(
                     //    0, plane, xd->tmp_conv_dst, tmp_dst_stride, is_compound, xd->bd);
                    conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, BLOCK_SIZE_64, is_compound, EB_8BIT);
                    conv_params.use_jnt_comp_avg = 0;

                    EbPictureBufferDesc_t                  *ref_pic = this_mbmi->ref_frame[0] == LAST_FRAME ? ref_pic_list0 : ref_pic_list1;
                    ASSERT(ref_pic != NULL);
                    src_ptr = (uint16_t*)ref_pic->bufferCb + (ref_pic->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic->strideCb;
                    dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;

                    src_stride = ref_pic->strideCb;
                    dst_stride = prediction_ptr->strideCb;
                    src_ptr = src_ptr + x + y * ref_pic->strideCb;
                    dst_ptr = dst_ptr + x + y * prediction_ptr->strideCb;

                    const MV mv = this_mbmi->mv[0].as_mv;
                    mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
                    subpel_x = mv_q4.col & SUBPEL_MASK;
                    subpel_y = mv_q4.row & SUBPEL_MASK;
                    src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);

                    EbPictureBufferDesc_t                  *local_buffer = this_mbmi->ref_frame[0] == LAST_FRAME ? context_ptr->mcp_context->local_reference_block8_bitl0 : context_ptr->mcp_context->local_reference_block8_bitl1;

                    av1_get_convolve_filter_params(cu_ptr->interp_filters, &filter_params_x,
                        &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

                    Av1UnPackReferenceBlock(
                        src_ptr,
                        src_stride,
                        b4_w,
                        b4_h,
                        local_buffer->bufferCb,
                        local_buffer->strideCb,
                        EB_FALSE,
                        asm_type,
                        8);

                    convolve[subpel_x != 0][subpel_y != 0][is_compound](
                        local_buffer->bufferCb + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl1->strideCb),
                        local_buffer->strideCb,
                        dst_ptr,
                        dst_stride,
                        b4_w,
                        b4_h,
                        &filter_params_x,
                        &filter_params_y,
                        subpel_x,
                        subpel_y,
                        &conv_params);


                    //Cr
                    conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, BLOCK_SIZE_64, is_compound, EB_8BIT);
                    conv_params.use_jnt_comp_avg = 0;

                    src_ptr = (uint16_t*)ref_pic->bufferCr + (ref_pic->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic->strideCr;
                    dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
                    src_stride = ref_pic->strideCr;
                    dst_stride = prediction_ptr->strideCr;
                    src_ptr = src_ptr + x + y * ref_pic->strideCr;
                    dst_ptr = dst_ptr + x + y * prediction_ptr->strideCr;

                    // const MV mv = this_mbmi->mv[0].as_mv;
                    mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
                    subpel_x = mv_q4.col & SUBPEL_MASK;
                    subpel_y = mv_q4.row & SUBPEL_MASK;
                    src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);

                    av1_get_convolve_filter_params(cu_ptr->interp_filters, &filter_params_x,
                        &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

                    //EbPictureBufferDesc_t                  *local_buffer = this_mbmi->ref_frame[0] == LAST_FRAME ? context_ptr->mcp_context->local_reference_block8_bitl0 : context_ptr->mcp_context->local_reference_block8_bitl1;

                    Av1UnPackReferenceBlock(
                        src_ptr,
                        src_stride,
                        b4_w,
                        b4_h,
                        local_buffer->bufferCr,
                        local_buffer->strideCr,
                        EB_FALSE,
                        asm_type,
                        8);

                    convolve[subpel_x != 0][subpel_y != 0][is_compound](
                        local_buffer->bufferCr + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl1->strideCr),
                        local_buffer->strideCr,
                        dst_ptr,
                        dst_stride,
                        b4_w,
                        b4_h,
                        &filter_params_x,
                        &filter_params_y,
                        subpel_x,
                        subpel_y,
                        &conv_params);

                    ++col;
                }
                ++row;
            }

            //for (ref = 0; ref < 2; ++ref) pd->pre[ref] = orig_pred_buf[ref];

            //return;
        }


    }



    if (mv_unit->predDirection == UNI_PRED_LIST_0 || mv_unit->predDirection == BI_PRED) {

        //List0-Y
        mv.col = mv_unit->mv[REF_LIST_0].x;
        mv.row = mv_unit->mv[REF_LIST_0].y;

        src_ptr = (uint16_t*)ref_pic_list0->buffer_y + ref_pic_list0->origin_x + pu_origin_x + (ref_pic_list0->origin_y + pu_origin_y) * ref_pic_list0->stride_y;
        dst_ptr = prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
        src_stride = ref_pic_list0->stride_y;
        dst_stride = prediction_ptr->stride_y;
        mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, bwidth, bheight, 0, 0);//mv_q4 has 1 extra bit for fractionnal to accomodate chroma when accessing filter coeffs.

        subpel_x = mv_q4.col & SUBPEL_MASK;
        subpel_y = mv_q4.row & SUBPEL_MASK;
        src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
        conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstY, 128, is_compound, EB_8BIT);

#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
        av1_get_convolve_filter_params(interp_filters, &filter_params_x,
            &filter_params_y, bwidth, bheight);
#endif

        Av1UnPackReferenceBlock(
            src_ptr,
            ref_pic_list0->stride_y,
            bwidth,
            bheight,
            context_ptr->mcp_context->local_reference_block8_bitl0->buffer_y,
            context_ptr->mcp_context->local_reference_block8_bitl0->stride_y,
            EB_FALSE,
            asm_type,
            16);

        convolve[subpel_x != 0][subpel_y != 0][is_compound](
            context_ptr->mcp_context->local_reference_block8_bitl0->buffer_y + 8 + (8 * context_ptr->mcp_context->local_reference_block8_bitl0->stride_y),
            context_ptr->mcp_context->local_reference_block8_bitl0->stride_y,
            dst_ptr,
            dst_stride,
            bwidth,
            bheight,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
            &filter_params_x,
            &filter_params_y,
#else
            &av1RegularFilter,
            &av1RegularFilter,
#endif
            subpel_x,
            subpel_y,
            &conv_params);

#if CHROMA_BLIND
        if (perform_chroma && blk_geom->has_uv && sub8x8_inter == 0) {
#else
        if (blk_geom->has_uv && sub8x8_inter == 0) {
#endif
            //List0-Cb
            src_ptr = (uint16_t*)ref_pic_list0->bufferCb + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCb;
            dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            src_stride = ref_pic_list0->strideCb;
            dst_stride = prediction_ptr->strideCb;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);

            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, 64, is_compound, EB_8BIT);
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
            av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

#endif

            Av1UnPackReferenceBlock(
                src_ptr,
                ref_pic_list0->strideCb,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                context_ptr->mcp_context->local_reference_block8_bitl0->bufferCb,
                context_ptr->mcp_context->local_reference_block8_bitl0->strideCb,
                EB_FALSE,
                asm_type,
                8);

            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                context_ptr->mcp_context->local_reference_block8_bitl0->bufferCb + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl0->strideCb),
                context_ptr->mcp_context->local_reference_block8_bitl0->strideCb,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
                &filter_params_x,
                &filter_params_y,
#else
                puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
#endif
                subpel_x,
                subpel_y,
                &conv_params);

            //List0-Cr
            src_ptr = (uint16_t*)ref_pic_list0->bufferCr + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCr;
            dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            src_stride = ref_pic_list0->strideCr;
            dst_stride = prediction_ptr->strideCr;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);

            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, 64, is_compound, EB_8BIT);

            Av1UnPackReferenceBlock(
                src_ptr,
                ref_pic_list0->strideCr,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                context_ptr->mcp_context->local_reference_block8_bitl0->bufferCr,
                context_ptr->mcp_context->local_reference_block8_bitl0->strideCr,
                EB_FALSE,
                asm_type,
                8);

            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                context_ptr->mcp_context->local_reference_block8_bitl0->bufferCr + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl1->strideCr),
                context_ptr->mcp_context->local_reference_block8_bitl0->strideCr,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
                &filter_params_x,
                &filter_params_y,
#else
                puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
#endif
                subpel_x,
                subpel_y,
                &conv_params);
        }

    }


    if (mv_unit->predDirection == UNI_PRED_LIST_1 || mv_unit->predDirection == BI_PRED) {

        //List1-Y
        mv.col = mv_unit->mv[REF_LIST_1].x;
        mv.row = mv_unit->mv[REF_LIST_1].y;
        ASSERT(ref_pic_list1 != NULL);
        src_ptr = (uint16_t*)ref_pic_list1->buffer_y + ref_pic_list1->origin_x + pu_origin_x + (ref_pic_list1->origin_y + pu_origin_y) * ref_pic_list1->stride_y;
        dst_ptr = prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
        src_stride = ref_pic_list1->stride_y;
        dst_stride = prediction_ptr->stride_y;

        mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, bwidth, bheight, 0, 0);//mv_q4 has 1 extra bit for fractionnal to accomodate chroma when accessing filter coeffs.

        subpel_x = mv_q4.col & SUBPEL_MASK;
        subpel_y = mv_q4.row & SUBPEL_MASK;

        src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
        conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstY, 128, is_compound, EB_8BIT);

#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
        av1_get_convolve_filter_params(interp_filters, &filter_params_x,
            &filter_params_y, bwidth, bheight);

#endif

        Av1UnPackReferenceBlock(
            src_ptr,
            ref_pic_list1->stride_y,
            bwidth,
            bheight,
            context_ptr->mcp_context->local_reference_block8_bitl1->buffer_y,
            context_ptr->mcp_context->local_reference_block8_bitl1->stride_y,
            EB_FALSE,
            asm_type,
            16);

        convolve[subpel_x != 0][subpel_y != 0][is_compound](
            context_ptr->mcp_context->local_reference_block8_bitl1->buffer_y + 8 + (8 * context_ptr->mcp_context->local_reference_block8_bitl1->stride_y),
            context_ptr->mcp_context->local_reference_block8_bitl1->stride_y,
            dst_ptr,
            dst_stride,
            bwidth,
            bheight,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
            &filter_params_x,
            &filter_params_y,
#else
            &av1RegularFilter,
            &av1RegularFilter,
#endif
            subpel_x,
            subpel_y,
            &conv_params);

#if CHROMA_BLIND
        if (perform_chroma && blk_geom->has_uv && sub8x8_inter == 0) {
#else
        if (blk_geom->has_uv && sub8x8_inter == 0) {
#endif
            //List1-Cb
            src_ptr = (uint16_t*)ref_pic_list1->bufferCb + (ref_pic_list1->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list1->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list1->strideCb;
            dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            src_stride = ref_pic_list1->strideCb;
            dst_stride = prediction_ptr->strideCb;
            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstCb, 64, is_compound, EB_8BIT);


#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
            av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);
#endif
            Av1UnPackReferenceBlock(
                src_ptr,
                ref_pic_list1->strideCb,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                context_ptr->mcp_context->local_reference_block8_bitl1->bufferCb,
                context_ptr->mcp_context->local_reference_block8_bitl1->strideCb,
                EB_FALSE,
                asm_type,
                8);

            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                context_ptr->mcp_context->local_reference_block8_bitl1->bufferCb + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl1->strideCb),
                context_ptr->mcp_context->local_reference_block8_bitl1->strideCb,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
                &filter_params_x,
                &filter_params_y,
#else
                puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
#endif
                subpel_x,
                subpel_y,
                &conv_params);

            //List1-Cr
            src_ptr = (uint16_t*)ref_pic_list1->bufferCr + (ref_pic_list1->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list1->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list1->strideCr;
            dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            src_stride = ref_pic_list1->strideCr;
            dst_stride = prediction_ptr->strideCr;
            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);

            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstCr, 64, is_compound, EB_8BIT);
            Av1UnPackReferenceBlock(
                src_ptr,
                ref_pic_list1->strideCr,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                context_ptr->mcp_context->local_reference_block8_bitl1->bufferCr,
                context_ptr->mcp_context->local_reference_block8_bitl1->strideCr,
                EB_FALSE,
                asm_type,
                8);

            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                context_ptr->mcp_context->local_reference_block8_bitl1->bufferCr + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl1->strideCr),
                context_ptr->mcp_context->local_reference_block8_bitl1->strideCr,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
                &filter_params_x,
                &filter_params_y,
#else
                puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
#endif
                subpel_x,
                subpel_y,
                &conv_params);
        }

    }


    return return_error;
}


EbErrorType av1_inter_prediction_hbd(
    PictureControlSet_t                    *picture_control_set_ptr,
    uint8_t                                   ref_frame_type,
    CodingUnit_t                           *cu_ptr,
    MvUnit_t                               *mv_unit,
    uint16_t                                  pu_origin_x,
    uint16_t                                  pu_origin_y,
    uint8_t                                   bwidth,
    uint8_t                                   bheight,
    EbPictureBufferDesc_t                  *ref_pic_list0,
    EbPictureBufferDesc_t                  *ref_pic_list1,
    EbPictureBufferDesc_t                  *prediction_ptr,
    uint16_t                                  dst_origin_x,
    uint16_t                                  dst_origin_y,
    uint8_t                                   bit_depth,
    EbAsm                                  asm_type)
{
    (void)asm_type;
    EbErrorType  return_error = EB_ErrorNone;
    uint8_t         is_compound = (mv_unit->predDirection == BI_PRED) ? 1 : 0;
    DECLARE_ALIGNED(32, uint16_t, tmp_dstY[128 * 128]);//move this to context if stack does not hold.
    DECLARE_ALIGNED(32, uint16_t, tmp_dstCb[64 * 64]);
    DECLARE_ALIGNED(32, uint16_t, tmp_dstCr[64 * 64]);
    MV  mv, mv_q4;

    int32_t subpel_x, subpel_y;
    uint16_t * src_ptr;
    uint16_t * dst_ptr;
    int32_t src_stride;
    int32_t dst_stride;
    ConvolveParams conv_params;
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
    InterpFilterParams filter_params_x, filter_params_y;
#endif

    const BlockGeom * blk_geom = get_blk_geom_mds(cu_ptr->mds_idx);

    //special treatment for chroma in 4XN/NX4 blocks
   //if one of the neighbour blocks of the parent square is intra the chroma prediction will follow the normal path using the luma MV of the current nsq block which is the latest sub8x8.
   //for this case: only uniPred is allowed.

    int32_t sub8x8_inter = 0;

    if (blk_geom->has_uv &&
        (blk_geom->bwidth == 4 || blk_geom->bheight == 4)
        )
    {

        //CHKN setup input param

        int32_t bw = blk_geom->bwidth_uv;
        int32_t bh = blk_geom->bheight_uv;
        UNUSED_VARIABLE(bw);
        UNUSED_VARIABLE(bh);

        uint32_t mi_x = pu_origin_x;       //these are luma picture wise
        uint32_t mi_y = pu_origin_y;

        MacroBlockD  *xd = cu_ptr->av1xd;
        xd->mi_stride = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->picture_width_in_sb*(BLOCK_SIZE_64 / 4);
        const int32_t offset = (mi_y >> MI_SIZE_LOG2) * xd->mi_stride + (mi_x >> MI_SIZE_LOG2);
        xd->mi = picture_control_set_ptr->mi_grid_base + offset;

        //CHKN fill current mi from current block
        {
            ModeInfo *miPtr = *xd->mi;
            uint8_t  miX, miY;
            MvReferenceFrame rf[2];
            av1_set_ref_frame(rf, ref_frame_type);
            for (miY = 0; miY < (blk_geom->bheight >> MI_SIZE_LOG2); miY++) {
                for (miX = 0; miX < (blk_geom->bwidth >> MI_SIZE_LOG2); miX++) {
                    miPtr[miX + miY * xd->mi_stride].mbmi.ref_frame[0] = rf[0];
                    if (mv_unit->predDirection == UNI_PRED_LIST_0) {
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_0].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_0].y;
                    }
                    else if (mv_unit->predDirection == UNI_PRED_LIST_1) {
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_1].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_1].y;
                    }
                    else {

                        // printf("ERRRRRRR");

                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.col = mv_unit->mv[REF_LIST_0].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[0].as_mv.row = mv_unit->mv[REF_LIST_0].y;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[1].as_mv.col = mv_unit->mv[REF_LIST_1].x;
                        miPtr[miX + miY * xd->mi_stride].mbmi.mv[1].as_mv.row = mv_unit->mv[REF_LIST_1].y;
                    }

                }
            }
        }

        int32_t build_for_obmc = 0;

        const block_size bsize = blk_geom->bsize;//mi->sb_type;
        ASSERT(bsize < BlockSizeS_ALL);
        const int32_t ss_x = 1;// pd->subsampling_x;
        const int32_t ss_y = 1;//pd->subsampling_y;
        sub8x8_inter = (block_size_wide[bsize] < 8 && ss_x) ||
            (block_size_high[bsize] < 8 && ss_y);

        //if (is_intrabc) sub8x8_inter = 0;

        // For sub8x8 chroma blocks, we may be covering more than one luma block's
        // worth of pixels. Thus (mi_x, mi_y) may not be the correct coordinates for
        // the top-left corner of the prediction source - the correct top-left corner
        // is at (pre_x, pre_y).
        const int32_t row_start =
            (block_size_high[bsize] == 4) && ss_y && !build_for_obmc ? -1 : 0;
        const int32_t col_start =
            (block_size_wide[bsize] == 4) && ss_x && !build_for_obmc ? -1 : 0;

        const int32_t pre_x = (mi_x + MI_SIZE * col_start) >> ss_x;
        const int32_t pre_y = (mi_y + MI_SIZE * row_start) >> ss_y;
        UNUSED_VARIABLE(pre_x);
        UNUSED_VARIABLE(pre_y);

        sub8x8_inter = sub8x8_inter && !build_for_obmc;
        if (sub8x8_inter) {
            for (int32_t row = row_start; row <= 0 && sub8x8_inter; ++row) {
                for (int32_t col = col_start; col <= 0; ++col) {

                    ModeInfo *miPtr = *xd->mi;
                    const MbModeInfo *this_mbmi = &miPtr[row * xd->mi_stride + col].mbmi;

                    if (!is_inter_block(this_mbmi)) sub8x8_inter = 0;
                    //if (is_intrabc_block(this_mbmi)) sub8x8_inter = 0;
                }
            }
        }

        if (sub8x8_inter) {
            // block size
            const int32_t b4_w = block_size_wide[bsize] >> ss_x;
            const int32_t b4_h = block_size_high[bsize] >> ss_y;
            const block_size plane_bsize = scale_chroma_bsize(bsize, ss_x, ss_y);
            ASSERT(plane_bsize < BlockSizeS_ALL);
            const int32_t b8_w = block_size_wide[plane_bsize] >> ss_x;
            const int32_t b8_h = block_size_high[plane_bsize] >> ss_y;

            assert(!is_compound);

            if (is_compound)
                printf("ETTTT");

            //const struct Buf2d orig_pred_buf[2] = { pd->pre[0], pd->pre[1] };

            int32_t row = row_start;
            int32_t src_stride;
            for (int32_t y = 0; y < b8_h; y += b4_h) {
                int32_t col = col_start;
                for (int32_t x = 0; x < b8_w; x += b4_w) {

                    ModeInfo *miPtr = *xd->mi;
                    const MbModeInfo *this_mbmi = &miPtr[row * xd->mi_stride + col].mbmi;

                    // MbModeInfo *this_mbmi = xd->mi[row * xd->mi_stride + col];
                     //is_compound = has_second_ref(this_mbmi);
                    int32_t tmp_dst_stride = 8;
                    UNUSED_VARIABLE(tmp_dst_stride);
                    assert(bw < 8 || bh < 8);

                    // ConvolveParams conv_params = get_conv_params_no_round(
                     //    0, plane, xd->tmp_conv_dst, tmp_dst_stride, is_compound, xd->bd);
                    conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, BLOCK_SIZE_64, is_compound, bit_depth);
                    conv_params.use_jnt_comp_avg = 0;
                    EbPictureBufferDesc_t                  *ref_pic = this_mbmi->ref_frame[0] == LAST_FRAME ? ref_pic_list0 : ref_pic_list1;
                    src_ptr = (uint16_t*)ref_pic->bufferCb + (ref_pic->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic->strideCb;
                    dst_ptr = (uint16_t*)prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
                    src_stride = ref_pic->strideCb;
                    dst_stride = prediction_ptr->strideCb;
                    src_ptr = src_ptr + x + y * ref_pic->strideCb;
                    dst_ptr = dst_ptr + x + y * prediction_ptr->strideCb;

                    const MV mv = this_mbmi->mv[0].as_mv;
                    mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
                    subpel_x = mv_q4.col & SUBPEL_MASK;
                    subpel_y = mv_q4.row & SUBPEL_MASK;
                    src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);

                    av1_get_convolve_filter_params(cu_ptr->interp_filters, &filter_params_x,
                        &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);


                    convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
                        src_ptr,
                        src_stride,
                        dst_ptr,
                        dst_stride,
                        b4_w,
                        b4_h,
                        &filter_params_x,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                        &filter_params_y,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                        subpel_x,
                        subpel_y,
                        &conv_params,
                        bit_depth);

                    //Cr
                    conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, BLOCK_SIZE_64, is_compound, bit_depth);
                    conv_params.use_jnt_comp_avg = 0;

                    src_ptr = (uint16_t*)ref_pic->bufferCr + (ref_pic->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic->strideCr;
                    dst_ptr = (uint16_t*)prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
                    src_stride = ref_pic->strideCr;
                    dst_stride = prediction_ptr->strideCr;
                    src_ptr = src_ptr + x + y * ref_pic->strideCr;
                    dst_ptr = dst_ptr + x + y * prediction_ptr->strideCr;

                    // const MV mv = this_mbmi->mv[0].as_mv;
                    mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
                    subpel_x = mv_q4.col & SUBPEL_MASK;
                    subpel_y = mv_q4.row & SUBPEL_MASK;
                    src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);

                    av1_get_convolve_filter_params(cu_ptr->interp_filters, &filter_params_x,
                        &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

                    convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
                        src_ptr,
                        src_stride,
                        dst_ptr,
                        dst_stride,
                        b4_w,
                        b4_h,
                        &filter_params_x,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                        &filter_params_y,//puSize > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                        subpel_x,
                        subpel_y,
                        &conv_params,
                        bit_depth);

                    ++col;
                }
                ++row;
            }

            //for (ref = 0; ref < 2; ++ref) pd->pre[ref] = orig_pred_buf[ref];

            //return;
        }


    }



    if (mv_unit->predDirection == UNI_PRED_LIST_0 || mv_unit->predDirection == BI_PRED) {

        //List0-Y
        mv.col = mv_unit->mv[REF_LIST_0].x;
        mv.row = mv_unit->mv[REF_LIST_0].y;

        src_ptr = (uint16_t*)ref_pic_list0->buffer_y + ref_pic_list0->origin_x + pu_origin_x + (ref_pic_list0->origin_y + pu_origin_y) * ref_pic_list0->stride_y;
        dst_ptr = (uint16_t*)prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
        src_stride = ref_pic_list0->stride_y;
        dst_stride = prediction_ptr->stride_y;

        mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, bwidth, bheight, 0, 0);//mv_q4 has 1 extra bit for fractionnal to accomodate chroma when accessing filter coeffs.
        subpel_x = mv_q4.col & SUBPEL_MASK;
        subpel_y = mv_q4.row & SUBPEL_MASK;
        src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
        conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstY, 128, is_compound, bit_depth);
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
        av1_get_convolve_filter_params(cu_ptr->interp_filters, &filter_params_x,
            &filter_params_y, bwidth, bheight);
#endif

        convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
            src_ptr,
            src_stride,
            dst_ptr,
            dst_stride,
            bwidth,
            bheight,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
            &filter_params_x,
            &filter_params_y,
#else
            &av1RegularFilter,
            &av1RegularFilter,
#endif
            subpel_x,
            subpel_y,
            &conv_params,
            bit_depth);

        if (blk_geom->has_uv && sub8x8_inter == 0) {

            //List0-Cb
            src_ptr = (uint16_t*)ref_pic_list0->bufferCb + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCb;
            dst_ptr = (uint16_t*)prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            src_stride = ref_pic_list0->strideCb;
            dst_stride = prediction_ptr->strideCb;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, 64, is_compound, bit_depth);


            av1_get_convolve_filter_params(cu_ptr->interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

            convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
                &filter_params_x,
                &filter_params_y,
#else
                bwidth > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                bheight > 8 ? &av1RegularFilter : &av1RegularFilterW4,
#endif
                subpel_x,
                subpel_y,
                &conv_params,
                bit_depth);

            //List0-Cr
            src_ptr = (uint16_t*)ref_pic_list0->bufferCr + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCr;
            dst_ptr = (uint16_t*)prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            src_stride = ref_pic_list0->strideCr;
            dst_stride = prediction_ptr->strideCr;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, 64, is_compound, bit_depth);
            convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
                &filter_params_x,
                &filter_params_y,
#else
                bwidth > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                bheight > 8 ? &av1RegularFilter : &av1RegularFilterW4,
#endif
                subpel_x,
                subpel_y,
                &conv_params,
                bit_depth);
        }

    }


    if (mv_unit->predDirection == UNI_PRED_LIST_1 || mv_unit->predDirection == BI_PRED) {

        //List0-Y
        mv.col = mv_unit->mv[REF_LIST_1].x;
        mv.row = mv_unit->mv[REF_LIST_1].y;

        src_ptr = (uint16_t*)ref_pic_list1->buffer_y + ref_pic_list1->origin_x + pu_origin_x + (ref_pic_list1->origin_y + pu_origin_y) * ref_pic_list1->stride_y;
        dst_ptr = (uint16_t*)prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
        src_stride = ref_pic_list1->stride_y;
        dst_stride = prediction_ptr->stride_y;

        mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, bwidth, bheight, 0, 0);//mv_q4 has 1 extra bit for fractionnal to accomodate chroma when accessing filter coeffs.
        subpel_x = mv_q4.col & SUBPEL_MASK;
        subpel_y = mv_q4.row & SUBPEL_MASK;

        src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
        conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstY, 128, is_compound, bit_depth);
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
        av1_get_convolve_filter_params(cu_ptr->interp_filters, &filter_params_x,
            &filter_params_y, bwidth, bheight);
#endif

        convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
            src_ptr,
            src_stride,
            dst_ptr,
            dst_stride,
            bwidth,
            bheight,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
            &filter_params_x,
            &filter_params_y,
#else
            &av1RegularFilter,
            &av1RegularFilter,
#endif
            subpel_x,
            subpel_y,
            &conv_params,
            bit_depth);

        if (blk_geom->has_uv && sub8x8_inter == 0) {

            //List0-Cb
            src_ptr = (uint16_t*)ref_pic_list1->bufferCb + (ref_pic_list1->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list1->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list1->strideCb;
            dst_ptr = (uint16_t*)prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            src_stride = ref_pic_list1->strideCb;
            dst_stride = prediction_ptr->strideCb;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstCb, 64, is_compound, bit_depth);
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
            av1_get_convolve_filter_params(cu_ptr->interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);
#endif

            convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
                &filter_params_x,
                &filter_params_y,
#else
                bwidth > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                bheight > 8 ? &av1RegularFilter : &av1RegularFilterW4,
#endif
                subpel_x,
                subpel_y,
                &conv_params,
                bit_depth);

            //List0-Cr
            src_ptr = (uint16_t*)ref_pic_list1->bufferCr + (ref_pic_list1->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list1->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list1->strideCr;
            dst_ptr = (uint16_t*)prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            src_stride = ref_pic_list1->strideCr;
            dst_stride = prediction_ptr->strideCr;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, (mv_unit->predDirection == BI_PRED) ? 1 : 0, 0, tmp_dstCr, 64, is_compound, bit_depth);
            convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
                &filter_params_x,
                &filter_params_y,
#else
                bwidth > 8 ? &av1RegularFilter : &av1RegularFilterW4,
                bheight > 8 ? &av1RegularFilter : &av1RegularFilterW4,
#endif
                subpel_x,
                subpel_y,
                &conv_params,
                bit_depth);
        }


    }


    return return_error;
}


EbErrorType warped_motion_prediction(
    MvUnit_t                               *mv_unit,
    uint16_t                                pu_origin_x,
    uint16_t                                pu_origin_y,
    CodingUnit_t                           *cu_ptr,
    const BlockGeom                        *blk_geom,
    EbPictureBufferDesc_t                  *ref_pic_list0,
    EbPictureBufferDesc_t                  *prediction_ptr,
    uint16_t                                dst_origin_x,
    uint16_t                                dst_origin_y,
    EbWarpedMotionParams                   *wm_params,
    uint8_t                                 bit_depth,
#if CHROMA_BLIND
    EbBool                                  perform_chroma,
#endif
    EbAsm                                   asm_type)
{
    (void)asm_type;

    EbErrorType  return_error = EB_ErrorNone;
    uint8_t is_compound = (mv_unit->predDirection == BI_PRED) ? 1 : 0;
    assert(!is_compound);
    EbBool  is16bit = (EbBool)(bit_depth > EB_8BIT);

    int32_t src_stride;
    int32_t dst_stride;
    uint16_t buf_width;
    uint16_t buf_height;
    ConvolveParams conv_params;
    uint8_t ss_x = 1; // subsamplings
    uint8_t ss_y = 1;

    if (!is16bit) {
        uint8_t *src_ptr;
        uint8_t *dst_ptr;

        // Y - UNI_PRED_LIST_0
        src_ptr = ref_pic_list0->buffer_y + ref_pic_list0->origin_x + ref_pic_list0->origin_y * ref_pic_list0->stride_y;
        src_stride = ref_pic_list0->stride_y;
        buf_width = ref_pic_list0->width;
        buf_height = ref_pic_list0->height;

        dst_ptr = prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
        dst_stride = prediction_ptr->stride_y;
        conv_params = get_conv_params_no_round(0, 0, 0, NULL, 128, is_compound, EB_8BIT);

        av1_warp_plane(
            wm_params,
            (int) is16bit,
            bit_depth,
            src_ptr,
            (int) buf_width,
            (int) buf_height,
            src_stride,
            dst_ptr,
            pu_origin_x,
            pu_origin_y,
            blk_geom->bwidth,
            blk_geom->bheight,
            dst_stride,
            0, //int subsampling_x,
            0, //int subsampling_y,
            &conv_params);

        if (!blk_geom->has_uv)
            return return_error;


#if CHROMA_BLIND
        if (perform_chroma) {
#endif
         if (blk_geom->bwidth >= 16  && blk_geom->bheight >= 16 ) {
            // Cb
            src_ptr = ref_pic_list0->bufferCb + ref_pic_list0->origin_x / 2 + (ref_pic_list0->origin_y / 2) * ref_pic_list0->strideCb;
            src_stride = ref_pic_list0->strideCb;

            dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            dst_stride = prediction_ptr->strideCb;
            conv_params = get_conv_params_no_round(0, 0, 0, NULL, 64, is_compound, EB_8BIT);

            av1_warp_plane(
                wm_params,
                (int) is16bit,
                bit_depth,
                src_ptr,
                buf_width >> ss_x,
                buf_height >> ss_y,
                src_stride,
                dst_ptr,
                pu_origin_x >> ss_x,
                pu_origin_y >> ss_y,
                blk_geom->bwidth >> ss_x,
                blk_geom->bheight >> ss_y,
                dst_stride,
                ss_x,
                ss_y,
                &conv_params);

            // Cr
            src_ptr = ref_pic_list0->bufferCr + ref_pic_list0->origin_x / 2 + (ref_pic_list0->origin_y / 2 ) * ref_pic_list0->strideCr;
            src_stride = ref_pic_list0->strideCr;

            dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            dst_stride = prediction_ptr->strideCr;

            conv_params = get_conv_params_no_round(0, 0, 0, NULL, 64, is_compound, EB_8BIT);

            av1_warp_plane(
                wm_params,
                (int) is16bit,
                bit_depth,
                src_ptr,
                (int) buf_width >> ss_x,
                (int) buf_height >> ss_y,
                src_stride,
                dst_ptr,
                pu_origin_x >> ss_x,
                pu_origin_y >> ss_y,
                blk_geom->bwidth >> ss_x,
                blk_geom->bheight >> ss_y,
                dst_stride,
                ss_x,
                ss_y,
                &conv_params);
        } else { // Translation prediction when chroma block is smaller than 8x8
            DECLARE_ALIGNED(32, uint16_t, tmp_dstCb[64 * 64]);
            DECLARE_ALIGNED(32, uint16_t, tmp_dstCr[64 * 64]);
            InterpFilterParams filter_params_x, filter_params_y;
            const uint32_t interp_filters = 0;
            MV  mv, mv_q4;
            int32_t subpel_x, subpel_y;

            mv.col = mv_unit->mv[REF_LIST_0].x;
            mv.row = mv_unit->mv[REF_LIST_0].y;

            //List0-Cb
            src_ptr = ref_pic_list0->bufferCb + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCb;
            dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            src_stride = ref_pic_list0->strideCb;
            dst_stride = prediction_ptr->strideCb;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, 64, is_compound, EB_8BIT);

            av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                &filter_params_x,
                &filter_params_y,
                subpel_x,
                subpel_y,
                &conv_params);

            //List0-Cr
            src_ptr = ref_pic_list0->bufferCr + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCr;
            dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            src_stride = ref_pic_list0->strideCr;
            dst_stride = prediction_ptr->strideCr;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, 64, is_compound, EB_8BIT);
            convolve[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                &filter_params_x,
                &filter_params_y,
                subpel_x,
                subpel_y,
                &conv_params);
            }
#if CHROMA_BLIND
        }
#endif
    } else { // HBD
        uint16_t *src_ptr;
        uint16_t *dst_ptr;

        // Y - UNI_PRED_LIST_0
        src_ptr = (uint16_t *)ref_pic_list0->buffer_y + ref_pic_list0->origin_x + ref_pic_list0->origin_y * ref_pic_list0->stride_y;
        src_stride = ref_pic_list0->stride_y;
        buf_width = ref_pic_list0->width;
        buf_height = ref_pic_list0->height;

        dst_ptr = (uint16_t *)prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
        dst_stride = prediction_ptr->stride_y;
        conv_params = get_conv_params_no_round(0, 0, 0, NULL, 128, is_compound, bit_depth);

        av1_warp_plane_hbd(
            wm_params,
            bit_depth,
            src_ptr,
            (int) buf_width,
            (int) buf_height,
            src_stride,
            dst_ptr,
            pu_origin_x,
            pu_origin_y,
            blk_geom->bwidth,
            blk_geom->bheight,
            dst_stride,
            0, //int subsampling_x,
            0, //int subsampling_y,
            &conv_params);

        if (!blk_geom->has_uv)
            return return_error;

#if CHROMA_BLIND
        if (perform_chroma) {
#endif
         if (blk_geom->bwidth >= 16  && blk_geom->bheight >= 16 ) {
            // Cb
            src_ptr = (uint16_t *)ref_pic_list0->bufferCb + ref_pic_list0->origin_x / 2 + (ref_pic_list0->origin_y / 2) * ref_pic_list0->strideCb;
            src_stride = ref_pic_list0->strideCb;

            dst_ptr = (uint16_t *)prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            dst_stride = prediction_ptr->strideCb;
            conv_params = get_conv_params_no_round(0, 0, 0, NULL, 64, is_compound, bit_depth);

            av1_warp_plane_hbd(
                wm_params,
                bit_depth,
                src_ptr,
                buf_width >> ss_x,
                buf_height >> ss_y,
                src_stride,
                dst_ptr,
                pu_origin_x >> ss_x,
                pu_origin_y >> ss_y,
                blk_geom->bwidth >> ss_x,
                blk_geom->bheight >> ss_y,
                dst_stride,
                ss_x,
                ss_y,
                &conv_params);

            // Cr
            src_ptr = (uint16_t *)ref_pic_list0->bufferCr + ref_pic_list0->origin_x / 2 + (ref_pic_list0->origin_y / 2 ) * ref_pic_list0->strideCr;
            src_stride = ref_pic_list0->strideCr;

            dst_ptr = (uint16_t *)prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            dst_stride = prediction_ptr->strideCr;

            conv_params = get_conv_params_no_round(0, 0, 0, NULL, 64, is_compound, bit_depth);

            av1_warp_plane_hbd(
                wm_params,
                bit_depth,
                src_ptr,
                (int) buf_width >> ss_x,
                (int) buf_height >> ss_y,
                src_stride,
                dst_ptr,
                pu_origin_x >> ss_x,
                pu_origin_y >> ss_y,
                blk_geom->bwidth >> ss_x,
                blk_geom->bheight >> ss_y,
                dst_stride,
                ss_x,
                ss_y,
                &conv_params);
        } else { // Simple translation prediction when chroma block is smaller than 8x8
            DECLARE_ALIGNED(32, uint16_t, tmp_dstCb[64 * 64]);
            DECLARE_ALIGNED(32, uint16_t, tmp_dstCr[64 * 64]);
            InterpFilterParams filter_params_x, filter_params_y;
            const uint32_t interp_filters = 0;
            MV  mv, mv_q4;
            int32_t subpel_x, subpel_y;

            mv.col = mv_unit->mv[REF_LIST_0].x;
            mv.row = mv_unit->mv[REF_LIST_0].y;

            //List0-Cb
            src_ptr = (uint16_t *)ref_pic_list0->bufferCb + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCb;
            dst_ptr = (uint16_t *)prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
            src_stride = ref_pic_list0->strideCb;
            dst_stride = prediction_ptr->strideCb;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, 64, is_compound, bit_depth);

            av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

            convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                &filter_params_x,
                &filter_params_y,
                subpel_x,
                subpel_y,
                &conv_params,
                bit_depth);

            //List0-Cr
            src_ptr = (uint16_t *)ref_pic_list0->bufferCr + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCr;
            dst_ptr = (uint16_t *)prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
            src_stride = ref_pic_list0->strideCr;
            dst_stride = prediction_ptr->strideCr;

            mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
            subpel_x = mv_q4.col & SUBPEL_MASK;
            subpel_y = mv_q4.row & SUBPEL_MASK;
            src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
            conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, 64, is_compound, bit_depth);
            convolveHbd[subpel_x != 0][subpel_y != 0][is_compound](
                src_ptr,
                src_stride,
                dst_ptr,
                dst_stride,
                blk_geom->bwidth_uv,
                blk_geom->bheight_uv,
                &filter_params_x,
                &filter_params_y,
                subpel_x,
                subpel_y,
                &conv_params,
                bit_depth);
            }
#if CHROMA_BLIND
        }
#endif
    }

    return return_error;
}


EbErrorType warped_motion_prediction_md(
    MvUnit_t                               *mv_unit,
    ModeDecisionContext_t                  *md_context_ptr,
    uint16_t                                pu_origin_x,
    uint16_t                                pu_origin_y,
    CodingUnit_t                           *cu_ptr,
    const BlockGeom                        *blk_geom,
    EbPictureBufferDesc_t                  *ref_pic_list0,
    EbPictureBufferDesc_t                  *prediction_ptr,
    uint16_t                                dst_origin_x,
    uint16_t                                dst_origin_y,
    EbWarpedMotionParams                   *wm_params,
    EbAsm                                   asm_type)
{
    EbErrorType  return_error = EB_ErrorNone;
    InterPredictionContext_t *context_ptr = (InterPredictionContext_t*)(md_context_ptr->inter_prediction_context);
    uint8_t is_compound = (mv_unit->predDirection == BI_PRED) ? 1 : 0;
    assert(!is_compound);

    int32_t src_stride;
    int32_t dst_stride;
    uint16_t buf_width;
    uint16_t buf_height;
    ConvolveParams conv_params;
    uint8_t ss_x = 1; // subsamplings
    uint8_t ss_y = 1;

    uint16_t *src_ptr;
    uint8_t *dst_ptr;

    // Y - UNI_PRED_LIST_0
    src_ptr = (uint16_t *)ref_pic_list0->buffer_y + ref_pic_list0->origin_x + ref_pic_list0->origin_y * ref_pic_list0->stride_y;
    src_stride = ref_pic_list0->stride_y;
    buf_width = ref_pic_list0->width;
    buf_height = ref_pic_list0->height;

    dst_ptr = prediction_ptr->buffer_y + prediction_ptr->origin_x + dst_origin_x + (prediction_ptr->origin_y + dst_origin_y) * prediction_ptr->stride_y;
    dst_stride = prediction_ptr->stride_y;
    conv_params = get_conv_params_no_round(0, 0, 0, NULL, 128, is_compound, EB_8BIT);

    Av1UnPackReferenceBlock(
        src_ptr,
        ref_pic_list0->stride_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        context_ptr->mcp_context->local_reference_block8_bitl0->buffer_y,
        context_ptr->mcp_context->local_reference_block8_bitl0->stride_y,
        EB_FALSE,
        asm_type,
        16);

    av1_warp_plane(
        wm_params,
        0,  // int use_hbd
        8,  // int bd
        context_ptr->mcp_context->local_reference_block8_bitl0->buffer_y + 8 + (8 * context_ptr->mcp_context->local_reference_block8_bitl0->stride_y),
        (int) buf_width,
        (int) buf_height,
        context_ptr->mcp_context->local_reference_block8_bitl0->stride_y,
        dst_ptr,
        pu_origin_x,
        pu_origin_y,
        blk_geom->bwidth,
        blk_geom->bheight,
        dst_stride,
        0, // int subsampling_x
        0, // int subsampling_y
        &conv_params);

    if (!blk_geom->has_uv)
        return return_error;
#if CHROMA_BLIND
    if (md_context_ptr->chroma_level == CHROMA_MODE_0) {
#endif

     if (blk_geom->bwidth >= 16  && blk_geom->bheight >= 16 ) {
        // Cb
        src_ptr = (uint16_t *)ref_pic_list0->bufferCb + ref_pic_list0->origin_x / 2 + (ref_pic_list0->origin_y / 2) * ref_pic_list0->strideCb;
        src_stride = ref_pic_list0->strideCb;

        dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
        dst_stride = prediction_ptr->strideCb;
        conv_params = get_conv_params_no_round(0, 0, 0, NULL, 64, is_compound, EB_8BIT);

        Av1UnPackReferenceBlock(
            src_ptr,
            ref_pic_list0->strideCb,
            blk_geom->bwidth_uv,
            blk_geom->bheight_uv,
            context_ptr->mcp_context->local_reference_block8_bitl0->bufferCb,
            context_ptr->mcp_context->local_reference_block8_bitl0->strideCb,
            EB_FALSE,
            asm_type,
            8);

        av1_warp_plane(
            wm_params,
            0, // int use_hbd
            8, // int bd
            context_ptr->mcp_context->local_reference_block8_bitl0->bufferCb + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl0->strideCb),
            buf_width >> ss_x,
            buf_height >> ss_y,
            context_ptr->mcp_context->local_reference_block8_bitl0->strideCb,
            dst_ptr,
            pu_origin_x >> ss_x,
            pu_origin_y >> ss_y,
            blk_geom->bwidth >> ss_x,
            blk_geom->bheight >> ss_y,
            dst_stride,
            ss_x,
            ss_y,
            &conv_params);

        // Cr
        src_ptr = (uint16_t *)ref_pic_list0->bufferCr + ref_pic_list0->origin_x / 2 + (ref_pic_list0->origin_y / 2 ) * ref_pic_list0->strideCr;
        src_stride = ref_pic_list0->strideCr;

        dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
        dst_stride = prediction_ptr->strideCr;


            conv_params = get_conv_params_no_round(0, 0, 0, NULL, 64, is_compound, EB_8BIT);

        Av1UnPackReferenceBlock(
            src_ptr,
            ref_pic_list0->strideCr,
            blk_geom->bwidth_uv,
            blk_geom->bheight_uv,
            context_ptr->mcp_context->local_reference_block8_bitl0->bufferCr,
            context_ptr->mcp_context->local_reference_block8_bitl0->strideCr,
            EB_FALSE,
            asm_type,
            8);

        av1_warp_plane(
            wm_params,
            0, // int use_hbd
            8, // int bd
            context_ptr->mcp_context->local_reference_block8_bitl0->bufferCr + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl1->strideCr),
            (int) buf_width >> ss_x,
            (int) buf_height >> ss_y,
            context_ptr->mcp_context->local_reference_block8_bitl0->strideCr,
            dst_ptr,
            pu_origin_x >> ss_x,
            pu_origin_y >> ss_y,
            blk_geom->bwidth >> ss_x,
            blk_geom->bheight >> ss_y,
            dst_stride,
            ss_x,
            ss_y,
            &conv_params);
    } else { // Simple translation prediction when chroma block is smaller than 8x8
        DECLARE_ALIGNED(32, uint16_t, tmp_dstCb[64 * 64]);
        DECLARE_ALIGNED(32, uint16_t, tmp_dstCr[64 * 64]);
        InterpFilterParams filter_params_x, filter_params_y;
        const uint32_t interp_filters = 0;
        MV  mv, mv_q4;
        int32_t subpel_x, subpel_y;

            mv.col = mv_unit->mv[REF_LIST_0].x;
            mv.row = mv_unit->mv[REF_LIST_0].y;
          
        //List0-Cb
        src_ptr = (uint16_t *)ref_pic_list0->bufferCb + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCb;
        dst_ptr = prediction_ptr->bufferCb + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCb;
        src_stride = ref_pic_list0->strideCb;
        dst_stride = prediction_ptr->strideCb;

        mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
        subpel_x = mv_q4.col & SUBPEL_MASK;
        subpel_y = mv_q4.row & SUBPEL_MASK;
        src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
        conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCb, 64, is_compound, EB_8BIT);


            av1_get_convolve_filter_params(interp_filters, &filter_params_x,
                &filter_params_y, blk_geom->bwidth_uv, blk_geom->bheight_uv);

        Av1UnPackReferenceBlock(
            src_ptr,
            ref_pic_list0->strideCb,
            blk_geom->bwidth_uv,
            blk_geom->bheight_uv,
            context_ptr->mcp_context->local_reference_block8_bitl0->bufferCb,
            context_ptr->mcp_context->local_reference_block8_bitl0->strideCb,
            EB_FALSE,
            asm_type,
            8);

        convolve[subpel_x != 0][subpel_y != 0][is_compound](
            context_ptr->mcp_context->local_reference_block8_bitl0->bufferCb + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl0->strideCb),
            context_ptr->mcp_context->local_reference_block8_bitl0->strideCb,
            dst_ptr,
            dst_stride,
            blk_geom->bwidth_uv,
            blk_geom->bheight_uv,
            &filter_params_x,
            &filter_params_y,
            subpel_x,
            subpel_y,
            &conv_params);

        //List0-Cr
        src_ptr = (uint16_t *)ref_pic_list0->bufferCr + (ref_pic_list0->origin_x + ((pu_origin_x >> 3) << 3)) / 2 + (ref_pic_list0->origin_y + ((pu_origin_y >> 3) << 3)) / 2 * ref_pic_list0->strideCr;
        dst_ptr = prediction_ptr->bufferCr + (prediction_ptr->origin_x + ((dst_origin_x >> 3) << 3)) / 2 + (prediction_ptr->origin_y + ((dst_origin_y >> 3) << 3)) / 2 * prediction_ptr->strideCr;
        src_stride = ref_pic_list0->strideCr;
        dst_stride = prediction_ptr->strideCr;

        mv_q4 = clamp_mv_to_umv_border_sb(cu_ptr->av1xd, &mv, blk_geom->bwidth_uv, blk_geom->bheight_uv, 1, 1);
        subpel_x = mv_q4.col & SUBPEL_MASK;
        subpel_y = mv_q4.row & SUBPEL_MASK;
        src_ptr = src_ptr + (mv_q4.row >> SUBPEL_BITS) * src_stride + (mv_q4.col >> SUBPEL_BITS);
        conv_params = get_conv_params_no_round(0, 0, 0, tmp_dstCr, 64, is_compound, EB_8BIT);

        Av1UnPackReferenceBlock(
            src_ptr,
            ref_pic_list0->strideCr,
            blk_geom->bwidth_uv,
            blk_geom->bheight_uv,
            context_ptr->mcp_context->local_reference_block8_bitl0->bufferCr,
            context_ptr->mcp_context->local_reference_block8_bitl0->strideCr,
            EB_FALSE,
            asm_type,
            8);

        convolve[subpel_x != 0][subpel_y != 0][is_compound](
            context_ptr->mcp_context->local_reference_block8_bitl0->bufferCr + 4 + (4 * context_ptr->mcp_context->local_reference_block8_bitl1->strideCr),
            context_ptr->mcp_context->local_reference_block8_bitl0->strideCr,
            dst_ptr,
            dst_stride,
            blk_geom->bwidth_uv,
            blk_geom->bheight_uv,
            &filter_params_x,
            &filter_params_y,
            subpel_x,
            subpel_y,
            &conv_params);
        }
#if CHROMA_BLIND
    }
#endif
    return return_error;
}


#define SWITCHABLE_INTERP_RATE_FACTOR 1
extern int32_t av1_get_pred_context_switchable_interp(
    NeighborArrayUnit_t     *ref_frame_type_neighbor_array,
    MvReferenceFrame rf0,
    MvReferenceFrame rf1,
    NeighborArrayUnit32_t     *interpolation_type_neighbor_array,
    uint32_t cu_origin_x,
    uint32_t cu_origin_y,
    int32_t dir
);


int32_t av1_get_switchable_rate(
    ModeDecisionCandidateBuffer_t *candidate_buffer_ptr,
    const Av1Common *const cm,
    ModeDecisionContext_t *md_context_ptr//,
    // Macroblock *x,
    // const MacroBlockD *xd
) {
    if (cm->interp_filter == SWITCHABLE) {
        // const MbModeInfo *const mbmi = xd->mi[0];
        int32_t inter_filter_cost = 0;
        int32_t dir;

        for (dir = 0; dir < 2; ++dir) {

            MvReferenceFrame rf[2];
            av1_set_ref_frame(rf, candidate_buffer_ptr->candidate_ptr->ref_frame_type);
            const int32_t ctx = av1_get_pred_context_switchable_interp(
                md_context_ptr->ref_frame_type_neighbor_array,
                rf[0],
                rf[1],
                md_context_ptr->interpolation_type_neighbor_array,
                md_context_ptr->cu_origin_x,
                md_context_ptr->cu_origin_y,
                //xd,
                dir
            );

            const InterpFilter filter = av1_extract_interp_filter(/*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters, dir);
            ASSERT(ctx < SWITCHABLE_FILTER_CONTEXTS);
            ASSERT(filter < SWITCHABLE_FILTERS);
            inter_filter_cost +=  /*x->switchable_interp_costs*/md_context_ptr->md_rate_estimation_ptr->switchable_interp_FacBitss[ctx][filter];
        }
        return SWITCHABLE_INTERP_RATE_FACTOR * inter_filter_cost;
    }
    else {
        return 0;
    }
}
//void model_rd_norm(int32_t xsq_q10, int32_t *r_q10, int32_t *d_q10) {
 // NOTE: The tables below must be of the same size.

 // The functions described below are sampled at the four most significant
 // bits of x^2 + 8 / 256.




void highbd_variance64_c(const uint8_t *a8, int32_t a_stride,
    const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h,
    uint64_t *sse) {
    const uint8_t *a = a8;//CONVERT_TO_SHORTPTR(a8);
    const uint8_t *b = b8;//CONVERT_TO_SHORTPTR(b8);
    uint64_t tsse = 0;
    for (int32_t i = 0; i < h; ++i) {
        for (int32_t j = 0; j < w; ++j) {
            const int32_t diff = a[j] - b[j];
            tsse += (uint32_t)(diff * diff);
        }
        a += a_stride;
        b += b_stride;
    }
    *sse = tsse;
}
void highbd_8_variance(const uint8_t *a8, int32_t a_stride,
    const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h,
    uint32_t *sse) {
    uint64_t sse_long = 0;
    highbd_variance64(a8, a_stride, b8, b_stride, w, h, &sse_long);
    *sse = (uint32_t)sse_long;
}
/*static*/ /*INLINE*/ void variance4x4_64_sse4_1(uint8_t *a8, int32_t a_stride,
    uint8_t *b8, int32_t b_stride,
    uint64_t *sse, int64_t *sum) {
    __m128i u0, u1, u2, u3;
    __m128i s0, s1, s2, s3;
    __m128i t0, t1, x0, y0;
    __m128i a0, a1, a2, a3;
    __m128i b0, b1, b2, b3;
    __m128i k_one_epi16 = _mm_set1_epi16((int16_t)1);

    uint8_t *a = a8;//CONVERT_TO_SHORTPTR(a8);
    uint8_t *b = b8;//CONVERT_TO_SHORTPTR(b8);

    a0 = _mm_loadl_epi64((__m128i const *)(a + 0 * a_stride));
    a1 = _mm_loadl_epi64((__m128i const *)(a + 1 * a_stride));
    a2 = _mm_loadl_epi64((__m128i const *)(a + 2 * a_stride));
    a3 = _mm_loadl_epi64((__m128i const *)(a + 3 * a_stride));

    b0 = _mm_loadl_epi64((__m128i const *)(b + 0 * b_stride));
    b1 = _mm_loadl_epi64((__m128i const *)(b + 1 * b_stride));
    b2 = _mm_loadl_epi64((__m128i const *)(b + 2 * b_stride));
    b3 = _mm_loadl_epi64((__m128i const *)(b + 3 * b_stride));

    u0 = _mm_unpacklo_epi16(a0, a1);
    u1 = _mm_unpacklo_epi16(a2, a3);
    u2 = _mm_unpacklo_epi16(b0, b1);
    u3 = _mm_unpacklo_epi16(b2, b3);

    s0 = _mm_sub_epi16(u0, u2);
    s1 = _mm_sub_epi16(u1, u3);

    t0 = _mm_madd_epi16(s0, k_one_epi16);
    t1 = _mm_madd_epi16(s1, k_one_epi16);

    s2 = _mm_hadd_epi32(t0, t1);
    s3 = _mm_hadd_epi32(s2, s2);
    y0 = _mm_hadd_epi32(s3, s3);

    t0 = _mm_madd_epi16(s0, s0);
    t1 = _mm_madd_epi16(s1, s1);

    s2 = _mm_hadd_epi32(t0, t1);
    s3 = _mm_hadd_epi32(s2, s2);
    x0 = _mm_hadd_epi32(s3, s3);

    *sse = (uint64_t)_mm_extract_epi32(x0, 0);
    *sum = (int64_t)_mm_extract_epi32(y0, 0);
}


#define RDDIV_BITS 7
#define RDCOST(RM, R, D)                                            \
  (ROUND_POWER_OF_TWO(((uint64_t)(R)) * (RM), AV1_PROB_COST_SHIFT) + \
   ((D) * (1 << RDDIV_BITS)))



static void model_rd_norm(int32_t xsq_q10, int32_t *r_q10, int32_t *d_q10) {
    // NOTE: The tables below must be of the same size.

    // The functions described below are sampled at the four most significant
    // bits of x^2 + 8 / 256.

    // Normalized rate:
    // This table models the rate for a Laplacian source with given variance
    // when quantized with a uniform quantizer with given stepsize. The
    // closed form expression is:
    // Rn(x) = H(sqrt(r)) + sqrt(r)*[1 + H(r)/(1 - r)],
    // where r = exp(-sqrt(2) * x) and x = qpstep / sqrt(variance),
    // and H(x) is the binary entropy function.
    static const int32_t rate_tab_q10[] = {
      65536, 6086, 5574, 5275, 5063, 4899, 4764, 4651, 4553, 4389, 4255, 4142,
      4044,  3958, 3881, 3811, 3748, 3635, 3538, 3453, 3376, 3307, 3244, 3186,
      3133,  3037, 2952, 2877, 2809, 2747, 2690, 2638, 2589, 2501, 2423, 2353,
      2290,  2232, 2179, 2130, 2084, 2001, 1928, 1862, 1802, 1748, 1698, 1651,
      1608,  1530, 1460, 1398, 1342, 1290, 1243, 1199, 1159, 1086, 1021, 963,
      911,   864,  821,  781,  745,  680,  623,  574,  530,  490,  455,  424,
      395,   345,  304,  269,  239,  213,  190,  171,  154,  126,  104,  87,
      73,    61,   52,   44,   38,   28,   21,   16,   12,   10,   8,    6,
      5,     3,    2,    1,    1,    1,    0,    0,
    };
    // Normalized distortion:
    // This table models the normalized distortion for a Laplacian source
    // with given variance when quantized with a uniform quantizer
    // with given stepsize. The closed form expression is:
    // Dn(x) = 1 - 1/sqrt(2) * x / sinh(x/sqrt(2))
    // where x = qpstep / sqrt(variance).
    // Note the actual distortion is Dn * variance.
    static const int32_t dist_tab_q10[] = {
      0,    0,    1,    1,    1,    2,    2,    2,    3,    3,    4,    5,
      5,    6,    7,    7,    8,    9,    11,   12,   13,   15,   16,   17,
      18,   21,   24,   26,   29,   31,   34,   36,   39,   44,   49,   54,
      59,   64,   69,   73,   78,   88,   97,   106,  115,  124,  133,  142,
      151,  167,  184,  200,  215,  231,  245,  260,  274,  301,  327,  351,
      375,  397,  418,  439,  458,  495,  528,  559,  587,  613,  637,  659,
      680,  717,  749,  777,  801,  823,  842,  859,  874,  899,  919,  936,
      949,  960,  969,  977,  983,  994,  1001, 1006, 1010, 1013, 1015, 1017,
      1018, 1020, 1022, 1022, 1023, 1023, 1023, 1024,
    };
    static const int32_t xsq_iq_q10[] = {
      0,      4,      8,      12,     16,     20,     24,     28,     32,
      40,     48,     56,     64,     72,     80,     88,     96,     112,
      128,    144,    160,    176,    192,    208,    224,    256,    288,
      320,    352,    384,    416,    448,    480,    544,    608,    672,
      736,    800,    864,    928,    992,    1120,   1248,   1376,   1504,
      1632,   1760,   1888,   2016,   2272,   2528,   2784,   3040,   3296,
      3552,   3808,   4064,   4576,   5088,   5600,   6112,   6624,   7136,
      7648,   8160,   9184,   10208,  11232,  12256,  13280,  14304,  15328,
      16352,  18400,  20448,  22496,  24544,  26592,  28640,  30688,  32736,
      36832,  40928,  45024,  49120,  53216,  57312,  61408,  65504,  73696,
      81888,  90080,  98272,  106464, 114656, 122848, 131040, 147424, 163808,
      180192, 196576, 212960, 229344, 245728,
    };
    const int32_t tmp = (xsq_q10 >> 2) + 8;
    const int32_t k = get_msb(tmp) - 3;
    const int32_t xq = (k << 3) + ((tmp >> k) & 0x7);
    const int32_t one_q10 = 1 << 10;
    const int32_t a_q10 = ((xsq_q10 - xsq_iq_q10[xq]) << 10) >> (2 + k);
    const int32_t b_q10 = one_q10 - a_q10;
    *r_q10 = (rate_tab_q10[xq] * b_q10 + rate_tab_q10[xq + 1] * a_q10) >> 10;
    *d_q10 = (dist_tab_q10[xq] * b_q10 + dist_tab_q10[xq + 1] * a_q10) >> 10;
}

void av1_model_rd_from_var_lapndz(int64_t var, uint32_t n_log2,
    uint32_t qstep, int32_t *rate,
    int64_t *dist) {
    // This function models the rate and distortion for a Laplacian
    // source with given variance when quantized with a uniform quantizer
    // with given stepsize. The closed form expressions are in:
    // Hang and Chen, "Source Model for transform video coder and its
    // application - Part I: Fundamental Theory", IEEE Trans. Circ.
    // Sys. for Video Tech., April 1997.
    if (var == 0) {
        *rate = 0;
        *dist = 0;
    }
    else {
        int32_t d_q10, r_q10;
        static const uint32_t MAX_XSQ_Q10 = 245727;
        const uint64_t xsq_q10_64 =
            (((uint64_t)qstep * qstep << (n_log2 + 10)) + (var >> 1)) / var;
        const int32_t xsq_q10 = (int32_t)MIN(xsq_q10_64, MAX_XSQ_Q10);
        model_rd_norm(xsq_q10, &r_q10, &d_q10);
        *rate = ROUND_POWER_OF_TWO(r_q10 << n_log2, 10 - AV1_PROB_COST_SHIFT);
        *dist = (var * (int64_t)d_q10 + 512) >> 10;
    }
}

/*static*/ void model_rd_from_sse(
    block_size bsize,
    int16_t quantizer,
    //const AV1_COMP *const cpi,
    //const MacroBlockD *const xd,
    //block_size bsize,
    //int32_t plane,
    int64_t sse,
    int32_t *rate,
    int64_t *dist
)
{
    // OMK
  //const struct MacroblockdPlane *const pd = &xd->plane[plane];
    int32_t dequant_shift = 3;
    /* OMK (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) ? xd->bd - 5 :3;*/

// Fast approximate the modelling function.
    if (0/*cpi->sf.simple_model_rd_from_var*/) {
        int64_t square_error = sse;
        quantizer = quantizer >> dequant_shift;

        if (quantizer < 120)
            *rate = (int32_t)((square_error * (280 - quantizer)) >>
            (16 - AV1_PROB_COST_SHIFT));
        else
            *rate = 0;
        *dist = (square_error * quantizer) >> 8;
    }
    else {

        av1_model_rd_from_var_lapndz(sse, num_pels_log2_lookup[bsize],
            quantizer >> dequant_shift, rate,
            dist);
    }

    *dist <<= 4;
}

/*static*/ void model_rd_for_sb(

    PictureControlSet_t *picture_control_set_ptr,
    EbPictureBufferDesc_t *prediction_ptr,
    ModeDecisionContext_t *md_context_ptr,
    //const AV1_COMP *const cpi,
    //block_size bsize,
    //Macroblock *x,
    //MacroBlockD *xd,
    int32_t plane_from,
    int32_t plane_to,
    int32_t *out_rate_sum,
    int64_t *out_dist_sum,
    int32_t *skip_txfm_sb,
    int64_t *skip_sse_sb,
    int32_t *plane_rate,
    int64_t *plane_sse,
    int64_t *plane_dist) {
    // Note our transform coeffs are 8 times an orthogonal transform.
    // Hence quantizer step is also 8 times. To get effective quantizer
    // we need to divide by 8 before sending to modeling function.
    int32_t plane;
    // const int32_t ref = xd->mi[0]->ref_frame[0];

    int64_t rate_sum = 0;
    int64_t dist_sum = 0;
    int64_t total_sse = 0;


    EbPictureBufferDesc_t                  *input_picture_ptr = picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;
    const uint32_t inputOriginIndex = (md_context_ptr->cu_origin_y + input_picture_ptr->origin_y) * input_picture_ptr->stride_y + (md_context_ptr->cu_origin_x + input_picture_ptr->origin_x);
    const uint32_t inputChromaOriginIndex = ((md_context_ptr->cu_origin_y + input_picture_ptr->origin_y) * input_picture_ptr->strideCb + (md_context_ptr->cu_origin_x + input_picture_ptr->origin_x)) / 2;

    for (plane = plane_from; plane <= plane_to; ++plane) {
        // struct MacroblockPlane *const p = &x->plane[plane];
         // struct MacroblockdPlane *const pd = &xd->plane[plane];
         // const block_size bs = get_plane_block_size(bsize, pd);
        uint32_t sse;
        int32_t rate;

        int64_t dist;

        // if (x->skip_chroma_rd && plane) continue;

         // TODO(geza): Write direct sse functions that do not compute
         // variance as well.
        uint32_t offset;


        if (plane)
            offset = (prediction_ptr->origin_x + md_context_ptr->blk_geom->origin_x + (prediction_ptr->origin_y + md_context_ptr->blk_geom->origin_y) * prediction_ptr->strideCb) / 2;
        else
            offset = prediction_ptr->origin_x + md_context_ptr->blk_geom->origin_x + (prediction_ptr->origin_y + md_context_ptr->blk_geom->origin_y) * prediction_ptr->stride_y;

        highbd_8_variance(
            plane == 0 ? (&(input_picture_ptr->buffer_y[inputOriginIndex])) : plane == 1 ? (&(input_picture_ptr->bufferCb[inputChromaOriginIndex])) : (&(input_picture_ptr->bufferCr[inputChromaOriginIndex])),
            plane == 0 ? input_picture_ptr->stride_y : plane == 1 ? input_picture_ptr->strideCb : input_picture_ptr->strideCr,
            plane == 0 ? (prediction_ptr->buffer_y + offset) : plane == 1 ? (prediction_ptr->bufferCb + offset) : (prediction_ptr->bufferCr + offset),
            plane == 0 ? prediction_ptr->stride_y : plane == 1 ? prediction_ptr->strideCb : prediction_ptr->strideCr,
            plane == 0 ? md_context_ptr->blk_geom->bwidth : md_context_ptr->blk_geom->bwidth_uv,
            plane == 0 ? md_context_ptr->blk_geom->bheight : md_context_ptr->blk_geom->bheight_uv,
            &sse
        );


        total_sse += sse;

        int32_t current_q_index = MAX(0, MIN(QINDEX_RANGE - 1, picture_control_set_ptr->parent_pcs_ptr->base_qindex));
        Dequants *const dequants = &picture_control_set_ptr->parent_pcs_ptr->deq;

        int16_t quantizer = dequants->y_dequant_Q3[current_q_index][1];
        model_rd_from_sse(
            plane == 0 ? md_context_ptr->blk_geom->bsize : md_context_ptr->blk_geom->bsize_uv,
            quantizer,
            sse,
            &rate,
            &dist);

        rate_sum += rate;
        dist_sum += dist;
        if (plane_rate) plane_rate[plane] = rate;
        if (plane_sse) plane_sse[plane] = sse;
        if (plane_dist) plane_dist[plane] = dist;
    }

    *skip_txfm_sb = total_sse == 0;
    *skip_sse_sb = total_sse << 4;
    *out_rate_sum = (int32_t)rate_sum;
    *out_dist_sum = dist_sum;
}

/*static*/ /*INLINE*/ int32_t is_nontrans_global_motion(
    block_size sb_type,
    ModeDecisionCandidateBuffer_t *candidate_buffer_ptr,
    PictureControlSet_t *picture_control_set_ptr
) {
    int32_t ref;

    // First check if all modes are GLOBALMV
    if (candidate_buffer_ptr->candidate_ptr->pred_mode != GLOBALMV && candidate_buffer_ptr->candidate_ptr->pred_mode != GLOBAL_GLOBALMV)
        return 0;

    if (MIN(mi_size_wide[sb_type], mi_size_high[sb_type]) < 2)
        return 0;
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, candidate_buffer_ptr->candidate_ptr->ref_frame_type);
    // Now check if all global motion is non translational
    for (ref = 0; ref < 1 + candidate_buffer_ptr->candidate_ptr->is_compound/*has_second_ref(mbmi)*/; ++ref) {
        if (picture_control_set_ptr->parent_pcs_ptr->global_motion[ref ? rf[1] : rf[0]].wmtype == TRANSLATION)
            //if (xd->global_motion[mbmi->ref_frame[ref]].wmtype == TRANSLATION)
            return 0;
    }
    return 1;
}
static INLINE int32_t av1_is_interp_needed(
    ModeDecisionCandidateBuffer_t *candidate_buffer_ptr,
    PictureControlSet_t *picture_control_set_ptr,
    block_size bsize)
{
    if (candidate_buffer_ptr->candidate_ptr->merge_flag)
        return 0;

    if (candidate_buffer_ptr->candidate_ptr->motion_mode == WARPED_CAUSAL)
        return 0;

    if (is_nontrans_global_motion( bsize,
        candidate_buffer_ptr, picture_control_set_ptr))
        return 0;

    return 1;
}

#define DUAL_FILTER_SET_SIZE (SWITCHABLE_FILTERS * SWITCHABLE_FILTERS)
static const int32_t filter_sets[DUAL_FILTER_SET_SIZE][2] = {
  { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 0 }, { 1, 1 },
  { 1, 2 }, { 2, 0 }, { 2, 1 }, { 2, 2 },
};

/*static*/ void interpolation_filter_search(
    PictureControlSet_t *picture_control_set_ptr,
    EbPictureBufferDesc_t *prediction_ptr,
    ModeDecisionContext_t *md_context_ptr,
    ModeDecisionCandidateBuffer_t *candidate_buffer_ptr,
    MvUnit_t mv_unit,
    EbPictureBufferDesc_t  *ref_pic_list0,
    EbPictureBufferDesc_t  *ref_pic_list1,
    EbAsm asm_type,
    //Macroblock *const xd,
    //const AV1_COMP *const cpi,
    //block_size bsize,
    //int32_t mi_row,
    //int32_t mi_col,
    //const BUFFER_SET *const tmp_dst,
    //BUFFER_SET *const orig_dst,
    /* InterpFilter (*const single_filter)[REF_FRAMES],*/
    int64_t *const rd,
    int32_t *const switchable_rate,
    int32_t *const skip_txfm_sb,
    int64_t *const skip_sse_sb) {


    const Av1Common *cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;//&cpi->common;

    const int32_t num_planes = MAX_MB_PLANE;

    int32_t i;
    int32_t tmp_rate;
    int64_t tmp_dist;

    //(void)single_filter;

    InterpFilter assign_filter = SWITCHABLE;

    if (cm->interp_filter != SWITCHABLE)
        assign_filter = cm->interp_filter;

    //set_default_interp_filters(mbmi, assign_filter);
    /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters =//EIGHTTAP_REGULAR ;
        av1_broadcast_interp_filter(av1_unswitchable_filter(assign_filter));

    *switchable_rate = av1_get_switchable_rate(
        candidate_buffer_ptr,
        cm,
        md_context_ptr//,
        //x,
        //xd
    );

    av1_inter_prediction(
        picture_control_set_ptr,
        candidate_buffer_ptr->candidate_ptr->interp_filters,
        md_context_ptr->cu_ptr,
        candidate_buffer_ptr->candidate_ptr->ref_frame_type,
        &mv_unit,
        md_context_ptr->cu_origin_x,
        md_context_ptr->cu_origin_y,
        md_context_ptr->blk_geom->bwidth,
        md_context_ptr->blk_geom->bheight,
        ref_pic_list0,
        ref_pic_list1,
        prediction_ptr,
        md_context_ptr->blk_geom->origin_x,
        md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
        md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
        asm_type);


    model_rd_for_sb(
        picture_control_set_ptr,
        prediction_ptr,
        md_context_ptr,
        0,
        num_planes - 1,
        &tmp_rate,
        &tmp_dist,
        skip_txfm_sb,
        skip_sse_sb,
        NULL, NULL, NULL);

    *rd = RDCOST(md_context_ptr->full_lambda/*x->rdmult*/, *switchable_rate + tmp_rate, tmp_dist);

    if (assign_filter == SWITCHABLE) {
        // do interp_filter search

        if (av1_is_interp_needed(candidate_buffer_ptr, picture_control_set_ptr, md_context_ptr->blk_geom->bsize) /*&& av1_is_interp_search_needed(xd)*/) {

            const int32_t filter_set_size = DUAL_FILTER_SET_SIZE;
            int32_t best_in_temp = 0;
            uint32_t best_filters = 0;// mbmi->interp_filters;

#if INTERPOLATION_SEARCH_LEVELS
            if (picture_control_set_ptr->parent_pcs_ptr->interpolation_search_level &&
#else
            if (picture_control_set_ptr->parent_pcs_ptr->interpolation_filter_search_mode == 1 &&
#endif
                picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->enable_dual_filter) {
                int32_t tmp_skip_sb = 0;
                int64_t tmp_skip_sse = INT64_MAX;
                int32_t tmp_rs;
                int64_t tmp_rd;

                // default to (R,R): EIGHTTAP_REGULARxEIGHTTAP_REGULAR
                int32_t best_dual_mode = 0;
                // Find best of {R}x{R,Sm,Sh}
                // EIGHTTAP_REGULAR mode is calculated beforehand
                for (i = 1; i < SWITCHABLE_FILTERS; ++i) {
                    tmp_skip_sb = 0;
                    tmp_skip_sse = INT64_MAX;

                    /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters = (InterpFilter)
                        av1_make_interp_filters((InterpFilter)filter_sets[i][0], (InterpFilter)filter_sets[i][1]);

                    tmp_rs = av1_get_switchable_rate(
                        candidate_buffer_ptr,
                        cm,
                        md_context_ptr//,
                        //x,
                        //xd
                    );

                    //av1_build_inter_predictors_sb(
                    //                              cm,
                    //                              xd,
                    //                              mi_row,
                    //                              mi_col,
                    //                              orig_dst,
                    //                              bsize);
                    av1_inter_prediction(
                        picture_control_set_ptr,
                        candidate_buffer_ptr->candidate_ptr->interp_filters,
                        md_context_ptr->cu_ptr,
                        candidate_buffer_ptr->candidate_ptr->ref_frame_type,
                        &mv_unit,
                        md_context_ptr->cu_origin_x,
                        md_context_ptr->cu_origin_y,
                        md_context_ptr->blk_geom->bwidth,
                        md_context_ptr->blk_geom->bheight,
                        ref_pic_list0,
                        ref_pic_list1,
                        prediction_ptr,
                        md_context_ptr->blk_geom->origin_x,
                        md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
                        md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
                        asm_type);

                    model_rd_for_sb(
                        picture_control_set_ptr,
                        prediction_ptr,
                        md_context_ptr,
                        0,
                        num_planes - 1,
                        &tmp_rate,
                        &tmp_dist,
                        &tmp_skip_sb,
                        &tmp_skip_sse,
                        NULL, NULL, NULL);
                    tmp_rd = RDCOST(md_context_ptr->full_lambda/*x->rdmult*/, tmp_rs + tmp_rate, tmp_dist);

                    if (tmp_rd < *rd) {
                        best_dual_mode = i;

                        *rd = tmp_rd;
                        *switchable_rate = tmp_rs;
                        best_filters = /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters;
                        *skip_txfm_sb = tmp_skip_sb;
                        *skip_sse_sb = tmp_skip_sse;
                        best_in_temp = !best_in_temp;
                        /*if (best_in_temp) {
                          restore_dst_buf(xd, *orig_dst, num_planes);
                        } else {
                          restore_dst_buf(xd, *tmp_dst, num_planes);
                        }*/
                    }
                }

                // From best of horizontal EIGHTTAP_REGULAR modes, check vertical modes
                for (i = best_dual_mode + SWITCHABLE_FILTERS; i < filter_set_size;
                    i += SWITCHABLE_FILTERS) {
                    tmp_skip_sb = 0;
                    tmp_skip_sse = INT64_MAX;

                    /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters =
                        av1_make_interp_filters((InterpFilter)filter_sets[i][0], (InterpFilter)filter_sets[i][1]);

                    tmp_rs = av1_get_switchable_rate(
                        candidate_buffer_ptr,
                        cm,
                        md_context_ptr//,
                        //x,
                        //xd
                    );
                    //av1_build_inter_predictors_sb(
                    //                              cm,
                    //                              xd,
                    //                              mi_row,
                    //                              mi_col,
                    //                              orig_dst,
                    //                              bsize);

                    av1_inter_prediction(
                        picture_control_set_ptr,
                        candidate_buffer_ptr->candidate_ptr->interp_filters,
                        md_context_ptr->cu_ptr,
                        candidate_buffer_ptr->candidate_ptr->ref_frame_type,
                        &mv_unit,
                        md_context_ptr->cu_origin_x,
                        md_context_ptr->cu_origin_y,
                        md_context_ptr->blk_geom->bwidth,
                        md_context_ptr->blk_geom->bheight,
                        ref_pic_list0,
                        ref_pic_list1,
                        prediction_ptr,
                        md_context_ptr->blk_geom->origin_x,
                        md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
                        md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
                        asm_type);


                    model_rd_for_sb(
                        picture_control_set_ptr,
                        prediction_ptr,
                        md_context_ptr,
                        0,
                        num_planes - 1,
                        &tmp_rate,
                        &tmp_dist,
                        &tmp_skip_sb,
                        &tmp_skip_sse,
                        NULL, NULL, NULL);
                    tmp_rd = RDCOST(md_context_ptr->full_lambda/*x->rdmult*/, tmp_rs + tmp_rate, tmp_dist);

                    if (tmp_rd < *rd) {
                        *rd = tmp_rd;
                        *switchable_rate = tmp_rs;
                        best_filters = /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters;
                        *skip_txfm_sb = tmp_skip_sb;
                        *skip_sse_sb = tmp_skip_sse;
                        best_in_temp = !best_in_temp;
                        /*if (best_in_temp) {
                          restore_dst_buf(xd, *orig_dst, num_planes);
                        } else {
                          restore_dst_buf(xd, *tmp_dst, num_planes);
                        }*/
                    }
                }
            }
            else {
                // EIGHTTAP_REGULAR mode is calculated beforehand
                for (i = 1; i < filter_set_size; ++i) {
                    int32_t tmp_skip_sb = 0;
                    int64_t tmp_skip_sse = INT64_MAX;
                    int32_t tmp_rs;
                    int64_t tmp_rd;

                    if (/*cm->seq_params.enable_dual_filter*/picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->enable_dual_filter == 0)
                        if (filter_sets[i][0] != filter_sets[i][1]) continue;

                    /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters = av1_make_interp_filters((InterpFilter)filter_sets[i][0], (InterpFilter)filter_sets[i][1]);

                    tmp_rs = av1_get_switchable_rate(
                        candidate_buffer_ptr,
                        cm,
                        md_context_ptr//,
                        //x,
                        //xd
                    );
                    //av1_build_inter_predictors_sb(
                    //                              cm,
                    //                              xd,
                    //                              mi_row,
                    //                              mi_col,
                    //                              orig_dst,
                    //                              bsize);

                    av1_inter_prediction(
                        picture_control_set_ptr,
                        candidate_buffer_ptr->candidate_ptr->interp_filters,
                        md_context_ptr->cu_ptr,
                        candidate_buffer_ptr->candidate_ptr->ref_frame_type,
                        &mv_unit,
                        md_context_ptr->cu_origin_x,
                        md_context_ptr->cu_origin_y,
                        md_context_ptr->blk_geom->bwidth,
                        md_context_ptr->blk_geom->bheight,
                        ref_pic_list0,
                        ref_pic_list1,
                        prediction_ptr,
                        md_context_ptr->blk_geom->origin_x,
                        md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
                        md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
                        asm_type);

                    model_rd_for_sb(
                        picture_control_set_ptr,
                        prediction_ptr,
                        md_context_ptr,
                        0,
                        num_planes - 1,
                        &tmp_rate,
                        &tmp_dist,
                        &tmp_skip_sb,
                        &tmp_skip_sse,
                        NULL, NULL, NULL);
                    tmp_rd = RDCOST(md_context_ptr->full_lambda/*x->rdmult*/, tmp_rs + tmp_rate, tmp_dist);

                    if (tmp_rd < *rd) {
                        *rd = tmp_rd;
                        *switchable_rate = tmp_rs;
                        best_filters = /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters;
                        *skip_txfm_sb = tmp_skip_sb;
                        *skip_sse_sb = tmp_skip_sse;
                        best_in_temp = !best_in_temp;
                        /*if (best_in_temp) {
                          restore_dst_buf(xd, *orig_dst, num_planes);
                        } else {
                          restore_dst_buf(xd, *tmp_dst, num_planes);
                        }*/
                    }
                }
            }

            /*if (best_in_temp) {
              restore_dst_buf(xd, *tmp_dst, num_planes);
            } else {
              restore_dst_buf(xd, *orig_dst, num_planes);
            }*/
            /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters = best_filters;
        }
        else {
            candidate_buffer_ptr->candidate_ptr->interp_filters = 0;

            /*assert(mbmi->cu_ptr->interp_filters ==
                   av1_broadcast_interp_filter(EIGHTTAP_REGULAR));*/
        }
    }
    //  return 0;
}

#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
/*static*/ void interpolation_filter_search_HBD(
    PictureControlSet_t *picture_control_set_ptr,
    EbPictureBufferDesc_t *prediction_ptr,
    ModeDecisionContext_t *md_context_ptr,
    ModeDecisionCandidateBuffer_t *candidate_buffer_ptr,
    MvUnit_t mv_unit,
    EbPictureBufferDesc_t  *ref_pic_list0,
    EbPictureBufferDesc_t  *ref_pic_list1,
    EbAsm asm_type,
    //Macroblock *const xd,
    //const AV1_COMP *const cpi,
    //block_size bsize,
    //int32_t mi_row,
    //int32_t mi_col,
    //const BUFFER_SET *const tmp_dst,
    //BUFFER_SET *const orig_dst,
    /* InterpFilter (*const single_filter)[REF_FRAMES],*/
    int64_t *const rd,
    int32_t *const switchable_rate,
    int32_t *const skip_txfm_sb,
    int64_t *const skip_sse_sb) {


    const Av1Common *cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;//&cpi->common;

#if CHROMA_BLIND
    const int32_t num_planes = (md_context_ptr->blk_geom->has_uv && md_context_ptr->chroma_level == CHROMA_MODE_0) ? MAX_MB_PLANE : 1;
#else
    const int32_t num_planes = md_context_ptr->blk_geom->has_uv ? MAX_MB_PLANE : 1;
#endif
    int32_t i;
    int32_t tmp_rate;
    int64_t tmp_dist;

    //(void)single_filter;

    InterpFilter assign_filter = SWITCHABLE;

    if (cm->interp_filter != SWITCHABLE)
        assign_filter = cm->interp_filter;

    //set_default_interp_filters(mbmi, assign_filter);
    /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters =//EIGHTTAP_REGULAR ;
        av1_broadcast_interp_filter(av1_unswitchable_filter(assign_filter));

    *switchable_rate = av1_get_switchable_rate(
        candidate_buffer_ptr,
        cm,
        md_context_ptr//,
        //x,
        //xd
    );

    AV1InterPrediction10BitMD(
        candidate_buffer_ptr->candidate_ptr->interp_filters,
        picture_control_set_ptr,
        candidate_buffer_ptr->candidate_ptr->ref_frame_type,
        md_context_ptr,
        md_context_ptr->cu_ptr,
        &mv_unit,
        md_context_ptr->cu_origin_x,
        md_context_ptr->cu_origin_y,
        md_context_ptr->blk_geom->bwidth,
        md_context_ptr->blk_geom->bheight,
        ref_pic_list0,
        ref_pic_list1,
        prediction_ptr,
        md_context_ptr->blk_geom->origin_x,
        md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
        md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
        asm_type);

    model_rd_for_sb(
        picture_control_set_ptr,
        prediction_ptr,
        md_context_ptr,
        0,
        num_planes - 1,
        &tmp_rate,
        &tmp_dist,
        skip_txfm_sb,
        skip_sse_sb,
        NULL, NULL, NULL);

    *rd = RDCOST(md_context_ptr->full_lambda/*x->rdmult*/, *switchable_rate + tmp_rate, tmp_dist);

    if (assign_filter == SWITCHABLE) {
        // do interp_filter search
        if (av1_is_interp_needed(candidate_buffer_ptr, picture_control_set_ptr, md_context_ptr->blk_geom->bsize) /*&& av1_is_interp_search_needed(xd)*/) {
            const int32_t filter_set_size = DUAL_FILTER_SET_SIZE;
            int32_t best_in_temp = 0;
            uint32_t best_filters = 0;// mbmi->interp_filters;

#if INTERPOLATION_SEARCH_LEVELS
            if (picture_control_set_ptr->parent_pcs_ptr->interpolation_search_level &&
#else
            if (picture_control_set_ptr->parent_pcs_ptr->interpolation_filter_search_mode == 1 &&
#endif
                picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->enable_dual_filter) {
                int32_t tmp_skip_sb = 0;
                int64_t tmp_skip_sse = INT64_MAX;
                int32_t tmp_rs;
                int64_t tmp_rd;

                // default to (R,R): EIGHTTAP_REGULARxEIGHTTAP_REGULAR
                int32_t best_dual_mode = 0;
                // Find best of {R}x{R,Sm,Sh}
                // EIGHTTAP_REGULAR mode is calculated beforehand
                for (i = 1; i < SWITCHABLE_FILTERS; ++i) {
                    tmp_skip_sb = 0;
                    tmp_skip_sse = INT64_MAX;

                    /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters = (InterpFilter)
                        av1_make_interp_filters((InterpFilter)filter_sets[i][0], (InterpFilter)filter_sets[i][1]);

                    tmp_rs = av1_get_switchable_rate(
                        candidate_buffer_ptr,
                        cm,
                        md_context_ptr//,
                        //x,
                        //xd
                    );

                    //av1_build_inter_predictors_sb(
                    //                              cm,
                    //                              xd,
                    //                              mi_row,
                    //                              mi_col,
                    //                              orig_dst,
                    //                              bsize);

                    AV1InterPrediction10BitMD(
                        candidate_buffer_ptr->candidate_ptr->interp_filters,
                        picture_control_set_ptr,
                        candidate_buffer_ptr->candidate_ptr->ref_frame_type,
                        md_context_ptr,
                        md_context_ptr->cu_ptr,
                        &mv_unit,
                        md_context_ptr->cu_origin_x,
                        md_context_ptr->cu_origin_y,
                        md_context_ptr->blk_geom->bwidth,
                        md_context_ptr->blk_geom->bheight,
                        ref_pic_list0,
                        ref_pic_list1,
                        prediction_ptr,
                        md_context_ptr->blk_geom->origin_x,
                        md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
                        md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
                        asm_type);

                    model_rd_for_sb(
                        picture_control_set_ptr,
                        prediction_ptr,
                        md_context_ptr,
                        0,
                        num_planes - 1,
                        &tmp_rate,
                        &tmp_dist,
                        &tmp_skip_sb,
                        &tmp_skip_sse,
                        NULL, NULL, NULL);
                    tmp_rd = RDCOST(md_context_ptr->full_lambda/*x->rdmult*/, tmp_rs + tmp_rate, tmp_dist);

                    if (tmp_rd < *rd) {
                        best_dual_mode = i;

                        *rd = tmp_rd;
                        *switchable_rate = tmp_rs;
                        best_filters = /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters;
                        *skip_txfm_sb = tmp_skip_sb;
                        *skip_sse_sb = tmp_skip_sse;
                        best_in_temp = !best_in_temp;
                        /*if (best_in_temp) {
                        restore_dst_buf(xd, *orig_dst, num_planes);
                        } else {
                        restore_dst_buf(xd, *tmp_dst, num_planes);
                        }*/
                    }
                }

                // From best of horizontal EIGHTTAP_REGULAR modes, check vertical modes
                for (i = best_dual_mode + SWITCHABLE_FILTERS; i < filter_set_size;
                    i += SWITCHABLE_FILTERS) {
                    tmp_skip_sb = 0;
                    tmp_skip_sse = INT64_MAX;

                    /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters =
                        av1_make_interp_filters((InterpFilter)filter_sets[i][0], (InterpFilter)filter_sets[i][1]);

                    tmp_rs = av1_get_switchable_rate(
                        candidate_buffer_ptr,
                        cm,
                        md_context_ptr//,
                        //x,
                        //xd
                    );
                    //av1_build_inter_predictors_sb(
                    //                              cm,
                    //                              xd,
                    //                              mi_row,
                    //                              mi_col,
                    //                              orig_dst,
                    //                              bsize);

                    AV1InterPrediction10BitMD(
                        candidate_buffer_ptr->candidate_ptr->interp_filters,
                        picture_control_set_ptr,
                        candidate_buffer_ptr->candidate_ptr->ref_frame_type,
                        md_context_ptr,
                        md_context_ptr->cu_ptr,
                        &mv_unit,
                        md_context_ptr->cu_origin_x,
                        md_context_ptr->cu_origin_y,
                        md_context_ptr->blk_geom->bwidth,
                        md_context_ptr->blk_geom->bheight,
                        ref_pic_list0,
                        ref_pic_list1,
                        prediction_ptr,
                        md_context_ptr->blk_geom->origin_x,
                        md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
                        md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
                        asm_type);


                    model_rd_for_sb(
                        picture_control_set_ptr,
                        prediction_ptr,
                        md_context_ptr,
                        0,
                        num_planes - 1,
                        &tmp_rate,
                        &tmp_dist,
                        &tmp_skip_sb,
                        &tmp_skip_sse,
                        NULL, NULL, NULL);
                    tmp_rd = RDCOST(md_context_ptr->full_lambda/*x->rdmult*/, tmp_rs + tmp_rate, tmp_dist);

                    if (tmp_rd < *rd) {
                        *rd = tmp_rd;
                        *switchable_rate = tmp_rs;
                        best_filters = /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters;
                        *skip_txfm_sb = tmp_skip_sb;
                        *skip_sse_sb = tmp_skip_sse;
                        best_in_temp = !best_in_temp;
                        /*if (best_in_temp) {
                        restore_dst_buf(xd, *orig_dst, num_planes);
                        } else {
                        restore_dst_buf(xd, *tmp_dst, num_planes);
                        }*/
                    }
                }
            }
            else {
                // EIGHTTAP_REGULAR mode is calculated beforehand
                for (i = 1; i < filter_set_size; ++i) {
                    int32_t tmp_skip_sb = 0;
                    int64_t tmp_skip_sse = INT64_MAX;
                    int32_t tmp_rs;
                    int64_t tmp_rd;

                    if (/*cm->seq_params.enable_dual_filter*/picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->enable_dual_filter == 0)
                        if (filter_sets[i][0] != filter_sets[i][1]) continue;

                    /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters = av1_make_interp_filters((InterpFilter)filter_sets[i][0], (InterpFilter)filter_sets[i][1]);

                    tmp_rs = av1_get_switchable_rate(
                        candidate_buffer_ptr,
                        cm,
                        md_context_ptr//,
                        //x,
                        //xd
                    );
                    //av1_build_inter_predictors_sb(
                    //                              cm,
                    //                              xd,
                    //                              mi_row,
                    //                              mi_col,
                    //                              orig_dst,
                    //                              bsize);

                    AV1InterPrediction10BitMD(
                        candidate_buffer_ptr->candidate_ptr->interp_filters,
                        picture_control_set_ptr,
                        candidate_buffer_ptr->candidate_ptr->ref_frame_type,
                        md_context_ptr,
                        md_context_ptr->cu_ptr,
                        &mv_unit,
                        md_context_ptr->cu_origin_x,
                        md_context_ptr->cu_origin_y,
                        md_context_ptr->blk_geom->bwidth,
                        md_context_ptr->blk_geom->bheight,
                        ref_pic_list0,
                        ref_pic_list1,
                        prediction_ptr,
                        md_context_ptr->blk_geom->origin_x,
                        md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
                        md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
                        asm_type);

                    model_rd_for_sb(
                        picture_control_set_ptr,
                        prediction_ptr,
                        md_context_ptr,
                        0,
                        num_planes - 1,
                        &tmp_rate,
                        &tmp_dist,
                        &tmp_skip_sb,
                        &tmp_skip_sse,
                        NULL, NULL, NULL);
                    tmp_rd = RDCOST(md_context_ptr->full_lambda/*x->rdmult*/, tmp_rs + tmp_rate, tmp_dist);

                    if (tmp_rd < *rd) {
                        *rd = tmp_rd;
                        *switchable_rate = tmp_rs;
                        best_filters = /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters;
                        *skip_txfm_sb = tmp_skip_sb;
                        *skip_sse_sb = tmp_skip_sse;
                        best_in_temp = !best_in_temp;
                        /*if (best_in_temp) {
                        restore_dst_buf(xd, *orig_dst, num_planes);
                        } else {
                        restore_dst_buf(xd, *tmp_dst, num_planes);
                        }*/
                    }
                }
            }

            /*if (best_in_temp) {
            restore_dst_buf(xd, *tmp_dst, num_planes);
            } else {
            restore_dst_buf(xd, *orig_dst, num_planes);
            }*/
            /*mbmi*/candidate_buffer_ptr->candidate_ptr->interp_filters = best_filters;
        }
        else {
            candidate_buffer_ptr->candidate_ptr->interp_filters = 0;

            /*assert(mbmi->cu_ptr->interp_filters ==
            av1_broadcast_interp_filter(EIGHTTAP_REGULAR));*/
        }
    }
    //  return 0;
}


#endif

EbErrorType inter_pu_prediction_av1(
    ModeDecisionContext_t                  *md_context_ptr,
#if !CHROMA_BLIND
    uint32_t                                  component_mask,
#endif
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
    EbAsm                                   asm_type)
{
#if !CHROMA_BLIND
    (void)component_mask;
#endif
    EbErrorType            return_error = EB_ErrorNone;
    EbPictureBufferDesc_t  *ref_pic_list0;
    EbPictureBufferDesc_t  *ref_pic_list1 = NULL;
    ModeDecisionCandidate_t *const candidate_ptr = candidate_buffer_ptr->candidate_ptr;

    MvUnit_t mv_unit;
    mv_unit.predDirection = candidate_buffer_ptr->candidate_ptr->prediction_direction[md_context_ptr->pu_itr];
    mv_unit.mv[0].x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
    mv_unit.mv[0].y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;
    mv_unit.mv[1].x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
    mv_unit.mv[1].y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

    SequenceControlSet_t* sequence_control_set_ptr = ((SequenceControlSet_t*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr));
    EbBool  is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    int64_t skip_sse_sb = INT64_MAX;
    int32_t skip_txfm_sb = 0;
    int32_t rs = 0;
    int64_t rd = INT64_MAX;

    if (is16bit) {
        ref_pic_list0 = ((EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr)->referencePicture16bit;
        if (picture_control_set_ptr->slice_type == B_SLICE)
            ref_pic_list1 = ((EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr)->referencePicture16bit;
    } else {
        ref_pic_list0 = ((EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr)->referencePicture;
        if (picture_control_set_ptr->slice_type == B_SLICE)
            ref_pic_list1 = ((EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr)->referencePicture;
    }

    if (picture_control_set_ptr->parent_pcs_ptr->allow_warped_motion
        && candidate_ptr->motion_mode != WARPED_CAUSAL)
            wm_count_samples(
                md_context_ptr->cu_ptr,
                md_context_ptr->blk_geom,
                md_context_ptr->cu_origin_x,
                md_context_ptr->cu_origin_y,
                candidate_ptr->ref_frame_type,
                picture_control_set_ptr,
                &candidate_ptr->num_proj_ref);

    if (candidate_ptr->motion_mode == WARPED_CAUSAL) {
        if (is16bit) {
            warped_motion_prediction_md(
                &mv_unit,
                md_context_ptr,
                md_context_ptr->cu_origin_x,
                md_context_ptr->cu_origin_y,
                md_context_ptr->cu_ptr,
                md_context_ptr->blk_geom,
                ref_pic_list0,
                candidate_buffer_ptr->prediction_ptr,
                md_context_ptr->blk_geom->origin_x,
                md_context_ptr->blk_geom->origin_y,
                &candidate_ptr->wm_params,
                asm_type);
        } else {
            warped_motion_prediction(
                &mv_unit,
                md_context_ptr->cu_origin_x,
                md_context_ptr->cu_origin_y,
                md_context_ptr->cu_ptr,
                md_context_ptr->blk_geom,
                ref_pic_list0,
                candidate_buffer_ptr->prediction_ptr,
                md_context_ptr->blk_geom->origin_x,
                md_context_ptr->blk_geom->origin_y,
                &candidate_ptr->wm_params,
                (uint8_t) sequence_control_set_ptr->static_config.encoder_bit_depth,
#if CHROMA_BLIND
                md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
                asm_type);
        }
        return return_error;
    }

    if (is16bit) {
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
        candidate_buffer_ptr->candidate_ptr->interp_filters = 0;
#if INTERPOLATION_SEARCH_LEVELS
        if (!md_context_ptr->skip_interpolation_search) {
#else
        if (picture_control_set_ptr->parent_pcs_ptr->interpolation_filter_search_mode > 0) {
#endif
            if (md_context_ptr->blk_geom->bwidth > 4 && md_context_ptr->blk_geom->bheight > 4)
                interpolation_filter_search_HBD(
                    picture_control_set_ptr,
                    candidate_buffer_ptr->predictionPtrTemp,
                    md_context_ptr,
                    candidate_buffer_ptr,
                    mv_unit,
                    ref_pic_list0,
                    ref_pic_list1,
                    asm_type,
                    &rd,
                    &rs,
                    &skip_txfm_sb,
                    &skip_sse_sb);
        }

#endif
#if INTERPOL_FILTER_SEARCH_10BIT_SUPPORT
        AV1InterPrediction10BitMD(
            candidate_buffer_ptr->candidate_ptr->interp_filters,
#else
        AV1MDInterPrediction(
#endif
            picture_control_set_ptr,
            candidate_buffer_ptr->candidate_ptr->ref_frame_type,
            md_context_ptr,
            md_context_ptr->cu_ptr,
            &mv_unit,
            md_context_ptr->cu_origin_x,
            md_context_ptr->cu_origin_y,
            md_context_ptr->blk_geom->bwidth,
            md_context_ptr->blk_geom->bheight,
            ref_pic_list0,
            ref_pic_list1,
            candidate_buffer_ptr->prediction_ptr,
            md_context_ptr->blk_geom->origin_x,
            md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
            md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
            asm_type);
    } else {
        candidate_buffer_ptr->candidate_ptr->interp_filters = 0;
#if INTERPOLATION_SEARCH_LEVELS
        if (!md_context_ptr->skip_interpolation_search) {
#else
        if (picture_control_set_ptr->parent_pcs_ptr->interpolation_filter_search_mode > 0) {
#endif
            if (md_context_ptr->blk_geom->bwidth > 4 && md_context_ptr->blk_geom->bheight > 4)
                interpolation_filter_search(
                    picture_control_set_ptr,
                    candidate_buffer_ptr->predictionPtrTemp,
                    md_context_ptr,
                    candidate_buffer_ptr,
                    mv_unit,
                    ref_pic_list0,
                    ref_pic_list1,
                    asm_type,
                    &rd,
                    &rs,
                    &skip_txfm_sb,
                    &skip_sse_sb);
        }

        av1_inter_prediction(
            picture_control_set_ptr,
            candidate_buffer_ptr->candidate_ptr->interp_filters,
            md_context_ptr->cu_ptr,
            candidate_buffer_ptr->candidate_ptr->ref_frame_type,
            &mv_unit,
            md_context_ptr->cu_origin_x,
            md_context_ptr->cu_origin_y,
            md_context_ptr->blk_geom->bwidth,
            md_context_ptr->blk_geom->bheight,
            ref_pic_list0,
            ref_pic_list1,
            candidate_buffer_ptr->prediction_ptr,
            md_context_ptr->blk_geom->origin_x,
            md_context_ptr->blk_geom->origin_y,
#if CHROMA_BLIND
            md_context_ptr->chroma_level == CHROMA_MODE_0,
#endif
            asm_type);
    }

    return return_error;
}


EbErrorType inter_prediction_context_ctor(
    InterPredictionContext_t **inter_prediction_context,
    uint16_t                     max_cu_width,
    uint16_t                     max_cu_height)

{
    EbErrorType              return_error = EB_ErrorNone;
    InterPredictionContext_t *context_ptr;
    EB_MALLOC(InterPredictionContext_t*, context_ptr, sizeof(InterPredictionContext_t), EB_N_PTR);

    (*inter_prediction_context) = context_ptr;

    return_error = motion_compensation_prediction_context_ctor(
        &context_ptr->mcp_context,
        max_cu_width,
        max_cu_height);

    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    return EB_ErrorNone;
}


void RoundMvOnTheFly(
    int16_t *motionVector_x,
    int16_t *motionVector_y)
{
    *motionVector_x = (*motionVector_x + 2)&~0x03;
    *motionVector_y = (*motionVector_y + 2)&~0x03;

    return;
}

#if !CHROMA_BLIND
EbErrorType inter2_nx2_n_pu_prediction_avc(
    ModeDecisionContext_t                  *md_context_ptr,
    uint32_t                                  component_mask,
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
    EbAsm                                  asm_type)
{
    EbErrorType            return_error = EB_ErrorNone;
    EbPictureBufferDesc_t  *ref_pic_list0;
    EbPictureBufferDesc_t  *ref_pic_list1;
    EbReferenceObject_t    *referenceObject;

    uint32_t                  refList0PosX;
    uint32_t                  refList0PosY;
    uint32_t                  refList1PosX;
    uint32_t                  refList1PosY;
    const uint32_t pu_origin_x = md_context_ptr->cu_origin_x;
    const uint32_t pu_origin_y = md_context_ptr->cu_origin_y;
    const uint32_t pu_width = md_context_ptr->blk_geom->bwidth;
    const uint32_t pu_height = md_context_ptr->blk_geom->bheight;
    const uint32_t pu_index = md_context_ptr->pu_itr;

    //InterPredictionContext_t  *context_ptr = (InterPredictionContext_t*)prediction_context_ptr;
    InterPredictionContext_t *context_ptr = (InterPredictionContext_t*)(md_context_ptr->inter_prediction_context);
    uint32_t                  puOriginIndex;
    uint32_t                  puChromaOriginIndex;
    EncodeContext_t        *encode_context_ptr = ((SequenceControlSet_t*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr))->encode_context_ptr;

    SequenceControlSet_t* sequence_control_set_ptr = ((SequenceControlSet_t*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr));
    EbBool  is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    int16_t motionVector_x;
    int16_t motionVector_y;
    (void)asm_type;

    if (md_context_ptr->cu_use_ref_src_flag)

        is16bit = EB_FALSE;

    switch (candidate_buffer_ptr->candidate_ptr->prediction_direction[pu_index]) {
    case UNI_PRED_LIST_0:

        if (is16bit) {

            //----NOTE------- In 10Bit mode, the below Prediction code assumes no fractionnal whatsover

            puOriginIndex = ((pu_origin_y  & (63)) * 64) + (pu_origin_x & (63));
            puChromaOriginIndex = (((pu_origin_y  & (63)) * 32) + (pu_origin_x & (63))) >> 1;
            referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;
            ref_pic_list0 = (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit;


            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;


            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            //all the input quantites are in Qpel resolution
            refList0PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );


            refList0PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );



            estimate_uni_pred_interpolation_avc_lumaRef10Bit(
                ref_pic_list0,
                refList0PosX,
                refList0PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                candidate_buffer_ptr->sub_sampled_pred,
                candidate_buffer_ptr->sub_sampled_pred_chroma,
                asm_type);


        }
        else {
            puOriginIndex = ((pu_origin_y  & (63)) * 64) + (pu_origin_x & (63));
            puChromaOriginIndex = (((pu_origin_y  & (63)) * 32) + (pu_origin_x & (63))) >> 1;
            referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;
            ref_pic_list0 = (EbPictureBufferDesc_t*)referenceObject->referencePicture;

            if (md_context_ptr->cu_use_ref_src_flag)
                ref_pic_list0 = referenceObject->refDenSrcPicture;
            else
                ref_pic_list0 = referenceObject->referencePicture;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            //all the input quantites are in Qpel resolution
            refList0PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );


            refList0PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            uni_pred_i_free_ref8_bit(
                ref_pic_list0,
                refList0PosX,
                refList0PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                candidate_buffer_ptr->sub_sampled_pred,
                candidate_buffer_ptr->sub_sampled_pred_chroma,
                asm_type); // temp until mismatch fixed
        }

        break;

    case UNI_PRED_LIST_1:
        puOriginIndex = ((pu_origin_y  & (63)) * 64) + (pu_origin_x & (63));
        puChromaOriginIndex = (((pu_origin_y  & (63)) * 32) + (pu_origin_x & (63))) >> 1;
        referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr;

        if (is16bit) {

            ref_pic_list1 = (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;


            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }
            refList1PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList1PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            estimate_uni_pred_interpolation_avc_lumaRef10Bit(
                ref_pic_list1,
                refList1PosX,
                refList1PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                candidate_buffer_ptr->sub_sampled_pred,
                candidate_buffer_ptr->sub_sampled_pred_chroma,
                asm_type);



        }
        else {
            ref_pic_list1 = (EbPictureBufferDesc_t*)referenceObject->referencePicture;

            if (md_context_ptr->cu_use_ref_src_flag)
                ref_pic_list1 = referenceObject->refDenSrcPicture;
            else
                ref_pic_list1 = referenceObject->referencePicture;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;


            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            refList1PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList1PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            uni_pred_i_free_ref8_bit(
                ref_pic_list1,
                refList1PosX,
                refList1PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                candidate_buffer_ptr->sub_sampled_pred,
                candidate_buffer_ptr->sub_sampled_pred_chroma,
                asm_type); // temp until mismatch fixed

        }

        break;

    case BI_PRED:
        puOriginIndex = ((pu_origin_y  & (63)) * 64) + (pu_origin_x & (63));
        puChromaOriginIndex = (((pu_origin_y  & (63)) * 32) + (pu_origin_x & (63))) >> 1;

        // List0
        referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;

        if (is16bit) {

            ref_pic_list0 = (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;
            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            refList0PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList0PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            // List1
            referenceObject = (EbReferenceObject_t *)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr;
            ref_pic_list1 = (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;


            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            refList1PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList1PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            //to test getting the preLoad fucntion inside to do it just prior to interpolation, also seperate chroma from luma (do all luma, all Cb, all Cr).
            estimate_bi_pred_interpolation_avc_luma_ref10_bit(
                ref_pic_list0,
                ref_pic_list1,
                refList0PosX,
                refList0PosY,
                refList1PosX,
                refList1PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf1,
                context_ptr->mcp_context->avc_style_mcp_two_d_interpolation_first_pass_filter_result_buf,
                candidate_buffer_ptr->sub_sampled_pred,
                candidate_buffer_ptr->sub_sampled_pred_chroma,
                asm_type);



        }
        else {
            ref_pic_list0 = (EbPictureBufferDesc_t*)referenceObject->referencePicture;

            if (md_context_ptr->cu_use_ref_src_flag)
                ref_pic_list0 = referenceObject->refDenSrcPicture;
            else
                ref_pic_list0 = referenceObject->referencePicture;


            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            refList0PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList0PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            // List1
            referenceObject = (EbReferenceObject_t *)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr;
            ref_pic_list1 = (EbPictureBufferDesc_t*)referenceObject->referencePicture;

            if (md_context_ptr->cu_use_ref_src_flag)
                ref_pic_list1 = referenceObject->refDenSrcPicture;
            else
                ref_pic_list1 = referenceObject->referencePicture;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;


            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            refList1PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList1PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            bi_pred_i_free_ref8_bit(
                ref_pic_list0,
                ref_pic_list1,
                refList0PosX,
                refList0PosY,
                refList1PosX,
                refList1PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf1,
                context_ptr->mcp_context->avc_style_mcp_two_d_interpolation_first_pass_filter_result_buf,
                candidate_buffer_ptr->sub_sampled_pred,
                candidate_buffer_ptr->sub_sampled_pred_chroma,
                asm_type); // temp until mismatch fixed
        }

        break;

    default:
        CHECK_REPORT_ERROR_NC(
            encode_context_ptr->app_callback_ptr,
            EB_ENC_INTER_PRED_ERROR0);
        break;
    }

    return return_error;
}
#endif


/***************************************************
*  PreLoad Reference Block  for 16bit mode
***************************************************/
void UnPackReferenceLumaBlock(
    EbPictureBufferDesc_t *refFramePic,
    uint32_t                 posX,
    uint32_t                 posY,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    EbPictureBufferDesc_t *dst,
    EbBool                sub_pred,
    EbAsm                 asm_type)
{

    pu_width += 4;
    pu_height += 4;
    uint32_t inPosx = (posX >> 2) - 2;
    uint32_t inPosy = (posY >> 2) - 2;
    uint16_t *ptr16 = (uint16_t *)refFramePic->buffer_y + inPosx + inPosy * refFramePic->stride_y;

    extract8_bitdata_safe_sub(
        ptr16,
        refFramePic->stride_y << sub_pred,
        dst->buffer_y,
        dst->stride_y << sub_pred,
        pu_width,
        pu_height >> sub_pred,
        sub_pred,
        asm_type
    );
}

#if !CHROMA_BLIND
EbErrorType inter2_nx2_n_pu_prediction_avc_style(
    ModeDecisionContext_t                  *md_context_ptr,
    uint32_t                                  component_mask,
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
    EbAsm                                  asm_type)
{
    EbErrorType            return_error = EB_ErrorNone;
    EbPictureBufferDesc_t  *ref_pic_list0;
    EbPictureBufferDesc_t  *ref_pic_list1;
    EbReferenceObject_t    *referenceObject;

    uint32_t                  refList0PosX;
    uint32_t                  refList0PosY;
    uint32_t                  refList1PosX;
    uint32_t                  refList1PosY;
    const uint32_t pu_origin_x = md_context_ptr->cu_origin_x;
    const uint32_t pu_origin_y = md_context_ptr->cu_origin_y;
    const uint32_t pu_width = md_context_ptr->blk_geom->bwidth;
    const uint32_t pu_height = md_context_ptr->blk_geom->bheight;
    const uint32_t pu_index = md_context_ptr->pu_itr;


    //InterPredictionContext_t  *context_ptr = (InterPredictionContext_t*)prediction_context_ptr;
    InterPredictionContext_t *context_ptr = (InterPredictionContext_t*)(md_context_ptr->inter_prediction_context);
    uint32_t                  puOriginIndex;
    uint32_t                  puChromaOriginIndex;
    EncodeContext_t        *encode_context_ptr = ((SequenceControlSet_t*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr))->encode_context_ptr;

    SequenceControlSet_t* sequence_control_set_ptr = ((SequenceControlSet_t*)(picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr));
    EbBool  is16bit = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    int16_t motionVector_x;
    int16_t motionVector_y;
    (void)asm_type;

    EbBool   subSamplePredFlag = candidate_buffer_ptr->sub_sampled_pred;

    if (md_context_ptr->cu_use_ref_src_flag)

        is16bit = EB_FALSE;

    switch (candidate_buffer_ptr->candidate_ptr->prediction_direction[pu_index]) {
    case UNI_PRED_LIST_0:

        if (is16bit) {

            //----NOTE------- In 10Bit mode, the below Prediction code assumes no fractionnal whatsover

            puOriginIndex = ((pu_origin_y  & (63)) * 64) + (pu_origin_x & (63));
            puChromaOriginIndex = (((pu_origin_y  & (63)) * 32) + (pu_origin_x & (63))) >> 1;
            referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;
            ref_pic_list0 = (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit;


            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }
            //all the input quantites are in Qpel resolution
            refList0PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );


            refList0PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );


            if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {

                UnPackReferenceLumaBlock(
                    ref_pic_list0,
                    refList0PosX,
                    refList0PosY,
                    pu_width,
                    pu_height,
                    context_ptr->mcp_context->local_reference_block8_bitl0,
                    candidate_buffer_ptr->sub_sampled_pred,
                    asm_type);

                estimate_uni_pred_interpolation_unpacked_avc_style(
                    context_ptr->mcp_context->local_reference_block8_bitl0,
                    refList0PosX,
                    refList0PosY,
                    pu_width,
                    pu_height,
                    candidate_buffer_ptr->prediction_ptr,
                    puOriginIndex,
                    puChromaOriginIndex,
                    component_mask,
                    context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                    subSamplePredFlag,
                    asm_type); // temp until mismatch fixed
            }

            if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {


                estimate_uni_pred_interpolation_avc_chroma_ref10_bit(
                    ref_pic_list0,
                    refList0PosX,
                    refList0PosY,
                    pu_width,
                    pu_height,
                    candidate_buffer_ptr->prediction_ptr,
                    puOriginIndex,
                    puChromaOriginIndex,
                    component_mask,
                    context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                    candidate_buffer_ptr->sub_sampled_pred,
                    asm_type);
            }


        }
        else {
            puOriginIndex = ((pu_origin_y  & (63)) * 64) + (pu_origin_x & (63));
            puChromaOriginIndex = (((pu_origin_y  & (63)) * 32) + (pu_origin_x & (63))) >> 1;
            referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;
            ref_pic_list0 = (EbPictureBufferDesc_t*)referenceObject->referencePicture;

            if (md_context_ptr->cu_use_ref_src_flag)
                ref_pic_list0 = referenceObject->refDenSrcPicture;
            else
                ref_pic_list0 = referenceObject->referencePicture;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            //all the input quantites are in Qpel resolution
            refList0PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );


            refList0PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            estimate_uni_pred_interpolation_avc_luma(
                ref_pic_list0,
                refList0PosX,
                refList0PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                subSamplePredFlag,
                asm_type); // temp until mismatch fixed
        }

        break;

    case UNI_PRED_LIST_1:
        puOriginIndex = ((pu_origin_y  & (63)) * 64) + (pu_origin_x & (63));
        puChromaOriginIndex = (((pu_origin_y  & (63)) * 32) + (pu_origin_x & (63))) >> 1;
        referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr;

        if (is16bit) {

            ref_pic_list1 = (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }
            refList1PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList1PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {

                UnPackReferenceLumaBlock(
                    ref_pic_list1,
                    refList1PosX,
                    refList1PosY,
                    pu_width,
                    pu_height,
                    context_ptr->mcp_context->local_reference_block8_bitl1,
                    candidate_buffer_ptr->sub_sampled_pred,
                    asm_type);

                estimate_uni_pred_interpolation_unpacked_avc_style(
                    context_ptr->mcp_context->local_reference_block8_bitl1,
                    refList1PosX,
                    refList1PosY,
                    pu_width,
                    pu_height,
                    candidate_buffer_ptr->prediction_ptr,
                    puOriginIndex,
                    puChromaOriginIndex,
                    component_mask,
                    context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf1,
                    subSamplePredFlag,
                    asm_type);
            }

            if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {

                estimate_uni_pred_interpolation_avc_chroma_ref10_bit(
                    ref_pic_list1,
                    refList1PosX,
                    refList1PosY,
                    pu_width,
                    pu_height,
                    candidate_buffer_ptr->prediction_ptr,
                    puOriginIndex,
                    puChromaOriginIndex,
                    component_mask,
                    context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                    candidate_buffer_ptr->sub_sampled_pred,
                    asm_type);
            }



        }
        else {
            ref_pic_list1 = (EbPictureBufferDesc_t*)referenceObject->referencePicture;

            if (md_context_ptr->cu_use_ref_src_flag)
                ref_pic_list1 = referenceObject->refDenSrcPicture;
            else
                ref_pic_list1 = referenceObject->referencePicture;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }


            refList1PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList1PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );


            estimate_uni_pred_interpolation_avc_luma(
                ref_pic_list1,
                refList1PosX,
                refList1PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                subSamplePredFlag,
                asm_type); // temp until mismatch fixed
        }

        break;

    case BI_PRED:
        puOriginIndex = ((pu_origin_y  & (63)) * 64) + (pu_origin_x & (63));
        puChromaOriginIndex = (((pu_origin_y  & (63)) * 32) + (pu_origin_x & (63))) >> 1;

        // List0
        referenceObject = (EbReferenceObject_t*)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_0]->object_ptr;

        if (is16bit) {

            ref_pic_list0 = (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            refList0PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList0PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            // List1
            referenceObject = (EbReferenceObject_t *)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr;
            ref_pic_list1 = (EbPictureBufferDesc_t*)referenceObject->referencePicture16bit;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;


            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }

            refList1PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList1PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            if (component_mask & PICTURE_BUFFER_DESC_LUMA_MASK) {

                UnPackReferenceLumaBlock(
                    ref_pic_list0,
                    refList0PosX,
                    refList0PosY,
                    pu_width,
                    pu_height,
                    context_ptr->mcp_context->local_reference_block8_bitl0,
                    candidate_buffer_ptr->sub_sampled_pred,
                    asm_type);
                UnPackReferenceLumaBlock(
                    ref_pic_list1,
                    refList1PosX,
                    refList1PosY,
                    pu_width,
                    pu_height,
                    context_ptr->mcp_context->local_reference_block8_bitl1,
                    candidate_buffer_ptr->sub_sampled_pred,
                    asm_type);

                estimate_bi_pred_interpolation_unpacked_avc_style(
                    context_ptr->mcp_context->local_reference_block8_bitl0,
                    context_ptr->mcp_context->local_reference_block8_bitl1,
                    refList0PosX,
                    refList0PosY,
                    refList1PosX,
                    refList1PosY,
                    pu_width,
                    pu_height,
                    candidate_buffer_ptr->prediction_ptr,
                    puOriginIndex,
                    puChromaOriginIndex,
                    component_mask,
                    context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                    context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf1,
                    context_ptr->mcp_context->avc_style_mcp_two_d_interpolation_first_pass_filter_result_buf,
                    subSamplePredFlag,
                    asm_type); // temp until mismatch fixed
            }
            if (component_mask & PICTURE_BUFFER_DESC_CHROMA_MASK) {
                //to test getting the preLoad fucntion inside to do it just prior to interpolation, also seperate chroma from luma (do all luma, all Cb, all Cr).
                estimate_bi_pred_interpolation_avc_chroma_ref10_bit(
                    ref_pic_list0,
                    ref_pic_list1,
                    refList0PosX,
                    refList0PosY,
                    refList1PosX,
                    refList1PosY,
                    pu_width,
                    pu_height,
                    candidate_buffer_ptr->prediction_ptr,
                    puOriginIndex,
                    puChromaOriginIndex,
                    component_mask,
                    context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                    context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf1,
                    context_ptr->mcp_context->avc_style_mcp_two_d_interpolation_first_pass_filter_result_buf,
                    candidate_buffer_ptr->sub_sampled_pred,
                    asm_type);
            }

        }
        else {
            ref_pic_list0 = (EbPictureBufferDesc_t*)referenceObject->referencePicture;

            if (md_context_ptr->cu_use_ref_src_flag)
                ref_pic_list0 = referenceObject->refDenSrcPicture;
            else
                ref_pic_list0 = referenceObject->referencePicture;

            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L0;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L0;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }
            refList0PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList0PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list0->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );



            // List1
            referenceObject = (EbReferenceObject_t *)picture_control_set_ptr->ref_pic_ptr_array[REF_LIST_1]->object_ptr;
            ref_pic_list1 = (EbPictureBufferDesc_t*)referenceObject->referencePicture;

            if (md_context_ptr->cu_use_ref_src_flag)
                ref_pic_list1 = referenceObject->refDenSrcPicture;
            else
                ref_pic_list1 = referenceObject->referencePicture;


            motionVector_x = candidate_buffer_ptr->candidate_ptr->motionVector_x_L1;
            motionVector_y = candidate_buffer_ptr->candidate_ptr->motionVector_y_L1;

            motionVector_x = motionVector_x >> 1;
            motionVector_y = motionVector_y >> 1;

            if (md_context_ptr->round_mv_to_integer) {
                RoundMvOnTheFly(
                    &motionVector_x,
                    &motionVector_y);
            }


            refList1PosX = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->width) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_x << 2) + REFPADD_QPEL + motionVector_x))
            );
            refList1PosY = (uint32_t)CLIP3((int32_t)MVBOUNDLOW,
                (int32_t)(((ref_pic_list1->height) << 2) + MVBOUNDHIGH),
                (int32_t)(((pu_origin_y << 2) + REFPADD_QPEL + motionVector_y))
            );

            estimate_bi_pred_interpolation_avc_luma(
                ref_pic_list0,
                ref_pic_list1,
                refList0PosX,
                refList0PosY,
                refList1PosX,
                refList1PosY,
                pu_width,
                pu_height,
                candidate_buffer_ptr->prediction_ptr,
                puOriginIndex,
                puChromaOriginIndex,
                component_mask,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf0,
                context_ptr->mcp_context->avc_style_mcp_intermediate_result_buf1,
                context_ptr->mcp_context->avc_style_mcp_two_d_interpolation_first_pass_filter_result_buf,
                subSamplePredFlag,
                asm_type); // temp until mismatch fixed
        }

        break;

    default:
        CHECK_REPORT_ERROR_NC(
            encode_context_ptr->app_callback_ptr,
            EB_ENC_INTER_PRED_ERROR0);
        break;
    }

    return return_error;
}
#endif
/** choose_mvp_idx_v2 function is used to choose the best AMVP candidate.
    @param *candidate_ptr(output)
        candidate_ptr points to the prediction result.
    @param cu_ptr(input)
        pointer to the CU where the target PU belongs to.
    @param *pu_index(input)
        the index of the PU inside a CU
    @param ref0AMVPCandArray(input)
    @param ref0_num_available_amvp_cand(input)
    @param ref1AMVPCandArray(input)
    @param ref1NumAvailableAMVPCand(input)
 */
EbErrorType choose_mvp_idx_v2(
    ModeDecisionCandidate_t  *candidate_ptr,
    uint32_t                    cu_origin_x,
    uint32_t                    cu_origin_y,
    uint32_t                    pu_index,
    uint32_t                    tbSize,
    int16_t                   *ref0_amvp_cand_array_x,
    int16_t                   *ref0_amvp_cand_array_y,
    uint32_t                    ref0_num_available_amvp_cand,
    int16_t                   *ref1_amvp_cand_array_x,
    int16_t                   *ref1_amvp_cand_array_y,
    uint32_t                    ref1NumAvailableAMVPCand,
    PictureControlSet_t      *picture_control_set_ptr)
{
    EbErrorType  return_error = EB_ErrorNone;
    uint8_t         mvpRef0Idx;
    uint8_t         mvpRef1Idx;


    uint32_t        picture_width = ((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr)->luma_width;
    uint32_t        picture_height = ((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr)->luma_height;

    uint32_t   mvd0, mvd1;

    switch (candidate_ptr->prediction_direction[pu_index]) {
    case UNI_PRED_LIST_0:
        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motionVector_x_L0,
            &candidate_ptr->motionVector_y_L0,
            picture_width,
            picture_height,
            tbSize);

        // Choose the AMVP candidate
        switch (ref0_num_available_amvp_cand) {
        case 0:
        case 1:
            //mvpRef0Idx = 0;
            candidate_ptr->motion_vector_pred_idx[REF_LIST_0] = 0;
            candidate_ptr->motion_vector_pred_x[REF_LIST_0] = ref0_amvp_cand_array_x[0];
            candidate_ptr->motion_vector_pred_y[REF_LIST_0] = ref0_amvp_cand_array_y[0];
            break;
        case 2:

            mvd0 = EB_ABS_DIFF(ref0_amvp_cand_array_x[0], candidate_ptr->motionVector_x_L0) +
                EB_ABS_DIFF(ref0_amvp_cand_array_y[0], candidate_ptr->motionVector_y_L0);

            mvd1 = EB_ABS_DIFF(ref0_amvp_cand_array_x[1], candidate_ptr->motionVector_x_L0) +
                EB_ABS_DIFF(ref0_amvp_cand_array_y[1], candidate_ptr->motionVector_y_L0);

            mvpRef0Idx = ((mvd0) <= (mvd1)) ? 0 : 1;

            candidate_ptr->motion_vector_pred_idx[REF_LIST_0] = mvpRef0Idx;
            candidate_ptr->motion_vector_pred_x[REF_LIST_0] = ref0_amvp_cand_array_x[mvpRef0Idx];
            candidate_ptr->motion_vector_pred_y[REF_LIST_0] = ref0_amvp_cand_array_y[mvpRef0Idx];
            break;
        default:
            break;
        }

        break;

    case UNI_PRED_LIST_1:

        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motionVector_x_L1,
            &candidate_ptr->motionVector_y_L1,
            picture_width,
            picture_height,
            tbSize);

        // Choose the AMVP candidate
        switch (ref1NumAvailableAMVPCand) {
        case 0:
        case 1:
            //mvpRef1Idx = 0;
            candidate_ptr->motion_vector_pred_idx[REF_LIST_1] = 0;
            candidate_ptr->motion_vector_pred_x[REF_LIST_1] = ref1_amvp_cand_array_x[0];
            candidate_ptr->motion_vector_pred_y[REF_LIST_1] = ref1_amvp_cand_array_y[0];
            break;
        case 2:

            mvd0 = EB_ABS_DIFF(ref1_amvp_cand_array_x[0], candidate_ptr->motionVector_x_L1) +
                EB_ABS_DIFF(ref1_amvp_cand_array_y[0], candidate_ptr->motionVector_y_L1);

            mvd1 = EB_ABS_DIFF(ref1_amvp_cand_array_x[1], candidate_ptr->motionVector_x_L1) +
                EB_ABS_DIFF(ref1_amvp_cand_array_y[1], candidate_ptr->motionVector_y_L1);

            mvpRef1Idx = ((mvd0) <= (mvd1)) ? 0 : 1;



            candidate_ptr->motion_vector_pred_idx[REF_LIST_1] = mvpRef1Idx;
            candidate_ptr->motion_vector_pred_x[REF_LIST_1] = ref1_amvp_cand_array_x[mvpRef1Idx];
            candidate_ptr->motion_vector_pred_y[REF_LIST_1] = ref1_amvp_cand_array_y[mvpRef1Idx];
            break;
        default:
            break;
        }

        // MVP in ref_pic_list0
        //mvpRef0Idx = 0;
        //candidate_ptr->motion_vector_pred_idx[REF_LIST_0][pu_index] = mvpRef0Idx;
        //candidate_ptr->motion_vector_pred_x[REF_LIST_0][pu_index]  = 0;
        //candidate_ptr->motion_vector_pred_y[REF_LIST_0][pu_index]  = 0;

        break;

    case BI_PRED:

        // Choose the MVP in list0
        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motionVector_x_L0,
            &candidate_ptr->motionVector_y_L0,
            picture_width,
            picture_height,
            tbSize);

        // Choose the AMVP candidate
        switch (ref0_num_available_amvp_cand) {
        case 0:
        case 1:
            //mvpRef0Idx = 0;
            candidate_ptr->motion_vector_pred_idx[REF_LIST_0] = 0;
            candidate_ptr->motion_vector_pred_x[REF_LIST_0] = ref0_amvp_cand_array_x[0];
            candidate_ptr->motion_vector_pred_y[REF_LIST_0] = ref0_amvp_cand_array_y[0];
            break;
        case 2:

            mvd0 = EB_ABS_DIFF(ref0_amvp_cand_array_x[0], candidate_ptr->motionVector_x_L0) +
                EB_ABS_DIFF(ref0_amvp_cand_array_y[0], candidate_ptr->motionVector_y_L0);

            mvd1 = EB_ABS_DIFF(ref0_amvp_cand_array_x[1], candidate_ptr->motionVector_x_L0) +
                EB_ABS_DIFF(ref0_amvp_cand_array_y[1], candidate_ptr->motionVector_y_L0);

            mvpRef0Idx = ((mvd0) <= (mvd1)) ? 0 : 1;



            candidate_ptr->motion_vector_pred_idx[REF_LIST_0] = mvpRef0Idx;
            candidate_ptr->motion_vector_pred_x[REF_LIST_0] = ref0_amvp_cand_array_x[mvpRef0Idx];
            candidate_ptr->motion_vector_pred_y[REF_LIST_0] = ref0_amvp_cand_array_y[mvpRef0Idx];
            break;
        default:
            break;
        }

        // Choose the MVP in list1
        // Clip the input MV
        clip_mv(
            cu_origin_x,
            cu_origin_y,
            &candidate_ptr->motionVector_x_L1,
            &candidate_ptr->motionVector_y_L1,
            picture_width,
            picture_height,
            tbSize);

        // Choose the AMVP candidate
        switch (ref1NumAvailableAMVPCand) {
        case 0:
        case 1:
            //mvpRef1Idx = 0;
            candidate_ptr->motion_vector_pred_idx[REF_LIST_1] = 0;
            candidate_ptr->motion_vector_pred_x[REF_LIST_1] = ref1_amvp_cand_array_x[0];
            candidate_ptr->motion_vector_pred_y[REF_LIST_1] = ref1_amvp_cand_array_y[0];
            break;
        case 2:

            mvd0 = EB_ABS_DIFF(ref1_amvp_cand_array_x[0], candidate_ptr->motionVector_x_L1) +
                EB_ABS_DIFF(ref1_amvp_cand_array_y[0], candidate_ptr->motionVector_y_L1);

            mvd1 = EB_ABS_DIFF(ref1_amvp_cand_array_x[1], candidate_ptr->motionVector_x_L1) +
                EB_ABS_DIFF(ref1_amvp_cand_array_y[1], candidate_ptr->motionVector_y_L1);

            mvpRef1Idx = ((mvd0) <= (mvd1)) ? 0 : 1;

            candidate_ptr->motion_vector_pred_idx[REF_LIST_1] = mvpRef1Idx;
            candidate_ptr->motion_vector_pred_x[REF_LIST_1] = ref1_amvp_cand_array_x[mvpRef1Idx];
            candidate_ptr->motion_vector_pred_y[REF_LIST_1] = ref1_amvp_cand_array_y[mvpRef1Idx];
            break;
        default:
            break;
        }

        break;

    default:
        break;
    }

    return return_error;
}
