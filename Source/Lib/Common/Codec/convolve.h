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

#ifndef AV1_COMMON_AV1_CONVOLVE_H_
#define AV1_COMMON_AV1_CONVOLVE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "EbDefinitions.h"

    //--filtyer.h
#define MAX_FILTER_TAP 8

    // With CONFIG_DUAL_FILTER, pack two InterpFilter's into a uint32_t: since
    // there are at most 10 filters, we can use 16 bits for each and have more than
    // enough space. This reduces argument passing and unifies the operation of
    // setting a (pair of) filters.
    //
    // Without CONFIG_DUAL_FILTER,
    typedef uint32_t InterpFilters;
    static INLINE InterpFilter av1_extract_interp_filter(InterpFilters filters,
        int32_t x_filter) {
        return (InterpFilter)((filters >> (x_filter ? 16 : 0)) & 0xffff);
    }

    static INLINE InterpFilters av1_make_interp_filters(InterpFilter y_filter,
        InterpFilter x_filter) {
        uint16_t y16 = y_filter & 0xffff;
        uint16_t x16 = x_filter & 0xffff;
        return y16 | ((uint32_t)x16 << 16);
    }

    static INLINE InterpFilters av1_broadcast_interp_filter(InterpFilter filter) {
        return av1_make_interp_filters(filter, filter);
    }

    static INLINE InterpFilter av1_unswitchable_filter(InterpFilter filter) {
        return filter == SWITCHABLE ? EIGHTTAP_REGULAR : filter;
    }

#define LOG_SWITCHABLE_FILTERS \
  2 /* (1 << LOG_SWITCHABLE_FILTERS) > SWITCHABLE_FILTERS */

#define MAX_SUBPEL_TAPS 12
#define SWITCHABLE_FILTER_CONTEXTS ((SWITCHABLE_FILTERS + 1) * 4)
#define INTER_FILTER_COMP_OFFSET (SWITCHABLE_FILTERS + 1)
#define INTER_FILTER_DIR_OFFSET ((SWITCHABLE_FILTERS + 1) * 2)
    /*
        typedef struct InterpFilterParams {
            const int16_t *filter_ptr;
            uint16_t taps;
            uint16_t subpel_shifts;
            InterpFilter interp_filter;
        } InterpFilterParams;*/
        //--filter.h

    //typedef uint16_t ConvBufType;
    //typedef struct ConvolveParams {
    //  int32_t ref;
    //  int32_t do_average;
    //  ConvBufType *dst;
    //  int32_t dst_stride;
    //  int32_t round_0;
    //  int32_t round_1;
    //  int32_t plane;
    //  int32_t is_compound;
    //  int32_t use_jnt_comp_avg;
    //  int32_t fwd_offset;
    //  int32_t bck_offset;
    //} ConvolveParams;

#define ROUND0_BITS 3
#define COMPOUND_ROUND1_BITS 7
#define WIENER_ROUND0_BITS 3

#define WIENER_CLAMP_LIMIT(r0, bd) (1 << ((bd) + 1 + FILTER_BITS - r0))

    typedef void(*aom_convolve_fn_t)(const uint8_t *src, int32_t src_stride,
        uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h,
        InterpFilterParams *filter_params_x,
        InterpFilterParams *filter_params_y,
        const int32_t subpel_x_q4, const int32_t subpel_y_q4,
        ConvolveParams *conv_params);

    typedef void(*aom_highbd_convolve_fn_t)(
        const uint16_t *src, int32_t src_stride, uint16_t *dst, int32_t dst_stride, int32_t w,
        int32_t h, const InterpFilterParams *filter_params_x,
        const InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
        const int32_t subpel_y_q4, ConvolveParams *conv_params, int32_t bd);

    //static INLINE void av1_get_convolve_filter_params(InterpFilters interp_filters,
    //                                                  InterpFilterParams *params_x,
    //                                                  InterpFilterParams *params_y,
    //                                                  int32_t w, int32_t h) {
    //  InterpFilter filter_x = av1_extract_interp_filter(interp_filters, 1);
    //  InterpFilter filter_y = av1_extract_interp_filter(interp_filters, 0);
    //  *params_x = av1_get_interp_filter_params_with_block_size(filter_x, w);
    //  *params_y = av1_get_interp_filter_params_with_block_size(filter_y, h);
    //}

    struct AV1Common;
    struct scale_factors;

    static INLINE ConvolveParams get_conv_params_no_round(int32_t ref, int32_t do_average,
        int32_t plane,
        ConvBufType *dst,
        int32_t dst_stride,
        int32_t is_compound, int32_t bd) {
        (void)plane;
        (void)ref;
        ConvolveParams conv_params;
        // conv_params.ref = ref;
        conv_params.do_average = do_average;
        assert(IMPLIES(do_average, is_compound));
        conv_params.is_compound = is_compound;
        conv_params.round_0 = ROUND0_BITS;
        conv_params.round_1 = is_compound ? COMPOUND_ROUND1_BITS
            : 2 * FILTER_BITS - conv_params.round_0;
        const int32_t intbufrange = bd + FILTER_BITS - conv_params.round_0 + 2;
        ASSERT(IMPLIES(bd < 12, intbufrange <= 16));
        if (intbufrange > 16) {
            conv_params.round_0 += intbufrange - 16;
            if (!is_compound) conv_params.round_1 -= intbufrange - 16;
        }
        // TODO(yunqing): The following dst should only be valid while
        // is_compound = 1;
        conv_params.dst = dst;
        conv_params.dst_stride = dst_stride;
        // conv_params.plane = plane;
        conv_params.use_jnt_comp_avg = 0;

        return conv_params;
    }

    static INLINE ConvolveParams get_conv_params(int32_t ref, int32_t do_average, int32_t plane,
        int32_t bd) {
        return get_conv_params_no_round(ref, do_average, plane, NULL, 0, 0, bd);
    }

    static INLINE ConvolveParams get_conv_params_wiener(int32_t bd) {
        ConvolveParams conv_params;
        (void)bd;
        conv_params.ref = 0;
        conv_params.do_average = 0;
        conv_params.is_compound = 0;
        conv_params.round_0 = WIENER_ROUND0_BITS;
        conv_params.round_1 = 2 * FILTER_BITS - conv_params.round_0;
        const int32_t intbufrange = bd + FILTER_BITS - conv_params.round_0 + 2;
        ASSERT(IMPLIES(bd < 12, intbufrange <= 16));
        if (intbufrange > 16) {
            conv_params.round_0 += intbufrange - 16;
            conv_params.round_1 -= intbufrange - 16;
        }
        conv_params.dst = NULL;
        conv_params.dst_stride = 0;
        conv_params.plane = 0;
        return conv_params;
    }

    void av1_highbd_convolve_2d_facade(const uint8_t *src8, int32_t src_stride,
        uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h,
        InterpFilters interp_filters,
        const int32_t subpel_x_q4, int32_t x_step_q4,
        const int32_t subpel_y_q4, int32_t y_step_q4,
        int32_t scaled, ConvolveParams *conv_params,
        const struct scale_factors *sf, int32_t bd);

    static INLINE const int16_t *av1_get_interp_filter_subpel_kernel(
        const InterpFilterParams filter_params, const int32_t subpel) {
        return filter_params.filter_ptr + filter_params.taps * subpel;
    }

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_AV1_CONVOLVE_H_
