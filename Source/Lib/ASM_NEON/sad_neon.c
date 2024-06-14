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

#include "mem_neon.h"
#include "aom_dsp_rtcd.h"
#include "sum_neon.h"

void svt_aom_compute8x4_sad_kernel_dual_neon(uint8_t *restrict src_ptr, uint32_t src_stride, uint8_t *restrict ref_ptr,
                                             uint32_t ref_stride, uint32_t *restrict result_0,
                                             uint32_t *restrict result_1) {
    uint16x8_t sum_0 = vdupq_n_u16(0);
    uint16x8_t sum_1 = vdupq_n_u16(0);

    for (unsigned int i = 0; i < 4; ++i) {
        uint8x8_t s0 = vld1_u8(src_ptr);
        uint8x8_t r0 = vld1_u8(ref_ptr);
        sum_0        = vabal_u8(sum_0, s0, r0);

        uint8x8_t s1 = vld1_u8(src_ptr + 8);
        uint8x8_t r1 = vld1_u8(ref_ptr + 8);
        sum_1        = vabal_u8(sum_1, s1, r1);

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    }

    *result_0 = vaddlvq_u16(sum_0);
    *result_1 = vaddlvq_u16(sum_1);
}

void compute8x8_sad_kernel_dual_neon(uint8_t *restrict src_ptr, uint32_t src_stride, uint8_t *restrict ref_ptr,
                                     uint32_t ref_stride, uint32_t *restrict result_0, uint32_t *restrict result_1) {
    uint16x8_t sum_0 = vdupq_n_u16(0);
    uint16x8_t sum_1 = vdupq_n_u16(0);

    for (unsigned int i = 0; i < 8; ++i) {
        uint8x8_t s1 = vld1_u8(src_ptr);
        uint8x8_t r1 = vld1_u8(ref_ptr);
        sum_0        = vabal_u8(sum_0, s1, r1);

        uint8x8_t s2 = vld1_u8(src_ptr + 8);
        uint8x8_t r2 = vld1_u8(ref_ptr + 8);
        sum_1        = vabal_u8(sum_1, s2, r2);

        src_ptr += src_stride;
        ref_ptr += ref_stride;
    }

    *result_0 = vaddlvq_u16(sum_0);
    *result_1 = vaddlvq_u16(sum_1);
}

/*******************************************
Calculate SAD for 16x16 and its 8x8 sublcoks
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void svt_ext_sad_calculation_8x8_16x16_neon_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                   uint32_t *p_best_sad_8x8, uint32_t *p_best_sad_16x16,
                                                   uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16, uint32_t mv,
                                                   uint32_t *p_sad16x16, uint32_t *p_sad8x8, Bool sub_sad) {
    uint32_t sad16x16;

    if (sub_sad) {
        svt_aom_compute8x4_sad_kernel_dual_neon(src + 0 * src_stride + 0,
                                                2 * src_stride,
                                                ref + 0 * ref_stride + 0,
                                                2 * ref_stride,
                                                &p_sad8x8[0],
                                                &p_sad8x8[1]);
        p_sad8x8[0] = p_sad8x8[0] << 1;
        p_sad8x8[1] = p_sad8x8[1] << 1;

        svt_aom_compute8x4_sad_kernel_dual_neon(src + 8 * src_stride + 0,
                                                2 * src_stride,
                                                ref + 8 * ref_stride + 0,
                                                2 * ref_stride,
                                                &p_sad8x8[2],
                                                &p_sad8x8[3]);
        p_sad8x8[2] = p_sad8x8[2] << 1;
        p_sad8x8[3] = p_sad8x8[3] << 1;
    } else {
        compute8x8_sad_kernel_dual_neon(
            src + 0 * src_stride + 0, src_stride, ref + 0 * ref_stride + 0, ref_stride, &p_sad8x8[0], &p_sad8x8[1]);
        compute8x8_sad_kernel_dual_neon(
            src + 8 * src_stride + 0, src_stride, ref + 8 * ref_stride + 0, ref_stride, &p_sad8x8[2], &p_sad8x8[3]);
    }

    if (p_sad8x8[0] < p_best_sad_8x8[0]) {
        p_best_sad_8x8[0] = (uint32_t)p_sad8x8[0];
        p_best_mv8x8[0]   = mv;
    }

    if (p_sad8x8[1] < p_best_sad_8x8[1]) {
        p_best_sad_8x8[1] = (uint32_t)p_sad8x8[1];
        p_best_mv8x8[1]   = mv;
    }

    if (p_sad8x8[2] < p_best_sad_8x8[2]) {
        p_best_sad_8x8[2] = (uint32_t)p_sad8x8[2];
        p_best_mv8x8[2]   = mv;
    }

    if (p_sad8x8[3] < p_best_sad_8x8[3]) {
        p_best_sad_8x8[3] = (uint32_t)p_sad8x8[3];
        p_best_mv8x8[3]   = mv;
    }

    sad16x16 = p_sad8x8[0] + p_sad8x8[1] + p_sad8x8[2] + p_sad8x8[3];
    if (sad16x16 < p_best_sad_16x16[0]) {
        p_best_sad_16x16[0] = (uint32_t)sad16x16;
        p_best_mv16x16[0]   = mv;
    }

    *p_sad16x16 = (uint32_t)sad16x16;
}

/*******************************************
 * svt_ext_eight_sad_calculation_8x8_16x16_neon
 *******************************************/
static void svt_ext_eight_sad_calculation_8x8_16x16_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                                         uint32_t ref_stride, uint32_t mv, uint32_t start_16x16_pos,
                                                         uint32_t *p_best_sad_8x8, uint32_t *p_best_sad_16x16,
                                                         uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16,
                                                         uint32_t p_eight_sad16x16[16][8],
                                                         uint32_t p_eight_sad8x8[64][8], Bool sub_sad) {
    const uint32_t start_8x8_pos = 4 * start_16x16_pos;
    int16_t        x_mv, y_mv;

    (void)p_eight_sad8x8;

    p_best_sad_8x8 += start_8x8_pos;
    p_best_mv8x8 += start_8x8_pos;
    p_best_sad_16x16 += start_16x16_pos;
    p_best_mv16x16 += start_16x16_pos;
    if (sub_sad) {
        uint32_t src_stride_sub = (src_stride << 1);
        uint32_t ref_stride_sub = (ref_stride << 1);
        for (int search_index = 0; search_index < 8; search_index++) {
            uint32_t sad8x8_0;
            uint32_t sad8x8_1;
            svt_aom_compute8x4_sad_kernel_dual_neon(
                src, src_stride_sub, ref + search_index, ref_stride_sub, &sad8x8_0, &sad8x8_1);

            sad8x8_0 = sad8x8_0 << 1;
            sad8x8_1 = sad8x8_1 << 1;

            if (sad8x8_0 < p_best_sad_8x8[0]) {
                p_best_sad_8x8[0] = (uint32_t)sad8x8_0;
                x_mv              = _MVXT(mv) + (int16_t)search_index;
                y_mv              = _MVYT(mv);
                p_best_mv8x8[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }

            if (sad8x8_1 < p_best_sad_8x8[1]) {
                p_best_sad_8x8[1] = (uint32_t)sad8x8_1;
                x_mv              = _MVXT(mv) + (int16_t)search_index;
                y_mv              = _MVYT(mv);
                p_best_mv8x8[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }

            uint32_t sad8x8_2;
            uint32_t sad8x8_3;
            svt_aom_compute8x4_sad_kernel_dual_neon(src + (src_stride << 3),
                                                    src_stride_sub,
                                                    ref + (ref_stride << 3) + search_index,
                                                    ref_stride_sub,
                                                    &sad8x8_2,
                                                    &sad8x8_3);
            sad8x8_2 = sad8x8_2 << 1;
            sad8x8_3 = sad8x8_3 << 1;

            if (sad8x8_2 < p_best_sad_8x8[2]) {
                p_best_sad_8x8[2] = (uint32_t)sad8x8_2;
                x_mv              = _MVXT(mv) + (int16_t)search_index;
                y_mv              = _MVYT(mv);
                p_best_mv8x8[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }

            if (sad8x8_3 < p_best_sad_8x8[3]) {
                p_best_sad_8x8[3] = (uint32_t)sad8x8_3;
                x_mv              = _MVXT(mv) + (int16_t)search_index;
                y_mv              = _MVYT(mv);
                p_best_mv8x8[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }
            uint32_t sad16x16 = p_eight_sad16x16[start_16x16_pos][search_index] = sad8x8_0 + sad8x8_1 + sad8x8_2 +
                sad8x8_3;
            if (sad16x16 < p_best_sad_16x16[0]) {
                p_best_sad_16x16[0] = (uint32_t)sad16x16;
                x_mv                = _MVXT(mv) + (int16_t)search_index;
                y_mv                = _MVYT(mv);
                p_best_mv16x16[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }
        }
    } else {
        for (int search_index = 0; search_index < 8; search_index++) {
            uint32_t sad8x8_0, sad8x8_1, sad8x8_2, sad8x8_3;
            compute8x8_sad_kernel_dual_neon(src, src_stride, ref + search_index, ref_stride, &sad8x8_0, &sad8x8_1);
            compute8x8_sad_kernel_dual_neon(src + (src_stride << 3),
                                            src_stride,
                                            ref + (ref_stride << 3) + search_index,
                                            ref_stride,
                                            &sad8x8_2,
                                            &sad8x8_3);

            if (sad8x8_0 < p_best_sad_8x8[0]) {
                p_best_sad_8x8[0] = (uint32_t)sad8x8_0;
                x_mv              = _MVXT(mv) + (int16_t)search_index;
                y_mv              = _MVYT(mv);
                p_best_mv8x8[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }

            if (sad8x8_1 < p_best_sad_8x8[1]) {
                p_best_sad_8x8[1] = (uint32_t)sad8x8_1;
                x_mv              = _MVXT(mv) + (int16_t)search_index;
                y_mv              = _MVYT(mv);
                p_best_mv8x8[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }

            if (sad8x8_2 < p_best_sad_8x8[2]) {
                p_best_sad_8x8[2] = (uint32_t)sad8x8_2;
                x_mv              = _MVXT(mv) + (int16_t)search_index;
                y_mv              = _MVYT(mv);
                p_best_mv8x8[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }

            if (sad8x8_3 < p_best_sad_8x8[3]) {
                p_best_sad_8x8[3] = (uint32_t)sad8x8_3;
                x_mv              = _MVXT(mv) + (int16_t)search_index;
                y_mv              = _MVYT(mv);
                p_best_mv8x8[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }

            uint32_t sad16x16 = p_eight_sad16x16[start_16x16_pos][search_index] = sad8x8_0 + sad8x8_1 + sad8x8_2 +
                sad8x8_3;
            if (sad16x16 < p_best_sad_16x16[0]) {
                p_best_sad_16x16[0] = (uint32_t)sad16x16;
                x_mv                = _MVXT(mv) + (int16_t)search_index;
                y_mv                = _MVYT(mv);
                p_best_mv16x16[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            }
        }
    }
}

void svt_ext_all_sad_calculation_8x8_16x16_neon(uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
                                                uint32_t mv, uint32_t *p_best_sad_8x8, uint32_t *p_best_sad_16x16,
                                                uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16,
                                                uint32_t p_eight_sad16x16[16][8], uint32_t p_eight_sad8x8[64][8],
                                                Bool sub_sad) {
    static const char offsets[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

    //---- 16x16 : 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            const uint32_t block_index           = 16 * y * src_stride + 16 * x;
            const uint32_t search_position_index = 16 * y * ref_stride + 16 * x;
            svt_ext_eight_sad_calculation_8x8_16x16_neon(src + block_index,
                                                         src_stride,
                                                         ref + search_position_index,
                                                         ref_stride,
                                                         mv,
                                                         offsets[4 * y + x],
                                                         p_best_sad_8x8,
                                                         p_best_sad_16x16,
                                                         p_best_mv8x8,
                                                         p_best_mv16x16,
                                                         p_eight_sad16x16,
                                                         p_eight_sad8x8,
                                                         sub_sad);
        }
    }
}

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
