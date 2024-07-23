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

#include <arm_neon.h>

#include "aom_dsp_rtcd.h"
#include "mem_neon.h"
#include "sum_neon.h"

int svt_aom_satd_neon(const TranLow *coeff, int length) {
    const int32x4_t zero  = vdupq_n_s32(0);
    int32x4_t       accum = zero;
    do {
        const int32x4_t src0  = vld1q_s32(&coeff[0]);
        const int32x4_t src8  = vld1q_s32(&coeff[4]);
        const int32x4_t src16 = vld1q_s32(&coeff[8]);
        const int32x4_t src24 = vld1q_s32(&coeff[12]);
        accum                 = vabaq_s32(accum, src0, zero);
        accum                 = vabaq_s32(accum, src8, zero);
        accum                 = vabaq_s32(accum, src16, zero);
        accum                 = vabaq_s32(accum, src24, zero);
        length -= 16;
        coeff += 16;
    } while (length != 0);

    // satd: 26 bits, dynamic range [-32640 * 1024, 32640 * 1024]
    return horizontal_add_s32x4(accum);
}

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
