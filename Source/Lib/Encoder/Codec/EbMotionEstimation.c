/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */
/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include <stdio.h>
#include <inttypes.h>

#include "aom_dsp_rtcd.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbMotionEstimation.h"
#include "EbUtility.h"

#include "EbComputeSAD.h"
#include "EbReferenceObject.h"
#include "EbAvcStyleMcp.h"

#include "EbEncIntraPrediction.h"
#include "EbLambdaRateTables.h"

#include "EbLog.h"
#include "EbResize.h"

#define AVCCODEL
/********************************************
 * Constants
 ********************************************/

#define MAX_INTRA_IN_MD 9
#define REFERENCE_PIC_LIST_0 0
#define REFERENCE_PIC_LIST_1 1
#define SC_HME_TH_STILL 1000
#define SC_HME_TH_EASY  100
#define SC_SR_DENOM_STILL 16
#define SC_SR_DENOM_EASY 8
/*******************************************
 * Compute8x4SAD_Default
 *   Unoptimized 8x4 SAD
 *******************************************/
uint32_t compute8x4_sad_kernel_c(uint8_t *src, // input parameter, source samples Ptr
                                 uint32_t src_stride, // input parameter, source stride
                                 uint8_t *ref, // input parameter, reference samples Ptr
                                 uint32_t ref_stride) // input parameter, reference stride
{
    uint32_t row_number_in_blocks_8x4;
    uint32_t sad_block_8x4 = 0;

    for (row_number_in_blocks_8x4 = 0; row_number_in_blocks_8x4 < 4; ++row_number_in_blocks_8x4) {
        sad_block_8x4 += EB_ABS_DIFF(src[0x00], ref[0x00]);
        sad_block_8x4 += EB_ABS_DIFF(src[0x01], ref[0x01]);
        sad_block_8x4 += EB_ABS_DIFF(src[0x02], ref[0x02]);
        sad_block_8x4 += EB_ABS_DIFF(src[0x03], ref[0x03]);
        sad_block_8x4 += EB_ABS_DIFF(src[0x04], ref[0x04]);
        sad_block_8x4 += EB_ABS_DIFF(src[0x05], ref[0x05]);
        sad_block_8x4 += EB_ABS_DIFF(src[0x06], ref[0x06]);
        sad_block_8x4 += EB_ABS_DIFF(src[0x07], ref[0x07]);
        src += src_stride;
        ref += ref_stride;
    }

    return sad_block_8x4;
}
/*******************************************
 * Compute8x8SAD_Default
 *   Unoptimized 8x8 SAD
 *******************************************/
uint32_t compute8x8_sad_kernel_c(uint8_t *src, // input parameter, source samples Ptr
                                 uint32_t src_stride, // input parameter, source stride
                                 uint8_t *ref, // input parameter, reference samples Ptr
                                 uint32_t ref_stride) // input parameter, reference stride
{
    uint32_t row_number_in_blocks_8x8;
    uint32_t sad_block_8x8 = 0;

    for (row_number_in_blocks_8x8 = 0; row_number_in_blocks_8x8 < 8; ++row_number_in_blocks_8x8) {
        sad_block_8x8 += EB_ABS_DIFF(src[0x00], ref[0x00]);
        sad_block_8x8 += EB_ABS_DIFF(src[0x01], ref[0x01]);
        sad_block_8x8 += EB_ABS_DIFF(src[0x02], ref[0x02]);
        sad_block_8x8 += EB_ABS_DIFF(src[0x03], ref[0x03]);
        sad_block_8x8 += EB_ABS_DIFF(src[0x04], ref[0x04]);
        sad_block_8x8 += EB_ABS_DIFF(src[0x05], ref[0x05]);
        sad_block_8x8 += EB_ABS_DIFF(src[0x06], ref[0x06]);
        sad_block_8x8 += EB_ABS_DIFF(src[0x07], ref[0x07]);
        src += src_stride;
        ref += ref_stride;
    }

    return sad_block_8x8;
}

/*******************************************
Calculate SAD for 16x16 and its 8x8 sublcoks
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void ext_sad_calculation_8x8_16x16_c(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                     uint32_t ref_stride, uint32_t *p_best_sad_8x8,
                                     uint32_t *p_best_sad_16x16, uint32_t *p_best_mv8x8,
                                     uint32_t *p_best_mv16x16, uint32_t mv, uint32_t *p_sad16x16,
                                     uint32_t *p_sad8x8, EbBool sub_sad) {
    uint32_t sad16x16;

    if (sub_sad) {
        p_sad8x8[0] =
            (compute8x4_sad_kernel_c(
                src + 0 * src_stride + 0, 2 * src_stride, ref + 0 * ref_stride + 0, 2 * ref_stride))
            << 1;
        p_sad8x8[1] =
            (compute8x4_sad_kernel_c(
                src + 0 * src_stride + 8, 2 * src_stride, ref + 0 * ref_stride + 8, 2 * ref_stride))
            << 1;
        p_sad8x8[2] =
            (compute8x4_sad_kernel_c(
                src + 8 * src_stride + 0, 2 * src_stride, ref + 8 * ref_stride + 0, 2 * ref_stride))
            << 1;
        p_sad8x8[3] =
            (compute8x4_sad_kernel_c(
                src + 8 * src_stride + 8, 2 * src_stride, ref + 8 * ref_stride + 8, 2 * ref_stride))
            << 1;
    } else {
        p_sad8x8[0] = compute8x8_sad_kernel_c(
            src + 0 * src_stride + 0, src_stride, ref + 0 * ref_stride + 0, ref_stride);
        p_sad8x8[1] = compute8x8_sad_kernel_c(
            src + 0 * src_stride + 8, src_stride, ref + 0 * ref_stride + 8, ref_stride);
        p_sad8x8[2] = compute8x8_sad_kernel_c(
            src + 8 * src_stride + 0, src_stride, ref + 8 * ref_stride + 0, ref_stride);
        p_sad8x8[3] = compute8x8_sad_kernel_c(
            src + 8 * src_stride + 8, src_stride, ref + 8 * ref_stride + 8, ref_stride);
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
Calculate SAD for 32x32,64x64 from 16x16
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void ext_sad_calculation_32x32_64x64_c(uint32_t *p_sad16x16, uint32_t *p_best_sad_32x32,
                                       uint32_t *p_best_sad_64x64, uint32_t *p_best_mv32x32,
                                       uint32_t *p_best_mv64x64, uint32_t mv,
                                       uint32_t *p_sad32x32) {
    uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;

    p_sad32x32[0] = sad32x32_0 = p_sad16x16[0] + p_sad16x16[1] + p_sad16x16[2] + p_sad16x16[3];
    if (sad32x32_0 < p_best_sad_32x32[0]) {
        p_best_sad_32x32[0] = sad32x32_0;
        p_best_mv32x32[0]   = mv;
    }

    p_sad32x32[1] = sad32x32_1 = p_sad16x16[4] + p_sad16x16[5] + p_sad16x16[6] + p_sad16x16[7];
    if (sad32x32_1 < p_best_sad_32x32[1]) {
        p_best_sad_32x32[1] = sad32x32_1;
        p_best_mv32x32[1]   = mv;
    }

    p_sad32x32[2] = sad32x32_2 = p_sad16x16[8] + p_sad16x16[9] + p_sad16x16[10] + p_sad16x16[11];
    if (sad32x32_2 < p_best_sad_32x32[2]) {
        p_best_sad_32x32[2] = sad32x32_2;
        p_best_mv32x32[2]   = mv;
    }

    p_sad32x32[3] = sad32x32_3 = p_sad16x16[12] + p_sad16x16[13] + p_sad16x16[14] + p_sad16x16[15];
    if (sad32x32_3 < p_best_sad_32x32[3]) {
        p_best_sad_32x32[3] = sad32x32_3;
        p_best_mv32x32[3]   = mv;
    }
    sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
    if (sad64x64 < p_best_sad_64x64[0]) {
        p_best_sad_64x64[0] = sad64x64;
        p_best_mv64x64[0]   = mv;
    }
}

/*******************************************
 * GetEightHorizontalSearchPointResults_8x8_16x16_PU
 *******************************************/
void get_eight_horizontal_search_point_results_8x8_16x16_pu_c(
    uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride, uint32_t *p_best_sad_8x8,
    uint32_t *p_best_mv8x8, uint32_t *p_best_sad_16x16, uint32_t *p_best_mv16x16, uint32_t mv,
    uint16_t *p_sad16x16, EbBool sub_sad) {
    uint32_t x_search_index;
    int16_t  x_mv, y_mv;
    uint32_t sad8x8[4];
    uint16_t sad16x16;

    /*
    -------------------------------------   -----------------------------------
    | 8x8_00 | 8x8_01 | 8x8_04 | 8x8_05 |   8x8_16 | 8x8_17 | 8x8_20 | 8x8_21 |
    -------------------------------------   -----------------------------------
    | 8x8_02 | 8x8_03 | 8x8_06 | 8x8_07 |   8x8_18 | 8x8_19 | 8x8_22 | 8x8_23 |
    -----------------------   -----------   ----------------------   ----------
    | 8x8_08 | 8x8_09 | 8x8_12 | 8x8_13 |   8x8_24 | 8x8_25 | 8x8_29 | 8x8_29 |
    ----------------------    -----------   ---------------------    ----------
    | 8x8_10 | 8x8_11 | 8x8_14 | 8x8_15 |   8x8_26 | 8x8_27 | 8x8_30 | 8x8_31 |
    -------------------------------------   -----------------------------------

    -------------------------------------   -----------------------------------
    | 8x8_32 | 8x8_33 | 8x8_36 | 8x8_37 |   8x8_48 | 8x8_49 | 8x8_52 | 8x8_53 |
    -------------------------------------   -----------------------------------
    | 8x8_34 | 8x8_35 | 8x8_38 | 8x8_39 |   8x8_50 | 8x8_51 | 8x8_54 | 8x8_55 |
    -----------------------   -----------   ----------------------   ----------
    | 8x8_40 | 8x8_41 | 8x8_44 | 8x8_45 |   8x8_56 | 8x8_57 | 8x8_60 | 8x8_61 |
    ----------------------    -----------   ---------------------    ----------
    | 8x8_42 | 8x8_43 | 8x8_46 | 8x8_48 |   8x8_58 | 8x8_59 | 8x8_62 | 8x8_63 |
    -------------------------------------   -----------------------------------
    */

    /*
    ----------------------    ----------------------
    |  16x16_0  |  16x16_1  |  16x16_4  |  16x16_5  |
    ----------------------    ----------------------
    |  16x16_2  |  16x16_3  |  16x16_6  |  16x16_7  |
    -----------------------   -----------------------
    |  16x16_8  |  16x16_9  |  16x16_12 |  16x16_13 |
    ----------------------    ----------------------
    |  16x16_10 |  16x16_11 |  16x16_14 |  16x16_15 |
    -----------------------   -----------------------
    */

    for (x_search_index = 0; x_search_index < 8; x_search_index++) {
        if (sub_sad) {
            sad8x8[0] = compute8x4_sad_kernel_c(src + 0 * src_stride + 0,
                                                2 * src_stride,
                                                ref + 0 * ref_stride + 0 + x_search_index,
                                                2 * ref_stride)
                        << 1;
            sad8x8[1] = compute8x4_sad_kernel_c(src + 0 * src_stride + 8,
                                                2 * src_stride,
                                                ref + 0 * ref_stride + 8 + x_search_index,
                                                2 * ref_stride)
                        << 1;
            sad8x8[2] = compute8x4_sad_kernel_c(src + 8 * src_stride + 0,
                                                2 * src_stride,
                                                ref + 8 * ref_stride + 0 + x_search_index,
                                                2 * ref_stride)
                        << 1;
            sad8x8[3] = compute8x4_sad_kernel_c(src + 8 * src_stride + 8,
                                                2 * src_stride,
                                                ref + 8 * ref_stride + 8 + x_search_index,
                                                2 * ref_stride)
                        << 1;
        } else {
            sad8x8[0] = compute8x8_sad_kernel_c(src + 0 * src_stride + 0,
                                                src_stride,
                                                ref + 0 * ref_stride + 0 + x_search_index,
                                                ref_stride);
            sad8x8[1] = compute8x8_sad_kernel_c(src + 0 * src_stride + 8,
                                                src_stride,
                                                ref + 0 * ref_stride + 8 + x_search_index,
                                                ref_stride);
            sad8x8[2] = compute8x8_sad_kernel_c(src + 8 * src_stride + 0,
                                                src_stride,
                                                ref + 8 * ref_stride + 0 + x_search_index,
                                                ref_stride);
            sad8x8[3] = compute8x8_sad_kernel_c(src + 8 * src_stride + 8,
                                                src_stride,
                                                ref + 8 * ref_stride + 8 + x_search_index,
                                                ref_stride);
        }

        // 8x8_0
        if (sad8x8[0] < p_best_sad_8x8[0]) {
            p_best_sad_8x8[0] = sad8x8[0];
            x_mv              = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv              = _MVYT(mv);
            p_best_mv8x8[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 8x8_1
        if (sad8x8[1] < p_best_sad_8x8[1]) {
            p_best_sad_8x8[1] = sad8x8[1];
            x_mv              = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv              = _MVYT(mv);
            p_best_mv8x8[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 8x8_2
        if (sad8x8[2] < p_best_sad_8x8[2]) {
            p_best_sad_8x8[2] = sad8x8[2];
            x_mv              = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv              = _MVYT(mv);
            p_best_mv8x8[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 8x8_3
        if (sad8x8[3] < p_best_sad_8x8[3]) {
            p_best_sad_8x8[3] = sad8x8[3];
            x_mv              = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv              = _MVYT(mv);
            p_best_mv8x8[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 16x16
        sad16x16                   = (uint16_t)(sad8x8[0] + sad8x8[1] + sad8x8[2] + sad8x8[3]);
        p_sad16x16[x_search_index] = sad16x16; // store the intermediate 16x16 SAD for 32x32.
        if ((uint32_t)(sad16x16) < p_best_sad_16x16[0]) {
            p_best_sad_16x16[0] = sad16x16;
            x_mv                = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x16[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }
}

/*******************************************
Calculate SAD for 32x32,64x64 from 16x16
and check if there is improvement, if yes keep
the best SAD+MV
*******************************************/
void get_eight_horizontal_search_point_results_32x32_64x64_pu_c(
    uint16_t *p_sad16x16, uint32_t *p_best_sad_32x32, uint32_t *p_best_sad_64x64,
    uint32_t *p_best_mv32x32, uint32_t *p_best_mv64x64, uint32_t mv) {
    int16_t  x_mv, y_mv;
    uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;
    uint32_t x_search_index;

    /*--------------------
    |  32x32_0  |  32x32_1
    ----------------------
    |  32x32_2  |  32x32_3
    ----------------------*/

    /*  data ordering in p_sad16x16 buffer

    Search    Search            Search
    Point 0   Point 1           Point 7
    ---------------------------------------
    16x16_0    |    x    |    x    | ...... |    x    |
    ---------------------------------------
    16x16_1    |    x    |    x    | ...... |    x    |

    16x16_n    |    x    |    x    | ...... |    x    |

    ---------------------------------------
    16x16_15   |    x    |    x    | ...... |    x    |
    ---------------------------------------
    */

    for (x_search_index = 0; x_search_index < 8; x_search_index++) {
        // 32x32_0
        sad32x32_0 = p_sad16x16[0 * 8 + x_search_index] + p_sad16x16[1 * 8 + x_search_index] +
                     p_sad16x16[2 * 8 + x_search_index] + p_sad16x16[3 * 8 + x_search_index];

        if (sad32x32_0 < p_best_sad_32x32[0]) {
            p_best_sad_32x32[0] = sad32x32_0;
            x_mv                = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x32[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x32_1
        sad32x32_1 = p_sad16x16[4 * 8 + x_search_index] + p_sad16x16[5 * 8 + x_search_index] +
                     p_sad16x16[6 * 8 + x_search_index] + p_sad16x16[7 * 8 + x_search_index];

        if (sad32x32_1 < p_best_sad_32x32[1]) {
            p_best_sad_32x32[1] = sad32x32_1;
            x_mv                = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x32[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x32_2
        sad32x32_2 = p_sad16x16[8 * 8 + x_search_index] + p_sad16x16[9 * 8 + x_search_index] +
                     p_sad16x16[10 * 8 + x_search_index] + p_sad16x16[11 * 8 + x_search_index];

        if (sad32x32_2 < p_best_sad_32x32[2]) {
            p_best_sad_32x32[2] = sad32x32_2;
            x_mv                = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x32[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x32_3
        sad32x32_3 = p_sad16x16[12 * 8 + x_search_index] + p_sad16x16[13 * 8 + x_search_index] +
                     p_sad16x16[14 * 8 + x_search_index] + p_sad16x16[15 * 8 + x_search_index];

        if (sad32x32_3 < p_best_sad_32x32[3]) {
            p_best_sad_32x32[3] = sad32x32_3;
            x_mv                = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x32[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 64x64
        sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
        if (sad64x64 < p_best_sad_64x64[0]) {
            p_best_sad_64x64[0] = sad64x64;
            x_mv                = _MVXT(mv) + (int16_t)x_search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv64x64[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }
}

/*******************************************
Calculate SAD for 16x16 and its 8x8 sublcoks
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void sad_calculation_8x8_16x16_c(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                 uint32_t ref_stride, uint32_t *p_best_sad_8x8,
                                 uint32_t *p_best_sad_16x16, uint32_t *p_best_mv8x8,
                                 uint32_t *p_best_mv16x16, uint32_t mv, uint32_t *p_sad16x16,
                                 EbBool sub_sad) {
    uint64_t sad8x8[4];
    uint64_t sad16x16;

    if (sub_sad) {
        sad8x8[0] =
            (compute8x4_sad_kernel_c(
                src + 0 * src_stride + 0, 2 * src_stride, ref + 0 * ref_stride + 0, 2 * ref_stride))
            << 1;
        sad8x8[1] =
            (compute8x4_sad_kernel_c(
                src + 0 * src_stride + 8, 2 * src_stride, ref + 0 * ref_stride + 8, 2 * ref_stride))
            << 1;
        sad8x8[2] =
            (compute8x4_sad_kernel_c(
                src + 8 * src_stride + 0, 2 * src_stride, ref + 8 * ref_stride + 0, 2 * ref_stride))
            << 1;
        sad8x8[3] =
            (compute8x4_sad_kernel_c(
                src + 8 * src_stride + 8, 2 * src_stride, ref + 8 * ref_stride + 8, 2 * ref_stride))
            << 1;
    } else {
        sad8x8[0] = compute8x8_sad_kernel_c(
            src + 0 * src_stride + 0, src_stride, ref + 0 * ref_stride + 0, ref_stride);
        sad8x8[1] = compute8x8_sad_kernel_c(
            src + 0 * src_stride + 8, src_stride, ref + 0 * ref_stride + 8, ref_stride);
        sad8x8[2] = compute8x8_sad_kernel_c(
            src + 8 * src_stride + 0, src_stride, ref + 8 * ref_stride + 0, ref_stride);
        sad8x8[3] = compute8x8_sad_kernel_c(
            src + 8 * src_stride + 8, src_stride, ref + 8 * ref_stride + 8, ref_stride);
    }

    if (sad8x8[0] < p_best_sad_8x8[0]) {
        p_best_sad_8x8[0] = (uint32_t)sad8x8[0];
        p_best_mv8x8[0]   = mv;
    }

    if (sad8x8[1] < p_best_sad_8x8[1]) {
        p_best_sad_8x8[1] = (uint32_t)sad8x8[1];
        p_best_mv8x8[1]   = mv;
    }

    if (sad8x8[2] < p_best_sad_8x8[2]) {
        p_best_sad_8x8[2] = (uint32_t)sad8x8[2];
        p_best_mv8x8[2]   = mv;
    }

    if (sad8x8[3] < p_best_sad_8x8[3]) {
        p_best_sad_8x8[3] = (uint32_t)sad8x8[3];
        p_best_mv8x8[3]   = mv;
    }

    sad16x16 = sad8x8[0] + sad8x8[1] + sad8x8[2] + sad8x8[3];
    if (sad16x16 < p_best_sad_16x16[0]) {
        p_best_sad_16x16[0] = (uint32_t)sad16x16;
        p_best_mv16x16[0]   = mv;
    }

    *p_sad16x16 = (uint32_t)sad16x16;
}

/*******************************************
Calculate SAD for 32x32,64x64 from 16x16
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void sad_calculation_32x32_64x64_c(uint32_t *p_sad16x16, uint32_t *p_best_sad_32x32,
                                   uint32_t *p_best_sad_64x64, uint32_t *p_best_mv32x32,
                                   uint32_t *p_best_mv64x64, uint32_t mv) {
    uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;

    sad32x32_0 = p_sad16x16[0] + p_sad16x16[1] + p_sad16x16[2] + p_sad16x16[3];
    if (sad32x32_0 < p_best_sad_32x32[0]) {
        p_best_sad_32x32[0] = sad32x32_0;
        p_best_mv32x32[0]   = mv;
    }

    sad32x32_1 = p_sad16x16[4] + p_sad16x16[5] + p_sad16x16[6] + p_sad16x16[7];
    if (sad32x32_1 < p_best_sad_32x32[1]) {
        p_best_sad_32x32[1] = sad32x32_1;
        p_best_mv32x32[1]   = mv;
    }

    sad32x32_2 = p_sad16x16[8] + p_sad16x16[9] + p_sad16x16[10] + p_sad16x16[11];
    if (sad32x32_2 < p_best_sad_32x32[2]) {
        p_best_sad_32x32[2] = sad32x32_2;
        p_best_mv32x32[2]   = mv;
    }

    sad32x32_3 = p_sad16x16[12] + p_sad16x16[13] + p_sad16x16[14] + p_sad16x16[15];
    if (sad32x32_3 < p_best_sad_32x32[3]) {
        p_best_sad_32x32[3] = sad32x32_3;
        p_best_mv32x32[3]   = mv;
    }

    sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
    if (sad64x64 < p_best_sad_64x64[0]) {
        p_best_sad_64x64[0] = sad64x64;
        p_best_mv64x64[0]   = mv;
    }
}

/****************************************************
Calculate SAD for Rect H, V and H4, V4 partitions

and update its Motion info if the result SAD is better
****************************************************/
void ext_sad_calculation(uint32_t *p_sad8x8, uint32_t *p_sad16x16, uint32_t *p_sad32x32,
                         uint32_t *p_best_sad_64x32, uint32_t *p_best_mv64x32,
                         uint32_t *p_best_sad_32x16, uint32_t *p_best_mv32x16,
                         uint32_t *p_best_sad_16x8, uint32_t *p_best_mv16x8,
                         uint32_t *p_best_sad_32x64, uint32_t *p_best_mv32x64,
                         uint32_t *p_best_sad_16x32, uint32_t *p_best_mv16x32,
                         uint32_t *p_best_sad_8x16, uint32_t *p_best_mv8x16,
                         uint32_t *p_best_sad_32x8, uint32_t *p_best_mv32x8,
                         uint32_t *p_best_sad_8x32, uint32_t *p_best_mv8x32,
                         uint32_t *p_best_sad_64x16, uint32_t *p_best_mv64x16,
                         uint32_t *p_best_sad_16x64, uint32_t *p_best_mv16x64, uint32_t mv) {
    uint32_t sad;

    uint32_t sad_16x8[32];
    uint32_t sad_8x16[32];
    uint32_t sad_32x16[8];
    uint32_t sad_16x32[8];

    // 64x32
    sad = p_sad32x32[0] + p_sad32x32[1];
    if (sad < p_best_sad_64x32[0]) {
        p_best_sad_64x32[0] = sad;
        p_best_mv64x32[0]   = mv;
    }

    sad = p_sad32x32[2] + p_sad32x32[3];
    if (sad < p_best_sad_64x32[1]) {
        p_best_sad_64x32[1] = sad;
        p_best_mv64x32[1]   = mv;
    }

    // 32x16
    sad_32x16[0] = p_sad16x16[0] + p_sad16x16[1];
    if (sad_32x16[0] < p_best_sad_32x16[0]) {
        p_best_sad_32x16[0] = sad_32x16[0];
        p_best_mv32x16[0]   = mv;
    }

    sad_32x16[1] = p_sad16x16[2] + p_sad16x16[3];
    if (sad_32x16[1] < p_best_sad_32x16[1]) {
        p_best_sad_32x16[1] = sad_32x16[1];
        p_best_mv32x16[1]   = mv;
    }

    sad_32x16[2] = p_sad16x16[4] + p_sad16x16[5];
    if (sad_32x16[2] < p_best_sad_32x16[2]) {
        p_best_sad_32x16[2] = sad_32x16[2];
        p_best_mv32x16[2]   = mv;
    }

    sad_32x16[3] = p_sad16x16[6] + p_sad16x16[7];
    if (sad_32x16[3] < p_best_sad_32x16[3]) {
        p_best_sad_32x16[3] = sad_32x16[3];
        p_best_mv32x16[3]   = mv;
    }

    sad_32x16[4] = p_sad16x16[8] + p_sad16x16[9];
    if (sad_32x16[4] < p_best_sad_32x16[4]) {
        p_best_sad_32x16[4] = sad_32x16[4];
        p_best_mv32x16[4]   = mv;
    }

    sad_32x16[5] = p_sad16x16[10] + p_sad16x16[11];
    if (sad_32x16[5] < p_best_sad_32x16[5]) {
        p_best_sad_32x16[5] = sad_32x16[5];
        p_best_mv32x16[5]   = mv;
    }

    sad_32x16[6] = p_sad16x16[12] + p_sad16x16[13];
    if (sad_32x16[6] < p_best_sad_32x16[6]) {
        p_best_sad_32x16[6] = sad_32x16[6];
        p_best_mv32x16[6]   = mv;
    }

    sad_32x16[7] = p_sad16x16[14] + p_sad16x16[15];
    if (sad_32x16[7] < p_best_sad_32x16[7]) {
        p_best_sad_32x16[7] = sad_32x16[7];
        p_best_mv32x16[7]   = mv;
    }

    // 64x16
    sad = sad_32x16[0] + sad_32x16[2];
    if (sad < p_best_sad_64x16[0]) {
        p_best_sad_64x16[0] = sad;
        p_best_mv64x16[0]   = mv;
    }
    sad = sad_32x16[1] + sad_32x16[3];
    if (sad < p_best_sad_64x16[1]) {
        p_best_sad_64x16[1] = sad;
        p_best_mv64x16[1]   = mv;
    }

    sad = sad_32x16[4] + sad_32x16[6];
    if (sad < p_best_sad_64x16[2]) {
        p_best_sad_64x16[2] = sad;
        p_best_mv64x16[2]   = mv;
    }
    sad = sad_32x16[5] + sad_32x16[7];
    if (sad < p_best_sad_64x16[3]) {
        p_best_sad_64x16[3] = sad;
        p_best_mv64x16[3]   = mv;
    }

    // 16x8
    sad_16x8[0] = p_sad8x8[0] + p_sad8x8[1];
    if (sad_16x8[0] < p_best_sad_16x8[0]) {
        p_best_sad_16x8[0] = sad_16x8[0];
        p_best_mv16x8[0]   = mv;
    }

    sad_16x8[1] = p_sad8x8[2] + p_sad8x8[3];
    if (sad_16x8[1] < p_best_sad_16x8[1]) {
        p_best_sad_16x8[1] = sad_16x8[1];
        p_best_mv16x8[1]   = mv;
    }

    sad_16x8[2] = p_sad8x8[4] + p_sad8x8[5];
    if (sad_16x8[2] < p_best_sad_16x8[2]) {
        p_best_sad_16x8[2] = sad_16x8[2];
        p_best_mv16x8[2]   = mv;
    }

    sad_16x8[3] = p_sad8x8[6] + p_sad8x8[7];
    if (sad_16x8[3] < p_best_sad_16x8[3]) {
        p_best_sad_16x8[3] = sad_16x8[3];
        p_best_mv16x8[3]   = mv;
    }

    sad_16x8[4] = p_sad8x8[8] + p_sad8x8[9];
    if (sad_16x8[4] < p_best_sad_16x8[4]) {
        p_best_sad_16x8[4] = sad_16x8[4];
        p_best_mv16x8[4]   = mv;
    }

    sad_16x8[5] = p_sad8x8[10] + p_sad8x8[11];
    if (sad_16x8[5] < p_best_sad_16x8[5]) {
        p_best_sad_16x8[5] = sad_16x8[5];
        p_best_mv16x8[5]   = mv;
    }

    sad_16x8[6] = p_sad8x8[12] + p_sad8x8[13];
    if (sad_16x8[6] < p_best_sad_16x8[6]) {
        p_best_sad_16x8[6] = sad_16x8[6];
        p_best_mv16x8[6]   = mv;
    }

    sad_16x8[7] = p_sad8x8[14] + p_sad8x8[15];
    if (sad_16x8[7] < p_best_sad_16x8[7]) {
        p_best_sad_16x8[7] = sad_16x8[7];
        p_best_mv16x8[7]   = mv;
    }

    sad_16x8[8] = p_sad8x8[16] + p_sad8x8[17];
    if (sad_16x8[8] < p_best_sad_16x8[8]) {
        p_best_sad_16x8[8] = sad_16x8[8];
        p_best_mv16x8[8]   = mv;
    }

    sad_16x8[9] = p_sad8x8[18] + p_sad8x8[19];
    if (sad_16x8[9] < p_best_sad_16x8[9]) {
        p_best_sad_16x8[9] = sad_16x8[9];
        p_best_mv16x8[9]   = mv;
    }

    sad_16x8[10] = p_sad8x8[20] + p_sad8x8[21];
    if (sad_16x8[10] < p_best_sad_16x8[10]) {
        p_best_sad_16x8[10] = sad_16x8[10];
        p_best_mv16x8[10]   = mv;
    }

    sad_16x8[11] = p_sad8x8[22] + p_sad8x8[23];
    if (sad_16x8[11] < p_best_sad_16x8[11]) {
        p_best_sad_16x8[11] = sad_16x8[11];
        p_best_mv16x8[11]   = mv;
    }

    sad_16x8[12] = p_sad8x8[24] + p_sad8x8[25];
    if (sad_16x8[12] < p_best_sad_16x8[12]) {
        p_best_sad_16x8[12] = sad_16x8[12];
        p_best_mv16x8[12]   = mv;
    }

    sad_16x8[13] = p_sad8x8[26] + p_sad8x8[27];
    if (sad_16x8[13] < p_best_sad_16x8[13]) {
        p_best_sad_16x8[13] = sad_16x8[13];
        p_best_mv16x8[13]   = mv;
    }

    sad_16x8[14] = p_sad8x8[28] + p_sad8x8[29];
    if (sad_16x8[14] < p_best_sad_16x8[14]) {
        p_best_sad_16x8[14] = sad_16x8[14];
        p_best_mv16x8[14]   = mv;
    }

    sad_16x8[15] = p_sad8x8[30] + p_sad8x8[31];
    if (sad_16x8[15] < p_best_sad_16x8[15]) {
        p_best_sad_16x8[15] = sad_16x8[15];
        p_best_mv16x8[15]   = mv;
    }

    sad_16x8[16] = p_sad8x8[32] + p_sad8x8[33];
    if (sad_16x8[16] < p_best_sad_16x8[16]) {
        p_best_sad_16x8[16] = sad_16x8[16];
        p_best_mv16x8[16]   = mv;
    }

    sad_16x8[17] = p_sad8x8[34] + p_sad8x8[35];
    if (sad_16x8[17] < p_best_sad_16x8[17]) {
        p_best_sad_16x8[17] = sad_16x8[17];
        p_best_mv16x8[17]   = mv;
    }

    sad_16x8[18] = p_sad8x8[36] + p_sad8x8[37];
    if (sad_16x8[18] < p_best_sad_16x8[18]) {
        p_best_sad_16x8[18] = sad_16x8[18];
        p_best_mv16x8[18]   = mv;
    }

    sad_16x8[19] = p_sad8x8[38] + p_sad8x8[39];
    if (sad_16x8[19] < p_best_sad_16x8[19]) {
        p_best_sad_16x8[19] = sad_16x8[19];
        p_best_mv16x8[19]   = mv;
    }

    sad_16x8[20] = p_sad8x8[40] + p_sad8x8[41];
    if (sad_16x8[20] < p_best_sad_16x8[20]) {
        p_best_sad_16x8[20] = sad_16x8[20];
        p_best_mv16x8[20]   = mv;
    }

    sad_16x8[21] = p_sad8x8[42] + p_sad8x8[43];
    if (sad_16x8[21] < p_best_sad_16x8[21]) {
        p_best_sad_16x8[21] = sad_16x8[21];
        p_best_mv16x8[21]   = mv;
    }

    sad_16x8[22] = p_sad8x8[44] + p_sad8x8[45];
    if (sad_16x8[22] < p_best_sad_16x8[22]) {
        p_best_sad_16x8[22] = sad_16x8[22];
        p_best_mv16x8[22]   = mv;
    }

    sad_16x8[23] = p_sad8x8[46] + p_sad8x8[47];
    if (sad_16x8[23] < p_best_sad_16x8[23]) {
        p_best_sad_16x8[23] = sad_16x8[23];
        p_best_mv16x8[23]   = mv;
    }

    sad_16x8[24] = p_sad8x8[48] + p_sad8x8[49];
    if (sad_16x8[24] < p_best_sad_16x8[24]) {
        p_best_sad_16x8[24] = sad_16x8[24];
        p_best_mv16x8[24]   = mv;
    }

    sad_16x8[25] = p_sad8x8[50] + p_sad8x8[51];
    if (sad_16x8[25] < p_best_sad_16x8[25]) {
        p_best_sad_16x8[25] = sad_16x8[25];
        p_best_mv16x8[25]   = mv;
    }

    sad_16x8[26] = p_sad8x8[52] + p_sad8x8[53];
    if (sad_16x8[26] < p_best_sad_16x8[26]) {
        p_best_sad_16x8[26] = sad_16x8[26];
        p_best_mv16x8[26]   = mv;
    }

    sad_16x8[27] = p_sad8x8[54] + p_sad8x8[55];
    if (sad_16x8[27] < p_best_sad_16x8[27]) {
        p_best_sad_16x8[27] = sad_16x8[27];
        p_best_mv16x8[27]   = mv;
    }

    sad_16x8[28] = p_sad8x8[56] + p_sad8x8[57];
    if (sad_16x8[28] < p_best_sad_16x8[28]) {
        p_best_sad_16x8[28] = sad_16x8[28];
        p_best_mv16x8[28]   = mv;
    }

    sad_16x8[29] = p_sad8x8[58] + p_sad8x8[59];
    if (sad_16x8[29] < p_best_sad_16x8[29]) {
        p_best_sad_16x8[29] = sad_16x8[29];
        p_best_mv16x8[29]   = mv;
    }

    sad_16x8[30] = p_sad8x8[60] + p_sad8x8[61];
    if (sad_16x8[30] < p_best_sad_16x8[30]) {
        p_best_sad_16x8[30] = sad_16x8[30];
        p_best_mv16x8[30]   = mv;
    }

    sad_16x8[31] = p_sad8x8[62] + p_sad8x8[63];
    if (sad_16x8[31] < p_best_sad_16x8[31]) {
        p_best_sad_16x8[31] = sad_16x8[31];
        p_best_mv16x8[31]   = mv;
    }

    // 32x64
    sad = p_sad32x32[0] + p_sad32x32[2];
    if (sad < p_best_sad_32x64[0]) {
        p_best_sad_32x64[0] = sad;
        p_best_mv32x64[0]   = mv;
    }

    sad = p_sad32x32[1] + p_sad32x32[3];
    if (sad < p_best_sad_32x64[1]) {
        p_best_sad_32x64[1] = sad;
        p_best_mv32x64[1]   = mv;
    }

    // 16x32
    sad_16x32[0] = p_sad16x16[0] + p_sad16x16[2];
    if (sad_16x32[0] < p_best_sad_16x32[0]) {
        p_best_sad_16x32[0] = sad_16x32[0];
        p_best_mv16x32[0]   = mv;
    }

    sad_16x32[1] = p_sad16x16[1] + p_sad16x16[3];
    if (sad_16x32[1] < p_best_sad_16x32[1]) {
        p_best_sad_16x32[1] = sad_16x32[1];
        p_best_mv16x32[1]   = mv;
    }

    sad_16x32[2] = p_sad16x16[4] + p_sad16x16[6];
    if (sad_16x32[2] < p_best_sad_16x32[2]) {
        p_best_sad_16x32[2] = sad_16x32[2];
        p_best_mv16x32[2]   = mv;
    }

    sad_16x32[3] = p_sad16x16[5] + p_sad16x16[7];
    if (sad_16x32[3] < p_best_sad_16x32[3]) {
        p_best_sad_16x32[3] = sad_16x32[3];
        p_best_mv16x32[3]   = mv;
    }

    sad_16x32[4] = p_sad16x16[8] + p_sad16x16[10];
    if (sad_16x32[4] < p_best_sad_16x32[4]) {
        p_best_sad_16x32[4] = sad_16x32[4];
        p_best_mv16x32[4]   = mv;
    }

    sad_16x32[5] = p_sad16x16[9] + p_sad16x16[11];
    if (sad_16x32[5] < p_best_sad_16x32[5]) {
        p_best_sad_16x32[5] = sad_16x32[5];
        p_best_mv16x32[5]   = mv;
    }

    sad_16x32[6] = p_sad16x16[12] + p_sad16x16[14];
    if (sad_16x32[6] < p_best_sad_16x32[6]) {
        p_best_sad_16x32[6] = sad_16x32[6];
        p_best_mv16x32[6]   = mv;
    }

    sad_16x32[7] = p_sad16x16[13] + p_sad16x16[15];
    if (sad_16x32[7] < p_best_sad_16x32[7]) {
        p_best_sad_16x32[7] = sad_16x32[7];
        p_best_mv16x32[7]   = mv;
    }

    sad = sad_16x32[0] + sad_16x32[4];
    if (sad < p_best_sad_16x64[0]) {
        p_best_sad_16x64[0] = sad;
        p_best_mv16x64[0]   = mv;
    }
    sad = sad_16x32[1] + sad_16x32[5];
    if (sad < p_best_sad_16x64[1]) {
        p_best_sad_16x64[1] = sad;
        p_best_mv16x64[1]   = mv;
    }

    sad = sad_16x32[2] + sad_16x32[6];
    if (sad < p_best_sad_16x64[2]) {
        p_best_sad_16x64[2] = sad;
        p_best_mv16x64[2]   = mv;
    }

    sad = sad_16x32[3] + sad_16x32[7];
    if (sad < p_best_sad_16x64[3]) {
        p_best_sad_16x64[3] = sad;
        p_best_mv16x64[3]   = mv;
    }

    // 8x16
    sad_8x16[0] = p_sad8x8[0] + p_sad8x8[2];
    if (sad_8x16[0] < p_best_sad_8x16[0]) {
        p_best_sad_8x16[0] = sad_8x16[0];
        p_best_mv8x16[0]   = mv;
    }

    sad_8x16[1] = p_sad8x8[1] + p_sad8x8[3];
    if (sad_8x16[1] < p_best_sad_8x16[1]) {
        p_best_sad_8x16[1] = sad_8x16[1];
        p_best_mv8x16[1]   = mv;
    }

    sad_8x16[2] = p_sad8x8[4] + p_sad8x8[6];
    if (sad_8x16[2] < p_best_sad_8x16[2]) {
        p_best_sad_8x16[2] = sad_8x16[2];
        p_best_mv8x16[2]   = mv;
    }

    sad_8x16[3] = p_sad8x8[5] + p_sad8x8[7];
    if (sad_8x16[3] < p_best_sad_8x16[3]) {
        p_best_sad_8x16[3] = sad_8x16[3];
        p_best_mv8x16[3]   = mv;
    }

    sad_8x16[4] = p_sad8x8[8] + p_sad8x8[10];
    if (sad_8x16[4] < p_best_sad_8x16[4]) {
        p_best_sad_8x16[4] = sad_8x16[4];
        p_best_mv8x16[4]   = mv;
    }

    sad_8x16[5] = p_sad8x8[9] + p_sad8x8[11];
    if (sad_8x16[5] < p_best_sad_8x16[5]) {
        p_best_sad_8x16[5] = sad_8x16[5];
        p_best_mv8x16[5]   = mv;
    }

    sad_8x16[6] = p_sad8x8[12] + p_sad8x8[14];
    if (sad_8x16[6] < p_best_sad_8x16[6]) {
        p_best_sad_8x16[6] = sad_8x16[6];
        p_best_mv8x16[6]   = mv;
    }

    sad_8x16[7] = p_sad8x8[13] + p_sad8x8[15];
    if (sad_8x16[7] < p_best_sad_8x16[7]) {
        p_best_sad_8x16[7] = sad_8x16[7];
        p_best_mv8x16[7]   = mv;
    }

    sad_8x16[8] = p_sad8x8[16] + p_sad8x8[18];
    if (sad_8x16[8] < p_best_sad_8x16[8]) {
        p_best_sad_8x16[8] = sad_8x16[8];
        p_best_mv8x16[8]   = mv;
    }

    sad_8x16[9] = p_sad8x8[17] + p_sad8x8[19];
    if (sad_8x16[9] < p_best_sad_8x16[9]) {
        p_best_sad_8x16[9] = sad_8x16[9];
        p_best_mv8x16[9]   = mv;
    }

    sad_8x16[10] = p_sad8x8[20] + p_sad8x8[22];
    if (sad_8x16[10] < p_best_sad_8x16[10]) {
        p_best_sad_8x16[10] = sad_8x16[10];
        p_best_mv8x16[10]   = mv;
    }

    sad_8x16[11] = p_sad8x8[21] + p_sad8x8[23];
    if (sad_8x16[11] < p_best_sad_8x16[11]) {
        p_best_sad_8x16[11] = sad_8x16[11];
        p_best_mv8x16[11]   = mv;
    }

    sad_8x16[12] = p_sad8x8[24] + p_sad8x8[26];
    if (sad_8x16[12] < p_best_sad_8x16[12]) {
        p_best_sad_8x16[12] = sad_8x16[12];
        p_best_mv8x16[12]   = mv;
    }

    sad_8x16[13] = p_sad8x8[25] + p_sad8x8[27];
    if (sad_8x16[13] < p_best_sad_8x16[13]) {
        p_best_sad_8x16[13] = sad_8x16[13];
        p_best_mv8x16[13]   = mv;
    }

    sad_8x16[14] = p_sad8x8[28] + p_sad8x8[30];
    if (sad_8x16[14] < p_best_sad_8x16[14]) {
        p_best_sad_8x16[14] = sad_8x16[14];
        p_best_mv8x16[14]   = mv;
    }

    sad_8x16[15] = p_sad8x8[29] + p_sad8x8[31];
    if (sad_8x16[15] < p_best_sad_8x16[15]) {
        p_best_sad_8x16[15] = sad_8x16[15];
        p_best_mv8x16[15]   = mv;
    }

    sad_8x16[16] = p_sad8x8[32] + p_sad8x8[34];
    if (sad_8x16[16] < p_best_sad_8x16[16]) {
        p_best_sad_8x16[16] = sad_8x16[16];
        p_best_mv8x16[16]   = mv;
    }

    sad_8x16[17] = p_sad8x8[33] + p_sad8x8[35];
    if (sad_8x16[17] < p_best_sad_8x16[17]) {
        p_best_sad_8x16[17] = sad_8x16[17];
        p_best_mv8x16[17]   = mv;
    }

    sad_8x16[18] = p_sad8x8[36] + p_sad8x8[38];
    if (sad_8x16[18] < p_best_sad_8x16[18]) {
        p_best_sad_8x16[18] = sad_8x16[18];
        p_best_mv8x16[18]   = mv;
    }

    sad_8x16[19] = p_sad8x8[37] + p_sad8x8[39];
    if (sad_8x16[19] < p_best_sad_8x16[19]) {
        p_best_sad_8x16[19] = sad_8x16[19];
        p_best_mv8x16[19]   = mv;
    }

    sad_8x16[20] = p_sad8x8[40] + p_sad8x8[42];
    if (sad_8x16[20] < p_best_sad_8x16[20]) {
        p_best_sad_8x16[20] = sad_8x16[20];
        p_best_mv8x16[20]   = mv;
    }

    sad_8x16[21] = p_sad8x8[41] + p_sad8x8[43];
    if (sad_8x16[21] < p_best_sad_8x16[21]) {
        p_best_sad_8x16[21] = sad_8x16[21];
        p_best_mv8x16[21]   = mv;
    }

    sad_8x16[22] = p_sad8x8[44] + p_sad8x8[46];
    if (sad_8x16[22] < p_best_sad_8x16[22]) {
        p_best_sad_8x16[22] = sad_8x16[22];
        p_best_mv8x16[22]   = mv;
    }

    sad_8x16[23] = p_sad8x8[45] + p_sad8x8[47];
    if (sad_8x16[23] < p_best_sad_8x16[23]) {
        p_best_sad_8x16[23] = sad_8x16[23];
        p_best_mv8x16[23]   = mv;
    }

    sad_8x16[24] = p_sad8x8[48] + p_sad8x8[50];
    if (sad_8x16[24] < p_best_sad_8x16[24]) {
        p_best_sad_8x16[24] = sad_8x16[24];
        p_best_mv8x16[24]   = mv;
    }

    sad_8x16[25] = p_sad8x8[49] + p_sad8x8[51];
    if (sad_8x16[25] < p_best_sad_8x16[25]) {
        p_best_sad_8x16[25] = sad_8x16[25];
        p_best_mv8x16[25]   = mv;
    }

    sad_8x16[26] = p_sad8x8[52] + p_sad8x8[54];
    if (sad_8x16[26] < p_best_sad_8x16[26]) {
        p_best_sad_8x16[26] = sad_8x16[26];
        p_best_mv8x16[26]   = mv;
    }

    sad_8x16[27] = p_sad8x8[53] + p_sad8x8[55];
    if (sad_8x16[27] < p_best_sad_8x16[27]) {
        p_best_sad_8x16[27] = sad_8x16[27];
        p_best_mv8x16[27]   = mv;
    }

    sad_8x16[28] = p_sad8x8[56] + p_sad8x8[58];
    if (sad_8x16[28] < p_best_sad_8x16[28]) {
        p_best_sad_8x16[28] = sad_8x16[28];
        p_best_mv8x16[28]   = mv;
    }

    sad_8x16[29] = p_sad8x8[57] + p_sad8x8[59];
    if (sad_8x16[29] < p_best_sad_8x16[29]) {
        p_best_sad_8x16[29] = sad_8x16[29];
        p_best_mv8x16[29]   = mv;
    }

    sad_8x16[30] = p_sad8x8[60] + p_sad8x8[62];
    if (sad_8x16[30] < p_best_sad_8x16[30]) {
        p_best_sad_8x16[30] = sad_8x16[30];
        p_best_mv8x16[30]   = mv;
    }

    sad_8x16[31] = p_sad8x8[61] + p_sad8x8[63];
    if (sad_8x16[31] < p_best_sad_8x16[31]) {
        p_best_sad_8x16[31] = sad_8x16[31];
        p_best_mv8x16[31]   = mv;
    }

    // 32x8
    sad = sad_16x8[0] + sad_16x8[2];
    if (sad < p_best_sad_32x8[0]) {
        p_best_sad_32x8[0] = sad;
        p_best_mv32x8[0]   = mv;
    }

    sad = sad_16x8[1] + sad_16x8[3];
    if (sad < p_best_sad_32x8[1]) {
        p_best_sad_32x8[1] = sad;
        p_best_mv32x8[1]   = mv;
    }

    sad = sad_16x8[4] + sad_16x8[6];
    if (sad < p_best_sad_32x8[2]) {
        p_best_sad_32x8[2] = sad;
        p_best_mv32x8[2]   = mv;
    }

    sad = sad_16x8[5] + sad_16x8[7];
    if (sad < p_best_sad_32x8[3]) {
        p_best_sad_32x8[3] = sad;
        p_best_mv32x8[3]   = mv;
    }

    sad = sad_16x8[8] + sad_16x8[10];
    if (sad < p_best_sad_32x8[4]) {
        p_best_sad_32x8[4] = sad;
        p_best_mv32x8[4]   = mv;
    }

    sad = sad_16x8[9] + sad_16x8[11];
    if (sad < p_best_sad_32x8[5]) {
        p_best_sad_32x8[5] = sad;
        p_best_mv32x8[5]   = mv;
    }

    sad = sad_16x8[12] + sad_16x8[14];
    if (sad < p_best_sad_32x8[6]) {
        p_best_sad_32x8[6] = sad;
        p_best_mv32x8[6]   = mv;
    }

    sad = sad_16x8[13] + sad_16x8[15];
    if (sad < p_best_sad_32x8[7]) {
        p_best_sad_32x8[7] = sad;
        p_best_mv32x8[7]   = mv;
    }

    sad = sad_16x8[16] + sad_16x8[18];
    if (sad < p_best_sad_32x8[8]) {
        p_best_sad_32x8[8] = sad;
        p_best_mv32x8[8]   = mv;
    }

    sad = sad_16x8[17] + sad_16x8[19];
    if (sad < p_best_sad_32x8[9]) {
        p_best_sad_32x8[9] = sad;
        p_best_mv32x8[9]   = mv;
    }

    sad = sad_16x8[20] + sad_16x8[22];
    if (sad < p_best_sad_32x8[10]) {
        p_best_sad_32x8[10] = sad;
        p_best_mv32x8[10]   = mv;
    }

    sad = sad_16x8[21] + sad_16x8[23];
    if (sad < p_best_sad_32x8[11]) {
        p_best_sad_32x8[11] = sad;
        p_best_mv32x8[11]   = mv;
    }

    sad = sad_16x8[24] + sad_16x8[26];
    if (sad < p_best_sad_32x8[12]) {
        p_best_sad_32x8[12] = sad;
        p_best_mv32x8[12]   = mv;
    }

    sad = sad_16x8[25] + sad_16x8[27];
    if (sad < p_best_sad_32x8[13]) {
        p_best_sad_32x8[13] = sad;
        p_best_mv32x8[13]   = mv;
    }

    sad = sad_16x8[28] + sad_16x8[30];
    if (sad < p_best_sad_32x8[14]) {
        p_best_sad_32x8[14] = sad;
        p_best_mv32x8[14]   = mv;
    }

    sad = sad_16x8[29] + sad_16x8[31];
    if (sad < p_best_sad_32x8[15]) {
        p_best_sad_32x8[15] = sad;
        p_best_mv32x8[15]   = mv;
    }

    // 8x32
    sad = sad_8x16[0] + sad_8x16[4];
    if (sad < p_best_sad_8x32[0]) {
        p_best_sad_8x32[0] = sad;
        p_best_mv8x32[0]   = mv;
    }

    sad = sad_8x16[1] + sad_8x16[5];
    if (sad < p_best_sad_8x32[1]) {
        p_best_sad_8x32[1] = sad;
        p_best_mv8x32[1]   = mv;
    }

    sad = sad_8x16[2] + sad_8x16[6];
    if (sad < p_best_sad_8x32[2]) {
        p_best_sad_8x32[2] = sad;
        p_best_mv8x32[2]   = mv;
    }

    sad = sad_8x16[3] + sad_8x16[7];
    if (sad < p_best_sad_8x32[3]) {
        p_best_sad_8x32[3] = sad;
        p_best_mv8x32[3]   = mv;
    }

    sad = sad_8x16[8] + sad_8x16[12];
    if (sad < p_best_sad_8x32[4]) {
        p_best_sad_8x32[4] = sad;
        p_best_mv8x32[4]   = mv;
    }

    sad = sad_8x16[9] + sad_8x16[13];
    if (sad < p_best_sad_8x32[5]) {
        p_best_sad_8x32[5] = sad;
        p_best_mv8x32[5]   = mv;
    }

    sad = sad_8x16[10] + sad_8x16[14];
    if (sad < p_best_sad_8x32[6]) {
        p_best_sad_8x32[6] = sad;
        p_best_mv8x32[6]   = mv;
    }

    sad = sad_8x16[11] + sad_8x16[15];
    if (sad < p_best_sad_8x32[7]) {
        p_best_sad_8x32[7] = sad;
        p_best_mv8x32[7]   = mv;
    }

    sad = sad_8x16[16] + sad_8x16[20];
    if (sad < p_best_sad_8x32[8]) {
        p_best_sad_8x32[8] = sad;
        p_best_mv8x32[8]   = mv;
    }

    sad = sad_8x16[17] + sad_8x16[21];
    if (sad < p_best_sad_8x32[9]) {
        p_best_sad_8x32[9] = sad;
        p_best_mv8x32[9]   = mv;
    }

    sad = sad_8x16[18] + sad_8x16[22];
    if (sad < p_best_sad_8x32[10]) {
        p_best_sad_8x32[10] = sad;
        p_best_mv8x32[10]   = mv;
    }

    sad = sad_8x16[19] + sad_8x16[23];
    if (sad < p_best_sad_8x32[11]) {
        p_best_sad_8x32[11] = sad;
        p_best_mv8x32[11]   = mv;
    }

    sad = sad_8x16[24] + sad_8x16[28];
    if (sad < p_best_sad_8x32[12]) {
        p_best_sad_8x32[12] = sad;
        p_best_mv8x32[12]   = mv;
    }

    sad = sad_8x16[25] + sad_8x16[29];
    if (sad < p_best_sad_8x32[13]) {
        p_best_sad_8x32[13] = sad;
        p_best_mv8x32[13]   = mv;
    }

    sad = sad_8x16[26] + sad_8x16[30];
    if (sad < p_best_sad_8x32[14]) {
        p_best_sad_8x32[14] = sad;
        p_best_mv8x32[14]   = mv;
    }

    sad = sad_8x16[27] + sad_8x16[31];
    if (sad < p_best_sad_8x32[15]) {
        p_best_sad_8x32[15] = sad;
        p_best_mv8x32[15]   = mv;
    }
}

/****************************************************
Calculate SAD for Rect H, V and H4, V4 partitions
and update its Motion info if the result SAD is better
****************************************************/
void ext_eigth_sad_calculation_nsq_c(
    uint32_t p_sad8x8[64][8], uint32_t p_sad16x16[16][8], uint32_t p_sad32x32[4][8],
    uint32_t *p_best_sad_64x32, uint32_t *p_best_mv64x32, uint32_t *p_best_sad_32x16,
    uint32_t *p_best_mv32x16, uint32_t *p_best_sad_16x8, uint32_t *p_best_mv16x8,
    uint32_t *p_best_sad_32x64, uint32_t *p_best_mv32x64, uint32_t *p_best_sad_16x32,
    uint32_t *p_best_mv16x32, uint32_t *p_best_sad_8x16, uint32_t *p_best_mv8x16,
    uint32_t *p_best_sad_32x8, uint32_t *p_best_mv32x8, uint32_t *p_best_sad_8x32,
    uint32_t *p_best_mv8x32, uint32_t *p_best_sad_64x16, uint32_t *p_best_mv64x16,
    uint32_t *p_best_sad_16x64, uint32_t *p_best_mv16x64, uint32_t mv) {
    uint8_t  search_index;
    uint32_t sad;
    uint32_t sad_16x8[32];
    uint32_t sad_8x16[32];
    uint32_t sad_32x16[8];
    uint32_t sad_16x32[8];

    int16_t x_mv, y_mv;

    for (search_index = 0; search_index < 8; search_index++) {
        // 64x32
        sad = p_sad32x32[0][search_index] + p_sad32x32[1][search_index];
        if (sad < p_best_sad_64x32[0]) {
            p_best_sad_64x32[0] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv64x32[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = p_sad32x32[2][search_index] + p_sad32x32[3][search_index];
        if (sad < p_best_sad_64x32[1]) {
            p_best_sad_64x32[1] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv64x32[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x16
        sad_32x16[0] = p_sad16x16[0][search_index] + p_sad16x16[1][search_index];
        if (sad_32x16[0] < p_best_sad_32x16[0]) {
            p_best_sad_32x16[0] = sad_32x16[0];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x16[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[1] = p_sad16x16[2][search_index] + p_sad16x16[3][search_index];
        if (sad_32x16[1] < p_best_sad_32x16[1]) {
            p_best_sad_32x16[1] = sad_32x16[1];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x16[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[2] = p_sad16x16[4][search_index] + p_sad16x16[5][search_index];
        if (sad_32x16[2] < p_best_sad_32x16[2]) {
            p_best_sad_32x16[2] = sad_32x16[2];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x16[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[3] = p_sad16x16[6][search_index] + p_sad16x16[7][search_index];
        if (sad_32x16[3] < p_best_sad_32x16[3]) {
            p_best_sad_32x16[3] = sad_32x16[3];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x16[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[4] = p_sad16x16[8][search_index] + p_sad16x16[9][search_index];
        if (sad_32x16[4] < p_best_sad_32x16[4]) {
            p_best_sad_32x16[4] = sad_32x16[4];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x16[4]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[5] = p_sad16x16[10][search_index] + p_sad16x16[11][search_index];
        if (sad_32x16[5] < p_best_sad_32x16[5]) {
            p_best_sad_32x16[5] = sad_32x16[5];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x16[5]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[6] = p_sad16x16[12][search_index] + p_sad16x16[13][search_index];
        if (sad_32x16[6] < p_best_sad_32x16[6]) {
            p_best_sad_32x16[6] = sad_32x16[6];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x16[6]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[7] = p_sad16x16[14][search_index] + p_sad16x16[15][search_index];
        if (sad_32x16[7] < p_best_sad_32x16[7]) {
            p_best_sad_32x16[7] = sad_32x16[7];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x16[7]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 64x16
        sad = sad_32x16[0] + sad_32x16[2];
        if (sad < p_best_sad_64x16[0]) {
            p_best_sad_64x16[0] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv64x16[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        sad = sad_32x16[1] + sad_32x16[3];
        if (sad < p_best_sad_64x16[1]) {
            p_best_sad_64x16[1] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv64x16[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_32x16[4] + sad_32x16[6];
        if (sad < p_best_sad_64x16[2]) {
            p_best_sad_64x16[2] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv64x16[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        sad = sad_32x16[5] + sad_32x16[7];
        if (sad < p_best_sad_64x16[3]) {
            p_best_sad_64x16[3] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv64x16[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        // 16x8
        sad_16x8[0] = p_sad8x8[0][search_index] + p_sad8x8[1][search_index];
        if (sad_16x8[0] < p_best_sad_16x8[0]) {
            p_best_sad_16x8[0] = sad_16x8[0];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[1] = p_sad8x8[2][search_index] + p_sad8x8[3][search_index];
        if (sad_16x8[1] < p_best_sad_16x8[1]) {
            p_best_sad_16x8[1] = sad_16x8[1];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[2] = p_sad8x8[4][search_index] + p_sad8x8[5][search_index];
        if (sad_16x8[2] < p_best_sad_16x8[2]) {
            p_best_sad_16x8[2] = sad_16x8[2];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[3] = p_sad8x8[6][search_index] + p_sad8x8[7][search_index];
        if (sad_16x8[3] < p_best_sad_16x8[3]) {
            p_best_sad_16x8[3] = sad_16x8[3];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[4] = p_sad8x8[8][search_index] + p_sad8x8[9][search_index];
        if (sad_16x8[4] < p_best_sad_16x8[4]) {
            p_best_sad_16x8[4] = sad_16x8[4];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[4]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[5] = p_sad8x8[10][search_index] + p_sad8x8[11][search_index];
        if (sad_16x8[5] < p_best_sad_16x8[5]) {
            p_best_sad_16x8[5] = sad_16x8[5];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[5]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[6] = p_sad8x8[12][search_index] + p_sad8x8[13][search_index];
        if (sad_16x8[6] < p_best_sad_16x8[6]) {
            p_best_sad_16x8[6] = sad_16x8[6];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[6]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[7] = p_sad8x8[14][search_index] + p_sad8x8[15][search_index];
        if (sad_16x8[7] < p_best_sad_16x8[7]) {
            p_best_sad_16x8[7] = sad_16x8[7];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[7]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[8] = p_sad8x8[16][search_index] + p_sad8x8[17][search_index];
        if (sad_16x8[8] < p_best_sad_16x8[8]) {
            p_best_sad_16x8[8] = sad_16x8[8];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[8]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[9] = p_sad8x8[18][search_index] + p_sad8x8[19][search_index];
        if (sad_16x8[9] < p_best_sad_16x8[9]) {
            p_best_sad_16x8[9] = sad_16x8[9];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv16x8[9]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[10] = p_sad8x8[20][search_index] + p_sad8x8[21][search_index];
        if (sad_16x8[10] < p_best_sad_16x8[10]) {
            p_best_sad_16x8[10] = sad_16x8[10];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[10]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[11] = p_sad8x8[22][search_index] + p_sad8x8[23][search_index];
        if (sad_16x8[11] < p_best_sad_16x8[11]) {
            p_best_sad_16x8[11] = sad_16x8[11];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[11]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[12] = p_sad8x8[24][search_index] + p_sad8x8[25][search_index];
        if (sad_16x8[12] < p_best_sad_16x8[12]) {
            p_best_sad_16x8[12] = sad_16x8[12];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[12]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[13] = p_sad8x8[26][search_index] + p_sad8x8[27][search_index];
        if (sad_16x8[13] < p_best_sad_16x8[13]) {
            p_best_sad_16x8[13] = sad_16x8[13];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[13]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[14] = p_sad8x8[28][search_index] + p_sad8x8[29][search_index];
        if (sad_16x8[14] < p_best_sad_16x8[14]) {
            p_best_sad_16x8[14] = sad_16x8[14];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[14]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[15] = p_sad8x8[30][search_index] + p_sad8x8[31][search_index];
        if (sad_16x8[15] < p_best_sad_16x8[15]) {
            p_best_sad_16x8[15] = sad_16x8[15];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[15]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[16] = p_sad8x8[32][search_index] + p_sad8x8[33][search_index];
        if (sad_16x8[16] < p_best_sad_16x8[16]) {
            p_best_sad_16x8[16] = sad_16x8[16];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[16]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[17] = p_sad8x8[34][search_index] + p_sad8x8[35][search_index];
        if (sad_16x8[17] < p_best_sad_16x8[17]) {
            p_best_sad_16x8[17] = sad_16x8[17];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[17]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[18] = p_sad8x8[36][search_index] + p_sad8x8[37][search_index];
        if (sad_16x8[18] < p_best_sad_16x8[18]) {
            p_best_sad_16x8[18] = sad_16x8[18];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[18]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[19] = p_sad8x8[38][search_index] + p_sad8x8[39][search_index];
        if (sad_16x8[19] < p_best_sad_16x8[19]) {
            p_best_sad_16x8[19] = sad_16x8[19];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[19]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[20] = p_sad8x8[40][search_index] + p_sad8x8[41][search_index];
        if (sad_16x8[20] < p_best_sad_16x8[20]) {
            p_best_sad_16x8[20] = sad_16x8[20];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[20]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[21] = p_sad8x8[42][search_index] + p_sad8x8[43][search_index];
        if (sad_16x8[21] < p_best_sad_16x8[21]) {
            p_best_sad_16x8[21] = sad_16x8[21];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[21]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[22] = p_sad8x8[44][search_index] + p_sad8x8[45][search_index];
        if (sad_16x8[22] < p_best_sad_16x8[22]) {
            p_best_sad_16x8[22] = sad_16x8[22];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[22]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[23] = p_sad8x8[46][search_index] + p_sad8x8[47][search_index];
        if (sad_16x8[23] < p_best_sad_16x8[23]) {
            p_best_sad_16x8[23] = sad_16x8[23];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[23]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[24] = p_sad8x8[48][search_index] + p_sad8x8[49][search_index];
        if (sad_16x8[24] < p_best_sad_16x8[24]) {
            p_best_sad_16x8[24] = sad_16x8[24];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[24]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[25] = p_sad8x8[50][search_index] + p_sad8x8[51][search_index];
        if (sad_16x8[25] < p_best_sad_16x8[25]) {
            p_best_sad_16x8[25] = sad_16x8[25];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[25]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[26] = p_sad8x8[52][search_index] + p_sad8x8[53][search_index];
        if (sad_16x8[26] < p_best_sad_16x8[26]) {
            p_best_sad_16x8[26] = sad_16x8[26];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[26]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[27] = p_sad8x8[54][search_index] + p_sad8x8[55][search_index];
        if (sad_16x8[27] < p_best_sad_16x8[27]) {
            p_best_sad_16x8[27] = sad_16x8[27];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[27]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[28] = p_sad8x8[56][search_index] + p_sad8x8[57][search_index];
        if (sad_16x8[28] < p_best_sad_16x8[28]) {
            p_best_sad_16x8[28] = sad_16x8[28];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[28]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[29] = p_sad8x8[58][search_index] + p_sad8x8[59][search_index];
        if (sad_16x8[29] < p_best_sad_16x8[29]) {
            p_best_sad_16x8[29] = sad_16x8[29];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[29]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[30] = p_sad8x8[60][search_index] + p_sad8x8[61][search_index];
        if (sad_16x8[30] < p_best_sad_16x8[30]) {
            p_best_sad_16x8[30] = sad_16x8[30];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[30]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[31] = p_sad8x8[62][search_index] + p_sad8x8[63][search_index];
        if (sad_16x8[31] < p_best_sad_16x8[31]) {
            p_best_sad_16x8[31] = sad_16x8[31];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x8[31]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x64
        sad = p_sad32x32[0][search_index] + p_sad32x32[2][search_index];
        if (sad < p_best_sad_32x64[0]) {
            p_best_sad_32x64[0] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x64[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = p_sad32x32[1][search_index] + p_sad32x32[3][search_index];
        if (sad < p_best_sad_32x64[1]) {
            p_best_sad_32x64[1] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x64[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 16x32
        sad_16x32[0] = p_sad16x16[0][search_index] + p_sad16x16[2][search_index];
        if (sad_16x32[0] < p_best_sad_16x32[0]) {
            p_best_sad_16x32[0] = sad_16x32[0];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x32[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[1] = p_sad16x16[1][search_index] + p_sad16x16[3][search_index];
        if (sad_16x32[1] < p_best_sad_16x32[1]) {
            p_best_sad_16x32[1] = sad_16x32[1];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x32[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[2] = p_sad16x16[4][search_index] + p_sad16x16[6][search_index];
        if (sad_16x32[2] < p_best_sad_16x32[2]) {
            p_best_sad_16x32[2] = sad_16x32[2];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x32[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[3] = p_sad16x16[5][search_index] + p_sad16x16[7][search_index];
        if (sad_16x32[3] < p_best_sad_16x32[3]) {
            p_best_sad_16x32[3] = sad_16x32[3];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x32[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[4] = p_sad16x16[8][search_index] + p_sad16x16[10][search_index];
        if (sad_16x32[4] < p_best_sad_16x32[4]) {
            p_best_sad_16x32[4] = sad_16x32[4];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x32[4]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[5] = p_sad16x16[9][search_index] + p_sad16x16[11][search_index];
        if (sad_16x32[5] < p_best_sad_16x32[5]) {
            p_best_sad_16x32[5] = sad_16x32[5];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x32[5]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[6] = p_sad16x16[12][search_index] + p_sad16x16[14][search_index];
        if (sad_16x32[6] < p_best_sad_16x32[6]) {
            p_best_sad_16x32[6] = sad_16x32[6];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x32[6]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[7] = p_sad16x16[13][search_index] + p_sad16x16[15][search_index];
        if (sad_16x32[7] < p_best_sad_16x32[7]) {
            p_best_sad_16x32[7] = sad_16x32[7];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x32[7]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x32[0] + sad_16x32[4];
        if (sad < p_best_sad_16x64[0]) {
            p_best_sad_16x64[0] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x64[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        sad = sad_16x32[1] + sad_16x32[5];
        if (sad < p_best_sad_16x64[1]) {
            p_best_sad_16x64[1] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x64[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x32[2] + sad_16x32[6];
        if (sad < p_best_sad_16x64[2]) {
            p_best_sad_16x64[2] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x64[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x32[3] + sad_16x32[7];
        if (sad < p_best_sad_16x64[3]) {
            p_best_sad_16x64[3] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x64[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        // 8x16
        sad_8x16[0] = p_sad8x8[0][search_index] + p_sad8x8[2][search_index];
        if (sad_8x16[0] < p_best_sad_8x16[0]) {
            p_best_sad_8x16[0] = sad_8x16[0];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[1] = p_sad8x8[1][search_index] + p_sad8x8[3][search_index];
        if (sad_8x16[1] < p_best_sad_8x16[1]) {
            p_best_sad_8x16[1] = sad_8x16[1];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[2] = p_sad8x8[4][search_index] + p_sad8x8[6][search_index];
        if (sad_8x16[2] < p_best_sad_8x16[2]) {
            p_best_sad_8x16[2] = sad_8x16[2];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[3] = p_sad8x8[5][search_index] + p_sad8x8[7][search_index];
        if (sad_8x16[3] < p_best_sad_8x16[3]) {
            p_best_sad_8x16[3] = sad_8x16[3];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[4] = p_sad8x8[8][search_index] + p_sad8x8[10][search_index];
        if (sad_8x16[4] < p_best_sad_8x16[4]) {
            p_best_sad_8x16[4] = sad_8x16[4];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[4]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[5] = p_sad8x8[9][search_index] + p_sad8x8[11][search_index];
        if (sad_8x16[5] < p_best_sad_8x16[5]) {
            p_best_sad_8x16[5] = sad_8x16[5];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[5]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[6] = p_sad8x8[12][search_index] + p_sad8x8[14][search_index];
        if (sad_8x16[6] < p_best_sad_8x16[6]) {
            p_best_sad_8x16[6] = sad_8x16[6];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[6]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[7] = p_sad8x8[13][search_index] + p_sad8x8[15][search_index];
        if (sad_8x16[7] < p_best_sad_8x16[7]) {
            p_best_sad_8x16[7] = sad_8x16[7];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[7]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[8] = p_sad8x8[16][search_index] + p_sad8x8[18][search_index];
        if (sad_8x16[8] < p_best_sad_8x16[8]) {
            p_best_sad_8x16[8] = sad_8x16[8];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[8]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[9] = p_sad8x8[17][search_index] + p_sad8x8[19][search_index];
        if (sad_8x16[9] < p_best_sad_8x16[9]) {
            p_best_sad_8x16[9] = sad_8x16[9];
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x16[9]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[10] = p_sad8x8[20][search_index] + p_sad8x8[22][search_index];
        if (sad_8x16[10] < p_best_sad_8x16[10]) {
            p_best_sad_8x16[10] = sad_8x16[10];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[10]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[11] = p_sad8x8[21][search_index] + p_sad8x8[23][search_index];
        if (sad_8x16[11] < p_best_sad_8x16[11]) {
            p_best_sad_8x16[11] = sad_8x16[11];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[11]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[12] = p_sad8x8[24][search_index] + p_sad8x8[26][search_index];
        if (sad_8x16[12] < p_best_sad_8x16[12]) {
            p_best_sad_8x16[12] = sad_8x16[12];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[12]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[13] = p_sad8x8[25][search_index] + p_sad8x8[27][search_index];
        if (sad_8x16[13] < p_best_sad_8x16[13]) {
            p_best_sad_8x16[13] = sad_8x16[13];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[13]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[14] = p_sad8x8[28][search_index] + p_sad8x8[30][search_index];
        if (sad_8x16[14] < p_best_sad_8x16[14]) {
            p_best_sad_8x16[14] = sad_8x16[14];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[14]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[15] = p_sad8x8[29][search_index] + p_sad8x8[31][search_index];
        if (sad_8x16[15] < p_best_sad_8x16[15]) {
            p_best_sad_8x16[15] = sad_8x16[15];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[15]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[16] = p_sad8x8[32][search_index] + p_sad8x8[34][search_index];
        if (sad_8x16[16] < p_best_sad_8x16[16]) {
            p_best_sad_8x16[16] = sad_8x16[16];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[16]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[17] = p_sad8x8[33][search_index] + p_sad8x8[35][search_index];
        if (sad_8x16[17] < p_best_sad_8x16[17]) {
            p_best_sad_8x16[17] = sad_8x16[17];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[17]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[18] = p_sad8x8[36][search_index] + p_sad8x8[38][search_index];
        if (sad_8x16[18] < p_best_sad_8x16[18]) {
            p_best_sad_8x16[18] = sad_8x16[18];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[18]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[19] = p_sad8x8[37][search_index] + p_sad8x8[39][search_index];
        if (sad_8x16[19] < p_best_sad_8x16[19]) {
            p_best_sad_8x16[19] = sad_8x16[19];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[19]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[20] = p_sad8x8[40][search_index] + p_sad8x8[42][search_index];
        if (sad_8x16[20] < p_best_sad_8x16[20]) {
            p_best_sad_8x16[20] = sad_8x16[20];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[20]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[21] = p_sad8x8[41][search_index] + p_sad8x8[43][search_index];
        if (sad_8x16[21] < p_best_sad_8x16[21]) {
            p_best_sad_8x16[21] = sad_8x16[21];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[21]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[22] = p_sad8x8[44][search_index] + p_sad8x8[46][search_index];
        if (sad_8x16[22] < p_best_sad_8x16[22]) {
            p_best_sad_8x16[22] = sad_8x16[22];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[22]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[23] = p_sad8x8[45][search_index] + p_sad8x8[47][search_index];
        if (sad_8x16[23] < p_best_sad_8x16[23]) {
            p_best_sad_8x16[23] = sad_8x16[23];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[23]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[24] = p_sad8x8[48][search_index] + p_sad8x8[50][search_index];
        if (sad_8x16[24] < p_best_sad_8x16[24]) {
            p_best_sad_8x16[24] = sad_8x16[24];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[24]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[25] = p_sad8x8[49][search_index] + p_sad8x8[51][search_index];
        if (sad_8x16[25] < p_best_sad_8x16[25]) {
            p_best_sad_8x16[25] = sad_8x16[25];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[25]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[26] = p_sad8x8[52][search_index] + p_sad8x8[54][search_index];
        if (sad_8x16[26] < p_best_sad_8x16[26]) {
            p_best_sad_8x16[26] = sad_8x16[26];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[26]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[27] = p_sad8x8[53][search_index] + p_sad8x8[55][search_index];
        if (sad_8x16[27] < p_best_sad_8x16[27]) {
            p_best_sad_8x16[27] = sad_8x16[27];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[27]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[28] = p_sad8x8[56][search_index] + p_sad8x8[58][search_index];
        if (sad_8x16[28] < p_best_sad_8x16[28]) {
            p_best_sad_8x16[28] = sad_8x16[28];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[28]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[29] = p_sad8x8[57][search_index] + p_sad8x8[59][search_index];
        if (sad_8x16[29] < p_best_sad_8x16[29]) {
            p_best_sad_8x16[29] = sad_8x16[29];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[29]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[30] = p_sad8x8[60][search_index] + p_sad8x8[62][search_index];
        if (sad_8x16[30] < p_best_sad_8x16[30]) {
            p_best_sad_8x16[30] = sad_8x16[30];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[30]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[31] = p_sad8x8[61][search_index] + p_sad8x8[63][search_index];
        if (sad_8x16[31] < p_best_sad_8x16[31]) {
            p_best_sad_8x16[31] = sad_8x16[31];
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x16[31]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x8
        sad = sad_16x8[0] + sad_16x8[2];
        if (sad < p_best_sad_32x8[0]) {
            p_best_sad_32x8[0] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[1] + sad_16x8[3];
        if (sad < p_best_sad_32x8[1]) {
            p_best_sad_32x8[1] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[4] + sad_16x8[6];
        if (sad < p_best_sad_32x8[2]) {
            p_best_sad_32x8[2] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[5] + sad_16x8[7];
        if (sad < p_best_sad_32x8[3]) {
            p_best_sad_32x8[3] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[8] + sad_16x8[10];
        if (sad < p_best_sad_32x8[4]) {
            p_best_sad_32x8[4] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[4]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[9] + sad_16x8[11];
        if (sad < p_best_sad_32x8[5]) {
            p_best_sad_32x8[5] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[5]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[12] + sad_16x8[14];
        if (sad < p_best_sad_32x8[6]) {
            p_best_sad_32x8[6] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[6]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[13] + sad_16x8[15];
        if (sad < p_best_sad_32x8[7]) {
            p_best_sad_32x8[7] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[7]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[16] + sad_16x8[18];
        if (sad < p_best_sad_32x8[8]) {
            p_best_sad_32x8[8] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[8]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[17] + sad_16x8[19];
        if (sad < p_best_sad_32x8[9]) {
            p_best_sad_32x8[9] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv32x8[9]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[20] + sad_16x8[22];
        if (sad < p_best_sad_32x8[10]) {
            p_best_sad_32x8[10] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x8[10]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[21] + sad_16x8[23];
        if (sad < p_best_sad_32x8[11]) {
            p_best_sad_32x8[11] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x8[11]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[24] + sad_16x8[26];
        if (sad < p_best_sad_32x8[12]) {
            p_best_sad_32x8[12] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x8[12]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[25] + sad_16x8[27];
        if (sad < p_best_sad_32x8[13]) {
            p_best_sad_32x8[13] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x8[13]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[28] + sad_16x8[30];
        if (sad < p_best_sad_32x8[14]) {
            p_best_sad_32x8[14] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x8[14]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[29] + sad_16x8[31];
        if (sad < p_best_sad_32x8[15]) {
            p_best_sad_32x8[15] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x8[15]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        // 8x32
        sad = sad_8x16[0] + sad_8x16[4];
        if (sad < p_best_sad_8x32[0]) {
            p_best_sad_8x32[0] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[1] + sad_8x16[5];
        if (sad < p_best_sad_8x32[1]) {
            p_best_sad_8x32[1] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[2] + sad_8x16[6];
        if (sad < p_best_sad_8x32[2]) {
            p_best_sad_8x32[2] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[3] + sad_8x16[7];
        if (sad < p_best_sad_8x32[3]) {
            p_best_sad_8x32[3] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[8] + sad_8x16[12];
        if (sad < p_best_sad_8x32[4]) {
            p_best_sad_8x32[4] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[4]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[9] + sad_8x16[13];
        if (sad < p_best_sad_8x32[5]) {
            p_best_sad_8x32[5] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[5]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[10] + sad_8x16[14];
        if (sad < p_best_sad_8x32[6]) {
            p_best_sad_8x32[6] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[6]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[11] + sad_8x16[15];
        if (sad < p_best_sad_8x32[7]) {
            p_best_sad_8x32[7] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[7]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[16] + sad_8x16[20];
        if (sad < p_best_sad_8x32[8]) {
            p_best_sad_8x32[8] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[8]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[17] + sad_8x16[21];
        if (sad < p_best_sad_8x32[9]) {
            p_best_sad_8x32[9] = sad;
            x_mv               = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv               = _MVYT(mv);
            p_best_mv8x32[9]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[18] + sad_8x16[22];
        if (sad < p_best_sad_8x32[10]) {
            p_best_sad_8x32[10] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x32[10]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[19] + sad_8x16[23];
        if (sad < p_best_sad_8x32[11]) {
            p_best_sad_8x32[11] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x32[11]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[24] + sad_8x16[28];
        if (sad < p_best_sad_8x32[12]) {
            p_best_sad_8x32[12] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x32[12]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[25] + sad_8x16[29];
        if (sad < p_best_sad_8x32[13]) {
            p_best_sad_8x32[13] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x32[13]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[26] + sad_8x16[30];
        if (sad < p_best_sad_8x32[14]) {
            p_best_sad_8x32[14] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x32[14]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[27] + sad_8x16[31];
        if (sad < p_best_sad_8x32[15]) {
            p_best_sad_8x32[15] = sad;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv8x32[15]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }
}

/*******************************************
 * ext_eight_sad_calculation_8x8_16x16
 *******************************************/
static void ext_eight_sad_calculation_8x8_16x16(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                                uint32_t ref_stride, uint32_t mv,
                                                uint32_t start_16x16_pos, uint32_t *p_best_sad_8x8,
                                                uint32_t *p_best_sad_16x16, uint32_t *p_best_mv8x8,
                                                uint32_t *p_best_mv16x16,
                                                uint32_t  p_eight_sad16x16[16][8],
                                                uint32_t  p_eight_sad8x8[64][8]) {
    const uint32_t start_8x8_pos = 4 * start_16x16_pos;
    uint32_t       sad8x8_0, sad8x8_1, sad8x8_2, sad8x8_3;
    uint32_t       sad16x16;
    uint32_t       search_index;
    int16_t        x_mv, y_mv;
    uint32_t       src_stride_sub = (src_stride << 1);
    uint32_t       ref_stride_sub = (ref_stride << 1);

    p_best_sad_8x8 += start_8x8_pos;
    p_best_mv8x8 += start_8x8_pos;
    p_best_sad_16x16 += start_16x16_pos;
    p_best_mv16x16 += start_16x16_pos;

    for (search_index = 0; search_index < 8; search_index++) {
        p_eight_sad8x8[0 + start_8x8_pos][search_index] = sad8x8_0 =
            (compute8x4_sad_kernel_c(src, src_stride_sub, ref + search_index, ref_stride_sub)) << 1;
        if (sad8x8_0 < p_best_sad_8x8[0]) {
            p_best_sad_8x8[0] = (uint32_t)sad8x8_0;
            x_mv              = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv              = _MVYT(mv);
            p_best_mv8x8[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_eight_sad8x8[1 + start_8x8_pos][search_index] = sad8x8_1 =
            (compute8x4_sad_kernel_c(
                src + 8, src_stride_sub, ref + 8 + search_index, ref_stride_sub))
            << 1;
        if (sad8x8_1 < p_best_sad_8x8[1]) {
            p_best_sad_8x8[1] = (uint32_t)sad8x8_1;
            x_mv              = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv              = _MVYT(mv);
            p_best_mv8x8[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_eight_sad8x8[2 + start_8x8_pos][search_index] = sad8x8_2 =
            (compute8x4_sad_kernel_c(src + (src_stride << 3),
                                     src_stride_sub,
                                     ref + (ref_stride << 3) + search_index,
                                     ref_stride_sub))
            << 1;
        if (sad8x8_2 < p_best_sad_8x8[2]) {
            p_best_sad_8x8[2] = (uint32_t)sad8x8_2;
            x_mv              = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv              = _MVYT(mv);
            p_best_mv8x8[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_eight_sad8x8[3 + start_8x8_pos][search_index] = sad8x8_3 =
            (compute8x4_sad_kernel_c(src + (src_stride << 3) + 8,
                                     src_stride_sub,
                                     ref + (ref_stride << 3) + 8 + search_index,
                                     ref_stride_sub))
            << 1;
        if (sad8x8_3 < p_best_sad_8x8[3]) {
            p_best_sad_8x8[3] = (uint32_t)sad8x8_3;
            x_mv              = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv              = _MVYT(mv);
            p_best_mv8x8[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_eight_sad16x16[start_16x16_pos][search_index] = sad16x16 =
            sad8x8_0 + sad8x8_1 + sad8x8_2 + sad8x8_3;
        if (sad16x16 < p_best_sad_16x16[0]) {
            p_best_sad_16x16[0] = (uint32_t)sad16x16;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv16x16[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }
}

void ext_all_sad_calculation_8x8_16x16_c(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                                         uint32_t ref_stride, uint32_t mv, uint32_t *p_best_sad_8x8,
                                         uint32_t *p_best_sad_16x16, uint32_t *p_best_mv8x8,
                                         uint32_t *p_best_mv16x16, uint32_t p_eight_sad16x16[16][8],
                                         uint32_t p_eight_sad8x8[64][8]) {
    static const char offsets[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

    //---- 16x16 : 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            const uint32_t block_index           = 16 * y * src_stride + 16 * x;
            const uint32_t search_position_index = 16 * y * ref_stride + 16 * x;
            ext_eight_sad_calculation_8x8_16x16(src + block_index,
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
                                                p_eight_sad8x8);
        }
    }
}

/*******************************************
Calculate SAD for 32x32,64x64 from 16x16
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void ext_eight_sad_calculation_32x32_64x64_c(uint32_t p_sad16x16[16][8], uint32_t *p_best_sad_32x32,
                                             uint32_t *p_best_sad_64x64, uint32_t *p_best_mv32x32,
                                             uint32_t *p_best_mv64x64, uint32_t mv,
                                             uint32_t p_sad32x32[4][8]) {
    uint32_t search_index;
    int16_t  x_mv, y_mv;
    for (search_index = 0; search_index < 8; search_index++) {
        uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;

        p_sad32x32[0][search_index] = sad32x32_0 =
            p_sad16x16[0][search_index] + p_sad16x16[1][search_index] +
            p_sad16x16[2][search_index] + p_sad16x16[3][search_index];
        if (sad32x32_0 < p_best_sad_32x32[0]) {
            p_best_sad_32x32[0] = sad32x32_0;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x32[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_sad32x32[1][search_index] = sad32x32_1 =
            p_sad16x16[4][search_index] + p_sad16x16[5][search_index] +
            p_sad16x16[6][search_index] + p_sad16x16[7][search_index];
        if (sad32x32_1 < p_best_sad_32x32[1]) {
            p_best_sad_32x32[1] = sad32x32_1;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x32[1]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_sad32x32[2][search_index] = sad32x32_2 =
            p_sad16x16[8][search_index] + p_sad16x16[9][search_index] +
            p_sad16x16[10][search_index] + p_sad16x16[11][search_index];
        if (sad32x32_2 < p_best_sad_32x32[2]) {
            p_best_sad_32x32[2] = sad32x32_2;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x32[2]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_sad32x32[3][search_index] = sad32x32_3 =
            p_sad16x16[12][search_index] + p_sad16x16[13][search_index] +
            p_sad16x16[14][search_index] + p_sad16x16[15][search_index];
        if (sad32x32_3 < p_best_sad_32x32[3]) {
            p_best_sad_32x32[3] = sad32x32_3;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv32x32[3]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
        if (sad64x64 < p_best_sad_64x64[0]) {
            p_best_sad_64x64[0] = sad64x64;
            x_mv                = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv                = _MVYT(mv);
            p_best_mv64x64[0]   = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }
}
/*******************************************
 * inherit nsq MVs from SQ MVs
 *******************************************/
void generate_nsq_mv(MeContext *context_ptr) {
    uint32_t *p_sad8x8 = context_ptr->p_best_sad_8x8;
    uint32_t *p_sad16x16 = context_ptr->p_best_sad_16x16;
    uint32_t *p_sad32x32 = context_ptr->p_best_sad_32x32;
    uint32_t *p_best_mv8x8 = context_ptr->p_best_mv8x8;
    uint32_t *p_best_mv16x16 = context_ptr->p_best_mv16x16;
    uint32_t *p_best_mv32x32 = context_ptr->p_best_mv32x32;
    uint32_t *p_best_mv64x32 = context_ptr->p_best_mv64x32;
    uint32_t *p_best_mv32x16 = context_ptr->p_best_mv32x16;
    uint32_t *p_best_mv16x8 = context_ptr->p_best_mv16x8;
    uint32_t *p_best_mv32x64 = context_ptr->p_best_mv32x64;
    uint32_t *p_best_mv16x32 = context_ptr->p_best_mv16x32;
    uint32_t *p_best_mv8x16 = context_ptr->p_best_mv8x16;
    uint32_t *p_best_mv32x8 = context_ptr->p_best_mv32x8;
    uint32_t *p_best_mv8x32 = context_ptr->p_best_mv8x32;
    uint32_t *p_best_mv64x16 = context_ptr->p_best_mv64x16;
    uint32_t *p_best_mv16x64 = context_ptr->p_best_mv16x64;
    // 64x32
    p_best_mv64x32[0] = p_sad32x32[0] < p_sad32x32[1] ? p_best_mv32x32[0] : p_best_mv32x32[1];
    p_best_mv64x32[1] = p_sad32x32[2] < p_sad32x32[3] ? p_best_mv32x32[2] : p_best_mv32x32[3];
    // 32x16
    p_best_mv32x16[0] =p_sad16x16[0] < p_sad16x16[1] ? p_best_mv16x16[0] : p_best_mv16x16[1];
    p_best_mv32x16[1] =p_sad16x16[2] < p_sad16x16[3] ? p_best_mv16x16[2] : p_best_mv16x16[3];
    p_best_mv32x16[2] =p_sad16x16[4] < p_sad16x16[5] ? p_best_mv16x16[4] : p_best_mv16x16[5];
    p_best_mv32x16[3] =p_sad16x16[6] < p_sad16x16[7] ? p_best_mv16x16[6] : p_best_mv16x16[7];
    p_best_mv32x16[4] =p_sad16x16[8] < p_sad16x16[9] ? p_best_mv16x16[8] : p_best_mv16x16[9];
    p_best_mv32x16[5] =p_sad16x16[10] < p_sad16x16[11] ? p_best_mv16x16[10] : p_best_mv16x16[11];
    p_best_mv32x16[6] =p_sad16x16[12] < p_sad16x16[13] ? p_best_mv16x16[12] : p_best_mv16x16[13];
    p_best_mv32x16[7] =p_sad16x16[14] < p_sad16x16[15] ? p_best_mv16x16[14] : p_best_mv16x16[15];
    // 64x16
    p_best_mv64x16[0] =(p_sad16x16[0] + p_sad16x16[1]) < (p_sad16x16[4] + p_sad16x16[5])  ? p_best_mv32x16[0] : p_best_mv32x16[2];
    p_best_mv64x16[1] =(p_sad16x16[2] + p_sad16x16[3]) < (p_sad16x16[6] + p_sad16x16[7])  ? p_best_mv32x16[1] : p_best_mv32x16[3];
    p_best_mv64x16[2] =(p_sad16x16[8] + p_sad16x16[9]) < (p_sad16x16[12] + p_sad16x16[13])  ? p_best_mv32x16[4] : p_best_mv32x16[6];
    p_best_mv64x16[3] =(p_sad16x16[10] + p_sad16x16[11]) < (p_sad16x16[14] + p_sad16x16[15])  ? p_best_mv32x16[5] : p_best_mv32x16[7];
    // 16x8
    p_best_mv16x8[0] =p_sad8x8[0] < p_sad8x8[1]  ? p_best_mv8x8[0] : p_best_mv8x8[1];
    p_best_mv16x8[1] =p_sad8x8[2] < p_sad8x8[3]  ? p_best_mv8x8[2] : p_best_mv8x8[3];
    p_best_mv16x8[2] =p_sad8x8[4] < p_sad8x8[5]  ? p_best_mv8x8[4] : p_best_mv8x8[5];
    p_best_mv16x8[3] =p_sad8x8[6] < p_sad8x8[7]  ? p_best_mv8x8[6] : p_best_mv8x8[7];
    p_best_mv16x8[4] =p_sad8x8[8] < p_sad8x8[9]  ? p_best_mv8x8[8] : p_best_mv8x8[9];
    p_best_mv16x8[5] =p_sad8x8[10] < p_sad8x8[11] ? p_best_mv8x8[10] : p_best_mv8x8[11];
    p_best_mv16x8[6] =p_sad8x8[12] < p_sad8x8[13] ? p_best_mv8x8[12] : p_best_mv8x8[13];
    p_best_mv16x8[7] =p_sad8x8[14] < p_sad8x8[15] ? p_best_mv8x8[14] : p_best_mv8x8[15];
    p_best_mv16x8[8] =p_sad8x8[16] < p_sad8x8[17] ? p_best_mv8x8[16] : p_best_mv8x8[17];
    p_best_mv16x8[9] =p_sad8x8[18] < p_sad8x8[19] ? p_best_mv8x8[18] : p_best_mv8x8[19];
    p_best_mv16x8[10] =p_sad8x8[20] < p_sad8x8[21] ? p_best_mv8x8[20] : p_best_mv8x8[21];
    p_best_mv16x8[11] =p_sad8x8[22] < p_sad8x8[23] ? p_best_mv8x8[22] : p_best_mv8x8[23];
    p_best_mv16x8[12] =p_sad8x8[24] < p_sad8x8[25] ? p_best_mv8x8[24] : p_best_mv8x8[25];
    p_best_mv16x8[13] =p_sad8x8[26] < p_sad8x8[27] ? p_best_mv8x8[26] : p_best_mv8x8[27];
    p_best_mv16x8[14] =p_sad8x8[28] < p_sad8x8[29] ? p_best_mv8x8[28] : p_best_mv8x8[29];
    p_best_mv16x8[15] =p_sad8x8[30] < p_sad8x8[31] ? p_best_mv8x8[30] : p_best_mv8x8[31];
    p_best_mv16x8[16] =p_sad8x8[32] < p_sad8x8[33] ? p_best_mv8x8[32] : p_best_mv8x8[33];
    p_best_mv16x8[17] =p_sad8x8[34] < p_sad8x8[35] ? p_best_mv8x8[34] : p_best_mv8x8[35];
    p_best_mv16x8[18] =p_sad8x8[36] < p_sad8x8[37] ? p_best_mv8x8[36] : p_best_mv8x8[37];
    p_best_mv16x8[19] =p_sad8x8[38] < p_sad8x8[39] ? p_best_mv8x8[38] : p_best_mv8x8[39];
    p_best_mv16x8[20] =p_sad8x8[40] < p_sad8x8[41] ? p_best_mv8x8[40] : p_best_mv8x8[41];
    p_best_mv16x8[21] =p_sad8x8[42] < p_sad8x8[43] ? p_best_mv8x8[42] : p_best_mv8x8[43];
    p_best_mv16x8[22] =p_sad8x8[44] < p_sad8x8[45] ? p_best_mv8x8[44] : p_best_mv8x8[45];
    p_best_mv16x8[23] =p_sad8x8[46] < p_sad8x8[47] ? p_best_mv8x8[46] : p_best_mv8x8[47];
    p_best_mv16x8[24] =p_sad8x8[48] < p_sad8x8[49] ? p_best_mv8x8[48] : p_best_mv8x8[49];
    p_best_mv16x8[25] =p_sad8x8[50] < p_sad8x8[51] ? p_best_mv8x8[50] : p_best_mv8x8[51];
    p_best_mv16x8[26] =p_sad8x8[52] < p_sad8x8[53] ? p_best_mv8x8[52] : p_best_mv8x8[53];
    p_best_mv16x8[27] =p_sad8x8[54] < p_sad8x8[55] ? p_best_mv8x8[54] : p_best_mv8x8[55];
    p_best_mv16x8[28] =p_sad8x8[56] < p_sad8x8[57] ? p_best_mv8x8[56] : p_best_mv8x8[57];
    p_best_mv16x8[29] =p_sad8x8[58] < p_sad8x8[59] ? p_best_mv8x8[58] : p_best_mv8x8[59];
    p_best_mv16x8[30] =p_sad8x8[60] < p_sad8x8[61] ? p_best_mv8x8[60] : p_best_mv8x8[61];
    p_best_mv16x8[31] =p_sad8x8[62] < p_sad8x8[63] ? p_best_mv8x8[62] : p_best_mv8x8[63];
    // 32x64
    p_best_mv32x64[0] = p_sad32x32[0] < p_sad32x32[2] ? p_best_mv32x32[0] : p_best_mv32x32[2];
    p_best_mv32x64[1] = p_sad32x32[1] < p_sad32x32[3] ? p_best_mv32x32[1] : p_best_mv32x32[3];
    // 16x32
    p_best_mv16x32[0] =p_sad16x16[0] < p_sad16x16[2] ? p_best_mv16x16[0] : p_best_mv16x16[2];
    p_best_mv16x32[1] =p_sad16x16[1] < p_sad16x16[3] ? p_best_mv16x16[1] : p_best_mv16x16[3];
    p_best_mv16x32[2] =p_sad16x16[4] < p_sad16x16[6] ? p_best_mv16x16[4] : p_best_mv16x16[6];
    p_best_mv16x32[3] =p_sad16x16[5] < p_sad16x16[7] ? p_best_mv16x16[5] : p_best_mv16x16[7];
    p_best_mv16x32[4] =p_sad16x16[8] < p_sad16x16[10] ? p_best_mv16x16[8] : p_best_mv16x16[10];
    p_best_mv16x32[5] =p_sad16x16[9] < p_sad16x16[11] ? p_best_mv16x16[9] : p_best_mv16x16[11];
    p_best_mv16x32[6] =p_sad16x16[12] < p_sad16x16[14] ? p_best_mv16x16[12] : p_best_mv16x16[14];
    p_best_mv16x32[7] =p_sad16x16[13] < p_sad16x16[15] ? p_best_mv16x16[13] : p_best_mv16x16[15];
    // 16x64
    p_best_mv16x64[0] =(p_sad16x16[0] + p_sad16x16[2]) < (p_sad16x16[8] + p_sad16x16[10])  ? p_best_mv32x16[0] : p_best_mv32x16[4];
    p_best_mv16x64[1] =(p_sad16x16[1] + p_sad16x16[3]) < (p_sad16x16[9] + p_sad16x16[11])  ? p_best_mv32x16[1] : p_best_mv32x16[5];
    p_best_mv16x64[2] =(p_sad16x16[4] + p_sad16x16[6]) < (p_sad16x16[12] + p_sad16x16[14])  ? p_best_mv32x16[2] : p_best_mv32x16[6];
    p_best_mv16x64[3] =(p_sad16x16[5] + p_sad16x16[7]) < (p_sad16x16[13] + p_sad16x16[15])  ? p_best_mv32x16[3] : p_best_mv32x16[7];
    // 8x16
    p_best_mv8x16[0] =p_sad8x8[0] < p_sad8x8[2]  ? p_best_mv8x8[0] : p_best_mv8x8[2];
    p_best_mv8x16[1] =p_sad8x8[1] < p_sad8x8[3]  ? p_best_mv8x8[1] : p_best_mv8x8[3];
    p_best_mv8x16[2] =p_sad8x8[4] < p_sad8x8[6]  ? p_best_mv8x8[4] : p_best_mv8x8[6];
    p_best_mv8x16[3] =p_sad8x8[5] < p_sad8x8[7]  ? p_best_mv8x8[5] : p_best_mv8x8[7];
    p_best_mv8x16[4] =p_sad8x8[8] < p_sad8x8[10]  ? p_best_mv8x8[8] : p_best_mv8x8[10];
    p_best_mv8x16[5] =p_sad8x8[9] < p_sad8x8[11]  ? p_best_mv8x8[9] : p_best_mv8x8[11];
    p_best_mv8x16[6] =p_sad8x8[12] < p_sad8x8[14]  ? p_best_mv8x8[12] : p_best_mv8x8[14];
    p_best_mv8x16[7] =p_sad8x8[13] < p_sad8x8[15]  ? p_best_mv8x8[13] : p_best_mv8x8[15];
    p_best_mv8x16[8] =p_sad8x8[16] < p_sad8x8[18]  ? p_best_mv8x8[16] : p_best_mv8x8[18];
    p_best_mv8x16[9] =p_sad8x8[17] < p_sad8x8[19]  ? p_best_mv8x8[17] : p_best_mv8x8[19];
    p_best_mv8x16[10] =p_sad8x8[20] < p_sad8x8[22]  ? p_best_mv8x8[20] : p_best_mv8x8[22];
    p_best_mv8x16[11] =p_sad8x8[21] < p_sad8x8[23]  ? p_best_mv8x8[21] : p_best_mv8x8[23];
    p_best_mv8x16[12] =p_sad8x8[24] < p_sad8x8[26]  ? p_best_mv8x8[24] : p_best_mv8x8[26];
    p_best_mv8x16[13] =p_sad8x8[25] < p_sad8x8[27]  ? p_best_mv8x8[25] : p_best_mv8x8[27];
    p_best_mv8x16[14] =p_sad8x8[28] < p_sad8x8[30]  ? p_best_mv8x8[28] : p_best_mv8x8[30];
    p_best_mv8x16[15] =p_sad8x8[29] < p_sad8x8[31]  ? p_best_mv8x8[29] : p_best_mv8x8[31];
    p_best_mv8x16[16] =p_sad8x8[32] < p_sad8x8[34]  ? p_best_mv8x8[32] : p_best_mv8x8[34];
    p_best_mv8x16[17] =p_sad8x8[33] < p_sad8x8[35]  ? p_best_mv8x8[33] : p_best_mv8x8[35];
    p_best_mv8x16[18] =p_sad8x8[36] < p_sad8x8[38]  ? p_best_mv8x8[36] : p_best_mv8x8[38];
    p_best_mv8x16[19] =p_sad8x8[37] < p_sad8x8[39]  ? p_best_mv8x8[37] : p_best_mv8x8[39];
    p_best_mv8x16[20] =p_sad8x8[40] < p_sad8x8[42]  ? p_best_mv8x8[40] : p_best_mv8x8[42];
    p_best_mv8x16[21] =p_sad8x8[41] < p_sad8x8[43]  ? p_best_mv8x8[41] : p_best_mv8x8[43];
    p_best_mv8x16[22] =p_sad8x8[44] < p_sad8x8[46]  ? p_best_mv8x8[44] : p_best_mv8x8[46];
    p_best_mv8x16[23] =p_sad8x8[45] < p_sad8x8[47]  ? p_best_mv8x8[45] : p_best_mv8x8[47];
    p_best_mv8x16[24] =p_sad8x8[48] < p_sad8x8[50]  ? p_best_mv8x8[48] : p_best_mv8x8[50];
    p_best_mv8x16[25] =p_sad8x8[49] < p_sad8x8[51]  ? p_best_mv8x8[49] : p_best_mv8x8[51];
    p_best_mv8x16[26] =p_sad8x8[52] < p_sad8x8[54]  ? p_best_mv8x8[52] : p_best_mv8x8[54];
    p_best_mv8x16[27] =p_sad8x8[53] < p_sad8x8[55]  ? p_best_mv8x8[53] : p_best_mv8x8[55];
    p_best_mv8x16[28] =p_sad8x8[56] < p_sad8x8[58]  ? p_best_mv8x8[56] : p_best_mv8x8[58];
    p_best_mv8x16[29] =p_sad8x8[57] < p_sad8x8[59]  ? p_best_mv8x8[57] : p_best_mv8x8[59];
    p_best_mv8x16[30] =p_sad8x8[60] < p_sad8x8[62]  ? p_best_mv8x8[60] : p_best_mv8x8[62];
    p_best_mv8x16[31] =p_sad8x8[61] < p_sad8x8[63]  ? p_best_mv8x8[61] : p_best_mv8x8[63];
    // 32x8
    p_best_mv32x8[0] =(p_sad8x8[0] + p_sad8x8[1]) < (p_sad16x16[4] + p_sad8x8[5])  ? p_best_mv16x8[0] : p_best_mv16x8[2];
    p_best_mv32x8[1] =(p_sad8x8[2] + p_sad8x8[3]) < (p_sad16x16[6] + p_sad8x8[7])  ? p_best_mv16x8[1] : p_best_mv16x8[3];
    p_best_mv32x8[2] =(p_sad8x8[8] + p_sad8x8[9]) < (p_sad16x16[12] + p_sad8x8[13])  ? p_best_mv16x8[4] : p_best_mv16x8[6];
    p_best_mv32x8[3] =(p_sad8x8[10] + p_sad8x8[11]) < (p_sad16x16[14] + p_sad8x8[15])  ? p_best_mv16x8[5] : p_best_mv16x8[7];
    p_best_mv32x8[4] =(p_sad8x8[16] + p_sad8x8[17]) < (p_sad16x16[20] + p_sad8x8[21])  ? p_best_mv16x8[8] : p_best_mv16x8[10];
    p_best_mv32x8[5] =(p_sad8x8[18] + p_sad8x8[19]) < (p_sad16x16[22] + p_sad8x8[23])  ? p_best_mv16x8[9] : p_best_mv16x8[11];
    p_best_mv32x8[6] =(p_sad8x8[24] + p_sad8x8[25]) < (p_sad16x16[28] + p_sad8x8[29])  ? p_best_mv16x8[12] : p_best_mv16x8[14];
    p_best_mv32x8[7] =(p_sad8x8[26] + p_sad8x8[27]) < (p_sad16x16[30] + p_sad8x8[31])  ? p_best_mv16x8[13] : p_best_mv16x8[15];
    p_best_mv32x8[8] =(p_sad8x8[32] + p_sad8x8[33]) < (p_sad16x16[36] + p_sad8x8[37])  ? p_best_mv16x8[16] : p_best_mv16x8[18];
    p_best_mv32x8[9] =(p_sad8x8[34] + p_sad8x8[36]) < (p_sad16x16[38] + p_sad8x8[39])  ? p_best_mv16x8[17] : p_best_mv16x8[19];
    p_best_mv32x8[10] =(p_sad8x8[40] + p_sad8x8[42]) < (p_sad16x16[44] + p_sad8x8[45])  ? p_best_mv16x8[20] : p_best_mv16x8[22];
    p_best_mv32x8[11] =(p_sad8x8[41] + p_sad8x8[43]) < (p_sad16x16[46] + p_sad8x8[47])  ? p_best_mv16x8[21] : p_best_mv16x8[23];
    p_best_mv32x8[12] =(p_sad8x8[48] + p_sad8x8[49]) < (p_sad16x16[52] + p_sad8x8[53])  ? p_best_mv16x8[24] : p_best_mv16x8[26];
    p_best_mv32x8[13] =(p_sad8x8[50] + p_sad8x8[51]) < (p_sad16x16[54] + p_sad8x8[55])  ? p_best_mv16x8[25] : p_best_mv16x8[27];
    p_best_mv32x8[14] =(p_sad8x8[56] + p_sad8x8[57]) < (p_sad16x16[60] + p_sad8x8[61])  ? p_best_mv16x8[28] : p_best_mv16x8[30];
    p_best_mv32x8[15] =(p_sad8x8[58] + p_sad8x8[59]) < (p_sad16x16[66] + p_sad8x8[63])  ? p_best_mv16x8[29] : p_best_mv16x8[31];
    // 8x32
    p_best_mv8x32[0] =(p_sad8x8[0] + p_sad8x8[2]) < (p_sad16x16[8] + p_sad8x8[10])  ? p_best_mv16x8[0] : p_best_mv16x8[4];
    p_best_mv8x32[1] =(p_sad8x8[1] + p_sad8x8[3]) < (p_sad16x16[9] + p_sad8x8[11])  ? p_best_mv16x8[1] : p_best_mv16x8[5];
    p_best_mv8x32[2] =(p_sad8x8[4] + p_sad8x8[6]) < (p_sad16x16[12] + p_sad8x8[14])  ? p_best_mv16x8[2] : p_best_mv16x8[6];
    p_best_mv8x32[3] =(p_sad8x8[5] + p_sad8x8[7]) < (p_sad16x16[13] + p_sad8x8[15])  ? p_best_mv16x8[3] : p_best_mv16x8[7];
    p_best_mv8x32[4] =(p_sad8x8[16] + p_sad8x8[18]) < (p_sad16x16[24] + p_sad8x8[26])  ? p_best_mv16x8[8] : p_best_mv16x8[12];
    p_best_mv8x32[5] =(p_sad8x8[17] + p_sad8x8[19]) < (p_sad16x16[25] + p_sad8x8[27])  ? p_best_mv16x8[9] : p_best_mv16x8[13];
    p_best_mv8x32[6] =(p_sad8x8[20] + p_sad8x8[22]) < (p_sad16x16[28] + p_sad8x8[30])  ? p_best_mv16x8[10] : p_best_mv16x8[14];
    p_best_mv8x32[7] =(p_sad8x8[21] + p_sad8x8[23]) < (p_sad16x16[29] + p_sad8x8[31])  ? p_best_mv16x8[11] : p_best_mv16x8[15];
    p_best_mv8x32[8] =(p_sad8x8[32] + p_sad8x8[34]) < (p_sad16x16[40] + p_sad8x8[42])  ? p_best_mv16x8[16] : p_best_mv16x8[20];
    p_best_mv8x32[9] =(p_sad8x8[33] + p_sad8x8[35]) < (p_sad16x16[41] + p_sad8x8[43])  ? p_best_mv16x8[17] : p_best_mv16x8[21];
    p_best_mv8x32[10] =(p_sad8x8[36] + p_sad8x8[38]) < (p_sad16x16[44] + p_sad8x8[46])  ? p_best_mv16x8[18] : p_best_mv16x8[22];
    p_best_mv8x32[11] =(p_sad8x8[37] + p_sad8x8[39]) < (p_sad16x16[45] + p_sad8x8[47])  ? p_best_mv16x8[19] : p_best_mv16x8[23];
    p_best_mv8x32[12] =(p_sad8x8[48] + p_sad8x8[50]) < (p_sad16x16[56] + p_sad8x8[58])  ? p_best_mv16x8[24] : p_best_mv16x8[28];
    p_best_mv8x32[13] =(p_sad8x8[49] + p_sad8x8[51]) < (p_sad16x16[57] + p_sad8x8[59])  ? p_best_mv16x8[25] : p_best_mv16x8[29];
    p_best_mv8x32[14] =(p_sad8x8[52] + p_sad8x8[54]) < (p_sad16x16[60] + p_sad8x8[62])  ? p_best_mv16x8[26] : p_best_mv16x8[30];
    p_best_mv8x32[15] =(p_sad8x8[53] + p_sad8x8[55]) < (p_sad16x16[61] + p_sad8x8[63])  ? p_best_mv16x8[27] : p_best_mv16x8[31];
}
/*******************************************
 * open_loop_me_get_search_point_results_block
 *******************************************/
static void open_loop_me_get_eight_search_point_results_block(
    MeContext *context_ptr, // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t   list_index, // input parameter, reference list index
    uint32_t   ref_pic_index,
    uint32_t   search_region_index, // input parameter, search area origin, used to
    // point to reference samples
    int32_t x_search_index, // input parameter, search region position in the
    // horizontal direction, used to derive xMV
    int32_t y_search_index // input parameter, search region position in the
    // vertical direction, used to derive yMV
) {
    // uint32_t ref_luma_stride = ref_pic_ptr->stride_y; // NADER
    // uint8_t  *ref_ptr = ref_pic_ptr->buffer_y; // NADER
    uint32_t ref_luma_stride = context_ptr->interpolated_full_stride[list_index][ref_pic_index];
    uint8_t *ref_ptr =
        context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
        ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[list_index][ref_pic_index]) +
        (ME_FILTER_TAP >> 1) + search_region_index;

    uint32_t curr_mv_1 = (((uint16_t)y_search_index) << 18);
    uint16_t curr_mv_2 = (((uint16_t)x_search_index << 2));
    uint32_t curr_mv   = curr_mv_1 | curr_mv_2;

    ext_all_sad_calculation_8x8_16x16(context_ptr->sb_src_ptr,
                                      context_ptr->sb_src_stride,
                                      ref_ptr,
                                      ref_luma_stride,
                                      curr_mv,
                                      context_ptr->p_best_sad_8x8,
                                      context_ptr->p_best_sad_16x16,
                                      context_ptr->p_best_mv8x8,
                                      context_ptr->p_best_mv16x16,
                                      context_ptr->p_eight_sad16x16,
                                      context_ptr->p_eight_sad8x8);

    ext_eight_sad_calculation_32x32_64x64(context_ptr->p_eight_sad16x16,
                                          context_ptr->p_best_sad_32x32,
                                          context_ptr->p_best_sad_64x64,
                                          context_ptr->p_best_mv32x32,
                                          context_ptr->p_best_mv64x64,
                                          curr_mv,
                                          context_ptr->p_eight_sad32x32);
    uint8_t perform_nsq_flag = 1;
    perform_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 1 && (list_index != context_ptr->best_list_idx || ref_pic_index != context_ptr->best_ref_idx)) ? 0 : perform_nsq_flag;
    perform_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 2 && ref_pic_index) ? 0: perform_nsq_flag;
    perform_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 3 ) ? 0 : perform_nsq_flag;
    if(perform_nsq_flag)
        ext_eigth_sad_calculation_nsq(context_ptr->p_eight_sad8x8,
                                  context_ptr->p_eight_sad16x16,
                                  context_ptr->p_eight_sad32x32,
                                  context_ptr->p_best_sad_64x32,
                                  context_ptr->p_best_mv64x32,
                                  context_ptr->p_best_sad_32x16,
                                  context_ptr->p_best_mv32x16,
                                  context_ptr->p_best_sad_16x8,
                                  context_ptr->p_best_mv16x8,
                                  context_ptr->p_best_sad_32x64,
                                  context_ptr->p_best_mv32x64,
                                  context_ptr->p_best_sad_16x32,
                                  context_ptr->p_best_mv16x32,
                                  context_ptr->p_best_sad_8x16,
                                  context_ptr->p_best_mv8x16,
                                  context_ptr->p_best_sad_32x8,
                                  context_ptr->p_best_mv32x8,
                                  context_ptr->p_best_sad_8x32,
                                  context_ptr->p_best_mv8x32,
                                  context_ptr->p_best_sad_64x16,
                                  context_ptr->p_best_mv64x16,
                                  context_ptr->p_best_sad_16x64,
                                  context_ptr->p_best_mv16x64,
                                  curr_mv);
}
/*******************************************
 * open_loop_me_get_search_point_results_block
 *******************************************/
static void open_loop_me_get_search_point_results_block(
    MeContext *context_ptr, // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t   list_index, // input parameter, reference list index
    uint32_t   ref_pic_index,
    uint32_t   search_region_index, // input parameter, search area origin, used to
    // point to reference samples
    int32_t x_search_index, // input parameter, search region position in the
    // horizontal direction, used to derive xMV
    int32_t y_search_index) // input parameter, search region position in the
// vertical direction, used to derive yMV
{
    const EbBool sub_sad = (context_ptr->me_search_method == SUB_SAD_SEARCH);
    uint8_t *    src_ptr = context_ptr->sb_src_ptr;

    // uint8_t  *ref_ptr = ref_pic_ptr->buffer_y; // NADER
    uint8_t *ref_ptr =
        context_ptr->integer_buffer_ptr[list_index][ref_pic_index] + (ME_FILTER_TAP >> 1) +
        ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[list_index][ref_pic_index]);
    // uint32_t ref_luma_stride = ref_pic_ptr->stride_y; // NADER
    uint32_t ref_luma_stride = context_ptr->interpolated_full_stride[list_index][ref_pic_index];
    uint32_t search_position_tl_index = search_region_index;
    uint32_t search_position_index;
    uint32_t block_index;
    uint32_t src_next_16x16_offset = (BLOCK_SIZE_64 << 4);
    // uint32_t ref_next_16x16_offset = (ref_pic_ptr->stride_y << 4); // NADER
    uint32_t  ref_next_16x16_offset = (ref_luma_stride << 4);
    uint32_t  curr_mv_1             = (((uint16_t)y_search_index) << 18);
    uint16_t  curr_mv_2             = (((uint16_t)x_search_index << 2));
    uint32_t  curr_mv               = curr_mv_1 | curr_mv_2;
    uint32_t *p_best_sad_8x8        = context_ptr->p_best_sad_8x8;
    uint32_t *p_best_sad_16x16      = context_ptr->p_best_sad_16x16;
    uint32_t *p_best_sad_32x32      = context_ptr->p_best_sad_32x32;
    uint32_t *p_best_sad_64x64      = context_ptr->p_best_sad_64x64;
    uint32_t *p_best_sad_64x32      = context_ptr->p_best_sad_64x32;
    uint32_t *p_best_sad_32x16      = context_ptr->p_best_sad_32x16;
    uint32_t *p_best_sad_16x8       = context_ptr->p_best_sad_16x8;
    uint32_t *p_best_sad_32x64      = context_ptr->p_best_sad_32x64;
    uint32_t *p_best_sad_16x32      = context_ptr->p_best_sad_16x32;
    uint32_t *p_best_sad_8x16       = context_ptr->p_best_sad_8x16;
    uint32_t *p_best_sad_32x8       = context_ptr->p_best_sad_32x8;
    uint32_t *p_best_sad_8x32       = context_ptr->p_best_sad_8x32;
    uint32_t *p_best_sad_64x16      = context_ptr->p_best_sad_64x16;
    uint32_t *p_best_sad_16x64      = context_ptr->p_best_sad_16x64;
    uint32_t *p_best_mv8x8          = context_ptr->p_best_mv8x8;
    uint32_t *p_best_mv16x16        = context_ptr->p_best_mv16x16;
    uint32_t *p_best_mv32x32        = context_ptr->p_best_mv32x32;
    uint32_t *p_best_mv64x64        = context_ptr->p_best_mv64x64;
    uint32_t *p_best_mv64x32        = context_ptr->p_best_mv64x32;
    uint32_t *p_best_mv32x16        = context_ptr->p_best_mv32x16;
    uint32_t *p_best_mv16x8         = context_ptr->p_best_mv16x8;
    uint32_t *p_best_mv32x64        = context_ptr->p_best_mv32x64;
    uint32_t *p_best_mv16x32        = context_ptr->p_best_mv16x32;
    uint32_t *p_best_mv8x16         = context_ptr->p_best_mv8x16;
    uint32_t *p_best_mv32x8         = context_ptr->p_best_mv32x8;
    uint32_t *p_best_mv8x32         = context_ptr->p_best_mv8x32;
    uint32_t *p_sad32x32            = context_ptr->p_sad32x32;
    uint32_t *p_sad16x16            = context_ptr->p_sad16x16;
    uint32_t *p_sad8x8              = context_ptr->p_sad8x8;
    uint32_t *p_best_mv64x16        = context_ptr->p_best_mv64x16;
    uint32_t *p_best_mv16x64        = context_ptr->p_best_mv16x64;

    // TODO: block_index search_position_index could be removed

    const uint32_t src_stride = context_ptr->sb_src_stride;
    src_next_16x16_offset     = src_stride << 4;

    //---- 16x16 : 0
    block_index           = 0;
    search_position_index = search_position_tl_index;

    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[0],
                                  &p_best_sad_16x16[0],
                                  &p_best_mv8x8[0],
                                  &p_best_mv16x16[0],
                                  curr_mv,
                                  &p_sad16x16[0],
                                  &p_sad8x8[0],
                                  sub_sad);

    //---- 16x16 : 1
    block_index           = block_index + 16;
    search_position_index = search_position_tl_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[4],
                                  &p_best_sad_16x16[1],
                                  &p_best_mv8x8[4],
                                  &p_best_mv16x16[1],
                                  curr_mv,
                                  &p_sad16x16[1],
                                  &p_sad8x8[4],
                                  sub_sad);
    //---- 16x16 : 4
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;

    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[16],
                                  &p_best_sad_16x16[4],
                                  &p_best_mv8x8[16],
                                  &p_best_mv16x16[4],
                                  curr_mv,
                                  &p_sad16x16[4],
                                  &p_sad8x8[16],
                                  sub_sad);

    //---- 16x16 : 5
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[20],
                                  &p_best_sad_16x16[5],
                                  &p_best_mv8x8[20],
                                  &p_best_mv16x16[5],
                                  curr_mv,
                                  &p_sad16x16[5],
                                  &p_sad8x8[20],
                                  sub_sad);

    //---- 16x16 : 2
    block_index           = src_next_16x16_offset;
    search_position_index = search_position_tl_index + ref_next_16x16_offset;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[8],
                                  &p_best_sad_16x16[2],
                                  &p_best_mv8x8[8],
                                  &p_best_mv16x16[2],
                                  curr_mv,
                                  &p_sad16x16[2],
                                  &p_sad8x8[8],
                                  sub_sad);
    //---- 16x16 : 3
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[12],
                                  &p_best_sad_16x16[3],
                                  &p_best_mv8x8[12],
                                  &p_best_mv16x16[3],
                                  curr_mv,
                                  &p_sad16x16[3],
                                  &p_sad8x8[12],
                                  sub_sad);
    //---- 16x16 : 6
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[24],
                                  &p_best_sad_16x16[6],
                                  &p_best_mv8x8[24],
                                  &p_best_mv16x16[6],
                                  curr_mv,
                                  &p_sad16x16[6],
                                  &p_sad8x8[24],
                                  sub_sad);
    //---- 16x16 : 7
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[28],
                                  &p_best_sad_16x16[7],
                                  &p_best_mv8x8[28],
                                  &p_best_mv16x16[7],
                                  curr_mv,
                                  &p_sad16x16[7],
                                  &p_sad8x8[28],
                                  sub_sad);

    //---- 16x16 : 8
    block_index           = (src_next_16x16_offset << 1);
    search_position_index = search_position_tl_index + (ref_next_16x16_offset << 1);
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[32],
                                  &p_best_sad_16x16[8],
                                  &p_best_mv8x8[32],
                                  &p_best_mv16x16[8],
                                  curr_mv,
                                  &p_sad16x16[8],
                                  &p_sad8x8[32],
                                  sub_sad);
    //---- 16x16 : 9
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[36],
                                  &p_best_sad_16x16[9],
                                  &p_best_mv8x8[36],
                                  &p_best_mv16x16[9],
                                  curr_mv,
                                  &p_sad16x16[9],
                                  &p_sad8x8[36],
                                  sub_sad);
    //---- 16x16 : 12
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[48],
                                  &p_best_sad_16x16[12],
                                  &p_best_mv8x8[48],
                                  &p_best_mv16x16[12],
                                  curr_mv,
                                  &p_sad16x16[12],
                                  &p_sad8x8[48],
                                  sub_sad);
    //---- 16x16 : 13
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[52],
                                  &p_best_sad_16x16[13],
                                  &p_best_mv8x8[52],
                                  &p_best_mv16x16[13],
                                  curr_mv,
                                  &p_sad16x16[13],
                                  &p_sad8x8[52],
                                  sub_sad);

    //---- 16x16 : 10
    block_index           = (src_next_16x16_offset * 3);
    search_position_index = search_position_tl_index + (ref_next_16x16_offset * 3);
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[40],
                                  &p_best_sad_16x16[10],
                                  &p_best_mv8x8[40],
                                  &p_best_mv16x16[10],
                                  curr_mv,
                                  &p_sad16x16[10],
                                  &p_sad8x8[40],
                                  sub_sad);
    //---- 16x16 : 11
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[44],
                                  &p_best_sad_16x16[11],
                                  &p_best_mv8x8[44],
                                  &p_best_mv16x16[11],
                                  curr_mv,
                                  &p_sad16x16[11],
                                  &p_sad8x8[44],
                                  sub_sad);
    //---- 16x16 : 14
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[56],
                                  &p_best_sad_16x16[14],
                                  &p_best_mv8x8[56],
                                  &p_best_mv16x16[14],
                                  curr_mv,
                                  &p_sad16x16[14],
                                  &p_sad8x8[56],
                                  sub_sad);
    //---- 16x16 : 15
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    ext_sad_calculation_8x8_16x16(src_ptr + block_index,
                                  src_stride,
                                  ref_ptr + search_position_index,
                                  ref_luma_stride,
                                  &p_best_sad_8x8[60],
                                  &p_best_sad_16x16[15],
                                  &p_best_mv8x8[60],
                                  &p_best_mv16x16[15],
                                  curr_mv,
                                  &p_sad16x16[15],
                                  &p_sad8x8[60],
                                  sub_sad);

    ext_sad_calculation_32x32_64x64(p_sad16x16,
                                    p_best_sad_32x32,
                                    p_best_sad_64x64,
                                    p_best_mv32x32,
                                    p_best_mv64x64,
                                    curr_mv,
                                    &p_sad32x32[0]);

    ext_sad_calculation(p_sad8x8,
                        p_sad16x16,
                        p_sad32x32,
                        p_best_sad_64x32,
                        p_best_mv64x32,
                        p_best_sad_32x16,
                        p_best_mv32x16,
                        p_best_sad_16x8,
                        p_best_mv16x8,
                        p_best_sad_32x64,
                        p_best_mv32x64,
                        p_best_sad_16x32,
                        p_best_mv16x32,
                        p_best_sad_8x16,
                        p_best_mv8x16,
                        p_best_sad_32x8,
                        p_best_mv32x8,
                        p_best_sad_8x32,
                        p_best_mv8x32,
                        p_best_sad_64x16,
                        p_best_mv64x16,
                        p_best_sad_16x64,
                        p_best_mv16x64,
                        curr_mv);
}

/*******************************************
 * get_search_point_results
 *******************************************/
static void get_search_point_results(
    MeContext *context_ptr, // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t   list_index, // input parameter, reference list index
    uint32_t   ref_pic_index,
    uint32_t   search_region_index, // input parameter, search area origin, used to
    // point to reference samples
    int32_t x_search_index, // input parameter, search region position in the
    // horizontal direction, used to derive xMV
    int32_t y_search_index) // input parameter, search region position in the
// vertical direction, used to derive yMV
{
    const EbBool sub_sad = (context_ptr->me_search_method == SUB_SAD_SEARCH);
    uint8_t *    src_ptr = context_ptr->sb_src_ptr;

    // uint8_t  *ref_ptr = ref_pic_ptr->buffer_y; // NADER
    uint8_t *ref_ptr =
        context_ptr->integer_buffer_ptr[list_index][ref_pic_index] + (ME_FILTER_TAP >> 1) +
        ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[list_index][ref_pic_index]);
    // uint32_t ref_luma_stride = ref_pic_ptr->stride_y; // NADER
    uint32_t ref_luma_stride = context_ptr->interpolated_full_stride[list_index][ref_pic_index];

    uint32_t search_position_tl_index = search_region_index;
    uint32_t search_position_index;
    uint32_t block_index;

    uint32_t src_next_16x16_offset = (BLOCK_SIZE_64 << 4);
    // uint32_t ref_next_16x16_offset = (ref_pic_ptr->stride_y << 4); // NADER
    uint32_t ref_next_16x16_offset = (ref_luma_stride << 4);

    uint32_t curr_mv_1 = (((uint16_t)y_search_index) << 18);
    uint16_t curr_mv_2 = (((uint16_t)x_search_index << 2));
    uint32_t curr_mv   = curr_mv_1 | curr_mv_2;

    uint32_t *p_best_sad_8x8   = context_ptr->p_best_sad_8x8;
    uint32_t *p_best_sad_16x16 = context_ptr->p_best_sad_16x16;
    uint32_t *p_best_sad_32x32 = context_ptr->p_best_sad_32x32;
    uint32_t *p_best_sad_64x64 = context_ptr->p_best_sad_64x64;

    uint32_t *p_best_mv8x8   = context_ptr->p_best_mv8x8;
    uint32_t *p_best_mv16x16 = context_ptr->p_best_mv16x16;
    uint32_t *p_best_mv32x32 = context_ptr->p_best_mv32x32;
    uint32_t *p_best_mv64x64 = context_ptr->p_best_mv64x64;
    uint32_t *p_sad16x16     = context_ptr->p_sad16x16;

    // TODO: block_index search_position_index could be removed

    const uint32_t src_stride = context_ptr->sb_src_stride;
    src_next_16x16_offset     = src_stride << 4;

    //---- 16x16 : 0
    block_index           = 0;
    search_position_index = search_position_tl_index;

    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[0],
                              &p_best_sad_16x16[0],
                              &p_best_mv8x8[0],
                              &p_best_mv16x16[0],
                              curr_mv,
                              &p_sad16x16[0],
                              sub_sad);

    //---- 16x16 : 1
    block_index           = block_index + 16;
    search_position_index = search_position_tl_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[4],
                              &p_best_sad_16x16[1],
                              &p_best_mv8x8[4],
                              &p_best_mv16x16[1],
                              curr_mv,
                              &p_sad16x16[1],
                              sub_sad);
    //---- 16x16 : 4
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;

    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[16],
                              &p_best_sad_16x16[4],
                              &p_best_mv8x8[16],
                              &p_best_mv16x16[4],
                              curr_mv,
                              &p_sad16x16[4],
                              sub_sad);

    //---- 16x16 : 5
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[20],
                              &p_best_sad_16x16[5],
                              &p_best_mv8x8[20],
                              &p_best_mv16x16[5],
                              curr_mv,
                              &p_sad16x16[5],
                              sub_sad);

    //---- 16x16 : 2
    block_index           = src_next_16x16_offset;
    search_position_index = search_position_tl_index + ref_next_16x16_offset;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[8],
                              &p_best_sad_16x16[2],
                              &p_best_mv8x8[8],
                              &p_best_mv16x16[2],
                              curr_mv,
                              &p_sad16x16[2],
                              sub_sad);
    //---- 16x16 : 3
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[12],
                              &p_best_sad_16x16[3],
                              &p_best_mv8x8[12],
                              &p_best_mv16x16[3],
                              curr_mv,
                              &p_sad16x16[3],
                              sub_sad);
    //---- 16x16 : 6
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[24],
                              &p_best_sad_16x16[6],
                              &p_best_mv8x8[24],
                              &p_best_mv16x16[6],
                              curr_mv,
                              &p_sad16x16[6],
                              sub_sad);
    //---- 16x16 : 7
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[28],
                              &p_best_sad_16x16[7],
                              &p_best_mv8x8[28],
                              &p_best_mv16x16[7],
                              curr_mv,
                              &p_sad16x16[7],
                              sub_sad);

    //---- 16x16 : 8
    block_index           = (src_next_16x16_offset << 1);
    search_position_index = search_position_tl_index + (ref_next_16x16_offset << 1);
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[32],
                              &p_best_sad_16x16[8],
                              &p_best_mv8x8[32],
                              &p_best_mv16x16[8],
                              curr_mv,
                              &p_sad16x16[8],
                              sub_sad);
    //---- 16x16 : 9
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[36],
                              &p_best_sad_16x16[9],
                              &p_best_mv8x8[36],
                              &p_best_mv16x16[9],
                              curr_mv,
                              &p_sad16x16[9],
                              sub_sad);
    //---- 16x16 : 12
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[48],
                              &p_best_sad_16x16[12],
                              &p_best_mv8x8[48],
                              &p_best_mv16x16[12],
                              curr_mv,
                              &p_sad16x16[12],
                              sub_sad);
    //---- 16x16 : 13
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[52],
                              &p_best_sad_16x16[13],
                              &p_best_mv8x8[52],
                              &p_best_mv16x16[13],
                              curr_mv,
                              &p_sad16x16[13],
                              sub_sad);

    //---- 16x16 : 10
    block_index           = (src_next_16x16_offset * 3);
    search_position_index = search_position_tl_index + (ref_next_16x16_offset * 3);
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[40],
                              &p_best_sad_16x16[10],
                              &p_best_mv8x8[40],
                              &p_best_mv16x16[10],
                              curr_mv,
                              &p_sad16x16[10],
                              sub_sad);
    //---- 16x16 : 11
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[44],
                              &p_best_sad_16x16[11],
                              &p_best_mv8x8[44],
                              &p_best_mv16x16[11],
                              curr_mv,
                              &p_sad16x16[11],
                              sub_sad);
    //---- 16x16 : 14
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[56],
                              &p_best_sad_16x16[14],
                              &p_best_mv8x8[56],
                              &p_best_mv16x16[14],
                              curr_mv,
                              &p_sad16x16[14],
                              sub_sad);
    //---- 16x16 : 15
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    sad_calculation_8x8_16x16(src_ptr + block_index,
                              src_stride,
                              ref_ptr + search_position_index,
                              ref_luma_stride,
                              &p_best_sad_8x8[60],
                              &p_best_sad_16x16[15],
                              &p_best_mv8x8[60],
                              &p_best_mv16x16[15],
                              curr_mv,
                              &p_sad16x16[15],
                              sub_sad);

    sad_calculation_32x32_64x64(
        p_sad16x16, p_best_sad_32x32, p_best_sad_64x64, p_best_mv32x32, p_best_mv64x64, curr_mv);
}

/*******************************************
 * GetEightHorizontalSearchPointResultsAll85CUs
 *******************************************/
static void get_eight_horizontal_search_point_results_all_85_pus(
    MeContext *context_ptr, uint32_t list_index, uint32_t ref_pic_index,
    uint32_t search_region_index,
    int32_t  x_search_index, // input parameter, search region position in the
    // horizontal direction, used to derive xMV
    int32_t y_search_index) { // input parameter, search region position in the
    // vertical direction, used to derive yMV
    const EbBool sub_sad = (context_ptr->me_search_method == SUB_SAD_SEARCH);
    uint8_t *    src_ptr = context_ptr->sb_src_ptr;
    uint8_t *    ref_ptr =
        context_ptr->integer_buffer_ptr[list_index][ref_pic_index] + (ME_FILTER_TAP >> 1) +
        ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[list_index][ref_pic_index]);
    uint32_t ref_luma_stride = context_ptr->interpolated_full_stride[list_index][ref_pic_index];

    uint32_t search_position_tl_index = search_region_index;
    uint32_t search_position_index;
    uint32_t block_index;

    uint32_t src_next_16x16_offset = (BLOCK_SIZE_64 << 4);
    uint32_t ref_next_16x16_offset = (ref_luma_stride << 4);

    uint32_t curr_mv_y = (((uint16_t)y_search_index) << 18);
    uint16_t curr_mv_x = (((uint16_t)x_search_index << 2));
    uint32_t curr_mv   = curr_mv_y | curr_mv_x;

    uint32_t *p_best_sad_8x8   = context_ptr->p_best_sad_8x8;
    uint32_t *p_best_sad_16x16 = context_ptr->p_best_sad_16x16;
    uint32_t *p_best_sad_32x32 = context_ptr->p_best_sad_32x32;
    uint32_t *p_best_sad_64x64 = context_ptr->p_best_sad_64x64;

    uint32_t *p_best_mv8x8   = context_ptr->p_best_mv8x8;
    uint32_t *p_best_mv16x16 = context_ptr->p_best_mv16x16;
    uint32_t *p_best_mv32x32 = context_ptr->p_best_mv32x32;
    uint32_t *p_best_mv64x64 = context_ptr->p_best_mv64x64;

    uint16_t *p_sad16x16 = context_ptr->p_eight_pos_sad16x16;

    /*
    ----------------------    ----------------------
    |  16x16_0  |  16x16_1  |  16x16_4  |  16x16_5  |
    ----------------------    ----------------------
    |  16x16_2  |  16x16_3  |  16x16_6  |  16x16_7  |
    -----------------------   -----------------------
    |  16x16_8  |  16x16_9  |  16x16_12 |  16x16_13 |
    ----------------------    ----------------------
    |  16x16_10 |  16x16_11 |  16x16_14 |  16x16_15 |
    -----------------------   -----------------------
    */

    const uint32_t src_stride = context_ptr->sb_src_stride;
    src_next_16x16_offset     = src_stride << 4;

    //---- 16x16_0
    block_index           = 0;
    search_position_index = search_position_tl_index;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[0],
                                                           &p_best_mv8x8[0],
                                                           &p_best_sad_16x16[0],
                                                           &p_best_mv16x16[0],
                                                           curr_mv,
                                                           &p_sad16x16[0 * 8],
                                                           sub_sad);
    //---- 16x16_1
    block_index           = block_index + 16;
    search_position_index = search_position_tl_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[4],
                                                           &p_best_mv8x8[4],
                                                           &p_best_sad_16x16[1],
                                                           &p_best_mv16x16[1],
                                                           curr_mv,
                                                           &p_sad16x16[1 * 8],
                                                           sub_sad);
    //---- 16x16_4
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[16],
                                                           &p_best_mv8x8[16],
                                                           &p_best_sad_16x16[4],
                                                           &p_best_mv16x16[4],
                                                           curr_mv,
                                                           &p_sad16x16[4 * 8],
                                                           sub_sad);
    //---- 16x16_5
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[20],
                                                           &p_best_mv8x8[20],
                                                           &p_best_sad_16x16[5],
                                                           &p_best_mv16x16[5],
                                                           curr_mv,
                                                           &p_sad16x16[5 * 8],
                                                           sub_sad);

    //---- 16x16_2
    block_index           = src_next_16x16_offset;
    search_position_index = search_position_tl_index + ref_next_16x16_offset;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[8],
                                                           &p_best_mv8x8[8],
                                                           &p_best_sad_16x16[2],
                                                           &p_best_mv16x16[2],
                                                           curr_mv,
                                                           &p_sad16x16[2 * 8],
                                                           sub_sad);
    //---- 16x16_3
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[12],
                                                           &p_best_mv8x8[12],
                                                           &p_best_sad_16x16[3],
                                                           &p_best_mv16x16[3],
                                                           curr_mv,
                                                           &p_sad16x16[3 * 8],
                                                           sub_sad);
    //---- 16x16_6
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[24],
                                                           &p_best_mv8x8[24],
                                                           &p_best_sad_16x16[6],
                                                           &p_best_mv16x16[6],
                                                           curr_mv,
                                                           &p_sad16x16[6 * 8],
                                                           sub_sad);
    //---- 16x16_7
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[28],
                                                           &p_best_mv8x8[28],
                                                           &p_best_sad_16x16[7],
                                                           &p_best_mv16x16[7],
                                                           curr_mv,
                                                           &p_sad16x16[7 * 8],
                                                           sub_sad);

    //---- 16x16_8
    block_index           = (src_next_16x16_offset << 1);
    search_position_index = search_position_tl_index + (ref_next_16x16_offset << 1);
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[32],
                                                           &p_best_mv8x8[32],
                                                           &p_best_sad_16x16[8],
                                                           &p_best_mv16x16[8],
                                                           curr_mv,
                                                           &p_sad16x16[8 * 8],
                                                           sub_sad);
    //---- 16x16_9
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[36],
                                                           &p_best_mv8x8[36],
                                                           &p_best_sad_16x16[9],
                                                           &p_best_mv16x16[9],
                                                           curr_mv,
                                                           &p_sad16x16[9 * 8],
                                                           sub_sad);
    //---- 16x16_12
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[48],
                                                           &p_best_mv8x8[48],
                                                           &p_best_sad_16x16[12],
                                                           &p_best_mv16x16[12],
                                                           curr_mv,
                                                           &p_sad16x16[12 * 8],
                                                           sub_sad);
    //---- 16x1_13
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[52],
                                                           &p_best_mv8x8[52],
                                                           &p_best_sad_16x16[13],
                                                           &p_best_mv16x16[13],
                                                           curr_mv,
                                                           &p_sad16x16[13 * 8],
                                                           sub_sad);

    //---- 16x16_10
    block_index           = (src_next_16x16_offset * 3);
    search_position_index = search_position_tl_index + (ref_next_16x16_offset * 3);
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[40],
                                                           &p_best_mv8x8[40],
                                                           &p_best_sad_16x16[10],
                                                           &p_best_mv16x16[10],
                                                           curr_mv,
                                                           &p_sad16x16[10 * 8],
                                                           sub_sad);
    //---- 16x16_11
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[44],
                                                           &p_best_mv8x8[44],
                                                           &p_best_sad_16x16[11],
                                                           &p_best_mv16x16[11],
                                                           curr_mv,
                                                           &p_sad16x16[11 * 8],
                                                           sub_sad);
    //---- 16x16_14
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[56],
                                                           &p_best_mv8x8[56],
                                                           &p_best_sad_16x16[14],
                                                           &p_best_mv16x16[14],
                                                           curr_mv,
                                                           &p_sad16x16[14 * 8],
                                                           sub_sad);
    //---- 16x16_15
    block_index           = block_index + 16;
    search_position_index = search_position_index + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(src_ptr + block_index,
                                                           context_ptr->sb_src_stride,
                                                           ref_ptr + search_position_index,
                                                           ref_luma_stride,
                                                           &p_best_sad_8x8[60],
                                                           &p_best_mv8x8[60],
                                                           &p_best_sad_16x16[15],
                                                           &p_best_mv16x16[15],
                                                           curr_mv,
                                                           &p_sad16x16[15 * 8],
                                                           sub_sad);
    // 32x32 and 64x64
    get_eight_horizontal_search_point_results_32x32_64x64_pu(
        p_sad16x16, p_best_sad_32x32, p_best_sad_64x64, p_best_mv32x32, p_best_mv64x64, curr_mv);
}

/*******************************************
 * full_pel_search_sb
 *******************************************/
static void full_pel_search_sb(MeContext *context_ptr, uint32_t list_index, uint32_t ref_pic_index,
                               int16_t x_search_area_origin, int16_t y_search_area_origin,
                               uint32_t search_area_width, uint32_t search_area_height) {
    uint32_t x_search_index, y_search_index;

    uint32_t search_area_width_rest_8 = search_area_width & 7;
    uint32_t search_area_width_mult_8 = search_area_width - search_area_width_rest_8;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index < search_area_width_mult_8; x_search_index += 8) {
            // this function will do:  x_search_index, +1, +2, ..., +7
            get_eight_horizontal_search_point_results_all_85_pus(
                context_ptr,
                list_index,
                ref_pic_index,
                x_search_index +
                    y_search_index *
                        context_ptr->interpolated_full_stride[list_index][ref_pic_index],
                (int32_t)x_search_index + x_search_area_origin,
                (int32_t)y_search_index + y_search_area_origin);
        }

        for (x_search_index = search_area_width_mult_8; x_search_index < search_area_width;
             x_search_index++) {
            get_search_point_results(
                context_ptr,
                list_index,
                ref_pic_index,
                x_search_index +
                    y_search_index *
                        context_ptr->interpolated_full_stride[list_index][ref_pic_index],
                (int32_t)x_search_index + x_search_area_origin,
                (int32_t)y_search_index + y_search_area_origin);
        }
    }
}
/*******************************************
 * pu_half_pel_refinement
 *   performs Half Pel refinement for one PU
 *******************************************/
static void half_pel_refinement_block(
    MeContext *context_ptr, // input parameter, ME context Ptr, used to get SB Ptr
    uint8_t *ref_buffer, uint32_t ref_stride, uint32_t *p_best_ssd,
    uint32_t src_block_index, // input parameter, PU origin, used to point to
    // source samples
    uint8_t *pos_b_buffer, // input parameter, position "b" interpolated search
    // area Ptr
    uint8_t *pos_h_buffer, // input parameter, position "h" interpolated search
    // area Ptr
    uint8_t *pos_j_buffer, // input parameter, position "j" interpolated search
    // area Ptr
    uint32_t pu_width, // input parameter, PU width
    uint32_t pu_height, // input parameter, PU height
    int16_t  x_search_area_origin, // input parameter, search area origin in the
    // horizontal direction, used to point to
    // reference samples
    int16_t y_search_area_origin, // input parameter, search area origin in the
    // vertical direction, used to point to
    // reference samples
    uint32_t  search_area_height, // input parameter, search area height
    uint32_t  search_area_width, // input parameter, search area width
    uint32_t *p_best_sad, uint32_t *p_best_mv, uint8_t *p_sub_pel_direction,
    uint32_t *best_pervious_stage_mv, uint32_t ineteger_mv) {
    int32_t  search_region_index;
    uint64_t distortion_left_position     = 0;
    uint64_t distortion_top_position      = 0;
    uint64_t distortion_topleft_position  = 0;
    uint64_t distortion_topright_position = 0;
    int16_t  half_mv_x[8];
    int16_t  half_mv_y[8];
    int16_t  x_best_mv;
    int16_t  y_best_mv;
    int16_t  x_mv;
    int16_t  y_mv;
    int16_t  search_index_x;
    int16_t  search_index_y;
    (void)p_sub_pel_direction;
    (void)ineteger_mv;
    // copute distance between best mv and the integer mv candidate
    int16_t offset_x, offset_y;
    int8_t h_pel_search_wind = context_ptr->h_pel_search_wind;

    for (offset_x = -h_pel_search_wind; offset_x <= h_pel_search_wind; offset_x++) {
        for (offset_y = -h_pel_search_wind; offset_y <= h_pel_search_wind; offset_y++) {
            x_best_mv            = _MVXT(*best_pervious_stage_mv);
            y_best_mv            = _MVYT(*best_pervious_stage_mv);
            x_mv                 = x_best_mv + (offset_x * 4);
            y_mv                 = y_best_mv + (offset_y * 4);
            search_index_x       = (x_mv >> 2) - x_search_area_origin;
            search_index_y       = (y_mv >> 2) - y_search_area_origin;
            uint32_t integer_mv1 = (((uint16_t)(y_mv >> 2)) << 18);
            uint16_t integer_mv2 = (((uint16_t)(x_mv >> 2) << 2));
            uint32_t integer_mv  = integer_mv1 | integer_mv2;
            if (search_index_x < 0 || search_index_x > (int16_t)(search_area_width - 1)) {
                continue;
            }
            if (search_index_y < 0 || search_index_y > (int16_t)(search_area_height - 1)) {
                continue;
            }
            half_mv_x[0] = x_mv - 2; // L  position
            half_mv_x[1] = x_mv + 2; // R  position
            half_mv_x[2] = x_mv; // T  position
            half_mv_x[3] = x_mv; // b  position
            half_mv_x[4] = x_mv - 2; // TL position
            half_mv_x[5] = x_mv + 2; // TR position
            half_mv_x[6] = x_mv + 2; // BR position
            half_mv_x[7] = x_mv - 2; // BL position
            half_mv_y[0] = y_mv; // L  position
            half_mv_y[1] = y_mv; // R  position
            half_mv_y[2] = y_mv - 2; // T  position
            half_mv_y[3] = y_mv + 2; // b  position
            half_mv_y[4] = y_mv - 2; // TL position
            half_mv_y[5] = y_mv - 2; // TR position
            half_mv_y[6] = y_mv + 2; // BR position
            half_mv_y[7] = y_mv + 2; // BL position
            // Compute SSD for the best full search candidate
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                uint32_t integer_sse = (uint32_t)spatial_full_distortion_kernel(
                    context_ptr->sb_src_ptr,
                    src_block_index,
                    context_ptr->sb_src_stride,
                    ref_buffer,
                    search_index_y * ref_stride + search_index_x,
                    ref_stride,
                    pu_width,
                    pu_height);
                if (integer_sse < *p_best_ssd) {
                    *p_best_ssd = integer_sse;
                    *p_best_mv  = integer_mv;
                }
            }
            // L position
            search_region_index =
                search_index_x + (int16_t)context_ptr->interpolated_stride * search_index_y;
            if (context_ptr->fractional_search_method == SSD_SEARCH)
                distortion_left_position =
                    spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                   src_block_index,
                                                   context_ptr->sb_src_stride,
                                                   pos_b_buffer,
                                                   search_region_index,
                                                   context_ptr->interpolated_stride,
                                                   pu_width,
                                                   pu_height);
            else if (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                distortion_left_position =
                    (nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                    context_ptr->sb_src_stride << 1,
                                    &(pos_b_buffer[search_region_index]),
                                    context_ptr->interpolated_stride << 1,
                                    pu_height >> 1,
                                    pu_width))
                    << 1;
            else
                distortion_left_position =
                    nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                   context_ptr->sb_src_stride,
                                   &(pos_b_buffer[search_region_index]),
                                   context_ptr->interpolated_stride,
                                   pu_height,
                                   pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (distortion_left_position < *p_best_ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                                 context_ptr->sb_src_stride,
                                                 &(pos_b_buffer[search_region_index]),
                                                 context_ptr->interpolated_stride,
                                                 pu_height,
                                                 pu_width);
                    *p_best_mv  = ((uint16_t)half_mv_y[0] << 16) | ((uint16_t)half_mv_x[0]);
                    *p_best_ssd = (uint32_t)distortion_left_position;
                }
            } else {
                if (distortion_left_position < *p_best_sad) {
                    *p_best_sad = (uint32_t)distortion_left_position;
                    *p_best_mv  = ((uint16_t)half_mv_y[0] << 16) | ((uint16_t)half_mv_x[0]);
                }
            }
            // T position
            search_region_index =
                search_index_x + (int16_t)context_ptr->interpolated_stride * search_index_y;
            if (context_ptr->fractional_search_method == SSD_SEARCH)
                distortion_top_position =
                    spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                   src_block_index,
                                                   context_ptr->sb_src_stride,
                                                   pos_h_buffer,
                                                   search_region_index,
                                                   context_ptr->interpolated_stride,
                                                   pu_width,
                                                   pu_height);
            else if (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                distortion_top_position =
                    (nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                    context_ptr->sb_src_stride << 1,
                                    &(pos_h_buffer[search_region_index]),
                                    context_ptr->interpolated_stride << 1,
                                    pu_height >> 1,
                                    pu_width))
                    << 1;
            else
                distortion_top_position =
                    nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                   context_ptr->sb_src_stride,
                                   &(pos_h_buffer[search_region_index]),
                                   context_ptr->interpolated_stride,
                                   pu_height,
                                   pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (distortion_top_position < *p_best_ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                                 context_ptr->sb_src_stride,
                                                 &(pos_h_buffer[search_region_index]),
                                                 context_ptr->interpolated_stride,
                                                 pu_height,
                                                 pu_width);
                    *p_best_mv  = ((uint16_t)half_mv_y[2] << 16) | ((uint16_t)half_mv_x[2]);
                    *p_best_ssd = (uint32_t)distortion_top_position;
                }
            } else {
                if (distortion_top_position < *p_best_sad) {
                    *p_best_sad = (uint32_t)distortion_top_position;
                    *p_best_mv  = ((uint16_t)half_mv_y[2] << 16) | ((uint16_t)half_mv_x[2]);
                }
            }
            // TL position
            search_region_index =
                search_index_x + (int16_t)context_ptr->interpolated_stride * search_index_y;
            if (context_ptr->fractional_search_method == SSD_SEARCH)
                distortion_topleft_position =
                    spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                   src_block_index,
                                                   context_ptr->sb_src_stride,
                                                   pos_j_buffer,
                                                   search_region_index,
                                                   context_ptr->interpolated_stride,
                                                   pu_width,
                                                   pu_height);
            else if (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                distortion_topleft_position =
                    (nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                    context_ptr->sb_src_stride << 1,
                                    &(pos_j_buffer[search_region_index]),
                                    context_ptr->interpolated_stride << 1,
                                    pu_height >> 1,
                                    pu_width))
                    << 1;
            else
                distortion_topleft_position =
                    nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                   context_ptr->sb_src_stride,
                                   &(pos_j_buffer[search_region_index]),
                                   context_ptr->interpolated_stride,
                                   pu_height,
                                   pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (distortion_topleft_position < *p_best_ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                                 context_ptr->sb_src_stride,
                                                 &(pos_j_buffer[search_region_index]),
                                                 context_ptr->interpolated_stride,
                                                 pu_height,
                                                 pu_width);
                    *p_best_mv  = ((uint16_t)half_mv_y[4] << 16) | ((uint16_t)half_mv_x[4]);
                    *p_best_ssd = (uint32_t)distortion_topleft_position;
                }
            } else {
                if (distortion_topleft_position < *p_best_sad) {
                    *p_best_sad = (uint32_t)distortion_topleft_position;
                    *p_best_mv  = ((uint16_t)half_mv_y[4] << 16) | ((uint16_t)half_mv_x[4]);
                }
            }
            // TR position
            search_region_index++;
            if (context_ptr->fractional_search_method == SSD_SEARCH)
                distortion_topright_position =
                    spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                   src_block_index,
                                                   context_ptr->sb_src_stride,
                                                   pos_j_buffer,
                                                   search_region_index,
                                                   context_ptr->interpolated_stride,
                                                   pu_width,
                                                   pu_height);
            else if (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                distortion_topright_position =
                    (nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                    context_ptr->sb_src_stride << 1,
                                    &(pos_j_buffer[search_region_index]),
                                    context_ptr->interpolated_stride << 1,
                                    pu_height >> 1,
                                    pu_width))
                    << 1;
            else
                distortion_topright_position =
                    nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                   context_ptr->sb_src_stride,
                                   &(pos_j_buffer[search_region_index]),
                                   context_ptr->interpolated_stride,
                                   pu_height,
                                   pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (distortion_topright_position < *p_best_ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[src_block_index]),
                                                 context_ptr->sb_src_stride,
                                                 &(pos_j_buffer[search_region_index]),
                                                 context_ptr->interpolated_stride,
                                                 pu_height,
                                                 pu_width);
                    *p_best_mv  = ((uint16_t)half_mv_y[5] << 16) | ((uint16_t)half_mv_x[5]);
                    *p_best_ssd = (uint32_t)distortion_topright_position;
                }
            } else {
                if (distortion_topright_position < *p_best_sad) {
                    *p_best_sad = (uint32_t)distortion_topright_position;
                    *p_best_mv  = ((uint16_t)half_mv_y[5] << 16) | ((uint16_t)half_mv_x[5]);
                }
            }
        }
    }
    return;
}

/*******************************************
 * half_pel_search_sb
 *   performs Half Pel refinement for the 85 PUs
 *******************************************/
void half_pel_refinement_sb(
    PictureParentControlSet *pcs_ptr,
    MeContext *              context_ptr, // input/output parameter, ME context Ptr, used to
    // get/update ME results
    uint8_t *refBuffer, uint32_t ref_stride,
    uint8_t *pos_b_buffer, // input parameter, position "b" interpolated search
    // area Ptr
    uint8_t *pos_h_buffer, // input parameter, position "h" interpolated search
    // area Ptr
    uint8_t *pos_j_buffer, // input parameter, position "j" interpolated search
    // area Ptr
    int16_t x_search_area_origin, // input parameter, search area origin in the
    // horizontal direction, used to point to
    // reference samples
    int16_t y_search_area_origin, // input parameter, search area origin in the
    // vertical direction, used to point to
    // reference samples
    uint32_t search_area_height, // input parameter, search area height
    uint32_t search_area_width, // input parameter, search area width
    uint8_t list_index, // reference picture list
    uint8_t ref_pic_index, // reference picture index
    uint32_t inetger_mv) {
    uint32_t idx;
    uint32_t pu_index;
    uint32_t block_index_shift_x;
    uint32_t block_index_shift_y;
    uint32_t src_block_index;
    uint32_t posb_buffer_index;
    uint32_t posh_buffer_index;
    uint32_t posj_buffer_index;
    // 64x64 [1 partition]
    half_pel_refinement_block(context_ptr,
                                &(refBuffer[0]),
                                ref_stride,
                                context_ptr->p_best_ssd64x64,
                                0,
                                &(pos_b_buffer[0]),
                                &(pos_h_buffer[0]),
                                &(pos_j_buffer[0]),
                                64,
                                64,
                                x_search_area_origin,
                                y_search_area_origin,
                                search_area_height,
                                search_area_width,
                                context_ptr->p_best_sad_64x64,
                                context_ptr->p_best_mv64x64,
                                &context_ptr->psub_pel_direction64x64,
                                context_ptr->p_best_full_pel_mv64x64,
                                inetger_mv);
    // 32x32 [4 partitions]
    for (pu_index = 0; pu_index < 4; ++pu_index) {
        block_index_shift_x = (pu_index & 0x01) << 5;
        block_index_shift_y = (pu_index >> 1) << 5;
        src_block_index = block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
        posb_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        posh_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        posj_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        half_pel_refinement_block(
            context_ptr,
            &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
            ref_stride,
            &context_ptr->p_best_ssd32x32[pu_index],
            src_block_index,
            &(pos_b_buffer[posb_buffer_index]),
            &(pos_h_buffer[posh_buffer_index]),
            &(pos_j_buffer[posj_buffer_index]),
            32,
            32,
            x_search_area_origin,
            y_search_area_origin,
            search_area_height,
            search_area_width,
            &context_ptr->p_best_sad_32x32[pu_index],
            &context_ptr->p_best_mv32x32[pu_index],
            &context_ptr->psub_pel_direction32x32[pu_index],
            &context_ptr->p_best_full_pel_mv32x32[pu_index],
            inetger_mv);
    }
    // 16x16 [16 partitions]
    for (pu_index = 0; pu_index < 16; ++pu_index) {
        idx                 = tab16x16[pu_index];
        block_index_shift_x = (pu_index & 0x03) << 4;
        block_index_shift_y = (pu_index >> 2) << 4;
        src_block_index = block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
        posb_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        posh_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        posj_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        half_pel_refinement_block(
            context_ptr,
            &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
            ref_stride,
            &context_ptr->p_best_ssd16x16[idx],
            src_block_index,
            &(pos_b_buffer[posb_buffer_index]),
            &(pos_h_buffer[posh_buffer_index]),
            &(pos_j_buffer[posj_buffer_index]),
            16,
            16,
            x_search_area_origin,
            y_search_area_origin,
            search_area_height,
            search_area_width,
            &context_ptr->p_best_sad_16x16[idx],
            &context_ptr->p_best_mv16x16[idx],
            &context_ptr->psub_pel_direction16x16[idx],
            &context_ptr->p_best_full_pel_mv16x16[idx],
            inetger_mv);
    }
    // 8x8   [64 partitions]
    for (pu_index = 0; pu_index < 64; ++pu_index) {
        idx                 = tab8x8[pu_index]; // TODO bitwise this
        block_index_shift_x = (pu_index & 0x07) << 3;
        block_index_shift_y = (pu_index >> 3) << 3;
        src_block_index = block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
        posb_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        posh_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        posj_buffer_index =
            block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
        half_pel_refinement_block(
            context_ptr,
            &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
            ref_stride,
            &context_ptr->p_best_ssd8x8[idx],
            src_block_index,
            &(pos_b_buffer[posb_buffer_index]),
            &(pos_h_buffer[posh_buffer_index]),
            &(pos_j_buffer[posj_buffer_index]),
            8,
            8,
            x_search_area_origin,
            y_search_area_origin,
            search_area_height,
            search_area_width,
            &context_ptr->p_best_sad_8x8[idx],
            &context_ptr->p_best_mv8x8[idx],
            &context_ptr->psub_pel_direction8x8[idx],
            &context_ptr->p_best_full_pel_mv8x8[idx],
            inetger_mv);
    }
    uint8_t perform_nsq_flag = 1;
    perform_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 1 && (list_index != context_ptr->best_list_idx || ref_pic_index != context_ptr->best_ref_idx)) ? 0 : perform_nsq_flag;
    perform_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 2 && ref_pic_index) ? 0: perform_nsq_flag;
    perform_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 3 ) ? 0 : perform_nsq_flag;
    if(perform_nsq_flag){
        if (pcs_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
            // 64x32
            for (pu_index = 0; pu_index < 2; ++pu_index) {
                block_index_shift_x = 0;
                block_index_shift_y = pu_index << 5;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd64x32[pu_index],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    64,
                    32,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_64x32[pu_index],
                    &context_ptr->p_best_mv64x32[pu_index],
                    &context_ptr->psub_pel_direction64x32[pu_index],
                    &context_ptr->p_best_full_pel_mv64x32[pu_index],
                    inetger_mv);
            }
            // 32x16
            for (pu_index = 0; pu_index < 8; ++pu_index) {
                idx = tab32x16[pu_index]; // TODO bitwise this
                block_index_shift_x = (pu_index & 0x01) << 5;
                block_index_shift_y = (pu_index >> 1) << 4;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd32x16[idx],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    32,
                    16,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_32x16[idx],
                    &context_ptr->p_best_mv32x16[idx],
                    &context_ptr->psub_pel_direction32x16[idx],
                    &context_ptr->p_best_full_pel_mv32x16[idx],
                    inetger_mv);
            }
            // 16x8
            for (pu_index = 0; pu_index < 32; ++pu_index) {
                idx = tab16x8[pu_index];
                block_index_shift_x = (pu_index & 0x03) << 4;
                block_index_shift_y = (pu_index >> 2) << 3;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd16x8[idx],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    16,
                    8,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_16x8[idx],
                    &context_ptr->p_best_mv16x8[idx],
                    &context_ptr->psub_pel_direction16x8[idx],
                    &context_ptr->p_best_full_pel_mv16x8[idx],
                    inetger_mv);
            }
            // 32x64
            for (pu_index = 0; pu_index < 2; ++pu_index) {
                block_index_shift_x = pu_index << 5;
                block_index_shift_y = 0;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd32x64[pu_index],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    32,
                    64,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_32x64[pu_index],
                    &context_ptr->p_best_mv32x64[pu_index],
                    &context_ptr->psub_pel_direction32x64[pu_index],
                    &context_ptr->p_best_full_pel_mv32x64[pu_index],
                    inetger_mv);
            }
            // 16x32
            for (pu_index = 0; pu_index < 8; ++pu_index) {
                idx = tab16x32[pu_index];
                block_index_shift_x = (pu_index & 0x03) << 4;
                block_index_shift_y = (pu_index >> 2) << 5;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd16x32[idx],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    16,
                    32,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_16x32[idx],
                    &context_ptr->p_best_mv16x32[idx],
                    &context_ptr->psub_pel_direction16x32[idx],
                    &context_ptr->p_best_full_pel_mv16x32[idx],
                    inetger_mv);
            }
            // 8x16
            for (pu_index = 0; pu_index < 32; ++pu_index) {
                idx = tab8x16[pu_index];
                block_index_shift_x = (pu_index & 0x07) << 3;
                block_index_shift_y = (pu_index >> 3) << 4;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd8x16[idx],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    8,
                    16,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_8x16[idx],
                    &context_ptr->p_best_mv8x16[idx],
                    &context_ptr->psub_pel_direction8x16[idx],
                    &context_ptr->p_best_full_pel_mv8x16[idx],
                    inetger_mv);
            }
            // 32x8
            for (pu_index = 0; pu_index < 16; ++pu_index) {
                idx = tab32x8[pu_index];
                block_index_shift_x = (pu_index & 0x01) << 5;
                block_index_shift_y = (pu_index >> 1) << 3;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd32x8[idx],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    32,
                    8,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_32x8[idx],
                    &context_ptr->p_best_mv32x8[idx],
                    &context_ptr->psub_pel_direction32x8[idx],
                    &context_ptr->p_best_full_pel_mv32x8[idx],
                    inetger_mv);
            }
            for (pu_index = 0; pu_index < 16; ++pu_index) {
                idx = tab8x32[pu_index];
                block_index_shift_x = (pu_index & 0x07) << 3;
                block_index_shift_y = (pu_index >> 3) << 5;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd8x32[idx],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    8,
                    32,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_8x32[idx],
                    &context_ptr->p_best_mv8x32[idx],
                    &context_ptr->psub_pel_direction8x32[idx],
                    &context_ptr->p_best_full_pel_mv8x32[idx],
                    inetger_mv);
            }
            for (pu_index = 0; pu_index < 4; ++pu_index) {
                idx = pu_index;
                block_index_shift_x = 0;
                block_index_shift_y = pu_index << 4;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd64x16[idx],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    64,
                    16,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_64x16[idx],
                    &context_ptr->p_best_mv64x16[idx],
                    &context_ptr->psub_pel_direction64x16[idx],
                    &context_ptr->p_best_full_pel_mv64x16[idx],
                    inetger_mv);
            }
            for (pu_index = 0; pu_index < 4; ++pu_index) {
                idx = pu_index;
                block_index_shift_x = pu_index << 4;
                block_index_shift_y = 0;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->sb_src_stride;
                posb_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posh_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                posj_buffer_index =
                    block_index_shift_x + block_index_shift_y * context_ptr->interpolated_stride;
                half_pel_refinement_block(
                    context_ptr,
                    &(refBuffer[block_index_shift_y * ref_stride + block_index_shift_x]),
                    ref_stride,
                    &context_ptr->p_best_ssd16x64[idx],
                    src_block_index,
                    &(pos_b_buffer[posb_buffer_index]),
                    &(pos_h_buffer[posh_buffer_index]),
                    &(pos_j_buffer[posj_buffer_index]),
                    16,
                    64,
                    x_search_area_origin,
                    y_search_area_origin,
                    search_area_height,
                    search_area_width,
                    &context_ptr->p_best_sad_16x64[idx],
                    &context_ptr->p_best_mv16x64[idx],
                    &context_ptr->psub_pel_direction16x64[idx],
                    &context_ptr->p_best_full_pel_mv16x64[idx],
                    inetger_mv);
            }
        }
    }
    return;
}
/*******************************************
 * open_loop_me_half_pel_search_sblock
 *******************************************/
static void open_loop_me_half_pel_search_sblock(
    PictureParentControlSet *pcs_ptr, MeContext *context_ptr, uint32_t list_index,
    uint32_t ref_pic_index, int16_t x_search_area_origin, int16_t y_search_area_origin,
    uint32_t search_area_width, uint32_t search_area_height) {
    half_pel_refinement_sb(
        pcs_ptr,
        context_ptr,
        context_ptr->integer_buffer_ptr[list_index][ref_pic_index] + (ME_FILTER_TAP >> 1) +
            ((ME_FILTER_TAP >> 1) *
             context_ptr->interpolated_full_stride[list_index][ref_pic_index]),
        context_ptr->interpolated_full_stride[list_index][ref_pic_index],
        &(context_ptr->pos_b_buffer[list_index][ref_pic_index]
                                   [(ME_FILTER_TAP >> 1) * context_ptr->interpolated_stride]),

        &(context_ptr->pos_h_buffer[list_index][ref_pic_index][1]),
        &(context_ptr->pos_j_buffer[list_index][ref_pic_index][0]),
        x_search_area_origin,
        y_search_area_origin,
        search_area_height,
        search_area_width,
        list_index,
        ref_pic_index,
        0);

    if (pcs_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
        uint8_t gather_nsq_flag = 0;
        gather_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 1 && (list_index != context_ptr->best_list_idx || ref_pic_index != context_ptr->best_ref_idx)) ? 1 : gather_nsq_flag;
        gather_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 2 && ref_pic_index) ? 1 : gather_nsq_flag;
        gather_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 3) ? 1 : gather_nsq_flag;
        if (gather_nsq_flag)
            generate_nsq_mv(context_ptr);
    }
}

/*******************************************
 * open_loop_me_fullpel_search_sblock
 *******************************************/
static void open_loop_me_fullpel_search_sblock(MeContext *context_ptr, uint32_t list_index,
                                               uint32_t ref_pic_index, int16_t x_search_area_origin,
                                               int16_t  y_search_area_origin,
                                               uint32_t search_area_width,
                                               uint32_t search_area_height,
                                               uint8_t pic_depth_mode) {
    uint32_t x_search_index, y_search_index;
    uint32_t search_area_width_rest_8 = search_area_width & 7;
    uint32_t search_area_width_mult_8 = search_area_width - search_area_width_rest_8;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index < search_area_width_mult_8; x_search_index += 8) {
            // this function will do:  x_search_index, +1, +2, ..., +7
            open_loop_me_get_eight_search_point_results_block(
                context_ptr,
                list_index,
                ref_pic_index,
                x_search_index +
                    y_search_index *
                        context_ptr->interpolated_full_stride[list_index][ref_pic_index],
                (int32_t)x_search_index + x_search_area_origin,
                (int32_t)y_search_index + y_search_area_origin);
        }

        for (x_search_index = search_area_width_mult_8; x_search_index < search_area_width;
             x_search_index++) {
            open_loop_me_get_search_point_results_block(
                context_ptr,
                list_index,
                ref_pic_index,
                x_search_index +
                    y_search_index *
                        context_ptr->interpolated_full_stride[list_index][ref_pic_index],
                (int32_t)x_search_index + x_search_area_origin,
                (int32_t)y_search_index + y_search_area_origin);
        }
    }
    if (pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
        uint8_t gather_nsq_flag = 0;
        gather_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 1 && (list_index != context_ptr->best_list_idx || ref_pic_index != context_ptr->best_ref_idx)) ? 1 : gather_nsq_flag;
        gather_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 2 && ref_pic_index) ? 1 : gather_nsq_flag;
        gather_nsq_flag = (context_ptr->inherit_rec_mv_from_sq_block == 3) ? 1 : gather_nsq_flag;
        if (gather_nsq_flag)
            generate_nsq_mv(context_ptr);
    }
}

#ifndef AVCCODEL
/*******************************************
 * HorizontalPelInterpolation
 *   interpolates the search region in the horizontal direction
 *******************************************/
static void HorizontalPelInterpolation(
    uint8_t *      src, // input parameter, input samples Ptr
    uint32_t       src_stride, // input parameter, input stride
    uint32_t       width, // input parameter, input area width
    uint32_t       height, // input parameter, input area height
    const int32_t *ifCoeff, // input parameter, interpolation filter coefficients Ptr
    uint32_t       inputBitDepth, // input parameter, input sample bit depth
    uint32_t       dst_stride, // input parameter, output stride
    uint8_t *      dst) // output parameter, interpolated samples Ptr
{
    uint32_t      x, y;
    const int32_t maxSampleValue = (1 << inputBitDepth) - 1;
    const int32_t ifOffset       = 1 << (if_shift - 1);
    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            dst[x] = (uint8_t)CLIP3(
                0,
                (int32_t)maxSampleValue,
                ((((int32_t)src[x] + (int32_t)src[x + 3]) * ifCoeff[0] +
                  ((int32_t)src[x + 1] + (int32_t)src[x + 2]) * ifCoeff[1] + ifOffset) >>
                 if_shift));
        }
        src += src_stride;
        dst += dst_stride;
    }

    return;
}

/*******************************************
 * VerticalPelInterpolation
 *   interpolates the serach region in the vertical direction
 *******************************************/
static void VerticalPelInterpolation(
    uint8_t *     src, // input parameter, input samples ptr
    uint32_t      src_stride, // input parameter, input stride
    uint32_t      width, // input parameter, input area width
    uint32_t      height, // input parameter, input area height
    const int32_t ifCoeff[4], // input parameter, interpolation filter coefficients Ptr
    uint32_t      inputBitDepth, // input parameter, input sample bit depth
    uint32_t      dst_stride, // input parameter, output stride
    uint8_t *     dst) // output parameter, interpolated samples Ptr
{
    uint32_t x, y;

    const int32_t maxSampleValue = (1 << inputBitDepth) - 1;
    const int32_t ifOffset       = 1 << (if_shift - 1);

    const uint32_t srcStride2 = src_stride << 1;
    const uint32_t srcStride3 = srcStride2 + src_stride;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (uint8_t)CLIP3(
                0,
                maxSampleValue,
                ((((int32_t)src[x] + (int32_t)src[x + srcStride3]) * ifCoeff[0] +
                  ((int32_t)src[x + src_stride] + (int32_t)src[x + srcStride2]) * ifCoeff[1] +
                  ifOffset) >>
                 if_shift));
        }
        src += src_stride;
        dst += dst_stride;
    }

    return;
}

/*******************************************
 * AvcStyleInterpolation
 *   interpolates the search region in the horizontal direction
 *******************************************/
static void AvcStyleInterpolation(uint8_t *srcOne, // input parameter, input samples Ptr
                                  uint32_t srcOneStride, // input parameter, input stride
                                  uint8_t *srcTwo, // input parameter, input samples Ptr
                                  uint32_t srcTwoStride, // input parameter, input stride
                                  uint32_t width, // input parameter, input area width
                                  uint32_t height, // input parameter, input area height
                                  uint32_t inputBitDepth, // input parameter, input sample bit depth
                                  uint32_t dst_stride, // input parameter, output stride
                                  uint8_t *dst) // output parameter, interpolated samples Ptr
{
    uint32_t x, y;
    int32_t  maxSampleValue = POW2(inputBitDepth) - 1;

    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            dst[x] =
                (uint8_t)CLIP3(0,
                               (int32_t)maxSampleValue,
                               (((int32_t)srcOne[x] + (int32_t)srcTwo[x] + 1) >> IFShiftAvcStyle));
        }
        srcOne += srcOneStride;
        srcTwo += srcTwoStride;
        dst += dst_stride;
    }

    return;
}
#endif
/*******************************************
 * InterpolateSearchRegion AVC
 *   interpolates the search area
 *   the whole search area is interpolated 15 times
 *   for each sub position an interpolation is done
 *   15 buffers are required for the storage of the interpolated samples.
 *   F0: {-4, 54, 16, -2}
 *   F1: {-4, 36, 36, -4}
 *   F2: {-2, 16, 54, -4}
 ********************************************/
void interpolate_search_region_avc(
    MeContext *context_ptr, // input/output parameter, ME context ptr, used to
    // get/set interpolated search area Ptr
    uint32_t list_index, // Refrence picture list index
    uint32_t ref_pic_index,
    uint8_t *searchRegionBuffer, // input parameter, search region index, used
    // to point to reference samples
    uint32_t luma_stride, // input parameter, reference Picture stride
    uint32_t search_area_width, // input parameter, search area width
    uint32_t search_area_height, // input parameter, search area height
    uint32_t inputBitDepth) // input parameter, input sample bit depth
{
    //      0    1    2    3
    // 0    A    a    b    c
    // 1    d    e    f    g
    // 2    h    i    j    k
    // 3    n    p    q    r

    // Position  Frac-pos Y  Frac-pos X  Horizontal filter  Vertical filter
    // A         0           0           -                  -
    // a         0           1           F0                 -
    // b         0           2           F1                 -
    // c         0           3           F2                 -
    // d         1           0           -                  F0
    // e         1           1           F0                 F0
    // f         1           2           F1                 F0
    // g         1           3           F2                 F0
    // h         2           0           -                  F1
    // i         2           1           F0                 F1
    // j         2           2           F1                 F1
    // k         2           3           F2                 F1
    // n         3           0           -                  F2
    // p         3           1           F0                 F2
    // q         3           2           F1                 F2
    // r         3           3           F2                 F2

    // Start a b c

    // The Search area needs to be a multiple of 8 to align with the ASM kernel
    // Also the search area must be oversized by 2 to account for edge
    // conditions
    uint32_t search_area_width_for_asm = ROUND_UP_MUL_8(search_area_width + 2);

#ifdef AVCCODEL

    (void)inputBitDepth;
    // Half pel interpolation of the search region using f1 -> pos_b_buffer
    if (search_area_width_for_asm) {
        avc_style_luma_interpolation_filter(
            searchRegionBuffer - (ME_FILTER_TAP >> 1) * luma_stride - (ME_FILTER_TAP >> 1) + 1,
            luma_stride,
            context_ptr->pos_b_buffer[list_index][ref_pic_index],
            context_ptr->interpolated_stride,
            search_area_width_for_asm,
            search_area_height + ME_FILTER_TAP,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            2);
    }

    // Half pel interpolation of the search region using f1 -> pos_h_buffer
    if (search_area_width_for_asm) {
        avc_style_luma_interpolation_filter(
            searchRegionBuffer - (ME_FILTER_TAP >> 1) * luma_stride - 1 + luma_stride,
            luma_stride,
            context_ptr->pos_h_buffer[list_index][ref_pic_index],
            context_ptr->interpolated_stride,
            search_area_width_for_asm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            8);
    }

    if (search_area_width_for_asm) {
        // Half pel interpolation of the search region using f1 -> pos_j_buffer
        avc_style_luma_interpolation_filter(
            context_ptr->pos_b_buffer[list_index][ref_pic_index] + context_ptr->interpolated_stride,
            context_ptr->interpolated_stride,
            context_ptr->pos_j_buffer[list_index][ref_pic_index],
            context_ptr->interpolated_stride,
            search_area_width_for_asm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            8);
    }

#else

    // Half pel interpolation of the search region using f1 -> pos_b_buffer
    HorizontalPelInterpolation(
        searchRegionBuffer - (ME_FILTER_TAP >> 1) * luma_stride - (ME_FILTER_TAP >> 1),
        luma_stride,
        search_area_width + 1,
        search_area_height + ME_FILTER_TAP,
        &(me_if_coeff[F1][0]),
        inputBitDepth,
        context_ptr->interpolated_stride,
        context_ptr->pos_b_buffer);

    // Half pel interpolation of the search region using f1 -> pos_h_buffer
    VerticalPelInterpolation(searchRegionBuffer - (ME_FILTER_TAP >> 1) * luma_stride - 1,
                             luma_stride,
                             search_area_width + 2,
                             search_area_height + 1,
                             &(me_if_coeff[F1][0]),
                             inputBitDepth,
                             context_ptr->interpolated_stride,
                             context_ptr->pos_h_buffer);

    // Half pel interpolation of the search region using f1 -> pos_j_buffer
    VerticalPelInterpolation(context_ptr->pos_b_buffer,
                             context_ptr->interpolated_stride,
                             search_area_width + 1,
                             search_area_height + 1,
                             &(me_if_coeff[F1][0]),
                             inputBitDepth,
                             context_ptr->interpolated_stride,
                             context_ptr->pos_j_buffer);

#endif

    return;
}

/*******************************************
 * pu_half_pel_refinement
 *   performs Half Pel refinement for one PU
 *******************************************/
static void pu_half_pel_refinement(
    SequenceControlSet *scs_ptr, // input parameter, Sequence control set Ptr
    MeContext *         context_ptr, // input parameter, ME context Ptr, used to get SB Ptr
    uint8_t *refBuffer, uint32_t ref_stride, uint32_t *p_best_Ssd,
    uint32_t blk_sb_buffer_index, // input parameter, PU origin, used to point to
    // source samples
    uint8_t *pos_b_buffer, // input parameter, position "b" interpolated search
    // area Ptr
    uint8_t *pos_h_buffer, // input parameter, position "h" interpolated search
    // area Ptr
    uint8_t *pos_j_buffer, // input parameter, position "j" interpolated search
    // area Ptr
    uint32_t pu_width, // input parameter, PU width
    uint32_t pu_height, // input parameter, PU height
    int16_t  x_search_area_origin, // input parameter, search area origin in the
    // horizontal direction, used to point to
    // reference samples
    int16_t y_search_area_origin, // input parameter, search area origin in the
    // vertical direction, used to point to
    // reference samples
    uint32_t *p_best_sad, uint32_t *p_best_MV, uint8_t *psubPelDirection) {
    EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;

    int32_t  search_region_index;
    uint64_t best_half_sad               = 0;
    uint64_t distortion_left_pos         = 0;
    uint64_t distortion_right_pos        = 0;
    uint64_t distortion_top_pos          = 0;
    uint64_t distortion_bottom_pos       = 0;
    uint64_t distortion_top_left_pos     = 0;
    uint64_t distortion_top_right_pos    = 0;
    uint64_t distortion_bottom_left_pos  = 0;
    uint64_t distortion_bottom_right_pos = 0;

    int16_t x_mv_half[8];
    int16_t y_mv_half[8];

    int16_t x_mv           = _MVXT(*p_best_MV);
    int16_t y_mv           = _MVYT(*p_best_MV);
    int32_t x_search_index = (x_mv >> 2) - x_search_area_origin;
    int32_t y_search_index = (y_mv >> 2) - y_search_area_origin;

    (void)scs_ptr;
    (void)encode_context_ptr;

    // TODO : remove these, and update the MV by just shifts

    x_mv_half[0] = x_mv - 2; // L  position
    x_mv_half[1] = x_mv + 2; // R  position
    x_mv_half[2] = x_mv; // T  position
    x_mv_half[3] = x_mv; // b  position
    x_mv_half[4] = x_mv - 2; // TL position
    x_mv_half[5] = x_mv + 2; // TR position
    x_mv_half[6] = x_mv + 2; // BR position
    x_mv_half[7] = x_mv - 2; // BL position

    y_mv_half[0] = y_mv; // L  position
    y_mv_half[1] = y_mv; // R  position
    y_mv_half[2] = y_mv - 2; // T  position
    y_mv_half[3] = y_mv + 2; // b  position
    y_mv_half[4] = y_mv - 2; // TL position
    y_mv_half[5] = y_mv - 2; // TR position
    y_mv_half[6] = y_mv + 2; // BR position
    y_mv_half[7] = y_mv + 2; // BL position

    // Compute SSD for the best full search candidate
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        *p_best_Ssd =
            (uint32_t)spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                     blk_sb_buffer_index,
                                                     context_ptr->sb_src_stride,
                                                     refBuffer,
                                                     y_search_index * (int32_t)ref_stride + x_search_index,
                                                     ref_stride,
                                                     pu_width,
                                                     pu_height);
    }
    // Use SATD only when QP mod, and RC are OFF
    // QP mod, and RC assume that ME distotion is always SAD.
    // This problem might be solved by computing SAD for the best position after
    // fractional search is done, or by considring the full pel resolution SAD.
    {
        // L position
        search_region_index =
            x_search_index + (int32_t)context_ptr->interpolated_stride * y_search_index;
        distortion_left_pos =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                 blk_sb_buffer_index,
                                                 context_ptr->sb_src_stride,
                                                 pos_b_buffer,
                                                 search_region_index,
                                                 context_ptr->interpolated_stride,
                                                 pu_width,
                                                 pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride << 1,
                                        &(pos_b_buffer[search_region_index]),
                                        context_ptr->interpolated_stride << 1,
                                        pu_height >> 1,
                                        pu_width))
                            << 1
                      : nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                       context_ptr->sb_src_stride,
                                       &(pos_b_buffer[search_region_index]),
                                       context_ptr->interpolated_stride,
                                       pu_height,
                                       pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortion_left_pos < *p_best_Ssd) {
                *p_best_sad =
                    (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                             context_ptr->sb_src_stride,
                                             &(pos_b_buffer[search_region_index]),
                                             context_ptr->interpolated_stride,
                                             pu_height,
                                             pu_width);
                *p_best_MV  = ((uint16_t)y_mv_half[0] << 16) | ((uint16_t)x_mv_half[0]);
                *p_best_Ssd = (uint32_t)distortion_left_pos;
            }
        } else {
            if (distortion_left_pos < *p_best_sad) {
                *p_best_sad = (uint32_t)distortion_left_pos;
                *p_best_MV  = ((uint16_t)y_mv_half[0] << 16) | ((uint16_t)x_mv_half[0]);
            }
        }
        // R position
        search_region_index++;
        distortion_right_pos =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                 blk_sb_buffer_index,
                                                 context_ptr->sb_src_stride,
                                                 pos_b_buffer,
                                                 search_region_index,
                                                 context_ptr->interpolated_stride,
                                                 pu_width,
                                                 pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride << 1,
                                        &(pos_b_buffer[search_region_index]),
                                        context_ptr->interpolated_stride << 1,
                                        pu_height >> 1,
                                        pu_width))
                            << 1
                      : nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                       context_ptr->sb_src_stride,
                                       &(pos_b_buffer[search_region_index]),
                                       context_ptr->interpolated_stride,
                                       pu_height,
                                       pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortion_right_pos < *p_best_Ssd) {
                *p_best_sad =
                    (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                             context_ptr->sb_src_stride,
                                             &(pos_b_buffer[search_region_index]),
                                             context_ptr->interpolated_stride,
                                             pu_height,
                                             pu_width);
                *p_best_MV  = ((uint16_t)y_mv_half[1] << 16) | ((uint16_t)x_mv_half[1]);
                *p_best_Ssd = (uint32_t)distortion_right_pos;
            }
        } else {
            if (distortion_right_pos < *p_best_sad) {
                *p_best_sad = (uint32_t)distortion_right_pos;
                *p_best_MV  = ((uint16_t)y_mv_half[1] << 16) | ((uint16_t)x_mv_half[1]);
            }
        }
        // T position
        search_region_index =
            x_search_index + (int32_t)context_ptr->interpolated_stride * y_search_index;
        distortion_top_pos =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                 blk_sb_buffer_index,
                                                 context_ptr->sb_src_stride,
                                                 pos_h_buffer,
                                                 search_region_index,
                                                 context_ptr->interpolated_stride,
                                                 pu_width,
                                                 pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride << 1,
                                        &(pos_h_buffer[search_region_index]),
                                        context_ptr->interpolated_stride << 1,
                                        pu_height >> 1,
                                        pu_width))
                            << 1
                      : nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                       context_ptr->sb_src_stride,
                                       &(pos_h_buffer[search_region_index]),
                                       context_ptr->interpolated_stride,
                                       pu_height,
                                       pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortion_top_pos < *p_best_Ssd) {
                *p_best_sad =
                    (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                             context_ptr->sb_src_stride,
                                             &(pos_h_buffer[search_region_index]),
                                             context_ptr->interpolated_stride,
                                             pu_height,
                                             pu_width);
                *p_best_MV  = ((uint16_t)y_mv_half[2] << 16) | ((uint16_t)x_mv_half[2]);
                *p_best_Ssd = (uint32_t)distortion_top_pos;
            }
        } else {
            if (distortion_top_pos < *p_best_sad) {
                *p_best_sad = (uint32_t)distortion_top_pos;
                *p_best_MV  = ((uint16_t)y_mv_half[2] << 16) | ((uint16_t)x_mv_half[2]);
            }
        }

        // b position
        search_region_index += (int32_t)context_ptr->interpolated_stride;
        distortion_bottom_pos =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                 blk_sb_buffer_index,
                                                 context_ptr->sb_src_stride,
                                                 pos_h_buffer,
                                                 search_region_index,
                                                 context_ptr->interpolated_stride,
                                                 pu_width,
                                                 pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride << 1,
                                        &(pos_h_buffer[search_region_index]),
                                        context_ptr->interpolated_stride << 1,
                                        pu_height >> 1,
                                        pu_width))
                            << 1
                      : nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                       context_ptr->sb_src_stride,
                                       &(pos_h_buffer[search_region_index]),
                                       context_ptr->interpolated_stride,
                                       pu_height,
                                       pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortion_bottom_pos < *p_best_Ssd) {
                *p_best_sad =
                    (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                             context_ptr->sb_src_stride,
                                             &(pos_h_buffer[search_region_index]),
                                             context_ptr->interpolated_stride,
                                             pu_height,
                                             pu_width);
                *p_best_MV  = ((uint16_t)y_mv_half[3] << 16) | ((uint16_t)x_mv_half[3]);
                *p_best_Ssd = (uint32_t)distortion_bottom_pos;
            }
        } else {
            if (distortion_bottom_pos < *p_best_sad) {
                *p_best_sad = (uint32_t)distortion_bottom_pos;
                *p_best_MV  = ((uint16_t)y_mv_half[3] << 16) | ((uint16_t)x_mv_half[3]);
            }
        }

        // TL position
        search_region_index =
            x_search_index + (int32_t)context_ptr->interpolated_stride * y_search_index;
        distortion_top_left_pos =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                 blk_sb_buffer_index,
                                                 context_ptr->sb_src_stride,
                                                 pos_j_buffer,
                                                 search_region_index,
                                                 context_ptr->interpolated_stride,
                                                 pu_width,
                                                 pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride << 1,
                                        &(pos_j_buffer[search_region_index]),
                                        context_ptr->interpolated_stride << 1,
                                        pu_height >> 1,
                                        pu_width))
                            << 1
                      : nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                       context_ptr->sb_src_stride,
                                       &(pos_j_buffer[search_region_index]),
                                       context_ptr->interpolated_stride,
                                       pu_height,
                                       pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortion_top_left_pos < *p_best_Ssd) {
                *p_best_sad =
                    (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                             context_ptr->sb_src_stride,
                                             &(pos_j_buffer[search_region_index]),
                                             context_ptr->interpolated_stride,
                                             pu_height,
                                             pu_width);
                *p_best_MV  = ((uint16_t)y_mv_half[4] << 16) | ((uint16_t)x_mv_half[4]);
                *p_best_Ssd = (uint32_t)distortion_top_left_pos;
            }
        } else {
            if (distortion_top_left_pos < *p_best_sad) {
                *p_best_sad = (uint32_t)distortion_top_left_pos;
                *p_best_MV  = ((uint16_t)y_mv_half[4] << 16) | ((uint16_t)x_mv_half[4]);
            }
        }

        // TR position
        search_region_index++;
        distortion_top_right_pos =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                 blk_sb_buffer_index,
                                                 context_ptr->sb_src_stride,
                                                 pos_j_buffer,
                                                 search_region_index,
                                                 context_ptr->interpolated_stride,
                                                 pu_width,
                                                 pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride << 1,
                                        &(pos_j_buffer[search_region_index]),
                                        context_ptr->interpolated_stride << 1,
                                        pu_height >> 1,
                                        pu_width))
                            << 1
                      : nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                       context_ptr->sb_src_stride,
                                       &(pos_j_buffer[search_region_index]),
                                       context_ptr->interpolated_stride,
                                       pu_height,
                                       pu_width);

        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortion_top_right_pos < *p_best_Ssd) {
                *p_best_sad =
                    (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                             context_ptr->sb_src_stride,
                                             &(pos_j_buffer[search_region_index]),
                                             context_ptr->interpolated_stride,
                                             pu_height,
                                             pu_width);
                *p_best_MV  = ((uint16_t)y_mv_half[5] << 16) | ((uint16_t)x_mv_half[5]);
                *p_best_Ssd = (uint32_t)distortion_top_right_pos;
            }
        } else {
            if (distortion_top_right_pos < *p_best_sad) {
                *p_best_sad = (uint32_t)distortion_top_right_pos;
                *p_best_MV  = ((uint16_t)y_mv_half[5] << 16) | ((uint16_t)x_mv_half[5]);
            }
        }

        // BR position
        search_region_index += (int32_t)context_ptr->interpolated_stride;
        distortion_bottom_right_pos =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                 blk_sb_buffer_index,
                                                 context_ptr->sb_src_stride,
                                                 pos_j_buffer,
                                                 search_region_index,
                                                 context_ptr->interpolated_stride,
                                                 pu_width,
                                                 pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride << 1,
                                        &(pos_j_buffer[search_region_index]),
                                        context_ptr->interpolated_stride << 1,
                                        pu_height >> 1,
                                        pu_width))
                            << 1
                      : nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                       context_ptr->sb_src_stride,
                                       &(pos_j_buffer[search_region_index]),
                                       context_ptr->interpolated_stride,
                                       pu_height,
                                       pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortion_bottom_right_pos < *p_best_Ssd) {
                *p_best_sad =
                    (uint32_t)nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                             context_ptr->sb_src_stride,
                                             &(pos_j_buffer[search_region_index]),
                                             context_ptr->interpolated_stride,
                                             pu_height,
                                             pu_width);
                *p_best_MV  = ((uint16_t)y_mv_half[6] << 16) | ((uint16_t)x_mv_half[6]);
                *p_best_Ssd = (uint32_t)distortion_bottom_right_pos;
            }
        } else {
            if (distortion_bottom_right_pos < *p_best_sad) {
                *p_best_sad = (uint32_t)distortion_bottom_right_pos;
                *p_best_MV  = ((uint16_t)y_mv_half[6] << 16) | ((uint16_t)x_mv_half[6]);
            }
        }

        // BL position
        search_region_index--;
        distortion_bottom_left_pos =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(context_ptr->sb_src_ptr,
                                                 blk_sb_buffer_index,
                                                 context_ptr->sb_src_stride,
                                                 pos_j_buffer,
                                                 search_region_index,
                                                 context_ptr->interpolated_stride,
                                                 pu_width,
                                                 pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride << 1,
                                        &(pos_j_buffer[search_region_index]),
                                        context_ptr->interpolated_stride << 1,
                                        pu_height >> 1,
                                        pu_width))
                            << 1
                      : (nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                        context_ptr->sb_src_stride,
                                        &(pos_j_buffer[search_region_index]),
                                        context_ptr->interpolated_stride,
                                        pu_height,
                                        pu_width));

        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortion_bottom_left_pos < *p_best_Ssd) {
                *p_best_sad =
                    (uint32_t)(nxm_sad_kernel(&(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
                                              context_ptr->sb_src_stride,
                                              &(pos_j_buffer[search_region_index]),
                                              context_ptr->interpolated_stride,
                                              pu_height,
                                              pu_width));
                *p_best_MV  = ((uint16_t)y_mv_half[7] << 16) | ((uint16_t)x_mv_half[7]);
                *p_best_Ssd = (uint32_t)distortion_bottom_left_pos;
            }
        } else {
            if (distortion_bottom_left_pos < *p_best_sad) {
                *p_best_sad = (uint32_t)distortion_bottom_left_pos;
                *p_best_MV  = ((uint16_t)y_mv_half[7] << 16) | ((uint16_t)x_mv_half[7]);
            }
        }
    }

    best_half_sad =
        MIN(distortion_left_pos,
            MIN(distortion_right_pos,
                MIN(distortion_top_pos,
                    MIN(distortion_bottom_pos,
                        MIN(distortion_top_left_pos,
                            MIN(distortion_top_right_pos,
                                MIN(distortion_bottom_left_pos, distortion_bottom_right_pos)))))));

    if (best_half_sad == distortion_left_pos)
        *psubPelDirection = LEFT_POSITION;
    else if (best_half_sad == distortion_right_pos)
        *psubPelDirection = RIGHT_POSITION;
    else if (best_half_sad == distortion_top_pos)
        *psubPelDirection = TOP_POSITION;
    else if (best_half_sad == distortion_bottom_pos)
        *psubPelDirection = BOTTOM_POSITION;
    else if (best_half_sad == distortion_top_left_pos)
        *psubPelDirection = TOP_LEFT_POSITION;
    else if (best_half_sad == distortion_top_right_pos)
        *psubPelDirection = TOP_RIGHT_POSITION;
    else if (best_half_sad == distortion_bottom_left_pos)
        *psubPelDirection = BOTTOM_LEFT_POSITION;
    else if (best_half_sad == distortion_bottom_right_pos)
        *psubPelDirection = BOTTOM_RIGHT_POSITION;
    return;
}

/*******************************************
 * half_pel_search_sb
 *   performs Half Pel refinement for the 85 PUs
 *******************************************/
void half_pel_search_sb(SequenceControlSet *scs_ptr, // input parameter, Sequence control set Ptr
                        PictureParentControlSet *pcs_ptr,
                        MeContext *context_ptr, // input/output parameter, ME context Ptr, used to
                        // get/update ME results
                        uint8_t *refBuffer, uint32_t ref_stride,
                        uint8_t *pos_b_buffer, // input parameter, position "b" interpolated search
                        // area Ptr
                        uint8_t *pos_h_buffer, // input parameter, position "h" interpolated search
                        // area Ptr
                        uint8_t *pos_j_buffer, // input parameter, position "j" interpolated search
                        // area Ptr
                        int16_t x_search_area_origin, // input parameter, search area origin in the
                        // horizontal direction, used to point to
                        // reference samples
                        int16_t y_search_area_origin, // input parameter, search area origin in the
                        // vertical direction, used to point to
                        // reference samples
                        EbBool enable_half_pel_32x32,
                        EbBool enable_half_pel_16x16, EbBool enable_half_pel_8x8) {
    uint32_t idx;
    uint32_t pu_index;
    uint32_t pu_shift_x_index;
    uint32_t pu_shift_y_index;
    uint32_t blk_sb_buffer_index;
    uint32_t pos_b_buffer_index;
    uint32_t pos_h_buffer_index;
    uint32_t pos_j_buffer_index;

    // 64x64 [1 partition]
    pu_half_pel_refinement(scs_ptr,
                            context_ptr,
                            &(refBuffer[0]),
                            ref_stride,
                            context_ptr->p_best_ssd64x64,
                            0,
                            &(pos_b_buffer[0]),
                            &(pos_h_buffer[0]),
                            &(pos_j_buffer[0]),
                            64,
                            64,
                            x_search_area_origin,
                            y_search_area_origin,
                            context_ptr->p_best_sad_64x64,
                            context_ptr->p_best_mv64x64,
                            &context_ptr->psub_pel_direction64x64);

    if (enable_half_pel_32x32) {
        // 32x32 [4 partitions]
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            pu_shift_x_index = (pu_index & 0x01) << 5;
            pu_shift_y_index = (pu_index >> 1) << 5;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;
            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd32x32[pu_index],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   32,
                                   32,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_32x32[pu_index],
                                   &context_ptr->p_best_mv32x32[pu_index],
                                   &context_ptr->psub_pel_direction32x32[pu_index]);
        }
    }
    if (enable_half_pel_16x16) {
        // 16x16 [16 partitions]
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            idx = tab16x16[pu_index];

            pu_shift_x_index = (pu_index & 0x03) << 4;
            pu_shift_y_index = (pu_index >> 2) << 4;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;
            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd16x16[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   16,
                                   16,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_16x16[idx],
                                   &context_ptr->p_best_mv16x16[idx],
                                   &context_ptr->psub_pel_direction16x16[idx]);
        }
    }
    if (enable_half_pel_8x8) {
        // 8x8   [64 partitions]
        for (pu_index = 0; pu_index < 64; ++pu_index) {
            idx = tab8x8[pu_index]; // TODO bitwise this

            pu_shift_x_index = (pu_index & 0x07) << 3;
            pu_shift_y_index = (pu_index >> 3) << 3;

            blk_sb_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(
                scs_ptr,
                context_ptr,
                &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                ref_stride,
                &context_ptr->p_best_ssd8x8[idx],
                blk_sb_buffer_index,
                &(pos_b_buffer[pos_b_buffer_index]),
                &(pos_h_buffer[pos_h_buffer_index]),
                &(pos_j_buffer[pos_j_buffer_index]),
                8,
                8,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad_8x8[idx],
                &context_ptr->p_best_mv8x8[idx],
                &context_ptr->psub_pel_direction8x8[idx]);
        }
    }
    if (pcs_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
        // 64x32
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            pu_shift_x_index = 0;
            pu_shift_y_index = pu_index << 5;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd64x32[pu_index],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   64,
                                   32,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_64x32[pu_index],
                                   &context_ptr->p_best_mv64x32[pu_index],
                                   &context_ptr->psub_pel_direction64x32[pu_index]);
        }

        // 32x16
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            idx = tab32x16[pu_index]; // TODO bitwise this

            pu_shift_x_index = (pu_index & 0x01) << 5;
            pu_shift_y_index = (pu_index >> 1) << 4;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd32x16[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   32,
                                   16,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_32x16[idx],
                                   &context_ptr->p_best_mv32x16[idx],
                                   &context_ptr->psub_pel_direction32x16[idx]);
        }

        // 16x8
        for (pu_index = 0; pu_index < 32; ++pu_index) {
            idx = tab16x8[pu_index];

            pu_shift_x_index = (pu_index & 0x03) << 4;
            pu_shift_y_index = (pu_index >> 2) << 3;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd16x8[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   16,
                                   8,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_16x8[idx],
                                   &context_ptr->p_best_mv16x8[idx],
                                   &context_ptr->psub_pel_direction16x8[idx]);
        }

        // 32x64
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            pu_shift_x_index = pu_index << 5;
            pu_shift_y_index = 0;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd32x64[pu_index],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   32,
                                   64,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_32x64[pu_index],
                                   &context_ptr->p_best_mv32x64[pu_index],
                                   &context_ptr->psub_pel_direction32x64[pu_index]);
        }

        // 16x32
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            idx = tab16x32[pu_index];

            pu_shift_x_index = (pu_index & 0x03) << 4;
            pu_shift_y_index = (pu_index >> 2) << 5;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd16x32[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   16,
                                   32,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_16x32[idx],
                                   &context_ptr->p_best_mv16x32[idx],
                                   &context_ptr->psub_pel_direction16x32[idx]);
        }

        // 8x16
        for (pu_index = 0; pu_index < 32; ++pu_index) {
            idx = tab8x16[pu_index];

            pu_shift_x_index = (pu_index & 0x07) << 3;
            pu_shift_y_index = (pu_index >> 3) << 4;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd8x16[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   8,
                                   16,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_8x16[idx],
                                   &context_ptr->p_best_mv8x16[idx],
                                   &context_ptr->psub_pel_direction8x16[idx]);
        }

        // 32x8
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            idx = tab32x8[pu_index];

            pu_shift_x_index = (pu_index & 0x01) << 5;
            pu_shift_y_index = (pu_index >> 1) << 3;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd32x8[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   32,
                                   8,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_32x8[idx],
                                   &context_ptr->p_best_mv32x8[idx],
                                   &context_ptr->psub_pel_direction32x8[idx]);
        }

        for (pu_index = 0; pu_index < 16; ++pu_index) {
            idx = tab8x32[pu_index];

            pu_shift_x_index = (pu_index & 0x07) << 3;
            pu_shift_y_index = (pu_index >> 3) << 5;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;
            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd8x32[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   8,
                                   32,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_8x32[idx],
                                   &context_ptr->p_best_mv8x32[idx],
                                   &context_ptr->psub_pel_direction8x32[idx]);
        }

        for (pu_index = 0; pu_index < 4; ++pu_index) {
            idx = pu_index;

            pu_shift_x_index = 0;
            pu_shift_y_index = pu_index << 4;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;
            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd64x16[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   64,
                                   16,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_64x16[idx],
                                   &context_ptr->p_best_mv64x16[idx],
                                   &context_ptr->psub_pel_direction64x16[idx]);
        }

        for (pu_index = 0; pu_index < 4; ++pu_index) {
            idx = pu_index;

            pu_shift_x_index = pu_index << 4;
            pu_shift_y_index = 0;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;
            pos_b_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_h_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;
            pos_j_buffer_index =
                pu_shift_x_index + pu_shift_y_index * context_ptr->interpolated_stride;

            pu_half_pel_refinement(scs_ptr,
                                   context_ptr,
                                   &(refBuffer[pu_shift_y_index * ref_stride + pu_shift_x_index]),
                                   ref_stride,
                                   &context_ptr->p_best_ssd16x64[idx],
                                   blk_sb_buffer_index,
                                   &(pos_b_buffer[pos_b_buffer_index]),
                                   &(pos_h_buffer[pos_h_buffer_index]),
                                   &(pos_j_buffer[pos_j_buffer_index]),
                                   16,
                                   64,
                                   x_search_area_origin,
                                   y_search_area_origin,
                                   &context_ptr->p_best_sad_16x64[idx],
                                   &context_ptr->p_best_mv16x64[idx],
                                   &context_ptr->psub_pel_direction16x64[idx]);
        }
    }

    return;
}
/*******************************************
 * combined_averaging_ssd
 *
 *******************************************/
uint32_t combined_averaging_ssd_c(uint8_t *src, ptrdiff_t src_stride, uint8_t *ref1,
                                  ptrdiff_t ref1_stride, uint8_t *ref2, ptrdiff_t ref2_stride,
                                  uint32_t height, uint32_t width) {
    uint32_t x, y;
    uint32_t ssd = 0;
    uint8_t  avgpel;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            avgpel = (ref1[x] + ref2[x] + 1) >> 1;
            ssd += SQR((int64_t)(src[x]) - (avgpel));
        }
        src += src_stride;
        ref1 += ref1_stride;
        ref2 += ref2_stride;
    }
    return ssd;
}
/*******************************************
 * pu_quarter_pel_refinement_on_the_fly
 *   performs Quarter Pel refinement for each PU
 *******************************************/
static void pu_quarter_pel_refinement_on_the_fly(
    MeContext *context_ptr, // [IN] ME context Ptr, used to get SB Ptr
    uint32_t * p_best_Ssd,
    uint32_t   blk_sb_buffer_index, // [IN] PU origin, used to point to source samples
    uint8_t ** buf1, // [IN]
    uint32_t * buf_1_stride,
    uint8_t ** buf2, // [IN]
    uint32_t * buf_2_stride,
    uint32_t   pu_width, // [IN]  PU width
    uint32_t   pu_height, // [IN]  PU height
    int16_t    x_search_area_origin, // [IN] search area origin in the horizontal
    // direction, used to point to reference samples
    int16_t y_search_area_origin, // [IN] search area origin in the vertical
    // direction, used to point to reference samples
    uint32_t *p_best_sad, uint32_t *p_best_MV, uint8_t sub_pel_direction) {
    int16_t x_mv = _MVXT(*p_best_MV);
    int16_t y_mv = _MVYT(*p_best_MV);

    int16_t x_search_index = ((x_mv + 2) >> 2) - x_search_area_origin;
    int16_t y_search_index = ((y_mv + 2) >> 2) - y_search_area_origin;

    uint64_t dist;

    EbBool valid_tl, valid_t, valid_tr, valid_r, valid_br, valid_b, valid_bl, valid_l;

    int16_t x_mv_quarter[8];
    int16_t y_mv_quarter[8];
    int32_t search_region_index_1 = 0;
    int32_t search_region_index_2 = 0;
    if (context_ptr->full_quarter_pel_refinement) {
        valid_tl = EB_TRUE;
        valid_t  = EB_TRUE;
        valid_tr = EB_TRUE;
        valid_r  = EB_TRUE;
        valid_br = EB_TRUE;
        valid_b  = EB_TRUE;
        valid_bl = EB_TRUE;
        valid_l  = EB_TRUE;
    } else {
        if ((y_mv & 2) + ((x_mv & 2) >> 1)) {
            valid_tl = (EbBool)(sub_pel_direction == RIGHT_POSITION ||
                                sub_pel_direction == BOTTOM_RIGHT_POSITION ||
                                sub_pel_direction == BOTTOM_POSITION);
            valid_t  = (EbBool)(sub_pel_direction == BOTTOM_RIGHT_POSITION ||
                               sub_pel_direction == BOTTOM_POSITION ||
                               sub_pel_direction == BOTTOM_LEFT_POSITION);
            valid_tr = (EbBool)(sub_pel_direction == BOTTOM_POSITION ||
                                sub_pel_direction == BOTTOM_LEFT_POSITION ||
                                sub_pel_direction == LEFT_POSITION);
            valid_r  = (EbBool)(sub_pel_direction == BOTTOM_LEFT_POSITION ||
                               sub_pel_direction == LEFT_POSITION ||
                               sub_pel_direction == TOP_LEFT_POSITION);
            valid_br = (EbBool)(sub_pel_direction == LEFT_POSITION ||
                                sub_pel_direction == TOP_LEFT_POSITION ||
                                sub_pel_direction == TOP_POSITION);
            valid_b  = (EbBool)(sub_pel_direction == TOP_LEFT_POSITION ||
                               sub_pel_direction == TOP_POSITION ||
                               sub_pel_direction == TOP_RIGHT_POSITION);
            valid_bl = (EbBool)(sub_pel_direction == TOP_POSITION ||
                                sub_pel_direction == TOP_RIGHT_POSITION ||
                                sub_pel_direction == RIGHT_POSITION);
            valid_l  = (EbBool)(sub_pel_direction == TOP_RIGHT_POSITION ||
                               sub_pel_direction == RIGHT_POSITION ||
                               sub_pel_direction == BOTTOM_RIGHT_POSITION);
        } else {
            valid_tl = (EbBool)(sub_pel_direction == LEFT_POSITION ||
                                sub_pel_direction == TOP_LEFT_POSITION ||
                                sub_pel_direction == TOP_POSITION);
            valid_t  = (EbBool)(sub_pel_direction == TOP_LEFT_POSITION ||
                               sub_pel_direction == TOP_POSITION ||
                               sub_pel_direction == TOP_RIGHT_POSITION);
            valid_tr = (EbBool)(sub_pel_direction == TOP_POSITION ||
                                sub_pel_direction == TOP_RIGHT_POSITION ||
                                sub_pel_direction == RIGHT_POSITION);
            valid_r  = (EbBool)(sub_pel_direction == TOP_RIGHT_POSITION ||
                               sub_pel_direction == RIGHT_POSITION ||
                               sub_pel_direction == BOTTOM_RIGHT_POSITION);
            valid_br = (EbBool)(sub_pel_direction == RIGHT_POSITION ||
                                sub_pel_direction == BOTTOM_RIGHT_POSITION ||
                                sub_pel_direction == BOTTOM_POSITION);
            valid_b  = (EbBool)(sub_pel_direction == BOTTOM_RIGHT_POSITION ||
                               sub_pel_direction == BOTTOM_POSITION ||
                               sub_pel_direction == BOTTOM_LEFT_POSITION);
            valid_bl = (EbBool)(sub_pel_direction == BOTTOM_POSITION ||
                                sub_pel_direction == BOTTOM_LEFT_POSITION ||
                                sub_pel_direction == LEFT_POSITION);
            valid_l  = (EbBool)(sub_pel_direction == BOTTOM_LEFT_POSITION ||
                               sub_pel_direction == LEFT_POSITION ||
                               sub_pel_direction == TOP_LEFT_POSITION);
        }
    }
    x_mv_quarter[0] = x_mv - 1; // L  position
    x_mv_quarter[1] = x_mv + 1; // R  position
    x_mv_quarter[2] = x_mv; // T  position
    x_mv_quarter[3] = x_mv; // b  position
    x_mv_quarter[4] = x_mv - 1; // TL position
    x_mv_quarter[5] = x_mv + 1; // TR position
    x_mv_quarter[6] = x_mv + 1; // BR position
    x_mv_quarter[7] = x_mv - 1; // BL position

    y_mv_quarter[0] = y_mv; // L  position
    y_mv_quarter[1] = y_mv; // R  position
    y_mv_quarter[2] = y_mv - 1; // T  position
    y_mv_quarter[3] = y_mv + 1; // b  position
    y_mv_quarter[4] = y_mv - 1; // TL position
    y_mv_quarter[5] = y_mv - 1; // TR position
    y_mv_quarter[6] = y_mv + 1; // BR position
    y_mv_quarter[7] = y_mv + 1; // BL position

    // Use SATD only when QP mod, and RC are OFF
    // QP mod, and RC assume that ME distotion is always SAD.
    // This problem might be solved by computing SAD for the best position after
    // fractional search is done, or by considring the full pel resolution SAD.

    {
        // L position
        if (valid_l) {
            search_region_index_1 =
                (int32_t)x_search_index + (int32_t)buf_1_stride[0] * (int32_t)y_search_index;
            search_region_index_2 =
                (int32_t)x_search_index + (int32_t)buf_2_stride[0] * (int32_t)y_search_index;

            dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                       ? combined_averaging_ssd(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                BLOCK_SIZE_64,
                                                buf1[0] + search_region_index_1,
                                                buf_1_stride[0],
                                                buf2[0] + search_region_index_2,
                                                buf_2_stride[0],
                                                pu_height,
                                                pu_width)
                       : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                             ? (nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                   BLOCK_SIZE_64 << 1,
                                                   buf1[0] + search_region_index_1,
                                                   buf_1_stride[0] << 1,
                                                   buf2[0] + search_region_index_2,
                                                   buf_2_stride[0] << 1,
                                                   pu_height >> 1,
                                                   pu_width))
                                   << 1
                             : nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                  BLOCK_SIZE_64,
                                                  buf1[0] + search_region_index_1,
                                                  buf_1_stride[0],
                                                  buf2[0] + search_region_index_2,
                                                  buf_2_stride[0],
                                                  pu_height,
                                                  pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *p_best_Ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                     BLOCK_SIZE_64,
                                                     buf1[0] + search_region_index_1,
                                                     buf_1_stride[0],
                                                     buf2[0] + search_region_index_2,
                                                     buf_2_stride[0],
                                                     pu_height,
                                                     pu_width);
                    *p_best_MV  = ((uint16_t)y_mv_quarter[0] << 16) | ((uint16_t)x_mv_quarter[0]);
                    *p_best_Ssd = (uint32_t)dist;
                }
            } else {
                if (dist < *p_best_sad) {
                    *p_best_sad = (uint32_t)dist;
                    *p_best_MV  = ((uint16_t)y_mv_quarter[0] << 16) | ((uint16_t)x_mv_quarter[0]);
                }
            }
        }

        // R positions
        if (valid_r) {
            search_region_index_1 =
                (int32_t)x_search_index + (int32_t)buf_1_stride[1] * (int32_t)y_search_index;
            search_region_index_2 =
                (int32_t)x_search_index + (int32_t)buf_2_stride[1] * (int32_t)y_search_index;
            dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                       ? combined_averaging_ssd(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                BLOCK_SIZE_64,
                                                buf1[1] + search_region_index_1,
                                                buf_1_stride[1],
                                                buf2[1] + search_region_index_2,
                                                buf_2_stride[1],
                                                pu_height,
                                                pu_width)
                       : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                             ? (nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                   BLOCK_SIZE_64 << 1,
                                                   buf1[1] + search_region_index_1,
                                                   buf_1_stride[1] << 1,
                                                   buf2[1] + search_region_index_2,
                                                   buf_2_stride[1] << 1,
                                                   pu_height >> 1,
                                                   pu_width))
                                   << 1
                             : nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                  BLOCK_SIZE_64,
                                                  buf1[1] + search_region_index_1,
                                                  buf_1_stride[1],
                                                  buf2[1] + search_region_index_2,
                                                  buf_2_stride[1],
                                                  pu_height,
                                                  pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *p_best_Ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                     BLOCK_SIZE_64,
                                                     buf1[1] + search_region_index_1,
                                                     buf_1_stride[1],
                                                     buf2[1] + search_region_index_2,
                                                     buf_2_stride[1],
                                                     pu_height,
                                                     pu_width);
                    *p_best_MV  = ((uint16_t)y_mv_quarter[1] << 16) | ((uint16_t)x_mv_quarter[1]);
                    *p_best_Ssd = (uint32_t)dist;
                }
            } else {
                if (dist < *p_best_sad) {
                    *p_best_sad = (uint32_t)dist;
                    *p_best_MV  = ((uint16_t)y_mv_quarter[1] << 16) | ((uint16_t)x_mv_quarter[1]);
                }
            }
        }

        // T position
        if (valid_t) {
            search_region_index_1 =
                (int32_t)x_search_index + (int32_t)buf_1_stride[2] * (int32_t)y_search_index;
            search_region_index_2 =
                (int32_t)x_search_index + (int32_t)buf_2_stride[2] * (int32_t)y_search_index;
            dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                       ? combined_averaging_ssd(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                BLOCK_SIZE_64,
                                                buf1[2] + search_region_index_1,
                                                buf_1_stride[2],
                                                buf2[2] + search_region_index_2,
                                                buf_2_stride[2],
                                                pu_height,
                                                pu_width)
                       : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                             ? (nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                   BLOCK_SIZE_64 << 1,
                                                   buf1[2] + search_region_index_1,
                                                   buf_1_stride[2] << 1,
                                                   buf2[2] + search_region_index_2,
                                                   buf_2_stride[2] << 1,
                                                   pu_height >> 1,
                                                   pu_width))
                                   << 1
                             : nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                  BLOCK_SIZE_64,
                                                  buf1[2] + search_region_index_1,
                                                  buf_1_stride[2],
                                                  buf2[2] + search_region_index_2,
                                                  buf_2_stride[2],
                                                  pu_height,
                                                  pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *p_best_Ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                     BLOCK_SIZE_64,
                                                     buf1[2] + search_region_index_1,
                                                     buf_1_stride[2],
                                                     buf2[2] + search_region_index_2,
                                                     buf_2_stride[2],
                                                     pu_height,
                                                     pu_width);
                    *p_best_MV  = ((uint16_t)y_mv_quarter[2] << 16) | ((uint16_t)x_mv_quarter[2]);
                    *p_best_Ssd = (uint32_t)dist;
                }
            } else {
                if (dist < *p_best_sad) {
                    *p_best_sad = (uint32_t)dist;
                    *p_best_MV  = ((uint16_t)y_mv_quarter[2] << 16) | ((uint16_t)x_mv_quarter[2]);
                }
            }
        }

        // b position
        if (valid_b) {
            search_region_index_1 =
                (int32_t)x_search_index + (int32_t)buf_1_stride[3] * (int32_t)y_search_index;
            search_region_index_2 =
                (int32_t)x_search_index + (int32_t)buf_2_stride[3] * (int32_t)y_search_index;
            dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                       ? combined_averaging_ssd(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                BLOCK_SIZE_64,
                                                buf1[3] + search_region_index_1,
                                                buf_1_stride[3],
                                                buf2[3] + search_region_index_2,
                                                buf_2_stride[3],
                                                pu_height,
                                                pu_width)
                       : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                             ? (nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                   BLOCK_SIZE_64 << 1,
                                                   buf1[3] + search_region_index_1,
                                                   buf_1_stride[3] << 1,
                                                   buf2[3] + search_region_index_2,
                                                   buf_2_stride[3] << 1,
                                                   pu_height >> 1,
                                                   pu_width))
                                   << 1
                             : nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                  BLOCK_SIZE_64,
                                                  buf1[3] + search_region_index_1,
                                                  buf_1_stride[3],
                                                  buf2[3] + search_region_index_2,
                                                  buf_2_stride[3],
                                                  pu_height,
                                                  pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *p_best_Ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                     BLOCK_SIZE_64,
                                                     buf1[3] + search_region_index_1,
                                                     buf_1_stride[3],
                                                     buf2[3] + search_region_index_2,
                                                     buf_2_stride[3],
                                                     pu_height,
                                                     pu_width);
                    *p_best_MV  = ((uint16_t)y_mv_quarter[3] << 16) | ((uint16_t)x_mv_quarter[3]);
                    *p_best_Ssd = (uint32_t)dist;
                }
            } else {
                if (dist < *p_best_sad) {
                    *p_best_sad = (uint32_t)dist;
                    *p_best_MV  = ((uint16_t)y_mv_quarter[3] << 16) | ((uint16_t)x_mv_quarter[3]);
                }
            }
        }

        // TL position
        if (valid_tl) {
            search_region_index_1 =
                (int32_t)x_search_index + (int32_t)buf_1_stride[4] * (int32_t)y_search_index;
            search_region_index_2 =
                (int32_t)x_search_index + (int32_t)buf_2_stride[4] * (int32_t)y_search_index;
            dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                       ? combined_averaging_ssd(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                BLOCK_SIZE_64,
                                                buf1[4] + search_region_index_1,
                                                buf_1_stride[4],
                                                buf2[4] + search_region_index_2,
                                                buf_2_stride[4],
                                                pu_height,
                                                pu_width)
                       : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                             ? (nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                   BLOCK_SIZE_64 << 1,
                                                   buf1[4] + search_region_index_1,
                                                   buf_1_stride[4] << 1,
                                                   buf2[4] + search_region_index_2,
                                                   buf_2_stride[4] << 1,
                                                   pu_height >> 1,
                                                   pu_width))
                                   << 1
                             : nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                  BLOCK_SIZE_64,
                                                  buf1[4] + search_region_index_1,
                                                  buf_1_stride[4],
                                                  buf2[4] + search_region_index_2,
                                                  buf_2_stride[4],
                                                  pu_height,
                                                  pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *p_best_Ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                     BLOCK_SIZE_64,
                                                     buf1[4] + search_region_index_1,
                                                     buf_1_stride[4],
                                                     buf2[4] + search_region_index_2,
                                                     buf_2_stride[4],
                                                     pu_height,
                                                     pu_width);
                    *p_best_MV  = ((uint16_t)y_mv_quarter[4] << 16) | ((uint16_t)x_mv_quarter[4]);
                    *p_best_Ssd = (uint32_t)dist;
                }
            } else {
                if (dist < *p_best_sad) {
                    *p_best_sad = (uint32_t)dist;
                    *p_best_MV  = ((uint16_t)y_mv_quarter[4] << 16) | ((uint16_t)x_mv_quarter[4]);
                }
            }
        }

        // TR position
        if (valid_tr) {
            search_region_index_1 =
                (int32_t)x_search_index + (int32_t)buf_1_stride[5] * (int32_t)y_search_index;
            search_region_index_2 =
                (int32_t)x_search_index + (int32_t)buf_2_stride[5] * (int32_t)y_search_index;
            dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                       ? combined_averaging_ssd(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                BLOCK_SIZE_64,
                                                buf1[5] + search_region_index_1,
                                                buf_1_stride[5],
                                                buf2[5] + search_region_index_2,
                                                buf_2_stride[5],
                                                pu_height,
                                                pu_width)
                       : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                             ? (nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                   BLOCK_SIZE_64 << 1,
                                                   buf1[5] + search_region_index_1,
                                                   buf_1_stride[5] << 1,
                                                   buf2[5] + search_region_index_2,
                                                   buf_2_stride[5] << 1,
                                                   pu_height >> 1,
                                                   pu_width))
                                   << 1
                             : nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                  BLOCK_SIZE_64,
                                                  buf1[5] + search_region_index_1,
                                                  buf_1_stride[5],
                                                  buf2[5] + search_region_index_2,
                                                  buf_2_stride[5],
                                                  pu_height,
                                                  pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *p_best_Ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                     BLOCK_SIZE_64,
                                                     buf1[5] + search_region_index_1,
                                                     buf_1_stride[5],
                                                     buf2[5] + search_region_index_2,
                                                     buf_2_stride[5],
                                                     pu_height,
                                                     pu_width);
                    *p_best_MV  = ((uint16_t)y_mv_quarter[5] << 16) | ((uint16_t)x_mv_quarter[5]);
                    *p_best_Ssd = (uint32_t)dist;
                }
            } else {
                if (dist < *p_best_sad) {
                    *p_best_sad = (uint32_t)dist;
                    *p_best_MV  = ((uint16_t)y_mv_quarter[5] << 16) | ((uint16_t)x_mv_quarter[5]);
                }
            }
        }

        // BR position
        if (valid_br) {
            search_region_index_1 =
                (int32_t)x_search_index + (int32_t)buf_1_stride[6] * (int32_t)y_search_index;
            search_region_index_2 =
                (int32_t)x_search_index + (int32_t)buf_2_stride[6] * (int32_t)y_search_index;
            dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                       ? combined_averaging_ssd(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                BLOCK_SIZE_64,
                                                buf1[6] + search_region_index_1,
                                                buf_1_stride[6],
                                                buf2[6] + search_region_index_2,
                                                buf_2_stride[6],
                                                pu_height,
                                                pu_width)
                       : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                             ? (nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                   BLOCK_SIZE_64 << 1,
                                                   buf1[6] + search_region_index_1,
                                                   buf_1_stride[6] << 1,
                                                   buf2[6] + search_region_index_2,
                                                   buf_2_stride[6] << 1,
                                                   pu_height >> 1,
                                                   pu_width))
                                   << 1
                             : nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                  BLOCK_SIZE_64,
                                                  buf1[6] + search_region_index_1,
                                                  buf_1_stride[6],
                                                  buf2[6] + search_region_index_2,
                                                  buf_2_stride[6],
                                                  pu_height,
                                                  pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *p_best_Ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                     BLOCK_SIZE_64,
                                                     buf1[6] + search_region_index_1,
                                                     buf_1_stride[6],
                                                     buf2[6] + search_region_index_2,
                                                     buf_2_stride[6],
                                                     pu_height,
                                                     pu_width);
                    *p_best_MV  = ((uint16_t)y_mv_quarter[6] << 16) | ((uint16_t)x_mv_quarter[6]);
                    *p_best_Ssd = (uint32_t)dist;
                }
            } else {
                if (dist < *p_best_sad) {
                    *p_best_sad = (uint32_t)dist;
                    *p_best_MV  = ((uint16_t)y_mv_quarter[6] << 16) | ((uint16_t)x_mv_quarter[6]);
                }
            }
        }

        // BL position
        if (valid_bl) {
            search_region_index_1 =
                (int32_t)x_search_index + (int32_t)buf_1_stride[7] * (int32_t)y_search_index;
            search_region_index_2 =
                (int32_t)x_search_index + (int32_t)buf_2_stride[7] * (int32_t)y_search_index;
            dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                       ? combined_averaging_ssd(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                BLOCK_SIZE_64,
                                                buf1[7] + search_region_index_1,
                                                buf_1_stride[7],
                                                buf2[7] + search_region_index_2,
                                                buf_2_stride[7],
                                                pu_height,
                                                pu_width)
                       : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                             ? (nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                   BLOCK_SIZE_64 << 1,
                                                   buf1[7] + search_region_index_1,
                                                   buf_1_stride[7] << 1,
                                                   buf2[7] + search_region_index_2,
                                                   buf_2_stride[7] << 1,
                                                   pu_height >> 1,
                                                   pu_width))
                                   << 1
                             : nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                  BLOCK_SIZE_64,
                                                  buf1[7] + search_region_index_1,
                                                  buf_1_stride[7],
                                                  buf2[7] + search_region_index_2,
                                                  buf_2_stride[7],
                                                  pu_height,
                                                  pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *p_best_Ssd) {
                    *p_best_sad =
                        (uint32_t)nxm_sad_avg_kernel(&(context_ptr->sb_buffer[blk_sb_buffer_index]),
                                                     BLOCK_SIZE_64,
                                                     buf1[7] + search_region_index_1,
                                                     buf_1_stride[7],
                                                     buf2[7] + search_region_index_2,
                                                     buf_2_stride[7],
                                                     pu_height,
                                                     pu_width);
                    *p_best_MV  = ((uint16_t)y_mv_quarter[7] << 16) | ((uint16_t)x_mv_quarter[7]);
                    *p_best_Ssd = (uint32_t)dist;
                }
            } else {
                if (dist < *p_best_sad) {
                    *p_best_sad = (uint32_t)dist;
                    *p_best_MV  = ((uint16_t)y_mv_quarter[7] << 16) | ((uint16_t)x_mv_quarter[7]);
                }
            }
        }
    }

    return;
}

/*******************************************
* set_quarter_pel_refinement_inputs_on_the_fly
*   determine the 2 half pel buffers to do
averaging for Quarter Pel Refinement
*******************************************/
static void set_quarter_pel_refinement_inputs_on_the_fly(uint8_t * pos_Full, //[IN] points to A
                                                         uint32_t  FullStride, //[IN]
                                                         uint8_t * pos_b, //[IN] points to b
                                                         uint8_t * pos_h, //[IN] points to h
                                                         uint8_t * pos_j, //[IN] points to j
                                                         uint32_t  Stride, //[IN]
                                                         int16_t   x_mv, //[IN]
                                                         int16_t   y_mv, //[IN]
                                                         uint8_t **buf1, //[OUT]
                                                         uint32_t *buf_1_stride, //[OUT]
                                                         uint8_t **buf2, //[OUT]
                                                         uint32_t *buf_2_stride //[OUT]
) {
    uint32_t quarter_pel_refinement_method = (y_mv & 2) + ((x_mv & 2) >> 1);

    // for each one of the 8 postions, we need to determine the 2 half pel
    // buffers to  do averaging

    //     A    a    b    c
    //     d    e    f    g
    //     h    i    j    k
    //     n    p    q    r

    switch (quarter_pel_refinement_method) {
    case EB_QUARTER_IN_FULL:

        /*c=b+A*/ buf1[0] = pos_b;
        buf_1_stride[0]   = Stride;
        buf2[0]           = pos_Full;
        buf_2_stride[0]   = FullStride;
        /*a=A+b*/ buf1[1] = pos_Full;
        buf_1_stride[1]   = FullStride;
        buf2[1]           = pos_b + 1;
        buf_2_stride[1]   = Stride;
        /*n=h+A*/ buf1[2] = pos_h;
        buf_1_stride[2]   = Stride;
        buf2[2]           = pos_Full;
        buf_2_stride[2]   = FullStride;
        /*d=A+h*/ buf1[3] = pos_Full;
        buf_1_stride[3]   = FullStride;
        buf2[3]           = pos_h + Stride;
        buf_2_stride[3]   = Stride;
        /*r=b+h*/ buf1[4] = pos_b;
        buf_1_stride[4]   = Stride;
        buf2[4]           = pos_h;
        buf_2_stride[4]   = Stride;
        /*p=h+b*/ buf1[5] = pos_h;
        buf_1_stride[5]   = Stride;
        buf2[5]           = pos_b + 1;
        buf_2_stride[5]   = Stride;
        /*e=h+b*/ buf1[6] = pos_h + Stride;
        buf_1_stride[6]   = Stride;
        buf2[6]           = pos_b + 1;
        buf_2_stride[6]   = Stride;
        /*g=b+h*/ buf1[7] = pos_b;
        buf_1_stride[7]   = Stride;
        buf2[7]           = pos_h + Stride;
        buf_2_stride[7]   = Stride;

        break;

    case EB_QUARTER_IN_HALF_HORIZONTAL:

        /*a=A+b*/ buf1[0] = pos_Full - 1;
        buf_1_stride[0]   = FullStride;
        buf2[0]           = pos_b;
        buf_2_stride[0]   = Stride;
        /*c=b+A*/ buf1[1] = pos_b;
        buf_1_stride[1]   = Stride;
        buf2[1]           = pos_Full;
        buf_2_stride[1]   = FullStride;
        /*q=j+b*/ buf1[2] = pos_j;
        buf_1_stride[2]   = Stride;
        buf2[2]           = pos_b;
        buf_2_stride[2]   = Stride;
        /*f=b+j*/ buf1[3] = pos_b;
        buf_1_stride[3]   = Stride;
        buf2[3]           = pos_j + Stride;
        buf_2_stride[3]   = Stride;
        /*p=h+b*/ buf1[4] = pos_h - 1;
        buf_1_stride[4]   = Stride;
        buf2[4]           = pos_b;
        buf_2_stride[4]   = Stride;
        /*r=b+h*/ buf1[5] = pos_b;
        buf_1_stride[5]   = Stride;
        buf2[5]           = pos_h;
        buf_2_stride[5]   = Stride;
        /*g=b+h*/ buf1[6] = pos_b;
        buf_1_stride[6]   = Stride;
        buf2[6]           = pos_h + Stride;
        buf_2_stride[6]   = Stride;
        /*e=h+b*/ buf1[7] = pos_h - 1 + Stride;
        buf_1_stride[7]   = Stride;
        buf2[7]           = pos_b;
        buf_2_stride[7]   = Stride;

        break;

    case EB_QUARTER_IN_HALF_VERTICAL:

        /*k=j+h*/ buf1[0] = pos_j;
        buf_1_stride[0]   = Stride;
        buf2[0]           = pos_h;
        buf_2_stride[0]   = Stride;
        /*i=h+j*/ buf1[1] = pos_h;
        buf_1_stride[1]   = Stride;
        buf2[1]           = pos_j + 1;
        buf_2_stride[1]   = Stride;
        /*d=A+h*/ buf1[2] = pos_Full - FullStride;
        buf_1_stride[2]   = FullStride;
        buf2[2]           = pos_h;
        buf_2_stride[2]   = Stride;
        /*n=h+A*/ buf1[3] = pos_h;
        buf_1_stride[3]   = Stride;
        buf2[3]           = pos_Full;
        buf_2_stride[3]   = FullStride;
        /*g=b+h*/ buf1[4] = pos_b - Stride;
        buf_1_stride[4]   = Stride;
        buf2[4]           = pos_h;
        buf_2_stride[4]   = Stride;
        /*e=h+b*/ buf1[5] = pos_h;
        buf_1_stride[5]   = Stride;
        buf2[5]           = pos_b + 1 - Stride;
        buf_2_stride[5]   = Stride;
        /*p=h+b*/ buf1[6] = pos_h;
        buf_1_stride[6]   = Stride;
        buf2[6]           = pos_b + 1;
        buf_2_stride[6]   = Stride;
        /*r=b+h*/ buf1[7] = pos_b;
        buf_1_stride[7]   = Stride;
        buf2[7]           = pos_h;
        buf_2_stride[7]   = Stride;

        break;

    case EB_QUARTER_IN_HALF_DIAGONAL:

        /*i=h+j*/ buf1[0] = pos_h - 1;
        buf_1_stride[0]   = Stride;
        buf2[0]           = pos_j;
        buf_2_stride[0]   = Stride;
        /*k=j+h*/ buf1[1] = pos_j;
        buf_1_stride[1]   = Stride;
        buf2[1]           = pos_h;
        buf_2_stride[1]   = Stride;
        /*f=b+j*/ buf1[2] = pos_b - Stride;
        buf_1_stride[2]   = Stride;
        buf2[2]           = pos_j;
        buf_2_stride[2]   = Stride;
        /*q=j+b*/ buf1[3] = pos_j;
        buf_1_stride[3]   = Stride;
        buf2[3]           = pos_b;
        buf_2_stride[3]   = Stride;
        /*e=h+b*/ buf1[4] = pos_h - 1;
        buf_1_stride[4]   = Stride;
        buf2[4]           = pos_b - Stride;
        buf_2_stride[4]   = Stride;
        /*g=b+h*/ buf1[5] = pos_b - Stride;
        buf_1_stride[5]   = Stride;
        buf2[5]           = pos_h;
        buf_2_stride[5]   = Stride;
        /*r=b+h*/ buf1[6] = pos_b;
        buf_1_stride[6]   = Stride;
        buf2[6]           = pos_h;
        buf_2_stride[6]   = Stride;
        /*p=h+b*/ buf1[7] = pos_h - 1;
        buf_1_stride[7]   = Stride;
        buf2[7]           = pos_b;
        buf_2_stride[7]   = Stride;

        break;

    default:
        SVT_ERROR("Unhandled quarter_pel_refinement_method: %" PRIu32 "\n",
            quarter_pel_refinement_method);
        assert(0);
        break;
    }

    return;
}

/*******************************************
 * quarter_pel_search_sb
 *   performs Quarter Pel refinement for the 85 PUs
 *******************************************/
static void quarter_pel_search_sb(
    MeContext *context_ptr, //[IN/OUT]  ME context Ptr, used to get/update ME results
    uint8_t *  pos_Full, //[IN]
    uint32_t   FullStride, //[IN]
    uint8_t *  pos_b, //[IN]
    uint8_t *  pos_h, //[IN]
    uint8_t *  pos_j, //[IN]
    int16_t    x_search_area_origin, //[IN] search area origin in the horizontal
    // direction, used to point to reference samples
    int16_t y_search_area_origin, //[IN] search area origin in the vertical
    // direction, used to point to reference samples
    EbBool enable_half_pel32x32, EbBool enable_half_pel16x16,
    EbBool enable_half_pel8x8, EbBool enable_quarter_pel, EbBool ext_block_flag) {
    uint32_t pu_index;

    uint32_t pu_shift_x_index;
    uint32_t pu_shift_y_index;

    uint32_t blk_sb_buffer_index;

    // for each one of the 8 positions, we need to determine the 2 buffers to do
    // averaging
    uint8_t *buf1[8];
    uint8_t *buf2[8];

    uint32_t buf_1_stride[8];
    uint32_t buf_2_stride[8];

    int16_t  x_mv, y_mv;
    uint32_t nidx;

    // 64x64 [1 partition]
    x_mv = _MVXT(*context_ptr->p_best_mv64x64);
    y_mv = _MVYT(*context_ptr->p_best_mv64x64);

    set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                    FullStride,
                                                    pos_b,
                                                    pos_h,
                                                    pos_j,
                                                    context_ptr->interpolated_stride,
                                                    x_mv,
                                                    y_mv,
                                                    buf1,
                                                    buf_1_stride,
                                                    buf2,
                                                    buf_2_stride);

    buf1[0] = buf1[0];
    buf2[0] = buf2[0];
    buf1[1] = buf1[1];
    buf2[1] = buf2[1];
    buf1[2] = buf1[2];
    buf2[2] = buf2[2];
    buf1[3] = buf1[3];
    buf2[3] = buf2[3];
    buf1[4] = buf1[4];
    buf2[4] = buf2[4];
    buf1[5] = buf1[5];
    buf2[5] = buf2[5];
    buf1[6] = buf1[6];
    buf2[6] = buf2[6];
    buf1[7] = buf1[7];
    buf2[7] = buf2[7];

    pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                            context_ptr->p_best_ssd64x64,
                                            0,
                                            buf1,
                                            buf_1_stride,
                                            buf2,
                                            buf_2_stride,
                                            64,
                                            64,
                                            x_search_area_origin,
                                            y_search_area_origin,
                                            context_ptr->p_best_sad_64x64,
                                            context_ptr->p_best_mv64x64,
                                            context_ptr->psub_pel_direction64x64);

    if (enable_quarter_pel && enable_half_pel32x32) {
        // 32x32 [4 partitions]
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            x_mv = _MVXT(context_ptr->p_best_mv32x32[pu_index]);
            y_mv = _MVYT(context_ptr->p_best_mv32x32[pu_index]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            pu_shift_x_index = (pu_index & 0x01) << 5;
            pu_shift_y_index = (pu_index >> 1) << 5;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd32x32[pu_index],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 32,
                                                 32,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_32x32[pu_index],
                                                 &context_ptr->p_best_mv32x32[pu_index],
                                                 context_ptr->psub_pel_direction32x32[pu_index]);
        }
    }

    if (enable_quarter_pel && enable_half_pel16x16) {
        // 16x16 [16 partitions]
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            nidx = tab16x16[pu_index];

            x_mv = _MVXT(context_ptr->p_best_mv16x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x16[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            pu_shift_x_index = (pu_index & 0x03) << 4;
            pu_shift_y_index = (pu_index >> 2) << 4;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd16x16[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 16,
                                                 16,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_16x16[nidx],
                                                 &context_ptr->p_best_mv16x16[nidx],
                                                 context_ptr->psub_pel_direction16x16[nidx]);
        }
    }

    if (enable_quarter_pel && enable_half_pel8x8) {
        // 8x8   [64 partitions]
        for (pu_index = 0; pu_index < 64; ++pu_index) {
            nidx = tab8x8[pu_index];

            x_mv = _MVXT(context_ptr->p_best_mv8x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x8[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                            FullStride,
                                                            pos_b,
                                                            pos_h,
                                                            pos_j,
                                                            context_ptr->interpolated_stride,
                                                            x_mv,
                                                            y_mv,
                                                            buf1,
                                                            buf_1_stride,
                                                            buf2,
                                                            buf_2_stride);

            pu_shift_x_index = (pu_index & 0x07) << 3;
            pu_shift_y_index = (pu_index >> 3) << 3;

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                    &context_ptr->p_best_ssd8x8[nidx],
                                                    blk_sb_buffer_index,
                                                    buf1,
                                                    buf_1_stride,
                                                    buf2,
                                                    buf_2_stride,
                                                    8,
                                                    8,
                                                    x_search_area_origin,
                                                    y_search_area_origin,
                                                    &context_ptr->p_best_sad_8x8[nidx],
                                                    &context_ptr->p_best_mv8x8[nidx],
                                                    context_ptr->psub_pel_direction8x8[nidx]);
        }
    }

    if (ext_block_flag) {
        // 64x32
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            pu_shift_x_index = 0;
            pu_shift_y_index = pu_index << 5;

            x_mv = _MVXT(context_ptr->p_best_mv64x32[pu_index]);
            y_mv = _MVYT(context_ptr->p_best_mv64x32[pu_index]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd64x32[pu_index],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 64,
                                                 32,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_64x32[pu_index],
                                                 &context_ptr->p_best_mv64x32[pu_index],
                                                 context_ptr->psub_pel_direction64x32[pu_index]);
        }

        // 32x16
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            nidx = tab32x16[pu_index]; // TODO bitwise this

            pu_shift_x_index = (pu_index & 0x01) << 5;
            pu_shift_y_index = (pu_index >> 1) << 4;

            x_mv = _MVXT(context_ptr->p_best_mv32x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv32x16[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd32x16[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 32,
                                                 16,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_32x16[nidx],
                                                 &context_ptr->p_best_mv32x16[nidx],
                                                 context_ptr->psub_pel_direction32x16[nidx]);
        }

        // 16x8
        for (pu_index = 0; pu_index < 32; ++pu_index) {
            nidx = tab16x8[pu_index];

            pu_shift_x_index = (pu_index & 0x03) << 4;
            pu_shift_y_index = (pu_index >> 2) << 3;

            x_mv = _MVXT(context_ptr->p_best_mv16x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x8[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd16x8[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 16,
                                                 8,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_16x8[nidx],
                                                 &context_ptr->p_best_mv16x8[nidx],
                                                 context_ptr->psub_pel_direction16x8[nidx]);
        }

        // 32x64
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            pu_shift_x_index = pu_index << 5;
            pu_shift_y_index = 0;

            x_mv = _MVXT(context_ptr->p_best_mv32x64[pu_index]);
            y_mv = _MVYT(context_ptr->p_best_mv32x64[pu_index]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd32x64[pu_index],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 32,
                                                 64,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_32x64[pu_index],
                                                 &context_ptr->p_best_mv32x64[pu_index],
                                                 context_ptr->psub_pel_direction32x64[pu_index]);
        }

        // 16x32
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            nidx = tab16x32[pu_index];

            pu_shift_x_index = (pu_index & 0x03) << 4;
            pu_shift_y_index = (pu_index >> 2) << 5;

            x_mv = _MVXT(context_ptr->p_best_mv16x32[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x32[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd16x32[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 16,
                                                 32,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_16x32[nidx],
                                                 &context_ptr->p_best_mv16x32[nidx],
                                                 context_ptr->psub_pel_direction16x32[nidx]);
        }

        // 8x16
        for (pu_index = 0; pu_index < 32; ++pu_index) {
            nidx = tab8x16[pu_index];

            pu_shift_x_index = (pu_index & 0x07) << 3;
            pu_shift_y_index = (pu_index >> 3) << 4;

            x_mv = _MVXT(context_ptr->p_best_mv8x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x16[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd8x16[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 8,
                                                 16,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_8x16[nidx],
                                                 &context_ptr->p_best_mv8x16[nidx],
                                                 context_ptr->psub_pel_direction8x16[nidx]);
        }

        // 32x8
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            nidx = tab32x8[pu_index];

            pu_shift_x_index = (pu_index & 0x01) << 5;
            pu_shift_y_index = (pu_index >> 1) << 3;

            x_mv = _MVXT(context_ptr->p_best_mv32x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv32x8[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd32x8[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 32,
                                                 8,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_32x8[nidx],
                                                 &context_ptr->p_best_mv32x8[nidx],
                                                 context_ptr->psub_pel_direction32x8[nidx]);
        }

        // 8x32
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            nidx = tab8x32[pu_index];

            pu_shift_x_index = (pu_index & 0x07) << 3;
            pu_shift_y_index = (pu_index >> 3) << 5;

            x_mv = _MVXT(context_ptr->p_best_mv8x32[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x32[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd8x32[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 8,
                                                 32,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_8x32[nidx],
                                                 &context_ptr->p_best_mv8x32[nidx],
                                                 context_ptr->psub_pel_direction8x32[nidx]);
        }

        // 64x16
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            nidx = pu_index;

            pu_shift_x_index = 0;
            pu_shift_y_index = pu_index << 4;

            x_mv = _MVXT(context_ptr->p_best_mv64x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv64x16[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd64x16[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 64,
                                                 16,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_64x16[nidx],
                                                 &context_ptr->p_best_mv64x16[nidx],
                                                 context_ptr->psub_pel_direction64x16[nidx]);
        }

        // 16x64
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            nidx = pu_index;

            pu_shift_x_index = pu_index << 4;
            pu_shift_y_index = 0;

            x_mv = _MVXT(context_ptr->p_best_mv16x64[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x64[nidx]);

            set_quarter_pel_refinement_inputs_on_the_fly(pos_Full,
                                                         FullStride,
                                                         pos_b,
                                                         pos_h,
                                                         pos_j,
                                                         context_ptr->interpolated_stride,
                                                         x_mv,
                                                         y_mv,
                                                         buf1,
                                                         buf_1_stride,
                                                         buf2,
                                                         buf_2_stride);

            blk_sb_buffer_index = pu_shift_x_index + pu_shift_y_index * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[0];
            buf2[0] = buf2[0] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[0];
            buf1[1] = buf1[1] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[1];
            buf2[1] = buf2[1] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[1];
            buf1[2] = buf1[2] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[2];
            buf2[2] = buf2[2] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[2];
            buf1[3] = buf1[3] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[3];
            buf2[3] = buf2[3] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[3];
            buf1[4] = buf1[4] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[4];
            buf2[4] = buf2[4] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[4];
            buf1[5] = buf1[5] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[5];
            buf2[5] = buf2[5] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[5];
            buf1[6] = buf1[6] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[6];
            buf2[6] = buf2[6] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[6];
            buf1[7] = buf1[7] + pu_shift_x_index + pu_shift_y_index * buf_1_stride[7];
            buf2[7] = buf2[7] + pu_shift_x_index + pu_shift_y_index * buf_2_stride[7];

            pu_quarter_pel_refinement_on_the_fly(context_ptr,
                                                 &context_ptr->p_best_ssd16x64[nidx],
                                                 blk_sb_buffer_index,
                                                 buf1,
                                                 buf_1_stride,
                                                 buf2,
                                                 buf_2_stride,
                                                 16,
                                                 64,
                                                 x_search_area_origin,
                                                 y_search_area_origin,
                                                 &context_ptr->p_best_sad_16x64[nidx],
                                                 &context_ptr->p_best_mv16x64[nidx],
                                                 context_ptr->psub_pel_direction16x64[nidx]);
        }
    }

    return;
}
void hme_one_quadrant_level_0(
    PictureParentControlSet *pcs_ptr,
    MeContext *              context_ptr, // input/output parameter, ME context Ptr, used to
    // get/update ME results
    int16_t origin_x, // input parameter, SB position in the horizontal
    // direction- sixteenth resolution
    int16_t origin_y, // input parameter, SB position in the vertical
    // direction- sixteenth resolution
    uint32_t sb_width, // input parameter, SB pwidth - sixteenth resolution
    uint32_t sb_height, // input parameter, SB height - sixteenth resolution
    int16_t  x_hme_search_center, // input parameter, HME search center in the
    // horizontal direction
    int16_t y_hme_search_center, // input parameter, HME search center in the
    // vertical direction
    EbPictureBufferDesc *sixteenth_ref_pic_ptr, // input parameter, sixteenth reference Picture Ptr
    uint64_t *           level0Bestsad_, // output parameter, Level0 SAD at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    int16_t *xLevel0SearchCenter, // output parameter, Level0 xMV at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    int16_t *yLevel0SearchCenter, // output parameter, Level0 yMV at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    uint32_t searchAreaMultiplierX, uint32_t searchAreaMultiplierY) {
    int16_t  x_top_left_search_region;
    int16_t  y_top_left_search_region;
    uint32_t search_region_index;
    int16_t  x_search_area_origin;
    int16_t  y_search_area_origin;
    int16_t  x_search_region_distance;
    int16_t  y_search_region_distance;

    int16_t pad_width;
    int16_t pad_height;

    (void)pcs_ptr;
    // Round up x_HME_L0 to be a multiple of 16
    int16_t search_area_width = (int16_t)(
        (((((context_ptr->hme_level0_total_search_area_width * searchAreaMultiplierX) / 100))) +
         15) &
        ~0x0F);
    int16_t search_area_height = (int16_t)(
        ((context_ptr->hme_level0_total_search_area_height * searchAreaMultiplierY) / 100));
    x_search_region_distance = x_hme_search_center;
    y_search_region_distance = y_hme_search_center;
    pad_width                = (int16_t)(sixteenth_ref_pic_ptr->origin_x) - 1;
    pad_height               = (int16_t)(sixteenth_ref_pic_ptr->origin_y) - 1;

    x_search_area_origin = -(int16_t)(search_area_width >> 1) + x_search_region_distance;
    y_search_area_origin = -(int16_t)(search_area_height >> 1) + y_search_region_distance;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -pad_width) ? -pad_width - origin_x
                                                                            : x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -pad_width)
                            ? search_area_width - (-pad_width - (origin_x + x_search_area_origin))
                            : search_area_width;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin =
        ((origin_x + x_search_area_origin) > (int16_t)sixteenth_ref_pic_ptr->width - 1)
            ? x_search_area_origin -
                  ((origin_x + x_search_area_origin) - ((int16_t)sixteenth_ref_pic_ptr->width - 1))
            : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) >
         (int16_t)sixteenth_ref_pic_ptr->width)
            ? MAX(1,
                  search_area_width - ((origin_x + x_search_area_origin + search_area_width) -
                                       (int16_t)sixteenth_ref_pic_ptr->width))
            : search_area_width;

    // Round down x_HME to be a multiple of 16 as cropping already performed
    search_area_width = (search_area_width < 16) ? search_area_width : search_area_width & ~0x0F;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -pad_height)
                               ? -pad_height - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -pad_height)
            ? search_area_height - (-pad_height - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin =
        ((origin_y + y_search_area_origin) > (int16_t)sixteenth_ref_pic_ptr->height - 1)
            ? y_search_area_origin -
                  ((origin_y + y_search_area_origin) - ((int16_t)sixteenth_ref_pic_ptr->height - 1))
            : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height >
         (int16_t)sixteenth_ref_pic_ptr->height)
            ? MAX(1,
                  search_area_height - ((origin_y + y_search_area_origin + search_area_height) -
                                        (int16_t)sixteenth_ref_pic_ptr->height))
            : search_area_height;

    x_top_left_search_region =
        ((int16_t)sixteenth_ref_pic_ptr->origin_x + origin_x) + x_search_area_origin;
    y_top_left_search_region =
        ((int16_t)sixteenth_ref_pic_ptr->origin_y + origin_y) + y_search_area_origin;
    search_region_index =
        x_top_left_search_region + y_top_left_search_region * sixteenth_ref_pic_ptr->stride_y;

    if (context_ptr->hme_search_type == HME_SPARSE) {
        sad_loop_kernel_sparse(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride,
            &sixteenth_ref_pic_ptr->buffer_y[search_region_index],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? sixteenth_ref_pic_ptr->stride_y
                : sixteenth_ref_pic_ptr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
            sb_width,
            /* results */
            level0Bestsad_,
            xLevel0SearchCenter,
            yLevel0SearchCenter,
            /* range */
            sixteenth_ref_pic_ptr->stride_y,
            search_area_width,
            search_area_height);
    } else {
        if ((search_area_width & 15) == 0) {
            // Only width equals 16 (SB equals 64) is updated
            // other width sizes work with the old code as the one
            // in"sad_loop_kernel_sse4_1_intrin"
            sad_loop_kernel_hme_l0(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenth_ref_pic_ptr->buffer_y[search_region_index],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenth_ref_pic_ptr->stride_y
                    : sixteenth_ref_pic_ptr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
                sb_width,
                /* results */
                level0Bestsad_,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenth_ref_pic_ptr->stride_y,
                search_area_width,
                search_area_height);
        } else {
            // Put the first search location into level0 results
            sad_loop_kernel(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenth_ref_pic_ptr->buffer_y[search_region_index],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenth_ref_pic_ptr->stride_y
                    : sixteenth_ref_pic_ptr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
                sb_width,
                /* results */
                level0Bestsad_,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenth_ref_pic_ptr->stride_y,
                search_area_width,
                search_area_height);
        }
    }

    *level0Bestsad_ =
        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
            ? *level0Bestsad_
            : *level0Bestsad_ * 2; // Multiply by 2 because considered only ever other line
    *xLevel0SearchCenter += x_search_area_origin;
    *xLevel0SearchCenter *= 4; // Multiply by 4 because operating on 1/4 resolution
    *yLevel0SearchCenter += y_search_area_origin;
    *yLevel0SearchCenter *= 4; // Multiply by 4 because operating on 1/4 resolution

    return;
}

void hme_level_0(
    PictureParentControlSet *pcs_ptr,
    MeContext *              context_ptr, // input/output parameter, ME context Ptr, used to
    // get/update ME results
    int16_t origin_x, // input parameter, SB position in the horizontal
    // direction- sixteenth resolution
    int16_t origin_y, // input parameter, SB position in the vertical
    // direction- sixteenth resolution
    uint32_t sb_width, // input parameter, SB pwidth - sixteenth resolution
    uint32_t sb_height, // input parameter, SB height - sixteenth resolution
    int16_t  x_hme_search_center, // input parameter, HME search center in the
    // horizontal direction
    int16_t y_hme_search_center, // input parameter, HME search center in the
    // vertical direction
    EbPictureBufferDesc *sixteenth_ref_pic_ptr, // input parameter, sixteenth reference Picture Ptr
    uint32_t             search_region_number_in_width, // input parameter, search region
    // number in the horizontal direction
    uint32_t search_region_number_in_height, // input parameter, search region
    // number in the vertical direction
    uint64_t *level0Bestsad_, // output parameter, Level0 SAD at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    int16_t *xLevel0SearchCenter, // output parameter, Level0 xMV at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    int16_t *yLevel0SearchCenter, // output parameter, Level0 yMV at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    uint32_t searchAreaMultiplierX, uint32_t searchAreaMultiplierY) {
    int16_t  x_top_left_search_region;
    int16_t  y_top_left_search_region;
    uint32_t search_region_index;
    int16_t  x_search_area_origin;
    int16_t  y_search_area_origin;
    int16_t  x_search_region_distance;
    int16_t  y_search_region_distance;

    int16_t pad_width;
    int16_t pad_height;

    // Adjust SR size based on the searchAreaShift
    (void)pcs_ptr;
    // Round up x_HME_L0 to be a multiple of 16
    int16_t search_area_width = (int16_t)(
        (((((context_ptr->hme_level0_search_area_in_width_array[search_region_number_in_width] *
             searchAreaMultiplierX) /
            100))) +
         15) &
        ~0x0F);
    int16_t search_area_height = (int16_t)(
        ((context_ptr->hme_level0_search_area_in_height_array[search_region_number_in_height] *
          searchAreaMultiplierY) /
         100));

    x_search_region_distance = x_hme_search_center;
    y_search_region_distance = y_hme_search_center;
    pad_width                = (int16_t)(sixteenth_ref_pic_ptr->origin_x) - 1;
    pad_height               = (int16_t)(sixteenth_ref_pic_ptr->origin_y) - 1;

    while (search_region_number_in_width) {
        search_region_number_in_width--;
        x_search_region_distance += (int16_t)(
            ((context_ptr->hme_level0_search_area_in_width_array[search_region_number_in_width] *
              searchAreaMultiplierX) /
             100));
    }

    while (search_region_number_in_height) {
        search_region_number_in_height--;
        y_search_region_distance += (int16_t)(
            ((context_ptr->hme_level0_search_area_in_height_array[search_region_number_in_height] *
              searchAreaMultiplierY) /
             100));
    }
    x_search_area_origin =
        -(int16_t)(
            (((context_ptr->hme_level0_total_search_area_width * searchAreaMultiplierX) / 100)) >>
            1) +
        x_search_region_distance;
    y_search_area_origin =
        -(int16_t)(
            (((context_ptr->hme_level0_total_search_area_height * searchAreaMultiplierY) / 100)) >>
            1) +
        y_search_region_distance;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -pad_width) ? -pad_width - origin_x
                                                                            : x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -pad_width)
                            ? search_area_width - (-pad_width - (origin_x + x_search_area_origin))
                            : search_area_width;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin =
        ((origin_x + x_search_area_origin) > (int16_t)sixteenth_ref_pic_ptr->width - 1)
            ? x_search_area_origin -
                  ((origin_x + x_search_area_origin) - ((int16_t)sixteenth_ref_pic_ptr->width - 1))
            : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) >
         (int16_t)sixteenth_ref_pic_ptr->width)
            ? MAX(1,
                  search_area_width - ((origin_x + x_search_area_origin + search_area_width) -
                                       (int16_t)sixteenth_ref_pic_ptr->width))
            : search_area_width;

    // Round down x_HME to be a multiple of 16 as cropping already performed
    search_area_width = (search_area_width < 16) ? search_area_width : search_area_width & ~0x0F;

    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -pad_height)
                               ? -pad_height - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -pad_height)
            ? search_area_height - (-pad_height - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin =
        ((origin_y + y_search_area_origin) > (int16_t)sixteenth_ref_pic_ptr->height - 1)
            ? y_search_area_origin -
                  ((origin_y + y_search_area_origin) - ((int16_t)sixteenth_ref_pic_ptr->height - 1))
            : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height >
         (int16_t)sixteenth_ref_pic_ptr->height)
            ? MAX(1,
                  search_area_height - ((origin_y + y_search_area_origin + search_area_height) -
                                        (int16_t)sixteenth_ref_pic_ptr->height))
            : search_area_height;

    x_top_left_search_region =
        ((int16_t)sixteenth_ref_pic_ptr->origin_x + origin_x) + x_search_area_origin;
    y_top_left_search_region =
        ((int16_t)sixteenth_ref_pic_ptr->origin_y + origin_y) + y_search_area_origin;
    search_region_index =
        x_top_left_search_region + y_top_left_search_region * sixteenth_ref_pic_ptr->stride_y;

    if (((sb_width & 7) == 0) || (sb_width == 4)) {
        if ((search_area_width & 15) == 0) {
            // Only width equals 16 (SB equals 64) is updated
            // other width sizes work with the old code as the one
            // in"sad_loop_kernel_sse4_1_intrin"
            sad_loop_kernel_hme_l0(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenth_ref_pic_ptr->buffer_y[search_region_index],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenth_ref_pic_ptr->stride_y
                    : sixteenth_ref_pic_ptr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
                sb_width,
                /* results */
                level0Bestsad_,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenth_ref_pic_ptr->stride_y,
                search_area_width,
                search_area_height);
        } else {
            // Put the first search location into level0 results
            sad_loop_kernel(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenth_ref_pic_ptr->buffer_y[search_region_index],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenth_ref_pic_ptr->stride_y
                    : sixteenth_ref_pic_ptr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
                sb_width,
                /* results */
                level0Bestsad_,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenth_ref_pic_ptr->stride_y,
                search_area_width,
                search_area_height);
        }
    } else {
        sad_loop_kernel_c(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride,
            &sixteenth_ref_pic_ptr->buffer_y[search_region_index],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? sixteenth_ref_pic_ptr->stride_y
                : sixteenth_ref_pic_ptr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
            sb_width,
            /* results */
            level0Bestsad_,
            xLevel0SearchCenter,
            yLevel0SearchCenter,
            /* range */
            sixteenth_ref_pic_ptr->stride_y,
            search_area_width,
            search_area_height);
    }

    *level0Bestsad_ =
        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
            ? *level0Bestsad_
            : *level0Bestsad_ * 2; // Multiply by 2 because considered only ever other line
    *xLevel0SearchCenter += x_search_area_origin;
    *xLevel0SearchCenter *= 4; // Multiply by 4 because operating on 1/4 resolution
    *yLevel0SearchCenter += y_search_area_origin;
    *yLevel0SearchCenter *= 4; // Multiply by 4 because operating on 1/4 resolution

    return;
}

void hme_level_1(
    MeContext *context_ptr, // input/output parameter, ME context Ptr, used to
    // get/update ME results
    int16_t origin_x, // input parameter, SB position in the horizontal
    // direction - quarter resolution
    int16_t origin_y, // input parameter, SB position in the vertical direction
    // - quarter resolution
    uint32_t             sb_width, // input parameter, SB pwidth - quarter resolution
    uint32_t             sb_height, // input parameter, SB height - quarter resolution
    EbPictureBufferDesc *quarter_ref_pic_ptr, // input parameter, quarter reference Picture Ptr
    int16_t              hme_level1_search_area_in_width, // input parameter, hme level 1 search
    // area in width
    int16_t hme_level1_search_area_in_height, // input parameter, hme level 1 search
    // area in height
    int16_t xLevel0SearchCenter, // input parameter, best Level0 xMV at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    int16_t yLevel0SearchCenter, // input parameter, best Level0 yMV at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    uint64_t *level1Bestsad_, // output parameter, Level1 SAD at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    int16_t *xLevel1SearchCenter, // output parameter, Level1 xMV at
    // (search_region_number_in_width,
    // search_region_number_in_height)
    int16_t *yLevel1SearchCenter // output parameter, Level1 yMV at
    // (search_region_number_in_width,
    // search_region_number_in_height)
) {
    int16_t  x_top_left_search_region;
    int16_t  y_top_left_search_region;
    uint32_t search_region_index;
    // Round up x_HME_L0 to be a multiple of 8
    int16_t search_area_width  = (int16_t)((hme_level1_search_area_in_width + 7) & ~0x07);
    int16_t search_area_height = hme_level1_search_area_in_height;

    int16_t x_search_area_origin;
    int16_t y_search_area_origin;

    int16_t pad_width  = (int16_t)(quarter_ref_pic_ptr->origin_x) - 1;
    int16_t pad_height = (int16_t)(quarter_ref_pic_ptr->origin_y) - 1;

    x_search_area_origin = -(search_area_width >> 1) + xLevel0SearchCenter;
    y_search_area_origin = -(search_area_height >> 1) + yLevel0SearchCenter;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -pad_width) ? -pad_width - origin_x
                                                                            : x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -pad_width)
                            ? search_area_width - (-pad_width - (origin_x + x_search_area_origin))
                            : search_area_width;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin =
        ((origin_x + x_search_area_origin) > (int16_t)quarter_ref_pic_ptr->width - 1)
            ? x_search_area_origin -
                  ((origin_x + x_search_area_origin) - ((int16_t)quarter_ref_pic_ptr->width - 1))
            : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) >
         (int16_t)quarter_ref_pic_ptr->width)
            ? MAX(1,
                  search_area_width - ((origin_x + x_search_area_origin + search_area_width) -
                                       (int16_t)quarter_ref_pic_ptr->width))
            : search_area_width;

    // Constrain x_HME_L1 to be a multiple of 8 (round down as cropping already
    // performed)
    search_area_width = (search_area_width < 8) ? search_area_width : search_area_width & ~0x07;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -pad_height)
                               ? -pad_height - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -pad_height)
            ? search_area_height - (-pad_height - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin =
        ((origin_y + y_search_area_origin) > (int16_t)quarter_ref_pic_ptr->height - 1)
            ? y_search_area_origin -
                  ((origin_y + y_search_area_origin) - ((int16_t)quarter_ref_pic_ptr->height - 1))
            : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height >
         (int16_t)quarter_ref_pic_ptr->height)
            ? MAX(1,
                  search_area_height - ((origin_y + y_search_area_origin + search_area_height) -
                                        (int16_t)quarter_ref_pic_ptr->height))
            : search_area_height;

    // Move to the top left of the search region
    x_top_left_search_region =
        ((int16_t)quarter_ref_pic_ptr->origin_x + origin_x) + x_search_area_origin;
    y_top_left_search_region =
        ((int16_t)quarter_ref_pic_ptr->origin_y + origin_y) + y_search_area_origin;
    search_region_index =
        x_top_left_search_region + y_top_left_search_region * quarter_ref_pic_ptr->stride_y;

    if (((sb_width & 7) == 0) || (sb_width == 4)) {
        // Put the first search location into level0 results
        sad_loop_kernel(
            &context_ptr->quarter_sb_buffer[0],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? context_ptr->quarter_sb_buffer_stride
                : context_ptr->quarter_sb_buffer_stride * 2,
            &quarter_ref_pic_ptr->buffer_y[search_region_index],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? quarter_ref_pic_ptr->stride_y
                                                                : quarter_ref_pic_ptr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
            sb_width,
            /* results */
            level1Bestsad_,
            xLevel1SearchCenter,
            yLevel1SearchCenter,
            /* range */
            quarter_ref_pic_ptr->stride_y,
            search_area_width,
            search_area_height);
    } else {
        sad_loop_kernel_c(
            &context_ptr->quarter_sb_buffer[0],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? context_ptr->quarter_sb_buffer_stride
                : context_ptr->quarter_sb_buffer_stride * 2,
            &quarter_ref_pic_ptr->buffer_y[search_region_index],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? quarter_ref_pic_ptr->stride_y
                                                                : quarter_ref_pic_ptr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
            sb_width,
            /* results */
            level1Bestsad_,
            xLevel1SearchCenter,
            yLevel1SearchCenter,
            /* range */
            quarter_ref_pic_ptr->stride_y,
            search_area_width,
            search_area_height);
    }

    *level1Bestsad_ =
        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
            ? *level1Bestsad_
            : *level1Bestsad_ * 2; // Multiply by 2 because considered only ever other line
    *xLevel1SearchCenter += x_search_area_origin;
    *xLevel1SearchCenter *= 2; // Multiply by 2 because operating on 1/2 resolution
    *yLevel1SearchCenter += y_search_area_origin;
    *yLevel1SearchCenter *= 2; // Multiply by 2 because operating on 1/2 resolution

    return;
}

void hme_level_2(PictureParentControlSet *pcs_ptr, // input parameter, Picture control set Ptr
                 MeContext *context_ptr, // input/output parameter, ME context Ptr, used to
                 // get/update ME results
                 int16_t  origin_x, // input parameter, SB position in the horizontal direction
                 int16_t  origin_y, // input parameter, SB position in the vertical direction
                 uint32_t sb_width, // input parameter, SB pwidth - full resolution
                 uint32_t sb_height, // input parameter, SB height - full resolution
                 EbPictureBufferDesc *ref_pic_ptr, // input parameter, reference Picture Ptr
                 uint32_t search_region_number_in_width, // input parameter, search region
                 // number in the horizontal direction
                 uint32_t search_region_number_in_height, // input parameter, search region
                 // number in the vertical direction
                 int16_t xLevel1SearchCenter, // input parameter, best Level1 xMV
                 // at(search_region_number_in_width,
                 // search_region_number_in_height)
                 int16_t yLevel1SearchCenter, // input parameter, best Level1 yMV
                 // at(search_region_number_in_width,
                 // search_region_number_in_height)
                 uint64_t *level2Bestsad_, // output parameter, Level2 SAD at
                 // (search_region_number_in_width,
                 // search_region_number_in_height)
                 int16_t *xLevel2SearchCenter, // output parameter, Level2 xMV at
                 // (search_region_number_in_width,
                 // search_region_number_in_height)
                 int16_t *yLevel2SearchCenter // output parameter, Level2 yMV at
                 // (search_region_number_in_width,
                 // search_region_number_in_height)
) {
    int16_t  x_top_left_search_region;
    int16_t  y_top_left_search_region;
    uint32_t search_region_index;

    // round the search region width to nearest multiple of 8 if it is less than
    // 8 or non multiple of 8 SAD calculation performance is the same for
    // searchregion width from 1 to 8
    (void)pcs_ptr;
    int16_t hme_level2_search_area_in_width =
        (int16_t)context_ptr->hme_level2_search_area_in_width_array[search_region_number_in_width];
    // Round up x_HME_L0 to be a multiple of 8
    int16_t search_area_width = (int16_t)((hme_level2_search_area_in_width + 7) & ~0x07);
    int16_t search_area_height =
        (int16_t)
            context_ptr->hme_level2_search_area_in_height_array[search_region_number_in_height];
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;

    int16_t pad_width  = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t pad_height = (int16_t)BLOCK_SIZE_64 - 1;

    x_search_area_origin = -(search_area_width >> 1) + xLevel1SearchCenter;
    y_search_area_origin = -(search_area_height >> 1) + yLevel1SearchCenter;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -pad_width) ? -pad_width - origin_x
                                                                            : x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -pad_width)
                            ? search_area_width - (-pad_width - (origin_x + x_search_area_origin))
                            : search_area_width;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) > (int16_t)ref_pic_ptr->width - 1)
                               ? x_search_area_origin - ((origin_x + x_search_area_origin) -
                                                         ((int16_t)ref_pic_ptr->width - 1))
                               : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) > (int16_t)ref_pic_ptr->width)
            ? MAX(1,
                  search_area_width - ((origin_x + x_search_area_origin + search_area_width) -
                                       (int16_t)ref_pic_ptr->width))
            : search_area_width;

    // Constrain x_HME_L1 to be a multiple of 8 (round down as cropping already
    // performed)
    search_area_width = (search_area_width < 8) ? search_area_width : search_area_width & ~0x07;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -pad_height)
                               ? -pad_height - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -pad_height)
            ? search_area_height - (-pad_height - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) > (int16_t)ref_pic_ptr->height - 1)
                               ? y_search_area_origin - ((origin_y + y_search_area_origin) -
                                                         ((int16_t)ref_pic_ptr->height - 1))
                               : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height > (int16_t)ref_pic_ptr->height)
            ? MAX(1,
                  search_area_height - ((origin_y + y_search_area_origin + search_area_height) -
                                        (int16_t)ref_pic_ptr->height))
            : search_area_height;

    // Move to the top left of the search region
    x_top_left_search_region = ((int16_t)ref_pic_ptr->origin_x + origin_x) + x_search_area_origin;
    y_top_left_search_region = ((int16_t)ref_pic_ptr->origin_y + origin_y) + y_search_area_origin;
    search_region_index =
        x_top_left_search_region + y_top_left_search_region * ref_pic_ptr->stride_y;
    if ((((sb_width & 7) == 0) && (sb_width != 40) && (sb_width != 56))) {
        // Put the first search location into level0 results
        sad_loop_kernel(
            context_ptr->sb_src_ptr,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? context_ptr->sb_src_stride
                                                                : context_ptr->sb_src_stride * 2,
            &ref_pic_ptr->buffer_y[search_region_index],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? ref_pic_ptr->stride_y
                                                                : ref_pic_ptr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
            sb_width,
            /* results */
            level2Bestsad_,
            xLevel2SearchCenter,
            yLevel2SearchCenter,
            /* range */
            ref_pic_ptr->stride_y,
            search_area_width,
            search_area_height);
    } else {
        // Put the first search location into level0 results
        sad_loop_kernel_c(
            context_ptr->sb_src_ptr,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? context_ptr->sb_src_stride
                                                                : context_ptr->sb_src_stride * 2,
            &ref_pic_ptr->buffer_y[search_region_index],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? ref_pic_ptr->stride_y
                                                                : ref_pic_ptr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH) ? sb_height : sb_height >> 1,
            sb_width,
            /* results */
            level2Bestsad_,
            xLevel2SearchCenter,
            yLevel2SearchCenter,
            /* range */
            ref_pic_ptr->stride_y,
            search_area_width,
            search_area_height);
    }

    *level2Bestsad_ =
        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
            ? *level2Bestsad_
            : *level2Bestsad_ * 2; // Multiply by 2 because considered only ever other line
    *xLevel2SearchCenter += x_search_area_origin;
    *yLevel2SearchCenter += y_search_area_origin;

    return;
}

static void select_buffer_luma(uint32_t  pu_index, //[IN]
                               uint8_t   fracPosition, //[IN]
                               uint32_t  pu_width, //[IN] Refrence picture list index
                               uint32_t  pu_height, //[IN] Refrence picture index in the list
                               uint8_t * pos_Full, //[IN]
                               uint8_t * pos_b, //[IN]
                               uint8_t * pos_h, //[IN]
                               uint8_t * pos_j, //[IN]
                               uint32_t  refHalfStride, //[IN]
                               uint32_t  refBufferFullStride,
                               uint8_t **dst_ptr, //[OUT]
                               uint32_t *DstPtrStride) //[OUT]
{
    (void)pu_width;
    (void)pu_height;

    uint32_t pu_shift_x_index = pu_search_index_map[pu_index][0];
    uint32_t pu_shift_y_index = pu_search_index_map[pu_index][1];
    uint32_t ref_stride       = refHalfStride;

    // for each one of the 8 positions, we need to determine the 2 buffers to do
    // averaging
    uint8_t *buf1 = pos_Full;

    switch (fracPosition) {
    case 0: // integer
        buf1       = pos_Full;
        ref_stride = refBufferFullStride;
        break;
    case 2: // b
        buf1 = pos_b;
        break;
    case 8: // h
        buf1 = pos_h;
        break;
    case 10: // j
        buf1 = pos_j;
        break;
    default: break;
    }

    buf1 = buf1 + pu_shift_x_index + pu_shift_y_index * ref_stride;

    *dst_ptr      = buf1;
    *DstPtrStride = ref_stride;

    return;
}

static void quarder_pel_compensation(uint32_t pu_index, //[IN]
                                     uint8_t  fracPosition, //[IN]
                                     uint32_t pu_width, //[IN] Refrence picture list index
                                     uint32_t pu_height, //[IN] Refrence picture index in the list
                                     uint8_t *pos_Full, //[IN]
                                     uint8_t *pos_b, //[IN]
                                     uint8_t *pos_h, //[IN]
                                     uint8_t *pos_j, //[IN]
                                     uint32_t refHalfStride, //[IN]
                                     uint32_t refBufferFullStride,
                                     uint8_t *Dst, //[IN]
                                     uint32_t DstStride) { //[IN]

    uint32_t pu_shift_x_index = pu_search_index_map[pu_index][0];
    uint32_t pu_shift_y_index = pu_search_index_map[pu_index][1];
    uint32_t ref_stride_1     = refHalfStride;
    uint32_t ref_stride_2     = refHalfStride;

    // for each one of the 8 positions, we need to determine the 2 buffers to do
    // averaging
    uint8_t *buf1 = pos_Full;
    uint8_t *buf2 = pos_Full;

    switch (fracPosition) {
    case 1: // a
        buf1         = pos_Full;
        buf2         = pos_b;
        ref_stride_1 = refBufferFullStride;
        break;

    case 3: // c
        buf1         = pos_b;
        buf2         = pos_Full + 1;
        ref_stride_2 = refBufferFullStride;
        break;

    case 4: // d
        buf1         = pos_Full;
        buf2         = pos_h;
        ref_stride_1 = refBufferFullStride;
        break;

    case 5: // e
        buf1 = pos_b;
        buf2 = pos_h;
        break;

    case 6: // f
        buf1 = pos_b;
        buf2 = pos_j;
        break;

    case 7: // g
        buf1 = pos_b;
        buf2 = pos_h + 1;
        break;

    case 9: // i
        buf1 = pos_h;
        buf2 = pos_j;
        break;

    case 11: // k
        buf1 = pos_j;
        buf2 = pos_h + 1;
        break;

    case 12: // L
        buf1         = pos_h;
        buf2         = pos_Full + refBufferFullStride;
        ref_stride_2 = refBufferFullStride;
        break;

    case 13: // m
        buf1 = pos_h;
        buf2 = pos_b + refHalfStride;
        break;

    case 14: // n
        buf1 = pos_j;
        buf2 = pos_b + refHalfStride;
        break;
    case 15: // 0
        buf1 = pos_h + 1;
        buf2 = pos_b + refHalfStride;
        break;
    default: break;
    }

    buf1 = buf1 + pu_shift_x_index + pu_shift_y_index * ref_stride_1;
    buf2 = buf2 + pu_shift_x_index + pu_shift_y_index * ref_stride_2;

    picture_average_kernel(
        buf1, ref_stride_1, buf2, ref_stride_2, Dst, DstStride, pu_width, pu_height);

    return;
}
/*******************************************************************************
 * Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
 * Requirement: pu_height % 2 = 0
 * Requirement: skip         = 0 or 1
 * Requirement (x86 only): temp_buf % 16 = 0
 * Requirement (x86 only): (dst->buffer_y  + dst_luma_index  ) % 16 = 0 when
 *pu_width %16 = 0 Requirement (x86 only): (dst->buffer_cb + dst_chroma_index) %
 *16 = 0 when pu_width %32 = 0 Requirement (x86 only): (dst->buffer_cr +
 *dst_chroma_index) % 16 = 0 when pu_width %32 = 0 Requirement (x86 only):
 *dst->stride_y   % 16 = 0 when pu_width %16 = 0 Requirement (x86 only):
 *dst->chroma_stride % 16 = 0 when pu_width %32 = 0
 *******************************************************************************/
uint32_t bi_pred_averaging(MeContext *context_ptr, MePredUnit *me_candidate, uint32_t pu_index,
                           uint8_t *sourcePic, uint32_t luma_stride, uint8_t firstFracPos,
                           uint8_t secondFracPos, uint32_t pu_width, uint32_t pu_height,
                           uint8_t *firstRefInteger, uint8_t *firstRefPosB, uint8_t *firstRefPosh_,
                           uint8_t *firstRefPosJ, uint8_t *secondRefInteger, uint8_t *secondRefPosB,
                           uint8_t *secondRefPosh_, uint8_t *secondRefPosJ,
                           uint32_t refBufferStride, uint32_t refBufferFullList0Stride,
                           uint32_t refBufferFullList1Stride, uint8_t *firstRefTempDst,
                           uint8_t *secondRefTempDst) {
    uint8_t *ptr_list0, *ptr_list1;
    uint32_t ptr_list0_stride, ptr_list1_stride;

    // Buffer Selection and quater-pel compensation on the fly
    if (sub_position_type[firstFracPos] != 2) {
        select_buffer_luma(pu_index,
                           firstFracPos,
                           pu_width,
                           pu_height,
                           firstRefInteger,
                           firstRefPosB,
                           firstRefPosh_,
                           firstRefPosJ,
                           refBufferStride,
                           refBufferFullList0Stride,
                           &ptr_list0,
                           &ptr_list0_stride);
    } else {
        quarder_pel_compensation(pu_index,
                                 firstFracPos,
                                 pu_width,
                                 pu_height,
                                 firstRefInteger,
                                 firstRefPosB,
                                 firstRefPosh_,
                                 firstRefPosJ,
                                 refBufferStride,
                                 refBufferFullList0Stride,
                                 firstRefTempDst,
                                 BLOCK_SIZE_64);

        ptr_list0        = firstRefTempDst;
        ptr_list0_stride = BLOCK_SIZE_64;
    }

    if (sub_position_type[secondFracPos] != 2) {
        select_buffer_luma(pu_index,
                           secondFracPos,
                           pu_width,
                           pu_height,
                           secondRefInteger,
                           secondRefPosB,
                           secondRefPosh_,
                           secondRefPosJ,
                           refBufferStride,
                           refBufferFullList1Stride,
                           &ptr_list1,
                           &ptr_list1_stride);
    } else {
        // uni-prediction List1 luma
        // doing the luma interpolation
        quarder_pel_compensation(pu_index,
                                 secondFracPos,
                                 pu_width,
                                 pu_height,
                                 secondRefInteger,
                                 secondRefPosB,
                                 secondRefPosh_,
                                 secondRefPosJ,
                                 refBufferStride,
                                 refBufferFullList1Stride,
                                 secondRefTempDst,
                                 BLOCK_SIZE_64);

        ptr_list1        = secondRefTempDst;
        ptr_list1_stride = BLOCK_SIZE_64;
    }

    // bi-pred luma
    me_candidate->distortion = (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                                   ? nxm_sad_avg_kernel(sourcePic,
                                                        luma_stride << 1,
                                                        ptr_list0,
                                                        ptr_list0_stride << 1,
                                                        ptr_list1,
                                                        ptr_list1_stride << 1,
                                                        pu_height >> 1,
                                                        pu_width)
                                         << 1
                                   : nxm_sad_avg_kernel(sourcePic,
                                                        luma_stride,
                                                        ptr_list0,
                                                        ptr_list0_stride,
                                                        ptr_list1,
                                                        ptr_list1_stride,
                                                        pu_height,
                                                        pu_width);

    return me_candidate->distortion;
}

/*******************************************
 * BiPredictionComponsation
 *   performs componsation fro List 0 and
 *   List1 Candidates and then compute the
 *   average
 *******************************************/
EbErrorType biprediction_compensation(MeContext *context_ptr, uint32_t pu_index,
                                      MePredUnit *me_candidate, uint32_t firstList,
                                      uint8_t first_list_ref_pic_idx, uint32_t firstRefMv,
                                      uint32_t secondList, uint8_t second_list_ref_pic_idx,
                                      uint32_t secondRefMv) {
    EbErrorType return_error = EB_ErrorNone;

    int16_t first_ref_pos_x;
    int16_t first_ref_pos_y;
    int16_t first_ref_integ_pos_x;
    int16_t first_ref_integ_pos_y;
    uint8_t first_ref_frac_pos_x;
    uint8_t first_ref_frac_pos_y;
    uint8_t first_ref_frac_pos;
    int32_t x_first_search_index;
    int32_t y_first_search_index;
    int32_t first_search_region_index_pos_integ;
    int32_t first_search_region_index_pos_b;
    int32_t first_search_region_index_pos_h;
    int32_t first_search_region_index_pos_j;

    int16_t second_ref_pos_x;
    int16_t second_ref_pos_y;
    int16_t second_ref_integ_pos_x;
    int16_t second_ref_integ_pos_y;
    uint8_t second_ref_frac_pos_x;
    uint8_t second_ref_frac_pos_y;
    uint8_t second_ref_frac_pos;
    int32_t x_second_search_index;
    int32_t y_second_search_index;
    int32_t second_search_region_index_pos_integ;
    int32_t second_search_region_index_pos_b;
    int32_t second_search_region_index_pos_h;
    int32_t second_search_region_index_pos_j;

    uint32_t pu_shift_x_index = pu_search_index_map[pu_index][0];
    uint32_t pu_shift_y_index = pu_search_index_map[pu_index][1];

    const uint32_t blk_sb_buffer_index =
        pu_shift_x_index + pu_shift_y_index * context_ptr->sb_src_stride;

    me_candidate->prediction_direction = BI_PRED;

    // First refrence
    // Set Candidate information
    first_ref_pos_x            = _MVXT(firstRefMv);
    first_ref_pos_y            = _MVYT(firstRefMv);
    me_candidate->ref_index[0] = (uint8_t)first_list_ref_pic_idx;
    me_candidate->ref0_list    = (uint8_t)firstList;

    first_ref_integ_pos_x = (first_ref_pos_x >> 2);
    first_ref_integ_pos_y = (first_ref_pos_y >> 2);
    first_ref_frac_pos_x  = (uint8_t)first_ref_pos_x & 0x03;
    first_ref_frac_pos_y  = (uint8_t)first_ref_pos_y & 0x03;

    first_ref_frac_pos   = (uint8_t)(first_ref_frac_pos_x + (first_ref_frac_pos_y << 2));
    x_first_search_index = (int32_t)first_ref_integ_pos_x -
                           context_ptr->x_search_area_origin[firstList][first_list_ref_pic_idx];
    y_first_search_index = (int32_t)first_ref_integ_pos_y -
                           context_ptr->y_search_area_origin[firstList][first_list_ref_pic_idx];
    first_search_region_index_pos_integ =
        (int32_t)(x_first_search_index + (ME_FILTER_TAP >> 1)) +
        (int32_t)context_ptr->interpolated_full_stride[firstList][first_list_ref_pic_idx] *
            (int32_t)(y_first_search_index + (ME_FILTER_TAP >> 1));

    first_search_region_index_pos_b = (int32_t)(x_first_search_index + (ME_FILTER_TAP >> 1) - 1) +
                                      (int32_t)context_ptr->interpolated_stride *
                                          (int32_t)(y_first_search_index + (ME_FILTER_TAP >> 1));
    first_search_region_index_pos_h =
        (int32_t)(x_first_search_index + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(y_first_search_index + (ME_FILTER_TAP >> 1) - 1);
    first_search_region_index_pos_j =
        (int32_t)(x_first_search_index + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(y_first_search_index + (ME_FILTER_TAP >> 1) - 1);

    // Second refrence

    // Set Candidate information
    second_ref_pos_x           = _MVXT(secondRefMv);
    second_ref_pos_y           = _MVYT(secondRefMv);
    me_candidate->ref_index[1] = (uint8_t)second_list_ref_pic_idx;
    me_candidate->ref1_list    = (uint8_t)secondList;
    second_ref_integ_pos_x     = (second_ref_pos_x >> 2);
    second_ref_integ_pos_y     = (second_ref_pos_y >> 2);
    second_ref_frac_pos_x      = (uint8_t)second_ref_pos_x & 0x03;
    second_ref_frac_pos_y      = (uint8_t)second_ref_pos_y & 0x03;

    second_ref_frac_pos   = (uint8_t)(second_ref_frac_pos_x + (second_ref_frac_pos_y << 2));
    x_second_search_index = second_ref_integ_pos_x -
                            context_ptr->x_search_area_origin[secondList][second_list_ref_pic_idx];
    y_second_search_index = second_ref_integ_pos_y -
                            context_ptr->y_search_area_origin[secondList][second_list_ref_pic_idx];
    second_search_region_index_pos_integ =
        (int32_t)(x_second_search_index + (ME_FILTER_TAP >> 1)) +
        (int32_t)context_ptr->interpolated_full_stride[secondList][second_list_ref_pic_idx] *
            (int32_t)(y_second_search_index + (ME_FILTER_TAP >> 1));
    second_search_region_index_pos_b = (int32_t)(x_second_search_index + (ME_FILTER_TAP >> 1) - 1) +
                                       (int32_t)context_ptr->interpolated_stride *
                                           (int32_t)(y_second_search_index + (ME_FILTER_TAP >> 1));
    second_search_region_index_pos_h =
        (int32_t)(x_second_search_index + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(y_second_search_index + (ME_FILTER_TAP >> 1) - 1);
    second_search_region_index_pos_j =
        (int32_t)(x_second_search_index + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(y_second_search_index + (ME_FILTER_TAP >> 1) - 1);

    uint32_t n_index;

    if (pu_index > 200)
        n_index = pu_index;
    else if (pu_index > 184)
        n_index = tab8x32[pu_index - 185] + 185;
    else if (pu_index > 168)
        n_index = tab32x8[pu_index - 169] + 169;
    else if (pu_index > 136)
        n_index = tab8x16[pu_index - 137] + 137;
    else if (pu_index > 128)
        n_index = tab16x32[pu_index - 129] + 129;
    else if (pu_index > 126)
        n_index = pu_index;
    else if (pu_index > 94)
        n_index = tab16x8[pu_index - 95] + 95;
    else if (pu_index > 86)
        n_index = tab32x16[pu_index - 87] + 87;
    else if (pu_index > 84)
        n_index = pu_index;
    else if (pu_index > 20)
        n_index = tab8x8[pu_index - 21] + 21;
    else if (pu_index > 4)
        n_index = tab16x16[pu_index - 5] + 5;
    else
        n_index = pu_index;
    context_ptr->p_sb_bipred_sad[n_index] =

        bi_pred_averaging(
            context_ptr,
            me_candidate,
            pu_index,
            &(context_ptr->sb_src_ptr[blk_sb_buffer_index]),
            context_ptr->sb_src_stride,
            first_ref_frac_pos,
            second_ref_frac_pos,
            partition_width[pu_index],
            partition_height[pu_index],
            &(context_ptr->integer_buffer_ptr[firstList][first_list_ref_pic_idx]
                                             [first_search_region_index_pos_integ]),
            &(context_ptr->pos_b_buffer[firstList][first_list_ref_pic_idx]
                                       [first_search_region_index_pos_b]),
            &(context_ptr->pos_h_buffer[firstList][first_list_ref_pic_idx]
                                       [first_search_region_index_pos_h]),
            &(context_ptr->pos_j_buffer[firstList][first_list_ref_pic_idx]
                                       [first_search_region_index_pos_j]),
            &(context_ptr->integer_buffer_ptr[secondList][second_list_ref_pic_idx]
                                             [second_search_region_index_pos_integ]),
            &(context_ptr->pos_b_buffer[secondList][second_list_ref_pic_idx]
                                       [second_search_region_index_pos_b]),
            &(context_ptr->pos_h_buffer[secondList][second_list_ref_pic_idx]
                                       [second_search_region_index_pos_h]),
            &(context_ptr->pos_j_buffer[secondList][second_list_ref_pic_idx]
                                       [second_search_region_index_pos_j]),
            context_ptr->interpolated_stride,
            context_ptr->interpolated_full_stride[firstList][first_list_ref_pic_idx],
            context_ptr->interpolated_full_stride[secondList][second_list_ref_pic_idx],
            &(context_ptr->one_d_intermediate_results_buf0[0]),
            &(context_ptr->one_d_intermediate_results_buf1[0]));

    return return_error;
}

uint8_t skip_bi_pred(PictureParentControlSet *pcs_ptr, uint8_t ref_type, uint8_t ref_type_table[7],
                     uint8_t size) {
    if (!pcs_ptr->prune_unipred_at_me) return 1;

    uint8_t allow_cand = 0;
    uint8_t ref_idx;
    for (ref_idx = 0; ref_idx < MIN(PRUNE_REF_ME_TH, size); ref_idx++) {
        if (ref_type == ref_type_table[ref_idx]) allow_cand = 1;
    }
    return allow_cand;
}

/*******************************************
 * bi_prediction_search
 *   performs Bi-Prediction Search (SB)
 *******************************************/
// This function enables all 16 Bipred candidates when MRP is ON
EbErrorType bi_prediction_search(SequenceControlSet *scs_ptr, MeContext *context_ptr,
                                 uint32_t pu_index, uint8_t cand_index,
                                 uint32_t activeRefPicFirstLisNum,
                                 uint32_t activeRefPicSecondLisNum,
                                 uint8_t *total_me_candidate_index, uint8_t ref_type_table[7],
                                 PictureParentControlSet *pcs_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    uint32_t first_list_ref_pict_idx;
    uint32_t second_list_ref_pict_idx;

    (void)pcs_ptr;

    uint32_t n_index;

    if (pu_index > 200)
        n_index = pu_index;
    else if (pu_index > 184)
        n_index = tab8x32[pu_index - 185] + 185;
    else if (pu_index > 168)
        n_index = tab32x8[pu_index - 169] + 169;
    else if (pu_index > 136)
        n_index = tab8x16[pu_index - 137] + 137;
    else if (pu_index > 128)
        n_index = tab16x32[pu_index - 129] + 129;
    else if (pu_index > 126)
        n_index = pu_index;
    else if (pu_index > 94)
        n_index = tab16x8[pu_index - 95] + 95;
    else if (pu_index > 86)
        n_index = tab32x16[pu_index - 87] + 87;
    else if (pu_index > 84)
        n_index = pu_index;
    else if (pu_index > 20)
        n_index = tab8x8[pu_index - 21] + 21;
    else if (pu_index > 4)
        n_index = tab16x16[pu_index - 5] + 5;
    else
        n_index = pu_index;
    // NM: Inter list bipred.
    //(LAST,BWD) , (LAST,ALT)  and (LAST,ALT2)
    //(LAST2,BWD), (LAST2,ALT) and (LAST2,ALT2)
    //(LAST3,BWD), (LAST3,ALT) and (LAST3,ALT2)
    //(GOLD,BWD) , (GOLD,ALT)  and (GOLD,ALT2)
    for (first_list_ref_pict_idx = 0; first_list_ref_pict_idx < activeRefPicFirstLisNum;
         first_list_ref_pict_idx++) {
        for (second_list_ref_pict_idx = 0; second_list_ref_pict_idx < activeRefPicSecondLisNum;
             second_list_ref_pict_idx++) {
            {
                uint8_t to_inject_ref_type_0 =
                    svt_get_ref_frame_type(REF_LIST_0, first_list_ref_pict_idx);
                uint8_t to_inject_ref_type_1 =
                    svt_get_ref_frame_type(REF_LIST_1, second_list_ref_pict_idx);
                uint8_t add_bi = skip_bi_pred(pcs_ptr, to_inject_ref_type_0, ref_type_table, *total_me_candidate_index);
                add_bi += skip_bi_pred(pcs_ptr, to_inject_ref_type_1, ref_type_table, *total_me_candidate_index);
                //if one of the references is skipped at ME, do not consider bi for this cand
                if (context_ptr->hme_results[REF_LIST_0][first_list_ref_pict_idx].do_ref == 0 ||
                    context_ptr->hme_results[REF_LIST_1][second_list_ref_pict_idx].do_ref == 0  )
                    add_bi = 0;
                if (add_bi) {
                    biprediction_compensation(
                        context_ptr,
                        pu_index,
                        &(context_ptr->me_candidate[cand_index].pu[pu_index]),
                        REFERENCE_PIC_LIST_0,
                        first_list_ref_pict_idx,
                        context_ptr
                            ->p_sb_best_mv[REFERENCE_PIC_LIST_0][first_list_ref_pict_idx][n_index],
                        REFERENCE_PIC_LIST_1,
                        second_list_ref_pict_idx,
                        context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_1][second_list_ref_pict_idx]
                                                 [n_index]);

                    cand_index++;
                }
            }
        }
    }

    if (scs_ptr->mrp_mode == 0) {
        // NM: Within list 0    bipred: (LAST,LAST2)    (LAST,LAST3) (LAST,GOLD)
        for (first_list_ref_pict_idx = 1; first_list_ref_pict_idx < activeRefPicFirstLisNum;
             first_list_ref_pict_idx++) {
            uint8_t to_inject_ref_type_0 =
                svt_get_ref_frame_type(REF_LIST_0, first_list_ref_pict_idx);
            uint8_t add_bi = skip_bi_pred(pcs_ptr, to_inject_ref_type_0, ref_type_table,*total_me_candidate_index);
            //if one of the references is skipped at ME, do not consider bi for this cand
            if (context_ptr->hme_results[REF_LIST_0][0].do_ref == 0 ||
                context_ptr->hme_results[REF_LIST_0][first_list_ref_pict_idx].do_ref == 0)
                add_bi = 0;
            if (add_bi) {
                biprediction_compensation(
                    context_ptr,
                    pu_index,
                    &(context_ptr->me_candidate[cand_index].pu[pu_index]),
                    REFERENCE_PIC_LIST_0,
                    0,
                    context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_0][0][n_index],
                    REFERENCE_PIC_LIST_0,
                    first_list_ref_pict_idx,
                    context_ptr
                        ->p_sb_best_mv[REFERENCE_PIC_LIST_0][first_list_ref_pict_idx][n_index]);

                cand_index++;
            }
        }
        // NM: Within list 1    bipred: (BWD, ALT)
        for (second_list_ref_pict_idx = 1;
             second_list_ref_pict_idx < MIN(activeRefPicSecondLisNum, 1);
             second_list_ref_pict_idx++) {
            uint8_t to_inject_ref_type_0 =
                svt_get_ref_frame_type(REF_LIST_0, first_list_ref_pict_idx);
            uint8_t add_bi = skip_bi_pred(pcs_ptr, to_inject_ref_type_0, ref_type_table,*total_me_candidate_index);
            //if one of the references is skipped at ME, do not consider bi for this cand
            if (context_ptr->hme_results[REF_LIST_1][0].do_ref == 0 ||
                context_ptr->hme_results[REF_LIST_1][first_list_ref_pict_idx].do_ref == 0)
                add_bi = 0;
            if (add_bi) {
                biprediction_compensation(
                    context_ptr,
                    pu_index,
                    &(context_ptr->me_candidate[cand_index].pu[pu_index]),
                    REFERENCE_PIC_LIST_1,
                    0,
                    context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_1][0][n_index],
                    REFERENCE_PIC_LIST_1,
                    second_list_ref_pict_idx,
                    context_ptr
                        ->p_sb_best_mv[REFERENCE_PIC_LIST_1][second_list_ref_pict_idx][n_index]);

                cand_index++;
            }
        }
    }
    *total_me_candidate_index = cand_index;

    return return_error;
}

// Nader - to be replaced by loock-up table
/*******************************************
 * get_me_info_index
 *   search the correct index of the motion
 *   info that corresponds to the input
 *   md candidate
 *******************************************/
uint32_t get_me_info_index(uint32_t max_me_block, const BlockGeom *blk_geom, uint32_t geom_offset_x,
                           uint32_t geom_offset_y) {
    // search for motion info
    uint32_t block_index;
    uint32_t me_info_index = 0xFFFFFFF;

    for (block_index = 0; block_index < max_me_block; block_index++) {
        if ((blk_geom->bwidth == partition_width[block_index]) &&
            (blk_geom->bheight == partition_height[block_index]) &&
            ((blk_geom->origin_x - geom_offset_x) == pu_search_index_map[block_index][0]) &&
            ((blk_geom->origin_y - geom_offset_y) == pu_search_index_map[block_index][1])) {
            me_info_index = block_index;
            break;
        }
    }
    return me_info_index;
}

#define NSET_CAND(me_pu_result, num, dist, dir)                      \
    (me_pu_result)->distortion_direction[(num)].distortion = (dist); \
    (me_pu_result)->distortion_direction[(num)].direction  = (dir);

EbErrorType check_00_center(EbPictureBufferDesc *ref_pic_ptr, MeContext *context_ptr,
                            uint32_t sb_origin_x, uint32_t sb_origin_y, uint32_t sb_width,
                            uint32_t sb_height, int16_t *x_search_center, int16_t *y_search_center)

{
    EbErrorType return_error = EB_ErrorNone;
    uint32_t    search_region_index, zero_mv_sad, hme_mv_sad, hme_mvd_rate;
    uint64_t    hme_mv_cost, zero_mv_cost, search_center_cost;
    int16_t     origin_x      = (int16_t)sb_origin_x;
    int16_t     origin_y      = (int16_t)sb_origin_y;
    uint32_t    subsample_sad = 1;
    int16_t     pad_width     = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t     pad_height    = (int16_t)BLOCK_SIZE_64 - 1;

    search_region_index = (int16_t)ref_pic_ptr->origin_x + origin_x +
                          ((int16_t)ref_pic_ptr->origin_y + origin_y) * ref_pic_ptr->stride_y;

    zero_mv_sad = nxm_sad_kernel(context_ptr->sb_src_ptr,
                                 context_ptr->sb_src_stride << subsample_sad,
                                 &(ref_pic_ptr->buffer_y[search_region_index]),
                                 ref_pic_ptr->stride_y << subsample_sad,
                                 sb_height >> subsample_sad,
                                 sb_width);

    zero_mv_sad = zero_mv_sad << subsample_sad;

    // FIX
    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    *x_search_center =
        ((origin_x + *x_search_center) < -pad_width) ? -pad_width - origin_x : *x_search_center;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    *x_search_center =
        ((origin_x + *x_search_center) > (int16_t)ref_pic_ptr->width - 1)
            ? *x_search_center - ((origin_x + *x_search_center) - ((int16_t)ref_pic_ptr->width - 1))
            : *x_search_center;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    *y_search_center =
        ((origin_y + *y_search_center) < -pad_height) ? -pad_height - origin_y : *y_search_center;
    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    *y_search_center = ((origin_y + *y_search_center) > (int16_t)ref_pic_ptr->height - 1)
                           ? *y_search_center - ((origin_y + *y_search_center) -
                                                 ((int16_t)ref_pic_ptr->height - 1))
                           : *y_search_center;
    ///

    zero_mv_cost = zero_mv_sad << COST_PRECISION;
    search_region_index =
        (int16_t)(ref_pic_ptr->origin_x + origin_x) + *x_search_center +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + *y_search_center) * ref_pic_ptr->stride_y;

    hme_mv_sad = nxm_sad_kernel(context_ptr->sb_src_ptr,
                               context_ptr->sb_src_stride << subsample_sad,
                               &(ref_pic_ptr->buffer_y[search_region_index]),
                               ref_pic_ptr->stride_y << subsample_sad,
                               sb_height >> subsample_sad,
                               sb_width);

    hme_mv_sad = hme_mv_sad << subsample_sad;

    hme_mvd_rate = 0;

    hme_mv_cost = (hme_mv_sad << COST_PRECISION) +
                  (((context_ptr->lambda * hme_mvd_rate) + MD_OFFSET) >> MD_SHIFT);
    search_center_cost = MIN(zero_mv_cost, hme_mv_cost);

    *x_search_center = (search_center_cost == zero_mv_cost) ? 0 : *x_search_center;
    *y_search_center = (search_center_cost == zero_mv_cost) ? 0 : *y_search_center;

    return return_error;
}

EbErrorType su_pel_enable(MeContext *context_ptr, PictureParentControlSet *pcs_ptr,
                          uint32_t list_index, uint32_t ref_pic_index,
                          EbBool *enable_half_pel_32x32, EbBool *enable_half_pel_16x16,
                          EbBool *enable_half_pel_8x8) {
    EbErrorType return_error = EB_ErrorNone;

    uint32_t mv_mag_32x32  = 0;
    uint32_t mv_mag_16x16  = 0;
    uint32_t mv_mag_8x8    = 0;
    uint32_t avgsad_32x32  = 0;
    uint32_t avgsad_16x16  = 0;
    uint32_t avgsad_8x8    = 0;
    uint32_t avg_mvx_32x32 = 0;
    uint32_t avg_mvy_32x32 = 0;
    uint32_t avg_mvx_16x16 = 0;
    uint32_t avg_mvy_16x16 = 0;
    uint32_t avg_mvx_8x8   = 0;
    uint32_t avg_mvy_8x8   = 0;

    avg_mvx_32x32 =
        (_MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_0]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_1]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_2]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_3])) >>
        2;
    avg_mvy_32x32 =
        (_MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_0]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_1]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_2]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_3])) >>
        2;
    mv_mag_32x32 = SQR(avg_mvx_32x32) + SQR(avg_mvy_32x32);

    avg_mvx_16x16 =
        (_MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_0]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_1]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_2]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_3]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_4]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_5]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_6]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_7]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_8]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_9]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_10]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_11]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_12]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_13]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_14]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_15])) >>
        4;
    avg_mvy_16x16 =
        (_MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_0]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_1]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_2]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_3]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_4]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_5]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_6]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_7]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_8]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_9]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_10]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_11]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_12]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_13]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_14]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_15])) >>
        4;
    mv_mag_16x16 = SQR(avg_mvx_16x16) + SQR(avg_mvy_16x16);

    avg_mvx_8x8 =
        (_MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_1]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_2]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_3]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_4]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_5]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_6]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_7]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_8]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_9]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_10]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_11]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_12]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_13]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_14]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_15]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_16]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_17]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_18]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_19]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_20]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_21]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_22]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_23]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_24]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_25]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_26]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_27]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_28]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_29]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_30]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_31]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_32]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_33]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_34]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_35]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_36]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_37]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_38]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_39]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_40]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_41]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_42]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_43]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_44]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_45]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_46]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_47]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_48]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_49]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_50]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_51]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_52]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_53]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_54]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_55]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_56]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_57]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_58]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_59]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_60]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_61]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_62]) +
         _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_63])) >>
        6;
    avg_mvy_8x8 =
        (_MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_1]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_2]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_3]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_4]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_5]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_6]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_7]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_8]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_9]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_10]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_11]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_12]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_13]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_14]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_15]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_16]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_17]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_18]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_19]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_20]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_21]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_22]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_23]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_24]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_25]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_26]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_27]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_28]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_29]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_30]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_31]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_32]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_33]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_34]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_35]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_36]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_37]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_38]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_39]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_40]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_41]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_42]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_43]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_44]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_45]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_46]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_47]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_48]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_49]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_50]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_51]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_52]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_53]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_54]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_55]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_56]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_57]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_58]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_59]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_60]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_61]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_62]) +
         _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_63])) >>
        6;
    mv_mag_8x8 = SQR(avg_mvx_8x8) + SQR(avg_mvy_8x8);

    avgsad_32x32 =
        (context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_0] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_1] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_2] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_3]) >>
        2;
    avgsad_16x16 =
        (context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_0] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_1] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_2] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_3] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_4] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_5] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_6] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_7] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_8] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_9] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_10] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_11] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_12] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_13] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_14] +
         context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_15]) >>
        4;
    avgsad_8x8 = (context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_1] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_2] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_3] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_4] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_5] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_6] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_7] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_8] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_9] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_10] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_11] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_12] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_13] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_14] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_15] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_16] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_17] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_18] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_19] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_20] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_21] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_22] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_23] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_24] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_25] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_26] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_27] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_28] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_29] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_30] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_31] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_32] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_33] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_34] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_35] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_36] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_37] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_38] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_39] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_40] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_41] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_42] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_43] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_44] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_45] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_46] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_47] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_48] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_49] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_50] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_51] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_52] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_53] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_54] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_55] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_56] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_57] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_58] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_59] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_60] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_61] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_62] +
                  context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_63]) >>
                 6;

    if (pcs_ptr->temporal_layer_index == 0) {
        // 32x32
        if ((mv_mag_32x32 < SQR(48)) && (avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_0
        else if ((mv_mag_32x32 < SQR(48)) && !(avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_FALSE; // CLASS_1
        else if (!(mv_mag_32x32 < SQR(48)) && (avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_2
        else
            *enable_half_pel_32x32 = EB_FALSE; // CLASS_3
        // 16x16
        if ((mv_mag_16x16 < SQR(48)) && (avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_FALSE; // CLASS_0
        else if ((mv_mag_16x16 < SQR(48)) && !(avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_16x16 < SQR(48)) && (avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_FALSE; // CLASS_2
        else
            *enable_half_pel_16x16 = EB_TRUE; // CLASS_3
        // 8x8
        if ((mv_mag_8x8 < SQR(48)) && (avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_FALSE; // CLASS_0
        else if ((mv_mag_8x8 < SQR(48)) && !(avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_8x8 < SQR(48)) && (avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_FALSE; // CLASS_2
        else
            *enable_half_pel_8x8 = EB_TRUE; // CLASS_3
    }

    else if (pcs_ptr->temporal_layer_index == 1) {
        // 32x32
        if ((mv_mag_32x32 < SQR(32)) && (avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_0
        else if ((mv_mag_32x32 < SQR(32)) && !(avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_FALSE; // CLASS_1
        else if (!(mv_mag_32x32 < SQR(32)) && (avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_2
        else
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_3
        // 16x16
        if ((mv_mag_16x16 < SQR(32)) && (avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_FALSE; // CLASS_0
        else if ((mv_mag_16x16 < SQR(32)) && !(avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_16x16 < SQR(32)) && (avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_FALSE; // CLASS_2
        else
            *enable_half_pel_16x16 = EB_TRUE; // CLASS_3
        // 8x8
        if ((mv_mag_8x8 < SQR(32)) && (avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_FALSE; // CLASS_0
        else if ((mv_mag_8x8 < SQR(32)) && !(avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_8x8 < SQR(32)) && (avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_FALSE; // CLASS_2
        else
            *enable_half_pel_8x8 = EB_TRUE; // CLASS_3
    } else if (pcs_ptr->temporal_layer_index == 2) {
        // 32x32
        if ((mv_mag_32x32 < SQR(80)) && (avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_0
        else if ((mv_mag_32x32 < SQR(80)) && !(avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_FALSE; // CLASS_1
        else if (!(mv_mag_32x32 < SQR(80)) && (avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_2
        else
            *enable_half_pel_32x32 = EB_FALSE; // CLASS_3
        // 16x16
        if ((mv_mag_16x16 < SQR(80)) && (avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_FALSE; // CLASS_0
        else if ((mv_mag_16x16 < SQR(80)) && !(avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_16x16 < SQR(80)) && (avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_FALSE; // CLASS_2
        else
            *enable_half_pel_16x16 = EB_TRUE; // CLASS_3
        // 8x8
        if ((mv_mag_8x8 < SQR(80)) && (avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_FALSE; // CLASS_0
        else if ((mv_mag_8x8 < SQR(80)) && !(avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_8x8 < SQR(80)) && (avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_FALSE; // CLASS_2
        else
            *enable_half_pel_8x8 = EB_TRUE; // CLASS_3
    } else {
        // 32x32
        if ((mv_mag_32x32 < SQR(48)) && (avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_0
        else if ((mv_mag_32x32 < SQR(48)) && !(avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_32x32 < SQR(48)) && (avgsad_32x32 < 32 * 32 * 6))
            *enable_half_pel_32x32 = EB_TRUE; // CLASS_2
        else
            *enable_half_pel_32x32 = EB_FALSE; // CLASS_3
        // 16x16
        if ((mv_mag_16x16 < SQR(48)) && (avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_FALSE; // CLASS_0
        else if ((mv_mag_16x16 < SQR(48)) && !(avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_16x16 < SQR(48)) && (avgsad_16x16 < 16 * 16 * 2))
            *enable_half_pel_16x16 = EB_FALSE; // CLASS_2
        else
            *enable_half_pel_16x16 = EB_TRUE; // CLASS_3
        // 8x8
        if ((mv_mag_8x8 < SQR(48)) && (avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_FALSE; // CLASS_0
        else if ((mv_mag_8x8 < SQR(48)) && !(avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_TRUE; // CLASS_1
        else if (!(mv_mag_8x8 < SQR(48)) && (avgsad_8x8 < 8 * 8 * 2))
            *enable_half_pel_8x8 = EB_FALSE; // CLASS_2
        else
            *enable_half_pel_8x8 = EB_FALSE; // EB_TRUE; //CLASS_3
    }

    return return_error;
}

static void hme_mv_center_check(EbPictureBufferDesc *ref_pic_ptr, MeContext *context_ptr,
                                int16_t *xsc, int16_t *ysc, uint32_t list_index, int16_t origin_x,
                                int16_t origin_y, uint32_t sb_width, uint32_t sb_height) {
    // Search for (-srx/2, 0),  (+srx/2, 0),  (0, -sry/2), (0, +sry/2),
    /*
    |------------C-------------|
    |--------------------------|
    |--------------------------|
    A            0             b
    |--------------------------|
    |--------------------------|
    |------------D-------------|
    */
    uint32_t search_region_index;
    int16_t  search_center_x = *xsc;
    int16_t  search_center_y = *ysc;
    uint64_t best_cost;
    uint64_t direct_mv_cost = 0xFFFFFFFFFFFFF;
    uint8_t  sparce_scale   = 1;
    int16_t  pad_width      = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t  pad_height     = (int16_t)BLOCK_SIZE_64 - 1;
    // O pos

    search_region_index = (int16_t)ref_pic_ptr->origin_x + origin_x +
                          ((int16_t)ref_pic_ptr->origin_y + origin_y) * ref_pic_ptr->stride_y;

    uint32_t sub_sampled_sad = 1;
    uint64_t zero_mv_sad     = nxm_sad_kernel(context_ptr->sb_src_ptr,
                                          context_ptr->sb_src_stride << sub_sampled_sad,
                                          &(ref_pic_ptr->buffer_y[search_region_index]),
                                          ref_pic_ptr->stride_y << sub_sampled_sad,
                                          sb_height >> sub_sampled_sad,
                                          sb_width);

    zero_mv_sad = zero_mv_sad << sub_sampled_sad;

    uint64_t zero_mv_cost = zero_mv_sad << COST_PRECISION;

    // A pos
    search_center_x = 0 - (context_ptr->hme_level0_total_search_area_width * sparce_scale);
    search_center_y = 0;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) < -pad_width) ? -pad_width - origin_x : search_center_x;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
            ? search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1))
            : search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) < -pad_height) ? -pad_height - origin_y : search_center_y;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
            ? search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1))
            : search_center_y;

    search_region_index =
        (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) * ref_pic_ptr->stride_y;

    uint64_t mv_a_sad = nxm_sad_kernel(context_ptr->sb_src_ptr,
                                       context_ptr->sb_src_stride << sub_sampled_sad,
                                       &(ref_pic_ptr->buffer_y[search_region_index]),
                                       ref_pic_ptr->stride_y << sub_sampled_sad,
                                       sb_height >> sub_sampled_sad,
                                       sb_width);

    mv_a_sad = mv_a_sad << sub_sampled_sad;

    uint64_t mv_a_cost = mv_a_sad << COST_PRECISION;

    // b pos
    search_center_x = (context_ptr->hme_level0_total_search_area_width * sparce_scale);
    search_center_y = 0;
    ///////////////// correct
    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) < -pad_width) ? -pad_width - origin_x : search_center_x;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
            ? search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1))
            : search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) < -pad_height) ? -pad_height - origin_y : search_center_y;
    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
            ? search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1))
            : search_center_y;

    search_region_index =
        (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) * ref_pic_ptr->stride_y;

    uint64_t mv_b_sad = nxm_sad_kernel(context_ptr->sb_src_ptr,
                                       context_ptr->sb_src_stride << sub_sampled_sad,
                                       &(ref_pic_ptr->buffer_y[search_region_index]),
                                       ref_pic_ptr->stride_y << sub_sampled_sad,
                                       sb_height >> sub_sampled_sad,
                                       sb_width);

    mv_b_sad = mv_b_sad << sub_sampled_sad;

    uint64_t mv_b_cost = mv_b_sad << COST_PRECISION;
    // C pos
    search_center_x = 0;
    search_center_y = 0 - (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    ///////////////// correct
    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) < -pad_width) ? -pad_width - origin_x : search_center_x;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
            ? search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1))
            : search_center_x;

    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) < -pad_height) ? -pad_height - origin_y : search_center_y;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
            ? search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1))
            : search_center_y;

    search_region_index =
        (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) * ref_pic_ptr->stride_y;

    uint64_t mv_c_sad = nxm_sad_kernel(context_ptr->sb_src_ptr,
                                       context_ptr->sb_src_stride << sub_sampled_sad,
                                       &(ref_pic_ptr->buffer_y[search_region_index]),
                                       ref_pic_ptr->stride_y << sub_sampled_sad,
                                       sb_height >> sub_sampled_sad,
                                       sb_width);

    mv_c_sad = mv_c_sad << sub_sampled_sad;

    uint64_t mv_c_cost = mv_c_sad << COST_PRECISION;

    // D pos
    search_center_x = 0;
    search_center_y = (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) < -pad_width) ? -pad_width - origin_x : search_center_x;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
            ? search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1))
            : search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) < -pad_height) ? -pad_height - origin_y : search_center_y;
    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
            ? search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1))
            : search_center_y;
    search_region_index =
        (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) * ref_pic_ptr->stride_y;
    uint64_t mv_d_sad = nxm_sad_kernel(context_ptr->sb_src_ptr,
                                       context_ptr->sb_src_stride << sub_sampled_sad,
                                       &(ref_pic_ptr->buffer_y[search_region_index]),
                                       ref_pic_ptr->stride_y << sub_sampled_sad,
                                       sb_height >> sub_sampled_sad,
                                       sb_width);

    mv_d_sad = mv_d_sad << sub_sampled_sad;

    uint64_t mv_d_cost = mv_d_sad << COST_PRECISION;
    best_cost = MIN(zero_mv_cost,
                    MIN(mv_a_cost, MIN(mv_b_cost, MIN(mv_c_cost, MIN(mv_d_cost, direct_mv_cost)))));

    if (best_cost == zero_mv_cost) {
        search_center_x = 0;
        search_center_y = 0;
    } else if (best_cost == mv_a_cost) {
        search_center_x = 0 - (context_ptr->hme_level0_total_search_area_width * sparce_scale);
        search_center_y = 0;
    } else if (best_cost == mv_b_cost) {
        search_center_x = (context_ptr->hme_level0_total_search_area_width * sparce_scale);
        search_center_y = 0;
    } else if (best_cost == mv_c_cost) {
        search_center_x = 0;
        search_center_y = 0 - (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    } else if (best_cost == direct_mv_cost) {
        search_center_x = list_index ? 0 - (_MVXT(context_ptr->p_sb_best_mv[0][0][0]) >> 2) : 0;
        search_center_y = list_index ? 0 - (_MVYT(context_ptr->p_sb_best_mv[0][0][0]) >> 2) : 0;
    } else if (best_cost == mv_d_cost) {
        search_center_x = 0;
        search_center_y = (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    }

    else
        SVT_ERROR("no center selected");
    *xsc = search_center_x;
    *ysc = search_center_y;
}
/*
 swap the content of two MePredUnit structures
*/
void swap_me_candidate(MePredUnit *a, MePredUnit *b) {
    MePredUnit temp_ptr;
    temp_ptr = *a;
    *a       = *b;
    *b       = temp_ptr;
}
/*******************************************
 *   performs integer search motion estimation for
 all avaiable references frames
 *******************************************/
void integer_search_sb(
    PictureParentControlSet   *pcs_ptr,
    uint32_t                   sb_index,
    uint32_t                   sb_origin_x,
    uint32_t                   sb_origin_y,
    MeContext                 *context_ptr,
    EbPictureBufferDesc       *input_ptr) {

    SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
    int16_t picture_width = pcs_ptr->aligned_width;
    int16_t picture_height = pcs_ptr->aligned_height;
    uint32_t sb_width = (input_ptr->width - sb_origin_x) < BLOCK_SIZE_64
                            ? input_ptr->width - sb_origin_x
                            : BLOCK_SIZE_64;
    uint32_t sb_height = (input_ptr->height - sb_origin_y) < BLOCK_SIZE_64
                             ? input_ptr->height - sb_origin_y
                             : BLOCK_SIZE_64;
    int16_t pad_width  = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t pad_height = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t origin_x = (int16_t)sb_origin_x;
    int16_t origin_y = (int16_t)sb_origin_y;
    int16_t search_area_width;
    int16_t search_area_height;
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;
    int16_t  x_top_left_search_region;
    int16_t  y_top_left_search_region;
    uint32_t search_region_index;
    uint32_t num_of_list_to_search;
    uint32_t list_index;
    uint8_t ref_pic_index;
    uint8_t num_of_ref_pic_to_search;
    EbPaReferenceObject *reference_object; // input parameter, reference Object Ptr
    // Final ME Search Center
    int16_t x_search_center = 0;
    int16_t y_search_center = 0;
    EbPictureBufferDesc *ref_pic_ptr;
    num_of_list_to_search =
        (pcs_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
    if (context_ptr->me_alt_ref == EB_TRUE) num_of_list_to_search = 0;

    // Uni-Prediction motion estimation loop
    // List Loop
    for (list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
        if (context_ptr->me_alt_ref == EB_TRUE) {
            num_of_ref_pic_to_search = 1;
        } else {
            num_of_ref_pic_to_search = (pcs_ptr->slice_type == P_SLICE)
                ? pcs_ptr->ref_list0_count_try
                : (list_index == REF_LIST_0) ? pcs_ptr->ref_list0_count_try
                : pcs_ptr->ref_list1_count_try;

            reference_object =
                (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[0][0]->object_ptr;
        }

        // Ref Picture Loop
        for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            if (context_ptr->me_alt_ref == EB_TRUE) {
                reference_object = (EbPaReferenceObject *)context_ptr->alt_ref_reference_ptr;
            } else {
                if (num_of_list_to_search) {
                    reference_object =
                        (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[1][0]->object_ptr;
                }

                reference_object =
                    (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]
                        ->object_ptr;
            }

            ref_pic_ptr = (EbPictureBufferDesc *)reference_object->input_padded_picture_ptr;
            // Get hme results
            if (context_ptr->hme_results[list_index][ref_pic_index].do_ref == 0)
                continue;  //so will not get ME results for those references.
            x_search_center = context_ptr->hme_results[list_index][ref_pic_index].hme_sc_x;
            y_search_center = context_ptr->hme_results[list_index][ref_pic_index].hme_sc_y;
            search_area_width = context_ptr->search_area_width;
            search_area_height = context_ptr->search_area_height;

            uint16_t dist = (context_ptr->me_alt_ref == EB_TRUE) ?
                ABS((int16_t)(context_ptr->tf_frame_index - context_ptr->tf_index_center)) :
                ABS((int16_t)(pcs_ptr->picture_number -
                    pcs_ptr->ref_pic_poc_array[list_index][ref_pic_index]));
            // factor to slowdown the ME search region growth to MAX
            if (!pcs_ptr->sc_content_detected && context_ptr->me_alt_ref == 0) {
                int8_t round_up = ((dist%8) == 0) ? 0 : 1;
                dist = ((dist * 5) / 8) + round_up;
            }
            search_area_width = MIN((search_area_width*dist),context_ptr->max_me_search_width);
            search_area_height = MIN((search_area_height*dist),context_ptr->max_me_search_height);

            // Constrain x_ME to be a multiple of 8 (round up)
            // Update ME search reagion size based on hme-data
            if (context_ptr->reduce_me_sr_flag[list_index][ref_pic_index] == SC_HME_TH_STILL) {
                search_area_width = ((search_area_width / SC_SR_DENOM_STILL) + 7) & ~0x07;
                search_area_height = (search_area_height / SC_SR_DENOM_STILL);
            }
            else if (context_ptr->reduce_me_sr_flag[list_index][ref_pic_index] == SC_HME_TH_EASY) {
                search_area_width = ((search_area_width / SC_SR_DENOM_EASY) + 7) & ~0x07;
                search_area_height = (search_area_height / SC_SR_DENOM_EASY);
            }
            else if (context_ptr->reduce_me_sr_flag[list_index][ref_pic_index]) {
                search_area_width = ((search_area_width / 8) + 7) & ~0x07;
                search_area_height = (search_area_height / 8);
            }
            else {
                search_area_width = (search_area_width + 7) & ~0x07;
            }

            if ((x_search_center != 0 || y_search_center != 0) &&
                (pcs_ptr->is_used_as_reference_flag == EB_TRUE)) {
                check_00_center(ref_pic_ptr,
                                context_ptr,
                                sb_origin_x,
                                sb_origin_y,
                                sb_width,
                                sb_height,
                                &x_search_center,
                                &y_search_center);
            }
            x_search_area_origin = x_search_center - (search_area_width >> 1);
            y_search_area_origin = y_search_center - (search_area_height >> 1);

            if (scs_ptr->static_config.unrestricted_motion_vector == 0) {
                int tile_start_x = pcs_ptr->sb_params_array[sb_index].tile_start_x;
                int tile_end_x   = pcs_ptr->sb_params_array[sb_index].tile_end_x;
                // Correct the left edge of the Search Area if it is not on the
                // reference Picture
                x_search_area_origin = ((origin_x + x_search_area_origin) < tile_start_x)
                                           ? tile_start_x - origin_x
                                           : x_search_area_origin;
                search_area_width =
                    ((origin_x + x_search_area_origin) < tile_start_x)
                        ? search_area_width - (tile_start_x - (origin_x + x_search_area_origin))
                        : search_area_width;
                // Correct the right edge of the Search Area if its not on the
                // reference Picture
                x_search_area_origin =
                    ((origin_x + x_search_area_origin) > tile_end_x - 1)
                        ? x_search_area_origin -
                              ((origin_x + x_search_area_origin) - (tile_end_x - 1))
                        : x_search_area_origin;
                search_area_width =
                    ((origin_x + x_search_area_origin + search_area_width) > tile_end_x)
                        ? MAX(1,
                              search_area_width -
                                  ((origin_x + x_search_area_origin + search_area_width) -
                                   tile_end_x))
                        : search_area_width;
                // Constrain x_ME to be a multiple of 8 (round down as cropping
                // already performed)
                search_area_width =
                    (search_area_width < 8) ? search_area_width : search_area_width & ~0x07;
            } else {
                // Correct the left edge of the Search Area if it is not on the
                // reference Picture
                x_search_area_origin = ((origin_x + x_search_area_origin) < -pad_width)
                                           ? -pad_width - origin_x
                                           : x_search_area_origin;
                search_area_width =
                    ((origin_x + x_search_area_origin) < -pad_width)
                        ? search_area_width - (-pad_width - (origin_x + x_search_area_origin))
                        : search_area_width;
                // Correct the right edge of the Search Area if its not on the
                // reference Picture
                x_search_area_origin =
                    ((origin_x + x_search_area_origin) > picture_width - 1)
                        ? x_search_area_origin -
                              ((origin_x + x_search_area_origin) - (picture_width - 1))
                        : x_search_area_origin;

                search_area_width =
                    ((origin_x + x_search_area_origin + search_area_width) > picture_width)
                        ? MAX(1,
                              search_area_width -
                                  ((origin_x + x_search_area_origin + search_area_width) -
                                   picture_width))
                        : search_area_width;

                // Constrain x_ME to be a multiple of 8 (round down as cropping
                // already performed)
                search_area_width =
                    (search_area_width < 8) ? search_area_width : search_area_width & ~0x07;
            }
            if (scs_ptr->static_config.unrestricted_motion_vector == 0) {
                int tile_start_y = pcs_ptr->sb_params_array[sb_index].tile_start_y;
                int tile_end_y   = pcs_ptr->sb_params_array[sb_index].tile_end_y;

                // Correct the top edge of the Search Area if it is not on the
                // reference Picture
                y_search_area_origin = ((origin_y + y_search_area_origin) < tile_start_y)
                                           ? tile_start_y - origin_y
                                           : y_search_area_origin;

                search_area_height =
                    ((origin_y + y_search_area_origin) < tile_start_y)
                        ? search_area_height - (tile_start_y - (origin_y + y_search_area_origin))
                        : search_area_height;

                // Correct the bottom edge of the Search Area if its not on the
                // reference Picture
                y_search_area_origin =
                    ((origin_y + y_search_area_origin) > tile_end_y - 1)
                        ? y_search_area_origin -
                              ((origin_y + y_search_area_origin) - (tile_end_y - 1))
                        : y_search_area_origin;

                search_area_height =
                    (origin_y + y_search_area_origin + search_area_height > tile_end_y)
                        ? MAX(1,
                              search_area_height -
                                  ((origin_y + y_search_area_origin + search_area_height) -
                                   tile_end_y))
                        : search_area_height;
            } else {
                // Correct the top edge of the Search Area if it is not on the
                // reference Picture
                y_search_area_origin = ((origin_y + y_search_area_origin) < -pad_height)
                                           ? -pad_height - origin_y
                                           : y_search_area_origin;
                search_area_height =
                    ((origin_y + y_search_area_origin) < -pad_height)
                        ? search_area_height - (-pad_height - (origin_y + y_search_area_origin))
                        : search_area_height;
                // Correct the bottom edge of the Search Area if its not on the
                // reference Picture
                y_search_area_origin =
                    ((origin_y + y_search_area_origin) > picture_height - 1)
                        ? y_search_area_origin -
                              ((origin_y + y_search_area_origin) - (picture_height - 1))
                        : y_search_area_origin;
                search_area_height =
                    (origin_y + y_search_area_origin + search_area_height > picture_height)
                        ? MAX(1,
                              search_area_height -
                                  ((origin_y + y_search_area_origin + search_area_height) -
                                   picture_height))
                        : search_area_height;
            }
            context_ptr->x_search_area_origin[list_index][ref_pic_index] = x_search_area_origin;
            context_ptr->y_search_area_origin[list_index][ref_pic_index] = y_search_area_origin;
            context_ptr->adj_search_area_width  = search_area_width;
            context_ptr->adj_search_area_height = search_area_height;
            x_top_left_search_region = (int16_t)(ref_pic_ptr->origin_x + sb_origin_x) -
                                       (ME_FILTER_TAP >> 1) + x_search_area_origin;
            y_top_left_search_region = (int16_t)(ref_pic_ptr->origin_y + sb_origin_y) -
                                       (ME_FILTER_TAP >> 1) + y_search_area_origin;
            search_region_index =
                (x_top_left_search_region) + (y_top_left_search_region)*ref_pic_ptr->stride_y;
            context_ptr->integer_buffer_ptr[list_index][ref_pic_index] =
                &(ref_pic_ptr->buffer_y[search_region_index]);
            context_ptr->interpolated_full_stride[list_index][ref_pic_index] =
                ref_pic_ptr->stride_y;
            // Move to the top left of the search region
            x_top_left_search_region =
                (int16_t)(ref_pic_ptr->origin_x + sb_origin_x) + x_search_area_origin;
            y_top_left_search_region =
                (int16_t)(ref_pic_ptr->origin_y + sb_origin_y) + y_search_area_origin;
            search_region_index =
                x_top_left_search_region + y_top_left_search_region * ref_pic_ptr->stride_y;
            if (pcs_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
                initialize_buffer_32bits(
                            context_ptr->p_sb_best_sad[list_index][ref_pic_index],
                            52,
                            1,
                            MAX_SAD_VALUE);
                context_ptr->p_best_sad_64x64 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_sad_32x32 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_sad_16x16 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_sad_8x8 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_sad_64x32 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_sad_32x16 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_sad_16x8 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_sad_32x64 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_sad_16x32 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_sad_8x16 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_sad_32x8 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_sad_8x32 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_sad_64x16 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_sad_16x64 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x64_0]);

                        context_ptr->p_best_mv64x64 =
                            &(context_ptr
                                  ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_mv32x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_mv16x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_mv8x8 =
                            &(context_ptr
                                  ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_mv64x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_mv32x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_mv16x8 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_mv32x64 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_mv16x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_mv8x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_mv32x8 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_mv8x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_mv64x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_mv16x64 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x64_0]);

                        context_ptr->p_best_ssd64x64 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_ssd32x32 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_ssd16x16 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_ssd8x8 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_ssd64x32 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_ssd32x16 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_ssd16x8 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_ssd32x64 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_ssd16x32 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_ssd8x16 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_ssd32x8 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_ssd8x32 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_ssd64x16 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_ssd16x64 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x64_0]);

                        open_loop_me_fullpel_search_sblock(context_ptr,
                                                           list_index,
                                                           ref_pic_index,
                                                           x_search_area_origin,
                                                           y_search_area_origin,
                                                           search_area_width,
                                                           search_area_height,
                                                           pcs_ptr->pic_depth_mode);

            }
            else {
                 initialize_buffer_32bits(
                            context_ptr->p_sb_best_sad[list_index][ref_pic_index],
                            21,
                            1,
                            MAX_SAD_VALUE);

                        context_ptr->full_quarter_pel_refinement = 0;

                        context_ptr->p_best_sad_64x64 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_sad_32x32 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_sad_16x16 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_sad_8x8 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);

                        context_ptr->p_best_mv64x64 =
                            &(context_ptr
                                  ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_mv32x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_mv16x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_mv8x8 =
                            &(context_ptr
                                  ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);

                        context_ptr->p_best_ssd64x64 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_ssd32x32 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_ssd16x16 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_ssd8x8 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
                        full_pel_search_sb(context_ptr,
                                           list_index,
                                           ref_pic_index,
                                           x_search_area_origin,
                                           y_search_area_origin,
                                           search_area_width,
                                           search_area_height);
            }
            context_ptr->x_search_area_origin[list_index][ref_pic_index] = x_search_area_origin;
            context_ptr->y_search_area_origin[list_index][ref_pic_index] = y_search_area_origin;
            context_ptr->sa_width[list_index][ref_pic_index] = search_area_width;
            context_ptr->sa_height[list_index][ref_pic_index] = search_area_height;
        }
    }
}

/*
  using previous stage ME results (Integer Search) for each reference
  frame. keep only the references that are close to the best reference.
*/
void prune_references_fp(
    PictureParentControlSet   *pcs_ptr,
    MeContext                 *context_ptr)
{
    HmeResults sorted[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    uint32_t num_of_cand_to_sort = MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH;
    uint8_t list_index, ref_pic_index;
    uint8_t num_of_ref_pic_to_search, num_of_list_to_search;
    uint32_t idx;
    uint32_t pu_index;
    num_of_list_to_search = (pcs_ptr->slice_type == P_SLICE)
        ? (uint32_t)REF_LIST_0
        : (uint32_t)REF_LIST_1;

    if (context_ptr->me_alt_ref == EB_TRUE)
        num_of_list_to_search = 0;

    for (list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {

        if (context_ptr->me_alt_ref == EB_TRUE) {
            num_of_ref_pic_to_search = 1;
        }
        else {
            num_of_ref_pic_to_search =
                (pcs_ptr->slice_type == P_SLICE)
                ? pcs_ptr->ref_list0_count_try
                : (list_index == REF_LIST_0)
                ? pcs_ptr->ref_list0_count_try
                : pcs_ptr->ref_list1_count_try;
        }
        // Ref Picture Loop
        for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            context_ptr->hme_results[list_index][ref_pic_index].hme_sad = 0;
            // Get hme results
            if (context_ptr->hme_results[list_index][ref_pic_index].do_ref == 0) {
                context_ptr->hme_results[list_index][ref_pic_index].hme_sad = MAX_SAD_VALUE * 64;
                continue;
            }
            context_ptr->p_best_sad_8x8 = &(context_ptr->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
            // 8x8   [64 partitions]
            for (pu_index = 0; pu_index < 64; ++pu_index) {
                idx = tab8x8[pu_index];
                context_ptr->hme_results[list_index][ref_pic_index].hme_sad += context_ptr->p_best_sad_8x8[idx];
            }
        }
    }
    memcpy(sorted, context_ptr->hme_results, sizeof(HmeResults)*MAX_NUM_OF_REF_PIC_LIST*REF_LIST_MAX_DEPTH);
    HmeResults     * res_p = sorted[0];
    uint32_t i, j;
    for (i = 0; i < num_of_cand_to_sort - 1; ++i) {
        for (j = i + 1; j < num_of_cand_to_sort; ++j) {
            if (res_p[j].hme_sad < res_p[i].hme_sad) {
                HmeResults temp = res_p[i];
                res_p[i] = res_p[j];
                res_p[j]= temp;
            }
        }
    }
    uint8_t  BIGGER_THAN_TH = 30;
    uint64_t best = sorted[0][0].hme_sad;//is this always the best?
    for (uint32_t li = 0; li < MAX_NUM_OF_REF_PIC_LIST; li++) {
        for (uint32_t ri = 0; ri < REF_LIST_MAX_DEPTH; ri++){
           // uint32_t dev = ((context_ptr->hme_results[li][ri].hme_sad - best) * 100) / best;
            if ((context_ptr->hme_results[li][ri].hme_sad - best) * 100  > BIGGER_THAN_TH*best)
                context_ptr->hme_results[li][ri].do_ref = 0;
            if (context_ptr->half_pel_mode == SWITCHABLE_HP_MODE)
                if (context_ptr->hme_results[li][ri].hme_sad > sorted[0][1].hme_sad)
                    context_ptr->local_hp_mode[li][ri] = REFINEMENT_HP_MODE;
        }
    }
}

/*******************************************
 *   performs hierarchical ME for every ref frame
 *******************************************/
void hme_sb(
    PictureParentControlSet   *pcs_ptr,
    uint32_t                   sb_origin_x,
    uint32_t                   sb_origin_y,
    MeContext                 *context_ptr,
    EbPictureBufferDesc       *input_ptr
){

    SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
    uint32_t sb_width = (input_ptr->width - sb_origin_x) < BLOCK_SIZE_64
                            ? input_ptr->width - sb_origin_x
                            : BLOCK_SIZE_64;
    uint32_t sb_height = (input_ptr->height - sb_origin_y) < BLOCK_SIZE_64
                             ? input_ptr->height - sb_origin_y
                             : BLOCK_SIZE_64;

    int16_t origin_x = (int16_t)sb_origin_x;
    int16_t origin_y = (int16_t)sb_origin_y;
    uint32_t num_of_list_to_search;
    uint32_t list_index;
    uint8_t ref_pic_index;
    uint8_t num_of_ref_pic_to_search;
    EbPaReferenceObject *reference_object; // input parameter, reference Object Ptr
    // HME
    uint32_t search_region_number_in_width  = 0;
    uint32_t search_region_number_in_height = 0;
    int16_t  x_hme_level_0_search_center[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                       [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t y_hme_level_0_search_center[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                       [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t hme_level0_sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t  x_hme_level_1_search_center[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                       [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t y_hme_level_1_search_center[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                       [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t hme_level1_sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t  x_hme_level_2_search_center[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                       [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t y_hme_level_2_search_center[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                       [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t hme_level2_sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    // Final ME Search Center
    int16_t x_search_center = 0;
    int16_t y_search_center = 0;
    // Hierarchical ME Search Center
    int16_t x_hme_search_center = 0;
    int16_t y_hme_search_center = 0;
    // Search Center SADs
    uint64_t hme_mv_sad = 0;
    EbPictureBufferDesc *ref_pic_ptr;
    EbPictureBufferDesc *quarter_ref_pic_ptr;
    EbPictureBufferDesc *sixteenth_ref_pic_ptr;
    int16_t temp_x_hme_search_center = 0;
    int16_t temp_y_hme_search_center = 0;
    uint32_t num_quad_in_width;
    uint32_t total_me_quad;
    uint32_t quad_index;
    uint32_t next_quad_index;
    uint64_t temp_x_hme_sad;
    uint64_t ref_0_poc = 0;
    uint64_t ref_1_poc = 0;
    int16_t hme_level1_search_area_in_width;
    int16_t hme_level1_search_area_in_height;
    // Configure HME level 0, level 1 and level 2 from static config parameters
    EbBool enable_hme_level0_flag = context_ptr->enable_hme_level0_flag;
    EbBool enable_hme_level1_flag = context_ptr->enable_hme_level1_flag;
    EbBool enable_hme_level2_flag = context_ptr->enable_hme_level2_flag;
    uint64_t best_cost = (uint64_t)~0;
    context_ptr->best_list_idx = 0;
    context_ptr->best_ref_idx = 0;
    EbBool one_quadrant_hme      = EB_FALSE;
    one_quadrant_hme = scs_ptr->input_resolution < INPUT_SIZE_4K_RANGE ? 0 : one_quadrant_hme;
    num_of_list_to_search =
        (pcs_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
    if (context_ptr->me_alt_ref == EB_TRUE) num_of_list_to_search = 0;
    // Uni-Prediction motion estimation loop
    // List Loop
    for (list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
        if (context_ptr->me_alt_ref == EB_TRUE)
            num_of_ref_pic_to_search = 1;
        else {
            num_of_ref_pic_to_search = (pcs_ptr->slice_type == P_SLICE)
                ? pcs_ptr->ref_list0_count_try
                : (list_index == REF_LIST_0) ? pcs_ptr->ref_list0_count_try
                : pcs_ptr->ref_list1_count_try;
            reference_object =
                (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[0][0]->object_ptr;
            ref_0_poc = pcs_ptr->ref_pic_poc_array[0][0];
        }
        // Ref Picture Loop
        for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            if (context_ptr->me_alt_ref == EB_TRUE)
                reference_object = (EbPaReferenceObject *)context_ptr->alt_ref_reference_ptr;
            else {
                if (num_of_list_to_search) {
                    reference_object =
                        (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[1][0]->object_ptr;
                    ref_1_poc = pcs_ptr->ref_pic_poc_array[1][0];
                }
                reference_object =
                    (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]
                        ->object_ptr;
            }
            ref_pic_ptr = (EbPictureBufferDesc *)reference_object->input_padded_picture_ptr;
            // Set 1/4 and 1/16 ME reference buffer(s); filtered or decimated
            quarter_ref_pic_ptr =
                (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED)
                    ? (EbPictureBufferDesc *)reference_object->quarter_filtered_picture_ptr
                    : (EbPictureBufferDesc *)reference_object->quarter_decimated_picture_ptr;

            sixteenth_ref_pic_ptr =
                (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED)
                    ? (EbPictureBufferDesc *)reference_object->sixteenth_filtered_picture_ptr
                    : (EbPictureBufferDesc *)reference_object->sixteenth_decimated_picture_ptr;

            EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_picture_ptr;

            // Use scaled source references if resolution of the reference is different that of the input
            use_scaled_source_refs_if_needed(pcs_ptr,
                                             input_picture_ptr,
                                             reference_object,
                                             &ref_pic_ptr,
                                             &quarter_ref_pic_ptr,
                                             &sixteenth_ref_pic_ptr);

            if (pcs_ptr->temporal_layer_index > 0 || list_index == 0) {
                if (context_ptr->update_hme_search_center_flag)
                    hme_mv_center_check(ref_pic_ptr,
                                        context_ptr,
                                        &x_search_center,
                                        &y_search_center,
                                        list_index,
                                        origin_x,
                                        origin_y,
                                        sb_width,
                                        sb_height);
                else {
                    x_search_center = 0;
                    y_search_center = 0;
                }
                if (context_ptr->enable_hme_flag && sb_height == BLOCK_SIZE_64){
                    while (search_region_number_in_height <
                           context_ptr->number_hme_search_region_in_height){
                        while (search_region_number_in_width <
                               context_ptr->number_hme_search_region_in_width){
                            x_hme_level_0_search_center[search_region_number_in_width]
                                                       [search_region_number_in_height] =
                                                           x_search_center;
                            y_hme_level_0_search_center[search_region_number_in_width]
                                                       [search_region_number_in_height] =
                                                           y_search_center;

                            x_hme_level_1_search_center[search_region_number_in_width]
                                                       [search_region_number_in_height] =
                                                           x_search_center;
                            y_hme_level_1_search_center[search_region_number_in_width]
                                                       [search_region_number_in_height] =
                                                           y_search_center;

                            x_hme_level_2_search_center[search_region_number_in_width]
                                                       [search_region_number_in_height] =
                                                           x_search_center;
                            y_hme_level_2_search_center[search_region_number_in_width]
                                                       [search_region_number_in_height] =
                                                           y_search_center;

                            search_region_number_in_width++;
                        }
                        search_region_number_in_width = 0;
                        search_region_number_in_height++;
                    }
                    // HME: Level0 search
                    if (enable_hme_level0_flag) {
                        if (one_quadrant_hme && !enable_hme_level1_flag &&
                            !enable_hme_level2_flag) {
                            search_region_number_in_height = 0;
                            search_region_number_in_width  = 0;
                            hme_one_quadrant_level_0(
                                pcs_ptr,
                                context_ptr,
                                origin_x >> 2,
                                origin_y >> 2,
                                sb_width >> 2,
                                sb_height >> 2,
                                x_search_center >> 2,
                                y_search_center >> 2,
                                sixteenth_ref_pic_ptr,
                                &(hme_level0_sad[search_region_number_in_width]
                                               [search_region_number_in_height]),
                                &(x_hme_level_0_search_center[search_region_number_in_width]
                                                             [search_region_number_in_height]),
                                &(y_hme_level_0_search_center[search_region_number_in_width]
                                                             [search_region_number_in_height]),
                                hme_level_0_search_area_multiplier_x[pcs_ptr->hierarchical_levels]
                                                                    [pcs_ptr->temporal_layer_index],
                                hme_level_0_search_area_multiplier_y
                                    [pcs_ptr->hierarchical_levels][pcs_ptr->temporal_layer_index]);
                        } else {
                            search_region_number_in_height = 0;
                            search_region_number_in_width  = 0;
                            while (search_region_number_in_height <
                                context_ptr->number_hme_search_region_in_height) {
                                while (search_region_number_in_width <
                                    context_ptr->number_hme_search_region_in_width) {
                                    hme_level_0(
                                        pcs_ptr,
                                        context_ptr,
                                        origin_x >> 2,
                                        origin_y >> 2,
                                        sb_width >> 2,
                                        sb_height >> 2,
                                        x_search_center >> 2,
                                        y_search_center >> 2,
                                        sixteenth_ref_pic_ptr,
                                        search_region_number_in_width,
                                        search_region_number_in_height,
                                        &(hme_level0_sad[search_region_number_in_width]
                                        [search_region_number_in_height]),
                                        &(x_hme_level_0_search_center
                                        [search_region_number_in_width]
                                        [search_region_number_in_height]),
                                        &(y_hme_level_0_search_center
                                        [search_region_number_in_width]
                                        [search_region_number_in_height]),
                                        hme_level_0_search_area_multiplier_x
                                        [pcs_ptr->hierarchical_levels]
                                        [pcs_ptr->temporal_layer_index],
                                        hme_level_0_search_area_multiplier_y
                                        [pcs_ptr->hierarchical_levels]
                                        [pcs_ptr->temporal_layer_index]);
                                    search_region_number_in_width++;
                                }
                                search_region_number_in_width = 0;
                                search_region_number_in_height++;
                            }
                        }
                    }
                    // HME: Level1 search
                    if (enable_hme_level1_flag) {
                        search_region_number_in_height = 0;
                        search_region_number_in_width  = 0;
                        while (search_region_number_in_height <
                            context_ptr->number_hme_search_region_in_height) {
                            while (search_region_number_in_width <
                                context_ptr->number_hme_search_region_in_width) {
                                // When HME level 0 has been disabled,
                                // increase the search area width and height
                                hme_level1_search_area_in_width =
                                    (int16_t)context_ptr->hme_level1_search_area_in_width_array
                                    [search_region_number_in_width];
                                hme_level1_search_area_in_height =
                                    (int16_t)context_ptr->hme_level1_search_area_in_height_array
                                    [search_region_number_in_height];
                                hme_level_1(context_ptr,
                                    origin_x >> 1,
                                    origin_y >> 1,
                                    sb_width >> 1,
                                    sb_height >> 1,
                                    quarter_ref_pic_ptr,
                                    hme_level1_search_area_in_width,
                                    hme_level1_search_area_in_height,
                                    x_hme_level_0_search_center
                                    [search_region_number_in_width]
                                    [search_region_number_in_height] >>
                                        1,
                                    y_hme_level_0_search_center
                                    [search_region_number_in_width]
                                    [search_region_number_in_height] >>
                                        1,
                                    &(hme_level1_sad[search_region_number_in_width]
                                    [search_region_number_in_height]),
                                    &(x_hme_level_1_search_center
                                    [search_region_number_in_width]
                                    [search_region_number_in_height]),
                                    &(y_hme_level_1_search_center
                                    [search_region_number_in_width]
                                    [search_region_number_in_height]));
                                    search_region_number_in_width++;
                            }
                            search_region_number_in_width = 0;
                            search_region_number_in_height++;
                        }
                    }
                    // HME: Level2 search
                    if (enable_hme_level2_flag) {
                        search_region_number_in_height = 0;
                        search_region_number_in_width  = 0;
                        {
                            while (search_region_number_in_height <
                                   context_ptr->number_hme_search_region_in_height) {
                                while (search_region_number_in_width <
                                       context_ptr->number_hme_search_region_in_width) {
                                    hme_level_2(
                                        pcs_ptr,
                                        context_ptr,
                                        origin_x,
                                        origin_y,
                                        sb_width,
                                        sb_height,
                                        ref_pic_ptr,
                                        search_region_number_in_width,
                                        search_region_number_in_height,
                                        x_hme_level_1_search_center[search_region_number_in_width]
                                                                   [search_region_number_in_height],
                                        y_hme_level_1_search_center[search_region_number_in_width]
                                                                   [search_region_number_in_height],
                                        &(hme_level2_sad[search_region_number_in_width]
                                                       [search_region_number_in_height]),
                                        &(x_hme_level_2_search_center
                                              [search_region_number_in_width]
                                              [search_region_number_in_height]),
                                        &(y_hme_level_2_search_center
                                              [search_region_number_in_width]
                                              [search_region_number_in_height]));

                                    search_region_number_in_width++;
                                }
                                search_region_number_in_width = 0;
                                search_region_number_in_height++;
                            }
                        }
                    }
                    // Hierarchical ME - Search Center
                    if (enable_hme_level0_flag && !enable_hme_level1_flag &&
                        !enable_hme_level2_flag) {
                        if (one_quadrant_hme) {
                            x_hme_search_center = x_hme_level_0_search_center[0][0];
                            y_hme_search_center = y_hme_level_0_search_center[0][0];
                            hme_mv_sad           = hme_level0_sad[0][0];
                        } else {
                            x_hme_search_center = x_hme_level_0_search_center[0][0];
                            y_hme_search_center = y_hme_level_0_search_center[0][0];
                            hme_mv_sad           = hme_level0_sad[0][0];
                            search_region_number_in_width  = 1;
                            search_region_number_in_height = 0;
                            while (search_region_number_in_height <
                                   context_ptr->number_hme_search_region_in_height) {
                                while (search_region_number_in_width <
                                       context_ptr->number_hme_search_region_in_width) {
                                    x_hme_search_center =
                                        (hme_level0_sad[search_region_number_in_width]
                                                      [search_region_number_in_height] < hme_mv_sad)
                                            ? x_hme_level_0_search_center
                                                  [search_region_number_in_width]
                                                  [search_region_number_in_height]
                                            : x_hme_search_center;
                                    y_hme_search_center =
                                        (hme_level0_sad[search_region_number_in_width]
                                                      [search_region_number_in_height] < hme_mv_sad)
                                            ? y_hme_level_0_search_center
                                                  [search_region_number_in_width]
                                                  [search_region_number_in_height]
                                            : y_hme_search_center;
                                    hme_mv_sad =
                                        (hme_level0_sad[search_region_number_in_width]
                                                      [search_region_number_in_height] < hme_mv_sad)
                                            ? hme_level0_sad[search_region_number_in_width]
                                                           [search_region_number_in_height]
                                            : hme_mv_sad;
                                    search_region_number_in_width++;
                                }
                                search_region_number_in_width = 0;
                                search_region_number_in_height++;
                            }
                        }
                    }
                    if (enable_hme_level1_flag && !enable_hme_level2_flag) {
                        x_hme_search_center = x_hme_level_1_search_center[0][0];
                        y_hme_search_center = y_hme_level_1_search_center[0][0];
                        hme_mv_sad           = hme_level1_sad[0][0];
                        search_region_number_in_width  = 1;
                        search_region_number_in_height = 0;
                        while (search_region_number_in_height <
                               context_ptr->number_hme_search_region_in_height) {
                            while (search_region_number_in_width <
                                   context_ptr->number_hme_search_region_in_width) {
                                x_hme_search_center =
                                    (hme_level1_sad[search_region_number_in_width]
                                                  [search_region_number_in_height] < hme_mv_sad)
                                        ? x_hme_level_1_search_center
                                              [search_region_number_in_width]
                                              [search_region_number_in_height]
                                        : x_hme_search_center;
                                y_hme_search_center =
                                    (hme_level1_sad[search_region_number_in_width]
                                                  [search_region_number_in_height] < hme_mv_sad)
                                        ? y_hme_level_1_search_center
                                              [search_region_number_in_width]
                                              [search_region_number_in_height]
                                        : y_hme_search_center;
                                hme_mv_sad =
                                    (hme_level1_sad[search_region_number_in_width]
                                                  [search_region_number_in_height] < hme_mv_sad)
                                        ? hme_level1_sad[search_region_number_in_width]
                                                       [search_region_number_in_height]
                                        : hme_mv_sad;
                                search_region_number_in_width++;
                            }
                            search_region_number_in_width = 0;
                            search_region_number_in_height++;
                        }
                    }
                    if (enable_hme_level2_flag) {
                        x_hme_search_center = x_hme_level_2_search_center[0][0];
                        y_hme_search_center = y_hme_level_2_search_center[0][0];
                        hme_mv_sad           = hme_level2_sad[0][0];
                        search_region_number_in_width  = 1;
                        search_region_number_in_height = 0;
                        while (search_region_number_in_height <
                               context_ptr->number_hme_search_region_in_height) {
                            while (search_region_number_in_width <
                                   context_ptr->number_hme_search_region_in_width) {
                                x_hme_search_center =
                                    (hme_level2_sad[search_region_number_in_width]
                                                  [search_region_number_in_height] < hme_mv_sad)
                                        ? x_hme_level_2_search_center
                                              [search_region_number_in_width]
                                              [search_region_number_in_height]
                                        : x_hme_search_center;
                                y_hme_search_center =
                                    (hme_level2_sad[search_region_number_in_width]
                                                  [search_region_number_in_height] < hme_mv_sad)
                                        ? y_hme_level_2_search_center
                                              [search_region_number_in_width]
                                              [search_region_number_in_height]
                                        : y_hme_search_center;
                                hme_mv_sad =
                                    (hme_level2_sad[search_region_number_in_width]
                                                  [search_region_number_in_height] < hme_mv_sad)
                                        ? hme_level2_sad[search_region_number_in_width]
                                                       [search_region_number_in_height]
                                        : hme_mv_sad;
                                search_region_number_in_width++;
                            }
                            search_region_number_in_width = 0;
                            search_region_number_in_height++;
                        }

                        num_quad_in_width = context_ptr->number_hme_search_region_in_width;
                        total_me_quad     = context_ptr->number_hme_search_region_in_height *
                                        context_ptr->number_hme_search_region_in_width;
                        if ((ref_0_poc == ref_1_poc) && (list_index == 1) && (total_me_quad > 1)) {
                            for (quad_index = 0; quad_index < total_me_quad - 1; ++quad_index) {
                                for (next_quad_index = quad_index + 1;
                                     next_quad_index < total_me_quad;
                                     ++next_quad_index) {
                                    if (hme_level2_sad[quad_index / num_quad_in_width]
                                                     [quad_index % num_quad_in_width] >
                                        hme_level2_sad[next_quad_index / num_quad_in_width]
                                                     [next_quad_index % num_quad_in_width]) {
                                        temp_x_hme_search_center =
                                            x_hme_level_2_search_center[quad_index /
                                                                        num_quad_in_width]
                                                                       [quad_index %
                                                                        num_quad_in_width];
                                        temp_y_hme_search_center =
                                            y_hme_level_2_search_center[quad_index /
                                                                        num_quad_in_width]
                                                                       [quad_index %
                                                                        num_quad_in_width];
                                        temp_x_hme_sad =
                                            hme_level2_sad[quad_index / num_quad_in_width]
                                                         [quad_index % num_quad_in_width];

                                        x_hme_level_2_search_center
                                            [quad_index / num_quad_in_width]
                                            [quad_index % num_quad_in_width] =
                                                x_hme_level_2_search_center[next_quad_index /
                                                                            num_quad_in_width]
                                                                           [next_quad_index %
                                                                            num_quad_in_width];
                                        y_hme_level_2_search_center
                                            [quad_index / num_quad_in_width]
                                            [quad_index % num_quad_in_width] =
                                                y_hme_level_2_search_center[next_quad_index /
                                                                            num_quad_in_width]
                                                                           [next_quad_index %
                                                                            num_quad_in_width];
                                        hme_level2_sad[quad_index /
                                                      num_quad_in_width][quad_index %
                                                                         num_quad_in_width] =
                                            hme_level2_sad[next_quad_index / num_quad_in_width]
                                                         [next_quad_index % num_quad_in_width];

                                        x_hme_level_2_search_center[next_quad_index /
                                                                    num_quad_in_width]
                                                                   [next_quad_index %
                                                                    num_quad_in_width] =
                                                                       temp_x_hme_search_center;
                                        y_hme_level_2_search_center[next_quad_index /
                                                                    num_quad_in_width]
                                                                   [next_quad_index %
                                                                    num_quad_in_width] =
                                                                       temp_y_hme_search_center;
                                        hme_level2_sad[next_quad_index / num_quad_in_width]
                                                     [next_quad_index % num_quad_in_width] =
                                                         temp_x_hme_sad;
                                    }
                                }
                            }
                            x_hme_search_center = x_hme_level_2_search_center[0][1];
                            y_hme_search_center = y_hme_level_2_search_center[0][1];
                        }
                    }
                    x_search_center = x_hme_search_center;
                    y_search_center = y_hme_search_center;
                }
            }else {
                x_search_center = 0;
                y_search_center = 0;
            }
            //sc valid for all cases. 0,0 if hme not done.
            context_ptr->hme_results[list_index][ref_pic_index].hme_sc_x = x_search_center;
            context_ptr->hme_results[list_index][ref_pic_index].hme_sc_y = y_search_center;
            context_ptr->hme_results[list_index][ref_pic_index].hme_sad = hme_mv_sad;//this is not valid in all cases. only when HME is done, and when HMELevel2 is done
            //also for base layer some references are redundant!!
            context_ptr->hme_results[list_index][ref_pic_index].do_ref = 1;
            if (hme_mv_sad < best_cost) {
                best_cost = hme_mv_sad;
                context_ptr->best_list_idx = list_index;
                context_ptr->best_ref_idx = ref_pic_index;
            }
        }
    }
}
void prune_references_sc(
    MeContext *context_ptr)
{
    for (uint32_t li = 0; li < MAX_NUM_OF_REF_PIC_LIST; li++) {
        for (uint32_t ri = 0; ri < REF_LIST_MAX_DEPTH; ri++){
            if (context_ptr->hme_results[li][ri].hme_sc_x == 0 &&
                context_ptr->hme_results[li][ri].hme_sc_y == 0 &&
                context_ptr->hme_results[li][ri].hme_sad < SC_HME_TH_EASY)

                context_ptr->reduce_me_sr_flag[li][ri] = SC_HME_TH_STILL;

            else if (context_ptr->hme_results[li][ri].hme_sad < SC_HME_TH_EASY)
                context_ptr->reduce_me_sr_flag[li][ri] = SC_HME_TH_EASY;
        }
    }
}
void prune_references(
    MeContext                 *context_ptr)
{
    HmeResults    sorted[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    uint32_t      num_of_cand_to_sort = MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH;
    memcpy(sorted, context_ptr->hme_results, sizeof(HmeResults)*MAX_NUM_OF_REF_PIC_LIST*REF_LIST_MAX_DEPTH);
    HmeResults     * res_p = sorted[0];
    uint32_t i, j;
    for (i = 0; i < num_of_cand_to_sort - 1; ++i) {
        for (j = i + 1; j < num_of_cand_to_sort; ++j) {
            if (res_p[j].hme_sad < res_p[i].hme_sad) {
                HmeResults temp = res_p[i];
                res_p[i] = res_p[j];
                res_p[j]= temp;
            }
        }
    }
    uint8_t  BIGGER_THAN_TH = 80;
    uint64_t best = sorted[0][0].hme_sad;//is this always the best?
    uint64_t REDUCE_SR_TH = 6000;
    int16_t  displacement_th = 4;
    for (uint32_t li = 0; li < MAX_NUM_OF_REF_PIC_LIST; li++) {
        for (uint32_t ri = 0; ri < REF_LIST_MAX_DEPTH; ri++){
            if ((context_ptr->hme_results[li][ri].hme_sad - best) * 100  > BIGGER_THAN_TH*best)
                context_ptr->hme_results[li][ri].do_ref = 0;
            if (context_ptr->hme_results[li][ri].hme_sad < REDUCE_SR_TH)
                context_ptr->reduce_me_sr_flag[li][ri] = 1;
            if (context_ptr->hme_results[li][ri].hme_sc_x <= displacement_th && context_ptr->hme_results[li][ri].hme_sc_y <= displacement_th && context_ptr->hme_results[li][ri].hme_sad < (2*REDUCE_SR_TH))
                context_ptr->reduce_me_sr_flag[li][ri] = 1;
        }
    }
}
/*******************************************
 * motion_estimate_sb
 *   performs ME (SB)
 *******************************************/
EbErrorType motion_estimate_sb(
    PictureParentControlSet *pcs_ptr, // input parameter, Picture Control Set Ptr
    uint32_t                 sb_index, // input parameter, SB Index
    uint32_t                 sb_origin_x, // input parameter, SB Origin X
    uint32_t                 sb_origin_y, // input parameter, SB Origin X
    MeContext               *context_ptr, // input parameter, ME Context Ptr, used to store decimated/interpolated SB/SR
    EbPictureBufferDesc     *input_ptr) // input parameter, source Picture Ptr

{
    EbErrorType return_error = EB_ErrorNone;

    SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
    uint32_t sb_height = (input_ptr->height - sb_origin_y) < BLOCK_SIZE_64
                             ? input_ptr->height - sb_origin_y
                             : BLOCK_SIZE_64;
    uint32_t pu_index;

    uint32_t max_number_of_pus_per_sb = pcs_ptr->max_number_of_pus_per_sb;

    uint32_t             num_of_list_to_search;
    uint32_t             list_index;
    uint8_t              cand_index               = 0;
    uint8_t              total_me_candidate_index = 0;
    EbPaReferenceObject *reference_object; // input parameter, reference Object Ptr

    uint8_t  ref_pic_index;
    uint8_t  num_of_ref_pic_to_search;
    uint8_t  candidate_index      = 0;
    uint32_t next_candidate_index = 0;

    MePredUnit *         me_candidate;
    EbPictureBufferDesc *ref_pic_ptr;
    uint64_t i;
    EbBool enable_half_pel_32x32 = EB_FALSE;
    EbBool enable_half_pel_16x16 = EB_FALSE;
    EbBool enable_half_pel_8x8   = EB_FALSE;
    EbBool enable_quarter_pel    = EB_FALSE;
    EbBool one_quadrant_hme      = EB_FALSE;

    one_quadrant_hme = scs_ptr->input_resolution < INPUT_SIZE_4K_RANGE ? 0 : one_quadrant_hme;

    num_of_list_to_search =
        (pcs_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
    //pruning of the references is not done for alt-ref / Base-Layer (HME not done for list1 refs) / non-complete-SBs when HMeLevel2 is done
    uint8_t prune_ref = (context_ptr->enable_hme_flag && context_ptr->enable_hme_level2_flag &&
        context_ptr->me_alt_ref == EB_FALSE && sb_height == BLOCK_SIZE_64 &&
        pcs_ptr->temporal_layer_index > 0) ? 1 : 0;
    //init hme results buffer
    for (uint32_t li = 0; li < MAX_NUM_OF_REF_PIC_LIST; li++) {
        for (uint32_t ri = 0; ri < REF_LIST_MAX_DEPTH; ri++) {
            context_ptr->hme_results[li][ri].list_i = li;
            context_ptr->hme_results[li][ri].ref_i = ri;
            context_ptr->hme_results[li][ri].do_ref = 1;
            context_ptr->hme_results[li][ri].hme_sad = 0xFFFFFFFF;
            context_ptr->reduce_me_sr_flag[li][ri] = 0;
            context_ptr->local_hp_mode[li][ri] =
                context_ptr->half_pel_mode == SWITCHABLE_HP_MODE ? EX_HP_MODE :
                context_ptr->half_pel_mode;
            // R2R FIX: no winner integer MV is set in special case like initial p_sb_best_mv for overlay case,
            // then it sends dirty p_sb_best_mv to MD, initializing it is necessary
            for(uint32_t pi = 0; pi < MAX_ME_PU_COUNT; pi++)
                context_ptr->p_sb_best_mv[li][ri][pi] = 0;
        }
    }
    // HME: Perform Hierachical Motion Estimation for all refrence frames.
    hme_sb(
        pcs_ptr,
        sb_origin_x,
        sb_origin_y,
        context_ptr,
        input_ptr);
    // prune the refrence frames based on the HME outputs.
    if (pcs_ptr->prune_ref_based_me && prune_ref)
        prune_references(
            context_ptr);
    else if (pcs_ptr->sc_content_detected && prune_ref)
        prune_references_sc(
            context_ptr);
    // Full pel: Perform the Integer Motion Estimation on the allowed refrence frames.
    integer_search_sb(
        pcs_ptr,
        sb_index,
        sb_origin_x,
        sb_origin_y,
        context_ptr,
        input_ptr);
    // prune the refrence frames based on the Full pel outputs.
    if (pcs_ptr->prune_ref_based_me && prune_ref)
        prune_references_fp(
            pcs_ptr,
            context_ptr );

    if (context_ptr->me_alt_ref == EB_TRUE) num_of_list_to_search = 0;

    // Uni-Prediction motion estimation loop
    // List Loop
    for (list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
        if (context_ptr->me_alt_ref == EB_TRUE) {
            num_of_ref_pic_to_search = 1;
        } else {
            num_of_ref_pic_to_search = (pcs_ptr->slice_type == P_SLICE)
                ? pcs_ptr->ref_list0_count_try
                : (list_index == REF_LIST_0) ? pcs_ptr->ref_list0_count_try
                : pcs_ptr->ref_list1_count_try;

            reference_object =
                (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[0][0]->object_ptr;
        }

        // Ref Picture Loop
        for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            if (context_ptr->me_alt_ref == EB_TRUE) {
                reference_object = (EbPaReferenceObject *)context_ptr->alt_ref_reference_ptr;
            } else {
                if (num_of_list_to_search) {
                    reference_object =
                        (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[1][0]->object_ptr;
                }

                reference_object =
                    (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]
                        ->object_ptr;
            }

            ref_pic_ptr = (EbPictureBufferDesc *)reference_object->input_padded_picture_ptr;
            if (ref_pic_ptr == NULL)
                printf("ERR NULL POINTER");
            if (context_ptr->hme_results[list_index][ref_pic_index].do_ref == 0)
                continue;  //so will not get ME results for those references. what will happen next, shall we just fill in max sads?
                           //we can also make the ME small and shut subpel
            {
                {
                    if (pcs_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {

                        context_ptr->p_best_sad_64x64 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_sad_32x32 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_sad_16x16 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_sad_8x8 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_sad_64x32 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_sad_32x16 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_sad_16x8 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_sad_32x64 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_sad_16x32 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_sad_8x16 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_sad_32x8 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_sad_8x32 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_sad_64x16 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_sad_16x64 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x64_0]);

                        context_ptr->p_best_mv64x64 =
                            &(context_ptr
                                  ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_mv32x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_mv16x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_mv8x8 =
                            &(context_ptr
                                  ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_mv64x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_mv32x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_mv16x8 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_mv32x64 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_mv16x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_mv8x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_mv32x8 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_mv8x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_mv64x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_mv16x64 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x64_0]);

                        context_ptr->p_best_ssd64x64 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_ssd32x32 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_ssd16x16 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_ssd8x8 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_ssd64x32 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_ssd32x16 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_ssd16x8 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_ssd32x64 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_ssd16x32 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_ssd8x16 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_ssd32x8 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_ssd8x32 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_ssd64x16 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_ssd16x64 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x64_0]);
                        context_ptr->full_quarter_pel_refinement = 0;
                        if (context_ptr->half_pel_mode == EX_HP_MODE ||
                            context_ptr->local_hp_mode[list_index][ref_pic_index] == EX_HP_MODE) {
                            // Move to the top left of the search region
                            //x_top_left_search_region =
                            //    (int16_t)(ref_pic_ptr->origin_x + sb_origin_x) +
                            //    context_ptr->x_search_area_origin[list_index][ref_pic_index];
                            //y_top_left_search_region =
                            //    (int16_t)(ref_pic_ptr->origin_y + sb_origin_y) +
                            //    context_ptr->y_search_area_origin[list_index][ref_pic_index];
                            // Interpolate the search region for Half-Pel
                            // Refinements H - AVC Style
                            interpolate_search_region_avc(
                                context_ptr,
                                list_index,
                                ref_pic_index,
                                context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
                                    (ME_FILTER_TAP >> 1) +
                                    ((ME_FILTER_TAP >> 1) *
                                     context_ptr
                                         ->interpolated_full_stride[list_index][ref_pic_index]),
                                context_ptr->interpolated_full_stride[list_index][ref_pic_index],
                                MAX(1, (uint32_t)context_ptr->sa_width[list_index][ref_pic_index]) + (BLOCK_SIZE_64 - 1),
                                MAX(1, (uint32_t)context_ptr->sa_height[list_index][ref_pic_index]) + (BLOCK_SIZE_64 - 1),
                                8);

                            initialize_buffer_32bits(
                                context_ptr->p_sb_best_ssd[list_index][ref_pic_index],
                                52,
                                1,
                                MAX_SSE_VALUE);
                            memcpy(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index],
                                   context_ptr->p_sb_best_mv[list_index][ref_pic_index],
                                   MAX_ME_PU_COUNT * sizeof(uint32_t));
                            context_ptr->full_quarter_pel_refinement = 1;
                            context_ptr->p_best_full_pel_mv64x64 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_64x64]);
                            context_ptr->p_best_full_pel_mv32x32 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_32x32_0]);
                            context_ptr->p_best_full_pel_mv16x16 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_16x16_0]);
                            context_ptr->p_best_full_pel_mv8x8 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_8x8_0]);
                            context_ptr->p_best_full_pel_mv64x32 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_64x32_0]);
                            context_ptr->p_best_full_pel_mv32x16 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_32x16_0]);
                            context_ptr->p_best_full_pel_mv16x8 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_16x8_0]);
                            context_ptr->p_best_full_pel_mv32x64 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_32x64_0]);
                            context_ptr->p_best_full_pel_mv16x32 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_16x32_0]);
                            context_ptr->p_best_full_pel_mv8x16 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_8x16_0]);
                            context_ptr->p_best_full_pel_mv32x8 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_32x8_0]);
                            context_ptr->p_best_full_pel_mv8x32 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_8x32_0]);
                            context_ptr->p_best_full_pel_mv64x16 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_64x16_0]);
                            context_ptr->p_best_full_pel_mv16x64 =
                                &(context_ptr->p_sb_best_full_pel_mv[list_index][ref_pic_index]
                                                                    [ME_TIER_ZERO_PU_16x64_0]);
                            // half-Pel search
                            open_loop_me_half_pel_search_sblock(pcs_ptr,
                                                                context_ptr,
                                                                list_index,
                                                                ref_pic_index,
                                                                context_ptr->x_search_area_origin[list_index][ref_pic_index],
                                                                context_ptr->y_search_area_origin[list_index][ref_pic_index],
                                                                context_ptr->sa_width[list_index][ref_pic_index],
                                                                context_ptr->sa_height[list_index][ref_pic_index]);
                        }


                    } else {
                        initialize_buffer_32bits(
                            context_ptr->p_sb_best_sad[list_index][ref_pic_index],
                            21,
                            1,
                            MAX_SAD_VALUE);

                        context_ptr->full_quarter_pel_refinement = 0;

                        context_ptr->p_best_sad_64x64 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_sad_32x32 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_sad_16x16 =
                            &(context_ptr->p_sb_best_sad[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_sad_8x8 = &(
                            context_ptr
                                ->p_sb_best_sad[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);

                        context_ptr->p_best_mv64x64 =
                            &(context_ptr
                                  ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_mv32x32 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_mv16x16 = &(
                            context_ptr
                                ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_mv8x8 =
                            &(context_ptr
                                  ->p_sb_best_mv[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);

                        context_ptr->p_best_ssd64x64 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_ssd32x32 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_ssd16x16 =
                            &(context_ptr->p_sb_best_ssd[list_index][ref_pic_index]
                                                        [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_ssd8x8 = &(
                            context_ptr
                                ->p_sb_best_ssd[list_index][ref_pic_index][ME_TIER_ZERO_PU_8x8_0]);
                    }
                }

                if (context_ptr->fractional_search_model == 0) {
                    enable_half_pel_32x32 = EB_TRUE;
                    enable_half_pel_16x16 = EB_TRUE;
                    enable_half_pel_8x8   = EB_TRUE;
                    enable_quarter_pel    = EB_TRUE;
                } else if (context_ptr->fractional_search_model == 1) {
                    su_pel_enable(context_ptr,
                                  pcs_ptr,
                                  list_index,
                                  0,
                                  &enable_half_pel_32x32,
                                  &enable_half_pel_16x16,
                                  &enable_half_pel_8x8);
                    enable_quarter_pel = EB_TRUE;
                } else {
                    enable_half_pel_32x32 = EB_FALSE;
                    enable_half_pel_16x16 = EB_FALSE;
                    enable_half_pel_8x8   = EB_FALSE;
                    enable_quarter_pel    = EB_FALSE;
                }
                if (enable_half_pel_32x32 || enable_half_pel_16x16 || enable_half_pel_8x8 ||
                    enable_quarter_pel) {
                    // if((pcs_ptr->is_used_as_reference_flag ==
                    // EB_TRUE)) {
                    // Move to the top left of the search region
                    //x_top_left_search_region =
                    //    (int16_t)(ref_pic_ptr->origin_x + sb_origin_x) + context_ptr->x_search_area_origin[list_index][ref_pic_index];
                    //y_top_left_search_region =
                    //    (int16_t)(ref_pic_ptr->origin_y + sb_origin_y) + context_ptr->y_search_area_origin[list_index][ref_pic_index];

                    // Interpolate the search region for Half-Pel Refinements
                    // H - AVC Style
                    if (context_ptr->half_pel_mode == REFINEMENT_HP_MODE ||
                        context_ptr->local_hp_mode[list_index][ref_pic_index] == REFINEMENT_HP_MODE) {
                        interpolate_search_region_avc(
                            context_ptr,
                            list_index,
                            ref_pic_index,
                            context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
                                (ME_FILTER_TAP >> 1) +
                                ((ME_FILTER_TAP >> 1) *
                                 context_ptr->interpolated_full_stride[list_index][ref_pic_index]),
                            context_ptr->interpolated_full_stride[list_index][ref_pic_index],
                            MAX(1, (uint32_t)context_ptr->sa_width[list_index][ref_pic_index]) + (BLOCK_SIZE_64 - 1),
                            MAX(1, (uint32_t)context_ptr->sa_height[list_index][ref_pic_index]) + (BLOCK_SIZE_64 - 1),
                            8);

                        // Half-Pel Refinement [8 search positions]
                        half_pel_search_sb(
                            scs_ptr,
                            pcs_ptr,
                            context_ptr,
                            context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
                                (ME_FILTER_TAP >> 1) +
                                ((ME_FILTER_TAP >> 1) *
                                 context_ptr->interpolated_full_stride[list_index][ref_pic_index]),
                            context_ptr->interpolated_full_stride[list_index][ref_pic_index],
                            &(context_ptr->pos_b_buffer[list_index][ref_pic_index]
                                                       [(ME_FILTER_TAP >> 1) *
                                                        context_ptr->interpolated_stride]),
                            &(context_ptr->pos_h_buffer[list_index][ref_pic_index][1]),
                            &(context_ptr->pos_j_buffer[list_index][ref_pic_index][0]),
                            context_ptr->x_search_area_origin[list_index][ref_pic_index],
                            context_ptr->y_search_area_origin[list_index][ref_pic_index],
                            enable_half_pel_32x32,
                            enable_half_pel_16x16,
                            enable_half_pel_8x8);
                    }

                    {
                        // Quarter-Pel Refinement [8 search positions]
                        quarter_pel_search_sb(
                            context_ptr,
                            context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
                                (ME_FILTER_TAP >> 1) +
                                ((ME_FILTER_TAP >> 1) *
                                 context_ptr->interpolated_full_stride[list_index][ref_pic_index]),
                            context_ptr->interpolated_full_stride[list_index][ref_pic_index],
                            &(context_ptr
                                  ->pos_b_buffer[list_index][ref_pic_index]
                                                [(ME_FILTER_TAP >> 1) *
                                                 context_ptr->interpolated_stride]), // points to b
                            // position of
                            // the figure
                            // above

                            &(context_ptr->pos_h_buffer[list_index][ref_pic_index]
                                                       [1]), // points to h position
                            // of the figure above
                            &(context_ptr->pos_j_buffer[list_index][ref_pic_index]
                                                       [0]), // points to j position
                            // of the figure above
                            context_ptr->x_search_area_origin[list_index][ref_pic_index],
                            context_ptr->y_search_area_origin[list_index][ref_pic_index],
                            enable_half_pel_32x32,
                            enable_half_pel_16x16,
                            enable_half_pel_8x8,
                            enable_quarter_pel,
                            pcs_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE);
                    }
                }
            }
        }
    }

    if (context_ptr->me_alt_ref == EB_FALSE) {
        // Bi-Prediction motion estimation loop
        for (pu_index = 0; pu_index < max_number_of_pus_per_sb; ++pu_index) {
            cand_index = 0;

            uint32_t n_idx;

            if (pu_index > 200)
                n_idx = pu_index;
            else if (pu_index > 184)
                n_idx = tab8x32[pu_index - 185] + 185;
            else if (pu_index > 168)
                n_idx = tab32x8[pu_index - 169] + 169;
            else if (pu_index > 136)
                n_idx = tab8x16[pu_index - 137] + 137;
            else if (pu_index > 128)
                n_idx = tab16x32[pu_index - 129] + 129;
            else if (pu_index > 126)
                n_idx = pu_index;
            else if (pu_index > 94)
                n_idx = tab16x8[pu_index - 95] + 95;
            else if (pu_index > 86)
                n_idx = tab32x16[pu_index - 87] + 87;
            else if (pu_index > 84)
                n_idx = pu_index;
            else if (pu_index > 20)
                n_idx = tab8x8[pu_index - 21] + 21;
            else if (pu_index > 4)
                n_idx = tab16x16[pu_index - 5] + 5;
            else
                n_idx = pu_index;
            for (list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
                num_of_ref_pic_to_search = (pcs_ptr->slice_type == P_SLICE)
                    ? pcs_ptr->ref_list0_count_try
                    : (list_index == REF_LIST_0)
                    ? pcs_ptr->ref_list0_count_try
                    : pcs_ptr->ref_list1_count_try;

                // Ref Picture Loop
                for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
                    //ME was skipped, so do not add this Unipred candidate
                    if (context_ptr->hme_results[list_index][ref_pic_index].do_ref == 0)
                        continue;
                    me_candidate = &(context_ptr->me_candidate[cand_index].pu[pu_index]);
                    me_candidate->prediction_direction  = list_index;
                    me_candidate->ref_index[list_index] = ref_pic_index;
                    me_candidate->ref0_list =
                        me_candidate->prediction_direction == 0 ? list_index : 24;
                    me_candidate->ref1_list =
                        me_candidate->prediction_direction == 1 ? list_index : 24;
                    me_candidate->distortion =
                        context_ptr->p_sb_best_sad[list_index][ref_pic_index][n_idx];
                    cand_index++;
                }
            }

            total_me_candidate_index = cand_index;
            uint8_t ref_type_table[7];
            if (pcs_ptr->prune_unipred_at_me) {
                // Sorting of the ME candidates
                for (candidate_index = 0; candidate_index < total_me_candidate_index - 1;
                     ++candidate_index) {
                    for (next_candidate_index = candidate_index + 1;
                         next_candidate_index < total_me_candidate_index;
                         ++next_candidate_index) {
                        if (context_ptr->me_candidate[candidate_index].pu[pu_index].distortion >
                            context_ptr->me_candidate[next_candidate_index]
                                .pu[pu_index]
                                .distortion) {
                            swap_me_candidate(
                                &(context_ptr->me_candidate[candidate_index].pu[pu_index]),
                                &(context_ptr->me_candidate[next_candidate_index].pu[pu_index]));
                        }
                    }
                }
                for (candidate_index = 0; candidate_index < total_me_candidate_index;
                     ++candidate_index) {
                    me_candidate = &(context_ptr->me_candidate[candidate_index].pu[pu_index]);

                    if (me_candidate->prediction_direction == 0)
                        ref_type_table[candidate_index] = svt_get_ref_frame_type(
                            me_candidate->ref0_list, me_candidate->ref_index[0]);
                    else
                        ref_type_table[candidate_index] = svt_get_ref_frame_type(
                            me_candidate->ref1_list, me_candidate->ref_index[1]);
                }
            }
            if (num_of_list_to_search) {
                bi_prediction_search(scs_ptr,
                                        context_ptr,
                                        pu_index,
                                        cand_index,
                                        pcs_ptr->ref_list0_count_try,
                                        pcs_ptr->ref_list1_count_try,
                                        &total_me_candidate_index,
                                        ref_type_table,
                                        pcs_ptr);
            }

            // Sorting of the ME candidates
            for (candidate_index = 0; candidate_index < total_me_candidate_index - 1;
                 ++candidate_index) {
                for (next_candidate_index = candidate_index + 1;
                     next_candidate_index < total_me_candidate_index;
                     ++next_candidate_index) {
                    if (context_ptr->me_candidate[candidate_index].pu[pu_index].distortion >
                        context_ptr->me_candidate[next_candidate_index].pu[pu_index].distortion) {
                        swap_me_candidate(
                            &(context_ptr->me_candidate[candidate_index].pu[pu_index]),
                            &(context_ptr->me_candidate[next_candidate_index].pu[pu_index]));
                    }
                }
            }

            MeSbResults *me_pu_result                        = pcs_ptr->me_results[sb_index];
            me_pu_result->total_me_candidate_index[pu_index] = total_me_candidate_index;
            me_pu_result->total_me_candidate_index[pu_index] =
                MIN(total_me_candidate_index, ME_RES_CAND_MRP_MODE_0);
            // Assining the ME candidates to the me Results buffer
            for (cand_index = 0; cand_index < total_me_candidate_index; ++cand_index) {
                me_candidate = &(context_ptr->me_candidate[cand_index].pu[pu_index]);
                pcs_ptr->me_results[sb_index]->me_candidate[pu_index][cand_index].direction =
                    me_candidate->prediction_direction;
                pcs_ptr->me_results[sb_index]->me_candidate[pu_index][cand_index].ref_idx_l0 =
                    me_candidate->ref_index[0];
                pcs_ptr->me_results[sb_index]->me_candidate[pu_index][cand_index].ref_idx_l1 =
                    me_candidate->ref_index[1];
                pcs_ptr->me_results[sb_index]->me_candidate[pu_index][cand_index].ref0_list =
                    me_candidate->ref0_list;
                pcs_ptr->me_results[sb_index]->me_candidate[pu_index][cand_index].ref1_list =
                    me_candidate->ref1_list;
            }

            for (list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
                num_of_ref_pic_to_search = (pcs_ptr->slice_type == P_SLICE)
                    ? pcs_ptr->ref_list0_count_try
                    : (list_index == REF_LIST_0)
                    ? pcs_ptr->ref_list0_count_try
                    : pcs_ptr->ref_list1_count_try;

                // Ref Picture Loop
                for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
                    pcs_ptr->me_results[sb_index]
                        ->me_mv_array[pu_index][((list_index && scs_ptr->mrp_mode == 0)
                                                     ? 4
                                                     : list_index ? 2 : 0) +
                                                ref_pic_index]
                        .x_mv = _MVXT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][n_idx]);
                    pcs_ptr->me_results[sb_index]
                        ->me_mv_array[pu_index][((list_index && scs_ptr->mrp_mode == 0)
                                                     ? 4
                                                     : list_index ? 2 : 0) +
                                                ref_pic_index]
                        .y_mv = _MVYT(context_ptr->p_sb_best_mv[list_index][ref_pic_index][n_idx]);
                }
            }
        }
        {
            // Compute the sum of the distortion of all 16 16x16 (best) blocks
            // in the SB
            pcs_ptr->rc_me_distortion[sb_index] = 0;
            for (i = 0; i < 16; i++) {
                me_candidate = &(context_ptr->me_candidate[0].pu[5 + i]);
                pcs_ptr->rc_me_distortion[sb_index] +=
                    me_candidate->distortion;
            }
        }
    }

    return return_error;
}

EbErrorType open_loop_intra_search_sb(PictureParentControlSet *pcs_ptr, uint32_t sb_index,
                                      MotionEstimationContext_t *context_ptr,
                                      EbPictureBufferDesc *      input_ptr) {
    EbErrorType         return_error = EB_ErrorNone;

    uint32_t      blk_origin_x;
    uint32_t      blk_origin_y;
    uint32_t      pa_blk_index       = 0;
    SbParams *    sb_params          = &pcs_ptr->sb_params_array[sb_index];
    OisSbResults *ois_sb_results_ptr = pcs_ptr->ois_sb_results[sb_index];
    uint8_t *     above_row;
    uint8_t *     left_col;

    DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);
    while (pa_blk_index < CU_MAX_COUNT) {
        const CodedBlockStats *blk_stats_ptr;
        blk_stats_ptr = get_coded_blk_stats(pa_blk_index);
        uint8_t bsize = blk_stats_ptr->size;
        TxSize  tx_size =
            bsize == 8 ? TX_8X8 : bsize == 16 ? TX_16X16 : bsize == 32 ? TX_32X32 : TX_64X64;
        if (sb_params->raster_scan_blk_validity[md_scan_to_raster_scan[pa_blk_index]]) {
            OisCandidate *ois_blk_ptr = ois_sb_results_ptr->ois_candidate_array[pa_blk_index];
            blk_origin_x              = sb_params->origin_x + blk_stats_ptr->origin_x;
            blk_origin_y              = sb_params->origin_y + blk_stats_ptr->origin_y;
            above_row                 = above_data + 16;
            left_col                  = left_data + 16;

            // Fill Neighbor Arrays
            update_neighbor_samples_array_open_loop(above_row - 1,
                                                    left_col - 1,
                                                    input_ptr,
                                                    input_ptr->stride_y,
                                                    blk_origin_x,
                                                    blk_origin_y,
                                                    bsize,
                                                    bsize);
            uint8_t  ois_intra_mode;
            uint8_t  ois_intra_count             = 0;
            uint8_t  best_intra_ois_index        = 0;
            uint32_t best_intra_ois_distortion   = 64 * 64 * 255;
            uint8_t  intra_mode_start            = DC_PRED;
            EbBool   enable_paeth                = pcs_ptr->scs_ptr->static_config.enable_paeth == DEFAULT ? EB_TRUE : (EbBool) pcs_ptr->scs_ptr->static_config.enable_paeth;
            EbBool   enable_smooth               = pcs_ptr->scs_ptr->static_config.enable_smooth == DEFAULT ? EB_TRUE : (EbBool) pcs_ptr->scs_ptr->static_config.enable_smooth;
            uint8_t  intra_mode_end              = enable_paeth ? PAETH_PRED : enable_smooth ? SMOOTH_H_PRED : D67_PRED;

            uint8_t  angle_delta_counter         = 0;
            uint8_t  angle_delta_shift           = 1;
            EbBool   use_angle_delta             = av1_use_angle_delta(bsize, pcs_ptr->scs_ptr->static_config.intra_angle_delta);
            uint8_t  angle_delta_candidate_count = use_angle_delta ? 7 : 1;
            uint8_t  disable_angular_prediction  = 0;
            if (pcs_ptr->intra_pred_mode == 5) {
                intra_mode_end =
                    (pcs_ptr->is_used_as_reference_flag == 0) ? DC_PRED : intra_mode_end;
                disable_angular_prediction =
                    pcs_ptr->temporal_layer_index > 0 ? 1 : (bsize > 16) ? 1 : 0;
                angle_delta_candidate_count =
                    disable_angular_prediction ? 1 : use_angle_delta ? 5 : 1;
                angle_delta_shift = 1;
            } else if (pcs_ptr->intra_pred_mode == 6) {
                intra_mode_end =
                    (pcs_ptr->is_used_as_reference_flag == 0) ? DC_PRED : intra_mode_end;
                disable_angular_prediction =
                    pcs_ptr->temporal_layer_index > 0 ? 1 : (bsize > 16) ? 1 : 0;
                angle_delta_candidate_count = 1;
                angle_delta_shift           = 1;
            } else {
                if (pcs_ptr->slice_type == I_SLICE) {
                    intra_mode_end              = /*is_16_bit ? SMOOTH_H_PRED :*/ enable_paeth ? PAETH_PRED : enable_smooth ? SMOOTH_H_PRED : D67_PRED;
                    angle_delta_candidate_count = use_angle_delta ? 5 : 1;
                    disable_angular_prediction  = 0;
                    angle_delta_shift           = 1;
                } else if (pcs_ptr->temporal_layer_index == 0) {
                    intra_mode_end              = /*is_16_bit ? SMOOTH_H_PRED :*/ enable_paeth ? PAETH_PRED : enable_smooth ? SMOOTH_H_PRED : D67_PRED;
                    angle_delta_candidate_count = (bsize > 16) ? 1 : use_angle_delta ? 2 : 1;
                    disable_angular_prediction  = 0;
                    angle_delta_shift           = 3;
                } else {
                    intra_mode_end              = DC_PRED;
                    disable_angular_prediction  = 1;
                    angle_delta_candidate_count = 1;
                    angle_delta_shift           = 1;
                }
            }
            for (ois_intra_mode = intra_mode_start; ois_intra_mode <= intra_mode_end;
                 ++ois_intra_mode) {
                if (av1_is_directional_mode((PredictionMode)ois_intra_mode)) {
                    if (!disable_angular_prediction) {
                        for (angle_delta_counter = 0;
                             angle_delta_counter < angle_delta_candidate_count;
                             ++angle_delta_counter) {
                            int32_t angle_delta =
                                angle_delta_shift *
                                (angle_delta_candidate_count == 1
                                     ? 0
                                     : angle_delta_counter - (angle_delta_candidate_count >> 1));
                            int32_t p_angle = mode_to_angle_map[(PredictionMode)ois_intra_mode] +
                                              angle_delta * ANGLE_STEP;
                            // PRED
                            intra_prediction_open_loop(p_angle,
                                                       ois_intra_mode,
                                                       blk_origin_x,
                                                       blk_origin_y,
                                                       tx_size,
                                                       above_row,
                                                       left_col,
                                                       context_ptr);
                            // Distortion
                            ois_blk_ptr[ois_intra_count].distortion =
                                (uint32_t)nxm_sad_kernel( // Always SAD without weighting
                                    &(input_ptr->buffer_y[(input_ptr->origin_y + blk_origin_y) *
                                                              input_ptr->stride_y +
                                                          (input_ptr->origin_x + blk_origin_x)]),
                                    input_ptr->stride_y,
                                    &(context_ptr->me_context_ptr->sb_buffer[0]),
                                    BLOCK_SIZE_64,
                                    bsize,
                                    bsize);
                            // kepp track of best SAD
                            if (ois_blk_ptr[ois_intra_count].distortion <
                                best_intra_ois_distortion) {
                                best_intra_ois_index      = ois_intra_count;
                                best_intra_ois_distortion = ois_blk_ptr[ois_intra_count].distortion;
                            }
                            ois_blk_ptr[ois_intra_count].intra_mode       = ois_intra_mode;
                            ois_blk_ptr[ois_intra_count].valid_distortion = EB_TRUE;
                            ois_blk_ptr[ois_intra_count++].angle_delta    = angle_delta;
                        }
                    }
                } else {
                    // PRED
                    intra_prediction_open_loop(0,
                                               ois_intra_mode,
                                               blk_origin_x,
                                               blk_origin_y,
                                               tx_size,
                                               above_row,
                                               left_col,
                                               context_ptr);
                    // Distortion
                    ois_blk_ptr[ois_intra_count].distortion =
                        (uint32_t)nxm_sad_kernel( // Always SAD without weighting
                            &(input_ptr->buffer_y[(input_ptr->origin_y + blk_origin_y) *
                                                      input_ptr->stride_y +
                                                  (input_ptr->origin_x + blk_origin_x)]),
                            input_ptr->stride_y,
                            &(context_ptr->me_context_ptr->sb_buffer[0]),
                            BLOCK_SIZE_64,
                            bsize,
                            bsize);
                    // kepp track of best SAD
                    if (ois_blk_ptr[ois_intra_count].distortion < best_intra_ois_distortion) {
                        best_intra_ois_index      = ois_intra_count;
                        best_intra_ois_distortion = ois_blk_ptr[ois_intra_count].distortion;
                    }
                    ois_blk_ptr[ois_intra_count].intra_mode       = ois_intra_mode;
                    ois_blk_ptr[ois_intra_count].valid_distortion = EB_TRUE;
                    ois_blk_ptr[ois_intra_count++].angle_delta    = 0;
                }
            }
            ois_sb_results_ptr->best_distortion_index[pa_blk_index]     = best_intra_ois_index;
            ois_sb_results_ptr->total_ois_intra_candidate[pa_blk_index] = ois_intra_count;
        }
        pa_blk_index++;
    }
    return return_error;
}
